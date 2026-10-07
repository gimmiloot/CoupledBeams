"""Source energy, dispersion and provenance checks; no graph-derived goldens."""
from dataclasses import replace
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from scripts.analysis import reproduce_mindlin_herrmann_timoshenko_literature as cli
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as timo


@pytest.fixture(scope="module")
def fixture():
    return cli.source_check()[0]


@pytest.fixture(scope="module", params=("rucka_2010", "jang_2014_bare_isotropic"))
def model(request, fixture):
    return cli.make_model(fixture, request.param, 5/6)


def test_exact_rectangle_integrals_and_centroidal_kinetic_cross_terms():
    b, h = F(1, 50), F(1, 500)
    def integral(power, length):
        return ((length/2)**(power+1)-(-length/2)**(power+1))/(power+1)
    integrated = {"A": b*h, "Qy": integral(1, b)*h, "Qz": b*integral(1, h),
        "Iyz": integral(1, b)*integral(1, h), "Iy": b*integral(2, h), "Iz": h*integral(2, b)}
    integrated["Ip"] = integrated["Iy"]+integrated["Iz"]
    assert mh.rectangle_moments(b, h) == integrated
    # Integrating (u_t-z*theta_t)^2+(w_t+z*c_t)^2 from source fields.
    rho = F(2700)
    assert -rho*integrated["Qz"] == rho*integrated["Qz"] == 0
    offset = F(1, 1000)
    shifted_qz = b*((offset+h/2)**2-(offset-h/2)**2)/2
    assert -rho*shifted_qz != 0 and rho*shifted_qz != 0


def test_dimensions_energy_per_length_and_pde_force():
    units = mh.COEFFICIENT_UNITS
    add = lambda a, b: tuple(x+y for x, y in zip(a, b))
    force = (1, 1, -2)
    assert units["C"] == units["S"] == force
    assert add(units["H"], (0, -2, 0)) == force
    assert add(units["B"], (0, -2, 0)) == force
    assert add(units["m"], (0, 2, -2)) == force
    assert add(units["j"], (0, 0, -2)) == force
    assert add(units["r"], (0, 0, -2)) == force
    assert add(units["m"], (0, 1, -2)) == add(units["C"], (0, -1, 0))


def test_source_field_energy_normal_reduction_exact():
    nu, ea = F(33, 100), F(400)
    cst = ea/(1-nu**2)
    # Quasistatic stationary contraction c=-nu*u_x: C*(1-nu^2)=EA.
    ux = F(7, 11)
    contraction = -nu*ux
    normal_energy = cst*(ux**2+contraction**2+2*nu*ux*contraction)/2
    assert normal_energy == ea*ux**2/2
    # Two separate stress reductions: MH keeps sigma_zz, Timo suppresses it.
    assert cst > ea  # hence C*I cannot silently replace the bending EI term
    compliance = ((1/ea, -nu/ea), (-nu/ea, 1/ea))  # sigma_yy=0, x/z retained
    reduced = ((cst, nu*cst), (nu*cst, cst))
    product = [[sum(reduced[i][k]*compliance[k][j] for k in (0, 1))
                for j in (0, 1)] for i in (0, 1)]
    assert product == [[1, 0], [0, 1]]
    # Additional sigma_zz=0 yields epsilon_xx=sigma_xx/E exactly.
    sx = F(13, 7)
    assert compliance[0][0]*sx == sx/ea


def test_source_matrices_symmetry_positive_mass_and_no_cross_block(model):
    elastic, mass = mh.energy_matrices(model)
    np.testing.assert_array_equal(elastic, elastic.T)
    assert np.all(np.linalg.eigvalsh(elastic) > 0)
    assert np.all(np.diag(mass) > 0)
    for k in (0, 1, 100, 1000):
        stiffness, actual_mass = mh.fourier_matrices(model, k)
        np.testing.assert_array_equal(stiffness, stiffness.conj().T)
        np.testing.assert_array_equal(stiffness[:2, 2:], np.zeros((2, 2)))
        np.testing.assert_array_equal(actual_mass[:2, 2:], np.zeros((2, 2)))
    # Adding a constitutive mixed entry creates coupling, it is not discarded.
    altered = elastic.copy()
    altered[0, 4] = altered[4, 0] = 1.
    d = mh.strain_operator(2.)
    assert (d.conj().T@altered@d)[0, 3] != 0


def test_exact_family_mapping_and_actual_published_non_equivalence(fixture):
    params = fixture["cases"]["rucka"]["parameters_si"]
    jang = cli.make_model(fixture, "jang_2014_bare_isotropic", 5/6, params)
    mapped = mh.source_model(params, mh_shear_factor=5/6, mh_inertia_factor=1,
        tim_shear_factor=5/6, tim_rotary_factor=1, variant="rucka_family_mapped_not_fitted")
    for k in (0, 10, 200):
        for a, b in zip(mh.fourier_matrices(jang, k), mh.fourier_matrices(mapped, k)):
            np.testing.assert_array_equal(a, b)
    fitted = cli.make_model(fixture, "rucka_2010")
    assert fitted.coefficients["j"] != jang.coefficients["j"]
    # Even an unstated Jang kappa cannot remove Rucka K_MH2=2.1 vs 1.
    assert mh.blocks(fitted)[0].cutoff_hz*math.sqrt(2.1) == pytest.approx(mh.blocks(jang)[0].cutoff_hz)


def test_low_frequency_acoustic_and_bending_limits(model, fixture):
    axial, bending = mh.blocks(model)
    p = model.coefficients
    for kh in fixture["numerical_contract"]["low_k_dimensionless_values"]:
        k = kh/model.section.thickness
        r = axial.temporal(k)[0]
        assert math.sqrt(r["omega_squared"])/k == pytest.approx(math.sqrt(model.section.E/model.section.rho), rel=2e-6)
        r = bending.temporal(k)[0]
        assert math.sqrt(r["omega_squared"])/k**2 == pytest.approx(math.sqrt(p["B"]/p["m"]), rel=2e-6)
    assert axial.temporal(0)[0]["group_velocity_m_s"] == pytest.approx(
        math.sqrt(model.section.E/model.section.rho),
        rel=fixture["numerical_contract"]["eigenvalue_relative_tol"])
    assert bending.temporal(0)[0]["group_velocity_m_s"] == 0


def test_contraction_and_shear_cutoff_and_branch_count(model):
    p = model.coefficients
    for block, expected in zip(mh.blocks(model), (math.sqrt(p["C"]/p["j"]), math.sqrt(p["S"]/p["r"]))):
        assert block.temporal(0)[1]["omega_squared"] == pytest.approx(expected**2)
        assert sum(r["state"] == "PROPAGATING" for r in block.spatial(.99*block.cutoff_hz)) == 1
        assert sum(r["state"] == "PROPAGATING" for r in block.spatial(1.01*block.cutoff_hz)) == 2
        assert block.spatial(block.cutoff_hz)[1]["wavenumber_per_m"] == 0


def test_independent_variational_eigenproblem_and_group_derivative(model, fixture):
    checks = cli.verify_model(model, fixture["numerical_contract"])
    assert checks["status"] == "PASS", checks["maxima"]


def test_project_timoshenko_exact_self_terms_and_spatial_roots(model, fixture):
    # Source Rucka rotary correction is source-specific; reset it only for this
    # explicit project limit, keeping the existing project's own section input.
    project = replace(model, tim_rotary_factor=1.)
    p = project.coefficients
    assert p["B"] == project.section.EI
    assert p["S"] == project.section.KGA
    assert p["r"] == project.section.rhoI
    block = mh.blocks(project)[1]
    for multiplier in (.1, .9, 1., 1.1, 2.):
        f = block.cutoff_hz*multiplier
        expected = timo.timoshenko_spatial_basis(2*math.pi*f, project.section)
        actual = sorted(-r["k_squared_per_m2"] for r in block.spatial(f))
        np.testing.assert_allclose(actual, sorted((expected.z_a, expected.z_b)),
            rtol=fixture["numerical_contract"]["project_spatial_root_relative_tol"], atol=1e-8)


def test_full_operator_determinant_factorization(model):
    k = 100.
    stiffness, mass = mh.fourier_matrices(model, k)
    inv = 1/np.sqrt(np.diag(mass))
    matrix = stiffness*inv[:, None]*inv[None, :]
    scale = np.linalg.norm(matrix, 2)
    matrix = matrix/scale-.37*np.eye(4)
    full = np.linalg.det(matrix)
    product = np.linalg.det(matrix[:2, :2])*np.linalg.det(matrix[2:, 2:])
    assert full == pytest.approx(product, rel=1e-13, abs=1e-16)


def test_bending_identical_when_only_mh_coefficients_change(model):
    variant = replace(model, mh_shear_factor=.7, mh_inertia_factor=1.6)
    for k in (0., 10., 500.):
        for a, b in zip(mh.fourier_matrices(model, k), mh.fourier_matrices(variant, k)):
            np.testing.assert_array_equal(a[2:, 2:], b[2:, 2:])
        assert mh.blocks(model)[1].temporal(k) == mh.blocks(variant)[1].temporal(k)


def test_rucka_printed_mode_statements(fixture):
    source = cli.make_model(fixture, "rucka_2010")
    for statement in fixture["sources"]["rucka"]["statements"]:
        block = mh.blocks(source)[statement["block"] == "tim"]
        assert block.cutoff_hz > statement["frequency_interval_hz"][1]
        for f in statement["frequency_interval_hz"]:
            assert sum(r["state"] == "PROPAGATING" for r in block.spatial(f)) == statement["propagating_count"]


def test_jang_kappa_is_never_guessed(fixture, tmp_path):
    assert fixture["sources"]["jang"]["factors"]["kappa_b_numeric"] is None
    with pytest.raises(ValueError, match="not established"):
        cli.make_model(fixture, "jang_2014_bare_isotropic")
    assert cli.main(["--compute", "--case", "jang", "--output-dir", str(tmp_path)]) == 2
    summary = next(tmp_path.glob("*/jang/summary.json"))
    assert json.loads(summary.read_text())["status"] == "SOURCE_NUMERIC_CONFIG_UNRESOLVED"


def test_cache_artifact_integrity_input_identity_and_plot_only(fixture, tmp_path, monkeypatch):
    argv = ["--case", "rucka", "--output-dir", str(tmp_path)]
    assert cli.main(["--compute", *argv]) == 0
    current = json.loads((tmp_path/"current.json").read_text())
    directory = __import__('pathlib').Path(current["directory"])/"rucka"
    fixture, checked = cli.source_check()
    assert cli.provenance(checked, "5/6")["fingerprint"] != cli.provenance(checked, "0.9")["fingerprint"]
    def forbidden(*a, **kw):
        raise AssertionError("Plot-only/cache reuse invoked computation")
    monkeypatch.setattr(cli, "compute", forbidden)
    monkeypatch.setattr(cli, "verify_model", forbidden)
    monkeypatch.setattr(mh, "source_model", forbidden)
    assert cli.main(["--compute", *argv]) == 0
    assert cli.main(["--plot-only", *argv]) == 0
    assert (directory/"rucka.png").exists()
    with pytest.raises(ValueError, match="Stale"):
        cli.validate_bundle(directory, "other-inputs")
    with (directory/"dispersion.csv").open("a") as f:
        f.write("corruption\n")
    with pytest.raises(ValueError, match="Artifact changed"):
        cli.validate_bundle(directory, current["fingerprint"])


def test_fixture_and_equation_changes_invalidate_provenance(fixture, monkeypatch):
    _, checked = cli.source_check()
    original = cli.provenance(checked, "5/6")["fingerprint"]
    original_sha = cli.sha
    monkeypatch.setattr(cli, "sha", lambda p: "changed-input-bytes" if p == cli.FIXTURE else original_sha(p))
    assert cli.provenance(checked, "5/6")["fingerprint"] != original
    monkeypatch.setattr(mh, "EQUATIONS_VERSION", "new-equations")
    with pytest.raises(ValueError, match="version mismatch"):
        cli.source_check()


def test_modified_si_input_cannot_keep_source_reproduction_label(fixture, tmp_path, monkeypatch):
    changed = json.loads(json.dumps(fixture))
    changed["cases"]["rucka"]["parameters_si"]["h"] *= 2
    path = tmp_path/"changed.json"
    path.write_text(json.dumps(changed), encoding="utf-8")
    monkeypatch.setattr(cli, "FIXTURE", path)
    with pytest.raises(ValueError, match="Source/SI parameter conflict"):
        cli.source_check()


def test_rectangular_candidate_cannot_bypass_two_source_and_convention_gate(fixture):
    prescription = fixture["production_prescription"]
    assert prescription["adopted"] is False
    assert prescription["code_preset"] is None
    assert prescription["status"] == "PRODUCTION_MH_COEFFICIENTS_UNRESOLVED"
    assert prescription["mapping_status"] == "RECTANGULAR_MH_PRESET_MAPPING_UNRESOLVED"
    assert prescription["direct_rectangular_factor_confirmations"] == ["ng"]
    assert fixture["source_variants"]["fernandes_2022_rectangular"]["source"] is None
    assert "NOT_A_SECOND_RECTANGULAR" in fixture["sources"]["elishakoff_tharu"]["prescription_status"]


def test_ng_printed_factors_map_to_stiffness_and_inertia_without_squaring(fixture):
    entry = fixture["sources"]["ng"]
    printed = entry["source_expressions"]
    assert printed["S1"] == "12/pi^2"
    assert printed["S2_j"] == "S1*((1+nu_j)/(0.87+1.12*nu_j))^2"
    assert "mu_j*I_j*S1*phi_j,xx" in printed["pde_contraction"]
    assert "rho_j*I_j*S2,j*phi_j,tt" in printed["pde_contraction"]
    assert entry["factor_mapping"]["S1"]["project_symbol"] == "K_MH1"
    assert entry["factor_mapping"]["S2_j"]["project_symbol"] == "K_MH2"
    # Arithmetic of a printed candidate, not an adopted project helper/default.
    for nu in (0., .3, .49):
        s1 = 12/math.pi**2
        ratio = (1+F(str(nu)))/(F('0.87')+F('1.12')*F(str(nu)))
        s2 = s1*((1+nu)/(.87+1.12*nu))**2
        assert s1 > 0 and math.isfinite(s2) and s2 > 0
        assert s2/s1 == pytest.approx(float(ratio**2), rel=5e-15)


def test_ng_lame_normal_block_is_not_current_reduced_energy(fixture):
    assert "NOT_EQUIVALENT" in fixture["sources"]["ng"]["normal_block"]["project_normal_mapping"]
    for nu in (F(0), F(3, 10), F(49, 100)):
        # E=A=1; direct reconstruction from Ng (2) and printed Lame definitions.
        mu = 1/(2*(1+nu))
        lame = nu/((1+nu)*(1-2*nu))
        d, cross = 2*mu+lame, lame
        current_c = 1/(1-nu**2)
        assert cross/d == nu/(1-nu)
        assert d-cross**2/d == current_c
        assert current_c-(nu*current_c)**2/current_c == 1
        assert ((d, cross) == (current_c, nu*current_c)) == (nu == 0)
        # Normalized coupling is invariant under independent DOF scalings.
        su, sc = F(7, 3), F(5, 2)
        scaled_cross_squared = (cross*su*sc)**2/((d*su**2)*(d*sc**2))
        assert scaled_cross_squared == (cross/d)**2
        if nu:
            assert scaled_cross_squared != nu**2


def test_new_source_identity_hashes_and_citation_chain(fixture):
    _, checked = cli.source_check()
    for name, count in (("ng", 41), ("elishakoff_tharu", 100)):
        entry = fixture["sources"][name]
        assert checked[name]["sha256"] == cli.sha(cli.ROOT/entry["path"])
        assert entry["pdf_page_count"] == count
    assert fixture["sources"]["ng"]["citation_chain"][0]["reference"] == 37
    preprint = fixture["sources"]["elishakoff_tharu"]
    assert preprint["year"] is None and preprint["doi"] is None
    assert "not peer reviewed" in preprint["version"]
    assert [r["reference"] for r in preprint["citation_chain"]] == [10, 15]


def test_unadopted_candidate_does_not_override_rucka_or_guess_jang(fixture):
    rucka = cli.make_model(fixture, "rucka_2010")
    assert (rucka.mh_shear_factor, rucka.mh_inertia_factor,
            rucka.section.K, rucka.tim_rotary_factor) == (1.1, 2.1, .95, 12*.95/math.pi**2)
    with pytest.raises(ValueError, match="not established"):
        cli.make_model(fixture, "jang_2014_bare_isotropic")
    jang = cli.make_model(fixture, "jang_2014_bare_isotropic", 5/6)
    assert (jang.mh_shear_factor, jang.mh_inertia_factor,
            jang.section.K, jang.tim_rotary_factor) == (5/6, 1., 5/6, 1.)
    assert jang.coefficients["r"] == jang.section.rhoI


@pytest.mark.parametrize("factor", [0, -1, float('nan')])
def test_invalid_correction_factors_rejected(fixture, factor):
    with pytest.raises(ValueError):
        cli.make_model(fixture, "jang_2014_bare_isotropic", factor)
