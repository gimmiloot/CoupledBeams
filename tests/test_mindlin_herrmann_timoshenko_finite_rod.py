"""Finite source-reduced Jang project gate; no plot-derived golden values."""
from dataclasses import replace
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from scripts.analysis import verify_mindlin_herrmann_timoshenko_single_rod as finite
from scripts.analysis import reproduce_mindlin_herrmann_timoshenko_literature as source
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as tim


@pytest.fixture(scope="module")
def setup():
    return finite.check_inputs()


@pytest.fixture(scope="module")
def computed(setup):
    config, model, length, _ = setup
    return finite.compute(config, model, length)


def test_selected_coefficients_preserve_project_section_exactly(setup):
    _, model, _, _ = setup
    s, p = model.section, model.coefficients
    assert model.variant == "project_jang_reduced_rectangular"
    assert s.K == model.mh_shear_factor == 5/6
    assert model.mh_inertia_factor == model.tim_rotary_factor == 1.
    assert p == {"C": s.EA/(1-s.nu**2), "H": s.K*s.G*s.inertia,
                 "m": s.rhoA, "j": s.rhoI, "B": s.EI, "S": s.KGA, "r": s.rhoI}


def test_project_kappa_gate_cannot_guess_or_accept_other_preset(setup):
    section = setup[1].section
    with pytest.raises(ValueError, match="accepted rectangular"):
        mh.project_jang_reduced_rectangular(replace(section, K=.9))
    original, _ = source.source_check()
    with pytest.raises(ValueError, match="not established"):
        source.make_model(original, "jang_2014_bare_isotropic")


def test_variation_signs_and_quasistatic_reduction_exact():
    # Exact central differentiation of the quadratic potential, no symbolic package.
    nu, c, h, b, shear = F(3, 10), F(100, 91), F(2, 7), F(3, 11), F(5, 13)
    jets = [F(2, 3), F(3, 7), F(4, 9), F(5, 11), F(6, 13), F(7, 15)]
    def potential(v):
        ux, contraction, cx, wx, theta, thetax = v
        return (c*(ux**2+2*nu*ux*contraction+contraction**2)+
                h*cx**2+b*thetax**2+shear*(wx-theta)**2)/2
    ux, contraction, cx, wx, theta, thetax = jets
    expected = [c*(ux+nu*contraction), c*(contraction+nu*ux), h*cx,
                shear*(wx-theta), -shear*(wx-theta), b*thetax]
    step = F(1, 17)
    for i, derivative in enumerate(expected):
        plus, minus = jets.copy(), jets.copy()
        plus[i] += step
        minus[i] -= step
        assert (potential(plus)-potential(minus))/(2*step) == derivative
    assert c*(ux**2-2*nu**2*ux**2+nu**2*ux**2) == ux**2  # EA=1


def test_full_clamp_is_four_field_dirichlet_in_selected_kinematics():
    # Ux=u-z*theta, Uz=w+z*c at both distinct face coordinates.
    z = F(1, 40)
    face = ((1, -z), (1, z))
    assert face[0][0]*face[1][1]-face[0][1]*face[1][0] != 0
    # Both 2x2 systems are invertible: vanishing point displacements force
    # u=theta=0 and w=c=0; this is a planar 4-field, not a full 3D assertion.


def test_low_k_acoustic_and_positive_contraction_cutoff(setup):
    _, model, _, _ = setup
    block = mh.blocks(model)[0]
    speed = math.sqrt(model.section.E/model.section.rho)
    assert block.temporal(0)[0]["group_velocity_m_s"] == pytest.approx(speed)
    assert block.cutoff_hz == math.sqrt(model.coefficients["C"]/model.coefficients["j"])/(2*math.pi)
    assert block.cutoff_hz > 0
    k = 1e-4/model.section.thickness
    assert math.sqrt(block.temporal(k)[0]["omega_squared"])/k == pytest.approx(speed, rel=2e-6)


@pytest.mark.parametrize("block", ["mh", "timoshenko"])
def test_state_basis_satisfies_independently_derived_state_and_resultants(setup, block):
    _, model, length, _ = setup
    omega, x = 2., np.array([0., .2, .7, 1.])*length
    values = mh.finite_state_basis(model, length, omega, x, block)
    gradients = mh.finite_state_basis(model, length, omega, x, block, 1)
    expected = np.einsum("ij,njk->nik", mh.harmonic_state_matrix(model, omega, block), values)
    np.testing.assert_allclose(gradients, expected, rtol=2e-12, atol=2e-12)
    p = model.coefficients
    np.testing.assert_allclose(gradients[:, 2], -p["m"]*omega**2*values[:, 0], atol=2e-14)
    if block == "mh":
        np.testing.assert_allclose(values[:, 2], p["C"]*(gradients[:, 0]+model.section.nu*values[:, 1]))
        np.testing.assert_allclose(gradients[:, 3], (p["C"]*(1-model.section.nu**2)-p["j"]*omega**2)*values[:, 1]+model.section.nu*values[:, 2], atol=2e-14)
    else:
        np.testing.assert_allclose(values[:, 2], p["S"]*(gradients[:, 0]-values[:, 1]))
        np.testing.assert_allclose(gradients[:, 3], -p["r"]*omega**2*values[:, 1]-values[:, 2], atol=2e-14)


def test_energy_block_structure_precedes_full_state_factorization(setup):
    model = setup[1]
    energy, mass = mh.energy_matrices(model)
    assert not np.any(energy[:3, 3:]) and not np.any(mass[:2, 2:])
    matrix = mh.full_harmonic_state_matrix(model, 2.)
    axial, bending = (0, 1, 4, 5), (2, 3, 6, 7)
    assert not np.any(matrix[np.ix_(axial, bending)])
    assert not np.any(matrix[np.ix_(bending, axial)])


@pytest.mark.parametrize("block", ["mh", "timoshenko"])
def test_primary_independent_frequencies_boundaries_energy_and_count(setup, computed, block):
    config = setup[0]
    result = computed[0][block]
    assert result["status"] == "PASS"
    assert len(result["roots"]) == result["completeness"]["upper_count"]
    assert result["below_search_count_bound"]["upper_count"] == 0
    assert result["roots"][-1]["guard"] is True
    for root in result["roots"]:
        assert root["converged"]
        assert root["bracket_determinants"][0]*root["bracket_determinants"][1] <= 0
        assert root["independent"]["relative_difference"] < config["policy"]["independent_frequency_relative_tol"]
        assert root["diagnostics"]["boundary_scaled_residual"] < config["policy"]["boundary_scaled_residual_tol"]
        assert root["diagnostics"]["energy_relative_error"] < config["policy"]["energy_relative_tol"]
    assert result["mass_orthogonality_max_error"] < config["policy"]["mass_orthogonality_tol"]


def test_young_completeness_inequality_exact_and_saturated(setup, computed):
    nu, eta = F(3, 10), F(18, 100)
    for ux, c in ((F(2, 3), F(-4, 7)), (F(0), F(5)), (F(2), F(3))):
        difference = ux**2+2*nu*ux*c+c**2-((1-eta)*ux**2+(1-nu**2/eta)*c**2)
        assert difference == eta*(ux+nu*c/eta)**2 >= 0
    assert computed[0]["mh"]["completeness"]["contraction_count"] == 0
    assert computed[0]["mh"]["completeness"]["upper_count"] == 7


def test_bending_matches_existing_project_basis_and_is_unchanged(setup, computed):
    _, model, length, _ = setup
    changed_axial = replace(model, mh_shear_factor=1.1, mh_inertia_factor=2.1)
    for root in computed[0]["timoshenko"]["roots"][:3]:
        omega = root["omega"]
        expected = tim.clamped_bending_columns(length, tim.timoshenko_spatial_basis(omega, model.section))
        matrix = np.array([expected["w"], length*expected["psi"]])
        matrix /= np.linalg.norm(matrix, axis=1)[:, None]
        assert abs(np.linalg.det(matrix)) < 1e-9
        np.testing.assert_array_equal(mh.finite_boundary_matrix(model, length, omega, "timoshenko"),
                                      mh.finite_boundary_matrix(changed_axial, length, omega, "timoshenko"))


def test_elementary_and_planar_rayleigh_love_are_exact_second_order(setup, computed):
    model, length = setup[1:3]
    section = model.section
    frequencies = computed[0]["hierarchy"]["axial_frequencies_hz"]
    for i, (wave, love) in enumerate(zip(frequencies["elementary"], frequencies["rayleigh_love_planar"]), 1):
        k = i*math.pi/length
        assert wave == pytest.approx(math.sqrt(section.E/section.rho)*k/(2*math.pi), rel=5e-15)
        assert love == pytest.approx(math.sqrt(section.EA*k*k/(section.rhoA+section.nu**2*section.rhoI*k*k))/(2*math.pi), rel=5e-15)
    assert all(r["H"] == 0 for r in computed[0]["hierarchy"]["reduced_boundary_checks"])


def test_combined_boundary_determinant_and_union_have_same_zeros(setup, computed):
    from scipy.linalg import block_diag
    model, length = setup[1:3]
    for omega in (2., 4., 7.):
        axial = mh.finite_boundary_matrix(model, length, omega)
        bending = mh.finite_boundary_matrix(model, length, omega, "timoshenko")
        assert np.linalg.det(block_diag(axial, bending)) == pytest.approx(np.linalg.det(axial)*np.linalg.det(bending), rel=2e-13)
    for block in ("mh", "timoshenko"):
        for root in computed[0][block]["roots"][:3]:
            matrices = [mh.finite_boundary_matrix(model, length, root["omega"], b) for b in ("mh", "timoshenko")]
            singular = np.linalg.svd(block_diag(*matrices), compute_uv=False)
            assert singular[-1]/singular[0] < setup[0]["policy"]["boundary_scaled_residual_tol"]
    for spectrum in computed[0]["hierarchy"]["combined_sorted_prefix"].values():
        assert [r["sorted_position"] for r in spectrum] == list(range(1, 13))
        assert all(r["family"] in ("bending", "axial_acoustic") for r in spectrum)


def test_no_optical_branch_is_claimed_in_the_low_inventory(setup, computed):
    result = computed[0]
    assert max(r["frequency_hz"] for r in result["mh"]["roots"]) < result["limits"]["contraction_cutoff_hz"]
    assert max(r["frequency_hz"] for r in result["timoshenko"]["roots"]) < result["limits"]["shear_cutoff_hz"]
    assert result["hierarchy_status"] == "PARTIAL_PASS"  # full resolved c-clamp not represented in reduced models


def test_nu_zero_separate_axial_basis_control(setup):
    model, length = setup[1:3]
    zero = mh.project_jang_reduced_rectangular(replace(model.section, nu=0.))
    omega = math.pi/length*math.sqrt(zero.section.E/zero.section.rho)
    singular = np.linalg.svd(mh.finite_boundary_matrix(zero, length, omega), compute_uv=False)
    assert singular[-1]/singular[0] < 1e-12


def test_source_alternative_does_not_select_production_factors(setup):
    config, model, _, checked = setup
    assert "fernandes" in checked
    fixture, _ = source.source_check()
    assert fixture["sources"]["fernandes"]["source_expressions"]["K_r1"] == "12/pi^2"
    assert model.mh_shear_factor != 12/math.pi**2
    assert fixture["production_prescription"]["adopted"] is False
    assert fixture["sources"]["jang"]["factors"]["kappa_b_numeric"] is None


def test_cache_identity_reuse_and_integrity(setup, computed, tmp_path, monkeypatch):
    monkeypatch.setattr(finite, "compute", lambda *args: computed)
    args = ["--compute", "--output-dir", str(tmp_path)]
    assert finite.main(args) == 0
    def forbidden(*args):
        raise AssertionError("Matching cache recomputed roots")
    monkeypatch.setattr(finite, "compute", forbidden)
    assert finite.main(args) == 0
    directory = next(p for p in tmp_path.iterdir() if p.is_dir())
    (directory/"mh_modes.csv").write_text("corrupted", encoding="utf-8")
    with pytest.raises(ValueError, match="Stale or changed"):
        finite.main(args)


def test_above_cutoff_basis_and_unbounded_transfer_are_rejected(setup):
    model, length = setup[1:3]
    with pytest.raises(ValueError, match="below optical"):
        mh.finite_state_basis(model, length, 2*math.pi*mh.blocks(model)[0].cutoff_hz, 0.)
    with pytest.raises(ArithmeticError, match="step budget"):
        mh.transfer_boundary_matrix(model, length, 2., max_steps=1)
