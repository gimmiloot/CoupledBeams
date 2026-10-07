"""Seven-field reduced action audit; no PDE trajectories or mode classification."""
from fractions import Fraction
import math

import numpy as np
import pytest

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import yartsev_ch2_monoclinic_rod as book
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section
from scripts.analysis import verify_weakly_nonlinear_spatial_rod as audit


@pytest.fixture(scope="module")
def model():
    # The lightweight exact derivation is generated once per test session.
    return rod.derive_polynomials()


@pytest.fixture(scope="module")
def coefficients():
    G = 1 / 2.6
    material = book.BookMaterial(E1_real=1., E2_real=1., G12_real=G,
        G13_real=G, G23_real=G, nu12=.3, rho=1., eta1=0., eta2=0.,
        eta12=0., eta13=0., eta23=0.)
    point = book.make_rod_point(0., material=material,
        geometry=book.Geometry(a=.20, b=.05, length=1., shear_factor=5/6))
    return rod.RodCoefficients.rectangular(1., 1., .3, .20, .05,
                                          float(point.torsion.C_T.real))


def test_field_and_jet_contract_are_explicit():
    assert rod.FIELD_ORDER == ("u", "w", "v", "Phi", "psi", "theta", "c")
    assert rod.JET_ORDER == ("q", "qs", "qt", "qss", "qst", "qtt")
    with pytest.raises(ValueError, match="seven finite"):
        rod.FieldJet(*(np.zeros(6) for _ in range(6)))


def test_exact_sparse_ring_keeps_material_constants_degree_zero():
    u, m = rod.Polynomial.symbol("u"), rod.Polynomial.symbol("m")
    polynomial = Fraction(1, 2) * m * u**2
    assert polynomial.homogeneous(2) == polynomial
    assert polynomial.derivative("u") == m * u
    assert polynomial.total_derivative("s") == m * u * rod.Polynomial.symbol("u_s")
    assert not (u**5)
    assert rod.Polynomial.deserialize(polynomial.serialize()) == polynomial


@pytest.mark.parametrize("field", range(7), ids=rod.FIELD_ORDER)
@pytest.mark.parametrize("degree", (1, 2, 3))
def test_twenty_one_exact_independent_action_balance_identities(model, field, degree):
    assert not (model.residual_a[field] - model.residual_b[field]).homogeneous(degree)


def test_linear_operator_recovers_both_bending_axes_and_generalized_torsion(model):
    p = model.symbols
    expected = (
        p["m"]*p["u_tt"]-p["C"]*(p["u_ss"]+p["nu"]*p["c_s"]),
        p["m"]*p["w_tt"]-p["S"]*(p["w_ss"]-p["theta_s"]),
        p["m"]*p["v_tt"]-p["S"]*(p["v_ss"]-p["psi_s"]),
        (p["jp"]+p["jb"])*p["Phi_tt"]-p["CT"]*p["Phi_ss"],
        p["jb"]*p["psi_tt"]-p["Bb"]*p["psi_ss"]-p["S"]*(p["v_s"]-p["psi"]),
        p["jp"]*p["theta_tt"]-p["Bp"]*p["theta_ss"]-p["S"]*(p["w_s"]-p["theta"]),
        p["jp"]*p["c_tt"]-p["H"]*p["c_ss"]+p["C"]*(p["c"]+p["nu"]*p["u_s"]),
    )
    assert all(not (a.homogeneous(1)-b) for a, b in zip(model.residual_a, expected))


def test_pure_axial_subspace_contains_no_quadratic_or_cubic_terms(model):
    zero = audit.restriction(("w", "v", "Phi", "psi", "theta"))
    for index, polynomial in enumerate(model.residual_a):
        restricted = polynomial.substitute(zero)
        assert restricted == restricted.homogeneous(1)
        if index not in (0, 6):
            assert not restricted


def test_planar_invariant_subspace_against_separate_trigonometric_balances(model):
    p = model.symbols
    theta, c = p["theta"], p["c"]
    cos, sin = 1-theta**2/2, theta-theta**3/6
    gamma1 = ((1+p["u_s"])*cos+p["w_s"]*sin-1).truncate(3)
    gamma2 = (-(1+p["u_s"])*sin+p["w_s"]*cos).truncate(3)
    N, Q = p["C"]*(gamma1+p["nu"]*c), p["S"]*gamma2
    expected = (
        p["m"]*p["u_tt"]-(N*cos-Q*sin).truncate(3).total_derivative("s"),
        p["m"]*p["w_tt"]-(N*sin+Q*cos).truncate(3).total_derivative("s"),
        rod.Polynomial(), rod.Polynomial(), rod.Polynomial(),
        (p["jp"]*(1+c)**2*p["theta_t"]).total_derivative("t")
            -p["Bp"]*p["theta_ss"]-((1+gamma1)*Q-gamma2*N).truncate(3),
        p["jp"]*p["c_tt"]-p["H"]*p["c_ss"]+p["C"]*(c+p["nu"]*gamma1)
            -p["jp"]*(1+c)*p["theta_t"]**2,
    )
    zero = audit.restriction(("v", "Phi", "psi"))
    assert all(not (a.substitute(zero)-b.truncate(3)) for a, b in zip(model.residual_a, expected))


def test_reflection_parity_in_action_residuals_and_boundary_covectors(model):
    signs = (1, 1, -1, -1, -1, 1, 1)
    substitution = {name+suffix: sign*model.symbols[name+suffix]
        for name, sign in zip(rod.FIELD_ORDER, signs)
        for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")}
    assert not (model.T4.substitute(substitution)-model.T4)
    assert not (model.V4.substitute(substitution)-model.V4)
    for outputs in (model.residual_a, model.flux_a):
        assert all(not (p.substitute(substitution)-sign*p) for p, sign in zip(outputs, signs))


def test_boundary_fluxes_and_coordinate_moment_transform_exact(model):
    assert all(not (a-b) for a, b in zip(model.flux_a, model.flux_b))
    p = model.symbols
    assert model.flux_a[6] == p["H"]*p["c_s"]
    # The coordinate moment differs from a bare body-vector moment away from zero.
    bare = (p["CT"]*model.chi_b[0], -p["Bb"]*model.chi_b[1], p["Bp"]*model.chi_b[2])
    assert any((flux-body).truncate(3) for flux, body in zip(model.flux_a[3:6], bare))


def test_polynomial_energy_identity_has_zero_exact_coefficients(model):
    p = model.symbols
    power = sum((p[name+"_t"]*flux for name, flux in zip(rod.FIELD_ORDER, model.flux_a)), rod.Polynomial())
    work = sum((p[name+"_t"]*residual for name, residual in zip(rod.FIELD_ORDER, model.residual_a)), rod.Polynomial())
    assert not ((model.T4+model.V4).total_derivative("t")-power.total_derivative("s")-work)


def test_negative_controls_detect_omitted_Jr_and_contraction_inertia(model):
    wrong_moments = (model.body_b[0], -model.body_b[1], model.body_b[2])
    assert all((a-b).truncate(3) for a, b in zip(model.residual_a[3:6], wrong_moments))
    p = model.symbols
    omitted = (p["jp"]*(1+p["c"])*(model.omega_b[0]**2+model.omega_b[2]**2)).truncate(3)
    assert omitted
    # nu=0 does not remove the independent contraction excitation by rotation.
    assert omitted.substitute({"nu": 0}) == omitted


def test_mass_matrix_symmetric_positive_only_in_declared_small_neighborhood(model, coefficients):
    for rotation in (np.zeros(3), np.array([.02, -.03, .04])):
        for contraction in (-.1, 0., .1):
            q = np.zeros(7); q[3:6] = rotation; q[6] = contraction
            jet = rod.FieldJet(q, *(np.zeros(7) for _ in range(5)))
            mass = rod.quartic_mass_matrix(jet, coefficients, model)
            np.testing.assert_allclose(mass, mass.T, rtol=0, atol=1e-18)
            assert np.linalg.eigvalsh(mass).min() > 0
    zero = rod.FieldJet(*(np.zeros(7) for _ in range(6)))
    expected = np.diag([coefficients.m]*3+[coefficients.jp+coefficients.jb,
        coefficients.jb, coefficients.jp, coefficients.jp])
    np.testing.assert_allclose(rod.quartic_mass_matrix(zero, coefficients, model), expected, rtol=1e-15, atol=0)


def test_exact_rotation_and_right_jacobian_directional_derivatives():
    a = np.array([.08, -.03, .05]); direction = np.array([.2, .4, -.1])
    R, J, dR, dJ = rod.rotation_and_right_jacobian(a, direction)
    np.testing.assert_allclose(R.T@R, np.eye(3), rtol=0, atol=3e-16)
    assert np.linalg.det(R) == pytest.approx(1., abs=3e-16)
    # Independent matrix exponential gives R; finite differences only verify
    # the evaluator's analytic derivative, never group velocity or PDE order.
    from scipy.linalg import expm
    np.testing.assert_allclose(R, expm(rod.skew(a)), rtol=0, atol=3e-16)
    step = 2e-6
    plus = rod.rotation_and_right_jacobian(a+step*direction)
    minus = rod.rotation_and_right_jacobian(a-step*direction)
    np.testing.assert_allclose(dR, (plus[0]-minus[0])/(2*step), rtol=2e-8, atol=2e-10)
    np.testing.assert_allclose(dJ, (plus[1]-minus[1])/(2*step), rtol=2e-8, atol=2e-10)
    np.testing.assert_allclose(R.T@dR, rod.skew(J@direction), rtol=0, atol=2e-16)


def test_rigid_motion_exact_zero_strain_and_order_correct_cubic_flux(model, coefficients):
    q = np.zeros(7); q[3:6] = (.07, -.02, .05)
    a = np.diag([1., -1., 1.])@q[3:6]
    R, _ = rod.rotation_and_right_jacobian(a)
    qs = np.zeros(7); qs[:3] = R[:, 0]-np.array([1., 0., 0.])
    jet = rod.FieldJet(q, qs, *(np.zeros(7) for _ in range(4)))
    full = rod.full_evaluate(jet, coefficients)
    assert np.linalg.norm(full["Gamma"]) < 2e-16
    assert abs(full["V"]) < 1e-30
    assert np.linalg.norm(full["flux"]) < 2e-17
    # Finite-amplitude cubic flux need not vanish exactly; the exact algebraic
    # order check substitutes the exponential series and truncates at degree3.
    assert audit.exact_checks(model)["rigid_flux"]["status"] == "PASS"


def test_quasistatic_effective_EA_and_low_frequency_speed(coefficients):
    E, rho, A, nu = 1., 1., .2*.05, .3
    assert coefficients.C*(1-nu**2) == pytest.approx(E*A, rel=2e-16)
    assert math.sqrt(coefficients.C*(1-nu**2)/coefficients.m) == pytest.approx(math.sqrt(E/rho), rel=2e-16)


def test_linear_state_maps_independent_old_operators_after_axis_permutation(model, coefficients):
    section = rectangular_section(E=1., nu=.3, rho=1., width=.2, thickness=.05, K=5/6)
    old = mh.project_jang_reduced_rectangular(section)
    state, _ = audit.linear_state_from_action(coefficients, model)
    ids = {"mh": (0, 6, 7, 13), "timoshenko": (1, 5, 8, 12)}
    for omega in (.5, 2., 9.):
        for block, indices in ids.items():
            np.testing.assert_allclose(state(omega)[np.ix_(indices, indices)],
                mh.harmonic_state_matrix(old, omega, block), rtol=2e-15, atol=1e-15)
    # The swapped Chapter-2 bending axis must be the width-cubed moment.
    assert coefficients.jb == pytest.approx(.05*.2**3/12, rel=1e-15)
    assert coefficients.jb != coefficients.jp
    assert coefficients.CT != pytest.approx(section.G*(coefficients.jp+coefficients.jb), rel=1e-3)


@pytest.mark.parametrize("beta", (0., 45., 90.))
def test_full_joint_fourteen_conditions_and_dual_physical_moments(beta):
    transforms = audit.frame_transforms(beta)
    joint = audit.joint_operator(beta)
    assert joint.shape == (14, 28)
    assert np.linalg.matrix_rank(joint) == 14
    np.testing.assert_allclose(joint@joint.T, 2*np.eye(14), rtol=0, atol=2e-15)
    rng = np.random.default_rng(9831)
    for T in transforms:
        np.testing.assert_allclose(T.T@T, np.eye(7), rtol=0, atol=5e-16)
        assert T[6, 6] == 1. and np.count_nonzero(T[6]) == 1
        for _ in range(12):
            f, variation = rng.normal(size=(2, 7))
            assert f@(T.T@variation) == pytest.approx((T@f)@variation, rel=2e-14, abs=2e-14)


def test_section_rotation_clamp_does_not_add_centerline_slope_constraint(coefficients):
    q = np.zeros(7); qs = np.array([.02, .03, -.04, .0, .0, .0, .0])
    jet = rod.FieldJet(q, qs, *(np.zeros(7) for _ in range(4)))
    full = rod.full_evaluate(jet, coefficients)
    assert full["N"][1] == pytest.approx(coefficients.S*.03)
    assert full["N"][2] == pytest.approx(coefficients.S*-.04)
    assert np.all(jet.q == 0) and np.any(jet.qs[:3] != 0)


@pytest.fixture(scope="module")
def completed_audit():
    # A validated saved audit is reused; a missing audit makes one bounded run.
    # Derivation is lightweight and is never executed separately for every test.
    result, bundle, performance = audit.compute()
    return result, bundle


def test_all_named_coefficients_have_exact_dimensional_consistency(model):
    assert audit.exact_checks(model)["dimensions"] == {
        "status": "PASS", "bad_monomials": [], "units": "kg,m,s exponents; densities per original length"}


def test_supplied_equations_compared_independently_and_exactly(model):
    comparison = audit.compare_supplied(model, audit.ROOT/"data/input/cubic_seven_field_expansion_supplied.md")
    assert comparison["status"] == "MATCH"
    assert len(comparison["comparisons"]) == 21
    assert all(not part["difference"]["terms"] for part in comparison["comparisons"])


@pytest.mark.parametrize("expression", (r"u^{5}", r"u u u u u", "unknown u", r"\sqrt{u}"))
def test_supplied_adapter_rejects_unknown_or_silently_truncated_inputs(expression):
    with pytest.raises(ValueError):
        audit.SuppliedPolynomialParser(expression).parse()


def test_supplied_adapter_detects_a_wrong_printed_coefficient(model, tmp_path):
    source = audit.ROOT/"data/input/cubic_seven_field_expansion_supplied.md"
    wrong = tmp_path/"wrong.md"
    wrong.write_text(source.read_text(encoding="utf8").replace("m u_{tt}", "2 m u_{tt}", 1), encoding="utf8")
    result = audit.compare_supplied(model, wrong)
    assert result["status"] == "DISCREPANCY"
    assert result["comparisons"][0]["difference"]["terms"]


def test_manufactured_amplitude_residuals_have_fourth_order(completed_audit):
    result, _ = completed_audit
    rows = result["amplitude"]["rows"]
    assert [row["epsilon_a"] for row in rows] == [.04,.02,.01,.005,.0025]
    assert all(3.7 < row["aggregate_order"] < 4.3 for row in rows[1:])
    assert result["amplitude"]["final_above_floor"]
    assert len(rows[-1]["component_errors"]) == 7
    assert rows[-1]["component_above_reporting_floor"][6] is False


def test_certified_linear_prefix_does_not_hide_axial_family(completed_audit):
    result, _ = completed_audit
    linear = result["linear"]
    assert len(linear["full_first6_guard7"]) == 7
    assert linear["completeness"]["found_count"] == linear["completeness"]["upper_count"] == 9
    assert linear["completeness"]["ceiling_omega"] == pytest.approx(2.7)
    assert linear["completeness"]["outplane"]["upper_count"] == 2
    assert len(linear["mh_family_first3"]) == 3
    assert all(row["relative_frequency_difference"] < 2e-8 for row in linear["mh_family_first3"])
    assert len(linear["profile_comparison"]) == 10


def test_straight_split_profiles_and_c_R_are_transparent(completed_audit):
    result, _ = completed_audit
    splits = result["linear"]["straight_split"]
    assert [case["split"] for case in splits] == [.5,.35]
    for case in splits:
        assert len(case["mh_family_first3"]) == 3
        assert len(case["profile_comparison"]) == 10
        assert all(max(profile["joint_scaled_residual"]) < 1e-9 for profile in case["profile_comparison"])
    assert max(case["action_additivity_error"] for case in result["joint"]["action_splits"]) < 1e-12
    assert max(abs(x) for x in result["joint"]["nonlinear_interface_balance"]) < 1e-12


def test_nonlinear_audit_does_not_claim_unavailable_same_clamp_angular_reference(completed_audit):
    result, _ = completed_audit
    assert result["linear"]["book_slope_clamp_reference_used"] is False
    assert result["linear"]["nonzero_angle_old_same_clamp_reference"].startswith("UNAVAILABLE")
    assert result["statuses"]["NLSP_CUBIC_MODEL_AUDIT"] == "PASS"


def test_matching_cache_and_report_only_make_zero_derivation_and_root_calls(completed_audit, monkeypatch, capsys):
    result, bundle = completed_audit
    def forbidden(*args, **kwargs):
        raise AssertionError("Saved report must not regenerate a derivation or roots")
    monkeypatch.setattr(audit, "run_audit", forbidden)
    monkeypatch.setattr(rod, "derive_polynomials", forbidden)
    reused, _, metrics = audit.compute()
    assert reused["fingerprint"] == result["fingerprint"]
    assert metrics == {"cache_reused": True, "derivation_calls": 0, "root_evaluations": 0}
    assert audit.main(["--report-only", str(bundle)]) == 0
    assert '"derivation_calls": 0' in capsys.readouterr().out


def test_cache_identity_includes_supplied_input_and_all_frozen_references(monkeypatch):
    _, original, fingerprint = audit.identity()
    reference_keys = [name for name in original["hashes"] if "results" in name]
    assert reference_keys and any("mh_modes.csv" in name for name in reference_keys)
    sha = audit.sha
    monkeypatch.setattr(audit, "sha", lambda path: "0"*64 if str(path).endswith("mh_modes.csv") else sha(path))
    assert audit.identity()[2] != fingerprint


def test_corrupted_cache_is_rejected_without_mutating_original(completed_audit, tmp_path):
    import json
    import shutil
    _, bundle = completed_audit
    clone = tmp_path/"cache"/bundle.name
    shutil.copytree(bundle, clone)
    (clone/"result.json").write_text("{}", encoding="utf8")
    with pytest.raises(ValueError, match="artifact mismatch"):
        audit.compute(output=clone.parent)


def test_zero_poisson_ratio_does_not_remove_rotation_driven_contraction(coefficients):
    no_poisson = rod.RodCoefficients.rectangular(1.,1.,0.,.20,.05,coefficients.CT)
    velocities = np.zeros(7); velocities[5] = .1
    jet = rod.FieldJet(np.zeros(7), np.zeros(7), velocities, *(np.zeros(7) for _ in range(3)))
    assert rod.full_evaluate(jet,no_poisson)["residual"][6] == pytest.approx(-no_poisson.jp*.1**2)


def test_pure_twist_with_c_zero_is_not_a_false_invariant_subspace(coefficients):
    velocities = np.zeros(7); velocities[3] = .1
    jet = rod.FieldJet(np.zeros(7), np.zeros(7), velocities, *(np.zeros(7) for _ in range(3)))
    assert rod.full_evaluate(jet,coefficients)["residual"][6] == pytest.approx(-coefficients.jp*.1**2)


@pytest.mark.parametrize("change", ("field_order", "kappa"))
def test_audit_config_cannot_silently_change_physics_or_unknown_order(tmp_path, change):
    import json
    config=json.loads(audit.CONFIG.read_text(encoding="utf8"))
    if change=="field_order":
        config["field_order"][0],config["field_order"][1]=config["field_order"][1],config["field_order"][0]
    else:
        config["material_geometry"]["kappa"]="1"
    path=tmp_path/"bad_config.json"
    path.write_text(json.dumps(config),encoding="utf8")
    with pytest.raises(ValueError):
        audit.identity(path)
