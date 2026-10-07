"""General frame gates; fixed cases, exact maps, independent energy count."""
from fractions import Fraction as F
import json
import math

import numpy as np
import pytest

from scripts.analysis import verify_mindlin_herrmann_timoshenko_general_beta_joint as audit
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint
from scripts.lib import mindlin_herrmann_longitudinal as mh


@pytest.fixture(scope="module")
def setup():
    return audit.check_inputs()


@pytest.fixture(scope="module")
def computed(setup):
    config, model, length, _, direct, frozen = setup
    return audit.compute(config, model, length, direct, frozen)


@pytest.mark.parametrize("beta", (0., 5., 45., 90., -45., 13.7))
def test_project_geometry_orthonormal_orientation_and_inverse(beta):
    for frame in joint.frames(beta):
        np.testing.assert_allclose(frame.translation.T@frame.translation, np.eye(2), atol=2e-16)
        assert np.linalg.det(frame.translation) == pytest.approx(-1., abs=2e-16)
        np.testing.assert_allclose(frame.nodal_transform.T@frame.nodal_transform, np.eye(4), atol=2e-16)
        # c/theta invariant under proper rotations, no artificial trig factors.
        np.testing.assert_array_equal(frame.nodal_transform@[0., 2., 0., -3.], [0., 2., 0., -3.])


def test_exact_beta0_and_axis_aligned_right_angle():
    zero, right = joint.frames(0.), joint.frames(90.)
    np.testing.assert_array_equal(zero[0].translation, [[1., 0.], [0., -1.]])
    np.testing.assert_array_equal(zero[1].translation, [[-1., 0.], [0., 1.]])
    np.testing.assert_allclose(right[1].translation, [[0., -1.], [-1., 0.]], rtol=0, atol=7e-17)


def test_rational_duality_before_equilibrium():
    t, n = (F(3, 5), F(4, 5)), (F(4, 5), -F(3, 5))
    p, delta = (F(5, 7), F(11, 13)), (F(2, 3), F(7, 9))
    lhs = p[0]*sum(a*b for a, b in zip(t, delta))+p[1]*sum(a*b for a, b in zip(n, delta))
    rhs = sum((p[0]*t[i]+p[1]*n[i])*delta[i] for i in range(2))
    assert lhs == rhs


def test_deterministic_arbitrary_state_virtual_work(setup):
    r = audit.structural_checks(setup[1], setup[2])
    assert max(x["virtual_work_max_error"] for x in r["geometry"]) < audit.POLICY["coordinate_duality_tol"]


@pytest.mark.parametrize("beta", (0., 5., 45., 90.))
def test_eight_conditions_rank_eight_and_exact_gram(beta):
    operator = joint.joint_matrix(joint.frames(beta))
    assert operator.shape == (8, 16) and np.linalg.matrix_rank(operator) == 8
    np.testing.assert_allclose(operator@operator.T, 2*np.eye(8), atol=5e-16)


def test_only_one_boundary_assembly_including_frozen_wrapper(setup):
    _, model, length, *_ = setup
    for split in (.5, .35, .65):
        np.testing.assert_array_equal(joint.boundary_matrix(model, length, split, 3.),
            joint.frame_boundary_matrix(model, (length*split, length*(1-split)), 3., beta_deg=0.))


def test_frozen_beta0_all_inventories_frequencies_modes_and_c_R(computed):
    r = computed["beta0_regression"]
    assert r["status"] == "PASS" and len(r["rows"]) == 54
    assert r["matrix_max_absolute_error"] == 0.
    assert r["fresh_arm_swap"]["status"] == "PASS"
    for row in r["rows"]:
        assert row["frequency_relative_difference"] < audit.POLICY["beta0_frequency_relative_tol"]
        assert max(row["component_L2_relative"]) < audit.POLICY["beta0_component_L2_tol"]
        assert row["joint_residuals"]["c"] < 1e-9 and row["joint_residuals"]["R_node"] < 1e-9


def test_zero_limit_geometry_operator_and_first_three_roots(computed):
    assert [x["beta_deg"] for x in computed["zero_limit"]] == list(audit.SMALL_ANGLES)
    for row in computed["zero_limit"]:
        assert max(row["first_three_relative_differences"]) < row["allowance"]
        assert row["joint_operator_norm_difference"] <= math.radians(row["beta_deg"])*(1+1e-10)
    assert max(computed["zero_limit"][0]["first_three_relative_differences"]) < 1e-10


def test_arm_label_swap_operator_frequency_shape_and_balance(computed):
    assert computed["structure"]["arm_swap_operator_error"] == 0
    rows = computed["symmetry_gate"]["arm_swap"]["rows"]
    assert len(rows) >= 3
    assert computed["symmetry_gate"]["swapped"]["arm_frames"] == computed["symmetry_gate"]["canonical"]["arm_frames"][::-1]
    for row in rows:
        assert row["frequency_relative_difference"] < audit.POLICY["symmetry_frequency_relative_tol"]
        assert max(row["component_L2_relative"]) < audit.POLICY["symmetry_component_L2_tol"]


def test_reflection_pseudoscalar_theta_M_and_scalar_c_R(computed):
    frames = joint.frames(45.)
    for mirrored, canonical in zip(joint.reflected_frames(frames), joint.frames(-45.)):
        np.testing.assert_array_equal(mirrored.translation, canonical.translation)
    state = np.arange(1., 9.)
    mirror_global = np.array([1., 1., -1., -1.])  # dX,c,dY,theta
    for a, b in zip(frames, joint.reflected_frames(frames)):
        np.testing.assert_allclose(joint.nodal_efforts(state*joint.MIRROR_STATE, b, "right"),
                                   joint.nodal_efforts(state, a, "right")*mirror_global)
    assert joint.MIRROR_STATE[1] == joint.MIRROR_STATE[5] == 1.
    assert joint.MIRROR_STATE[3] == joint.MIRROR_STATE[7] == -1.
    assert computed["symmetry_gate"]["reflection"]["status"] == "PASS"
    for saved, frame in zip(computed["symmetry_gate"]["mirrored"]["arm_frames"], joint.frames(-45.)):
        np.testing.assert_array_equal(saved["local_to_global"], frame.nodal_transform)


def test_local_operators_unchanged_but_global_translations_mix(setup):
    model = setup[1]
    local = mh.full_harmonic_state_matrix(model, 3.)
    assert not np.any(local[np.ix_(joint.BLOCK_INDICES["mh"], joint.BLOCK_INDICES["timoshenko"])])
    assert joint.joint_matrix(joint.frames(45.))[0, 10] != 0.
    assert joint.joint_matrix(joint.frames(0.))[0, 10] == 0.
    assert model.section.K == model.mh_shear_factor == 5/6
    assert model.mh_inertia_factor == model.tim_rotary_factor == 1.


def test_energy_schur_symmetry_and_pole_certificates(setup, computed):
    _, model, length, *_ = setup
    for beta in (0., *audit.ANGLES):
        matrix, diagnostic = joint.nodal_schur_matrix(model, (length/2, length/2), 3., beta_deg=beta)
        np.testing.assert_array_equal(matrix, matrix.T)
        assert diagnostic["scaled_skew_residual"] < audit.POLICY["schur_symmetry_tol"]
    for block in computed["fixed_arm_poles"].values():
        assert len(block["roots"]) == block["count_bound"]["upper_count"]
        assert all(r["independent"]["relative_difference"] < 1e-10 for r in block["roots"])


def test_schur_inertia_count_matches_known_straight_roots(setup, computed):
    _, model, length, _, _, frozen = setup
    count = audit.count_function(model, length, computed["fixed_arm_poles"])
    known = [r["omega"] for block in ("mh", "timoshenko") for r in frozen["splits"][0][block]["roots"]]
    for frequency in (.04, .10, .20, .40, .60, .70):
        omega = frequency*2*math.pi
        actual, _ = count(omega, 0.)
        assert actual == sum(w < omega for w in known)


def test_count_guard_handles_dirichlet_pole_without_fake_global_root(setup, computed):
    _, model, length, *_ = setup
    count = audit.count_function(model, length, computed["fixed_arm_poles"])
    pole = computed["fixed_arm_poles"]["timoshenko"]["roots"][0]["omega"]
    value, diagnostic = count(pole, 45.)
    assert diagnostic["both_counts"] == [value, value]
    assert "pole_exclusion" in diagnostic


@pytest.mark.parametrize("beta", audit.ANGLES)
def test_pilot_count_certified_inventory_and_guard(computed, beta):
    case = next(c for c in computed["pilot"] if c["beta_deg"] == beta)
    assert len(case["roots"]) == case["search"]["upper_count"]
    assert case["search"]["lower_count"] == 0 and len(case["roots"]) >= 13
    assert not case["search"]["failed_intervals"]
    for root in case["roots"]:
        assert root["bracket_counts"][1]-root["bracket_counts"][0] == 1
        assert root["bracket_determinants"][0]*root["bracket_determinants"][1] <= 0
    frequencies = [r["frequency_hz"] for r in case["roots"]]
    assert all(a < b for a, b in zip(frequencies[:-1], frequencies[1:]))


@pytest.mark.parametrize("row", joint.JOINT_ROWS)
def test_each_joint_compatibility_force_R_moment_residual(computed, row):
    for case in computed["pilot"]:
        for root in case["roots"]:
            assert root["diagnostics"]["joint_residual_scaled"][row] < audit.POLICY["boundary_joint_scaled_tol"]


def test_clamps_equations_mass_energy_and_conditioning(computed):
    for case in computed["pilot"]:
        assert case["mass_gram_max_error"] < audit.POLICY["mass_gram_tol"]
        for root in case["roots"]:
            d = root["diagnostics"]
            assert d["mass_norm"] == pytest.approx(1., abs=1e-12)
            assert d["clamp_scaled_residual"] < 1e-9 and d["equation_scaled_residual"] < 1e-9
            assert d["energy_relative_error"] < 5e-8
            assert d["nonzero_singular_condition"] < 1e8


def test_structural_failure_prevents_any_nonzero_spectrum(monkeypatch, setup):
    def fail(*args):
        raise ArithmeticError("duality gate failed")
    monkeypatch.setattr(audit, "structural_checks", fail)
    monkeypatch.setattr(audit, "solve_case", lambda *args, **kwargs: pytest.fail("Forbidden spectrum"))
    with pytest.raises(ArithmeticError, match="duality"):
        audit.compute(setup[0], setup[1], setup[2], setup[4], setup[5])


def test_beta0_failure_prevents_nonzero_angles(monkeypatch, setup):
    def fail(*args):
        raise ArithmeticError("beta0 regression failed")
    monkeypatch.setattr(audit, "beta0_regression", fail)
    monkeypatch.setattr(audit, "solve_case", lambda *args, **kwargs: pytest.fail("Forbidden spectrum"))
    with pytest.raises(ArithmeticError, match="beta0"):
        audit.compute(setup[0], setup[1], setup[2], setup[4], setup[5])


def test_cache_reuse_and_corruption(tmp_path, monkeypatch, computed):
    monkeypatch.setattr(audit, "compute", lambda *args: computed)
    assert audit.main(["--compute", "--output-dir", str(tmp_path)]) == 0
    monkeypatch.setattr(audit, "compute", lambda *args: pytest.fail("Matching cache must not calculate roots"))
    assert audit.main(["--compute", "--output-dir", str(tmp_path)]) == 0
    pointer = json.loads((tmp_path/"current.json").read_text(encoding="utf-8"))
    (tmp_path/pointer["fingerprint"]/"result.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="stale"):
        audit.main(["--compute", "--output-dir", str(tmp_path)])
