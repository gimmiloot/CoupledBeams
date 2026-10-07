"""Reduced joint closure/transparency; beta0 only, no source-fit goldens."""
from fractions import Fraction as F
import json

import numpy as np
import pytest

from scripts.analysis import verify_mindlin_herrmann_timoshenko_beta0_joint as audit
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint
from scripts.lib import mindlin_herrmann_longitudinal as mh


@pytest.fixture(scope="module")
def setup():
    return audit.check_reference()


@pytest.fixture(scope="module")
def computed(setup):
    config, model, length, _, _, _, ref = setup
    return audit.compute(config, model, length, ref)[0]


def test_eight_independent_invariant_joint_conditions():
    matrix = joint.joint_matrix()
    assert matrix.shape == (8, 16)
    assert np.linalg.matrix_rank(matrix) == len(joint.JOINT_ROWS) == 8


def test_translation_transform_orthogonal_and_dual():
    for frame in (*joint.beta0_frames(), joint.Frame((.6, .8), (.8, -.6))):
        np.testing.assert_allclose(frame.translation.T@frame.translation, np.eye(2), atol=1e-16)
        np.testing.assert_allclose(frame.nodal_transform.T@frame.nodal_transform, np.eye(4), atol=1e-16)
        assert frame.nodal_transform[1, 1] == frame.nodal_transform[3, 3] == 1
        assert np.count_nonzero(frame.nodal_transform[1]) == np.count_nonzero(frame.nodal_transform[3]) == 1
    with pytest.raises(ValueError, match="orthonormal"):
        joint.Frame((1., 1.), (0., -1.))


def test_scalar_c_and_theta_under_proper_planar_rotation():
    q = np.array([0., 2., 0., -3.])
    for frame in (*joint.beta0_frames(), joint.Frame((.8, -.6), (-.6, -.8))):
        np.testing.assert_array_equal(frame.nodal_transform@q, q)


def test_boundary_signs_from_exact_integration_by_parts():
    # p=x^2+2x+3, deltaq=x^3+x: integral (p*dq' + p'*dq)=[p*dq].
    p, q = [F(3), F(2), F(1)], [F(0), F(1), F(0), F(1)]
    derivative = lambda a: [i*a[i] for i in range(1, len(a))]
    product = lambda a, b: [sum((a[i]*b[n-i] for i in range(len(a)) if 0 <= n-i < len(b)), F(0)) for n in range(len(a)+len(b)-1)]
    a, b = product(p, derivative(q)), product(derivative(p), q)
    integral = sum((v/F(i+1) for i, v in enumerate(a)), F(0))+sum((v/F(i+1) for i, v in enumerate(b)), F(0))
    assert integral == sum(p)*sum(q)-p[0]*q[0]
    assert joint.endpoint_sign("left") == -1 and joint.endpoint_sign("right") == 1
    with pytest.raises(ValueError):
        joint.endpoint_sign("joint")


def test_exact_rational_duality_and_arbitrary_nodal_virtual_work():
    t, n, force, variation = (F(3, 5), F(4, 5)), (F(4, 5), -F(3, 5)), (F(7, 3), F(5, 7)), (F(11, 13), F(3, 2))
    local = sum(force[i]*sum(a[j]*variation[j] for j in range(2)) for i, a in enumerate((t, n)))
    global_value = sum(sum(force[i]*a[j] for i, a in enumerate((t, n)))*variation[j] for j in range(2))
    assert local == global_value
    assert audit.virtual_work_check()["max_scaled_error"] < audit.POLICY["virtual_work_relative_tol"]


def test_c_R_transparency_signs_are_not_naive_sum_or_difference():
    state = np.arange(1., 9.)
    right = state*joint.REFLECTION
    np.testing.assert_array_equal(joint.joint_residual(state, right), np.zeros(8))
    assert right[1] == state[1] and right[5] == -state[5]
    assert right[3] == state[3] and right[7] == -state[7]


def test_joint_force_transform_dual_for_non_equilibrium_state():
    state, variation = np.arange(1., 9.), np.array([3., -2., 4., 1.])
    for frame in joint.beta0_frames():
        for end in ("left", "right"):
            assert joint.nodal_efforts(state, frame, end)@variation == pytest.approx(
                joint.endpoint_sign(end)*(state[4:]@(frame.nodal_transform.T@variation)), abs=1e-14)


@pytest.mark.parametrize("split", audit.SPLITS)
def test_primary_family_frequencies_match_direct_and_count_guard(computed, split):
    case = next(c for c in computed["splits"] if c["split"] == split)
    for block, count in (("mh", 7), ("timoshenko", 11)):
        assert len(case[block]["roots"]) == case[block]["completeness"]["upper_count"] == count
        assert case[block]["below_search_count_bound"]["upper_count"] == 0
        assert len(case[block]["search_attempts"]) == 1
        assert not case[block]["search_attempts"][0]["failed_intervals"]
        for root in case[block]["roots"]:
            assert root["relative_difference"] < audit.POLICY["frequency_relative_tol"]
            assert root["bracket_determinants"][0]*root["bracket_determinants"][1] <= 0
            assert root["bracket_omega"][0] <= root["omega"] <= root["bracket_omega"][1]
    assert len(case["combined"]) == 13
    assert case["combined"][-1]["guard"]
    assert not any(r["guard"] for r in case["combined"][:-1])


def test_split_location_invariance(computed):
    spectra = [[r["frequency_hz"] for r in case["combined"]] for case in computed["splits"]]
    for values in spectra[1:]:
        np.testing.assert_allclose(values, spectra[0], rtol=audit.POLICY["frequency_relative_tol"], atol=0)


def test_arm_swap_reflection_with_all_fields(computed):
    assert computed["arm_swap"]["status"] == "PASS"
    assert len(computed["arm_swap"]["rows"]) == 18
    for row in computed["arm_swap"]["rows"]:
        assert 1-row["mass_MAC"] < audit.POLICY["MAC_loss_tol"]
        assert max(row["component_L2_relative"]) < audit.POLICY["component_L2_relative_tol"]


@pytest.mark.parametrize("block", ["mh", "timoshenko"])
def test_direct_coupled_mass_overlap_including_c(computed, block):
    for case in computed["splits"]:
        for root in case[block]["roots"]:
            comparison = root["mode_comparison"]
            assert 1-comparison["mass_MAC"] < audit.POLICY["MAC_loss_tol"]
            assert max(comparison["component_L2_relative"].values()) < audit.POLICY["component_L2_relative_tol"]
            assert comparison["mass_norm"] == pytest.approx(1., abs=1e-12)
            if block == "mh":
                assert "c" in comparison["component_L2_relative"]


@pytest.mark.parametrize("row", joint.JOINT_ROWS)
def test_every_interface_compatibility_and_equilibrium_row(computed, row):
    for case in computed["splits"]:
        for block in ("mh", "timoshenko"):
            for root in case[block]["roots"]:
                assert root["mode_comparison"]["joint_residual_scaled"][row] < audit.POLICY["interface_scaled_residual_tol"]
                assert root["mode_comparison"]["clamp_scaled_residual"] < audit.POLICY["clamp_scaled_residual_tol"]


def test_full_matrix_permutation_has_exact_zero_mixed_blocks(setup):
    _, model, length, *_ = setup
    rows = (0, 1, 4, 5, 8, 10, 12, 14, 2, 3, 6, 7, 9, 11, 13, 15)
    cols = (*range(4), *range(8, 12), *range(4, 8), *range(12, 16))
    for split in audit.SPLITS:
        full = joint.boundary_matrix(model, length, split, 3.)
        ordered = full[np.ix_(rows, cols)]
        assert not np.any(ordered[:8, 8:]) and not np.any(ordered[8:, :8])
        np.testing.assert_array_equal(ordered[:8, :8], joint.block_matrix(model, length, split, 3., "mh"))
        np.testing.assert_array_equal(ordered[8:, 8:], joint.block_matrix(model, length, split, 3., "timoshenko"))


def test_full_matrix_zeros_equal_union_of_verified_blocks(setup, computed):
    model, length = setup[1:3]
    for case in computed["splits"]:
        for block in ("mh", "timoshenko"):
            for root in case[block]["roots"]:
                matrix = joint.boundary_matrix(model, length, case["split"], root["omega"])
                matrix /= np.linalg.norm(matrix, axis=1)[:, None]
                singular = np.linalg.svd(matrix, compute_uv=False)
                assert singular[-1]/singular[0] < 1e-9


def test_segmented_state_semigroup_is_independent_of_joint_basis(computed):
    for case in computed["splits"]:
        for block in ("mh", "timoshenko"):
            for root in case[block]["roots"]:
                assert root["segmented_transfer"]["direct_projected_boundary_difference"] < audit.POLICY["segmented_transfer_tol"]


def test_no_optical_cutoff_mislabeled_as_finite_mode(computed):
    for case in computed["splits"]:
        assert case["mh"]["roots"][-1]["frequency_hz"] < computed["limits"]["contraction_cutoff_hz"]
        assert case["timoshenko"]["roots"][-1]["frequency_hz"] < computed["limits"]["shear_cutoff_hz"]
        assert {r["family"] for r in case["combined"]} == {"bending", "axial_acoustic"}


def test_nonzero_angle_spectrum_is_rejected(setup, computed):
    with pytest.raises(ValueError, match="Only beta=0"):
        joint.boundary_matrix(setup[1], setup[2], .5, 2., beta_deg=1.)
    assert all(case["beta_deg"] == 0 for case in computed["splits"])
    assert computed["coordinate_contract"]["beta_deg"] == 0


def test_arm_physics_unchanged_and_no_reverse_import_dependency(setup):
    model = setup[1]
    assert model.variant == mh.PROJECT_VARIANT
    assert model.mh_shear_factor == model.section.K == 5/6
    assert model.mh_inertia_factor == model.tim_rotary_factor == 1.
    for path in (audit.single.CONFIG, audit.ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py",
                 audit.ROOT/"scripts/analysis/reproduce_mindlin_herrmann_timoshenko_literature.py"):
        assert "mindlin_herrmann_timoshenko_joint" not in path.read_text(encoding="utf-8")
    reference_manifest = setup[5]
    for path, digest in reference_manifest["identity"]["files"].items():
        assert audit.sha(audit.ROOT/path) == digest


def test_reference_reuse_and_cache_corruption_gate(tmp_path, monkeypatch, computed):
    monkeypatch.setattr(audit, "compute", lambda *args: (computed, [{"x": 0., "c": 0.}]))
    assert audit.main(["--compute", "--output-dir", str(tmp_path)]) == 0
    monkeypatch.setattr(audit, "compute", lambda *args: pytest.fail("Matching cache must not recompute roots"))
    assert audit.main(["--compute", "--output-dir", str(tmp_path)]) == 0
    pointer = json.loads((tmp_path/"current.json").read_text(encoding="utf-8"))
    (tmp_path/pointer["fingerprint"]/"result.json").write_text("{}", encoding="utf-8")
    with pytest.raises(ValueError, match="stale"):
        audit.main(["--compute", "--output-dir", str(tmp_path)])


def test_hard_variational_gate_precedes_any_spectrum(monkeypatch, setup):
    def failed():
        raise ArithmeticError("Virtual-work gate failed")
    monkeypatch.setattr(audit, "virtual_work_check", failed)
    monkeypatch.setattr(audit, "solve_split", lambda *args: pytest.fail("Hard gate must stop root search"))
    with pytest.raises(ArithmeticError, match="Virtual-work"):
        audit.compute(setup[0], setup[1], setup[2], setup[-1])


def test_cli_gate_precedes_optional_direct_reference_generation(monkeypatch):
    def failed():
        raise ArithmeticError("Virtual-work gate failed")
    monkeypatch.setattr(audit, "virtual_work_check", failed)
    monkeypatch.setattr(audit, "check_reference", lambda **kwargs: pytest.fail("No reference generation before gate"))
    with pytest.raises(ArithmeticError, match="Virtual-work"):
        audit.main(["--compute"])
