"""Postprocessing contracts; no FEM, ODE, equilibrium or eigen solves."""
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.polynomial.legendre import legval

from scripts.lib import nlsp_profile_1d_diagnostics as profile


def coefficients():
    return SimpleNamespace(C=2., H=.02, S=.75, nu=.3)


def test_canonical_planar_definitions_and_axial_resultant():
    fields = np.array([[.1, .2, .05, -.03]])
    gradients = np.array([[.02, .04, .01, .001]])
    result = profile.planar_measures(fields, gradients, coefficients())
    expected = 1.02 * np.cos(.05) + .04 * np.sin(.05) - 1
    assert result["Gamma1"][0] == pytest.approx(expected)
    assert result["N"][0] == pytest.approx(2 * (expected - .009))
    assert result["F_axial_global"][0] == pytest.approx(result["N"][0] * np.cos(.05)
        - result["Q"][0] * np.sin(.05))
    assert result["Gamma2"][0] == pytest.approx(-1.02 * np.sin(.05) + .04 * np.cos(.05))
    assert result["minus_nu_Gamma1"][0] == pytest.approx(-.3 * expected)
    assert result["theta_minus_w_s"][0] == pytest.approx(.01)


def test_rigid_rotation_has_no_full_planar_strain():
    theta = .2
    f = np.array([[0., 0., theta, 0.]])
    d = np.array([[np.cos(theta) - 1, np.sin(theta), 0., 0.]])
    result = profile.planar_measures(f, d, coefficients())
    assert result["Gamma1"][0] == pytest.approx(0., abs=3e-16)
    assert result["Gamma2"][0] == pytest.approx(0., abs=3e-16)
    assert result["N"][0] == pytest.approx(0., abs=6e-16)
    assert result["F_axial_global"][0] == pytest.approx(0., abs=1e-15)


def test_linear_static_reference_has_zero_first_order_axial_resultant():
    from scripts.analysis.verify_nlsp_nonlinear_static_3d_fem import fem2_tim_uniform
    x = np.linspace(0., 1., 401)
    f, d = fem2_tim_uniform(x, .002, 1., .3, .75)
    # At first order u_s=c=0; geometric terms remain a second-order diagnostic.
    first_order_N = coefficients().C * (d[:, 0] + coefficients().nu * f[:, 3])
    np.testing.assert_array_equal(first_order_N, np.zeros(len(x)))
    assert np.max(abs(d[:, 1] - f[:, 2])) > 0
    assert profile.symmetry_metrics(x, f)["theta"]["max_absolute_error"] < 1e-17


def test_retained_cubic_measures_are_not_silently_replaced():
    f = np.array([[0., 0., .15, -.002]])
    d = np.array([[.04, .2, 0., 0.]])
    r = profile.planar_measures(f, d, coefficients())
    assert r["Gamma1_retained"][0] == pytest.approx(.04 + .15*.2 - .15**2/2 - .04*.15**2/2)
    assert r["Gamma2_retained"][0] == pytest.approx(.2 - .15 - .15*.04 - .2*.15**2/2 + .15**3/6)
    assert abs(r["Gamma1"][0] - r["Gamma1_retained"][0]) > 1e-6


def test_uniform_axial_stretch_and_independent_contraction():
    f = np.zeros((3, 4)); d = np.zeros_like(f)
    d[:, 0] = .01
    f[:, 3] = -.003
    r = profile.planar_measures(f, d, coefficients())
    np.testing.assert_allclose(r["Gamma1"], .01, atol=3e-16)
    np.testing.assert_allclose(r["c_plus_nu_Gamma1"], 0., atol=3e-16)
    np.testing.assert_allclose(r["N"], 2*(.01-.0009))


@pytest.mark.parametrize("actual,requested", [([0., 1.], [1.0001]), ([0., 1.], [-.1]), ([0., 1.], [.5])])
def test_unsaved_times_are_not_available(actual, requested):
    with pytest.raises(ValueError, match="NOT_AVAILABLE"):
        profile.exact_saved_indices(actual, requested)


def test_exact_saved_times_are_preserved_without_interpolation():
    np.testing.assert_array_equal(profile.exact_saved_indices([0., .25, .5, .75, 1.], [0., .5, 1.]), [0, 2, 4])


@pytest.mark.parametrize("times", [[0., .5, .5], [0., .5, np.nan], [1., 0.]])
def test_invalid_saved_time_axes_rejected(times):
    with pytest.raises(ValueError):
        profile.exact_saved_indices(times, [0.])


def test_shen_conversion_preserves_exact_polynomial_and_endpoints():
    raw = np.array([.2, -.3, .01, .05, -.08])
    xi = np.linspace(-1, 1, 51)
    independent = sum(a * (np.polynomial.legendre.Legendre.basis(i)(xi)
        - np.polynomial.legendre.Legendre.basis(i+2)(xi)) for i, a in enumerate(raw))
    result = profile.shen_to_legendre(raw)
    np.testing.assert_allclose(legval(xi, result), independent, atol=2e-16)
    np.testing.assert_allclose(legval(np.array([-1., 1.]), result), 0., atol=2e-16)
    np.testing.assert_array_equal(raw, [.2, -.3, .01, .05, -.08])


def test_legendre_l2_has_physical_length_weight():
    assert profile.polynomial_l2(np.array([3.]), 2.) == pytest.approx(3*np.sqrt(2))
    assert profile.polynomial_l2(np.array([0., 3.]), 2.) == pytest.approx(np.sqrt(6))


def test_tail_is_explanatory_and_never_changes_input():
    co = np.array([1., .1, .01, .001, .002, .003, .004, .005, .006])
    original = co.copy(); x = np.linspace(0, 1, 101)
    result = profile.legendre_tail_diagnostic(co, 1., x)
    assert result["degree_cutoff"] == 6
    assert result["tail_L2_fraction"] < 1
    assert result["tail_derivative_L2_fraction"] > result["tail_L2_fraction"]
    np.testing.assert_array_equal(co, original)


def test_classical_polynomial_benchmark_is_clamped_and_antisymmetric():
    # A symmetric quartic fixed-fixed bending shape, represented in Legendre.
    w = np.polynomial.Polynomial([1., 0., -2., 0., 1.]).convert(kind=np.polynomial.Legendre)
    x = np.linspace(0, 2., 201)
    result = profile.classical_axial_benchmark(.01*w.coef, 2., x)
    assert result["mean_bending_extension"] > 0
    assert result["u_endpoint_error"] < 1e-18
    np.testing.assert_allclose(result["u"], -result["u"][::-1], atol=2e-19)
    np.testing.assert_allclose(result["u_s"], result["u_s"][::-1], atol=2e-18)


def test_symmetry_and_four_independent_essential_clamps():
    x = np.linspace(0, 1, 501)
    even = x * (1-x)
    odd = even * (x-.5)
    f = np.column_stack((odd, even, 2*odd, -.1*even))
    metrics = profile.symmetry_metrics(x, f)
    assert all(v["max_absolute_error"] < 1e-16 for v in metrics.values())
    assert all(v["endpoint_max"] == 0 for v in metrics.values())
    assert metrics["u"]["parity"] == metrics["theta"]["parity"] == "odd"
    assert metrics["w"]["parity"] == metrics["c"]["parity"] == "even"


def test_asymmetry_remains_visible():
    x = np.linspace(0, 1, 101)
    f = np.column_stack((x*(1-x), x*(1-x), x*(1-x), -x*(1-x)))
    metrics = profile.symmetry_metrics(x, f)
    assert metrics["theta"]["max_absolute_error"] == pytest.approx(.5)
    assert metrics["u"]["max_absolute_error"] == pytest.approx(.5)


def test_zero_crossing_coordinates_are_physical_and_signed():
    assert profile.zero_crossings([0., .3, .5, 1.], [-1., 1., 0., -1.]) == pytest.approx([.15, .5])


def test_region_partition_keeps_full_and_boundary_results():
    x = np.linspace(0, 1., 1001)
    value = np.ones_like(x); value[x < .05] = 10.
    r = profile._regional_measure(value, x, profile._grid_weights(x), .02)
    assert r["full"]["max_abs"] == r["boundary"]["max_abs"] == 10
    assert r["interior"]["max_abs"] == 1
    assert r["interior"]["sample_coordinate_range"][0] == pytest.approx(.06)


def test_helper_source_has_no_execution_entrypoint_or_calls():
    import inspect
    source = inspect.getsource(profile)
    for forbidden in ("integrate_case(", "solve_ivp(", "fem2_static_newton(",
                      "linear_eigenpairs(", "derive_polynomials(", "subprocess."):
        assert forbidden not in source


def test_audit_cannot_write_inside_immutable_source(tmp_path):
    source = tmp_path / "historical_bundle"
    source.mkdir()
    with pytest.raises(ValueError, match="must not overwrite"):
        profile.audit_saved_one_d(source, source / "new_child")
    assert list(source.iterdir()) == []


def test_actual_fem_time_extension_reads_each_unique_time_once(tmp_path, monkeypatch):
    import json
    from scripts.lib import nlsp_fem3b_continuation as previous

    source, output = tmp_path / "source", tmp_path / "audit"
    source.mkdir(); output.mkdir()
    (source / "dense_p64.npz").write_bytes(b"mocked accepted dense-polynomial fixture")
    q0 = np.array([.1, .2, .03, -.001])
    np.savez(source / "one_d_p64_nonlinear.npz", times=np.array([0., 1.]),
        q=np.array([q0, q0]), velocity=np.zeros((2, 4)))
    hashes = {name: profile._sha(source / name) for name in ("dense_p64.npz", "one_d_p64_nonlinear.npz")}
    (source / "manifest.json").write_text(json.dumps({"artifact_hashes": hashes}))
    (output / "one_d_profile_audit.json").write_bytes(b"immutable previous five-time evidence")
    old_bytes = (output / "one_d_profile_audit.json").read_bytes()

    class SavedDisc:
        ndof, nq, length = 4, 129, 1.
        coefficients = coefficients()

        def counters(self):
            return {name: 0 for name in ("rhs_calls", "jacobian_calls", "mass_factorizations",
                "force_evaluations", "linear_eigendecompositions")}

        def reconstruct_series(self, coordinates, x, derivative=0):
            if derivative:
                return np.zeros((len(coordinates), len(x), 4))
            return np.broadcast_to(coordinates[:, None, :], (len(coordinates), len(x), 4)).copy()

    monkeypatch.setattr(profile, "_load_discretizations", lambda source: ({64: SavedDisc()}, {"T1": 1.}))
    evaluations = []

    def evaluate_saved_dense(path, times):
        evaluations.append(np.asarray(times).copy())
        # A deliberately different zero-time interpolation result must not
        # replace the exactly saved physical initial state.
        return np.column_stack([np.tile(99. + times[:, None], (1, 4)), np.zeros((len(times), 4))])

    monkeypatch.setattr(previous, "evaluate_dense_records", evaluate_saved_dense)
    selected = [dict(name=name, actual_time=t, requested_time=t + offset, requested_tau=t + offset)
        for name, t, offset in (("static_a", 0., 0.), ("medium", .25, .01),
                              ("fine", .25, .01), ("later", .5, 0.))]
    result = profile.audit_actual_one_d_states(source, output, selected, spatial_points=21)
    assert result["new_scientific_solver_calls"] == 0
    assert result["states"] == 4 and result["distinct_physical_times"] == 3
    assert len(evaluations) == 1
    np.testing.assert_array_equal(evaluations[0], [0., .25, .5])
    with np.load(output / "actual_time_one_d.npz") as saved:
        np.testing.assert_array_equal(saved["q"][0], q0)
        np.testing.assert_array_equal(saved["q"][1], saved["q"][2])
        np.testing.assert_array_equal(saved["times"], [0., .25, .25, .5])
        np.testing.assert_array_equal(saved["requested_times"], [0., .26, .26, .5])
    assert result["source_hashes_unchanged_after_processing"]
    assert (output / "one_d_profile_audit.json").read_bytes() == old_bytes


def test_actual_fem_time_extension_rejects_corrupt_saved_source(tmp_path):
    import json
    source = tmp_path / "source"
    source.mkdir()
    (source / "dense_p64.npz").write_bytes(b"changed source")
    (source / "manifest.json").write_text(json.dumps({"artifact_hashes": {"dense_p64.npz": "wrong"}}))
    with pytest.raises(ValueError, match="absent/corrupted"):
        profile.audit_actual_one_d_states(source, tmp_path / "output",
            [dict(name="zero", actual_time=0., requested_time=0., requested_tau=0.)])
    assert not (tmp_path / "output").exists()
