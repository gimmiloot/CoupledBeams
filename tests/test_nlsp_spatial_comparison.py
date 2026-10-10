"""Saved-array/synthetic comparisons; no real FEM, Newton, eigen, or IVP jobs."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import nlsp_spatial_comparison as comparison
from scripts.lib import weakly_nonlinear_spatial_rod as rod


class SyntheticDisc:
    def __init__(self, p):
        self.p, self.ndof, self.length = p, 7, 1.
        self.coefficients = SimpleNamespace(C=1., jp=1.)

    def reconstruct_series(self, rows, x, derivative=0):
        x = np.asarray(x)
        shape = np.sin(np.pi*x) if derivative == 0 else np.pi*np.cos(np.pi*x)
        return np.asarray(rows)[:, None, :]*shape[None, :, None]


def save_case(path, times, nonlinear, linear):
    path.mkdir()
    np.savez_compressed(path/"trajectory.npz", times=times, q=nonlinear, velocity=nonlinear*.03)
    np.savez_compressed(path/"linear_trajectory.npz", times=times, q=linear, velocity=linear*.03)
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in path.iterdir() if p.is_file()}
    (path/"case.json").write_text(json.dumps({"status": "PASS", "artifact_hashes": files}), encoding="utf8")


@pytest.fixture
def stage_b(tmp_path):
    times = np.linspace(0., comparison.DEFAULT_T1/4, 5)
    linear = np.zeros((len(times), 7)); linear[:, 1] = .004*np.cos(.7*times); linear[:, 2] = .001*np.cos(1.1*times)
    delta = np.zeros_like(linear); delta[:, 1] = -1e-6*(1+times)
    delta[:, 2] = -1.5e-6*(1+times); delta[:, 3] = 2e-5*(1+times/2)
    delta[:, 0] = 3e-5; delta[:, 6] = -2e-5
    joint64 = linear+delta; joint48 = joint64.copy(); joint48[:, 1] += 1e-10*(1+times)
    isolated_w, isolated_v = np.zeros_like(linear), np.zeros_like(linear)
    isolated_w[:, 1] = linear[:, 1]+delta[:, 1]/2
    isolated_v[:, 2] = linear[:, 2]+delta[:, 2]/2
    isolated_w[:, 0] = isolated_v[:, 0] = delta[:, 0]/3
    arrays = {"joint_p64": (joint64, linear), "joint_p48": (joint48, linear),
              "isolated_w_p64": (isolated_w, isolated_w), "isolated_v_p64": (isolated_v, isolated_v)}
    cases = {}
    for name, (nl, li) in arrays.items():
        cases[name] = tmp_path/name; save_case(cases[name], times, nl, li)
    return cases, {48: SyntheticDisc(48), 64: SyntheticDisc(64)}, tmp_path/"comparison"


def test_full_horizon_denominators_and_signed_metrics():
    x, times = np.linspace(0., 1., 101), np.array((0., 1.))
    second = np.stack((np.sin(np.pi*x), np.sin(np.pi*x)*.1))
    first = second-.01*np.sin(np.pi*x)
    row = comparison.field_metrics(first, second, x, times, fixed_scale=.1, tolerance=.1)
    assert row["relative_max"] == pytest.approx(.01)
    assert row["signed_difference_at_absolute_max"] < 0
    assert row["reference_max_abs"] == 1. and row["status"] == "PASS"
    assert row["relative_max_fixed"] == pytest.approx(.1)


def test_planning_is_not_model_accuracy_threshold():
    row = comparison.planning_signal(8e-7, 1e-9)
    assert row["resolved_for_pre_FEM_planning"]
    assert row["not_physical_accuracy_gate_or_strict_error_bound"]
    row = comparison.planning_signal(2e-8, 1e-8)
    assert not row["resolved_for_pre_FEM_planning"]
    with pytest.raises(ValueError):
        comparison.planning_signal(-1., 0.)


def test_rotation_mapping_and_finite_orientation_composition():
    fields = np.zeros((2, 7)); fields[0, 3:6] = (.01, .02, .03)
    rotations = comparison.physical_rotations(fields)
    expected = rod.rotation_and_right_jacobian(np.array((.01, -.02, .03)))[0]
    np.testing.assert_allclose(rotations[0], expected, rtol=1e-15, atol=1e-15)
    np.testing.assert_allclose(rotations@rotations.transpose(0, 2, 1), np.broadcast_to(np.eye(3), rotations.shape), atol=3e-15)
    assert comparison.orientation_difference(rotations, rotations).max() < 1e-15
    first = np.zeros((1, 7)); second = first.copy(); second[0, 3] = 1e-10
    assert comparison.orientation_difference(comparison.physical_rotations(first), comparison.physical_rotations(second))[0] == pytest.approx(1e-10, rel=1e-12)


def test_retained_curvature_contains_mixed_rotation_source():
    fields = np.zeros((1, 7)); gradient = fields.copy()
    fields[0, 4:6] = (.02, .03); gradient[0, 4:6] = (.04, .02)
    result = comparison.retained_curvatures(fields, gradient)
    assert result[0, 0] != 0. and gradient[0, 3] == 0.
    assert result[0, 0] == pytest.approx(-.5*((-.02)*.02-.03*(-.04)))


def test_axis_nonplanarity_and_plane_change_are_distinct_from_coupling():
    x = np.linspace(0., 1., 101)
    fields = np.zeros((2, len(x), 7)); fields[:, :, 1] = np.sin(np.pi*x)*.01
    fields[:, :, 2] = fields[:, :, 1]*.4
    first = comparison.axis_plane_diagnostic(fields, x)
    assert first["RMS_distance"].max() < 1e-15
    fields[1, :, 2] = np.sin(2*np.pi*x)*.004
    second = comparison.axis_plane_diagnostic(fields, x)
    assert second["RMS_distance"][1] > 1e-4
    assert second["nonplanarity_is_not_by_itself_nonlinear_coupling"]


def test_stage_b_saved_controls_all14_and_no_fitting(stage_b, monkeypatch):
    cases, discs, output = stage_b
    monkeypatch.setattr(rod, "derive_polynomials", lambda: pytest.fail("No new symbolic derivation"))
    summary = comparison.analyze_stage_b(cases, discs, output, chunk_size=2)
    assert summary["scientific_calls"] == 0
    assert sum(len(summary["all14_spatial"][part]["fields"]) for part in ("q", "velocity")) == 14
    assert summary["all14_spatial"]["q"]["fields"]["c"]["tolerance"] == 1e-3
    assert summary["all14_spatial"]["q"]["fields"]["Phi"]["tolerance"] == 1e-4
    assert summary["isolated_p64_controls_have_no_independent_p_control"]
    assert summary["historical_p64_reference_gate_denominators_preserved"]
    assert summary["responses"]["mix_evolution"]["v"]["absolute_max"] > 0.
    with np.load(output/"stage_b_comparison.npz", allow_pickle=False) as z:
        assert z["evolution_midspan"][0].max() == 0.
        assert z["physical_rotation_snapshots"].shape[-2:] == (3, 3)
        assert z["spatial_q_max"].shape[-1] == 7


def test_source_hash_mismatch_and_timestamp_guard(stage_b):
    cases, discs, output = stage_b
    with (cases["joint_p48"]/"trajectory.npz").open("ab") as stream:
        stream.write(b"changed")
    with pytest.raises(ValueError, match="artifact changed"):
        comparison.analyze_stage_b(cases, discs, output)


def test_interpolation_never_extrapolates_and_methods_are_distinct():
    source = np.array((0., .2, .6, 1.)); values = source**2
    target = np.linspace(0., 1., 201)
    linear = comparison._interpolate(values, source, target, "linear")
    cubic = comparison._interpolate(values, source, target, "pchip")
    assert np.max(abs(linear-cubic)) > 1e-3
    with pytest.raises(ValueError, match="coverage"):
        comparison._interpolate(values, source, np.array((0., 1.1)), "linear")
    with pytest.raises(ValueError, match="predeclared"):
        comparison._interpolate(values, source, target, "phase_fit")


def test_full_matrix_space_transfer_preserves_so3_and_known_face_values():
    fields = np.zeros((2, 3, 7)); fields[:, :, 3] = .003
    history = {"time": np.array((0., 1.)), "x": np.linspace(0., 1., 5),
        "raw_rotation_x": np.array((.2, .5, .8)), "raw_rotation_matrices": comparison.physical_rotations(fields)}
    result = comparison._raw_rotations_on_sections(history)
    np.testing.assert_allclose(result@result.transpose(0, 1, 3, 2), np.broadcast_to(np.eye(3), result.shape), atol=2e-15)
    np.testing.assert_allclose(result[:, (0, -1)], np.broadcast_to(np.eye(3), result[:, (0, -1)].shape), atol=1e-15)


def save_dense_synthetic_case(path, disc):
    path.mkdir()
    H = comparison.DEFAULT_T1/4
    q0 = np.zeros(disc.ndof); q0[1:3] = (.004, .001)
    nonlinear0 = q0.copy(); nonlinear0[1] -= 1e-6
    initial = np.r_[nonlinear0, np.zeros(disc.ndof)]
    Q = np.zeros((1, 2*disc.ndof, 3)); Q[0, 1, 0] = -2e-6
    np.savez_compressed(path/"accepted_dense.npz", t_old=np.array((0.,)), t=np.array((H,)),
                        y_old=initial[None, :], Q=Q)
    np.savez_compressed(path/"static_states.npz", q_linear=q0, q_nonlinear=nonlinear0)
    np.savez_compressed(path/"linear_modes.npz", omega=np.ones(disc.ndof), vectors=np.eye(disc.ndof), M0=np.eye(disc.ndof))
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in path.iterdir()}
    (path/"case.json").write_text(json.dumps({"status": "PASS", "ndof": disc.ndof,
        "fields": list(comparison.FIELDS), "artifact_hashes": files}), encoding="utf8")
    return path


def test_saved_linear_and_dense_reader_without_eigensolve(tmp_path, monkeypatch):
    disc = SyntheticDisc(64)
    case = save_dense_synthetic_case(tmp_path/"case", disc)
    monkeypatch.setattr(np.linalg, "eigh", lambda *a, **kw: pytest.fail("No eigen solve in saved-state postprocessing"))
    times = np.array((0., comparison.DEFAULT_T1/8, comparison.DEFAULT_T1/4))
    result = comparison.evaluate_saved_one_d(case, disc, times)
    assert result["scientific_calls"] == 0
    assert result["q_linear"][0, 1] == .004
    assert result["q_linear"][-1, 1] == pytest.approx(.004*np.cos(times[-1]))
    assert result["q_nonlinear"][-1, 1] == pytest.approx(.004-3e-6)
    assert np.max(abs(result["v_nonlinear"])) == 0.


def test_stage_c_arrays_and_interpolation_are_labelled(tmp_path):
    disc = SyntheticDisc(64); case = save_dense_synthetic_case(tmp_path/"one_d", disc)
    H = comparison.DEFAULT_T1/4
    time, x = np.linspace(0., H, 9), np.linspace(0., 1., 5)
    saved = comparison.evaluate_saved_one_d(case, disc, time)
    linear = disc.reconstruct_series(saved["q_linear"], x)
    nonlinear = disc.reconstruct_series(saved["q_nonlinear"], x)
    raw_x = np.array((.2, .5, .8))
    rotations = np.broadcast_to(np.eye(3), (len(time), len(raw_x), 3, 3)).copy()
    first = {"time": time, "x": x, "fields": linear, "raw_rotation_x": raw_x, "raw_rotation_matrices": rotations}
    second = {**first, "fields": nonlinear}
    output = tmp_path/"comparison"
    summary = comparison.compare_fem_pair(first, second, case, disc, output)
    assert summary["scientific_calls"] == 0 and summary["time_grid"]["count"] == 201
    assert summary["native_samples_and_interpolated_values_are_distinct"]
    assert "not a native" in summary["zero_state_source"]
    assert "proxy" in summary["metrics"]["nonlinear"]["c"]["qualification"]
    assert summary["physical_model_accuracy_status"] == "NOT_AUTOMATICALLY_ASSIGNED"
    with np.load(output/"one_d_three_d_comparison.npz", allow_pickle=False) as z:
        assert np.max(abs(z["one_d_evolution"][0])) == 0.
        assert z["three_d_R_nonlinear"].shape == (201, 5, 3, 3)


def test_native_zero_source_is_static_not_dynamic_frame(tmp_path):
    x = np.linspace(0., 1., 5); raw_x = np.array((.2, .5, .8))
    initial = np.zeros((5, 7)); initial[:, 1] = .004*np.sin(np.pi*x)
    np.savez_compressed(tmp_path/"initial_sections.npz", x=x, fields=initial,
                        raw_rotation_x=raw_x, raw_rotation_matrices=np.broadcast_to(np.eye(3), (3, 3, 3)))
    np.savez_compressed(tmp_path/"section_history.npz", time=np.array((.01, .02)), x=x,
                        fields=np.stack((initial*.999, initial*.998)), raw_rotation_x=raw_x,
                        raw_rotation_matrices=np.broadcast_to(np.eye(3), (2, 3, 3, 3)))
    history = comparison.read_fem_trajectory(tmp_path)
    assert np.array_equal(history["time"], (0., .01, .02))
    assert not history["zero_frame_is_native_dynamic"]
    np.testing.assert_array_equal(history["fields"][0], initial)


def test_four_figures_use_saved_arrays_and_no_new_scientific_calls(stage_b, monkeypatch):
    cases, discs, output = stage_b
    comparison.analyze_stage_b(cases, discs, output, chunk_size=2)
    monkeypatch.setattr(rod, "derive_polynomials", lambda: pytest.fail("Plotting must not derive an action"))
    figures = comparison.render_figures(output)
    assert len(figures) == 8
    assert len({Path(name).stem for name in figures}) == 4
    assert all((output/name).exists() for name in figures)
    evidence = json.loads((output/"figure_data_provenance.json").read_text(encoding="utf8"))
    assert evidence["scientific_calls"] == 0 and evidence["3D_level"] is None
    assert evidence["three_dimensional_view_uses_equal_physical_scale"]
    assert not evidence["source_3D_samples_are_explicitly_postprocessed"]


def test_independent_mixed_p_control_does_not_hide_velocity_partial(stage_b):
    cases, discs, output = stage_b
    for field in ("w", "v"):
        source = cases["isolated_"+field+"_p64"]
        with np.load(source/"trajectory.npz", allow_pickle=False) as z:
            times, nonlinear = z["times"], z["q"].copy()
        with np.load(source/"linear_trajectory.npz", allow_pickle=False) as z:
            linear = z["q"]
        nonlinear[:, 1 if field == "w" else 2] += 1e-11*times
        name = "isolated_"+field+"_p48"; cases[name] = source.parent/name
        save_case(cases[name], times, nonlinear, linear)
    # A synthetic velocity discrepancy leaves the displacement signal and
    # mixed-response p control unchanged. Its inherited status must survive.
    source = cases["joint_p48"]
    with np.load(source/"trajectory.npz", allow_pickle=False) as z:
        rows = {name: z[name].copy() for name in z.files}
    rows["velocity"][:, 6] += 1e-4*rows["times"]
    np.savez_compressed(source/"trajectory.npz", **rows)
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir() if p.name != "case.json"}
    (source/"case.json").write_text(json.dumps({"status": "PASS", "artifact_hashes": files}), encoding="utf8")
    summary = comparison.analyze_stage_b(cases, discs, output, chunk_size=2)
    assert summary["status"] == "NUMERICAL_PARTIAL"
    assert summary["all14_spatial"]["velocity"]["status"] == "PARTIAL"
    assert summary["stage_c_allowed"]
    assert summary["independent_mixed_p_control_available"]
    assert not summary["mixed_response_p_sensitivity_uses_joint_control_only"]
    assert summary["spatial_mixed_evolution"]["fields"]["w"]["absolute_max"] > 0.
    assert summary["spatial_velocity_correction"]["fields"]["c"]["absolute_max"] > 0.
    assert len(summary["source_cases"]) == 6


def test_one_extra_isolated_p_control_is_not_sufficient(stage_b):
    cases, discs, output = stage_b
    cases["isolated_w_p48"] = cases["isolated_w_p64"]
    with pytest.raises(ValueError, match="Both isolated p48"):
        comparison.analyze_stage_b(cases, discs, output)


def test_physical_curvature_distinguishes_Phi_from_its_gradient():
    x = np.linspace(0., 1., 101)
    fields = np.zeros((1, len(x), 7)); gradients = fields.copy()
    fields[0, :, 3] = .002*x; gradients[0, :, 3] = .002
    exact = comparison.exact_diagnostic_curvatures(fields, gradients)
    recovered = comparison.physical_curvature_from_rotations(comparison.physical_rotations(fields), x)
    np.testing.assert_allclose(exact[..., 0], .002, rtol=3e-15, atol=1e-17)
    np.testing.assert_allclose(recovered["interval_chi_geodesic"][..., 0], .002, rtol=2e-12, atol=1e-16)
    np.testing.assert_allclose(recovered["chi_gradient"][..., 0], .002, rtol=2e-9, atol=1e-15)
    assert recovered["Phi_is_not_chi1"]
    assert fields[0, 0, 3] == 0. and exact[0, 0, 0] != 0.


def test_uniform_rigid_rotation_has_no_orientation_curvature():
    x = np.linspace(0., 1., 21)
    fields = np.zeros((2, len(x), 7)); fields[..., 3:6] = (.1, -.03, .04)
    rotations = comparison.physical_rotations(fields)
    result = comparison.physical_curvature_from_rotations(rotations, x)
    assert np.max(abs(result["chi_gradient"])) < 2e-14
    assert np.max(abs(result["interval_chi_geodesic"])) < 1e-14
    assert np.max(abs(comparison.exact_diagnostic_curvatures(fields, np.zeros_like(fields)))) == 0.


def test_mesh_comparison_reads_two_saved_common_grids_only(tmp_path):
    disc = SyntheticDisc(64); case = save_dense_synthetic_case(tmp_path/"one_d", disc)
    H = comparison.DEFAULT_T1/4
    time, x = np.linspace(0., H, 9), np.linspace(0., 1., 5)
    saved = comparison.evaluate_saved_one_d(case, disc, time)
    raw_x = np.array((.2, .5, .8))
    rotations = np.broadcast_to(np.eye(3), (len(time), len(raw_x), 3, 3)).copy()
    linear = {"time": time, "x": x, "fields": disc.reconstruct_series(saved["q_linear"], x),
              "raw_rotation_x": raw_x, "raw_rotation_matrices": rotations}
    nonlinear = {**linear, "fields": disc.reconstruct_series(saved["q_nonlinear"], x)}
    assert comparison.analyze_mesh_comparison(tmp_path)["status"] == "NOT_RUN"
    for level in ("medium", "fine"):
        comparison.compare_fem_pair(linear, nonlinear, case, disc, tmp_path/("comparison_"+level))
    result = comparison.analyze_mesh_comparison(tmp_path)
    assert result["scientific_calls"] == 0
    assert result["mesh_metrics"]["evolution"]["w"]["absolute_max"] == 0.
    assert result["physical_orientation_mesh_changes"]["nonlinear"]["maximum_principal_angle_radians"] == 0.
    assert result["section_orientation_curvature"]["nonlinear"]["levels"]["fine"]["twist_proxy_max"] == 0.
    assert result["mesh_change_is_not_a_continuum_error_bound"]
    assert result["single_3D_time_level_temporal_certification"] == "PARTIAL"
    assert (tmp_path/"mesh_comparison.npz").exists()


def test_saved_operator_initial_state_decomposition_identity(tmp_path, monkeypatch):
    discs = {48: SyntheticDisc(48), 64: SyntheticDisc(64)}
    cases = {}
    for degree in (48, 64):
        for kind in ("joint", "isolated_w", "isolated_v"):
            name = kind+"_p"+str(degree)
            cases[name] = save_dense_synthetic_case(tmp_path/name, discs[degree])
    monkeypatch.setattr(np.linalg, "eigh", lambda *a, **kw: pytest.fail("No new eigenanalysis"))
    result = comparison.write_initial_state_decomposition(cases, tmp_path/"decomposition", discs=discs)
    assert result["scientific_calls"] == 0
    assert result["operation_counts"]["saved_full_linear_operator_evaluations"] == 12
    assert result["pre_FEM_decision_unchanged"] and result["diagnostic_not_a_new_stage_C_gate"]
    for name in ("joint_p48", "joint_p64", "mixed_p48", "mixed_p64"):
        assert result["components"][name]["physical_identity_max_absolute"] < 1e-15
    assert result["components"]["joint_p64"]["fields"]["w"]["same_initial_nonlinear_max"] > 0.
    with np.load(tmp_path/"decomposition"/"initial_state_decomposition.npz", allow_pickle=False) as z:
        np.testing.assert_allclose(z["mixed_p64_total_q"],
            z["mixed_p64_propagated_initial_q"]+z["mixed_p64_same_initial_nonlinear_q"], atol=2e-18)


def test_basis_only_saved_shen_reconstruction_without_action_derivation(tmp_path, monkeypatch):
    from scripts.lib.weakly_nonlinear_spatial_dynamics import SpatialGalerkin
    from types import SimpleNamespace
    root = Path(__file__).resolve().parents[1]
    action = root/"results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json"
    if not action.exists():
        pytest.skip("Frozen action unavailable")
    pol = json.loads(action.read_text(encoding="utf8"))["polynomials"]
    model = SimpleNamespace(T4=rod.Polynomial.deserialize(pol["T4"]), V4=rod.Polynomial.deserialize(pol["V4"]),
                            residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol["residuals_A"]))
    coefficients = rod.RodCoefficients.rectangular(1., 1., .3, .2, .1, 1.759089824002232e-5)
    disc = SpatialGalerkin(coefficients, 6, model=model)
    q = disc.from_raw_coefficients(np.arange(disc.ndof)*1e-7)
    np.savez_compressed(tmp_path/"static_states.npz", q_linear=q, q_nonlinear=q,
                        raw_linear=disc.raw_coefficients(q), raw_nonlinear=disc.raw_coefficients(q))
    files = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in tmp_path.iterdir()}
    (tmp_path/"case.json").write_text(json.dumps({"fields": list(comparison.FIELDS), "ndof": disc.ndof,
        "request": {"p": disc.p, "nq": disc.nq, "length": 1., "whiten": True,
                    "coefficients": coefficients.values()}, "artifact_hashes": files}), encoding="utf8")
    monkeypatch.setattr(rod, "derive_polynomials", lambda: pytest.fail("Basis reconstruction must not derive an action"))
    adapter = comparison.SavedShenReconstruction(tmp_path)
    x = np.linspace(0., 1., 101)
    np.testing.assert_allclose(adapter.reconstruct_series(q[None, :], x), disc.reconstruct_series(q[None, :], x), rtol=5e-13, atol=1e-16)


def test_support_moment_origin_transfer_matches_direct_force_cloud():
    points = np.array(((0., -.05, .03), (0., .02, -.07), (0., .03, .04)))
    forces = np.array(((1e-8, 2e-6, 3e-6), (-2e-8, -1e-6, 4e-6), (1e-8, 3e-6, -2e-6)))
    old_origin = points.mean(axis=0)
    new_origin = np.array((0., 0., 0.))
    old_moment = np.cross(points-old_origin, forces).sum(axis=0)
    direct = np.cross(points-new_origin, forces).sum(axis=0)
    corrected = comparison.transfer_moment_origin(old_moment, forces.sum(axis=0), old_origin, new_origin)
    np.testing.assert_allclose(corrected, direct, rtol=2e-15, atol=1e-21)
    np.testing.assert_allclose(
        comparison.transfer_moment_origin(corrected, forces.sum(axis=0), new_origin, old_origin),
        old_moment, rtol=2e-15, atol=1e-21)


def test_off_axis_force_centroid_does_not_imply_beam_axis_torque():
    force = np.array((0., 2e-5, 3e-5))
    node_mean = np.array((0., .002, -.001))
    axis_center = np.zeros(3)
    # A pure resultant applied on the beam axis has zero moment there.
    moment_about_node_mean = np.cross(axis_center-node_mean, force)
    assert abs(moment_about_node_mean[0]) > 1e-9
    np.testing.assert_allclose(
        comparison.transfer_moment_origin(moment_about_node_mean, force, node_mean, axis_center),
        np.zeros(3), atol=1e-22)

