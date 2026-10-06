"""Composition, target-prefix regression and separate historical qualification.

Reuse this diagnostic's saved Gate-A pair only when its source hashes,
parameters and settings match; missing data never trigger a root search.
No physical model, root search or solver internals are reproduced here.
"""
from dataclasses import asdict
from copy import deepcopy
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch
import shutil

import numpy as np
import pytest
from numpy.testing import assert_allclose

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.analysis.joint_review import check_circular_eb_spring_spectrum as check


@pytest.fixture(autouse=True)
def no_new_pilot_root_search(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("This saved-spectrum test must not search for roots")
    monkeypatch.setattr(check, "resolve_matrix_spectrum", forbidden)


@pytest.fixture
def params():
    return check.BeamParams(E=2.1e11, rho=7800., r=.005, L_total=2.)


def test_baseline_epsilon(params):
    assert_allclose(params.eps, params.r / (2 * params.L_base), rtol=1e-14)


@pytest.mark.parametrize("state", [1, 10, 100, "RIGID"])
def test_composition_uses_baseline_lengths_frequency_and_exact_joint(params, state):
    marker = np.eye(6)
    with patch.object(check, "segment_lengths", wraps=check.segment_lengths) as lengths, \
         patch.object(check, "lambdas_to_frequencies", wraps=check.lambdas_to_frequencies) as conversion, \
         patch.object(check, "boundary_assembly", return_value=SimpleNamespace(dimensionless=marker)) as assembly:
        _, provider = check.providers(params, .30, 15., state)
        assert provider(3.25) is marker
    lengths.assert_called_once_with(params, .30)
    conversion.assert_called_once()
    assert_allclose(conversion.call_args.args[0], [3.25], rtol=0, atol=0)
    assert conversion.call_args.args[1] is params
    omega, left, right, beta, joint, reference = assembly.call_args.args
    assert (left.L, right.L) == check.segment_lengths(params, .30)
    for arm in (left, right, reference):
        assert (arm.A, arm.D, arm.m) == (params.E * params.S, params.E * params.I, params.rho * params.S)
    assert reference.L == params.L_base
    assert_allclose(beta, np.deg2rad(15.))
    assert_allclose(omega, 3.25**2 * np.sqrt(reference.D / reference.m) / reference.L**2, rtol=1e-14)
    if state == "RIGID":
        assert joint.mode == "RIGID" and joint.k_theta is None
    else:
        assert joint.mode == "SPRING"
        assert joint.k_theta == state * reference.D / reference.L


def test_baseline_provider_forwards_frozen_matrix(params):
    marker = np.eye(6)
    with patch.object(check, "assemble_clamped_coupled_matrix", return_value=marker) as matrix:
        baseline, _ = check.providers(params, .30, 15.)
        assert baseline(3.25) is marker
    matrix.assert_called_once_with(3.25, np.deg2rad(15.), .30, params.eps)


@pytest.fixture(scope="module")
def rigid_results():
    params = check.BeamParams(E=2.1e11, rho=7800., r=.005, L_total=2.)
    settings = check.SearchSettings()
    checkpoint = check.OUTPUT / "diagnostics.json"
    payload = json.loads(checkpoint.read_text(encoding="utf-8")) if checkpoint.exists() else {}
    runs = payload.get("runs", {})
    reusable = (
        all(payload.get("source_hashes", {}).get(p) == check.source_hashes()[p] for p in check.SOURCE_PATHS[:-1])
        and payload.get("provider_source_sha256") == check.provider_source_hash()
        and payload.get("params") == asdict(params)
        and payload.get("settings") == asdict(settings)
        and all(key in runs for key in ("main_baseline", "main_RIGID"))
    )
    if not reusable:
        pytest.fail("Source-matched saved Gate A required; root recalculation is forbidden")
    results = [runs[key]["result"] for key in ("main_baseline", "main_RIGID")]
    for result in results:
        assert result["geometry"] == asdict(check.Geometry(params.eps, 15., .30, 0.))
        assert len(result["roots"]) >= check.K_GUARD
    return results


def test_rigid_target_prefix_regression(rigid_results):
    baseline, rigid = rigid_results
    comparison, rows = check.rigid_comparison(baseline, rigid)
    assert comparison["target_prefix_status"] == "RIGID_TARGET_PREFIX_PASS"
    assert all(r["within_root_match_tol"] for r in rows[:check.K_GUARD])
    assert all(check.target_prefix(r)["target_prefix_status"] == "TARGET_PREFIX_PASS" for r in rigid_results)


def test_higher_spectrum_qualification_is_retained_metadata():
    checkpoint = check.OUTPUT / "diagnostics.json"
    if not checkpoint.exists():
        pytest.skip("Historical qualification requires the saved diagnostic artifact")
    data = json.loads(checkpoint.read_text(encoding="utf-8"))
    history = data.get("previous_full_spectrum_gate")
    if history is None:
        pytest.skip("This checkpoint does not contain the historical full-spectrum audit")
    assert history["status"] == "STOP_GATE_A"
    assert history["audit"]["test_result"]["failed"] == 1
    assert not history["gates"]["A"][0]["passed"]
    baseline, rigid = (data["runs"][key]["result"] for key in ("main_baseline", "main_RIGID"))
    comparison, _ = check.rigid_comparison(baseline, rigid)
    assert not comparison["full12_equivalent"]  # This asserts the qualification, NOT equivalence.
    assert comparison["higher_mismatched_positions"] == [11, 12]
    assert rigid["roots"][10]["detected_nullity"] == 2
    assert rigid["roots"][10]["Lambda"] == rigid["roots"][11]["Lambda"]
    assert baseline["roots"][10]["detected_nullity"] == 1


@pytest.mark.parametrize("above", [False, True])
def test_native_unresolved_interval_is_scoped_to_guard_margin(rigid_results, above):
    result = deepcopy(rigid_results[1])
    limit = check.target_prefix(result)["guard_audit_limit"]
    left = limit + (1 if above else -1) * result["settings"]["root_match_tol"]
    result["primary"]["unresolved_intervals"].append(f"{left}:{left + .01}:native_test_record")
    result["spectrum_status"] = "unresolved"
    result["exclusion_reason"] = "unresolved_low_sigma_interval"
    prefix = check.target_prefix(result)
    assert prefix["solver_spectrum_status"] == "unresolved"
    assert result["spectrum_status"] == "unresolved"
    assert prefix["target_prefix_status"] == ("TARGET_PREFIX_PASS" if above else "TARGET_PREFIX_FAIL")


def test_missing_high_spectrum_reserve_does_not_replace_prefix_gate(rigid_results):
    result = deepcopy(rigid_results[1])
    result["roots"] = result["roots"][:8]
    for name in ("primary", "verification"):
        result[name]["roots"] = result[name]["roots"][:8]
    for row in result["primary_vs_verification"][8:]:
        row["status"] = "missing_root"
        row["multiplicity_agreement"] = False
    result.update(spectrum_status="unresolved", independent_agreement=False,
        root11_available=False, root12_available=False, root12_boundary_warning=True,
        exclusion_reason="found_only_8_of_12;unresolved_independent_search_disagreement;candidate_boundary_warning")
    prefix = check.target_prefix(result)
    assert prefix["solver_spectrum_status"] == "unresolved"
    assert prefix["target_prefix_status"] == "TARGET_PREFIX_PASS"
    assert prefix["higher_spectrum_status"] == "HIGHER_SPECTRUM_QUALIFICATION"


@pytest.mark.parametrize("problem", ["missing_guard", "disagreement", "nullity", "quality", "order", "coverage", "candidate", "unlocalized_interval"])
def test_prefix_rejects_each_native_target_problem(rigid_results, problem):
    result = deepcopy(rigid_results[1])
    primary = result["primary"]
    if problem == "missing_guard":
        primary["roots"] = primary["roots"][:6]
    elif problem == "disagreement":
        result["primary_vs_verification"][0]["status"] = "disagreement"
    elif problem == "nullity":
        primary["roots"][0]["detected_nullity"] = 2
    elif problem == "quality":
        primary["roots"][0]["sigma_1"] = 2 * result["settings"]["sigma_accept"]
    elif problem == "order":
        primary["roots"][0], primary["roots"][1] = primary["roots"][1], primary["roots"][0]
    elif problem == "coverage":
        primary["lambda_upper"] = check.target_prefix(result)["guard_audit_limit"]
    elif problem == "candidate":
        primary["candidates"][0]["acceptance_status"] = "unresolved_test_record"
    else:
        primary["unresolved_intervals"].append("unlocalized")
    assert check.target_prefix(result)["target_prefix_status"] == "TARGET_PREFIX_FAIL"


@pytest.fixture
def recovered_saved_mode():
    params, saved, _ = check.saved_shape_inputs(check.OUTPUT)
    _, provider = check.providers(params, .30, 15., "1")
    row, shape, diagnostic = check.reconstruct_saved_mode(params, "1", saved["1"][0], provider.case)
    return row, shape, diagnostic, provider.case


def test_saved_simple_root_recovers_without_root_search(recovered_saved_mode):
    row, shape, diagnostic, _ = recovered_saved_mode
    assert row["reconstruction_status"] == "CONFIRMED"
    assert row["detected_nullity"] == 1
    assert not diagnostic["failures"]
    assert shape["states"].shape == (2, 129, 6)
    assert row["sigma_ratio"] <= check.SHAPE_GATES["sigma_ratio"]


def test_normalization_uses_both_actual_arm_masses(recovered_saved_mode):
    row, shape, diagnostic, case = recovered_saved_mode
    arms = case[:2]
    assert [a.L for a in arms] == [.7, 1.3]
    _, weights = check.modes.quadrature(129)
    physical = np.concatenate([check.modes.mass_vector(y[None, ...], arm, weights)
                               for y, arm in zip(shape["states"], arms)])
    assert_allclose(shape["physical_vector"], physical, rtol=1e-14, atol=1e-14)
    assert np.vdot(physical, physical).real == pytest.approx(1., abs=1e-14)
    assert row["mass_before_normalization"] == pytest.approx(sum(diagnostic["arm_masses_before_normalization"]))
    assert min(diagnostic["arm_masses_before_normalization"]) > 0


def test_sensitivity_uses_existing_formula_and_normalized_mass(recovered_saved_mode):
    row, _, _, case = recovered_saved_mode
    reference = case[-1]
    expected = reference.D/reference.L * row["Delta_psi_mass_normalized"]**2 / (
        (2*np.pi*row["frequency_hz"])**2 * row["mass_after_normalization"])
    assert row["s"] == pytest.approx(expected, rel=1e-14)


def test_all_saved_rigid_modes_meet_constraint_and_have_no_finite_s():
    params, saved, _ = check.saved_shape_inputs(check.OUTPUT)
    _, provider = check.providers(params, .30, 15., "RIGID")
    for root in saved["RIGID"]:
        row, shape, diagnostic = check.reconstruct_saved_mode(params, "RIGID", root, provider.case)
        assert row["reconstruction_status"] == "CONFIRMED"
        assert abs(diagnostic["normalized_reconstructed_physical_residuals"][2]) <= check.SHAPE_GATES["physical_residual"]
        assert row["Delta_psi_mass_normalized"] == shape["states"][0, -1, 2] - shape["states"][1, -1, 2]
        assert row["s"] == "NOT_APPLICABLE_RIGID"


def test_mass_mac_permutation_can_select_guard_and_ignores_frequency():
    vectors = np.eye(7)
    permutation = [6, 1, 2, 3, 4, 5, 0]
    left = [dict(physical_vector=vectors[i], frequency_hz=i+1) for i in range(6)]
    right = [dict(physical_vector=vectors[i], frequency_hz=j+1) for j, i in enumerate(permutation)]
    columns, mac, margin = check.assign_physical_modes(left, right)
    assert columns.tolist() == [6, 1, 2, 3, 4, 5]
    assert_allclose(mac[np.arange(6), columns], 1)
    assert_allclose(margin, 1)
    for row in right:
        row["frequency_hz"] = 1e12  # No frequency enters the assignment cost.
    changed, _, _ = check.assign_physical_modes(left, right)
    assert_allclose(changed, columns)


def test_missing_guard_aborts_shapes_instead_of_recalculating(tmp_path):
    data = json.loads((check.OUTPUT / "diagnostics.json").read_text(encoding="utf-8"))
    data["runs"]["main_kappa_1"]["result"]["roots"] = data["runs"]["main_kappa_1"]["result"]["roots"][:6]
    (tmp_path / "diagnostics.json").write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError, match="MISSING_SAVED_TARGET_PREFIX"):
        check.saved_shape_inputs(tmp_path)


def test_shapes_orchestration_preserves_frequency_inputs_guard_pool_and_phase(tmp_path):
    names = ("diagnostics.json", "rigid_equivalence.csv", "spring_spectrum.csv", "trend_summary.csv", "report.md")
    for name in names:
        shutil.copyfile(check.OUTPUT / name, tmp_path / name)
    original = {name: (tmp_path / name).read_bytes() for name in names[:-1]}
    result = check.shapes_run(tmp_path)
    assert result["summary"]["root_search_calls"] == 0
    assert result["summary"]["confirmed_modes"] == 28
    assert all((tmp_path / name).read_bytes() == contents for name, contents in original.items())
    assert all(a["candidate_sorted_positions"] == list(range(1, 8)) for a in result["assignments"])
    with np.load(tmp_path / "shapes.npz") as archive:
        records = {(f"{r['kappa_theta']:g}" if r["joint_mode"] == "SPRING" else "RIGID", r["sorted_position"]): r
                   for r in result["reconstruction"]}
        for row in result["reconstruction"]:
            states = archive[row["shape_key"] + "__states"]
            assert row["Delta_psi_mass_normalized"] == states[0, -1, 2] - states[1, -1, 2]
            vector = archive[row["shape_key"] + "__physical_vector"]
            assert np.vdot(vector, vector).real == pytest.approx(1., abs=1e-14)
        for row in result["mapping"]:
            if row["mapping_status"] != "CONFIRMED":
                continue
            source = records[row["source_kappa"], row["source_sorted_position"]]["shape_key"]
            target = records[row["target_kappa"], row["target_sorted_position"]]["shape_key"]
            assert np.vdot(archive[source + "__physical_vector"], archive[target + "__physical_vector"]).real >= 0
    for seed in result["summary"]["unresolved_descendants"]:
        assert next(r for r in result["trends"] if r["branch_id"] == seed)["rotation_trend_status"] == "UNRESOLVED"


def continuation_pool(angle, permutation=None):
    """Synthetic orthonormal mass vectors only; no physical/root implementation."""
    vectors = np.eye(7)
    cosine, sine = np.cos(np.deg2rad(angle)), np.sin(np.deg2rad(angle))
    vectors[:2, :2] = [[cosine, sine], [-sine, cosine]]
    order = list(range(7)) if permutation is None else permutation
    return [dict(physical_vector=vectors[i].copy(), states=np.ones((2, 2, 6)), reactions=np.ones((2, 3)),
        psi1_mass_normalized=1., psi2_mass_normalized=0., Delta_psi_mass_normalized=1.,
        sorted_position=j+1, branch_id=f"k1_seed_{j+1:02d}", frequency_hz=1e9/(j+1), s=99.)
        for j, i in enumerate(order)]


def test_geometric_midpoint_and_rigid_is_not_numeric():
    assert check.geometric_midpoint(1, 10) == np.sqrt(10)
    assert check.geometric_midpoint(10, 100) == np.sqrt(1000)
    with pytest.raises(ValueError, match="RIGID_IS_NOT_A_FINITE_KAPPA"):
        check.geometric_midpoint(100, "RIGID")
    with pytest.raises(ValueError, match="RIGID_IS_NOT_A_FINITE_KAPPA"):
        check.continue_finite_interval([], 100, "RIGID", None, [])


def test_continuation_delegates_mass_only_assignment_and_keeps_guard(monkeypatch):
    assert check.SHAPE_GATES["MAC"] == .95
    assert check.SHAPE_GATES["margin"] == .20
    existing = check.modes.assign
    calls = []

    def delegate(left, right):
        calls.append((left, right))
        assert all(isinstance(v, np.ndarray) for v in left + right)
        return existing(left, right)

    monkeypatch.setattr(check.modes, "assign", delegate)
    source = continuation_pool(0)[:6]
    target = continuation_pool(0, [6, 1, 2, 3, 4, 5, 0])
    before = deepcopy(target)
    selected, rows = check.continuation_assignment(source, target, 1, 10, 1, [])
    assert len(calls[0][0]) == 6 and len(calls[0][1]) == 7
    assert [r["target_sorted_position"] for r in rows] == [7, 2, 3, 4, 5, 6]
    for mode in target:
        mode.update(frequency_hz=-9e15, s=-9e20, Delta_psi_mass_normalized=0.)
    changed, _ = check.continuation_assignment(source, target, 1, 10, 1, [])
    assert [r["sorted_position"] for r in changed] == [r["sorted_position"] for r in selected]
    for old, now in zip(before, target):
        assert_allclose(old["physical_vector"], now["physical_vector"])
    with pytest.raises(ValueError, match="GUARD_7"):
        check.continuation_assignment(source, target[:6], 1, 10, 1, [])


def test_direct_failure_remains_historical_after_two_local_passes():
    calls = []

    def point(kappa):
        calls.append(kappa)
        return continuation_pool(22*np.log10(kappa))

    source = continuation_pool(0)[:6]
    _, direct_rows = check.continuation_assignment(source, point(10), 1, 10, 0, [])
    direct = dict(direct_rows[0], mapping_status=direct_rows[0]["step_status"])
    unchanged = deepcopy(direct)
    assert direct["MAC"] < .95 and direct["mapping_status"] == "UNRESOLVED"
    attempts = []
    final, steps, path = check.continue_finite_interval(source, 1, 10, point, attempts, force_midpoint=True)
    assert path == [1, np.sqrt(10), 10]
    assert all(r["step_status"] == "CONFIRMED" for step in steps for r in step)
    status, agreement = check.continuation_status(direct, final[0], [step[0] for step in steps])
    assert status == "DIRECT_UNRESOLVED_CONTINUATION_CONFIRMED"
    assert agreement == "DIRECT_ENDPOINT_AGREEMENT"
    assert direct == unchanged  # Continuation cannot make direct MAC pass.


def test_continuation_phase_alignment_changes_only_copies():
    source, target = continuation_pool(0)[:6], continuation_pool(0)
    for key in ("physical_vector", "states", "reactions"):
        target[0][key] *= -1
    for key in ("psi1_mass_normalized", "psi2_mass_normalized", "Delta_psi_mass_normalized"):
        target[0][key] *= -1
    snapshot = deepcopy(target)
    selected, rows = check.continuation_assignment(source, target, 1, 10, 1, [])
    assert rows[0]["phase_sign"] == -1
    assert np.vdot(source[0]["physical_vector"], selected[0]["physical_vector"]).real > 0
    assert selected[0]["Delta_psi_mass_normalized"] == 1.
    for before, after in zip(snapshot, target):
        for key in ("physical_vector", "states", "reactions", "Delta_psi_mass_normalized"):
            assert_allclose(before[key], after[key])


def test_endpoint_disagreement_is_a_path_conflict():
    direct = dict(mapping_status="UNRESOLVED", target_sorted_position=1)
    status, agreement = check.continuation_status(direct, dict(sorted_position=2), [dict(step_status="CONFIRMED")])
    assert status == "PATH_MAPPING_CONFLICT"
    assert agreement == "DIRECT_ENDPOINT_DISAGREEMENT"


def test_refinement_stops_at_quarters_even_when_still_unresolved():
    calls, attempts = [], []

    def point(kappa):
        calls.append(kappa)
        return continuation_pool(120*np.log10(kappa))

    _, steps, path = check.continue_finite_interval(continuation_pool(0)[:6], 1, 10, point, attempts, force_midpoint=True)
    new = sorted(set(k for k in calls if k != 10))
    assert_allclose(new, [10**.25, 10**.5, 10**.75], rtol=1e-15)
    assert len(new) == 3
    assert max(r["refinement_depth"] for a in attempts for r in a["rows"]) == 2
    assert any(r["step_status"] == "UNRESOLVED" for step in steps for r in step)
    assert len(path) == 5
    with pytest.raises(ValueError, match="DEPTH_BUDGET"):
        check.continue_finite_interval([], 1, 10, point, attempts, depth=3)


def test_only_failed_half_is_refined():
    calls = []

    def point(kappa):
        calls.append(kappa)
        return continuation_pool(60*max(0, np.log10(kappa)-.5))

    check.continue_finite_interval(continuation_pool(0)[:6], 1, 10, point, [], force_midpoint=True)
    new = sorted(set(k for k in calls if k != 10))
    assert_allclose(new, [10**.5, 10**.75], rtol=1e-15)


@pytest.mark.parametrize("final_angle,middle_angle,expected_bridge", [(10, 5, False), (20, 10, True), (60, 30, True)])
def test_exact_rigid_direct_or_only_one_1000_bridge(final_angle, middle_angle, expected_bridge):
    calls, attempts = [], []

    def point(kappa):
        calls.append(kappa)
        assert kappa in ("RIGID", 1000.)
        return continuation_pool(final_angle if kappa == "RIGID" else middle_angle)

    _, steps, path = check.continue_rigid_endpoint(continuation_pool(0)[:6], point, attempts)
    assert (1000. in calls) == expected_bridge
    assert path == ([100., 1000., "RIGID"] if expected_bridge else [100., "RIGID"])
    if final_angle == 60:
        assert any(r["step_status"] == "UNRESOLVED" for step in steps for r in step)
    else:
        assert all(r["step_status"] == "CONFIRMED" for step in steps for r in step)


def copy_continuation_inputs(destination):
    for name in (*check.CONTINUATION_INPUTS, "report.md", "kappa_continuation_diagnostics.json"):
        shutil.copyfile(check.OUTPUT / name, destination / name)


@pytest.mark.parametrize("failure", ["prefix", "reconstruction"])
def test_new_point_failure_stops_before_any_further_point(tmp_path, monkeypatch, failure):
    copy_continuation_inputs(tmp_path)
    path = tmp_path / "kappa_continuation_diagnostics.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    key = repr(float(np.sqrt(10)))
    data["points"] = {key: data["points"][key]}
    if failure == "prefix":
        data["points"][key]["prefix"]["target_prefix_status"] = "TARGET_PREFIX_FAIL"
    else:
        data["points"][key]["reconstruction"][0]["reconstruction_status"] = "FAILED"
    path.write_text(json.dumps(data), encoding="utf-8")

    def forbidden(*args, **kwargs):
        raise AssertionError("A failed cached point must not be repaired or skipped")

    monkeypatch.setattr(check, "reconstruct_saved_mode", forbidden)
    result = check.continuation_run(tmp_path)
    assert result["status"] == "STOP_CONTINUATION"
    assert result["summary"]["full_confirmed_count"] == 0
    assert result["solver_calls_this_invocation"] == 0
    assert list(result["points"]) == [key]
    assert result["mapping"] == []


def test_completed_continuation_replays_without_roots_or_shapes_and_preserves_inputs(tmp_path, monkeypatch):
    copy_continuation_inputs(tmp_path)
    cached = json.loads((tmp_path / "kappa_continuation_diagnostics.json").read_text(encoding="utf-8"))
    assert cached["status"] == "COMPLETED_BOUNDED_CONTINUATION"
    before = {name: (tmp_path / name).read_bytes() for name in check.CONTINUATION_INPUTS}
    report = (tmp_path / "report.md").read_bytes()
    for marker in (b"\n## Bounded kappa continuation\n", b"\r\n## Bounded kappa continuation\r\n"):
        report = report.split(marker)[0]

    def forbidden(*args, **kwargs):
        raise AssertionError("Saved points must not be reconstructed again")

    monkeypatch.setattr(check, "reconstruct_saved_mode", forbidden)
    result = check.continuation_run(tmp_path)
    assert result["solver_calls_this_invocation"] == 0
    assert result["summary"]["full_confirmed_count"] == cached["summary"]["full_confirmed_count"]
    assert all((tmp_path / name).read_bytes() == content for name, content in before.items())
    assert (tmp_path / "report.md").read_bytes().startswith(report)
    assert result["DIRECT_1_TO_10"] == cached["DIRECT_1_TO_10"]
    for attempt in result["attempts"]:
        assert attempt["candidate_sorted_positions"] == list(range(1, 8))
        assert len({r["target_sorted_position"] for r in attempt["rows"]}) == len(attempt["rows"])
    for row in result["trends"]:
        if row["full_path_status"] == "FULL_KAPPA_PATH_CONFIRMED":
            assert row["minimum_step_MAC"] >= .95
            assert row["minimum_step_margin"] >= .20
    for point in result["points"].values():
        assert not float(point["kappa_theta"]).is_integer()  # Current run needed finite log midpoints only.
        assert all(row["kappa_theta"] == point["kappa_theta"] for row in point["reconstruction"])
