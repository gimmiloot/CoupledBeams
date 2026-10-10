"""FEM-3B algebra, scheduling and saved-data contracts; no native/ODE jobs."""
from __future__ import annotations

import copy
import hashlib
import importlib
import json
from types import SimpleNamespace
import time
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate._ivp.radau import RadauDenseOutput

from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
PARENT_SHA = "187870ece572dfe1b999b83d99016d01c1fdbe10e8835ed45543d1119d6f36ea"
NEW_BUNDLE = ROOT / "results/nlsp_nonlinear_dynamic_long_horizon/7d2b499e6a1eb990"


@pytest.fixture
def config():
    return json.loads((ROOT / "data/input/nlsp_nonlinear_dynamic_long_horizon.json").read_text(encoding="utf8"))


@pytest.fixture
def continuation():
    return importlib.import_module("scripts.lib.nlsp_fem3b_continuation")


class _FakeDisc:
    """The observer test concerns storage, so it never evaluates a model RHS."""
    ndof = 2
    p = 3
    nq = 7

    def reset_counters(self):
        pass

    def counters(self):
        return {"fake_RHS_calls": 0}

    def rhs(self, *args):
        pytest.fail("The fake observer test must not integrate an actual RHS")

    def jacobian(self, *args):
        pytest.fail("The fake observer test must not evaluate a Jacobian")


class _FakeRadau:
    """Two deterministic polynomial steps, not a numerical integration."""
    nfev = njev = nlu = 0

    def __init__(self, rhs, start, initial, end, **kwargs):
        self.t = start
        self.t_bound = end
        self.y = initial.copy()
        self.status = "running"
        self._step = 0
        self.settings = kwargs

    def step(self):
        old = self.t
        self.t = (self._step + 1) * self.t_bound / 2
        Q = np.array([[1., -.1, .02], [-.5, .2, .01], [.3, 0., -.01],
                      [-.2, .04, 0.]]) * (self._step + 1) * 1e-3
        self._dense = RadauDenseOutput(old, self.t, self.y.copy(), Q)
        self.y = self._dense(self.t)
        self._step += 1
        if self._step == 2:
            self.status = "finished"

    def dense_output(self):
        return self._dense


def _fake_run(monkeypatch, observer=None):
    monkeypatch.setattr(runner, "np", np, raising=False)
    monkeypatch.setattr(runner, "Radau", _FakeRadau, raising=False)
    settings = {"rtol": 1e-10, "atol": np.full(4, 1e-14), "max_step": .005,
                "coordinate_scales": [1., 1.], "velocity_scale_multiplier": 1.}
    monkeypatch.setattr(runner, "time_settings", lambda *a: settings)
    config = {"material_geometry": {"h": .1}, "safety": {}}
    return runner.integrate_case(_FakeDisc(), None, None, config, .05, "tight",
                                 np.array([0., .1, .2, .3, .4]),
                                 time.perf_counter() + 1,
                                 initial_coordinates=np.array([.001, -.002]),
                                 dense_output_observer=observer)


def test_dense_observer_leaves_default_sampled_history_unchanged(monkeypatch):
    plain, old = _fake_run(monkeypatch)
    records = []
    observed, new = _fake_run(monkeypatch, records.append)
    np.testing.assert_array_equal(plain, observed)
    assert len(records) == 2
    for key in ("status", "p", "ndof", "nq", "rtol", "atol", "max_step",
                "accepted_internal_steps", "internal_time_steps", "nfev", "njev",
                "nlu", "counters", "time_end", "target_time_end"):
        assert old[key] == new[key]
    assert records[0]["t_old"] == 0 and records[-1]["t"] == .4
    assert records[0]["Q"].shape == (4, 3)


def test_dense_observer_receives_copies_of_native_polynomial_arrays(monkeypatch):
    plain, _ = _fake_run(monkeypatch)
    def damage_only_observer_copies(record):
        record["Q"][:] = np.nan
        record["y_old"][:] = np.inf
    unchanged, _ = _fake_run(monkeypatch, damage_only_observer_copies)
    np.testing.assert_array_equal(plain, unchanged)


def test_actual_radau_dense_polynomial_has_no_extra_step_multiplier():
    q = np.arange(12, dtype=float).reshape(4, 3) / 100
    initial = np.array([.001, -.002, .003, -.004])
    times = np.array([.2, .21, .245, .3])
    local = (times - .2) / .1
    manual = initial[:, None] + q @ np.vstack((local, local**2, local**3))
    native = RadauDenseOutput(.2, .3, initial, q)
    np.testing.assert_allclose(native(times), manual, rtol=0, atol=3e-16)


def test_new_authorization_and_frozen_scope(continuation, config):
    assert continuation.validate_config(config) is config
    assert config["authorization"]["id"] == "explicit_user_FEM3B_2026_10_09"
    assert config["authorization"]["maximum_production_CCX_jobs"] == 2
    assert config["authorization"]["maximum_nonlinear_1D_integrations"] == 2
    assert config["authorization"]["case_order"] == ["linear", "nonlinear"]
    assert config["authorization"]["automatic_retry"] is False
    assert config["parent_completed"]["manifest_sha256"] == PARENT_SHA
    assert config["one_d"]["main_degree"] == 64
    assert config["one_d"]["optional_spatial_degree"] == 48
    assert config["one_d"]["target_T1_fraction"] == .5
    assert config["one_d"]["time_level"] == "tight"
    assert config["admitted"] is False
    assert config["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED"
    assert config["threads"] == 1 and config["job_memory_limit_bytes"] == 4 * 1024**3
    assert config["numerical_budget_seconds"] == 6000
    assert config["job_timeout_seconds_by_horizon"] == {"0.25": 1200, "0.5": 2400}


@pytest.mark.parametrize("key,value", [
    ("schema", "old-continuation"), ("threads", 2), ("numerical_budget_seconds", 6001),
    ("job_memory_limit_bytes", 8 * 1024**3), ("execution_mode", "STRICT_ADMITTED"),
    ("admitted", True), ("new_meshes", True), ("new_modal_jobs", True),
    ("new_static_only_jobs", True), ("new_time_or_space_levels", True)])
def test_forbidden_scope_extension_rejected(continuation, config, key, value):
    config[key] = value
    with pytest.raises(ValueError):
        continuation.validate_config(config)


@pytest.mark.parametrize("key,value", [("id", "explicit_user_FEM3AR_2026_10_09"),
    ("maximum_production_CCX_jobs", 3), ("maximum_nonlinear_1D_integrations", 3),
    ("automatic_retry", True), ("case_order", ["nonlinear", "linear"])])
def test_permission_does_not_expand_or_reuse_old_authorization(continuation, config, key, value):
    config["authorization"][key] = value
    with pytest.raises(ValueError):
        continuation.validate_config(config)


def test_parent_success_manifest_and_old_ledgers_remain_unchanged(config):
    if not (PARENT / "manifest.json").exists():
        pytest.skip("Immutable local FEM-3AR bundle unavailable")
    assert hashlib.sha256((PARENT / "manifest.json").read_bytes()).hexdigest() == PARENT_SHA
    old = json.loads((PARENT / "summary.json").read_text(encoding="utf8"))
    assert old["overall"] == "PILOT_COMPLETE_WITH_QUALIFICATIONS"
    assert old["job_calls"]["CCX_production"] == 2
    assert old["job_calls"]["1D_nonlinear_ODE"] == 1
    assert old["resume_statuses"]["NLSP_FEM3AR_ENERGY_DIAGNOSTICS"] == "PARTIAL"
    assert old["authorization"]["id"] != config["authorization"]["id"]
    failed = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7"
    historical = json.loads((failed / "summary.json").read_text(encoding="utf8"))
    assert historical["overall"] == "BLOCKED_BY_SOLVER"
    assert historical["job_calls"]["CCX_production"] == 1


def test_static_preload_is_not_reported_as_native_dynamic_zero():
    if not (PARENT / "manifest.json").exists():
        pytest.skip("Immutable local FEM-3AR bundle unavailable")
    for kind in ("linear", "nonlinear"):
        with np.load(PARENT / "cases" / kind / "section_history.npz") as z:
            assert z["time"][0] > 0
            assert len(z["time"]) == 102
            assert np.all(np.diff(z["time"]) > 0)


def test_correction_decomposition_preserves_initial_state_and_evolution(continuation):
    # Physical values on a common x grid, including explicitly known t=0.
    linear = np.array([[0., .005, 0.], [0., .004, 0.], [0., .002, 0.]])
    inherited = np.array([[0., -9e-6, 0.], [0., -7e-6, 0.], [0., -4e-6, 0.]])
    dynamic = np.array([[0., 0., 0.], [0., -2e-6, 0.], [0., -6e-6, 0.]])
    common = linear + inherited
    nonlinear = common + dynamic
    d = continuation.correction_decomposition(nonlinear, linear, common)
    np.testing.assert_allclose(d["initial_state_component"], inherited, rtol=0, atol=5e-19)
    np.testing.assert_allclose(d["same_ic_nonlinear_component"], dynamic, rtol=0, atol=5e-19)
    np.testing.assert_allclose(d["total_correction"], inherited + dynamic, rtol=0, atol=5e-19)
    np.testing.assert_array_equal(d["evolving_correction"][0], 0)
    np.testing.assert_allclose(d["evolving_correction"], nonlinear-linear-(nonlinear-linear)[0], rtol=0, atol=0)
    assert d["identity_max_abs"] < 1e-18
    assert d["physical_initial_states_not_amplitude_aligned"] is True


@pytest.mark.parametrize("shape", [(3, 4), (2, 3, 1), (3,)])
def test_decomposition_incompatible_shapes_rejected(continuation, shape):
    with pytest.raises(ValueError):
        continuation.correction_decomposition(np.zeros(shape), np.zeros((2, 4)), np.zeros((2, 4)))


def _candidate(fraction, signal, p_difference=1e-9):
    return {"horizon_T1": fraction, "midspan_evolution_signal": signal,
            "p48_p64_evolution_difference": p_difference,
            "safety_passed": True, "resource_preflight_passed": True}


def test_preferred_horizon_uses_tenfold_prespecified_indicator_before_fem(continuation, config):
    threshold = 10 * 7.78842e-9
    rows = [_candidate(.25, threshold), _candidate(.5, 2*threshold)]
    selected = continuation.choose_horizon(rows, config)
    assert selected["horizon_T1"] == .25
    assert selected["selected_before_any_new_3D_result"] is True
    assert selected["chosen_by_future_1D_3D_agreement"] is False
    assert selected["planning_diagnostics"]["0.25"]["physical_validation_threshold"] is False
    assert selected["planning_diagnostics"]["0.25"]["planning_requirement"] == threshold


def test_preferred_horizon_must_also_exceed_observed_p_sensitivity(continuation, config):
    rows = [_candidate(.25, 1e-7, p_difference=2e-8), _candidate(.5, 4e-7, p_difference=2e-8)]
    selected = continuation.choose_horizon(rows, config)
    assert selected["horizon_T1"] == .5
    assert selected["planning_diagnostics"]["0.25"]["planning_requirement"] == 2e-7
    assert selected["exploratory_signal_resolution_qualification"] is False


def test_both_unresolved_candidates_remain_exploratory(continuation, config):
    selected = continuation.choose_horizon([_candidate(.25, 1e-10), _candidate(.5, 1e-9)], config)
    assert selected["horizon_T1"] == .5
    assert selected["exploratory_signal_resolution_qualification"] is True


def test_missing_p48_is_qualified_not_replaced_by_new_static_solve(continuation, config):
    decision = continuation.choose_horizon([_candidate(.25, 2e-7, None), _candidate(.5, 1e-6, None)], config)
    assert decision["horizon_T1"] == .25
    assert decision["planning_diagnostics"]["0.25"]["p48_p64_available"] is False


@pytest.mark.parametrize("key", ["safety_passed", "resource_preflight_passed"])
def test_horizon_needs_safety_and_resource_preflight(continuation, config, key):
    rows = [_candidate(.25, 2e-7), _candidate(.5, 1e-6)]
    rows[0][key] = False
    with pytest.raises(ValueError):
        continuation.choose_horizon(rows, config)


@pytest.mark.parametrize("signal", [-1., np.nan, np.inf])
def test_invalid_planning_signal_rejected(continuation, config, signal):
    with pytest.raises(ValueError):
        continuation.choose_horizon([_candidate(.25, signal), _candidate(.5, 1e-6)], config)


def test_only_two_predeclared_horizon_rows_are_accepted(continuation, config):
    with pytest.raises(ValueError):
        continuation.choose_horizon([_candidate(.25, 1e-7), _candidate(.25, 2e-7), _candidate(.5, 1e-6)], config)


def test_horizon_cannot_be_chosen_after_any_native_attempt(continuation, config, tmp_path):
    item = {"long_horizon_config": config}
    for calls, attempts in [(1, []), (0, [{"status": "STARTED"}])]:
        summary = {"job_calls": {"CCX_production": calls}, "attempts": attempts}
        with pytest.raises(ValueError, match="after native"):
            continuation.freeze_horizon(tmp_path, item, summary, [])


def test_dense_archive_matches_native_cubic_and_preserves_actual_prefix(continuation, tmp_path):
    y = np.array([.001, -.002])
    q1 = np.array([[.01, .02, -.001], [.03, -.01, .002]])
    step1 = RadauDenseOutput(0., .3, y, q1)
    q2 = q1 * -.7
    step2 = RadauDenseOutput(.3, .7, step1(.3), q2)
    records = [{"t_old": 0., "t": .3, "y_old": y, "Q": q1},
               {"t_old": .3, "t": .7, "y_old": step1(.3), "Q": q2}]
    path = tmp_path / "dense.npz"
    meta = continuation.save_dense_records(path, records)
    times = np.array([0., .17, .3, .4, .7])
    expected = np.vstack([step1(t) if t <= .3 else step2(t) for t in times])
    np.testing.assert_allclose(continuation.evaluate_dense_records(path, times), expected, rtol=0, atol=1e-17)
    assert meta["actual_end"] == .7 and meta["accepted_steps"] == 2
    assert meta["new_integrations"] == 0 and meta["no_trajectory_interpolation"] is True
    for requested in ([-.001, .1], [.1, .700001], [.2, .1], [np.nan]):
        with pytest.raises(ValueError, match="actual Radau prefix"):
            continuation.evaluate_dense_records(path, requested)


def test_dense_archive_rejects_noncontiguous_or_nonfinite_records(continuation, tmp_path):
    basic = {"t_old": 0., "t": .2, "y_old": np.zeros(2), "Q": np.zeros((2, 3))}
    invalid = copy.deepcopy(basic); invalid["Q"][0, 0] = np.nan
    with pytest.raises(ValueError, match="Nonfinite"):
        continuation.save_dense_records(tmp_path / "bad.npz", [invalid])
    second = copy.deepcopy(basic); second.update(t_old=.21, t=.4)
    with pytest.raises(ValueError, match="contiguous"):
        continuation.save_dense_records(tmp_path / "gap.npz", [basic, second])


def test_equal_native_schedules_are_used_without_interpolation(continuation):
    times = np.array([.01, .02, .04])
    linear = np.arange(9, dtype=float).reshape(3, 3)
    nonlinear = linear - 1e-7
    result = continuation.pair_native_histories(times, linear, times.copy(), nonlinear)
    np.testing.assert_array_equal(result["time"], times)
    np.testing.assert_array_equal(result["correction"], nonlinear-linear)
    assert result["sample_origin"] == "native" and result["interpolation_used"] is False
    assert result["interpolation_correction_difference_max_abs"] == 0


def test_divergent_schedules_need_explicit_common_grid(continuation):
    with pytest.raises(ValueError, match="predeclared explicit"):
        continuation.pair_native_histories([.01, .03], np.zeros((2, 3)), [.02, .04], np.zeros((2, 3)))


def test_interpolation_keeps_linear_and_pchip_difference_visible(continuation):
    tl = np.array([.01, .1, .3, .5])
    tn = np.array([.02, .13, .31, .52])
    grid = np.array([.05, .15, .25, .45])
    linear = (tl**2)[:, None]
    nonlinear = (tn**2 + 1e-6*tn**3)[:, None]
    result = continuation.pair_native_histories(tl, linear, tn, nonlinear, common_grid=grid)
    assert result["sample_origin"] == "interpolated"
    assert result["interpolation_used"] is True and result["native_values_not_relabelled"] is True
    assert result["primary_interpolation"] == "piecewise_linear"
    assert result["diagnostic_interpolation"] == "shape_preserving_cubic_PCHIP"
    assert result["no_extrapolation"] is True
    np.testing.assert_array_equal(result["correction"], result["nonlinear"]-result["linear"])
    difference = np.max(abs(result["correction_PCHIP"]-result["correction"]))
    assert result["interpolation_correction_difference_max_abs"] == difference
    assert difference > 1e-8


@pytest.mark.parametrize("grid", [[0., .1], [.1, .6], [.2, .1], [.1]])
def test_interpolation_cannot_extrapolate_or_invent_missing_prefix(continuation, grid):
    with pytest.raises(ValueError, match="actual native overlap"):
        continuation.pair_native_histories([.01, .3], np.zeros((2, 1)), [.02, .4], np.zeros((2, 1)), common_grid=grid)


@pytest.mark.parametrize("bad", [np.nan, np.inf])
def test_native_pairing_rejects_nonfinite_time_and_grid(continuation, bad):
    with pytest.raises(ValueError):
        continuation.pair_native_histories([.01, bad], np.zeros((2, 1)), [.01, bad], np.zeros((2, 1)))
    with pytest.raises(ValueError):
        continuation.pair_native_histories([.01, .3], np.zeros((2, 1)), [.02, .4], np.zeros((2, 1)), common_grid=[.1, bad])


def _forbid_science(monkeypatch, module):
    def forbidden(*args, **kwargs):
        pytest.fail("Cached/synthetic test attempted a scientific solver call")
    monkeypatch.setattr(module.base.one, "integrate_nonlinear_reference", forbidden)
    monkeypatch.setattr(module.base.one, "exact_linear_reference", forbidden)
    monkeypatch.setattr(module.base.one.runner, "integrate_case", forbidden)
    monkeypatch.setattr(module.base.base.fem1, "run_job", forbidden)
    monkeypatch.setattr(module.base.base.fem1.single, "generate_mesh_with_gmsh_python", forbidden)
    monkeypatch.setattr(module.base.base.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)


def _summary():
    return {"overall": "NOT_RUN", "cases": {}, "attempts": [],
            "job_calls": {"CCX_production": 0}, "numerical_seconds": 0.,
            "long_horizon_statuses": {}}


def _frozen(continuation, tmp_path, config):
    continuation.write_json(tmp_path / "horizon_decision.json", {"horizon_T1": .25})
    science = {"job_timeout_seconds": 1200}
    continuation.write_json(tmp_path / "frozen_3d_config.json", science)
    summary = _summary()
    summary["horizon_decision_sha256"] = continuation.sha(tmp_path / "horizon_decision.json")
    return {"long_horizon_config": config}, summary


def test_freeze_sets_horizon_timeout_and_preserves_all_dynamic_settings(continuation, config, tmp_path, monkeypatch):
    original = {"dynamic": {"initial_increment_T1_fraction": 1/4000,
                           "max_increment_T1_fraction": 1/2000, "ALPHA": 0},
                "geometry": {"L": 1., "b": .2, "h": .1},
                "job_timeout_seconds": 1200, "horizon_T1": .05}
    item = {"config": original, "long_horizon_config": config, "authorization": config["authorization"]}
    summary = _summary(); summary["preflight"] = {"T1": 10.37828159055014}
    monkeypatch.setattr(continuation, "save", lambda *a: None)
    rows = [_candidate(.25, 1e-10), _candidate(.5, 1e-6)]
    decision = continuation.freeze_horizon(tmp_path, item, summary, rows)
    saved = continuation.read_json(tmp_path / "frozen_3d_config.json")
    assert decision["horizon_T1"] == .5 and saved["job_timeout_seconds"] == 2400
    assert decision["target_end"] == .5 * summary["preflight"]["T1"]
    assert saved["dynamic"] == original["dynamic"]
    assert saved["geometry"] == original["geometry"]
    assert original["horizon_T1"] == .05 and original["job_timeout_seconds"] == 1200
    with pytest.raises(ValueError, match="already frozen"):
        continuation.freeze_horizon(tmp_path, item, summary, rows)


def test_no_native_job_before_frozen_horizon(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    with pytest.raises(ValueError, match="Freeze"):
        continuation.run_3d_case(tmp_path, {"long_horizon_config": config}, _summary(), "linear")


def test_horizon_hash_cannot_change_before_native_job(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    item, summary = _frozen(continuation, tmp_path, config)
    continuation.write_json(tmp_path / "horizon_decision.json", {"horizon_T1": .5})
    with pytest.raises(ValueError, match="decision changed"):
        continuation.run_3d_case(tmp_path, item, summary, "linear")


@pytest.mark.parametrize("linear", [{}, {"status": "FAIL"}, {"status": "PASS"},
    {"status": "PASS", "continuation_audit": {"status": "FAIL"}}])
def test_nonlinear_requires_actual_linear_execution_and_audit(continuation, config, tmp_path, monkeypatch, linear):
    _forbid_science(monkeypatch, continuation)
    item, summary = _frozen(continuation, tmp_path, config)
    if linear:
        summary["cases"]["linear"] = linear
    with pytest.raises(ValueError, match="Actual linear"):
        continuation.run_3d_case(tmp_path, item, summary, "nonlinear")


@pytest.mark.parametrize("status", ["FAIL", "STARTED"])
def test_prior_attempt_not_retried_even_with_new_code_hash(continuation, config, tmp_path, monkeypatch, status):
    _forbid_science(monkeypatch, continuation)
    item, summary = _frozen(continuation, tmp_path, config)
    item["helper_sha256"] = {"updated-helper": "new-code-hash"}
    summary["cases"]["linear"] = {"status": status}
    assert continuation.run_3d_case(tmp_path, item, summary, "linear") is False
    assert summary["job_calls"]["CCX_production"] == 0


def test_successful_case_replays_without_native_job_or_recovery(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    item, summary = _frozen(continuation, tmp_path, config)
    summary["cases"]["linear"] = {"status": "PASS"}
    monkeypatch.setattr(continuation.base, "verify_sources", lambda *a: pytest.fail("Successful case should replay"))
    assert continuation.run_3d_case(tmp_path, item, summary, "linear") is True


def test_native_failure_stops_once_and_preserves_attempt_evidence(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    item, summary = _frozen(continuation, tmp_path, config)
    monkeypatch.setattr(continuation, "save", lambda *a: None)
    monkeypatch.setattr(continuation.base, "verify_sources", lambda *a: ({"science_config": {"ccx_exe": "fake-ccx"}}, {}, {}, {}))
    def input_only(path, *a):
        path.write_text("synthetic unchanged input", encoding="utf8")
        return {"no_actual_mesh_or_deck_generation": True}
    monkeypatch.setattr(continuation.base, "write_input", input_only)
    monkeypatch.setattr(continuation.resume, "output_safety", lambda *a: {"status": "PASS"})
    calls = []
    def fake_failure(command, case, timeout, memory, prefix, env):
        calls.append(command)
        assert timeout == 1200 and memory == 4 * 1024**3
        assert all(env[k] == "1" for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS"))
        (case / "motion.stdout.txt").write_text("*ERROR input rejected", encoding="utf8")
        return SimpleNamespace(returncode=1), {"failure": None, "returncode": 1}
    monkeypatch.setattr(continuation.base.base.fem1, "run_job", fake_failure)
    assert continuation.run_3d_case(tmp_path, item, summary, "linear") is False
    assert calls == [["fake-ccx", "motion"]]
    assert summary["job_calls"]["CCX_production"] == 1
    assert summary["hard_stop"] is True and summary["overall"] == "BLOCKED_BY_SOLVER"
    assert summary["attempts"][0]["authorization_id"] == "explicit_user_FEM3B_2026_10_09"
    assert summary["attempts"][0]["status"] == "FAIL"
    assert (tmp_path / "cases/linear/motion.inp").exists()
    assert (tmp_path / "cases/linear/motion.stdout.txt").exists()
    assert continuation.run_3d_case(tmp_path, item, summary, "nonlinear") is False
    assert continuation.run_3d_case(tmp_path, item, summary, "linear") is False
    assert len(calls) == 1


def test_cache_finds_authorization_attempt_despite_new_helper_hash(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    bundle = tmp_path / "old-fingerprint"; bundle.mkdir()
    identity = {"authorization": config["authorization"], "long_horizon_config": config, "helper_sha256": {"helper": "oldhash"}}
    continuation.write_json(bundle / "provenance.json", identity)
    continuation.write_json(bundle / "manifest.json", {})
    summary = {"overall": "BLOCKED_BY_SOLVER"}
    monkeypatch.setattr(continuation, "OUTPUT", tmp_path)
    monkeypatch.setattr(continuation, "validate_cache", lambda *a: summary)
    found = continuation.existing_attempt(config)
    assert found[0] == bundle and found[1] == identity and found[2] is summary


def test_unmanifested_attempt_cannot_be_recreated_as_hidden_retry(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    bundle = tmp_path / "interrupted"; bundle.mkdir()
    continuation.write_json(bundle / "provenance.json", {"authorization": config["authorization"], "long_horizon_config": config})
    monkeypatch.setattr(continuation, "OUTPUT", tmp_path)
    with pytest.raises(RuntimeError, match="no hidden retry"):
        continuation.existing_attempt(config)


def test_matching_completed_compute_runs_zero_science(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary["overall"] = "PILOT_COMPLETE_WITH_QUALIFICATIONS"; summary["completed"] = True
    monkeypatch.setattr(continuation, "prepare_stage", lambda *a: (tmp_path, {}, summary))
    monkeypatch.setattr(continuation, "run_3d_case", lambda *a: pytest.fail("Completed compute attempted job routing"))
    assert continuation.main(["--run-3d"]) is summary


def test_report_only_performs_zero_science(continuation, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary["overall"] = "PARTIAL"
    monkeypatch.setattr(continuation, "validate_cache", lambda *a: summary)
    monkeypatch.setattr(continuation, "prepare_stage", lambda *a: pytest.fail("Report attempted preparation"))
    assert continuation.main(["--report-only", str(tmp_path)]) is summary


def test_old_diagnostic_must_precede_any_new_one_d_integration(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary["long_horizon_statuses"]["NLSP_FEM3B_OLD_SIGNAL_DIAGNOSTIC"] = "NOT_RUN"
    monkeypatch.setattr(continuation, "load_reference", lambda *a: pytest.fail("1D starts before old-data gate"))
    with pytest.raises(ValueError, match="old-signal"):
        continuation.run_1d_stage(tmp_path, {"long_horizon_config": config, "config": {}}, summary)


def test_cached_one_d_candidate_evidence_performs_zero_scientific_calls(continuation, config, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary["long_horizon_statuses"]["NLSP_FEM3B_OLD_SIGNAL_DIAGNOSTIC"] = "PASS"
    continuation.write_json(tmp_path / "candidate_signals.json", [{"frozen": True}])
    monkeypatch.setattr(continuation, "load_reference", lambda *a: pytest.fail("Matching 1D cache rebuilt reference/modes"))
    assert continuation.run_1d_stage(tmp_path, {"long_horizon_config": config, "config": {}}, summary) is summary


@pytest.mark.parametrize("status", ["STARTED", "PARTIAL", "FAIL"])
def test_interrupted_one_d_attempt_stops_before_reference_or_eigensolve(continuation, config, tmp_path, monkeypatch, status):
    _forbid_science(monkeypatch, continuation)
    summary = _summary()
    summary["long_horizon_statuses"]["NLSP_FEM3B_OLD_SIGNAL_DIAGNOSTIC"] = "PASS"
    summary["one_d_attempts"] = [{"p": 64, "status": status}]
    monkeypatch.setattr(continuation, "load_reference", lambda *a: pytest.fail("Interrupted attempt constructed another reference/modes"))
    monkeypatch.setattr(continuation, "save", lambda *a: None)
    assert continuation.run_1d_stage(tmp_path, {"long_horizon_config": config, "config": {}}, summary) is summary
    assert summary["hard_stop"] is True


def test_matching_completed_full_compute_skips_all_stage_routing(continuation, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary["completed"] = True
    monkeypatch.setattr(continuation, "prepare_stage", lambda *a: (tmp_path, {}, summary))
    monkeypatch.setattr(continuation, "run_1d_stage", lambda *a: pytest.fail("Cached compute reran 1D planning"))
    monkeypatch.setattr(continuation, "run_3d_case", lambda *a: pytest.fail("Cached compute reran CCX"))
    assert continuation.main(["--compute"]) is summary


def test_long_horizon_dispatch_does_not_change_old_resume_mode(continuation, monkeypatch):
    calls = []
    monkeypatch.setattr(continuation, "main", lambda args: calls.append(args) or {"cached": True})
    result = continuation.resume.main(["--long-horizon", "--report-only", "synthetic-bundle"])
    assert result == {"cached": True}
    assert calls == [["--report-only", "synthetic-bundle"]]


@pytest.fixture
def diagnostics():
    return importlib.import_module("scripts.lib.nlsp_fem3b_diagnostics")


def test_sampled_evolution_norms_keep_sign_location_and_fixed_domain(diagnostics):
    x = np.array([0., .5, 1.]); times = np.array([.1, .2, .3])
    values = np.array([[0., -1e-6, 0.], [0., 2e-6, 0.], [0., -1.5e-6, 0.]])
    metrics, history = diagnostics.sampled_norms(values, x, times)
    assert metrics["absolute_max"] == 2e-6
    assert metrics["signed_at_max"] == 2e-6 and metrics["x_at_max"] == .5
    assert metrics["time_at_max"] == .2
    assert metrics["midspan_final"] == -1.5e-6
    assert metrics["max_time_L2"] == pytest.approx(np.sqrt(.5)*2e-6)
    np.testing.assert_array_equal(history["cumulative_max"], [1e-6, 2e-6, 2e-6])
    assert metrics["sampled_maxima_only"] is True


def test_evolving_correction_uses_confirmed_static_initial_profile(diagnostics):
    initial = np.array([0., -2e-6, 0.])
    positive_time_correction = np.array([[0., -2.5e-6, 0.], [0., -3e-6, 0.]])
    evolved = diagnostics.evolving_correction(positive_time_correction, initial)
    np.testing.assert_allclose(evolved[:, 1], [-.5e-6, -1e-6], rtol=0, atol=5e-22)
    assert evolved[0, 1] != 0  # The first dynamic sample is not a fabricated t=0.


def test_signal_resolution_is_observed_scale_diagnostic_not_validation(diagnostics):
    result = diagnostics.signal_resolution(2e-6, {"output": 5e-9, "recovery": 1e-8, "interpolation": 0., "time": None})
    assert result["signal_to_largest_observed_difference"] == 200
    assert result["status"] == "SIGNAL_EXCEEDS_OBSERVED_DIAGNOSTIC_DIFFERENCES"
    assert result["physical_validation_threshold"] is False
    assert result["continuum_error_bound_claimed"] is False
    assert result["independent_3D_temporal_spatial_certification"] is False
    assert "time" not in result["observed_diagnostic_differences"]
    unresolved = diagnostics.signal_resolution(5e-9, {"output": 5e-9, "recovery": 1e-8})
    assert unresolved["status"] == "DYNAMIC_NONLINEAR_EVOLUTION_UNRESOLVED"
    zero = diagnostics.signal_resolution(0., {"output": 0.})
    assert zero["signal_to_largest_observed_difference"] is None
    assert zero["status"] == "DYNAMIC_NONLINEAR_EVOLUTION_UNRESOLVED"


@pytest.mark.parametrize("signal,differences", [(-1., {"d": 0.}), (np.nan, {}), (1., {"d": -1.}), (1., {"d": np.inf})])
def test_invalid_signal_resolution_inputs_rejected(diagnostics, signal, differences):
    with pytest.raises(ValueError):
        diagnostics.signal_resolution(signal, differences)


def test_one_d_qualifications_measure_all8_and_additive_components_on_fixed_scale(diagnostics, tmp_path):
    times = np.array([0., .1, .2]); x = np.array([0., .5, 1.])
    linear = np.zeros((3, 3, 4)); linear[:, 1, :] = np.array([1., .9, .7])[:, None] * 1e-3
    initial = np.zeros_like(linear); initial[:, 1, :] = np.array([-9., -7., -4.])[:, None] * 1e-6
    same = np.zeros_like(linear); same[:, 1, :] = np.array([0., -2., -6.])[:, None] * 1e-6
    total = initial + same
    np.savez(tmp_path / "one_d_p64_decomposition.npz", times=times, x=x,
             total_correction=total, initial_state_component=initial,
             same_ic_nonlinear_component=same, evolving_correction=total-total[0])
    for p in (48, 64):
        for kind, values in (("linear", linear), ("nonlinear", linear+total)):
            field = values + (0 if p == 64 else 1e-8)
            np.savez(tmp_path / ("one_d_p"+str(p)+"_"+kind+".npz"), times=times, x=x,
                     fields=field, physical_velocities=np.zeros_like(field)+(0 if p == 64 else 1e-7),
                     energy=np.array([1., 1.+1e-12, 1.+2e-12]), energy_relative_drift=np.array([0., 1e-12, 2e-12]))
    result = diagnostics.one_d_qualifications(tmp_path, .2)
    assert result["decomposition"]["identity_max_abs"] < 1e-20
    control = result["spatial_degree_control"]
    assert control["status"] == "COMPLETED_DIAGNOSTIC"
    assert len(control["fields"]) == 8
    assert set(control["correction"]) == {"u", "w", "theta", "c_eff_diagnostic"}
    assert control["fields"]["w"]["absolute_max"] == pytest.approx(1e-8)
    assert control["fields"]["w_t"]["absolute_max"] == 1e-7
    assert control["evolution"]["w"]["absolute_max"] < 1e-18
    row = result["decomposition"]["fields"]["w"]
    assert row["total_correction"]["midspan_final"] == pytest.approx(-1e-5, rel=0, abs=1e-20)
    assert row["initial_state_component"]["midspan_final"] == -4e-6
    assert row["same_ic_nonlinear_component"]["midspan_final"] == -6e-6
    assert result["temporal_control"]["status"] == "NOT_CERTIFIED"


def test_plot_only_uses_saved_arrays_and_selected_horizon_zero_science(continuation, diagnostics, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary.update(preflight={"T1": 1.}, selected_horizon_T1=.25)
    x = np.linspace(0., 1., 41); times = np.array([.1, .25])
    one = np.zeros((2, 41, 4)); three = np.zeros((2, 41, 7))
    np.savez(tmp_path / "dynamic_comparison.npz", time=times, x=x,
             linear_1D=one, nonlinear_1D=one, linear_3D=three, nonlinear_3D=three,
             one_d_initial_correction=np.zeros((41, 4)), three_d_initial_correction=np.zeros((41, 7)),
             one_d_correction=one, three_d_correction=three, one_d_evolution=one, three_d_evolution=three)
    np.savez(tmp_path / "one_d_p64_decomposition.npz", times=np.r_[0., times], x=x,
             total_correction=np.zeros((3, 41, 4)), initial_state_component=np.zeros((3, 41, 4)),
             same_ic_nonlinear_component=np.zeros((3, 41, 4)))
    continuation.write_json(tmp_path / "dynamic_comparison.json", {"pairing": {"interpolation_used": False}})
    continuation.write_json(tmp_path / "provenance.json", {"frozen": True})
    calls = []
    monkeypatch.setattr(continuation, "validate_cache", lambda *a: calls.append("validate") or summary)
    monkeypatch.setattr(continuation, "save", lambda *a: calls.append("save"))
    from matplotlib.figure import Figure
    def inspect_figure(self, path, **kwargs):
        assert all(tuple(ax.get_xlim()) == (0., .25) for ax in self.axes if ax.get_xlabel() == "t/T1")
        calls.append(Path(path).suffix)
    monkeypatch.setattr(Figure, "savefig", inspect_figure)
    result = diagnostics.plot_bundle(tmp_path)
    assert result == {"figures": 3, "new_scientific_calls": 0}
    assert calls.count(".pdf") == calls.count(".png") == 3
    assert calls[0] == "validate" and calls[-1] == "save"


def test_saved_new_one_d_runs_reuse_each_frozen_static_state_exactly():
    if not (NEW_BUNDLE / "candidate_signals.json").exists():
        pytest.skip("Completed local FEM-3B Stage B unavailable")
    identity = json.loads((NEW_BUNDLE / "provenance.json").read_text(encoding="utf8"))
    source = ROOT / identity["config"]["source_static"]["bundle"]
    previous_times = None
    for p in (64, 48):
        with np.load(source / ("one_d_p"+str(p)+".npz")) as static:
            for kind in ("linear", "nonlinear"):
                with np.load(NEW_BUNDLE / ("one_d_p"+str(p)+"_"+kind+".npz")) as data:
                    assert data["q"].shape[1] == data["velocity"].shape[1] == 4*(p-1)
                    np.testing.assert_array_equal(data["q"][0], static["q_"+kind])
                    np.testing.assert_array_equal(data["velocity"][0], np.zeros(4*(p-1)))
                    assert data["times"][0] == 0.
                    assert data["times"][-1] == .5 * 10.37828159055014
                    assert .25 * 10.37828159055014 in data["times"]
                    if previous_times is not None:
                        np.testing.assert_array_equal(data["times"], previous_times)
                    previous_times = data["times"].copy()
        stats = json.loads((NEW_BUNDLE / ("one_d_p"+str(p)+"_nonlinear.json")).read_text(encoding="utf8"))["execution"]
        assert stats["status"] == "PASS" and stats["new_ODE_integrations"] == 1
        assert stats["initial_coordinates_reused_exactly"] is True
        assert stats["no_dynamic_derivative_constraints"] is True
        assert stats["rtol"] == 1e-10 and len(stats["atol"]) == 8*(p-1)
        assert stats["max_step"] == .0072093929877006385
        assert stats["authorization_id"] == "explicit_user_FEM3B_2026_10_09"


def test_actual_stage_b_horizon_and_qualifications_not_rewritten_as_all8_pass():
    if not (NEW_BUNDLE / "candidate_signals.json").exists():
        pytest.skip("Completed local FEM-3B Stage B unavailable")
    choice = json.loads((NEW_BUNDLE / "horizon_decision.json").read_text(encoding="utf8"))
    assert choice["horizon_T1"] == .25 and choice["job_timeout_seconds"] == 1200
    assert choice["created_before_native_attempts"] == 0
    assert choice["selected_before_any_new_3D_result"] is True
    assert choice["chosen_by_future_1D_3D_agreement"] is False
    short = choice["planning_diagnostics"]["0.25"]
    assert short["signal_exceeds_planning_requirement"] is True
    assert short["physical_validation_threshold"] is False
    all8 = json.loads((NEW_BUNDLE / "one_d_all8_spatial.json").read_text(encoding="utf8"))
    assert all8["status"] == "PARTIAL"
    assert len(all8["fields"]) == 8
    assert sum(row["pass"] for row in all8["fields"].values()) == 2
    assert all8["fields"]["velocity_theta"]["tolerance"] == 1e-4
    assert all8["fields"]["velocity_c"]["tolerance"] == 1e-3


def test_saved_linear_native_job_has_completed_preload_release_and_real_frames():
    case = NEW_BUNDLE / "cases/linear"
    if not (case / "continuation_audit.json").exists():
        pytest.skip("Completed local FEM-3B linear native audit unavailable")
    job = json.loads((case / "job.json").read_text(encoding="utf8"))
    audit = json.loads((case / "continuation_audit.json").read_text(encoding="utf8"))
    result = json.loads((case / "recovery.json").read_text(encoding="utf8"))
    assert job["returncode"] == 0 and job["failure"] is None
    assert job["peak_working_set_bytes"] < 4*1024**3
    assert audit["status"] == "PASS" and audit["warning_lines"] == []
    assert audit["preload_force_imbalance_relative"] < audit["equilibrium_gate_unchanged"] == 1e-5
    assert audit["preload_moment_imbalance_relative"] < 1e-5
    assert result["status"] == "PASS"
    assert result["dynamic_output_frames"] == result["dynamic_increments"] == 502
    assert result["dynamic_time_start"] > 0
    assert result["dynamic_time_end"] == .25*10.37828159055014
    assert result["cutbacks"] == 0
    assert result["maximum_native_external_work_after_release"] == 0
    assert result["maximum_native_damping_work_after_release"] == 0
    assert result["energy_status"] == "PARTIAL"
    transfer = result["preload_transfer"]
    assert transfer["status"] == "PASS"
    for key in ("node_displacement_max_difference", "support_RF_max_difference",
                "section_profile_max_difference", "stress_difference", "strain_difference"):
        assert transfer[key] == 0
    assert result["final_strain_diagnostics"]["minimum_det_deformation_gradient"] > 0


def test_fresh_compute_finalizes_and_reports_actual_call_deltas_then_replays_cache(continuation, diagnostics, tmp_path, monkeypatch, capsys):
    """Every numeric stage is a fixture; no real native/ODE/eigen job occurs."""
    _forbid_science(monkeypatch, continuation)
    summary = _summary(); summary["job_calls"]["1D_nonlinear_ODE"] = 0
    item = {"frozen": True}
    order = []
    monkeypatch.setattr(continuation, "prepare_stage", lambda *a: (tmp_path, item, summary))
    def fake_one_d(*args):
        order.append("1D")
        summary["job_calls"]["1D_nonlinear_ODE"] = 2
    def fake_native(bundle, identity, state, kind):
        order.append(kind)
        if kind == "nonlinear":
            assert state["cases"]["linear"]["continuation_audit"]["status"] == "PASS"
        state["cases"][kind] = {"status": "PASS", "continuation_audit": {"status": "PASS"}}
        state["job_calls"]["CCX_production"] += 1
        return True
    monkeypatch.setattr(continuation, "run_1d_stage", fake_one_d)
    monkeypatch.setattr(continuation, "run_3d_case", fake_native)
    monkeypatch.setattr(diagnostics, "actual_pairing_times", lambda *a: order.append("times") or np.array([0., .25]))
    monkeypatch.setattr(continuation, "evaluate_saved_one_d", lambda *a: order.append("saved-state evaluation") or {"times": np.array([0., .25])})
    monkeypatch.setattr(diagnostics, "complete_comparison", lambda *a: order.append("comparison"))
    def save(*args):
        assert summary["completed"] is True
        assert summary["overall"] == "FEM3B_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
        order.append("save")
    def plot(*args):
        assert summary["completed"] is True
        assert (tmp_path / "comparison_one_d.npz").exists()
        order.append("plot")
    monkeypatch.setattr(continuation, "save", save)
    monkeypatch.setattr(diagnostics, "plot_bundle", plot)
    assert continuation.main(["--compute"]) is summary
    assert order == ["1D", "linear", "nonlinear", "times", "saved-state evaluation", "comparison", "save", "plot", "save"]
    stdout = json.loads(capsys.readouterr().out)
    assert stdout["new_scientific_calls"] == {"CCX_production": 2, "1D_nonlinear_ODE": 2}
    order.clear()
    assert continuation.main(["--compute"]) is summary
    assert order == []
    assert json.loads(capsys.readouterr().out)["new_scientific_calls"] == 0


def test_read_only_reference_completion_reuses_saved_modes_after_successful_nonlinear_run(continuation, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    reference = {"disc": SimpleNamespace(p=64, _linear_modes={})}
    times = np.array([0., .1]); x = np.array([0., .5, 1.])
    values = np.zeros((2, 3, 4)); values[:, 1, 1] = [.005, .004]
    nonlinear = values.copy(); nonlinear[:, 1, 1] -= [9e-6, 8e-6]
    np.savez(tmp_path / "linear_modes_p64.npz", omega=np.array([1., 2.]), vectors=np.eye(2))
    np.savez(tmp_path / "one_d_p64_nonlinear.npz", fields=nonlinear, x=x)
    calls = []
    def fixture_linear(bundle, name, ref, requested, initial_kind):
        assert None in ref["disc"]._linear_modes
        np.testing.assert_array_equal(ref["disc"]._linear_modes[None]["omega"], [1., 2.])
        calls.append(initial_kind)
        field = values.copy()
        if initial_kind == "nonlinear":
            field[:, 1, 1] -= [9e-6, 7e-6]
        np.savez(bundle / (name+".npz"), fields=field, x=x)
        continuation.write_json(bundle / (name+".json"), {"source": "synthetic saved factor evaluation"})
    monkeypatch.setattr(continuation, "_linear_history", fixture_linear)
    continuation._complete_degree_references(tmp_path, "one_d_p64", reference, times)
    assert calls == ["linear", "nonlinear"]
    proof = continuation.read_json(tmp_path / "one_d_p64_decomposition.json")
    assert proof["saved_complete_linear_factors_reused"] is True
    assert proof["no_additional_nonlinear_integration"] is True
    calls.clear()
    continuation._complete_degree_references(tmp_path, "one_d_p64", reference, times)
    assert calls == []


def test_missing_saved_linear_factors_do_not_trigger_eigenanalysis_retry(continuation, tmp_path, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    reference = {"disc": SimpleNamespace(p=64, _linear_modes={})}
    with pytest.raises(ValueError, match="no eigenanalysis retry"):
        continuation._complete_degree_references(tmp_path, "one_d_p64", reference, np.array([0., .1]))
