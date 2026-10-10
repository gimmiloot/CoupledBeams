"""FEM-3C orchestration contracts using fixtures and saved data, never jobs."""
from __future__ import annotations

import copy
import hashlib
import importlib
import json
import os
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "results/nlsp_nonlinear_dynamic_long_horizon/7d2b499e6a1eb990"
PARENT_SHA = "e4bcc291fab5f04a2fb73103c35a2fa6481aecf5f058b8ce6ae4a5661efa0c4e"
OLD_FAILURE = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7"
OLD_FAILURE_SHA = "bf217d9c77d42170d1021b55c61873f8b2a94ed1fa67950d5ff9ed82c5cd257f"
NEW_BUNDLE = ROOT / "results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef"


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


@pytest.fixture
def continuation():
    return importlib.import_module("scripts.lib.nlsp_fem3c_validation")


@pytest.fixture
def config(continuation):
    return copy.deepcopy(_read(continuation.CONFIG))


@pytest.fixture(scope="module")
def science():
    return _read(PARENT / "provenance.json")["config"]


def _forbid_science(monkeypatch, module):
    def forbidden(*args, **kwargs):
        pytest.fail("FEM-3C tests must not run FEM, meshes, nonlinear ODE, static, eigen or symbolic work")
    base = module.base
    monkeypatch.setattr(base, "run_case", forbidden)
    monkeypatch.setattr(base, "finish_references_and_comparison", forbidden)
    monkeypatch.setattr(base.one, "integrate_nonlinear_reference", forbidden)
    monkeypatch.setattr(base.one, "exact_linear_reference", forbidden)
    monkeypatch.setattr(base.one.runner, "integrate_case", forbidden)
    monkeypatch.setattr(base.one.dynamics.PlanarGalerkin, "linear_eigenpairs", forbidden)
    monkeypatch.setattr(base.one.rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(base.base, "fem2_static_newton", forbidden)
    monkeypatch.setattr(base.base.fem1, "run_job", forbidden)
    monkeypatch.setattr(base.base.fem1.single, "generate_mesh_with_gmsh_python", forbidden)
    monkeypatch.setattr(base.base.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)


def test_historical_success_and_failure_manifests_are_unchanged():
    assert _sha(PARENT / "manifest.json") == PARENT_SHA
    assert _sha(OLD_FAILURE / "manifest.json") == OLD_FAILURE_SHA
    previous = _read(PARENT / "summary.json")
    assert previous["overall"] == "FEM3B_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
    assert previous["job_calls"]["CCX_production"] == 2
    assert previous["job_calls"]["1D_nonlinear_ODE"] == 2
    failed = _read(OLD_FAILURE / "summary.json")
    assert failed["overall"] == "BLOCKED_BY_SOLVER"
    assert failed["cases"]["linear"]["status"] == "FAIL"


def test_frozen_physics_and_original_generator_hashes_remain_unchanged():
    expected = {
        "scripts/analysis/pilot_nlsp_nonlinear_dynamic_3d_fem.py":
            "8f1f7c8782543546b317cd3ba7cac30162c7fa9b136a101a9e198fdc1f1d2a2c",
        "scripts/lib/weakly_nonlinear_planar_dynamics.py":
            "eea98b77babcb4ed840e325efb270f37a68c32cff1d20ededf43c9989defc548",
    }
    for path, value in expected.items():
        assert _sha(ROOT / path) == value


def test_saved_static_initial_states_have_all_independent_coordinates(science):
    source = ROOT / science["source_static"]["bundle"]
    for p in (48, 64):
        with np.load(source / f"one_d_p{p}.npz") as static:
            for kind in ("linear", "nonlinear"):
                state = static["q_" + kind]
                assert state.shape == (4 * (p - 1),)
                assert np.isfinite(state).all()
                with np.load(PARENT / f"one_d_p{p}_{kind}.npz") as prior:
                    np.testing.assert_array_equal(prior["q"][0], state)
                    np.testing.assert_array_equal(prior["velocity"][0], np.zeros_like(state))
    # Own static equilibria remain distinct; no amplitude matching is allowed.
    with np.load(source / "one_d_p64.npz") as static:
        assert not np.array_equal(static["q_linear"], static["q_nonlinear"])


def test_saved_linear_factors_restored_before_evaluation_without_eigenanalysis(tmp_path, monkeypatch):
    one = importlib.import_module("scripts.lib.nlsp_fem3c_1d")
    old = {"long_horizon_config": {"saved": True}}
    (tmp_path / "provenance.json").write_text(json.dumps(old), encoding="utf8")
    np.savez(tmp_path / "linear_modes_p64.npz", eigenvalues=np.array([1., 4.]),
             omega=np.array([1., 2.]), frequency_hz=np.array([1., 2.])/(2*np.pi),
             vectors=np.eye(2))
    disc = SimpleNamespace(ndof=2, _linear_modes={}, linear_eigendecompositions=0)
    reference = {"disc": disc}
    monkeypatch.setattr(one.previous, "load_reference", lambda cfg, p: reference)
    monkeypatch.setattr(one.previous.base.one.dynamics.PlanarGalerkin, "linear_eigenpairs",
                        lambda *args, **kwargs: pytest.fail("Saved-factor replay starts eigenanalysis"))
    item = {"validation_config": {"parent_completed": {"bundle": str(tmp_path)}}}
    restored, parent, factors = one._reference(item, 64)
    assert restored is reference and parent == tmp_path
    assert disc.linear_eigendecompositions == 0
    assert disc._linear_modes[None] is factors
    np.testing.assert_array_equal(factors["vectors"], np.eye(2))
    assert factors["vectors"].shape == (disc.ndof, disc.ndof)


def test_missing_saved_factors_never_falls_back_to_eigenanalysis(tmp_path, monkeypatch):
    one = importlib.import_module("scripts.lib.nlsp_fem3c_1d")
    (tmp_path / "provenance.json").write_text('{"long_horizon_config": {}}', encoding="utf8")
    reference = {"disc": SimpleNamespace(ndof=2, _linear_modes={}, linear_eigendecompositions=0)}
    monkeypatch.setattr(one.previous, "load_reference", lambda *args: reference)
    monkeypatch.setattr(one.previous.base.one.dynamics.PlanarGalerkin, "linear_eigenpairs",
                        lambda *args, **kwargs: pytest.fail("Missing factors trigger forbidden eigen retry"))
    item = {"validation_config": {"parent_completed": {"bundle": str(tmp_path)}}}
    with pytest.raises(FileNotFoundError):
        one._reference(item, 64)


def test_auxiliary_linear_solution_keeps_own_q0_and_zero_velocities(tmp_path, monkeypatch):
    one = importlib.import_module("scripts.lib.nlsp_fem3c_1d")
    initial = np.array([1e-3, -2e-3])
    calls = []
    def exact(q0, v0, times):
        np.testing.assert_array_equal(q0, initial)
        np.testing.assert_array_equal(v0, np.zeros(2))
        calls.append(times.copy())
        # The initial factor reconstruction has harmless roundoff; q0 must
        # still be stored exactly from the immutable source.
        return {"q": np.broadcast_to(q0+1e-17, (len(times), 2)).copy(),
                "velocity": np.ones((len(times), 2))*1e-17}
    disc = SimpleNamespace(ndof=2, linear_reference=exact)
    reference = {"disc": disc, "saved": {"q_nonlinear": initial}}
    path = tmp_path / "scratch.npy"; path.touch()
    history = np.zeros((3, 4))
    monkeypatch.setattr(one.previous, "_working_history", lambda *args: (path, history))
    def save(path, ref, data, times, **kwargs):
        np.testing.assert_array_equal(data[0, :2], initial)
        np.testing.assert_array_equal(data[0, 2:], np.zeros(2))
        assert kwargs["linear"] is True
        assert kwargs["execution"]["initial_kind"] == "nonlinear"
        assert kwargs["execution"]["ODE_integrations"] == 0
        assert kwargs["execution"]["no_modal_truncation"] is True
    monkeypatch.setattr(one.previous, "_save_trajectory", save)
    one._linear(tmp_path, "same_initial_state", reference, np.array([0., .5, 1.]),
                "nonlinear", "explicit_user_FEM3C_2026_10_09")
    assert len(calls) == 1 and not path.exists()


@pytest.mark.parametrize("flag", ["one_d_completed", "hard_stop"])
def test_full_period_cached_or_stopped_does_not_construct_reference(flag, tmp_path, monkeypatch):
    one = importlib.import_module("scripts.lib.nlsp_fem3c_1d")
    monkeypatch.setattr(one, "_reference", lambda *args: pytest.fail("Cached/stopped stage constructs reference"))
    monkeypatch.setattr(one.previous.base.one.runner, "integrate_case",
                        lambda *args, **kwargs: pytest.fail("Cached/stopped stage starts ODE"))
    summary = {flag: True}
    assert one.run_full_period(tmp_path, {}, summary) is summary


def test_separate_authorization_and_frozen_continuation_scope(continuation, config):
    assert continuation.validate_config(config) is config
    assert config["authorization"]["id"] == "explicit_user_FEM3C_2026_10_09"
    assert config["authorization"]["maximum_production_CCX_jobs"] == 6
    assert config["authorization"]["maximum_nonlinear_1D_integrations"] == 2
    assert config["authorization"]["automatic_retry"] is False
    assert config["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED"
    assert config["admitted"] is False
    assert config["parent_completed"]["manifest_sha256"] == PARENT_SHA
    assert continuation.CASE_ORDER == (
        ("medium_refined_time", "linear"), ("medium_refined_time", "nonlinear"),
        ("fine_refined_time", "linear"), ("fine_refined_time", "nonlinear"),
        ("full_period_medium", "linear"), ("full_period_medium", "nonlinear"))
    assert config["threads"] == 1 and config["job_memory_limit_bytes"] == 4*1024**3
    assert config["job_timeout_seconds"] <= 5400
    assert config["numerical_budget_seconds"] <= 24000


@pytest.mark.parametrize("path,value", [
    (("authorization", "id"), "explicit_user_FEM3B_2026_10_09"),
    (("authorization", "automatic_retry"), True),
    (("authorization", "maximum_production_CCX_jobs"), 7),
    (("authorization", "maximum_nonlinear_1D_integrations"), 3),
    (("stages", "medium_refined_time", "initial_T1_fraction"), 1/16000),
    (("stages", "fine_refined_time", "maximum_T1_fraction"), 1/2000),
    (("stages", "fine_refined_time", "mesh_level"), "refined"),
    (("stages", "full_period_medium", "maximum_increments"), 2000),
    (("stages", "full_period_medium", "conditional_on_robustness"), False),
    (("dynamic", "alpha"), -.05),
    (("comparison", "points"), 200),
    (("comparison", "baseline_model_discrepancy"), 3e-7),
    (("comparison", "temporal_ratio_limit"), .3),
    (("comparison", "spatial_ratio_limit"), .3),
    (("comparison", "phase_amplitude_fitting"), True),
    (("one_d", "time_level"), "loose"),
    (("one_d", "main_degree"), 96),
    (("one_d", "target_T1_fraction"), 5.),
    (("execution_mode",), "STRICT_ADMITTED"),
    (("admitted",), True),
    (("threads",), 2),
    (("new_meshes",), True),
    (("new_modal_jobs",), True),
    (("new_static_only_jobs",), True),
    (("seven_field_observations", "three_d_c"), "identical_MH_coordinate"),
])
def test_frozen_settings_cannot_be_changed_for_agreement(continuation, config, path, value):
    current = config
    for key in path[:-1]:
        current = current[key]
    current[path[-1]] = value
    with pytest.raises(ValueError):
        continuation.validate_config(config)


@pytest.mark.parametrize("stage", ["medium_refined_time", "fine_refined_time", "full_period_medium"])
def test_derived_case_preserves_physics_and_refines_only_declared_policy(continuation, config, science, stage):
    original = copy.deepcopy(science)
    derived = continuation.case_science(science, config, stage)
    assert science == original
    for key in ("geometry", "material", "g", "q", "omega1", "preload_reproduction",
                "source_resume", "source_static", "source_fem1", "source_action"):
        assert derived[key] == original[key]
    policy = config["stages"][stage]
    assert derived["mesh_level"] == policy["mesh_level"]
    assert derived["horizon_T1"] == policy["horizon_T1"]
    settings = continuation.base.dynamic_settings(derived)
    assert settings["minimum_increment"] == settings["initial_increment"]*1e-4
    assert settings["alpha"] == 0
    assert settings["maximum_increments"] == policy["maximum_increments"]
    if stage in continuation.CONTROL_STAGES:
        assert settings["initial_increment"] == settings["T1"]/8000
        assert settings["maximum_increment"] == settings["T1"]/4000
        assert settings["duration"] == .25*settings["T1"]
    else:
        assert settings["duration"] == settings["T1"]
        assert settings["maximum_increments"] > 2000


def test_corrupt_parent_stops_before_any_scientific_call(continuation, config, monkeypatch):
    _forbid_science(monkeypatch, continuation)
    config["parent_completed"]["manifest_sha256"] = "0"*64
    monkeypatch.setattr(continuation.previous, "validate_cache",
                        lambda *args: pytest.fail("Corrupt source accepted for replay"))
    with pytest.raises(ValueError, match="manifest changed"):
        continuation.load_parent(config)


@pytest.mark.parametrize("level,expected", [("medium", (5649, 3120)), ("fine", (11553, 6670))])
def test_each_source_mesh_and_preload_keeps_its_own_level(continuation, science, monkeypatch, level, expected):
    _forbid_science(monkeypatch, continuation)
    medium = SimpleNamespace(nodes=range(5649), solid_elements=range(3120))
    fine = SimpleNamespace(nodes=range(11553), solid_elements=range(6670))
    old = {"science_config": {"frozen": True}}
    calls = []
    monkeypatch.setattr(continuation.base, "verify_sources", lambda cfg: (old, ROOT/"medium", medium, {"level": "medium"}))
    monkeypatch.setattr(continuation.base.base, "load_fem2_sources", lambda cfg: {"saved": True})
    def saved(config, requested, sources):
        assert config == old["science_config"] and requested == "fine" and sources == {"saved": True}
        calls.append(requested)
        return ROOT/"fine", fine, {"level": "fine"}
    monkeypatch.setattr(continuation.base.base, "fem2_source_mesh", saved)
    _, source, mesh, audit = continuation.source_for_level(science, level)
    assert (len(mesh.nodes), len(mesh.solid_elements)) == expected
    assert source.name == level and audit["level"] == level
    assert calls == (["fine"] if level == "fine" else [])


@pytest.mark.parametrize("level", ["coarse", "refined", "new", "fine2"])
def test_no_extra_mesh_level_can_be_routed(continuation, science, monkeypatch, level):
    _forbid_science(monkeypatch, continuation)
    monkeypatch.setattr(continuation.base, "verify_sources", lambda *args: pytest.fail("Unauthorized mesh source read"))
    with pytest.raises(ValueError, match="medium/fine"):
        continuation.source_for_level(science, level)


@pytest.mark.parametrize("stage", ["medium_refined_time", "full_period_medium"])
@pytest.mark.parametrize("kind", ["linear", "nonlinear"])
def test_output_adapter_changes_only_dynamic_cadence_and_keeps_safe_static(continuation, config, science, monkeypatch, tmp_path, stage, kind):
    _forbid_science(monkeypatch, continuation)
    derived = continuation.case_science(science, config, stage)
    original = (PARENT / "cases" / kind / "motion.inp").read_text(encoding="utf8")
    lines = original.splitlines()
    d = continuation.base.dynamic_settings(derived)
    for i, line in enumerate(lines):
        if line.startswith("*DYNAMIC"):
            lines[i+1] = ",".join(format(d[k], ".15g") for k in
                ("initial_increment", "duration", "minimum_increment", "maximum_increment"))
    prepared = "\n".join(lines)+"\n"
    def generator(path, *args):
        Path(path).write_text(prepared, encoding="utf8")
        return {"status": "PASS"}
    monkeypatch.setattr(continuation.base, "write_input", generator)
    path = tmp_path / "motion.inp"
    gate = continuation.write_case_input(path, derived, {}, None, None, None, kind == "nonlinear")
    text = path.read_text(encoding="utf8")
    static_before, dynamic_before = prepared.split("*END STEP\n", 1)
    static_after, dynamic_after = text.split("*END STEP\n", 1)
    assert static_after == static_before
    frequency = config["stages"][stage]["output_frequency"]
    expected_dynamic = dynamic_before.replace("FREQUENCY=1", f"FREQUENCY={frequency}")
    marker = f"*NODE FILE, GLOBAL=YES, FREQUENCY={frequency}"
    expected_dynamic = expected_dynamic.replace(marker,
        f"*EL FILE, GLOBAL=YES, FREQUENCY={frequency}\nS,E,ENER\n"+marker, 1)
    assert dynamic_after == expected_dynamic
    assert f"*EL FILE, GLOBAL=YES, FREQUENCY={frequency}\nS,E,ENER\n" in dynamic_after
    assert all(f"FREQUENCY={frequency}" in line for line in dynamic_after.splitlines()
               if line.startswith(("*NODE FILE", "*EL FILE", "*NODE PRINT", "*EL PRINT")))
    assert gate["static_output_frequency"] == 1
    assert gate["dynamic_output_frequency"] == frequency
    assert gate["final_frame_required"] is True
    assert gate["physical_generation_unchanged"] is True
    if kind == "linear":
        assert "ELKE" not in static_after
    assert "*DLOAD, OP=NEW\nSOLID,GRAV,0.,0.,-1.,0." in dynamic_after
    assert "*DYNAMIC, ALPHA=0" in dynamic_after


@pytest.mark.parametrize("kind", ["linear", "nonlinear"])
def test_support_recovery_subtracts_independent_bodyload_for_own_mesh_level(continuation, monkeypatch, tmp_path, kind):
    _forbid_science(monkeypatch, continuation)
    (tmp_path / "frames").mkdir(); (tmp_path / "dat_fields").mkdir()
    (tmp_path / "motion.stdout.txt").write_text("Job finished\n", encoding="utf8")
    continuation.write_json(tmp_path / "increments.json", {
        "accepted_increments": [{"step": 1, "increment": 10, "total_time": 1.}]})
    ids = np.array([10, 20]); xyz = np.array([[0., 0., 0.], [1., 0., 0.]])
    gravity = np.array([[0., -.5, 0.], [0., -.5, 0.]])
    np.savez(tmp_path / "frames/step1_inc00010.npz", U=np.zeros((2, 3)))
    # Raw RF includes consistent applied bodyload at the restrained nodes.
    # Subtraction, rather than blind RF summation, recovers the support forces.
    for name in ("LEFT_FIXED", "RIGHT_FIXED"):
        np.savez(tmp_path / f"dat_fields/1_10_{name}_FORC.npz", values=np.zeros((1, 3)))
    calls = []
    def independent(mesh, rho, g):
        calls.append((mesh, rho, g))
        return ids, xyz, gravity, 1.
    monkeypatch.setattr(continuation.base.base, "consistent_gravity_loads", independent)
    science = {"material": {"rho": 1.}, "g": .0014224751066856333, "mesh_level": "fine"}
    old = {"science_config": {"gates": {"equilibrium_relative": 1e-5}}}
    record = {"maximum_native_external_work_after_release": 0.,
              "maximum_native_damping_work_after_release": 0.}
    mesh = object()
    audit = {"fixed_left_ids": [10], "fixed_right_ids": [20]}
    result = continuation.audit_saved_case(tmp_path, science, old, mesh, audit, kind, record)
    assert calls == [(mesh, 1., .0014224751066856333)]
    assert result["source_mesh_level"] == "fine" and result["status"] == "PASS"
    assert result["preload_force_imbalance_relative"] == 0
    assert result["preload_moment_imbalance_relative"] == 0
    for support in result["support_resultants"].values():
        np.testing.assert_array_equal(support["force"], [0., .5, 0.])
    assert result["equilibrium_gate_unchanged"] == 1e-5
    assert result["bodyload_recovered_independently"] is True
    assert record["continuation_audit"] is result


def _state():
    return {"cases": {}, "attempts": [], "statuses": {}, "numerical_seconds": 0.,
            "job_calls": {"CCX_production": 0}, "overall": "NOT_RUN"}


def _decision(continuation, bundle, state):
    continuation.write_json(bundle / "pre_fem_decision.json", {"before_results": True})
    state["pre_fem_decision_sha256"] = continuation.sha(bundle / "pre_fem_decision.json")


@pytest.mark.parametrize("stage,kind", [
    ("medium_refined_time", "nonlinear"), ("fine_refined_time", "linear"),
    ("fine_refined_time", "nonlinear"), ("full_period_medium", "linear"),
])
def test_case_order_requires_actual_preceding_success_and_audit(continuation, config, science, monkeypatch, tmp_path, stage, kind):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    monkeypatch.setattr(continuation, "source_for_level", lambda *args: pytest.fail("Missing prior gate loads next case"))
    item = {"config": science, "validation_config": config}
    with pytest.raises(ValueError, match="Sequential production gate"):
        continuation.run_case(tmp_path, item, state, stage, kind)
    assert state["attempts"] == [] and state["job_calls"]["CCX_production"] == 0


@pytest.mark.parametrize("prior_audit", [None, "FAIL", "PARTIAL"])
def test_linear_exit_success_without_recovery_audit_does_not_admit_nonlinear(continuation, config, science, monkeypatch, tmp_path, prior_audit):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    state["cases"]["medium_refined_time/linear"] = {"status": "PASS"}
    if prior_audit:
        state["cases"]["medium_refined_time/linear"]["continuation_audit"] = {"status": prior_audit}
    item = {"config": science, "validation_config": config}
    with pytest.raises(ValueError, match="Sequential production gate"):
        continuation.run_case(tmp_path, item, state, "medium_refined_time", "nonlinear")


def test_full_period_3d_requires_frozen_positive_robustness_decision(continuation, config, science, monkeypatch, tmp_path):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    for stage in continuation.CONTROL_STAGES:
        for kind in ("linear", "nonlinear"):
            state["cases"][stage+"/"+kind] = {"status": "PASS", "continuation_audit": {"status": "PASS"}}
    monkeypatch.setattr(continuation, "source_for_level", lambda *args: pytest.fail("Full-period gate not checked before source"))
    item = {"config": science, "validation_config": config}
    with pytest.raises(ValueError, match="frozen successful numerical-robustness"):
        continuation.run_case(tmp_path, item, state, continuation.FULL_STAGE, "linear")
    continuation.write_json(tmp_path / "full_period_decision.json", {"full_period_3D_authorized_by_actual_gates": False})
    state["full_period_decision_sha256"] = continuation.sha(tmp_path / "full_period_decision.json")
    with pytest.raises(ValueError, match="frozen successful numerical-robustness"):
        continuation.run_case(tmp_path, item, state, continuation.FULL_STAGE, "linear")
    assert state["job_calls"]["CCX_production"] == 0


@pytest.mark.parametrize("status", ["STARTED", "FAIL", "PARTIAL"])
def test_interrupted_or_failed_attempt_never_auto_retries(continuation, config, science, monkeypatch, tmp_path, status):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    state["cases"]["medium_refined_time/linear"] = {"status": status}
    state["job_calls"]["CCX_production"] = 1
    monkeypatch.setattr(continuation, "source_for_level", lambda *args: pytest.fail("Interrupted case loads new work"))
    assert continuation.run_case(tmp_path, {"config": science, "validation_config": config}, state,
                                 "medium_refined_time", "linear") is False
    assert state["job_calls"]["CCX_production"] == 1


@pytest.mark.parametrize("returncode,log", [(1, "Job finished\n"), (0, "*ERROR Invalid input\n"), (0, "Interrupted\n")])
def test_real_failure_and_incomplete_native_log_stop_series_without_retry(continuation, config, science, monkeypatch, tmp_path, returncode, log):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    source = tmp_path / "source"; source.mkdir()
    (source / "solid_mesh.inp").write_text("immutable existing mesh", encoding="utf8")
    old = {"science_config": {"ccx_exe": "saved_ccx.exe"}}
    monkeypatch.setattr(continuation, "source_for_level", lambda *args: (old, source, object(), {}))
    def input_fixture(path, *args):
        Path(path).write_text("fixture input, never executed", encoding="utf8")
        return {"status": "PASS"}
    monkeypatch.setattr(continuation, "write_case_input", input_fixture)
    snapshots = []
    monkeypatch.setattr(continuation, "save", lambda *args: snapshots.append(copy.deepcopy(state)))
    calls = []
    def fixture_job(cmd, cwd, timeout, memory, prefix, env):
        calls.append(cmd)
        assert snapshots[-1]["attempts"][0]["status"] == "STARTED"
        assert snapshots[-1]["attempts"][0]["case"] == "medium_refined_time/linear"
        assert snapshots[-1]["cases"]["medium_refined_time/linear"]["status"] == "STARTED"
        assert snapshots[-1]["job_calls"]["CCX_production"] == 1
        assert timeout <= 5400 and memory == 4*1024**3
        assert all(env[key] == "1" for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS"))
        (Path(cwd) / "motion.stdout.txt").write_text(log, encoding="utf8")
        return SimpleNamespace(returncode=returncode), {"failure": None, "seconds": .001}
    monkeypatch.setattr(continuation.base.base.fem1, "run_job", fixture_job)
    monkeypatch.setattr(continuation, "recover_saved_case", lambda *args: pytest.fail("Invalid native result enters recovery"))
    frozen = (Path(continuation.base.__file__), Path(continuation.base.base.__file__),
              Path(continuation.base.base.fem1.__file__),
              Path(continuation.base.one.dynamics.__file__), Path(continuation.base.one.rod.__file__))
    # Different helper names exercise the snapshot loops while the original
    # scientific case identity must remain unchanged in the actual ledger.
    item = {"config": science, "validation_config": config,
            "helper_sha256": {p.relative_to(ROOT).as_posix(): _sha(p) for p in frozen}}
    (tmp_path / "execution_code").mkdir()
    assert continuation.run_case(tmp_path, item, state, "medium_refined_time", "linear") is False
    assert state["hard_stop"] is True and state["overall"] == "BLOCKED_BY_SOLVER"
    assert state["attempts"][0]["status"] == "FAIL"
    assert state["attempts"][0]["automatic_retry"] is False
    assert continuation.run_case(tmp_path, item, state, "medium_refined_time", "linear") is False
    assert len(calls) == 1 and state["job_calls"]["CCX_production"] == 1


def test_finished_native_case_can_be_reparsed_without_another_solver_call(continuation, config, science, monkeypatch, tmp_path):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    name = "medium_refined_time/linear"
    state["cases"][name] = {"status": "OUTPUT_RECOVERY_PENDING", "job": {"returncode": 0}}
    state["attempts"] = [{"case": name, "status": "SOLVER_FINISHED"}]
    state["job_calls"]["CCX_production"] = 1
    monkeypatch.setattr(continuation, "save", lambda *args: None)
    monkeypatch.setattr(continuation, "source_for_level", lambda *args: ({}, None, None, None))
    calls = []
    def reparse(*args):
        calls.append("reparse")
        return {"status": "PASS", "native_fields_preserved": True}
    def audit(*args):
        args[-1]["continuation_audit"] = {"status": "PASS"}
    monkeypatch.setattr(continuation, "recover_saved_case", reparse)
    monkeypatch.setattr(continuation, "audit_saved_case", audit)
    item = {"config": science, "validation_config": config}
    assert continuation.run_case(tmp_path, item, state, "medium_refined_time", "linear") is True
    assert calls == ["reparse"] and state["job_calls"]["CCX_production"] == 1
    assert state["attempts"][0]["status"] == "PASS"
    assert state["cases"][name]["continuation_audit"]["status"] == "PASS"


@pytest.mark.parametrize("condition,allowed", [("all", True), ("robustness", False), ("audit", False), ("resource", False)])
def test_full_period_decision_requires_all_actual_controls_and_resources(continuation, config, monkeypatch, tmp_path, condition, allowed):
    _forbid_science(monkeypatch, continuation)
    state = _state()
    for stage in continuation.CONTROL_STAGES:
        for kind in ("linear", "nonlinear"):
            state["cases"][stage+"/"+kind] = {"status": "PASS",
                "continuation_audit": {"status": "PASS"}, "job": {"seconds": 10.}}
    if condition == "audit":
        state["cases"]["fine_refined_time/nonlinear"]["continuation_audit"]["status"] = "PARTIAL"
    if condition == "resource":
        state["numerical_seconds"] = config["numerical_budget_seconds"]-5.
    continuation.write_json(tmp_path / "robustness_comparison.json", {
        "full_period_numerical_robustness_gate": condition != "robustness"})
    continuation.write_json(tmp_path / "pre_fem_decision.json", {
        "planning_CCX_estimates_seconds": {continuation.FULL_STAGE: 40.}})
    monkeypatch.setattr(continuation, "save", lambda *args: None)
    item = {"validation_config": config}
    assert continuation.freeze_full_period(tmp_path, item, state) is allowed
    decision = continuation.read_json(tmp_path / "full_period_decision.json")
    assert decision["full_period_3D_authorized_by_actual_gates"] is allowed
    assert decision["before_full_period_results"] is True
    assert decision["illustration_policy"] == config["stages"][continuation.FULL_STAGE]
    assert decision["no_choice_by_future_1D_3D_agreement"] is True
    assert decision["full_period_convergence_not_certified_by_quarter_period_control"] is True
    assert continuation.freeze_full_period(tmp_path, item, state) is allowed
    if not allowed:
        assert "no full-period 3D jobs" in state["full_period_stop_reason"]
    # A result must not be overwritten by a later reclassification.
    decision["full_period_3D_authorized_by_actual_gates"] = not allowed
    continuation.write_json(tmp_path / "full_period_decision.json", decision)
    with pytest.raises(ValueError, match="decision changed"):
        continuation.freeze_full_period(tmp_path, item, state)


@pytest.mark.parametrize("mode", ["--report-only", "--plot-only"])
def test_report_plot_are_saved_data_routes_with_zero_scientific_calls(continuation, monkeypatch, tmp_path, mode, capsys):
    _forbid_science(monkeypatch, continuation)
    state = _state(); state["overall"] = "PARTIAL"
    calls = []
    monkeypatch.setattr(continuation, "validate_cache", lambda path: calls.append("validate") or state)
    monkeypatch.setattr(continuation, "prepare_stage", lambda *args: pytest.fail("Report/plot prepare a new scientific stage"))
    monkeypatch.setattr(continuation, "save", lambda *args: calls.append("save"))
    continuation.write_json(tmp_path / "provenance.json", {"fixture": True})
    diagnostics = importlib.import_module("scripts.lib.nlsp_fem3c_diagnostics")
    monkeypatch.setattr(continuation, "plot_saved", lambda *args: calls.append("plot") or [])
    assert continuation.main([mode, str(tmp_path)]) is state
    assert calls == (["validate", "plot", "save"] if mode == "--plot-only" else ["validate"])
    assert json.loads(capsys.readouterr().out)["new_scientific_calls"] == 0


def test_completed_compute_never_routes_native_or_one_d_stages(continuation, monkeypatch, tmp_path, capsys):
    _forbid_science(monkeypatch, continuation)
    state = _state(); state["completed"] = True
    monkeypatch.setattr(continuation, "prepare_stage", lambda *args: (tmp_path, {}, state))
    monkeypatch.setattr(continuation, "run_controls", lambda *args: pytest.fail("Completed cache repeats controls"))
    monkeypatch.setattr(continuation, "run_case", lambda *args: pytest.fail("Completed cache repeats native job"))
    one = importlib.import_module("scripts.lib.nlsp_fem3c_1d")
    monkeypatch.setattr(one, "run_full_period", lambda *args: pytest.fail("Completed cache repeats nonlinear 1D solve"))
    assert continuation.main(["--compute"]) is state
    assert json.loads(capsys.readouterr().out)["new_scientific_calls"] == 0


def test_unmanifested_attempt_cannot_be_replaced_after_code_hash_change(continuation, config, monkeypatch, tmp_path):
    _forbid_science(monkeypatch, continuation)
    output = tmp_path / "namespace"; bundle = output / "old_fingerprint"
    bundle.mkdir(parents=True)
    continuation.write_json(bundle / "provenance.json", {
        "authorization": config["authorization"], "validation_config": config,
        "helper_sha256": {"different_historical_code": "0"*64}})
    monkeypatch.setattr(continuation, "OUTPUT", output)
    with pytest.raises(RuntimeError, match="no automatic retry"):
        continuation.existing_attempt(config)


def test_same_authorization_rejects_a_changed_cache_configuration(continuation, config, monkeypatch, tmp_path):
    _forbid_science(monkeypatch, continuation)
    output = tmp_path / "namespace"; bundle = output / "old_fingerprint"
    bundle.mkdir(parents=True)
    old = copy.deepcopy(config); old["comparison"]["points"] = 203
    continuation.write_json(bundle / "provenance.json", {
        "authorization": config["authorization"], "validation_config": old})
    monkeypatch.setattr(continuation, "OUTPUT", output)
    with pytest.raises(ValueError, match="already used with another config"):
        continuation.existing_attempt(config)


def test_validation_dispatch_preserves_the_existing_entrypoint(continuation, monkeypatch):
    calls = []
    monkeypatch.setattr(continuation, "main", lambda args: calls.append(args) or {"cached": True})
    result = continuation.resume.main(["--validation", "--report-only", "saved_bundle"])
    assert result == {"cached": True}
    assert calls == [["--report-only", "saved_bundle"]]


@pytest.mark.parametrize("same_mesh", [True, False])
def test_preflight_and_production_depths_resolve_to_the_same_frozen_mesh(continuation, config, science, monkeypatch, tmp_path, same_mesh):
    _forbid_science(monkeypatch, continuation)
    state = _state(); _decision(continuation, tmp_path, state)
    source = tmp_path / "source"; source.mkdir()
    mesh_path = source / "solid_mesh.inp"
    mesh_path.write_text("immutable source", encoding="utf8")
    other = source / "other_mesh.inp"; other.write_text("different mesh", encoding="utf8")
    old = {"science_config": {"ccx_exe": "fixture.exe"}}
    monkeypatch.setattr(continuation, "source_for_level", lambda *args: (old, source, None, None))
    def text(path, actual_mesh):
        include = os.path.relpath(actual_mesh, Path(path).parent)
        return (f"*INCLUDE, INPUT={include}\n*STEP\n*STATIC\n*END STEP\n"
                "*STEP\n*DYNAMIC, ALPHA=0\n*NODE FILE, GLOBAL=YES, FREQUENCY=2\nU,V,RF\n*END STEP\n")
    preview = tmp_path / "input_gate/medium_refined_time/linear.inp"
    preview.parent.mkdir(parents=True)
    preview.write_text(text(preview, mesh_path), encoding="utf8")
    def generated(path, *args):
        Path(path).write_text(text(path, mesh_path if same_mesh else other), encoding="utf8")
        return {"status": "PASS"}
    monkeypatch.setattr(continuation, "write_case_input", generated)
    monkeypatch.setattr(continuation, "save", lambda *args: None)
    frozen = (Path(continuation.base.__file__), Path(continuation.base.base.__file__),
              Path(continuation.base.base.fem1.__file__), Path(continuation.base.one.dynamics.__file__),
              Path(continuation.base.one.rod.__file__))
    item = {"config": science, "validation_config": config,
            "helper_sha256": {p.relative_to(ROOT).as_posix(): _sha(p) for p in frozen}}
    (tmp_path / "execution_code").mkdir()
    calls = []
    def native_fixture(cmd, cwd, *args):
        calls.append(cmd)
        (Path(cwd) / "motion.stdout.txt").write_text("fixture scientific failure", encoding="utf8")
        return SimpleNamespace(returncode=1), {"failure": None, "seconds": .001}
    monkeypatch.setattr(continuation.base.base.fem1, "run_job", native_fixture)
    if same_mesh:
        assert continuation.run_case(tmp_path, item, state, "medium_refined_time", "linear") is False
        revision = continuation.read_json(tmp_path / "cases/medium_refined_time/linear/input_adapter_revision.json")
        assert revision["physical_preflight_equivalence"] is True
        proof = revision["resolved_include_evidence"]
        assert proof["preflight"]["written"] != proof["production"]["written"]
        assert proof["preflight"]["resolved"] == proof["production"]["resolved"] == str(mesh_path.resolve())
        assert proof["preflight"]["sha256"] == proof["production"]["sha256"] == _sha(mesh_path)
        assert len(calls) == state["job_calls"]["CCX_production"] == 1
    else:
        with pytest.raises(ValueError, match="frozen source mesh"):
            continuation.run_case(tmp_path, item, state, "medium_refined_time", "linear")
        assert calls == [] and state["job_calls"]["CCX_production"] == 0
        assert state["attempts"] == []


@pytest.mark.parametrize("problem", [None, "missing", "audit", "started", "hard_stop"])
def test_full_period_one_d_waits_for_the_bounded_c1_outcome(continuation, problem):
    state = _state()
    for stage in continuation.CONTROL_STAGES:
        for kind in ("linear", "nonlinear"):
            state["cases"][stage+"/"+kind] = {"status": "PASS", "continuation_audit": {"status": "PASS"}}
    key = "fine_refined_time/nonlinear"
    if problem == "missing":
        state["cases"].pop(key)
    elif problem == "audit":
        state["cases"][key]["continuation_audit"]["status"] = "PARTIAL"
    elif problem == "started":
        state["cases"][key]["status"] = "STARTED"
    elif problem == "hard_stop":
        state["hard_stop"] = True
    assert continuation.one_d_stage_allowed(state) is (problem is None)
    if problem in ("missing", "audit"):
        # Only an explicit, non-solver C1 stop qualification permits the
        # optional saved-IC 1D illustration after an incomplete C1 outcome.
        state["c1_terminated_without_solver_failure"] = True
        assert continuation.one_d_stage_allowed(state) is True


def test_explicit_one_d_cli_before_c1_is_rejected_without_ode(continuation, monkeypatch, tmp_path):
    _forbid_science(monkeypatch, continuation)
    state = _state()
    monkeypatch.setattr(continuation, "prepare_stage", lambda *args: (tmp_path, {}, state))
    monkeypatch.setattr(continuation, "record_postprocessing_phase", lambda *args: pytest.fail("Premature 1D phase begins"))
    with pytest.raises(ValueError, match="bounded C1 outcome"):
        continuation.main(["--run-1d"])
    assert state["job_calls"]["CCX_production"] == 0


@pytest.mark.parametrize("robust", [True, False])
def test_fresh_compute_records_actual_call_deltas_then_closes_qualified_cache(continuation, config, science, monkeypatch, tmp_path, robust, capsys):
    """Every numeric stage is a synthetic fixture; no real solver runs."""
    _forbid_science(monkeypatch, continuation)
    state = _state(); state["job_calls"]["1D_nonlinear_ODE"] = 0
    item = {"config": science, "validation_config": config}
    calls = []
    monkeypatch.setattr(continuation, "prepare_stage", lambda *args: (tmp_path, item, state))
    def controls(*args):
        calls.append("controls")
        for stage in continuation.CONTROL_STAGES:
            for kind in ("linear", "nonlinear"):
                state["cases"][stage+"/"+kind] = {"status": "PASS", "continuation_audit": {"status": "PASS"}}
        state["job_calls"]["CCX_production"] = 4
    def robustness(*args):
        assert continuation.controls_ready(state)
        calls.append("robustness")
        state["robustness_gate"] = robust
        continuation.write_json(tmp_path / "robustness_comparison.json", {"full_period_numerical_robustness_gate": robust})
    def one_d(*args):
        assert continuation.controls_ready(state)
        calls.append("one_d")
        state["one_d_completed"] = True
        state["job_calls"]["1D_nonlinear_ODE"] = 2
    def freeze(*args):
        assert state["one_d_completed"] is True
        calls.append("full_period_decision")
        state["full_period_3D_allowed"] = robust
        if not robust:
            state["full_period_stop_reason"] = "Numerical robustness gate fails; no full-period 3D jobs"
        return robust
    def native(bundle, identity, summary, stage, kind):
        assert robust and stage == continuation.FULL_STAGE
        if kind == "nonlinear":
            assert summary["cases"][stage+"/linear"]["continuation_audit"]["status"] == "PASS"
        calls.append("full_"+kind)
        state["cases"][stage+"/"+kind] = {"status": "PASS", "continuation_audit": {"status": "PASS"}}
        state["job_calls"]["CCX_production"] += 1
        return True
    monkeypatch.setattr(continuation, "run_controls", controls)
    monkeypatch.setattr(continuation, "complete_robustness", robustness)
    one = importlib.import_module("scripts.lib.nlsp_fem3c_1d")
    monkeypatch.setattr(one, "run_full_period", one_d)
    monkeypatch.setattr(continuation, "freeze_full_period", freeze)
    monkeypatch.setattr(continuation, "run_case", native)
    monkeypatch.setattr(continuation, "save", lambda *args: None)
    monkeypatch.setattr(continuation, "record_postprocessing_phase", lambda *args: calls.append("phase:"+args[-1]))
    diagnostics = importlib.import_module("scripts.lib.nlsp_fem3c_diagnostics")
    def compared(*args):
        assert robust and state["job_calls"]["CCX_production"] == 6
        calls.append("full_compare")
        return {"one_d_spatial_status": "PARTIAL"}
    monkeypatch.setattr(diagnostics, "full_period", compared)
    monkeypatch.setattr(continuation, "plot_saved", lambda *args: calls.append("plot_saved") or [])
    assert continuation.main(["--compute"]) is state
    printed = json.loads(capsys.readouterr().out)
    assert printed["new_scientific_calls"] == {"CCX_production": 6 if robust else 4, "1D_nonlinear_ODE": 2}
    assert state["completed"] is True
    assert state["universal_nonlinear_validation"] is False
    assert state["experimental_validation"] is False
    assert calls.index("controls") < calls.index("robustness") < calls.index("one_d")
    if robust:
        assert calls.index("one_d") < calls.index("full_period_decision") < calls.index("full_linear") < calls.index("full_nonlinear") < calls.index("full_compare")
        assert state["overall"] == "STRAIGHT_ROD_NONLINEAR_3D_FEM_VERIFICATION_COMPLETE_WITH_QUALIFICATIONS"
        assert state["statuses"]["NLSP_FEM3C_ENERGY_DIAGNOSTICS"] == "PARTIAL"
        assert state["full_period_1d_spatial_qualification"] == "PARTIAL"
        assert state["full_period_illustrative_only"] is True
    else:
        assert "full_linear" not in calls and "full_nonlinear" not in calls and "full_compare" not in calls
        assert state["overall"] == "PARTIAL"
        assert state["statuses"]["NLSP_FEM3C_VERIFICATION_SUMMARY"] == "PARTIAL"
    calls.clear()
    assert continuation.main(["--compute"]) is state
    assert calls == []
    assert json.loads(capsys.readouterr().out)["new_scientific_calls"] == 0


def test_actual_completed_c1_cases_keep_own_preloads_and_real_sparse_output():
    if not (NEW_BUNDLE / "robustness_comparison.json").exists():
        pytest.skip("Completed local FEM-3C1 evidence unavailable")
    for level, stage in (("medium", "medium_refined_time"), ("fine", "fine_refined_time")):
        for kind in ("linear", "nonlinear"):
            case = NEW_BUNDLE / "cases" / stage / kind
            record = _read(case / "recovery.json")
            audit = _read(case / "continuation_audit.json")
            assert record["status"] == audit["status"] == "PASS"
            assert record["dynamic_increments"] == 1002
            assert record["dynamic_output_frames"] == 501
            assert record["dynamic_output_frames"] != record["dynamic_increments"]
            assert record["output_frequency"] == 2
            assert record["cutbacks"] == 0 and audit["warning_lines"] == []
            assert audit["source_mesh_level"] == level
            assert audit["equilibrium_gate_unchanged"] == 1e-5
            assert audit["preload_force_imbalance_relative"] <= 1e-5
            assert audit["preload_moment_imbalance_relative"] <= 1e-5
            assert record["final_native_frame_reached"] is True
            assert record["dynamic_time_end"] == .25*10.37828159055014
            assert record["maximum_native_external_work_after_release"] == 0
            assert record["maximum_native_damping_work_after_release"] == 0
            assert record["energy_status"] == "PARTIAL"
            raw_jump = (record["native_dynamic_bookkeeping_initial_energy"]
                        / record["native_initial_internal_energy"] - 1.)
            np.testing.assert_allclose(record["native_energy_static_to_dynamic_reference_jump_relative"],
                                       raw_jump, rtol=2*np.finfo(float).eps, atol=0)
            assert raw_jump != 0. and record["native_energy_history_not_offset_corrected"] is True
            assert record["all_saved_frames_strain_diagnostics"]["minimum_det_deformation_gradient"] > 0
            preload = record["preload_transfer"]
            assert f"/cases/{level}/{kind}" in preload["source_static_case"]
            for metric in ("node_displacement_max_difference", "support_RF_max_difference",
                           "stress_difference", "strain_difference", "section_profile_max_difference"):
                assert preload[metric] == 0
            with np.load(case / "section_history.npz") as history, np.load(case / "initial_sections.npz") as initial:
                assert history["time"][0] > 0
                assert history["time"][-1] == record["dynamic_time_end"]
                assert np.all(np.diff(history["time"]) > 0)
                assert history["fields"].shape == (501, 41, 7)
                assert bool(initial["not_a_native_dynamic_zero_frame"]) is True


def test_actual_c1_effects_use_fixed_baseline_and_explicit_model_denominators():
    if not (NEW_BUNDLE / "robustness_comparison.json").exists():
        pytest.skip("Completed local FEM-3C1 evidence unavailable")
    result = _read(NEW_BUNDLE / "robustness_comparison.json")
    policy = result["comparison_policy"]
    assert policy["baseline_model_discrepancy"] == 2.824717e-7
    assert policy["temporal_ratio_limit"] == policy["spatial_ratio_limit"] == .25
    assert policy["phase_amplitude_fitting"] is False
    assert result["comparison_points"] == 201 and result["section_points"] == 41
    assert result["scientific_calls"] == {"CCX": 0, "Gmsh": 0, "nonlinear_ODE": 0,
                                          "eigenanalysis": 0, "static_equilibrium": 0}
    with np.load(NEW_BUNDLE / "robustness_comparison.npz") as data:
        for control, first, second in (("temporal", "medium_refined_time", "existing_medium"),
                                        ("spatial", "fine_refined_time", "medium_refined_time")):
            delta = data[first+"_evolution_fields"][..., 1]-data[second+"_evolution_fields"][..., 1]
            observed = float(np.max(abs(delta)))
            row = result["refinement"][control]
            assert row["evolution_w_max"] == observed
            assert row["ratio_to_fixed_baseline_model_discrepancy"] == observed/2.824717e-7
            assert row["ratio_limit"] == .25 and row["status"] == "PASS"
            assert row["error_upper_bound_claimed"] is False
    for record in result["resolutions"].values():
        metric = record["model_comparisons"]["evolution"]["w"]
        assert metric["relative_max"] == metric["absolute_max"]/metric["characteristic_scale"]
        assert metric["relative_max_on_historical_characteristic_scale"] == metric["absolute_max"]/3.86809158318442e-6
        assert metric["phase_amplitude_alignment"] is False
    assert result["interpolation_comparability"]["status"] == "PASS"
    assert result["numerical_robustness_status"] == "PASS"
    assert result["full_period_numerical_robustness_gate"] is True


def _saved_first_row(path, key):
    """Inspect the initial state without loading a whole generated history."""
    with zipfile.ZipFile(path) as archive, archive.open(key+".npy") as stream:
        version = np.lib.format.read_magic(stream)
        header = (np.lib.format.read_array_header_1_0 if version == (1, 0)
                  else np.lib.format.read_array_header_2_0)
        shape, fortran_order, dtype = header(stream)
        assert not fortran_order and len(shape) >= 2
        count = int(np.prod(shape[1:]))
        row = np.frombuffer(stream.read(count*dtype.itemsize), dtype=dtype).reshape(shape[1:])
    return shape, row


def test_actual_full_period_one_d_reuses_all_saved_coordinates_and_zero_velocities(science):
    if not (NEW_BUNDLE / "one_d_all8_spatial.json").exists():
        pytest.skip("Completed local full-period 1D data unavailable")
    static_bundle = ROOT / science["source_static"]["bundle"]
    previous_times = None
    for p in (64, 48):
        with np.load(static_bundle / f"one_d_p{p}.npz") as static:
            for kind in ("linear", "nonlinear"):
                path = NEW_BUNDLE / f"one_d_p{p}_{kind}.npz"
                q_shape, q0 = _saved_first_row(path, "q")
                v_shape, v0 = _saved_first_row(path, "velocity")
                assert q_shape == v_shape == (28499, 4*(p-1))
                np.testing.assert_array_equal(q0, static["q_"+kind])
                np.testing.assert_array_equal(v0, np.zeros_like(q0))
                with np.load(path) as data:
                    times = data["times"]
                    assert times[0] == 0. and times[-1] == 10.37828159055014
                    assert np.all(np.diff(times) > 0)
                    for fraction in (0., .05, .25, .5, .75, 1.):
                        assert fraction*10.37828159055014 in times
                    if previous_times is not None:
                        np.testing.assert_array_equal(times, previous_times)
                    previous_times = times.copy()
        row = _read(NEW_BUNDLE / f"one_d_p{p}_nonlinear.json")
        execution = row["execution"]
        assert execution["status"] == "PASS" and execution["new_ODE_integrations"] == 1
        assert execution["authorization_id"] == "explicit_user_FEM3C_2026_10_09"
        assert execution["initial_coordinates_reused_exactly"] is True
        assert execution["no_dynamic_derivative_constraints"] is True
        assert execution["external_force_after_release"] == 0
        assert execution["rtol"] == 1e-10 and execution["max_step"] == .0072093929877006385
        assert len(execution["atol"]) == 8*(p-1)
        assert execution["counters"]["linear_eigendecompositions"] == 0
        assert execution["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED" and execution["admitted"] is False
        assert execution["strict_float64_strong_weak"] == "PARTIAL"
        diagnostics = row["diagnostics"]
        assert diagnostics["max_relative_energy_drift"] <= 1e-6
        assert diagnostics["removed_GRAV_potential_included"] is False
        assert diagnostics["temporal_convergence_claimed"] is False
        assert diagnostics["spatial_convergence_claimed"] is False
        assert diagnostics["sampled_safety"]["relative_mass_lower_bound"] > 0
        assert diagnostics["sampled_safety"]["min_one_plus_c"] > 0


def test_actual_all8_full_period_qualification_keeps_every_failed_component_visible():
    if not (NEW_BUNDLE / "one_d_all8_spatial.json").exists():
        pytest.skip("Completed local full-period 1D data unavailable")
    result = _read(NEW_BUNDLE / "one_d_all8_spatial.json")
    assert result["status"] == "PARTIAL" and len(result["fields"]) == 8
    assert sum(row["pass"] for row in result["fields"].values()) == 2
    assert result["fields"]["q_u"]["pass"] is True and result["fields"]["q_w"]["pass"] is True
    for name, row in result["fields"].items():
        expected = 1e-4 if name in ("q_w", "q_theta", "velocity_w", "velocity_theta") else 1e-3
        assert row["tolerance"] == expected
        assert row["floor_limited"] is False
        assert row["pass"] is (row["relative_L2"] <= expected and row["relative_max"] <= expected)
        if name not in ("q_u", "q_w"):
            assert row["pass"] is False


def test_actual_one_d_correction_sensitivity_does_not_certify_all8_fields():
    path = NEW_BUNDLE / "one_d_correction_spatial_diagnostic.json"
    if not path.exists():
        pytest.skip("Completed local saved-correction diagnostic unavailable")
    result = _read(path)
    assert result["scientific_calls"] == 0 and result["full_four_field_spatial_status"] == "PARTIAL"
    assert result["actual_samples"] == 28499 and result["degrees"] == [48, 64]
    for label, window in result["windows"].items():
        for component in ("Delta_w", "delta_evol_w"):
            row = window["components"][component]
            assert row["sampled_maxima_only"] is True
            assert row["observed_difference_over_fixed_C1_baseline_model_discrepancy"] == row["absolute_max"]/2.824717e-7
    prefix = result["prefix_reproduction"]
    for kind in ("linear", "nonlinear"):
        assert prefix[kind]["shared_samples"] == 14301
        assert prefix[kind]["exact_old_timestamps_preserved"] is True
        assert prefix[kind]["q0_bitwise_unchanged"] is True and prefix[kind]["v0_bitwise_unchanged"] is True
    assert "correlated" in result["qualification"]


def test_actual_completed_validation_has_exactly_six_native_attempts_and_qualified_status():
    if not (NEW_BUNDLE / "full_period_comparison.json").exists():
        pytest.skip("Completed local full-period FEM-3C evidence unavailable")
    summary = _read(NEW_BUNDLE / "summary.json")
    assert summary["overall"] == "STRAIGHT_ROD_NONLINEAR_3D_FEM_VERIFICATION_COMPLETE_WITH_QUALIFICATIONS"
    assert summary["completed"] is True and summary["universal_nonlinear_validation"] is False
    assert summary["experimental_validation"] is False
    assert len(summary["attempts"]) == summary["job_calls"]["CCX_production"] == 6
    assert summary["job_calls"]["1D_nonlinear_ODE"] == 2
    for key in ("Gmsh", "1D_static", "physical_root_search", "symbolic_derivations"):
        assert summary["job_calls"][key] == 0
    assert summary["statuses"]["NLSP_FEM3C_ROBUSTNESS"] == "PASS"
    assert summary["statuses"]["NLSP_FEM3C_FULL_PERIOD_3D"] == "PASS"
    assert summary["statuses"]["NLSP_FEM3C_ENERGY_DIAGNOSTICS"] == "PARTIAL"
    assert summary["full_period_1d_spatial_qualification"] == "PARTIAL"
    assert summary["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED" and summary["admitted"] is False
    assert summary["strict_float64_qualification"] == "PARTIAL"
    assert summary["numerical_seconds"] <= 24000
    for ordinal, attempt in enumerate(summary["attempts"], 1):
        assert attempt["ordinal"] == ordinal and attempt["status"] == "PASS"
        assert attempt["automatic_retry"] is False
        assert attempt["authorization_id"] == "explicit_user_FEM3C_2026_10_09"


def test_actual_full_medium_cases_have_final_native_output_and_preload_release_evidence():
    if not (NEW_BUNDLE / "full_period_comparison.json").exists():
        pytest.skip("Completed local full-period FEM-3C evidence unavailable")
    for kind in ("linear", "nonlinear"):
        case = NEW_BUNDLE / "cases/full_period_medium" / kind
        result = _read(case / "recovery.json")
        audit = _read(case / "continuation_audit.json")
        job = _read(case / "job.json")
        assert result["status"] == audit["status"] == "PASS"
        assert job["returncode"] == 0 and job["failure"] is None
        assert job["seconds"] < 5400 and job["peak_working_set_bytes"] < 4*1024**3
        assert result["dynamic_time_end"] == 10.37828159055014
        assert result["dynamic_increments"] == 2002 and result["dynamic_output_frames"] == 401
        assert result["output_frequency"] == 5 and result["final_native_frame_reached"] is True
        assert result["cutbacks"] == 0 and audit["warning_lines"] == []
        assert result["maximum_native_external_work_after_release"] == 0
        assert result["maximum_native_damping_work_after_release"] == 0
        assert result["all_saved_frames_strain_diagnostics"]["minimum_det_deformation_gradient"] > 0
        for key in ("node_displacement_max_difference", "support_RF_max_difference",
                    "stress_difference", "strain_difference", "section_profile_max_difference"):
            assert result["preload_transfer"][key] == 0
        with np.load(case / "section_history.npz") as native, np.load(case / "initial_sections.npz") as initial:
            assert native["time"][0] > 0 and native["time"][-1] == 10.37828159055014
            assert len(native["time"]) == 401 and np.all(np.diff(native["time"]) > 0)
            assert bool(initial["not_a_native_dynamic_zero_frame"]) is True


def test_actual_full_comparison_preserves_seven_field_and_uncertainty_qualifications():
    path = NEW_BUNDLE / "full_period_comparison.json"
    if not path.exists():
        pytest.skip("Completed local full-period comparison unavailable")
    result = _read(path)
    assert result["comparison_points"] == 401 and result["section_points"] == 41
    assert set(result["inactive_one_d_fields"]) == {"v", "Phi", "psi"}
    for field in result["inactive_one_d_fields"].values():
        assert field["one_d_identically_zero_by_planar_subspace"] is True
        assert field["out_of_plane_stability_verified"] is False
        assert field["remainder_over_active_physical_scale"] == field["three_d_observed_max"]/field["physical_scale"]
    assert result["one_d_spatial_status"] == "PARTIAL"
    assert result["quarter_period_robustness_not_extended_to_full_period"] is True
    assert result["no_phase_amplitude_frequency_time_fitting"] is True
    for key in ("nonlinear_periodic_orbit_assumed_or_found", "p64_exact_continuum_truth",
                "experimental_validation", "general_seven_field_nonlinear_validation",
                "coupled_rod_joint_validation"):
        assert result[key] is False
    assert result["observations"]["c_eff_diagnostic"]["qualification"] == "3D effective contraction proxy, not identical generalized coordinate"
    for name in ("v", "Phi", "psi"):
        assert result["observations"][name]["one_d_max_at_observation"] == 0
        assert "not stability evidence" in result["observations"][name]["qualification"]
    energy = result["energy"]
    assert energy["status"] == "PARTIAL" and energy["raw_native_data_unchanged"] is True
    assert energy["independent_internal_StVK_energy"] == "NOT_RUN"
    for case in energy["cases"].values():
        assert case["native_energy_static_to_dynamic_reference_jump_relative"] != 0
        kinetic = case["independent_kinetic_energy"]
        assert kinetic["samples"] == 401 and kinetic["status"] == "COMPLETED_DIAGNOSTIC"
        assert "rho" in kinetic["definition"] and "consistent mass" in kinetic["quadrature"]
        assert kinetic["native_internal_energy_reference_corrected"] is False
    for record in result["model_comparisons"].values():
        metric = record["w"]
        assert metric["relative_max"] == metric["absolute_max"]/metric["characteristic_scale"]
        assert metric["sampled_maxima_only"] is True and metric["phase_amplitude_alignment"] is False

