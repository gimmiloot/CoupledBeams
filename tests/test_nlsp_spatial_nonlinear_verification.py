"""Scoped orchestration/source tests; no real numerical or native jobs."""
import copy
from pathlib import Path

import numpy as np
import pytest

from scripts.lib import nlsp_spatial_nonlinear_verification as workflow

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def config():
    return workflow.read(REPO/"data/input/nlsp_spatial_nonlinear_3d_fem_verification.json")


@pytest.mark.parametrize("path,value", [
    (("authorization",), "historical_FEM3C_permission"),
    (("geometry", "h"), .12), (("material", "nu"), .31),
    (("fields",), ["u", "w", "v", "Phi", "theta", "psi", "c"]),
    (("degrees",), [64, 96]), (("horizon_T1",), 1.), (("omega1",), .6),
    (("rotation_vector_local",), ["Phi", "psi", "theta"]),
    (("new_meshes",), True), (("full_period_FEM",), True),
    (("admitted",), True), (("execution_mode",), "STRICT_ADMITTED"),
    (("load_policy", "primary_w_over_h"), .05),
    (("load_policy", "g_k_over_g_n"), 2.),
    (("FEM_dynamic", "alpha"), -.05),
    (("FEM_budget", "automatic_retry"), True),
    (("comparison", "u_c_relative"), .01),
    (("comparison", "bending_rotation_relative"), .001),
    (("comparison", "energy_relative_drift"), 1e-5),
    (("comparison", "planning_signal_ratio"), 1.),
    (("comparison", "phase_amplitude_time_fitting"), True),
])
def test_frozen_config_rejects_unapproved_changes(config, path, value):
    changed = copy.deepcopy(config); target = changed
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        workflow.validate_config(changed)


def test_primary_config_keeps_four_original_cases_and_all14_thresholds(config):
    assert workflow.validate_config(config) is config
    assert config["one_d_cases"] == [["joint_p48",48,"joint"],["joint_p64",64,"joint"],
        ["isolated_w_p64",64,"w"],["isolated_v_p64",64,"v"]]
    assert config["one_d_maximum_nonlinear_ODE_calls"] == 4
    assert config["comparison"]["u_c_relative"] == 1e-3
    assert config["comparison"]["bending_rotation_relative"] == 1e-4
    assert config["comparison"]["relative_numerical_floor"] == 1e-10


@pytest.fixture
def source_fixture(tmp_path, monkeypatch, config):
    monkeypatch.setattr(workflow, "ROOT", tmp_path)
    files = {"action": ["result.json"], "FEM1": ["preflight.json",
        "meshes/medium/solid_mesh.inp", "meshes/medium/mesh_audit.json",
        "meshes/fine/solid_mesh.inp", "meshes/fine/mesh_audit.json"],
        "FEM2": ["one_d_p64.npz", "one_d_preflight.json"],
        "FEM3AR": ["one_d_nonlinear.npz"],
        "FEM3C": ["provenance.json", "one_d_p48_nonlinear.npz", "one_d_p64_nonlinear.npz",
            "cases/medium_refined_time/nonlinear/section_history.npz",
            "cases/fine_refined_time/nonlinear/section_history.npz"],
        "profile_audit": ["summary.json", "numbers.json"]}
    sources = {}
    for name, names in files.items():
        source = tmp_path/config["source_bundles"][name]["path"]
        source.mkdir(parents=True); hashes = {}
        for filename in names:
            target = source/filename; target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text("immutable synthetic " + name + "/" + filename)
            hashes[filename] = workflow.sha(target)
        workflow.write(source/"manifest.json", {"artifact_hashes": hashes})
        config["source_bundles"][name]["manifest_sha256"] = workflow.sha(source/"manifest.json")
        sources[name] = source
    binary = tmp_path/"ccx.exe"; binary.write_bytes(b"not an executable")
    config.update(ccx_exe=str(binary), ccx_sha256=workflow.sha(binary))
    return config, sources


def test_selected_source_registry_uses_hashes_without_replaying_studies(source_fixture):
    config, sources = source_fixture
    found, registry = workflow.source_registry(config)
    assert found == sources and len(registry) == 16
    for name, digest in registry.items():
        assert workflow.sha(workflow.ROOT/name) == digest


def test_corrupt_source_manifest_stops_before_science(source_fixture):
    config, sources = source_fixture
    (sources["FEM3C"]/"manifest.json").write_text("corrupt")
    with pytest.raises(ValueError, match="manifest"):
        workflow.source_registry(config)


def test_corrupt_selected_mesh_or_saved_state_is_not_recomputed(source_fixture):
    config, sources = source_fixture
    (sources["FEM1"]/"meshes/medium/solid_mesh.inp").write_text("changed")
    with pytest.raises(ValueError, match="corrupt immutable source"):
        workflow.source_registry(config)


def cached_fixture(tmp_path, config, sources):
    bundle = tmp_path/"bundle"; bundle.mkdir()
    summary = {"completed": True, "hard_stop": False, "overall": "NUMERICAL_PARTIAL",
        "stage_A": "PASS_WITH_QUALIFICATIONS", "stage_B": "NUMERICAL_PARTIAL", "stage_C": "NOT_RUN",
        "calls": {"nonlinear_1D_ODE": 6, "CCX_production": 0, "Gmsh": 0, "FEM_modal": 0}}
    item = {"config": config, "source_artifacts": {}}
    workflow.write(bundle/"provenance.json", item)
    workflow.save(bundle, item, summary)
    return bundle, item, summary


@pytest.mark.parametrize("mode", ["--report-only", "--plot-only"])
def test_cached_modes_do_not_construct_or_solve(source_fixture, tmp_path, monkeypatch, mode):
    config, sources = source_fixture
    bundle, _, original = cached_fixture(tmp_path, config, sources)
    def forbidden(*args, **kwargs):
        raise AssertionError("Cached output must not execute science")
    monkeypatch.setattr(workflow, "prepare", forbidden)
    monkeypatch.setattr(workflow, "load_discretizations", forbidden)
    monkeypatch.setattr(workflow, "compute_stage_ab", forbidden)
    monkeypatch.setattr(workflow, "compute_additional_controls", forbidden)
    from scripts.lib import nlsp_spatial_comparison as comparison
    monkeypatch.setattr(comparison, "render_figures", lambda path: [])
    result = workflow.main([mode, str(bundle)])
    assert result == original


def test_cache_rejects_changed_own_artifact(source_fixture, tmp_path):
    config, sources = source_fixture
    bundle, _, _ = cached_fixture(tmp_path, config, sources)
    (bundle/"summary.json").write_text("modified")
    with pytest.raises(ValueError, match="cache artifact"):
        workflow.validate_cache(bundle)


def test_root_manifest_hashes_nested_immutable_manifest_and_detects_tamper(source_fixture, tmp_path):
    config,sources = source_fixture
    bundle,item,summary = cached_fixture(tmp_path,config,sources)
    child = bundle/"recovery_sensitivity_endpoints/medium"; child.mkdir(parents=True)
    data = child/"profiles.npz"; data.write_bytes(b"saved synthetic profile array")
    nested = child/"manifest.json"
    workflow.write(nested,{"artifact_hashes":{"profiles.npz":workflow.sha(data)}})
    original = nested.read_bytes()
    workflow.save(bundle,item,summary)
    hashes = workflow.read(bundle/"manifest.json")["artifact_hashes"]
    assert "manifest.json" not in hashes
    assert hashes["recovery_sensitivity_endpoints/medium/manifest.json"] == workflow.sha(nested)
    assert nested.read_bytes() == original
    assert workflow.validate_cache(bundle) == summary
    workflow.write(nested,{"artifact_hashes":{"profiles.npz":"changed metadata only"}})
    assert data.read_bytes() == b"saved synthetic profile array"
    with pytest.raises(ValueError,match="cache artifact.*manifest.json"):
        workflow.validate_cache(bundle)


def test_existing_cli_dispatches_scoped_flag_only(monkeypatch):
    from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as cli
    called = []
    monkeypatch.setattr(workflow, "main", lambda argv: called.append(argv) or "dispatched")
    assert cli.main(["--spatial-verification", "--preflight"]) == "dispatched"
    assert called == [["--preflight"]]


def test_existing_cli_historical_routes_remain_separate(monkeypatch, tmp_path):
    from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as cli
    monkeypatch.setattr(workflow, "main", lambda *args: pytest.fail("Old route entered new experiment"))
    result = {"resume_statuses": {}, "overall": "PILOT_COMPLETE_WITH_QUALIFICATIONS"}
    monkeypatch.setattr(cli, "validate_cache", lambda bundle: result)
    assert cli.main(["--report-only", str(tmp_path)]) == result


@pytest.fixture
def additional_fixture(tmp_path, monkeypatch, config):
    bundle = tmp_path/"additional"; bundle.mkdir()
    item = {"config": config}
    summary = {"completed": False, "hard_stop": False, "stage_B": "NUMERICAL_PARTIAL",
        "overall": "ONE_D_PREFLIGHT_COMPLETE", "one_d_seconds": 437.,
        "one_d_cases": {name: {"status": "PASS", "p": p} for name,p,_ in config["one_d_cases"]},
        "calls": {"nonlinear_1D_ODE": 4}}
    workflow.write(bundle/"config.json", config)
    workflow.write(bundle/"frozen_load.json", {"q_w": 2.3e-5, "q_v": 2.875e-5})
    workflow.write(bundle/"summary.json", summary)
    workflow.write(bundle/"stage_b_comparison.json", {"status": "NUMERICAL_PARTIAL", "stage_c_allowed": False})
    workflow.write(bundle/"manifest.json", {"preserved_primary": True})
    np.save(bundle/"one_d_common_times.npy", np.array([0., 1., 2.594570397637535]))
    original = {name: (bundle/name).read_bytes() for name in ("config.json", "frozen_load.json",
        "summary.json", "stage_b_comparison.json", "manifest.json")}
    for name,_,_ in config["one_d_cases"]:
        target = bundle/"one_d"/name; target.mkdir(parents=True)
        workflow.write(target/"case.json", {"calls": {"nonlinear_ODE": 1}})
    monkeypatch.setattr(workflow, "load_discretizations", lambda c: ({48: "synthetic p48"}, None, None))
    monkeypatch.setattr(workflow, "save", lambda b, i, s: workflow.write(b/"summary.json", s))
    return bundle, item, summary, original


def test_additional_controls_keep_original_phase_and_only_two_same_excitations(additional_fixture, monkeypatch):
    bundle, item, summary, original = additional_fixture
    from scripts.lib import nlsp_spatial_1d_program as program
    from scripts.lib import nlsp_spatial_comparison as comparison
    calls = []
    def run(disc, loads, times, deadline, target, **kwargs):
        calls.append((disc, loads, times.copy(), kwargs))
        target.mkdir(parents=True)
        workflow.write(target/"case.json", {"calls": {"nonlinear_ODE": 1}})
        return {"status": "PASS"}
    monkeypatch.setattr(program, "run_case", run)
    monkeypatch.setattr(comparison, "analyze_stage_b", lambda *a, **k: {"status": "NUMERICAL_PARTIAL"})
    result = workflow.compute_additional_controls(bundle,item,summary)
    assert result["stage_B"] == "NUMERICAL_PARTIAL" and result["calls"]["nonlinear_1D_ODE"] == 6
    assert len(calls) == 2 and calls[0][1] == (2.3e-5, 0.) and calls[1][1] == (0., 2.875e-5)
    assert all(call[0] == "synthetic p48" for call in calls)
    assert (bundle/"config.json").read_bytes() == original["config.json"]
    assert (bundle/"frozen_load.json").read_bytes() == original["frozen_load.json"]
    baseline = bundle/"primary_four_case_evidence"
    assert (baseline/"phase_manifest.json").read_bytes() == original["manifest.json"]
    assert (baseline/"summary.json").read_bytes() == original["summary.json"]
    assert (baseline/"stage_b_comparison.json").read_bytes() == original["stage_b_comparison.json"]
    evidence = workflow.read(bundle/"additional_1d_controls_config.json")
    assert evidence["maximum_additional_nonlinear_ODE_calls"] == 2 and not evidence["automatic_retry"]
    assert evidence["full_field_and_velocity_thresholds_unchanged"]


def test_additional_failure_stops_before_second_control_and_no_retry(additional_fixture, monkeypatch):
    bundle,item,summary,_ = additional_fixture
    from scripts.lib import nlsp_spatial_1d_program as program
    calls = []
    monkeypatch.setattr(program, "run_case", lambda *a, **k: calls.append(1) or {"status": "FAIL", "failure": "synthetic failure"})
    result = workflow.compute_additional_controls(bundle,item,summary)
    assert result["hard_stop"] and len(calls) == 1
    workflow.compute_additional_controls(bundle,item,summary)
    assert len(calls) == 1 and "isolated_v_p48" not in result["one_d_cases"]


def test_completed_additional_cache_has_zero_new_calls(monkeypatch):
    summary = {"completed": False, "hard_stop": False, "stage_B": "NUMERICAL_PARTIAL", "additional_controls": "PASS"}
    monkeypatch.setattr(workflow, "load_discretizations", lambda *a: pytest.fail("Cached controls constructed solver"))
    assert workflow.compute_additional_controls(Path("unused"),{},summary) is summary


def frozen_decision_fixture(bundle, item, summary):
    baseline = bundle/"primary_four_case_evidence"; baseline.mkdir(exist_ok=True)
    workflow.write(baseline/"stage_b_comparison.json", {"stage_c_allowed": False, "status": "NUMERICAL_PARTIAL"})
    workflow.write(bundle/"stage_a_checks.json", {"execution_gate": "PASS"})
    workflow.write(bundle/"stage_b_comparison.json", {"stage_c_allowed": True,
        "stage_c_decision": {"limited_evolving_bending_and_mixed_displacement_resolved": True},
        "all14_spatial": {"status": "PARTIAL", "velocity": {"status": "PARTIAL"}}})
    workflow.write(bundle/"additional_1d_controls_config.json", {"maximum_additional_nonlinear_ODE_calls": 2})
    item["source_artifacts"] = {}
    summary["additional_controls"] = "PASS"


def test_frozen_fem_decision_keeps_full14_partial_and_original_rule(additional_fixture, monkeypatch):
    bundle,item,summary,_ = additional_fixture
    frozen_decision_fixture(bundle,item,summary)
    monkeypatch.setattr(workflow.shutil, "disk_usage", lambda p: type("Disk", (), {"free": 20*1024**3})())
    result = workflow.freeze_fem_decision(bundle,item,summary)
    assert result["stage_C_allowed"] and result["selected_before_any_new_3D_result"]
    assert result["full14_spatial_status"] == "PARTIAL"
    assert result["full_velocity_spatial_status"] == "PARTIAL"
    assert result["stage_B_overall_status_unchanged"] == "NUMERICAL_PARTIAL"
    assert result["maximum_new_production_CCX_jobs"] == 4
    original = (bundle/"pre_fem_decision.json").read_bytes()
    workflow.write(bundle/"stage_b_comparison.json", {"stage_c_allowed": False})
    assert workflow.freeze_fem_decision(bundle,item,summary) == result
    assert (bundle/"pre_fem_decision.json").read_bytes() == original


def test_insufficient_output_disk_blocks_before_native(additional_fixture, monkeypatch):
    bundle,item,summary,_ = additional_fixture
    frozen_decision_fixture(bundle,item,summary)
    monkeypatch.setattr(workflow.shutil, "disk_usage", lambda p: type("Disk", (), {"free": 1024**3})())
    result = workflow.freeze_fem_decision(bundle,item,summary)
    assert not result["stage_C_allowed"] and "disk" in result["stop_reason"]


@pytest.fixture
def native_stage_fixture(additional_fixture, monkeypatch, tmp_path):
    bundle,item,summary,_ = additional_fixture
    item["source_artifacts"] = {}
    summary.update(additional_controls="PASS", stage_C="NOT_RUN", FEM_cases={}, attempts=[], CCX_seconds=0.)
    summary["calls"].update(CCX_production=0, Gmsh=0, FEM_modal=0)
    workflow.write(bundle/"stage_b_comparison.json", {"status": "NUMERICAL_PARTIAL"})
    decision = {"stage_C_allowed": True, "planning_runtime_seconds":
        {"medium_linear": 1122., "medium_nonlinear": 1185., "fine_linear": 2922., "fine_nonlinear": 3071.}}
    workflow.write(bundle/"pre_fem_decision.json", decision)
    monkeypatch.setattr(workflow, "freeze_fem_decision", lambda *args: decision)
    source = tmp_path/"FEM1"
    for level in ("medium", "fine"):
        directory = source/"meshes"/level; directory.mkdir(parents=True)
        (directory/"solid_mesh.inp").write_text("existing immutable " + level)
        workflow.write(directory/"mesh_audit.json", {"status": "PASS"})
    monkeypatch.setattr(workflow, "source_registry", lambda config: ({"FEM1": source}, {}))
    monkeypatch.setattr(workflow, "load_discretizations", lambda config: ({64: "synthetic p64"}, None, None))
    from scripts.lib import nlsp_spatial_fem_protocol as protocol
    from scripts.lib import nlsp_spatial_comparison as comparison
    monkeypatch.setattr(protocol.fem1.single, "read_gmsh_inp_mesh_data", lambda path: "existing synthetic mesh")
    informative = {"value": True}
    def compare(first, second, one, disc, output, **kwargs):
        output.mkdir()
        workflow.write(output/"one_d_three_d_comparison.json", {"signal_planning":
            {field: {"resolved_for_pre_FEM_planning": informative["value"]} for field in ("w", "v")}})
        return {}
    monkeypatch.setattr(comparison, "compare_fem_pair", compare)
    return bundle,item,summary,informative


def test_native_order_exact_four_and_cached_completion_zero_calls(native_stage_fixture, monkeypatch):
    bundle,item,summary,_ = native_stage_fixture
    from scripts.lib import nlsp_spatial_native_program as native
    calls = []
    def run(case, source, mesh, audit, config, load, *, authorization, remainingBudget):
        calls.append(authorization.copy()); case.mkdir(parents=True)
        workflow.write(case/"native_attempt.json", {"solver_calls": 1})
        assert source.name == authorization["mesh_level"]
        assert authorization["stage_A_pass"] and authorization["stage_B_pass"]
        assert (bundle/"pre_fem_decision.json").is_file()
        assert remainingBudget > 0
        return {"status": "PASS", "new_solver_calls": 1, "native_seconds": 10., "recovery_seconds": 1.}
    monkeypatch.setattr(native, "run_attempt", run)
    result = workflow.compute_stage_c(bundle,item,summary)
    assert [(call["mesh_level"], call["nonlinear"], call["ordinal"]) for call in calls] == [
        ("medium",False,1),("medium",True,2),("fine",False,3),("fine",True,4)]
    assert result["calls"]["CCX_production"] == 4 and result["completed"]
    assert result["overall"] == "NUMERICAL_PARTIAL" and result["stage_B"] == "NUMERICAL_PARTIAL"
    monkeypatch.setattr(workflow, "load_discretizations", lambda *args: pytest.fail("Cached FEM constructed 1D solver"))
    assert workflow.compute_stage_c(bundle,item,summary) is result and len(calls) == 4


@pytest.mark.parametrize("status", ["FAIL", "OUTPUT_RECOVERY_PENDING"])
def test_first_native_failure_or_pending_parser_stops_and_never_retries(native_stage_fixture, monkeypatch, status):
    bundle,item,summary,_ = native_stage_fixture
    from scripts.lib import nlsp_spatial_native_program as native
    calls = []
    def run(case, *args, **kwargs):
        calls.append(1); case.mkdir(parents=True)
        workflow.write(case/"native_attempt.json", {"solver_calls": 1})
        return {"status": status, "new_solver_calls": 1, "native_seconds": 10., "recovery_seconds": 0.,
            "native_failure": "synthetic stopped prefix"}
    monkeypatch.setattr(native, "run_attempt", run)
    result = workflow.compute_stage_c(bundle,item,summary)
    assert result["hard_stop"] and len(calls) == 1
    assert set(result["FEM_cases"]) == {"medium_linear"}
    workflow.compute_stage_c(bundle,item,summary)
    assert len(calls) == 1


def test_uninformative_medium_pair_forbids_fine(native_stage_fixture, monkeypatch):
    bundle,item,summary,informative = native_stage_fixture
    informative["value"] = False
    from scripts.lib import nlsp_spatial_native_program as native
    calls = []
    def run(case, *args, **kwargs):
        calls.append(1); case.mkdir(parents=True)
        workflow.write(case/"native_attempt.json", {"solver_calls": 1})
        return {"status": "PASS", "new_solver_calls": 1, "native_seconds": 10., "recovery_seconds": 1.}
    monkeypatch.setattr(native, "run_attempt", run)
    result = workflow.compute_stage_c(bundle,item,summary)
    assert len(calls) == 2 and result["completed"]
    assert result["overall"] == "SPATIAL_NONLINEAR_SIGNAL_NOT_RESOLVED"
    assert set(result["FEM_cases"]) == {"medium_linear", "medium_nonlinear"}


def test_failing_prefem_signal_gate_prevents_all_native_calls(native_stage_fixture, monkeypatch):
    bundle,item,summary,_ = native_stage_fixture
    monkeypatch.setattr(workflow, "freeze_fem_decision", lambda *a: {"stage_C_allowed": False})
    from scripts.lib import nlsp_spatial_native_program as native
    monkeypatch.setattr(native, "run_attempt", lambda *a, **k: pytest.fail("Unresolved signal entered FEM"))
    result = workflow.compute_stage_c(bundle,item,summary)
    assert result["stage_C"] == "NOT_RUN" and result["calls"]["CCX_production"] == 0


def test_postprocess_unfinished_prefix_rejected_before_model_load(tmp_path, monkeypatch):
    summary = {"completed": False, "hard_stop": False}
    monkeypatch.setattr(workflow, "load_discretizations", lambda *args: pytest.fail("Active prefix constructed model"))
    with pytest.raises(ValueError, match="prefix"):
        workflow.postprocess_saved(tmp_path, {}, summary)
    assert not (tmp_path/"execution_code/final_postprocessing").exists()


def test_finished_postprocess_preserves_first_metrics_source_code_and_zero_science(additional_fixture, monkeypatch):
    bundle,item,summary,_ = additional_fixture
    from scripts.lib import nlsp_spatial_comparison as comparison
    from scripts.lib import nlsp_spatial_native_program as native
    from scripts.lib import nlsp_spatial_1d_program as program
    def forbidden(*args, **kwargs):
        pytest.fail("Saved postprocessing launched physical science")
    monkeypatch.setattr(native, "run_attempt", forbidden)
    monkeypatch.setattr(program, "run_case", forbidden)
    monkeypatch.setattr(workflow, "compute_stage_ab", forbidden)
    monkeypatch.setattr(workflow, "compute_stage_c", forbidden)
    monkeypatch.setattr(workflow, "load_discretizations", lambda *args: ({64: "saved reference evaluator"}, None, None))
    summary.update(completed=True, FEM_cases={"medium_linear": {"status": "PASS"},
        "medium_nonlinear": {"status": "PASS"}}, calls={"nonlinear_1D_ODE": 6, "CCX_production": 2})
    output = bundle/"comparison_medium"; output.mkdir()
    old = b'{"first_comparison": true}\n'
    (output/"one_d_three_d_comparison.json").write_bytes(old)
    old_csv = b"time,old_value\n0,1\n"; (output/"old_metrics.csv").write_bytes(old_csv)
    calls = []
    def compare(*args, **kwargs):
        destination = args[4]
        calls.append(destination)
        workflow.write(destination/"one_d_three_d_comparison.json", {"added_chi_diagnostics": True})
    monkeypatch.setattr(comparison, "compare_fem_pair", compare)
    monkeypatch.setattr(comparison, "analyze_mesh_comparison", lambda *args: {})
    monkeypatch.setattr(comparison, "write_initial_state_decomposition", lambda *args, **kwargs: {})
    support_audits = []
    monkeypatch.setattr(comparison, "audit_static_torsional_supports",
        lambda *args, **kwargs: support_audits.append((args, kwargs)) or {})
    output_tables = []
    monkeypatch.setattr(workflow, "write_output_tables",
        lambda *args: output_tables.append(args) or {})
    from scripts.lib import nlsp_spatial_recovery_endpoints as endpoints
    endpoint_audits = []
    def audit_endpoint(*args):
        endpoint_audits.append(args)
        assert (output/"first_processing_evidence/one_d_three_d_comparison.json").read_bytes() == old
        return {"status":"DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS","recovery_calls_81":4,
                "scientific_calls":{"CCX":0,"Radau":0}}
    monkeypatch.setattr(endpoints,"audit_recovery_endpoints",audit_endpoint)
    monkeypatch.setattr(comparison, "render_figures", lambda *args: ["saved_figure.png"])
    original_calls = copy.deepcopy(summary["calls"])
    result = workflow.postprocess_saved(bundle,item,summary)
    assert result["postprocessing_complete"] and result["calls"] == original_calls
    assert endpoint_audits == [(bundle,"medium")]
    assert set(result["recovery_endpoint_sensitivity"]) == {"medium"}
    assert result["recovery_endpoint_sensitivity"]["medium"]["recorded_81_section_recoveries"] == 4
    assert len(output_tables) == 1 and output_tables[0] == (bundle, summary)
    assert len(support_audits) == 1
    assert support_audits[0][0][0] == bundle
    assert support_audits[0][1]["discs"] == {64: "saved reference evaluator"}
    assert support_audits[0][1]["action_result"] == REPO/item["config"]["source_bundles"]["action"]["path"]/"result.json"
    assert all(count == 0 for count in result["postprocessing_scientific_calls"].values())
    backup = output/"first_processing_evidence"
    assert (backup/"one_d_three_d_comparison.json").read_bytes() == old
    assert (backup/"old_metrics.csv").read_bytes() == old_csv
    phase = bundle/"execution_code/final_postprocessing"
    for name in ("nlsp_spatial_nonlinear_verification", "nlsp_spatial_comparison"):
        assert (phase/(name+".py")).read_bytes() == (REPO/"scripts/lib"/(name+".py")).read_bytes()
    # Matching replay cannot overwrite the preserved first evidence or snapshot.
    snapshot = (phase/"nlsp_spatial_comparison.py").read_bytes()
    monkeypatch.setattr(workflow, "load_discretizations", forbidden)
    monkeypatch.setattr(comparison, "render_figures", forbidden)
    assert workflow.postprocess_saved(bundle,item,summary) is result and len(calls) == 1
    assert len(support_audits) == 1
    assert len(output_tables) == 1
    assert endpoint_audits == [(bundle,"medium")]
    assert (backup/"one_d_three_d_comparison.json").read_bytes() == old
    assert (phase/"nlsp_spatial_comparison.py").read_bytes() == snapshot


def test_completed_matching_plot_preserves_manifest_and_image_bytes(source_fixture, tmp_path, monkeypatch):
    config,sources = source_fixture
    bundle,item,summary = cached_fixture(tmp_path,config,sources)
    workflow.write(bundle/"figure_data_provenance.json", {"completed_render": True})
    image = bundle/"figure.png"; image.write_bytes(b"saved synthetic PNG")
    workflow.save(bundle,item,summary)
    manifest = (bundle/"manifest.json").read_bytes(); raw_image = image.read_bytes()
    from scripts.lib import nlsp_spatial_comparison as comparison
    def forbidden(*args, **kwargs):
        pytest.fail("Matching plot rewrote completed evidence")
    monkeypatch.setattr(comparison, "render_figures", forbidden)
    monkeypatch.setattr(workflow, "save", forbidden)
    monkeypatch.setattr(workflow, "load_discretizations", forbidden)
    assert workflow.main(["--plot-only",str(bundle)]) == summary
    assert (bundle/"manifest.json").read_bytes() == manifest
    assert image.read_bytes() == raw_image


def test_matching_postprocess_cli_cache_does_not_construct_or_regenerate(source_fixture, tmp_path, monkeypatch):
    config,sources = source_fixture
    bundle,item,summary = cached_fixture(tmp_path,config,sources)
    summary["postprocessing_complete"] = True
    workflow.save(bundle,item,summary)
    manifest = (bundle/"manifest.json").read_bytes()
    monkeypatch.setattr(workflow, "load_discretizations", lambda *args: pytest.fail("Cached postprocess constructed model"))
    monkeypatch.setattr(workflow, "save", lambda *args: pytest.fail("Cached postprocess rewrote manifest"))
    assert workflow.main(["--postprocess-only",str(bundle)]) == summary
    assert (bundle/"manifest.json").read_bytes() == manifest


def test_saved_initial_decomposition_is_not_reevaluated(additional_fixture, monkeypatch):
    bundle,item,summary,_ = additional_fixture
    summary.update(completed=True, FEM_cases={})
    from scripts.lib import nlsp_spatial_comparison as comparison
    evidence = bundle/"initial_state_decomposition.json"
    evidence.write_bytes(b'{"saved_decomposition": true}\n')
    original = evidence.read_bytes()
    monkeypatch.setattr(workflow, "load_discretizations", lambda *args: ({64: "saved evaluator"}, None, None))
    monkeypatch.setattr(comparison, "analyze_mesh_comparison", lambda *args: {})
    monkeypatch.setattr(comparison, "render_figures", lambda *args: [])
    monkeypatch.setattr(comparison, "write_initial_state_decomposition",
        lambda *args, **kwargs: pytest.fail("Saved decomposition was reevaluated"))
    monkeypatch.setattr(comparison, "audit_static_torsional_supports", lambda *args, **kwargs: {})
    result = workflow.postprocess_saved(bundle,item,summary)
    assert result["postprocessing_complete"]
    assert evidence.read_bytes() == original
    assert all(value == 0 for value in result["postprocessing_scientific_calls"].values())


def kinetic_table_fixture(tmp_path, *, native_K=(.4, .1), integrated_K=(.1001, .3999),
                          independent_times=(.1, .2), section_times=(.1, .2), initial_energy=.5):
    case = tmp_path/"FEM/medium_linear"; case.mkdir(parents=True)
    workflow.write(case/"recovery.json", {"native_initial_internal_energy": initial_energy})
    records = [{"step": 1, "increment": 1, "kinetic_energy": 99.}]
    for increment, value in zip((14, 2), native_K):
        row = {"step": 2, "increment": increment}
        if value is not None:
            row["kinetic_energy"] = value
        records.append(row)
    workflow.write(case/"energy.json", {"records": records})
    np.savez_compressed(case/"independent_kinetic_energy.npz", time=independent_times,
                        kinetic_energy=integrated_K)
    np.savez_compressed(case/"section_history.npz", time=section_times, increments=(2, 14))
    summary = {"FEM_cases": {"medium_linear": {"status": "PASS"}}}
    return case, summary


def test_saved_kinetic_table_matches_increment_ids_not_record_order(tmp_path):
    case, summary = kinetic_table_fixture(tmp_path)
    result = workflow.write_output_tables(tmp_path,summary)["medium_linear"]
    assert result["max_absolute_kinetic_difference"] == pytest.approx(.0001)
    assert result["common_full_horizon_kinetic_scale"] == .4
    assert result["relative_kinetic_difference_on_common_scale"] == pytest.approx(.00025)
    assert result["relative_difference_to_own_STATIC_internal_energy"] == pytest.approx(.0002)
    assert result["energy_status"] == "PARTIAL" and result["independent_internal_energy"] == "NOT_RUN"
    assert result["native_reference_jump_not_corrected"]
    import csv
    with (case/"kinetic_comparison.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert [int(row["increment"]) for row in rows] == [2, 14]
    assert [float(row["native_K"]) for row in rows] == [.1, .4]


def test_desynchronized_saved_kinetic_timestamps_are_rejected(tmp_path):
    case, summary = kinetic_table_fixture(tmp_path, independent_times=(.1, .21))
    with pytest.raises(ValueError, match="timestamps differ"):
        workflow.write_output_tables(tmp_path,summary)
    assert not (case/"kinetic_comparison.csv").exists()


def test_missing_native_kinetic_increment_stays_partial(tmp_path):
    case, summary = kinetic_table_fixture(tmp_path, native_K=(None, .1))
    result = workflow.write_output_tables(tmp_path,summary)["medium_linear"]
    assert result["energy_status"] == "PARTIAL"
    assert result["missing_native_kinetic_increment_ids"] == [14]
    assert result["independent_internal_energy"] == "NOT_RUN"
    assert not (case/"kinetic_comparison.csv").exists()


@pytest.mark.parametrize("initial", [0., None])
def test_zero_kinetic_scale_and_missing_or_zero_initial_energy_never_divide_by_zero(tmp_path, initial):
    _, summary = kinetic_table_fixture(tmp_path, native_K=(0., 0.), integrated_K=(0., 0.),
                                      initial_energy=initial)
    result = workflow.write_output_tables(tmp_path,summary)["medium_linear"]
    assert result["max_absolute_kinetic_difference"] == 0.
    assert result["relative_kinetic_difference_on_common_scale"] is None
    assert result["relative_difference_to_own_STATIC_internal_energy"] is None
    assert result["energy_status"] == "PARTIAL"
