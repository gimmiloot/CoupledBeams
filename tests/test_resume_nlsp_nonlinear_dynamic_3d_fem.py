"""FEM-3AR authorization, output safety and saved-data contracts; zero real jobs."""
import copy
from pathlib import Path

import numpy as np
import pytest

from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as cli

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_pilot/a69310e3bb30bab7"
PARENT_SHA = "bf217d9c77d42170d1021b55c61873f8b2a94ed1fa67950d5ff9ed82c5cd257f"


@pytest.fixture
def config():
    return cli.read_json(cli.CONFIG)


@pytest.fixture(scope="module")
def science():
    return cli.read_json(PARENT / "provenance.json")["config"]


@pytest.fixture(scope="module")
def decks():
    return {
        "linear": (PARENT / "remediation_preview/linear_corrected_NOT_RUN.inp").read_text(encoding="utf8"),
        "nonlinear": (PARENT / "input_gate/nonlinear.inp").read_text(encoding="utf8"),
    }


def _forbid_science(monkeypatch):
    def forbidden(*a, **k):
        pytest.fail("No production FEM, ODE, equilibrium, mesh, eigen or symbolic work in this test")
    monkeypatch.setattr(cli.base, "run_case", forbidden)
    monkeypatch.setattr(cli.base, "finish_references_and_comparison", forbidden)
    monkeypatch.setattr(cli.base.one, "integrate_nonlinear_reference", forbidden)
    monkeypatch.setattr(cli.base.one, "exact_linear_reference", forbidden)
    monkeypatch.setattr(cli.base.one.dynamics.PlanarGalerkin, "linear_eigenpairs", forbidden)
    monkeypatch.setattr(cli.base.one.rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(cli.base.base, "fem2_static_newton", forbidden)
    monkeypatch.setattr(cli.base.base.fem1, "run_job", forbidden)
    monkeypatch.setattr(cli.base.base.fem1.single, "generate_mesh_with_gmsh_python", forbidden)
    monkeypatch.setattr(cli.base.base.fem1.single, "generate_mesh_with_gmsh_cli", forbidden)


def _summary():
    return {"overall": "NOT_RUN", "cases": {}, "attempts": [], "statuses": {},
            "resume_statuses": {}, "job_calls": {}, "numerical_seconds": 0.}


def _cached(monkeypatch, tmp_path, science, summary):
    monkeypatch.setattr(cli, "existing_attempt", lambda c: (tmp_path, {"config": science}, summary))
    monkeypatch.setattr(cli, "save", lambda *a: None)
    monkeypatch.setattr(cli, "plot_bundle", lambda *a: {"figures": 0, "new_scientific_calls": 0})


def test_separate_authorization_and_frozen_resource_contract(config):
    assert cli.validate_config(config) is config
    assert config["authorization"]["id"] == "explicit_user_FEM3AR_2026_10_09"
    assert config["authorization"]["id"] != cli.read_json(PARENT / "summary.json")["authorization"]["id"]
    assert config["authorization"]["case_order"] == ["linear", "nonlinear"]
    assert config["authorization"]["maximum_production_CCX_jobs"] == 2
    assert config["authorization"]["maximum_nonlinear_1D_integrations"] == 1
    assert config["authorization"]["automatic_retry"] is False
    assert config["numerical_budget_seconds"] == 3600
    assert config["job_timeout_seconds"] == 1200
    assert config["job_memory_limit_bytes"] == 4 * 1024**3
    assert config["threads"] == 1


@pytest.mark.parametrize("key,value", [
    ("schema", "historical-pilot"), ("threads", 2), ("job_timeout_seconds", 1201),
    ("job_memory_limit_bytes", 8 * 1024**3), ("numerical_budget_seconds", 3601),
    ("science_policy", "rebuild_initial_conditions"), ("new_meshes", True),
    ("new_modal_jobs", True), ("new_static_only_jobs", True), ("new_time_or_space_levels", True)])
def test_unauthorized_extensions_rejected(config, key, value):
    config[key] = value
    with pytest.raises(ValueError):
        cli.validate_config(config)


@pytest.mark.parametrize("key,value", [("id", "explicit_user_FEM3A_2026_10_09"),
    ("maximum_production_CCX_jobs", 3), ("maximum_nonlinear_1D_integrations", 2),
    ("case_order", ["nonlinear", "linear"]), ("automatic_retry", True),
    ("basis", "renamed old attempt")])
def test_authorization_cannot_be_reused_or_expanded(config, key, value):
    config["authorization"][key] = value
    with pytest.raises(ValueError, match="authorization"):
        cli.validate_config(config)


def test_parent_failed_evidence_immutable_and_separate(config):
    assert cli.sha(PARENT / "manifest.json") == PARENT_SHA == config["parent_failed"]["manifest_sha256"]
    assert cli.sha(PARENT / "cases/linear/motion.inp") == "a38fb741b7638103bedee3893280266d32053cdae1b50d66804dd99369650998"
    old = cli.read_json(PARENT / "summary.json")
    assert old["overall"] == "BLOCKED_BY_SOLVER"
    assert old["cases"]["linear"]["status"] == "FAIL"
    assert old["job_calls"]["CCX_production"] == 1
    assert old["job_calls"]["1D_nonlinear_ODE"] == 0
    assert not (PARENT / "one_d_nonlinear.npz").exists()
    assert not (PARENT / "figures").exists()
    assert (PARENT / "cases/linear/motion.dat").stat().st_size == 0
    assert (PARENT / "cases/linear/motion.sta").stat().st_size == 0


def test_parent_load_verifies_original_artifacts_and_science(config, science, monkeypatch):
    _forbid_science(monkeypatch)
    before = cli.sha(PARENT / "manifest.json")
    path, identity, failed, inherited = cli.load_parent(config)
    assert path == PARENT
    assert inherited == science == identity["config"]
    assert failed["overall"] == "BLOCKED_BY_SOLVER"
    assert cli.sha(PARENT / "manifest.json") == before
    assert cli.sha(cli.base.__file__) == config["corrected_generator_sha256"]


@pytest.mark.parametrize("key", ["parent", "generator"])
def test_corrupt_source_stops_before_scientific_calls(config, monkeypatch, key):
    _forbid_science(monkeypatch)
    if key == "parent":
        config["parent_failed"]["manifest_sha256"] = "0" * 64
    else:
        config["corrected_generator_sha256"] = "0" * 64
    with pytest.raises(ValueError, match="changed"):
        cli.load_parent(config)


def test_saved_p64_initial_states_reused_without_new_projection(science):
    pre = cli.read_json(PARENT / "one_d_preflight.json")
    with np.load(ROOT / science["source_static"]["bundle"] / "one_d_p64.npz") as z:
        for kind in ("linear", "nonlinear"):
            row = pre["initial_states"][kind]
            np.testing.assert_array_equal(row["q0"], z["q_" + kind])
            assert len(row["q0"]) == 4 * 63
            assert np.count_nonzero(row["v0"]) == 0
            assert row["midspan_w_acceleration"] == pytest.approx(-.001418192734439, rel=1e-12)
            assert row["midspan_w_acceleration"] < 0
    assert pre["admitted"] is False
    assert pre["strict_float64_strong_weak"] == "PARTIAL"
    assert pre["slope_constraints"] is False


@pytest.mark.parametrize("kind", ["linear", "nonlinear"])
def test_existing_decks_preserve_physics_release_and_time(science, decks, kind):
    text = decks[kind]
    gate = cli.base.input_contract(text, science)
    assert gate["status"] == "PASS"
    assert cli.output_safety(text)["status"] == "PASS"
    assert text.count("*STATIC") == text.count("*DYNAMIC, ALPHA=0") == 1
    assert text.count("*END STEP") == 2
    assert "*DLOAD, OP=NEW\nSOLID,GRAV,0.,0.,-1.,0." in text
    assert "AMPLITUDE=STEP" in text
    assert "ALL_NODES,1,0.\nALL_NODES,2,0.\nALL_NODES,3,0." in text
    assert "LEFT_FIXED,1,3" in text and "RIGHT_FIXED,1,3" in text
    assert "*INCLUDE" in text and "*RESTART" not in text
    assert all(len(token) <= 20 for token in gate["dynamic_native_numeric_fields"])
    assert science["geometry"] == {"L": 1., "b": .2, "h": .1}
    assert science["g"] == .0014224751066856333
    assert science["q"] == 2.844950213371267e-5
    d = cli.base.dynamic_settings(science)
    assert d["duration"] == .05 * d["T1"]
    assert d["T1"] == pytest.approx(10.37828159055014, rel=1e-15)
    assert d["initial_increment"] == d["T1"] / 4000
    assert d["maximum_increment"] == d["T1"] / 2000


def test_only_documented_linear_preload_output_changed(decks):
    old = (PARENT / "cases/linear/motion.inp").read_text(encoding="utf8")
    corrected = decks["linear"]
    strip = lambda text: "\n".join(row for row in text.splitlines() if not row.startswith("*INCLUDE"))
    assert strip(old.replace("ELSE,ELKE\n*END STEP", "ELSE\n*END STEP", 1)) == strip(corrected)
    first, dynamic = corrected.split("*END STEP", 1)
    assert "ELKE" not in first and "S,E,ENER" in first and "\nELSE\n" in first
    assert "ELSE,ELKE" in dynamic


@pytest.mark.parametrize("routing", ["", ", NLGEOM=NO"])
def test_unsafe_linear_static_elke_rejected_in_both_routes(decks, routing):
    text = decks["linear"].replace("*STEP", "*STEP" + routing, 1)
    text = text.replace("\nELSE\n*END STEP", "\nELSE,ELKE\n*END STEP", 1)
    with pytest.raises(ValueError, match="Unsafe linear STATIC ELKE"):
        cli.output_safety(text)


@pytest.mark.parametrize("routing", [", NLGEOM", ", NLGEOM=YES"])
def test_nonlinear_static_route_distinguished_explicitly(decks, routing):
    text = decks["nonlinear"].replace(", NLGEOM", routing, 1)
    assert cli.output_safety(text)["nonlinear_static"] is True


def test_missing_dynamic_kinetic_energy_is_not_hidden(decks):
    first, dynamic = decks["linear"].split("*END STEP", 1)
    with pytest.raises(ValueError, match="output missing"):
        cli.output_safety(first + "*END STEP" + dynamic.replace("ELSE,ELKE", "ELSE"))


def test_lin_nl_material_load_clamps_equivalent_except_routing_and_static_output(decks):
    def normalized(text):
        first, dynamic = text.split("*END STEP", 1)
        return (first.replace("ELSE,ELKE", "ELSE") + "*END STEP" + dynamic).replace(", NLGEOM=NO", "").replace(", NLGEOM", "")
    assert normalized(decks["linear"]) == normalized(decks["nonlinear"])


@pytest.mark.parametrize("key", ["hard_stop", "blocked", "complete"])
def test_cached_stops_or_complete_never_execute_again(science, tmp_path, monkeypatch, key):
    _forbid_science(monkeypatch)
    s = _summary()
    if key == "hard_stop":
        s["hard_stop"] = True
    elif key == "blocked":
        s["overall"] = "BLOCKED_BY_SOLVER"
    else:
        s["overall"] = "PILOT_COMPLETE_WITH_QUALIFICATIONS"
    _cached(monkeypatch, tmp_path, science, s)
    assert cli.main(["--run-pilot"]) is s


def test_linear_failure_stops_before_nonlinear_and_one_d(science, tmp_path, monkeypatch):
    _forbid_science(monkeypatch)
    s = _summary(); _cached(monkeypatch, tmp_path, science, s); calls = []
    def fail(b, c, item, summary, kind):
        calls.append(kind); summary["cases"][kind] = {"status": "FAIL"}
        summary["overall"] = "BLOCKED_BY_SOLVER"; return False
    monkeypatch.setattr(cli.base, "run_case", fail)
    assert cli.main(["--run-pilot"]) is s
    assert calls == ["linear"] and "nonlinear" not in s["cases"]


def test_linear_actual_audit_failure_hard_stops_without_nonlinear(science, tmp_path, monkeypatch):
    _forbid_science(monkeypatch)
    s = _summary(); _cached(monkeypatch, tmp_path, science, s); calls = []
    def success(b, c, item, summary, kind):
        calls.append(kind); summary["cases"][kind] = {"status": "PASS"}; return True
    monkeypatch.setattr(cli.base, "run_case", success)
    monkeypatch.setattr(cli, "audit_case", lambda *a: (_ for _ in ()).throw(ValueError("actual transfer mismatch")))
    cli.main(["--run-pilot"])
    assert calls == ["linear"]
    assert s["hard_stop"] is True and s["cases"]["linear"]["status"] == "FAIL"
    assert s["cases"]["linear"]["audit_failure"] == "actual transfer mismatch"


def test_through_linear_preserves_medium_first_no_one_d(science, tmp_path, monkeypatch):
    _forbid_science(monkeypatch)
    s = _summary(); _cached(monkeypatch, tmp_path, science, s); calls = []
    def success(b, c, item, summary, kind):
        calls.append(kind); summary["cases"][kind] = {"status": "PASS"}; return True
    def audit(b, c, summary, kind):
        summary["cases"][kind]["continuation_audit"] = {"status": "PASS"}
    monkeypatch.setattr(cli.base, "run_case", success); monkeypatch.setattr(cli, "audit_case", audit)
    cli.main(["--run-pilot", "--through-case", "linear"])
    assert calls == ["linear"] and "nonlinear" not in s["cases"]


def test_one_d_and_comparison_only_after_both_actual_audits(science, tmp_path, monkeypatch):
    _forbid_science(monkeypatch)
    s = _summary(); _cached(monkeypatch, tmp_path, science, s); calls = []
    def run(b, c, item, summary, kind):
        calls.append("job:" + kind); summary["cases"][kind] = {"status": "PASS"}; return True
    def audit(b, c, summary, kind):
        calls.append("audit:" + kind); summary["cases"][kind]["continuation_audit"] = {"status": "PASS"}
    def finish(b, c, item, summary):
        assert all(summary["cases"][k]["continuation_audit"]["status"] == "PASS" for k in ("linear", "nonlinear"))
        calls.append("references"); summary["overall"] = "PILOT_COMPLETE_WITH_QUALIFICATIONS"
    monkeypatch.setattr(cli.base, "run_case", run); monkeypatch.setattr(cli, "audit_case", audit)
    monkeypatch.setattr(cli.base, "finish_references_and_comparison", finish)
    cli.main(["--run-pilot"])
    assert calls == ["job:linear", "audit:linear", "job:nonlinear", "audit:nonlinear", "references"]
    assert s["one_d_authorization"]["id"] == "explicit_user_FEM3AR_2026_10_09"


def test_missing_energy_remains_partial_and_execution_is_separate():
    s = _summary()
    for kind in ("linear", "nonlinear"):
        s["cases"][kind] = {"status": "PASS", "continuation_audit": {"status": "PASS"}, "energy_status": "PARTIAL"}
    s["overall"] = "PILOT_COMPLETE_WITH_QUALIFICATIONS"; s["comparison"] = {"saved": True}
    cli.update_statuses(s)
    assert s["resume_statuses"]["NLSP_FEM3AR_LINEAR_DYNAMIC"] == "PASS"
    assert s["resume_statuses"]["NLSP_FEM3AR_NONLINEAR_DYNAMIC"] == "PASS"
    assert s["resume_statuses"]["NLSP_FEM3AR_ENERGY_DIAGNOSTICS"] == "PARTIAL"


@pytest.mark.parametrize("mode", ["--report-only", "--plot-only"])
def test_report_plot_zero_science_from_cache(mode, tmp_path, monkeypatch):
    _forbid_science(monkeypatch)
    s = _summary(); s["overall"] = "PARTIAL"
    monkeypatch.setattr(cli, "validate_cache", lambda b: s)
    plots = []; monkeypatch.setattr(cli.base, "plot_bundle", lambda b: plots.append(b))
    assert cli.main([mode, str(tmp_path)]) is s
    assert plots == ([tmp_path] if mode == "--plot-only" else [])


def test_preflight_cache_does_not_repeat_preparation(science, tmp_path, monkeypatch):
    _forbid_science(monkeypatch)
    s = _summary(); _cached(monkeypatch, tmp_path, science, s)
    monkeypatch.setattr(cli, "prepare", lambda *a: pytest.fail("matching preflight re-prepares"))
    assert cli.main(["--preflight"]) is s


def test_interrupted_attempt_cannot_be_silently_restarted(config, tmp_path, monkeypatch):
    b = tmp_path / "changed_code_hash"; b.mkdir()
    cli.write_json(b / "provenance.json", {"authorization": config["authorization"], "continuation_config": config})
    monkeypatch.setattr(cli, "OUTPUT", tmp_path)
    with pytest.raises(RuntimeError, match="no hidden retry"):
        cli.existing_attempt(config)


def test_changed_code_hash_does_not_bypass_authorized_attempt(config, tmp_path, monkeypatch):
    b = tmp_path / "original_fingerprint"; b.mkdir()
    item = {"authorization": config["authorization"], "continuation_config": config, "continuation_code_sha256": "old-hash"}
    cli.write_json(b / "provenance.json", item); cli.write_json(b / "manifest.json", {})
    monkeypatch.setattr(cli, "OUTPUT", tmp_path)
    s = {"overall": "BLOCKED_BY_SOLVER"}; monkeypatch.setattr(cli, "validate_cache", lambda path: s)
    found = cli.existing_attempt(config)
    assert found[0] == b and found[1]["continuation_code_sha256"] == "old-hash" and found[2] is s


def test_actual_times_must_match_without_interpolation():
    np.testing.assert_array_equal(cli.base._matches(np.array([0., .01, .025]), np.array([.01, .025])), [1, 2])
    with pytest.raises(ValueError, match="Actual-time"):
        cli.base._matches(np.array([0., .01, .025]), np.array([.02]))


def test_small_nonlinear_signal_uses_common_scale_without_alignment():
    x = np.linspace(0., 1., 5)
    first = np.array([[0., -1e-6, -2e-6, -1e-6, 0.]])
    second = first * 1.1
    result = cli.base.field_difference(first, second, x, scale=2.2e-6)
    assert result["absolute_max"] == pytest.approx(2e-7)
    assert result["relative_max"] == pytest.approx(2e-7 / 2.2e-6)
    assert result["phase_amplitude_alignment"] is False
    assert result["sampled_maxima_only"] is True



def test_actual_successful_linear_sta_print_precision_read_only(science):
    """Seven significant STA digits must not be mistaken for solver failure."""
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    job = cli.read_json(b / "cases/linear/job.json")
    assert job["returncode"] == 0 and job["failure"] is None
    result = cli.base.io.read_transient_sta(b / "cases/linear/motion.sta")
    static = [r for r in result["accepted_increments"] if r["step"] == 1]
    dynamic = [r for r in result["accepted_increments"] if r["step"] == 2]
    assert len(static) == 1 and len(dynamic) == 102
    assert static[-1]["step_time"] == 1.
    assert dynamic[0]["step_time"] > 0
    assert dynamic[-1]["step_time"] == pytest.approx(cli.base.dynamic_settings(science)["duration"], abs=5e-7)
    assert result["reported_cutbacks"] == 0
    assert all(r["total_time"] > 1. for r in dynamic)


def test_successful_solver_outputs_retained_when_recovery_pending():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    job = cli.read_json(b / "cases/linear/job.json")
    assert job["returncode"] == 0
    assert all((b / "cases/linear" / ("motion." + suffix)).stat().st_size > 0
               for suffix in ("inp", "dat", "frd", "sta", "stdout.txt"))
    assert not (PARENT / "cases/linear/motion.dat").stat().st_size



def test_actual_linear_preload_reproduces_saved_full_state():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    r = cli.read_json(b / "cases/linear/recovery.json")
    transfer = r["preload_transfer"]
    assert transfer["status"] == "PASS"
    for name in ("node_displacement_max_difference", "support_RF_max_difference",
                 "section_profile_max_difference", "stress_difference", "strain_difference"):
        assert transfer[name] == 0.
    assert transfer["static_end_total_time"] == 1.
    assert transfer["actual_metadata"]["dynamic_time"] is None


def test_actual_linear_motion_is_restoring_after_zero_work_release():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    r = cli.read_json(b / "cases/linear/recovery.json")
    assert r["dynamic_output_frames"] == r["dynamic_increments"] == 102
    assert 0 < r["dynamic_time_start"] < r["dynamic_time_end"]
    assert r["final_midspan_w"] < r["first_midspan_w"] < r["initial_midspan_w"]
    assert r["first_midspan_w_velocity"] < 0
    assert r["maximum_native_external_work_after_release"] == 0
    assert r["maximum_native_damping_work_after_release"] == 0
    assert r["cutbacks"] == 0
    audit = cli.read_json(b / "cases/linear/continuation_audit.json")
    assert audit["preload_force_imbalance_relative"] <= audit["equilibrium_gate_unchanged"] == 1e-5
    assert audit["preload_moment_imbalance_relative"] <= 1e-5


def test_actual_native_energy_offset_is_retained_as_qualified_partial():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    e = cli.read_json(b / "energy_transfer_diagnostic.json")
    assert e["status"] == "LOCALIZED_NATIVE_INITIALIZATION_ENERGY_BOOKKEEPING_OFFSET"
    assert e["actual_dynamic_initial_energy_stdout"] == 2 * e["actual_static_internal_energy_DAT"]
    assert e["max_abs_mechanical_change_relative_to_preload_ELSE"] > 1
    assert e["separate_max_abs_native_step_reference_drift"] < 1e-6
    assert e["energy_status_recommendation"] == "PARTIAL"
    assert e["native_data_not_renormalized_or_corrected"] is True
    assert e["new_solver_calls"] == 0
    assert e["external_work_after_release"] == e["damping_work_after_release"] == 0



def test_actual_pilot_attempt_accounting_and_all_gates():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    s = cli.read_json(b / "summary.json")
    assert s["overall"] == "PILOT_COMPLETE_WITH_QUALIFICATIONS"
    assert s["job_calls"]["CCX_production"] == 2
    assert s["job_calls"]["1D_nonlinear_ODE"] == 1
    for key in ("CCX_fixture", "Gmsh", "1D_static", "physical_root_search", "symbolic_derivations"):
        assert s["job_calls"][key] == 0
    assert [r["kind"] for r in s["attempts"]] == ["linear", "nonlinear"]
    assert all(r["status"] == "PASS" for r in s["attempts"])
    assert s["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED" and s["admitted"] is False
    assert s["strict_float64_qualification"] == "PARTIAL"
    for name, result in s["resume_statuses"].items():
        assert result == ("PARTIAL" if name.endswith("ENERGY_DIAGNOSTICS") else "PASS")


def test_actual_nonlinear_preload_release_and_energy_qualification():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    s = cli.read_json(b / "summary.json"); r = s["cases"]["nonlinear"]
    assert r["status"] == "PASS" and r["energy_status"] == "PARTIAL"
    assert r["static_increments"] == 10 and r["dynamic_increments"] == 102
    assert r["preload_transfer"]["node_displacement_max_difference"] == 0
    assert r["preload_transfer"]["support_RF_max_difference"] == 0
    assert r["preload_transfer"]["stress_difference"] == r["preload_transfer"]["strain_difference"] == 0
    assert r["maximum_native_external_work_after_release"] == r["maximum_native_damping_work_after_release"] == 0
    assert r["first_midspan_w_velocity"] < 0
    assert r["final_midspan_w"] < r["initial_midspan_w"]
    assert r["native_energy_static_to_dynamic_reference_jump_relative"] == 1.
    assert r["native_energy_history_not_offset_corrected"] is True


def test_actual_reference_initial_states_unfiltered_and_zero_velocities(science):
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    with np.load(ROOT / science["source_static"]["bundle"] / "one_d_p64.npz") as source:
        for kind in ("linear", "nonlinear"):
            with np.load(b / ("one_d_" + kind + ".npz")) as z:
                assert z["q"].shape[1] == z["velocity"].shape[1] == 252
                np.testing.assert_allclose(z["q"][0], source["q_" + kind], atol=5e-18, rtol=0)
                np.testing.assert_array_equal(z["velocity"][0], 0 * source["q_" + kind])
                assert z["times"][0] == 0 and z["times"][-1] == cli.base.dynamic_settings(science)["duration"]
    nl = cli.read_json(b / "one_d_nonlinear.json")["execution"]
    assert nl["rtol"] == 1e-10 and len(nl["atol"]) == 504
    assert nl["initial_coordinates_reused_exactly"] is True
    assert nl["no_dynamic_derivative_constraints"] is True
    assert nl["new_ODE_integrations"] == 1
    linear = cli.read_json(b / "one_d_linear.json")["execution"]
    assert linear["exact_in_time"] is True and linear["ODE_integrations"] == 0
    assert linear["linear_eigendecompositions"] == 1


def test_actual_correction_uses_all_native_shared_times_no_interpolation():
    b = ROOT / "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893"
    with np.load(b / "cases/linear/section_history.npz") as linear, np.load(b / "cases/nonlinear/section_history.npz") as nl, np.load(b / "dynamic_comparison.npz") as diff:
        np.testing.assert_array_equal(linear["time"], nl["time"])
        np.testing.assert_array_equal(diff["time"], linear["time"])
        assert len(diff["time"]) == 102
        np.testing.assert_array_equal(diff["three_d_correction"], nl["fields"] - linear["fields"])
    c = cli.read_json(b / "dynamic_comparison.json")
    assert c["dynamic_nonlinear_signal_certified"] is False
    assert c["c_eff_not_identical_to_MH_DOF"] is True
    assert c["initial_fields_not_amplitude_aligned"] is True
    for metric in c["nonlinear_corrections"].values():
        assert metric["phase_amplitude_alignment"] is False
        assert metric["sampled_maxima_only"] is True



def test_reference_authorization_annotation_preserves_legacy_execution(config, tmp_path):
    original = {"execution": {"authorization": {"user_authorized_FEM3A": True}, "new_ODE_integrations": 1},
                "diagnostics": {"max_relative_energy_drift": 1e-12}, "source": {"p": 64}}
    path = tmp_path / "one_d_nonlinear.json"
    cli.write_json(path, original)
    cli.record_reference_authorization(tmp_path, config["authorization"])
    annotated = cli.read_json(path)
    assert annotated["execution"] == original["execution"]
    assert annotated["diagnostics"] == original["diagnostics"]
    assert annotated["source"] == original["source"]
    assert annotated["continuation_authorization"]["id"] == "explicit_user_FEM3AR_2026_10_09"
    assert "not treated as a second authorization" in annotated["authorization_qualification"]
    assert not (tmp_path / "one_d_linear.json").exists()



def test_plot_validates_before_renderer_mutation_and_saves_after(tmp_path, monkeypatch):
    """Rendering may update a tracked figure; validate its old hash only first."""
    _forbid_science(monkeypatch)
    artifact = tmp_path / "figure.txt"; artifact.write_text("original", encoding="utf8")
    np.savez(tmp_path / "dynamic_comparison.npz", present=np.array(True))
    summary = {"preflight": {"T1": 1.}, "cases": {}}
    for kind in ("linear", "nonlinear"):
        np.savez(tmp_path / ("one_d_" + kind + ".npz"), times=[0., .05], energy_relative_drift=[0., 0.])
        case = tmp_path / "cases" / kind; case.mkdir(parents=True)
        cli.write_json(case / "energy.json", {"records": [{"step": 2, "dynamic_time": .05, "mechanical_energy": 2.}]})
        summary["cases"][kind] = {"native_initial_internal_energy": 1., "native_dynamic_bookkeeping_initial_energy": 2.}
    cli.write_json(tmp_path / "provenance.json", {"frozen": True})
    calls = []
    def validate(path):
        assert artifact.read_text(encoding="utf8") == "original", "validation incorrectly repeated after renderer changed figure"
        calls.append("validate"); return summary
    def render(path):
        calls.append("render"); artifact.write_text("rendered", encoding="utf8")
        return {"figures": 3, "new_scientific_calls": 0}
    def save(path, identity, data):
        assert artifact.read_text(encoding="utf8") == "rendered"
        assert identity == {"frozen": True} and data is summary
        calls.append("save")
    from matplotlib.figure import Figure
    monkeypatch.setattr(cli, "validate_cache", validate)
    monkeypatch.setattr(cli.base, "plot_bundle", render)
    monkeypatch.setattr(cli, "save", save)
    monkeypatch.setattr(Figure, "savefig", lambda *a, **k: None)
    result = cli.plot_bundle(tmp_path)
    assert calls == ["validate", "render", "save"]
    assert result["new_scientific_calls"] == 0
