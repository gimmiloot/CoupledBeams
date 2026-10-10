"""Attempt-accounting fixtures monkeypatch the only native execution call."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts.lib import nlsp_spatial_native_program as native


@pytest.fixture
def fixture(tmp_path, monkeypatch):
    root = Path(__file__).resolve().parents[1]
    config = json.loads((root/"data/input/nlsp_spatial_nonlinear_3d_fem_verification.json").read_text())
    binary = tmp_path/"ccx.exe"; binary.write_bytes(b"synthetic-not-executable")
    config.update(ccx_exe=str(binary), ccx_sha256=native.sha(binary))
    source = tmp_path/"source"; source.mkdir(); include = source/"solid_mesh.inp"; include.write_text("immutable fixture")
    decision = tmp_path/"pre_fem_decision.json"; decision.write_text('{"Stage_A":"PASS","Stage_B":"PASS"}')
    authorization = {"id": config["authorization"], "explicit_user_authorization": True,
        "stage_A_pass": True, "stage_B_pass": True, "source_manifest_verified": True,
        "frozen_decision_path": str(decision), "frozen_decision_sha256": native.sha(decision),
        "source_mesh_sha256": native.sha(include), "mesh_level": "medium", "nonlinear": False, "ordinal": 1}
    load = native.protocol.load_contract(.001, .00125)
    calls = {"native": 0, "recovery": 0}
    def generate(path, *args, **kwargs):
        Path(path).write_text("corrected generated fixture input")
        return {"load": load}, {"status": "PASS"}
    def execute(command, cwd, timeout, memory, prefix, env):
        calls["native"] += 1
        assert timeout == 1800 and memory == 4*1024**3
        assert all(env[name] == "1" for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS"))
        assert (cwd/"native_attempt.json").exists()
        assert json.loads((cwd/"native_attempt.json").read_text())["status"] == "STARTED"
        assert (cwd/"execution_code/before_native_execution").exists()
        (cwd/"motion.stdout.txt").write_text("JOB FINISHED\n")
        (cwd/"motion.stderr.txt").write_text("")
        return SimpleNamespace(returncode=0), {"returncode": 0, "failure": None, "seconds": .01}
    def recover(*args, **kwargs):
        calls["recovery"] += 1
        return {"status": "PASS", "energy_status": "PARTIAL"}
    monkeypatch.setattr(native.protocol, "write_spatial_input", generate)
    monkeypatch.setattr(native.protocol.fem1, "run_job", execute)
    monkeypatch.setattr(native.protocol, "recover_saved_outputs", recover)
    return tmp_path/"case", source, config, load, authorization, calls


def run(fixture, **auth_changes):
    case, source, config, load, authorization, _ = fixture
    authorization = {**authorization, **auth_changes}
    return native.run_attempt(case, source, None, {}, config, load,
                              authorization=authorization, remainingBudget=14400)


def test_one_call_then_cached_compute_zero_calls(fixture):
    first = run(fixture)
    assert first["status"] == "PASS" and first["new_solver_calls"] == 1
    second = run(fixture)
    assert second["status"] == "PASS" and second["new_solver_calls"] == 0
    assert fixture[-1] == {"native": 1, "recovery": 1}
    assert first["native_seconds"] >= 0 and first["recovery_seconds"] >= 0


def test_completed_cache_read_allowed_after_budget_exhaustion(fixture):
    run(fixture)
    case, source, config, load, auth, _ = fixture
    result = native.run_attempt(case, source, None, {}, config, load,
                                authorization=auth, remainingBudget=0)
    assert result["new_solver_calls"] == 0 and result["status"] == "PASS"
    assert fixture[-1]["native"] == 1


@pytest.mark.parametrize("key", ["explicit_user_authorization", "stage_A_pass", "stage_B_pass", "source_manifest_verified"])
def test_unpassed_prerequisite_forbids_native_call(fixture, key):
    with pytest.raises(ValueError, match="gate"):
        run(fixture, **{key: False})
    assert fixture[-1]["native"] == 0


def test_source_tamper_forbids_execution(fixture):
    (fixture[1]/"solid_mesh.inp").write_text("changed")
    with pytest.raises(ValueError, match="INP"):
        run(fixture)
    assert fixture[-1]["native"] == 0


def test_decision_tamper_forbids_execution(fixture):
    Path(fixture[4]["frozen_decision_path"]).write_text("changed")
    with pytest.raises(ValueError, match="decision"):
        run(fixture)


def test_fine_requires_informative_completed_medium_pair(fixture):
    with pytest.raises(ValueError, match="Fine pair"):
        run(fixture, ordinal=3, mesh_level="fine")
    assert fixture[-1]["native"] == 0


def test_invalid_case_order_and_fifth_job_are_forbidden(fixture):
    with pytest.raises(ValueError, match="order"):
        run(fixture, nonlinear=True)
    with pytest.raises(ValueError, match="four"):
        run(fixture, ordinal=5)


def test_actual_native_failure_is_immutable_and_no_retry(fixture, monkeypatch):
    def fail(command, case, *args):
        fixture[-1]["native"] += 1
        (case/"motion.stdout.txt").write_text("*ERROR actual failure\n")
        return SimpleNamespace(returncode=-1), {"failure": "native crash"}
    monkeypatch.setattr(native.protocol.fem1, "run_job", fail)
    first = run(fixture)
    assert first["status"] == "FAIL" and first["hard_stop"]
    second = run(fixture)
    assert second["status"] == "FAIL" and second["new_solver_calls"] == 0
    assert fixture[-1]["native"] == 1 and fixture[-1]["recovery"] == 0


def test_finished_output_parser_failure_reparses_without_ccx(fixture, monkeypatch):
    def fail(*args, **kwargs):
        fixture[-1]["recovery"] += 1
        raise ValueError("output parser fixture failure")
    monkeypatch.setattr(native.protocol, "recover_saved_outputs", fail)
    first = run(fixture)
    assert first["status"] == "OUTPUT_RECOVERY_PENDING" and first["hard_stop"]
    second = run(fixture)
    assert second["new_solver_calls"] == 0 and fixture[-1]["recovery"] == 1
    monkeypatch.setattr(native.protocol, "recover_saved_outputs", lambda *a, **k: {"status": "PASS"})
    third = run(fixture, reparse_pending=True)
    assert third["status"] == "PASS" and third["new_solver_calls"] == 0
    assert fixture[-1]["native"] == 1
    assert (fixture[0]/"execution_code/recovery_2").exists()


def test_started_attempt_refused_even_after_helper_code_change(fixture):
    run(fixture)
    case = fixture[0]; record = json.loads((case/native.RECORD).read_text())
    record.update(status="STARTED", hard_stop=True)
    native._save(case, record)
    next_record = run(fixture)
    assert next_record["status"] == "STARTED" and next_record["new_solver_calls"] == 0
    assert fixture[-1]["native"] == 1


def test_unledgered_files_never_overwritten(fixture):
    fixture[0].mkdir(); (fixture[0]/"motion.dat").write_text("old output")
    with pytest.raises(RuntimeError, match="unledgered"):
        run(fixture)
    assert (fixture[0]/"motion.dat").read_text() == "old output"


def test_cached_native_artifact_tamper_not_hidden(fixture):
    run(fixture); (fixture[0]/"motion.stdout.txt").write_text("modified")
    with pytest.raises(ValueError, match="artifact"):
        run(fixture)
    assert fixture[-1]["native"] == 1


def test_frozen_load_changed_not_reused(fixture):
    run(fixture)
    fixture[3].update(native.protocol.load_contract(.002, .0025))
    with pytest.raises(ValueError, match="immutable scientific attempt"):
        run(fixture)


def test_scope_and_remaining_budget_not_relaxed(fixture):
    case, source, config, load, auth, _ = fixture
    changed = copy.deepcopy(config); changed["FEM_budget"]["memory_limit_bytes"] = 8*1024**3
    with pytest.raises(ValueError, match="scope"):
        native.run_attempt(case, source, None, {}, changed, load, authorization=auth, remainingBudget=14400)
    with pytest.raises(ValueError, match="budget"):
        native.run_attempt(case, source, None, {}, config, load, authorization=auth, remainingBudget=0)
