"""One explicitly gated native spatial attempt; no CLI or mechanical solver.

The enclosing seven-field study owns the global four-attempt ledger and budget.
This helper prevents physical replay of an individual started/failed attempt.
Only an explicit read-only reparse may recover already finished native outputs.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
import time
from pathlib import Path

from scripts.lib import nlsp_spatial_fem_protocol as protocol

ROOT = Path(__file__).resolve().parents[2]
RECORD = "native_attempt.json"
MANIFEST = "native_attempt_manifest.json"
AUTHORIZATION = "explicit_user_NLSP_spatial_seven_field_verification_2026_10_10"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _save(case, record):
    protocol.static.write_json(case/RECORD, record)
    names = (RECORD, "motion.inp", "motion.stdout.txt", "motion.stderr.txt",
        "motion.dat", "motion.frd", "motion.sta", "job.json", "science_config.json",
        "input_contract.json", "execution_environment.json", "authorization_gate.json",
        "section_history.npz", "initial_sections.npz", "increments.json", "recovery.json",
        "preload_equilibrium.json", "energy.json", "stdout_energy.json")
    protocol.static.write_json(case/MANIFEST,
        {"schema": "nlsp-spatial-native-attempt-v1", "artifact_hashes":
         {name: sha(case/name) for name in names if (case/name).is_file()}})


def _verify_cache(case):
    manifest = protocol.static.read_json(case/MANIFEST)
    for name, digest in manifest["artifact_hashes"].items():
        if sha(case/name) != digest:
            raise ValueError("Saved native attempt artifact changed: " + name)
    return protocol.static.read_json(case/RECORD)


def _gate(config, frozen_load, authorization, source_mesh, remaining_budget):
    if config["authorization"] != AUTHORIZATION or authorization.get("id") != AUTHORIZATION:
        raise ValueError("Independent explicit seven-field authorization required")
    for key in ("explicit_user_authorization", "stage_A_pass", "stage_B_pass", "source_manifest_verified"):
        if authorization.get(key) is not True:
            raise ValueError("Unpassed prerequisite gate: " + key)
    decision = Path(authorization["frozen_decision_path"])
    if sha(decision) != authorization["frozen_decision_sha256"]:
        raise ValueError("Frozen pre-FEM decision hash mismatch")
    if sha(source_mesh) != authorization["source_mesh_sha256"]:
        raise ValueError("Immutable source INP hash mismatch")
    if (config["geometry"] != {"L": 1., "b": .2, "h": .1}
        or config["material"] != {"E": 1., "rho": 1., "nu": .3, "kappa": 5/6}
        or config["omega1"] != .6054167303477958 or config["horizon_T1"] != .25
        or config["FEM_dynamic"] != protocol.DEFAULT_DYNAMIC
        or config["FEM_budget"] != protocol.RESOURCE_POLICY
        or config["execution_mode"] != "EXPLORATORY_NOT_CERTIFIED" or config["admitted"] is not False
        or any(config[name] for name in ("new_meshes", "new_modal_FEM_jobs", "full_period_FEM"))):
        raise ValueError("Frozen scientific/resource scope changed")
    ordinal = authorization["ordinal"]
    if ordinal not in (1, 2, 3, 4) or type(authorization["nonlinear"]) is not bool:
        raise ValueError("One of the four predeclared cases required")
    level = "medium" if ordinal < 3 else "fine"
    if authorization["mesh_level"] != level or authorization["nonlinear"] != (ordinal % 2 == 0):
        raise ValueError("Medium linear/NL then conditional fine linear/NL order required")
    if level == "fine" and (authorization.get("completed_medium_pair") is not True
                           or authorization.get("medium_signal_resolved") is not True):
        raise ValueError("Fine pair requires completed and informative medium pair")
    if not math.isfinite(remaining_budget) or not 0 <= remaining_budget <= 14400:
        raise ValueError("Finite nonnegative remaining CCX budget required")
    local = frozen_load["local_acceleration"]
    physical = protocol.load_contract(local[1], local[2], rho=config["material"]["rho"])
    for key in ("global_acceleration", "global_direction", "magnitude", "line_load_local", "total_force_global"):
        if frozen_load[key] != physical[key]:
            raise ValueError("Frozen load component/resultant contract differs: " + key)
    binary = Path(config["ccx_exe"])
    if sha(binary) != config["ccx_sha256"]:
        raise ValueError("Previously verified native binary changed or unavailable")
    return level, binary


def _snapshot(case, phase):
    directory = case/"execution_code"/phase
    directory.mkdir(parents=True, exist_ok=False)
    paths = (Path(__file__), Path(protocol.__file__), Path(protocol.pilot.__file__),
        Path(protocol.static.__file__), Path(protocol.fem1.__file__), Path(protocol.io.__file__))
    hashes = {}
    for source in paths:
        name = source.relative_to(ROOT).as_posix()
        destination = directory/name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        hashes[name] = sha(source)
    protocol.static.write_json(directory/"code_identity.json", hashes)
    return hashes


def _recover(case, mesh, audit, record, authorization):
    phase = "recovery_" + str(record.get("recovery_attempts", 0)+1)
    record["recovery_helper_sha256"] = _snapshot(case, phase)
    record["recovery_attempts"] = record.get("recovery_attempts", 0)+1
    _save(case, record)
    start = time.perf_counter()
    try:
        recovered = protocol.recover_saved_outputs(case, mesh, audit, record["science_config"],
            nonlinear=authorization["nonlinear"])
        record.update(status="PASS", recovery=recovered, hard_stop=False)
        record.pop("recovery_failure", None)
    except Exception as error:
        record.update(status="OUTPUT_RECOVERY_PENDING", recovery_failure=str(error), hard_stop=True)
    finally:
        record["recovery_seconds"] = record.get("recovery_seconds", 0.)+time.perf_counter()-start
        _save(case, record)
    return record


def run_attempt(caseDir, sourceMeshDir, mesh, audit, config, frozenLoad, *,
                authorization, remainingBudget):
    """Return status and actual new-call count; never restart a native attempt.

Authorization keys: id, explicit_user_authorization, stage_A_pass,
stage_B_pass, frozen_decision_path/sha256, source_manifest_verified,
source_mesh_sha256, mesh_level, nonlinear, ordinal. The conditional fine pair
also requires completed_medium_pair and medium_signal_resolved. A pending
finished output is reparsed only with reparse_pending=True; scientific calls=0.
"""
    case = Path(caseDir); source = Path(sourceMeshDir)
    level, binary = _gate(config, frozenLoad, authorization, source/"solid_mesh.inp", float(remainingBudget))
    identity = {"authorization_id": authorization["id"], "ordinal": authorization["ordinal"],
        "mesh_level": level, "nonlinear": authorization["nonlinear"],
        "frozen_decision_sha256": authorization["frozen_decision_sha256"],
        "source_mesh_sha256": authorization["source_mesh_sha256"],
        "frozen_load_sha256": fingerprint(frozenLoad), "config_sha256": fingerprint(config)}
    if case.exists() and (case/RECORD).exists():
        record = _verify_cache(case)
        if record["identity"] != identity:
            raise ValueError("Existing case uses a different immutable scientific attempt")
        if record["status"] == "OUTPUT_RECOVERY_PENDING" and authorization.get("reparse_pending") is True:
            record = _recover(case, mesh, audit, record, authorization)
        return {**record, "new_solver_calls": 0, "cache_or_reparse": True}
    if case.exists() and any(case.iterdir()):
        raise RuntimeError("Existing unledgered native files: no automatic overwrite or retry")
    if remainingBudget <= 0:
        raise ValueError("Positive remaining CCX budget required for a new attempt")
    case.mkdir(parents=True, exist_ok=True)
    static_settings = config.get("FEM_static_settings", protocol.static.STATIC_CONTROL_DEFAULTS)
    science, input_gate = protocol.write_spatial_input(case/"motion.inp", source, mesh, audit,
        material=config["material"], g_n=frozenLoad["local_acceleration"][1],
        g_k=frozenLoad["local_acceleration"][2], omega1=config["omega1"],
        static_settings=static_settings, nonlinear=authorization["nonlinear"],
        dynamic=config["FEM_dynamic"], horizon_T1=config["horizon_T1"])
    protocol.static.write_json(case/"science_config.json", science)
    protocol.static.write_json(case/"input_contract.json", input_gate)
    protocol.static.write_json(case/"authorization_gate.json", authorization)
    helpers = _snapshot(case, "before_native_execution")
    env = dict(os.environ)
    thread_names = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS")
    env.update({name: "1" for name in thread_names})
    timeout = min(config["FEM_budget"]["job_timeouts_seconds"][level], float(remainingBudget))
    command = [str(binary), "motion"]
    protocol.static.write_json(case/"execution_environment.json", {"command": command, "cwd": str(case.resolve()),
        "thread_environment": {name: env[name] for name in thread_names}, "PATH_unchanged": True,
        "timeout_seconds": timeout, "memory_limit_bytes": 4*1024**3,
        "binary_sha256": sha(binary), "runtime_DLLs": {path.name: sha(path) for path in binary.parent.glob("*.dll")}})
    record = {"identity": identity, "status": "STARTED", "solver_calls": 1, "hard_stop": True,
        "automatic_retry": False, "science_config": science, "input_sha256": sha(case/"motion.inp"),
        "actual_helper_sha256": helpers, "native_seconds": 0., "recovery_seconds": 0.}
    _save(case, record)
    start = time.perf_counter()
    try:
        result, stats = protocol.fem1.run_job(command, case, timeout, 4*1024**3, case/"motion", env)
        protocol.static.write_json(case/"job.json", stats)
        record["job"] = stats
        text = (case/"motion.stdout.txt").read_text(encoding="utf8", errors="strict")
        if result.returncode != 0 or stats.get("failure") or "JOB FINISHED" not in text.upper() or "*ERROR" in text.upper():
            raise RuntimeError("Actual native failure: " + str(stats))
        record.update(status="OUTPUT_RECOVERY_PENDING", hard_stop=True)
    except Exception as error:
        record.update(status="FAIL", native_failure=str(error), hard_stop=True)
    finally:
        record["native_seconds"] = time.perf_counter()-start
        _save(case, record)
    if record["status"] == "OUTPUT_RECOVERY_PENDING":
        record = _recover(case, mesh, audit, record, authorization)
    return {**record, "new_solver_calls": 1, "cache_or_reparse": False}
