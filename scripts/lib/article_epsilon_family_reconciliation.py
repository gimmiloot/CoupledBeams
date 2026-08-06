"""Zero-solve promotion of verified family-repair shadow results.

The reconciliation layer reads immutable point-cache/shadow artifacts and
writes a versioned article-facing overlay.  It intentionally imports no matrix
evaluator, root solver, continuation helper, or local-repair implementation.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
from typing import Mapping, Sequence


RECONCILIATION_VERSION = "article_epsilon_family_reconciliation_v1"
PROMOTION_POLICY = "verified-only"
SCIENTIFIC_SCOPE = "isotropic_circular_coupled_rods_eb_timoshenko"
OUTPUT_DIRECTORY_NAME = "family_local_repair_reconciliation"
SHADOW_DIRECTORY_NAME = "family_local_repair_shadow"
K_MAX = 10
DELTA_TOL = 0.10

REQUIRED_SHADOW_GATES = (
    "scope_isolation_gate",
    "source_cache_preservation_gate",
    "resolved_case_preservation_gate",
    "shadow_repair_gate",
    "strict_avoidance_gate",
    "resume_determinism_gate",
    "production_resume_readiness_gate",
)

RECONCILIATION_MANIFEST_FIELDS = (
    "case_id", "case_identity", "grid_group", "epsilon_0", "beta", "mu", "eta",
    "shadow_present", "point_cache_present", "source_class", "source_cache_hash",
    "scientific_scope",
)
ELIGIBILITY_FIELDS = (
    "case_id", "candidate_type", "eligible", "scientific_scope_pass",
    "shadow_gates_pass", "cache_fingerprint_pass", "manifest_sha_pass",
    "shadow_status_pass", "n_true_pass", "required_guard_pass",
    "cluster_guard_pass", "matrix_confirmation_pass", "inventory_stability_pass",
    "force_strict_pass", "original_result_conflict_pass", "provenance_complete",
    "eligibility_reason", "source_cache_hash", "repair_cache_path",
)
PROMOTION_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta", "execution_status", "N_true",
    "first_failed_mode", "required_guard", "EB_roots_json", "Timoshenko_roots_json",
    "delta_f_json", *(f"delta_f_{index}" for index in range(1, K_MAX + 1)),
    "detector_provenance_json", "local_repair_provenance_json", "inferred_window_json",
    "accepted_roots_json", "source_cache_hash", "repair_cache_hash",
    "detector_version", "repair_version", "general_evaluator_version",
    "EB_evaluator_version", "Timoshenko_evaluator_version", "tolerance_hash",
    "result_origin", "promotion_status", "promotion_version", "promotion_timestamp_utc",
    "scientific_scope",
)
RECONCILED_FIELDS = (
    "case_id", "case_identity", "grid_group", "claim_eligible", "epsilon_0", "beta",
    "mu", "eta", "s_max", "thin_0p1_flag", "original_execution_status",
    "reconciled_execution_status", "original_N_true", "reconciled_N_true",
    "first_failed_mode", "required_guard", "result_origin", "promotion_status",
    "reconciliation_status", "detector_status", "local_repair_status",
    "required_prefix_guard_status", "upper_spectrum_audit_status",
    "full_K10_control_status", "unresolved_reason", "point_cache_present",
    "source_cache_hash", "EB_roots_json", "Timoshenko_roots_json", "delta_f_json",
    *(f"delta_f_{index}" for index in range(1, K_MAX + 1)),
    "promotion_provenance_json", "scientific_scope",
)
COMPARISON_FIELDS = (
    "case_id", "original_execution_status", "reconciled_execution_status",
    "original_N_true", "reconciled_N_true", "N_true_changed",
    "first_failed_mode_changed", "root_inventory_changed", "reconciliation_status",
    "result_origin", "source_cache_hash", "scientific_scope",
)
DEFERRED_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta", "execution_status", "N_true",
    "first_apparent_failed_mode", "required_guard", "unresolved_reason",
    "inferred_window_json", "local_repair_status", "force_strict_requested",
    "force_strict_executed", "source_cache_hash", "scientific_scope",
)
UNRESOLVED_FIELDS = (
    "case_id", "original_execution_status", "reconciled_execution_status",
    "first_apparent_failed_mode", "required_guard", "detector_status",
    "local_repair_status", "inferred_window_json", "candidate_outcome_json",
    "unresolved_reason", "future_audit_category", "source_cache_hash",
    "scientific_scope",
)
RESUME_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta", "resume_status", "skip_reason",
    "original_execution_status", "reconciled_execution_status", "point_cache_present",
    "promotion_status", "scientific_scope",
)
GATE_FIELDS = ("gate", "status", "metric", "value", "explanation")
OPERATION_FIELDS = ("operation", "count", "scope")

DETERMINISTIC_DATA_CSVS = (
    "reconciliation_manifest.csv", "eligibility_audit.csv", "promotion_overlay.csv",
    "original_vs_reconciled_cases.csv", "reconciled_case_results.csv",
    "promoted_cases.csv", "deferred_cases.csv", "unresolved_provenance.csv",
    "resume_plan.csv", "operation_counts.csv",
)
ALL_CSVS = (*DETERMINISTIC_DATA_CSVS, "gate_summary.csv")


class ReconciliationIntegrityError(RuntimeError):
    """Raised before promotion when an immutable input contract is violated."""


def validate_scope(scientific_scope: str) -> str:
    if scientific_scope != SCIENTIFIC_SCOPE:
        raise ValueError(
            f"reconciliation supports only {SCIENTIFIC_SCOPE!r}; received {scientific_scope!r}"
        )
    return scientific_scope


def shadow_input_contract_valid(
    *,
    required_gates_pass: bool,
    actual_cache_fingerprint: str,
    expected_cache_fingerprint: str,
    actual_manifest_sha256: str,
    expected_manifest_sha256: str,
    artifact_mismatches: Sequence[str] = (),
) -> bool:
    """Pure integrity predicate used before any promotion is considered."""

    return (
        required_gates_pass
        and actual_cache_fingerprint == expected_cache_fingerprint
        and actual_manifest_sha256 == expected_manifest_sha256
        and not artifact_mismatches
    )


def confirmed_n_true_value(value: object) -> bool:
    """Distinguish the valid scientific value zero from an absent result."""

    if value is None:
        return False
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return False
    parsed = float(text)
    return math.isfinite(parsed) and parsed.is_integer() and 0 <= parsed <= K_MAX


def forbid_scientific_operation(operation: str) -> None:
    raise ReconciliationIntegrityError(
        f"reconcile-only workflow forbids scientific operation: {operation}"
    )


def pending_family_status(
    *,
    primary_complete: bool,
    family_context_complete: bool,
    poststage_outcome: str = "",
) -> str:
    """Return a parent-owned transient/final status independent of worker order."""

    if not primary_complete:
        return "not_attempted"
    if not family_context_complete and not poststage_outcome:
        return "pending_family_inventory_check"
    outcomes = {
        "resolved": "resolved_after_local_repair",
        "deferred": "deferred_expensive_strict",
        "context_incomplete": "deferred_family_context_incomplete",
    }
    if poststage_outcome not in outcomes:
        raise ValueError("completed family post-stage requires a declared outcome")
    return outcomes[poststage_outcome]


def deterministic_family_poststage(
    rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Sort primary completions by family/beta before assigning final statuses."""

    ordered = sorted(
        (dict(row) for row in rows),
        key=lambda row: (
            float(row["epsilon_0"]), float(row["mu"]), float(row["eta"]),
            float(row["beta"]), str(row["case_id"]),
        ),
    )
    for row in ordered:
        row["execution_status"] = pending_family_status(
            primary_complete=bool(row.get("primary_complete", True)),
            family_context_complete=bool(row.get("family_context_complete", False)),
            poststage_outcome=str(row.get("poststage_outcome", "context_incomplete")),
        )
    return ordered


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})
    os.replace(temporary, path)


def _atomic_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _directory_fingerprint(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    count = 0
    size = 0
    for item in sorted((value for value in path.rglob("*") if value.is_file()), key=lambda value: value.relative_to(path).as_posix().lower()):
        relative = item.relative_to(path).as_posix()
        file_hash = _sha256(item)
        stat = item.stat()
        count += 1
        size += stat.st_size
        digest.update(relative.encode("utf-8"))
        digest.update(str(stat.st_size).encode("ascii"))
        digest.update(file_hash.encode("ascii"))
    return {"file_count": count, "size_bytes": size, "sha256": digest.hexdigest()}


def _point_cache_fingerprint(
    cache_dir: Path,
    inventory_rows: Sequence[Mapping[str, str]],
) -> tuple[str, list[str]]:
    digest = hashlib.sha256()
    mismatches: list[str] = []
    expected_paths = {str(row["relative_path"]) for row in inventory_rows}
    actual_paths = {
        item.relative_to(cache_dir).as_posix()
        for item in cache_dir.rglob("*") if item.is_file()
    }
    if expected_paths != actual_paths:
        mismatches.append("cache_path_inventory_mismatch")
    for row in sorted(inventory_rows, key=lambda item: str(item["relative_path"]).lower()):
        relative = str(row["relative_path"])
        path = (cache_dir / relative).resolve()
        if not path.is_relative_to(cache_dir.resolve()) or not path.is_file():
            mismatches.append(relative + ":missing_or_outside")
            continue
        actual_hash = _sha256(path)
        actual_size = path.stat().st_size
        if actual_hash != str(row["sha256"]):
            mismatches.append(relative + ":sha256")
        if actual_size != int(row["size_bytes"]):
            mismatches.append(relative + ":size")
        digest.update(relative.encode("utf-8"))
        digest.update(str(actual_size).encode("ascii"))
        digest.update(actual_hash.encode("ascii"))
    return digest.hexdigest(), mismatches


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _format_float(value: float) -> str:
    return format(float(value), ".17g")


def _case_sort(row: Mapping[str, object]) -> tuple[float, float, float, float, str]:
    return (
        float(row.get("epsilon_0", 0.0)), float(row.get("beta", row.get("beta_deg", 0.0))),
        float(row.get("mu", 0.0)), float(row.get("eta", 0.0)), str(row.get("case_id", "")),
    )


def _spectra_map(
    rows: Sequence[Mapping[str, str]],
) -> dict[tuple[str, str], dict[str, list[float]]]:
    result: dict[tuple[str, str], dict[str, list[float]]] = defaultdict(
        lambda: {"before": [], "shadow": []}
    )
    for row in sorted(rows, key=lambda item: (item["case_id"], item["theory"], int(item["sorted_rank"]))):
        key = (row["case_id"], row["theory"])
        if row.get("lambda_before", ""):
            result[key]["before"].append(float(row["lambda_before"]))
        if row.get("lambda_shadow", ""):
            result[key]["shadow"].append(float(row["lambda_shadow"]))
    return result


def _delta_payload(
    eb_roots: Sequence[float],
    timo_roots: Sequence[float],
) -> tuple[list[float], int | None, int | None]:
    deltas = [
        abs(float(eb_roots[index]) ** 2 - float(timo_roots[index]) ** 2)
        / float(timo_roots[index]) ** 2
        for index in range(min(K_MAX, len(eb_roots), len(timo_roots)))
    ]
    first_failure = next(
        (index for index, value in enumerate(deltas, start=1) if value > DELTA_TOL),
        None,
    )
    if first_failure is None and len(deltas) < K_MAX:
        return deltas, None, None
    n_true = K_MAX if first_failure is None else first_failure - 1
    required_guard = 11 if n_true == K_MAX else first_failure + 1
    return deltas, n_true, required_guard


def _stable_promotion_timestamp(shadow_dir: Path) -> str:
    candidates: list[str] = []
    for path in sorted((shadow_dir / "logs").glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        gates = payload.get("gates", {})
        if isinstance(gates, Mapping) and gates.get("production_resume_readiness_gate") == "PASS":
            candidates.append(str(payload.get("timestamp", "")))
    if not candidates:
        raise ReconciliationIntegrityError("shadow PASS log timestamp is missing")
    return max(candidates)


def _shadow_artifact_hashes(shadow_dir: Path, metadata: Mapping[str, object]) -> list[str]:
    mismatches: list[str] = []
    expected = metadata.get("deterministic_csv_hashes", {})
    if not isinstance(expected, Mapping) or not expected:
        return ["missing_deterministic_csv_hashes"]
    for name, value in expected.items():
        path = shadow_dir / str(name)
        if not path.is_file() or _sha256(path) != str(value):
            mismatches.append(str(name))
    return mismatches


def _repair_cache_payload(shadow_dir: Path, case_id: str) -> tuple[Path | None, dict[str, object]]:
    matches = sorted((shadow_dir / "cache" / "coarse_cases").glob(f"{case_id}_*.json"))
    if len(matches) != 1:
        return None, {}
    return matches[0], json.loads(matches[0].read_text(encoding="utf-8"))


def _candidate_consistency(
    entries: Sequence[Mapping[str, object]],
    candidate_rows: Sequence[Mapping[str, str]],
) -> bool:
    accepted = [
        float(row["lambda_candidate"])
        for row in candidate_rows
        if _truthy(row.get("accepted", "")) and row.get("lambda_candidate", "")
    ]
    return bool(entries) and all(
        any(abs(float(entry["value"]) - value) <= 1.0e-8 * max(1.0, abs(value)) for value in accepted)
        for entry in entries
    )


def _promotion_candidate(
    row: Mapping[str, str],
    *,
    manifest_row: Mapping[str, str],
    shadow_dir: Path,
    shadow_integrity: Mapping[str, bool],
    spectra: Mapping[tuple[str, str], Mapping[str, list[float]]],
    detector_rows: Sequence[Mapping[str, str]],
    candidate_rows: Sequence[Mapping[str, str]],
    window_rows: Sequence[Mapping[str, str]],
    force_strict_executed: int,
    promotion_timestamp: str,
) -> tuple[dict[str, object], dict[str, object] | None]:
    case_id = str(row["case_id"])
    repair_path, repair_payload = _repair_cache_payload(shadow_dir, case_id)
    identity = repair_payload.get("identity", {}) if isinstance(repair_payload, Mapping) else {}
    result = repair_payload.get("result", {}) if isinstance(repair_payload, Mapping) else {}
    identity = identity if isinstance(identity, Mapping) else {}
    result = result if isinstance(result, Mapping) else {}
    entries = result.get("entries", ()) if isinstance(result.get("entries", ()), Sequence) else ()
    case_candidates = [item for item in candidate_rows if item["case_id"] == case_id]
    case_detectors = [item for item in detector_rows if item["case_id"] == case_id]
    case_windows = [item for item in window_rows if item["case_id"] == case_id]
    eb_roots = list(spectra.get((case_id, "Euler-Bernoulli"), {}).get("shadow", ()))
    timo_roots = list(spectra.get((case_id, "Timoshenko"), {}).get("shadow", ()))
    deltas, computed_n, computed_guard = _delta_payload(eb_roots, timo_roots)
    shadow_n_text = str(row.get("shadow_N_true", "")).strip()
    try:
        shadow_n = int(float(shadow_n_text))
        n_finite = math.isfinite(float(shadow_n_text)) and 0 <= shadow_n <= K_MAX
    except (TypeError, ValueError):
        shadow_n = -1
        n_finite = False
    original_n = str(row.get("original_N_true", "")).strip()
    scope_pass = row.get("scientific_scope") == SCIENTIFIC_SCOPE and identity.get("scientific_scope") == SCIENTIFIC_SCOPE
    status_pass = row.get("shadow_execution_status") == "shadow_resolved_after_local_repair"
    n_pass = n_finite and computed_n == shadow_n
    guard = int(float(row.get("required_guard", "0") or 0))
    guard_pass = (
        row.get("required_prefix_guard_status") == "resolved_after_local_repair"
        and computed_guard == guard
        and len(eb_roots) >= guard
        and len(timo_roots) >= guard
    )
    cluster_pass = row.get("family_inventory_status") == "required_defect_matrix_confirmed"
    matrix_pass = (
        result.get("status") == "resolved_after_local_repair"
        and _candidate_consistency(entries, case_candidates)
    )
    stable_pass = result.get("status") == "resolved_after_local_repair" and bool(result.get("stage_roots"))
    force_pass = force_strict_executed == 0 and identity.get("force_strict_enabled") is False
    conflict_pass = not original_n or original_n == shadow_n_text
    provenance_keys = (
        "source_hash", "detector_version", "repair_algorithm_version",
        "general_spectrum_evaluator_version", "eb_evaluator_version",
        "timo_evaluator_version", "tolerance_hash", "inferred_window",
    )
    provenance_pass = (
        repair_path is not None
        and all(identity.get(key) not in (None, "", {}) for key in provenance_keys)
        and identity.get("source_hash") == row.get("source_cache_hash")
        and bool(entries) and bool(case_windows) and bool(case_detectors)
    )
    checks = {
        "scientific_scope_pass": scope_pass,
        "shadow_gates_pass": bool(shadow_integrity["gates"]),
        "cache_fingerprint_pass": bool(shadow_integrity["cache"]),
        "manifest_sha_pass": bool(shadow_integrity["manifest"]),
        "shadow_status_pass": status_pass,
        "n_true_pass": n_pass,
        "required_guard_pass": guard_pass,
        "cluster_guard_pass": cluster_pass,
        "matrix_confirmation_pass": matrix_pass,
        "inventory_stability_pass": stable_pass,
        "force_strict_pass": force_pass,
        "original_result_conflict_pass": conflict_pass,
        "provenance_complete": provenance_pass,
    }
    eligible = all(checks.values())
    failed = [key for key, passed in checks.items() if not passed]
    audit = {
        "case_id": case_id, "candidate_type": "shadow_local_repair_resolution",
        "eligible": eligible, **checks,
        "eligibility_reason": "verified_shadow_promotable" if eligible else ";".join(failed),
        "source_cache_hash": row.get("source_cache_hash", ""),
        "repair_cache_path": repair_path.relative_to(shadow_dir).as_posix() if repair_path else "",
    }
    if not eligible:
        return audit, None
    first_failure = "" if computed_n == K_MAX else int(computed_n) + 1
    delta_fields = {
        f"delta_f_{index}": _format_float(deltas[index - 1]) if index <= len(deltas) else ""
        for index in range(1, K_MAX + 1)
    }
    detector_provenance = [dict(item) for item in case_detectors]
    repair_provenance = {
        "repair_stage": result.get("repair_stage"),
        "status": result.get("status"),
        "matrix_evaluations_saved_reference": result.get("matrix_evaluations"),
        "stage_matrix_evaluations": result.get("stage_matrix_evaluations"),
        "stage_roots": result.get("stage_roots"),
        "beta_probes": result.get("beta_probes"),
        "block_classification": result.get("block_classification"),
    }
    promotion = {
        "case_id": case_id, "epsilon_0": manifest_row["epsilon_0"],
        "beta": manifest_row["beta_deg"], "mu": manifest_row["mu"], "eta": manifest_row["eta"],
        "execution_status": "resolved_after_local_repair", "N_true": computed_n,
        "first_failed_mode": first_failure, "required_guard": computed_guard,
        "EB_roots_json": _json(eb_roots), "Timoshenko_roots_json": _json(timo_roots),
        "delta_f_json": _json(deltas), **delta_fields,
        "detector_provenance_json": _json(detector_provenance),
        "local_repair_provenance_json": _json(repair_provenance),
        "inferred_window_json": _json(identity["inferred_window"]),
        "accepted_roots_json": _json(entries),
        "source_cache_hash": row["source_cache_hash"],
        "repair_cache_hash": _sha256(repair_path),
        "detector_version": identity["detector_version"],
        "repair_version": identity["repair_algorithm_version"],
        "general_evaluator_version": identity["general_spectrum_evaluator_version"],
        "EB_evaluator_version": identity["eb_evaluator_version"],
        "Timoshenko_evaluator_version": identity["timo_evaluator_version"],
        "tolerance_hash": identity["tolerance_hash"],
        "result_origin": "family_local_repair_promotion",
        "promotion_status": "verified_shadow_promoted",
        "promotion_version": RECONCILIATION_VERSION,
        "promotion_timestamp_utc": promotion_timestamp,
        "scientific_scope": SCIENTIFIC_SCOPE,
    }
    return audit, promotion


def _deferred_audit_row(row: Mapping[str, str]) -> dict[str, object]:
    checks = {
        "scientific_scope_pass": row.get("scientific_scope") == SCIENTIFIC_SCOPE,
        "shadow_gates_pass": True, "cache_fingerprint_pass": True,
        "manifest_sha_pass": True, "shadow_status_pass": False,
        "n_true_pass": row.get("shadow_N_true", "") in ("", "NaN"),
        "required_guard_pass": row.get("required_prefix_guard_status") == "unresolved_without_expensive_strict",
        "cluster_guard_pass": False, "matrix_confirmation_pass": False,
        "inventory_stability_pass": False, "force_strict_pass": True,
        "original_result_conflict_pass": row.get("original_N_true", "") == "",
        "provenance_complete": bool(row.get("source_cache_hash")),
    }
    return {
        "case_id": row["case_id"], "candidate_type": "shadow_deferred",
        "eligible": False, **checks, "eligibility_reason": "shadow_not_resolved",
        "source_cache_hash": row.get("source_cache_hash", ""), "repair_cache_path": "",
    }


def _reference_validation(coarse_dir: Path, shadow_metadata: Mapping[str, object]) -> dict[str, object]:
    readiness_dir = coarse_dir.parent / "solver_readiness_v2"
    validation = _read_csv(readiness_dir / "validation_24_cases.csv")
    readiness_pass = len(validation) == 24 and all(
        row["solver_readiness_case_status"] == "READY"
        and int(row["full_N_true"]) == int(row["auto_N_true"])
        for row in validation
    )
    former = _read_csv(readiness_dir / "seven_case_manifest.csv")
    validation_by_id = {row["validation_id"]: row for row in validation}
    former_pass = len(former) == 7 and all(
        row["validation_id"] in validation_by_id
        and int(validation_by_id[row["validation_id"]]["full_N_true"])
        == int(validation_by_id[row["validation_id"]]["auto_N_true"])
        for row in former
    )
    s3 = shadow_metadata.get("S3_regressions", {})
    s3_pass = (
        isinstance(s3, Mapping)
        and int(s3.get("S3_12_N_true", -1)) == 4
        and int(s3.get("S3_14_N_true", -1)) == 4
        and abs(float(s3.get("S3_12_delta_f_5", math.nan)) - 0.11739469908796035) <= 1.0e-15
        and abs(float(s3.get("S3_14_delta_f_5", math.nan)) - 0.10050934855181458) <= 1.0e-15
    )
    return {
        "readiness_pass": readiness_pass, "readiness_count": len(validation),
        "former_pass": former_pass, "former_count": len(former),
        "S3_pass": s3_pass, "S3": dict(s3) if isinstance(s3, Mapping) else {},
    }


def _csv_hashes(output_dir: Path, names: Sequence[str]) -> dict[str, str]:
    return {name: _sha256(output_dir / name) for name in names}


def load_ready_resume_ids(coarse_dir: Path) -> set[str]:
    output_dir = coarse_dir / OUTPUT_DIRECTORY_NAME
    gates = {row["gate"]: row["status"] for row in _read_csv(output_dir / "gate_summary.csv")}
    if gates.get("production_resume_readiness_gate") != "PASS":
        raise ReconciliationIntegrityError("production resume requires a PASS reconciliation gate")
    return {
        row["case_id"] for row in _read_csv(output_dir / "resume_plan.csv")
        if row["resume_status"] == "ready_not_attempted"
    }


def reconcile(
    coarse_dir: Path,
    *,
    shadow_dir: Path | None = None,
    promotion_policy: str = PROMOTION_POLICY,
    scientific_scope: str = SCIENTIFIC_SCOPE,
) -> dict[str, object]:
    """Promote only verified shadow results and build a 1554-row overlay."""

    started = time.perf_counter()
    validate_scope(scientific_scope)
    if promotion_policy != PROMOTION_POLICY:
        raise ValueError(f"unsupported promotion policy: {promotion_policy}")
    coarse_dir = coarse_dir.resolve()
    shadow_dir = (shadow_dir or coarse_dir / SHADOW_DIRECTORY_NAME).resolve()
    output_dir = coarse_dir / OUTPUT_DIRECTORY_NAME
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)

    point_cache_dir = coarse_dir / "cache"
    source_inventory = _read_csv(shadow_dir / "source_cache_inventory.csv")
    cache_fingerprint_before, cache_mismatches_before = _point_cache_fingerprint(
        point_cache_dir, source_inventory
    )
    shadow_tree_before = _directory_fingerprint(shadow_dir)
    immutable_paths = (
        coarse_dir / "grid_manifest.csv", coarse_dir / "partial_case_summary.csv",
        coarse_dir / "partial_unresolved_audit.csv",
    )
    immutable_hashes_before = {path.name: _sha256(path) for path in immutable_paths}

    shadow_metadata = json.loads((shadow_dir / "run_metadata.json").read_text(encoding="utf-8"))
    shadow_gates = {row["gate"]: row["status"] for row in _read_csv(shadow_dir / "gate_summary.csv")}
    shadow_gate_pass = all(shadow_gates.get(gate) == "PASS" for gate in REQUIRED_SHADOW_GATES)
    artifact_mismatches = _shadow_artifact_hashes(shadow_dir, shadow_metadata)
    manifest_path = coarse_dir / "grid_manifest.csv"
    manifest_sha = _sha256(manifest_path)
    expected_manifest_sha = str(
        shadow_metadata.get("original_csv_hashes_before", {}).get("grid_manifest.csv", "")
    ) if isinstance(shadow_metadata.get("original_csv_hashes_before", {}), Mapping) else ""
    scope_pass = shadow_metadata.get("scientific_scope") == SCIENTIFIC_SCOPE
    expected_cache_fingerprint = str(shadow_metadata.get("cache_fingerprint_before", ""))
    cache_pass = (
        not cache_mismatches_before
        and cache_fingerprint_before == expected_cache_fingerprint
        and cache_fingerprint_before == shadow_metadata.get("cache_fingerprint_after")
    )
    manifest_pass = manifest_sha == expected_manifest_sha
    if not shadow_input_contract_valid(
        required_gates_pass=scope_pass and shadow_gate_pass and cache_pass,
        actual_cache_fingerprint=cache_fingerprint_before,
        expected_cache_fingerprint=expected_cache_fingerprint,
        actual_manifest_sha256=manifest_sha,
        expected_manifest_sha256=expected_manifest_sha,
        artifact_mismatches=artifact_mismatches,
    ):
        raise ReconciliationIntegrityError(
            "shadow input integrity failed: "
            + _json({
                "scope": scope_pass, "gates": shadow_gate_pass, "cache": cache_pass,
                "manifest": manifest_pass, "artifact_mismatches": artifact_mismatches,
                "cache_mismatches": cache_mismatches_before,
            })
        )

    manifest_rows = _read_csv(manifest_path)
    if len(manifest_rows) != 1554 or len({row["case_id"] for row in manifest_rows}) != 1554:
        raise ReconciliationIntegrityError("approved manifest must contain 1554 unique case IDs")
    manifest_by_id = {row["case_id"]: row for row in manifest_rows}
    shadow_manifest = _read_csv(shadow_dir / "shadow_manifest.csv")
    comparisons = _read_csv(shadow_dir / "original_vs_shadow_cases.csv")
    comparison_by_id = {row["case_id"]: row for row in comparisons}
    if len(comparison_by_id) != len(comparisons) or not set(comparison_by_id) <= set(manifest_by_id):
        raise ReconciliationIntegrityError("shadow case identities are duplicated or outside manifest")
    shadow_resolved = [
        row for row in _read_csv(shadow_dir / "shadow_resolved_cases.csv")
        if row["shadow_execution_status"] == "shadow_resolved_after_local_repair"
    ]
    shadow_deferred = _read_csv(shadow_dir / "shadow_deferred_cases.csv")
    if set(row["case_id"] for row in shadow_resolved) & set(row["case_id"] for row in shadow_deferred):
        raise ReconciliationIntegrityError("a shadow case cannot be both resolved and deferred")
    spectra_rows = _read_csv(shadow_dir / "repaired_spectra.csv")
    spectra = _spectra_map(spectra_rows)
    detector_rows = _read_csv(shadow_dir / "detector_events.csv")
    candidate_rows = _read_csv(shadow_dir / "local_root_candidates.csv")
    window_rows = _read_csv(shadow_dir / "inferred_repair_windows.csv")
    promotion_timestamp = _stable_promotion_timestamp(shadow_dir)
    operation_counts = shadow_metadata.get("operation_counts", {})
    force_strict_executed = int(operation_counts.get("force_strict_executed", -1)) if isinstance(operation_counts, Mapping) else -1
    shadow_integrity = {"gates": shadow_gate_pass, "cache": cache_pass, "manifest": manifest_pass}

    eligibility_rows: list[dict[str, object]] = []
    promotions: list[dict[str, object]] = []
    for row in sorted(shadow_resolved, key=_case_sort):
        audit, promotion = _promotion_candidate(
            row, manifest_row=manifest_by_id[row["case_id"]], shadow_dir=shadow_dir,
            shadow_integrity=shadow_integrity, spectra=spectra,
            detector_rows=detector_rows, candidate_rows=candidate_rows,
            window_rows=window_rows, force_strict_executed=force_strict_executed,
            promotion_timestamp=promotion_timestamp,
        )
        eligibility_rows.append(audit)
        if promotion is not None:
            promotions.append(promotion)
    for row in sorted(shadow_deferred, key=_case_sort):
        eligibility_rows.append(_deferred_audit_row(row))
    if len({row["case_id"] for row in promotions}) != len(promotions):
        raise ReconciliationIntegrityError("duplicate promotion is forbidden")
    promotion_by_id = {str(row["case_id"]): row for row in promotions}
    deferred_by_id = {row["case_id"]: row for row in shadow_deferred}

    cache_case_ids = {
        row["case_id"] for row in source_inventory
        if row.get("case_id", "") and _truthy(row.get("selected_for_shadow", ""))
    }
    manifest_output: list[dict[str, object]] = []
    reconciled: list[dict[str, object]] = []
    comparisons_output: list[dict[str, object]] = []
    resume_rows: list[dict[str, object]] = []
    for manifest_row in manifest_rows:
        case_id = manifest_row["case_id"]
        original = comparison_by_id.get(case_id)
        promotion = promotion_by_id.get(case_id)
        deferred = deferred_by_id.get(case_id)
        point_present = case_id in cache_case_ids
        source_hash = str(original.get("source_cache_hash", "")) if original else ""
        original_status = str(original.get("original_execution_status", "not_attempted")) if original else "not_attempted"
        original_n = str(original.get("original_N_true", "")) if original else ""
        original_n_csv = original_n if original_n else "NaN"
        first_failure = str(original.get("original_first_failure", "")) if original else ""
        required_guard = str(original.get("required_guard", "")) if original else ""
        detector_status = str(original.get("family_inventory_status", "not_attempted")) if original else "not_attempted"
        local_status = str(original.get("local_repair_status", "not_attempted")) if original else "not_attempted"
        prefix_status = str(original.get("required_prefix_guard_status", "not_attempted")) if original else "not_attempted"
        upper_status = str(original.get("upper_spectrum_audit_status", "not_attempted")) if original else "not_attempted"
        full_status = str(original.get("full_K10_control_status", "not_attempted")) if original else "not_attempted"
        unresolved_reason = str(original.get("defer_reason", "")) if original else ""
        roots_eb: list[float] = []
        roots_timo: list[float] = []
        deltas: list[float] = []
        promotion_provenance = ""
        if promotion:
            reconciled_status = "resolved_after_local_repair"
            reconciled_n = str(promotion["N_true"])
            first_failure = str(promotion["first_failed_mode"])
            required_guard = str(promotion["required_guard"])
            result_origin = "family_local_repair_promotion"
            promotion_status = "verified_shadow_promoted"
            reconciliation_status = "promoted_verified_shadow"
            detector_status = "required_defect_matrix_confirmed"
            local_status = "resolved_after_local_repair"
            prefix_status = "resolved_after_local_repair"
            upper_status = str(original.get("upper_spectrum_audit_status", "not_required")) if original else "not_required"
            roots_eb = json.loads(str(promotion["EB_roots_json"]))
            roots_timo = json.loads(str(promotion["Timoshenko_roots_json"]))
            deltas = json.loads(str(promotion["delta_f_json"]))
            promotion_provenance = _json({
                "promotion_version": promotion["promotion_version"],
                "promotion_timestamp_utc": promotion["promotion_timestamp_utc"],
                "detector_version": promotion["detector_version"],
                "repair_version": promotion["repair_version"],
                "tolerance_hash": promotion["tolerance_hash"],
                "repair_cache_hash": promotion["repair_cache_hash"],
            })
        elif deferred:
            reconciled_status = "deferred_expensive_strict"
            reconciled_n = "NaN"
            result_origin = "family_local_repair_shadow_deferred"
            promotion_status = "not_promoted"
            reconciliation_status = "deferred_after_failed_local_repair"
            detector_status = str(deferred.get("family_inventory_status", detector_status))
            local_status = str(deferred.get("local_repair_status", local_status))
            prefix_status = "unresolved_without_expensive_strict"
            upper_status = "incomplete_above_required_guard"
            full_status = "not_attempted"
            unresolved_reason = str(deferred.get("defer_reason", "local_repair_did_not_resolve_required_prefix"))
        elif original and original_n:
            reconciled_status = (
                "resolved_full_K10" if original.get("shadow_execution_status") == "resolved_full_K10"
                else "resolved_primary"
            )
            reconciled_n = original_n
            result_origin = "immutable_point_cache"
            promotion_status = "not_required"
            reconciliation_status = "unchanged_original_resolved"
            roots_eb = list(spectra.get((case_id, "Euler-Bernoulli"), {}).get("before", ()))
            roots_timo = list(spectra.get((case_id, "Timoshenko"), {}).get("before", ()))
            deltas, computed_n, computed_guard = _delta_payload(roots_eb, roots_timo)
            if computed_n is not None and computed_n != int(original_n):
                raise ReconciliationIntegrityError(f"original resolved N_true mismatch: {case_id}")
            # Guard conventions are carried verbatim from the accepted point
            # cache/shadow row.  They differ intentionally between an early
            # failed-mode guard and the historical complete-K10 inventory, so
            # reconciliation validates N_true/delta but does not reinterpret
            # an already accepted guard number.
        else:
            reconciled_status = "not_attempted"
            reconciled_n = "NaN"
            result_origin = "none"
            promotion_status = "not_applicable"
            reconciliation_status = "not_attempted"
            detector_status = "not_attempted"
            local_status = "not_attempted"
            prefix_status = "not_attempted"
            upper_status = "not_attempted"
            full_status = "not_attempted"
        delta_fields = {
            f"delta_f_{index}": _format_float(deltas[index - 1]) if index <= len(deltas) else ""
            for index in range(1, K_MAX + 1)
        }
        row_out = {
            "case_id": case_id, "case_identity": manifest_row["case_identity"],
            "grid_group": manifest_row["grid_group"], "claim_eligible": manifest_row["claim_eligible"],
            "epsilon_0": manifest_row["epsilon_0"], "beta": manifest_row["beta_deg"],
            "mu": manifest_row["mu"], "eta": manifest_row["eta"],
            "s_max": manifest_row["s_max"], "thin_0p1_flag": manifest_row["thin_0p1_flag"],
            "original_execution_status": original_status,
            "reconciled_execution_status": reconciled_status,
            "original_N_true": original_n_csv, "reconciled_N_true": reconciled_n,
            "first_failed_mode": first_failure, "required_guard": required_guard,
            "result_origin": result_origin, "promotion_status": promotion_status,
            "reconciliation_status": reconciliation_status, "detector_status": detector_status,
            "local_repair_status": local_status, "required_prefix_guard_status": prefix_status,
            "upper_spectrum_audit_status": upper_status, "full_K10_control_status": full_status,
            "unresolved_reason": unresolved_reason, "point_cache_present": point_present,
            "source_cache_hash": source_hash,
            "EB_roots_json": _json(roots_eb) if roots_eb else "",
            "Timoshenko_roots_json": _json(roots_timo) if roots_timo else "",
            "delta_f_json": _json(deltas) if deltas else "", **delta_fields,
            "promotion_provenance_json": promotion_provenance,
            "scientific_scope": SCIENTIFIC_SCOPE,
        }
        reconciled.append(row_out)
        original_resolved = bool(original_n)
        comparisons_output.append({
            "case_id": case_id, "original_execution_status": original_status,
            "reconciled_execution_status": reconciled_status,
            "original_N_true": original_n_csv, "reconciled_N_true": reconciled_n,
            "N_true_changed": original_resolved and original_n != reconciled_n,
            "first_failed_mode_changed": original_resolved and first_failure != str(original.get("original_first_failure", "")),
            "root_inventory_changed": promotion is not None,
            "reconciliation_status": reconciliation_status, "result_origin": result_origin,
            "source_cache_hash": source_hash, "scientific_scope": SCIENTIFIC_SCOPE,
        })
        source_class = (
            "promoted_shadow_resolved" if promotion else "shadow_deferred" if deferred
            else "original_resolved" if original_resolved else "not_attempted"
        )
        manifest_output.append({
            "case_id": case_id, "case_identity": manifest_row["case_identity"],
            "grid_group": manifest_row["grid_group"], "epsilon_0": manifest_row["epsilon_0"],
            "beta": manifest_row["beta_deg"], "mu": manifest_row["mu"], "eta": manifest_row["eta"],
            "shadow_present": original is not None, "point_cache_present": point_present,
            "source_class": source_class, "source_cache_hash": source_hash,
            "scientific_scope": SCIENTIFIC_SCOPE,
        })
        if original_resolved:
            resume_status = "skipped_existing_resolved"
            skip_reason = "immutable point-cache result already accepted"
        elif promotion:
            resume_status = "skipped_promoted_resolved"
            skip_reason = "verified promotion overlay supplies exact result"
        elif deferred:
            resume_status = "skipped_deferred"
            skip_reason = "complex-case pass required; main pass skips deferred"
        elif original_status == "interrupted_incomplete":
            resume_status = "skipped_interrupted"
            skip_reason = "explicit --skip-interrupted policy"
        elif original is None:
            resume_status = "ready_not_attempted"
            skip_reason = ""
        else:
            resume_status = "blocked_invalid_metadata"
            skip_reason = "calculated case has neither accepted nor deferred reconciliation"
        resume_rows.append({
            "case_id": case_id, "epsilon_0": manifest_row["epsilon_0"],
            "beta": manifest_row["beta_deg"], "mu": manifest_row["mu"], "eta": manifest_row["eta"],
            "resume_status": resume_status, "skip_reason": skip_reason,
            "original_execution_status": original_status,
            "reconciled_execution_status": reconciled_status,
            "point_cache_present": point_present, "promotion_status": promotion_status,
            "scientific_scope": SCIENTIFIC_SCOPE,
        })

    reconciled.sort(key=_case_sort)
    manifest_output.sort(key=_case_sort)
    comparisons_output.sort(key=lambda row: str(row["case_id"]))
    resume_rows.sort(key=_case_sort)
    promotions.sort(key=_case_sort)
    eligibility_rows.sort(key=lambda row: (str(row["candidate_type"]), str(row["case_id"])))
    deferred_output: list[dict[str, object]] = []
    unresolved_output: list[dict[str, object]] = []
    window_by_case: dict[str, list[dict[str, str]]] = defaultdict(list)
    for item in window_rows:
        window_by_case[item["case_id"]].append(dict(item))
    candidates_by_case: dict[str, list[dict[str, str]]] = defaultdict(list)
    for item in candidate_rows:
        candidates_by_case[item["case_id"]].append(dict(item))
    for row in sorted(shadow_deferred, key=_case_sort):
        case_id = row["case_id"]
        windows_json = _json(window_by_case.get(case_id, ()))
        reason = row.get("defer_reason", "local_repair_did_not_resolve_required_prefix")
        deferred_output.append({
            "case_id": case_id, "epsilon_0": row["epsilon_0"], "beta": row["beta"],
            "mu": row["mu"], "eta": row["eta"],
            "execution_status": "deferred_expensive_strict", "N_true": "NaN",
            "first_apparent_failed_mode": row.get("original_first_failure", ""),
            "required_guard": row["required_guard"], "unresolved_reason": reason,
            "inferred_window_json": windows_json, "local_repair_status": row["local_repair_status"],
            "force_strict_requested": 1, "force_strict_executed": 0,
            "source_cache_hash": row["source_cache_hash"], "scientific_scope": SCIENTIFIC_SCOPE,
        })
        unresolved_output.append({
            "case_id": case_id, "original_execution_status": row["original_execution_status"],
            "reconciled_execution_status": "deferred_expensive_strict",
            "first_apparent_failed_mode": row.get("original_first_failure", ""),
            "required_guard": row["required_guard"],
            "detector_status": row["family_inventory_status"],
            "local_repair_status": row["local_repair_status"],
            "inferred_window_json": windows_json,
            "candidate_outcome_json": _json(candidates_by_case.get(case_id, ())),
            "unresolved_reason": reason, "future_audit_category": "targeted_complex_case_audit",
            "source_cache_hash": row["source_cache_hash"], "scientific_scope": SCIENTIFIC_SCOPE,
        })

    operation_counts_reconciliation = {
        "point_solver_calls": 0, "matrix_evaluator_calls": 0, "local_repair_calls": 0,
        "detector_calls": 0, "force_strict_requested": len(deferred_output),
        "force_strict_executed": 0, "processed_not_attempted_points": 0,
        "branch_tracking_calls": 0, "MAC_calls": 0, "shape_calls": 0,
        "FEM_calls": 0,
    }
    operation_rows = [
        {"operation": key, "count": value, "scope": "reconcile_only"}
        for key, value in sorted(operation_counts_reconciliation.items())
    ]

    _write_csv(output_dir / "reconciliation_manifest.csv", manifest_output, RECONCILIATION_MANIFEST_FIELDS)
    _write_csv(output_dir / "eligibility_audit.csv", eligibility_rows, ELIGIBILITY_FIELDS)
    _write_csv(output_dir / "promotion_overlay.csv", promotions, PROMOTION_FIELDS)
    _write_csv(output_dir / "original_vs_reconciled_cases.csv", comparisons_output, COMPARISON_FIELDS)
    _write_csv(output_dir / "reconciled_case_results.csv", reconciled, RECONCILED_FIELDS)
    _write_csv(output_dir / "promoted_cases.csv", promotions, PROMOTION_FIELDS)
    _write_csv(output_dir / "deferred_cases.csv", deferred_output, DEFERRED_FIELDS)
    _write_csv(output_dir / "unresolved_provenance.csv", unresolved_output, UNRESOLVED_FIELDS)
    _write_csv(output_dir / "resume_plan.csv", resume_rows, RESUME_FIELDS)
    _write_csv(output_dir / "operation_counts.csv", operation_rows, OPERATION_FIELDS)

    reference = _reference_validation(coarse_dir, shadow_metadata)
    changed_original_n = sum(_truthy(row["N_true_changed"]) for row in comparisons_output)
    changed_original_failure = sum(_truthy(row["first_failed_mode_changed"]) for row in comparisons_output)
    false_n = sum(row["N_true"] != "NaN" for row in deferred_output)
    eligible_count = sum(_truthy(row["eligible"]) for row in eligibility_rows)
    promotion_ids = {row["case_id"] for row in promotions}
    deferred_ids = {row["case_id"] for row in deferred_output}
    resume_counts = Counter(row["resume_status"] for row in resume_rows)
    original_resolved_count = sum(row["reconciliation_status"] == "unchanged_original_resolved" for row in reconciled)
    point_cache_missing_count = sum(not _truthy(row["point_cache_present"]) for row in reconciled)
    original_interrupted_count = sum(row["original_execution_status"] == "interrupted_incomplete" for row in reconciled)
    final_pending_count = sum(row["reconciled_execution_status"] == "pending_family_inventory_check" for row in reconciled)

    previous_metadata_path = output_dir / "run_metadata.json"
    previous_metadata = json.loads(previous_metadata_path.read_text(encoding="utf-8")) if previous_metadata_path.exists() else {}
    data_hashes = _csv_hashes(output_dir, DETERMINISTIC_DATA_CSVS)
    prior_data_hashes = previous_metadata.get("deterministic_csv_hashes", {})
    idempotence_pass = bool(prior_data_hashes) and prior_data_hashes == data_hashes

    prohibited_loaded = sorted(
        name for name in sys.modules
        if "anisotropic_rods" in name or "yartsev" in name.lower()
    )
    gate_conditions = {
        "scope_isolation_gate": scope_pass and not prohibited_loaded,
        "shadow_input_integrity_gate": shadow_gate_pass and cache_pass and manifest_pass and not artifact_mismatches,
        "promotion_eligibility_gate": eligible_count == len(promotions) and not (promotion_ids & deferred_ids) and all(_truthy(row["provenance_complete"]) for row in eligibility_rows if _truthy(row["eligible"])),
        "original_resolved_preservation_gate": changed_original_n == 0 and changed_original_failure == 0 and reference["S3_pass"] and reference["readiness_pass"] and reference["former_pass"],
        "deferred_safety_gate": false_n == 0 and operation_counts_reconciliation["force_strict_executed"] == 0,
        "no_solve_reconciliation_gate": all(operation_counts_reconciliation[key] == 0 for key in ("point_solver_calls", "matrix_evaluator_calls", "local_repair_calls", "processed_not_attempted_points")),
        "idempotence_gate": idempotence_pass,
        "production_resume_plan_gate": (
            len(resume_rows) == len(manifest_rows)
            and not (promotion_ids & {row["case_id"] for row in resume_rows if row["resume_status"] == "ready_not_attempted"})
            and not (deferred_ids & {row["case_id"] for row in resume_rows if row["resume_status"] == "ready_not_attempted"})
            and resume_counts["blocked_invalid_metadata"] == 0
            and final_pending_count == 0
        ),
    }
    gate_rows: list[dict[str, object]] = []
    gate_metrics = {
        "scope_isolation_gate": ("scope/prohibited modules", f"{scientific_scope}/{len(prohibited_loaded)}"),
        "shadow_input_integrity_gate": ("shadow gates/cache/manifest/artifacts", f"{shadow_gate_pass}/{cache_pass}/{manifest_pass}/{len(artifact_mismatches)}"),
        "promotion_eligibility_gate": ("eligible/promoted/deferred", f"{eligible_count}/{len(promotions)}/{len(deferred_ids)}"),
        "original_resolved_preservation_gate": ("changed N/failure; S3/readiness/former", f"{changed_original_n}/{changed_original_failure};{reference['S3_pass']}/{reference['readiness_count']}/{reference['former_count']}"),
        "deferred_safety_gate": ("deferred/false N/strict executed", f"{len(deferred_ids)}/{false_n}/0"),
        "no_solve_reconciliation_gate": ("point/matrix/repair/not-attempted processed", "0/0/0/0"),
        "idempotence_gate": ("prior deterministic CSV hash match", idempotence_pass),
        "production_resume_plan_gate": ("ready/point-cache-missing/pending/blocked", f"{resume_counts['ready_not_attempted']}/{point_cache_missing_count}/{final_pending_count}/{resume_counts['blocked_invalid_metadata']}"),
    }
    explanations = {
        "scope_isolation_gate": "isotropic circular metadata only; no prohibited research module loaded",
        "shadow_input_integrity_gate": "all verified shadow gates and immutable hashes agree",
        "promotion_eligibility_gate": "only matrix-confirmed shadow local-repair resolutions are promoted",
        "original_resolved_preservation_gate": "original accepted results and immutable references are unchanged",
        "deferred_safety_gate": "deferred results carry NaN, never a provisional article N_true",
        "no_solve_reconciliation_gate": "reconciliation performs file/array operations only",
        "idempotence_gate": "PASS requires a prior identical reconciliation data set",
        "production_resume_plan_gate": "ready queue excludes accepted, promoted, deferred, and interrupted results",
    }
    for gate in (
        "scope_isolation_gate", "shadow_input_integrity_gate", "promotion_eligibility_gate",
        "original_resolved_preservation_gate", "deferred_safety_gate",
        "no_solve_reconciliation_gate", "idempotence_gate", "production_resume_plan_gate",
    ):
        metric, value = gate_metrics[gate]
        gate_rows.append({
            "gate": gate, "status": "PASS" if gate_conditions[gate] else "PENDING" if gate == "idempotence_gate" else "FAIL",
            "metric": metric, "value": value, "explanation": explanations[gate],
        })
    production_pass = all(row["status"] == "PASS" for row in gate_rows)
    gate_rows.append({
        "gate": "production_resume_readiness_gate",
        "status": "PASS" if production_pass else "PENDING",
        "metric": "component gates A-H", "value": production_pass,
        "explanation": "PASS prepares but never executes the explicit future production command",
    })
    _write_csv(output_dir / "gate_summary.csv", gate_rows, GATE_FIELDS)

    cache_fingerprint_after, cache_mismatches_after = _point_cache_fingerprint(
        point_cache_dir, source_inventory
    )
    shadow_tree_after = _directory_fingerprint(shadow_dir)
    immutable_hashes_after = {path.name: _sha256(path) for path in immutable_paths}
    if (
        cache_fingerprint_after != cache_fingerprint_before
        or cache_mismatches_after
        or shadow_tree_after != shadow_tree_before
        or immutable_hashes_after != immutable_hashes_before
    ):
        raise ReconciliationIntegrityError("an immutable source changed during reconciliation")

    final_all_hashes = _csv_hashes(output_dir, ALL_CSVS)
    previous_all_hashes = previous_metadata.get("all_csv_sha256", {})
    all_csv_repeat_match = bool(previous_all_hashes) and previous_all_hashes == final_all_hashes
    future_command = (
        "D:\\python\\Pycharm\\pythonProject\\.venv\\Scripts\\python.exe -B "
        "scripts/analysis/thickness_mismatch/audits/run_article_epsilon_upper_envelope_grid.py "
        "--output-dir results/article_epsilon_upper_envelope/coarse_grid_v1 "
        "--prefix-until-failure --prefix-strategy paired --strict-policy main-pass "
        "--family-inventory-policy local-repair --defer-expensive-strict --workers 4 "
        "--reuse-cache --main-pass-only --skip-existing-unresolved --skip-deferred --skip-interrupted"
    )
    metadata = {
        "reconciliation_version": RECONCILIATION_VERSION,
        "promotion_policy": promotion_policy, "scientific_scope": scientific_scope,
        "source_shadow_directory": shadow_dir.relative_to(coarse_dir.parents[2]).as_posix(),
        "canonical_article_facing_table": (output_dir / "reconciled_case_results.csv").relative_to(coarse_dir.parents[2]).as_posix(),
        "promotion_overlay": (output_dir / "promotion_overlay.csv").relative_to(coarse_dir.parents[2]).as_posix(),
        "source_backup_metadata": {"immutable_hashes_before": immutable_hashes_before, "shadow_tree_before": shadow_tree_before},
        "source_cache_fingerprint_before": cache_fingerprint_before,
        "source_cache_fingerprint_after": cache_fingerprint_after,
        "manifest_sha256_before": manifest_sha, "manifest_sha256_after": _sha256(manifest_path),
        "shadow_tree_fingerprint_before": shadow_tree_before,
        "shadow_tree_fingerprint_after": shadow_tree_after,
        "immutable_hashes_before": immutable_hashes_before,
        "immutable_hashes_after": immutable_hashes_after,
        "counts": {
            "manifest": len(manifest_rows), "shadow": len(comparisons),
            "original_resolved": original_resolved_count, "eligible_promotions": eligible_count,
            "promoted": len(promotions), "deferred": len(deferred_ids),
            "point_cache_missing": point_cache_missing_count,
            "original_interrupted": original_interrupted_count,
            "remaining_interrupted": resume_counts["skipped_interrupted"],
            "ready_not_attempted": resume_counts["ready_not_attempted"],
        },
        "operation_counts": operation_counts_reconciliation,
        "reference_validation": reference,
        "gates": {row["gate"]: row["status"] for row in gate_rows},
        "deterministic_csv_hashes": data_hashes,
        "all_csv_sha256": final_all_hashes,
        "all_csv_repeat_match": all_csv_repeat_match,
        "promotion_timestamp_utc": promotion_timestamp,
        "current_invocation": {"wall_seconds": time.perf_counter() - started, "root_calculations": 0},
        "future_production_command": future_command,
        "future_production_command_executed": False,
    }
    _atomic_json(output_dir / "run_metadata.json", metadata)
    report_lines = [
        "# Family local-repair reconciliation", "",
        f"Scientific scope: `{scientific_scope}`.",
        "The point cache and verified shadow directory are immutable inputs. No solver, detector, or local repair was called.",
        "The reconciled table compares sorted spectral positions and does not use descendant branches.", "",
        "## Counts", "",
        f"- manifest/shadow: {len(manifest_rows)}/{len(comparisons)}",
        f"- unchanged original resolved/promoted/deferred: {original_resolved_count}/{len(promotions)}/{len(deferred_ids)}",
        f"- point-cache missing/ready queue: {point_cache_missing_count}/{resume_counts['ready_not_attempted']}",
        f"- original/remaining interrupted: {original_interrupted_count}/{resume_counts['skipped_interrupted']}", "",
        "## Gates", "",
    ]
    report_lines.extend(f"- {row['gate']}: {row['status']}" for row in gate_rows)
    report_lines.extend(["", "## Future production command (not executed)", "", f"`{future_command}`", ""])
    (output_dir / "report.md").write_text("\n".join(report_lines), encoding="utf-8", newline="\n")
    log_path = output_dir / "logs" / f"reconcile_{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')}.json"
    _atomic_json(log_path, {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "wall_seconds": time.perf_counter() - started, "root_calculations": 0,
        "data_hash_repeat_match": idempotence_pass,
        "all_csv_repeat_match": all_csv_repeat_match,
        "gates": metadata["gates"],
    })
    return {
        "output_dir": output_dir, "root_calculations": 0,
        "promoted_count": len(promotions), "deferred_count": len(deferred_ids),
        "ready_count": resume_counts["ready_not_attempted"],
        "gates": metadata["gates"], "metadata": metadata,
    }


__all__ = [
    "OUTPUT_DIRECTORY_NAME", "PROMOTION_POLICY", "RECONCILIATION_VERSION",
    "ReconciliationIntegrityError", "SCIENTIFIC_SCOPE", "deterministic_family_poststage",
    "confirmed_n_true_value", "forbid_scientific_operation", "load_ready_resume_ids",
    "pending_family_status", "reconcile", "shadow_input_contract_valid", "validate_scope",
]
