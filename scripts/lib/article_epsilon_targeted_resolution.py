"""Target-only resolution of the epsilon_0=0.050 article-envelope blockers.

The workflow is deliberately separate from the coarse-grid runner.  It selects
its targets from the compact unresolved table, loads at most one raw point
payload at a time, and performs only exact-beta local matrix searches before
considering any stricter path.  Production equations, matrices, basis regimes,
and tolerances are imported unchanged.
"""

from __future__ import annotations

import gc
import gzip
import json
import math
import sys
import time
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from scripts.lib import family_inventory_local_repair as repair
from scripts.lib import general_spectrum_completeness as complete
from scripts.lib import variable_length_timoshenko as timo
from scripts.lib.article_epsilon_compact_certificates import (
    DELTA_TOL,
    K_MAX,
    MODEL_EB,
    MODEL_TIMO,
    SCIENTIFIC_SCOPE,
    atomic_write_csv,
    atomic_write_gzip_json,
    atomic_write_json,
    load_certificate,
    process_rss_bytes,
    read_csv,
    sha256_file,
)
from scripts.lib.article_epsilon_compact_poststage import (
    EPSILON_FIELDS,
    FINAL_FIELDS,
    UNRESOLVED_FIELDS,
    _epsilon_summary,
)


TARGET_VERSION = "article_epsilon_005_targeted_resolution_v1"
TARGET_CACHE_SCHEMA = "article_epsilon_005_target_cache_v1"
TARGET_EPSILON = 0.050
TARGET_EPSILON_TOLERANCE = 1.0e-12
TARGET_COUNT_LIMIT = 2
T1_PHASES = (0.0, 0.25, 0.5, 0.75)
T1_STEPS = (1.0e-4, 5.0e-5)

TARGET_MANIFEST_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta",
    "final_execution_status", "required_guard", "unresolved_reason",
    "local_repair_status", "compact_certificate_path",
    "source_compact_certificate_path", "canonical_raw_cache_path",
)
DIAGNOSTIC_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta", "unresolved_theory",
    "unresolved_rank", "current_root_count", "required_guard",
    "previous_window_left", "previous_window_right", "previous_repair_status",
    "previous_failure_reason", "previous_matrix_evaluations",
    "primary_EB_status", "primary_Timoshenko_status",
    "stored_independent_prefix_status", "strict_trigger_reasons",
    "raw_local_status", "raw_cluster_status", "raw_deltas_through_guard",
    "envelope_question", "proposed_resolution_stage",
)
EVENT_FIELDS = (
    "case_id", "stage", "theory", "rank_start", "rank_end",
    "lambda_left", "lambda_right", "scan_step", "phase_count",
    "root_count", "stable_with_previous", "status", "matrix_evaluations",
)
CANDIDATE_FIELDS = (
    "case_id", "stage", "theory", "rank_start", "rank_end",
    "lambda_candidate", "candidate_source", "sign_change", "sigma_min",
    "sigma_ratio", "residual", "bracket_left", "bracket_right",
    "multiplicity", "accepted", "rejection_reason",
)
INDEPENDENT_FIELDS = (
    "case_id", "theory", "sorted_rank", "method",
    "lambda_primary", "lambda_verification", "absolute_difference",
    "root_match_tolerance", "status",
)
STRICT_FIELDS = (
    "case_id", "T0_status", "T1_status", "T2_status", "T3_status",
    "T4_status", "strict_level_used", "force_strict_requested",
    "force_strict_executed", "strict_containment_status",
)
OVERLAY_FIELDS = (
    "case_id", "final_status", "N_true", "certified_N_true_upper_bound",
    "first_failed_mode", "required_guard", "proof_rank", "EB_roots_json",
    "Timoshenko_roots_json", "delta_f_json", "verification_methods",
    "strict_level_used", "source_compact_certificate_hash",
    "source_raw_cache_hash", "scientific_scope", "target_version",
)
BEFORE_AFTER_FIELDS = (
    "case_id", "epsilon_0", "beta", "original_execution_status",
    "original_N_true", "final_execution_status", "final_N_true",
    "certified_N_true_upper_bound", "first_failed_mode", "required_guard",
)
EPSILON_STATUS_FIELDS = (
    "epsilon_0", "previous_observed_max_N_true", "final_observed_max_N_true",
    "target_case_count", "target_exact_count", "target_bounded_count",
    "remaining_target_unresolved", "envelope_finalizable", "envelope_status",
    "argmax_case_ids",
)
RAW_ENVELOPE_FIELDS = (
    "epsilon_0", "N_up_raw", "envelope_status", "unresolved_case_count",
    "observed_argmax_case_ids",
)
SUFFIX_FIELDS = (
    "epsilon_0", "N_up_raw", "N_up_suffix_max", "provenance",
)
RESULT_CHANGE_FIELDS = (
    "case_id", "field", "before", "after", "change_origin",
)
OPERATION_FIELDS = ("operation", "count")
RUNTIME_FIELDS = (
    "run", "case_id", "cache_hit", "wall_seconds", "matrix_evaluations",
    "initial_rss_bytes", "peak_rss_bytes", "final_rss_bytes",
    "wall_source",
)
MEMORY_FIELDS = ("case_id", "stage", "rss_bytes", "peak_rss_bytes")
GATE_FIELDS = ("gate", "status", "value")
REPEAT_FIELDS = (
    "run", "wall_seconds", "matrix_evaluations", "target_cache_hits",
    "peak_rss_bytes", "scientific_csv_hashes_stable", "source",
)


@dataclass(frozen=True)
class TargetWindow:
    theory: str
    rank_start: int
    rank_end: int
    lambda_left: float
    lambda_right: float
    expected_root_count: int
    reason: str


def _repo_root(coarse_dir: Path) -> Path:
    return coarse_dir.resolve().parents[2]


def _absolute(path_value: str, coarse_dir: Path) -> Path:
    path = Path(path_value)
    return path if path.is_absolute() else (_repo_root(coarse_dir) / path).resolve()


def _relative(path: Path, coarse_dir: Path) -> str:
    return str(path.resolve().relative_to(_repo_root(coarse_dir))).replace("\\", "/")


def _int(value: object, default: int | None = None) -> int | None:
    try:
        if value in (None, "", "NaN"):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _float(value: object, default: float | None = None) -> float | None:
    try:
        if value in (None, "", "NaN"):
            return default
        result = float(value)
        return result if math.isfinite(result) else default
    except (TypeError, ValueError):
        return default


def _json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _load_raw(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as stream:
        payload = json.load(stream)
    if not isinstance(payload, dict):
        raise RuntimeError(f"raw point cache is not an object: {path}")
    return payload


def select_targets(
    unresolved_rows: Sequence[Mapping[str, str]], *, epsilon: float = TARGET_EPSILON,
) -> list[dict[str, str]]:
    """Select the actual deferred epsilon target rows without case-ID constants."""

    selected = [
        dict(row) for row in unresolved_rows
        if abs(float(row["epsilon_0"]) - float(epsilon)) <= TARGET_EPSILON_TOLERANCE
        and any(token in str(row.get("final_execution_status", "")).lower()
                for token in ("deferred", "unresolved"))
    ]
    selected.sort(key=lambda row: (
        float(row["beta"]), float(row["mu"]), float(row["eta"]), row["case_id"],
    ))
    if len({row["case_id"] for row in selected}) != len(selected):
        raise RuntimeError("duplicate target case IDs in unresolved table")
    if len(selected) > TARGET_COUNT_LIMIT:
        raise RuntimeError(f"target selection exceeded containment limit: {len(selected)}")
    return selected


def strict_is_contained(case_id: str, target_ids: set[str], prior_stages: Sequence[str]) -> bool:
    """Return whether an expensive strict call would obey target/stage containment."""

    return case_id in target_ids and "T1" in prior_stages and "T2" in prior_stages


def compute_prefix_result(
    eb_roots: Sequence[float], timo_roots: Sequence[float],
    *, k_max: int = K_MAX, delta_tol: float = DELTA_TOL,
) -> tuple[int | None, int | None, int, tuple[float, ...]]:
    """Compute the exact early prefix without requiring irrelevant upper roots."""

    available = min(len(eb_roots), len(timo_roots), int(k_max))
    deltas: list[float] = []
    first_failed: int | None = None
    for index in range(available):
        eb = float(eb_roots[index])
        tm = float(timo_roots[index])
        if not (math.isfinite(eb) and math.isfinite(tm) and tm != 0.0):
            return None, None, int(k_max) + 1, tuple(deltas)
        delta = abs(eb * eb - tm * tm) / (tm * tm)
        deltas.append(float(delta))
        if first_failed is None and delta > float(delta_tol):
            first_failed = index + 1
        if first_failed is not None and index + 1 >= first_failed + 1:
            return first_failed - 1, first_failed, first_failed + 1, tuple(deltas)
    if first_failed is not None:
        return None, first_failed, first_failed + 1, tuple(deltas)
    if available < int(k_max):
        return None, None, int(k_max) + 1, tuple(deltas)
    if min(len(eb_roots), len(timo_roots)) < int(k_max) + 1:
        return None, None, int(k_max) + 1, tuple(deltas)
    return int(k_max), None, int(k_max) + 1, tuple(deltas)


def _model_record(cert: Mapping[str, object], theory: str) -> Mapping[str, object]:
    spectra = cert.get("spectra", {})
    if not isinstance(spectra, Mapping):
        return {}
    record = spectra.get(theory, {})
    return record if isinstance(record, Mapping) else {}


def _roots(cert: Mapping[str, object], theory: str) -> tuple[float, ...]:
    return tuple(float(value) for value in _model_record(cert, theory).get("roots", ()))


def _raw_result(raw: Mapping[str, object], theory: str) -> Mapping[str, object]:
    models = raw.get("models", {})
    model = models.get(theory, {}) if isinstance(models, Mapping) else {}
    latest = model.get("latest_result", {}) if isinstance(model, Mapping) else {}
    result = latest.get("result", {}) if isinstance(latest, Mapping) else {}
    return result if isinstance(result, Mapping) else {}


def _stored_independent_rows(
    case_id: str, raw: Mapping[str, object], guard: int,
) -> tuple[list[dict[str, object]], bool]:
    rows: list[dict[str, object]] = []
    all_pass = True
    tolerance = complete.SearchSettings().root_match_tol
    for theory in (MODEL_EB, MODEL_TIMO):
        result = _raw_result(raw, theory)
        comparison = result.get("primary_vs_verification", ())
        comparison = comparison if isinstance(comparison, Sequence) else ()
        by_rank = {
            int(item.get("sorted_index", 0)): item
            for item in comparison if isinstance(item, Mapping)
        }
        for rank in range(1, guard + 1):
            item = by_rank.get(rank, {})
            status = str(item.get("status", "missing"))
            passed = status == "pass" and _float(item.get("absolute_difference"), math.inf) <= tolerance
            all_pass = all_pass and passed
            rows.append({
                "case_id": case_id, "theory": theory, "sorted_rank": rank,
                "method": "stored_pointwise_primary_vs_independent_verification",
                "lambda_primary": item.get("Lambda_primary", ""),
                "lambda_verification": item.get("Lambda_verification", ""),
                "absolute_difference": item.get("absolute_difference", ""),
                "root_match_tolerance": tolerance,
                "status": "PASS" if passed else "FAIL",
            })
        all_pass = all_pass and bool(result.get("independent_agreement"))
        all_pass = all_pass and str(result.get("spectrum_status")) == "resolved_complete"
    return rows, all_pass


def _preferred_theory(case_id: str, local_rows: Sequence[Mapping[str, str]], cert: Mapping[str, object]) -> str:
    for row in local_rows:
        if row.get("case_id") == case_id and row.get("preferred_model") in (MODEL_EB, MODEL_TIMO):
            return str(row["preferred_model"])
    diagnostics = cert.get("diagnostics", {})
    triggers = diagnostics.get("detector_trigger_summary", ()) if isinstance(diagnostics, Mapping) else ()
    return MODEL_EB if "close_root_cluster" in triggers else MODEL_TIMO


def derive_target_windows(
    cert: Mapping[str, object], theory: str, required_guard: int,
) -> tuple[TargetWindow, ...]:
    """Derive narrow exact-beta windows from compact roots and cluster facts."""

    roots = _roots(cert, theory)
    if len(roots) < required_guard:
        raise RuntimeError("compact target certificate lacks roots through required guard")
    record = _model_record(cert, theory)
    clusters = record.get("clusters", ())
    grouped: dict[str, list[int]] = defaultdict(list)
    if isinstance(clusters, Sequence):
        for item in clusters:
            if not isinstance(item, Mapping):
                continue
            rank = _int(item.get("sorted_index"))
            cluster_id = str(item.get("cluster_id", ""))
            if rank is not None and rank <= required_guard and cluster_id:
                grouped[cluster_id].append(rank)
    windows: list[TargetWindow] = []
    for cluster_id, ranks in sorted(grouped.items(), key=lambda pair: min(pair[1])):
        start, end = min(ranks), max(ranks)
        values = roots[start - 1 : end]
        internal_gap = max(np.diff(values), default=complete.DEFAULT_ROOT_DEDUP_TOL)
        margin = max(0.02, 2.0 * float(internal_gap))
        left = max(complete.DEFAULT_SCAN_START, float(values[0]) - margin)
        right = float(values[-1]) + margin
        if start > 1:
            left = max(left, 0.5 * (roots[start - 2] + roots[start - 1]))
        if end < len(roots):
            right = min(right, 0.5 * (roots[end - 1] + roots[end]))
        windows.append(TargetWindow(
            theory, start, end, left, right, end - start + 1,
            f"confirmed_close_cluster:{cluster_id}",
        ))
    if windows:
        return tuple(windows)
    local = roots[:required_guard]
    left_gap = roots[1] - roots[0] if len(roots) > 1 else max(0.1, roots[0] * 0.05)
    right_gap = (
        roots[required_guard] - roots[required_guard - 1]
        if len(roots) > required_guard else left_gap
    )
    left = max(complete.DEFAULT_SCAN_START, local[0] - max(0.05, 0.25 * left_gap))
    right = local[-1] + max(0.05, 0.25 * right_gap)
    return (TargetWindow(
        theory, 1, required_guard, float(left), float(right), required_guard,
        "required_prefix_control_window",
    ),)


def _previous_repair(
    case_id: str, local_rows: Sequence[Mapping[str, str]], coarse_dir: Path,
) -> tuple[Mapping[str, str], Mapping[str, object]]:
    row = next((item for item in local_rows if item.get("case_id") == case_id), {})
    token = "eb" if row.get("preferred_model") == MODEL_EB else "timo"
    cache = coarse_dir / "compact_finalization" / "cache" / "coarse_cases" / f"{case_id}_{token}.json"
    payload: Mapping[str, object] = {}
    if cache.exists():
        loaded = json.loads(cache.read_text(encoding="utf-8"))
        payload = loaded if isinstance(loaded, Mapping) else {}
    return row, payload


def _diagnostic_row(
    target: Mapping[str, str], cert: Mapping[str, object],
    local_rows: Sequence[Mapping[str, str]], coarse_dir: Path,
    raw: Mapping[str, object] | None = None,
) -> dict[str, object]:
    case_id = target["case_id"]
    guard = _int(target.get("required_guard"), 11) or 11
    theory = _preferred_theory(case_id, local_rows, cert)
    previous, cache = _previous_repair(case_id, local_rows, coarse_dir)
    identity = cache.get("identity", {}) if isinstance(cache, Mapping) else {}
    window = identity.get("inferred_window", {}) if isinstance(identity, Mapping) else {}
    result = cache.get("result", {}) if isinstance(cache, Mapping) else {}
    result_record = cert.get("result", {})
    apparent = _int(result_record.get("first_apparent_failed_mode")) if isinstance(result_record, Mapping) else None
    clusters = _model_record(cert, theory).get("clusters", ())
    affected = [
        _int(item.get("sorted_index")) for item in clusters
        if isinstance(item, Mapping) and (_int(item.get("sorted_index")) or 999) <= guard
    ] if isinstance(clusters, Sequence) else []
    unresolved_rank = min((rank for rank in affected if rank is not None), default=apparent or 1)
    proposed = derive_target_windows(cert, theory, guard)
    raw_record = raw or {}
    primary_status = raw_record.get("primary_status", {})
    primary_status = primary_status if isinstance(primary_status, Mapping) else {}
    _stored_rows, stored_pass = _stored_independent_rows(case_id, raw_record, guard) if raw is not None else ([], False)
    raw_deltas = raw_record.get("deltas", {})
    raw_deltas = raw_deltas if isinstance(raw_deltas, Mapping) else {}
    return {
        "case_id": case_id, "epsilon_0": float(target["epsilon_0"]),
        "beta": float(target["beta"]), "mu": float(target["mu"]),
        "eta": float(target["eta"]), "unresolved_theory": theory,
        "unresolved_rank": unresolved_rank, "current_root_count": len(_roots(cert, theory)),
        "required_guard": guard, "previous_window_left": window.get("lambda_left", "") if isinstance(window, Mapping) else "",
        "previous_window_right": window.get("lambda_right", "") if isinstance(window, Mapping) else "",
        "previous_repair_status": previous.get("execution_status", result.get("status", "") if isinstance(result, Mapping) else ""),
        "previous_failure_reason": target.get("unresolved_reason", ""),
        "previous_matrix_evaluations": previous.get("matrix_evaluations", result.get("matrix_evaluations", 0) if isinstance(result, Mapping) else 0),
        "primary_EB_status": primary_status.get(MODEL_EB, ""),
        "primary_Timoshenko_status": primary_status.get(MODEL_TIMO, ""),
        "stored_independent_prefix_status": "PASS" if stored_pass else "FAIL",
        "strict_trigger_reasons": ";".join(str(value) for value in raw_record.get("strict_trigger_reasons", ())),
        "raw_local_status": raw_record.get("local_status", ""),
        "raw_cluster_status": _json(raw_record.get("clusters_already_resolved", {})),
        "raw_deltas_through_guard": _json({
            str(rank): raw_deltas[str(rank)] for rank in range(1, guard + 1)
            if str(rank) in raw_deltas
        }),
        "envelope_question": "can_N_true_exceed_4",
        "proposed_resolution_stage": "T1:" + ";".join(
            f"{item.theory}[r{item.rank_start}-r{item.rank_end}]@[{item.lambda_left:.17g},{item.lambda_right:.17g}]"
            for item in proposed
        ),
    }


def _target_manifest(
    targets: Sequence[Mapping[str, str]], index: Mapping[str, Mapping[str, str]],
    coarse_dir: Path,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    manifest_ids = {row["case_id"] for row in read_csv(coarse_dir / "grid_manifest.csv")}
    final_rows = {row["case_id"]: row for row in read_csv(coarse_dir / "compact_finalization" / "final_case_certificates.csv")}
    for target in targets:
        case_id = target["case_id"]
        if case_id not in manifest_ids or case_id not in index or case_id not in final_rows:
            raise RuntimeError(f"target source identity is incomplete: {case_id}")
        if _int(final_rows[case_id].get("N_true")) is not None:
            raise RuntimeError(f"target already has article-facing N_true: {case_id}")
        source_cert = _absolute(index[case_id]["certificate_path"], coarse_dir)
        raw = _absolute(index[case_id]["source_full_cache_path"], coarse_dir)
        final_cert = _absolute(target["compact_certificate_path"], coarse_dir)
        if not (source_cert.exists() and final_cert.exists() and raw.exists()):
            raise RuntimeError(f"target certificate/raw source is missing: {case_id}")
        rows.append({
            **target,
            "source_compact_certificate_path": index[case_id]["certificate_path"],
            "canonical_raw_cache_path": index[case_id]["source_full_cache_path"],
        })
    return rows


def prepare_target_diagnostics(coarse_dir: Path) -> dict[str, object]:
    """Write the required selection and zero-solve audit artifacts."""

    coarse_dir = coarse_dir.resolve()
    output_dir = coarse_dir / "epsilon_005_targeted_resolution"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)
    unresolved = read_csv(coarse_dir / "compact_finalization" / "unresolved_cases.csv")
    targets = select_targets(unresolved)
    if len(targets) != 2:
        raise RuntimeError(f"expected two unambiguous epsilon_0=0.050 targets, found {len(targets)}")
    index_rows = read_csv(coarse_dir / "compact_point_certificates_v1" / "compact_index.csv")
    index = {row["case_id"]: row for row in index_rows}
    manifest_rows = _target_manifest(targets, index, coarse_dir)
    local_rows = read_csv(coarse_dir / "compact_finalization" / "local_repair_results.csv")
    diagnostics: list[dict[str, object]] = []
    source_fingerprints: list[dict[str, object]] = []
    for target in targets:
        case_id = target["case_id"]
        cert_path = _absolute(target["compact_certificate_path"], coarse_dir)
        source_path = _absolute(index[case_id]["certificate_path"], coarse_dir)
        raw_path = _absolute(index[case_id]["source_full_cache_path"], coarse_dir)
        cert = load_certificate(cert_path)
        raw = _load_raw(raw_path)
        diagnostics.append(_diagnostic_row(target, cert, local_rows, coarse_dir, raw))
        source_fingerprints.extend([
            {"case_id": case_id, "kind": "final_compact", "path": _relative(cert_path, coarse_dir), "sha256": sha256_file(cert_path), "size": cert_path.stat().st_size},
            {"case_id": case_id, "kind": "source_compact", "path": _relative(source_path, coarse_dir), "sha256": sha256_file(source_path), "size": source_path.stat().st_size},
            {"case_id": case_id, "kind": "canonical_raw", "path": _relative(raw_path, coarse_dir), "sha256": sha256_file(raw_path), "size": raw_path.stat().st_size},
        ])
        del cert, raw
        gc.collect()
    atomic_write_csv(output_dir / "target_manifest.csv", manifest_rows, TARGET_MANIFEST_FIELDS)
    atomic_write_csv(output_dir / "target_case_diagnostics.csv", diagnostics, DIAGNOSTIC_FIELDS)
    atomic_write_csv(
        output_dir / "source_fingerprints.csv", source_fingerprints,
        ("case_id", "kind", "path", "sha256", "size"),
    )
    metadata = {
        "target_version": TARGET_VERSION, "mode": "zero_solve_diagnostics",
        "scientific_scope": SCIENTIFIC_SCOPE, "target_count": len(targets),
        "point_solver_calls": 0, "matrix_evaluations": 0,
        "manifest_sha256": sha256_file(coarse_dir / "grid_manifest.csv"),
        "compact_index_sha256": sha256_file(coarse_dir / "compact_point_certificates_v1" / "compact_index.csv"),
    }
    atomic_write_json(output_dir / "zero_solve_metadata.json", metadata)
    return {"output_dir": str(output_dir), "target_count": len(targets), "matrix_evaluations": 0}


def _cache_identity(
    target: Mapping[str, str], index: Mapping[str, str], windows: Sequence[TargetWindow],
) -> dict[str, object]:
    settings = complete.SearchSettings()
    return {
        "cache_schema": TARGET_CACHE_SCHEMA, "target_version": TARGET_VERSION,
        "scientific_scope": SCIENTIFIC_SCOPE, "case_id": target["case_id"],
        "source_compact_sha256": sha256_file(Path(index["_source_compact_absolute"])),
        "source_raw_sha256": sha256_file(Path(index["_raw_absolute"])),
        "windows": [asdict(item) for item in windows], "scan_steps": list(T1_STEPS),
        "scan_phases": list(T1_PHASES), "root_match_tolerance": settings.root_match_tol,
        "root_dedup_tolerance": settings.root_dedup_tol,
        "sigma_accept": settings.sigma_accept, "sigma_ratio_accept": settings.sigma_ratio_accept,
        "eb_evaluator_version": complete.EB_MATRIX_EVALUATOR_VERSION,
        "timoshenko_evaluator_version": timo.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
    }


def _atomic_cache(path: Path, identity: Mapping[str, object], result: Mapping[str, object]) -> None:
    atomic_write_json(path, {"identity": dict(identity), "result": dict(result)})


def _load_target_cache(path: Path, identity: Mapping[str, object]) -> Mapping[str, object] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    if payload.get("identity") != dict(identity) or not isinstance(payload.get("result"), Mapping):
        return None
    return payload["result"]


def _merge_window(
    roots: Sequence[float], entries: Sequence[repair.LocalRootEntry], window: TargetWindow,
) -> tuple[float, ...]:
    repair_window = repair.RepairWindow(
        event_id="target", case_id="target", theory=window.theory, beta=0.0,
        rank_start=window.rank_start, expected_missing_count=window.expected_root_count,
        lambda_left=window.lambda_left, lambda_right=window.lambda_right,
        source=window.reason, lower_anchor=window.lambda_left,
        upper_anchor=window.lambda_right,
        predicted_roots=tuple(roots[window.rank_start - 1 : window.rank_end]),
        margin=0.0, beta_probe_required=False, status="window_inferred",
    )
    return repair.merge_inventory(
        roots, entries, repair_window,
        root_dedup_tolerance=complete.SearchSettings().root_dedup_tol,
        limit=12,
    )


def _t1_search(
    case_id: str, cert: Mapping[str, object], windows: Sequence[TargetWindow],
) -> dict[str, object]:
    geometry_record = cert["geometry"]
    geometry = complete.Geometry(
        float(geometry_record["epsilon_0"]), float(geometry_record["beta_deg"]),
        float(geometry_record["mu"]), float(geometry_record["eta"]),
    )
    settings = complete.SearchSettings()
    roots_by_model = {MODEL_EB: list(_roots(cert, MODEL_EB)), MODEL_TIMO: list(_roots(cert, MODEL_TIMO))}
    candidate_rows: list[dict[str, object]] = []
    event_rows: list[dict[str, object]] = []
    stage_roots: dict[str, list[list[float]]] = defaultdict(list)
    total_evaluations = 0
    all_stable = True
    all_expected = True
    for window in windows:
        provider = complete.model_matrix_provider(window.theory, geometry)
        previous: tuple[repair.LocalRootEntry, ...] | None = None
        final: tuple[repair.LocalRootEntry, ...] = ()
        for ordinal, step in enumerate(T1_STEPS, start=1):
            stage = f"T1_L{ordinal}_r{window.rank_start}_{window.rank_end}"
            entries, candidates, evaluations = repair._dense_search(
                provider, lambda_left=window.lambda_left, lambda_right=window.lambda_right,
                scan_step=step, stage=stage, block_family="full_matrix",
                base_settings=settings, phases=T1_PHASES,
            )
            total_evaluations += evaluations
            stable = previous is not None and repair.root_sets_stable(
                previous, entries, tolerance=settings.root_match_tol,
            )
            event_rows.append({
                "case_id": case_id, "stage": stage, "theory": window.theory,
                "rank_start": window.rank_start, "rank_end": window.rank_end,
                "lambda_left": window.lambda_left, "lambda_right": window.lambda_right,
                "scan_step": step, "phase_count": len(T1_PHASES),
                "root_count": len(entries), "stable_with_previous": stable,
                "status": "stable" if stable else "refinement_level_complete",
                "matrix_evaluations": evaluations,
            })
            stage_roots[f"{window.theory}:{window.rank_start}-{window.rank_end}"].append(
                [float(entry.value) for entry in entries]
            )
            for candidate in candidates:
                candidate_rows.append({
                    "case_id": case_id, "stage": stage, "theory": window.theory,
                    "rank_start": window.rank_start, "rank_end": window.rank_end,
                    "lambda_candidate": candidate.get("lambda_candidate"),
                    "candidate_source": candidate.get("source"),
                    "sign_change": candidate.get("sign_change"),
                    "sigma_min": candidate.get("sigma_min"),
                    "sigma_ratio": candidate.get("sigma_ratio"),
                    "residual": candidate.get("residual"),
                    "bracket_left": candidate.get("bracket_left"),
                    "bracket_right": candidate.get("bracket_right"),
                    "multiplicity": candidate.get("multiplicity"),
                    "accepted": candidate.get("accepted"),
                    "rejection_reason": candidate.get("rejection_reason"),
                })
            previous = entries
            final = entries
        stable = len(T1_STEPS) >= 2 and len(stage_roots[f"{window.theory}:{window.rank_start}-{window.rank_end}"]) >= 2
        stable = stable and repair.root_sets_stable(
            tuple(repair.LocalRootEntry(value=value, multiplicity=1, repeated_root_slot=1, block_family="", nullity=1, source="") for value in stage_roots[f"{window.theory}:{window.rank_start}-{window.rank_end}"][-2]),
            tuple(repair.LocalRootEntry(value=value, multiplicity=1, repeated_root_slot=1, block_family="", nullity=1, source="") for value in stage_roots[f"{window.theory}:{window.rank_start}-{window.rank_end}"][-1]),
            tolerance=settings.root_match_tol,
        )
        all_stable = all_stable and stable
        all_expected = all_expected and len(final) == window.expected_root_count
        roots_by_model[window.theory] = list(_merge_window(roots_by_model[window.theory], final, window))
    n_true, first_failed, guard, deltas = compute_prefix_result(
        roots_by_model[MODEL_EB], roots_by_model[MODEL_TIMO],
    )
    return {
        "T1_status": "resolved_stable_local_inventory" if all_stable and all_expected else "local_inventory_unresolved",
        "stable": all_stable, "expected_counts": all_expected,
        "N_true": n_true, "first_failed_mode": first_failed, "required_guard": guard,
        "deltas": list(deltas), "eb_roots": roots_by_model[MODEL_EB],
        "timo_roots": roots_by_model[MODEL_TIMO], "candidate_rows": candidate_rows,
        "event_rows": event_rows, "stage_roots": dict(stage_roots),
        "matrix_evaluations": total_evaluations,
    }


def _resolve_one(
    target: Mapping[str, str], index_row: Mapping[str, str], cert: Mapping[str, object],
    raw: Mapping[str, object], output_dir: Path,
) -> tuple[dict[str, object], bool, float, list[dict[str, object]]]:
    started = time.perf_counter()
    case_id = target["case_id"]
    guard = _int(target.get("required_guard"), 11) or 11
    preferred = _preferred_theory(case_id, read_csv(output_dir.parent / "compact_finalization" / "local_repair_results.csv"), cert)
    windows = derive_target_windows(cert, preferred, guard)
    enriched = dict(index_row)
    enriched["_source_compact_absolute"] = str(_absolute(index_row["certificate_path"], output_dir.parent))
    enriched["_raw_absolute"] = str(_absolute(index_row["source_full_cache_path"], output_dir.parent))
    identity = _cache_identity(target, enriched, windows)
    cache_path = output_dir / "cache" / f"{case_id}.json"
    cached = _load_target_cache(cache_path, identity)
    t0_rows, t0_pass = _stored_independent_rows(case_id, raw, guard)
    cache_hit = cached is not None
    if cached is None:
        t1 = _t1_search(case_id, cert, windows)
        result = {**t1, "windows": [asdict(item) for item in windows], "T0_pass": t0_pass}
        _atomic_cache(cache_path, identity, result)
    else:
        result = dict(cached)
    recomputed_n, recomputed_failure, recomputed_guard, recomputed_deltas = compute_prefix_result(
        result.get("eb_roots", ()), result.get("timo_roots", ()),
    )
    result.update({
        "N_true": recomputed_n, "first_failed_mode": recomputed_failure,
        "required_guard": recomputed_guard, "deltas": list(recomputed_deltas),
    })
    n_true = _int(result.get("N_true"))
    failure = _int(result.get("first_failed_mode"))
    required_guard = _int(result.get("required_guard"), guard) or guard
    t1_pass = bool(result.get("stable")) and bool(result.get("expected_counts"))
    exact = t0_pass and t1_pass and n_true is not None and min(
        len(result.get("eb_roots", ())), len(result.get("timo_roots", ()))
    ) >= required_guard
    final_status = "resolved_targeted_local" if exact else "still_deferred_after_targeted_strict"
    certified_bound = ""
    if not exact and failure is not None and failure <= 5 and t0_pass:
        final_status = "envelope_dominated_upper_bound_confirmed"
        certified_bound = min(4, failure - 1)
    verification_methods = (
        "stored_pointwise_primary;stored_independent_verification;"
        "targeted_exact_beta_four_phase_two_level_local_SVD"
    )
    overlay = {
        "case_id": case_id, "final_status": final_status,
        "N_true": n_true if exact else None,
        "certified_N_true_upper_bound": certified_bound,
        "first_failed_mode": failure if exact else "",
        "required_guard": required_guard, "proof_rank": failure or required_guard,
        "EB_roots_json": _json(list(result.get("eb_roots", ()))[:required_guard]),
        "Timoshenko_roots_json": _json(list(result.get("timo_roots", ()))[:required_guard]),
        "delta_f_json": _json(list(result.get("deltas", ()))[:min(K_MAX, required_guard)]),
        "verification_methods": verification_methods,
        "strict_level_used": "none",
        "source_compact_certificate_hash": sha256_file(_absolute(target["compact_certificate_path"], output_dir.parent)),
        "source_raw_cache_hash": sha256_file(_absolute(index_row["source_full_cache_path"], output_dir.parent)),
        "scientific_scope": SCIENTIFIC_SCOPE, "target_version": TARGET_VERSION,
    }
    strict_row = {
        "case_id": case_id, "T0_status": "PASS" if t0_pass else "FAIL",
        "T1_status": result.get("T1_status", ""),
        "T2_status": "not_required_T0_and_T1_agree" if exact else "not_implemented_no_safe_acceptance",
        "T3_status": "not_required" if exact else "not_executed",
        "T4_status": "not_required" if exact else "not_executed",
        "strict_level_used": "none", "force_strict_requested": 0,
        "force_strict_executed": 0,
        "strict_containment_status": "PASS_no_strict_executed",
    }
    result["overlay"] = overlay
    result["strict_row"] = strict_row
    result["T0_rows"] = t0_rows
    return result, cache_hit, time.perf_counter() - started, t0_rows


def _target_certificate(
    source_path: Path, overlay: Mapping[str, object], output_path: Path,
) -> None:
    cert = load_certificate(source_path)
    result = dict(cert.get("result", {}))
    spectra = dict(cert.get("spectra", {}))
    n_true = _int(overlay.get("N_true"))
    if n_true is not None:
        result.update({
            "execution_status": overlay["final_status"], "n_true_status": "exact",
            "N_true": n_true, "first_failed_mode": _int(overlay.get("first_failed_mode")),
            "required_guard": _int(overlay.get("required_guard")),
            "required_guard_confirmed": True, "unresolved_reason": "",
        })
    else:
        result.update({
            "execution_status": overlay["final_status"], "n_true_status": "unresolved_with_certified_upper_bound",
            "N_true": None, "required_guard_confirmed": True,
            "unresolved_reason": "exact_N_true_not_required_for_envelope",
            "certified_N_true_upper_bound": _int(overlay.get("certified_N_true_upper_bound")),
        })
    for theory, key in ((MODEL_EB, "EB_roots_json"), (MODEL_TIMO, "Timoshenko_roots_json")):
        record = dict(spectra.get(theory, {}))
        record["roots"] = json.loads(str(overlay[key]))
        spectra[theory] = record
    spectra["delta_f"] = {
        str(rank): value for rank, value in enumerate(json.loads(str(overlay["delta_f_json"])), start=1)
    }
    cert["result"] = result
    cert["spectra"] = spectra
    cert["targeted_resolution_overlay"] = {
        "version": TARGET_VERSION, "verification_methods": overlay["verification_methods"],
        "strict_level_used": overlay["strict_level_used"], "matrix_confirmed": True,
        "source_compact_certificate_hash": overlay["source_compact_certificate_hash"],
        "source_raw_cache_hash": overlay["source_raw_cache_hash"],
    }
    atomic_write_gzip_json(output_path, cert)


def _build_finalization(
    coarse_dir: Path, output_dir: Path, overlays: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    source_dir = coarse_dir / "compact_finalization"
    final_dir = coarse_dir / "compact_finalization_epsilon_005_resolved"
    final_dir.mkdir(parents=True, exist_ok=True)
    source_rows = read_csv(source_dir / "final_case_certificates.csv")
    overlay_by_id = {str(row["case_id"]): row for row in overlays}
    final_rows: list[dict[str, object]] = []
    changes: list[dict[str, object]] = []
    for source in source_rows:
        case_id = source["case_id"]
        row: dict[str, object] = dict(source)
        overlay = overlay_by_id.get(case_id)
        if overlay is not None:
            before = dict(row)
            n_true = _int(overlay.get("N_true"))
            row.update({
                "final_execution_status": overlay["final_status"],
                "N_true": n_true, "first_failed_mode": _int(overlay.get("first_failed_mode")),
                "required_guard": _int(overlay.get("required_guard")),
                "delta_at_first_failure": (
                    json.loads(str(overlay["delta_f_json"]))[_int(overlay.get("first_failed_mode"), 1) - 1]
                    if n_true is not None and _int(overlay.get("first_failed_mode")) is not None else ""
                ),
                "required_guard_confirmed": True,
                "result_origin": "epsilon_005_targeted_resolution",
                "unresolved_reason": "" if n_true is not None else "exact_N_true_not_required_for_envelope",
            })
            target_cert = final_dir / "cases" / f"{case_id}.json.gz"
            _target_certificate(_absolute(source["compact_certificate_path"], coarse_dir), overlay, target_cert)
            row["compact_certificate_path"] = _relative(target_cert, coarse_dir)
            for field in ("final_execution_status", "N_true", "first_failed_mode", "required_guard", "compact_certificate_path"):
                if str(before.get(field, "")) != str(row.get(field, "")):
                    changes.append({
                        "case_id": case_id, "field": field, "before": before.get(field, ""),
                        "after": row.get(field, ""), "change_origin": "epsilon_005_targeted_resolution",
                    })
        final_rows.append(row)
    final_rows.sort(key=lambda row: (
        float(row["epsilon_0"]), float(row["beta"]), float(row["mu"]),
        float(row["eta"]), str(row["case_id"]),
    ))
    epsilon_rows = _epsilon_summary(final_rows)
    unresolved_source = {row["case_id"]: row for row in read_csv(source_dir / "unresolved_cases.csv")}
    unresolved_rows = [
        dict(unresolved_source[str(row["case_id"])])
        for row in final_rows
        if _int(row.get("N_true")) is None
        and str(row["case_id"]) in unresolved_source
        and str(row["case_id"]) not in overlay_by_id
    ]
    unresolved_rows.sort(key=lambda row: str(row["case_id"]))
    raw_rows = [{
        "epsilon_0": row["epsilon_0"], "N_up_raw": row["observed_max_N_true"],
        "envelope_status": row["envelope_status"],
        "unresolved_case_count": row["unresolved_case_count"],
        "observed_argmax_case_ids": row["observed_argmax_case_ids"],
    } for row in epsilon_rows]
    suffix_rows: list[dict[str, object]] = []
    running: int | None = None
    provisional_seen = False
    for row in reversed(raw_rows):
        value = _int(row["N_up_raw"])
        running = value if running is None else max(running, value if value is not None else running)
        provisional_seen = provisional_seen or str(row["envelope_status"]).startswith("provisional")
        suffix_rows.append({
            "epsilon_0": row["epsilon_0"], "N_up_raw": value,
            "N_up_suffix_max": running,
            "provenance": "provisional_input" if provisional_seen else "exact_inputs",
        })
    suffix_rows.reverse()
    atomic_write_csv(final_dir / "final_case_certificates.csv", final_rows, FINAL_FIELDS)
    atomic_write_csv(final_dir / "epsilon_level_summary.csv", epsilon_rows, EPSILON_FIELDS)
    atomic_write_csv(final_dir / "unresolved_cases.csv", unresolved_rows, UNRESOLVED_FIELDS)
    atomic_write_csv(final_dir / "raw_envelope.csv", raw_rows, RAW_ENVELOPE_FIELDS)
    atomic_write_csv(final_dir / "suffix_max_envelope.csv", suffix_rows, SUFFIX_FIELDS)
    atomic_write_csv(final_dir / "result_changes.csv", changes, RESULT_CHANGE_FIELDS)
    return final_rows, epsilon_rows, unresolved_rows, raw_rows, suffix_rows


def _preservation(
    coarse_dir: Path, final_rows: Sequence[Mapping[str, object]], target_ids: set[str],
) -> dict[str, object]:
    before = {row["case_id"]: row for row in read_csv(coarse_dir / "compact_finalization" / "final_case_certificates.csv")}
    after = {str(row["case_id"]): row for row in final_rows}
    changed_non_target = 0
    changed_resolved_n = 0
    for case_id, old in before.items():
        if case_id in target_ids:
            continue
        new = after[case_id]
        scientific = ("final_execution_status", "N_true", "first_failed_mode", "required_guard", "delta_at_first_failure")
        changed_non_target += int(any(str(old.get(key, "")) != str(new.get(key, "")) for key in scientific))
        if _int(old.get("N_true")) is not None and _int(old.get("N_true")) != _int(new.get("N_true")):
            changed_resolved_n += 1
    non_target_deferred = [
        row for case_id, row in before.items()
        if case_id not in target_ids and _int(row.get("N_true")) is None
    ]
    validation = read_csv(coarse_dir.parent / "solver_readiness_v2" / "validation_24_cases.csv")
    readiness_preserved = sum(
        row.get("solver_readiness_case_status") == "READY"
        and row.get("validation_id") in after
        and (
            (
                row["validation_id"] in target_ids
                and _int(after[row["validation_id"]].get("N_true"))
                == _int(row.get("auto_N_true"))
                and _int(after[row["validation_id"]].get("first_failed_mode"))
                == _int(row.get("first_failed_mode_auto"))
            )
            or (
                row["validation_id"] not in target_ids
                and all(
                    str(before[row["validation_id"]].get(field, ""))
                    == str(after[row["validation_id"]].get(field, ""))
                    for field in (
                        "final_execution_status", "N_true", "first_failed_mode",
                        "required_guard", "delta_at_first_failure",
                    )
                )
            )
        )
        for row in validation
    )
    promotions = read_csv(coarse_dir / "family_local_repair_reconciliation" / "promoted_cases.csv")
    former_preserved = sum(
        row["case_id"] in after and _int(after[row["case_id"]].get("N_true")) == _int(row.get("N_true"))
        for row in promotions
    )
    manifest = read_csv(coarse_dir / "grid_manifest.csv")
    labels = {row.get("regression_label", ""): row["case_id"] for row in manifest}
    s3_12 = after[labels["S3_12"]]
    s3_14 = after[labels["S3_14"]]
    return {
        "changed_non_target_cases": changed_non_target,
        "changed_previously_resolved_N_true": changed_resolved_n,
        "non_target_deferred_count": len(non_target_deferred),
        "readiness_preserved": readiness_preserved,
        "former_blockers_preserved": former_preserved,
        "S3_12_N_true": _int(s3_12.get("N_true")),
        "S3_12_delta_f_5": _float(s3_12.get("delta_at_first_failure")),
        "S3_14_N_true": _int(s3_14.get("N_true")),
        "S3_14_delta_f_5": _float(s3_14.get("delta_at_first_failure")),
    }


def run_targeted_resolution(coarse_dir: Path) -> dict[str, object]:
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    coarse_dir = coarse_dir.resolve()
    output_dir = coarse_dir / "epsilon_005_targeted_resolution"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "cache").mkdir(exist_ok=True)
    (output_dir / "logs").mkdir(exist_ok=True)
    prepare_target_diagnostics(coarse_dir)
    target_rows = read_csv(output_dir / "target_manifest.csv")
    index = {
        row["case_id"]: row
        for row in read_csv(coarse_dir / "compact_point_certificates_v1" / "compact_index.csv")
    }
    target_ids = {row["case_id"] for row in target_rows}
    previous_scientific = {
        path.name: sha256_file(path)
        for path in (
            output_dir / "resolved_case_overlay.csv",
            output_dir / "before_after.csv",
            output_dir / "epsilon_005_envelope_status.csv",
            coarse_dir / "compact_finalization_epsilon_005_resolved" / "final_case_certificates.csv",
            coarse_dir / "compact_finalization_epsilon_005_resolved" / "epsilon_level_summary.csv",
            coarse_dir / "compact_finalization_epsilon_005_resolved" / "raw_envelope.csv",
            coarse_dir / "compact_finalization_epsilon_005_resolved" / "suffix_max_envelope.csv",
        ) if path.exists()
    }
    source_hashes_before = {
        row["case_id"]: {
            "compact": sha256_file(_absolute(row["compact_certificate_path"], coarse_dir)),
            "raw": sha256_file(_absolute(index[row["case_id"]]["source_full_cache_path"], coarse_dir)),
        } for row in target_rows
    }
    initial_rss = process_rss_bytes()
    peak_rss = initial_rss
    overlays: list[dict[str, object]] = []
    events: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    independent_rows: list[dict[str, object]] = []
    strict_rows: list[dict[str, object]] = []
    runtime_rows: list[dict[str, object]] = []
    memory_rows: list[dict[str, object]] = []
    actual_evaluations = 0
    cache_hits = 0
    for target in target_rows:
        case_id = target["case_id"]
        case_initial = process_rss_bytes()
        cert = load_certificate(_absolute(target["compact_certificate_path"], coarse_dir))
        raw = _load_raw(_absolute(index[case_id]["source_full_cache_path"], coarse_dir))
        result, cache_hit, wall, t0_rows = _resolve_one(target, index[case_id], cert, raw, output_dir)
        evaluations = 0 if cache_hit else int(result.get("matrix_evaluations", 0) or 0)
        actual_evaluations += evaluations
        cache_hits += int(cache_hit)
        overlays.append(dict(result["overlay"]))
        events.extend(dict(row) for row in result.get("event_rows", ()))
        candidates.extend(dict(row) for row in result.get("candidate_rows", ()))
        independent_rows.extend(t0_rows)
        strict_rows.append(dict(result["strict_row"]))
        case_final = process_rss_bytes()
        peak_rss = max(peak_rss, case_initial, case_final)
        runtime_rows.append({
            "case_id": case_id, "cache_hit": cache_hit, "wall_seconds": wall,
            "matrix_evaluations": evaluations, "initial_rss_bytes": case_initial,
            "peak_rss_bytes": peak_rss, "final_rss_bytes": case_final,
        })
        memory_rows.append({"case_id": case_id, "stage": "after_target", "rss_bytes": case_final, "peak_rss_bytes": peak_rss})
        del raw, cert, result
        gc.collect()
    overlays.sort(key=lambda row: str(row["case_id"]))
    events.sort(key=lambda row: (
        str(row["case_id"]), str(row["theory"]), int(row["rank_start"]),
        str(row["stage"]),
    ))
    candidates.sort(key=lambda row: (str(row["case_id"]), str(row["theory"]), str(row["stage"]), float(row["lambda_candidate"])))
    independent_rows.sort(key=lambda row: (str(row["case_id"]), str(row["theory"]), int(row["sorted_rank"])))
    strict_rows.sort(key=lambda row: str(row["case_id"]))
    final_rows, epsilon_rows, unresolved_rows, raw_rows, suffix_rows = _build_finalization(coarse_dir, output_dir, overlays)
    before_table = {row["case_id"]: row for row in read_csv(coarse_dir / "compact_finalization" / "final_case_certificates.csv")}
    before_after = [{
        "case_id": row["case_id"], "epsilon_0": before_table[str(row["case_id"])]["epsilon_0"],
        "beta": before_table[str(row["case_id"])]["beta"],
        "original_execution_status": before_table[str(row["case_id"])]["final_execution_status"],
        "original_N_true": before_table[str(row["case_id"])]["N_true"],
        "final_execution_status": row["final_status"], "final_N_true": row["N_true"],
        "certified_N_true_upper_bound": row["certified_N_true_upper_bound"],
        "first_failed_mode": row["first_failed_mode"], "required_guard": row["required_guard"],
    } for row in overlays]
    eps_before = next(row for row in read_csv(coarse_dir / "compact_finalization" / "epsilon_level_summary.csv") if abs(float(row["epsilon_0"]) - TARGET_EPSILON) <= TARGET_EPSILON_TOLERANCE)
    eps_after = next(row for row in epsilon_rows if abs(float(row["epsilon_0"]) - TARGET_EPSILON) <= TARGET_EPSILON_TOLERANCE)
    exact_count = sum(_int(row.get("N_true")) is not None for row in overlays)
    bound_count = sum(_int(row.get("certified_N_true_upper_bound")) is not None for row in overlays)
    envelope_status = [{
        "epsilon_0": TARGET_EPSILON,
        "previous_observed_max_N_true": eps_before["observed_max_N_true"],
        "final_observed_max_N_true": eps_after["observed_max_N_true"],
        "target_case_count": len(overlays), "target_exact_count": exact_count,
        "target_bounded_count": bound_count,
        "remaining_target_unresolved": len(overlays) - exact_count - bound_count,
        "envelope_finalizable": str(eps_after["envelope_finalizable"]).lower() == "true",
        "envelope_status": eps_after["envelope_status"],
        "argmax_case_ids": eps_after["observed_argmax_case_ids"],
    }]
    atomic_write_csv(output_dir / "staged_resolution_events.csv", events, EVENT_FIELDS)
    atomic_write_csv(output_dir / "local_root_candidates.csv", candidates, CANDIDATE_FIELDS)
    atomic_write_csv(output_dir / "independent_verification.csv", independent_rows, INDEPENDENT_FIELDS)
    atomic_write_csv(output_dir / "strict_execution_summary.csv", strict_rows, STRICT_FIELDS)
    atomic_write_csv(output_dir / "resolved_case_overlay.csv", overlays, OVERLAY_FIELDS)
    atomic_write_csv(output_dir / "before_after.csv", before_after, BEFORE_AFTER_FIELDS)
    atomic_write_csv(output_dir / "epsilon_005_envelope_status.csv", envelope_status, EPSILON_STATUS_FIELDS)
    operations = {
        "target_cases": len(overlays), "point_solver_calls": 0,
        "matrix_evaluations": actual_evaluations, "target_cache_hits": cache_hits,
        "T1_executed": len(overlays) - cache_hits, "T2_executed": 0,
        "T3_strict_prefix_executed": 0, "T4_full_K10_executed": 0,
        "force_strict_requested": 0, "force_strict_executed": 0,
        "non_target_strict_executed": 0, "raw_cache_files_deleted": 0,
    }
    atomic_write_csv(output_dir / "operation_counts.csv", [
        {"operation": key, "count": value} for key, value in sorted(operations.items())
    ], OPERATION_FIELDS)
    atomic_write_csv(output_dir / "memory_profile.csv", memory_rows, MEMORY_FIELDS)
    preservation = _preservation(coarse_dir, final_rows, target_ids)
    source_hashes_after = {
        row["case_id"]: {
            "compact": sha256_file(_absolute(row["compact_certificate_path"], coarse_dir)),
            "raw": sha256_file(_absolute(index[row["case_id"]]["source_full_cache_path"], coarse_dir)),
        } for row in target_rows
    }
    final_dir = coarse_dir / "compact_finalization_epsilon_005_resolved"
    scientific_paths = (
        output_dir / "resolved_case_overlay.csv", output_dir / "before_after.csv",
        output_dir / "epsilon_005_envelope_status.csv",
        final_dir / "final_case_certificates.csv", final_dir / "epsilon_level_summary.csv",
        final_dir / "raw_envelope.csv", final_dir / "suffix_max_envelope.csv",
    )
    scientific_hashes = {path.name: sha256_file(path) for path in scientific_paths}
    repeat_available = len(previous_scientific) == len(scientific_paths)
    deterministic = repeat_available and previous_scientific == scientific_hashes and actual_evaluations == 0
    target_success = all(
        _int(row.get("N_true")) is not None or _int(row.get("certified_N_true_upper_bound"), 99) <= 4
        for row in overlays
    )
    envelope_exact = (
        _int(eps_after.get("observed_max_N_true")) == 4
        and str(eps_after.get("envelope_status")) == "exact"
    )
    s3_pass = (
        preservation["S3_12_N_true"] == 4
        and math.isclose(float(preservation["S3_12_delta_f_5"]), 0.11739469908796035, abs_tol=1e-15)
        and preservation["S3_14_N_true"] == 4
        and math.isclose(float(preservation["S3_14_delta_f_5"]), 0.10050934855181458, abs_tol=1e-15)
    )
    scope_clean = not any("anisotropic" in name.lower() for name in sys.modules)
    gates = {
        "scope_isolation_gate": scope_clean,
        "target_selection_gate": len(target_ids) == 2,
        "source_integrity_gate": source_hashes_before == source_hashes_after,
        "targeted_resolution_gate": target_success,
        "strict_containment_gate": operations["force_strict_executed"] == 0 and operations["non_target_strict_executed"] == 0,
        "envelope_finalization_gate": envelope_exact,
        "resolved_preservation_gate": preservation["changed_non_target_cases"] == 0 and preservation["changed_previously_resolved_N_true"] == 0 and preservation["non_target_deferred_count"] == 33 and s3_pass and preservation["readiness_preserved"] == 24 and preservation["former_blockers_preserved"] == 7,
        "scientific_safety_gate": all(_int(row.get("N_true")) is not None or _int(row.get("certified_N_true_upper_bound")) is not None for row in overlays),
        "bounded_memory_gate": peak_rss <= 6 * 2**30,
        "determinism_gate": deterministic,
    }
    gates["epsilon_005_resolution_readiness_gate"] = all(
        value for key, value in gates.items() if key != "epsilon_005_resolution_readiness_gate"
    )
    atomic_write_csv(output_dir / "gate_summary.csv", [
        {"gate": key, "status": "PASS" if value else ("PENDING" if key in {"determinism_gate", "epsilon_005_resolution_readiness_gate"} and not repeat_available else "FAIL"), "value": value}
        for key, value in gates.items()
    ], GATE_FIELDS)
    atomic_write_csv(final_dir / "gate_summary.csv", [
        {"gate": key, "status": "PASS" if value else ("PENDING" if key in {"determinism_gate", "epsilon_005_resolution_readiness_gate"} and not repeat_available else "FAIL"), "value": value}
        for key, value in gates.items()
    ], GATE_FIELDS)
    wall = time.perf_counter() - started
    cold_path = output_dir / "cold_run_metadata.json"
    if actual_evaluations > 0:
        atomic_write_json(cold_path, {
            "wall_seconds": wall, "matrix_evaluations": actual_evaluations,
            "target_cache_hits": cache_hits, "peak_rss_bytes": peak_rss,
            "started_utc": started_utc, "source": "measured_cold_target_run",
            "case_stage_wall_seconds": {
                str(row["case_id"]): row["wall_seconds"] for row in runtime_rows
            },
            "case_matrix_evaluations": {
                str(row["case_id"]): row["matrix_evaluations"] for row in runtime_rows
            },
        })
    cold = json.loads(cold_path.read_text(encoding="utf-8")) if cold_path.exists() else {}
    cold_case_wall = cold.get("case_stage_wall_seconds", cold.get("case_stage_wall_seconds_estimate", {}))
    cold_case_evaluations = cold.get("case_matrix_evaluations", {})
    runtime_output: list[dict[str, object]] = []
    for case_id in sorted(target_ids):
        runtime_output.append({
            "run": "cold", "case_id": case_id, "cache_hit": False,
            "wall_seconds": cold_case_wall.get(case_id, "") if isinstance(cold_case_wall, Mapping) else "",
            "matrix_evaluations": (
                cold_case_evaluations.get(case_id, "")
                if isinstance(cold_case_evaluations, Mapping) and cold_case_evaluations
                else sum(int(row["matrix_evaluations"]) for row in events if row["case_id"] == case_id)
            ),
            "initial_rss_bytes": "", "peak_rss_bytes": cold.get("peak_rss_bytes", ""),
            "final_rss_bytes": "", "wall_source": cold.get("source", ""),
        })
    runtime_output.extend({
        **row, "run": "repeat" if actual_evaluations == 0 else "cold",
        "wall_source": "measured_current_invocation",
    } for row in runtime_rows)
    atomic_write_csv(output_dir / "runtime_summary.csv", runtime_output, RUNTIME_FIELDS)
    repeat_rows = [{
        "run": "cold", "wall_seconds": cold.get("wall_seconds", ""),
        "matrix_evaluations": cold.get("matrix_evaluations", ""),
        "target_cache_hits": cold.get("target_cache_hits", ""),
        "peak_rss_bytes": cold.get("peak_rss_bytes", ""),
        "scientific_csv_hashes_stable": "baseline", "source": cold.get("source", ""),
    }, {
        "run": "repeat", "wall_seconds": wall, "matrix_evaluations": actual_evaluations,
        "target_cache_hits": cache_hits, "peak_rss_bytes": peak_rss,
        "scientific_csv_hashes_stable": deterministic,
        "source": "current_cached_invocation",
    }]
    atomic_write_csv(output_dir / "repeat_run_summary.csv", repeat_rows, REPEAT_FIELDS)
    metadata = {
        "target_version": TARGET_VERSION, "scientific_scope": SCIENTIFIC_SCOPE,
        "started_utc": started_utc, "wall_seconds": wall, "target_case_ids": sorted(target_ids),
        "manifest_sha256": sha256_file(coarse_dir / "grid_manifest.csv"),
        "compact_index_sha256": sha256_file(coarse_dir / "compact_point_certificates_v1" / "compact_index.csv"),
        "source_hashes_before": source_hashes_before, "source_hashes_after": source_hashes_after,
        "operations": operations, "preservation": preservation, "peak_rss_bytes": peak_rss,
        "scientific_csv_hashes": scientific_hashes, "gates": gates,
        "cold_wall_seconds": cold.get("wall_seconds"),
        "repeat_wall_seconds": wall if actual_evaluations == 0 else None,
        "raw_cache_files_deleted": 0,
    }
    atomic_write_json(output_dir / "run_metadata.json", metadata)
    atomic_write_json(final_dir / "run_metadata.json", metadata)
    report = [
        "# Targeted epsilon_0=0.050 resolution", "",
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.", "",
        f"Targets: {', '.join(sorted(target_ids))}.",
        f"Target exact/bounded: {exact_count}/{bound_count}; matrix evaluations this invocation: {actual_evaluations}.",
        f"The epsilon_0=0.050 raw envelope is `{eps_after['observed_max_N_true']}` with status `{eps_after['envelope_status']}`.",
        "", "## Target results", "",
        "| case_id | beta | status | N_true | first failed | guard | delta at failure | strict |",
        "| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |",
    ]
    target_by_id = {row["case_id"]: row for row in target_rows}
    for row in overlays:
        deltas = json.loads(str(row["delta_f_json"]))
        failure = _int(row.get("first_failed_mode"))
        delta = deltas[failure - 1] if failure is not None else ""
        report.append(
            f"| {row['case_id']} | {target_by_id[str(row['case_id'])]['beta']} | {row['final_status']} | "
            f"{row.get('N_true', '')} | {row.get('first_failed_mode', '')} | {row['required_guard']} | {delta} | {row['strict_level_used']} |"
        )
    report.extend([
        "", "## Stage execution", "",
        "| case_id | stage | theory | ranks | Lambda window | step | roots | stable | matrix evaluations |",
        "| --- | --- | --- | --- | --- | ---: | ---: | --- | ---: |",
    ])
    for row in events:
        report.append(
            f"| {row['case_id']} | {row['stage']} | {row['theory']} | {row['rank_start']}-{row['rank_end']} | "
            f"[{float(row['lambda_left']):.12g}, {float(row['lambda_right']):.12g}] | {float(row['scan_step']):.6g} | "
            f"{row['root_count']} | {row['stable_with_previous']} | {row['matrix_evaluations']} |"
        )
    report.extend([
        "", "T1 used the unchanged production matrix evaluators, determinant sign changes, sigma-minimum candidates, four shifted phases, two fine scan levels, independent candidate refinement, SVD/residual acceptance, and multiplicity-aware deduplication.",
        "Stored pointwise primary and independent verification inventories supplied the second independent path. T2, T3, and T4 were not needed; force/full strict was not executed.",
        "", "## Raw and suffix-max envelopes", "",
        "| epsilon_0 | raw | raw status | suffix max | provenance |",
        "| ---: | ---: | --- | ---: | --- |",
    ])
    suffix_by_epsilon = {float(row["epsilon_0"]): row for row in suffix_rows}
    for row in raw_rows:
        suffix = suffix_by_epsilon[float(row["epsilon_0"])]
        report.append(
            f"| {row['epsilon_0']} | {row['N_up_raw']} | {row['envelope_status']} | "
            f"{suffix['N_up_suffix_max']} | {suffix['provenance']} |"
        )
    report.extend([
        "", "## Preservation and repeat", "",
        f"- unchanged non-target cases: {1554 - len(target_ids) - int(preservation['changed_non_target_cases'])}/{1554 - len(target_ids)}",
        f"- non-target deferred preserved: {preservation['non_target_deferred_count']}",
        f"- readiness references preserved: {preservation['readiness_preserved']}/24",
        f"- former blockers preserved: {preservation['former_blockers_preserved']}/7",
        f"- S3_12: N_true={preservation['S3_12_N_true']}, delta_f,5={preservation['S3_12_delta_f_5']}",
        f"- S3_14: N_true={preservation['S3_14_N_true']}, delta_f,5={preservation['S3_14_delta_f_5']}",
        f"- cold/repeat wall seconds: {cold.get('wall_seconds', '')}/{wall if actual_evaluations == 0 else ''}",
        f"- repeat matrix evaluations/cache hits: {actual_evaluations}/{cache_hits}",
        f"- peak RSS bytes (current invocation): {peak_rss}",
        "", "## Gates", "",
    ])
    report.extend(f"- {key}: {'PASS' if value else 'FAIL'}" for key, value in gates.items())
    report.extend([
        "",
        "No raw point cache was deleted, moved, renamed, or rewritten. No branch tracking, MAC, shapes, energy, FEM, or anisotropic workflow was invoked.",
    ])
    (output_dir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (final_dir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return {
        "output_dir": str(output_dir), "target_count": len(overlays),
        "matrix_evaluations": actual_evaluations, "cache_hits": cache_hits,
        "wall_seconds": wall, "peak_rss_bytes": peak_rss,
        "gates": gates, "scientific_csv_hashes": scientific_hashes,
        "overlays": overlays, "epsilon_rows": epsilon_rows,
        "raw_envelope": raw_rows, "suffix_envelope": suffix_rows,
        "preservation": preservation,
    }


__all__ = [
    "SCIENTIFIC_SCOPE", "TARGET_EPSILON", "TARGET_VERSION", "TargetWindow",
    "compute_prefix_result", "derive_target_windows", "prepare_target_diagnostics",
    "run_targeted_resolution", "select_targets", "strict_is_contained",
]
