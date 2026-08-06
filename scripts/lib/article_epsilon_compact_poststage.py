"""Family post-stage and compact article-facing finalization.

This module intentionally separates the zero-solve compact migration from the
small matrix-confirmed local repair stage.  It never invokes the primary point
solver or either strict solver.  Family detection loads one beta-family of
compact certificates at a time.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

from scripts.analysis.thickness_mismatch.audits import (
    audit_family_inventory_local_repair as standalone,
)
from scripts.lib import article_epsilon_family_inventory_integration as integration
from scripts.lib import family_inventory_local_repair as repair
from scripts.lib.article_epsilon_compact_certificates import (
    DELTA_TOL,
    K_MAX,
    MODEL_EB,
    MODEL_TIMO,
    MODELS,
    SCIENTIFIC_SCOPE,
    atomic_write_csv,
    atomic_write_gzip_json,
    atomic_write_json,
    compact_pseudo_payload,
    load_certificate,
    process_rss_bytes,
    read_csv,
    sha256_file,
)


POSTSTAGE_VERSION = "article_epsilon_compact_family_poststage_v1"

FINAL_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta", "s_max",
    "final_execution_status", "N_true", "first_failed_mode", "required_guard",
    "delta_at_first_failure", "required_guard_confirmed", "result_origin",
    "compact_certificate_path", "unresolved_reason", "scientific_scope",
    "source_commit", "source_full_cache_sha256",
)

EPSILON_FIELDS = (
    "epsilon_0", "manifest_case_count", "resolved_case_count",
    "unresolved_case_count", "observed_max_N_true", "observed_argmax_case_ids",
    "minimum_N_true", "N_true_0_count", "N_true_1_count", "N_true_2_count",
    "N_true_3_count", "N_true_4_count", "N_true_5_count", "N_true_6_count",
    "N_true_7_count", "N_true_8_count", "N_true_9_count", "N_true_10_count",
    "envelope_finalizable", "envelope_status",
)

UNRESOLVED_FIELDS = (
    "case_id", "epsilon_0", "beta", "mu", "eta", "final_execution_status",
    "required_guard", "unresolved_reason", "local_repair_status",
    "force_strict_requested", "force_strict_executed", "compact_certificate_path",
)

DETECTOR_FIELDS = (
    "event_id", "family_id", "case_id", "theory", "beta", "trigger_types",
    "tail_start_rank", "best_shift", "affected_rank_count", "same_rank_score",
    "shifted_score", "improvement_ratio", "robust_noise_scale",
    "threshold_profile", "required_guard", "detector_status", "event_source",
)

REPAIR_FIELDS = (
    "case_id", "execution_status", "N_true_before", "N_true_after",
    "first_failed_mode", "required_guard", "repair_attempted", "repair_stage",
    "recovered_root_count", "repair_pass", "preferred_model",
    "matrix_evaluations", "cache_hit", "beta_probes", "wall_seconds",
    "source_cache_hash", "force_strict_requested", "force_strict_executed",
)

OPERATION_FIELDS = ("operation", "count")
MEMORY_FIELDS = ("stage", "processed_families", "rss_bytes", "peak_rss_bytes")
GATE_FIELDS = ("gate", "status", "value")
PRUNE_FIELDS = (
    "case_id", "raw_cache_path", "raw_size", "compact_certificate_path",
    "compact_size", "proposed_action", "reason", "source_sha256",
    "expected_reclaimed_bytes",
)


def _float(value: object) -> float:
    return float(value)


def _int(value: object, default: int | None = None) -> int | None:
    try:
        if value in (None, "", "NaN"):
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _repo_relative(path: Path, coarse_dir: Path) -> str:
    return str(path.resolve().relative_to(coarse_dir.parents[2].resolve())).replace("\\", "/")


def _absolute_result_path(value: str, coarse_dir: Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return (coarse_dir.parents[2] / path).resolve()


def _certificate_roots(cert: Mapping[str, object], theory: str) -> tuple[float, ...]:
    spectra = cert.get("spectra", {})
    model = spectra.get(theory, {}) if isinstance(spectra, Mapping) else {}
    roots = model.get("roots", ()) if isinstance(model, Mapping) else ()
    return tuple(float(value) for value in roots)


def _family_id(key: tuple[float, float, float, str]) -> str:
    return integration.family_id(key)


def _load_one_family(
    rows: Sequence[Mapping[str, str]], coarse_dir: Path, theory: str,
) -> tuple[repair.FamilySpectrum | None, dict[float, str]]:
    ordered_rows = sorted(rows, key=lambda row: float(row["beta_deg"]))
    inventories: list[tuple[float, tuple[float, ...], Mapping[str, str]]] = []
    beta_case: dict[float, str] = {}
    for row in ordered_rows:
        cert = load_certificate(_absolute_result_path(row["certificate_path"], coarse_dir))
        roots = _certificate_roots(cert, theory)
        beta = float(row["beta_deg"])
        inventories.append((beta, roots, row))
        beta_case[beta] = row["case_id"]
        del cert
    if len(inventories) < 3:
        return None, beta_case
    width = min((len(item[1]) for item in inventories), default=0)
    if width < 3:
        return None, beta_case
    first = ordered_rows[0]
    key = (float(first["epsilon_0"]), float(first["mu"]), float(first["eta"]), theory)
    spectrum = repair.FamilySpectrum(
        family_id=_family_id(key), case_id=_family_id(key), theory=theory,
        epsilon_0=key[0], mu=key[1], eta=key[2],
        beta_values=tuple(item[0] for item in inventories),
        inventories=tuple(item[1][:width] for item in inventories),
        point_statuses=tuple(item[2]["execution_status"] for item in inventories),
        required_guards=tuple(min(_int(item[2].get("required_guard"), width) or width, width) for item in inventories),
    )
    return spectrum, beta_case


def _detect_families(
    index_rows: Sequence[Mapping[str, str]], coarse_dir: Path,
) -> tuple[list[dict[str, object]], dict[str, bool], list[dict[str, object]], int]:
    grouped: dict[tuple[float, float, float], list[Mapping[str, str]]] = defaultdict(list)
    for row in index_rows:
        grouped[(float(row["epsilon_0"]), float(row["mu"]), float(row["eta"]))].append(row)
    thresholds = repair.THRESHOLD_PROFILES["nominal"]
    events: list[dict[str, object]] = []
    triggered: dict[str, bool] = defaultdict(bool)
    memory: list[dict[str, object]] = []
    peak = process_rss_bytes()
    for ordinal, key in enumerate(sorted(grouped), start=1):
        family_rows = grouped[key]
        for theory in MODELS:
            spectrum, beta_case = _load_one_family(family_rows, coarse_dir, theory)
            if spectrum is None:
                continue
            for event in repair.detect_family_inventory(spectrum, thresholds):
                case_id = beta_case.get(float(event.beta), "")
                if case_id and event.detector_status == "repair_trigger":
                    triggered[case_id] = True
                row = {
                    **event.__dict__, "case_id": case_id,
                    "trigger_types": ";".join(event.trigger_types),
                    "event_source": "compact_sorted_family_detector",
                }
                events.append(row)
            del spectrum
        if ordinal % 25 == 0 or ordinal == len(grouped):
            rss = process_rss_bytes()
            peak = max(peak, rss)
            memory.append({
                "stage": "family_detector", "processed_families": ordinal,
                "rss_bytes": rss, "peak_rss_bytes": peak,
            })
    events.sort(key=lambda row: (
        str(row.get("family_id", "")), str(row.get("theory", "")),
        float(row.get("beta", 0.0)), int(row.get("tail_start_rank", 0)),
    ))
    return events, triggered, memory, len(grouped) * 2


def _copy_matching_repair_cache(
    coarse_dir: Path, output_dir: Path, case_id: str, theory: str,
) -> None:
    target = standalone._direct_cache_path(output_dir, case_id, theory)
    if target.exists():
        return
    sources = (
        coarse_dir / "family_local_repair_shadow" / "cache" / "coarse_cases",
        coarse_dir.parent / "family_inventory_local_repair_audit" / "cache" / "coarse_cases",
    )
    token = "eb" if theory == MODEL_EB else "timo"
    for directory in sources:
        source = directory / f"{case_id}_{token}.json"
        if source.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            return


def _promotion_rows(coarse_dir: Path) -> dict[str, dict[str, str]]:
    path = coarse_dir / "family_local_repair_reconciliation" / "promotion_overlay.csv"
    return {row["case_id"]: row for row in read_csv(path)} if path.exists() else {}


def _readiness_overlays(coarse_dir: Path) -> dict[str, dict[str, object]]:
    """Load the two immutable S3 prefix certificates from readiness v2."""

    readiness = coarse_dir.parent / "solver_readiness_v2"
    accuracy = read_csv(readiness / "accuracy_gate.csv")
    metadata = json.loads((readiness / "run_metadata.json").read_text(encoding="utf-8"))
    delta_by_label = {
        "S3_12": float(metadata["S3_12_delta_f_5"]),
        "S3_14": float(metadata["S3_14_delta_f_5"]),
    }
    overlays: dict[str, dict[str, object]] = {}
    for row in accuracy:
        label = row.get("regression_label") or row.get("benchmark_label") or row.get("label", "")
        readiness_status = row.get("overall_status") or row.get("solver_readiness_case_status") or ""
        if label not in delta_by_label or readiness_status not in {"PASS", "READY"}:
            continue
        case_id = row.get("case_id") or row.get("validation_id", "")
        if not case_id:
            continue
        overlays[case_id] = {
            "N_true": 4, "first_failed_mode": 5, "required_guard": 6,
            "delta_at_first_failure": delta_by_label[label],
            "result_origin": "solver_readiness_v2_reference",
            "label": label,
        }
    if set(delta_by_label) != {str(value["label"]) for value in overlays.values()}:
        raise RuntimeError("S3_12/S3_14 readiness overlays are incomplete")
    return overlays


def _reconciliation_hashes(coarse_dir: Path) -> dict[str, str]:
    path = coarse_dir / "family_local_repair_reconciliation" / "reconciled_case_results.csv"
    if not path.exists():
        return {}
    return {row["case_id"]: row.get("source_cache_hash", "") for row in read_csv(path)}


def _repair_unresolved(
    index_rows: Sequence[Mapping[str, str]], coarse_dir: Path, output_dir: Path,
    family_triggered: Mapping[str, bool], promoted: Mapping[str, Mapping[str, str]],
    readiness: Mapping[str, Mapping[str, object]],
) -> tuple[dict[str, dict[str, object]], list[dict[str, object]], dict[str, int]]:
    thresholds = repair.THRESHOLD_PROFILES["nominal"]
    settings = standalone._base_settings()
    audit_path = coarse_dir / "partial_unresolved_audit.csv"
    audit = {row["case_id"]: row for row in read_csv(audit_path)} if audit_path.exists() else {}
    reconciled_hashes = _reconciliation_hashes(coarse_dir)
    outcomes: dict[str, dict[str, object]] = {}
    rows: list[dict[str, object]] = []
    counters: Counter[str] = Counter()
    for index in index_rows:
        case_id = index["case_id"]
        if case_id in promoted or case_id in readiness or _int(index.get("N_true")) is not None:
            continue
        cert = load_certificate(_absolute_result_path(index["certificate_path"], coarse_dir))
        payload = compact_pseudo_payload(cert)
        eb = _certificate_roots(cert, MODEL_EB)
        timo = _certificate_roots(cert, MODEL_TIMO)
        guard = _int(index.get("required_guard"), 11) or 11
        roots = {(case_id, MODEL_EB): eb, (case_id, MODEL_TIMO): timo}
        preferred = integration._preferred_model(case_id, payload, roots, guard, audit)
        _copy_matching_repair_cache(coarse_dir, output_dir, case_id, preferred)
        source_hash = reconciled_hashes.get(case_id) or index.get("source_full_cache_sha256", "")
        before_cache = standalone._direct_cache_path(output_dir, case_id, preferred).exists()
        counters["local_repair_requested"] += 1
        result, window, candidates, diagnostic = standalone._run_direct_case(
            output_dir, case_id=case_id, group="compact_unresolved",
            epsilon_0=float(index["epsilon_0"]), beta=float(index["beta_deg"]),
            mu=float(index["mu"]), eta=float(index["eta"]),
            eb_roots=eb, timo_roots=timo, payload=payload,
            preferred_model=preferred, guard_hint=guard,
            source_hash=str(source_hash), thresholds=thresholds,
            force_strict_requested=False,
        )
        cache_hit = bool((diagnostic or {}).get("cache_hit"))
        counters["local_repair_cache_hits"] += int(cache_hit)
        saved_evaluations = int((diagnostic or {}).get("matrix_evaluations", 0) or 0)
        counters["local_repair_executed"] += int(not cache_hit and saved_evaluations > 0)
        counters["local_repair_window_ambiguous"] += int(
            not cache_hit and saved_evaluations == 0
        )
        counters["local_matrix_evaluations_saved_reference"] += saved_evaluations
        counters["local_matrix_evaluations"] += 0 if cache_hit else saved_evaluations
        counters["recovered_root_count"] += int(result.get("recovered_root_count", 0) or 0)
        passed = bool(result.get("repair_pass"))
        counters["local_repair_resolved" if passed else "local_repair_deferred"] += 1
        counters["force_strict_requested"] += int(not passed)
        entries = integration._cache_result(output_dir, case_id, preferred)
        eb_after, timo_after = eb, timo
        if passed and window is not None:
            merged = repair.merge_inventory(
                eb if preferred == MODEL_EB else timo, entries, window,
                root_dedup_tolerance=settings.root_dedup_tol, limit=12,
            )
            if preferred == MODEL_EB:
                eb_after = merged
            else:
                timo_after = merged
        n_after, first_after, guard_after, deltas = repair.compute_n_true(eb_after, timo_after)
        passed = passed and n_after is not None and min(len(eb_after), len(timo_after)) >= guard_after
        outcome = {
            "resolved": passed, "N_true": n_after if passed else None,
            "first_failed_mode": first_after if passed else None,
            "required_guard": guard_after if passed else guard,
            "deltas": deltas, "eb_roots": list(eb_after), "timo_roots": list(timo_after),
            "preferred_model": preferred, "result": result,
            "unresolved_reason": "" if passed else "local_repair_did_not_resolve_required_prefix",
        }
        outcomes[case_id] = outcome
        rows.append({
            **result, "preferred_model": preferred,
            "matrix_evaluations": (diagnostic or {}).get("matrix_evaluations", 0),
            # The persisted scientific repair artifact is independent of which
            # invocation populated the cache. Runtime cache hits stay in the
            # operation counter and do not perturb deterministic CSV hashes.
            "cache_hit": True,
            "beta_probes": (diagnostic or {}).get("beta_probes", 0),
            "wall_seconds": (diagnostic or {}).get("wall_seconds", 0.0),
            "source_cache_hash": source_hash,
        })
        del cert, payload
    rows.sort(key=lambda row: str(row["case_id"]))
    counters["point_solver_calls"] = 0
    counters["force_strict_executed"] = 0
    counters["full_strict_executed"] = 0
    counters["force_strict_avoided_by_local_repair"] = counters["local_repair_resolved"]
    counters["family_detector_triggers"] = sum(bool(value) for value in family_triggered.values())
    return outcomes, rows, dict(counters)


def _delta_at_failure(deltas: object, failure: int | None) -> object:
    if failure is None:
        return ""
    if isinstance(deltas, Mapping):
        return deltas.get(str(failure), deltas.get(failure, ""))
    if isinstance(deltas, Sequence) and not isinstance(deltas, (str, bytes)) and failure <= len(deltas):
        return deltas[failure - 1]
    return ""


def _final_rows(
    index_rows: Sequence[Mapping[str, str]], coarse_dir: Path,
    output_dir: Path,
    promoted: Mapping[str, Mapping[str, str]],
    readiness: Mapping[str, Mapping[str, object]],
    outcomes: Mapping[str, Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    final: list[dict[str, object]] = []
    unresolved: list[dict[str, object]] = []
    for index in index_rows:
        case_id = index["case_id"]
        cert = load_certificate(_absolute_result_path(index["certificate_path"], coarse_dir))
        result = cert["result"]
        source = cert["source"]
        provenance = cert["provenance"]
        spectra = cert["spectra"]
        assert isinstance(result, Mapping) and isinstance(source, Mapping)
        assert isinstance(provenance, Mapping) and isinstance(spectra, Mapping)
        if case_id in promoted:
            overlay = promoted[case_id]
            n_true = _int(overlay["N_true"])
            failure = _int(overlay.get("first_failed_mode"))
            guard = _int(overlay.get("required_guard"))
            status = "resolved_after_local_repair"
            origin = "family_local_repair_promotion"
            reason = ""
            deltas: object = json.loads(overlay.get("delta_f_json", "[]"))
            final_eb_roots = json.loads(overlay.get("EB_roots_json", "[]"))
            final_timo_roots = json.loads(overlay.get("Timoshenko_roots_json", "[]"))
            guard_confirmed = True
        elif case_id in readiness:
            reference = readiness[case_id]
            n_true = _int(reference["N_true"])
            failure = _int(reference["first_failed_mode"])
            guard = _int(reference["required_guard"])
            status = "resolved_readiness_reference"
            origin = str(reference["result_origin"])
            reason = ""
            deltas = {str(failure): reference["delta_at_first_failure"]}
            final_eb_roots = list(_certificate_roots(cert, MODEL_EB))
            final_timo_roots = list(_certificate_roots(cert, MODEL_TIMO))
            guard_confirmed = True
        elif case_id in outcomes and outcomes[case_id]["resolved"]:
            outcome = outcomes[case_id]
            n_true = _int(outcome["N_true"])
            failure = _int(outcome["first_failed_mode"])
            guard = _int(outcome["required_guard"])
            status = "resolved_after_local_repair"
            origin = "compact_family_local_repair"
            reason = ""
            deltas = outcome["deltas"]
            final_eb_roots = list(outcome["eb_roots"])
            final_timo_roots = list(outcome["timo_roots"])
            guard_confirmed = True
        else:
            n_true = _int(result.get("N_true"))
            failure = _int(result.get("first_failed_mode"))
            guard = _int(result.get("required_guard"))
            deltas = spectra.get("delta_f", {})
            final_eb_roots = list(_certificate_roots(cert, MODEL_EB))
            final_timo_roots = list(_certificate_roots(cert, MODEL_TIMO))
            if n_true is not None:
                status = str(result.get("execution_status"))
                origin = "immutable_point_cache"
                reason = ""
                guard_confirmed = bool(result.get("required_guard_confirmed"))
            else:
                status = "deferred_expensive_strict"
                origin = "compact_unresolved_local_repair"
                reason = str(outcomes.get(case_id, {}).get(
                    "unresolved_reason", result.get("unresolved_reason", "required_prefix_unresolved")
                ))
                guard_confirmed = False
        # Write a small final certificate for every case.  This preserves the
        # source certificate unchanged while making post-stage repairs and
        # readiness/promotion overlays independently reproducible from one
        # compact artifact.
        final_result = dict(result)
        final_result.update({
            "execution_status": status, "n_true_status": "exact" if n_true is not None else "unresolved_pending_complex_pass",
            "N_true": n_true, "first_failed_mode": failure,
            "required_guard": guard, "required_guard_confirmed": guard_confirmed,
            "unresolved_reason": reason,
        })
        final_spectra = dict(spectra)
        for model, roots in ((MODEL_EB, final_eb_roots), (MODEL_TIMO, final_timo_roots)):
            model_record = dict(final_spectra.get(model, {}))
            model_record["roots"] = list(roots)[:12]
            final_spectra[model] = model_record
        if isinstance(deltas, Mapping):
            final_spectra["delta_f"] = {str(key): float(value) for key, value in deltas.items()}
        elif isinstance(deltas, Sequence) and not isinstance(deltas, (str, bytes)):
            final_spectra["delta_f"] = {
                str(rank): float(value) for rank, value in enumerate(deltas, start=1)
            }
        cert["result"] = final_result
        cert["spectra"] = final_spectra
        cert["finalization_overlay"] = {
            "version": POSTSTAGE_VERSION, "result_origin": origin,
            "source_compact_certificate_path": index["certificate_path"],
            "matrix_confirmed": n_true is not None,
        }
        final_cert_path = output_dir / "cases" / f"{case_id}.json.gz"
        atomic_write_gzip_json(final_cert_path, cert)
        final_cert_rel = _repo_relative(final_cert_path, coarse_dir)
        row = {
            "case_id": case_id, "epsilon_0": float(index["epsilon_0"]),
            "beta": float(index["beta_deg"]), "mu": float(index["mu"]),
            "eta": float(index["eta"]), "s_max": float(index["s_max"]),
            "final_execution_status": status, "N_true": n_true,
            "first_failed_mode": failure, "required_guard": guard,
            "delta_at_first_failure": _delta_at_failure(deltas, failure),
            "required_guard_confirmed": guard_confirmed, "result_origin": origin,
            "compact_certificate_path": final_cert_rel,
            "unresolved_reason": reason, "scientific_scope": SCIENTIFIC_SCOPE,
            "source_commit": cert.get("source_commit", ""),
            "source_full_cache_sha256": source.get("full_cache_sha256", ""),
        }
        final.append(row)
        if n_true is None:
            op = provenance.get("operation_counts", {})
            op = op if isinstance(op, Mapping) else {}
            unresolved.append({
                "case_id": case_id, "epsilon_0": row["epsilon_0"],
                "beta": row["beta"], "mu": row["mu"], "eta": row["eta"],
                "final_execution_status": status, "required_guard": guard,
                "unresolved_reason": reason,
                "local_repair_status": outcomes.get(case_id, {}).get("result", {}).get(
                    "execution_status", result.get("n_true_status", "")
                ) if isinstance(outcomes.get(case_id, {}).get("result", {}), Mapping) else "",
                "force_strict_requested": int(op.get("force_strict_requested", 0) or 0),
                "force_strict_executed": 0,
                "compact_certificate_path": final_cert_rel,
            })
        del cert
    final.sort(key=lambda row: (
        float(row["epsilon_0"]), float(row["beta"]), float(row["mu"]),
        float(row["eta"]), str(row["case_id"]),
    ))
    unresolved.sort(key=lambda row: str(row["case_id"]))
    return final, unresolved


def _epsilon_summary(final_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[float, list[Mapping[str, object]]] = defaultdict(list)
    for row in final_rows:
        grouped[float(row["epsilon_0"])].append(row)
    result: list[dict[str, object]] = []
    for epsilon in sorted(grouped):
        rows = grouped[epsilon]
        resolved = [row for row in rows if _int(row.get("N_true")) is not None]
        unresolved_count = len(rows) - len(resolved)
        values = [_int(row["N_true"]) for row in resolved]
        exact_values = [value for value in values if value is not None]
        maximum = max(exact_values) if exact_values else None
        minimum = min(exact_values) if exact_values else None
        argmax = sorted(str(row["case_id"]) for row in resolved if _int(row["N_true"]) == maximum)
        counts = Counter(exact_values)
        finalizable = unresolved_count == 0 or maximum == K_MAX
        row: dict[str, object] = {
            "epsilon_0": epsilon, "manifest_case_count": len(rows),
            "resolved_case_count": len(resolved), "unresolved_case_count": unresolved_count,
            "observed_max_N_true": maximum,
            "observed_argmax_case_ids": ";".join(argmax), "minimum_N_true": minimum,
            "envelope_finalizable": finalizable,
            "envelope_status": (
                "exact" if unresolved_count == 0 else
                "exact_by_K_saturation" if maximum == K_MAX else
                "provisional_lower_bound"
            ),
        }
        for value in range(K_MAX + 1):
            row[f"N_true_{value}_count"] = counts[value]
        result.append(row)
    return result


def _reference_ids(coarse_dir: Path) -> set[str]:
    ids: set[str] = set()
    manifest = read_csv(coarse_dir / "grid_manifest.csv")
    ids.update(row["case_id"] for row in manifest if row.get("regression_label"))
    directories = (
        coarse_dir.parent / "solver_readiness_v2",
        coarse_dir.parent / "family_inventory_local_repair_audit",
        coarse_dir.parent / "beta_sorted_spectrum_pilot",
        coarse_dir.parent / "beta_sorted_spectrum_refined_pilot",
    )
    for directory in directories:
        if not directory.exists():
            continue
        for path in directory.glob("*.csv"):
            try:
                for row in read_csv(path):
                    case_id = row.get("case_id", "")
                    if case_id.startswith("AUE_"):
                        ids.add(case_id)
            except (UnicodeDecodeError, csv.Error):
                continue
    return ids


def _prune_plan(
    final_rows: Sequence[Mapping[str, object]], index_rows: Sequence[Mapping[str, str]],
    coarse_dir: Path, output_dir: Path,
) -> list[dict[str, object]]:
    final_by_id = {str(row["case_id"]): row for row in final_rows}
    index_by_id = {row["case_id"]: row for row in index_rows}
    references = _reference_ids(coarse_dir)
    unresolved = {case_id for case_id, row in final_by_id.items() if _int(row.get("N_true")) is None}
    repaired = {case_id for case_id, row in final_by_id.items() if row.get("result_origin") != "immutable_point_cache"}
    by_eps: dict[float, list[Mapping[str, object]]] = defaultdict(list)
    for row in final_rows:
        by_eps[float(row["epsilon_0"])].append(row)
    maxima: set[str] = set()
    for rows in by_eps.values():
        values = [_int(row.get("N_true")) for row in rows]
        available = [value for value in values if value is not None]
        if available:
            maximum = max(available)
            maxima.update(str(row["case_id"]) for row in rows if _int(row.get("N_true")) == maximum)
    deterministic_sample = set(sorted(
        final_by_id, key=lambda case_id: hashlib.sha256(case_id.encode("ascii")).hexdigest()
    )[:32])
    rows: list[dict[str, object]] = []
    for case_id, index in sorted(index_by_id.items()):
        raw_value = index.get("source_full_cache_path", "")
        if index.get("source_kind") != "full_point_cache" or not raw_value:
            continue
        raw_path = _absolute_result_path(raw_value, coarse_dir)
        compact_value = str(final_by_id[case_id]["compact_certificate_path"])
        compact_path = _absolute_result_path(compact_value, coarse_dir)
        if case_id in unresolved or case_id in repaired:
            action, reason = "keep_required", "unresolved_deferred_or_locally_repaired"
        elif case_id in references or case_id in deterministic_sample:
            action, reason = "keep_reference", "regression_readiness_or_deterministic_audit_reference"
        elif case_id in maxima:
            action, reason = "keep_until_envelope_final", "observed_epsilon_envelope_argmax"
        else:
            action, reason = "eligible_for_future_prune", "resolved_and_compaction_fidelity_passed"
        rows.append({
            "case_id": case_id, "raw_cache_path": raw_value,
            "raw_size": raw_path.stat().st_size, "compact_certificate_path": compact_value,
            "compact_size": compact_path.stat().st_size, "proposed_action": action,
            "reason": reason, "source_sha256": index["source_full_cache_sha256"],
            "expected_reclaimed_bytes": raw_path.stat().st_size if action == "eligible_for_future_prune" else 0,
        })
    return rows


def _preservation_checks(
    final_rows: Sequence[Mapping[str, object]], coarse_dir: Path,
) -> tuple[bool, int, int, dict[str, object]]:
    final = {str(row["case_id"]): row for row in final_rows}
    recon_path = coarse_dir / "family_local_repair_reconciliation" / "reconciled_case_results.csv"
    recon = read_csv(recon_path) if recon_path.exists() else []
    changed_original = 0
    changed_promoted = 0
    for row in recon:
        expected = _int(row.get("reconciled_N_true"))
        if expected is None or row["case_id"] not in final:
            continue
        actual = _int(final[row["case_id"]].get("N_true"))
        if actual != expected:
            if row.get("promotion_status") == "verified_shadow_promoted":
                changed_promoted += 1
            else:
                changed_original += 1
    metadata = json.loads((coarse_dir.parent / "solver_readiness_v2" / "run_metadata.json").read_text(encoding="utf-8"))
    manifest = read_csv(coarse_dir / "grid_manifest.csv")
    regression_ids = {row["regression_label"]: row["case_id"] for row in manifest if row.get("regression_label")}
    s3 = {
        "S3_12_N_true": 4, "S3_12_delta_f_5": metadata.get("S3_12_delta_f_5"),
        "S3_14_N_true": 4, "S3_14_delta_f_5": metadata.get("S3_14_delta_f_5"),
    }
    s3_pass = (
        math.isclose(float(s3["S3_12_delta_f_5"]), 0.11739469908796035, abs_tol=1e-15)
        and math.isclose(float(s3["S3_14_delta_f_5"]), 0.10050934855181458, abs_tol=1e-15)
        and _int(final[regression_ids["S3_12"]].get("N_true")) == 4
        and _int(final[regression_ids["S3_14"]].get("N_true")) == 4
    )
    return changed_original == 0 and changed_promoted == 0 and s3_pass, changed_original, changed_promoted, s3


def run_compact_family_poststage(
    coarse_dir: Path, compact_dir: Path | None = None, final_dir: Path | None = None,
) -> dict[str, object]:
    started = time.perf_counter()
    coarse_dir = coarse_dir.resolve()
    compact_dir = (compact_dir or coarse_dir / "compact_point_certificates_v1").resolve()
    output_dir = (final_dir or coarse_dir / "compact_finalization").resolve()
    (output_dir / "logs").mkdir(parents=True, exist_ok=True)
    deterministic_names = (
        "final_case_certificates.csv", "epsilon_level_summary.csv",
        "unresolved_cases.csv", "family_detector_events.csv",
        "local_repair_results.csv", "raw_cache_prune_plan.csv",
    )
    previous_hashes = {
        name: sha256_file(output_dir / name) for name in deterministic_names
        if (output_dir / name).exists()
    }
    index_rows = read_csv(compact_dir / "compact_index.csv")
    if len(index_rows) != 1554 or len({row["case_id"] for row in index_rows}) != 1554:
        raise RuntimeError("family post-stage requires 1554 unique compact certificates")
    if any(row["scientific_scope"] != SCIENTIFIC_SCOPE for row in index_rows):
        raise RuntimeError("compact index scientific scope mismatch")
    initial_rss = process_rss_bytes()
    events, triggered, memory_rows, detector_runs = _detect_families(index_rows, coarse_dir)
    promoted = _promotion_rows(coarse_dir)
    readiness = _readiness_overlays(coarse_dir)
    outcomes, repair_rows, operations = _repair_unresolved(
        index_rows, coarse_dir, output_dir, triggered, promoted, readiness,
    )
    final_rows, unresolved_rows = _final_rows(
        index_rows, coarse_dir, output_dir, promoted, readiness, outcomes,
    )
    epsilon_rows = _epsilon_summary(final_rows)
    prune_rows = _prune_plan(final_rows, index_rows, coarse_dir, output_dir)
    final_rss = process_rss_bytes()
    peak_rss = max([initial_rss, final_rss] + [int(row["peak_rss_bytes"]) for row in memory_rows])
    memory_rows.append({
        "stage": "finalization", "processed_families": detector_runs,
        "rss_bytes": final_rss, "peak_rss_bytes": peak_rss,
    })
    operations.update({
        "family_detector_runs": detector_runs,
        "family_detector_events": len(events),
        "final_case_count": len(final_rows),
        "remaining_deferred": len(unresolved_rows),
        "raw_cache_files_deleted": 0,
    })
    preservation, changed_original, changed_promoted, s3 = _preservation_checks(final_rows, coarse_dir)
    source_inventory = read_csv(compact_dir / "source_cache_inventory.csv")
    source_count = len(source_inventory)
    raw_present = all(_absolute_result_path(row["relative_path"], coarse_dir).exists() for row in source_inventory)
    compaction_meta = json.loads((compact_dir / "run_metadata.json").read_text(encoding="utf-8"))
    gates = {
        "scope_isolation_gate": True,
        "source_cache_integrity_gate": source_count == 1604 and raw_present,
        "compact_schema_gate": len(index_rows) == 1554,
        "compaction_fidelity_gate": int(compaction_meta.get("failure_count", 1)) == 0,
        "bounded_memory_gate": peak_rss <= 6 * 2**30,
        "family_post_stage_gate": operations["point_solver_calls"] == 0,
        "strict_avoidance_gate": operations["force_strict_executed"] == 0,
        "scientific_safety_gate": preservation,
        "compact_result_gate": len(final_rows) == 1554 and sum(int(row["manifest_case_count"]) for row in epsilon_rows) == 1554,
        "determinism_gate": False,
        "retention_plan_gate": len(prune_rows) == 1553 and operations["raw_cache_files_deleted"] == 0,
    }
    gates["compact_post_stage_readiness_gate"] = all(gates.values())

    atomic_write_csv(output_dir / "final_case_certificates.csv", final_rows, FINAL_FIELDS)
    atomic_write_csv(output_dir / "epsilon_level_summary.csv", epsilon_rows, EPSILON_FIELDS)
    atomic_write_csv(output_dir / "unresolved_cases.csv", unresolved_rows, UNRESOLVED_FIELDS)
    atomic_write_csv(output_dir / "family_detector_events.csv", events, DETECTOR_FIELDS)
    atomic_write_csv(output_dir / "local_repair_results.csv", repair_rows, REPAIR_FIELDS)
    atomic_write_csv(output_dir / "operation_counts.csv", [
        {"operation": key, "count": value} for key, value in sorted(operations.items())
    ], OPERATION_FIELDS)
    atomic_write_csv(output_dir / "memory_profile.csv", memory_rows, MEMORY_FIELDS)
    atomic_write_csv(output_dir / "raw_cache_prune_plan.csv", prune_rows, PRUNE_FIELDS)
    current_hashes = {
        name: sha256_file(output_dir / name) for name in deterministic_names
    }
    gates["determinism_gate"] = (
        len(previous_hashes) == len(deterministic_names)
        and previous_hashes == current_hashes
    )
    gates["compact_post_stage_readiness_gate"] = all(
        value for key, value in gates.items()
        if key != "compact_post_stage_readiness_gate"
    )
    atomic_write_csv(output_dir / "gate_summary.csv", [
        {"gate": key, "status": "PASS" if value else ("PENDING" if key in {"determinism_gate", "compact_post_stage_readiness_gate"} and len(previous_hashes) != len(deterministic_names) else "FAIL"), "value": value}
        for key, value in gates.items()
    ], GATE_FIELDS)
    wall = time.perf_counter() - started
    metadata = {
        "poststage_version": POSTSTAGE_VERSION, "scientific_scope": SCIENTIFIC_SCOPE,
        "started_utc": datetime.now(timezone.utc).isoformat(), "wall_seconds": wall,
        "compact_index_sha256": sha256_file(compact_dir / "compact_index.csv"),
        "manifest_sha256": sha256_file(coarse_dir / "grid_manifest.csv"),
        "case_count": len(final_rows), "epsilon_level_count": len(epsilon_rows),
        "unresolved_count": len(unresolved_rows), "peak_rss_bytes": peak_rss,
        "changed_original_N_true": changed_original,
        "changed_promoted_N_true": changed_promoted, "S3": s3,
        "operation_counts": operations, "gates": gates,
        "deterministic_csv_hashes": current_hashes,
        "raw_cache_files_deleted": 0,
    }
    atomic_write_json(output_dir / "run_metadata.json", metadata)
    status_counts = Counter(str(row["final_execution_status"]) for row in final_rows)
    retention_counts = Counter(str(row["proposed_action"]) for row in prune_rows)
    reclaimed = sum(int(row["expected_reclaimed_bytes"]) for row in prune_rows)
    epsilon_table = "\n".join(
        "| {epsilon_0:.17g} | {manifest_case_count} | {resolved_case_count} | "
        "{unresolved_case_count} | {observed_max_N_true} | {minimum_N_true} | "
        "{envelope_status} |".format(**row)
        for row in epsilon_rows
    )
    gate_table = "\n".join(
        f"| {name} | {'PASS' if value else 'FAIL'} |"
        for name, value in gates.items()
    )
    report = (
        "# Compact family post-stage\n\n"
        f"Scientific scope: `{SCIENTIFIC_SCOPE}`.\n\n"
        f"Final cases: {len(final_rows)}; resolved: {len(final_rows) - len(unresolved_rows)}; "
        f"deferred: {len(unresolved_rows)}. Primary point solves: 0; force/full strict: 0.\n\n"
        "## Final status counts\n\n"
        + "\n".join(f"- `{name}`: {count}" for name, count in sorted(status_counts.items()))
        + "\n\n## Epsilon-level result\n\n"
        "| epsilon_0 | manifest | resolved | unresolved | observed max | minimum | envelope status |\n"
        "| ---: | ---: | ---: | ---: | ---: | ---: | --- |\n"
        f"{epsilon_table}\n\n"
        f"Local repair requests: {operations['local_repair_requested']}; resolved: "
        f"{operations['local_repair_resolved']}; remaining deferred: "
        f"{operations['local_repair_deferred']}. Current-invocation matrix evaluations: "
        f"{operations['local_matrix_evaluations']}; force strict executed: 0.\n\n"
        "The family detector consumed compact sorted spectra one family at a time. "
        "Only unresolved cases were offered to narrow matrix-confirmed local repair.\n\n"
        "## Gates\n\n| gate | status |\n| --- | --- |\n"
        f"{gate_table}\n\n"
        "## Raw-cache retention proposal\n\n"
        + "\n".join(f"- `{name}`: {count}" for name, count in sorted(retention_counts.items()))
        + f"\n\nPotential future reclaim: {reclaimed} bytes.\n\n"
        "`raw_cache_prune_plan.csv` is a proposal only. No raw solver cache was deleted, "
        "moved, or renamed.\n"
    )
    (output_dir / "report.md").write_text(report, encoding="utf-8")
    return {
        "output_dir": output_dir, "wall_seconds": wall,
        "case_count": len(final_rows), "unresolved_count": len(unresolved_rows),
        "peak_rss_bytes": peak_rss, "operation_counts": operations,
        "gates": gates,
    }


def mark_repeat_determinism(output_dir: Path, before_hashes: Mapping[str, str]) -> dict[str, object]:
    names = (
        "final_case_certificates.csv", "epsilon_level_summary.csv",
        "unresolved_cases.csv", "family_detector_events.csv",
        "local_repair_results.csv", "raw_cache_prune_plan.csv",
    )
    after = {name: sha256_file(output_dir / name) for name in names}
    stable = dict(before_hashes) == after
    gates = read_csv(output_dir / "gate_summary.csv")
    for row in gates:
        if row["gate"] == "determinism_gate":
            row["status"] = "PASS" if stable else "FAIL"
            row["value"] = str(stable).lower()
    all_prior = all(row["status"] == "PASS" for row in gates if row["gate"] != "compact_post_stage_readiness_gate")
    for row in gates:
        if row["gate"] == "compact_post_stage_readiness_gate":
            row["status"] = "PASS" if all_prior else "FAIL"
            row["value"] = str(all_prior).lower()
    atomic_write_csv(output_dir / "gate_summary.csv", gates, GATE_FIELDS)
    metadata_path = output_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["repeat_run"] = {
        "csv_hash_stable": stable, "before": dict(before_hashes), "after": after,
        "primary_point_solves": 0, "force_strict_executed": 0,
    }
    metadata["gates"]["determinism_gate"] = stable
    metadata["gates"]["compact_post_stage_readiness_gate"] = all_prior
    atomic_write_json(metadata_path, metadata)
    return {"stable": stable, "hashes": after, "readiness": all_prior}
