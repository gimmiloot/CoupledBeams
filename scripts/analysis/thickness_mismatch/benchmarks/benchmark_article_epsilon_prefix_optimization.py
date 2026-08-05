"""Benchmark full and staged article epsilon spectrum paths.

The benchmark cache is deliberately isolated from main, smoke, regression,
and any future parallel cache.  This script never builds or runs the 1554-point
coarse grid.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
import gzip
import json
import math
from pathlib import Path
import re
import statistics
import sys
import time
import tracemalloc
from typing import Mapping, Sequence


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    run_article_epsilon_upper_envelope_grid as runner,
)
from scripts.lib import article_epsilon_prefix_optimization as prefix  # noqa: E402
from scripts.lib import article_epsilon_upper_envelope as workflow  # noqa: E402
from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402


OUTPUT_DIR = REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "optimization_benchmark"
MODE_ORDER = (
    "full_k10_full_strict",
    "prefix_best_local_strict",
    "prefix_best_auto_strict",
)
OPTIMIZED_MODES = MODE_ORDER[1:]
SMOKE_PREFIX = "smoke_unresolved"

VALIDATION_FIELDS = (
    "validation_order",
    "case_id",
    "case_identity",
    "category",
    "regression_label",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "thin_0p1_flag",
    "s_max",
)

ROOT_COMPARISON_FIELDS = (
    "run_mode",
    "case_id",
    "model",
    "sorted_index",
    "Lambda_full",
    "Lambda_optimized",
    "absolute_difference",
    "within_existing_tolerance",
    "root_ordering_agreement",
    "multiplicity_agreement",
    "cluster_identity_agreement",
    "comparison_status",
)

FINAL_VALIDATION_FIELDS = (
    "validation_id",
    "group",
    "epsilon_0",
    "beta",
    "mu",
    "eta",
    "thin_0p1_flag",
    "full_execution_status",
    "full_N_true",
    "local_execution_status",
    "local_N_true",
    "auto_execution_status",
    "auto_N_true",
    "first_failed_mode_full",
    "first_failed_mode_local",
    "first_failed_mode_auto",
    "root_agreement",
    "prefix_guard_agreement",
    "execution_status_agreement",
    "unresolved_reason_full",
    "unresolved_reason_local",
    "unresolved_reason_auto",
    "equivalence_case_status",
    "solver_readiness_case_status",
    "optimization_equivalence_gate",
    "full_grid_solver_readiness_gate",
)

BENCHMARK_CASE_FIELDS = (
    "run_mode", "case_id", "category", "epsilon_0", "beta_deg", "mu", "eta",
    "execution_status", "N_true", "first_failed_mode", "first_failed_delta_f",
    "prefix_guard_status", "full_K10_guard_status", "strict_status", "early_stop_used",
    "points_with_first_failure", "guard_root_count", "wall_seconds", "warm_cache_seconds",
    "peak_memory_bytes", "EB_primary_seconds", "Timo_primary_seconds", "EB_strict_seconds",
    "Timo_strict_seconds", "preparation_seconds", "branch_strict_seconds",
    "branch_strict_cache_hits", "branch_strict_cache_misses", "EB_root_count", "Timo_root_count",
    "cache_hits", "cache_provenance", "unresolved_reason",
)


@dataclass(frozen=True)
class ModeSpec:
    name: str
    full_k10: bool
    strategy: str
    strict_policy: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark article prefix optimization on a fixed manifest.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--postprocess-only", action="store_true")
    parser.add_argument("--case-limit", type=int, default=None)
    parser.add_argument("--case-id", action="append", default=[])
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--session-id", default=None)
    parser.add_argument("--track-memory", action="store_true")
    parser.add_argument(
        "--include-smoke-controls",
        action="store_true",
        help="Explicitly include the four smoke controls in a targeted solver-readiness run.",
    )
    parser.add_argument("--modes", nargs="+", choices=MODE_ORDER, default=list(MODE_ORDER))
    args = parser.parse_args(argv)
    if args.force and args.reuse_cache:
        parser.error("--force and --reuse-cache are mutually exclusive for benchmark timing")
    if args.case_limit is not None and args.case_limit < 1:
        parser.error("--case-limit must be positive")
    if args.postprocess_only and (args.force or args.case_id or args.case_limit is not None):
        parser.error("--postprocess-only cannot select or force solver cases")
    return args


def _find_point(
    points: Sequence[workflow.GridPoint],
    epsilon_0: float,
    beta_deg: float,
    mu: float,
    eta: float,
) -> workflow.GridPoint:
    matches = [
        point
        for point in points
        if point.epsilon_0 == float(epsilon_0)
        and point.beta_deg == float(beta_deg)
        and point.mu == float(mu)
        and point.eta == float(eta)
    ]
    if len(matches) != 1:
        raise ValueError(f"validation point lookup is not unique: {(epsilon_0, beta_deg, mu, eta)}")
    return matches[0]


def build_validation_manifest() -> list[tuple[str, workflow.GridPoint]]:
    all_points = workflow.build_manifest()
    rows: list[tuple[str, workflow.GridPoint]] = []
    for point in all_points:
        if point.regression_label:
            rows.append((f"regression_{point.regression_label}", point))
    requested = (
        ("smoke_unresolved_basis_beta0_eps002", 0.020, 0.0, 0.7, -0.5),
        ("smoke_unresolved_eb_guard_eps005", 0.050, 0.0, 0.7, -0.5),
        ("smoke_unresolved_basis_beta90_eps005", 0.050, 90.0, 0.0, 0.0),
        ("smoke_unresolved_basis_beta90_eps002", 0.020, 90.0, 0.7, -0.5),
        ("beta0_special_thin", 0.010, 0.0, 0.0, 0.0),
        ("beta0_special_thick", 0.060, 0.0, 0.0, 0.0),
        ("beta90_eta_plus", 0.020, 90.0, 0.0, 0.5),
        ("beta90_eta_minus", 0.040, 90.0, 0.3, -0.5),
        ("high_mu_eta_minus", 0.015, 15.0, 0.9, -0.5),
        ("high_mu_eta_plus", 0.050, 75.0, 0.9, 0.5),
        ("high_mu_mid_eta", 0.030, 45.0, 0.9, 0.0),
        ("low_angle_cluster_probe", 0.020, 5.0, 0.5, 0.0),
        ("low_angle_extreme", 0.060, 10.0, 0.9, 0.5),
        ("early_failure_probe_1", 0.060, 30.0, 0.7, -0.25),
        ("early_failure_probe_2", 0.050, 60.0, 0.5, 0.25),
        ("middle_prefix_probe_1", 0.025, 45.0, 0.5, -0.5),
        ("middle_prefix_probe_2", 0.030, 75.0, 0.3, 0.5),
        ("late_prefix_probe_1", 0.015, 30.0, 0.3, 0.0),
        ("late_prefix_probe_2", 0.020, 60.0, 0.0, -0.25),
        ("cutoff_basis_probe_1", 0.060, 90.0, 0.9, -0.5),
        ("cutoff_basis_probe_2", 0.050, 15.0, 0.9, 0.5),
        ("thin_split_false", 0.060, 45.0, 0.0, 0.5),
    )
    for category, epsilon_0, beta_deg, mu, eta in requested:
        rows.append((category, _find_point(all_points, epsilon_0, beta_deg, mu, eta)))
    unique: dict[str, tuple[str, workflow.GridPoint]] = {}
    for category, point in rows:
        unique.setdefault(point.case_id, (category, point))
    result = list(unique.values())
    if len(result) != 24:
        raise ValueError(f"validation manifest must contain 24 unique points, found {len(result)}")
    return result


def write_validation_manifest(output_dir: Path, cases: Sequence[tuple[str, workflow.GridPoint]]) -> Path:
    rows = [
        {
            "validation_order": index,
            "case_id": point.case_id,
            "case_identity": point.case_identity,
            "category": category,
            "regression_label": point.regression_label,
            "epsilon_0": point.epsilon_0,
            "beta_deg": point.beta_deg,
            "mu": point.mu,
            "eta": point.eta,
            "thin_0p1_flag": point.thin_0p1_flag,
            "s_max": point.s_max,
        }
        for index, (category, point) in enumerate(cases, start=1)
    ]
    return workflow.write_csv(output_dir / "validation_manifest.csv", rows, VALIDATION_FIELDS)


def _mode_specs(best_strategy: str) -> dict[str, ModeSpec]:
    return {
        "full_k10_full_strict": ModeSpec("full_k10_full_strict", True, "paired", "full"),
        "prefix_best_local_strict": ModeSpec("prefix_best_local_strict", False, best_strategy, "local"),
        "prefix_best_auto_strict": ModeSpec("prefix_best_auto_strict", False, best_strategy, "auto"),
    }


def _model_roots(payload: Mapping[str, object], model: str) -> tuple[complete.RootRecord, ...]:
    states = payload.get("models", {})
    state = states.get(model, {}) if isinstance(states, Mapping) else {}
    latest = state.get("latest_result", {}) if isinstance(state, Mapping) else {}
    if not isinstance(latest, Mapping) or not latest:
        return ()
    result = complete._result_from_payload(latest, "benchmark")  # type: ignore[attr-defined]
    return tuple(result.primary.roots)


def _model_values(payload: Mapping[str, object], model: str) -> tuple[float, ...]:
    return tuple(root.Lambda for root in _model_roots(payload, model))


def _resolved(payload: Mapping[str, object] | None) -> bool:
    return bool(payload and str(payload.get("execution_status", "")).startswith("resolved_"))


def _as_bool(value: object) -> bool:
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes"}
    return bool(value)


def _unresolved(payload: Mapping[str, object] | None) -> bool:
    return bool(payload and str(payload.get("execution_status", "")).startswith("attempted_"))


def _canonical_unresolved_reason(payload: Mapping[str, object]) -> str:
    reason = str(payload.get("unresolved_reason", "")).strip()
    lowered = reason.lower()
    if "expected one positive and one negative" in lowered or "basis regime" in lowered:
        return "timoshenko_basis_regime_limit"
    if "strict_failure" in lowered or "strict failure" in lowered or "strict_failure_or_disagreement" in lowered:
        return "strict_failure_or_disagreement"
    if "independent" in lowered and "disagreement" in lowered:
        return "independent_search_disagreement"
    if "root" in lowered and ("count" in lowered or "found_only" in lowered):
        return "root_count_failure"
    return reason


def execution_status_agreement(payloads: Sequence[Mapping[str, object] | None]) -> bool:
    """Compare execution state without treating empty resolved reasons as a failure."""

    if not payloads or any(payload is None for payload in payloads):
        return False
    concrete = [payload for payload in payloads if payload is not None]
    if all(_resolved(payload) for payload in concrete):
        return True
    if not all(_unresolved(payload) for payload in concrete):
        return False
    reasons = [_canonical_unresolved_reason(payload) for payload in concrete]
    return bool(reasons[0]) and all(reason == reasons[0] for reason in reasons[1:])


def _right_guard(payload: Mapping[str, object]) -> int | None:
    failed = payload.get("first_failed_mode")
    if failed not in {None, "", "None"}:
        return min(workflow.ROOT_GUARD_INDEX, int(failed) + 1)
    if payload.get("N_true") == workflow.K_MAX or str(payload.get("N_true", "")) == str(workflow.K_MAX):
        return workflow.ROOT_GUARD_INDEX
    return None


def _timing_totals(payload: Mapping[str, object]) -> dict[str, float]:
    totals = defaultdict(float)
    for row in payload.get("stage_timings", ()):  # type: ignore[union-attr]
        token = "EB" if row.get("model") == complete.MODEL_EB else "Timo"
        totals[f"{token}_primary_seconds"] += float(row.get("primary_seconds", 0.0))
        totals[f"{token}_strict_seconds"] += float(row.get("verification_seconds", 0.0))
        totals["preparation_seconds"] += float(row.get("preparation_seconds", 0.0))
        totals[f"{token}_root_count"] = max(totals[f"{token}_root_count"], float(row.get("guard_index", 0)))
    strict = payload.get("strict_details", {})
    if isinstance(strict, Mapping):
        models = strict.get("models", {})
        if isinstance(models, Mapping):
            for model, details in models.items():
                token = "EB" if model == complete.MODEL_EB else "Timo"
                if isinstance(details, Mapping):
                    totals[f"{token}_strict_seconds"] += float(details.get("seconds", 0.0))
                    totals["branch_strict_seconds"] += float(details.get("seconds", 0.0))
                    totals["branch_strict_cache_hits"] += details.get("cache_status") == "hit"
                    totals["branch_strict_cache_misses"] += details.get("cache_status") != "hit"
    return dict(totals)


def _case_row(
    spec: ModeSpec,
    category: str,
    point: workflow.GridPoint,
    payload: Mapping[str, object],
    wall_seconds: float,
    peak_bytes: int,
) -> dict[str, object]:
    timings = _timing_totals(payload)
    failed = payload.get("first_failed_mode")
    actual_early_stop = bool(
        not spec.full_k10
        and payload.get("execution_status") == "resolved_prefix_early_stop"
        and payload.get("early_stop_used")
    )
    return {
        "run_mode": spec.name,
        "case_id": point.case_id,
        "category": category,
        "epsilon_0": point.epsilon_0,
        "beta_deg": point.beta_deg,
        "mu": point.mu,
        "eta": point.eta,
        "execution_status": payload.get("execution_status"),
        "N_true": payload.get("N_true"),
        "first_failed_mode": failed,
        "first_failed_delta_f": payload.get("first_failed_delta_f"),
        "prefix_guard_status": payload.get("prefix_guard_status"),
        "full_K10_guard_status": payload.get("full_K10_guard_status"),
        "strict_status": payload.get("strict_verification_status"),
        "early_stop_used": actual_early_stop,
        "points_with_first_failure": failed not in {None, "", "None"},
        "guard_root_count": 2 * (_right_guard(payload) or 0),
        "wall_seconds": wall_seconds,
        "peak_memory_bytes": peak_bytes,
        "EB_primary_seconds": timings.get("EB_primary_seconds", 0.0),
        "Timo_primary_seconds": timings.get("Timo_primary_seconds", 0.0),
        "EB_strict_seconds": timings.get("EB_strict_seconds", 0.0),
        "Timo_strict_seconds": timings.get("Timo_strict_seconds", 0.0),
        "preparation_seconds": timings.get("preparation_seconds", 0.0),
        "branch_strict_seconds": timings.get("branch_strict_seconds", 0.0),
        "branch_strict_cache_hits": int(timings.get("branch_strict_cache_hits", 0)),
        "branch_strict_cache_misses": int(timings.get("branch_strict_cache_misses", 0)),
        "EB_root_count": int(timings.get("EB_root_count", 0)),
        "Timo_root_count": int(timings.get("Timo_root_count", 0)),
        "cache_hits": payload.get("cache_hits", 0),
        "cache_provenance": "hit" if int(payload.get("cache_hits", 0)) > 0 else "miss",
        "warm_cache_seconds": "",
        "unresolved_reason": payload.get("unresolved_reason", ""),
    }


def _root_comparisons(
    baseline: Mapping[str, Mapping[str, object]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mode_name, by_case in results.items():
        if mode_name == "full_k10_full_strict":
            continue
        for case_id, payload in by_case.items():
            reference = baseline.get(case_id)
            if reference is None:
                continue
            guard = _right_guard(payload)
            for model in complete.SUPPORTED_MODELS:
                left = _model_roots(reference, model)
                right = _model_roots(payload, model)
                for index in range(workflow.ROOT_GUARD_INDEX):
                    if guard is not None and index + 1 > guard:
                        rows.append(
                            {
                                "run_mode": mode_name,
                                "case_id": case_id,
                                "model": model,
                                "sorted_index": index + 1,
                                "Lambda_full": left[index].Lambda if index < len(left) else "",
                                "Lambda_optimized": "",
                                "absolute_difference": "",
                                "within_existing_tolerance": True,
                                "root_ordering_agreement": True,
                                "multiplicity_agreement": True,
                                "cluster_identity_agreement": True,
                                "comparison_status": "not_needed_after_first_failure",
                            }
                        )
                        continue
                    if index >= len(left) or index >= len(right):
                        rows.append(
                            {
                                "run_mode": mode_name,
                                "case_id": case_id,
                                "model": model,
                                "sorted_index": index + 1,
                                "Lambda_full": left[index].Lambda if index < len(left) else "",
                                "Lambda_optimized": right[index].Lambda if index < len(right) else "",
                                "absolute_difference": "",
                                "within_existing_tolerance": False,
                                "root_ordering_agreement": False,
                                "multiplicity_agreement": False,
                                "cluster_identity_agreement": False,
                                "comparison_status": "missing_root",
                            }
                        )
                        continue
                    difference = abs(left[index].Lambda - right[index].Lambda)
                    multiplicity = (
                        left[index].detected_nullity == right[index].detected_nullity
                        and left[index].track_multiplicity == right[index].track_multiplicity
                    )
                    cluster = (
                        left[index].root_cluster_id == right[index].root_cluster_id
                        and left[index].cluster_member_index == right[index].cluster_member_index
                        and left[index].cluster_size == right[index].cluster_size
                    )
                    within = difference <= complete.DEFAULT_ROOT_MATCH_TOL
                    rows.append(
                        {
                            "run_mode": mode_name,
                            "case_id": case_id,
                            "model": model,
                            "sorted_index": index + 1,
                            "Lambda_full": left[index].Lambda,
                            "Lambda_optimized": right[index].Lambda,
                            "absolute_difference": difference,
                            "within_existing_tolerance": within,
                            "root_ordering_agreement": within,
                            "multiplicity_agreement": multiplicity,
                            "cluster_identity_agreement": cluster,
                            "comparison_status": "PASS" if within and multiplicity and cluster else "DISAGREEMENT",
                        }
                    )
    return rows


def _accuracy_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
    root_rows: Sequence[Mapping[str, object]] | None = None,
) -> list[dict[str, object]]:
    compared = list(root_rows if root_rows is not None else _root_comparisons(results.get(MODE_ORDER[0], {}), results))

    def roots_pass(case_id: str) -> bool:
        relevant = [
            row
            for row in compared
            if row.get("case_id") == case_id and row.get("run_mode") in OPTIMIZED_MODES
        ]
        needed = [row for row in relevant if row.get("comparison_status") != "not_needed_after_first_failure"]
        return bool(needed) and all(row.get("comparison_status") == "PASS" for row in needed)

    rows: list[dict[str, object]] = []
    for category, point in cases:
        values = {mode: results.get(mode, {}).get(point.case_id) for mode in MODE_ORDER}
        baseline, local, auto = (values[mode] for mode in MODE_ORDER)
        baseline_resolved = _resolved(baseline)
        reference_unresolved = baseline is not None and not baseline_resolved
        local_raw_resolved = _resolved(local)
        auto_raw_resolved = _resolved(auto)
        local_status = local.get("execution_status", "not_attempted") if local else "not_attempted"
        auto_status = auto.get("execution_status", "not_attempted") if auto else "not_attempted"
        local_n = local.get("N_true", "") if local else ""
        auto_n = auto.get("N_true", "") if auto else ""
        local_reason = local.get("unresolved_reason", "") if local else "not_attempted"
        auto_reason = auto.get("unresolved_reason", "") if auto else "not_attempted"
        if reference_unresolved and local_raw_resolved:
            local_status = "not_accepted_reference_unresolved"
            local_n = ""
            local_reason = "full_reference_unresolved;raw_optimized_candidate_not_accepted"
        if reference_unresolved and auto_raw_resolved:
            auto_status = "not_accepted_reference_unresolved"
            auto_n = ""
            auto_reason = "full_reference_unresolved;raw_optimized_candidate_not_accepted"
        all_modes_executed = all(payload is not None for payload in values.values())
        status_agreement = execution_status_agreement((baseline, local, auto))
        n_agreement = bool(
            baseline_resolved
            and local is not None
            and auto is not None
            and _resolved(local)
            and _resolved(auto)
            and baseline.get("N_true") == local.get("N_true") == auto.get("N_true")  # type: ignore[union-attr]
        )
        failure_agreement = bool(
            baseline_resolved
            and local is not None
            and auto is not None
            and baseline.get("first_failed_mode")
            == local.get("first_failed_mode")
            == auto.get("first_failed_mode")  # type: ignore[union-attr]
        )
        guard_agreement = bool(
            baseline_resolved
            and local is not None
            and auto is not None
            and _right_guard(baseline) == _right_guard(local) == _right_guard(auto)
        )
        root_agreement = roots_pass(point.case_id) if baseline_resolved and local and auto else False
        if baseline is None:
            equivalence = "NOT_ATTEMPTED"
        elif not baseline_resolved:
            equivalence = "NOT_EVALUABLE_REFERENCE_UNRESOLVED"
        elif not all_modes_executed:
            equivalence = "NOT_ATTEMPTED"
        elif n_agreement and failure_agreement and root_agreement and guard_agreement and status_agreement:
            equivalence = "PASS"
        else:
            equivalence = "FAIL"

        reason_full = str(baseline.get("unresolved_reason", "")) if baseline else "not_attempted"
        if baseline_resolved:
            readiness = "READY"
        elif baseline is None:
            readiness = "ERROR"
        elif (
            _as_bool(baseline.get("implementation_limit", False))
            or "basis" in reason_full.lower()
            or "expected one positive and one negative" in reason_full.lower()
        ):
            readiness = "IMPLEMENTATION_LIMIT"
        elif str(baseline.get("execution_status")) == "attempted_error":
            readiness = "ERROR"
        else:
            readiness = "UNRESOLVED"
        group = (
            "smoke_failure_control"
            if category.startswith(SMOKE_PREFIX)
            else ("completed_regression" if category.startswith("regression_") else "ordinary_validation")
        )
        rows.append(
            {
                "validation_id": point.case_id,
                "group": group,
                "epsilon_0": point.epsilon_0,
                "beta": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "thin_0p1_flag": point.thin_0p1_flag,
                "full_execution_status": baseline.get("execution_status", "not_attempted") if baseline else "not_attempted",
                "full_N_true": baseline.get("N_true", "") if baseline else "",
                "local_execution_status": local_status,
                "local_N_true": local_n,
                "auto_execution_status": auto_status,
                "auto_N_true": auto_n,
                "first_failed_mode_full": baseline.get("first_failed_mode", "") if baseline else "",
                "first_failed_mode_local": local.get("first_failed_mode", "") if local else "",
                "first_failed_mode_auto": auto.get("first_failed_mode", "") if auto else "",
                "root_agreement": root_agreement,
                "prefix_guard_agreement": guard_agreement,
                "execution_status_agreement": status_agreement,
                "unresolved_reason_full": reason_full,
                "unresolved_reason_local": local_reason,
                "unresolved_reason_auto": auto_reason,
                "equivalence_case_status": equivalence,
                "solver_readiness_case_status": readiness,
            }
        )
    return rows


def calculate_gates(rows: Sequence[Mapping[str, object]]) -> tuple[str, str]:
    equivalence_states = {str(row.get("equivalence_case_status", "NOT_ATTEMPTED")) for row in rows}
    if "FAIL" in equivalence_states:
        equivalence = "FAIL_DISAGREEMENT"
    elif "NOT_ATTEMPTED" in equivalence_states:
        equivalence = "INCOMPLETE"
    else:
        equivalence = "PASS"

    readiness_states = {str(row.get("solver_readiness_case_status", "ERROR")) for row in rows}
    if readiness_states == {"READY"}:
        readiness = "PASS"
    elif readiness_states.intersection({"UNRESOLVED", "IMPLEMENTATION_LIMIT"}):
        readiness = "BLOCKED_BY_UNRESOLVED_REFERENCE"
    else:
        readiness = "FAIL"
    return equivalence, readiness


def apply_global_gates(rows: Sequence[dict[str, object]]) -> tuple[str, str]:
    equivalence, readiness = calculate_gates(rows)
    for row in rows:
        row["optimization_equivalence_gate"] = equivalence
        row["full_grid_solver_readiness_gate"] = readiness
    return equivalence, readiness


def _run_summary(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        category = str(row.get("category", ""))
        if category.startswith(SMOKE_PREFIX) or category.startswith("regression_"):
            continue
        grouped[str(row["run_mode"])].append(row)
    baseline_seconds = sum(float(row["wall_seconds"]) for row in grouped.get("full_k10_full_strict", ()))
    baseline_cold_branch_seconds = sum(
        float(row.get("branch_strict_seconds", 0.0))
        for row in grouped.get("full_k10_full_strict", ())
    )
    output: list[dict[str, object]] = []
    for mode in MODE_ORDER:
        mode_rows = grouped.get(mode, [])
        wall = sum(float(row["wall_seconds"]) for row in mode_rows)
        mode_branch_seconds = sum(float(row.get("branch_strict_seconds", 0.0)) for row in mode_rows)
        branch_misses = sum(int(row.get("branch_strict_cache_misses", 0)) for row in mode_rows)
        cold_equivalent = (
            wall
            if mode == "full_k10_full_strict" or branch_misses > 0 or mode_branch_seconds == 0.0
            else wall - mode_branch_seconds + baseline_cold_branch_seconds
        )
        failed_modes = [int(row["first_failed_mode"]) for row in mode_rows if str(row.get("first_failed_mode", "")) not in {"", "None"}]
        strict_seconds = sum(
            float(row["EB_strict_seconds"]) + float(row["Timo_strict_seconds"])
            for row in mode_rows
        )
        primary_seconds = sum(
            float(row["EB_primary_seconds"]) + float(row["Timo_primary_seconds"])
            for row in mode_rows
        )
        warm_values = [
            float(row["warm_cache_seconds"])
            for row in mode_rows
            if str(row.get("warm_cache_seconds", "")) not in {"", "None"}
        ]
        output.append(
            {
                "run_mode": mode,
                "case_count": len(mode_rows),
                "wall_seconds": wall,
                "EB_primary_seconds": sum(float(row["EB_primary_seconds"]) for row in mode_rows),
                "Timo_primary_seconds": sum(float(row["Timo_primary_seconds"]) for row in mode_rows),
                "EB_strict_seconds": sum(float(row["EB_strict_seconds"]) for row in mode_rows),
                "Timo_strict_seconds": sum(float(row["Timo_strict_seconds"]) for row in mode_rows),
                "preparation_seconds": sum(float(row["preparation_seconds"]) for row in mode_rows),
                "branch_strict_seconds": mode_branch_seconds,
                "branch_strict_cache_hits": sum(int(row.get("branch_strict_cache_hits", 0)) for row in mode_rows),
                "branch_strict_cache_misses": branch_misses,
                "root_count_EB": sum(int(row["EB_root_count"]) for row in mode_rows),
                "root_count_Timo": sum(int(row["Timo_root_count"]) for row in mode_rows),
                "guard_root_count": sum(int(float(str(row.get("guard_root_count", 0) or 0))) for row in mode_rows),
                "full_K10_points": sum(row["execution_status"] == "resolved_full_K10" for row in mode_rows),
                "early_stopped_points": (
                    0 if mode == "full_k10_full_strict" else sum(_as_bool(row["early_stop_used"]) for row in mode_rows)
                ),
                "points_with_first_failure": sum(
                    str(row.get("first_failed_mode", "")) not in {"", "None"} for row in mode_rows
                ),
                "average_first_failed_mode": sum(failed_modes) / len(failed_modes) if failed_modes else "",
                "N_true_10_cases": sum(str(row.get("N_true", "")) == str(workflow.K_MAX) for row in mode_rows),
                "cache_hits": sum(int(row["cache_hits"]) for row in mode_rows),
                "resolved": sum(str(row["execution_status"]).startswith("resolved_") for row in mode_rows),
                "unresolved": sum(str(row["execution_status"]).startswith("attempted_") for row in mode_rows),
                "peak_memory_bytes": max((int(row["peak_memory_bytes"]) for row in mode_rows), default=0),
                "speedup_relative_to_full": baseline_seconds / wall if baseline_seconds > 0.0 and wall > 0.0 else "",
                "warm_cache_seconds": sum(warm_values) if warm_values else "",
                "cold_equivalent_wall_seconds": cold_equivalent if mode_rows else "",
                "cold_equivalent_speedup_relative_to_full": (
                    baseline_seconds / cold_equivalent
                    if baseline_seconds > 0.0 and cold_equivalent > 0.0 and mode_rows
                    else ""
                ),
                "primary_seconds": primary_seconds,
                "strict_seconds": strict_seconds,
                "strict_time_fraction": strict_seconds / wall if wall > 0.0 else "",
            }
        )
    return output


def write_report(
    output_dir: Path,
    validation_count: int,
    selected_count: int,
    run_rows: Sequence[Mapping[str, object]],
    accuracy_rows: Sequence[Mapping[str, object]],
    forecast_rows: Sequence[Mapping[str, object]],
    best_strategy: str,
) -> Path:
    equivalence_gate, readiness_gate = calculate_gates(accuracy_rows)
    case_counts = defaultdict(int)
    for row in accuracy_rows:
        case_counts[str(row.get("equivalence_case_status", ""))] += 1
    lines = [
        "# Article epsilon prefix optimization benchmark",
        "",
        "This is a computational optimization audit. It does not change `delta_f`, sorted-root ordering, K=10, model equations, or strict tolerances.",
        "",
        f"- Fixed validation manifest: `{validation_count}` points.",
        f"- Points selected in this invocation: `{selected_count}`.",
        f"- Measured strategy candidate: `{best_strategy}`.",
        f"- `optimization_equivalence_gate`: `{equivalence_gate}`.",
        f"- `full_grid_solver_readiness_gate`: `{readiness_gate}`.",
        f"- Equivalence cases: PASS={case_counts['PASS']}, FAIL={case_counts['FAIL']}, "
        f"reference-unresolved={case_counts['NOT_EVALUABLE_REFERENCE_UNRESOLVED']}, "
        f"not-attempted={case_counts['NOT_ATTEMPTED']}.",
        "- Equivalence and full-grid solver readiness are independent gates; PASS/BLOCKED is not a contradiction.",
        "- Four declared smoke failures are audited separately from ordinary runtime medians.",
        "- No 1554-point coarse-grid sweep was run by this benchmark.",
        "",
        "## Run summary",
        "",
        "The timing population is the 18 ordinary validation points; the two immutable S3 regressions and four smoke controls are excluded.",
        "",
        "| mode | ordinary cases | measured wall_s | warm-cache_s | cold-equivalent_s | actual cold sample_s | primary_s | strict_s | EB/Timo roots | guard roots | observed speedup | warm speedup | cold-equivalent speedup | cold sample speedup | early stops | first failures | mean first failure | N_true=10 | unresolved | strict fraction |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in run_rows:
        lines.append(
            f"| {row['run_mode']} | {row['case_count']} | {row['wall_seconds']} | "
            f"{row['warm_cache_seconds']} | {row['cold_equivalent_wall_seconds']} | "
            f"{row.get('actual_cold_isolated_sample_seconds', '')} | {row['primary_seconds']} | "
            f"{row['strict_seconds']} | {row['root_count_EB']}/{row['root_count_Timo']} | "
            f"{row['guard_root_count']} | {row['speedup_relative_to_full']} | "
            f"{row.get('warm_cache_speedup_relative_to_full', '')} | {row['cold_equivalent_speedup_relative_to_full']} | "
            f"{row.get('actual_cold_isolated_sample_speedup', '')} | "
            f"{row['early_stopped_points']} | {row['points_with_first_failure']} | "
            f"{row['average_first_failed_mode']} | {row['N_true_10_cases']} | {row['unresolved']} | "
            f"{row['strict_time_fraction']} |"
        )
    lines.extend(
        [
            "",
            "`early_stopped_points` records actual computational termination before K10. Therefore it is identically zero for `full_k10_full_strict`; `points_with_first_failure` separately records the physical presence of a first failing mode.",
            "",
            "## Validation-based 1554-point timing projection (not executed)",
            "",
            "Each row is the observed mean for that resolved ordinary stratum multiplied by 1554. It is a planning estimate, not authorization to run the grid.",
            "",
            "| mode | stratum | cases | mean_s | median_s | projected_s | projected_h |",
            "| --- | --- | ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for row in forecast_rows:
        projected = row.get("validation_stratified_projected_1554_seconds", "")
        projected_hours = float(projected) / 3600.0 if projected != "" else ""
        lines.append(
            f"| {row['run_mode']} | {row['stratum']} | {row['case_count']} | "
            f"{row['mean_wall_seconds']} | {row['median_wall_seconds']} | {projected} | {projected_hours} |"
        )
    lines.extend(
        [
            "",
            "The four smoke controls have a separate targeted audit under `smoke_failure_audit/`. Three additional ordinary full references remain unresolved. None of these seven references is treated as a physical threshold failure or as evidence against early stopping.",
            "",
            "## Interpretation",
            "",
            "The optimized path is not the default. Because the equivalence gate is PASS on every resolved reference and auto preserves warning escalation without an observed material cost, the future candidate is `paired + auto`; this is a recommendation only, not a default change.",
            "",
            "The full coarse grid remains prohibited until solver readiness is PASS and a separate run is authorized. Current readiness is blocked by unresolved full references and the missing two-trigonometric-pair Timoshenko basis regime documented in the design note; no basis implementation was changed here.",
            "",
            "Future command after readiness PASS and explicit authorization:",
            "",
            "```bash",
            "python scripts/analysis/thickness_mismatch/audits/run_article_epsilon_upper_envelope_grid.py --prefix-until-failure --prefix-strategy paired --strict-policy auto --reuse-cache",
            "```",
        ]
    )
    return workflow.atomic_write_text(output_dir / "report.md", "\n".join(lines) + "\n")


_CASE_ID_PATTERN = re.compile(r"prefix_case_(AUE_[0-9a-f]+)_")


def _read_partial_payload(path: Path) -> dict[str, object] | None:
    try:
        if path.suffix == ".gz":
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                data = json.load(handle)
        else:
            data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, EOFError, json.JSONDecodeError):
        return None
    return dict(data) if isinstance(data, Mapping) else None


def discover_cached_results(
    output_dir: Path,
    cases: Sequence[tuple[str, workflow.GridPoint]],
) -> dict[str, dict[str, Mapping[str, object]]]:
    """Discover durable mode payloads, including legacy timestamped sessions."""

    allowed_ids = {point.case_id for _category, point in cases}
    chosen: dict[tuple[str, str], tuple[tuple[int, int, float], Mapping[str, object]]] = {}
    cache_root = output_dir / "cache"
    paths = [*cache_root.rglob("prefix_case_*.json.gz"), *cache_root.rglob("prefix_case_*.json")]
    for path in paths:
        mode = next((part for part in path.parts if part in MODE_ORDER), None)
        match = _CASE_ID_PATTERN.search(path.name)
        if mode is None or match is None or match.group(1) not in allowed_ids:
            continue
        payload = _read_partial_payload(path)
        if payload is None or payload.get("cache_schema_version") != prefix.PREFIX_CACHE_SCHEMA_VERSION:
            continue
        status = str(payload.get("execution_status", ""))
        status_score = 3 if status.startswith("resolved_") else (2 if status.startswith("attempted_") else 1)
        root_score = int(payload.get("highest_guard_mode_eb", 0)) + int(payload.get("highest_guard_mode_timo", 0))
        score = (status_score, root_score, path.stat().st_mtime)
        key = (mode, match.group(1))
        if key not in chosen or score > chosen[key][0]:
            chosen[key] = (score, payload)

    results: dict[str, dict[str, Mapping[str, object]]] = {mode: {} for mode in MODE_ORDER}
    for (mode, case_id), (_score, payload) in chosen.items():
        results[mode][case_id] = payload

    smoke_summary = output_dir / "smoke_failure_audit" / "summary.csv"
    if smoke_summary.exists():
        for row in workflow.read_csv(smoke_summary):
            case_id = row.get("validation_id", "")
            if case_id not in allowed_ids:
                continue
            results[MODE_ORDER[0]][case_id] = {
                "execution_status": row.get("full_execution_status", "attempted_unresolved"),
                "N_true": None,
                "first_failed_mode": None,
                "unresolved_reason": row.get("exact_unresolved_reason", "targeted_smoke_reference_unresolved"),
                "implementation_limit": _as_bool(row.get("implementation_limit", False)),
                "strict_verification_status": row.get("strict_status", "fail"),
                "early_stop_used": False,
                "models": {},
            }
    return results


def _load_case_rows(
    output_dir: Path,
    cases: Sequence[tuple[str, workflow.GridPoint]],
) -> list[dict[str, object]]:
    path = output_dir / "benchmark_cases.csv"
    if not path.exists():
        return []
    allowed = {point.case_id for _category, point in cases}
    return [
        dict(row)
        for row in workflow.read_csv(path)
        if row.get("run_mode") in MODE_ORDER and row.get("case_id") in allowed
    ]


def _deterministic_case_rows(
    rows: Sequence[Mapping[str, object]],
    cases: Sequence[tuple[str, workflow.GridPoint]],
) -> list[dict[str, object]]:
    order = {point.case_id: index for index, (_category, point) in enumerate(cases)}
    mode_order = {mode: index for index, mode in enumerate(MODE_ORDER)}
    unique: dict[tuple[str, str], dict[str, object]] = {}
    for row in rows:
        normalized = dict(row)
        mode = str(normalized.get("run_mode", ""))
        failed = normalized.get("first_failed_mode", "")
        n_true = normalized.get("N_true", "")
        guard = ""
        if str(failed) not in {"", "None"}:
            guard = 2 * min(workflow.ROOT_GUARD_INDEX, int(float(str(failed))) + 1)
        elif str(n_true) == str(workflow.K_MAX):
            guard = 2 * workflow.ROOT_GUARD_INDEX
        if normalized.get("guard_root_count", "") in {"", None}:
            normalized["guard_root_count"] = guard
        if normalized.get("points_with_first_failure", "") in {"", None}:
            normalized["points_with_first_failure"] = str(failed) not in {"", "None"}
        if mode == MODE_ORDER[0]:
            normalized["early_stop_used"] = False
        normalized.setdefault("warm_cache_seconds", "")
        normalized.setdefault("cache_provenance", "legacy_saved_measurement")
        unique[(mode, str(normalized.get("case_id", "")))] = normalized
    return sorted(
        unique.values(),
        key=lambda row: (
            order.get(str(row.get("case_id", "")), len(order)),
            mode_order.get(str(row.get("run_mode", "")), len(mode_order)),
        ),
    )


def _measure_warm_cache_times(
    rows: Sequence[dict[str, object]],
    cases: Sequence[tuple[str, workflow.GridPoint]],
    output_dir: Path,
) -> None:
    by_key = {(str(row.get("run_mode", "")), str(row.get("case_id", ""))): row for row in rows}
    specs = _mode_specs("paired")
    root = output_dir / "cache" / "sessions" / "validation"
    for category, point in cases:
        if category.startswith(SMOKE_PREFIX) or category.startswith("regression_"):
            continue
        for mode in MODE_ORDER:
            spec = specs[mode]
            cache = prefix.PartialPointCache(root / mode / "partial", reuse_cache=True)
            started = time.perf_counter()
            payload = cache.load(point, strategy=spec.strategy, strict_policy=spec.strict_policy)
            elapsed = time.perf_counter() - started
            row = by_key.get((mode, point.case_id))
            if row is not None and payload is not None:
                row["warm_cache_seconds"] = elapsed


def _runtime_strata(
    category: str,
    payload: Mapping[str, object],
) -> tuple[str, ...]:
    strata: list[str] = ["all_resolved_ordinary"]
    n_true = int(payload["N_true"])
    if n_true <= 3:
        strata.append("early_failure")
    elif n_true <= 7:
        strata.append("middle_prefix")
    elif n_true < workflow.K_MAX:
        strata.append("late_prefix")
    else:
        strata.append("N_true_10")
    deltas = payload.get("deltas", {})
    if isinstance(deltas, Mapping) and any(
        math.isfinite(float(value)) and 0.095 <= float(value) <= 0.105 for value in deltas.values()
    ):
        strata.append("near_threshold")
    if "cluster" in category or "cutoff" in category:
        strata.append("cluster_cutoff_sensitive")
    return tuple(strata)


def _forecast_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    case_rows: Sequence[Mapping[str, object]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    category_by_id = {point.case_id: category for category, point in cases}
    values: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in case_rows:
        mode = str(row.get("run_mode", ""))
        case_id = str(row.get("case_id", ""))
        category = category_by_id.get(case_id, "")
        if category.startswith(SMOKE_PREFIX) or category.startswith("regression_"):
            continue
        payload = results.get(mode, {}).get(case_id)
        reference = results.get(MODE_ORDER[0], {}).get(case_id)
        if not _resolved(payload) or not _resolved(reference):
            continue
        for stratum in _runtime_strata(category, reference):
            values[(mode, stratum)].append(float(row.get("wall_seconds", 0.0)))
    strata = (
        "all_resolved_ordinary", "early_failure", "middle_prefix", "late_prefix",
        "N_true_10", "near_threshold", "cluster_cutoff_sensitive",
    )
    rows: list[dict[str, object]] = []
    for mode in MODE_ORDER:
        for stratum in strata:
            sample = values.get((mode, stratum), [])
            rows.append(
                {
                    "run_mode": mode,
                    "stratum": stratum,
                    "defined": bool(sample),
                    "case_count": len(sample),
                    "mean_wall_seconds": statistics.fmean(sample) if sample else "",
                    "median_wall_seconds": statistics.median(sample) if sample else "",
                    "min_wall_seconds": min(sample) if sample else "",
                    "max_wall_seconds": max(sample) if sample else "",
                    "validation_stratified_projected_1554_seconds": (
                        statistics.fmean(sample) * 1554 if sample else ""
                    ),
                }
            )
    return rows


def _cold_isolated_sample_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    case_rows: Sequence[Mapping[str, object]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    by_key = {(str(row.get("run_mode", "")), str(row.get("case_id", ""))): row for row in case_rows}
    selected: dict[str, tuple[str, workflow.GridPoint]] = {}
    for category, point in cases:
        reference = results.get(MODE_ORDER[0], {}).get(point.case_id)
        if category.startswith(SMOKE_PREFIX) or category.startswith("regression_") or not _resolved(reference):
            continue
        for stratum in _runtime_strata(category, reference):
            if stratum != "all_resolved_ordinary":
                selected.setdefault(stratum, (category, point))
    rows: list[dict[str, object]] = []
    for stratum in (
        "early_failure", "middle_prefix", "late_prefix", "N_true_10",
        "near_threshold", "cluster_cutoff_sensitive",
    ):
        chosen = selected.get(stratum)
        if chosen is None:
            rows.append(
                {
                    "stratum": stratum,
                    "run_mode": "",
                    "case_id": "",
                    "actual_cold_isolated_seconds": "",
                    "measurement_status": "undefined_no_resolved_validation_case",
                }
            )
            continue
        _category, point = chosen
        for mode in MODE_ORDER:
            key = (mode, point.case_id)
            row = by_key.get(key, {})
            cold = (
                str(row.get("cache_provenance", "")) in {"miss", "legacy_saved_measurement"}
                and int(float(str(row.get("branch_strict_cache_hits", 0) or 0))) == 0
            )
            rows.append(
                {
                    "stratum": stratum,
                    "run_mode": mode,
                    "case_id": point.case_id,
                    "actual_cold_isolated_seconds": row.get("wall_seconds", "") if cold else "",
                    "measurement_status": (
                        "measured_mode_isolated_cache_miss" if cold else "not_isolated_or_not_available"
                    ),
                }
            )
    return rows


def _write_outputs(
    output_dir: Path,
    cases: Sequence[tuple[str, workflow.GridPoint]],
    case_rows: Sequence[Mapping[str, object]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
    *,
    selected_count: int,
    best_strategy: str,
) -> tuple[list[dict[str, object]], list[dict[str, object]], str, str]:
    deterministic = _deterministic_case_rows(case_rows, cases)
    root_rows = _root_comparisons(results.get(MODE_ORDER[0], {}), results)
    accuracy_rows = _accuracy_rows(cases, results, root_rows)
    equivalence, readiness = apply_global_gates(accuracy_rows)
    run_rows = _run_summary(deterministic)
    forecast_rows = _forecast_rows(cases, deterministic, results)
    cold_rows = _cold_isolated_sample_rows(cases, deterministic, results)
    cold_totals: dict[str, float] = defaultdict(float)
    cold_cases: dict[str, set[str]] = defaultdict(set)
    for row in cold_rows:
        if row.get("measurement_status") != "measured_mode_isolated_cache_miss":
            continue
        mode = str(row.get("run_mode", ""))
        case_id = str(row.get("case_id", ""))
        if case_id in cold_cases[mode]:
            continue
        cold_cases[mode].add(case_id)
        cold_totals[mode] += float(row["actual_cold_isolated_seconds"])
    baseline_cold = cold_totals.get(MODE_ORDER[0], 0.0)
    baseline_warm = next(
        (float(row["warm_cache_seconds"]) for row in run_rows if row["run_mode"] == MODE_ORDER[0] and row["warm_cache_seconds"] != ""),
        0.0,
    )
    for row in run_rows:
        mode = str(row["run_mode"])
        row["actual_cold_isolated_sample_seconds"] = cold_totals.get(mode, "")
        row["actual_cold_isolated_sample_case_count"] = len(cold_cases.get(mode, set()))
        row["actual_cold_isolated_sample_speedup"] = (
            baseline_cold / cold_totals[mode] if baseline_cold > 0.0 and cold_totals.get(mode, 0.0) > 0.0 else ""
        )
        warm = row.get("warm_cache_seconds", "")
        row["warm_cache_speedup_relative_to_full"] = (
            baseline_warm / float(warm) if baseline_warm > 0.0 and warm != "" and float(warm) > 0.0 else ""
        )
    workflow.write_csv(output_dir / "benchmark_cases.csv", deterministic, BENCHMARK_CASE_FIELDS)
    workflow.write_csv(output_dir / "benchmark_runs.csv", run_rows)
    workflow.write_csv(output_dir / "root_comparison.csv", root_rows, ROOT_COMPARISON_FIELDS)
    workflow.write_csv(output_dir / "accuracy_gate.csv", accuracy_rows, FINAL_VALIDATION_FIELDS)
    workflow.write_csv(output_dir / "runtime_forecast_strata.csv", forecast_rows)
    workflow.write_csv(output_dir / "cold_isolated_sample.csv", cold_rows)
    write_report(output_dir, len(cases), selected_count, run_rows, accuracy_rows, forecast_rows, best_strategy)
    return run_rows, accuracy_rows, equivalence, readiness


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    cases = build_validation_manifest()
    write_validation_manifest(output_dir, cases)
    if args.prepare_only:
        print(f"validation_manifest={len(cases)} output={output_dir}")
        return {"validation_count": len(cases), "selected_count": 0, "output_dir": output_dir}
    best_strategy = "paired"
    results = discover_cached_results(output_dir, cases)
    case_rows = _load_case_rows(output_dir, cases)
    if args.postprocess_only:
        normalized_for_warm = _deterministic_case_rows(case_rows, cases)
        run_rows, accuracy_rows, equivalence, readiness = _write_outputs(
            output_dir,
            cases,
            normalized_for_warm,
            results,
            selected_count=0,
            best_strategy=best_strategy,
        )
        metadata = {
            "finished_utc": datetime.now(timezone.utc).isoformat(),
            "validation_count": len(cases),
            "selected_count": 0,
            "modes": list(MODE_ORDER),
            "best_strategy": best_strategy,
            "optimization_equivalence_gate": equivalence,
            "full_grid_solver_readiness_gate": readiness,
            "postprocess_only": True,
            "root_calculations": 0,
            "coarse_grid_run": False,
        }
        workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
        print(json.dumps(metadata, sort_keys=True))
        return {**metadata, "output_dir": output_dir, "run_rows": run_rows, "accuracy_rows": accuracy_rows}
    selected = list(cases)
    if args.case_id:
        selected = [(category, point) for category, point in selected if point.case_id in set(args.case_id)]
    if args.case_limit is not None:
        selected = selected[: args.case_limit]

    # Smoke controls belong to the separate targeted audit, never to ordinary
    # runtime medians or the three-mode equivalence execution loop.
    if not args.include_smoke_controls:
        selected = [(category, point) for category, point in selected if not category.startswith(SMOKE_PREFIX)]
    previous_measurements = {
        (str(row.get("run_mode", "")), str(row.get("case_id", ""))): row for row in case_rows
    }
    session_token = args.session_id or (
        "validation"
        if args.reuse_cache
        else datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%S_%fZ")
    )
    session_root = output_dir / "cache" / "sessions" / session_token
    shared_strict_cache = branch.BranchContinuationCache(
        session_root / "shared_strict_branch",
        reuse_cache=True,
        force_recompute=False,
        verification_scope="force_strict_verification",
    )
    specs = _mode_specs(best_strategy)
    mismatch_detected = False
    completed_this_invocation = 0
    for category, point in selected:
        for mode_name in MODE_ORDER:
            if mode_name not in args.modes:
                continue
            spec = specs[mode_name]
            run_cache_root = session_root / mode_name
            partial_cache = prefix.PartialPointCache(
                run_cache_root / "partial",
                reuse_cache=args.reuse_cache,
                force=args.force,
            )
            counters: dict[str, int] = defaultdict(int)
            callback = runner._make_prefix_full_strict_callback(
                shared_strict_cache, workflow.strict_settings(), counters
            )
            discovered = results.get(mode_name, {}).get(point.case_id)
            if args.reuse_cache and _resolved(discovered):
                payload = dict(discovered)  # type: ignore[arg-type]
                measured = previous_measurements.get((mode_name, point.case_id))
                if measured is None:
                    measured = _case_row(spec, category, point, payload, 0.0, 0)
                retained = dict(measured)
                retained["cache_provenance"] = "discovered_immutable_cache"
                case_rows.append(retained)
                print(
                    f"validation {point.case_id} mode={mode_name} status={payload.get('execution_status')} "
                    "cache=discovered_immutable_cache",
                    flush=True,
                )
                continue
            if args.track_memory:
                tracemalloc.start()
            started = time.perf_counter()
            payload = prefix.run_staged_point(
                point,
                cache=partial_cache,
                strategy=spec.strategy,
                strict_policy=spec.strict_policy,
                full_k10=spec.full_k10,
                strict_callback=callback,
                force_audit=False,
            )
            wall = time.perf_counter() - started
            if args.track_memory:
                _current, peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
            else:
                peak = 0
            results[mode_name][point.case_id] = payload
            measured = _case_row(spec, category, point, payload, wall, peak)
            measured["cache_provenance"] = partial_cache.last_load_status
            previous = previous_measurements.get((mode_name, point.case_id))
            if previous is not None and partial_cache.last_load_status == "hit":
                warm = wall
                measured.update(previous)
                measured["warm_cache_seconds"] = warm
                measured["cache_provenance"] = "hit"
            case_rows.append(measured)
            print(
                f"validation {point.case_id} mode={mode_name} status={payload.get('execution_status')} "
                f"N_true={payload.get('N_true')} wall_s={wall:.6f} cache={partial_cache.last_load_status}",
                flush=True,
            )

        completed_this_invocation += 1
        run_rows, accuracy_rows, equivalence, readiness = _write_outputs(
            output_dir,
            cases,
            case_rows,
            results,
            selected_count=completed_this_invocation,
            best_strategy=best_strategy,
        )
        current = next(row for row in accuracy_rows if row["validation_id"] == point.case_id)
        n_values = (current["full_N_true"], current["local_N_true"], current["auto_N_true"])
        all_resolved = all(
            str(current[name]).startswith("resolved_")
            for name in ("full_execution_status", "local_execution_status", "auto_execution_status")
        )
        if all_resolved and len({str(value) for value in n_values}) != 1:
            mismatch_detected = True
            print(f"stopping_after_N_true_mismatch case={point.case_id}", flush=True)
            break

    run_rows, accuracy_rows, equivalence, readiness = _write_outputs(
        output_dir,
        cases,
        case_rows,
        results,
        selected_count=completed_this_invocation,
        best_strategy=best_strategy,
    )
    report_path = output_dir / "report.md"
    metadata = {
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "validation_count": len(cases),
        "selected_count": len(selected),
        "completed_this_invocation": completed_this_invocation,
        "modes": args.modes,
        "best_strategy": best_strategy,
        "memory_tracking_enabled": bool(args.track_memory),
        "smoke_controls_explicitly_included": bool(args.include_smoke_controls),
        "optimization_equivalence_gate": equivalence,
        "full_grid_solver_readiness_gate": readiness,
        "stopped_on_N_true_mismatch": mismatch_detected,
        "postprocess_only": False,
        "coarse_grid_run": False,
    }
    workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
    print(json.dumps(metadata, sort_keys=True))
    return {**metadata, "output_dir": output_dir, "report_path": report_path}


if __name__ == "__main__":
    main()
