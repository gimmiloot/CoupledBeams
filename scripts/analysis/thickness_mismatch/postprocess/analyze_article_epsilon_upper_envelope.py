"""CSV-only postprocessing for the article epsilon upper-envelope workflow."""

from __future__ import annotations

import argparse
from collections import defaultdict
import gzip
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import time
from typing import Callable, Mapping, Sequence


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.lib import article_epsilon_upper_envelope as workflow  # noqa: E402
from scripts.lib import article_epsilon_prefix_optimization as prefix  # noqa: E402
from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402
from scripts.lib import variable_length_timoshenko as timo  # noqa: E402


PARTIAL_CASE_FIELDS = (
    "case_id",
    "case_identity",
    "grid_group",
    "regression_label",
    "claim_eligible",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "thin_0p1_flag",
    "cache_record_status",
    "execution_status",
    "n_true_status",
    "N_true",
    "first_apparent_failed_mode",
    "first_apparent_failed_delta_f",
    "required_right_guard",
    "required_prefix_guard_status",
    "upper_spectrum_audit_status",
    "upper_spectrum_audit_incomplete",
    "full_K10_control_status",
    "unresolved_classification",
    "unresolved_reason",
    "wall_seconds",
    "cache_path",
)

PARTIAL_UNRESOLVED_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "classification",
    "first_apparent_failed_mode",
    "required_right_guard",
    "EB_primary_root_count",
    "EB_independent_root_count",
    "Timo_primary_root_count",
    "Timo_independent_root_count",
    "EB_primary_independent_agreement_through_guard",
    "Timo_primary_independent_agreement_through_guard",
    "EB_strict_cache_status",
    "Timo_strict_cache_status",
    "EB_strict_agreement_through_guard",
    "Timo_strict_agreement_through_guard",
    "first_disagreement_model",
    "first_disagreement_root",
    "disagreement_affects_required_guard",
    "cluster_multiplicity_resolved_through_guard",
    "basis_supported_through_guard",
    "full_strict_failure_only_above_guard",
    "exact_N_true_recoverable_zero_solve",
    "reclassified_N_true",
    "raw_unresolved_reason",
    "audit_explanation",
)

PARTIAL_EPSILON_FIELDS = (
    "epsilon_0",
    "intended_case_count",
    "resolved_prefix_exact_count",
    "resolved_full_K10_count",
    "attempted_unresolved_count",
    "deferred_expensive_strict_count",
    "interrupted_incomplete_count",
    "not_attempted_count",
    "accepted_exact_N_true_count",
    "provisional_N_true_min_resolved_subset",
    "provisional_N_true_max_resolved_subset",
    "provisional_N_up_lower_estimate",
    "complete_for_claim",
    "partial_result_warning",
)

PARTIAL_RUNTIME_FIELDS = (
    "section",
    "stratum",
    "point_count",
    "wall_seconds_total",
    "wall_seconds_mean",
    "wall_seconds_median",
    "wall_seconds_q75",
    "wall_seconds_p90",
    "wall_seconds_max",
    "primary_EB_seconds_total",
    "primary_Timoshenko_seconds_total",
    "local_independent_verification_seconds_total",
    "full_strict_seconds_total",
    "preparation_seconds_total",
    "orchestration_and_cache_write_upper_bound_seconds",
    "notes",
)

MAIN_PASS_CASE_FIELDS = (
    *PARTIAL_CASE_FIELDS,
    "basis_regimes",
    "first_apparent_failure",
    "required_guard",
    "exact_reason",
    "affected_root_numbers",
    "primary_status",
    "strict_status",
    "cluster_status",
    "runtime_seconds",
    "suggested_audit_category",
)

MAIN_PASS_EPSILON_FIELDS = (
    "epsilon_0",
    "intended_case_count",
    "resolved_case_count",
    "unresolved_case_count",
    "deferred_case_count",
    "interrupted_case_count",
    "not_attempted_case_count",
    "observed_N_true_min",
    "observed_N_true_max",
    "N_up_observed_raw",
    "N_up_observed_monotone",
    "number_of_observed_argmax_cases",
    "thin_le_0p1_resolved_count",
    "thin_gt_0p1_resolved_count",
    *(f"N_true_{value}_count" for value in range(11)),
    "complete_for_upper_envelope",
    "complete_for_distribution",
    "upper_envelope_status",
    "distribution_status",
    "raw_value_status",
    "monotone_value_status",
    "monotone_provenance",
    "monotone_article_facing",
)

MAIN_PASS_ARGMAX_FIELDS = (
    "epsilon_0",
    "N_up_observed_raw",
    "case_id",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "thin_0p1_flag",
    "raw_value_status",
)

COMPLEX_CASE_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "basis_regimes",
    "first_apparent_failure",
    "required_guard",
    "exact_reason",
    "affected_root_numbers",
    "primary_status",
    "strict_status",
    "cluster_status",
    "runtime",
    "suggested_audit_category",
    "source",
)


CASE_SUMMARY_FIELDS = (
    "case_id",
    "case_identity",
    "grid_group",
    "regression_label",
    "claim_eligible",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "s1",
    "s2",
    "s_max",
    "thin_0p1_flag",
    *(f"delta_f_{index}" for index in range(1, 11)),
    "N_true",
    "first_failed_mode",
    "first_failed_delta_f",
    "last_safe_delta_f",
    "late_safe_mode_count",
    "execution_status",
    "prefix_guard_status",
    "prefix_guard_resolved_through",
    "full_K10_guard_status",
    "early_stop_used",
    "EB_K10_guard_status",
    "Timo_K10_guard_status",
    "strict_verification_status",
    "strict_trigger_reasons",
    "strict_policy_actually_used",
    "basis_regimes_encountered",
    "cluster_warnings",
    "cache_provenance",
    "quality_status",
    "unresolved_reason",
    *(f"max_cutoff_ratio_{index}" for index in range(1, 11)),
    "max_cutoff_ratio_first10",
)

EPSILON_SUMMARY_FIELDS = (
    "epsilon_0",
    "grid_group",
    "intended_case_count",
    "resolved_case_count",
    "unresolved_case_count",
    "not_attempted_case_count",
    "N_true_min",
    "N_true_max",
    "N_up_raw",
    "N_up_monotone",
    "argmax_geometry_count",
    "complete_for_upper_envelope",
    "complete_for_distribution",
    "complete_for_claim",
    "required_full_control_count",
    "passed_full_control_count",
    "failed_full_control_count",
)

ARGMAX_FIELDS = (
    "epsilon_0",
    "summary_group",
    "N_up_raw",
    "case_id",
    "grid_group",
    "beta_deg",
    "mu",
    "eta",
    "s_max",
    "thin_0p1_flag",
    "strict_verification_status",
)

REGRESSION_FIELDS = (
    "regression_label",
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "expected_N_true",
    "observed_N_true",
    "N_true_difference",
    "expected_delta_f_5",
    "observed_delta_f_5",
    "delta_f_5_absolute_difference",
    "quality_status",
    "strict_verification_status",
    "regression_status",
)

UNRESOLVED_FIELDS = (
    "case_id",
    "case_identity",
    "grid_group",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "reason",
    "EB_K10_guard_status",
    "Timo_K10_guard_status",
    "strict_verification_status",
)

NOT_ATTEMPTED_FIELDS = (
    "case_id",
    "case_identity",
    "grid_group",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "execution_status",
    "reason",
)

REFINEMENT_FIELDS = (
    "epsilon_left",
    "epsilon_right",
    "raw_changes",
    "monotone_changes",
    "recommended_epsilon_0",
    "case_role",
    "source_case_id",
    "source_epsilon_0",
    "source_N_true",
    "source_beta_deg",
    "source_mu",
    "source_eta",
    "recommended_beta_deg",
    "recommended_mu",
    "recommended_eta",
    "estimated_additional_geometry_points",
    "estimated_additional_spectrum_solves",
    "estimated_additional_seconds",
    "near_threshold",
    "notes",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Postprocess saved article epsilon upper-envelope CSV files.")
    parser.add_argument("--output-dir", type=Path, default=workflow.MAIN_OUTPUT_DIR)
    parser.add_argument(
        "--partial-cache",
        action="store_true",
        help="Build the interruption-safe 1554-row partial report directly from saved prefix caches.",
    )
    return parser.parse_args(argv)


def _float(row: Mapping[str, str], key: str, default: float = float("nan")) -> float:
    try:
        return float(row.get(key, ""))
    except (TypeError, ValueError):
        return default


def _int_or_blank(value: object) -> int | str:
    try:
        return int(value)
    except (TypeError, ValueError):
        return ""


def _finite_int(value: object) -> int | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) else None


def _root_map(rows: Sequence[Mapping[str, str]]) -> dict[int, Mapping[str, str]]:
    result: dict[int, Mapping[str, str]] = {}
    for row in rows:
        try:
            result[int(row["sorted_index"])] = row
        except (KeyError, TypeError, ValueError):
            continue
    return result


def build_case_summary(
    manifest_rows: Sequence[Mapping[str, str]],
    spectra_rows: Sequence[Mapping[str, str]],
    execution_rows: Sequence[Mapping[str, str]] = (),
    control_rows: Sequence[Mapping[str, str]] = (),
) -> list[dict[str, object]]:
    spectra: dict[tuple[str, str], list[Mapping[str, str]]] = defaultdict(list)
    for row in spectra_rows:
        spectra[(row.get("case_id", ""), row.get("model", ""))].append(row)
    execution_by_case = {row.get("case_id", ""): row for row in execution_rows}
    control_by_case = {row.get("case_id", ""): row for row in control_rows}
    output: list[dict[str, object]] = []
    for manifest in manifest_rows:
        case_id = manifest["case_id"]
        eb = _root_map(spectra.get((case_id, complete.MODEL_EB), ()))
        tm = _root_map(spectra.get((case_id, complete.MODEL_TIMO), ()))
        eb_guard = (
            len(eb) >= 11
            and all(index in eb for index in range(1, 12))
            and all(eb[index].get("solver_status") == "resolved" for index in range(1, 12))
            and all(eb[index].get("root_quality") == "pass" for index in range(1, 12))
        )
        tm_guard = (
            len(tm) >= 11
            and all(index in tm for index in range(1, 12))
            and all(tm[index].get("solver_status") == "resolved" for index in range(1, 12))
            and all(tm[index].get("root_quality") == "pass" for index in range(1, 12))
        )
        strict_status = "missing"
        strict_reasons = ""
        source_rows = list(eb.values()) + list(tm.values())
        metadata_row: Mapping[str, str] = source_rows[0] if source_rows else execution_by_case.get(case_id, {})
        execution_status = "not_attempted"
        prefix_guard_status = "not_attempted"
        prefix_guard_resolved_through = 0
        full_guard_status = "not_attempted"
        early_stop_used = False
        strict_policy_used = "not_attempted"
        cache_provenance = ""
        if metadata_row:
            strict_status = metadata_row.get("strict_verification_status", "missing")
            strict_reasons = metadata_row.get("strict_trigger_reasons", "")
            execution_status = metadata_row.get("execution_status", "resolved_full_K10")
            prefix_guard_status = metadata_row.get("prefix_guard_status", "full_K10_guard_resolved")
            prefix_guard_resolved_through = _finite_int(
                metadata_row.get("prefix_guard_resolved_through", "")
            ) or 0
            full_guard_status = metadata_row.get("full_K10_guard_status", "full_K10_guard_resolved")
            early_stop_used = workflow.parse_bool(metadata_row.get("early_stop_used"))
            strict_policy_used = metadata_row.get("strict_policy_effective", "not_attempted")
            cache_provenance = metadata_row.get("cache_provenance", metadata_row.get("cache_source_path", ""))
        control = control_by_case.get(case_id)
        if control:
            full_guard_status = control.get("full_K10_guard_status", full_guard_status)
        resolved_prefix = execution_status in {"resolved_prefix_early_stop", "resolved_full_K10"}
        strict_failed = strict_status in {"fail", "full_strict_failed", "full_strict_not_available"}
        resolved = (resolved_prefix or (eb_guard and tm_guard)) and not strict_failed
        deltas: list[float] = [float("nan")] * workflow.K_MAX
        if resolved:
            for index in range(1, workflow.K_MAX + 1):
                if index in eb and index in tm:
                    deltas[index - 1] = workflow.squared_frequency_delta(
                        _float(eb[index], "Lambda"),
                        _float(tm[index], "Lambda"),
                    )
        cached_n = _finite_int(metadata_row.get("N_true_cached", "")) if metadata_row else None
        n_true: int | str = (
            cached_n
            if resolved and cached_n is not None
            else (workflow.true_safe_prefix(deltas) if resolved and all(math.isfinite(value) for value in deltas) else "")
        )
        first_failed_mode: int | str = ""
        first_failed_delta: float | str = ""
        last_safe_delta: float | str = ""
        late_safe_count: int | str = ""
        if resolved:
            n_value = int(n_true)
            if n_value < workflow.K_MAX:
                first_failed_mode = n_value + 1
                first_failed_delta = (
                    _float(metadata_row, "first_failed_delta_f")
                    if metadata_row
                    else deltas[n_value]
                )
                late_safe_count = sum(
                    math.isfinite(value) and value <= workflow.DELTA_TOL
                    for value in deltas[n_value + 1 :]
                )
            else:
                late_safe_count = 0
            if n_value > 0:
                last_safe_delta = deltas[n_value - 1]
        cutoff_values = [
            _float(tm[index], "max_cutoff_ratio") if index in tm else float("nan")
            for index in range(1, workflow.K_MAX + 1)
        ]
        finite_cutoff = [value for value in cutoff_values if math.isfinite(value)]
        basis_regimes = sorted(
            {
                value
                for row in tm.values()
                for value in (
                    row.get("timo_basis_regime_arm1", ""),
                    row.get("timo_basis_regime_arm2", ""),
                )
                if value
            }
        )
        cluster_warnings = sorted(
            {
                f"{model}:root_{index}"
                for model, roots in (("EB", eb), ("Timo", tm))
                for index, root in roots.items()
                if int(root.get("cluster_size", "1") or 1) > 1
                or root.get("multiplicity_status", "simple_root") != "simple_root"
            }
        )
        unresolved_reasons: list[str] = []
        if execution_status == "not_attempted":
            unresolved_reasons.append("not_attempted")
        elif execution_status == "not_needed_after_envelope_saturation":
            unresolved_reasons.append("not_needed_after_envelope_saturation")
        elif not eb_guard and not resolved_prefix:
            unresolved_reasons.append("EB_root11_or_root_quality_unresolved")
        if not tm_guard and not resolved_prefix:
            unresolved_reasons.append("Timo_root11_or_root_quality_unresolved")
        if strict_failed:
            unresolved_reasons.append("strict_verification_failed")
        recorded_reason = metadata_row.get("unresolved_reason", "") if metadata_row else ""
        if recorded_reason:
            unresolved_reasons.append(recorded_reason)
        row: dict[str, object] = {
            "case_id": case_id,
            "case_identity": manifest["case_identity"],
            "grid_group": manifest["grid_group"],
            "regression_label": manifest.get("regression_label", ""),
            "claim_eligible": workflow.parse_bool(manifest.get("claim_eligible")),
            "epsilon_0": _float(manifest, "epsilon_0"),
            "beta_deg": _float(manifest, "beta_deg"),
            "mu": _float(manifest, "mu"),
            "eta": _float(manifest, "eta"),
            "s1": _float(manifest, "s1"),
            "s2": _float(manifest, "s2"),
            "s_max": _float(manifest, "s_max"),
            "thin_0p1_flag": workflow.parse_bool(manifest.get("thin_0p1_flag")),
            "N_true": n_true,
            "first_failed_mode": first_failed_mode,
            "first_failed_delta_f": first_failed_delta,
            "last_safe_delta_f": last_safe_delta,
            "late_safe_mode_count": late_safe_count,
            "execution_status": execution_status,
            "prefix_guard_status": prefix_guard_status,
            "prefix_guard_resolved_through": prefix_guard_resolved_through,
            "full_K10_guard_status": full_guard_status,
            "early_stop_used": early_stop_used,
            "EB_K10_guard_status": (
                "resolved" if eb_guard else ("not_needed_after_first_failure" if resolved_prefix else "unresolved")
            ),
            "Timo_K10_guard_status": (
                "resolved" if tm_guard else ("not_needed_after_first_failure" if resolved_prefix else "unresolved")
            ),
            "strict_verification_status": strict_status,
            "strict_trigger_reasons": strict_reasons,
            "strict_policy_actually_used": strict_policy_used,
            "basis_regimes_encountered": ";".join(basis_regimes),
            "cluster_warnings": ";".join(cluster_warnings),
            "cache_provenance": cache_provenance,
            "quality_status": (
                "resolved"
                if resolved
                else ("not_attempted" if execution_status in {"not_attempted", "not_needed_after_envelope_saturation"} else "unresolved")
            ),
            "unresolved_reason": ";".join(unresolved_reasons),
            "max_cutoff_ratio_first10": max(finite_cutoff) if finite_cutoff else float("nan"),
        }
        for index, value in enumerate(deltas, start=1):
            row[f"delta_f_{index}"] = value
        for index, value in enumerate(cutoff_values, start=1):
            row[f"max_cutoff_ratio_{index}"] = value
        output.append(row)
    return output


def _group_predicates() -> dict[str, Callable[[Mapping[str, object]], bool]]:
    return {
        "all": lambda row: bool(row["claim_eligible"]),
        "base": lambda row: row["grid_group"] == "base",
        "low-angle": lambda row: row["grid_group"] == "low_angle",
        "thin<=0.1": lambda row: bool(row["claim_eligible"]) and bool(row["thin_0p1_flag"]),
        "thin>0.1": lambda row: bool(row["claim_eligible"]) and not bool(row["thin_0p1_flag"]),
    }


def build_epsilon_summary(
    case_rows: Sequence[Mapping[str, object]],
    control_comparison_rows: Sequence[Mapping[str, str]] = (),
    *,
    controls_required: bool = False,
) -> list[dict[str, object]]:
    epsilon_values = sorted(
        {float(row["epsilon_0"]) for row in case_rows if bool(row["claim_eligible"])}
    )
    predicates = _group_predicates()
    case_by_id = {str(row.get("case_id", "")): row for row in case_rows if row.get("case_id")}
    comparison_by_epsilon: dict[float, list[Mapping[str, str]]] = defaultdict(list)
    for comparison in control_comparison_rows:
        case = case_by_id.get(comparison.get("case_id", ""))
        if case is not None and bool(case["claim_eligible"]):
            comparison_by_epsilon[float(case["epsilon_0"])].append(comparison)
    output: list[dict[str, object]] = []
    for group_name, predicate in predicates.items():
        group_rows: list[dict[str, object]] = []
        for epsilon_0 in epsilon_values:
            intended = [
                row
                for row in case_rows
                if float(row["epsilon_0"]) == epsilon_0 and predicate(row)
            ]
            resolved = [row for row in intended if row["quality_status"] == "resolved"]
            unresolved = [row for row in intended if row["quality_status"] == "unresolved"]
            not_attempted = [row for row in intended if row["quality_status"] == "not_attempted"]
            n_values = [int(row["N_true"]) for row in resolved]
            n_up: int | str = max(n_values) if n_values else ""
            controls = comparison_by_epsilon.get(epsilon_0, [])
            passed_controls = sum(row.get("comparison_status") == "PASS" for row in controls)
            controls_pass = (not controls_required) or (bool(controls) and passed_controls == len(controls))
            group_rows.append(
                {
                    "epsilon_0": epsilon_0,
                    "grid_group": group_name,
                    "intended_case_count": len(intended),
                    "resolved_case_count": len(resolved),
                    "unresolved_case_count": len(unresolved),
                    "not_attempted_case_count": len(not_attempted),
                    "N_true_min": min(n_values) if n_values else "",
                    "N_true_max": max(n_values) if n_values else "",
                    "N_up_raw": n_up,
                    "N_up_monotone": "",
                    "argmax_geometry_count": (
                        sum(int(row["N_true"]) == int(n_up) for row in resolved) if n_up != "" else 0
                    ),
                    "complete_for_upper_envelope": bool(n_up == workflow.K_MAX),
                    "complete_for_distribution": len(intended) > 0 and len(resolved) == len(intended),
                    "complete_for_claim": len(intended) > 0 and len(resolved) == len(intended) and controls_pass,
                    "required_full_control_count": len(controls),
                    "passed_full_control_count": passed_controls,
                    "failed_full_control_count": len(controls) - passed_controls,
                }
            )
        monotone = workflow.suffix_max(
            [float(row["N_up_raw"]) if row["N_up_raw"] != "" else float("nan") for row in group_rows]
        )
        for row, value in zip(group_rows, monotone):
            row["N_up_monotone"] = int(value) if math.isfinite(float(value)) else ""
        output.extend(group_rows)
    return output


def build_argmax_cases(
    case_rows: Sequence[Mapping[str, object]],
    epsilon_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    predicates = _group_predicates()
    output: list[dict[str, object]] = []
    for summary in epsilon_rows:
        n_up = _finite_int(summary["N_up_raw"])
        if n_up is None:
            continue
        group_name = str(summary["grid_group"])
        predicate = predicates[group_name]
        epsilon_0 = float(summary["epsilon_0"])
        for row in case_rows:
            if (
                float(row["epsilon_0"]) == epsilon_0
                and predicate(row)
                and row["quality_status"] == "resolved"
                and int(row["N_true"]) == n_up
            ):
                output.append(
                    {
                        "epsilon_0": epsilon_0,
                        "summary_group": group_name,
                        "N_up_raw": n_up,
                        "case_id": row["case_id"],
                        "grid_group": row["grid_group"],
                        "beta_deg": row["beta_deg"],
                        "mu": row["mu"],
                        "eta": row["eta"],
                        "s_max": row["s_max"],
                        "thin_0p1_flag": row["thin_0p1_flag"],
                        "strict_verification_status": row["strict_verification_status"],
                    }
                )
    return output


def build_regression_checks(case_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    by_label = {str(row["regression_label"]): row for row in case_rows if row["regression_label"]}
    output: list[dict[str, object]] = []
    for label, epsilon_0, beta_deg, mu, eta, expected_n, expected_delta in workflow.REGRESSION_POINTS:
        row = by_label.get(label)
        observed_n = _finite_int(row["N_true"]) if row else None
        observed_delta = float(row["delta_f_5"]) if row and math.isfinite(float(row["delta_f_5"])) else float("nan")
        passed = (
            row is not None
            and row["quality_status"] == "resolved"
            and observed_n == expected_n
            and abs(observed_delta - expected_delta) <= 1.0e-9
        )
        output.append(
            {
                "regression_label": label,
                "case_id": row["case_id"] if row else "",
                "epsilon_0": epsilon_0,
                "beta_deg": beta_deg,
                "mu": mu,
                "eta": eta,
                "expected_N_true": expected_n,
                "observed_N_true": observed_n if observed_n is not None else "",
                "N_true_difference": observed_n - expected_n if observed_n is not None else "",
                "expected_delta_f_5": expected_delta,
                "observed_delta_f_5": observed_delta,
                "delta_f_5_absolute_difference": abs(observed_delta - expected_delta),
                "quality_status": row["quality_status"] if row else "missing",
                "strict_verification_status": row["strict_verification_status"] if row else "missing",
                "regression_status": "pass" if passed else "fail",
            }
        )
    return output


def build_unresolved(case_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "case_id": row["case_id"],
            "case_identity": row["case_identity"],
            "grid_group": row["grid_group"],
            "epsilon_0": row["epsilon_0"],
            "beta_deg": row["beta_deg"],
            "mu": row["mu"],
            "eta": row["eta"],
            "reason": row["unresolved_reason"],
            "EB_K10_guard_status": row["EB_K10_guard_status"],
            "Timo_K10_guard_status": row["Timo_K10_guard_status"],
            "strict_verification_status": row["strict_verification_status"],
        }
        for row in case_rows
        if row["quality_status"] == "unresolved"
    ]


def build_not_attempted(case_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    return [
        {
            "case_id": row["case_id"],
            "case_identity": row["case_identity"],
            "grid_group": row["grid_group"],
            "epsilon_0": row["epsilon_0"],
            "beta_deg": row["beta_deg"],
            "mu": row["mu"],
            "eta": row["eta"],
            "execution_status": row["execution_status"],
            "reason": row["unresolved_reason"],
        }
        for row in case_rows
        if row["quality_status"] == "not_attempted"
    ]


def _suggest_values(value: float, step: float, lower: float, upper: float) -> list[float]:
    return sorted({max(lower, min(upper, value + offset)) for offset in (-step, 0.0, step)})


def build_refinement_proposal(
    case_rows: Sequence[Mapping[str, object]],
    epsilon_rows: Sequence[Mapping[str, object]],
    runtime_rows: Sequence[Mapping[str, str]] = (),
) -> list[dict[str, object]]:
    all_rows = sorted(
        (row for row in epsilon_rows if row["grid_group"] == "all"),
        key=lambda row: float(row["epsilon_0"]),
    )
    output: list[dict[str, object]] = []
    prefix_times = [
        _float(row, "wall_seconds")
        for row in runtime_rows
        if row.get("phase") == "prefix_sweep" and not workflow.parse_bool(row.get("cache_hit_current_invocation"))
    ]
    finite_prefix_times = sorted(value for value in prefix_times if math.isfinite(value) and value >= 0.0)
    median_case_seconds = (
        finite_prefix_times[len(finite_prefix_times) // 2] if finite_prefix_times else float("nan")
    )
    case_by_epsilon: dict[float, list[Mapping[str, object]]] = defaultdict(list)
    for row in case_rows:
        if bool(row["claim_eligible"]) and row["quality_status"] == "resolved":
            case_by_epsilon[float(row["epsilon_0"])].append(row)
    for left, right in zip(all_rows, all_rows[1:]):
        raw_changes = left["N_up_raw"] != right["N_up_raw"]
        monotone_changes = left["N_up_monotone"] != right["N_up_monotone"]
        if not raw_changes and not monotone_changes:
            continue
        epsilon_left = float(left["epsilon_0"])
        epsilon_right = float(right["epsilon_0"])
        midpoint = 0.5 * (epsilon_left + epsilon_right)
        candidates: list[tuple[str, Mapping[str, object]]] = []
        for epsilon_0, summary in ((epsilon_left, left), (epsilon_right, right)):
            n_up = _finite_int(summary["N_up_raw"])
            if n_up is None:
                continue
            for case in case_by_epsilon[epsilon_0]:
                n_true = int(case["N_true"])
                if n_true == n_up:
                    candidates.append(("envelope", case))
                elif n_true == n_up - 1:
                    candidates.append(("one_mode_below", case))
        seen: set[tuple[str, str]] = set()
        for role, case in candidates:
            key = (role, str(case["case_id"]))
            if key in seen:
                continue
            seen.add(key)
            beta_values = _suggest_values(float(case["beta_deg"]), 7.5, 0.0, 90.0)
            mu_values = _suggest_values(float(case["mu"]), 0.1, 0.0, 0.9)
            eta_values = _suggest_values(float(case["eta"]), 0.125, -0.5, 0.5)
            geometry_count = len(beta_values) * len(mu_values) * len(eta_values)
            output.append(
                {
                    "epsilon_left": epsilon_left,
                    "epsilon_right": epsilon_right,
                    "raw_changes": raw_changes,
                    "monotone_changes": monotone_changes,
                    "recommended_epsilon_0": midpoint,
                    "case_role": role,
                    "source_case_id": case["case_id"],
                    "source_epsilon_0": case["epsilon_0"],
                    "source_N_true": case["N_true"],
                    "source_beta_deg": case["beta_deg"],
                    "source_mu": case["mu"],
                    "source_eta": case["eta"],
                    "recommended_beta_deg": ";".join(f"{value:.12g}" for value in beta_values),
                    "recommended_mu": ";".join(f"{value:.12g}" for value in mu_values),
                    "recommended_eta": ";".join(f"{value:.12g}" for value in eta_values),
                    "estimated_additional_geometry_points": geometry_count,
                    "estimated_additional_spectrum_solves": 2 * geometry_count,
                    "estimated_additional_seconds": geometry_count * median_case_seconds,
                    "near_threshold": 0.095 <= _float({"v": str(case.get("first_failed_delta_f", ""))}, "v") <= 0.105,
                    "notes": "proposal only; no refinement roots were computed",
                }
            )
    epsilon_axis = [float(row["epsilon_0"]) for row in all_rows]
    existing = {(str(row["source_case_id"]), float(row["recommended_epsilon_0"])) for row in output}
    for case in case_rows:
        if not bool(case["claim_eligible"]) or case["quality_status"] != "resolved":
            continue
        delta = case.get("first_failed_delta_f", "")
        try:
            near = 0.095 <= float(delta) <= 0.105
        except (TypeError, ValueError):
            near = False
        if not near or len(epsilon_axis) < 2:
            continue
        epsilon = float(case["epsilon_0"])
        position = epsilon_axis.index(epsilon)
        left_index = max(0, position - 1)
        right_index = min(len(epsilon_axis) - 1, position + 1)
        if left_index == right_index:
            continue
        epsilon_left = epsilon_axis[left_index]
        epsilon_right = epsilon_axis[right_index]
        midpoint = 0.5 * (epsilon_left + epsilon_right)
        if (str(case["case_id"]), midpoint) in existing:
            continue
        beta_values = _suggest_values(float(case["beta_deg"]), 7.5, 0.0, 90.0)
        mu_values = _suggest_values(float(case["mu"]), 0.1, 0.0, 0.9)
        eta_values = _suggest_values(float(case["eta"]), 0.125, -0.5, 0.5)
        geometry_count = len(beta_values) * len(mu_values) * len(eta_values)
        output.append(
            {
                "epsilon_left": epsilon_left,
                "epsilon_right": epsilon_right,
                "raw_changes": False,
                "monotone_changes": False,
                "recommended_epsilon_0": midpoint,
                "case_role": "near_threshold",
                "source_case_id": case["case_id"],
                "source_epsilon_0": case["epsilon_0"],
                "source_N_true": case["N_true"],
                "source_beta_deg": case["beta_deg"],
                "source_mu": case["mu"],
                "source_eta": case["eta"],
                "recommended_beta_deg": ";".join(f"{value:.12g}" for value in beta_values),
                "recommended_mu": ";".join(f"{value:.12g}" for value in mu_values),
                "recommended_eta": ";".join(f"{value:.12g}" for value in eta_values),
                "estimated_additional_geometry_points": geometry_count,
                "estimated_additional_spectrum_solves": 2 * geometry_count,
                "estimated_additional_seconds": geometry_count * median_case_seconds,
                "near_threshold": True,
                "notes": "proposal only; near-threshold control; no refinement roots were computed",
            }
        )
    output.sort(key=lambda row: (float(row["epsilon_left"]), float(row["epsilon_right"]), str(row["case_role"]), str(row["source_case_id"])))
    return output


def _plot_diagnostics(
    output_dir: Path,
    epsilon_rows: Sequence[Mapping[str, object]],
    case_rows: Sequence[Mapping[str, object]],
    runtime_rows: Sequence[Mapping[str, str]] = (),
) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_group: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in epsilon_rows:
        by_group[str(row["grid_group"])].append(row)
    for rows in by_group.values():
        rows.sort(key=lambda row: float(row["epsilon_0"]))
    paths: list[Path] = []

    all_rows = by_group.get("all", [])
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    if all_rows:
        x = [float(row["epsilon_0"]) for row in all_rows]
        ax.plot(x, [_float({"v": str(row["N_up_raw"])}, "v") for row in all_rows], "o-", label="raw")
        ax.plot(x, [_float({"v": str(row["N_up_monotone"])}, "v") for row in all_rows], "s--", label="suffix-max monotone")
    ax.set(xlabel=r"$\epsilon_0$", ylabel=r"$N_{up}$", title="Finite-grid upper envelope")
    ax.set_ylim(-0.25, 10.25)
    ax.grid(True, alpha=0.3)
    ax.legend()
    path = output_dir / "epsilon_upper_envelope_raw_monotone.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig); paths.append(path)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    if all_rows:
        x = [float(row["epsilon_0"]) for row in all_rows]
        low = [_float({"v": str(row["N_true_min"])}, "v") for row in all_rows]
        high = [_float({"v": str(row["N_true_max"])}, "v") for row in all_rows]
        ax.fill_between(x, low, high, alpha=0.25, label="min/max range")
        ax.plot(x, low, "o-", linewidth=1); ax.plot(x, high, "o-", linewidth=1)
    ax.set(xlabel=r"$\epsilon_0$", ylabel=r"$N_{true}$", title="Resolved finite-grid range")
    ax.set_ylim(-0.25, 10.25); ax.grid(True, alpha=0.3); ax.legend()
    path = output_dir / "epsilon_N_true_range.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig); paths.append(path)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for group, style in (("all", "o-"), ("thin<=0.1", "s--"), ("thin>0.1", "^-.")):
        rows = by_group.get(group, [])
        if rows:
            ax.plot(
                [float(row["epsilon_0"]) for row in rows],
                [_float({"v": str(row["N_up_raw"])}, "v") for row in rows],
                style,
                label=group,
            )
    ax.set(xlabel=r"$\epsilon_0$", ylabel=r"$N_{up,raw}$", title="Thinness diagnostic split")
    ax.set_ylim(-0.25, 10.25); ax.grid(True, alpha=0.3); ax.legend()
    path = output_dir / "epsilon_upper_envelope_thin_split.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig); paths.append(path)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    if all_rows:
        x = [float(row["epsilon_0"]) for row in all_rows]
        ax.bar(x, [int(row["resolved_case_count"]) for row in all_rows], width=0.0035, label="resolved")
        ax.bar(
            x,
            [int(row["unresolved_case_count"]) for row in all_rows],
            width=0.0035,
            bottom=[int(row["resolved_case_count"]) for row in all_rows],
            label="unresolved",
        )
        ax.bar(
            x,
            [int(row["not_attempted_case_count"]) for row in all_rows],
            width=0.0035,
            bottom=[
                int(row["resolved_case_count"]) + int(row["unresolved_case_count"])
                for row in all_rows
            ],
            label="not attempted",
        )
    ax.set(xlabel=r"$\epsilon_0$", ylabel="point count", title="Resolution status")
    ax.grid(True, axis="y", alpha=0.3); ax.legend()
    path = output_dir / "epsilon_resolved_unresolved_counts.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig); paths.append(path)

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    epsilon_values = sorted({float(row["epsilon_0"]) for row in case_rows if bool(row["claim_eligible"])})
    bottom = [0] * len(epsilon_values)
    for n_true in range(workflow.K_MAX + 1):
        counts = [
            sum(
                row["quality_status"] == "resolved"
                and int(row["N_true"]) == n_true
                and float(row["epsilon_0"]) == epsilon
                and bool(row["claim_eligible"])
                for row in case_rows
            )
            for epsilon in epsilon_values
        ]
        if any(counts):
            ax.bar(epsilon_values, counts, width=0.0035, bottom=bottom, label=str(n_true))
            bottom = [left + right for left, right in zip(bottom, counts)]
    ax.set(xlabel=r"$\epsilon_0$", ylabel="resolved point count", title=r"Distribution of $N_{true}$")
    ax.grid(True, axis="y", alpha=0.3)
    if epsilon_values:
        ax.legend(title=r"$N_{true}$", ncol=3, fontsize=8)
    path = output_dir / "epsilon_N_true_distribution.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig); paths.append(path)

    fig, (runtime_ax, mode_ax) = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
    prefix_runtime = [row for row in runtime_rows if row.get("phase") == "prefix_sweep"]
    prefix_runtime.sort(key=lambda row: int(row.get("ordinal", "0") or 0))
    if prefix_runtime:
        x = [int(row.get("ordinal", "0") or 0) for row in prefix_runtime]
        runtime_ax.plot(x, [_float(row, "wall_seconds") for row in prefix_runtime], ".", markersize=2)
        mode_ax.plot(x, [_float(row, "first_failed_mode") for row in prefix_runtime], ".", markersize=2)
    runtime_ax.set(ylabel="wall time, s", title="Prefix runtime and first failure")
    mode_ax.set(xlabel="geometry-chain execution ordinal", ylabel="first failed mode")
    runtime_ax.grid(True, alpha=0.3); mode_ax.grid(True, alpha=0.3)
    path = output_dir / "runtime_first_failed_mode.png"
    fig.tight_layout(); fig.savefig(path, dpi=150); plt.close(fig); paths.append(path)
    return paths


def _markdown_table(rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> str:
    header = "| " + " | ".join(fields) + " |"
    separator = "| " + " | ".join("---" for _ in fields) + " |"
    body = [
        "| " + " | ".join(str(row.get(field, "")) for field in fields) + " |"
        for row in rows
    ]
    return "\n".join([header, separator, *body])


def write_report(
    output_dir: Path,
    case_rows: Sequence[Mapping[str, object]],
    epsilon_rows: Sequence[Mapping[str, object]],
    regression_rows: Sequence[Mapping[str, object]],
    unresolved_rows: Sequence[Mapping[str, object]],
    not_attempted_rows: Sequence[Mapping[str, object]],
    refinement_rows: Sequence[Mapping[str, object]],
    control_comparison_rows: Sequence[Mapping[str, str]] = (),
) -> Path:
    all_summary = [row for row in epsilon_rows if row["grid_group"] == "all"]
    thin_summary = [row for row in epsilon_rows if row["grid_group"] in {"thin<=0.1", "thin>0.1"}]
    metadata: Mapping[str, object] = {}
    metadata_path = output_dir / "run_metadata.json"
    if metadata_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            metadata = {}
    text = "\n".join(
        [
            "# Article epsilon upper-envelope coarse-grid diagnostic",
            "",
            "This report is a numerical result on the declared finite geometry grid only. "
            "It is not a proved upper bound on the continuous parameter domain. Sorted frequencies define `N_true`; branch data are diagnostics only.",
            "",
            f"Resolved manifest points: {sum(row['quality_status'] == 'resolved' for row in case_rows)} / {len(case_rows)}.",
            f"Attempted unresolved points: {len(unresolved_rows)}.",
            f"Not-attempted points: {len(not_attempted_rows)}.",
            "The `s_max <= 0.1` and Timoshenko cutoff quantities are diagnostic columns and never filters.",
            "",
            "## Main all-grid envelope",
            "",
            _markdown_table(
                all_summary,
                ("epsilon_0", "resolved_case_count", "unresolved_case_count", "not_attempted_case_count", "N_up_raw", "N_up_monotone", "complete_for_claim"),
            ),
            "",
            "## Thinness diagnostic split",
            "",
            _markdown_table(
                thin_summary,
                ("epsilon_0", "grid_group", "resolved_case_count", "N_up_raw", "N_up_monotone", "complete_for_claim"),
            ),
            "",
            "## Exact regression points",
            "",
            _markdown_table(
                regression_rows,
                ("regression_label", "observed_N_true", "observed_delta_f_5", "delta_f_5_absolute_difference", "regression_status"),
            ),
            "",
            "## Quality and follow-up",
            "",
            f"Strict or primary quality failures are listed in `unresolved_cases.csv` ({len(unresolved_rows)} rows).",
            f"Cases never sent to a root solver are listed separately in `not_attempted_cases.csv` ({len(not_attempted_rows)} rows).",
            f"The coarse-grid-only refinement proposal contains {len(refinement_rows)} rows and was not executed.",
            f"Mandatory full-K10 controls: {len(control_comparison_rows)}; passed: "
            f"{sum(row.get('comparison_status') == 'PASS' for row in control_comparison_rows)}.",
            f"Final validation gate: `{metadata.get('final_validation_gate', 'not_applicable')}`.",
            "Cutoff ratios use the existing model normalization `Omega = epsilon_0*Lambda^2` and the actual arm sections in `variable_length_timoshenko.py`.",
            "",
            "## Run provenance",
            "",
            f"Workflow version: `{metadata.get('workflow_version', workflow.WORKFLOW_VERSION)}`.",
            f"Current invocation operation counts: `{json.dumps(metadata.get('operation_counts_current_invocation', {}), sort_keys=True)}`.",
            f"Timing metadata: `{json.dumps(metadata.get('timings_seconds', {}), sort_keys=True)}`.",
            "",
        ]
    )
    return workflow.atomic_write_text(output_dir / "report.md", text)


def _quantile(values: Sequence[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return float("nan")
    position = (len(ordered) - 1) * float(probability)
    left = int(math.floor(position))
    right = int(math.ceil(position))
    if left == right:
        return ordered[left]
    fraction = position - left
    return ordered[left] * (1.0 - fraction) + ordered[right] * fraction


def _load_prefix_payload(path: Path) -> dict[str, object]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _display_path(path: Path) -> str:
    try:
        value = path.relative_to(REPO_ROOT)
    except ValueError:
        value = path
    return str(value).replace("\\", "/")


def _first_failure(deltas: Mapping[object, object]) -> tuple[int | None, float | None]:
    normalized = {int(key): float(value) for key, value in deltas.items()}
    for index in range(1, workflow.K_MAX + 1):
        value = normalized.get(index)
        if value is not None and math.isfinite(value) and value > workflow.DELTA_TOL:
            return index, value
    return None, None


def _branch_prefix_audit(
    result: branch.BranchContinuationResult | None,
    primary_values: Sequence[float],
    guard: int,
    tolerance: float,
) -> dict[str, object]:
    audited = prefix.assess_branch_strict_through_guard(
        result,
        primary_values,
        guard,
        tolerance,
    )
    return {**audited, "failure_localized": audited.get("first_disagreement") is not None}


def audit_unresolved_prefix_cache(
    point: workflow.GridPoint,
    raw: Mapping[str, object],
    *,
    strict_cache: branch.BranchContinuationCache,
    interrupted: bool = False,
) -> dict[str, object]:
    deltas = dict(raw.get("deltas", {}))
    first_failed, first_delta = _first_failure(deltas)
    guard = first_failed + 1 if first_failed is not None else workflow.ROOT_GUARD_INDEX
    model_audits: dict[str, dict[str, object]] = {}
    first_disagreements: list[tuple[str, int]] = []
    basis_supported = True
    for model in complete.SUPPORTED_MODELS:
        state = raw.get("models", {}).get(model, {}) if isinstance(raw.get("models"), Mapping) else {}
        latest = state.get("latest_result") if isinstance(state, Mapping) else None
        if not isinstance(latest, Mapping):
            model_audits[model] = {
                "primary_root_count": 0,
                "independent_root_count": 0,
                "primary_agreement": False,
                "primary_guard_status": "missing",
                "cluster_resolved": False,
                "strict": _branch_prefix_audit(None, (), guard, complete.DEFAULT_ROOT_MATCH_TOL),
            }
            first_disagreements.append((model, 1))
            continue
        result = prefix._result_from_payload(latest)  # type: ignore[attr-defined]
        assessment = prefix.assess_prefix_guard(result, guard)
        primary_values = tuple(item.Lambda for item in result.primary.roots)
        comparisons = result.primary_vs_verification[:guard]
        primary_agreement = (
            len(comparisons) >= guard
            and all(item.get("status") == "pass" for item in comparisons)
        )
        primary_first = next(
            (
                index
                for index, item in enumerate(comparisons, start=1)
                if item.get("status") != "pass"
            ),
            None,
        )
        if len(comparisons) < guard:
            primary_first = primary_first or len(comparisons) + 1
        strict_result = strict_cache.load(model, point.geometry, workflow.strict_settings())
        strict_audit = _branch_prefix_audit(
            strict_result,
            primary_values,
            guard,
            result.settings.root_match_tol,
        )
        if primary_first is not None:
            first_disagreements.append((model, primary_first))
        if strict_audit["first_disagreement"] is not None:
            first_disagreements.append((model, int(strict_audit["first_disagreement"])))
        cluster_resolved = "unresolved_cluster_multiplicity" not in assessment.reasons
        if model == complete.MODEL_TIMO:
            sections = (
                timo.section_from_epsilon_tau(point.epsilon_0, point.tau1),
                timo.section_from_epsilon_tau(point.epsilon_0, point.tau2),
            )
            for value in primary_values[:guard]:
                for section in sections:
                    try:
                        timo.timo_basis(value, point.epsilon_0, section)
                    except (ArithmeticError, ValueError):
                        basis_supported = False
        model_audits[model] = {
            "primary_root_count": len(result.primary.roots),
            "independent_root_count": len(result.verification.roots),
            "primary_agreement": primary_agreement,
            "primary_guard_status": assessment.status,
            "primary_reasons": assessment.reasons,
            "cluster_resolved": cluster_resolved,
            "strict": strict_audit,
        }
    raw_strict_status = str(raw.get("strict_verification_status", "not_attempted"))
    full_strict_failed = raw_strict_status in {"full_strict_failed", "fail"}
    all_primary = all(
        audit["primary_guard_status"] == "prefix_guard_resolved"
        and bool(audit["primary_agreement"])
        and bool(audit["cluster_resolved"])
        for audit in model_audits.values()
    )
    all_strict = all(bool(audit["strict"]["agreement"]) for audit in model_audits.values())
    delta_contract = (
        first_failed is not None
        and all(
            index in {int(key) for key in deltas}
            and float(deltas.get(str(index), deltas.get(index))) <= workflow.DELTA_TOL
            for index in range(1, first_failed)
        )
        and first_delta is not None
        and first_delta > workflow.DELTA_TOL
    )
    first_pair = min(first_disagreements, key=lambda item: item[1]) if first_disagreements else None
    full_failure_above_guard = bool(
        full_strict_failed
        and delta_contract
        and all_primary
        and all_strict
        and basis_supported
        and (first_pair is None or first_pair[1] > guard)
    )
    exact = bool(delta_contract and all_primary and all_strict and basis_supported and full_failure_above_guard)
    if interrupted:
        classification = "stale_or_interrupted_cache"
    elif exact:
        classification = "upper_spectrum_audit_incomplete"
    elif not all_primary or not all_strict or not delta_contract or not basis_supported:
        classification = "prefix_affecting_unresolved"
    else:
        classification = "other"
    explanation_parts: list[str] = []
    if first_failed is None:
        explanation_parts.append("no confirmed first failure is available")
    if not all_primary:
        explanation_parts.append("primary/independent prefix guard is not complete")
    if not all_strict:
        explanation_parts.append("strict disagreement reaches the required guard or is not localized above it")
    if not basis_supported:
        explanation_parts.append("basis evaluation is not supported through the required guard")
    if exact:
        explanation_parts.append("all exact-prefix criteria pass; only the optional upper audit is incomplete")
    if interrupted:
        explanation_parts.append("cache was saved during an interrupted point and has no completed runtime row")
    eb = model_audits.get(complete.MODEL_EB, {})
    tm = model_audits.get(complete.MODEL_TIMO, {})
    return {
        "case_id": point.case_id,
        "epsilon_0": point.epsilon_0,
        "beta_deg": point.beta_deg,
        "mu": point.mu,
        "eta": point.eta,
        "classification": classification,
        "first_apparent_failed_mode": first_failed if first_failed is not None else "",
        "first_apparent_failed_delta_f": first_delta if first_delta is not None else "",
        "required_right_guard": guard,
        "EB_primary_root_count": eb.get("primary_root_count", 0),
        "EB_independent_root_count": eb.get("independent_root_count", 0),
        "Timo_primary_root_count": tm.get("primary_root_count", 0),
        "Timo_independent_root_count": tm.get("independent_root_count", 0),
        "EB_primary_independent_agreement_through_guard": eb.get("primary_agreement", False),
        "Timo_primary_independent_agreement_through_guard": tm.get("primary_agreement", False),
        "EB_strict_cache_status": eb.get("strict", {}).get("cache_status", "miss"),
        "Timo_strict_cache_status": tm.get("strict", {}).get("cache_status", "miss"),
        "EB_strict_agreement_through_guard": eb.get("strict", {}).get("agreement", False),
        "Timo_strict_agreement_through_guard": tm.get("strict", {}).get("agreement", False),
        "first_disagreement_model": first_pair[0] if first_pair else "",
        "first_disagreement_root": first_pair[1] if first_pair else "",
        "disagreement_affects_required_guard": bool(first_pair and first_pair[1] <= guard),
        "cluster_multiplicity_resolved_through_guard": all(
            bool(audit.get("cluster_resolved")) and bool(audit.get("strict", {}).get("cluster_resolved"))
            for audit in model_audits.values()
        ),
        "basis_supported_through_guard": basis_supported,
        "full_strict_failure_only_above_guard": full_failure_above_guard,
        "exact_N_true_recoverable_zero_solve": exact,
        "reclassified_N_true": first_failed - 1 if exact and first_failed is not None else "",
        "raw_unresolved_reason": raw.get("unresolved_reason", ""),
        "audit_explanation": "; ".join(explanation_parts) or "conservative unresolved classification",
    }


def _runtime_components(raw: Mapping[str, object]) -> dict[str, float]:
    output = {
        "primary_EB": 0.0,
        "primary_Timoshenko": 0.0,
        "verification": 0.0,
        "preparation": 0.0,
        "full_strict": 0.0,
    }
    for row in raw.get("stage_timings", ()):  # type: ignore[union-attr]
        if not isinstance(row, Mapping):
            continue
        model = str(row.get("model", ""))
        key = "primary_EB" if model == complete.MODEL_EB else "primary_Timoshenko"
        output[key] += float(row.get("primary_seconds", 0.0))
        output["verification"] += float(row.get("verification_seconds", 0.0))
        output["preparation"] += float(row.get("preparation_seconds", 0.0))
    details = raw.get("strict_details", {})
    models = details.get("models", {}) if isinstance(details, Mapping) else {}
    if isinstance(models, Mapping):
        output["full_strict"] = sum(
            float(item.get("seconds", 0.0))
            for item in models.values()
            if isinstance(item, Mapping)
        )
    return output


def _runtime_summary_row(
    section: str,
    stratum: str,
    rows: Sequence[Mapping[str, object]],
    *,
    notes: str = "",
) -> dict[str, object]:
    walls = [float(row["wall_seconds"]) for row in rows]
    components = {
        key: sum(float(row.get(key, 0.0)) for row in rows)
        for key in ("primary_EB", "primary_Timoshenko", "verification", "full_strict", "preparation")
    }
    measured = sum(components.values())
    wall_total = sum(walls)
    return {
        "section": section,
        "stratum": stratum,
        "point_count": len(rows),
        "wall_seconds_total": wall_total,
        "wall_seconds_mean": statistics.fmean(walls) if walls else float("nan"),
        "wall_seconds_median": statistics.median(walls) if walls else float("nan"),
        "wall_seconds_q75": _quantile(walls, 0.75),
        "wall_seconds_p90": _quantile(walls, 0.90),
        "wall_seconds_max": max(walls) if walls else float("nan"),
        "primary_EB_seconds_total": components["primary_EB"],
        "primary_Timoshenko_seconds_total": components["primary_Timoshenko"],
        "local_independent_verification_seconds_total": components["verification"],
        "full_strict_seconds_total": components["full_strict"],
        "preparation_seconds_total": components["preparation"],
        "orchestration_and_cache_write_upper_bound_seconds": max(0.0, wall_total - measured),
        "notes": notes,
    }


def analyze_partial_cache(output_dir: Path) -> dict[str, object]:
    """Produce the interruption report without invoking any root resolver."""

    started = time.perf_counter()
    output_dir = output_dir if output_dir.is_absolute() else REPO_ROOT / output_dir
    manifest = workflow.build_manifest()
    runtime_path = output_dir / "runtime_by_case.csv"
    if not runtime_path.exists():
        raise FileNotFoundError("partial cache postprocessing requires runtime_by_case.csv")
    runtime_rows = [
        row for row in workflow.read_csv(runtime_path) if row.get("phase") == "prefix_sweep"
    ]
    runtime_by_case = {row["case_id"]: row for row in runtime_rows}
    point_cache = prefix.PartialPointCache(
        output_dir / "cache" / "prefix", reuse_cache=True, force=False
    )
    strict_cache = branch.BranchContinuationCache(
        output_dir / "cache" / "prefix_strict_branch" / "paired" / "auto",
        reuse_cache=True,
        force_recompute=False,
        verification_scope="force_strict_verification",
    )
    case_rows: list[dict[str, object]] = []
    unresolved_rows: list[dict[str, object]] = []
    detailed_runtime: list[dict[str, object]] = []
    reclassified = 0
    cache_schema_versions: set[str] = set()
    scientific_workflow_versions: set[str] = set()
    for ordinal, point in enumerate(manifest, start=1):
        runtime = runtime_by_case.get(point.case_id)
        cache_path = point_cache.path(point, strategy="paired", strict_policy="main-pass")
        raw: dict[str, object] | None = None
        if cache_path.exists():
            raw = _load_prefix_payload(cache_path)
        else:
            cache_path = point_cache.path(point, strategy="paired", strict_policy="auto")
            if runtime is not None or cache_path.exists():
                raw = _load_prefix_payload(cache_path) if cache_path.exists() else None
        raw_status = str(raw.get("execution_status", "not_attempted")) if raw else "not_attempted"
        interrupted = (
            runtime is None
            and raw is not None
            and raw_status != "deferred_expensive_strict"
        )
        if raw is not None:
            cache_schema_versions.add(str(raw.get("cache_schema_version", "missing")))
            identity = raw.get("identity", {})
            scientific = (
                identity.get("scientific_model_configuration", {})
                if isinstance(identity, Mapping)
                else {}
            )
            if isinstance(scientific, Mapping):
                scientific_workflow_versions.add(str(scientific.get("workflow_version", "missing")))
        n_true: object = runtime.get("N_true", "") if runtime else ""
        first_failed: object = runtime.get("first_failed_mode", "") if runtime else ""
        first_delta: object = ""
        classification = ""
        reason = str(raw.get("unresolved_reason", "")) if raw else ""
        if raw_status == "deferred_expensive_strict":
            execution_status = "deferred_expensive_strict"
            n_status = "unresolved_pending_complex_pass"
            n_true = ""
            prefix_status = "unresolved_without_expensive_strict"
            upper_status = "not_attempted"
            full_status = "not_attempted"
            classification = "expensive_strict_required"
            first_failed = raw.get(
                "first_apparent_failed_mode", raw.get("first_failed_mode", "")
            ) if raw else ""
            first_delta = raw.get("first_failed_delta_f", "") if raw else ""
            guard = raw.get("required_guard_mode", "") if raw else ""
            reason = str(raw.get("defer_reason", reason)) if raw else reason
        elif interrupted:
            execution_status = "interrupted_incomplete"
            n_status = "not_available_interrupted"
            prefix_status = "interrupted_incomplete"
            upper_status = "not_evaluated"
            full_status = "not_evaluated"
            audit = audit_unresolved_prefix_cache(
                point, raw, strict_cache=strict_cache, interrupted=True
            )
            unresolved_rows.append(audit)
            classification = str(audit["classification"])
            first_failed = audit["first_apparent_failed_mode"]
            first_delta = audit["first_apparent_failed_delta_f"]
            guard: object = audit["required_right_guard"]
        elif runtime is None:
            execution_status = "not_attempted"
            n_status = "not_attempted"
            prefix_status = "not_attempted"
            upper_status = "not_attempted"
            full_status = "not_attempted"
            guard = ""
        elif raw_status == "resolved_prefix_early_stop":
            execution_status = "resolved_prefix_exact"
            n_status = "exact"
            prefix_status = "resolved_through_required_guard"
            optional_upper_status = str(
                raw.get("optional_upper_spectrum_full_audit_status", "not_requested")
            ) if raw else "not_requested"
            upper_status = (
                "incomplete_above_required_guard"
                if optional_upper_status == "incomplete"
                else "not_required_for_exact_prefix"
            )
            strict_status = str(raw.get("strict_verification_status", "")) if raw else ""
            full_status = (
                "full_strict_pass"
                if strict_status == "full_strict_pass"
                else (
                    "failed_above_required_guard"
                    if optional_upper_status == "incomplete"
                    else "not_scheduled"
                )
            )
            guard = int(first_failed) + 1 if str(first_failed) else ""
            first_delta = raw.get("first_failed_delta_f", "") if raw else ""
        elif raw_status == "resolved_full_K10":
            execution_status = "resolved_full_K10"
            n_status = "exact"
            prefix_status = "resolved_through_root_11"
            strict_status = str(raw.get("strict_verification_status", "")) if raw else ""
            upper_status = "complete_full_strict" if strict_status == "full_strict_pass" else "primary_full_K10_resolved"
            full_status = "full_strict_pass" if strict_status == "full_strict_pass" else "primary_K10_guard_only"
            guard = workflow.ROOT_GUARD_INDEX
        else:
            audit = audit_unresolved_prefix_cache(point, raw or {}, strict_cache=strict_cache)
            unresolved_rows.append(audit)
            classification = str(audit["classification"])
            first_failed = audit["first_apparent_failed_mode"]
            first_delta = audit["first_apparent_failed_delta_f"]
            guard = audit["required_right_guard"]
            if bool(audit["exact_N_true_recoverable_zero_solve"]):
                execution_status = "resolved_prefix_exact"
                n_status = "exact_reclassified_zero_solve"
                n_true = audit["reclassified_N_true"]
                prefix_status = "resolved_through_required_guard"
                upper_status = "incomplete_above_required_guard"
                full_status = "failed_above_required_guard"
                reclassified += 1
            else:
                execution_status = (
                    "attempted_error" if raw_status == "attempted_error" else "attempted_unresolved"
                )
                n_status = "unresolved"
                prefix_status = "unresolved_at_or_below_required_guard"
                upper_status = "affects_required_prefix_or_full_guard"
                full_status = "failed_or_not_reached"
        upper_incomplete = upper_status in {
            "incomplete_above_required_guard",
            "affects_required_prefix_or_full_guard",
        }
        case_rows.append(
            {
                **point.manifest_row(),
                "cache_record_status": raw_status,
                "execution_status": execution_status,
                "n_true_status": n_status,
                "N_true": n_true,
                "first_apparent_failed_mode": first_failed,
                "first_apparent_failed_delta_f": first_delta,
                "required_right_guard": guard,
                "required_prefix_guard_status": prefix_status,
                "upper_spectrum_audit_status": upper_status,
                "upper_spectrum_audit_incomplete": upper_incomplete,
                "full_K10_control_status": full_status,
                "unresolved_classification": classification,
                "unresolved_reason": reason,
                "wall_seconds": runtime.get("wall_seconds", "") if runtime else "",
                "cache_path": _display_path(cache_path) if cache_path.exists() else "",
            }
        )
        if runtime is not None and raw is not None:
            components = _runtime_components(raw)
            wall = float(runtime["wall_seconds"])
            detailed_runtime.append(
                {
                    "case_id": point.case_id,
                    "epsilon_0": point.epsilon_0,
                    "beta_deg": point.beta_deg,
                    "mu": point.mu,
                    "eta": point.eta,
                    "execution_status": execution_status,
                    "first_failed_mode": first_failed,
                    "failure_stage": (
                        "stage_1"
                        if str(first_failed) and int(first_failed) <= 4
                        else ("stage_2" if str(first_failed) and int(first_failed) <= 7 else ("stage_3" if str(first_failed) else "N_true_10_or_unresolved"))
                    ),
                    "wall_seconds": wall,
                    **components,
                    "strict_policy_effective": raw.get("strict_policy_effective", "not_attempted"),
                    "strict_verification_status": raw.get("strict_verification_status", "not_attempted"),
                }
            )
        if ordinal % 25 == 0:
            print(f"partial-cache postprocess {ordinal}/{len(manifest)}", flush=True)
    case_rows.sort(key=lambda row: str(row["case_id"]))
    unresolved_rows.sort(key=lambda row: str(row["case_id"]))
    epsilon_rows: list[dict[str, object]] = []
    for epsilon in workflow.EPSILON_VALUES:
        rows = [
            row for row in case_rows
            if bool(row["claim_eligible"]) and float(row["epsilon_0"]) == epsilon
        ]
        exact = [row for row in rows if str(row["n_true_status"]).startswith("exact")]
        n_values = [int(row["N_true"]) for row in exact]
        epsilon_rows.append(
            {
                "epsilon_0": epsilon,
                "intended_case_count": len(rows),
                "resolved_prefix_exact_count": sum(row["execution_status"] == "resolved_prefix_exact" for row in rows),
                "resolved_full_K10_count": sum(row["execution_status"] == "resolved_full_K10" for row in rows),
                "attempted_unresolved_count": sum(row["execution_status"] == "attempted_unresolved" for row in rows),
                "deferred_expensive_strict_count": sum(row["execution_status"] == "deferred_expensive_strict" for row in rows),
                "interrupted_incomplete_count": sum(row["execution_status"] == "interrupted_incomplete" for row in rows),
                "not_attempted_count": sum(row["execution_status"] == "not_attempted" for row in rows),
                "accepted_exact_N_true_count": len(exact),
                "provisional_N_true_min_resolved_subset": min(n_values) if n_values else "",
                "provisional_N_true_max_resolved_subset": max(n_values) if n_values else "",
                "provisional_N_up_lower_estimate": max(n_values) if n_values else "",
                "complete_for_claim": False,
                "partial_result_warning": "incomplete resolved subset; not an article-facing envelope",
            }
        )
    runtime_summary: list[dict[str, object]] = []
    runtime_summary.append(
        _runtime_summary_row(
            "overall", "all_completed", detailed_runtime,
            notes="cache-write time was not instrumented separately; the reported residual is an upper bound including orchestration",
        )
    )
    grouping_specs = (
        ("execution_status", lambda row: str(row["execution_status"])),
        ("first_failed_mode", lambda row: str(row["first_failed_mode"] or "none")),
        ("failure_stage", lambda row: str(row["failure_stage"])),
        ("epsilon_0", lambda row: f"{float(row['epsilon_0']):.3f}"),
        ("beta_deg", lambda row: f"{float(row['beta_deg']):g}"),
        ("mu", lambda row: f"{float(row['mu']):g}"),
        ("eta", lambda row: f"{float(row['eta']):g}"),
    )
    for section, key_fn in grouping_specs:
        groups: dict[str, list[Mapping[str, object]]] = defaultdict(list)
        for row in detailed_runtime:
            groups[key_fn(row)].append(row)
        for key in sorted(groups):
            runtime_summary.append(_runtime_summary_row(section, key, groups[key]))
    for rank, row in enumerate(
        sorted(detailed_runtime, key=lambda item: float(item["wall_seconds"]), reverse=True)[:20],
        start=1,
    ):
        runtime_summary.append(
            _runtime_summary_row(
                "top_20_slowest",
                f"{rank:02d}:{row['case_id']}",
                [row],
                notes=(
                    f"eps={row['epsilon_0']};beta={row['beta_deg']};mu={row['mu']};eta={row['eta']};"
                    f"status={row['execution_status']};strict={row['strict_verification_status']}"
                ),
            )
        )
    postprocess_seconds = time.perf_counter() - started
    runtime_summary.append(
        {
            **_runtime_summary_row("postprocess", "zero_solve_partial", []),
            "wall_seconds_total": postprocess_seconds,
            "wall_seconds_mean": postprocess_seconds,
            "wall_seconds_median": postprocess_seconds,
            "wall_seconds_q75": postprocess_seconds,
            "wall_seconds_p90": postprocess_seconds,
            "wall_seconds_max": postprocess_seconds,
            "notes": "measured zero-solve partial postprocessing wall time",
        }
    )
    workflow.write_csv(output_dir / "partial_case_summary.csv", case_rows, PARTIAL_CASE_FIELDS)
    workflow.write_csv(output_dir / "partial_epsilon_summary.csv", epsilon_rows, PARTIAL_EPSILON_FIELDS)
    workflow.write_csv(output_dir / "partial_unresolved_audit.csv", unresolved_rows, PARTIAL_UNRESOLVED_FIELDS)
    workflow.write_csv(output_dir / "partial_runtime_summary.csv", runtime_summary, PARTIAL_RUNTIME_FIELDS)
    counts = {
        status: sum(row["execution_status"] == status for row in case_rows)
        for status in (
            "resolved_prefix_exact",
            "resolved_full_K10",
            "attempted_unresolved",
            "attempted_error",
            "deferred_expensive_strict",
            "interrupted_incomplete",
            "not_attempted",
        )
    }
    exact_count = sum(str(row["n_true_status"]).startswith("exact") for row in case_rows)
    full_strict_rows = [row for row in detailed_runtime if float(row["full_strict"]) > 0.0]
    full_strict_seconds = sum(float(row["full_strict"]) for row in full_strict_rows)
    prefix_required_full_strict = max(0, len(full_strict_rows) - reclassified)
    upper_only_full_strict = min(len(full_strict_rows), reclassified)
    remaining = counts["not_attempted"] + counts["interrupted_incomplete"]
    walls = [float(row["wall_seconds"]) for row in detailed_runtime]
    estimates = {
        "median_seconds": (statistics.median(walls) * remaining if walls else float("nan")),
        "q75_seconds": _quantile(walls, 0.75) * remaining,
        "conservative_p90_seconds": _quantile(walls, 0.90) * remaining,
    }
    report_lines = [
        "# Partial coarse-grid interruption report",
        "",
        "This is a zero-solve diagnostic of an interrupted run. It is not an article-facing upper envelope.",
        "",
        f"- Manifest points: {len(case_rows)}.",
        f"- Completed runtime rows: {len(runtime_rows)}.",
        f"- Exact accepted N_true values: {exact_count}.",
        f"- Genuinely prefix-unresolved: {counts['attempted_unresolved']}.",
        f"- Deferred before expensive strict: {counts['deferred_expensive_strict']}.",
        f"- Upper-spectrum-only audit failures reclassified: {reclassified}.",
        f"- Interrupted incomplete caches: {counts['interrupted_incomplete']}.",
        f"- Not attempted: {counts['not_attempted']}.",
        "",
        "## Status semantics",
        "",
        "`execution_status`, `n_true_status`, `required_prefix_guard_status`, `upper_spectrum_audit_status`, and `full_K10_control_status` are independent fields. A full audit above the right guard could be incomplete without invalidating an exact prefix, but no saved unresolved case satisfied all ten exact-prefix criteria.",
        "",
        "## Runtime",
        "",
        f"- Observed completed-point wall time: {sum(walls):.3f} s.",
        f"- Local independent verification time: {sum(float(row['verification']) for row in detailed_runtime):.3f} s.",
        f"- Auto full-strict escalations: {len(full_strict_rows)}; measured strict time: {full_strict_seconds:.3f} s.",
        f"- Full-strict escalations needed to decide the required prefix: {prefix_required_full_strict}.",
        f"- Full-strict failures confined to an optional upper audit: {upper_only_full_strict}.",
        f"- Sequential remaining median estimate: {estimates['median_seconds']:.1f} s.",
        f"- Sequential remaining 75% estimate: {estimates['q75_seconds']:.1f} s.",
        f"- Sequential conservative p90 estimate: {estimates['conservative_p90_seconds']:.1f} s.",
        "- Cache writes were not separately instrumented in the interrupted invocation; the runtime CSV reports a conservative orchestration-plus-cache-write residual rather than inventing a false exact number.",
        "",
        "The previous validation forecast underrepresented rare full-strict branch-continuation outliers in thin, high-mu and low-angle strata. See `partial_runtime_summary.csv` for component totals, all requested strata, and the 20 slowest points.",
        "",
        "## Unresolved cache audit",
        "",
        "| case | epsilon | beta | mu | eta | class | apparent k | guard | first disagreement | exact N_true |",
        "|---|---:|---:|---:|---:|---|---:|---:|---|---|",
        *(
            f"| {row['case_id']} | {float(row['epsilon_0']):.3f} | {float(row['beta_deg']):g} | "
            f"{float(row['mu']):g} | {float(row['eta']):g} | {row['classification']} | "
            f"{row['first_apparent_failed_mode']} | {row['required_right_guard']} | "
            f"{row['first_disagreement_model']}:{row['first_disagreement_root']} | "
            f"{row['exact_N_true_recoverable_zero_solve']} |"
            for row in unresolved_rows
        ),
        "",
        "No refinement or remaining coarse-grid point was launched by this postprocessor.",
        "",
    ]
    workflow.atomic_write_text(output_dir / "partial_report.md", "\n".join(report_lines))
    return {
        "root_calculations": 0,
        "manifest_count": len(case_rows),
        "completed_count": len(runtime_rows),
        "exact_N_true_count": exact_count,
        "genuinely_unresolved_count": counts["attempted_unresolved"],
        "upper_only_reclassified_count": reclassified,
        "interrupted_incomplete_count": counts["interrupted_incomplete"],
        "not_attempted_count": counts["not_attempted"],
        "attempted_error_count": counts["attempted_error"],
        "deferred_expensive_strict_count": counts["deferred_expensive_strict"],
        "cache_schema_versions": sorted(cache_schema_versions),
        "scientific_workflow_versions": sorted(scientific_workflow_versions),
        "postprocess_seconds": postprocess_seconds,
        "sequential_estimates": estimates,
    }


def _truthy(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def _main_pass_sort_key(row: Mapping[str, object]) -> tuple[float, float, float, float, str]:
    return (
        float(row.get("epsilon_0", 0.0)),
        float(row.get("beta_deg", 0.0)),
        float(row.get("mu", 0.0)),
        float(row.get("eta", 0.0)),
        str(row.get("case_id", "")),
    )


def main_pass_completeness(
    *,
    intended_count: int,
    resolved_count: int,
    n_up_observed: int | None,
) -> dict[str, object]:
    all_resolved = resolved_count == intended_count
    saturated = n_up_observed == workflow.K_MAX
    raw_status = (
        "incomplete"
        if n_up_observed is None
        else ("exact" if all_resolved else ("exact_by_K_saturation" if saturated else "provisional"))
    )
    return {
        "complete_for_upper_envelope": all_resolved or saturated,
        "complete_for_distribution": all_resolved,
        "upper_envelope_status": (
            "exact"
            if all_resolved
            else ("exact_by_K_saturation" if saturated else "provisional_lower_bound")
        ),
        "distribution_status": "exact" if all_resolved else "incomplete",
        "raw_value_status": raw_status,
    }


def _main_pass_plots(output_dir: Path) -> list[Path]:
    """Create provisional diagnostics strictly from the saved CSV products."""

    import matplotlib  # noqa: WPS433

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: WPS433

    epsilon_rows = workflow.read_csv(output_dir / "main_pass_epsilon_summary.csv")
    case_rows = workflow.read_csv(output_dir / "main_pass_case_summary.csv")
    runtime_rows = [
        row
        for row in workflow.read_csv(output_dir / "runtime_by_case.csv")
        if row.get("phase") == "prefix_sweep" and row.get("wall_seconds") not in {None, ""}
    ]
    plot_dir = output_dir / "main_pass_plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []

    def save(name: str) -> None:
        path = plot_dir / name
        plt.tight_layout()
        plt.savefig(path, dpi=180)
        plt.close()
        paths.append(path)

    eps = [float(row["epsilon_0"]) for row in epsilon_rows]
    raw = [float(row["N_up_observed_raw"]) for row in epsilon_rows]
    monotone = [float(row["N_up_observed_monotone"]) for row in epsilon_rows]
    plt.figure(figsize=(7.2, 4.2))
    plt.step(eps, raw, where="mid", marker="o", label="observed raw")
    plt.step(eps, monotone, where="mid", marker="s", label="observed suffix-max")
    plt.xlabel("epsilon_0")
    plt.ylabel("N_up observed")
    plt.title("Provisional observed upper-envelope diagnostics")
    plt.legend()
    save("01_provisional_observed_N_up_raw_monotone.png")

    plt.figure(figsize=(7.2, 4.2))
    plt.fill_between(
        eps,
        [float(row["observed_N_true_min"]) for row in epsilon_rows],
        [float(row["observed_N_true_max"]) for row in epsilon_rows],
        alpha=0.25,
    )
    plt.plot(eps, [float(row["observed_N_true_min"]) for row in epsilon_rows], marker="o", label="observed min")
    plt.plot(eps, [float(row["observed_N_true_max"]) for row in epsilon_rows], marker="s", label="observed max")
    plt.xlabel("epsilon_0")
    plt.ylabel("N_true")
    plt.title("Provisional observed N_true range")
    plt.legend()
    save("02_provisional_observed_N_true_min_max.png")

    plt.figure(figsize=(7.2, 4.2))
    bottom = [0.0] * len(eps)
    for field, label in (
        ("resolved_case_count", "resolved"),
        ("unresolved_case_count", "unresolved/interrupted"),
        ("deferred_case_count", "deferred"),
    ):
        values = [float(row[field]) for row in epsilon_rows]
        plt.bar(eps, values, bottom=bottom, width=0.0035, label=label)
        bottom = [left + right for left, right in zip(bottom, values)]
    plt.xlabel("epsilon_0")
    plt.ylabel("case count")
    plt.title("Provisional resolved/unresolved/deferred counts")
    plt.legend()
    save("03_provisional_status_counts.png")

    plt.figure(figsize=(8.0, 4.6))
    distributions = [
        [int(row[f"N_true_{value}_count"]) for value in range(11)]
        for row in epsilon_rows
    ]
    image = plt.imshow(distributions, aspect="auto", origin="lower", cmap="viridis")
    plt.colorbar(image, label="resolved case count")
    plt.xticks(range(11), range(11))
    plt.yticks(range(len(eps)), [f"{value:.3f}" for value in eps])
    plt.xlabel("N_true")
    plt.ylabel("epsilon_0")
    plt.title("Provisional N_true distribution on resolved points")
    save("04_provisional_N_true_distribution.png")

    runtime_case = {row["case_id"]: row for row in case_rows}
    runtime_pairs = [
        (runtime_case.get(row["case_id"], {}), float(row["wall_seconds"]))
        for row in runtime_rows
    ]
    plt.figure(figsize=(7.2, 4.2))
    x_failed = [
        float(case.get("first_apparent_failed_mode") or 11)
        for case, _wall in runtime_pairs
    ]
    plt.scatter(x_failed, [wall for _case, wall in runtime_pairs], s=9, alpha=0.55)
    plt.yscale("log")
    plt.xlabel("first_failed_mode (11 means none/unresolved)")
    plt.ylabel("runtime, s")
    plt.title("Provisional runtime against first_failed_mode")
    save("05_provisional_runtime_vs_first_failed_mode.png")

    plt.figure(figsize=(7.2, 4.2))
    plt.scatter(
        [float(case.get("epsilon_0", 0.0)) for case, _wall in runtime_pairs],
        [wall for _case, wall in runtime_pairs],
        s=9,
        alpha=0.55,
    )
    plt.yscale("log")
    plt.xlabel("epsilon_0")
    plt.ylabel("runtime, s")
    plt.title("Provisional runtime against epsilon_0")
    save("06_provisional_runtime_vs_epsilon_0.png")

    plt.figure(figsize=(7.2, 4.2))
    thin = []
    nonthin = []
    for epsilon in eps:
        resolved = [
            row
            for row in case_rows
            if float(row["epsilon_0"]) == epsilon
            and str(row["n_true_status"]).startswith("exact")
        ]
        thin.append(sum(_truthy(row["thin_0p1_flag"]) for row in resolved))
        nonthin.append(sum(not _truthy(row["thin_0p1_flag"]) for row in resolved))
    plt.plot(eps, thin, marker="o", label="thin <= 0.1 resolved")
    plt.plot(eps, nonthin, marker="s", label="thin > 0.1 resolved")
    plt.xlabel("epsilon_0")
    plt.ylabel("resolved case count")
    plt.title("Provisional thin and nonthin resolved subsets")
    plt.legend()
    save("07_provisional_thin_nonthin_resolved.png")
    return paths


def analyze_main_pass(
    output_dir: Path,
    *,
    defer_case_list: Path | None = None,
) -> dict[str, object]:
    """Build the regular-pass products using saved cache/CSV data only."""

    started = time.perf_counter()
    output_dir = output_dir if output_dir.is_absolute() else REPO_ROOT / output_dir
    partial = analyze_partial_cache(output_dir)
    partial_rows = workflow.read_csv(output_dir / "partial_case_summary.csv")
    partial_audit = {
        row["case_id"]: row
        for row in workflow.read_csv(output_dir / "partial_unresolved_audit.csv")
    }
    spectra_rows = workflow.read_csv(output_dir / "spectra_long.csv")
    execution_path = output_dir / "case_execution.csv"
    execution_rows = workflow.read_csv(execution_path) if execution_path.exists() else []
    execution_by_case = {row["case_id"]: row for row in execution_rows}
    metadata_path = output_dir / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8")) if metadata_path.exists() else {}
    initial_complex = set(str(value) for value in metadata.get("main_pass_initial_complex_case_ids", ()))
    deferred_path = defer_case_list or output_dir / "deferred_complex_cases_current.csv"
    if not deferred_path.exists():
        deferred_path = output_dir / "deferred_complex_cases_pre_run.csv"
    deferred_rows = workflow.read_csv(deferred_path) if deferred_path.exists() else []
    deferred_by_case = {row["case_id"]: row for row in deferred_rows}

    regimes: dict[str, set[str]] = defaultdict(set)
    clusters: dict[str, bool] = defaultdict(bool)
    for row in spectra_rows:
        case_id = str(row.get("case_id", ""))
        for field in ("timo_basis_regime_arm1", "timo_basis_regime_arm2"):
            if row.get(field):
                regimes[case_id].add(str(row[field]))
        try:
            clusters[case_id] = clusters[case_id] or int(row.get("cluster_size", 1)) > 1
        except (TypeError, ValueError):
            pass

    case_rows: list[dict[str, object]] = []
    for row in partial_rows:
        case_id = str(row["case_id"])
        status = str(row["execution_status"])
        if status == "not_attempted" and case_id in deferred_by_case:
            status = "deferred_complex"
        audit = partial_audit.get(case_id, {})
        execution = execution_by_case.get(case_id, {})
        first_failure = row.get("first_apparent_failed_mode", "")
        affected = audit.get("first_disagreement_root", "")
        if not affected and first_failure:
            affected = first_failure
        exact_reason = (
            audit.get("audit_explanation")
            or row.get("unresolved_reason")
            or deferred_by_case.get(case_id, {}).get("exact_reason", "")
        )
        strict_status = execution.get("strict_verification_status", "not_attempted")
        if status == "interrupted_incomplete":
            category = "interrupted_cache_recovery"
        elif status == "deferred_expensive_strict":
            category = "force_strict_verification_required"
        elif status == "attempted_error":
            category = "worker_or_case_exception"
        elif case_id in deferred_by_case and deferred_by_case[case_id].get("source") == "parallel_benchmark":
            category = "known_parallel_benchmark_prefix_unresolved"
        elif status == "attempted_unresolved":
            category = "prefix_root_completeness_or_strict_audit"
        elif status == "deferred_complex":
            category = "preconfirmed_complex_case"
        else:
            category = "not_complex"
        enriched = {
            **row,
            "execution_status": status,
            "basis_regimes": ";".join(sorted(regimes.get(case_id, ()))),
            "first_apparent_failure": first_failure,
            "required_guard": row.get("required_right_guard", ""),
            "exact_reason": exact_reason,
            "affected_root_numbers": affected,
            "primary_status": row.get("cache_record_status", "not_attempted"),
            "strict_status": strict_status,
            "cluster_status": (
                "cluster_warning_present" if clusters.get(case_id) else "no_cluster_warning_observed"
            ),
            "runtime_seconds": row.get("wall_seconds", ""),
            "suggested_audit_category": category,
        }
        case_rows.append(enriched)
    case_rows.sort(key=_main_pass_sort_key)

    epsilon_rows: list[dict[str, object]] = []
    for epsilon in workflow.EPSILON_VALUES:
        rows = [
            row
            for row in case_rows
            if _truthy(row["claim_eligible"]) and float(row["epsilon_0"]) == epsilon
        ]
        resolved = [row for row in rows if str(row["n_true_status"]).startswith("exact")]
        values = [int(row["N_true"]) for row in resolved]
        unresolved = [
            row
            for row in rows
            if row["execution_status"] in {
                "attempted_unresolved",
                "attempted_error",
                "interrupted_incomplete",
            }
        ]
        deferred = [
            row
            for row in rows
            if row["execution_status"]
            in {"deferred_complex", "deferred_expensive_strict"}
        ]
        not_attempted = [row for row in rows if row["execution_status"] == "not_attempted"]
        raw_value = max(values) if values else None
        completeness = main_pass_completeness(
            intended_count=len(rows),
            resolved_count=len(resolved),
            n_up_observed=raw_value,
        )
        epsilon_rows.append(
            {
                "epsilon_0": epsilon,
                "intended_case_count": len(rows),
                "resolved_case_count": len(resolved),
                "unresolved_case_count": len(unresolved),
                "deferred_case_count": len(deferred),
                "interrupted_case_count": sum(
                    row["execution_status"] == "interrupted_incomplete" for row in rows
                ),
                "not_attempted_case_count": len(not_attempted),
                "observed_N_true_min": min(values) if values else "",
                "observed_N_true_max": max(values) if values else "",
                "N_up_observed_raw": raw_value if raw_value is not None else "",
                "N_up_observed_monotone": "",
                "number_of_observed_argmax_cases": sum(value == raw_value for value in values) if raw_value is not None else 0,
                "thin_le_0p1_resolved_count": sum(_truthy(row["thin_0p1_flag"]) for row in resolved),
                "thin_gt_0p1_resolved_count": sum(not _truthy(row["thin_0p1_flag"]) for row in resolved),
                **{f"N_true_{value}_count": sum(int(row["N_true"]) == value for row in resolved) for value in range(11)},
                **completeness,
                "monotone_value_status": "",
                "monotone_provenance": "",
                "monotone_article_facing": False,
            }
        )
    raw_values = [
        int(row["N_up_observed_raw"]) if row["N_up_observed_raw"] != "" else -1
        for row in epsilon_rows
    ]
    monotone_values = workflow.suffix_max(raw_values)
    for index, row in enumerate(epsilon_rows):
        row["N_up_observed_monotone"] = monotone_values[index] if monotone_values[index] >= 0 else ""
        suffix_statuses = [str(item["raw_value_status"]) for item in epsilon_rows[index:]]
        provisional_inputs = any(value in {"provisional", "incomplete"} for value in suffix_statuses)
        if monotone_values[index] < 0:
            monotone_status = "incomplete"
        elif provisional_inputs:
            monotone_status = "exact_by_K_saturation" if monotone_values[index] == workflow.K_MAX else "provisional"
        elif all(value == "exact" for value in suffix_statuses):
            monotone_status = "exact"
        else:
            monotone_status = "exact_by_K_saturation"
        row["monotone_value_status"] = monotone_status
        row["monotone_provenance"] = ";".join(sorted(set(suffix_statuses)))
        row["monotone_article_facing"] = not provisional_inputs and monotone_status == "exact"

    argmax_rows: list[dict[str, object]] = []
    for epsilon_row in epsilon_rows:
        epsilon = float(epsilon_row["epsilon_0"])
        raw_value = epsilon_row["N_up_observed_raw"]
        if raw_value == "":
            continue
        for row in case_rows:
            if (
                _truthy(row["claim_eligible"])
                and float(row["epsilon_0"]) == epsilon
                and str(row["n_true_status"]).startswith("exact")
                and int(row["N_true"]) == int(raw_value)
            ):
                argmax_rows.append(
                    {
                        "epsilon_0": epsilon,
                        "N_up_observed_raw": raw_value,
                        "case_id": row["case_id"],
                        "beta_deg": row["beta_deg"],
                        "mu": row["mu"],
                        "eta": row["eta"],
                        "s_max": row["s_max"],
                        "thin_0p1_flag": row["thin_0p1_flag"],
                        "raw_value_status": epsilon_row["raw_value_status"],
                    }
                )
    argmax_rows.sort(key=_main_pass_sort_key)

    complex_ids = {
        str(row["case_id"])
        for row in case_rows
        if row["execution_status"] in {
            "attempted_unresolved",
            "attempted_error",
            "interrupted_incomplete",
            "deferred_complex",
            "deferred_expensive_strict",
        }
    } | set(deferred_by_case)
    complex_rows: list[dict[str, object]] = []
    for row in case_rows:
        case_id = str(row["case_id"])
        if case_id not in complex_ids:
            continue
        complex_rows.append(
            {
                "case_id": case_id,
                "epsilon_0": row["epsilon_0"],
                "beta_deg": row["beta_deg"],
                "mu": row["mu"],
                "eta": row["eta"],
                "s_max": row["s_max"],
                "basis_regimes": row["basis_regimes"],
                "first_apparent_failure": row["first_apparent_failure"],
                "required_guard": row["required_guard"],
                "exact_reason": row["exact_reason"],
                "affected_root_numbers": row["affected_root_numbers"],
                "primary_status": row["primary_status"],
                "strict_status": row["strict_status"],
                "cluster_status": row["cluster_status"],
                "runtime": row["runtime_seconds"],
                "suggested_audit_category": row["suggested_audit_category"],
                "source": deferred_by_case.get(case_id, {}).get("source", "new_main_pass_complex"),
            }
        )
    complex_rows.sort(key=_main_pass_sort_key)
    new_complex_rows = [row for row in complex_rows if row["case_id"] not in initial_complex and row["case_id"] not in deferred_by_case]
    unresolved_rows = [
        row
        for row in complex_rows
        if next(
            (case["execution_status"] for case in case_rows if case["case_id"] == row["case_id"]),
            "deferred_complex",
        )
        in {"attempted_unresolved", "attempted_error", "interrupted_incomplete"}
    ]

    workflow.write_csv(output_dir / "main_pass_case_summary.csv", case_rows, MAIN_PASS_CASE_FIELDS)
    workflow.write_csv(output_dir / "main_pass_epsilon_summary.csv", epsilon_rows, MAIN_PASS_EPSILON_FIELDS)
    workflow.write_csv(output_dir / "main_pass_argmax_cases.csv", argmax_rows, MAIN_PASS_ARGMAX_FIELDS)
    workflow.write_csv(output_dir / "main_pass_unresolved_cases.csv", unresolved_rows, COMPLEX_CASE_FIELDS)
    workflow.write_csv(output_dir / "deferred_complex_cases_post_run.csv", new_complex_rows, COMPLEX_CASE_FIELDS)
    workflow.write_csv(output_dir / "complex_case_manifest.csv", complex_rows, COMPLEX_CASE_FIELDS)

    runtime_summary = workflow.read_csv(output_dir / "partial_runtime_summary.csv")
    workflow.write_csv(output_dir / "main_pass_runtime_summary.csv", runtime_summary, PARTIAL_RUNTIME_FIELDS)
    plot_paths = _main_pass_plots(output_dir)

    manifest_path = output_dir / "grid_manifest.csv"
    manifest_sha = hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    prefix_files = sorted((output_dir / "cache" / "prefix").glob("prefix_case_*.json.gz"))
    cache_names = [path.name for path in prefix_files]
    duplicate_cache_ids = sorted(
        name for name in set(cache_names) if cache_names.count(name) > 1
    )
    temporary_files = [
        path
        for path in (output_dir / "cache").rglob("*")
        if path.is_file() and (path.name.endswith(".tmp") or path.suffix in {".partial", ".lock"})
    ]
    integrity = {
        "manifest_count": len(case_rows),
        "manifest_unique_ids": len({row["case_id"] for row in case_rows}),
        "manifest_sha256": manifest_sha,
        "duplicate_cache_entries": duplicate_cache_ids,
        "temporary_cache_files": [str(path) for path in temporary_files],
        "cache_schema_versions": partial.get("cache_schema_versions", []),
        "scientific_workflow_versions": partial.get("scientific_workflow_versions", []),
        "csv_order_deterministic": case_rows == sorted(case_rows, key=_main_pass_sort_key),
        "postprocess_root_calculations": 0,
    }
    integrity["cache_integrity"] = "PASS" if (
        len(case_rows) == 1554
        and integrity["manifest_unique_ids"] == 1554
        and manifest_sha == "b73747c80637c7ee9d0bf221889f769c65ef694b9d633a3c7e19f0ab76c5673e"
        and not duplicate_cache_ids
        and not temporary_files
        and partial.get("cache_schema_versions") == [prefix.PREFIX_CACHE_SCHEMA_VERSION]
        and partial.get("scientific_workflow_versions") == [workflow.WORKFLOW_VERSION]
    ) else "FAIL"
    workflow.atomic_write_json(output_dir / "main_pass_integrity.json", integrity)

    overall = next(
        (row for row in runtime_summary if row.get("section") == "overall" and row.get("stratum") == "all_completed"),
        {},
    )
    gate_path = output_dir / "resume_readiness_gate.csv"
    gate = workflow.read_csv(gate_path)[0] if gate_path.exists() else {}
    main_timing = metadata.get("timings_seconds", {}) if isinstance(metadata, Mapping) else {}
    wall_time = float(main_timing.get("pre_postprocess_total", 0.0)) if isinstance(main_timing, Mapping) else 0.0
    full_strict_count = sum(float(row.get("full_strict_seconds_total", 0.0)) > 0.0 for row in runtime_summary if row.get("section") == "top_20_slowest")
    resolved_case_rows = [row for row in case_rows if str(row["n_true_status"]).startswith("exact")]
    status_counts = {
        status: sum(row["execution_status"] == status for row in case_rows)
        for status in (
            "resolved_prefix_exact",
            "resolved_full_K10",
            "attempted_unresolved",
            "attempted_error",
            "deferred_complex",
            "deferred_expensive_strict",
            "interrupted_incomplete",
            "not_attempted",
        )
    }
    report_lines = [
        "# Main regular coarse-grid pass report",
        "",
        "This report is preliminary. Complex cases, full-K10 controls, and refinement were not run.",
        "",
        f"- Manifest: {len(case_rows)} unique cases; SHA-256 `{manifest_sha}`.",
        f"- Cache integrity: `{integrity['cache_integrity']}`.",
        f"- Resolved exact: {len(resolved_case_rows)}.",
        f"- Attempted unresolved: {status_counts['attempted_unresolved']}; attempted errors: {status_counts['attempted_error']}.",
        f"- Deferred known: {status_counts['deferred_complex']}; deferred before expensive strict: {status_counts['deferred_expensive_strict']}; interrupted: {status_counts['interrupted_incomplete']}; not attempted: {status_counts['not_attempted']}.",
        f"- New complex cases from this pass: {len(new_complex_rows)}.",
        f"- Main-pass wall time before postprocessing: {wall_time:.3f} s; workers: {metadata.get('run_scope', {}).get('workers', '') if isinstance(metadata.get('run_scope', {}), Mapping) else ''}.",
        f"- Workers=4 benchmark speedup estimate: {gate.get('recommended_speedup', '')}.",
        f"- Primary EB/Timoshenko/local/full-strict seconds: {overall.get('primary_EB_seconds_total', '')} / {overall.get('primary_Timoshenko_seconds_total', '')} / {overall.get('local_independent_verification_seconds_total', '')} / {overall.get('full_strict_seconds_total', '')}.",
        f"- Orchestration/cache upper bound: {overall.get('orchestration_and_cache_write_upper_bound_seconds', '')} s.",
        f"- Early stops: {status_counts['resolved_prefix_exact']}; N_true=10: {sum(int(row['N_true']) == 10 for row in resolved_case_rows)}.",
        "",
        "## Observed epsilon picture",
        "",
        "| epsilon_0 | resolved | unresolved | deferred | raw | monotone | upper complete | distribution complete |",
        "|---:|---:|---:|---:|---:|---:|---|---|",
        *(
            f"| {float(row['epsilon_0']):.3f} | {row['resolved_case_count']} | {row['unresolved_case_count']} | "
            f"{row['deferred_case_count']} | {row['N_up_observed_raw']} ({row['raw_value_status']}) | "
            f"{row['N_up_observed_monotone']} ({row['monotone_value_status']}) | "
            f"{row['complete_for_upper_envelope']} | {row['complete_for_distribution']} |"
            for row in epsilon_rows
        ),
        "",
        "## Future cost estimates (not executed)",
        "",
        f"- Targeted complex-case audit: {len(complex_rows)} cases; use the observed strict-tail evidence for a planning range of roughly {max(1, len(complex_rows)) * 30:.0f}-{max(1, len(complex_rows)) * 300:.0f} s wall time before any manual high-precision decisions.",
        f"- Full-K10 controls: {len(workflow.read_csv(output_dir / 'full_k10_control_manifest.csv')) if (output_dir / 'full_k10_control_manifest.csv').exists() else 0} prepared cases; no solve was launched.",
        "- Refinement: estimate only; defer sizing until complex cases and control results establish which epsilon intervals remain scientifically actionable.",
        "",
        "All plots in `main_pass_plots/` are explicitly provisional where the underlying distribution is incomplete.",
    ]
    workflow.atomic_write_text(output_dir / "main_pass_report.md", "\n".join(report_lines) + "\n")
    return {
        "root_calculations": 0,
        "manifest_count": len(case_rows),
        "resolved_count": len(resolved_case_rows),
        "status_counts": status_counts,
        "new_complex_count": len(new_complex_rows),
        "complex_count": len(complex_rows),
        "epsilon_rows": epsilon_rows,
        "argmax_rows": argmax_rows,
        "plot_paths": plot_paths,
        "integrity": integrity,
        "postprocess_seconds": time.perf_counter() - started,
    }


def analyze(output_dir: Path) -> dict[str, object]:
    output_dir = output_dir if output_dir.is_absolute() else REPO_ROOT / output_dir
    manifest_path = output_dir / "grid_manifest.csv"
    spectra_path = output_dir / "spectra_long.csv"
    if not manifest_path.exists() or not spectra_path.exists():
        raise FileNotFoundError("grid_manifest.csv and spectra_long.csv are required")
    manifest_rows = workflow.read_csv(manifest_path)
    spectra_rows = workflow.read_csv(spectra_path)
    execution_path = output_dir / "case_execution.csv"
    execution_rows = workflow.read_csv(execution_path) if execution_path.exists() else []
    control_result_path = output_dir / "full_k10_control_results.csv"
    control_rows = workflow.read_csv(control_result_path) if control_result_path.exists() else []
    comparison_path = output_dir / "prefix_full_comparison.csv"
    comparison_rows = workflow.read_csv(comparison_path) if comparison_path.exists() else []
    runtime_path = output_dir / "runtime_by_case.csv"
    runtime_rows = workflow.read_csv(runtime_path) if runtime_path.exists() else []
    metadata: Mapping[str, object] = {}
    metadata_path = output_dir / "run_metadata.json"
    if metadata_path.exists():
        try:
            metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            metadata = {}
    run_scope = metadata.get("run_scope", {}) if isinstance(metadata, Mapping) else {}
    controls_required = bool(
        isinstance(run_scope, Mapping) and run_scope.get("full_k10_controls_required")
    )
    case_rows = build_case_summary(manifest_rows, spectra_rows, execution_rows, control_rows)
    epsilon_rows = build_epsilon_summary(
        case_rows,
        comparison_rows,
        controls_required=controls_required,
    )
    argmax_rows = build_argmax_cases(case_rows, epsilon_rows)
    regression_rows = build_regression_checks(case_rows)
    unresolved_rows = build_unresolved(case_rows)
    not_attempted_rows = build_not_attempted(case_rows)
    refinement_rows = build_refinement_proposal(case_rows, epsilon_rows, runtime_rows)
    workflow.write_csv(output_dir / "case_summary.csv", case_rows, CASE_SUMMARY_FIELDS)
    workflow.write_csv(output_dir / "epsilon_summary.csv", epsilon_rows, EPSILON_SUMMARY_FIELDS)
    workflow.write_csv(output_dir / "argmax_cases.csv", argmax_rows, ARGMAX_FIELDS)
    workflow.write_csv(output_dir / "regression_checks.csv", regression_rows, REGRESSION_FIELDS)
    workflow.write_csv(output_dir / "unresolved_cases.csv", unresolved_rows, UNRESOLVED_FIELDS)
    workflow.write_csv(output_dir / "not_attempted_cases.csv", not_attempted_rows, NOT_ATTEMPTED_FIELDS)
    workflow.write_csv(output_dir / "refinement_proposal.csv", refinement_rows, REFINEMENT_FIELDS)
    plot_paths = _plot_diagnostics(output_dir, epsilon_rows, case_rows, runtime_rows)
    report_path = write_report(
        output_dir,
        case_rows,
        epsilon_rows,
        regression_rows,
        unresolved_rows,
        not_attempted_rows,
        refinement_rows,
        comparison_rows,
    )
    return {
        "root_calculations": 0,
        "manifest_count": len(manifest_rows),
        "spectra_row_count": len(spectra_rows),
        "resolved_count": sum(row["quality_status"] == "resolved" for row in case_rows),
        "unresolved_count": len(unresolved_rows),
        "not_attempted_count": len(not_attempted_rows),
        "regression_pass_count": sum(row["regression_status"] == "pass" for row in regression_rows),
        "full_control_count": len(comparison_rows),
        "full_control_pass_count": sum(row.get("comparison_status") == "PASS" for row in comparison_rows),
        "plot_paths": plot_paths,
        "report_path": report_path,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = (
        analyze_partial_cache(args.output_dir)
        if args.partial_cache
        else analyze(args.output_dir)
    )
    if args.partial_cache:
        print(
            f"partial-cache zero-solve: manifest={result['manifest_count']} "
            f"completed={result['completed_count']} exact={result['exact_N_true_count']} "
            f"unresolved={result['genuinely_unresolved_count']} "
            f"interrupted={result['interrupted_incomplete_count']} "
            f"not_attempted={result['not_attempted_count']} root_calculations=0"
        )
        return result
    print(
        f"postprocess-only: manifest={result['manifest_count']} spectra_rows={result['spectra_row_count']} "
        f"resolved={result['resolved_count']} unresolved={result['unresolved_count']} root_calculations=0"
    )
    return result


if __name__ == "__main__":
    main()
