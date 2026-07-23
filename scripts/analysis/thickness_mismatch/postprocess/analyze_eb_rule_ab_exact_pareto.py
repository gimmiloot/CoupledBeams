from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass, field
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_safe_prefix_certification as certification,
)


ALGORITHM_VERSION = "eb_rule_ab_exact_pareto_v2"
K_MAX = 10
FREQUENCY_ERROR_THRESHOLD = 0.10
GEOMETRY_TOLERANCE = 1.0e-12
DEFAULT_N_SHAPE_POINTS = 401
DEFAULT_DEVELOPMENT_DIR = REPO_ROOT / "results" / "eb_epsilon_apriori_pilot_branch_continuation_v1"
DEFAULT_PILOT_METADATA_DIR = REPO_ROOT / "results" / "eb_epsilon_apriori_pilot"
DEFAULT_STEP3A_DIR = REPO_ROOT / "results" / "eb_epsilon_lower_envelope_step3a"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "eb_rule_ab_exact_pareto"
LOCKED_HOLDOUT_IDS = ("S3_12", "S3_14")

PARTITION_DEVELOPMENT = "development_in_sample"
PARTITION_BASELINE = "step3a_baseline_controls"
PARTITION_DIRECTED = "step3a_nonbaseline_directed_validation"
PARTITION_HOLDOUT = "locked_adversarial_holdout"
PARTITION_EXCLUDED = "excluded"

NUMERICAL_OUTPUT_NAMES = (
    "data_partition_audit.csv",
    "included_geometry_audit.csv",
    "excluded_geometry_audit.csv",
    "prefix_predictor_metrics.csv",
    "rule_A_exact_search.csv",
    "rule_A_predictions.csv",
    "rule_B_exact_search_summary.csv",
    "rule_B_optimal_threshold_pairs.csv",
    "rule_B_pareto_frontier.csv",
    "rule_B_dominance_witnesses.csv",
    "rule_B_predictions.csv",
    "rule_B_equal_optimum_prediction_equivalence.csv",
    "rule_S_exact_search.csv",
    "rule_S_predictions.csv",
    "rule_S_validation_summary.csv",
    "rule_S_equivalence_audit.csv",
    "rule_S_rule_B_equivalence_audit.csv",
    "rule_S_threshold_margin_audit.csv",
    "rule_S_leave_one_geometry_out.csv",
    "rule_S_predictor_implementation_equivalence.csv",
    "partition_validation_summary.csv",
    "threshold_pair_validation_sensitivity.csv",
    "operation_counts.csv",
)
PLOT_NAMES = (
    "rule_A_prefix_points.png",
    "rule_B_prefix_plane.png",
    "partition_N_hat_vs_N_true.png",
)
REPORT_NAME = "eb_rule_ab_exact_pareto_report.md"
COST_PROPOSAL_NAME = "rule_B_cost_break_even_proposal.csv"


PARTITION_AUDIT_FIELDS = (
    "source",
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "canonical_geometry_key",
    "assigned_partition",
    "duplicate_status",
    "exclusion_reason",
    "quality_status",
    "root_source",
    "predictor_source",
    "algorithm_version",
)
PREFIX_FIELDS = (
    "source",
    "case_id",
    "assigned_partition",
    "canonical_geometry_key",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "prefix_n",
    "N_true",
    "prefix_label",
    "P_prefix_max_Pi_EB",
    "S_prefix_max_Pi_shear",
    "R_prefix_max_Pi_rotary",
    "Pi_EB_mode",
    "Pi_shear_mode",
    "Pi_rotary_mode",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_f",
    "root_source",
    "predictor_source",
)
PREDICTION_FIELDS = (
    "rule",
    "source",
    "case_id",
    "assigned_partition",
    "canonical_geometry_key",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "N_true",
    "N_hat",
    "false_safe",
    "false_safe_frequency_count",
    "conservative_loss",
    "exact_match",
    "first_rejected_mode",
    "trigger_component",
    "T_A",
    "T_s",
    "T_r",
    "threshold_pair_id",
    "threshold_provenance",
)
SUMMARY_FIELDS = (
    "summary_partition",
    "rule",
    "geometry_count",
    "false_safe_geometry_count",
    "false_safe_frequency_count",
    "worst_overprediction",
    "mean_N_true",
    "mean_N_hat",
    "exact_match_rate",
    "mean_conservative_loss",
    "maximum_conservative_loss",
    "usable_frequency_retention",
    "N_hat_zero_count",
    "N_hat_ten_count",
    "selected_thresholds",
    "threshold_provenance",
)
OPERATION_FIELDS = (
    "scope",
    "partition",
    "rule",
    "case_id",
    "existing_CSV_rows_read",
    "existing_root_records_consumed",
    "EB_mode_reconstructions",
    "EB_characteristic_matrix_evaluations",
    "SVD_6x6_calls",
    "shape_grid_point_evaluations",
    "quadrature_calls",
    "quadrature_point_evaluations",
    "rule_A_threshold_candidates",
    "rule_B_threshold_pairs_evaluated",
    "rule_S_threshold_candidates",
    "scalar_comparisons",
    "dominance_comparisons",
    "predictor_reconstruction_cache_hits",
    "specialized_shear_predictor_calls",
    "specialized_shear_quadrature_calls",
    "specialized_shear_quadrature_point_evaluations",
    "estimated_online_mode_reconstructions",
    "estimated_online_SVD_6x6_calls",
    "estimated_online_quadrature_calls",
    "estimated_online_scalar_comparisons",
    "Timoshenko_suffix_frequencies",
    "input_files",
    "notes",
)


@dataclass(frozen=True)
class Args:
    development_dir: Path
    pilot_metadata_dir: Path
    step3a_dir: Path
    output_dir: Path
    n_shape_points: int
    force: bool
    plot_only: bool


@dataclass
class GeometryData:
    source: str
    case_id: str
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float
    modes: list[dict[str, object]]
    N_true: int
    quality_status: str
    root_source: str
    predictor_source: str
    assigned_partition: str = ""
    duplicate_status: str = "unique"
    exclusion_reason: str = ""

    @property
    def key(self) -> tuple[int, int, int, int]:
        return canonical_geometry_key(self.epsilon_0, self.beta_deg, self.mu, self.eta)


@dataclass
class OfflineCounts:
    existing_csv_rows_read: int = 0
    existing_root_records_consumed: int = 0
    input_files: list[str] = field(default_factory=list)
    reconstruction: certification.OperationCounts = field(default_factory=certification.OperationCounts)
    rule_A_threshold_candidates: int = 0
    rule_B_threshold_pairs_evaluated: int = 0
    rule_S_threshold_candidates: int = 0
    scalar_comparisons: int = 0
    dominance_comparisons: int = 0
    predictor_reconstruction_cache_hits: int = 0
    specialized_shear_predictor_calls: int = 0
    specialized_shear_quadrature_calls: int = 0
    specialized_shear_quadrature_point_evaluations: int = 0


@dataclass(frozen=True)
class RuleASearchResult:
    candidates: tuple[float, ...]
    selected_threshold: float
    objective: int
    rows: tuple[dict[str, object], ...]
    scalar_comparisons: int


@dataclass(frozen=True)
class RuleBSearchResult:
    shear_candidates: tuple[float, ...]
    rotary_candidates: tuple[float, ...]
    evaluated_pair_count: int
    rejected_false_safe_count: int
    objective: int
    all_optimal_pairs: tuple[tuple[float, float], ...]
    pareto_optimal_pairs: tuple[tuple[float, float], ...]
    representative: tuple[float, float]
    scalar_comparisons: int


@dataclass(frozen=True)
class RuleSSearchResult:
    candidates: tuple[float, ...]
    selected_threshold: float
    objective: int
    rows: tuple[dict[str, object], ...]
    scalar_comparisons: int


class DataQualityError(RuntimeError):
    pass


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Run the exact CSV-first Rule-A/Rule-B safe-prefix experiment. "
            "Saved roots are consumed; no eigenvalue problem is solved."
        ),
    )
    parser.add_argument("--development-dir", type=Path, default=DEFAULT_DEVELOPMENT_DIR)
    parser.add_argument("--pilot-metadata-dir", type=Path, default=DEFAULT_PILOT_METADATA_DIR)
    parser.add_argument("--step3a-dir", type=Path, default=DEFAULT_STEP3A_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_N_SHAPE_POINTS)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if int(ns.n_shape_points) < 51:
        raise ValueError("--n-shape-points must be at least 51")
    return Args(
        development_dir=repo_path(Path(ns.development_dir)),
        pilot_metadata_dir=repo_path(Path(ns.pilot_metadata_dir)),
        step3a_dir=repo_path(Path(ns.step3a_dir)),
        output_dir=repo_path(Path(ns.output_dir)),
        n_shape_points=int(ns.n_shape_points),
        force=bool(ns.force),
        plot_only=bool(ns.plot_only),
    )


def fmt(value: object) -> object:
    return certification.fmt(value)


def read_csv_counted(path: Path, counts: OfflineCounts) -> list[dict[str, str]]:
    rows = certification.read_csv(path)
    counts.existing_csv_rows_read += len(rows)
    counts.input_files.append(rel(path))
    return rows


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def finite_float(row: Mapping[str, object], key: str) -> float:
    return certification.finite_float(row, key)


def bool_value(value: object) -> bool:
    return certification.bool_value(value)


def canonical_geometry_key(
    epsilon_0: float,
    beta_deg: float,
    mu: float,
    eta: float,
    *,
    tolerance: float = GEOMETRY_TOLERANCE,
) -> tuple[int, int, int, int]:
    values = (epsilon_0, beta_deg, mu, eta)
    if any(not math.isfinite(float(value)) for value in values):
        raise ValueError("geometry coordinates must be finite")
    return tuple(int(round(float(value) / float(tolerance))) for value in values)  # type: ignore[return-value]


def canonical_geometry_key_text(key: Sequence[int]) -> str:
    return "|".join(str(int(value)) for value in key)


def same_geometry(
    left: GeometryData,
    right: GeometryData,
    *,
    tolerance: float = GEOMETRY_TOLERANCE,
) -> bool:
    return all(
        abs(float(left_value) - float(right_value)) <= float(tolerance)
        for left_value, right_value in zip(
            (left.epsilon_0, left.beta_deg, left.mu, left.eta),
            (right.epsilon_0, right.beta_deg, right.mu, right.eta),
            strict=True,
        )
    )


def group_by_case(rows: Iterable[Mapping[str, object]]) -> dict[str, list[dict[str, object]]]:
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get("case_id", ""))].append(dict(row))
    return grouped


def exact_mode_rows(rows: Sequence[Mapping[str, object]], index_field: str = "sorted_index") -> list[dict[str, object]]:
    selected = sorted(
        (dict(row) for row in rows if 1 <= int(float(row.get(index_field, 0) or 0)) <= K_MAX),
        key=lambda row: int(float(row[index_field])),
    )
    indices = [int(float(row[index_field])) for row in selected]
    if indices != list(range(1, K_MAX + 1)):
        raise DataQualityError(f"expected exactly sorted modes 1..10, observed {indices}")
    return selected


def validate_true_prefix(modes: Sequence[Mapping[str, object]], expected: int) -> int:
    deltas = [finite_float(row, "delta_f") for row in modes]
    observed = certification.true_safe_prefix(
        deltas,
        threshold=FREQUENCY_ERROR_THRESHOLD,
        k_max=K_MAX,
        strict=True,
    )
    if observed != int(expected):
        raise DataQualityError(f"stored N_true={expected} disagrees with saved roots/deltas: {observed}")
    return observed


def load_development(args: Args, counts: OfflineCounts) -> list[GeometryData]:
    manifest_rows = read_csv_counted(args.development_dir / "epsilon_pilot_case_manifest_resolved.csv", counts)
    geometry_rows = read_csv_counted(args.development_dir / "epsilon_pilot_geometry_metrics.csv", counts)
    mode_rows = read_csv_counted(args.development_dir / "epsilon_pilot_mode_metrics.csv", counts)
    manifest = {row["case_id"]: row for row in manifest_rows}
    geometry = {row["case_id"]: row for row in geometry_rows}
    modes_by_case = group_by_case(mode_rows)
    if set(manifest) != set(geometry) or len(manifest) != 21:
        raise DataQualityError("development set must contain the same 21 manifest and geometry case IDs")
    records: list[GeometryData] = []
    for case_id in sorted(manifest):
        meta = manifest[case_id]
        summary = geometry[case_id]
        if str(meta.get("spectrum_method", "")) != "branch_informed_continuation_v1":
            raise DataQualityError(f"{case_id}: branch-informed roots are required; legacy substitution is forbidden")
        required = (
            str(summary.get("quality_status", "")) == "included"
            and bool_value(summary.get("complete_K10", False))
            and bool_value(summary.get("root11_guard_available", False))
            and str(summary.get("spectrum_method", "")) == "branch_informed_continuation_v1"
            and str(summary.get("EB_spectrum_status", "")) in {"K10_guard_resolved", "resolved_complete"}
            and str(summary.get("Timo_spectrum_status", "")) in {"K10_guard_resolved", "resolved_complete"}
            and int(float(summary.get("root_warning_count", 0) or 0)) == 0
            and int(float(summary.get("candidate_boundary_count", 0) or 0)) == 0
            and not str(summary.get("exclusion_reason", "")).strip()
        )
        if not required:
            raise DataQualityError(f"{case_id}: development K10 quality contract is not satisfied")
        selected = exact_mode_rows(modes_by_case.get(case_id, ()))
        normalized_modes: list[dict[str, object]] = []
        for row in selected:
            values = {
                "sorted_index": int(float(row["sorted_index"])),
                "Lambda_EB": finite_float(row, "Lambda_EB"),
                "Lambda_Timo": finite_float(row, "Lambda_Timo"),
                "delta_f": certification.squared_lambda_delta(
                    finite_float(row, "Lambda_EB"), finite_float(row, "Lambda_Timo")
                ),
                "Pi_shear_EB": finite_float(row, "Pi_shear_EB"),
                "Pi_rotary_EB": finite_float(row, "Pi_rotary_EB"),
                "Pi_EB": finite_float(row, "Pi_EB"),
            }
            if not all(math.isfinite(float(value)) for key, value in values.items() if key != "sorted_index"):
                raise DataQualityError(f"{case_id}: non-finite saved frequency or Pi value")
            if not math.isclose(
                float(values["Pi_EB"]),
                float(values["Pi_shear_EB"]) + float(values["Pi_rotary_EB"]),
                rel_tol=1.0e-12,
                abs_tol=1.0e-14,
            ):
                raise DataQualityError(f"{case_id}: stored Pi_EB is not the saved component sum")
            if not math.isclose(
                float(values["delta_f"]),
                finite_float(row, "delta_f"),
                rel_tol=1.0e-8,
                abs_tol=1.0e-10,
            ):
                raise DataQualityError(f"{case_id}: saved delta_f disagrees with squared Lambda values")
            normalized_modes.append(values)
        n_true = validate_true_prefix(normalized_modes, int(float(summary["N_true"])))
        records.append(
            GeometryData(
                source="branch_informed_pilot",
                case_id=case_id,
                epsilon_0=finite_float(meta, "epsilon_0"),
                beta_deg=finite_float(meta, "beta_deg"),
                mu=finite_float(meta, "mu"),
                eta=finite_float(meta, "eta"),
                modes=normalized_modes,
                N_true=n_true,
                quality_status="quality_approved",
                root_source="branch_informed_saved_roots",
                predictor_source="stored_branch_informed_pilot_Pi",
                assigned_partition=PARTITION_DEVELOPMENT,
            )
        )
        counts.existing_root_records_consumed += 2 * K_MAX
    if any(same_geometry(left, right) for index, left in enumerate(records) for right in records[index + 1 :]):
        raise DataQualityError("development set contains duplicate geometries at tolerance 1e-12")
    return records


def strict_quality_by_case(rows: Sequence[Mapping[str, object]]) -> dict[str, bool]:
    grouped = group_by_case(rows)
    quality: dict[str, bool] = {}
    for case_id, items in grouped.items():
        indices = {int(float(row["prefix_n"])) for row in items}
        quality[case_id] = indices == set(range(1, 11)) and all(
            str(row.get("verification_status", "")) == "pass"
            and bool_value(row.get("primary_K10_resolved", False))
            and bool_value(row.get("verification_K10_resolved", False))
            and bool_value(row.get("roots_agree", False))
            and bool_value(row.get("cluster_agreement", False))
            for row in items
        )
    return quality


def load_step3a(args: Args, counts: OfflineCounts) -> list[GeometryData]:
    manifest_rows = read_csv_counted(args.step3a_dir / "step3a_manifest_resolved.csv", counts)
    case_rows = read_csv_counted(args.step3a_dir / "step3a_case_summary.csv", counts)
    spectrum_rows = read_csv_counted(args.step3a_dir / "step3a_case_spectrum_summary.csv", counts)
    mode_rows = read_csv_counted(args.step3a_dir / "step3a_mode_metrics.csv", counts)
    strict_rows = read_csv_counted(args.step3a_dir / "step3a_strict_verification_audit.csv", counts)
    baseline_rows = read_csv_counted(args.step3a_dir / "step3a_baseline_control_audit.csv", counts)
    manifest = {row["case_id"]: row for row in manifest_rows}
    cases = {row["case_id"]: row for row in case_rows}
    modes_by_case = group_by_case(mode_rows)
    spectra_by_case = group_by_case(spectrum_rows)
    strict_quality = strict_quality_by_case(strict_rows)
    if set(manifest) != set(cases) or len(cases) != 28:
        raise DataQualityError("Step-3A must contain the same 28 manifest and case-summary IDs")
    if not baseline_rows or any(str(row.get("baseline_control_status", "")) != "pass" for row in baseline_rows):
        raise DataQualityError("Step-3A baseline-control oracle audit is incomplete or failed")
    records: list[GeometryData] = []
    for case_id in sorted(cases):
        meta = manifest[case_id]
        summary = cases[case_id]
        spectra = spectra_by_case.get(case_id, [])
        models = {str(row.get("model", "")): row for row in spectra}
        spectrum_ok = len(models) == 2 and all(
            str(row.get("algorithm_version", "")) == "branch_informed_continuation_v1"
            and str(row.get("run_scope", "")) == "primary"
            and str(row.get("roots_1_11_status", "")) == "resolved"
            and bool_value(row.get("K10_guard_resolved", False))
            and bool_value(row.get("guard_passed", False))
            and not str(row.get("warnings", "")).strip()
            for row in models.values()
        )
        if (
            str(summary.get("K10_status", "")) != "included"
            or not spectrum_ok
            or str(meta.get("manifest_validation_status", "")) != "pass"
            or str(meta.get("required_spectrum_method", "")) != "branch_informed_continuation_v1"
            or not bool_value(meta.get("required_K10_guard", False))
        ):
            raise DataQualityError(f"{case_id}: Step-3A primary K10 quality contract is not satisfied")
        verification_required = bool_value(summary.get("strict_verification_triggered", False))
        if verification_required and not strict_quality.get(case_id, False):
            raise DataQualityError(f"{case_id}: required independent verification did not pass")
        selected = exact_mode_rows(modes_by_case.get(case_id, ()))
        normalized_modes: list[dict[str, object]] = []
        root_source = "independently_verified_step3a_roots" if verification_required else "quality_approved_step3a_primary_roots"
        for row in selected:
            if str(row.get("EB_root_quality_status", "")) != "pass" or str(row.get("Timo_root_quality_status", "")) != "pass":
                raise DataQualityError(f"{case_id}: primary mode residual quality failed")
            if verification_required:
                if str(row.get("verification_residual_status", "")) != "pass":
                    raise DataQualityError(f"{case_id}: verification mode residual quality failed")
                lambda_eb = finite_float(row, "Lambda_EB_verification")
                lambda_timo = finite_float(row, "Lambda_Timo_verification")
                stored_delta = finite_float(row, "delta_f_verification")
            else:
                lambda_eb = finite_float(row, "Lambda_EB_primary")
                lambda_timo = finite_float(row, "Lambda_Timo_primary")
                stored_delta = finite_float(row, "delta_f_primary")
            delta = certification.squared_lambda_delta(lambda_eb, lambda_timo)
            if not all(math.isfinite(value) for value in (lambda_eb, lambda_timo, delta, stored_delta)):
                raise DataQualityError(f"{case_id}: selected existing roots are missing or non-finite")
            if not math.isclose(delta, stored_delta, rel_tol=1.0e-8, abs_tol=1.0e-10):
                raise DataQualityError(f"{case_id}: selected delta disagrees with squared roots")
            normalized_modes.append(
                {
                    "sorted_index": int(float(row["sorted_index"])),
                    "Lambda_EB": lambda_eb,
                    "Lambda_Timo": lambda_timo,
                    "delta_f": delta,
                }
            )
        n_true = validate_true_prefix(normalized_modes, int(float(summary["N_true"])))
        records.append(
            GeometryData(
                source="step3a",
                case_id=case_id,
                epsilon_0=finite_float(meta, "epsilon"),
                beta_deg=finite_float(meta, "beta_deg"),
                mu=finite_float(meta, "mu"),
                eta=finite_float(meta, "eta"),
                modes=normalized_modes,
                N_true=n_true,
                quality_status="quality_approved",
                root_source=root_source,
                predictor_source="reconstructed_from_selected_saved_EB_root",
            )
        )
        counts.existing_root_records_consumed += 2 * K_MAX
    return records


def assign_partitions(development: Sequence[GeometryData], step3a: Sequence[GeometryData]) -> None:
    if any(record.case_id in LOCKED_HOLDOUT_IDS for record in development):
        raise AssertionError("S3_12/S3_14 must not occur in development")
    for record in step3a:
        case_group = "baseline_control" if record.beta_deg == 0.0 and record.mu == 0.0 and record.eta == 0.0 else "nonbaseline"
        if record.case_id in LOCKED_HOLDOUT_IDS:
            record.assigned_partition = PARTITION_HOLDOUT
        elif any(same_geometry(record, development_record) for development_record in development):
            record.assigned_partition = PARTITION_EXCLUDED
            record.duplicate_status = "duplicate_of_development"
            record.exclusion_reason = "exact_geometry_duplicate_of_development_within_1e-12"
        elif case_group == "baseline_control":
            record.assigned_partition = PARTITION_BASELINE
        else:
            record.assigned_partition = PARTITION_DIRECTED

    directed = [record for record in step3a if record.assigned_partition in {PARTITION_BASELINE, PARTITION_DIRECTED}]
    holdout = [record for record in step3a if record.assigned_partition == PARTITION_HOLDOUT]
    if {record.case_id for record in holdout} != set(LOCKED_HOLDOUT_IDS):
        raise AssertionError("locked adversarial holdout must contain exactly S3_12 and S3_14")
    if any(record.case_id in LOCKED_HOLDOUT_IDS for record in directed):
        raise AssertionError("S3_12/S3_14 must not occur in directed validation")
    if any(same_geometry(development_record, validation_record) for development_record in development for validation_record in directed):
        raise AssertionError("development and directed validation contain a geometry duplicate")
    included = [*development, *directed, *holdout]
    if any(same_geometry(left, right) for index, left in enumerate(included) for right in included[index + 1 :]):
        raise AssertionError("one geometry belongs to more than one included partition")


def reconstruct_step3a_predictors(
    records: Sequence[GeometryData],
    *,
    n_shape_points: int,
    counts: OfflineCounts,
) -> list[dict[str, object]]:
    targets = [record for record in records if record.assigned_partition != PARTITION_EXCLUDED]
    equivalence_rows: list[dict[str, object]] = []
    for record in targets:
        for mode in record.modes:
            try:
                result = certification.fixed_scan.eb_mode_result(
                    epsilon=record.epsilon_0,
                    beta_deg=record.beta_deg,
                    mu=record.mu,
                    eta=record.eta,
                    sorted_index=int(mode["sorted_index"]),
                    Lambda=finite_float(mode, "Lambda_EB"),
                    n_points=int(n_shape_points),
                    root_warnings=(),
                )
                predictors = certification.fixed_scan.eb_predictors(result)
                specialized = certification.fixed_scan.eb_shear_predictor(result)
            except (FloatingPointError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
                raise DataQualityError(
                    f"{record.case_id}: EB predictor reconstruction failed at mode "
                    f"{mode.get('sorted_index')}: {exc}"
                ) from exc
            counts.reconstruction.EB_mode_reconstructions += 1
            counts.reconstruction.EB_characteristic_matrix_evaluations += 1
            counts.reconstruction.SVD_6x6_calls += 1
            counts.reconstruction.EB_shape_grid_points += 2 * int(n_shape_points)
            counts.reconstruction.quadrature_calls += 14
            counts.reconstruction.quadrature_point_evaluations += 14 * int(n_shape_points)
            counts.specialized_shear_predictor_calls += 1
            counts.specialized_shear_quadrature_calls += 2
            counts.specialized_shear_quadrature_point_evaluations += 2 * int(n_shape_points)
            existing_shear = float(predictors["Pi_shear_EB"])
            specialized_shear = float(specialized["Pi_shear_EB"])
            absolute_difference = abs(existing_shear - specialized_shear)
            equivalent = math.isclose(
                existing_shear,
                specialized_shear,
                rel_tol=1.0e-13,
                abs_tol=1.0e-15,
            )
            equivalence_rows.append(
                {
                    "source": record.source,
                    "case_id": record.case_id,
                    "partition": record.assigned_partition,
                    "sorted_index": int(mode["sorted_index"]),
                    "Lambda_EB": finite_float(mode, "Lambda_EB"),
                    "Pi_shear_existing": existing_shear,
                    "Pi_shear_specialized": specialized_shear,
                    "absolute_difference": absolute_difference,
                    "relative_difference": absolute_difference / max(abs(existing_shear), 1.0e-300),
                    "relative_tolerance": 1.0e-13,
                    "absolute_tolerance": 1.0e-15,
                    "equivalent": equivalent,
                    "existing_predictor_quadrature_calls": 12,
                    "specialized_rule_S_total_quadrature_calls": 6,
                    "specialized_rule_S_shear_quadrature_calls": 2,
                }
            )
            if not equivalent:
                raise DataQualityError(
                    f"{record.case_id}: specialized Pi_shear mismatch at mode {mode.get('sorted_index')}"
                )
            mode.update(
                {
                    "Pi_shear_EB": predictors["Pi_shear_EB"],
                    "Pi_rotary_EB": predictors["Pi_rotary_EB"],
                    "Pi_EB": predictors["Pi_EB"],
                }
            )
            if not all(
                math.isfinite(finite_float(mode, field))
                for field in ("Pi_shear_EB", "Pi_rotary_EB", "Pi_EB")
            ):
                raise DataQualityError(f"{record.case_id}: reconstructed Pi is non-finite")
            if not math.isclose(
                finite_float(mode, "Pi_EB"),
                finite_float(mode, "Pi_shear_EB") + finite_float(mode, "Pi_rotary_EB"),
                rel_tol=1.0e-12,
                abs_tol=1.0e-14,
            ):
                raise DataQualityError(f"{record.case_id}: reconstructed Pi component sum failed")
    return equivalence_rows


def prefix_extrema(record: GeometryData) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    p_value = float("-inf")
    s_value = float("-inf")
    r_value = float("-inf")
    for mode in sorted(record.modes, key=lambda row: int(row["sorted_index"])):
        p_value = max(p_value, finite_float(mode, "Pi_EB"))
        s_value = max(s_value, finite_float(mode, "Pi_shear_EB"))
        r_value = max(r_value, finite_float(mode, "Pi_rotary_EB"))
        prefix_n = int(mode["sorted_index"])
        rows.append(
            {
                "source": record.source,
                "case_id": record.case_id,
                "assigned_partition": record.assigned_partition,
                "canonical_geometry_key": canonical_geometry_key_text(record.key),
                "epsilon_0": record.epsilon_0,
                "beta_deg": record.beta_deg,
                "mu": record.mu,
                "eta": record.eta,
                "prefix_n": prefix_n,
                "N_true": record.N_true,
                "prefix_label": "safe_prefix" if prefix_n <= record.N_true else "unsafe_prefix",
                "P_prefix_max_Pi_EB": p_value,
                "S_prefix_max_Pi_shear": s_value,
                "R_prefix_max_Pi_rotary": r_value,
                "Pi_EB_mode": finite_float(mode, "Pi_EB"),
                "Pi_shear_mode": finite_float(mode, "Pi_shear_EB"),
                "Pi_rotary_mode": finite_float(mode, "Pi_rotary_EB"),
                "Lambda_EB": finite_float(mode, "Lambda_EB"),
                "Lambda_Timo": finite_float(mode, "Lambda_Timo"),
                "delta_f": finite_float(mode, "delta_f"),
                "root_source": record.root_source,
                "predictor_source": record.predictor_source,
            }
        )
    return rows


def predict_rule_a(record: GeometryData, threshold: float) -> tuple[int, int, str]:
    comparisons = 0
    for row in prefix_extrema(record):
        comparisons += 1
        if float(row["P_prefix_max_Pi_EB"]) > float(threshold):
            return int(row["prefix_n"]) - 1, comparisons, "Pi_EB"
    return K_MAX, comparisons, ""


def predict_rule_b(record: GeometryData, threshold_s: float, threshold_r: float) -> tuple[int, int, str]:
    comparisons = 0
    for row in prefix_extrema(record):
        comparisons += 2
        shear_pass = float(row["S_prefix_max_Pi_shear"]) <= float(threshold_s)
        rotary_pass = float(row["R_prefix_max_Pi_rotary"]) <= float(threshold_r)
        if not (shear_pass and rotary_pass):
            return int(row["prefix_n"]) - 1, comparisons, "Pi_shear" if not shear_pass else "Pi_rotary"
    return K_MAX, comparisons, ""


def predict_rule_s(record: GeometryData, threshold_s: float) -> tuple[int, int, str]:
    comparisons = 0
    for row in prefix_extrema(record):
        comparisons += 1
        if float(row["S_prefix_max_Pi_shear"]) > float(threshold_s):
            return int(row["prefix_n"]) - 1, comparisons, "Pi_shear"
    return K_MAX, comparisons, ""


def rule_a_candidate_axis(records: Sequence[GeometryData]) -> tuple[float, ...]:
    values = {
        float(row["P_prefix_max_Pi_EB"])
        for record in records
        for row in prefix_extrema(record)
    }
    return (float("-inf"), *tuple(sorted(values)))


def rule_b_candidate_axes(records: Sequence[GeometryData]) -> tuple[tuple[float, ...], tuple[float, ...]]:
    points = [row for record in records for row in prefix_extrema(record)]
    shear = (float("-inf"), *tuple(sorted({float(row["S_prefix_max_Pi_shear"]) for row in points})))
    rotary = (float("-inf"), *tuple(sorted({float(row["R_prefix_max_Pi_rotary"]) for row in points})))
    return shear, rotary


def rule_s_candidate_axis(records: Sequence[GeometryData]) -> tuple[float, ...]:
    values = {
        float(row["S_prefix_max_Pi_shear"])
        for record in records
        for row in prefix_extrema(record)
    }
    return (float("-inf"), *tuple(sorted(values)))


def exact_rule_a_search(records: Sequence[GeometryData]) -> RuleASearchResult:
    candidates = rule_a_candidate_axis(records)
    rows: list[dict[str, object]] = []
    best_threshold = float("-inf")
    best_objective = -1
    comparisons = 0
    for threshold in candidates:
        predictions: list[int] = []
        for record in records:
            n_hat, used, _trigger = predict_rule_a(record, threshold)
            comparisons += used
            predictions.append(n_hat)
        false_safe = [max(n_hat - record.N_true, 0) for record, n_hat in zip(records, predictions, strict=True)]
        admissible = not any(false_safe)
        objective = sum(predictions)
        if admissible and (objective > best_objective or (objective == best_objective and threshold < best_threshold)):
            best_objective = objective
            best_threshold = threshold
        rows.append(
            {
                "candidate_index": len(rows) + 1,
                "T_A": threshold,
                "zero_development_false_safe": admissible,
                "false_safe_geometry_count": sum(value > 0 for value in false_safe),
                "false_safe_frequency_count": sum(false_safe),
                "retained_frequencies": objective,
                "conservative_loss": sum(max(record.N_true - n_hat, 0) for record, n_hat in zip(records, predictions, strict=True)),
                "selected": False,
                "search_is_exact": True,
                "threshold_provenance": "observed_development_prefix_extrema_plus_reject_all",
            }
        )
    for row in rows:
        row["selected"] = float(row["T_A"]) == best_threshold
    return RuleASearchResult(candidates, best_threshold, best_objective, tuple(rows), comparisons)


def pareto_nondominated_threshold_pairs(
    pairs: Iterable[tuple[float, float]],
) -> tuple[tuple[float, float], ...]:
    unique = sorted(set((float(left), float(right)) for left, right in pairs))
    frontier = [
        pair
        for pair in unique
        if not any(
            other != pair
            and other[0] >= pair[0]
            and other[1] >= pair[1]
            and (other[0] > pair[0] or other[1] > pair[1])
            for other in unique
        )
    ]
    return tuple(sorted(frontier))


def brute_force_rule_b_search(records: Sequence[GeometryData]) -> RuleBSearchResult:
    shear, rotary = rule_b_candidate_axes(records)
    evaluated = 0
    rejected = 0
    comparisons = 0
    best_objective = -1
    optimal: list[tuple[float, float]] = []
    for threshold_s in shear:
        for threshold_r in rotary:
            evaluated += 1
            objective = 0
            admissible = True
            for record in records:
                n_hat, used, _trigger = predict_rule_b(record, threshold_s, threshold_r)
                comparisons += used
                if n_hat > record.N_true:
                    admissible = False
                    break
                objective += n_hat
            if not admissible:
                rejected += 1
                continue
            if objective > best_objective:
                best_objective = objective
                optimal = [(threshold_s, threshold_r)]
            elif objective == best_objective:
                optimal.append((threshold_s, threshold_r))
    if not optimal:
        raise RuntimeError("Rule B exact search found no admissible candidate, including reject-all")
    frontier = pareto_nondominated_threshold_pairs(optimal)
    representative = min(optimal, key=lambda pair: (pair[0], pair[1]))
    return RuleBSearchResult(
        shear_candidates=shear,
        rotary_candidates=rotary,
        evaluated_pair_count=evaluated,
        rejected_false_safe_count=rejected,
        objective=best_objective,
        all_optimal_pairs=tuple(sorted(set(optimal))),
        pareto_optimal_pairs=frontier,
        representative=representative,
        scalar_comparisons=comparisons,
    )


def exact_rule_b_pareto_search(records: Sequence[GeometryData]) -> RuleBSearchResult:
    """Exact Cartesian search followed by exact threshold-space Pareto filtering."""

    return brute_force_rule_b_search(records)


def exact_rule_s_search(records: Sequence[GeometryData]) -> RuleSSearchResult:
    candidates = rule_s_candidate_axis(records)
    rows: list[dict[str, object]] = []
    best_threshold = float("-inf")
    best_objective = -1
    comparisons = 0
    for threshold in candidates:
        predictions: list[int] = []
        for record in records:
            n_hat, used, _trigger = predict_rule_s(record, threshold)
            comparisons += used
            predictions.append(n_hat)
        false_safe = [max(n_hat - record.N_true, 0) for record, n_hat in zip(records, predictions, strict=True)]
        admissible = not any(false_safe)
        objective = sum(predictions)
        if admissible and (objective > best_objective or (objective == best_objective and threshold < best_threshold)):
            best_objective = objective
            best_threshold = threshold
        rows.append(
            {
                "candidate_index": len(rows) + 1,
                "T_s": threshold,
                "zero_development_false_safe": admissible,
                "false_safe_geometry_count": sum(value > 0 for value in false_safe),
                "false_safe_frequency_count": sum(false_safe),
                "retained_frequencies": objective,
                "conservative_loss": sum(
                    max(record.N_true - n_hat, 0)
                    for record, n_hat in zip(records, predictions, strict=True)
                ),
                "selected": False,
                "search_is_exact": True,
                "threshold_provenance": "independent_rule_S_development_prefix_extrema_plus_reject_all",
            }
        )
    for row in rows:
        row["selected"] = float(row["T_s"]) == best_threshold
    return RuleSSearchResult(candidates, best_threshold, best_objective, tuple(rows), comparisons)


def prediction_row(
    record: GeometryData,
    *,
    rule: str,
    n_hat: int,
    trigger: str,
    threshold_a: float | str = "",
    threshold_s: float | str = "",
    threshold_r: float | str = "",
    pair_id: str = "",
) -> dict[str, object]:
    return {
        "rule": rule,
        "source": record.source,
        "case_id": record.case_id,
        "assigned_partition": record.assigned_partition,
        "canonical_geometry_key": canonical_geometry_key_text(record.key),
        "epsilon_0": record.epsilon_0,
        "beta_deg": record.beta_deg,
        "mu": record.mu,
        "eta": record.eta,
        "N_true": record.N_true,
        "N_hat": n_hat,
        "false_safe": n_hat > record.N_true,
        "false_safe_frequency_count": max(n_hat - record.N_true, 0),
        "conservative_loss": max(record.N_true - n_hat, 0),
        "exact_match": n_hat == record.N_true,
        "first_rejected_mode": n_hat + 1 if n_hat < K_MAX else "",
        "trigger_component": trigger,
        "T_A": threshold_a,
        "T_s": threshold_s,
        "T_r": threshold_r,
        "threshold_pair_id": pair_id,
        "threshold_provenance": "exact_development_prefix_extrema_search_frozen_before_validation",
    }


def predictions_for_rule_a(
    records: Sequence[GeometryData],
    threshold: float,
    *,
    counts: OfflineCounts | None = None,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for record in records:
        n_hat, comparisons, trigger = predict_rule_a(record, threshold)
        if counts is not None:
            counts.scalar_comparisons += comparisons
        output.append(prediction_row(record, rule="Rule_A", n_hat=n_hat, trigger=trigger, threshold_a=threshold))
    return output


def predictions_for_rule_b(
    records: Sequence[GeometryData],
    pair: tuple[float, float],
    *,
    pair_id: str,
    counts: OfflineCounts | None = None,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for record in records:
        n_hat, comparisons, trigger = predict_rule_b(record, *pair)
        if counts is not None:
            counts.scalar_comparisons += comparisons
        output.append(
            prediction_row(
                record,
                rule="Rule_B",
                n_hat=n_hat,
                trigger=trigger,
                threshold_s=pair[0],
                threshold_r=pair[1],
                pair_id=pair_id,
            )
        )
    return output


def predictions_for_rule_s(
    records: Sequence[GeometryData],
    threshold_s: float,
    *,
    counts: OfflineCounts | None = None,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for record in records:
        n_hat, comparisons, trigger = predict_rule_s(record, threshold_s)
        if counts is not None:
            counts.scalar_comparisons += comparisons
        output.append(
            prediction_row(
                record,
                rule="Rule_S",
                n_hat=n_hat,
                trigger=trigger,
                threshold_s=threshold_s,
            )
        )
    return output


def metrics_for_predictions(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    if not rows:
        return {
            "geometry_count": 0,
            "false_safe_geometry_count": 0,
            "false_safe_frequency_count": 0,
            "worst_overprediction": 0,
            "mean_N_true": float("nan"),
            "mean_N_hat": float("nan"),
            "exact_match_rate": float("nan"),
            "mean_conservative_loss": float("nan"),
            "maximum_conservative_loss": 0,
            "usable_frequency_retention": float("nan"),
            "N_hat_zero_count": 0,
            "N_hat_ten_count": 0,
        }
    n_true = np.asarray([int(float(row["N_true"])) for row in rows], dtype=int)
    n_hat = np.asarray([int(float(row["N_hat"])) for row in rows], dtype=int)
    over = np.maximum(n_hat - n_true, 0)
    loss = np.maximum(n_true - n_hat, 0)
    usable = int(np.sum(n_true))
    retained = int(np.sum(np.minimum(n_true, n_hat)))
    return {
        "geometry_count": len(rows),
        "false_safe_geometry_count": int(np.count_nonzero(over)),
        "false_safe_frequency_count": int(np.sum(over)),
        "worst_overprediction": int(np.max(over)),
        "mean_N_true": float(np.mean(n_true)),
        "mean_N_hat": float(np.mean(n_hat)),
        "exact_match_rate": float(np.mean(n_true == n_hat)),
        "mean_conservative_loss": float(np.mean(loss)),
        "maximum_conservative_loss": int(np.max(loss)),
        "usable_frequency_retention": float(retained / usable) if usable else float("nan"),
        "N_hat_zero_count": int(np.count_nonzero(n_hat == 0)),
        "N_hat_ten_count": int(np.count_nonzero(n_hat == K_MAX)),
    }


def summary_partition_rows(
    rule_a_predictions: Sequence[dict[str, object]],
    rule_b_predictions: Sequence[dict[str, object]],
    *,
    threshold_a: float,
    threshold_b: tuple[float, float],
) -> list[dict[str, object]]:
    definitions = (
        (PARTITION_DEVELOPMENT, lambda row: row["assigned_partition"] == PARTITION_DEVELOPMENT),
        (PARTITION_BASELINE, lambda row: row["assigned_partition"] == PARTITION_BASELINE),
        (PARTITION_DIRECTED, lambda row: row["assigned_partition"] == PARTITION_DIRECTED),
        ("S3_12", lambda row: row["case_id"] == "S3_12"),
        ("S3_14", lambda row: row["case_id"] == "S3_14"),
        ("S3_combined", lambda row: row["case_id"] in LOCKED_HOLDOUT_IDS),
    )
    output: list[dict[str, object]] = []
    for rule, predictions in (("Rule_A", rule_a_predictions), ("Rule_B", rule_b_predictions)):
        for partition, predicate in definitions:
            selected = [row for row in predictions if predicate(row)]
            thresholds = {"T_A": threshold_a} if rule == "Rule_A" else {"T_s": threshold_b[0], "T_r": threshold_b[1]}
            output.append(
                {
                    "summary_partition": partition,
                    "rule": rule,
                    **metrics_for_predictions(selected),
                    "selected_thresholds": thresholds,
                    "threshold_provenance": "development_only_exact_search",
                }
            )
    return output


def rule_s_summary_rows(
    predictions: Sequence[dict[str, object]],
    *,
    threshold_s: float,
) -> list[dict[str, object]]:
    definitions = (
        (PARTITION_DEVELOPMENT, lambda row: row["assigned_partition"] == PARTITION_DEVELOPMENT),
        (PARTITION_BASELINE, lambda row: row["assigned_partition"] == PARTITION_BASELINE),
        (PARTITION_DIRECTED, lambda row: row["assigned_partition"] == PARTITION_DIRECTED),
        ("S3_12", lambda row: row["case_id"] == "S3_12"),
        ("S3_14", lambda row: row["case_id"] == "S3_14"),
        ("S3_combined", lambda row: row["case_id"] in LOCKED_HOLDOUT_IDS),
    )
    return [
        {
            "summary_partition": partition,
            "rule": "Rule_S",
            **metrics_for_predictions([row for row in predictions if predicate(row)]),
            "selected_thresholds": {"T_s": threshold_s},
            "threshold_provenance": "independent_development_only_exact_rule_S_search",
        }
        for partition, predicate in definitions
    ]


def dominance_audit(
    development_prefix_rows: Sequence[Mapping[str, object]],
) -> tuple[list[dict[str, object]], dict[str, object]]:
    safe = [row for row in development_prefix_rows if row["prefix_label"] == "safe_prefix"]
    unsafe = [row for row in development_prefix_rows if row["prefix_label"] == "unsafe_prefix"]
    witnesses: list[dict[str, object]] = []
    comparisons = 0
    for safe_row in safe:
        candidates: list[Mapping[str, object]] = []
        for unsafe_row in unsafe:
            comparisons += 1
            if (
                float(unsafe_row["S_prefix_max_Pi_shear"]) <= float(safe_row["S_prefix_max_Pi_shear"])
                and float(unsafe_row["R_prefix_max_Pi_rotary"]) <= float(safe_row["R_prefix_max_Pi_rotary"])
            ):
                candidates.append(unsafe_row)
        if candidates:
            witness = min(
                candidates,
                key=lambda row: (
                    float(row["S_prefix_max_Pi_shear"]) + float(row["R_prefix_max_Pi_rotary"]),
                    str(row["case_id"]),
                    int(row["prefix_n"]),
                ),
            )
            witnesses.append(
                {
                    "safe_case_id": safe_row["case_id"],
                    "safe_prefix_n": safe_row["prefix_n"],
                    "safe_S": safe_row["S_prefix_max_Pi_shear"],
                    "safe_R": safe_row["R_prefix_max_Pi_rotary"],
                    "unsafe_case_id": witness["case_id"],
                    "unsafe_prefix_n": witness["prefix_n"],
                    "unsafe_S": witness["S_prefix_max_Pi_shear"],
                    "unsafe_R": witness["R_prefix_max_Pi_rotary"],
                    "dominance_relation": "unsafe_S<=safe_S_and_unsafe_R<=safe_R",
                }
            )
    summary = {
        "safe_prefix_point_count": len(safe),
        "unsafe_prefix_point_count": len(unsafe),
        "unavoidably_rejected_safe_prefix_count": len(witnesses),
        "unavoidably_rejected_safe_prefix_fraction": len(witnesses) / len(safe) if safe else float("nan"),
        "dominance_comparisons": comparisons,
    }
    return witnesses, summary


def threshold_pair_rows(result: RuleBSearchResult) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    frontier = set(result.pareto_optimal_pairs)
    all_rows: list[dict[str, object]] = []
    for index, pair in enumerate(result.all_optimal_pairs, start=1):
        all_rows.append(
            {
                "threshold_pair_id": f"BOPT_{index:04d}",
                "T_s": pair[0],
                "T_r": pair[1],
                "objective_retained_frequencies": result.objective,
                "pareto_nondominated": pair in frontier,
                "representative": pair == result.representative,
                "tie_break_order": "minimum_T_s_then_minimum_T_r",
            }
        )
    pareto_rows = [row for row in all_rows if bool(row["pareto_nondominated"])]
    return all_rows, pareto_rows


def sensitivity_rows(
    records: Sequence[GeometryData],
    result: RuleBSearchResult,
    optimal_rows: Sequence[Mapping[str, object]],
    *,
    counts: OfflineCounts | None = None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object], bool]:
    definitions = (
        (PARTITION_BASELINE, lambda record: record.assigned_partition == PARTITION_BASELINE),
        (PARTITION_DIRECTED, lambda record: record.assigned_partition == PARTITION_DIRECTED),
        ("S3_12", lambda record: record.case_id == "S3_12"),
        ("S3_14", lambda record: record.case_id == "S3_14"),
        ("S3_combined", lambda record: record.case_id in LOCKED_HOLDOUT_IDS),
    )
    validation_records = [
        record
        for record in records
        if record.assigned_partition in {PARTITION_BASELINE, PARTITION_DIRECTED, PARTITION_HOLDOUT}
    ]
    output: list[dict[str, object]] = []
    equivalence: list[dict[str, object]] = []
    pair_predictions: dict[str, list[dict[str, object]]] = {}
    for pair_row in optimal_rows:
        pair = (finite_float(pair_row, "T_s"), finite_float(pair_row, "T_r"))
        pair_predictions[str(pair_row["threshold_pair_id"])] = predictions_for_rule_b(
            validation_records,
            pair,
            pair_id=str(pair_row["threshold_pair_id"]),
            counts=counts,
        )
    representative_row = next(row for row in optimal_rows if bool(row["representative"]))
    representative_id = str(representative_row["threshold_pair_id"])
    representative_predictions = {
        str(row["case_id"]): row for row in pair_predictions[representative_id]
    }

    def vector_hash(rows: Sequence[Mapping[str, object]], *, safety_only: bool) -> str:
        payload = [
            (
                str(row["case_id"]),
                bool(row["false_safe"]),
            )
            if safety_only
            else (
                str(row["case_id"]),
                int(row["N_true"]),
                int(row["N_hat"]),
            )
            for row in sorted(rows, key=lambda item: str(item["case_id"]))
        ]
        return hashlib.sha256(
            json.dumps(payload, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
        ).hexdigest()

    pair_vector_equal: list[bool] = []
    pair_safety_equal: list[bool] = []
    for pair_row in optimal_rows:
        pair_id = str(pair_row["threshold_pair_id"])
        pair = (finite_float(pair_row, "T_s"), finite_float(pair_row, "T_r"))
        predictions = pair_predictions[pair_id]
        prediction_equal_for_pair = all(
            int(row["N_hat"]) == int(representative_predictions[str(row["case_id"])]["N_hat"])
            for row in predictions
        )
        safety_equal_for_pair = all(
            bool(row["false_safe"])
            == bool(representative_predictions[str(row["case_id"])]["false_safe"])
            for row in predictions
        )
        pair_vector_equal.append(prediction_equal_for_pair)
        pair_safety_equal.append(safety_equal_for_pair)
        for row in predictions:
            representative = representative_predictions[str(row["case_id"])]
            equivalence.append(
                {
                    "threshold_pair_id": pair_id,
                    "representative": pair_id == representative_id,
                    "partition": row["assigned_partition"],
                    "case_id": row["case_id"],
                    "N_true": row["N_true"],
                    "N_hat": row["N_hat"],
                    "representative_N_hat": representative["N_hat"],
                    "prediction_equal": int(row["N_hat"]) == int(representative["N_hat"]),
                    "safety_equal": bool(row["false_safe"]) == bool(representative["false_safe"]),
                }
            )
        for partition, predicate in definitions:
            selected = [
                row
                for row in predictions
                if predicate(next(record for record in validation_records if record.case_id == row["case_id"]))
            ]
            representative_selected = [
                row
                for row in pair_predictions[representative_id]
                if predicate(next(record for record in validation_records if record.case_id == row["case_id"]))
            ]
            metrics = metrics_for_predictions(selected)
            prediction_hash = vector_hash(selected, safety_only=False)
            safety_hash = vector_hash(selected, safety_only=True)
            output.append(
                {
                    "threshold_pair_id": pair_id,
                    "representative": pair_id == representative_id,
                    "T_s": pair[0],
                    "T_r": pair[1],
                    "partition": partition,
                    **metrics,
                    "prediction_vector_hash": prediction_hash,
                    "safety_signature_hash": safety_hash,
                    "prediction_vector_equals_representative": prediction_hash
                    == vector_hash(representative_selected, safety_only=False),
                    "safety_signature_equals_representative": safety_hash
                    == vector_hash(representative_selected, safety_only=True),
                    "safety_result": "false_safe_observed" if int(metrics["false_safe_geometry_count"]) else "no_false_safe_observed",
                }
            )
    invariants = {
        "equal_optimum_pair_count": len(optimal_rows),
        "equal_optimum_pairs_checked": len(pair_predictions),
        "validation_summary_partition_count": len(definitions),
        "expected_sensitivity_row_count": len(optimal_rows) * len(definitions),
        "actual_sensitivity_row_count": len(output),
        "all_equal_optimum_pairs_checked": len(pair_predictions) == len(optimal_rows),
        "all_equal_optimum_prediction_vectors_identical": all(pair_vector_equal),
        "all_equal_optimum_safety_signatures_identical": all(pair_safety_equal),
    }
    if int(invariants["actual_sensitivity_row_count"]) != int(invariants["expected_sensitivity_row_count"]):
        raise AssertionError("equal-optimum sensitivity row count is incomplete")
    unstable = not (
        bool(invariants["all_equal_optimum_pairs_checked"])
        and bool(invariants["all_equal_optimum_prediction_vectors_identical"])
        and bool(invariants["all_equal_optimum_safety_signatures_identical"])
    )
    return output, equivalence, invariants, unstable


def rule_s_rule_b_equivalence_rows(
    records: Sequence[GeometryData],
    rule_s_predictions: Sequence[Mapping[str, object]],
    rule_b_predictions: Sequence[Mapping[str, object]],
    *,
    threshold_s: float,
    threshold_r: float,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    s_by_case = {str(row["case_id"]): row for row in rule_s_predictions}
    b_by_case = {str(row["case_id"]): row for row in rule_b_predictions}
    rows: list[dict[str, object]] = []
    for record in records:
        s_prediction = s_by_case[record.case_id]
        b_prediction = b_by_case[record.case_id]
        first_rejected = int(s_prediction["N_hat"]) + 1 if int(s_prediction["N_hat"]) < K_MAX else None
        rejected_row = (
            next(row for row in prefix_extrema(record) if int(row["prefix_n"]) == first_rejected)
            if first_rejected is not None
            else None
        )
        shear_value = (
            float(rejected_row["S_prefix_max_Pi_shear"]) if rejected_row is not None else float("nan")
        )
        rotary_value = (
            float(rejected_row["R_prefix_max_Pi_rotary"]) if rejected_row is not None else float("nan")
        )
        shear_threshold_exceeded = rejected_row is not None and shear_value > float(threshold_s)
        rotary_threshold_exceeded = rejected_row is not None and rotary_value > float(threshold_r)
        shear_binding = shear_threshold_exceeded
        rotary_binding = rotary_threshold_exceeded and not shear_threshold_exceeded
        trigger = (
            "Pi_shear"
            if shear_binding
            else "Pi_rotary"
            if rotary_binding
            else "none"
        )
        rows.append(
            {
                "partition": record.assigned_partition,
                "case_id": record.case_id,
                "N_true": record.N_true,
                "N_S": int(s_prediction["N_hat"]),
                "N_B": int(b_prediction["N_hat"]),
                "prediction_equal": int(s_prediction["N_hat"]) == int(b_prediction["N_hat"]),
                "first_rejected_mode": first_rejected if first_rejected is not None else "",
                "shear_value": shear_value,
                "rotary_value": rotary_value,
                "T_s": threshold_s,
                "representative_T_r": threshold_r,
                "shear_threshold_exceeded": shear_threshold_exceeded,
                "rotary_threshold_exceeded": rotary_threshold_exceeded,
                "shear_binding": shear_binding,
                "rotary_binding": rotary_binding,
                "trigger": trigger,
            }
        )
    summary = {
        "geometry_count": len(rows),
        "Rule_S_equals_Rule_B_on_all_checked_geometries": all(bool(row["prediction_equal"]) for row in rows),
        "Pi_rotary_nonbinding_on_all_checked_predictions": not any(bool(row["rotary_binding"]) for row in rows),
        "Pi_shear_binding_rejection_count": sum(bool(row["shear_binding"]) for row in rows),
        "Pi_rotary_binding_rejection_count": sum(bool(row["rotary_binding"]) for row in rows),
        "Pi_rotary_coexceedance_count": sum(
            bool(row["rotary_threshold_exceeded"]) and bool(row["shear_threshold_exceeded"])
            for row in rows
        ),
        "T_s": threshold_s,
        "representative_T_r": threshold_r,
        "scope": "49_quality_approved_checked_geometries",
    }
    return rows, summary


def rule_s_threshold_margin_rows(
    prefix_rows: Sequence[Mapping[str, object]],
    *,
    threshold_s: float,
) -> list[dict[str, object]]:
    definitions = (
        (PARTITION_DEVELOPMENT, lambda row: row["assigned_partition"] == PARTITION_DEVELOPMENT),
        (PARTITION_DIRECTED, lambda row: row["assigned_partition"] == PARTITION_DIRECTED),
        ("S3_combined", lambda row: row["case_id"] in LOCKED_HOLDOUT_IDS),
    )
    rows: list[dict[str, object]] = []
    for partition, predicate in definitions:
        selected = [row for row in prefix_rows if predicate(row)]
        accepted = [row for row in selected if float(row["S_prefix_max_Pi_shear"]) <= threshold_s]
        unsafe_above = [
            row
            for row in selected
            if row["prefix_label"] == "unsafe_prefix"
            and float(row["S_prefix_max_Pi_shear"]) > threshold_s
        ]
        accepted_nearest = max(accepted, key=lambda row: float(row["S_prefix_max_Pi_shear"]), default=None)
        unsafe_nearest = min(unsafe_above, key=lambda row: float(row["S_prefix_max_Pi_shear"]), default=None)
        accepted_value = (
            float(accepted_nearest["S_prefix_max_Pi_shear"])
            if accepted_nearest is not None
            else float("nan")
        )
        unsafe_value = (
            float(unsafe_nearest["S_prefix_max_Pi_shear"])
            if unsafe_nearest is not None
            else float("nan")
        )
        accepted_margin = threshold_s - accepted_value if math.isfinite(accepted_value) else float("nan")
        unsafe_margin = unsafe_value - threshold_s if math.isfinite(unsafe_value) else float("nan")
        denominator = max(abs(threshold_s), 1.0e-300)
        rows.append(
            {
                "partition": partition,
                "T_s": threshold_s,
                "nearest_accepted_case_id": accepted_nearest["case_id"] if accepted_nearest else "",
                "nearest_accepted_prefix_n": accepted_nearest["prefix_n"] if accepted_nearest else "",
                "nearest_accepted_S": accepted_value,
                "accepted_absolute_margin": accepted_margin,
                "accepted_relative_margin": accepted_margin / denominator,
                "nearest_unsafe_above_case_id": unsafe_nearest["case_id"] if unsafe_nearest else "",
                "nearest_unsafe_above_prefix_n": unsafe_nearest["prefix_n"] if unsafe_nearest else "",
                "nearest_unsafe_above_S": unsafe_value,
                "unsafe_absolute_margin": unsafe_margin,
                "unsafe_relative_margin": unsafe_margin / denominator,
            }
        )
    return rows


def rule_s_leave_one_geometry_out_rows(
    development: Sequence[GeometryData],
) -> tuple[list[dict[str, object]], int]:
    rows: list[dict[str, object]] = []
    comparisons = 0
    for held_out in sorted(development, key=lambda record: record.case_id):
        training = [record for record in development if record.case_id != held_out.case_id]
        search = exact_rule_s_search(training)
        n_hat, used, trigger = predict_rule_s(held_out, search.selected_threshold)
        comparisons += search.scalar_comparisons + used
        rows.append(
            {
                "held_out_case_id": held_out.case_id,
                "training_geometry_count": len(training),
                "candidate_count": len(search.candidates),
                "selected_T_s": search.selected_threshold,
                "training_objective": search.objective,
                "N_true": held_out.N_true,
                "N_hat": n_hat,
                "false_safe": n_hat > held_out.N_true,
                "false_safe_frequency_count": max(n_hat - held_out.N_true, 0),
                "conservative_loss": max(held_out.N_true - n_hat, 0),
                "first_rejected_mode": n_hat + 1 if n_hat < K_MAX else "",
                "trigger_component": trigger,
                "frozen_full_development_threshold_changed": False,
            }
        )
    return rows, comparisons


def partition_audit_rows(records: Sequence[GeometryData]) -> list[dict[str, object]]:
    return [
        {
            "source": record.source,
            "case_id": record.case_id,
            "epsilon_0": record.epsilon_0,
            "beta_deg": record.beta_deg,
            "mu": record.mu,
            "eta": record.eta,
            "canonical_geometry_key": canonical_geometry_key_text(record.key),
            "assigned_partition": record.assigned_partition,
            "duplicate_status": record.duplicate_status,
            "exclusion_reason": record.exclusion_reason,
            "quality_status": record.quality_status,
            "root_source": record.root_source,
            "predictor_source": record.predictor_source,
            "algorithm_version": ALGORITHM_VERSION,
        }
        for record in records
    ]


def offline_operation_rows(counts: OfflineCounts) -> list[dict[str, object]]:
    reconstruction = counts.reconstruction
    return [
        {
            "scope": "offline_research_calibration_total",
            "existing_CSV_rows_read": counts.existing_csv_rows_read,
            "existing_root_records_consumed": counts.existing_root_records_consumed,
            "EB_mode_reconstructions": reconstruction.EB_mode_reconstructions,
            "EB_characteristic_matrix_evaluations": reconstruction.EB_characteristic_matrix_evaluations,
            "SVD_6x6_calls": reconstruction.SVD_6x6_calls,
            "shape_grid_point_evaluations": reconstruction.EB_shape_grid_points,
            "quadrature_calls": reconstruction.quadrature_calls,
            "quadrature_point_evaluations": reconstruction.quadrature_point_evaluations,
            "rule_A_threshold_candidates": counts.rule_A_threshold_candidates,
            "rule_B_threshold_pairs_evaluated": counts.rule_B_threshold_pairs_evaluated,
            "rule_S_threshold_candidates": counts.rule_S_threshold_candidates,
            "scalar_comparisons": counts.scalar_comparisons,
            "dominance_comparisons": counts.dominance_comparisons,
            "predictor_reconstruction_cache_hits": counts.predictor_reconstruction_cache_hits,
            "specialized_shear_predictor_calls": counts.specialized_shear_predictor_calls,
            "specialized_shear_quadrature_calls": counts.specialized_shear_quadrature_calls,
            "specialized_shear_quadrature_point_evaluations": counts.specialized_shear_quadrature_point_evaluations,
            "input_files": ";".join(counts.input_files),
            "notes": "offline CSV/reconstruction/search costs only; new root calculations=0",
        }
    ]


def online_operation_rows(
    predictions: Sequence[Mapping[str, object]],
    *,
    n_shape_points: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    grouped: dict[tuple[str, str], list[Mapping[str, object]]] = defaultdict(list)
    for prediction in predictions:
        grouped[(str(prediction["assigned_partition"]), str(prediction["rule"]))].append(prediction)
        if prediction["assigned_partition"] == PARTITION_EXCLUDED:
            continue
        n_hat = int(prediction["N_hat"])
        mode_count = K_MAX if n_hat == K_MAX else n_hat + 1
        rule = str(prediction["rule"])
        quadrature_per_mode = 6 if rule == "Rule_S" else 10
        rows.append(
            {
                "scope": "estimated_online_selector_cost_geometry",
                "partition": prediction["assigned_partition"],
                "rule": rule,
                "case_id": prediction["case_id"],
                "estimated_online_mode_reconstructions": mode_count,
                "estimated_online_SVD_6x6_calls": mode_count,
                "estimated_online_quadrature_calls": quadrature_per_mode * mode_count,
                "estimated_online_scalar_comparisons": mode_count if rule == "Rule_A" else 2 * mode_count,
                "shape_grid_point_evaluations": 2 * int(n_shape_points) * mode_count,
                "quadrature_point_evaluations": quadrature_per_mode * int(n_shape_points) * mode_count,
                "Timoshenko_suffix_frequencies": K_MAX - n_hat,
                "notes": (
                    "per-geometry sequential Rule-S estimate; EB mode energy plus shear-only quadratures"
                    if rule == "Rule_S"
                    else "per-geometry sequential estimate; EB mode energy plus both Pi components"
                ),
            }
        )
    for (partition, rule), items in sorted(grouped.items()):
        if partition == PARTITION_EXCLUDED:
            continue
        mode_count = sum(K_MAX if int(row["N_hat"]) == K_MAX else int(row["N_hat"]) + 1 for row in items)
        scalar = mode_count if rule == "Rule_A" else 2 * mode_count
        if rule == "Rule_S":
            scalar = mode_count
        quadrature_per_mode = 6 if rule == "Rule_S" else 10
        rows.append(
            {
                "scope": "estimated_online_selector_cost_partition_total",
                "partition": partition,
                "rule": rule,
                "case_id": "ALL",
                "estimated_online_mode_reconstructions": mode_count,
                "estimated_online_SVD_6x6_calls": mode_count,
                "estimated_online_quadrature_calls": quadrature_per_mode * mode_count,
                "estimated_online_scalar_comparisons": scalar,
                "shape_grid_point_evaluations": 2 * int(n_shape_points) * mode_count,
                "quadrature_point_evaluations": quadrature_per_mode * int(n_shape_points) * mode_count,
                "Timoshenko_suffix_frequencies": sum(K_MAX - int(row["N_hat"]) for row in items),
                "notes": (
                    "EB roots already exist; specialized shear-only predictor for each sequentially considered mode"
                    if rule == "Rule_S"
                    else "EB roots already exist; both Pi components are formed for each sequentially considered mode"
                ),
            }
        )
    return rows


def create_plots(output_dir: Path) -> tuple[Path, ...]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prefix_rows = certification.read_csv(output_dir / "prefix_predictor_metrics.csv")
    rule_a_search = certification.read_csv(output_dir / "rule_A_exact_search.csv")
    optimal_rows = certification.read_csv(output_dir / "rule_B_optimal_threshold_pairs.csv")
    pareto_rows = certification.read_csv(output_dir / "rule_B_pareto_frontier.csv")
    predictions = [
        *certification.read_csv(output_dir / "rule_A_predictions.csv"),
        *certification.read_csv(output_dir / "rule_B_predictions.csv"),
    ]
    development = [row for row in prefix_rows if row["assigned_partition"] == PARTITION_DEVELOPMENT]
    selected_a = next(row for row in rule_a_search if bool_value(row["selected"]))

    paths: list[Path] = []
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for label, marker, color in (("safe_prefix", "o", "tab:blue"), ("unsafe_prefix", "x", "tab:red")):
        rows = [row for row in development if row["prefix_label"] == label]
        ax.scatter(
            [int(row["prefix_n"]) for row in rows],
            [finite_float(row, "P_prefix_max_Pi_EB") for row in rows],
            marker=marker,
            color=color,
            alpha=0.65,
            label=label,
        )
    threshold_a = finite_float(selected_a, "T_A")
    if math.isfinite(threshold_a):
        ax.axhline(threshold_a, color="black", linestyle="--", label="selected T_A")
    ax.set(xlabel="prefix n", ylabel="P_g,n = max Pi_EB", title="Rule A development prefix points")
    ax.legend()
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[0]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(7.0, 5.8))
    for label, marker, color in (("safe_prefix", "o", "tab:blue"), ("unsafe_prefix", "x", "tab:red")):
        rows = [row for row in development if row["prefix_label"] == label]
        ax.scatter(
            [finite_float(row, "S_prefix_max_Pi_shear") for row in rows],
            [finite_float(row, "R_prefix_max_Pi_rotary") for row in rows],
            marker=marker,
            color=color,
            alpha=0.6,
            label=label,
        )
    representative = next(row for row in optimal_rows if bool_value(row["representative"]))
    t_s, t_r = finite_float(representative, "T_s"), finite_float(representative, "T_r")
    ax.axvline(t_s, color="black", linestyle="--")
    ax.axhline(t_r, color="black", linestyle="--", label="selected rectangle")
    if len(pareto_rows) > 1:
        ordered = sorted(pareto_rows, key=lambda row: finite_float(row, "T_s"))
        ax.plot(
            [finite_float(row, "T_s") for row in ordered],
            [finite_float(row, "T_r") for row in ordered],
            color="tab:green",
            marker="s",
            linewidth=1.0,
            label="optimal threshold Pareto frontier",
        )
    ax.set(xlabel="S_g,n = max Pi_shear", ylabel="R_g,n = max Pi_rotary", title="Rule B development prefix plane")
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[1]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    markers = {PARTITION_DEVELOPMENT: "o", PARTITION_BASELINE: "s", PARTITION_DIRECTED: "^", PARTITION_HOLDOUT: "X"}
    colors = {"Rule_A": "tab:orange", "Rule_B": "tab:blue"}
    for (partition, rule), rows in _group_rows(predictions, ("assigned_partition", "rule")).items():
        if partition == PARTITION_EXCLUDED:
            continue
        ax.scatter(
            [int(float(row["N_true"])) for row in rows],
            [int(float(row["N_hat"])) for row in rows],
            marker=markers.get(partition, "."),
            color=colors.get(rule, "black"),
            alpha=0.65,
            label=f"{rule}: {partition}",
        )
    ax.plot([0, K_MAX], [0, K_MAX], color="black", linewidth=0.8)
    ax.set(xlabel="N_true", ylabel="N_hat", xlim=(-0.4, 10.4), ylim=(-0.4, 10.4), title="Frozen Rule A/B predictions")
    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=6)
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[2]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)
    return tuple(paths)


def _group_rows(
    rows: Sequence[Mapping[str, object]], fields: Sequence[str]
) -> dict[tuple[str, ...], list[Mapping[str, object]]]:
    grouped: dict[tuple[str, ...], list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[tuple(str(row[field]) for field in fields)].append(row)
    return grouped


def stop_go_status(
    rule_b_predictions: Sequence[Mapping[str, object]],
    *,
    unstable: bool,
) -> str:
    decisive = [
        row
        for row in rule_b_predictions
        if row["assigned_partition"] == PARTITION_DIRECTED or row["case_id"] in LOCKED_HOLDOUT_IDS
    ]
    if any(bool(row["false_safe"]) for row in decisive):
        return "rule_B_rejected_false_safe"
    if unstable:
        return "rule_B_rejected_unstable_optimum"
    return "rule_B_safety_survives_cost_test_required"


def rule_a_status(rule_a_predictions: Sequence[Mapping[str, object]]) -> str:
    checked = [row for row in rule_a_predictions if row["assigned_partition"] != PARTITION_DEVELOPMENT]
    return (
        "rule_A_false_safe_observed"
        if any(bool(row["false_safe"]) for row in checked)
        else "rule_A_no_false_safe_on_checked_sets_benchmark_only"
    )


def write_cost_proposal(output_dir: Path, records: Sequence[GeometryData]) -> Path:
    selected_ids = ("B01", "G04", "S3_06", "S3_12", "S3_14")
    lookup = {record.case_id: record for record in records}
    rows = [
        {
            "proposal_status": "proposal_only",
            "execution_status": "not_executed",
            "case_id": case_id,
            "epsilon_0": lookup[case_id].epsilon_0,
            "beta_deg": lookup[case_id].beta_deg,
            "mu": lookup[case_id].mu,
            "eta": lookup[case_id].eta,
            "purpose": "future_direct_Timoshenko_K10_cost_break_even_benchmark",
        }
        for case_id in selected_ids
        if case_id in lookup
    ]
    return write_csv(
        output_dir / COST_PROPOSAL_NAME,
        rows,
        ("proposal_status", "execution_status", "case_id", "epsilon_0", "beta_deg", "mu", "eta", "purpose"),
    )


def write_report(
    output_dir: Path,
    *,
    records: Sequence[GeometryData],
    rule_a: RuleASearchResult,
    rule_b: RuleBSearchResult,
    rule_s: RuleSSearchResult,
    dominance: Mapping[str, object],
    summaries: Sequence[Mapping[str, object]],
    summaries_s: Sequence[Mapping[str, object]],
    sensitivity_invariants: Mapping[str, object],
    rule_s_equivalence: Mapping[str, object],
    margin_rows: Sequence[Mapping[str, object]],
    leave_one_rows: Sequence[Mapping[str, object]],
    phase_i_gate: Mapping[str, object],
    phase_i_gate_passed: bool,
    rule_a_decision: str,
    rule_b_decision: str,
    unstable: bool,
    counts: OfflineCounts,
) -> Path:
    lookup = {
        (str(row["summary_partition"]), str(row["rule"])): row
        for row in [*summaries, *summaries_s]
    }
    reconstruction = counts.reconstruction
    lines = [
        "# Exact Rule A/B/S safe-prefix experiment",
        "",
        "## Scope and data lock",
        "",
        f"Algorithm: `{ALGORITHM_VERSION}`. Development uses exactly the 21-case branch-informed pilot. Directed validation uses the existing Step-3A cases after removing exact development duplicates, while `S3_12` and `S3_14` are a locked adversarial holdout. Geometry deduplication uses an absolute `{GEOMETRY_TOLERANCE:.0e}` tolerance independently on epsilon, beta, mu, and eta.",
        "",
        "No EB or Timoshenko eigenvalue problem, general certified search, branch continuation, Step 3B, FEM, or root solver was run. Step-3A predictors were reconstructed only from selected saved EB roots. Existing A-gap/C/D and epsilon certificates were not recalibrated or continued.",
        "",
        "## Exact calibration",
        "",
        f"Rule A evaluated `{len(rule_a.candidates)}` exact candidates and selected `T_A={rule_a.selected_threshold:.16e}` with development retained-frequency objective `{rule_a.objective}` and zero development false-safe.",
        "",
        f"Rule B evaluated `{rule_b.evaluated_pair_count}` exact Cartesian pairs on `{len(rule_b.shear_candidates)}` shear by `{len(rule_b.rotary_candidates)}` rotary candidates. The maximum zero-false-safe development objective is `{rule_b.objective}`. There are `{len(rule_b.all_optimal_pairs)}` equally optimal pairs and `{len(rule_b.pareto_optimal_pairs)}` Pareto-nondominated equally optimal pairs. The frozen representative is `T_s={rule_b.representative[0]:.16e}`, `T_r={rule_b.representative[1]:.16e}` by the predeclared minimum-T_s then minimum-T_r tie-break.",
        "",
        f"The independent exact one-dimensional Rule-S search evaluated `{len(rule_s.candidates)}` observed shear-prefix extrema plus reject-all, selected `T_s={rule_s.selected_threshold:.16e}`, and retained `{rule_s.objective}` of `{sum(record.N_true for record in records if record.assigned_partition == PARTITION_DEVELOPMENT)}` development frequencies with zero development false-safe.",
        "",
        "## Equal-optimum Rule-B audit",
        "",
        f"All `{sensitivity_invariants['equal_optimum_pairs_checked']}` of `{sensitivity_invariants['equal_optimum_pair_count']}` equally optimal pairs were applied to five validation summaries, producing `{sensitivity_invariants['actual_sensitivity_row_count']}` rows. All prediction vectors identical: `{str(bool(sensitivity_invariants['all_equal_optimum_prediction_vectors_identical'])).lower()}`. All safety signatures identical: `{str(bool(sensitivity_invariants['all_equal_optimum_safety_signatures_identical'])).lower()}`.",
        "",
        "## Rule-S extraction",
        "",
        f"Rule S equals frozen Rule B on all `{rule_s_equivalence['geometry_count']}` checked geometries: `{str(bool(rule_s_equivalence['Rule_S_equals_Rule_B_on_all_checked_geometries'])).lower()}`. Pi_rotary is nonbinding on all checked predictions: `{str(bool(rule_s_equivalence['Pi_rotary_nonbinding_on_all_checked_predictions'])).lower()}`. Thus Rule B degenerates to a shear-only Rule S on all checked data; this is an empirical finite-sample statement only.",
        "",
        "## Finite-sample separability audit",
        "",
        f"Development contains `{dominance['safe_prefix_point_count']}` safe and `{dominance['unsafe_prefix_point_count']}` unsafe prefix points. `{dominance['unavoidably_rejected_safe_prefix_count']}` safe points (`{float(dominance['unavoidably_rejected_safe_prefix_fraction']):.6f}`) have at least one componentwise-lower unsafe witness. The maximum theoretically attainable zero-false-safe Rule-B retention on this development sample is `{rule_b.objective}` frequencies.",
        "",
        "This is a **finite-sample limitation of the monotone rectangular Rule-B class**, not an impossibility proof over the continuous parameter domain.",
        "",
        "## Frozen validation",
        "",
        "| partition | rule | geometries | false-safe geometries | false-safe frequencies | worst overprediction | mean N_true | mean N_hat | retention | mean loss |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    order = (PARTITION_DEVELOPMENT, PARTITION_BASELINE, PARTITION_DIRECTED, "S3_12", "S3_14", "S3_combined")
    for partition in order:
        for rule in ("Rule_A", "Rule_B", "Rule_S"):
            row = lookup[(partition, rule)]
            lines.append(
                f"| {partition} | {rule} | {row['geometry_count']} | {row['false_safe_geometry_count']} | {row['false_safe_frequency_count']} | {row['worst_overprediction']} | {float(row['mean_N_true']):.4g} | {float(row['mean_N_hat']):.4g} | {float(row['usable_frequency_retention']):.4g} | {float(row['mean_conservative_loss']):.4g} |"
            )
    lines.extend(
        [
            "",
            f"Equal-optimum Rule-B validation sensitivity changes the checked safety signature: `{str(unstable).lower()}`. Thresholds were not changed after directed validation or locked-holdout inspection.",
            "",
            "## Rule-S margins and leave-one-geometry stability",
            "",
            "| partition | nearest accepted | accepted margin | nearest unsafe above | unsafe margin |",
            "|---|---|---:|---|---:|",
            *[
                f"| {row['partition']} | {row['nearest_accepted_case_id']} n={row['nearest_accepted_prefix_n']} | {float(row['accepted_absolute_margin']):.6e} | {row['nearest_unsafe_above_case_id']} n={row['nearest_unsafe_above_prefix_n']} | {float(row['unsafe_absolute_margin']):.6e} |"
                for row in margin_rows
            ],
            "",
            f"Leave-one-development-geometry-out recalibrated 21 train-20 thresholds without changing the frozen full-development threshold. False-safe held-out geometries: `{sum(bool(row['false_safe']) for row in leave_one_rows)}`; threshold range: `{min(float(row['selected_T_s']) for row in leave_one_rows):.16e}` to `{max(float(row['selected_T_s']) for row in leave_one_rows):.16e}`.",
            "",
            "## Decisions",
            "",
            f"- Rule A diagnostic status: `{rule_a_decision}`.",
            f"- Rule B primary status: `{rule_b_decision}`.",
            f"- Phase-I gate passed: `{str(phase_i_gate_passed).lower()}`.",
            *[f"- `{key}`: `{fmt(value)}`." for key, value in phase_i_gate.items()],
            "",
            "Zero observed false-safe on these finite checked sets would not be a guarantee over the continuous parameter domain. Rule A remains a benchmark and cannot become a production certificate from this experiment alone.",
            "",
            "## Offline operation counts",
            "",
            f"The postprocessor read `{counts.existing_csv_rows_read}` existing CSV rows and consumed `{counts.existing_root_records_consumed}` saved EB/Timoshenko root records. It performed `{reconstruction.EB_mode_reconstructions}` EB mode reconstructions, `{reconstruction.EB_characteristic_matrix_evaluations}` EB characteristic-matrix evaluations, `{reconstruction.SVD_6x6_calls}` 6x6 SVD calls, `{reconstruction.EB_shape_grid_points}` shape-grid point evaluations, `{reconstruction.quadrature_calls}` quadrature calls, and `{reconstruction.quadrature_point_evaluations}` quadrature-point evaluations. Predictor reconstruction cache hits: `{counts.predictor_reconstruction_cache_hits}`.",
            "",
            f"The specialized shear helper matched the existing predictor for all `{counts.specialized_shear_predictor_calls}` reconstructed Step-3A modes. It preserves the EB matrix/SVD and shape-grid work but reduces the online predictor path from 10 quadratures per mode for both Pi components (or 12 in the historical full predictor helper) to 6 total quadratures per mode: 4 mode-energy plus 2 shear integrals. Rule S is therefore logically simpler and saves four quadratures per reconstructed online mode relative to the minimal two-component Pi path; it does not save the 6x6 SVD or shape reconstruction.",
            "",
            f"Threshold search plus frozen and equal-optimum sensitivity applications used `{counts.scalar_comparisons}` scalar threshold comparisons. The dominance audit used `{counts.dominance_comparisons}` prefix-pair comparisons.",
            "",
            "Online selector costs are reported separately in `operation_counts.csv`; Timoshenko suffix counts are structural frequency counts and are not converted into a full Timoshenko operation cost.",
            "",
            "## Scientific limits",
            "",
            "`epsilon_0` is not a certificate. The historical A-gap/C/D paths are not continued here. This experiment corrects exact Rule B and formally extracts exact Rule S. No subgroup threshold, spectral-gap guard, modal-character guard, regression fit, ML predictor, or new composite scalar was introduced. Rule S remains a finite-sample validated selector, not a proven continuous-domain certificate.",
        ]
    )
    path = output_dir / REPORT_NAME
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def ensure_output_policy(args: Args) -> None:
    expected = (*NUMERICAL_OUTPUT_NAMES, *PLOT_NAMES, REPORT_NAME, COST_PROPOSAL_NAME)
    existing = [args.output_dir / name for name in expected if (args.output_dir / name).exists()]
    if existing and not args.force:
        raise FileExistsError("Rule A/B outputs already exist; pass --force to replace them")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def plot_only(args: Args) -> dict[str, object]:
    required = [args.output_dir / name for name in NUMERICAL_OUTPUT_NAMES]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"--plot-only requires saved numerical CSV files; missing: {missing}")
    plots = create_plots(args.output_dir)
    return {"args": args, "plot_only": True, "root_calculations": 0, "mode_reconstructions": 0, "plots": plots}


def write_inconclusive_quality_report(args: Args, error: Exception) -> Path:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    path = args.output_dir / REPORT_NAME
    path.write_text(
        "\n".join(
            (
                "# Exact Rule A/B safe-prefix experiment",
                "",
                "## Data-quality decision",
                "",
                "- Rule A status: `rule_A_inconclusive_data_quality`.",
                "- Rule B status: `rule_B_inconclusive_data_quality`.",
                "- New root calculations: `0`.",
                f"- Blocking existing-data condition: `{error}`",
                "",
                "The missing or insufficient saved data were not replaced with legacy roots, and no root solver was run.",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    return path


def run_analysis(args: Args) -> dict[str, object]:
    ensure_output_policy(args)
    counts = OfflineCounts()
    development = load_development(args, counts)
    step3a = load_step3a(args, counts)
    assign_partitions(development, step3a)
    predictor_equivalence = reconstruct_step3a_predictors(
        step3a,
        n_shape_points=args.n_shape_points,
        counts=counts,
    )
    all_records = [*development, *step3a]
    included_records = [record for record in all_records if record.assigned_partition != PARTITION_EXCLUDED]
    development_prefix = [row for record in development for row in prefix_extrema(record)]
    all_prefix = [row for record in included_records for row in prefix_extrema(record)]

    rule_a = exact_rule_a_search(development)
    rule_b = exact_rule_b_pareto_search(development)
    rule_s = exact_rule_s_search(development)
    counts.rule_A_threshold_candidates = len(rule_a.candidates)
    counts.rule_B_threshold_pairs_evaluated = rule_b.evaluated_pair_count
    counts.rule_S_threshold_candidates = len(rule_s.candidates)
    counts.scalar_comparisons = (
        rule_a.scalar_comparisons + rule_b.scalar_comparisons + rule_s.scalar_comparisons
    )
    witnesses, dominance = dominance_audit(development_prefix)
    counts.dominance_comparisons = int(dominance["dominance_comparisons"])

    optimal_rows, pareto_rows = threshold_pair_rows(rule_b)
    representative_row = next(row for row in optimal_rows if bool(row["representative"]))
    representative_id = str(representative_row["threshold_pair_id"])
    predictions_a = predictions_for_rule_a(
        included_records,
        rule_a.selected_threshold,
        counts=counts,
    )
    predictions_b = predictions_for_rule_b(
        included_records,
        rule_b.representative,
        pair_id=representative_id,
        counts=counts,
    )
    predictions_s = predictions_for_rule_s(
        included_records,
        rule_s.selected_threshold,
        counts=counts,
    )
    summaries = summary_partition_rows(
        predictions_a,
        predictions_b,
        threshold_a=rule_a.selected_threshold,
        threshold_b=rule_b.representative,
    )
    summaries_s = rule_s_summary_rows(predictions_s, threshold_s=rule_s.selected_threshold)
    sensitivity, equal_optimum_equivalence, sensitivity_invariants, unstable = sensitivity_rows(
        included_records,
        rule_b,
        optimal_rows,
        counts=counts,
    )
    rule_s_rule_b_rows, rule_s_equivalence = rule_s_rule_b_equivalence_rows(
        included_records,
        predictions_s,
        predictions_b,
        threshold_s=rule_s.selected_threshold,
        threshold_r=rule_b.representative[1],
    )
    margin_rows = rule_s_threshold_margin_rows(all_prefix, threshold_s=rule_s.selected_threshold)
    leave_one_rows, leave_one_comparisons = rule_s_leave_one_geometry_out_rows(development)
    counts.scalar_comparisons += leave_one_comparisons
    decision_b = stop_go_status(predictions_b, unstable=unstable)
    decision_a = rule_a_status(predictions_a)
    corrected_representative = rule_b.representative == min(
        rule_b.all_optimal_pairs,
        key=lambda pair: (pair[0], pair[1]),
    )
    rule_s_exact_completed = (
        math.isfinite(rule_s.selected_threshold)
        and rule_s.objective == sum(int(row["N_hat"]) for row in predictions_s if row["assigned_partition"] == PARTITION_DEVELOPMENT)
        and not any(
            bool(row["false_safe"])
            for row in predictions_s
            if row["assigned_partition"] == PARTITION_DEVELOPMENT
        )
    )
    phase_i_gate = {
        "corrected_representative_pair": corrected_representative,
        "all_equal_optimum_pairs_checked": bool(sensitivity_invariants["all_equal_optimum_pairs_checked"]),
        "all_equal_optimum_prediction_vectors_identical": bool(
            sensitivity_invariants["all_equal_optimum_prediction_vectors_identical"]
        ),
        "all_equal_optimum_safety_signatures_identical": bool(
            sensitivity_invariants["all_equal_optimum_safety_signatures_identical"]
        ),
        "Rule_S_exact_search_completed": rule_s_exact_completed,
        "Rule_S_threshold_equals_Rule_B_T_s": math.isclose(
            rule_s.selected_threshold,
            rule_b.representative[0],
            rel_tol=0.0,
            abs_tol=0.0,
        ),
        "Rule_S_equals_Rule_B_on_all_checked_geometries": bool(
            rule_s_equivalence["Rule_S_equals_Rule_B_on_all_checked_geometries"]
        ),
        "Pi_rotary_nonbinding_on_all_checked_predictions": bool(
            rule_s_equivalence["Pi_rotary_nonbinding_on_all_checked_predictions"]
        ),
        "specialized_Pi_shear_matches_existing_for_all_280_modes": len(predictor_equivalence) == 280
        and all(bool(row["equivalent"]) for row in predictor_equivalence),
        "new_root_calculations_in_Phase_I": 0,
    }
    phase_i_gate_passed = all(
        value == 0 if key == "new_root_calculations_in_Phase_I" else bool(value)
        for key, value in phase_i_gate.items()
    )

    audit_rows = partition_audit_rows(all_records)
    included_audit = [row for row in audit_rows if row["assigned_partition"] != PARTITION_EXCLUDED]
    excluded_audit = [row for row in audit_rows if row["assigned_partition"] == PARTITION_EXCLUDED]
    summary_b = {
        "algorithm_version": ALGORITHM_VERSION,
        "search_is_exact": True,
        "shear_candidate_count": len(rule_b.shear_candidates),
        "rotary_candidate_count": len(rule_b.rotary_candidates),
        "candidate_pair_count": rule_b.evaluated_pair_count,
        "evaluated_pair_count": rule_b.evaluated_pair_count,
        "rejected_pairs_with_false_safe": rule_b.rejected_false_safe_count,
        "maximum_objective_retained_frequencies": rule_b.objective,
        "all_equally_optimal_pair_count": len(rule_b.all_optimal_pairs),
        "pareto_nondominated_optimal_pair_count": len(rule_b.pareto_optimal_pairs),
        "selected_T_s": rule_b.representative[0],
        "selected_T_r": rule_b.representative[1],
        "representative_pair_id": representative_id,
        "tie_break": "minimum_T_s_then_minimum_T_r",
        "threshold_provenance": "development_only_observed_prefix_extrema_plus_reject_all",
        "rule_B_status": decision_b,
        **sensitivity_invariants,
        **phase_i_gate,
        "phase_I_gate_passed": phase_i_gate_passed,
        **dominance,
        "maximum_theoretically_attainable_zero_false_safe_retention": rule_b.objective,
    }
    summary_s = {
        "algorithm_version": ALGORITHM_VERSION,
        "search_is_exact": True,
        "candidate_count": len(rule_s.candidates),
        "selected_T_s": rule_s.selected_threshold,
        "maximum_objective_retained_frequencies": rule_s.objective,
        "development_true_frequency_count": sum(record.N_true for record in development),
        **rule_s_equivalence,
        **phase_i_gate,
        "phase_I_gate_passed": phase_i_gate_passed,
        "scientific_status": (
            "rule_S_finite_sample_validated_selector"
            if phase_i_gate_passed
            else "rule_S_phase_I_gate_failed"
        ),
    }

    write_csv(args.output_dir / "data_partition_audit.csv", audit_rows, PARTITION_AUDIT_FIELDS)
    write_csv(args.output_dir / "included_geometry_audit.csv", included_audit, PARTITION_AUDIT_FIELDS)
    write_csv(args.output_dir / "excluded_geometry_audit.csv", excluded_audit, PARTITION_AUDIT_FIELDS)
    write_csv(args.output_dir / "prefix_predictor_metrics.csv", all_prefix, PREFIX_FIELDS)
    write_csv(
        args.output_dir / "rule_A_exact_search.csv",
        rule_a.rows,
        (
            "candidate_index", "T_A", "zero_development_false_safe", "false_safe_geometry_count",
            "false_safe_frequency_count", "retained_frequencies", "conservative_loss", "selected",
            "search_is_exact", "threshold_provenance",
        ),
    )
    write_csv(args.output_dir / "rule_A_predictions.csv", predictions_a, PREDICTION_FIELDS)
    write_csv(args.output_dir / "rule_B_exact_search_summary.csv", [summary_b], tuple(summary_b))
    pair_fields = (
        "threshold_pair_id", "T_s", "T_r", "objective_retained_frequencies",
        "pareto_nondominated", "representative", "tie_break_order",
    )
    write_csv(args.output_dir / "rule_B_optimal_threshold_pairs.csv", optimal_rows, pair_fields)
    write_csv(args.output_dir / "rule_B_pareto_frontier.csv", pareto_rows, pair_fields)
    write_csv(
        args.output_dir / "rule_B_dominance_witnesses.csv",
        witnesses,
        (
            "safe_case_id", "safe_prefix_n", "safe_S", "safe_R", "unsafe_case_id",
            "unsafe_prefix_n", "unsafe_S", "unsafe_R", "dominance_relation",
        ),
    )
    write_csv(args.output_dir / "rule_B_predictions.csv", predictions_b, PREDICTION_FIELDS)
    write_csv(
        args.output_dir / "rule_B_equal_optimum_prediction_equivalence.csv",
        equal_optimum_equivalence,
        (
            "threshold_pair_id", "representative", "partition", "case_id", "N_true", "N_hat",
            "representative_N_hat", "prediction_equal", "safety_equal",
        ),
    )
    write_csv(
        args.output_dir / "rule_S_exact_search.csv",
        rule_s.rows,
        (
            "candidate_index", "T_s", "zero_development_false_safe", "false_safe_geometry_count",
            "false_safe_frequency_count", "retained_frequencies", "conservative_loss", "selected",
            "search_is_exact", "threshold_provenance",
        ),
    )
    write_csv(args.output_dir / "rule_S_predictions.csv", predictions_s, PREDICTION_FIELDS)
    write_csv(args.output_dir / "rule_S_validation_summary.csv", summaries_s, SUMMARY_FIELDS)
    write_csv(args.output_dir / "rule_S_equivalence_audit.csv", [summary_s], tuple(summary_s))
    write_csv(
        args.output_dir / "rule_S_rule_B_equivalence_audit.csv",
        rule_s_rule_b_rows,
        (
            "partition", "case_id", "N_true", "N_S", "N_B", "prediction_equal",
            "first_rejected_mode", "shear_value", "rotary_value", "T_s", "representative_T_r",
            "shear_threshold_exceeded", "rotary_threshold_exceeded",
            "shear_binding", "rotary_binding", "trigger",
        ),
    )
    write_csv(
        args.output_dir / "rule_S_threshold_margin_audit.csv",
        margin_rows,
        (
            "partition", "T_s", "nearest_accepted_case_id", "nearest_accepted_prefix_n",
            "nearest_accepted_S", "accepted_absolute_margin", "accepted_relative_margin",
            "nearest_unsafe_above_case_id", "nearest_unsafe_above_prefix_n", "nearest_unsafe_above_S",
            "unsafe_absolute_margin", "unsafe_relative_margin",
        ),
    )
    write_csv(
        args.output_dir / "rule_S_leave_one_geometry_out.csv",
        leave_one_rows,
        (
            "held_out_case_id", "training_geometry_count", "candidate_count", "selected_T_s",
            "training_objective", "N_true", "N_hat", "false_safe", "false_safe_frequency_count",
            "conservative_loss", "first_rejected_mode", "trigger_component",
            "frozen_full_development_threshold_changed",
        ),
    )
    write_csv(
        args.output_dir / "rule_S_predictor_implementation_equivalence.csv",
        predictor_equivalence,
        (
            "source", "case_id", "partition", "sorted_index", "Lambda_EB",
            "Pi_shear_existing", "Pi_shear_specialized", "absolute_difference", "relative_difference",
            "relative_tolerance", "absolute_tolerance", "equivalent",
            "existing_predictor_quadrature_calls", "specialized_rule_S_total_quadrature_calls",
            "specialized_rule_S_shear_quadrature_calls",
        ),
    )
    write_csv(args.output_dir / "partition_validation_summary.csv", summaries, SUMMARY_FIELDS)
    sensitivity_fields = (
        "threshold_pair_id", "representative", "T_s", "T_r", "partition", "geometry_count",
        "false_safe_geometry_count", "false_safe_frequency_count", "worst_overprediction",
        "mean_N_true", "mean_N_hat", "exact_match_rate", "mean_conservative_loss",
        "maximum_conservative_loss", "usable_frequency_retention", "N_hat_zero_count",
        "N_hat_ten_count", "prediction_vector_hash", "safety_signature_hash",
        "prediction_vector_equals_representative", "safety_signature_equals_representative",
        "safety_result",
    )
    write_csv(args.output_dir / "threshold_pair_validation_sensitivity.csv", sensitivity, sensitivity_fields)
    operation_rows = [
        *offline_operation_rows(counts),
        *online_operation_rows(
            [*predictions_a, *predictions_b, *predictions_s],
            n_shape_points=args.n_shape_points,
        ),
    ]
    write_csv(args.output_dir / "operation_counts.csv", operation_rows, OPERATION_FIELDS)
    plots = create_plots(args.output_dir)
    report = write_report(
        args.output_dir,
        records=included_records,
        rule_a=rule_a,
        rule_b=rule_b,
        rule_s=rule_s,
        dominance=dominance,
        summaries=summaries,
        summaries_s=summaries_s,
        sensitivity_invariants=sensitivity_invariants,
        rule_s_equivalence=rule_s_equivalence,
        margin_rows=margin_rows,
        leave_one_rows=leave_one_rows,
        phase_i_gate=phase_i_gate,
        phase_i_gate_passed=phase_i_gate_passed,
        rule_a_decision=decision_a,
        rule_b_decision=decision_b,
        unstable=unstable,
        counts=counts,
    )
    proposal = (
        write_cost_proposal(args.output_dir, included_records)
        if decision_b == "rule_B_safety_survives_cost_test_required" and phase_i_gate_passed
        else None
    )
    if proposal is None and (args.output_dir / COST_PROPOSAL_NAME).exists():
        (args.output_dir / COST_PROPOSAL_NAME).unlink()
    return {
        "args": args,
        "records": included_records,
        "rule_a": rule_a,
        "rule_b": rule_b,
        "rule_s": rule_s,
        "dominance": dominance,
        "summaries": summaries,
        "rule_A_status": decision_a,
        "rule_B_status": decision_b,
        "unstable_optimum": unstable,
        "phase_I_gate": phase_i_gate,
        "phase_I_gate_passed": phase_i_gate_passed,
        "counts": counts,
        "plots": plots,
        "report": report,
        "cost_proposal": proposal,
        "root_calculations": 0,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    if args.plot_only:
        result = plot_only(args)
    else:
        try:
            result = run_analysis(args)
        except (DataQualityError, FileNotFoundError) as exc:
            report = write_inconclusive_quality_report(args, exc)
            result = {
                "args": args,
                "rule_A_status": "rule_A_inconclusive_data_quality",
                "rule_B_status": "rule_B_inconclusive_data_quality",
                "data_quality_error": str(exc),
                "report": report,
                "root_calculations": 0,
            }
    print(f"algorithm: {ALGORITHM_VERSION}")
    print(f"output: {args.output_dir}")
    print("new root calculations: 0")
    if not args.plot_only:
        print(f"Rule A status: {result['rule_A_status']}")
        print(f"Rule B status: {result['rule_B_status']}")
    return result


if __name__ == "__main__":
    try:
        main()
    except (FileExistsError, FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from None
