from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from dataclasses import dataclass
import json
import math
from math import isfinite
from pathlib import Path
import statistics
import sys
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
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


DEFAULT_PILOT_DIR = Path("results") / "eb_epsilon_apriori_pilot"
DEFAULT_OUTPUT_DIR = DEFAULT_PILOT_DIR / "analysis"
SMOKE_PILOT_DIR = Path("results") / "_smoke" / "eb_epsilon_apriori_pilot"
SMOKE_OUTPUT_DIR = SMOKE_PILOT_DIR / "analysis"
DEFAULT_K_MAX = 10
DEFAULT_FREQUENCY_ERROR_THRESHOLD = 0.10
DEFAULT_CANDIDATE_GRID_SIZES = (8, 16, 32)

E_RULES = ("Rule_E0_ref", "Rule_Emax_ref", "Rule_E0_cal", "Rule_Emax_cal")
EB_RULES = certification.RULES
ALL_RULES = (*E_RULES, *EB_RULES)

OUTPUT_NAMES = (
    "epsilon_rule_thresholds.csv",
    "epsilon_rule_fold_predictions.csv",
    "epsilon_rule_validation_summary.csv",
    "epsilon_reference_curve.csv",
    "matched_epsilon_max_group_audit.csv",
    "same_epsilon0_group_audit.csv",
    "epsilon_vs_existing_rules_geometry_comparison.csv",
    "epsilon_rule_operation_counts.csv",
    "epsilon_calibration_search_convergence.csv",
    "epsilon_apriori_pilot_report.md",
    "N_true_vs_epsilon0.png",
    "N_true_vs_epsilon_max.png",
    "baseline_delta_f_vs_epsilon.png",
    "matched_epsilon_max_groups.png",
    "rule_false_safe_retention_comparison.png",
    "geometry_rule_vs_EB_rule_prefixes.png",
)

THRESHOLD_FIELDS = [
    "split_kind",
    "fold_id",
    "fold_status",
    "rule",
    "predictor",
    "mode_n",
    "threshold",
    "train_X_min",
    "train_X_max",
    "exact_calibration",
    "monotonicity_status",
    "train_false_safe_geometry_count",
    "train_retained_frequency_objective",
    "notes",
]

PREDICTION_FIELDS = [
    "split_kind",
    "fold_id",
    "fold_status",
    "evaluation_scope",
    "rule",
    "case_id",
    "geometry_id",
    "epsilon_0",
    "epsilon_max",
    "beta_deg",
    "mu",
    "eta",
    "N_true",
    "N_hat",
    "N_hat_guarded",
    "N_hat_unguarded",
    "false_safe",
    "false_safe_frequency_count",
    "conservative_loss",
    "out_of_X_calibration_domain",
    "first_predicted_failure",
    "trigger_reason",
    "cluster_guard_triggered",
    "mixed_mode_guard_triggered",
    "thresholds_json",
    "notes",
]

SUMMARY_FIELDS = [
    "split_kind",
    "rule",
    "fold_count",
    "fold_statuses",
    "validation_geometry_count",
    "false_safe_geometry_count",
    "false_safe_frequency_count",
    "false_safe_geometry_rate",
    "false_safe_frequency_rate",
    "mean_N_true",
    "mean_N_hat",
    "exact_N_match_rate",
    "mean_conservative_loss",
    "median_conservative_loss",
    "max_conservative_loss",
    "usable_EB_frequency_count",
    "retained_EB_frequency_count",
    "usable_frequency_retention",
    "fully_EB_geometry_count",
    "zero_EB_geometry_count",
    "out_of_X_calibration_domain_count",
    "cluster_guard_trigger_count",
    "mixed_mode_guard_trigger_count",
    "calibration_search_converged",
    "notes",
]

REFERENCE_FIELDS = [
    "case_id",
    "epsilon",
    "epsilon_max",
    "N_true_raw",
    "N_true_monotone_reference",
    "reference_N_from_thresholds",
    "first_failed_mode",
    "late_pass_after_failure",
    *[f"delta_f_{index}" for index in range(1, 11)],
    "thresholds_json",
    "notes",
]

MATCHED_FIELDS = [
    "group",
    "epsilon_max_target",
    "case_ids",
    "N_true_values",
    "range_N_true",
    "std_N_true",
    "maximum_per_mode_delta_f_spread",
    "mode_of_maximum_spread",
    "first_failed_mode_identical",
    "first_failed_modes",
    "reordering_warning",
    "cluster_warning",
    "notes",
]

SAME_EPSILON_FIELDS = [
    "case_id",
    "epsilon_0",
    "epsilon_max",
    "beta_deg",
    "mu",
    "eta",
    "N_true",
    "first_failed_mode",
    "group_range_N_true",
    "group_std_N_true",
    "epsilon_max_N_true_rank_correlation",
    "same_epsilon_material_variation",
    "notes",
]

COMPARISON_FIELDS = [
    "case_id",
    "epsilon_0",
    "epsilon_max",
    "beta_deg",
    "mu",
    "eta",
    "N_true",
    "N_hat_E0_ref",
    "N_hat_Emax_ref",
    "N_hat_E0_cal",
    "N_hat_Emax_cal",
    "N_hat_Rule_A",
    "N_hat_Rule_A_gap",
    "N_hat_Rule_B",
    "N_hat_Rule_C",
    "N_hat_Rule_D",
    "notes",
]

OPERATION_FIELDS = [
    "cost_scope",
    "cost_status",
    "split_kind",
    "fold_id",
    "rule",
    "geometry_count",
    "frequency_count",
    "geometry_rule_X_evaluations",
    "tau_factor_evaluations",
    "local_epsilon_evaluations",
    "max_operations",
    "threshold_comparisons",
    "EB_root_operations_if_available",
    "EB_characteristic_matrix_evaluations",
    "EB_determinant_evaluations",
    "Timoshenko_characteristic_matrix_evaluations",
    "Timoshenko_determinant_evaluations",
    "root_bracketing_intervals",
    "Brent_function_evaluations",
    "EB_mode_reconstructions",
    "SVD_6x6_calls",
    "EB_shape_grid_points",
    "quadrature_calls",
    "quadrature_point_evaluations",
    "scalar_predictor_comparisons",
    "spectral_gap_operations",
    "local_Timoshenko_root_attempts",
    "full_Timoshenko_fallbacks",
    "Timoshenko_frequencies_required",
    "Timoshenko_frequencies_avoided",
    "threshold_combinations_evaluated",
    "source_root_cache_hits",
    "source_root_cache_misses",
    "source_boundary_retries",
    "notes",
]

CONVERGENCE_FIELDS = [
    "split_kind",
    "fold_id",
    "fold_status",
    "rule",
    "max_candidates_per_axis",
    "retained_frequency_objective",
    "thresholds_json",
    "calibration_false_safe_geometry_count",
    "candidate_count_per_axis",
    "evaluated_combinations",
    "search_is_exact",
    "prediction_changes_from_previous_grid",
    "objective_changed_from_previous_grid",
    "calibration_search_converged",
    "notes",
]


@dataclass(frozen=True)
class Args:
    pilot_dir: Path
    output_dir: Path
    k_max: int
    frequency_error_threshold: float
    candidate_grid_sizes: tuple[int, ...]
    force: bool
    smoke: bool


@dataclass(frozen=True)
class PilotFold:
    split_kind: str
    fold_id: str
    train: tuple[certification.GeometryRecord, ...]
    validation: tuple[certification.GeometryRecord, ...]
    optimistic_in_sample: bool = False


def repo_path(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def fmt(value: object) -> object:
    return certification.fmt(value)


def read_csv(path: Path) -> list[dict[str, str]]:
    return certification.read_csv(path)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def finite_float(row: Mapping[str, object], key: str) -> float:
    return certification.finite_float(row, key)


def bool_value(value: object) -> bool:
    return certification.bool_value(value)


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Post-process the EB epsilon a-priori pilot CSV files without solving roots.",
    )
    parser.add_argument("--pilot-dir", type=Path, default=DEFAULT_PILOT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--k-max", type=int, default=DEFAULT_K_MAX)
    parser.add_argument("--frequency-error-threshold", type=float, default=DEFAULT_FREQUENCY_ERROR_THRESHOLD)
    parser.add_argument("--candidate-grid-sizes", type=int, nargs="+", default=DEFAULT_CANDIDATE_GRID_SIZES)
    parser.add_argument("--force", action="store_true", help="Allow replacement of existing analysis outputs.")
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    pilot_dir = repo_path(Path(ns.pilot_dir))
    output_dir = repo_path(Path(ns.output_dir))
    if bool(ns.smoke):
        if Path(ns.pilot_dir) == DEFAULT_PILOT_DIR:
            pilot_dir = repo_path(SMOKE_PILOT_DIR)
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            output_dir = repo_path(SMOKE_OUTPUT_DIR)
    grids = tuple(sorted(set(int(value) for value in ns.candidate_grid_sizes)))
    args = Args(
        pilot_dir=pilot_dir,
        output_dir=output_dir,
        k_max=int(ns.k_max),
        frequency_error_threshold=float(ns.frequency_error_threshold),
        candidate_grid_sizes=grids,
        force=bool(ns.force),
        smoke=bool(ns.smoke),
    )
    if args.k_max != 10:
        raise ValueError("this pilot contract requires --k-max 10")
    if not math.isclose(args.frequency_error_threshold, 0.10, rel_tol=0.0, abs_tol=1.0e-12):
        raise ValueError("this pilot contract requires --frequency-error-threshold 0.10")
    if not grids or any(value < 3 for value in grids):
        raise ValueError("candidate grid sizes must be integers >= 3")
    if not args.smoke and not {8, 16, 32}.issubset(grids):
        raise ValueError("normal pilot analysis requires candidate grids 8, 16, and 32")
    return args


def ensure_output_policy(args: Args) -> None:
    existing = [args.output_dir / name for name in OUTPUT_NAMES if (args.output_dir / name).exists()]
    if existing and not args.force:
        names = ", ".join(path.name for path in existing)
        raise FileExistsError(f"pilot analysis outputs already exist ({names}); pass --force to replace them")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def load_pilot_records(
    args: Args,
) -> tuple[list[certification.GeometryRecord], dict[str, dict[str, str]], dict[str, list[dict[str, str]]], list[dict[str, str]]]:
    manifest_rows = read_csv(args.pilot_dir / "epsilon_pilot_case_manifest_resolved.csv")
    mode_rows = read_csv(args.pilot_dir / "epsilon_pilot_mode_metrics.csv")
    geometry_rows = read_csv(args.pilot_dir / "epsilon_pilot_geometry_metrics.csv")
    root_operation_rows = read_csv(args.pilot_dir / "epsilon_pilot_root_operation_counts.csv")
    manifest_by_case = {row["case_id"]: row for row in manifest_rows}
    geometry_by_case = {row["case_id"]: row for row in geometry_rows}
    modes_by_case: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in mode_rows:
        modes_by_case[row["case_id"]].append(row)

    records: list[certification.GeometryRecord] = []
    for case_id in sorted(manifest_by_case):
        manifest = manifest_by_case[case_id]
        geometry = geometry_by_case.get(case_id)
        if geometry is None or str(geometry.get("quality_status", "")) != "included":
            continue
        modes = sorted(modes_by_case.get(case_id, []), key=lambda row: int(float(row["sorted_index"])))
        if len(modes) != args.k_max:
            continue
        key = certification.geometry_key(
            finite_float(manifest, "epsilon_0"),
            finite_float(manifest, "beta_deg"),
            finite_float(manifest, "mu"),
            finite_float(manifest, "eta"),
        )
        n_true = int(float(geometry["N_true"]))
        records.append(
            certification.GeometryRecord(
                source_dataset="epsilon_apriori_pilot",
                geometry_id=str(manifest["geometry_id"]),
                geometry_key=key,
                epsilon_0=key[0],
                beta_deg=key[1],
                mu=key[2],
                eta=key[3],
                modes=[dict(row) for row in modes],
                source_file=args.pilot_dir / "epsilon_pilot_mode_metrics.csv",
                complete=True,
                included=True,
                explicit_warning=False,
                N_true=n_true,
            )
        )
    return records, manifest_by_case, modes_by_case, root_operation_rows


def case_id_for(record: certification.GeometryRecord, manifest_by_case: Mapping[str, Mapping[str, object]]) -> str:
    for case_id, row in manifest_by_case.items():
        if str(row.get("geometry_id", "")) == record.geometry_id:
            return case_id
    raise KeyError(record.geometry_id)


def groups_for(record: certification.GeometryRecord, manifest_by_case: Mapping[str, Mapping[str, object]]) -> set[str]:
    case_id = case_id_for(record, manifest_by_case)
    return {part for part in str(manifest_by_case[case_id].get("case_groups", "")).split(";") if part}


def epsilon_max(record: certification.GeometryRecord) -> float:
    return finite_float(record.modes[0], "epsilon_max")


def build_folds(
    records: Sequence[certification.GeometryRecord],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> list[PilotFold]:
    ordered = tuple(sorted(records, key=lambda item: case_id_for(item, manifest_by_case)))
    baseline = tuple(record for record in ordered if "baseline_reference" in groups_for(record, manifest_by_case))
    nonbaseline = tuple(record for record in ordered if record not in baseline)
    folds = [PilotFold("baseline_reference_transfer", "baseline_to_nonbaseline", baseline, nonbaseline)]

    def leave_one(split_kind: str, attr: str, *, minimum_validation_size: int) -> None:
        values = sorted({float(getattr(record, attr)) for record in ordered})
        for value in values:
            validation = tuple(record for record in ordered if math.isclose(float(getattr(record, attr)), value, abs_tol=1.0e-12))
            if len(validation) < minimum_validation_size:
                continue
            train = tuple(record for record in ordered if record not in validation)
            token = f"{value:.12g}".replace("-", "m").replace(".", "p")
            folds.append(PilotFold(split_kind, f"holdout_{attr}_{token}", train, validation))

    leave_one("leave_one_beta_out", "beta_deg", minimum_validation_size=2)
    leave_one("leave_one_mu_out", "mu", minimum_validation_size=2)
    leave_one("leave_one_epsilon_out", "epsilon_0", minimum_validation_size=2)
    leave_one("leave_one_eta_out", "eta", minimum_validation_size=1)
    folds.append(PilotFold("combined_in_sample", "optimistic", ordered, ordered, optimistic_in_sample=True))
    return folds


def fold_status(fold: PilotFold) -> str:
    if not fold.train:
        return "no_training_geometries"
    if not fold.validation:
        return "no_validation_geometries"
    if len(fold.train) < 3:
        return "limited_training_size"
    targets = {int(record.N_true or 0) for record in fold.train}
    if len(targets) < 2:
        return "limited_train_target_variation"
    return "ok"


def predictor_name(rule: str) -> str:
    if rule in {"Rule_E0_ref", "Rule_E0_cal"}:
        return "epsilon_0"
    if rule in {"Rule_Emax_ref", "Rule_Emax_cal"}:
        return "epsilon_max"
    raise ValueError(rule)


def predictor_value(record: certification.GeometryRecord, rule: str) -> float:
    return record.epsilon_0 if predictor_name(rule) == "epsilon_0" else epsilon_max(record)


def calibrate_epsilon_thresholds(
    records: Sequence[certification.GeometryRecord],
    *,
    predictor: str,
    k_max: int,
) -> tuple[float, ...]:
    if not records:
        return tuple(float("-inf") for _ in range(k_max))
    pairs = [
        (record.epsilon_0 if predictor == "epsilon_0" else epsilon_max(record), int(record.N_true or 0))
        for record in records
    ]
    maximum_x = max(value for value, _target in pairs)
    thresholds: list[float] = []
    for mode_n in range(1, int(k_max) + 1):
        unsafe = [value for value, target in pairs if target < mode_n]
        threshold = maximum_x if not unsafe else float(np.nextafter(min(unsafe), -math.inf))
        thresholds.append(threshold)
    if any(right > left for left, right in zip(thresholds, thresholds[1:])):
        raise AssertionError("epsilon threshold vector must be monotone non-increasing")
    return tuple(thresholds)


def epsilon_prediction(
    record: certification.GeometryRecord,
    rule: str,
    thresholds: Sequence[float],
    train: Sequence[certification.GeometryRecord],
    *,
    guarded: bool,
) -> dict[str, object]:
    x_value = predictor_value(record, rule)
    unguarded = max((index for index, threshold in enumerate(thresholds, start=1) if x_value <= threshold), default=0)
    train_x = [predictor_value(item, rule) for item in train]
    outside = not train_x or x_value < min(train_x) - 1.0e-12 or x_value > max(train_x) + 1.0e-12
    n_hat = 0 if guarded and outside else unguarded
    return {
        "N_hat": n_hat,
        "N_hat_guarded": n_hat,
        "N_hat_unguarded": unguarded,
        "out_of_X_calibration_domain": outside,
        "first_predicted_failure": n_hat + 1 if n_hat < len(thresholds) else "",
        "trigger_reason": "out_of_X_calibration_domain" if guarded and outside else ("epsilon_threshold" if n_hat < len(thresholds) else ""),
        "cluster_guard_triggered": False,
        "mixed_mode_guard_triggered": False,
    }


def epsilon_threshold_train_metrics(
    records: Sequence[certification.GeometryRecord],
    rule: str,
    thresholds: Sequence[float],
) -> dict[str, int]:
    predictions = [epsilon_prediction(record, rule, thresholds, records, guarded=False) for record in records]
    false_safe = sum(int(item["N_hat"]) > int(record.N_true or 0) for item, record in zip(predictions, records, strict=True))
    objective = sum(int(item["N_hat"]) for item in predictions)
    return {"false_safe": false_safe, "objective": objective}


def thresholds_json(values: Sequence[float] | Mapping[str, float]) -> str:
    if isinstance(values, Mapping):
        payload = {key: fmt(value) for key, value in values.items()}
    else:
        payload = {f"T_{index}": fmt(value) for index, value in enumerate(values, start=1)}
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def prediction_row(
    fold: PilotFold,
    status: str,
    rule: str,
    record: certification.GeometryRecord,
    detail: Mapping[str, object],
    thresholds: Sequence[float] | Mapping[str, float],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    n_true = int(record.N_true or 0)
    n_hat = int(detail["N_hat"])
    return {
        "split_kind": fold.split_kind,
        "fold_id": fold.fold_id,
        "fold_status": status,
        "evaluation_scope": "optimistic_in_sample" if fold.optimistic_in_sample else "held_out_complete_geometry",
        "rule": rule,
        "case_id": case_id_for(record, manifest_by_case),
        "geometry_id": record.geometry_id,
        "epsilon_0": record.epsilon_0,
        "epsilon_max": epsilon_max(record),
        "beta_deg": record.beta_deg,
        "mu": record.mu,
        "eta": record.eta,
        "N_true": n_true,
        "N_hat": n_hat,
        "N_hat_guarded": detail.get("N_hat_guarded", n_hat),
        "N_hat_unguarded": detail.get("N_hat_unguarded", n_hat),
        "false_safe": n_hat > n_true,
        "false_safe_frequency_count": max(n_hat - n_true, 0),
        "conservative_loss": max(n_true - n_hat, 0),
        "out_of_X_calibration_domain": detail.get("out_of_X_calibration_domain", False),
        "first_predicted_failure": detail.get("first_predicted_failure", detail.get("first_predicted_failed_mode", "")),
        "trigger_reason": detail.get("trigger_reason", ""),
        "cluster_guard_triggered": detail.get("cluster_guard_triggered", False),
        "mixed_mode_guard_triggered": detail.get("mixed_mode_guard_triggered", False),
        "thresholds_json": thresholds_json(thresholds),
        "notes": "same-geometry mode rows remain in one fold",
    }


def run_epsilon_rules(
    folds: Sequence[PilotFold],
    manifest_by_case: Mapping[str, Mapping[str, object]],
    *,
    k_max: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[tuple[str, str, str], tuple[float, ...]]]:
    threshold_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    saved: dict[tuple[str, str, str], tuple[float, ...]] = {}
    baseline_fold = next((fold for fold in folds if fold.split_kind == "baseline_reference_transfer"), None)
    baseline_thresholds: tuple[float, ...] | None = None
    if baseline_fold is not None:
        baseline_thresholds = calibrate_epsilon_thresholds(baseline_fold.train, predictor="epsilon_0", k_max=k_max)

    for fold in folds:
        status = fold_status(fold)
        rules = ["Rule_E0_cal", "Rule_Emax_cal"]
        if fold.split_kind == "baseline_reference_transfer":
            rules = ["Rule_E0_ref", "Rule_Emax_ref", *rules]
        for rule in rules:
            if rule in {"Rule_E0_ref", "Rule_Emax_ref"}:
                if baseline_thresholds is None:
                    continue
                thresholds = baseline_thresholds
            else:
                thresholds = calibrate_epsilon_thresholds(
                    fold.train,
                    predictor=predictor_name(rule),
                    k_max=k_max,
                )
            saved[(fold.split_kind, fold.fold_id, rule)] = thresholds
            train_metrics = epsilon_threshold_train_metrics(fold.train, rule, thresholds)
            if train_metrics["false_safe"] != 0:
                raise AssertionError(f"epsilon calibration produced false-safe for {rule}")
            train_x = [predictor_value(record, rule) for record in fold.train]
            monotone = not any(right > left for left, right in zip(thresholds, thresholds[1:]))
            for mode_n, threshold in enumerate(thresholds, start=1):
                threshold_rows.append(
                    {
                        "split_kind": fold.split_kind,
                        "fold_id": fold.fold_id,
                        "fold_status": status,
                        "rule": rule,
                        "predictor": predictor_name(rule),
                        "mode_n": mode_n,
                        "threshold": threshold,
                        "train_X_min": min(train_x) if train_x else float("nan"),
                        "train_X_max": max(train_x) if train_x else float("nan"),
                        "exact_calibration": True,
                        "monotonicity_status": "monotone_nonincreasing" if monotone else "failed",
                        "train_false_safe_geometry_count": train_metrics["false_safe"],
                        "train_retained_frequency_objective": train_metrics["objective"],
                        "notes": "exact finite-observation threshold; equality is accepted",
                    }
                )
            for record in fold.validation:
                detail = epsilon_prediction(record, rule, thresholds, fold.train, guarded=not fold.optimistic_in_sample)
                prediction_rows.append(
                    prediction_row(fold, status, rule, record, detail, thresholds, manifest_by_case)
                )
    return threshold_rows, prediction_rows, saved


def run_eb_rules(
    folds: Sequence[PilotFold],
    manifest_by_case: Mapping[str, Mapping[str, object]],
    *,
    k_max: int,
    grid_sizes: Sequence[int],
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[tuple[str, str, str], dict[str, object]],
]:
    prediction_rows: list[dict[str, object]] = []
    convergence_rows: list[dict[str, object]] = []
    offline_operation_rows: list[dict[str, object]] = []
    primary: dict[tuple[str, str, str], dict[str, object]] = {}
    primary_grid = max(grid_sizes)
    for fold in folds:
        status = fold_status(fold)
        if not fold.train:
            continue
        for rule in EB_RULES:
            previous_objective: int | None = None
            previous_predictions: dict[str, int] | None = None
            grid_results: dict[int, dict[str, object]] = {}
            for grid_size in grid_sizes:
                operations = certification.OperationCounts()
                result = certification.calibrate_rule(
                    fold.train,
                    rule,
                    k_max=k_max,
                    max_candidates_per_axis=grid_size,
                    operations=operations,
                )
                thresholds = result["thresholds"]
                assert isinstance(thresholds, Mapping)
                current_predictions = {
                    case_id_for(record, manifest_by_case): int(
                        certification.predict_safe_prefix(record.modes, rule, thresholds)["N_hat"]
                    )
                    for record in fold.validation
                }
                objective = int(result["objective_retained_frequencies"])
                changes = (
                    sum(current_predictions.get(key) != previous_predictions.get(key) for key in current_predictions)
                    if previous_predictions is not None
                    else 0
                )
                objective_changed = previous_objective is not None and objective != previous_objective
                converged = not objective_changed and changes == 0 if previous_objective is not None else False
                metrics = result.get("metrics", {})
                assert isinstance(metrics, Mapping)
                convergence_rows.append(
                    {
                        "split_kind": fold.split_kind,
                        "fold_id": fold.fold_id,
                        "fold_status": status,
                        "rule": rule,
                        "max_candidates_per_axis": grid_size,
                        "retained_frequency_objective": objective,
                        "thresholds_json": thresholds_json(thresholds),
                        "calibration_false_safe_geometry_count": metrics.get("false_safe_geometry_count", ""),
                        "candidate_count_per_axis": result.get("candidate_count_per_axis", {}),
                        "evaluated_combinations": result.get("evaluated_combination_count", 0),
                        "search_is_exact": result.get("search_is_exact", False),
                        "prediction_changes_from_previous_grid": changes,
                        "objective_changed_from_previous_grid": objective_changed,
                        "calibration_search_converged": converged,
                        "notes": "convergence compares held-out N_hat and train retained-frequency objective",
                    }
                )
                offline_operation_rows.append(
                    {
                        "cost_scope": "offline_threshold_calibration",
                        "cost_status": "measured_postprocessor_counter",
                        "split_kind": fold.split_kind,
                        "fold_id": fold.fold_id,
                        "rule": rule,
                        "geometry_count": len(fold.train),
                        "frequency_count": len(fold.train) * k_max,
                        "threshold_comparisons": operations.scalar_predictor_comparisons,
                        "scalar_predictor_comparisons": operations.scalar_predictor_comparisons,
                        "threshold_combinations_evaluated": operations.threshold_combinations_evaluated,
                        "Timoshenko_frequencies_required": 0,
                        "Timoshenko_frequencies_avoided": 0,
                        "notes": f"offline candidate grid {grid_size}; not mixed with online inference",
                    }
                )
                grid_results[grid_size] = result
                previous_objective = objective
                previous_predictions = current_predictions
            result = grid_results[primary_grid]
            primary[(fold.split_kind, fold.fold_id, rule)] = result
            thresholds = result["thresholds"]
            assert isinstance(thresholds, Mapping)
            for record in fold.validation:
                detail = certification.predict_safe_prefix(record.modes, rule, thresholds)
                prediction_rows.append(
                    prediction_row(fold, status, rule, record, detail, thresholds, manifest_by_case)
                )
    return prediction_rows, convergence_rows, offline_operation_rows, primary


def prediction_metrics(rows: Sequence[Mapping[str, object]], *, k_max: int) -> dict[str, object]:
    normalized = [
        {
            "N_true": int(float(row["N_true"])),
            "N_hat": int(float(row["N_hat_guarded"])),
            "cluster_guard_triggered": row.get("cluster_guard_triggered", False),
            "mixed_mode_guard_triggered": row.get("mixed_mode_guard_triggered", False),
            "out_of_calibration_domain": row.get("out_of_X_calibration_domain", False),
        }
        for row in rows
    ]
    return certification.prediction_metrics(normalized, k_max=k_max)


def validation_summaries(
    prediction_rows: Sequence[dict[str, object]],
    convergence_rows: Sequence[dict[str, object]],
    *,
    k_max: int,
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    fold_ids: dict[tuple[str, str], set[str]] = defaultdict(set)
    statuses: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in prediction_rows:
        key = (str(row["split_kind"]), str(row["rule"]))
        grouped[key].append(row)
        fold_ids[key].add(str(row["fold_id"]))
        statuses[key].add(str(row["fold_status"]))
        if row["split_kind"] != "combined_in_sample":
            aggregate = ("all_held_out", str(row["rule"]))
            grouped[aggregate].append(row)
            fold_ids[aggregate].add(f"{row['split_kind']}:{row['fold_id']}")
            statuses[aggregate].add(str(row["fold_status"]))
    convergence_lookup: dict[tuple[str, str], list[bool]] = defaultdict(list)
    maximum_grid = max((int(row["max_candidates_per_axis"]) for row in convergence_rows), default=0)
    for row in convergence_rows:
        if int(row["max_candidates_per_axis"]) != maximum_grid:
            continue
        key = (str(row["split_kind"]), str(row["rule"]))
        convergence_lookup[key].append(bool_value(row["calibration_search_converged"]))
        if row["split_kind"] != "combined_in_sample":
            convergence_lookup[("all_held_out", str(row["rule"]))].append(bool_value(row["calibration_search_converged"]))
    rows: list[dict[str, object]] = []
    for key in sorted(grouped):
        split_kind, rule = key
        metrics = prediction_metrics(grouped[key], k_max=k_max)
        convergence = convergence_lookup.get(key, [])
        rows.append(
            {
                "split_kind": split_kind,
                "rule": rule,
                "fold_count": len(fold_ids[key]),
                "fold_statuses": ";".join(sorted(statuses[key])),
                **metrics,
                "out_of_X_calibration_domain_count": metrics.get("out_of_calibration_domain_count", 0),
                "calibration_search_converged": all(convergence) if convergence else (True if rule in E_RULES else ""),
                "notes": "optimistic in-sample diagnostic" if split_kind == "combined_in_sample" else "held-out complete-geometry fold evaluations",
            }
        )
    return rows


def comparison_table(
    records: Sequence[certification.GeometryRecord],
    folds: Sequence[PilotFold],
    epsilon_thresholds: Mapping[tuple[str, str, str], tuple[float, ...]],
    eb_primary: Mapping[tuple[str, str, str], Mapping[str, object]],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    baseline = next(fold for fold in folds if fold.split_kind == "baseline_reference_transfer")
    combined = next(fold for fold in folds if fold.split_kind == "combined_in_sample")
    rows: list[dict[str, object]] = []
    for record in sorted(records, key=lambda item: case_id_for(item, manifest_by_case)):
        row: dict[str, object] = {
            "case_id": case_id_for(record, manifest_by_case),
            "epsilon_0": record.epsilon_0,
            "epsilon_max": epsilon_max(record),
            "beta_deg": record.beta_deg,
            "mu": record.mu,
            "eta": record.eta,
            "N_true": int(record.N_true or 0),
        }
        for rule, field in (
            ("Rule_E0_ref", "N_hat_E0_ref"),
            ("Rule_Emax_ref", "N_hat_Emax_ref"),
        ):
            thresholds = epsilon_thresholds[(baseline.split_kind, baseline.fold_id, rule)]
            row[field] = epsilon_prediction(record, rule, thresholds, baseline.train, guarded=True)["N_hat"]
        for rule, field in (
            ("Rule_E0_cal", "N_hat_E0_cal"),
            ("Rule_Emax_cal", "N_hat_Emax_cal"),
        ):
            thresholds = epsilon_thresholds[(combined.split_kind, combined.fold_id, rule)]
            row[field] = epsilon_prediction(record, rule, thresholds, combined.train, guarded=False)["N_hat"]
        for rule in EB_RULES:
            result = eb_primary[(combined.split_kind, combined.fold_id, rule)]
            thresholds = result["thresholds"]
            assert isinstance(thresholds, Mapping)
            row[f"N_hat_{rule}"] = certification.predict_safe_prefix(record.modes, rule, thresholds)["N_hat"]
        row["notes"] = "reference rules use the baseline X-domain guard; calibrated and EB rules use optimistic combined calibration"
        rows.append(row)
    return rows


def reference_curve_rows(
    records: Sequence[certification.GeometryRecord],
    baseline_fold: PilotFold,
    epsilon_thresholds: Mapping[tuple[str, str, str], tuple[float, ...]],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    thresholds = epsilon_thresholds[(baseline_fold.split_kind, baseline_fold.fold_id, "Rule_E0_ref")]
    ordered = sorted(baseline_fold.train, key=lambda item: item.epsilon_0)
    running_min = 10
    rows: list[dict[str, object]] = []
    for record in ordered:
        running_min = min(running_min, int(record.N_true or 0))
        detail = epsilon_prediction(record, "Rule_E0_ref", thresholds, baseline_fold.train, guarded=False)
        deltas = {int(float(mode["sorted_index"])): finite_float(mode, "delta_f") for mode in record.modes}
        geometry = manifest_by_case[case_id_for(record, manifest_by_case)]
        del geometry
        first_failed = int(record.N_true or 0) + 1 if int(record.N_true or 0) < 10 else ""
        passes = [deltas[index] <= 0.10 for index in range(1, 11)]
        first_failure_position = next((index for index, value in enumerate(passes) if not value), None)
        rows.append(
            {
                "case_id": case_id_for(record, manifest_by_case),
                "epsilon": record.epsilon_0,
                "epsilon_max": epsilon_max(record),
                "N_true_raw": int(record.N_true or 0),
                "N_true_monotone_reference": running_min,
                "reference_N_from_thresholds": detail["N_hat"],
                "first_failed_mode": first_failed,
                "late_pass_after_failure": bool(first_failure_position is not None and any(passes[first_failure_position + 1 :])),
                **{f"delta_f_{index}": deltas[index] for index in range(1, 11)},
                "thresholds_json": thresholds_json(thresholds),
                "notes": "monotone reference is the cumulative minimum of observed N_true",
            }
        )
    return rows


def matched_group_rows(
    records: Sequence[certification.GeometryRecord],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for group, target in (("matched_epsilon_max_0p02", 0.02), ("matched_epsilon_max_0p04", 0.04)):
        members = [record for record in records if group in groups_for(record, manifest_by_case)]
        members.sort(key=lambda item: case_id_for(item, manifest_by_case))
        if not members:
            continue
        n_values = [int(record.N_true or 0) for record in members]
        spreads: dict[int, float] = {}
        for index in range(1, 11):
            values = [finite_float(record.modes[index - 1], "delta_f") for record in members]
            spreads[index] = max(values) - min(values)
        max_mode = max(spreads, key=spreads.get)
        first_failed = [value + 1 if value < 10 else 0 for value in n_values]
        reordering = any(
            int(float(mode.get("matched_timo_sorted_index", mode["sorted_index"]))) != int(float(mode["sorted_index"]))
            for record in members
            for mode in record.modes
            if str(mode.get("matched_timo_sorted_index", "")).strip()
        )
        cluster = any(
            str(mode.get("cluster_id", "")).strip()
            or str(mode.get("matching_status", "")) in {"ambiguous", "reliable_cluster"}
            for record in members
            for mode in record.modes
        )
        rows.append(
            {
                "group": group,
                "epsilon_max_target": target,
                "case_ids": ",".join(case_id_for(record, manifest_by_case) for record in members),
                "N_true_values": ",".join(str(value) for value in n_values),
                "range_N_true": max(n_values) - min(n_values),
                "std_N_true": statistics.pstdev(n_values),
                "maximum_per_mode_delta_f_spread": spreads[max_mode],
                "mode_of_maximum_spread": max_mode,
                "first_failed_mode_identical": len(set(first_failed)) == 1,
                "first_failed_modes": ",".join(str(value or "none") for value in first_failed),
                "reordering_warning": reordering,
                "cluster_warning": cluster,
                "notes": "equal epsilon_max is a local-thickness diagnostic and does not imply equal spectra",
            }
        )
    return rows


def rank_correlation(x_values: Sequence[float], y_values: Sequence[float]) -> float:
    if len(x_values) < 2 or len(set(x_values)) < 2 or len(set(y_values)) < 2:
        return float("nan")
    x_ranks = np.argsort(np.argsort(np.asarray(x_values, dtype=float))).astype(float)
    y_ranks = np.argsort(np.argsort(np.asarray(y_values, dtype=float))).astype(float)
    return float(np.corrcoef(x_ranks, y_ranks)[0, 1])


def same_epsilon_rows(
    records: Sequence[certification.GeometryRecord],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    members = [record for record in records if "fixed_epsilon_0p02_block" in groups_for(record, manifest_by_case)]
    members.sort(key=lambda item: case_id_for(item, manifest_by_case))
    n_values = [int(record.N_true or 0) for record in members]
    eps_values = [epsilon_max(record) for record in members]
    variation = max(n_values) - min(n_values) if n_values else 0
    correlation = rank_correlation(eps_values, n_values)
    return [
        {
            "case_id": case_id_for(record, manifest_by_case),
            "epsilon_0": record.epsilon_0,
            "epsilon_max": epsilon_max(record),
            "beta_deg": record.beta_deg,
            "mu": record.mu,
            "eta": record.eta,
            "N_true": int(record.N_true or 0),
            "first_failed_mode": int(record.N_true or 0) + 1 if int(record.N_true or 0) < 10 else "",
            "group_range_N_true": variation,
            "group_std_N_true": statistics.pstdev(n_values) if n_values else float("nan"),
            "epsilon_max_N_true_rank_correlation": correlation,
            "same_epsilon_material_variation": variation > 0,
            "notes": "beta and mu vary while epsilon_0 remains 0.02",
        }
        for record in members
    ]


def online_comparisons(n_hat_values: Sequence[int], k_max: int) -> int:
    return sum(k_max if value >= k_max else value + 1 for value in n_hat_values)


def operation_rows(
    prediction_rows: Sequence[dict[str, object]],
    offline_rows: Sequence[dict[str, object]],
    root_rows: Sequence[Mapping[str, object]],
    *,
    k_max: int,
) -> list[dict[str, object]]:
    rows = list(offline_rows)
    rows.append(
        {
            "cost_scope": "offline_reference_generation",
            "cost_status": "measured_existing_wrapper_counters",
            "geometry_count": len(root_rows),
            "frequency_count": len(root_rows) * k_max,
            "EB_root_operations_if_available": sum(int(float(row.get("EB_root_solver_calls", 0) or 0)) for row in root_rows),
            "EB_characteristic_matrix_evaluations": "unavailable",
            "EB_determinant_evaluations": "unavailable",
            "Timoshenko_characteristic_matrix_evaluations": "unavailable",
            "Timoshenko_determinant_evaluations": "unavailable",
            "root_bracketing_intervals": "unavailable",
            "Brent_function_evaluations": "unavailable",
            "source_root_cache_hits": sum(int(float(row.get("root_cache_hits", 0) or 0)) for row in root_rows),
            "source_root_cache_misses": sum(int(float(row.get("root_cache_misses", 0) or 0)) for row in root_rows),
            "source_boundary_retries": sum(int(float(row.get("boundary_retries", 0) or 0)) for row in root_rows),
            "Timoshenko_frequencies_required": len(root_rows) * k_max,
            "Timoshenko_frequencies_avoided": 0,
            "notes": "reference target generation; determinant-level counters remain unavailable in the source wrapper",
        }
    )
    grouped: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in prediction_rows:
        grouped[(str(row["split_kind"]), str(row["fold_id"]), str(row["rule"]))].append(row)
    n_shape_points = int(float(root_rows[0].get("n_shape_points", 0))) if root_rows else 0
    for (split_kind, fold_id, rule), items in sorted(grouped.items()):
        n_hats = [int(float(row["N_hat_guarded"])) for row in items]
        geometry_count = len(items)
        threshold_count = online_comparisons(n_hats, k_max)
        timo_required = sum(k_max - value for value in n_hats)
        if rule in E_RULES:
            emax = rule in {"Rule_Emax_ref", "Rule_Emax_cal"}
            out_count = sum(bool_value(row.get("out_of_X_calibration_domain", False)) for row in items)
            rows.append(
                {
                    "cost_scope": "online_geometry_only_inference",
                    "cost_status": "exact_algorithmic_primitive_count",
                    "split_kind": split_kind,
                    "fold_id": fold_id,
                    "rule": rule,
                    "geometry_count": geometry_count,
                    "frequency_count": geometry_count * k_max,
                    "geometry_rule_X_evaluations": geometry_count,
                    "tau_factor_evaluations": geometry_count if emax else 0,
                    "local_epsilon_evaluations": 2 * geometry_count if emax else 0,
                    "max_operations": geometry_count if emax else 0,
                    "threshold_comparisons": max(threshold_count - out_count, 0),
                    "EB_root_operations_if_available": 0,
                    "EB_characteristic_matrix_evaluations": 0,
                    "EB_determinant_evaluations": 0,
                    "Timoshenko_characteristic_matrix_evaluations": 0,
                    "Timoshenko_determinant_evaluations": 0,
                    "root_bracketing_intervals": 0,
                    "Brent_function_evaluations": 0,
                    "EB_mode_reconstructions": 0,
                    "SVD_6x6_calls": 0,
                    "EB_shape_grid_points": 0,
                    "quadrature_calls": 0,
                    "quadrature_point_evaluations": 0,
                    "scalar_predictor_comparisons": max(threshold_count - out_count, 0),
                    "spectral_gap_operations": 0,
                    "local_Timoshenko_root_attempts": 0,
                    "full_Timoshenko_fallbacks": 0,
                    "Timoshenko_frequencies_required": timo_required,
                    "Timoshenko_frequencies_avoided": sum(n_hats),
                    "notes": "domain checks precede threshold comparisons; no eigenvalue or shape solve",
                }
            )
        else:
            evaluated_modes = online_comparisons(n_hats, k_max)
            gap_rule = rule in {"Rule_A_gap", "Rule_C", "Rule_D"}
            rows.append(
                {
                    "cost_scope": "online_EB_based_inference",
                    "cost_status": "structural_primitive_count; EB root internals unavailable",
                    "split_kind": split_kind,
                    "fold_id": fold_id,
                    "rule": rule,
                    "geometry_count": geometry_count,
                    "frequency_count": geometry_count * k_max,
                    "geometry_rule_X_evaluations": 0,
                    "tau_factor_evaluations": 0,
                    "local_epsilon_evaluations": 0,
                    "max_operations": 0,
                    "threshold_comparisons": threshold_count * (2 if rule in {"Rule_B", "Rule_C", "Rule_D"} else 1),
                    "EB_root_operations_if_available": "unavailable",
                    "EB_characteristic_matrix_evaluations": "unavailable",
                    "EB_determinant_evaluations": "unavailable",
                    "Timoshenko_characteristic_matrix_evaluations": 0,
                    "Timoshenko_determinant_evaluations": 0,
                    "root_bracketing_intervals": "unavailable",
                    "Brent_function_evaluations": "unavailable",
                    "EB_mode_reconstructions": evaluated_modes,
                    "SVD_6x6_calls": evaluated_modes,
                    "EB_shape_grid_points": 2 * n_shape_points * evaluated_modes,
                    "quadrature_calls": 12 * evaluated_modes,
                    "quadrature_point_evaluations": 12 * n_shape_points * evaluated_modes,
                    "scalar_predictor_comparisons": threshold_count * (2 if rule in {"Rule_B", "Rule_C", "Rule_D"} else 1),
                    "spectral_gap_operations": geometry_count * k_max if gap_rule else 0,
                    "local_Timoshenko_root_attempts": 0,
                    "full_Timoshenko_fallbacks": 0,
                    "Timoshenko_frequencies_required": timo_required,
                    "Timoshenko_frequencies_avoided": sum(n_hats),
                    "notes": "Timoshenko fallback count is structural and not a FLOP estimate",
                }
            )
    return rows


def plot_n_true(
    output: Path,
    comparison_rows: Sequence[Mapping[str, object]],
    *,
    x_field: str,
    filename: str,
) -> Path:
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    colors = np.asarray([finite_float(row, "beta_deg") for row in comparison_rows], dtype=float)
    scatter = ax.scatter(
        [finite_float(row, x_field) for row in comparison_rows],
        [finite_float(row, "N_true") for row in comparison_rows],
        c=colors,
        cmap="viridis",
        edgecolors="black",
        linewidths=0.5,
    )
    for row in comparison_rows:
        ax.annotate(str(row["case_id"]), (finite_float(row, x_field), finite_float(row, "N_true")), xytext=(3, 3), textcoords="offset points", fontsize=7)
    ax.set_xlabel(x_field)
    ax.set_ylabel("N_true")
    ax.set_ylim(-0.4, 10.6)
    ax.grid(True, alpha=0.25)
    fig.colorbar(scatter, ax=ax, label="beta [deg]")
    fig.tight_layout()
    path = output / filename
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def generate_plots(
    output: Path,
    comparison_rows: Sequence[Mapping[str, object]],
    reference_rows: Sequence[Mapping[str, object]],
    matched_rows: Sequence[Mapping[str, object]],
    summary_rows: Sequence[Mapping[str, object]],
    records: Sequence[certification.GeometryRecord],
    manifest_by_case: Mapping[str, Mapping[str, object]],
) -> list[Path]:
    paths = [
        plot_n_true(output, comparison_rows, x_field="epsilon_0", filename="N_true_vs_epsilon0.png"),
        plot_n_true(output, comparison_rows, x_field="epsilon_max", filename="N_true_vs_epsilon_max.png"),
    ]
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for index in range(1, 11):
        ax.plot(
            [finite_float(row, "epsilon") for row in reference_rows],
            [finite_float(row, f"delta_f_{index}") for row in reference_rows],
            marker="o",
            linewidth=1.0,
            markersize=3,
            label=f"k={index}",
        )
    ax.axhline(0.10, color="black", linestyle="--", linewidth=1.0)
    ax.set_xlabel("epsilon")
    ax.set_ylabel("delta_f")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=2, fontsize=7)
    fig.tight_layout()
    path = output / "baseline_delta_f_vs_epsilon.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), sharey=True)
    for axis, group in zip(axes, ("matched_epsilon_max_0p02", "matched_epsilon_max_0p04"), strict=True):
        members = [record for record in records if group in groups_for(record, manifest_by_case)]
        for record in members:
            axis.plot(
                range(1, 11),
                [finite_float(mode, "delta_f") for mode in record.modes],
                marker="o",
                markersize=3,
                label=case_id_for(record, manifest_by_case),
            )
        axis.axhline(0.10, color="black", linestyle="--", linewidth=1.0)
        axis.set_title(group.replace("matched_epsilon_max_", "epsilon_max="))
        axis.set_xlabel("sorted mode")
        axis.grid(True, alpha=0.25)
        if members:
            axis.legend(fontsize=7)
    axes[0].set_ylabel("delta_f")
    fig.tight_layout()
    path = output / "matched_epsilon_max_groups.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    baseline_summary = [row for row in summary_rows if row.get("split_kind") == "baseline_reference_transfer"]
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for row in baseline_summary:
        x = finite_float(row, "false_safe_geometry_count")
        y = finite_float(row, "usable_frequency_retention")
        ax.scatter(x, y, s=45)
        ax.annotate(str(row["rule"]), (x, y), xytext=(4, 3), textcoords="offset points", fontsize=8)
    ax.set_xlabel("false-safe geometry count")
    ax.set_ylabel("usable EB frequency retention")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    path = output / "rule_false_safe_retention_comparison.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(11.0, 5.2))
    positions = np.arange(len(comparison_rows))
    ax.plot(positions, [finite_float(row, "N_true") for row in comparison_rows], color="black", marker="o", linewidth=1.8, label="N_true")
    for field, label in (("N_hat_Emax_ref", "Emax-ref"), ("N_hat_Emax_cal", "Emax-cal"), ("N_hat_Rule_A_gap", "Rule A-gap"), ("N_hat_Rule_D", "Rule D")):
        ax.plot(positions, [finite_float(row, field) for row in comparison_rows], marker=".", linewidth=1.0, label=label)
    ax.set_xticks(positions, [str(row["case_id"]) for row in comparison_rows], rotation=60)
    ax.set_ylabel("safe-prefix length")
    ax.set_ylim(-0.4, 10.6)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(ncol=3, fontsize=8)
    fig.tight_layout()
    path = output / "geometry_rule_vs_EB_rule_prefixes.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)
    del matched_rows
    return paths


def summary_lookup(rows: Sequence[Mapping[str, object]], split: str, rule: str) -> Mapping[str, object] | None:
    return next((row for row in rows if row.get("split_kind") == split and row.get("rule") == rule), None)


def metric_text(row: Mapping[str, object] | None, key: str, digits: int = 4) -> str:
    if row is None:
        return "n/a"
    value = finite_float(row, key)
    return f"{value:.{digits}g}" if isfinite(value) else "n/a"


def write_report(
    args: Args,
    records: Sequence[certification.GeometryRecord],
    comparison_rows: Sequence[Mapping[str, object]],
    reference_rows: Sequence[Mapping[str, object]],
    matched_rows: Sequence[Mapping[str, object]],
    same_rows: Sequence[Mapping[str, object]],
    summary_rows: Sequence[Mapping[str, object]],
    convergence_rows: Sequence[Mapping[str, object]],
) -> Path:
    report = args.output_dir / "epsilon_apriori_pilot_report.md"
    baseline_rows = [row for row in summary_rows if row.get("split_kind") == "baseline_reference_transfer"]
    zero_safe = [row for row in baseline_rows if finite_float(row, "false_safe_geometry_count") == 0.0]
    best = sorted(
        zero_safe or baseline_rows,
        key=lambda row: (
            finite_float(row, "false_safe_geometry_count"),
            -finite_float(row, "usable_frequency_retention"),
            finite_float(row, "mean_conservative_loss"),
        ),
    )[0] if baseline_rows else None
    e0 = summary_lookup(summary_rows, "baseline_reference_transfer", "Rule_E0_ref")
    emax = summary_lookup(summary_rows, "baseline_reference_transfer", "Rule_Emax_ref")
    emax_better = (
        emax is not None
        and e0 is not None
        and (
            finite_float(emax, "false_safe_geometry_count"),
            -finite_float(emax, "usable_frequency_retention"),
        )
        < (
            finite_float(e0, "false_safe_geometry_count"),
            -finite_float(e0, "usable_frequency_retention"),
        )
    )
    convergence32 = [row for row in convergence_rows if int(float(row["max_candidates_per_axis"])) == max(args.candidate_grid_sizes)]
    all_converged = all(bool_value(row["calibration_search_converged"]) for row in convergence32) if convergence32 else False
    lines = [
        "# EB Epsilon A-Priori Pilot Report",
        "",
        "## Scope",
        "",
        f"This diagnostic uses `K = {args.k_max}`, a `10%` dimensional-frequency discrepancy threshold, and `{len(records)}` included manifest geometries. The target is the first continuous prefix of sorted frequencies. Timoshenko is the one-dimensional reference; no FEM was rerun, and no analytic formula, determinant, root solver, tolerance, or shear coefficient was modified.",
        "",
        "## Pilot Cases",
        "",
        "| case | epsilon_0 | epsilon_max | beta | mu | eta | N_true |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in comparison_rows:
        lines.append(
            f"| {row['case_id']} | {finite_float(row, 'epsilon_0'):.4g} | {finite_float(row, 'epsilon_max'):.4g} | {finite_float(row, 'beta_deg'):.4g} | {finite_float(row, 'mu'):.4g} | {finite_float(row, 'eta'):.4g} | {int(finite_float(row, 'N_true'))} |"
        )
    lines.extend([
        "",
        "## Baseline Reference Line",
        "",
        "| epsilon | N_true | monotone reference | first failed | late pass |",
        "|---:|---:|---:|---:|---|",
    ])
    for row in reference_rows:
        lines.append(
            f"| {finite_float(row, 'epsilon'):.4g} | {int(finite_float(row, 'N_true_raw'))} | {int(finite_float(row, 'N_true_monotone_reference'))} | {row.get('first_failed_mode', '')} | {row.get('late_pass_after_failure', '')} |"
        )
    lines.extend([
        "",
        "The monotone reference is a conservative envelope of the seven observed baseline cases. Any non-monotonic raw prefix or late mode re-entry remains visible in the CSV and is not used to repair an earlier sorted-mode failure.",
        "",
        "## Direct Epsilon Comparison",
        "",
        "| rule | false-safe geometries | retention | mean loss | X-domain fallbacks |",
        "|---|---:|---:|---:|---:|",
    ])
    for rule in E_RULES:
        row = summary_lookup(summary_rows, "baseline_reference_transfer", rule)
        if row is None:
            continue
        lines.append(
            f"| {rule} | {metric_text(row, 'false_safe_geometry_count')} | {metric_text(row, 'usable_frequency_retention')} | {metric_text(row, 'mean_conservative_loss')} | {metric_text(row, 'out_of_calibration_domain_count')} |"
        )
    lines.extend([
        "",
        f"On this selected transfer set, `{'epsilon_max' if emax_better else 'epsilon_0'}` is the better of the two baseline-reference indicators under the false-safe-first ordering. This is a pilot observation, not a universal criterion.",
        "",
        "## Matched Epsilon-Max Groups",
        "",
        "| group | N_true | range | max mode-level delta spread | mode | first failure identical | reorder | cluster |",
        "|---|---|---:|---:|---:|---|---|---|",
    ])
    for row in matched_rows:
        lines.append(
            f"| {row['group']} | {row['N_true_values']} | {row['range_N_true']} | {metric_text(row, 'maximum_per_mode_delta_f_spread')} | {row['mode_of_maximum_spread']} | {row['first_failed_mode_identical']} | {row['reordering_warning']} | {row['cluster_warning']} |"
        )
    same_range = int(finite_float(same_rows[0], "group_range_N_true")) if same_rows else 0
    same_corr = metric_text(same_rows[0] if same_rows else None, "epsilon_max_N_true_rank_correlation")
    lines.extend([
        "",
        "Equal `epsilon_max` does not imply equal individual discrepancies or equal spectra. The matched groups test only whether local maximum thickness is a useful coarse ordering variable.",
        "",
        "## Same Epsilon-0 Block",
        "",
        f"At `epsilon_0 = 0.02`, the observed `N_true` range is `{same_range}` modes. The diagnostic rank correlation between `epsilon_max` and `N_true` is `{same_corr}`. Residual differences among equal or nearby `epsilon_max` cases expose beta and modal-order effects that a scalar local thickness cannot represent.",
        "",
        "## Comparison With Existing Rules",
        "",
        "Baseline-reference transfer is the common direct comparison below.",
        "",
        "| rule | false-safe | retention | mean loss | converged |",
        "|---|---:|---:|---:|---|",
    ])
    for rule in ALL_RULES:
        row = summary_lookup(summary_rows, "baseline_reference_transfer", rule)
        if row is None:
            continue
        lines.append(
            f"| {rule} | {metric_text(row, 'false_safe_geometry_count')} | {metric_text(row, 'usable_frequency_retention')} | {metric_text(row, 'mean_conservative_loss')} | {row.get('calibration_search_converged', '')} |"
        )
    lines.extend([
        "",
        f"The false-safe-first ranking selects `{best.get('rule') if best else 'none'}` on this finite transfer set. Candidate-grid convergence between the two largest grids is `{'satisfied for every audited fold/rule' if all_converged else 'not satisfied for every fold/rule; rankings require a convergence caveat'}`.",
        "",
        "## Operation Counts",
        "",
        "Geometry-only online evaluation uses no eigenvalue solve, shape reconstruction, SVD, or quadrature. `epsilon_0` needs one X read plus threshold comparisons; `epsilon_max` additionally needs tau factors, two local epsilon evaluations, and one maximum. EB rules require EB roots and modal predictor construction before their comparisons. The Timoshenko fallback count is `K - N_hat`; it is a structural frequency count, not a FLOP estimate. Offline reference generation and threshold-search operations are reported in separate scopes.",
        "",
        "## Limitations",
        "",
        "- Only 21 deliberately selected geometries were requested; the set is not space-filling.",
        "- Zero observed false-safe is not a mathematical guarantee on the continuous parameter domain.",
        "- The beta=0 baseline-transfer hypothesis remains unproved.",
        "- Equal epsilon_max does not imply equal spectra.",
        "- Timoshenko and existing 3D validation remain diagnostic rather than exact 3D truth.",
        "- The production K=10 Stage-1 and fixed-epsilon grids were not run.",
        "",
        "## Recommendation",
        "",
        f"The geometry-only result does `{'support' if best and str(best.get('rule', '')).startswith('Rule_E') else 'not support'}` replacing the EB-based certificate on this pilot. Retain a cascade in which geometry-only screening may provide a cheap early fallback, followed by EB frequencies/gaps, Pi-based guards, and Timoshenko for the uncertified suffix. A calibrated geometry-only threshold should be promoted only after zero observed false-safe and useful retention transfer to a larger held-out design.",
        "",
        "The next 10-30 points should concentrate on: intermediate beta values inside both matched-epsilon_max groups; eta sign pairs at epsilon_max near 0.02 and 0.04; mu values between 0.5 and 0.7 at fixed epsilon_0; and epsilon values immediately around baseline prefix drops. These points directly test the observed failure structure without launching a full tensor grid.",
    ])
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    ensure_output_policy(args)
    records, manifest_by_case, _modes_by_case, root_rows = load_pilot_records(args)
    if not records:
        raise RuntimeError("no included complete K=10 pilot geometries are available")
    folds = build_folds(records, manifest_by_case)
    threshold_rows, epsilon_predictions, epsilon_thresholds = run_epsilon_rules(
        folds,
        manifest_by_case,
        k_max=args.k_max,
    )
    eb_predictions, convergence_rows, offline_rows, eb_primary = run_eb_rules(
        folds,
        manifest_by_case,
        k_max=args.k_max,
        grid_sizes=args.candidate_grid_sizes,
    )
    prediction_rows = [*epsilon_predictions, *eb_predictions]
    summary_rows = validation_summaries(
        prediction_rows,
        convergence_rows,
        k_max=args.k_max,
    )
    comparison_rows = comparison_table(
        records,
        folds,
        epsilon_thresholds,
        eb_primary,
        manifest_by_case,
    )
    baseline_fold = next(fold for fold in folds if fold.split_kind == "baseline_reference_transfer")
    reference_rows = reference_curve_rows(records, baseline_fold, epsilon_thresholds, manifest_by_case)
    matched_rows = matched_group_rows(records, manifest_by_case)
    same_rows = same_epsilon_rows(records, manifest_by_case)
    cost_rows = operation_rows(
        prediction_rows,
        offline_rows,
        root_rows,
        k_max=args.k_max,
    )

    write_csv(args.output_dir / "epsilon_rule_thresholds.csv", threshold_rows, THRESHOLD_FIELDS)
    write_csv(args.output_dir / "epsilon_rule_fold_predictions.csv", prediction_rows, PREDICTION_FIELDS)
    write_csv(args.output_dir / "epsilon_rule_validation_summary.csv", summary_rows, SUMMARY_FIELDS)
    write_csv(args.output_dir / "epsilon_reference_curve.csv", reference_rows, REFERENCE_FIELDS)
    write_csv(args.output_dir / "matched_epsilon_max_group_audit.csv", matched_rows, MATCHED_FIELDS)
    write_csv(args.output_dir / "same_epsilon0_group_audit.csv", same_rows, SAME_EPSILON_FIELDS)
    write_csv(args.output_dir / "epsilon_vs_existing_rules_geometry_comparison.csv", comparison_rows, COMPARISON_FIELDS)
    write_csv(args.output_dir / "epsilon_rule_operation_counts.csv", cost_rows, OPERATION_FIELDS)
    write_csv(args.output_dir / "epsilon_calibration_search_convergence.csv", convergence_rows, CONVERGENCE_FIELDS)
    plot_paths = generate_plots(
        args.output_dir,
        comparison_rows,
        reference_rows,
        matched_rows,
        summary_rows,
        records,
        manifest_by_case,
    )
    report = write_report(
        args,
        records,
        comparison_rows,
        reference_rows,
        matched_rows,
        same_rows,
        summary_rows,
        convergence_rows,
    )
    print(f"included pilot geometries: {len(records)}")
    print(f"fold predictions: {len(prediction_rows)}")
    print(f"output: {args.output_dir}")
    return {
        "args": args,
        "records": records,
        "folds": folds,
        "threshold_rows": threshold_rows,
        "prediction_rows": prediction_rows,
        "summary_rows": summary_rows,
        "reference_rows": reference_rows,
        "matched_rows": matched_rows,
        "same_rows": same_rows,
        "comparison_rows": comparison_rows,
        "operation_rows": cost_rows,
        "convergence_rows": convergence_rows,
        "plot_paths": plot_paths,
        "report": report,
    }


if __name__ == "__main__":
    main()
