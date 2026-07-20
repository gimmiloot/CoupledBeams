from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field
import hashlib
import itertools
import json
import math
from math import isfinite
from pathlib import Path
import sys
from typing import Mapping, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_fixed_epsilon_geometry_scan as fixed_scan,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_vs_timoshenko_stage1 as stage1,
)


DEFAULT_STAGE1_DIR = Path("results") / "eb_validity_vs_timoshenko_stage1"
DEFAULT_FIXED_EPSILON_DIR = Path("results") / "eb_validity_fixed_epsilon_geometry_scan"
DEFAULT_OUTPUT_DIR = Path("results") / "eb_safe_prefix_certification"
DEFAULT_K_MAX = 10
DEFAULT_FREQUENCY_ERROR_THRESHOLD = 0.10
DEFAULT_N_SHAPE_POINTS = fixed_scan.DEFAULT_SHAPE_POINTS
DEFAULT_MAX_CANDIDATES_PER_AXIS = 16
GEOMETRY_ROUND_DIGITS = 12
DELTA_ABS_TOLERANCE = 1.0e-10
DELTA_REL_TOLERANCE = 1.0e-8
PREDICTOR_CONSISTENCY_TOLERANCE = 1.0e-8

SOURCE_STAGE1 = "stage1"
SOURCE_FIXED = "fixed_epsilon"
RULES = ("Rule_A", "Rule_A_gap", "Rule_B", "Rule_C", "Rule_D")

MODE_AUDIT_FIELDS = [
    "source_dataset",
    "geometry_id",
    "geometry_key",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "sorted_index",
    "Lambda_EB",
    "Lambda_Timo",
    "stored_delta_f",
    "delta_f",
    "true_pass_10",
    "epsilon_max",
    "diameter_to_length_1",
    "diameter_to_length_2",
    "chi_max_EB",
    "chi_eff_EB",
    "Theta_max_EB",
    "Pi_shear_EB",
    "Pi_rotary_EB",
    "Pi_EB",
    "EB_axial_fraction",
    "EB_bending_fraction",
    "EB_classification",
    "gap_left_EB",
    "gap_right_EB",
    "gap_EB",
    "gap_sidedness",
    "mode10_upper_guard_source",
    "source_warning",
    "included_in_certification",
    "exclusion_reason",
    "source_notes",
]

GEOMETRY_FIELDS = [
    "source_dataset",
    "geometry_id",
    "geometry_key",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "N_true",
    "first_true_failed_mode",
    "target_mode_count",
    "complete_K10",
    "data_quality_status",
    "epsilon_1",
    "epsilon_2",
    "epsilon_max",
    "diameter_to_length_1",
    "diameter_to_length_2",
    "max_diameter_to_length",
    "thin_rod_criterion_pass",
    "max_Pi_EB_true_prefix",
    "max_Pi_shear_EB_true_prefix",
    "max_Pi_rotary_EB_true_prefix",
    "max_chi_max_EB_true_prefix",
    "max_chi_eff_EB_true_prefix",
    "max_Theta_max_EB_true_prefix",
    "minimum_gap_EB_true_prefix",
    "minimum_gap_EB_firstK",
    "mean_gap_EB_firstK",
    "mode10_gap_sidedness",
    "explicit_root_warning",
    "source_warnings",
    "exclusion_reason",
    "notes",
]

CALIBRATION_FIELDS = [
    "split_kind",
    "fold_id",
    "fold_status",
    "rule",
    "train_count",
    "T_pi",
    "T_shear",
    "T_rotary",
    "gap_min",
    "thresholds_json",
    "objective_retained_frequencies",
    "calibration_geometry_count",
    "calibration_false_safe_geometry_count",
    "calibration_false_safe_frequency_count",
    "mean_N_true",
    "mean_N_hat",
    "exact_N_match_rate",
    "mean_conservative_loss",
    "max_conservative_loss",
    "usable_frequency_retention",
    "search_is_exact",
    "candidate_count_per_axis",
    "evaluated_combination_count",
    "tie_break_notes",
    "overlap_geometry_count_excluded",
    "notes",
]

PREDICTION_FIELDS = [
    "split_kind",
    "fold_id",
    "fold_status",
    "rule",
    "source_dataset",
    "geometry_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "N_true",
    "N_hat",
    "false_safe",
    "false_safe_frequency_count",
    "conservative_loss",
    "first_predicted_failed_mode",
    "first_trigger",
    "trigger_reason",
    "cluster_guard_triggered",
    "mixed_mode_guard_triggered",
    "out_of_calibration_domain",
    "thresholds_json",
    "notes",
]

VALIDATION_FIELDS = [
    "split_kind",
    "rule",
    "fold_count",
    "fold_statuses",
    "validation_geometry_count",
    "validation_frequency_count",
    "false_safe_geometry_count",
    "false_safe_frequency_count",
    "false_safe_geometry_rate",
    "false_safe_frequency_rate",
    "mean_N_true",
    "mean_N_hat",
    "exact_N_match_count",
    "exact_N_match_rate",
    "mean_conservative_loss",
    "median_conservative_loss",
    "max_conservative_loss",
    "usable_EB_frequency_count",
    "retained_EB_frequency_count",
    "usable_frequency_retention",
    "accepted_geometry_count",
    "fully_EB_geometry_count",
    "zero_EB_geometry_count",
    "cluster_guard_trigger_count",
    "mixed_mode_guard_trigger_count",
    "incomplete_geometry_count",
    "explicit_warning_geometry_count",
    "out_of_calibration_domain_count",
    "notes",
]

EXCLUSION_FIELDS = [
    "source_dataset",
    "geometry_id",
    "geometry_key",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "observed_sorted_indices",
    "target_mode_count",
    "complete_K10",
    "explicit_warning",
    "exclusion_reason",
    "source_file",
    "notes",
]

OVERLAP_FIELDS = [
    "geometry_id",
    "geometry_key",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "stage1_included",
    "fixed_epsilon_included",
    "max_abs_Lambda_EB_difference",
    "max_abs_Lambda_Timo_difference",
    "max_abs_delta_f_difference",
    "N_true_stage1",
    "N_true_fixed_epsilon",
    "N_true_difference",
    "max_abs_Pi_EB_difference",
    "max_abs_chi_max_EB_difference",
    "duplicate_use_policy",
    "notes",
]

CONSISTENCY_FIELDS = [
    "source_dataset",
    "quantity",
    "compared_row_count",
    "max_abs_difference",
    "max_rel_difference",
    "mismatch_count",
    "tolerance",
    "status",
    "notes",
]

OPERATION_FIELDS = [
    "scope",
    "split_kind",
    "fold_id",
    "rule",
    "cost_component",
    "measurement_status",
    "geometry_count",
    "frequency_count",
    "EB_characteristic_matrix_evaluations",
    "EB_determinant_evaluations",
    "Timoshenko_characteristic_matrix_evaluations",
    "Timoshenko_determinant_evaluations",
    "SVD_6x6_calls",
    "root_bracketing_intervals",
    "Brent_function_evaluations",
    "EB_mode_reconstructions",
    "EB_shape_grid_points",
    "quadrature_calls",
    "quadrature_point_evaluations",
    "scalar_predictor_comparisons",
    "spectral_gap_operations",
    "threshold_combinations_evaluated",
    "local_Timoshenko_root_attempts",
    "full_Timoshenko_fallbacks",
    "EB_frequencies_accepted",
    "Timoshenko_frequencies_required",
    "Timoshenko_frequencies_avoided",
    "fraction_firstK_avoided_in_Timoshenko",
    "notes",
]


@dataclass(frozen=True)
class Args:
    stage1_dir: Path
    fixed_epsilon_dir: Path
    output_dir: Path
    k_max: int
    frequency_error_threshold: float
    n_shape_points: int
    max_candidates_per_axis: int
    strict_complete: bool


@dataclass
class OperationCounts:
    EB_characteristic_matrix_evaluations: int = 0
    EB_determinant_evaluations: int = 0
    Timoshenko_characteristic_matrix_evaluations: int = 0
    Timoshenko_determinant_evaluations: int = 0
    SVD_6x6_calls: int = 0
    root_bracketing_intervals: int = 0
    Brent_function_evaluations: int = 0
    EB_mode_reconstructions: int = 0
    EB_shape_grid_points: int = 0
    quadrature_calls: int = 0
    quadrature_point_evaluations: int = 0
    scalar_predictor_comparisons: int = 0
    spectral_gap_operations: int = 0
    threshold_combinations_evaluated: int = 0
    local_Timoshenko_root_attempts: int = 0
    full_Timoshenko_fallbacks: int = 0

    def as_dict(self) -> dict[str, int]:
        return {name: int(getattr(self, name)) for name in self.__dataclass_fields__}


@dataclass
class GeometryRecord:
    source_dataset: str
    geometry_id: str
    geometry_key: tuple[float, float, float, float]
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float
    modes: list[dict[str, object]]
    source_file: Path
    complete: bool = False
    included: bool = False
    explicit_warning: bool = False
    exclusion_reasons: list[str] = field(default_factory=list)
    N_true: int | None = None


@dataclass(frozen=True)
class Fold:
    split_kind: str
    fold_id: str
    train: tuple[GeometryRecord, ...]
    validation: tuple[GeometryRecord, ...]
    overlap_excluded: int = 0
    optimistic_in_sample: bool = False


def repo_path(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Calibrate and validate conservative EB-only safe-prefix rules from existing "
            "sorted EB/Timoshenko applicability CSV files. No roots are solved."
        ),
    )
    parser.add_argument("--stage1-dir", type=Path, default=DEFAULT_STAGE1_DIR)
    parser.add_argument("--fixed-epsilon-dir", type=Path, default=DEFAULT_FIXED_EPSILON_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--k-max", type=int, default=DEFAULT_K_MAX)
    parser.add_argument(
        "--frequency-error-threshold",
        type=float,
        default=DEFAULT_FREQUENCY_ERROR_THRESHOLD,
    )
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_N_SHAPE_POINTS)
    parser.add_argument(
        "--max-candidates-per-axis",
        type=int,
        default=DEFAULT_MAX_CANDIDATES_PER_AXIS,
    )
    parser.add_argument("--strict-complete", dest="strict_complete", action="store_true", default=True)
    parser.add_argument("--no-strict-complete", dest="strict_complete", action="store_false")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if int(ns.k_max) < 1:
        raise ValueError("--k-max must be positive")
    if not (0.0 < float(ns.frequency_error_threshold) < 1.0):
        raise ValueError("--frequency-error-threshold must lie in (0, 1)")
    if int(ns.n_shape_points) < 51:
        raise ValueError("--n-shape-points must be at least 51")
    if int(ns.max_candidates_per_axis) < 3:
        raise ValueError("--max-candidates-per-axis must be at least 3")
    return Args(
        stage1_dir=repo_path(Path(ns.stage1_dir)),
        fixed_epsilon_dir=repo_path(Path(ns.fixed_epsilon_dir)),
        output_dir=repo_path(Path(ns.output_dir)),
        k_max=int(ns.k_max),
        frequency_error_threshold=float(ns.frequency_error_threshold),
        n_shape_points=int(ns.n_shape_points),
        max_candidates_per_axis=int(ns.max_candidates_per_axis),
        strict_complete=bool(ns.strict_complete),
    )


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if math.isnan(value_f):
            return "nan"
        if math.isinf(value_f):
            return "inf" if value_f > 0.0 else "-inf"
        return f"{value_f:.16e}"
    if isinstance(value, (bool, np.bool_)):
        return "true" if bool(value) else "false"
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required input CSV does not exist: {path}")
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: fmt(row.get(name, "")) for name in fields})
    return path


def finite_float(row: Mapping[str, object], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))
    except (TypeError, ValueError):
        return float("nan")


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def squared_lambda_delta(lambda_eb: float, lambda_timo: float) -> float:
    eb = float(lambda_eb)
    timo = float(lambda_timo)
    if not (isfinite(eb) and isfinite(timo)) or abs(timo) <= 1.0e-30:
        return float("nan")
    return abs(eb * eb - timo * timo) / (timo * timo)


def true_safe_prefix(
    deltas: Sequence[float] | Mapping[int, float],
    threshold: float = DEFAULT_FREQUENCY_ERROR_THRESHOLD,
    k_max: int = DEFAULT_K_MAX,
    *,
    strict: bool = True,
) -> int:
    if isinstance(deltas, Mapping):
        missing = [index for index in range(1, int(k_max) + 1) if index not in deltas]
        if strict and missing:
            raise ValueError(f"missing sorted mode indices: {missing}")
        values = [float(deltas[index]) for index in range(1, int(k_max) + 1) if index in deltas]
    else:
        values = [float(value) for value in deltas]
        if strict and len(values) < int(k_max):
            raise ValueError(f"expected {int(k_max)} delta values, found {len(values)}")
        values = values[: int(k_max)]
    if strict and len(values) != int(k_max):
        raise ValueError(f"expected {int(k_max)} delta values, found {len(values)}")
    count = 0
    for value in values:
        if not isfinite(value):
            if strict:
                raise ValueError("delta_f values must be finite")
            break
        if value <= float(threshold):
            count += 1
        else:
            break
    return count


def deltas_from_mode_rows(
    rows: Sequence[Mapping[str, object]],
    *,
    k_max: int = DEFAULT_K_MAX,
    strict: bool = True,
) -> list[float]:
    indexed: dict[int, float] = {}
    duplicates: list[int] = []
    for row in rows:
        try:
            index = int(float(row.get("sorted_index", row.get("eb_sorted_index", float("nan")))))
        except (TypeError, ValueError, OverflowError):
            continue
        if not 1 <= index <= int(k_max):
            continue
        if index in indexed:
            duplicates.append(index)
        indexed[index] = finite_float(row, "delta_f")
    if strict and duplicates:
        raise ValueError(f"duplicate sorted mode indices: {sorted(set(duplicates))}")
    missing = [index for index in range(1, int(k_max) + 1) if index not in indexed]
    if strict and missing:
        raise ValueError(f"missing sorted mode indices: {missing}")
    return [indexed[index] for index in range(1, int(k_max) + 1) if index in indexed]


def geometry_key(epsilon_0: float, beta_deg: float, mu: float, eta: float) -> tuple[float, float, float, float]:
    return tuple(
        round(float(value), GEOMETRY_ROUND_DIGITS)
        for value in (epsilon_0, beta_deg, mu, eta)
    )  # type: ignore[return-value]


def geometry_key_text(key: tuple[float, float, float, float]) -> str:
    return "|".join(f"{value:.12f}" for value in key)


def geometry_identifier(key: tuple[float, float, float, float]) -> str:
    digest = hashlib.sha256(geometry_key_text(key).encode("ascii")).hexdigest()[:16]
    return f"geometry_{digest}"


def source_index(row: Mapping[str, object]) -> int | None:
    for key in ("eb_sorted_index", "sorted_index"):
        try:
            value = int(float(row.get(key, float("nan"))))
        except (TypeError, ValueError, OverflowError):
            continue
        return value
    return None


def source_warning(row: Mapping[str, object]) -> tuple[bool, str]:
    def warning_present(value: object) -> bool:
        text_value = str(value).strip().lower()
        return text_value not in {"", "0", "false", "no", "n", "none", "nan", "[]", "()"}

    parts: list[str] = []
    if warning_present(row.get("root_warning", "false")):
        parts.append("root_warning")
    if warning_present(row.get("candidate_boundary", "false")):
        parts.append("candidate_boundary")
    tracking = str(row.get("tracking_warning", "")).strip()
    if warning_present(tracking):
        parts.append(tracking)
    explicit = str(row.get("explicit_root_warning", "")).strip()
    if warning_present(explicit):
        parts.append("explicit_root_warning")
    unique = list(dict.fromkeys(part for part in parts if part))
    return bool(unique), "; ".join(unique)


def normalize_source_row(source_dataset: str, row: Mapping[str, object]) -> dict[str, object]:
    index = source_index(row)
    warning, warning_text = source_warning(row)
    return {
        "source_dataset": source_dataset,
        "epsilon_0": finite_float(row, "epsilon_0"),
        "beta_deg": finite_float(row, "beta_deg"),
        "mu": finite_float(row, "mu"),
        "eta": finite_float(row, "eta"),
        "sorted_index": index if index is not None else "",
        "Lambda_EB": finite_float(row, "Lambda_EB"),
        "Lambda_Timo": finite_float(row, "Lambda_Timo"),
        "stored_delta_f": finite_float(row, "delta_f"),
        "delta_f": squared_lambda_delta(finite_float(row, "Lambda_EB"), finite_float(row, "Lambda_Timo")),
        "explicit_root_warning": warning,
        "source_warning": warning_text,
        "source_notes": str(row.get("notes", "")),
        "_raw": dict(row),
    }


def validate_geometry_mode_rows(
    rows: Sequence[Mapping[str, object]],
    *,
    k_max: int = DEFAULT_K_MAX,
) -> list[str]:
    reasons: list[str] = []
    indices = [source_index(row) for row in rows]
    target_indices = [index for index in indices if index is not None and 1 <= index <= int(k_max)]
    counts = Counter(target_indices)
    missing = [index for index in range(1, int(k_max) + 1) if counts[index] == 0]
    duplicates = [index for index in range(1, int(k_max) + 1) if counts[index] > 1]
    if missing:
        reasons.append("missing_sorted_indices=" + ",".join(str(index) for index in missing))
    if duplicates:
        reasons.append("duplicate_sorted_indices=" + ",".join(str(index) for index in duplicates))
    selected = sorted(
        (row for row in rows if source_index(row) is not None and 1 <= int(source_index(row) or 0) <= int(k_max)),
        key=lambda row: int(source_index(row) or 0),
    )
    if any(
        not all(isfinite(finite_float(row, key)) for key in ("Lambda_EB", "Lambda_Timo", "delta_f"))
        for row in selected
    ):
        reasons.append("non_finite_frequency_or_delta")
    inconsistent: list[int] = []
    for row in selected:
        recomputed = squared_lambda_delta(finite_float(row, "Lambda_EB"), finite_float(row, "Lambda_Timo"))
        stored = finite_float(row, "stored_delta_f") if "stored_delta_f" in row else finite_float(row, "delta_f")
        tolerance = DELTA_ABS_TOLERANCE + DELTA_REL_TOLERANCE * max(abs(recomputed), abs(stored), 1.0e-12)
        if not (isfinite(recomputed) and isfinite(stored)) or abs(recomputed - stored) > tolerance:
            index = source_index(row)
            if index is not None:
                inconsistent.append(index)
    if inconsistent:
        reasons.append("stored_delta_f_mismatch_at=" + ",".join(str(index) for index in sorted(set(inconsistent))))
    if len(selected) == int(k_max):
        eb_values = [finite_float(row, "Lambda_EB") for row in selected]
        timo_values = [finite_float(row, "Lambda_Timo") for row in selected]
        if any(right <= left for left, right in zip(eb_values, eb_values[1:])):
            reasons.append("Lambda_EB_not_strictly_increasing")
        if any(right <= left for left, right in zip(timo_values, timo_values[1:])):
            reasons.append("Lambda_Timo_not_strictly_increasing")
    if any(bool(row.get("explicit_root_warning", False)) for row in selected):
        reasons.append("explicit_root_warning")
    return reasons


def load_source_records(
    source_dataset: str,
    source_path: Path,
    *,
    k_max: int,
) -> list[GeometryRecord]:
    source_rows = [
        row
        for row in read_csv(source_path)
        if str(row.get("comparison_type", "")).strip() == "sorted_index"
    ]
    grouped: dict[tuple[float, float, float, float], list[dict[str, object]]] = defaultdict(list)
    invalid_parameter_rows: list[dict[str, object]] = []
    for raw in source_rows:
        row = normalize_source_row(source_dataset, raw)
        params = [finite_float(row, key) for key in ("epsilon_0", "beta_deg", "mu", "eta")]
        if not all(isfinite(value) for value in params):
            invalid_parameter_rows.append(row)
            continue
        grouped[geometry_key(*params)].append(row)

    records: list[GeometryRecord] = []
    for key in sorted(grouped):
        rows = grouped[key]
        selected_rows = [
            row
            for row in rows
            if isinstance(row.get("sorted_index"), int) and 1 <= int(row["sorted_index"]) <= int(k_max)
        ]
        selected_rows.sort(key=lambda row: int(row["sorted_index"]))
        reasons = validate_geometry_mode_rows(selected_rows, k_max=int(k_max))
        complete = not any(reason.startswith(("missing_sorted_indices", "duplicate_sorted_indices")) for reason in reasons)
        record = GeometryRecord(
            source_dataset=source_dataset,
            geometry_id=geometry_identifier(key),
            geometry_key=key,
            epsilon_0=key[0],
            beta_deg=key[1],
            mu=key[2],
            eta=key[3],
            modes=selected_rows,
            source_file=source_path,
            complete=complete,
            included=not reasons,
            explicit_warning=any(bool(row.get("explicit_root_warning", False)) for row in selected_rows),
            exclusion_reasons=list(reasons),
        )
        records.append(record)

    for index, row in enumerate(invalid_parameter_rows, start=1):
        key = (float("nan"), float(index), float("nan"), float("nan"))
        records.append(
            GeometryRecord(
                source_dataset=source_dataset,
                geometry_id=f"invalid_geometry_{source_dataset}_{index:04d}",
                geometry_key=key,
                epsilon_0=finite_float(row, "epsilon_0"),
                beta_deg=finite_float(row, "beta_deg"),
                mu=finite_float(row, "mu"),
                eta=finite_float(row, "eta"),
                modes=[row],
                source_file=source_path,
                complete=False,
                included=False,
                exclusion_reasons=["non_finite_geometry_parameters"],
            )
        )
    return records


def incomplete_records(records: Sequence[GeometryRecord]) -> list[GeometryRecord]:
    return [
        record
        for record in records
        if any(reason.startswith(("missing_sorted_indices", "duplicate_sorted_indices")) for reason in record.exclusion_reasons)
    ]


def k10_rerun_message(k_max: int, incomplete: Sequence[GeometryRecord]) -> str:
    by_source = Counter(record.source_dataset for record in incomplete)
    details = ", ".join(f"{source}={count}" for source, count in sorted(by_source.items()))
    stage_candidates = max(int(k_max) + 6, 16)
    fixed_candidates = max(int(k_max) + 10, 20)
    return "\n".join(
        [
            f"K={int(k_max)} complete sorted rows are required for every geometry; incomplete geometries: {details}.",
            "Regenerate the source CSV files with:",
            f"  python scripts/analysis/thickness_mismatch/audits/audit_eb_validity_vs_timoshenko_stage1.py --n-reported-modes {int(k_max)} --n-candidate-roots {stage_candidates}",
            f"  python scripts/analysis/thickness_mismatch/audits/audit_eb_validity_fixed_epsilon_geometry_scan.py --n-reported-modes {int(k_max)} --n-candidate-roots {fixed_candidates}",
            "No K=8 fallback was used." if int(k_max) == 10 else "No lower-K fallback was used.",
        ]
    )


def exclusion_rows(records: Sequence[GeometryRecord], *, k_max: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in records:
        if record.included:
            continue
        indices = sorted(
            int(mode["sorted_index"])
            for mode in record.modes
            if isinstance(mode.get("sorted_index"), int)
        )
        rows.append(
            {
                "source_dataset": record.source_dataset,
                "geometry_id": record.geometry_id,
                "geometry_key": geometry_key_text(record.geometry_key) if all(isfinite(value) for value in record.geometry_key) else "invalid",
                "epsilon_0": record.epsilon_0,
                "beta_deg": record.beta_deg,
                "mu": record.mu,
                "eta": record.eta,
                "observed_sorted_indices": ",".join(str(index) for index in indices),
                "target_mode_count": int(k_max),
                "complete_K10": record.complete,
                "explicit_warning": record.explicit_warning,
                "exclusion_reason": "; ".join(record.exclusion_reasons),
                "source_file": rel(record.source_file),
                "notes": "excluded explicitly; no silent row removal",
            }
        )
    return rows


def reconstruct_predictors(
    records: Sequence[GeometryRecord],
    *,
    n_shape_points: int,
    k_max: int,
    threshold: float,
    operations: OperationCounts,
) -> None:
    for record in records:
        if not record.included:
            continue
        local = stage1.local_thickness_parameters(record.epsilon_0, record.mu, record.eta)
        failed = False
        for mode in record.modes:
            try:
                result = fixed_scan.eb_mode_result(
                    epsilon=record.epsilon_0,
                    beta_deg=record.beta_deg,
                    mu=record.mu,
                    eta=record.eta,
                    sorted_index=int(mode["sorted_index"]),
                    Lambda=finite_float(mode, "Lambda_EB"),
                    n_points=int(n_shape_points),
                    root_warnings=(),
                )
                predictors = fixed_scan.eb_predictors(result)
            except (FloatingPointError, ValueError, OverflowError, np.linalg.LinAlgError) as exc:
                record.exclusion_reasons.append(
                    f"EB_predictor_reconstruction_failed_mode_{mode.get('sorted_index')}={exc}"
                )
                failed = True
                break
            operations.EB_mode_reconstructions += 1
            operations.EB_characteristic_matrix_evaluations += 1
            operations.SVD_6x6_calls += 1
            operations.EB_shape_grid_points += 2 * int(n_shape_points)
            operations.quadrature_calls += 12
            operations.quadrature_point_evaluations += 12 * int(n_shape_points)
            mode.update(predictors)
            mode.update(
                {
                    "epsilon_max": local["epsilon_max"],
                    "diameter_to_length_1": local["diameter_to_length_1"],
                    "diameter_to_length_2": local["diameter_to_length_2"],
                    "EB_axial_fraction": result.energy["axial_fraction"],
                    "EB_bending_fraction": result.energy["bending_fraction"],
                    "EB_classification": result.energy["classification"],
                    "true_pass_10": finite_float(mode, "delta_f") <= float(threshold),
                }
            )
        if failed:
            record.included = False
            continue
        required = (
            "epsilon_max",
            "chi_max_EB",
            "chi_eff_EB",
            "Theta_max_EB",
            "Pi_shear_EB",
            "Pi_rotary_EB",
            "Pi_EB",
            "EB_axial_fraction",
            "EB_bending_fraction",
        )
        if any(not isfinite(finite_float(mode, key)) for mode in record.modes for key in required):
            record.exclusion_reasons.append("non_finite_recomputed_EB_predictor")
            record.included = False


def add_spectral_gaps(
    records: Sequence[GeometryRecord],
    *,
    k_max: int,
    operations: OperationCounts,
) -> None:
    for record in records:
        if not record.included:
            continue
        modes = sorted(record.modes, key=lambda row: int(row["sorted_index"]))
        pair_gaps: list[float] = []
        for left, right in zip(modes, modes[1:]):
            lambda_left = finite_float(left, "Lambda_EB")
            lambda_right = finite_float(right, "Lambda_EB")
            gap = 2.0 * abs(lambda_right - lambda_left) / (lambda_right + lambda_left)
            pair_gaps.append(float(gap))
            operations.spectral_gap_operations += 1
        for position, mode in enumerate(modes):
            left_gap = pair_gaps[position - 1] if position > 0 else float("nan")
            right_gap = pair_gaps[position] if position < len(pair_gaps) else float("nan")
            available = [value for value in (left_gap, right_gap) if isfinite(value)]
            mode["gap_left_EB"] = left_gap
            mode["gap_right_EB"] = right_gap
            mode["gap_EB"] = min(available) if available else float("nan")
            if position == 0:
                mode["gap_sidedness"] = "right_only"
            elif position == len(modes) - 1:
                mode["gap_sidedness"] = "left_only"
            else:
                mode["gap_sidedness"] = "two_sided"
            mode["mode10_upper_guard_source"] = "unavailable" if int(mode["sorted_index"]) == int(k_max) else "not_applicable"


def finalize_geometry_targets(
    records: Sequence[GeometryRecord],
    *,
    k_max: int,
    threshold: float,
) -> None:
    for record in records:
        if not record.included:
            continue
        deltas = [finite_float(mode, "delta_f") for mode in sorted(record.modes, key=lambda row: int(row["sorted_index"]))]
        record.N_true = true_safe_prefix(deltas, threshold=float(threshold), k_max=int(k_max), strict=True)


def fixed_predictor_consistency_rows(records: Sequence[GeometryRecord]) -> list[dict[str, object]]:
    fixed_records = [record for record in records if record.source_dataset == SOURCE_FIXED and record.included]
    numeric_fields = (
        "epsilon_max",
        "chi_max_EB",
        "chi_eff_EB",
        "Theta_max_EB",
        "Pi_shear_EB",
        "Pi_rotary_EB",
        "Pi_EB",
        "EB_axial_fraction",
        "EB_bending_fraction",
    )
    rows: list[dict[str, object]] = []
    for field_name in numeric_fields:
        differences: list[float] = []
        relative: list[float] = []
        mismatch_count = 0
        for record in fixed_records:
            for mode in record.modes:
                raw = mode.get("_raw", {})
                if not isinstance(raw, Mapping):
                    continue
                stored = finite_float(raw, field_name)
                recomputed = finite_float(mode, field_name)
                if not (isfinite(stored) and isfinite(recomputed)):
                    continue
                diff = abs(recomputed - stored)
                rel_diff = diff / max(abs(stored), abs(recomputed), 1.0e-30)
                differences.append(diff)
                relative.append(rel_diff)
                if diff > PREDICTOR_CONSISTENCY_TOLERANCE and rel_diff > PREDICTOR_CONSISTENCY_TOLERANCE:
                    mismatch_count += 1
        rows.append(
            {
                "source_dataset": SOURCE_FIXED,
                "quantity": field_name,
                "compared_row_count": len(differences),
                "max_abs_difference": max(differences) if differences else float("nan"),
                "max_rel_difference": max(relative) if relative else float("nan"),
                "mismatch_count": mismatch_count,
                "tolerance": PREDICTOR_CONSISTENCY_TOLERANCE,
                "status": "consistent" if differences and mismatch_count == 0 else ("no_stored_values" if not differences else "mismatch"),
                "notes": "recomputed from Lambda_EB and EB mode data only",
            }
        )
    classification_compared = 0
    classification_mismatches = 0
    for record in fixed_records:
        for mode in record.modes:
            raw = mode.get("_raw", {})
            if not isinstance(raw, Mapping):
                continue
            stored = str(raw.get("EB_classification", "")).strip()
            if not stored:
                continue
            classification_compared += 1
            if stored != str(mode.get("EB_classification", "")):
                classification_mismatches += 1
    rows.append(
        {
            "source_dataset": SOURCE_FIXED,
            "quantity": "EB_classification",
            "compared_row_count": classification_compared,
            "max_abs_difference": "",
            "max_rel_difference": "",
            "mismatch_count": classification_mismatches,
            "tolerance": "exact_string_match",
            "status": "consistent" if classification_compared and classification_mismatches == 0 else ("no_stored_values" if not classification_compared else "mismatch"),
            "notes": "Timoshenko classification is not used",
        }
    )
    return rows


def mode_audit_rows(records: Sequence[GeometryRecord], *, threshold: float) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in records:
        reason = "; ".join(record.exclusion_reasons)
        key_text = geometry_key_text(record.geometry_key) if all(isfinite(value) for value in record.geometry_key) else "invalid"
        for mode in record.modes:
            row = {
                **{key: mode.get(key, "") for key in MODE_AUDIT_FIELDS},
                "source_dataset": record.source_dataset,
                "geometry_id": record.geometry_id,
                "geometry_key": key_text,
                "epsilon_0": record.epsilon_0,
                "beta_deg": record.beta_deg,
                "mu": record.mu,
                "eta": record.eta,
                "true_pass_10": finite_float(mode, "delta_f") <= float(threshold) if isfinite(finite_float(mode, "delta_f")) else "",
                "included_in_certification": record.included,
                "exclusion_reason": reason,
            }
            rows.append(row)
    return sorted(rows, key=lambda row: (str(row["source_dataset"]), str(row["geometry_id"]), int(row.get("sorted_index", 0) or 0)))


def maximum_over_prefix(record: GeometryRecord, field_name: str) -> float:
    if record.N_true is None or record.N_true <= 0:
        return float("nan")
    values = [
        finite_float(mode, field_name)
        for mode in record.modes
        if int(mode["sorted_index"]) <= int(record.N_true) and isfinite(finite_float(mode, field_name))
    ]
    return max(values) if values else float("nan")


def geometry_metric_rows(records: Sequence[GeometryRecord], *, k_max: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for record in records:
        local = (
            stage1.local_thickness_parameters(record.epsilon_0, record.mu, record.eta)
            if all(isfinite(value) for value in (record.epsilon_0, record.mu, record.eta))
            else {}
        )
        gaps = [finite_float(mode, "gap_EB") for mode in record.modes if isfinite(finite_float(mode, "gap_EB"))]
        prefix_gaps = [
            finite_float(mode, "gap_EB")
            for mode in record.modes
            if record.N_true is not None
            and int(mode.get("sorted_index", 0) or 0) <= record.N_true
            and isfinite(finite_float(mode, "gap_EB"))
        ]
        d1 = float(local.get("diameter_to_length_1", float("nan")))
        d2 = float(local.get("diameter_to_length_2", float("nan")))
        first_failed = record.N_true + 1 if record.N_true is not None and record.N_true < int(k_max) else ""
        rows.append(
            {
                "source_dataset": record.source_dataset,
                "geometry_id": record.geometry_id,
                "geometry_key": geometry_key_text(record.geometry_key) if all(isfinite(value) for value in record.geometry_key) else "invalid",
                "epsilon_0": record.epsilon_0,
                "beta_deg": record.beta_deg,
                "mu": record.mu,
                "eta": record.eta,
                "N_true": record.N_true if record.N_true is not None else "",
                "first_true_failed_mode": first_failed,
                "target_mode_count": int(k_max),
                "complete_K10": record.complete,
                "data_quality_status": "included" if record.included else "excluded",
                "epsilon_1": local.get("epsilon_1", ""),
                "epsilon_2": local.get("epsilon_2", ""),
                "epsilon_max": local.get("epsilon_max", ""),
                "diameter_to_length_1": d1,
                "diameter_to_length_2": d2,
                "max_diameter_to_length": max(d1, d2) if isfinite(d1) and isfinite(d2) else float("nan"),
                "thin_rod_criterion_pass": d1 <= 0.1 and d2 <= 0.1 if isfinite(d1) and isfinite(d2) else "",
                "max_Pi_EB_true_prefix": maximum_over_prefix(record, "Pi_EB"),
                "max_Pi_shear_EB_true_prefix": maximum_over_prefix(record, "Pi_shear_EB"),
                "max_Pi_rotary_EB_true_prefix": maximum_over_prefix(record, "Pi_rotary_EB"),
                "max_chi_max_EB_true_prefix": maximum_over_prefix(record, "chi_max_EB"),
                "max_chi_eff_EB_true_prefix": maximum_over_prefix(record, "chi_eff_EB"),
                "max_Theta_max_EB_true_prefix": maximum_over_prefix(record, "Theta_max_EB"),
                "minimum_gap_EB_true_prefix": min(prefix_gaps) if prefix_gaps else float("nan"),
                "minimum_gap_EB_firstK": min(gaps) if gaps else float("nan"),
                "mean_gap_EB_firstK": float(np.mean(gaps)) if gaps else float("nan"),
                "mode10_gap_sidedness": next((mode.get("gap_sidedness", "") for mode in record.modes if int(mode.get("sorted_index", 0) or 0) == int(k_max)), ""),
                "explicit_root_warning": record.explicit_warning,
                "source_warnings": "; ".join(str(mode.get("source_warning", "")) for mode in record.modes if str(mode.get("source_warning", ""))),
                "exclusion_reason": "; ".join(record.exclusion_reasons),
                "notes": "thin-rod ratios are diagnostic and are not exclusion filters",
            }
        )
    return rows


def close_cluster_members(
    modes: Sequence[Mapping[str, object]],
    gap_min: float,
    *,
    operations: OperationCounts | None = None,
) -> set[int]:
    blocked: set[int] = set()
    ordered = sorted(modes, key=lambda row: int(row["sorted_index"]))
    for left, right in zip(ordered, ordered[1:]):
        if operations is not None:
            operations.scalar_predictor_comparisons += 1
        gap = finite_float(left, "gap_right_EB")
        if isfinite(gap) and gap < float(gap_min):
            blocked.add(int(left["sorted_index"]))
            blocked.add(int(right["sorted_index"]))
    if ordered:
        if operations is not None:
            operations.scalar_predictor_comparisons += 1
        right_guard = finite_float(ordered[-1], "gap_right_EB")
        if isfinite(right_guard) and right_guard < float(gap_min):
            blocked.add(int(ordered[-1]["sorted_index"]))
    return blocked


def predict_safe_prefix(
    modes: Sequence[Mapping[str, object]],
    rule: str,
    thresholds: Mapping[str, float],
    *,
    operations: OperationCounts | None = None,
) -> dict[str, object]:
    if rule not in RULES:
        raise ValueError(f"unknown rule: {rule}")
    ordered = sorted(modes, key=lambda row: int(row["sorted_index"]))
    blocked = (
        close_cluster_members(ordered, float(thresholds.get("gap_min", 0.0)), operations=operations)
        if rule in {"Rule_A_gap", "Rule_C", "Rule_D"}
        else set()
    )
    count = 0
    first_trigger = ""
    trigger_reason = ""
    cluster_triggered = False
    mixed_triggered = False
    for mode in ordered:
        index = int(mode["sorted_index"])
        passes = True
        reason = ""
        if rule in {"Rule_A", "Rule_A_gap"}:
            if operations is not None:
                operations.scalar_predictor_comparisons += 1
            passes = finite_float(mode, "Pi_EB") <= float(thresholds["T_pi"])
            reason = "Pi_EB_limit"
            if passes and rule == "Rule_A_gap" and index in blocked:
                passes = False
                reason = "close_EB_cluster_guard"
                cluster_triggered = True
        else:
            if operations is not None:
                operations.scalar_predictor_comparisons += 2
            shear_pass = finite_float(mode, "Pi_shear_EB") <= float(thresholds["T_shear"])
            rotary_pass = finite_float(mode, "Pi_rotary_EB") <= float(thresholds["T_rotary"])
            passes = shear_pass and rotary_pass
            reason = "Pi_shear_limit" if not shear_pass else ("Pi_rotary_limit" if not rotary_pass else "")
            if passes and rule in {"Rule_C", "Rule_D"} and index in blocked:
                passes = False
                reason = "close_EB_cluster_guard"
                cluster_triggered = True
            if passes and rule == "Rule_D":
                if operations is not None:
                    operations.scalar_predictor_comparisons += 1
                if str(mode.get("EB_classification", "")) == "mixed":
                    passes = False
                    reason = "mixed_EB_mode_guard"
                    mixed_triggered = True
        if not passes:
            first_trigger = index
            trigger_reason = reason
            break
        count += 1
    return {
        "N_hat": count,
        "first_predicted_failed_mode": count + 1 if count < len(ordered) else "",
        "first_trigger": first_trigger,
        "trigger_reason": trigger_reason,
        "cluster_guard_triggered": cluster_triggered,
        "mixed_mode_guard_triggered": mixed_triggered,
    }


def predictor_boundary_values(records: Sequence[GeometryRecord], field_name: str) -> list[float]:
    values: list[float] = []
    for record in records:
        if record.N_true is None:
            continue
        for index in (record.N_true, record.N_true + 1):
            if not 1 <= index <= len(record.modes):
                continue
            value = finite_float(record.modes[index - 1], field_name)
            if isfinite(value):
                values.append(value)
    return values


def prefix_extrema_values(
    records: Sequence[GeometryRecord],
    field_name: str,
    *,
    use_minimum: bool = False,
) -> list[float]:
    values: set[float] = set()
    for record in records:
        running: list[float] = []
        for mode in sorted(record.modes, key=lambda row: int(row["sorted_index"])):
            value = finite_float(mode, field_name)
            if isfinite(value):
                running.append(value)
            if running:
                values.add(min(running) if use_minimum else max(running))
    return sorted(values)


def deterministic_candidate_axis(
    values: Sequence[float],
    *,
    max_candidates: int,
    include_negative_infinity: bool,
    priority_values: Sequence[float] = (),
) -> tuple[list[float], bool]:
    finite_values = sorted({float(value) for value in values if isfinite(float(value))})
    all_values = ([float("-inf")] if include_negative_infinity else []) + finite_values
    if len(all_values) <= int(max_candidates):
        return all_values, True
    selected: set[float] = {all_values[0], all_values[-1]}
    finite_set = set(finite_values)
    if finite_values:
        selected.add(finite_values[0])
        selected.add(finite_values[-1])
    for value in sorted({float(item) for item in priority_values if isfinite(float(item))}):
        value_f = float(value)
        if value_f in finite_set:
            selected.add(value_f)
        if len(selected) >= int(max_candidates):
            break
    remaining = int(max_candidates) - len(selected)
    if remaining > 0:
        indices = np.linspace(0, len(all_values) - 1, remaining + 2)[1:-1]
        for index in indices:
            selected.add(all_values[int(round(float(index)))])
    if len(selected) < int(max_candidates):
        for value in all_values:
            selected.add(value)
            if len(selected) >= int(max_candidates):
                break
    return sorted(selected)[: int(max_candidates)], False


def threshold_candidate_axes(
    records: Sequence[GeometryRecord],
    rule: str,
    *,
    max_candidates_per_axis: int,
) -> tuple[dict[str, list[float]], bool]:
    modes = [mode for record in records for mode in record.modes]
    if rule == "Rule_A":
        values = sorted({finite_float(mode, "Pi_EB") for mode in modes if isfinite(finite_float(mode, "Pi_EB"))})
        return {"T_pi": [float("-inf"), *values]}, True
    if rule == "Rule_A_gap":
        total, total_exact = deterministic_candidate_axis(
            [finite_float(mode, "Pi_EB") for mode in modes],
            max_candidates=int(max_candidates_per_axis),
            include_negative_infinity=True,
            priority_values=[
                *predictor_boundary_values(records, "Pi_EB"),
                *prefix_extrema_values(records, "Pi_EB"),
            ],
        )
        axes = {"T_pi": total}
        exact = total_exact
    else:
        shear, shear_exact = deterministic_candidate_axis(
            [finite_float(mode, "Pi_shear_EB") for mode in modes],
            max_candidates=int(max_candidates_per_axis),
            include_negative_infinity=True,
            priority_values=[
                *predictor_boundary_values(records, "Pi_shear_EB"),
                *prefix_extrema_values(records, "Pi_shear_EB"),
            ],
        )
        rotary, rotary_exact = deterministic_candidate_axis(
            [finite_float(mode, "Pi_rotary_EB") for mode in modes],
            max_candidates=int(max_candidates_per_axis),
            include_negative_infinity=True,
            priority_values=[
                *predictor_boundary_values(records, "Pi_rotary_EB"),
                *prefix_extrema_values(records, "Pi_rotary_EB"),
            ],
        )
        axes = {"T_shear": shear, "T_rotary": rotary}
        exact = shear_exact and rotary_exact
    if rule in {"Rule_A_gap", "Rule_C", "Rule_D"}:
        raw_gaps = [
            finite_float(mode, "gap_right_EB")
            for mode in modes
            if isfinite(finite_float(mode, "gap_right_EB"))
        ]
        gap_candidates = [0.0, *[float(np.nextafter(value, math.inf)) for value in raw_gaps]]
        priority_gaps = [
            float(np.nextafter(value, math.inf))
            for value in prefix_extrema_values(records, "gap_EB", use_minimum=True)
        ]
        gap_axis, gap_exact = deterministic_candidate_axis(
            gap_candidates,
            max_candidates=int(max_candidates_per_axis),
            include_negative_infinity=False,
            priority_values=priority_gaps,
        )
        axes["gap_min"] = gap_axis
        exact = exact and gap_exact
    return axes, exact


def conservative_threshold_key(rule: str, thresholds: Mapping[str, float]) -> tuple[float, ...]:
    if rule == "Rule_A":
        return (float(thresholds["T_pi"]),)
    if rule == "Rule_A_gap":
        return (float(thresholds["T_pi"]), -float(thresholds.get("gap_min", 0.0)))
    if rule == "Rule_B":
        return (float(thresholds["T_shear"]), float(thresholds["T_rotary"]))
    return (
        float(thresholds["T_shear"]),
        float(thresholds["T_rotary"]),
        -float(thresholds.get("gap_min", 0.0)),
    )


def threshold_is_more_conservative(
    rule: str,
    candidate: Mapping[str, float],
    reference: Mapping[str, float],
) -> bool:
    if rule == "Rule_A":
        return float(candidate["T_pi"]) < float(reference["T_pi"])
    upper_names = ("T_pi",) if rule == "Rule_A_gap" else ("T_shear", "T_rotary")
    weakly_tighter = all(
        float(candidate[name]) <= float(reference[name])
        for name in upper_names
    )
    strictly_tighter = any(
        float(candidate[name]) < float(reference[name])
        for name in upper_names
    )
    if rule in {"Rule_A_gap", "Rule_C", "Rule_D"}:
        weakly_tighter = weakly_tighter and float(candidate["gap_min"]) >= float(reference["gap_min"])
        strictly_tighter = strictly_tighter or float(candidate["gap_min"]) > float(reference["gap_min"])
    return weakly_tighter and strictly_tighter


def prediction_metrics(predictions: Sequence[Mapping[str, object]], *, k_max: int) -> dict[str, object]:
    count = len(predictions)
    if count == 0:
        return {name: float("nan") for name in VALIDATION_FIELDS if name not in {"split_kind", "rule", "fold_count", "fold_statuses", "notes"}}
    n_true = np.array([float(row["N_true"]) for row in predictions], dtype=float)
    n_hat = np.array([float(row["N_hat"]) for row in predictions], dtype=float)
    false_safe = n_hat > n_true
    false_safe_frequency = np.maximum(n_hat - n_true, 0.0)
    conservative = np.maximum(n_true - n_hat, 0.0)
    retained = np.minimum(n_hat, n_true)
    usable = float(np.sum(n_true))
    return {
        "validation_geometry_count": count,
        "validation_frequency_count": count * int(k_max),
        "false_safe_geometry_count": int(np.sum(false_safe)),
        "false_safe_frequency_count": int(np.sum(false_safe_frequency)),
        "false_safe_geometry_rate": float(np.mean(false_safe)),
        "false_safe_frequency_rate": safe_ratio(float(np.sum(false_safe_frequency)), count * int(k_max)),
        "mean_N_true": float(np.mean(n_true)),
        "mean_N_hat": float(np.mean(n_hat)),
        "exact_N_match_count": int(np.sum(n_true == n_hat)),
        "exact_N_match_rate": float(np.mean(n_true == n_hat)),
        "mean_conservative_loss": float(np.mean(conservative)),
        "median_conservative_loss": float(np.median(conservative)),
        "max_conservative_loss": int(np.max(conservative)),
        "usable_EB_frequency_count": int(usable),
        "retained_EB_frequency_count": int(np.sum(retained)),
        "usable_frequency_retention": safe_ratio(float(np.sum(retained)), usable),
        "accepted_geometry_count": int(np.sum(n_hat > 0)),
        "fully_EB_geometry_count": int(np.sum(n_hat == int(k_max))),
        "zero_EB_geometry_count": int(np.sum(n_hat == 0)),
        "cluster_guard_trigger_count": sum(bool_value(row.get("cluster_guard_triggered", False)) for row in predictions),
        "mixed_mode_guard_trigger_count": sum(bool_value(row.get("mixed_mode_guard_triggered", False)) for row in predictions),
        "incomplete_geometry_count": sum(int(float(row.get("incomplete_geometry", 0) or 0)) for row in predictions),
        "explicit_warning_geometry_count": sum(int(float(row.get("explicit_warning_geometry", 0) or 0)) for row in predictions),
        "out_of_calibration_domain_count": sum(bool_value(row.get("out_of_calibration_domain", False)) for row in predictions),
    }


def calibration_predictions(
    records: Sequence[GeometryRecord],
    rule: str,
    thresholds: Mapping[str, float],
    *,
    operations: OperationCounts | None = None,
) -> list[dict[str, object]]:
    predictions: list[dict[str, object]] = []
    for record in records:
        detail = predict_safe_prefix(record.modes, rule, thresholds, operations=operations)
        predictions.append(
            {
                **detail,
                "N_true": int(record.N_true or 0),
                "cluster_guard_triggered": detail["cluster_guard_triggered"],
                "mixed_mode_guard_triggered": detail["mixed_mode_guard_triggered"],
                "out_of_calibration_domain": False,
                "incomplete_geometry": 0,
                "explicit_warning_geometry": int(record.explicit_warning),
            }
        )
    return predictions


def calibrate_rule(
    records: Sequence[GeometryRecord],
    rule: str,
    *,
    k_max: int,
    max_candidates_per_axis: int,
    operations: OperationCounts | None = None,
) -> dict[str, object]:
    if not records:
        return {
            "rule": rule,
            "status": "no_training_geometries",
            "thresholds": {},
            "search_is_exact": False,
            "candidate_count_per_axis": {},
            "evaluated_combination_count": 0,
            "metrics": {},
        }
    axes, search_exact = threshold_candidate_axes(
        records,
        rule,
        max_candidates_per_axis=int(max_candidates_per_axis),
    )
    axis_names = list(axes)
    best_thresholds: dict[str, float] | None = None
    best_objective = -1
    evaluated = 0
    for values in itertools.product(*(axes[name] for name in axis_names)):
        thresholds = dict(zip(axis_names, (float(value) for value in values), strict=True))
        evaluated += 1
        if operations is not None:
            operations.threshold_combinations_evaluated += 1
        objective = 0
        valid = True
        for record in records:
            detail = predict_safe_prefix(record.modes, rule, thresholds, operations=operations)
            n_hat = int(detail["N_hat"])
            if n_hat > int(record.N_true or 0):
                valid = False
                break
            objective += n_hat
        if not valid:
            continue
        if objective > best_objective:
            best_objective = objective
            best_thresholds = thresholds
        elif objective == best_objective and best_thresholds is not None:
            if threshold_is_more_conservative(rule, thresholds, best_thresholds):
                best_thresholds = thresholds
            elif not threshold_is_more_conservative(rule, best_thresholds, thresholds) and (
                conservative_threshold_key(rule, thresholds)
                < conservative_threshold_key(rule, best_thresholds)
            ):
                best_thresholds = thresholds
    if best_thresholds is None:
        raise RuntimeError(f"no zero-false-safe calibration candidate found for {rule}")
    best_predictions = calibration_predictions(records, rule, best_thresholds, operations=operations)
    metrics = prediction_metrics(best_predictions, k_max=int(k_max))
    return {
        "rule": rule,
        "status": "ok",
        "thresholds": best_thresholds,
        "objective_retained_frequencies": best_objective,
        "search_is_exact": search_exact,
        "candidate_count_per_axis": {name: len(values) for name, values in axes.items()},
        "evaluated_combination_count": evaluated,
        "metrics": metrics,
        "tie_break_notes": "equal objectives prefer Pareto-tighter Pi limits and, for Rules C/D, a higher gap guard; incomparable ties use deterministic axis order",
    }


def exclude_transfer_overlaps(
    train: Sequence[GeometryRecord],
    validation: Sequence[GeometryRecord],
) -> tuple[list[GeometryRecord], int]:
    train_keys = {record.geometry_key for record in train}
    kept = [record for record in validation if record.geometry_key not in train_keys]
    return kept, len(validation) - len(kept)


def fold_value_token(value: float) -> str:
    return f"{float(value):.12g}".replace("-", "m").replace(".", "p")


def build_folds(records: Sequence[GeometryRecord]) -> list[Fold]:
    stage = tuple(record for record in records if record.included and record.source_dataset == SOURCE_STAGE1)
    fixed = tuple(record for record in records if record.included and record.source_dataset == SOURCE_FIXED)
    folds: list[Fold] = []
    fixed_validation, overlap = exclude_transfer_overlaps(stage, fixed)
    folds.append(Fold("stage1_to_fixed_epsilon", "transfer", stage, tuple(fixed_validation), overlap_excluded=overlap))
    stage_validation, overlap_reverse = exclude_transfer_overlaps(fixed, stage)
    folds.append(Fold("fixed_epsilon_to_stage1", "transfer", fixed, tuple(stage_validation), overlap_excluded=overlap_reverse))

    def leave_one(source: tuple[GeometryRecord, ...], split_kind: str, attr: str) -> None:
        values = sorted({float(getattr(record, attr)) for record in source})
        for value in values:
            validation = tuple(record for record in source if abs(float(getattr(record, attr)) - value) <= 1.0e-12)
            train = tuple(record for record in source if abs(float(getattr(record, attr)) - value) > 1.0e-12)
            folds.append(Fold(split_kind, f"holdout_{attr}_{fold_value_token(value)}", train, validation))

    leave_one(stage, "stage1_leave_one_mu_out", "mu")
    leave_one(stage, "stage1_leave_one_epsilon_out", "epsilon_0")
    leave_one(fixed, "fixed_epsilon_leave_one_beta_out", "beta_deg")
    leave_one(fixed, "fixed_epsilon_leave_one_mu_out", "mu")
    leave_one(fixed, "fixed_epsilon_leave_one_eta_out", "eta")

    combined_by_key: dict[tuple[float, float, float, float], GeometryRecord] = {
        record.geometry_key: record for record in stage
    }
    for record in fixed:
        combined_by_key[record.geometry_key] = record
    combined = tuple(combined_by_key[key] for key in sorted(combined_by_key))
    folds.append(Fold("combined_in_sample", "optimistic", combined, combined, optimistic_in_sample=True))
    return folds


def fold_status(fold: Fold) -> str:
    if not fold.train:
        return "no_training_geometries"
    if not fold.validation:
        return "no_validation_geometries"
    if len(fold.train) < 2:
        return "limited_training_size"
    safe_frequency_count = sum(int(record.N_true or 0) for record in fold.train)
    unsafe_frequency_count = sum(
        len(record.modes) - int(record.N_true or 0)
        for record in fold.train
    )
    if safe_frequency_count == 0 or unsafe_frequency_count == 0:
        return "limited_train_target_variation"
    return "ok"


def out_of_calibration_domain(record: GeometryRecord, train: Sequence[GeometryRecord]) -> bool:
    if not train:
        return True
    for attr in ("epsilon_0", "beta_deg", "mu", "eta"):
        train_values = [float(getattr(item, attr)) for item in train]
        value = float(getattr(record, attr))
        if value < min(train_values) - 1.0e-12 or value > max(train_values) + 1.0e-12:
            return True
    return False


def thresholds_json(thresholds: Mapping[str, float]) -> str:
    return json.dumps({key: fmt(float(value)) for key, value in thresholds.items()}, sort_keys=True, separators=(",", ":"))


def calibration_row(
    fold: Fold,
    status: str,
    result: Mapping[str, object],
) -> dict[str, object]:
    thresholds = result.get("thresholds", {})
    thresholds_map = thresholds if isinstance(thresholds, Mapping) else {}
    metrics = result.get("metrics", {})
    metrics_map = metrics if isinstance(metrics, Mapping) else {}
    return {
        "split_kind": fold.split_kind,
        "fold_id": fold.fold_id,
        "fold_status": status,
        "rule": result.get("rule", ""),
        "train_count": len(fold.train),
        "T_pi": thresholds_map.get("T_pi", ""),
        "T_shear": thresholds_map.get("T_shear", ""),
        "T_rotary": thresholds_map.get("T_rotary", ""),
        "gap_min": thresholds_map.get("gap_min", ""),
        "thresholds_json": thresholds_json(thresholds_map) if thresholds_map else "{}",
        "objective_retained_frequencies": result.get("objective_retained_frequencies", ""),
        "calibration_geometry_count": metrics_map.get("validation_geometry_count", len(fold.train)),
        "calibration_false_safe_geometry_count": metrics_map.get("false_safe_geometry_count", ""),
        "calibration_false_safe_frequency_count": metrics_map.get("false_safe_frequency_count", ""),
        "mean_N_true": metrics_map.get("mean_N_true", ""),
        "mean_N_hat": metrics_map.get("mean_N_hat", ""),
        "exact_N_match_rate": metrics_map.get("exact_N_match_rate", ""),
        "mean_conservative_loss": metrics_map.get("mean_conservative_loss", ""),
        "max_conservative_loss": metrics_map.get("max_conservative_loss", ""),
        "usable_frequency_retention": metrics_map.get("usable_frequency_retention", ""),
        "search_is_exact": result.get("search_is_exact", ""),
        "candidate_count_per_axis": result.get("candidate_count_per_axis", {}),
        "evaluated_combination_count": result.get("evaluated_combination_count", 0),
        "tie_break_notes": result.get("tie_break_notes", ""),
        "overlap_geometry_count_excluded": fold.overlap_excluded,
        "notes": "optimistic in-sample diagnostic" if fold.optimistic_in_sample else "thresholds use train geometries only",
    }


def prediction_row(
    fold: Fold,
    status: str,
    rule: str,
    thresholds: Mapping[str, float],
    record: GeometryRecord,
    detail: Mapping[str, object],
) -> dict[str, object]:
    n_true = int(record.N_true or 0)
    n_hat = int(detail["N_hat"])
    return {
        "split_kind": fold.split_kind,
        "fold_id": fold.fold_id,
        "fold_status": status,
        "rule": rule,
        "source_dataset": record.source_dataset,
        "geometry_id": record.geometry_id,
        "epsilon_0": record.epsilon_0,
        "beta_deg": record.beta_deg,
        "mu": record.mu,
        "eta": record.eta,
        "N_true": n_true,
        "N_hat": n_hat,
        "false_safe": n_hat > n_true,
        "false_safe_frequency_count": max(n_hat - n_true, 0),
        "conservative_loss": max(n_true - n_hat, 0),
        "first_predicted_failed_mode": detail.get("first_predicted_failed_mode", ""),
        "first_trigger": detail.get("first_trigger", ""),
        "trigger_reason": detail.get("trigger_reason", ""),
        "cluster_guard_triggered": detail.get("cluster_guard_triggered", False),
        "mixed_mode_guard_triggered": detail.get("mixed_mode_guard_triggered", False),
        "out_of_calibration_domain": out_of_calibration_domain(record, fold.train),
        "thresholds_json": thresholds_json(thresholds),
        "incomplete_geometry": 0,
        "explicit_warning_geometry": int(record.explicit_warning),
        "notes": "same-geometry mode rows were kept in one fold",
    }


def operation_row(
    *,
    scope: str,
    cost_component: str,
    measurement_status: str,
    operations: OperationCounts | None = None,
    split_kind: str = "",
    fold_id: str = "",
    rule: str = "",
    geometry_count: int = 0,
    frequency_count: int = 0,
    eb_accepted: int | str = "",
    timo_required: int | str = "",
    notes: str = "",
) -> dict[str, object]:
    row: dict[str, object] = {
        "scope": scope,
        "split_kind": split_kind,
        "fold_id": fold_id,
        "rule": rule,
        "cost_component": cost_component,
        "measurement_status": measurement_status,
        "geometry_count": geometry_count,
        "frequency_count": frequency_count,
        "EB_frequencies_accepted": eb_accepted,
        "Timoshenko_frequencies_required": timo_required,
        "Timoshenko_frequencies_avoided": eb_accepted,
        "fraction_firstK_avoided_in_Timoshenko": safe_ratio(float(eb_accepted), frequency_count) if isinstance(eb_accepted, int) else "",
        "notes": notes,
    }
    if operations is not None:
        row.update(operations.as_dict())
    return row


def run_calibration_and_validation(
    records: Sequence[GeometryRecord],
    *,
    k_max: int,
    max_candidates_per_axis: int,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    calibration_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    operation_rows: list[dict[str, object]] = []
    for fold in build_folds(records):
        status = fold_status(fold)
        for rule in RULES:
            rule_operations = OperationCounts()
            if status == "ok":
                result = calibrate_rule(
                    fold.train,
                    rule,
                    k_max=int(k_max),
                    max_candidates_per_axis=int(max_candidates_per_axis),
                    operations=rule_operations,
                )
            else:
                result = {
                    "rule": rule,
                    "status": status,
                    "thresholds": {},
                    "search_is_exact": False,
                    "candidate_count_per_axis": {},
                    "evaluated_combination_count": 0,
                    "metrics": {},
                    "tie_break_notes": "not calibrated because the fold is not evaluable",
                }
            calibration_rows.append(calibration_row(fold, status, result))
            thresholds = result.get("thresholds", {})
            thresholds_map = thresholds if isinstance(thresholds, Mapping) else {}
            fold_predictions: list[dict[str, object]] = []
            if thresholds_map:
                for record in fold.validation:
                    detail = predict_safe_prefix(record.modes, rule, thresholds_map, operations=rule_operations)
                    row = prediction_row(fold, status, rule, thresholds_map, record, detail)
                    prediction_rows.append(row)
                    fold_predictions.append(row)
            eb_accepted = sum(int(row["N_hat"]) for row in fold_predictions)
            frequency_count = len(fold_predictions) * int(k_max)
            operation_rows.append(
                operation_row(
                    scope="fold_rule",
                    split_kind=fold.split_kind,
                    fold_id=fold.fold_id,
                    rule=rule,
                    cost_component="projected_hybrid_cost",
                    measurement_status=(
                        "measured_postprocessor_and_structural_projection"
                        if status == "ok"
                        else f"not_evaluated_{status}"
                    ),
                    operations=rule_operations,
                    geometry_count=len(fold_predictions),
                    frequency_count=frequency_count,
                    eb_accepted=eb_accepted,
                    timo_required=frequency_count - eb_accepted,
                    notes="avoided frequency count is structural and is not a FLOP count",
                )
            )
    return calibration_rows, prediction_rows, operation_rows


def validation_summary_rows(
    prediction_rows: Sequence[dict[str, object]],
    *,
    k_max: int,
    calibration_rows: Sequence[dict[str, object]] = (),
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    fold_metadata: dict[tuple[str, str], set[tuple[str, str, str]]] = defaultdict(set)
    for row in calibration_rows:
        split_kind = str(row["split_kind"])
        rule = str(row["rule"])
        item = (split_kind, str(row["fold_id"]), str(row["fold_status"]))
        fold_metadata[(split_kind, rule)].add(item)
        if split_kind != "combined_in_sample":
            fold_metadata[("all_held_out", rule)].add(item)
    for row in prediction_rows:
        split_kind = str(row["split_kind"])
        rule = str(row["rule"])
        grouped[(split_kind, rule)].append(row)
        fold_metadata[(split_kind, rule)].add(
            (split_kind, str(row["fold_id"]), str(row["fold_status"]))
        )
        if split_kind != "combined_in_sample":
            grouped[("all_held_out", str(row["rule"]))].append(row)
            fold_metadata[("all_held_out", rule)].add(
                (split_kind, str(row["fold_id"]), str(row["fold_status"]))
            )
    rows: list[dict[str, object]] = []
    for split_kind, rule in sorted(set(grouped) | set(fold_metadata)):
        predictions = grouped.get((split_kind, rule), [])
        metadata = fold_metadata.get((split_kind, rule), set())
        metrics = prediction_metrics(predictions, k_max=int(k_max))
        statuses = sorted({status for _source_split, _fold_id, status in metadata})
        rows.append(
            {
                "split_kind": split_kind,
                "rule": rule,
                "fold_count": len({(source_split, fold_id) for source_split, fold_id, _status in metadata}),
                "fold_statuses": ";".join(statuses),
                **metrics,
                "notes": (
                    "optimistic in-sample diagnostic"
                    if split_kind == "combined_in_sample"
                    else (
                        "held-out complete geometries"
                        if predictions
                        else "not evaluated; see fold_statuses and calibration audit"
                    )
                ),
            }
        )
    return rows


def overlap_audit_rows(records: Sequence[GeometryRecord]) -> list[dict[str, object]]:
    by_source: dict[str, dict[tuple[float, float, float, float], GeometryRecord]] = {
        SOURCE_STAGE1: {},
        SOURCE_FIXED: {},
    }
    for record in records:
        by_source.setdefault(record.source_dataset, {})[record.geometry_key] = record
    shared = sorted(set(by_source[SOURCE_STAGE1]) & set(by_source[SOURCE_FIXED]))
    rows: list[dict[str, object]] = []
    for key in shared:
        left = by_source[SOURCE_STAGE1][key]
        right = by_source[SOURCE_FIXED][key]

        def max_mode_difference(field_name: str) -> float:
            left_by_index = {int(mode["sorted_index"]): mode for mode in left.modes if isinstance(mode.get("sorted_index"), int)}
            right_by_index = {int(mode["sorted_index"]): mode for mode in right.modes if isinstance(mode.get("sorted_index"), int)}
            values = [
                abs(finite_float(left_by_index[index], field_name) - finite_float(right_by_index[index], field_name))
                for index in sorted(set(left_by_index) & set(right_by_index))
                if isfinite(finite_float(left_by_index[index], field_name)) and isfinite(finite_float(right_by_index[index], field_name))
            ]
            return max(values) if values else float("nan")

        rows.append(
            {
                "geometry_id": geometry_identifier(key),
                "geometry_key": geometry_key_text(key),
                "epsilon_0": key[0],
                "beta_deg": key[1],
                "mu": key[2],
                "eta": key[3],
                "stage1_included": left.included,
                "fixed_epsilon_included": right.included,
                "max_abs_Lambda_EB_difference": max_mode_difference("Lambda_EB"),
                "max_abs_Lambda_Timo_difference": max_mode_difference("Lambda_Timo"),
                "max_abs_delta_f_difference": max_mode_difference("delta_f"),
                "N_true_stage1": left.N_true if left.N_true is not None else "",
                "N_true_fixed_epsilon": right.N_true if right.N_true is not None else "",
                "N_true_difference": (left.N_true - right.N_true) if left.N_true is not None and right.N_true is not None else "",
                "max_abs_Pi_EB_difference": max_mode_difference("Pi_EB"),
                "max_abs_chi_max_EB_difference": max_mode_difference("chi_max_EB"),
                "duplicate_use_policy": "excluded from cross-source held-out validation; fixed-epsilon row preferred in combined in-sample data",
                "notes": "overlap retained only for consistency audit",
            }
        )
    return rows


def best_current_rule(
    validation_rows: Sequence[dict[str, object]],
    operation_rows: Sequence[dict[str, object]] = (),
) -> tuple[str, str]:
    rows = [
        row
        for row in validation_rows
        if row.get("split_kind") == "all_held_out"
        and isfinite(finite_float(row, "validation_geometry_count"))
        and finite_float(row, "validation_geometry_count") > 0.0
    ]
    if not rows:
        return "none", "no held-out predictions are available"
    zero_false_safe = [row for row in rows if finite_float(row, "false_safe_geometry_count") == 0.0]
    candidates = zero_false_safe if zero_false_safe else rows
    complexity = {"Rule_A": 1, "Rule_A_gap": 2, "Rule_B": 2, "Rule_C": 3, "Rule_D": 4}
    operation_cost: Counter[str] = Counter()
    for row in operation_rows:
        if row.get("scope") != "fold_rule" or row.get("split_kind") == "combined_in_sample":
            continue
        if not str(row.get("measurement_status", "")).startswith("measured_"):
            continue
        operation_cost[str(row.get("rule", ""))] += int(
            finite_float(row, "scalar_predictor_comparisons")
        )
    best = sorted(
        candidates,
        key=lambda row: (
            0.0 if zero_false_safe else finite_float(row, "false_safe_geometry_count"),
            0.0 if zero_false_safe else finite_float(row, "false_safe_frequency_count"),
            -finite_float(row, "usable_frequency_retention"),
            finite_float(row, "mean_conservative_loss"),
            operation_cost.get(str(row.get("rule", "")), math.inf),
            complexity.get(str(row.get("rule", "")), 99),
        ),
    )[0]
    if zero_false_safe:
        reason = "selected among rules with zero observed held-out false-safe"
    else:
        reason = "no rule achieved zero observed held-out false-safe; ranking first minimizes observed false-safe"
    return str(best["rule"]), reason


def markdown_number(row: Mapping[str, object], key: str, digits: int = 4) -> str:
    value = finite_float(row, key)
    return f"{value:.{digits}g}" if isfinite(value) else "n/a"


def write_report(
    output_dir: Path,
    *,
    args: Args,
    records: Sequence[GeometryRecord],
    calibration_rows: Sequence[dict[str, object]],
    prediction_rows: Sequence[dict[str, object]],
    validation_rows: Sequence[dict[str, object]],
    overlap_rows: Sequence[dict[str, object]],
    consistency_rows: Sequence[dict[str, object]],
    operation_rows: Sequence[dict[str, object]],
) -> Path:
    report = output_dir / "eb_safe_prefix_certification_report.md"
    included = [record for record in records if record.included]
    excluded = [record for record in records if not record.included]
    best_rule, best_reason = best_current_rule(validation_rows, operation_rows)
    heldout = [
        row
        for row in validation_rows
        if row.get("split_kind") == "all_held_out"
        and isfinite(finite_float(row, "validation_geometry_count"))
        and finite_float(row, "validation_geometry_count") > 0.0
    ]
    cross_source = [
        row
        for row in validation_rows
        if row.get("split_kind") in {"stage1_to_fixed_epsilon", "fixed_epsilon_to_stage1"}
    ]
    consistency_mismatches = sum(int(float(row.get("mismatch_count", 0) or 0)) for row in consistency_rows)
    exclusion_counts = Counter(
        reason.split("=", 1)[0]
        for record in excluded
        for reason in record.exclusion_reasons
    )
    postprocessor_operations = next(
        (row for row in operation_rows if row.get("scope") == "postprocessor"),
        {},
    )
    heldout_operations = [
        row
        for row in operation_rows
        if row.get("scope") == "fold_rule"
        and row.get("split_kind") != "combined_in_sample"
        and str(row.get("measurement_status", "")).startswith("measured_")
    ]
    next_points = (
        "Generate the complete K=10 source slices first, then prioritize parameter-domain boundaries "
        "flagged as out of calibration and intersections between varying epsilon and varying beta/eta."
    )
    lines = [
        "# EB Safe-Spectrum-Prefix Certification",
        "",
        "## Scope",
        "",
        f"- target sorted frequencies: `K={args.k_max}`",
        f"- dimensional-frequency discrepancy threshold: `{args.frequency_error_threshold:g}`",
        "- target: first continuous sorted-frequency prefix",
        "- Timoshenko is the one-dimensional reference model, not exact 3D elasticity",
        "- the postprocessor performs no new EB or Timoshenko root calculations",
        "- no 3D FEM, Gmsh, or CalculiX calculation was run",
        "",
        "## Data",
        "",
        "| source | geometries | included | excluded |",
        "|---|---:|---:|---:|",
    ]
    for source in (SOURCE_STAGE1, SOURCE_FIXED):
        source_records = [record for record in records if record.source_dataset == source]
        lines.append(
            f"| {source} | {len(source_records)} | {sum(record.included for record in source_records)} | {sum(not record.included for record in source_records)} |"
        )
    lines.extend(
        [
            "",
            f"- complete included K={args.k_max} geometries: `{len(included)}`",
            f"- excluded geometries: `{len(excluded)}`",
            f"- exact cross-source geometry overlaps: `{len(overlap_rows)}`",
            f"- geometries carrying explicit source warnings: `{sum(record.explicit_warning for record in records)}` (excluded from certification)",
            f"- fixed-source stored/recomputed predictor mismatches: `{consistency_mismatches}`",
            "",
            "The classical thin-rod diameter ratio is retained as diagnostic metadata and is not an exclusion filter.",
            "Source warning fields are honored when present; a zero count does not imply that unavailable source diagnostics were reconstructed.",
            "",
            "Exclusion reasons:",
        ]
    )
    if exclusion_counts:
        for reason, count in sorted(exclusion_counts.items()):
            lines.append(f"- `{reason}`: `{count}` geometries")
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "## Candidate Rules",
            "",
            "- Rule A applies one total `Pi_EB` upper limit.",
            "- Rule A-gap applies the same total `Pi_EB` limit plus the adjacent sorted EB gap guard.",
            "- Rule B applies separate EB shear and rotary upper limits.",
            "- Rule C adds an adjacent sorted EB gap guard and blocks complete close clusters.",
            "- Rule D adds the EB-only `mixed` modal-character guard to Rule C.",
            "",
            "No Timoshenko frequency, Timoshenko energy fraction, Timoshenko classification, cross-model MAC, or observed `delta_f` is used as an inference feature.",
            "",
            "## Calibration",
            "",
            "| split | fold | status | rule | thresholds | objective | calibration false-safe | retention | exact search |",
            "|---|---|---|---|---|---:|---:|---:|---|",
        ]
    )
    selected_calibration = [
        row
        for row in calibration_rows
        if row.get("split_kind") in {"stage1_to_fixed_epsilon", "fixed_epsilon_to_stage1", "combined_in_sample"}
    ]
    for row in selected_calibration:
        lines.append(
            "| "
            f"{row.get('split_kind', '')} | {row.get('fold_id', '')} | {row.get('fold_status', '')} | "
            f"{row.get('rule', '')} | `{row.get('thresholds_json', '{}')}` | "
            f"{row.get('objective_retained_frequencies', '')} | "
            f"{row.get('calibration_false_safe_geometry_count', '')} | "
            f"{markdown_number(row, 'usable_frequency_retention')} | {row.get('search_is_exact', '')} |"
        )
    lines.extend(
        [
            "",
            "Thresholds are calibrated separately on each train fold and are never modified after held-out targets are evaluated.",
            "",
            "## Held-Out Validation",
            "",
            "| split | rule | folds | statuses | geometries | false-safe geometries | false-safe frequencies | mean N_true | mean N_hat | retention | mean loss |",
            "|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in validation_rows:
        if row.get("split_kind") == "combined_in_sample":
            continue
        lines.append(
            "| "
            f"{row.get('split_kind', '')} | {row.get('rule', '')} | {row.get('fold_count', '')} | "
            f"{row.get('fold_statuses', '')} | {markdown_number(row, 'validation_geometry_count')} | "
            f"{markdown_number(row, 'false_safe_geometry_count')} | "
            f"{markdown_number(row, 'false_safe_frequency_count')} | {markdown_number(row, 'mean_N_true')} | "
            f"{markdown_number(row, 'mean_N_hat')} | {markdown_number(row, 'usable_frequency_retention')} | "
            f"{markdown_number(row, 'mean_conservative_loss')} |"
        )
    lines.extend(
        [
            "",
            "`all_held_out` aggregates prediction events across deterministic split families; a geometry held out in different transfer checks may therefore appear more than once. Per-split rows remain the primary transfer diagnostics.",
            "",
            "## Cross-Source Transfer",
            "",
            "Exact geometry duplicates are removed from the held-out side of both transfer directions. They remain in the overlap audit only.",
            "",
        ]
    )
    for row in cross_source:
        lines.append(
            f"- {row.get('split_kind')} / {row.get('rule')}: false-safe geometries "
            f"`{markdown_number(row, 'false_safe_geometry_count')}`, retention `{markdown_number(row, 'usable_frequency_retention')}`."
        )
    lines.extend(
        [
            "",
            "## Best Current Rule",
            "",
            f"Current diagnostic selection: **{best_rule}**; {best_reason.rstrip('.')}.",
            "",
            "The ranking first requires zero observed held-out false-safe when available; otherwise it minimizes observed false-safe. It then prefers higher safe-frequency retention, lower conservative loss, fewer measured scalar indicator comparisons, and simpler interpretation.",
            "",
            "## Operation Counts",
            "",
            "Raw primitive counts are stored in `eb_safe_prefix_operation_counts.csv`. Source root-operation counts that are absent from the input workflows are marked unavailable and are not inferred from wall-clock time.",
            "",
            "| measured postprocessor primitive | count |",
            "|---|---:|",
        ]
    )
    for field_name in (
        "EB_characteristic_matrix_evaluations",
        "SVD_6x6_calls",
        "EB_mode_reconstructions",
        "EB_shape_grid_points",
        "quadrature_calls",
        "quadrature_point_evaluations",
        "spectral_gap_operations",
    ):
        lines.append(f"| `{field_name}` | {postprocessor_operations.get(field_name, 0)} |")
    lines.extend(
        [
            "",
            "| rule | threshold combinations | scalar comparisons | EB frequencies accepted | Timoshenko frequencies required | avoided fraction |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for rule in RULES:
        rule_rows = [row for row in heldout_operations if row.get("rule") == rule]
        combinations = sum(int(finite_float(row, "threshold_combinations_evaluated")) for row in rule_rows)
        comparisons = sum(int(finite_float(row, "scalar_predictor_comparisons")) for row in rule_rows)
        accepted = sum(int(finite_float(row, "EB_frequencies_accepted")) for row in rule_rows)
        required = sum(int(finite_float(row, "Timoshenko_frequencies_required")) for row in rule_rows)
        avoided_fraction = safe_ratio(accepted, accepted + required)
        lines.append(
            f"| {rule} | {combinations} | {comparisons} | {accepted} | {required} | "
            f"{f'{avoided_fraction:.4g}' if isfinite(avoided_fraction) else 'n/a'} |"
        )
    lines.extend(
        [
            "",
            "The rule-level comparison counts include deterministic threshold calibration and held-out application on evaluable folds. Structural avoided-frequency counts are not a complete FLOP count. End-to-end EB/Timoshenko root-cost ratios remain unavailable until source workflows expose primitive counters.",
            "",
            "## Limitations",
            "",
            "- Evidence covers two limited parameter slices: varying epsilon/mu at beta=45, eta=0, and varying beta/mu/eta at epsilon=0.02.",
            "- Zero observed false-safe on finite held-out data is not a mathematical guarantee over the continuous parameter domain.",
            "- Timoshenko and 3D FEM validation status remains diagnostic.",
            f"- Mode {args.k_max} uses a one-sided left EB gap because root {args.k_max + 1} is not required and was not invented.",
            "- No geometry-only pre-solution criterion is provided.",
            "- No local Timoshenko fallback is implemented.",
            "",
            "## Recommendation",
            "",
        ]
    )
    by_rule = {str(row.get("rule")): row for row in heldout}
    for rule in RULES:
        row = by_rule.get(rule)
        if row is None:
            continue
        lines.append(
            f"- {rule}: held-out false-safe geometries `{row.get('false_safe_geometry_count')}`, "
            f"retention `{markdown_number(row, 'usable_frequency_retention')}`, mean conservative loss "
            f"`{markdown_number(row, 'mean_conservative_loss')}`."
        )

    def comparison_sentence(left_rule: str, right_rule: str, description: str) -> str:
        left = by_rule.get(left_rule)
        right = by_rule.get(right_rule)
        if left is None or right is None:
            return f"- {description}: unavailable because the required held-out folds were not evaluable."
        return (
            f"- {description}: false-safe geometries `{left.get('false_safe_geometry_count')}` to "
            f"`{right.get('false_safe_geometry_count')}`, retention "
            f"`{markdown_number(left, 'usable_frequency_retention')}` to "
            f"`{markdown_number(right, 'usable_frequency_retention')}`, and mean conservative loss "
            f"`{markdown_number(left, 'mean_conservative_loss')}` to "
            f"`{markdown_number(right, 'mean_conservative_loss')}`."
        )

    rule_a = by_rule.get("Rule_A")
    if rule_a is None:
        pi_assessment = "- `Pi_EB` sufficiency cannot yet be assessed because no held-out Rule A result is available."
    elif finite_float(rule_a, "false_safe_geometry_count") == 0.0:
        pi_assessment = (
            "- Rule A (`Pi_EB`) achieved zero observed held-out false-safe on the evaluated folds, "
            f"with retention `{markdown_number(rule_a, 'usable_frequency_retention')}`; this is finite-data evidence, not a continuous-domain guarantee."
        )
    else:
        pi_assessment = (
            "- Rule A (`Pi_EB`) did not achieve zero observed held-out false-safe: "
            f"`{rule_a.get('false_safe_geometry_count')}` geometries were false-safe."
        )

    if any(finite_float(row, "false_safe_geometry_count") == 0.0 for row in heldout):
        new_criterion_assessment = (
            "- At least one current rule achieved zero observed held-out false-safe. A new physical criterion should therefore be deferred until retention and the remaining conservative losses are judged on the complete K=10 datasets."
        )
    elif heldout:
        new_criterion_assessment = (
            "- Every evaluated rule retained held-out false-safe cases. Inspect those geometries and their cluster/character triggers before proposing a new physically motivated dependence."
        )
    else:
        new_criterion_assessment = (
            "- No held-out rule comparison is available, so there is no evidence yet for introducing a new physical criterion."
        )
    lines.extend(
        [
            "",
            pi_assessment,
            comparison_sentence("Rule_A", "Rule_A_gap", "Adjacent EB gap ablation (Rule A to Rule A-gap)"),
            comparison_sentence("Rule_A", "Rule_B", "Separate shear/rotary limits (Rule A to Rule B)"),
            comparison_sentence("Rule_B", "Rule_C", "Adjacent EB gap guard (Rule B to Rule C)"),
            comparison_sentence("Rule_C", "Rule_D", "EB-only modal-character guard (Rule C to Rule D)"),
            f"- EB shape/Pi construction used `{postprocessor_operations.get('EB_mode_reconstructions', 0)}` reconstructions and `{postprocessor_operations.get('quadrature_point_evaluations', 0)}` quadrature-point evaluations; a root-cost ratio cannot be claimed from the current source CSV files.",
            new_criterion_assessment,
            "",
            f"Next computation recommendation: {next_points}",
        ]
    )
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    stage_path = args.stage1_dir / "eb_timo_mode_level_metrics.csv"
    fixed_path = args.fixed_epsilon_dir / "fixed_epsilon_mode_metrics.csv"
    records = [
        *load_source_records(SOURCE_STAGE1, stage_path, k_max=args.k_max),
        *load_source_records(SOURCE_FIXED, fixed_path, k_max=args.k_max),
    ]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    incomplete = incomplete_records(records)
    if args.strict_complete and incomplete:
        exclusion_path = write_csv(
            args.output_dir / "eb_safe_prefix_exclusion_audit.csv",
            exclusion_rows(records, k_max=args.k_max),
            EXCLUSION_FIELDS,
        )
        raise ValueError(k10_rerun_message(args.k_max, incomplete) + f"\nExclusion audit: {rel(exclusion_path)}")

    operations = OperationCounts()
    reconstruct_predictors(
        records,
        n_shape_points=args.n_shape_points,
        k_max=args.k_max,
        threshold=args.frequency_error_threshold,
        operations=operations,
    )
    add_spectral_gaps(records, k_max=args.k_max, operations=operations)
    finalize_geometry_targets(
        records,
        k_max=args.k_max,
        threshold=args.frequency_error_threshold,
    )
    consistency = fixed_predictor_consistency_rows(records)
    overlaps = overlap_audit_rows(records)
    calibration, predictions, fold_operations = run_calibration_and_validation(
        records,
        k_max=args.k_max,
        max_candidates_per_axis=args.max_candidates_per_axis,
    )
    validation = validation_summary_rows(
        predictions,
        k_max=args.k_max,
        calibration_rows=calibration,
    )
    included_count = sum(record.included for record in records)
    operation_rows = [
        operation_row(
            scope="source_inputs",
            cost_component="common_EB_solution_cost",
            measurement_status="root_primitive_counts_unavailable_from_source_CSV",
            geometry_count=included_count,
            frequency_count=included_count * args.k_max,
            notes="EB roots already exist in source CSV files; no root counts inferred from timing",
        ),
        operation_row(
            scope="postprocessor",
            cost_component="incremental_indicator_cost",
            measurement_status="measured_exact_postprocessor_counts",
            operations=operations,
            geometry_count=included_count,
            frequency_count=included_count * args.k_max,
            notes="12 trapezoidal quadratures per reconstructed EB mode: four EB energy and eight predictor integrals",
        ),
        operation_row(
            scope="source_inputs",
            cost_component="reference_full_Timoshenko_cost",
            measurement_status="root_primitive_counts_unavailable_from_source_CSV",
            geometry_count=included_count,
            frequency_count=included_count * args.k_max,
            notes="full Timoshenko frequencies are reference targets already stored in the source CSV files",
        ),
        *fold_operations,
    ]
    output_paths = [
        write_csv(args.output_dir / "eb_safe_prefix_mode_audit.csv", mode_audit_rows(records, threshold=args.frequency_error_threshold), MODE_AUDIT_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_geometry_metrics.csv", geometry_metric_rows(records, k_max=args.k_max), GEOMETRY_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_rule_calibration.csv", calibration, CALIBRATION_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_fold_predictions.csv", predictions, PREDICTION_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_validation_summary.csv", validation, VALIDATION_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_exclusion_audit.csv", exclusion_rows(records, k_max=args.k_max), EXCLUSION_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_source_overlap_audit.csv", overlaps, OVERLAP_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_predictor_consistency_audit.csv", consistency, CONSISTENCY_FIELDS),
        write_csv(args.output_dir / "eb_safe_prefix_operation_counts.csv", operation_rows, OPERATION_FIELDS),
    ]
    report = write_report(
        args.output_dir,
        args=args,
        records=records,
        calibration_rows=calibration,
        prediction_rows=predictions,
        validation_rows=validation,
        overlap_rows=overlaps,
        consistency_rows=consistency,
        operation_rows=operation_rows,
    )
    print("generated EB safe-prefix certification outputs:")
    for path in [report, *output_paths]:
        print(f"  {rel(path)}")
    print(f"included/excluded geometries: {included_count}/{len(records) - included_count}")
    heldout_prediction_count = sum(
        row.get("split_kind") != "combined_in_sample"
        for row in predictions
    )
    print(f"prediction rows (total/held-out): {len(predictions)}/{heldout_prediction_count}")
    print("root calculations performed by postprocessor: 0")
    return {
        "args": args,
        "records": records,
        "calibration_rows": calibration,
        "prediction_rows": predictions,
        "validation_rows": validation,
        "operation_rows": operation_rows,
        "output_paths": output_paths,
        "report": report,
    }


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from None
