from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import dataclass, field, replace
from decimal import Decimal, ROUND_FLOOR
import json
import math
from math import isfinite
from pathlib import Path
import sys
import time
from typing import Callable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.compare_single_rod_eb_timoshenko import (  # noqa: E402
    fixed_fixed_eb_roots,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_fixed_epsilon_geometry_scan as fixed_scan,
)
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as root_workflow,
)
from scripts.lib import straight_rod_factorized_spectrum as factorized  # noqa: E402


DEFAULT_OUTPUT_DIR = Path("results") / "eb_epsilon_baseline_thresholds"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_epsilon_baseline_thresholds"
DEFAULT_EPSILON_MIN = 0.005
DEFAULT_EPSILON_MAX = 0.060
DEFAULT_PRIMARY_EPSILON_MIN = 0.010
DEFAULT_PRIMARY_EPSILON_MAX = 0.050
DEFAULT_COARSE_STEP = 0.0005
DEFAULT_K_MAX = 10
DEFAULT_N_SPECTRUM_ROOTS = 12
DEFAULT_N_CANDIDATE_ROOTS = 20
DEFAULT_EPSILON_TOLERANCE = 1.0e-6
DEFAULT_DELTA_TOLERANCE = 1.0e-6
DEFAULT_THRESHOLD = 0.10

# These constants affect only event/family diagnostics. They do not change any
# determinant, root scan, root tolerance, or recovery setting.
EVENT_PROXIMITY = 0.01
CLOSE_GAP_THRESHOLD = fixed_scan.DEFAULT_CLUSTER_GAP_THRESHOLD
FAMILY_AMBIGUITY_GAP = 2.0e-5
REFERENCE_ABS_TOLERANCE = 2.0e-5
REFERENCE_REL_TOLERANCE = 2.0e-6
MU_EB_ABS_TOLERANCE = 1.0e-8
MU_EB_REL_TOLERANCE = 1.0e-8
MU_TIMO_ABS_TOLERANCE = 1.0e-6
MU_TIMO_REL_TOLERANCE = 1.0e-7
SINGULAR_MULTIPLICITY_THRESHOLD = root_workflow.SVD_RECOVERY_MULTIPLICITY
MANDATORY_EPSILONS = (0.006, 0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060)
MU_AUDIT_VALUES = (0.0, 0.3, 0.7, 0.9)
MU_AUDIT_EPSILONS = (0.005, 0.010, 0.015, 0.020, 0.030, 0.040, 0.050, 0.060)
VERIFICATION_CANDIDATE_ROOTS = 24
ROUND_DIGITS = 12
GENERAL_MATCH_ABS_TOLERANCE = 2.0e-5
GENERAL_MATCH_REL_TOLERANCE = 2.0e-6
LOCAL_FULL_SVD_RADIUS = 0.004
REGRESSION_EPSILON_OFFSET = 1.0e-6
REGRESSION_CASES = (
    ("R1", 0.0228017578125, 8.2932, 8.3000),
    ("R2", 0.0159306640625, 9.9287, 9.9299),
    ("R3", 0.015625, 14.1773, 14.1796),
)

OUTPUT_NAMES = (
    "baseline_epsilon_point_summary.csv",
    "baseline_mode_metrics.csv",
    "baseline_critical_prefix_thresholds.csv",
    "baseline_threshold_refinement_audit.csv",
    "baseline_transition_audit.csv",
    "baseline_spectral_family_audit.csv",
    "baseline_mu_invariance_audit.csv",
    "baseline_solver_quality_audit.csv",
    "baseline_operation_counts.csv",
    "baseline_factorized_spectrum_audit.csv",
    "baseline_axial_completeness_audit.csv",
    "baseline_general_vs_factorized_spectrum_audit.csv",
    "baseline_legacy_threshold_comparison.csv",
    "baseline_regression_close_root_pairs.csv",
    "eb_epsilon_baseline_thresholds_report.md",
    "baseline_N_true_vs_epsilon.png",
    "baseline_delta_f_by_mode.png",
    "baseline_Delta_prefixes.png",
    "baseline_certified_prefix_staircase.png",
    "baseline_sorted_family_map.png",
    "baseline_adjacent_gap_map.png",
)

POINT_FIELDS = [
    "epsilon",
    "in_primary_range",
    "N_true",
    "first_failed_mode",
    "late_pass_count",
    "late_pass_indices",
    "max_delta_first10",
    "minimum_gap_first10",
    "root11_available",
    "family_signature_EB",
    "family_signature_Timo",
    "reorder_detected",
    "cluster_detected",
    "possible_multiplicity",
    "root_quality_status",
    "cache_status",
    "point_type",
    "notes",
]

MODE_FIELDS = [
    "epsilon",
    "sorted_index",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_f",
    "Delta_n",
    "true_pass_10",
    "EB_family",
    "Timo_family",
    "family_ambiguity",
    "gap_left_EB",
    "gap_right_EB",
    "gap_EB",
    "gap_left_Timo",
    "gap_right_Timo",
    "gap_Timo",
    "axial_reference_Lambda",
    "EB_reference_Lambda",
    "reference_match_errors",
    "possible_multiplicity",
    "smallest_singular_value_EB",
    "second_smallest_singular_value_EB",
    "smallest_singular_value_Timo",
    "second_smallest_singular_value_Timo",
    "source_warnings",
    "point_type",
]

THRESHOLD_FIELDS = [
    "prefix_n",
    "threshold_status",
    "epsilon_certified_n",
    "epsilon_safe_lower",
    "epsilon_unsafe_upper",
    "epsilon_star_estimate",
    "epsilon_star_display_5dp",
    "epsilon_certified_floor_5dp",
    "epsilon_certified_floor_4dp",
    "bracket_width",
    "Delta_n_safe",
    "Delta_n_unsafe",
    "triggering_sorted_indices",
    "triggering_delta",
    "triggering_EB_family",
    "triggering_Timo_family",
    "reorder_near_threshold",
    "cluster_near_threshold",
    "simultaneous_transition_group",
    "first_loss_iteration_count",
    "reentry_count_after_first_loss",
    "inside_primary_range",
    "monotonicity_status",
    "continuity_status",
    "threshold_continuity_status",
    "numerical_verification_status",
    "verification_estimate",
    "verification_abs_difference",
    "verification_bracket_agreement",
    "verification_root_order_agreement",
    "epsilon_near_n",
    "epsilon_buffer_n",
    "notes",
]

FACTORIZED_FIELDS = [
    "epsilon",
    "mu",
    "model",
    "sorted_index",
    "sorted_position",
    "Lambda_factorized",
    "family",
    "family_index",
    "detection_source",
    "block_det",
    "block_sigma_min",
    "full_matrix_sigma_min",
    "full_matrix_second_sigma",
    "off_block_max_abs",
    "multiplicity_group",
    "multiplicity_group_id",
    "within_family_duplicate",
    "cross_family_cluster",
    "cluster_gap",
    "cross_family_gap",
    "raw_general_match",
    "raw_general_missing",
    "local_general_recovery_status",
    "cache_status",
    "algorithm_version",
    "root_quality_status",
    "quality_status",
    "notes",
]

AXIAL_COMPLETENESS_FIELDS = [
    "epsilon",
    "mu",
    "axial_index",
    "Lambda_exact",
    "Lambda_axial_exact",
    "axial_block_det",
    "axial_block_sigma_min",
    "axial_block_second_sigma",
    "axial_block_sigma_ratio",
    "full_matrix_sigma_min",
    "full_matrix_sigma_ratio",
    "expected_within_first12",
    "present_in_factorized_first12",
    "factorized_present",
    "raw_general_present",
    "full_matrix_verified",
    "factorized_sorted_index",
    "sorted_position_factorized",
    "nearest_raw_general_root",
    "difference",
    "verification_status",
    "status",
    "algorithm_version",
    "notes",
]

GENERAL_VS_FACTORIZED_FIELDS = [
    "scope",
    "epsilon",
    "mu",
    "model",
    "factorized_sorted_index",
    "Lambda_factorized",
    "factorized_family",
    "factorized_family_index",
    "present_in_raw_general",
    "nearest_general_root",
    "Lambda_general_6x6",
    "general_sorted_index",
    "abs_difference",
    "match_tolerance",
    "match_status",
    "missing_general_root",
    "raw_general_missing_root",
    "same_family_when_known",
    "general_sigma_min_at_factorized_root",
    "local_full_SVD_root",
    "local_full_SVD_sigma_min",
    "full_matrix_SVD_confirmation",
    "general_cache_status",
    "general_root_warnings",
    "notes",
]

LEGACY_COMPARISON_FIELDS = [
    "prefix_n",
    "legacy_status",
    "legacy_threshold_status",
    "corrected_status",
    "corrected_threshold_status",
    "legacy_epsilon_safe_lower",
    "legacy_epsilon_certified_n",
    "corrected_epsilon_certified_n",
    "corrected_epsilon_safe_lower",
    "legacy_epsilon_star_estimate",
    "corrected_epsilon_star_estimate",
    "absolute_estimate_change",
    "absolute_difference",
    "relative_difference",
    "legacy_trigger_mode",
    "corrected_trigger_mode",
    "missing_root_in_legacy_bracket",
    "missing_root_family",
    "correction_reason",
    "scientific_status",
    "comparison_status",
    "legacy_source",
    "algorithm_version",
    "notes",
]

REGRESSION_FIELDS = [
    "case_id",
    "epsilon_nominal",
    "epsilon_offset",
    "epsilon",
    "bending_expected_window_center",
    "axial_expected_window_center",
    "Lambda_bending_factorized",
    "Lambda_axial_exact",
    "pair_gap",
    "bending_sorted_index",
    "axial_sorted_index",
    "both_roots_present",
    "cross_family_cluster",
    "general_bending_present",
    "general_axial_present",
    "general_missing_root_count",
    "full_matrix_SVD_bending",
    "full_matrix_SVD_axial",
    "status",
    "algorithm_version",
    "notes",
]

REFINEMENT_FIELDS = [
    "refinement_kind",
    "prefix_n",
    "iteration",
    "epsilon",
    "Delta_n",
    "N_true",
    "side",
    "maximum_triggering_sorted_mode",
    "root_warnings",
    "family_signature",
    "adjacent_gap",
    "cache_status",
    "root_quality_status",
    "interval_left",
    "interval_right",
    "event_reasons",
    "notes",
]

TRANSITION_FIELDS = [
    "event_type",
    "prefix_n",
    "sorted_index",
    "epsilon_left",
    "epsilon_right",
    "epsilon_event_estimate",
    "state_left",
    "state_right",
    "family_signature_left",
    "family_signature_right",
    "related_to_reorder",
    "related_to_cluster",
    "is_first_loss",
    "notes",
]

FAMILY_FIELDS = [
    "epsilon",
    "sorted_index",
    "Lambda_EB_general",
    "Lambda_EB_analytic_union",
    "EB_reference_family",
    "EB_reference_family_index",
    "EB_reference_abs_error",
    "EB_reference_rel_error",
    "Lambda_axial_reference",
    "Lambda_Timo_factorized",
    "Lambda_Timo_general_6x6",
    "Timo_family",
    "Timo_factorized_family_exact",
    "Timo_axial_abs_error",
    "Timo_axial_rel_error",
    "family_ambiguity",
    "adjacent_gap_EB",
    "adjacent_gap_Timo",
    "possible_multiplicity",
    "point_type",
    "notes",
]

MU_FIELDS = [
    "epsilon",
    "mu",
    "sorted_index",
    "Lambda_EB",
    "Lambda_EB_mu0",
    "EB_abs_difference",
    "EB_rel_difference",
    "Lambda_Timo",
    "Lambda_Timo_mu0",
    "Timo_abs_difference",
    "Timo_rel_difference",
    "N_true",
    "N_true_mu0",
    "N_true_equal",
    "family_order_EB_equal",
    "family_order_Timo_equal",
    "multiplicity_flags_equal",
    "multiplicity_groups_Timo_equal",
    "root_quality_status",
    "status",
    "notes",
]

QUALITY_FIELDS = [
    "scope",
    "epsilon",
    "mu",
    "n_candidate_roots",
    "EB_requested_root_count",
    "EB_found_root_count",
    "Timo_requested_root_count",
    "Timo_found_root_count",
    "EB_root_warning",
    "Timo_root_warning",
    "candidate_boundary_warning",
    "cache_status",
    "general_cache_status",
    "algorithm_version",
    "EB_retry_attempted",
    "Timo_retry_attempted",
    "SVD_recovery_attempted",
    "SVD_recovery_changed_roots",
    "root11_available",
    "strictly_increasing_EB",
    "strictly_increasing_Timo",
    "multiplicity_preserving_order_EB",
    "multiplicity_preserving_order_Timo",
    "analytic_EB_union_max_abs_error",
    "analytic_EB_union_max_rel_error",
    "analytic_EB_union_mismatch_count",
    "possible_multiplicity_count",
    "stored_recomputed_delta_mismatch_count",
    "root_quality_status",
    "quality_reasons",
    "point_type",
    "notes",
]

OPERATION_FIELDS = ["scope", "metric", "value", "measurement_status", "notes"]


@dataclass(frozen=True)
class Args:
    epsilon_min: float
    epsilon_max: float
    primary_epsilon_min: float
    primary_epsilon_max: float
    coarse_step: float
    k_max: int
    n_spectrum_roots: int
    n_candidate_roots: int
    epsilon_tolerance: float
    delta_tolerance: float
    threshold: float
    output_dir: Path
    cache_dir: Path
    reuse_cache: bool
    force_recompute: bool
    plot_only: bool
    smoke: bool


@dataclass(frozen=True)
class ReferenceRoot:
    value: float
    family: str
    family_index: int
    ambiguous: bool = False


@dataclass
class Evaluation:
    epsilon: float
    mu: float
    eb_roots: tuple[float, ...]
    timo_roots: tuple[float, ...]
    factorized_timo: factorized.FactorizedSpectrum
    raw_eb_roots: tuple[float, ...]
    raw_timo_roots: tuple[float, ...]
    deltas: tuple[float, ...]
    prefix_deltas: tuple[float, ...]
    n_true: int
    eb_families: tuple[str, ...]
    timo_families: tuple[str, ...]
    family_ambiguity: tuple[bool, ...]
    eb_references: tuple[ReferenceRoot, ...]
    axial_references: tuple[float, ...]
    eb_reference_abs_errors: tuple[float, ...]
    eb_reference_rel_errors: tuple[float, ...]
    timo_axial_abs_errors: tuple[float, ...]
    timo_axial_rel_errors: tuple[float, ...]
    eb_gaps: tuple[float, ...]
    timo_gaps: tuple[float, ...]
    eb_singular_values: tuple[tuple[float, float], ...]
    timo_singular_values: tuple[tuple[float, float], ...]
    possible_multiplicity: tuple[bool, ...]
    quality_status: str
    quality_reasons: tuple[str, ...]
    cache_status: str
    general_cache_status: str
    eb_found: int
    timo_found: int
    eb_warnings: tuple[str, ...]
    timo_warnings: tuple[str, ...]
    raw_eb_warnings: tuple[str, ...]
    raw_timo_warnings: tuple[str, ...]
    eb_retry_attempted: bool
    timo_retry_attempted: bool
    candidate_boundary_warning: bool
    svd_recovery_attempted: bool
    svd_recovery_changed_roots: bool
    svd_recovery_call_count: int
    point_types: set[str] = field(default_factory=set)

    @property
    def family_signature_eb(self) -> str:
        return ",".join(self.eb_families)

    @property
    def family_signature_timo(self) -> str:
        return ",".join(self.timo_families)

    @property
    def point_type(self) -> str:
        return ";".join(sorted(self.point_types))

    @property
    def resolved(self) -> bool:
        return self.quality_status == "resolved"


@dataclass
class EvaluationManager:
    args: Args
    n_candidate_roots: int | None = None
    reuse_cache: bool | None = None
    force_recompute: bool | None = None
    scope: str = "primary"
    runtime: fixed_scan.RuntimeTracker = field(default_factory=fixed_scan.RuntimeTracker)
    evaluations: dict[tuple[float, float], Evaluation] = field(default_factory=dict)
    analytic_reference_evaluations: int = 0
    svd_recovery_calls: int = 0
    epsilon_requests: int = 0

    def evaluate(self, epsilon: float, point_type: str, *, mu: float = 0.0) -> Evaluation:
        key = (round(float(epsilon), ROUND_DIGITS), round(float(mu), ROUND_DIGITS))
        self.epsilon_requests += 1
        if key in self.evaluations:
            result = self.evaluations[key]
            result.point_types.add(str(point_type))
            return result
        result = _solve_evaluation(
            self.args,
            float(epsilon),
            float(mu),
            n_candidate_roots=int(self.n_candidate_roots or self.args.n_candidate_roots),
            reuse_cache=self.args.reuse_cache if self.reuse_cache is None else bool(self.reuse_cache),
            force_recompute=self.args.force_recompute if self.force_recompute is None else bool(self.force_recompute),
            runtime=self.runtime,
            scope=self.scope,
        )
        result.point_types.add(str(point_type))
        self.evaluations[key] = result
        self.analytic_reference_evaluations += 1
        self.svd_recovery_calls += int(result.svd_recovery_call_count)
        return result

    def existing(self, epsilon: float, *, mu: float = 0.0) -> Evaluation:
        key = (round(float(epsilon), ROUND_DIGITS), round(float(mu), ROUND_DIGITS))
        if key not in self.evaluations:
            raise KeyError(f"evaluation not available for epsilon={epsilon}, mu={mu}")
        return self.evaluations[key]


def repo_path(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


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
    if isinstance(value, (dict, list, tuple, set)):
        payload = sorted(value) if isinstance(value, set) else value
        return json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return value


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field_name: fmt(row.get(field_name, "")) for field_name in fields})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"required baseline output does not exist: {path}")
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def finite_float(row: Mapping[str, object], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))
    except (TypeError, ValueError):
        return float("nan")


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Audit first-loss K=10 EB applicability thresholds for the straight homogeneous baseline.",
    )
    parser.add_argument("--epsilon-min", type=float, default=DEFAULT_EPSILON_MIN)
    parser.add_argument("--epsilon-max", type=float, default=DEFAULT_EPSILON_MAX)
    parser.add_argument("--primary-epsilon-min", type=float, default=DEFAULT_PRIMARY_EPSILON_MIN)
    parser.add_argument("--primary-epsilon-max", type=float, default=DEFAULT_PRIMARY_EPSILON_MAX)
    parser.add_argument("--coarse-step", type=float, default=DEFAULT_COARSE_STEP)
    parser.add_argument("--k-max", type=int, default=DEFAULT_K_MAX)
    parser.add_argument("--n-spectrum-roots", type=int, default=DEFAULT_N_SPECTRUM_ROOTS)
    parser.add_argument("--n-candidate-roots", type=int, default=DEFAULT_N_CANDIDATE_ROOTS)
    parser.add_argument("--epsilon-tolerance", type=float, default=DEFAULT_EPSILON_TOLERANCE)
    parser.add_argument("--delta-tolerance", type=float, default=DEFAULT_DELTA_TOLERANCE)
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    output_dir = repo_path(Path(ns.output_dir))
    if bool(ns.smoke) and Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
        output_dir = repo_path(SMOKE_OUTPUT_DIR)
    cache_dir = repo_path(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / "cache"
    args = Args(
        epsilon_min=float(ns.epsilon_min),
        epsilon_max=float(ns.epsilon_max),
        primary_epsilon_min=float(ns.primary_epsilon_min),
        primary_epsilon_max=float(ns.primary_epsilon_max),
        coarse_step=float(ns.coarse_step),
        k_max=int(ns.k_max),
        n_spectrum_roots=int(ns.n_spectrum_roots),
        n_candidate_roots=int(ns.n_candidate_roots),
        epsilon_tolerance=float(ns.epsilon_tolerance),
        delta_tolerance=float(ns.delta_tolerance),
        threshold=float(ns.threshold),
        output_dir=output_dir,
        cache_dir=cache_dir,
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        plot_only=bool(ns.plot_only),
        smoke=bool(ns.smoke),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if args.epsilon_min <= 0.0:
        raise ValueError("--epsilon-min must be positive")
    if args.epsilon_max <= args.epsilon_min:
        raise ValueError("--epsilon-max must be greater than --epsilon-min")
    if not (
        args.epsilon_min <= args.primary_epsilon_min
        <= args.primary_epsilon_max <= args.epsilon_max
    ):
        raise ValueError("the primary epsilon range must lie inside the full scan range")
    if args.coarse_step <= args.epsilon_tolerance:
        raise ValueError("--coarse-step must be greater than --epsilon-tolerance")
    if args.k_max != DEFAULT_K_MAX:
        raise ValueError("this baseline audit contract requires --k-max 10")
    if args.n_spectrum_roots < args.k_max + 1:
        raise ValueError("--n-spectrum-roots must be at least --k-max + 1")
    if args.n_candidate_roots <= args.n_spectrum_roots:
        raise ValueError("--n-candidate-roots must be greater than --n-spectrum-roots")
    if args.epsilon_tolerance <= 0.0 or args.delta_tolerance <= 0.0:
        raise ValueError("epsilon and delta tolerances must be positive")
    if not 0.0 < args.threshold < 1.0:
        raise ValueError("--threshold must lie in (0, 1)")


def running_prefix_maxima(deltas: Sequence[float]) -> tuple[float, ...]:
    maxima: list[float] = []
    current = float("-inf")
    for value in deltas:
        value_f = float(value)
        current = max(current, value_f) if isfinite(value_f) else float("nan")
        maxima.append(current)
        if not isfinite(current):
            maxima.extend([float("nan")] * (len(deltas) - len(maxima)))
            break
    return tuple(maxima)


def true_safe_prefix(deltas: Sequence[float], threshold: float = DEFAULT_THRESHOLD) -> int:
    count = 0
    for value in deltas:
        value_f = float(value)
        if not isfinite(value_f) or value_f > float(threshold):
            break
        count += 1
    return count


def relative_frequency_delta(lambda_eb: float, lambda_timo: float) -> float:
    eb = float(lambda_eb)
    timo = float(lambda_timo)
    if not (isfinite(eb) and isfinite(timo)) or abs(timo) <= 1.0e-30:
        return float("nan")
    return abs(eb * eb - timo * timo) / (timo * timo)


def normalized_adjacent_gaps(roots: Sequence[float]) -> tuple[float, ...]:
    values = [float(value) for value in roots]
    gaps: list[float] = []
    for left, right in zip(values, values[1:]):
        if not (isfinite(left) and isfinite(right)) or left + right <= 0.0:
            gaps.append(float("nan"))
        else:
            gaps.append(2.0 * abs(right - left) / (right + left))
    return tuple(gaps)


def sided_gap(gaps: Sequence[float], index_zero_based: int) -> tuple[float, float, float]:
    left = float(gaps[index_zero_based - 1]) if index_zero_based > 0 else float("nan")
    right = float(gaps[index_zero_based]) if index_zero_based < len(gaps) else float("nan")
    finite = [value for value in (left, right) if isfinite(value)]
    return left, right, min(finite) if finite else float("nan")


def analytic_axial_roots(epsilon: float, count: int) -> np.ndarray:
    epsilon_f = float(epsilon)
    if epsilon_f <= 0.0:
        raise ValueError("epsilon must be positive")
    indices = np.arange(1, int(count) + 1, dtype=float)
    # Project normalization: theta = epsilon*Lambda^2 and a fixed-fixed total
    # length of 2 gives theta*2=m*pi, hence Lambda=sqrt(m*pi/(2*epsilon)).
    return np.sqrt(indices * np.pi / (2.0 * epsilon_f))


_BENDING_CACHE: dict[int, np.ndarray] = {}


def analytic_bending_roots(count: int) -> np.ndarray:
    count_i = int(count)
    if count_i not in _BENDING_CACHE:
        # fixed_fixed_eb_roots returns alpha for cosh(alpha)cos(alpha)=1.
        # The project total-length-two coordinate has alpha=2*Lambda.
        _BENDING_CACHE[count_i] = np.asarray(fixed_fixed_eb_roots(count_i), dtype=float) / 2.0
    return np.array(_BENDING_CACHE[count_i], copy=True)


def analytic_eb_union(epsilon: float, count: int) -> tuple[ReferenceRoot, ...]:
    count_i = int(count)
    # Each component family contributes an increasing positive sequence, so
    # the first ``count`` roots of their union cannot require an index larger
    # than ``count`` from either family. Avoid evaluating unused high-alpha
    # mpmath roots, which are deliberately expensive at 60-digit precision.
    axial_count = count_i
    bending_count = count_i
    components = [
        *(ReferenceRoot(float(value), "axial", index) for index, value in enumerate(analytic_axial_roots(epsilon, axial_count), start=1)),
        *(ReferenceRoot(float(value), "bending_EB", index) for index, value in enumerate(analytic_bending_roots(bending_count), start=1)),
    ]
    components.sort(key=lambda item: (item.value, item.family, item.family_index))
    selected = components[:count_i]
    out: list[ReferenceRoot] = []
    for index, item in enumerate(selected):
        ambiguous = False
        for neighbor_index in (index - 1, index + 1):
            if not 0 <= neighbor_index < len(selected):
                continue
            neighbor = selected[neighbor_index]
            gap = 2.0 * abs(item.value - neighbor.value) / (item.value + neighbor.value)
            if neighbor.family != item.family and gap <= FAMILY_AMBIGUITY_GAP:
                ambiguous = True
        out.append(replace(item, ambiguous=ambiguous))
    return tuple(out)


def abs_rel_difference(value: float, reference: float) -> tuple[float, float]:
    if not (isfinite(float(value)) and isfinite(float(reference))):
        return float("nan"), float("nan")
    absolute = abs(float(value) - float(reference))
    relative = absolute / abs(float(reference)) if abs(float(reference)) > 0.0 else float("nan")
    return absolute, relative


def reference_match_ok(value: float, reference: float) -> bool:
    absolute, _relative = abs_rel_difference(value, reference)
    return isfinite(absolute) and absolute <= REFERENCE_ABS_TOLERANCE + REFERENCE_REL_TOLERANCE * abs(reference)


def nearest_axial_match(value: float, axial_roots: Sequence[float]) -> tuple[float, float, float]:
    finite = [float(item) for item in axial_roots if isfinite(float(item))]
    if not finite or not isfinite(float(value)):
        return float("nan"), float("nan"), float("nan")
    reference = min(finite, key=lambda item: abs(float(value) - item))
    absolute, relative = abs_rel_difference(float(value), reference)
    return reference, absolute, relative


def classify_timo_families(
    roots: Sequence[float],
    axial_roots: Sequence[float],
    gaps: Sequence[float],
) -> tuple[tuple[str, ...], tuple[bool, ...], tuple[float, ...], tuple[float, ...]]:
    families: list[str] = []
    ambiguous_values: list[bool] = []
    abs_errors: list[float] = []
    rel_errors: list[float] = []
    for index, root in enumerate(roots):
        reference, absolute, relative = nearest_axial_match(float(root), axial_roots)
        axial = isfinite(reference) and reference_match_ok(float(root), reference)
        left, right, local_gap = sided_gap(gaps, index)
        del left, right
        ambiguity = axial and isfinite(local_gap) and local_gap <= FAMILY_AMBIGUITY_GAP
        families.append("axial_bending_cluster" if ambiguity else ("axial" if axial else "bending_Timoshenko"))
        ambiguous_values.append(ambiguity)
        abs_errors.append(absolute)
        rel_errors.append(relative)
    return tuple(families), tuple(ambiguous_values), tuple(abs_errors), tuple(rel_errors)


def fixed_args(
    args: Args,
    epsilon: float,
    *,
    n_candidate_roots: int,
    reuse_cache: bool,
    force_recompute: bool,
) -> fixed_scan.Args:
    return fixed_scan.Args(
        epsilon=float(epsilon),
        beta_values=(0.0,),
        mu_values=(0.0,),
        eta_values=(0.0,),
        n_reported_modes=args.n_spectrum_roots,
        n_candidate_roots=int(n_candidate_roots),
        n_shape_points=fixed_scan.DEFAULT_SHAPE_POINTS,
        cluster_gap_threshold=CLOSE_GAP_THRESHOLD,
        benchmark_local_timo=False,
        workers=1,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        reuse_cache=bool(reuse_cache),
        force_recompute=bool(force_recompute),
        plot_only=False,
        smoke=args.smoke,
    )


def roots_strictly_increasing(roots: Sequence[float], count: int) -> bool:
    values = [float(value) for value in roots[: int(count)]]
    return len(values) == int(count) and all(isfinite(value) for value in values) and all(
        right > left for left, right in zip(values, values[1:])
    )


def roots_multiplicity_preserving_order(roots: Sequence[float], count: int) -> bool:
    values = [float(value) for value in roots[: int(count)]]
    return len(values) == int(count) and all(isfinite(value) for value in values) and all(
        right >= left for left, right in zip(values, values[1:])
    )


def factorized_cache_dir(args: Args) -> Path:
    return args.output_dir / f"cache_{factorized.ALGORITHM_VERSION}"


def conservative_decimal_floor(value: object, digits: int) -> object:
    try:
        number = Decimal(str(value))
    except Exception:
        return ""
    if not number.is_finite():
        return ""
    quantum = Decimal(1).scaleb(-int(digits))
    return format(number.quantize(quantum, rounding=ROUND_FLOOR), f".{int(digits)}f")


def presentation_round(value: object, digits: int) -> object:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return ""
    return f"{number:.{int(digits)}f}" if isfinite(number) else ""


def cache_status(eb_entry: fixed_scan.RootCacheEntry, timo_entry: fixed_scan.RootCacheEntry) -> str:
    if eb_entry.cache_hit and timo_entry.cache_hit:
        return "hit"
    if not eb_entry.cache_hit and not timo_entry.cache_hit:
        return "miss"
    return "mixed"


def singular_values_for_roots(
    epsilon: float,
    mu: float,
    roots: Sequence[float],
    model: str,
    count: int,
) -> tuple[tuple[float, float], ...]:
    case = root_workflow.CaseSpec(mu=float(mu), eta=0.0, epsilon=float(epsilon))
    values: list[tuple[float, float]] = []
    for root in roots[: int(count)]:
        if not isfinite(float(root)):
            values.append((float("nan"), float("nan")))
        else:
            values.append(root_workflow.model_singular_values(case, 0.0, float(root), model))
    while len(values) < int(count):
        values.append((float("nan"), float("nan")))
    return tuple(values)


def eb_full_matrix_singular_values(
    Lambda: float,
    mu: float,
    epsilon: float,
) -> tuple[float, float]:
    try:
        matrix = root_workflow.assemble_clamped_coupled_matrix_eta(
            float(Lambda), 0.0, float(mu), float(epsilon), 0.0
        )
        reference = root_workflow.assemble_clamped_coupled_matrix_eta(
            float(Lambda) + factorized.FULL_MATRIX_SVD_ROW_SCALE_OFFSET,
            0.0,
            float(mu),
            float(epsilon),
            0.0,
        )
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan"), float("nan")
    if not (np.all(np.isfinite(matrix)) and np.all(np.isfinite(reference))):
        return float("nan"), float("nan")
    return factorized.fixed_row_scale_singular_values(matrix, reference)


def _entry_from_root_result(result: root_workflow.RootResult) -> fixed_scan.RootCacheEntry:
    return fixed_scan.RootCacheEntry(
        roots=tuple(float(value) for value in result.roots),
        warnings=tuple(result.warnings),
        root_count_found=int(result.root_count_found),
        lambda_max_used=float(result.lambda_max_used),
        scan_step_used=float(result.scan_step_used),
        retry_attempted=bool(result.retry_attempted),
        retry_changed_value=bool(result.retry_changed_value),
        notes=tuple(result.notes),
        cache_hit=False,
    )


def _analytic_reference_errors(
    roots: Sequence[float], references: Sequence[ReferenceRoot], count: int
) -> tuple[tuple[float, ...], tuple[float, ...], int]:
    absolute: list[float] = []
    relative: list[float] = []
    mismatches = 0
    for index in range(int(count)):
        value = float(roots[index]) if index < len(roots) else float("nan")
        reference = references[index].value if index < len(references) else float("nan")
        abs_diff, rel_diff = abs_rel_difference(value, reference)
        absolute.append(abs_diff)
        relative.append(rel_diff)
        if not reference_match_ok(value, reference):
            mismatches += 1
    return tuple(absolute), tuple(relative), mismatches


def _try_existing_svd_recovery(
    epsilon: float,
    mu: float,
    entry: fixed_scan.RootCacheEntry,
    model: str,
    n_candidate_roots: int,
) -> fixed_scan.RootCacheEntry:
    case = root_workflow.CaseSpec(mu=float(mu), eta=0.0, epsilon=float(epsilon))
    finite = [float(value) for value in entry.roots if isfinite(float(value))]
    upper_hint = max(finite) if finite else None
    recovered = root_workflow.solve_model(
        case,
        0.0,
        int(n_candidate_roots),
        model,
        retry=True,
        upper_hint=upper_hint,
        sv_recovery=True,
        sv_target_indices=tuple(range(min(12, int(n_candidate_roots)))),
    )
    return _entry_from_root_result(recovered)


def _try_reference_guided_eb_recovery(
    epsilon: float,
    mu: float,
    entry: fixed_scan.RootCacheEntry,
    references: Sequence[ReferenceRoot],
    n_candidate_roots: int,
) -> tuple[fixed_scan.RootCacheEntry, int, int]:
    """Recover only analytically detected omissions through existing 6x6 helpers.

    The analytic union supplies disjoint search windows, not production root
    values.  Every added value is returned either by the frozen EB sign-scan
    helper or by the existing row-normalized 6x6 SVD minimizer, and it must
    satisfy both the general-matrix singular-value acceptance criterion and
    the independent-reference tolerance.  Existing global roots are retained
    whenever they already match a reference.  This local orchestration is what
    keeps a near axial/bending crossing from shifting the rest of the sorted
    spectrum without changing any determinant, scan step, or solver tolerance.
    """

    finite = [float(value) for value in entry.roots if isfinite(float(value))]
    unused = set(range(len(finite)))
    selected: list[float] = []
    sign_scan_calls = 0
    svd_calls = 0
    case = root_workflow.CaseSpec(mu=float(mu), eta=0.0, epsilon=float(epsilon))
    reference_values = [float(item.value) for item in references]

    for index, reference in enumerate(reference_values):
        existing_matches = [
            candidate_index
            for candidate_index in unused
            if reference_match_ok(finite[candidate_index], reference)
        ]
        if existing_matches:
            chosen_index = min(existing_matches, key=lambda candidate_index: abs(finite[candidate_index] - reference))
            selected.append(finite[chosen_index])
            unused.remove(chosen_index)
            continue

        left_gap = reference - reference_values[index - 1] if index > 0 else float("inf")
        right_gap = reference_values[index + 1] - reference if index + 1 < len(reference_values) else float("inf")
        finite_gaps = [gap for gap in (left_gap, right_gap) if isfinite(gap) and gap > 0.0]
        local_radius = min(0.02, 0.24 * min(finite_gaps)) if finite_gaps else 0.02
        local_radius = max(local_radius, 1.0e-8)
        left = max(root_workflow.ROOT_SCAN_START, reference - local_radius)
        right = reference + local_radius

        sign_scan_calls += 1
        raw = root_workflow.find_roots_scan_bisect_eta(
            beta=0.0,
            mu=float(mu),
            epsilon=float(epsilon),
            eta=0.0,
            n_roots=1,
            Lmin=float(left),
            Lmax=float(right),
            scan_step=float(root_workflow.EB_RETRY_SCAN_STEP),
        )
        candidates = [float(value) for value in raw if isfinite(float(value))]

        svd_calls += 1
        svd_root, svd_minimum, _svd_second = root_workflow.refine_svd_minimum(
            case,
            0.0,
            root_workflow.MODEL_EB,
            float(left),
            float(right),
        )
        if isfinite(svd_root) and svd_minimum <= root_workflow.SVD_RECOVERY_ACCEPT:
            candidates.append(float(svd_root))

        accepted: list[float] = []
        for candidate in candidates:
            singular_minimum, _singular_second = root_workflow.model_singular_values(
                case, 0.0, candidate, root_workflow.MODEL_EB
            )
            if (
                singular_minimum <= root_workflow.SVD_RECOVERY_ACCEPT
                and reference_match_ok(candidate, reference)
            ):
                accepted.append(candidate)
        if not accepted:
            return entry, sign_scan_calls, svd_calls
        selected.append(min(accepted, key=lambda candidate: abs(candidate - reference)))

    remaining = [finite[index] for index in sorted(unused)]
    recovered_roots = sorted([*selected, *remaining])[: int(n_candidate_roots)]
    if len(recovered_roots) < int(n_candidate_roots):
        return entry, sign_scan_calls, svd_calls
    recovered = replace(
        entry,
        roots=tuple(recovered_roots),
        root_count_found=min(len(recovered_roots), int(n_candidate_roots)),
        retry_attempted=True,
        retry_changed_value=not root_workflow.roots_close(entry.roots, recovered_roots),
        notes=tuple(
            [
                *entry.notes,
                "analytic-guided local recovery used existing general-6x6 sign/SVD helpers",
            ]
        ),
        cache_hit=False,
    )
    return recovered, sign_scan_calls, svd_calls


def _solve_evaluation(
    args: Args,
    epsilon: float,
    mu: float,
    *,
    n_candidate_roots: int,
    reuse_cache: bool,
    force_recompute: bool,
    runtime: fixed_scan.RuntimeTracker,
    scope: str,
) -> Evaluation:
    point_args = fixed_args(
        args,
        epsilon,
        n_candidate_roots=n_candidate_roots,
        reuse_cache=reuse_cache,
        force_recompute=force_recompute,
    )
    cache = fixed_scan.RootCache(point_args, runtime)
    # The target requests the first 12 roots from each theory.  The larger
    # candidate margin is retained for EB because it guards sorted-root loss
    # against the independent analytic union.  Asking the current low-spectrum
    # Timoshenko wrapper for unused roots beyond 12 can cross its documented
    # high-frequency basis range at epsilon=0.06 and trigger repeated upper-
    # bound growth; root 12 already supplies the required boundary margin.
    eb_entry = cache.roots(
        model=fixed_scan.MODEL_EB,
        beta_deg=0.0,
        mu=float(mu),
        eta=0.0,
        n_roots=int(n_candidate_roots),
    )
    raw_general_eb_entry = eb_entry
    finite_eb = [float(value) for value in eb_entry.roots if isfinite(float(value))]
    timo_entry = cache.roots(
        model=fixed_scan.MODEL_TIMO,
        beta_deg=0.0,
        mu=float(mu),
        eta=0.0,
        n_roots=int(args.n_spectrum_roots),
        upper_hint=max(finite_eb) if finite_eb else None,
    )
    references = analytic_eb_union(epsilon, args.n_spectrum_roots)
    eb_abs, eb_rel, reference_mismatches = _analytic_reference_errors(
        eb_entry.roots, references, args.n_spectrum_roots
    )
    svd_attempted = False
    svd_changed = False
    svd_call_count = 0
    if reference_mismatches and float(mu) == 0.0:
        svd_attempted = True
        svd_call_count += 1
        recovered = _try_existing_svd_recovery(
            epsilon,
            mu,
            eb_entry,
            root_workflow.MODEL_EB,
            n_candidate_roots,
        )
        runtime.eb_root_calls += 1
        recovered_abs, recovered_rel, recovered_mismatches = _analytic_reference_errors(
            recovered.roots, references, args.n_spectrum_roots
        )
        current_error = max((value for value in eb_abs if isfinite(value)), default=float("inf"))
        recovered_error = max((value for value in recovered_abs if isfinite(value)), default=float("inf"))
        if (recovered_mismatches, recovered_error) < (reference_mismatches, current_error):
            svd_changed = not root_workflow.roots_close(eb_entry.roots, recovered.roots)
            eb_entry = recovered
            eb_abs, eb_rel, reference_mismatches = recovered_abs, recovered_rel, recovered_mismatches

    if reference_mismatches:
        svd_attempted = True
        locally_recovered, sign_calls, local_svd_calls = _try_reference_guided_eb_recovery(
            epsilon,
            mu,
            eb_entry,
            references,
            n_candidate_roots,
        )
        runtime.eb_root_calls += sign_calls
        svd_call_count += local_svd_calls
        local_abs, local_rel, local_mismatches = _analytic_reference_errors(
            locally_recovered.roots, references, args.n_spectrum_roots
        )
        current_error = max((value for value in eb_abs if isfinite(value)), default=float("inf"))
        local_error = max((value for value in local_abs if isfinite(value)), default=float("inf"))
        if (local_mismatches, local_error) < (reference_mismatches, current_error):
            svd_changed = svd_changed or not root_workflow.roots_close(eb_entry.roots, locally_recovered.roots)
            eb_entry = locally_recovered
            eb_abs, eb_rel, reference_mismatches = local_abs, local_rel, local_mismatches

    spectrum_count = args.n_spectrum_roots
    factorized_cache = factorized.FactorizedSpectrumCache(
        factorized_cache_dir(args),
        reuse=bool(reuse_cache),
        force_recompute=bool(force_recompute),
    )
    factorized_timo = factorized_cache.get(
        float(epsilon),
        float(mu),
        int(spectrum_count),
        int(n_candidate_roots),
    )
    raw_eb_roots = tuple(
        float(value) for value in raw_general_eb_entry.roots[:spectrum_count]
    )
    raw_timo_roots = tuple(float(value) for value in timo_entry.roots[:spectrum_count])
    # The straight baseline source of truth is the exact family union.  The EB
    # family values are analytic; the Timoshenko bending values come only from
    # the exact 4x4 block extracted from the unchanged 6x6 matrix.  The raw
    # general roots above remain available for the independent completeness
    # audit and are never substituted into the corrected target.
    eb_roots = tuple(float(item.value) for item in references)
    timo_roots = tuple(float(item.value) for item in factorized_timo.roots)
    eb_abs = tuple(0.0 for _item in references)
    eb_rel = tuple(0.0 for _item in references)
    eb_gaps = normalized_adjacent_gaps(eb_roots)
    timo_gaps = normalized_adjacent_gaps(timo_roots)
    axial = analytic_axial_roots(epsilon, max(spectrum_count + 12, 24))
    eb_families = tuple(
        "axial_bending_cluster" if item.ambiguous else item.family for item in references
    )
    timo_families = tuple(
        "axial_bending_cluster" if item.cross_family_cluster else item.family
        for item in factorized_timo.roots
    )
    timo_ambiguity = tuple(item.cross_family_cluster for item in factorized_timo.roots)
    timo_ax_abs_list: list[float] = []
    timo_ax_rel_list: list[float] = []
    for item in factorized_timo.roots:
        _reference, absolute, relative = nearest_axial_match(item.value, axial)
        timo_ax_abs_list.append(absolute)
        timo_ax_rel_list.append(relative)
    timo_ax_abs = tuple(timo_ax_abs_list)
    timo_ax_rel = tuple(timo_ax_rel_list)
    family_ambiguity = tuple(
        bool(references[index].ambiguous or timo_ambiguity[index])
        for index in range(spectrum_count)
    )
    eb_singular = tuple(
        eb_full_matrix_singular_values(root, mu, epsilon)
        for root in eb_roots
    )
    timo_singular = tuple(
        factorized.full_matrix_singular_values(root, mu, epsilon)
        for root in timo_roots
    )
    multiplicity: list[bool] = []
    for index in range(spectrum_count):
        _left_eb, _right_eb, gap_eb = sided_gap(eb_gaps, index)
        _left_t, _right_t, gap_t = sided_gap(timo_gaps, index)
        second_eb = eb_singular[index][1]
        second_timo = timo_singular[index][1]
        multiplicity.append(
            family_ambiguity[index]
            or factorized_timo.roots[index].cross_family_cluster
            or (isfinite(gap_eb) and gap_eb <= FAMILY_AMBIGUITY_GAP)
            or (isfinite(gap_t) and gap_t <= FAMILY_AMBIGUITY_GAP)
            or (isfinite(second_eb) and second_eb <= SINGULAR_MULTIPLICITY_THRESHOLD)
            or (isfinite(second_timo) and second_timo <= SINGULAR_MULTIPLICITY_THRESHOLD)
        )

    deltas = tuple(
        relative_frequency_delta(eb_roots[index], timo_roots[index])
        for index in range(args.k_max)
    )
    prefix_deltas = running_prefix_maxima(deltas)
    n_true = true_safe_prefix(deltas, args.threshold)
    candidate_boundary = len(eb_roots) <= args.k_max or len(timo_roots) <= args.k_max
    reasons: list[str] = []
    if len(factorized_timo.roots) < spectrum_count:
        reasons.append(f"missing_factorized_Timo_roots={len(factorized_timo.roots)}/{spectrum_count}")
    if len(eb_roots) < spectrum_count or not all(isfinite(value) for value in eb_roots):
        reasons.append("nonfinite_EB_first12")
    if len(timo_roots) < spectrum_count or not all(isfinite(value) for value in timo_roots):
        reasons.append("nonfinite_Timo_first12")
    if not roots_multiplicity_preserving_order(eb_roots, spectrum_count):
        reasons.append("EB_roots_not_multiplicity_preserving_order")
    if not roots_multiplicity_preserving_order(timo_roots, spectrum_count):
        reasons.append("Timo_roots_not_multiplicity_preserving_order")
    if factorized_timo.bending_search.warnings:
        reasons.append("unresolved_factorized_bending_root_warning")
    if candidate_boundary:
        reasons.append("mode10_or_root11_at_candidate_boundary")
    off_block_values: list[float] = []
    for root in timo_roots:
        matrix, _bending, _axial, _warnings = factorized.straight_blocks(root, mu, epsilon)
        off_block_values.append(factorized.off_block_max_abs(matrix))
    if any(value > factorized.BLOCK_OFFDIAGONAL_ABS_TOL for value in off_block_values):
        reasons.append("nonzero_straight_off_block_entry")
    if any(
        not isfinite(values[0]) or values[0] > factorized.FULL_MATRIX_SVD_ACCEPT
        for values in timo_singular
    ):
        reasons.append("factorized_root_failed_full_matrix_SVD_confirmation")
    for item in factorized_timo.roots:
        if item.family == "axial" and (
            not isfinite(item.block_sigma_min)
            or item.block_sigma_min > factorized.AXIAL_BLOCK_SVD_ACCEPT
        ):
            reasons.append("exact_axial_root_failed_axial_block_confirmation")
            break
    if not all(isfinite(value) for value in deltas):
        reasons.append("nonfinite_target_delta")
    quality_status = "resolved" if not reasons else "unresolved"
    del scope
    return Evaluation(
        epsilon=float(epsilon),
        mu=float(mu),
        eb_roots=eb_roots,
        timo_roots=timo_roots,
        factorized_timo=factorized_timo,
        raw_eb_roots=raw_eb_roots,
        raw_timo_roots=raw_timo_roots,
        deltas=deltas,
        prefix_deltas=prefix_deltas,
        n_true=n_true,
        eb_families=eb_families,
        timo_families=timo_families,
        family_ambiguity=family_ambiguity,
        eb_references=references,
        axial_references=tuple(float(value) for value in axial),
        eb_reference_abs_errors=eb_abs,
        eb_reference_rel_errors=eb_rel,
        timo_axial_abs_errors=timo_ax_abs,
        timo_axial_rel_errors=timo_ax_rel,
        eb_gaps=eb_gaps,
        timo_gaps=timo_gaps,
        eb_singular_values=eb_singular,
        timo_singular_values=timo_singular,
        possible_multiplicity=tuple(multiplicity),
        quality_status=quality_status,
        quality_reasons=tuple(dict.fromkeys(reasons)),
        cache_status=factorized_timo.cache_status,
        general_cache_status=cache_status(raw_general_eb_entry, timo_entry),
        eb_found=len(eb_roots),
        timo_found=len(timo_roots),
        eb_warnings=(),
        timo_warnings=tuple(factorized_timo.bending_search.warnings),
        raw_eb_warnings=tuple(raw_general_eb_entry.warnings),
        raw_timo_warnings=tuple(timo_entry.warnings),
        eb_retry_attempted=bool(eb_entry.retry_attempted),
        timo_retry_attempted=bool(timo_entry.retry_attempted),
        candidate_boundary_warning=candidate_boundary,
        svd_recovery_attempted=svd_attempted,
        svd_recovery_changed_roots=svd_changed,
        svd_recovery_call_count=svd_call_count,
    )


def epsilon_grid(args: Args) -> list[float]:
    if args.smoke:
        return [float(value) for value in np.linspace(args.epsilon_min, args.epsilon_max, 4)]
    count = int(math.floor((args.epsilon_max - args.epsilon_min) / args.coarse_step + 1.0e-10))
    values = [args.epsilon_min + index * args.coarse_step for index in range(count + 1)]
    if not math.isclose(values[-1], args.epsilon_max, rel_tol=0.0, abs_tol=1.0e-12):
        values.append(args.epsilon_max)
    values.extend(
        value for value in MANDATORY_EPSILONS
        if args.epsilon_min - 1.0e-12 <= value <= args.epsilon_max + 1.0e-12
    )
    return sorted({round(float(value), ROUND_DIGITS) for value in values})


def interval_event_reasons(
    left: Evaluation,
    midpoint: Evaluation,
    right: Evaluation,
    *,
    threshold: float,
    event_proximity: float = EVENT_PROXIMITY,
    close_gap_threshold: float = CLOSE_GAP_THRESHOLD,
) -> tuple[str, ...]:
    reasons: list[str] = []
    if left.n_true != right.n_true:
        reasons.append("N_true_transition")
    for prefix_index in range(min(len(left.prefix_deltas), len(right.prefix_deltas))):
        a = left.prefix_deltas[prefix_index] - threshold
        b = right.prefix_deltas[prefix_index] - threshold
        m = midpoint.prefix_deltas[prefix_index]
        if isfinite(a) and isfinite(b) and a * b <= 0.0:
            reasons.append(f"prefix_{prefix_index + 1}_threshold_endpoint")
        if (
            isfinite(m)
            and isfinite(left.prefix_deltas[prefix_index])
            and isfinite(right.prefix_deltas[prefix_index])
            and m > max(left.prefix_deltas[prefix_index], right.prefix_deltas[prefix_index])
            and abs(m - threshold) <= event_proximity
        ):
            reasons.append(f"prefix_{prefix_index + 1}_midpoint_peak")
        if any(
            isfinite(value) and abs(value - threshold) <= event_proximity
            for value in (left.prefix_deltas[prefix_index], m, right.prefix_deltas[prefix_index])
        ):
            reasons.append(f"prefix_{prefix_index + 1}_near_threshold")
    for mode_index in range(min(len(left.deltas), len(right.deltas))):
        a = left.deltas[mode_index] - threshold
        b = right.deltas[mode_index] - threshold
        if isfinite(a) and isfinite(b) and a * b <= 0.0:
            reasons.append(f"mode_{mode_index + 1}_threshold_endpoint")
    if (
        left.family_signature_eb != right.family_signature_eb
        or left.family_signature_timo != right.family_signature_timo
        or midpoint.family_signature_eb not in {left.family_signature_eb, right.family_signature_eb}
        or midpoint.family_signature_timo not in {left.family_signature_timo, right.family_signature_timo}
    ):
        reasons.append("family_reorder")
    gaps = [*left.eb_gaps, *left.timo_gaps, *midpoint.eb_gaps, *midpoint.timo_gaps, *right.eb_gaps, *right.timo_gaps]
    if any(isfinite(value) and value <= close_gap_threshold for value in gaps):
        reasons.append("close_cluster")
    if any((*left.possible_multiplicity, *midpoint.possible_multiplicity, *right.possible_multiplicity)):
        reasons.append("possible_multiplicity")
    if left.svd_recovery_attempted or midpoint.svd_recovery_attempted or right.svd_recovery_attempted:
        reasons.append("SVD_recovery")
    if not (left.resolved and midpoint.resolved and right.resolved):
        reasons.append("root_quality_warning")
    return tuple(dict.fromkeys(reasons))


def refinement_row(
    evaluation: Evaluation,
    *,
    kind: str,
    prefix_n: int | str,
    iteration: int,
    side: str,
    interval_left: float,
    interval_right: float,
    event_reasons: Sequence[str] = (),
    notes: str = "",
) -> dict[str, object]:
    prefix_index = int(prefix_n) - 1 if str(prefix_n).isdigit() and int(prefix_n) > 0 else None
    delta = evaluation.prefix_deltas[prefix_index] if prefix_index is not None else float("nan")
    trigger = ""
    if prefix_index is not None and prefix_index < len(evaluation.deltas):
        prefix_values = evaluation.deltas[: prefix_index + 1]
        maximum = max(prefix_values)
        trigger = ";".join(str(index + 1) for index, value in enumerate(prefix_values) if math.isclose(value, maximum, rel_tol=1.0e-10, abs_tol=1.0e-12))
    gaps = [value for value in (*evaluation.eb_gaps, *evaluation.timo_gaps) if isfinite(value)]
    warnings = [*evaluation.eb_warnings, *evaluation.timo_warnings, *evaluation.quality_reasons]
    return {
        "refinement_kind": kind,
        "prefix_n": prefix_n,
        "iteration": iteration,
        "epsilon": evaluation.epsilon,
        "Delta_n": delta,
        "N_true": evaluation.n_true,
        "side": side,
        "maximum_triggering_sorted_mode": trigger,
        "root_warnings": ";".join(warnings),
        "family_signature": f"EB={evaluation.family_signature_eb}|Timo={evaluation.family_signature_timo}",
        "adjacent_gap": min(gaps) if gaps else float("nan"),
        "cache_status": evaluation.cache_status,
        "root_quality_status": evaluation.quality_status,
        "interval_left": interval_left,
        "interval_right": interval_right,
        "event_reasons": ";".join(event_reasons),
        "notes": notes,
    }


def adaptive_screen(
    args: Args,
    manager: EvaluationManager,
    coarse_values: Sequence[float],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    audit_rows: list[dict[str, object]] = []
    interval_rows: list[dict[str, object]] = []
    candidates: list[tuple[float, float, tuple[str, ...]]] = []
    for left_value, right_value in zip(coarse_values, coarse_values[1:]):
        left = manager.evaluate(left_value, "coarse")
        right = manager.evaluate(right_value, "coarse")
        midpoint_value = 0.5 * (left_value + right_value)
        midpoint = manager.evaluate(midpoint_value, "midpoint")
        reasons = interval_event_reasons(left, midpoint, right, threshold=args.threshold)
        if reasons:
            candidates.append((left_value, right_value, reasons))
            interval_rows.append(
                {
                    "event_type": "adaptive_event_interval",
                    "prefix_n": "",
                    "sorted_index": "",
                    "epsilon_left": left_value,
                    "epsilon_right": right_value,
                    "epsilon_event_estimate": midpoint_value,
                    "state_left": left.n_true,
                    "state_right": right.n_true,
                    "family_signature_left": f"{left.family_signature_eb}|{left.family_signature_timo}",
                    "family_signature_right": f"{right.family_signature_eb}|{right.family_signature_timo}",
                    "related_to_reorder": "family_reorder" in reasons,
                    "related_to_cluster": "close_cluster" in reasons,
                    "is_first_loss": False,
                    "notes": ";".join(reasons),
                }
            )
    if args.smoke and candidates:
        candidates = [min(candidates, key=lambda item: abs(0.5 * (item[0] + item[1]) - 0.5 * (args.epsilon_min + args.epsilon_max)))]
    max_depth = 1 if args.smoke else 2
    queue = [(left, right, reasons, 0) for left, right, reasons in candidates]
    while queue:
        left_value, right_value, reasons, depth = queue.pop(0)
        if depth >= max_depth or right_value - left_value <= max(8.0 * args.epsilon_tolerance, args.coarse_step / 8.0):
            continue
        midpoint = 0.5 * (left_value + right_value)
        quarter_values = (0.5 * (left_value + midpoint), 0.5 * (midpoint + right_value))
        for quarter in quarter_values:
            evaluation = manager.evaluate(quarter, "event_refinement")
            audit_rows.append(
                refinement_row(
                    evaluation,
                    kind="event_adaptive",
                    prefix_n="",
                    iteration=depth + 1,
                    side="screen",
                    interval_left=left_value,
                    interval_right=right_value,
                    event_reasons=reasons,
                    notes="deterministic quarter-point event refinement",
                )
            )
        if not args.smoke:
            # All five values already exist.  Looking them up without adding a
            # new label preserves the actual coarse/midpoint/refinement counts
            # when recursively screening sub-intervals.
            sub_points = [
                manager.existing(left_value),
                manager.existing(quarter_values[0]),
                manager.existing(midpoint),
                manager.existing(quarter_values[1]),
                manager.existing(right_value),
            ]
            for first, middle, last in ((sub_points[0], sub_points[1], sub_points[2]), (sub_points[2], sub_points[3], sub_points[4])):
                sub_reasons = interval_event_reasons(first, middle, last, threshold=args.threshold)
                if sub_reasons:
                    queue.append((first.epsilon, last.epsilon, sub_reasons, depth + 1))
    return audit_rows, interval_rows


def ordered_baseline_evaluations(manager: EvaluationManager, args: Args) -> list[Evaluation]:
    return sorted(
        (
            evaluation for (_key_eps, key_mu), evaluation in manager.evaluations.items()
            if abs(key_mu) <= 1.0e-12
            and args.epsilon_min - 1.0e-12 <= evaluation.epsilon <= args.epsilon_max + 1.0e-12
            and not evaluation.point_types.intersection({"verification_safe", "verification_unsafe", "verification_refinement"})
        ),
        key=lambda item: item.epsilon,
    )


def first_loss_bracket(
    evaluations: Sequence[Evaluation], prefix_n: int, threshold: float
) -> tuple[str, Evaluation | None, Evaluation | None, str]:
    ordered = sorted(evaluations, key=lambda item: item.epsilon)
    if not ordered:
        return "not_reached_in_scan_range", None, None, "no evaluated epsilon points"
    index = int(prefix_n) - 1
    first = ordered[0]
    if not first.resolved:
        return "unresolved_root_quality", None, None, "minimum scan point is unresolved"
    if first.prefix_deltas[index] > threshold:
        return "below_scan_range", None, first, "first verified scan point is already unsafe"
    previous = first
    for current in ordered[1:]:
        if not current.resolved:
            if previous.prefix_deltas[index] <= threshold:
                return "unresolved_root_quality", previous, current, "unresolved point precedes a possible first loss"
            previous = current
            continue
        if current.prefix_deltas[index] > threshold:
            if not previous.resolved:
                return "unresolved_root_quality", previous, current, "unsafe point follows unresolved root quality"
            if previous.prefix_deltas[index] <= threshold:
                return "resolved", previous, current, "first safe-to-unsafe bracket"
        previous = current
    return "not_reached_in_scan_range", ordered[-1], None, "no unsafe point was reached"


def refine_first_loss(
    args: Args,
    manager: EvaluationManager,
    prefix_n: int,
    safe: Evaluation,
    unsafe: Evaluation,
    *,
    kind: str = "threshold_bisection",
) -> tuple[Evaluation, Evaluation, list[dict[str, object]], int]:
    safe_point = safe
    unsafe_point = unsafe
    rows = [
        refinement_row(
            safe_point,
            kind=kind,
            prefix_n=prefix_n,
            iteration=0,
            side="safe",
            interval_left=safe.epsilon,
            interval_right=unsafe.epsilon,
            notes="initial safe bracket endpoint",
        ),
        refinement_row(
            unsafe_point,
            kind=kind,
            prefix_n=prefix_n,
            iteration=0,
            side="unsafe",
            interval_left=safe.epsilon,
            interval_right=unsafe.epsilon,
            notes="initial unsafe bracket endpoint",
        ),
    ]
    iteration = 0
    max_iterations = 80
    while unsafe_point.epsilon - safe_point.epsilon > args.epsilon_tolerance and iteration < max_iterations:
        iteration += 1
        midpoint = 0.5 * (safe_point.epsilon + unsafe_point.epsilon)
        evaluation = manager.evaluate(
            midpoint,
            "verification_refinement" if kind.startswith("verification") else "threshold_refinement",
        )
        if not evaluation.resolved:
            rows.append(
                refinement_row(
                    evaluation,
                    kind=kind,
                    prefix_n=prefix_n,
                    iteration=iteration,
                    side="unresolved",
                    interval_left=safe_point.epsilon,
                    interval_right=unsafe_point.epsilon,
                    notes="root quality blocks further threshold refinement",
                )
            )
            break
        delta = evaluation.prefix_deltas[prefix_n - 1]
        side = "safe" if delta <= args.threshold else "unsafe"
        rows.append(
            refinement_row(
                evaluation,
                kind=kind,
                prefix_n=prefix_n,
                iteration=iteration,
                side=side,
                interval_left=safe_point.epsilon,
                interval_right=unsafe_point.epsilon,
                notes=(
                    "Delta_n is within delta_tolerance"
                    if abs(delta - args.threshold) <= args.delta_tolerance
                    else "epsilon bracket refinement"
                ),
            )
        )
        if delta <= args.threshold:
            safe_point = evaluation
        else:
            unsafe_point = evaluation
    return safe_point, unsafe_point, rows, iteration


def triggering_modes(evaluation: Evaluation, prefix_n: int) -> list[int]:
    values = evaluation.deltas[: int(prefix_n)]
    maximum = max(values)
    return [
        index + 1 for index, value in enumerate(values)
        if math.isclose(value, maximum, rel_tol=1.0e-9, abs_tol=1.0e-11)
    ]


def transition_rows(evaluations: Sequence[Evaluation], args: Args) -> list[dict[str, object]]:
    ordered = sorted(evaluations, key=lambda item: item.epsilon)
    rows: list[dict[str, object]] = []
    first_loss_seen = [False] * args.k_max
    for left, right in zip(ordered, ordered[1:]):
        related_reorder = (
            left.family_signature_eb != right.family_signature_eb
            or left.family_signature_timo != right.family_signature_timo
        )
        related_cluster = any(
            isfinite(value) and value <= CLOSE_GAP_THRESHOLD
            for value in (*left.eb_gaps, *left.timo_gaps, *right.eb_gaps, *right.timo_gaps)
        )
        common = {
            "epsilon_left": left.epsilon,
            "epsilon_right": right.epsilon,
            "epsilon_event_estimate": 0.5 * (left.epsilon + right.epsilon),
            "family_signature_left": f"{left.family_signature_eb}|{left.family_signature_timo}",
            "family_signature_right": f"{right.family_signature_eb}|{right.family_signature_timo}",
            "related_to_reorder": related_reorder,
            "related_to_cluster": related_cluster,
        }
        if left.n_true != right.n_true:
            rows.append({"event_type": "N_true_change", "prefix_n": "", "sorted_index": "", "state_left": left.n_true, "state_right": right.n_true, "is_first_loss": False, "notes": "sampled/refined N_true transition", **common})
        if related_reorder:
            rows.append({"event_type": "family_reorder", "prefix_n": "", "sorted_index": "", "state_left": left.family_signature_eb, "state_right": right.family_signature_eb, "is_first_loss": False, "notes": "sorted analytic family signature changed", **common})
        for prefix_index in range(args.k_max):
            left_safe = left.prefix_deltas[prefix_index] <= args.threshold
            right_safe = right.prefix_deltas[prefix_index] <= args.threshold
            if left_safe == right_safe:
                continue
            first_loss = left_safe and not right_safe and not first_loss_seen[prefix_index]
            if first_loss:
                first_loss_seen[prefix_index] = True
                event_type = "prefix_threshold_crossing"
            elif not left_safe and right_safe:
                event_type = "safe_reentry"
            else:
                event_type = "unsafe_reentry"
            rows.append({"event_type": event_type, "prefix_n": prefix_index + 1, "sorted_index": "", "state_left": "safe" if left_safe else "unsafe", "state_right": "safe" if right_safe else "unsafe", "is_first_loss": first_loss, "notes": "equality is classified on the safe side", **common})
        for mode_index in range(args.k_max):
            left_safe = left.deltas[mode_index] <= args.threshold
            right_safe = right.deltas[mode_index] <= args.threshold
            if left_safe != right_safe:
                rows.append({"event_type": "individual_mode_crossing", "prefix_n": "", "sorted_index": mode_index + 1, "state_left": "safe" if left_safe else "unsafe", "state_right": "safe" if right_safe else "unsafe", "is_first_loss": False, "notes": "individual pass does not redefine prefix safety", **common})
    for evaluation in ordered:
        if any(evaluation.possible_multiplicity):
            rows.append(
                {
                    "event_type": "multiplicity_warning",
                    "prefix_n": "",
                    "sorted_index": ";".join(str(index + 1) for index, value in enumerate(evaluation.possible_multiplicity) if value),
                    "epsilon_left": evaluation.epsilon,
                    "epsilon_right": evaluation.epsilon,
                    "epsilon_event_estimate": evaluation.epsilon,
                    "state_left": "possible",
                    "state_right": "possible",
                    "family_signature_left": f"{evaluation.family_signature_eb}|{evaluation.family_signature_timo}",
                    "family_signature_right": f"{evaluation.family_signature_eb}|{evaluation.family_signature_timo}",
                    "related_to_reorder": False,
                    "related_to_cluster": True,
                    "is_first_loss": False,
                    "notes": "small gap, ambiguous family, or small second singular value",
                }
            )
    return rows


def reentry_count(evaluations: Sequence[Evaluation], prefix_n: int, threshold: float) -> int:
    ordered = sorted(evaluations, key=lambda item: item.epsilon)
    states = [evaluation.prefix_deltas[prefix_n - 1] <= threshold for evaluation in ordered if evaluation.resolved]
    first_unsafe = next((index for index, state in enumerate(states) if not state), None)
    if first_unsafe is None:
        return 0
    return sum(not states[index - 1] and states[index] for index in range(first_unsafe + 1, len(states)))


def build_threshold_rows(
    args: Args,
    manager: EvaluationManager,
    evaluations: Sequence[Evaluation],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    rows: list[dict[str, object]] = []
    refinement_rows: list[dict[str, object]] = []
    raw_brackets = {
        prefix_n: first_loss_bracket(evaluations, prefix_n, args.threshold)
        for prefix_n in range(1, args.k_max + 1)
    }
    smoke_interval: tuple[float, float] | None = None
    smoke_prefix: int | None = None
    if args.smoke:
        resolved_intervals = [
            (safe.epsilon, unsafe.epsilon)
            for status, safe, unsafe, _note in raw_brackets.values()
            if status == "resolved" and safe is not None and unsafe is not None
        ]
        if resolved_intervals:
            scan_midpoint = 0.5 * (args.epsilon_min + args.epsilon_max)
            smoke_interval = min(
                resolved_intervals,
                key=lambda interval: abs(0.5 * (interval[0] + interval[1]) - scan_midpoint),
            )
            smoke_prefix = max(
                prefix_n
                for prefix_n, (status, safe, unsafe, _note) in raw_brackets.items()
                if status == "resolved"
                and safe is not None
                and unsafe is not None
                and (safe.epsilon, unsafe.epsilon) == smoke_interval
            )
    for prefix_n in range(1, args.k_max + 1):
        status, safe, unsafe, note = raw_brackets[prefix_n]
        iteration_count = 0
        if (
            args.smoke
            and status == "resolved"
            and safe is not None
            and unsafe is not None
            and smoke_interval is not None
            and (
                (safe.epsilon, unsafe.epsilon) != smoke_interval
                or prefix_n != smoke_prefix
            )
        ):
            status = "smoke_not_refined"
            note += "; smoke mode refines one local event interval only"
        if status == "resolved" and safe is not None and unsafe is not None:
            safe, unsafe, local_rows, iteration_count = refine_first_loss(
                args, manager, prefix_n, safe, unsafe
            )
            refinement_rows.extend(local_rows)
            if not (safe.resolved and unsafe.resolved):
                status = "unresolved_root_quality"
                note += "; unresolved refinement point"
        estimate = 0.5 * (safe.epsilon + unsafe.epsilon) if status == "resolved" and safe is not None and unsafe is not None else float("nan")
        trigger_indices = triggering_modes(unsafe, prefix_n) if unsafe is not None and unsafe.resolved else []
        trigger_delta = max((unsafe.deltas[index - 1] for index in trigger_indices), default=float("nan")) if unsafe is not None else float("nan")
        eb_families = [unsafe.eb_families[index - 1] for index in trigger_indices] if unsafe is not None else []
        timo_families = [unsafe.timo_families[index - 1] for index in trigger_indices] if unsafe is not None else []
        nearby = [
            evaluation for evaluation in manager.evaluations.values()
            if isfinite(estimate) and abs(evaluation.epsilon - estimate) <= max(args.coarse_step, 4.0 * args.epsilon_tolerance)
            and abs(evaluation.mu) <= 1.0e-12
        ]
        reorder_near = len({(item.family_signature_eb, item.family_signature_timo) for item in nearby}) > 1
        cluster_near = any(
            isfinite(value) and value <= CLOSE_GAP_THRESHOLD
            for item in nearby for value in (*item.eb_gaps, *item.timo_gaps)
        )
        inside_primary = (
            args.primary_epsilon_min <= estimate <= args.primary_epsilon_max
            if isfinite(estimate) else False
        )
        certified_safe = (
            safe.epsilon
            if safe is not None
            and safe.resolved
            and status in {"resolved", "not_reached_in_scan_range"}
            else ""
        )
        certified_delta = (
            safe.prefix_deltas[prefix_n - 1]
            if safe is not None
            and safe.resolved
            and status in {"resolved", "not_reached_in_scan_range"}
            else ""
        )
        if status == "not_reached_in_scan_range" and safe is not None and safe.resolved:
            note += "; epsilon_certified_n is the right-censored verified scan endpoint, not a first-loss estimate"
        continuity_status = {
            "resolved": (
                "sorted_family_reorder_crossing"
                if reorder_near
                else "smooth_threshold_crossing"
            ),
            "not_reached_in_scan_range": "right_censored_safe_over_scan",
            "below_scan_range": "unsafe_at_left_scan_boundary",
            "unresolved_root_quality": "unresolved_root_quality",
            "smoke_not_refined": "smoke_not_refined",
        }.get(status, status)
        rows.append(
            {
                "prefix_n": prefix_n,
                "threshold_status": status,
                "epsilon_certified_n": certified_safe,
                "epsilon_safe_lower": certified_safe,
                "epsilon_unsafe_upper": unsafe.epsilon if status == "resolved" and unsafe is not None else "",
                "epsilon_star_estimate": estimate,
                "epsilon_star_display_5dp": presentation_round(estimate, 5),
                "epsilon_certified_floor_5dp": conservative_decimal_floor(certified_safe, 5),
                "epsilon_certified_floor_4dp": conservative_decimal_floor(certified_safe, 4),
                "bracket_width": unsafe.epsilon - safe.epsilon if status == "resolved" and safe is not None and unsafe is not None else "",
                "Delta_n_safe": certified_delta,
                "Delta_n_unsafe": unsafe.prefix_deltas[prefix_n - 1] if status == "resolved" and unsafe is not None else "",
                "triggering_sorted_indices": ";".join(str(value) for value in trigger_indices),
                "triggering_delta": trigger_delta,
                "triggering_EB_family": ";".join(eb_families),
                "triggering_Timo_family": ";".join(timo_families),
                "reorder_near_threshold": reorder_near,
                "cluster_near_threshold": cluster_near,
                "simultaneous_transition_group": "",
                "first_loss_iteration_count": iteration_count,
                "reentry_count_after_first_loss": reentry_count(evaluations, prefix_n, args.threshold),
                "inside_primary_range": inside_primary,
                "monotonicity_status": "pending",
                "continuity_status": continuity_status,
                "threshold_continuity_status": continuity_status,
                "numerical_verification_status": "pending" if status == "resolved" else "not_applicable",
                "verification_estimate": "",
                "verification_abs_difference": "",
                "verification_bracket_agreement": "",
                "verification_root_order_agreement": "",
                "epsilon_near_n": 0.999 * estimate if isfinite(estimate) else "",
                "epsilon_buffer_n": 0.99 * estimate if isfinite(estimate) else "",
                "notes": note,
            }
        )
    assign_simultaneous_groups(rows, args.epsilon_tolerance)
    check_threshold_ordering(rows, args.epsilon_tolerance)
    return rows, refinement_rows


def assign_simultaneous_groups(rows: Sequence[dict[str, object]], tolerance: float) -> None:
    resolved = [row for row in rows if isfinite(finite_float(row, "epsilon_star_estimate"))]
    groups: list[list[dict[str, object]]] = []
    for row in sorted(resolved, key=lambda item: finite_float(item, "epsilon_star_estimate")):
        estimate = finite_float(row, "epsilon_star_estimate")
        if not groups or abs(estimate - finite_float(groups[-1][0], "epsilon_star_estimate")) > tolerance:
            groups.append([row])
        else:
            groups[-1].append(row)
    group_index = 1
    for group in groups:
        if len(group) <= 1:
            continue
        label = f"simultaneous_{group_index:02d}"
        group_index += 1
        for row in group:
            row["simultaneous_transition_group"] = label


def check_threshold_ordering(rows: Sequence[dict[str, object]], tolerance: float) -> bool:
    violation = False
    previous: float | None = None
    for row in sorted(rows, key=lambda item: int(item["prefix_n"])):
        estimate = finite_float(row, "epsilon_star_estimate")
        if not isfinite(estimate):
            row["monotonicity_status"] = "not_applicable"
            continue
        if previous is not None and estimate > previous + tolerance:
            violation = True
        previous = estimate
    for row in rows:
        if isfinite(finite_float(row, "epsilon_star_estimate")):
            row["monotonicity_status"] = "violation" if violation else "pass"
            if violation:
                row["threshold_status"] = "unresolved_threshold_ordering"
                row["numerical_verification_status"] = "fail"
    return not violation


def verify_thresholds(
    args: Args,
    threshold_rows: Sequence[dict[str, object]],
) -> tuple[EvaluationManager, list[dict[str, object]]]:
    verification = EvaluationManager(
        args,
        n_candidate_roots=max(VERIFICATION_CANDIDATE_ROOTS, args.n_candidate_roots + 4),
        reuse_cache=False,
        force_recompute=True,
        scope="independent_verification",
    )
    audit_rows: list[dict[str, object]] = []
    cache: dict[tuple[float, float, int], tuple[Evaluation, Evaluation, float, bool, bool, list[dict[str, object]]]] = {}
    for row in threshold_rows:
        if row.get("threshold_status") != "resolved":
            continue
        safe_epsilon = finite_float(row, "epsilon_certified_n")
        unsafe_epsilon = finite_float(row, "epsilon_unsafe_upper")
        prefix_n = int(row["prefix_n"])
        bracket_key = (round(safe_epsilon, ROUND_DIGITS), round(unsafe_epsilon, ROUND_DIGITS), prefix_n)
        if bracket_key not in cache:
            safe = verification.evaluate(safe_epsilon, "verification_safe")
            unsafe = verification.evaluate(unsafe_epsilon, "verification_unsafe")
            verified_safe, verified_unsafe, local_rows, _iterations = refine_first_loss(
                args,
                verification,
                prefix_n,
                safe,
                unsafe,
                kind="verification_bisection",
            )
            estimate = 0.5 * (verified_safe.epsilon + verified_unsafe.epsilon)
            side_agreement = (
                verified_safe.resolved
                and verified_unsafe.resolved
                and verified_safe.prefix_deltas[prefix_n - 1] <= args.threshold
                and verified_unsafe.prefix_deltas[prefix_n - 1] > args.threshold
            )
            primary_mid = 0.5 * (safe_epsilon + unsafe_epsilon)
            verification_mid = verification.evaluate(estimate, "verification_refinement")
            primary_reference = analytic_eb_union(primary_mid, args.n_spectrum_roots)
            root_order_agreement = (
                verification_mid.resolved
                and tuple(item.family for item in primary_reference)
                == tuple(item.family for item in verification_mid.eb_references)
            )
            cache[bracket_key] = (
                verified_safe,
                verified_unsafe,
                estimate,
                side_agreement,
                root_order_agreement,
                local_rows,
            )
        verified_safe, verified_unsafe, estimate, side_agreement, root_order_agreement, local_rows = cache[bracket_key]
        audit_rows.extend(
            {**item, "prefix_n": prefix_n}
            for item in local_rows
        )
        difference = abs(estimate - finite_float(row, "epsilon_star_estimate"))
        passed = (
            difference <= args.epsilon_tolerance
            and side_agreement
            and root_order_agreement
            and verified_safe.resolved
            and verified_unsafe.resolved
        )
        row["verification_estimate"] = estimate
        row["verification_abs_difference"] = difference
        row["verification_bracket_agreement"] = side_agreement
        row["verification_root_order_agreement"] = root_order_agreement
        row["numerical_verification_status"] = "pass" if passed else "fail"
        if not passed:
            row["threshold_status"] = "unresolved_numerical_verification"
            row["notes"] = str(row.get("notes", "")) + "; independent verification failed"
    return verification, audit_rows


def mu_invariance_rows(
    args: Args,
    manager: EvaluationManager,
    threshold_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    if args.smoke:
        representative = min(
            MU_AUDIT_EPSILONS,
            key=lambda value: abs(value - 0.5 * (args.epsilon_min + args.epsilon_max)),
        )
        epsilon_values = {min(max(representative, args.epsilon_min), args.epsilon_max)}
    else:
        epsilon_values = {
            value for value in MU_AUDIT_EPSILONS
            if args.epsilon_min - 1.0e-12 <= value <= args.epsilon_max + 1.0e-12
        }
    threshold_estimates = [
        finite_float(row, "epsilon_star_estimate")
        for row in threshold_rows
        if row.get("threshold_status") == "resolved" and isfinite(finite_float(row, "epsilon_star_estimate"))
    ]
    if args.smoke:
        if threshold_estimates:
            epsilon_values.add(threshold_estimates[0])
    else:
        epsilon_values.update(threshold_estimates)
    deduplicated: list[float] = []
    for epsilon in sorted(epsilon_values):
        if not deduplicated or abs(epsilon - deduplicated[-1]) > args.epsilon_tolerance:
            deduplicated.append(epsilon)
    rows: list[dict[str, object]] = []
    for epsilon in deduplicated:
        baseline = manager.evaluate(epsilon, "mu_invariance")
        for mu in MU_AUDIT_VALUES:
            current = manager.evaluate(epsilon, "mu_invariance", mu=mu)
            for index in range(args.n_spectrum_roots):
                eb_abs, eb_rel = abs_rel_difference(current.eb_roots[index], baseline.eb_roots[index])
                timo_abs, timo_rel = abs_rel_difference(current.timo_roots[index], baseline.timo_roots[index])
                eb_pass = eb_abs <= MU_EB_ABS_TOLERANCE and eb_rel <= MU_EB_REL_TOLERANCE
                timo_pass = timo_abs <= MU_TIMO_ABS_TOLERANCE and timo_rel <= MU_TIMO_REL_TOLERANCE
                family_eb_equal = current.eb_families == baseline.eb_families
                family_timo_equal = current.timo_families == baseline.timo_families
                # Physical cross-family multiplicity comes from the retained
                # family records.  `possible_multiplicity` also contains a
                # conditioning-sensitive second-SVD diagnostic and is not a
                # canonical multiplicity identity across artificial splits.
                multiplicity_flags_equal = (
                    tuple(item.ambiguous for item in current.eb_references)
                    == tuple(item.ambiguous for item in baseline.eb_references)
                    and tuple(
                        item.cross_family_cluster
                        for item in current.factorized_timo.roots
                    )
                    == tuple(
                        item.cross_family_cluster
                        for item in baseline.factorized_timo.roots
                    )
                )
                multiplicity_groups_equal = tuple(
                    item.multiplicity_group for item in current.factorized_timo.roots
                ) == tuple(
                    item.multiplicity_group for item in baseline.factorized_timo.roots
                )
                n_equal = current.n_true == baseline.n_true
                status = (
                    "pass"
                    if current.resolved
                    and eb_pass
                    and timo_pass
                    and family_eb_equal
                    and family_timo_equal
                    and multiplicity_flags_equal
                    and multiplicity_groups_equal
                    and n_equal
                    else "fail"
                )
                rows.append(
                    {
                        "epsilon": epsilon,
                        "mu": mu,
                        "sorted_index": index + 1,
                        "Lambda_EB": current.eb_roots[index],
                        "Lambda_EB_mu0": baseline.eb_roots[index],
                        "EB_abs_difference": eb_abs,
                        "EB_rel_difference": eb_rel,
                        "Lambda_Timo": current.timo_roots[index],
                        "Lambda_Timo_mu0": baseline.timo_roots[index],
                        "Timo_abs_difference": timo_abs,
                        "Timo_rel_difference": timo_rel,
                        "N_true": current.n_true,
                        "N_true_mu0": baseline.n_true,
                        "N_true_equal": n_equal,
                        "family_order_EB_equal": family_eb_equal,
                        "family_order_Timo_equal": family_timo_equal,
                        "multiplicity_flags_equal": multiplicity_flags_equal,
                        "multiplicity_groups_Timo_equal": multiplicity_groups_equal,
                        "root_quality_status": current.quality_status,
                        "status": status,
                        "notes": "mu changes only the artificial internal joint location in the straight homogeneous rod",
                    }
                )
    return rows


def late_pass_indices(evaluation: Evaluation, threshold: float) -> list[int]:
    first_failure = next((index for index, value in enumerate(evaluation.deltas) if value > threshold), None)
    if first_failure is None:
        return []
    return [index + 1 for index, value in enumerate(evaluation.deltas) if index > first_failure and value <= threshold]


def point_rows(evaluations: Sequence[Evaluation], args: Args) -> list[dict[str, object]]:
    ordered = sorted(evaluations, key=lambda item: item.epsilon)
    signatures = [(item.family_signature_eb, item.family_signature_timo) for item in ordered]
    rows: list[dict[str, object]] = []
    for position, evaluation in enumerate(ordered):
        late = late_pass_indices(evaluation, args.threshold)
        gaps = [value for value in (*evaluation.eb_gaps[: args.k_max], *evaluation.timo_gaps[: args.k_max]) if isfinite(value)]
        reorder = (
            (position > 0 and signatures[position] != signatures[position - 1])
            or (position + 1 < len(signatures) and signatures[position] != signatures[position + 1])
        )
        rows.append(
            {
                "epsilon": evaluation.epsilon,
                "in_primary_range": args.primary_epsilon_min <= evaluation.epsilon <= args.primary_epsilon_max,
                "N_true": evaluation.n_true,
                "first_failed_mode": evaluation.n_true + 1 if evaluation.n_true < args.k_max else "",
                "late_pass_count": len(late),
                "late_pass_indices": ";".join(str(value) for value in late),
                "max_delta_first10": max(evaluation.deltas),
                "minimum_gap_first10": min(gaps) if gaps else float("nan"),
                "root11_available": len(evaluation.eb_roots) > args.k_max and isfinite(evaluation.eb_roots[args.k_max]),
                "family_signature_EB": evaluation.family_signature_eb,
                "family_signature_Timo": evaluation.family_signature_timo,
                "reorder_detected": reorder,
                "cluster_detected": any(isfinite(value) and value <= CLOSE_GAP_THRESHOLD for value in gaps),
                "possible_multiplicity": any(evaluation.possible_multiplicity),
                "root_quality_status": evaluation.quality_status,
                "cache_status": evaluation.cache_status,
                "point_type": evaluation.point_type,
                "notes": ";".join(evaluation.quality_reasons),
            }
        )
    return rows


def mode_rows(evaluations: Sequence[Evaluation], args: Args) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in sorted(evaluations, key=lambda item: item.epsilon):
        for index in range(args.k_max):
            eb_left, eb_right, eb_gap = sided_gap(evaluation.eb_gaps, index)
            t_left, t_right, t_gap = sided_gap(evaluation.timo_gaps, index)
            axial_ref, _ax_abs, _ax_rel = nearest_axial_match(evaluation.eb_roots[index], evaluation.axial_references)
            rows.append(
                {
                    "epsilon": evaluation.epsilon,
                    "sorted_index": index + 1,
                    "Lambda_EB": evaluation.eb_roots[index],
                    "Lambda_Timo": evaluation.timo_roots[index],
                    "delta_f": evaluation.deltas[index],
                    "Delta_n": evaluation.prefix_deltas[index],
                    "true_pass_10": evaluation.deltas[index] <= args.threshold,
                    "EB_family": evaluation.eb_families[index],
                    "Timo_family": evaluation.timo_families[index],
                    "family_ambiguity": evaluation.family_ambiguity[index],
                    "gap_left_EB": eb_left,
                    "gap_right_EB": eb_right,
                    "gap_EB": eb_gap,
                    "gap_left_Timo": t_left,
                    "gap_right_Timo": t_right,
                    "gap_Timo": t_gap,
                    "axial_reference_Lambda": axial_ref,
                    "EB_reference_Lambda": evaluation.eb_references[index].value,
                    "reference_match_errors": json.dumps(
                        {
                            "EB_abs": evaluation.eb_reference_abs_errors[index],
                            "EB_rel": evaluation.eb_reference_rel_errors[index],
                            "Timo_axial_abs": evaluation.timo_axial_abs_errors[index],
                            "Timo_axial_rel": evaluation.timo_axial_rel_errors[index],
                        },
                        sort_keys=True,
                    ),
                    "possible_multiplicity": evaluation.possible_multiplicity[index],
                    "smallest_singular_value_EB": evaluation.eb_singular_values[index][0],
                    "second_smallest_singular_value_EB": evaluation.eb_singular_values[index][1],
                    "smallest_singular_value_Timo": evaluation.timo_singular_values[index][0],
                    "second_smallest_singular_value_Timo": evaluation.timo_singular_values[index][1],
                    "source_warnings": ";".join((*evaluation.eb_warnings, *evaluation.timo_warnings, *evaluation.quality_reasons)),
                    "point_type": evaluation.point_type,
                }
            )
    return rows


def family_rows(evaluations: Sequence[Evaluation], args: Args) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in sorted(evaluations, key=lambda item: item.epsilon):
        for index in range(args.n_spectrum_roots):
            axial_ref, timo_ax_abs, timo_ax_rel = nearest_axial_match(
                evaluation.timo_roots[index], evaluation.axial_references
            )
            _eb_left, _eb_right, eb_gap = sided_gap(evaluation.eb_gaps, index)
            _t_left, _t_right, timo_gap = sided_gap(evaluation.timo_gaps, index)
            reference = evaluation.eb_references[index]
            rows.append(
                {
                    "epsilon": evaluation.epsilon,
                    "sorted_index": index + 1,
                    "Lambda_EB_general": evaluation.eb_roots[index],
                    "Lambda_EB_analytic_union": reference.value,
                    "EB_reference_family": reference.family,
                    "EB_reference_family_index": reference.family_index,
                    "EB_reference_abs_error": evaluation.eb_reference_abs_errors[index],
                    "EB_reference_rel_error": evaluation.eb_reference_rel_errors[index],
                    "Lambda_axial_reference": axial_ref,
                    "Lambda_Timo_factorized": evaluation.timo_roots[index],
                    "Lambda_Timo_general_6x6": (
                        evaluation.raw_timo_roots[index]
                        if index < len(evaluation.raw_timo_roots)
                        else float("nan")
                    ),
                    "Timo_family": evaluation.timo_families[index],
                    "Timo_factorized_family_exact": evaluation.factorized_timo.roots[index].family,
                    "Timo_axial_abs_error": timo_ax_abs,
                    "Timo_axial_rel_error": timo_ax_rel,
                    "family_ambiguity": evaluation.family_ambiguity[index],
                    "adjacent_gap_EB": eb_gap,
                    "adjacent_gap_Timo": timo_gap,
                    "possible_multiplicity": evaluation.possible_multiplicity[index],
                    "point_type": evaluation.point_type,
                    "notes": "corrected baseline uses the exact family union; the raw general 6x6 value is retained only for completeness comparison",
                }
            )
    return rows


def local_eb_svd_confirmation(
    evaluation: Evaluation,
    index: int,
) -> tuple[float, float, bool]:
    root = evaluation.eb_roots[index]
    neighbor_gaps = []
    if index > 0:
        neighbor_gaps.append(root - evaluation.eb_roots[index - 1])
    if index + 1 < len(evaluation.eb_roots):
        neighbor_gaps.append(evaluation.eb_roots[index + 1] - root)
    positive = [gap for gap in neighbor_gaps if gap > 0.0]
    radius = min(
        LOCAL_FULL_SVD_RADIUS,
        0.45 * min(positive) if positive else LOCAL_FULL_SVD_RADIUS,
    )
    radius = max(radius, 5.0e-7)
    exact_sigma, _second = eb_full_matrix_singular_values(
        root, evaluation.mu, evaluation.epsilon
    )
    result = minimize_scalar(
        lambda value: eb_full_matrix_singular_values(
            float(value), evaluation.mu, evaluation.epsilon
        )[0],
        bounds=(max(root_workflow.ROOT_SCAN_START, root - radius), root + radius),
        method="bounded",
        options={"xatol": 1.0e-12, "maxiter": 160},
    )
    candidates = [(root, exact_sigma)]
    if result.success and isfinite(float(result.fun)):
        candidates.append((float(result.x), float(result.fun)))
    best_root, best_sigma = min(candidates, key=lambda pair: pair[1])
    return best_root, best_sigma, best_sigma <= factorized.FULL_MATRIX_SVD_ACCEPT


def factorized_spectrum_rows(
    evaluations: Sequence[Evaluation],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in sorted(evaluations, key=lambda item: (item.epsilon, item.mu)):
        eb_components = tuple(
            factorized.FamilyRoot(
                value=item.value,
                family=item.family,
                family_index=item.family_index,
                detection_source=(
                    "exact_axial_formula"
                    if item.family == "axial"
                    else "analytic_fixed_fixed_EB"
                ),
                block_sigma_min=evaluation.eb_singular_values[index][0],
                multiplicity_group=(
                    f"eb_cross_family_{index + 1:03d}" if item.ambiguous else ""
                ),
                cross_family_cluster=item.ambiguous,
                cluster_gap=(
                    sided_gap(evaluation.eb_gaps, index)[2]
                    if item.ambiguous
                    else float("nan")
                ),
            )
            for index, item in enumerate(evaluation.eb_references)
        )
        eb_matches = match_general_roots(eb_components, evaluation.raw_eb_roots)
        for index, (item, match) in enumerate(zip(eb_components, eb_matches)):
            general_index, general_value, _difference = match
            missing = general_index is None
            if missing:
                _local_root, _local_sigma, recovered = local_eb_svd_confirmation(
                    evaluation, index
                )
            else:
                recovered = True
            full_smallest, full_second = eb_full_matrix_singular_values(
                item.value, evaluation.mu, evaluation.epsilon
            )
            matrix = root_workflow.assemble_clamped_coupled_matrix_eta(
                item.value, 0.0, evaluation.mu, evaluation.epsilon, 0.0
            )
            rows.append(
                {
                    "epsilon": evaluation.epsilon,
                    "mu": evaluation.mu,
                    "model": "Euler-Bernoulli",
                    "sorted_index": index + 1,
                    "sorted_position": index + 1,
                    "Lambda_factorized": item.value,
                    "family": item.family,
                    "family_index": item.family_index,
                    "detection_source": item.detection_source,
                    "block_det": float(np.linalg.det(matrix)),
                    "block_sigma_min": item.block_sigma_min,
                    "full_matrix_sigma_min": full_smallest,
                    "full_matrix_second_sigma": full_second,
                    "off_block_max_abs": float("nan"),
                    "multiplicity_group": item.multiplicity_group,
                    "multiplicity_group_id": item.multiplicity_group,
                    "within_family_duplicate": item.within_family_duplicate,
                    "cross_family_cluster": item.cross_family_cluster,
                    "cluster_gap": item.cluster_gap,
                    "cross_family_gap": item.cluster_gap,
                    "raw_general_match": general_value,
                    "raw_general_missing": missing,
                    "local_general_recovery_status": (
                        "confirmed_by_local_full_matrix_SVD" if missing and recovered else "not_needed"
                    ),
                    "cache_status": evaluation.cache_status,
                    "algorithm_version": factorized.ALGORITHM_VERSION,
                    "root_quality_status": (
                        "pass"
                        if full_smallest <= factorized.FULL_MATRIX_SVD_ACCEPT and recovered
                        else "fail"
                    ),
                    "quality_status": (
                        "missing_required_axial_root_in_general_scan"
                        if missing and item.family == "axial"
                        else ("missing_in_general_sign_scan" if missing else "pass")
                    ),
                    "notes": "exact analytic EB family union; general 6x6 retained as completeness audit",
                }
            )

        timo_matches = match_general_roots(
            evaluation.factorized_timo.roots, evaluation.raw_timo_roots
        )
        for index, item in enumerate(evaluation.factorized_timo.roots, start=1):
            general_index, general_value, _difference = timo_matches[index - 1]
            missing = general_index is None
            if missing:
                _local_root, _local_sigma, recovered = local_full_svd_confirmation(
                    evaluation, index - 1
                )
            else:
                recovered = True
            matrix, _bending, _axial, _warnings = factorized.straight_blocks(
                item.value, evaluation.mu, evaluation.epsilon
            )
            local_off_block = factorized.off_block_max_abs(matrix)
            full_smallest, full_second = factorized.full_matrix_singular_values(
                item.value, evaluation.mu, evaluation.epsilon
            )
            status = (
                "pass"
                if local_off_block <= factorized.BLOCK_OFFDIAGONAL_ABS_TOL
                and item.block_sigma_min <= (
                    factorized.AXIAL_BLOCK_SVD_ACCEPT
                    if item.family == "axial"
                    else factorized.BENDING_SVD_ACCEPT
                )
                and full_smallest <= factorized.FULL_MATRIX_SVD_ACCEPT
                and recovered
                else "fail"
            )
            block_det = (
                factorized.axial_block_det(item.value, evaluation.mu, evaluation.epsilon)
                if item.family == "axial"
                else factorized.bending_block_det(item.value, evaluation.mu, evaluation.epsilon)
            )
            rows.append(
                {
                    "epsilon": evaluation.epsilon,
                    "mu": evaluation.mu,
                    "model": "Timoshenko",
                    "sorted_index": index,
                    "sorted_position": index,
                    "Lambda_factorized": item.value,
                    "family": item.family,
                    "family_index": item.family_index,
                    "detection_source": item.detection_source,
                    "block_det": block_det,
                    "block_sigma_min": item.block_sigma_min,
                    "full_matrix_sigma_min": full_smallest,
                    "full_matrix_second_sigma": full_second,
                    "off_block_max_abs": local_off_block,
                    "multiplicity_group": item.multiplicity_group,
                    "multiplicity_group_id": item.multiplicity_group,
                    "within_family_duplicate": item.within_family_duplicate,
                    "cross_family_cluster": item.cross_family_cluster,
                    "cluster_gap": item.cluster_gap,
                    "cross_family_gap": item.cluster_gap,
                    "raw_general_match": general_value,
                    "raw_general_missing": missing,
                    "local_general_recovery_status": (
                        "confirmed_by_local_full_matrix_SVD" if missing and recovered else "not_needed"
                    ),
                    "cache_status": evaluation.cache_status,
                    "algorithm_version": factorized.ALGORITHM_VERSION,
                    "root_quality_status": status,
                    "quality_status": (
                        "missing_required_axial_root_in_general_scan"
                        if missing and item.family == "axial"
                        else ("missing_in_general_sign_scan" if missing else "pass")
                    ),
                    "notes": "exact beta=0 eta=0 block extraction from the unchanged 6x6 Timoshenko matrix",
                }
            )
    return rows


def axial_completeness_rows(
    evaluations: Sequence[Evaluation],
    count: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for evaluation in sorted(evaluations, key=lambda item: (item.epsilon, item.mu)):
        matches = match_general_roots(
            evaluation.factorized_timo.roots, evaluation.raw_timo_roots
        )
        for sorted_index, item in enumerate(evaluation.factorized_timo.roots, start=1):
            if item.family != "axial":
                continue
            axial_index = item.family_index
            value = item.value
            block_smallest, block_second = factorized.axial_block_singular_values(
                float(value), evaluation.mu, evaluation.epsilon
            )
            full_smallest, full_second = factorized.full_matrix_singular_values(
                float(value), evaluation.mu, evaluation.epsilon
            )
            _general_index, nearest_general, difference = matches[sorted_index - 1]
            raw_present = isfinite(nearest_general)
            block_det = factorized.axial_block_det(
                float(value), evaluation.mu, evaluation.epsilon
            )
            passed = (
                block_smallest <= factorized.AXIAL_BLOCK_SVD_ACCEPT
                and full_smallest <= factorized.FULL_MATRIX_SVD_ACCEPT
            )
            status = (
                "fail"
                if not passed
                else (
                    "pass"
                    if raw_present
                    else "missing_required_axial_root_in_general_scan"
                )
            )
            rows.append(
                {
                    "epsilon": evaluation.epsilon,
                    "mu": evaluation.mu,
                    "axial_index": axial_index,
                    "Lambda_exact": value,
                    "Lambda_axial_exact": value,
                    "axial_block_det": block_det,
                    "axial_block_sigma_min": block_smallest,
                    "axial_block_second_sigma": block_second,
                    "axial_block_sigma_ratio": (
                        block_smallest / block_second
                        if isfinite(block_second) and block_second > 0.0
                        else float("nan")
                    ),
                    "full_matrix_sigma_min": full_smallest,
                    "full_matrix_sigma_ratio": (
                        full_smallest / full_second
                        if isfinite(full_second) and full_second > 0.0
                        else float("nan")
                    ),
                    "expected_within_first12": True,
                    "present_in_factorized_first12": True,
                    "factorized_present": True,
                    "raw_general_present": raw_present,
                    "full_matrix_verified": passed,
                    "factorized_sorted_index": sorted_index,
                    "sorted_position_factorized": sorted_index,
                    "nearest_raw_general_root": nearest_general,
                    "difference": difference,
                    "verification_status": "pass" if passed else "fail",
                    "status": status,
                    "algorithm_version": factorized.ALGORITHM_VERSION,
                    "notes": "exact axial root retained in the requested sorted prefix; raw scan status is independent evidence",
                }
            )
    del count
    return rows


def general_match_tolerance(value: float) -> float:
    return GENERAL_MATCH_ABS_TOLERANCE + GENERAL_MATCH_REL_TOLERANCE * abs(float(value))


def match_general_roots(
    factorized_roots: Sequence[factorized.FamilyRoot],
    general_roots: Sequence[float],
) -> list[tuple[int | None, float, float]]:
    finite = [
        (index, float(value))
        for index, value in enumerate(general_roots, start=1)
        if isfinite(float(value))
    ]
    unused = set(range(len(finite)))
    matches: list[tuple[int | None, float, float]] = []
    for item in factorized_roots:
        tolerance = general_match_tolerance(item.value)
        candidates = [
            position
            for position in unused
            if abs(finite[position][1] - item.value) <= tolerance
        ]
        if not candidates:
            matches.append((None, float("nan"), float("nan")))
            continue
        position = min(candidates, key=lambda value: abs(finite[value][1] - item.value))
        unused.remove(position)
        general_index, general_value = finite[position]
        matches.append((general_index, general_value, abs(general_value - item.value)))
    return matches


def local_full_svd_confirmation(
    evaluation: Evaluation,
    index: int,
) -> tuple[float, float, bool]:
    item = evaluation.factorized_timo.roots[index]
    neighbor_gaps = []
    if index > 0:
        neighbor_gaps.append(item.value - evaluation.factorized_timo.roots[index - 1].value)
    if index + 1 < len(evaluation.factorized_timo.roots):
        neighbor_gaps.append(evaluation.factorized_timo.roots[index + 1].value - item.value)
    positive_gaps = [value for value in neighbor_gaps if value > 0.0]
    radius = min(
        LOCAL_FULL_SVD_RADIUS,
        0.45 * min(positive_gaps) if positive_gaps else LOCAL_FULL_SVD_RADIUS,
    )
    radius = max(radius, 5.0e-7)
    left = max(root_workflow.ROOT_SCAN_START, item.value - radius)
    right = item.value + radius
    exact_sigma, _exact_second = factorized.full_matrix_singular_values(
        item.value, evaluation.mu, evaluation.epsilon
    )
    result = minimize_scalar(
        lambda value: factorized.full_matrix_singular_values(
            float(value), evaluation.mu, evaluation.epsilon
        )[0],
        bounds=(left, right),
        method="bounded",
        options={"xatol": 1.0e-12, "maxiter": 160},
    )
    candidates = [(item.value, exact_sigma)]
    if result.success and isfinite(float(result.fun)):
        candidates.append((float(result.x), float(result.fun)))
    best_root, best_sigma = min(candidates, key=lambda pair: pair[1])
    return best_root, best_sigma, best_sigma <= factorized.FULL_MATRIX_SVD_ACCEPT


def general_vs_factorized_rows(
    managers: Sequence[EvaluationManager],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for manager in managers:
        for evaluation in sorted(manager.evaluations.values(), key=lambda item: (item.epsilon, item.mu)):
            eb_components = tuple(
                factorized.FamilyRoot(
                    value=reference.value,
                    family=reference.family,
                    family_index=reference.family_index,
                    detection_source="exact_analytic_union",
                    block_sigma_min=evaluation.eb_singular_values[index][0],
                )
                for index, reference in enumerate(evaluation.eb_references)
            )
            model_payloads = (
                (
                    "Euler-Bernoulli",
                    eb_components,
                    evaluation.raw_eb_roots,
                    evaluation.raw_eb_warnings,
                ),
                (
                    "Timoshenko",
                    evaluation.factorized_timo.roots,
                    evaluation.raw_timo_roots,
                    evaluation.raw_timo_warnings,
                ),
            )
            for model, components, general_roots, general_warnings in model_payloads:
                matches = match_general_roots(components, general_roots)
                for index, (item, match) in enumerate(zip(components, matches)):
                    general_index, general_value, difference = match
                    missing = general_index is None
                    if model == "Timoshenko":
                        if missing:
                            local_root, local_sigma, confirmed = local_full_svd_confirmation(
                                evaluation, index
                            )
                        else:
                            local_root = item.value
                            local_sigma, _second = factorized.full_matrix_singular_values(
                                item.value, evaluation.mu, evaluation.epsilon
                            )
                            confirmed = local_sigma <= factorized.FULL_MATRIX_SVD_ACCEPT
                    else:
                        if missing:
                            local_root, local_sigma, confirmed = local_eb_svd_confirmation(
                                evaluation, index
                            )
                        else:
                            local_root = item.value
                            local_sigma, _second = eb_full_matrix_singular_values(
                                item.value, evaluation.mu, evaluation.epsilon
                            )
                            confirmed = local_sigma <= factorized.FULL_MATRIX_SVD_ACCEPT
                    missing_status = (
                        "missing_required_axial_root_in_general_scan"
                        if missing and item.family == "axial"
                        else "missing_in_general_sign_scan"
                    )
                    rows.append(
                        {
                            "scope": manager.scope,
                            "epsilon": evaluation.epsilon,
                            "mu": evaluation.mu,
                            "model": model,
                            "factorized_sorted_index": index + 1,
                            "Lambda_factorized": item.value,
                            "factorized_family": item.family,
                            "factorized_family_index": item.family_index,
                            "present_in_raw_general": not missing,
                            "nearest_general_root": general_value,
                            "Lambda_general_6x6": general_value,
                            "general_sorted_index": general_index or "",
                            "abs_difference": difference,
                            "match_tolerance": general_match_tolerance(item.value),
                            "match_status": missing_status if missing else "matched",
                            "missing_general_root": missing,
                            "raw_general_missing_root": missing,
                            "same_family_when_known": True if not missing else "unknown",
                            "general_sigma_min_at_factorized_root": (
                                factorized.full_matrix_singular_values(
                                    item.value, evaluation.mu, evaluation.epsilon
                                )[0]
                                if model == "Timoshenko"
                                else eb_full_matrix_singular_values(
                                    item.value, evaluation.mu, evaluation.epsilon
                                )[0]
                            ),
                            "local_full_SVD_root": local_root,
                            "local_full_SVD_sigma_min": local_sigma,
                            "full_matrix_SVD_confirmation": confirmed,
                            "general_cache_status": evaluation.general_cache_status,
                            "general_root_warnings": ";".join(general_warnings),
                            "notes": (
                                "missing general root independently confirmed by local full-6x6 SVD refinement"
                                if missing and confirmed
                                else "one-to-one multiplicity-preserving general/factorized comparison"
                            ),
                        }
                    )
    return rows


def legacy_threshold_comparison_rows(
    args: Args,
    corrected_rows: Sequence[Mapping[str, object]],
    general_comparison: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    legacy_path = (
        args.output_dir
        / "legacy_pre_factorized_root_fix"
        / "baseline_critical_prefix_thresholds.csv"
    )
    if not legacy_path.exists():
        return []
    try:
        legacy_source = str(legacy_path.relative_to(REPO_ROOT))
    except ValueError:
        legacy_source = str(legacy_path)
    legacy_by_prefix = {
        int(float(row["prefix_n"])): row for row in read_csv(legacy_path)
    }
    rows: list[dict[str, object]] = []
    for corrected in sorted(corrected_rows, key=lambda row: int(row["prefix_n"])):
        prefix_n = int(corrected["prefix_n"])
        legacy = legacy_by_prefix.get(prefix_n, {})
        legacy_status = str(legacy.get("threshold_status", ""))
        corrected_status = str(corrected.get("threshold_status", ""))
        legacy_estimate = finite_float(legacy, "epsilon_star_estimate")
        corrected_estimate = finite_float(corrected, "epsilon_star_estimate")
        difference = (
            abs(legacy_estimate - corrected_estimate)
            if isfinite(legacy_estimate) and isfinite(corrected_estimate)
            else float("nan")
        )
        if "unresolved" in corrected_status:
            comparison_status = "still_unresolved"
        elif "unresolved" in legacy_status and "unresolved" not in corrected_status:
            comparison_status = "previously_unresolved_now_resolved"
        elif legacy_status == "not_reached_in_scan_range" and corrected_status == legacy_status:
            comparison_status = "not_reached"
        elif isfinite(difference) and difference <= args.epsilon_tolerance:
            comparison_status = "unchanged_within_tolerance"
        elif corrected_status == "not_reached_in_scan_range":
            comparison_status = "not_reached"
        else:
            comparison_status = "corrected_due_to_missing_root"
        nearby_missing = [
            audit_row
            for audit_row in general_comparison
            if bool_value(audit_row.get("missing_general_root"))
            and isfinite(legacy_estimate)
            and abs(finite_float(audit_row, "epsilon") - legacy_estimate)
            <= args.coarse_step + args.epsilon_tolerance
        ]
        missing_families = sorted(
            {str(item.get("factorized_family")) for item in nearby_missing}
        )
        relative_difference = (
            difference / abs(legacy_estimate)
            if isfinite(difference) and abs(legacy_estimate) > 0.0
            else float("nan")
        )
        if comparison_status == "corrected_due_to_missing_root":
            correction_reason = (
                "raw general sign scan omitted a factorized root near the legacy bracket"
            )
        elif comparison_status == "not_reached":
            correction_reason = "first loss was not reached in either scan; no finite bracket exists"
        else:
            correction_reason = "corrected and legacy estimates agree within epsilon tolerance"
        rows.append(
            {
                "prefix_n": prefix_n,
                "legacy_status": legacy_status,
                "legacy_threshold_status": legacy_status,
                "corrected_status": corrected_status,
                "corrected_threshold_status": corrected_status,
                "legacy_epsilon_safe_lower": legacy.get("epsilon_certified_n", ""),
                "legacy_epsilon_certified_n": legacy.get("epsilon_certified_n", ""),
                "corrected_epsilon_certified_n": corrected.get("epsilon_certified_n", ""),
                "corrected_epsilon_safe_lower": corrected.get("epsilon_certified_n", ""),
                "legacy_epsilon_star_estimate": legacy.get("epsilon_star_estimate", ""),
                "corrected_epsilon_star_estimate": corrected.get("epsilon_star_estimate", ""),
                "absolute_estimate_change": difference,
                "absolute_difference": difference,
                "relative_difference": relative_difference,
                "legacy_trigger_mode": legacy.get("triggering_sorted_indices", ""),
                "corrected_trigger_mode": corrected.get("triggering_sorted_indices", ""),
                "missing_root_in_legacy_bracket": bool(nearby_missing),
                "missing_root_family": ";".join(missing_families),
                "correction_reason": correction_reason,
                "scientific_status": comparison_status,
                "comparison_status": comparison_status,
                "legacy_source": legacy_source,
                "algorithm_version": factorized.ALGORITHM_VERSION,
                "notes": "legacy values are preserved for provenance and are not used by the corrected certificate",
            }
        )
    return rows


def regression_close_root_rows(
    manager: EvaluationManager,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_id, nominal, expected_bending, expected_axial in REGRESSION_CASES:
        for offset in (-REGRESSION_EPSILON_OFFSET, 0.0, REGRESSION_EPSILON_OFFSET):
            epsilon = nominal + offset
            evaluation = manager.evaluate(epsilon, "close_root_regression")
            bending_candidates = [
                (index, item)
                for index, item in enumerate(evaluation.factorized_timo.roots, start=1)
                if item.family == "bending_Timoshenko"
            ]
            axial_candidates = [
                (index, item)
                for index, item in enumerate(evaluation.factorized_timo.roots, start=1)
                if item.family == "axial"
            ]
            if not bending_candidates or not axial_candidates:
                rows.append(
                    {
                        "case_id": case_id,
                        "epsilon_nominal": nominal,
                        "epsilon_offset": offset,
                        "epsilon": epsilon,
                        "bending_expected_window_center": expected_bending,
                        "axial_expected_window_center": expected_axial,
                        "Lambda_bending_factorized": "",
                        "Lambda_axial_exact": "",
                        "pair_gap": "",
                        "bending_sorted_index": "",
                        "axial_sorted_index": "",
                        "both_roots_present": False,
                        "cross_family_cluster": False,
                        "general_bending_present": False,
                        "general_axial_present": False,
                        "general_missing_root_count": "",
                        "full_matrix_SVD_bending": "",
                        "full_matrix_SVD_axial": "",
                        "status": "unavailable_synthetic_evaluation",
                        "algorithm_version": factorized.ALGORITHM_VERSION,
                        "notes": "regression families unavailable in injected synthetic evaluation",
                    }
                )
                continue
            bending_index, bending_item = min(
                bending_candidates,
                key=lambda pair: abs(pair[1].value - expected_bending),
            )
            axial_index, axial_item = min(
                axial_candidates,
                key=lambda pair: abs(pair[1].value - expected_axial),
            )
            matches = match_general_roots(
                evaluation.factorized_timo.roots,
                evaluation.raw_timo_roots,
            )
            general_bending = matches[bending_index - 1][0] is not None
            general_axial = matches[axial_index - 1][0] is not None
            missing_count = sum(match[0] is None for match in matches)
            bending_sigma, _ = factorized.full_matrix_singular_values(
                bending_item.value, evaluation.mu, evaluation.epsilon
            )
            axial_sigma, _ = factorized.full_matrix_singular_values(
                axial_item.value, evaluation.mu, evaluation.epsilon
            )
            both_present = (
                abs(bending_item.value - expected_bending) <= 0.05
                and abs(axial_item.value - expected_axial) <= 0.05
                and bending_index != axial_index
            )
            passed = (
                both_present
                and bending_sigma <= factorized.FULL_MATRIX_SVD_ACCEPT
                and axial_sigma <= factorized.FULL_MATRIX_SVD_ACCEPT
            )
            rows.append(
                {
                    "case_id": case_id,
                    "epsilon_nominal": nominal,
                    "epsilon_offset": offset,
                    "epsilon": epsilon,
                    "bending_expected_window_center": expected_bending,
                    "axial_expected_window_center": expected_axial,
                    "Lambda_bending_factorized": bending_item.value,
                    "Lambda_axial_exact": axial_item.value,
                    "pair_gap": abs(axial_item.value - bending_item.value),
                    "bending_sorted_index": bending_index,
                    "axial_sorted_index": axial_index,
                    "both_roots_present": both_present,
                    "cross_family_cluster": bending_item.cross_family_cluster and axial_item.cross_family_cluster,
                    "general_bending_present": general_bending,
                    "general_axial_present": general_axial,
                    "general_missing_root_count": missing_count,
                    "full_matrix_SVD_bending": bending_sigma,
                    "full_matrix_SVD_axial": axial_sigma,
                    "status": "pass" if passed else "fail",
                    "algorithm_version": factorized.ALGORITHM_VERSION,
                    "notes": "expected values define regression windows only; no root is hard-coded",
                }
            )
    return rows


def quality_rows(
    managers: Sequence[EvaluationManager], args: Args
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for manager in managers:
        for evaluation in sorted(manager.evaluations.values(), key=lambda item: (item.epsilon, item.mu)):
            max_abs = max((value for value in evaluation.eb_reference_abs_errors if isfinite(value)), default=float("nan"))
            max_rel = max((value for value in evaluation.eb_reference_rel_errors if isfinite(value)), default=float("nan"))
            mismatch_count = sum(
                not reference_match_ok(evaluation.eb_roots[index], evaluation.eb_references[index].value)
                for index in range(args.n_spectrum_roots)
            )
            rows.append(
                {
                    "scope": manager.scope,
                    "epsilon": evaluation.epsilon,
                    "mu": evaluation.mu,
                    "n_candidate_roots": manager.n_candidate_roots or args.n_candidate_roots,
                    "EB_requested_root_count": manager.n_candidate_roots or args.n_candidate_roots,
                    "EB_found_root_count": evaluation.eb_found,
                    "Timo_requested_root_count": args.n_spectrum_roots,
                    "Timo_found_root_count": evaluation.timo_found,
                    "EB_root_warning": ";".join(evaluation.eb_warnings),
                    "Timo_root_warning": ";".join(evaluation.timo_warnings),
                    "candidate_boundary_warning": evaluation.candidate_boundary_warning,
                    "cache_status": evaluation.cache_status,
                    "EB_retry_attempted": evaluation.eb_retry_attempted,
                    "Timo_retry_attempted": evaluation.timo_retry_attempted,
                    "SVD_recovery_attempted": evaluation.svd_recovery_attempted,
                    "SVD_recovery_changed_roots": evaluation.svd_recovery_changed_roots,
                    "root11_available": len(evaluation.eb_roots) > args.k_max and isfinite(evaluation.eb_roots[args.k_max]),
                    "strictly_increasing_EB": roots_strictly_increasing(evaluation.eb_roots, args.n_spectrum_roots),
                    "strictly_increasing_Timo": roots_strictly_increasing(evaluation.timo_roots, args.n_spectrum_roots),
                    "multiplicity_preserving_order_EB": roots_multiplicity_preserving_order(evaluation.eb_roots, args.n_spectrum_roots),
                    "multiplicity_preserving_order_Timo": roots_multiplicity_preserving_order(evaluation.timo_roots, args.n_spectrum_roots),
                    "analytic_EB_union_max_abs_error": max_abs,
                    "analytic_EB_union_max_rel_error": max_rel,
                    "analytic_EB_union_mismatch_count": mismatch_count,
                    "possible_multiplicity_count": sum(evaluation.possible_multiplicity),
                    "stored_recomputed_delta_mismatch_count": 0,
                    "root_quality_status": evaluation.quality_status,
                    "quality_reasons": ";".join(evaluation.quality_reasons),
                    "point_type": evaluation.point_type,
                    "notes": "corrected source is the exact straight-family union; raw general 6x6 roots remain independent completeness evidence",
                }
            )
    return rows


def operation_rows(
    args: Args,
    primary: EvaluationManager,
    verification: EvaluationManager,
    baseline_evaluations: Sequence[Evaluation],
    refinement_rows: Sequence[Mapping[str, object]],
    elapsed_seconds: float,
) -> list[dict[str, object]]:
    managers = (primary, verification)
    root_calls_eb = sum(manager.runtime.eb_root_calls for manager in managers)
    root_calls_timo = sum(manager.runtime.timo_root_calls for manager in managers)
    hits = sum(manager.runtime.root_cache_hits for manager in managers)
    misses = sum(manager.runtime.root_cache_misses for manager in managers)
    svd_calls = sum(manager.svd_recovery_calls for manager in managers)
    analytic_calls = sum(manager.analytic_reference_evaluations for manager in managers)
    factorized_evaluations = [
        evaluation for manager in managers for evaluation in manager.evaluations.values()
    ]
    factorized_hits = sum(evaluation.cache_status == "hit" for evaluation in factorized_evaluations)
    factorized_misses = sum(evaluation.cache_status != "hit" for evaluation in factorized_evaluations)
    factorized_det_evaluations = sum(
        evaluation.factorized_timo.bending_search.determinant_evaluations
        for evaluation in factorized_evaluations
        if evaluation.cache_status != "hit"
    )
    factorized_svd_refinements = sum(
        evaluation.factorized_timo.bending_search.svd_refinements
        for evaluation in factorized_evaluations
        if evaluation.cache_status != "hit"
    )
    unique_eps = len({round(item.epsilon, ROUND_DIGITS) for item in baseline_evaluations})
    primary_mu0 = [
        evaluation
        for (_epsilon, mu), evaluation in primary.evaluations.items()
        if abs(mu) <= 1.0e-12
    ]
    verification_mu0 = [
        evaluation
        for (_epsilon, mu), evaluation in verification.evaluations.items()
        if abs(mu) <= 1.0e-12
    ]
    coarse_count = sum("coarse" in evaluation.point_types for evaluation in primary_mu0)
    midpoint_count = sum("midpoint" in evaluation.point_types for evaluation in primary_mu0)
    refinement_labels = {"event_refinement", "threshold_refinement"}
    primary_refinement_count = sum(bool(evaluation.point_types.intersection(refinement_labels)) for evaluation in primary_mu0)
    verification_refinement_count = len(verification_mu0)
    epsilon_bisections = sum(str(row.get("refinement_kind", "")).endswith("bisection") and int(float(row.get("iteration", 0))) > 0 for row in refinement_rows)
    measured = [
        ("epsilon points evaluated", unique_eps, "measured", "unique baseline mu=0 epsilon points"),
        ("coarse epsilon evaluations", coarse_count, "measured", "unique initial coarse-grid root evaluations"),
        ("midpoint epsilon evaluations", midpoint_count, "measured", "unique mandatory half-step root evaluations"),
        ("epsilon refinement evaluations", primary_refinement_count + verification_refinement_count, "measured", "unique event, threshold, and independent-verification root evaluations"),
        ("EB root solves", root_calls_eb, "measured_wrapper_counter", "includes explicit existing SVD-recovery reruns"),
        ("Timoshenko root solves", root_calls_timo, "measured_wrapper_counter", "general 6x6 workflow"),
        ("cache hits", hits, "measured_wrapper_counter", "both theories counted separately"),
        ("cache misses", misses, "measured_wrapper_counter", "both theories counted separately"),
        ("SVD recovery calls", svd_calls, "measured", "existing global helper plus audit-local calls to the existing 6x6 SVD minimizer"),
        ("epsilon bisection evaluations", epsilon_bisections, "measured", "adaptive epsilon threshold refinement"),
        ("analytic reference evaluations", analytic_calls, "measured", "axial plus fixed-fixed EB union audits"),
        ("factorized spectrum cache hits", factorized_hits, "measured", factorized.ALGORITHM_VERSION),
        ("factorized spectrum cache misses/recomputes", factorized_misses, "measured", factorized.ALGORITHM_VERSION),
        ("factorized bending determinant evaluations", factorized_det_evaluations, "measured", "exact extracted 4x4 block"),
        ("factorized bending SVD refinements", factorized_svd_refinements, "measured", "audit-local complementary minima"),
        ("wall-clock seconds", elapsed_seconds, "auxiliary", "not used as the primary cost metric"),
    ]
    unavailable = [
        ("determinant/matrix evaluations", "unavailable"),
        ("root brackets", "unavailable"),
        ("root Brent function evaluations", "unavailable"),
    ]
    rows = [
        {"scope": "complete_run", "metric": metric, "value": value, "measurement_status": status, "notes": notes}
        for metric, value, status, notes in measured
    ]
    rows.extend(
        {"scope": "complete_run", "metric": metric, "value": value, "measurement_status": "unavailable", "notes": "existing wrappers do not expose this primitive count"}
        for metric, value in unavailable
    )
    del args
    return rows


def make_plots(
    output_dir: Path,
    points: Sequence[Mapping[str, object]],
    modes: Sequence[Mapping[str, object]],
    thresholds: Sequence[Mapping[str, object]],
    families: Sequence[Mapping[str, object]],
    *,
    threshold: float,
    primary_min: float,
    primary_max: float,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    point_order = sorted(points, key=lambda row: finite_float(row, "epsilon"))
    epsilon = np.array([finite_float(row, "epsilon") for row in point_order], dtype=float)
    n_true = np.array([finite_float(row, "N_true") for row in point_order], dtype=float)
    estimates = [finite_float(row, "epsilon_star_estimate") for row in thresholds if isfinite(finite_float(row, "epsilon_star_estimate"))]
    paths: list[Path] = []

    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.axvspan(primary_min, primary_max, color="#dcefdc", alpha=0.5, label="primary range")
    ax.step(epsilon, n_true, where="post", color="#1f4e79", label="sampled N_true")
    ax.scatter(epsilon, n_true, s=10, color="#1f4e79")
    for value in sorted(set(estimates)):
        ax.axvline(value, color="#aa3333", alpha=0.35, linewidth=0.8)
    ax.set(xlabel="epsilon", ylabel="N_true", ylim=(-0.3, 10.5))
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    path = output_dir / "baseline_N_true_vs_epsilon.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for index in range(1, DEFAULT_K_MAX + 1):
        rows = sorted((row for row in modes if int(float(row.get("sorted_index", 0))) == index), key=lambda row: finite_float(row, "epsilon"))
        ax.plot([finite_float(row, "epsilon") for row in rows], [finite_float(row, "delta_f") for row in rows], linewidth=1.0, label=f"k={index}")
    ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0)
    ax.axvspan(primary_min, primary_max, color="#dcefdc", alpha=0.25)
    ax.set(xlabel="epsilon", ylabel="delta_f")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=7)
    fig.tight_layout()
    path = output_dir / "baseline_delta_f_by_mode.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for index in range(1, DEFAULT_K_MAX + 1):
        rows = sorted((row for row in modes if int(float(row.get("sorted_index", 0))) == index), key=lambda row: finite_float(row, "epsilon"))
        ax.plot([finite_float(row, "epsilon") for row in rows], [finite_float(row, "Delta_n") for row in rows], linewidth=1.0, label=f"n={index}")
    ax.axhline(threshold, color="black", linestyle="--", linewidth=1.0)
    for value in estimates:
        ax.scatter([value], [threshold], s=18, color="#aa3333", zorder=4)
    ax.set(xlabel="epsilon", ylabel="Delta_n")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=7)
    fig.tight_layout()
    path = output_dir / "baseline_Delta_prefixes.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    certified = []
    for value in epsilon:
        certified.append(max((int(row["prefix_n"]) for row in thresholds if isfinite(finite_float(row, "epsilon_certified_n")) and value <= finite_float(row, "epsilon_certified_n")), default=0))
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    ax.step(epsilon, certified, where="post", color="#8b4513", label="conservative N_certified,0")
    ax.scatter(epsilon, n_true, s=8, color="#1f4e79", alpha=0.45, label="sampled N_true")
    ax.axvspan(primary_min, primary_max, color="#dcefdc", alpha=0.35)
    ax.set(xlabel="epsilon", ylabel="prefix length", ylim=(-0.3, 10.5))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    path = output_dir / "baseline_certified_prefix_staircase.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    family_colors = {"axial": "#d62728", "bending_EB": "#1f77b4", "bending_Timoshenko": "#2ca02c", "axial_bending_cluster": "#9467bd"}
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for family, color in family_colors.items():
        selected = [row for row in families if str(row.get("EB_reference_family")) == family or str(row.get("Timo_family")) == family]
        ax.scatter([finite_float(row, "epsilon") for row in selected], [finite_float(row, "sorted_index") for row in selected], s=9, color=color, alpha=0.65, label=family)
    ax.set(xlabel="epsilon", ylabel="sorted index", ylim=(0.5, DEFAULT_N_SPECTRUM_ROOTS + 0.5))
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.tight_layout()
    path = output_dir / "baseline_sorted_family_map.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for index in range(1, DEFAULT_K_MAX + 1):
        rows = sorted((row for row in families if int(float(row.get("sorted_index", 0))) == index), key=lambda row: finite_float(row, "epsilon"))
        ax.plot([finite_float(row, "epsilon") for row in rows], [finite_float(row, "adjacent_gap_EB") for row in rows], linewidth=0.9, label=f"EB gap near {index}")
    ax.axhline(CLOSE_GAP_THRESHOLD, color="black", linestyle="--", linewidth=1.0)
    ax.set(xlabel="epsilon", ylabel="normalized adjacent gap", yscale="log")
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=6)
    fig.tight_layout()
    path = output_dir / "baseline_adjacent_gap_map.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)
    return paths


def markdown_value(value: object, digits: int = 8) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    return f"{number:.{digits}g}" if isfinite(number) else "--"


def write_report(
    args: Args,
    points: Sequence[Mapping[str, object]],
    modes: Sequence[Mapping[str, object]],
    thresholds: Sequence[Mapping[str, object]],
    transitions: Sequence[Mapping[str, object]],
    families: Sequence[Mapping[str, object]],
    mu_rows: Sequence[Mapping[str, object]],
    quality: Sequence[Mapping[str, object]],
    operations: Sequence[Mapping[str, object]],
    factorized_rows: Sequence[Mapping[str, object]],
    axial_rows: Sequence[Mapping[str, object]],
    general_comparison: Sequence[Mapping[str, object]],
    legacy_comparison: Sequence[Mapping[str, object]],
    regression_rows: Sequence[Mapping[str, object]],
) -> Path:
    report = args.output_dir / "eb_epsilon_baseline_thresholds_report.md"
    resolved_points = [row for row in points if str(row.get("root_quality_status")) == "resolved"]
    quality_counts = Counter(str(row.get("root_quality_status")) for row in quality)
    cache_counts = Counter(str(row.get("cache_status")) for row in quality)
    max_eb_abs = max((finite_float(row, "EB_reference_abs_error") for row in families), default=float("nan"))
    axial_eb_errors = [finite_float(row, "EB_reference_abs_error") for row in families if str(row.get("EB_reference_family")) == "axial"]
    axial_timo_errors = [
        finite_float(row, "Timo_axial_abs_error")
        for row in families
        if str(row.get("Timo_factorized_family_exact")) == "axial"
    ]
    mu_eb = max((finite_float(row, "EB_abs_difference") for row in mu_rows), default=float("nan"))
    mu_timo = max((finite_float(row, "Timo_abs_difference") for row in mu_rows), default=float("nan"))
    mu_failures = sum(str(row.get("status")) != "pass" for row in mu_rows)
    mu_n_true_failures = sum(not bool_value(row.get("N_true_equal")) for row in mu_rows)
    mu_eb_order_failures = sum(not bool_value(row.get("family_order_EB_equal")) for row in mu_rows)
    mu_timo_order_failures = sum(not bool_value(row.get("family_order_Timo_equal")) for row in mu_rows)
    mu_multiplicity_failures = sum(
        not bool_value(row.get("multiplicity_flags_equal"))
        or not bool_value(row.get("multiplicity_groups_Timo_equal"))
        for row in mu_rows
    )
    eb_warning_rows = sum(bool(str(row.get("EB_root_warning", "")).strip()) for row in quality)
    timo_warning_rows = sum(bool(str(row.get("Timo_root_warning", "")).strip()) for row in quality)
    candidate_boundary_rows = sum(str(row.get("candidate_boundary_warning")) == "true" for row in quality)
    multiplicity_point_rows = sum(int(float(row.get("possible_multiplicity_count", 0) or 0)) > 0 for row in quality)
    multiplicity_flags = sum(int(float(row.get("possible_multiplicity_count", 0) or 0)) for row in quality)
    reentries = [row for row in transitions if str(row.get("event_type")) in {"safe_reentry", "unsafe_reentry"}]
    reentry_counts = Counter(str(row.get("event_type")) for row in reentries)
    late_pass_points = [row for row in points if int(float(row.get("late_pass_count", 0) or 0)) > 0]
    factorized_failures = sum(str(row.get("root_quality_status")) != "pass" for row in factorized_rows)
    axial_failures = sum(str(row.get("status")) == "fail" for row in axial_rows)
    missing_general = [row for row in general_comparison if bool_value(row.get("missing_general_root"))]
    missing_general_confirmed = sum(
        bool_value(row.get("full_matrix_SVD_confirmation")) for row in missing_general
    )
    missing_general_by_model = Counter(str(row.get("model")) for row in missing_general)
    missing_required_axial = [
        row
        for row in missing_general
        if str(row.get("model")) == "Timoshenko"
        and str(row.get("factorized_family")) == "axial"
    ]
    full_matrix_verified = sum(
        str(row.get("root_quality_status")) == "pass" for row in factorized_rows
    )
    max_off_block = max(
        (
            finite_float(row, "off_block_max_abs")
            for row in factorized_rows
            if str(row.get("model")) == "Timoshenko"
        ),
        default=float("nan"),
    )
    legacy_counts = Counter(str(row.get("comparison_status")) for row in legacy_comparison)
    regression_failures = sum(str(row.get("status")) != "pass" for row in regression_rows)
    finite_bracket_widths = [
        finite_float(row, "bracket_width")
        for row in thresholds
        if isfinite(finite_float(row, "bracket_width"))
    ]
    regression_table_rows = [
        "| "
        f"{row.get('case_id')} | {markdown_value(row.get('epsilon_offset'))} | "
        f"{markdown_value(row.get('Lambda_bending_factorized'))} | "
        f"{markdown_value(row.get('Lambda_axial_exact'))} | "
        f"{markdown_value(row.get('pair_gap'))} | {row.get('status')} |"
        for row in regression_rows
    ]
    legacy_table_rows = [
        "| "
        f"{row.get('prefix_n')} | {row.get('legacy_threshold_status')} | "
        f"{markdown_value(row.get('legacy_epsilon_star_estimate'))} | "
        f"{row.get('corrected_threshold_status')} | "
        f"{markdown_value(row.get('corrected_epsilon_star_estimate'))} | "
        f"{row.get('correction_reason')} | {row.get('scientific_status')} |"
        for row in sorted(legacy_comparison, key=lambda item: int(item["prefix_n"]))
    ]
    lines = [
        "# EB Epsilon Baseline Thresholds Report",
        "",
        "## Correction Notice",
        "",
        "The earlier step-2 thresholds for prefixes 5--10 were produced by a general 6x6 determinant sign scan that can miss two close simple roots when both lie inside one Lambda step. The earlier independent verification repeated that same root-detection mechanism, so it could reproduce rather than expose the omission. Affected legacy thresholds are superseded.",
        "Legacy results and cache are preserved under `legacy_pre_factorized_root_fix/` and `cache_legacy_pre_factorized_root_fix/`; they are provenance only and do not contribute to this corrected certificate.",
        f"`spectrum_algorithm_version = {factorized.ALGORITHM_VERSION}`.",
        "Cache settings include algorithm version, epsilon, beta, mu, eta, requested/candidate root counts, and audit-local solver settings; an old-version entry is rejected as `stale_cache_algorithm_version`.",
        "",
        "## Scope",
        "",
        "- Baseline geometry: `beta=0 deg`, `mu=0`, `eta=0`.",
        f"- Sorted-spectrum target: `K={args.k_max}` with a `{100.0 * args.threshold:g}%` dimensional-frequency threshold.",
        f"- Primary engineering range: `{args.primary_epsilon_min:g}..{args.primary_epsilon_max:g}`.",
        f"- Computational buffer: `{args.epsilon_min:g}..{args.epsilon_max:g}`.",
        "- Timoshenko is a one-dimensional reference, not an exact 3D solution.",
        "- No FEM rerun, formula/determinant/shared-root-solver modification, or beta/mu/eta adversarial search was performed; research step 3 was not started.",
        "",
        "## Exact Straight-System Factorization",
        "",
        "At `beta=0`, `eta=0`, the unchanged Timoshenko 6x6 matrix with unknown ordering `(A1,B1,P1,A2,B2,P2)` is exactly block diagonal after selecting bending rows/columns `[0,2,3,4]/[0,1,3,4]` and axial rows/columns `[1,5]/[2,5]`. This is a matrix factorization, not a new physical model or a new characteristic equation.",
        "Changing `mu` moves only the artificial internal joint inside one straight homogeneous fixed-fixed rod; it does not change the physical rod or its spectrum.",
        "The independent axial reference follows `Lambda_axial,m=sqrt(m*pi/(2*epsilon))` from `theta=epsilon*Lambda^2` over total length 2.",
        "The independent EB bending reference uses `cosh(2*Lambda) cos(2*Lambda)=1`; `fixed_fixed_eb_roots` supplies the existing alpha convention and the audit uses `Lambda=alpha/2`.",
        "Within-family numerical duplicates are merged; axial/bending records are never merged, so exact and near cross-family crossings retain multiplicity two with deterministic stable ordering.",
        f"- Maximum off-block matrix norm: `{markdown_value(max_off_block)}`.",
        f"- Maximum corrected analytic EB union absolute discrepancy: `{markdown_value(max_eb_abs)}`.",
        f"- Maximum EB axial-family absolute discrepancy: `{markdown_value(max(axial_eb_errors, default=float('nan')))}`.",
        f"- Maximum Timoshenko axial-family absolute discrepancy: `{markdown_value(max(axial_timo_errors, default=float('nan')))}`.",
        f"- Mu-invariance maximum EB/Timoshenko absolute drift: `{markdown_value(mu_eb)}` / `{markdown_value(mu_timo)}`.",
        f"- Corrected factorized mu-invariance failing rows: `{mu_failures}`; `N_true` mismatches: `{mu_n_true_failures}`; EB/Timoshenko family-order mismatches: `{mu_eb_order_failures}/{mu_timo_order_failures}`; multiplicity mismatches: `{mu_multiplicity_failures}`.",
        "",
        "## Root Completeness",
        "",
        f"- Total factorized root records evaluated: `{len(factorized_rows)}`; unresolved: `{factorized_failures}`.",
        f"- Full-matrix SVD verified factorized records: `{full_matrix_verified}` / `{len(factorized_rows)}`.",
        f"- Exact axial roots required in first 12: `{len(axial_rows)}`; block/full-matrix failures: `{axial_failures}`.",
        f"- Raw general-6x6 missing factorized roots: `{len(missing_general)}` (EB `{missing_general_by_model['Euler-Bernoulli']}`, Timoshenko `{missing_general_by_model['Timoshenko']}`); missing required axial records: `{len(missing_required_axial)}`.",
        f"- Missing roots locally confirmed by full-matrix SVD: `{missing_general_confirmed}`; unresolved confirmations: `{len(missing_general) - missing_general_confirmed}`.",
        f"- R1--R3 and +/-1e-6 regression failures: `{regression_failures}` / `{len(regression_rows)}` rows.",
        "",
        "| case | epsilon offset | bending Lambda | axial Lambda | gap | status |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
        *regression_table_rows,
        "",
        "## Mu-Invariance",
        "",
        "The physical factorized straight-rod spectrum is independent of the artificial internal split `mu`; the raw general solver is retained only as a diagnostic comparison.",
        f"- Maximum corrected EB/Timoshenko absolute drift: `{markdown_value(mu_eb)}` / `{markdown_value(mu_timo)}`.",
        f"- Failing rows: `{mu_failures}` / `{len(mu_rows)}`; `N_true` mismatches: `{mu_n_true_failures}`; EB/Timoshenko family-order mismatches: `{mu_eb_order_failures}/{mu_timo_order_failures}`; multiplicity mismatches: `{mu_multiplicity_failures}`.",
        f"- The raw general comparison reports `{len(missing_general)}` missing records; those diagnostic omissions are not interpreted as physical mu-dependence.",
        "",
        "## Coarse and Adaptive Scan",
        "",
        f"- Unique evaluated baseline epsilon points: `{len(points)}` (`{len(resolved_points)}` resolved).",
        f"- Root-quality resolved/unresolved rows across all scopes: `{quality_counts['resolved']}/{quality_counts['unresolved']}`.",
        f"- Cache hit/miss/mixed/force-recomputed rows: `{cache_counts['hit']}/{cache_counts['miss']}/{cache_counts['mixed']}/{cache_counts['force_recomputed']}`.",
        f"- EB/Timoshenko root-warning rows: `{eb_warning_rows}/{timo_warning_rows}`; candidate-boundary rows: `{candidate_boundary_rows}`.",
        f"- Possible-multiplicity point rows/first-12 flags: `{multiplicity_point_rows}/{multiplicity_flags}`.",
        f"- Recorded transition/event rows: `{len(transitions)}`.",
        "- Every coarse interval received a deterministic midpoint screening pass; event intervals received local deterministic refinement.",
        "",
        "## Corrected Thresholds",
        "",
        "| n | status | safe lower | unsafe upper | epsilon_star_estimate (display 5dp) | floor 5dp / 4dp | trigger mode/family | continuity | verification |",
        "| ---: | --- | ---: | ---: | ---: | --- | --- | --- | --- |",
    ]
    for row in sorted(thresholds, key=lambda item: int(item["prefix_n"])):
        lines.append(
            f"| {row['prefix_n']} | {row.get('threshold_status')} | {markdown_value(row.get('epsilon_safe_lower'))} | {markdown_value(row.get('epsilon_unsafe_upper'))} | {markdown_value(row.get('epsilon_star_estimate'))} ({row.get('epsilon_star_display_5dp') or '--'}) | {row.get('epsilon_certified_floor_5dp') or '--'} / {row.get('epsilon_certified_floor_4dp') or '--'} | "
            f"{row.get('triggering_sorted_indices') or '--'} / {row.get('triggering_EB_family') or '--'} -> {row.get('triggering_Timo_family') or '--'} | {row.get('threshold_continuity_status')} | "
            f"{row.get('numerical_verification_status')} |"
        )
    lines.extend(
        [
            "",
            "## Baseline Certificate",
            "",
            "The conservative straight-system numerical certificate is",
            "",
            "```text",
            "N_certified_0(epsilon) = max { n : epsilon <= epsilon_certified_n }.",
            "```",
            "",
            "For a resolved first loss, `epsilon_certified_n` is the verified safe lower endpoint, not the midpoint estimate. For `not_reached_in_scan_range`, it is the right-censored verified endpoint of the computational buffer; no first-loss estimate is claimed. This certificate applies only to the straight baseline system and only over the scanned epsilon range.",
            "",
            "## Precision Policy",
            "",
            f"Internal first-loss brackets are refined to approximately `{args.epsilon_tolerance:g}` or tighter. Publication uses five decimal places; practical certificate fields use downward decimal flooring at five and four places. Presentation fields never replace full-precision values, and any future step-3 points use the full-precision verified bracket/estimate.",
            f"- Maximum resolved bracket width: `{markdown_value(max(finite_bracket_widths, default=float('nan')))}`.",
            "",
            "## Legacy Comparison",
            "",
            f"- Unchanged within tolerance: `{legacy_counts['unchanged_within_tolerance']}`.",
            f"- Corrected due to a missing root: `{legacy_counts['corrected_due_to_missing_root']}`.",
            f"- Previously unresolved now resolved / still unresolved / not reached: `{legacy_counts['previously_unresolved_now_resolved']}` / `{legacy_counts['still_unresolved']}` / `{legacy_counts['not_reached']}`.",
            "",
            "| n | legacy status | legacy estimate | corrected status | corrected estimate | reason | scientific status |",
            "| ---: | --- | ---: | --- | ---: | --- | --- |",
            *legacy_table_rows,
            "",
            "## Re-Entry",
            "",
            f"- Prefix safe/unsafe re-entry events: `{reentry_counts['safe_reentry']}` / `{reentry_counts['unsafe_reentry']}` (`{len(reentries)}` total).",
            f"- Evaluated points containing late individual passes: `{len(late_pass_points)}`.",
            f"- Family reorder events: `{sum(str(row.get('event_type')) == 'family_reorder' for row in transitions)}`.",
            "- First loss defines every threshold; a later safe interval or individual late pass never increases the certificate.",
            "",
            "## Numerical Verification",
            "",
            f"- Thresholds passing independent force-recompute verification: `{sum(str(row.get('numerical_verification_status')) == 'pass' for row in thresholds)}` / `{sum(str(row.get('threshold_status')) not in {'below_scan_range', 'not_reached_in_scan_range'} for row in thresholds)}` applicable rows.",
            f"- Independent verification candidate margin: `{max(VERIFICATION_CANDIDATE_ROOTS, args.n_candidate_roots + 4)}` for the raw audit; the corrected target force-recomputed `{factorized.ALGORITHM_VERSION}` without primary cache reuse.",
            "- Verification repeated the first-12 multiplicity-preserving family order, exact axial block checks, and full-6x6 SVD confirmation.",
            "",
            "## Operation Counts",
            "",
            "| metric | value | status |",
            "| --- | ---: | --- |",
        ]
    )
    for row in operations:
        lines.append(f"| {row.get('metric')} | {row.get('value')} | {row.get('measurement_status')} |")
    lines.extend(
        [
            "",
            "No online epsilon-to-Pi cascade cost is claimed by this baseline audit.",
            "",
            "## Limitations",
            "",
            "- This is a straight-system baseline only; it is not a proof of a lower-envelope property over `beta`, `mu`, or `eta`.",
            "- No adversarial geometry search has been run yet.",
            "- Timoshenko remains a one-dimensional reference.",
            "- Existing 3D support is strongest in the primary range `0.01..0.05` and does not automatically extend over the whole computational buffer.",
            "- Exact crossings and multiplicity require special numerical care; cross-family duplicates are retained and only within-family numerical duplicates are removed.",
            "- The raw general 6x6 sign scan remains a useful independent completeness diagnostic, but its missing-root rows are not a production source for this special-limit certificate.",
            "- The general production solver outside the straight limit was not changed; future general-geometry work still requires a dedicated root-completeness recovery audit.",
            "- Thresholds outside the primary range have weaker physical validation status.",
            "",
            "## Boundary of This Correction",
            "",
            "| n | epsilon_star_n | epsilon_near_n | epsilon_buffer_n |",
            "| ---: | ---: | ---: | ---: |",
        ]
    )
    for row in sorted(thresholds, key=lambda item: int(item["prefix_n"])):
        lines.append(
            f"| {row['prefix_n']} | {markdown_value(row.get('epsilon_star_estimate'))} | {markdown_value(row.get('epsilon_near_n'))} | {markdown_value(row.get('epsilon_buffer_n'))} |"
        )
    lines.extend(
        [
            "",
            "These values are retained as possible inputs to a future targeted counterexample search only. No geometry search was implemented or run, and this task does not transition to research step 3.",
        ]
    )
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    del modes
    return report


def run_plot_only(args: Args) -> dict[str, object]:
    points = read_csv(args.output_dir / "baseline_epsilon_point_summary.csv")
    modes = read_csv(args.output_dir / "baseline_mode_metrics.csv")
    thresholds = read_csv(args.output_dir / "baseline_critical_prefix_thresholds.csv")
    transitions = read_csv(args.output_dir / "baseline_transition_audit.csv")
    families = read_csv(args.output_dir / "baseline_spectral_family_audit.csv")
    mu_rows = read_csv(args.output_dir / "baseline_mu_invariance_audit.csv")
    quality = read_csv(args.output_dir / "baseline_solver_quality_audit.csv")
    operations = read_csv(args.output_dir / "baseline_operation_counts.csv")
    factorized_rows = read_csv(args.output_dir / "baseline_factorized_spectrum_audit.csv")
    axial_rows = read_csv(args.output_dir / "baseline_axial_completeness_audit.csv")
    general_comparison = read_csv(args.output_dir / "baseline_general_vs_factorized_spectrum_audit.csv")
    legacy_comparison = read_csv(args.output_dir / "baseline_legacy_threshold_comparison.csv")
    regression_rows = read_csv(args.output_dir / "baseline_regression_close_root_pairs.csv")
    plots = make_plots(
        args.output_dir,
        points,
        modes,
        thresholds,
        families,
        threshold=args.threshold,
        primary_min=args.primary_epsilon_min,
        primary_max=args.primary_epsilon_max,
    )
    report = write_report(
        args,
        points,
        modes,
        thresholds,
        transitions,
        families,
        mu_rows,
        quality,
        operations,
        factorized_rows,
        axial_rows,
        general_comparison,
        legacy_comparison,
        regression_rows,
    )
    print(f"plot-only regenerated {len(plots)} plots; root calculations performed: 0")
    return {"args": args, "plot_paths": plots, "report": report, "root_solves": 0}


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    if args.plot_only:
        return run_plot_only(args)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    primary = EvaluationManager(args)
    coarse_values = epsilon_grid(args)
    for epsilon in coarse_values:
        primary.evaluate(epsilon, "coarse")
    regression_rows = regression_close_root_rows(primary)
    adaptive_rows, adaptive_intervals = adaptive_screen(args, primary, coarse_values)
    baseline_before_thresholds = ordered_baseline_evaluations(primary, args)
    threshold_rows, threshold_refinement_rows = build_threshold_rows(
        args, primary, baseline_before_thresholds
    )
    verification, verification_rows = verify_thresholds(args, threshold_rows)
    baseline_evaluations = ordered_baseline_evaluations(primary, args)
    mu_rows = mu_invariance_rows(args, primary, threshold_rows)
    transitions = [*adaptive_intervals, *transition_rows(baseline_evaluations, args)]
    points = point_rows(baseline_evaluations, args)
    modes = mode_rows(baseline_evaluations, args)
    families = family_rows(baseline_evaluations, args)
    refinement = [*adaptive_rows, *threshold_refinement_rows, *verification_rows]
    quality = quality_rows((primary, verification), args)
    primary_all_evaluations = list(primary.evaluations.values())
    factorized_rows = factorized_spectrum_rows(primary_all_evaluations)
    axial_rows = axial_completeness_rows(primary_all_evaluations, args.n_spectrum_roots)
    general_comparison = general_vs_factorized_rows((primary, verification))
    legacy_comparison = legacy_threshold_comparison_rows(
        args, threshold_rows, general_comparison
    )
    elapsed = time.perf_counter() - started
    operations = operation_rows(
        args,
        primary,
        verification,
        baseline_evaluations,
        refinement,
        elapsed,
    )

    write_csv(args.output_dir / "baseline_epsilon_point_summary.csv", points, POINT_FIELDS)
    write_csv(args.output_dir / "baseline_mode_metrics.csv", modes, MODE_FIELDS)
    write_csv(args.output_dir / "baseline_critical_prefix_thresholds.csv", threshold_rows, THRESHOLD_FIELDS)
    write_csv(args.output_dir / "baseline_threshold_refinement_audit.csv", refinement, REFINEMENT_FIELDS)
    write_csv(args.output_dir / "baseline_transition_audit.csv", transitions, TRANSITION_FIELDS)
    write_csv(args.output_dir / "baseline_spectral_family_audit.csv", families, FAMILY_FIELDS)
    write_csv(args.output_dir / "baseline_mu_invariance_audit.csv", mu_rows, MU_FIELDS)
    write_csv(args.output_dir / "baseline_solver_quality_audit.csv", quality, QUALITY_FIELDS)
    write_csv(args.output_dir / "baseline_operation_counts.csv", operations, OPERATION_FIELDS)
    write_csv(args.output_dir / "baseline_factorized_spectrum_audit.csv", factorized_rows, FACTORIZED_FIELDS)
    write_csv(args.output_dir / "baseline_axial_completeness_audit.csv", axial_rows, AXIAL_COMPLETENESS_FIELDS)
    write_csv(args.output_dir / "baseline_general_vs_factorized_spectrum_audit.csv", general_comparison, GENERAL_VS_FACTORIZED_FIELDS)
    write_csv(args.output_dir / "baseline_legacy_threshold_comparison.csv", legacy_comparison, LEGACY_COMPARISON_FIELDS)
    write_csv(args.output_dir / "baseline_regression_close_root_pairs.csv", regression_rows, REGRESSION_FIELDS)
    plots = make_plots(
        args.output_dir,
        points,
        modes,
        threshold_rows,
        families,
        threshold=args.threshold,
        primary_min=args.primary_epsilon_min,
        primary_max=args.primary_epsilon_max,
    )
    report = write_report(
        args,
        points,
        modes,
        threshold_rows,
        transitions,
        families,
        mu_rows,
        quality,
        operations,
        factorized_rows,
        axial_rows,
        general_comparison,
        legacy_comparison,
        regression_rows,
    )
    unresolved = sum(str(row.get("root_quality_status")) != "resolved" for row in quality)
    print(f"baseline epsilon points: {len(points)}; unresolved quality rows: {unresolved}")
    print(f"root cache hits/misses: {primary.runtime.root_cache_hits}/{primary.runtime.root_cache_misses}")
    print(f"output: {args.output_dir}")
    return {
        "args": args,
        "primary_manager": primary,
        "verification_manager": verification,
        "point_rows": points,
        "mode_rows": modes,
        "threshold_rows": threshold_rows,
        "refinement_rows": refinement,
        "transition_rows": transitions,
        "family_rows": families,
        "mu_rows": mu_rows,
        "quality_rows": quality,
        "operation_rows": operations,
        "factorized_rows": factorized_rows,
        "axial_rows": axial_rows,
        "general_comparison_rows": general_comparison,
        "legacy_comparison_rows": legacy_comparison,
        "regression_rows": regression_rows,
        "plot_paths": plots,
        "report": report,
        "elapsed_seconds": elapsed,
    }


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from None
