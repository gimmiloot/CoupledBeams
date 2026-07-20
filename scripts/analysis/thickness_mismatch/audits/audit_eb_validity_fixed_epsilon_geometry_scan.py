from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import math
from math import isfinite
from pathlib import Path
import sys
import time
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq, linear_sum_assignment, minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    thickness_mismatch_factors,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_timoshenko_shape_construction as shape_audit,
)
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow,
)
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


CACHE_VERSION = "eb_validity_fixed_epsilon_geometry_scan_v1"
DEFAULT_OUTPUT_DIR = Path("results") / "eb_validity_fixed_epsilon_geometry_scan"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_validity_fixed_epsilon_geometry_scan"
DEFAULT_CACHE_SUBDIR = "cache"

DEFAULT_EPSILON = 0.02
DEFAULT_BETA_VALUES = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
DEFAULT_MU_VALUES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
DEFAULT_ETA_VALUES = (-0.1, 0.0, 0.1)
SMOKE_BETA_VALUES = (0.0, 45.0, 90.0)
SMOKE_MU_VALUES = (0.0, 0.7)
SMOKE_ETA_VALUES = (-0.1, 0.1)
DEFAULT_N_REPORTED_MODES = 8
DEFAULT_N_CANDIDATE_ROOTS = 16
BOUNDARY_RETRY_ROOTS = 20
SMOKE_N_REPORTED_MODES = 4
SMOKE_N_CANDIDATE_ROOTS = 8
DEFAULT_SHAPE_POINTS = 401
SMOKE_SHAPE_POINTS = 201
DEFAULT_CLUSTER_GAP_THRESHOLD = 0.02

MODEL_EB = beta_workflow.MODEL_EB
MODEL_TIMO = beta_workflow.MODEL_TIMO
MAIN_THRESHOLD = 0.10
MAC_MIN_RELIABLE = 0.80
MAC_MARGIN_MIN_RELIABLE = 0.05
NUMERICAL_FLOOR = 1.0e-12

PREDICTORS = (
    "epsilon_max",
    "chi_max_EB",
    "chi_eff_EB",
    "Theta_max_EB",
    "Pi_EB",
    "raw_parameter_baseline",
)
SCALAR_PREDICTORS = (
    "epsilon_max",
    "chi_max_EB",
    "chi_eff_EB",
    "Theta_max_EB",
    "Pi_EB",
)
SUBGROUPS = ("all_reliable", "bending", "mixed", "longitudinal", "cluster_level")
SPLITS = ("leave_one_beta_out", "leave_one_mu_out", "leave_one_eta_out")
THRESHOLD_KINDS = (
    ("central", None),
    ("safe_0pct", 0.0),
    ("safe_1pct", 0.01),
    ("safe_5pct", 0.05),
)

MODE_FIELDS = [
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "tau_1",
    "tau_2",
    "epsilon_1",
    "epsilon_2",
    "epsilon_max",
    "l1_over_r1",
    "l2_over_r2",
    "l1_over_d1",
    "l2_over_d2",
    "comparison_type",
    "eb_sorted_index",
    "timo_sorted_index",
    "cluster_id",
    "cluster_size_EB",
    "cluster_size_Timo",
    "subspace_MAC",
    "cluster_centroid_delta_f",
    "cluster_max_delta_f",
    "Pi_cluster_max",
    "Theta_cluster_max",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_Lambda",
    "delta_f",
    "delta_f_symmetric",
    "bias_f",
    "MAC_u",
    "MAC_w",
    "MAC_uw",
    "MAC_psi_vs_EB_slope",
    "MAC_margin",
    "matching_status",
    "include_in_calibration",
    "EB_axial_fraction",
    "EB_bending_fraction",
    "EB_classification",
    "Timo_axial_fraction",
    "Timo_bending_fraction",
    "Timo_shear_fraction",
    "Timo_bending_shear_fraction",
    "Timo_classification",
    "chi_1_EB",
    "chi_2_EB",
    "chi_max_EB",
    "chi_eff_EB",
    "Theta_max_EB",
    "Pi_shear_EB",
    "Pi_rotary_EB",
    "Pi_EB",
    "delta_f_from_Pi_half",
    "residual_Pi_half",
    "pass_10pct",
    "root_warning",
    "candidate_boundary",
    "notes",
]

MATCHING_FIELDS = [
    "beta_deg",
    "mu",
    "eta",
    "eb_sorted_index",
    "timo_candidate_index",
    "Lambda_EB",
    "Lambda_Timo",
    "MAC_u",
    "MAC_w",
    "MAC_uw",
    "MAC_margin",
    "gap_EB",
    "gap_Timo",
    "cluster_EB",
    "cluster_Timo",
    "subspace_MAC",
    "assignment_selected",
    "matching_status",
    "warning",
    "notes",
]

POINT_SUMMARY_FIELDS = [
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "target_mode_count",
    "N10_sorted_consecutive",
    "pass_pattern_10",
    "n_pass_10_among_first8",
    "n_pass_10_among_firstK",
    "N10_current_bending",
    "n_bending_pass_10",
    "reliable_individual_count",
    "reliable_cluster_count",
    "ambiguous_count",
    "root_warning_count",
    "candidate_boundary_count",
    "cluster_count_EB",
    "cluster_count_Timo",
    "matched_cluster_count",
    "max_delta_f_first8",
    "mean_delta_f_first8",
    "max_delta_f_firstK",
    "mean_delta_f_firstK",
    "notes",
]

FIT_FIELDS = [
    "predictor",
    "subgroup",
    "sample_filter",
    "point_count",
    "C",
    "p",
    "R2",
    "RMSE",
    "MAE",
    "threshold_10_central",
    "threshold_safe_0pct",
    "threshold_safe_1pct",
    "threshold_safe_5pct",
    "false_safe_rate_central",
    "false_unsafe_rate_central",
    "balanced_accuracy_central",
    "safe_coverage_central",
    "geometry_extended_R2",
    "geometry_extended_RMSE",
    "geometry_extended_MAE",
    "notes",
]

CV_FIELDS = [
    "predictor",
    "subgroup",
    "split_kind",
    "held_out_value",
    "train_count",
    "test_count",
    "C",
    "p",
    "threshold_10",
    "test_RMSE",
    "test_MAE",
    "false_safe_rate",
    "false_unsafe_rate",
    "balanced_accuracy",
    "safe_coverage",
    "notes",
]

GEOMETRY_FIELDS = [
    "predictor",
    "grouping_variable",
    "group_value",
    "row_count",
    "mean_residual",
    "median_residual",
    "max_abs_residual",
    "p95_abs_residual",
    "false_safe_rate",
    "false_unsafe_rate",
    "notes",
]

THRESHOLD_CLASS_FIELDS = [
    "predictor",
    "threshold_kind",
    "threshold_value",
    "subgroup",
    "beta_deg",
    "mu",
    "eta",
    "true_safe",
    "predicted_safe",
    "false_safe",
    "false_unsafe",
    "notes",
]

RUNTIME_FIELDS = [
    "stage",
    "parameter_point_count",
    "mode_count",
    "total_seconds",
    "mean_seconds_per_point",
    "mean_seconds_per_mode",
    "matrix_evaluation_count",
    "relative_cost_vs_full_reference",
    "notes",
]

LOCAL_TIMO_FIELDS = [
    "beta_deg",
    "mu",
    "eta",
    "mode",
    "Lambda_EB",
    "Lambda_Timo_full",
    "Lambda_Timo_local",
    "delta_f_full",
    "delta_f_local",
    "local_timo_matrix_evaluations",
    "full_timo_matrix_evaluations",
    "local_runtime",
    "full_runtime",
    "speedup",
    "local_delta_error",
    "local_root_failure",
    "cluster_warning",
    "notes",
]


@dataclass(frozen=True)
class Args:
    epsilon: float
    beta_values: tuple[float, ...]
    mu_values: tuple[float, ...]
    eta_values: tuple[float, ...]
    n_reported_modes: int
    n_candidate_roots: int
    n_shape_points: int
    cluster_gap_threshold: float
    benchmark_local_timo: bool
    workers: int
    output_dir: Path
    cache_dir: Path
    reuse_cache: bool
    force_recompute: bool
    plot_only: bool
    smoke: bool


@dataclass(frozen=True)
class RootCacheEntry:
    roots: tuple[float, ...]
    warnings: tuple[str, ...]
    root_count_found: int
    lambda_max_used: float
    scan_step_used: float
    retry_attempted: bool
    retry_changed_value: bool
    notes: tuple[str, ...]
    cache_hit: bool


@dataclass(frozen=True)
class RodFields:
    x: np.ndarray
    u: np.ndarray
    w: np.ndarray
    u_prime: np.ndarray
    w_prime: np.ndarray
    w_second: np.ndarray | None
    w_third: np.ndarray | None
    psi: np.ndarray | None
    psi_prime: np.ndarray | None
    gamma: np.ndarray | None


@dataclass(frozen=True)
class ModeResult:
    model: str
    epsilon: float
    beta_deg: float
    mu: float
    eta: float
    sorted_index: int
    Lambda: float
    coeff: np.ndarray
    rod1: RodFields
    rod2: RodFields
    energy: dict[str, float | str]
    warnings: tuple[str, ...]


@dataclass
class RuntimeTracker:
    eb_roots_seconds: float = 0.0
    theta_seconds: float = 0.0
    eb_shape_pi_seconds: float = 0.0
    timo_roots_seconds: float = 0.0
    matching_seconds: float = 0.0
    local_timo_seconds: float = 0.0
    eb_root_calls: int = 0
    timo_root_calls: int = 0
    root_cache_hits: int = 0
    root_cache_misses: int = 0
    boundary_retries: int = 0


@dataclass(frozen=True)
class PointProducts:
    mode_rows: list[dict[str, object]]
    matching_rows: list[dict[str, object]]
    summary_row: dict[str, object]
    local_timo_rows: list[dict[str, object]]


def repo_path(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def rel(path: Path) -> str:
    try:
        return str(Path(path).resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if not isfinite(value_f):
            return "nan"
        return f"{value_f:.16e}"
    if isinstance(value, (bool, np.bool_)):
        return "true" if bool(value) else "false"
    return value


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def finite_float(row: dict[str, object] | dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def trapz(values: np.ndarray, coordinates: np.ndarray) -> float:
    return shape_audit.trapz(np.asarray(values, dtype=float), np.asarray(coordinates, dtype=float))


def parse_float_list(values: Sequence[str] | None, default: Sequence[float]) -> tuple[float, ...]:
    if values is None:
        return tuple(float(value) for value in default)
    return tuple(float(value) for value in values)


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Fixed-epsilon diagnostic EB validity scan relative to Timoshenko.",
    )
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--beta-values", type=str, nargs="+", default=None)
    parser.add_argument("--mu-values", type=str, nargs="+", default=None)
    parser.add_argument("--eta-values", type=str, nargs="+", default=None)
    parser.add_argument("--n-reported-modes", type=int, default=DEFAULT_N_REPORTED_MODES)
    parser.add_argument("--n-candidate-roots", type=int, default=DEFAULT_N_CANDIDATE_ROOTS)
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_SHAPE_POINTS)
    parser.add_argument("--cluster-gap-threshold", type=float, default=DEFAULT_CLUSTER_GAP_THRESHOLD)
    parser.add_argument("--benchmark-local-timo", action="store_true")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    beta_values = parse_float_list(ns.beta_values, DEFAULT_BETA_VALUES)
    mu_values = parse_float_list(ns.mu_values, DEFAULT_MU_VALUES)
    eta_values = parse_float_list(ns.eta_values, DEFAULT_ETA_VALUES)
    output_dir = repo_path(Path(ns.output_dir))
    n_reported = int(ns.n_reported_modes)
    n_candidate = int(ns.n_candidate_roots)
    n_shape = int(ns.n_shape_points)
    if ns.smoke:
        beta_values = SMOKE_BETA_VALUES if ns.beta_values is None else beta_values
        mu_values = SMOKE_MU_VALUES if ns.mu_values is None else mu_values
        eta_values = SMOKE_ETA_VALUES if ns.eta_values is None else eta_values
        n_reported = min(n_reported, SMOKE_N_REPORTED_MODES)
        n_candidate = min(n_candidate, SMOKE_N_CANDIDATE_ROOTS)
        n_shape = min(n_shape, SMOKE_SHAPE_POINTS)
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            output_dir = repo_path(SMOKE_OUTPUT_DIR)
        ns.benchmark_local_timo = False
    cache_dir = repo_path(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / DEFAULT_CACHE_SUBDIR
    args = Args(
        epsilon=float(ns.epsilon),
        beta_values=tuple(float(value) for value in beta_values),
        mu_values=tuple(float(value) for value in mu_values),
        eta_values=tuple(float(value) for value in eta_values),
        n_reported_modes=n_reported,
        n_candidate_roots=n_candidate,
        n_shape_points=n_shape,
        cluster_gap_threshold=float(ns.cluster_gap_threshold),
        benchmark_local_timo=bool(ns.benchmark_local_timo),
        workers=max(1, int(ns.workers)),
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
    if args.epsilon <= 0.0:
        raise ValueError("--epsilon must be positive")
    if args.n_reported_modes < 1:
        raise ValueError("--n-reported-modes must be positive")
    if args.n_candidate_roots < args.n_reported_modes:
        raise ValueError("--n-candidate-roots must be at least --n-reported-modes")
    if args.n_shape_points < 51:
        raise ValueError("--n-shape-points must be at least 51")
    if args.cluster_gap_threshold <= 0.0:
        raise ValueError("--cluster-gap-threshold must be positive")
    if not args.beta_values or not args.mu_values or not args.eta_values:
        raise ValueError("beta, mu, and eta grids must be nonempty")
    for beta in args.beta_values:
        if not (0.0 <= float(beta) <= 90.0):
            raise ValueError("beta values must lie in [0, 90] degrees")
    for mu in args.mu_values:
        for eta in args.eta_values:
            thickness_mismatch_factors(float(mu), float(eta))


def local_geometry(epsilon: float, mu: float, eta: float) -> dict[str, float]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    eps1 = float(epsilon) * factors.tau1 / (1.0 - float(mu))
    eps2 = float(epsilon) * factors.tau2 / (1.0 + float(mu))
    return {
        "tau_1": float(factors.tau1),
        "tau_2": float(factors.tau2),
        "epsilon_1": float(eps1),
        "epsilon_2": float(eps2),
        "epsilon_max": max(float(eps1), float(eps2)),
        "l1_over_r1": 1.0 / (2.0 * float(eps1)),
        "l2_over_r2": 1.0 / (2.0 * float(eps2)),
        "l1_over_d1": 1.0 / (4.0 * float(eps1)),
        "l2_over_d2": 1.0 / (4.0 * float(eps2)),
    }


def frequency_metrics(lambda_eb: float, lambda_timo: float) -> dict[str, float]:
    eb = float(lambda_eb)
    timo = float(lambda_timo)
    if not (isfinite(eb) and isfinite(timo)) or abs(timo) <= 1.0e-30:
        return {
            "delta_Lambda": float("nan"),
            "delta_f": float("nan"),
            "delta_f_symmetric": float("nan"),
            "bias_f": float("nan"),
        }
    eb_sq = eb * eb
    timo_sq = timo * timo
    denom_sym = eb_sq + timo_sq
    return {
        "delta_Lambda": abs(eb - timo) / abs(timo),
        "delta_f": abs(eb_sq - timo_sq) / timo_sq,
        "delta_f_symmetric": 2.0 * abs(eb_sq - timo_sq) / denom_sym if denom_sym > 0.0 else float("nan"),
        "bias_f": (eb_sq - timo_sq) / timo_sq,
    }


def pass_bool(delta_f: float, threshold: float = MAIN_THRESHOLD) -> bool:
    return isfinite(float(delta_f)) and float(delta_f) <= float(threshold)


def consecutive_pass_count(deltas: Sequence[float], threshold: float = MAIN_THRESHOLD) -> int:
    count = 0
    for value in deltas:
        if pass_bool(float(value), threshold):
            count += 1
        else:
            break
    return count


def pass_pattern(deltas: Sequence[float], threshold: float = MAIN_THRESHOLD) -> str:
    return ",".join("Y" if pass_bool(float(value), threshold) else "N" for value in deltas)


def count_passes(deltas: Sequence[float], threshold: float = MAIN_THRESHOLD) -> int:
    return sum(1 for value in deltas if pass_bool(float(value), threshold))


def k_aware_point_summary_fields(
    deltas: Sequence[float],
    *,
    target_mode_count: int,
) -> dict[str, object]:
    target = [float(value) for value in deltas[: int(target_mode_count)]]
    legacy_first8 = target[:8]
    return {
        "target_mode_count": int(target_mode_count),
        "N10_sorted_consecutive": consecutive_pass_count(target),
        "pass_pattern_10": pass_pattern(target),
        "n_pass_10_among_first8": count_passes(legacy_first8),
        "n_pass_10_among_firstK": count_passes(target),
        "max_delta_f_first8": max(legacy_first8) if legacy_first8 else float("nan"),
        "mean_delta_f_first8": float(np.mean(legacy_first8)) if legacy_first8 else float("nan"),
        "max_delta_f_firstK": max(target) if target else float("nan"),
        "mean_delta_f_firstK": float(np.mean(target)) if target else float("nan"),
    }


def classify_energy(model: str, axial_fraction: float, bending_fraction: float, shear_fraction: float) -> str:
    if isfinite(float(axial_fraction)) and float(axial_fraction) >= 0.70:
        return "longitudinal_dominated"
    if model == MODEL_TIMO:
        bending_shear = float(bending_fraction) + float(shear_fraction)
        if isfinite(bending_shear) and bending_shear >= 0.70:
            return "bending_dominated"
    elif isfinite(float(bending_fraction)) and float(bending_fraction) >= 0.70:
        return "bending_dominated"
    return "mixed"


def energy_from_terms(model: str, U_axial: float, U_bending: float, U_shear: float) -> dict[str, float | str]:
    total = float(U_axial) + float(U_bending) + float(U_shear)
    axial = safe_ratio(float(U_axial), total)
    bending = safe_ratio(float(U_bending), total)
    shear = safe_ratio(float(U_shear), total)
    return {
        "U_axial": max(float(U_axial), 0.0),
        "U_bending": max(float(U_bending), 0.0),
        "U_shear": max(float(U_shear), 0.0),
        "U_total": max(float(total), 0.0),
        "axial_fraction": axial,
        "bending_fraction": bending,
        "shear_fraction": shear,
        "classification": classify_energy(model, axial, bending, shear),
    }


def theta_factor() -> float:
    kappa = TIMO.circular_shear_coefficient(TIMO.NU)
    shear_modulus = TIMO.E / (2.0 * (1.0 + TIMO.NU))
    return 1.0 + TIMO.E / (kappa * shear_modulus)


def chi_values_eb(lambda_eb: float, epsilon: float, mu: float, eta: float) -> dict[str, float]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    chi1 = float(epsilon) * float(lambda_eb) * math.sqrt(factors.tau1)
    chi2 = float(epsilon) * float(lambda_eb) * math.sqrt(factors.tau2)
    return {"chi_1_EB": chi1, "chi_2_EB": chi2, "chi_max_EB": max(chi1, chi2)}


def theta_max_eb(lambda_eb: float, epsilon: float, mu: float, eta: float) -> float:
    chis = chi_values_eb(float(lambda_eb), float(epsilon), float(mu), float(eta))
    return theta_factor() * chis["chi_max_EB"] ** 2


def chi_eff_from_eb_bending_energies(chi1: float, chi2: float, bending1: float, bending2: float) -> float:
    denom = float(bending1) + float(bending2)
    if not isfinite(denom) or denom <= 1.0e-30:
        return float("nan")
    value = (float(bending1) * float(chi1) ** 2 + float(bending2) * float(chi2) ** 2) / denom
    return math.sqrt(max(float(value), 0.0))


class RootCache:
    def __init__(self, args: Args, runtime: RuntimeTracker | None = None) -> None:
        self.args = args
        self.runtime = runtime

    def settings(self, *, model: str, beta_deg: float, mu: float, eta: float, n_roots: int) -> dict[str, object]:
        return {
            "cache_version": CACHE_VERSION,
            "model": str(model),
            "epsilon_0": round(float(self.args.epsilon), 12),
            "beta_deg": round(float(beta_deg), 12),
            "mu": round(float(mu), 12),
            "eta": round(float(eta), 12),
            "n_roots": int(n_roots),
            "shape_sample_count": int(self.args.n_shape_points),
            "k_prime": float(TIMO.circular_shear_coefficient(TIMO.NU)),
            "root_settings": {
                "root_scan_start": beta_workflow.ROOT_SCAN_START,
                "eb_lambda_max": beta_workflow.EB_LAMBDA_MAX,
                "eb_retry_lambda_max": beta_workflow.EB_RETRY_LAMBDA_MAX,
                "eb_scan_step": beta_workflow.EB_SCAN_STEP,
                "eb_retry_scan_step": beta_workflow.EB_RETRY_SCAN_STEP,
                "timo_scan_step": beta_workflow.TIMO_SCAN_STEP,
                "timo_retry_scan_step": beta_workflow.TIMO_RETRY_SCAN_STEP,
                "retry_missing_roots": True,
            },
        }

    def path(self, *, model: str, beta_deg: float, mu: float, eta: float, n_roots: int) -> Path:
        settings_json = json.dumps(
            self.settings(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots),
            sort_keys=True,
        )
        digest = hashlib.sha256(settings_json.encode("utf-8")).hexdigest()[:20]
        token = "eb" if model == MODEL_EB else "timo"
        return self.args.cache_dir / f"root_{token}_{digest}.npz"

    def load(self, *, model: str, beta_deg: float, mu: float, eta: float, n_roots: int) -> RootCacheEntry | None:
        if self.args.force_recompute or not self.args.reuse_cache:
            return None
        path = self.path(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots)
        if not path.exists():
            return None
        data = np.load(path, allow_pickle=True)
        try:
            expected = json.dumps(
                self.settings(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots),
                sort_keys=True,
            )
            if str(data["settings_json"].item()) != expected:
                return None
            if self.runtime is not None:
                self.runtime.root_cache_hits += 1
            return RootCacheEntry(
                roots=tuple(float(value) for value in np.asarray(data["roots"], dtype=float)),
                warnings=tuple(json.loads(str(data["warnings_json"].item()))),
                root_count_found=int(data["root_count_found"].item()),
                lambda_max_used=float(data["lambda_max_used"].item()),
                scan_step_used=float(data["scan_step_used"].item()),
                retry_attempted=bool(data["retry_attempted"].item()),
                retry_changed_value=bool(data["retry_changed_value"].item()),
                notes=tuple(json.loads(str(data["notes_json"].item()))),
                cache_hit=True,
            )
        finally:
            data.close()

    def save(
        self,
        *,
        model: str,
        beta_deg: float,
        mu: float,
        eta: float,
        n_roots: int,
        result: beta_workflow.RootResult,
    ) -> None:
        path = self.path(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            settings_json=np.array(
                json.dumps(
                    self.settings(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots),
                    sort_keys=True,
                ),
                dtype=object,
            ),
            roots=np.asarray(result.roots, dtype=float),
            warnings_json=np.array(json.dumps(list(result.warnings)), dtype=object),
            notes_json=np.array(json.dumps(list(result.notes)), dtype=object),
            root_count_found=np.array(int(result.root_count_found), dtype=int),
            lambda_max_used=np.array(float(result.lambda_max_used), dtype=float),
            scan_step_used=np.array(float(result.scan_step_used), dtype=float),
            retry_attempted=np.array(bool(result.retry_attempted), dtype=bool),
            retry_changed_value=np.array(bool(result.retry_changed_value), dtype=bool),
        )

    def solve_uncached(
        self,
        *,
        model: str,
        beta_deg: float,
        mu: float,
        eta: float,
        n_roots: int,
        upper_hint: float | None,
    ) -> beta_workflow.RootResult:
        case = beta_workflow.CaseSpec(mu=float(mu), eta=float(eta), epsilon=float(self.args.epsilon))
        start = time.perf_counter()
        result = beta_workflow.solve_model(
            case,
            float(beta_deg),
            int(n_roots),
            model,
            upper_hint=upper_hint,
        )
        if result.root_count_found < int(n_roots):
            result = beta_workflow.retry_missing_roots(
                result,
                case,
                float(beta_deg),
                int(n_roots),
                model,
                upper_hint=upper_hint,
            )
        elapsed = time.perf_counter() - start
        if self.runtime is not None:
            if model == MODEL_EB:
                self.runtime.eb_roots_seconds += elapsed
                self.runtime.eb_root_calls += 1
            else:
                self.runtime.timo_roots_seconds += elapsed
                self.runtime.timo_root_calls += 1
        return result

    def roots(
        self,
        *,
        model: str,
        beta_deg: float,
        mu: float,
        eta: float,
        n_roots: int,
        upper_hint: float | None = None,
    ) -> RootCacheEntry:
        cached = self.load(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots)
        if cached is not None:
            return cached
        result = self.solve_uncached(
            model=model,
            beta_deg=beta_deg,
            mu=mu,
            eta=eta,
            n_roots=n_roots,
            upper_hint=upper_hint,
        )
        self.save(model=model, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots, result=result)
        if self.runtime is not None:
            self.runtime.root_cache_misses += 1
        return RootCacheEntry(
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


def eb_mode_result(
    *,
    epsilon: float,
    beta_deg: float,
    mu: float,
    eta: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    root_warnings: Sequence[str],
) -> ModeResult:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    l1, l2 = TIMO.segment_lengths(float(mu))
    section1 = TIMO.section_from_epsilon_tau(float(epsilon), factors.tau1)
    section2 = TIMO.section_from_epsilon_tau(float(epsilon), factors.tau2)
    matrix = assemble_clamped_coupled_matrix_eta(
        float(Lambda),
        math.radians(float(beta_deg)),
        float(mu),
        float(epsilon),
        float(eta),
    )
    coeff, _smallest, _ratio = analytic_null_vector(matrix)
    A1, B1, A2, B2, P1, P2 = [float(value) for value in np.asarray(coeff, dtype=float)]

    def rod_fields(x: np.ndarray, tau: float, A: float, B: float, P: float) -> RodFields:
        alpha = float(Lambda) / math.sqrt(float(tau))
        z = alpha * x
        omega = TIMO.project_omega(float(Lambda), float(epsilon))
        theta = omega
        u = P * np.sin(theta * x)
        u_prime = P * theta * np.cos(theta * x)
        w = A * (np.cos(z) - np.cosh(z)) + B * (np.sin(z) - np.sinh(z))
        w_prime = alpha * (A * (-np.sin(z) - np.sinh(z)) + B * (np.cos(z) - np.cosh(z)))
        w_second = alpha**2 * (A * (-np.cos(z) - np.cosh(z)) + B * (-np.sin(z) - np.sinh(z)))
        w_third = alpha**3 * (A * (np.sin(z) - np.sinh(z)) + B * (-np.cos(z) - np.cosh(z)))
        return RodFields(
            x=np.asarray(x, dtype=float),
            u=u,
            w=w,
            u_prime=u_prime,
            w_prime=w_prime,
            w_second=w_second,
            w_third=w_third,
            psi=None,
            psi_prime=None,
            gamma=None,
        )

    rod1 = rod_fields(np.linspace(0.0, l1, int(n_points), dtype=float), factors.tau1, A1, B1, P1)
    rod2 = rod_fields(np.linspace(0.0, -l2, int(n_points), dtype=float), factors.tau2, A2, B2, P2)

    x1_abs = np.abs(rod1.x)
    x2_abs = np.abs(rod2.x)
    U_axial = 0.5 * TIMO.E * section1.area * trapz(rod1.u_prime**2, x1_abs)
    U_axial += 0.5 * TIMO.E * section2.area * trapz(rod2.u_prime**2, x2_abs)
    U_bending = 0.5 * section1.bending_stiffness * trapz(np.asarray(rod1.w_second) ** 2, x1_abs)
    U_bending += 0.5 * section2.bending_stiffness * trapz(np.asarray(rod2.w_second) ** 2, x2_abs)
    return ModeResult(
        model=MODEL_EB,
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        mu=float(mu),
        eta=float(eta),
        sorted_index=int(sorted_index),
        Lambda=float(Lambda),
        coeff=np.asarray(coeff, dtype=float),
        rod1=rod1,
        rod2=rod2,
        energy=energy_from_terms(MODEL_EB, max(U_axial, 0.0), max(U_bending, 0.0), 0.0),
        warnings=tuple(root_warnings),
    )


def timo_mode_result(
    *,
    epsilon: float,
    beta_deg: float,
    mu: float,
    eta: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    root_warnings: Sequence[str],
) -> ModeResult:
    mode = TIMO.timo_mode_coefficients(float(Lambda), float(beta_deg), float(mu), float(epsilon), float(eta))
    fields = TIMO.timo_mode_fields(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=mode.coeff,
        n_points=int(n_points),
    )
    rod1_raw = fields["rod1"]
    rod2_raw = fields["rod2"]
    rod1 = RodFields(
        x=np.asarray(rod1_raw["x"], dtype=float),
        u=np.asarray(rod1_raw["u"], dtype=float),
        w=np.asarray(rod1_raw["w"], dtype=float),
        u_prime=np.asarray(rod1_raw["u_prime"], dtype=float),
        w_prime=np.asarray(rod1_raw["w_prime"], dtype=float),
        w_second=None,
        w_third=None,
        psi=np.asarray(rod1_raw["psi"], dtype=float),
        psi_prime=np.asarray(rod1_raw["psi_prime"], dtype=float),
        gamma=np.asarray(rod1_raw["w_prime"], dtype=float) - np.asarray(rod1_raw["psi"], dtype=float),
    )
    rod2 = RodFields(
        x=np.asarray(rod2_raw["x"], dtype=float),
        u=np.asarray(rod2_raw["u"], dtype=float),
        w=np.asarray(rod2_raw["w"], dtype=float),
        u_prime=np.asarray(rod2_raw["u_prime"], dtype=float),
        w_prime=np.asarray(rod2_raw["w_prime"], dtype=float),
        w_second=None,
        w_third=None,
        psi=np.asarray(rod2_raw["psi"], dtype=float),
        psi_prime=np.asarray(rod2_raw["psi_prime"], dtype=float),
        gamma=np.asarray(rod2_raw["w_prime"], dtype=float) - np.asarray(rod2_raw["psi"], dtype=float),
    )
    partition = TIMO.timo_energy_partition(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=mode.coeff,
        n_points=int(n_points),
    )
    warnings = list(root_warnings)
    warnings.extend(mode.warnings)
    warnings.extend(fields.get("warnings", ()))
    warnings = list(dict.fromkeys(item for item in warnings if item))
    return ModeResult(
        model=MODEL_TIMO,
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        mu=float(mu),
        eta=float(eta),
        sorted_index=int(sorted_index),
        Lambda=float(Lambda),
        coeff=np.asarray(mode.coeff, dtype=float),
        rod1=rod1,
        rod2=rod2,
        energy=energy_from_terms(
            MODEL_TIMO,
            max(float(partition["U_a_total"]), 0.0),
            max(float(partition["U_b_total"]), 0.0),
            max(float(partition["U_s_total"]), 0.0),
        ),
        warnings=tuple(warnings),
    )


def bending_energies_by_rod_eb(result: ModeResult) -> tuple[float, float]:
    factors = thickness_mismatch_factors(result.mu, result.eta)
    section1 = TIMO.section_from_epsilon_tau(result.epsilon, factors.tau1)
    section2 = TIMO.section_from_epsilon_tau(result.epsilon, factors.tau2)
    if result.rod1.w_second is None or result.rod2.w_second is None:
        return float("nan"), float("nan")
    e1 = 0.5 * section1.bending_stiffness * trapz(np.asarray(result.rod1.w_second) ** 2, np.abs(result.rod1.x))
    e2 = 0.5 * section2.bending_stiffness * trapz(np.asarray(result.rod2.w_second) ** 2, np.abs(result.rod2.x))
    return max(float(e1), 0.0), max(float(e2), 0.0)


def pi_eb_metrics(result: ModeResult) -> dict[str, float]:
    if result.model != MODEL_EB:
        raise ValueError("Pi_EB is defined here for EB modes only")
    factors = thickness_mismatch_factors(result.mu, result.eta)
    sections = {
        1: TIMO.section_from_epsilon_tau(result.epsilon, factors.tau1),
        2: TIMO.section_from_epsilon_tau(result.epsilon, factors.tau2),
    }
    omega = TIMO.project_omega(result.Lambda, result.epsilon)
    U_shear_star = 0.0
    T_trans = 0.0
    T_rot = 0.0
    U_total = float(result.energy.get("U_total", float("nan")))
    for rod_index, rod in ((1, result.rod1), (2, result.rod2)):
        if rod.w_third is None:
            raise ValueError("EB third derivative is required for Pi_EB")
        section = sections[rod_index]
        coordinate = np.abs(rod.x)
        q_eb = -section.bending_stiffness * np.asarray(rod.w_third, dtype=float)
        U_shear_star += 0.5 * trapz(q_eb**2 / section.shear_stiffness, coordinate)
        T_trans += 0.5 * omega**2 * section.mass_per_length * trapz(rod.u**2 + rod.w**2, coordinate)
        T_rot += 0.5 * omega**2 * section.rotary_inertia_per_length * trapz(rod.w_prime**2, coordinate)
    pi_shear = safe_ratio(U_shear_star, U_total)
    pi_rot = safe_ratio(T_rot, T_trans)
    pi = pi_shear + pi_rot if isfinite(pi_shear) and isfinite(pi_rot) else float("nan")
    return {
        "Pi_shear_EB": pi_shear,
        "Pi_rotary_EB": pi_rot,
        "Pi_EB": pi,
        "delta_f_from_Pi_half": 0.5 * pi if isfinite(pi) else float("nan"),
    }


def eb_predictors(result: ModeResult) -> dict[str, float]:
    start = time.perf_counter()
    chis = chi_values_eb(result.Lambda, result.epsilon, result.mu, result.eta)
    theta = theta_factor() * chis["chi_max_EB"] ** 2
    del start
    b1, b2 = bending_energies_by_rod_eb(result)
    chi_eff = chi_eff_from_eb_bending_energies(chis["chi_1_EB"], chis["chi_2_EB"], b1, b2)
    pi = pi_eb_metrics(result)
    return {
        **chis,
        "chi_eff_EB": chi_eff,
        "Theta_max_EB": theta,
        **pi,
    }


def trapezoid_weights(x: np.ndarray) -> np.ndarray:
    coordinates = np.abs(np.asarray(x, dtype=float))
    if coordinates.size == 1:
        return np.ones(1, dtype=float)
    dx = np.diff(coordinates)
    weights = np.zeros_like(coordinates, dtype=float)
    weights[0] = 0.5 * abs(float(dx[0]))
    weights[-1] = 0.5 * abs(float(dx[-1]))
    if coordinates.size > 2:
        weights[1:-1] = 0.5 * (np.abs(dx[:-1]) + np.abs(dx[1:]))
    return weights


def weighted_components(result: ModeResult, component: str) -> np.ndarray:
    factors = thickness_mismatch_factors(result.mu, result.eta)
    sections = {
        1: TIMO.section_from_epsilon_tau(result.epsilon, factors.tau1),
        2: TIMO.section_from_epsilon_tau(result.epsilon, factors.tau2),
    }
    out: list[np.ndarray] = []
    for rod_index, rod in ((1, result.rod1), (2, result.rod2)):
        section = sections[rod_index]
        if component == "uw":
            arrays = (rod.u, rod.w)
            density = section.mass_per_length
        elif component == "u":
            arrays = (rod.u,)
            density = section.mass_per_length
        elif component == "w":
            arrays = (rod.w,)
            density = section.mass_per_length
        elif component == "psi_slope":
            if result.model == MODEL_TIMO:
                if rod.psi is None:
                    arrays = (np.zeros_like(rod.u),)
                else:
                    arrays = (rod.psi,)
            else:
                arrays = (rod.w_prime,)
            density = section.rotary_inertia_per_length
        else:
            raise ValueError(f"unknown component: {component}")
        weights = np.sqrt(np.maximum(trapezoid_weights(rod.x) * density, 0.0))
        for array in arrays:
            out.append(weights * np.asarray(array, dtype=float))
    vector = np.concatenate(out) if out else np.array([], dtype=float)
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-30 or not np.all(np.isfinite(vector)):
        return np.full_like(vector, np.nan, dtype=float)
    return vector / norm


def mac_value(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return 0.0
    denom = float(np.dot(a, a) * np.dot(b, b))
    if denom <= 1.0e-30:
        return 0.0
    return abs(float(np.dot(a, b))) ** 2 / denom


def mac_metrics(eb: ModeResult, timo: ModeResult) -> dict[str, float]:
    return {
        "MAC_u": mac_value(weighted_components(eb, "u"), weighted_components(timo, "u")),
        "MAC_w": mac_value(weighted_components(eb, "w"), weighted_components(timo, "w")),
        "MAC_uw": mac_value(weighted_components(eb, "uw"), weighted_components(timo, "uw")),
        "MAC_psi_vs_EB_slope": mac_value(
            weighted_components(eb, "psi_slope"),
            weighted_components(timo, "psi_slope"),
        ),
    }


def weighted_uw_vector(result: ModeResult) -> np.ndarray:
    return weighted_components(result, "uw")


def orthonormal_columns(vectors: Sequence[np.ndarray]) -> np.ndarray:
    columns = [np.asarray(vector, dtype=float) for vector in vectors if np.all(np.isfinite(vector))]
    if not columns:
        return np.empty((0, 0), dtype=float)
    matrix = np.column_stack(columns)
    q, _r = np.linalg.qr(matrix)
    return q


def subspace_mac(left_vectors: Sequence[np.ndarray], right_vectors: Sequence[np.ndarray]) -> float:
    q_left = orthonormal_columns(left_vectors)
    q_right = orthonormal_columns(right_vectors)
    if q_left.size == 0 or q_right.size == 0 or q_left.shape[0] != q_right.shape[0]:
        return float("nan")
    singular_values = np.linalg.svd(q_left.T @ q_right, compute_uv=False)
    if singular_values.size == 0:
        return float("nan")
    return float(np.mean(np.clip(singular_values, 0.0, 1.0) ** 2))


def sorted_gap(roots: Sequence[float], index_one_based: int) -> float:
    values = [float(value) for value in roots if isfinite(float(value))]
    idx = int(index_one_based) - 1
    gaps: list[float] = []
    if idx > 0:
        left = values[idx - 1]
        current = values[idx]
        gaps.append(2.0 * abs(current - left) / (current + left))
    if idx < len(values) - 1:
        current = values[idx]
        right = values[idx + 1]
        gaps.append(2.0 * abs(right - current) / (right + current))
    return min(gaps) if gaps else float("nan")


def cluster_ids(roots: Sequence[float], threshold: float, max_index: int) -> dict[int, str]:
    ids: dict[int, str] = {}
    cluster_index = 1
    current: list[int] = [1]
    values = [float(value) for value in roots]
    for idx in range(1, min(len(values), int(max_index))):
        left = values[idx - 1]
        right = values[idx]
        gap = 2.0 * abs(right - left) / (right + left) if isfinite(left) and isfinite(right) else float("inf")
        if gap <= float(threshold):
            current.append(idx + 1)
        else:
            if len(current) > 1:
                cid = f"cluster_{cluster_index:02d}"
                for item in current:
                    ids[item] = cid
                cluster_index += 1
            current = [idx + 1]
    if len(current) > 1:
        cid = f"cluster_{cluster_index:02d}"
        for item in current:
            ids[item] = cid
    return ids


def mode_result_for_root(
    *,
    model: str,
    epsilon: float,
    beta_deg: float,
    mu: float,
    eta: float,
    sorted_index: int,
    Lambda: float,
    n_points: int,
    root_warnings: Sequence[str],
) -> ModeResult:
    if model == MODEL_EB:
        return eb_mode_result(
            epsilon=epsilon,
            beta_deg=beta_deg,
            mu=mu,
            eta=eta,
            sorted_index=sorted_index,
            Lambda=Lambda,
            n_points=n_points,
            root_warnings=root_warnings,
        )
    if model == MODEL_TIMO:
        return timo_mode_result(
            epsilon=epsilon,
            beta_deg=beta_deg,
            mu=mu,
            eta=eta,
            sorted_index=sorted_index,
            Lambda=Lambda,
            n_points=n_points,
            root_warnings=root_warnings,
        )
    raise ValueError(f"unknown model: {model}")


def solve_point_roots(
    args: Args,
    cache: RootCache,
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    n_roots: int,
) -> tuple[RootCacheEntry, RootCacheEntry]:
    eb = cache.roots(model=MODEL_EB, beta_deg=beta_deg, mu=mu, eta=eta, n_roots=n_roots)
    finite_eb = [float(value) for value in eb.roots if isfinite(float(value))]
    upper_hint = max(finite_eb) if finite_eb else None
    timo = cache.roots(
        model=MODEL_TIMO,
        beta_deg=beta_deg,
        mu=mu,
        eta=eta,
        n_roots=n_roots,
        upper_hint=upper_hint,
    )
    return eb, timo


def build_modes(
    args: Args,
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    eb_roots: Sequence[float],
    timo_roots: Sequence[float],
    eb_warnings: Sequence[str],
    timo_warnings: Sequence[str],
    n_timo_candidates: int,
    runtime: RuntimeTracker,
) -> tuple[list[ModeResult], list[ModeResult]]:
    eb_modes: list[ModeResult] = []
    timo_modes: list[ModeResult] = []
    start = time.perf_counter()
    for index, Lambda in enumerate(eb_roots[: args.n_reported_modes], start=1):
        if isfinite(float(Lambda)):
            eb_modes.append(
                mode_result_for_root(
                    model=MODEL_EB,
                    epsilon=args.epsilon,
                    beta_deg=beta_deg,
                    mu=mu,
                    eta=eta,
                    sorted_index=index,
                    Lambda=float(Lambda),
                    n_points=args.n_shape_points,
                    root_warnings=eb_warnings,
                )
            )
    for index, Lambda in enumerate(timo_roots[:n_timo_candidates], start=1):
        if isfinite(float(Lambda)):
            timo_modes.append(
                mode_result_for_root(
                    model=MODEL_TIMO,
                    epsilon=args.epsilon,
                    beta_deg=beta_deg,
                    mu=mu,
                    eta=eta,
                    sorted_index=index,
                    Lambda=float(Lambda),
                    n_points=args.n_shape_points,
                    root_warnings=timo_warnings,
                )
            )
    runtime.eb_shape_pi_seconds += time.perf_counter() - start
    return eb_modes, timo_modes


def assignment_for_modes(eb_modes: Sequence[ModeResult], timo_modes: Sequence[ModeResult]) -> tuple[np.ndarray, np.ndarray]:
    mac = np.array(
        [[mac_metrics(eb, timo)["MAC_uw"] for timo in timo_modes] for eb in eb_modes],
        dtype=float,
    )
    rows, cols = linear_sum_assignment(1.0 - mac)
    assignment = np.full(len(eb_modes), -1, dtype=int)
    for row, col in zip(rows, cols, strict=True):
        assignment[int(row)] = int(col)
    return assignment, mac


def matching_status(
    *,
    mac_uw: float,
    mac_margin: float,
    candidate_boundary: bool,
    root_warning: bool,
    cluster_member: bool,
    subspace_value: float,
) -> str:
    if candidate_boundary:
        return "candidate_boundary"
    if root_warning:
        return "root_warning"
    if cluster_member:
        if isfinite(float(subspace_value)) and float(subspace_value) >= MAC_MIN_RELIABLE:
            return "reliable_cluster"
        return "ambiguous"
    if not isfinite(float(mac_uw)):
        return "unmatched"
    if float(mac_uw) >= MAC_MIN_RELIABLE and isfinite(float(mac_margin)) and float(mac_margin) >= MAC_MARGIN_MIN_RELIABLE:
        return "reliable_individual"
    return "ambiguous"


def common_row_fields(epsilon: float, beta_deg: float, mu: float, eta: float) -> dict[str, float]:
    return {
        "epsilon_0": float(epsilon),
        "beta_deg": float(beta_deg),
        "mu": float(mu),
        "eta": float(eta),
        **local_geometry(float(epsilon), float(mu), float(eta)),
    }


def comparison_row(
    *,
    args: Args,
    beta_deg: float,
    mu: float,
    eta: float,
    comparison_type: str,
    eb: ModeResult,
    timo: ModeResult,
    macs: dict[str, float],
    mac_margin: float,
    status: str,
    include: bool,
    cluster_id: str,
    cluster_size_eb: int,
    cluster_size_timo: int,
    subspace_value: float,
    root_warning: bool,
    candidate_boundary: bool,
    notes: str,
) -> dict[str, object]:
    metrics = frequency_metrics(eb.Lambda, timo.Lambda)
    predictors = eb_predictors(eb)
    chis = chi_values_eb(eb.Lambda, args.epsilon, mu, eta)
    b1, b2 = bending_energies_by_rod_eb(eb)
    chi_eff = chi_eff_from_eb_bending_energies(chis["chi_1_EB"], chis["chi_2_EB"], b1, b2)
    predictors["chi_eff_EB"] = chi_eff
    row: dict[str, object] = {
        **common_row_fields(args.epsilon, beta_deg, mu, eta),
        "comparison_type": comparison_type,
        "eb_sorted_index": eb.sorted_index,
        "timo_sorted_index": timo.sorted_index,
        "cluster_id": cluster_id,
        "cluster_size_EB": int(cluster_size_eb),
        "cluster_size_Timo": int(cluster_size_timo),
        "subspace_MAC": subspace_value,
        "cluster_centroid_delta_f": float("nan"),
        "cluster_max_delta_f": float("nan"),
        "Pi_cluster_max": float("nan"),
        "Theta_cluster_max": float("nan"),
        "Lambda_EB": eb.Lambda,
        "Lambda_Timo": timo.Lambda,
        **metrics,
        **macs,
        "MAC_margin": mac_margin,
        "matching_status": status,
        "include_in_calibration": "true" if include else "false",
        "EB_axial_fraction": eb.energy["axial_fraction"],
        "EB_bending_fraction": eb.energy["bending_fraction"],
        "EB_classification": eb.energy["classification"],
        "Timo_axial_fraction": timo.energy["axial_fraction"],
        "Timo_bending_fraction": timo.energy["bending_fraction"],
        "Timo_shear_fraction": timo.energy["shear_fraction"],
        "Timo_bending_shear_fraction": float(timo.energy["bending_fraction"]) + float(timo.energy["shear_fraction"]),
        "Timo_classification": timo.energy["classification"],
        **predictors,
        "residual_Pi_half": metrics["delta_f"] - predictors["delta_f_from_Pi_half"]
        if isfinite(metrics["delta_f"]) and isfinite(predictors["delta_f_from_Pi_half"])
        else float("nan"),
        "pass_10pct": "true" if pass_bool(metrics["delta_f"]) else "false",
        "root_warning": "true" if root_warning else "false",
        "candidate_boundary": "true" if candidate_boundary else "false",
        "notes": notes,
    }
    return row


def cluster_summary_rows(
    *,
    args: Args,
    beta_deg: float,
    mu: float,
    eta: float,
    eb_modes: Sequence[ModeResult],
    timo_modes: Sequence[ModeResult],
    assignment: Sequence[int],
    individual_rows: Sequence[dict[str, object]],
    eb_clusters: dict[int, str],
    timo_clusters: dict[int, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    cluster_names = sorted(set(eb_clusters.values()))
    for name in cluster_names:
        eb_indices = sorted(index for index, cid in eb_clusters.items() if cid == name)
        timo_cols = [int(assignment[index - 1]) for index in eb_indices if 0 <= index - 1 < len(assignment)]
        timo_indices = sorted({timo_modes[col].sorted_index for col in timo_cols if 0 <= col < len(timo_modes)})
        if not eb_indices or not timo_indices:
            continue
        eb_cluster_modes = [eb_modes[index - 1] for index in eb_indices if 0 <= index - 1 < len(eb_modes)]
        timo_cluster_modes = [timo_modes[index - 1] for index in timo_indices if 0 <= index - 1 < len(timo_modes)]
        if not eb_cluster_modes or not timo_cluster_modes:
            continue
        sub = subspace_mac([weighted_uw_vector(item) for item in eb_cluster_modes], [weighted_uw_vector(item) for item in timo_cluster_modes])
        lambda_eb = float(np.mean([mode.Lambda for mode in eb_cluster_modes]))
        lambda_timo = float(np.mean([mode.Lambda for mode in timo_cluster_modes]))
        metrics = frequency_metrics(lambda_eb, lambda_timo)
        related_rows = [row for row in individual_rows if int(row["eb_sorted_index"]) in eb_indices]
        max_delta = max((finite_float(row, "delta_f") for row in related_rows), default=float("nan"))
        pi_max = max((finite_float(row, "Pi_EB") for row in related_rows), default=float("nan"))
        theta_max = max((finite_float(row, "Theta_max_EB") for row in related_rows), default=float("nan"))
        representative_eb = eb_cluster_modes[0]
        representative_timo = timo_cluster_modes[0]
        row = {
            **common_row_fields(args.epsilon, beta_deg, mu, eta),
            "comparison_type": "matched_cluster",
            "eb_sorted_index": ",".join(str(index) for index in eb_indices),
            "timo_sorted_index": ",".join(str(index) for index in timo_indices),
            "cluster_id": name,
            "cluster_size_EB": len(eb_indices),
            "cluster_size_Timo": len(timo_indices),
            "subspace_MAC": sub,
            "cluster_centroid_delta_f": metrics["delta_f"],
            "cluster_max_delta_f": max_delta,
            "Pi_cluster_max": pi_max,
            "Theta_cluster_max": theta_max,
            "Lambda_EB": lambda_eb,
            "Lambda_Timo": lambda_timo,
            **metrics,
            "MAC_u": float("nan"),
            "MAC_w": float("nan"),
            "MAC_uw": float("nan"),
            "MAC_psi_vs_EB_slope": float("nan"),
            "MAC_margin": float("nan"),
            "matching_status": "reliable_cluster" if isfinite(sub) and sub >= MAC_MIN_RELIABLE else "ambiguous",
            "include_in_calibration": "true" if isfinite(sub) and sub >= MAC_MIN_RELIABLE else "false",
            "EB_axial_fraction": representative_eb.energy["axial_fraction"],
            "EB_bending_fraction": representative_eb.energy["bending_fraction"],
            "EB_classification": representative_eb.energy["classification"],
            "Timo_axial_fraction": representative_timo.energy["axial_fraction"],
            "Timo_bending_fraction": representative_timo.energy["bending_fraction"],
            "Timo_shear_fraction": representative_timo.energy["shear_fraction"],
            "Timo_bending_shear_fraction": float(representative_timo.energy["bending_fraction"]) + float(representative_timo.energy["shear_fraction"]),
            "Timo_classification": "cluster_level",
            "chi_1_EB": max((finite_float(row, "chi_1_EB") for row in related_rows), default=float("nan")),
            "chi_2_EB": max((finite_float(row, "chi_2_EB") for row in related_rows), default=float("nan")),
            "chi_max_EB": max((finite_float(row, "chi_max_EB") for row in related_rows), default=float("nan")),
            "chi_eff_EB": max((finite_float(row, "chi_eff_EB") for row in related_rows), default=float("nan")),
            "Theta_max_EB": theta_max,
            "Pi_shear_EB": max((finite_float(row, "Pi_shear_EB") for row in related_rows), default=float("nan")),
            "Pi_rotary_EB": max((finite_float(row, "Pi_rotary_EB") for row in related_rows), default=float("nan")),
            "Pi_EB": pi_max,
            "delta_f_from_Pi_half": 0.5 * pi_max if isfinite(pi_max) else float("nan"),
            "residual_Pi_half": max_delta - 0.5 * pi_max if isfinite(max_delta) and isfinite(pi_max) else float("nan"),
            "pass_10pct": "true" if pass_bool(max_delta) else "false",
            "root_warning": "false",
            "candidate_boundary": "false",
            "notes": "cluster-level conservative row; predictors are maxima over EB cluster members",
        }
        rows.append(row)
    del timo_clusters
    return rows


def build_matching_rows(
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    eb_roots: Sequence[float],
    timo_roots: Sequence[float],
    eb_modes: Sequence[ModeResult],
    timo_modes: Sequence[ModeResult],
    assignment: Sequence[int],
    mac_matrix: np.ndarray,
    eb_clusters: dict[int, str],
    timo_clusters: dict[int, str],
    status_by_eb: dict[int, str],
    selected_margin: dict[int, float],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for eb_pos, eb in enumerate(eb_modes):
        mac_row = mac_matrix[eb_pos] if eb_pos < mac_matrix.shape[0] else np.array([], dtype=float)
        selected_col = int(assignment[eb_pos]) if eb_pos < len(assignment) else -1
        for col, timo in enumerate(timo_modes):
            macs = mac_metrics(eb, timo)
            selected = col == selected_col
            rows.append(
                {
                    "beta_deg": beta_deg,
                    "mu": mu,
                    "eta": eta,
                    "eb_sorted_index": eb.sorted_index,
                    "timo_candidate_index": timo.sorted_index,
                    "Lambda_EB": eb.Lambda,
                    "Lambda_Timo": timo.Lambda,
                    **macs,
                    "MAC_margin": selected_margin.get(eb.sorted_index, float("nan")) if selected else float("nan"),
                    "gap_EB": sorted_gap(eb_roots, eb.sorted_index),
                    "gap_Timo": sorted_gap(timo_roots, timo.sorted_index),
                    "cluster_EB": eb_clusters.get(eb.sorted_index, ""),
                    "cluster_Timo": timo_clusters.get(timo.sorted_index, ""),
                    "subspace_MAC": float("nan"),
                    "assignment_selected": "true" if selected else "false",
                    "matching_status": status_by_eb.get(eb.sorted_index, "") if selected else "",
                    "warning": "" if selected else "",
                    "notes": "Hungarian one-to-one assignment over mass-weighted local u,w MAC",
                }
            )
        if selected_col < 0 or selected_col >= len(timo_modes) or mac_row.size == 0:
            rows.append(
                {
                    "beta_deg": beta_deg,
                    "mu": mu,
                    "eta": eta,
                    "eb_sorted_index": eb.sorted_index,
                    "timo_candidate_index": "",
                    "Lambda_EB": eb.Lambda,
                    "Lambda_Timo": float("nan"),
                    "MAC_u": float("nan"),
                    "MAC_w": float("nan"),
                    "MAC_uw": float("nan"),
                    "MAC_margin": float("nan"),
                    "gap_EB": sorted_gap(eb_roots, eb.sorted_index),
                    "gap_Timo": float("nan"),
                    "cluster_EB": eb_clusters.get(eb.sorted_index, ""),
                    "cluster_Timo": "",
                    "subspace_MAC": float("nan"),
                    "assignment_selected": "false",
                    "matching_status": "unmatched",
                    "warning": "unmatched EB mode",
                    "notes": "assignment failed",
                }
            )
    return rows


def point_products(
    args: Args,
    cache: RootCache,
    runtime: RuntimeTracker,
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    n_roots: int | None = None,
) -> PointProducts:
    n_candidate = int(args.n_candidate_roots if n_roots is None else n_roots)
    eb_entry, timo_entry = solve_point_roots(
        args,
        cache,
        beta_deg=beta_deg,
        mu=mu,
        eta=eta,
        n_roots=n_candidate,
    )
    eb_modes, timo_modes = build_modes(
        args,
        beta_deg=beta_deg,
        mu=mu,
        eta=eta,
        eb_roots=eb_entry.roots,
        timo_roots=timo_entry.roots,
        eb_warnings=eb_entry.warnings,
        timo_warnings=timo_entry.warnings,
        n_timo_candidates=n_candidate,
        runtime=runtime,
    )
    start_matching = time.perf_counter()
    assignment, mac_matrix = assignment_for_modes(eb_modes, timo_modes)
    runtime.matching_seconds += time.perf_counter() - start_matching

    boundary_selected = any(
        int(col) == n_candidate - 1 for col in assignment if int(col) >= 0
    )
    if boundary_selected and n_candidate < BOUNDARY_RETRY_ROOTS:
        runtime.boundary_retries += 1
        return point_products(
            args,
            cache,
            runtime,
            beta_deg=beta_deg,
            mu=mu,
            eta=eta,
            n_roots=BOUNDARY_RETRY_ROOTS,
        )

    eb_clusters = cluster_ids(eb_entry.roots, args.cluster_gap_threshold, args.n_reported_modes)
    timo_clusters = cluster_ids(timo_entry.roots, args.cluster_gap_threshold, n_candidate)
    root_warning = bool(eb_entry.warnings or timo_entry.warnings)
    rows: list[dict[str, object]] = []
    status_by_eb: dict[int, str] = {}
    selected_margin: dict[int, float] = {}

    for eb_pos, eb in enumerate(eb_modes):
        selected_col = int(assignment[eb_pos]) if eb_pos < len(assignment) else -1
        if selected_col < 0 or selected_col >= len(timo_modes):
            continue
        timo = timo_modes[selected_col]
        macs = mac_metrics(eb, timo)
        row_values = np.sort(mac_matrix[eb_pos])[::-1]
        selected_mac = float(mac_matrix[eb_pos, selected_col])
        second_best = float(row_values[1]) if row_values.size > 1 else float("nan")
        margin = selected_mac - second_best if isfinite(second_best) else float("nan")
        cluster_member = eb.sorted_index in eb_clusters or timo.sorted_index in timo_clusters
        if cluster_member:
            eb_members = [idx for idx, cid in eb_clusters.items() if cid == eb_clusters.get(eb.sorted_index, "")]
            if not eb_members:
                eb_members = [eb.sorted_index]
            timo_members = [idx for idx, cid in timo_clusters.items() if cid == timo_clusters.get(timo.sorted_index, "")]
            if not timo_members:
                timo_members = [timo.sorted_index]
            sub = subspace_mac(
                [weighted_uw_vector(eb_modes[idx - 1]) for idx in eb_members if 0 <= idx - 1 < len(eb_modes)],
                [weighted_uw_vector(timo_modes[idx - 1]) for idx in timo_members if 0 <= idx - 1 < len(timo_modes)],
            )
        else:
            sub = float("nan")
        candidate_boundary = selected_col == n_candidate - 1
        status = matching_status(
            mac_uw=macs["MAC_uw"],
            mac_margin=margin,
            candidate_boundary=candidate_boundary,
            root_warning=root_warning,
            cluster_member=cluster_member,
            subspace_value=sub,
        )
        include = status == "reliable_individual"
        status_by_eb[eb.sorted_index] = status
        selected_margin[eb.sorted_index] = margin
        rows.append(
            comparison_row(
                args=args,
                beta_deg=beta_deg,
                mu=mu,
                eta=eta,
                comparison_type="homologous_mode",
                eb=eb,
                timo=timo,
                macs=macs,
                mac_margin=margin,
                status=status,
                include=include,
                cluster_id=eb_clusters.get(eb.sorted_index, timo_clusters.get(timo.sorted_index, "")),
                cluster_size_eb=sum(1 for cid in eb_clusters.values() if cid == eb_clusters.get(eb.sorted_index, "")) if eb.sorted_index in eb_clusters else 1,
                cluster_size_timo=sum(1 for cid in timo_clusters.values() if cid == timo_clusters.get(timo.sorted_index, "")) if timo.sorted_index in timo_clusters else 1,
                subspace_value=sub,
                root_warning=root_warning,
                candidate_boundary=candidate_boundary,
                notes="homologous EB/Timoshenko mode from per-point mass-weighted MAC assignment",
            )
        )

    sorted_rows: list[dict[str, object]] = []
    for index in range(1, min(args.n_reported_modes, len(eb_modes), len(timo_modes)) + 1):
        eb = eb_modes[index - 1]
        timo = timo_modes[index - 1]
        sorted_rows.append(
            comparison_row(
                args=args,
                beta_deg=beta_deg,
                mu=mu,
                eta=eta,
                comparison_type="sorted_index",
                eb=eb,
                timo=timo,
                macs=mac_metrics(eb, timo),
                mac_margin=float("nan"),
                status="sorted_comparison",
                include=False,
                cluster_id="",
                cluster_size_eb=1,
                cluster_size_timo=1,
                subspace_value=float("nan"),
                root_warning=root_warning,
                candidate_boundary=False,
                notes="sorted-spectrum diagnostic; no mode identity is implied",
            )
        )

    cluster_rows = cluster_summary_rows(
        args=args,
        beta_deg=beta_deg,
        mu=mu,
        eta=eta,
        eb_modes=eb_modes,
        timo_modes=timo_modes,
        assignment=assignment,
        individual_rows=rows,
        eb_clusters=eb_clusters,
        timo_clusters=timo_clusters,
    )
    matching_rows = build_matching_rows(
        beta_deg=beta_deg,
        mu=mu,
        eta=eta,
        eb_roots=eb_entry.roots,
        timo_roots=timo_entry.roots,
        eb_modes=eb_modes,
        timo_modes=timo_modes,
        assignment=assignment,
        mac_matrix=mac_matrix,
        eb_clusters=eb_clusters,
        timo_clusters=timo_clusters,
        status_by_eb=status_by_eb,
        selected_margin=selected_margin,
    )
    mode_rows = [*sorted_rows, *rows, *cluster_rows]
    sorted_deltas = [finite_float(row, "delta_f") for row in sorted_rows]
    bending_sorted = [
        row for row in sorted_rows if str(row.get("Timo_classification", "")) == "bending_dominated"
    ]
    bending_deltas = [finite_float(row, "delta_f") for row in bending_sorted]
    summary = {
        "epsilon_0": args.epsilon,
        "beta_deg": beta_deg,
        "mu": mu,
        "eta": eta,
        **k_aware_point_summary_fields(
            sorted_deltas,
            target_mode_count=args.n_reported_modes,
        ),
        "N10_current_bending": consecutive_pass_count(bending_deltas),
        "n_bending_pass_10": count_passes(bending_deltas),
        "reliable_individual_count": sum(1 for row in rows if row["matching_status"] == "reliable_individual"),
        "reliable_cluster_count": sum(1 for row in [*rows, *cluster_rows] if row["matching_status"] == "reliable_cluster"),
        "ambiguous_count": sum(1 for row in rows if row["matching_status"] == "ambiguous"),
        "root_warning_count": sum(1 for row in rows if row["matching_status"] == "root_warning"),
        "candidate_boundary_count": sum(1 for row in rows if row["matching_status"] == "candidate_boundary"),
        "cluster_count_EB": len(set(eb_clusters.values())),
        "cluster_count_Timo": len(set(timo_clusters.values())),
        "matched_cluster_count": len(cluster_rows),
        "notes": "current bending metric uses current Timoshenko sorted energy classification",
    }
    local_rows = (
        local_timo_benchmark_rows(args, cache, beta_deg=beta_deg, mu=mu, eta=eta, eb_roots=eb_entry.roots, full_rows=rows, runtime=runtime)
        if args.benchmark_local_timo
        else []
    )
    return PointProducts(mode_rows=mode_rows, matching_rows=matching_rows, summary_row=summary, local_timo_rows=local_rows)


def subgroup_filter(rows: Sequence[dict[str, object]], subgroup: str) -> list[dict[str, object]]:
    if subgroup == "cluster_level":
        return [
            row
            for row in rows
            if row.get("comparison_type") == "matched_cluster"
            and str(row.get("include_in_calibration", "")).lower() == "true"
        ]
    base = [
        row
        for row in rows
        if row.get("comparison_type") == "homologous_mode"
        and str(row.get("include_in_calibration", "")).lower() == "true"
    ]
    if subgroup == "all_reliable":
        return base
    if subgroup == "bending":
        return [row for row in base if row.get("Timo_classification") == "bending_dominated"]
    if subgroup == "mixed":
        return [row for row in base if row.get("Timo_classification") == "mixed"]
    if subgroup == "longitudinal":
        return [row for row in base if row.get("Timo_classification") == "longitudinal_dominated"]
    raise ValueError(f"unknown subgroup {subgroup}")


def finite_xy(rows: Sequence[dict[str, object]], predictor: str) -> tuple[np.ndarray, np.ndarray]:
    x_values: list[float] = []
    y_values: list[float] = []
    for row in rows:
        x = finite_float(row, predictor)
        y = finite_float(row, "delta_f")
        if isfinite(x) and isfinite(y):
            x_values.append(x)
            y_values.append(y)
    return np.asarray(x_values, dtype=float), np.asarray(y_values, dtype=float)


def fit_power_law_arrays(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    mask = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > NUMERICAL_FLOOR)
    xf = np.asarray(x[mask], dtype=float)
    yf = np.asarray(y[mask], dtype=float)
    if xf.size < 3 or float(np.std(np.log(xf))) <= 1.0e-14:
        return {
            "C": float("nan"),
            "p": float("nan"),
            "R2": float("nan"),
            "RMSE": float("nan"),
            "MAE": float("nan"),
            "threshold_10": float("nan"),
        }
    p, logc = np.polyfit(np.log(xf), np.log(yf), 1)
    c = float(np.exp(logc))
    pred_log = float(p) * np.log(xf) + float(logc)
    pred = c * xf ** float(p)
    ss_res = float(np.sum((np.log(yf) - pred_log) ** 2))
    ss_tot = float(np.sum((np.log(yf) - float(np.mean(np.log(yf)))) ** 2))
    threshold = float((MAIN_THRESHOLD / c) ** (1.0 / float(p))) if c > 0.0 and p > 0.0 else float("nan")
    return {
        "C": c,
        "p": float(p),
        "R2": 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan"),
        "RMSE": float(np.sqrt(np.mean((pred - yf) ** 2))),
        "MAE": float(np.mean(np.abs(pred - yf))),
        "threshold_10": threshold,
    }


def scalar_prediction(row: dict[str, object], fit: dict[str, float], predictor: str) -> float:
    x = finite_float(row, predictor)
    c = float(fit.get("C", float("nan")))
    p = float(fit.get("p", float("nan")))
    if not (isfinite(x) and x > 0.0 and isfinite(c) and c > 0.0 and isfinite(p)):
        return float("nan")
    return c * x**p


def classification_metrics(
    rows: Sequence[dict[str, object]],
    *,
    predictor: str,
    threshold: float | None = None,
    predictions: Sequence[float] | None = None,
) -> dict[str, float | int]:
    finite_rows: list[tuple[dict[str, object], bool]] = []
    if predictions is None:
        if threshold is None or not isfinite(float(threshold)):
            return {
                "row_count": 0,
                "false_safe_count": 0,
                "false_safe_rate": float("nan"),
                "false_unsafe_count": 0,
                "false_unsafe_rate": float("nan"),
                "balanced_accuracy": float("nan"),
                "safe_coverage": float("nan"),
            }
        for row in rows:
            x = finite_float(row, predictor)
            y = finite_float(row, "delta_f")
            if isfinite(x) and isfinite(y):
                finite_rows.append((row, x <= float(threshold)))
    else:
        for row, pred in zip(rows, predictions, strict=False):
            y = finite_float(row, "delta_f")
            if isfinite(float(pred)) and isfinite(y):
                finite_rows.append((row, float(pred) <= MAIN_THRESHOLD))
    actual_safe = [(row, pred_safe) for row, pred_safe in finite_rows if finite_float(row, "delta_f") <= MAIN_THRESHOLD]
    actual_unsafe = [(row, pred_safe) for row, pred_safe in finite_rows if finite_float(row, "delta_f") > MAIN_THRESHOLD]
    true_safe = [(row, pred_safe) for row, pred_safe in actual_safe if pred_safe]
    true_unsafe = [(row, pred_safe) for row, pred_safe in actual_unsafe if not pred_safe]
    false_safe = [(row, pred_safe) for row, pred_safe in actual_unsafe if pred_safe]
    false_unsafe = [(row, pred_safe) for row, pred_safe in actual_safe if not pred_safe]
    tpr_safe = safe_ratio(len(true_safe), len(actual_safe))
    tpr_unsafe = safe_ratio(len(true_unsafe), len(actual_unsafe))
    return {
        "row_count": len(finite_rows),
        "false_safe_count": len(false_safe),
        "false_safe_rate": safe_ratio(len(false_safe), len(actual_unsafe)),
        "false_unsafe_count": len(false_unsafe),
        "false_unsafe_rate": safe_ratio(len(false_unsafe), len(actual_safe)),
        "balanced_accuracy": 0.5 * (tpr_safe + tpr_unsafe) if isfinite(tpr_safe) and isfinite(tpr_unsafe) else float("nan"),
        "safe_coverage": safe_ratio(len(true_safe), len(actual_safe)),
    }


def conservative_threshold(rows: Sequence[dict[str, object]], predictor: str, max_false_safe_rate: float) -> float:
    values = sorted({finite_float(row, predictor) for row in rows if isfinite(finite_float(row, predictor))})
    if not values:
        return float("nan")
    best = float("nan")
    for threshold in values:
        cls = classification_metrics(rows, predictor=predictor, threshold=threshold)
        rate = float(cls["false_safe_rate"])
        if isfinite(rate) and rate <= float(max_false_safe_rate) + 1.0e-15:
            best = float(threshold)
    return best


def raw_features(rows: Sequence[dict[str, object]]) -> np.ndarray:
    data: list[list[float]] = []
    for row in rows:
        beta = math.radians(finite_float(row, "beta_deg"))
        mu = finite_float(row, "mu")
        eta = finite_float(row, "eta")
        lam = finite_float(row, "Lambda_EB")
        data.append(
            [
                1.0,
                mu,
                eta,
                math.sin(beta),
                math.cos(beta),
                lam,
                lam * lam,
                mu * eta,
                mu * math.sin(beta),
                eta * math.sin(beta),
            ]
        )
    return np.asarray(data, dtype=float)


def fit_raw_baseline(rows: Sequence[dict[str, object]]) -> dict[str, object]:
    usable = [row for row in rows if finite_float(row, "delta_f") > NUMERICAL_FLOOR and isfinite(finite_float(row, "Lambda_EB"))]
    if len(usable) < 6:
        return {"coef": np.array([], dtype=float), "R2": float("nan"), "RMSE": float("nan"), "MAE": float("nan")}
    x = raw_features(usable)
    y = np.log(np.asarray([finite_float(row, "delta_f") for row in usable], dtype=float))
    coef = np.linalg.pinv(x.T @ x + 1.0e-10 * np.eye(x.shape[1])) @ x.T @ y
    pred = np.exp(x @ coef)
    actual = np.exp(y)
    ss_res = float(np.sum((y - x @ coef) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    return {
        "coef": coef,
        "R2": 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan"),
        "RMSE": float(np.sqrt(np.mean((pred - actual) ** 2))),
        "MAE": float(np.mean(np.abs(pred - actual))),
    }


def raw_predict(rows: Sequence[dict[str, object]], fit: dict[str, object]) -> np.ndarray:
    coef = np.asarray(fit.get("coef", np.array([], dtype=float)), dtype=float)
    if coef.size == 0:
        return np.full(len(rows), np.nan, dtype=float)
    return np.exp(raw_features(rows) @ coef)


def fit_summary_rows(mode_rows: Sequence[dict[str, object]]) -> tuple[list[dict[str, object]], dict[tuple[str, str], dict[str, float]]]:
    rows_out: list[dict[str, object]] = []
    scalar_fits: dict[tuple[str, str], dict[str, float]] = {}
    for subgroup in SUBGROUPS:
        rows = subgroup_filter(mode_rows, subgroup)
        for predictor in PREDICTORS:
            if predictor == "raw_parameter_baseline":
                raw_fit = fit_raw_baseline(rows)
                preds = raw_predict(rows, raw_fit)
                cls = classification_metrics(rows, predictor=predictor, predictions=preds)
                rows_out.append(
                    {
                        "predictor": predictor,
                        "subgroup": subgroup,
                        "sample_filter": "reliable homologous rows; multivariate log baseline",
                        "point_count": len(rows),
                        "C": float("nan"),
                        "p": float("nan"),
                        "R2": raw_fit["R2"],
                        "RMSE": raw_fit["RMSE"],
                        "MAE": raw_fit["MAE"],
                        "threshold_10_central": float("nan"),
                        "threshold_safe_0pct": float("nan"),
                        "threshold_safe_1pct": float("nan"),
                        "threshold_safe_5pct": float("nan"),
                        "false_safe_rate_central": cls["false_safe_rate"],
                        "false_unsafe_rate_central": cls["false_unsafe_rate"],
                        "balanced_accuracy_central": cls["balanced_accuracy"],
                        "safe_coverage_central": cls["safe_coverage"],
                        "geometry_extended_R2": raw_fit["R2"],
                        "geometry_extended_RMSE": raw_fit["RMSE"],
                        "geometry_extended_MAE": raw_fit["MAE"],
                        "notes": "raw-parameter OLS/ridge-like baseline predicts delta_f directly; no scalar threshold",
                    }
                )
                continue
            x, y = finite_xy(rows, predictor)
            fit = fit_power_law_arrays(x, y)
            scalar_fits[(subgroup, predictor)] = fit
            cls = classification_metrics(rows, predictor=predictor, threshold=fit["threshold_10"])
            thresholds = {
                "threshold_safe_0pct": conservative_threshold(rows, predictor, 0.0),
                "threshold_safe_1pct": conservative_threshold(rows, predictor, 0.01),
                "threshold_safe_5pct": conservative_threshold(rows, predictor, 0.05),
            }
            extended = fit_geometry_extended(rows, predictor)
            rows_out.append(
                {
                    "predictor": predictor,
                    "subgroup": subgroup,
                    "sample_filter": "reliable individual rows" if subgroup != "cluster_level" else "reliable cluster rows",
                    "point_count": len(x),
                    "C": fit["C"],
                    "p": fit["p"],
                    "R2": fit["R2"],
                    "RMSE": fit["RMSE"],
                    "MAE": fit["MAE"],
                    "threshold_10_central": fit["threshold_10"],
                    **thresholds,
                    "false_safe_rate_central": cls["false_safe_rate"],
                    "false_unsafe_rate_central": cls["false_unsafe_rate"],
                    "balanced_accuracy_central": cls["balanced_accuracy"],
                    "safe_coverage_central": cls["safe_coverage"],
                    "geometry_extended_R2": extended["R2"],
                    "geometry_extended_RMSE": extended["RMSE"],
                    "geometry_extended_MAE": extended["MAE"],
                    "notes": "R2 is log-space; RMSE/MAE are in delta_f space",
                }
            )
    return rows_out, scalar_fits


def fit_geometry_extended(rows: Sequence[dict[str, object]], predictor: str) -> dict[str, float]:
    usable = [
        row
        for row in rows
        if finite_float(row, predictor) > 0.0 and finite_float(row, "delta_f") > NUMERICAL_FLOOR
    ]
    if len(usable) < 8:
        return {"R2": float("nan"), "RMSE": float("nan"), "MAE": float("nan")}
    design: list[list[float]] = []
    y_values: list[float] = []
    for row in usable:
        beta = math.radians(finite_float(row, "beta_deg"))
        mu = finite_float(row, "mu")
        eta = finite_float(row, "eta")
        x = finite_float(row, predictor)
        design.append([1.0, math.log(x), math.sin(beta), math.cos(beta), mu, eta, mu * math.sin(beta), eta * math.sin(beta), mu * eta])
        y_values.append(math.log(finite_float(row, "delta_f")))
    xmat = np.asarray(design, dtype=float)
    y = np.asarray(y_values, dtype=float)
    coef = np.linalg.pinv(xmat.T @ xmat + 1.0e-10 * np.eye(xmat.shape[1])) @ xmat.T @ y
    pred_log = xmat @ coef
    pred = np.exp(pred_log)
    actual = np.exp(y)
    ss_res = float(np.sum((y - pred_log) ** 2))
    ss_tot = float(np.sum((y - float(np.mean(y))) ** 2))
    return {
        "R2": 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan"),
        "RMSE": float(np.sqrt(np.mean((pred - actual) ** 2))),
        "MAE": float(np.mean(np.abs(pred - actual))),
    }


def cross_validation_rows(mode_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for subgroup in SUBGROUPS:
        rows = subgroup_filter(mode_rows, subgroup)
        for split in SPLITS:
            key = {"leave_one_beta_out": "beta_deg", "leave_one_mu_out": "mu", "leave_one_eta_out": "eta"}[split]
            held_values = sorted({finite_float(row, key) for row in rows if isfinite(finite_float(row, key))})
            for predictor in PREDICTORS:
                for held in held_values:
                    train = [row for row in rows if abs(finite_float(row, key) - held) > 1.0e-12]
                    test = [row for row in rows if abs(finite_float(row, key) - held) <= 1.0e-12]
                    if predictor == "raw_parameter_baseline":
                        fit = fit_raw_baseline(train)
                        pred = raw_predict(test, fit)
                        actual = np.asarray([finite_float(row, "delta_f") for row in test], dtype=float)
                        mask = np.isfinite(pred) & np.isfinite(actual)
                        rmse = float(np.sqrt(np.mean((pred[mask] - actual[mask]) ** 2))) if np.any(mask) else float("nan")
                        mae = float(np.mean(np.abs(pred[mask] - actual[mask]))) if np.any(mask) else float("nan")
                        cls = classification_metrics(test, predictor=predictor, predictions=pred)
                        out.append(
                            {
                                "predictor": predictor,
                                "subgroup": subgroup,
                                "split_kind": split,
                                "held_out_value": held,
                                "train_count": len(train),
                                "test_count": len(test),
                                "C": float("nan"),
                                "p": float("nan"),
                                "threshold_10": float("nan"),
                                "test_RMSE": rmse,
                                "test_MAE": mae,
                                "false_safe_rate": cls["false_safe_rate"],
                                "false_unsafe_rate": cls["false_unsafe_rate"],
                                "balanced_accuracy": cls["balanced_accuracy"],
                                "safe_coverage": cls["safe_coverage"],
                                "notes": "raw-parameter baseline leave-one-group-out",
                            }
                        )
                        continue
                    x_train, y_train = finite_xy(train, predictor)
                    fit = fit_power_law_arrays(x_train, y_train)
                    pred_values = np.asarray([scalar_prediction(row, fit, predictor) for row in test], dtype=float)
                    y_test = np.asarray([finite_float(row, "delta_f") for row in test], dtype=float)
                    mask = np.isfinite(pred_values) & np.isfinite(y_test)
                    rmse = float(np.sqrt(np.mean((pred_values[mask] - y_test[mask]) ** 2))) if np.any(mask) else float("nan")
                    mae = float(np.mean(np.abs(pred_values[mask] - y_test[mask]))) if np.any(mask) else float("nan")
                    cls = classification_metrics(test, predictor=predictor, threshold=fit["threshold_10"])
                    out.append(
                        {
                            "predictor": predictor,
                            "subgroup": subgroup,
                            "split_kind": split,
                            "held_out_value": held,
                            "train_count": len(train),
                            "test_count": len(test),
                            "C": fit["C"],
                            "p": fit["p"],
                            "threshold_10": fit["threshold_10"],
                            "test_RMSE": rmse,
                            "test_MAE": mae,
                            "false_safe_rate": cls["false_safe_rate"],
                            "false_unsafe_rate": cls["false_unsafe_rate"],
                            "balanced_accuracy": cls["balanced_accuracy"],
                            "safe_coverage": cls["safe_coverage"],
                            "notes": "scalar predictor leave-one-group-out",
                        }
                    )
    return out


def threshold_classification_rows(mode_rows: Sequence[dict[str, object]], fit_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    fit_map = {(row["subgroup"], row["predictor"]): row for row in fit_rows}
    for subgroup in SUBGROUPS:
        rows = subgroup_filter(mode_rows, subgroup)
        for predictor in PREDICTORS:
            fit = fit_map.get((subgroup, predictor), {})
            if predictor == "raw_parameter_baseline":
                raw_fit = fit_raw_baseline(rows)
                predictions = raw_predict(rows, raw_fit)
                for row, pred in zip(rows, predictions, strict=False):
                    y = finite_float(row, "delta_f")
                    pred_safe = isfinite(float(pred)) and float(pred) <= MAIN_THRESHOLD
                    true_safe = isfinite(y) and y <= MAIN_THRESHOLD
                    out.append(
                        {
                            "predictor": predictor,
                            "threshold_kind": "central",
                            "threshold_value": float("nan"),
                            "subgroup": subgroup,
                            "beta_deg": finite_float(row, "beta_deg"),
                            "mu": finite_float(row, "mu"),
                            "eta": finite_float(row, "eta"),
                            "true_safe": "true" if true_safe else "false",
                            "predicted_safe": "true" if pred_safe else "false",
                            "false_safe": "true" if pred_safe and not true_safe else "false",
                            "false_unsafe": "true" if true_safe and not pred_safe else "false",
                            "notes": "raw baseline predicted delta_f <= 0.10",
                        }
                    )
                continue
            for kind, _limit in THRESHOLD_KINDS:
                threshold_key = "threshold_10_central" if kind == "central" else f"threshold_{kind}"
                threshold = finite_float(fit, threshold_key)
                for row in rows:
                    x = finite_float(row, predictor)
                    y = finite_float(row, "delta_f")
                    if not (isfinite(x) and isfinite(y) and isfinite(threshold)):
                        continue
                    pred_safe = x <= threshold
                    true_safe = y <= MAIN_THRESHOLD
                    out.append(
                        {
                            "predictor": predictor,
                            "threshold_kind": kind,
                            "threshold_value": threshold,
                            "subgroup": subgroup,
                            "beta_deg": finite_float(row, "beta_deg"),
                            "mu": finite_float(row, "mu"),
                            "eta": finite_float(row, "eta"),
                            "true_safe": "true" if true_safe else "false",
                            "predicted_safe": "true" if pred_safe else "false",
                            "false_safe": "true" if pred_safe and not true_safe else "false",
                            "false_unsafe": "true" if true_safe and not pred_safe else "false",
                            "notes": "safe means predictor <= threshold",
                        }
                    )
    return out


def geometry_dependence_rows(
    mode_rows: Sequence[dict[str, object]],
    scalar_fits: dict[tuple[str, str], dict[str, float]],
) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    rows = subgroup_filter(mode_rows, "all_reliable")
    for predictor in PREDICTORS:
        if predictor == "raw_parameter_baseline":
            raw_fit = fit_raw_baseline(rows)
            pred_map = {
                id(row): pred for row, pred in zip(rows, raw_predict(rows, raw_fit), strict=False)
            }
        else:
            fit = scalar_fits.get(("all_reliable", predictor), {})
            pred_map = {id(row): scalar_prediction(row, fit, predictor) for row in rows}
        for group_key, group_name in (("beta_deg", "beta"), ("mu", "mu"), ("eta", "eta")):
            for value in sorted({finite_float(row, group_key) for row in rows if isfinite(finite_float(row, group_key))}):
                subset = [row for row in rows if abs(finite_float(row, group_key) - value) <= 1.0e-12]
                residuals = np.asarray(
                    [
                        finite_float(row, "delta_f") - pred_map.get(id(row), float("nan"))
                        for row in subset
                        if isfinite(pred_map.get(id(row), float("nan"))) and isfinite(finite_float(row, "delta_f"))
                    ],
                    dtype=float,
                )
                if residuals.size == 0:
                    continue
                if predictor == "raw_parameter_baseline":
                    cls = classification_metrics(subset, predictor=predictor, predictions=[pred_map.get(id(row), float("nan")) for row in subset])
                else:
                    fit = scalar_fits.get(("all_reliable", predictor), {})
                    cls = classification_metrics(subset, predictor=predictor, threshold=fit.get("threshold_10", float("nan")))
                out.append(
                    {
                        "predictor": predictor,
                        "grouping_variable": group_name,
                        "group_value": value,
                        "row_count": len(subset),
                        "mean_residual": float(np.mean(residuals)),
                        "median_residual": float(np.median(residuals)),
                        "max_abs_residual": float(np.max(np.abs(residuals))),
                        "p95_abs_residual": float(np.percentile(np.abs(residuals), 95.0)),
                        "false_safe_rate": cls["false_safe_rate"],
                        "false_unsafe_rate": cls["false_unsafe_rate"],
                        "notes": "residual is delta_f minus predictor-only predicted delta_f",
                    }
                )
    return out


def timo_normalized_det_counter(beta_deg: float, mu: float, eta: float, epsilon: float) -> tuple[object, dict[str, int]]:
    counter = {"n": 0}

    def func(lam: float) -> float:
        counter["n"] += 1
        return TIMO.timo_det_for_scan(float(lam), float(beta_deg), float(mu), float(epsilon), float(eta))

    return func, counter


def local_timo_root(
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    epsilon: float,
    lambda_eb: float,
    left: float,
    right: float,
) -> tuple[float, int, bool, str]:
    func, counter = timo_normalized_det_counter(beta_deg, mu, eta, epsilon)
    grid = np.linspace(max(beta_workflow.ROOT_SCAN_START, left), max(right, left + 1.0e-3), 80)
    values = np.asarray([func(float(value)) for value in grid], dtype=float)
    best_root = float("nan")
    failure = True
    notes = ""
    for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
        if not (isfinite(float(fa)) and isfinite(float(fb))):
            continue
        if fa == 0.0 or fa * fb < 0.0:
            try:
                root = float(brentq(func, float(a), float(b), xtol=1.0e-10, rtol=1.0e-10, maxiter=80))
            except ValueError:
                continue
            if not isfinite(best_root) or abs(root - float(lambda_eb)) < abs(best_root - float(lambda_eb)):
                best_root = root
                failure = False
                notes = "local sign bracket"
    if failure:
        finite_mask = np.isfinite(values)
        if np.any(finite_mask):
            best_index = int(np.argmin(np.abs(values[finite_mask])))
            finite_grid = grid[finite_mask]
            center = float(finite_grid[best_index])
            width = max((right - left) / 20.0, 0.02)
            a = max(float(left), center - width)
            b = min(float(right), center + width)
            try:
                res = minimize_scalar(lambda x: abs(func(float(x))), bounds=(a, b), method="bounded", options={"xatol": 1.0e-8})
                if res.success:
                    best_root = float(res.x)
                    failure = False
                    notes = "local smallest-singular determinant minimum"
            except (ValueError, FloatingPointError, OverflowError):
                pass
    return best_root, int(counter["n"]), failure, notes


def local_timo_benchmark_rows(
    args: Args,
    cache: RootCache,
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    eb_roots: Sequence[float],
    full_rows: Sequence[dict[str, object]],
    runtime: RuntimeTracker,
) -> list[dict[str, object]]:
    if not (
        any(abs(beta_deg - item) <= 1.0e-12 for item in (0.0, 45.0, 90.0))
        and any(abs(mu - item) <= 1.0e-12 for item in (0.0, 0.3, 0.7))
        and any(abs(eta - item) <= 1.0e-12 for item in (-0.1, 0.0, 0.1))
    ):
        return []
    rows: list[dict[str, object]] = []
    roots = [float(value) for value in eb_roots if isfinite(float(value))]
    for mode_index, lam in enumerate(roots[: args.n_reported_modes], start=1):
        left_neighbor = roots[mode_index - 2] if mode_index >= 2 else beta_workflow.ROOT_SCAN_START
        right_neighbor = roots[mode_index] if mode_index < len(roots) else lam + max(0.4, lam - left_neighbor)
        left = max(beta_workflow.ROOT_SCAN_START, lam - 0.45 * (lam - left_neighbor))
        right = lam + 0.45 * (right_neighbor - lam)
        start = time.perf_counter()
        local_root, evals, failure, notes = local_timo_root(
            beta_deg=beta_deg,
            mu=mu,
            eta=eta,
            epsilon=args.epsilon,
            lambda_eb=lam,
            left=left,
            right=right,
        )
        local_seconds = time.perf_counter() - start
        runtime.local_timo_seconds += local_seconds
        full = next((row for row in full_rows if int(row["eb_sorted_index"]) == mode_index), None)
        full_lam = finite_float(full, "Lambda_Timo") if full is not None else float("nan")
        full_delta = finite_float(full, "delta_f") if full is not None else float("nan")
        local_delta = frequency_metrics(lam, local_root)["delta_f"] if isfinite(local_root) else float("nan")
        rows.append(
            {
                "beta_deg": beta_deg,
                "mu": mu,
                "eta": eta,
                "mode": mode_index,
                "Lambda_EB": lam,
                "Lambda_Timo_full": full_lam,
                "Lambda_Timo_local": local_root,
                "delta_f_full": full_delta,
                "delta_f_local": local_delta,
                "local_timo_matrix_evaluations": evals,
                "full_timo_matrix_evaluations": "",
                "local_runtime": local_seconds,
                "full_runtime": "",
                "speedup": "",
                "local_delta_error": abs(local_delta - full_delta) if isfinite(local_delta) and isfinite(full_delta) else float("nan"),
                "local_root_failure": "true" if failure else "false",
                "cluster_warning": "true" if sorted_gap(roots, mode_index) <= args.cluster_gap_threshold else "false",
                "notes": notes,
            }
        )
    return rows


def run_scan(args: Args) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]], RuntimeTracker]:
    runtime = RuntimeTracker()
    cache = RootCache(args, runtime)
    mode_rows: list[dict[str, object]] = []
    matching_rows: list[dict[str, object]] = []
    point_rows: list[dict[str, object]] = []
    local_rows: list[dict[str, object]] = []
    total_points = len(args.beta_values) * len(args.mu_values) * len(args.eta_values)
    point_number = 0
    for beta_deg in args.beta_values:
        for mu in args.mu_values:
            for eta in args.eta_values:
                point_number += 1
                print(
                    f"[{point_number}/{total_points}] beta={beta_deg:g}, mu={mu:g}, eta={eta:g}",
                    flush=True,
                )
                products = point_products(args, cache, runtime, beta_deg=float(beta_deg), mu=float(mu), eta=float(eta))
                mode_rows.extend(products.mode_rows)
                matching_rows.extend(products.matching_rows)
                point_rows.append(products.summary_row)
                local_rows.extend(products.local_timo_rows)
    return mode_rows, matching_rows, point_rows, local_rows, runtime


def runtime_rows(runtime: RuntimeTracker, *, point_count: int, mode_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    mode_count = sum(1 for row in mode_rows if row.get("comparison_type") == "homologous_mode")
    full_reference = max(runtime.eb_roots_seconds + runtime.timo_roots_seconds + runtime.matching_seconds, 1.0e-30)

    def row(stage: str, seconds: float, modes: int, notes: str) -> dict[str, object]:
        return {
            "stage": stage,
            "parameter_point_count": point_count,
            "mode_count": modes,
            "total_seconds": seconds,
            "mean_seconds_per_point": safe_ratio(seconds, point_count),
            "mean_seconds_per_mode": safe_ratio(seconds, modes),
            "matrix_evaluation_count": "",
            "relative_cost_vs_full_reference": safe_ratio(seconds, full_reference),
            "notes": notes,
        }

    return [
        row("EB roots only", runtime.eb_roots_seconds, point_count, f"root calls={runtime.eb_root_calls}"),
        row("EB roots + Theta_max_EB", runtime.eb_roots_seconds + runtime.theta_seconds, mode_count, "Theta uses EB Lambda only"),
        row("EB roots + EB shape + Pi_EB", runtime.eb_roots_seconds + runtime.eb_shape_pi_seconds, mode_count, "Pi uses EB shape reconstruction"),
        row("Full Timoshenko roots", runtime.timo_roots_seconds, point_count, f"root calls={runtime.timo_root_calls}"),
        row("Cross-model matching", runtime.matching_seconds, mode_count, "Hungarian assignment and mass-weighted MAC"),
        row("Optional local Timoshenko correction", runtime.local_timo_seconds, mode_count, "only populated when --benchmark-local-timo is enabled"),
    ]


def plot_n10_heatmaps(output_dir: Path, point_rows: Sequence[dict[str, object]]) -> list[Path]:
    paths: list[Path] = []
    target_counts = [
        int(finite_float(row, "target_mode_count"))
        for row in point_rows
        if isfinite(finite_float(row, "target_mode_count"))
    ]
    color_max = max(target_counts, default=DEFAULT_N_REPORTED_MODES)
    for eta in sorted({finite_float(row, "eta") for row in point_rows}):
        subset = [row for row in point_rows if abs(finite_float(row, "eta") - eta) <= 1.0e-12]
        betas = sorted({finite_float(row, "beta_deg") for row in subset})
        mus = sorted({finite_float(row, "mu") for row in subset})
        data = np.full((len(mus), len(betas)), np.nan, dtype=float)
        for row in subset:
            i = mus.index(finite_float(row, "mu"))
            j = betas.index(finite_float(row, "beta_deg"))
            data[i, j] = finite_float(row, "N10_sorted_consecutive")
        fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
        im = ax.imshow(data, origin="lower", aspect="auto", vmin=0, vmax=color_max)
        ax.set_xticks(range(len(betas)))
        ax.set_xticklabels([f"{value:g}" for value in betas])
        ax.set_yticks(range(len(mus)))
        ax.set_yticklabels([f"{value:g}" for value in mus])
        ax.set_xlabel("beta, deg")
        ax.set_ylabel("mu")
        ax.set_title(f"N10 sorted, eta={eta:g}")
        fig.colorbar(im, ax=ax, label="N10")
        token = "m0p1" if eta < 0 else ("p0p1" if eta > 0 else "0")
        path = output_dir / f"N10_sorted_heatmap_eta_{token}.png"
        fig.savefig(path, dpi=180)
        plt.close(fig)
        paths.append(path)
    return paths


def plot_predictor_collapse(output_dir: Path, mode_rows: Sequence[dict[str, object]], fit_rows: Sequence[dict[str, object]]) -> Path:
    rows = subgroup_filter(mode_rows, "all_reliable")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), constrained_layout=True)
    for ax, predictor in zip(axes, ("Theta_max_EB", "Pi_EB"), strict=True):
        for character, color in (
            ("bending_dominated", "#1f77b4"),
            ("mixed", "#d55e00"),
            ("longitudinal_dominated", "#009e73"),
        ):
            subset = [row for row in rows if row.get("Timo_classification") == character]
            ax.scatter(
                [finite_float(row, predictor) for row in subset],
                [finite_float(row, "delta_f") for row in subset],
                s=18,
                alpha=0.6,
                label=character.replace("_", " "),
                color=color,
            )
        fit = next((row for row in fit_rows if row["subgroup"] == "all_reliable" and row["predictor"] == predictor), None)
        if fit is not None and isfinite(finite_float(fit, "threshold_10_central")):
            ax.axvline(finite_float(fit, "threshold_10_central"), color="black", ls=":", lw=1.0)
        ax.axhline(MAIN_THRESHOLD, color="black", ls="--", lw=1.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(predictor)
        ax.set_ylabel("delta_f")
        ax.grid(True, which="both", alpha=0.25)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labels, fontsize=8)
    path = output_dir / "predictor_collapse_Theta_vs_Pi.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_residual_beta(output_dir: Path, mode_rows: Sequence[dict[str, object]], scalar_fits: dict[tuple[str, str], dict[str, float]]) -> Path:
    rows = subgroup_filter(mode_rows, "all_reliable")
    theta_fit = scalar_fits.get(("all_reliable", "Theta_max_EB"), {})
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8), constrained_layout=True)
    betas = [finite_float(row, "beta_deg") for row in rows]
    theta_res = [finite_float(row, "delta_f") - scalar_prediction(row, theta_fit, "Theta_max_EB") for row in rows]
    pi_res = [finite_float(row, "residual_Pi_half") for row in rows]
    axes[0].scatter(betas, theta_res, s=16, alpha=0.55)
    axes[1].scatter(betas, pi_res, s=16, alpha=0.55, color="#d55e00")
    axes[0].axhline(0.0, color="black", lw=1.0)
    axes[1].axhline(0.0, color="black", lw=1.0)
    axes[0].set_title("Theta fit residual")
    axes[1].set_title("delta_f - 0.5 Pi_EB")
    for ax in axes:
        ax.set_xlabel("beta, deg")
        ax.set_ylabel("residual")
        ax.grid(True, alpha=0.25)
    path = output_dir / "predictor_residual_vs_beta.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_threshold_beta_eta(output_dir: Path, threshold_rows: Sequence[dict[str, object]]) -> Path:
    rows = [
        row for row in threshold_rows
        if row["subgroup"] == "all_reliable" and row["threshold_kind"] == "central" and row["predictor"] in {"Theta_max_EB", "Pi_EB"}
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), constrained_layout=True)
    for ax, predictor in zip(axes, ("Theta_max_EB", "Pi_EB"), strict=True):
        subset = [row for row in rows if row["predictor"] == predictor]
        for eta in sorted({finite_float(row, "eta") for row in subset}):
            eta_rows = [row for row in subset if abs(finite_float(row, "eta") - eta) <= 1.0e-12]
            grouped: dict[float, list[int]] = {}
            for row in eta_rows:
                grouped.setdefault(finite_float(row, "beta_deg"), []).append(1 if str(row["false_safe"]).lower() == "true" else 0)
            xs = sorted(grouped)
            ys = [float(np.mean(grouped[x])) for x in xs]
            ax.plot(xs, ys, marker="o", label=f"eta={eta:g}")
        ax.set_xlabel("beta, deg")
        ax.set_ylabel("false-safe fraction at central threshold")
        ax.set_title(predictor)
        ax.grid(True, alpha=0.25)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(handles, labels, fontsize=8)
    path = output_dir / "predictor_threshold_by_beta_eta.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_false_safe_coverage(output_dir: Path, fit_rows: Sequence[dict[str, object]]) -> Path:
    rows = [row for row in fit_rows if row["subgroup"] == "all_reliable" and row["predictor"] in SCALAR_PREDICTORS]
    x = np.arange(len(rows), dtype=float)
    false_safe = [finite_float(row, "false_safe_rate_central") for row in rows]
    coverage = [finite_float(row, "safe_coverage_central") for row in rows]
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    width = 0.38
    ax.bar(x - width / 2.0, false_safe, width=width, label="false-safe")
    ax.bar(x + width / 2.0, coverage, width=width, label="safe coverage")
    ax.set_xticks(x)
    ax.set_xticklabels([str(row["predictor"]) for row in rows], rotation=25, ha="right")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    path = output_dir / "predictor_false_safe_coverage.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_delta_beta_selected(
    output_dir: Path,
    mode_rows: Sequence[dict[str, object]],
    *,
    n_reported_modes: int,
) -> Path:
    rows = [
        row for row in mode_rows
        if row.get("comparison_type") == "homologous_mode"
        and int(finite_float(row, "eb_sorted_index")) <= int(n_reported_modes)
        and any(abs(finite_float(row, "mu") - value) <= 1.0e-12 for value in (0.0, 0.3, 0.7))
        and any(abs(finite_float(row, "eta") - value) <= 1.0e-12 for value in (-0.1, 0.0, 0.1))
    ]
    fig, axes = plt.subplots(3, 3, figsize=(12.0, 9.0), sharex=True, sharey=True, constrained_layout=True)
    for i, mu in enumerate((0.0, 0.3, 0.7)):
        for j, eta in enumerate((-0.1, 0.0, 0.1)):
            ax = axes[i, j]
            subset = [row for row in rows if abs(finite_float(row, "mu") - mu) <= 1.0e-12 and abs(finite_float(row, "eta") - eta) <= 1.0e-12]
            for mode in range(1, int(n_reported_modes) + 1):
                mrows = sorted([row for row in subset if int(finite_float(row, "eb_sorted_index")) == mode], key=lambda row: finite_float(row, "beta_deg"))
                if mrows:
                    ax.plot([finite_float(row, "beta_deg") for row in mrows], [finite_float(row, "delta_f") for row in mrows], lw=1.0)
            ax.axhline(MAIN_THRESHOLD, color="black", ls="--", lw=0.8)
            ax.set_title(f"mu={mu:g}, eta={eta:g}", fontsize=9)
            ax.grid(True, alpha=0.2)
    for ax in axes[-1, :]:
        ax.set_xlabel("beta, deg")
    for ax in axes[:, 0]:
        ax.set_ylabel("delta_f")
    path = output_dir / "delta_f_vs_beta_selected_mu_eta.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_reduced_vs_raw_cv(output_dir: Path, cv_rows: Sequence[dict[str, object]]) -> Path:
    rows = [row for row in cv_rows if row["subgroup"] == "all_reliable" and row["split_kind"] in {"leave_one_beta_out", "leave_one_eta_out"}]
    predictors = ["Theta_max_EB", "Pi_EB", "raw_parameter_baseline"]
    splits = ["leave_one_beta_out", "leave_one_eta_out"]
    data = np.full((len(splits), len(predictors)), np.nan)
    for i, split in enumerate(splits):
        for j, predictor in enumerate(predictors):
            values = [finite_float(row, "test_RMSE") for row in rows if row["split_kind"] == split and row["predictor"] == predictor and isfinite(finite_float(row, "test_RMSE"))]
            if values:
                data[i, j] = float(np.mean(values))
    fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=True)
    x = np.arange(len(predictors))
    width = 0.36
    ax.bar(x - width / 2.0, data[0], width=width, label="leave beta")
    ax.bar(x + width / 2.0, data[1], width=width, label="leave eta")
    ax.set_xticks(x)
    ax.set_xticklabels(predictors, rotation=20, ha="right")
    ax.set_ylabel("mean test RMSE")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    path = output_dir / "reduced_vs_raw_parameter_model_cv.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def generate_plots(
    output_dir: Path,
    *,
    mode_rows: Sequence[dict[str, object]],
    point_rows: Sequence[dict[str, object]],
    fit_rows: Sequence[dict[str, object]],
    scalar_fits: dict[tuple[str, str], dict[str, float]],
    threshold_rows: Sequence[dict[str, object]],
    cv_rows: Sequence[dict[str, object]],
    n_reported_modes: int,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_mode_counts = [
        int(finite_float(row, "target_mode_count"))
        for row in point_rows
        if isfinite(finite_float(row, "target_mode_count"))
    ]
    plot_mode_count = max(summary_mode_counts, default=int(n_reported_modes))
    paths = plot_n10_heatmaps(output_dir, point_rows)
    paths.extend(
        [
            plot_predictor_collapse(output_dir, mode_rows, fit_rows),
            plot_residual_beta(output_dir, mode_rows, scalar_fits),
            plot_threshold_beta_eta(output_dir, threshold_rows),
            plot_false_safe_coverage(output_dir, fit_rows),
            plot_delta_beta_selected(
                output_dir,
                mode_rows,
                n_reported_modes=plot_mode_count,
            ),
            plot_reduced_vs_raw_cv(output_dir, cv_rows),
        ]
    )
    return paths


def best_cv_row(cv_rows: Sequence[dict[str, object]], *, split: str, metric: str, minimize: bool = True) -> dict[str, object] | None:
    rows = [
        row for row in cv_rows
        if row["subgroup"] == "all_reliable" and row["split_kind"] == split and isfinite(finite_float(row, metric))
    ]
    if not rows:
        return None
    grouped: list[dict[str, object]] = []
    for predictor in sorted({str(row["predictor"]) for row in rows}):
        values = [finite_float(row, metric) for row in rows if row["predictor"] == predictor and isfinite(finite_float(row, metric))]
        if values:
            grouped.append({"predictor": predictor, metric: float(np.mean(values))})
    if not grouped:
        return None
    return sorted(grouped, key=lambda row: finite_float(row, metric), reverse=not minimize)[0]


def report_table_thresholds(fit_rows: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| predictor | central | safe 0% | safe 1% | safe 5% | false-safe central | safe coverage central |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for predictor in ("Theta_max_EB", "Pi_EB"):
        row = next((item for item in fit_rows if item["subgroup"] == "all_reliable" and item["predictor"] == predictor), None)
        if row is None:
            continue
        lines.append(
            "| "
            f"{predictor} | "
            f"{finite_float(row, 'threshold_10_central'):.6g} | "
            f"{finite_float(row, 'threshold_safe_0pct'):.6g} | "
            f"{finite_float(row, 'threshold_safe_1pct'):.6g} | "
            f"{finite_float(row, 'threshold_safe_5pct'):.6g} | "
            f"{finite_float(row, 'false_safe_rate_central'):.4g} | "
            f"{finite_float(row, 'safe_coverage_central'):.4g} |"
        )
    return lines


def write_report(
    output_dir: Path,
    *,
    args: Args,
    mode_rows: Sequence[dict[str, object]],
    point_rows: Sequence[dict[str, object]],
    matching_rows: Sequence[dict[str, object]],
    fit_rows: Sequence[dict[str, object]],
    cv_rows: Sequence[dict[str, object]],
    runtime: RuntimeTracker,
    plot_paths: Sequence[Path],
    local_timo_rows: Sequence[dict[str, object]],
) -> Path:
    path = output_dir / "eb_validity_fixed_epsilon_geometry_scan_report.md"
    reliable = subgroup_filter(mode_rows, "all_reliable")
    beta_groups = {
        beta: [row for row in reliable if abs(finite_float(row, "beta_deg") - beta) <= 1.0e-12]
        for beta in sorted({finite_float(row, "beta_deg") for row in reliable})
    }
    beta_mean_delta = {beta: float(np.mean([finite_float(row, "delta_f") for row in rows])) for beta, rows in beta_groups.items() if rows}
    beta_effect_span = max(beta_mean_delta.values()) - min(beta_mean_delta.values()) if beta_mean_delta else float("nan")
    theta_row = next((row for row in fit_rows if row["subgroup"] == "all_reliable" and row["predictor"] == "Theta_max_EB"), {})
    pi_row = next((row for row in fit_rows if row["subgroup"] == "all_reliable" and row["predictor"] == "Pi_EB"), {})
    raw_row = next((row for row in fit_rows if row["subgroup"] == "all_reliable" and row["predictor"] == "raw_parameter_baseline"), {})
    best_beta_rmse = best_cv_row(cv_rows, split="leave_one_beta_out", metric="test_RMSE", minimize=True)
    best_eta_rmse = best_cv_row(cv_rows, split="leave_one_eta_out", metric="test_RMSE", minimize=True)
    lowest_fs = min(
        (row for row in fit_rows if row["subgroup"] == "all_reliable" and isfinite(finite_float(row, "false_safe_rate_central"))),
        key=lambda row: finite_float(row, "false_safe_rate_central"),
        default=None,
    )
    highest_cov = max(
        (row for row in fit_rows if row["subgroup"] == "all_reliable" and isfinite(finite_float(row, "safe_coverage_central"))),
        key=lambda row: finite_float(row, "safe_coverage_central"),
        default=None,
    )
    cluster_rows = [row for row in mode_rows if row.get("comparison_type") == "matched_cluster"]
    lines = [
        "# Fixed-Epsilon EB Validity Geometry Scan",
        "",
        "## Scope",
        "",
        "This diagnostic compares Euler--Bernoulli and Timoshenko in-plane beam theories at fixed `epsilon_0`. Timoshenko is used as a beam-theory reference model for calibration, not as exact 3D elasticity.",
        "",
        "No FEM, 3D FEM, Gmsh, CalculiX, article workspace, `main.tex`, or article figures are used.",
        "",
        "## Parameter Grid",
        "",
        f"- epsilon_0: `{args.epsilon:g}`",
        f"- beta values: `{', '.join(f'{value:g}' for value in args.beta_values)}`",
        f"- mu values: `{', '.join(f'{value:g}' for value in args.mu_values)}`",
        f"- eta values: `{', '.join(f'{value:g}' for value in args.eta_values)}`",
        f"- reported modes: `{args.n_reported_modes}`",
        f"- candidate roots: `{args.n_candidate_roots}`; boundary retry roots: `{BOUNDARY_RETRY_ROOTS}`",
        f"- shape points: `{args.n_shape_points}`",
        "",
        "## Definitions",
        "",
        "- `delta_f = abs(Lambda_EB^2 - Lambda_Timo^2) / Lambda_Timo^2`; the square is used because dimensional frequency has the common factor form `f=C*Lambda^2`.",
        "- Sorted comparison keeps current sorted index only.",
        "- Homologous-mode comparison is independent at each geometry point and uses mass-weighted local `[u,w]` MAC plus Hungarian assignment.",
        "- `Theta_max_EB=(1+E/(k'G))*chi_max_EB^2` uses EB frequency and geometry only.",
        "- `Pi_EB=Pi_shear_EB+Pi_rotary_EB` uses EB mode shape, analytic EB derivatives including `w'''`, and the current repository `k'`.",
        "",
        "## Data Quality",
        "",
        f"- Mode rows: `{len(mode_rows)}`",
        f"- Homologous reliable individual rows: `{len(reliable)}`",
        f"- Cluster rows: `{len(cluster_rows)}`",
        f"- Matching audit rows: `{len(matching_rows)}`",
        f"- Candidate-boundary retries: `{runtime.boundary_retries}`",
        f"- Root cache hits/misses: `{runtime.root_cache_hits}/{runtime.root_cache_misses}`",
        "",
        "## Thresholds",
        "",
        *report_table_thresholds(fit_rows),
        "",
        "## Answers",
        "",
        f"1. Beta effect at fixed epsilon: the mean reliable-row `delta_f` span over beta groups is `{beta_effect_span:.6g}`. See `fixed_epsilon_geometry_dependence.csv` for residual group statistics.",
        f"2. Pi_EB absorption: predictor-only R2 is `{finite_float(pi_row, 'R2'):.4g}` and geometry-extended R2 is `{finite_float(pi_row, 'geometry_extended_R2'):.4g}`.",
        f"3. Theta_max_EB conservative value: central false-safe rate is `{finite_float(theta_row, 'false_safe_rate_central'):.4g}`; zero/one/five-percent thresholds are reported above.",
        f"4. Lowest all-reliable false-safe predictor: `{lowest_fs.get('predictor') if lowest_fs else 'n/a'}`. Highest safe coverage: `{highest_cov.get('predictor') if highest_cov else 'n/a'}`.",
        f"5. Best leave-one-beta-out transfer by mean RMSE: `{best_beta_rmse.get('predictor') if best_beta_rmse else 'n/a'}`.",
        f"6. Best leave-one-eta-out transfer by mean RMSE: `{best_eta_rmse.get('predictor') if best_eta_rmse else 'n/a'}`.",
        "7. Longitudinal and mixed modes are kept as separate subgroups; ambiguous or cluster-member individual matches are excluded from individual calibration.",
        "8. Beta=0 and beta=90 are treated with close-cluster diagnostics; ordinary low MAC inside a close cluster is not by itself treated as a failed theory comparison.",
        f"9. One scalar predictor sufficiency: compare predictor-only R2 with geometry-extended R2. For Pi_EB the improvement is `{finite_float(pi_row, 'geometry_extended_R2') - finite_float(pi_row, 'R2'):.4g}` when both are finite.",
        f"10. Raw-parameter baseline all-reliable R2 is `{finite_float(raw_row, 'R2'):.4g}`; it is a diagnostic baseline, not a physical law.",
        f"11. Criterion cost: `Pi_EB` relative cost is saved in `fixed_epsilon_runtime_costs.csv`; EB roots plus Pi cost was `{runtime.eb_roots_seconds + runtime.eb_shape_pi_seconds:.4g}` s in this run.",
        "12. The fixed-epsilon scan cannot establish the final dependence on epsilon_0. It only tests transfer of the boundary near epsilon_0=0.02 across beta, mu, and eta.",
    ]
    if local_timo_rows:
        failures = sum(1 for row in local_timo_rows if str(row.get("local_root_failure", "")).lower() == "true")
        errors = [finite_float(row, "local_delta_error") for row in local_timo_rows if isfinite(finite_float(row, "local_delta_error"))]
        lines.extend(
            [
                "",
                "## Local Timoshenko Benchmark",
                "",
                f"- rows: `{len(local_timo_rows)}`",
                f"- local root failures: `{failures}`",
                f"- median absolute delta_f error: `{float(np.median(errors)) if errors else float('nan'):.6g}`",
            ]
        )
    else:
        lines.extend(["", "## Local Timoshenko Benchmark", "", "Not run; `--benchmark-local-timo` was not enabled."])
    lines.extend(
        [
            "",
            "## Plots",
            "",
        ]
    )
    for plot in plot_paths:
        lines.append(f"- `{rel(plot)}`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_outputs(
    args: Args,
    *,
    mode_rows: list[dict[str, object]],
    matching_rows: list[dict[str, object]],
    point_rows: list[dict[str, object]],
    local_rows: list[dict[str, object]],
    runtime: RuntimeTracker,
    elapsed: float,
) -> tuple[list[Path], Path, list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    fit_rows, scalar_fits = fit_summary_rows(mode_rows)
    cv_rows = cross_validation_rows(mode_rows)
    geometry_rows = geometry_dependence_rows(mode_rows, scalar_fits)
    threshold_rows = threshold_classification_rows(mode_rows, fit_rows)
    runtime_table = runtime_rows(
        runtime,
        point_count=len(args.beta_values) * len(args.mu_values) * len(args.eta_values),
        mode_rows=mode_rows,
    )
    runtime_table.append(
        {
            "stage": "Total wall time",
            "parameter_point_count": len(args.beta_values) * len(args.mu_values) * len(args.eta_values),
            "mode_count": len(mode_rows),
            "total_seconds": elapsed,
            "mean_seconds_per_point": safe_ratio(elapsed, len(args.beta_values) * len(args.mu_values) * len(args.eta_values)),
            "mean_seconds_per_mode": safe_ratio(elapsed, len(mode_rows)),
            "matrix_evaluation_count": "",
            "relative_cost_vs_full_reference": 1.0,
            "notes": f"timestamp={datetime.now().isoformat(timespec='seconds')}; workers requested={args.workers}; serial execution",
        }
    )
    paths = [
        write_csv(args.output_dir / "fixed_epsilon_mode_metrics.csv", mode_rows, MODE_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_matching_audit.csv", matching_rows, MATCHING_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_point_summary.csv", point_rows, POINT_SUMMARY_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_predictor_fit_summary.csv", fit_rows, FIT_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_predictor_cross_validation.csv", cv_rows, CV_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_geometry_dependence.csv", geometry_rows, GEOMETRY_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_threshold_classification.csv", threshold_rows, THRESHOLD_CLASS_FIELDS),
        write_csv(args.output_dir / "fixed_epsilon_runtime_costs.csv", runtime_table, RUNTIME_FIELDS),
    ]
    if args.benchmark_local_timo:
        paths.append(write_csv(args.output_dir / "local_timo_correction_benchmark.csv", local_rows, LOCAL_TIMO_FIELDS))
    plot_paths = generate_plots(
        args.output_dir,
        mode_rows=mode_rows,
        point_rows=point_rows,
        fit_rows=fit_rows,
        scalar_fits=scalar_fits,
        threshold_rows=threshold_rows,
        cv_rows=cv_rows,
        n_reported_modes=args.n_reported_modes,
    )
    report = write_report(
        args.output_dir,
        args=args,
        mode_rows=mode_rows,
        point_rows=point_rows,
        matching_rows=matching_rows,
        fit_rows=fit_rows,
        cv_rows=cv_rows,
        runtime=runtime,
        plot_paths=plot_paths,
        local_timo_rows=local_rows,
    )
    return [*paths, *plot_paths], report, fit_rows, cv_rows, threshold_rows


def plot_only(args: Args) -> tuple[list[Path], Path]:
    mode_rows = [dict(row) for row in read_csv(args.output_dir / "fixed_epsilon_mode_metrics.csv")]
    point_rows = [dict(row) for row in read_csv(args.output_dir / "fixed_epsilon_point_summary.csv")]
    matching_rows = [dict(row) for row in read_csv(args.output_dir / "fixed_epsilon_matching_audit.csv")]
    fit_rows = [dict(row) for row in read_csv(args.output_dir / "fixed_epsilon_predictor_fit_summary.csv")]
    cv_rows = [dict(row) for row in read_csv(args.output_dir / "fixed_epsilon_predictor_cross_validation.csv")]
    threshold_rows = [dict(row) for row in read_csv(args.output_dir / "fixed_epsilon_threshold_classification.csv")]
    _fit, scalar_fits = fit_summary_rows([dict(row) for row in mode_rows])
    plot_paths = generate_plots(
        args.output_dir,
        mode_rows=mode_rows,
        point_rows=point_rows,
        fit_rows=fit_rows,
        scalar_fits=scalar_fits,
        threshold_rows=threshold_rows,
        cv_rows=cv_rows,
        n_reported_modes=args.n_reported_modes,
    )
    runtime = runtime_from_cost_csv(args.output_dir / "fixed_epsilon_runtime_costs.csv", args.cache_dir)
    local_path = args.output_dir / "local_timo_correction_benchmark.csv"
    local_rows = read_csv(local_path) if local_path.exists() else []
    report = write_report(
        args.output_dir,
        args=args,
        mode_rows=mode_rows,
        point_rows=point_rows,
        matching_rows=matching_rows,
        fit_rows=fit_rows,
        cv_rows=cv_rows,
        runtime=runtime,
        plot_paths=plot_paths,
        local_timo_rows=local_rows,
    )
    return plot_paths, report


def runtime_from_cost_csv(path: Path, cache_dir: Path | None = None) -> RuntimeTracker:
    runtime = RuntimeTracker()
    if not path.exists():
        return runtime
    rows = read_csv(path)
    by_stage = {row.get("stage", ""): row for row in rows}
    eb = finite_float(by_stage.get("EB roots only", {}), "total_seconds")
    eb_plus_pi = finite_float(by_stage.get("EB roots + EB shape + Pi_EB", {}), "total_seconds")
    timo = finite_float(by_stage.get("Full Timoshenko roots", {}), "total_seconds")
    matching = finite_float(by_stage.get("Cross-model matching", {}), "total_seconds")
    local = finite_float(by_stage.get("Optional local Timoshenko correction", {}), "total_seconds")
    runtime.eb_roots_seconds = eb if isfinite(eb) else 0.0
    runtime.eb_shape_pi_seconds = max(eb_plus_pi - runtime.eb_roots_seconds, 0.0) if isfinite(eb_plus_pi) else 0.0
    runtime.timo_roots_seconds = timo if isfinite(timo) else 0.0
    runtime.matching_seconds = matching if isfinite(matching) else 0.0
    runtime.local_timo_seconds = local if isfinite(local) else 0.0
    if cache_dir is not None and cache_dir.exists():
        runtime.root_cache_hits = sum(1 for path_obj in cache_dir.glob("root_*.npz"))
        runtime.root_cache_misses = 0
    return runtime


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    start = time.perf_counter()
    if args.plot_only:
        paths, report = plot_only(args)
        elapsed = time.perf_counter() - start
        print("plot-only regeneration complete:")
        for path in [report, *paths]:
            print(f"  {rel(path)}")
        print(f"plot-only seconds: {elapsed:.3f}")
        return {"paths": paths, "report": report, "seconds": elapsed}
    mode_rows, matching_rows, point_rows, local_rows, runtime = run_scan(args)
    elapsed = time.perf_counter() - start
    paths, report, fit_rows, cv_rows, threshold_rows = write_outputs(
        args,
        mode_rows=mode_rows,
        matching_rows=matching_rows,
        point_rows=point_rows,
        local_rows=local_rows,
        runtime=runtime,
        elapsed=elapsed,
    )
    print("generated fixed-epsilon EB validity geometry scan outputs:")
    for path in [report, *paths]:
        print(f"  {rel(path)}")
    print(f"parameter points: {len(point_rows)}")
    print(f"mode rows: {len(mode_rows)}")
    print(f"cache hits/misses: {runtime.root_cache_hits}/{runtime.root_cache_misses}")
    print(f"candidate-boundary retries: {runtime.boundary_retries}")
    print(f"elapsed seconds: {elapsed:.3f}")
    return {
        "paths": paths,
        "report": report,
        "mode_rows": mode_rows,
        "matching_rows": matching_rows,
        "point_rows": point_rows,
        "fit_rows": fit_rows,
        "cv_rows": cv_rows,
        "threshold_rows": threshold_rows,
        "runtime": runtime,
        "seconds": elapsed,
    }


if __name__ == "__main__":
    main()
