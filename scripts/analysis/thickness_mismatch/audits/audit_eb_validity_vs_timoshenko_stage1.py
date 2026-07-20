from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from datetime import datetime
import hashlib
import json
from math import isfinite
from pathlib import Path
import sys
import time
from typing import Callable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    thickness_mismatch_factors,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_timoshenko_shape_construction as shape_audit,
)
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow,
)
from scripts.lib import in_plane_shape_geometry as DISPLAY  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


CACHE_VERSION = "eb_validity_vs_timoshenko_stage1_v1"
DEFAULT_OUTPUT_DIR = Path("results") / "eb_validity_vs_timoshenko_stage1"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_validity_vs_timoshenko_stage1"
DEFAULT_CACHE_SUBDIR = "cache"

DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
DEFAULT_MU_VALUES = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
DEFAULT_EPSILON_VALUES = (
    0.0025,
    0.005,
    0.010,
    0.015,
    0.020,
    0.025,
    0.030,
    0.035,
    0.040,
    0.045,
    0.050,
    0.055,
    0.060,
)
SMOKE_MU_VALUES = (0.0, 0.7)
SMOKE_EPSILON_VALUES = (0.0025, 0.03, 0.06)
DEFAULT_N_REPORTED_MODES = 8
DEFAULT_N_CANDIDATE_ROOTS = 12
SMOKE_N_REPORTED_MODES = 4
SMOKE_N_CANDIDATE_ROOTS = 6
DEFAULT_SHAPE_SAMPLE_COUNT = 401
SMOKE_SHAPE_SAMPLE_COUNT = 201
THRESHOLDS = (0.05, 0.10, 0.15)
MAIN_THRESHOLD = 0.10
MAC_WARNING_THRESHOLD = 0.80
MAC_MARGIN_WARNING_THRESHOLD = 0.05
CRITICAL_EPS_TOL = 1.0e-4
CRITICAL_DELTA_TOL = 1.0e-3
DEFAULT_MAX_REFINEMENT_ITERATIONS = 15

MODEL_EB = shape_audit.MODEL_EB
MODEL_TIMO = shape_audit.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)

MODE_LEVEL_FIELDS = [
    "beta_deg",
    "mu",
    "eta",
    "epsilon_0",
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
    "branch_id",
    "eb_sorted_index",
    "timo_sorted_index",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_Lambda",
    "delta_f",
    "delta_f_symmetric",
    "bias_f",
    "EB_axial_fraction",
    "EB_bending_fraction",
    "EB_classification",
    "Timo_axial_fraction",
    "Timo_bending_fraction",
    "Timo_shear_fraction",
    "Timo_bending_shear_fraction",
    "Timo_classification",
    "chi_1",
    "chi_2",
    "chi_max",
    "chi_eff",
    "pass_5pct",
    "pass_10pct",
    "pass_15pct",
    "shape_mac_EB_Timo",
    "tracking_warning",
    "notes",
]

SUMMARY_FIELDS = [
    "beta_deg",
    "mu",
    "eta",
    "epsilon_0",
    "epsilon_1",
    "epsilon_2",
    "epsilon_max",
    "target_mode_count",
    "N5_sorted_consecutive",
    "N10_sorted_consecutive",
    "N15_sorted_consecutive",
    "n_pass_5_among_first8",
    "n_pass_10_among_first8",
    "n_pass_15_among_first8",
    "n_pass_5_among_firstK",
    "n_pass_10_among_firstK",
    "n_pass_15_among_firstK",
    "max_delta_f_firstK",
    "mean_delta_f_firstK",
    "pass_pattern_5",
    "pass_pattern_10",
    "pass_pattern_15",
    "N5_bending_consecutive",
    "N10_bending_consecutive",
    "N15_bending_consecutive",
    "n_bending_pass_10_among_first8",
    "n_bending_pass_10_among_firstK",
    "first_failed_sorted_index_10",
    "first_failed_branch_id_10",
    "first_failed_delta_f_10",
    "first_failed_character_10",
    "warnings",
    "notes",
]

CRITICAL_FIELDS = [
    "beta_deg",
    "mu",
    "eta",
    "branch_id",
    "base_sorted_index",
    "epsilon_0_crit_10",
    "critical_status",
    "crossing_count",
    "nonmonotonic",
    "epsilon_low_bracket",
    "epsilon_high_bracket",
    "delta_low",
    "delta_high",
    "refinement_iterations",
    "epsilon_1_at_crit",
    "epsilon_2_at_crit",
    "epsilon_max_at_crit",
    "r0_over_l0_at_crit",
    "d0_over_l0_at_crit",
    "l1_over_r1_at_crit",
    "l2_over_r2_at_crit",
    "l1_over_d1_at_crit",
    "l2_over_d2_at_crit",
    "chi_max_at_crit",
    "chi_eff_at_crit",
    "EB_sorted_index_at_crit",
    "Timo_sorted_index_at_crit",
    "mode_character_at_crit",
    "axial_fraction_at_crit",
    "bending_fraction_at_crit",
    "shear_fraction_at_crit",
    "warning",
    "notes",
]

TRACKING_AUDIT_FIELDS = [
    "model",
    "mu",
    "eta",
    "epsilon_prev",
    "epsilon_current",
    "branch_id",
    "previous_sorted_index",
    "current_sorted_index",
    "Lambda",
    "MAC_to_previous",
    "second_best_MAC",
    "MAC_margin",
    "ambiguous",
    "warning",
    "notes",
]

THRESHOLD_AUDIT_FIELDS = [
    "mu",
    "branch_id",
    "iteration",
    "epsilon_low",
    "epsilon_high",
    "epsilon_trial",
    "delta_low",
    "delta_high",
    "delta_trial",
    "EB_sorted_index",
    "Timo_sorted_index",
    "shape_MAC",
    "status",
    "notes",
]

CHI_FIT_FIELDS = [
    "metric",
    "point_count",
    "C",
    "p",
    "R2",
    "chi_crit_for_delta_f_0p10",
    "notes",
]

TIMING_FIELDS = [
    "run_kind",
    "seconds",
    "cache_hits",
    "cache_misses",
    "solved_parameter_points",
    "threshold_refinement_iterations",
    "notes",
]


@dataclass(frozen=True)
class Args:
    beta_deg: float
    eta: float
    mu_values: tuple[float, ...]
    epsilon_values: tuple[float, ...]
    n_reported_modes: int
    n_candidate_roots: int
    n_shape_points: int
    output_dir: Path
    cache_dir: Path
    reuse_cache: bool
    force_recompute: bool
    plot_only: bool
    skip_critical_refinement: bool
    smoke: bool
    max_refinement_iterations: int


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
class BranchState:
    branch_id: str
    base_sorted_index: int
    current_sorted_index: int
    Lambda: float
    result: shape_audit.ModeResult
    vector: np.ndarray
    mac_to_previous: float
    second_best_mac: float
    mac_margin: float
    warning: str


@dataclass(frozen=True)
class BranchComparison:
    branch_id: str
    base_sorted_index: int
    eb: BranchState
    timo: BranchState
    shape_mac: float
    delta_f: float


@dataclass(frozen=True)
class RunProducts:
    mode_rows: list[dict[str, object]]
    summary_rows: list[dict[str, object]]
    critical_rows: list[dict[str, object]]
    tracking_rows: list[dict[str, object]]
    threshold_rows: list[dict[str, object]]
    chi_fit_rows: list[dict[str, object]]
    warnings: list[str]
    solved_parameter_points: int
    threshold_refinement_iterations: int


class RootCache:
    def __init__(self, args: Args) -> None:
        self.args = args
        self.hits = 0
        self.misses = 0

    def settings(self, *, model: str, epsilon: float, mu: float) -> dict[str, object]:
        kappa = TIMO.circular_shear_coefficient(TIMO.NU)
        return {
            "cache_version": CACHE_VERSION,
            "model": str(model),
            "beta_deg": round(float(self.args.beta_deg), 12),
            "mu": round(float(mu), 12),
            "eta": round(float(self.args.eta), 12),
            "epsilon_0": round(float(epsilon), 12),
            "n_roots": int(self.args.n_candidate_roots),
            "shape_sample_count": int(self.args.n_shape_points),
            "k_prime": float(kappa),
            "root_settings": {
                "root_scan_start": beta_workflow.ROOT_SCAN_START,
                "eb_lambda_max": beta_workflow.EB_LAMBDA_MAX,
                "eb_retry_lambda_max": beta_workflow.EB_RETRY_LAMBDA_MAX,
                "eb_scan_step": beta_workflow.EB_SCAN_STEP,
                "eb_retry_scan_step": beta_workflow.EB_RETRY_SCAN_STEP,
                "timo_scan_step": beta_workflow.TIMO_SCAN_STEP,
                "timo_retry_scan_step": beta_workflow.TIMO_RETRY_SCAN_STEP,
                "timo_find_roots_grow_factor": 1.35,
                "timo_find_roots_max_tries": 8,
                "retry_missing_roots": True,
            },
        }

    def path(self, *, model: str, epsilon: float, mu: float) -> Path:
        settings_json = json.dumps(self.settings(model=model, epsilon=epsilon, mu=mu), sort_keys=True)
        digest = hashlib.sha256(settings_json.encode("utf-8")).hexdigest()[:20]
        model_token = "eb" if model == MODEL_EB else "timo"
        return self.args.cache_dir / f"root_{model_token}_{digest}.npz"

    def load(self, *, model: str, epsilon: float, mu: float) -> RootCacheEntry | None:
        if self.args.force_recompute or not self.args.reuse_cache:
            return None
        path = self.path(model=model, epsilon=epsilon, mu=mu)
        if not path.exists():
            return None
        data = np.load(path, allow_pickle=True)
        try:
            settings_json = str(data["settings_json"].item())
            expected = json.dumps(self.settings(model=model, epsilon=epsilon, mu=mu), sort_keys=True)
            if settings_json != expected:
                return None
            self.hits += 1
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

    def save(self, *, model: str, epsilon: float, mu: float, result: beta_workflow.RootResult) -> None:
        path = self.path(model=model, epsilon=epsilon, mu=mu)
        path.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            path,
            settings_json=np.array(
                json.dumps(self.settings(model=model, epsilon=epsilon, mu=mu), sort_keys=True),
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

    def solve_uncached(self, *, model: str, epsilon: float, mu: float) -> beta_workflow.RootResult:
        case = beta_workflow.CaseSpec(mu=float(mu), eta=self.args.eta, epsilon=float(epsilon))
        result = beta_workflow.solve_model(
            case,
            self.args.beta_deg,
            self.args.n_candidate_roots,
            model,
        )
        if result.root_count_found < self.args.n_candidate_roots:
            result = beta_workflow.retry_missing_roots(
                result,
                case,
                self.args.beta_deg,
                self.args.n_candidate_roots,
                model,
            )
        return result

    def roots(self, *, model: str, epsilon: float, mu: float) -> RootCacheEntry:
        cached = self.load(model=model, epsilon=epsilon, mu=mu)
        if cached is not None:
            return cached
        result = self.solve_uncached(model=model, epsilon=epsilon, mu=mu)
        self.save(model=model, epsilon=epsilon, mu=mu, result=result)
        self.misses += 1
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
        return "yes" if bool(value) else "no"
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


def append_timing_row(path: Path, row: dict[str, object]) -> Path:
    rows: list[dict[str, object]] = []
    if path.exists():
        rows.extend(read_csv(path))
    rows.append(row)
    return write_csv(path, rows, TIMING_FIELDS)


def bool_text(value: bool) -> str:
    return "yes" if bool(value) else "no"


def finite_float(row: dict[str, object] | dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def number_token(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def mu_filename_token(value: float) -> str:
    return f"{float(value):.1f}".replace("-", "m").replace(".", "p")


def parse_float_list(values: Sequence[str] | None, default: Sequence[float]) -> tuple[float, ...]:
    if values is None:
        return tuple(float(value) for value in default)
    return tuple(float(value) for value in values)


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Diagnostic-only EB applicability scan relative to the in-plane Timoshenko model.",
    )
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--skip-critical-refinement", action="store_true")
    parser.add_argument("--epsilon-values", type=str, nargs="+", default=None)
    parser.add_argument("--mu-values", type=str, nargs="+", default=None)
    parser.add_argument("--n-reported-modes", type=int, default=DEFAULT_N_REPORTED_MODES)
    parser.add_argument("--n-candidate-roots", type=int, default=DEFAULT_N_CANDIDATE_ROOTS)
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_SHAPE_SAMPLE_COUNT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--max-refinement-iterations", type=int, default=DEFAULT_MAX_REFINEMENT_ITERATIONS)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    mu_values = parse_float_list(ns.mu_values, DEFAULT_MU_VALUES)
    epsilon_values = parse_float_list(ns.epsilon_values, DEFAULT_EPSILON_VALUES)
    n_reported = int(ns.n_reported_modes)
    n_candidate = int(ns.n_candidate_roots)
    n_shape_points = int(ns.n_shape_points)
    output_dir = repo_path(Path(ns.output_dir))
    if ns.smoke:
        mu_values = SMOKE_MU_VALUES if ns.mu_values is None else mu_values
        epsilon_values = SMOKE_EPSILON_VALUES if ns.epsilon_values is None else epsilon_values
        n_reported = min(n_reported, SMOKE_N_REPORTED_MODES)
        n_candidate = min(n_candidate, SMOKE_N_CANDIDATE_ROOTS)
        n_shape_points = min(n_shape_points, SMOKE_SHAPE_SAMPLE_COUNT)
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            output_dir = repo_path(SMOKE_OUTPUT_DIR)
        ns.skip_critical_refinement = True
    cache_dir = repo_path(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / DEFAULT_CACHE_SUBDIR
    args = Args(
        beta_deg=DEFAULT_BETA_DEG,
        eta=DEFAULT_ETA,
        mu_values=tuple(float(value) for value in mu_values),
        epsilon_values=tuple(float(value) for value in epsilon_values),
        n_reported_modes=n_reported,
        n_candidate_roots=n_candidate,
        n_shape_points=n_shape_points,
        output_dir=output_dir,
        cache_dir=cache_dir,
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        plot_only=bool(ns.plot_only),
        skip_critical_refinement=bool(ns.skip_critical_refinement),
        smoke=bool(ns.smoke),
        max_refinement_iterations=int(ns.max_refinement_iterations),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if args.n_reported_modes < 1:
        raise ValueError("--n-reported-modes must be positive")
    if args.n_candidate_roots < args.n_reported_modes:
        raise ValueError("--n-candidate-roots must be at least --n-reported-modes")
    if args.n_shape_points < 51:
        raise ValueError("--n-shape-points must be at least 51")
    if not args.mu_values:
        raise ValueError("at least one mu value is required")
    if not args.epsilon_values:
        raise ValueError("at least one epsilon value is required")
    if tuple(sorted(args.epsilon_values)) != args.epsilon_values:
        raise ValueError("epsilon values must be sorted increasingly")
    if abs(args.beta_deg - 45.0) > 1.0e-14 or abs(args.eta) > 1.0e-14:
        raise ValueError("stage 1 is fixed at beta=45 deg and eta=0")
    for mu in args.mu_values:
        thickness_mismatch_factors(float(mu), args.eta)
    for epsilon in args.epsilon_values:
        if float(epsilon) <= 0.0:
            raise ValueError("epsilon values must be positive")
    if args.max_refinement_iterations < 0:
        raise ValueError("--max-refinement-iterations must be non-negative")


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


def pass_bool(delta_f: float, threshold: float) -> bool:
    return isfinite(float(delta_f)) and float(delta_f) <= float(threshold)


def consecutive_pass_count(deltas: Sequence[float], threshold: float) -> int:
    count = 0
    for value in deltas:
        if pass_bool(float(value), float(threshold)):
            count += 1
        else:
            break
    return count


def pass_pattern(deltas: Sequence[float], threshold: float) -> str:
    return ",".join("Y" if pass_bool(float(value), float(threshold)) else "N" for value in deltas)


def count_passes(deltas: Sequence[float], threshold: float) -> int:
    return sum(1 for value in deltas if pass_bool(float(value), float(threshold)))


def k_aware_sorted_summary_fields(
    deltas: Sequence[float],
    *,
    target_mode_count: int,
) -> dict[str, object]:
    target = [float(value) for value in deltas[: int(target_mode_count)]]
    legacy_first8 = target[:8]
    return {
        "target_mode_count": int(target_mode_count),
        "N5_sorted_consecutive": consecutive_pass_count(target, 0.05),
        "N10_sorted_consecutive": consecutive_pass_count(target, 0.10),
        "N15_sorted_consecutive": consecutive_pass_count(target, 0.15),
        "n_pass_5_among_first8": count_passes(legacy_first8, 0.05),
        "n_pass_10_among_first8": count_passes(legacy_first8, 0.10),
        "n_pass_15_among_first8": count_passes(legacy_first8, 0.15),
        "n_pass_5_among_firstK": count_passes(target, 0.05),
        "n_pass_10_among_firstK": count_passes(target, 0.10),
        "n_pass_15_among_firstK": count_passes(target, 0.15),
        "max_delta_f_firstK": max(target) if target else float("nan"),
        "mean_delta_f_firstK": float(np.mean(target)) if target else float("nan"),
        "pass_pattern_5": pass_pattern(target, 0.05),
        "pass_pattern_10": pass_pattern(target, 0.10),
        "pass_pattern_15": pass_pattern(target, 0.15),
    }


def local_thickness_parameters(epsilon_0: float, mu: float, eta: float) -> dict[str, float]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    eps1 = float(epsilon_0) * factors.tau1 / (1.0 - float(mu))
    eps2 = float(epsilon_0) * factors.tau2 / (1.0 + float(mu))
    return {
        "tau_1": float(factors.tau1),
        "tau_2": float(factors.tau2),
        "epsilon_1": float(eps1),
        "epsilon_2": float(eps2),
        "epsilon_max": max(float(eps1), float(eps2)),
        "radius_to_length_1": 2.0 * float(eps1),
        "radius_to_length_2": 2.0 * float(eps2),
        "diameter_to_length_1": 4.0 * float(eps1),
        "diameter_to_length_2": 4.0 * float(eps2),
        "l1_over_r1": 1.0 / (2.0 * float(eps1)),
        "l2_over_r2": 1.0 / (2.0 * float(eps2)),
        "l1_over_d1": 1.0 / (4.0 * float(eps1)),
        "l2_over_d2": 1.0 / (4.0 * float(eps2)),
    }


def chi_values(lambda_global: float, epsilon_0: float, mu: float, eta: float) -> dict[str, float]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    chi1 = float(epsilon_0) * float(lambda_global) * np.sqrt(factors.tau1)
    chi2 = float(epsilon_0) * float(lambda_global) * np.sqrt(factors.tau2)
    return {
        "chi_1": float(chi1),
        "chi_2": float(chi2),
        "chi_max": max(float(chi1), float(chi2)),
    }


def trapz(values: np.ndarray, coordinates: np.ndarray) -> float:
    return shape_audit.trapz(np.asarray(values, dtype=float), np.asarray(coordinates, dtype=float))


def bending_energies_by_rod(result: shape_audit.ModeResult) -> tuple[float, float]:
    factors = TIMO.tau_factors(float(result.mu), float(result.eta))
    section1 = TIMO.section_from_epsilon_tau(float(result.epsilon), factors.tau1)
    section2 = TIMO.section_from_epsilon_tau(float(result.epsilon), factors.tau2)
    if result.model == MODEL_TIMO:
        if result.rod1.psi_prime is None or result.rod2.psi_prime is None:
            return float("nan"), float("nan")
        e1 = 0.5 * section1.bending_stiffness * trapz(result.rod1.psi_prime**2, np.abs(result.rod1.x_local))
        e2 = 0.5 * section2.bending_stiffness * trapz(result.rod2.psi_prime**2, np.abs(result.rod2.x_local))
    else:
        if result.rod1.w_second is None or result.rod2.w_second is None:
            return float("nan"), float("nan")
        e1 = 0.5 * section1.bending_stiffness * trapz(result.rod1.w_second**2, np.abs(result.rod1.x_local))
        e2 = 0.5 * section2.bending_stiffness * trapz(result.rod2.w_second**2, np.abs(result.rod2.x_local))
    return max(float(e1), 0.0), max(float(e2), 0.0)


def chi_eff(lambda_global: float, epsilon_0: float, mu: float, eta: float, timo_result: shape_audit.ModeResult) -> float:
    chis = chi_values(float(lambda_global), float(epsilon_0), float(mu), float(eta))
    b1, b2 = bending_energies_by_rod(timo_result)
    denom = b1 + b2
    if not isfinite(denom) or denom <= 1.0e-30:
        return float("nan")
    value = (b1 * chis["chi_1"] ** 2 + b2 * chis["chi_2"] ** 2) / denom
    return float(np.sqrt(max(value, 0.0)))


def mode_vector(result: shape_audit.ModeResult) -> np.ndarray:
    sign = float(result.sign_factor)
    if result.model == MODEL_EB:
        dx1, dy1 = DISPLAY.eb_rod1_local_displacement_to_display(sign * result.rod1.u, sign * result.rod1.w)
        dx2, dy2 = DISPLAY.eb_rod2_local_displacement_to_display(
            sign * result.rod2.u,
            sign * result.rod2.w,
            beta_deg=float(result.beta_deg),
        )
    else:
        dx1, dy1 = DISPLAY.rod1_local_displacement_to_display(sign * result.rod1.u, sign * result.rod1.w)
        dx2, dy2 = DISPLAY.rod2_local_displacement_to_display(
            sign * result.rod2.u,
            sign * result.rod2.w,
            beta_deg=float(result.beta_deg),
        )
    vector = np.concatenate([dx1, dy1, dx2, dy2])
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-30:
        return np.full_like(vector, np.nan, dtype=float)
    return vector / norm


def mac_value(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return 0.0
    denom = float(np.dot(a, a) * np.dot(b, b))
    if denom <= 1.0e-30:
        return 0.0
    return abs(float(np.dot(a, b))) ** 2 / denom


def mac_assignment(previous: Sequence[np.ndarray], candidates: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mac = np.array(
        [[mac_value(prev, cand) for cand in candidates] for prev in previous],
        dtype=float,
    )
    rows, cols = linear_sum_assignment(1.0 - mac)
    assignment = np.full(len(previous), -1, dtype=int)
    for row, col in zip(rows, cols, strict=True):
        assignment[int(row)] = int(col)
    if np.any(assignment < 0):
        raise RuntimeError("MAC assignment failed")
    return assignment, mac


def mode_result_for_root(
    *,
    model: str,
    epsilon: float,
    mu: float,
    sorted_index: int,
    Lambda: float,
    args: Args,
    root_warnings: Sequence[str],
) -> shape_audit.ModeResult:
    common = {
        "epsilon": float(epsilon),
        "beta_deg": float(args.beta_deg),
        "eta": float(args.eta),
        "mu": float(mu),
        "sorted_index": int(sorted_index),
        "Lambda": float(Lambda),
        "n_points": int(args.n_shape_points),
        "case_role": "eb_validity_vs_timoshenko_stage1",
        "root_source": "root_cache_or_solver",
        "root_warnings": tuple(root_warnings),
    }
    if model == MODEL_EB:
        return shape_audit.eb_mode_result(**common)
    if model == MODEL_TIMO:
        return shape_audit.timo_mode_result(**common)
    raise ValueError(f"unknown model: {model}")


def build_candidate_states(
    *,
    model: str,
    epsilon: float,
    mu: float,
    args: Args,
    root_cache: RootCache,
) -> list[BranchState]:
    root_entry = root_cache.roots(model=model, epsilon=float(epsilon), mu=float(mu))
    states: list[BranchState] = []
    for index, Lambda in enumerate(root_entry.roots[: args.n_candidate_roots], start=1):
        if not isfinite(float(Lambda)):
            continue
        result = mode_result_for_root(
            model=model,
            epsilon=float(epsilon),
            mu=float(mu),
            sorted_index=index,
            Lambda=float(Lambda),
            args=args,
            root_warnings=root_entry.warnings,
        )
        states.append(
            BranchState(
                branch_id="candidate",
                base_sorted_index=index,
                current_sorted_index=index,
                Lambda=float(Lambda),
                result=result,
                vector=mode_vector(result),
                mac_to_previous=float("nan"),
                second_best_mac=float("nan"),
                mac_margin=float("nan"),
                warning="",
            )
        )
    return states


def sorted_state(
    *,
    model: str,
    epsilon: float,
    mu: float,
    sorted_index: int,
    args: Args,
    root_cache: RootCache,
) -> BranchState:
    root_entry = root_cache.roots(model=model, epsilon=float(epsilon), mu=float(mu))
    roots = root_entry.roots
    if int(sorted_index) < 1 or int(sorted_index) > len(roots):
        raise RuntimeError(f"sorted index {sorted_index} is outside root array")
    Lambda = float(roots[int(sorted_index) - 1])
    result = mode_result_for_root(
        model=model,
        epsilon=float(epsilon),
        mu=float(mu),
        sorted_index=int(sorted_index),
        Lambda=Lambda,
        args=args,
        root_warnings=root_entry.warnings,
    )
    return BranchState(
        branch_id=f"sorted_{int(sorted_index):02d}",
        base_sorted_index=int(sorted_index),
        current_sorted_index=int(sorted_index),
        Lambda=Lambda,
        result=result,
        vector=mode_vector(result),
        mac_to_previous=float("nan"),
        second_best_mac=float("nan"),
        mac_margin=float("nan"),
        warning="",
    )


def update_state_identity(
    state: BranchState,
    *,
    branch_id: str,
    base_sorted_index: int,
    mac_to_previous: float,
    second_best_mac: float,
) -> BranchState:
    margin = float(mac_to_previous) - float(second_best_mac) if isfinite(second_best_mac) else float("nan")
    warning_parts: list[str] = []
    if isfinite(float(mac_to_previous)) and float(mac_to_previous) < MAC_WARNING_THRESHOLD:
        warning_parts.append("low_mac")
    if isfinite(margin) and margin < MAC_MARGIN_WARNING_THRESHOLD:
        warning_parts.append("low_mac_margin")
    return replace(
        state,
        branch_id=str(branch_id),
        base_sorted_index=int(base_sorted_index),
        mac_to_previous=float(mac_to_previous),
        second_best_mac=float(second_best_mac),
        mac_margin=float(margin),
        warning="; ".join(warning_parts),
    )


def initial_branch_pairing(
    *,
    mu: float,
    args: Args,
    root_cache: RootCache,
) -> tuple[dict[str, BranchState], dict[str, BranchState], list[str]]:
    epsilon0 = float(args.epsilon_values[0])
    eb_candidates = build_candidate_states(model=MODEL_EB, epsilon=epsilon0, mu=float(mu), args=args, root_cache=root_cache)
    timo_candidates = build_candidate_states(model=MODEL_TIMO, epsilon=epsilon0, mu=float(mu), args=args, root_cache=root_cache)
    eb_seed = eb_candidates[: args.n_reported_modes]
    timo_seed_candidates = timo_candidates[: args.n_candidate_roots]
    assignment, mac = mac_assignment([state.vector for state in eb_seed], [state.vector for state in timo_seed_candidates])
    warnings: list[str] = []
    eb_states: dict[str, BranchState] = {}
    timo_states: dict[str, BranchState] = {}
    for branch_zero, timo_col in enumerate(assignment):
        branch_index = branch_zero + 1
        branch_id = f"base_branch_{branch_index}"
        eb_state = update_state_identity(
            eb_seed[branch_zero],
            branch_id=branch_id,
            base_sorted_index=branch_index,
            mac_to_previous=float("nan"),
            second_best_mac=float("nan"),
        )
        timo_row_macs = np.sort(mac[branch_zero])[::-1]
        timo_mac = float(mac[branch_zero, int(timo_col)])
        second = float(timo_row_macs[1]) if len(timo_row_macs) > 1 else float("nan")
        timo_state = update_state_identity(
            timo_seed_candidates[int(timo_col)],
            branch_id=branch_id,
            base_sorted_index=branch_index,
            mac_to_previous=timo_mac,
            second_best_mac=second,
        )
        if int(timo_col) + 1 != branch_index:
            warnings.append(
                f"base EB/Timoshenko branch pairing for mu={float(mu):g}, {branch_id}: "
                f"EB sorted {branch_index} matched Timoshenko sorted {int(timo_col) + 1}"
            )
        if timo_state.warning:
            warnings.append(f"base EB/Timoshenko shape MAC warning for mu={float(mu):g}, {branch_id}: {timo_state.warning}")
        eb_states[branch_id] = eb_state
        timo_states[branch_id] = timo_state
    return eb_states, timo_states, warnings


def track_one_epsilon_step(
    *,
    model: str,
    mu: float,
    epsilon_prev: float,
    epsilon_current: float,
    previous_states: dict[str, BranchState],
    args: Args,
    root_cache: RootCache,
) -> tuple[dict[str, BranchState], list[dict[str, object]], list[str]]:
    candidates = build_candidate_states(
        model=model,
        epsilon=float(epsilon_current),
        mu=float(mu),
        args=args,
        root_cache=root_cache,
    )
    branch_ids = list(previous_states)
    assignment, mac = mac_assignment([previous_states[branch_id].vector for branch_id in branch_ids], [state.vector for state in candidates])
    next_states: dict[str, BranchState] = {}
    audit_rows: list[dict[str, object]] = []
    warnings: list[str] = []
    for row_index, branch_id in enumerate(branch_ids):
        col = int(assignment[row_index])
        row_macs = np.sort(mac[row_index])[::-1]
        best = float(mac[row_index, col])
        second = float(row_macs[1]) if len(row_macs) > 1 else float("nan")
        previous = previous_states[branch_id]
        state = update_state_identity(
            candidates[col],
            branch_id=branch_id,
            base_sorted_index=previous.base_sorted_index,
            mac_to_previous=best,
            second_best_mac=second,
        )
        ambiguous = bool(state.warning)
        if ambiguous:
            warnings.append(
                f"{model} tracking warning mu={float(mu):g}, {branch_id}, "
                f"eps {float(epsilon_prev):g}->{float(epsilon_current):g}: {state.warning}"
            )
        audit_rows.append(
            {
                "model": model,
                "mu": float(mu),
                "eta": args.eta,
                "epsilon_prev": float(epsilon_prev),
                "epsilon_current": float(epsilon_current),
                "branch_id": branch_id,
                "previous_sorted_index": int(previous.current_sorted_index),
                "current_sorted_index": int(state.current_sorted_index),
                "Lambda": float(state.Lambda),
                "MAC_to_previous": best,
                "second_best_MAC": second,
                "MAC_margin": state.mac_margin,
                "ambiguous": bool_text(ambiguous),
                "warning": state.warning,
                "notes": "shape-MAC continuation over candidate roots",
            }
        )
        next_states[branch_id] = state
    return next_states, audit_rows, warnings


def base_tracking_audit_rows(
    states: dict[str, BranchState],
    *,
    model: str,
    mu: float,
    eta: float,
    epsilon: float,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for branch_id, state in states.items():
        rows.append(
            {
                "model": model,
                "mu": float(mu),
                "eta": float(eta),
                "epsilon_prev": float("nan"),
                "epsilon_current": float(epsilon),
                "branch_id": branch_id,
                "previous_sorted_index": "",
                "current_sorted_index": int(state.current_sorted_index),
                "Lambda": float(state.Lambda),
                "MAC_to_previous": float("nan"),
                "second_best_MAC": float("nan"),
                "MAC_margin": float("nan"),
                "ambiguous": "no",
                "warning": "",
                "notes": "seed at epsilon_0 minimum",
            }
        )
    return rows


def comparison_row(
    *,
    args: Args,
    comparison_type: str,
    branch_id: str,
    eb_state: BranchState,
    timo_state: BranchState,
    shape_mac: float,
    tracking_warning: str,
    notes: str,
) -> dict[str, object]:
    mu = float(eb_state.result.mu)
    epsilon = float(eb_state.result.epsilon)
    local = local_thickness_parameters(epsilon, mu, args.eta)
    metrics = frequency_metrics(eb_state.Lambda, timo_state.Lambda)
    eb_energy = eb_state.result.energy
    timo_energy = timo_state.result.energy
    chis = chi_values(timo_state.Lambda, epsilon, mu, args.eta)
    chi_eff_value = chi_eff(timo_state.Lambda, epsilon, mu, args.eta, timo_state.result)
    return {
        "beta_deg": args.beta_deg,
        "mu": mu,
        "eta": args.eta,
        "epsilon_0": epsilon,
        "tau_1": local["tau_1"],
        "tau_2": local["tau_2"],
        "epsilon_1": local["epsilon_1"],
        "epsilon_2": local["epsilon_2"],
        "epsilon_max": local["epsilon_max"],
        "l1_over_r1": local["l1_over_r1"],
        "l2_over_r2": local["l2_over_r2"],
        "l1_over_d1": local["l1_over_d1"],
        "l2_over_d2": local["l2_over_d2"],
        "comparison_type": comparison_type,
        "branch_id": branch_id,
        "eb_sorted_index": int(eb_state.current_sorted_index),
        "timo_sorted_index": int(timo_state.current_sorted_index),
        "Lambda_EB": float(eb_state.Lambda),
        "Lambda_Timo": float(timo_state.Lambda),
        "delta_Lambda": metrics["delta_Lambda"],
        "delta_f": metrics["delta_f"],
        "delta_f_symmetric": metrics["delta_f_symmetric"],
        "bias_f": metrics["bias_f"],
        "EB_axial_fraction": eb_energy["axial_energy_fraction"],
        "EB_bending_fraction": eb_energy["bending_energy_fraction"],
        "EB_classification": eb_energy["classification"],
        "Timo_axial_fraction": timo_energy["axial_energy_fraction"],
        "Timo_bending_fraction": timo_energy["bending_energy_fraction"],
        "Timo_shear_fraction": timo_energy["shear_energy_fraction"],
        "Timo_bending_shear_fraction": float(timo_energy["bending_energy_fraction"]) + float(timo_energy["shear_energy_fraction"]),
        "Timo_classification": timo_energy["classification"],
        "chi_1": chis["chi_1"],
        "chi_2": chis["chi_2"],
        "chi_max": chis["chi_max"],
        "chi_eff": chi_eff_value,
        "pass_5pct": bool_text(pass_bool(metrics["delta_f"], 0.05)),
        "pass_10pct": bool_text(pass_bool(metrics["delta_f"], 0.10)),
        "pass_15pct": bool_text(pass_bool(metrics["delta_f"], 0.15)),
        "shape_mac_EB_Timo": float(shape_mac),
        "tracking_warning": tracking_warning,
        "notes": notes,
    }


def sorted_comparison_rows(args: Args, root_cache: RootCache, *, epsilon: float, mu: float) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for sorted_index in range(1, args.n_reported_modes + 1):
        eb = sorted_state(model=MODEL_EB, epsilon=epsilon, mu=mu, sorted_index=sorted_index, args=args, root_cache=root_cache)
        timo = sorted_state(model=MODEL_TIMO, epsilon=epsilon, mu=mu, sorted_index=sorted_index, args=args, root_cache=root_cache)
        rows.append(
            comparison_row(
                args=args,
                comparison_type="sorted_index",
                branch_id=f"sorted_{sorted_index:02d}",
                eb_state=eb,
                timo_state=timo,
                shape_mac=mac_value(eb.vector, timo.vector),
                tracking_warning="",
                notes="sorted-spectrum diagnostic; branch_id is not physical identity",
            )
        )
    return rows


def physical_comparison_rows(
    args: Args,
    *,
    eb_states: dict[str, BranchState],
    timo_states: dict[str, BranchState],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for branch_id in eb_states:
        eb = eb_states[branch_id]
        timo = timo_states[branch_id]
        tracking = "; ".join(item for item in (eb.warning, timo.warning) if item)
        rows.append(
            comparison_row(
                args=args,
                comparison_type="physical_branch",
                branch_id=branch_id,
                eb_state=eb,
                timo_state=timo,
                shape_mac=mac_value(eb.vector, timo.vector),
                tracking_warning=tracking,
                notes="physical branch seeded at epsilon_0 minimum and continued by shape MAC",
            )
        )
    return rows


def summary_for_point(
    *,
    args: Args,
    epsilon: float,
    mu: float,
    sorted_rows: Sequence[dict[str, object]],
    physical_rows: Sequence[dict[str, object]],
) -> dict[str, object]:
    local = local_thickness_parameters(float(epsilon), float(mu), args.eta)
    sorted_deltas = [finite_float(row, "delta_f") for row in sorted(sorted_rows, key=lambda row: int(row["eb_sorted_index"]))]
    summary: dict[str, object] = {
        "beta_deg": args.beta_deg,
        "mu": float(mu),
        "eta": args.eta,
        "epsilon_0": float(epsilon),
        "epsilon_1": local["epsilon_1"],
        "epsilon_2": local["epsilon_2"],
        "epsilon_max": local["epsilon_max"],
        **k_aware_sorted_summary_fields(
            sorted_deltas,
            target_mode_count=args.n_reported_modes,
        ),
    }

    bending_rows = [
        row
        for row in physical_rows
        if str(row.get("Timo_classification", "")) == "bending_dominated"
    ]
    bending_rows.sort(key=lambda row: finite_float(row, "Lambda_Timo"))
    bending_deltas = [finite_float(row, "delta_f") for row in bending_rows]
    legacy_bending_rows = [
        row
        for row in physical_rows
        if str(row.get("Timo_classification", "")) == "bending_dominated"
        and str(row.get("branch_id", "")).startswith("base_branch_")
        and int(str(row["branch_id"]).rsplit("_", 1)[-1]) <= 8
    ]
    legacy_bending_rows.sort(key=lambda row: finite_float(row, "Lambda_Timo"))
    legacy_bending_deltas = [finite_float(row, "delta_f") for row in legacy_bending_rows]
    summary.update(
        {
            "N5_bending_consecutive": consecutive_pass_count(bending_deltas, 0.05),
            "N10_bending_consecutive": consecutive_pass_count(bending_deltas, 0.10),
            "N15_bending_consecutive": consecutive_pass_count(bending_deltas, 0.15),
            "n_bending_pass_10_among_first8": count_passes(legacy_bending_deltas, 0.10),
            "n_bending_pass_10_among_firstK": count_passes(bending_deltas, 0.10),
        }
    )

    first_failed_index = ""
    first_failed_delta = float("nan")
    for row in sorted(sorted_rows, key=lambda item: int(item["eb_sorted_index"])):
        delta = finite_float(row, "delta_f")
        if not pass_bool(delta, 0.10):
            first_failed_index = int(row["eb_sorted_index"])
            first_failed_delta = delta
            break

    first_failed_branch = ""
    first_failed_character = ""
    if first_failed_index != "":
        matches = [
            row
            for row in physical_rows
            if int(row["timo_sorted_index"]) == int(first_failed_index)
            or int(row["eb_sorted_index"]) == int(first_failed_index)
        ]
        if matches:
            match = sorted(matches, key=lambda row: abs(finite_float(row, "Lambda_Timo") - finite_float(sorted_rows[int(first_failed_index) - 1], "Lambda_Timo")))[0]
            first_failed_branch = str(match["branch_id"])
            first_failed_character = str(match["Timo_classification"])

    warnings = [
        str(row["tracking_warning"])
        for row in physical_rows
        if str(row.get("tracking_warning", ""))
    ]
    summary.update(
        {
            "first_failed_sorted_index_10": first_failed_index,
            "first_failed_branch_id_10": first_failed_branch,
            "first_failed_delta_f_10": first_failed_delta,
            "first_failed_character_10": first_failed_character,
            "warnings": " | ".join(warnings),
            "notes": "bending metrics use Timoshenko energy classification and physical branch ids",
        }
    )
    return summary


def find_threshold_brackets(
    x_values: Sequence[float],
    y_values: Sequence[float],
    threshold: float,
) -> list[tuple[float, float, float, float]]:
    xs = [float(value) for value in x_values]
    ys = [float(value) for value in y_values]
    brackets: list[tuple[float, float, float, float]] = []
    for left_index in range(len(xs) - 1):
        y_left = ys[left_index]
        y_right = ys[left_index + 1]
        if not (isfinite(y_left) and isfinite(y_right)):
            continue
        if y_left < float(threshold) <= y_right:
            brackets.append((xs[left_index], xs[left_index + 1], y_left, y_right))
    return brackets


def is_nonmonotonic(values: Sequence[float], *, tolerance: float = 1.0e-6) -> bool:
    finite_values = [float(value) for value in values if isfinite(float(value))]
    return any((right - left) < -float(tolerance) for left, right in zip(finite_values, finite_values[1:]))


def bisect_threshold(
    func: Callable[[float], float],
    low: float,
    high: float,
    *,
    threshold: float,
    eps_tol: float = CRITICAL_EPS_TOL,
    delta_tol: float = CRITICAL_DELTA_TOL,
    max_iterations: int = DEFAULT_MAX_REFINEMENT_ITERATIONS,
) -> tuple[float, float, int]:
    left = float(low)
    right = float(high)
    y_left = float(func(left))
    y_right = float(func(right))
    if not (y_left < float(threshold) <= y_right):
        raise ValueError("threshold bisection requires low<threshold<=high")
    iterations = 0
    for iteration in range(int(max_iterations)):
        iterations = iteration + 1
        mid = 0.5 * (left + right)
        y_mid = float(func(mid))
        if y_mid >= float(threshold):
            right = mid
            y_right = y_mid
        else:
            left = mid
            y_left = y_mid
        if abs(right - left) <= float(eps_tol) or abs(y_mid - float(threshold)) <= float(delta_tol):
            break
    del y_left
    return right, y_right, iterations


def select_candidate_by_endpoint_mac(
    *,
    model: str,
    epsilon: float,
    mu: float,
    low_state: BranchState,
    high_state: BranchState,
    args: Args,
    root_cache: RootCache,
) -> BranchState:
    candidates = build_candidate_states(
        model=model,
        epsilon=float(epsilon),
        mu=float(mu),
        args=args,
        root_cache=root_cache,
    )
    if not candidates:
        raise RuntimeError(f"no {model} candidates at epsilon={float(epsilon):g}, mu={float(mu):g}")
    scores = []
    for candidate in candidates:
        mac_low = mac_value(low_state.vector, candidate.vector)
        mac_high = mac_value(high_state.vector, candidate.vector)
        scores.append((0.5 * (mac_low + mac_high), mac_low, mac_high))
    best_index = int(np.argmax([score[0] for score in scores]))
    sorted_best = sorted((score[0] for score in scores), reverse=True)
    best_score, mac_low, mac_high = scores[best_index]
    second = float(sorted_best[1]) if len(sorted_best) > 1 else float("nan")
    state = update_state_identity(
        candidates[best_index],
        branch_id=low_state.branch_id,
        base_sorted_index=low_state.base_sorted_index,
        mac_to_previous=float(best_score),
        second_best_mac=second,
    )
    extra_warning = []
    if mac_low < MAC_WARNING_THRESHOLD:
        extra_warning.append("low_mac_to_low_bracket")
    if mac_high < MAC_WARNING_THRESHOLD:
        extra_warning.append("low_mac_to_high_bracket")
    if extra_warning:
        state = replace(state, warning="; ".join(item for item in (state.warning, *extra_warning) if item))
    return state


def refine_critical_branch(
    *,
    args: Args,
    root_cache: RootCache,
    mu: float,
    branch_id: str,
    low_eb: BranchState,
    low_timo: BranchState,
    high_eb: BranchState,
    high_timo: BranchState,
    delta_low: float,
    delta_high: float,
) -> tuple[float, float, float, float, BranchState, BranchState, list[dict[str, object]], int]:
    left = float(low_eb.result.epsilon)
    right = float(high_eb.result.epsilon)
    y_left = float(delta_low)
    y_right = float(delta_high)
    current_low_eb = low_eb
    current_low_timo = low_timo
    current_high_eb = high_eb
    current_high_timo = high_timo
    audit_rows: list[dict[str, object]] = []
    iterations = 0
    for iteration in range(1, int(args.max_refinement_iterations) + 1):
        iterations = iteration
        trial = 0.5 * (left + right)
        eb_trial = select_candidate_by_endpoint_mac(
            model=MODEL_EB,
            epsilon=trial,
            mu=float(mu),
            low_state=current_low_eb,
            high_state=current_high_eb,
            args=args,
            root_cache=root_cache,
        )
        timo_trial = select_candidate_by_endpoint_mac(
            model=MODEL_TIMO,
            epsilon=trial,
            mu=float(mu),
            low_state=current_low_timo,
            high_state=current_high_timo,
            args=args,
            root_cache=root_cache,
        )
        metrics = frequency_metrics(eb_trial.Lambda, timo_trial.Lambda)
        y_trial = metrics["delta_f"]
        shape_mac = mac_value(eb_trial.vector, timo_trial.vector)
        if y_trial >= MAIN_THRESHOLD:
            status = "high_updated"
            right = trial
            y_right = y_trial
            current_high_eb = eb_trial
            current_high_timo = timo_trial
        else:
            status = "low_updated"
            left = trial
            y_left = y_trial
            current_low_eb = eb_trial
            current_low_timo = timo_trial
        audit_rows.append(
            {
                "mu": float(mu),
                "branch_id": branch_id,
                "iteration": int(iteration),
                "epsilon_low": left,
                "epsilon_high": right,
                "epsilon_trial": trial,
                "delta_low": y_left,
                "delta_high": y_right,
                "delta_trial": y_trial,
                "EB_sorted_index": int(eb_trial.current_sorted_index),
                "Timo_sorted_index": int(timo_trial.current_sorted_index),
                "shape_MAC": shape_mac,
                "status": status,
                "notes": "; ".join(item for item in (eb_trial.warning, timo_trial.warning) if item),
            }
        )
        if abs(right - left) <= CRITICAL_EPS_TOL or abs(y_trial - MAIN_THRESHOLD) <= CRITICAL_DELTA_TOL:
            break
    return left, right, y_left, y_right, current_high_eb, current_high_timo, audit_rows, iterations


def critical_local_values(epsilon: float, mu: float, eta: float) -> dict[str, float]:
    local = local_thickness_parameters(float(epsilon), float(mu), float(eta))
    return {
        "epsilon_1_at_crit": local["epsilon_1"],
        "epsilon_2_at_crit": local["epsilon_2"],
        "epsilon_max_at_crit": local["epsilon_max"],
        "r0_over_l0_at_crit": 2.0 * float(epsilon),
        "d0_over_l0_at_crit": 4.0 * float(epsilon),
        "l1_over_r1_at_crit": local["l1_over_r1"],
        "l2_over_r2_at_crit": local["l2_over_r2"],
        "l1_over_d1_at_crit": local["l1_over_d1"],
        "l2_over_d2_at_crit": local["l2_over_d2"],
    }


def critical_rows_for_all_branches(
    *,
    args: Args,
    root_cache: RootCache,
    physical_by_mu_eps: dict[tuple[float, float], list[BranchComparison]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], int, list[str]]:
    critical_rows: list[dict[str, object]] = []
    threshold_audit_rows: list[dict[str, object]] = []
    warnings: list[str] = []
    refinement_iterations_total = 0
    for mu in args.mu_values:
        for branch_index in range(1, args.n_reported_modes + 1):
            branch_id = f"base_branch_{branch_index}"
            comparisons = [
                next(
                    item
                    for item in physical_by_mu_eps[(round(float(mu), 12), round(float(eps), 12))]
                    if item.branch_id == branch_id
                )
                for eps in args.epsilon_values
            ]
            deltas = [float(item.delta_f) for item in comparisons]
            nonmono = is_nonmonotonic(deltas)
            brackets = find_threshold_brackets(args.epsilon_values, deltas, MAIN_THRESHOLD)
            status = "bracketed" if brackets else "above_scan_max"
            epsilon_crit = float("nan")
            delta_low = float("nan")
            delta_high = float("nan")
            epsilon_low = float("nan")
            epsilon_high = float("nan")
            crit_eb: BranchState | None = None
            crit_timo: BranchState | None = None
            iterations = 0
            note_parts: list[str] = []
            warning = ""

            if isfinite(deltas[0]) and deltas[0] >= MAIN_THRESHOLD:
                status = "below_scan_min"
                epsilon_crit = float(args.epsilon_values[0])
                epsilon_low = float("nan")
                epsilon_high = float(args.epsilon_values[0])
                delta_low = float("nan")
                delta_high = float(deltas[0])
                crit_eb = comparisons[0].eb
                crit_timo = comparisons[0].timo
            elif brackets:
                epsilon_low, epsilon_high, delta_low, delta_high = brackets[0]
                low_index = list(args.epsilon_values).index(epsilon_low)
                high_index = list(args.epsilon_values).index(epsilon_high)
                low_comp = comparisons[low_index]
                high_comp = comparisons[high_index]
                if args.skip_critical_refinement:
                    status = "coarse_bracket_only"
                    epsilon_crit = float(epsilon_high)
                    crit_eb = high_comp.eb
                    crit_timo = high_comp.timo
                    note_parts.append("critical refinement skipped")
                else:
                    (
                        epsilon_low,
                        epsilon_crit,
                        delta_low,
                        delta_high,
                        crit_eb,
                        crit_timo,
                        audit,
                        iterations,
                    ) = refine_critical_branch(
                        args=args,
                        root_cache=root_cache,
                        mu=float(mu),
                        branch_id=branch_id,
                        low_eb=low_comp.eb,
                        low_timo=low_comp.timo,
                        high_eb=high_comp.eb,
                        high_timo=high_comp.timo,
                        delta_low=delta_low,
                        delta_high=delta_high,
                    )
                    epsilon_high = epsilon_crit
                    threshold_audit_rows.extend(audit)
                    refinement_iterations_total += iterations
            else:
                note_parts.append(f"lower_bound={float(args.epsilon_values[-1]):g}")
                crit_eb = comparisons[-1].eb
                crit_timo = comparisons[-1].timo

            if nonmono:
                note_parts.append("delta_f is nonmonotonic on the coarse grid")
            if len(brackets) > 1:
                note_parts.append(
                    "all_crossing_brackets="
                    + ";".join(f"[{left:g},{right:g}]" for left, right, _dl, _dh in brackets)
                )
            if crit_eb is not None and crit_eb.warning:
                warning = "; ".join(item for item in (warning, f"EB: {crit_eb.warning}") if item)
            if crit_timo is not None and crit_timo.warning:
                warning = "; ".join(item for item in (warning, f"Timoshenko: {crit_timo.warning}") if item)
            if warning:
                warnings.append(f"critical {branch_id}, mu={float(mu):g}: {warning}")

            local_at_crit = (
                critical_local_values(epsilon_crit, float(mu), args.eta)
                if isfinite(float(epsilon_crit))
                else {field: float("nan") for field in CRITICAL_FIELDS if field.endswith("_at_crit")}
            )
            chi_max_value = float("nan")
            chi_eff_value = float("nan")
            character = ""
            axial = float("nan")
            bending = float("nan")
            shear = float("nan")
            eb_sorted = ""
            timo_sorted = ""
            if crit_timo is not None and crit_eb is not None and isfinite(float(epsilon_crit)):
                chis = chi_values(crit_timo.Lambda, epsilon_crit, float(mu), args.eta)
                chi_max_value = chis["chi_max"]
                chi_eff_value = chi_eff(crit_timo.Lambda, epsilon_crit, float(mu), args.eta, crit_timo.result)
                character = str(crit_timo.result.energy["classification"])
                axial = float(crit_timo.result.energy["axial_energy_fraction"])
                bending = float(crit_timo.result.energy["bending_energy_fraction"])
                shear = float(crit_timo.result.energy["shear_energy_fraction"])
                eb_sorted = int(crit_eb.current_sorted_index)
                timo_sorted = int(crit_timo.current_sorted_index)

            row = {
                "beta_deg": args.beta_deg,
                "mu": float(mu),
                "eta": args.eta,
                "branch_id": branch_id,
                "base_sorted_index": int(branch_index),
                "epsilon_0_crit_10": epsilon_crit,
                "critical_status": status,
                "crossing_count": len(brackets),
                "nonmonotonic": bool_text(nonmono),
                "epsilon_low_bracket": epsilon_low,
                "epsilon_high_bracket": epsilon_high,
                "delta_low": delta_low,
                "delta_high": delta_high,
                "refinement_iterations": int(iterations),
                "epsilon_1_at_crit": local_at_crit.get("epsilon_1_at_crit", float("nan")),
                "epsilon_2_at_crit": local_at_crit.get("epsilon_2_at_crit", float("nan")),
                "epsilon_max_at_crit": local_at_crit.get("epsilon_max_at_crit", float("nan")),
                "r0_over_l0_at_crit": local_at_crit.get("r0_over_l0_at_crit", float("nan")),
                "d0_over_l0_at_crit": local_at_crit.get("d0_over_l0_at_crit", float("nan")),
                "l1_over_r1_at_crit": local_at_crit.get("l1_over_r1_at_crit", float("nan")),
                "l2_over_r2_at_crit": local_at_crit.get("l2_over_r2_at_crit", float("nan")),
                "l1_over_d1_at_crit": local_at_crit.get("l1_over_d1_at_crit", float("nan")),
                "l2_over_d2_at_crit": local_at_crit.get("l2_over_d2_at_crit", float("nan")),
                "chi_max_at_crit": chi_max_value,
                "chi_eff_at_crit": chi_eff_value,
                "EB_sorted_index_at_crit": eb_sorted,
                "Timo_sorted_index_at_crit": timo_sorted,
                "mode_character_at_crit": character,
                "axial_fraction_at_crit": axial,
                "bending_fraction_at_crit": bending,
                "shear_fraction_at_crit": shear,
                "warning": warning,
                "notes": "; ".join(note_parts),
            }
            critical_rows.append(row)
    return critical_rows, threshold_audit_rows, refinement_iterations_total, warnings


def fit_power_law(rows: Sequence[dict[str, object]], metric: str) -> dict[str, object]:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        if str(row.get("comparison_type", "")) != "physical_branch":
            continue
        if str(row.get("Timo_classification", "")) != "bending_dominated":
            continue
        x = finite_float(row, metric)
        y = finite_float(row, "delta_f")
        if isfinite(x) and isfinite(y) and x > 0.0 and y > 1.0e-12:
            xs.append(x)
            ys.append(y)
    if len(xs) < 3:
        return {
            "metric": metric,
            "point_count": len(xs),
            "C": float("nan"),
            "p": float("nan"),
            "R2": float("nan"),
            "chi_crit_for_delta_f_0p10": float("nan"),
            "notes": "not enough bending-dominated points",
        }
    logx = np.log(np.asarray(xs, dtype=float))
    logy = np.log(np.asarray(ys, dtype=float))
    p, logc = np.polyfit(logx, logy, 1)
    pred = p * logx + logc
    ss_res = float(np.sum((logy - pred) ** 2))
    ss_tot = float(np.sum((logy - float(np.mean(logy))) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else float("nan")
    c = float(np.exp(logc))
    chi_crit = float((MAIN_THRESHOLD / c) ** (1.0 / p)) if c > 0.0 and p != 0.0 else float("nan")
    return {
        "metric": metric,
        "point_count": len(xs),
        "C": c,
        "p": float(p),
        "R2": r2,
        "chi_crit_for_delta_f_0p10": chi_crit,
        "notes": "exploratory log-log fit for bending-dominated branch rows only",
    }


def build_chi_fit_rows(mode_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    return [fit_power_law(mode_rows, "chi_max"), fit_power_law(mode_rows, "chi_eff")]


def solve_products(args: Args, root_cache: RootCache) -> RunProducts:
    mode_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    tracking_rows: list[dict[str, object]] = []
    warnings: list[str] = []
    physical_by_mu_eps: dict[tuple[float, float], list[BranchComparison]] = {}

    for mu in args.mu_values:
        print(f"tracking physical branches for mu={float(mu):g}")
        eb_states, timo_states, base_warnings = initial_branch_pairing(mu=float(mu), args=args, root_cache=root_cache)
        warnings.extend(base_warnings)
        tracking_rows.extend(
            base_tracking_audit_rows(
                eb_states,
                model=MODEL_EB,
                mu=float(mu),
                eta=args.eta,
                epsilon=float(args.epsilon_values[0]),
            )
        )
        tracking_rows.extend(
            base_tracking_audit_rows(
                timo_states,
                model=MODEL_TIMO,
                mu=float(mu),
                eta=args.eta,
                epsilon=float(args.epsilon_values[0]),
            )
        )
        for eps_index, epsilon in enumerate(args.epsilon_values):
            if eps_index > 0:
                previous_epsilon = float(args.epsilon_values[eps_index - 1])
                eb_states, eb_audit, eb_warnings = track_one_epsilon_step(
                    model=MODEL_EB,
                    mu=float(mu),
                    epsilon_prev=previous_epsilon,
                    epsilon_current=float(epsilon),
                    previous_states=eb_states,
                    args=args,
                    root_cache=root_cache,
                )
                timo_states, timo_audit, timo_warnings = track_one_epsilon_step(
                    model=MODEL_TIMO,
                    mu=float(mu),
                    epsilon_prev=previous_epsilon,
                    epsilon_current=float(epsilon),
                    previous_states=timo_states,
                    args=args,
                    root_cache=root_cache,
                )
                tracking_rows.extend(eb_audit)
                tracking_rows.extend(timo_audit)
                warnings.extend(eb_warnings)
                warnings.extend(timo_warnings)

            sorted_rows = sorted_comparison_rows(args, root_cache, epsilon=float(epsilon), mu=float(mu))
            physical_rows = physical_comparison_rows(args, eb_states=eb_states, timo_states=timo_states)
            mode_rows.extend(sorted_rows)
            mode_rows.extend(physical_rows)
            comparisons = [
                BranchComparison(
                    branch_id=str(row["branch_id"]),
                    base_sorted_index=int(str(row["branch_id"]).split("_")[-1]),
                    eb=eb_states[str(row["branch_id"])],
                    timo=timo_states[str(row["branch_id"])],
                    shape_mac=finite_float(row, "shape_mac_EB_Timo"),
                    delta_f=finite_float(row, "delta_f"),
                )
                for row in physical_rows
            ]
            physical_by_mu_eps[(round(float(mu), 12), round(float(epsilon), 12))] = comparisons
            summary_rows.append(
                summary_for_point(
                    args=args,
                    epsilon=float(epsilon),
                    mu=float(mu),
                    sorted_rows=sorted_rows,
                    physical_rows=physical_rows,
                )
            )

    critical_rows, threshold_rows, refinement_iterations, critical_warnings = critical_rows_for_all_branches(
        args=args,
        root_cache=root_cache,
        physical_by_mu_eps=physical_by_mu_eps,
    )
    warnings.extend(critical_warnings)
    chi_fit_rows = build_chi_fit_rows(mode_rows)
    return RunProducts(
        mode_rows=mode_rows,
        summary_rows=summary_rows,
        critical_rows=critical_rows,
        tracking_rows=tracking_rows,
        threshold_rows=threshold_rows,
        chi_fit_rows=chi_fit_rows,
        warnings=warnings,
        solved_parameter_points=len(args.mu_values) * len(args.epsilon_values),
        threshold_refinement_iterations=refinement_iterations,
    )


def heatmap_from_summary(
    summary_rows: Sequence[dict[str, object] | dict[str, str]],
    *,
    value_key: str,
    output_path: Path,
    title: str,
) -> Path:
    mus = sorted({finite_float(row, "mu") for row in summary_rows})
    epsilons = sorted({finite_float(row, "epsilon_0") for row in summary_rows})
    values = np.full((len(epsilons), len(mus)), np.nan, dtype=float)
    for row in summary_rows:
        mu = finite_float(row, "mu")
        eps = finite_float(row, "epsilon_0")
        if not (isfinite(mu) and isfinite(eps)):
            continue
        i = epsilons.index(eps)
        j = mus.index(mu)
        values[i, j] = finite_float(row, value_key)
    fig, ax = plt.subplots(figsize=(9.0, 5.8), constrained_layout=True)
    im = ax.imshow(values, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(mus)))
    ax.set_xticklabels([f"{value:g}" for value in mus])
    ax.set_yticks(range(len(epsilons)))
    ax.set_yticklabels([f"{value:g}" for value in epsilons])
    ax.set_xlabel("mu")
    ax.set_ylabel("epsilon_0")
    ax.set_title(title)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            if isfinite(float(values[i, j])):
                ax.text(j, i, f"{int(values[i, j])}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=ax, label=value_key)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_critical_by_branch(
    critical_rows: Sequence[dict[str, object] | dict[str, str]],
    *,
    value_key: str,
    output_path: Path,
    ylabel: str,
) -> Path:
    mus = sorted({finite_float(row, "mu") for row in critical_rows})
    fig, ax = plt.subplots(figsize=(9.0, 5.4), constrained_layout=True)
    for mu in mus:
        rows = [row for row in critical_rows if abs(finite_float(row, "mu") - mu) <= 1.0e-12]
        rows.sort(key=lambda row: finite_float(row, "base_sorted_index"))
        x = [finite_float(row, "base_sorted_index") for row in rows]
        y = [finite_float(row, value_key) for row in rows]
        ax.plot(x, y, marker="o", lw=1.5, label=f"mu={mu:g}")
    ax.set_xlabel("base branch index")
    ax.set_ylabel(ylabel)
    ax.set_title(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=2, fontsize=8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_delta_curves(
    mode_rows: Sequence[dict[str, object] | dict[str, str]],
    *,
    mu: float,
    output_path: Path,
) -> Path:
    rows = [
        row
        for row in mode_rows
        if str(row.get("comparison_type", "")) == "physical_branch"
        and abs(finite_float(row, "mu") - float(mu)) <= 1.0e-12
    ]
    fig, ax = plt.subplots(figsize=(9.0, 5.4), constrained_layout=True)
    branch_indices = sorted(
        {
            int(str(row.get("branch_id", "")).rsplit("_", 1)[-1])
            for row in rows
            if str(row.get("branch_id", "")).startswith("base_branch_")
        }
    )
    for branch_index in branch_indices:
        branch_id = f"base_branch_{branch_index}"
        branch_rows = [row for row in rows if str(row.get("branch_id", "")) == branch_id]
        if not branch_rows:
            continue
        branch_rows.sort(key=lambda row: finite_float(row, "epsilon_0"))
        ax.plot(
            [finite_float(row, "epsilon_0") for row in branch_rows],
            [finite_float(row, "delta_f") for row in branch_rows],
            marker="o",
            lw=1.2,
            label=branch_id.replace("base_branch_", "b"),
        )
    ax.axhline(MAIN_THRESHOLD, color="black", lw=1.0, ls="--", label="10%")
    ax.set_xlabel("epsilon_0")
    ax.set_ylabel("delta_f")
    ax.set_title(f"EB/Timoshenko frequency divergence, mu={float(mu):g}")
    ax.grid(True, alpha=0.3)
    ax.legend(ncol=4, fontsize=8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def plot_chi(mode_rows: Sequence[dict[str, object] | dict[str, str]], output_path: Path) -> Path:
    rows = [
        row
        for row in mode_rows
        if str(row.get("comparison_type", "")) == "physical_branch"
        and str(row.get("Timo_classification", "")) == "bending_dominated"
        and isfinite(finite_float(row, "chi_max"))
        and isfinite(finite_float(row, "delta_f"))
        and finite_float(row, "chi_max") > 0.0
        and finite_float(row, "delta_f") > 0.0
    ]
    fig, ax = plt.subplots(figsize=(8.0, 5.6), constrained_layout=True)
    branches = sorted({str(row.get("branch_id", "")) for row in rows})
    for branch in branches:
        branch_rows = [row for row in rows if str(row.get("branch_id", "")) == branch]
        ax.scatter(
            [finite_float(row, "chi_max") for row in branch_rows],
            [finite_float(row, "delta_f") for row in branch_rows],
            s=28,
            label=branch.replace("base_branch_", "b"),
        )
    ax.axhline(MAIN_THRESHOLD, color="black", lw=1.0, ls="--")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("chi_max")
    ax.set_ylabel("delta_f")
    ax.set_title("Bending-dominated branch rows")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(ncol=4, fontsize=8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)
    return output_path


def generate_plots(
    *,
    output_dir: Path,
    mode_rows: Sequence[dict[str, object] | dict[str, str]],
    summary_rows: Sequence[dict[str, object] | dict[str, str]],
    critical_rows: Sequence[dict[str, object] | dict[str, str]],
) -> list[Path]:
    paths = [
        heatmap_from_summary(
            summary_rows,
            value_key="N10_sorted_consecutive",
            output_path=output_dir / "N10_sorted_consecutive_heatmap.png",
            title="N10 sorted consecutive",
        ),
        heatmap_from_summary(
            summary_rows,
            value_key="N10_bending_consecutive",
            output_path=output_dir / "N10_bending_consecutive_heatmap.png",
            title="N10 bending consecutive",
        ),
        plot_critical_by_branch(
            critical_rows,
            value_key="epsilon_0_crit_10",
            output_path=output_dir / "critical_epsilon10_by_branch_and_mu.png",
            ylabel="epsilon_0_crit_10",
        ),
        plot_critical_by_branch(
            critical_rows,
            value_key="epsilon_max_at_crit",
            output_path=output_dir / "critical_epsilon_max10_by_branch_and_mu.png",
            ylabel="epsilon_max_at_crit",
        ),
        plot_chi(mode_rows, output_dir / "delta_f_vs_chi_max_bending_modes.png"),
    ]
    for mu in (0.0, 0.3, 0.7):
        if any(abs(finite_float(row, "mu") - mu) <= 1.0e-12 for row in mode_rows):
            paths.append(
                plot_delta_curves(
                    mode_rows,
                    mu=mu,
                    output_path=output_dir / f"delta_f_vs_epsilon_mu{mu_filename_token(mu)}.png",
                )
            )
    return paths


def main_applicability_markdown(summary_rows: Sequence[dict[str, object] | dict[str, str]]) -> list[str]:
    selected_mus = [0.0, 0.3, 0.7]
    lines = [
        "| epsilon_0 | mu | N10 sorted consecutive | N10 bending consecutive | pass pattern 10 |",
        "|---:|---:|---:|---:|---|",
    ]
    for row in sorted(summary_rows, key=lambda item: (finite_float(item, "epsilon_0"), finite_float(item, "mu"))):
        mu = finite_float(row, "mu")
        if not any(abs(mu - target) <= 1.0e-12 for target in selected_mus):
            continue
        lines.append(
            "| "
            f"{finite_float(row, 'epsilon_0'):.5g} | "
            f"{mu:.3g} | "
            f"{int(finite_float(row, 'N10_sorted_consecutive'))} | "
            f"{int(finite_float(row, 'N10_bending_consecutive'))} | "
            f"{row.get('pass_pattern_10', '')} |"
        )
    return lines


def critical_markdown(critical_rows: Sequence[dict[str, object] | dict[str, str]], *, limit: int = 24) -> list[str]:
    lines = [
        "| mu | branch | status | epsilon_0_crit_10 | epsilon_max_crit | character | l1/r1 | l2/r2 |",
        "|---:|---|---|---:|---:|---|---:|---:|",
    ]
    sorted_rows = sorted(critical_rows, key=lambda row: (finite_float(row, "mu"), finite_float(row, "base_sorted_index")))
    for row in sorted_rows[:limit]:
        lines.append(
            "| "
            f"{finite_float(row, 'mu'):.3g} | "
            f"{row.get('branch_id', '')} | "
            f"{row.get('critical_status', '')} | "
            f"{finite_float(row, 'epsilon_0_crit_10'):.5g} | "
            f"{finite_float(row, 'epsilon_max_at_crit'):.5g} | "
            f"{row.get('mode_character_at_crit', '')} | "
            f"{finite_float(row, 'l1_over_r1_at_crit'):.5g} | "
            f"{finite_float(row, 'l2_over_r2_at_crit'):.5g} |"
        )
    if len(sorted_rows) > limit:
        lines.append(f"| ... | ... | ... | ... | ... | ... | ... | ... |")
    return lines


def first_fail_lines(summary_rows: Sequence[dict[str, object] | dict[str, str]]) -> list[str]:
    lines: list[str] = []
    for mu in sorted({finite_float(row, "mu") for row in summary_rows}):
        rows = [row for row in summary_rows if abs(finite_float(row, "mu") - mu) <= 1.0e-12]
        rows.sort(key=lambda row: finite_float(row, "epsilon_0"))
        first = next((row for row in rows if finite_float(row, "first_failed_delta_f_10") >= MAIN_THRESHOLD), None)
        if first is None:
            lines.append(f"- mu={mu:g}: no sorted-mode 10% failure on the scanned epsilon_0 grid.")
        else:
            lines.append(
                "- "
                f"mu={mu:g}: first sorted failure at epsilon_0={finite_float(first, 'epsilon_0'):.5g}, "
                f"sorted k={first.get('first_failed_sorted_index_10', '')}, "
                f"branch={first.get('first_failed_branch_id_10', '')}, "
                f"character={first.get('first_failed_character_10', '')}, "
                f"delta_f={finite_float(first, 'first_failed_delta_f_10'):.4g}."
            )
    return lines


def write_report(
    *,
    output_dir: Path,
    args: Args,
    mode_rows: Sequence[dict[str, object] | dict[str, str]],
    summary_rows: Sequence[dict[str, object] | dict[str, str]],
    critical_rows: Sequence[dict[str, object] | dict[str, str]],
    tracking_rows: Sequence[dict[str, object] | dict[str, str]],
    threshold_rows: Sequence[dict[str, object] | dict[str, str]],
    chi_fit_rows: Sequence[dict[str, object] | dict[str, str]],
    warnings: Sequence[str],
    plot_paths: Sequence[Path],
) -> Path:
    report = output_dir / "eb_validity_vs_timoshenko_stage1_report.md"
    nonmonotonic_count = sum(1 for row in critical_rows if str(row.get("nonmonotonic", "")).lower() in {"yes", "true", "1"})
    tracking_warning_count = sum(1 for row in tracking_rows if str(row.get("warning", "")))
    threshold_iteration_count = len(threshold_rows)
    long_or_mixed_pass_after_fail = [
        row
        for row in summary_rows
        if "N,Y" in str(row.get("pass_pattern_10", "")) or "N,N,Y" in str(row.get("pass_pattern_10", ""))
    ]
    lines = [
        "# EB Validity Relative To Timoshenko, Stage 1",
        "",
        "## Scope",
        "",
        "This diagnostic estimates the boundary of applicability of Euler--Bernoulli theory relative to the in-plane Timoshenko beam model. It reports divergence between Euler--Bernoulli and Timoshenko theories, not an exact error relative to a real 3D elastic body.",
        "",
        "No FEM, 3D FEM, Gmsh, or CalculiX calculation is part of this stage.",
        "",
        "## Scientific Definitions",
        "",
        "- For a fixed physical parameter set, dimensional frequency has the common factor form `f = C * Lambda^2`. The shared factor `C` cancels, so the main metric is `delta_f = abs(Lambda_EB^2 - Lambda_Timo^2) / Lambda_Timo^2`.",
        "- `delta_Lambda = abs(Lambda_EB - Lambda_Timo) / Lambda_Timo` is saved as a secondary metric.",
        "- `bias_f = (Lambda_EB^2 - Lambda_Timo^2) / Lambda_Timo^2`; positive values mean EB is higher than Timoshenko.",
        "- `delta_f_symmetric` is saved for diagnostics, while the 10% criterion uses the Timoshenko denominator.",
        "- Sorted-spectrum comparison compares EB sorted index `k` to Timoshenko sorted index `k`.",
        f"- Physical-branch comparison seeds `base_branch_1..base_branch_{args.n_reported_modes}` at `epsilon_0=0.0025` and continues EB and Timoshenko branches separately by shape MAC over candidate roots.",
        "",
        "## Parameter Grid",
        "",
        f"- beta: `{args.beta_deg:g} deg`",
        f"- eta: `{args.eta:g}`",
        f"- mu values: `{', '.join(f'{value:g}' for value in args.mu_values)}`",
        f"- epsilon_0 values: `{', '.join(f'{value:g}' for value in args.epsilon_values)}`",
        f"- reported modes: `{args.n_reported_modes}`",
        f"- candidate roots for matching: `{args.n_candidate_roots}`",
        "",
        "## Main Applicability Table",
        "",
        *main_applicability_markdown(summary_rows),
        "",
        "## Critical Thickness Summary",
        "",
        *critical_markdown(critical_rows),
        "",
        "First sorted failures on the scanned grid:",
        "",
        *first_fail_lines(summary_rows),
        "",
        "## Longitudinal-Mode Caveat",
        "",
        "A longitudinal-dominated mode can pass the 10% criterion after a neighboring bending-dominated mode has already failed because Timoshenko shear and rotary-inertia corrections mainly affect bending/shear kinematics. Therefore the report stores `N10_sorted_consecutive`, the legacy `n_pass_10_among_first8`, the K-aware `n_pass_10_among_firstK`, and a pass pattern such as `Y,Y,N,Y`.",
        "",
        f"Nonconsecutive pass patterns with a later pass after a failure: `{len(long_or_mixed_pass_after_fail)}` rows.",
        "",
        "## Nonmonotonic And Reordering Diagnostics",
        "",
        f"- Tracking warning rows: `{tracking_warning_count}`",
        f"- Critical branches with nonmonotonic coarse `delta_f`: `{nonmonotonic_count}`",
        f"- Threshold refinement iterations: `{threshold_iteration_count}`",
        "",
        "Branch ids and current sorted indices are stored separately. Current sorted index is metadata only.",
        "",
        "## Exploratory Chi Scaling",
        "",
        "| metric | points | C | p | R2 | chi at delta_f=0.10 |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in chi_fit_rows:
        lines.append(
            "| "
            f"{row.get('metric', '')} | "
            f"{int(finite_float(row, 'point_count')) if isfinite(finite_float(row, 'point_count')) else 0} | "
            f"{finite_float(row, 'C'):.6g} | "
            f"{finite_float(row, 'p'):.6g} | "
            f"{finite_float(row, 'R2'):.6g} | "
            f"{finite_float(row, 'chi_crit_for_delta_f_0p10'):.6g} |"
        )
    lines.extend(
        [
            "",
            "This is exploratory diagnostic scaling for bending-dominated branch rows only. It is not presented as a universal law without additional eta scans.",
            "",
            "## Conclusion",
            "",
            "For stage 1, `chi_max` and `chi_eff` are the most directly modal predictors because they combine frequency level and local thickness. `epsilon_0` is convenient for input scans, while `epsilon_max` is better for engineering slenderness conversion, especially when `mu` changes local length. A further eta scan is justified before using any fitted chi scaling as a general criterion.",
            "",
            "## Limitations",
            "",
            "- This is EB vs Timoshenko, not EB vs exact 3D elasticity.",
            "- No 3D FEM was run.",
            "- Stage 1 uses only `beta=45 deg` and `eta=0`.",
            "- The Timoshenko model is a fuller beam reference model, not an absolute truth.",
            "",
            "## Outputs",
            "",
            "- `eb_timo_mode_level_metrics.csv`",
            "- `eb_timo_validity_summary.csv`",
            "- `eb_timo_critical_thickness_by_branch.csv`",
            "- `epsilon_branch_tracking_audit.csv`",
            "- `critical_threshold_refinement_audit.csv`",
            "- `chi_scaling_fit_summary.csv`",
        ]
    )
    for path in plot_paths:
        lines.append(f"- `{rel(path)}`")
    if warnings:
        lines.extend(["", "## Warnings", ""])
        for warning in warnings[:80]:
            lines.append(f"- {warning}")
        if len(warnings) > 80:
            lines.append(f"- ... {len(warnings) - 80} additional warnings in audit CSVs.")
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def write_outputs(
    *,
    args: Args,
    products: RunProducts,
    run_seconds: float,
    root_cache: RootCache,
    run_kind: str,
) -> tuple[list[Path], Path]:
    output_dir = args.output_dir
    csv_paths = [
        write_csv(output_dir / "eb_timo_mode_level_metrics.csv", products.mode_rows, MODE_LEVEL_FIELDS),
        write_csv(output_dir / "eb_timo_validity_summary.csv", products.summary_rows, SUMMARY_FIELDS),
        write_csv(output_dir / "eb_timo_critical_thickness_by_branch.csv", products.critical_rows, CRITICAL_FIELDS),
        write_csv(output_dir / "epsilon_branch_tracking_audit.csv", products.tracking_rows, TRACKING_AUDIT_FIELDS),
        write_csv(output_dir / "critical_threshold_refinement_audit.csv", products.threshold_rows, THRESHOLD_AUDIT_FIELDS),
        write_csv(output_dir / "chi_scaling_fit_summary.csv", products.chi_fit_rows, CHI_FIT_FIELDS),
    ]
    timing_row = {
        "run_kind": run_kind,
        "seconds": float(run_seconds),
        "cache_hits": int(root_cache.hits),
        "cache_misses": int(root_cache.misses),
        "solved_parameter_points": int(products.solved_parameter_points),
        "threshold_refinement_iterations": int(products.threshold_refinement_iterations),
        "notes": f"timestamp={datetime.now().isoformat(timespec='seconds')}",
    }
    csv_paths.append(append_timing_row(output_dir / "timing_report.csv", timing_row))
    plot_paths = generate_plots(
        output_dir=output_dir,
        mode_rows=products.mode_rows,
        summary_rows=products.summary_rows,
        critical_rows=products.critical_rows,
    )
    report = write_report(
        output_dir=output_dir,
        args=args,
        mode_rows=products.mode_rows,
        summary_rows=products.summary_rows,
        critical_rows=products.critical_rows,
        tracking_rows=products.tracking_rows,
        threshold_rows=products.threshold_rows,
        chi_fit_rows=products.chi_fit_rows,
        warnings=products.warnings,
        plot_paths=plot_paths,
    )
    return [*csv_paths, *plot_paths], report


def warnings_from_saved_audits(
    mode_rows: Sequence[dict[str, object] | dict[str, str]],
    tracking_rows: Sequence[dict[str, object] | dict[str, str]],
) -> list[str]:
    warnings: list[str] = []
    physical_rows = [
        row
        for row in mode_rows
        if str(row.get("comparison_type", "")) == "physical_branch"
    ]
    for mu in sorted({finite_float(row, "mu") for row in physical_rows}):
        mu_rows = [
            row
            for row in physical_rows
            if abs(finite_float(row, "mu") - mu) <= 1.0e-12
        ]
        epsilon_values = [
            finite_float(row, "epsilon_0")
            for row in mu_rows
            if isfinite(finite_float(row, "epsilon_0"))
        ]
        if not epsilon_values:
            continue
        epsilon_min = min(epsilon_values)
        for row in mu_rows:
            if abs(finite_float(row, "epsilon_0") - epsilon_min) > 1.0e-12:
                continue
            eb_index = int(finite_float(row, "eb_sorted_index"))
            timo_index = int(finite_float(row, "timo_sorted_index"))
            if eb_index != timo_index:
                warnings.append(
                    f"base EB/Timoshenko branch pairing for mu={mu:g}, {row.get('branch_id', '')}: "
                    f"EB sorted {eb_index} matched Timoshenko sorted {timo_index}"
                )
            mode_warning = str(row.get("tracking_warning", "")).strip()
            if mode_warning:
                warnings.append(
                    f"saved physical-branch warning for mu={mu:g}, {row.get('branch_id', '')}: "
                    f"{mode_warning}"
                )
    for row in tracking_rows:
        warning = str(row.get("warning", "")).strip()
        if not warning:
            continue
        warnings.append(
            f"{row.get('model', '')} tracking warning mu={finite_float(row, 'mu'):g}, "
            f"{row.get('branch_id', '')}, eps {finite_float(row, 'epsilon_prev'):g}->"
            f"{finite_float(row, 'epsilon_current'):g}: {warning}"
        )
    return list(dict.fromkeys(warnings))


def plot_only_outputs(args: Args) -> tuple[list[Path], Path]:
    output_dir = args.output_dir
    mode_rows = read_csv(output_dir / "eb_timo_mode_level_metrics.csv")
    summary_rows = read_csv(output_dir / "eb_timo_validity_summary.csv")
    critical_rows = read_csv(output_dir / "eb_timo_critical_thickness_by_branch.csv")
    tracking_rows = read_csv(output_dir / "epsilon_branch_tracking_audit.csv")
    threshold_rows = read_csv(output_dir / "critical_threshold_refinement_audit.csv")
    chi_fit_path = output_dir / "chi_scaling_fit_summary.csv"
    chi_fit_rows = read_csv(chi_fit_path) if chi_fit_path.exists() else build_chi_fit_rows(mode_rows)
    plot_paths = generate_plots(
        output_dir=output_dir,
        mode_rows=mode_rows,
        summary_rows=summary_rows,
        critical_rows=critical_rows,
    )
    report = write_report(
        output_dir=output_dir,
        args=args,
        mode_rows=mode_rows,
        summary_rows=summary_rows,
        critical_rows=critical_rows,
        tracking_rows=tracking_rows,
        threshold_rows=threshold_rows,
        chi_fit_rows=chi_fit_rows,
        warnings=warnings_from_saved_audits(mode_rows, tracking_rows),
        plot_paths=plot_paths,
    )
    timing_path = write_csv(
        output_dir / "timing_report_plot_only.csv",
        [
            {
                "run_kind": "plot_only",
                "seconds": float("nan"),
                "cache_hits": 0,
                "cache_misses": 0,
                "solved_parameter_points": 0,
                "threshold_refinement_iterations": 0,
                "notes": f"timestamp={datetime.now().isoformat(timespec='seconds')}",
            }
        ],
        TIMING_FIELDS,
    )
    return [*plot_paths, timing_path], report


def thin_limit_low_frequency_deltas(
    *,
    beta_deg: float = DEFAULT_BETA_DEG,
    mu: float = 0.0,
    eta: float = DEFAULT_ETA,
    epsilon: float = 0.0025,
    n_roots: int = 4,
) -> list[float]:
    case = beta_workflow.CaseSpec(mu=float(mu), eta=float(eta), epsilon=float(epsilon))
    eb = beta_workflow.solve_model(case, float(beta_deg), int(n_roots), MODEL_EB)
    if eb.root_count_found < int(n_roots):
        eb = beta_workflow.retry_missing_roots(eb, case, float(beta_deg), int(n_roots), MODEL_EB)
    timo = beta_workflow.solve_model(case, float(beta_deg), int(n_roots), MODEL_TIMO)
    if timo.root_count_found < int(n_roots):
        timo = beta_workflow.retry_missing_roots(timo, case, float(beta_deg), int(n_roots), MODEL_TIMO)
    return [frequency_metrics(float(a), float(b))["delta_f"] for a, b in zip(eb.roots, timo.roots)]


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    start = time.perf_counter()
    if args.plot_only:
        paths, report = plot_only_outputs(args)
        elapsed = time.perf_counter() - start
        timing_row = {
            "run_kind": "plot_only",
            "seconds": float(elapsed),
            "cache_hits": 0,
            "cache_misses": 0,
            "solved_parameter_points": 0,
            "threshold_refinement_iterations": 0,
            "notes": f"timestamp={datetime.now().isoformat(timespec='seconds')}",
        }
        timing_path = write_csv(args.output_dir / "timing_report_plot_only.csv", [timing_row], TIMING_FIELDS)
        append_timing_row(args.output_dir / "timing_report.csv", timing_row)
        if timing_path not in paths:
            paths.append(timing_path)
        print("plot-only regeneration complete:")
        for path in [report, *paths]:
            print(f"  {rel(path)}")
        print(f"plot-only seconds: {elapsed:.3f}")
        return {"paths": paths, "report": report, "seconds": elapsed}

    root_cache = RootCache(args)
    products = solve_products(args, root_cache)
    elapsed = time.perf_counter() - start
    paths, report = write_outputs(
        args=args,
        products=products,
        run_seconds=elapsed,
        root_cache=root_cache,
        run_kind="smoke" if args.smoke else ("cached" if root_cache.misses == 0 else "full"),
    )
    print("generated EB validity vs Timoshenko stage 1 outputs:")
    for path in [report, *paths]:
        print(f"  {rel(path)}")
    print(f"solved parameter points: {products.solved_parameter_points}")
    print(f"cache hits/misses: {root_cache.hits}/{root_cache.misses}")
    print(f"threshold refinement iterations: {products.threshold_refinement_iterations}")
    print(f"warnings: {len(products.warnings)}")
    print(f"elapsed seconds: {elapsed:.3f}")
    return {
        "paths": paths,
        "report": report,
        "products": products,
        "cache_hits": root_cache.hits,
        "cache_misses": root_cache.misses,
        "seconds": elapsed,
    }


if __name__ == "__main__":
    main()
