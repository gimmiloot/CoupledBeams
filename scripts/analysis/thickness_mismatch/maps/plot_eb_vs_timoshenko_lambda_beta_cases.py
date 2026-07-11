from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from math import isfinite
from pathlib import Path
import sys
import time
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    find_roots_scan_bisect_eta,
    thickness_mismatch_factors,
)
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


DEFAULT_OUTPUT_DIR = Path("results") / "eb_vs_timoshenko_lambda_beta_cases"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_vs_timoshenko_lambda_beta_cases"
DEFAULT_CACHE_SUBDIR = "cache"

DEFAULT_MU_ETA_CASES = ((0.0, 0.0), (0.0, 0.1), (0.3, 0.1))
DEFAULT_EPSILON_VALUES = (0.0025, 0.05)
DEFAULT_BETA_MIN = 0.0
DEFAULT_BETA_MAX = 90.0
DEFAULT_BETA_STEP = 0.5
DEFAULT_REFINED_BETA_STEP = 0.05
DEFAULT_N_ROOTS = 6

ROOT_SCAN_START = 0.2
EB_LAMBDA_MAX = 35.0
EB_RETRY_LAMBDA_MAX = 50.0
EB_SCAN_STEP = 0.02
EB_RETRY_SCAN_STEP = 0.001
TIMO_SCAN_STEP = 0.01
TIMO_RETRY_SCAN_STEP = 0.001

DEFAULT_JUMP_ABS_THRESHOLD = 0.75
DEFAULT_JUMP_REL_THRESHOLD = 0.08
DEFAULT_REPAIR_WINDOW_DEG = 2.0
MAX_RECOVERY_INTERVAL_DEG = 10.0
SVD_RECOVERY_GAP_THRESHOLD = 0.75
SVD_RECOVERY_PREFILTER = 1.0e-4
SVD_RECOVERY_ACCEPT = 5.0e-6
CONTINUATION_SVD_ACCEPT = 1.0e-4
SVD_RECOVERY_MULTIPLICITY = 2.0e-5
SVD_RECOVERY_CLUSTER_TOL = 2.0e-3
CACHE_VERSION = "eb_vs_timo_lambda_beta_cases_v2"
TIMING_REPORT_NAME = "eb_vs_timo_lambda_beta_timing_report.csv"

MODEL_EB = "Euler-Bernoulli"
MODEL_TIMO = "Timoshenko"
MODELS = (MODEL_EB, MODEL_TIMO)

CASE_CSV_FIELDS = [
    "beta_deg",
    "sorted_index",
    "Lambda_EB_raw",
    "Lambda_Timoshenko_raw",
    "Lambda_EB_plot",
    "Lambda_Timoshenko_plot",
    "rel_diff_abs_Timoshenko_vs_EB",
    "suspicious_EB",
    "suspicious_Timoshenko",
    "retry_attempted",
    "retry_fixed",
    "plotted_as_nan",
    "notes",
    "status_EB",
    "status_Timoshenko",
    "root_warning_EB",
    "root_warning_Timoshenko",
    "root_count_EB",
    "root_count_Timoshenko",
    "retry_attempted_EB",
    "retry_attempted_Timoshenko",
    "retry_fixed_EB",
    "retry_fixed_Timoshenko",
    "plotted_as_nan_EB",
    "plotted_as_nan_Timoshenko",
    "notes_EB",
    "notes_Timoshenko",
]

SPIKE_AUDIT_FIELDS = [
    "case_index",
    "mu",
    "eta",
    "epsilon",
    "model",
    "sorted_index",
    "beta_deg",
    "Lambda_raw",
    "Lambda_plot",
    "jump_prev_raw",
    "jump_next_raw",
    "jump_rel_prev_raw",
    "jump_rel_next_raw",
    "suspicious_raw",
    "retry_attempted",
    "retry_fixed",
    "plotted_as_nan",
    "jump_prev_plot",
    "jump_next_plot",
    "jump_rel_prev_plot",
    "jump_rel_next_plot",
    "suspicious_plot",
    "notes",
]

SUMMARY_FIELDS = [
    "mu",
    "eta",
    "epsilon",
    "sorted_index",
    "max_rel_diff_over_beta",
    "mean_rel_diff_over_beta",
    "beta_at_max_rel_diff",
    "raw_suspicious_point_count",
    "suspicious_point_count",
    "retry_fixed_count",
    "nan_count_after_cleanup",
]

TIMING_FIELDS = [
    "run_id",
    "model",
    "mu",
    "eta",
    "epsilon",
    "n_beta_points",
    "n_roots",
    "root_mode",
    "cache_hit",
    "ordinary_compute_seconds",
    "spike_audit_seconds",
    "repair_seconds",
    "plotting_seconds",
    "total_seconds",
    "warnings",
    "fallback_count",
    "repair_count",
    "sv_recovery_calls",
    "notes",
]


@dataclass(frozen=True)
class CaseSpec:
    mu: float
    eta: float
    epsilon: float


@dataclass(frozen=True)
class RefinementWindow:
    mu: float
    eta: float
    epsilon: float
    beta_min: float
    beta_max: float
    reason: str


@dataclass(frozen=True)
class Args:
    beta_min: float
    beta_max: float
    beta_step: float
    refined_beta_step: float
    jump_abs_threshold: float
    jump_rel_threshold: float
    n_roots: int
    output_dir: Path
    cache_dir: Path
    reuse_cache: bool
    force_recompute: bool
    repair_spikes: bool
    repair_window_deg: float
    sv_recovery_only_on_spikes: bool
    use_known_timo_spike_windows: bool
    timo_root_mode: str
    plot_only: bool
    mu_eta_cases: tuple[tuple[float, float], ...]
    epsilon_values: tuple[float, ...]
    smoke: bool


@dataclass(frozen=True)
class RootResult:
    roots: tuple[float, ...]
    warnings: tuple[str, ...]
    root_count_found: int
    lambda_max_used: float
    scan_step_used: float
    retry_attempted: bool
    retry_changed_value: bool
    notes: tuple[str, ...]


@dataclass(frozen=True)
class TimingStats:
    run_id: str
    cache_hit: bool
    ordinary_compute_seconds: float
    spike_audit_seconds: float
    repair_seconds: float
    plotting_seconds: float
    total_seconds: float
    fallback_count: int
    repair_count: int
    sv_recovery_calls: int
    notes: tuple[str, ...]


@dataclass(frozen=True)
class RootTables:
    raw: dict[tuple[int, str], np.ndarray]
    clean: dict[tuple[int, str], np.ndarray]
    status: dict[tuple[int, str], np.ndarray]
    notes: dict[tuple[int, str], np.ndarray]
    retry_attempted: dict[tuple[int, str], np.ndarray]
    retry_fixed: dict[tuple[int, str], np.ndarray]
    plotted_as_nan: dict[tuple[int, str], np.ndarray]
    raw_suspicious: dict[tuple[int, str], np.ndarray]
    plot_suspicious: dict[tuple[int, str], np.ndarray]


FORCED_REFINEMENT_WINDOWS = (
    RefinementWindow(
        mu=0.0,
        eta=0.0,
        epsilon=0.05,
        beta_min=50.0,
        beta_max=60.0,
        reason="known Timoshenko missed-root spike around beta ~55-57 deg",
    ),
    RefinementWindow(
        mu=0.0,
        eta=0.0,
        epsilon=0.05,
        beta_min=80.0,
        beta_max=88.0,
        reason="known Timoshenko missed-root spike around beta ~85 deg",
    ),
    RefinementWindow(
        mu=0.0,
        eta=0.1,
        epsilon=0.05,
        beta_min=56.0,
        beta_max=66.0,
        reason="known Timoshenko spike window around beta ~60 deg",
    ),
)


def _repo_output_dir(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def _fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if not isfinite(value_f):
            return "nan"
        return f"{value_f:.16e}"
    return value


def _float_label(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def _rel(path: Path) -> str:
    return str(Path(path).relative_to(REPO_ROOT))


def _write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row.get(key, "")) for key in fields})


def default_cases(
    mu_eta_cases: Sequence[tuple[float, float]] = DEFAULT_MU_ETA_CASES,
    epsilon_values: Sequence[float] = DEFAULT_EPSILON_VALUES,
) -> tuple[CaseSpec, ...]:
    return tuple(
        CaseSpec(mu=float(mu), eta=float(eta), epsilon=float(epsilon))
        for mu, eta in mu_eta_cases
        for epsilon in epsilon_values
    )


def parse_mu_eta_case(value: str) -> tuple[float, float]:
    parts = [part.strip() for part in str(value).split(",")]
    if len(parts) != 2:
        raise argparse.ArgumentTypeError(f"mu-eta case must look like 'mu,eta', got {value!r}")
    try:
        return float(parts[0]), float(parts[1])
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"mu-eta case must contain numeric values, got {value!r}") from exc


def cases_are_default(mu_eta_cases: Sequence[tuple[float, float]], epsilon_values: Sequence[float]) -> bool:
    if len(mu_eta_cases) != len(DEFAULT_MU_ETA_CASES) or len(epsilon_values) != len(DEFAULT_EPSILON_VALUES):
        return False
    return all(
        abs(float(a_mu) - float(b_mu)) <= 1.0e-12 and abs(float(a_eta) - float(b_eta)) <= 1.0e-12
        for (a_mu, a_eta), (b_mu, b_eta) in zip(mu_eta_cases, DEFAULT_MU_ETA_CASES)
    ) and all(abs(float(a) - float(b)) <= 1.0e-12 for a, b in zip(epsilon_values, DEFAULT_EPSILON_VALUES))


def parse_args(argv: list[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic-only sorted in-plane Lambda(beta) comparison between "
            "Euler-Bernoulli and Timoshenko theories."
        )
    )
    parser.add_argument("--beta-min", type=float, default=DEFAULT_BETA_MIN)
    parser.add_argument("--beta-max", type=float, default=DEFAULT_BETA_MAX)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--refined-beta-step", type=float, default=DEFAULT_REFINED_BETA_STEP)
    parser.add_argument("--local-beta-step", type=float, default=None)
    parser.add_argument("--jump-abs-threshold", type=float, default=DEFAULT_JUMP_ABS_THRESHOLD)
    parser.add_argument("--jump-rel-threshold", type=float, default=DEFAULT_JUMP_REL_THRESHOLD)
    parser.add_argument("--n-roots", type=int, default=DEFAULT_N_ROOTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--repair-spikes", dest="repair_spikes", action="store_true", default=True)
    parser.add_argument("--no-repair-spikes", dest="repair_spikes", action="store_false")
    parser.add_argument("--repair-window-deg", type=float, default=DEFAULT_REPAIR_WINDOW_DEG)
    parser.add_argument("--sv-recovery-only-on-spikes", dest="sv_recovery_only_on_spikes", action="store_true", default=True)
    parser.add_argument("--no-sv-recovery-only-on-spikes", dest="sv_recovery_only_on_spikes", action="store_false")
    parser.add_argument("--use-known-timo-spike-windows", dest="use_known_timo_spike_windows", action="store_true", default=None)
    parser.add_argument("--no-use-known-timo-spike-windows", dest="use_known_timo_spike_windows", action="store_false")
    parser.add_argument("--timo-root-mode", choices=("global", "continuation"), default="continuation")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--epsilon-values", type=float, nargs="+", default=None)
    parser.add_argument("--mu-eta-cases", type=parse_mu_eta_case, nargs="+", default=None)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)

    beta_min = float(ns.beta_min)
    beta_max = float(ns.beta_max)
    beta_step = float(ns.beta_step)
    refined_beta_step = float(ns.local_beta_step if ns.local_beta_step is not None else ns.refined_beta_step)
    jump_abs_threshold = float(ns.jump_abs_threshold)
    jump_rel_threshold = float(ns.jump_rel_threshold)
    n_roots = int(ns.n_roots)
    output_dir = _repo_output_dir(Path(ns.output_dir))
    cache_dir = _repo_output_dir(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / DEFAULT_CACHE_SUBDIR
    mu_eta_cases = tuple(ns.mu_eta_cases) if ns.mu_eta_cases is not None else tuple(DEFAULT_MU_ETA_CASES)
    epsilon_values = tuple(float(value) for value in ns.epsilon_values) if ns.epsilon_values is not None else tuple(DEFAULT_EPSILON_VALUES)
    repair_spikes = bool(ns.repair_spikes)
    use_known_windows = (
        cases_are_default(mu_eta_cases, epsilon_values)
        if ns.use_known_timo_spike_windows is None
        else bool(ns.use_known_timo_spike_windows)
    )
    if ns.smoke:
        beta_min = 0.0
        beta_max = 90.0
        beta_step = 45.0
        refined_beta_step = 45.0
        n_roots = min(3, n_roots)
        output_dir = _repo_output_dir(SMOKE_OUTPUT_DIR)
        cache_dir = output_dir / DEFAULT_CACHE_SUBDIR
        repair_spikes = False

    args = Args(
        beta_min=beta_min,
        beta_max=beta_max,
        beta_step=beta_step,
        refined_beta_step=refined_beta_step,
        jump_abs_threshold=jump_abs_threshold,
        jump_rel_threshold=jump_rel_threshold,
        n_roots=n_roots,
        output_dir=output_dir,
        cache_dir=cache_dir,
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        repair_spikes=repair_spikes,
        repair_window_deg=float(ns.repair_window_deg),
        sv_recovery_only_on_spikes=bool(ns.sv_recovery_only_on_spikes),
        use_known_timo_spike_windows=bool(use_known_windows),
        timo_root_mode=str(ns.timo_root_mode),
        plot_only=bool(ns.plot_only),
        mu_eta_cases=mu_eta_cases,
        epsilon_values=epsilon_values,
        smoke=bool(ns.smoke),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if not (isfinite(args.beta_min) and isfinite(args.beta_max) and isfinite(args.beta_step)):
        raise ValueError("beta grid values must be finite.")
    if args.beta_step <= 0.0:
        raise ValueError("beta-step must be positive.")
    if args.refined_beta_step <= 0.0:
        raise ValueError("refined-beta-step must be positive.")
    if args.jump_abs_threshold <= 0.0 or args.jump_rel_threshold <= 0.0:
        raise ValueError("jump thresholds must be positive.")
    if args.repair_window_deg <= 0.0:
        raise ValueError("repair-window-deg must be positive.")
    if args.beta_max < args.beta_min:
        raise ValueError("beta-max must be greater than or equal to beta-min.")
    if args.n_roots <= 0:
        raise ValueError("n-roots must be positive.")
    if args.timo_root_mode not in {"global", "continuation"}:
        raise ValueError("timo-root-mode must be 'global' or 'continuation'.")
    if not args.mu_eta_cases:
        raise ValueError("at least one mu-eta case is required.")
    if not args.epsilon_values:
        raise ValueError("at least one epsilon value is required.")
    for case in default_cases(args.mu_eta_cases, args.epsilon_values):
        if case.epsilon <= 0.0:
            raise ValueError("epsilon values must be positive.")
        thickness_mismatch_factors(case.mu, case.eta)
        TIMO.tau_factors(case.mu, case.eta)


def regular_grid(start: float, end: float, step: float) -> np.ndarray:
    count = int(np.floor((float(end) - float(start)) / float(step) + 0.5)) + 1
    values = float(start) + float(step) * np.arange(count, dtype=float)
    if values.size == 0 or values[-1] < float(end) - 1.0e-10:
        values = np.append(values, float(end))
    values[0] = float(start)
    values[-1] = float(end)
    return np.unique(np.round(values, 12))


def case_matches_window(case: CaseSpec, window: RefinementWindow) -> bool:
    return (
        abs(case.mu - window.mu) <= 1.0e-12
        and abs(case.eta - window.eta) <= 1.0e-12
        and abs(case.epsilon - window.epsilon) <= 1.0e-12
    )


def refinement_windows_for_case(case: CaseSpec, args: Args) -> tuple[RefinementWindow, ...]:
    if args.smoke:
        return ()
    windows = [
        window
        for window in FORCED_REFINEMENT_WINDOWS
        if args.use_known_timo_spike_windows and case_matches_window(case, window)
    ]
    windows.append(
        RefinementWindow(
            mu=case.mu,
            eta=case.eta,
            epsilon=case.epsilon,
            beta_min=0.0,
            beta_max=2.5,
            reason="low-beta guard for jump-relative audit of the first sorted root",
        )
    )
    return tuple(windows)


def in_known_timo_spike_window(case: CaseSpec, beta_deg: float, args: Args) -> bool:
    if not args.use_known_timo_spike_windows:
        return False
    for window in FORCED_REFINEMENT_WINDOWS:
        if case_matches_window(case, window) and window.beta_min - 1.0e-12 <= float(beta_deg) <= window.beta_max + 1.0e-12:
            return True
    return False


def beta_grid_for_case(case: CaseSpec, args: Args) -> np.ndarray:
    values = list(regular_grid(args.beta_min, args.beta_max, args.beta_step))
    for window in refinement_windows_for_case(case, args):
        left = max(args.beta_min, window.beta_min)
        right = min(args.beta_max, window.beta_max)
        if right >= left:
            values.extend(float(value) for value in regular_grid(left, right, args.refined_beta_step))
    return np.unique(np.round(np.asarray(values, dtype=float), 12))


def beta_grids_by_case(cases: Sequence[CaseSpec], args: Args) -> list[np.ndarray]:
    return [beta_grid_for_case(case, args) for case in cases]


def is_base_beta_point(beta_deg: float, args: Args) -> bool:
    if args.beta_step <= 0.0:
        return False
    position = (float(beta_deg) - float(args.beta_min)) / float(args.beta_step)
    return abs(position - round(position)) <= 1.0e-8


def cache_settings(cases: Sequence[CaseSpec], beta_values_by_case: Sequence[np.ndarray], args: Args) -> dict[str, object]:
    return {
        "cache_version": CACHE_VERSION,
        "cases": [{"mu": case.mu, "eta": case.eta, "epsilon": case.epsilon} for case in cases],
        "beta_grids": [[float(value) for value in beta_values] for beta_values in beta_values_by_case],
        "beta_min": args.beta_min,
        "beta_max": args.beta_max,
        "beta_step": args.beta_step,
        "refined_beta_step": args.refined_beta_step,
        "jump_abs_threshold": args.jump_abs_threshold,
        "jump_rel_threshold": args.jump_rel_threshold,
        "n_roots": args.n_roots,
        "repair_spikes": args.repair_spikes,
        "repair_window_deg": args.repair_window_deg,
        "sv_recovery_only_on_spikes": args.sv_recovery_only_on_spikes,
        "use_known_timo_spike_windows": args.use_known_timo_spike_windows,
        "timo_root_mode": args.timo_root_mode,
        "continuation_uses_base_grid_global_anchors": True,
        "continuation_disabled_in_known_timo_spike_windows": True,
        "mu_eta_cases": [[float(mu), float(eta)] for mu, eta in args.mu_eta_cases],
        "epsilon_values": [float(value) for value in args.epsilon_values],
        "root_scan_start": ROOT_SCAN_START,
        "eb_lambda_max": EB_LAMBDA_MAX,
        "eb_retry_lambda_max": EB_RETRY_LAMBDA_MAX,
        "eb_scan_step": EB_SCAN_STEP,
        "eb_retry_scan_step": EB_RETRY_SCAN_STEP,
        "timo_scan_step": TIMO_SCAN_STEP,
        "timo_retry_scan_step": TIMO_RETRY_SCAN_STEP,
        "sv_recovery_gap_threshold": SVD_RECOVERY_GAP_THRESHOLD,
        "sv_recovery_prefilter": SVD_RECOVERY_PREFILTER,
        "sv_recovery_accept": SVD_RECOVERY_ACCEPT,
        "continuation_svd_accept": CONTINUATION_SVD_ACCEPT,
        "sv_recovery_multiplicity": SVD_RECOVERY_MULTIPLICITY,
        "sv_recovery_cluster_tol": SVD_RECOVERY_CLUSTER_TOL,
    }


def cache_settings_json(cases: Sequence[CaseSpec], beta_values_by_case: Sequence[np.ndarray], args: Args) -> str:
    return json.dumps(cache_settings(cases, beta_values_by_case, args), sort_keys=True, separators=(",", ":"))


def cache_file_path(cases: Sequence[CaseSpec], beta_values_by_case: Sequence[np.ndarray], args: Args) -> Path:
    digest = hashlib.sha256(cache_settings_json(cases, beta_values_by_case, args).encode("utf-8")).hexdigest()[:16]
    return args.cache_dir / f"roots_all_cases_{CACHE_VERSION}_{digest}.npz"


def root_map_arrays(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
) -> dict[str, np.ndarray]:
    n_cases = len(cases)
    max_beta = max(len(values) for values in beta_values_by_case) if beta_values_by_case else 0
    n_models = len(MODELS)
    roots = np.full((n_cases, n_models, max_beta, args.n_roots), np.nan, dtype=float)
    root_count = np.zeros((n_cases, n_models, max_beta), dtype=int)
    lambda_max = np.full((n_cases, n_models, max_beta), np.nan, dtype=float)
    scan_step = np.full((n_cases, n_models, max_beta), np.nan, dtype=float)
    retry_attempted = np.full((n_cases, n_models, max_beta), False, dtype=bool)
    retry_changed = np.full((n_cases, n_models, max_beta), False, dtype=bool)
    warnings = np.empty((n_cases, n_models, max_beta), dtype=object)
    notes = np.empty((n_cases, n_models, max_beta), dtype=object)
    warnings.fill("")
    notes.fill("")
    for case_index, beta_values in enumerate(beta_values_by_case):
        for model_index, model in enumerate(MODELS):
            for beta_index in range(len(beta_values)):
                result = root_map[(case_index, model, beta_index)]
                roots[case_index, model_index, beta_index, :] = np.asarray(result.roots, dtype=float)
                root_count[case_index, model_index, beta_index] = int(result.root_count_found)
                lambda_max[case_index, model_index, beta_index] = float(result.lambda_max_used)
                scan_step[case_index, model_index, beta_index] = float(result.scan_step_used)
                retry_attempted[case_index, model_index, beta_index] = bool(result.retry_attempted)
                retry_changed[case_index, model_index, beta_index] = bool(result.retry_changed_value)
                warnings[case_index, model_index, beta_index] = "\n".join(result.warnings)
                notes[case_index, model_index, beta_index] = "\n".join(result.notes)
    return {
        "root_roots": roots,
        "root_count": root_count,
        "root_lambda_max": lambda_max,
        "root_scan_step": scan_step,
        "root_retry_attempted": retry_attempted,
        "root_retry_changed": retry_changed,
        "root_warnings": warnings,
        "root_notes": notes,
    }


def tables_arrays(cases: Sequence[CaseSpec], beta_values_by_case: Sequence[np.ndarray], args: Args, tables: RootTables) -> dict[str, np.ndarray]:
    n_cases = len(cases)
    max_beta = max(len(values) for values in beta_values_by_case) if beta_values_by_case else 0
    n_models = len(MODELS)
    shape = (n_cases, n_models, args.n_roots, max_beta)
    arrays: dict[str, np.ndarray] = {
        "table_raw": np.full(shape, np.nan, dtype=float),
        "table_clean": np.full(shape, np.nan, dtype=float),
        "table_status": np.empty(shape, dtype=object),
        "table_notes": np.empty(shape, dtype=object),
        "table_retry_attempted": np.full(shape, False, dtype=bool),
        "table_retry_fixed": np.full(shape, False, dtype=bool),
        "table_plotted_as_nan": np.full(shape, False, dtype=bool),
        "table_raw_suspicious": np.full(shape, False, dtype=bool),
        "table_plot_suspicious": np.full(shape, False, dtype=bool),
    }
    arrays["table_status"].fill("")
    arrays["table_notes"].fill("")
    table_map = {
        "table_raw": tables.raw,
        "table_clean": tables.clean,
        "table_status": tables.status,
        "table_notes": tables.notes,
        "table_retry_attempted": tables.retry_attempted,
        "table_retry_fixed": tables.retry_fixed,
        "table_plotted_as_nan": tables.plotted_as_nan,
        "table_raw_suspicious": tables.raw_suspicious,
        "table_plot_suspicious": tables.plot_suspicious,
    }
    for case_index, beta_values in enumerate(beta_values_by_case):
        stop = len(beta_values)
        for model_index, model in enumerate(MODELS):
            key = (case_index, model)
            for array_name, table_dict in table_map.items():
                arrays[array_name][case_index, model_index, :, :stop] = table_dict[key]
    return arrays


def save_cache(
    path: Path,
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    tables: RootTables,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, object] = {
        "settings_json": np.array(cache_settings_json(cases, beta_values_by_case, args), dtype=object),
        "timestamp_utc": np.array(datetime.now(timezone.utc).isoformat(), dtype=object),
        "models": np.array(MODELS, dtype=object),
        "beta_lengths": np.array([len(values) for values in beta_values_by_case], dtype=int),
    }
    max_beta = max(len(values) for values in beta_values_by_case) if beta_values_by_case else 0
    beta_grid = np.full((len(cases), max_beta), np.nan, dtype=float)
    for case_index, beta_values in enumerate(beta_values_by_case):
        beta_grid[case_index, : len(beta_values)] = beta_values
    payload["beta_grid"] = beta_grid
    payload.update(root_map_arrays(cases, beta_values_by_case, args, root_map))
    payload.update(tables_arrays(cases, beta_values_by_case, args, tables))
    np.savez_compressed(path, **payload)


def load_cache(
    path: Path,
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
) -> tuple[dict[tuple[int, str, int], RootResult], RootTables] | None:
    if not path.exists():
        return None
    try:
        data = np.load(path, allow_pickle=True)
    except (OSError, ValueError):
        return None
    try:
        if str(data["settings_json"].item()) != cache_settings_json(cases, beta_values_by_case, args):
            return None
        beta_lengths = np.asarray(data["beta_lengths"], dtype=int)
        root_map: dict[tuple[int, str, int], RootResult] = {}
        for case_index, beta_values in enumerate(beta_values_by_case):
            if int(beta_lengths[case_index]) != len(beta_values):
                return None
            for model_index, model in enumerate(MODELS):
                for beta_index in range(len(beta_values)):
                    warnings_text = str(data["root_warnings"][case_index, model_index, beta_index])
                    notes_text = str(data["root_notes"][case_index, model_index, beta_index])
                    root_map[(case_index, model, beta_index)] = RootResult(
                        roots=tuple(float(value) for value in data["root_roots"][case_index, model_index, beta_index, : args.n_roots]),
                        warnings=tuple(item for item in warnings_text.split("\n") if item),
                        root_count_found=int(data["root_count"][case_index, model_index, beta_index]),
                        lambda_max_used=float(data["root_lambda_max"][case_index, model_index, beta_index]),
                        scan_step_used=float(data["root_scan_step"][case_index, model_index, beta_index]),
                        retry_attempted=bool(data["root_retry_attempted"][case_index, model_index, beta_index]),
                        retry_changed_value=bool(data["root_retry_changed"][case_index, model_index, beta_index]),
                        notes=tuple(item for item in notes_text.split("\n") if item) or ("ok",),
                    )

        def table_dict(name: str) -> dict[tuple[int, str], np.ndarray]:
            out: dict[tuple[int, str], np.ndarray] = {}
            source = data[name]
            for case_index, beta_values in enumerate(beta_values_by_case):
                stop = len(beta_values)
                for model_index, model in enumerate(MODELS):
                    out[(case_index, model)] = np.array(source[case_index, model_index, :, :stop], copy=True)
            return out

        tables = RootTables(
            raw=table_dict("table_raw"),
            clean=table_dict("table_clean"),
            status=table_dict("table_status"),
            notes=table_dict("table_notes"),
            retry_attempted=table_dict("table_retry_attempted"),
            retry_fixed=table_dict("table_retry_fixed"),
            plotted_as_nan=table_dict("table_plotted_as_nan"),
            raw_suspicious=table_dict("table_raw_suspicious"),
            plot_suspicious=table_dict("table_plot_suspicious"),
        )
    finally:
        data.close()
    return root_map, tables


def normalize_roots(raw_roots: Sequence[float], n_roots: int) -> tuple[tuple[float, ...], list[str], int]:
    notes: list[str] = []
    raw = [float(root) for root in raw_roots]
    finite_roots = [root for root in raw if isfinite(root)]
    if len(finite_roots) != len(raw):
        notes.append("non-finite roots removed")
    if any(finite_roots[index + 1] < finite_roots[index] for index in range(len(finite_roots) - 1)):
        finite_roots = sorted(finite_roots)
        notes.append("roots sorted after solve")
    root_count = min(len(finite_roots), int(n_roots))
    roots = finite_roots[: int(n_roots)]
    if len(roots) < int(n_roots):
        roots.extend([float("nan")] * (int(n_roots) - len(roots)))
        notes.append("missing roots filled with NaN")
    if not notes:
        notes.append("ok")
    return tuple(roots), notes, root_count


def roots_close(left: Sequence[float], right: Sequence[float]) -> bool:
    if len(left) != len(right):
        return False
    for a, b in zip(left, right):
        a_f = float(a)
        b_f = float(b)
        if not (isfinite(a_f) or isfinite(b_f)):
            if isfinite(a_f) != isfinite(b_f):
                return False
            continue
        scale = max(abs(a_f), abs(b_f), 1.0)
        if abs(a_f - b_f) > 1.0e-8 + 1.0e-8 * scale:
            return False
    return True


def row_normalized_matrix(matrix: np.ndarray) -> np.ndarray:
    out = np.array(matrix, dtype=float, copy=True)
    row_norms = np.linalg.norm(out, axis=1)
    row_norms[row_norms == 0.0] = 1.0
    return out / row_norms[:, None]


def model_singular_values(case: CaseSpec, beta_deg: float, lambda_value: float, model: str) -> tuple[float, float]:
    try:
        if model == MODEL_EB:
            matrix = assemble_clamped_coupled_matrix_eta(
                float(lambda_value),
                float(np.deg2rad(beta_deg)),
                case.mu,
                case.epsilon,
                case.eta,
            )
            scaled = row_normalized_matrix(matrix)
        elif model == MODEL_TIMO:
            matrix, _warnings = TIMO.timo_coupling_matrix(
                float(lambda_value),
                float(beta_deg),
                case.mu,
                case.epsilon,
                case.eta,
            )
            scaled = TIMO.row_normalized(matrix)
        else:
            raise ValueError(f"unknown model: {model}")
        if not np.all(np.isfinite(scaled)):
            return float("inf"), float("inf")
        singular_values = np.linalg.svd(scaled, compute_uv=False)
    except (FloatingPointError, ValueError, OverflowError, np.linalg.LinAlgError):
        return float("inf"), float("inf")
    if len(singular_values) < 2:
        return float("inf"), float("inf")
    return float(singular_values[-1]), float(singular_values[-2])


def refine_svd_minimum(case: CaseSpec, beta_deg: float, model: str, left: float, right: float) -> tuple[float, float, float]:
    a = float(left)
    b = float(right)
    if not (isfinite(a) and isfinite(b)) or b <= a:
        middle = 0.5 * (a + b)
        s_min, s_second = model_singular_values(case, beta_deg, middle, model)
        return middle, s_min, s_second

    inv_phi = (np.sqrt(5.0) - 1.0) / 2.0
    inv_phi_sq = (3.0 - np.sqrt(5.0)) / 2.0
    c = a + inv_phi_sq * (b - a)
    d = a + inv_phi * (b - a)
    fc, _ = model_singular_values(case, beta_deg, c, model)
    fd, _ = model_singular_values(case, beta_deg, d, model)
    for _iteration in range(44):
        if fc <= fd:
            b = d
            d = c
            fd = fc
            c = a + inv_phi_sq * (b - a)
            fc, _ = model_singular_values(case, beta_deg, c, model)
        else:
            a = c
            c = d
            fc = fd
            d = a + inv_phi * (b - a)
            fd, _ = model_singular_values(case, beta_deg, d, model)
    root = 0.5 * (a + b)
    s_min, s_second = model_singular_values(case, beta_deg, root, model)
    return root, s_min, s_second


def svd_candidates_in_interval(
    case: CaseSpec,
    beta_deg: float,
    model: str,
    left: float,
    right: float,
    scan_step: float,
) -> list[tuple[float, int, float]]:
    margin = max(1.0e-5, 0.25 * float(scan_step))
    scan_left = float(left) + margin
    scan_right = float(right) - margin
    if scan_right <= scan_left + float(scan_step):
        return []

    grid = regular_grid(scan_left, scan_right, float(scan_step))
    if grid.size < 3:
        return []

    sigma = np.array(
        [model_singular_values(case, beta_deg, float(value), model)[0] for value in grid],
        dtype=float,
    )
    raw_candidates: list[tuple[float, int, float, float]] = []
    for index in range(1, grid.size - 1):
        current = float(sigma[index])
        if not isfinite(current) or current > SVD_RECOVERY_PREFILTER:
            continue
        if current <= float(sigma[index - 1]) and current <= float(sigma[index + 1]):
            root, s_min, s_second = refine_svd_minimum(
                case,
                beta_deg,
                model,
                float(grid[index - 1]),
                float(grid[index + 1]),
            )
            if s_min <= SVD_RECOVERY_ACCEPT:
                multiplicity = 2 if s_second <= SVD_RECOVERY_MULTIPLICITY else 1
                raw_candidates.append((root, multiplicity, s_min, s_second))

    if not raw_candidates:
        return []

    raw_candidates.sort(key=lambda item: item[0])
    grouped: list[list[tuple[float, int, float, float]]] = []
    for candidate in raw_candidates:
        if not grouped or abs(candidate[0] - grouped[-1][-1][0]) > SVD_RECOVERY_CLUSTER_TOL:
            grouped.append([candidate])
        else:
            grouped[-1].append(candidate)

    candidates: list[tuple[float, int, float]] = []
    for group in grouped:
        best = min(group, key=lambda item: item[2])
        roots = [item[0] for item in group]
        multiplicity = max(len(group), max(item[1] for item in group))
        for root in roots[:multiplicity]:
            candidates.append((float(root), 1, float(best[2])))
        for _extra in range(max(0, multiplicity - len(roots))):
            candidates.append((float(best[0]), 1, float(best[2])))
    return candidates


def augment_roots_with_svd_recovery(
    raw_roots: Sequence[float],
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    model: str,
    scan_step: float,
    target_indices: Sequence[int] | None = None,
) -> tuple[list[float], list[str]]:
    finite_roots = sorted(float(root) for root in raw_roots if isfinite(float(root)))
    if len(finite_roots) < 2:
        return list(raw_roots), []

    additions: list[float] = []
    inspected_roots = finite_roots[: max(int(n_roots), 2)]
    intervals: list[tuple[float, float]] = []
    if target_indices is None:
        if inspected_roots[0] - ROOT_SCAN_START > SVD_RECOVERY_GAP_THRESHOLD:
            intervals.append((ROOT_SCAN_START, inspected_roots[0]))
        intervals.extend(zip(inspected_roots, inspected_roots[1:]))
    else:
        for target_index in sorted({int(index) for index in target_indices if int(index) >= 0}):
            if target_index == 0 and inspected_roots[0] - ROOT_SCAN_START > SVD_RECOVERY_GAP_THRESHOLD:
                intervals.append((ROOT_SCAN_START, inspected_roots[0]))
            if 0 < target_index < len(inspected_roots):
                intervals.append((inspected_roots[target_index - 1], inspected_roots[target_index]))
            if target_index + 1 < len(inspected_roots):
                intervals.append((inspected_roots[target_index], inspected_roots[target_index + 1]))

    unique_intervals: list[tuple[float, float]] = []
    for left, right in intervals:
        key_interval = (round(float(left), 12), round(float(right), 12))
        if key_interval not in [(round(a, 12), round(b, 12)) for a, b in unique_intervals]:
            unique_intervals.append((float(left), float(right)))

    for left, right in unique_intervals:
        if right - left <= SVD_RECOVERY_GAP_THRESHOLD:
            continue
        candidates = svd_candidates_in_interval(case, beta_deg, model, left, right, scan_step)
        for root, _multiplicity, _residual in candidates:
            if any(abs(root - existing) <= SVD_RECOVERY_CLUSTER_TOL for existing in finite_roots):
                continue
            additions.append(float(root))

    if not additions:
        return list(raw_roots), []

    augmented = sorted(finite_roots + additions)
    note = (
        "SVD recovery added "
        f"{len(additions)} candidate root(s) at beta={float(beta_deg):.8g} deg for {model}"
    )
    return augmented, [note]


def timo_initial_upper(mu: float, eta: float, n_roots: int) -> float:
    factors = TIMO.tau_factors(float(mu), float(eta))
    l1, l2 = TIMO.segment_lengths(factors.mu)
    rod_scales = [
        float(np.sqrt(factors.tau1)) / max(l1, 0.08),
        float(np.sqrt(factors.tau2)) / max(l2, 0.08),
    ]
    return max(18.0, 1.35 * np.pi * (float(n_roots) + 4.0) * max(rod_scales))


def solve_eb(
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    *,
    retry: bool = False,
    sv_recovery: bool = False,
    sv_target_indices: Sequence[int] | None = None,
) -> RootResult:
    lambda_max = EB_RETRY_LAMBDA_MAX if retry else EB_LAMBDA_MAX
    scan_step = EB_RETRY_SCAN_STEP if retry else EB_SCAN_STEP
    warnings: list[str] = []
    try:
        raw = find_roots_scan_bisect_eta(
            beta=float(np.deg2rad(beta_deg)),
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            n_roots=int(n_roots),
            Lmin=ROOT_SCAN_START,
            Lmax=float(lambda_max),
            scan_step=float(scan_step),
        )
    except (FloatingPointError, ValueError, OverflowError) as exc:
        raw = [float("nan")] * int(n_roots)
        warnings.append(f"EB solve failed: {exc}")
    extra_notes: list[str] = []
    if retry and sv_recovery:
        raw, extra_notes = augment_roots_with_svd_recovery(
            raw,
            case,
            beta_deg,
            int(n_roots),
            MODEL_EB,
            float(scan_step),
            target_indices=sv_target_indices,
        )
    roots, notes, root_count = normalize_roots(raw, int(n_roots))
    notes.extend(extra_notes)
    if root_count < int(n_roots):
        warnings.append(f"EB found {root_count} roots, fewer than requested {int(n_roots)}.")
    return RootResult(
        roots=roots,
        warnings=tuple(warnings),
        root_count_found=root_count,
        lambda_max_used=float(lambda_max),
        scan_step_used=float(scan_step),
        retry_attempted=retry,
        retry_changed_value=False,
        notes=tuple(notes),
    )


def solve_timo(
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    *,
    retry: bool = False,
    upper_hint: float | None = None,
    sv_recovery: bool = False,
    sv_target_indices: Sequence[int] | None = None,
) -> RootResult:
    scan_step = TIMO_RETRY_SCAN_STEP if retry else TIMO_SCAN_STEP
    fallback_upper = timo_initial_upper(case.mu, case.eta, int(n_roots))
    if upper_hint is not None and isfinite(float(upper_hint)):
        hinted_upper = max(12.0, 1.35 * float(upper_hint))
    else:
        hinted_upper = fallback_upper
    upper = max(fallback_upper, hinted_upper) if retry else hinted_upper
    warnings: list[str] = []
    try:
        raw, root_warnings = TIMO.find_roots_by_sign_scan(
            lambda value: TIMO.timo_det_for_scan(
                float(value),
                float(beta_deg),
                case.mu,
                case.epsilon,
                case.eta,
            ),
            int(n_roots),
            start=TIMO.ROOT_SCAN_START,
            upper=float(upper),
            scan_step=float(scan_step),
            grow_factor=1.35,
            max_tries=8,
        )
        warnings.extend(root_warnings)
    except (FloatingPointError, ValueError, OverflowError) as exc:
        raw = [float("nan")] * int(n_roots)
        warnings.append(f"Timoshenko solve failed: {exc}")
    extra_notes: list[str] = []
    if retry and sv_recovery:
        raw, extra_notes = augment_roots_with_svd_recovery(
            raw,
            case,
            beta_deg,
            int(n_roots),
            MODEL_TIMO,
            float(scan_step),
            target_indices=sv_target_indices,
        )
    roots, notes, root_count = normalize_roots(raw, int(n_roots))
    notes.extend(extra_notes)
    if root_count < int(n_roots):
        warnings.append(f"Timoshenko found {root_count} roots, fewer than requested {int(n_roots)}.")
    return RootResult(
        roots=roots,
        warnings=tuple(warnings),
        root_count_found=root_count,
        lambda_max_used=float(upper),
        scan_step_used=float(scan_step),
        retry_attempted=retry,
        retry_changed_value=False,
        notes=tuple(notes),
    )


def local_eb_root_candidates(
    case: CaseSpec,
    beta_deg: float,
    previous_roots: Sequence[float],
    n_roots: int,
) -> list[float]:
    previous = [float(root) for root in previous_roots[: int(n_roots)] if isfinite(float(root))]
    if len(previous) < int(n_roots):
        return []
    previous = sorted(previous[: int(n_roots)])
    candidates: list[float] = []
    beta_rad = float(np.deg2rad(beta_deg))
    for index, root in enumerate(previous):
        if index == 0:
            right = 0.5 * (previous[0] + previous[1])
            left = max(ROOT_SCAN_START, root - max(0.35, right - root))
        elif index == len(previous) - 1:
            left = 0.5 * (previous[-2] + previous[-1])
            right = root + max(0.35, root - left)
        else:
            left = 0.5 * (previous[index - 1] + root)
            right = 0.5 * (root + previous[index + 1])
        left = max(ROOT_SCAN_START, float(left) - 0.03)
        right = float(right) + 0.03
        if right <= left:
            continue
        interval_step = min(EB_RETRY_SCAN_STEP, max((right - left) / 18.0, 0.0005))
        roots = find_roots_scan_bisect_eta(
            beta=beta_rad,
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            n_roots=1,
            Lmin=float(left),
            Lmax=float(right),
            scan_step=float(interval_step),
        )
        added_finite_candidate = False
        for candidate in roots:
            candidate_f = float(candidate)
            if not isfinite(candidate_f):
                continue
            if any(abs(candidate_f - existing) <= 1.0e-8 for existing in candidates):
                continue
            candidates.append(candidate_f)
            added_finite_candidate = True
        if not added_finite_candidate:
            sv_root, sv_min, _sv_second = refine_svd_minimum(case, beta_deg, MODEL_EB, float(left), float(right))
            if sv_min <= CONTINUATION_SVD_ACCEPT and not any(abs(sv_root - existing) <= 1.0e-8 for existing in candidates):
                candidates.append(float(sv_root))
    return sorted(candidates)


def solve_eb_continuation(
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    previous_roots: Sequence[float] | None,
) -> RootResult:
    if previous_roots is None:
        return solve_eb(case, beta_deg, n_roots, retry=True, sv_recovery=False)
    local_roots = local_eb_root_candidates(case, beta_deg, previous_roots, n_roots)
    if len(local_roots) >= int(n_roots):
        roots, notes, root_count = normalize_roots(local_roots, int(n_roots))
        notes.append("EB continuation local brackets")
        return RootResult(
            roots=roots,
            warnings=(),
            root_count_found=root_count,
            lambda_max_used=max(local_roots[: int(n_roots)]),
            scan_step_used=EB_RETRY_SCAN_STEP,
            retry_attempted=True,
            retry_changed_value=False,
            notes=tuple(notes),
        )
    fallback = solve_eb(case, beta_deg, n_roots, retry=True, sv_recovery=False)
    notes = list(fallback.notes)
    notes.append(f"EB continuation fallback to dense scan; local found {len(local_roots)} roots")
    return RootResult(
        roots=fallback.roots,
        warnings=fallback.warnings,
        root_count_found=fallback.root_count_found,
        lambda_max_used=fallback.lambda_max_used,
        scan_step_used=fallback.scan_step_used,
        retry_attempted=True,
        retry_changed_value=fallback.retry_changed_value,
        notes=tuple(notes),
    )


def local_timo_root_candidates(
    case: CaseSpec,
    beta_deg: float,
    previous_roots: Sequence[float],
    n_roots: int,
) -> list[float]:
    previous = [float(root) for root in previous_roots[: int(n_roots)] if isfinite(float(root))]
    if len(previous) < int(n_roots):
        return []
    previous = sorted(previous[: int(n_roots)])
    candidates: list[float] = []
    func = lambda value: TIMO.timo_det_for_scan(
        float(value),
        float(beta_deg),
        case.mu,
        case.epsilon,
        case.eta,
    )
    for index, root in enumerate(previous):
        if index == 0:
            right = 0.5 * (previous[0] + previous[1])
            left = max(ROOT_SCAN_START, root - max(0.35, right - root))
        elif index == len(previous) - 1:
            left = 0.5 * (previous[-2] + previous[-1])
            right = root + max(0.35, root - left)
        else:
            left = 0.5 * (previous[index - 1] + root)
            right = 0.5 * (root + previous[index + 1])

        left = max(ROOT_SCAN_START, float(left) - 0.03)
        right = float(right) + 0.03
        if right <= left:
            continue
        interval_step = min(TIMO_SCAN_STEP, max((right - left) / 18.0, 0.0005))
        roots, _warnings = TIMO.find_roots_by_sign_scan(
            func,
            1,
            start=float(left),
            upper=float(right),
            scan_step=float(interval_step),
            grow_factor=1.0,
            max_tries=1,
        )
        added_finite_candidate = False
        for candidate in roots:
            candidate_f = float(candidate)
            if not isfinite(candidate_f):
                continue
            if any(abs(candidate_f - existing) <= 1.0e-8 for existing in candidates):
                continue
            candidates.append(candidate_f)
            added_finite_candidate = True
        if not added_finite_candidate:
            sv_root, sv_min, _sv_second = refine_svd_minimum(case, beta_deg, MODEL_TIMO, float(left), float(right))
            if sv_min <= CONTINUATION_SVD_ACCEPT and not any(abs(sv_root - existing) <= 1.0e-8 for existing in candidates):
                candidates.append(float(sv_root))
    return sorted(candidates)


def solve_timo_continuation(
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    previous_roots: Sequence[float] | None,
    *,
    upper_hint: float | None = None,
) -> RootResult:
    if previous_roots is None:
        return solve_timo(case, beta_deg, n_roots, upper_hint=upper_hint)

    local_roots = local_timo_root_candidates(case, beta_deg, previous_roots, n_roots)
    if len(local_roots) >= int(n_roots):
        roots, notes, root_count = normalize_roots(local_roots, int(n_roots))
        notes.append("continuation local brackets")
        return RootResult(
            roots=roots,
            warnings=(),
            root_count_found=root_count,
            lambda_max_used=max(local_roots[: int(n_roots)]),
            scan_step_used=TIMO_SCAN_STEP,
            retry_attempted=False,
            retry_changed_value=False,
            notes=tuple(notes),
        )

    fallback = solve_timo(case, beta_deg, n_roots, upper_hint=upper_hint)
    notes = list(fallback.notes)
    notes.append(f"continuation fallback to global scan; local found {len(local_roots)} roots")
    return RootResult(
        roots=fallback.roots,
        warnings=fallback.warnings,
        root_count_found=fallback.root_count_found,
        lambda_max_used=fallback.lambda_max_used,
        scan_step_used=fallback.scan_step_used,
        retry_attempted=fallback.retry_attempted,
        retry_changed_value=fallback.retry_changed_value,
        notes=tuple(notes),
    )


def solve_model(
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    model: str,
    *,
    retry: bool = False,
    upper_hint: float | None = None,
    sv_recovery: bool = False,
    sv_target_indices: Sequence[int] | None = None,
) -> RootResult:
    if model == MODEL_EB:
        return solve_eb(
            case,
            beta_deg,
            n_roots,
            retry=retry,
            sv_recovery=sv_recovery,
            sv_target_indices=sv_target_indices,
        )
    if model == MODEL_TIMO:
        return solve_timo(
            case,
            beta_deg,
            n_roots,
            retry=retry,
            upper_hint=upper_hint,
            sv_recovery=sv_recovery,
            sv_target_indices=sv_target_indices,
        )
    raise ValueError(f"unknown model: {model}")


def retry_missing_roots(
    initial: RootResult,
    case: CaseSpec,
    beta_deg: float,
    n_roots: int,
    model: str,
    *,
    upper_hint: float | None = None,
) -> RootResult:
    retry = solve_model(case, beta_deg, n_roots, model, retry=True, upper_hint=upper_hint)
    changed = not roots_close(initial.roots, retry.roots)
    chosen = retry if retry.root_count_found > initial.root_count_found else initial
    notes = list(chosen.notes)
    notes.append("retry used for missing roots" if chosen is retry else "retry attempted for missing roots; initial kept")
    return RootResult(
        roots=chosen.roots,
        warnings=chosen.warnings,
        root_count_found=chosen.root_count_found,
        lambda_max_used=chosen.lambda_max_used,
        scan_step_used=chosen.scan_step_used,
        retry_attempted=True,
        retry_changed_value=changed,
        notes=tuple(notes),
    )


def compute_root_map(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
) -> dict[tuple[int, str, int], RootResult]:
    root_map: dict[tuple[int, str, int], RootResult] = {}
    for case_index, case in enumerate(cases):
        beta_values = np.asarray(beta_values_by_case[case_index], dtype=float)
        print(
            "solving case "
            f"{case_index + 1}/{len(cases)}: mu={case.mu:g}, eta={case.eta:g}, "
            f"epsilon={case.epsilon:g}, beta points={len(beta_values)}"
        )
        previous_timo_roots: tuple[float, ...] | None = None
        for beta_index, beta_deg in enumerate(beta_values):
            eb_result = solve_model(case, float(beta_deg), args.n_roots, MODEL_EB)
            if eb_result.root_count_found < args.n_roots:
                eb_result = retry_missing_roots(
                    eb_result,
                    case,
                    float(beta_deg),
                    args.n_roots,
                    MODEL_EB,
                )
            root_map[(case_index, MODEL_EB, beta_index)] = eb_result

            finite_eb_roots = [float(root) for root in eb_result.roots if isfinite(float(root))]
            upper_hint = max(finite_eb_roots) if finite_eb_roots else None
            use_global_timo = (
                args.timo_root_mode == "global"
                or is_base_beta_point(float(beta_deg), args)
                or in_known_timo_spike_window(case, float(beta_deg), args)
            )
            if not use_global_timo:
                timo_result = solve_timo_continuation(
                    case,
                    float(beta_deg),
                    args.n_roots,
                    previous_timo_roots,
                    upper_hint=upper_hint,
                )
            else:
                timo_result = solve_model(
                    case,
                    float(beta_deg),
                    args.n_roots,
                    MODEL_TIMO,
                    upper_hint=upper_hint,
                )
            if timo_result.root_count_found < args.n_roots:
                timo_result = retry_missing_roots(
                    timo_result,
                    case,
                    float(beta_deg),
                    args.n_roots,
                    MODEL_TIMO,
                    upper_hint=upper_hint,
                )
            root_map[(case_index, MODEL_TIMO, beta_index)] = timo_result
            previous_timo_roots = timo_result.roots
    return root_map


def jump_abs_rel(left: float, right: float) -> tuple[float, float]:
    if not (isfinite(float(left)) and isfinite(float(right))):
        return float("nan"), float("nan")
    jump_abs = abs(float(right) - float(left))
    scale = max(abs(float(left)), abs(float(right)), 1.0)
    return jump_abs, jump_abs / scale


def row_jump_metrics(values: np.ndarray, sorted_zero: int, beta_index: int) -> tuple[float, float, float, float]:
    row = np.asarray(values[sorted_zero, :], dtype=float)
    current = float(row[beta_index])
    if beta_index > 0:
        jump_prev, rel_prev = jump_abs_rel(float(row[beta_index - 1]), current)
    else:
        jump_prev, rel_prev = float("nan"), float("nan")
    if beta_index + 1 < row.size:
        jump_next, rel_next = jump_abs_rel(current, float(row[beta_index + 1]))
    else:
        jump_next, rel_next = float("nan"), float("nan")
    return jump_prev, jump_next, rel_prev, rel_next


def jump_is_suspicious(jump_abs: float, jump_rel: float, args: Args) -> bool:
    return (
        isfinite(float(jump_abs))
        and isfinite(float(jump_rel))
        and (float(jump_abs) > args.jump_abs_threshold or float(jump_rel) > args.jump_rel_threshold)
    )


def point_is_suspicious(values: np.ndarray, sorted_zero: int, beta_index: int, args: Args) -> bool:
    jump_prev, jump_next, rel_prev, rel_next = row_jump_metrics(values, sorted_zero, beta_index)
    return jump_is_suspicious(jump_prev, rel_prev, args) or jump_is_suspicious(jump_next, rel_next, args)


def suspicious_mask(values: np.ndarray, args: Args) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    out = np.full(array.shape, False, dtype=bool)
    for sorted_zero in range(array.shape[0]):
        for beta_index in range(array.shape[1]):
            if isfinite(float(array[sorted_zero, beta_index])):
                out[sorted_zero, beta_index] = point_is_suspicious(array, sorted_zero, beta_index, args)
    return out


def recovery_mask(values: np.ndarray, beta_values: np.ndarray, args: Args) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    beta_array = np.asarray(beta_values, dtype=float)
    out = suspicious_mask(array, args)
    if not args.repair_spikes:
        return np.full(array.shape, False, dtype=bool)
    for sorted_zero in range(array.shape[0]):
        bad_abs_edges: list[int] = []
        for beta_index in range(array.shape[1] - 1):
            jump_abs, _jump_rel = jump_abs_rel(array[sorted_zero, beta_index], array[sorted_zero, beta_index + 1])
            if isfinite(jump_abs) and jump_abs > args.jump_abs_threshold:
                bad_abs_edges.append(beta_index)

        if not bad_abs_edges:
            continue

        for left_edge, right_edge in zip(bad_abs_edges, bad_abs_edges[1:]):
            interval_width = abs(float(beta_array[right_edge]) - float(beta_array[left_edge]))
            if interval_width > MAX_RECOVERY_INTERVAL_DEG + 1.0e-12:
                continue
            out[sorted_zero, left_edge + 1 : right_edge + 1] = True
    return out


def make_root_tables(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    timing_parts: dict[str, float] | None = None,
) -> RootTables:
    audit_start = time.perf_counter()
    raw: dict[tuple[int, str], np.ndarray] = {}
    clean: dict[tuple[int, str], np.ndarray] = {}
    status: dict[tuple[int, str], np.ndarray] = {}
    notes: dict[tuple[int, str], np.ndarray] = {}
    retry_attempted: dict[tuple[int, str], np.ndarray] = {}
    retry_fixed: dict[tuple[int, str], np.ndarray] = {}
    plotted_as_nan: dict[tuple[int, str], np.ndarray] = {}
    raw_suspicious: dict[tuple[int, str], np.ndarray] = {}
    recovery_suspicious: dict[tuple[int, str], np.ndarray] = {}
    plot_suspicious: dict[tuple[int, str], np.ndarray] = {}

    for case_index, _case in enumerate(cases):
        beta_values = np.asarray(beta_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            raw_array = np.full((args.n_roots, len(beta_values)), np.nan, dtype=float)
            status_array = np.full(raw_array.shape, "ok", dtype=object)
            notes_array = np.full(raw_array.shape, "", dtype=object)
            retry_array = np.full(raw_array.shape, False, dtype=bool)
            retry_fixed_array = np.full(raw_array.shape, False, dtype=bool)
            nan_array = np.full(raw_array.shape, False, dtype=bool)
            for beta_index, _beta_deg in enumerate(beta_values):
                result = root_map[(case_index, model, beta_index)]
                retry_array[:, beta_index] = bool(result.retry_attempted)
                for sorted_index in range(args.n_roots):
                    value = float(result.roots[sorted_index])
                    raw_array[sorted_index, beta_index] = value
                    status_array[sorted_index, beta_index] = "ok" if isfinite(value) else "missing_root_nan"
                    notes_array[sorted_index, beta_index] = "; ".join(result.notes)
            raw[key] = raw_array
            clean[key] = raw_array.copy()
            status[key] = status_array
            notes[key] = notes_array
            retry_attempted[key] = retry_array
            retry_fixed[key] = retry_fixed_array
            plotted_as_nan[key] = nan_array
            raw_suspicious[key] = suspicious_mask(raw_array, args)
            recovery_suspicious[key] = recovery_mask(raw_array, beta_values, args)

    audit_seconds = time.perf_counter() - audit_start
    repair_start = time.perf_counter()
    for case_index, case in enumerate(cases):
        beta_values = np.asarray(beta_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            values = clean[key]
            spike_retry_cache: dict[tuple[int, str, int], RootResult] = {}
            suspicious_beta_indices = sorted(
                {
                    beta_index
                    for sorted_zero in range(args.n_roots)
                    for beta_index in range(len(beta_values))
                    if bool(recovery_suspicious[key][sorted_zero, beta_index])
                }
            )
            for beta_index in suspicious_beta_indices:
                upper_hint = None
                if model == MODEL_TIMO:
                    eb_roots = root_map[(case_index, MODEL_EB, beta_index)].roots
                    finite_eb_roots = [float(root) for root in eb_roots if isfinite(float(root))]
                    upper_hint = max(finite_eb_roots) if finite_eb_roots else None
                retry_key = (case_index, model, beta_index)
                if retry_key not in spike_retry_cache:
                    if model == MODEL_EB and beta_index > 0:
                        spike_retry_cache[retry_key] = solve_eb_continuation(
                            case,
                            float(beta_values[beta_index]),
                            args.n_roots,
                            tuple(float(value) for value in values[:, beta_index - 1]),
                        )
                    elif (
                        model == MODEL_TIMO
                        and beta_index > 0
                        and not in_known_timo_spike_window(case, float(beta_values[beta_index]), args)
                    ):
                        spike_retry_cache[retry_key] = solve_timo_continuation(
                            case,
                            float(beta_values[beta_index]),
                            args.n_roots,
                            tuple(float(value) for value in values[:, beta_index - 1]),
                            upper_hint=upper_hint,
                        )
                    else:
                        spike_retry_cache[retry_key] = solve_model(
                            case,
                            float(beta_values[beta_index]),
                            args.n_roots,
                            model,
                            retry=True,
                            upper_hint=upper_hint,
                            sv_recovery=False,
                        )
                retry = spike_retry_cache[retry_key]
                for sorted_zero in range(args.n_roots):
                    if not bool(recovery_suspicious[key][sorted_zero, beta_index]):
                        continue
                    is_raw_suspicious = bool(raw_suspicious[key][sorted_zero, beta_index])
                    retry_attempted[key][sorted_zero, beta_index] = True
                    retry_value = float(retry.roots[sorted_zero])
                    current = float(values[sorted_zero, beta_index])
                    if isfinite(retry_value) and abs(retry_value - current) > 1.0e-8:
                        values[sorted_zero, beta_index] = retry_value
                        retry_fixed[key][sorted_zero, beta_index] = True
                        status[key][sorted_zero, beta_index] = (
                            "retry_fixed_spike" if is_raw_suspicious else "retry_fixed_spike_interval"
                        )
                        notes[key][sorted_zero, beta_index] = (
                            "suspicious jump fixed by dense retry"
                            if is_raw_suspicious
                            else "dense retry fixed point inside suspicious jump interval"
                        )
                    elif not is_raw_suspicious and notes[key][sorted_zero, beta_index] == "ok":
                        notes[key][sorted_zero, beta_index] = (
                            "dense retry attempted inside suspicious jump interval; initial kept"
                        )

            sv_retry_cache: dict[tuple[int, str, int, tuple[int, ...]], RootResult] = {}
            for _sv_pass in range(4):
                plot_mask_after_dense = suspicious_mask(values, args)
                sv_beta_indices = sorted(
                    {
                        beta_index
                        for sorted_zero in range(args.n_roots)
                        for beta_index in range(len(beta_values))
                        if bool(plot_mask_after_dense[sorted_zero, beta_index])
                        and bool(args.repair_spikes)
                        and bool(args.sv_recovery_only_on_spikes)
                    }
                )
                if not sv_beta_indices:
                    break
                changed_in_pass = False
                for beta_index in sv_beta_indices:
                    target_indices = tuple(
                        sorted_zero
                        for sorted_zero in range(args.n_roots)
                        if bool(plot_mask_after_dense[sorted_zero, beta_index])
                    )
                    upper_hint = None
                    if model == MODEL_TIMO:
                        eb_roots = root_map[(case_index, MODEL_EB, beta_index)].roots
                        finite_eb_roots = [float(root) for root in eb_roots if isfinite(float(root))]
                        upper_hint = max(finite_eb_roots) if finite_eb_roots else None
                    retry_key = (case_index, model, beta_index, target_indices)
                    if retry_key not in sv_retry_cache:
                        sv_retry_cache[retry_key] = solve_model(
                            case,
                            float(beta_values[beta_index]),
                            args.n_roots,
                            model,
                            retry=True,
                            upper_hint=upper_hint,
                            sv_recovery=True,
                            sv_target_indices=target_indices,
                        )
                    retry = sv_retry_cache[retry_key]
                    for sorted_zero in target_indices:
                        retry_attempted[key][sorted_zero, beta_index] = True
                        retry_value = float(retry.roots[sorted_zero])
                        current = float(values[sorted_zero, beta_index])
                        if isfinite(retry_value) and abs(retry_value - current) > 1.0e-8:
                            values[sorted_zero, beta_index] = retry_value
                            retry_fixed[key][sorted_zero, beta_index] = True
                            changed_in_pass = True
                            status[key][sorted_zero, beta_index] = (
                                "retry_fixed_spike_svd"
                                if bool(raw_suspicious[key][sorted_zero, beta_index])
                                else "retry_fixed_spike_interval_svd"
                            )
                            notes[key][sorted_zero, beta_index] = "suspicious jump fixed by local SVD recovery"
                if not changed_in_pass:
                    break

            plot_mask_before_nan = suspicious_mask(values, args)
            for sorted_zero in range(args.n_roots):
                for beta_index in range(len(beta_values)):
                    if not bool(plot_mask_before_nan[sorted_zero, beta_index]):
                        continue
                    if (
                        bool(args.repair_spikes)
                        and bool(raw_suspicious[key][sorted_zero, beta_index])
                        and not bool(retry_fixed[key][sorted_zero, beta_index])
                    ):
                        values[sorted_zero, beta_index] = float("nan")
                        plotted_as_nan[key][sorted_zero, beta_index] = True
                        status[key][sorted_zero, beta_index] = "unresolved_spike_nan"
                        notes[key][sorted_zero, beta_index] = "unresolved suspicious jump hidden from plot"
            plot_suspicious[key] = suspicious_mask(values, args)

    repair_seconds = time.perf_counter() - repair_start
    if timing_parts is not None:
        timing_parts["spike_audit_seconds"] = timing_parts.get("spike_audit_seconds", 0.0) + audit_seconds
        timing_parts["repair_seconds"] = timing_parts.get("repair_seconds", 0.0) + repair_seconds

    return RootTables(
        raw=raw,
        clean=clean,
        status=status,
        notes=notes,
        retry_attempted=retry_attempted,
        retry_fixed=retry_fixed,
        plotted_as_nan=plotted_as_nan,
        raw_suspicious=raw_suspicious,
        plot_suspicious=plot_suspicious,
    )


def warning_text(result: RootResult) -> str:
    return " | ".join(result.warnings)


def rel_diff_abs(timo: float, eb: float) -> float:
    if not (isfinite(timo) and isfinite(eb)) or abs(eb) <= 1.0e-30:
        return float("nan")
    return abs(float(timo) - float(eb)) / abs(float(eb))


def case_csv_path(output_dir: Path, case: CaseSpec) -> Path:
    return (
        output_dir
        / "lambda_beta_eb_vs_timo_"
        f"mu{_float_label(case.mu)}_eta{_float_label(case.eta)}_eps{_float_label(case.epsilon)}.csv"
    )


def case_png_path(output_dir: Path, case: CaseSpec) -> Path:
    return case_csv_path(output_dir, case).with_suffix(".png")


def build_case_rows(
    case_index: int,
    case: CaseSpec,
    beta_values: np.ndarray,
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    tables: RootTables,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    eb_key = (case_index, MODEL_EB)
    timo_key = (case_index, MODEL_TIMO)
    for beta_index, beta_deg in enumerate(beta_values):
        eb_result = root_map[(case_index, MODEL_EB, beta_index)]
        timo_result = root_map[(case_index, MODEL_TIMO, beta_index)]
        for sorted_zero in range(args.n_roots):
            lambda_eb = float(tables.clean[eb_key][sorted_zero, beta_index])
            lambda_timo = float(tables.clean[timo_key][sorted_zero, beta_index])
            suspicious_eb = bool(tables.raw_suspicious[eb_key][sorted_zero, beta_index])
            suspicious_timo = bool(tables.raw_suspicious[timo_key][sorted_zero, beta_index])
            retry_attempted_eb = bool(tables.retry_attempted[eb_key][sorted_zero, beta_index])
            retry_attempted_timo = bool(tables.retry_attempted[timo_key][sorted_zero, beta_index])
            retry_fixed_eb = bool(tables.retry_fixed[eb_key][sorted_zero, beta_index])
            retry_fixed_timo = bool(tables.retry_fixed[timo_key][sorted_zero, beta_index])
            plotted_as_nan_eb = bool(tables.plotted_as_nan[eb_key][sorted_zero, beta_index])
            plotted_as_nan_timo = bool(tables.plotted_as_nan[timo_key][sorted_zero, beta_index])
            notes_eb = str(tables.notes[eb_key][sorted_zero, beta_index])
            notes_timo = str(tables.notes[timo_key][sorted_zero, beta_index])
            rows.append(
                {
                    "beta_deg": float(beta_deg),
                    "sorted_index": sorted_zero + 1,
                    "Lambda_EB_raw": float(tables.raw[eb_key][sorted_zero, beta_index]),
                    "Lambda_Timoshenko_raw": float(tables.raw[timo_key][sorted_zero, beta_index]),
                    "Lambda_EB_plot": lambda_eb,
                    "Lambda_Timoshenko_plot": lambda_timo,
                    "rel_diff_abs_Timoshenko_vs_EB": rel_diff_abs(lambda_timo, lambda_eb),
                    "suspicious_EB": suspicious_eb,
                    "suspicious_Timoshenko": suspicious_timo,
                    "retry_attempted": retry_attempted_eb or retry_attempted_timo,
                    "retry_fixed": retry_fixed_eb or retry_fixed_timo,
                    "plotted_as_nan": plotted_as_nan_eb or plotted_as_nan_timo,
                    "notes": "; ".join(
                        item
                        for item in (f"EB: {notes_eb}" if notes_eb != "ok" else "", f"Timoshenko: {notes_timo}" if notes_timo != "ok" else "")
                        if item
                    )
                    or "ok",
                    "status_EB": str(tables.status[eb_key][sorted_zero, beta_index]),
                    "status_Timoshenko": str(tables.status[timo_key][sorted_zero, beta_index]),
                    "root_warning_EB": warning_text(eb_result),
                    "root_warning_Timoshenko": warning_text(timo_result),
                    "root_count_EB": int(eb_result.root_count_found),
                    "root_count_Timoshenko": int(timo_result.root_count_found),
                    "retry_attempted_EB": retry_attempted_eb,
                    "retry_attempted_Timoshenko": retry_attempted_timo,
                    "retry_fixed_EB": retry_fixed_eb,
                    "retry_fixed_Timoshenko": retry_fixed_timo,
                    "plotted_as_nan_EB": plotted_as_nan_eb,
                    "plotted_as_nan_Timoshenko": plotted_as_nan_timo,
                    "notes_EB": notes_eb,
                    "notes_Timoshenko": notes_timo,
                }
            )
    return rows


def build_summary_rows(case_rows_by_index: dict[int, list[dict[str, object]]], cases: Sequence[CaseSpec]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_index, case in enumerate(cases):
        case_rows = case_rows_by_index[case_index]
        n_roots = max(int(row["sorted_index"]) for row in case_rows)
        for sorted_index in range(1, n_roots + 1):
            selected = [row for row in case_rows if int(row["sorted_index"]) == sorted_index]
            finite_rel = [
                (float(row["beta_deg"]), float(row["rel_diff_abs_Timoshenko_vs_EB"]))
                for row in selected
                if isfinite(float(row["rel_diff_abs_Timoshenko_vs_EB"]))
            ]
            if finite_rel:
                beta_at_max, max_rel = max(finite_rel, key=lambda item: item[1])
                mean_rel = float(np.mean([item[1] for item in finite_rel]))
            else:
                beta_at_max = float("nan")
                max_rel = float("nan")
                mean_rel = float("nan")

            statuses = [str(row["status_EB"]) for row in selected] + [
                str(row["status_Timoshenko"]) for row in selected
            ]
            nan_count = sum(
                1
                for row in selected
                for key in ("Lambda_EB_plot", "Lambda_Timoshenko_plot")
                if not isfinite(float(row[key]))
            )
            rows.append(
                {
                    "mu": case.mu,
                    "eta": case.eta,
                    "epsilon": case.epsilon,
                    "sorted_index": sorted_index,
                    "max_rel_diff_over_beta": max_rel,
                    "mean_rel_diff_over_beta": mean_rel,
                    "beta_at_max_rel_diff": beta_at_max,
                    "raw_suspicious_point_count": sum(
                        1
                        for row in selected
                        for key in ("suspicious_EB", "suspicious_Timoshenko")
                        if str(row[key]) == "True" or row[key] is True
                    ),
                    "suspicious_point_count": sum(1 for status in statuses if status == "unresolved_spike_nan"),
                    "retry_fixed_count": sum(1 for status in statuses if status.startswith("retry_fixed")),
                    "nan_count_after_cleanup": nan_count,
                }
            )
    return rows


def plotted_series(
    tables: RootTables,
    *,
    case_index: int,
    model: str,
    sorted_index: int,
) -> np.ndarray:
    return tables.clean[(case_index, model)][int(sorted_index) - 1, :]


def plot_case(
    case_index: int,
    case: CaseSpec,
    beta_values: np.ndarray,
    args: Args,
    tables: RootTables,
) -> Path:
    fig, ax = plt.subplots(figsize=(10.6, 6.4), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(args.n_roots, 2)))
    for sorted_index in range(1, args.n_roots + 1):
        color = colors[sorted_index - 1]
        ax.plot(
            beta_values,
            plotted_series(tables, case_index=case_index, model=MODEL_TIMO, sorted_index=sorted_index),
            color=color,
            linestyle="-",
            linewidth=1.65,
            alpha=0.95,
        )
        ax.plot(
            beta_values,
            plotted_series(tables, case_index=case_index, model=MODEL_EB, sorted_index=sorted_index),
            color=color,
            linestyle="--",
            linewidth=1.3,
            alpha=0.85,
        )
    ax.set_xlabel("beta, degrees")
    ax.set_ylabel("Lambda")
    ax.set_title(
        "Sorted in-plane Lambda(beta): EB dashed vs Timoshenko solid\n"
        f"mu={case.mu:g}, eta={case.eta:g}, epsilon={case.epsilon:g}"
    )
    ax.grid(True, alpha=0.25)
    style_legend = ax.legend(
        handles=[
            Line2D([0], [0], color="black", linestyle="-", linewidth=1.8, label="Timoshenko"),
            Line2D([0], [0], color="black", linestyle="--", linewidth=1.6, label="Euler-Bernoulli"),
        ],
        loc="upper left",
        frameon=False,
    )
    ax.add_artist(style_legend)
    mode_handles = [
        Line2D([0], [0], color=colors[index - 1], linestyle="-", linewidth=1.8, label=f"sorted {index}")
        for index in range(1, args.n_roots + 1)
    ]
    ax.legend(handles=mode_handles, loc="upper right", ncols=2, fontsize=8, frameon=False)
    path = case_png_path(args.output_dir, case)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_overview(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    tables: RootTables,
) -> Path:
    epsilons = []
    for case in cases:
        if not any(abs(case.epsilon - existing) <= 1.0e-12 for existing in epsilons):
            epsilons.append(float(case.epsilon))
    mu_eta_cases = []
    for case in cases:
        pair = (float(case.mu), float(case.eta))
        if not any(abs(pair[0] - existing[0]) <= 1.0e-12 and abs(pair[1] - existing[1]) <= 1.0e-12 for existing in mu_eta_cases):
            mu_eta_cases.append(pair)
    fig, axes = plt.subplots(
        len(epsilons),
        len(mu_eta_cases),
        figsize=(4.7 * len(mu_eta_cases), 3.55 * len(epsilons)),
        squeeze=False,
        constrained_layout=True,
    )
    case_lookup = {(case.mu, case.eta, case.epsilon): index for index, case in enumerate(cases)}
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(args.n_roots, 2)))
    for row_index, epsilon in enumerate(epsilons):
        for col_index, (mu, eta) in enumerate(mu_eta_cases):
            ax = axes[row_index][col_index]
            case_index = case_lookup.get((float(mu), float(eta), float(epsilon)))
            if case_index is None:
                ax.set_axis_off()
                continue
            beta_values = np.asarray(beta_values_by_case[case_index], dtype=float)
            for sorted_index in range(1, args.n_roots + 1):
                color = colors[sorted_index - 1]
                ax.plot(
                    beta_values,
                    plotted_series(tables, case_index=case_index, model=MODEL_TIMO, sorted_index=sorted_index),
                    color=color,
                    linestyle="-",
                    linewidth=1.05,
                    alpha=0.95,
                )
                ax.plot(
                    beta_values,
                    plotted_series(tables, case_index=case_index, model=MODEL_EB, sorted_index=sorted_index),
                    color=color,
                    linestyle="--",
                    linewidth=0.95,
                    alpha=0.8,
                )
            ax.set_title(f"mu={float(mu):g}, eta={float(eta):g}, eps={float(epsilon):g}", fontsize=10)
            ax.grid(True, alpha=0.2)
            if row_index == len(epsilons) - 1:
                ax.set_xlabel("beta, degrees")
            if col_index == 0:
                ax.set_ylabel("Lambda")
    fig.suptitle("EB dashed vs Timoshenko solid, sorted in-plane Lambda(beta)", fontsize=12)
    path = args.output_dir / "eb_vs_timo_lambda_beta_overview.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def warning_counter(root_map: dict[tuple[int, str, int], RootResult]) -> Counter[str]:
    counter: Counter[str] = Counter()
    for result in root_map.values():
        for warning in result.warnings:
            counter[warning] += 1
    return counter


def missing_root_points(root_map: dict[tuple[int, str, int], RootResult], n_roots: int) -> int:
    return sum(max(0, int(n_roots) - int(result.root_count_found)) for result in root_map.values())


def status_counter(case_rows_by_index: dict[int, list[dict[str, object]]]) -> Counter[str]:
    counter: Counter[str] = Counter()
    for rows in case_rows_by_index.values():
        for row in rows:
            counter[str(row["status_EB"])] += 1
            counter[str(row["status_Timoshenko"])] += 1
    return counter


def build_spike_audit_rows(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    tables: RootTables,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_index, case in enumerate(cases):
        beta_values = np.asarray(beta_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            for sorted_zero in range(args.n_roots):
                for beta_index, beta_deg in enumerate(beta_values):
                    raw_prev, raw_next, raw_rel_prev, raw_rel_next = row_jump_metrics(
                        tables.raw[key],
                        sorted_zero,
                        beta_index,
                    )
                    plot_prev, plot_next, plot_rel_prev, plot_rel_next = row_jump_metrics(
                        tables.clean[key],
                        sorted_zero,
                        beta_index,
                    )
                    rows.append(
                        {
                            "case_index": case_index + 1,
                            "mu": case.mu,
                            "eta": case.eta,
                            "epsilon": case.epsilon,
                            "model": model,
                            "sorted_index": sorted_zero + 1,
                            "beta_deg": float(beta_deg),
                            "Lambda_raw": float(tables.raw[key][sorted_zero, beta_index]),
                            "Lambda_plot": float(tables.clean[key][sorted_zero, beta_index]),
                            "jump_prev_raw": raw_prev,
                            "jump_next_raw": raw_next,
                            "jump_rel_prev_raw": raw_rel_prev,
                            "jump_rel_next_raw": raw_rel_next,
                            "suspicious_raw": bool(tables.raw_suspicious[key][sorted_zero, beta_index]),
                            "retry_attempted": bool(tables.retry_attempted[key][sorted_zero, beta_index]),
                            "retry_fixed": bool(tables.retry_fixed[key][sorted_zero, beta_index]),
                            "plotted_as_nan": bool(tables.plotted_as_nan[key][sorted_zero, beta_index]),
                            "jump_prev_plot": plot_prev,
                            "jump_next_plot": plot_next,
                            "jump_rel_prev_plot": plot_rel_prev,
                            "jump_rel_next_plot": plot_rel_next,
                            "suspicious_plot": bool(tables.plot_suspicious[key][sorted_zero, beta_index]),
                            "notes": str(tables.notes[key][sorted_zero, beta_index]),
                        }
                    )
    return rows


def markdown_summary_table(summary_rows: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| mu | eta | epsilon | sorted | max rel diff | mean rel diff | beta at max | raw suspicious | unresolved suspicious | retry fixed | NaN after cleanup |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{float(row['mu']):.6g} | "
            f"{float(row['eta']):.6g} | "
            f"{float(row['epsilon']):.6g} | "
            f"{int(row['sorted_index'])} | "
            f"{float(row['max_rel_diff_over_beta']):.6g} | "
            f"{float(row['mean_rel_diff_over_beta']):.6g} | "
            f"{float(row['beta_at_max_rel_diff']):.6g} | "
            f"{int(row['raw_suspicious_point_count'])} | "
            f"{int(row['suspicious_point_count'])} | "
            f"{int(row['retry_fixed_count'])} | "
            f"{int(row['nan_count_after_cleanup'])} |"
        )
    return lines


def case_level_maxima(summary_rows: Sequence[dict[str, object]]) -> dict[tuple[float, float, float], float]:
    out: dict[tuple[float, float, float], float] = {}
    for row in summary_rows:
        key = (float(row["mu"]), float(row["eta"]), float(row["epsilon"]))
        value = float(row["max_rel_diff_over_beta"])
        if not isfinite(value):
            continue
        out[key] = max(out.get(key, 0.0), value)
    return out


def overall_observation_lines(summary_rows: Sequence[dict[str, object]]) -> list[str]:
    lines: list[str] = []
    by_epsilon: dict[float, list[float]] = {}
    for row in summary_rows:
        value = float(row["max_rel_diff_over_beta"])
        if isfinite(value):
            by_epsilon.setdefault(float(row["epsilon"]), []).append(value)
    for epsilon in DEFAULT_EPSILON_VALUES:
        values = by_epsilon.get(float(epsilon), [])
        if values:
            lines.append(
                f"- epsilon={float(epsilon):g}: largest case/index max relative difference is "
                f"{max(values):.6g}; mean over case/index maxima is {float(np.mean(values)):.6g}."
            )

    small_values = by_epsilon.get(float(DEFAULT_EPSILON_VALUES[0]), [])
    large_values = by_epsilon.get(float(DEFAULT_EPSILON_VALUES[1]), [])
    if small_values and large_values:
        relation = "larger" if max(large_values) > max(small_values) else "not larger"
        lines.append(
            f"- The epsilon={DEFAULT_EPSILON_VALUES[1]:g} differences are {relation} than the "
            f"epsilon={DEFAULT_EPSILON_VALUES[0]:g} differences in this diagnostic set."
        )

    maxima = case_level_maxima(summary_rows)
    for epsilon in DEFAULT_EPSILON_VALUES:
        pieces = []
        for mu, eta in DEFAULT_MU_ETA_CASES:
            value = maxima.get((float(mu), float(eta), float(epsilon)), float("nan"))
            pieces.append(f"(mu={float(mu):g}, eta={float(eta):g}) {value:.6g}")
        lines.append(f"- At epsilon={float(epsilon):g}, case-level max relative differences: " + "; ".join(pieces) + ".")
    lines.append(
        "- The mu/eta changes alter the magnitude and beta location of the EB/Timoshenko gap; "
        "this report does not promote those sorted-frequency differences to a modal-character claim."
    )
    return lines


def write_report(
    path: Path,
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    case_rows_by_index: dict[int, list[dict[str, object]]],
    summary_rows: Sequence[dict[str, object]],
    spike_audit_rows: Sequence[dict[str, object]],
    output_paths: Sequence[Path],
) -> None:
    warnings = warning_counter(root_map)
    statuses = status_counter(case_rows_by_index)
    nan_count = statuses.get("missing_root_nan", 0) + statuses.get("unresolved_spike_nan", 0)
    raw_suspicious = sum(1 for row in spike_audit_rows if bool(row["suspicious_raw"]))
    final_suspicious = sum(1 for row in spike_audit_rows if bool(row["suspicious_plot"]))
    retry_fixed = sum(count for status, count in statuses.items() if status.startswith("retry_fixed"))
    refined_windows = [
        (case_index, case, window)
        for case_index, case in enumerate(cases)
        for window in refinement_windows_for_case(case, args)
    ]
    raw_suspicious_cases = sorted(
        {
            (float(row["mu"]), float(row["eta"]), float(row["epsilon"]))
            for row in spike_audit_rows
            if bool(row["suspicious_raw"])
        }
    )

    lines: list[str] = [
        "# EB vs Timoshenko Lambda(beta) Diagnostic Cases",
        "",
        "Diagnostic only. These plots use sorted in-plane frequencies at each beta, not descendant branch tracking.",
        "The workflow reuses the existing Euler-Bernoulli thickness-mismatch root helper and the existing variable-length Timoshenko helper.",
        "No analytic formulas, determinants, root solvers, Timoshenko shear coefficient, FEM workflows, article files, article figures, or baseline results were changed.",
        "",
        "## Cases",
        "",
        "| # | mu | eta | epsilon |",
        "| ---: | ---: | ---: | ---: |",
    ]
    for index, case in enumerate(cases, start=1):
        lines.append(f"| {index} | {case.mu:g} | {case.eta:g} | {case.epsilon:g} |")

    lines.extend(
        [
            "",
            "## Beta Grid",
            "",
            f"- base beta range: {args.beta_min:g} to {args.beta_max:g} degrees",
            f"- base beta step: {args.beta_step:g} degrees",
            f"- refined beta step: {args.refined_beta_step:g} degrees",
            f"- sorted frequencies per theory: {args.n_roots}",
            f"- spike thresholds: jump_abs > {args.jump_abs_threshold:g} or jump_rel > {args.jump_rel_threshold:g}",
            "",
            "Per-case beta point counts:",
            "",
            "| case | mu | eta | epsilon | beta points |",
            "| ---: | ---: | ---: | ---: | ---: |",
        ]
    )
    for case_index, case in enumerate(cases):
        lines.append(
            f"| {case_index + 1} | {case.mu:g} | {case.eta:g} | {case.epsilon:g} | "
            f"{len(beta_values_by_case[case_index])} |"
        )
    lines.extend(
        [
            "",
            "Forced refined windows:",
            "",
        ]
    )
    if refined_windows:
        for _case_index, case, window in refined_windows:
            lines.append(
                f"- mu={case.mu:g}, eta={case.eta:g}, epsilon={case.epsilon:g}: "
                f"beta {window.beta_min:g}..{window.beta_max:g} deg at step {args.refined_beta_step:g} deg "
                f"({window.reason})."
            )
    else:
        lines.append("- none")
    lines.extend(
        [
            "",
            "## Line Style",
            "",
            "- Timoshenko curves are solid.",
            "- Euler-Bernoulli curves are dashed.",
            "- Within each plot, the same sorted index uses the same color for both theories.",
            "",
            "## Root Warnings And Cleanup",
            "",
            f"- root warning message classes: {len(warnings)}",
            f"- root warning occurrences: {sum(warnings.values())}",
            f"- missing root slots after initial/retry solving: {missing_root_points(root_map, args.n_roots)}",
            f"- raw suspicious points detected: {raw_suspicious}",
            f"- retry-fixed suspicious points: {retry_fixed}",
            f"- unresolved NaN statuses after cleanup: {nan_count}",
            f"- final plotted suspicious points: {final_suspicious}",
            "",
        ]
    )
    if raw_suspicious_cases:
        lines.append("Cases with raw spike artifacts before recovery:")
        for mu, eta, epsilon in raw_suspicious_cases:
            lines.append(f"- mu={mu:g}, eta={eta:g}, epsilon={epsilon:g}")
        lines.append("")
    else:
        lines.append("No raw spike artifacts were detected by the configured jump thresholds.")
        lines.append("")
    if warnings:
        lines.append("Most frequent root warnings:")
        for warning, count in warnings.most_common(12):
            lines.append(f"- {count}x: {warning}")
        lines.append("")
    else:
        lines.append("No root-search or Timoshenko basis warnings were recorded.")
        lines.append("")

    lines.extend(["## Numeric Summary", ""])
    lines.extend(markdown_summary_table(summary_rows))
    lines.append("")
    lines.extend(["## General Diagnostic Observations", ""])
    lines.extend(overall_observation_lines(summary_rows))
    lines.append("")
    lines.append(
        "Visual artifact status: "
        + ("no artificial jump artifacts remain by the configured plotted-curve check." if final_suspicious == 0 else "inspect remaining suspicious plotted points in the spike audit CSV.")
    )
    lines.append("")
    lines.extend(["## Outputs", ""])
    for output_path in output_paths:
        lines.append(f"- `{_rel(output_path)}`")
    lines.append("")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_outputs(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    tables: RootTables,
) -> tuple[list[Path], list[dict[str, object]]]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_paths: list[Path] = []
    case_rows_by_index: dict[int, list[dict[str, object]]] = {}
    for case_index, case in enumerate(cases):
        beta_values = np.asarray(beta_values_by_case[case_index], dtype=float)
        rows = build_case_rows(case_index, case, beta_values, args, root_map, tables)
        case_rows_by_index[case_index] = rows
        csv_path = case_csv_path(args.output_dir, case)
        _write_csv(csv_path, rows, CASE_CSV_FIELDS)
        output_paths.append(csv_path)
        output_paths.append(plot_case(case_index, case, beta_values, args, tables))

    summary_rows = build_summary_rows(case_rows_by_index, cases)
    summary_path = args.output_dir / "eb_vs_timo_lambda_beta_case_summary.csv"
    _write_csv(summary_path, summary_rows, SUMMARY_FIELDS)
    output_paths.append(summary_path)

    spike_audit_rows = build_spike_audit_rows(cases, beta_values_by_case, args, tables)
    spike_audit_path = args.output_dir / "eb_vs_timo_lambda_beta_spike_audit.csv"
    _write_csv(spike_audit_path, spike_audit_rows, SPIKE_AUDIT_FIELDS)
    output_paths.append(spike_audit_path)

    overview_path = plot_overview(cases, beta_values_by_case, args, tables)
    output_paths.append(overview_path)

    report_path = args.output_dir / "eb_vs_timo_lambda_beta_cases_report.md"
    write_report(
        report_path,
        cases,
        beta_values_by_case,
        args,
        root_map,
        case_rows_by_index,
        summary_rows,
        spike_audit_rows,
        [*output_paths, report_path],
    )
    output_paths.append(report_path)
    return output_paths, summary_rows


def print_run_summary(output_paths: Sequence[Path], summary_rows: Sequence[dict[str, object]]) -> None:
    print("saved outputs:")
    for path in output_paths:
        print(f"- {_rel(path)}")
    finite_maxima = [float(row["max_rel_diff_over_beta"]) for row in summary_rows if isfinite(float(row["max_rel_diff_over_beta"]))]
    if finite_maxima:
        print(f"largest max relative difference over all rows: {max(finite_maxima):.8g}")
    print("sorted frequencies only; no descendant tracking; no FEM/Gmsh/CalculiX")


def note_contains(text: object, needle: str) -> bool:
    return needle in str(text)


def count_case_model_warnings(
    root_map: dict[tuple[int, str, int], RootResult],
    case_index: int,
    model: str,
    n_beta: int,
) -> int:
    return sum(len(root_map[(case_index, model, beta_index)].warnings) for beta_index in range(n_beta))


def count_case_model_fallbacks(
    root_map: dict[tuple[int, str, int], RootResult],
    case_index: int,
    model: str,
    n_beta: int,
) -> int:
    return sum(
        1
        for beta_index in range(n_beta)
        if any("continuation fallback" in note for note in root_map[(case_index, model, beta_index)].notes)
    )


def count_case_model_table_flags(
    tables: RootTables,
    case_index: int,
    model: str,
    name: str,
) -> int:
    table = getattr(tables, name)[(case_index, model)]
    return int(np.count_nonzero(np.asarray(table, dtype=bool)))


def count_case_model_svd_recovery(tables: RootTables, case_index: int, model: str) -> int:
    notes = np.asarray(tables.notes[(case_index, model)], dtype=object)
    return sum(1 for value in notes.ravel() if note_contains(value, "SVD recovery added"))


def build_timing_rows(
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    tables: RootTables,
    timing: TimingStats,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_index, case in enumerate(cases):
        n_beta = len(beta_values_by_case[case_index])
        for model in MODELS:
            rows.append(
                {
                    "run_id": timing.run_id,
                    "model": model,
                    "mu": case.mu,
                    "eta": case.eta,
                    "epsilon": case.epsilon,
                    "n_beta_points": n_beta,
                    "n_roots": args.n_roots,
                    "root_mode": args.timo_root_mode if model == MODEL_TIMO else "global",
                    "cache_hit": timing.cache_hit,
                    "ordinary_compute_seconds": timing.ordinary_compute_seconds,
                    "spike_audit_seconds": timing.spike_audit_seconds,
                    "repair_seconds": timing.repair_seconds,
                    "plotting_seconds": timing.plotting_seconds,
                    "total_seconds": timing.total_seconds,
                    "warnings": count_case_model_warnings(root_map, case_index, model, n_beta),
                    "fallback_count": count_case_model_fallbacks(root_map, case_index, model, n_beta),
                    "repair_count": count_case_model_table_flags(tables, case_index, model, "retry_fixed"),
                    "sv_recovery_calls": count_case_model_svd_recovery(tables, case_index, model),
                    "notes": "; ".join(timing.notes) or "ok",
                }
            )
    return rows


def write_timing_report(
    path: Path,
    cases: Sequence[CaseSpec],
    beta_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, int], RootResult],
    tables: RootTables,
    timing: TimingStats,
) -> list[dict[str, object]]:
    rows = build_timing_rows(cases, beta_values_by_case, args, root_map, tables, timing)
    _write_csv(path, rows, TIMING_FIELDS)
    return rows


def append_timing_section(report_path: Path, timing_path: Path, timing_rows: Sequence[dict[str, object]]) -> None:
    total_runtime = max((float(row["total_seconds"]) for row in timing_rows), default=0.0)
    cache_hits = sum(1 for row in timing_rows if bool(row["cache_hit"]))
    repair_count = sum(int(row["repair_count"]) for row in timing_rows)
    sv_calls = sum(int(row["sv_recovery_calls"]) for row in timing_rows)
    fallback_count = sum(int(row["fallback_count"]) for row in timing_rows)
    cache_hit_text = f"{cache_hits}/{len(timing_rows)} timing rows"
    lines = [
        "",
        "## Timing And Cache",
        "",
        f"- timing report: `{_rel(timing_path)}`",
        f"- total runtime: {total_runtime:.6g} s",
        f"- cache hits: {cache_hit_text}",
        f"- spike repair fixed values: {repair_count}",
        f"- singular-value recovery additions recorded in notes: {sv_calls}",
        f"- Timoshenko continuation global fallbacks: {fallback_count}",
        "- Compared with the previous slow behavior, cache and plot-only modes avoid root recomputation entirely when settings match.",
    ]
    with report_path.open("a", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main(argv: list[str] | None = None) -> None:
    run_start = time.perf_counter()
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    args = parse_args(argv)
    cases = default_cases(args.mu_eta_cases, args.epsilon_values)
    beta_values_by_case = beta_grids_by_case(cases, args)
    cache_path = cache_file_path(cases, beta_values_by_case, args)
    cache_hit = False
    timing_parts: dict[str, float] = {"spike_audit_seconds": 0.0, "repair_seconds": 0.0}
    notes: list[str] = []

    cached = None
    if args.reuse_cache and not args.force_recompute:
        cached = load_cache(cache_path, cases, beta_values_by_case, args)
        if cached is not None:
            cache_hit = True
            notes.append(f"loaded cache {_rel(cache_path)}")

    if args.plot_only and cached is None:
        raise FileNotFoundError(
            "plot-only mode requires a matching cache. Run without --plot-only or with --force-recompute first."
        )

    if cached is not None:
        root_map, tables = cached
        ordinary_compute_seconds = 0.0
    else:
        compute_start = time.perf_counter()
        root_map = compute_root_map(cases, beta_values_by_case, args)
        ordinary_compute_seconds = time.perf_counter() - compute_start
        tables = make_root_tables(cases, beta_values_by_case, args, root_map, timing_parts=timing_parts)
        save_cache(cache_path, cases, beta_values_by_case, args, root_map, tables)
        notes.append(f"saved cache {_rel(cache_path)}")

    plotting_start = time.perf_counter()
    output_paths, summary_rows = write_outputs(cases, beta_values_by_case, args, root_map, tables)
    plotting_seconds = time.perf_counter() - plotting_start
    total_seconds = time.perf_counter() - run_start
    timing = TimingStats(
        run_id=run_id,
        cache_hit=cache_hit,
        ordinary_compute_seconds=ordinary_compute_seconds,
        spike_audit_seconds=float(timing_parts.get("spike_audit_seconds", 0.0)),
        repair_seconds=float(timing_parts.get("repair_seconds", 0.0)),
        plotting_seconds=plotting_seconds,
        total_seconds=total_seconds,
        fallback_count=sum(
            count_case_model_fallbacks(root_map, case_index, MODEL_TIMO, len(beta_values_by_case[case_index]))
            for case_index in range(len(cases))
        ),
        repair_count=sum(
            count_case_model_table_flags(tables, case_index, model, "retry_fixed")
            for case_index in range(len(cases))
            for model in MODELS
        ),
        sv_recovery_calls=sum(
            count_case_model_svd_recovery(tables, case_index, model)
            for case_index in range(len(cases))
            for model in MODELS
        ),
        notes=tuple(notes) or ("ok",),
    )
    timing_path = args.output_dir / TIMING_REPORT_NAME
    timing_rows = write_timing_report(timing_path, cases, beta_values_by_case, args, root_map, tables, timing)
    append_timing_section(args.output_dir / "eb_vs_timo_lambda_beta_cases_report.md", timing_path, timing_rows)
    output_paths.append(timing_path)
    print_run_summary(output_paths, summary_rows)
    print(
        "timing: "
        f"total={timing.total_seconds:.3f}s, compute={timing.ordinary_compute_seconds:.3f}s, "
        f"audit={timing.spike_audit_seconds:.3f}s, repair={timing.repair_seconds:.3f}s, "
        f"plot={timing.plotting_seconds:.3f}s, cache_hit={timing.cache_hit}, "
        f"fallbacks={timing.fallback_count}, repairs={timing.repair_count}, "
        f"sv_recovery={timing.sv_recovery_calls}"
    )


if __name__ == "__main__":
    main()
