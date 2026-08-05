from __future__ import annotations

import argparse
import csv
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
from math import isfinite
import os
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
    find_first_n_roots_eta,
    find_roots_scan_bisect_eta,
    thickness_mismatch_factors,
)
from scripts.lib import analytic_branch_tracking as TRACK  # noqa: E402
from scripts.lib import general_spectrum_completeness as COMPLETE  # noqa: E402
from scripts.lib import straight_rod_factorized_spectrum as STRAIGHT  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import (  # noqa: E402
    analytic_shape_vector_eta,
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

PILOT_OUTPUT_DIR = (
    Path("results") / "article_epsilon_upper_envelope" / "beta_branch_pilot"
)
PILOT_CASES = (
    ("P1", 0.010, 0.0, 0.0),
    ("P2", 0.010, 0.9, -0.25),
    ("P3", 0.024798906738281248, 0.5, -0.1),
    ("P4", 0.050, 0.0, 0.0),
)
PILOT_N_TRACK = 12
PILOT_N_SOLVE = 18
PILOT_SHAPE_SAMPLES = 201
PILOT_MAC_THRESHOLD = 0.8
PILOT_MAC_MARGIN_THRESHOLD = 0.05
PILOT_MAX_SORTED_JUMP = 1
PILOT_MU_STEP = 0.025
PILOT_CSV_FIELDS = (
    "case_id",
    "epsilon_0",
    "mu",
    "eta",
    "beta_deg",
    "branch_index",
    "lambda_eb",
    "lambda_timo",
    "tracking_status",
)

SORTED_PILOT_OUTPUT_DIR = (
    Path("results") / "article_epsilon_upper_envelope" / "beta_sorted_spectrum_pilot"
)
SORTED_PILOT_CASES = PILOT_CASES
SORTED_PILOT_N_STORE = 12
SORTED_PILOT_N_PLOT = 10
SORTED_PILOT_CACHE_VERSION = "beta_sorted_spectrum_pointwise_primary_v1"
SORTED_PILOT_NOTICEABLE_JUMP_REL = 0.03
SORTED_PILOT_NOTICEABLE_JUMP_ABS = 0.5
SORTED_PILOT_SINGLE_JUMP_REL = 0.08
SORTED_PILOT_SIMULTANEOUS_COUNT = 3
SORTED_PILOT_SUGGESTED_BETA_STEP = 0.25
SORTED_PILOT_CSV_FIELDS = (
    "case_id",
    "epsilon_0",
    "mu",
    "eta",
    "beta_deg",
    "sorted_mode_index",
    "lambda_eb",
    "lambda_timo",
    "eb_root_count",
    "timo_root_count",
    "eb_candidate_root_count",
    "timo_candidate_root_count",
    "eb_point_status",
    "timo_point_status",
    "eb_root_source",
    "timo_root_source",
    "eb_lambda_min",
    "eb_lambda_max",
    "timo_lambda_min",
    "timo_lambda_max",
    "scan_step",
    "eb_matrix_evaluations",
    "timo_matrix_evaluations",
)
SORTED_PILOT_STEP_FIELDS = (
    "case_id",
    "theory",
    "sorted_mode_index",
    "beta_left",
    "beta_right",
    "lambda_left",
    "lambda_right",
    "absolute_jump",
    "relative_jump",
    "root_count_left",
    "root_count_right",
    "simultaneous_shift_count",
    "suspect_interval",
    "suspect_reason",
)
SORTED_PILOT_SUSPECT_FIELDS = (
    "case_id",
    "theory",
    "beta_left",
    "beta_right",
    "affected_sorted_modes",
    "reason",
    "max_relative_jump",
    "root_count_change",
    "suggested_local_beta_step",
    "suggested_lambda_interval_left",
    "suggested_lambda_interval_right",
)

REFINED_PILOT_OUTPUT_DIR = (
    Path("results") / "article_epsilon_upper_envelope" / "beta_sorted_spectrum_refined_pilot"
)
REFINED_PILOT_SOURCE_DIR = SORTED_PILOT_OUTPUT_DIR
REFINED_PILOT_CACHE_VERSION = "beta_sorted_spectrum_local_refinement_v2"
REFINED_PILOT_ALGORITHM_VERSION = "dense_local_two_phase_primary_candidates_v2"
REFINED_PILOT_CSV_FIELDS = (
    "case_id",
    "epsilon_0",
    "mu",
    "eta",
    "beta_deg",
    "beta_source",
    "sorted_mode_index",
    "lambda_eb",
    "lambda_timo",
    "eb_status",
    "timo_status",
    "eb_root_count",
    "timo_root_count",
    "eb_inventory_source",
    "timo_inventory_source",
    "eb_multiplicity",
    "timo_multiplicity",
    "eb_multiplicity_source",
    "timo_multiplicity_source",
    "eb_block_family",
    "timo_block_family",
    "eb_nullity",
    "timo_nullity",
    "eb_repeated_root_slot",
    "timo_repeated_root_slot",
    "local_region_id",
)
LOCAL_CANDIDATE_FIELDS = (
    "region_id",
    "case_id",
    "theory",
    "beta_deg",
    "lambda_candidate",
    "candidate_source",
    "sign_change",
    "sigma_min",
    "sigma_ratio",
    "residual",
    "bracket_left",
    "bracket_right",
    "multiplicity",
    "accepted",
    "rejection_reason",
    "block_family",
    "nullity",
    "multiplicity_source",
)
BEFORE_AFTER_FIELDS = (
    "case_id",
    "theory",
    "beta_deg",
    "sorted_mode_index",
    "lambda_before",
    "lambda_after",
    "absolute_difference",
    "status_before",
    "status_after",
    "root_inventory_shift",
)
LOCAL_REFINEMENT_SUMMARY_FIELDS = (
    "region_id",
    "beta_points",
    "resolved_points",
    "unresolved_points",
    "recovered_root_count",
    "multiplicity_two_count",
    "max_before_after_difference",
    "remaining_spike",
    "status",
)

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
    label: str = ""


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
    beta_branch_pilot: bool
    beta_sorted_spectrum_pilot: bool
    beta_sorted_spectrum_refined_pilot: bool


@dataclass(frozen=True)
class LocalRefinementRegion:
    region_id: str
    case_id: str
    model: str
    beta_min: float
    beta_max: float
    beta_step: float
    lambda_min: float
    lambda_max: float
    lambda_step: float
    expected_min_roots: int
    purpose: str


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
    parser.add_argument(
        "--beta-branch-pilot",
        action="store_true",
        help=(
            "run the fixed four-case, 12-descendant EB/Timoshenko beta pilot; "
            "ambiguous tracking is saved as NaN and no strict/refinement workflow is used"
        ),
    )
    parser.add_argument(
        "--beta-sorted-spectrum-pilot",
        action="store_true",
        help=(
            "run the fixed four-case independent pointwise sorted-spectrum pilot; "
            "no branch tracking, strict verification, or beta refinement is used"
        ),
    )
    parser.add_argument(
        "--beta-sorted-spectrum-refined-pilot",
        action="store_true",
        help=(
            "locally refine only the fixed R1-R7 windows of the independent "
            "sorted-spectrum pilot; no tracking or strict workflow is used"
        ),
    )
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
        beta_branch_pilot=bool(ns.beta_branch_pilot),
        beta_sorted_spectrum_pilot=bool(ns.beta_sorted_spectrum_pilot),
        beta_sorted_spectrum_refined_pilot=bool(ns.beta_sorted_spectrum_refined_pilot),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    pilot_count = sum(
        (
            args.beta_branch_pilot,
            args.beta_sorted_spectrum_pilot,
            args.beta_sorted_spectrum_refined_pilot,
        )
    )
    if pilot_count > 1:
        raise ValueError("beta pilot presets are mutually exclusive.")
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


def _pilot_normalized(vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(values))
    if not isfinite(norm) or norm <= 1.0e-28:
        raise ValueError("non-finite or zero tracking vector")
    return values / norm


def _pilot_components(vector: np.ndarray) -> dict[str, np.ndarray]:
    values = _pilot_normalized(vector)
    split = len(values) // 2
    return {
        "w_left": values[:split],
        "w_right": values[split:],
    }


def _pilot_roots(
    case: CaseSpec,
    beta_deg: float,
    model: str,
    previous_roots: Sequence[float] | None = None,
) -> tuple[float, ...]:
    if model == MODEL_EB:
        roots = find_first_n_roots_eta(
            beta=float(np.deg2rad(beta_deg)),
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            n_roots=PILOT_N_SOLVE,
            Lmin=ROOT_SCAN_START,
            Lmax0=45.0,
            scan_step=0.01,
            grow_factor=1.35,
            max_tries=8,
        )
    elif model == MODEL_TIMO:
        result = solve_timo_continuation(
            case,
            beta_deg,
            PILOT_N_SOLVE,
            previous_roots,
            upper_hint=(max(previous_roots) if previous_roots else None),
        )
        roots = result.roots
    else:  # pragma: no cover - protected by the fixed model tuple
        raise ValueError(f"unknown pilot model: {model}")
    return tuple(sorted(float(value) for value in roots if isfinite(float(value))))


def _pilot_state(
    case: CaseSpec,
    beta_deg: float,
    model: str,
    sorted_index: int,
    Lambda: float,
) -> TRACK.AnalyticModeState:
    if model == MODEL_EB:
        matrix = assemble_clamped_coupled_matrix_eta(
            Lambda, float(np.deg2rad(beta_deg)), case.mu, case.epsilon, case.eta
        )
        _u, singular, vh = np.linalg.svd(np.asarray(matrix, dtype=float))
        coeff = np.asarray(vh[-1], dtype=float)
        vector = analytic_shape_vector_eta(
            Lambda,
            beta_rad=float(np.deg2rad(beta_deg)),
            mu=case.mu,
            epsilon=case.epsilon,
            eta=case.eta,
            s_norm=np.linspace(0.0, 1.0, PILOT_SHAPE_SAMPLES),
        )
    else:
        mode = TIMO.timo_mode_coefficients(
            Lambda, beta_deg, case.mu, case.epsilon, case.eta
        )
        fields = TIMO.timo_mode_fields(
            Lambda,
            beta_deg,
            case.mu,
            case.epsilon,
            case.eta,
            coeff=mode.coeff,
            n_points=PILOT_SHAPE_SAMPLES,
        )
        rod1 = fields["rod1"]
        rod2 = fields["rod2"]
        vector = np.concatenate(
            [
                np.asarray(rod1["u"], dtype=float),  # type: ignore[index]
                np.asarray(rod1["w"], dtype=float),  # type: ignore[index]
                np.asarray(rod2["u"], dtype=float)[::-1],  # type: ignore[index]
                np.asarray(rod2["w"], dtype=float)[::-1],  # type: ignore[index]
            ]
        )
        coeff = np.asarray(mode.coeff, dtype=float)
        singular = np.array(
            [1.0, mode.smallest_singular_value], dtype=float
        )
    normalized = _pilot_normalized(vector)
    components = _pilot_components(normalized)
    smallest = float(singular[-1]) if len(singular) else float("nan")
    ratio = (
        float(singular[-1] / singular[-2])
        if len(singular) >= 2 and abs(float(singular[-2])) > 1.0e-28
        else float("nan")
    )
    return TRACK.AnalyticModeState(
        epsilon=case.epsilon,
        beta=beta_deg,
        mu=case.mu,
        current_sorted_index=sorted_index,
        Lambda=Lambda,
        coeff=coeff,
        components=components,
        shape_vector=normalized,
        smallest_singular_value=smallest,
        singular_value_ratio=ratio,
    )


def _pilot_states(
    case: CaseSpec,
    beta_deg: float,
    model: str,
    previous_roots: Sequence[float] | None = None,
) -> tuple[list[TRACK.AnalyticModeState], tuple[float, ...]]:
    roots = _pilot_roots(case, beta_deg, model, previous_roots)
    states: list[TRACK.AnalyticModeState] = []
    for sorted_index, root in enumerate(roots, start=1):
        try:
            states.append(
                _pilot_state(case, beta_deg, model, sorted_index, float(root))
            )
        except (FloatingPointError, ValueError, OverflowError, np.linalg.LinAlgError):
            continue
    return states, roots


def _pilot_seed_points(
    case: CaseSpec, model: str
) -> tuple[list[TRACK.BranchPoint], tuple[float, ...]]:
    states, roots = _pilot_states(case, 0.0, model)
    if len(states) < PILOT_N_TRACK:
        raise RuntimeError(
            f"{case.label} {model}: only {len(states)} usable seed states"
        )
    return [
        TRACK.make_branch_point(
            state=states[index],
            branch_id=TRACK.branch_id_from_base_sorted_index(
                index + 1, prefix="beta_desc"
            ),
            base_sorted_index=index + 1,
            mac_to_previous=1.0,
            relative_lambda_jump=0.0,
            sorted_lambdas=roots,
            step_type="base",
            step_index=0,
            mac_warning_threshold=PILOT_MAC_THRESHOLD,
        )
        for index in range(PILOT_N_TRACK)
    ], roots


def _pilot_transition(
    previous: Sequence[TRACK.BranchPoint],
    case: CaseSpec,
    beta_deg: float,
    model: str,
    *,
    step_type: str,
    step_index: int,
    previous_roots: Sequence[float] | None,
) -> tuple[list[TRACK.BranchPoint], list[str], list[float], tuple[float, ...]]:
    try:
        states, roots = _pilot_states(case, beta_deg, model, previous_roots)
        if len(states) < len(previous):
            raise RuntimeError(f"only {len(states)} usable candidate states")
        diagnostics = TRACK.assignment_diagnostics(
            previous_points=previous,
            current_states=states,
            sorted_lambdas=roots,
            freq_weight=TRACK.DEFAULT_FREQ_WEIGHT,
            shape_metric="transverse",
        )
        proposed = TRACK.build_next_points_from_assignment(
            source_points=previous,
            current_states=states,
            sorted_lambdas=roots,
            diagnostics=diagnostics,
            step_type=step_type,
            step_index=step_index,
            mac_warning_threshold=PILOT_MAC_THRESHOLD,
            shape_metric="transverse",
        )
    except (FloatingPointError, ValueError, OverflowError, RuntimeError, np.linalg.LinAlgError) as exc:
        return (
            list(previous),
            [f"gap_{model}:solve_or_assignment:{type(exc).__name__}"] * len(previous),
            [float("nan")] * len(previous),
            tuple(float(value) for value in (previous_roots or ())),
        )

    accepted: list[TRACK.BranchPoint] = []
    statuses: list[str] = []
    values: list[float] = []
    for row, candidate in enumerate(proposed):
        selected_col = int(diagnostics.assignment[row])
        selected_mac = float(diagnostics.mac[row, selected_col])
        other_macs = np.delete(np.asarray(diagnostics.mac[row], dtype=float), selected_col)
        second_best = float(np.nanmax(other_macs)) if other_macs.size else float("nan")
        margin = selected_mac - second_best if isfinite(second_best) else float("inf")
        sorted_jump = abs(
            int(candidate.current_sorted_index) - int(previous[row].current_sorted_index)
        )
        reasons: list[str] = []
        if not isfinite(selected_mac) or selected_mac < PILOT_MAC_THRESHOLD:
            reasons.append("low_continuity")
        if isfinite(margin) and margin < PILOT_MAC_MARGIN_THRESHOLD:
            reasons.append("ambiguous_margin")
        if sorted_jump > PILOT_MAX_SORTED_JUMP:
            reasons.append("large_sorted_jump")
        if reasons:
            retained_col = next(
                (
                    index
                    for index, state in enumerate(diagnostics.current_states)
                    if int(state.current_sorted_index)
                    == int(previous[row].current_sorted_index)
                ),
                None,
            )
            if retained_col is None:
                retained_point = previous[row]
            else:
                retained_state = TRACK.align_state_to_previous(
                    diagnostics.current_states[retained_col],
                    previous[row],
                    shape_metric="transverse",
                )
                relative_jump = abs(
                    float(retained_state.Lambda) - float(previous[row].Lambda)
                ) / max(abs(float(previous[row].Lambda)), 1.0e-28)
                retained_point = TRACK.make_branch_point(
                    state=retained_state,
                    branch_id=previous[row].branch_id,
                    base_sorted_index=previous[row].base_sorted_index,
                    mac_to_previous=float(diagnostics.mac[row, retained_col]),
                    relative_lambda_jump=relative_jump,
                    sorted_lambdas=roots,
                    step_type=step_type,
                    step_index=step_index,
                    mac_warning_threshold=PILOT_MAC_THRESHOLD,
                )
            accepted.append(retained_point)
            statuses.append(f"gap_{model}:" + ";".join(reasons))
            values.append(float("nan"))
        else:
            accepted.append(candidate)
            statuses.append("ok")
            values.append(float(candidate.Lambda))
    return accepted, statuses, values, roots


def _pilot_track_model(case: CaseSpec, model: str) -> dict[str, object]:
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = "1"
    seed_case = CaseSpec(mu=0.0, eta=case.eta, epsilon=case.epsilon, label=case.label)
    previous, root_seeds = _pilot_seed_points(seed_case, model)
    beta0_status = ["ok"] * PILOT_N_TRACK
    beta0_values = [float(point.Lambda) for point in previous]
    if abs(case.mu) > 1.0e-14:
        mu_count = max(2, int(np.ceil(abs(case.mu) / PILOT_MU_STEP)) + 1)
        for mu_index, mu_value in enumerate(np.linspace(0.0, case.mu, mu_count)[1:], start=1):
            current_case = CaseSpec(
                mu=float(mu_value), eta=case.eta, epsilon=case.epsilon, label=case.label
            )
            previous, beta0_status, beta0_values, root_seeds = _pilot_transition(
                previous,
                current_case,
                0.0,
                model,
                step_type="mu_seed",
                step_index=mu_index,
                previous_roots=root_seeds,
            )
    values = np.full((PILOT_N_TRACK, 91), np.nan, dtype=float)
    statuses = np.empty((PILOT_N_TRACK, 91), dtype=object)
    values[:, 0] = np.asarray(beta0_values, dtype=float)
    statuses[:, 0] = np.asarray(beta0_status, dtype=object)
    for beta_index in range(1, 91):
        previous, current_status, current_values, root_seeds = _pilot_transition(
            previous,
            case,
            float(beta_index),
            model,
            step_type="beta",
            step_index=beta_index,
            previous_roots=root_seeds,
        )
        values[:, beta_index] = np.asarray(current_values, dtype=float)
        statuses[:, beta_index] = np.asarray(current_status, dtype=object)
    return {"model": model, "values": values, "statuses": statuses}


def _pilot_case_worker(case_tuple: tuple[str, float, float, float]) -> dict[str, object]:
    label, epsilon, mu, eta = case_tuple
    case = CaseSpec(mu=mu, eta=eta, epsilon=epsilon, label=label)
    print(
        f"[{label}] tracking epsilon_0={epsilon:.17g}, mu={mu:g}, eta={eta:g}",
        flush=True,
    )
    models = {
        model: _pilot_track_model(case, model)
        for model in MODELS
    }
    return {"case": case_tuple, "models": models}


def _pilot_plot(
    output_dir: Path,
    case: CaseSpec,
    eb_values: np.ndarray,
    timo_values: np.ndarray,
) -> Path:
    beta = np.arange(91, dtype=float)
    colors = plt.get_cmap("tab20")(np.linspace(0.0, 1.0, PILOT_N_TRACK))
    fig, ax = plt.subplots(figsize=(12.0, 7.4))
    for branch_zero in range(PILOT_N_TRACK):
        color = colors[branch_zero]
        ax.plot(beta, eb_values[branch_zero], color=color, lw=1.35, ls="-")
        ax.plot(beta, timo_values[branch_zero], color=color, lw=1.35, ls="--")
    ax.set_xlabel(r"$\beta$, deg")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_xlim(0.0, 90.0)
    ax.set_xticks(np.arange(0.0, 91.0, 10.0))
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_title(
        rf"{case.label}: $\epsilon_0={case.epsilon:.16g}$, "
        rf"$\mu={case.mu:g}$, $\eta={case.eta:g}$"
    )
    branch_handles = [
        Line2D([0], [0], color=colors[index], lw=2.0, label=f"branch {index + 1}")
        for index in range(PILOT_N_TRACK)
    ]
    theory_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="-", label="Euler–Bernoulli"),
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="Timoshenko"),
    ]
    theory_legend = ax.legend(handles=theory_handles, loc="upper left", fontsize=9, frameon=True)
    ax.add_artist(theory_legend)
    ax.legend(
        handles=branch_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=6,
        fontsize=8,
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    path = output_dir / f"{case.label}_lambda_beta_comparison.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def sorted_pilot_search_settings() -> COMPLETE.SearchSettings:
    """Fixed one-pass settings for the independent sorted-spectrum pilot."""

    return COMPLETE.SearchSettings(
        requested_roots=SORTED_PILOT_N_STORE,
        candidate_roots=SORTED_PILOT_N_STORE,
        verification_candidate_roots=SORTED_PILOT_N_STORE + 1,
        max_upper_growth_tries=1,
    )


def _sorted_pilot_identity(case_tuple: tuple[str, float, float, float], beta_deg: int) -> dict[str, object]:
    return {
        "cache_version": SORTED_PILOT_CACHE_VERSION,
        "general_spectrum_algorithm_version": COMPLETE.GENERAL_SPECTRUM_ALGORITHM_VERSION,
        "eb_matrix_evaluator_version": COMPLETE.EB_MATRIX_EVALUATOR_VERSION,
        "timoshenko_basis_evaluator_version": COMPLETE.TIMO.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
        "case": list(case_tuple),
        "beta_deg": int(beta_deg),
        "settings": asdict(sorted_pilot_search_settings()),
        "independent_pointwise": True,
        "seed_roots": [],
        "verification_enabled": False,
    }


def _sorted_pilot_cache_path(output_dir: Path, case_id: str, beta_deg: int) -> Path:
    return output_dir / "cache" / f"{case_id}_beta_{int(beta_deg):03d}.json"


def _load_sorted_pilot_cache(
    path: Path,
    identity: dict[str, object],
) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return None
    if payload.get("identity") != identity:
        return None
    point = payload.get("point")
    return point if isinstance(point, dict) else None


def _save_sorted_pilot_cache(
    path: Path,
    identity: dict[str, object],
    point: dict[str, object],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps({"identity": identity, "point": point}, sort_keys=True, indent=2),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _sorted_pilot_model_search(
    model: str,
    geometry: COMPLETE.Geometry,
) -> dict[str, object]:
    settings = sorted_pilot_search_settings()
    result = COMPLETE.resolve_primary_spectrum(model, geometry, settings=settings)
    roots = list(result.roots[:SORTED_PILOT_N_STORE])
    accepted_count = len(result.roots)
    relevant_upper = (
        roots[-1].Lambda + max(settings.seed_half_width, 2.0 * settings.scan_step)
        if roots
        else float("inf")
    )
    relevant_unresolved = []
    for entry in result.unresolved_intervals:
        try:
            interval_left = float(str(entry).split(":", 1)[0])
        except (TypeError, ValueError):
            interval_left = float("-inf")
        if interval_left <= relevant_upper:
            relevant_unresolved.append(str(entry))
    resolved_cluster = any(root.cluster_size > 1 for root in roots)
    if accepted_count < SORTED_PILOT_N_STORE:
        status = "incomplete_root_inventory"
    elif relevant_unresolved:
        status = "diagnostic_unresolved_interval"
    elif resolved_cluster:
        status = "complete_with_resolved_cluster"
    else:
        status = "complete_root_inventory"
    return {
        "values": [float(root.Lambda) for root in roots],
        "root_sources": ["+".join(root.detection_sources) for root in roots],
        "root_count": int(accepted_count),
        "candidate_root_count": int(len(result.candidates)),
        "point_status": status,
        "lambda_min": float(settings.lambda_min),
        "lambda_max": float(result.lambda_upper),
        "scan_step": float(result.scan_step),
        "matrix_evaluations": int(result.operations.characteristic_matrix_evaluations),
        "unresolved_intervals": relevant_unresolved,
        "resolved_cluster_count": int(sum(root.cluster_size > 1 for root in roots)),
        "independent_verification_runs": int(result.operations.independent_verification_runs),
    }


def _sorted_pilot_case_worker(
    case_tuple: tuple[str, float, float, float],
    output_dir_text: str,
) -> dict[str, object]:
    label, epsilon, mu, eta = case_tuple
    output_dir = Path(output_dir_text)
    points: list[dict[str, object]] = []
    cache_hits = 0
    for beta_deg in range(91):
        identity = _sorted_pilot_identity(case_tuple, beta_deg)
        cache_path = _sorted_pilot_cache_path(output_dir, label, beta_deg)
        point = _load_sorted_pilot_cache(cache_path, identity)
        if point is None:
            geometry = COMPLETE.Geometry(
                epsilon_0=float(epsilon),
                beta_deg=float(beta_deg),
                mu=float(mu),
                eta=float(eta),
            )
            point = {
                "case_id": label,
                "epsilon_0": float(epsilon),
                "mu": float(mu),
                "eta": float(eta),
                "beta_deg": int(beta_deg),
                "models": {
                    MODEL_EB: _sorted_pilot_model_search(MODEL_EB, geometry),
                    MODEL_TIMO: _sorted_pilot_model_search(MODEL_TIMO, geometry),
                },
            }
            _save_sorted_pilot_cache(cache_path, identity, point)
        else:
            cache_hits += 1
        points.append(point)
        if beta_deg % 15 == 0 or beta_deg == 90:
            print(f"[{label}] beta={beta_deg:02d}/90, cache_hits={cache_hits}", flush=True)
    return {"case": case_tuple, "points": points, "cache_hits": cache_hits}


def _sorted_pilot_values(
    case_result: dict[str, object],
    model: str,
) -> np.ndarray:
    values = np.full((SORTED_PILOT_N_STORE, 91), np.nan, dtype=float)
    for beta_index, point in enumerate(case_result["points"]):  # type: ignore[index]
        roots = point["models"][model]["values"]  # type: ignore[index]
        values[: min(len(roots), SORTED_PILOT_N_STORE), beta_index] = roots[:SORTED_PILOT_N_STORE]
    return values


def _sorted_pilot_plot(
    output_dir: Path,
    case: CaseSpec,
    eb_values: np.ndarray,
    timo_values: np.ndarray,
) -> Path:
    beta = np.arange(91, dtype=float)
    colors = plt.get_cmap("tab10")(np.arange(SORTED_PILOT_N_PLOT))
    fig, ax = plt.subplots(figsize=(12.4, 7.6))
    for mode_zero in range(SORTED_PILOT_N_PLOT):
        color = colors[mode_zero]
        ax.plot(beta, eb_values[mode_zero], color=color, lw=1.35, ls="-")
        ax.plot(beta, timo_values[mode_zero], color=color, lw=1.35, ls="--")
    ax.set_xlabel(r"$\beta$, deg")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_xlim(0.0, 90.0)
    ax.set_xticks(np.arange(0.0, 91.0, 10.0))
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_title(
        rf"{case.label}: sorted spectra, $\epsilon_0={case.epsilon:.16g}$, "
        rf"$\mu={case.mu:g}$, $\eta={case.eta:g}$"
    )
    theory_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="-", label="Euler-Bernoulli"),
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="Timoshenko"),
    ]
    mode_handles = [
        Line2D([0], [0], color=colors[index], lw=2.0, label=f"sorted mode {index + 1}")
        for index in range(SORTED_PILOT_N_PLOT)
    ]
    theory_legend = ax.legend(handles=theory_handles, loc="upper left", fontsize=9, frameon=True)
    ax.add_artist(theory_legend)
    ax.legend(
        handles=mode_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=5,
        fontsize=8,
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    path = output_dir / f"{case.label}_sorted_lambda_beta_comparison.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def _sorted_pilot_csv_rows(
    results: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label, epsilon, mu, eta in SORTED_PILOT_CASES:
        for point in results[label]["points"]:  # type: ignore[index]
            eb = point["models"][MODEL_EB]  # type: ignore[index]
            tm = point["models"][MODEL_TIMO]  # type: ignore[index]
            for mode_zero in range(SORTED_PILOT_N_STORE):
                eb_value = eb["values"][mode_zero] if mode_zero < len(eb["values"]) else float("nan")
                tm_value = tm["values"][mode_zero] if mode_zero < len(tm["values"]) else float("nan")
                eb_source = eb["root_sources"][mode_zero] if mode_zero < len(eb["root_sources"]) else "missing"
                tm_source = tm["root_sources"][mode_zero] if mode_zero < len(tm["root_sources"]) else "missing"
                rows.append(
                    {
                        "case_id": label,
                        "epsilon_0": epsilon,
                        "mu": mu,
                        "eta": eta,
                        "beta_deg": point["beta_deg"],
                        "sorted_mode_index": mode_zero + 1,
                        "lambda_eb": eb_value,
                        "lambda_timo": tm_value,
                        "eb_root_count": eb["root_count"],
                        "timo_root_count": tm["root_count"],
                        "eb_candidate_root_count": eb["candidate_root_count"],
                        "timo_candidate_root_count": tm["candidate_root_count"],
                        "eb_point_status": eb["point_status"],
                        "timo_point_status": tm["point_status"],
                        "eb_root_source": eb_source,
                        "timo_root_source": tm_source,
                        "eb_lambda_min": eb["lambda_min"],
                        "eb_lambda_max": eb["lambda_max"],
                        "timo_lambda_min": tm["lambda_min"],
                        "timo_lambda_max": tm["lambda_max"],
                        "scan_step": eb["scan_step"],
                        "eb_matrix_evaluations": eb["matrix_evaluations"],
                        "timo_matrix_evaluations": tm["matrix_evaluations"],
                    }
                )
    return rows


def _sorted_pilot_step_rows(
    results: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for label, _epsilon, _mu, _eta in SORTED_PILOT_CASES:
        case_result = results[label]
        for model in MODELS:
            theory = "EB" if model == MODEL_EB else "Timoshenko"
            values = _sorted_pilot_values(case_result, model)
            counts = [
                int(point["models"][model]["root_count"])  # type: ignore[index]
                for point in case_result["points"]  # type: ignore[index]
            ]
            statuses = [
                str(point["models"][model]["point_status"])  # type: ignore[index]
                for point in case_result["points"]  # type: ignore[index]
            ]
            unresolved_by_beta = [
                list(point["models"][model]["unresolved_intervals"])  # type: ignore[index]
                for point in case_result["points"]  # type: ignore[index]
            ]
            for beta_left in range(90):
                left_values = values[:, beta_left]
                right_values = values[:, beta_left + 1]
                jumps = np.abs(right_values - left_values)
                scale = np.maximum(np.maximum(np.abs(left_values), np.abs(right_values)), 1.0e-12)
                relative = jumps / scale
                directions = np.sign(right_values - left_values)
                noticeable = np.isfinite(jumps) & (
                    (relative >= SORTED_PILOT_NOTICEABLE_JUMP_REL)
                    | (jumps >= SORTED_PILOT_NOTICEABLE_JUMP_ABS)
                )
                positive_count = int(np.count_nonzero(noticeable & (directions > 0.0)))
                negative_count = int(np.count_nonzero(noticeable & (directions < 0.0)))
                incomplete = counts[beta_left] < SORTED_PILOT_N_STORE or counts[beta_left + 1] < SORTED_PILOT_N_STORE
                simultaneous = max(positive_count, negative_count)
                diagnostic_ranges: list[tuple[float, float]] = []
                for entry in (*unresolved_by_beta[beta_left], *unresolved_by_beta[beta_left + 1]):
                    try:
                        range_left, range_right, *_rest = str(entry).split(":")
                        diagnostic_ranges.append((float(range_left), float(range_right)))
                    except (TypeError, ValueError):
                        continue
                diagnostic_status = (
                    statuses[beta_left] == "diagnostic_unresolved_interval"
                    or statuses[beta_left + 1] == "diagnostic_unresolved_interval"
                )
                diagnostic_affected = {
                    mode_zero
                    for mode_zero in range(SORTED_PILOT_N_STORE)
                    if any(
                        interval_left - 0.5
                        <= value
                        <= interval_right + 0.5
                        for interval_left, interval_right in diagnostic_ranges
                        for value in (left_values[mode_zero], right_values[mode_zero])
                        if isfinite(float(value))
                    )
                }
                if diagnostic_status and not diagnostic_affected:
                    diagnostic_affected = {SORTED_PILOT_N_STORE - 2, SORTED_PILOT_N_STORE - 1}
                for mode_zero in range(SORTED_PILOT_N_STORE):
                    same_direction_count = 0
                    if noticeable[mode_zero] and directions[mode_zero] > 0.0:
                        same_direction_count = positive_count
                    elif noticeable[mode_zero] and directions[mode_zero] < 0.0:
                        same_direction_count = negative_count
                    reasons: list[str] = []
                    if incomplete:
                        reasons.append("incomplete_root_inventory")
                    if diagnostic_status and mode_zero in diagnostic_affected:
                        reasons.append("pointwise_unresolved_interval")
                    if simultaneous >= SORTED_PILOT_SIMULTANEOUS_COUNT:
                        reasons.append("simultaneous_sorted_mode_shift")
                    if np.isfinite(relative[mode_zero]) and relative[mode_zero] >= SORTED_PILOT_SINGLE_JUMP_REL:
                        reasons.append("sharp_single_mode_jump")
                    rows.append(
                        {
                            "case_id": label,
                            "theory": theory,
                            "sorted_mode_index": mode_zero + 1,
                            "beta_left": beta_left,
                            "beta_right": beta_left + 1,
                            "lambda_left": left_values[mode_zero],
                            "lambda_right": right_values[mode_zero],
                            "absolute_jump": jumps[mode_zero],
                            "relative_jump": relative[mode_zero],
                            "root_count_left": counts[beta_left],
                            "root_count_right": counts[beta_left + 1],
                            "simultaneous_shift_count": same_direction_count,
                            "suspect_interval": bool(reasons),
                            "suspect_reason": "+".join(reasons),
                        }
                    )
    return rows


def _sorted_pilot_suspect_rows(
    step_rows: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    grouped: dict[tuple[str, str, int, int], list[dict[str, object]]] = {}
    for row in step_rows:
        if not bool(row["suspect_interval"]):
            continue
        key = (str(row["case_id"]), str(row["theory"]), int(row["beta_left"]), int(row["beta_right"]))
        grouped.setdefault(key, []).append(row)

    interval_records: list[dict[str, object]] = []
    for (case_id, theory, beta_left, beta_right), rows in sorted(grouped.items()):
        affected: set[int] = set()
        reasons: set[str] = set()
        relative_values: list[float] = []
        lambda_values: list[float] = []
        count_changes: set[str] = set()
        for row in rows:
            reasons.update(filter(None, str(row["suspect_reason"]).split("+")))
            relative = float(row["relative_jump"])
            if isfinite(relative):
                relative_values.append(relative)
            left_value = float(row["lambda_left"])
            right_value = float(row["lambda_right"])
            if isfinite(left_value):
                lambda_values.append(left_value)
            if isfinite(right_value):
                lambda_values.append(right_value)
            if (
                not isfinite(left_value)
                or not isfinite(right_value)
                or relative >= SORTED_PILOT_SINGLE_JUMP_REL
                or int(row["simultaneous_shift_count"]) >= SORTED_PILOT_SIMULTANEOUS_COUNT
                or "pointwise_unresolved_interval" in str(row["suspect_reason"])
            ):
                affected.add(int(row["sorted_mode_index"]))
            left_count = int(row["root_count_left"])
            right_count = int(row["root_count_right"])
            if left_count != right_count or left_count < SORTED_PILOT_N_STORE:
                count_changes.add(f"{left_count}->{right_count}")
        interval_records.append(
            {
                "case_id": case_id,
                "theory": theory,
                "beta_left": beta_left,
                "beta_right": beta_right,
                "affected": affected,
                "reasons": reasons,
                "relative_values": relative_values,
                "lambda_values": lambda_values,
                "count_changes": count_changes,
            }
        )

    merged: list[dict[str, object]] = []
    for record in interval_records:
        if (
            merged
            and merged[-1]["case_id"] == record["case_id"]
            and merged[-1]["theory"] == record["theory"]
            and int(record["beta_left"]) <= int(merged[-1]["beta_right"])
        ):
            merged[-1]["beta_right"] = record["beta_right"]
            for key in ("affected", "reasons", "relative_values", "lambda_values", "count_changes"):
                if isinstance(merged[-1][key], set):
                    merged[-1][key].update(record[key])
                else:
                    merged[-1][key].extend(record[key])
        else:
            merged.append(
                {
                    **record,
                    "affected": set(record["affected"]),
                    "reasons": set(record["reasons"]),
                    "relative_values": list(record["relative_values"]),
                    "lambda_values": list(record["lambda_values"]),
                    "count_changes": set(record["count_changes"]),
                }
            )

    output: list[dict[str, object]] = []
    for record in merged:
        lambdas = list(record["lambda_values"])
        output.append(
            {
                "case_id": record["case_id"],
                "theory": record["theory"],
                "beta_left": record["beta_left"],
                "beta_right": record["beta_right"],
                "affected_sorted_modes": ",".join(str(item) for item in sorted(record["affected"])),
                "reason": "+".join(sorted(record["reasons"])),
                "max_relative_jump": max(record["relative_values"], default=float("nan")),
                "root_count_change": ",".join(sorted(record["count_changes"])) or "none",
                "suggested_local_beta_step": SORTED_PILOT_SUGGESTED_BETA_STEP,
                "suggested_lambda_interval_left": min(lambdas) if lambdas else float("nan"),
                "suggested_lambda_interval_right": max(lambdas) if lambdas else float("nan"),
            }
        )
    return output


def validate_sorted_pilot_artifacts(
    csv_rows: Sequence[dict[str, object]],
    png_paths: Sequence[Path],
    output_dir: Path,
) -> dict[str, object]:
    expected_rows = len(SORTED_PILOT_CASES) * 91 * SORTED_PILOT_N_STORE
    if len(csv_rows) != expected_rows:
        raise AssertionError(f"expected {expected_rows} CSV rows, found {len(csv_rows)}")
    for label, _epsilon, _mu, _eta in SORTED_PILOT_CASES:
        case_rows = [row for row in csv_rows if row["case_id"] == label]
        betas = sorted({int(row["beta_deg"]) for row in case_rows})
        if betas != list(range(91)):
            raise AssertionError(f"{label}: beta grid is not exactly 0..90")
        for beta_deg in betas:
            point_rows = [row for row in case_rows if int(row["beta_deg"]) == beta_deg]
            for field in ("lambda_eb", "lambda_timo"):
                finite = [float(row[field]) for row in point_rows if isfinite(float(row[field]))]
                if any(right <= left for left, right in zip(finite, finite[1:])):
                    raise AssertionError(f"{label}, beta={beta_deg}, {field}: roots are not strictly sorted")
                if any(right - left <= COMPLETE.DEFAULT_ROOT_DEDUP_TOL for left, right in zip(finite, finite[1:])):
                    raise AssertionError(f"{label}, beta={beta_deg}, {field}: duplicate roots exceed dedup policy")
    if len(png_paths) != 4 or len(list(output_dir.glob("*.png"))) != 4:
        raise AssertionError("sorted pilot must produce exactly four PNG files")
    forbidden = [*output_dir.glob("*.pdf"), *output_dir.glob("*.svg"), *output_dir.glob("*.eps")]
    if forbidden:
        raise AssertionError(f"forbidden vector/PDF outputs found: {forbidden}")
    return {
        "csv_rows": len(csv_rows),
        "png_count": len(png_paths),
        "force_strict_verification_calls": 0,
        "tracking_function_calls": 0,
        "mac_function_calls": 0,
        "shape_reconstruction_calls": 0,
        "adaptive_beta_refinement_calls": 0,
    }


def run_beta_sorted_spectrum_pilot() -> dict[str, object]:
    output_dir = _repo_output_dir(SORTED_PILOT_OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, object]] = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(_sorted_pilot_case_worker, case, str(output_dir)): case[0]
            for case in SORTED_PILOT_CASES
        }
        for future in as_completed(futures):
            result = future.result()
            label = str(result["case"][0])  # type: ignore[index]
            results[label] = result
            print(f"[{label}] independent sorted spectrum complete", flush=True)

    csv_rows = _sorted_pilot_csv_rows(results)
    csv_path = output_dir / "beta_sorted_spectrum_pilot.csv"
    _write_csv(csv_path, csv_rows, SORTED_PILOT_CSV_FIELDS)
    step_rows = _sorted_pilot_step_rows(results)
    diagnostics_path = output_dir / "beta_step_diagnostics.csv"
    _write_csv(diagnostics_path, step_rows, SORTED_PILOT_STEP_FIELDS)
    suspect_rows = _sorted_pilot_suspect_rows(step_rows)
    suspect_path = output_dir / "suspect_beta_intervals.csv"
    _write_csv(suspect_path, suspect_rows, SORTED_PILOT_SUSPECT_FIELDS)

    png_paths: list[Path] = []
    for label, epsilon, mu, eta in SORTED_PILOT_CASES:
        case = CaseSpec(mu=mu, eta=eta, epsilon=epsilon, label=label)
        png_paths.append(
            _sorted_pilot_plot(
                output_dir,
                case,
                _sorted_pilot_values(results[label], MODEL_EB),
                _sorted_pilot_values(results[label], MODEL_TIMO),
            )
        )

    full_eb = 0
    full_timo = 0
    fewer_than_12_eb = 0
    fewer_than_12_timo = 0
    nan_eb = 0
    nan_timo = 0
    eb_evaluations = 0
    timo_evaluations = 0
    verification_runs = 0
    status_counts: Counter[str] = Counter()
    for label, _epsilon, _mu, _eta in SORTED_PILOT_CASES:
        for point in results[label]["points"]:  # type: ignore[index]
            eb = point["models"][MODEL_EB]  # type: ignore[index]
            tm = point["models"][MODEL_TIMO]  # type: ignore[index]
            full_eb += int(
                int(eb["root_count"]) >= SORTED_PILOT_N_STORE
                and str(eb["point_status"]).startswith("complete_")
            )
            full_timo += int(
                int(tm["root_count"]) >= SORTED_PILOT_N_STORE
                and str(tm["point_status"]).startswith("complete_")
            )
            fewer_than_12_eb += int(int(eb["root_count"]) < SORTED_PILOT_N_STORE)
            fewer_than_12_timo += int(int(tm["root_count"]) < SORTED_PILOT_N_STORE)
            nan_eb += max(0, SORTED_PILOT_N_STORE - len(eb["values"]))
            nan_timo += max(0, SORTED_PILOT_N_STORE - len(tm["values"]))
            eb_evaluations += int(eb["matrix_evaluations"])
            timo_evaluations += int(tm["matrix_evaluations"])
            verification_runs += int(eb["independent_verification_runs"]) + int(tm["independent_verification_runs"])
            status_counts[f"EB:{eb['point_status']}"] += 1
            status_counts[f"Timoshenko:{tm['point_status']}"] += 1

    validation = validate_sorted_pilot_artifacts(csv_rows, png_paths, output_dir)
    if verification_runs != 0:
        raise AssertionError(f"primary-only pilot unexpectedly executed {verification_runs} verification runs")
    settings = sorted_pilot_search_settings()
    report_lines = [
        "# Independent beta sorted-spectrum pilot",
        "",
        "Curves show sorted spectral positions, not physical descendant branches.",
        "",
        "Each beta point was solved independently with the production EB/Timoshenko matrix evaluators and the ordinary two-phase pointwise primary root search. No neighboring-beta seeds or assignments were used.",
        "",
        "## Search configuration",
        "",
        f"- stored sorted modes: {SORTED_PILOT_N_STORE}",
        f"- plotted sorted modes: {SORTED_PILOT_N_PLOT}",
        f"- lambda_min: {settings.lambda_min:g}",
        "- initial and only lambda_max: 22",
        f"- scan_step: {settings.scan_step:g}",
        f"- grid phases: 0 and {settings.shifted_grid_phase:g}",
        "- upper-range growth attempts: disabled (one primary range only)",
        "- independent verification/full strict/force strict: disabled",
        "",
        "## Inventory summary",
        "",
        f"- complete EB beta-points: {full_eb}/364",
        f"- complete Timoshenko beta-points: {full_timo}/364",
        f"- EB beta-points with fewer than 12 roots: {fewer_than_12_eb}",
        f"- Timoshenko beta-points with fewer than 12 roots: {fewer_than_12_timo}",
        f"- EB missing-value NaNs: {nan_eb}",
        f"- Timoshenko missing-value NaNs: {nan_timo}",
        f"- EB matrix evaluations: {eb_evaluations}",
        f"- Timoshenko matrix evaluations: {timo_evaluations}",
        f"- cache hits: {sum(int(results[label]['cache_hits']) for label, *_ in SORTED_PILOT_CASES)}",
        f"- point statuses: {dict(sorted(status_counts.items()))}",
        "",
        "## Suspect-interval diagnostic",
        "",
        f"A jump is noticeable when relative jump >= {SORTED_PILOT_NOTICEABLE_JUMP_REL:g} or absolute jump >= {SORTED_PILOT_NOTICEABLE_JUMP_ABS:g}. An interval is flagged for an incomplete endpoint inventory, at least {SORTED_PILOT_SIMULTANEOUS_COUNT} noticeable same-direction sorted-mode shifts, or a single relative jump >= {SORTED_PILOT_SINGLE_JUMP_REL:g}. Adjacent flagged one-degree intervals are merged. These flags are diagnostic only; they do not assert crossing or veering and trigger no recalculation.",
        f"- merged suspect intervals: {len(suspect_rows)}",
        f"- suggested manual follow-up beta step: {SORTED_PILOT_SUGGESTED_BETA_STEP:g} deg",
        "",
        "## Prohibited-operation counters",
        "",
        "- branch tracking calls: 0",
        "- MAC calls: 0",
        "- shape reconstruction calls: 0",
        "- force_strict_verification calls: 0",
        "- independent verification calls: 0",
        "- adaptive beta refinement calls: 0",
        "",
        "## Outputs",
        "",
    ]
    report_lines.extend(f"- `{_rel(path)}`" for path in png_paths)
    report_lines.extend(
        [
            f"- `{_rel(csv_path)}`",
            f"- `{_rel(diagnostics_path)}`",
            f"- `{_rel(suspect_path)}`",
            f"- `{_rel(output_dir / 'report.md')}`",
            "",
        ]
    )
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    return {
        "png_paths": png_paths,
        "csv_path": csv_path,
        "diagnostics_path": diagnostics_path,
        "suspect_path": suspect_path,
        "report_path": report_path,
        "suspect_rows": suspect_rows,
        "full_eb": full_eb,
        "full_timo": full_timo,
        "fewer_than_12_eb": fewer_than_12_eb,
        "fewer_than_12_timo": fewer_than_12_timo,
        "nan_eb": nan_eb,
        "nan_timo": nan_timo,
        "eb_evaluations": eb_evaluations,
        "timo_evaluations": timo_evaluations,
        "validation": validation,
    }


def refined_pilot_regions() -> tuple[LocalRefinementRegion, ...]:
    return (
        LocalRefinementRegion("R1", "P1", MODEL_TIMO, 58.0, 60.0, 0.25, 9.8, 10.3, 1.0e-4, 2, "recover close pair near beta=59"),
        LocalRefinementRegion("R2", "P1", MODEL_EB, 87.0, 89.0, 0.25, 12.0, 12.4, 1.0e-4, 2, "recover close pair near beta=88"),
        LocalRefinementRegion("R3", "P3", MODEL_TIMO, 0.0, 3.0, 0.1, 7.9, 8.5, 1.0e-4, 2, "resolve beta=0 block-family pair"),
        LocalRefinementRegion("R4", "P4", MODEL_TIMO, 56.0, 58.0, 0.25, 3.4, 4.1, 1.0e-4, 2, "recover low pair near beta=57"),
        LocalRefinementRegion("R5", "P4", MODEL_EB, 58.0, 60.0, 0.25, 5.2, 5.8, 1.0e-4, 2, "recover close pair near beta=59"),
        LocalRefinementRegion("R6", "P4", MODEL_TIMO, 71.0, 74.0, 0.25, 9.2, 9.8, 1.0e-4, 1, "recover root or close pair near beta=72-73"),
        LocalRefinementRegion("R7", "P4", MODEL_TIMO, 84.0, 86.0, 0.25, 5.0, 5.7, 1.0e-4, 2, "recover low pair near beta=85"),
    )


def _refined_case_tuple(case_id: str) -> tuple[str, float, float, float]:
    for case in SORTED_PILOT_CASES:
        if case[0] == case_id:
            return case
    raise KeyError(case_id)


def _region_beta_values(region: LocalRefinementRegion) -> tuple[float, ...]:
    return tuple(float(value) for value in regular_grid(region.beta_min, region.beta_max, region.beta_step))


def _is_integer_beta(beta_deg: float) -> bool:
    return abs(float(beta_deg) - round(float(beta_deg))) <= 1.0e-10


def _model_csv_prefix(model: str) -> str:
    return "eb" if model == MODEL_EB else "timo"


def _load_original_sorted_inventory(source_dir: Path) -> dict[tuple[str, int, str], dict[str, object]]:
    path = source_dir / "beta_sorted_spectrum_pilot.csv"
    grouped: dict[tuple[str, int, str], list[dict[str, str]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            case_id = str(row["case_id"])
            beta_deg = int(round(float(row["beta_deg"])))
            for model in MODELS:
                grouped.setdefault((case_id, beta_deg, model), []).append(row)
    result: dict[tuple[str, int, str], dict[str, object]] = {}
    for key, rows in grouped.items():
        model = key[2]
        prefix = _model_csv_prefix(model)
        ordered = sorted(rows, key=lambda row: int(row["sorted_mode_index"]))
        entries: list[dict[str, object]] = []
        for row in ordered:
            value = float(row[f"lambda_{prefix}"])
            entries.append(
                {
                    "value": value,
                    "source": str(row[f"{prefix}_root_source"]),
                    "multiplicity": 1,
                    "multiplicity_source": "original_scalar_inventory",
                    "block_family": "full_matrix",
                    "nullity": 1,
                    "repeated_root_slot": 1,
                }
            )
        result[key] = {
            "entries": entries,
            "status": str(ordered[0][f"{prefix}_point_status"]),
            "root_count": int(ordered[0][f"{prefix}_root_count"]),
            "inventory_source": "original_pointwise_primary",
        }
    return result


def _local_refinement_settings(region: LocalRefinementRegion) -> COMPLETE.SearchSettings:
    return replace(
        sorted_pilot_search_settings(),
        lambda_min=float(region.lambda_min),
        lambda_max=float(region.lambda_max),
        scan_step=float(region.lambda_step),
        max_upper_growth_tries=1,
    )


def _candidate_nullity(candidate: COMPLETE.RootCandidate, settings: COMPLETE.SearchSettings) -> int:
    diagnostics = candidate.diagnostics
    return 2 if (
        diagnostics.sigma_2 <= settings.nullity_sigma
        and isfinite(diagnostics.sigma_3)
        and diagnostics.sigma_3 > 0.0
        and diagnostics.sigma_2 / diagnostics.sigma_3 <= settings.sigma_ratio_accept
    ) else 1


def _root_entries_from_records(
    roots: Sequence[COMPLETE.RootRecord],
    *,
    source: str,
    block_family: str = "full_6x6",
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    slot_counter: dict[tuple[float, int], int] = {}
    for root in roots:
        multiplicity = max(1, int(root.detected_nullity))
        key = (round(float(root.Lambda), 10), multiplicity)
        slot_counter[key] = slot_counter.get(key, 0) + 1
        entries.append(
            {
                "value": float(root.Lambda),
                "source": source + "+" + "+".join(root.detection_sources),
                "multiplicity": multiplicity,
                "multiplicity_source": (
                    "full_matrix_nullity_2" if multiplicity > 1 else "simple_full_matrix_root"
                ),
                "block_family": block_family,
                "nullity": multiplicity,
                "repeated_root_slot": slot_counter[key] if multiplicity > 1 else 1,
            }
        )
    return entries


def _dense_local_candidate_search(
    provider,
    settings: COMPLETE.SearchSettings,
    *,
    candidate_source_prefix: str,
    block_family: str,
) -> dict[str, object]:
    operations = COMPLETE.OperationCounts()
    evaluator = COMPLETE._MatrixEvaluator(provider, operations)
    candidates, interval_rows, unresolved = COMPLETE._global_candidates(
        evaluator,
        settings,
        configuration="primary",
        scan_step=settings.scan_step,
        phases=(0.0, settings.shifted_grid_phase),
        upper=float(settings.lambda_max),
        seed_roots=(),
        seed_source="local_refinement_no_seeds",
    )
    merged = COMPLETE._merge_candidates(candidates, settings)
    roots = COMPLETE._root_records(merged, settings)
    entries = _root_entries_from_records(
        roots,
        source=candidate_source_prefix,
        block_family=block_family,
    )
    candidate_rows: list[dict[str, object]] = []
    for candidate in candidates:
        nullity = _candidate_nullity(candidate, settings)
        accepted = candidate.acceptance_status == "accepted_full_matrix_svd"
        candidate_rows.append(
            {
                "lambda_candidate": float(candidate.Lambda),
                "candidate_source": candidate_source_prefix + ":" + "+".join(candidate.detection_sources),
                "sign_change": any("sign_change" in item or "grid_zero" in item for item in candidate.detection_sources),
                "sigma_min": float(candidate.diagnostics.sigma_1),
                "sigma_ratio": float(candidate.diagnostics.sigma_ratio),
                "residual": float(candidate.diagnostics.sigma_1),
                "bracket_left": float(candidate.interval_left),
                "bracket_right": float(candidate.interval_right),
                "multiplicity": nullity,
                "accepted": accepted,
                "rejection_reason": "" if accepted else candidate.acceptance_status,
                "block_family": block_family,
                "nullity": nullity,
                "multiplicity_source": "full_matrix_SVD" if block_family == "full_6x6" else "block_SVD",
            }
        )
    return {
        "entries": entries,
        "candidate_rows": candidate_rows,
        "interval_rows": list(interval_rows),
        "unresolved_intervals": list(unresolved),
        "matrix_evaluations": int(operations.characteristic_matrix_evaluations),
        "operations": asdict(operations),
    }


def _annotate_p3_beta0_block_provenance(
    full_result: dict[str, object],
    bending_result: dict[str, object],
    axial_result: dict[str, object],
    settings: COMPLETE.SearchSettings,
) -> tuple[list[dict[str, object]], str]:
    bending_values = [float(entry["value"]) for entry in bending_result["entries"]]  # type: ignore[index]
    axial_values = [float(entry["value"]) for entry in axial_result["entries"]]  # type: ignore[index]
    annotated: list[dict[str, object]] = []
    for raw_entry in full_result["entries"]:  # type: ignore[index]
        entry = dict(raw_entry)
        value = float(entry["value"])
        families: list[str] = []
        if any(abs(value - item) <= settings.root_match_tol for item in bending_values):
            families.append("bending_block")
        if any(abs(value - item) <= settings.root_match_tol for item in axial_values):
            families.append("axial_block")
        entry["block_family"] = "+".join(families) if families else "full_6x6_only"
        if int(entry["multiplicity"]) > 1 and len(families) == 2:
            entry["multiplicity_source"] = "full_6x6_nullity_2+cross_block_coincidence"
        elif families:
            entry["multiplicity_source"] = "simple_full_matrix_root+" + "+".join(families)
        annotated.append(entry)
    has_confirmed_double = any(
        int(entry["multiplicity"]) > 1
        and entry["block_family"] == "bending_block+axial_block"
        for entry in annotated
    )
    family_values = sorted(
        float(entry["value"])
        for entry in annotated
        if entry["block_family"] in {"bending_block", "axial_block", "bending_block+axial_block"}
    )
    has_distinct_pair = any(
        right - left > settings.root_match_tol
        for left, right in zip(family_values, family_values[1:])
    )
    if has_confirmed_double:
        status = "resolved_double_root"
    elif has_distinct_pair:
        status = "resolved_distinct_pair"
    else:
        status = "local_refinement_unresolved"
    return annotated, status


def merge_primary_and_local_inventory(
    primary_entries: Sequence[dict[str, object]],
    local_entries: Sequence[dict[str, object]],
    *,
    lambda_min: float,
    lambda_max: float,
    settings: COMPLETE.SearchSettings,
) -> list[dict[str, object]]:
    kept = [
        dict(entry)
        for entry in primary_entries
        if not (
            isfinite(float(entry["value"]))
            and float(lambda_min) - settings.root_match_tol
            <= float(entry["value"])
            <= float(lambda_max) + settings.root_match_tol
        )
    ]
    combined = kept + [dict(entry) for entry in local_entries]
    combined.sort(key=lambda entry: (float(entry["value"]), int(entry.get("repeated_root_slot", 1))))
    deduplicated: list[dict[str, object]] = []
    for entry in combined:
        if not deduplicated:
            deduplicated.append(entry)
            continue
        previous = deduplicated[-1]
        gap = abs(float(entry["value"]) - float(previous["value"]))
        preserve_slots = (
            int(entry.get("multiplicity", 1)) > 1
            and int(previous.get("multiplicity", 1)) > 1
            and int(entry.get("repeated_root_slot", 1)) != int(previous.get("repeated_root_slot", 1))
        )
        if gap <= settings.root_dedup_tol and not preserve_slots:
            previous_local = "local_dense" in str(previous.get("source", ""))
            entry_local = "local_dense" in str(entry.get("source", ""))
            if entry_local and not previous_local:
                deduplicated[-1] = entry
            continue
        deduplicated.append(entry)
    return deduplicated[:SORTED_PILOT_N_STORE]


def merge_confirmed_inventory_union(
    primary_entries: Sequence[dict[str, object]],
    recovered_entries: Sequence[dict[str, object]],
    *,
    settings: COMPLETE.SearchSettings,
) -> list[dict[str, object]]:
    combined = [dict(entry) for entry in primary_entries] + [dict(entry) for entry in recovered_entries]
    combined.sort(key=lambda entry: (float(entry["value"]), int(entry.get("repeated_root_slot", 1))))
    result: list[dict[str, object]] = []
    for entry in combined:
        if not result:
            result.append(entry)
            continue
        previous = result[-1]
        gap = abs(float(entry["value"]) - float(previous["value"]))
        preserve_slots = (
            int(entry.get("multiplicity", 1)) > 1
            and int(previous.get("multiplicity", 1)) > 1
            and int(entry.get("repeated_root_slot", 1)) != int(previous.get("repeated_root_slot", 1))
        )
        if gap <= settings.root_dedup_tol and not preserve_slots:
            if "base_candidate_local_dense" in str(entry.get("source", "")):
                result[-1] = entry
            continue
        result.append(entry)
    return result


def unresolved_inventory_from_position(
    primary_entries: Sequence[dict[str, object]],
    *,
    lambda_min: float,
) -> list[dict[str, object]]:
    entries = [dict(entry) for entry in primary_entries[:SORTED_PILOT_N_STORE]]
    first = next(
        (
            index
            for index, entry in enumerate(entries)
            if isfinite(float(entry["value"])) and float(entry["value"]) >= float(lambda_min)
        ),
        len(entries),
    )
    for index in range(first, len(entries)):
        entries[index] = {
            "value": float("nan"),
            "source": "local_refinement_unresolved",
            "multiplicity": 0,
            "multiplicity_source": "unresolved",
            "block_family": "unresolved",
            "nullity": 0,
            "repeated_root_slot": 0,
        }
    return entries


def _refined_cache_path(output_dir: Path, region_id: str, beta_deg: float) -> Path:
    beta_token = f"{float(beta_deg):08.3f}".replace(".", "p")
    return output_dir / "cache" / region_id / f"beta_{beta_token}.json"


def _refined_cache_identity(
    region: LocalRefinementRegion,
    beta_deg: float,
    source_csv_sha256: str,
) -> dict[str, object]:
    settings = _local_refinement_settings(region)
    return {
        "cache_version": REFINED_PILOT_CACHE_VERSION,
        "algorithm_version": REFINED_PILOT_ALGORITHM_VERSION,
        "region": asdict(region),
        "beta_deg": float(beta_deg),
        "source_csv_sha256": source_csv_sha256,
        "general_spectrum_algorithm_version": COMPLETE.GENERAL_SPECTRUM_ALGORITHM_VERSION,
        "eb_matrix_evaluator_version": COMPLETE.EB_MATRIX_EVALUATOR_VERSION,
        "timoshenko_basis_evaluator_version": COMPLETE.TIMO.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
        "settings": asdict(settings),
        "force_strict_enabled": False,
        "tracking_enabled": False,
    }


def _load_refined_cache(path: Path, identity: dict[str, object]) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, TypeError, ValueError):
        return None
    point = payload.get("point")
    return point if payload.get("identity") == identity and isinstance(point, dict) else None


def _save_refined_cache(path: Path, identity: dict[str, object], point: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(
        json.dumps({"identity": identity, "point": point}, sort_keys=True, indent=2),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _primary_entries_for_fractional_beta(
    region: LocalRefinementRegion,
    beta_deg: float,
    source_inventory: dict[tuple[str, int, str], dict[str, object]],
) -> tuple[list[dict[str, object]], str, int, int, list[dict[str, object]]]:
    _label, epsilon, mu, eta = _refined_case_tuple(region.case_id)
    geometry = COMPLETE.Geometry(epsilon_0=epsilon, beta_deg=beta_deg, mu=mu, eta=eta)
    primary = COMPLETE.resolve_primary_spectrum(
        region.model,
        geometry,
        settings=sorted_pilot_search_settings(),
    )
    entries = _root_entries_from_records(
        primary.roots[:SORTED_PILOT_N_STORE],
        source="fractional_pointwise_primary",
    )
    base_candidate_rows: list[dict[str, object]] = []
    base_candidate_evaluations = 0
    recovered_entries: list[dict[str, object]] = []
    beta_left = int(np.floor(float(beta_deg)))
    beta_right = int(np.ceil(float(beta_deg)))
    endpoint_entries = [
        source_inventory[(region.case_id, beta_left, region.model)]["entries"],
        source_inventory[(region.case_id, beta_right, region.model)]["entries"],
    ]
    windows: list[tuple[float, float]] = []
    base_settings = sorted_pilot_search_settings()
    for mode_zero in range(SORTED_PILOT_N_STORE):
        endpoint_values = [
            float(endpoint[mode_zero]["value"])
            for endpoint in endpoint_entries
            if mode_zero < len(endpoint) and isfinite(float(endpoint[mode_zero]["value"]))
        ]
        if not endpoint_values:
            continue
        center_left = min(endpoint_values)
        center_right = max(endpoint_values)
        if center_right >= region.lambda_min - base_settings.root_match_tol and center_left <= region.lambda_max + base_settings.root_match_tol:
            continue
        windows.append(
            (
                max(base_settings.lambda_min, center_left - base_settings.seed_half_width),
                min(22.0, center_right + base_settings.seed_half_width),
            )
        )
    merged_windows: list[list[float]] = []
    for left, right in sorted(windows):
        if merged_windows and left <= merged_windows[-1][1] + base_settings.root_match_tol:
            merged_windows[-1][1] = max(merged_windows[-1][1], right)
        else:
            merged_windows.append([left, right])
    for window_left, window_right in merged_windows:
        window_settings = replace(
            base_settings,
            lambda_min=float(window_left),
            lambda_max=float(window_right),
            scan_step=1.0e-3,
            max_upper_growth_tries=1,
        )
        recovered = _dense_local_candidate_search(
            COMPLETE.model_matrix_provider(region.model, geometry),
            window_settings,
            candidate_source_prefix="fractional_base_candidate_local_dense",
            block_family="full_6x6",
        )
        recovered_entries.extend(recovered["entries"])
        base_candidate_rows.extend(recovered["candidate_rows"])
        base_candidate_evaluations += int(recovered["matrix_evaluations"])
    entries = merge_confirmed_inventory_union(entries, recovered_entries, settings=base_settings)
    status = "fractional_primary_complete" if len(entries) >= SORTED_PILOT_N_STORE else "fractional_primary_incomplete"
    return (
        entries,
        status,
        len(entries),
        int(primary.operations.characteristic_matrix_evaluations) + base_candidate_evaluations,
        base_candidate_rows,
    )


def _run_refinement_point(
    region: LocalRefinementRegion,
    beta_deg: float,
    source_inventory: dict[tuple[str, int, str], dict[str, object]],
) -> dict[str, object]:
    label, epsilon, mu, eta = _refined_case_tuple(region.case_id)
    settings = _local_refinement_settings(region)
    if _is_integer_beta(beta_deg):
        original = source_inventory[(region.case_id, int(round(beta_deg)), region.model)]
        primary_entries = [dict(entry) for entry in original["entries"]]  # type: ignore[index]
        primary_status = str(original["status"])
        primary_root_count = int(original["root_count"])
        base_primary_evaluations = 0
        base_candidate_rows: list[dict[str, object]] = []
        primary_source = "original_pointwise_primary"
    else:
        primary_entries, primary_status, primary_root_count, base_primary_evaluations, base_candidate_rows = (
            _primary_entries_for_fractional_beta(region, beta_deg, source_inventory)
        )
        primary_source = "fractional_pointwise_primary+base_candidate_matrix_recovery"

    geometry = COMPLETE.Geometry(epsilon_0=epsilon, beta_deg=beta_deg, mu=mu, eta=eta)
    provider = COMPLETE.model_matrix_provider(region.model, geometry)
    full = _dense_local_candidate_search(
        provider,
        settings,
        candidate_source_prefix=f"{region.region_id}_local_dense",
        block_family="full_6x6",
    )
    all_candidate_rows = list(base_candidate_rows) + list(full["candidate_rows"])
    local_matrix_evaluations = int(full["matrix_evaluations"])
    local_entries = [dict(entry) for entry in full["entries"]]
    special_status = ""
    block_details: dict[str, object] = {}
    if region.region_id == "R3" and abs(beta_deg) <= 1.0e-12:
        def bending_provider(value: float) -> np.ndarray:
            return provider(value)[np.ix_(STRAIGHT.BENDING_ROWS, STRAIGHT.BENDING_COLUMNS)]

        def axial_provider(value: float) -> np.ndarray:
            return provider(value)[np.ix_(STRAIGHT.AXIAL_ROWS, STRAIGHT.AXIAL_COLUMNS)]

        bending = _dense_local_candidate_search(
            bending_provider,
            settings,
            candidate_source_prefix="R3_bending_block_local_dense",
            block_family="bending_block",
        )
        axial = _dense_local_candidate_search(
            axial_provider,
            settings,
            candidate_source_prefix="R3_axial_block_local_dense",
            block_family="axial_block",
        )
        local_entries, special_status = _annotate_p3_beta0_block_provenance(full, bending, axial, settings)
        all_candidate_rows.extend(bending["candidate_rows"])
        all_candidate_rows.extend(axial["candidate_rows"])
        local_matrix_evaluations += int(bending["matrix_evaluations"]) + int(axial["matrix_evaluations"])
        block_details = {
            "bending_entries": bending["entries"],
            "axial_entries": axial["entries"],
            "bending_unresolved": bending["unresolved_intervals"],
            "axial_unresolved": axial["unresolved_intervals"],
        }

    local_slot_count = len(local_entries)
    unique_local_values: list[float] = []
    for entry in local_entries:
        value = float(entry["value"])
        if not unique_local_values or abs(value - unique_local_values[-1]) > settings.root_dedup_tol:
            unique_local_values.append(value)
    unresolved = bool(full["unresolved_intervals"])
    resolved = (
        len(primary_entries) >= SORTED_PILOT_N_STORE
        and local_slot_count >= region.expected_min_roots
        and not unresolved
        and special_status != "local_refinement_unresolved"
    )
    if resolved:
        merged = merge_primary_and_local_inventory(
            primary_entries,
            local_entries,
            lambda_min=region.lambda_min,
            lambda_max=region.lambda_max,
            settings=settings,
        )
        resolved = len(merged) >= SORTED_PILOT_N_STORE
    else:
        merged = []
    if not resolved:
        merged = unresolved_inventory_from_position(primary_entries, lambda_min=region.lambda_min)
        point_status = "local_refinement_unresolved"
    elif special_status:
        point_status = special_status
    elif any(int(entry["multiplicity"]) > 1 for entry in local_entries):
        point_status = "resolved_double_root"
    elif len(unique_local_values) >= 2:
        point_status = "resolved_distinct_pair"
    else:
        point_status = "resolved_local_root"

    base_window_count = sum(
        region.lambda_min - settings.root_match_tol
        <= float(entry["value"])
        <= region.lambda_max + settings.root_match_tol
        for entry in primary_entries
        if isfinite(float(entry["value"]))
    )
    multiplicity_two_groups = len(
        {
            round(float(entry["value"]), 10)
            for entry in local_entries
            if int(entry["multiplicity"]) > 1
        }
    )
    return {
        "region_id": region.region_id,
        "case_id": label,
        "model": region.model,
        "beta_deg": float(beta_deg),
        "beta_source": "original_integer" if _is_integer_beta(beta_deg) else "local_refinement",
        "primary_entries": primary_entries,
        "primary_status": primary_status,
        "primary_root_count": primary_root_count,
        "primary_source": primary_source,
        "entries": merged,
        "status": point_status,
        "root_count": sum(isfinite(float(entry["value"])) for entry in merged),
        "inventory_source": primary_source + f"+{region.region_id}_local_dense",
        "candidate_rows": all_candidate_rows,
        "unresolved_intervals": full["unresolved_intervals"],
        "local_matrix_evaluations": local_matrix_evaluations,
        "base_primary_matrix_evaluations": base_primary_evaluations,
        "recovered_root_count": max(0, local_slot_count - base_window_count) if resolved else 0,
        "multiplicity_two_count": multiplicity_two_groups,
        "block_details": block_details,
        "force_strict_calls": 0,
        "tracking_calls": 0,
        "mac_calls": 0,
        "shape_calls": 0,
        "continuation_calls": 0,
    }


def _refinement_region_worker(
    region: LocalRefinementRegion,
    source_dir_text: str,
    output_dir_text: str,
    source_csv_sha256: str,
) -> dict[str, object]:
    source_dir = Path(source_dir_text)
    output_dir = Path(output_dir_text)
    source_inventory = _load_original_sorted_inventory(source_dir)
    points: list[dict[str, object]] = []
    cache_hits = 0
    for beta_deg in _region_beta_values(region):
        identity = _refined_cache_identity(region, beta_deg, source_csv_sha256)
        path = _refined_cache_path(output_dir, region.region_id, beta_deg)
        point = _load_refined_cache(path, identity)
        if point is None:
            point = _run_refinement_point(region, beta_deg, source_inventory)
            _save_refined_cache(path, identity, point)
        else:
            cache_hits += 1
        points.append(point)
        print(
            f"[{region.region_id}] beta={beta_deg:g}, status={point['status']}, cache_hits={cache_hits}",
            flush=True,
        )
    return {
        "region": asdict(region),
        "points": points,
        "cache_hits": cache_hits,
    }


def _refinement_point_index(
    region_results: dict[str, dict[str, object]],
) -> dict[tuple[str, str, float], dict[str, object]]:
    index: dict[tuple[str, str, float], dict[str, object]] = {}
    for region in refined_pilot_regions():
        for point in region_results[region.region_id]["points"]:  # type: ignore[index]
            key = (region.case_id, region.model, round(float(point["beta_deg"]), 10))
            if key in index:
                raise AssertionError(f"duplicate local refinement target: {key}")
            index[key] = point
    return index


def _regions_at_case_beta(case_id: str, beta_deg: float) -> list[str]:
    result = []
    beta_key = round(float(beta_deg), 10)
    for region in refined_pilot_regions():
        if region.case_id == case_id and beta_key in {round(value, 10) for value in _region_beta_values(region)}:
            result.append(region.region_id)
    return result


def _assemble_refined_rows(
    source_inventory: dict[tuple[str, int, str], dict[str, object]],
    region_results: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    point_index = _refinement_point_index(region_results)
    rows: list[dict[str, object]] = []
    for case_id, epsilon, mu, eta in SORTED_PILOT_CASES:
        beta_values = {float(value) for value in range(91)}
        for region in refined_pilot_regions():
            if region.case_id == case_id:
                beta_values.update(_region_beta_values(region))
        for beta_deg in sorted(beta_values):
            model_inventories: dict[str, dict[str, object] | None] = {}
            for model in MODELS:
                local = point_index.get((case_id, model, round(beta_deg, 10)))
                if local is not None:
                    model_inventories[model] = {
                        "entries": local["entries"],
                        "status": local["status"],
                        "root_count": local["root_count"],
                        "inventory_source": local["inventory_source"],
                    }
                elif _is_integer_beta(beta_deg):
                    model_inventories[model] = source_inventory[(case_id, int(round(beta_deg)), model)]
                else:
                    model_inventories[model] = None
            region_ids = "+".join(_regions_at_case_beta(case_id, beta_deg))
            for mode_zero in range(SORTED_PILOT_N_STORE):
                row: dict[str, object] = {
                    "case_id": case_id,
                    "epsilon_0": epsilon,
                    "mu": mu,
                    "eta": eta,
                    "beta_deg": beta_deg,
                    "beta_source": "original_integer" if _is_integer_beta(beta_deg) else "local_refinement",
                    "sorted_mode_index": mode_zero + 1,
                    "local_region_id": region_ids,
                }
                for model in MODELS:
                    prefix = _model_csv_prefix(model)
                    inventory = model_inventories[model]
                    if inventory is None:
                        entry = None
                        status = "not_evaluated_at_fractional_beta"
                        root_count = 0
                        inventory_source = "original_integer_series_only"
                    else:
                        entries = inventory["entries"]  # type: ignore[index]
                        entry = entries[mode_zero] if mode_zero < len(entries) else None
                        status = str(inventory["status"])
                        root_count = int(inventory["root_count"])
                        inventory_source = str(inventory["inventory_source"])
                    row[f"lambda_{prefix}"] = float(entry["value"]) if entry is not None else float("nan")
                    row[f"{prefix}_status"] = status
                    row[f"{prefix}_root_count"] = root_count
                    row[f"{prefix}_inventory_source"] = inventory_source
                    row[f"{prefix}_multiplicity"] = int(entry["multiplicity"]) if entry is not None else 0
                    row[f"{prefix}_multiplicity_source"] = (
                        str(entry["multiplicity_source"]) if entry is not None else "not_evaluated"
                    )
                    row[f"{prefix}_block_family"] = str(entry["block_family"]) if entry is not None else "not_evaluated"
                    row[f"{prefix}_nullity"] = int(entry["nullity"]) if entry is not None else 0
                    row[f"{prefix}_repeated_root_slot"] = (
                        int(entry["repeated_root_slot"]) if entry is not None else 0
                    )
                rows.append(row)
    return rows


def _candidate_output_rows(
    region_results: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for region in refined_pilot_regions():
        theory = "Euler-Bernoulli" if region.model == MODEL_EB else "Timoshenko"
        points = sorted(
            region_results[region.region_id]["points"],  # type: ignore[index]
            key=lambda point: float(point["beta_deg"]),
        )
        for point in points:
            candidates = sorted(
                point["candidate_rows"],
                key=lambda row: (float(row["lambda_candidate"]), str(row["candidate_source"])),
            )
            for candidate in candidates:
                rows.append(
                    {
                        "region_id": region.region_id,
                        "case_id": region.case_id,
                        "theory": theory,
                        "beta_deg": point["beta_deg"],
                        **candidate,
                    }
                )
    return rows


def _before_after_rows(
    source_inventory: dict[tuple[str, int, str], dict[str, object]],
    region_results: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    settings_by_region = {region.region_id: _local_refinement_settings(region) for region in refined_pilot_regions()}
    for region in refined_pilot_regions():
        settings = settings_by_region[region.region_id]
        for point in sorted(
            region_results[region.region_id]["points"],  # type: ignore[index]
            key=lambda item: float(item["beta_deg"]),
        ):
            beta_deg = float(point["beta_deg"])
            if not _is_integer_beta(beta_deg):
                continue
            before = source_inventory[(region.case_id, int(round(beta_deg)), region.model)]
            before_entries = before["entries"]  # type: ignore[index]
            after_entries = point["entries"]
            changed_indices: list[int] = []
            for index in range(SORTED_PILOT_N_STORE):
                before_value = float(before_entries[index]["value"])
                after_value = float(after_entries[index]["value"])
                if not (isfinite(before_value) and isfinite(after_value)) or abs(after_value - before_value) > settings.root_match_tol:
                    changed_indices.append(index)
            first_changed = min(changed_indices) if changed_indices else None
            for index in range(SORTED_PILOT_N_STORE):
                before_value = float(before_entries[index]["value"])
                after_value = float(after_entries[index]["value"])
                difference = (
                    abs(after_value - before_value)
                    if isfinite(before_value) and isfinite(after_value)
                    else float("nan")
                )
                if first_changed is None:
                    shift = "unchanged"
                elif index < first_changed:
                    shift = "unchanged_prefix"
                elif isfinite(difference) and difference <= settings.root_match_tol:
                    shift = "retained_after_reindex"
                else:
                    shift = "shifted_after_local_recovery"
                rows.append(
                    {
                        "case_id": region.case_id,
                        "theory": "Euler-Bernoulli" if region.model == MODEL_EB else "Timoshenko",
                        "beta_deg": beta_deg,
                        "sorted_mode_index": index + 1,
                        "lambda_before": before_value,
                        "lambda_after": after_value,
                        "absolute_difference": difference,
                        "status_before": before["status"],
                        "status_after": point["status"],
                        "root_inventory_shift": shift,
                    }
                )
    return rows


def _inventory_is_sorted(entries: Sequence[dict[str, object]]) -> bool:
    finite = [float(entry["value"]) for entry in entries if isfinite(float(entry["value"]))]
    return all(right >= left for left, right in zip(finite, finite[1:]))


def _inventory_has_invalid_duplicate(
    entries: Sequence[dict[str, object]],
    settings: COMPLETE.SearchSettings,
) -> bool:
    finite = [entry for entry in entries if isfinite(float(entry["value"]))]
    for left, right in zip(finite, finite[1:]):
        if float(right["value"]) - float(left["value"]) <= settings.root_dedup_tol:
            valid_multiplicity = (
                int(left.get("multiplicity", 1)) > 1
                and int(right.get("multiplicity", 1)) > 1
                and int(left.get("repeated_root_slot", 1)) != int(right.get("repeated_root_slot", 1))
            )
            if not valid_multiplicity:
                return True
    return False


def _remaining_simultaneous_spike(
    points: Sequence[dict[str, object]],
) -> bool:
    ordered = sorted(points, key=lambda item: float(item["beta_deg"]))
    for left, right in zip(ordered, ordered[1:]):
        left_values = np.asarray([float(entry["value"]) for entry in left["entries"]], dtype=float)
        right_values = np.asarray([float(entry["value"]) for entry in right["entries"]], dtype=float)
        if left_values.size < SORTED_PILOT_N_STORE or right_values.size < SORTED_PILOT_N_STORE:
            return True
        jump = np.abs(right_values - left_values)
        scale = np.maximum(np.maximum(np.abs(left_values), np.abs(right_values)), 1.0e-12)
        relative = jump / scale
        direction = np.sign(right_values - left_values)
        noticeable = np.isfinite(jump) & (
            (relative >= SORTED_PILOT_NOTICEABLE_JUMP_REL)
            | (jump >= SORTED_PILOT_NOTICEABLE_JUMP_ABS)
        )
        positive = int(np.count_nonzero(noticeable & (direction > 0.0)))
        negative = int(np.count_nonzero(noticeable & (direction < 0.0)))
        if max(positive, negative) >= SORTED_PILOT_SIMULTANEOUS_COUNT:
            return True
    return False


def _refinement_summary_rows(
    source_inventory: dict[tuple[str, int, str], dict[str, object]],
    region_results: dict[str, dict[str, object]],
    before_after: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for region in refined_pilot_regions():
        points = sorted(
            region_results[region.region_id]["points"],  # type: ignore[index]
            key=lambda item: float(item["beta_deg"]),
        )
        unresolved = sum(point["status"] == "local_refinement_unresolved" for point in points)
        settings = _local_refinement_settings(region)
        sorted_ok = all(_inventory_is_sorted(point["entries"]) for point in points)
        duplicate_ok = not any(_inventory_has_invalid_duplicate(point["entries"], settings) for point in points)
        endpoints_match = True
        for beta_deg in (region.beta_min, region.beta_max):
            point = next(item for item in points if abs(float(item["beta_deg"]) - beta_deg) <= 1.0e-10)
            original = source_inventory[(region.case_id, int(round(beta_deg)), region.model)]["entries"]
            for before_entry, after_entry in zip(original, point["entries"]):
                before_value = float(before_entry["value"])
                after_value = float(after_entry["value"])
                if not (isfinite(before_value) and isfinite(after_value)) or abs(before_value - after_value) > settings.root_match_tol:
                    endpoints_match = False
                    break
        remaining_spike = _remaining_simultaneous_spike(points)
        relevant_before_after = [
            row
            for row in before_after
            if row["case_id"] == region.case_id
            and row["theory"] == ("Euler-Bernoulli" if region.model == MODEL_EB else "Timoshenko")
            and region.beta_min <= float(row["beta_deg"]) <= region.beta_max
        ]
        finite_differences = [
            float(row["absolute_difference"])
            for row in relevant_before_after
            if isfinite(float(row["absolute_difference"]))
        ]
        if unresolved:
            status = "local_refinement_unresolved"
        elif not sorted_ok or not duplicate_ok:
            status = "integrity_failure"
        elif remaining_spike:
            status = "resolved_with_remaining_spike"
        elif not endpoints_match:
            status = "resolved_edge_changed"
        else:
            status = "resolved"
        rows.append(
            {
                "region_id": region.region_id,
                "beta_points": len(points),
                "resolved_points": len(points) - unresolved,
                "unresolved_points": unresolved,
                "recovered_root_count": sum(int(point["recovered_root_count"]) for point in points),
                "multiplicity_two_count": sum(int(point["multiplicity_two_count"]) for point in points),
                "max_before_after_difference": max(finite_differences, default=0.0),
                "remaining_spike": remaining_spike,
                "status": status,
            }
        )
    return rows


def _plot_refined_case(
    output_dir: Path,
    case: CaseSpec,
    rows: Sequence[dict[str, object]],
) -> Path:
    case_rows = [row for row in rows if row["case_id"] == case.label]
    colors = plt.get_cmap("tab10")(np.arange(SORTED_PILOT_N_PLOT))
    fig, ax = plt.subplots(figsize=(12.4, 7.6))
    for mode_zero in range(SORTED_PILOT_N_PLOT):
        mode_rows = [row for row in case_rows if int(row["sorted_mode_index"]) == mode_zero + 1]
        for model, linestyle in ((MODEL_EB, "-"), (MODEL_TIMO, "--")):
            prefix = _model_csv_prefix(model)
            finite_points = sorted(
                (
                    (float(row["beta_deg"]), float(row[f"lambda_{prefix}"]))
                    for row in mode_rows
                    if isfinite(float(row[f"lambda_{prefix}"]))
                ),
                key=lambda item: item[0],
            )
            ax.plot(
                [item[0] for item in finite_points],
                [item[1] for item in finite_points],
                color=colors[mode_zero],
                lw=1.35,
                ls=linestyle,
            )
    ax.set_xlabel(r"$\beta$, deg")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_xlim(0.0, 90.0)
    ax.set_xticks(np.arange(0.0, 91.0, 10.0))
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_title(
        rf"{case.label}: locally refined sorted spectra, $\epsilon_0={case.epsilon:.16g}$, "
        rf"$\mu={case.mu:g}$, $\eta={case.eta:g}$"
    )
    theory_handles = [
        Line2D([0], [0], color="black", lw=1.8, ls="-", label="Euler-Bernoulli"),
        Line2D([0], [0], color="black", lw=1.8, ls="--", label="Timoshenko"),
    ]
    mode_handles = [
        Line2D([0], [0], color=colors[index], lw=2.0, label=f"sorted mode {index + 1}")
        for index in range(SORTED_PILOT_N_PLOT)
    ]
    theory_legend = ax.legend(handles=theory_handles, loc="upper left", fontsize=9, frameon=True)
    ax.add_artist(theory_legend)
    ax.legend(
        handles=mode_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.13),
        ncol=5,
        fontsize=8,
        frameon=False,
    )
    fig.tight_layout(rect=(0.0, 0.07, 1.0, 1.0))
    path = output_dir / f"{case.label}_refined_sorted_lambda_beta_comparison.png"
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def validate_refined_pilot_outputs(
    refined_rows: Sequence[dict[str, object]],
    source_inventory: dict[tuple[str, int, str], dict[str, object]],
    region_results: dict[str, dict[str, object]],
    png_paths: Sequence[Path],
    output_dir: Path,
) -> dict[str, object]:
    point_index = _refinement_point_index(region_results)
    unaffected_mismatches = 0
    p2_mismatches = 0
    for case_id, _epsilon, _mu, _eta in SORTED_PILOT_CASES:
        for beta_deg in range(91):
            point_rows = [
                row
                for row in refined_rows
                if row["case_id"] == case_id and abs(float(row["beta_deg"]) - beta_deg) <= 1.0e-12
            ]
            for model in MODELS:
                targeted = (case_id, model, round(float(beta_deg), 10)) in point_index
                prefix = _model_csv_prefix(model)
                original = source_inventory[(case_id, beta_deg, model)]["entries"]
                for mode_zero, row in enumerate(point_rows):
                    value = float(row[f"lambda_{prefix}"])
                    original_value = float(original[mode_zero]["value"])
                    if not targeted and value != original_value:
                        unaffected_mismatches += 1
                    if case_id == "P2" and value != original_value:
                        p2_mismatches += 1
    if unaffected_mismatches:
        raise AssertionError(f"{unaffected_mismatches} unaffected integer values changed")
    if p2_mismatches:
        raise AssertionError(f"P2 changed in {p2_mismatches} values")
    if len(png_paths) != 4 or len(list(output_dir.glob("*.png"))) != 4:
        raise AssertionError("refined pilot must contain exactly four PNG files")
    forbidden = [p for p in output_dir.rglob("*") if p.suffix.lower() in {".pdf", ".svg", ".eps"}]
    if forbidden:
        raise AssertionError(f"forbidden output formats: {forbidden}")
    for region in refined_pilot_regions():
        settings = _local_refinement_settings(region)
        for point in region_results[region.region_id]["points"]:  # type: ignore[index]
            entries = point["entries"]
            if not _inventory_is_sorted(entries):
                raise AssertionError(f"{region.region_id} beta={point['beta_deg']}: inventory is not sorted")
            if _inventory_has_invalid_duplicate(entries, settings):
                raise AssertionError(f"{region.region_id} beta={point['beta_deg']}: invalid duplicate")
            if point["status"] != "local_refinement_unresolved" and (
                len(entries) < SORTED_PILOT_N_STORE
                or any(not isfinite(float(entry["value"])) for entry in entries)
            ):
                raise AssertionError(f"{region.region_id} beta={point['beta_deg']}: resolved inventory is incomplete")
    return {
        "unaffected_integer_mismatches": unaffected_mismatches,
        "p2_mismatches": p2_mismatches,
        "png_count": len(png_paths),
        "force_strict_calls": 0,
        "tracking_calls": 0,
        "mac_calls": 0,
        "shape_calls": 0,
        "continuation_calls": 0,
    }


def _accepted_target_local_slots(
    point: dict[str, object],
    region: LocalRefinementRegion,
) -> int:
    settings = _local_refinement_settings(region)
    candidates = sorted(
        (
            row
            for row in point["candidate_rows"]  # type: ignore[index]
            if bool(row["accepted"])
            and str(row["candidate_source"]).startswith(region.region_id + "_local_dense:")
            and str(row["block_family"]) == "full_6x6"
        ),
        key=lambda row: float(row["lambda_candidate"]),
    )
    groups: list[dict[str, object]] = []
    for candidate in candidates:
        if groups and abs(float(candidate["lambda_candidate"]) - float(groups[-1]["lambda_candidate"])) <= settings.root_dedup_tol:
            if int(candidate["multiplicity"]) > int(groups[-1]["multiplicity"]):
                groups[-1] = candidate
        else:
            groups.append(candidate)
    return sum(max(1, int(candidate["multiplicity"])) for candidate in groups)


def run_beta_sorted_spectrum_refined_pilot() -> dict[str, object]:
    source_dir = _repo_output_dir(REFINED_PILOT_SOURCE_DIR)
    output_dir = _repo_output_dir(REFINED_PILOT_OUTPUT_DIR)
    source_csv = source_dir / "beta_sorted_spectrum_pilot.csv"
    if not source_csv.exists():
        raise FileNotFoundError(f"source sorted pilot is missing: {source_csv}")
    source_csv_sha256 = hashlib.sha256(source_csv.read_bytes()).hexdigest()
    output_dir.mkdir(parents=True, exist_ok=True)
    region_results: dict[str, dict[str, object]] = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = {
            pool.submit(
                _refinement_region_worker,
                region,
                str(source_dir),
                str(output_dir),
                source_csv_sha256,
            ): region.region_id
            for region in refined_pilot_regions()
        }
        for future in as_completed(futures):
            result = future.result()
            region_id = str(result["region"]["region_id"])  # type: ignore[index]
            region_results[region_id] = result
            print(f"[{region_id}] local refinement complete", flush=True)

    source_inventory = _load_original_sorted_inventory(source_dir)
    refined_rows = _assemble_refined_rows(source_inventory, region_results)
    refined_csv = output_dir / "refined_beta_sorted_spectrum.csv"
    _write_csv(refined_csv, refined_rows, REFINED_PILOT_CSV_FIELDS)
    candidate_rows = _candidate_output_rows(region_results)
    candidates_csv = output_dir / "local_root_candidates.csv"
    _write_csv(candidates_csv, candidate_rows, LOCAL_CANDIDATE_FIELDS)
    before_after = _before_after_rows(source_inventory, region_results)
    before_after_csv = output_dir / "before_after_integer_beta.csv"
    _write_csv(before_after_csv, before_after, BEFORE_AFTER_FIELDS)
    summary_rows = _refinement_summary_rows(source_inventory, region_results, before_after)
    summary_csv = output_dir / "local_refinement_summary.csv"
    _write_csv(summary_csv, summary_rows, LOCAL_REFINEMENT_SUMMARY_FIELDS)

    png_paths: list[Path] = []
    for label, epsilon, mu, eta in SORTED_PILOT_CASES:
        png_paths.append(
            _plot_refined_case(
                output_dir,
                CaseSpec(mu=mu, eta=eta, epsilon=epsilon, label=label),
                refined_rows,
            )
        )
    validation = validate_refined_pilot_outputs(
        refined_rows,
        source_inventory,
        region_results,
        png_paths,
        output_dir,
    )

    local_matrix_evaluations = sum(
        int(point["local_matrix_evaluations"])
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    base_primary_evaluations = sum(
        int(point["base_primary_matrix_evaluations"])
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    sign_change_candidates = sum(bool(row["sign_change"]) for row in candidate_rows)
    sigma_candidates = sum(
        not bool(row["sign_change"]) and "sigma" in str(row["candidate_source"])
        for row in candidate_rows
    )
    rejected_candidates = sum(not bool(row["accepted"]) for row in candidate_rows)
    accepted_local_roots = sum(
        _accepted_target_local_slots(point, region)
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    distinct_pairs = sum(
        point["status"] == "resolved_distinct_pair"
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    double_roots = sum(
        int(point["multiplicity_two_count"])
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    unresolved_points = sum(
        point["status"] == "local_refinement_unresolved"
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    additional_betas = sum(
        not _is_integer_beta(float(point["beta_deg"]))
        for region in refined_pilot_regions()
        for point in region_results[region.region_id]["points"]  # type: ignore[index]
    )
    p3_beta0 = next(
        point
        for point in region_results["R3"]["points"]  # type: ignore[index]
        if abs(float(point["beta_deg"])) <= 1.0e-12
    )
    cache_hits = sum(int(region_results[region.region_id]["cache_hits"]) for region in refined_pilot_regions())
    report_lines = [
        "# Locally refined independent beta sorted-spectrum pilot",
        "",
        "Curves show sorted spectral positions, not physical descendant branches.",
        "",
        "Only R1-R7 were refined. Integer data outside those theory-specific windows come directly from the original sorted-spectrum pilot. Fractional-beta primary inventories and all accepted local candidates are matrix-confirmed at their own beta; no frequency interpolation is used.",
        "",
        "## Operation summary",
        "",
        f"- additional fractional beta-points: {additional_betas}",
        f"- local dense-window matrix evaluations: {local_matrix_evaluations}",
        f"- fractional base-primary matrix evaluations: {base_primary_evaluations}",
        f"- raw sign-change candidates: {sign_change_candidates}",
        f"- raw sigma-minimum-only candidates: {sigma_candidates}",
        f"- rejected candidates: {rejected_candidates}",
        f"- accepted deduplicated target-local root slots: {accepted_local_roots}",
        f"- resolved distinct-pair points: {distinct_pairs}",
        f"- confirmed double-root groups: {double_roots}",
        f"- unresolved local points: {unresolved_points}",
        f"- P3 beta=0 status: {p3_beta0['status']}",
        f"- cache hits in this invocation: {cache_hits}",
        "",
        "## Region summary",
        "",
    ]
    for row in summary_rows:
        report_lines.append(
            f"- {row['region_id']}: beta_points={row['beta_points']}, resolved={row['resolved_points']}, "
            f"unresolved={row['unresolved_points']}, recovered={row['recovered_root_count']}, "
            f"multiplicity_two={row['multiplicity_two_count']}, remaining_spike={row['remaining_spike']}, "
            f"status={row['status']}"
        )
    report_lines.extend(
        [
            "",
            "## Prohibited-operation counters",
            "",
            "- branch tracking calls: 0",
            "- MAC calls: 0",
            "- shape reconstruction calls: 0",
            "- continuation calls: 0",
            "- force/full strict calls: 0",
            "- global K12 strict scans: 0",
            "",
            "## Outputs",
            "",
        ]
    )
    report_lines.extend(f"- `{_rel(path)}`" for path in png_paths)
    report_lines.extend(
        [
            f"- `{_rel(refined_csv)}`",
            f"- `{_rel(candidates_csv)}`",
            f"- `{_rel(before_after_csv)}`",
            f"- `{_rel(summary_csv)}`",
            f"- `{_rel(output_dir / 'report.md')}`",
            "",
        ]
    )
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    return {
        "png_paths": png_paths,
        "refined_csv": refined_csv,
        "candidates_csv": candidates_csv,
        "before_after_csv": before_after_csv,
        "summary_csv": summary_csv,
        "report_path": report_path,
        "summary_rows": summary_rows,
        "additional_betas": additional_betas,
        "local_matrix_evaluations": local_matrix_evaluations,
        "base_primary_evaluations": base_primary_evaluations,
        "sign_change_candidates": sign_change_candidates,
        "sigma_candidates": sigma_candidates,
        "accepted_local_roots": accepted_local_roots,
        "rejected_candidates": rejected_candidates,
        "distinct_pairs": distinct_pairs,
        "double_roots": double_roots,
        "unresolved_points": unresolved_points,
        "p3_beta0_status": p3_beta0["status"],
        "cache_hits": cache_hits,
        "validation": validation,
    }


def run_beta_branch_pilot() -> dict[str, object]:
    output_dir = _repo_output_dir(PILOT_OUTPUT_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)
    results: dict[str, dict[str, object]] = {}
    with ProcessPoolExecutor(max_workers=4) as pool:
        futures = {pool.submit(_pilot_case_worker, case): case[0] for case in PILOT_CASES}
        for future in as_completed(futures):
            result = future.result()
            label = str(result["case"][0])  # type: ignore[index]
            results[label] = result
            print(f"[{label}] complete", flush=True)

    csv_rows: list[dict[str, object]] = []
    png_paths: list[Path] = []
    full_by_case: dict[str, list[int]] = {}
    gaps_by_case: dict[str, list[int]] = {}
    for label, epsilon, mu, eta in PILOT_CASES:
        case = CaseSpec(mu=mu, eta=eta, epsilon=epsilon, label=label)
        models = results[label]["models"]
        eb = np.asarray(models[MODEL_EB]["values"], dtype=float)  # type: ignore[index]
        tm = np.asarray(models[MODEL_TIMO]["values"], dtype=float)  # type: ignore[index]
        eb_status = np.asarray(models[MODEL_EB]["statuses"], dtype=object)  # type: ignore[index]
        tm_status = np.asarray(models[MODEL_TIMO]["statuses"], dtype=object)  # type: ignore[index]
        png_paths.append(_pilot_plot(output_dir, case, eb, tm))
        full: list[int] = []
        gaps: list[int] = []
        for branch_zero in range(PILOT_N_TRACK):
            if np.all(np.isfinite(eb[branch_zero])) and np.all(np.isfinite(tm[branch_zero])):
                full.append(branch_zero + 1)
            else:
                gaps.append(branch_zero + 1)
            for beta_index in range(91):
                status_parts = []
                if str(eb_status[branch_zero, beta_index]) != "ok":
                    status_parts.append(str(eb_status[branch_zero, beta_index]))
                if str(tm_status[branch_zero, beta_index]) != "ok":
                    status_parts.append(str(tm_status[branch_zero, beta_index]))
                csv_rows.append(
                    {
                        "case_id": label,
                        "epsilon_0": epsilon,
                        "mu": mu,
                        "eta": eta,
                        "beta_deg": beta_index,
                        "branch_index": branch_zero + 1,
                        "lambda_eb": eb[branch_zero, beta_index],
                        "lambda_timo": tm[branch_zero, beta_index],
                        "tracking_status": ";".join(status_parts) if status_parts else "ok",
                    }
                )
        full_by_case[label] = full
        gaps_by_case[label] = gaps
    csv_path = output_dir / "beta_branch_pilot.csv"
    _write_csv(csv_path, csv_rows, PILOT_CSV_FIELDS)
    gap_rows = [row for row in csv_rows if row["tracking_status"] != "ok"]
    report_lines = [
        "# Beta branch pilot",
        "",
        "Diagnostic-only EB/Timoshenko descendant-frequency comparison.",
        "Branches are seeded at beta=0, mu=0 for each epsilon and fixed eta, continued to the requested mu, then followed over beta=0..90 deg in one-degree steps.",
        "Internal shape continuity is used only for branch assignment; no shape, MAC, energy, gap, veering, FEM, strict, refinement, or article figure is emitted.",
        "",
        "- plotted branches per theory: 12",
        "- candidate roots per tracking step: 18",
        f"- ambiguous CSV rows saved as NaN gaps: {len(gap_rows)}",
        "- force_strict_verification calls: 0",
        "- adaptive refinement calls: 0",
        "",
        "## Coverage",
        "",
    ]
    for label, epsilon, mu, eta in PILOT_CASES:
        report_lines.append(
            f"- {label} (epsilon_0={epsilon:.17g}, mu={mu:g}, eta={eta:g}): "
            f"full branches={full_by_case[label] or 'none'}; branches with gaps={gaps_by_case[label] or 'none'}"
        )
    report_lines.extend(["", "## Outputs", ""])
    report_lines.extend(f"- `{_rel(path)}`" for path in png_paths)
    report_lines.extend([f"- `{_rel(csv_path)}`", f"- `{_rel(output_dir / 'report.md')}`", ""])
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    return {
        "png_paths": png_paths,
        "csv_path": csv_path,
        "report_path": report_path,
        "gap_row_count": len(gap_rows),
        "full_by_case": full_by_case,
        "gaps_by_case": gaps_by_case,
    }


def main(argv: list[str] | None = None) -> None:
    run_start = time.perf_counter()
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    args = parse_args(argv)
    if args.beta_branch_pilot:
        result = run_beta_branch_pilot()
        print(f"saved {len(result['png_paths'])} PNG files")
        print(f"CSV: {_rel(result['csv_path'])}")
        print(f"tracking gaps: {result['gap_row_count']}")
        print("force_strict=0; refinement=0; FEM=0; PDF=0")
        return
    if args.beta_sorted_spectrum_pilot:
        result = run_beta_sorted_spectrum_pilot()
        print(f"saved {len(result['png_paths'])} PNG files")
        print(f"CSV: {_rel(result['csv_path'])}")
        print(f"suspect intervals: {len(result['suspect_rows'])}")
        print(
            f"complete inventories EB/Timoshenko: {result['full_eb']}/364, "
            f"{result['full_timo']}/364"
        )
        print("tracking=0; MAC=0; shapes=0; force_strict=0; refinement=0; FEM=0; PDF=0")
        return
    if args.beta_sorted_spectrum_refined_pilot:
        result = run_beta_sorted_spectrum_refined_pilot()
        print(f"saved {len(result['png_paths'])} refined PNG files")
        print(f"CSV: {_rel(result['refined_csv'])}")
        print(
            f"additional_beta={result['additional_betas']}; local_matrix_evaluations="
            f"{result['local_matrix_evaluations']}; unresolved={result['unresolved_points']}"
        )
        print(
            f"distinct_pairs={result['distinct_pairs']}; double_roots={result['double_roots']}; "
            f"P3_beta0={result['p3_beta0_status']}"
        )
        print("tracking=0; MAC=0; shapes=0; continuation=0; force_strict=0; FEM=0; PDF=0")
        return
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
