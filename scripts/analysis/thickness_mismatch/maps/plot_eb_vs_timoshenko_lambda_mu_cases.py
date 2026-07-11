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

try:  # Package import for tests; direct import for script execution.
    from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
        plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow,
    )
except ImportError:  # pragma: no cover - exercised when run as a file.
    import plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow  # type: ignore[no-redef]  # noqa: E402

from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


DEFAULT_OUTPUT_STEM = "eb_vs_timo_lambda_mu_beta45_eta0_eps_scan"
GENERIC_OUTPUT_STEM = "eb_vs_timo_lambda_mu_cases"
DEFAULT_OUTPUT_DIR = Path("results") / "eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_vs_timoshenko_lambda_mu_beta45_eta0_eps_scan"
DEFAULT_CACHE_SUBDIR = "cache"

DEFAULT_BETA_VALUES = (45.0,)
DEFAULT_ETA = 0.0
DEFAULT_EPSILON_VALUES = (0.005, 0.01, 0.02, 0.03, 0.04)
SMOKE_EPSILON_VALUES = (0.005, 0.04)
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.70
DEFAULT_MU_STEP = 0.005
DEFAULT_LOCAL_MU_STEP = 0.001
DEFAULT_REPAIR_WINDOW_MU = 0.02
DEFAULT_N_ROOTS = 6
SMOKE_N_ROOTS = 4

DEFAULT_JUMP_ABS_THRESHOLD = 0.75
DEFAULT_JUMP_REL_THRESHOLD = 0.08
MAX_RECOVERY_INTERVAL_MU = 0.08
CACHE_VERSION = "eb_vs_timo_lambda_mu_cases_v2"
TIMING_REPORT_NAME = "timing_report.csv"

MODEL_EB = beta_workflow.MODEL_EB
MODEL_TIMO = beta_workflow.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)


CASE_CSV_FIELDS = [
    "mu",
    "beta_deg",
    "eta",
    "epsilon",
    "sorted_index",
    "Lambda_EB_raw",
    "Lambda_Timoshenko_raw",
    "Lambda_EB_plot",
    "Lambda_Timoshenko_plot",
    "rel_diff_Timo_vs_EB",
    "suspicious_EB",
    "suspicious_Timoshenko",
    "retry_attempted",
    "retry_fixed",
    "plotted_as_nan",
    "root_solver_warning",
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
    "beta_deg",
    "eta",
    "epsilon",
    "model",
    "sorted_index",
    "mu",
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
    "beta_deg",
    "eta",
    "epsilon",
    "sorted_index",
    "max_rel_diff_over_mu",
    "mean_rel_diff_over_mu",
    "mu_at_max_rel_diff",
    "suspicious_point_count",
    "retry_fixed_count",
    "nan_count_after_cleanup",
    "notes",
    "initial_raw_suspicious_point_count",
    "raw_suspicious_point_count",
]

TIMING_FIELDS = [
    "run_id",
    "model",
    "beta_deg",
    "eta",
    "epsilon",
    "n_mu_points",
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
    beta_deg: float
    eta: float
    epsilon: float


@dataclass(frozen=True)
class RefinementWindow:
    case_index: int
    mu_min: float
    mu_max: float
    reason: str


@dataclass(frozen=True)
class Args:
    mu_min: float
    mu_max: float
    mu_step: float
    local_mu_step: float
    repair_window_mu: float
    jump_abs_threshold: float
    jump_rel_threshold: float
    n_roots: int
    output_dir: Path
    cache_dir: Path
    reuse_cache: bool
    force_recompute: bool
    plot_only: bool
    repair_spikes: bool
    sv_recovery_only_on_spikes: bool
    timo_root_mode: str
    beta_values: tuple[float, ...]
    eta: float
    epsilon_values: tuple[float, ...]
    smoke: bool


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
    root_warning: dict[tuple[int, str], np.ndarray]
    root_count: dict[tuple[int, str], np.ndarray]


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


def output_stem(cases: Sequence[CaseSpec]) -> str:
    betas = sorted({round(float(case.beta_deg), 12) for case in cases})
    etas = sorted({round(float(case.eta), 12) for case in cases})
    if len(betas) == 1 and len(etas) == 1:
        return (
            "eb_vs_timo_lambda_mu_"
            f"beta{_float_label(betas[0])}_eta{_float_label(etas[0])}_eps_scan"
        )
    return GENERIC_OUTPUT_STEM


def summary_csv_path(output_dir: Path, cases: Sequence[CaseSpec]) -> Path:
    return output_dir / f"{output_stem(cases)}_summary.csv"


def spike_audit_csv_path(output_dir: Path, cases: Sequence[CaseSpec]) -> Path:
    return output_dir / f"{output_stem(cases)}_spike_audit.csv"


def overview_png_path(output_dir: Path, cases: Sequence[CaseSpec]) -> Path:
    return output_dir / f"{output_stem(cases)}_overview.png"


def report_md_path(output_dir: Path, cases: Sequence[CaseSpec]) -> Path:
    return output_dir / f"{output_stem(cases)}_report.md"


def regular_grid(start: float, end: float, step: float) -> np.ndarray:
    return beta_workflow.regular_grid(float(start), float(end), float(step))


def default_cases(
    beta_values: Sequence[float] = DEFAULT_BETA_VALUES,
    eta: float = DEFAULT_ETA,
    epsilon_values: Sequence[float] = DEFAULT_EPSILON_VALUES,
) -> tuple[CaseSpec, ...]:
    return tuple(
        CaseSpec(beta_deg=float(beta_deg), eta=float(eta), epsilon=float(epsilon))
        for epsilon in epsilon_values
        for beta_deg in beta_values
    )


def parse_args(argv: list[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic-only sorted in-plane Lambda(mu) comparison between "
            "Euler-Bernoulli and Timoshenko theories."
        )
    )
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--local-mu-step", type=float, default=DEFAULT_LOCAL_MU_STEP)
    parser.add_argument("--repair-window-mu", type=float, default=DEFAULT_REPAIR_WINDOW_MU)
    parser.add_argument("--jump-abs-threshold", type=float, default=DEFAULT_JUMP_ABS_THRESHOLD)
    parser.add_argument("--jump-rel-threshold", type=float, default=DEFAULT_JUMP_REL_THRESHOLD)
    parser.add_argument("--n-roots", type=int, default=DEFAULT_N_ROOTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--repair-spikes", dest="repair_spikes", action="store_true", default=True)
    parser.add_argument("--no-repair-spikes", dest="repair_spikes", action="store_false")
    parser.add_argument("--sv-recovery-only-on-spikes", dest="sv_recovery_only_on_spikes", action="store_true", default=True)
    parser.add_argument("--no-sv-recovery-only-on-spikes", dest="sv_recovery_only_on_spikes", action="store_false")
    parser.add_argument("--timo-root-mode", choices=("global", "continuation"), default="continuation")
    parser.add_argument("--beta-deg", type=float, default=None)
    parser.add_argument("--beta-values", type=float, nargs="+", default=None)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--epsilon-values", type=float, nargs="+", default=None)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)

    mu_min = float(ns.mu_min)
    mu_max = float(ns.mu_max)
    mu_step = float(ns.mu_step)
    local_mu_step = float(ns.local_mu_step)
    repair_window_mu = float(ns.repair_window_mu)
    n_roots = int(ns.n_roots)
    if ns.beta_values is not None and ns.beta_deg is not None:
        raise ValueError("use either --beta-deg or --beta-values, not both.")
    if ns.beta_values is not None:
        beta_values = tuple(float(value) for value in ns.beta_values)
    elif ns.beta_deg is not None:
        beta_values = (float(ns.beta_deg),)
    else:
        beta_values = tuple(DEFAULT_BETA_VALUES)
    epsilon_values = (
        tuple(float(value) for value in ns.epsilon_values)
        if ns.epsilon_values is not None
        else tuple(DEFAULT_EPSILON_VALUES)
    )
    output_dir = _repo_output_dir(Path(ns.output_dir))
    cache_dir = _repo_output_dir(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / DEFAULT_CACHE_SUBDIR
    repair_spikes = bool(ns.repair_spikes)
    if ns.smoke:
        mu_min = 0.0
        mu_max = 0.70
        mu_step = 0.35
        local_mu_step = 0.07
        repair_window_mu = 0.02
        n_roots = SMOKE_N_ROOTS
        epsilon_values = tuple(SMOKE_EPSILON_VALUES)
        output_dir = _repo_output_dir(SMOKE_OUTPUT_DIR)
        cache_dir = output_dir / DEFAULT_CACHE_SUBDIR
        repair_spikes = False

    args = Args(
        mu_min=mu_min,
        mu_max=mu_max,
        mu_step=mu_step,
        local_mu_step=local_mu_step,
        repair_window_mu=repair_window_mu,
        jump_abs_threshold=float(ns.jump_abs_threshold),
        jump_rel_threshold=float(ns.jump_rel_threshold),
        n_roots=n_roots,
        output_dir=output_dir,
        cache_dir=cache_dir,
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        plot_only=bool(ns.plot_only),
        repair_spikes=repair_spikes,
        sv_recovery_only_on_spikes=bool(ns.sv_recovery_only_on_spikes),
        timo_root_mode=str(ns.timo_root_mode),
        beta_values=beta_values,
        eta=float(ns.eta),
        epsilon_values=epsilon_values,
        smoke=bool(ns.smoke),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if not (isfinite(args.mu_min) and isfinite(args.mu_max) and isfinite(args.mu_step)):
        raise ValueError("mu grid values must be finite.")
    if args.mu_step <= 0.0:
        raise ValueError("mu-step must be positive.")
    if args.local_mu_step <= 0.0:
        raise ValueError("local-mu-step must be positive.")
    if args.local_mu_step > args.mu_step / 5.0 + 1.0e-12:
        raise ValueError("local-mu-step must be at least 5 times smaller than mu-step.")
    if args.repair_window_mu <= 0.0:
        raise ValueError("repair-window-mu must be positive.")
    if args.mu_max < args.mu_min:
        raise ValueError("mu-max must be greater than or equal to mu-min.")
    if args.jump_abs_threshold <= 0.0 or args.jump_rel_threshold <= 0.0:
        raise ValueError("jump thresholds must be positive.")
    if args.n_roots <= 0:
        raise ValueError("n-roots must be positive.")
    if not args.beta_values:
        raise ValueError("at least one beta value is required.")
    if not args.epsilon_values:
        raise ValueError("at least one epsilon value is required.")
    for beta_deg in args.beta_values:
        if not isfinite(beta_deg):
            raise ValueError("beta values must be finite.")
    for epsilon in args.epsilon_values:
        if epsilon <= 0.0:
            raise ValueError("epsilon values must be positive.")
    for mu in (args.mu_min, args.mu_max):
        thickness_mismatch_factors(mu, args.eta)
        TIMO.tau_factors(mu, args.eta)


def base_mu_grid(args: Args) -> np.ndarray:
    return regular_grid(args.mu_min, args.mu_max, args.mu_step)


def cache_settings(cases: Sequence[CaseSpec], args: Args) -> dict[str, object]:
    return {
        "cache_version": CACHE_VERSION,
        "cases": [
            {"beta_deg": case.beta_deg, "eta": case.eta, "epsilon": case.epsilon}
            for case in cases
        ],
        "mu_min": args.mu_min,
        "mu_max": args.mu_max,
        "mu_step": args.mu_step,
        "local_mu_step": args.local_mu_step,
        "repair_window_mu": args.repair_window_mu,
        "jump_abs_threshold": args.jump_abs_threshold,
        "jump_rel_threshold": args.jump_rel_threshold,
        "n_roots": args.n_roots,
        "repair_spikes": args.repair_spikes,
        "sv_recovery_only_on_spikes": args.sv_recovery_only_on_spikes,
        "timo_root_mode": args.timo_root_mode,
        "beta_values": [float(value) for value in args.beta_values],
        "eta": args.eta,
        "epsilon_values": [float(value) for value in args.epsilon_values],
        "ordinary_pass": "base mu grid followed by spike-triggered local mu refinement",
        "timo_continuation": "determinant solves in local brackets, sorted roots saved at each mu",
    }


def cache_settings_json(cases: Sequence[CaseSpec], args: Args) -> str:
    return json.dumps(cache_settings(cases, args), sort_keys=True, separators=(",", ":"))


def cache_file_path(cases: Sequence[CaseSpec], args: Args) -> Path:
    digest = hashlib.sha256(cache_settings_json(cases, args).encode("utf-8")).hexdigest()[:16]
    return args.cache_dir / f"roots_all_cases_{CACHE_VERSION}_{digest}.npz"


def root_key(case_index: int, model: str, mu: float) -> tuple[int, str, float]:
    return case_index, model, round(float(mu), 12)


def solve_point(
    case: CaseSpec,
    mu: float,
    args: Args,
    model: str,
    *,
    previous_roots: Sequence[float] | None = None,
    retry: bool = False,
    upper_hint: float | None = None,
    sv_recovery: bool = False,
    sv_target_indices: Sequence[int] | None = None,
) -> beta_workflow.RootResult:
    beta_case = beta_workflow.CaseSpec(mu=float(mu), eta=case.eta, epsilon=case.epsilon)
    if model == MODEL_TIMO and args.timo_root_mode == "continuation" and not retry and previous_roots is not None:
        return beta_workflow.solve_timo_continuation(
            beta_case,
            case.beta_deg,
            args.n_roots,
            previous_roots,
            upper_hint=upper_hint,
        )
    return beta_workflow.solve_model(
        beta_case,
        case.beta_deg,
        args.n_roots,
        model,
        retry=retry,
        upper_hint=upper_hint,
        sv_recovery=sv_recovery,
        sv_target_indices=sv_target_indices,
    )


def compute_roots_for_grids(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    existing: dict[tuple[int, str, float], beta_workflow.RootResult] | None = None,
) -> dict[tuple[int, str, float], beta_workflow.RootResult]:
    root_map = dict(existing or {})
    for case_index, case in enumerate(cases):
        mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
        print(
            "solving Lambda(mu) case "
            f"{case_index + 1}/{len(cases)}: beta={case.beta_deg:g} deg, eta={case.eta:g}, "
            f"epsilon={case.epsilon:g}, mu points={len(mu_values)}"
        )
        previous_timo_roots: tuple[float, ...] | None = None
        for mu in mu_values:
            eb_key = root_key(case_index, MODEL_EB, float(mu))
            if eb_key not in root_map:
                eb_result = solve_point(case, float(mu), args, MODEL_EB)
                if eb_result.root_count_found < args.n_roots:
                    eb_result = solve_point(case, float(mu), args, MODEL_EB, retry=True)
                root_map[eb_key] = eb_result

            timo_key = root_key(case_index, MODEL_TIMO, float(mu))
            if timo_key not in root_map:
                eb_roots = root_map[eb_key].roots
                finite_eb = [float(root) for root in eb_roots if isfinite(float(root))]
                upper_hint = max(finite_eb) if finite_eb else None
                timo_result = solve_point(
                    case,
                    float(mu),
                    args,
                    MODEL_TIMO,
                    previous_roots=previous_timo_roots,
                    upper_hint=upper_hint,
                )
                if timo_result.root_count_found < args.n_roots:
                    timo_result = solve_point(
                        case,
                        float(mu),
                        args,
                        MODEL_TIMO,
                        retry=True,
                        upper_hint=upper_hint,
                    )
                root_map[timo_key] = timo_result
            previous_timo_roots = root_map[timo_key].roots
    return root_map


def jump_abs_rel(left: float, right: float) -> tuple[float, float]:
    return beta_workflow.jump_abs_rel(float(left), float(right))


def row_jump_metrics(values: np.ndarray, sorted_zero: int, mu_index: int) -> tuple[float, float, float, float]:
    row = np.asarray(values[sorted_zero, :], dtype=float)
    current = float(row[mu_index])
    if mu_index > 0:
        jump_prev, rel_prev = jump_abs_rel(float(row[mu_index - 1]), current)
    else:
        jump_prev, rel_prev = float("nan"), float("nan")
    if mu_index + 1 < row.size:
        jump_next, rel_next = jump_abs_rel(current, float(row[mu_index + 1]))
    else:
        jump_next, rel_next = float("nan"), float("nan")
    return jump_prev, jump_next, rel_prev, rel_next


def jump_is_suspicious(jump_abs: float, jump_rel: float, args: Args) -> bool:
    return (
        isfinite(float(jump_abs))
        and isfinite(float(jump_rel))
        and (float(jump_abs) > args.jump_abs_threshold or float(jump_rel) > args.jump_rel_threshold)
    )


def point_is_suspicious(values: np.ndarray, sorted_zero: int, mu_index: int, args: Args) -> bool:
    jump_prev, jump_next, rel_prev, rel_next = row_jump_metrics(values, sorted_zero, mu_index)
    return jump_is_suspicious(jump_prev, rel_prev, args) or jump_is_suspicious(jump_next, rel_next, args)


def suspicious_mask(values: np.ndarray, args: Args) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    out = np.full(array.shape, False, dtype=bool)
    for sorted_zero in range(array.shape[0]):
        for mu_index in range(array.shape[1]):
            if isfinite(float(array[sorted_zero, mu_index])):
                out[sorted_zero, mu_index] = point_is_suspicious(array, sorted_zero, mu_index, args)
    return out


def recovery_mask(values: np.ndarray, mu_values: np.ndarray, args: Args) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    mu_array = np.asarray(mu_values, dtype=float)
    if not args.repair_spikes:
        return np.full(array.shape, False, dtype=bool)
    out = suspicious_mask(array, args)
    for sorted_zero in range(array.shape[0]):
        bad_edges: list[int] = []
        for mu_index in range(array.shape[1] - 1):
            jump_abs, _jump_rel = jump_abs_rel(array[sorted_zero, mu_index], array[sorted_zero, mu_index + 1])
            if isfinite(jump_abs) and jump_abs > args.jump_abs_threshold:
                bad_edges.append(mu_index)
        for left_edge, right_edge in zip(bad_edges, bad_edges[1:]):
            interval_width = abs(float(mu_array[right_edge]) - float(mu_array[left_edge]))
            if interval_width <= MAX_RECOVERY_INTERVAL_MU + 1.0e-12:
                out[sorted_zero, left_edge + 1 : right_edge + 1] = True
    return out


def raw_arrays_from_root_map(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
) -> tuple[
    dict[tuple[int, str], np.ndarray],
    dict[tuple[int, str], np.ndarray],
    dict[tuple[int, str], np.ndarray],
    dict[tuple[int, str], np.ndarray],
]:
    raw: dict[tuple[int, str], np.ndarray] = {}
    root_warning: dict[tuple[int, str], np.ndarray] = {}
    root_count: dict[tuple[int, str], np.ndarray] = {}
    solve_notes: dict[tuple[int, str], np.ndarray] = {}
    for case_index, _case in enumerate(cases):
        mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            raw_array = np.full((args.n_roots, len(mu_values)), np.nan, dtype=float)
            warning_array = np.full((args.n_roots, len(mu_values)), "", dtype=object)
            count_array = np.full((args.n_roots, len(mu_values)), 0, dtype=int)
            notes_array = np.full((args.n_roots, len(mu_values)), "", dtype=object)
            for mu_index, mu in enumerate(mu_values):
                result = root_map[root_key(case_index, model, float(mu))]
                for sorted_zero in range(args.n_roots):
                    raw_array[sorted_zero, mu_index] = float(result.roots[sorted_zero])
                    warning_array[sorted_zero, mu_index] = beta_workflow.warning_text(result)
                    count_array[sorted_zero, mu_index] = int(result.root_count_found)
                    notes_array[sorted_zero, mu_index] = "; ".join(result.notes)
            raw[key] = raw_array
            root_warning[key] = warning_array
            root_count[key] = count_array
            solve_notes[key] = notes_array
    return raw, root_warning, root_count, solve_notes


def refinement_windows_from_raw(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    raw: dict[tuple[int, str], np.ndarray],
) -> tuple[RefinementWindow, ...]:
    if not args.repair_spikes:
        return ()
    windows: list[RefinementWindow] = []
    seen: set[tuple[int, float, float]] = set()
    for case_index, _case in enumerate(cases):
        mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
        for model in MODELS:
            mask = suspicious_mask(raw[(case_index, model)], args)
            for _sorted_zero, mu_index in zip(*np.where(mask)):
                center = float(mu_values[int(mu_index)])
                left = max(args.mu_min, center - args.repair_window_mu)
                right = min(args.mu_max, center + args.repair_window_mu)
                key = (case_index, round(left, 12), round(right, 12))
                if key in seen:
                    continue
                seen.add(key)
                windows.append(
                    RefinementWindow(
                        case_index=case_index,
                        mu_min=left,
                        mu_max=right,
                        reason=f"automatic {model} jump near mu={center:.8g}",
                    )
                )
    return tuple(windows)


def refine_mu_grids(
    cases: Sequence[CaseSpec],
    base_mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    windows: Sequence[RefinementWindow],
) -> list[np.ndarray]:
    grids = [list(np.asarray(values, dtype=float)) for values in base_mu_values_by_case]
    for window in windows:
        if window.case_index < 0 or window.case_index >= len(cases):
            continue
        grids[window.case_index].extend(
            float(value) for value in regular_grid(window.mu_min, window.mu_max, args.local_mu_step)
        )
    return [np.unique(np.round(np.asarray(values, dtype=float), 12)) for values in grids]


def make_root_tables(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
    *,
    initial_raw_suspicious_count: int,
) -> tuple[RootTables, dict[str, float], int]:
    audit_start = time.perf_counter()
    raw, root_warning, root_count, solve_notes = raw_arrays_from_root_map(cases, mu_values_by_case, args, root_map)
    clean: dict[tuple[int, str], np.ndarray] = {}
    status: dict[tuple[int, str], np.ndarray] = {}
    notes: dict[tuple[int, str], np.ndarray] = {}
    retry_attempted: dict[tuple[int, str], np.ndarray] = {}
    retry_fixed: dict[tuple[int, str], np.ndarray] = {}
    plotted_as_nan: dict[tuple[int, str], np.ndarray] = {}
    raw_suspicious: dict[tuple[int, str], np.ndarray] = {}
    plot_suspicious: dict[tuple[int, str], np.ndarray] = {}
    recovery_suspicious: dict[tuple[int, str], np.ndarray] = {}

    for case_index, _case in enumerate(cases):
        mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            raw_array = raw[key]
            clean[key] = raw_array.copy()
            status_array = np.full(raw_array.shape, "ok", dtype=object)
            notes_array = np.array(solve_notes[key], dtype=object, copy=True)
            retry_attempted[key] = np.full(raw_array.shape, False, dtype=bool)
            retry_fixed[key] = np.full(raw_array.shape, False, dtype=bool)
            plotted_as_nan[key] = np.full(raw_array.shape, False, dtype=bool)
            for sorted_zero in range(args.n_roots):
                for mu_index in range(len(mu_values)):
                    if not isfinite(float(raw_array[sorted_zero, mu_index])):
                        status_array[sorted_zero, mu_index] = "missing_root_nan"
            status[key] = status_array
            notes[key] = notes_array
            raw_suspicious[key] = suspicious_mask(raw_array, args)
            recovery_suspicious[key] = recovery_mask(raw_array, mu_values, args)

    audit_seconds = time.perf_counter() - audit_start
    repair_start = time.perf_counter()
    sv_recovery_calls = 0

    for case_index, case in enumerate(cases):
        mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            values = clean[key]
            retry_cache: dict[tuple[int, str, int], beta_workflow.RootResult] = {}
            suspicious_mu_indices = sorted(
                {
                    mu_index
                    for sorted_zero in range(args.n_roots)
                    for mu_index in range(len(mu_values))
                    if bool(recovery_suspicious[key][sorted_zero, mu_index])
                }
            )
            for mu_index in suspicious_mu_indices:
                mu = float(mu_values[mu_index])
                upper_hint = None
                if model == MODEL_TIMO:
                    eb_values = clean[(case_index, MODEL_EB)][:, mu_index]
                    finite_eb = [float(root) for root in eb_values if isfinite(float(root))]
                    upper_hint = max(finite_eb) if finite_eb else None
                retry_key = (case_index, model, mu_index)
                if retry_key not in retry_cache:
                    retry_cache[retry_key] = solve_point(
                        case,
                        mu,
                        args,
                        model,
                        retry=True,
                        upper_hint=upper_hint,
                        sv_recovery=False,
                    )
                retry = retry_cache[retry_key]
                for sorted_zero in range(args.n_roots):
                    if not bool(recovery_suspicious[key][sorted_zero, mu_index]):
                        continue
                    retry_attempted[key][sorted_zero, mu_index] = True
                    retry_value = float(retry.roots[sorted_zero])
                    current = float(values[sorted_zero, mu_index])
                    if isfinite(retry_value) and abs(retry_value - current) > 1.0e-8:
                        values[sorted_zero, mu_index] = retry_value
                        retry_fixed[key][sorted_zero, mu_index] = True
                        status[key][sorted_zero, mu_index] = (
                            "retry_fixed_spike"
                            if bool(raw_suspicious[key][sorted_zero, mu_index])
                            else "retry_fixed_spike_interval"
                        )
                        notes[key][sorted_zero, mu_index] = (
                            "suspicious jump fixed by dense retry"
                            if bool(raw_suspicious[key][sorted_zero, mu_index])
                            else "dense retry fixed point inside suspicious jump interval"
                        )

            sv_retry_cache: dict[tuple[int, str, int, tuple[int, ...]], beta_workflow.RootResult] = {}
            for _sv_pass in range(3):
                plot_mask_after_dense = suspicious_mask(values, args)
                sv_mu_indices = sorted(
                    {
                        mu_index
                        for sorted_zero in range(args.n_roots)
                        for mu_index in range(len(mu_values))
                        if bool(plot_mask_after_dense[sorted_zero, mu_index])
                        and bool(args.repair_spikes)
                        and bool(args.sv_recovery_only_on_spikes)
                    }
                )
                if not sv_mu_indices:
                    break
                changed = False
                for mu_index in sv_mu_indices:
                    target_indices = tuple(
                        sorted_zero
                        for sorted_zero in range(args.n_roots)
                        if bool(plot_mask_after_dense[sorted_zero, mu_index])
                    )
                    mu = float(mu_values[mu_index])
                    upper_hint = None
                    if model == MODEL_TIMO:
                        eb_values = clean[(case_index, MODEL_EB)][:, mu_index]
                        finite_eb = [float(root) for root in eb_values if isfinite(float(root))]
                        upper_hint = max(finite_eb) if finite_eb else None
                    retry_key = (case_index, model, mu_index, target_indices)
                    if retry_key not in sv_retry_cache:
                        sv_recovery_calls += 1
                        sv_retry_cache[retry_key] = solve_point(
                            case,
                            mu,
                            args,
                            model,
                            retry=True,
                            upper_hint=upper_hint,
                            sv_recovery=True,
                            sv_target_indices=target_indices,
                        )
                    retry = sv_retry_cache[retry_key]
                    for sorted_zero in target_indices:
                        retry_attempted[key][sorted_zero, mu_index] = True
                        retry_value = float(retry.roots[sorted_zero])
                        current = float(values[sorted_zero, mu_index])
                        if isfinite(retry_value) and abs(retry_value - current) > 1.0e-8:
                            values[sorted_zero, mu_index] = retry_value
                            retry_fixed[key][sorted_zero, mu_index] = True
                            changed = True
                            status[key][sorted_zero, mu_index] = (
                                "retry_fixed_spike_svd"
                                if bool(raw_suspicious[key][sorted_zero, mu_index])
                                else "retry_fixed_spike_interval_svd"
                            )
                            notes[key][sorted_zero, mu_index] = "suspicious jump fixed by local SVD recovery"
                if not changed:
                    break

            plot_mask_before_nan = suspicious_mask(values, args)
            for sorted_zero in range(args.n_roots):
                for mu_index in range(len(mu_values)):
                    if not bool(plot_mask_before_nan[sorted_zero, mu_index]):
                        continue
                    if bool(args.repair_spikes) and not bool(retry_fixed[key][sorted_zero, mu_index]):
                        values[sorted_zero, mu_index] = float("nan")
                        plotted_as_nan[key][sorted_zero, mu_index] = True
                        status[key][sorted_zero, mu_index] = "unresolved_spike_nan"
                        notes[key][sorted_zero, mu_index] = "unresolved suspicious jump hidden from plot"
            plot_suspicious[key] = suspicious_mask(values, args)

    repair_seconds = time.perf_counter() - repair_start
    tables = RootTables(
        raw=raw,
        clean=clean,
        status=status,
        notes=notes,
        retry_attempted=retry_attempted,
        retry_fixed=retry_fixed,
        plotted_as_nan=plotted_as_nan,
        raw_suspicious=raw_suspicious,
        plot_suspicious=plot_suspicious,
        root_warning=root_warning,
        root_count=root_count,
    )
    timing_parts = {"spike_audit_seconds": audit_seconds, "repair_seconds": repair_seconds}
    return tables, timing_parts, int(sv_recovery_calls)


def root_map_arrays(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
) -> dict[str, np.ndarray]:
    n_cases = len(cases)
    max_mu = max(len(values) for values in mu_values_by_case) if mu_values_by_case else 0
    n_models = len(MODELS)
    roots = np.full((n_cases, n_models, max_mu, args.n_roots), np.nan, dtype=float)
    root_count = np.zeros((n_cases, n_models, max_mu), dtype=int)
    warnings = np.full((n_cases, n_models, max_mu), "", dtype=object)
    notes = np.full((n_cases, n_models, max_mu), "", dtype=object)
    lambda_max = np.full((n_cases, n_models, max_mu), np.nan, dtype=float)
    scan_step = np.full((n_cases, n_models, max_mu), np.nan, dtype=float)
    retry_attempted = np.full((n_cases, n_models, max_mu), False, dtype=bool)
    retry_changed = np.full((n_cases, n_models, max_mu), False, dtype=bool)
    for case_index, _case in enumerate(cases):
        for mu_index, mu in enumerate(np.asarray(mu_values_by_case[case_index], dtype=float)):
            for model_index, model in enumerate(MODELS):
                result = root_map[root_key(case_index, model, float(mu))]
                roots[case_index, model_index, mu_index, :] = np.asarray(result.roots, dtype=float)
                root_count[case_index, model_index, mu_index] = int(result.root_count_found)
                warnings[case_index, model_index, mu_index] = " | ".join(result.warnings)
                notes[case_index, model_index, mu_index] = "; ".join(result.notes)
                lambda_max[case_index, model_index, mu_index] = float(result.lambda_max_used)
                scan_step[case_index, model_index, mu_index] = float(result.scan_step_used)
                retry_attempted[case_index, model_index, mu_index] = bool(result.retry_attempted)
                retry_changed[case_index, model_index, mu_index] = bool(result.retry_changed_value)
    return {
        "root_roots": roots,
        "root_count": root_count,
        "root_warnings": warnings,
        "root_notes": notes,
        "root_lambda_max": lambda_max,
        "root_scan_step": scan_step,
        "root_retry_attempted": retry_attempted,
        "root_retry_changed": retry_changed,
    }


def tables_arrays(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    tables: RootTables,
) -> dict[str, np.ndarray]:
    n_cases = len(cases)
    max_mu = max(len(values) for values in mu_values_by_case) if mu_values_by_case else 0
    n_models = len(MODELS)
    shape = (n_cases, n_models, args.n_roots, max_mu)
    arrays: dict[str, np.ndarray] = {
        "table_raw": np.full(shape, np.nan, dtype=float),
        "table_clean": np.full(shape, np.nan, dtype=float),
        "table_status": np.full(shape, "", dtype=object),
        "table_notes": np.full(shape, "", dtype=object),
        "table_retry_attempted": np.full(shape, False, dtype=bool),
        "table_retry_fixed": np.full(shape, False, dtype=bool),
        "table_plotted_as_nan": np.full(shape, False, dtype=bool),
        "table_raw_suspicious": np.full(shape, False, dtype=bool),
        "table_plot_suspicious": np.full(shape, False, dtype=bool),
        "table_root_warning": np.full(shape, "", dtype=object),
        "table_root_count": np.zeros(shape, dtype=int),
    }
    for case_index, _case in enumerate(cases):
        mu_count = len(mu_values_by_case[case_index])
        for model_index, model in enumerate(MODELS):
            key = (case_index, model)
            arrays["table_raw"][case_index, model_index, :, :mu_count] = tables.raw[key]
            arrays["table_clean"][case_index, model_index, :, :mu_count] = tables.clean[key]
            arrays["table_status"][case_index, model_index, :, :mu_count] = tables.status[key]
            arrays["table_notes"][case_index, model_index, :, :mu_count] = tables.notes[key]
            arrays["table_retry_attempted"][case_index, model_index, :, :mu_count] = tables.retry_attempted[key]
            arrays["table_retry_fixed"][case_index, model_index, :, :mu_count] = tables.retry_fixed[key]
            arrays["table_plotted_as_nan"][case_index, model_index, :, :mu_count] = tables.plotted_as_nan[key]
            arrays["table_raw_suspicious"][case_index, model_index, :, :mu_count] = tables.raw_suspicious[key]
            arrays["table_plot_suspicious"][case_index, model_index, :, :mu_count] = tables.plot_suspicious[key]
            arrays["table_root_warning"][case_index, model_index, :, :mu_count] = tables.root_warning[key]
            arrays["table_root_count"][case_index, model_index, :, :mu_count] = tables.root_count[key]
    return arrays


def save_cache(
    path: Path,
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
    tables: RootTables,
    metadata: dict[str, object],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {
        "settings_json": np.array(cache_settings_json(cases, args), dtype=object),
        "metadata_json": np.array(json.dumps(metadata, sort_keys=True), dtype=object),
        "case_beta_deg": np.asarray([case.beta_deg for case in cases], dtype=float),
        "case_eta": np.asarray([case.eta for case in cases], dtype=float),
        "case_epsilon": np.asarray([case.epsilon for case in cases], dtype=float),
    }
    for case_index, values in enumerate(mu_values_by_case):
        arrays[f"mu_grid_{case_index}"] = np.asarray(values, dtype=float)
    arrays.update(root_map_arrays(cases, mu_values_by_case, args, root_map))
    arrays.update(tables_arrays(cases, mu_values_by_case, args, tables))
    np.savez_compressed(path, **arrays)


def load_cache(
    path: Path,
    cases: Sequence[CaseSpec],
    args: Args,
) -> tuple[
    list[np.ndarray],
    dict[tuple[int, str, float], beta_workflow.RootResult],
    RootTables,
    dict[str, object],
] | None:
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    try:
        settings_json = str(data["settings_json"].item())
        if settings_json != cache_settings_json(cases, args):
            return None
        metadata = json.loads(str(data["metadata_json"].item()))
        mu_values_by_case = [
            np.asarray(data[f"mu_grid_{case_index}"], dtype=float)
            for case_index in range(len(cases))
        ]
        root_map: dict[tuple[int, str, float], beta_workflow.RootResult] = {}
        for case_index, _case in enumerate(cases):
            for mu_index, mu in enumerate(mu_values_by_case[case_index]):
                for model_index, model in enumerate(MODELS):
                    root_map[root_key(case_index, model, float(mu))] = beta_workflow.RootResult(
                        roots=tuple(
                            float(value)
                            for value in np.asarray(data["root_roots"][case_index, model_index, mu_index, :], dtype=float)
                        ),
                        warnings=tuple(
                            item
                            for item in str(data["root_warnings"][case_index, model_index, mu_index]).split(" | ")
                            if item
                        ),
                        root_count_found=int(data["root_count"][case_index, model_index, mu_index]),
                        lambda_max_used=float(data["root_lambda_max"][case_index, model_index, mu_index]),
                        scan_step_used=float(data["root_scan_step"][case_index, model_index, mu_index]),
                        retry_attempted=bool(data["root_retry_attempted"][case_index, model_index, mu_index]),
                        retry_changed_value=bool(data["root_retry_changed"][case_index, model_index, mu_index]),
                        notes=tuple(
                            item
                            for item in str(data["root_notes"][case_index, model_index, mu_index]).split("; ")
                            if item
                        ),
                    )
        tables = arrays_to_tables(cases, mu_values_by_case, args, data)
    finally:
        data.close()
    return mu_values_by_case, root_map, tables, metadata


def arrays_to_tables(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    data: object,
) -> RootTables:
    raw: dict[tuple[int, str], np.ndarray] = {}
    clean: dict[tuple[int, str], np.ndarray] = {}
    status: dict[tuple[int, str], np.ndarray] = {}
    notes: dict[tuple[int, str], np.ndarray] = {}
    retry_attempted: dict[tuple[int, str], np.ndarray] = {}
    retry_fixed: dict[tuple[int, str], np.ndarray] = {}
    plotted_as_nan: dict[tuple[int, str], np.ndarray] = {}
    raw_suspicious: dict[tuple[int, str], np.ndarray] = {}
    plot_suspicious: dict[tuple[int, str], np.ndarray] = {}
    root_warning: dict[tuple[int, str], np.ndarray] = {}
    root_count: dict[tuple[int, str], np.ndarray] = {}
    for case_index, _case in enumerate(cases):
        mu_count = len(mu_values_by_case[case_index])
        for model_index, model in enumerate(MODELS):
            key = (case_index, model)
            raw[key] = np.asarray(data["table_raw"][case_index, model_index, :, :mu_count], dtype=float)
            clean[key] = np.asarray(data["table_clean"][case_index, model_index, :, :mu_count], dtype=float)
            status[key] = np.asarray(data["table_status"][case_index, model_index, :, :mu_count], dtype=object)
            notes[key] = np.asarray(data["table_notes"][case_index, model_index, :, :mu_count], dtype=object)
            retry_attempted[key] = np.asarray(
                data["table_retry_attempted"][case_index, model_index, :, :mu_count],
                dtype=bool,
            )
            retry_fixed[key] = np.asarray(data["table_retry_fixed"][case_index, model_index, :, :mu_count], dtype=bool)
            plotted_as_nan[key] = np.asarray(
                data["table_plotted_as_nan"][case_index, model_index, :, :mu_count],
                dtype=bool,
            )
            raw_suspicious[key] = np.asarray(
                data["table_raw_suspicious"][case_index, model_index, :, :mu_count],
                dtype=bool,
            )
            plot_suspicious[key] = np.asarray(
                data["table_plot_suspicious"][case_index, model_index, :, :mu_count],
                dtype=bool,
            )
            root_warning[key] = np.asarray(
                data["table_root_warning"][case_index, model_index, :, :mu_count],
                dtype=object,
            )
            root_count[key] = np.asarray(data["table_root_count"][case_index, model_index, :, :mu_count], dtype=int)
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
        root_warning=root_warning,
        root_count=root_count,
    )


def rel_diff_abs(timo: float, eb: float) -> float:
    return beta_workflow.rel_diff_abs(float(timo), float(eb))


def case_csv_path(output_dir: Path, case: CaseSpec) -> Path:
    return (
        output_dir
        / "lambda_mu_eb_vs_timo_"
        f"beta{_float_label(case.beta_deg)}_eta{_float_label(case.eta)}_eps{_float_label(case.epsilon)}.csv"
    )


def case_png_path(output_dir: Path, case: CaseSpec) -> Path:
    return case_csv_path(output_dir, case).with_suffix(".png")


def build_case_rows(
    case_index: int,
    case: CaseSpec,
    mu_values: np.ndarray,
    args: Args,
    tables: RootTables,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    eb_key = (case_index, MODEL_EB)
    timo_key = (case_index, MODEL_TIMO)
    for mu_index, mu in enumerate(mu_values):
        for sorted_zero in range(args.n_roots):
            eb_raw = float(tables.raw[eb_key][sorted_zero, mu_index])
            timo_raw = float(tables.raw[timo_key][sorted_zero, mu_index])
            eb_plot = float(tables.clean[eb_key][sorted_zero, mu_index])
            timo_plot = float(tables.clean[timo_key][sorted_zero, mu_index])
            suspicious_eb = bool(tables.plot_suspicious[eb_key][sorted_zero, mu_index])
            suspicious_timo = bool(tables.plot_suspicious[timo_key][sorted_zero, mu_index])
            retry_attempted_eb = bool(tables.retry_attempted[eb_key][sorted_zero, mu_index])
            retry_attempted_timo = bool(tables.retry_attempted[timo_key][sorted_zero, mu_index])
            retry_fixed_eb = bool(tables.retry_fixed[eb_key][sorted_zero, mu_index])
            retry_fixed_timo = bool(tables.retry_fixed[timo_key][sorted_zero, mu_index])
            plotted_nan_eb = bool(tables.plotted_as_nan[eb_key][sorted_zero, mu_index])
            plotted_nan_timo = bool(tables.plotted_as_nan[timo_key][sorted_zero, mu_index])
            notes_eb = str(tables.notes[eb_key][sorted_zero, mu_index])
            notes_timo = str(tables.notes[timo_key][sorted_zero, mu_index])
            rel_diff = rel_diff_abs(timo_plot, eb_plot)
            root_solver_warning = "; ".join(
                item
                for item in (
                    f"EB: {tables.root_warning[eb_key][sorted_zero, mu_index]}"
                    if str(tables.root_warning[eb_key][sorted_zero, mu_index])
                    else "",
                    f"Timoshenko: {tables.root_warning[timo_key][sorted_zero, mu_index]}"
                    if str(tables.root_warning[timo_key][sorted_zero, mu_index])
                    else "",
                )
                if item
            )
            rows.append(
                {
                    "mu": float(mu),
                    "beta_deg": case.beta_deg,
                    "eta": case.eta,
                    "epsilon": case.epsilon,
                    "sorted_index": sorted_zero + 1,
                    "Lambda_EB_raw": eb_raw,
                    "Lambda_Timoshenko_raw": timo_raw,
                    "Lambda_EB_plot": eb_plot,
                    "Lambda_Timoshenko_plot": timo_plot,
                    "rel_diff_Timo_vs_EB": rel_diff,
                    "rel_diff_abs_Timoshenko_vs_EB": rel_diff,
                    "suspicious_EB": suspicious_eb,
                    "suspicious_Timoshenko": suspicious_timo,
                    "retry_attempted": retry_attempted_eb or retry_attempted_timo,
                    "retry_fixed": retry_fixed_eb or retry_fixed_timo,
                    "plotted_as_nan": plotted_nan_eb or plotted_nan_timo,
                    "root_solver_warning": root_solver_warning,
                    "notes": "; ".join(
                        item
                        for item in (
                            f"EB: {notes_eb}" if notes_eb and notes_eb != "ok" else "",
                            f"Timoshenko: {notes_timo}" if notes_timo and notes_timo != "ok" else "",
                        )
                        if item
                    )
                    or "ok",
                    "status_EB": str(tables.status[eb_key][sorted_zero, mu_index]),
                    "status_Timoshenko": str(tables.status[timo_key][sorted_zero, mu_index]),
                    "root_warning_EB": str(tables.root_warning[eb_key][sorted_zero, mu_index]),
                    "root_warning_Timoshenko": str(tables.root_warning[timo_key][sorted_zero, mu_index]),
                    "root_count_EB": int(tables.root_count[eb_key][sorted_zero, mu_index]),
                    "root_count_Timoshenko": int(tables.root_count[timo_key][sorted_zero, mu_index]),
                    "retry_attempted_EB": retry_attempted_eb,
                    "retry_attempted_Timoshenko": retry_attempted_timo,
                    "retry_fixed_EB": retry_fixed_eb,
                    "retry_fixed_Timoshenko": retry_fixed_timo,
                    "plotted_as_nan_EB": plotted_nan_eb,
                    "plotted_as_nan_Timoshenko": plotted_nan_timo,
                    "notes_EB": notes_eb,
                    "notes_Timoshenko": notes_timo,
                }
            )
    return rows


def build_summary_rows(
    case_rows_by_index: dict[int, list[dict[str, object]]],
    cases: Sequence[CaseSpec],
    tables: RootTables,
    initial_raw_suspicious_by_case: dict[int, int],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_index, case in enumerate(cases):
        case_rows = case_rows_by_index[case_index]
        sorted_indices = sorted({int(row["sorted_index"]) for row in case_rows})
        for sorted_index in sorted_indices:
            selected = [row for row in case_rows if int(row["sorted_index"]) == sorted_index]
            rel_values = [
                (float(row["mu"]), float(row["rel_diff_Timo_vs_EB"]))
                for row in selected
                if isfinite(float(row["rel_diff_Timo_vs_EB"]))
            ]
            if rel_values:
                mu_at_max, max_rel = max(rel_values, key=lambda item: item[1])
                mean_rel = float(np.mean([item[1] for item in rel_values]))
            else:
                mu_at_max, max_rel, mean_rel = float("nan"), float("nan"), float("nan")
            sorted_zero = sorted_index - 1
            raw_suspicious = sum(
                int(np.count_nonzero(tables.raw_suspicious[(case_index, model)][sorted_zero, :]))
                for model in MODELS
            )
            final_suspicious = sum(
                int(np.count_nonzero(tables.plot_suspicious[(case_index, model)][sorted_zero, :]))
                for model in MODELS
            )
            retry_fixed = sum(1 for row in selected if bool(row["retry_fixed"]))
            nan_count = sum(1 for row in selected if bool(row["plotted_as_nan"]))
            notes = (
                "ok"
                if raw_suspicious == 0 and final_suspicious == 0 and nan_count == 0
                else (
                    f"initial raw suspicious points in case={int(initial_raw_suspicious_by_case.get(case_index, 0))}; "
                    f"raw suspicious after refinement for sorted index={raw_suspicious}; "
                    f"final plotted suspicious={final_suspicious}; plotted NaN rows={nan_count}"
                )
            )
            rows.append(
                {
                    "beta_deg": case.beta_deg,
                    "eta": case.eta,
                    "epsilon": case.epsilon,
                    "sorted_index": sorted_index,
                    "max_rel_diff_over_mu": max_rel,
                    "mean_rel_diff_over_mu": mean_rel,
                    "mu_at_max_rel_diff": mu_at_max,
                    "initial_raw_suspicious_point_count": int(initial_raw_suspicious_by_case.get(case_index, 0)),
                    "raw_suspicious_point_count": raw_suspicious,
                    "suspicious_point_count": final_suspicious,
                    "retry_fixed_count": retry_fixed,
                    "nan_count_after_cleanup": nan_count,
                    "notes": notes,
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
    mu_values: np.ndarray,
    args: Args,
    tables: RootTables,
) -> Path:
    fig, ax = plt.subplots(figsize=(10.6, 6.4), constrained_layout=True)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(args.n_roots, 2)))
    for sorted_index in range(1, args.n_roots + 1):
        color = colors[sorted_index - 1]
        ax.plot(
            mu_values,
            plotted_series(tables, case_index=case_index, model=MODEL_TIMO, sorted_index=sorted_index),
            color=color,
            linestyle="-",
            linewidth=1.65,
            alpha=0.95,
        )
        ax.plot(
            mu_values,
            plotted_series(tables, case_index=case_index, model=MODEL_EB, sorted_index=sorted_index),
            color=color,
            linestyle="--",
            linewidth=1.3,
            alpha=0.85,
        )
    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    ax.set_title(
        "Sorted in-plane Lambda(mu): EB dashed vs Timoshenko solid\n"
        f"beta={case.beta_deg:g} deg, eta={case.eta:g}, epsilon={case.epsilon:g}"
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
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    tables: RootTables,
) -> Path:
    epsilons: list[float] = []
    betas: list[float] = []
    for case in cases:
        if not any(abs(case.epsilon - existing) <= 1.0e-12 for existing in epsilons):
            epsilons.append(float(case.epsilon))
        if not any(abs(case.beta_deg - existing) <= 1.0e-12 for existing in betas):
            betas.append(float(case.beta_deg))
    fig, axes = plt.subplots(
        len(epsilons),
        len(betas),
        figsize=(5.0 * len(betas), 3.8 * len(epsilons)),
        squeeze=False,
        constrained_layout=True,
    )
    case_lookup = {(case.beta_deg, case.epsilon): index for index, case in enumerate(cases)}
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(args.n_roots, 2)))
    for row_index, epsilon in enumerate(epsilons):
        for col_index, beta_deg in enumerate(betas):
            ax = axes[row_index][col_index]
            case_index = case_lookup.get((float(beta_deg), float(epsilon)))
            if case_index is None:
                ax.set_axis_off()
                continue
            mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
            for sorted_index in range(1, args.n_roots + 1):
                color = colors[sorted_index - 1]
                ax.plot(
                    mu_values,
                    plotted_series(tables, case_index=case_index, model=MODEL_TIMO, sorted_index=sorted_index),
                    color=color,
                    linestyle="-",
                    linewidth=1.05,
                    alpha=0.95,
                )
                ax.plot(
                    mu_values,
                    plotted_series(tables, case_index=case_index, model=MODEL_EB, sorted_index=sorted_index),
                    color=color,
                    linestyle="--",
                    linewidth=0.95,
                    alpha=0.8,
                )
            ax.set_title(f"beta={float(beta_deg):g} deg, eta={args.eta:g}, eps={float(epsilon):g}", fontsize=10)
            ax.grid(True, alpha=0.2)
            if row_index == len(epsilons) - 1:
                ax.set_xlabel("mu")
            if col_index == 0:
                ax.set_ylabel("Lambda")
    fig.suptitle("EB dashed vs Timoshenko solid, sorted in-plane Lambda(mu)", fontsize=12)
    path = overview_png_path(args.output_dir, cases)
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def build_spike_audit_rows(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    tables: RootTables,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_index, case in enumerate(cases):
        mu_values = np.asarray(mu_values_by_case[case_index], dtype=float)
        for model in MODELS:
            key = (case_index, model)
            for sorted_zero in range(args.n_roots):
                for mu_index, mu in enumerate(mu_values):
                    raw_prev, raw_next, raw_rel_prev, raw_rel_next = row_jump_metrics(
                        tables.raw[key],
                        sorted_zero,
                        mu_index,
                    )
                    plot_prev, plot_next, plot_rel_prev, plot_rel_next = row_jump_metrics(
                        tables.clean[key],
                        sorted_zero,
                        mu_index,
                    )
                    rows.append(
                        {
                            "case_index": case_index + 1,
                            "beta_deg": case.beta_deg,
                            "eta": case.eta,
                            "epsilon": case.epsilon,
                            "model": model,
                            "sorted_index": sorted_zero + 1,
                            "mu": float(mu),
                            "Lambda_raw": float(tables.raw[key][sorted_zero, mu_index]),
                            "Lambda_plot": float(tables.clean[key][sorted_zero, mu_index]),
                            "jump_prev_raw": raw_prev,
                            "jump_next_raw": raw_next,
                            "jump_rel_prev_raw": raw_rel_prev,
                            "jump_rel_next_raw": raw_rel_next,
                            "suspicious_raw": bool(tables.raw_suspicious[key][sorted_zero, mu_index]),
                            "retry_attempted": bool(tables.retry_attempted[key][sorted_zero, mu_index]),
                            "retry_fixed": bool(tables.retry_fixed[key][sorted_zero, mu_index]),
                            "plotted_as_nan": bool(tables.plotted_as_nan[key][sorted_zero, mu_index]),
                            "jump_prev_plot": plot_prev,
                            "jump_next_plot": plot_next,
                            "jump_rel_prev_plot": plot_rel_prev,
                            "jump_rel_next_plot": plot_rel_next,
                            "suspicious_plot": bool(tables.plot_suspicious[key][sorted_zero, mu_index]),
                            "notes": str(tables.notes[key][sorted_zero, mu_index]),
                        }
                    )
    return rows


def warning_counter(root_map: dict[tuple[int, str, float], beta_workflow.RootResult]) -> Counter[str]:
    counter: Counter[str] = Counter()
    for result in root_map.values():
        for warning in result.warnings:
            counter[warning] += 1
    return counter


def missing_root_points(root_map: dict[tuple[int, str, float], beta_workflow.RootResult], n_roots: int) -> int:
    return sum(max(0, int(n_roots) - int(result.root_count_found)) for result in root_map.values())


def note_contains(text: object, needle: str) -> bool:
    return needle.lower() in str(text).lower()


def count_model_fallbacks(tables: RootTables, case_index: int, model: str) -> int:
    return sum(
        1
        for value in np.ravel(tables.notes[(case_index, model)])
        if note_contains(value, "continuation fallback")
    )


def count_table_flags(tables: RootTables, case_index: int, model: str, attribute: str) -> int:
    array = getattr(tables, attribute)[(case_index, model)]
    return int(np.count_nonzero(np.asarray(array, dtype=bool)))


def count_svd_notes(tables: RootTables, case_index: int, model: str) -> int:
    return sum(1 for value in np.ravel(tables.notes[(case_index, model)]) if note_contains(value, "SVD recovery"))


def markdown_summary_table(summary_rows: Sequence[dict[str, object]]) -> list[str]:
    lines = [
        "| beta | eta | epsilon | sorted | max rel diff | mean rel diff | mu at max | initial suspicious | final suspicious | retry fixed | NaN after cleanup |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            f"{float(row['beta_deg']):.6g} | "
            f"{float(row['eta']):.6g} | "
            f"{float(row['epsilon']):.6g} | "
            f"{int(row['sorted_index'])} | "
            f"{float(row['max_rel_diff_over_mu']):.6g} | "
            f"{float(row['mean_rel_diff_over_mu']):.6g} | "
            f"{float(row['mu_at_max_rel_diff']):.6g} | "
            f"{int(row['initial_raw_suspicious_point_count'])} | "
            f"{int(row['suspicious_point_count'])} | "
            f"{int(row['retry_fixed_count'])} | "
            f"{int(row['nan_count_after_cleanup'])} |"
        )
    return lines


def overall_observation_lines(summary_rows: Sequence[dict[str, object]]) -> list[str]:
    lines: list[str] = []
    by_epsilon: dict[float, list[float]] = {}
    for row in summary_rows:
        value = float(row["max_rel_diff_over_mu"])
        if isfinite(value):
            by_epsilon.setdefault(float(row["epsilon"]), []).append(value)
    for epsilon in sorted(by_epsilon):
        values = by_epsilon[epsilon]
        lines.append(
            f"- epsilon={epsilon:g}: largest case/index max relative difference is "
            f"{max(values):.6g}; mean over case/index maxima is {float(np.mean(values)):.6g}."
        )
    epsilon_keys = sorted(by_epsilon)
    if len(epsilon_keys) >= 2:
        small_epsilon = epsilon_keys[0]
        large_epsilon = epsilon_keys[-1]
        small = by_epsilon[small_epsilon]
        large = by_epsilon[large_epsilon]
        relation = "larger" if max(large) > max(small) else "not larger"
        lines.append(
            f"- The epsilon={large_epsilon:g} differences are {relation} than the "
            f"epsilon={small_epsilon:g} differences in this diagnostic set."
        )
    lines.append(
        "- This fixed-beta sorted-frequency report does not assign modal-character branches."
    )
    lines.extend(
        [
            "- Smaller epsilon values are expected to keep Euler-Bernoulli and Timoshenko closer.",
            "- As epsilon grows across the scan, EB/Timoshenko differences should become more visible.",
            "- Higher sorted frequencies are expected to show larger Timoshenko corrections.",
        ]
    )
    return lines


def write_report(
    path: Path,
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
    tables: RootTables,
    summary_rows: Sequence[dict[str, object]],
    output_paths: Sequence[Path],
    timing: TimingStats,
    refinement_windows: Sequence[RefinementWindow],
    initial_raw_suspicious_count: int,
) -> Path:
    warnings = warning_counter(root_map)
    total_raw = sum(int(np.count_nonzero(tables.raw_suspicious[key])) for key in tables.raw_suspicious)
    total_plot = sum(int(np.count_nonzero(tables.plot_suspicious[key])) for key in tables.plot_suspicious)
    total_retry_fixed = sum(int(np.count_nonzero(tables.retry_fixed[key])) for key in tables.retry_fixed)
    total_nan = sum(int(np.count_nonzero(tables.plotted_as_nan[key])) for key in tables.plotted_as_nan)
    missing = missing_root_points(root_map, args.n_roots)

    lines = [
        "# EB vs Timoshenko Lambda(mu) Diagnostic Cases",
        "",
        "Diagnostic only. These plots use sorted in-plane frequencies at each mu, not descendant branch tracking.",
        "The workflow reuses the existing Euler-Bernoulli thickness-mismatch root helper and the existing variable-length Timoshenko helper.",
        "No analytic formulas, determinants, root solvers, Timoshenko shear coefficient/k', FEM workflows, article files, article figures, or baseline results were changed.",
        "",
        "## Cases",
        "",
        "| # | beta deg | eta | epsilon | mu points |",
        "| ---: | ---: | ---: | ---: | ---: |",
    ]
    for case_index, case in enumerate(cases):
        lines.append(
            f"| {case_index + 1} | {case.beta_deg:g} | {case.eta:g} | {case.epsilon:g} | "
            f"{len(mu_values_by_case[case_index])} |"
        )

    lines.extend(
        [
            "",
            "## Mu Grid",
            "",
            f"- base mu range: {args.mu_min:g} to {args.mu_max:g}",
            f"- base mu step: {args.mu_step:g}",
            f"- local refined mu step: {args.local_mu_step:g}",
            f"- sorted frequencies per theory: {args.n_roots}",
            f"- spike thresholds: jump_abs > {args.jump_abs_threshold:g} or jump_rel > {args.jump_rel_threshold:g}",
            "",
            "Refined windows:",
            "",
        ]
    )
    if refinement_windows:
        for window in refinement_windows:
            case = cases[window.case_index]
            lines.append(
                f"- beta={case.beta_deg:g}, eta={case.eta:g}, epsilon={case.epsilon:g}: "
                f"mu {window.mu_min:.8g}..{window.mu_max:.8g} at step {args.local_mu_step:g} ({window.reason})."
            )
    else:
        lines.append("- none; the base mu grid did not trigger local refinement.")

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
            f"- missing root slots after initial/retry solving: {missing}",
            f"- initial raw suspicious points detected on base grid: {initial_raw_suspicious_count}",
            f"- raw suspicious points detected after local refinement: {total_raw}",
            f"- retry-fixed points: {total_retry_fixed}",
            f"- unresolved NaN statuses after cleanup: {total_nan}",
            f"- final plotted suspicious points: {total_plot}",
            "",
        ]
    )
    if warnings:
        lines.append("Warning classes:")
        for warning, count in warnings.most_common():
            lines.append(f"- {count} x {warning}")
    else:
        lines.append("No root-search or Timoshenko basis warnings were recorded.")

    lines.extend(["", "## Numeric Summary", "", *markdown_summary_table(summary_rows), ""])
    lines.extend(["## General Diagnostic Observations", "", *overall_observation_lines(summary_rows), ""])
    lines.append(
        "Visual artifact status: "
        + (
            "no artificial jump artifacts remain in plotted curves by the configured check; "
            "unresolved points, if any, are plotted as NaN rather than false spikes."
            if total_plot == 0
            else "inspect spike audit CSV; plotted suspicious points remain."
        )
    )
    lines.extend(
        [
            "",
            "## Explicit Confirmations",
            "",
            "- Sorted frequencies were used.",
            "- Descendant tracking was not used.",
            "- No single-rod reference curves were plotted.",
            "- FEM, 3D FEM, Gmsh, and CalculiX were not run or used.",
            "- Formulas, determinants, unknown ordering, root solvers, Timoshenko shear coefficient/k', and protected article/baseline areas were not changed.",
        ]
    )
    lines.extend(["", "## Outputs", ""])
    for output_path in output_paths:
        lines.append(f"- `{_rel(output_path)}`")
    lines.extend(
        [
            "",
            "## Timing And Cache",
            "",
            f"- timing report: `{_rel(args.output_dir / TIMING_REPORT_NAME)}`",
            f"- total runtime: {timing.total_seconds:.6g} s",
            f"- cache hit: {timing.cache_hit}",
            f"- ordinary compute seconds: {timing.ordinary_compute_seconds:.6g}",
            f"- spike audit seconds: {timing.spike_audit_seconds:.6g}",
            f"- repair seconds: {timing.repair_seconds:.6g}",
            f"- plotting seconds: {timing.plotting_seconds:.6g}",
            f"- retry-fixed values: {timing.repair_count}",
            f"- explicit local SVD recovery calls: {timing.sv_recovery_calls}",
            f"- Timoshenko continuation global fallbacks: {timing.fallback_count}",
            "- Cache and plot-only modes avoid root recomputation when settings match.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_outputs(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
    tables: RootTables,
    timing: TimingStats,
    refinement_windows: Sequence[RefinementWindow],
    initial_raw_suspicious_count: int,
    initial_raw_suspicious_by_case: dict[int, int],
) -> tuple[list[Path], list[dict[str, object]]]:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    case_rows_by_index: dict[int, list[dict[str, object]]] = {}
    output_paths: list[Path] = []
    for case_index, case in enumerate(cases):
        rows = build_case_rows(case_index, case, np.asarray(mu_values_by_case[case_index], dtype=float), args, tables)
        case_rows_by_index[case_index] = rows
        csv_path = case_csv_path(args.output_dir, case)
        _write_csv(csv_path, rows, CASE_CSV_FIELDS)
        output_paths.append(csv_path)
        output_paths.append(plot_case(case_index, case, np.asarray(mu_values_by_case[case_index], dtype=float), args, tables))
    summary_rows = build_summary_rows(case_rows_by_index, cases, tables, initial_raw_suspicious_by_case)
    summary_path = summary_csv_path(args.output_dir, cases)
    _write_csv(summary_path, summary_rows, SUMMARY_FIELDS)
    output_paths.append(summary_path)
    spike_path = spike_audit_csv_path(args.output_dir, cases)
    _write_csv(spike_path, build_spike_audit_rows(cases, mu_values_by_case, args, tables), SPIKE_AUDIT_FIELDS)
    output_paths.append(spike_path)
    overview_path = plot_overview(cases, mu_values_by_case, args, tables)
    output_paths.append(overview_path)
    report_path = write_report(
        report_md_path(args.output_dir, cases),
        cases,
        mu_values_by_case,
        args,
        root_map,
        tables,
        summary_rows,
        output_paths,
        timing,
        refinement_windows,
        initial_raw_suspicious_count,
    )
    output_paths.append(report_path)
    return output_paths, summary_rows


def build_timing_rows(
    cases: Sequence[CaseSpec],
    mu_values_by_case: Sequence[np.ndarray],
    args: Args,
    root_map: dict[tuple[int, str, float], beta_workflow.RootResult],
    tables: RootTables,
    timing: TimingStats,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case_index, case in enumerate(cases):
        for model in MODELS:
            warnings_count = sum(
                1
                for mu in np.asarray(mu_values_by_case[case_index], dtype=float)
                if root_map[root_key(case_index, model, float(mu))].warnings
            )
            root_mode = "global" if model == MODEL_EB else args.timo_root_mode
            rows.append(
                {
                    "run_id": timing.run_id,
                    "model": model,
                    "beta_deg": case.beta_deg,
                    "eta": case.eta,
                    "epsilon": case.epsilon,
                    "n_mu_points": len(mu_values_by_case[case_index]),
                    "n_roots": args.n_roots,
                    "root_mode": root_mode,
                    "cache_hit": timing.cache_hit,
                    "ordinary_compute_seconds": timing.ordinary_compute_seconds,
                    "spike_audit_seconds": timing.spike_audit_seconds,
                    "repair_seconds": timing.repair_seconds,
                    "plotting_seconds": timing.plotting_seconds,
                    "total_seconds": timing.total_seconds,
                    "warnings": warnings_count,
                    "fallback_count": count_model_fallbacks(tables, case_index, model),
                    "repair_count": count_table_flags(tables, case_index, model, "retry_fixed"),
                    "sv_recovery_calls": count_svd_notes(tables, case_index, model),
                    "notes": "; ".join(timing.notes),
                }
            )
    return rows


def write_timing_report(path: Path, rows: Sequence[dict[str, object]]) -> Path:
    _write_csv(path, rows, TIMING_FIELDS)
    return path


def print_run_summary(
    output_paths: Sequence[Path],
    summary_rows: Sequence[dict[str, object]],
    timing: TimingStats,
) -> None:
    print("generated EB-vs-Timoshenko Lambda(mu) diagnostic outputs:")
    for path in output_paths:
        print(f"  {_rel(path)}")
    finite_max = [
        float(row["max_rel_diff_over_mu"])
        for row in summary_rows
        if isfinite(float(row["max_rel_diff_over_mu"]))
    ]
    if finite_max:
        print(f"largest max relative difference over all rows: {max(finite_max):.8g}")
    print(
        "timing: "
        f"total={timing.total_seconds:.3f}s, compute={timing.ordinary_compute_seconds:.3f}s, "
        f"audit={timing.spike_audit_seconds:.3f}s, repair={timing.repair_seconds:.3f}s, "
        f"plot={timing.plotting_seconds:.3f}s, cache_hit={timing.cache_hit}, "
        f"fallbacks={timing.fallback_count}, repairs={timing.repair_count}, "
        f"sv_recovery={timing.sv_recovery_calls}"
    )


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    cases = default_cases(args.beta_values, args.eta, args.epsilon_values)
    cache_path = cache_file_path(cases, args)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    total_start = time.perf_counter()

    cache_loaded = None
    if args.reuse_cache and not args.force_recompute:
        cache_loaded = load_cache(cache_path, cases, args)
    if args.plot_only and cache_loaded is None:
        raise FileNotFoundError(
            f"--plot-only requires a matching cache, but none was found at {cache_path}"
        )

    refinement_windows: tuple[RefinementWindow, ...] = ()
    initial_raw_suspicious_count = 0
    initial_raw_suspicious_by_case = {case_index: 0 for case_index in range(len(cases))}
    compute_seconds = 0.0
    spike_audit_seconds = 0.0
    repair_seconds = 0.0
    sv_recovery_calls = 0

    if cache_loaded is not None:
        mu_values_by_case, root_map, tables, metadata = cache_loaded
        refinement_windows = tuple(
            RefinementWindow(
                case_index=int(item["case_index"]),
                mu_min=float(item["mu_min"]),
                mu_max=float(item["mu_max"]),
                reason=str(item["reason"]),
            )
            for item in metadata.get("refinement_windows", [])
        )
        initial_raw_suspicious_count = int(metadata.get("initial_raw_suspicious_count", 0))
        initial_raw_suspicious_by_case = {
            int(key): int(value)
            for key, value in dict(metadata.get("initial_raw_suspicious_by_case", {})).items()
        }
        sv_recovery_calls = int(metadata.get("sv_recovery_calls", 0))
        notes = (f"loaded cache {cache_path}",)
        cache_hit = True
    else:
        compute_start = time.perf_counter()
        base_grids = [base_mu_grid(args) for _case in cases]
        root_map = compute_roots_for_grids(cases, base_grids, args)
        compute_seconds = time.perf_counter() - compute_start

        audit_start = time.perf_counter()
        base_raw, _base_warnings, _base_counts, _base_notes = raw_arrays_from_root_map(cases, base_grids, args, root_map)
        initial_raw_suspicious_by_case = {}
        for case_index, _case in enumerate(cases):
            count = sum(
                int(np.count_nonzero(suspicious_mask(base_raw[(case_index, model)], args)))
                for model in MODELS
            )
            initial_raw_suspicious_by_case[case_index] = count
        initial_raw_suspicious_count = sum(initial_raw_suspicious_by_case.values())
        refinement_windows = refinement_windows_from_raw(cases, base_grids, args, base_raw)
        refined_grids = refine_mu_grids(cases, base_grids, args, refinement_windows)
        spike_audit_seconds += time.perf_counter() - audit_start

        if any(len(refined_grids[index]) > len(base_grids[index]) for index in range(len(cases))):
            compute_start = time.perf_counter()
            root_map = compute_roots_for_grids(cases, refined_grids, args, existing=root_map)
            compute_seconds += time.perf_counter() - compute_start
        tables, timing_parts, sv_recovery_calls = make_root_tables(
            cases,
            refined_grids,
            args,
            root_map,
            initial_raw_suspicious_count=initial_raw_suspicious_count,
        )
        spike_audit_seconds += float(timing_parts["spike_audit_seconds"])
        repair_seconds += float(timing_parts["repair_seconds"])
        mu_values_by_case = refined_grids
        metadata = {
            "cache_version": CACHE_VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "initial_raw_suspicious_count": initial_raw_suspicious_count,
            "initial_raw_suspicious_by_case": {
                str(key): int(value) for key, value in initial_raw_suspicious_by_case.items()
            },
            "refinement_windows": [
                {
                    "case_index": window.case_index,
                    "mu_min": window.mu_min,
                    "mu_max": window.mu_max,
                    "reason": window.reason,
                }
                for window in refinement_windows
            ],
            "sv_recovery_calls": sv_recovery_calls,
        }
        save_cache(cache_path, cases, mu_values_by_case, args, root_map, tables, metadata)
        notes = (f"wrote cache {cache_path}",)
        cache_hit = False

    plotting_start = time.perf_counter()
    provisional_timing = TimingStats(
        run_id=run_id,
        cache_hit=cache_hit,
        ordinary_compute_seconds=compute_seconds,
        spike_audit_seconds=spike_audit_seconds,
        repair_seconds=repair_seconds,
        plotting_seconds=0.0,
        total_seconds=0.0,
        fallback_count=sum(count_model_fallbacks(tables, case_index, MODEL_TIMO) for case_index in range(len(cases))),
        repair_count=sum(int(np.count_nonzero(array)) for array in tables.retry_fixed.values()),
        sv_recovery_calls=sv_recovery_calls,
        notes=notes,
    )
    output_paths, summary_rows = write_outputs(
        cases,
        mu_values_by_case,
        args,
        root_map,
        tables,
        provisional_timing,
        refinement_windows,
        initial_raw_suspicious_count,
        initial_raw_suspicious_by_case,
    )
    plotting_seconds = time.perf_counter() - plotting_start
    timing = TimingStats(
        run_id=run_id,
        cache_hit=cache_hit,
        ordinary_compute_seconds=compute_seconds,
        spike_audit_seconds=spike_audit_seconds,
        repair_seconds=repair_seconds,
        plotting_seconds=plotting_seconds,
        total_seconds=time.perf_counter() - total_start,
        fallback_count=sum(count_model_fallbacks(tables, case_index, MODEL_TIMO) for case_index in range(len(cases))),
        repair_count=sum(int(np.count_nonzero(array)) for array in tables.retry_fixed.values()),
        sv_recovery_calls=sv_recovery_calls,
        notes=notes,
    )
    timing_path = write_timing_report(
        args.output_dir / TIMING_REPORT_NAME,
        build_timing_rows(cases, mu_values_by_case, args, root_map, tables, timing),
    )
    if timing_path not in output_paths:
        output_paths.append(timing_path)
    write_report(
        report_md_path(args.output_dir, cases),
        cases,
        mu_values_by_case,
        args,
        root_map,
        tables,
        summary_rows,
        output_paths,
        timing,
        refinement_windows,
        initial_raw_suspicious_count,
    )
    print_run_summary(output_paths, summary_rows, timing)


if __name__ == "__main__":
    main()
