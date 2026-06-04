from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    find_roots_scan_bisect_eta,
    thickness_mismatch_factors,
)
from my_project.analytic.solvers_out_of_plane import (  # noqa: E402
    find_first_n_roots_out_of_plane_with_warnings,
)


DEFAULT_EPSILON = 0.0025
DEFAULT_POISSON = 0.3
DEFAULT_MU_VALUES = (0.0, 0.3, 0.6)
DEFAULT_ETA_VALUES = (-0.5, 0.0, 0.5)
DEFAULT_BETA_MIN = 0.0
DEFAULT_BETA_MAX = 90.0
DEFAULT_BETA_STEP = 0.25
DEFAULT_N_ROOTS = 16
DEFAULT_LAMBDA_MAX = 35.0
DEFAULT_LAMBDA_MAX_RETRY = 45.0
DEFAULT_ROOT_SCAN_DENSITY = 0.02
DEFAULT_RETRY_SCAN_DENSITY = 0.01
DEFAULT_JUMP_ABS_THRESHOLD = 1.0
DEFAULT_JUMP_REL_THRESHOLD = 0.15
DEFAULT_OUTPUT_DIR = Path("results") / "in_plane_vs_out_of_plane_lambda_beta"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "in_plane_vs_out_of_plane_lambda_beta"
SUBSYSTEMS = ("in_plane", "out_of_plane")


@dataclass(frozen=True)
class Args:
    epsilon: float
    poisson: float
    mu_values: tuple[float, ...]
    eta_values: tuple[float, ...]
    beta_min: float
    beta_max: float
    beta_step: float
    n_roots: int
    lambda_max: float
    lambda_max_retry: float
    root_scan_density: float
    retry_scan_density: float
    jump_abs_threshold: float
    jump_rel_threshold: float
    output_dir: Path
    smoke: bool
    force: bool


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


def _repo_output_dir(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def _param_label(value: float) -> str:
    return f"{float(value):.1f}".replace("-", "m").replace(".", "p")


def _fmt(value: object) -> object:
    if isinstance(value, float):
        value_f = float(value)
        if not isfinite(value_f):
            return "nan"
        return f"{value_f:.16e}"
    return value


def rel(path: Path) -> str:
    return str(Path(path).relative_to(REPO_ROOT))


def write_csv(path: Path, rows: Sequence[dict[str, object]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _fmt(row.get(key, "")) for key in fieldnames})


def _parse_float_tuple(values: list[str], name: str) -> tuple[float, ...]:
    parsed = tuple(float(value) for value in values)
    if not parsed:
        raise ValueError(f"{name} must contain at least one value.")
    if not all(isfinite(value) for value in parsed):
        raise ValueError(f"{name} values must be finite.")
    return parsed


def beta_grid(beta_min: float, beta_max: float, beta_step: float) -> np.ndarray:
    beta_min_f = float(beta_min)
    beta_max_f = float(beta_max)
    beta_step_f = float(beta_step)
    if not (isfinite(beta_min_f) and isfinite(beta_max_f) and isfinite(beta_step_f)):
        raise ValueError("beta grid values must be finite.")
    if beta_step_f <= 0.0:
        raise ValueError("beta_step must be positive.")
    if beta_max_f < beta_min_f:
        raise ValueError("beta_max must be greater than or equal to beta_min.")
    count = int(np.floor((beta_max_f - beta_min_f) / beta_step_f + 0.5)) + 1
    values = beta_min_f + beta_step_f * np.arange(count, dtype=float)
    if values[-1] < beta_max_f - 1.0e-10:
        values = np.append(values, beta_max_f)
    values[-1] = min(values[-1], beta_max_f)
    return values


def validate_args(args: Args) -> None:
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not isfinite(args.poisson) or args.poisson <= -1.0:
        raise ValueError("poisson must be finite and greater than -1.")
    if args.n_roots <= 0:
        raise ValueError("n_roots must be positive.")
    if args.lambda_max <= 0.0 or args.lambda_max_retry <= 0.0:
        raise ValueError("lambda maxima must be positive.")
    if args.lambda_max_retry < args.lambda_max:
        raise ValueError("lambda_max_retry must be at least lambda_max.")
    if args.root_scan_density <= 0.0 or args.retry_scan_density <= 0.0:
        raise ValueError("scan densities are scan steps and must be positive.")
    if args.jump_abs_threshold <= 0.0 or args.jump_rel_threshold <= 0.0:
        raise ValueError("jump thresholds must be positive.")
    beta_grid(args.beta_min, args.beta_max, args.beta_step)
    for mu in args.mu_values:
        if not (-1.0 < float(mu) < 1.0):
            raise ValueError("all mu values must lie inside (-1, 1).")
    for mu in args.mu_values:
        for eta in args.eta_values:
            thickness_mismatch_factors(float(mu), float(eta))


def parse_args(argv: list[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        description="Plot sorted in-plane vs out-of-plane analytic Lambda(beta) diagnostics."
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--mu-values", nargs="+", default=[str(value) for value in DEFAULT_MU_VALUES])
    parser.add_argument("--eta-values", nargs="+", default=[str(value) for value in DEFAULT_ETA_VALUES])
    parser.add_argument("--beta-min", type=float, default=DEFAULT_BETA_MIN)
    parser.add_argument("--beta-max", type=float, default=DEFAULT_BETA_MAX)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--n-roots", type=int, default=DEFAULT_N_ROOTS)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--lambda-max-retry", type=float, default=DEFAULT_LAMBDA_MAX_RETRY)
    parser.add_argument(
        "--root-scan-density",
        type=float,
        default=DEFAULT_ROOT_SCAN_DENSITY,
        help="Root scan step; smaller values are denser.",
    )
    parser.add_argument(
        "--retry-scan-density",
        type=float,
        default=DEFAULT_RETRY_SCAN_DENSITY,
        help="Retry root scan step; smaller values are denser.",
    )
    parser.add_argument("--jump-abs-threshold", type=float, default=DEFAULT_JUMP_ABS_THRESHOLD)
    parser.add_argument("--jump-rel-threshold", type=float, default=DEFAULT_JUMP_REL_THRESHOLD)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--force", action="store_true")
    ns = parser.parse_args(argv)

    mu_values = _parse_float_tuple(ns.mu_values, "mu-values")
    eta_values = _parse_float_tuple(ns.eta_values, "eta-values")
    beta_min = float(ns.beta_min)
    beta_max = float(ns.beta_max)
    beta_step = float(ns.beta_step)
    n_roots = int(ns.n_roots)
    output_dir = _repo_output_dir(Path(ns.output_dir))
    if ns.smoke:
        mu_values = (0.0, 0.3)
        eta_values = (0.0,)
        beta_min = 0.0
        beta_max = 90.0
        beta_step = 30.0
        n_roots = 6
        output_dir = _repo_output_dir(SMOKE_OUTPUT_DIR)

    args = Args(
        epsilon=float(ns.epsilon),
        poisson=float(ns.poisson),
        mu_values=mu_values,
        eta_values=eta_values,
        beta_min=beta_min,
        beta_max=beta_max,
        beta_step=beta_step,
        n_roots=n_roots,
        lambda_max=float(ns.lambda_max),
        lambda_max_retry=float(ns.lambda_max_retry),
        root_scan_density=float(ns.root_scan_density),
        retry_scan_density=float(ns.retry_scan_density),
        jump_abs_threshold=float(ns.jump_abs_threshold),
        jump_rel_threshold=float(ns.jump_rel_threshold),
        output_dir=output_dir,
        smoke=bool(ns.smoke),
        force=bool(ns.force),
    )
    validate_args(args)
    return args


def expected_outputs(output_dir: Path) -> list[Path]:
    return [
        output_dir / "in_plane_vs_out_of_plane_lambda_beta_sorted_roots.csv",
        output_dir / "in_plane_vs_out_of_plane_lambda_beta_case_summary.csv",
        output_dir / "in_plane_vs_out_of_plane_lambda_beta_spike_audit.csv",
        output_dir / "full_analytic_union_lambda_beta.csv",
        output_dir / "in_plane_vs_out_of_plane_lambda_beta_overview.png",
        output_dir / "in_plane_vs_out_of_plane_lambda_beta_report.md",
    ]


def prepare_output_dir(args: Args) -> None:
    existing = [path for path in expected_outputs(args.output_dir) if path.exists()]
    if existing and not args.force:
        listed = ", ".join(rel(path) for path in existing[:4])
        raise FileExistsError(f"output files already exist ({listed}); rerun with --force to overwrite them.")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def normalize_roots(raw_roots: Sequence[float], n_roots: int) -> tuple[tuple[float, ...], list[str], int]:
    notes: list[str] = []
    finite_roots = [float(root) for root in raw_roots if isfinite(float(root))]
    if len(finite_roots) != len(raw_roots):
        notes.append("non-finite roots removed")
    if any(finite_roots[index + 1] < finite_roots[index] for index in range(len(finite_roots) - 1)):
        notes.append("roots were not nondecreasing and were sorted")
        finite_roots = sorted(finite_roots)
    count_found = min(len(finite_roots), int(n_roots))
    roots = finite_roots[: int(n_roots)]
    if len(roots) < int(n_roots):
        notes.append("missing roots filled with NaN")
        roots.extend([float("nan")] * (int(n_roots) - len(roots)))
    if not notes:
        notes.append("ok")
    return tuple(roots), notes, count_found


def roots_close(left: Sequence[float], right: Sequence[float], *, atol: float = 1.0e-8, rtol: float = 1.0e-8) -> bool:
    if len(left) != len(right):
        return False
    for a, b in zip(left, right):
        if not (isfinite(float(a)) or isfinite(float(b))):
            if isfinite(float(a)) != isfinite(float(b)):
                return False
            continue
        scale = max(abs(float(a)), abs(float(b)), 1.0)
        if abs(float(a) - float(b)) > atol + rtol * scale:
            return False
    return True


def solve_roots(
    *,
    subsystem: str,
    beta_rad: float,
    mu: float,
    eta: float,
    epsilon: float,
    poisson: float,
    n_roots: int,
    lambda_max: float,
    scan_step: float,
) -> RootResult:
    warnings: list[str] = []
    if subsystem == "in_plane":
        raw_roots = find_roots_scan_bisect_eta(
            beta=float(beta_rad),
            mu=float(mu),
            epsilon=float(epsilon),
            eta=float(eta),
            n_roots=int(n_roots),
            Lmin=0.2,
            Lmax=float(lambda_max),
            scan_step=float(scan_step),
        )
    elif subsystem == "out_of_plane":
        result = find_first_n_roots_out_of_plane_with_warnings(
            beta=float(beta_rad),
            mu=float(mu),
            epsilon=float(epsilon),
            eta=float(eta),
            poisson=float(poisson),
            n_roots=int(n_roots),
            lambda_max=float(lambda_max),
            scan_step=float(scan_step),
        )
        raw_roots = result.roots
        warnings.extend(result.warnings)
    else:
        raise ValueError(f"unknown subsystem: {subsystem}")
    roots, notes, count_found = normalize_roots(raw_roots, int(n_roots))
    if count_found < int(n_roots):
        warnings.append(f"found {count_found} roots, fewer than requested {int(n_roots)}")
    return RootResult(
        roots=roots,
        warnings=tuple(warnings),
        root_count_found=count_found,
        lambda_max_used=float(lambda_max),
        scan_step_used=float(scan_step),
        retry_attempted=False,
        retry_changed_value=False,
        notes=tuple(notes),
    )


def retry_roots(
    initial: RootResult,
    *,
    subsystem: str,
    beta_rad: float,
    mu: float,
    eta: float,
    args: Args,
) -> RootResult:
    retry = solve_roots(
        subsystem=subsystem,
        beta_rad=beta_rad,
        mu=mu,
        eta=eta,
        epsilon=args.epsilon,
        poisson=args.poisson,
        n_roots=args.n_roots,
        lambda_max=args.lambda_max_retry,
        scan_step=args.retry_scan_density,
    )
    changed = not roots_close(initial.roots, retry.roots)
    use_retry = retry.root_count_found > initial.root_count_found
    notes = list(retry.notes if use_retry else initial.notes)
    notes.append(
        "retry used" if use_retry else "retry attempted but initial roots kept because count did not improve"
    )
    chosen = retry if use_retry else initial
    return RootResult(
        roots=chosen.roots,
        warnings=tuple(chosen.warnings),
        root_count_found=chosen.root_count_found,
        lambda_max_used=chosen.lambda_max_used,
        scan_step_used=chosen.scan_step_used,
        retry_attempted=True,
        retry_changed_value=changed,
        notes=tuple(notes),
    )


def compute_root_map(args: Args, beta_values: np.ndarray) -> dict[tuple[float, float, str, float], RootResult]:
    root_map: dict[tuple[float, float, str, float], RootResult] = {}
    for eta in args.eta_values:
        for mu in args.mu_values:
            for beta_deg in beta_values:
                beta_rad = float(np.deg2rad(beta_deg))
                for subsystem in SUBSYSTEMS:
                    result = solve_roots(
                        subsystem=subsystem,
                        beta_rad=beta_rad,
                        mu=float(mu),
                        eta=float(eta),
                        epsilon=args.epsilon,
                        poisson=args.poisson,
                        n_roots=args.n_roots,
                        lambda_max=args.lambda_max,
                        scan_step=args.root_scan_density,
                    )
                    if result.root_count_found < args.n_roots:
                        result = retry_roots(
                            result,
                            subsystem=subsystem,
                            beta_rad=beta_rad,
                            mu=float(mu),
                            eta=float(eta),
                            args=args,
                        )
                    root_map[(float(mu), float(eta), subsystem, float(beta_deg))] = result
    return root_map


def jump_value(current: float, neighbor: float) -> tuple[float, float]:
    if not (isfinite(float(current)) and isfinite(float(neighbor))):
        return float("nan"), float("nan")
    jump_abs = abs(float(current) - float(neighbor))
    scale = max(abs(float(neighbor)), 1.0e-12)
    return jump_abs, jump_abs / scale


def is_suspicious_jump(current: float, neighbor: float, args: Args) -> bool:
    jump_abs, jump_rel = jump_value(current, neighbor)
    if not (isfinite(jump_abs) and isfinite(jump_rel)):
        return False
    return jump_abs > args.jump_abs_threshold or jump_rel > args.jump_rel_threshold


def build_rows_and_spike_audit(
    args: Args,
    beta_values: np.ndarray,
    root_map: dict[tuple[float, float, str, float], RootResult],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    row_lookup: dict[tuple[float, float, str, int, float], dict[str, object]] = {}
    spike_rows: list[dict[str, object]] = []

    for mu in args.mu_values:
        for eta in args.eta_values:
            for subsystem in SUBSYSTEMS:
                for sorted_index in range(1, args.n_roots + 1):
                    for beta_pos, beta_deg in enumerate(beta_values):
                        beta_deg_f = float(beta_deg)
                        result = root_map[(float(mu), float(eta), subsystem, beta_deg_f)]
                        Lambda = float(result.roots[sorted_index - 1])
                        missing = not isfinite(Lambda)
                        plotted_as_nan = missing
                        notes = list(result.notes)
                        if missing:
                            notes.append("missing root not connected in plot")

                        prev_lambda = float("nan")
                        next_lambda = float("nan")
                        if beta_pos > 0:
                            prev_result = root_map[(float(mu), float(eta), subsystem, float(beta_values[beta_pos - 1]))]
                            prev_lambda = float(prev_result.roots[sorted_index - 1])
                        if beta_pos + 1 < len(beta_values):
                            next_result = root_map[(float(mu), float(eta), subsystem, float(beta_values[beta_pos + 1]))]
                            next_lambda = float(next_result.roots[sorted_index - 1])

                        suspicious_prev = is_suspicious_jump(Lambda, prev_lambda, args)
                        suspicious_next = is_suspicious_jump(Lambda, next_lambda, args)
                        suspicious = suspicious_prev or suspicious_next
                        retry_attempted = bool(result.retry_attempted)
                        retry_changed = bool(result.retry_changed_value)
                        if suspicious and isfinite(Lambda):
                            retry_attempted = True
                            retry_result = retry_roots(
                                result,
                                subsystem=subsystem,
                                beta_rad=float(np.deg2rad(beta_deg_f)),
                                mu=float(mu),
                                eta=float(eta),
                                args=args,
                            )
                            retry_lambda = float(retry_result.roots[sorted_index - 1])
                            retry_changed = retry_changed or abs(retry_lambda - Lambda) > 1.0e-8
                            retry_suspicious = is_suspicious_jump(retry_lambda, prev_lambda, args) or is_suspicious_jump(
                                retry_lambda, next_lambda, args
                            )
                            if retry_changed and not retry_suspicious:
                                result = retry_result
                                Lambda = retry_lambda
                                root_map[(float(mu), float(eta), subsystem, beta_deg_f)] = retry_result
                                suspicious = False
                                notes = list(retry_result.notes) + ["suspicious jump fixed by retry"]
                            else:
                                plotted_as_nan = True
                                notes.append("unresolved suspicious jump hidden from line plot")

                            jump_prev, _jump_prev_rel = jump_value(Lambda, prev_lambda)
                            jump_next, _jump_next_rel = jump_value(Lambda, next_lambda)
                            spike_rows.append(
                                {
                                    "mu": float(mu),
                                    "eta": float(eta),
                                    "subsystem": subsystem,
                                    "sorted_index": sorted_index,
                                    "beta_deg": beta_deg_f,
                                    "Lambda_raw": Lambda,
                                    "Lambda_prev": prev_lambda,
                                    "Lambda_next": next_lambda,
                                    "jump_prev": jump_prev,
                                    "jump_next": jump_next,
                                    "retry_attempted": retry_attempted,
                                    "retry_changed_value": retry_changed,
                                    "plotted_as_nan": plotted_as_nan,
                                    "notes": "; ".join(notes),
                                }
                            )

                        row = {
                            "beta_deg": beta_deg_f,
                            "beta_rad": float(np.deg2rad(beta_deg_f)),
                            "mu": float(mu),
                            "eta": float(eta),
                            "epsilon": args.epsilon,
                            "poisson": args.poisson,
                            "subsystem": subsystem,
                            "sorted_index": sorted_index,
                            "Lambda": Lambda,
                            "root_solver_warning": " | ".join(result.warnings),
                            "root_count_found": result.root_count_found,
                            "root_count_requested": args.n_roots,
                            "lambda_max_used": result.lambda_max_used,
                            "suspicious_jump": suspicious,
                            "plotted_as_nan": plotted_as_nan,
                            "notes": "; ".join(notes),
                        }
                        row_lookup[(float(mu), float(eta), subsystem, sorted_index, beta_deg_f)] = row

    rows = sorted(
        row_lookup.values(),
        key=lambda row: (
            float(row["mu"]),
            float(row["eta"]),
            str(row["subsystem"]),
            float(row["beta_deg"]),
            int(row["sorted_index"]),
        ),
    )
    return rows, spike_rows


def finite_lambdas(rows: Sequence[dict[str, object]]) -> list[float]:
    values: list[float] = []
    for row in rows:
        try:
            value = float(row["Lambda"])
        except (KeyError, TypeError, ValueError):
            continue
        if isfinite(value):
            values.append(value)
    return values


def build_case_summary(rows: Sequence[dict[str, object]], args: Args, beta_values: np.ndarray) -> list[dict[str, object]]:
    summaries: list[dict[str, object]] = []
    for mu in args.mu_values:
        for eta in args.eta_values:
            for subsystem in SUBSYSTEMS:
                selected = [
                    row
                    for row in rows
                    if float(row["mu"]) == float(mu)
                    and float(row["eta"]) == float(eta)
                    and str(row["subsystem"]) == subsystem
                ]
                values = finite_lambdas(selected)
                warning_count = sum(1 for row in selected if str(row.get("root_solver_warning", "")))
                missing = sum(1 for row in selected if not isfinite(float(row["Lambda"])))
                suspicious = sum(1 for row in selected if bool(row["suspicious_jump"]))
                max_step_jump = 0.0
                for sorted_index in range(1, args.n_roots + 1):
                    branch = [
                        row
                        for row in selected
                        if int(row["sorted_index"]) == sorted_index and isfinite(float(row["Lambda"]))
                    ]
                    branch.sort(key=lambda row: float(row["beta_deg"]))
                    for left, right in zip(branch, branch[1:]):
                        max_step_jump = max(max_step_jump, abs(float(right["Lambda"]) - float(left["Lambda"])))
                summaries.append(
                    {
                        "mu": float(mu),
                        "eta": float(eta),
                        "subsystem": subsystem,
                        "n_roots_requested": args.n_roots,
                        "n_beta_points": len(beta_values),
                        "n_missing_roots": missing,
                        "n_suspicious_jumps": suspicious,
                        "min_lambda": min(values) if values else float("nan"),
                        "max_lambda": max(values) if values else float("nan"),
                        "max_step_jump": max_step_jump if values else float("nan"),
                        "warning_count": warning_count,
                        "notes": "ok" if not (missing or suspicious or warning_count) else "inspect warnings/spike audit",
                    }
                )
    return summaries


def build_union_rows(rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[float, float, float], list[dict[str, object]]] = {}
    for row in rows:
        if isfinite(float(row["Lambda"])):
            grouped.setdefault((float(row["beta_deg"]), float(row["mu"]), float(row["eta"])), []).append(row)

    union_rows: list[dict[str, object]] = []
    for (beta_deg, mu, eta), items in grouped.items():
        sorted_items = sorted(items, key=lambda row: (float(row["Lambda"]), str(row["subsystem"]), int(row["sorted_index"])))
        for full_index, item in enumerate(sorted_items, start=1):
            union_rows.append(
                {
                    "beta_deg": beta_deg,
                    "mu": mu,
                    "eta": eta,
                    "full_sorted_index": full_index,
                    "Lambda": float(item["Lambda"]),
                    "source_subsystem": item["subsystem"],
                    "source_sorted_index_within_subsystem": item["sorted_index"],
                }
            )
    union_rows.sort(
        key=lambda row: (float(row["mu"]), float(row["eta"]), float(row["beta_deg"]), int(row["full_sorted_index"]))
    )
    return union_rows


def plotted_series(
    rows: Sequence[dict[str, object]],
    *,
    mu: float,
    eta: float,
    subsystem: str,
    sorted_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    selected = [
        row
        for row in rows
        if float(row["mu"]) == float(mu)
        and float(row["eta"]) == float(eta)
        and str(row["subsystem"]) == subsystem
        and int(row["sorted_index"]) == int(sorted_index)
    ]
    selected.sort(key=lambda row: float(row["beta_deg"]))
    x_values = np.asarray([float(row["beta_deg"]) for row in selected], dtype=float)
    y_values = []
    for row in selected:
        if bool(row["plotted_as_nan"]):
            y_values.append(float("nan"))
        else:
            y_values.append(float(row["Lambda"]))
    return x_values, np.asarray(y_values, dtype=float)


def plot_case(rows: Sequence[dict[str, object]], args: Args, *, mu: float, eta: float) -> Path:
    fig, ax = plt.subplots(figsize=(10.5, 6.4), constrained_layout=True)
    colors = plt.cm.tab20(np.linspace(0.0, 1.0, max(args.n_roots, 2)))
    for index in range(1, args.n_roots + 1):
        color = colors[index - 1]
        for subsystem, linestyle, alpha, linewidth in (
            ("in_plane", "-", 0.9, 1.35),
            ("out_of_plane", "--", 0.75, 1.2),
        ):
            x_values, y_values = plotted_series(rows, mu=mu, eta=eta, subsystem=subsystem, sorted_index=index)
            ax.plot(x_values, y_values, linestyle=linestyle, color=color, alpha=alpha, linewidth=linewidth)
    ax.set_xlabel("beta, degrees")
    ax.set_ylabel("Lambda")
    ax.set_title(
        f"epsilon={args.epsilon:g}, poisson={args.poisson:g}, mu={float(mu):g}, eta={float(eta):g}\n"
        "solid: in-plane; dashed: out-of-plane; sorted frequencies"
    )
    ax.grid(True, alpha=0.25)
    ax.legend(
        handles=[
            Line2D([0], [0], color="black", linestyle="-", label="in-plane sorted roots"),
            Line2D([0], [0], color="black", linestyle="--", label="out-of-plane sorted roots"),
        ],
        loc="best",
        frameon=False,
    )
    path = (
        args.output_dir
        / f"lambda_beta_mu_{_param_label(mu)}_eta_{_param_label(eta)}_in_plane_vs_out_of_plane.png"
    )
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_overview(rows: Sequence[dict[str, object]], args: Args) -> Path:
    n_rows = len(args.eta_values)
    n_cols = len(args.mu_values)
    overview_roots = min(8, args.n_roots)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(4.3 * n_cols, 3.2 * n_rows),
        squeeze=False,
        constrained_layout=True,
    )
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(overview_roots, 2)))
    for row_index, eta in enumerate(args.eta_values):
        for col_index, mu in enumerate(args.mu_values):
            ax = axes[row_index][col_index]
            for index in range(1, overview_roots + 1):
                color = colors[index - 1]
                for subsystem, linestyle, alpha in (("in_plane", "-", 0.9), ("out_of_plane", "--", 0.75)):
                    x_values, y_values = plotted_series(
                        rows,
                        mu=float(mu),
                        eta=float(eta),
                        subsystem=subsystem,
                        sorted_index=index,
                    )
                    ax.plot(x_values, y_values, linestyle=linestyle, color=color, alpha=alpha, linewidth=1.0)
            ax.set_title(f"mu={float(mu):g}, eta={float(eta):g}", fontsize=10)
            ax.grid(True, alpha=0.2)
            if row_index == n_rows - 1:
                ax.set_xlabel("beta, degrees")
            if col_index == 0:
                ax.set_ylabel("Lambda")
    fig.suptitle(
        f"In-plane solid vs out-of-plane dashed sorted Lambda(beta), first {overview_roots} roots; "
        f"epsilon={args.epsilon:g}, poisson={args.poisson:g}",
        fontsize=12,
    )
    path = args.output_dir / "in_plane_vs_out_of_plane_lambda_beta_overview.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def count_retry_fixed(spike_rows: Sequence[dict[str, object]]) -> int:
    return sum(
        1
        for row in spike_rows
        if bool(row["retry_attempted"]) and bool(row["retry_changed_value"]) and not bool(row["plotted_as_nan"])
    )


def mu0_eta0_spike_rows(spike_rows: Sequence[dict[str, object]]) -> list[dict[str, object]]:
    return [
        row
        for row in spike_rows
        if abs(float(row["mu"])) <= 1.0e-15 and abs(float(row["eta"])) <= 1.0e-15
    ]


def coincidence_observation(rows: Sequence[dict[str, object]]) -> str:
    selected = [
        row
        for row in rows
        if abs(float(row["mu"])) <= 1.0e-15
        and abs(float(row["eta"])) <= 1.0e-15
        and abs(float(row["beta_deg"])) <= 1.0e-15
        and int(row["sorted_index"]) <= 8
        and isfinite(float(row["Lambda"]))
    ]
    in_values = [float(row["Lambda"]) for row in selected if row["subsystem"] == "in_plane"]
    out_values = [float(row["Lambda"]) for row in selected if row["subsystem"] == "out_of_plane"]
    if not in_values or not out_values:
        return "mu=0, eta=0 coincidence check unavailable because roots are missing."
    matches = 0
    for value in in_values:
        if any(abs(value - other) <= 1.0e-6 * max(abs(value), 1.0) for other in out_values):
            matches += 1
    return (
        f"At beta=0, mu=0, eta=0, {matches} of the first {len(in_values)} in-plane roots "
        "have an out-of-plane root at the same Lambda within 1e-6 relative tolerance; overlapping solid/dashed curves are expected."
    )


def torsion_observation(rows: Sequence[dict[str, object]]) -> str:
    selected = [
        row
        for row in rows
        if abs(float(row["mu"])) <= 1.0e-15
        and abs(float(row["eta"])) <= 1.0e-15
        and abs(float(row["beta_deg"])) <= 1.0e-15
        and str(row["subsystem"]) == "out_of_plane"
        and isfinite(float(row["Lambda"]))
    ]
    near_torsion = [row for row in selected if 18.0 <= float(row["Lambda"]) <= 21.0]
    if not near_torsion:
        return "The beta=0, mu=0, eta=0 out-of-plane roots did not include a Lambda in the 18-21 torsion-reference window."
    root = min(near_torsion, key=lambda row: abs(float(row["Lambda"]) - 19.7399749487))
    return (
        "For beta=0, mu=0, eta=0, an out-of-plane root in the torsion-reference window "
        f"appears at sorted index {int(root['sorted_index'])}, Lambda={float(root['Lambda']):.8g}."
    )


def beta_effect_observation(rows: Sequence[dict[str, object]], args: Args) -> str:
    selected = [
        row
        for row in rows
        if float(row["mu"]) == 0.3
        and abs(float(row["eta"])) <= 1.0e-15
        and int(row["sorted_index"]) <= min(4, args.n_roots)
        and isfinite(float(row["Lambda"]))
    ]
    if not selected:
        return "Beta-dependence observation unavailable for mu=0.3, eta=0."
    ranges: list[str] = []
    for subsystem in SUBSYSTEMS:
        values = [float(row["Lambda"]) for row in selected if row["subsystem"] == subsystem]
        if values:
            ranges.append(f"{subsystem}: range {max(values) - min(values):.4g} over first low sorted roots")
    return "Beta affects the two sorted subsystems differently in the diagnostic grid (" + "; ".join(ranges) + ")."


def write_report(
    *,
    path: Path,
    args: Args,
    beta_values: np.ndarray,
    rows: Sequence[dict[str, object]],
    summary_rows: Sequence[dict[str, object]],
    spike_rows: Sequence[dict[str, object]],
    generated_files: Sequence[Path],
) -> None:
    warning_count = sum(1 for row in rows if str(row.get("root_solver_warning", "")))
    missing_count = sum(1 for row in rows if not isfinite(float(row["Lambda"])))
    suspicious_count = sum(1 for row in rows if bool(row["suspicious_jump"]))
    plotted_nan_count = sum(1 for row in rows if bool(row["plotted_as_nan"]))
    suspicious_plotted_nan_count = sum(1 for row in spike_rows if bool(row["plotted_as_nan"]))
    retry_fixed = count_retry_fixed(spike_rows)
    mu0_spikes = mu0_eta0_spike_rows(spike_rows)
    total_possible = len(beta_values) * len(args.mu_values) * len(args.eta_values) * len(SUBSYSTEMS) * args.n_roots
    root_counts = []
    for subsystem in SUBSYSTEMS:
        counts = [
            int(row["root_count_found"])
            for row in rows
            if str(row["subsystem"]) == subsystem and int(row["sorted_index"]) == 1
        ]
        if counts:
            root_counts.append(f"{subsystem}: min {min(counts)}, max {max(counts)} roots found per beta point")
    lines = [
        "# In-Plane vs Out-of-Plane Lambda(beta) Diagnostic",
        "",
        "Analytic-only sorted-frequency comparison of the two independent linear subsystems.",
        "",
        "## Parameters",
        "",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        f"- beta grid: {float(beta_values[0]):g}..{float(beta_values[-1]):g} deg, step {args.beta_step:g} deg, {len(beta_values)} points",
        f"- mu values: {', '.join(f'{value:g}' for value in args.mu_values)}",
        f"- eta values: {', '.join(f'{value:g}' for value in args.eta_values)}",
        f"- roots requested per subsystem: {args.n_roots}",
        f"- lambda_max: {args.lambda_max:g}",
        f"- lambda_max_retry: {args.lambda_max_retry:g}",
        f"- root scan step: {args.root_scan_density:g}",
        f"- retry scan step: {args.retry_scan_density:g}",
        f"- jump thresholds: abs>{args.jump_abs_threshold:g} or rel>{args.jump_rel_threshold:g}",
        f"- smoke mode: {'yes' if args.smoke else 'no'}",
        "",
        "## Root Counts And Warnings",
        "",
        f"- expected CSV rows: {total_possible}",
        f"- missing root rows: {missing_count}",
        f"- root-finder warning rows: {warning_count}",
        f"- root-count summary: {'; '.join(root_counts) if root_counts else 'unavailable'}",
        "",
        "## Spike Audit",
        "",
        f"- suspicious jumps found: {suspicious_count}",
        f"- suspicious jump audit rows: {len(spike_rows)}",
        f"- fixed by retry: {retry_fixed}",
        f"- suspicious jumps hidden by NaN segmentation: {suspicious_plotted_nan_count}",
        f"- missing-root rows also plotted as NaN: {missing_count}",
        f"- total plotted-as-NaN rows: {plotted_nan_count}",
        f"- mu=0, eta=0 suspicious jumps: {len(mu0_spikes)}",
        f"- mu=0, eta=0 plotted as NaN: {sum(1 for row in mu0_spikes if bool(row['plotted_as_nan']))}",
        "",
        "## Qualitative Observations",
        "",
        f"- {coincidence_observation(rows)}",
        f"- {torsion_observation(rows)}",
        f"- {beta_effect_observation(rows, args)}",
        "",
        "## Interpretation Limits",
        "",
        "- These are sorted-frequency diagnostic plots.",
        "- No descendant tracking was used.",
        "- No crossing or no-crossing claim is made.",
        "- No FEM, 1D FEM diagnostic, Gmsh, CalculiX, or 3D FEM workflow was run.",
        "- In-plane curves are plotted as solid lines; out-of-plane curves are plotted as dashed lines.",
        "- Duplicate or overlapping subsystem roots are not removed or perturbed.",
        "",
        "## Case Summary",
        "",
    ]
    for summary in summary_rows:
        lines.append(
            f"- mu={float(summary['mu']):g}, eta={float(summary['eta']):g}, {summary['subsystem']}: "
            f"missing={summary['n_missing_roots']}, suspicious={summary['n_suspicious_jumps']}, "
            f"warnings={summary['warning_count']}, Lambda range=[{float(summary['min_lambda']):.6g}, {float(summary['max_lambda']):.6g}]"
        )
    lines.extend(["", "## Generated Files", ""])
    for item in generated_files:
        lines.append(f"- `{rel(item)}`")
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def run(args: Args) -> dict[str, object]:
    prepare_output_dir(args)
    beta_values = beta_grid(args.beta_min, args.beta_max, args.beta_step)
    root_map = compute_root_map(args, beta_values)
    rows, spike_rows = build_rows_and_spike_audit(args, beta_values, root_map)
    summary_rows = build_case_summary(rows, args, beta_values)
    union_rows = build_union_rows(rows)

    sorted_path = args.output_dir / "in_plane_vs_out_of_plane_lambda_beta_sorted_roots.csv"
    summary_path = args.output_dir / "in_plane_vs_out_of_plane_lambda_beta_case_summary.csv"
    spike_path = args.output_dir / "in_plane_vs_out_of_plane_lambda_beta_spike_audit.csv"
    union_path = args.output_dir / "full_analytic_union_lambda_beta.csv"
    write_csv(
        sorted_path,
        rows,
        [
            "beta_deg",
            "beta_rad",
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "subsystem",
            "sorted_index",
            "Lambda",
            "root_solver_warning",
            "root_count_found",
            "root_count_requested",
            "lambda_max_used",
            "suspicious_jump",
            "plotted_as_nan",
            "notes",
        ],
    )
    write_csv(
        summary_path,
        summary_rows,
        [
            "mu",
            "eta",
            "subsystem",
            "n_roots_requested",
            "n_beta_points",
            "n_missing_roots",
            "n_suspicious_jumps",
            "min_lambda",
            "max_lambda",
            "max_step_jump",
            "warning_count",
            "notes",
        ],
    )
    write_csv(
        spike_path,
        spike_rows,
        [
            "mu",
            "eta",
            "subsystem",
            "sorted_index",
            "beta_deg",
            "Lambda_raw",
            "Lambda_prev",
            "Lambda_next",
            "jump_prev",
            "jump_next",
            "retry_attempted",
            "retry_changed_value",
            "plotted_as_nan",
            "notes",
        ],
    )
    write_csv(
        union_path,
        union_rows,
        [
            "beta_deg",
            "mu",
            "eta",
            "full_sorted_index",
            "Lambda",
            "source_subsystem",
            "source_sorted_index_within_subsystem",
        ],
    )

    figure_paths: list[Path] = []
    for mu in args.mu_values:
        for eta in args.eta_values:
            figure_paths.append(plot_case(rows, args, mu=float(mu), eta=float(eta)))
    overview_path = plot_overview(rows, args)
    report_path = args.output_dir / "in_plane_vs_out_of_plane_lambda_beta_report.md"
    generated_files = [sorted_path, summary_path, spike_path, union_path, *figure_paths, overview_path, report_path]
    write_report(
        path=report_path,
        args=args,
        beta_values=beta_values,
        rows=rows,
        summary_rows=summary_rows,
        spike_rows=spike_rows,
        generated_files=generated_files,
    )
    print(f"saved sorted roots CSV: {sorted_path}")
    print(f"saved case summary CSV: {summary_path}")
    print(f"saved spike audit CSV: {spike_path}")
    print(f"saved union CSV: {union_path}")
    print(f"saved {len(figure_paths)} per-case PNGs")
    print(f"saved overview PNG: {overview_path}")
    print(f"saved report: {report_path}")
    print(f"beta grid: {float(beta_values[0]):g}..{float(beta_values[-1]):g} step {args.beta_step:g} ({len(beta_values)} points)")
    print(f"mu values: {', '.join(f'{value:g}' for value in args.mu_values)}")
    print(f"eta values: {', '.join(f'{value:g}' for value in args.eta_values)}")
    print(f"roots requested: {args.n_roots}")
    print(f"warning rows: {sum(1 for row in rows if str(row.get('root_solver_warning', '')))}")
    print(f"suspicious jumps: {sum(1 for row in rows if bool(row['suspicious_jump']))}")
    return {
        "sorted_path": sorted_path,
        "summary_path": summary_path,
        "spike_path": spike_path,
        "union_path": union_path,
        "figure_paths": figure_paths,
        "overview_path": overview_path,
        "report_path": report_path,
    }


def main(argv: list[str] | None = None) -> int:
    run(parse_args(argv))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
