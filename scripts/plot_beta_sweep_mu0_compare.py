from __future__ import annotations

import csv
from contextlib import contextmanager
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import linear_sum_assignment


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import BeamParams, lambdas_to_frequencies  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots, find_roots_scan_bisect  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.sweep_grid_policy import (  # noqa: E402
    LOCAL_BETA_REFINEMENT_STEP,
    PRESENTATION_BETA_STEP,
    nominal_step,
    presentation_beta_grid as build_presentation_beta_grid,
)


MU_VALUE = 0.0
RADIUS_VALUES = (0.005, 0.01)
N_BRANCHES = 10
N_ANALYTIC_CANDIDATES = 16
N_FEM_SOLVE = 16
ANALYTIC_LMIN = 0.2
ANALYTIC_LMAX0 = 45.0
ANALYTIC_SCAN_STEP = 0.02
ANALYTIC_FINE_SCAN_STEP = 0.0015
ANALYTIC_CONTINUATION_SCAN_STEP = 0.0002
ANALYTIC_GROW_FACTOR = 1.35
ANALYTIC_MAX_TRIES = 8
LOCAL_CONTINUATION_WINDOW = 0.08

PROBLEM_WINDOWS_BY_RADIUS = {
    0.005: ((28.6, 29.3),),
    0.01: ((43.5, 45.5), (75.5, 77.5)),
    0.015: ((22.5, 23.5), (77.5, 78.5)),
    0.02: ((29.5, 31.5), (87.5, 88.5)),
}

BASE_E = 2.1e11
BASE_RHO = 7800.0
BASE_L_TOTAL = 2.0

FEM_STATE_KEYS = (
    "r",
    "E",
    "rho",
    "L_tot",
    "ell",
    "A",
    "I",
    "EI",
    "EA",
    "rhoA",
    "eps",
    "scale",
    "EI_nd",
    "rhoA_nd",
    "EA_nd",
)

RESULTS_DIR = REPO_ROOT / "results"


def build_params(radius: float) -> BeamParams:
    return BeamParams(E=BASE_E, rho=BASE_RHO, r=radius, L_total=BASE_L_TOTAL)


def radius_key(radius: float) -> float:
    return round(float(radius), 6)


def radius_tag(radius: float) -> str:
    return f"r{int(round(radius * 1000.0)):03d}"


def analytic_plot_path(radius: float) -> Path:
    return RESULTS_DIR / f"beta_sweep_mu0_{radius_tag(radius)}_compare_smoothed.png"


def analytic_csv_path(radius: float) -> Path:
    return RESULTS_DIR / f"beta_sweep_mu0_{radius_tag(radius)}_compare_smoothed.csv"


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def relative_error(analytic_hz: float, fem_hz: float) -> float:
    return abs(analytic_hz - fem_hz) / abs(fem_hz) if fem_hz else np.nan


def problem_windows(radius: float) -> tuple[tuple[float, float], ...]:
    return PROBLEM_WINDOWS_BY_RADIUS.get(radius_key(radius), ())


def beta_in_problem_window(radius: float, beta_deg: float) -> bool:
    beta = float(beta_deg)
    return any(left <= beta <= right for left, right in problem_windows(radius))


def presentation_beta_grid(radius: float) -> np.ndarray:
    return build_presentation_beta_grid(local_windows=problem_windows(radius))


def analytic_scan_step(radius: float, beta_deg: float) -> float:
    return ANALYTIC_FINE_SCAN_STEP if beta_in_problem_window(radius, beta_deg) else ANALYTIC_SCAN_STEP


def analytic_global_roots(params: BeamParams, beta_deg: float) -> np.ndarray:
    roots = find_first_n_roots(
        beta=np.deg2rad(float(beta_deg)),
        mu=MU_VALUE,
        eps=params.eps,
        n_roots=N_ANALYTIC_CANDIDATES,
        Lmin=ANALYTIC_LMIN,
        Lmax0=ANALYTIC_LMAX0,
        scan_step=analytic_scan_step(params.r, float(beta_deg)),
        grow_factor=ANALYTIC_GROW_FACTOR,
        max_tries=ANALYTIC_MAX_TRIES,
    )
    roots = np.asarray(roots, dtype=float)
    return np.sort(roots[np.isfinite(roots)])


def unique_sorted_roots(candidates: list[float] | np.ndarray) -> np.ndarray:
    values = np.asarray(candidates, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return values
    values = np.sort(values)
    merged = [float(values[0])]
    for value in values[1:]:
        if abs(float(value) - merged[-1]) > 1e-8:
            merged.append(float(value))
    return np.asarray(merged, dtype=float)


def local_continuation_roots(
    params: BeamParams,
    beta_deg: float,
    previous_lambdas: np.ndarray,
) -> np.ndarray:
    candidates: list[float] = []
    for lambda_prev in previous_lambdas:
        if not np.isfinite(lambda_prev):
            continue
        roots_local = find_roots_scan_bisect(
            beta=np.deg2rad(float(beta_deg)),
            mu=MU_VALUE,
            eps=params.eps,
            n_roots=6,
            Lmin=max(ANALYTIC_LMIN, float(lambda_prev) - LOCAL_CONTINUATION_WINDOW),
            Lmax=float(lambda_prev) + LOCAL_CONTINUATION_WINDOW,
            scan_step=ANALYTIC_CONTINUATION_SCAN_STEP,
            bisect_iters=100,
        )
        candidates.extend(float(root) for root in roots_local)
    return unique_sorted_roots(candidates)


def assign_by_frequency(previous_lambdas: np.ndarray, candidate_lambdas: np.ndarray) -> np.ndarray:
    cost = np.abs(previous_lambdas[:, None] - candidate_lambdas[None, :])
    rows, cols = linear_sum_assignment(cost)
    if len(rows) != len(previous_lambdas):
        raise RuntimeError("Failed to assign all analytic branches.")
    assigned = np.full(len(previous_lambdas), np.nan, dtype=float)
    assigned[rows] = candidate_lambdas[cols]
    if np.any(~np.isfinite(assigned)):
        raise RuntimeError("Some analytic branches remained unassigned.")
    return assigned


def _apply_fem_params(params: BeamParams) -> None:
    fem.r = params.r
    fem.E = params.E
    fem.rho = params.rho
    fem.L_tot = params.L_total
    fem.ell = params.L_total / 2.0
    fem.A = np.pi * params.r**2
    fem.I = np.pi * params.r**4 / 4.0
    fem.EI = fem.E * fem.I
    fem.EA = fem.E * fem.A
    fem.rhoA = fem.rho * fem.A
    fem.eps = np.sqrt(fem.I / fem.A) / fem.ell
    fem.scale = np.sqrt(fem.EI / fem.rhoA) / (2.0 * np.pi * fem.ell**2)
    fem.EI_nd = 1.0
    fem.rhoA_nd = 1.0
    fem.EA_nd = 1.0 / fem.eps**2


@contextmanager
def fem_parameter_override(params: BeamParams):
    saved = {name: getattr(fem, name) for name in FEM_STATE_KEYS}
    try:
        _apply_fem_params(params)
        yield
    finally:
        for name, value in saved.items():
            setattr(fem, name, value)


def solve_fem_modes(params: BeamParams, beta_deg: float, n_modes: int) -> tuple[np.ndarray, np.ndarray]:
    with fem_parameter_override(params):
        omega, _ = fem.fem_solve(mu=MU_VALUE, beta_deg=beta_deg, n_modes=n_modes)
        lambda_vals = np.sqrt(np.maximum(omega, 0.0))
        freq_hz = omega * fem.scale
        return lambda_vals, freq_hz


def resolve_beta_grid(params: BeamParams, beta_values: np.ndarray) -> tuple[np.ndarray, float]:
    requested_beta_max = float(beta_values[-1])
    try:
        find_first_n_roots(
            beta=np.deg2rad(requested_beta_max),
            mu=MU_VALUE,
            eps=params.eps,
            n_roots=N_ANALYTIC_CANDIDATES,
            Lmin=ANALYTIC_LMIN,
            Lmax0=ANALYTIC_LMAX0,
            scan_step=analytic_scan_step(params.r, requested_beta_max),
            grow_factor=ANALYTIC_GROW_FACTOR,
            max_tries=ANALYTIC_MAX_TRIES,
        )
        solve_fem_modes(params=params, beta_deg=requested_beta_max, n_modes=N_BRANCHES)
        return beta_values, requested_beta_max
    except Exception:
        fallback_beta = requested_beta_max - 0.1
        beta_adjusted = beta_values.copy()
        beta_adjusted[-1] = fallback_beta
        return beta_adjusted, fallback_beta


def compute_analytic_branches(params: BeamParams, beta_values: np.ndarray, n_branches: int) -> tuple[np.ndarray, np.ndarray]:
    lambdas_tracked = np.full((N_ANALYTIC_CANDIDATES, len(beta_values)), np.nan, dtype=float)
    initial_roots = analytic_global_roots(params=params, beta_deg=float(beta_values[0]))
    if len(initial_roots) < N_ANALYTIC_CANDIDATES:
        raise RuntimeError("Not enough analytic candidate roots at beta=0 for presentation sweep.")
    lambdas_tracked[:, 0] = initial_roots[:N_ANALYTIC_CANDIDATES]

    for col, beta_deg in enumerate(beta_values[1:], start=1):
        previous_lambdas = lambdas_tracked[:, col - 1]
        candidate_roots = analytic_global_roots(params=params, beta_deg=float(beta_deg))
        if beta_in_problem_window(params.r, float(beta_deg)) or beta_in_problem_window(params.r, float(beta_values[col - 1])):
            local_roots = local_continuation_roots(
                params=params,
                beta_deg=float(beta_deg),
                previous_lambdas=previous_lambdas,
            )
            candidate_roots = unique_sorted_roots(np.concatenate((candidate_roots, local_roots)))
        if len(candidate_roots) < N_ANALYTIC_CANDIDATES:
            raise RuntimeError(f"Not enough analytic candidate roots at beta={beta_deg}")
        lambdas_tracked[:, col] = assign_by_frequency(
            previous_lambdas=previous_lambdas,
            candidate_lambdas=candidate_roots,
        )

    lambdas_presentation = lambdas_tracked[:n_branches]
    freq_hz_tracked = lambdas_to_frequencies(lambdas_presentation, params)
    return lambdas_presentation, freq_hz_tracked


def compute_fem_branches(
    params: BeamParams,
    beta_values: np.ndarray,
    analytic_lambda: np.ndarray,
    n_modes_solve: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_branches = analytic_lambda.shape[0]
    lambda_matched = np.full((n_branches, len(beta_values)), np.nan, dtype=float)
    freq_hz_matched = np.full((n_branches, len(beta_values)), np.nan, dtype=float)
    fem_mode_ids = np.full((n_branches, len(beta_values)), -1, dtype=int)

    for col, beta_deg in enumerate(beta_values):
        cur_lambda, cur_freq_hz = solve_fem_modes(params=params, beta_deg=float(beta_deg), n_modes=n_modes_solve)
        cost = np.abs(analytic_lambda[:, col][:, None] - cur_lambda[None, :])
        rows, cols = linear_sum_assignment(cost)
        if len(rows) != n_branches:
            raise RuntimeError(f"Failed to match all FEM modes at beta={beta_deg}")
        lambda_matched[rows, col] = cur_lambda[cols]
        freq_hz_matched[rows, col] = cur_freq_hz[cols]
        fem_mode_ids[rows, col] = cols + 1

    return lambda_matched, freq_hz_matched, fem_mode_ids


def build_rows(
    beta_values: np.ndarray,
    analytic_lambda: np.ndarray,
    analytic_hz: np.ndarray,
    fem_lambda: np.ndarray,
    fem_hz: np.ndarray,
    fem_mode_ids: np.ndarray,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for branch_pos in range(analytic_lambda.shape[0]):
        branch_id = f"branch_{branch_pos + 1:02d}"
        for col, beta_deg in enumerate(beta_values):
            rows.append(
                {
                    "beta": float(beta_deg),
                    "branch_id": branch_id,
                    "fem_mode_id": int(fem_mode_ids[branch_pos, col]),
                    "analytic_lambda": float(analytic_lambda[branch_pos, col]),
                    "fem_lambda": float(fem_lambda[branch_pos, col]),
                    "analytic_hz": float(analytic_hz[branch_pos, col]),
                    "fem_hz": float(fem_hz[branch_pos, col]),
                    "relative_error": float(relative_error(analytic_hz[branch_pos, col], fem_hz[branch_pos, col])),
                }
            )
    return rows


def plot_comparison(
    params: BeamParams,
    beta_values: np.ndarray,
    analytic_lambda: np.ndarray,
    fem_lambda: np.ndarray,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(11.2, 6.5))
    for branch_pos in range(analytic_lambda.shape[0]):
        (line,) = ax.plot(beta_values, analytic_lambda[branch_pos], linewidth=2.0)
        color = line.get_color()
        ax.plot(
            beta_values,
            fem_lambda[branch_pos],
            linestyle="None",
            marker="o",
            markersize=3.5,
            color=color,
            markerfacecolor=color,
            markeredgecolor=color,
        )

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="solid — analytic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.0, label="markers — FEM"),
    ]
    ax.legend(handles=style_handles, fontsize=9, ncols=2, loc="upper left")
    ax.set_xlabel("beta, deg")
    ax.set_ylabel("Безразмерный частотный параметр Λ")
    ax.set_title(f"Analytic vs FEM, mu=0, r={params.r:.3f} m, first {analytic_lambda.shape[0]} branches")
    ax.set_xlim(beta_values[0], beta_values[-1])
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def run_radius_case(radius: float) -> dict[str, object]:
    params = build_params(radius)
    beta_values_requested = presentation_beta_grid(radius)
    beta_values, beta_endpoint_used = resolve_beta_grid(params=params, beta_values=beta_values_requested)
    analytic_lambda, analytic_hz = compute_analytic_branches(params=params, beta_values=beta_values, n_branches=N_BRANCHES)
    fem_lambda, fem_hz, fem_mode_ids = compute_fem_branches(
        params=params,
        beta_values=beta_values,
        analytic_lambda=analytic_lambda,
        n_modes_solve=N_FEM_SOLVE,
    )

    rows = build_rows(
        beta_values=beta_values,
        analytic_lambda=analytic_lambda,
        analytic_hz=analytic_hz,
        fem_lambda=fem_lambda,
        fem_hz=fem_hz,
        fem_mode_ids=fem_mode_ids,
    )
    csv_path = analytic_csv_path(radius)
    plot_path = analytic_plot_path(radius)
    write_csv_rows(csv_path, rows)
    plot_comparison(
        params=params,
        beta_values=beta_values,
        analytic_lambda=analytic_lambda,
        fem_lambda=fem_lambda,
        output_path=plot_path,
    )

    max_rel_error = max(float(row["relative_error"]) for row in rows)
    return {
        "radius": radius,
        "params": params,
        "beta_values": beta_values,
        "beta_count": len(beta_values),
        "beta_endpoint_used": beta_endpoint_used,
        "csv_path": csv_path,
        "plot_path": plot_path,
        "rows": rows,
        "max_rel_error": max_rel_error,
        "shown_branches": [f"branch_{idx + 1:02d}" for idx in range(N_BRANCHES)],
        "problem_windows": problem_windows(radius),
        "analytic_candidates": N_ANALYTIC_CANDIDATES,
        "fine_scan_step": ANALYTIC_FINE_SCAN_STEP,
        "continuation_scan_step": ANALYTIC_CONTINUATION_SCAN_STEP,
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    summaries = [run_radius_case(radius) for radius in RADIUS_VALUES]

    for summary in summaries:
        print(f"radius={summary['radius']:.3f} m")
        print(f"  beta endpoint used: {summary['beta_endpoint_used']:.1f} deg")
        print(f"  beta points: {summary['beta_count']}")
        print(f"  presentation beta base step: {PRESENTATION_BETA_STEP:.1f} deg")
        print(f"  effective min beta step: {nominal_step(summary['beta_values']):.3f} deg")
        print(f"  problem windows: {summary['problem_windows']}")
        print(f"  analytic candidates: {summary['analytic_candidates']}")
        print(f"  local beta refinement step: {LOCAL_BETA_REFINEMENT_STEP:.1f} deg")
        print(f"  local fine scan_step: {summary['fine_scan_step']}")
        print(f"  local continuation scan_step: {summary['continuation_scan_step']}")
        print(f"  shown branches: {', '.join(summary['shown_branches'])}")
        print(f"  saved plot: {summary['plot_path']}")
        print(f"  saved csv: {summary['csv_path']}")
        print(f"  max relative error: {summary['max_rel_error']:.6e}")


if __name__ == "__main__":
    main()
