from __future__ import annotations

import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path
import os
import sys
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

DEFAULT_EPSILON = 0.0025
DEFAULT_ROOT_SCAN_STEP = 0.02
DEFAULT_ROOT_LMAX0 = 35.0

LAMBDA_ETA_STEP = 0.01
LAMBDA_ETA_MIN = -0.5
LAMBDA_ETA_MAX = 0.5

SLICE_BETA_STEP = 0.25
SLICE_MU_VALUES = (0.1, 0.2, 0.3)
SLICE_ETA_VALUES = (-0.3, -0.2, -0.1, 0.1, 0.2, 0.3)

HEATMAP_ETA_STEP = 0.025
HEATMAP_MU_STEP = 0.05
HEATMAP_BETA_STEP = 1.0
HEATMAP_ROOTS = 10

DISPLAY_ROOTS = 6
GAP_ROOTS = 7
RELATIVE_DENOMINATOR_TOL = 1.0e-12

LAMBDA_ETA_PNG = RESULTS_DIR / "diagnostic_lambda_eta_beta0_mu0_eps0p0025.png"
LAMBDA_ETA_CSV = RESULTS_DIR / "diagnostic_lambda_eta_beta0_mu0_eps0p0025.csv"
LAMBDA_BETA_CSV = RESULTS_DIR / "diagnostic_lambda_beta_eps0p0025_mu_eta_slices.csv"
LAMBDA_BETA_OVERVIEW_PNG = RESULTS_DIR / "diagnostic_lambda_beta_eps0p0025_mu_eta_grid_overview.png"
HEATMAP_CSV = RESULTS_DIR / "diagnostic_eta_mu_beta_heatmap_metrics_eps0p0025.csv"
HEATMAP_REPORT = RESULTS_DIR / "diagnostic_eta_mu_beta_heatmap_metrics_eps0p0025_report.md"

HEATMAP_G_MIN_PNG = RESULTS_DIR / "diagnostic_heatmap_g_min_eta_mu_eps0p0025.png"
HEATMAP_BETA_G_MIN_PNG = RESULTS_DIR / "diagnostic_heatmap_beta_at_g_min_eta_mu_eps0p0025.png"
HEATMAP_PAIR_G_MIN_PNG = RESULTS_DIR / "diagnostic_heatmap_pair_at_g_min_eta_mu_eps0p0025.png"
HEATMAP_S_MEAN_PNG = RESULTS_DIR / "diagnostic_heatmap_S_mean_eta_mu_eps0p0025.png"
HEATMAP_S_MAX_PNG = RESULTS_DIR / "diagnostic_heatmap_S_max_eta_mu_eps0p0025.png"
HEATMAP_S_MEAN_REL_PNG = RESULTS_DIR / "diagnostic_heatmap_S_mean_rel_eta_mu_eps0p0025.png"

LAMBDA_ETA_FIELDNAMES = [
    "eta",
    "beta_deg",
    "mu",
    "epsilon",
    "sorted_index",
    "Lambda",
]

LAMBDA_BETA_FIELDNAMES = [
    "mu",
    "eta",
    "epsilon",
    "beta_deg",
    "sorted_index",
    "Lambda",
    "adjacent_gap_to_next_if_available",
]

HEATMAP_FIELDNAMES = [
    "eta",
    "mu",
    "epsilon",
    "g_min",
    "beta_at_g_min",
    "pair_at_g_min",
    "lambda_lower_at_g_min",
    "lambda_upper_at_g_min",
    "S_mean",
    "S_max",
    "S_mean_rel",
    "S_branch1",
    "S_branch2",
    "S_branch3",
    "S_branch4",
    "S_branch5",
    "S_branch6",
    "S_rel_branch1",
    "S_rel_branch2",
    "S_rel_branch3",
    "S_rel_branch4",
    "S_rel_branch5",
    "S_rel_branch6",
    "solver_warning_count",
    "missing_root_count",
    "notes",
]


@dataclass(frozen=True)
class RootSolve:
    roots: np.ndarray
    warning_count: int
    missing_root_count: int
    note: str


def number_token(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def number_text(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    value_f = float(value)
    if not np.isfinite(value_f):
        return ""
    return f"{value_f:.12g}"


def pair_label(pair_index: int) -> str:
    return f"{int(pair_index)}-{int(pair_index) + 1}"


def inclusive_grid(start: float, stop: float, step: float) -> np.ndarray:
    start_f = float(start)
    stop_f = float(stop)
    step_f = float(step)
    if step_f <= 0.0:
        raise ValueError("grid step must be positive")
    values = np.arange(start_f, stop_f + 0.5 * step_f, step_f, dtype=float)
    if values.size == 0:
        values = np.array([start_f, stop_f], dtype=float)
    if abs(float(values[0]) - start_f) > 1.0e-12:
        values = np.insert(values, 0, start_f)
    if abs(float(values[-1]) - stop_f) > 1.0e-12:
        values = np.append(values, stop_f)
    values[0] = start_f
    values[-1] = stop_f
    return np.unique(np.round(values, 12))


def solve_roots(
    *,
    beta_deg: float,
    mu: float,
    eta: float,
    epsilon: float,
    n_roots: int,
    root_scan_step: float,
    root_lmax0: float,
) -> RootSolve:
    notes: list[str] = []
    warning_count = 0
    try:
        roots = find_first_n_roots_eta(
            float(np.deg2rad(beta_deg)),
            float(mu),
            float(epsilon),
            float(eta),
            int(n_roots),
            Lmax0=float(root_lmax0),
            scan_step=float(root_scan_step),
        )
    except Exception as exc:  # pragma: no cover - diagnostic failure path
        roots = np.full(int(n_roots), np.nan, dtype=float)
        warning_count += 1
        notes.append(f"root solve failed at beta={float(beta_deg):g}: {exc}")

    roots = np.asarray(roots, dtype=float)
    if roots.shape[0] < int(n_roots):
        padded = np.full(int(n_roots), np.nan, dtype=float)
        padded[: roots.shape[0]] = roots
        roots = padded
    missing = int(np.count_nonzero(~np.isfinite(roots[: int(n_roots)])))
    if missing:
        warning_count += 1
        notes.append(f"missing {missing} roots at beta={float(beta_deg):g}")
    return RootSolve(roots=roots[: int(n_roots)], warning_count=warning_count, missing_root_count=missing, note="; ".join(notes))


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def roots_grid_for_parameter_scan(
    *,
    beta_values: Sequence[float],
    mu: float,
    eta: float,
    epsilon: float,
    n_roots: int,
    root_scan_step: float,
    root_lmax0: float,
) -> tuple[np.ndarray, int, int, list[str]]:
    beta_array = np.asarray(beta_values, dtype=float)
    roots = np.full((len(beta_array), int(n_roots)), np.nan, dtype=float)
    warning_count = 0
    missing_count = 0
    notes: list[str] = []
    for idx, beta_deg in enumerate(beta_array):
        solved = solve_roots(
            beta_deg=float(beta_deg),
            mu=float(mu),
            eta=float(eta),
            epsilon=float(epsilon),
            n_roots=int(n_roots),
            root_scan_step=float(root_scan_step),
            root_lmax0=float(root_lmax0),
        )
        roots[idx] = solved.roots
        warning_count += solved.warning_count
        missing_count += solved.missing_root_count
        if solved.note:
            notes.append(solved.note)
    return roots, warning_count, missing_count, notes


def compute_lambda_eta(
    *,
    epsilon: float,
    eta_values: np.ndarray,
    root_scan_step: float,
    root_lmax0: float,
) -> tuple[np.ndarray, list[dict[str, object]], int, int, list[str]]:
    roots = np.full((len(eta_values), DISPLAY_ROOTS), np.nan, dtype=float)
    rows: list[dict[str, object]] = []
    warning_count = 0
    missing_count = 0
    notes: list[str] = []
    for eta_idx, eta in enumerate(eta_values):
        solved = solve_roots(
            beta_deg=0.0,
            mu=0.0,
            eta=float(eta),
            epsilon=float(epsilon),
            n_roots=DISPLAY_ROOTS,
            root_scan_step=float(root_scan_step),
            root_lmax0=float(root_lmax0),
        )
        roots[eta_idx] = solved.roots[:DISPLAY_ROOTS]
        warning_count += solved.warning_count
        missing_count += solved.missing_root_count
        if solved.note:
            notes.append(solved.note)
        for sorted_idx, value in enumerate(solved.roots[:DISPLAY_ROOTS], start=1):
            rows.append(
                {
                    "eta": number_text(float(eta)),
                    "beta_deg": number_text(0.0),
                    "mu": number_text(0.0),
                    "epsilon": number_text(float(epsilon)),
                    "sorted_index": sorted_idx,
                    "Lambda": number_text(float(value)),
                }
            )
    return roots, rows, warning_count, missing_count, notes


def plot_lambda_eta(eta_values: np.ndarray, roots: np.ndarray, output: Path, *, epsilon: float) -> None:
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    for root_idx in range(DISPLAY_ROOTS):
        ax.plot(
            eta_values,
            roots[:, root_idx],
            lw=1.45,
            color=colors[root_idx % len(colors)],
            label=f"sorted {root_idx + 1}",
        )
    ax.axvline(0.0, color="0.68", lw=0.8, ls=":")
    ax.set_xlabel("eta")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Diagnostic sorted Lambda(eta), beta=0 deg, mu=0, epsilon={epsilon:g}")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def lambda_beta_png_path(mu: float, eta: float) -> Path:
    return RESULTS_DIR / f"diagnostic_lambda_beta_eps0p0025_mu{number_token(mu)}_eta_{number_token(eta)}.png"


def compute_lambda_beta_slices(
    *,
    epsilon: float,
    beta_values: np.ndarray,
    root_scan_step: float,
    root_lmax0: float,
) -> tuple[dict[tuple[float, float], np.ndarray], list[dict[str, object]], int, int, list[str]]:
    rows: list[dict[str, object]] = []
    grids: dict[tuple[float, float], np.ndarray] = {}
    warning_count = 0
    missing_count = 0
    notes: list[str] = []
    for mu in SLICE_MU_VALUES:
        for eta in SLICE_ETA_VALUES:
            roots, warnings, missing, local_notes = roots_grid_for_parameter_scan(
                beta_values=beta_values,
                mu=float(mu),
                eta=float(eta),
                epsilon=float(epsilon),
                n_roots=GAP_ROOTS,
                root_scan_step=float(root_scan_step),
                root_lmax0=float(root_lmax0),
            )
            grids[(float(mu), float(eta))] = roots
            warning_count += warnings
            missing_count += missing
            notes.extend(local_notes)
            for beta_idx, beta_deg in enumerate(beta_values):
                for sorted_idx in range(1, DISPLAY_ROOTS + 1):
                    value = float(roots[beta_idx, sorted_idx - 1])
                    next_value = float(roots[beta_idx, sorted_idx])
                    gap = next_value - value if np.isfinite(value) and np.isfinite(next_value) else float("nan")
                    rows.append(
                        {
                            "mu": number_text(mu),
                            "eta": number_text(eta),
                            "epsilon": number_text(epsilon),
                            "beta_deg": number_text(float(beta_deg)),
                            "sorted_index": sorted_idx,
                            "Lambda": number_text(value),
                            "adjacent_gap_to_next_if_available": number_text(gap),
                        }
                    )
    return grids, rows, warning_count, missing_count, notes


def plot_lambda_beta_case(
    beta_values: np.ndarray,
    roots: np.ndarray,
    output: Path,
    *,
    epsilon: float,
    mu: float,
    eta: float,
) -> None:
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(figsize=(9.4, 5.4))
    for root_idx in range(DISPLAY_ROOTS):
        ax.plot(
            beta_values,
            roots[:, root_idx],
            lw=1.35,
            color=colors[root_idx % len(colors)],
            label=f"sorted {root_idx + 1}",
        )
    gaps = roots[:, 1:GAP_ROOTS] - roots[:, : GAP_ROOTS - 1]
    if np.any(np.isfinite(gaps)):
        min_gap = float(np.nanmin(gaps))
        if min_gap < 1.0e-3:
            beta_idx, gap_idx = np.unravel_index(int(np.nanargmin(gaps)), gaps.shape)
            ax.text(
                0.02,
                0.02,
                f"candidate close approach: sorted {gap_idx + 1}-{gap_idx + 2}, "
                f"gap {min_gap:.3g} at beta={float(beta_values[beta_idx]):g} deg",
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=8,
                color="0.24",
            )
    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Diagnostic sorted Lambda(beta), mu={mu:g}, eta={eta:g}, epsilon={epsilon:g}")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", ncol=2, fontsize=8, frameon=False)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_lambda_beta_overview(
    beta_values: np.ndarray,
    grids: dict[tuple[float, float], np.ndarray],
    output: Path,
    *,
    epsilon: float,
) -> None:
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    nrows = len(SLICE_ETA_VALUES)
    ncols = len(SLICE_MU_VALUES)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(12.2, 2.2 * nrows),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes_array = np.asarray(axes)
    for row_idx, eta in enumerate(SLICE_ETA_VALUES):
        for col_idx, mu in enumerate(SLICE_MU_VALUES):
            ax = axes_array[row_idx, col_idx]
            roots = grids[(float(mu), float(eta))]
            for root_idx in range(DISPLAY_ROOTS):
                ax.plot(beta_values, roots[:, root_idx], lw=0.95, color=colors[root_idx % len(colors)])
            ax.set_title(f"mu={mu:g}, eta={eta:g}", fontsize=9)
            ax.grid(True, color="0.9", linewidth=0.5)
            if row_idx == nrows - 1:
                ax.set_xlabel("beta (deg)")
            if col_idx == 0:
                ax.set_ylabel("Lambda")
    handles = [
        plt.Line2D([0], [0], color=colors[idx % len(colors)], lw=1.4, label=f"sorted {idx + 1}")
        for idx in range(DISPLAY_ROOTS)
    ]
    fig.legend(handles=handles, loc="upper center", ncol=DISPLAY_ROOTS, frameon=False, fontsize=8)
    fig.suptitle(f"Diagnostic sorted Lambda(beta) slices, epsilon={epsilon:g}", fontsize=12)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def finite_range(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return float("nan"), float("nan")
    return float(np.min(finite)), float(np.max(finite))


def compute_heatmap_metric(task: tuple[float, float, tuple[float, ...], float, int, float, float]) -> dict[str, object]:
    eta, mu, beta_values_tuple, epsilon, n_roots, root_scan_step, root_lmax0 = task
    beta_values = np.asarray(beta_values_tuple, dtype=float)
    roots, solver_warning_count, missing_root_count, solve_notes = roots_grid_for_parameter_scan(
        beta_values=beta_values,
        mu=float(mu),
        eta=float(eta),
        epsilon=float(epsilon),
        n_roots=int(n_roots),
        root_scan_step=float(root_scan_step),
        root_lmax0=float(root_lmax0),
    )

    notes: list[str] = []
    if solve_notes:
        notes.append(f"root solve notes: {len(solve_notes)}")

    gap_matrix = roots[:, 1:GAP_ROOTS] - roots[:, : GAP_ROOTS - 1]
    finite_gaps = np.isfinite(gap_matrix)
    if np.any(finite_gaps):
        safe_gaps = np.where(finite_gaps, gap_matrix, np.nan)
        beta_idx, gap_idx = np.unravel_index(int(np.nanargmin(safe_gaps)), safe_gaps.shape)
        g_min = float(safe_gaps[beta_idx, gap_idx])
        beta_at_g_min = float(beta_values[beta_idx])
        pair_index = int(gap_idx) + 1
        lambda_lower = float(roots[beta_idx, gap_idx])
        lambda_upper = float(roots[beta_idx, gap_idx + 1])
    else:
        g_min = float("nan")
        beta_at_g_min = float("nan")
        pair_index = -1
        lambda_lower = float("nan")
        lambda_upper = float("nan")
        notes.append("no finite adjacent gaps among roots 1..7")

    branch_ranges: list[float] = []
    branch_rel_ranges: list[float] = []
    for root_idx in range(DISPLAY_ROOTS):
        branch_values = roots[:, root_idx]
        finite = branch_values[np.isfinite(branch_values)]
        if finite.size:
            sensitivity = float(np.max(finite) - np.min(finite))
        else:
            sensitivity = float("nan")
            notes.append(f"branch {root_idx + 1} has no finite values")
        branch_ranges.append(sensitivity)

        base = float(roots[0, root_idx])
        if np.isfinite(base) and abs(base) > RELATIVE_DENOMINATOR_TOL and np.isfinite(sensitivity):
            branch_rel_ranges.append(float(sensitivity / base))
        else:
            branch_rel_ranges.append(float("nan"))
            notes.append(f"relative sensitivity skipped for branch {root_idx + 1}")

    branch_array = np.asarray(branch_ranges, dtype=float)
    branch_rel_array = np.asarray(branch_rel_ranges, dtype=float)
    s_mean = float(np.nanmean(branch_array)) if np.any(np.isfinite(branch_array)) else float("nan")
    s_max = float(np.nanmax(branch_array)) if np.any(np.isfinite(branch_array)) else float("nan")
    s_mean_rel = float(np.nanmean(branch_rel_array)) if np.any(np.isfinite(branch_rel_array)) else float("nan")

    row: dict[str, object] = {
        "eta": number_text(eta),
        "mu": number_text(mu),
        "epsilon": number_text(epsilon),
        "g_min": number_text(g_min),
        "beta_at_g_min": number_text(beta_at_g_min),
        "pair_at_g_min": pair_label(pair_index) if pair_index > 0 else "",
        "lambda_lower_at_g_min": number_text(lambda_lower),
        "lambda_upper_at_g_min": number_text(lambda_upper),
        "S_mean": number_text(s_mean),
        "S_max": number_text(s_max),
        "S_mean_rel": number_text(s_mean_rel),
        "solver_warning_count": solver_warning_count,
        "missing_root_count": missing_root_count,
        "notes": "; ".join(notes),
        "_pair_index_at_g_min": pair_index,
        "_eta_float": float(eta),
        "_mu_float": float(mu),
        "_g_min_float": g_min,
        "_beta_at_g_min_float": beta_at_g_min,
        "_lambda_lower_float": lambda_lower,
        "_lambda_upper_float": lambda_upper,
        "_S_mean_float": s_mean,
        "_S_max_float": s_max,
        "_S_mean_rel_float": s_mean_rel,
    }
    for idx, value in enumerate(branch_ranges, start=1):
        row[f"S_branch{idx}"] = number_text(value)
        row[f"_S_branch{idx}_float"] = value
    for idx, value in enumerate(branch_rel_ranges, start=1):
        row[f"S_rel_branch{idx}"] = number_text(value)
        row[f"_S_rel_branch{idx}_float"] = value
    return row


def run_heatmap_metrics(
    *,
    eta_values: np.ndarray,
    mu_values: np.ndarray,
    beta_values: np.ndarray,
    epsilon: float,
    n_roots: int,
    root_scan_step: float,
    root_lmax0: float,
    workers: int,
) -> list[dict[str, object]]:
    tasks = [
        (
            float(eta),
            float(mu),
            tuple(float(value) for value in beta_values),
            float(epsilon),
            int(n_roots),
            float(root_scan_step),
            float(root_lmax0),
        )
        for mu in mu_values
        for eta in eta_values
    ]
    rows: list[dict[str, object]] = []
    total = len(tasks)
    if int(workers) <= 1:
        for idx, task in enumerate(tasks, start=1):
            rows.append(compute_heatmap_metric(task))
            if idx == 1 or idx % 25 == 0 or idx == total:
                print(f"heatmap metrics: {idx}/{total} eta-mu points")
        return rows

    with ProcessPoolExecutor(max_workers=int(workers)) as executor:
        for idx, row in enumerate(executor.map(compute_heatmap_metric, tasks, chunksize=4), start=1):
            rows.append(row)
            if idx == 1 or idx % 25 == 0 or idx == total:
                print(f"heatmap metrics: {idx}/{total} eta-mu points")
    return rows


def row_lookup(rows: Sequence[dict[str, object]]) -> dict[tuple[float, float], dict[str, object]]:
    return {
        (round(float(row["_eta_float"]), 12), round(float(row["_mu_float"]), 12)): row
        for row in rows
    }


def metric_matrix(
    rows: Sequence[dict[str, object]],
    eta_values: np.ndarray,
    mu_values: np.ndarray,
    key: str,
) -> np.ndarray:
    lookup = row_lookup(rows)
    matrix = np.full((len(mu_values), len(eta_values)), np.nan, dtype=float)
    for mu_idx, mu in enumerate(mu_values):
        for eta_idx, eta in enumerate(eta_values):
            row = lookup.get((round(float(eta), 12), round(float(mu), 12)))
            if row is not None:
                matrix[mu_idx, eta_idx] = float(row[key])
    return matrix


def save_heatmap(
    matrix: np.ndarray,
    eta_values: np.ndarray,
    mu_values: np.ndarray,
    output: Path,
    *,
    title: str,
    colorbar_label: str,
    log_if_dynamic: bool = False,
    cmap: str | mcolors.Colormap = "viridis",
    norm: mcolors.Normalize | None = None,
    colorbar_ticks: Sequence[float] | None = None,
    colorbar_ticklabels: Sequence[str] | None = None,
) -> None:
    extent = [float(eta_values[0]), float(eta_values[-1]), float(mu_values[0]), float(mu_values[-1])]
    data = np.asarray(matrix, dtype=float)
    local_norm = norm
    if local_norm is None and log_if_dynamic:
        finite = data[np.isfinite(data) & (data > 0.0)]
        if finite.size and float(np.max(finite) / np.min(finite)) >= 100.0:
            local_norm = mcolors.LogNorm(vmin=float(np.min(finite)), vmax=float(np.max(finite)))

    fig, ax = plt.subplots(figsize=(8.7, 5.5))
    image = ax.imshow(
        data,
        origin="lower",
        aspect="auto",
        extent=extent,
        interpolation="nearest",
        cmap=cmap,
        norm=local_norm,
    )
    ax.set_xlabel("eta")
    ax.set_ylabel("mu")
    ax.set_title(title)
    cbar = fig.colorbar(image, ax=ax)
    cbar.set_label(colorbar_label)
    if colorbar_ticks is not None:
        cbar.set_ticks(list(colorbar_ticks))
    if colorbar_ticklabels is not None:
        cbar.ax.set_yticklabels(list(colorbar_ticklabels))
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_heatmaps(rows: Sequence[dict[str, object]], eta_values: np.ndarray, mu_values: np.ndarray) -> None:
    save_heatmap(
        metric_matrix(rows, eta_values, mu_values, "_g_min_float"),
        eta_values,
        mu_values,
        HEATMAP_G_MIN_PNG,
        title="diagnostic minimum gap over beta",
        colorbar_label="min adjacent gap",
        log_if_dynamic=True,
    )
    save_heatmap(
        metric_matrix(rows, eta_values, mu_values, "_beta_at_g_min_float"),
        eta_values,
        mu_values,
        HEATMAP_BETA_G_MIN_PNG,
        title="diagnostic beta at minimum gap",
        colorbar_label="beta at min gap (deg)",
        cmap="plasma",
    )
    pair_cmap = plt.get_cmap("tab10", 6)
    pair_norm = mcolors.BoundaryNorm(np.arange(0.5, 7.5, 1.0), pair_cmap.N)
    save_heatmap(
        metric_matrix(rows, eta_values, mu_values, "_pair_index_at_g_min"),
        eta_values,
        mu_values,
        HEATMAP_PAIR_G_MIN_PNG,
        title="diagnostic pair at minimum gap",
        colorbar_label="pair at min gap",
        cmap=pair_cmap,
        norm=pair_norm,
        colorbar_ticks=[1, 2, 3, 4, 5, 6],
        colorbar_ticklabels=["1-2", "2-3", "3-4", "4-5", "5-6", "6-7"],
    )
    save_heatmap(
        metric_matrix(rows, eta_values, mu_values, "_S_mean_float"),
        eta_values,
        mu_values,
        HEATMAP_S_MEAN_PNG,
        title="diagnostic mean beta sensitivity",
        colorbar_label="S_mean",
        cmap="magma",
    )
    save_heatmap(
        metric_matrix(rows, eta_values, mu_values, "_S_max_float"),
        eta_values,
        mu_values,
        HEATMAP_S_MAX_PNG,
        title="diagnostic max beta sensitivity",
        colorbar_label="S_max",
        cmap="magma",
    )
    save_heatmap(
        metric_matrix(rows, eta_values, mu_values, "_S_mean_rel_float"),
        eta_values,
        mu_values,
        HEATMAP_S_MEAN_REL_PNG,
        title="diagnostic mean relative beta sensitivity",
        colorbar_label="S_mean_rel",
        cmap="cividis",
    )


def finite_rows(rows: Sequence[dict[str, object]], key: str) -> list[dict[str, object]]:
    return [row for row in rows if np.isfinite(float(row[key]))]


def best_row(rows: Sequence[dict[str, object]], key: str, *, largest: bool) -> dict[str, object]:
    local = finite_rows(rows, key)
    if not local:
        raise RuntimeError(f"No finite rows for {key}")
    return max(local, key=lambda row: float(row[key])) if largest else min(local, key=lambda row: float(row[key]))


def row_summary(row: dict[str, object], *, include_gap: bool = True) -> str:
    parts = [
        f"eta={float(row['_eta_float']):.12g}",
        f"mu={float(row['_mu_float']):.12g}",
    ]
    if include_gap:
        parts.extend(
            [
                f"beta={float(row['_beta_at_g_min_float']):.12g} deg",
                f"pair={row['pair_at_g_min']}",
                f"gap={float(row['_g_min_float']):.12g}",
            ]
        )
    return ", ".join(parts)


def top_close_rows(rows: Sequence[dict[str, object]], *, limit: int = 8) -> list[dict[str, object]]:
    candidates = sorted(finite_rows(rows, "_g_min_float"), key=lambda row: float(row["_g_min_float"]))
    return candidates[: int(limit)]


def paired_eta_asymmetry(
    rows: Sequence[dict[str, object]],
    eta_values: np.ndarray,
    mu_values: np.ndarray,
    key: str,
) -> tuple[float, float, tuple[float, float, float] | None]:
    lookup = row_lookup(rows)
    diffs: list[tuple[float, float, float]] = []
    for mu in mu_values:
        for eta in eta_values:
            eta_f = float(eta)
            if eta_f <= 0.0:
                continue
            plus = lookup.get((round(eta_f, 12), round(float(mu), 12)))
            minus = lookup.get((round(-eta_f, 12), round(float(mu), 12)))
            if plus is None or minus is None:
                continue
            plus_value = float(plus[key])
            minus_value = float(minus[key])
            if np.isfinite(plus_value) and np.isfinite(minus_value):
                diffs.append((abs(plus_value - minus_value), eta_f, float(mu)))
    if not diffs:
        return float("nan"), float("nan"), None
    values = np.array([item[0] for item in diffs], dtype=float)
    max_item = max(diffs, key=lambda item: item[0])
    return float(np.mean(values)), float(max_item[0]), max_item


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_report(
    *,
    heatmap_rows: Sequence[dict[str, object]],
    eta_values: np.ndarray,
    mu_values: np.ndarray,
    beta_values: np.ndarray,
    epsilon: float,
    heatmap_eta_step: float,
    heatmap_mu_step: float,
    heatmap_beta_step: float,
    heatmap_roots: int,
    root_scan_step: float,
    root_lmax0: float,
    lambda_eta_warning_count: int,
    lambda_eta_missing_count: int,
    lambda_beta_warning_count: int,
    lambda_beta_missing_count: int,
) -> None:
    min_gap = best_row(heatmap_rows, "_g_min_float", largest=False)
    max_s_mean = best_row(heatmap_rows, "_S_mean_float", largest=True)
    max_s_max = best_row(heatmap_rows, "_S_max_float", largest=True)
    max_s_mean_rel = best_row(heatmap_rows, "_S_mean_rel_float", largest=True)
    total_solver_warnings = sum(int(row["solver_warning_count"]) for row in heatmap_rows)
    total_missing_roots = sum(int(row["missing_root_count"]) for row in heatmap_rows)
    near_symmetric = (
        abs(float(min_gap["_eta_float"])) <= float(heatmap_eta_step) + 1.0e-12
        and float(min_gap["_mu_float"]) <= float(heatmap_mu_step) + 1.0e-12
    )
    mean_s_asym, max_s_asym, max_s_asym_item = paired_eta_asymmetry(
        heatmap_rows,
        eta_values,
        mu_values,
        "_S_mean_float",
    )
    mean_gap_asym, max_gap_asym, max_gap_asym_item = paired_eta_asymmetry(
        heatmap_rows,
        eta_values,
        mu_values,
        "_g_min_float",
    )

    close_rows = top_close_rows(heatmap_rows, limit=8)
    close_table = [
        [
            number_text(row["_eta_float"]),
            number_text(row["_mu_float"]),
            number_text(row["_beta_at_g_min_float"]),
            str(row["pair_at_g_min"]),
            number_text(row["_g_min_float"]),
        ]
        for row in close_rows
    ]
    recommendation_rows = close_rows[:4]
    for row in (max_s_mean, max_s_max):
        if row not in recommendation_rows:
            recommendation_rows.append(row)

    lines = [
        "# Diagnostic Eta-Mu-Beta Frequency Maps",
        "",
        "## 1. What Was Computed",
        "",
        "Diagnostic-only Euler-Bernoulli analytic thickness-mismatch maps using",
        "sorted frequencies as the primary quantity. The run generated:",
        "",
        f"- `Lambda(eta)` at `beta=0 deg`, `mu=0`, `epsilon={epsilon:g}`",
        "- sorted `Lambda(beta)` slices for the requested mu/eta combinations",
        "- eta-mu heatmaps of minimum adjacent sorted gaps over beta",
        "- eta-mu heatmaps of beta sensitivity metrics for the first six sorted roots",
        "",
        "No descendant branch tracking is used as a source of truth in these maps.",
        "No FEM, Gmsh, CalculiX, Timoshenko model, article workspace, article",
        "figures, old determinant, old solvers, or baseline results were touched.",
        "",
        "## 2. Grid Used",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- epsilon: `{epsilon:g}`",
        f"- Lambda(eta): eta `{LAMBDA_ETA_MIN:g}..{LAMBDA_ETA_MAX:g}` step `{LAMBDA_ETA_STEP:g}`, beta `0 deg`, mu `0`, first `{DISPLAY_ROOTS}` sorted roots",
        f"- Lambda(beta) slices: beta `0..90 deg` step `{SLICE_BETA_STEP:g}`, mu `{', '.join(str(value) for value in SLICE_MU_VALUES)}`, eta `{', '.join(str(value) for value in SLICE_ETA_VALUES)}`, first `{DISPLAY_ROOTS}` plotted roots",
        f"- heatmap eta grid: `{float(eta_values[0]):g}..{float(eta_values[-1]):g}` step `{heatmap_eta_step:g}` ({len(eta_values)} samples)",
        f"- heatmap mu grid: `{float(mu_values[0]):g}..{float(mu_values[-1]):g}` step `{heatmap_mu_step:g}` ({len(mu_values)} samples)",
        f"- heatmap beta grid: `{float(beta_values[0]):g}..{float(beta_values[-1]):g}` step `{heatmap_beta_step:g}` ({len(beta_values)} samples)",
        f"- heatmap sorted roots solved per beta: first `{heatmap_roots}`",
        f"- adjacent gaps audited: sorted pairs `1-2` through `6-7`",
        f"- root scan step: `{root_scan_step:g}`; root Lmax0: `{root_lmax0:g}`",
        "",
        "The heatmap uses the allowed fallback grid rather than the denser suggested",
        "`eta/mu/beta = 0.02/0.02/0.5` grid because that denser run would require",
        "roughly 424k determinant root solves.",
        "",
        "## 3. Definition Of g_min",
        "",
        "`g_min(eta, mu)` is the minimum over `beta in [0, 90 deg]` and adjacent",
        "sorted pairs `k-(k+1)` for `k=1..6` of",
        "`Lambda_sorted_{k+1}(beta; eta, mu) - Lambda_sorted_k(beta; eta, mu)`.",
        "",
        "## 4. Definition Of S Metrics",
        "",
        "`S_k(eta, mu)` is `max_beta Lambda_sorted_k - min_beta Lambda_sorted_k`",
        "for sorted branch `k`. `S_mean` is the mean of `S_k` over `k=1..6`.",
        "`S_max` is the maximum over `k=1..6`. `S_mean_rel` is the mean of",
        "`S_k / Lambda_sorted_k(beta=0)` for valid nonzero denominators.",
        "",
        "## 5. Smallest g_min",
        "",
        f"- location: {row_summary(min_gap)}",
        f"- lower Lambda: `{number_text(min_gap['_lambda_lower_float'])}`",
        f"- upper Lambda: `{number_text(min_gap['_lambda_upper_float'])}`",
        f"- near symmetric eta~0, mu~0 by one-grid-step criterion: `{'yes' if near_symmetric else 'no'}`",
        "",
        "Smallest gap candidates:",
        "",
    ]
    lines.extend(markdown_table(["eta", "mu", "beta deg", "pair", "g_min"], close_table))
    lines.extend(
        [
            "",
            "## 6. Symmetric-Case Check",
            "",
            (
                "The smallest heatmap gap occurs near the symmetric point by the one-step criterion."
                if near_symmetric
                else "The smallest heatmap gap does not occur near eta=0, mu=0 by the one-step criterion."
            ),
            "This is a diagnostic grid statement only; symmetric or near-symmetric",
            "close approaches still need strict local verification before article claims.",
            "",
            "## 7. Largest Angle Sensitivity",
            "",
            f"- max `S_mean`: `{number_text(max_s_mean['_S_mean_float'])}` at eta=`{number_text(max_s_mean['_eta_float'])}`, mu=`{number_text(max_s_mean['_mu_float'])}`",
            f"- max `S_max`: `{number_text(max_s_max['_S_max_float'])}` at eta=`{number_text(max_s_max['_eta_float'])}`, mu=`{number_text(max_s_max['_mu_float'])}`",
            f"- max `S_mean_rel`: `{number_text(max_s_mean_rel['_S_mean_rel_float'])}` at eta=`{number_text(max_s_mean_rel['_eta_float'])}`, mu=`{number_text(max_s_mean_rel['_mu_float'])}`",
            "",
            "## 8. Eta Sign Symmetry",
            "",
            f"- paired same-mu `S_mean(+eta)-S_mean(-eta)` mean absolute difference: `{number_text(mean_s_asym)}`",
            f"- paired same-mu `S_mean` max absolute difference: `{number_text(max_s_asym)}`"
            + (
                f" at eta=`{number_text(max_s_asym_item[1])}`, mu=`{number_text(max_s_asym_item[2])}`"
                if max_s_asym_item is not None
                else ""
            ),
            f"- paired same-mu `g_min(+eta)-g_min(-eta)` mean absolute difference: `{number_text(mean_gap_asym)}`",
            f"- paired same-mu `g_min` max absolute difference: `{number_text(max_gap_asym)}`"
            + (
                f" at eta=`{number_text(max_gap_asym_item[1])}`, mu=`{number_text(max_gap_asym_item[2])}`"
                if max_gap_asym_item is not None
                else ""
            ),
            "",
            "At fixed positive `mu`, positive and negative `eta` are generally",
            "asymmetric in this diagnostic map. The model's rod-exchange symmetry",
            "is a different comparison, `(mu, eta) -> (-mu, -eta)`, while this",
            "heatmap only scans nonnegative `mu`.",
            "",
            "## 9. Solver Warnings Or Missing Roots",
            "",
            f"- Lambda(eta) warning count: `{lambda_eta_warning_count}`",
            f"- Lambda(eta) missing root count: `{lambda_eta_missing_count}`",
            f"- Lambda(beta) slice warning count: `{lambda_beta_warning_count}`",
            f"- Lambda(beta) slice missing root count: `{lambda_beta_missing_count}`",
            f"- heatmap solver warning count: `{total_solver_warnings}`",
            f"- heatmap missing root count: `{total_missing_roots}`",
            "",
            "## 10. Diagnostic Conclusion",
            "",
            "Close approaches concentrate in the rows listed above. Strong beta",
            "sensitivity is represented by the max `S_mean` and max `S_max` rows.",
            "Small g_min regions should be checked by strict local gap verification",
            "before article claims.",
            "",
            "Recommended representative slices for future article-figure exploration:",
            "",
        ]
    )
    for row in recommendation_rows:
        lines.append(
            f"- eta=`{number_text(row['_eta_float'])}`, mu=`{number_text(row['_mu_float'])}`, "
            f"beta-near-gap=`{number_text(row['_beta_at_g_min_float'])}` deg, "
            f"pair `{row['pair_at_g_min']}`, g_min `{number_text(row['_g_min_float'])}`, "
            f"S_mean `{number_text(row['_S_mean_float'])}`, S_max `{number_text(row['_S_max_float'])}`"
        )
    lines.extend(
        [
            "",
            "These recommendations are diagnostic candidates, not crossing/no-crossing",
            "claims. Use the previous strict positive-gap verification machinery, or",
            "a new strict local check, before promoting any close approach.",
            "",
            "## Outputs",
            "",
            f"- Lambda(eta) PNG: `{LAMBDA_ETA_PNG.relative_to(REPO_ROOT)}`",
            f"- Lambda(eta) CSV: `{LAMBDA_ETA_CSV.relative_to(REPO_ROOT)}`",
            f"- Lambda(beta) slice CSV: `{LAMBDA_BETA_CSV.relative_to(REPO_ROOT)}`",
            f"- Lambda(beta) overview PNG: `{LAMBDA_BETA_OVERVIEW_PNG.relative_to(REPO_ROOT)}`",
            f"- heatmap CSV: `{HEATMAP_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{HEATMAP_REPORT.relative_to(REPO_ROOT)}`",
            f"- heatmaps: `{HEATMAP_G_MIN_PNG.relative_to(REPO_ROOT)}`, `{HEATMAP_BETA_G_MIN_PNG.relative_to(REPO_ROOT)}`, `{HEATMAP_PAIR_G_MIN_PNG.relative_to(REPO_ROOT)}`, `{HEATMAP_S_MEAN_PNG.relative_to(REPO_ROOT)}`, `{HEATMAP_S_MAX_PNG.relative_to(REPO_ROOT)}`, `{HEATMAP_S_MEAN_REL_PNG.relative_to(REPO_ROOT)}`",
            "",
        ]
    )
    HEATMAP_REPORT.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--root-lmax0", type=float, default=DEFAULT_ROOT_LMAX0)
    parser.add_argument("--heatmap-eta-step", type=float, default=HEATMAP_ETA_STEP)
    parser.add_argument("--heatmap-mu-step", type=float, default=HEATMAP_MU_STEP)
    parser.add_argument("--heatmap-beta-step", type=float, default=HEATMAP_BETA_STEP)
    parser.add_argument("--heatmap-roots", type=int, default=HEATMAP_ROOTS)
    parser.add_argument(
        "--workers",
        type=int,
        default=max(1, min(8, os.cpu_count() or 1)),
        help="Parallel workers for eta-mu heatmap rows.",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Use a tiny grid for wiring checks. Smoke outputs use the same file names.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    epsilon = float(args.epsilon)
    root_scan_step = float(args.root_scan_step)
    root_lmax0 = float(args.root_lmax0)

    lambda_eta_values = inclusive_grid(LAMBDA_ETA_MIN, LAMBDA_ETA_MAX, LAMBDA_ETA_STEP)
    slice_beta_values = inclusive_grid(0.0, 90.0, SLICE_BETA_STEP)
    heatmap_eta_step = float(args.heatmap_eta_step)
    heatmap_mu_step = float(args.heatmap_mu_step)
    heatmap_beta_step = float(args.heatmap_beta_step)
    heatmap_roots = int(args.heatmap_roots)
    workers = int(args.workers)

    if bool(args.smoke):
        lambda_eta_values = inclusive_grid(-0.05, 0.05, 0.05)
        slice_beta_values = inclusive_grid(0.0, 10.0, 5.0)
        heatmap_eta_step = 0.5
        heatmap_mu_step = 0.45
        heatmap_beta_step = 45.0
        heatmap_roots = GAP_ROOTS
        workers = 1

    heatmap_eta_values = inclusive_grid(-0.5, 0.5, heatmap_eta_step)
    heatmap_mu_values = inclusive_grid(0.0, 0.9, heatmap_mu_step)
    heatmap_beta_values = inclusive_grid(0.0, 90.0, heatmap_beta_step)

    print("diagnostic eta-mu-beta frequency maps")
    print(f"epsilon={epsilon:g}, root_scan_step={root_scan_step:g}, root_lmax0={root_lmax0:g}")

    lambda_eta_roots, lambda_eta_rows, lambda_eta_warnings, lambda_eta_missing, _lambda_eta_notes = compute_lambda_eta(
        epsilon=epsilon,
        eta_values=lambda_eta_values,
        root_scan_step=root_scan_step,
        root_lmax0=root_lmax0,
    )
    write_csv(LAMBDA_ETA_CSV, LAMBDA_ETA_FIELDNAMES, lambda_eta_rows)
    plot_lambda_eta(lambda_eta_values, lambda_eta_roots, LAMBDA_ETA_PNG, epsilon=epsilon)
    print(f"saved Lambda(eta): {LAMBDA_ETA_PNG}")

    lambda_beta_grids, lambda_beta_rows, lambda_beta_warnings, lambda_beta_missing, _lambda_beta_notes = (
        compute_lambda_beta_slices(
            epsilon=epsilon,
            beta_values=slice_beta_values,
            root_scan_step=root_scan_step,
            root_lmax0=root_lmax0,
        )
    )
    write_csv(LAMBDA_BETA_CSV, LAMBDA_BETA_FIELDNAMES, lambda_beta_rows)
    for (mu, eta), roots in lambda_beta_grids.items():
        output = lambda_beta_png_path(mu, eta)
        plot_lambda_beta_case(slice_beta_values, roots, output, epsilon=epsilon, mu=mu, eta=eta)
    plot_lambda_beta_overview(slice_beta_values, lambda_beta_grids, LAMBDA_BETA_OVERVIEW_PNG, epsilon=epsilon)
    print(f"saved Lambda(beta) slices: {len(lambda_beta_grids)} PNGs plus {LAMBDA_BETA_OVERVIEW_PNG}")

    print(
        "heatmap grid: "
        f"eta {heatmap_eta_values[0]:g}..{heatmap_eta_values[-1]:g} step {heatmap_eta_step:g} "
        f"({len(heatmap_eta_values)}), "
        f"mu {heatmap_mu_values[0]:g}..{heatmap_mu_values[-1]:g} step {heatmap_mu_step:g} "
        f"({len(heatmap_mu_values)}), "
        f"beta {heatmap_beta_values[0]:g}..{heatmap_beta_values[-1]:g} step {heatmap_beta_step:g} "
        f"({len(heatmap_beta_values)}), workers={workers}"
    )
    heatmap_rows = run_heatmap_metrics(
        eta_values=heatmap_eta_values,
        mu_values=heatmap_mu_values,
        beta_values=heatmap_beta_values,
        epsilon=epsilon,
        n_roots=heatmap_roots,
        root_scan_step=root_scan_step,
        root_lmax0=root_lmax0,
        workers=workers,
    )
    write_csv(HEATMAP_CSV, HEATMAP_FIELDNAMES, heatmap_rows)
    plot_heatmaps(heatmap_rows, heatmap_eta_values, heatmap_mu_values)
    write_report(
        heatmap_rows=heatmap_rows,
        eta_values=heatmap_eta_values,
        mu_values=heatmap_mu_values,
        beta_values=heatmap_beta_values,
        epsilon=epsilon,
        heatmap_eta_step=heatmap_eta_step,
        heatmap_mu_step=heatmap_mu_step,
        heatmap_beta_step=heatmap_beta_step,
        heatmap_roots=heatmap_roots,
        root_scan_step=root_scan_step,
        root_lmax0=root_lmax0,
        lambda_eta_warning_count=lambda_eta_warnings,
        lambda_eta_missing_count=lambda_eta_missing,
        lambda_beta_warning_count=lambda_beta_warnings,
        lambda_beta_missing_count=lambda_beta_missing,
    )
    print(f"saved heatmap CSV: {HEATMAP_CSV}")
    print(f"saved heatmap report: {HEATMAP_REPORT}")
    min_gap = best_row(heatmap_rows, "_g_min_float", largest=False)
    max_s_mean = best_row(heatmap_rows, "_S_mean_float", largest=True)
    max_s_max = best_row(heatmap_rows, "_S_max_float", largest=True)
    print(f"min g_min: {row_summary(min_gap)}")
    print(
        "max S_mean: "
        f"eta={float(max_s_mean['_eta_float']):.12g}, mu={float(max_s_mean['_mu_float']):.12g}, "
        f"S_mean={float(max_s_mean['_S_mean_float']):.12g}"
    )
    print(
        "max S_max: "
        f"eta={float(max_s_max['_eta_float']):.12g}, mu={float(max_s_max['_mu_float']):.12g}, "
        f"S_max={float(max_s_max['_S_max_float']):.12g}"
    )
    return {
        "lambda_eta_png": LAMBDA_ETA_PNG,
        "lambda_beta_png_count": len(lambda_beta_grids),
        "heatmap_rows": len(heatmap_rows),
        "report": HEATMAP_REPORT,
    }


if __name__ == "__main__":
    main()
