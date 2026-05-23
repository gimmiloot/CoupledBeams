from __future__ import annotations

import csv
from pathlib import Path
import sys
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
LAMBDA_ETA_PNG = RESULTS_DIR / "thickness_mismatch_lambda_eta_beta15_eps0p0025.png"
LAMBDA_ETA_CSV = RESULTS_DIR / "thickness_mismatch_lambda_eta_beta15_eps0p0025.csv"
LAMBDA_MU_PNG = RESULTS_DIR / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_sweep.png"
LAMBDA_MU_CSV = RESULTS_DIR / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_sweep.csv"

BETA_DEG = 15.0
EPSILON = 0.0025
NUM_ROOTS = 6
MU_VALUES_FOR_ETA_SWEEP = (0.0, 0.3, 0.6)
ETA_VALUES = np.round(np.arange(-0.15, 0.1500001, 0.01), 10)
ETA_VALUES_FOR_MU_SWEEP = (-0.1, 0.0, 0.1)
MU_VALUES = np.round(np.arange(0.0, 0.9000001, 0.02), 10)
ROOT_SCAN_STEP = 0.015
ROOT_LMAX0 = 35.0

FIELDNAMES = ["beta_deg", "epsilon", "mu", "eta", "root_index", "Lambda"]


def write_csv(path: Path, rows: Iterable[dict[str, float | int]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def roots_for(*, beta_deg: float, mu: float, eta: float, epsilon: float, num_roots: int) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(beta_deg)),
        float(mu),
        float(epsilon),
        float(eta),
        int(num_roots),
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots)):
        raise RuntimeError(f"Missing roots for beta={beta_deg:g}, mu={mu:g}, eta={eta:g}.")
    return roots


def lambda_eta_rows() -> tuple[list[dict[str, float | int]], dict[float, np.ndarray]]:
    rows: list[dict[str, float | int]] = []
    grids: dict[float, np.ndarray] = {}
    for mu in MU_VALUES_FOR_ETA_SWEEP:
        grid = np.full((NUM_ROOTS, len(ETA_VALUES)), np.nan, dtype=float)
        for col, eta in enumerate(ETA_VALUES):
            roots = roots_for(beta_deg=BETA_DEG, mu=float(mu), eta=float(eta), epsilon=EPSILON, num_roots=NUM_ROOTS)
            grid[:, col] = roots
            for idx, value in enumerate(roots, start=1):
                rows.append(
                    {
                        "beta_deg": BETA_DEG,
                        "epsilon": EPSILON,
                        "mu": float(mu),
                        "eta": float(eta),
                        "root_index": int(idx),
                        "Lambda": float(value),
                    }
                )
        grids[float(mu)] = grid
    return rows, grids


def lambda_mu_rows() -> tuple[list[dict[str, float | int]], dict[float, np.ndarray]]:
    rows: list[dict[str, float | int]] = []
    grids: dict[float, np.ndarray] = {}
    for eta in ETA_VALUES_FOR_MU_SWEEP:
        grid = np.full((NUM_ROOTS, len(MU_VALUES)), np.nan, dtype=float)
        for col, mu in enumerate(MU_VALUES):
            roots = roots_for(beta_deg=BETA_DEG, mu=float(mu), eta=float(eta), epsilon=EPSILON, num_roots=NUM_ROOTS)
            grid[:, col] = roots
            for idx, value in enumerate(roots, start=1):
                rows.append(
                    {
                        "beta_deg": BETA_DEG,
                        "epsilon": EPSILON,
                        "mu": float(mu),
                        "eta": float(eta),
                        "root_index": int(idx),
                        "Lambda": float(value),
                    }
                )
        grids[float(eta)] = grid
    return rows, grids


def plot_lambda_eta(grids: dict[float, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(1, len(MU_VALUES_FOR_ETA_SWEEP), figsize=(10.5, 3.6), sharey=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, mu in zip(axes, MU_VALUES_FOR_ETA_SWEEP, strict=True):
        grid = grids[float(mu)]
        for root_idx in range(NUM_ROOTS):
            ax.plot(
                ETA_VALUES,
                grid[root_idx],
                lw=1.4,
                color=colors[root_idx % len(colors)],
                label=f"root {root_idx + 1}" if ax is axes[0] else "_nolegend_",
            )
        ax.axvline(0.0, color="0.55", lw=0.8, ls=":")
        ax.set_title(rf"$\mu={mu:g}$", fontsize=11)
        ax.set_xlabel(r"$\eta$")
        ax.grid(True, color="0.88", linewidth=0.6)
    axes[0].set_ylabel(r"$\Lambda$")
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    fig.suptitle(rf"Thickness mismatch diagnostic: $\beta={BETA_DEG:g}^\circ$, $\varepsilon={EPSILON:g}$", fontsize=12)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_lambda_mu(grids: dict[float, np.ndarray], output: Path) -> None:
    fig, axes = plt.subplots(1, len(ETA_VALUES_FOR_MU_SWEEP), figsize=(10.5, 3.6), sharey=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, eta in zip(axes, ETA_VALUES_FOR_MU_SWEEP, strict=True):
        grid = grids[float(eta)]
        for root_idx in range(NUM_ROOTS):
            ax.plot(
                MU_VALUES,
                grid[root_idx],
                lw=1.4,
                color=colors[root_idx % len(colors)],
                label=f"root {root_idx + 1}" if ax is axes[0] else "_nolegend_",
            )
        ax.set_title(rf"$\eta={eta:g}$", fontsize=11)
        ax.set_xlabel(r"$\mu$")
        ax.grid(True, color="0.88", linewidth=0.6)
    axes[0].set_ylabel(r"$\Lambda$")
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(rf"Thickness mismatch diagnostic: $\beta={BETA_DEG:g}^\circ$, $\varepsilon={EPSILON:g}$", fontsize=12)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def smoothness_observation(lambda_eta_grids: dict[float, np.ndarray]) -> tuple[float, float, int]:
    small_eta = np.array([-0.1, -0.05, 0.0, 0.05, 0.1], dtype=float)
    indices = [int(np.where(np.isclose(ETA_VALUES, value, atol=1e-12))[0][0]) for value in small_eta]
    max_jump = -np.inf
    max_mu = np.nan
    max_root = -1
    for mu, grid in lambda_eta_grids.items():
        small_grid = grid[:, indices]
        jumps = np.abs(np.diff(small_grid, axis=1))
        local = float(np.max(jumps))
        if local > max_jump:
            root_idx, _ = np.unravel_index(int(np.argmax(jumps)), jumps.shape)
            max_jump = local
            max_mu = float(mu)
            max_root = int(root_idx) + 1
    return float(max_jump), float(max_mu), int(max_root)


def main() -> dict[str, Path | float | int]:
    eta_rows, eta_grids = lambda_eta_rows()
    mu_rows, mu_grids = lambda_mu_rows()

    write_csv(LAMBDA_ETA_CSV, eta_rows, FIELDNAMES)
    write_csv(LAMBDA_MU_CSV, mu_rows, FIELDNAMES)
    plot_lambda_eta(eta_grids, LAMBDA_ETA_PNG)
    plot_lambda_mu(mu_grids, LAMBDA_MU_PNG)

    max_jump, max_jump_mu, max_jump_root = smoothness_observation(eta_grids)
    print(f"saved Lambda(eta) PNG: {LAMBDA_ETA_PNG}")
    print(f"saved Lambda(eta) CSV: {LAMBDA_ETA_CSV}")
    print(f"saved Lambda(mu) eta-sweep PNG: {LAMBDA_MU_PNG}")
    print(f"saved Lambda(mu) eta-sweep CSV: {LAMBDA_MU_CSV}")
    print(
        "small-eta smoothness scan: "
        f"max adjacent root jump={max_jump:.6e} at mu={max_jump_mu:g}, root={max_jump_root}"
    )
    return {
        "lambda_eta_png": LAMBDA_ETA_PNG,
        "lambda_eta_csv": LAMBDA_ETA_CSV,
        "lambda_mu_png": LAMBDA_MU_PNG,
        "lambda_mu_csv": LAMBDA_MU_CSV,
        "small_eta_max_adjacent_jump": max_jump,
        "small_eta_max_jump_mu": max_jump_mu,
        "small_eta_max_jump_root": max_jump_root,
    }


if __name__ == "__main__":
    main()
