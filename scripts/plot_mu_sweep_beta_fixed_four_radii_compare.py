from __future__ import annotations

import numpy as np
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.compare_beta0_analytic_vs_fem import (  # noqa: E402
    N_BENDING,
    N_FEM_TYPED,
    build_params,
    fem_axial_fractions,
    fem_parameter_override,
    relative_error,
)
from scripts.plot_mu_sweep_beta0_four_radii_compare import (  # noqa: E402
    RADIUS_VALUES,
    build_combined_rows,
    frequencies_to_lambdas,
    load_csv_rows,
    plot_reference_curves,
    radius_key,
    rows_to_case_data,
    write_csv_rows,
)
from my_project.analytic.formulas import BeamParams, lambdas_to_frequencies  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.sweep_grid_policy import PRESENTATION_MU_STEP, presentation_mu_grid  # noqa: E402


BETA_VALUES = (7.5, 15.0)
MU_VALUES = presentation_mu_grid()
N_PLOTTED_BRANCHES = N_BENDING

ANALYTIC_LMIN = 0.2
ANALYTIC_LMAX0 = 70.0
ANALYTIC_SCAN_STEP = 0.01
ANALYTIC_GROW_FACTOR = 1.35
ANALYTIC_MAX_TRIES = 8

RESULTS_DIR = REPO_ROOT / "results"
FIGURE_PATH_BY_BETA = {
    7.5: RESULTS_DIR / "mu_sweep_beta7p5_four_radii_compare.png",
    15.0: RESULTS_DIR / "mu_sweep_beta15_four_radii_compare.png",
}
COMBINED_CSV_PATH_BY_BETA = {
    7.5: RESULTS_DIR / "mu_sweep_beta7p5_four_radii_compare.csv",
    15.0: RESULTS_DIR / "mu_sweep_beta15_four_radii_compare.csv",
}
CSV_PATHS_BY_BETA = {
    7.5: {
        0.005: RESULTS_DIR / "mu_sweep_beta7p5_r005_compare.csv",
        0.01: RESULTS_DIR / "mu_sweep_beta7p5_r010_compare.csv",
        0.015: RESULTS_DIR / "mu_sweep_beta7p5_r015_compare.csv",
        0.02: RESULTS_DIR / "mu_sweep_beta7p5_r020_compare.csv",
    },
    15.0: {
        0.005: RESULTS_DIR / "mu_sweep_beta15_r005_compare.csv",
        0.01: RESULTS_DIR / "mu_sweep_beta15_r010_compare.csv",
        0.015: RESULTS_DIR / "mu_sweep_beta15_r015_compare.csv",
        0.02: RESULTS_DIR / "mu_sweep_beta15_r020_compare.csv",
    },
}


def beta_key(beta_deg: float) -> float:
    return round(float(beta_deg), 6)


def beta_title(beta_deg: float) -> str:
    return f"{int(beta_deg)}" if float(beta_deg).is_integer() else f"{beta_deg:.1f}"


def solve_fem_modes_beta(
    params: BeamParams,
    mu: float,
    beta_deg: float,
    n_modes: int,
) -> tuple[np.ndarray, np.ndarray]:
    with fem_parameter_override(params):
        omega, vecs = fem.fem_solve(mu=mu, beta_deg=beta_deg, n_modes=n_modes)
        return omega * fem.scale, vecs


def analytic_sorted_lambdas(
    params: BeamParams,
    beta_deg: float,
    mu: float,
    n_modes: int,
) -> np.ndarray:
    roots = find_first_n_roots(
        beta=np.deg2rad(beta_deg),
        mu=float(mu),
        eps=params.eps,
        n_roots=n_modes,
        Lmin=ANALYTIC_LMIN,
        Lmax0=ANALYTIC_LMAX0,
        scan_step=ANALYTIC_SCAN_STEP,
        grow_factor=ANALYTIC_GROW_FACTOR,
        max_tries=ANALYTIC_MAX_TRIES,
    )
    return np.asarray(roots, dtype=float)


def compute_case_rows(beta_deg: float, radius: float) -> list[dict[str, float | int | str]]:
    params = build_params(radius)

    analytic_lambda = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)
    analytic_hz = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)
    fem_hz = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)
    fem_lambda = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)
    fem_mode_ids = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), -1, dtype=int)
    fem_axial_fraction = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)

    for col, mu in enumerate(MU_VALUES):
        roots = analytic_sorted_lambdas(
            params=params,
            beta_deg=beta_deg,
            mu=float(mu),
            n_modes=N_PLOTTED_BRANCHES,
        )
        analytic_lambda[:, col] = roots[:N_PLOTTED_BRANCHES]
        analytic_hz[:, col] = lambdas_to_frequencies(roots[:N_PLOTTED_BRANCHES], params)

        cur_fem_hz, cur_vecs = solve_fem_modes_beta(
            params=params,
            mu=float(mu),
            beta_deg=beta_deg,
            n_modes=N_FEM_TYPED,
        )
        cur_fem_hz = np.asarray(cur_fem_hz, dtype=float)
        cur_fem_lambda = frequencies_to_lambdas(cur_fem_hz, params)
        cur_axial_fraction = fem_axial_fractions(cur_vecs)

        fem_hz[:, col] = cur_fem_hz[:N_PLOTTED_BRANCHES]
        fem_lambda[:, col] = cur_fem_lambda[:N_PLOTTED_BRANCHES]
        fem_mode_ids[:, col] = np.arange(1, N_PLOTTED_BRANCHES + 1, dtype=int)
        fem_axial_fraction[:, col] = cur_axial_fraction[:N_PLOTTED_BRANCHES]

    rows: list[dict[str, float | int | str]] = []
    for branch_idx in range(N_PLOTTED_BRANCHES):
        branch_id = f"bending_{branch_idx + 1}"
        for col, mu in enumerate(MU_VALUES):
            analytic_hz_value = float(analytic_hz[branch_idx, col])
            fem_hz_value = float(fem_hz[branch_idx, col])
            rows.append(
                {
                    "beta": float(beta_deg),
                    "radius": float(radius),
                    "radius_mm": float(radius * 1e3),
                    "mu": float(mu),
                    "analytic_branch_id": branch_id,
                    "mode_type": "bending",
                    "fem_mode_id": int(fem_mode_ids[branch_idx, col]),
                    "analytic_lambda": float(analytic_lambda[branch_idx, col]),
                    "fem_lambda": float(fem_lambda[branch_idx, col]),
                    "analytic_hz": analytic_hz_value,
                    "fem_hz": fem_hz_value,
                    "relative_error": float(relative_error(analytic_hz_value, fem_hz_value)),
                    "fem_axial_fraction": float(fem_axial_fraction[branch_idx, col]),
                }
            )
    return rows


def get_case_data(beta_deg: float, radius: float) -> dict[str, object]:
    csv_path = CSV_PATHS_BY_BETA[beta_key(beta_deg)][radius_key(radius)]
    if csv_path.exists():
        rows = load_csv_rows(csv_path)
        return rows_to_case_data(radius=radius, rows=rows, source_csv=csv_path)

    rows = compute_case_rows(beta_deg=beta_deg, radius=radius)
    write_csv_rows(csv_path, rows)
    return rows_to_case_data(radius=radius, rows=rows, source_csv=csv_path)


def plot_four_radii_for_beta(beta_deg: float, cases: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(17.5, 13.2), squeeze=False)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, case in zip(axes.ravel(), cases):
        mu_values = np.asarray(case["mu_values"], dtype=float)
        analytic_lambda = np.asarray(case["analytic_lambda"], dtype=float)
        fem_lambda = np.asarray(case["fem_lambda"], dtype=float)
        radius = float(case["radius"])

        plot_reference_curves(ax=ax, case=case, colors=colors)

        for branch_pos in range(min(N_PLOTTED_BRANCHES, analytic_lambda.shape[0])):
            color = colors[branch_pos % len(colors)]
            ax.plot(mu_values, analytic_lambda[branch_pos], linewidth=2.0, color=color, zorder=3)
            ax.plot(
                mu_values,
                fem_lambda[branch_pos],
                linestyle="None",
                marker="o",
                markersize=3.3,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
                zorder=4,
            )

        y_main = float(np.nanmax(np.concatenate((analytic_lambda.ravel(), fem_lambda.ravel()))))
        ax.set_xlim(float(mu_values[0]), float(mu_values[-1]))
        ax.set_ylim(0.0, 1.08 * y_main)
        ax.set_title(f"beta = {beta_title(beta_deg)} deg, r = {radius:.3f}")
        ax.set_xlabel("mu")
        ax.set_ylabel(r"Dimensionless frequency parameter $\Lambda$")
        ax.grid(True, alpha=0.3)

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="solid - analytic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.0, label="markers - FEM"),
        Line2D(
            [0],
            [0],
            color=(0.0, 0.0, 0.0, 0.24),
            lw=1.1,
            linestyle="--",
            label=r"dashed - single-rod CS, $L^+ = \ell(1+\mu)$",
        ),
        Line2D(
            [0],
            [0],
            color=(0.0, 0.0, 0.0, 0.18),
            lw=1.0,
            linestyle=(0, (2.0, 2.8)),
            label=r"dashed - single-rod CS, $L^- = \ell(1-\mu)$",
        ),
    ]

    fig.suptitle(
        f"Analytic vs FEM at beta = {beta_title(beta_deg)} deg in Lambda\n"
        "same four-radii style, first 6 low branches across mu",
        y=0.985,
    )
    fig.legend(handles=style_handles, loc="upper center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(FIGURE_PATH_BY_BETA[beta_key(beta_deg)], dpi=220, bbox_inches="tight")
    plt.close(fig)


def run_beta_case(beta_deg: float) -> dict[str, object]:
    cases = [get_case_data(beta_deg=beta_deg, radius=radius) for radius in RADIUS_VALUES]
    combined_rows = build_combined_rows(cases)
    write_csv_rows(COMBINED_CSV_PATH_BY_BETA[beta_key(beta_deg)], combined_rows)
    plot_four_radii_for_beta(beta_deg=beta_deg, cases=cases)

    return {
        "beta": float(beta_deg),
        "figure_path": FIGURE_PATH_BY_BETA[beta_key(beta_deg)],
        "combined_csv_path": COMBINED_CSV_PATH_BY_BETA[beta_key(beta_deg)],
        "cases": cases,
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    summaries = [run_beta_case(beta_deg) for beta_deg in BETA_VALUES]

    for summary in summaries:
        print(f"beta={float(summary['beta']):.1f} deg")
        print(f"  presentation mu base step: {PRESENTATION_MU_STEP:.3f}")
        print("  local mu refinement windows: none")
        print(f"  saved figure: {summary['figure_path']}")
        print(f"  saved combined csv: {summary['combined_csv_path']}")
        for case in summary["cases"]:
            print(f"  radius={float(case['radius']):.3f} m")
            print(f"    source csv: {case['source_csv']}")
            print(f"    shown branches: {', '.join(case['branch_ids'])}")
            print(f"    max relative error: {float(case['max_relative_error']):.6e}")


if __name__ == "__main__":
    main()
