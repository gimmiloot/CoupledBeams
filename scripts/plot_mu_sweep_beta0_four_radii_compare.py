from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import to_rgba
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
    analytic_bending_mu_sweep,
    build_params,
    fem_mode_rows,
    relative_error,
    select_mode_rows,
)
from my_project.analytic.FreqMuNet import roots_clamped_supported, single_lambda  # noqa: E402
from my_project.analytic.formulas import BeamParams, frequency_scale  # noqa: E402
from scripts.sweep_grid_policy import PRESENTATION_MU_STEP, presentation_mu_grid  # noqa: E402


RADIUS_VALUES = (0.005, 0.01, 0.015, 0.02)
MU_VALUES = presentation_mu_grid()
N_PLOTTED_BRANCHES = N_BENDING
N_REFERENCE_LPLUS = 6
N_REFERENCE_LMINUS = 3

RESULTS_DIR = REPO_ROOT / "results"
FIGURE_PATH = RESULTS_DIR / "mu_sweep_beta0_four_radii_compare.png"
COMBINED_CSV_PATH = RESULTS_DIR / "mu_sweep_beta0_four_radii_compare.csv"
CSV_PATH_BY_RADIUS = {
    0.005: RESULTS_DIR / "mu_sweep_beta0_r005_compare.csv",
    0.01: RESULTS_DIR / "mu_sweep_beta0_r010_compare.csv",
    0.015: RESULTS_DIR / "mu_sweep_beta0_r015_compare.csv",
    0.02: RESULTS_DIR / "mu_sweep_beta0_r020_compare.csv",
}


def radius_key(radius: float) -> float:
    return round(float(radius), 6)


def branch_sort_key(branch_id: str) -> tuple[int, str]:
    try:
        return int(branch_id.split("_")[-1]), branch_id
    except Exception:
        return (10**9, branch_id)


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def load_csv_rows(path: Path) -> list[dict[str, float | int | str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def frequencies_to_lambdas(freq_hz: np.ndarray, params: BeamParams) -> np.ndarray:
    return np.sqrt(np.maximum(np.asarray(freq_hz, dtype=float) / frequency_scale(params), 0.0))


def build_reference_curves(params: BeamParams, mu_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    roots_cs = roots_clamped_supported(max(N_REFERENCE_LPLUS, N_REFERENCE_LMINUS))
    l_base = params.L_base
    mu_values = np.asarray(mu_values, dtype=float)
    l_plus = l_base * (1.0 + mu_values)
    l_minus = l_base * (1.0 - mu_values)
    lambda_lplus = single_lambda(roots_cs[:N_REFERENCE_LPLUS], l_base, l_plus)
    lambda_lminus = single_lambda(roots_cs[:N_REFERENCE_LMINUS], l_base, l_minus)
    return lambda_lplus, lambda_lminus


def rows_to_case_data(
    radius: float,
    rows: list[dict[str, float | int | str]],
    source_csv: Path,
) -> dict[str, object]:
    params = build_params(radius)
    branch_ids = sorted({str(row["analytic_branch_id"]) for row in rows}, key=branch_sort_key)
    mu_values = np.asarray(sorted({float(row["mu"]) for row in rows}), dtype=float)

    analytic_lambda = np.full((len(branch_ids), len(mu_values)), np.nan, dtype=float)
    fem_lambda = np.full((len(branch_ids), len(mu_values)), np.nan, dtype=float)
    rel_error = np.full((len(branch_ids), len(mu_values)), np.nan, dtype=float)

    branch_lookup = {branch_id: idx for idx, branch_id in enumerate(branch_ids)}
    mu_lookup = {float(mu): idx for idx, mu in enumerate(mu_values)}

    for row in rows:
        i = branch_lookup[str(row["analytic_branch_id"])]
        j = mu_lookup[float(row["mu"])]
        analytic_lambda[i, j] = float(row["analytic_lambda"])
        fem_lambda[i, j] = float(row["fem_lambda"])
        rel_error[i, j] = float(row["relative_error"])

    lambda_lplus, lambda_lminus = build_reference_curves(params=params, mu_values=mu_values)
    return {
        "radius": radius,
        "params": params,
        "rows": rows,
        "mu_values": mu_values,
        "branch_ids": branch_ids,
        "analytic_lambda": analytic_lambda,
        "fem_lambda": fem_lambda,
        "reference_lplus": lambda_lplus,
        "reference_lminus": lambda_lminus,
        "max_relative_error": float(np.nanmax(rel_error)),
        "source_csv": source_csv,
    }


def compute_case_rows(radius: float) -> list[dict[str, float | int | str]]:
    params = build_params(radius)
    analytic_hz = analytic_bending_mu_sweep(
        params=params,
        mu_values=MU_VALUES,
        n_modes=N_PLOTTED_BRANCHES,
    )
    analytic_lambda = frequencies_to_lambdas(analytic_hz, params)

    fem_hz = np.full_like(analytic_hz, np.nan, dtype=float)
    fem_mode_ids = np.full(analytic_hz.shape, -1, dtype=int)
    fem_axial_fractions = np.full(analytic_hz.shape, np.nan, dtype=float)

    for col, mu in enumerate(MU_VALUES):
        mode_rows = fem_mode_rows(params, mu=float(mu), n_modes=N_FEM_TYPED)
        bending_rows = select_mode_rows(mode_rows, mode_type="bending", count=N_PLOTTED_BRANCHES)
        for branch_idx, fem_row in enumerate(bending_rows):
            fem_hz[branch_idx, col] = float(fem_row["frequency_hz"])
            fem_mode_ids[branch_idx, col] = int(fem_row["mode_id"])
            fem_axial_fractions[branch_idx, col] = float(fem_row["axial_fraction"])

    fem_lambda = frequencies_to_lambdas(fem_hz, params)

    rows: list[dict[str, float | int | str]] = []
    for branch_idx in range(N_PLOTTED_BRANCHES):
        branch_id = f"bending_{branch_idx + 1}"
        for col, mu in enumerate(MU_VALUES):
            analytic_hz_value = float(analytic_hz[branch_idx, col])
            fem_hz_value = float(fem_hz[branch_idx, col])
            rows.append(
                {
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
                    "fem_axial_fraction": float(fem_axial_fractions[branch_idx, col]),
                }
            )
    return rows


def get_case_data(radius: float) -> dict[str, object]:
    csv_path = CSV_PATH_BY_RADIUS[radius_key(radius)]
    if csv_path.exists():
        rows = load_csv_rows(csv_path)
        return rows_to_case_data(radius=radius, rows=rows, source_csv=csv_path)

    rows = compute_case_rows(radius=radius)
    write_csv_rows(csv_path, rows)
    return rows_to_case_data(radius=radius, rows=rows, source_csv=csv_path)


def build_combined_rows(cases: list[dict[str, object]]) -> list[dict[str, float | int | str]]:
    combined_rows: list[dict[str, float | int | str]] = []
    for case in cases:
        for row in case["rows"]:
            combined_rows.append(
                {
                    "source_csv": str(case["source_csv"]),
                    **row,
                }
            )
    return combined_rows


def plot_reference_curves(ax: plt.Axes, case: dict[str, object], colors: list[str]) -> None:
    mu_values = np.asarray(case["mu_values"], dtype=float)
    reference_lplus = np.asarray(case["reference_lplus"], dtype=float)
    reference_lminus = np.asarray(case["reference_lminus"], dtype=float)

    for idx in range(reference_lplus.shape[0]):
        ax.plot(
            mu_values,
            reference_lplus[idx],
            linewidth=1.1,
            linestyle="--",
            color=to_rgba(colors[idx % len(colors)], alpha=0.24),
            zorder=1,
        )

    for idx in range(reference_lminus.shape[0]):
        ax.plot(
            mu_values,
            reference_lminus[idx],
            linewidth=1.0,
            linestyle=(0, (2.0, 2.8)),
            color=to_rgba(colors[idx % len(colors)], alpha=0.18),
            zorder=1,
        )


def plot_four_radii(cases: list[dict[str, object]]) -> None:
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
        ax.set_title(f"beta = 0, r = {radius:.3f}")
        ax.set_xlabel("mu")
        ax.set_ylabel(r"Dimensionless frequency parameter $\Lambda$")
        ax.grid(True, alpha=0.3)

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="solid - analytic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.0, label="markers - FEM"),
        Line2D(
            [0],
            [0],
            color=to_rgba("black", alpha=0.24),
            lw=1.1,
            linestyle="--",
            label=r"dashed - single-rod CS, $L^+ = \ell(1+\mu)$",
        ),
        Line2D(
            [0],
            [0],
            color=to_rgba("black", alpha=0.18),
            lw=1.0,
            linestyle=(0, (2.0, 2.8)),
            label=r"dashed - single-rod CS, $L^- = \ell(1-\mu)$",
        ),
    ]

    fig.suptitle(
        "Analytic vs FEM at beta = 0 in Lambda\n"
        "type-aware comparison for the first 6 bending branches across four radii",
        y=0.985,
    )
    fig.legend(handles=style_handles, loc="upper center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(FIGURE_PATH, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    cases = [get_case_data(radius) for radius in RADIUS_VALUES]
    combined_rows = build_combined_rows(cases)
    write_csv_rows(COMBINED_CSV_PATH, combined_rows)
    plot_four_radii(cases)

    print(f"saved figure: {FIGURE_PATH}")
    print(f"saved combined csv: {COMBINED_CSV_PATH}")
    print(f"presentation mu base step: {PRESENTATION_MU_STEP:.3f}")
    print("local mu refinement windows: none")
    for case in cases:
        print(f"radius={float(case['radius']):.3f} m")
        print(f"  source csv: {case['source_csv']}")
        print(f"  shown branches: {', '.join(case['branch_ids'])}")
        print(f"  max relative error: {float(case['max_relative_error']):.6e}")


if __name__ == "__main__":
    main()
