from __future__ import annotations

import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.plot_beta_sweep_mu0_compare import (  # noqa: E402
    N_BRANCHES,
    N_FEM_SOLVE,
    build_params,
    build_rows,
    compute_analytic_branches,
    compute_fem_branches,
    presentation_beta_grid,
    resolve_beta_grid,
    write_csv_rows,
)


RADIUS_VALUES = (0.005, 0.01, 0.015, 0.02)
RESULTS_DIR = REPO_ROOT / "results"
FIGURE_PATH = RESULTS_DIR / "beta_sweep_mu0_four_radii_compare.png"
COMBINED_CSV_PATH = RESULTS_DIR / "beta_sweep_mu0_four_radii_compare.csv"
REUSED_CSV_BY_RADIUS = {
    0.005: RESULTS_DIR / "beta_sweep_mu0_r005_compare_smoothed.csv",
    0.01: RESULTS_DIR / "beta_sweep_mu0_r010_compare_smoothed.csv",
}
NEW_CSV_BY_RADIUS = {
    0.015: RESULTS_DIR / "beta_sweep_mu0_r015_compare.csv",
    0.02: RESULTS_DIR / "beta_sweep_mu0_r020_compare.csv",
}


def radius_key(radius: float) -> float:
    return round(float(radius), 6)


def branch_sort_key(branch_id: str) -> tuple[int, str]:
    try:
        return int(branch_id.split("_")[-1]), branch_id
    except Exception:
        return (10**9, branch_id)


def load_csv_rows(path: Path) -> list[dict[str, float | int | str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def rows_to_case_data(radius: float, rows: list[dict[str, float | int | str]], source_csv: Path) -> dict[str, object]:
    branch_ids = sorted({str(row["branch_id"]) for row in rows}, key=branch_sort_key)
    beta_values = np.asarray(sorted({float(row["beta"]) for row in rows}), dtype=float)

    analytic_lambda = np.full((len(branch_ids), len(beta_values)), np.nan, dtype=float)
    fem_lambda = np.full((len(branch_ids), len(beta_values)), np.nan, dtype=float)
    rel_error = np.full((len(branch_ids), len(beta_values)), np.nan, dtype=float)

    branch_lookup = {branch_id: idx for idx, branch_id in enumerate(branch_ids)}
    beta_lookup = {float(beta): idx for idx, beta in enumerate(beta_values)}

    for row in rows:
        i = branch_lookup[str(row["branch_id"])]
        j = beta_lookup[float(row["beta"])]
        analytic_lambda[i, j] = float(row["analytic_lambda"])
        fem_lambda[i, j] = float(row["fem_lambda"])
        rel_error[i, j] = float(row["relative_error"])

    return {
        "radius": radius,
        "rows": rows,
        "beta_values": beta_values,
        "branch_ids": branch_ids,
        "analytic_lambda": analytic_lambda,
        "fem_lambda": fem_lambda,
        "max_relative_error": float(np.nanmax(rel_error)),
        "source_csv": source_csv,
    }


def compute_case_rows(radius: float) -> list[dict[str, float | int | str]]:
    params = build_params(radius)
    beta_requested = presentation_beta_grid(radius)
    beta_values, _ = resolve_beta_grid(params=params, beta_values=beta_requested)
    analytic_lambda, analytic_hz = compute_analytic_branches(
        params=params,
        beta_values=beta_values,
        n_branches=N_BRANCHES,
    )
    fem_lambda, fem_hz, fem_mode_ids = compute_fem_branches(
        params=params,
        beta_values=beta_values,
        analytic_lambda=analytic_lambda,
        n_modes_solve=N_FEM_SOLVE,
    )
    return build_rows(
        beta_values=beta_values,
        analytic_lambda=analytic_lambda,
        analytic_hz=analytic_hz,
        fem_lambda=fem_lambda,
        fem_hz=fem_hz,
        fem_mode_ids=fem_mode_ids,
    )


def get_case_data(radius: float) -> dict[str, object]:
    reuse_path = REUSED_CSV_BY_RADIUS.get(radius_key(radius))
    if reuse_path is not None and reuse_path.exists():
        rows = load_csv_rows(reuse_path)
        return rows_to_case_data(radius=radius, rows=rows, source_csv=reuse_path)

    rows = compute_case_rows(radius=radius)
    save_path = NEW_CSV_BY_RADIUS[radius_key(radius)]
    write_csv_rows(save_path, rows)
    return rows_to_case_data(radius=radius, rows=rows, source_csv=save_path)


def build_combined_rows(cases: list[dict[str, object]]) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for case in cases:
        radius = float(case["radius"])
        source_csv = str(case["source_csv"])
        for row in case["rows"]:
            rows.append(
                {
                    "radius": radius,
                    "source_csv": source_csv,
                    **row,
                }
            )
    return rows


def plot_four_radii(cases: list[dict[str, object]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16.5, 12.5), squeeze=False)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, case in zip(axes.ravel(), cases):
        beta_values = np.asarray(case["beta_values"], dtype=float)
        analytic_lambda = np.asarray(case["analytic_lambda"], dtype=float)
        fem_lambda = np.asarray(case["fem_lambda"], dtype=float)
        radius = float(case["radius"])

        for branch_pos in range(min(N_BRANCHES, analytic_lambda.shape[0])):
            color = colors[branch_pos % len(colors)]
            ax.plot(beta_values, analytic_lambda[branch_pos], linewidth=2.0, color=color)
            ax.plot(
                beta_values,
                fem_lambda[branch_pos],
                linestyle="None",
                marker="o",
                markersize=3.2,
                color=color,
                markerfacecolor=color,
                markeredgecolor=color,
            )

        ax.set_title(f"mu = 0, r = {radius:.3f}")
        ax.set_xlabel("beta, degrees")
        ax.set_ylabel("Безразмерный частотный параметр Λ")
        ax.set_xlim(float(beta_values[0]), float(beta_values[-1]))
        ax.grid(True, alpha=0.3)

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="solid — analytic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=5.0, label="markers — FEM"),
    ]
    fig.suptitle("Analytic vs FEM beta-sweep comparison, mu = 0, first 10 branches", y=0.985)
    fig.legend(handles=style_handles, loc="upper center", ncol=2, fontsize=10, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=(0, 0, 1, 0.93))
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
    for case in cases:
        print(f"radius={float(case['radius']):.3f} m")
        print(f"  source csv: {case['source_csv']}")
        print(f"  shown branches: {', '.join(case['branch_ids'][:N_BRANCHES])}")
        print(f"  max relative error: {float(case['max_relative_error']):.6e}")


if __name__ == "__main__":
    main()
