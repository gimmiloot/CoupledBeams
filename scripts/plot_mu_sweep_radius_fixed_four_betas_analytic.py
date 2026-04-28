from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.compare_beta0_analytic_vs_fem import N_BENDING, build_params  # noqa: E402
from scripts.plot_mu_sweep_beta0_four_radii_compare import (  # noqa: E402
    branch_sort_key,
    build_combined_rows,
    build_reference_curves,
    load_csv_rows,
    plot_reference_curves,
    write_csv_rows,
)
from my_project.analytic.formulas import BeamParams, lambdas_to_frequencies  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402
from scripts.sweep_grid_policy import PRESENTATION_MU_STEP, presentation_mu_grid  # noqa: E402


FIXED_RADIUS = 0.015
DEFAULT_BETA_VALUES = (15.0, 30.0, 45.0, 60.0)
MU_VALUES = presentation_mu_grid()
N_PLOTTED_BRANCHES = N_BENDING

ANALYTIC_LMIN = 0.2
ANALYTIC_LMAX0 = 70.0
ANALYTIC_SCAN_STEP = 0.01
ANALYTIC_GROW_FACTOR = 1.35
ANALYTIC_MAX_TRIES = 8

RESULTS_DIR = REPO_ROOT / "results"


def beta_key(beta_deg: float) -> float:
    return round(float(beta_deg), 6)


def radius_tag(radius: float) -> str:
    return f"r{int(round(float(radius) * 1e3)):03d}"


def beta_slug(beta_deg: float) -> str:
    beta = beta_key(beta_deg)
    if float(beta).is_integer():
        return f"{int(beta)}"
    return f"{beta:.6f}".rstrip("0").rstrip(".").replace(".", "p")


def beta_title(beta_deg: float) -> str:
    return f"{int(beta_deg)}" if float(beta_deg).is_integer() else f"{beta_deg:.1f}"


def format_beta_list(beta_values: tuple[float, ...]) -> str:
    return ", ".join(beta_title(beta) for beta in beta_values)


def csv_path_for_beta(beta_deg: float) -> Path:
    return RESULTS_DIR / f"mu_sweep_{radius_tag(FIXED_RADIUS)}_beta{beta_slug(beta_deg)}_analytic.csv"


def figure_path_for_run() -> Path:
    return RESULTS_DIR / f"mu_sweep_{radius_tag(FIXED_RADIUS)}_selected_betas_analytic.png"


def combined_csv_path_for_run() -> Path:
    return RESULTS_DIR / f"mu_sweep_{radius_tag(FIXED_RADIUS)}_selected_betas_analytic.csv"


def normalize_beta_values(beta_values: list[float]) -> tuple[float, ...]:
    normalized: list[float] = []
    seen: set[float] = set()
    for beta in beta_values:
        beta_norm = beta_key(beta)
        if beta_norm not in seen:
            seen.add(beta_norm)
            normalized.append(beta_norm)
    return tuple(normalized)


def parse_beta_values(argv: list[str] | None = None) -> tuple[float, ...]:
    parser = argparse.ArgumentParser(
        description="Plot analytic mu-sweeps at fixed radius for selected beta values."
    )
    parser.add_argument(
        "--betas",
        nargs="*",
        default=None,
        help="Beta values in degrees, for example --betas 7.5 15 30 45 or --betas 7.5,15,30,45",
    )
    args = parser.parse_args(argv)

    if args.betas is None:
        return DEFAULT_BETA_VALUES
    if not args.betas:
        parser.error("--betas was provided without any values")

    beta_values: list[float] = []
    for token in args.betas:
        for part in token.split(","):
            text = part.strip()
            if not text:
                continue
            try:
                beta_values.append(float(text))
            except ValueError as exc:
                parser.error(f"could not parse beta value '{text}': {exc}")

    if not beta_values:
        parser.error("no beta values were parsed from --betas")
    return normalize_beta_values(beta_values)


def subplot_shape_for_case_count(case_count: int) -> tuple[int, int]:
    if case_count <= 0:
        raise ValueError("Need at least one beta value to plot.")
    if case_count == 1:
        return 1, 1
    if case_count == 2:
        return 1, 2
    ncols = int(math.ceil(math.sqrt(case_count)))
    nrows = int(math.ceil(case_count / ncols))
    return nrows, ncols


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


def compute_case_rows(beta_deg: float) -> list[dict[str, float | int | str]]:
    params = build_params(FIXED_RADIUS)
    analytic_lambda = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)
    analytic_hz = np.full((N_PLOTTED_BRANCHES, len(MU_VALUES)), np.nan, dtype=float)

    for col, mu in enumerate(MU_VALUES):
        roots = analytic_sorted_lambdas(
            params=params,
            beta_deg=beta_deg,
            mu=float(mu),
            n_modes=N_PLOTTED_BRANCHES,
        )
        analytic_lambda[:, col] = roots[:N_PLOTTED_BRANCHES]
        analytic_hz[:, col] = lambdas_to_frequencies(roots[:N_PLOTTED_BRANCHES], params)

    rows: list[dict[str, float | int | str]] = []
    for branch_idx in range(N_PLOTTED_BRANCHES):
        branch_id = f"bending_{branch_idx + 1}"
        for col, mu in enumerate(MU_VALUES):
            rows.append(
                {
                    "beta": float(beta_deg),
                    "radius": float(FIXED_RADIUS),
                    "radius_mm": float(FIXED_RADIUS * 1e3),
                    "mu": float(mu),
                    "analytic_branch_id": branch_id,
                    "mode_type": "bending",
                    "analytic_lambda": float(analytic_lambda[branch_idx, col]),
                    "analytic_hz": float(analytic_hz[branch_idx, col]),
                }
            )
    return rows


def rows_to_case_data(
    beta_deg: float,
    rows: list[dict[str, float | int | str]],
    source_csv: Path,
) -> dict[str, object]:
    params = build_params(FIXED_RADIUS)
    branch_ids = sorted({str(row["analytic_branch_id"]) for row in rows}, key=branch_sort_key)
    mu_values = np.asarray(sorted({float(row["mu"]) for row in rows}), dtype=float)

    analytic_lambda = np.full((len(branch_ids), len(mu_values)), np.nan, dtype=float)
    branch_lookup = {branch_id: idx for idx, branch_id in enumerate(branch_ids)}
    mu_lookup = {float(mu): idx for idx, mu in enumerate(mu_values)}

    for row in rows:
        i = branch_lookup[str(row["analytic_branch_id"])]
        j = mu_lookup[float(row["mu"])]
        analytic_lambda[i, j] = float(row["analytic_lambda"])

    reference_lplus, reference_lminus = build_reference_curves(params=params, mu_values=mu_values)
    return {
        "beta": float(beta_deg),
        "radius": float(FIXED_RADIUS),
        "params": params,
        "rows": rows,
        "mu_values": mu_values,
        "branch_ids": branch_ids,
        "analytic_lambda": analytic_lambda,
        "reference_lplus": reference_lplus,
        "reference_lminus": reference_lminus,
        "source_csv": source_csv,
    }


def get_case_data(beta_deg: float) -> dict[str, object]:
    csv_path = csv_path_for_beta(beta_deg)
    if csv_path.exists():
        rows = load_csv_rows(csv_path)
        return rows_to_case_data(beta_deg=beta_deg, rows=rows, source_csv=csv_path)

    rows = compute_case_rows(beta_deg=beta_deg)
    write_csv_rows(csv_path, rows)
    return rows_to_case_data(beta_deg=beta_deg, rows=rows, source_csv=csv_path)


def plot_beta_cases(cases: list[dict[str, object]], beta_values: tuple[float, ...], figure_path: Path) -> None:
    nrows, ncols = subplot_shape_for_case_count(len(cases))
    fig, axes = plt.subplots(  # keep the same per-panel scale as the existing 2x2 figure family
        nrows,
        ncols,
        figsize=(8.75 * ncols, 6.6 * nrows),
        squeeze=False,
    )
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    flat_axes = list(axes.ravel())

    for ax, case in zip(flat_axes, cases):
        mu_values = np.asarray(case["mu_values"], dtype=float)
        analytic_lambda = np.asarray(case["analytic_lambda"], dtype=float)
        reference_lplus = np.asarray(case["reference_lplus"], dtype=float)
        reference_lminus = np.asarray(case["reference_lminus"], dtype=float)
        beta_deg = float(case["beta"])

        plot_reference_curves(ax=ax, case=case, colors=colors)

        for branch_pos in range(min(N_PLOTTED_BRANCHES, analytic_lambda.shape[0])):
            color = colors[branch_pos % len(colors)]
            ax.plot(mu_values, analytic_lambda[branch_pos], linewidth=2.0, color=color, zorder=3)

        y_main = float(
            np.nanmax(
                np.concatenate(
                    (
                        analytic_lambda.ravel(),
                        reference_lplus.ravel(),
                        reference_lminus.ravel(),
                    )
                )
            )
        )
        ax.set_xlim(float(mu_values[0]), float(mu_values[-1]))
        ax.set_ylim(0.0, min(15.0, 1.08 * y_main))
        ax.set_title(f"beta = {beta_title(beta_deg)} deg, r = {FIXED_RADIUS:.3f}")
        ax.set_xlabel("mu")
        ax.set_ylabel(r"Dimensionless frequency parameter $\Lambda$")
        ax.grid(True, alpha=0.3)

    for ax in flat_axes[len(cases) :]:
        ax.set_visible(False)

    style_handles = [
        Line2D([0], [0], color="black", lw=2.0, label="solid - analytic"),
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
        "Analytic mu-sweeps at fixed radius r = 5 mm in Lambda\n"
        f"selected beta values: {format_beta_list(beta_values)} deg",
        y=0.985,
    )
    fig.legend(handles=style_handles, loc="upper center", ncol=3, fontsize=10, bbox_to_anchor=(0.5, 0.955))
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(figure_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> None:
    beta_values = parse_beta_values(argv)
    figure_path = figure_path_for_run()
    combined_csv_path = combined_csv_path_for_run()

    RESULTS_DIR.mkdir(exist_ok=True)
    cases = [get_case_data(beta_deg) for beta_deg in beta_values]
    combined_rows = build_combined_rows(cases)
    write_csv_rows(combined_csv_path, combined_rows)
    plot_beta_cases(cases=cases, beta_values=beta_values, figure_path=figure_path)

    print(f"fixed radius: {FIXED_RADIUS:.3f} m")
    print(f"requested betas (deg): {format_beta_list(beta_values)}")
    print(f"presentation mu base step: {PRESENTATION_MU_STEP:.3f}")
    print("local mu refinement windows: none")
    print(f"saved figure: {figure_path}")
    print(f"saved combined csv: {combined_csv_path}")
    for case in cases:
        print(f"beta={float(case['beta']):.1f} deg")
        print(f"  source csv: {case['source_csv']}")
        print(f"  shown branches: {', '.join(case['branch_ids'])}")


if __name__ == "__main__":
    main()
