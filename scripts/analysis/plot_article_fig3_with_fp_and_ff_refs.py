from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.solvers import fixed_fixed_lambdas  # noqa: E402
from paper_dorofeev_style import generate_article_spectral_figures as article_figures  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
DEFAULT_OUTPUT = RESULTS_DIR / "article_fig3_with_fp_and_ff_refs.png"
DEFAULT_CSV_OUTPUT = RESULTS_DIR / "article_fig3_with_fp_and_ff_refs.csv"

FF_REFERENCE_LPLUS_WIDTH = 1.55
FF_REFERENCE_LPLUS_LINESTYLE = (0, (6.0, 2.0, 1.4, 2.0))
FF_REFERENCE_LPLUS_BLEND = 0.55
FF_REFERENCE_LMINUS_WIDTH = 1.45
FF_REFERENCE_LMINUS_LINESTYLE = ":"
FF_REFERENCE_LMINUS_BLEND = 0.55

LEGEND_FONT_SIZE = 10
LEGEND_LOC = "upper left"

CSV_FIELDNAMES = [
    "curve_group",
    "boundary_condition",
    "member_length_family",
    "branch_or_mode_index",
    "mu",
    "Lambda",
]


def resolve_repo_path(path: str | Path | None) -> Path | None:
    if path is None:
        return None
    resolved = Path(path)
    if resolved.is_absolute():
        return resolved
    return REPO_ROOT / resolved


def build_fixed_fixed_reference_curves(mu_values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    params = article_figures.build_params(article_figures.RADIUS)
    roots_ff = fixed_fixed_lambdas(max(article_figures.N_REFERENCE_LPLUS, article_figures.N_REFERENCE_LMINUS))
    l_base = params.L_base
    mu_values = np.asarray(mu_values, dtype=float)
    l_plus = l_base * (1.0 + mu_values)
    l_minus = l_base * (1.0 - mu_values)
    lambda_lplus = article_figures.single_lambda(roots_ff[: article_figures.N_REFERENCE_LPLUS], l_base, l_plus)
    lambda_lminus = article_figures.single_lambda(roots_ff[: article_figures.N_REFERENCE_LMINUS], l_base, l_minus)
    return lambda_lplus, lambda_lminus


def add_curve_rows(
    rows: list[dict[str, float | int | str]],
    *,
    curve_group: str,
    boundary_condition: str,
    member_length_family: str,
    mu_values: np.ndarray,
    values: np.ndarray,
) -> None:
    values = np.asarray(values, dtype=float)
    for row_idx in range(values.shape[0]):
        for col, mu in enumerate(mu_values):
            rows.append(
                {
                    "curve_group": curve_group,
                    "boundary_condition": boundary_condition,
                    "member_length_family": member_length_family,
                    "branch_or_mode_index": int(row_idx + 1),
                    "mu": float(mu),
                    "Lambda": float(values[row_idx, col]),
                }
            )


def write_curve_csv(
    path: Path,
    *,
    mu_values: np.ndarray,
    analytic: np.ndarray,
    fem: np.ndarray,
    reference_cs_lplus: np.ndarray,
    reference_cs_lminus: np.ndarray,
    reference_ff_lplus: np.ndarray,
    reference_ff_lminus: np.ndarray,
) -> None:
    rows: list[dict[str, float | int | str]] = []
    add_curve_rows(
        rows,
        curve_group="coupled_analytic",
        boundary_condition="coupled_clamped_arms",
        member_length_family="coupled_system",
        mu_values=mu_values,
        values=analytic,
    )
    add_curve_rows(
        rows,
        curve_group="coupled_fem",
        boundary_condition="coupled_clamped_arms",
        member_length_family="coupled_system",
        mu_values=mu_values,
        values=fem,
    )
    add_curve_rows(
        rows,
        curve_group="single_rod_reference",
        boundary_condition="clamped_pinned",
        member_length_family="L_plus",
        mu_values=mu_values,
        values=reference_cs_lplus,
    )
    add_curve_rows(
        rows,
        curve_group="single_rod_reference",
        boundary_condition="clamped_pinned",
        member_length_family="L_minus",
        mu_values=mu_values,
        values=reference_cs_lminus,
    )
    add_curve_rows(
        rows,
        curve_group="single_rod_reference",
        boundary_condition="clamped_clamped",
        member_length_family="L_plus",
        mu_values=mu_values,
        values=reference_ff_lplus,
    )
    add_curve_rows(
        rows,
        curve_group="single_rod_reference",
        boundary_condition="clamped_clamped",
        member_length_family="L_minus",
        mu_values=mu_values,
        values=reference_ff_lminus,
    )

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def add_reference_legend(ax: plt.Axes) -> None:
    handles = [
        Line2D([0], [0], color="black", lw=article_figures.FIGURE_3_ANALYTIC_LINE_WIDTH, ls="-", label="coupled branches"),
        Line2D(
            [0],
            [0],
            color="0.35",
            lw=article_figures.FIGURE_3_REFERENCE_LPLUS_WIDTH,
            ls=article_figures.FIGURE_3_REFERENCE_LPLUS_LINESTYLE,
            label="FP single rods, L+",
        ),
        Line2D(
            [0],
            [0],
            color="0.35",
            lw=article_figures.FIGURE_3_REFERENCE_LMINUS_WIDTH,
            ls=article_figures.FIGURE_3_REFERENCE_LMINUS_LINESTYLE,
            label="FP single rods, L-",
        ),
        Line2D(
            [0],
            [0],
            color="0.2",
            lw=FF_REFERENCE_LPLUS_WIDTH,
            ls=FF_REFERENCE_LPLUS_LINESTYLE,
            label="FF single rods, L+",
        ),
        Line2D(
            [0],
            [0],
            color="0.2",
            lw=FF_REFERENCE_LMINUS_WIDTH,
            ls=FF_REFERENCE_LMINUS_LINESTYLE,
            label="FF single rods, L-",
        ),
    ]
    ax.legend(handles=handles, loc=LEGEND_LOC, frameon=False, fontsize=LEGEND_FONT_SIZE, handlelength=3.6)


def render_diagnostic_figure(output_path: Path, csv_output_path: Path | None) -> dict[str, Path | int]:
    mu_values, analytic, fem, reference_cs_lplus, reference_cs_lminus = article_figures.build_mu_tracking_case()
    reference_ff_lplus, reference_ff_lminus = build_fixed_fixed_reference_curves(mu_values)

    article_figures.configure_matplotlib()
    fig, ax = plt.subplots(figsize=article_figures.FIGURE_3_FIGSIZE)

    colors: list[str] = []
    for branch_idx in range(article_figures.N_BRANCHES_FIGURE_3):
        (line,) = ax.plot(
            mu_values,
            analytic[branch_idx],
            linewidth=article_figures.FIGURE_3_ANALYTIC_LINE_WIDTH,
            zorder=3,
            rasterized=False,
        )
        colors.append(line.get_color())

    for idx in range(reference_cs_lplus.shape[0]):
        ax.plot(
            mu_values,
            reference_cs_lplus[idx],
            linewidth=article_figures.FIGURE_3_REFERENCE_LPLUS_WIDTH,
            linestyle=article_figures.FIGURE_3_REFERENCE_LPLUS_LINESTYLE,
            color=article_figures.blend_with_white(colors[idx % len(colors)], article_figures.FIGURE_3_REFERENCE_LPLUS_BLEND),
            zorder=1,
            rasterized=False,
        )
    for idx in range(reference_cs_lminus.shape[0]):
        ax.plot(
            mu_values,
            reference_cs_lminus[idx],
            linewidth=article_figures.FIGURE_3_REFERENCE_LMINUS_WIDTH,
            linestyle=article_figures.FIGURE_3_REFERENCE_LMINUS_LINESTYLE,
            color=article_figures.blend_with_white(colors[idx % len(colors)], article_figures.FIGURE_3_REFERENCE_LMINUS_BLEND),
            zorder=1,
            rasterized=False,
        )

    for idx in range(reference_ff_lplus.shape[0]):
        ax.plot(
            mu_values,
            reference_ff_lplus[idx],
            linewidth=FF_REFERENCE_LPLUS_WIDTH,
            linestyle=FF_REFERENCE_LPLUS_LINESTYLE,
            color=article_figures.blend_with_white(colors[idx % len(colors)], FF_REFERENCE_LPLUS_BLEND),
            zorder=2,
            rasterized=False,
        )
    for idx in range(reference_ff_lminus.shape[0]):
        ax.plot(
            mu_values,
            reference_ff_lminus[idx],
            linewidth=FF_REFERENCE_LMINUS_WIDTH,
            linestyle=FF_REFERENCE_LMINUS_LINESTYLE,
            color=article_figures.blend_with_white(colors[idx % len(colors)], FF_REFERENCE_LMINUS_BLEND),
            zorder=2,
            rasterized=False,
        )

    for branch_idx, color in enumerate(colors):
        ax.plot(
            mu_values,
            fem[branch_idx],
            linestyle="None",
            marker="o",
            markersize=article_figures.FIGURE_3_FEM_MARKER_SIZE,
            color=color,
            markerfacecolor=color,
            markeredgecolor="none",
            markeredgewidth=article_figures.FIGURE_3_FEM_MARKER_EDGE_WIDTH,
            alpha=1.0,
            zorder=4,
            rasterized=False,
        )

    article_figures.style_axes(
        ax,
        xlabel=article_figures.FIGURE_3_XLABEL,
        ylabel=article_figures.FIGURE_3_YLABEL,
        xlabel_font_size=article_figures.FIGURE_3_XLABEL_FONT_SIZE,
        ylabel_font_size=article_figures.FIGURE_3_YLABEL_FONT_SIZE,
        tick_label_size=article_figures.FIGURE_3_TICK_LABEL_SIZE,
        tick_width=article_figures.FIGURE_3_TICK_WIDTH,
        tick_length=article_figures.FIGURE_3_TICK_LENGTH,
        tick_pad=article_figures.FIGURE_3_TICK_PAD,
        title_text=article_figures.FIGURE_3_TITLE_TEXT,
        title_font_size=article_figures.FIGURE_3_TITLE_FONT_SIZE,
        spine_width=article_figures.FIGURE_3_SPINE_WIDTH,
    )
    ax.set_xlim(float(mu_values[0]), float(mu_values[-1]))
    ax.set_xticks(np.arange(0.0, 1.0, 0.1))
    y_main = float(np.max(np.concatenate((analytic.ravel(), fem.ravel()))))
    ax.set_ylim(0.0, 1.08 * y_main)

    article_figures.add_figure_3_internal_mode_labels(
        ax=ax,
        annotations=article_figures.FIGURE_3_INTERNAL_MODE_LABELS,
    )

    article_figures.add_right_labels(
        ax=ax,
        x_value=float(mu_values[-1]),
        y_values=analytic[:, -1],
        labels=[rf"$\mathit{{{index}}}$" for index in range(1, article_figures.N_BRANCHES_FIGURE_3 + 1)],
        colors=["black"] * article_figures.N_BRANCHES_FIGURE_3,
        x_offset=article_figures.FIGURE_3_RIGHT_LABEL_X_OFFSET,
        font_size=article_figures.FIGURE_3_RIGHT_LABEL_FONT_SIZE,
        font_weight=article_figures.FIGURE_3_RIGHT_LABEL_FONT_WEIGHT,
    )

    add_reference_legend(ax)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    save_format = output_path.suffix.lstrip(".").lower() or None
    fig.savefig(output_path, format=save_format, dpi=article_figures.FIGURE_3_DPI, bbox_inches="tight")
    plt.close(fig)

    if csv_output_path is not None:
        write_curve_csv(
            csv_output_path,
            mu_values=mu_values,
            analytic=analytic,
            fem=fem,
            reference_cs_lplus=reference_cs_lplus,
            reference_cs_lminus=reference_cs_lminus,
            reference_ff_lplus=reference_ff_lplus,
            reference_ff_lminus=reference_ff_lminus,
        )

    return {
        "output": output_path,
        "csv_output": csv_output_path or Path(),
        "mu_points": int(len(mu_values)),
        "fp_lplus_curves": int(reference_cs_lplus.shape[0]),
        "fp_lminus_curves": int(reference_cs_lminus.shape[0]),
        "ff_lplus_curves": int(reference_ff_lplus.shape[0]),
        "ff_lminus_curves": int(reference_ff_lminus.shape[0]),
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a results-only diagnostic version of article Figure 3 with both "
            "clamped-pinned and clamped-clamped single-rod dashed references."
        ),
    )
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Figure output path.")
    parser.add_argument("--csv-output", default=DEFAULT_CSV_OUTPUT, help="Curve CSV output path; use empty string to skip.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> dict[str, Path | int]:
    args = parse_args(argv)
    output_path = resolve_repo_path(args.output)
    csv_output_path = resolve_repo_path(args.csv_output) if str(args.csv_output).strip() else None
    if output_path is None:
        raise RuntimeError("Output path could not be resolved.")

    result = render_diagnostic_figure(output_path=output_path, csv_output_path=csv_output_path)
    print(f"saved figure: {result['output']}")
    if csv_output_path is not None:
        print(f"saved CSV: {result['csv_output']}")
    print(f"mu points: {result['mu_points']}")
    print(f"FP references: L+={result['fp_lplus_curves']}, L-={result['fp_lminus_curves']}")
    print(f"FF references: L+={result['ff_lplus_curves']}, L-={result['ff_lminus_curves']}")
    print("article figure files were not written")
    return result


if __name__ == "__main__":
    main()
