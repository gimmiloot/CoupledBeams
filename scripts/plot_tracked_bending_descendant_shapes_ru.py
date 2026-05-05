from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib.tracked_bending_descendant_shapes import (  # noqa: E402
    ANALYSIS_BETA_STEP,
    ANALYSIS_MU_STEP,
    DEFAULT_BETA_DEG,
    DEFAULT_BRANCH_ID,
    DEFAULT_TARGET_MUS,
    DEFAULT_TARGET_RADII,
    build_single_shape_plot,
    collect_branch_rows,
    collect_branch_shape_cases,
    collect_single_branch_shape,
    default_output_path,
    default_single_output_path,
    draw_shape_case,
    filename_number_token,
    infer_title_label,
    lambda_from_frequency_hz,
    mu_label,
    radius_from_epsilon,
    resolve_output_path,
    resolve_single_output_path,
    shape_axis_limits,
    shape_case_from_tracking,
    shape_legend_handles,
    target_mus_token,
    track_branch_family,
)


def build_shape_plot(
    branch_id: str,
    beta_deg: float,
    target_mus: Sequence[float],
    radii: Sequence[float],
    output_path: Path,
    title_label: str,
    show_spectral_index: bool,
    spectral_index_label: str,
    axis_label_style: str,
) -> Path:
    shape_cases = collect_branch_shape_cases(
        branch_id=branch_id,
        beta_deg=beta_deg,
        target_mus=target_mus,
        radii=radii,
    )
    x_limits, y_limits = shape_axis_limits(shape_cases)
    fig, axes = plt.subplots(len(radii), len(target_mus), figsize=(15.2, 7.7), squeeze=False)
    radius_positions = {float(radius): idx for idx, radius in enumerate(radii)}
    mu_positions = {float(mu_value): idx for idx, mu_value in enumerate(target_mus)}
    x_axis_label = "координата x" if axis_label_style == "ru" else "x"
    y_axis_label = "координата y" if axis_label_style == "ru" else "y"

    for shape_case in shape_cases:
        radius = float(shape_case["r"])
        mu_value = float(shape_case["mu"])
        row_pos = radius_positions[radius]
        col_pos = mu_positions[mu_value]
        ax = axes[row_pos, col_pos]
        draw_shape_case(ax=ax, shape_case=shape_case, x_limits=x_limits, y_limits=y_limits)
        title = (
            f"β = {beta_deg:g}°, r = {radius * 1e3:.0f} мм, "
            f"μ = {mu_label(mu_value)}, Λ = {float(shape_case['lambda']):.3f}"
        )
        if show_spectral_index:
            title = f"{title}, {spectral_index_label}: {int(shape_case['current_sorted_index'])}"
        ax.set_title(title, fontsize=9.6 if show_spectral_index else 10.2)
        if row_pos == len(radii) - 1:
            ax.set_xlabel(x_axis_label)
        if col_pos == 0:
            ax.set_ylabel(y_axis_label)

    fig.legend(handles=shape_legend_handles(), loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965), fontsize=10)
    fig.suptitle(f"Формы колебаний: β = {beta_deg:g}°, {title_label}", y=0.992, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return output_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a Russian-labeled mode-shape plot for a tracked bending descendant branch."
    )
    parser.add_argument("--branch-id", default=DEFAULT_BRANCH_ID, help="Tracked branch id, e.g. bending_desc_01.")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA_DEG, help="Target beta angle in degrees.")
    parser.add_argument("--target-mus", nargs="+", type=float, default=DEFAULT_TARGET_MUS, help="Target mu values.")
    parser.add_argument(
        "--radii",
        nargs="+",
        type=float,
        default=tuple(float(value) for value in DEFAULT_TARGET_RADII),
        help="Target beam radii in meters.",
    )
    parser.add_argument("--output", default=None, help="Output PNG path. Relative paths are resolved from the repo root.")
    parser.add_argument("--title-label", default=None, help="Russian branch label for the figure title.")
    parser.add_argument("--no-spectral-index", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--spectral-index-label", default="место в спектре", help=argparse.SUPPRESS)
    parser.add_argument("--axis-label-style", choices=("ru", "xy"), default="ru", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    args.target_mus = tuple(float(value) for value in args.target_mus)
    args.radii = tuple(float(value) for value in args.radii)
    if not args.target_mus:
        parser.error("--target-mus must contain at least one value.")
    if not args.radii:
        parser.error("--radii must contain at least one value.")
    return args


def main(argv: Sequence[str] | None = None) -> Path:
    args = parse_args(argv)
    title_label = args.title_label if args.title_label is not None else infer_title_label(args.branch_id)
    output_path = resolve_output_path(
        value=args.output,
        branch_id=args.branch_id,
        beta_deg=args.beta,
        target_mus=args.target_mus,
    )
    out_path = build_shape_plot(
        branch_id=args.branch_id,
        beta_deg=args.beta,
        target_mus=args.target_mus,
        radii=args.radii,
        output_path=output_path,
        title_label=title_label,
        show_spectral_index=not args.no_spectral_index,
        spectral_index_label=args.spectral_index_label,
        axis_label_style=args.axis_label_style,
    )
    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    print(f"saved PNG: {out_path}")
    return out_path


if __name__ == "__main__":
    main()
