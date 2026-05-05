from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib.tracked_bending_descendant_shapes import (  # noqa: E402
    ANALYSIS_BETA_STEP,
    ANALYSIS_MU_STEP,
    MODE_SHAPE_SCALE,
    NEAR_ZERO_NORM,
    NORMALIZE_KINDS,
    PLOT_KINDS,
    build_single_shape_plot,
    collect_single_branch_shape,
    infer_title_label,
    radius_from_epsilon,
    resolve_normalization,
    resolve_single_output_path,
    single_shape_diagnostics,
)


DEFAULT_BRANCH_NUMBER = 5
DEFAULT_MU = 0.0
DEFAULT_EPSILON = 0.0025
DEFAULT_BETA_DEG = 90.0


def branch_id_from_number(branch_number: int) -> str:
    return f"bending_desc_{branch_number:02d}"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build one Russian-labeled mode-shape plot for a tracked bending descendant branch."
    )
    parser.add_argument(
        "--branch-number",
        type=int,
        default=None,
        help=f"Descendant branch number, e.g. 1, 2, 4, or 6. Default: {DEFAULT_BRANCH_NUMBER}.",
    )
    parser.add_argument(
        "--branch-id",
        default=None,
        help="Tracked branch id alternative, e.g. bending_desc_04.",
    )
    parser.add_argument("--mu", type=float, default=DEFAULT_MU, help="Target mu value.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Beam slenderness epsilon.")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA_DEG, help="Target beta angle in degrees.")
    parser.add_argument(
        "--plot-kind",
        choices=PLOT_KINDS,
        default="full",
        help="Plot mode: full modal displacement, local transverse deflection, or local component diagnostics.",
    )
    parser.add_argument(
        "--mode-scale",
        type=float,
        default=MODE_SHAPE_SCALE,
        help="Visual deformation scale for full/transverse geometry plots.",
    )
    parser.add_argument(
        "--normalize",
        choices=NORMALIZE_KINDS,
        default=None,
        help="Normalization for modal postprocessing. Defaults depend on --plot-kind.",
    )
    parser.add_argument("--output", default=None, help="Output PNG path. Relative paths are resolved from the repo root.")
    parser.add_argument("--title-label", default=None, help="Russian branch label for the figure title.")
    args = parser.parse_args(argv)

    if args.branch_number is not None and args.branch_id is not None:
        parser.error("--branch-number and --branch-id are alternatives; provide only one.")
    if args.branch_number is not None and args.branch_number <= 0:
        parser.error("--branch-number must be positive.")
    if args.mu < 0.0:
        parser.error("--mu must be non-negative; tracking starts at mu=0.")
    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.beta < 0.0:
        parser.error("--beta must be non-negative; tracking starts at beta=0.")
    if args.mode_scale <= 0.0:
        parser.error("--mode-scale must be positive.")
    return args


def main(argv: Sequence[str] | None = None) -> Path:
    args = parse_args(argv)
    branch_id = args.branch_id
    if branch_id is None:
        branch_number = args.branch_number if args.branch_number is not None else DEFAULT_BRANCH_NUMBER
        branch_id = branch_id_from_number(branch_number)

    radius = radius_from_epsilon(args.epsilon)
    title_label = args.title_label if args.title_label is not None else infer_title_label(branch_id)
    normalization = resolve_normalization(plot_kind=args.plot_kind, normalize=args.normalize)
    output_path = resolve_single_output_path(
        value=args.output,
        branch_id=branch_id,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
        plot_kind=args.plot_kind,
    )
    try:
        shape_case = collect_single_branch_shape(
            branch_id=branch_id,
            beta_deg=args.beta,
            mu_value=args.mu,
            radius=radius,
        )
        out_path = build_single_shape_plot(
            branch_id=branch_id,
            beta_deg=args.beta,
            mu_value=args.mu,
            radius=radius,
            epsilon=args.epsilon,
            output_path=output_path,
            title_label=title_label,
            plot_kind=args.plot_kind,
            mode_scale=args.mode_scale,
            normalize=normalization,
            shape_case=shape_case,
        )
    except KeyError as exc:
        raise SystemExit(str(exc)) from exc
    diagnostics = single_shape_diagnostics(
        shape_case,
        plot_kind=args.plot_kind,
        mode_scale=args.mode_scale,
        normalize=normalization,
    )
    print(f"branch id: {branch_id}")
    print(f"beta: {args.beta:g} deg")
    print(f"mu: {args.mu:g}")
    print(f"epsilon: {args.epsilon:g}")
    print(f"radius: {radius:g} m ({radius * 1e3:g} mm)")
    print(f"Lambda: {float(diagnostics['lambda']):.6g}")
    print(f"current sorted index: {int(diagnostics['current_sorted_index'])}")
    print(f"plot kind: {args.plot_kind}")
    print(f"mode scale: {args.mode_scale:g}")
    print(f"normalization: {normalization}")
    print(f"max full displacement norm before normalization: {float(diagnostics['max_full_displacement']):.6e}")
    print(f"max transverse displacement before normalization: {float(diagnostics['max_abs_transverse']):.6e}")
    print(f"max axial displacement before normalization: {float(diagnostics['max_abs_axial']):.6e}")
    print(f"max_abs_transverse / max_full: {float(diagnostics['transverse_to_full_ratio']):.6g}")
    print(f"max_abs_axial / max_full: {float(diagnostics['axial_to_full_ratio']):.6g}")
    axial_fraction = float(diagnostics["axial_fraction"])
    if axial_fraction == axial_fraction:
        print(f"axial_fraction: {axial_fraction:.6g}")
    else:
        print("axial_fraction: not available")
    if normalization == "max-transverse" and float(diagnostics["max_abs_transverse"]) <= NEAR_ZERO_NORM:
        print("warning: max transverse displacement is near zero; normalization denominator was clamped")
    print("tracking MAC diagnostics: not exported by current helper")
    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    print(f"saved PNG: {out_path}")
    return out_path


if __name__ == "__main__":
    main()
