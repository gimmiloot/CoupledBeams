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
    build_single_shape_plot,
    infer_title_label,
    radius_from_epsilon,
    resolve_single_output_path,
)


DEFAULT_BRANCH_NUMBER = 1
DEFAULT_MU = 0.0
DEFAULT_EPSILON = 0.0025
DEFAULT_BETA_DEG = 15.0


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
        help="Descendant branch number, e.g. 1, 2, 4, or 6. Default: 1.",
    )
    parser.add_argument(
        "--branch-id",
        default=None,
        help="Tracked branch id alternative, e.g. bending_desc_04.",
    )
    parser.add_argument("--mu", type=float, default=DEFAULT_MU, help="Target mu value.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Beam slenderness epsilon.")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA_DEG, help="Target beta angle in degrees.")
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
    return args


def main(argv: Sequence[str] | None = None) -> Path:
    args = parse_args(argv)
    branch_id = args.branch_id
    if branch_id is None:
        branch_number = args.branch_number if args.branch_number is not None else DEFAULT_BRANCH_NUMBER
        branch_id = branch_id_from_number(branch_number)

    radius = radius_from_epsilon(args.epsilon)
    title_label = args.title_label if args.title_label is not None else infer_title_label(branch_id)
    output_path = resolve_single_output_path(
        value=args.output,
        branch_id=branch_id,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
    )
    try:
        out_path = build_single_shape_plot(
            branch_id=branch_id,
            beta_deg=args.beta,
            mu_value=args.mu,
            radius=radius,
            epsilon=args.epsilon,
            output_path=output_path,
            title_label=title_label,
        )
    except KeyError as exc:
        raise SystemExit(str(exc)) from exc
    print(f"branch id: {branch_id}")
    print(f"epsilon: {args.epsilon:g}")
    print(f"radius: {radius:g} m ({radius * 1e3:g} mm)")
    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    print(f"saved PNG: {out_path}")
    return out_path


if __name__ == "__main__":
    main()
