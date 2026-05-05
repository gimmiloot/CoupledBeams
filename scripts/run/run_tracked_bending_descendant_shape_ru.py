from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence
import warnings


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


DEFAULT_BRANCH_NUMBER = 1
DEFAULT_MU = 0.0
DEFAULT_EPSILON = 0.0025
DEFAULT_BETA_DEG = 15.0
DEFAULT_L_TOTAL = 2.0
DEFAULT_DPI = 240


def branch_id_from_number(branch_number: int) -> str:
    return f"bending_desc_{branch_number:02d}"


def branch_number_from_id(branch_id: str) -> int | None:
    prefix = "bending_desc_"
    if branch_id.startswith(prefix):
        suffix = branch_id.removeprefix(prefix)
        if suffix.isdecimal():
            return int(suffix)
    return None


def resolve_repo_path(value: str) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build one Russian-labeled mode-shape plot for a tracked bending descendant branch."
    )
    parser.add_argument(
        "--branch-number",
        type=int,
        default=None,
        help=(
            f"Descendant branch number used to form bending_desc_NN, not the current sorted spectral index. "
            f"Default: {DEFAULT_BRANCH_NUMBER}."
        ),
    )
    parser.add_argument(
        "--branch-id",
        default=None,
        help="Tracked branch id alternative to --branch-number, e.g. bending_desc_04.",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=DEFAULT_MU,
        help="Target non-negative length-asymmetry parameter mu; values above 0.9 print a warning.",
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Thickness/slenderness parameter epsilon.")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA_DEG, help="Target coupling angle in degrees.")
    parser.add_argument(
        "--l-total",
        type=float,
        default=DEFAULT_L_TOTAL,
        help="Total beam length used only for epsilon -> radius conversion. Default: 2.0.",
    )
    parser.add_argument(
        "--plot-kind",
        choices=PLOT_KINDS,
        default="full",
        help=(
            "Plot mode: full = global modal displacement; transverse = local bending deflection only; "
            "components = local axial/transverse diagnostic curves."
        ),
    )
    parser.add_argument(
        "--mode-scale",
        type=float,
        default=MODE_SHAPE_SCALE,
        help="Visual-only deformation scale for full/transverse geometry plots; does not change eigenvectors or frequencies.",
    )
    parser.add_argument(
        "--normalize",
        choices=NORMALIZE_KINDS,
        default="auto",
        help="Normalization mode. auto resolves to max-full for full/components and max-transverse for transverse.",
    )
    parser.add_argument("--output", default=None, help="Output PNG path. Relative paths are resolved from the repo root.")
    parser.add_argument("--title-label", default=None, help="Russian branch label for the figure title.")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="Output PNG DPI.")
    parser.add_argument(
        "--figsize",
        nargs=2,
        type=float,
        metavar=("WIDTH", "HEIGHT"),
        default=None,
        help="Optional matplotlib figure size in inches.",
    )
    parser.add_argument("--show", action="store_true", help="Display the figure after saving, using the active backend.")
    parser.add_argument(
        "--print-diagnostics",
        dest="print_diagnostics",
        action="store_true",
        default=True,
        help="Print one-case diagnostics to stdout. Default: enabled.",
    )
    parser.add_argument(
        "--no-print-diagnostics",
        dest="print_diagnostics",
        action="store_false",
        help="Suppress diagnostic stdout, except for saved output paths.",
    )
    parser.add_argument(
        "--save-diagnostics-csv",
        default=None,
        help="Optional one-row diagnostics CSV path. Relative paths are resolved from the repo root.",
    )
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
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if args.mode_scale <= 0.0:
        parser.error("--mode-scale must be positive.")
    if args.dpi <= 0:
        parser.error("--dpi must be positive.")
    if args.figsize is not None and any(value <= 0.0 for value in args.figsize):
        parser.error("--figsize values must be positive.")
    return args


def build_diagnostics_row(
    *,
    diagnostics: dict[str, float | int | str],
    branch_number: int | None,
    epsilon: float,
    l_total: float,
    radius: float,
    normalize_requested: str,
    normalize_resolved: str,
    output_path: Path,
) -> dict[str, float | int | str]:
    axial_fraction = float(diagnostics["axial_fraction"])
    return {
        "branch_id": str(diagnostics["branch_id"]),
        "branch_number": "" if branch_number is None else int(branch_number),
        "beta": float(diagnostics["beta"]),
        "mu": float(diagnostics["mu"]),
        "epsilon": float(epsilon),
        "l_total": float(l_total),
        "radius": float(radius),
        "lambda": float(diagnostics["lambda"]),
        "current_sorted_index": int(diagnostics["current_sorted_index"]),
        "plot_kind": str(diagnostics["plot_kind"]),
        "mode_scale": float(diagnostics["mode_scale"]),
        "normalize_requested": normalize_requested,
        "normalize_resolved": normalize_resolved,
        "output_path": str(output_path),
        "max_full_displacement": float(diagnostics["max_full_displacement"]),
        "max_transverse_displacement": float(diagnostics["max_abs_transverse"]),
        "max_axial_displacement": float(diagnostics["max_abs_axial"]),
        "max_abs_transverse_over_max_full": float(diagnostics["transverse_to_full_ratio"]),
        "max_abs_axial_over_max_full": float(diagnostics["axial_to_full_ratio"]),
        "axial_fraction": "" if axial_fraction != axial_fraction else axial_fraction,
        "tracking_mac_diagnostics": "not exported",
    }


def print_diagnostics(row: dict[str, float | int | str]) -> None:
    print(f"branch id: {row['branch_id']}")
    print(f"branch number: {row['branch_number'] if row['branch_number'] != '' else 'not available'}")
    print(f"beta: {float(row['beta']):g} deg")
    print(f"mu: {float(row['mu']):g}")
    print(f"epsilon: {float(row['epsilon']):g}")
    print(f"l_total: {float(row['l_total']):g}")
    print(f"radius: {float(row['radius']):g} m ({float(row['radius']) * 1e3:g} mm)")
    print(f"Lambda: {float(row['lambda']):.6g}")
    print(f"current sorted index: {int(row['current_sorted_index'])}")
    print(f"plot kind: {row['plot_kind']}")
    print(f"mode scale: {float(row['mode_scale']):g}")
    print(f"normalize requested: {row['normalize_requested']}")
    print(f"normalize resolved: {row['normalize_resolved']}")
    print(f"output path: {row['output_path']}")
    print(f"max full displacement norm before normalization: {float(row['max_full_displacement']):.6e}")
    print(f"max transverse displacement before normalization: {float(row['max_transverse_displacement']):.6e}")
    print(f"max axial displacement before normalization: {float(row['max_axial_displacement']):.6e}")
    print(f"max_abs_transverse / max_full: {float(row['max_abs_transverse_over_max_full']):.6g}")
    print(f"max_abs_axial / max_full: {float(row['max_abs_axial_over_max_full']):.6g}")
    print(f"axial_fraction: {row['axial_fraction'] if row['axial_fraction'] != '' else 'not available'}")
    print(f"tracking MAC diagnostics: {row['tracking_mac_diagnostics']}")


def write_diagnostics_csv(path: Path, row: dict[str, float | int | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def main(argv: Sequence[str] | None = None) -> Path:
    args = parse_args(argv)
    if args.mu > 0.9:
        warnings.warn("--mu is outside the usual [0, 0.9] analysis window; continuing with requested value.")

    branch_id = args.branch_id
    if branch_id is None:
        branch_number = args.branch_number if args.branch_number is not None else DEFAULT_BRANCH_NUMBER
        branch_id = branch_id_from_number(branch_number)
    else:
        branch_number = branch_number_from_id(branch_id)

    radius = radius_from_epsilon(args.epsilon, l_total=args.l_total)
    title_label = args.title_label if args.title_label is not None else infer_title_label(branch_id)
    normalization = resolve_normalization(plot_kind=args.plot_kind, normalize=args.normalize)
    output_path = resolve_single_output_path(
        value=args.output,
        branch_id=branch_id,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
        plot_kind=args.plot_kind,
        mode_scale=args.mode_scale,
        normalize=args.normalize,
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
            dpi=args.dpi,
            figsize=None if args.figsize is None else tuple(args.figsize),
            show=args.show,
        )
    except KeyError as exc:
        raise SystemExit(str(exc)) from exc

    diagnostics = single_shape_diagnostics(
        shape_case,
        plot_kind=args.plot_kind,
        mode_scale=args.mode_scale,
        normalize=normalization,
    )
    row = build_diagnostics_row(
        diagnostics=diagnostics,
        branch_number=branch_number,
        epsilon=args.epsilon,
        l_total=args.l_total,
        radius=radius,
        normalize_requested=args.normalize,
        normalize_resolved=normalization,
        output_path=out_path,
    )

    if args.print_diagnostics:
        print_diagnostics(row)
        if normalization == "max-transverse" and float(row["max_transverse_displacement"]) <= NEAR_ZERO_NORM:
            print("warning: max transverse displacement is near zero; normalization denominator was clamped")
        print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
        print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
        print("local refinement windows: none")
    else:
        print(f"saved PNG: {out_path}")

    if args.save_diagnostics_csv is not None:
        diagnostics_path = resolve_repo_path(args.save_diagnostics_csv)
        write_diagnostics_csv(diagnostics_path, row)
        print(f"saved diagnostics CSV: {diagnostics_path}")

    if args.print_diagnostics:
        print(f"saved PNG: {out_path}")
    return out_path


if __name__ == "__main__":
    main()
