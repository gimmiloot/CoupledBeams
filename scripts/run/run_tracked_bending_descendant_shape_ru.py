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
    DIAGNOSTICS_LEVELS,
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


# ============================================================
# USER PARAMETERS
# Change these values for ordinary runs.
# CLI arguments, if provided, override these defaults.
# ============================================================

DEFAULT_BRANCH_NUMBER = 5
DEFAULT_BRANCH_ID = None

DEFAULT_MU = 0.7
DEFAULT_EPSILON = 0.01
DEFAULT_BETA = 15.0
DEFAULT_L_TOTAL = 2.0

DEFAULT_PLOT_KIND = "full"
DEFAULT_MODE_SCALE = 0.05
DEFAULT_NORMALIZE = "auto"

DEFAULT_OUTPUT = None
DEFAULT_TITLE_LABEL = None
DEFAULT_DPI = 240
DEFAULT_FIGSIZE = None
DEFAULT_SHOW = False
DEFAULT_PRINT_DIAGNOSTICS = True
DEFAULT_SAVE_DIAGNOSTICS_CSV = None
DEFAULT_DIAGNOSTICS_LEVEL = "basic"

PARAMETER_SOURCE_LABEL = "USER PARAMETERS + CLI overrides"
PROJECTION_RELATIVE_WARNING_THRESHOLD = 1e-10

EXTENDED_DIAGNOSTIC_KEYS = (
    "left_projection_max_abs_u_error",
    "left_projection_max_abs_v_error",
    "left_projection_max_norm_error",
    "left_projection_relative_norm_error",
    "right_projection_max_abs_u_error",
    "right_projection_max_abs_v_error",
    "right_projection_max_norm_error",
    "right_projection_relative_norm_error",
    "max_abs_du_left_ds",
    "max_abs_du_right_ds",
    "rms_du_left_ds",
    "rms_du_right_ds",
    "max_abs_dw_left_ds",
    "max_abs_dw_right_ds",
    "rms_dw_left_ds",
    "rms_dw_right_ds",
    "left_axial_energy",
    "right_axial_energy",
    "left_bending_energy",
    "right_bending_energy",
    "total_axial_energy",
    "total_bending_energy",
    "total_energy",
    "axial_energy_fraction_from_arm_energies",
    "right_axial_share",
    "right_bending_share",
)


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


def option_was_provided(argv: Sequence[str], option_name: str) -> bool:
    return any(arg == option_name or arg.startswith(f"{option_name}=") for arg in argv)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    raw_args = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Build one Russian-labeled mode-shape plot for a tracked bending descendant branch. "
            "Defaults are defined in the USER PARAMETERS block at the top of this file; "
            "CLI arguments override them."
        )
    )
    branch_number_explicit = option_was_provided(raw_args, "--branch-number")
    branch_id_explicit = option_was_provided(raw_args, "--branch-id")
    parser.add_argument(
        "--branch-number",
        type=int,
        default=DEFAULT_BRANCH_NUMBER,
        help=(
            f"Descendant branch number used to form bending_desc_NN, not the current sorted spectral index. "
            f"Default: {DEFAULT_BRANCH_NUMBER}."
        ),
    )
    parser.add_argument(
        "--branch-id",
        default=DEFAULT_BRANCH_ID,
        help=f"Tracked branch id alternative to --branch-number, e.g. bending_desc_04. Default: {DEFAULT_BRANCH_ID}.",
    )
    parser.add_argument(
        "--mu",
        type=float,
        default=DEFAULT_MU,
        help="Target non-negative length-asymmetry parameter mu; values above 0.9 print a warning.",
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Thickness/slenderness parameter epsilon.")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA, help="Target coupling angle in degrees.")
    parser.add_argument(
        "--l-total",
        type=float,
        default=DEFAULT_L_TOTAL,
        help="Total beam length used only for epsilon -> radius conversion. Default: 2.0.",
    )
    parser.add_argument(
        "--plot-kind",
        choices=PLOT_KINDS,
        default=DEFAULT_PLOT_KIND,
        help=(
            "Plot mode: full = global modal displacement; transverse = local bending deflection only; "
            "components = local axial/transverse diagnostic curves."
        ),
    )
    parser.add_argument(
        "--mode-scale",
        type=float,
        default=DEFAULT_MODE_SCALE,
        help="Visual-only deformation scale for full/transverse geometry plots; does not change eigenvectors or frequencies.",
    )
    parser.add_argument(
        "--normalize",
        choices=NORMALIZE_KINDS,
        default=DEFAULT_NORMALIZE,
        help="Normalization mode. auto resolves to max-full for full/components and max-transverse for transverse.",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help="Output PNG path. Relative paths are resolved from the repo root.",
    )
    parser.add_argument("--title-label", default=DEFAULT_TITLE_LABEL, help="Russian branch label for the figure title.")
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI, help="Output PNG DPI.")
    parser.add_argument(
        "--figsize",
        nargs=2,
        type=float,
        metavar=("WIDTH", "HEIGHT"),
        default=DEFAULT_FIGSIZE,
        help="Optional matplotlib figure size in inches.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        default=DEFAULT_SHOW,
        help="Display the figure after saving, using the active backend.",
    )
    parser.add_argument(
        "--print-diagnostics",
        dest="print_diagnostics",
        action="store_true",
        default=DEFAULT_PRINT_DIAGNOSTICS,
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
        default=DEFAULT_SAVE_DIAGNOSTICS_CSV,
        help="Optional one-row diagnostics CSV path. Relative paths are resolved from the repo root.",
    )
    parser.add_argument(
        "--diagnostics-level",
        choices=DIAGNOSTICS_LEVELS,
        default=DEFAULT_DIAGNOSTICS_LEVEL,
        help="Diagnostic detail level: basic keeps the standard output; all adds projection, derivative, and arm-energy checks.",
    )
    args = parser.parse_args(raw_args)

    if branch_number_explicit and branch_id_explicit:
        parser.error("--branch-number and --branch-id are alternatives; provide only one.")
    if branch_number_explicit and not branch_id_explicit:
        args.branch_id = None
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
    parameter_source: str,
    diagnostics_level: str,
) -> dict[str, float | int | str]:
    axial_fraction = float(diagnostics["axial_fraction"])
    row: dict[str, float | int | str] = {
        "parameter_source": parameter_source,
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
    if diagnostics_level == "all":
        row["diagnostics_level"] = diagnostics_level
        for key in EXTENDED_DIAGNOSTIC_KEYS:
            row[key] = diagnostics[key]
    return row


def print_diagnostics(row: dict[str, float | int | str]) -> None:
    print(f"parameter_source: {row['parameter_source']}")
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


def print_extended_diagnostics(row: dict[str, float | int | str]) -> None:
    print("diagnostics level: all")
    print("projection reconstruction check:")
    print(f"  left max |u error|: {float(row['left_projection_max_abs_u_error']):.6e}")
    print(f"  left max |v error|: {float(row['left_projection_max_abs_v_error']):.6e}")
    print(f"  left max norm error: {float(row['left_projection_max_norm_error']):.6e}")
    print(f"  left relative norm error: {float(row['left_projection_relative_norm_error']):.6e}")
    print(f"  right max |u error|: {float(row['right_projection_max_abs_u_error']):.6e}")
    print(f"  right max |v error|: {float(row['right_projection_max_abs_v_error']):.6e}")
    print(f"  right max norm error: {float(row['right_projection_max_norm_error']):.6e}")
    print(f"  right relative norm error: {float(row['right_projection_relative_norm_error']):.6e}")
    if (
        float(row["left_projection_relative_norm_error"]) > PROJECTION_RELATIVE_WARNING_THRESHOLD
        or float(row["right_projection_relative_norm_error"]) > PROJECTION_RELATIVE_WARNING_THRESHOLD
    ):
        print("warning: local projection reconstruction relative error exceeds 1e-10")

    print("local derivative diagnostics (nodal amplitude slopes, not exact FEM strain fields):")
    print(f"  max |du_left/ds|: {float(row['max_abs_du_left_ds']):.6e}")
    print(f"  max |du_right/ds|: {float(row['max_abs_du_right_ds']):.6e}")
    print(f"  rms du_left/ds: {float(row['rms_du_left_ds']):.6e}")
    print(f"  rms du_right/ds: {float(row['rms_du_right_ds']):.6e}")
    print(f"  max |dw_left/ds|: {float(row['max_abs_dw_left_ds']):.6e}")
    print(f"  max |dw_right/ds|: {float(row['max_abs_dw_right_ds']):.6e}")
    print(f"  rms dw_left/ds: {float(row['rms_dw_left_ds']):.6e}")
    print(f"  rms dw_right/ds: {float(row['rms_dw_right_ds']):.6e}")

    print("arm-wise diagnostic energies:")
    print(f"  left axial energy: {float(row['left_axial_energy']):.6e}")
    print(f"  right axial energy: {float(row['right_axial_energy']):.6e}")
    print(f"  left bending energy: {float(row['left_bending_energy']):.6e}")
    print(f"  right bending energy: {float(row['right_bending_energy']):.6e}")
    print(f"  total axial energy: {float(row['total_axial_energy']):.6e}")
    print(f"  total bending energy: {float(row['total_bending_energy']):.6e}")
    print(f"  total energy: {float(row['total_energy']):.6e}")
    print(f"  axial energy fraction from arm energies: {float(row['axial_energy_fraction_from_arm_energies']):.6e}")
    print(f"  right axial share: {float(row['right_axial_share']):.6e}")
    print(f"  right bending share: {float(row['right_bending_share']):.6e}")


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
        diagnostics_level=args.diagnostics_level,
        epsilon=args.epsilon,
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
        parameter_source=PARAMETER_SOURCE_LABEL,
        diagnostics_level=args.diagnostics_level,
    )

    if args.print_diagnostics:
        print_diagnostics(row)
        if args.diagnostics_level == "all":
            print_extended_diagnostics(row)
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
