from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.analysis.compare_analytic_fem_tracked_descendant_shape import (  # noqa: E402
    filename_number_token,
    main as compare_shape_main,
)


DEFAULT_BRANCH_IDS = ("bending_desc_01", "bending_desc_02")
DEFAULT_BETA = 15.0
DEFAULT_EPSILON = 0.0025
DEFAULT_MUS = (0.0, 0.1, 0.2)
DEFAULT_L_TOTAL = 2.0
DEFAULT_NUM_ROOTS = 20
DEFAULT_ROOT_WINDOW = 1.0
DEFAULT_OUTPUT = REPO_ROOT / "results" / "article_shape_validation_beta15_summary.csv"
DEFAULT_OUTPUT_PREFIX_ROOT = REPO_ROOT / "results" / "article_shape_validation"

SUMMARY_FIELDNAMES = [
    "branch_id",
    "mu",
    "epsilon",
    "beta",
    "fem_lambda",
    "analytic_lambda",
    "rel_lambda_diff",
    "analytic_root_index",
    "component_rel_l2_all",
    "rel_l2_w_left",
    "rel_l2_w_right",
    "rel_l2_u_left",
    "rel_l2_u_right",
    "analytic_axial_energy_fraction",
    "fem_axial_energy_fraction",
    "analytic_right_axial_share",
    "fem_right_axial_share",
    "conclusion_flag",
]

LAMBDA_MATCH_TOL = 1e-5
COMPONENT_REL_L2_TOL = 0.25
COMPONENT_PANEL_REL_L2_TOL = 0.35
ENERGY_FRACTION_TOL = 0.10
RIGHT_AXIAL_SHARE_TOL = 0.25


def resolve_repo_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Run analytic-vs-FEM tracked-descendant shape validation for the beta=15 article cases. "
            "This is diagnostic postprocessing only."
        ),
    )
    parser.add_argument(
        "--branch-ids",
        nargs="+",
        default=list(DEFAULT_BRANCH_IDS),
        help="Tracked FEM descendant branch ids. Default: bending_desc_01 bending_desc_02.",
    )
    parser.add_argument(
        "--mus",
        nargs="+",
        type=float,
        default=list(DEFAULT_MUS),
        help="Mu values to validate. Default: 0 0.1 0.2.",
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA, help="Coupling angle in degrees.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Slenderness epsilon.")
    parser.add_argument(
        "--l-total",
        type=float,
        default=DEFAULT_L_TOTAL,
        help="Total beam length used only for epsilon-to-radius conversion.",
    )
    parser.add_argument(
        "--num-roots",
        type=int,
        default=DEFAULT_NUM_ROOTS,
        help="Number of analytic roots acquired by the underlying comparison script.",
    )
    parser.add_argument(
        "--root-window",
        type=float,
        default=DEFAULT_ROOT_WINDOW,
        help="Local determinant scan half-window around the FEM Lambda.",
    )
    parser.add_argument(
        "--right-coordinate",
        choices=("external-to-joint", "joint-to-external"),
        default="external-to-joint",
        help="Right-arm coordinate convention passed to the underlying comparison script.",
    )
    parser.add_argument(
        "--use-best-orientation",
        action="store_true",
        help="Use the best right-coordinate/sign variant reported by the underlying comparison script.",
    )
    parser.add_argument(
        "--output",
        default=str(DEFAULT_OUTPUT),
        help="Summary CSV path. Default: results/article_shape_validation_beta15_summary.csv.",
    )
    parser.add_argument(
        "--output-prefix-root",
        default=str(DEFAULT_OUTPUT_PREFIX_ROOT),
        help=(
            "Root path prefix for per-case overlay/sample/diagnostic files. "
            "Default produces results/article_shape_validation_beta15_<branch>_mu... files."
        ),
    )
    args = parser.parse_args(argv)

    if args.beta < 0.0:
        parser.error("--beta must be non-negative.")
    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if args.num_roots <= 0:
        parser.error("--num-roots must be positive.")
    if args.root_window <= 0.0:
        parser.error("--root-window must be positive.")
    if not args.branch_ids:
        parser.error("--branch-ids must contain at least one branch id.")
    if not args.mus:
        parser.error("--mus must contain at least one value.")
    return args


def case_output_prefix(prefix_root: Path, *, beta_deg: float, branch_id: str, mu_value: float, epsilon: float) -> Path:
    return Path(
        f"{prefix_root}_beta{filename_number_token(beta_deg)}_{branch_id}"
        f"_mu{filename_number_token(mu_value)}_eps{filename_number_token(epsilon)}"
    )


def as_float(row: dict[str, float | int | str], key: str) -> float:
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return float("nan")


def as_int(row: dict[str, float | int | str], key: str) -> int:
    try:
        return int(row[key])
    except (KeyError, TypeError, ValueError):
        return -1


def is_finite(value: float) -> bool:
    return bool(np.isfinite(float(value)))


def classify_case(summary_row: dict[str, float | int | str]) -> str:
    rel_lambda_diff = abs(float(summary_row["rel_lambda_diff"]))
    component_all = abs(float(summary_row["component_rel_l2_all"]))
    component_panels = [
        abs(float(summary_row["rel_l2_w_left"])),
        abs(float(summary_row["rel_l2_w_right"])),
        abs(float(summary_row["rel_l2_u_left"])),
        abs(float(summary_row["rel_l2_u_right"])),
    ]
    energy_fraction_diff = abs(
        float(summary_row["analytic_axial_energy_fraction"]) - float(summary_row["fem_axial_energy_fraction"])
    )
    right_axial_share_diff = abs(
        float(summary_row["analytic_right_axial_share"]) - float(summary_row["fem_right_axial_share"])
    )

    required_values = [rel_lambda_diff, component_all, energy_fraction_diff, right_axial_share_diff, *component_panels]
    if not all(is_finite(value) for value in required_values):
        return "needs_review"

    lambda_ok = rel_lambda_diff <= LAMBDA_MATCH_TOL
    components_ok = component_all <= COMPONENT_REL_L2_TOL and all(
        value <= COMPONENT_PANEL_REL_L2_TOL for value in component_panels
    )
    energy_ok = energy_fraction_diff <= ENERGY_FRACTION_TOL and right_axial_share_diff <= RIGHT_AXIAL_SHARE_TOL

    if not lambda_ok:
        return "needs_review"
    if components_ok and energy_ok:
        return "good_match"
    if components_ok and not energy_ok:
        return "energy_mismatch"
    if not components_ok and energy_ok:
        return "component_mismatch"
    return "acceptable_frequency_only"


def build_summary_row(compare_row: dict[str, float | int | str]) -> dict[str, float | int | str]:
    summary_row: dict[str, float | int | str] = {
        "branch_id": str(compare_row["branch_id"]),
        "mu": as_float(compare_row, "mu"),
        "epsilon": as_float(compare_row, "epsilon"),
        "beta": as_float(compare_row, "beta"),
        "fem_lambda": as_float(compare_row, "fem_lambda"),
        "analytic_lambda": as_float(compare_row, "analytic_lambda"),
        "rel_lambda_diff": as_float(compare_row, "relative_lambda_difference"),
        "analytic_root_index": as_int(compare_row, "analytic_root_index"),
        "component_rel_l2_all": as_float(compare_row, "relative_L2_error_all_components"),
        "rel_l2_w_left": as_float(compare_row, "rel_l2_w_left"),
        "rel_l2_w_right": as_float(compare_row, "rel_l2_w_right"),
        "rel_l2_u_left": as_float(compare_row, "rel_l2_u_left"),
        "rel_l2_u_right": as_float(compare_row, "rel_l2_u_right"),
        "analytic_axial_energy_fraction": as_float(compare_row, "analytic_axial_energy_fraction"),
        "fem_axial_energy_fraction": as_float(compare_row, "fem_axial_energy_fraction"),
        "analytic_right_axial_share": as_float(compare_row, "analytic_right_axial_share"),
        "fem_right_axial_share": as_float(compare_row, "fem_right_axial_share"),
        "conclusion_flag": "",
    }
    summary_row["conclusion_flag"] = classify_case(summary_row)
    return summary_row


def write_summary_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def compare_one_case(
    *,
    branch_id: str,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    l_total: float,
    num_roots: int,
    root_window: float,
    right_coordinate: str,
    use_best_orientation: bool,
    output_prefix: Path,
) -> dict[str, float | int | str]:
    compare_args = [
        "--branch-id",
        branch_id,
        "--beta",
        f"{float(beta_deg):g}",
        "--mu",
        f"{float(mu_value):g}",
        "--epsilon",
        f"{float(epsilon):g}",
        "--l-total",
        f"{float(l_total):g}",
        "--num-roots",
        str(int(num_roots)),
        "--root-window",
        f"{float(root_window):g}",
        "--right-coordinate",
        right_coordinate,
        "--output-prefix",
        str(output_prefix),
        "--compare-energies",
    ]
    if use_best_orientation:
        compare_args.append("--use-best-orientation")
    return compare_shape_main(compare_args)


def print_batch_summary(rows: list[dict[str, float | int | str]], output_path: Path) -> None:
    print("Article beta=15 analytic/FEM shape validation")
    print(f"summary_csv: {output_path}")
    print("branch_id, mu, analytic_root_index, rel_lambda_diff, component_rel_l2_all, energy_fraction_diff, flag")
    for row in rows:
        energy_fraction_diff = abs(
            float(row["analytic_axial_energy_fraction"]) - float(row["fem_axial_energy_fraction"])
        )
        print(
            f"{row['branch_id']}, {float(row['mu']):g}, {int(row['analytic_root_index'])}, "
            f"{float(row['rel_lambda_diff']):.3e}, {float(row['component_rel_l2_all']):.3e}, "
            f"{energy_fraction_diff:.3e}, {row['conclusion_flag']}"
        )


def main(argv: Sequence[str] | None = None) -> list[dict[str, float | int | str]]:
    args = parse_args(argv)
    output_path = resolve_repo_path(args.output)
    output_prefix_root = resolve_repo_path(args.output_prefix_root)

    summary_rows: list[dict[str, float | int | str]] = []
    for branch_id in args.branch_ids:
        for mu_value in args.mus:
            prefix = case_output_prefix(
                output_prefix_root,
                beta_deg=float(args.beta),
                branch_id=str(branch_id),
                mu_value=float(mu_value),
                epsilon=float(args.epsilon),
            )
            print(
                f"Running case: branch_id={branch_id}, beta={float(args.beta):g} deg, "
                f"mu={float(mu_value):g}, epsilon={float(args.epsilon):g}"
            )
            compare_row = compare_one_case(
                branch_id=str(branch_id),
                beta_deg=float(args.beta),
                mu_value=float(mu_value),
                epsilon=float(args.epsilon),
                l_total=float(args.l_total),
                num_roots=int(args.num_roots),
                root_window=float(args.root_window),
                right_coordinate=str(args.right_coordinate),
                use_best_orientation=bool(args.use_best_orientation),
                output_prefix=prefix,
            )
            summary_rows.append(build_summary_row(compare_row))

    write_summary_csv(output_path, summary_rows)
    print_batch_summary(summary_rows, output_path)
    return summary_rows


if __name__ == "__main__":
    main()
