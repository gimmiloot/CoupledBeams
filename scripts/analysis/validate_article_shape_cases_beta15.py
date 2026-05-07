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
DEFAULT_REPORT = REPO_ROOT / "results" / "article_shape_validation_beta15_report.md"
DEFAULT_OUTPUT_PREFIX_ROOT = REPO_ROOT / "results" / "article_shape_validation"

SUMMARY_FIELDNAMES = [
    "branch_id",
    "mu",
    "epsilon",
    "beta",
    "analytic_root_index",
    "fem_lambda",
    "analytic_lambda",
    "rel_lambda_diff",
    "component_rel_l2_all",
    "max_abs_error_all_components",
    "rel_l2_u_left",
    "rel_l2_w_left",
    "rel_l2_u_right",
    "rel_l2_w_right",
    "max_abs_error_u_left",
    "max_abs_error_w_left",
    "max_abs_error_u_right",
    "max_abs_error_w_right",
    "transverse_rel_l2_all",
    "transverse_shape_dot",
    "max_abs_error_transverse",
    "axial_rel_l2_all",
    "axial_shape_dot",
    "analytic_axial_energy_fraction",
    "fem_direct_axial_energy_fraction",
    "analytic_right_axial_share",
    "fem_direct_right_axial_share",
    "analytic_right_bending_share",
    "fem_direct_right_bending_share",
    "energy_source",
    "frequency_flag",
    "transverse_shape_flag",
    "full_component_flag",
    "energy_flag",
    "article_use_flag",
]

FREQUENCY_GOOD_TOL = 1e-5
SHAPE_GOOD_TOL = 0.05
SHAPE_REVIEW_TOL = 0.15
ENERGY_GOOD_TOL = 0.05
ENERGY_REVIEW_TOL = 0.15
NEAR_ZERO_NORM = 1e-12


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
        "--report",
        default=str(DEFAULT_REPORT),
        help="Compact markdown report path. Default: results/article_shape_validation_beta15_report.md.",
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


def safe_relative_l2(analytic: np.ndarray, fem_values: np.ndarray) -> float:
    denominator = float(np.linalg.norm(fem_values))
    if denominator <= NEAR_ZERO_NORM:
        return float("nan")
    return float(np.linalg.norm(analytic - fem_values) / denominator)


def safe_shape_dot(analytic: np.ndarray, fem_values: np.ndarray) -> float:
    denominator = float(np.linalg.norm(analytic) * np.linalg.norm(fem_values))
    if denominator <= NEAR_ZERO_NORM:
        return float("nan")
    return float(np.dot(analytic, fem_values) / denominator)


def channel_metrics(analytic: np.ndarray, fem_values: np.ndarray) -> dict[str, float]:
    diff = analytic - fem_values
    return {
        "relative_l2": safe_relative_l2(analytic, fem_values),
        "shape_dot": safe_shape_dot(analytic, fem_values),
        "max_abs_error": float(np.max(np.abs(diff))) if diff.size else float("nan"),
    }


def load_sample_metrics(path: Path) -> dict[str, float]:
    values: dict[str, list[float]] = {
        "fem_u_left": [],
        "analytic_u_left": [],
        "fem_w_left": [],
        "analytic_w_left": [],
        "fem_u_right": [],
        "analytic_u_right": [],
        "fem_w_right": [],
        "analytic_w_right": [],
    }
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            for key in values:
                values[key].append(float(row[key]))

    fem_transverse = np.concatenate(
        [
            np.asarray(values["fem_w_left"], dtype=float),
            np.asarray(values["fem_w_right"], dtype=float),
        ]
    )
    analytic_transverse = np.concatenate(
        [
            np.asarray(values["analytic_w_left"], dtype=float),
            np.asarray(values["analytic_w_right"], dtype=float),
        ]
    )
    fem_axial = np.concatenate(
        [
            np.asarray(values["fem_u_left"], dtype=float),
            np.asarray(values["fem_u_right"], dtype=float),
        ]
    )
    analytic_axial = np.concatenate(
        [
            np.asarray(values["analytic_u_left"], dtype=float),
            np.asarray(values["analytic_u_right"], dtype=float),
        ]
    )

    transverse = channel_metrics(analytic_transverse, fem_transverse)
    axial = channel_metrics(analytic_axial, fem_axial)
    return {
        "transverse_rel_l2_all": transverse["relative_l2"],
        "transverse_shape_dot": transverse["shape_dot"],
        "max_abs_error_transverse": transverse["max_abs_error"],
        "axial_rel_l2_all": axial["relative_l2"],
        "axial_shape_dot": axial["shape_dot"],
    }


def frequency_flag(rel_lambda_diff: float) -> str:
    return "good" if is_finite(rel_lambda_diff) and abs(float(rel_lambda_diff)) < FREQUENCY_GOOD_TOL else "bad"


def threshold_flag(value: float, *, good_tol: float, review_tol: float) -> str:
    if not is_finite(value):
        return "bad"
    magnitude = abs(float(value))
    if magnitude < good_tol:
        return "good"
    if magnitude < review_tol:
        return "review"
    return "bad"


def energy_flag(analytic_fraction: float, fem_direct_fraction: float) -> str:
    if not is_finite(analytic_fraction) or not is_finite(fem_direct_fraction):
        return "bad"
    return threshold_flag(
        abs(float(analytic_fraction) - float(fem_direct_fraction)),
        good_tol=ENERGY_GOOD_TOL,
        review_tol=ENERGY_REVIEW_TOL,
    )


def article_use_flag(*, freq_flag: str, transverse_flag: str) -> str:
    if freq_flag != "good":
        return "needs_review"
    if transverse_flag == "good":
        return "safe_transverse_plot"
    if transverse_flag == "review":
        return "frequency_only"
    return "needs_review"


def build_summary_row(compare_row: dict[str, float | int | str]) -> dict[str, float | int | str]:
    sample_path = resolve_repo_path(str(compare_row["components_samples_csv"]))
    sample_metrics = load_sample_metrics(sample_path)

    rel_lambda_diff = as_float(compare_row, "relative_lambda_difference")
    component_rel_l2_all = as_float(compare_row, "relative_L2_error_all_components")
    analytic_axial_fraction = as_float(compare_row, "analytic_axial_energy_fraction")
    fem_direct_axial_fraction = as_float(compare_row, "fem_direct_axial_energy_fraction")

    freq = frequency_flag(rel_lambda_diff)
    transverse = threshold_flag(
        sample_metrics["transverse_rel_l2_all"],
        good_tol=SHAPE_GOOD_TOL,
        review_tol=SHAPE_REVIEW_TOL,
    )
    full_component = threshold_flag(
        component_rel_l2_all,
        good_tol=SHAPE_GOOD_TOL,
        review_tol=SHAPE_REVIEW_TOL,
    )
    energy = energy_flag(analytic_axial_fraction, fem_direct_axial_fraction)

    return {
        "branch_id": str(compare_row["branch_id"]),
        "mu": as_float(compare_row, "mu"),
        "epsilon": as_float(compare_row, "epsilon"),
        "beta": as_float(compare_row, "beta"),
        "analytic_root_index": as_int(compare_row, "analytic_root_index"),
        "fem_lambda": as_float(compare_row, "fem_lambda"),
        "analytic_lambda": as_float(compare_row, "analytic_lambda"),
        "rel_lambda_diff": rel_lambda_diff,
        "component_rel_l2_all": component_rel_l2_all,
        "max_abs_error_all_components": as_float(compare_row, "max_abs_error_all_components"),
        "rel_l2_u_left": as_float(compare_row, "rel_l2_u_left"),
        "rel_l2_w_left": as_float(compare_row, "rel_l2_w_left"),
        "rel_l2_u_right": as_float(compare_row, "rel_l2_u_right"),
        "rel_l2_w_right": as_float(compare_row, "rel_l2_w_right"),
        "max_abs_error_u_left": as_float(compare_row, "max_abs_error_u_left"),
        "max_abs_error_w_left": as_float(compare_row, "max_abs_error_w_left"),
        "max_abs_error_u_right": as_float(compare_row, "max_abs_error_u_right"),
        "max_abs_error_w_right": as_float(compare_row, "max_abs_error_w_right"),
        "transverse_rel_l2_all": sample_metrics["transverse_rel_l2_all"],
        "transverse_shape_dot": sample_metrics["transverse_shape_dot"],
        "max_abs_error_transverse": sample_metrics["max_abs_error_transverse"],
        "axial_rel_l2_all": sample_metrics["axial_rel_l2_all"],
        "axial_shape_dot": sample_metrics["axial_shape_dot"],
        "analytic_axial_energy_fraction": analytic_axial_fraction,
        "fem_direct_axial_energy_fraction": fem_direct_axial_fraction,
        "analytic_right_axial_share": as_float(compare_row, "analytic_right_axial_share"),
        "fem_direct_right_axial_share": as_float(compare_row, "fem_direct_right_axial_share"),
        "analytic_right_bending_share": as_float(compare_row, "analytic_right_bending_share"),
        "fem_direct_right_bending_share": as_float(compare_row, "fem_direct_right_bending_share"),
        "energy_source": "direct_element_stiffness",
        "frequency_flag": freq,
        "transverse_shape_flag": transverse,
        "full_component_flag": full_component,
        "energy_flag": energy,
        "article_use_flag": article_use_flag(freq_flag=freq, transverse_flag=transverse),
    }


def write_summary_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def markdown_float(value: float, precision: int = 3) -> str:
    if not is_finite(value):
        return "nan"
    return f"{float(value):.{precision}e}"


def write_markdown_report(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Article shape validation beta=15",
        "",
        "Energy source: `direct_element_stiffness`.",
        "",
        "| branch_id | mu | rel Lambda diff | transverse rel L2 | full rel L2 | analytic axial frac | FEM direct axial frac | article use |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            f"{row['branch_id']} | "
            f"{float(row['mu']):g} | "
            f"{markdown_float(float(row['rel_lambda_diff']))} | "
            f"{markdown_float(float(row['transverse_rel_l2_all']))} | "
            f"{markdown_float(float(row['component_rel_l2_all']))} | "
            f"{markdown_float(float(row['analytic_axial_energy_fraction']))} | "
            f"{markdown_float(float(row['fem_direct_axial_energy_fraction']))} | "
            f"{row['article_use_flag']} |"
        )
    lines.extend(
        [
            "",
            "Flags:",
            "",
            "- `frequency_flag`: good when `rel_lambda_diff < 1e-5`.",
            "- `transverse_shape_flag`: good below 0.05, review below 0.15, bad otherwise.",
            "- `full_component_flag`: good below 0.05, review below 0.15, bad otherwise.",
            "- `energy_flag`: good when axial-fraction difference is below 0.05, review below 0.15, bad otherwise.",
            "- `article_use_flag` evaluates transverse plotting separately from axial-component mismatch.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


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
        "--check-fem-direct-energy",
    ]
    if use_best_orientation:
        compare_args.append("--use-best-orientation")
    return compare_shape_main(compare_args)


def print_batch_summary(
    rows: list[dict[str, float | int | str]],
    *,
    output_path: Path,
    report_path: Path,
) -> None:
    print("Article beta=15 analytic/FEM shape validation")
    print(f"summary_csv: {output_path}")
    print(f"markdown_report: {report_path}")
    print(
        "branch_id, mu, rel_lambda_diff, transverse_rel_l2_all, component_rel_l2_all, "
        "analytic_axial_fraction, fem_direct_axial_fraction, article_use_flag"
    )
    for row in rows:
        print(
            f"{row['branch_id']}, {float(row['mu']):g}, "
            f"{float(row['rel_lambda_diff']):.3e}, {float(row['transverse_rel_l2_all']):.3e}, "
            f"{float(row['component_rel_l2_all']):.3e}, "
            f"{float(row['analytic_axial_energy_fraction']):.3e}, "
            f"{float(row['fem_direct_axial_energy_fraction']):.3e}, {row['article_use_flag']}"
        )


def main(argv: Sequence[str] | None = None) -> list[dict[str, float | int | str]]:
    args = parse_args(argv)
    output_path = resolve_repo_path(args.output)
    report_path = resolve_repo_path(args.report)
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
    write_markdown_report(report_path, summary_rows)
    print_batch_summary(summary_rows, output_path=output_path, report_path=report_path)
    return summary_rows


if __name__ == "__main__":
    main()
