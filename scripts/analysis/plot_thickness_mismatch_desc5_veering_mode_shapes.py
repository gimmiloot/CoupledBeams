from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis import plot_thickness_mismatch_branch_shapes_vs_eta as shape_plotter  # noqa: E402


BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5
DESCENDANT_ID = 5

MU_STEP = 0.001
NUM_SORTED_ROOTS = 8
NUM_TRACKED_BRANCHES = 8
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401
NORMALIZE = "max-full"

SOURCE_DESC6_DIR = REPO_ROOT / "results" / "thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5"
SOURCE_DESC6_POINTS_CSV = SOURCE_DESC6_DIR / "desc6_veering_points.csv"

OUTPUT_DIR = REPO_ROOT / "results" / "thickness_mismatch_desc5_veering_mode_shapes_beta15_eps0p0025_eta0p5"
OUTPUT_CSV = OUTPUT_DIR / "desc5_veering_points.csv"
OUTPUT_REPORT = OUTPUT_DIR / "desc5_veering_points_report.md"

CSV_FIELDNAMES = [
    "veering_case",
    "rank_within_case",
    "mu",
    "eta",
    "beta_deg",
    "epsilon",
    "descendant_id",
    "lambda_desc5",
    "sorted_position",
    "nearest_sorted_pair",
    "local_gap",
    "source_mu_from_desc6_table",
    "tracking_mac",
    "tracking_warning",
    "plot_file",
]

EXPECTED_POINTS = {
    ("near_4_5", 1): 0.379,
    ("near_4_5", 2): 0.037,
    ("near_4_5", 3): 0.660,
    ("near_5_6", 1): 0.716,
    ("near_5_6", 2): 0.488,
    ("near_5_6", 3): 0.217,
}


@dataclass(frozen=True)
class VeeringPoint:
    veering_case: str
    rank_within_case: int
    mu: float
    nearest_sorted_pair: str
    local_gap: float
    source_file: Path

    @property
    def plot_stem(self) -> str:
        if self.veering_case == "near_4_5":
            return f"desc5_near45_rank{self.rank_within_case}"
        if self.veering_case == "near_5_6":
            return f"desc5_near56_rank{self.rank_within_case}"
        raise ValueError(f"Unexpected veering case: {self.veering_case}")

    @property
    def title_pair(self) -> str:
        if self.nearest_sorted_pair == "4-5":
            return "4-5"
        if self.nearest_sorted_pair == "5-6":
            return "5-6"
        return self.nearest_sorted_pair


def display_path(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_float(row: dict[str, str], key: str) -> float:
    value = str(row.get(key, "")).strip()
    if value == "":
        return float("nan")
    return float(value)


def read_source_points() -> list[VeeringPoint]:
    if not SOURCE_DESC6_POINTS_CSV.exists():
        raise FileNotFoundError(f"Missing required source CSV: {SOURCE_DESC6_POINTS_CSV}")

    points: list[VeeringPoint] = []
    with SOURCE_DESC6_POINTS_CSV.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            veering_case = str(row["veering_case"])
            rank = int(row["rank_within_case"])
            key = (veering_case, rank)
            if key not in EXPECTED_POINTS:
                continue
            mu = parse_float(row, "mu")
            expected_mu = EXPECTED_POINTS[key]
            if not np.isclose(mu, expected_mu, rtol=0.0, atol=5e-13):
                raise RuntimeError(
                    f"Source point {veering_case} rank {rank} has mu={mu:g}, "
                    f"expected {expected_mu:g}."
                )
            points.append(
                VeeringPoint(
                    veering_case=veering_case,
                    rank_within_case=rank,
                    mu=mu,
                    nearest_sorted_pair=str(row["nearest_sorted_pair"]),
                    local_gap=parse_float(row, "local_gap"),
                    source_file=SOURCE_DESC6_POINTS_CSV,
                )
            )

    points.sort(key=lambda point: (point.veering_case, point.rank_within_case))
    if len(points) != len(EXPECTED_POINTS):
        raise RuntimeError(f"Expected {len(EXPECTED_POINTS)} source points, found {len(points)}.")
    return points


def build_title(point: VeeringPoint) -> str:
    return (
        f"Thickness mismatch descendant 5 near {point.title_pair}: "
        f"beta = {BETA_DEG:g} deg, epsilon = {EPSILON:g}, "
        f"eta = {ETA:g}, mu = {point.mu:g}"
    )


def shape_args(points: list[VeeringPoint]):
    args = shape_plotter.parse_args(
        [
            "--branch-index",
            str(DESCENDANT_ID),
            "--beta-deg",
            str(BETA_DEG),
            "--epsilon",
            str(EPSILON),
            "--mus",
            *[f"{point.mu:.12g}" for point in points],
            "--etas",
            str(ETA),
            "--mu-step",
            str(MU_STEP),
            "--num-tracked-branches",
            str(NUM_TRACKED_BRANCHES),
            "--num-sorted-roots",
            str(NUM_SORTED_ROOTS),
            "--root-scan-step",
            str(ROOT_SCAN_STEP),
            "--root-lmax0",
            str(ROOT_LMAX0),
            "--num-shape-samples",
            str(NUM_SHAPE_SAMPLES),
            "--normalize",
            NORMALIZE,
            "--output-dir",
            str(OUTPUT_DIR),
            "--plot-stems",
            *[point.plot_stem for point in points],
            "--plot-titles",
            *[build_title(point) for point in points],
        ]
    )
    shape_plotter.validate_args(args)
    args.mus = [float(point.mu) for point in points]
    args.etas = [ETA]
    args.output_dir = OUTPUT_DIR
    return args


def tracking_mac_text(row: dict[str, float | int | str]) -> str:
    for key in ("accepted_mac_to_previous", "diagnostic_candidate_mac_to_previous", "mac_to_previous"):
        value = row.get(key, "")
        try:
            value_f = float(value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value_f):
            return f"{value_f:.10g}"
    return ""


def tracking_warning_text(row: dict[str, float | int | str], *, sorted_position: int) -> str:
    flags = []
    if sorted_position != 5:
        flags.append("sorted_position_not_5")
    for key, label in (
        ("requires_refined_check", "requires_refined_check"),
        ("suspicious_assignment", "suspicious_assignment"),
        ("unresolved_assignment", "unresolved_assignment"),
        ("frequency_mac_disagreement", "frequency_mac_disagreement"),
        ("low_mac", "low_mac"),
        ("low_margin", "low_margin"),
    ):
        if str(row.get(key, "no")) == "yes":
            flags.append(label)
    return ";".join(dict.fromkeys(flags)) if flags else "no"


def fmt(value: float) -> str:
    if not np.isfinite(float(value)):
        return "nan"
    return f"{float(value):.12g}"


def build_shape_rows(points: list[VeeringPoint]) -> list[dict[str, str]]:
    args = shape_args(points)
    tracking_mu_values = shape_plotter.tracking_grid_for_targets([point.mu for point in points], MU_STEP)
    tracking_cases, _tracking_rows_by_eta = shape_plotter.compute_tracking_cases(args, tracking_mu_values)
    shape_cases_by_mu = shape_plotter.build_cases_by_mu(args, tracking_cases)
    s_norm = np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float)

    rows: list[dict[str, str]] = []
    for plot_index, point in enumerate(points):
        cases = shape_cases_by_mu[float(point.mu)]
        if len(cases) != 1:
            raise RuntimeError(f"Expected one eta case for mu={point.mu:g}, found {len(cases)}.")
        output_png = shape_plotter.plot_cases_for_mu(
            args,
            float(point.mu),
            cases,
            s_norm,
            plot_index=plot_index,
            title_override=build_title(point),
        )
        shape_case = cases[0]
        tracking_case = tracking_cases[ETA][float(point.mu)]
        tracking_row = tracking_case.tracking_row
        sorted_position = int(shape_case.current_sorted_index)
        rows.append(
            {
                "veering_case": point.veering_case,
                "rank_within_case": str(point.rank_within_case),
                "mu": fmt(point.mu),
                "eta": fmt(ETA),
                "beta_deg": fmt(BETA_DEG),
                "epsilon": fmt(EPSILON),
                "descendant_id": str(DESCENDANT_ID),
                "lambda_desc5": fmt(shape_case.Lambda),
                "sorted_position": str(sorted_position),
                "nearest_sorted_pair": point.nearest_sorted_pair,
                "local_gap": fmt(point.local_gap),
                "source_mu_from_desc6_table": display_path(point.source_file),
                "tracking_mac": tracking_mac_text(tracking_row),
                "tracking_warning": tracking_warning_text(tracking_row, sorted_position=sorted_position),
                "plot_file": display_path(output_png),
            }
        )
        print(
            f"saved descendant-5 shape PNG for {point.veering_case} rank {point.rank_within_case}: "
            f"{output_png}"
        )
    return rows


def write_csv(rows: list[dict[str, str]]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| veering case | rank | mu | Lambda desc5 | sorted position | local gap | plot | tracking warning |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["veering_case"],
                    row["rank_within_case"],
                    row["mu"],
                    row["lambda_desc5"],
                    row["sorted_position"],
                    row["local_gap"],
                    f"`{row['plot_file']}`",
                    row["tracking_warning"],
                ]
            )
            + " |"
        )
    return lines


def write_report(rows: list[dict[str, str]]) -> None:
    warning_rows = [row for row in rows if row["tracking_warning"] != "no"]
    sorted_not_5 = [row for row in rows if row["sorted_position"] != "5"]
    all_sorted5 = not sorted_not_5

    lines = [
        "# Descendant 5 Veering Mode-Shape Points",
        "",
        "## Purpose",
        "",
        "This diagnostic corrects the previous descendant/sorted identity confusion",
        "for the eta=0.5 near-veering shape plots. The recent `Lambda(beta)` audit",
        "showed that sorted 5 is descendant 5 at the target eta=0.5 point, so this",
        "run builds descendant 5 shapes at the six already selected near-approach",
        "points.",
        "",
        "## Data Source",
        "",
        f"- source table: `{display_path(SOURCE_DESC6_POINTS_CSV)}`",
        "- the source table is used only for the six selected `mu` values, near-approach labels, ranks, and local gaps",
        "- mode shapes are reconstructed for descendant 5, not descendant 6",
        "- no new veering-point selection is performed",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- descendant id: {DESCENDANT_ID}",
        f"- tracking grid step in mu: {MU_STEP:g}",
        f"- sorted roots used for tracking: first {NUM_SORTED_ROOTS}",
        f"- shape normalization: `{NORMALIZE}`",
        "",
        "## Six Plots",
        "",
    ]
    lines.extend(markdown_table(rows))

    lines.extend(
        [
            "",
            "## Explicit Answer",
            "",
            f"- descendant 5 equals sorted 5 in all six selected points: {'yes' if all_sorted5 else 'no'}",
        ]
    )
    if sorted_not_5:
        lines.append("- sorted position differs from 5 at:")
        for row in sorted_not_5:
            lines.append(
                f"  - {row['veering_case']} rank {row['rank_within_case']}, "
                f"mu={row['mu']}: sorted position {row['sorted_position']}"
            )
    else:
        lines.append("- no selected point has descendant 5 away from sorted position 5.")

    lines.extend(
        [
            "",
            "## Relation To Previous Descendant-6 Plots",
            "",
            "The previous descendant-6 plots are retained as useful descendant-6",
            "diagnostics. They are not the intended sorted-5 shapes for eta=0.5.",
            "This new output directory is the corrected descendant-5 set for the",
            "original sorted-5 identity question.",
            "",
            "## Tracking Warnings And Ambiguities",
            "",
            f"- rows with warnings: {len(warning_rows)}",
        ]
    )
    if warning_rows:
        for row in warning_rows:
            lines.append(
                f"  - {row['veering_case']} rank {row['rank_within_case']}, "
                f"mu={row['mu']}: {row['tracking_warning']}"
            )
    else:
        lines.append("- no tracking warnings were found for the six descendant-5 plotted rows.")

    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- Euler-Bernoulli analytic-only thickness-mismatch diagnostic",
            "- no Timoshenko curves",
            "- no FEM, Gmsh, CalculiX, or heavy 3D FEM",
            "- no new veering-point selection",
            "- shape sign and absolute scale are arbitrary",
            "",
            "This diagnostic writes only the descendant-5 output directory. It does not",
            "delete or overwrite descendant-6 outputs, article files, article figures,",
            "`paper_dorofeev_style`, the old determinant, `src/my_project/analytic/formulas.py`,",
            "old solvers, baseline results, or the FEM physical model.",
            "",
            "## Outputs",
            "",
            f"- CSV: `{display_path(OUTPUT_CSV)}`",
            f"- report: `{display_path(OUTPUT_REPORT)}`",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, object]:
    print("Descendant-5 corrected veering mode-shape diagnostic")
    print(f"reading selected points from: {SOURCE_DESC6_POINTS_CSV}")
    points = read_source_points()
    rows = build_shape_rows(points)
    write_csv(rows)
    write_report(rows)
    print(f"saved descendant-5 veering CSV: {OUTPUT_CSV}")
    print(f"saved descendant-5 veering report: {OUTPUT_REPORT}")
    return {"csv": OUTPUT_CSV, "report": OUTPUT_REPORT, "rows": rows}


if __name__ == "__main__":
    main()
