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

from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
from scripts.analysis import plot_thickness_mismatch_branch_shapes_vs_eta as shape_plotter  # noqa: E402


BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5
DESCENDANT_ID = 6
MU_STEP = 0.001
NUM_SORTED_ROOTS = 8
NUM_TRACKED_BRANCHES = 8
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_TRACKING_SHAPE_SAMPLES = 401
NUM_INTEGRATION_SAMPLES = 2001
NORMALIZE = "max-full"
LOCALIZATION_THRESHOLD = 0.75

OUTPUT_DIR = REPO_ROOT / "results" / "thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5"
INPUT_POINTS_CSV = OUTPUT_DIR / "desc6_veering_points.csv"
OUTPUT_CSV = OUTPUT_DIR / "desc6_localization_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "desc6_localization_audit_report.md"

CSV_FIELDNAMES = [
    "veering_case",
    "rank_within_case",
    "mu",
    "lambda_desc6",
    "sorted_position",
    "local_gap",
    "tau1",
    "tau2",
    "l1",
    "l2",
    "U_w_rod1",
    "U_w_rod2",
    "rod1_displacement_fraction",
    "rod2_displacement_fraction",
    "A_rod1",
    "A_rod2",
    "max_abs_w_ratio_rod1_over_rod2",
    "U_slope_rod1",
    "U_slope_rod2",
    "rod1_slope_fraction",
    "rod2_slope_fraction",
    "U_curv_rod1",
    "U_curv_rod2",
    "rod1_curvature_fraction",
    "rod2_curvature_fraction",
    "Eb_rod1",
    "Eb_rod2",
    "rod1_bending_energy_fraction",
    "rod2_bending_energy_fraction",
    "displacement_classification",
    "energy_classification",
    "tracking_warning",
]


@dataclass(frozen=True)
class VeeringPoint:
    veering_case: str
    rank_within_case: int
    mu: float
    local_gap: float
    tracking_warning: str


def display_path(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value == "":
        return float("nan")
    return float(value)


def read_points() -> list[VeeringPoint]:
    if not INPUT_POINTS_CSV.exists():
        raise FileNotFoundError(f"Missing required input CSV: {INPUT_POINTS_CSV}")
    points: list[VeeringPoint] = []
    with INPUT_POINTS_CSV.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            points.append(
                VeeringPoint(
                    veering_case=str(row["veering_case"]),
                    rank_within_case=int(row["rank_within_case"]),
                    mu=parse_float(row, "mu"),
                    local_gap=parse_float(row, "local_gap"),
                    tracking_warning=str(row.get("tracking_warning", "")),
                )
            )
    if len(points) != 6:
        raise RuntimeError(f"Expected 6 descendant-6 veering points, found {len(points)}.")
    return points


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
            str(NUM_TRACKING_SHAPE_SAMPLES),
            "--normalize",
            NORMALIZE,
        ]
    )
    shape_plotter.validate_args(args)
    args.mus = [float(point.mu) for point in points]
    args.etas = [ETA]
    return args


def fraction(left: float, right: float) -> tuple[float, float]:
    total = float(left) + float(right)
    if total <= 0.0:
        return float("nan"), float("nan")
    return float(left) / total, float(right) / total


def ratio(left: float, right: float) -> float:
    if abs(float(right)) <= 1e-30:
        return float("inf")
    return float(left) / float(right)


def classify_displacement(rod1_fraction: float, rod2_fraction: float) -> str:
    if rod1_fraction >= LOCALIZATION_THRESHOLD:
        return "displacement_localized_rod1"
    if rod2_fraction >= LOCALIZATION_THRESHOLD:
        return "displacement_localized_rod2"
    return "displacement_mixed"


def classify_energy(rod1_fraction: float, rod2_fraction: float) -> str:
    if rod1_fraction >= LOCALIZATION_THRESHOLD:
        return "energy_localized_rod1"
    if rod2_fraction >= LOCALIZATION_THRESHOLD:
        return "energy_localized_rod2"
    return "energy_mixed"


def integrate(values: np.ndarray, coordinates: np.ndarray) -> float:
    integrate_func = getattr(np, "trapezoid", None)
    if integrate_func is not None:
        return float(integrate_func(values, coordinates))
    return float(np.trapz(values, coordinates))


def transverse_metrics(w_values: np.ndarray, length: float) -> dict[str, float]:
    xi = np.linspace(0.0, 1.0, len(w_values), dtype=float)
    s_values = float(length) * xi
    w = np.asarray(w_values, dtype=float)
    dw_ds = np.gradient(w, s_values, edge_order=2)
    d2w_ds2 = np.gradient(dw_ds, s_values, edge_order=2)
    return {
        "U_w": max(integrate(w * w, s_values), 0.0),
        "A": float(np.max(np.abs(w))),
        "U_slope": max(integrate(dw_ds * dw_ds, s_values), 0.0),
        "U_curv": max(integrate(d2w_ds2 * d2w_ds2, s_values), 0.0),
    }


def fmt(value: float) -> str:
    if not np.isfinite(float(value)):
        return str(float(value))
    return f"{float(value):.12g}"


def compute_rows(points: list[VeeringPoint]) -> list[dict[str, str]]:
    args = shape_args(points)
    target_mus = [float(point.mu) for point in points]
    tracking_mu_values = shape_plotter.tracking_grid_for_targets(target_mus, MU_STEP)
    tracking_cases, _rows_by_eta = shape_plotter.compute_tracking_cases(args, tracking_mu_values)

    s_norm = np.linspace(0.0, 1.0, NUM_INTEGRATION_SAMPLES, dtype=float)
    out_rows: list[dict[str, str]] = []
    for point in points:
        tracking_case = tracking_cases[ETA][float(point.mu)]
        shape_case = shape_plotter.build_shape_case(args, tracking_case, s_norm)
        factors = thickness_mismatch_factors(float(point.mu), ETA)
        l1 = 1.0 - float(point.mu)
        l2 = 1.0 + float(point.mu)

        rod1 = transverse_metrics(shape_case.components["w_left"], l1)
        rod2 = transverse_metrics(shape_case.components["w_right"], l2)
        disp1, disp2 = fraction(rod1["U_w"], rod2["U_w"])
        slope1, slope2 = fraction(rod1["U_slope"], rod2["U_slope"])
        curv1, curv2 = fraction(rod1["U_curv"], rod2["U_curv"])

        eb1 = factors.tau1**4 * rod1["U_curv"]
        eb2 = factors.tau2**4 * rod2["U_curv"]
        energy1, energy2 = fraction(eb1, eb2)

        out_rows.append(
            {
                "veering_case": point.veering_case,
                "rank_within_case": str(point.rank_within_case),
                "mu": fmt(point.mu),
                "lambda_desc6": fmt(shape_case.Lambda),
                "sorted_position": str(shape_case.current_sorted_index),
                "local_gap": fmt(point.local_gap),
                "tau1": fmt(factors.tau1),
                "tau2": fmt(factors.tau2),
                "l1": fmt(l1),
                "l2": fmt(l2),
                "U_w_rod1": fmt(rod1["U_w"]),
                "U_w_rod2": fmt(rod2["U_w"]),
                "rod1_displacement_fraction": fmt(disp1),
                "rod2_displacement_fraction": fmt(disp2),
                "A_rod1": fmt(rod1["A"]),
                "A_rod2": fmt(rod2["A"]),
                "max_abs_w_ratio_rod1_over_rod2": fmt(ratio(rod1["A"], rod2["A"])),
                "U_slope_rod1": fmt(rod1["U_slope"]),
                "U_slope_rod2": fmt(rod2["U_slope"]),
                "rod1_slope_fraction": fmt(slope1),
                "rod2_slope_fraction": fmt(slope2),
                "U_curv_rod1": fmt(rod1["U_curv"]),
                "U_curv_rod2": fmt(rod2["U_curv"]),
                "rod1_curvature_fraction": fmt(curv1),
                "rod2_curvature_fraction": fmt(curv2),
                "Eb_rod1": fmt(eb1),
                "Eb_rod2": fmt(eb2),
                "rod1_bending_energy_fraction": fmt(energy1),
                "rod2_bending_energy_fraction": fmt(energy2),
                "displacement_classification": classify_displacement(disp1, disp2),
                "energy_classification": classify_energy(energy1, energy2),
                "tracking_warning": point.tracking_warning,
            }
        )
    return out_rows


def write_csv(rows: list[dict[str, str]]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def class_side(label: str) -> str:
    if label.endswith("rod1"):
        return "rod1"
    if label.endswith("rod2"):
        return "rod2"
    return "mixed"


def markdown_table(rows: list[dict[str, str]]) -> list[str]:
    lines = [
        "| case | rank | mu | disp rod1 | disp rod2 | energy rod1 | energy rod2 | disp class | energy class |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["veering_case"],
                    row["rank_within_case"],
                    row["mu"],
                    row["rod1_displacement_fraction"],
                    row["rod2_displacement_fraction"],
                    row["rod1_bending_energy_fraction"],
                    row["rod2_bending_energy_fraction"],
                    row["displacement_classification"],
                    row["energy_classification"],
                ]
            )
            + " |"
        )
    return lines


def write_report(rows: list[dict[str, str]]) -> None:
    displacement_classes = sorted({row["displacement_classification"] for row in rows})
    energy_classes = sorted({row["energy_classification"] for row in rows})
    warning_rows = [row for row in rows if row["tracking_warning"].lower() == "yes"]
    disagreements = [
        row
        for row in rows
        if class_side(row["displacement_classification"]) != class_side(row["energy_classification"])
    ]

    lines = [
        "# Descendant 6 Localization Audit",
        "",
        "## What Was Checked",
        "",
        "This audit uses the six already selected descendant-6 veering points from",
        f"`{display_path(INPUT_POINTS_CSV)}`. No new mu values are selected and the",
        "veering-point selection is not recomputed.",
        "",
        "For each point, the Euler-Bernoulli thickness-mismatch descendant-6 shape is",
        "reconstructed with the same shape-MAC branch identity workflow used by the",
        "mode-shape plotter. Sorted position is diagnostic metadata only.",
        "",
        "## Parameters And Normalization",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- descendant id: {DESCENDANT_ID}",
        f"- integration samples per rod: {NUM_INTEGRATION_SAMPLES}",
        f"- shape normalization: `{NORMALIZE}` from the existing shape plotter",
        "- integrals are computed over each rod physical local coordinate",
        "- bending-energy proxy uses relative circular-section stiffness `J_i/J0 = tau_i^4`",
        "",
        "The absolute eigenvector scale is arbitrary, but the reported fractions and",
        "rod1/rod2 ratios are invariant under a common scaling of the full mode shape.",
        "",
        "## Six Points",
        "",
    ]
    for row in rows:
        lines.append(
            f"- {row['veering_case']} rank {row['rank_within_case']}: "
            f"mu={row['mu']}, Lambda={row['lambda_desc6']}, sorted={row['sorted_position']}, "
            f"local_gap={row['local_gap']}"
        )

    lines.extend(["", "## Fraction Table", ""])
    lines.extend(markdown_table(rows))

    lines.extend(
        [
            "",
            "## Conclusions",
            "",
            f"- displacement classifications present: {', '.join(displacement_classes)}",
            f"- bending-energy classifications present: {', '.join(energy_classes)}",
        ]
    )
    if all(row["displacement_classification"] == "displacement_mixed" for row in rows):
        lines.append("- by displacement L2 norm, none of the six points is localized at the 0.75 threshold.")
    else:
        lines.append("- by displacement L2 norm, at least one point crosses the 0.75 localization threshold.")
    if all(row["energy_classification"] == "energy_mixed" for row in rows):
        lines.append("- by tau-weighted EB bending energy, none of the six points is localized at the 0.75 threshold.")
    else:
        lines.append("- by tau-weighted EB bending energy, at least one point crosses the 0.75 localization threshold.")

    lines.extend(["", "## Visual-Amplitude Versus Bending-Energy Check", ""])
    if disagreements:
        lines.append("The displacement and bending-energy classifications differ in these rows:")
        for row in disagreements:
            lines.append(
                f"- {row['veering_case']} rank {row['rank_within_case']}, mu={row['mu']}: "
                f"{row['displacement_classification']} vs {row['energy_classification']}"
            )
    else:
        lines.append("No row has different localization sides by displacement fraction and bending-energy fraction.")
    lines.append(
        "A larger smooth displacement amplitude need not imply larger bending energy, because "
        "bending energy depends on curvature and on the tau-weighted section stiffness."
    )

    lines.extend(["", "## Transition From near_4_5 To near_5_6", ""])
    for case in ("near_4_5", "near_5_6"):
        case_rows = [row for row in rows if row["veering_case"] == case]
        disp_values = ", ".join(
            f"rank {row['rank_within_case']}: rod1 {row['rod1_displacement_fraction']}"
            for row in case_rows
        )
        energy_values = ", ".join(
            f"rank {row['rank_within_case']}: rod1 {row['rod1_bending_energy_fraction']}"
            for row in case_rows
        )
        lines.append(f"- {case} displacement rod1 fractions: {disp_values}")
        lines.append(f"- {case} bending-energy rod1 fractions: {energy_values}")

    lines.extend(["", "## Tracking Warnings", ""])
    if warning_rows:
        lines.append(f"- WARNING: {len(warning_rows)} row(s) carry tracking warnings:")
        for row in warning_rows:
            lines.append(f"  - {row['veering_case']} rank {row['rank_within_case']}, mu={row['mu']}")
    else:
        lines.append("- no tracking warnings are present in the input six-point table.")

    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- this is an analytic-only Euler-Bernoulli thickness-mismatch shape audit",
            "- no Timoshenko curves, FEM markers, Gmsh, CalculiX, or 3D FEM tasks are used",
            "- shape sign and absolute scale are arbitrary; fractions are common-scale invariant",
            "- if future tracking warnings appear, localization should be interpreted cautiously",
            "",
            "This diagnostic writes only localization audit outputs. It does not modify",
            "article files, article figures, `paper_dorofeev_style`, the old determinant,",
            "`src/my_project/analytic/formulas.py`, old solvers, Gmsh/CalculiX workflows,",
            "baseline results, or the FEM physical model.",
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
    print("Descendant-6 localization audit")
    print(f"reading selected points: {INPUT_POINTS_CSV}")
    points = read_points()
    rows = compute_rows(points)
    write_csv(rows)
    write_report(rows)
    print(f"saved CSV: {OUTPUT_CSV}")
    print(f"saved report: {OUTPUT_REPORT}")
    return {"csv": OUTPUT_CSV, "report": OUTPUT_REPORT, "rows": rows}


if __name__ == "__main__":
    main()
