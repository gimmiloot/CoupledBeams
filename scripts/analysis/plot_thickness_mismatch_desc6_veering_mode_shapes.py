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

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402
from scripts.analysis import plot_thickness_mismatch_branch_shapes_vs_eta as shape_plotter  # noqa: E402


BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5
DESCENDANT_ID = 6

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.001
NUM_SORTED_ROOTS = 8
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401
MIN_DISTINCT_MU_SEPARATION = 0.02

PAIRS = ((4, 5), (5, 6))
OUTPUT_DIR = REPO_ROOT / "results" / "thickness_mismatch_desc6_veering_mode_shapes_beta15_eps0p0025_eta0p5"
OUTPUT_CSV = OUTPUT_DIR / "desc6_veering_points.csv"
OUTPUT_REPORT = OUTPUT_DIR / "desc6_veering_points_report.md"
AUDIT_CSV = REPO_ROOT / "results" / "eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.csv"
AUDIT_REPORT = REPO_ROOT / "results" / "eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.md"
AUDIT_PNG = REPO_ROOT / "results" / "eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.png"

LEGACY_PLOTS = {
    "near_4_5": OUTPUT_DIR / "desc6_near_4_5_veering.png",
    "near_5_6": OUTPUT_DIR / "desc6_near_5_6_veering.png",
}

CSV_FIELDNAMES = [
    "veering_case",
    "rank_within_case",
    "mu",
    "lambda_desc6",
    "descendant_id",
    "sorted_position",
    "nearest_sorted_pair",
    "local_gap",
    "tracking_mac",
    "tracking_warning",
    "plot_file",
    "tracking_status",
    "source_audit_candidate_type",
    "source_audit_pair_tracking_status",
    "legacy_plot_file",
    "legacy_plot_matched",
]


@dataclass(frozen=True)
class NearApproachPoint:
    veering_case: str
    pair_label: str
    rank_within_case: int
    mu: float
    local_gap: float
    plot_stem: str


def display_path(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def mu_values() -> np.ndarray:
    values = np.arange(MU_MIN, MU_MAX + 0.5 * MU_STEP, MU_STEP, dtype=float)
    values[0] = MU_MIN
    values[-1] = MU_MAX
    return np.unique(np.round(values, 12))


def solve_sorted_roots(mu_grid: np.ndarray) -> np.ndarray:
    beta_rad = float(np.deg2rad(BETA_DEG))
    roots = np.full((len(mu_grid), NUM_SORTED_ROOTS), np.nan, dtype=float)
    for idx, mu in enumerate(mu_grid):
        roots[idx] = find_first_n_roots_eta(
            beta_rad,
            float(mu),
            EPSILON,
            ETA,
            NUM_SORTED_ROOTS,
            Lmax0=ROOT_LMAX0,
            scan_step=ROOT_SCAN_STEP,
        )
        if np.any(~np.isfinite(roots[idx, : max(max(pair) for pair in PAIRS)])):
            raise RuntimeError(f"Missing sorted roots at mu={float(mu):g}.")
        if idx == 0 or (idx + 1) % 100 == 0 or idx + 1 == len(mu_grid):
            print(f"computed sorted roots for {idx + 1}/{len(mu_grid)} mu values")
    return roots


def strict_local_minima(mu_grid: np.ndarray, gaps: np.ndarray) -> list[tuple[float, float]]:
    minima: list[tuple[float, float]] = []
    for idx in range(1, len(mu_grid) - 1):
        left = float(gaps[idx - 1])
        mid = float(gaps[idx])
        right = float(gaps[idx + 1])
        if mid <= left and mid <= right and (mid < left or mid < right):
            minima.append((float(mu_grid[idx]), mid))
    return minima


def select_distinct_minima(
    mu_grid: np.ndarray,
    roots: np.ndarray,
    pair: tuple[int, int],
) -> tuple[list[NearApproachPoint], int, int]:
    left_sorted, right_sorted = pair
    gaps = roots[:, right_sorted - 1] - roots[:, left_sorted - 1]
    local_minima = strict_local_minima(mu_grid, gaps)
    selected: list[tuple[float, float]] = []
    filtered_duplicates = 0

    for mu, gap in sorted(local_minima, key=lambda item: item[1]):
        if any(abs(mu - existing_mu) < MIN_DISTINCT_MU_SEPARATION for existing_mu, _gap in selected):
            filtered_duplicates += 1
            continue
        selected.append((mu, gap))
        if len(selected) == 3:
            break

    case = f"near_{left_sorted}_{right_sorted}"
    stem_pair = f"{left_sorted}{right_sorted}"
    points = [
        NearApproachPoint(
            veering_case=case,
            pair_label=f"{left_sorted}-{right_sorted}",
            rank_within_case=rank,
            mu=float(mu),
            local_gap=float(gap),
            plot_stem=f"desc6_near{stem_pair}_rank{rank}",
        )
        for rank, (mu, gap) in enumerate(selected, start=1)
    ]
    return points, len(local_minima), filtered_duplicates


def read_audit_rows() -> dict[str, dict[str, str]]:
    if not AUDIT_CSV.exists():
        return {}
    with AUDIT_CSV.open("r", newline="", encoding="utf-8") as handle:
        return {row["pair"]: row for row in csv.DictReader(handle)}


def build_title(point: NearApproachPoint) -> str:
    return (
        f"Thickness mismatch descendant 6 near {point.pair_label}: "
        f"beta = {BETA_DEG:g} deg, epsilon = {EPSILON:g}, "
        f"eta = {ETA:g}, mu = {point.mu:g}"
    )


def run_shape_plotter(points: list[NearApproachPoint]) -> list[dict[str, str]]:
    plotter_args = [
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
        "--num-sorted-roots",
        str(NUM_SORTED_ROOTS),
        "--num-tracked-branches",
        str(NUM_SORTED_ROOTS),
        "--root-scan-step",
        str(ROOT_SCAN_STEP),
        "--root-lmax0",
        str(ROOT_LMAX0),
        "--num-shape-samples",
        str(NUM_SHAPE_SAMPLES),
        "--output-dir",
        str(OUTPUT_DIR),
        "--plot-stems",
        *[point.plot_stem for point in points],
        "--plot-titles",
        *[build_title(point) for point in points],
        "--skip-per-mu-reports",
        "--summary-kind",
        "veering",
        "--veering-pairs",
        *[point.pair_label for point in points],
        "--summary-csv",
        OUTPUT_CSV.name,
        "--summary-report",
        OUTPUT_REPORT.name,
    ]
    shape_plotter.main(plotter_args)
    with OUTPUT_CSV.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def legacy_match(point: NearApproachPoint, audit_rows: dict[str, dict[str, str]]) -> tuple[Path | None, bool]:
    legacy = LEGACY_PLOTS.get(point.veering_case)
    if legacy is None:
        return None, False
    audit_row = audit_rows.get(point.pair_label, {})
    try:
        audit_mu = float(audit_row.get("mu_candidate", "nan"))
    except ValueError:
        audit_mu = float("nan")
    matched = (
        bool(legacy.exists())
        and point.rank_within_case == 1
        and np.isfinite(audit_mu)
        and abs(float(point.mu) - audit_mu) <= 0.5 * MU_STEP + 1e-12
    )
    return legacy, matched


def build_final_rows(
    points: list[NearApproachPoint],
    plotter_rows: list[dict[str, str]],
    audit_rows: dict[str, dict[str, str]],
) -> list[dict[str, str]]:
    by_mu = {round(float(row["mu"]), 12): row for row in plotter_rows}
    rows: list[dict[str, str]] = []
    for point in points:
        plotter_row = by_mu[round(float(point.mu), 12)]
        audit_row = audit_rows.get(point.pair_label, {})
        legacy, matched = legacy_match(point, audit_rows)
        rows.append(
            {
                "veering_case": point.veering_case,
                "rank_within_case": str(point.rank_within_case),
                "mu": f"{point.mu:.12g}",
                "lambda_desc6": plotter_row["lambda_desc6"],
                "descendant_id": str(DESCENDANT_ID),
                "sorted_position": plotter_row["sorted_position"],
                "nearest_sorted_pair": point.pair_label,
                "local_gap": f"{point.local_gap:.12g}",
                "tracking_mac": plotter_row.get("tracking_mac", ""),
                "tracking_warning": plotter_row.get("tracking_warning", ""),
                "plot_file": display_path(OUTPUT_DIR / f"{point.plot_stem}.png"),
                "tracking_status": plotter_row.get("tracking_status", ""),
                "source_audit_candidate_type": str(audit_row.get("candidate_type", "")),
                "source_audit_pair_tracking_status": str(audit_row.get("tracking_status", "")),
                "legacy_plot_file": display_path(legacy) if legacy is not None else "",
                "legacy_plot_matched": "yes" if matched else "no",
            }
        )
    return rows


def write_final_csv(rows: list[dict[str, str]]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def write_report(
    rows: list[dict[str, str]],
    *,
    local_minima_counts: dict[str, int],
    filtered_duplicates: dict[str, int],
    legacy_existing_before: dict[str, bool],
    audit_rows: dict[str, dict[str, str]],
) -> None:
    warning_rows = [row for row in rows if row["tracking_warning"] == "yes"]
    existing_lines: list[str] = []
    for case, existed in legacy_existing_before.items():
        legacy = LEGACY_PLOTS[case]
        matches = [row for row in rows if row["veering_case"] == case and row["rank_within_case"] == "1"]
        match_text = "matched rank 1" if matches and matches[0]["legacy_plot_matched"] == "yes" else "not matched"
        status = "found before this run" if existed else "not found before this run"
        existing_lines.append(f"- `{display_path(legacy)}`: {status}; {match_text}.")

    added_missing = [
        row
        for row in rows
        if row["veering_case"] in {"near_4_5", "near_5_6"} and row["rank_within_case"] in {"2", "3"}
    ]

    lines = [
        "# Descendant 6 Veering Mode-Shape Points",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- descendant id: {DESCENDANT_ID}",
        f"- dense sorted-gap search: mu={MU_MIN:g}..{MU_MAX:g}, step={MU_STEP:g}",
        f"- distinct-local-minimum separation: {MIN_DISTINCT_MU_SEPARATION:g}",
        f"- source audit CSV: `{display_path(AUDIT_CSV)}`",
        f"- source audit report: `{display_path(AUDIT_REPORT)}`",
        f"- source audit PNG: `{display_path(AUDIT_PNG)}`",
        "",
        "This is an analytic-only Euler-Bernoulli thickness-mismatch diagnostic.",
        "It does not add Timoshenko curves, FEM markers, Gmsh, CalculiX, or 3D FEM runs.",
        "",
        "## Selection Method",
        "",
        "For each adjacent sorted pair, the script recomputes sorted EB roots on the dense",
        "audit-style grid, forms the sorted adjacent gap, keeps strict local minima only,",
        "then ranks distinct minima by increasing local gap. Descendant 6 shapes are then",
        "tracked from `mu=0`; sorted position is recorded only as metadata.",
        "",
        "## Existing And Added Plots",
        "",
    ]
    lines.extend(existing_lines)
    lines.append("- Standardized rank-named output files were written for all six selected points.")
    lines.append("- The four missing rank-2/rank-3 plots added by this run are:")
    for row in added_missing:
        lines.append(f"  - `{row['plot_file']}`")

    for case, title in (("near_4_5", "Pair 4-5"), ("near_5_6", "Pair 5-6")):
        case_rows = [row for row in rows if row["veering_case"] == case]
        lines.extend(["", f"## Selected Points: {title}", ""])
        lines.append(f"- strict local minima found before distinct filtering: {local_minima_counts.get(case, 0)}")
        lines.append(f"- near-duplicate local minima filtered: {filtered_duplicates.get(case, 0)}")
        if len(case_rows) < 3:
            lines.append(f"- WARNING: only {len(case_rows)} distinct local minima were found.")
        lines.append("")
        lines.append("| rank | mu | local gap | Lambda desc 6 | sorted position | tracking warning | plot |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for row in case_rows:
            lines.append(
                "| "
                + " | ".join(
                    [
                        row["rank_within_case"],
                        row["mu"],
                        row["local_gap"],
                        row["lambda_desc6"],
                        row["sorted_position"],
                        row["tracking_warning"],
                        f"`{row['plot_file']}`",
                    ]
                )
                + " |"
            )

    lines.extend(["", "## Source Audit Context", ""])
    for pair_label in ("4-5", "5-6"):
        audit_row = audit_rows.get(pair_label, {})
        if not audit_row:
            lines.append(f"- pair {pair_label}: no matching row found in the source audit CSV.")
            continue
        lines.append(
            f"- pair {pair_label}: audit candidate mu={audit_row.get('mu_candidate', '')}, "
            f"candidate_type={audit_row.get('candidate_type', '')}, "
            f"tracking_status={audit_row.get('tracking_status', '')}, "
            f"abs_gap={audit_row.get('abs_gap', '')}."
        )
    lines.append(
        "These are pair-level crossing-audit statuses; the six descendant-6 shape rows above "
        "use their own adjacent-step shape-MAC tracking metadata."
    )

    lines.extend(
        [
            "",
            "## Tracking Summary",
            "",
            f"- tracking warning rows: {len(warning_rows)}",
        ]
    )
    if warning_rows:
        for row in warning_rows:
            lines.append(
                f"  - {row['veering_case']} rank {row['rank_within_case']}, "
                f"mu={row['mu']}: {row['tracking_status']}"
            )
    else:
        lines.append("- no tracking warnings in the six descendant-6 plotted rows.")
    lines.extend(
        [
            "",
            "## Output Files",
            "",
            f"- CSV: `{display_path(OUTPUT_CSV)}`",
            f"- report: `{display_path(OUTPUT_REPORT)}`",
            "",
            "This diagnostic writes only results files. It does not modify article files,",
            "article figures, the old determinant, `src/my_project/analytic/formulas.py`,",
            "old solvers, Gmsh/CalculiX workflows, or the FEM physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, object]:
    print("Selecting descendant-6 veering mode-shape points")
    print(f"beta={BETA_DEG:g} deg, epsilon={EPSILON:g}, eta={ETA:g}")
    legacy_existing_before = {case: path.exists() for case, path in LEGACY_PLOTS.items()}
    audit_rows = read_audit_rows()
    mu_grid = mu_values()
    roots = solve_sorted_roots(mu_grid)

    points: list[NearApproachPoint] = []
    local_minima_counts: dict[str, int] = {}
    filtered_duplicates: dict[str, int] = {}
    for pair in PAIRS:
        pair_points, n_local_minima, n_filtered = select_distinct_minima(mu_grid, roots, pair)
        case = f"near_{pair[0]}_{pair[1]}"
        local_minima_counts[case] = n_local_minima
        filtered_duplicates[case] = n_filtered
        points.extend(pair_points)
        print(f"{case}: selected {len(pair_points)} of {n_local_minima} local minima")
        for point in pair_points:
            print(f"  rank {point.rank_within_case}: mu={point.mu:.12g}, gap={point.local_gap:.12g}")

    plotter_rows = run_shape_plotter(points)
    final_rows = build_final_rows(points, plotter_rows, audit_rows)
    write_final_csv(final_rows)
    write_report(
        final_rows,
        local_minima_counts=local_minima_counts,
        filtered_duplicates=filtered_duplicates,
        legacy_existing_before=legacy_existing_before,
        audit_rows=audit_rows,
    )
    print(f"saved ranked CSV: {OUTPUT_CSV}")
    print(f"saved ranked report: {OUTPUT_REPORT}")
    return {"points": points, "csv": OUTPUT_CSV, "report": OUTPUT_REPORT}


if __name__ == "__main__":
    main()
