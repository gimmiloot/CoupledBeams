from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Iterable, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
TRACKED_CSV = RESULTS_DIR / "thickness_mismatch_lambda_eta_beta15_eps0p0025_tracked.csv"
SORTED_CSV = RESULTS_DIR / "thickness_mismatch_lambda_eta_beta15_eps0p0025_sorted_for_tracking.csv"
TRACKED_PNG = RESULTS_DIR / "thickness_mismatch_lambda_eta_beta15_eps0p0025_tracked.png"
REPORT_MD = RESULTS_DIR / "thickness_mismatch_lambda_eta_tracking_report.md"

BETA_DEG = 15.0
EPSILON = 0.0025
MU_VALUES = (0.0, 0.3, 0.6)
ETA_VALUES = np.round(np.arange(-0.15, 0.1500001, 0.01), 10)
NUM_TRACKED_BRANCHES = 6
NUM_SORTED_ROOTS = 10
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0

TRACKED_FIELDNAMES = [
    "beta_deg",
    "epsilon",
    "mu",
    "eta",
    "branch_index_from_eta0",
    "Lambda_tracked",
    "nearest_sorted_root_index",
    "nearest_sorted_Lambda",
    "abs_diff_to_nearest_sorted",
    "tracking_direction",
    "tracking_step_status",
    "min_gap_to_other_sorted_roots",
    "distance_from_previous",
    "second_nearest_distance",
    "assignment_margin",
]
SORTED_FIELDNAMES = [
    "beta_deg",
    "epsilon",
    "mu",
    "eta",
    "sorted_root_index",
    "Lambda_sorted",
]


@dataclass(frozen=True)
class AssignmentInfo:
    branch_index: int
    root_index: int
    distance_from_previous: float
    second_nearest_distance: float

    @property
    def assignment_margin(self) -> float:
        return self.second_nearest_distance - self.distance_from_previous


def write_csv(path: Path, rows: Iterable[dict[str, float | int | str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def roots_for(mu: float, eta: float) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        float(eta),
        NUM_SORTED_ROOTS,
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots)):
        raise RuntimeError(f"Missing roots for mu={mu:g}, eta={eta:g}.")
    return roots


def sorted_roots_by_mu_eta() -> dict[float, dict[float, np.ndarray]]:
    by_mu: dict[float, dict[float, np.ndarray]] = {}
    for mu in MU_VALUES:
        by_mu[float(mu)] = {}
        for eta in ETA_VALUES:
            by_mu[float(mu)][float(eta)] = roots_for(float(mu), float(eta))
    return by_mu


def unique_nearest_assignment(previous: np.ndarray, candidates: np.ndarray) -> list[AssignmentInfo]:
    costs = np.abs(previous[:, None] - candidates[None, :])
    pairs: list[tuple[float, int, int]] = [
        (float(costs[row, col]), int(row), int(col))
        for row in range(costs.shape[0])
        for col in range(costs.shape[1])
    ]
    pairs.sort(key=lambda item: item[0])

    assigned_rows: set[int] = set()
    assigned_cols: set[int] = set()
    assignment: dict[int, int] = {}
    for _distance, row, col in pairs:
        if row in assigned_rows or col in assigned_cols:
            continue
        assignment[row] = col
        assigned_rows.add(row)
        assigned_cols.add(col)
        if len(assignment) == len(previous):
            break

    if len(assignment) != len(previous):
        raise RuntimeError("Could not assign all tracked branches to current sorted roots.")

    info: list[AssignmentInfo] = []
    for row in range(len(previous)):
        col = assignment[row]
        sorted_distances = np.sort(costs[row])
        second = float(sorted_distances[1]) if len(sorted_distances) > 1 else np.inf
        info.append(
            AssignmentInfo(
                branch_index=row + 1,
                root_index=col + 1,
                distance_from_previous=float(costs[row, col]),
                second_nearest_distance=second,
            )
        )
    return info


def row_for_tracked(
    *,
    mu: float,
    eta: float,
    branch_index: int,
    value: float,
    roots: np.ndarray,
    direction: str,
    status: str,
    assignment: AssignmentInfo | None,
) -> dict[str, float | int | str]:
    distances = np.abs(roots - float(value))
    nearest_idx = int(np.argmin(distances))
    other_distances = np.delete(distances, nearest_idx)
    min_gap = float(np.min(other_distances)) if len(other_distances) else np.inf
    distance_from_previous = np.nan if assignment is None else float(assignment.distance_from_previous)
    second_nearest = np.nan if assignment is None else float(assignment.second_nearest_distance)
    assignment_margin = np.nan if assignment is None else float(assignment.assignment_margin)
    return {
        "beta_deg": BETA_DEG,
        "epsilon": EPSILON,
        "mu": float(mu),
        "eta": float(eta),
        "branch_index_from_eta0": int(branch_index),
        "Lambda_tracked": float(value),
        "nearest_sorted_root_index": int(nearest_idx) + 1,
        "nearest_sorted_Lambda": float(roots[nearest_idx]),
        "abs_diff_to_nearest_sorted": float(distances[nearest_idx]),
        "tracking_direction": direction,
        "tracking_step_status": status,
        "min_gap_to_other_sorted_roots": min_gap,
        "distance_from_previous": distance_from_previous,
        "second_nearest_distance": second_nearest,
        "assignment_margin": assignment_margin,
    }


def track_one_mu(mu: float, roots_by_eta: dict[float, np.ndarray]) -> tuple[list[dict[str, float | int | str]], np.ndarray]:
    zero_idx = int(np.where(np.isclose(ETA_VALUES, 0.0, rtol=0.0, atol=1e-12))[0][0])
    tracked = np.full((NUM_TRACKED_BRANCHES, len(ETA_VALUES)), np.nan, dtype=float)
    tracked[:, zero_idx] = roots_by_eta[0.0][:NUM_TRACKED_BRANCHES]

    rows: list[dict[str, float | int | str]] = []
    for branch_idx, value in enumerate(tracked[:, zero_idx], start=1):
        rows.append(
            row_for_tracked(
                mu=mu,
                eta=0.0,
                branch_index=branch_idx,
                value=float(value),
                roots=roots_by_eta[0.0],
                direction="seed",
                status="seed_eta0",
                assignment=None,
            )
        )

    previous = tracked[:, zero_idx].copy()
    for col in range(zero_idx + 1, len(ETA_VALUES)):
        eta = float(ETA_VALUES[col])
        roots = roots_by_eta[eta]
        assignments = unique_nearest_assignment(previous, roots)
        current = np.full(NUM_TRACKED_BRANCHES, np.nan, dtype=float)
        for assignment in assignments:
            value = float(roots[assignment.root_index - 1])
            current[assignment.branch_index - 1] = value
            status = "ok" if assignment.assignment_margin > 1e-5 else "ambiguous"
            rows.append(
                row_for_tracked(
                    mu=mu,
                    eta=eta,
                    branch_index=assignment.branch_index,
                    value=value,
                    roots=roots,
                    direction="positive",
                    status=status,
                    assignment=assignment,
                )
            )
        tracked[:, col] = current
        previous = current

    previous = tracked[:, zero_idx].copy()
    for col in range(zero_idx - 1, -1, -1):
        eta = float(ETA_VALUES[col])
        roots = roots_by_eta[eta]
        assignments = unique_nearest_assignment(previous, roots)
        current = np.full(NUM_TRACKED_BRANCHES, np.nan, dtype=float)
        for assignment in assignments:
            value = float(roots[assignment.root_index - 1])
            current[assignment.branch_index - 1] = value
            status = "ok" if assignment.assignment_margin > 1e-5 else "ambiguous"
            rows.append(
                row_for_tracked(
                    mu=mu,
                    eta=eta,
                    branch_index=assignment.branch_index,
                    value=value,
                    roots=roots,
                    direction="negative",
                    status=status,
                    assignment=assignment,
                )
            )
        tracked[:, col] = current
        previous = current

    rows.sort(key=lambda row: (float(row["eta"]), int(row["branch_index_from_eta0"])))
    return rows, tracked


def sorted_rows(sorted_roots: dict[float, dict[float, np.ndarray]]) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for mu in MU_VALUES:
        for eta in ETA_VALUES:
            roots = sorted_roots[float(mu)][float(eta)]
            for idx, value in enumerate(roots[:NUM_TRACKED_BRANCHES], start=1):
                rows.append(
                    {
                        "beta_deg": BETA_DEG,
                        "epsilon": EPSILON,
                        "mu": float(mu),
                        "eta": float(eta),
                        "sorted_root_index": int(idx),
                        "Lambda_sorted": float(value),
                    }
                )
    return rows


def plot_tracked(
    *,
    sorted_roots: dict[float, dict[float, np.ndarray]],
    tracked_by_mu: dict[float, np.ndarray],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, len(MU_VALUES), figsize=(10.8, 3.8), sharey=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, mu in zip(axes, MU_VALUES, strict=True):
        tracked = tracked_by_mu[float(mu)]
        sorted_grid = np.column_stack([sorted_roots[float(mu)][float(eta)][:NUM_TRACKED_BRANCHES] for eta in ETA_VALUES])
        for branch_idx in range(NUM_TRACKED_BRANCHES):
            color = colors[branch_idx % len(colors)]
            ax.plot(ETA_VALUES, sorted_grid[branch_idx], color=color, lw=0.8, ls=":", alpha=0.55)
            ax.plot(
                ETA_VALUES,
                tracked[branch_idx],
                color=color,
                lw=1.7,
                label=f"branch {branch_idx + 1}" if ax is axes[0] else "_nolegend_",
            )
        ax.axvline(0.0, color="0.55", lw=0.8, ls="--")
        ax.set_title(rf"$\mu={mu:g}$", fontsize=11)
        ax.set_xlabel(r"$\eta$")
        ax.grid(True, color="0.88", linewidth=0.6)
    axes[0].set_ylabel(r"$\Lambda$")
    axes[0].legend(loc="upper left", fontsize=8, frameon=False)
    fig.suptitle(
        rf"Tracked Lambda(eta), beta={BETA_DEG:g}$^\circ$, epsilon={EPSILON:g}; dotted = same sorted index",
        fontsize=12,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def compute_same_index_diffs(
    sorted_roots: dict[float, dict[float, np.ndarray]],
    tracked_by_mu: dict[float, np.ndarray],
) -> list[dict[str, float | int]]:
    diffs: list[dict[str, float | int]] = []
    for mu in MU_VALUES:
        tracked = tracked_by_mu[float(mu)]
        for branch_idx in range(NUM_TRACKED_BRANCHES):
            for col, eta in enumerate(ETA_VALUES):
                sorted_value = sorted_roots[float(mu)][float(eta)][branch_idx]
                tracked_value = tracked[branch_idx, col]
                diffs.append(
                    {
                        "mu": float(mu),
                        "eta": float(eta),
                        "branch_index": int(branch_idx) + 1,
                        "abs_diff": abs(float(sorted_value) - float(tracked_value)),
                    }
                )
    return diffs


def max_mu0_symmetry_error(tracked_by_mu: dict[float, np.ndarray]) -> float:
    tracked = tracked_by_mu[0.0]
    error = 0.0
    for col, eta in enumerate(ETA_VALUES):
        mirror_cols = np.where(np.isclose(ETA_VALUES, -float(eta), rtol=0.0, atol=1e-12))[0]
        if len(mirror_cols) != 1:
            continue
        mirror_col = int(mirror_cols[0])
        error = max(error, float(np.max(np.abs(tracked[:, col] - tracked[:, mirror_col]))))
    return error


def small_eta_jump(kind_grid: np.ndarray, branch_index: int) -> tuple[float, float, float]:
    small_eta = np.array([-0.1, -0.05, 0.0, 0.05, 0.1], dtype=float)
    indices = [int(np.where(np.isclose(ETA_VALUES, value, atol=1e-12))[0][0]) for value in small_eta]
    values = kind_grid[int(branch_index) - 1, indices]
    jumps = np.abs(np.diff(values))
    jump_idx = int(np.argmax(jumps))
    return float(jumps[jump_idx]), float(small_eta[jump_idx]), float(small_eta[jump_idx + 1])


def build_report(
    *,
    sorted_roots: dict[float, dict[float, np.ndarray]],
    tracked_rows: Sequence[dict[str, float | int | str]],
    tracked_by_mu: dict[float, np.ndarray],
    same_index_diffs: Sequence[dict[str, float | int]],
    output: Path,
) -> dict[str, float | int | str]:
    max_diff_row = max(same_index_diffs, key=lambda row: float(row["abs_diff"]))
    switches = [
        row
        for row in tracked_rows
        if int(row["nearest_sorted_root_index"]) != int(row["branch_index_from_eta0"])
    ]
    ambiguous = [row for row in tracked_rows if str(row["tracking_step_status"]) == "ambiguous"]
    symmetry_error = max_mu0_symmetry_error(tracked_by_mu)

    sorted_grid_mu03 = np.column_stack([sorted_roots[0.3][float(eta)][:NUM_TRACKED_BRANCHES] for eta in ETA_VALUES])
    sorted_jump, sorted_jump_left, sorted_jump_right = small_eta_jump(sorted_grid_mu03, branch_index=5)
    tracked_jump, tracked_jump_left, tracked_jump_right = small_eta_jump(tracked_by_mu[0.3], branch_index=5)
    jump_is_artifact = (
        sorted_jump > 5.0 * max(tracked_jump, 1e-12)
        or any(
            int(row["branch_index_from_eta0"]) == 5
            and float(row["mu"]) == 0.3
            and sorted_jump_left - 1e-12 <= float(row["eta"]) <= sorted_jump_right + 1e-12
            and int(row["nearest_sorted_root_index"]) != 5
            for row in tracked_rows
        )
    )

    lines = [
        "# Thickness-Mismatch Lambda(eta) Tracking Report",
        "",
        "## Setup",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- mu values: {', '.join(f'{value:g}' for value in MU_VALUES)}",
        f"- eta range: {float(ETA_VALUES[0]):g} .. {float(ETA_VALUES[-1]):g}",
        f"- tracked branches: first {NUM_TRACKED_BRANCHES} sorted roots at eta=0",
        "",
        "## Method",
        "",
        "Branches are seeded at eta=0. The script continues from eta=0 to",
        "+0.15 and from eta=0 to -0.15. Each step assigns previous branch values",
        "to the current sorted candidate roots by unique nearest-root matching.",
        "The CSV records the nearest sorted index and the assignment margin.",
        "",
        "## Outputs",
        "",
        f"- tracked CSV: `{TRACKED_CSV.relative_to(REPO_ROOT)}`",
        f"- sorted comparison CSV: `{SORTED_CSV.relative_to(REPO_ROOT)}`",
        f"- tracked plot: `{TRACKED_PNG.relative_to(REPO_ROOT)}`",
        "",
        "## Checks",
        "",
        f"- max mu=0 tracked symmetry error Lambda(eta)-Lambda(-eta): {symmetry_error:.6e}",
        f"- max same-index sorted-vs-tracked difference: {float(max_diff_row['abs_diff']):.6e}",
        (
            "- max difference location: "
            f"mu={float(max_diff_row['mu']):g}, "
            f"branch/root={int(max_diff_row['branch_index'])}, "
            f"eta={float(max_diff_row['eta']):g}"
        ),
        f"- rows where nearest sorted index differs from eta=0 branch index: {len(switches)}",
        f"- ambiguous nearest-root assignments: {len(ambiguous)}",
        "",
        "## Previous Sorted-Root Jump",
        "",
        (
            "For mu=0.3, sorted root 5 over eta = [-0.1, -0.05, 0, 0.05, 0.1] "
            f"has max adjacent jump {sorted_jump:.6e} between eta={sorted_jump_left:g} and {sorted_jump_right:g}."
        ),
        (
            "The tracked branch seeded as eta=0 branch 5 over the same eta set "
            f"has max adjacent jump {tracked_jump:.6e} between eta={tracked_jump_left:g} and {tracked_jump_right:g}."
        ),
        (
            "Conclusion: the previous root-5 jump is "
            + ("likely a sorted-root artifact." if jump_is_artifact else "not explained by sorted-index switching in this nearest-root diagnostic.")
        ),
        "",
        "## Suspected Sorted-Index Switches",
        "",
    ]

    if switches:
        lines.append("| mu | eta | branch from eta=0 | nearest sorted index | Lambda tracked | sorted Lambda |")
        lines.append("| --- | --- | --- | --- | --- | --- |")
        for row in switches[:24]:
            lines.append(
                "| "
                f"{float(row['mu']):g} | "
                f"{float(row['eta']):g} | "
                f"{int(row['branch_index_from_eta0'])} | "
                f"{int(row['nearest_sorted_root_index'])} | "
                f"{float(row['Lambda_tracked']):.8g} | "
                f"{float(row['nearest_sorted_Lambda']):.8g} |"
            )
        if len(switches) > 24:
            lines.append(f"| ... | ... | ... | ... | ... | ... |")
    else:
        lines.append("No sorted-index switches were detected for the tracked branches on this eta grid.")

    lines.extend(
        [
            "",
            "These are diagnostic observations only. The script does not change the",
            "baseline determinant, solvers, FEM model, article files, or article figures.",
            "",
        ]
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    return {
        "max_mu0_symmetry_error": symmetry_error,
        "max_same_index_diff": float(max_diff_row["abs_diff"]),
        "max_same_index_diff_mu": float(max_diff_row["mu"]),
        "max_same_index_diff_branch": int(max_diff_row["branch_index"]),
        "max_same_index_diff_eta": float(max_diff_row["eta"]),
        "switch_count": len(switches),
        "ambiguous_count": len(ambiguous),
        "sorted_mu03_root5_jump": sorted_jump,
        "tracked_mu03_branch5_jump": tracked_jump,
        "previous_jump_is_artifact": "yes" if jump_is_artifact else "no",
    }


def main() -> dict[str, Path | float | int | str]:
    sorted_roots = sorted_roots_by_mu_eta()
    sorted_comparison_rows = sorted_rows(sorted_roots)
    tracked_rows: list[dict[str, float | int | str]] = []
    tracked_by_mu: dict[float, np.ndarray] = {}

    for mu in MU_VALUES:
        rows, tracked = track_one_mu(float(mu), sorted_roots[float(mu)])
        tracked_rows.extend(rows)
        tracked_by_mu[float(mu)] = tracked

    write_csv(SORTED_CSV, sorted_comparison_rows, SORTED_FIELDNAMES)
    write_csv(TRACKED_CSV, tracked_rows, TRACKED_FIELDNAMES)
    plot_tracked(sorted_roots=sorted_roots, tracked_by_mu=tracked_by_mu, output=TRACKED_PNG)

    same_index_diffs = compute_same_index_diffs(sorted_roots, tracked_by_mu)
    report = build_report(
        sorted_roots=sorted_roots,
        tracked_rows=tracked_rows,
        tracked_by_mu=tracked_by_mu,
        same_index_diffs=same_index_diffs,
        output=REPORT_MD,
    )

    print(f"saved tracked CSV: {TRACKED_CSV}")
    print(f"saved sorted comparison CSV: {SORTED_CSV}")
    print(f"saved tracked plot: {TRACKED_PNG}")
    print(f"saved tracking report: {REPORT_MD}")
    print(f"max mu=0 tracked symmetry error: {float(report['max_mu0_symmetry_error']):.6e}")
    print(f"max same-index sorted-vs-tracked diff: {float(report['max_same_index_diff']):.6e}")
    print(f"sorted-index switch rows: {int(report['switch_count'])}")
    print(
        "mu=0.3 root/branch 5 small-eta jump: "
        f"sorted={float(report['sorted_mu03_root5_jump']):.6e}, "
        f"tracked={float(report['tracked_mu03_branch5_jump']):.6e}, "
        f"artifact={report['previous_jump_is_artifact']}"
    )

    return {
        "tracked_csv": TRACKED_CSV,
        "sorted_csv": SORTED_CSV,
        "tracked_png": TRACKED_PNG,
        "report_md": REPORT_MD,
        **report,
    }


if __name__ == "__main__":
    main()
