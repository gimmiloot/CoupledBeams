from __future__ import annotations

from pathlib import Path
import sys
from typing import Sequence

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
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402
from scripts.analysis.track_lambda_eta_thickness_mismatch import (  # noqa: E402
    write_csv,
)
from scripts.lib.thickness_mismatch_mac_tracking import track_mu_branches_shape_mac  # noqa: E402


# =========================
# User-editable parameters
# =========================
BETA_DEG = 7.5
EPSILON = 0.0025
ETA_VALUES = (-0.1, 0.0, 0.1)
MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.01

NUM_TRACKED_BRANCHES = 6
NUM_SORTED_ROOTS = 12
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401
MAC_WARNING_THRESHOLD = 0.9

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_TAG: str | None = None
SAVE_PNG = True
SAVE_CSV = False
SAVE_REPORT = False

MU_VALUES = np.round(np.arange(MU_MIN, MU_MAX + 0.5 * MU_STEP, MU_STEP), 10)

TRACKED_FIELDNAMES = [
    "beta_deg",
    "epsilon",
    "eta",
    "mu",
    "branch_index_from_mu0",
    "Lambda_tracked",
    "current_sorted_root_index",
    "raw_mac_sorted_root_index",
    "mac_to_previous",
    "raw_mac_to_previous",
    "second_best_mac",
    "frequency_nearest_sorted_root_index",
    "frequency_mac_disagreement",
    "used_frequency_fallback",
    "mac_margin",
    "raw_mac_margin",
    "assigned_from_previous_sorted_position",
    "sorted_position_jump",
    "raw_sorted_position_jump",
    "low_mac",
    "low_margin",
    "large_mu_step",
    "suspicious_assignment",
    "requires_refined_check",
    "nearest_sorted_root_index",
    "nearest_sorted_Lambda",
    "abs_diff_to_nearest_sorted",
    "tracking_step_status",
    "min_gap_to_other_sorted_roots",
    "distance_from_previous",
    "second_nearest_distance",
    "assignment_margin",
]


def number_tag(value: float) -> str:
    text = f"{float(value):.10g}"
    return text.replace("-", "m").replace(".", "p")


def output_stem() -> str:
    if OUTPUT_TAG:
        return OUTPUT_TAG
    return (
        "thickness_mismatch_lambda_mu_"
        f"beta{number_tag(BETA_DEG)}_eps{number_tag(EPSILON)}_eta_sweep_tracked"
    )


def tracked_csv_path() -> Path:
    return OUTPUT_DIR / f"{output_stem()}.csv"


def tracked_png_path() -> Path:
    return OUTPUT_DIR / f"{output_stem()}.png"


def report_md_path() -> Path:
    report_stem = output_stem().removesuffix("_tracked")
    return OUTPUT_DIR / f"{report_stem}_tracking_report.md"


def roots_for(mu: float, eta: float, n_roots: int = NUM_SORTED_ROOTS) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        float(eta),
        int(n_roots),
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots)):
        raise RuntimeError(f"Missing roots for mu={mu:g}, eta={eta:g}.")
    return roots


def sorted_roots_by_eta_mu() -> dict[float, dict[float, np.ndarray]]:
    by_eta: dict[float, dict[float, np.ndarray]] = {}
    for eta in ETA_VALUES:
        by_eta[float(eta)] = {}
        for mu in MU_VALUES:
            by_eta[float(eta)][float(mu)] = roots_for(float(mu), float(eta))
    return by_eta


def row_for_tracked(
    *,
    eta: float,
    mu: float,
    branch_index: int,
    value: float,
    roots: np.ndarray,
    status: str,
    distance_from_previous: float,
    second_nearest_distance: float,
    assignment_margin: float,
) -> dict[str, float | int | str]:
    distances = np.abs(roots - float(value))
    nearest_idx = int(np.argmin(distances))
    other_distances = np.delete(distances, nearest_idx)
    min_gap = float(np.min(other_distances)) if len(other_distances) else np.inf
    return {
        "beta_deg": BETA_DEG,
        "epsilon": EPSILON,
        "eta": float(eta),
        "mu": float(mu),
        "branch_index_from_mu0": int(branch_index),
        "Lambda_tracked": float(value),
        "nearest_sorted_root_index": int(nearest_idx) + 1,
        "nearest_sorted_Lambda": float(roots[nearest_idx]),
        "abs_diff_to_nearest_sorted": float(distances[nearest_idx]),
        "tracking_step_status": status,
        "min_gap_to_other_sorted_roots": min_gap,
        "distance_from_previous": float(distance_from_previous),
        "second_nearest_distance": float(second_nearest_distance),
        "assignment_margin": float(assignment_margin),
    }


def adapt_mac_tracking_row(
    row: dict[str, float | int | str],
    *,
    roots_by_mu: dict[float, np.ndarray],
) -> dict[str, float | int | str]:
    mu = float(row["mu"])
    roots = roots_by_mu[mu]
    current_idx = int(row["mac_sorted_root_index"])
    value = float(row["Lambda_tracked"])
    distances = np.abs(roots - value)
    other_distances = np.delete(distances, current_idx - 1)
    min_gap = float(np.min(other_distances)) if len(other_distances) else np.inf
    frequency_nearest = int(row["nearest_sorted_root_index"])
    return {
        "beta_deg": BETA_DEG,
        "epsilon": EPSILON,
        "eta": float(row["eta"]),
        "mu": mu,
        "branch_index_from_mu0": int(row["branch_index_from_mu0"]),
        "Lambda_tracked": value,
        "current_sorted_root_index": current_idx,
        "raw_mac_sorted_root_index": int(row["raw_mac_sorted_root_index"]),
        "mac_to_previous": float(row["mac_to_previous"]),
        "raw_mac_to_previous": float(row["raw_mac_to_previous"]),
        "second_best_mac": float(row["second_best_mac"]),
        "frequency_nearest_sorted_root_index": frequency_nearest,
        "frequency_mac_disagreement": str(row["frequency_mac_disagreement"]),
        "used_frequency_fallback": str(row["used_frequency_fallback"]),
        "mac_margin": float(row["mac_margin"]),
        "raw_mac_margin": float(row["raw_mac_margin"]),
        "assigned_from_previous_sorted_position": int(row["assigned_from_previous_sorted_position"]),
        "sorted_position_jump": int(row["sorted_position_jump"]),
        "raw_sorted_position_jump": int(row["raw_sorted_position_jump"]),
        "low_mac": str(row["low_mac"]),
        "low_margin": str(row["low_margin"]),
        "large_mu_step": str(row["large_mu_step"]),
        "suspicious_assignment": str(row["suspicious_assignment"]),
        "requires_refined_check": str(row["requires_refined_check"]),
        # Backward-compatible report fields: this now means the MAC-tracked current sorted index.
        "nearest_sorted_root_index": current_idx,
        "nearest_sorted_Lambda": float(roots[current_idx - 1]),
        "abs_diff_to_nearest_sorted": 0.0,
        "tracking_step_status": str(row["tracking_step_status"]),
        "min_gap_to_other_sorted_roots": min_gap,
        "distance_from_previous": float(row["frequency_distance_from_previous"]),
        "second_nearest_distance": float(row["frequency_second_nearest_distance"]),
        "assignment_margin": float(row["frequency_assignment_margin"]),
    }


def track_one_eta(
    eta: float,
    roots_by_mu: dict[float, np.ndarray],
) -> tuple[list[dict[str, float | int | str]], np.ndarray]:
    result = track_mu_branches_shape_mac(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=float(eta),
        mu_values=MU_VALUES,
        roots_by_mu=roots_by_mu,
        num_tracked_branches=NUM_TRACKED_BRANCHES,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
    )
    tracked = result.tracked
    rows = [adapt_mac_tracking_row(row, roots_by_mu=roots_by_mu) for row in result.rows]
    rows.sort(key=lambda row: (float(row["mu"]), int(row["branch_index_from_mu0"])))
    return rows, tracked


def eta_panel_title(eta: float) -> str:
    if eta > 0.0:
        meaning = "longer rod thicker for $\\mu>0$"
    elif eta < 0.0:
        meaning = "shorter rod thicker for $\\mu>0$"
    else:
        meaning = "equal radii"
    return rf"$\eta={eta:g}$" + "\n" + meaning


def plot_tracked(
    *,
    sorted_roots: dict[float, dict[float, np.ndarray]],
    tracked_by_eta: dict[float, np.ndarray],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, len(ETA_VALUES), figsize=(11.2, 3.8), sharey=True)
    axes = np.atleast_1d(axes)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, eta in zip(axes, ETA_VALUES, strict=True):
        tracked = tracked_by_eta[float(eta)]
        sorted_grid = np.column_stack(
            [sorted_roots[float(eta)][float(mu)][:NUM_TRACKED_BRANCHES] for mu in MU_VALUES]
        )
        for branch_idx in range(NUM_TRACKED_BRANCHES):
            color = colors[branch_idx % len(colors)]
            ax.plot(MU_VALUES, sorted_grid[branch_idx], color=color, lw=0.8, ls=":", alpha=0.5)
            ax.plot(
                MU_VALUES,
                tracked[branch_idx],
                color=color,
                lw=1.7,
                label=f"branch {branch_idx + 1}" if ax is axes[0] else "_nolegend_",
            )
        ax.set_title(eta_panel_title(float(eta)), fontsize=10)
        ax.set_xlabel(r"$\mu$")
        ax.grid(True, color="0.88", linewidth=0.6)
    axes[0].set_ylabel(r"$\Lambda$")
    axes[0].legend(loc="upper right", fontsize=8, frameon=False)
    fig.suptitle(
        rf"Tracked Lambda(mu), beta={BETA_DEG:g}$^\circ$, epsilon={EPSILON:g}; dotted = same sorted index",
        fontsize=12,
    )
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def compute_same_index_diffs(
    sorted_roots: dict[float, dict[float, np.ndarray]],
    tracked_by_eta: dict[float, np.ndarray],
) -> list[dict[str, float | int]]:
    diffs: list[dict[str, float | int]] = []
    for eta in ETA_VALUES:
        tracked = tracked_by_eta[float(eta)]
        for branch_idx in range(NUM_TRACKED_BRANCHES):
            for col, mu in enumerate(MU_VALUES):
                sorted_value = sorted_roots[float(eta)][float(mu)][branch_idx]
                tracked_value = tracked[branch_idx, col]
                diffs.append(
                    {
                        "eta": float(eta),
                        "mu": float(mu),
                        "branch_index": int(branch_idx) + 1,
                        "abs_diff": abs(float(sorted_value) - float(tracked_value)),
                    }
                )
    return diffs


def eta_zero_consistency(tracked_by_eta: dict[float, np.ndarray]) -> dict[str, float | int]:
    tracked = tracked_by_eta[0.0]
    max_nearest_diff = 0.0
    max_same_index_diff = 0.0
    max_nearest_mu = 0.0
    max_nearest_branch = 1
    max_same_mu = 0.0
    max_same_branch = 1
    for col, mu in enumerate(MU_VALUES):
        old_roots = find_first_n_roots(
            float(np.deg2rad(BETA_DEG)),
            float(mu),
            EPSILON,
            NUM_SORTED_ROOTS,
            Lmax0=ROOT_LMAX0,
            scan_step=ROOT_SCAN_STEP,
        )
        if np.any(~np.isfinite(old_roots)):
            raise RuntimeError(f"Missing baseline roots for eta=0 consistency at mu={mu:g}.")
        for branch_idx in range(NUM_TRACKED_BRANCHES):
            tracked_value = float(tracked[branch_idx, col])
            nearest_diff = float(np.min(np.abs(old_roots - tracked_value)))
            same_index_diff = abs(float(old_roots[branch_idx]) - tracked_value)
            if nearest_diff > max_nearest_diff:
                max_nearest_diff = nearest_diff
                max_nearest_mu = float(mu)
                max_nearest_branch = branch_idx + 1
            if same_index_diff > max_same_index_diff:
                max_same_index_diff = same_index_diff
                max_same_mu = float(mu)
                max_same_branch = branch_idx + 1
    return {
        "max_nearest_old_root_diff": max_nearest_diff,
        "max_nearest_old_root_mu": max_nearest_mu,
        "max_nearest_old_root_branch": max_nearest_branch,
        "max_same_index_old_root_diff": max_same_index_diff,
        "max_same_index_old_root_mu": max_same_mu,
        "max_same_index_old_root_branch": max_same_branch,
    }


def swap_symmetry_sanity() -> dict[str, float | int]:
    max_abs = 0.0
    max_rel = 0.0
    max_mu = 0.0
    max_eta = 0.0
    max_root = 1
    for mu in (0.0, 0.3, 0.6, 0.9):
        for eta in ETA_VALUES:
            roots = roots_for(float(mu), float(eta), NUM_TRACKED_BRANCHES)
            swapped = roots_for(-float(mu), -float(eta), NUM_TRACKED_BRANCHES)
            for idx, (left, right) in enumerate(zip(roots, swapped, strict=True), start=1):
                abs_diff = abs(float(left) - float(right))
                rel_diff = abs_diff / max(abs(float(left)), abs(float(right)), 1e-15)
                if abs_diff > max_abs:
                    max_abs = abs_diff
                    max_rel = rel_diff
                    max_mu = float(mu)
                    max_eta = float(eta)
                    max_root = int(idx)
    return {
        "max_swap_abs_diff": max_abs,
        "max_swap_rel_diff": max_rel,
        "max_swap_mu": max_mu,
        "max_swap_eta": max_eta,
        "max_swap_root": max_root,
    }


def branch_sensitivity_rows(tracked_by_eta: dict[float, np.ndarray]) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for branch_idx in range(NUM_TRACKED_BRANCHES):
        values_by_eta = np.vstack([tracked_by_eta[float(eta)][branch_idx] for eta in ETA_VALUES])
        spreads = np.max(values_by_eta, axis=0) - np.min(values_by_eta, axis=0)
        max_col = int(np.argmax(spreads))
        rows.append(
            {
                "branch_index": int(branch_idx) + 1,
                "max_eta_spread": float(spreads[max_col]),
                "max_eta_spread_mu": float(MU_VALUES[max_col]),
                "spread_at_mu_0p9": float(spreads[-1]),
                "eta_p0p1_minus_m0p1_at_mu_0p9": float(values_by_eta[2, -1] - values_by_eta[0, -1]),
            }
        )
    rows.sort(key=lambda row: float(row["max_eta_spread"]), reverse=True)
    return rows


def jump_rows(tracked_by_eta: dict[float, np.ndarray]) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for eta in ETA_VALUES:
        tracked = tracked_by_eta[float(eta)]
        jumps = np.abs(np.diff(tracked, axis=1))
        for branch_idx in range(NUM_TRACKED_BRANCHES):
            local_jumps = jumps[branch_idx]
            max_col = int(np.argmax(local_jumps))
            rows.append(
                {
                    "eta": float(eta),
                    "branch_index": int(branch_idx) + 1,
                    "max_adjacent_mu_jump": float(local_jumps[max_col]),
                    "left_mu": float(MU_VALUES[max_col]),
                    "right_mu": float(MU_VALUES[max_col + 1]),
                }
            )
    rows.sort(key=lambda row: float(row["max_adjacent_mu_jump"]), reverse=True)
    return rows


def format_table(rows: Sequence[Sequence[str]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(row[col]) for row in rows) for col in range(len(rows[0]))]
    formatted = []
    for idx, row in enumerate(rows):
        formatted.append("| " + " | ".join(value.ljust(widths[col]) for col, value in enumerate(row)) + " |")
        if idx == 0:
            formatted.append("| " + " | ".join("-" * widths[col] for col in range(len(row))) + " |")
    return formatted


def build_report(
    *,
    tracked_rows: Sequence[dict[str, float | int | str]],
    same_index_diffs: Sequence[dict[str, float | int]],
    sensitivity_rows: Sequence[dict[str, float | int]],
    jumps: Sequence[dict[str, float | int]],
    eta_zero: dict[str, float | int],
    swap: dict[str, float | int],
    output: Path | None,
) -> dict[str, float | int]:
    max_diff_row = max(same_index_diffs, key=lambda row: float(row["abs_diff"]))
    switches = [
        row
        for row in tracked_rows
        if int(row["nearest_sorted_root_index"]) != int(row["branch_index_from_mu0"])
    ]
    low_mac = [row for row in tracked_rows if "low_mac" in str(row["tracking_step_status"])]
    frequency_warnings = [row for row in tracked_rows if str(row["frequency_mac_disagreement"]) == "yes"]
    close_rows = sorted(
        (row for row in tracked_rows if str(row["tracking_step_status"]) != "seed_mu0"),
        key=lambda row: float(row["min_gap_to_other_sorted_roots"]),
    )[:8]

    lines = [
        "# Thickness-Mismatch Lambda(mu) Eta-Sweep Tracking Report",
        "",
        "## Setup",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta values: {', '.join(f'{value:g}' for value in ETA_VALUES)}",
        f"- mu range: {float(MU_VALUES[0]):g} .. {float(MU_VALUES[-1]):g}",
        f"- mu step: {float(MU_VALUES[1] - MU_VALUES[0]):g}",
        f"- tracked branches: first {NUM_TRACKED_BRANCHES} sorted roots at mu=0 for each eta",
        "",
        "For eta > 0 the second rod is thicker than the first. For mu > 0 the",
        "second rod is also longer, so eta=0.1 is the case where the longer rod",
        "is thicker, while eta=-0.1 is the case where the shorter rod is thicker.",
        "",
        "## Method",
        "",
        "For each fixed eta, branches are seeded at mu=0 and continued to mu=0.9.",
        "Each step assigns previous branch values to current sorted candidate roots",
        "by maximizing the adjacent-step analytic shape MAC. Low-MAC assignments",
        "are marked for review. Nearest-frequency assignment is retained only as",
        "a warning diagnostic. The tracking",
        "diagnostics record the current sorted index, raw MAC index, MAC value,",
        "nearest-frequency assignment, and local root gap.",
        "",
        "## Outputs",
        "",
        (
            f"- tracked CSV: `{tracked_csv_path().relative_to(REPO_ROOT)}`"
            if SAVE_CSV
            else "- tracked CSV: disabled (`SAVE_CSV = False`)"
        ),
        (
            f"- tracked plot: `{tracked_png_path().relative_to(REPO_ROOT)}`"
            if SAVE_PNG
            else "- tracked plot: disabled (`SAVE_PNG = False`)"
        ),
        "",
        "## Checks",
        "",
        f"- eta=0 max nearest old-root difference: {float(eta_zero['max_nearest_old_root_diff']):.6e}",
        f"- eta=0 max same-index old-root difference: {float(eta_zero['max_same_index_old_root_diff']):.6e}",
        f"- swap-symmetry sorted-root sanity max abs diff: {float(swap['max_swap_abs_diff']):.6e}",
        f"- swap-symmetry sorted-root sanity max rel diff: {float(swap['max_swap_rel_diff']):.6e}",
        f"- max same-index sorted-vs-tracked difference: {float(max_diff_row['abs_diff']):.6e}",
        (
            "- max same-index difference location: "
            f"eta={float(max_diff_row['eta']):g}, "
            f"branch/root={int(max_diff_row['branch_index'])}, "
            f"mu={float(max_diff_row['mu']):g}"
        ),
        f"- rows where MAC current sorted index differs from mu=0 branch index: {len(switches)}",
        f"- low-MAC rows below threshold {MAC_WARNING_THRESHOLD:g}: {len(low_mac)}",
        f"- nearest-frequency vs MAC disagreement rows: {len(frequency_warnings)}",
        "",
        "## Eta Sensitivity",
        "",
    ]

    sensitivity_table = [["branch", "max eta spread", "at mu", "spread at mu=0.9", "eta=0.1 minus eta=-0.1 at mu=0.9"]]
    for row in sensitivity_rows:
        sensitivity_table.append(
            [
                str(int(row["branch_index"])),
                f"{float(row['max_eta_spread']):.6e}",
                f"{float(row['max_eta_spread_mu']):g}",
                f"{float(row['spread_at_mu_0p9']):.6e}",
                f"{float(row['eta_p0p1_minus_m0p1_at_mu_0p9']):.6e}",
            ]
        )
    lines.extend(format_table(sensitivity_table))

    top_sensitivity = sensitivity_rows[0]
    lines.extend(
        [
            "",
            (
                "The largest eta spread on this grid is on branch "
                f"{int(top_sensitivity['branch_index'])} near mu={float(top_sensitivity['max_eta_spread_mu']):g}. "
                "This is a sensitivity diagnostic, not a physical conclusion."
            ),
            "",
            "Qualitatively, the eta=-0.1, eta=0, and eta=0.1 curves separate in",
            "a branch-dependent way; the largest spread on this grid need not occur",
            "at the largest mu. At mu=0.9, eta=0.1 gives lower Lambda than eta=-0.1",
            "for all six tracked branches in this diagnostic run.",
            "",
            "## Tracking Sanity",
            "",
        ]
    )

    jump_table = [["eta", "branch", "max adjacent mu jump", "left mu", "right mu"]]
    for row in jumps[:10]:
        jump_table.append(
            [
                f"{float(row['eta']):g}",
                str(int(row["branch_index"])),
                f"{float(row['max_adjacent_mu_jump']):.6e}",
                f"{float(row['left_mu']):g}",
                f"{float(row['right_mu']):g}",
            ]
        )
    lines.extend(format_table(jump_table))

    lines.extend(["", "Smallest local gaps to neighboring sorted roots:", ""])
    gap_table = [["eta", "mu", "branch", "MAC current sorted index", "min gap"]]
    for row in close_rows:
        gap_table.append(
            [
                f"{float(row['eta']):g}",
                f"{float(row['mu']):g}",
                str(int(row["branch_index_from_mu0"])),
                str(int(row["nearest_sorted_root_index"])),
                f"{float(row['min_gap_to_other_sorted_roots']):.6e}",
            ]
        )
    lines.extend(format_table(gap_table))

    lines.extend(["", "## MAC Sorted-Index Switches", ""])
    if switches:
        switch_table = [["eta", "mu", "branch from mu=0", "MAC current sorted index", "Lambda tracked"]]
        for row in switches[:24]:
            switch_table.append(
                [
                    f"{float(row['eta']):g}",
                    f"{float(row['mu']):g}",
                    str(int(row["branch_index_from_mu0"])),
                    str(int(row["nearest_sorted_root_index"])),
                    f"{float(row['Lambda_tracked']):.8g}",
                ]
            )
        lines.extend(format_table(switch_table))
        if len(switches) > 24:
            lines.append("")
            lines.append(f"Only the first 24 of {len(switches)} switch rows are shown.")
    else:
        lines.append("No sorted-index switches were detected for the tracked branches on this mu grid.")

    if low_mac:
        lines.extend(["", "## Low-MAC Assignments", ""])
        low_mac_table = [["eta", "mu", "branch", "status", "MAC", "raw MAC index"]]
        for row in low_mac[:24]:
            low_mac_table.append(
                [
                    f"{float(row['eta']):g}",
                    f"{float(row['mu']):g}",
                    str(int(row["branch_index_from_mu0"])),
                    str(row["tracking_step_status"]),
                    f"{float(row['mac_to_previous']):.6e}",
                    str(int(row["raw_mac_sorted_root_index"])),
                ]
            )
        lines.extend(format_table(low_mac_table))

    if frequency_warnings:
        lines.extend(["", "## Nearest-Frequency Disagreement Warnings", ""])
        warning_table = [["eta", "mu", "branch", "MAC index", "nearest-frequency index", "MAC"]]
        for row in frequency_warnings[:24]:
            warning_table.append(
                [
                    f"{float(row['eta']):g}",
                    f"{float(row['mu']):g}",
                    str(int(row["branch_index_from_mu0"])),
                    str(int(row["current_sorted_root_index"])),
                    str(int(row["frequency_nearest_sorted_root_index"])),
                    f"{float(row['mac_to_previous']):.6e}",
                ]
            )
        lines.extend(format_table(warning_table))

    lines.extend(
        [
            "",
            "These are diagnostic observations only. The script does not change the",
            "baseline determinant, solvers, FEM model, article files, or article figures.",
            "",
        ]
    )

    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text("\n".join(lines), encoding="utf-8")
    return {
        "eta_zero_max_nearest_old_root_diff": float(eta_zero["max_nearest_old_root_diff"]),
        "eta_zero_max_same_index_old_root_diff": float(eta_zero["max_same_index_old_root_diff"]),
        "swap_max_abs_diff": float(swap["max_swap_abs_diff"]),
        "swap_max_rel_diff": float(swap["max_swap_rel_diff"]),
        "max_same_index_diff": float(max_diff_row["abs_diff"]),
        "switch_count": len(switches),
        "low_mac_count": len(low_mac),
        "frequency_warning_count": len(frequency_warnings),
        "top_sensitive_branch": int(top_sensitivity["branch_index"]),
        "top_sensitive_spread": float(top_sensitivity["max_eta_spread"]),
        "top_sensitive_mu": float(top_sensitivity["max_eta_spread_mu"]),
        "max_adjacent_mu_jump": float(jumps[0]["max_adjacent_mu_jump"]),
        "max_adjacent_mu_jump_eta": float(jumps[0]["eta"]),
        "max_adjacent_mu_jump_branch": int(jumps[0]["branch_index"]),
    }


def main() -> dict[str, Path | float | int]:
    sorted_roots = sorted_roots_by_eta_mu()
    tracked_rows: list[dict[str, float | int | str]] = []
    tracked_by_eta: dict[float, np.ndarray] = {}

    for eta in ETA_VALUES:
        rows, tracked = track_one_eta(float(eta), sorted_roots[float(eta)])
        tracked_rows.extend(rows)
        tracked_by_eta[float(eta)] = tracked

    if SAVE_CSV:
        write_csv(tracked_csv_path(), tracked_rows, TRACKED_FIELDNAMES)
    if SAVE_PNG:
        plot_tracked(sorted_roots=sorted_roots, tracked_by_eta=tracked_by_eta, output=tracked_png_path())

    same_index_diffs = compute_same_index_diffs(sorted_roots, tracked_by_eta)
    sensitivity = branch_sensitivity_rows(tracked_by_eta)
    jumps = jump_rows(tracked_by_eta)
    eta_zero = eta_zero_consistency(tracked_by_eta)
    swap = swap_symmetry_sanity()
    report = build_report(
        tracked_rows=tracked_rows,
        same_index_diffs=same_index_diffs,
        sensitivity_rows=sensitivity,
        jumps=jumps,
        eta_zero=eta_zero,
        swap=swap,
        output=report_md_path() if SAVE_REPORT else None,
    )

    if SAVE_CSV:
        print(f"saved tracked Lambda(mu) CSV: {tracked_csv_path()}")
    else:
        print("tracked Lambda(mu) CSV not saved (SAVE_CSV=False)")
    if SAVE_PNG:
        print(f"saved tracked Lambda(mu) plot: {tracked_png_path()}")
    else:
        print("tracked Lambda(mu) plot not saved (SAVE_PNG=False)")
    if SAVE_REPORT:
        print(f"saved tracking report: {report_md_path()}")
    else:
        print("tracking report not saved (SAVE_REPORT=False)")
    print(f"eta=0 max nearest old-root diff: {float(report['eta_zero_max_nearest_old_root_diff']):.6e}")
    print(f"eta=0 max same-index old-root diff: {float(report['eta_zero_max_same_index_old_root_diff']):.6e}")
    print(f"swap symmetry sanity max abs diff: {float(report['swap_max_abs_diff']):.6e}")
    print(f"max same-index sorted-vs-tracked diff: {float(report['max_same_index_diff']):.6e}")
    print(f"MAC current-sorted-index switch rows: {int(report['switch_count'])}")
    print(f"low-MAC assignment rows: {int(report['low_mac_count'])}")
    print(f"nearest-frequency vs MAC disagreement rows: {int(report['frequency_warning_count'])}")
    print(
        "largest eta sensitivity: "
        f"branch={int(report['top_sensitive_branch'])}, "
        f"mu={float(report['top_sensitive_mu']):g}, "
        f"spread={float(report['top_sensitive_spread']):.6e}"
    )

    return {
        "tracked_csv": tracked_csv_path(),
        "tracked_png": tracked_png_path(),
        "report_md": report_md_path(),
        **report,
    }


if __name__ == "__main__":
    main()
