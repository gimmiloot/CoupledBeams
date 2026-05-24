from __future__ import annotations

import csv
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

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    find_first_n_roots_eta,
    thickness_to_length_ratios,
    thin_rod_validity,
)
from scripts.lib.thickness_mismatch_mac_tracking import track_mu_branches_shape_mac  # noqa: E402


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.005

NUM_ROOTS = 8
NUM_DESCENDANTS = 8
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
THICKNESS_RATIO_LIMIT = 0.1

ZOOM_MU_MIN = 0.25
ZOOM_MU_MAX = 0.9
ZOOM_MODE_MIN = 4
ZOOM_MODE_MAX = 8

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.png"
OUTPUT_ZOOM_PNG = OUTPUT_DIR / "thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes_zoom_4_8.png"
OUTPUT_CSV = OUTPUT_DIR / "thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes.csv"
OUTPUT_REPORT = OUTPUT_DIR / "thickness_mismatch_eta_p0p5_beta15_global_spectrum_8modes_report.md"

MU_VALUES = np.round(np.arange(MU_MIN, MU_MAX + 0.5 * MU_STEP, MU_STEP), 10)

CSV_FIELDNAMES = [
    "mu",
    "sorted_mode_index",
    "Lambda_sorted",
    "descendant_branch",
    "Lambda_descendant",
    "current_sorted_position",
    "MAC_to_previous",
    "accepted_MAC_to_previous",
    "diagnostic_candidate_position",
    "diagnostic_candidate_Lambda",
    "diagnostic_candidate_MAC",
    "diagnostic_candidate_position_jump",
    "raw_mac_sorted_position",
    "low_mac",
    "low_margin",
    "sorted_position_jump",
    "raw_sorted_position_jump",
    "unresolved_assignment",
    "blocked_by_unresolved_neighbor",
    "requires_refined_check",
    "thickness_ratio_1",
    "thickness_ratio_2",
    "thin_rod_valid",
    "suspicious_assignment",
    "mac_margin",
    "frequency_nearest_sorted_position",
    "frequency_mac_disagreement",
    "tracking_step_status",
]


def roots_for(mu: float) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        ETA,
        NUM_ROOTS,
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots)):
        raise RuntimeError(f"Missing roots for mu={mu:g}, eta={ETA:g}.")
    return roots


def compute_roots_by_mu() -> dict[float, np.ndarray]:
    return {float(mu): roots_for(float(mu)) for mu in MU_VALUES}


def ratios_for(mu: float) -> tuple[float, float, bool]:
    ratio1, ratio2 = thickness_to_length_ratios(EPSILON, float(mu), ETA)
    valid = thin_rod_validity(EPSILON, float(mu), ETA, limit=THICKNESS_RATIO_LIMIT)
    return ratio1, ratio2, bool(valid)


def ratio_rows() -> dict[float, tuple[float, float, bool]]:
    return {float(mu): ratios_for(float(mu)) for mu in MU_VALUES}


def tracking_rows_by_mu_branch(
    rows: Sequence[dict[str, float | int | str]],
) -> dict[tuple[float, int], dict[str, float | int | str]]:
    return {
        (round(float(row["mu"]), 10), int(row["branch_index_from_mu0"])): row
        for row in rows
    }


def occupants_by_mu_sorted_position(
    rows: Sequence[dict[str, float | int | str]],
) -> dict[tuple[float, int], list[dict[str, float | int | str]]]:
    occupants: dict[tuple[float, int], list[dict[str, float | int | str]]] = {}
    for row in rows:
        key = (round(float(row["mu"]), 10), int(row["mac_sorted_root_index"]))
        occupants.setdefault(key, []).append(row)
    return occupants


def write_csv(
    *,
    path: Path,
    roots_by_mu: dict[float, np.ndarray],
    tracking_rows: Sequence[dict[str, float | int | str]],
    ratios_by_mu: dict[float, tuple[float, float, bool]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    occupants = occupants_by_mu_sorted_position(tracking_rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        for mu in MU_VALUES:
            mu_f = float(mu)
            ratio1, ratio2, valid = ratios_by_mu[mu_f]
            roots = roots_by_mu[mu_f]
            for sorted_idx in range(1, NUM_ROOTS + 1):
                local_occupants = occupants.get((round(mu_f, 10), sorted_idx), [])
                if local_occupants:
                    for row in local_occupants:
                        writer.writerow(
                            {
                                "mu": mu_f,
                                "sorted_mode_index": sorted_idx,
                                "Lambda_sorted": float(roots[sorted_idx - 1]),
                                "descendant_branch": int(row["branch_index_from_mu0"]),
                                "Lambda_descendant": float(row["Lambda_tracked"]),
                                "current_sorted_position": int(row["mac_sorted_root_index"]),
                                "MAC_to_previous": float(row["mac_to_previous"]),
                                "accepted_MAC_to_previous": float(row["accepted_mac_to_previous"]),
                                "diagnostic_candidate_position": int(row["diagnostic_candidate_sorted_position"]),
                                "diagnostic_candidate_Lambda": float(row["diagnostic_candidate_Lambda"]),
                                "diagnostic_candidate_MAC": float(row["diagnostic_candidate_mac_to_previous"]),
                                "diagnostic_candidate_position_jump": int(row["diagnostic_candidate_sorted_position_jump"]),
                                "raw_mac_sorted_position": int(row["raw_mac_sorted_root_index"]),
                                "low_mac": str(row["low_mac"]),
                                "low_margin": str(row["low_margin"]),
                                "sorted_position_jump": int(row["sorted_position_jump"]),
                                "raw_sorted_position_jump": int(row["raw_sorted_position_jump"]),
                                "unresolved_assignment": str(row["unresolved_assignment"]),
                                "blocked_by_unresolved_neighbor": str(row["blocked_by_unresolved_neighbor"]),
                                "requires_refined_check": str(row["requires_refined_check"]),
                                "thickness_ratio_1": ratio1,
                                "thickness_ratio_2": ratio2,
                                "thin_rod_valid": "yes" if valid else "no",
                                "suspicious_assignment": str(row["suspicious_assignment"]),
                                "mac_margin": float(row["mac_margin"]),
                                "frequency_nearest_sorted_position": int(row["nearest_sorted_root_index"]),
                                "frequency_mac_disagreement": str(row["frequency_mac_disagreement"]),
                                "tracking_step_status": str(row["tracking_step_status"]),
                            }
                        )
                else:
                    writer.writerow(
                        {
                            "mu": mu_f,
                            "sorted_mode_index": sorted_idx,
                            "Lambda_sorted": float(roots[sorted_idx - 1]),
                            "descendant_branch": "",
                            "Lambda_descendant": "",
                            "current_sorted_position": "",
                            "MAC_to_previous": "",
                            "accepted_MAC_to_previous": "",
                            "diagnostic_candidate_position": "",
                            "diagnostic_candidate_Lambda": "",
                            "diagnostic_candidate_MAC": "",
                            "diagnostic_candidate_position_jump": "",
                            "raw_mac_sorted_position": "",
                            "low_mac": "",
                            "low_margin": "",
                            "sorted_position_jump": "",
                            "raw_sorted_position_jump": "",
                            "unresolved_assignment": "",
                            "blocked_by_unresolved_neighbor": "",
                            "requires_refined_check": "",
                            "thickness_ratio_1": ratio1,
                            "thickness_ratio_2": ratio2,
                            "thin_rod_valid": "yes" if valid else "no",
                            "suspicious_assignment": "",
                            "mac_margin": "",
                            "frequency_nearest_sorted_position": "",
                            "frequency_mac_disagreement": "",
                            "tracking_step_status": "",
                        }
                    )


def contiguous_true_runs(mask: np.ndarray) -> list[tuple[int, int]]:
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(mask):
        if bool(value) and start is None:
            start = idx
        elif not bool(value) and start is not None:
            runs.append((start, idx - 1))
            start = None
    if start is not None:
        runs.append((start, len(mask) - 1))
    return runs


def plot_segmented(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    valid: np.ndarray,
    *,
    color: str,
    lw: float,
    label: str | None,
    alpha: float = 1.0,
) -> None:
    for idx, mask_value in enumerate((True, False)):
        mask = valid == mask_value
        linestyle = "-" if mask_value else "--"
        local_label = label if idx == 0 else None
        for start, end in contiguous_true_runs(mask):
            sl = slice(start, end + 1)
            ax.plot(x[sl], y[sl], color=color, lw=lw, ls=linestyle, alpha=alpha, label=local_label)
            local_label = None


def plot_overview(
    *,
    output: Path,
    roots_by_mu: dict[float, np.ndarray],
    tracking_rows: Sequence[dict[str, float | int | str]],
    ratios_by_mu: dict[float, tuple[float, float, bool]],
    mu_min: float,
    mu_max: float,
    mode_min: int,
    mode_max: int,
    title_suffix: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    mu = MU_VALUES[(MU_VALUES >= mu_min - 1e-12) & (MU_VALUES <= mu_max + 1e-12)]
    valid = np.array([bool(ratios_by_mu[float(value)][2]) for value in mu], dtype=bool)
    tracking_by_branch = tracking_rows_by_mu_branch(tracking_rows)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    modes = tuple(range(int(mode_min), int(mode_max) + 1))

    fig, axes = plt.subplots(3, 1, figsize=(11.0, 9.2), sharex=True, height_ratios=[1.35, 1.65, 1.0])
    ax_sorted, ax_desc, ax_pos = axes

    for mode in modes:
        y_sorted = np.array([float(roots_by_mu[float(value)][mode - 1]) for value in mu], dtype=float)
        plot_segmented(
            ax_sorted,
            mu,
            y_sorted,
            valid,
            color="0.35",
            lw=1.0,
            alpha=0.72,
            label=f"sorted {mode}" if len(modes) <= 5 else None,
        )
        plot_segmented(
            ax_desc,
            mu,
            y_sorted,
            valid,
            color="0.72",
            lw=0.8,
            alpha=0.55,
            label="sorted roots" if mode == modes[0] else None,
        )

    suspicious_mu: list[float] = []
    suspicious_lambda: list[float] = []
    suspicious_pos_mu: list[float] = []
    suspicious_pos_y: list[int] = []
    for branch in modes:
        y_desc: list[float] = []
        y_pos: list[int] = []
        for mu_value in mu:
            row = tracking_by_branch[(round(float(mu_value), 10), branch)]
            y_desc.append(float(row["Lambda_tracked"]))
            y_pos.append(int(row["mac_sorted_root_index"]))
            if str(row["requires_refined_check"]) == "yes":
                suspicious_mu.append(float(mu_value))
                suspicious_lambda.append(float(row["Lambda_tracked"]))
                suspicious_pos_mu.append(float(mu_value))
                suspicious_pos_y.append(int(row["mac_sorted_root_index"]))

        color = colors[(branch - 1) % len(colors)]
        plot_segmented(
            ax_desc,
            mu,
            np.asarray(y_desc, dtype=float),
            valid,
            color=color,
            lw=1.9,
            label=f"desc {branch}",
        )
        ax_pos.step(mu, y_pos, where="post", color=color, lw=1.5, label=f"desc {branch}")

    if suspicious_mu:
        ax_desc.scatter(
            suspicious_mu,
            suspicious_lambda,
            marker="x",
            s=28,
            color="black",
            linewidths=0.9,
            label="requires refined check",
            zorder=5,
        )
        ax_pos.scatter(
            suspicious_pos_mu,
            suspicious_pos_y,
            marker="x",
            s=24,
            color="black",
            linewidths=0.9,
            zorder=5,
        )

    ax_sorted.set_ylabel(r"sorted $\Lambda$")
    ax_sorted.set_title(
        rf"eta={ETA:g}, beta={BETA_DEG:g}$^\circ$, epsilon={EPSILON:g}: sorted spectrum {title_suffix}"
    )
    ax_sorted.grid(True, color="0.88", linewidth=0.6)

    ax_desc.set_ylabel(r"descendant $\Lambda$")
    ax_desc.set_title("Descendant branches seeded at mu=0; gray lines are sorted roots")
    ax_desc.grid(True, color="0.88", linewidth=0.6)
    ax_desc.legend(loc="best", ncol=4, fontsize=8, frameon=False)

    ax_pos.set_xlabel(r"$\mu$")
    ax_pos.set_ylabel("accepted sorted position")
    ax_pos.set_yticks(range(max(1, mode_min - 1), mode_max + 2))
    ax_pos.set_ylim(mode_min - 1.25, mode_max + 1.25)
    ax_pos.grid(True, color="0.88", linewidth=0.6)
    ax_pos.legend(loc="best", ncol=4, fontsize=8, frameon=False)

    fig.text(
        0.5,
        0.01,
        "Solid/dashed segments show 2r_i/l_i <= 0.1; accepted positions retain previous value when candidate tracking is unresolved.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 1.0))
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def sorted_gap_rows(roots_by_mu: dict[float, np.ndarray], top_n: int = 12) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for mu in MU_VALUES:
        roots = roots_by_mu[float(mu)]
        for mode in range(1, NUM_ROOTS):
            gap = float(roots[mode] - roots[mode - 1])
            rows.append(
                {
                    "mu": float(mu),
                    "mode_left": int(mode),
                    "mode_right": int(mode) + 1,
                    "gap": gap,
                    "Lambda_left": float(roots[mode - 1]),
                    "Lambda_right": float(roots[mode]),
                }
            )
    rows.sort(key=lambda row: float(row["gap"]))
    return rows[:top_n]


def grouped_mu_ranges(rows: Sequence[dict[str, float | int | str]]) -> dict[int, tuple[float, float, int]]:
    grouped: dict[int, list[float]] = {}
    for row in rows:
        grouped.setdefault(int(row["branch_index_from_mu0"]), []).append(float(row["mu"]))
    return {
        branch: (float(min(values)), float(max(values)), len(values))
        for branch, values in sorted(grouped.items())
    }


def first_nonself_rows(
    tracking_rows: Sequence[dict[str, float | int | str]],
) -> dict[int, dict[str, float | int | str] | None]:
    out: dict[int, dict[str, float | int | str] | None] = {}
    for branch in range(1, NUM_DESCENDANTS + 1):
        branch_rows = sorted(
            (row for row in tracking_rows if int(row["branch_index_from_mu0"]) == branch),
            key=lambda row: float(row["mu"]),
        )
        out[branch] = next(
            (row for row in branch_rows if int(row["mac_sorted_root_index"]) != branch),
            None,
        )
    return out


def positions_by_descendant_after(
    tracking_rows: Sequence[dict[str, float | int | str]],
    *,
    mu_threshold: float,
) -> dict[int, list[int]]:
    return {
        branch: sorted(
            {
                int(row["mac_sorted_root_index"])
                for row in tracking_rows
                if int(row["branch_index_from_mu0"]) == branch and float(row["mu"]) >= mu_threshold - 1e-12
            }
        )
        for branch in range(1, NUM_DESCENDANTS + 1)
    }


def occupants_by_position_after(
    tracking_rows: Sequence[dict[str, float | int | str]],
    *,
    mu_threshold: float,
) -> dict[int, list[int]]:
    out: dict[int, set[int]] = {position: set() for position in range(1, NUM_ROOTS + 1)}
    for row in tracking_rows:
        if float(row["mu"]) < mu_threshold - 1e-12:
            continue
        position = int(row["mac_sorted_root_index"])
        if 1 <= position <= NUM_ROOTS:
            out[position].add(int(row["branch_index_from_mu0"]))
    return {position: sorted(values) for position, values in out.items()}


def missing_duplicate_rows(
    tracking_rows: Sequence[dict[str, float | int | str]],
) -> list[dict[str, object]]:
    by_mu: dict[float, list[int]] = {}
    for row in tracking_rows:
        by_mu.setdefault(float(row["mu"]), []).append(int(row["mac_sorted_root_index"]))
    rows: list[dict[str, object]] = []
    expected = set(range(1, NUM_ROOTS + 1))
    for mu, positions in sorted(by_mu.items()):
        present = set(positions)
        duplicates = sorted({position for position in positions if positions.count(position) > 1})
        missing = sorted(expected - present)
        if duplicates or missing:
            rows.append({"mu": mu, "missing": missing, "duplicates": duplicates})
    return rows


def validity_segments(ratios_by_mu: dict[float, tuple[float, float, bool]], *, valid: bool) -> list[tuple[float, float]]:
    mask = np.array([bool(ratios_by_mu[float(mu)][2]) == bool(valid) for mu in MU_VALUES], dtype=bool)
    return [(float(MU_VALUES[start]), float(MU_VALUES[end])) for start, end in contiguous_true_runs(mask)]


def format_table(rows: Sequence[Sequence[str]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(row[col]) for row in rows) for col in range(len(rows[0]))]
    lines = []
    for idx, row in enumerate(rows):
        lines.append("| " + " | ".join(value.ljust(widths[col]) for col, value in enumerate(row)) + " |")
        if idx == 0:
            lines.append("| " + " | ".join("-" * widths[col] for col in range(len(row))) + " |")
    return lines


def write_report(
    *,
    path: Path,
    roots_by_mu: dict[float, np.ndarray],
    tracking_rows: Sequence[dict[str, float | int | str]],
    ratios_by_mu: dict[float, tuple[float, float, bool]],
) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    low_mac = [row for row in tracking_rows if str(row["low_mac"]) == "yes"]
    low_margin = [row for row in tracking_rows if str(row["low_margin"]) == "yes"]
    unresolved = [row for row in tracking_rows if str(row["unresolved_assignment"]) == "yes"]
    blocked = [row for row in tracking_rows if str(row["blocked_by_unresolved_neighbor"]) == "yes"]
    suspicious = [row for row in tracking_rows if str(row["suspicious_assignment"]) == "yes"]
    refined = [row for row in tracking_rows if str(row["requires_refined_check"]) == "yes"]
    jumps = [row for row in tracking_rows if abs(int(row["sorted_position_jump"])) > MAX_SORTED_POSITION_JUMP]
    disagreements = [row for row in tracking_rows if str(row["frequency_mac_disagreement"]) == "yes"]
    missing_duplicate = missing_duplicate_rows(tracking_rows)
    positions_after = positions_by_descendant_after(tracking_rows, mu_threshold=0.7)
    occupants_after = occupants_by_position_after(tracking_rows, mu_threshold=0.7)
    first_nonself = first_nonself_rows(tracking_rows)
    gap_rows = sorted_gap_rows(roots_by_mu)
    invalid_segments = validity_segments(ratios_by_mu, valid=False)

    ratios = np.array([[ratios_by_mu[float(mu)][0], ratios_by_mu[float(mu)][1]] for mu in MU_VALUES], dtype=float)
    max_flat = int(np.argmax(ratios))
    max_mu_idx, max_rod_idx = np.unravel_index(max_flat, ratios.shape)
    max_ratio = float(ratios[max_mu_idx, max_rod_idx])
    max_ratio_mu = float(MU_VALUES[max_mu_idx])
    max_ratio_rod = int(max_rod_idx) + 1

    lines = [
        "# Eta=0.5 Global Spectrum Diagnostic",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu range: {MU_MIN:g} .. {MU_MAX:g}",
        f"- mu step: {MU_STEP:g}",
        f"- sorted roots: first {NUM_ROOTS}",
        f"- descendant branches: first {NUM_DESCENDANTS}",
        f"- MAC threshold: {MAC_WARNING_THRESHOLD:g}",
        f"- MAC margin threshold: {MAC_MARGIN_WARNING_THRESHOLD:g}",
        "",
        "## What The Figures Show",
        "",
        f"- overview PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
        f"- zoom PNG: `{OUTPUT_ZOOM_PNG.relative_to(REPO_ROOT)}`",
        "- top panel: sorted spectrum roots",
        "- middle panel: descendant branches seeded at `mu=0`; gray lines repeat sorted roots",
        "- bottom panel: accepted canonical sorted position of each descendant branch",
        "- black x markers: unresolved, low-MAC, low-margin, jump, or other assignments requiring refined review",
        "",
        "This is an overview diagnostic only. It does not resolve branch identity",
        "by itself and does not rewrite the large-eta plot.",
        "",
        "## Thin-Rod Criterion",
        "",
        f"- criterion: `2*r_i/l_i <= {THICKNESS_RATIO_LIMIT:g}` for both rods",
        f"- max ratio on grid: {max_ratio:.6g} at mu={max_ratio_mu:g}, rod {max_ratio_rod}",
    ]
    if invalid_segments:
        lines.append("- WARNING: criterion violations were found:")
        for start, end in invalid_segments:
            lines.append(f"  - mu={start:g}..{end:g}")
    else:
        lines.append("- No violations of `2*r_i/l_i <= 0.1` were found on the grid.")

    lines.extend(["", "## Closest Sorted-Spectrum Gaps", ""])
    gap_table = [["rank", "mu", "sorted modes", "gap", "Lambda left", "Lambda right"]]
    for idx, row in enumerate(gap_rows[:10], start=1):
        gap_table.append(
            [
                str(idx),
                f"{float(row['mu']):.3f}",
                f"{int(row['mode_left'])}-{int(row['mode_right'])}",
                f"{float(row['gap']):.6g}",
                f"{float(row['Lambda_left']):.6g}",
                f"{float(row['Lambda_right']):.6g}",
            ]
        )
    lines.extend(format_table(gap_table))

    lines.extend(["", "## Sorted Positions After mu about 0.7", ""])
    desc_table = [["descendant", "accepted sorted positions for mu>=0.7"]]
    for branch in range(1, NUM_DESCENDANTS + 1):
        desc_table.append([str(branch), ", ".join(str(value) for value in positions_after[branch])])
    lines.extend(format_table(desc_table))
    lines.append("")
    occupant_table = [["sorted position", "occupying descendants for mu>=0.7"]]
    for position in range(1, NUM_ROOTS + 1):
        values = occupants_after[position]
        occupant_table.append([str(position), ", ".join(str(value) for value in values) if values else "none"])
    lines.extend(format_table(occupant_table))

    lines.extend(["", "## Tracking Warnings", ""])
    lines.extend(
        [
            f"- low-MAC rows: {len(low_mac)}",
            f"- low-margin rows: {len(low_margin)}",
            f"- unresolved rows: {len(unresolved)}",
            f"- rows blocked by an unresolved neighbor: {len(blocked)}",
            f"- suspicious rows: {len(suspicious)}",
            f"- rows requiring refined check: {len(refined)}",
            f"- sorted-position jumps larger than {MAX_SORTED_POSITION_JUMP}: {len(jumps)}",
            f"- nearest-frequency vs MAC disagreement rows: {len(disagreements)}",
        ]
    )

    grouped_refined = grouped_mu_ranges(refined)
    if grouped_refined:
        lines.extend(["", "Rows requiring refined check by descendant:", ""])
        refined_table = [["descendant", "mu min", "mu max", "rows"]]
        for branch, (start, end, count) in grouped_refined.items():
            refined_table.append([str(branch), f"{start:.3f}", f"{end:.3f}", str(count)])
        lines.extend(format_table(refined_table))
    else:
        lines.append("")
        lines.append("No rows require refined check on this grid.")

    lines.extend(["", "First non-seed sorted position by descendant:", ""])
    first_table = [["descendant", "first accepted non-seed position", "mu", "candidate position", "candidate MAC", "status"]]
    for branch in range(1, NUM_DESCENDANTS + 1):
        row = first_nonself[branch]
        if row is None:
            first_table.append([str(branch), "none", "-", "-", "-", "-"])
        else:
            first_table.append(
                [
                    str(branch),
                    str(int(row["mac_sorted_root_index"])),
                    f"{float(row['mu']):.3f}",
                    str(int(row["diagnostic_candidate_sorted_position"])),
                    f"{float(row['diagnostic_candidate_mac_to_previous']):.6g}",
                    str(row["tracking_step_status"]),
                ]
            )
    lines.extend(format_table(first_table))

    lines.extend(["", "## Around mu about 0.38", ""])
    local_rows = [
        row
        for row in tracking_rows
        if 0.35 <= float(row["mu"]) <= 0.41 and str(row["requires_refined_check"]) == "yes"
    ]
    if local_rows:
        local_table = [["mu", "descendant", "accepted position", "candidate position", "candidate MAC", "flags"]]
        for row in local_rows[:16]:
            flags = ",".join(
                name
                for name in (
                    "low_mac",
                    "low_margin",
                    "blocked_by_unresolved_neighbor",
                    "unresolved_assignment",
                    "requires_refined_check",
                )
                if str(row[name]) == "yes"
            )
            local_table.append(
                [
                    f"{float(row['mu']):.3f}",
                    str(int(row["branch_index_from_mu0"])),
                    str(int(row["mac_sorted_root_index"])),
                    str(int(row["diagnostic_candidate_sorted_position"])),
                    f"{float(row['diagnostic_candidate_mac_to_previous']):.6g}",
                    flags,
                ]
            )
        lines.extend(format_table(local_table))
    else:
        lines.append("No rows requiring refined check were found in `0.35 <= mu <= 0.41`.")

    lines.extend(["", "## Missing Or Duplicate Sorted Positions", ""])
    if missing_duplicate:
        lines.append(
            f"Found {len(missing_duplicate)} mu points with missing or duplicate sorted positions among descendants 1..{NUM_DESCENDANTS}."
        )
        miss_table = [["mu", "missing", "duplicates"]]
        for row in missing_duplicate[:20]:
            miss_table.append(
                [
                    f"{float(row['mu']):.3f}",
                    ",".join(str(value) for value in row["missing"]) or "none",
                    ",".join(str(value) for value in row["duplicates"]) or "none",
                ]
            )
        lines.extend(format_table(miss_table))
    else:
        lines.append(f"No missing or duplicate sorted positions were found among descendants 1..{NUM_DESCENDANTS}.")

    lines.extend(["", "## Overview Interpretation", ""])
    lines.extend(
        [
            "The closest sorted-spectrum gaps concentrate in the lower and middle",
            "part of the spectrum, especially the pairs listed above. For the",
            "branch-identity question, the key diagnostic is that accepted",
            "canonical sorted positions are not changed by low-MAC or",
            "low-margin candidate assignments.",
            "",
            f"- descendant 5 accepted positions for mu>=0.7: {', '.join(str(value) for value in positions_after[5])}",
            f"- descendant 6 accepted positions for mu>=0.7: {', '.join(str(value) for value in positions_after[6])}",
            f"- descendant 7 accepted positions for mu>=0.7: {', '.join(str(value) for value in positions_after[7])}",
        ]
    )
    if positions_after[5] == [5] and positions_after[6] == [6] and positions_after[7] == [7]:
        lines.append("- On this overview grid, descendants 5, 6, and 7 stay at sorted positions 5, 6, and 7 after mu about 0.7.")
    else:
        lines.append(
            "- On this overview grid, descendants 5, 6, and 7 do not all stay at sorted positions 5, 6, and 7 after mu about 0.7."
        )
    first_branch5 = first_nonself[5]
    if first_branch5 is not None:
        lines.append(
            "- The first accepted non-seed sorted position for descendant 5 occurs at "
            f"mu={float(first_branch5['mu']):g}, accepted position {int(first_branch5['mac_sorted_root_index'])}, "
            f"candidate position {int(first_branch5['diagnostic_candidate_sorted_position'])}, "
            f"candidate MAC={float(first_branch5['diagnostic_candidate_mac_to_previous']):.6g}. "
            "This is a warning event, not a physical conclusion."
        )
    else:
        local_branch5 = [
            row
            for row in tracking_rows
            if (
                int(row["branch_index_from_mu0"]) == 5
                and str(row["unresolved_assignment"]) == "yes"
                and int(row["diagnostic_candidate_sorted_position"]) != int(row["mac_sorted_root_index"])
            )
        ]
        if local_branch5:
            first_unresolved = local_branch5[0]
            lines.append(
                "- Descendant 5 has no accepted non-seed sorted-position change. "
                f"The raw candidate first suggests position {int(first_unresolved['diagnostic_candidate_sorted_position'])} "
                f"at mu={float(first_unresolved['mu']):g}, but this step is unresolved "
                f"(candidate MAC={float(first_unresolved['diagnostic_candidate_mac_to_previous']):.6g}), "
                "so the previous canonical position is retained."
            )

    lines.extend(
        [
            "",
            "The article files, article figures, old determinant, old solvers,",
            "and FEM physical model are not modified by this diagnostic.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")

    return {
        "low_mac": len(low_mac),
        "low_margin": len(low_margin),
        "unresolved": len(unresolved),
        "blocked": len(blocked),
        "suspicious": len(suspicious),
        "requires_refined": len(refined),
        "jumps": len(jumps),
        "disagreements": len(disagreements),
        "missing_duplicate": len(missing_duplicate),
        "invalid_segments": len(invalid_segments),
        "max_ratio": max_ratio,
        "max_ratio_mu": max_ratio_mu,
        "max_ratio_rod": max_ratio_rod,
        "positions_after": positions_after,
        "occupants_after": occupants_after,
    }


def main() -> dict[str, object]:
    if not np.isclose(float(MU_VALUES[0]), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Descendant branches are seeded at mu=0; set MU_MIN=0.0.")

    roots_by_mu = compute_roots_by_mu()
    ratios_by_mu = ratio_rows()
    tracking = track_mu_branches_shape_mac(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        roots_by_mu=roots_by_mu,
        num_tracked_branches=NUM_DESCENDANTS,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
        max_mu_step_for_confidence=MU_STEP,
    )

    write_csv(path=OUTPUT_CSV, roots_by_mu=roots_by_mu, tracking_rows=tracking.rows, ratios_by_mu=ratios_by_mu)
    plot_overview(
        output=OUTPUT_PNG,
        roots_by_mu=roots_by_mu,
        tracking_rows=tracking.rows,
        ratios_by_mu=ratios_by_mu,
        mu_min=MU_MIN,
        mu_max=MU_MAX,
        mode_min=1,
        mode_max=NUM_ROOTS,
        title_suffix="1..8",
    )
    plot_overview(
        output=OUTPUT_ZOOM_PNG,
        roots_by_mu=roots_by_mu,
        tracking_rows=tracking.rows,
        ratios_by_mu=ratios_by_mu,
        mu_min=ZOOM_MU_MIN,
        mu_max=ZOOM_MU_MAX,
        mode_min=ZOOM_MODE_MIN,
        mode_max=ZOOM_MODE_MAX,
        title_suffix="zoom 4..8",
    )
    report = write_report(
        path=OUTPUT_REPORT,
        roots_by_mu=roots_by_mu,
        tracking_rows=tracking.rows,
        ratios_by_mu=ratios_by_mu,
    )

    print(f"saved global spectrum CSV: {OUTPUT_CSV}")
    print(f"saved global spectrum PNG: {OUTPUT_PNG}")
    print(f"saved global spectrum zoom PNG: {OUTPUT_ZOOM_PNG}")
    print(f"saved global spectrum report: {OUTPUT_REPORT}")
    print(f"thin-rod violations: {report['invalid_segments']}")
    print(
        "max thickness ratio: "
        f"{float(report['max_ratio']):.6g} at mu={float(report['max_ratio_mu']):g}, "
        f"rod {int(report['max_ratio_rod'])}"
    )
    print(
        "tracking warnings: "
        f"low_mac={int(report['low_mac'])}, low_margin={int(report['low_margin'])}, "
        f"unresolved={int(report['unresolved'])}, blocked={int(report['blocked'])}, "
        f"suspicious={int(report['suspicious'])}, requires_refined={int(report['requires_refined'])}"
    )
    print(f"missing/duplicate sorted-position mu points: {int(report['missing_duplicate'])}")
    for branch in (5, 6, 7):
        positions = ", ".join(str(value) for value in report["positions_after"][branch])
        print(f"descendant {branch} accepted sorted positions for mu>=0.7: {positions}")

    return {
        "csv": OUTPUT_CSV,
        "png": OUTPUT_PNG,
        "zoom_png": OUTPUT_ZOOM_PNG,
        "report": OUTPUT_REPORT,
        **report,
    }


if __name__ == "__main__":
    main()
