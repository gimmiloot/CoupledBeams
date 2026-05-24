from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
)
from scripts.lib.thickness_mismatch_mac_tracking import track_mu_branches_shape_mac  # noqa: E402


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA_VALUES = (-0.5, 0.0, 0.5)

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.01

NUM_TRACKED_BRANCHES = 10
NUM_PLOTTED_BRANCHES = 7
NUM_SORTED_ROOTS = 14
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401
MAC_WARNING_THRESHOLD = 0.9

THICKNESS_RATIO_LIMIT = 0.1
OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_PNG = (
    OUTPUT_DIR
    / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_thickness_ratio_split.png"
)
OUTPUT_WARNING_CSV = (
    OUTPUT_DIR
    / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_mac_tracking_warnings.csv"
)
OUTPUT_REPORT = (
    OUTPUT_DIR
    / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_mac_tracking_report.md"
)

VALID_TOL = 1e-12
MU_VALUES = np.round(np.arange(MU_MIN, MU_MAX + 0.5 * MU_STEP, MU_STEP), 10)
BRANCH_CHECK_MUS = (0.70, 0.75, 0.80, 0.85, 0.90)
BRANCH_CHECK_INDICES = (5, 6, 7)

WARNING_FIELDNAMES = [
    "eta",
    "mu_prev",
    "mu",
    "branch_index_from_mu0",
    "Lambda_tracked",
    "mac_sorted_root_index",
    "raw_mac_sorted_root_index",
    "nearest_sorted_root_index",
    "nearest_sorted_Lambda",
    "tracking_step_status",
    "mac_to_previous",
    "raw_mac_to_previous",
    "second_best_mac",
    "frequency_distance_from_previous",
    "frequency_second_nearest_distance",
    "frequency_assignment_margin",
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
]


@dataclass(frozen=True)
class ViolationSegment:
    start_mu: float
    end_mu: float
    rods: tuple[int, ...]


@dataclass(frozen=True)
class ThicknessRatioSummary:
    eta: float
    segments: tuple[ViolationSegment, ...]
    max_ratio: float
    max_ratio_rod: int
    max_ratio_mu: float
    max_ratio_1: float
    max_ratio_2: float


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


def track_one_eta(
    eta: float,
    roots_by_mu: dict[float, np.ndarray],
) -> tuple[np.ndarray, list[dict[str, float | int | str]]]:
    if not np.isclose(float(MU_VALUES[0]), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("This diagnostic seeds tracked branches at mu=0; set MU_MIN = 0.0.")

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
    return result.tracked, result.rows


def thickness_ratios(mu: float, eta: float) -> tuple[float, float]:
    return thickness_to_length_ratios(EPSILON, float(mu), float(eta))


def violating_rods(mu: float, eta: float) -> tuple[int, ...]:
    ratio1, ratio2 = thickness_ratios(mu, eta)
    rods: list[int] = []
    if ratio1 > THICKNESS_RATIO_LIMIT + VALID_TOL:
        rods.append(1)
    if ratio2 > THICKNESS_RATIO_LIMIT + VALID_TOL:
        rods.append(2)
    return tuple(rods)


def valid_mask_for_eta(eta: float) -> np.ndarray:
    return np.array([len(violating_rods(float(mu), float(eta))) == 0 for mu in MU_VALUES], dtype=bool)


def violation_segments_for_eta(eta: float) -> tuple[ViolationSegment, ...]:
    segments: list[ViolationSegment] = []
    start_idx: int | None = None
    current_rods: tuple[int, ...] = ()

    for idx, mu in enumerate(MU_VALUES):
        rods = violating_rods(float(mu), float(eta))
        if rods:
            if start_idx is None:
                start_idx = idx
                current_rods = rods
            elif rods != current_rods:
                segments.append(
                    ViolationSegment(
                        start_mu=float(MU_VALUES[start_idx]),
                        end_mu=float(MU_VALUES[idx - 1]),
                        rods=current_rods,
                    )
                )
                start_idx = idx
                current_rods = rods
        elif start_idx is not None:
            segments.append(
                ViolationSegment(
                    start_mu=float(MU_VALUES[start_idx]),
                    end_mu=float(MU_VALUES[idx - 1]),
                    rods=current_rods,
                )
            )
            start_idx = None
            current_rods = ()

    if start_idx is not None:
        segments.append(
            ViolationSegment(
                start_mu=float(MU_VALUES[start_idx]),
                end_mu=float(MU_VALUES[-1]),
                rods=current_rods,
            )
        )
    return tuple(segments)


def thickness_ratio_summary_for_eta(eta: float) -> ThicknessRatioSummary:
    ratios = np.array([thickness_ratios(float(mu), float(eta)) for mu in MU_VALUES], dtype=float)
    flat_idx = int(np.argmax(ratios))
    mu_idx, rod_idx = np.unravel_index(flat_idx, ratios.shape)
    return ThicknessRatioSummary(
        eta=float(eta),
        segments=violation_segments_for_eta(float(eta)),
        max_ratio=float(ratios[mu_idx, rod_idx]),
        max_ratio_rod=int(rod_idx) + 1,
        max_ratio_mu=float(MU_VALUES[mu_idx]),
        max_ratio_1=float(np.max(ratios[:, 0])),
        max_ratio_2=float(np.max(ratios[:, 1])),
    )


def eta_panel_title(eta: float) -> str:
    if eta > 0.0:
        meaning = "rod 2 thicker"
    elif eta < 0.0:
        meaning = "rod 1 thicker"
    else:
        meaning = "equal radii"
    return rf"$\eta={eta:g}$" + "\n" + meaning


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


def plot_style_split(
    ax: plt.Axes,
    x_values: np.ndarray,
    y_values: np.ndarray,
    valid_mask: np.ndarray,
    *,
    color: str,
    linewidth: float,
) -> None:
    for start, end in contiguous_true_runs(valid_mask):
        if end > start:
            ax.plot(x_values[start : end + 1], y_values[start : end + 1], color=color, lw=linewidth, ls="-")

    invalid_mask = ~valid_mask
    for start, end in contiguous_true_runs(invalid_mask):
        plot_start = max(0, start - 1)
        plot_end = min(len(x_values) - 1, end + 1)
        if plot_end > plot_start:
            ax.plot(
                x_values[plot_start : plot_end + 1],
                y_values[plot_start : plot_end + 1],
                color=color,
                lw=linewidth,
                ls="--",
            )


def plot_tracked(
    *,
    tracked_by_eta: dict[float, np.ndarray],
    output: Path,
) -> None:
    fig, axes = plt.subplots(1, len(ETA_VALUES), figsize=(11.4, 4.45), sharey=True)
    axes = np.atleast_1d(axes)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, eta in zip(axes, ETA_VALUES, strict=True):
        tracked = tracked_by_eta[float(eta)]
        valid_mask = valid_mask_for_eta(float(eta))
        for branch_idx in range(NUM_PLOTTED_BRANCHES):
            color = colors[branch_idx % len(colors)]
            plot_style_split(
                ax,
                MU_VALUES,
                tracked[branch_idx],
                valid_mask,
                color=color,
                linewidth=1.75,
            )
        ax.set_title(eta_panel_title(float(eta)), fontsize=10)
        ax.set_xlabel(r"$\mu$")
        ax.grid(True, color="0.88", linewidth=0.6)

    branch_handles = [
        Line2D([0], [0], color=colors[idx % len(colors)], lw=1.75, label=f"branch {idx + 1}")
        for idx in range(NUM_PLOTTED_BRANCHES)
    ]
    style_handles = [
        Line2D([0], [0], color="0.2", lw=1.75, ls="-", label=r"solid: $2r_i/l_i \leq 0.1$"),
        Line2D([0], [0], color="0.2", lw=1.75, ls="--", label="dashed: criterion violated"),
    ]

    axes[0].set_ylabel(r"$\Lambda$")
    fig.legend(
        handles=branch_handles + style_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=4,
        fontsize=8,
        frameon=False,
    )
    fig.suptitle(
        rf"Tracked Lambda(mu), beta={BETA_DEG:g}$^\circ$, epsilon={EPSILON:g}; "
        rf"thickness split by $2r_i/l_i$",
        fontsize=12,
    )
    fig.tight_layout(rect=(0.0, 0.16, 1.0, 0.92))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def print_tracking_summary(rows: list[dict[str, float | int | str]]) -> None:
    low_mac = [row for row in rows if "low_mac" in str(row["tracking_step_status"])]
    switches = [
        row
        for row in rows
        if int(row["mac_sorted_root_index"]) != int(row["branch_index_from_mu0"])
    ]
    disagreements = [row for row in rows if str(row["frequency_mac_disagreement"]) == "yes"]
    suspicious = [row for row in rows if str(row["suspicious_assignment"]) == "yes"]
    refined = [row for row in rows if str(row["requires_refined_check"]) == "yes"]
    print("tracking method: analytic shape-MAC between adjacent mu steps")
    print(f"MAC warning threshold: {MAC_WARNING_THRESHOLD:g}")
    print(f"MAC current-sorted-index switch rows: {len(switches)}")
    print(f"low-MAC assignment rows: {len(low_mac)}")
    print(f"nearest-frequency vs MAC disagreement rows: {len(disagreements)}")
    print(f"suspicious assignment rows: {len(suspicious)}")
    print(f"rows requiring refined check: {len(refined)}")


def rods_label(rods: tuple[int, ...]) -> str:
    return "both" if rods == (1, 2) else ", ".join(str(rod) for rod in rods)


def print_thickness_ratio_summary(summaries: list[ThicknessRatioSummary]) -> None:
    print(f"thickness-ratio criterion: 2*r_i/l_i <= {THICKNESS_RATIO_LIMIT:g}")
    for summary in summaries:
        if not summary.segments:
            print(
                f"eta={summary.eta:g}: no violations; "
                f"max ratio={summary.max_ratio:.6g} "
                f"(rod {summary.max_ratio_rod}, mu={summary.max_ratio_mu:g})"
            )
            continue

        print(f"eta={summary.eta:g}: violations detected")
        for segment in summary.segments:
            print(
                f"  mu={segment.start_mu:g}..{segment.end_mu:g}, "
                f"rod(s)={rods_label(segment.rods)}"
            )
        print(
            f"  max ratio={summary.max_ratio:.6g} "
            f"(rod {summary.max_ratio_rod}, mu={summary.max_ratio_mu:g}); "
            f"max rod1={summary.max_ratio_1:.6g}, max rod2={summary.max_ratio_2:.6g}"
        )


def write_warning_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=WARNING_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def eta_branch_rows(
    rows: list[dict[str, float | int | str]],
    *,
    eta: float,
    mus: tuple[float, ...],
    branches: tuple[int, ...],
) -> list[dict[str, float | int | str]]:
    selected = []
    for row in rows:
        if not np.isclose(float(row["eta"]), float(eta), rtol=0.0, atol=1e-12):
            continue
        if not any(np.isclose(float(row["mu"]), float(mu), rtol=0.0, atol=1e-12) for mu in mus):
            continue
        if int(row["branch_index_from_mu0"]) not in branches:
            continue
        selected.append(row)
    selected.sort(key=lambda row: (float(row["mu"]), int(row["branch_index_from_mu0"])))
    return selected


def branch_mapping_lines(rows: list[dict[str, float | int | str]]) -> list[str]:
    selected = eta_branch_rows(rows, eta=0.5, mus=BRANCH_CHECK_MUS, branches=BRANCH_CHECK_INDICES)
    by_mu: dict[float, list[dict[str, float | int | str]]] = {}
    for row in selected:
        by_mu.setdefault(float(row["mu"]), []).append(row)

    lines: list[str] = []
    for mu in BRANCH_CHECK_MUS:
        local = sorted(by_mu.get(float(mu), []), key=lambda row: int(row["branch_index_from_mu0"]))
        if not local:
            lines.append(f"mu={mu:g}: no rows")
            continue
        mapping = ", ".join(
            (
                f"branch {int(row['branch_index_from_mu0'])}"
                f" current sorted position {int(row['mac_sorted_root_index'])}"
                f" (MAC={float(row['mac_to_previous']):.4f}, "
                f"jump={int(row['sorted_position_jump'])}, "
                f"suspicious={row['suspicious_assignment']})"
            )
            for row in local
        )
        lines.append(f"mu={mu:g}: {mapping}")
    return lines


def write_tracking_report(
    *,
    path: Path,
    rows: list[dict[str, float | int | str]],
    warning_rows: list[dict[str, float | int | str]],
    summaries: list[ThicknessRatioSummary],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    low_mac = [row for row in rows if "low_mac" in str(row["tracking_step_status"])]
    suspicious = [row for row in rows if str(row["suspicious_assignment"]) == "yes"]
    refined = [row for row in rows if str(row["requires_refined_check"]) == "yes"]
    switches = [
        row
        for row in rows
        if int(row["mac_sorted_root_index"]) != int(row["branch_index_from_mu0"])
    ]
    lines = [
        "# Large-Eta Thickness-Mismatch MAC Tracking Diagnostic",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta values: {', '.join(f'{float(value):g}' for value in ETA_VALUES)}",
        f"- mu range: {float(MU_VALUES[0]):g} .. {float(MU_VALUES[-1]):g}",
        f"- tracked branches: first {NUM_TRACKED_BRANCHES}",
        f"- plotted branches: first {NUM_PLOTTED_BRANCHES}",
        f"- sorted root candidates per step: {NUM_SORTED_ROOTS}",
        f"- shape samples per analytic arm: {NUM_SHAPE_SAMPLES}",
        "",
        "## Method",
        "",
        "Branches are seeded at `mu=0` for each fixed eta. Each next mu step",
        "reconstructs analytic null-vector mode shapes for all sorted-root",
        "candidates and assigns branches by maximizing the adjacent-step MAC",
        "matrix. Low-MAC assignments are marked for review. Nearest-frequency",
        "assignment is computed only as a comparator; MAC-vs-frequency",
        "disagreements are stored as warnings.",
        "",
        "## Warning Summary",
        "",
        f"- MAC current-sorted-index switch rows: {len(switches)}",
        f"- low-MAC rows below threshold {MAC_WARNING_THRESHOLD:g}: {len(low_mac)}",
        f"- nearest-frequency vs MAC disagreement rows: {len(warning_rows)}",
        f"- suspicious assignment rows: {len(suspicious)}",
        f"- rows requiring refined check: {len(refined)}",
        f"- warning CSV: `{OUTPUT_WARNING_CSV.relative_to(REPO_ROOT)}`",
        "",
        "## Eta=0.5 Branches 5--7",
        "",
    ]
    lines.extend(f"- {line}" for line in branch_mapping_lines(rows))
    lines.extend(
        [
            "",
            "Interpretation: these are descendant branches, and the listed sorted",
            "positions are diagnostic metadata. Any branch-5--7 sorted-position",
            "jump in this coarse large-eta plot is a suspicious assignment until",
            "a refined local mu-grid audit confirms it. Do not interpret this",
            "table as automatic renumbering of the branches.",
        ]
    )

    lines.extend(["", "## Thin-Rod Criterion", ""])
    for summary in summaries:
        if not summary.segments:
            lines.append(
                f"- eta={summary.eta:g}: no violations; max ratio={summary.max_ratio:.6g} "
                f"(rod {summary.max_ratio_rod}, mu={summary.max_ratio_mu:g})"
            )
        else:
            segment_text = "; ".join(
                f"mu={segment.start_mu:g}..{segment.end_mu:g}, rod(s)={rods_label(segment.rods)}"
                for segment in summary.segments
            )
            lines.append(
                f"- eta={summary.eta:g}: violations: {segment_text}; "
                f"max ratio={summary.max_ratio:.6g} "
                f"(rod {summary.max_ratio_rod}, mu={summary.max_ratio_mu:g})"
            )

    lines.extend(
        [
            "",
            "This is diagnostic-only. It does not change the article files, article",
            "figures, baseline equal-radius determinant, old solvers, or FEM model.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def print_branch_5_7_summary(rows: list[dict[str, float | int | str]]) -> None:
    print("eta=0.5 branch 5--7 diagnostic current sorted positions:")
    for line in branch_mapping_lines(rows):
        print(f"  {line}")
    print("warning: branch 5--7 sorted-position jumps require the refined branch-identity audit.")


def main() -> dict[str, object]:
    sorted_roots = sorted_roots_by_eta_mu()
    tracked_rows: list[dict[str, float | int | str]] = []
    tracked_by_eta: dict[float, np.ndarray] = {}

    for eta in ETA_VALUES:
        tracked, rows = track_one_eta(float(eta), sorted_roots[float(eta)])
        tracked_by_eta[float(eta)] = tracked
        tracked_rows.extend(rows)

    summaries = [thickness_ratio_summary_for_eta(float(eta)) for eta in ETA_VALUES]
    warning_rows = [
        row
        for row in tracked_rows
        if str(row["frequency_mac_disagreement"]) == "yes" or str(row["requires_refined_check"]) == "yes"
    ]
    write_warning_csv(OUTPUT_WARNING_CSV, warning_rows)
    write_tracking_report(path=OUTPUT_REPORT, rows=tracked_rows, warning_rows=warning_rows, summaries=summaries)
    plot_tracked(tracked_by_eta=tracked_by_eta, output=OUTPUT_PNG)

    print(f"saved tracked Lambda(mu) thickness-ratio plot: {OUTPUT_PNG}")
    print(f"saved MAC/frequency disagreement warnings: {OUTPUT_WARNING_CSV}")
    print(f"saved MAC tracking report: {OUTPUT_REPORT}")
    print_tracking_summary(tracked_rows)
    print_branch_5_7_summary(tracked_rows)
    print_thickness_ratio_summary(summaries)

    return {
        "output_png": OUTPUT_PNG,
        "warning_csv": OUTPUT_WARNING_CSV,
        "report": OUTPUT_REPORT,
        "tracked_by_eta": tracked_by_eta,
        "tracking_rows": tracked_rows,
        "warning_rows": warning_rows,
        "thickness_ratio_summaries": summaries,
    }


if __name__ == "__main__":
    main()
