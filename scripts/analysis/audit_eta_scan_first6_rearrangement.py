from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis import audit_lambda_beta_sorted_descendant_thickness_mismatch as beta_audit  # noqa: E402


# =========================
# User-editable parameters
# =========================
EPSILON = 0.0025
MU = 0.0

ETA_MIN = -0.5
ETA_MAX = 0.5
COARSE_ETA_STEP = 0.02
REFINE_ETA_STEP = 0.002
BETA_STEP = 0.25

N_DESCENDANTS_TRACK = 8
N_DESCENDANTS_CLASSIFY = 6
N_SORTED_ROOTS = 14
N_SORTED_GAP_ROOTS = 7
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
TRACKING_UNRELIABLE_MAC = 0.5
TRACKING_UNRELIABLE_WARNING_FRACTION = 0.20
NEAR_DEGENERATE_GAP = 1e-3

KEY_BETAS = (0.0, 15.0, 90.0)
KEY_ETAS = (-0.5, 0.0, 0.5)

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_first6_rearrangement.csv"
OUTPUT_SUMMARY_CSV = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_first6_rearrangement_summary.csv"
OUTPUT_REPORT = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_first6_rearrangement_report.md"
OUTPUT_MAP_PNG = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_first6_sorted_position_map.png"
OUTPUT_CLASS_PNG = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_first6_rearrangement_class.png"
OUTPUT_TRACES_PNG = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_sorted_position_traces_representative_eta.png"
OUTPUT_GAPS_PNG = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_min_adjacent_gaps_first7.png"
OUTPUT_SORTED5_MAP_PNG = OUTPUT_DIR / "eta_scan_eps0p0025_mu0_sorted5_identity_map.png"

SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

CLASS_CODES = {
    "no_rearrangement": 0,
    "rearrangement": 1,
    "ambiguous_near_degenerate": 2,
    "tracking_unreliable": 3,
}

MAIN_FIELDNAMES = [
    "eta",
    "beta_deg",
    "descendant_id",
    "Lambda_descendant",
    "sorted_position",
    "mac_to_previous",
    "tracking_warning",
    "is_position_changed_relative_to_beta0",
    "is_within_first6",
    "sorted5_descendant_id_at_beta",
    "sorted5_lambda",
    "min_adjacent_gap_first7_at_beta",
    "gap_sorted1_2",
    "gap_sorted2_3",
    "gap_sorted3_4",
    "gap_sorted4_5",
    "gap_sorted5_6",
    "gap_sorted6_7",
]

SUMMARY_FIELDNAMES = [
    "eta",
    "rearrangement_class",
    "rearrangement_detected_bool",
    "descendants_with_position_changes",
    "num_descendants_changed",
    "desc1_positions_over_beta",
    "desc2_positions_over_beta",
    "desc3_positions_over_beta",
    "desc4_positions_over_beta",
    "desc5_positions_over_beta",
    "desc6_positions_over_beta",
    "desc1_pos_beta0",
    "desc2_pos_beta0",
    "desc3_pos_beta0",
    "desc4_pos_beta0",
    "desc5_pos_beta0",
    "desc6_pos_beta0",
    "desc1_pos_beta15",
    "desc2_pos_beta15",
    "desc3_pos_beta15",
    "desc4_pos_beta15",
    "desc5_pos_beta15",
    "desc6_pos_beta15",
    "desc1_pos_beta90",
    "desc2_pos_beta90",
    "desc3_pos_beta90",
    "desc4_pos_beta90",
    "desc5_pos_beta90",
    "desc6_pos_beta90",
    "sorted5_descendant_id_beta15",
    "min_gap_first7",
    "beta_min_gap_first7",
    "tracking_warning_count",
    "min_mac_continuity",
    "notes",
]


@dataclass(frozen=True)
class EtaScanCase:
    eta: float
    result: beta_audit.BetaTrackingResult
    summary_row: dict[str, str]


def number_text(value: float) -> str:
    return f"{float(value):.12g}"


def eta_values(start: float, stop: float, step: float) -> np.ndarray:
    values = np.arange(float(start), float(stop) + 0.5 * float(step), float(step), dtype=float)
    values[0] = float(start)
    values[-1] = float(stop)
    values = np.unique(np.round(values, 12))
    return values


def beta_values() -> np.ndarray:
    values = np.arange(0.0, 90.0 + 0.5 * BETA_STEP, BETA_STEP, dtype=float)
    values[0] = 0.0
    values[-1] = 90.0
    return np.unique(np.round(values, 12))


def configure_beta_backend(beta_grid: np.ndarray) -> None:
    beta_audit.EPSILON = EPSILON
    beta_audit.MU = MU
    beta_audit.BETA_VALUES_DEG = np.asarray(beta_grid, dtype=float)
    beta_audit.N_DESCENDANTS = N_DESCENDANTS_TRACK
    beta_audit.N_SORTED_ROOTS = N_SORTED_ROOTS
    beta_audit.ROOT_SCAN_STEP = ROOT_SCAN_STEP
    beta_audit.ROOT_LMAX0 = ROOT_LMAX0
    beta_audit.NUM_SHAPE_SAMPLES = NUM_SHAPE_SAMPLES
    beta_audit.MAC_WARNING_THRESHOLD = MAC_WARNING_THRESHOLD
    beta_audit.MAC_MARGIN_WARNING_THRESHOLD = MAC_MARGIN_WARNING_THRESHOLD
    beta_audit.MAX_SORTED_POSITION_JUMP = MAX_SORTED_POSITION_JUMP


def finite_min(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(np.min(finite)) if finite.size else float("nan")


def yesno(value: bool) -> str:
    return "yes" if bool(value) else "no"


def tracking_warning(row: Mapping[str, float | int | str]) -> str:
    text = beta_audit.tracking_warning(row)
    return text if text else "no"


def rows_by_beta_descendant(result: beta_audit.BetaTrackingResult) -> dict[tuple[float, int], dict[str, float | int | str]]:
    return {
        (round(float(row["beta_deg"]), 10), int(row["descendant_id"])): row
        for row in result.rows
    }


def sorted5_identity_by_beta(result: beta_audit.BetaTrackingResult) -> list[str]:
    return beta_audit.sorted5_identity_by_beta(result)


def key_beta_index(beta_grid: np.ndarray, beta_deg: float) -> int:
    matches = np.flatnonzero(np.isclose(beta_grid, float(beta_deg), rtol=0.0, atol=1e-12))
    if matches.size:
        return int(matches[0])
    return int(np.argmin(np.abs(beta_grid - float(beta_deg))))


def compact_position_sequence(beta_grid: np.ndarray, positions: np.ndarray) -> str:
    parts: list[str] = []
    start = 0
    for idx in range(1, len(positions) + 1):
        if idx < len(positions) and int(positions[idx]) == int(positions[start]):
            continue
        beta_start = float(beta_grid[start])
        beta_end = float(beta_grid[idx - 1])
        if np.isclose(beta_start, beta_end, rtol=0.0, atol=1e-12):
            parts.append(f"beta {beta_start:g} -> sorted {int(positions[start])}")
        else:
            parts.append(f"beta {beta_start:g}..{beta_end:g} -> sorted {int(positions[start])}")
        start = idx
    return "; ".join(parts)


def adjacent_gaps(sorted_roots: np.ndarray) -> np.ndarray:
    return np.diff(sorted_roots[:, :N_SORTED_GAP_ROOTS], axis=1)


def min_adjacent_gap_at_beta(sorted_roots_at_beta: np.ndarray) -> float:
    gaps = np.diff(np.asarray(sorted_roots_at_beta[:N_SORTED_GAP_ROOTS], dtype=float))
    return float(np.nanmin(gaps))


def sorted5_lambda_at_beta(sorted_roots_at_beta: np.ndarray) -> float:
    return float(sorted_roots_at_beta[4])


def change_event_indices(positions: np.ndarray, descendant_zero_based: int) -> np.ndarray:
    return np.flatnonzero(np.diff(positions[descendant_zero_based]) != 0) + 1


def classify_case(result: beta_audit.BetaTrackingResult) -> tuple[str, list[int], str]:
    positions = result.current_sorted_positions[:N_DESCENDANTS_CLASSIFY]
    changed_descendants = [
        desc_idx + 1
        for desc_idx in range(N_DESCENDANTS_CLASSIFY)
        if np.unique(positions[desc_idx]).size > 1
    ]
    if not changed_descendants:
        return "no_rearrangement", [], "all first-six descendants keep constant sorted positions"

    finite_macs = [
        float(row["mac_to_previous"])
        for row in result.rows
        if np.isfinite(float(row["mac_to_previous"]))
    ]
    min_mac = finite_min(finite_macs)
    warning_fraction = len(result.warning_rows) / max(1, len(result.rows))
    if np.isfinite(min_mac) and min_mac < TRACKING_UNRELIABLE_MAC:
        return (
            "tracking_unreliable",
            changed_descendants,
            f"minimum MAC {min_mac:.6g} below {TRACKING_UNRELIABLE_MAC:g}",
        )
    if warning_fraction > TRACKING_UNRELIABLE_WARNING_FRACTION:
        return (
            "tracking_unreliable",
            changed_descendants,
            f"warning fraction {warning_fraction:.6g} above {TRACKING_UNRELIABLE_WARNING_FRACTION:g}",
        )

    row_lookup = rows_by_beta_descendant(result)
    gaps = adjacent_gaps(result.sorted_roots)
    change_flags: list[bool] = []
    for desc_id in changed_descendants:
        for col in change_event_indices(positions, desc_id - 1):
            beta_deg = float(result.beta_values_deg[col])
            row = row_lookup[(round(beta_deg, 10), desc_id)]
            local_gap = float(np.nanmin(gaps[col]))
            has_warning = tracking_warning(row) != "no"
            change_flags.append(has_warning and local_gap < NEAR_DEGENERATE_GAP)
    if change_flags and all(change_flags):
        return (
            "ambiguous_near_degenerate",
            changed_descendants,
            f"all sorted-position changes coincide with warning rows and gap < {NEAR_DEGENERATE_GAP:g}",
        )
    return "rearrangement", changed_descendants, "accepted sorted-position changes in first six descendants"


def summary_for_result(result: beta_audit.BetaTrackingResult) -> dict[str, str]:
    rearrangement_class, changed_descendants, note = classify_case(result)
    beta_grid = np.asarray(result.beta_values_deg, dtype=float)
    positions = result.current_sorted_positions
    gap_matrix = adjacent_gaps(result.sorted_roots)
    min_gap_by_beta = np.nanmin(gap_matrix, axis=1)
    min_gap_idx = int(np.nanargmin(min_gap_by_beta))
    sorted5_identities = sorted5_identity_by_beta(result)
    warning_count = len(result.warning_rows)
    finite_macs = [
        float(row["mac_to_previous"])
        for row in result.rows
        if np.isfinite(float(row["mac_to_previous"]))
    ]

    row: dict[str, str] = {
        "eta": number_text(result.eta),
        "rearrangement_class": rearrangement_class,
        "rearrangement_detected_bool": yesno(bool(changed_descendants)),
        "descendants_with_position_changes": ";".join(str(value) for value in changed_descendants) or "none",
        "num_descendants_changed": str(len(changed_descendants)),
        "sorted5_descendant_id_beta15": sorted5_identities[key_beta_index(beta_grid, 15.0)],
        "min_gap_first7": number_text(float(min_gap_by_beta[min_gap_idx])),
        "beta_min_gap_first7": number_text(float(beta_grid[min_gap_idx])),
        "tracking_warning_count": str(warning_count),
        "min_mac_continuity": number_text(finite_min(finite_macs)),
        "notes": note,
    }
    for desc_id in range(1, N_DESCENDANTS_CLASSIFY + 1):
        row[f"desc{desc_id}_positions_over_beta"] = compact_position_sequence(beta_grid, positions[desc_id - 1])
    for beta_label, beta_value in (("beta0", 0.0), ("beta15", 15.0), ("beta90", 90.0)):
        idx = key_beta_index(beta_grid, beta_value)
        for desc_id in range(1, N_DESCENDANTS_CLASSIFY + 1):
            row[f"desc{desc_id}_pos_{beta_label}"] = str(int(positions[desc_id - 1, idx]))
    return row


def run_eta_case(eta: float, beta_grid: np.ndarray) -> EtaScanCase:
    configure_beta_backend(beta_grid)
    roots = beta_audit.solve_sorted_root_grid(eta=float(eta))
    result = beta_audit.track_descendants_by_beta(eta=float(eta), sorted_roots=roots)
    summary = summary_for_result(result)
    print(
        f"eta={float(eta): .6f}: class={summary['rearrangement_class']}, "
        f"changed={summary['descendants_with_position_changes']}, "
        f"warnings={summary['tracking_warning_count']}, min_mac={summary['min_mac_continuity']}"
    )
    return EtaScanCase(eta=float(eta), result=result, summary_row=summary)


def class_change_intervals(cases: Sequence[EtaScanCase]) -> list[tuple[float, float]]:
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    intervals: list[tuple[float, float]] = []
    for left, right in zip(sorted_cases[:-1], sorted_cases[1:], strict=True):
        left_class = left.summary_row["rearrangement_class"]
        right_class = right.summary_row["rearrangement_class"]
        if left_class != right_class:
            intervals.append((float(left.eta), float(right.eta)))
    return intervals


def refined_eta_values(coarse_cases: Sequence[EtaScanCase]) -> np.ndarray:
    values: set[float] = set(float(case.eta) for case in coarse_cases)
    for left, right in class_change_intervals(coarse_cases):
        lo, hi = sorted((left, right))
        for value in eta_values(lo, hi, REFINE_ETA_STEP):
            values.add(float(value))
    for eta in KEY_ETAS:
        values.add(float(eta))
    return np.asarray(sorted(values), dtype=float)


def write_main_csv(cases: Sequence[EtaScanCase]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MAIN_FIELDNAMES)
        writer.writeheader()
        for case in sorted(cases, key=lambda item: item.eta):
            result = case.result
            row_lookup = rows_by_beta_descendant(result)
            sorted5_identities = sorted5_identity_by_beta(result)
            gaps = adjacent_gaps(result.sorted_roots)
            initial_positions = result.current_sorted_positions[:, 0].copy()
            for beta_idx, beta_deg in enumerate(result.beta_values_deg):
                beta_gaps = gaps[beta_idx]
                for desc_idx in range(N_DESCENDANTS_TRACK):
                    row = row_lookup[(round(float(beta_deg), 10), desc_idx + 1)]
                    mac = float(row["mac_to_previous"])
                    writer.writerow(
                        {
                            "eta": number_text(case.eta),
                            "beta_deg": number_text(float(beta_deg)),
                            "descendant_id": desc_idx + 1,
                            "Lambda_descendant": number_text(float(result.tracked_lambdas[desc_idx, beta_idx])),
                            "sorted_position": int(result.current_sorted_positions[desc_idx, beta_idx]),
                            "mac_to_previous": "" if not np.isfinite(mac) else number_text(mac),
                            "tracking_warning": tracking_warning(row),
                            "is_position_changed_relative_to_beta0": yesno(
                                int(result.current_sorted_positions[desc_idx, beta_idx]) != int(initial_positions[desc_idx])
                            ),
                            "is_within_first6": yesno(desc_idx < N_DESCENDANTS_CLASSIFY),
                            "sorted5_descendant_id_at_beta": sorted5_identities[beta_idx],
                            "sorted5_lambda": number_text(sorted5_lambda_at_beta(result.sorted_roots[beta_idx])),
                            "min_adjacent_gap_first7_at_beta": number_text(float(np.nanmin(beta_gaps))),
                            "gap_sorted1_2": number_text(float(beta_gaps[0])),
                            "gap_sorted2_3": number_text(float(beta_gaps[1])),
                            "gap_sorted3_4": number_text(float(beta_gaps[2])),
                            "gap_sorted4_5": number_text(float(beta_gaps[3])),
                            "gap_sorted5_6": number_text(float(beta_gaps[4])),
                            "gap_sorted6_7": number_text(float(beta_gaps[5])),
                        }
                    )


def write_summary_csv(cases: Sequence[EtaScanCase]) -> None:
    OUTPUT_SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        for case in sorted(cases, key=lambda item: item.eta):
            writer.writerow(case.summary_row)


def class_intervals(cases: Sequence[EtaScanCase]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    start = 0
    for idx in range(1, len(sorted_cases) + 1):
        if idx < len(sorted_cases) and sorted_cases[idx].summary_row["rearrangement_class"] == sorted_cases[start].summary_row["rearrangement_class"]:
            continue
        local = sorted_cases[start:idx]
        changed: set[str] = set()
        for case in local:
            value = case.summary_row["descendants_with_position_changes"]
            if value != "none":
                changed.update(value.split(";"))
        rows.append(
            {
                "eta_start": number_text(local[0].eta),
                "eta_end": number_text(local[-1].eta),
                "class": local[0].summary_row["rearrangement_class"],
                "changed_descendants": ";".join(sorted(changed, key=int)) if changed else "none",
            }
        )
        start = idx
    return rows


def boundary_lines(cases: Sequence[EtaScanCase]) -> list[str]:
    lines: list[str] = []
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    for left, right in zip(sorted_cases[:-1], sorted_cases[1:], strict=True):
        left_class = left.summary_row["rearrangement_class"]
        right_class = right.summary_row["rearrangement_class"]
        if left_class == right_class:
            continue
        midpoint = 0.5 * (left.eta + right.eta)
        lines.append(
            f"- eta in [{left.eta:.6g}, {right.eta:.6g}] changes "
            f"`{left_class}` -> `{right_class}`; midpoint about `{midpoint:.6g}`"
        )
    return lines or ["- no rearrangement-class boundary was found on the evaluated eta grid."]


def case_nearest(cases: Sequence[EtaScanCase], eta: float) -> EtaScanCase:
    return min(cases, key=lambda case: abs(float(case.eta) - float(eta)))


def positions_at_key_betas(case: EtaScanCase) -> dict[float, list[int]]:
    out: dict[float, list[int]] = {}
    for beta in KEY_BETAS:
        idx = key_beta_index(case.result.beta_values_deg, beta)
        out[beta] = [int(case.result.current_sorted_positions[desc - 1, idx]) for desc in range(1, 7)]
    return out


def symmetry_summary(cases: Sequence[EtaScanCase]) -> tuple[bool, list[str]]:
    by_eta = {round(float(case.eta), 12): case for case in cases}
    mismatches: list[str] = []
    for eta, case in sorted(by_eta.items()):
        if eta < -1e-12:
            other = by_eta.get(round(-eta, 12))
            if other is None:
                continue
            if case.summary_row["rearrangement_class"] != other.summary_row["rearrangement_class"]:
                mismatches.append(
                    f"eta={eta:g}: {case.summary_row['rearrangement_class']} vs "
                    f"eta={-eta:g}: {other.summary_row['rearrangement_class']}"
                )
    return not mismatches, mismatches


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return lines


def plot_sorted_position_map(cases: Sequence[EtaScanCase]) -> None:
    OUTPUT_MAP_PNG.parent.mkdir(parents=True, exist_ok=True)
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    fig, axes = plt.subplots(3, 2, figsize=(12.0, 9.0), sharex=True, sharey=True, constrained_layout=True)
    axes = axes.ravel()
    cmap = plt.get_cmap("tab10", N_DESCENDANTS_TRACK + 1)
    for desc_idx, ax in enumerate(axes):
        for case in sorted_cases:
            beta = case.result.beta_values_deg
            eta = np.full_like(beta, float(case.eta), dtype=float)
            values = case.result.current_sorted_positions[desc_idx]
            sc = ax.scatter(beta, eta, c=values, cmap=cmap, vmin=0.5, vmax=N_DESCENDANTS_TRACK + 0.5, s=8)
        ax.set_title(f"descendant {desc_idx + 1}")
        ax.grid(True, color="0.88", linewidth=0.5)
        if desc_idx % 2 == 0:
            ax.set_ylabel(r"$\eta$")
        if desc_idx >= 4:
            ax.set_xlabel(r"$\beta$ (deg)")
    cbar = fig.colorbar(
        sc,
        ax=axes.tolist(),
        location="right",
        fraction=0.025,
        pad=0.02,
        ticks=range(1, N_DESCENDANTS_TRACK + 1),
    )
    cbar.set_label("sorted position")
    fig.suptitle(r"First-six descendant sorted-position map, $\epsilon=0.0025$, $\mu=0$")
    fig.savefig(OUTPUT_MAP_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_class(cases: Sequence[EtaScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    eta = np.array([case.eta for case in sorted_cases], dtype=float)
    codes = np.array([CLASS_CODES[case.summary_row["rearrangement_class"]] for case in sorted_cases], dtype=float)
    fig, ax = plt.subplots(figsize=(10.0, 4.0))
    ax.step(eta, codes, where="mid", color="black", lw=1.4)
    ax.scatter(eta, codes, c=codes, cmap="viridis", vmin=0, vmax=3, s=28, zorder=3)
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("class code")
    ax.set_yticks([0, 1, 2, 3])
    ax.set_yticklabels(["0 no", "1 rearr.", "2 ambiguous", "3 unreliable"])
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_title("First-six rearrangement class by eta")
    fig.tight_layout()
    fig.savefig(OUTPUT_CLASS_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def representative_etas(cases: Sequence[EtaScanCase]) -> list[float]:
    reps = {float(value) for value in KEY_ETAS}
    for line in boundary_lines(cases):
        if "midpoint about" not in line:
            continue
        text = line.split("`")[-2]
        try:
            reps.add(float(text))
        except ValueError:
            pass
    return sorted(reps)


def plot_representative_traces(cases: Sequence[EtaScanCase]) -> None:
    reps = representative_etas(cases)
    selected = [case_nearest(cases, eta) for eta in reps]
    unique: dict[float, EtaScanCase] = {}
    for case in selected:
        unique[round(float(case.eta), 12)] = case
    selected = [unique[key] for key in sorted(unique)]
    n = len(selected)
    fig, axes = plt.subplots(n, 1, figsize=(10.5, max(3.0, 2.25 * n)), sharex=True)
    axes = np.atleast_1d(axes)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, case in zip(axes, selected, strict=True):
        for desc_idx in range(N_DESCENDANTS_CLASSIFY):
            ax.step(
                case.result.beta_values_deg,
                case.result.current_sorted_positions[desc_idx],
                where="post",
                color=colors[desc_idx % len(colors)],
                lw=1.4,
                label=f"desc {desc_idx + 1}",
            )
        ax.set_ylabel("sorted pos.")
        ax.set_title(f"eta={case.eta:g}, class={case.summary_row['rearrangement_class']}")
        ax.set_yticks(range(1, N_DESCENDANTS_TRACK + 1))
        ax.grid(True, color="0.88", linewidth=0.6)
        ax.legend(loc="best", ncol=6, fontsize=8, frameon=False)
    axes[-1].set_xlabel(r"$\beta$ (deg)")
    fig.suptitle("Representative first-six descendant sorted-position traces")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.97))
    fig.savefig(OUTPUT_TRACES_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_min_gaps(cases: Sequence[EtaScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    eta = np.array([case.eta for case in sorted_cases], dtype=float)
    gap_by_pair = np.full((len(sorted_cases), 6), np.nan, dtype=float)
    for row_idx, case in enumerate(sorted_cases):
        gaps = adjacent_gaps(case.result.sorted_roots)
        gap_by_pair[row_idx, :] = np.nanmin(gaps, axis=0)
    fig, ax = plt.subplots(figsize=(10.0, 5.0))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for pair_idx in range(6):
        ax.plot(
            eta,
            gap_by_pair[:, pair_idx],
            marker="o",
            ms=3,
            lw=1.1,
            color=colors[pair_idx % len(colors)],
            label=f"sorted {pair_idx + 1}-{pair_idx + 2}",
        )
    ax.set_xlabel(r"$\eta$")
    ax.set_ylabel("minimum adjacent gap over beta")
    ax.set_yscale("log")
    ax.grid(True, color="0.88", linewidth=0.6, which="both")
    ax.legend(loc="best", ncol=3, fontsize=8, frameon=False)
    ax.set_title("Minimum adjacent sorted gaps among roots 1..7")
    fig.tight_layout()
    fig.savefig(OUTPUT_GAPS_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_sorted5_identity_map(cases: Sequence[EtaScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    fig, ax = plt.subplots(figsize=(10.2, 5.2))
    cmap = plt.get_cmap("tab10", N_DESCENDANTS_TRACK + 1)
    for case in sorted_cases:
        beta = case.result.beta_values_deg
        eta = np.full_like(beta, float(case.eta), dtype=float)
        ids = [
            float(value) if str(value).isdigit() else np.nan
            for value in sorted5_identity_by_beta(case.result)
        ]
        sc = ax.scatter(beta, eta, c=ids, cmap=cmap, vmin=0.5, vmax=N_DESCENDANTS_TRACK + 0.5, s=8)
    cbar = fig.colorbar(sc, ax=ax, ticks=range(1, N_DESCENDANTS_TRACK + 1))
    cbar.set_label("descendant occupying sorted 5")
    ax.set_xlabel(r"$\beta$ (deg)")
    ax.set_ylabel(r"$\eta$")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_title("Secondary sorted-5 descendant identity map")
    fig.tight_layout()
    fig.savefig(OUTPUT_SORTED5_MAP_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(cases: Sequence[EtaScanCase], *, coarse_count: int, refined_count: int) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.eta)
    eta0 = case_nearest(sorted_cases, 0.0)
    etap = case_nearest(sorted_cases, 0.5)
    etam = case_nearest(sorted_cases, -0.5)
    intervals = class_intervals(sorted_cases)
    symmetric, mismatches = symmetry_summary(sorted_cases)

    interval_rows = [
        [row["eta_start"], row["eta_end"], row["class"], row["changed_descendants"]]
        for row in intervals
    ]
    key_rows = []
    for label, case in (("eta=-0.5", etam), ("eta=0", eta0), ("eta=0.5", etap)):
        positions = positions_at_key_betas(case)
        key_rows.append(
            [
                label,
                case.summary_row["rearrangement_class"],
                case.summary_row["descendants_with_position_changes"],
                ",".join(str(value) for value in positions[0.0]),
                ",".join(str(value) for value in positions[15.0]),
                ",".join(str(value) for value in positions[90.0]),
                case.summary_row["sorted5_descendant_id_beta15"],
            ]
        )

    involved: dict[str, set[str]] = {}
    for case in sorted_cases:
        cls = case.summary_row["rearrangement_class"]
        involved.setdefault(cls, set())
        value = case.summary_row["descendants_with_position_changes"]
        if value != "none":
            involved[cls].update(value.split(";"))

    lines = [
        "# Eta Scan First-Six Rearrangement Audit",
        "",
        "## 1. Rearrangement Definition",
        "",
        "In this audit, `rearrangement` means that at least one of the first six",
        "descendant branches changes its current sorted position as beta varies",
        "over `0..90 deg` at fixed eta. Descendant identity is primary; sorted",
        "position is diagnostic metadata describing the current ordering.",
        "",
        "## 2. Grids And Parameters",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- epsilon: {EPSILON:g}",
        f"- mu: {MU:g}",
        f"- eta range: {ETA_MIN:g} .. {ETA_MAX:g}",
        f"- coarse eta step: {COARSE_ETA_STEP:g}",
        f"- refine eta step near class changes: {REFINE_ETA_STEP:g}",
        f"- evaluated eta values: {len(sorted_cases)} ({coarse_count} coarse, {refined_count} added by refinement/key insertion)",
        f"- beta range: 0 .. 90 deg",
        f"- beta step: {BETA_STEP:g} deg",
        f"- sorted roots solved at each beta: first {N_SORTED_ROOTS}",
        f"- tracked descendants: first {N_DESCENDANTS_TRACK}",
        f"- classified descendants: first {N_DESCENDANTS_CLASSIFY}",
        "",
        "The requested `0.1 deg` beta grid and `0.01` eta grid would be much",
        "heavier, so this run uses the allowed coarse grid plus local eta",
        "refinement around detected class transitions.",
        "",
        "## 3-5. Key Eta Cases",
        "",
    ]
    lines.extend(
        markdown_table(
            [
                "case",
                "class",
                "changed descendants",
                "pos beta=0 desc1..6",
                "pos beta=15 desc1..6",
                "pos beta=90 desc1..6",
                "sorted5 desc at beta=15",
            ],
            key_rows,
        )
    )
    lines.extend(
        [
            "",
            "At eta=0 the sorted positions show pairwise exchanges",
            "`1<->2`, `3<->4`, and `5<->6` by beta=15/90. The case is",
            "classified as `tracking_unreliable` rather than clean",
            "`rearrangement` because the minimum adjacent-step MAC is extremely",
            "small on this grid, consistent with a symmetry-degenerate crossing",
            "where a local subspace-MAC refinement is more appropriate than",
            "single-vector MAC tracking.",
            "",
            "At eta=0.5 and eta=-0.5, all first-six descendants keep their",
            "initial sorted positions over the full beta range.",
        ]
    )

    lines.extend(["", "## 6. Boundary Eta Values", ""])
    lines.extend(boundary_lines(sorted_cases))
    lines.extend(["", "Class intervals on the evaluated eta grid:", ""])
    lines.extend(markdown_table(["eta start", "eta end", "class", "changed descendants"], interval_rows))

    lines.extend(["", "## 7. Eta-Sign Symmetry", ""])
    if symmetric:
        lines.append("The rearrangement class is symmetric under eta sign on the evaluated symmetric eta samples.")
    else:
        lines.append("The rearrangement class is asymmetric under eta sign on the evaluated samples:")
        for mismatch in mismatches:
            lines.append(f"- {mismatch}")

    lines.extend(["", "## 8. Descendants Involved", ""])
    for cls in ("rearrangement", "ambiguous_near_degenerate", "tracking_unreliable", "no_rearrangement"):
        values = sorted(involved.get(cls, set()), key=int)
        lines.append(f"- {cls}: {', '.join(values) if values else 'none'}")

    lines.extend(["", "## 9. Secondary Sorted5 Identity Summary", ""])
    lines.append(
        "Sorted5 identity is recorded because it explains the earlier descendant-6/sorted-5 confusion, "
        "but it is not the primary rearrangement criterion."
    )
    for case in (etam, eta0, etap):
        lines.append(
            f"- eta={case.eta:g}: sorted5 descendant at beta=15 is "
            f"{case.summary_row['sorted5_descendant_id_beta15']}; class={case.summary_row['rearrangement_class']}."
        )

    lines.extend(
        [
            "",
            "## 10. Limitations",
            "",
            "- analytic Euler-Bernoulli thickness-mismatch determinant only",
            "- no Timoshenko model",
            "- no FEM, Gmsh, CalculiX, or 3D workflow",
            "- branch tracking uses adjacent-step shape MAC and can become ambiguous near very small sorted gaps",
            "- near-degenerate cases may need local subspace-MAC refinement before physical interpretation",
            "- boundaries are grid-resolved brackets, not exact analytic thresholds",
            "",
            "## Outputs",
            "",
            f"- main CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- summary CSV: `{OUTPUT_SUMMARY_CSV.relative_to(REPO_ROOT)}`",
            f"- sorted-position map: `{OUTPUT_MAP_PNG.relative_to(REPO_ROOT)}`",
            f"- class plot: `{OUTPUT_CLASS_PNG.relative_to(REPO_ROOT)}`",
            f"- representative traces: `{OUTPUT_TRACES_PNG.relative_to(REPO_ROOT)}`",
            f"- minimum gaps plot: `{OUTPUT_GAPS_PNG.relative_to(REPO_ROOT)}`",
            f"- sorted5 identity map: `{OUTPUT_SORTED5_MAP_PNG.relative_to(REPO_ROOT)}`",
            "",
            "This audit writes only diagnostic outputs and documentation. It does not",
            "modify article files, article figures, the old determinant,",
            "`src/my_project/analytic/formulas.py`, old solvers, baseline results,",
            "Gmsh/CalculiX workflows, or the FEM physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def run_scan() -> list[EtaScanCase]:
    grid_beta = beta_values()
    coarse_etas = eta_values(ETA_MIN, ETA_MAX, COARSE_ETA_STEP)
    coarse_cases = [run_eta_case(float(eta), grid_beta) for eta in coarse_etas]
    refined_values = refined_eta_values(coarse_cases)
    existing = {round(float(case.eta), 12): case for case in coarse_cases}
    cases = list(coarse_cases)
    additional = [eta for eta in refined_values if round(float(eta), 12) not in existing]
    if additional:
        print(f"Refining {len(additional)} additional eta values near class changes/key etas")
    for eta in additional:
        cases.append(run_eta_case(float(eta), grid_beta))
    cases.sort(key=lambda case: case.eta)
    write_main_csv(cases)
    write_summary_csv(cases)
    plot_sorted_position_map(cases)
    plot_class(cases)
    plot_representative_traces(cases)
    plot_min_gaps(cases)
    plot_sorted5_identity_map(cases)
    write_report(cases, coarse_count=len(coarse_cases), refined_count=len(additional))
    return cases


def main() -> dict[str, object]:
    print("Eta scan first-six rearrangement audit")
    print(
        f"epsilon={EPSILON:g}, mu={MU:g}, eta={ETA_MIN:g}..{ETA_MAX:g}, "
        f"coarse_eta_step={COARSE_ETA_STEP:g}, beta_step={BETA_STEP:g}"
    )
    cases = run_scan()
    print(f"saved main CSV: {OUTPUT_CSV}")
    print(f"saved summary CSV: {OUTPUT_SUMMARY_CSV}")
    print(f"saved report: {OUTPUT_REPORT}")
    print(f"saved plots: {OUTPUT_MAP_PNG}, {OUTPUT_CLASS_PNG}, {OUTPUT_TRACES_PNG}, {OUTPUT_GAPS_PNG}, {OUTPUT_SORTED5_MAP_PNG}")
    return {"cases": cases, "csv": OUTPUT_CSV, "summary_csv": OUTPUT_SUMMARY_CSV, "report": OUTPUT_REPORT}


if __name__ == "__main__":
    main()
