from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Mapping, Sequence
import argparse

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
ETA = 0.0

MU_MIN = 0.0
MU_MAX = 0.9
COARSE_MU_STEP = 0.01
REFINE_MU_STEP = 0.001
BETA_STEP = 0.25

BETA_MIN_CHECK = 1.0
ROBUST_BETA_MIN_CHECKS = (0.5, 1.0, 2.0)
KEY_BETAS = (0.0, 1.0, 5.0, 15.0, 90.0)
KEY_MUS = (0.0, 0.9)

N_DESCENDANTS_TRACK = 8
N_DESCENDANTS_CLASSIFY = 6
N_SORTED_ROOTS = 14
N_TRACKING_CANDIDATE_ROOTS = 10
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
FINITE_INTERVAL_MIN_WIDTH_DEG = 0.5

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_first6_rearrangement.csv"
OUTPUT_SUMMARY_CSV = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_first6_rearrangement_summary.csv"
OUTPUT_REPORT = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_first6_rearrangement_report.md"
OUTPUT_MAP_PNG = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_first6_sorted_position_map.png"
OUTPUT_CLASS_PNG = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_first6_rearrangement_class.png"
OUTPUT_GAPS_PNG = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_min_adjacent_gaps_first7.png"
OUTPUT_TRACES_PNG = OUTPUT_DIR / "mu_scan_eta0_eps0p0025_sorted_position_traces_representative_mu.png"

SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

CLASS_CODES = {
    "no_rearrangement": 0,
    "beta0_seed_degeneracy_or_label_split": 1,
    "true_beta_rearrangement": 2,
    "ambiguous_near_degenerate": 3,
    "tracking_unreliable": 4,
}

MAIN_FIELDNAMES = [
    "mu",
    "eta",
    "epsilon",
    "beta_deg",
    "descendant_id",
    "Lambda_descendant",
    "sorted_position",
    "mac_to_previous",
    "tracking_warning",
    "is_position_changed_relative_to_beta0",
    "is_position_changed_relative_to_beta_min_check",
    "sorted5_descendant_id_at_beta",
    "sorted5_lambda",
    "gap_sorted1_2",
    "gap_sorted2_3",
    "gap_sorted3_4",
    "gap_sorted4_5",
    "gap_sorted5_6",
    "gap_sorted6_7",
    "min_adjacent_gap_first7_at_beta",
    "is_within_first6",
]

SUMMARY_FIELDNAMES = [
    "mu",
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
    "desc1_pos_beta1",
    "desc2_pos_beta1",
    "desc3_pos_beta1",
    "desc4_pos_beta1",
    "desc5_pos_beta1",
    "desc6_pos_beta1",
    "desc1_pos_beta5",
    "desc2_pos_beta5",
    "desc3_pos_beta5",
    "desc4_pos_beta5",
    "desc5_pos_beta5",
    "desc6_pos_beta5",
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
    "sorted5_identity_over_beta",
    "sorted5_identity_changes_over_beta",
    "desc5_desc6_swaps_over_beta",
    "min_gap_first7",
    "beta_min_gap_first7",
    "min_gap_sorted1_2",
    "beta_min_gap_sorted1_2",
    "min_gap_sorted2_3",
    "beta_min_gap_sorted2_3",
    "min_gap_sorted3_4",
    "beta_min_gap_sorted3_4",
    "min_gap_sorted4_5",
    "beta_min_gap_sorted4_5",
    "min_gap_sorted5_6",
    "beta_min_gap_sorted5_6",
    "min_gap_sorted6_7",
    "beta_min_gap_sorted6_7",
    "tracking_warning_count",
    "min_mac_continuity",
    "rearrangement_class_beta_min_0p5",
    "descendants_changed_beta_min_0p5",
    "rearrangement_class_beta_min_2",
    "descendants_changed_beta_min_2",
    "notes",
]


@dataclass(frozen=True)
class Classification:
    name: str
    changed_descendants: list[int]
    note: str


@dataclass(frozen=True)
class MuScanCase:
    mu: float
    result: beta_audit.BetaTrackingResult
    classifications: dict[float, Classification]
    summary_row: dict[str, str]


def number_text(value: float) -> str:
    return f"{float(value):.12g}"


def beta_min_token(value: float) -> str:
    text = f"{float(value):g}".replace(".", "p")
    return text


def yesno(value: bool) -> str:
    return "yes" if bool(value) else "no"


def mu_values(start: float, stop: float, step: float) -> np.ndarray:
    values = np.arange(float(start), float(stop) + 0.5 * float(step), float(step), dtype=float)
    values[0] = float(start)
    values[-1] = float(stop)
    return np.unique(np.round(values, 12))


def beta_values() -> np.ndarray:
    values = np.arange(0.0, 90.0 + 0.5 * BETA_STEP, BETA_STEP, dtype=float)
    values[0] = 0.0
    values[-1] = 90.0
    return np.unique(np.round(values, 12))


def configure_beta_backend(mu: float, beta_grid: np.ndarray) -> None:
    beta_audit.EPSILON = EPSILON
    beta_audit.MU = float(mu)
    beta_audit.BETA_VALUES_DEG = np.asarray(beta_grid, dtype=float)
    beta_audit.N_DESCENDANTS = N_DESCENDANTS_TRACK
    beta_audit.N_SORTED_ROOTS = N_SORTED_ROOTS
    beta_audit.ROOT_SCAN_STEP = ROOT_SCAN_STEP
    beta_audit.ROOT_LMAX0 = ROOT_LMAX0
    beta_audit.NUM_SHAPE_SAMPLES = NUM_SHAPE_SAMPLES
    beta_audit.MAC_WARNING_THRESHOLD = MAC_WARNING_THRESHOLD
    beta_audit.MAC_MARGIN_WARNING_THRESHOLD = MAC_MARGIN_WARNING_THRESHOLD
    beta_audit.MAX_SORTED_POSITION_JUMP = MAX_SORTED_POSITION_JUMP
    beta_audit.shape_vectors_for_beta = limited_shape_vectors_for_beta


def limited_shape_vectors_for_beta(
    sorted_roots: np.ndarray,
    *,
    eta: float,
    beta_deg: float,
) -> list[np.ndarray]:
    roots = np.asarray(sorted_roots, dtype=float)[:N_TRACKING_CANDIDATE_ROOTS]
    if len(roots) < N_DESCENDANTS_TRACK:
        raise RuntimeError(
            f"Need at least {N_DESCENDANTS_TRACK} roots for tracking candidates; got {len(roots)}."
        )
    return beta_audit.analytic_shape_vectors_for_roots(
        roots,
        beta_rad=float(np.deg2rad(beta_deg)),
        mu=float(beta_audit.MU),
        epsilon=float(beta_audit.EPSILON),
        eta=float(eta),
        s_norm=np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float),
    )


def finite_min(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(np.min(finite)) if finite.size else float("nan")


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


def first_beta_at_or_after(beta_grid: np.ndarray, beta_deg: float) -> int:
    matches = np.flatnonzero(beta_grid >= float(beta_deg) - 1e-12)
    if matches.size:
        return int(matches[0])
    return len(beta_grid) - 1


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


def compact_string_sequence(beta_grid: np.ndarray, values: Sequence[str], *, label: str) -> str:
    parts: list[str] = []
    start = 0
    values_list = [str(value) for value in values]
    for idx in range(1, len(values_list) + 1):
        if idx < len(values_list) and values_list[idx] == values_list[start]:
            continue
        beta_start = float(beta_grid[start])
        beta_end = float(beta_grid[idx - 1])
        if np.isclose(beta_start, beta_end, rtol=0.0, atol=1e-12):
            parts.append(f"beta {beta_start:g} -> {label} {values_list[start]}")
        else:
            parts.append(f"beta {beta_start:g}..{beta_end:g} -> {label} {values_list[start]}")
        start = idx
    return "; ".join(parts)


def adjacent_gaps(sorted_roots: np.ndarray) -> np.ndarray:
    return np.diff(sorted_roots[:, :N_SORTED_GAP_ROOTS], axis=1)


def sorted5_lambda_at_beta(sorted_roots_at_beta: np.ndarray) -> float:
    return float(sorted_roots_at_beta[4])


def change_event_indices(positions: np.ndarray, descendant_zero_based: int, *, start_col: int = 0) -> np.ndarray:
    changes = np.flatnonzero(np.diff(positions[descendant_zero_based]) != 0) + 1
    return changes[changes >= int(start_col)]


def run_width_for_position(beta_grid: np.ndarray, positions: np.ndarray, col: int) -> float:
    value = int(positions[col])
    start = int(col)
    while start > 0 and int(positions[start - 1]) == value:
        start -= 1
    end = int(col)
    while end + 1 < len(positions) and int(positions[end + 1]) == value:
        end += 1
    return float(beta_grid[end] - beta_grid[start])


def classify_case(result: beta_audit.BetaTrackingResult, *, beta_min_check: float) -> Classification:
    beta_grid = np.asarray(result.beta_values_deg, dtype=float)
    positions = result.current_sorted_positions[:N_DESCENDANTS_CLASSIFY]
    start_col = first_beta_at_or_after(beta_grid, float(beta_min_check))
    active_positions = positions[:, start_col:]
    active_changed_descendants = [
        desc_idx + 1
        for desc_idx in range(N_DESCENDANTS_CLASSIFY)
        if np.unique(active_positions[desc_idx]).size > 1
    ]
    all_changed_descendants = [
        desc_idx + 1
        for desc_idx in range(N_DESCENDANTS_CLASSIFY)
        if np.unique(positions[desc_idx]).size > 1
    ]

    if not active_changed_descendants:
        if all_changed_descendants:
            return Classification(
                "beta0_seed_degeneracy_or_label_split",
                all_changed_descendants,
                f"changes occur before beta_min_check={float(beta_min_check):g} deg only",
            )
        return Classification(
            "no_rearrangement",
            [],
            f"all first-six descendants keep constant sorted positions for beta >= {float(beta_min_check):g} deg",
        )

    active_rows = [
        row for row in result.rows if float(row["beta_deg"]) >= float(beta_min_check) - 1e-12
    ]
    finite_macs = [
        float(row["mac_to_previous"])
        for row in active_rows
        if np.isfinite(float(row["mac_to_previous"]))
    ]
    min_mac = finite_min(finite_macs)
    active_warning_rows = [row for row in active_rows if tracking_warning(row) != "no"]
    warning_fraction = len(active_warning_rows) / max(1, len(active_rows))

    if np.isfinite(min_mac) and min_mac < TRACKING_UNRELIABLE_MAC:
        return Classification(
            "tracking_unreliable",
            active_changed_descendants,
            f"minimum active-beta MAC {min_mac:.6g} below {TRACKING_UNRELIABLE_MAC:g}",
        )
    if warning_fraction > TRACKING_UNRELIABLE_WARNING_FRACTION:
        return Classification(
            "tracking_unreliable",
            active_changed_descendants,
            f"active-beta warning fraction {warning_fraction:.6g} above {TRACKING_UNRELIABLE_WARNING_FRACTION:g}",
        )

    row_lookup = rows_by_beta_descendant(result)
    gaps = adjacent_gaps(result.sorted_roots)
    suspicious_flags: list[bool] = []
    suspicious_reasons: list[str] = []
    for desc_id in active_changed_descendants:
        desc_positions = positions[desc_id - 1]
        for col in change_event_indices(positions, desc_id - 1, start_col=start_col):
            beta_deg = float(beta_grid[col])
            row = row_lookup[(round(beta_deg, 10), desc_id)]
            local_gap = float(np.nanmin(gaps[col]))
            mac = float(row["mac_to_previous"])
            run_width = run_width_for_position(beta_grid, desc_positions, col)
            warning = tracking_warning(row)
            suspicious = (
                warning != "no"
                or (np.isfinite(mac) and mac < MAC_WARNING_THRESHOLD)
                or local_gap < NEAR_DEGENERATE_GAP
                or run_width < FINITE_INTERVAL_MIN_WIDTH_DEG
            )
            suspicious_flags.append(bool(suspicious))
            if suspicious:
                reason_bits = []
                if warning != "no":
                    reason_bits.append(warning)
                if np.isfinite(mac) and mac < MAC_WARNING_THRESHOLD:
                    reason_bits.append(f"MAC {mac:.4g}")
                if local_gap < NEAR_DEGENERATE_GAP:
                    reason_bits.append(f"gap {local_gap:.4g}")
                if run_width < FINITE_INTERVAL_MIN_WIDTH_DEG:
                    reason_bits.append(f"run width {run_width:.4g} deg")
                suspicious_reasons.append(
                    f"desc {desc_id} beta {beta_deg:g}: " + ", ".join(reason_bits)
                )

    if suspicious_flags and all(suspicious_flags):
        return Classification(
            "ambiguous_near_degenerate",
            active_changed_descendants,
            "all active sorted-position changes are local warning/near-degenerate events: "
            + "; ".join(suspicious_reasons[:4]),
        )

    return Classification(
        "true_beta_rearrangement",
        active_changed_descendants,
        f"accepted finite-beta sorted-position changes for beta >= {float(beta_min_check):g} deg",
    )


def classification_for(result: beta_audit.BetaTrackingResult) -> dict[float, Classification]:
    return {
        float(beta_min): classify_case(result, beta_min_check=float(beta_min))
        for beta_min in ROBUST_BETA_MIN_CHECKS
    }


def sorted5_changes_text(beta_grid: np.ndarray, identities: Sequence[str]) -> str:
    if len(set(str(value) for value in identities)) <= 1:
        return "no"
    return "yes: " + compact_string_sequence(beta_grid, identities, label="desc")


def desc5_desc6_swap_text(beta_grid: np.ndarray, positions: np.ndarray) -> str:
    desc5 = positions[4]
    desc6 = positions[5]
    swap_mask = (desc5 == 6) & (desc6 == 5)
    if not np.any(swap_mask):
        return "no"
    values = ["swap" if bool(flag) else "no_swap" for flag in swap_mask]
    return "yes: " + compact_string_sequence(beta_grid, values, label="state")


def min_gap_pair_columns(beta_grid: np.ndarray, gaps: np.ndarray) -> dict[str, str]:
    row: dict[str, str] = {}
    for pair_idx in range(gaps.shape[1]):
        gap_values = np.asarray(gaps[:, pair_idx], dtype=float)
        min_idx = int(np.nanargmin(gap_values))
        pair_label = f"sorted{pair_idx + 1}_{pair_idx + 2}"
        row[f"min_gap_{pair_label}"] = number_text(float(gap_values[min_idx]))
        row[f"beta_min_gap_{pair_label}"] = number_text(float(beta_grid[min_idx]))
    return row


def summary_for_result(
    result: beta_audit.BetaTrackingResult,
    classifications: Mapping[float, Classification],
) -> dict[str, str]:
    primary = classifications[float(BETA_MIN_CHECK)]
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
        "mu": number_text(float(beta_audit.MU)),
        "rearrangement_class": primary.name,
        "rearrangement_detected_bool": yesno(primary.name == "true_beta_rearrangement"),
        "descendants_with_position_changes": ";".join(str(value) for value in primary.changed_descendants) or "none",
        "num_descendants_changed": str(len(primary.changed_descendants)),
        "sorted5_descendant_id_beta15": sorted5_identities[key_beta_index(beta_grid, 15.0)],
        "sorted5_identity_over_beta": compact_string_sequence(beta_grid, sorted5_identities, label="desc"),
        "sorted5_identity_changes_over_beta": sorted5_changes_text(beta_grid, sorted5_identities),
        "desc5_desc6_swaps_over_beta": desc5_desc6_swap_text(beta_grid, positions),
        "min_gap_first7": number_text(float(min_gap_by_beta[min_gap_idx])),
        "beta_min_gap_first7": number_text(float(beta_grid[min_gap_idx])),
        "tracking_warning_count": str(warning_count),
        "min_mac_continuity": number_text(finite_min(finite_macs)),
        "notes": primary.note,
    }
    row.update(min_gap_pair_columns(beta_grid, gap_matrix))

    for desc_id in range(1, N_DESCENDANTS_CLASSIFY + 1):
        row[f"desc{desc_id}_positions_over_beta"] = compact_position_sequence(beta_grid, positions[desc_id - 1])

    beta_labels = (
        ("beta0", 0.0),
        ("beta1", 1.0),
        ("beta5", 5.0),
        ("beta15", 15.0),
        ("beta90", 90.0),
    )
    for beta_label, beta_value in beta_labels:
        idx = key_beta_index(beta_grid, beta_value)
        for desc_id in range(1, N_DESCENDANTS_CLASSIFY + 1):
            row[f"desc{desc_id}_pos_{beta_label}"] = str(int(positions[desc_id - 1, idx]))

    for beta_min in (0.5, 2.0):
        cls = classifications[float(beta_min)]
        token = beta_min_token(beta_min)
        row[f"rearrangement_class_beta_min_{token}"] = cls.name
        row[f"descendants_changed_beta_min_{token}"] = ";".join(str(value) for value in cls.changed_descendants) or "none"

    return row


def run_mu_case(mu: float, beta_grid: np.ndarray) -> MuScanCase:
    configure_beta_backend(float(mu), beta_grid)
    roots = beta_audit.solve_sorted_root_grid(eta=ETA)
    result = beta_audit.track_descendants_by_beta(eta=ETA, sorted_roots=roots)
    classifications = classification_for(result)
    summary = summary_for_result(result, classifications)
    print(
        f"mu={float(mu): .6f}: class={summary['rearrangement_class']}, "
        f"changed={summary['descendants_with_position_changes']}, "
        f"warnings={summary['tracking_warning_count']}, min_mac={summary['min_mac_continuity']}"
    )
    return MuScanCase(mu=float(mu), result=result, classifications=classifications, summary_row=summary)


def class_change_intervals(cases: Sequence[MuScanCase]) -> list[tuple[float, float]]:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    intervals: list[tuple[float, float]] = []
    for left, right in zip(sorted_cases[:-1], sorted_cases[1:], strict=True):
        left_class = left.summary_row["rearrangement_class"]
        right_class = right.summary_row["rearrangement_class"]
        if left_class != right_class:
            intervals.append((float(left.mu), float(right.mu)))
    return intervals


def refined_mu_values(coarse_cases: Sequence[MuScanCase]) -> np.ndarray:
    values: set[float] = set(float(case.mu) for case in coarse_cases)
    for left, right in class_change_intervals(coarse_cases):
        lo, hi = sorted((left, right))
        for value in mu_values(lo, hi, REFINE_MU_STEP):
            values.add(float(value))
    for mu in KEY_MUS:
        values.add(float(mu))
    return np.asarray(sorted(values), dtype=float)


def write_main_csv(cases: Sequence[MuScanCase]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=MAIN_FIELDNAMES)
        writer.writeheader()
        for case in sorted(cases, key=lambda item: item.mu):
            result = case.result
            row_lookup = rows_by_beta_descendant(result)
            sorted5_identities = sorted5_identity_by_beta(result)
            gaps = adjacent_gaps(result.sorted_roots)
            initial_positions = result.current_sorted_positions[:, 0].copy()
            beta_min_idx = first_beta_at_or_after(np.asarray(result.beta_values_deg, dtype=float), BETA_MIN_CHECK)
            beta_min_positions = result.current_sorted_positions[:, beta_min_idx].copy()
            for beta_idx, beta_deg in enumerate(result.beta_values_deg):
                beta_gaps = gaps[beta_idx]
                for desc_idx in range(N_DESCENDANTS_TRACK):
                    row = row_lookup[(round(float(beta_deg), 10), desc_idx + 1)]
                    mac = float(row["mac_to_previous"])
                    writer.writerow(
                        {
                            "mu": number_text(case.mu),
                            "eta": number_text(ETA),
                            "epsilon": number_text(EPSILON),
                            "beta_deg": number_text(float(beta_deg)),
                            "descendant_id": desc_idx + 1,
                            "Lambda_descendant": number_text(float(result.tracked_lambdas[desc_idx, beta_idx])),
                            "sorted_position": int(result.current_sorted_positions[desc_idx, beta_idx]),
                            "mac_to_previous": "" if not np.isfinite(mac) else number_text(mac),
                            "tracking_warning": tracking_warning(row),
                            "is_position_changed_relative_to_beta0": yesno(
                                int(result.current_sorted_positions[desc_idx, beta_idx]) != int(initial_positions[desc_idx])
                            ),
                            "is_position_changed_relative_to_beta_min_check": yesno(
                                int(result.current_sorted_positions[desc_idx, beta_idx])
                                != int(beta_min_positions[desc_idx])
                            ),
                            "sorted5_descendant_id_at_beta": sorted5_identities[beta_idx],
                            "sorted5_lambda": number_text(sorted5_lambda_at_beta(result.sorted_roots[beta_idx])),
                            "gap_sorted1_2": number_text(float(beta_gaps[0])),
                            "gap_sorted2_3": number_text(float(beta_gaps[1])),
                            "gap_sorted3_4": number_text(float(beta_gaps[2])),
                            "gap_sorted4_5": number_text(float(beta_gaps[3])),
                            "gap_sorted5_6": number_text(float(beta_gaps[4])),
                            "gap_sorted6_7": number_text(float(beta_gaps[5])),
                            "min_adjacent_gap_first7_at_beta": number_text(float(np.nanmin(beta_gaps))),
                            "is_within_first6": yesno(desc_idx < N_DESCENDANTS_CLASSIFY),
                        }
                    )


def write_summary_csv(cases: Sequence[MuScanCase]) -> None:
    OUTPUT_SUMMARY_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_SUMMARY_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        for case in sorted(cases, key=lambda item: item.mu):
            writer.writerow(case.summary_row)


def class_intervals(cases: Sequence[MuScanCase]) -> list[dict[str, str]]:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    if not sorted_cases:
        return []
    intervals: list[dict[str, str]] = []
    start = 0
    for idx in range(1, len(sorted_cases) + 1):
        if idx < len(sorted_cases) and (
            sorted_cases[idx].summary_row["rearrangement_class"]
            == sorted_cases[start].summary_row["rearrangement_class"]
            and sorted_cases[idx].summary_row["descendants_with_position_changes"]
            == sorted_cases[start].summary_row["descendants_with_position_changes"]
        ):
            continue
        interval_cases = sorted_cases[start:idx]
        changed: set[str] = set()
        for case in interval_cases:
            value = case.summary_row["descendants_with_position_changes"]
            if value != "none":
                changed.update(value.split(";"))
        intervals.append(
            {
                "mu_start": number_text(float(interval_cases[0].mu)),
                "mu_end": number_text(float(interval_cases[-1].mu)),
                "class": interval_cases[0].summary_row["rearrangement_class"],
                "changed_descendants": ",".join(sorted(changed, key=int)) if changed else "none",
            }
        )
        start = idx
    return intervals


def boundary_lines(cases: Sequence[MuScanCase]) -> list[str]:
    lines: list[str] = []
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    for left, right in zip(sorted_cases[:-1], sorted_cases[1:], strict=True):
        left_class = left.summary_row["rearrangement_class"]
        right_class = right.summary_row["rearrangement_class"]
        if left_class == right_class:
            continue
        mid = 0.5 * (float(left.mu) + float(right.mu))
        lines.append(
            f"- mu in [{float(left.mu):.12g}, {float(right.mu):.12g}] changes "
            f"`{left_class}` -> `{right_class}`; midpoint about `{mid:.12g}`"
        )
    return lines or ["- no class changes were detected on the evaluated mu grid"]


def case_nearest(cases: Sequence[MuScanCase], mu: float) -> MuScanCase:
    return min(cases, key=lambda case: abs(float(case.mu) - float(mu)))


def first_positive_case(cases: Sequence[MuScanCase]) -> MuScanCase:
    positives = [case for case in cases if float(case.mu) > 0.0]
    return min(positives, key=lambda case: float(case.mu)) if positives else case_nearest(cases, 0.0)


def positions_at_key_betas(case: MuScanCase) -> dict[float, list[int]]:
    beta_grid = np.asarray(case.result.beta_values_deg, dtype=float)
    out: dict[float, list[int]] = {}
    for beta in KEY_BETAS:
        idx = key_beta_index(beta_grid, beta)
        out[beta] = [int(case.result.current_sorted_positions[desc - 1, idx]) for desc in range(1, 7)]
    return out


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def plot_sorted_position_map(cases: Sequence[MuScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    fig, axes = plt.subplots(3, 2, figsize=(11.0, 9.5), sharex=True, sharey=True)
    cmap = plt.get_cmap("tab10", N_DESCENDANTS_TRACK + 1)
    last_scatter = None
    for desc_idx, ax in enumerate(axes.flat):
        for case in sorted_cases:
            beta = case.result.beta_values_deg
            mu_line = np.full_like(beta, float(case.mu), dtype=float)
            values = case.result.current_sorted_positions[desc_idx]
            last_scatter = ax.scatter(
                beta,
                mu_line,
                c=values,
                cmap=cmap,
                vmin=0.5,
                vmax=N_DESCENDANTS_TRACK + 0.5,
                s=7,
                linewidths=0,
            )
        ax.set_title(f"descendant {desc_idx + 1}")
        ax.grid(True, color="0.88", linewidth=0.5)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\beta$ (deg)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\mu$")
    if last_scatter is not None:
        cbar = fig.colorbar(last_scatter, ax=axes.ravel().tolist(), ticks=range(1, N_DESCENDANTS_TRACK + 1), shrink=0.92)
        cbar.set_label("sorted position")
    fig.suptitle(r"First-six descendant sorted positions over $\beta$ and $\mu$ at $\eta=0$")
    fig.savefig(OUTPUT_MAP_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_class(cases: Sequence[MuScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    fig, ax = plt.subplots(figsize=(10.0, 3.8))
    x = [float(case.mu) for case in sorted_cases]
    y = [CLASS_CODES[case.summary_row["rearrangement_class"]] for case in sorted_cases]
    ax.step(x, y, where="mid", color="#226f54", linewidth=1.6)
    ax.scatter(x, y, s=16, color="#d1495b", zorder=3)
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel("class")
    ax.set_yticks(list(CLASS_CODES.values()))
    ax.set_yticklabels(list(CLASS_CODES.keys()))
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.set_title(r"First-six rearrangement class over $\mu$ at $\eta=0$")
    fig.tight_layout()
    fig.savefig(OUTPUT_CLASS_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_min_gaps(cases: Sequence[MuScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    fig, ax = plt.subplots(figsize=(10.2, 5.2))
    mu = np.array([float(case.mu) for case in sorted_cases], dtype=float)
    pair_colors = ["#2a6fbb", "#b6465f", "#378b5b", "#8661c1", "#c8792a", "#497174"]
    for pair_idx, color in enumerate(pair_colors):
        values = []
        for case in sorted_cases:
            gaps = adjacent_gaps(case.result.sorted_roots)
            values.append(float(np.nanmin(gaps[:, pair_idx])))
        ax.plot(mu, values, marker="o", markersize=2.8, linewidth=1.0, color=color, label=f"sorted {pair_idx + 1}-{pair_idx + 2}")
    global_values = [
        float(np.nanmin(adjacent_gaps(case.result.sorted_roots)))
        for case in sorted_cases
    ]
    ax.plot(mu, global_values, color="black", linewidth=1.4, label="global min 1..7")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel("minimum adjacent gap over beta")
    ax.set_yscale("log")
    ax.grid(True, color="0.88", linewidth=0.6, which="both")
    ax.legend(loc="best", ncol=2, fontsize=8, frameon=False)
    ax.set_title("Minimum adjacent sorted gaps among roots 1..7")
    fig.tight_layout()
    fig.savefig(OUTPUT_GAPS_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def representative_mus(cases: Sequence[MuScanCase]) -> list[float]:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    reps = {0.0, 0.9}
    intervals = class_intervals(sorted_cases)
    for interval in intervals:
        mu_start = float(interval["mu_start"])
        mu_end = float(interval["mu_end"])
        reps.add(0.5 * (mu_start + mu_end))
    for left, right in class_change_intervals(sorted_cases):
        reps.add(float(left))
        reps.add(float(right))
    nearest = sorted({case_nearest(sorted_cases, value).mu for value in reps})
    if len(nearest) > 10:
        keep = {nearest[0], nearest[-1]}
        for value in nearest[1:-1]:
            if len(keep) >= 10:
                break
            keep.add(value)
        nearest = sorted(keep)
    return nearest


def plot_representative_traces(cases: Sequence[MuScanCase]) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    reps = representative_mus(sorted_cases)
    ncols = 2
    nrows = int(np.ceil(len(reps) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(11.0, max(3.0, 2.55 * nrows)), sharex=True, sharey=True)
    axes_array = np.atleast_1d(axes).ravel()
    colors = ["#2a6fbb", "#b6465f", "#378b5b", "#8661c1", "#c8792a", "#497174"]
    for ax, mu in zip(axes_array, reps, strict=False):
        case = case_nearest(sorted_cases, mu)
        beta = case.result.beta_values_deg
        for desc_idx in range(N_DESCENDANTS_CLASSIFY):
            ax.step(
                beta,
                case.result.current_sorted_positions[desc_idx],
                where="post",
                color=colors[desc_idx],
                linewidth=1.2,
                label=f"desc {desc_idx + 1}" if ax is axes_array[0] else None,
            )
        ax.set_title(f"mu={float(case.mu):g}; {case.summary_row['rearrangement_class']}")
        ax.grid(True, color="0.88", linewidth=0.5)
    for ax in axes_array[len(reps):]:
        ax.axis("off")
    for ax in axes_array:
        ax.set_xlabel(r"$\beta$ (deg)")
        ax.set_ylabel("sorted position")
        ax.set_yticks(range(1, N_DESCENDANTS_TRACK + 1))
    axes_array[0].legend(loc="best", ncol=3, fontsize=8, frameon=False)
    fig.suptitle(r"Representative first-six sorted-position traces at $\eta=0$")
    fig.tight_layout()
    fig.savefig(OUTPUT_TRACES_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def selected_case_rows(cases: Sequence[MuScanCase]) -> list[list[str]]:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    selected: list[tuple[str, MuScanCase]] = [
        ("mu=0", case_nearest(sorted_cases, 0.0)),
        ("small positive mu", first_positive_case(sorted_cases)),
        ("mu=0.9", case_nearest(sorted_cases, 0.9)),
    ]
    for left, right in class_change_intervals(sorted_cases):
        selected.append((f"boundary left {float(left):g}", case_nearest(sorted_cases, left)))
        selected.append((f"boundary right {float(right):g}", case_nearest(sorted_cases, right)))

    seen: set[tuple[str, float]] = set()
    rows: list[list[str]] = []
    for label, case in selected:
        key = (label, round(float(case.mu), 12))
        if key in seen:
            continue
        seen.add(key)
        positions = positions_at_key_betas(case)
        rows.append(
            [
                label,
                number_text(case.mu),
                case.summary_row["rearrangement_class"],
                case.summary_row["descendants_with_position_changes"],
                ",".join(str(value) for value in positions[0.0]),
                ",".join(str(value) for value in positions[1.0]),
                ",".join(str(value) for value in positions[5.0]),
                ",".join(str(value) for value in positions[15.0]),
                ",".join(str(value) for value in positions[90.0]),
                case.summary_row["sorted5_descendant_id_beta15"],
                case.summary_row["min_gap_first7"],
                case.summary_row["tracking_warning_count"],
                case.summary_row["min_mac_continuity"],
            ]
        )
    return rows


def write_report(cases: Sequence[MuScanCase], *, coarse_count: int, refined_count: int) -> None:
    sorted_cases = sorted(cases, key=lambda case: case.mu)
    intervals = class_intervals(sorted_cases)
    interval_rows = [
        [row["mu_start"], row["mu_end"], row["class"], row["changed_descendants"]]
        for row in intervals
    ]

    involved: dict[str, set[str]] = {}
    for case in sorted_cases:
        cls = case.summary_row["rearrangement_class"]
        involved.setdefault(cls, set())
        value = case.summary_row["descendants_with_position_changes"]
        if value != "none":
            involved[cls].update(value.split(";"))

    sorted5_beta15_counts: dict[str, int] = {}
    sorted5_change_count = 0
    desc5_desc6_swap_count = 0
    for case in sorted_cases:
        sorted5_beta15_counts[case.summary_row["sorted5_descendant_id_beta15"]] = (
            sorted5_beta15_counts.get(case.summary_row["sorted5_descendant_id_beta15"], 0) + 1
        )
        if case.summary_row["sorted5_identity_changes_over_beta"].startswith("yes"):
            sorted5_change_count += 1
        if case.summary_row["desc5_desc6_swaps_over_beta"].startswith("yes"):
            desc5_desc6_swap_count += 1

    true_intervals = [row for row in intervals if row["class"] == "true_beta_rearrangement"]
    no_intervals = [row for row in intervals if row["class"] == "no_rearrangement"]
    ambiguous_intervals = [
        row for row in intervals if row["class"] in {"ambiguous_near_degenerate", "tracking_unreliable"}
    ]

    lines = [
        "# Mu Scan Eta=0 First-Six Rearrangement Audit",
        "",
        "## 1. What Is Being Checked",
        "",
        "This analytic-only Euler-Bernoulli diagnostic checks whether the first six",
        "descendant branches change sorted position as beta varies from 0 to 90 deg",
        "at fixed `eta=0`, `epsilon=0.0025`, and scanned `mu`.",
        "",
        "The calculation is diagnostic-only. It does not generate article figures",
        "and does not modify article workspaces, determinant formulas, old solvers,",
        "FEM models, Gmsh/CalculiX workflows, or baseline results.",
        "",
        "## 2. Definition Of Rearrangement",
        "",
        "Primary identity is the descendant branch; sorted position is diagnostic",
        "metadata. For each fixed mu, descendants 1..6 are tracked over beta.",
        f"The primary classification uses `beta_min_check = {BETA_MIN_CHECK:g} deg`",
        "to avoid interpreting an exact beta=0 symmetry seed split as a finite-beta",
        "rearrangement.",
        "",
        "Classes:",
        "",
        "- `no_rearrangement`: first-six sorted positions are constant for beta >= beta_min_check",
        "- `true_beta_rearrangement`: at least one first-six descendant changes sorted position over a finite beta interval after beta_min_check",
        "- `beta0_seed_degeneracy_or_label_split`: changes occur only before beta_min_check",
        "- `ambiguous_near_degenerate`: apparent changes are confined to warning rows, tiny gaps, low MAC, or single-grid-point runs",
        "- `tracking_unreliable`: active-beta tracking has too-low MAC or too many warnings",
        "",
        "## 3. Grids Used",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu range: {MU_MIN:g} .. {MU_MAX:g}",
        f"- coarse mu step: {COARSE_MU_STEP:g}",
        f"- refine mu step near class changes: {REFINE_MU_STEP:g}",
        f"- coarse mu cases: {coarse_count}",
        f"- added refined/key mu cases: {refined_count}",
        f"- total evaluated mu cases: {len(sorted_cases)}",
        f"- beta range: 0 .. 90 deg",
        f"- beta step: {BETA_STEP:g} deg",
        f"- sorted roots solved per beta: first {N_SORTED_ROOTS}",
        f"- descendants tracked: first {N_DESCENDANTS_TRACK}",
        f"- sorted roots used as MAC tracking candidates: first {N_TRACKING_CANDIDATE_ROOTS}",
        f"- descendants classified: first {N_DESCENDANTS_CLASSIFY}",
        f"- shape samples per arm component: {NUM_SHAPE_SAMPLES}",
        f"- beta_min_check robustness: {', '.join(str(value) for value in ROBUST_BETA_MIN_CHECKS)} deg",
        "",
        "## 4. Beta=0 Degeneracy Handling",
        "",
        "The primary class ignores changes that disappear before beta_min_check.",
        "If all first-six descendants are constant for beta >= beta_min_check but",
        "differ only at beta=0 or the immediate near-zero beta samples, the case is",
        "classified as `beta0_seed_degeneracy_or_label_split`. Active-beta low-MAC",
        "or warning-heavy cases are not promoted to true rearrangement.",
        "",
        "## 5. Classification Summary By Mu Intervals",
        "",
    ]
    lines.extend(markdown_table(["mu start", "mu end", "class", "changed descendants"], interval_rows))
    lines.extend(["", "True beta rearrangement intervals:", ""])
    if true_intervals:
        for row in true_intervals:
            lines.append(
                f"- mu in [{row['mu_start']}, {row['mu_end']}]: descendants {row['changed_descendants']}"
            )
    else:
        lines.append("- none on the evaluated grid")
    lines.extend(["", "No rearrangement intervals:", ""])
    if no_intervals:
        for row in no_intervals:
            lines.append(f"- mu in [{row['mu_start']}, {row['mu_end']}]")
    else:
        lines.append("- none on the evaluated grid")
    lines.extend(["", "Ambiguous or tracking-unreliable intervals:", ""])
    if ambiguous_intervals:
        for row in ambiguous_intervals:
            lines.append(
                f"- mu in [{row['mu_start']}, {row['mu_end']}]: {row['class']}, descendants {row['changed_descendants']}"
            )
    else:
        lines.append("- none on the evaluated grid")

    lines.extend(["", "## 6. Approximate Boundary Mu Values", ""])
    lines.extend(boundary_lines(sorted_cases))

    lines.extend(["", "## 7. Descendants Involved", ""])
    for cls in CLASS_CODES:
        values = involved.get(cls, set())
        pretty = ", ".join(sorted(values, key=int)) if values else "none"
        lines.append(f"- {cls}: {pretty}")

    lines.extend(["", "## 8. Key Mu Cases", ""])
    lines.extend(
        markdown_table(
            [
                "case",
                "mu",
                "class",
                "changed",
                "pos beta=0",
                "pos beta=1",
                "pos beta=5",
                "pos beta=15",
                "pos beta=90",
                "sorted5 desc beta=15",
                "min gap",
                "warnings",
                "min MAC",
            ],
            selected_case_rows(sorted_cases),
        )
    )

    lines.extend(
        [
            "",
            "## 9. Sorted5 / Descendant6 Diagnostic",
            "",
            "Sorted5 is a secondary diagnostic; the primary rearrangement criterion is",
            "whether any of descendants 1..6 changes sorted position. Counts for the",
            "descendant occupying sorted position 5 at beta=15 are:",
        ]
    )
    for identity, count in sorted(sorted5_beta15_counts.items(), key=lambda item: item[0]):
        lines.append(f"- sorted5 descendant `{identity}` at beta=15: {count} mu cases")
    lines.extend(
        [
            f"- sorted5 identity changes over beta in {sorted5_change_count} of {len(sorted_cases)} evaluated mu cases",
            f"- desc5/desc6 sorted-position swaps occur in {desc5_desc6_swap_count} of {len(sorted_cases)} evaluated mu cases",
        ]
    )

    lines.extend(["", "## 10. Tracking / Unreliable Regions", ""])
    unreliable_cases = [
        case for case in sorted_cases if case.summary_row["rearrangement_class"] == "tracking_unreliable"
    ]
    ambiguous_cases = [
        case for case in sorted_cases if case.summary_row["rearrangement_class"] == "ambiguous_near_degenerate"
    ]
    if unreliable_cases:
        lines.append(
            "- tracking_unreliable mu values: "
            + ", ".join(number_text(case.mu) for case in unreliable_cases)
        )
    else:
        lines.append("- tracking_unreliable mu values: none")
    if ambiguous_cases:
        lines.append(
            "- ambiguous_near_degenerate mu values: "
            + ", ".join(number_text(case.mu) for case in ambiguous_cases)
        )
    else:
        lines.append("- ambiguous_near_degenerate mu values: none")
    lines.append("")
    lines.append("Representative notes for non-clean cases:")
    for case in (unreliable_cases + ambiguous_cases)[:12]:
        lines.append(
            f"- mu={case.mu:g}: {case.summary_row['rearrangement_class']}; "
            f"{case.summary_row['notes']}"
        )

    lines.extend(
        [
            "",
            "## 11. Limitations",
            "",
            "- EB analytic thickness-mismatch determinant only, evaluated at eta=0",
            "- no Timoshenko model",
            "- no FEM, Gmsh, CalculiX, or heavy 3D workflow",
            "- single-vector MAC tracking near multiple roots may require subspace-MAC refinement",
            "- boundaries are grid-resolved brackets, not exact analytic thresholds",
            "- the beta=0 seed is symmetry-sensitive and is excluded from the primary true-rearrangement decision",
            "",
            "## Outputs",
            "",
            f"- main CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- summary CSV: `{OUTPUT_SUMMARY_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
            f"- sorted-position map: `{OUTPUT_MAP_PNG.relative_to(REPO_ROOT)}`",
            f"- class plot: `{OUTPUT_CLASS_PNG.relative_to(REPO_ROOT)}`",
            f"- minimum gaps plot: `{OUTPUT_GAPS_PNG.relative_to(REPO_ROOT)}`",
            f"- representative traces: `{OUTPUT_TRACES_PNG.relative_to(REPO_ROOT)}`",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def run_scan() -> list[MuScanCase]:
    grid_beta = beta_values()
    coarse_mus = mu_values(MU_MIN, MU_MAX, COARSE_MU_STEP)
    coarse_cases = [run_mu_case(float(mu), grid_beta) for mu in coarse_mus]
    refined_values = refined_mu_values(coarse_cases)
    existing = {round(float(case.mu), 12): case for case in coarse_cases}
    cases = list(coarse_cases)
    additional = [mu for mu in refined_values if round(float(mu), 12) not in existing]
    if additional:
        print(f"Refining {len(additional)} additional mu values near class changes/key mus")
    for mu in additional:
        cases.append(run_mu_case(float(mu), grid_beta))
    cases.sort(key=lambda case: case.mu)
    write_main_csv(cases)
    write_summary_csv(cases)
    plot_sorted_position_map(cases)
    plot_class(cases)
    plot_min_gaps(cases)
    plot_representative_traces(cases)
    write_report(cases, coarse_count=len(coarse_cases), refined_count=len(additional))
    return cases


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run a tiny no-output smoke scan with coarse beta/mu grids.",
    )
    return parser.parse_args(argv)


def run_smoke() -> list[MuScanCase]:
    global BETA_STEP, COARSE_MU_STEP, REFINE_MU_STEP, MU_MIN, MU_MAX
    BETA_STEP = 30.0
    COARSE_MU_STEP = 0.02
    REFINE_MU_STEP = 0.01
    MU_MIN = 0.0
    MU_MAX = 0.02
    grid_beta = beta_values()
    cases = [run_mu_case(float(mu), grid_beta) for mu in mu_values(MU_MIN, MU_MAX, COARSE_MU_STEP)]
    print(
        "smoke summary: "
        + str([(case.mu, case.summary_row["rearrangement_class"]) for case in cases])
    )
    return cases


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    print("Mu scan eta=0 first-six rearrangement audit")
    if args.smoke:
        cases = run_smoke()
        return {"cases": cases}
    print(
        f"epsilon={EPSILON:g}, eta={ETA:g}, mu={MU_MIN:g}..{MU_MAX:g}, "
        f"coarse_mu_step={COARSE_MU_STEP:g}, beta_step={BETA_STEP:g}"
    )
    cases = run_scan()
    print(f"saved main CSV: {OUTPUT_CSV}")
    print(f"saved summary CSV: {OUTPUT_SUMMARY_CSV}")
    print(f"saved report: {OUTPUT_REPORT}")
    print(f"saved plots: {OUTPUT_MAP_PNG}, {OUTPUT_CLASS_PNG}, {OUTPUT_GAPS_PNG}, {OUTPUT_TRACES_PNG}")
    return {"cases": cases, "csv": OUTPUT_CSV, "summary_csv": OUTPUT_SUMMARY_CSV, "report": OUTPUT_REPORT}


if __name__ == "__main__":
    main()
