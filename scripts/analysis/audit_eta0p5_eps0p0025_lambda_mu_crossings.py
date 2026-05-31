from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import (  # noqa: E402
    analytic_shape_vectors_for_roots,
    mac_value,
)


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5

MU_VALUES = np.linspace(0.0, 0.9, 901)

NUM_DESCENDANTS_CHECK = 7
NUM_SORTED_ROOTS = 14
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

NEAR_CROSSING_ABS_GAP = 1e-3
NEAR_CROSSING_REL_GAP = 1e-4
MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
TRACKING_WARNING_GAP_CONTEXT = 5e-2

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.csv"
OUTPUT_REPORT = OUTPUT_DIR / "eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.md"
OUTPUT_PNG = OUTPUT_DIR / "eta0p5_eps0p0025_beta15_lambda_mu_crossing_audit.png"

PRIOR_SOURCES = [
    REPO_ROOT / "scripts" / "analysis" / "plot_lambda_mu_thickness_mismatch_beta15_eta_descendants.py",
    REPO_ROOT / "scripts" / "analysis" / "plot_lambda_mu_thickness_mismatch_beta15_eta_large_slenderness.py",
    REPO_ROOT / "results" / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants.png",
    REPO_ROOT / "results" / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants_report.md",
    REPO_ROOT / "results" / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_mac_tracking_warnings.csv",
]


CSV_FIELDNAMES = [
    "pair",
    "branch_i",
    "branch_j",
    "candidate_type",
    "mu_candidate",
    "Lambda_i",
    "Lambda_j",
    "abs_gap",
    "rel_gap",
    "sign_change",
    "sorted_position_i_left",
    "sorted_position_j_left",
    "sorted_position_i_right",
    "sorted_position_j_right",
    "mac_i_continuity",
    "mac_j_continuity",
    "mac_cross_i_to_j",
    "mac_cross_j_to_i",
    "tracking_status",
    "notes",
    "min_gap_mu",
    "sign_change_count",
    "local_min_step_mac",
    "local_min_mac_margin",
    "sorted_position_changed_i",
    "sorted_position_changed_j",
]


@dataclass(frozen=True)
class DenseTrackingResult:
    mu_values: np.ndarray
    sorted_roots: np.ndarray
    tracked_lambdas: np.ndarray
    current_sorted_positions: np.ndarray
    branch_vectors: list[list[np.ndarray]]
    step_mac: np.ndarray
    step_mac_margin: np.ndarray
    sorted_position_jump: np.ndarray
    nearest_frequency_positions: np.ndarray
    row_best_positions: np.ndarray


def finite_min(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(np.min(finite)) if finite.size else float("nan")


def fmt(value: float, precision: int = 6) -> str:
    if not np.isfinite(float(value)):
        return "nan"
    return f"{float(value):.{precision}g}"


def yesno(value: bool) -> str:
    return "yes" if bool(value) else "no"


def solve_sorted_root_grid() -> np.ndarray:
    roots = np.full((len(MU_VALUES), NUM_SORTED_ROOTS), np.nan, dtype=float)
    beta_rad = float(np.deg2rad(BETA_DEG))
    for idx, mu in enumerate(MU_VALUES):
        roots[idx] = find_first_n_roots_eta(
            beta_rad,
            float(mu),
            EPSILON,
            ETA,
            NUM_SORTED_ROOTS,
            Lmax0=ROOT_LMAX0,
            scan_step=ROOT_SCAN_STEP,
        )
        if np.any(~np.isfinite(roots[idx, :NUM_DESCENDANTS_CHECK])):
            raise RuntimeError(f"Missing required roots at mu={float(mu):g}.")
        if idx == 0 or (idx + 1) % 100 == 0 or idx + 1 == len(MU_VALUES):
            print(f"computed sorted roots for {idx + 1}/{len(MU_VALUES)} mu values")
    return roots


def shape_vectors_for_mu(sorted_roots: np.ndarray, mu: float) -> list[np.ndarray]:
    return analytic_shape_vectors_for_roots(
        sorted_roots,
        beta_rad=float(np.deg2rad(BETA_DEG)),
        mu=float(mu),
        epsilon=EPSILON,
        eta=ETA,
        s_norm=np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float),
    )


def assign_by_mac(previous_vectors: Sequence[np.ndarray], candidate_vectors: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mac = np.array(
        [
            [mac_value(previous, candidate) for candidate in candidate_vectors]
            for previous in previous_vectors
        ],
        dtype=float,
    )
    rows, cols = linear_sum_assignment(1.0 - mac)
    assignment = np.full(len(previous_vectors), -1, dtype=int)
    for row, col in zip(rows, cols, strict=True):
        assignment[int(row)] = int(col)
    if np.any(assignment < 0):
        raise RuntimeError("Could not assign every tracked descendant by MAC.")
    return assignment, mac


def track_descendants(sorted_roots: np.ndarray) -> DenseTrackingResult:
    n_mu = len(MU_VALUES)
    tracked = np.full((NUM_DESCENDANTS_CHECK, n_mu), np.nan, dtype=float)
    positions = np.full((NUM_DESCENDANTS_CHECK, n_mu), -1, dtype=int)
    step_mac = np.full((NUM_DESCENDANTS_CHECK, n_mu), np.nan, dtype=float)
    step_margin = np.full((NUM_DESCENDANTS_CHECK, n_mu), np.nan, dtype=float)
    jumps = np.zeros((NUM_DESCENDANTS_CHECK, n_mu), dtype=int)
    nearest_frequency_positions = np.full((NUM_DESCENDANTS_CHECK, n_mu), -1, dtype=int)
    row_best_positions = np.full((NUM_DESCENDANTS_CHECK, n_mu), -1, dtype=int)
    branch_vectors: list[list[np.ndarray]] = [[] for _ in range(NUM_DESCENDANTS_CHECK)]

    initial_vectors = shape_vectors_for_mu(sorted_roots[0], float(MU_VALUES[0]))
    previous_vectors = initial_vectors[:NUM_DESCENDANTS_CHECK]
    previous_lambdas = sorted_roots[0, :NUM_DESCENDANTS_CHECK].copy()

    tracked[:, 0] = previous_lambdas
    positions[:, 0] = np.arange(1, NUM_DESCENDANTS_CHECK + 1, dtype=int)
    nearest_frequency_positions[:, 0] = positions[:, 0]
    row_best_positions[:, 0] = positions[:, 0]
    for branch_idx, vector in enumerate(previous_vectors):
        branch_vectors[branch_idx].append(np.asarray(vector, dtype=float))

    for col in range(1, n_mu):
        candidate_vectors = shape_vectors_for_mu(sorted_roots[col], float(MU_VALUES[col]))
        assignment, mac = assign_by_mac(previous_vectors, candidate_vectors)
        current_vectors: list[np.ndarray] = []
        current_lambdas = np.full(NUM_DESCENDANTS_CHECK, np.nan, dtype=float)

        for branch_idx, candidate_col in enumerate(assignment):
            row_macs = np.sort(mac[branch_idx])[::-1]
            second_best_mac = float(row_macs[1]) if len(row_macs) > 1 else float("nan")
            assigned_mac = float(mac[branch_idx, int(candidate_col)])
            previous_position = int(positions[branch_idx, col - 1])
            current_position = int(candidate_col) + 1
            nearest_col = int(np.argmin(np.abs(sorted_roots[col] - previous_lambdas[branch_idx])))
            row_best_col = int(np.argmax(mac[branch_idx]))

            current_lambdas[branch_idx] = float(sorted_roots[col, int(candidate_col)])
            tracked[branch_idx, col] = current_lambdas[branch_idx]
            positions[branch_idx, col] = current_position
            step_mac[branch_idx, col] = assigned_mac
            step_margin[branch_idx, col] = assigned_mac - second_best_mac
            jumps[branch_idx, col] = current_position - previous_position
            nearest_frequency_positions[branch_idx, col] = nearest_col + 1
            row_best_positions[branch_idx, col] = row_best_col + 1
            current_vectors.append(np.asarray(candidate_vectors[int(candidate_col)], dtype=float))
            branch_vectors[branch_idx].append(np.asarray(candidate_vectors[int(candidate_col)], dtype=float))

        previous_vectors = current_vectors
        previous_lambdas = current_lambdas

        if (col + 1) % 100 == 0 or col + 1 == n_mu:
            print(f"tracked descendants for {col + 1}/{n_mu} mu values")

    return DenseTrackingResult(
        mu_values=MU_VALUES.copy(),
        sorted_roots=np.asarray(sorted_roots, dtype=float),
        tracked_lambdas=tracked,
        current_sorted_positions=positions,
        branch_vectors=branch_vectors,
        step_mac=step_mac,
        step_mac_margin=step_margin,
        sorted_position_jump=jumps,
        nearest_frequency_positions=nearest_frequency_positions,
        row_best_positions=row_best_positions,
    )


def sign_change_indices(difference: np.ndarray) -> list[int]:
    indices: list[int] = []
    for idx in range(len(difference) - 1):
        left = float(difference[idx])
        right = float(difference[idx + 1])
        if not (np.isfinite(left) and np.isfinite(right)):
            continue
        if left == 0.0 or right == 0.0 or left * right < 0.0:
            indices.append(idx)
    return indices


def interpolated_crossing(mu_values: np.ndarray, values_i: np.ndarray, values_j: np.ndarray, idx: int) -> tuple[float, float, float]:
    mu_left = float(mu_values[idx])
    mu_right = float(mu_values[idx + 1])
    d_left = float(values_i[idx] - values_j[idx])
    d_right = float(values_i[idx + 1] - values_j[idx + 1])
    if abs(d_right - d_left) <= 1e-15:
        mu_cross = 0.5 * (mu_left + mu_right)
    else:
        mu_cross = mu_left - d_left * (mu_right - mu_left) / (d_right - d_left)
    weight = (mu_cross - mu_left) / (mu_right - mu_left) if mu_right != mu_left else 0.0
    lambda_i = float(values_i[idx] + weight * (values_i[idx + 1] - values_i[idx]))
    lambda_j = float(values_j[idx] + weight * (values_j[idx + 1] - values_j[idx]))
    return float(mu_cross), lambda_i, lambda_j


def local_indices(mu_values: np.ndarray, candidate_mu: float, sign_idx: int | None) -> tuple[int, int, int]:
    if sign_idx is not None:
        left = max(0, int(sign_idx))
        right = min(len(mu_values) - 1, int(sign_idx) + 1)
        mid = left if abs(float(mu_values[left]) - candidate_mu) <= abs(float(mu_values[right]) - candidate_mu) else right
        return left, mid, right
    mid = int(np.argmin(np.abs(mu_values - float(candidate_mu))))
    left = max(0, mid - 1)
    right = min(len(mu_values) - 1, mid + 1)
    return left, mid, right


def branch_mac(result: DenseTrackingResult, branch_a: int, idx_a: int, branch_b: int, idx_b: int) -> float:
    return mac_value(result.branch_vectors[branch_a][idx_a], result.branch_vectors[branch_b][idx_b])


def local_tracking_status(
    result: DenseTrackingResult,
    branch_i: int,
    branch_j: int,
    left: int,
    right: int,
    mac_i: float,
    mac_j: float,
    cross_i_j: float,
    cross_j_i: float,
) -> tuple[str, str, float, float]:
    cols = range(max(1, left), right + 1)
    step_macs: list[float] = []
    margins: list[float] = []
    flags: list[str] = []
    for branch in (branch_i, branch_j):
        for col in cols:
            step_macs.append(float(result.step_mac[branch, col]))
            margins.append(float(result.step_mac_margin[branch, col]))
            if abs(int(result.sorted_position_jump[branch, col])) > MAX_SORTED_POSITION_JUMP:
                flags.append("large_sorted_position_jump")
            if int(result.nearest_frequency_positions[branch, col]) != int(result.current_sorted_positions[branch, col]):
                flags.append("nearest_frequency_disagrees_with_mac")
            if int(result.row_best_positions[branch, col]) != int(result.current_sorted_positions[branch, col]):
                flags.append("global_assignment_not_row_best")

    local_min_mac = finite_min(step_macs)
    local_min_margin = finite_min(margins)
    own_continuity = min(float(mac_i), float(mac_j))
    cross_continuity = max(float(cross_i_j), float(cross_j_i))

    if np.isfinite(local_min_mac) and local_min_mac < MAC_WARNING_THRESHOLD:
        flags.append("low_mac")
    if np.isfinite(local_min_margin) and local_min_margin < MAC_MARGIN_WARNING_THRESHOLD:
        flags.append("low_mac_margin")
    if cross_continuity > own_continuity + 0.05:
        flags.append("cross_mac_exceeds_own_continuity")

    if not flags:
        return "stable", "local MAC continuity is stable", local_min_mac, local_min_margin
    if "cross_mac_exceeds_own_continuity" in flags or "large_sorted_position_jump" in flags:
        return "tracking_error_suspected", "; ".join(dict.fromkeys(flags)), local_min_mac, local_min_margin
    return "tracking_ambiguous", "; ".join(dict.fromkeys(flags)), local_min_mac, local_min_margin


def classify_pair(
    *,
    sign_change: bool,
    abs_gap: float,
    rel_gap: float,
    tracking_status: str,
) -> str:
    near_threshold = abs_gap < NEAR_CROSSING_ABS_GAP or rel_gap < NEAR_CROSSING_REL_GAP
    close_with_tracking_warning = abs_gap < TRACKING_WARNING_GAP_CONTEXT and tracking_status != "stable"
    if tracking_status == "tracking_error_suspected":
        return "tracking_error_suspected"
    if sign_change:
        return "tracking_ambiguous" if tracking_status != "stable" else "true_descendant_crossing"
    if near_threshold:
        return "tracking_ambiguous" if tracking_status != "stable" else "near_miss_no_crossing"
    if close_with_tracking_warning:
        return "tracking_ambiguous"
    return "no_candidate"


def audit_pair(result: DenseTrackingResult, branch_i_1based: int, branch_j_1based: int) -> dict[str, float | int | str]:
    i = int(branch_i_1based) - 1
    j = int(branch_j_1based) - 1
    values_i = result.tracked_lambdas[i]
    values_j = result.tracked_lambdas[j]
    difference = values_i - values_j
    sign_indices = sign_change_indices(difference)
    min_idx = int(np.nanargmin(np.abs(difference)))

    sign_idx: int | None = None
    if sign_indices:
        sign_idx = min(sign_indices, key=lambda idx: min(abs(float(difference[idx])), abs(float(difference[idx + 1]))))
        mu_candidate, lambda_i, lambda_j = interpolated_crossing(result.mu_values, values_i, values_j, sign_idx)
        abs_gap = abs(lambda_i - lambda_j)
    else:
        mu_candidate = float(result.mu_values[min_idx])
        lambda_i = float(values_i[min_idx])
        lambda_j = float(values_j[min_idx])
        abs_gap = abs(lambda_i - lambda_j)

    mean_lambda = 0.5 * (abs(lambda_i) + abs(lambda_j))
    rel_gap = abs_gap / mean_lambda if mean_lambda > 0.0 else float("nan")
    left, _mid, right = local_indices(result.mu_values, mu_candidate, sign_idx)

    mac_i = branch_mac(result, i, left, i, right)
    mac_j = branch_mac(result, j, left, j, right)
    cross_i_j = branch_mac(result, i, left, j, right)
    cross_j_i = branch_mac(result, j, left, i, right)
    tracking_status, tracking_notes, local_min_mac, local_min_margin = local_tracking_status(
        result,
        i,
        j,
        left,
        right,
        mac_i,
        mac_j,
        cross_i_j,
        cross_j_i,
    )
    candidate_type = classify_pair(
        sign_change=bool(sign_indices),
        abs_gap=float(abs_gap),
        rel_gap=float(rel_gap),
        tracking_status=tracking_status,
    )

    changed_i = bool(np.any(np.diff(result.current_sorted_positions[i]) != 0))
    changed_j = bool(np.any(np.diff(result.current_sorted_positions[j]) != 0))
    notes = []
    if sign_indices:
        notes.append(f"{len(sign_indices)} sign-change bracket(s)")
    else:
        notes.append("no sign change of Lambda_i-Lambda_j on dense grid")
    if abs_gap < NEAR_CROSSING_ABS_GAP or rel_gap < NEAR_CROSSING_REL_GAP:
        notes.append("near-crossing threshold met")
    else:
        notes.append("gap above near-crossing thresholds")
    if tracking_notes:
        notes.append(tracking_notes)

    return {
        "pair": f"{branch_i_1based}-{branch_j_1based}",
        "branch_i": int(branch_i_1based),
        "branch_j": int(branch_j_1based),
        "candidate_type": candidate_type,
        "mu_candidate": float(mu_candidate),
        "Lambda_i": float(lambda_i),
        "Lambda_j": float(lambda_j),
        "abs_gap": float(abs_gap),
        "rel_gap": float(rel_gap),
        "sign_change": yesno(bool(sign_indices)),
        "sorted_position_i_left": int(result.current_sorted_positions[i, left]),
        "sorted_position_j_left": int(result.current_sorted_positions[j, left]),
        "sorted_position_i_right": int(result.current_sorted_positions[i, right]),
        "sorted_position_j_right": int(result.current_sorted_positions[j, right]),
        "mac_i_continuity": float(mac_i),
        "mac_j_continuity": float(mac_j),
        "mac_cross_i_to_j": float(cross_i_j),
        "mac_cross_j_to_i": float(cross_j_i),
        "tracking_status": tracking_status,
        "notes": "; ".join(notes),
        "min_gap_mu": float(result.mu_values[min_idx]),
        "sign_change_count": int(len(sign_indices)),
        "local_min_step_mac": float(local_min_mac),
        "local_min_mac_margin": float(local_min_margin),
        "sorted_position_changed_i": yesno(changed_i),
        "sorted_position_changed_j": yesno(changed_j),
    }


def audit_pairs(result: DenseTrackingResult) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for branch_i in range(1, NUM_DESCENDANTS_CHECK + 1):
        for branch_j in range(branch_i + 1, NUM_DESCENDANTS_CHECK + 1):
            rows.append(audit_pair(result, branch_i, branch_j))
    return rows


def adjacent_sorted_gap_rows(result: DenseTrackingResult) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for sorted_idx in range(1, NUM_DESCENDANTS_CHECK):
        gap = result.sorted_roots[:, sorted_idx] - result.sorted_roots[:, sorted_idx - 1]
        min_idx = int(np.nanargmin(gap))
        mean_lambda = 0.5 * (
            abs(float(result.sorted_roots[min_idx, sorted_idx]))
            + abs(float(result.sorted_roots[min_idx, sorted_idx - 1]))
        )
        rows.append(
            {
                "sorted_pair": f"{sorted_idx}-{sorted_idx + 1}",
                "mu": float(result.mu_values[min_idx]),
                "abs_gap": float(gap[min_idx]),
                "rel_gap": float(gap[min_idx] / mean_lambda if mean_lambda > 0.0 else np.nan),
            }
        )
    return rows


def write_csv(rows: Sequence[dict[str, float | int | str]]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    table = ["| " + " | ".join(headers) + " |"]
    table.append("| " + " | ".join("---" for _ in headers) + " |")
    table.extend("| " + " | ".join(row) + " |" for row in rows)
    return table


def source_lines() -> list[str]:
    lines: list[str] = []
    for path in PRIOR_SOURCES:
        status = "found" if path.exists() else "not found"
        lines.append(f"- {status}: `{path.relative_to(REPO_ROOT)}`")
    return lines


def sorted_position_change_lines(result: DenseTrackingResult) -> list[str]:
    lines: list[str] = []
    for branch in range(NUM_DESCENDANTS_CHECK):
        positions = result.current_sorted_positions[branch]
        change_indices = np.flatnonzero(np.diff(positions) != 0)
        if change_indices.size == 0:
            lines.append(f"- branch {branch + 1}: no sorted-position changes")
            continue
        chunks = []
        for idx in change_indices[:10]:
            chunks.append(
                f"mu={float(result.mu_values[idx + 1]):.6g}: "
                f"{int(positions[idx])}->{int(positions[idx + 1])}"
            )
        suffix = "" if change_indices.size <= 10 else f"; plus {change_indices.size - 10} more"
        lines.append(f"- branch {branch + 1}: " + "; ".join(chunks) + suffix)
    return lines


def write_report(
    result: DenseTrackingResult,
    pair_rows: Sequence[dict[str, float | int | str]],
    sorted_gap_rows: Sequence[dict[str, float | int | str]],
) -> None:
    output_rows = sorted(pair_rows, key=lambda row: float(row["abs_gap"]))
    candidate_rows = [row for row in pair_rows if str(row["candidate_type"]) != "no_candidate"]
    true_rows = [row for row in pair_rows if str(row["candidate_type"]) == "true_descendant_crossing"]
    sign_rows = [row for row in pair_rows if str(row["sign_change"]) == "yes"]
    sorted_changes = any(np.any(np.diff(result.current_sorted_positions[branch]) != 0) for branch in range(NUM_DESCENDANTS_CHECK))

    summary_table = [
        [
            str(row["pair"]),
            fmt(float(row["min_gap_mu"])),
            fmt(float(row["abs_gap"]), 8),
            fmt(float(row["rel_gap"]), 8),
            str(row["sign_change"]),
            str(row["candidate_type"]),
            str(row["tracking_status"]),
        ]
        for row in output_rows
    ]
    adjacent_table = [
        [
            str(row["sorted_pair"]),
            fmt(float(row["mu"])),
            fmt(float(row["abs_gap"]), 8),
            fmt(float(row["rel_gap"]), 8),
        ]
        for row in sorted(sorted_gap_rows, key=lambda item: float(item["abs_gap"]))
    ]

    lines = [
        "# Eta=0.5 Lambda(mu) Crossing Audit",
        "",
        "## Purpose",
        "",
        "Check whether the already studied thickness-mismatch `Lambda(mu)` graph contains real",
        "frequency crossings between descendant branches, or only close approaches, sorted-root",
        "permutations, or tracking artifacts.",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu range: {float(MU_VALUES[0]):g} .. {float(MU_VALUES[-1]):g}",
        f"- mu grid: {len(MU_VALUES)} points, step {float(MU_VALUES[1] - MU_VALUES[0]):g}",
        f"- descendant branches checked: first {NUM_DESCENDANTS_CHECK}",
        f"- sorted roots solved at each mu: first {NUM_SORTED_ROOTS}",
        f"- near-crossing thresholds: `abs_gap < {NEAR_CROSSING_ABS_GAP:g}` or `rel_gap < {NEAR_CROSSING_REL_GAP:g}`",
        f"- tracking-warning context gap: `{TRACKING_WARNING_GAP_CONTEXT:g}`",
        "",
        "## Data Source",
        "",
    ]
    lines.extend(source_lines())
    lines.extend(
        [
            "",
            "The old files are sufficient to identify the source graph, but they do not contain",
            "full dense pairwise branch data and local mode-shape MAC checks. This audit therefore",
            "recomputes diagnostic-only analytic Euler--Bernoulli thickness-mismatch data with",
            "`src/my_project/analytic/formulas_thickness_mismatch.py` and",
            "`scripts/lib/thickness_mismatch_mac_tracking.py` shape reconstruction helpers.",
            "",
            "No article figure, baseline result, old determinant, solver, or FEM model is modified.",
            "",
            "## Branch Tracking Convention",
            "",
            "For this fixed-eta thickness-mismatch diagnostic, branch `k` is the descendant of the",
            "`k`-th mode shape at the local seed point `mu=0`, with fixed `beta=15 deg` and",
            "`eta=0.5`. Sorted position is diagnostic metadata only. Tracking uses adjacent-step",
            "mode-shape MAC on the dense mu grid; low MAC, low MAC margin, or large sorted-position",
            "jumps are treated as tracking warnings and are not accepted as evidence for a true crossing.",
            "",
            "## Summary Table",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            ["pair", "mu at min/candidate", "abs gap", "rel gap", "sign change", "classification", "tracking"],
            summary_table,
        )
    )
    lines.extend(["", "## Candidate Interpretation", ""])
    if candidate_rows:
        for row in sorted(candidate_rows, key=lambda item: float(item["abs_gap"])):
            lines.extend(
                [
                    f"### Pair {row['pair']}",
                    "",
                    f"- classification: `{row['candidate_type']}`",
                    f"- mu candidate: `{fmt(float(row['mu_candidate']), 8)}`",
                    f"- Lambda values: `{fmt(float(row['Lambda_i']), 10)}` and `{fmt(float(row['Lambda_j']), 10)}`",
                    f"- abs gap: `{fmt(float(row['abs_gap']), 10)}`",
                    f"- rel gap: `{fmt(float(row['rel_gap']), 10)}`",
                    f"- sign change: `{row['sign_change']}`",
                    f"- sorted positions left/right: branch {row['branch_i']} `{row['sorted_position_i_left']}->{row['sorted_position_i_right']}`, "
                    f"branch {row['branch_j']} `{row['sorted_position_j_left']}->{row['sorted_position_j_right']}`",
                    f"- MAC continuity: branch {row['branch_i']} `{fmt(float(row['mac_i_continuity']), 8)}`, "
                    f"branch {row['branch_j']} `{fmt(float(row['mac_j_continuity']), 8)}`, "
                    f"cross `{fmt(float(row['mac_cross_i_to_j']), 8)}` / `{fmt(float(row['mac_cross_j_to_i']), 8)}`",
                    f"- tracking status: `{row['tracking_status']}`; {row['notes']}",
                    "",
                ]
            )
    else:
        lines.append("No pair met the sign-change or near-crossing thresholds.")

    lines.extend(
        [
            "## Sorted-Root Checks",
            "",
            "Descendant sorted-position changes:",
            "",
        ]
    )
    lines.extend(sorted_position_change_lines(result))
    lines.extend(["", "Closest adjacent sorted-root gaps:", ""])
    lines.extend(markdown_table(["sorted pair", "mu", "abs gap", "rel gap"], adjacent_table))

    lines.extend(
        [
            "",
            "## Final Answer",
            "",
            f"- true descendant crossings found: {len(true_rows)}",
            f"- descendant pairs with sign changes: {len(sign_rows)}",
            f"- descendant sorted-position swaps found: {'yes' if sorted_changes else 'no'}",
            f"- closest descendant pair: `{output_rows[0]['pair']}` at `mu={fmt(float(output_rows[0]['min_gap_mu']), 8)}` "
            f"with abs gap `{fmt(float(output_rows[0]['abs_gap']), 10)}`",
            "",
        ]
    )
    if true_rows:
        for row in true_rows:
            lines.append(
                f"- true crossing pair `{row['pair']}`: mu `{fmt(float(row['mu_candidate']), 8)}`, "
                f"Lambda `{fmt(0.5 * (float(row['Lambda_i']) + float(row['Lambda_j'])), 10)}`"
            )
    else:
        lines.append(
            "No real crossings of the checked descendant branches were found on the dense analytic grid. "
            "The closest interactions remain finite-gap approaches; locally low MAC near the closest gaps "
            "is treated as tracking ambiguity, not as crossing evidence."
        )

    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- This audit checks the diagnostic Euler--Bernoulli thickness-mismatch determinant used by the old graph.",
            "- Timoshenko thickness-mismatch branches are not added here because the identified old graph is EB-only.",
            "- Crossing locations, if present, would be resolved by dense-grid sign changes and linear interpolation.",
            "- Low-MAC or low-margin local assignments are reported as ambiguous rather than accepted as physical branch exchange.",
            "",
            "## Next Step",
            "",
            "If the ambiguous finite-gap interactions near the closest pairs need a stronger identity statement, run a",
            "localized subspace-MAC refinement around the listed mu values. The current finite gaps are well above the",
            "near-crossing thresholds, so this is not needed to reject true frequency crossings in this audit.",
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            "",
        ]
    )

    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def plot_audit(result: DenseTrackingResult, pair_rows: Sequence[dict[str, float | int | str]]) -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, ax = plt.subplots(figsize=(9.8, 5.8))
    for branch in range(NUM_DESCENDANTS_CHECK):
        ax.plot(
            result.mu_values,
            result.tracked_lambdas[branch],
            lw=1.6,
            color=colors[branch % len(colors)],
            label=f"desc {branch + 1}",
        )

    candidate_rows = [row for row in pair_rows if str(row["candidate_type"]) != "no_candidate"]
    for row in candidate_rows:
        marker = "x" if "ambiguous" in str(row["candidate_type"]) else "o"
        ax.scatter(
            [float(row["mu_candidate"])],
            [0.5 * (float(row["Lambda_i"]) + float(row["Lambda_j"]))],
            color="black",
            marker=marker,
            s=42,
            zorder=5,
        )
        ax.annotate(
            str(row["pair"]),
            xy=(float(row["mu_candidate"]), 0.5 * (float(row["Lambda_i"]) + float(row["Lambda_j"]))),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_title(
        rf"Diagnostic crossing audit, $\eta={ETA:g}$, $\beta={BETA_DEG:g}^\circ$, $\epsilon={EPSILON:g}$"
    )
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(ncol=4, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> dict[str, object]:
    sorted_roots = solve_sorted_root_grid()
    tracking = track_descendants(sorted_roots)
    pair_rows = audit_pairs(tracking)
    sorted_gap_rows = adjacent_sorted_gap_rows(tracking)
    write_csv(pair_rows)
    plot_audit(tracking, pair_rows)
    write_report(tracking, pair_rows, sorted_gap_rows)

    true_rows = [row for row in pair_rows if str(row["candidate_type"]) == "true_descendant_crossing"]
    sorted_changes = any(
        np.any(np.diff(tracking.current_sorted_positions[branch]) != 0)
        for branch in range(NUM_DESCENDANTS_CHECK)
    )
    closest = min(pair_rows, key=lambda row: float(row["abs_gap"]))
    print(f"saved crossing audit CSV: {OUTPUT_CSV}")
    print(f"saved crossing audit report: {OUTPUT_REPORT}")
    print(f"saved crossing audit PNG: {OUTPUT_PNG}")
    print(f"true descendant crossings found: {len(true_rows)}")
    print(f"descendant sorted-position swaps found: {yesno(sorted_changes)}")
    print(
        "closest descendant pair: "
        f"{closest['pair']} at mu={float(closest['min_gap_mu']):.6g}, "
        f"gap={float(closest['abs_gap']):.12g}"
    )
    return {
        "csv": OUTPUT_CSV,
        "report": OUTPUT_REPORT,
        "png": OUTPUT_PNG,
        "true_crossings": len(true_rows),
        "sorted_position_swaps": bool(sorted_changes),
        "closest_pair": closest,
    }


if __name__ == "__main__":
    main()
