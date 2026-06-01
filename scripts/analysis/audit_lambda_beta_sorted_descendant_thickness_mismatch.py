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

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import (  # noqa: E402
    analytic_shape_vectors_for_roots,
    mac_assignment,
    unique_nearest_frequency_assignment,
)


# =========================
# User-editable parameters
# =========================
ETA_VALUES = (0.5, 0.0)
EPSILON = 0.0025
MU = 0.0
BETA_VALUES_DEG = np.linspace(0.0, 90.0, 901)

N_DESCENDANTS = 8
N_SORTED_ROOTS = 14
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
TARGET_BETA_DEG = 15.0

OUTPUT_DIR = REPO_ROOT / "results"
SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

AUDIT_FIELDNAMES = [
    "eta",
    "epsilon",
    "mu",
    "beta_deg",
    "descendant_id",
    "Lambda_descendant",
    "sorted_position",
    "sorted_index_value_at_same_lambda",
    "mac_to_previous",
    "tracking_warning",
    "is_sorted5_at_this_beta",
    "sorted5_descendant_id_at_beta",
    "sorted5_lambda",
    "gap_to_sorted4",
    "gap_to_sorted6",
]

SUMMARY_FIELDNAMES = [
    "eta",
    "beta_interval_start",
    "beta_interval_end",
    "sorted5_descendant_id",
    "min_gap_to_sorted4",
    "min_gap_to_sorted6",
    "notes",
]


@dataclass(frozen=True)
class BetaTrackingResult:
    eta: float
    beta_values_deg: np.ndarray
    sorted_roots: np.ndarray
    tracked_lambdas: np.ndarray
    current_sorted_positions: np.ndarray
    rows: list[dict[str, float | int | str]]
    warning_rows: list[dict[str, float | int | str]]


def number_token(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def eta_stem(value: float) -> str:
    value_f = float(value)
    if np.isclose(value_f, 0.0, rtol=0.0, atol=1e-12):
        return "eta0"
    prefix = "eta" if value_f > 0.0 else "eta_m"
    return prefix + number_token(abs(value_f))


def output_stem(eta: float) -> str:
    return (
        f"{eta_stem(float(eta))}_eps{number_token(EPSILON)}_mu{number_token(MU)}"
        "_lambda_beta_sorted_descendant_audit"
    )


def summary_stem(eta: float) -> str:
    return (
        f"{eta_stem(float(eta))}_eps{number_token(EPSILON)}_mu{number_token(MU)}"
        "_sorted5_identity_summary"
    )


def output_paths(eta: float) -> dict[str, Path]:
    stem = output_stem(float(eta))
    focus_stem = (
        f"{eta_stem(float(eta))}_eps{number_token(EPSILON)}_mu{number_token(MU)}"
        "_sorted5_identity_vs_beta"
    )
    return {
        "csv": OUTPUT_DIR / f"{stem}.csv",
        "png": OUTPUT_DIR / f"{stem}.png",
        "summary_csv": OUTPUT_DIR / f"{summary_stem(float(eta))}.csv",
        "report": OUTPUT_DIR / f"{stem}_report.md",
        "focus_png": OUTPUT_DIR / f"{focus_stem}.png",
    }


def solve_sorted_root_grid(*, eta: float) -> np.ndarray:
    roots = np.full((len(BETA_VALUES_DEG), N_SORTED_ROOTS), np.nan, dtype=float)
    for idx, beta_deg in enumerate(BETA_VALUES_DEG):
        values = find_first_n_roots_eta(
            float(np.deg2rad(beta_deg)),
            MU,
            EPSILON,
            float(eta),
            N_SORTED_ROOTS,
            Lmax0=ROOT_LMAX0,
            scan_step=ROOT_SCAN_STEP,
        )
        if np.any(~np.isfinite(values[:N_DESCENDANTS])):
            raise RuntimeError(f"Missing roots at eta={float(eta):g}, beta={float(beta_deg):g} deg.")
        roots[idx] = values
        if idx == 0 or (idx + 1) % 100 == 0 or idx + 1 == len(BETA_VALUES_DEG):
            print(
                f"eta={float(eta):g}: computed sorted roots for "
                f"{idx + 1}/{len(BETA_VALUES_DEG)} beta values"
            )
    return roots


def shape_vectors_for_beta(sorted_roots: np.ndarray, *, eta: float, beta_deg: float) -> list[np.ndarray]:
    return analytic_shape_vectors_for_roots(
        sorted_roots,
        beta_rad=float(np.deg2rad(beta_deg)),
        mu=MU,
        epsilon=EPSILON,
        eta=float(eta),
        s_norm=np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float),
    )


def yesno(value: bool) -> str:
    return "yes" if bool(value) else "no"


def tracking_warning(row: Mapping[str, float | int | str]) -> str:
    flags = []
    for key, label in (
        ("low_mac", "low_mac"),
        ("low_margin", "low_margin"),
        ("large_beta_step", "large_beta_step"),
        ("blocked_by_unresolved_neighbor", "blocked_by_unresolved_neighbor"),
        ("unresolved_assignment", "unresolved_assignment"),
        ("frequency_mac_disagreement", "frequency_mac_disagreement"),
    ):
        if str(row.get(key, "no")) == "yes":
            flags.append(label)
    if str(row.get("tracking_step_status", "")) not in {"", "seed_beta0", "mac_ok"}:
        status = str(row["tracking_step_status"])
        if status not in flags:
            flags.append(status)
    return ";".join(flags)


def track_descendants_by_beta(*, eta: float, sorted_roots: np.ndarray) -> BetaTrackingResult:
    beta_grid = np.asarray(BETA_VALUES_DEG, dtype=float)
    if not np.isclose(float(beta_grid[0]), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Beta descendants must be seeded at beta=0.")
    if sorted_roots.shape[1] < N_DESCENDANTS:
        raise ValueError("sorted_roots does not contain enough roots to seed descendants.")

    tracked = np.full((N_DESCENDANTS, len(beta_grid)), np.nan, dtype=float)
    positions = np.full((N_DESCENDANTS, len(beta_grid)), -1, dtype=int)
    rows: list[dict[str, float | int | str]] = []
    warning_rows: list[dict[str, float | int | str]] = []

    initial_vectors = shape_vectors_for_beta(sorted_roots[0], eta=float(eta), beta_deg=float(beta_grid[0]))
    previous_vectors = initial_vectors[:N_DESCENDANTS]
    previous_lambdas = sorted_roots[0, :N_DESCENDANTS].copy()

    tracked[:, 0] = previous_lambdas
    positions[:, 0] = np.arange(1, N_DESCENDANTS + 1, dtype=int)
    max_beta_step = float(beta_grid[1] - beta_grid[0]) if len(beta_grid) > 1 else None

    for branch_idx, value in enumerate(previous_lambdas, start=1):
        rows.append(
            {
                "eta": float(eta),
                "beta_prev_deg": np.nan,
                "beta_deg": float(beta_grid[0]),
                "descendant_id": int(branch_idx),
                "Lambda_tracked": float(value),
                "current_sorted_position": int(branch_idx),
                "diagnostic_candidate_sorted_position": int(branch_idx),
                "diagnostic_candidate_Lambda": float(value),
                "nearest_sorted_root_index": int(branch_idx),
                "tracking_step_status": "seed_beta0",
                "mac_to_previous": np.nan,
                "accepted_mac_to_previous": np.nan,
                "diagnostic_candidate_mac_to_previous": np.nan,
                "second_best_mac": np.nan,
                "frequency_distance_from_previous": np.nan,
                "frequency_second_nearest_distance": np.nan,
                "frequency_assignment_margin": np.nan,
                "frequency_mac_disagreement": "no",
                "mac_margin": np.nan,
                "assigned_from_previous_sorted_position": int(branch_idx),
                "sorted_position_jump": 0,
                "diagnostic_candidate_sorted_position_jump": 0,
                "low_mac": "no",
                "low_margin": "no",
                "large_beta_step": "no",
                "blocked_by_unresolved_neighbor": "no",
                "unresolved_assignment": "no",
                "suspicious_assignment": "no",
                "requires_refined_check": "no",
            }
        )

    for col in range(1, len(beta_grid)):
        beta_prev = float(beta_grid[col - 1])
        beta = float(beta_grid[col])
        roots = np.asarray(sorted_roots[col], dtype=float)
        candidate_vectors = shape_vectors_for_beta(roots, eta=float(eta), beta_deg=beta)
        raw_cols, mac = mac_assignment(previous_vectors, candidate_vectors)
        freq_assignments = unique_nearest_frequency_assignment(previous_lambdas, roots)
        previous_positions = positions[:, col - 1].astype(int)

        second_best_values = np.full(N_DESCENDANTS, np.nan, dtype=float)
        candidate_macs = np.full(N_DESCENDANTS, np.nan, dtype=float)
        mac_margins = np.full(N_DESCENDANTS, np.nan, dtype=float)
        candidate_jumps = np.zeros(N_DESCENDANTS, dtype=int)
        low_mac_flags = np.zeros(N_DESCENDANTS, dtype=bool)
        low_margin_flags = np.zeros(N_DESCENDANTS, dtype=bool)
        large_step_flags = np.zeros(N_DESCENDANTS, dtype=bool)
        unresolved = np.zeros(N_DESCENDANTS, dtype=bool)

        for branch_row in range(N_DESCENDANTS):
            candidate_col = int(raw_cols[branch_row])
            row_macs = np.sort(mac[branch_row])[::-1]
            second_best = float(row_macs[1]) if len(row_macs) > 1 else np.nan
            candidate_mac = float(mac[branch_row, candidate_col])
            mac_margin = candidate_mac - second_best if np.isfinite(second_best) else np.nan
            candidate_jump = int(candidate_col) + 1 - int(previous_positions[branch_row])
            low_mac = candidate_mac < MAC_WARNING_THRESHOLD
            low_margin = np.isfinite(mac_margin) and mac_margin < MAC_MARGIN_WARNING_THRESHOLD
            large_step = (
                max_beta_step is not None and abs(beta - beta_prev) > max_beta_step + 1e-12
            )

            second_best_values[branch_row] = second_best
            candidate_macs[branch_row] = candidate_mac
            mac_margins[branch_row] = mac_margin
            candidate_jumps[branch_row] = candidate_jump
            low_mac_flags[branch_row] = bool(low_mac)
            low_margin_flags[branch_row] = bool(low_margin)
            large_step_flags[branch_row] = bool(large_step)
            unresolved[branch_row] = (
                abs(candidate_jump) > MAX_SORTED_POSITION_JUMP
                or low_mac
                or low_margin
                or large_step
            )

        blocked_by_unresolved_neighbor = np.zeros(N_DESCENDANTS, dtype=bool)
        while True:
            retained_positions = {
                int(previous_positions[row])
                for row in range(N_DESCENDANTS)
                if bool(unresolved[row])
            }
            updated = unresolved.copy()
            for branch_row in range(N_DESCENDANTS):
                if bool(updated[branch_row]):
                    continue
                candidate_position = int(raw_cols[branch_row]) + 1
                if candidate_position in retained_positions and candidate_position != int(previous_positions[branch_row]):
                    updated[branch_row] = True
                    blocked_by_unresolved_neighbor[branch_row] = True
            if np.array_equal(updated, unresolved):
                break
            unresolved = updated

        accepted_cols = np.array(
            [
                int(previous_positions[row]) - 1 if bool(unresolved[row]) else int(raw_cols[row])
                for row in range(N_DESCENDANTS)
            ],
            dtype=int,
        )
        current_lambdas = np.full(N_DESCENDANTS, np.nan, dtype=float)
        current_vectors: list[np.ndarray] = []

        for branch_row in range(N_DESCENDANTS):
            candidate_col = int(raw_cols[branch_row])
            root_col = int(accepted_cols[branch_row])
            frequency = freq_assignments[branch_row]
            previous_sorted_position = int(previous_positions[branch_row])
            candidate_mac = float(candidate_macs[branch_row])
            second_best = float(second_best_values[branch_row])
            mac_margin = float(mac_margins[branch_row])
            diagnostic_candidate_jump = int(candidate_jumps[branch_row])
            unresolved_assignment = bool(unresolved[branch_row])
            if root_col < 0 or root_col >= len(roots):
                raise RuntimeError(
                    "Cannot retain previous canonical sorted position "
                    f"{previous_sorted_position} at beta={beta:g}; only {len(roots)} roots are available."
                )

            sorted_position_jump = int(root_col) + 1 - previous_sorted_position
            accepted_mac = float(mac[branch_row, root_col])
            disagreement = int(frequency.root_index) != int(candidate_col) + 1
            status = "unresolved_assignment" if unresolved_assignment else "mac_ok"
            row = {
                "eta": float(eta),
                "beta_prev_deg": beta_prev,
                "beta_deg": beta,
                "descendant_id": int(branch_row) + 1,
                "Lambda_tracked": float(roots[root_col]),
                "current_sorted_position": int(root_col) + 1,
                "diagnostic_candidate_sorted_position": int(candidate_col) + 1,
                "diagnostic_candidate_Lambda": float(roots[candidate_col]),
                "nearest_sorted_root_index": int(frequency.root_index),
                "nearest_sorted_Lambda": float(roots[int(frequency.root_index) - 1]),
                "tracking_step_status": status,
                "mac_to_previous": candidate_mac,
                "accepted_mac_to_previous": accepted_mac,
                "diagnostic_candidate_mac_to_previous": candidate_mac,
                "second_best_mac": second_best,
                "frequency_distance_from_previous": float(frequency.distance_from_previous),
                "frequency_second_nearest_distance": float(frequency.second_nearest_distance),
                "frequency_assignment_margin": float(frequency.assignment_margin),
                "frequency_mac_disagreement": yesno(disagreement),
                "mac_margin": mac_margin,
                "assigned_from_previous_sorted_position": previous_sorted_position,
                "sorted_position_jump": int(sorted_position_jump),
                "diagnostic_candidate_sorted_position_jump": int(diagnostic_candidate_jump),
                "low_mac": yesno(bool(low_mac_flags[branch_row])),
                "low_margin": yesno(bool(low_margin_flags[branch_row])),
                "large_beta_step": yesno(bool(large_step_flags[branch_row])),
                "blocked_by_unresolved_neighbor": yesno(bool(blocked_by_unresolved_neighbor[branch_row])),
                "unresolved_assignment": yesno(unresolved_assignment),
                "suspicious_assignment": yesno(unresolved_assignment),
                "requires_refined_check": yesno(unresolved_assignment),
            }
            rows.append(row)
            if tracking_warning(row):
                warning_rows.append(row)
            current_lambdas[branch_row] = float(roots[root_col])
            positions[branch_row, col] = int(root_col) + 1
            current_vectors.append(candidate_vectors[root_col])

        tracked[:, col] = current_lambdas
        previous_lambdas = current_lambdas
        previous_vectors = current_vectors

        if (col + 1) % 100 == 0 or col + 1 == len(beta_grid):
            print(f"eta={float(eta):g}: tracked descendants for {col + 1}/{len(beta_grid)} beta values")

    return BetaTrackingResult(
        eta=float(eta),
        beta_values_deg=beta_grid,
        sorted_roots=np.asarray(sorted_roots, dtype=float),
        tracked_lambdas=tracked,
        current_sorted_positions=positions,
        rows=rows,
        warning_rows=warning_rows,
    )


def sorted5_identity_by_beta(result: BetaTrackingResult) -> list[str]:
    identities: list[str] = []
    for col in range(len(result.beta_values_deg)):
        occupants = [
            int(branch_idx) + 1
            for branch_idx in range(N_DESCENDANTS)
            if int(result.current_sorted_positions[branch_idx, col]) == 5
        ]
        if len(occupants) == 1:
            identities.append(str(occupants[0]))
        elif not occupants:
            identities.append("none")
        else:
            identities.append("ambiguous:" + ";".join(str(value) for value in occupants))
    return identities


def finite_float_text(value: float) -> str:
    return "" if not np.isfinite(float(value)) else f"{float(value):.12g}"


def write_audit_csv(result: BetaTrackingResult, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    identities = sorted5_identity_by_beta(result)
    sorted5 = result.sorted_roots[:, 4]
    gap45 = result.sorted_roots[:, 4] - result.sorted_roots[:, 3]
    gap56 = result.sorted_roots[:, 5] - result.sorted_roots[:, 4]
    rows_by_key = {
        (round(float(row["beta_deg"]), 10), int(row["descendant_id"])): row
        for row in result.rows
    }
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=AUDIT_FIELDNAMES)
        writer.writeheader()
        for col, beta_deg in enumerate(result.beta_values_deg):
            roots = result.sorted_roots[col]
            for branch_idx in range(N_DESCENDANTS):
                row = rows_by_key[(round(float(beta_deg), 10), branch_idx + 1)]
                lambda_desc = float(result.tracked_lambdas[branch_idx, col])
                nearest_sorted = int(np.nanargmin(np.abs(roots - lambda_desc))) + 1
                is_sorted5 = int(result.current_sorted_positions[branch_idx, col]) == 5
                mac_value = float(row["mac_to_previous"])
                writer.writerow(
                    {
                        "eta": float(result.eta),
                        "epsilon": EPSILON,
                        "mu": MU,
                        "beta_deg": float(beta_deg),
                        "descendant_id": branch_idx + 1,
                        "Lambda_descendant": lambda_desc,
                        "sorted_position": int(result.current_sorted_positions[branch_idx, col]),
                        "sorted_index_value_at_same_lambda": nearest_sorted,
                        "mac_to_previous": finite_float_text(mac_value),
                        "tracking_warning": tracking_warning(row),
                        "is_sorted5_at_this_beta": yesno(is_sorted5),
                        "sorted5_descendant_id_at_beta": identities[col],
                        "sorted5_lambda": float(sorted5[col]),
                        "gap_to_sorted4": float(gap45[col]),
                        "gap_to_sorted6": float(gap56[col]),
                    }
                )


def contiguous_identity_intervals(result: BetaTrackingResult) -> list[dict[str, float | str]]:
    identities = sorted5_identity_by_beta(result)
    gap45 = result.sorted_roots[:, 4] - result.sorted_roots[:, 3]
    gap56 = result.sorted_roots[:, 5] - result.sorted_roots[:, 4]
    warning_betas = {
        round(float(row["beta_deg"]), 10)
        for row in result.warning_rows
    }
    rows: list[dict[str, float | str]] = []
    start = 0
    for idx in range(1, len(identities) + 1):
        if idx < len(identities) and identities[idx] == identities[start]:
            continue
        interval_slice = slice(start, idx)
        beta_start = float(result.beta_values_deg[start])
        beta_end = float(result.beta_values_deg[idx - 1])
        interval_betas = result.beta_values_deg[interval_slice]
        notes = []
        if beta_start <= TARGET_BETA_DEG <= beta_end:
            notes.append(f"contains beta={TARGET_BETA_DEG:g} deg")
        if any(round(float(beta), 10) in warning_betas for beta in interval_betas):
            notes.append("tracking warning present")
        rows.append(
            {
                "eta": f"{result.eta:g}",
                "beta_interval_start": f"{beta_start:.10g}",
                "beta_interval_end": f"{beta_end:.10g}",
                "sorted5_descendant_id": identities[start],
                "min_gap_to_sorted4": f"{float(np.nanmin(gap45[interval_slice])):.12g}",
                "min_gap_to_sorted6": f"{float(np.nanmin(gap56[interval_slice])):.12g}",
                "notes": "; ".join(notes),
            }
        )
        start = idx
    return rows


def write_summary_csv(result: BetaTrackingResult, path: Path) -> list[dict[str, float | str]]:
    rows = contiguous_identity_intervals(result)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return rows


def target_col(result: BetaTrackingResult, target_beta: float = TARGET_BETA_DEG) -> int:
    matches = np.flatnonzero(np.isclose(result.beta_values_deg, float(target_beta), rtol=0.0, atol=1e-12))
    if matches.size:
        return int(matches[0])
    return int(np.argmin(np.abs(result.beta_values_deg - float(target_beta))))


def local_minima_indices(values: np.ndarray, *, top_n: int = 6) -> list[int]:
    candidates: list[int] = []
    for idx in range(len(values)):
        left = values[idx - 1] if idx > 0 else np.inf
        right = values[idx + 1] if idx + 1 < len(values) else np.inf
        if float(values[idx]) <= float(left) and float(values[idx]) <= float(right):
            candidates.append(idx)
    candidates.sort(key=lambda idx: float(values[idx]))
    return candidates[: int(top_n)]


def nearest_gap_summary(result: BetaTrackingResult) -> dict[str, object]:
    gap45 = result.sorted_roots[:, 4] - result.sorted_roots[:, 3]
    gap56 = result.sorted_roots[:, 5] - result.sorted_roots[:, 4]
    idx45 = int(np.nanargmin(gap45))
    idx56 = int(np.nanargmin(gap56))
    return {
        "gap45_min_beta": float(result.beta_values_deg[idx45]),
        "gap45_min_value": float(gap45[idx45]),
        "gap56_min_beta": float(result.beta_values_deg[idx56]),
        "gap56_min_value": float(gap56[idx56]),
        "gap45_local_minima": local_minima_indices(gap45),
        "gap56_local_minima": local_minima_indices(gap56),
    }


def plot_audit(result: BetaTrackingResult, output_png: Path) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    identities = sorted5_identity_by_beta(result)
    target_idx = target_col(result)
    target_identity = identities[target_idx]

    fig, axes = plt.subplots(2, 1, figsize=(10.2, 7.2), sharex=True, height_ratios=[2.4, 1.0])
    ax_lambda, ax_identity = axes
    for branch_idx in range(N_DESCENDANTS):
        ax_lambda.plot(
            result.beta_values_deg,
            result.tracked_lambdas[branch_idx],
            color=colors[branch_idx % len(colors)],
            lw=1.45,
            label=f"desc {branch_idx + 1}",
        )
    ax_lambda.plot(
        result.beta_values_deg,
        result.sorted_roots[:, 4],
        color="black",
        lw=2.0,
        ls="--",
        label="sorted 5",
        zorder=5,
    )
    ax_lambda.axvline(TARGET_BETA_DEG, color="0.35", lw=0.9, ls=":", zorder=1)
    ax_lambda.scatter(
        [float(result.beta_values_deg[target_idx])],
        [float(result.sorted_roots[target_idx, 4])],
        color="black",
        s=34,
        zorder=6,
    )
    ax_lambda.annotate(
        f"beta={float(result.beta_values_deg[target_idx]):g} deg: sorted 5 = desc {target_identity}",
        xy=(float(result.beta_values_deg[target_idx]), float(result.sorted_roots[target_idx, 4])),
        xytext=(8, 8),
        textcoords="offset points",
        fontsize=8,
    )
    ax_lambda.set_ylabel(r"$\Lambda$")
    ax_lambda.set_title(
        rf"EB thickness-mismatch beta audit, $\eta={result.eta:g}$, "
        rf"$\epsilon={EPSILON:g}$, $\mu={MU:g}$"
    )
    ax_lambda.grid(True, color="0.88", linewidth=0.6)
    ax_lambda.legend(loc="best", ncol=5, fontsize=8, frameon=False)

    identity_numeric = np.array(
        [float(value) if value.isdigit() else np.nan for value in identities],
        dtype=float,
    )
    ax_identity.step(result.beta_values_deg, identity_numeric, where="post", color="black", lw=1.4)
    ax_identity.scatter(result.beta_values_deg, identity_numeric, c=identity_numeric, cmap="tab10", s=8, zorder=3)
    ax_identity.axvline(TARGET_BETA_DEG, color="0.35", lw=0.9, ls=":")
    ax_identity.set_xlabel(r"$\beta$ (deg)")
    ax_identity.set_ylabel("sorted 5 descendant")
    ax_identity.set_yticks(range(1, N_DESCENDANTS + 1))
    ax_identity.set_ylim(0.5, N_DESCENDANTS + 0.5)
    ax_identity.grid(True, color="0.88", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(output_png, dpi=240, bbox_inches="tight")
    plt.close(fig)


def plot_sorted5_focus(result: BetaTrackingResult, output_png: Path) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    identities = sorted5_identity_by_beta(result)
    fig, ax = plt.subplots(figsize=(10.0, 4.4))
    for descendant_id in range(1, N_DESCENDANTS + 1):
        mask = np.array([identity == str(descendant_id) for identity in identities], dtype=bool)
        for start, end in contiguous_true_runs(mask):
            sl = slice(start, end + 1)
            ax.plot(
                result.beta_values_deg[sl],
                result.sorted_roots[sl, 4],
                color=colors[(descendant_id - 1) % len(colors)],
                lw=2.0,
                label=f"desc {descendant_id}" if start == 0 or not any(
                    identities[idx] == str(descendant_id) for idx in range(start)
                ) else None,
            )
    target_idx = target_col(result)
    ax.axvline(TARGET_BETA_DEG, color="0.35", lw=0.9, ls=":")
    ax.scatter(
        [float(result.beta_values_deg[target_idx])],
        [float(result.sorted_roots[target_idx, 4])],
        color="black",
        s=34,
        zorder=5,
    )
    ax.set_xlabel(r"$\beta$ (deg)")
    ax.set_ylabel(r"sorted 5 $\Lambda$")
    ax.set_title(
        rf"Sorted 5 identity vs beta, $\eta={result.eta:g}$, "
        rf"$\epsilon={EPSILON:g}$, $\mu={MU:g}$"
    )
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", ncol=4, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=240, bbox_inches="tight")
    plt.close(fig)


def contiguous_true_runs(mask: Sequence[bool]) -> list[tuple[int, int]]:
    values = np.asarray(mask, dtype=bool)
    runs: list[tuple[int, int]] = []
    start: int | None = None
    for idx, value in enumerate(values):
        if bool(value) and start is None:
            start = idx
        elif not bool(value) and start is not None:
            runs.append((start, idx - 1))
            start = None
    if start is not None:
        runs.append((start, len(values) - 1))
    return runs


def sorted_position_change_lines(result: BetaTrackingResult) -> list[str]:
    lines: list[str] = []
    for branch_idx in range(N_DESCENDANTS):
        positions = result.current_sorted_positions[branch_idx]
        change_indices = np.flatnonzero(np.diff(positions) != 0)
        if change_indices.size == 0:
            lines.append(f"- descendant {branch_idx + 1}: no sorted-position changes")
            continue
        chunks = []
        for idx in change_indices[:14]:
            chunks.append(
                f"beta={float(result.beta_values_deg[idx + 1]):.6g}: "
                f"{int(positions[idx])}->{int(positions[idx + 1])}"
            )
        suffix = "" if change_indices.size <= 14 else f"; plus {change_indices.size - 14} more"
        lines.append(f"- descendant {branch_idx + 1}: " + "; ".join(chunks) + suffix)
    return lines


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    table = ["| " + " | ".join(headers) + " |"]
    table.append("| " + " | ".join("---" for _ in headers) + " |")
    table.extend("| " + " | ".join(row) + " |" for row in rows)
    return table


def result_at_beta(result: BetaTrackingResult, *, beta_deg: float = TARGET_BETA_DEG) -> dict[str, object]:
    col = target_col(result, beta_deg)
    identities = sorted5_identity_by_beta(result)
    desc6_idx = 5
    return {
        "beta_deg": float(result.beta_values_deg[col]),
        "sorted5_lambda": float(result.sorted_roots[col, 4]),
        "sorted5_descendant_id": identities[col],
        "desc6_sorted_position": int(result.current_sorted_positions[desc6_idx, col]),
        "desc6_lambda": float(result.tracked_lambdas[desc6_idx, col]),
        "sorted5_equals_desc6": identities[col] == "6",
    }


def warning_summary(result: BetaTrackingResult) -> dict[str, int]:
    rows = result.rows
    return {
        "warning_rows": len(result.warning_rows),
        "low_mac": sum(1 for row in rows if str(row.get("low_mac", "no")) == "yes"),
        "low_margin": sum(1 for row in rows if str(row.get("low_margin", "no")) == "yes"),
        "unresolved": sum(1 for row in rows if str(row.get("unresolved_assignment", "no")) == "yes"),
        "frequency_mac_disagreement": sum(
            1 for row in rows if str(row.get("frequency_mac_disagreement", "no")) == "yes"
        ),
    }


def write_report(
    *,
    main_result: BetaTrackingResult,
    comparison_result: BetaTrackingResult | None,
    main_summary_rows: Sequence[dict[str, float | str]],
    comparison_summary_rows: Sequence[dict[str, float | str]] | None,
    paths: Mapping[str, Path],
) -> None:
    target_main = result_at_beta(main_result)
    target_comparison = result_at_beta(comparison_result) if comparison_result is not None else None
    gaps = nearest_gap_summary(main_result)
    warnings = warning_summary(main_result)
    comparison_warnings = warning_summary(comparison_result) if comparison_result is not None else None
    intervals = [
        [
            str(row["beta_interval_start"]),
            str(row["beta_interval_end"]),
            str(row["sorted5_descendant_id"]),
            str(row["min_gap_to_sorted4"]),
            str(row["min_gap_to_sorted6"]),
            str(row["notes"]),
        ]
        for row in main_summary_rows
    ]
    comparison_intervals = [
        [
            str(row["beta_interval_start"]),
            str(row["beta_interval_end"]),
            str(row["sorted5_descendant_id"]),
            str(row["min_gap_to_sorted4"]),
            str(row["min_gap_to_sorted6"]),
            str(row["notes"]),
        ]
        for row in (comparison_summary_rows or [])
    ]
    gap45_local = [
        [
            f"{float(main_result.beta_values_deg[idx]):.6g}",
            f"{float(main_result.sorted_roots[idx, 4] - main_result.sorted_roots[idx, 3]):.10g}",
        ]
        for idx in gaps["gap45_local_minima"]
    ]
    gap56_local = [
        [
            f"{float(main_result.beta_values_deg[idx]):.6g}",
            f"{float(main_result.sorted_roots[idx, 5] - main_result.sorted_roots[idx, 4]):.10g}",
        ]
        for idx in gaps["gap56_local_minima"]
    ]

    correction_needed = not bool(target_main["sorted5_equals_desc6"])
    lines = [
        "# Eta=0.5 Lambda(beta) Sorted/Descendant Audit",
        "",
        "## 1. What Was Checked",
        "",
        "This analytic-only Euler-Bernoulli thickness-mismatch diagnostic tracks",
        "descendant identities over beta and compares them with independently",
        "sorted roots at each beta. The target question is whether sorted",
        "frequency 5 at `beta=15 deg`, `eta=0.5`, `epsilon=0.0025`, `mu=0`",
        "is descendant 6 or a different descendant.",
        "",
        "No FEM workflow, Gmsh, CalculiX, article file, article figure, old",
        "determinant, old solver, or baseline result is modified by this audit.",
        "",
        "## 2. Parameters",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- eta main case: {main_result.eta:g}",
        f"- eta comparison case: {comparison_result.eta:g}" if comparison_result else "- eta comparison case: not run",
        f"- epsilon: {EPSILON:g}",
        f"- mu: {MU:g}",
        f"- beta range: {float(main_result.beta_values_deg[0]):g} .. {float(main_result.beta_values_deg[-1]):g} deg",
        f"- beta grid points: {len(main_result.beta_values_deg)}",
        f"- beta step: {float(main_result.beta_values_deg[1] - main_result.beta_values_deg[0]):g} deg",
        f"- tracked descendants: first {N_DESCENDANTS}",
        f"- sorted roots solved at each beta: first {N_SORTED_ROOTS}",
        f"- target beta: {float(target_main['beta_deg']):g} deg",
        "",
        "## 3. Branch Tracking Seed Convention",
        "",
        "`descendant k` is the continuation of the `k`-th sorted mode shape at",
        "`beta=0`, `mu=0`, for each fixed eta independently. Tracking uses",
        "adjacent-step analytic shape MAC from",
        "`scripts/lib/thickness_mismatch_mac_tracking.py`. Sorted position is",
        "metadata only.",
        "",
        "Low-MAC, low-margin, or large sorted-position-jump assignments are",
        "reported as tracking warnings and are not treated as branch renaming.",
        "",
        "## 4. Target Result: eta=0.5, beta=15 deg",
        "",
        f"- sorted 5 Lambda: `{float(target_main['sorted5_lambda']):.12g}`",
        f"- sorted 5 descendant identity: `{target_main['sorted5_descendant_id']}`",
        f"- descendant 6 sorted position: `{target_main['desc6_sorted_position']}`",
        f"- descendant 6 Lambda: `{float(target_main['desc6_lambda']):.12g}`",
        f"- sorted 5 equals descendant 6: `{'yes' if target_main['sorted5_equals_desc6'] else 'no'}`",
        "",
        "## 5. Eta=0 Comparison",
        "",
    ]
    if target_comparison is None:
        lines.append("The eta=0 comparison was not run.")
    else:
        lines.extend(
            [
                f"- sorted 5 Lambda at beta={float(target_comparison['beta_deg']):g} deg: "
                f"`{float(target_comparison['sorted5_lambda']):.12g}`",
                f"- sorted 5 descendant identity: `{target_comparison['sorted5_descendant_id']}`",
                f"- descendant 6 sorted position: `{target_comparison['desc6_sorted_position']}`",
                f"- descendant 6 Lambda: `{float(target_comparison['desc6_lambda']):.12g}`",
                f"- sorted 5 equals descendant 6: `{'yes' if target_comparison['sorted5_equals_desc6'] else 'no'}`",
                f"- eta=0 tracking warning rows: `{comparison_warnings['warning_rows']}`",
                "",
                "Thus the equal-thickness relation `sorted 5 = descendant 6` is",
                "confirmed at beta=15 in the eta=0 comparison, but it does not",
                "carry over to eta=0.5.",
            ]
        )

    lines.extend(
        [
            "",
            "## 6. Sorted 5 Identity Intervals",
            "",
        ]
    )
    lines.extend(
        markdown_table(
            ["beta start", "beta end", "sorted 5 descendant", "min gap 4-5", "min gap 5-6", "notes"],
            intervals,
        )
    )
    if comparison_intervals:
        lines.extend(["", "Eta=0 comparison intervals:", ""])
        lines.extend(
            markdown_table(
                ["beta start", "beta end", "sorted 5 descendant", "min gap 4-5", "min gap 5-6", "notes"],
                comparison_intervals,
            )
        )

    lines.extend(
        [
            "",
            "## 7. Nearest Sorted Gaps",
            "",
            f"- nearest sorted 4-5 gap: `{float(gaps['gap45_min_value']):.12g}` "
            f"at beta `{float(gaps['gap45_min_beta']):.10g}` deg",
            f"- nearest sorted 5-6 gap: `{float(gaps['gap56_min_value']):.12g}` "
            f"at beta `{float(gaps['gap56_min_beta']):.10g}` deg",
            "",
            "Closest local minima for sorted 4-5:",
            "",
        ]
    )
    lines.extend(markdown_table(["beta", "gap"], gap45_local))
    lines.extend(["", "Closest local minima for sorted 5-6:", ""])
    lines.extend(markdown_table(["beta", "gap"], gap56_local))

    lines.extend(
        [
            "",
            "Sorted-position changes among tracked descendants:",
            "",
        ]
    )
    lines.extend(sorted_position_change_lines(main_result))

    lines.extend(
        [
            "",
            "## 8. Previous Mode-Shape Task",
            "",
        ]
    )
    if correction_needed:
        lines.extend(
            [
                "At the target eta=0.5 point, sorted 5 is not descendant 6. Therefore",
                "the previous descendant-6 veering-shape plots are not the requested",
                "sorted-5 shapes for this parameter point.",
                "",
                "Recommended follow-up prompt:",
                "",
                "> Build sorted-5 shape plots at the six selected veering points.",
            ]
        )
    else:
        lines.append(
            "At the target eta=0.5 point, sorted 5 is descendant 6. The previous "
            "descendant-6 veering-shape plots are therefore not identity-wrong for "
            "the sorted-5 question at beta=15, mu=0 under this audit convention."
        )

    lines.extend(
        [
            "",
            "## 9. Tracking Warnings",
            "",
            f"- warning rows: {warnings['warning_rows']}",
            f"- low-MAC rows: {warnings['low_mac']}",
            f"- low-margin rows: {warnings['low_margin']}",
            f"- unresolved rows: {warnings['unresolved']}",
            f"- nearest-frequency/MAC disagreement rows: {warnings['frequency_mac_disagreement']}",
        ]
    )
    if comparison_warnings is not None:
        lines.extend(
            [
                "",
                "Eta=0 comparison warnings:",
                "",
                f"- warning rows: {comparison_warnings['warning_rows']}",
                f"- low-MAC rows: {comparison_warnings['low_mac']}",
                f"- low-margin rows: {comparison_warnings['low_margin']}",
                f"- unresolved rows: {comparison_warnings['unresolved']}",
                "- these warnings occur in the comparison case, not in the eta=0.5 target case.",
            ]
        )
    if main_result.warning_rows:
        lines.append("")
        lines.append("First warning rows:")
        for row in main_result.warning_rows[:12]:
            lines.append(
                f"- beta={float(row['beta_deg']):.6g}, desc={int(row['descendant_id'])}, "
                f"accepted position={int(row['current_sorted_position'])}, "
                f"candidate position={int(row['diagnostic_candidate_sorted_position'])}, "
                f"MAC={float(row['diagnostic_candidate_mac_to_previous']):.6g}, "
                f"warning={tracking_warning(row)}"
            )
        if len(main_result.warning_rows) > 12:
            lines.append(f"- plus {len(main_result.warning_rows) - 12} additional warning rows.")

    lines.extend(
        [
            "",
            "## 10. Limitations",
            "",
            "- This is an Euler-Bernoulli analytic thickness-mismatch audit only.",
            "- The root finder is the existing sign-scan/bisection path in",
            "  `formulas_thickness_mismatch.py`; no determinant entries or",
            "  unknown ordering were changed.",
            "- Tracking is adjacent-step shape MAC in the existing local-component",
            "  reconstruction convention.",
            "- If a target interval contains tracking warnings, use a localized",
            "  refined/subspace audit before turning that local assignment into a",
            "  physical statement.",
            "",
            "## Outputs",
            "",
            f"- main CSV: `{paths['csv'].relative_to(REPO_ROOT)}`",
            f"- summary CSV: `{paths['summary_csv'].relative_to(REPO_ROOT)}`",
            f"- main PNG: `{paths['png'].relative_to(REPO_ROOT)}`",
            f"- sorted-5 focus PNG: `{paths['focus_png'].relative_to(REPO_ROOT)}`",
        ]
    )
    if comparison_result is not None:
        comparison_paths = output_paths(comparison_result.eta)
        lines.extend(
            [
                f"- eta=0 comparison CSV: `{comparison_paths['csv'].relative_to(REPO_ROOT)}`",
                f"- eta=0 comparison summary CSV: `{comparison_paths['summary_csv'].relative_to(REPO_ROOT)}`",
                f"- eta=0 comparison PNG: `{comparison_paths['png'].relative_to(REPO_ROOT)}`",
                f"- eta=0 sorted-5 focus PNG: `{comparison_paths['focus_png'].relative_to(REPO_ROOT)}`",
            ]
        )
    lines.append("")

    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].write_text("\n".join(lines), encoding="utf-8")


def run_case(eta: float) -> tuple[BetaTrackingResult, list[dict[str, float | str]], dict[str, Path]]:
    paths = output_paths(float(eta))
    sorted_roots = solve_sorted_root_grid(eta=float(eta))
    tracking = track_descendants_by_beta(eta=float(eta), sorted_roots=sorted_roots)
    write_audit_csv(tracking, paths["csv"])
    summary_rows = write_summary_csv(tracking, paths["summary_csv"])
    plot_audit(tracking, paths["png"])
    plot_sorted5_focus(tracking, paths["focus_png"])
    target = result_at_beta(tracking)
    print(f"saved audit CSV: {paths['csv']}")
    print(f"saved summary CSV: {paths['summary_csv']}")
    print(f"saved audit PNG: {paths['png']}")
    print(f"saved sorted-5 focus PNG: {paths['focus_png']}")
    print(
        f"eta={float(eta):g}, beta={float(target['beta_deg']):g}: "
        f"sorted5 descendant={target['sorted5_descendant_id']}, "
        f"desc6 sorted position={target['desc6_sorted_position']}, "
        f"same={'yes' if target['sorted5_equals_desc6'] else 'no'}"
    )
    warnings = warning_summary(tracking)
    print(
        f"eta={float(eta):g}: warning rows={warnings['warning_rows']}, "
        f"low_mac={warnings['low_mac']}, low_margin={warnings['low_margin']}, "
        f"unresolved={warnings['unresolved']}"
    )
    return tracking, summary_rows, paths


def main() -> dict[str, object]:
    results: dict[float, BetaTrackingResult] = {}
    summaries: dict[float, list[dict[str, float | str]]] = {}
    paths_by_eta: dict[float, dict[str, Path]] = {}
    for eta in ETA_VALUES:
        tracking, summary_rows, paths = run_case(float(eta))
        results[float(eta)] = tracking
        summaries[float(eta)] = summary_rows
        paths_by_eta[float(eta)] = paths

    main_eta = float(ETA_VALUES[0])
    comparison_eta = float(ETA_VALUES[1]) if len(ETA_VALUES) > 1 else None
    write_report(
        main_result=results[main_eta],
        comparison_result=results.get(comparison_eta) if comparison_eta is not None else None,
        main_summary_rows=summaries[main_eta],
        comparison_summary_rows=summaries.get(comparison_eta) if comparison_eta is not None else None,
        paths=paths_by_eta[main_eta],
    )
    print(f"saved main report: {paths_by_eta[main_eta]['report']}")
    return {
        "main": results[main_eta],
        "comparison": results.get(comparison_eta) if comparison_eta is not None else None,
        "paths": paths_by_eta,
    }


if __name__ == "__main__":
    main()
