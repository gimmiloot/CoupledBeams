from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.compare_beta0_analytic_vs_fem import (  # noqa: E402
    build_params,
    fem_axial_fractions,
    fem_parameter_override,
)
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.sweep_grid_policy import (  # noqa: E402
    ANALYSIS_MU_STEP,
    analysis_mu_grid,
)


EXISTING_SWEEP_PATHS = {
    0.0: REPO_ROOT / "results" / "mu_sweep_beta0_four_radii_compare.csv",
    7.5: REPO_ROOT / "results" / "mu_sweep_beta7p5_four_radii_compare.csv",
    15.0: REPO_ROOT / "results" / "mu_sweep_beta15_four_radii_compare.csv",
}
RESULTS_DIR = REPO_ROOT / "results"
FLATTENING_METRICS_PATH = RESULTS_DIR / "flat_mu_branch_flattening_metrics.csv"
CANDIDATE_SUMMARY_PATH = RESULTS_DIR / "flat_mu_branch_candidates.csv"
OCCUPANCY_PATH = RESULTS_DIR / "flat_mu_slot_occupancy.csv"
REPRESENTATIVE_CASES_PATH = RESULTS_DIR / "flat_mu_representative_cases.csv"
AXIAL_RADIUS_PLOT_PATH = RESULTS_DIR / "flat_mu_axial_fraction_vs_radius.png"

RADIUS_VALUES = (0.005, 0.01, 0.015, 0.02)
TARGET_BETA_VALUES = (0.0, 7.5, 15.0)
BETA_CONTINUATION_VALUES = (0.0, 2.5, 5.0, 7.5, 10.0, 12.5, 15.0)
LOW_SLOTS = (1, 2, 3, 4, 5, 6)
SELECTED_MU_TARGETS = (0.0, 0.45, 0.9)
REPRESENTATIVE_RADII = (0.005, 0.02)
MU_VALUES = analysis_mu_grid()
N_TRACK = 20
N_SOLVE = 36
FREQ_WEIGHT = 0.03
AXIAL_THRESHOLD = 0.9
BENDING_THRESHOLD = 0.1
FLATTENING_RATIO_THRESHOLD = 0.65
MAX_REPRESENTATIVE_SLOTS_PER_BETA = 8
MODE_SHAPE_SCALE_FRACTION = 0.12


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def branch_slot(branch_id: str) -> int:
    return int(str(branch_id).split("_")[-1])


def reference_mode_type(axial_fraction: float) -> str:
    if axial_fraction >= AXIAL_THRESHOLD:
        return "axial"
    if axial_fraction <= BENDING_THRESHOLD:
        return "bending"
    return "mixed"


def content_label(axial_fraction: float) -> str:
    if axial_fraction >= AXIAL_THRESHOLD:
        return "axial-like"
    if axial_fraction <= BENDING_THRESHOLD:
        return "bending-dominated"
    return "mixed"


def solve_fem_modes(radius: float, mu: float, beta_deg: float, n_modes: int) -> tuple[np.ndarray, np.ndarray]:
    params = build_params(radius)
    with fem_parameter_override(params):
        omega, vecs = fem.fem_solve(mu=mu, beta_deg=beta_deg, n_modes=n_modes)
        return omega * fem.scale, vecs


def mac_matrix(prev_vecs: np.ndarray, cur_vecs: np.ndarray) -> np.ndarray:
    out = np.zeros((prev_vecs.shape[1], cur_vecs.shape[1]), dtype=float)
    for row in range(prev_vecs.shape[1]):
        for col in range(cur_vecs.shape[1]):
            out[row, col] = fem.mac_value(prev_vecs[:, row], cur_vecs[:, col])
    return out


def assign_by_mac_and_frequency(
    prev_vecs: np.ndarray,
    prev_freqs: np.ndarray,
    cur_vecs: np.ndarray,
    cur_freqs: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    mac = mac_matrix(prev_vecs, cur_vecs)
    freq_penalty = np.abs(prev_freqs[:, None] - cur_freqs[None, :]) / np.maximum(prev_freqs[:, None], 1.0)
    cost = (1.0 - mac) + FREQ_WEIGHT * freq_penalty
    rows, cols = linear_sum_assignment(cost)
    assignment = np.full(prev_vecs.shape[1], -1, dtype=int)
    assigned_mac = np.full(prev_vecs.shape[1], np.nan, dtype=float)
    for row, col in zip(rows, cols):
        assignment[row] = col
        assigned_mac[row] = mac[row, col]
    if np.any(assignment < 0):
        raise RuntimeError("Failed to assign all tracked descendants.")
    return assignment, assigned_mac


def free_dof_indices() -> np.ndarray:
    n1 = n2 = fem.N_ELEM
    ndof = 3 * (n1 + n2 + 1)
    bc = {0, 1, 2, 3 * (n1 + n2), 3 * (n1 + n2) + 1, 3 * (n1 + n2) + 2}
    return np.array(sorted(set(range(ndof)) - bc), dtype=int)


def expand_full_vector(vec: np.ndarray) -> np.ndarray:
    free = free_dof_indices()
    full = np.zeros(3 * (2 * fem.N_ELEM + 1), dtype=float)
    full[free] = np.asarray(vec, dtype=float)
    return full


def build_node_coordinates(mu: float, beta_deg: float) -> tuple[np.ndarray, np.ndarray]:
    n1 = n2 = fem.N_ELEM
    l1 = fem.ell * (1.0 - mu)
    l2 = fem.ell * (1.0 + mu)
    beta_rad = np.deg2rad(beta_deg)

    x = np.zeros(n1 + n2 + 1, dtype=float)
    y = np.zeros(n1 + n2 + 1, dtype=float)

    x[: n1 + 1] = np.linspace(0.0, l1, n1 + 1)
    s2 = np.linspace(0.0, l2, n2 + 1)
    x[n1:] = l1 + s2 * np.cos(beta_rad)
    y[n1:] = s2 * np.sin(beta_rad)
    return x, y


def local_mode_metrics(full_vec: np.ndarray, beta_deg: float) -> dict[str, float]:
    n1 = n2 = fem.N_ELEM
    beta_rad = np.deg2rad(beta_deg)
    cos_beta = float(np.cos(beta_rad))
    sin_beta = float(np.sin(beta_rad))

    total_sq = 0.0
    local_axial_sq = 0.0
    arm1_sq = 0.0
    arm2_sq = 0.0

    for node_idx in range(n1 + n2 + 1):
        u, v, theta = full_vec[3 * node_idx : 3 * node_idx + 3]
        total_sq += u * u + v * v + theta * theta
        if node_idx <= n1:
            u_local = u
        else:
            u_local = cos_beta * u + sin_beta * v
        local_axial_sq += u_local * u_local

        if node_idx < n1:
            arm1_sq += u * u + v * v + theta * theta
        elif node_idx > n1:
            arm2_sq += u * u + v * v + theta * theta

    localization_den = arm1_sq + arm2_sq
    return {
        "local_axial_fraction": local_axial_sq / total_sq if total_sq else np.nan,
        "arm1_fraction_excl_joint": arm1_sq / localization_den if localization_den else np.nan,
        "arm2_fraction_excl_joint": arm2_sq / localization_den if localization_den else np.nan,
    }


def nearest_mu_value(target: float) -> float:
    idx = int(np.argmin(np.abs(MU_VALUES - target)))
    return float(MU_VALUES[idx])


SELECTED_MU_VALUES = tuple(nearest_mu_value(mu) for mu in SELECTED_MU_TARGETS)


def load_existing_mu_sweep_rows() -> list[dict[str, float | int | str]]:
    out: list[dict[str, float | int | str]] = []
    for beta_deg, path in EXISTING_SWEEP_PATHS.items():
        rows = load_csv_rows(path)
        for row in rows:
            out.append(
                {
                    "beta": float(row.get("beta", beta_deg)),
                    "radius": float(row["radius"]),
                    "radius_mm": float(row["radius_mm"]),
                    "mu": float(row["mu"]),
                    "low_slot": branch_slot(row["analytic_branch_id"]),
                    "branch_id": str(row["analytic_branch_id"]),
                    "mode_type": str(row["mode_type"]),
                    "fem_hz": float(row["fem_hz"]),
                    "fem_axial_fraction": float(row["fem_axial_fraction"]),
                    "source_csv": str(row.get("source_csv", path)),
                }
            )
    return out


def compute_flattening_metrics(
    existing_rows: list[dict[str, float | int | str]],
) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]], dict[float, list[int]]]:
    grouped: dict[tuple[float, int, float], list[dict[str, float | int | str]]] = defaultdict(list)
    for row in existing_rows:
        grouped[(float(row["beta"]), int(row["low_slot"]), float(row["radius"]))].append(row)

    metric_rows: list[dict[str, float | int | str]] = []
    by_branch: dict[tuple[float, int], list[dict[str, float | int | str]]] = defaultdict(list)

    for (beta_deg, low_slot, radius), rows in sorted(grouped.items()):
        rows_sorted = sorted(rows, key=lambda item: float(item["mu"]))
        mu_values = np.asarray([float(item["mu"]) for item in rows_sorted], dtype=float)
        freq_hz = np.asarray([float(item["fem_hz"]) for item in rows_sorted], dtype=float)
        axial = np.asarray([float(item["fem_axial_fraction"]) for item in rows_sorted], dtype=float)
        spread_hz = float(np.max(freq_hz) - np.min(freq_hz))
        mean_abs_hz = float(np.mean(np.abs(freq_hz)))
        relative_spread = spread_hz / mean_abs_hz if mean_abs_hz else np.nan
        max_abs_df_dmu = float(np.max(np.abs(np.gradient(freq_hz, mu_values))))

        metric_row = {
            "beta": float(beta_deg),
            "low_slot": int(low_slot),
            "radius": float(radius),
            "radius_mm": float(radius * 1e3),
            "frequency_hz_min": float(np.min(freq_hz)),
            "frequency_hz_max": float(np.max(freq_hz)),
            "spread_hz": spread_hz,
            "relative_spread": float(relative_spread),
            "max_abs_df_dmu_hz_per_unit": max_abs_df_dmu,
            "axial_fraction_min": float(np.min(axial)),
            "axial_fraction_max": float(np.max(axial)),
            "axial_fraction_mu0": float(axial[0]),
            "axial_fraction_mu_mid": float(axial[len(axial) // 2]),
            "axial_fraction_mu_last": float(axial[-1]),
        }
        metric_rows.append(metric_row)
        by_branch[(float(beta_deg), int(low_slot))].append(metric_row)

    candidate_rows: list[dict[str, float | int | str]] = []
    selected_by_beta: dict[float, list[int]] = defaultdict(list)
    for (beta_deg, low_slot), rows in sorted(by_branch.items()):
        rows_sorted = sorted(rows, key=lambda item: float(item["radius"]))
        first = rows_sorted[0]
        last = rows_sorted[-1]
        spread_ratio = float(last["spread_hz"]) / float(first["spread_hz"]) if float(first["spread_hz"]) else np.nan
        relative_ratio = (
            float(last["relative_spread"]) / float(first["relative_spread"])
            if float(first["relative_spread"])
            else np.nan
        )
        slope_ratio = (
            float(last["max_abs_df_dmu_hz_per_unit"]) / float(first["max_abs_df_dmu_hz_per_unit"])
            if float(first["max_abs_df_dmu_hz_per_unit"])
            else np.nan
        )
        is_flattening = beta_deg > 0.0 and np.isfinite(spread_ratio) and spread_ratio < FLATTENING_RATIO_THRESHOLD

        candidate_rows.append(
            {
                "beta": float(beta_deg),
                "low_slot": int(low_slot),
                "spread_ratio_rmax_over_rmin": spread_ratio,
                "relative_spread_ratio_rmax_over_rmin": relative_ratio,
                "slope_ratio_rmax_over_rmin": slope_ratio,
                "selected_as_flattening_candidate": int(is_flattening),
            }
        )
        if is_flattening:
            selected_by_beta[float(beta_deg)].append(int(low_slot))

    trimmed_selected: dict[float, list[int]] = {}
    for beta_deg, slots in selected_by_beta.items():
        ranked = sorted(
            (row for row in candidate_rows if float(row["beta"]) == beta_deg and int(row["selected_as_flattening_candidate"]) == 1),
            key=lambda row: float(row["spread_ratio_rmax_over_rmin"]),
        )
        trimmed_selected[beta_deg] = [int(row["low_slot"]) for row in ranked[:MAX_REPRESENTATIVE_SLOTS_PER_BETA]]

    return metric_rows, candidate_rows, trimmed_selected


def build_branch_metadata(seed_freqs: np.ndarray, seed_fracs: np.ndarray) -> list[dict[str, float | int | str]]:
    counters = {"bending": 0, "mixed": 0, "axial": 0}
    metadata: list[dict[str, float | int | str]] = []
    for ref_index in range(N_TRACK):
        ref_type = reference_mode_type(float(seed_fracs[ref_index]))
        counters[ref_type] += 1
        branch_id = f"{ref_type}_desc_{counters[ref_type]:02d}"
        metadata.append(
            {
                "branch_id": branch_id,
                "reference_mode_index_at_beta0": ref_index + 1,
                "reference_mode_type_at_beta0": ref_type,
                "reference_frequency_hz_at_beta0": float(seed_freqs[ref_index]),
                "reference_axial_fraction_at_beta0": float(seed_fracs[ref_index]),
            }
        )
    return metadata


def track_descendants_for_radius(radius: float) -> dict[str, object]:
    seed_freqs, seed_vecs = solve_fem_modes(radius=radius, mu=0.0, beta_deg=0.0, n_modes=N_SOLVE)
    seed_fracs = fem_axial_fractions(seed_vecs)
    metadata = build_branch_metadata(seed_freqs=seed_freqs[:N_TRACK], seed_fracs=seed_fracs[:N_TRACK])

    beta_seed_vecs: dict[float, np.ndarray] = {0.0: seed_vecs[:, :N_TRACK].copy()}
    beta_seed_freqs: dict[float, np.ndarray] = {0.0: seed_freqs[:N_TRACK].copy()}
    beta_seed_index: dict[float, np.ndarray] = {0.0: np.arange(1, N_TRACK + 1, dtype=int)}
    beta_seed_frac: dict[float, np.ndarray] = {0.0: seed_fracs[:N_TRACK].copy()}

    prev_vecs = beta_seed_vecs[0.0]
    prev_freqs = beta_seed_freqs[0.0]

    for beta_deg in BETA_CONTINUATION_VALUES[1:]:
        cur_freqs, cur_vecs = solve_fem_modes(radius=radius, mu=0.0, beta_deg=beta_deg, n_modes=N_SOLVE)
        cur_fracs = fem_axial_fractions(cur_vecs)
        assign, _ = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        beta_seed_vecs[float(beta_deg)] = cur_vecs[:, assign].copy()
        beta_seed_freqs[float(beta_deg)] = cur_freqs[assign].copy()
        beta_seed_index[float(beta_deg)] = assign + 1
        beta_seed_frac[float(beta_deg)] = cur_fracs[assign].copy()
        prev_vecs = beta_seed_vecs[float(beta_deg)]
        prev_freqs = beta_seed_freqs[float(beta_deg)]

    beta_tracks: dict[float, dict[str, np.ndarray]] = {}
    for beta_deg in TARGET_BETA_VALUES:
        tracked_freqs = np.full((len(MU_VALUES), N_TRACK), np.nan, dtype=float)
        tracked_index = np.full((len(MU_VALUES), N_TRACK), -1, dtype=int)
        tracked_frac = np.full((len(MU_VALUES), N_TRACK), np.nan, dtype=float)
        tracked_mac = np.full((len(MU_VALUES), N_TRACK), np.nan, dtype=float)

        tracked_freqs[0] = beta_seed_freqs[float(beta_deg)]
        tracked_index[0] = beta_seed_index[float(beta_deg)]
        tracked_frac[0] = beta_seed_frac[float(beta_deg)]

        prev_vecs = beta_seed_vecs[float(beta_deg)].copy()
        prev_freqs = beta_seed_freqs[float(beta_deg)].copy()

        for mu_pos in range(1, len(MU_VALUES)):
            mu_value = float(MU_VALUES[mu_pos])
            cur_freqs, cur_vecs = solve_fem_modes(radius=radius, mu=mu_value, beta_deg=beta_deg, n_modes=N_SOLVE)
            cur_fracs = fem_axial_fractions(cur_vecs)
            assign, mac_values = assign_by_mac_and_frequency(
                prev_vecs=prev_vecs,
                prev_freqs=prev_freqs,
                cur_vecs=cur_vecs,
                cur_freqs=cur_freqs,
            )
            tracked_freqs[mu_pos] = cur_freqs[assign]
            tracked_index[mu_pos] = assign + 1
            tracked_frac[mu_pos] = cur_fracs[assign]
            tracked_mac[mu_pos] = mac_values
            prev_vecs = cur_vecs[:, assign]
            prev_freqs = cur_freqs[assign]

        beta_tracks[float(beta_deg)] = {
            "frequency_hz": tracked_freqs,
            "current_sorted_index": tracked_index,
            "axial_fraction": tracked_frac,
            "assignment_mac": tracked_mac,
        }

    return {
        "radius": float(radius),
        "metadata": metadata,
        "tracks": beta_tracks,
    }


def build_slot_occupancy_rows(
    tracked_by_radius: dict[float, dict[str, object]],
    selected_slots_by_beta: dict[float, list[int]],
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []

    for radius, tracking in tracked_by_radius.items():
        metadata = list(tracking["metadata"])
        tracks = dict(tracking["tracks"])

        for beta_deg in TARGET_BETA_VALUES:
            slots = LOW_SLOTS if beta_deg == 0.0 else tuple(selected_slots_by_beta.get(beta_deg, []))
            if not slots:
                continue
            sorted_index = np.asarray(tracks[float(beta_deg)]["current_sorted_index"], dtype=int)
            frequency_hz = np.asarray(tracks[float(beta_deg)]["frequency_hz"], dtype=float)
            axial_fraction = np.asarray(tracks[float(beta_deg)]["axial_fraction"], dtype=float)
            assignment_mac = np.asarray(tracks[float(beta_deg)]["assignment_mac"], dtype=float)

            for mu_pos, mu_value in enumerate(MU_VALUES):
                for low_slot in slots:
                    matches = np.where(sorted_index[mu_pos] == int(low_slot))[0]
                    if len(matches) == 1:
                        branch_pos = int(matches[0])
                        branch_meta = metadata[branch_pos]
                        rows.append(
                            {
                                "beta": float(beta_deg),
                                "radius": float(radius),
                                "radius_mm": float(radius * 1e3),
                                "mu": float(mu_value),
                                "low_slot": int(low_slot),
                                "occupant_branch_id": str(branch_meta["branch_id"]),
                                "occupant_reference_type_at_beta0": str(branch_meta["reference_mode_type_at_beta0"]),
                                "occupant_reference_mode_index_at_beta0": int(branch_meta["reference_mode_index_at_beta0"]),
                                "occupant_frequency_hz": float(frequency_hz[mu_pos, branch_pos]),
                                "occupant_axial_fraction": float(axial_fraction[mu_pos, branch_pos]),
                                "occupant_assignment_mac": float(assignment_mac[mu_pos, branch_pos]),
                            }
                        )
                    else:
                        rows.append(
                            {
                                "beta": float(beta_deg),
                                "radius": float(radius),
                                "radius_mm": float(radius * 1e3),
                                "mu": float(mu_value),
                                "low_slot": int(low_slot),
                                "occupant_branch_id": "unresolved",
                                "occupant_reference_type_at_beta0": "unresolved",
                                "occupant_reference_mode_index_at_beta0": -1,
                                "occupant_frequency_hz": np.nan,
                                "occupant_axial_fraction": np.nan,
                                "occupant_assignment_mac": np.nan,
                            }
                        )
    return rows


def summarize_slot_identity(
    occupancy_rows: list[dict[str, float | int | str]],
) -> dict[tuple[float, float, int], dict[str, float | int | str]]:
    grouped: dict[tuple[float, float, int], list[dict[str, float | int | str]]] = defaultdict(list)
    for row in occupancy_rows:
        grouped[(float(row["beta"]), float(row["radius"]), int(row["low_slot"]))].append(row)

    summary: dict[tuple[float, float, int], dict[str, float | int | str]] = {}
    for key, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda item: float(item["mu"]))
        seq = [str(row["occupant_branch_id"]) for row in rows_sorted if str(row["occupant_branch_id"]) != "unresolved"]
        switches = 0
        if len(seq) >= 2:
            switches = sum(1 for left, right in zip(seq[:-1], seq[1:]) if left != right)
        counts = Counter(seq)
        dominant_branch = counts.most_common(1)[0][0] if counts else "unresolved"
        dominant_fraction = counts[dominant_branch] / len(seq) if seq else np.nan
        summary[key] = {
            "identity_switches_over_mu": int(switches),
            "dominant_occupant_branch_id": dominant_branch,
            "dominant_occupant_fraction_over_mu": float(dominant_fraction) if np.isfinite(dominant_fraction) else np.nan,
        }
    return summary


def build_representative_case_rows(
    occupancy_rows: list[dict[str, float | int | str]],
    identity_summary: dict[tuple[float, float, int], dict[str, float | int | str]],
    selected_slots_by_beta: dict[float, list[int]],
) -> list[dict[str, float | int | str]]:
    occupancy_lookup = {
        (float(row["beta"]), float(row["radius"]), float(row["mu"]), int(row["low_slot"])): row
        for row in occupancy_rows
    }

    rows: list[dict[str, float | int | str]] = []
    for beta_deg in (7.5, 15.0):
        for low_slot in selected_slots_by_beta.get(beta_deg, []):
            for radius in REPRESENTATIVE_RADII:
                identity_info = identity_summary[(float(beta_deg), float(radius), int(low_slot))]
                for mu_value in SELECTED_MU_VALUES:
                    occupancy = occupancy_lookup[(float(beta_deg), float(radius), float(mu_value), int(low_slot))]
                    freq_hz, vecs = solve_fem_modes(radius=radius, mu=mu_value, beta_deg=beta_deg, n_modes=max(N_SOLVE, low_slot))
                    sorted_pos = int(low_slot) - 1
                    phi = vecs[:, sorted_pos]
                    global_axial_fraction = float(fem_axial_fractions(vecs[:, [sorted_pos]])[0])
                    full_vec = expand_full_vector(phi)
                    local_metrics = local_mode_metrics(full_vec=full_vec, beta_deg=beta_deg)
                    long_arm_fraction = float(local_metrics["arm2_fraction_excl_joint"]) if mu_value >= 0.0 else np.nan

                    rows.append(
                        {
                            "beta": float(beta_deg),
                            "radius": float(radius),
                            "radius_mm": float(radius * 1e3),
                            "mu": float(mu_value),
                            "low_slot": int(low_slot),
                            "occupant_branch_id": str(occupancy["occupant_branch_id"]),
                            "occupant_reference_type_at_beta0": str(occupancy["occupant_reference_type_at_beta0"]),
                            "occupant_reference_mode_index_at_beta0": int(occupancy["occupant_reference_mode_index_at_beta0"]),
                            "frequency_hz": float(freq_hz[sorted_pos]),
                            "current_sorted_index": int(low_slot),
                            "global_axial_fraction": global_axial_fraction,
                            "local_axial_fraction": float(local_metrics["local_axial_fraction"]),
                            "content_label_by_local_axial_fraction": content_label(float(local_metrics["local_axial_fraction"])),
                            "arm1_fraction_excl_joint": float(local_metrics["arm1_fraction_excl_joint"]),
                            "arm2_fraction_excl_joint": float(local_metrics["arm2_fraction_excl_joint"]),
                            "long_arm_fraction_excl_joint": long_arm_fraction,
                            "identity_switches_over_mu": int(identity_info["identity_switches_over_mu"]),
                            "dominant_occupant_branch_id_over_mu": str(identity_info["dominant_occupant_branch_id"]),
                            "dominant_occupant_fraction_over_mu": float(identity_info["dominant_occupant_fraction_over_mu"]),
                        }
                    )
    return rows


def build_axial_fraction_radius_plot(
    existing_rows: list[dict[str, float | int | str]],
    selected_slots_by_beta: dict[float, list[int]],
) -> None:
    grouped: dict[tuple[float, float, int, float], float] = {}
    for row in existing_rows:
        grouped[(float(row["beta"]), float(row["radius"]), int(row["low_slot"]), float(row["mu"]))] = float(row["fem_axial_fraction"])

    fig, axes = plt.subplots(2, len(SELECTED_MU_VALUES), figsize=(16.5, 7.6), sharex=True, sharey=True, squeeze=False)
    color_map = {3: "#4c78a8", 4: "#f58518", 5: "#54a24b", 6: "#e45756"}

    for row_pos, beta_deg in enumerate((7.5, 15.0)):
        selected_slots = selected_slots_by_beta.get(beta_deg, [])
        for col_pos, mu_value in enumerate(SELECTED_MU_VALUES):
            ax = axes[row_pos, col_pos]
            for low_slot in selected_slots:
                y_values = np.asarray(
                    [grouped[(beta_deg, radius, low_slot, mu_value)] for radius in RADIUS_VALUES],
                    dtype=float,
                )
                ax.plot(
                    np.asarray(RADIUS_VALUES, dtype=float) * 1e3,
                    y_values,
                    marker="o",
                    linewidth=2.0,
                    color=color_map.get(low_slot, None),
                    label=f"slot {low_slot}",
                )

            ax.set_yscale("log")
            ax.grid(True, alpha=0.3)
            ax.set_title(f"beta = {beta_deg:g} deg, mu = {mu_value:.6f}")
            if row_pos == 1:
                ax.set_xlabel("radius [mm]")
            if col_pos == 0:
                ax.set_ylabel("existing fem_axial_fraction")
            if selected_slots:
                ax.legend(loc="best", fontsize=8)

    fig.suptitle("Low-slot axial_fraction vs radius for the selected flattening branches", y=0.98)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(AXIAL_RADIUS_PLOT_PATH, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_mode_shape_grid(
    beta_deg: float,
    low_slot: int,
    representative_rows: list[dict[str, float | int | str]],
) -> None:
    rows = [row for row in representative_rows if float(row["beta"]) == beta_deg and int(row["low_slot"]) == low_slot]
    if not rows:
        return

    rows_sorted = sorted(rows, key=lambda item: (float(item["radius"]), float(item["mu"])))
    fig, axes = plt.subplots(2, len(SELECTED_MU_VALUES), figsize=(16.5, 8.6), squeeze=False)

    for row in rows_sorted:
        radius = float(row["radius"])
        mu_value = float(row["mu"])
        row_pos = REPRESENTATIVE_RADII.index(radius)
        col_pos = SELECTED_MU_VALUES.index(mu_value)
        ax = axes[row_pos, col_pos]

        _, vecs = solve_fem_modes(radius=radius, mu=mu_value, beta_deg=beta_deg, n_modes=max(N_SOLVE, low_slot))
        phi = vecs[:, int(low_slot) - 1]
        full_vec = expand_full_vector(phi)
        x_base, y_base = build_node_coordinates(mu=mu_value, beta_deg=beta_deg)
        u = full_vec[0::3]
        v = full_vec[1::3]
        max_disp = float(np.max(np.sqrt(u * u + v * v)))
        total_length = float(fem.ell * (1.0 - mu_value) + fem.ell * (1.0 + mu_value))
        scale = MODE_SHAPE_SCALE_FRACTION * total_length / max(max_disp, 1e-12)
        x_def = x_base + scale * u
        y_def = y_base + scale * v
        joint = fem.N_ELEM

        ax.plot(x_base[: joint + 1], y_base[: joint + 1], color="0.7", linestyle="--", linewidth=1.1)
        ax.plot(x_base[joint:], y_base[joint:], color="0.7", linestyle="--", linewidth=1.1)
        ax.plot(x_def[: joint + 1], y_def[: joint + 1], color="#1f77b4", linewidth=2.0)
        ax.plot(x_def[joint:], y_def[joint:], color="#ff7f0e", linewidth=2.0)
        ax.scatter([x_base[joint]], [y_base[joint]], color="black", s=14, zorder=5)

        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.25)
        ax.set_title(
            f"r = {radius * 1e3:.1f} mm, mu = {mu_value:.6f}\n"
            f"{row['occupant_branch_id']} | {row['content_label_by_local_axial_fraction']}"
        )
        ax.text(
            0.02,
            0.02,
            (
                f"f = {float(row['frequency_hz']):.2f} Hz\n"
                f"global ax = {float(row['global_axial_fraction']):.4f}\n"
                f"local ax = {float(row['local_axial_fraction']):.4f}\n"
                f"long-arm frac = {float(row['long_arm_fraction_excl_joint']):.3f}\n"
                f"switches = {int(row['identity_switches_over_mu'])}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.92},
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    fig.suptitle(
        f"Representative sorted-slot mode shapes: beta = {beta_deg:g} deg, low slot = {low_slot}",
        y=0.98,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out_path = RESULTS_DIR / f"flat_mu_mode_shapes_beta{str(beta_deg).replace('.', 'p')}_slot{low_slot}.png"
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    existing_rows = load_existing_mu_sweep_rows()
    metric_rows, candidate_rows, selected_slots_by_beta = compute_flattening_metrics(existing_rows=existing_rows)
    write_csv_rows(FLATTENING_METRICS_PATH, metric_rows)
    write_csv_rows(CANDIDATE_SUMMARY_PATH, candidate_rows)

    tracked_by_radius = {
        float(radius): track_descendants_for_radius(radius=float(radius))
        for radius in REPRESENTATIVE_RADII
    }
    occupancy_rows = build_slot_occupancy_rows(
        tracked_by_radius=tracked_by_radius,
        selected_slots_by_beta=selected_slots_by_beta,
    )
    write_csv_rows(OCCUPANCY_PATH, occupancy_rows)

    identity_summary = summarize_slot_identity(occupancy_rows=occupancy_rows)
    representative_rows = build_representative_case_rows(
        occupancy_rows=occupancy_rows,
        identity_summary=identity_summary,
        selected_slots_by_beta=selected_slots_by_beta,
    )
    write_csv_rows(REPRESENTATIVE_CASES_PATH, representative_rows)

    build_axial_fraction_radius_plot(
        existing_rows=existing_rows,
        selected_slots_by_beta=selected_slots_by_beta,
    )
    for beta_deg, slots in selected_slots_by_beta.items():
        for low_slot in slots:
            plot_mode_shape_grid(
                beta_deg=float(beta_deg),
                low_slot=int(low_slot),
                representative_rows=representative_rows,
            )

    print(f"saved: {FLATTENING_METRICS_PATH}")
    print(f"saved: {CANDIDATE_SUMMARY_PATH}")
    print(f"saved: {OCCUPANCY_PATH}")
    print(f"saved: {REPRESENTATIVE_CASES_PATH}")
    print(f"saved: {AXIAL_RADIUS_PLOT_PATH}")
    print(f"analysis mu base step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    for beta_deg, slots in sorted(selected_slots_by_beta.items()):
        print(f"beta={beta_deg:g} deg selected slots: {', '.join(str(slot) for slot in slots)}")


if __name__ == "__main__":
    main()
