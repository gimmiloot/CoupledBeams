from __future__ import annotations

import csv
from math import ceil
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.optimize import linear_sum_assignment


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.compare_beta0_analytic_vs_fem import (  # noqa: E402
    BASE_PARAMS,
    fem_axial_fractions,
    fem_parameter_override,
)
from my_project.fem import python_fem as fem  # noqa: E402


BETA_VALUES = np.arange(0.0, 15.0 + 0.25, 0.5)
MU_VALUES = np.linspace(0.0, 0.9, 61)
N_TRACK = 30
N_SOLVE = 50
FREQ_WEIGHT = 0.03
AXIAL_THRESHOLD = 0.9
BENDING_THRESHOLD = 0.1
LOW_SPECTRUM_CUTOFF = 15
SELECTED_BETAS = (0.0, 2.0, 5.0, 10.0, 15.0)
SELECTED_MUS = (0.0, 0.15, 0.3, 0.45, 0.6, 0.9)
SELECTED_AXIAL_IDS = ("axial_desc_01", "axial_desc_02", "axial_desc_03")

RESULTS_DIR = REPO_ROOT / "results"
BRANCH_GRID_PATH = RESULTS_DIR / "branchwise_fem_grid.csv"
BRANCH_SUMMARY_PATH = RESULTS_DIR / "branchwise_fem_branch_summary.csv"
EXTREMA_PATH = RESULTS_DIR / "branchwise_fem_extrema.csv"
PAIR_GAPS_PATH = RESULTS_DIR / "branchwise_fem_near_crossings.csv"
POINT_GAPS_PATH = RESULTS_DIR / "branchwise_fem_min_gap_map.csv"
FREQ_MU_PLOT_PATH = RESULTS_DIR / "branchwise_freq_vs_mu_by_beta.png"
FREQ_BETA_PLOT_PATH = RESULTS_DIR / "branchwise_freq_vs_beta_by_mu.png"
SORTED_SPECTRUM_MU_PLOT_PATH = RESULTS_DIR / "branchwise_sorted_spectrum_vs_mu.png"
SORTED_SPECTRUM_BETA_PLOT_PATH = RESULTS_DIR / "branchwise_sorted_spectrum_vs_beta.png"
AXIAL_HEATMAP_PATH = RESULTS_DIR / "branchwise_axial_fraction_heatmaps.png"
MIN_GAP_PLOT_PATH = RESULTS_DIR / "branchwise_min_gap_map.png"


def solve_fem_modes(mu: float, beta_deg: float, n_modes: int) -> tuple[np.ndarray, np.ndarray]:
    with fem_parameter_override(BASE_PARAMS):
        omega, vecs = fem.fem_solve(mu=mu, beta_deg=beta_deg, n_modes=n_modes)
        return omega * fem.scale, vecs


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def reference_mode_type(axial_fraction: float) -> str:
    if axial_fraction >= AXIAL_THRESHOLD:
        return "axial"
    if axial_fraction <= BENDING_THRESHOLD:
        return "bending"
    return "mixed"


def descendant_label(reference_type: str) -> str:
    return f"{reference_type}-descendant"


def branch_palette(branch_ids: list[str]) -> dict[str, tuple[float, float, float, float]]:
    cmap = plt.get_cmap("turbo")
    colors = cmap(np.linspace(0.02, 0.98, len(branch_ids)))
    return {branch_id: tuple(color) for branch_id, color in zip(branch_ids, colors)}


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
        raise RuntimeError("Failed to assign all tracked branches.")
    return assignment, assigned_mac


def make_branch_metadata(seed_freqs: np.ndarray, seed_fracs: np.ndarray) -> list[dict[str, float | int | str]]:
    counters = {"bending": 0, "axial": 0, "mixed": 0}
    metadata = []
    for ref_index in range(N_TRACK):
        ref_type = reference_mode_type(float(seed_fracs[ref_index]))
        counters[ref_type] += 1
        branch_id = f"{ref_type}_desc_{counters[ref_type]:02d}"
        metadata.append(
            {
                "branch_id": branch_id,
                "reference_mode_index_at_beta0": ref_index + 1,
                "reference_mode_type_at_beta0": ref_type,
                "reference_descendant_label": descendant_label(ref_type),
                "reference_frequency_hz_at_beta0": float(seed_freqs[ref_index]),
                "reference_axial_fraction_at_beta0": float(seed_fracs[ref_index]),
            }
        )
    return metadata


def track_branch_grid() -> dict[str, object]:
    beta_count = len(BETA_VALUES)
    mu_count = len(MU_VALUES)

    full_freqs = np.full((beta_count, mu_count, N_SOLVE), np.nan, dtype=float)
    tracked_freqs = np.full((beta_count, mu_count, N_TRACK), np.nan, dtype=float)
    tracked_index = np.full((beta_count, mu_count, N_TRACK), -1, dtype=int)
    tracked_frac = np.full((beta_count, mu_count, N_TRACK), np.nan, dtype=float)
    tracked_mac = np.full((beta_count, mu_count, N_TRACK), np.nan, dtype=float)

    beta_seed_vecs: list[np.ndarray | None] = [None] * beta_count
    beta_seed_freqs: list[np.ndarray | None] = [None] * beta_count

    seed_freqs, seed_vecs = solve_fem_modes(mu=float(MU_VALUES[0]), beta_deg=float(BETA_VALUES[0]), n_modes=N_SOLVE)
    seed_fracs = fem_axial_fractions(seed_vecs)
    metadata = make_branch_metadata(seed_freqs=seed_freqs[:N_TRACK], seed_fracs=seed_fracs[:N_TRACK])

    full_freqs[0, 0] = seed_freqs
    tracked_freqs[0, 0] = seed_freqs[:N_TRACK]
    tracked_index[0, 0] = np.arange(1, N_TRACK + 1, dtype=int)
    tracked_frac[0, 0] = seed_fracs[:N_TRACK]
    beta_seed_vecs[0] = seed_vecs[:, :N_TRACK].copy()
    beta_seed_freqs[0] = seed_freqs[:N_TRACK].copy()

    prev_beta_vecs = beta_seed_vecs[0]
    prev_beta_freqs = beta_seed_freqs[0]

    for beta_pos in range(1, beta_count):
        beta_deg = float(BETA_VALUES[beta_pos])
        cur_freqs, cur_vecs = solve_fem_modes(mu=float(MU_VALUES[0]), beta_deg=beta_deg, n_modes=N_SOLVE)
        cur_fracs = fem_axial_fractions(cur_vecs)
        assign, mac_values = assign_by_mac_and_frequency(
            prev_vecs=prev_beta_vecs,
            prev_freqs=prev_beta_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        full_freqs[beta_pos, 0] = cur_freqs
        tracked_freqs[beta_pos, 0] = cur_freqs[assign]
        tracked_index[beta_pos, 0] = assign + 1
        tracked_frac[beta_pos, 0] = cur_fracs[assign]
        tracked_mac[beta_pos, 0] = mac_values
        beta_seed_vecs[beta_pos] = cur_vecs[:, assign].copy()
        beta_seed_freqs[beta_pos] = cur_freqs[assign].copy()
        prev_beta_vecs = beta_seed_vecs[beta_pos]
        prev_beta_freqs = beta_seed_freqs[beta_pos]

    for beta_pos, beta_deg in enumerate(BETA_VALUES):
        prev_vecs = np.asarray(beta_seed_vecs[beta_pos], dtype=float)
        prev_freqs = np.asarray(beta_seed_freqs[beta_pos], dtype=float)
        for mu_pos in range(1, mu_count):
            mu = float(MU_VALUES[mu_pos])
            cur_freqs, cur_vecs = solve_fem_modes(mu=mu, beta_deg=float(beta_deg), n_modes=N_SOLVE)
            cur_fracs = fem_axial_fractions(cur_vecs)
            assign, mac_values = assign_by_mac_and_frequency(
                prev_vecs=prev_vecs,
                prev_freqs=prev_freqs,
                cur_vecs=cur_vecs,
                cur_freqs=cur_freqs,
            )
            full_freqs[beta_pos, mu_pos] = cur_freqs
            tracked_freqs[beta_pos, mu_pos] = cur_freqs[assign]
            tracked_index[beta_pos, mu_pos] = assign + 1
            tracked_frac[beta_pos, mu_pos] = cur_fracs[assign]
            tracked_mac[beta_pos, mu_pos] = mac_values
            prev_vecs = cur_vecs[:, assign]
            prev_freqs = cur_freqs[assign]

    return {
        "metadata": metadata,
        "full_freqs": full_freqs,
        "tracked_freqs": tracked_freqs,
        "tracked_index": tracked_index,
        "tracked_frac": tracked_frac,
        "tracked_mac": tracked_mac,
    }


def compute_neighbor_gaps(
    full_freqs: np.ndarray,
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lower_gap = np.full_like(tracked_freqs, np.nan)
    upper_gap = np.full_like(tracked_freqs, np.nan)
    min_gap = np.full_like(tracked_freqs, np.nan)

    for beta_pos in range(tracked_freqs.shape[0]):
        for mu_pos in range(tracked_freqs.shape[1]):
            spectrum = full_freqs[beta_pos, mu_pos]
            for branch_pos in range(tracked_freqs.shape[2]):
                order = int(tracked_index[beta_pos, mu_pos, branch_pos])
                freq = float(tracked_freqs[beta_pos, mu_pos, branch_pos])
                if order > 1:
                    lower = (freq - spectrum[order - 2]) / freq
                    lower_gap[beta_pos, mu_pos, branch_pos] = lower
                if order < len(spectrum):
                    upper = (spectrum[order] - freq) / freq
                    upper_gap[beta_pos, mu_pos, branch_pos] = upper
                gaps = [value for value in (lower_gap[beta_pos, mu_pos, branch_pos], upper_gap[beta_pos, mu_pos, branch_pos]) if np.isfinite(value)]
                if gaps:
                    min_gap[beta_pos, mu_pos, branch_pos] = min(gaps)

    return lower_gap, upper_gap, min_gap


def branch_slice_sign_changes(values: np.ndarray) -> int:
    scale = max(float(np.nanmax(np.abs(values))), 1.0)
    tol = 1e-3 * scale
    signs = np.where(values > tol, 1, np.where(values < -tol, -1, 0))
    compact = [int(sign) for sign in signs if sign != 0]
    if len(compact) < 2:
        return 0
    return int(sum(1 for left, right in zip(compact[:-1], compact[1:]) if left != right))


def monotonicity_label(values: np.ndarray) -> str:
    scale = max(float(np.nanmax(np.abs(values))), 1.0)
    tol = 1e-3 * scale
    if np.all(values >= -tol):
        return "increasing"
    if np.all(values <= tol):
        return "decreasing"
    return "mixed"


def collect_extrema_rows(
    metadata: list[dict[str, float | int | str]],
    tracked_freqs: np.ndarray,
    df_dmu: np.ndarray,
    df_dbeta: np.ndarray,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []

    for branch_pos, branch_meta in enumerate(metadata):
        branch_id = str(branch_meta["branch_id"])
        reference_type = str(branch_meta["reference_mode_type_at_beta0"])

        for beta_pos, beta_deg in enumerate(BETA_VALUES):
            deriv = df_dmu[beta_pos, :, branch_pos]
            freq = tracked_freqs[beta_pos, :, branch_pos]
            for mu_pos in range(1, len(MU_VALUES) - 1):
                left = deriv[mu_pos - 1]
                right = deriv[mu_pos + 1]
                scale = max(abs(left), abs(right), 1.0)
                tol = 1e-3 * scale
                if left > tol and right < -tol:
                    rows.append(
                        {
                            "branch_id": branch_id,
                            "reference_mode_type_at_beta0": reference_type,
                            "sweep_parameter": "mu",
                            "fixed_parameter_name": "beta",
                            "fixed_parameter_value": float(beta_deg),
                            "extremum_parameter_value": float(MU_VALUES[mu_pos]),
                            "frequency_hz": float(freq[mu_pos]),
                            "extremum_kind": "local_max",
                        }
                    )
                if left < -tol and right > tol:
                    rows.append(
                        {
                            "branch_id": branch_id,
                            "reference_mode_type_at_beta0": reference_type,
                            "sweep_parameter": "mu",
                            "fixed_parameter_name": "beta",
                            "fixed_parameter_value": float(beta_deg),
                            "extremum_parameter_value": float(MU_VALUES[mu_pos]),
                            "frequency_hz": float(freq[mu_pos]),
                            "extremum_kind": "local_min",
                        }
                    )

        for mu_pos, mu in enumerate(MU_VALUES):
            deriv = df_dbeta[:, mu_pos, branch_pos]
            freq = tracked_freqs[:, mu_pos, branch_pos]
            for beta_pos in range(1, len(BETA_VALUES) - 1):
                left = deriv[beta_pos - 1]
                right = deriv[beta_pos + 1]
                scale = max(abs(left), abs(right), 1.0)
                tol = 1e-3 * scale
                if left > tol and right < -tol:
                    rows.append(
                        {
                            "branch_id": branch_id,
                            "reference_mode_type_at_beta0": reference_type,
                            "sweep_parameter": "beta",
                            "fixed_parameter_name": "mu",
                            "fixed_parameter_value": float(mu),
                            "extremum_parameter_value": float(BETA_VALUES[beta_pos]),
                            "frequency_hz": float(freq[beta_pos]),
                            "extremum_kind": "local_max",
                        }
                    )
                if left < -tol and right > tol:
                    rows.append(
                        {
                            "branch_id": branch_id,
                            "reference_mode_type_at_beta0": reference_type,
                            "sweep_parameter": "beta",
                            "fixed_parameter_name": "mu",
                            "fixed_parameter_value": float(mu),
                            "extremum_parameter_value": float(BETA_VALUES[beta_pos]),
                            "frequency_hz": float(freq[beta_pos]),
                            "extremum_kind": "local_min",
                        }
                    )

    return rows


def build_grid_rows(
    metadata: list[dict[str, float | int | str]],
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
    tracked_frac: np.ndarray,
    tracked_mac: np.ndarray,
    lower_gap: np.ndarray,
    upper_gap: np.ndarray,
    min_gap: np.ndarray,
    df_dmu: np.ndarray,
    df_dbeta: np.ndarray,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for beta_pos, beta_deg in enumerate(BETA_VALUES):
        for mu_pos, mu in enumerate(MU_VALUES):
            for branch_pos, branch_meta in enumerate(metadata):
                rows.append(
                    {
                        "beta": float(beta_deg),
                        "mu": float(mu),
                        "branch_id": str(branch_meta["branch_id"]),
                        "reference_mode_index_at_beta0": int(branch_meta["reference_mode_index_at_beta0"]),
                        "reference_mode_type_at_beta0": str(branch_meta["reference_mode_type_at_beta0"]),
                        "reference_descendant_label": str(branch_meta["reference_descendant_label"]),
                        "current_sorted_index": int(tracked_index[beta_pos, mu_pos, branch_pos]),
                        "frequency_hz": float(tracked_freqs[beta_pos, mu_pos, branch_pos]),
                        "axial_fraction": float(tracked_frac[beta_pos, mu_pos, branch_pos]),
                        "relative_gap_lower": float(lower_gap[beta_pos, mu_pos, branch_pos]),
                        "relative_gap_upper": float(upper_gap[beta_pos, mu_pos, branch_pos]),
                        "relative_gap_min_neighbor": float(min_gap[beta_pos, mu_pos, branch_pos]),
                        "df_dmu_hz_per_unit": float(df_dmu[beta_pos, mu_pos, branch_pos]),
                        "df_dbeta_hz_per_deg": float(df_dbeta[beta_pos, mu_pos, branch_pos]),
                        "assignment_mac": float(tracked_mac[beta_pos, mu_pos, branch_pos]),
                    }
                )
    return rows


def build_branch_summary_rows(
    metadata: list[dict[str, float | int | str]],
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
    tracked_frac: np.ndarray,
    min_gap: np.ndarray,
    df_dmu: np.ndarray,
    df_dbeta: np.ndarray,
    extrema_rows: list[dict[str, float | int | str]],
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    extrema_by_branch: dict[str, list[dict[str, float | int | str]]] = {}
    for row in extrema_rows:
        extrema_by_branch.setdefault(str(row["branch_id"]), []).append(row)

    for branch_pos, branch_meta in enumerate(metadata):
        branch_id = str(branch_meta["branch_id"])
        branch_extrema = extrema_by_branch.get(branch_id, [])
        mu_extrema_count = sum(1 for row in branch_extrema if row["sweep_parameter"] == "mu")
        beta_extrema_count = sum(1 for row in branch_extrema if row["sweep_parameter"] == "beta")
        mu_sign_changes = [
            branch_slice_sign_changes(df_dmu[beta_pos, :, branch_pos])
            for beta_pos in range(len(BETA_VALUES))
        ]
        beta_sign_changes = [
            branch_slice_sign_changes(df_dbeta[:, mu_pos, branch_pos])
            for mu_pos in range(len(MU_VALUES))
        ]
        rows.append(
            {
                "branch_id": branch_id,
                "reference_mode_index_at_beta0": int(branch_meta["reference_mode_index_at_beta0"]),
                "reference_mode_type_at_beta0": str(branch_meta["reference_mode_type_at_beta0"]),
                "reference_descendant_label": str(branch_meta["reference_descendant_label"]),
                "frequency_hz_min": float(np.nanmin(tracked_freqs[:, :, branch_pos])),
                "frequency_hz_max": float(np.nanmax(tracked_freqs[:, :, branch_pos])),
                "current_sorted_index_min": int(np.nanmin(tracked_index[:, :, branch_pos])),
                "current_sorted_index_max": int(np.nanmax(tracked_index[:, :, branch_pos])),
                "axial_fraction_min": float(np.nanmin(tracked_frac[:, :, branch_pos])),
                "axial_fraction_max": float(np.nanmax(tracked_frac[:, :, branch_pos])),
                "relative_gap_min_neighbor_min": float(np.nanmin(min_gap[:, :, branch_pos])),
                "df_dmu_min": float(np.nanmin(df_dmu[:, :, branch_pos])),
                "df_dmu_max": float(np.nanmax(df_dmu[:, :, branch_pos])),
                "df_dbeta_min": float(np.nanmin(df_dbeta[:, :, branch_pos])),
                "df_dbeta_max": float(np.nanmax(df_dbeta[:, :, branch_pos])),
                "mu_monotonicity_global": monotonicity_label(df_dmu[:, :, branch_pos]),
                "beta_monotonicity_global": monotonicity_label(df_dbeta[:, :, branch_pos]),
                "mu_sign_changes_max_over_beta": int(max(mu_sign_changes)),
                "beta_sign_changes_max_over_mu": int(max(beta_sign_changes)),
                "mu_extrema_count": int(mu_extrema_count),
                "beta_extrema_count": int(beta_extrema_count),
                "points_in_first_15": int(np.sum(tracked_index[:, :, branch_pos] <= LOW_SPECTRUM_CUTOFF)),
            }
        )
    return rows


def build_pair_gap_rows(
    metadata: list[dict[str, float | int | str]],
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
    tracked_frac: np.ndarray,
) -> list[dict[str, float | int | str]]:
    branch_ids = [str(item["branch_id"]) for item in metadata]
    ref_types = [str(item["reference_mode_type_at_beta0"]) for item in metadata]
    pair_best: dict[tuple[str, str], dict[str, float | int | str]] = {}

    for beta_pos, beta_deg in enumerate(BETA_VALUES):
        for mu_pos, mu in enumerate(MU_VALUES):
            order = np.argsort(tracked_index[beta_pos, mu_pos])
            for left_pos, right_pos in zip(order[:-1], order[1:]):
                left_idx = int(tracked_index[beta_pos, mu_pos, left_pos])
                right_idx = int(tracked_index[beta_pos, mu_pos, right_pos])
                if right_idx - left_idx != 1:
                    continue
                left_freq = float(tracked_freqs[beta_pos, mu_pos, left_pos])
                right_freq = float(tracked_freqs[beta_pos, mu_pos, right_pos])
                mean_freq = max(0.5 * (left_freq + right_freq), 1.0)
                gap = (right_freq - left_freq) / mean_freq
                pair_key = tuple(sorted((branch_ids[left_pos], branch_ids[right_pos])))
                candidate = {
                    "branch_id_left": branch_ids[left_pos],
                    "branch_id_right": branch_ids[right_pos],
                    "reference_mode_type_left_at_beta0": ref_types[left_pos],
                    "reference_mode_type_right_at_beta0": ref_types[right_pos],
                    "beta": float(beta_deg),
                    "mu": float(mu),
                    "current_sorted_index_left": left_idx,
                    "current_sorted_index_right": right_idx,
                    "frequency_left_hz": left_freq,
                    "frequency_right_hz": right_freq,
                    "axial_fraction_left": float(tracked_frac[beta_pos, mu_pos, left_pos]),
                    "axial_fraction_right": float(tracked_frac[beta_pos, mu_pos, right_pos]),
                    "relative_gap_pair": float(gap),
                }
                if pair_key not in pair_best or gap < float(pair_best[pair_key]["relative_gap_pair"]):
                    pair_best[pair_key] = candidate

    rows = list(pair_best.values())
    rows.sort(key=lambda row: float(row["relative_gap_pair"]))
    return rows


def build_point_gap_rows(
    metadata: list[dict[str, float | int | str]],
    min_gap: np.ndarray,
) -> list[dict[str, float | int | str]]:
    branch_ids = [str(item["branch_id"]) for item in metadata]
    ref_types = [str(item["reference_mode_type_at_beta0"]) for item in metadata]
    rows: list[dict[str, float | int | str]] = []

    point_branch = np.nanargmin(min_gap, axis=2)
    point_gap = np.take_along_axis(min_gap, point_branch[:, :, None], axis=2)[:, :, 0]

    for beta_pos, beta_deg in enumerate(BETA_VALUES):
        for mu_pos, mu in enumerate(MU_VALUES):
            branch_pos = int(point_branch[beta_pos, mu_pos])
            rows.append(
                {
                    "beta": float(beta_deg),
                    "mu": float(mu),
                    "min_relative_gap": float(point_gap[beta_pos, mu_pos]),
                    "responsible_branch_id": branch_ids[branch_pos],
                    "responsible_reference_mode_type_at_beta0": ref_types[branch_pos],
                }
            )
    return rows


def branches_for_low_spectrum_plot(
    metadata: list[dict[str, float | int | str]],
    tracked_index: np.ndarray,
) -> list[str]:
    branch_ids = [str(item["branch_id"]) for item in metadata]
    selected = []
    for branch_pos, branch_id in enumerate(branch_ids):
        if int(np.nanmin(tracked_index[:, :, branch_pos])) <= LOW_SPECTRUM_CUTOFF:
            selected.append(branch_id)
    return selected


def find_beta_position(beta_deg: float) -> int:
    return int(np.argmin(np.abs(BETA_VALUES - beta_deg)))


def find_mu_position(mu: float) -> int:
    return int(np.argmin(np.abs(MU_VALUES - mu)))


def plot_freq_vs_mu(
    metadata: list[dict[str, float | int | str]],
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
) -> None:
    selected_branch_ids = branches_for_low_spectrum_plot(metadata, tracked_index)
    color_map = branch_palette(selected_branch_ids)
    branch_ids = [str(item["branch_id"]) for item in metadata]
    plot_indices = [branch_ids.index(branch_id) for branch_id in selected_branch_ids]

    ncols = 2
    nrows = ceil(len(SELECTED_BETAS) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 4.8 * nrows), squeeze=False)

    for subplot_pos, beta_deg in enumerate(SELECTED_BETAS):
        beta_pos = find_beta_position(beta_deg)
        ax = axes[subplot_pos // ncols][subplot_pos % ncols]
        for branch_pos in plot_indices:
            branch_id = branch_ids[branch_pos]
            ax.plot(
                MU_VALUES,
                tracked_freqs[beta_pos, :, branch_pos],
                linewidth=1.8,
                color=color_map[branch_id],
                label=branch_id,
            )
        ax.set_title(f"beta = {BETA_VALUES[beta_pos]:.1f} deg")
        ax.set_xlabel("mu")
        ax.set_ylabel("Frequency (Hz)")
        ax.grid(True, alpha=0.25)

    for subplot_pos in range(len(SELECTED_BETAS), nrows * ncols):
        axes[subplot_pos // ncols][subplot_pos % ncols].axis("off")

    handles = [plt.Line2D([0], [0], color=color_map[branch_id], lw=2.0, label=branch_id) for branch_id in selected_branch_ids]
    fig.suptitle("Tracked branch frequencies vs mu for fixed beta", y=0.995)
    fig.legend(handles=handles, loc="upper center", ncols=4, bbox_to_anchor=(0.5, 0.99))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FREQ_MU_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_freq_vs_beta(
    metadata: list[dict[str, float | int | str]],
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
) -> None:
    selected_branch_ids = branches_for_low_spectrum_plot(metadata, tracked_index)
    color_map = branch_palette(selected_branch_ids)
    branch_ids = [str(item["branch_id"]) for item in metadata]
    plot_indices = [branch_ids.index(branch_id) for branch_id in selected_branch_ids]

    ncols = 2
    nrows = ceil(len(SELECTED_MUS) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 4.8 * nrows), squeeze=False)

    for subplot_pos, mu in enumerate(SELECTED_MUS):
        mu_pos = find_mu_position(mu)
        ax = axes[subplot_pos // ncols][subplot_pos % ncols]
        for branch_pos in plot_indices:
            branch_id = branch_ids[branch_pos]
            ax.plot(
                BETA_VALUES,
                tracked_freqs[:, mu_pos, branch_pos],
                linewidth=1.8,
                color=color_map[branch_id],
                label=branch_id,
            )
        ax.set_title(f"mu = {MU_VALUES[mu_pos]:.3f}")
        ax.set_xlabel("beta (deg)")
        ax.set_ylabel("Frequency (Hz)")
        ax.grid(True, alpha=0.25)

    for subplot_pos in range(len(SELECTED_MUS), nrows * ncols):
        axes[subplot_pos // ncols][subplot_pos % ncols].axis("off")

    handles = [plt.Line2D([0], [0], color=color_map[branch_id], lw=2.0, label=branch_id) for branch_id in selected_branch_ids]
    fig.suptitle("Tracked branch frequencies vs beta for fixed mu", y=0.995)
    fig.legend(handles=handles, loc="upper center", ncols=4, bbox_to_anchor=(0.5, 0.99))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(FREQ_BETA_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_sorted_spectrum_vs_mu(
    metadata: list[dict[str, float | int | str]],
    full_freqs: np.ndarray,
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
) -> None:
    selected_branch_ids = branches_for_low_spectrum_plot(metadata, tracked_index)
    color_map = branch_palette(selected_branch_ids)
    branch_ids = [str(item["branch_id"]) for item in metadata]
    plot_indices = [branch_ids.index(branch_id) for branch_id in selected_branch_ids]

    ncols = 2
    nrows = ceil(len(SELECTED_BETAS) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 4.8 * nrows), squeeze=False)

    for subplot_pos, beta_deg in enumerate(SELECTED_BETAS):
        beta_pos = find_beta_position(beta_deg)
        ax = axes[subplot_pos // ncols][subplot_pos % ncols]
        for mode_pos in range(min(24, full_freqs.shape[2])):
            ax.plot(MU_VALUES, full_freqs[beta_pos, :, mode_pos], color="0.85", linewidth=0.9, zorder=1)
        for branch_pos in plot_indices:
            branch_id = branch_ids[branch_pos]
            ax.scatter(
                MU_VALUES,
                tracked_freqs[beta_pos, :, branch_pos],
                s=12,
                color=color_map[branch_id],
                alpha=0.95,
                zorder=2,
                label=branch_id,
            )
        ax.set_title(f"beta = {BETA_VALUES[beta_pos]:.1f} deg")
        ax.set_xlabel("mu")
        ax.set_ylabel("Frequency (Hz)")
        ax.grid(True, alpha=0.25)

    for subplot_pos in range(len(SELECTED_BETAS), nrows * ncols):
        axes[subplot_pos // ncols][subplot_pos % ncols].axis("off")

    handles = [plt.Line2D([0], [0], color=color_map[branch_id], marker="o", linestyle="None", markersize=5, label=branch_id) for branch_id in selected_branch_ids]
    fig.suptitle("Sorted spectrum vs mu, colored by branch identity", y=0.995)
    fig.legend(handles=handles, loc="upper center", ncols=4, bbox_to_anchor=(0.5, 0.99))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(SORTED_SPECTRUM_MU_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_sorted_spectrum_vs_beta(
    metadata: list[dict[str, float | int | str]],
    full_freqs: np.ndarray,
    tracked_freqs: np.ndarray,
    tracked_index: np.ndarray,
) -> None:
    selected_branch_ids = branches_for_low_spectrum_plot(metadata, tracked_index)
    color_map = branch_palette(selected_branch_ids)
    branch_ids = [str(item["branch_id"]) for item in metadata]
    plot_indices = [branch_ids.index(branch_id) for branch_id in selected_branch_ids]

    ncols = 2
    nrows = ceil(len(SELECTED_MUS) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(14.0, 4.8 * nrows), squeeze=False)

    for subplot_pos, mu in enumerate(SELECTED_MUS):
        mu_pos = find_mu_position(mu)
        ax = axes[subplot_pos // ncols][subplot_pos % ncols]
        for mode_pos in range(min(24, full_freqs.shape[2])):
            ax.plot(BETA_VALUES, full_freqs[:, mu_pos, mode_pos], color="0.85", linewidth=0.9, zorder=1)
        for branch_pos in plot_indices:
            branch_id = branch_ids[branch_pos]
            ax.scatter(
                BETA_VALUES,
                tracked_freqs[:, mu_pos, branch_pos],
                s=12,
                color=color_map[branch_id],
                alpha=0.95,
                zorder=2,
                label=branch_id,
            )
        ax.set_title(f"mu = {MU_VALUES[mu_pos]:.3f}")
        ax.set_xlabel("beta (deg)")
        ax.set_ylabel("Frequency (Hz)")
        ax.grid(True, alpha=0.25)

    for subplot_pos in range(len(SELECTED_MUS), nrows * ncols):
        axes[subplot_pos // ncols][subplot_pos % ncols].axis("off")

    handles = [plt.Line2D([0], [0], color=color_map[branch_id], marker="o", linestyle="None", markersize=5, label=branch_id) for branch_id in selected_branch_ids]
    fig.suptitle("Sorted spectrum vs beta, colored by branch identity", y=0.995)
    fig.legend(handles=handles, loc="upper center", ncols=4, bbox_to_anchor=(0.5, 0.99))
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(SORTED_SPECTRUM_BETA_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_axial_fraction_heatmaps(
    metadata: list[dict[str, float | int | str]],
    tracked_frac: np.ndarray,
) -> None:
    branch_lookup = {str(item["branch_id"]): idx for idx, item in enumerate(metadata)}
    selected = [branch_id for branch_id in SELECTED_AXIAL_IDS if branch_id in branch_lookup]

    ncols = 1
    nrows = len(selected)
    fig, axes = plt.subplots(nrows, ncols, figsize=(10.0, 3.8 * nrows), squeeze=False)
    last_image = None

    for row_pos, branch_id in enumerate(selected):
        branch_pos = branch_lookup[branch_id]
        ax = axes[row_pos][0]
        last_image = ax.imshow(
            tracked_frac[:, :, branch_pos],
            origin="lower",
            aspect="auto",
            extent=(MU_VALUES[0], MU_VALUES[-1], BETA_VALUES[0], BETA_VALUES[-1]),
            vmin=0.0,
            vmax=1.0,
            cmap="viridis",
        )
        ax.contour(
            MU_VALUES,
            BETA_VALUES,
            tracked_frac[:, :, branch_pos],
            levels=[0.1, 0.5, 0.9],
            colors=["white", "orange", "red"],
            linewidths=0.9,
        )
        ax.set_title(f"{branch_id}: axial_fraction(beta, mu)")
        ax.set_xlabel("mu")
        ax.set_ylabel("beta (deg)")

    if last_image is not None:
        fig.colorbar(last_image, ax=axes.ravel().tolist(), label="axial_fraction")

    fig.tight_layout()
    fig.savefig(AXIAL_HEATMAP_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_min_gap_map(point_gap_rows: list[dict[str, float | int | str]]) -> None:
    gap_grid = np.full((len(BETA_VALUES), len(MU_VALUES)), np.nan, dtype=float)
    for row in point_gap_rows:
        beta_pos = find_beta_position(float(row["beta"]))
        mu_pos = find_mu_position(float(row["mu"]))
        gap_grid[beta_pos, mu_pos] = float(row["min_relative_gap"])

    fig, ax = plt.subplots(figsize=(10.0, 5.2))
    image = ax.imshow(
        gap_grid,
        origin="lower",
        aspect="auto",
        extent=(MU_VALUES[0], MU_VALUES[-1], BETA_VALUES[0], BETA_VALUES[-1]),
        cmap="magma_r",
        norm=LogNorm(vmin=max(float(np.nanmin(gap_grid)), 1e-6), vmax=float(np.nanmax(gap_grid))),
    )
    ax.set_xlabel("mu")
    ax.set_ylabel("beta (deg)")
    ax.set_title("Minimal neighbor gap over tracked branches")
    fig.colorbar(image, ax=ax, label="min relative gap")
    fig.tight_layout()
    fig.savefig(MIN_GAP_PLOT_PATH, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    tracked = track_branch_grid()
    metadata = tracked["metadata"]
    full_freqs = np.asarray(tracked["full_freqs"], dtype=float)
    tracked_freqs = np.asarray(tracked["tracked_freqs"], dtype=float)
    tracked_index = np.asarray(tracked["tracked_index"], dtype=int)
    tracked_frac = np.asarray(tracked["tracked_frac"], dtype=float)
    tracked_mac = np.asarray(tracked["tracked_mac"], dtype=float)

    lower_gap, upper_gap, min_gap = compute_neighbor_gaps(
        full_freqs=full_freqs,
        tracked_freqs=tracked_freqs,
        tracked_index=tracked_index,
    )
    df_dbeta = np.gradient(tracked_freqs, BETA_VALUES, axis=0, edge_order=2)
    df_dmu = np.gradient(tracked_freqs, MU_VALUES, axis=1, edge_order=2)

    extrema_rows = collect_extrema_rows(
        metadata=metadata,
        tracked_freqs=tracked_freqs,
        df_dmu=df_dmu,
        df_dbeta=df_dbeta,
    )
    grid_rows = build_grid_rows(
        metadata=metadata,
        tracked_freqs=tracked_freqs,
        tracked_index=tracked_index,
        tracked_frac=tracked_frac,
        tracked_mac=tracked_mac,
        lower_gap=lower_gap,
        upper_gap=upper_gap,
        min_gap=min_gap,
        df_dmu=df_dmu,
        df_dbeta=df_dbeta,
    )
    branch_summary_rows = build_branch_summary_rows(
        metadata=metadata,
        tracked_freqs=tracked_freqs,
        tracked_index=tracked_index,
        tracked_frac=tracked_frac,
        min_gap=min_gap,
        df_dmu=df_dmu,
        df_dbeta=df_dbeta,
        extrema_rows=extrema_rows,
    )
    pair_gap_rows = build_pair_gap_rows(
        metadata=metadata,
        tracked_freqs=tracked_freqs,
        tracked_index=tracked_index,
        tracked_frac=tracked_frac,
    )
    point_gap_rows = build_point_gap_rows(metadata=metadata, min_gap=min_gap)

    write_csv_rows(BRANCH_GRID_PATH, grid_rows)
    write_csv_rows(BRANCH_SUMMARY_PATH, branch_summary_rows)
    write_csv_rows(EXTREMA_PATH, extrema_rows)
    write_csv_rows(PAIR_GAPS_PATH, pair_gap_rows)
    write_csv_rows(POINT_GAPS_PATH, point_gap_rows)

    plot_freq_vs_mu(metadata=metadata, tracked_freqs=tracked_freqs, tracked_index=tracked_index)
    plot_freq_vs_beta(metadata=metadata, tracked_freqs=tracked_freqs, tracked_index=tracked_index)
    plot_sorted_spectrum_vs_mu(
        metadata=metadata,
        full_freqs=full_freqs,
        tracked_freqs=tracked_freqs,
        tracked_index=tracked_index,
    )
    plot_sorted_spectrum_vs_beta(
        metadata=metadata,
        full_freqs=full_freqs,
        tracked_freqs=tracked_freqs,
        tracked_index=tracked_index,
    )
    plot_axial_fraction_heatmaps(metadata=metadata, tracked_frac=tracked_frac)
    plot_min_gap_map(point_gap_rows=point_gap_rows)

    print(f"tracked branches: {len(metadata)}")
    print(f"beta grid: {BETA_VALUES[0]:.1f} .. {BETA_VALUES[-1]:.1f} ({len(BETA_VALUES)} points)")
    print(f"mu grid: {MU_VALUES[0]:.3f} .. {MU_VALUES[-1]:.3f} ({len(MU_VALUES)} points)")
    print(f"saved table: {BRANCH_GRID_PATH}")
    print(f"saved table: {BRANCH_SUMMARY_PATH}")
    print(f"saved table: {EXTREMA_PATH}")
    print(f"saved table: {PAIR_GAPS_PATH}")
    print(f"saved table: {POINT_GAPS_PATH}")
    print(f"saved plot: {FREQ_MU_PLOT_PATH}")
    print(f"saved plot: {FREQ_BETA_PLOT_PATH}")
    print(f"saved plot: {SORTED_SPECTRUM_MU_PLOT_PATH}")
    print(f"saved plot: {SORTED_SPECTRUM_BETA_PLOT_PATH}")
    print(f"saved plot: {AXIAL_HEATMAP_PATH}")
    print(f"saved plot: {MIN_GAP_PLOT_PATH}")


if __name__ == "__main__":
    main()
