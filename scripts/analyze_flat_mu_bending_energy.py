from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
from scripts.sweep_grid_policy import (  # noqa: E402
    ANALYSIS_BETA_STEP,
    ANALYSIS_MU_STEP,
    analysis_beta_grid,
    analysis_mu_grid,
)
from my_project.fem import python_fem as fem  # noqa: E402


RESULTS_DIR = REPO_ROOT / "results"
CANDIDATE_PATH = RESULTS_DIR / "flat_mu_branch_candidates.csv"
REPRESENTATIVE_PATH = RESULTS_DIR / "flat_mu_representative_cases.csv"
ENERGY_CSV_PATH = RESULTS_DIR / "flat_mu_bending_energy_cases.csv"
SUMMARY_PLOT_PATH = RESULTS_DIR / "flat_mu_bending_energy_summary.png"

PRIMARY_BETA = 15.0
SECONDARY_BETA = 7.5
TARGET_RADII = (0.005, 0.02)
TARGET_MUS = (0.0, 0.5, 0.9)
TARGET_BRANCH_COUNT = 3
MIN_BRANCH_COUNT = 2
N_TRACK = 20
N_SOLVE = 36
FREQ_WEIGHT = 0.03
AXIAL_THRESHOLD = 0.9
BENDING_THRESHOLD = 0.1
MODE_SHAPE_SCALE = 0.22

MU_TRACK_VALUES = analysis_mu_grid(anchor_points=TARGET_MUS)


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


def mu_key(mu: float) -> float:
    return round(float(mu), 10)


def branch_sort_key(branch_id: str) -> tuple[int, str]:
    try:
        return int(branch_id.split("_")[-1]), branch_id
    except Exception:
        return (10**9, branch_id)


def reference_mode_type(axial_fraction: float) -> str:
    if axial_fraction >= AXIAL_THRESHOLD:
        return "axial"
    if axial_fraction <= BENDING_THRESHOLD:
        return "bending"
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
        raise RuntimeError("Failed to assign all tracked branches.")
    return assignment, assigned_mac


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


def build_node_coordinates(mu: float, beta_deg: float, ell: float) -> tuple[np.ndarray, np.ndarray]:
    n1 = n2 = fem.N_ELEM
    l1 = ell * (1.0 - mu)
    l2 = ell * (1.0 + mu)
    beta_rad = np.deg2rad(beta_deg)

    x = np.zeros(n1 + n2 + 1, dtype=float)
    y = np.zeros(n1 + n2 + 1, dtype=float)
    x[: n1 + 1] = np.linspace(0.0, l1, n1 + 1)
    s2 = np.linspace(0.0, l2, n2 + 1)
    x[n1:] = l1 + s2 * np.cos(beta_rad)
    y[n1:] = s2 * np.sin(beta_rad)
    return x, y


def element_bending_stiffness(le: float) -> np.ndarray:
    c = 1.0 / le**3
    return c * np.array(
        [
            [12.0, 6.0 * le, -12.0, 6.0 * le],
            [6.0 * le, 4.0 * le**2, -6.0 * le, 2.0 * le**2],
            [-12.0, -6.0 * le, 12.0, -6.0 * le],
            [6.0 * le, 2.0 * le**2, -6.0 * le, 4.0 * le**2],
        ],
        dtype=float,
    )


def bending_energy_fractions(full_vec: np.ndarray, mu: float, beta_deg: float, ell: float) -> tuple[float, float]:
    n1 = n2 = fem.N_ELEM
    l1 = ell * (1.0 - mu)
    l2 = ell * (1.0 + mu)
    le1 = l1 / n1
    le2 = l2 / n2
    kb1 = element_bending_stiffness(le1)
    kb2 = element_bending_stiffness(le2)
    t6 = fem.rotation_matrix_6x6(np.deg2rad(beta_deg))

    left_energy = 0.0
    for i in range(n1):
        dofs = [3 * i, 3 * i + 1, 3 * i + 2, 3 * (i + 1), 3 * (i + 1) + 1, 3 * (i + 1) + 2]
        q_global = full_vec[dofs]
        q_b = q_global[[1, 2, 4, 5]]
        left_energy += 0.5 * float(q_b @ kb1 @ q_b)

    right_energy = 0.0
    for j in range(n2):
        dofs = [
            3 * (n1 + j),
            3 * (n1 + j) + 1,
            3 * (n1 + j) + 2,
            3 * (n1 + j + 1),
            3 * (n1 + j + 1) + 1,
            3 * (n1 + j + 1) + 2,
        ]
        q_global = full_vec[dofs]
        q_local = t6 @ q_global
        q_b = q_local[[1, 2, 4, 5]]
        right_energy += 0.5 * float(q_b @ kb2 @ q_b)

    left_energy = max(left_energy, 0.0)
    right_energy = max(right_energy, 0.0)
    total = left_energy + right_energy
    if total <= 0.0:
        return np.nan, np.nan
    return left_energy / total, right_energy / total


def choose_selected_branches() -> list[dict[str, float | str | int]]:
    candidate_rows = load_csv_rows(CANDIDATE_PATH)
    representative_rows = load_csv_rows(REPRESENTATIVE_PATH)

    def collect_for_beta(beta_deg: float) -> list[dict[str, float | str | int]]:
        selected_slots = sorted(
            (
                row
                for row in candidate_rows
                if abs(float(row["beta"]) - beta_deg) < 1e-9 and int(row["selected_as_flattening_candidate"]) == 1
            ),
            key=lambda row: float(row["spread_ratio_rmax_over_rmin"]),
        )
        picked: list[dict[str, float | str | int]] = []
        used = set()
        for slot_row in selected_slots:
            low_slot = int(slot_row["low_slot"])
            rep_row = next(
                (
                    row
                    for row in representative_rows
                    if abs(float(row["beta"]) - beta_deg) < 1e-9
                    and abs(float(row["radius_mm"]) - 20.0) < 1e-9
                    and abs(float(row["mu"]) - 0.0) < 1e-9
                    and int(row["low_slot"]) == low_slot
                ),
                None,
            )
            if rep_row is None:
                continue
            branch_id = str(rep_row["dominant_occupant_branch_id_over_mu"])
            if branch_id in used:
                continue
            used.add(branch_id)
            picked.append(
                {
                    "beta": float(beta_deg),
                    "branch_id": branch_id,
                    "source_low_slot": low_slot,
                    "spread_ratio_rmax_over_rmin": float(slot_row["spread_ratio_rmax_over_rmin"]),
                }
            )
            if len(picked) >= TARGET_BRANCH_COUNT:
                break
        return picked

    selected = collect_for_beta(PRIMARY_BETA)
    if len(selected) < MIN_BRANCH_COUNT:
        supplemental = collect_for_beta(SECONDARY_BETA)
        for item in supplemental:
            if item["branch_id"] not in {row["branch_id"] for row in selected}:
                selected.append(item)
            if len(selected) >= TARGET_BRANCH_COUNT:
                break
    return selected


def track_branch_family(radius: float, beta_deg: float) -> dict[str, object]:
    seed_freqs, seed_vecs = solve_fem_modes(radius=radius, mu=0.0, beta_deg=0.0, n_modes=N_SOLVE)
    seed_fracs = fem_axial_fractions(seed_vecs)
    metadata = build_branch_metadata(seed_freqs=seed_freqs[:N_TRACK], seed_fracs=seed_fracs[:N_TRACK])

    prev_vecs = seed_vecs[:, :N_TRACK].copy()
    prev_freqs = seed_freqs[:N_TRACK].copy()
    prev_indices = np.arange(1, N_TRACK + 1, dtype=int)
    prev_fracs = seed_fracs[:N_TRACK].copy()

    beta_values = analysis_beta_grid(start=0.0, stop=beta_deg, anchor_points=(beta_deg,))
    for beta_value in beta_values[1:]:
        cur_freqs, cur_vecs = solve_fem_modes(radius=radius, mu=0.0, beta_deg=float(beta_value), n_modes=N_SOLVE)
        cur_fracs = fem_axial_fractions(cur_vecs)
        assign, _ = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        prev_vecs = cur_vecs[:, assign]
        prev_freqs = cur_freqs[assign]
        prev_indices = assign + 1
        prev_fracs = cur_fracs[assign]

    snapshots: dict[float, dict[str, np.ndarray]] = {
        mu_key(0.0): {
            "frequency_hz": prev_freqs.copy(),
            "current_sorted_index": prev_indices.copy(),
            "axial_fraction": prev_fracs.copy(),
            "vecs": prev_vecs.copy(),
        }
    }

    for mu_value in MU_TRACK_VALUES[1:]:
        cur_freqs, cur_vecs = solve_fem_modes(radius=radius, mu=float(mu_value), beta_deg=beta_deg, n_modes=N_SOLVE)
        cur_fracs = fem_axial_fractions(cur_vecs)
        assign, _ = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        prev_vecs = cur_vecs[:, assign]
        prev_freqs = cur_freqs[assign]
        prev_indices = assign + 1
        prev_fracs = cur_fracs[assign]
        snapshots[mu_key(mu_value)] = {
            "frequency_hz": prev_freqs.copy(),
            "current_sorted_index": prev_indices.copy(),
            "axial_fraction": prev_fracs.copy(),
            "vecs": prev_vecs.copy(),
        }

    return {
        "radius": float(radius),
        "beta": float(beta_deg),
        "metadata": metadata,
        "snapshots": snapshots,
    }


def collect_case_rows(
    selected_branches: list[dict[str, float | str | int]],
) -> tuple[list[dict[str, float | int | str]], dict[tuple[float, float], dict[str, object]]]:
    rows: list[dict[str, float | int | str]] = []
    cache: dict[tuple[float, float], dict[str, object]] = {}

    for item in selected_branches:
        beta_deg = float(item["beta"])
        branch_id = str(item["branch_id"])

        for radius in TARGET_RADII:
            key = (float(radius), beta_deg)
            if key not in cache:
                cache[key] = track_branch_family(radius=float(radius), beta_deg=beta_deg)

            tracking = cache[key]
            metadata = list(tracking["metadata"])
            branch_lookup = {str(meta["branch_id"]): idx for idx, meta in enumerate(metadata)}
            if branch_id not in branch_lookup:
                raise KeyError(f"Selected branch {branch_id} is not present in tracked metadata.")
            branch_pos = branch_lookup[branch_id]
            ell = build_params(radius).L_total / 2.0

            for mu_value in TARGET_MUS:
                snapshot = tracking["snapshots"][mu_key(mu_value)]
                phi = np.asarray(snapshot["vecs"][:, branch_pos], dtype=float)
                full_vec = expand_full_vector(phi)
                left_frac, right_frac = bending_energy_fractions(
                    full_vec=full_vec,
                    mu=float(mu_value),
                    beta_deg=beta_deg,
                    ell=ell,
                )
                rows.append(
                    {
                        "beta": beta_deg,
                        "branch_id": branch_id,
                        "r": float(radius),
                        "radius_mm": float(radius * 1e3),
                        "mu": float(mu_value),
                        "current_sorted_index": int(snapshot["current_sorted_index"][branch_pos]),
                        "frequency_hz": float(snapshot["frequency_hz"][branch_pos]),
                        "left_bending_energy_fraction": float(left_frac),
                        "right_bending_energy_fraction": float(right_frac),
                        "axial_fraction": float(snapshot["axial_fraction"][branch_pos]),
                    }
                )

    rows.sort(key=lambda row: (float(row["beta"]), branch_sort_key(str(row["branch_id"])), float(row["r"]), float(row["mu"])))
    return rows, cache


def build_shape_plot_for_branch(
    branch_id: str,
    beta_deg: float,
    case_rows: list[dict[str, float | int | str]],
    tracking_cache: dict[tuple[float, float], dict[str, object]],
) -> Path:
    rows = [row for row in case_rows if str(row["branch_id"]) == branch_id and abs(float(row["beta"]) - beta_deg) < 1e-9]
    rows.sort(key=lambda row: (float(row["r"]), float(row["mu"])))
    if not rows:
        raise ValueError(f"No rows for branch {branch_id} at beta={beta_deg}")

    fig, axes = plt.subplots(len(TARGET_RADII), len(TARGET_MUS), figsize=(15.5, 8.1), squeeze=False)
    normalized_cases = []
    ell = build_params(TARGET_RADII[0]).L_total / 2.0

    for row in rows:
        radius = float(row["r"])
        mu_value = float(row["mu"])
        tracking = tracking_cache[(radius, beta_deg)]
        metadata = list(tracking["metadata"])
        branch_lookup = {str(meta["branch_id"]): idx for idx, meta in enumerate(metadata)}
        branch_pos = branch_lookup[branch_id]
        snapshot = tracking["snapshots"][mu_key(mu_value)]
        full_vec = expand_full_vector(np.asarray(snapshot["vecs"][:, branch_pos], dtype=float))
        x_base, y_base = build_node_coordinates(mu=mu_value, beta_deg=beta_deg, ell=ell)
        u = full_vec[0::3]
        v = full_vec[1::3]
        max_disp = float(np.max(np.sqrt(u * u + v * v)))
        norm = max(max_disp, 1e-12)
        normalized_cases.append(
            {
                "row": row,
                "x_base": x_base,
                "y_base": y_base,
                "u": u / norm,
                "v": v / norm,
            }
        )

    x_all = []
    y_all = []
    for item in normalized_cases:
        x_def = item["x_base"] + MODE_SHAPE_SCALE * item["u"]
        y_def = item["y_base"] + MODE_SHAPE_SCALE * item["v"]
        x_all.extend([float(np.min(item["x_base"])), float(np.max(item["x_base"])), float(np.min(x_def)), float(np.max(x_def))])
        y_all.extend([float(np.min(item["y_base"])), float(np.max(item["y_base"])), float(np.min(y_def)), float(np.max(y_def))])

    x_margin = 0.06 * max(max(x_all) - min(x_all), 1.0)
    y_margin = 0.10 * max(max(y_all) - min(y_all), 0.5)
    x_limits = (min(x_all) - x_margin, max(x_all) + x_margin)
    y_limits = (min(y_all) - y_margin, max(y_all) + y_margin)

    for item in normalized_cases:
        row = item["row"]
        radius = float(row["r"])
        mu_value = float(row["mu"])
        row_pos = TARGET_RADII.index(radius)
        col_pos = TARGET_MUS.index(mu_value)
        ax = axes[row_pos, col_pos]
        joint = fem.N_ELEM
        x_def = item["x_base"] + MODE_SHAPE_SCALE * item["u"]
        y_def = item["y_base"] + MODE_SHAPE_SCALE * item["v"]

        ax.plot(item["x_base"][: joint + 1], item["y_base"][: joint + 1], color="0.75", linestyle="--", linewidth=1.0)
        ax.plot(item["x_base"][joint:], item["y_base"][joint:], color="0.75", linestyle="--", linewidth=1.0)
        ax.plot(x_def[: joint + 1], y_def[: joint + 1], color="#1f77b4", linewidth=2.0)
        ax.plot(x_def[joint:], y_def[joint:], color="#ff7f0e", linewidth=2.0)
        ax.scatter([item["x_base"][joint]], [item["y_base"][joint]], color="black", s=12, zorder=5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.grid(True, alpha=0.25)
        ax.set_title(
            f"beta={beta_deg:g} deg, r={radius * 1e3:.1f} mm, mu={mu_value:.1f}\n"
            f"{branch_id}, sorted={int(row['current_sorted_index'])}"
        )
        ax.text(
            0.02,
            0.02,
            (
                f"f={float(row['frequency_hz']):.2f} Hz\n"
                f"L bend={float(row['left_bending_energy_fraction']):.3f}\n"
                f"R bend={float(row['right_bending_energy_fraction']):.3f}\n"
                f"axial={float(row['axial_fraction']):.4f}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.5,
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.92},
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    handles = [
        Line2D([0], [0], color="0.75", linestyle="--", linewidth=1.0, label="undeformed geometry"),
        Line2D([0], [0], color="#1f77b4", linewidth=2.0, label="left arm"),
        Line2D([0], [0], color="#ff7f0e", linewidth=2.0, label="right arm"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965), fontsize=9)
    fig.suptitle(f"Bending-energy follow-up mode shapes: beta={beta_deg:g} deg, {branch_id}", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_path = RESULTS_DIR / f"flat_mu_bending_shapes_beta{str(beta_deg).replace('.', 'p')}_{branch_id}.png"
    fig.savefig(out_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out_path


def build_summary_plot(case_rows: list[dict[str, float | int | str]], selected_branches: list[dict[str, float | str | int]]) -> Path:
    branch_keys = [(float(item["beta"]), str(item["branch_id"])) for item in selected_branches]
    fig, axes = plt.subplots(1, len(branch_keys), figsize=(5.6 * len(branch_keys), 4.8), squeeze=False)
    axes_row = axes[0]

    colors = {0.005: "#1f77b4", 0.02: "#d62728"}
    mu_values = np.asarray(TARGET_MUS, dtype=float)

    for ax, (beta_deg, branch_id) in zip(axes_row, branch_keys):
        rows = [
            row
            for row in case_rows
            if abs(float(row["beta"]) - beta_deg) < 1e-9 and str(row["branch_id"]) == branch_id
        ]
        grouped: dict[float, list[dict[str, float | int | str]]] = defaultdict(list)
        for row in rows:
            grouped[float(row["r"])].append(row)

        for radius in TARGET_RADII:
            r_rows = sorted(grouped[float(radius)], key=lambda row: float(row["mu"]))
            left = np.asarray([float(row["left_bending_energy_fraction"]) for row in r_rows], dtype=float)
            right = np.asarray([float(row["right_bending_energy_fraction"]) for row in r_rows], dtype=float)
            ax.plot(mu_values, left, color=colors[float(radius)], linewidth=2.0, marker="o", label=f"r={radius*1e3:.0f} mm left")
            ax.plot(
                mu_values,
                right,
                color=colors[float(radius)],
                linewidth=2.0,
                linestyle="--",
                marker="s",
                label=f"r={radius*1e3:.0f} mm right",
            )

        ax.set_title(f"beta={beta_deg:g} deg, {branch_id}")
        ax.set_xlabel("mu")
        ax.set_ylabel("bending energy fraction")
        ax.set_ylim(0.0, 1.0)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle("Left/right bending-energy redistribution for the selected flattening branches", y=0.99)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(SUMMARY_PLOT_PATH, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return SUMMARY_PLOT_PATH


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    selected_branches = choose_selected_branches()
    if len(selected_branches) < MIN_BRANCH_COUNT:
        raise RuntimeError("Could not select enough representative branches from the existing flattening analysis.")

    case_rows, tracking_cache = collect_case_rows(selected_branches=selected_branches)
    write_csv_rows(ENERGY_CSV_PATH, case_rows)
    shape_paths = []
    for item in selected_branches:
        shape_paths.append(
            build_shape_plot_for_branch(
                branch_id=str(item["branch_id"]),
                beta_deg=float(item["beta"]),
                case_rows=case_rows,
                tracking_cache=tracking_cache,
            )
        )
    build_summary_plot(case_rows=case_rows, selected_branches=selected_branches)

    print("selected branches:")
    for item in selected_branches:
        print(
            f"  beta={float(item['beta']):g} deg, branch_id={item['branch_id']}, "
            f"source low slot={int(item['source_low_slot'])}, "
            f"spread_ratio_rmax_over_rmin={float(item['spread_ratio_rmax_over_rmin']):.3f}"
        )
    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    print(f"saved CSV: {ENERGY_CSV_PATH}")
    for path in shape_paths:
        print(f"saved PNG: {path}")
    print(f"saved PNG: {SUMMARY_PLOT_PATH}")


if __name__ == "__main__":
    main()
