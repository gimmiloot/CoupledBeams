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


REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analyze_flat_mu_bending_energy import (  # noqa: E402
    N_SOLVE,
    N_TRACK,
    MODE_SHAPE_SCALE,
    RESULTS_DIR,
    assign_by_mac_and_frequency,
    build_branch_metadata,
    build_params,
    expand_full_vector,
    solve_fem_modes,
)
from scripts.compare_beta0_analytic_vs_fem import fem_axial_fractions  # noqa: E402
from scripts.sweep_grid_policy import (  # noqa: E402
    ANALYSIS_BETA_STEP,
    ANALYSIS_MU_STEP,
    LOCAL_MU_REFINEMENT_STEP,
    analysis_beta_grid,
    analysis_mu_grid,
)
from my_project.analytic.FreqMuNet import roots_clamped_supported  # noqa: E402
from my_project.analytic.formulas import frequency_scale  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402


BETA_DEG = 15.0
RADIUS = 0.005
BRANCH_LABELS = {
    "bending_desc_01": "потомок 1-й изгибной ветви",
    "bending_desc_02": "потомок 2-й изгибной ветви",
    "bending_desc_04": "потомок 4-й изгибной ветви",
}
TARGET_BRANCH_IDS = tuple(BRANCH_LABELS)

FULL_CSV_PATH = RESULTS_DIR / "target_descendants_beta15_r5_analysis.csv"
SELECTED_CSV_PATH = RESULTS_DIR / "target_descendants_beta15_r5_selected_points.csv"
SHAPE_PATH_TEMPLATE = "target_descendant_beta15_r5_{branch_id}_selected_shapes.png"

REFERENCE_ROOTS = roots_clamped_supported(12)


def write_csv_rows(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def branch_sort_key(branch_id: str) -> tuple[int, str]:
    try:
        return int(branch_id.split("_")[-1]), branch_id
    except Exception:
        return (10**9, branch_id)


def lambda_from_frequency_hz(freq_hz: float) -> float:
    params = build_params(RADIUS)
    return float(np.sqrt(max(float(freq_hz) / frequency_scale(params), 0.0)))


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


def nearest_grid_value(values: np.ndarray, target: float) -> float:
    idx = int(np.argmin(np.abs(np.asarray(values, dtype=float) - float(target))))
    return float(np.asarray(values, dtype=float)[idx])


def track_branch_family(radius: float, beta_deg: float, mu_values: np.ndarray) -> dict[str, object]:
    seed_freqs, seed_vecs = solve_fem_modes(radius=radius, mu=0.0, beta_deg=0.0, n_modes=N_SOLVE)
    seed_fracs = fem_axial_fractions(seed_vecs)
    metadata = build_branch_metadata(seed_freqs=seed_freqs[:N_TRACK], seed_fracs=seed_fracs[:N_TRACK])

    prev_vecs = seed_vecs[:, :N_TRACK].copy()
    prev_freqs = seed_freqs[:N_TRACK].copy()
    prev_indices = np.arange(1, N_TRACK + 1, dtype=int)

    for beta_value in analysis_beta_grid(start=0.0, stop=beta_deg, anchor_points=(beta_deg,))[1:]:
        cur_freqs, cur_vecs = solve_fem_modes(
            radius=radius,
            mu=0.0,
            beta_deg=float(beta_value),
            n_modes=N_SOLVE,
        )
        assign, _ = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        prev_vecs = cur_vecs[:, assign]
        prev_freqs = cur_freqs[assign]
        prev_indices = assign + 1

    snapshots: list[dict[str, object]] = [
        {
            "mu": 0.0,
            "frequency_hz": prev_freqs.copy(),
            "current_sorted_index": prev_indices.copy(),
            "vecs": prev_vecs.copy(),
            "assignment_mac": np.full(N_TRACK, np.nan, dtype=float),
        }
    ]

    for mu_value in np.asarray(mu_values[1:], dtype=float):
        cur_freqs, cur_vecs = solve_fem_modes(
            radius=radius,
            mu=float(mu_value),
            beta_deg=beta_deg,
            n_modes=N_SOLVE,
        )
        assign, mac_values = assign_by_mac_and_frequency(
            prev_vecs=prev_vecs,
            prev_freqs=prev_freqs,
            cur_vecs=cur_vecs,
            cur_freqs=cur_freqs,
        )
        prev_vecs = cur_vecs[:, assign]
        prev_freqs = cur_freqs[assign]
        prev_indices = assign + 1
        snapshots.append(
            {
                "mu": float(mu_value),
                "frequency_hz": prev_freqs.copy(),
                "current_sorted_index": prev_indices.copy(),
                "vecs": prev_vecs.copy(),
                "assignment_mac": mac_values.copy(),
            }
        )

    return {"metadata": metadata, "snapshots": snapshots}


def sample_local_bending_profile(
    full_vec: np.ndarray,
    mu: float,
    beta_deg: float,
    arm: str,
    pts_per_elem: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    n = fem.N_ELEM
    beta_rad = np.deg2rad(beta_deg)
    cos_beta = float(np.cos(beta_rad))
    sin_beta = float(np.sin(beta_rad))

    if arm == "left":
        node_ids = list(range(n + 1))
        length = fem.ell * (1.0 - mu)
        le = length / n

        def local_node_state(node: int) -> tuple[float, float, float]:
            u, v, theta = full_vec[3 * node : 3 * node + 3]
            return float(u), float(v), float(theta)

    elif arm == "right":
        node_ids = list(range(n, 2 * n + 1))
        length = fem.ell * (1.0 + mu)
        le = length / n

        def local_node_state(node: int) -> tuple[float, float, float]:
            u, v, theta = full_vec[3 * node : 3 * node + 3]
            u_local = cos_beta * u + sin_beta * v
            w_local = -sin_beta * u + cos_beta * v
            return float(u_local), float(w_local), float(theta)

    else:
        raise ValueError(f"Unknown arm: {arm}")

    coord = []
    w_values = []
    for elem_idx in range(n):
        node_left = node_ids[elem_idx]
        node_right = node_ids[elem_idx + 1]
        _, w1, theta1 = local_node_state(node_left)
        _, w2, theta2 = local_node_state(node_right)
        xi_values = np.linspace(0.0, 1.0, pts_per_elem, endpoint=False)
        if elem_idx == n - 1:
            xi_values = np.append(xi_values, 1.0)
        s0 = elem_idx * le

        for xi in xi_values:
            n1 = 1.0 - 3.0 * xi**2 + 2.0 * xi**3
            n2 = le * (xi - 2.0 * xi**2 + xi**3)
            n3 = 3.0 * xi**2 - 2.0 * xi**3
            n4 = le * (-xi**2 + xi**3)
            coord.append(s0 + xi * le)
            w_values.append(n1 * w1 + n2 * theta1 + n3 * w2 + n4 * theta2)

    coord_arr = np.asarray(coord, dtype=float)
    w_arr = np.asarray(w_values, dtype=float)
    linear_trend = np.linspace(w_arr[0], w_arr[-1], len(w_arr))
    return coord_arr, w_arr - linear_trend


def count_halfwaves(w_values: np.ndarray) -> int:
    amplitude = float(np.max(np.abs(w_values)))
    if amplitude <= 1e-12:
        return 0

    tol = max(1e-8, 0.04 * amplitude)
    signs = np.zeros(len(w_values), dtype=int)
    signs[w_values > tol] = 1
    signs[w_values < -tol] = -1

    for idx in range(1, len(signs)):
        if signs[idx] == 0:
            signs[idx] = signs[idx - 1]
    for idx in range(len(signs) - 2, -1, -1):
        if signs[idx] == 0:
            signs[idx] = signs[idx + 1]

    min_segment_length = max(4, int(0.03 * len(signs)))
    halfwaves = 0
    current = signs[0]
    start = 0
    for idx in range(1, len(signs)):
        if signs[idx] != current:
            if current != 0 and idx - start >= min_segment_length:
                halfwaves += 1
            current = signs[idx]
            start = idx
    if current != 0 and len(signs) - start >= min_segment_length:
        halfwaves += 1
    return halfwaves


def nearest_reference_family(lambda_value: float, mu: float, arm: str) -> int:
    if arm == "left":
        length = fem.ell * (1.0 - mu)
    elif arm == "right":
        length = fem.ell * (1.0 + mu)
    else:
        raise ValueError(f"Unknown arm: {arm}")
    reference_lambda = REFERENCE_ROOTS * fem.ell / length
    return int(np.argmin(np.abs(reference_lambda - lambda_value))) + 1


def detect_change_windows(rows_by_branch: dict[str, list[dict[str, float | int | str]]]) -> tuple[tuple[float, float], ...]:
    centers: list[float] = []
    half_width = 0.03

    for branch_id, rows in rows_by_branch.items():
        mu_values = np.asarray([float(row["mu"]) for row in rows], dtype=float)
        lambda_values = np.asarray([float(row["Lambda"]) for row in rows], dtype=float)
        derivative = np.diff(lambda_values)
        sign = np.sign(derivative)
        for idx in range(1, len(sign)):
            if sign[idx - 1] != 0 and sign[idx] != 0 and sign[idx - 1] != sign[idx]:
                centers.append(float(0.5 * (mu_values[idx] + mu_values[idx + 1])))

        for field in (
            "sorted_position",
            "n_halfwaves_left",
            "n_halfwaves_right",
            "nearest_left_reference_family",
            "nearest_right_reference_family",
        ):
            values = np.asarray([int(row[field]) for row in rows], dtype=int)
            change_idx = np.where(values[1:] != values[:-1])[0]
            for idx in change_idx:
                centers.append(float(0.5 * (mu_values[idx] + mu_values[idx + 1])))

    windows: list[tuple[float, float]] = []
    for center in sorted(centers):
        left = max(0.0, center - half_width)
        right = min(0.9, center + half_width)
        if not windows:
            windows.append((left, right))
            continue
        prev_left, prev_right = windows[-1]
        if left <= prev_right + 1e-12:
            windows[-1] = (prev_left, max(prev_right, right))
        else:
            windows.append((left, right))
    return tuple(windows)


def build_branch_rows(tracking: dict[str, object], branch_id: str) -> list[dict[str, float | int | str]]:
    metadata = list(tracking["metadata"])
    branch_lookup = {str(item["branch_id"]): idx for idx, item in enumerate(metadata)}
    branch_pos = branch_lookup[branch_id]

    rows: list[dict[str, float | int | str]] = []
    for snapshot in tracking["snapshots"]:
        mu_value = float(snapshot["mu"])
        freq_hz = float(np.asarray(snapshot["frequency_hz"], dtype=float)[branch_pos])
        lambda_value = lambda_from_frequency_hz(freq_hz)
        slot = int(np.asarray(snapshot["current_sorted_index"], dtype=int)[branch_pos])
        mac_value = float(np.asarray(snapshot["assignment_mac"], dtype=float)[branch_pos])
        full_vec = expand_full_vector(np.asarray(snapshot["vecs"], dtype=float)[:, branch_pos])

        _, left_w = sample_local_bending_profile(full_vec=full_vec, mu=mu_value, beta_deg=BETA_DEG, arm="left")
        _, right_w = sample_local_bending_profile(full_vec=full_vec, mu=mu_value, beta_deg=BETA_DEG, arm="right")

        rows.append(
            {
                "branch_id": branch_id,
                "branch_label": BRANCH_LABELS[branch_id],
                "mu": mu_value,
                "Lambda": lambda_value,
                "sorted_position": slot,
                "n_halfwaves_left": count_halfwaves(left_w),
                "n_halfwaves_right": count_halfwaves(right_w),
                "nearest_left_reference_family": nearest_reference_family(lambda_value=lambda_value, mu=mu_value, arm="left"),
                "nearest_right_reference_family": nearest_reference_family(lambda_value=lambda_value, mu=mu_value, arm="right"),
                "tracking_mac": mac_value,
            }
        )

    rows.sort(key=lambda row: float(row["mu"]))
    return rows


def characteristic_rows(branch_id: str, rows: list[dict[str, float | int | str]]) -> list[dict[str, float | int | str]]:
    mu_values = np.asarray([float(row["mu"]) for row in rows], dtype=float)
    lambda_values = np.asarray([float(row["Lambda"]) for row in rows], dtype=float)
    row_by_mu = {float(row["mu"]): row for row in rows}

    def pick(mu_target: float) -> dict[str, float | int | str]:
        return dict(row_by_mu[nearest_grid_value(mu_values, mu_target)])

    selected: list[dict[str, float | int | str]] = []
    if branch_id == "bending_desc_02":
        plan = [
            ("малая μ", 0.0),
            ("средняя μ", 0.45),
            ("большая μ", 0.9),
        ]
    elif branch_id == "bending_desc_01":
        interior = np.arange(1, len(lambda_values) - 1, dtype=int)
        max_idx = None
        for idx in interior:
            if lambda_values[idx - 1] < lambda_values[idx] and lambda_values[idx] > lambda_values[idx + 1]:
                max_idx = int(idx)
                break
        if max_idx is None:
            plan = [
                ("малая μ", 0.0),
                ("средняя μ", 0.3),
                ("большая μ", 0.9),
            ]
        else:
            center = float(mu_values[max_idx])
            plan = [
                ("до перехода на правом плече", 0.0),
                ("вблизи максимума", center),
                ("после максимума", min(0.9, center + 0.25)),
            ]
    elif branch_id == "bending_desc_04":
        change_centers: list[float] = []
        for field in ("n_halfwaves_left", "n_halfwaves_right"):
            values = np.asarray([int(row[field]) for row in rows], dtype=int)
            idx_arr = np.where(values[1:] != values[:-1])[0]
            for idx in idx_arr:
                change_centers.append(float(0.5 * (mu_values[idx] + mu_values[idx + 1])))
        if change_centers:
            first = min(change_centers)
            last = max(change_centers)
            center = 0.5 * (first + last)
            plan = [
                ("до перестройки", max(0.0, first - 0.12)),
                ("вблизи перестройки", center),
                ("после перестройки", min(0.9, last + 0.12)),
            ]
        else:
            plan = [
                ("малая μ", 0.0),
                ("средняя μ", 0.3),
                ("большая μ", 0.9),
            ]
    else:
        raise ValueError(f"Unknown branch: {branch_id}")

    used_mu: set[float] = set()
    for label, mu_target in plan:
        row = pick(mu_target)
        if float(row["mu"]) in used_mu:
            continue
        used_mu.add(float(row["mu"]))
        row["qualitative_behavior"] = label
        selected.append(row)

    return selected


@dataclass
class ShapeCase:
    row: dict[str, float | int | str]
    x_base: np.ndarray
    y_base: np.ndarray
    u_norm: np.ndarray
    v_norm: np.ndarray


def build_shape_plot(
    branch_id: str,
    selected_rows: list[dict[str, float | int | str]],
    tracking: dict[str, object],
) -> Path:
    metadata = list(tracking["metadata"])
    branch_lookup = {str(item["branch_id"]): idx for idx, item in enumerate(metadata)}
    branch_pos = branch_lookup[branch_id]
    ell = build_params(RADIUS).L_total / 2.0
    snapshots = {float(item["mu"]): item for item in tracking["snapshots"]}

    cases: list[ShapeCase] = []
    for row in selected_rows:
        mu_value = float(row["mu"])
        snapshot = snapshots[mu_value]
        full_vec = expand_full_vector(np.asarray(snapshot["vecs"], dtype=float)[:, branch_pos])
        x_base, y_base = build_node_coordinates(mu=mu_value, beta_deg=BETA_DEG, ell=ell)
        u = full_vec[0::3]
        v = full_vec[1::3]
        norm = max(float(np.max(np.sqrt(u * u + v * v))), 1e-12)
        cases.append(
            ShapeCase(
                row=row,
                x_base=x_base,
                y_base=y_base,
                u_norm=u / norm,
                v_norm=v / norm,
            )
        )

    fig, axes = plt.subplots(1, len(cases), figsize=(5.6 * len(cases), 4.7), squeeze=False)
    axes_row = axes[0]

    x_all = []
    y_all = []
    for case in cases:
        x_def = case.x_base + MODE_SHAPE_SCALE * case.u_norm
        y_def = case.y_base + MODE_SHAPE_SCALE * case.v_norm
        x_all.extend([float(np.min(case.x_base)), float(np.max(case.x_base)), float(np.min(x_def)), float(np.max(x_def))])
        y_all.extend([float(np.min(case.y_base)), float(np.max(case.y_base)), float(np.min(y_def)), float(np.max(y_def))])
    x_margin = 0.06 * max(max(x_all) - min(x_all), 1.0)
    y_margin = 0.10 * max(max(y_all) - min(y_all), 0.5)
    x_limits = (min(x_all) - x_margin, max(x_all) + x_margin)
    y_limits = (min(y_all) - y_margin, max(y_all) + y_margin)

    for ax, case in zip(axes_row, cases):
        joint = fem.N_ELEM
        x_def = case.x_base + MODE_SHAPE_SCALE * case.u_norm
        y_def = case.y_base + MODE_SHAPE_SCALE * case.v_norm
        ax.plot(case.x_base[: joint + 1], case.y_base[: joint + 1], color="0.78", linestyle="--", linewidth=1.1)
        ax.plot(case.x_base[joint:], case.y_base[joint:], color="0.78", linestyle="--", linewidth=1.1)
        ax.plot(x_def[: joint + 1], y_def[: joint + 1], color="#1f77b4", linewidth=2.2)
        ax.plot(x_def[joint:], y_def[joint:], color="#ff7f0e", linewidth=2.2)
        ax.scatter([case.x_base[joint]], [case.y_base[joint]], color="black", s=14, zorder=5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.grid(True, alpha=0.18)
        ax.set_title(
            f"μ = {float(case.row['mu']):.4f}, Λ = {float(case.row['Lambda']):.3f}\n"
            f"{str(case.row['qualitative_behavior'])}, место в спектре: {int(case.row['sorted_position'])}",
            fontsize=10.0,
        )
        ax.text(
            0.02,
            0.02,
            (
                f"полуволн слева: {int(case.row['n_halfwaves_left'])}\n"
                f"полуволн справа: {int(case.row['n_halfwaves_right'])}\n"
                f"семья слева: {int(case.row['nearest_left_reference_family'])}\n"
                f"семья справа: {int(case.row['nearest_right_reference_family'])}"
            ),
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.7,
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.95},
        )
        ax.set_xlabel("x")
        ax.set_ylabel("y")

    handles = [
        Line2D([0], [0], color="0.78", linestyle="--", linewidth=1.1, label="недеформированная геометрия"),
        Line2D([0], [0], color="#1f77b4", linewidth=2.2, label="левое плечо"),
        Line2D([0], [0], color="#ff7f0e", linewidth=2.2, label="правое плечо"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.98), fontsize=10)
    fig.suptitle(f"β = 15°, r = 5 мм, {BRANCH_LABELS[branch_id]}", y=0.995, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.92))

    out_path = RESULTS_DIR / SHAPE_PATH_TEMPLATE.format(branch_id=branch_id)
    fig.savefig(out_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return out_path


def attach_branch_behavior(rows_by_branch: dict[str, list[dict[str, float | int | str]]]) -> None:
    for branch_id, rows in rows_by_branch.items():
        lambda_values = np.asarray([float(row["Lambda"]) for row in rows], dtype=float)
        if np.all(np.diff(lambda_values) <= 1e-10):
            behavior = "монотонно убывает"
        else:
            diffs = np.diff(lambda_values)
            signs = np.sign(diffs)
            change_points: list[tuple[int, int]] = []
            for idx in range(1, len(signs)):
                if signs[idx - 1] != 0 and signs[idx] != 0 and signs[idx - 1] != signs[idx]:
                    change_points.append((idx, int(signs[idx - 1]), int(signs[idx])))

            if len(change_points) == 1:
                idx, left_sign, right_sign = change_points[0]
                mu_center = 0.5 * (float(rows[idx]["mu"]) + float(rows[idx + 1]["mu"]))
                if left_sign > 0 and right_sign < 0:
                    behavior = f"сначала возрастает, затем убывает; максимум около μ={mu_center:.4f}"
                elif left_sign < 0 and right_sign > 0:
                    behavior = f"сначала убывает, затем возрастает; минимум около μ={mu_center:.4f}"
                else:
                    behavior = "немонотонное поведение"
            elif len(change_points) >= 2:
                first_idx, first_left, first_right = change_points[0]
                second_idx, second_left, second_right = change_points[1]
                first_mu = 0.5 * (float(rows[first_idx]["mu"]) + float(rows[first_idx + 1]["mu"]))
                second_mu = 0.5 * (float(rows[second_idx]["mu"]) + float(rows[second_idx + 1]["mu"]))
                if first_left < 0 and first_right > 0 and second_left > 0 and second_right < 0:
                    behavior = (
                        f"сначала убывает, затем возрастает; минимум около μ={first_mu:.4f}, "
                        f"последующий максимум около μ={second_mu:.4f}"
                    )
                else:
                    behavior = "немонотонное поведение с несколькими перестройками"
            else:
                behavior = "без внутреннего экстремума на 0≤μ≤0.9"
        for row in rows:
            row["branch_behavior"] = behavior


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)

    coarse_mu = analysis_mu_grid(start=0.0, stop=0.9, anchor_points=(0.45, 0.9))
    coarse_tracking = track_branch_family(radius=RADIUS, beta_deg=BETA_DEG, mu_values=coarse_mu)
    coarse_rows = {branch_id: build_branch_rows(coarse_tracking, branch_id) for branch_id in TARGET_BRANCH_IDS}
    refinement_windows = detect_change_windows(coarse_rows)

    refined_mu = analysis_mu_grid(
        start=0.0,
        stop=0.9,
        local_windows=refinement_windows,
        anchor_points=(0.45, 0.9),
    )
    tracking = track_branch_family(radius=RADIUS, beta_deg=BETA_DEG, mu_values=refined_mu)
    rows_by_branch = {branch_id: build_branch_rows(tracking, branch_id) for branch_id in TARGET_BRANCH_IDS}
    attach_branch_behavior(rows_by_branch)

    full_rows: list[dict[str, float | int | str]] = []
    selected_rows: list[dict[str, float | int | str]] = []
    shape_paths: list[Path] = []

    for branch_id in sorted(TARGET_BRANCH_IDS, key=branch_sort_key):
        branch_rows = rows_by_branch[branch_id]
        full_rows.extend(branch_rows)
        chosen_rows = characteristic_rows(branch_id=branch_id, rows=branch_rows)
        selected_rows.extend(chosen_rows)
        shape_paths.append(build_shape_plot(branch_id=branch_id, selected_rows=chosen_rows, tracking=tracking))

    write_csv_rows(FULL_CSV_PATH, full_rows)
    write_csv_rows(SELECTED_CSV_PATH, selected_rows)

    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu base step: {ANALYSIS_MU_STEP:.4f}")
    print(f"analysis mu local step: {LOCAL_MU_REFINEMENT_STEP:.4f}")
    if refinement_windows:
        print("local refinement windows:")
        for left, right in refinement_windows:
            print(f"  {left:.4f} .. {right:.4f}")
    else:
        print("local refinement windows: none")
    print(f"saved CSV: {FULL_CSV_PATH}")
    print(f"saved CSV: {SELECTED_CSV_PATH}")
    for path in shape_paths:
        print(f"saved PNG: {path}")


if __name__ == "__main__":
    main()
