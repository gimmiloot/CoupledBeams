from __future__ import annotations

from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import frequency_scale  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.analyze_flat_mu_bending_energy import (  # noqa: E402
    MODE_SHAPE_SCALE,
    N_SOLVE,
    N_TRACK,
    RESULTS_DIR,
    TARGET_RADII as DEFAULT_TARGET_RADII,
    assign_by_mac_and_frequency,
    build_branch_metadata,
    build_node_coordinates,
    build_params,
    expand_full_vector,
    mu_key,
    solve_fem_modes,
)
from scripts.compare_beta0_analytic_vs_fem import fem_axial_fractions  # noqa: E402
from scripts.sweep_grid_policy import ANALYSIS_BETA_STEP, ANALYSIS_MU_STEP, analysis_beta_grid, analysis_mu_grid  # noqa: E402


DEFAULT_BETA_DEG = 15.0
DEFAULT_BRANCH_ID = "bending_desc_01"
DEFAULT_TARGET_MUS = (0.0, 0.1, 0.2)
TITLE_LABELS_RU = {
    "bending_desc_01": "потомок 1-й изгибной ветви",
    "bending_desc_02": "потомок 2-й изгибной ветви",
    "bending_desc_04": "потомок 4-й изгибной ветви",
    "bending_desc_06": "потомок 6-й изгибной ветви",
}


def mu_label(mu: float) -> str:
    rounded = round(float(mu))
    if abs(float(mu) - rounded) < 1e-9:
        return str(int(rounded))
    return f"{float(mu):.1f}".rstrip("0").rstrip(".")


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def target_mus_token(target_mus: Sequence[float]) -> str:
    return "mu" + "_".join(filename_number_token(mu) for mu in target_mus)


def default_output_path(branch_id: str, beta_deg: float, target_mus: Sequence[float]) -> Path:
    beta_token = filename_number_token(beta_deg)
    return RESULTS_DIR / f"tracked_bending_descendant_shapes_beta{beta_token}_{branch_id}_{target_mus_token(target_mus)}_ru.png"


def default_single_output_path(branch_id: str, beta_deg: float, mu_value: float, epsilon: float) -> Path:
    beta_token = filename_number_token(beta_deg)
    mu_token = filename_number_token(mu_value)
    epsilon_token = filename_number_token(epsilon)
    return RESULTS_DIR / (
        f"tracked_bending_descendant_shape_beta{beta_token}_{branch_id}_mu{mu_token}_eps{epsilon_token}_ru.png"
    )


def resolve_output_path(value: str | None, branch_id: str, beta_deg: float, target_mus: Sequence[float]) -> Path:
    if value is None:
        return default_output_path(branch_id=branch_id, beta_deg=beta_deg, target_mus=target_mus)
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def resolve_single_output_path(
    value: str | None,
    branch_id: str,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
) -> Path:
    if value is None:
        return default_single_output_path(
            branch_id=branch_id,
            beta_deg=beta_deg,
            mu_value=mu_value,
            epsilon=epsilon,
        )
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def infer_title_label(branch_id: str) -> str:
    return TITLE_LABELS_RU.get(branch_id, f"ветвь {branch_id}")


def lambda_from_frequency_hz(freq_hz: float, radius: float) -> float:
    params = build_params(float(radius))
    return float(np.sqrt(max(float(freq_hz) / frequency_scale(params), 0.0)))


def radius_from_epsilon(epsilon: float) -> float:
    params = build_params(float(DEFAULT_TARGET_RADII[0]))
    ell = params.L_total / 2.0
    return 2.0 * ell * float(epsilon)


def track_branch_family(radius: float, beta_deg: float, target_mus: Sequence[float]) -> dict[str, object]:
    mu_track_values = analysis_mu_grid(start=0.0, stop=max(target_mus), anchor_points=tuple(target_mus))
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

    for mu_value in mu_track_values[1:]:
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


def branch_position(tracking: dict[str, object], branch_id: str) -> int:
    metadata = list(tracking["metadata"])
    branch_lookup = {str(meta["branch_id"]): idx for idx, meta in enumerate(metadata)}
    if branch_id not in branch_lookup:
        raise KeyError(f"Tracked metadata does not contain {branch_id}.")
    return branch_lookup[branch_id]


def shape_case_from_tracking(
    tracking: dict[str, object],
    branch_id: str,
    mu_value: float,
    *,
    radius: float | None = None,
    beta_deg: float | None = None,
) -> dict[str, object]:
    case_radius = float(radius if radius is not None else tracking["radius"])
    case_beta = float(beta_deg if beta_deg is not None else tracking["beta"])
    branch_pos = branch_position(tracking=tracking, branch_id=branch_id)
    snapshot = tracking["snapshots"][mu_key(mu_value)]
    full_vec = expand_full_vector(np.asarray(snapshot["vecs"][:, branch_pos], dtype=float))
    ell = build_params(case_radius).L_total / 2.0
    x_base, y_base = build_node_coordinates(mu=float(mu_value), beta_deg=case_beta, ell=ell)
    u = full_vec[0::3]
    v = full_vec[1::3]
    max_disp = float(np.max(np.sqrt(u * u + v * v)))
    norm = max(max_disp, 1e-12)
    frequency_hz = float(snapshot["frequency_hz"][branch_pos])

    return {
        "beta": case_beta,
        "branch_id": branch_id,
        "r": case_radius,
        "mu": float(mu_value),
        "current_sorted_index": int(snapshot["current_sorted_index"][branch_pos]),
        "frequency_hz": frequency_hz,
        "lambda": lambda_from_frequency_hz(freq_hz=frequency_hz, radius=case_radius),
        "x_base": x_base,
        "y_base": y_base,
        "u": u / norm,
        "v": v / norm,
    }


def collect_single_branch_shape(
    branch_id: str,
    beta_deg: float,
    mu_value: float,
    radius: float,
) -> dict[str, object]:
    tracking = track_branch_family(radius=float(radius), beta_deg=beta_deg, target_mus=(float(mu_value),))
    return shape_case_from_tracking(
        tracking=tracking,
        branch_id=branch_id,
        mu_value=float(mu_value),
        radius=float(radius),
        beta_deg=beta_deg,
    )


def collect_branch_rows(
    branch_id: str,
    beta_deg: float,
    target_mus: Sequence[float],
    radii: Sequence[float],
) -> tuple[list[dict[str, float | int | str]], dict[tuple[float, float], dict[str, object]]]:
    rows: list[dict[str, float | int | str]] = []
    tracking_cache: dict[tuple[float, float], dict[str, object]] = {}

    for radius in radii:
        key = (float(radius), beta_deg)
        tracking_cache[key] = track_branch_family(radius=float(radius), beta_deg=beta_deg, target_mus=target_mus)
        tracking = tracking_cache[key]
        branch_pos = branch_position(tracking=tracking, branch_id=branch_id)

        for mu_value in target_mus:
            snapshot = tracking["snapshots"][mu_key(mu_value)]
            rows.append(
                {
                    "beta": beta_deg,
                    "branch_id": branch_id,
                    "r": float(radius),
                    "mu": float(mu_value),
                    "current_sorted_index": int(snapshot["current_sorted_index"][branch_pos]),
                    "frequency_hz": float(snapshot["frequency_hz"][branch_pos]),
                }
            )

    rows.sort(key=lambda row: (float(row["r"]), float(row["mu"])))
    return rows, tracking_cache


def collect_branch_shape_cases(
    branch_id: str,
    beta_deg: float,
    target_mus: Sequence[float],
    radii: Sequence[float],
) -> list[dict[str, object]]:
    cases: list[dict[str, object]] = []
    for radius in radii:
        tracking = track_branch_family(radius=float(radius), beta_deg=beta_deg, target_mus=target_mus)
        for mu_value in target_mus:
            cases.append(
                shape_case_from_tracking(
                    tracking=tracking,
                    branch_id=branch_id,
                    mu_value=float(mu_value),
                    radius=float(radius),
                    beta_deg=beta_deg,
                )
            )
    cases.sort(key=lambda case: (float(case["r"]), float(case["mu"])))
    return cases


def deformed_coordinates(shape_case: dict[str, object]) -> tuple[np.ndarray, np.ndarray]:
    x_base = np.asarray(shape_case["x_base"], dtype=float)
    y_base = np.asarray(shape_case["y_base"], dtype=float)
    u = np.asarray(shape_case["u"], dtype=float)
    v = np.asarray(shape_case["v"], dtype=float)
    return x_base + MODE_SHAPE_SCALE * u, y_base + MODE_SHAPE_SCALE * v


def shape_axis_limits(shape_cases: Sequence[dict[str, object]]) -> tuple[tuple[float, float], tuple[float, float]]:
    x_all: list[float] = []
    y_all: list[float] = []
    for shape_case in shape_cases:
        x_base = np.asarray(shape_case["x_base"], dtype=float)
        y_base = np.asarray(shape_case["y_base"], dtype=float)
        x_def, y_def = deformed_coordinates(shape_case)
        x_all.extend([float(np.min(x_base)), float(np.max(x_base)), float(np.min(x_def)), float(np.max(x_def))])
        y_all.extend([float(np.min(y_base)), float(np.max(y_base)), float(np.min(y_def)), float(np.max(y_def))])

    x_margin = 0.06 * max(max(x_all) - min(x_all), 1.0)
    y_margin = 0.10 * max(max(y_all) - min(y_all), 0.5)
    return (min(x_all) - x_margin, max(x_all) + x_margin), (min(y_all) - y_margin, max(y_all) + y_margin)


def draw_shape_case(
    ax: Axes,
    shape_case: dict[str, object],
    *,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
) -> None:
    joint = fem.N_ELEM
    x_base = np.asarray(shape_case["x_base"], dtype=float)
    y_base = np.asarray(shape_case["y_base"], dtype=float)
    x_def, y_def = deformed_coordinates(shape_case)

    ax.plot(x_base[: joint + 1], y_base[: joint + 1], color="0.78", linestyle="--", linewidth=1.1)
    ax.plot(x_base[joint:], y_base[joint:], color="0.78", linestyle="--", linewidth=1.1)
    ax.plot(x_def[: joint + 1], y_def[: joint + 1], color="#1f77b4", linewidth=2.2)
    ax.plot(x_def[joint:], y_def[joint:], color="#ff7f0e", linewidth=2.2)
    ax.scatter([x_base[joint]], [y_base[joint]], color="black", s=12, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    if x_limits is not None:
        ax.set_xlim(*x_limits)
    if y_limits is not None:
        ax.set_ylim(*y_limits)
    ax.grid(True, alpha=0.18)


def shape_legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color="0.78", linestyle="--", linewidth=1.1, label="недеформированная геометрия"),
        Line2D([0], [0], color="#1f77b4", linewidth=2.2, label="левое плечо"),
        Line2D([0], [0], color="#ff7f0e", linewidth=2.2, label="правое плечо"),
    ]


def build_single_shape_plot(
    branch_id: str,
    beta_deg: float,
    mu_value: float,
    radius: float,
    epsilon: float,
    output_path: Path,
    title_label: str,
) -> Path:
    shape_case = collect_single_branch_shape(
        branch_id=branch_id,
        beta_deg=beta_deg,
        mu_value=mu_value,
        radius=radius,
    )
    x_limits, y_limits = shape_axis_limits([shape_case])
    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    draw_shape_case(ax=ax, shape_case=shape_case, x_limits=x_limits, y_limits=y_limits)
    ax.set_xlabel("координата x")
    ax.set_ylabel("координата y")
    ax.set_title(
        (
            f"β = {beta_deg:g}°, μ = {mu_label(mu_value)}, ε = {epsilon:g}, r = {radius * 1e3:g} мм\n"
            f"{branch_id}, {title_label}, Λ = {float(shape_case['lambda']):.3f}, "
            f"место в спектре: {int(shape_case['current_sorted_index'])}"
        ),
        fontsize=10.6,
    )

    fig.legend(handles=shape_legend_handles(), loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.985), fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return output_path
