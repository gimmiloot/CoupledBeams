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
BENDING_DESCENDANT_PREFIX = "bending_desc_"
PLOT_KINDS = ("full", "transverse", "components")
NORMALIZE_KINDS = ("auto", "max-full", "max-transverse", "none")
NEAR_ZERO_NORM = 1e-12
DEFAULT_SINGLE_GEOMETRY_FIGSIZE = (8.0, 4.2)
DEFAULT_SINGLE_COMPONENTS_FIGSIZE = (8.0, 5.4)


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


def default_single_output_path(
    branch_id: str,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    plot_kind: str = "full",
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "auto",
) -> Path:
    beta_token = filename_number_token(beta_deg)
    mu_token = filename_number_token(mu_value)
    epsilon_token = filename_number_token(epsilon)
    scale_token = filename_number_token(mode_scale)
    normalize_token = "" if normalize == "auto" else f"_norm{normalize.replace('-', '_')}"
    return RESULTS_DIR / (
        f"tracked_bending_descendant_shape_{plot_kind}_beta{beta_token}_{branch_id}_mu{mu_token}"
        f"_eps{epsilon_token}_scale{scale_token}{normalize_token}_ru.png"
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
    plot_kind: str = "full",
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "auto",
) -> Path:
    if value is None:
        return default_single_output_path(
            branch_id=branch_id,
            beta_deg=beta_deg,
            mu_value=mu_value,
            epsilon=epsilon,
            plot_kind=plot_kind,
            mode_scale=mode_scale,
            normalize=normalize,
        )
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def infer_bending_descendant_title_label(branch_id: str) -> str:
    if branch_id.startswith(BENDING_DESCENDANT_PREFIX):
        number_part = branch_id.removeprefix(BENDING_DESCENDANT_PREFIX)
        if number_part.isdecimal():
            number = int(number_part)
            return f"потомок {number}-й изгибной ветви"
    return f"ветвь {branch_id}"


def infer_title_label(branch_id: str) -> str:
    return infer_bending_descendant_title_label(branch_id)


def lambda_from_frequency_hz(freq_hz: float, radius: float) -> float:
    params = build_params(float(radius))
    return float(np.sqrt(max(float(freq_hz) / frequency_scale(params), 0.0)))


def radius_from_epsilon(epsilon: float, l_total: float | None = None) -> float:
    ell = (build_params(float(DEFAULT_TARGET_RADII[0])).L_total if l_total is None else float(l_total)) / 2.0
    return 2.0 * ell * float(epsilon)


def resolve_normalization(plot_kind: str, normalize: str | None) -> str:
    if normalize is not None and normalize != "auto":
        return normalize
    if plot_kind == "transverse":
        return "max-transverse"
    return "max-full"


def normalization_factor(shape_case: dict[str, object], normalize: str) -> float:
    if normalize == "none":
        return 1.0
    if normalize == "max-transverse":
        value = float(shape_case["max_abs_transverse"])
    elif normalize == "max-full":
        value = float(shape_case["max_full_displacement"])
    else:
        raise ValueError(f"Unknown normalization mode: {normalize}")
    return max(value, NEAR_ZERO_NORM)


def compute_local_components(
    u_global: np.ndarray,
    v_global: np.ndarray,
    beta_deg: float,
) -> dict[str, np.ndarray]:
    joint = fem.N_ELEM
    beta_rad = np.deg2rad(beta_deg)
    cos_beta = float(np.cos(beta_rad))
    sin_beta = float(np.sin(beta_rad))
    u_right_global = u_global[joint:]
    v_right_global = v_global[joint:]
    return {
        "u_left_local": u_global[: joint + 1],
        "w_left_local": v_global[: joint + 1],
        "u_right_local": cos_beta * u_right_global + sin_beta * v_right_global,
        "w_right_local": -sin_beta * u_right_global + cos_beta * v_right_global,
    }


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
    local = compute_local_components(u_global=u, v_global=v, beta_deg=case_beta)
    full_displacement = np.sqrt(u * u + v * v)
    max_disp = float(np.max(full_displacement))
    max_abs_transverse = float(
        max(
            np.max(np.abs(local["w_left_local"])),
            np.max(np.abs(local["w_right_local"])),
        )
    )
    max_abs_axial = float(
        max(
            np.max(np.abs(local["u_left_local"])),
            np.max(np.abs(local["u_right_local"])),
        )
    )
    norm = max(max_disp, NEAR_ZERO_NORM)
    frequency_hz = float(snapshot["frequency_hz"][branch_pos])
    axial_fraction_values = snapshot.get("axial_fraction")
    axial_fraction = float(axial_fraction_values[branch_pos]) if axial_fraction_values is not None else np.nan
    full_for_ratio = max(max_disp, NEAR_ZERO_NORM)

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
        "u_raw": u,
        "v_raw": v,
        "u": u / norm,
        "v": v / norm,
        "u_left_local": local["u_left_local"],
        "w_left_local": local["w_left_local"],
        "u_right_local": local["u_right_local"],
        "w_right_local": local["w_right_local"],
        "s_left_normalized": np.linspace(0.0, 1.0, fem.N_ELEM + 1),
        "s_right_normalized": np.linspace(0.0, 1.0, fem.N_ELEM + 1),
        "max_full_displacement": max_disp,
        "max_abs_transverse": max_abs_transverse,
        "max_abs_axial": max_abs_axial,
        "transverse_to_full_ratio": max_abs_transverse / full_for_ratio,
        "axial_to_full_ratio": max_abs_axial / full_for_ratio,
        "axial_fraction": axial_fraction,
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


def full_deformed_coordinates(
    shape_case: dict[str, object],
    *,
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "max-full",
) -> tuple[np.ndarray, np.ndarray]:
    x_base = np.asarray(shape_case["x_base"], dtype=float)
    y_base = np.asarray(shape_case["y_base"], dtype=float)
    norm = normalization_factor(shape_case=shape_case, normalize=normalize)
    u = np.asarray(shape_case["u_raw"], dtype=float) / norm
    v = np.asarray(shape_case["v_raw"], dtype=float) / norm
    return x_base + mode_scale * u, y_base + mode_scale * v


def transverse_deformed_arm_coordinates(
    shape_case: dict[str, object],
    *,
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "max-transverse",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    joint = fem.N_ELEM
    x_base = np.asarray(shape_case["x_base"], dtype=float)
    y_base = np.asarray(shape_case["y_base"], dtype=float)
    beta_rad = np.deg2rad(float(shape_case["beta"]))
    normal_x = -float(np.sin(beta_rad))
    normal_y = float(np.cos(beta_rad))
    norm = normalization_factor(shape_case=shape_case, normalize=normalize)
    w_left = np.asarray(shape_case["w_left_local"], dtype=float) / norm
    w_right = np.asarray(shape_case["w_right_local"], dtype=float) / norm

    x_left = x_base[: joint + 1].copy()
    y_left = y_base[: joint + 1] + mode_scale * w_left
    x_right = x_base[joint:] + mode_scale * w_right * normal_x
    y_right = y_base[joint:] + mode_scale * w_right * normal_y
    return x_left, y_left, x_right, y_right


def deformed_arm_coordinates(
    shape_case: dict[str, object],
    *,
    plot_kind: str = "full",
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "max-full",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    joint = fem.N_ELEM
    if plot_kind == "full":
        x_def, y_def = full_deformed_coordinates(
            shape_case=shape_case,
            mode_scale=mode_scale,
            normalize=normalize,
        )
        return x_def[: joint + 1], y_def[: joint + 1], x_def[joint:], y_def[joint:]
    if plot_kind == "transverse":
        return transverse_deformed_arm_coordinates(
            shape_case=shape_case,
            mode_scale=mode_scale,
            normalize=normalize,
        )
    raise ValueError(f"Geometry deformation is not defined for plot kind: {plot_kind}")


def deformed_coordinates(
    shape_case: dict[str, object],
    *,
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "max-full",
) -> tuple[np.ndarray, np.ndarray]:
    return full_deformed_coordinates(shape_case=shape_case, mode_scale=mode_scale, normalize=normalize)


def shape_axis_limits(
    shape_cases: Sequence[dict[str, object]],
    *,
    plot_kind: str = "full",
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "max-full",
) -> tuple[tuple[float, float], tuple[float, float]]:
    x_all: list[float] = []
    y_all: list[float] = []
    for shape_case in shape_cases:
        x_base = np.asarray(shape_case["x_base"], dtype=float)
        y_base = np.asarray(shape_case["y_base"], dtype=float)
        x_left, y_left, x_right, y_right = deformed_arm_coordinates(
            shape_case=shape_case,
            plot_kind=plot_kind,
            mode_scale=mode_scale,
            normalize=normalize,
        )
        x_all.extend(
            [
                float(np.min(x_base)),
                float(np.max(x_base)),
                float(np.min(x_left)),
                float(np.max(x_left)),
                float(np.min(x_right)),
                float(np.max(x_right)),
            ]
        )
        y_all.extend(
            [
                float(np.min(y_base)),
                float(np.max(y_base)),
                float(np.min(y_left)),
                float(np.max(y_left)),
                float(np.min(y_right)),
                float(np.max(y_right)),
            ]
        )

    x_margin = 0.06 * max(max(x_all) - min(x_all), 1.0)
    y_margin = 0.10 * max(max(y_all) - min(y_all), 0.5)
    return (min(x_all) - x_margin, max(x_all) + x_margin), (min(y_all) - y_margin, max(y_all) + y_margin)


def draw_shape_case(
    ax: Axes,
    shape_case: dict[str, object],
    *,
    x_limits: tuple[float, float] | None = None,
    y_limits: tuple[float, float] | None = None,
    plot_kind: str = "full",
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str = "max-full",
) -> None:
    joint = fem.N_ELEM
    x_base = np.asarray(shape_case["x_base"], dtype=float)
    y_base = np.asarray(shape_case["y_base"], dtype=float)
    x_left, y_left, x_right, y_right = deformed_arm_coordinates(
        shape_case=shape_case,
        plot_kind=plot_kind,
        mode_scale=mode_scale,
        normalize=normalize,
    )

    ax.plot(x_base[: joint + 1], y_base[: joint + 1], color="0.78", linestyle="--", linewidth=1.1)
    ax.plot(x_base[joint:], y_base[joint:], color="0.78", linestyle="--", linewidth=1.1)
    ax.plot(x_left, y_left, color="#1f77b4", linewidth=2.2)
    ax.plot(x_right, y_right, color="#ff7f0e", linewidth=2.2)
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


def single_shape_title(
    shape_case: dict[str, object],
    *,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    title_label: str,
    plot_kind: str,
) -> str:
    prefix = "Локальные компоненты: " if plot_kind == "components" else ""
    return (
        f"{prefix}β = {beta_deg:g}°, μ = {mu_label(mu_value)}, ε = {epsilon:g}\n"
        f"{title_label}, Λ = {float(shape_case['lambda']):.3f}, "
        f"место в спектре: {int(shape_case['current_sorted_index'])}"
    )


def normalized_component_curves(
    shape_case: dict[str, object],
    *,
    normalize: str,
) -> dict[str, np.ndarray]:
    norm = normalization_factor(shape_case=shape_case, normalize=normalize)
    return {
        "left_axial": np.asarray(shape_case["u_left_local"], dtype=float) / norm,
        "left_transverse": np.asarray(shape_case["w_left_local"], dtype=float) / norm,
        "right_axial": np.asarray(shape_case["u_right_local"], dtype=float) / norm,
        "right_transverse": np.asarray(shape_case["w_right_local"], dtype=float) / norm,
    }


def build_component_diagnostic_plot(
    shape_case: dict[str, object],
    *,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    output_path: Path,
    title_label: str,
    normalize: str,
    dpi: int = 240,
    figsize: tuple[float, float] | None = None,
    show: bool = False,
) -> Path:
    curves = normalized_component_curves(shape_case=shape_case, normalize=normalize)
    s_left = np.asarray(shape_case["s_left_normalized"], dtype=float)
    s_right = np.asarray(shape_case["s_right_normalized"], dtype=float)

    fig, axes = plt.subplots(2, 1, figsize=figsize or DEFAULT_SINGLE_COMPONENTS_FIGSIZE, sharex=True)
    axes[0].plot(s_left, curves["left_axial"], color="#1f77b4", linewidth=2.0, label="левое плечо")
    axes[0].plot(s_right, curves["right_axial"], color="#ff7f0e", linewidth=2.0, label="правое плечо")
    axes[0].set_ylabel("локальная продольная")
    axes[0].grid(True, alpha=0.22)
    axes[0].legend(fontsize=9, loc="best")

    axes[1].plot(s_left, curves["left_transverse"], color="#1f77b4", linewidth=2.0, label="левое плечо")
    axes[1].plot(s_right, curves["right_transverse"], color="#ff7f0e", linewidth=2.0, label="правое плечо")
    axes[1].set_xlabel("локальная координата s / L")
    axes[1].set_ylabel("локальная поперечная")
    axes[1].grid(True, alpha=0.22)
    axes[1].legend(fontsize=9, loc="best")

    fig.suptitle(
        single_shape_title(
            shape_case=shape_case,
            beta_deg=beta_deg,
            mu_value=mu_value,
            epsilon=epsilon,
            title_label=title_label,
            plot_kind="components",
        ),
        fontsize=10.6,
    )
    fig.tight_layout(pad=0.8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return output_path


def single_shape_diagnostics(
    shape_case: dict[str, object],
    *,
    plot_kind: str,
    mode_scale: float,
    normalize: str,
) -> dict[str, float | int | str]:
    return {
        "branch_id": str(shape_case["branch_id"]),
        "beta": float(shape_case["beta"]),
        "mu": float(shape_case["mu"]),
        "radius": float(shape_case["r"]),
        "lambda": float(shape_case["lambda"]),
        "current_sorted_index": int(shape_case["current_sorted_index"]),
        "plot_kind": plot_kind,
        "mode_scale": float(mode_scale),
        "normalization": normalize,
        "max_full_displacement": float(shape_case["max_full_displacement"]),
        "max_abs_transverse": float(shape_case["max_abs_transverse"]),
        "max_abs_axial": float(shape_case["max_abs_axial"]),
        "transverse_to_full_ratio": float(shape_case["transverse_to_full_ratio"]),
        "axial_to_full_ratio": float(shape_case["axial_to_full_ratio"]),
        "axial_fraction": float(shape_case["axial_fraction"]),
    }


def build_single_shape_plot(
    branch_id: str,
    beta_deg: float,
    mu_value: float,
    radius: float,
    epsilon: float,
    output_path: Path,
    title_label: str,
    plot_kind: str = "full",
    mode_scale: float = MODE_SHAPE_SCALE,
    normalize: str | None = None,
    shape_case: dict[str, object] | None = None,
    dpi: int = 240,
    figsize: tuple[float, float] | None = None,
    show: bool = False,
) -> Path:
    if plot_kind not in PLOT_KINDS:
        raise ValueError(f"Unknown plot kind: {plot_kind}")
    resolved_normalize = resolve_normalization(plot_kind=plot_kind, normalize=normalize)
    if shape_case is None:
        shape_case = collect_single_branch_shape(
            branch_id=branch_id,
            beta_deg=beta_deg,
            mu_value=mu_value,
            radius=radius,
        )

    if plot_kind == "components":
        return build_component_diagnostic_plot(
            shape_case=shape_case,
            beta_deg=beta_deg,
            mu_value=mu_value,
            epsilon=epsilon,
            output_path=output_path,
            title_label=title_label,
            normalize=resolved_normalize,
            dpi=dpi,
            figsize=figsize,
            show=show,
        )

    x_limits, y_limits = shape_axis_limits(
        [shape_case],
        plot_kind=plot_kind,
        mode_scale=mode_scale,
        normalize=resolved_normalize,
    )
    fig, ax = plt.subplots(figsize=figsize or DEFAULT_SINGLE_GEOMETRY_FIGSIZE)
    draw_shape_case(
        ax=ax,
        shape_case=shape_case,
        x_limits=x_limits,
        y_limits=y_limits,
        plot_kind=plot_kind,
        mode_scale=mode_scale,
        normalize=resolved_normalize,
    )
    ax.set_xlabel("координата x")
    ax.set_ylabel("координата y")
    ax.set_title(
        single_shape_title(
            shape_case=shape_case,
            beta_deg=beta_deg,
            mu_value=mu_value,
            epsilon=epsilon,
            title_label=title_label,
            plot_kind=plot_kind,
        ),
        fontsize=10.6,
    )

    fig.tight_layout(pad=0.8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return output_path
