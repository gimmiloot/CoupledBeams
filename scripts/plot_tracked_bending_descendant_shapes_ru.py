from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Sequence

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


def resolve_output_path(value: str | None, branch_id: str, beta_deg: float, target_mus: Sequence[float]) -> Path:
    if value is None:
        return default_output_path(branch_id=branch_id, beta_deg=beta_deg, target_mus=target_mus)
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def infer_title_label(branch_id: str) -> str:
    return TITLE_LABELS_RU.get(branch_id, f"ветвь {branch_id}")


def lambda_from_frequency_hz(freq_hz: float, radius: float) -> float:
    params = build_params(float(radius))
    return float(np.sqrt(max(float(freq_hz) / frequency_scale(params), 0.0)))


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
        metadata = list(tracking["metadata"])
        branch_lookup = {str(meta["branch_id"]): idx for idx, meta in enumerate(metadata)}
        if branch_id not in branch_lookup:
            raise KeyError(f"Tracked metadata does not contain {branch_id}.")
        branch_pos = branch_lookup[branch_id]

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


def build_shape_plot(
    branch_id: str,
    beta_deg: float,
    target_mus: Sequence[float],
    radii: Sequence[float],
    output_path: Path,
    title_label: str,
    show_spectral_index: bool,
    spectral_index_label: str,
    axis_label_style: str,
) -> Path:
    rows, tracking_cache = collect_branch_rows(
        branch_id=branch_id,
        beta_deg=beta_deg,
        target_mus=target_mus,
        radii=radii,
    )
    fig, axes = plt.subplots(len(radii), len(target_mus), figsize=(15.2, 7.7), squeeze=False)
    normalized_cases = []
    ell = build_params(float(radii[0])).L_total / 2.0

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
    x_axis_label = "координата x" if axis_label_style == "ru" else "x"
    y_axis_label = "координата y" if axis_label_style == "ru" else "y"

    for item in normalized_cases:
        row = item["row"]
        radius = float(row["r"])
        mu_value = float(row["mu"])
        row_pos = tuple(radii).index(radius)
        col_pos = tuple(target_mus).index(mu_value)
        ax = axes[row_pos, col_pos]
        joint = fem.N_ELEM
        x_def = item["x_base"] + MODE_SHAPE_SCALE * item["u"]
        y_def = item["y_base"] + MODE_SHAPE_SCALE * item["v"]
        lambda_value = lambda_from_frequency_hz(freq_hz=float(row["frequency_hz"]), radius=radius)
        spectral_index = int(row["current_sorted_index"])

        ax.plot(item["x_base"][: joint + 1], item["y_base"][: joint + 1], color="0.78", linestyle="--", linewidth=1.1)
        ax.plot(item["x_base"][joint:], item["y_base"][joint:], color="0.78", linestyle="--", linewidth=1.1)
        ax.plot(x_def[: joint + 1], y_def[: joint + 1], color="#1f77b4", linewidth=2.2)
        ax.plot(x_def[joint:], y_def[joint:], color="#ff7f0e", linewidth=2.2)
        ax.scatter([item["x_base"][joint]], [item["y_base"][joint]], color="black", s=12, zorder=5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.grid(True, alpha=0.18)
        title = (
            f"β = {beta_deg:g}°, r = {radius * 1e3:.0f} мм, "
            f"μ = {mu_label(mu_value)}, Λ = {lambda_value:.3f}"
        )
        if show_spectral_index:
            title = f"{title}, {spectral_index_label}: {spectral_index}"
        ax.set_title(title, fontsize=9.6 if show_spectral_index else 10.2)
        if row_pos == len(radii) - 1:
            ax.set_xlabel(x_axis_label)
        if col_pos == 0:
            ax.set_ylabel(y_axis_label)

    handles = [
        Line2D([0], [0], color="0.78", linestyle="--", linewidth=1.1, label="недеформированная геометрия"),
        Line2D([0], [0], color="#1f77b4", linewidth=2.2, label="левое плечо"),
        Line2D([0], [0], color="#ff7f0e", linewidth=2.2, label="правое плечо"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965), fontsize=10)
    fig.suptitle(f"Формы колебаний: β = {beta_deg:g}°, {title_label}", y=0.992, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return output_path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a Russian-labeled mode-shape plot for a tracked bending descendant branch."
    )
    parser.add_argument("--branch-id", default=DEFAULT_BRANCH_ID, help="Tracked branch id, e.g. bending_desc_01.")
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA_DEG, help="Target beta angle in degrees.")
    parser.add_argument("--target-mus", nargs="+", type=float, default=DEFAULT_TARGET_MUS, help="Target mu values.")
    parser.add_argument(
        "--radii",
        nargs="+",
        type=float,
        default=tuple(float(value) for value in DEFAULT_TARGET_RADII),
        help="Target beam radii in meters.",
    )
    parser.add_argument("--output", default=None, help="Output PNG path. Relative paths are resolved from the repo root.")
    parser.add_argument("--title-label", default=None, help="Russian branch label for the figure title.")
    parser.add_argument("--no-spectral-index", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--spectral-index-label", default="место в спектре", help=argparse.SUPPRESS)
    parser.add_argument("--axis-label-style", choices=("ru", "xy"), default="ru", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    args.target_mus = tuple(float(value) for value in args.target_mus)
    args.radii = tuple(float(value) for value in args.radii)
    if not args.target_mus:
        parser.error("--target-mus must contain at least one value.")
    if not args.radii:
        parser.error("--radii must contain at least one value.")
    return args


def main(argv: Sequence[str] | None = None) -> Path:
    args = parse_args(argv)
    title_label = args.title_label if args.title_label is not None else infer_title_label(args.branch_id)
    output_path = resolve_output_path(
        value=args.output,
        branch_id=args.branch_id,
        beta_deg=args.beta,
        target_mus=args.target_mus,
    )
    out_path = build_shape_plot(
        branch_id=args.branch_id,
        beta_deg=args.beta,
        target_mus=args.target_mus,
        radii=args.radii,
        output_path=output_path,
        title_label=title_label,
        show_spectral_index=not args.no_spectral_index,
        spectral_index_label=args.spectral_index_label,
        axis_label_style=args.axis_label_style,
    )
    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    print(f"saved PNG: {out_path}")
    return out_path


if __name__ == "__main__":
    main()
