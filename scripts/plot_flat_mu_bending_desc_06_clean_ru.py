from __future__ import annotations

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

from my_project.analytic.formulas import frequency_scale  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.analyze_flat_mu_bending_energy import (  # noqa: E402
    MODE_SHAPE_SCALE,
    RESULTS_DIR,
    TARGET_MUS,
    TARGET_RADII,
    build_node_coordinates,
    build_params,
    collect_case_rows,
    expand_full_vector,
    mu_key,
)
from scripts.sweep_grid_policy import ANALYSIS_BETA_STEP, ANALYSIS_MU_STEP  # noqa: E402


BETA_DEG = 15.0
BRANCH_ID = "bending_desc_06"
OUTPUT_PATH = RESULTS_DIR / "flat_mu_bending_shapes_beta15_bending_desc_06_clean_ru.png"


def mu_label(mu: float) -> str:
    rounded = round(float(mu))
    if abs(float(mu) - rounded) < 1e-9:
        return str(int(rounded))
    return f"{float(mu):.1f}".rstrip("0").rstrip(".")


def lambda_from_frequency_hz(freq_hz: float, radius: float) -> float:
    params = build_params(float(radius))
    return float(np.sqrt(max(float(freq_hz) / frequency_scale(params), 0.0)))


def build_clean_shape_plot() -> Path:
    selected_branches = [{"beta": BETA_DEG, "branch_id": BRANCH_ID}]
    case_rows, tracking_cache = collect_case_rows(selected_branches=selected_branches)
    rows = [
        row
        for row in case_rows
        if str(row["branch_id"]) == BRANCH_ID and abs(float(row["beta"]) - BETA_DEG) < 1e-9
    ]
    rows.sort(key=lambda row: (float(row["r"]), float(row["mu"])))
    if not rows:
        raise RuntimeError(f"No tracked rows found for beta={BETA_DEG:g}, branch={BRANCH_ID}.")

    fig, axes = plt.subplots(len(TARGET_RADII), len(TARGET_MUS), figsize=(15.2, 7.7), squeeze=False)
    normalized_cases = []
    ell = build_params(TARGET_RADII[0]).L_total / 2.0

    for row in rows:
        radius = float(row["r"])
        mu_value = float(row["mu"])
        tracking = tracking_cache[(radius, BETA_DEG)]
        metadata = list(tracking["metadata"])
        branch_lookup = {str(meta["branch_id"]): idx for idx, meta in enumerate(metadata)}
        branch_pos = branch_lookup[BRANCH_ID]
        snapshot = tracking["snapshots"][mu_key(mu_value)]
        full_vec = expand_full_vector(np.asarray(snapshot["vecs"][:, branch_pos], dtype=float))
        x_base, y_base = build_node_coordinates(mu=mu_value, beta_deg=BETA_DEG, ell=ell)
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
        lambda_value = lambda_from_frequency_hz(freq_hz=float(row["frequency_hz"]), radius=radius)

        ax.plot(item["x_base"][: joint + 1], item["y_base"][: joint + 1], color="0.78", linestyle="--", linewidth=1.1)
        ax.plot(item["x_base"][joint:], item["y_base"][joint:], color="0.78", linestyle="--", linewidth=1.1)
        ax.plot(x_def[: joint + 1], y_def[: joint + 1], color="#1f77b4", linewidth=2.2)
        ax.plot(x_def[joint:], y_def[joint:], color="#ff7f0e", linewidth=2.2)
        ax.scatter([item["x_base"][joint]], [item["y_base"][joint]], color="black", s=12, zorder=5)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(*x_limits)
        ax.set_ylim(*y_limits)
        ax.grid(True, alpha=0.18)
        ax.set_title(
            f"β = {BETA_DEG:g}°, r = {radius * 1e3:.0f} мм, μ = {mu_label(mu_value)}, Λ = {lambda_value:.3f}",
            fontsize=10.2,
        )
        if row_pos == len(TARGET_RADII) - 1:
            ax.set_xlabel("x")
        if col_pos == 0:
            ax.set_ylabel("y")

    handles = [
        Line2D([0], [0], color="0.78", linestyle="--", linewidth=1.1, label="недеформированная геометрия"),
        Line2D([0], [0], color="#1f77b4", linewidth=2.2, label="левое плечо"),
        Line2D([0], [0], color="#ff7f0e", linewidth=2.2, label="правое плечо"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.965), fontsize=10)
    fig.suptitle("Формы колебаний: β = 15°, потомок 6-й изгибной ветви", y=0.992, fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUTPUT_PATH, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return OUTPUT_PATH


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    out_path = build_clean_shape_plot()
    print(f"analysis beta step: {ANALYSIS_BETA_STEP:.1f} deg")
    print(f"analysis mu step: {ANALYSIS_MU_STEP:.4f}")
    print("local refinement windows: none")
    print(f"saved PNG: {out_path}")


if __name__ == "__main__":
    main()
