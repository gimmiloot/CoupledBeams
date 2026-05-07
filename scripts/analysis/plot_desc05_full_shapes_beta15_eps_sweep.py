from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


# ============================================================
# USER PARAMETERS
# Change these values for ordinary runs.
# CLI arguments, if provided, override these defaults.
# ============================================================

DEFAULT_BRANCH_ID = "bending_desc_05"
DEFAULT_BRANCH_NUMBER = 5
DEFAULT_BETA = 45.0
DEFAULT_MUS = (0.6, 0.8)
DEFAULT_EPSILONS = (0.0025, 0.005, 0.01)
DEFAULT_PLOT_KIND = "full"
DEFAULT_MODE_SCALE = 0.12
DEFAULT_NORMALIZE = "max-full"
DEFAULT_DPI = 240
DEFAULT_OUTPUT_DIR = "results"
DEFAULT_ALLOW_LOW_MAC = False
DEFAULT_BETA_STEPS = 200
DEFAULT_MU_STEPS = 260
DEFAULT_MAX_REFINEMENT_DEPTH = 8
DEFAULT_MIN_BETA_STEP = 1e-3
DEFAULT_MIN_MU_STEP = 1e-4


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    BranchPoint,
    DEFAULT_MAC_WARNING_THRESHOLD,
    DEFAULT_N_SOLVE,
    DEFAULT_N_TRACK,
    base_sorted_index_from_branch_id,
    branch_id_from_base_sorted_index,
    dense_mu_values_for_targets,
    nearby_sorted_lambdas,
    track_mu_sweep,
)
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    analytic_arm_energy_diagnostics,
    endpoint_consistency_diagnostics,
)
from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402


DEFAULT_L_TOTAL = 2.0
DEFAULT_SHAPE_METRIC = "full"
NUM_SAMPLES = 401
CONSISTENCY_TOL = 1e-7
CLAMP_TOL = 1e-9
FIGSIZE = (8.6, 4.8)
SUMMARY_FILENAME = "analytic_full_shapes_desc05_beta15_eps_sweep_summary.csv"
KNOWN_DIAGNOSTIC_BETA = 15.0
KNOWN_DIAGNOSTIC_MU = 0.8
KNOWN_DIAGNOSTIC_EPSILON = 0.0025
REPORTED_OLD_FULL_SHAPE_INDEX = 7

SUMMARY_FIELDNAMES = [
    "branch_id",
    "branch_number",
    "epsilon",
    "mu",
    "beta",
    "base_sorted_index",
    "current_sorted_index",
    "lambda",
    "mac_to_previous",
    "min_tracking_mac",
    "mean_tracking_mac",
    "relative_lambda_jump",
    "relative_gap_lower",
    "relative_gap_upper",
    "analytic_axial_energy_fraction",
    "analytic_right_axial_share",
    "warning_flag",
    "output_png",
]


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def resolve_repo_path(value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot analytic full mode shapes for a branch identity defined at beta=0, mu=0. "
            "Branch selection is delegated to scripts/lib/analytic_branch_tracking.py; no FEM Lambda is used."
        ),
    )
    parser.add_argument("--branch-id", default=DEFAULT_BRANCH_ID)
    parser.add_argument("--branch-number", type=int, default=DEFAULT_BRANCH_NUMBER)
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA)
    parser.add_argument("--mus", type=float, nargs="+", default=list(DEFAULT_MUS))
    parser.add_argument("--epsilons", type=float, nargs="+", default=list(DEFAULT_EPSILONS))
    parser.add_argument("--plot-kind", default=DEFAULT_PLOT_KIND, choices=("full",))
    parser.add_argument("--mode-scale", type=float, default=DEFAULT_MODE_SCALE)
    parser.add_argument("--normalize", default=DEFAULT_NORMALIZE, choices=("max-full",))
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--n-track", type=int, default=DEFAULT_N_TRACK)
    parser.add_argument("--n-solve", type=int, default=DEFAULT_N_SOLVE)
    parser.add_argument(
        "--beta-steps",
        type=int,
        default=DEFAULT_BETA_STEPS,
        help="Number of continuation samples from beta=0 to the target beta.",
    )
    parser.add_argument(
        "--mu-steps",
        type=int,
        default=DEFAULT_MU_STEPS,
        help="Number of dense mu continuation samples; increase --mu-steps if low-MAC tracking appears at larger beta.",
    )
    parser.add_argument(
        "--max-refinement-depth",
        type=int,
        default=DEFAULT_MAX_REFINEMENT_DEPTH,
        help="Maximum adaptive bisection depth for a low-MAC tracking step.",
    )
    parser.add_argument(
        "--min-beta-step",
        type=float,
        default=DEFAULT_MIN_BETA_STEP,
        help="Minimum beta step, in degrees, allowed during adaptive tracking refinement.",
    )
    parser.add_argument(
        "--min-mu-step",
        type=float,
        default=DEFAULT_MIN_MU_STEP,
        help="Minimum mu step allowed during adaptive tracking refinement.",
    )
    parser.add_argument("--shape-metric", default=DEFAULT_SHAPE_METRIC, choices=("full", "transverse"))
    parser.add_argument("--save-tracking-debug", action="store_true", help="Write disposable tracking CSVs under results/debug/.")
    parser.add_argument(
        "--check-known-contradiction",
        action="store_true",
        help="Run the old beta=15, epsilon=0.0025, mu=0.8 diagnostic. Off by default.",
    )
    parser.add_argument(
        "--allow-low-mac",
        action="store_true",
        default=DEFAULT_ALLOW_LOW_MAC,
        help="Allow exploratory figures even if canonical tracking falls below the MAC warning threshold.",
    )
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    if args.branch_number <= 0:
        parser.error("--branch-number must be positive.")
    if args.beta < 0.0:
        parser.error("--beta must be non-negative.")
    if any(float(mu) < -1e-12 for mu in args.mus):
        parser.error("--mus must be non-negative for this beta-then-mu plotting script.")
    if not args.mus:
        parser.error("--mus must contain at least one value.")
    if not args.epsilons:
        parser.error("--epsilons must contain at least one value.")
    if any(epsilon <= 0.0 for epsilon in args.epsilons):
        parser.error("--epsilons must be positive.")
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if args.dpi <= 0:
        parser.error("--dpi must be positive.")
    if args.n_track < args.branch_number:
        parser.error("--n-track must be at least --branch-number.")
    if args.n_solve < args.n_track:
        parser.error("--n-solve must be at least --n-track.")
    if args.beta_steps < 2:
        parser.error("--beta-steps must be at least 2.")
    if args.mu_steps < 2:
        parser.error("--mu-steps must be at least 2.")
    if args.max_refinement_depth < 0:
        parser.error("--max-refinement-depth must be non-negative.")
    if args.min_beta_step <= 0.0:
        parser.error("--min-beta-step must be positive.")
    if args.min_mu_step <= 0.0:
        parser.error("--min-mu-step must be positive.")

    try:
        branch_number_from_id = base_sorted_index_from_branch_id(args.branch_id)
    except ValueError as exc:
        parser.error(str(exc))
    if branch_number_from_id != args.branch_number:
        parser.error("--branch-id suffix and --branch-number must describe the same base branch.")
    return args


def mu_label(mu: float) -> str:
    rounded = round(float(mu))
    if abs(float(mu) - rounded) < 1e-9:
        return str(int(rounded))
    return f"{float(mu):g}"


def output_png_path(output_dir: Path, *, beta_deg: float, mu_value: float) -> Path:
    return output_dir / (
        f"analytic_full_shapes_desc05_beta{filename_number_token(beta_deg)}"
        f"_mu{filename_number_token(mu_value)}_eps_sweep.png"
    )


def debug_csv_path(
    *,
    branch_id: str,
    beta_deg: float,
    epsilon: float,
    max_mu: float,
) -> Path:
    return REPO_ROOT / "results" / "debug" / (
        f"analytic_tracking_{branch_id}_beta{filename_number_token(beta_deg)}"
        f"_eps{filename_number_token(epsilon)}_mu{filename_number_token(max_mu)}.csv"
    )


def arm_lengths(*, mu_value: float, l_total: float) -> tuple[float, float]:
    ell = float(l_total) / 2.0
    return ell * (1.0 - float(mu_value)), ell * (1.0 + float(mu_value))


def base_coordinates(
    *,
    s_norm: np.ndarray,
    beta_rad: float,
    mu_value: float,
    l_total: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    left_length, right_length = arm_lengths(mu_value=mu_value, l_total=l_total)
    xi = np.asarray(s_norm, dtype=float)
    x_left = left_length * xi
    y_left = np.zeros_like(x_left)
    x_right = left_length + right_length * xi * float(np.cos(beta_rad))
    y_right = right_length * xi * float(np.sin(beta_rad))
    return x_left, y_left, x_right, y_right


def deformed_coordinates(
    *,
    components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    beta_rad: float,
    mu_value: float,
    l_total: float,
    mode_scale: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_left, y_left, x_right, y_right = base_coordinates(
        s_norm=s_norm,
        beta_rad=beta_rad,
        mu_value=mu_value,
        l_total=l_total,
    )
    tangent_right = np.array([float(np.cos(beta_rad)), float(np.sin(beta_rad))], dtype=float)
    normal_right = np.array([-float(np.sin(beta_rad)), float(np.cos(beta_rad))], dtype=float)
    u_left = np.asarray(components["u_left"], dtype=float)
    w_left = np.asarray(components["w_left"], dtype=float)
    u_right = np.asarray(components["u_right"], dtype=float)
    w_right = np.asarray(components["w_right"], dtype=float)

    x_left_def = x_left + float(mode_scale) * u_left
    y_left_def = y_left + float(mode_scale) * w_left
    x_right_def = x_right + float(mode_scale) * (u_right * tangent_right[0] + w_right * normal_right[0])
    y_right_def = y_right + float(mode_scale) * (u_right * tangent_right[1] + w_right * normal_right[1])
    return x_left_def, y_left_def, x_right_def, y_right_def


def align_components_to_reference(
    components: dict[str, np.ndarray],
    reference: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    vector = np.concatenate([np.asarray(components[key], dtype=float) for key in ("u_left", "w_left", "u_right", "w_right")])
    reference_vector = np.concatenate([np.asarray(reference[key], dtype=float) for key in ("u_left", "w_left", "u_right", "w_right")])
    if float(np.dot(vector, reference_vector)) >= 0.0:
        return components
    return {key: -np.asarray(value, dtype=float) for key, value in components.items()}


def target_consistency(point: BranchPoint) -> dict[str, float]:
    matrix = assemble_clamped_coupled_matrix(
        point.Lambda,
        float(np.deg2rad(point.beta)),
        point.mu,
        point.epsilon,
    )
    return endpoint_consistency_diagnostics(
        matrix,
        point.coeff,
        Lambda=point.Lambda,
        beta_rad=float(np.deg2rad(point.beta)),
        mu_value=point.mu,
        epsilon=point.epsilon,
    )


def consistency_warning(consistency: dict[str, float]) -> str:
    if float(consistency["max_matrix_residual"]) > CONSISTENCY_TOL:
        return "needs_review"
    if float(consistency["max_field_residual"]) > CONSISTENCY_TOL:
        return "needs_review"
    if float(consistency["matrix_field_residual_max_abs_difference"]) > CONSISTENCY_TOL:
        return "needs_review"
    if float(consistency["external_clamp_residual"]) > CLAMP_TOL:
        return "needs_review"
    return "ok"


def combine_warning_flags(*flags: str) -> str:
    return "needs_review" if any(flag != "ok" for flag in flags) else "ok"


def branch_mac_stats(points: Sequence[BranchPoint]) -> tuple[float, float]:
    finite = [point.mac_to_previous for point in points if point.step_type != "base" and np.isfinite(point.mac_to_previous)]
    if not finite:
        return np.nan, np.nan
    return float(min(finite)), float(np.mean(finite))


def branch_points_up_to_target(points: Sequence[BranchPoint], *, beta: float, target_mu: float) -> list[BranchPoint]:
    out: list[BranchPoint] = []
    for point in points:
        if point.step_type in {"base", "beta"}:
            out.append(point)
        elif abs(float(point.beta) - float(beta)) <= 1e-10 and float(point.mu) <= float(target_mu) + 1e-10:
            out.append(point)
    return out


def build_summary_row(
    *,
    output_png: Path,
    target_point: BranchPoint,
    branch_points: Sequence[BranchPoint],
) -> dict[str, float | int | str]:
    consistency = target_consistency(target_point)
    energy = analytic_arm_energy_diagnostics(
        target_point.components,
        mu_value=target_point.mu,
        epsilon=target_point.epsilon,
        s_norm=np.linspace(0.0, 1.0, NUM_SAMPLES),
    )
    min_mac, mean_mac = branch_mac_stats(branch_points)
    path_warning = "needs_review" if any(point.warning_flag != "ok" for point in branch_points) else "ok"
    return {
        "branch_id": target_point.branch_id,
        "branch_number": int(target_point.base_sorted_index),
        "epsilon": float(target_point.epsilon),
        "mu": float(target_point.mu),
        "beta": float(target_point.beta),
        "base_sorted_index": int(target_point.base_sorted_index),
        "current_sorted_index": int(target_point.current_sorted_index),
        "lambda": float(target_point.Lambda),
        "mac_to_previous": float(target_point.mac_to_previous),
        "min_tracking_mac": float(min_mac),
        "mean_tracking_mac": float(mean_mac),
        "relative_lambda_jump": float(target_point.relative_lambda_jump),
        "relative_gap_lower": float(target_point.relative_gap_lower),
        "relative_gap_upper": float(target_point.relative_gap_upper),
        "analytic_axial_energy_fraction": float(energy["axial_energy_fraction"]),
        "analytic_right_axial_share": float(energy["right_axial_share"]),
        "warning_flag": combine_warning_flags(path_warning, target_point.warning_flag, consistency_warning(consistency)),
        "output_png": str(output_png),
    }


def axis_limits_for_mu(
    *,
    all_components: Sequence[dict[str, np.ndarray]],
    s_norm: np.ndarray,
    beta_rad: float,
    mu_value: float,
    l_total: float,
    mode_scale: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    x_base_left, y_base_left, x_base_right, y_base_right = base_coordinates(
        s_norm=s_norm,
        beta_rad=beta_rad,
        mu_value=mu_value,
        l_total=l_total,
    )
    x_values = [x_base_left, x_base_right]
    y_values = [y_base_left, y_base_right]
    for components in all_components:
        x_left, y_left, x_right, y_right = deformed_coordinates(
            components=components,
            s_norm=s_norm,
            beta_rad=beta_rad,
            mu_value=mu_value,
            l_total=l_total,
            mode_scale=mode_scale,
        )
        x_values.extend([x_left, x_right])
        y_values.extend([y_left, y_right])
    x_all = np.concatenate([np.asarray(value, dtype=float) for value in x_values])
    y_all = np.concatenate([np.asarray(value, dtype=float) for value in y_values])
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.06 * x_span), float(np.max(x_all) + 0.06 * x_span)),
        (float(np.min(y_all) - 0.12 * y_span), float(np.max(y_all) + 0.12 * y_span)),
    )


def plot_mu_panel(
    *,
    output_png: Path,
    beta_deg: float,
    mu_value: float,
    epsilons: Sequence[float],
    current_indices: Sequence[int],
    components_by_epsilon: Sequence[dict[str, np.ndarray]],
    s_norm: np.ndarray,
    l_total: float,
    mode_scale: float,
    dpi: int,
) -> None:
    beta_rad = float(np.deg2rad(beta_deg))
    x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(
        s_norm=s_norm,
        beta_rad=beta_rad,
        mu_value=mu_value,
        l_total=l_total,
    )
    x_limits, y_limits = axis_limits_for_mu(
        all_components=components_by_epsilon,
        s_norm=s_norm,
        beta_rad=beta_rad,
        mu_value=mu_value,
        l_total=l_total,
        mode_scale=mode_scale,
    )
    colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.plot(x_left_base, y_left_base, color="0.68", linestyle="--", linewidth=1.1)
    ax.plot(x_right_base, y_right_base, color="0.68", linestyle="--", linewidth=1.1)

    legend_handles: list[Line2D] = [
        Line2D([0], [0], color="0.68", linestyle="--", linewidth=1.1, label="недеформированная геометрия")
    ]
    for idx, (epsilon, current_index, components) in enumerate(
        zip(epsilons, current_indices, components_by_epsilon)
    ):
        color = colors[idx % len(colors)]
        x_left, y_left, x_right, y_right = deformed_coordinates(
            components=components,
            s_norm=s_norm,
            beta_rad=beta_rad,
            mu_value=mu_value,
            l_total=l_total,
            mode_scale=mode_scale,
        )
        label = f"ε = {float(epsilon):g}, current index = {int(current_index)}"
        ax.plot(x_left, y_left, color=color, linewidth=2.1)
        ax.plot(x_right, y_right, color=color, linewidth=2.1)
        legend_handles.append(Line2D([0], [0], color=color, linewidth=2.1, label=label))

    ax.scatter([x_left_base[-1]], [y_left_base[-1]], color="black", s=14, zorder=5)
    ax.set_title(
        "Полные аналитические формы: потомок 5-й изгибной ветви\n"
        f"β = {float(beta_deg):g}°, μ = {mu_label(mu_value)}",
        fontsize=11,
    )
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.grid(True, alpha=0.20)
    ax.legend(handles=legend_handles, fontsize=9, loc="best")
    fig.tight_layout(pad=0.9)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def write_summary_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)


def print_summary_table(rows: list[dict[str, float | int | str]], *, beta: float) -> None:
    print(f"Desc05 beta={float(beta):g} analytic full-shape epsilon sweep")
    print("mu, epsilon, branch id, base index, current index, Lambda, min MAC, axial fraction, warning, output PNG")
    for row in rows:
        print(
            f"{float(row['mu']):g}, {float(row['epsilon']):g}, {row['branch_id']}, "
            f"{int(row['base_sorted_index'])}, {int(row['current_sorted_index'])}, "
            f"{float(row['lambda']):.10g}, {float(row['min_tracking_mac']):.6f}, "
            f"{float(row['analytic_axial_energy_fraction']):.3e}, {row['warning_flag']}, {row['output_png']}"
        )


def print_known_contradiction_diagnostic(
    *,
    branch_id: str,
    branch_number: int,
    result_by_epsilon: dict[float, object],
    current_beta: float,
) -> None:
    if branch_id != DEFAULT_BRANCH_ID or abs(float(current_beta) - KNOWN_DIAGNOSTIC_BETA) > 1e-10:
        print("Known contradiction diagnostic is only relevant for beta=15, epsilon=0.0025.")
        print("Known contradiction diagnostic skipped: current tracking result does not contain beta=15, mu=0.8.")
        return

    result = result_by_epsilon.get(KNOWN_DIAGNOSTIC_EPSILON)
    if result is None:
        print("Known contradiction diagnostic is only relevant for beta=15, epsilon=0.0025.")
        print("Known contradiction diagnostic skipped: current tracking result does not contain beta=15, mu=0.8.")
        return

    try:
        point = result.point_at(branch_id, beta=KNOWN_DIAGNOSTIC_BETA, mu=KNOWN_DIAGNOSTIC_MU)
    except KeyError:
        print("Known contradiction diagnostic skipped: current tracking result does not contain beta=15, mu=0.8.")
        return

    branch_points = branch_points_up_to_target(
        result.points_for_branch(branch_id),
        beta=KNOWN_DIAGNOSTIC_BETA,
        target_mu=KNOWN_DIAGNOSTIC_MU,
    )
    min_mac, mean_mac = branch_mac_stats(branch_points)
    nearby = nearby_sorted_lambdas(point, (5, 6, 7))
    nearby_text = ", ".join(f"{index}: {value:.10g}" for index, value in nearby.items())
    if int(point.current_sorted_index) == REPORTED_OLD_FULL_SHAPE_INDEX:
        cause = (
            "needs review: the canonical helper currently agrees with the reported old full-shape index"
        )
    else:
        cause = (
            "old current_sorted_index=7 is invalidated by the canonical helper; it came from a retired "
            "single-branch/coarse low-MAC tracking path rather than the regression-tested branch identity rule"
        )
    print("Known contradiction diagnostic")
    print(f"branch_id: {branch_id} (base_sorted_index={branch_number})")
    print(f"beta: {KNOWN_DIAGNOSTIC_BETA:g} deg, epsilon: {KNOWN_DIAGNOSTIC_EPSILON:g}, mu: {KNOWN_DIAGNOSTIC_MU:g}")
    print(f"canonical current_sorted_index: {int(point.current_sorted_index)}")
    print(f"canonical Lambda: {float(point.Lambda):.10g}")
    print(f"MAC to previous step: {float(point.mac_to_previous):.6f}")
    print(f"min/mean tracking MAC on canonical path: {float(min_mac):.6f} / {float(mean_mac):.6f}")
    print(f"nearby sorted roots: {nearby_text}")
    print(f"reported old full-shape index: {REPORTED_OLD_FULL_SHAPE_INDEX}")
    print(f"diagnosis: {cause}")


def main(argv: Sequence[str] | None = None) -> list[dict[str, float | int | str]]:
    args = parse_args(argv)
    output_dir = resolve_repo_path(args.output_dir)
    branch_id = branch_id_from_base_sorted_index(int(args.branch_number))
    s_norm = np.linspace(0.0, 1.0, NUM_SAMPLES)
    summary_rows: list[dict[str, float | int | str]] = []
    result_by_epsilon: dict[float, object] = {}
    debug_paths: list[Path] = []
    mu_tracking_values = dense_mu_values_for_targets(args.mus, mu_steps=int(args.mu_steps))
    if args.allow_low_mac:
        print(
            "WARNING: --allow-low-mac is enabled; this run is exploratory, "
            "and low-MAC cases remain marked needs_review in CSV diagnostics."
        )

    for epsilon in args.epsilons:
        result = track_mu_sweep(
            epsilon=float(epsilon),
            beta=float(args.beta),
            mu_values=mu_tracking_values,
            n_track=int(args.n_track),
            n_solve=int(args.n_solve),
            mac_warning_threshold=DEFAULT_MAC_WARNING_THRESHOLD,
            shape_metric=str(args.shape_metric),
            beta_steps=int(args.beta_steps),
            num_samples=NUM_SAMPLES,
            allow_low_mac=bool(args.allow_low_mac),
            required_branch_ids=[branch_id],
            max_refinement_depth=int(args.max_refinement_depth),
            min_beta_step=float(args.min_beta_step),
            min_mu_step=float(args.min_mu_step),
        )
        result_by_epsilon[float(epsilon)] = result
        if args.save_tracking_debug:
            path = debug_csv_path(
                branch_id=branch_id,
                beta_deg=float(args.beta),
                epsilon=float(epsilon),
                max_mu=float(np.max(mu_tracking_values)),
            )
            debug_paths.append(result.write_debug_csv(path, branch_id=branch_id))

    for mu_value in args.mus:
        output_png = output_png_path(output_dir, beta_deg=float(args.beta), mu_value=float(mu_value))
        components_for_mu: list[dict[str, np.ndarray]] = []
        current_indices_for_mu: list[int] = []
        reference_components: dict[str, np.ndarray] | None = None

        for epsilon in args.epsilons:
            result = result_by_epsilon[float(epsilon)]
            target_point = result.point_at(branch_id, beta=float(args.beta), mu=float(mu_value))
            branch_points = branch_points_up_to_target(
                result.points_for_branch(branch_id),
                beta=float(args.beta),
                target_mu=float(mu_value),
            )
            components_plot = target_point.components
            if reference_components is None:
                reference_components = components_plot
            else:
                components_plot = align_components_to_reference(components_plot, reference_components)
            components_for_mu.append(components_plot)
            current_indices_for_mu.append(int(target_point.current_sorted_index))
            summary_rows.append(
                build_summary_row(
                    output_png=output_png,
                    target_point=target_point,
                    branch_points=branch_points,
                )
            )

        plot_mu_panel(
            output_png=output_png,
            beta_deg=float(args.beta),
            mu_value=float(mu_value),
            epsilons=[float(value) for value in args.epsilons],
            current_indices=current_indices_for_mu,
            components_by_epsilon=components_for_mu,
            s_norm=s_norm,
            l_total=float(args.l_total),
            mode_scale=float(args.mode_scale),
            dpi=int(args.dpi),
        )

    summary_path = output_dir / SUMMARY_FILENAME
    write_summary_csv(summary_path, summary_rows)
    print_summary_table(summary_rows, beta=float(args.beta))
    print(f"saved summary CSV: {summary_path}")
    if debug_paths:
        for path in debug_paths:
            print(f"saved tracking debug CSV: {path}")
    else:
        print("tracking debug CSV: not saved (pass --save-tracking-debug to write disposable files under results/debug/)")
    if args.check_known_contradiction:
        print_known_contradiction_diagnostic(
            branch_id=branch_id,
            branch_number=int(args.branch_number),
            result_by_epsilon=result_by_epsilon,
            current_beta=float(args.beta),
        )
    return summary_rows


if __name__ == "__main__":
    main()
