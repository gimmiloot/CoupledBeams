from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_timoshenko_shape_construction as shape_audit,
)
from scripts.lib import in_plane_shape_geometry as DISPLAY  # noqa: E402


MODEL_EB = shape_audit.MODEL_EB
MODEL_TIMO = shape_audit.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)
MODEL_PREFIX = {MODEL_EB: "eb", MODEL_TIMO: "timo"}
MODEL_LABEL = {MODEL_EB: "Euler--Bernoulli", MODEL_TIMO: "Timoshenko"}
MODEL_COLOR = {MODEL_EB: "#1f77b4", MODEL_TIMO: "#d55e00"}

DEFAULT_EPSILON = 0.03
DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
DEFAULT_MU_VALUES = (0.0, 0.2, 0.4, 0.6)
DEFAULT_SORTED_INDICES = (4, 5, 6)
SMOKE_MU_VALUES = (0.0,)
SMOKE_SORTED_INDICES = (4, 5)
DEFAULT_N_POINTS = 801
SMOKE_N_POINTS = 301
DEFAULT_DEFORMATION_SCALE_FRACTION = 0.08
FALLBACK_DEFORMATION_SCALE_FRACTION = 0.05
MINIMUM_DEFORMATION_SCALE_FRACTION = 0.02
DEFAULT_OUTPUT_DIR = Path("results") / "eb_timo_clean_mode_shapes"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_timo_clean_mode_shapes"
SUMMARY_CSV_NAME = "clean_mode_shape_summary.csv"
REPORT_NAME = "clean_mode_shape_report.md"
JOINT_TOL = 1.0e-6

SUMMARY_FIELDS = [
    "model",
    "epsilon",
    "beta_deg",
    "eta",
    "mu",
    "sorted_index",
    "Lambda",
    "sign_factor",
    "common_scale",
    "fold_flag",
    "near_overlap_flag",
    "plotted_joint_gap",
    "figure_file",
    "notes",
]


@dataclass(frozen=True)
class ScaleCheck:
    passed: bool
    fold_flag: bool
    near_overlap_flag: bool
    point_order_ok: bool
    plotted_joint_gap: float
    max_abs_kinematic_gap: float
    transform_ok: bool
    notes: tuple[str, ...]


@dataclass(frozen=True)
class ShapeRecord:
    result: shape_audit.ModeResult
    common_scale: float
    fold_flag: bool
    near_overlap_flag: bool
    point_order_ok: bool
    plotted_joint_gap: float
    max_abs_kinematic_gap: float
    transform_ok: bool
    figure_file: Path
    notes: str


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def value_token(value: float, *, min_decimals: int, max_decimals: int) -> str:
    text = f"{float(value):.{max_decimals}f}".rstrip("0").rstrip(".")
    if "." not in text and min_decimals > 0:
        text += "." + "0" * min_decimals
    if "." in text:
        decimals = len(text.split(".", 1)[1])
        if decimals < min_decimals:
            text += "0" * (min_decimals - decimals)
    return text.replace("-", "m").replace(".", "p")


def mu_text(value: float) -> str:
    return f"{float(value):g}"


def eps_token(value: float) -> str:
    return value_token(float(value), min_decimals=2, max_decimals=4)


def beta_token(value: float) -> str:
    return value_token(float(value), min_decimals=0, max_decimals=2)


def eta_token(value: float) -> str:
    return value_token(float(value), min_decimals=0, max_decimals=3)


def mode_indices_token(mode_indices: Sequence[int]) -> str:
    values = tuple(int(value) for value in mode_indices)
    if not values:
        raise ValueError("at least one sorted mode index is required")
    if values == tuple(range(min(values), max(values) + 1)):
        return f"{min(values)}_{max(values)}"
    return "_".join(str(value) for value in values)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Clean diagnostic EB/Timoshenko full in-plane mode-shape grids "
            "from sorted analytic frequencies."
        ),
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--mu-values", type=float, nargs="+", default=None)
    parser.add_argument("--mode-indices", type=int, nargs="+", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument(
        "--deformation-scale-fraction",
        type=float,
        default=DEFAULT_DEFORMATION_SCALE_FRACTION,
    )
    parser.add_argument(
        "--fallback-deformation-scale-fraction",
        type=float,
        default=FALLBACK_DEFORMATION_SCALE_FRACTION,
    )
    parser.add_argument(
        "--minimum-deformation-scale-fraction",
        type=float,
        default=MINIMUM_DEFORMATION_SCALE_FRACTION,
    )
    parser.add_argument(
        "--clean-style",
        action="store_true",
        help="Use uncluttered scientific panel titles and one figure-level legend.",
    )
    parser.add_argument("--no-combined", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    if bool(args.smoke):
        if args.mode_indices is None:
            args.mode_indices = list(SMOKE_SORTED_INDICES)
        if args.mu_values is None:
            args.mu_values = list(SMOKE_MU_VALUES)
        if Path(args.output_dir) == DEFAULT_OUTPUT_DIR:
            args.output_dir = SMOKE_OUTPUT_DIR
        args.n_points = min(int(args.n_points), SMOKE_N_POINTS)
    else:
        if args.mode_indices is None:
            args.mode_indices = list(DEFAULT_SORTED_INDICES)
        if args.mu_values is None:
            args.mu_values = list(DEFAULT_MU_VALUES)

    args.output_dir = repo_path(Path(args.output_dir))
    args.mode_indices = tuple(int(value) for value in args.mode_indices)
    args.mu_values = tuple(float(value) for value in args.mu_values)

    if int(args.n_points) < 51:
        raise ValueError("--n-points must be at least 51")
    if abs(float(args.eta)) > 1.0e-14:
        raise ValueError("this clean EB/Timoshenko comparison path is restricted to eta=0")
    if any(int(value) < 1 for value in args.mode_indices):
        raise ValueError("--mode-indices must be positive sorted-frequency indices")
    if float(args.deformation_scale_fraction) <= 0.0:
        raise ValueError("--deformation-scale-fraction must be positive")
    if float(args.fallback_deformation_scale_fraction) <= 0.0:
        raise ValueError("--fallback-deformation-scale-fraction must be positive")
    if float(args.minimum_deformation_scale_fraction) <= 0.0:
        raise ValueError("--minimum-deformation-scale-fraction must be positive")
    if not bool(args.clean_style):
        print("clean-style was not requested explicitly; writing clean-style figures anyway")
    return args


def build_results(args: argparse.Namespace) -> list[shape_audit.ModeResult]:
    n_roots = max(int(value) for value in args.mode_indices)
    provider = shape_audit.RootProvider(float(args.beta_deg), float(args.eta), n_roots)
    results: list[shape_audit.ModeResult] = []
    for mu in args.mu_values:
        eb_roots, timo_roots, root_warnings = provider.roots(float(args.epsilon), float(mu))
        for sorted_index in args.mode_indices:
            root_index = int(sorted_index) - 1
            eb_lambda = float(eb_roots[root_index])
            timo_lambda = float(timo_roots[root_index])
            if not (isfinite(eb_lambda) and isfinite(timo_lambda)):
                raise RuntimeError(
                    f"non-finite root for epsilon={args.epsilon:g}, "
                    f"mu={mu:g}, sorted={sorted_index}"
                )
            results.append(
                shape_audit.eb_mode_result(
                    epsilon=float(args.epsilon),
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=eb_lambda,
                    n_points=int(args.n_points),
                    case_role="clean_grid",
                    root_source="sorted_global_solver",
                    root_warnings=root_warnings,
                )
            )
            results.append(
                shape_audit.timo_mode_result(
                    epsilon=float(args.epsilon),
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=timo_lambda,
                    n_points=int(args.n_points),
                    case_role="clean_grid",
                    root_source="sorted_global_solver",
                    root_warnings=root_warnings,
                )
            )
    return align_cross_model_signs(results)


def local_displacement_vector(result: shape_audit.ModeResult) -> np.ndarray:
    sign = float(result.sign_factor)
    u1 = sign * np.asarray(result.rod1.u, dtype=float)
    w1 = sign * np.asarray(result.rod1.w, dtype=float)
    u2 = sign * np.asarray(result.rod2.u, dtype=float)
    w2 = sign * np.asarray(result.rod2.w, dtype=float)
    return np.concatenate([u1, w1, u2, w2])


def sign_alignment_key(result: shape_audit.ModeResult) -> tuple[float, float, float, float, int]:
    return (
        round(float(result.epsilon), 12),
        round(float(result.beta_deg), 12),
        round(float(result.eta), 12),
        round(float(result.mu), 12),
        int(result.sorted_index),
    )


def align_cross_model_signs(
    results: Sequence[shape_audit.ModeResult],
) -> list[shape_audit.ModeResult]:
    by_key: dict[tuple[float, float, float, float, int], dict[str, shape_audit.ModeResult]] = {}
    for result in results:
        by_key.setdefault(sign_alignment_key(result), {})[result.model] = result

    replacements: dict[int, shape_audit.ModeResult] = {}
    for grouped in by_key.values():
        eb = grouped.get(MODEL_EB)
        timo = grouped.get(MODEL_TIMO)
        if eb is None or timo is None:
            continue
        inner = float(np.dot(local_displacement_vector(eb), local_displacement_vector(timo)))
        if inner < 0.0:
            replacements[id(timo)] = replace(timo, sign_factor=-float(timo.sign_factor))

    return [replacements.get(id(result), result) for result in results]


def point_order_ok(result: shape_audit.ModeResult) -> bool:
    for rod_id in (1, 2):
        local = shape_audit.plotted_local_coordinate(result, rod_id)
        diffs = np.diff(local)
        if not (np.all(diffs >= -1.0e-12) or np.all(diffs <= 1.0e-12)):
            return False
    return True


def corrected_timoshenko_basis_ok(beta_deg: float) -> bool:
    beta_rad = np.deg2rad(float(beta_deg))
    expected_t1 = np.array([1.0, 0.0], dtype=float)
    expected_n1 = np.array([0.0, -1.0], dtype=float)
    expected_t2 = np.array([float(np.cos(beta_rad)), float(np.sin(beta_rad))], dtype=float)
    expected_n2 = np.array([float(np.sin(beta_rad)), -float(np.cos(beta_rad))], dtype=float)
    basis1 = DISPLAY.rod1_display_basis()
    basis2 = DISPLAY.rod2_display_basis(float(beta_deg))
    checks = [
        np.allclose(basis1.tangent, expected_t1, rtol=0.0, atol=1.0e-14),
        np.allclose(basis1.normal, expected_n1, rtol=0.0, atol=1.0e-14),
        np.allclose(basis2.tangent, expected_t2, rtol=0.0, atol=1.0e-14),
        np.allclose(basis2.normal, expected_n2, rtol=0.0, atol=1.0e-14),
        abs(float(np.dot(basis1.tangent, basis1.normal))) <= 1.0e-14,
        abs(float(np.dot(basis2.tangent, basis2.normal))) <= 1.0e-14,
        abs(float(np.linalg.norm(basis1.tangent)) - 1.0) <= 1.0e-14,
        abs(float(np.linalg.norm(basis1.normal)) - 1.0) <= 1.0e-14,
        abs(float(np.linalg.norm(basis2.tangent)) - 1.0) <= 1.0e-14,
        abs(float(np.linalg.norm(basis2.normal)) - 1.0) <= 1.0e-14,
    ]
    return all(bool(value) for value in checks)


def plotted_joint_gap(result: shape_audit.ModeResult, scale_fraction: float) -> float:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(scale_fraction))
    xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
        result,
        transverse_only=False,
        scale=scale,
    )
    return float(np.linalg.norm(np.array([xd1[-1] - xd2[0], yd1[-1] - yd2[0]], dtype=float)))


def eb_kinematic_gap(result: shape_audit.ModeResult) -> float:
    sign = float(result.sign_factor)
    u1 = sign * float(result.rod1.u[-1])
    w1 = sign * float(result.rod1.w[-1])
    u2 = sign * float(result.rod2.u[-1])
    w2 = sign * float(result.rod2.w[-1])
    dx1, dy1 = DISPLAY.eb_rod1_local_displacement_to_display(
        np.array([u1], dtype=float),
        np.array([w1], dtype=float),
    )
    dx2, dy2 = DISPLAY.eb_rod2_local_displacement_to_display(
        np.array([u2], dtype=float),
        np.array([w2], dtype=float),
        beta_deg=result.beta_deg,
    )
    return max(abs(float(dx1[0] - dx2[0])), abs(float(dy1[0] - dy2[0])))


def max_abs_kinematic_gap(result: shape_audit.ModeResult) -> float:
    if result.model == MODEL_TIMO and result.joint_row is not None:
        return float(result.joint_row["max_abs_kinematic_gap"])
    if result.model == MODEL_EB:
        return eb_kinematic_gap(result)
    return float("nan")


def regularity_rows_for_scale(
    result: shape_audit.ModeResult,
    *,
    scale_fraction: float,
) -> list[dict[str, object]]:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(scale_fraction))
    return [
        shape_audit.regularity_row(
            result,
            rod_id=rod_id,
            scale_kind="fixed_fraction",
            deformation_scale_fraction=float(scale_fraction),
            deformation_scale_value=scale,
        )
        for rod_id in (1, 2)
    ]


def scale_check(result: shape_audit.ModeResult, *, scale_fraction: float) -> ScaleCheck:
    rows = regularity_rows_for_scale(result, scale_fraction=float(scale_fraction))
    fold_flag = any(
        (
            isfinite(float(row["min_abs_dr_plot_ds"]))
            and float(row["min_abs_dr_plot_ds"]) < shape_audit.CURVE_SPEED_TOL
        )
        or not bool(row["projected_x_is_monotone"])
        for row in rows
    )
    near_overlap_flag = any(int(row["near_overlap_pair_count"]) > 0 for row in rows)
    order_ok = point_order_ok(result)
    plot_gap = plotted_joint_gap(result, float(scale_fraction))
    kin_gap = max_abs_kinematic_gap(result)
    transform_ok = True
    if result.model == MODEL_TIMO:
        transform_ok = corrected_timoshenko_basis_ok(result.beta_deg)
    notes: list[str] = []
    if not order_ok:
        notes.append("point_order_warning")
    if fold_flag:
        notes.append("fold_warning")
    if near_overlap_flag:
        notes.append("near_overlap_warning")
    if plot_gap > JOINT_TOL:
        notes.append("plotted_joint_gap_warning")
    if isfinite(kin_gap) and kin_gap > JOINT_TOL:
        notes.append("kinematic_gap_warning")
    if not transform_ok:
        notes.append("display_basis_warning")
    passed = bool(
        not fold_flag
        and not near_overlap_flag
        and order_ok
        and transform_ok
        and plot_gap <= JOINT_TOL
        and (not isfinite(kin_gap) or kin_gap <= JOINT_TOL)
    )
    if passed:
        notes.append("regularity_checks_passed")
    return ScaleCheck(
        passed=passed,
        fold_flag=bool(fold_flag),
        near_overlap_flag=bool(near_overlap_flag),
        point_order_ok=bool(order_ok),
        plotted_joint_gap=float(plot_gap),
        max_abs_kinematic_gap=float(kin_gap),
        transform_ok=bool(transform_ok),
        notes=tuple(notes),
    )


def deformed_arrays(
    result: shape_audit.ModeResult,
    *,
    scale_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(scale_fraction))
    return shape_audit.deformed_coordinates_with_scale(result, transverse_only=False, scale=scale)


def axis_limits(records: Sequence[ShapeRecord]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for record in records:
        result = record.result
        x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        xd1, yd1, xd2, yd2 = deformed_arrays(result, scale_fraction=record.common_scale)
        xs.extend([x1, x2, xd1, xd2])
        ys.extend([y1, y2, yd1, yd2])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.08 * x_span), float(np.max(x_all) + 0.08 * x_span)),
        (float(np.min(y_all) - 0.18 * y_span), float(np.max(y_all) + 0.18 * y_span)),
    )


def draw_shape(
    ax: plt.Axes,
    record: ShapeRecord,
    *,
    limits: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    result = record.result
    x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
    xd1, yd1, xd2, yd2 = deformed_arrays(result, scale_fraction=record.common_scale)
    color = MODEL_COLOR[result.model]
    ax.plot(x1, y1, color="0.78", linestyle="--", linewidth=0.8)
    ax.plot(x2, y2, color="0.78", linestyle="--", linewidth=0.8)
    ax.plot(xd1, yd1, color=color, linewidth=1.65)
    ax.plot(xd2, yd2, color=color, linewidth=1.65)
    ax.scatter([x1[-1]], [y1[-1]], color="black", s=12, zorder=5)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.91", linewidth=0.45)
    ax.tick_params(labelsize=7, pad=1.5)
    ax.set_title(
        rf"$\mu$={mu_text(result.mu)}, k={result.sorted_index}, $\Lambda$={result.Lambda:.4f}",
        fontsize=8.6,
        pad=7.0,
    )


def model_grid_path(
    output_dir: Path,
    *,
    model: str,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mode_indices: Sequence[int],
) -> Path:
    return output_dir / (
        f"{MODEL_PREFIX[model]}_full_shapes_eps{eps_token(epsilon)}"
        f"_beta{beta_token(beta_deg)}_eta{eta_token(eta)}"
        f"_modes{mode_indices_token(mode_indices)}_clean.png"
    )


def combined_grid_path(
    output_dir: Path,
    *,
    epsilon: float,
    beta_deg: float,
    eta: float,
    mode_indices: Sequence[int],
) -> Path:
    return output_dir / (
        f"eb_vs_timo_full_shapes_eps{eps_token(epsilon)}"
        f"_beta{beta_token(beta_deg)}_eta{eta_token(eta)}"
        f"_modes{mode_indices_token(mode_indices)}_clean.png"
    )


def write_model_grid(
    output_dir: Path,
    records: Sequence[ShapeRecord],
    *,
    model: str,
    mu_values: Sequence[float],
    mode_indices: Sequence[int],
    epsilon: float,
    beta_deg: float,
    eta: float,
) -> Path:
    model_records = [record for record in records if record.result.model == model]
    by_key = {
        (int(record.result.sorted_index), round(float(record.result.mu), 12)): record
        for record in model_records
    }
    output = model_grid_path(
        output_dir,
        model=model,
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        eta=float(eta),
        mode_indices=mode_indices,
    )
    limits = axis_limits(model_records)
    fig, axes = plt.subplots(
        len(mode_indices),
        len(mu_values),
        figsize=(3.85 * len(mu_values), 2.75 * len(mode_indices) + 1.15),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(mode_indices):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            record = by_key[(int(sorted_index), round(float(mu), 12))]
            draw_shape(ax, record, limits=limits)
            if row_index == len(mode_indices) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=8)

    handles = [
        Line2D([0], [0], color="0.78", linestyle="--", linewidth=0.8, label="undeformed rods"),
        Line2D([0], [0], color=MODEL_COLOR[model], linewidth=1.65, label="full deformed centerline"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.4)
    fig.suptitle(
        f"{MODEL_LABEL[model]}: epsilon={epsilon:g}, beta={beta_deg:g} deg, eta={eta:g}",
        fontsize=12.0,
        y=0.992,
    )
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.935), h_pad=1.15, w_pad=0.85)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return output


def write_combined_grid(
    output_dir: Path,
    records: Sequence[ShapeRecord],
    *,
    mu_values: Sequence[float],
    mode_indices: Sequence[int],
    epsilon: float,
    beta_deg: float,
    eta: float,
) -> Path:
    by_key = {
        (record.result.model, int(record.result.sorted_index), round(float(record.result.mu), 12)): record
        for record in records
    }
    output = combined_grid_path(
        output_dir,
        epsilon=float(epsilon),
        beta_deg=float(beta_deg),
        eta=float(eta),
        mode_indices=mode_indices,
    )
    limits = axis_limits(records)
    row_specs = [(mode_index, model) for mode_index in mode_indices for model in MODELS]
    fig, axes = plt.subplots(
        len(row_specs),
        len(mu_values),
        figsize=(3.85 * len(mu_values), 2.35 * len(row_specs) + 1.25),
        squeeze=False,
    )
    for row_index, (sorted_index, model) in enumerate(row_specs):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            record = by_key[(model, int(sorted_index), round(float(mu), 12))]
            draw_shape(ax, record, limits=limits)
            if row_index == len(row_specs) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel(f"{MODEL_PREFIX[model].upper()}\nglobal y", fontsize=8)

    handles = [
        Line2D([0], [0], color="0.78", linestyle="--", linewidth=0.8, label="undeformed rods"),
        Line2D([0], [0], color=MODEL_COLOR[MODEL_EB], linewidth=1.65, label=MODEL_LABEL[MODEL_EB]),
        Line2D([0], [0], color=MODEL_COLOR[MODEL_TIMO], linewidth=1.65, label=MODEL_LABEL[MODEL_TIMO]),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8.4)
    fig.suptitle(
        (
            f"Euler--Bernoulli vs Timoshenko: epsilon={epsilon:g}, "
            f"beta={beta_deg:g} deg, eta={eta:g}"
        ),
        fontsize=12.0,
        y=0.992,
    )
    fig.tight_layout(rect=(0.0, 0.055, 1.0, 0.945), h_pad=1.0, w_pad=0.85)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return output


def build_records(
    results: Sequence[shape_audit.ModeResult],
    *,
    output_dir: Path,
    scale_candidates: Sequence[float],
    mode_indices: Sequence[int],
) -> list[ShapeRecord]:
    candidates = tuple(float(value) for value in scale_candidates)
    if not candidates:
        raise ValueError("at least one scale candidate is required")
    checks_by_scale = {
        scale: {id(result): scale_check(result, scale_fraction=scale) for result in results}
        for scale in candidates
    }
    common_scale = candidates[-1]
    for scale in candidates:
        if all(check.passed for check in checks_by_scale[scale].values()):
            common_scale = scale
            break

    final_checks = checks_by_scale[common_scale]
    records: list[ShapeRecord] = []
    for result in results:
        figure_file = model_grid_path(
            output_dir,
            model=result.model,
            epsilon=result.epsilon,
            beta_deg=result.beta_deg,
            eta=result.eta,
            mode_indices=mode_indices,
        )
        final = final_checks[id(result)]
        note_parts = [
            "sorted frequency index used directly",
            "determinant-consistent local fields",
            "sign aligned by EB/Timoshenko local displacement inner product",
            "point_order_ok=true" if final.point_order_ok else "point_order_ok=false",
            f"max_abs_kinematic_gap={final.max_abs_kinematic_gap:.6e}",
        ]
        if result.model == MODEL_TIMO:
            note_parts.append(
                "corrected Timoshenko display bases t1=(1,0), n1=(0,-1), "
                "t2=(cos beta,sin beta), n2=(sin beta,-cos beta)"
            )
            note_parts.append("display_basis_ok=true" if final.transform_ok else "display_basis_ok=false")
        if abs(common_scale - candidates[0]) <= 1.0e-12:
            note_parts.append(f"scale {candidates[0]:.2f} passed all requested checks")
        else:
            note_parts.append(
                "scale ladder "
                f"{' -> '.join(f'{value:.2f}' for value in candidates)} "
                f"selected common scale {common_scale:.2f}"
            )
        if final.notes:
            note_parts.append(",".join(final.notes))
        if not checks_by_scale[candidates[0]][id(result)].passed:
            note_parts.append("requested_scale_check_failed")
        if not final.passed:
            note_parts.append("final_scale_check_failed")
        records.append(
            ShapeRecord(
                result=result,
                common_scale=common_scale,
                fold_flag=final.fold_flag,
                near_overlap_flag=final.near_overlap_flag,
                point_order_ok=final.point_order_ok,
                plotted_joint_gap=final.plotted_joint_gap,
                max_abs_kinematic_gap=final.max_abs_kinematic_gap,
                transform_ok=final.transform_ok,
                figure_file=figure_file,
                notes="; ".join(note_parts),
            )
        )
    return records


def csv_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"true", "1", "yes"}


def summary_key(row: dict[str, object]) -> tuple[str, float, float, float, float, int]:
    return (
        str(row["model"]),
        round(float(row["epsilon"]), 12),
        round(float(row["beta_deg"]), 12),
        round(float(row["eta"]), 12),
        round(float(row["mu"]), 12),
        int(row["sorted_index"]),
    )


def summary_sort_key(row: dict[str, object]) -> tuple[float, int, int, float]:
    model_rank = 0 if str(row["model"]) == MODEL_EB else 1
    return (
        round(float(row["epsilon"]), 12),
        model_rank,
        int(row["sorted_index"]),
        round(float(row["mu"]), 12),
    )


def summary_row(record: ShapeRecord) -> dict[str, object]:
    result = record.result
    return {
        "model": result.model,
        "epsilon": result.epsilon,
        "beta_deg": result.beta_deg,
        "eta": result.eta,
        "mu": result.mu,
        "sorted_index": result.sorted_index,
        "Lambda": result.Lambda,
        "sign_factor": float(result.sign_factor),
        "common_scale": record.common_scale,
        "fold_flag": record.fold_flag,
        "near_overlap_flag": record.near_overlap_flag,
        "plotted_joint_gap": record.plotted_joint_gap,
        "figure_file": shape_audit.rel(record.figure_file),
        "notes": record.notes,
    }


def read_existing_summary(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != SUMMARY_FIELDS:
            return []
        return [dict(row) for row in reader]


def write_merged_summary(path: Path, current_rows: Sequence[dict[str, object]]) -> tuple[Path, list[dict[str, object]]]:
    merged: dict[tuple[str, float, float, float, float, int], dict[str, object]] = {}
    for row in read_existing_summary(path):
        merged[summary_key(row)] = row
    for row in current_rows:
        merged[summary_key(row)] = row
    rows = sorted(merged.values(), key=summary_sort_key)
    shape_audit.write_csv(path, rows, SUMMARY_FIELDS)
    return path, rows


def rows_by_epsilon(rows: Sequence[dict[str, object]]) -> dict[float, list[dict[str, object]]]:
    grouped: dict[float, list[dict[str, object]]] = {}
    for row in rows:
        grouped.setdefault(round(float(row["epsilon"]), 12), []).append(row)
    return dict(sorted(grouped.items()))


def fmt_float(value: object) -> str:
    return f"{float(value):g}"


def write_report(
    path: Path,
    *,
    summary_rows: Sequence[dict[str, object]],
    figure_paths: Sequence[Path],
) -> Path:
    grouped = rows_by_epsilon(summary_rows)
    fold_count = sum(1 for row in summary_rows if csv_bool(row["fold_flag"]))
    near_count = sum(1 for row in summary_rows if csv_bool(row["near_overlap_flag"]))
    max_plot_gap = max((float(row["plotted_joint_gap"]) for row in summary_rows), default=float("nan"))
    transform_ok = all(
        "display_basis_ok=false" not in str(row.get("notes", ""))
        for row in summary_rows
        if str(row["model"]) == MODEL_TIMO
    )
    point_order_ok_all = all("point_order_ok=false" not in str(row.get("notes", "")) for row in summary_rows)
    all_pass = bool(fold_count == 0 and near_count == 0 and max_plot_gap <= JOINT_TOL and point_order_ok_all and transform_ok)

    lines = [
        "# Clean EB/Timoshenko Mode Shape Figures",
        "",
        "## Parameters",
        "",
    ]
    for epsilon, rows in grouped.items():
        mu_values = sorted({round(float(row["mu"]), 12) for row in rows})
        mode_indices = sorted({int(row["sorted_index"]) for row in rows})
        beta_values = sorted({round(float(row["beta_deg"]), 12) for row in rows})
        eta_values = sorted({round(float(row["eta"]), 12) for row in rows})
        lines.append(
            (
                f"- epsilon `{epsilon:g}`: beta_deg `{', '.join(f'{value:g}' for value in beta_values)}`, "
                f"eta `{', '.join(f'{value:g}' for value in eta_values)}`, "
                f"mu `{', '.join(f'{value:g}' for value in mu_values)}`, "
                f"sorted modes `{', '.join(str(value) for value in mode_indices)}`"
            )
        )

    lines.extend(
        [
            "",
            "## Scale And Checks",
            "",
        ]
    )
    for epsilon, rows in grouped.items():
        scales = sorted({round(float(row["common_scale"]), 12) for row in rows})
        lines.append(f"- Common scale for epsilon `{epsilon:g}`: `{', '.join(f'{value:.2f}' for value in scales)}`")
    lines.extend(
        [
            f"- All final panels pass regularity checks: `{str(all_pass).lower()}`",
            f"- Final fold flags: `{fold_count}`",
            f"- Final near-overlap flags: `{near_count}`",
            f"- Maximum final plotted joint gap: `{max_plot_gap:.6e}`",
            f"- Point order check passed for all panels: `{str(point_order_ok_all).lower()}`",
            f"- Corrected Timoshenko display transform check passed: `{str(transform_ok).lower()}`",
            "- Timoshenko display bases used: `t1=(1,0)`, `n1=(0,-1)`, `t2=(cos beta,sin beta)`, `n2=(sin beta,-cos beta)`.",
            "",
            "The workflow uses sorted frequencies only and does not use descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article workspaces, article figures, old determinant code, old solvers, or baseline-result edits.",
            "",
            "No analytic formulas, determinant entries, root solvers, local mode fields, frequencies, or Timoshenko shear coefficient k' are modified by this plotting script.",
            "",
            "## Generated Files",
            "",
            f"- Summary CSV: `{shape_audit.rel(path.parent / SUMMARY_CSV_NAME)}`",
            f"- Report: `{shape_audit.rel(path)}`",
        ]
    )
    for figure in sorted(figure_paths, key=lambda item: shape_audit.rel(item)):
        lines.append(f"- `{shape_audit.rel(figure)}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def output_pngs(output_dir: Path) -> list[Path]:
    if not output_dir.exists():
        return []
    return sorted(output_dir.glob("*_clean.png"), key=lambda item: item.name)


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    results = build_results(args)
    records = build_records(
        results,
        output_dir=output_dir,
        scale_candidates=(
            float(args.deformation_scale_fraction),
            float(args.fallback_deformation_scale_fraction),
            float(args.minimum_deformation_scale_fraction),
        ),
        mode_indices=args.mode_indices,
    )

    figure_paths: list[Path] = []
    figure_paths.append(
        write_model_grid(
            output_dir,
            records,
            model=MODEL_EB,
            mu_values=args.mu_values,
            mode_indices=args.mode_indices,
            epsilon=float(args.epsilon),
            beta_deg=float(args.beta_deg),
            eta=float(args.eta),
        )
    )
    figure_paths.append(
        write_model_grid(
            output_dir,
            records,
            model=MODEL_TIMO,
            mu_values=args.mu_values,
            mode_indices=args.mode_indices,
            epsilon=float(args.epsilon),
            beta_deg=float(args.beta_deg),
            eta=float(args.eta),
        )
    )
    if not bool(args.no_combined):
        figure_paths.append(
            write_combined_grid(
                output_dir,
                records,
                mu_values=args.mu_values,
                mode_indices=args.mode_indices,
                epsilon=float(args.epsilon),
                beta_deg=float(args.beta_deg),
                eta=float(args.eta),
            )
        )

    current_summary_rows = [summary_row(record) for record in records]
    summary_csv, merged_rows = write_merged_summary(output_dir / SUMMARY_CSV_NAME, current_summary_rows)
    report = write_report(
        output_dir / REPORT_NAME,
        summary_rows=merged_rows,
        figure_paths=output_pngs(output_dir),
    )

    common_scales = sorted({float(record.common_scale) for record in records})
    final_fold = sum(1 for record in records if record.fold_flag)
    final_near = sum(1 for record in records if record.near_overlap_flag)
    max_plot_gap = max((float(record.plotted_joint_gap) for record in records), default=float("nan"))
    print("generated clean EB/Timoshenko full mode-shape outputs:")
    for path in [summary_csv, report, *figure_paths]:
        print(f"  {shape_audit.rel(path)}")
    print(f"common scale for epsilon {float(args.epsilon):g}: {', '.join(f'{value:.2f}' for value in common_scales)}")
    print(f"final fold flags: {final_fold}")
    print(f"final near-overlap flags: {final_near}")
    print(f"max plotted joint gap: {max_plot_gap:.6e}")
    return {
        "report": report,
        "summary_csv": summary_csv,
        "figure_paths": figure_paths,
        "records": records,
        "summary_rows": merged_rows,
    }


if __name__ == "__main__":
    main()
