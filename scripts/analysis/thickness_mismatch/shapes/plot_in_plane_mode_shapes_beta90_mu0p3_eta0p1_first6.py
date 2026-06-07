from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
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

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    find_first_n_roots_eta,
    thickness_to_length_ratios,
)
from scripts.analysis import plot_mode_shapes_eta_beta_scan as shape_scan  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import reconstruct_components_eta  # noqa: E402


DEFAULT_BETA_DEG = 90.0
DEFAULT_EPSILON = 0.0025
DEFAULT_MU = 0.3
DEFAULT_ETA = 0.1
DEFAULT_MODE_INDICES = (1, 2, 3, 4, 5, 6)
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "in_plane_mode_shapes_beta90_mu0p3_eta0p1_first6"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "in_plane_mode_shapes_beta90_mu0p3_eta0p1_first6"

DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
RETRY_ROOT_SCAN_STEP = 0.005
RETRY_ROOT_LMAX0 = 45.0
DEFAULT_NUM_SHAPE_SAMPLES = 401
DEFAULT_L_TOTAL = 2.0
DEFAULT_MODE_SCALE = 0.14
DEFAULT_DPI = 240
SIGN_AMBIGUOUS_TOL = 1.0e-12

SUMMARY_FIELDS = [
    "sorted_index",
    "Lambda",
    "sign_flipped_for_visualization",
    "plot_sign_factor",
    "sign_normalization_basis",
    "short_rod_id",
    "sign_reference_value_before",
    "sign_reference_value_after",
    "normalization_factor",
    "shape_max_amplitude_raw",
    "shape_max_amplitude_normalized",
    "residual_norm",
    "smallest_singular_value",
    "singular_value_ratio",
    "thickness_ratio_1",
    "thickness_ratio_2",
    "warning",
    "notes",
]


@dataclass(frozen=True)
class Args:
    beta_deg: float
    epsilon: float
    mu: float
    eta: float
    mode_indices: tuple[int, ...]
    output_dir: Path
    num_shape_samples: int
    l_total: float
    mode_scale: float
    dpi: int
    smoke: bool


@dataclass(frozen=True)
class ModeShapeCase:
    sorted_index: int
    Lambda: float
    raw_coeff: np.ndarray
    plot_coeff: np.ndarray
    raw_components: dict[str, np.ndarray]
    signed_components: dict[str, np.ndarray]
    components: dict[str, np.ndarray]
    sign_flipped_for_visualization: bool
    plot_sign_factor: float
    sign_normalization_basis: str
    short_rod_id: int
    sign_reference_value_before: float
    sign_reference_value_after: float
    normalization_factor: float
    shape_max_amplitude_raw: float
    shape_max_amplitude_normalized: float
    residual_norm: float
    smallest_singular_value: float
    singular_value_ratio: float
    thickness_ratio_1: float
    thickness_ratio_2: float
    warning: str
    notes: str
    output_png: Path


def _fmt(value: object) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16g}"


def number_text(value: float | int | str | bool) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, str):
        return value
    value_f = float(value)
    return "" if not np.isfinite(value_f) else f"{value_f:.12g}"


def decimal_token(value: float) -> str:
    value_f = float(value)
    if abs(value_f) < 5.0e-15:
        return "0p0"
    prefix = "m" if value_f < 0.0 else ""
    text = f"{abs(value_f):.10g}"
    if "e" in text or "E" in text:
        text = f"{abs(value_f):.10f}".rstrip("0").rstrip(".")
    if "." not in text:
        text = f"{text}.0"
    return prefix + text.replace(".", "p")


def angle_token(value: float) -> str:
    value_f = float(value)
    prefix = "m" if value_f < 0.0 else ""
    text = f"{abs(value_f):.10g}"
    if text.endswith(".0"):
        text = text[:-2]
    return prefix + text.replace(".", "p")


def param_token(args: Args) -> str:
    return f"beta{angle_token(args.beta_deg)}_mu{decimal_token(args.mu)}_eta{decimal_token(args.eta)}"


def mode_range_token(args: Args) -> str:
    values = tuple(int(value) for value in args.mode_indices)
    if values == tuple(range(values[0], values[-1] + 1)) and len(values) > 1:
        return f"sorted{values[0]}_to_sorted{values[-1]}"
    return "_".join(f"sorted{value}" for value in values)


def individual_png_path(args: Args, sorted_index: int) -> Path:
    return Path(args.output_dir) / f"mode_shape_sorted{int(sorted_index)}_{param_token(args)}.png"


def grid_png_path(args: Args) -> Path:
    return Path(args.output_dir) / f"mode_shapes_{mode_range_token(args)}_{param_token(args)}_grid.png"


def summary_csv_path(args: Args) -> Path:
    return Path(args.output_dir) / f"mode_shapes_{param_token(args)}_first{len(args.mode_indices)}_summary.csv"


def report_path(args: Args) -> Path:
    return Path(args.output_dir) / f"mode_shapes_{param_token(args)}_first{len(args.mode_indices)}_report.md"


def resolve_repo_path(path: str | Path) -> Path:
    value = Path(path)
    if value.is_absolute():
        return value
    return REPO_ROOT / value


def validate_args(args: Args) -> None:
    if not isfinite(float(args.beta_deg)):
        raise ValueError("beta-deg must be finite.")
    if not isfinite(float(args.epsilon)) or float(args.epsilon) <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not (-1.0 < float(args.mu) < 1.0):
        raise ValueError("mu must lie inside (-1, 1).")
    if not (-1.0 < float(args.eta) < 1.0):
        raise ValueError("eta must lie inside (-1, 1).")
    if not args.mode_indices:
        raise ValueError("mode-indices must not be empty.")
    if any(int(value) <= 0 for value in args.mode_indices):
        raise ValueError("mode-indices must be positive one-based sorted indices.")
    if len(set(args.mode_indices)) != len(args.mode_indices):
        raise ValueError("mode-indices must be unique.")
    if int(args.num_shape_samples) < 5:
        raise ValueError("num-shape-samples must be at least 5.")
    if float(args.l_total) <= 0.0 or float(args.mode_scale) <= 0.0:
        raise ValueError("l-total and mode-scale must be positive.")
    if int(args.dpi) <= 0:
        raise ValueError("dpi must be positive.")


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot diagnostic in-plane Euler-Bernoulli thickness-mismatch mode shapes "
            "for fixed sorted roots at one parameter point. No descendant tracking is used."
        ),
    )
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--mode-indices", type=int, nargs="+", default=list(DEFAULT_MODE_INDICES))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--num-shape-samples", type=int, default=DEFAULT_NUM_SHAPE_SAMPLES)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--mode-scale", type=float, default=DEFAULT_MODE_SCALE)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--smoke", action="store_true", help="Generate only the first 3 requested sorted mode shapes.")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    mode_indices = tuple(int(value) for value in ns.mode_indices)
    output_dir = resolve_repo_path(ns.output_dir)
    num_shape_samples = int(ns.num_shape_samples)
    if bool(ns.smoke):
        mode_indices = tuple(mode_indices[: min(3, len(mode_indices))])
        num_shape_samples = min(num_shape_samples, 121)
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            output_dir = SMOKE_OUTPUT_DIR

    args = Args(
        beta_deg=float(ns.beta_deg),
        epsilon=float(ns.epsilon),
        mu=float(ns.mu),
        eta=float(ns.eta),
        mode_indices=mode_indices,
        output_dir=output_dir,
        num_shape_samples=num_shape_samples,
        l_total=float(ns.l_total),
        mode_scale=float(ns.mode_scale),
        dpi=int(ns.dpi),
        smoke=bool(ns.smoke),
    )
    validate_args(args)
    return args


def solve_sorted_roots(args: Args) -> tuple[np.ndarray, list[str]]:
    n_roots = max(int(value) for value in args.mode_indices)
    beta_rad = float(np.deg2rad(float(args.beta_deg)))
    roots = np.asarray(
        find_first_n_roots_eta(
            beta_rad,
            float(args.mu),
            float(args.epsilon),
            float(args.eta),
            int(n_roots),
            scan_step=DEFAULT_ROOT_SCAN_STEP,
            Lmax0=DEFAULT_ROOT_LMAX0,
        ),
        dtype=float,
    )
    warnings: list[str] = []
    if roots.size < n_roots:
        padded = np.full(n_roots, np.nan, dtype=float)
        padded[: roots.size] = roots
        roots = padded
    if np.count_nonzero(np.isfinite(roots[:n_roots])) < n_roots:
        retry = np.asarray(
            find_first_n_roots_eta(
                beta_rad,
                float(args.mu),
                float(args.epsilon),
                float(args.eta),
                int(n_roots),
                scan_step=RETRY_ROOT_SCAN_STEP,
                Lmax0=RETRY_ROOT_LMAX0,
            ),
            dtype=float,
        )
        if retry.size < n_roots:
            padded = np.full(n_roots, np.nan, dtype=float)
            padded[: retry.size] = retry
            retry = padded
        if np.count_nonzero(np.isfinite(retry[:n_roots])) > np.count_nonzero(np.isfinite(roots[:n_roots])):
            roots = retry
            warnings.append("retry_root_search_used_after_missing_default_roots")
    roots = np.sort(roots[:n_roots])
    missing_indices = [index + 1 for index, value in enumerate(roots) if not (isfinite(float(value)) and value > 0.0)]
    if missing_indices:
        warnings.append("missing_sorted_roots:" + ",".join(str(value) for value in missing_indices))
    return roots, warnings


def full_amplitude(components: dict[str, np.ndarray]) -> float:
    left = np.sqrt(np.asarray(components["u_left"], dtype=float) ** 2 + np.asarray(components["w_left"], dtype=float) ** 2)
    right = np.sqrt(np.asarray(components["u_right"], dtype=float) ** 2 + np.asarray(components["w_right"], dtype=float) ** 2)
    return max(float(np.nanmax(left)), float(np.nanmax(right)))


def short_rod_id(mu: float) -> int:
    left_length = 1.0 - float(mu)
    right_length = 1.0 + float(mu)
    return 1 if left_length <= right_length else 2


def global_vertical_components(
    components: dict[str, np.ndarray],
    *,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    beta_rad = float(np.deg2rad(float(beta_deg)))
    left_vertical = np.asarray(components["w_left"], dtype=float)
    right_vertical = (
        np.asarray(components["u_right"], dtype=float) * np.sin(beta_rad)
        + np.asarray(components["w_right"], dtype=float) * np.cos(beta_rad)
    )
    return left_vertical, right_vertical


def dominant_signed_value(values: np.ndarray) -> tuple[float, float]:
    finite_values = np.asarray(values, dtype=float)
    finite_values = finite_values[np.isfinite(finite_values)]
    if finite_values.size == 0:
        return float("nan"), float("nan")
    index = int(np.nanargmax(np.abs(finite_values)))
    value = float(finite_values[index])
    return value, abs(value)


def sign_factor_for_visualization(
    args: Args,
    components: dict[str, np.ndarray],
) -> tuple[float, str, int, float, float, str]:
    short_id = short_rod_id(float(args.mu))
    left_vertical, right_vertical = global_vertical_components(components, beta_deg=float(args.beta_deg))
    short_vertical = left_vertical if short_id == 1 else right_vertical
    ref_before, ref_abs = dominant_signed_value(short_vertical)
    basis = "short_rod_dominant_vertical_lobe"
    warning = ""
    if not np.isfinite(ref_before) or ref_abs <= SIGN_AMBIGUOUS_TOL:
        ref_before, ref_abs = dominant_signed_value(np.concatenate([left_vertical, right_vertical]))
        basis = "fallback_global_vertical_lobe"
    if not np.isfinite(ref_before) or ref_abs <= SIGN_AMBIGUOUS_TOL:
        ref_before, ref_abs = dominant_signed_value(
            np.concatenate(
                [
                    np.asarray(components["u_left"], dtype=float),
                    np.asarray(components["w_left"], dtype=float),
                    np.asarray(components["u_right"], dtype=float),
                    np.asarray(components["w_right"], dtype=float),
                ]
            )
        )
        basis = "fallback_global_component"
    if not np.isfinite(ref_before) or ref_abs <= SIGN_AMBIGUOUS_TOL:
        warning = "sign_reference_ambiguous"
        return 1.0, basis, short_id, ref_before, ref_before, warning
    sign_factor = -1.0 if ref_before < 0.0 else 1.0
    return sign_factor, basis, short_id, ref_before, sign_factor * ref_before, warning


def sign_components(components: dict[str, np.ndarray], sign_factor: float) -> dict[str, np.ndarray]:
    return {key: float(sign_factor) * np.asarray(value, dtype=float) for key, value in components.items()}


def normalize_signed_components(
    signed_components: dict[str, np.ndarray],
    normalization_factor: float,
) -> dict[str, np.ndarray]:
    scale = max(float(normalization_factor), SIGN_AMBIGUOUS_TOL)
    return {key: np.asarray(value, dtype=float) / scale for key, value in signed_components.items()}


def reconstruct_raw_cases(args: Args, roots: np.ndarray, s_norm: np.ndarray) -> list[dict[str, object]]:
    beta_rad = float(np.deg2rad(float(args.beta_deg)))
    raw_cases: list[dict[str, object]] = []
    ratio1, ratio2 = thickness_to_length_ratios(float(args.epsilon), float(args.mu), float(args.eta))
    for sorted_index in args.mode_indices:
        Lambda = float(roots[int(sorted_index) - 1])
        matrix = assemble_clamped_coupled_matrix_eta(
            Lambda,
            beta_rad,
            float(args.mu),
            float(args.epsilon),
            float(args.eta),
        )
        coeff, smallest, ratio = analytic_null_vector(matrix)
        components = reconstruct_components_eta(
            Lambda,
            mu=float(args.mu),
            eta=float(args.eta),
            epsilon=float(args.epsilon),
            coeff=coeff,
            s_norm=s_norm,
        )
        sign_factor, basis, short_id, ref_before, ref_after, sign_warning = sign_factor_for_visualization(args, components)
        signed_components = sign_components(components, sign_factor)
        residual = np.asarray(matrix, dtype=float) @ np.asarray(coeff, dtype=float)
        raw_cases.append(
            {
                "sorted_index": int(sorted_index),
                "Lambda": Lambda,
                "raw_coeff": np.asarray(coeff, dtype=float),
                "plot_coeff": float(sign_factor) * np.asarray(coeff, dtype=float),
                "raw_components": {key: np.asarray(value, dtype=float) for key, value in components.items()},
                "signed_components": signed_components,
                "sign_factor": float(sign_factor),
                "sign_basis": basis,
                "short_rod_id": int(short_id),
                "ref_before": float(ref_before),
                "ref_after": float(ref_after),
                "residual_norm": float(np.linalg.norm(residual)),
                "smallest_singular_value": float(smallest),
                "singular_value_ratio": float(ratio),
                "thickness_ratio_1": float(ratio1),
                "thickness_ratio_2": float(ratio2),
                "shape_max_amplitude_raw": full_amplitude(signed_components),
                "sign_warning": sign_warning,
            }
        )
    return raw_cases


def build_cases(args: Args, roots: np.ndarray, root_warnings: Sequence[str]) -> list[ModeShapeCase]:
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)
    raw_cases = reconstruct_raw_cases(args, roots, s_norm)
    normalization_factor = max(
        max(float(item["shape_max_amplitude_raw"]) for item in raw_cases),
        SIGN_AMBIGUOUS_TOL,
    )
    cases: list[ModeShapeCase] = []
    for item in raw_cases:
        components = normalize_signed_components(dict(item["signed_components"]), normalization_factor)
        normalized_amplitude = full_amplitude(components)
        warning_parts = list(root_warnings)
        if not np.isfinite(float(item["Lambda"])):
            warning_parts.append("nonfinite_lambda")
        if float(item["residual_norm"]) > 1.0e-7:
            warning_parts.append("residual_norm_gt_1e-7")
        if str(item["sign_warning"]):
            warning_parts.append(str(item["sign_warning"]))
        if float(item["thickness_ratio_1"]) > 0.1 or float(item["thickness_ratio_2"]) > 0.1:
            warning_parts.append("thin_rod_ratio_gt_0p1")
        cases.append(
            ModeShapeCase(
                sorted_index=int(item["sorted_index"]),
                Lambda=float(item["Lambda"]),
                raw_coeff=np.asarray(item["raw_coeff"], dtype=float),
                plot_coeff=np.asarray(item["plot_coeff"], dtype=float),
                raw_components={key: np.asarray(value, dtype=float) for key, value in dict(item["raw_components"]).items()},
                signed_components={
                    key: np.asarray(value, dtype=float) for key, value in dict(item["signed_components"]).items()
                },
                components=components,
                sign_flipped_for_visualization=float(item["sign_factor"]) < 0.0,
                plot_sign_factor=float(item["sign_factor"]),
                sign_normalization_basis=str(item["sign_basis"]),
                short_rod_id=int(item["short_rod_id"]),
                sign_reference_value_before=float(item["ref_before"]),
                sign_reference_value_after=float(item["ref_after"]),
                normalization_factor=float(normalization_factor),
                shape_max_amplitude_raw=float(item["shape_max_amplitude_raw"]),
                shape_max_amplitude_normalized=float(normalized_amplitude),
                residual_norm=float(item["residual_norm"]),
                smallest_singular_value=float(item["smallest_singular_value"]),
                singular_value_ratio=float(item["singular_value_ratio"]),
                thickness_ratio_1=float(item["thickness_ratio_1"]),
                thickness_ratio_2=float(item["thickness_ratio_2"]),
                warning=";".join(dict.fromkeys(warning_parts)) if warning_parts else "no",
                notes=(
                    "fixed sorted in-plane mode at one parameter point; no descendant tracking; "
                    "plot_coeff and plotted components are sign-normalized for visualization only"
                ),
                output_png=individual_png_path(args, int(item["sorted_index"])),
            )
        )
    return cases


def plot_namespace(args: Args) -> argparse.Namespace:
    return argparse.Namespace(l_total=float(args.l_total), mu=float(args.mu), mode_scale=float(args.mode_scale))


def axis_limits_for_cases(
    args: Args,
    cases: Sequence[ModeShapeCase],
    s_norm: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    plot_args = plot_namespace(args)
    x_values: list[np.ndarray] = []
    y_values: list[np.ndarray] = []
    for case in cases:
        base = shape_scan.base_coordinates(plot_args, s_norm, beta_deg=float(args.beta_deg))
        deformed = shape_scan.deformed_coordinates(plot_args, case.components, s_norm, beta_deg=float(args.beta_deg))
        x_values.extend([base[0], base[2], deformed[0], deformed[2]])
        y_values.extend([base[1], base[3], deformed[1], deformed[3]])
    x_all = np.concatenate([np.asarray(value, dtype=float) for value in x_values])
    y_all = np.concatenate([np.asarray(value, dtype=float) for value in y_values])
    x_span = max(float(np.nanmax(x_all) - np.nanmin(x_all)), 1.0)
    y_span = max(float(np.nanmax(y_all) - np.nanmin(y_all)), 0.4)
    return (
        (float(np.nanmin(x_all) - 0.06 * x_span), float(np.nanmax(x_all) + 0.06 * x_span)),
        (float(np.nanmin(y_all) - 0.18 * y_span), float(np.nanmax(y_all) + 0.18 * y_span)),
    )


def draw_case(
    ax: plt.Axes,
    args: Args,
    case: ModeShapeCase,
    s_norm: np.ndarray,
    *,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
    include_labels: bool,
) -> None:
    plot_args = plot_namespace(args)
    x_left_base, y_left_base, x_right_base, y_right_base = shape_scan.base_coordinates(
        plot_args,
        s_norm,
        beta_deg=float(args.beta_deg),
    )
    x_left_def, y_left_def, x_right_def, y_right_def = shape_scan.deformed_coordinates(
        plot_args,
        case.components,
        s_norm,
        beta_deg=float(args.beta_deg),
    )
    ax.plot(x_left_base, y_left_base, color="0.70", linestyle="--", linewidth=1.0, label="undeformed" if include_labels else None)
    ax.plot(x_right_base, y_right_base, color="0.70", linestyle="--", linewidth=1.0)
    ax.plot(x_left_def, y_left_def, color="#006BA4", linewidth=2.0, label="deformed" if include_labels else None)
    ax.plot(x_right_def, y_right_def, color="#006BA4", linewidth=2.0)
    ax.scatter([x_left_base[-1]], [y_left_base[-1]], color="black", s=18, zorder=5, label="joint" if include_labels else None)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*axis_limits[0])
    ax.set_ylim(*axis_limits[1])
    ax.grid(True, color="0.88", linewidth=0.6)


def title_for_case(args: Args, case: ModeShapeCase) -> str:
    return (
        f"In-plane EB sorted mode {int(case.sorted_index)}; Lambda={case.Lambda:.8g}\n"
        f"beta={args.beta_deg:g} deg, epsilon={args.epsilon:g}, mu={args.mu:g}, eta={args.eta:g}"
    )


def plot_individuals(
    args: Args,
    cases: Sequence[ModeShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
) -> list[Path]:
    paths: list[Path] = []
    for case in cases:
        path = case.output_png
        path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        draw_case(ax, args, case, s_norm, axis_limits=axis_limits, include_labels=True)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_title(title_for_case(args, case), fontsize=10)
        ax.legend(loc="best", fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(path, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
    return paths


def legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color="0.70", linestyle="--", linewidth=1.0, label="undeformed"),
        Line2D([0], [0], color="#006BA4", linewidth=2.0, label="deformed"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]


def plot_grid(
    args: Args,
    cases: Sequence[ModeShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
) -> Path:
    path = grid_png_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    n_cases = len(cases)
    n_cols = min(3, n_cases)
    n_rows = int(np.ceil(n_cases / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.35 * n_cols, 2.75 * n_rows + 0.9), squeeze=False)
    for idx, ax in enumerate(axes.ravel()):
        if idx >= n_cases:
            ax.axis("off")
            continue
        case = cases[idx]
        draw_case(ax, args, case, s_norm, axis_limits=axis_limits, include_labels=False)
        ax.set_title(f"sorted {case.sorted_index}, Lambda={case.Lambda:.6g}", fontsize=9)
        ax.set_xlabel("x", fontsize=8)
        ax.set_ylabel("y", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.legend(handles=legend_handles(), loc="lower center", ncol=3, frameon=False, fontsize=9)
    fig.suptitle(
        (
            f"Diagnostic fixed sorted in-plane mode shapes: beta={args.beta_deg:g} deg, "
            f"epsilon={args.epsilon:g}, mu={args.mu:g}, eta={args.eta:g}"
        ),
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.08, 1.0, 0.94), h_pad=0.9, w_pad=0.8)
    fig.savefig(path, dpi=int(args.dpi), bbox_inches="tight")
    plt.close(fig)
    return path


def write_summary_csv(args: Args, cases: Sequence[ModeShapeCase]) -> Path:
    path = summary_csv_path(args)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        for case in cases:
            writer.writerow(
                {
                    "sorted_index": int(case.sorted_index),
                    "Lambda": number_text(case.Lambda),
                    "sign_flipped_for_visualization": number_text(case.sign_flipped_for_visualization),
                    "plot_sign_factor": number_text(case.plot_sign_factor),
                    "sign_normalization_basis": case.sign_normalization_basis,
                    "short_rod_id": int(case.short_rod_id),
                    "sign_reference_value_before": number_text(case.sign_reference_value_before),
                    "sign_reference_value_after": number_text(case.sign_reference_value_after),
                    "normalization_factor": number_text(case.normalization_factor),
                    "shape_max_amplitude_raw": number_text(case.shape_max_amplitude_raw),
                    "shape_max_amplitude_normalized": number_text(case.shape_max_amplitude_normalized),
                    "residual_norm": number_text(case.residual_norm),
                    "smallest_singular_value": number_text(case.smallest_singular_value),
                    "singular_value_ratio": number_text(case.singular_value_ratio),
                    "thickness_ratio_1": number_text(case.thickness_ratio_1),
                    "thickness_ratio_2": number_text(case.thickness_ratio_2),
                    "warning": case.warning,
                    "notes": case.notes,
                }
            )
    return path


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_report(args: Args, cases: Sequence[ModeShapeCase], outputs: Sequence[Path]) -> Path:
    path = report_path(args)
    warning_rows = [case for case in cases if case.warning != "no"]
    lines = [
        "# In-Plane Sorted Mode Shapes",
        "",
        "This is a diagnostic-only analytic Euler-Bernoulli thickness-mismatch mode-shape plot set.",
        "",
        "## Fixed Parameters",
        "",
        f"- beta_deg: {_fmt(args.beta_deg)}",
        f"- epsilon: {_fmt(args.epsilon)}",
        f"- mu: {_fmt(args.mu)}",
        f"- eta: {_fmt(args.eta)}",
        f"- sorted mode indices: {', '.join(str(value) for value in args.mode_indices)}",
        f"- shape samples per rod: {int(args.num_shape_samples)}",
        "",
        "## First Sorted Frequencies",
        "",
    ]
    lines.extend(
        markdown_table(
            ["sorted index", "Lambda", "sign flipped", "warning"],
            [
                [
                    int(case.sorted_index),
                    f"{case.Lambda:.12g}",
                    "yes" if case.sign_flipped_for_visualization else "no",
                    case.warning,
                ]
                for case in cases
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Sign Normalization",
            "",
            "- The plotted eigenvector sign is chosen by the dominant global vertical displacement on the shorter rod.",
            "- If the shorter-rod reference is ambiguous, the fallback is the dominant global vertical lobe over both rods, then the dominant raw displacement component.",
            "- The sign choice is for visualization only and does not change frequencies or physical conclusions.",
            "",
            "## Shape Reconstruction Checks",
            "",
            f"- warning rows: {len(warning_rows)}",
            f"- maximum residual norm: {max(float(case.residual_norm) for case in cases):.12g}",
            f"- thickness ratios: rod 1 = {cases[0].thickness_ratio_1:.12g}, rod 2 = {cases[0].thickness_ratio_2:.12g}",
            "",
            "## Outputs",
            "",
        ]
    )
    for output in outputs:
        lines.append(f"- {display_path(output)}")
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "Sorted in-plane frequencies were used directly at this fixed parameter point. No descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article artifacts, determinant edits, old-solver edits, baseline-result edits, or analytic formula changes are part of this workflow.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def display_path(path: Path | None) -> str:
    if path is None:
        return "not generated"
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def run(args: Args) -> dict[str, object]:
    roots, root_warnings = solve_sorted_roots(args)
    requested_roots = roots[np.asarray([int(value) - 1 for value in args.mode_indices], dtype=int)]
    if np.any(~np.isfinite(requested_roots)):
        raise RuntimeError("Missing one or more requested sorted roots; see root warnings.")
    cases = build_cases(args, roots, root_warnings)
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)
    axis_limits = axis_limits_for_cases(args, cases, s_norm)
    individual_paths = plot_individuals(args, cases, s_norm, axis_limits)
    grid_path = plot_grid(args, cases, s_norm, axis_limits)
    summary_path = write_summary_csv(args, cases)
    outputs = [*individual_paths, grid_path, summary_path]
    report = write_report(args, cases, [*outputs, report_path(args)])
    return {
        "roots": roots,
        "root_warnings": root_warnings,
        "cases": cases,
        "individual_paths": individual_paths,
        "grid_path": grid_path,
        "summary_path": summary_path,
        "report_path": report,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    result = run(args)
    cases: Sequence[ModeShapeCase] = result["cases"]
    print("diagnostic-only fixed sorted in-plane mode shapes")
    print(
        f"beta_deg={args.beta_deg:g}, epsilon={args.epsilon:g}, "
        f"mu={args.mu:g}, eta={args.eta:g}"
    )
    print("mode selection: sorted in-plane frequencies only; no descendant tracking")
    for case in cases:
        print(
            f"sorted {case.sorted_index}: Lambda={case.Lambda:.12g}; "
            f"sign_flipped={'yes' if case.sign_flipped_for_visualization else 'no'}; "
            f"basis={case.sign_normalization_basis}; warning={case.warning}"
        )
    print(f"individual PNG count: {len(result['individual_paths'])}")
    for path in result["individual_paths"]:
        print(f"  saved: {display_path(path)}")
    print(f"grid PNG: {display_path(result['grid_path'])}")
    print(f"summary CSV: {display_path(result['summary_path'])}")
    print(f"report: {display_path(result['report_path'])}")
    print("no FEM/Gmsh/CalculiX; no article/formula/baseline changes")
    return result


if __name__ == "__main__":
    main()
