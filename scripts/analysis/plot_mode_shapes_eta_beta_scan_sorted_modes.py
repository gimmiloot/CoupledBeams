from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import assemble_clamped_coupled_matrix_eta  # noqa: E402
from scripts.analysis import plot_mode_shapes_eta_beta_scan as shape_scan  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import mac_value, reconstruct_components_eta  # noqa: E402


DEFAULT_ETA = 0.1
DEFAULT_MU = 0.3
DEFAULT_EPSILON = 0.0025
DEFAULT_BETA_START = 0.0
DEFAULT_BETA_END = 15.0
DEFAULT_BETA_STEP = 1.0
DEFAULT_SORTED_MODES = (2, 3)
DEFAULT_NUM_SORTED_ROOTS = 8
DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
DEFAULT_NUM_SHAPE_SAMPLES = 401
DEFAULT_L_TOTAL = 2.0
DEFAULT_MODE_SCALE = 0.14
DEFAULT_DPI = 240
DEFAULT_GRID_COLUMNS = 4
DEFAULT_SIGN_NORMALIZATION = "short-rod-up"
SIGN_REFERENCE_AMBIGUOUS_TOL = 1e-12

COEFF_NAMES = ("A1", "B1", "A2", "B2", "P1", "P2")
COMPONENT_KEYS = ("u_left", "w_left", "u_right", "w_right")
SUMMARY_FIELDNAMES = [
    "beta_deg",
    "sorted_index",
    "Lambda",
    "sign_normalization",
    "plot_sign_factor",
    "short_rod_id",
    "short_rod_reference_value_before_sign",
    "short_rod_reference_value_after_sign",
    "coeff_A1",
    "coeff_B1",
    "coeff_A2",
    "coeff_B2",
    "coeff_P1",
    "coeff_P2",
    "raw_coeff_A1",
    "raw_coeff_B1",
    "raw_coeff_A2",
    "raw_coeff_B2",
    "raw_coeff_P1",
    "raw_coeff_P2",
    "plot_coeff_A1",
    "plot_coeff_B1",
    "plot_coeff_A2",
    "plot_coeff_B2",
    "plot_coeff_P1",
    "plot_coeff_P2",
    "normalization_factor",
    "shape_max_amplitude_raw",
    "shape_max_amplitude_normalized",
    "mac_to_previous_same_sorted",
    "mac_to_beta0_same_sorted",
    "mac_sorted2_to_sorted3_at_beta",
    "residual_norm",
    "smallest_singular_value",
    "singular_value_ratio",
    "root_solver_warning",
    "warning",
    "notes",
]
GAP_FIELDNAMES = [
    "beta_deg",
    "Lambda_sorted2",
    "Lambda_sorted3",
    "gap_23",
    "mac_sorted2_to_sorted3",
    "mac_sorted2_to_previous",
    "mac_sorted3_to_previous",
]


@dataclass(frozen=True)
class SortedShapeCase:
    beta_deg: float
    sorted_index: int
    Lambda: float
    raw_coeff: np.ndarray
    coeff: np.ndarray
    raw_components: dict[str, np.ndarray]
    components_raw: dict[str, np.ndarray]
    components: dict[str, np.ndarray]
    sign_normalization: str
    plot_sign_factor: float
    short_rod_id: int
    short_rod_reference_value_before_sign: float
    short_rod_reference_value_after_sign: float
    residual_norm: float
    smallest_singular_value: float
    singular_value_ratio: float
    shape_max_amplitude_raw: float
    shape_max_amplitude_normalized: float
    mac_to_previous_same_sorted: float
    mac_to_beta0_same_sorted: float
    mac_sorted2_to_sorted3_at_beta: float
    output_png: Path
    warning: str
    notes: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot diagnostic analytic Euler-Bernoulli thickness-mismatch mode shapes "
            "for fixed sorted indices at each beta. No descendant branch tracking is used "
            "for mode selection."
        ),
    )
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--beta-start", type=float, default=DEFAULT_BETA_START)
    parser.add_argument("--beta-end", type=float, default=DEFAULT_BETA_END)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--sorted-modes", type=int, nargs="+", default=list(DEFAULT_SORTED_MODES))
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--num-sorted-roots", type=int, default=DEFAULT_NUM_SORTED_ROOTS)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--root-lmax0", type=float, default=DEFAULT_ROOT_LMAX0)
    parser.add_argument("--num-shape-samples", type=int, default=DEFAULT_NUM_SHAPE_SAMPLES)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--mode-scale", type=float, default=DEFAULT_MODE_SCALE)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--grid-columns", type=int, default=DEFAULT_GRID_COLUMNS)
    parser.add_argument(
        "--sign-normalization",
        choices=("none", "short-rod-up"),
        default=DEFAULT_SIGN_NORMALIZATION,
        help=(
            "Visual eigenvector sign convention for plotting. 'short-rod-up' flips the whole "
            "mode shape if the dominant plotted vertical deformation on the shorter rod is negative."
        ),
    )
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    validate_args(args)
    if args.output_dir is None:
        args.output_dir = default_output_dir(args)
    else:
        args.output_dir = resolve_repo_path(args.output_dir)
    args.sorted_modes = tuple(int(value) for value in args.sorted_modes)
    return args


def validate_args(args: argparse.Namespace) -> None:
    if not (-1.0 < float(args.mu) < 1.0):
        raise ValueError("--mu must lie inside (-1, 1).")
    if not (-1.0 < float(args.eta) < 1.0):
        raise ValueError("--eta must lie inside (-1, 1).")
    if float(args.epsilon) <= 0.0:
        raise ValueError("--epsilon must be positive.")
    if float(args.beta_start) < -1e-12:
        raise ValueError("--beta-start must be non-negative.")
    if float(args.beta_end) < float(args.beta_start):
        raise ValueError("--beta-end must be greater than or equal to --beta-start.")
    if float(args.beta_step) <= 0.0:
        raise ValueError("--beta-step must be positive.")
    if not args.sorted_modes:
        raise ValueError("--sorted-modes must contain at least one mode index.")
    if any(int(value) <= 0 for value in args.sorted_modes):
        raise ValueError("--sorted-modes values must be positive.")
    if int(args.num_sorted_roots) < max(int(value) for value in args.sorted_modes):
        raise ValueError("--num-sorted-roots must include the requested sorted modes.")
    if int(args.num_shape_samples) < 5:
        raise ValueError("--num-shape-samples must be at least 5.")
    if float(args.root_scan_step) <= 0.0 or float(args.root_lmax0) <= 0.0:
        raise ValueError("root scan parameters must be positive.")
    if float(args.l_total) <= 0.0 or float(args.mode_scale) <= 0.0:
        raise ValueError("--l-total and --mode-scale must be positive.")
    if int(args.dpi) <= 0 or int(args.grid_columns) <= 0:
        raise ValueError("--dpi and --grid-columns must be positive.")


def resolve_repo_path(path: str | Path) -> Path:
    value = Path(path)
    if value.is_absolute():
        return value
    return REPO_ROOT / value


def number_token(value: float) -> str:
    return shape_scan.number_token(float(value))


def number_text(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    value_f = float(value)
    return "" if not np.isfinite(value_f) else f"{value_f:.12g}"


def default_output_dir(args: argparse.Namespace) -> Path:
    sorted_part = "_".join(f"sorted{int(value)}" for value in args.sorted_modes)
    return REPO_ROOT / "results" / (
        f"mode_shapes_{sorted_part}_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}"
    )


def file_stem(args: argparse.Namespace, sorted_index: int) -> str:
    return (
        f"mode_shape_sorted{int(sorted_index)}_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}"
    )


def individual_png_path(args: argparse.Namespace, *, sorted_index: int, beta_deg: float) -> Path:
    return Path(args.output_dir) / f"{file_stem(args, int(sorted_index))}_beta_{number_token(float(beta_deg))}deg.png"


def grid_png_path(args: argparse.Namespace, *, sorted_index: int, beta_values: Sequence[float]) -> Path:
    return Path(args.output_dir) / (
        f"mode_shapes_sorted{int(sorted_index)}_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}_"
        f"beta_{number_token(float(beta_values[0]))}_to_{number_token(float(beta_values[-1]))}deg_grid.png"
    )


def comparison_png_path(args: argparse.Namespace, beta_values: Sequence[float]) -> Path:
    sorted_part = "_vs_".join(f"sorted{int(value)}" for value in args.sorted_modes)
    return Path(args.output_dir) / (
        f"mode_shapes_{sorted_part}_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}_"
        f"beta_{number_token(float(beta_values[0]))}_to_{number_token(float(beta_values[-1]))}deg_comparison_grid.png"
    )


def summary_csv_path(args: argparse.Namespace) -> Path:
    sorted_part = "_".join(f"sorted{int(value)}" for value in args.sorted_modes)
    return Path(args.output_dir) / (
        f"mode_shapes_{sorted_part}_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}_summary.csv"
    )


def gap_csv_path(args: argparse.Namespace, beta_values: Sequence[float]) -> Path:
    return Path(args.output_dir) / (
        f"sorted_gap_23_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}_"
        f"beta_{number_token(float(beta_values[0]))}_to_{number_token(float(beta_values[-1]))}deg.csv"
    )


def gap_png_path(args: argparse.Namespace, beta_values: Sequence[float]) -> Path:
    return Path(args.output_dir) / (
        f"sorted_gap_23_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}_"
        f"beta_{number_token(float(beta_values[0]))}_to_{number_token(float(beta_values[-1]))}deg.png"
    )


def root_solve_namespace(args: argparse.Namespace) -> argparse.Namespace:
    return argparse.Namespace(
        eta=float(args.eta),
        mu=float(args.mu),
        epsilon=float(args.epsilon),
        num_sorted_roots=int(args.num_sorted_roots),
        num_tracked_branches=max(int(value) for value in args.sorted_modes),
        root_scan_step=float(args.root_scan_step),
        root_lmax0=float(args.root_lmax0),
    )


def plot_namespace(args: argparse.Namespace) -> argparse.Namespace:
    return argparse.Namespace(
        l_total=float(args.l_total),
        mu=float(args.mu),
        mode_scale=float(args.mode_scale),
    )


def full_amplitude(components: dict[str, np.ndarray]) -> float:
    left = np.sqrt(np.asarray(components["u_left"], dtype=float) ** 2 + np.asarray(components["w_left"], dtype=float) ** 2)
    right = np.sqrt(np.asarray(components["u_right"], dtype=float) ** 2 + np.asarray(components["w_right"], dtype=float) ** 2)
    return max(float(np.max(left)), float(np.max(right)))


def component_vector(components: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate([np.asarray(components[key], dtype=float) for key in COMPONENT_KEYS])


def short_rod_id(mu: float) -> int:
    left_length_factor = 1.0 - float(mu)
    right_length_factor = 1.0 + float(mu)
    return 1 if left_length_factor <= right_length_factor else 2


def short_rod_global_vertical_deformation(
    components: dict[str, np.ndarray],
    *,
    beta_deg: float,
    short_rod: int,
) -> np.ndarray:
    if int(short_rod) == 1:
        return np.asarray(components["w_left"], dtype=float)
    beta_rad = float(np.deg2rad(float(beta_deg)))
    u_right = np.asarray(components["u_right"], dtype=float)
    w_right = np.asarray(components["w_right"], dtype=float)
    return u_right * np.sin(beta_rad) + w_right * np.cos(beta_rad)


def short_rod_reference_value(
    components: dict[str, np.ndarray],
    *,
    beta_deg: float,
    short_rod: int,
) -> tuple[float, float]:
    values = short_rod_global_vertical_deformation(
        components,
        beta_deg=float(beta_deg),
        short_rod=int(short_rod),
    )
    if values.size == 0:
        return float("nan"), float("nan")
    index = int(np.nanargmax(np.abs(values)))
    value = float(values[index])
    return value, float(abs(value))


def sign_factor_for_plot(
    args: argparse.Namespace,
    components: dict[str, np.ndarray],
    *,
    beta_deg: float,
    short_rod: int,
) -> tuple[float, float, float, str]:
    before_value, before_abs = short_rod_reference_value(
        components,
        beta_deg=float(beta_deg),
        short_rod=int(short_rod),
    )
    if str(args.sign_normalization) == "none":
        return 1.0, before_value, before_value, ""
    if not np.isfinite(before_value) or before_abs <= SIGN_REFERENCE_AMBIGUOUS_TOL:
        return 1.0, before_value, before_value, "short_rod_sign_reference_ambiguous"
    sign_factor = -1.0 if before_value < 0.0 else 1.0
    return sign_factor, before_value, sign_factor * before_value, ""


def sign_components(components: dict[str, np.ndarray], sign_factor: float) -> dict[str, np.ndarray]:
    return {key: float(sign_factor) * np.asarray(value, dtype=float) for key, value in components.items()}


def reconstruct_raw_case(
    args: argparse.Namespace,
    *,
    beta_deg: float,
    sorted_index: int,
    Lambda: float,
    s_norm: np.ndarray,
) -> tuple[np.ndarray, dict[str, np.ndarray], float, float, float]:
    matrix = assemble_clamped_coupled_matrix_eta(
        float(Lambda),
        float(np.deg2rad(float(beta_deg))),
        float(args.mu),
        float(args.epsilon),
        float(args.eta),
    )
    coeff, smallest, ratio = analytic_null_vector(matrix)
    components = reconstruct_components_eta(
        float(Lambda),
        mu=float(args.mu),
        eta=float(args.eta),
        epsilon=float(args.epsilon),
        coeff=coeff,
        s_norm=s_norm,
    )
    residual = np.asarray(matrix, dtype=float) @ np.asarray(coeff, dtype=float)
    return coeff, components, float(np.linalg.norm(residual)), float(smallest), float(ratio)


def solve_and_reconstruct(args: argparse.Namespace, beta_values: np.ndarray, sorted_roots: np.ndarray) -> list[SortedShapeCase]:
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)
    raw_by_mode: dict[int, list[dict[str, object]]] = {int(mode): [] for mode in args.sorted_modes}
    anchor_components_by_mode: dict[int, dict[str, np.ndarray]] = {}
    short_rod = short_rod_id(float(args.mu))

    for beta_idx, beta_deg in enumerate(beta_values):
        roots = np.asarray(sorted_roots[beta_idx], dtype=float)
        for sorted_index in args.sorted_modes:
            Lambda = float(roots[int(sorted_index) - 1])
            coeff, components, residual_norm, smallest, ratio = reconstruct_raw_case(
                args,
                beta_deg=float(beta_deg),
                sorted_index=int(sorted_index),
                Lambda=Lambda,
                s_norm=s_norm,
            )
            if int(sorted_index) not in anchor_components_by_mode:
                anchor_components_by_mode[int(sorted_index)] = components
            sign_factor, ref_before, ref_after, sign_warning = sign_factor_for_plot(
                args,
                components,
                beta_deg=float(beta_deg),
                short_rod=short_rod,
            )
            plot_coeff = float(sign_factor) * np.asarray(coeff, dtype=float)
            plot_components = sign_components(components, float(sign_factor))
            raw_by_mode[int(sorted_index)].append(
                {
                    "beta_deg": float(beta_deg),
                    "sorted_index": int(sorted_index),
                    "Lambda": Lambda,
                    "raw_coeff": np.asarray(coeff, dtype=float),
                    "coeff": np.asarray(plot_coeff, dtype=float),
                    "raw_components": {key: np.asarray(value, dtype=float) for key, value in components.items()},
                    "plot_sign_factor": float(sign_factor),
                    "short_rod_id": int(short_rod),
                    "short_rod_reference_value_before_sign": float(ref_before),
                    "short_rod_reference_value_after_sign": float(ref_after),
                    "sign_warning": sign_warning,
                    "components_raw": {key: np.asarray(value, dtype=float) for key, value in components.items()},
                    "components_plot_raw": {key: np.asarray(value, dtype=float) for key, value in plot_components.items()},
                    "residual_norm": residual_norm,
                    "smallest_singular_value": smallest,
                    "singular_value_ratio": ratio,
                    "shape_max_amplitude_raw": full_amplitude(plot_components),
                }
            )

    normalization = max(
        float(entry["shape_max_amplitude_raw"])
        for entries in raw_by_mode.values()
        for entry in entries
    )
    normalization = max(normalization, shape_scan.NEAR_ZERO_NORM if hasattr(shape_scan, "NEAR_ZERO_NORM") else 1e-12)

    cases: list[SortedShapeCase] = []
    by_beta_mode_raw = {
        (round(float(entry["beta_deg"]), 10), int(entry["sorted_index"])): entry
        for entries in raw_by_mode.values()
        for entry in entries
    }
    for sorted_index in args.sorted_modes:
        entries = raw_by_mode[int(sorted_index)]
        anchor_vector = component_vector(anchor_components_by_mode[int(sorted_index)])
        for idx, entry in enumerate(entries):
            beta_deg = float(entry["beta_deg"])
            raw_components = {key: np.asarray(value, dtype=float) for key, value in dict(entry["raw_components"]).items()}
            components_raw = {
                key: np.asarray(value, dtype=float)
                for key, value in dict(entry["components_plot_raw"]).items()
            }
            components = {key: value / normalization for key, value in components_raw.items()}
            current_vector = component_vector(raw_components)
            previous_mac = np.nan
            if idx > 0:
                previous_vector = component_vector(dict(entries[idx - 1]["raw_components"]))
                previous_mac = mac_value(previous_vector, current_vector)
            anchor_mac = mac_value(anchor_vector, current_vector)
            cross_mac = np.nan
            if set(args.sorted_modes) >= {2, 3}:
                other_index = 3 if int(sorted_index) == 2 else 2
                other_entry = by_beta_mode_raw[(round(beta_deg, 10), other_index)]
                cross_mac = mac_value(current_vector, component_vector(dict(other_entry["raw_components"])))
            normalized_amplitude = full_amplitude(components)
            warning_parts = []
            if not np.isfinite(float(entry["Lambda"])):
                warning_parts.append("nonfinite_lambda")
            if float(entry["residual_norm"]) > 1e-7:
                warning_parts.append("residual_norm_gt_1e-7")
            if str(entry["sign_warning"]):
                warning_parts.append(str(entry["sign_warning"]))
            if float(entry["short_rod_reference_value_after_sign"]) < -SIGN_REFERENCE_AMBIGUOUS_TOL:
                warning_parts.append("short_rod_reference_after_sign_negative")
            cases.append(
                SortedShapeCase(
                    beta_deg=beta_deg,
                    sorted_index=int(sorted_index),
                    Lambda=float(entry["Lambda"]),
                    raw_coeff=np.asarray(entry["raw_coeff"], dtype=float),
                    coeff=np.asarray(entry["coeff"], dtype=float),
                    raw_components=raw_components,
                    components_raw=components_raw,
                    components=components,
                    sign_normalization=str(args.sign_normalization),
                    plot_sign_factor=float(entry["plot_sign_factor"]),
                    short_rod_id=int(entry["short_rod_id"]),
                    short_rod_reference_value_before_sign=float(entry["short_rod_reference_value_before_sign"]),
                    short_rod_reference_value_after_sign=float(entry["short_rod_reference_value_after_sign"]),
                    residual_norm=float(entry["residual_norm"]),
                    smallest_singular_value=float(entry["smallest_singular_value"]),
                    singular_value_ratio=float(entry["singular_value_ratio"]),
                    shape_max_amplitude_raw=float(entry["shape_max_amplitude_raw"]),
                    shape_max_amplitude_normalized=normalized_amplitude,
                    mac_to_previous_same_sorted=float(previous_mac),
                    mac_to_beta0_same_sorted=float(anchor_mac),
                    mac_sorted2_to_sorted3_at_beta=float(cross_mac),
                    output_png=individual_png_path(args, sorted_index=int(sorted_index), beta_deg=beta_deg),
                    warning=";".join(warning_parts) if warning_parts else "no",
                    notes=(
                        "fixed sorted-index selection at this beta; no descendant branch tracking used for selection; "
                        "coeff_* and plot_coeff_* are sign-normalized for plotting; raw_coeff_* preserve the "
                        "determinant-nullspace orientation before plot sign normalization"
                    ),
                )
            )
    return cases


def cases_for_mode(cases: Sequence[SortedShapeCase], sorted_index: int) -> list[SortedShapeCase]:
    return [case for case in cases if int(case.sorted_index) == int(sorted_index)]


def axis_limits_for_cases(
    args: argparse.Namespace,
    cases: Sequence[SortedShapeCase],
    s_norm: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    plot_args = plot_namespace(args)
    x_values: list[np.ndarray] = []
    y_values: list[np.ndarray] = []
    for case in cases:
        base = shape_scan.base_coordinates(plot_args, s_norm, beta_deg=float(case.beta_deg))
        deformed = shape_scan.deformed_coordinates(plot_args, case.components, s_norm, beta_deg=float(case.beta_deg))
        x_values.extend([base[0], base[2], deformed[0], deformed[2]])
        y_values.extend([base[1], base[3], deformed[1], deformed[3]])
    x_all = np.concatenate([np.asarray(value, dtype=float) for value in x_values])
    y_all = np.concatenate([np.asarray(value, dtype=float) for value in y_values])
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.4)
    return (
        (float(np.min(x_all) - 0.06 * x_span), float(np.max(x_all) + 0.06 * x_span)),
        (float(np.min(y_all) - 0.18 * y_span), float(np.max(y_all) + 0.18 * y_span)),
    )


def draw_case(
    ax: plt.Axes,
    args: argparse.Namespace,
    case: SortedShapeCase,
    s_norm: np.ndarray,
    *,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
    labels: bool,
) -> None:
    plot_args = plot_namespace(args)
    x_left_base, y_left_base, x_right_base, y_right_base = shape_scan.base_coordinates(
        plot_args,
        s_norm,
        beta_deg=float(case.beta_deg),
    )
    x_left_def, y_left_def, x_right_def, y_right_def = shape_scan.deformed_coordinates(
        plot_args,
        case.components,
        s_norm,
        beta_deg=float(case.beta_deg),
    )
    ax.plot(x_left_base, y_left_base, color="0.70", linestyle="--", linewidth=1.0, label="undeformed" if labels else None)
    ax.plot(x_right_base, y_right_base, color="0.70", linestyle="--", linewidth=1.0)
    ax.plot(x_left_def, y_left_def, color="#006BA4", linewidth=2.0, label="deformed sorted mode" if labels else None)
    ax.plot(x_right_def, y_right_def, color="#006BA4", linewidth=2.0)
    ax.scatter([x_left_base[-1]], [y_left_base[-1]], color="black", s=18, zorder=5, label="joint" if labels else None)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*axis_limits[0])
    ax.set_ylim(*axis_limits[1])
    ax.grid(True, color="0.88", linewidth=0.6)


def sign_label(args: argparse.Namespace) -> str:
    if str(args.sign_normalization) == "short-rod-up":
        return "sign: short-rod-up"
    return "sign: none"


def single_title(args: argparse.Namespace, case: SortedShapeCase) -> str:
    return (
        f"Euler-Bernoulli mode shape, sorted {int(case.sorted_index)}, "
        f"eta={float(args.eta):g}, mu={float(args.mu):g}, epsilon={float(args.epsilon):g}, "
        f"beta={case.beta_deg:g} deg, Lambda={case.Lambda:.8g}; {sign_label(args)}"
    )


def plot_individuals(
    args: argparse.Namespace,
    cases: Sequence[SortedShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
) -> list[Path]:
    paths: list[Path] = []
    for case in cases:
        output = case.output_png
        output.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        draw_case(ax, args, case, s_norm, axis_limits=axis_limits, labels=True)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(loc="best", fontsize=8, frameon=False)
        ax.set_title(single_title(args, case), fontsize=10)
        fig.tight_layout()
        fig.savefig(output, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig)
        paths.append(output)
    return paths


def legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color="0.70", linestyle="--", linewidth=1.0, label="undeformed"),
        Line2D([0], [0], color="#006BA4", linewidth=2.0, label="deformed sorted mode"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]


def plot_mode_grid(
    args: argparse.Namespace,
    cases: Sequence[SortedShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
    *,
    sorted_index: int,
    beta_values: Sequence[float],
) -> Path:
    mode_cases = cases_for_mode(cases, int(sorted_index))
    n_cases = len(mode_cases)
    n_cols = min(int(args.grid_columns), n_cases)
    n_rows = int(np.ceil(n_cases / n_cols))
    output = grid_png_path(args, sorted_index=int(sorted_index), beta_values=beta_values)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.8 * n_cols, 1.95 * n_rows + 0.6), squeeze=False)
    for idx, ax in enumerate(axes.ravel()):
        if idx >= n_cases:
            ax.axis("off")
            continue
        case = mode_cases[idx]
        draw_case(ax, args, case, s_norm, axis_limits=axis_limits, labels=False)
        ax.set_title(f"beta={case.beta_deg:g} deg, Lambda={case.Lambda:.5g}", fontsize=9)
        ax.set_xlabel("x", fontsize=8)
        ax.set_ylabel("y", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.suptitle(
        (
            f"Diagnostic fixed sorted {int(sorted_index)} mode shapes: eta={float(args.eta):g}, "
            f"mu={float(args.mu):g}, epsilon={float(args.epsilon):g}; {sign_label(args)}"
        ),
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93), h_pad=0.9, w_pad=0.6)
    fig.savefig(output, dpi=int(args.dpi), bbox_inches="tight")
    plt.close(fig)
    return output


def plot_comparison_grid(
    args: argparse.Namespace,
    cases: Sequence[SortedShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
    *,
    beta_values: Sequence[float],
) -> Path | None:
    if tuple(int(value) for value in args.sorted_modes) != (2, 3):
        return None
    output = comparison_png_path(args, beta_values)
    output.parent.mkdir(parents=True, exist_ok=True)
    by_key = {(round(float(case.beta_deg), 10), int(case.sorted_index)): case for case in cases}
    fig, axes = plt.subplots(len(beta_values), 2, figsize=(7.8, 1.55 * len(beta_values) + 0.8), squeeze=False)
    for row, beta in enumerate(beta_values):
        for col, sorted_index in enumerate((2, 3)):
            ax = axes[row, col]
            case = by_key[(round(float(beta), 10), sorted_index)]
            draw_case(ax, args, case, s_norm, axis_limits=axis_limits, labels=False)
            ax.set_title(f"beta={float(beta):g} deg, sorted {sorted_index}, Lambda={case.Lambda:.5g}", fontsize=8.5)
            ax.tick_params(labelsize=7)
            if row == len(beta_values) - 1:
                ax.set_xlabel("x", fontsize=8)
            if col == 0:
                ax.set_ylabel("y", fontsize=8)
    fig.suptitle(
        (
            f"Fixed sorted mode comparison: eta={float(args.eta):g}, "
            f"mu={float(args.mu):g}, epsilon={float(args.epsilon):g}; {sign_label(args)}"
        ),
        fontsize=13,
        y=0.997,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.985), h_pad=0.65, w_pad=0.5)
    fig.savefig(output, dpi=int(args.dpi), bbox_inches="tight")
    plt.close(fig)
    return output


def write_summary_csv(args: argparse.Namespace, cases: Sequence[SortedShapeCase]) -> Path:
    output = summary_csv_path(args)
    output.parent.mkdir(parents=True, exist_ok=True)
    normalization_factor = global_normalization(cases)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        for case in sorted(cases, key=lambda item: (float(item.beta_deg), int(item.sorted_index))):
            raw_coeff = {name: float(value) for name, value in zip(COEFF_NAMES, case.raw_coeff, strict=True)}
            plot_coeff = {name: float(value) for name, value in zip(COEFF_NAMES, case.coeff, strict=True)}
            writer.writerow(
                {
                    "beta_deg": number_text(case.beta_deg),
                    "sorted_index": int(case.sorted_index),
                    "Lambda": number_text(case.Lambda),
                    "sign_normalization": case.sign_normalization,
                    "plot_sign_factor": number_text(case.plot_sign_factor),
                    "short_rod_id": int(case.short_rod_id),
                    "short_rod_reference_value_before_sign": number_text(case.short_rod_reference_value_before_sign),
                    "short_rod_reference_value_after_sign": number_text(case.short_rod_reference_value_after_sign),
                    "coeff_A1": number_text(plot_coeff["A1"]),
                    "coeff_B1": number_text(plot_coeff["B1"]),
                    "coeff_A2": number_text(plot_coeff["A2"]),
                    "coeff_B2": number_text(plot_coeff["B2"]),
                    "coeff_P1": number_text(plot_coeff["P1"]),
                    "coeff_P2": number_text(plot_coeff["P2"]),
                    "raw_coeff_A1": number_text(raw_coeff["A1"]),
                    "raw_coeff_B1": number_text(raw_coeff["B1"]),
                    "raw_coeff_A2": number_text(raw_coeff["A2"]),
                    "raw_coeff_B2": number_text(raw_coeff["B2"]),
                    "raw_coeff_P1": number_text(raw_coeff["P1"]),
                    "raw_coeff_P2": number_text(raw_coeff["P2"]),
                    "plot_coeff_A1": number_text(plot_coeff["A1"]),
                    "plot_coeff_B1": number_text(plot_coeff["B1"]),
                    "plot_coeff_A2": number_text(plot_coeff["A2"]),
                    "plot_coeff_B2": number_text(plot_coeff["B2"]),
                    "plot_coeff_P1": number_text(plot_coeff["P1"]),
                    "plot_coeff_P2": number_text(plot_coeff["P2"]),
                    "normalization_factor": number_text(normalization_factor),
                    "shape_max_amplitude_raw": number_text(case.shape_max_amplitude_raw),
                    "shape_max_amplitude_normalized": number_text(case.shape_max_amplitude_normalized),
                    "mac_to_previous_same_sorted": number_text(case.mac_to_previous_same_sorted),
                    "mac_to_beta0_same_sorted": number_text(case.mac_to_beta0_same_sorted),
                    "mac_sorted2_to_sorted3_at_beta": number_text(case.mac_sorted2_to_sorted3_at_beta),
                    "residual_norm": number_text(case.residual_norm),
                    "smallest_singular_value": number_text(case.smallest_singular_value),
                    "singular_value_ratio": number_text(case.singular_value_ratio),
                    "root_solver_warning": "no",
                    "warning": case.warning,
                    "notes": case.notes,
                }
            )
    return output


def global_normalization(cases: Sequence[SortedShapeCase]) -> float:
    values = [
        float(case.shape_max_amplitude_raw) / max(float(case.shape_max_amplitude_normalized), 1e-12)
        for case in cases
        if np.isfinite(case.shape_max_amplitude_normalized)
    ]
    return float(np.median(values)) if values else float("nan")


def write_gap_csv(args: argparse.Namespace, cases: Sequence[SortedShapeCase], beta_values: Sequence[float]) -> Path | None:
    if set(int(value) for value in args.sorted_modes) < {2, 3}:
        return None
    output = gap_csv_path(args, beta_values)
    output.parent.mkdir(parents=True, exist_ok=True)
    by_key = {(round(float(case.beta_deg), 10), int(case.sorted_index)): case for case in cases}
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=GAP_FIELDNAMES)
        writer.writeheader()
        for beta in beta_values:
            sorted2 = by_key[(round(float(beta), 10), 2)]
            sorted3 = by_key[(round(float(beta), 10), 3)]
            writer.writerow(
                {
                    "beta_deg": number_text(float(beta)),
                    "Lambda_sorted2": number_text(sorted2.Lambda),
                    "Lambda_sorted3": number_text(sorted3.Lambda),
                    "gap_23": number_text(sorted3.Lambda - sorted2.Lambda),
                    "mac_sorted2_to_sorted3": number_text(sorted2.mac_sorted2_to_sorted3_at_beta),
                    "mac_sorted2_to_previous": number_text(sorted2.mac_to_previous_same_sorted),
                    "mac_sorted3_to_previous": number_text(sorted3.mac_to_previous_same_sorted),
                }
            )
    return output


def plot_gap_png(args: argparse.Namespace, cases: Sequence[SortedShapeCase], beta_values: Sequence[float]) -> Path | None:
    if set(int(value) for value in args.sorted_modes) < {2, 3}:
        return None
    output = gap_png_path(args, beta_values)
    output.parent.mkdir(parents=True, exist_ok=True)
    by_key = {(round(float(case.beta_deg), 10), int(case.sorted_index)): case for case in cases}
    gaps = np.array(
        [by_key[(round(float(beta), 10), 3)].Lambda - by_key[(round(float(beta), 10), 2)].Lambda for beta in beta_values],
        dtype=float,
    )
    min_idx = int(np.nanargmin(gaps))
    fig, ax = plt.subplots(figsize=(8.4, 4.5))
    ax.plot(beta_values, gaps, color="#C85200", linewidth=1.9, marker="o", markersize=3)
    ax.scatter([float(beta_values[min_idx])], [float(gaps[min_idx])], color="black", s=28, zorder=5)
    ax.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
    ax.annotate(
        f"min gap={float(gaps[min_idx]):.6g}\nat beta={float(beta_values[min_idx]):g} deg",
        xy=(float(beta_values[min_idx]), float(gaps[min_idx])),
        xytext=(0.05, 0.76),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": "0.25", "linewidth": 0.8},
        fontsize=9,
    )
    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("Lambda_3^sorted - Lambda_2^sorted")
    ax.set_title("Sorted gap between fixed sorted modes 2 and 3")
    ax.grid(True, color="0.88", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)
    return output


def display_path(path: Path | None) -> str:
    if path is None:
        return "not generated"
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def case_key(beta_deg: float, sorted_index: int) -> tuple[float, int]:
    return round(float(beta_deg), 10), int(sorted_index)


def read_existing_summary_lambdas(path: Path) -> dict[tuple[float, int], float]:
    if not Path(path).exists():
        return {}
    values: dict[tuple[float, int], float] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                values[case_key(float(row["beta_deg"]), int(row["sorted_index"]))] = float(row["Lambda"])
            except (KeyError, TypeError, ValueError):
                continue
    return values


def read_existing_gap_values(path: Path) -> dict[float, float]:
    if not Path(path).exists():
        return {}
    values: dict[float, float] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                values[round(float(row["beta_deg"]), 10)] = float(row["gap_23"])
            except (KeyError, TypeError, ValueError):
                continue
    return values


def lambda_change_status(
    cases: Sequence[SortedShapeCase],
    previous_lambdas: dict[tuple[float, int], float],
    *,
    tolerance: float = 1e-10,
) -> tuple[str, float]:
    if not previous_lambdas:
        return "previous summary CSV not available", float("nan")
    deltas = []
    for case in cases:
        previous = previous_lambdas.get(case_key(float(case.beta_deg), int(case.sorted_index)))
        if previous is not None:
            deltas.append(abs(float(case.Lambda) - float(previous)))
    if not deltas:
        return "no matching previous Lambda rows", float("nan")
    max_delta = float(max(deltas))
    status = "unchanged" if max_delta <= tolerance else "changed"
    return status, max_delta


def gap_change_status(
    cases: Sequence[SortedShapeCase],
    beta_values: Sequence[float],
    previous_gaps: dict[float, float],
    *,
    tolerance: float = 1e-10,
) -> tuple[str, float]:
    if set(int(value) for value in {case.sorted_index for case in cases}) < {2, 3}:
        return "not applicable", float("nan")
    if not previous_gaps:
        return "previous gap CSV not available", float("nan")
    by_key = {(round(float(case.beta_deg), 10), int(case.sorted_index)): case for case in cases}
    deltas = []
    for beta in beta_values:
        key = round(float(beta), 10)
        previous = previous_gaps.get(key)
        if previous is None:
            continue
        current = by_key[(key, 3)].Lambda - by_key[(key, 2)].Lambda
        deltas.append(abs(float(current) - float(previous)))
    if not deltas:
        return "no matching previous gap rows", float("nan")
    max_delta = float(max(deltas))
    status = "unchanged" if max_delta <= tolerance else "changed"
    return status, max_delta


def print_sanity_summary(
    cases: Sequence[SortedShapeCase],
    beta_values: Sequence[float],
    *,
    previous_lambdas: dict[tuple[float, int], float],
    previous_gaps: dict[float, float],
) -> None:
    print("sign-normalization sanity:")
    after_values = [float(case.short_rod_reference_value_after_sign) for case in cases]
    min_after = float(np.nanmin(after_values)) if after_values else float("nan")
    print(f"short-rod reference min after sign: {number_text(min_after)}")
    print(f"after-sign >= -1e-12: {'yes' if min_after >= -SIGN_REFERENCE_AMBIGUOUS_TOL else 'no'}")
    lambda_status, lambda_delta = lambda_change_status(cases, previous_lambdas)
    gap_status, gap_delta = gap_change_status(cases, beta_values, previous_gaps)
    print(f"Lambda comparison with existing summary: {lambda_status}; max_abs_delta={number_text(lambda_delta)}")
    print(f"gap_23 comparison with existing gap CSV: {gap_status}; max_abs_delta={number_text(gap_delta)}")
    for sorted_index in sorted({int(case.sorted_index) for case in cases}):
        signs = [
            int(np.sign(float(case.plot_sign_factor)))
            for case in cases
            if int(case.sorted_index) == sorted_index
        ]
        negative_count = sum(1 for value in signs if value < 0)
        positive_count = sum(1 for value in signs if value > 0)
        print(f"sorted {sorted_index} plot_sign_factor summary: +1={positive_count}, -1={negative_count}")
    print("per-case sign-normalization rows:")
    for case in sorted(cases, key=lambda item: (float(item.beta_deg), int(item.sorted_index))):
        print(
            "  "
            f"beta={case.beta_deg:g}, sorted={int(case.sorted_index)}, "
            f"Lambda={case.Lambda:.12g}, sign={case.plot_sign_factor:+.0f}, "
            f"short_ref_before={case.short_rod_reference_value_before_sign:.12g}, "
            f"short_ref_after={case.short_rod_reference_value_after_sign:.12g}, "
            f"warning={case.warning}"
        )


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    beta_values = shape_scan.inclusive_grid(float(args.beta_start), float(args.beta_end), float(args.beta_step))
    sorted_roots = shape_scan.solve_sorted_root_grid(root_solve_namespace(args), beta_values)
    cases = solve_and_reconstruct(args, beta_values, sorted_roots)
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)
    axis_limits = axis_limits_for_cases(args, cases, s_norm)
    previous_lambdas = read_existing_summary_lambdas(summary_csv_path(args))
    previous_gaps = read_existing_gap_values(gap_csv_path(args, beta_values))

    individual_paths = plot_individuals(args, cases, s_norm, axis_limits)
    grid_paths = [
        plot_mode_grid(args, cases, s_norm, axis_limits, sorted_index=int(sorted_index), beta_values=beta_values)
        for sorted_index in args.sorted_modes
    ]
    comparison_path = plot_comparison_grid(args, cases, s_norm, axis_limits, beta_values=beta_values)
    summary_path = write_summary_csv(args, cases)
    gap_csv = write_gap_csv(args, cases, beta_values)
    gap_png = plot_gap_png(args, cases, beta_values)

    print("fixed sorted-index mode shape scan:")
    print(f"output directory: {display_path(Path(args.output_dir))}")
    print(f"individual PNG count: {len(individual_paths)}")
    for path in grid_paths:
        print(f"grid PNG: {display_path(path)}")
    print(f"comparison PNG: {display_path(comparison_path)}")
    print(f"summary CSV: {display_path(summary_path)}")
    print(f"gap CSV: {display_path(gap_csv)}")
    print(f"gap PNG: {display_path(gap_png)}")
    print("mode selection: fixed sorted index at each beta; descendant tracking not used")
    print_sanity_summary(
        cases,
        beta_values,
        previous_lambdas=previous_lambdas,
        previous_gaps=previous_gaps,
    )
    return {
        "cases": cases,
        "individual_paths": individual_paths,
        "grid_paths": grid_paths,
        "comparison_path": comparison_path,
        "summary_path": summary_path,
        "gap_csv": gap_csv,
        "gap_png": gap_png,
    }


if __name__ == "__main__":
    main()
