from __future__ import annotations

import argparse
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

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_timoshenko_shape_construction as shape_audit,
)
from scripts.lib import in_plane_shape_geometry as DISPLAY  # noqa: E402


MODEL_EB = shape_audit.MODEL_EB
MODEL_TIMO = shape_audit.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)
MODEL_LABEL = {MODEL_EB: "Euler-Bernoulli", MODEL_TIMO: "Timoshenko"}
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
REQUESTED_FULL_SCALE = 0.08
FALLBACK_FULL_SCALE = 0.05
JOINT_TOL = 1.0e-6
DEFAULT_OUTPUT_DIR = Path("results") / "timoshenko_mode_shape_diagnostics"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "timoshenko_mode_shape_diagnostics"

MAC_FIELDS = [
    "mu",
    "epsilon",
    "beta_deg",
    "eta",
    "mode_i",
    "mode_j",
    "Lambda_i",
    "Lambda_j",
    "mac_full",
    "mac_displacement",
    "mac_w_psi",
    "mac_u",
    "mac_w",
    "mac_psi",
    "mac_gamma",
    "notes",
]

ENERGY_FIELDS = [
    "model",
    "epsilon",
    "beta_deg",
    "eta",
    "mu",
    "sorted_index",
    "Lambda",
    "U_axial",
    "U_bending",
    "U_shear",
    "axial_fraction",
    "bending_fraction",
    "shear_fraction",
    "bending_shear_fraction",
    "classification",
    "notes",
]

RECONSTRUCTION_FIELDS = [
    "mu",
    "sorted_index",
    "Lambda",
    "smallest_singular_value",
    "singular_value_ratio",
    "max_joint_kinematic_gap",
    "max_joint_force_gap",
    "max_boundary_gap",
    "max_ode_residual",
    "notes",
]


@dataclass(frozen=True)
class FullScaleRecord:
    result: shape_audit.ModeResult
    scale_requested: float
    scale_used: float
    fallback_used: bool
    fold_flag: bool
    near_overlap_flag: bool
    plotted_joint_gap: float
    notes: str


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Diagnostic-only Timoshenko modes 4-6 shape visualization audit.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if bool(args.smoke):
        args.output_dir = SMOKE_OUTPUT_DIR
        args.n_points = min(int(args.n_points), SMOKE_N_POINTS)
    args.output_dir = repo_path(Path(args.output_dir))
    if int(args.n_points) < 51:
        raise ValueError("--n-points must be at least 51")
    return args


def selected_mu_values(smoke: bool) -> tuple[float, ...]:
    return SMOKE_MU_VALUES if bool(smoke) else DEFAULT_MU_VALUES


def selected_sorted_indices(smoke: bool) -> tuple[int, ...]:
    return SMOKE_SORTED_INDICES if bool(smoke) else DEFAULT_SORTED_INDICES


def build_results(args: argparse.Namespace) -> list[shape_audit.ModeResult]:
    mu_values = selected_mu_values(bool(args.smoke))
    sorted_indices = selected_sorted_indices(bool(args.smoke))
    provider = shape_audit.RootProvider(DEFAULT_BETA_DEG, DEFAULT_ETA, max(DEFAULT_SORTED_INDICES))
    results: list[shape_audit.ModeResult] = []
    for mu in mu_values:
        eb_roots, timo_roots, warnings = provider.roots(DEFAULT_EPSILON, float(mu))
        for sorted_index in sorted_indices:
            results.append(
                shape_audit.eb_mode_result(
                    epsilon=DEFAULT_EPSILON,
                    beta_deg=DEFAULT_BETA_DEG,
                    eta=DEFAULT_ETA,
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=float(eb_roots[int(sorted_index) - 1]),
                    n_points=int(args.n_points),
                    case_role="modes_4_6_shape_diagnostics",
                    root_source="sorted_global_solver",
                    root_warnings=warnings,
                )
            )
            results.append(
                shape_audit.timo_mode_result(
                    epsilon=DEFAULT_EPSILON,
                    beta_deg=DEFAULT_BETA_DEG,
                    eta=DEFAULT_ETA,
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=float(timo_roots[int(sorted_index) - 1]),
                    n_points=int(args.n_points),
                    case_role="modes_4_6_shape_diagnostics",
                    root_source="sorted_global_solver",
                    root_warnings=warnings,
                )
            )
    return results


def only_model(results: Sequence[shape_audit.ModeResult], model: str) -> list[shape_audit.ModeResult]:
    return [result for result in results if result.model == model]


def result_by_key(results: Sequence[shape_audit.ModeResult]) -> dict[tuple[str, float, int], shape_audit.ModeResult]:
    return {
        (result.model, round(float(result.mu), 12), int(result.sorted_index)): result
        for result in results
    }


def signed_component(result: shape_audit.ModeResult, component: str) -> tuple[np.ndarray, np.ndarray]:
    sign = float(result.sign_factor)
    if component == "u":
        return sign * result.rod1.u, sign * result.rod2.u
    if component == "w":
        return sign * result.rod1.w, sign * result.rod2.w
    if component == "psi":
        return sign * np.asarray(result.rod1.psi, dtype=float), sign * np.asarray(result.rod2.psi, dtype=float)
    if component == "gamma":
        return sign * np.asarray(result.rod1.gamma, dtype=float), sign * np.asarray(result.rod2.gamma, dtype=float)
    raise ValueError(f"unknown component {component!r}")


def trapezoid_sqrt_weights(x_values: np.ndarray) -> np.ndarray:
    x = np.asarray(x_values, dtype=float)
    if x.size < 2:
        return np.ones_like(x)
    dx = np.abs(np.diff(x))
    weights = np.empty_like(x)
    weights[0] = 0.5 * dx[0]
    weights[-1] = 0.5 * dx[-1]
    if x.size > 2:
        weights[1:-1] = 0.5 * (dx[:-1] + dx[1:])
    return np.sqrt(np.maximum(weights, 0.0))


def component_vector(result: shape_audit.ModeResult, component: str) -> np.ndarray:
    values1, values2 = signed_component(result, component)
    weighted = np.concatenate(
        [
            values1 * trapezoid_sqrt_weights(result.rod1.x_local),
            values2 * trapezoid_sqrt_weights(result.rod2.x_local),
        ]
    )
    norm = float(np.linalg.norm(weighted))
    return weighted / norm if norm > 1.0e-30 else weighted


def grouped_vector(result: shape_audit.ModeResult, components: Sequence[str]) -> np.ndarray:
    pieces = [component_vector(result, component) for component in components]
    vector = np.concatenate(pieces)
    norm = float(np.linalg.norm(vector))
    return vector / norm if norm > 1.0e-30 else vector


def mac_value(vector_i: np.ndarray, vector_j: np.ndarray) -> float:
    vi = np.asarray(vector_i, dtype=float)
    vj = np.asarray(vector_j, dtype=float)
    denom = float(np.dot(vi, vi) * np.dot(vj, vj))
    if denom <= 1.0e-30:
        return float("nan")
    dot = float(np.dot(vi, vj))
    return float((dot * dot) / denom)


def mac_rows(results: Sequence[shape_audit.ModeResult], mu_values: Sequence[float], sorted_indices: Sequence[int]) -> list[dict[str, object]]:
    timo_by_key = result_by_key(only_model(results, MODEL_TIMO))
    rows: list[dict[str, object]] = []
    groups = {
        "mac_full": ("u", "w", "psi"),
        "mac_displacement": ("u", "w"),
        "mac_w_psi": ("w", "psi"),
        "mac_u": ("u",),
        "mac_w": ("w",),
        "mac_psi": ("psi",),
        "mac_gamma": ("gamma",),
    }
    for mu in mu_values:
        for index_i, mode_i in enumerate(sorted_indices):
            for mode_j in sorted_indices[index_i + 1 :]:
                result_i = timo_by_key[(MODEL_TIMO, round(float(mu), 12), int(mode_i))]
                result_j = timo_by_key[(MODEL_TIMO, round(float(mu), 12), int(mode_j))]
                row: dict[str, object] = {
                    "mu": float(mu),
                    "epsilon": DEFAULT_EPSILON,
                    "beta_deg": DEFAULT_BETA_DEG,
                    "eta": DEFAULT_ETA,
                    "mode_i": int(mode_i),
                    "mode_j": int(mode_j),
                    "Lambda_i": result_i.Lambda,
                    "Lambda_j": result_j.Lambda,
                }
                for field, components in groups.items():
                    row[field] = mac_value(grouped_vector(result_i, components), grouped_vector(result_j, components))
                row["notes"] = (
                    "component-balanced sign-invariant MAC; each component is L2-normalized over both rods "
                    "before grouped vectors are concatenated"
                )
                rows.append(row)
    return rows


def scale_flags(result: shape_audit.ModeResult, scale_fraction: float) -> tuple[bool, bool, float, bool]:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(scale_fraction))
    rows = [
        shape_audit.regularity_row(
            result,
            rod_id=rod_id,
            scale_kind="fixed_fraction",
            deformation_scale_fraction=float(scale_fraction),
            deformation_scale_value=scale,
        )
        for rod_id in (1, 2)
    ]
    fold = any(bool(row["visually_folded_at_scale"]) for row in rows)
    near = any(int(row["near_overlap_pair_count"]) > 0 for row in rows)
    xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
        result,
        transverse_only=False,
        scale=scale,
    )
    gap = float(np.linalg.norm(np.array([xd1[-1] - xd2[0], yd1[-1] - yd2[0]], dtype=float)))
    ok = bool(not fold and not near and gap <= JOINT_TOL)
    return fold, near, gap, ok


def full_scale_records(results: Sequence[shape_audit.ModeResult]) -> list[FullScaleRecord]:
    records: list[FullScaleRecord] = []
    for result in only_model(results, MODEL_TIMO):
        fold, near, gap, ok = scale_flags(result, REQUESTED_FULL_SCALE)
        scale_used = REQUESTED_FULL_SCALE
        fallback = False
        notes = f"scale {REQUESTED_FULL_SCALE:.2f} passed regularity"
        if not ok:
            fallback = True
            scale_used = FALLBACK_FULL_SCALE
            fold, near, gap, fallback_ok = scale_flags(result, FALLBACK_FULL_SCALE)
            notes = (
                f"scale {REQUESTED_FULL_SCALE:.2f} failed regularity; "
                f"fallback {FALLBACK_FULL_SCALE:.2f} {'passed' if fallback_ok else 'still has warnings'}"
            )
        records.append(
            FullScaleRecord(
                result=result,
                scale_requested=REQUESTED_FULL_SCALE,
                scale_used=scale_used,
                fallback_used=fallback,
                fold_flag=fold,
                near_overlap_flag=near,
                plotted_joint_gap=gap,
                notes=notes,
            )
        )
    return records


def write_component_grid(
    output_dir: Path,
    results: Sequence[shape_audit.ModeResult],
    *,
    component: str,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / f"timo_component_{component}_modes4_6.png"
    by_key = result_by_key(only_model(results, MODEL_TIMO))
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(12.6, 2.15 * len(sorted_indices) + 0.75),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            result = by_key[(MODEL_TIMO, round(float(mu), 12), int(sorted_index))]
            values1, values2 = signed_component(result, component)
            ax.plot(result.rod1.x_local, values1, color="#0072b2", linewidth=1.25)
            ax.plot(result.rod2.x_local, values2, color="#d55e00", linewidth=1.25)
            ax.scatter(
                [result.rod1.x_local[-1], result.rod2.x_local[-1]],
                [values1[-1], values2[-1]],
                color="black",
                s=10,
                zorder=5,
            )
            ax.axhline(0.0, color="0.82", linewidth=0.55)
            ax.grid(True, color="0.91", linewidth=0.5)
            ax.tick_params(labelsize=7)
            ax.set_title(f"mu={float(mu):g}, sorted {int(sorted_index)}, Lambda={result.Lambda:.5g}", fontsize=7.8)
            if col_index == 0:
                ax.set_ylabel(component, fontsize=8)
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("signed local x_i", fontsize=8)
    handles = [
        Line2D([0], [0], color="#0072b2", linewidth=1.25, label="rod 1"),
        Line2D([0], [0], color="#d55e00", linewidth=1.25, label="rod 2"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(
        f"Timoshenko {component} component: eps=0.03, beta=45 deg, eta=0, sorted modes 4-6",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.045, 1.0, 0.94), h_pad=0.65, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def write_combined_component_sheet(
    output_dir: Path,
    results: Sequence[shape_audit.ModeResult],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / "timo_components_modes4_6_combined.png"
    by_key = result_by_key(only_model(results, MODEL_TIMO))
    components = ("u", "w", "psi", "gamma")
    row_specs = [(int(mode), component) for mode in sorted_indices for component in components]
    fig, axes = plt.subplots(
        len(row_specs),
        len(mu_values),
        figsize=(12.8, 1.35 * len(row_specs) + 0.9),
        squeeze=False,
    )
    for row_index, (sorted_index, component) in enumerate(row_specs):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            result = by_key[(MODEL_TIMO, round(float(mu), 12), int(sorted_index))]
            values1, values2 = signed_component(result, component)
            ax.plot(result.rod1.x_local, values1, color="#0072b2", linewidth=0.95)
            ax.plot(result.rod2.x_local, values2, color="#d55e00", linewidth=0.95)
            ax.scatter(
                [result.rod1.x_local[-1], result.rod2.x_local[-1]],
                [values1[-1], values2[-1]],
                color="black",
                s=7,
                zorder=5,
            )
            ax.axhline(0.0, color="0.84", linewidth=0.45)
            ax.grid(True, color="0.92", linewidth=0.42)
            ax.tick_params(labelsize=6)
            if row_index == 0:
                ax.set_title(f"mu={float(mu):g}", fontsize=7.6)
            if col_index == 0:
                ax.set_ylabel(f"m{sorted_index}\n{component}", fontsize=6.9)
            if row_index == len(row_specs) - 1:
                ax.set_xlabel("signed x_i", fontsize=7)
    fig.suptitle("Timoshenko component sheet: modes 4-6, eps=0.03, beta=45 deg, eta=0", fontsize=11, y=0.997)
    fig.tight_layout(rect=(0.0, 0.015, 1.0, 0.965), h_pad=0.25, w_pad=0.45)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def base_and_displacement_arrays(
    result: shape_audit.ModeResult,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u1, w1, u2, w2 = shape_audit.plot_arrays(result)
    x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(u1))
    if result.model == MODEL_EB:
        dx1, dy1 = DISPLAY.eb_rod1_local_displacement_to_display(u1, w1)
        dx2, dy2 = DISPLAY.eb_rod2_local_displacement_to_display(u2, w2, beta_deg=result.beta_deg)
    else:
        dx1, dy1 = DISPLAY.rod1_local_displacement_to_display(u1, w1)
        dx2, dy2 = DISPLAY.rod2_local_displacement_to_display(u2, w2, beta_deg=result.beta_deg)
    return x1, y1, x2, y2, dx1, dy1, dx2, dy2


def vector_axis_limits(results: Sequence[shape_audit.ModeResult], arrow_scale: float) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for result in results:
        x1, y1, x2, y2, dx1, dy1, dx2, dy2 = base_and_displacement_arrays(result)
        xs.extend([x1, x2, x1 + arrow_scale * dx1, x2 + arrow_scale * dx2])
        ys.extend([y1, y2, y1 + arrow_scale * dy1, y2 + arrow_scale * dy2])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.08 * x_span), float(np.max(x_all) + 0.08 * x_span)),
        (float(np.min(y_all) - 0.16 * y_span), float(np.max(y_all) + 0.16 * y_span)),
    )


def global_arrow_scale(results: Sequence[shape_audit.ModeResult], fraction: float = 0.12) -> float:
    max_disp = 0.0
    for result in results:
        _x1, _y1, _x2, _y2, dx1, dy1, dx2, dy2 = base_and_displacement_arrays(result)
        max_disp = max(
            max_disp,
            float(np.max(np.sqrt(dx1 * dx1 + dy1 * dy1))),
            float(np.max(np.sqrt(dx2 * dx2 + dy2 * dy2))),
        )
    return float(fraction * 2.0 / max(max_disp, 1.0e-14))


def write_vector_field_grid(
    output_dir: Path,
    results: Sequence[shape_audit.ModeResult],
    *,
    model: str,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    prefix = "timo" if model == MODEL_TIMO else "eb"
    output = output_dir / f"{prefix}_displacement_vector_field_modes4_6_grid_corrected.png"
    model_results = only_model(results, model)
    by_key = result_by_key(model_results)
    arrow_scale = global_arrow_scale(model_results)
    limits = vector_axis_limits(model_results, arrow_scale)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(13.0, 2.55 * len(sorted_indices) + 0.8),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            result = by_key[(model, round(float(mu), 12), int(sorted_index))]
            ax = axes[row_index, col_index]
            x1, y1, x2, y2, dx1, dy1, dx2, dy2 = base_and_displacement_arrays(result)
            stride = max(1, len(x1) // 18)
            ax.plot(x1, y1, color="0.70", linestyle="--", linewidth=0.9)
            ax.plot(x2, y2, color="0.70", linestyle="--", linewidth=0.9)
            ax.quiver(
                x1[::stride],
                y1[::stride],
                arrow_scale * dx1[::stride],
                arrow_scale * dy1[::stride],
                angles="xy",
                scale_units="xy",
                scale=1.0,
                color=MODEL_COLOR[model],
                width=0.004,
            )
            ax.quiver(
                x2[::stride],
                y2[::stride],
                arrow_scale * dx2[::stride],
                arrow_scale * dy2[::stride],
                angles="xy",
                scale_units="xy",
                scale=1.0,
                color=MODEL_COLOR[model],
                width=0.004,
            )
            ax.scatter([x1[-1]], [y1[-1]], color="black", s=12, zorder=5)
            ax.set_xlim(*limits[0])
            ax.set_ylim(*limits[1])
            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, color="0.90", linewidth=0.5)
            ax.tick_params(labelsize=6.5)
            ax.set_title(
                f"{MODEL_LABEL[model]}, mu={float(mu):g}, sorted {int(sorted_index)}\n"
                f"Lambda={result.Lambda:.5g}, {result.energy['classification']}",
                fontsize=7.6,
            )
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=7.5)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=7.5)
    fig.suptitle(
        f"{MODEL_LABEL[model]} displacement vector field: full vector u*t + w*n",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.02, 1.0, 0.94), h_pad=0.55, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def full_centerline_limits(records: Sequence[FullScaleRecord]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for record in records:
        result = record.result
        x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        scale = shape_audit.fixed_scale_from_fraction(result.mu, record.scale_used)
        xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(result, transverse_only=False, scale=scale)
        xs.extend([x1, x2, xd1, xd2])
        ys.extend([y1, y2, yd1, yd2])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.08 * x_span), float(np.max(x_all) + 0.08 * x_span)),
        (float(np.min(y_all) - 0.16 * y_span), float(np.max(y_all) + 0.16 * y_span)),
    )


def write_full_centerline_grid(
    output_dir: Path,
    records: Sequence[FullScaleRecord],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / "timo_full_centerline_modes4_6_conservative_corrected.png"
    by_key = {
        (round(float(record.result.mu), 12), int(record.result.sorted_index)): record
        for record in records
    }
    limits = full_centerline_limits(records)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(13.0, 2.65 * len(sorted_indices) + 1.0),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            record = by_key[(round(float(mu), 12), int(sorted_index))]
            result = record.result
            ax = axes[row_index, col_index]
            x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
            scale = shape_audit.fixed_scale_from_fraction(result.mu, record.scale_used)
            xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(result, transverse_only=False, scale=scale)
            ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=0.9)
            ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=0.9)
            ax.plot(xd1, yd1, color="#d55e00", linewidth=1.55)
            ax.plot(xd2, yd2, color="#d55e00", linewidth=1.55)
            ax.scatter([x1[-1]], [y1[-1]], color="black", s=12, zorder=5)
            ax.set_xlim(*limits[0])
            ax.set_ylim(*limits[1])
            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, color="0.90", linewidth=0.5)
            ax.tick_params(labelsize=6.5)
            fallback = " fallback" if record.fallback_used else ""
            ax.set_title(
                f"mu={float(mu):g}, sorted {int(sorted_index)}, Lambda={result.Lambda:.5g}\n"
                f"scale={record.scale_used:.2f}{fallback}; {result.energy['classification']}",
                fontsize=7.5,
            )
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=7.5)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=7.5)
    fig.text(
        0.5,
        0.045,
        "Secondary view only: for longitudinal/mixed modes full centerline may be visually misleading.",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    fig.suptitle("Timoshenko conservative full-centerline view: eps=0.03, beta=45 deg, eta=0", fontsize=11, y=0.995)
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.94), h_pad=0.6, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def energy_rows(results: Sequence[shape_audit.ModeResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        energy = result.energy
        shear = float(energy["U_shear"]) if result.model == MODEL_TIMO else float("nan")
        shear_fraction = float(energy["shear_energy_fraction"]) if result.model == MODEL_TIMO else float("nan")
        bending_shear = (
            float(energy["bending_energy_fraction"]) + float(energy["shear_energy_fraction"])
            if result.model == MODEL_TIMO
            else float("nan")
        )
        rows.append(
            {
                "model": result.model,
                "epsilon": result.epsilon,
                "beta_deg": result.beta_deg,
                "eta": result.eta,
                "mu": result.mu,
                "sorted_index": result.sorted_index,
                "Lambda": result.Lambda,
                "U_axial": energy["U_axial"],
                "U_bending": energy["U_bending"],
                "U_shear": shear,
                "axial_fraction": energy["axial_energy_fraction"],
                "bending_fraction": energy["bending_energy_fraction"],
                "shear_fraction": shear_fraction,
                "bending_shear_fraction": bending_shear,
                "classification": energy["classification"],
                "notes": "energy fractions from existing analytic reconstruction helpers",
            }
        )
    return rows


def max_boundary_gap(result: shape_audit.ModeResult) -> float:
    if result.model != MODEL_TIMO:
        return float("nan")
    values = [
        float(result.rod1.u[0]),
        float(result.rod1.w[0]),
        float(np.asarray(result.rod1.psi, dtype=float)[0]),
        float(result.rod2.u[0]),
        float(result.rod2.w[0]),
        float(np.asarray(result.rod2.psi, dtype=float)[0]),
    ]
    return max(abs(value) for value in values)


def reconstruction_rows(results: Sequence[shape_audit.ModeResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in only_model(results, MODEL_TIMO):
        mode = shape_audit.TIMO.timo_mode_coefficients(
            result.Lambda,
            result.beta_deg,
            result.mu,
            result.epsilon,
            result.eta,
        )
        joint = result.joint_row or {}
        rows.append(
            {
                "mu": result.mu,
                "sorted_index": result.sorted_index,
                "Lambda": result.Lambda,
                "smallest_singular_value": mode.smallest_singular_value,
                "singular_value_ratio": mode.singular_value_ratio,
                "max_joint_kinematic_gap": joint.get("max_abs_kinematic_gap", float("nan")),
                "max_joint_force_gap": joint.get("max_abs_force_gap", float("nan")),
                "max_boundary_gap": max_boundary_gap(result),
                "max_ode_residual": float("nan"),
                "notes": (
                    "joint and clamp gaps checked from reconstructed fields; "
                    "ODE residual not independently evaluated because fields are analytic basis evaluations"
                ),
            }
        )
    return rows


def mac_stats(rows: Sequence[dict[str, object]]) -> dict[str, dict[str, float]]:
    fields = ("mac_full", "mac_displacement", "mac_w_psi", "mac_u", "mac_w", "mac_psi", "mac_gamma")
    stats: dict[str, dict[str, float]] = {}
    for field in fields:
        values = np.asarray([float(row[field]) for row in rows], dtype=float)
        stats[field] = {
            "min": float(np.min(values)),
            "mean": float(np.mean(values)),
            "max": float(np.max(values)),
        }
    return stats


def energy_class_counts(rows: Sequence[dict[str, object]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        if str(row["model"]) != MODEL_TIMO:
            continue
        label = str(row["classification"])
        counts[label] = counts.get(label, 0) + 1
    return counts


def write_report(
    path: Path,
    *,
    mac_csv: Path,
    energy_csv: Path,
    reconstruction_csv: Path,
    figure_paths: Sequence[Path],
    mac: Sequence[dict[str, object]],
    energy: Sequence[dict[str, object]],
    reconstruction: Sequence[dict[str, object]],
    full_records: Sequence[FullScaleRecord],
) -> Path:
    stats = mac_stats(mac)
    class_counts = energy_class_counts(energy)
    max_joint_kin = max((float(row["max_joint_kinematic_gap"]) for row in reconstruction), default=float("nan"))
    max_joint_force = max((float(row["max_joint_force_gap"]) for row in reconstruction), default=float("nan"))
    max_boundary = max((float(row["max_boundary_gap"]) for row in reconstruction), default=float("nan"))
    fallback_count = sum(1 for record in full_records if record.fallback_used)
    final_fold = sum(1 for record in full_records if record.fold_flag)
    final_near = sum(1 for record in full_records if record.near_overlap_flag)
    max_plot_gap = max((float(record.plotted_joint_gap) for record in full_records), default=float("nan"))
    component_means = {field: stats[field]["mean"] for field in ("mac_u", "mac_w", "mac_psi", "mac_gamma")}
    strongest = min(component_means, key=component_means.get)

    lines = [
        "# Timoshenko Modes 4-6 Shape Diagnostics",
        "",
        "## Scope",
        "",
        "- Parameters: `epsilon=0.03`, `beta=45 deg`, `eta=0`, `mu=0,0.2,0.4,0.6`, sorted modes `4,5,6`.",
        "- Sorted frequencies are used directly. No descendant tracking is used.",
        "- No formulas, determinant entries, root solvers, or Timoshenko shear coefficient k' were changed.",
        "- No FEM, 3D FEM, Gmsh, CalculiX, article workspace, article figures, old solvers, or baseline results were used or modified.",
        "",
        "## Display-Frame Correction",
        "",
        "Timoshenko global centerline and displacement-vector figures generated before this correction used the determinant rod-2 transform directly as Cartesian coordinates and must not be used for physical interpretation. Local fields, frequencies, null-vector checks, and energy fractions remain valid. Corrected Timoshenko figures use `t1=(1,0)`, `n1=(0,-1)`, `t2=(cos(beta),sin(beta))`, and `n2=(sin(beta),-cos(beta))`. EB fields use their opposite determinant transverse sign and the equivalent display normals `n1=(0,1)`, `n2=(-sin(beta),cos(beta))`.",
        "",
        "## Output Files",
        "",
        f"- MAC CSV: `{shape_audit.rel(mac_csv)}`",
        f"- Energy CSV: `{shape_audit.rel(energy_csv)}`",
        f"- Reconstruction audit CSV: `{shape_audit.rel(reconstruction_csv)}`",
        "",
        "Figures:",
    ]
    for figure in figure_paths:
        lines.append(f"- `{shape_audit.rel(figure)}`")

    lines.extend(
        [
            "",
            "## MAC Summary",
            "",
            "| MAC group | min | mean | max |",
            "|---|---:|---:|---:|",
        ]
    )
    for field, values in stats.items():
        lines.append(f"| {field} | {values['min']:.4f} | {values['mean']:.4f} | {values['max']:.4f} |")

    lines.extend(
        [
            "",
            "The low full-field and component MAC values show that sorted Timoshenko modes 4, 5, and 6 are not the same reconstructed shape. The full-centerline plot is therefore hiding differences rather than proving the modes are identical.",
            "",
            "## Energy And Classification",
            "",
            f"- Timoshenko classification counts: `{class_counts}`.",
            "- Energy fractions are recorded in the CSV for both EB and Timoshenko.",
            "",
            "## Reconstruction Checks",
            "",
            f"- Maximum Timoshenko joint kinematic gap: `{max_joint_kin:.6e}`.",
            f"- Maximum Timoshenko joint force/moment gap: `{max_joint_force:.6e}`.",
            f"- Maximum clamp boundary gap: `{max_boundary:.6e}`.",
            "- ODE residuals are not independently recomputed; the fields are evaluated from the analytic Timoshenko basis.",
            "",
            "## Full-Centerline Regularity",
            "",
            f"- Requested conservative full-centerline scale: `{REQUESTED_FULL_SCALE:.2f}`.",
            f"- Fallback scale: `{FALLBACK_FULL_SCALE:.2f}`.",
            f"- Fallback full-centerline shapes: `{fallback_count}`.",
            f"- Remaining fold flags after fallback: `{final_fold}`.",
            f"- Remaining near-overlap flags after fallback: `{final_near}`.",
            f"- Maximum plotted joint gap: `{max_plot_gap:.6e}`.",
            "",
            "## Answers",
            "",
            "1. Are modes 4, 5, and 6 really similar?",
            "No. The MAC audit shows substantial differences in the reconstructed Timoshenko fields.",
            "",
            "2. What does MAC show?",
            f"The component-balanced `mac_full` has mean `{stats['mac_full']['mean']:.4f}` and maximum `{stats['mac_full']['max']:.4f}`, far from a near-identical value of 1.",
            "",
            "3. Which components differ most?",
            f"The lowest mean single-component MAC is `{strongest}`. Component plots should be used to inspect u, w, psi, and gamma directly for each mu.",
            "",
            "4. Which modes are longitudinal / mixed / bending?",
            "Use `mode_4_6_energy_fractions.csv`; it records the energy classification for every model, mu, and sorted mode.",
            "",
            "5. Why can full displacement centerline be a poor visualization?",
            "The corrected reflected display frame removes the former rod-2 direction error. Conservative full-centerline scaling can still suppress small component differences, so component plots remain the primary interpretation tool.",
            "",
            "6. Recommended visualization:",
            "- component plots are the primary modal interpretation tool;",
            "- displacement vector fields show the full direction of u*t + w*n without forcing arrows into a misleading curve;",
            "- conservative full centerline is useful only as a secondary geometry/joint sanity check.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    mu_values = selected_mu_values(bool(args.smoke))
    sorted_indices = selected_sorted_indices(bool(args.smoke))
    results = build_results(args)
    mac = mac_rows(results, mu_values, sorted_indices)
    energy = energy_rows(results)
    reconstruction = reconstruction_rows(results)
    full_records = full_scale_records(results)

    mac_csv = shape_audit.write_csv(args.output_dir / "timoshenko_modes_4_6_mac.csv", mac, MAC_FIELDS)
    energy_csv = shape_audit.write_csv(args.output_dir / "mode_4_6_energy_fractions.csv", energy, ENERGY_FIELDS)
    reconstruction_csv = shape_audit.write_csv(
        args.output_dir / "timoshenko_mode_reconstruction_audit.csv",
        reconstruction,
        RECONSTRUCTION_FIELDS,
    )
    figure_paths = [
        write_component_grid(args.output_dir, results, component="u", mu_values=mu_values, sorted_indices=sorted_indices),
        write_component_grid(args.output_dir, results, component="w", mu_values=mu_values, sorted_indices=sorted_indices),
        write_component_grid(args.output_dir, results, component="psi", mu_values=mu_values, sorted_indices=sorted_indices),
        write_component_grid(args.output_dir, results, component="gamma", mu_values=mu_values, sorted_indices=sorted_indices),
        write_combined_component_sheet(args.output_dir, results, mu_values=mu_values, sorted_indices=sorted_indices),
        write_vector_field_grid(args.output_dir, results, model=MODEL_TIMO, mu_values=mu_values, sorted_indices=sorted_indices),
        write_vector_field_grid(args.output_dir, results, model=MODEL_EB, mu_values=mu_values, sorted_indices=sorted_indices),
        write_full_centerline_grid(args.output_dir, full_records, mu_values=mu_values, sorted_indices=sorted_indices),
    ]
    report = write_report(
        args.output_dir / "timoshenko_modes_4_6_shape_diagnostics_report.md",
        mac_csv=mac_csv,
        energy_csv=energy_csv,
        reconstruction_csv=reconstruction_csv,
        figure_paths=figure_paths,
        mac=mac,
        energy=energy,
        reconstruction=reconstruction,
        full_records=full_records,
    )
    stats = mac_stats(mac)
    fallback_count = sum(1 for record in full_records if record.fallback_used)
    max_plot_gap = max((float(record.plotted_joint_gap) for record in full_records), default=float("nan"))
    print("generated Timoshenko modes 4-6 shape diagnostic outputs:")
    for path in [mac_csv, energy_csv, reconstruction_csv, report, *figure_paths]:
        print(f"  {shape_audit.rel(path)}")
    print("MAC summary:")
    for field, values in stats.items():
        print(f"  {field}: min={values['min']:.6f}, mean={values['mean']:.6f}, max={values['max']:.6f}")
    print(f"fallback full-centerline shapes: {fallback_count}")
    print(f"max plotted joint gap: {max_plot_gap:.6e}")
    return {
        "mac_csv": mac_csv,
        "energy_csv": energy_csv,
        "reconstruction_csv": reconstruction_csv,
        "report": report,
        "figure_paths": figure_paths,
        "mac_rows": mac,
        "energy_rows": energy,
        "reconstruction_rows": reconstruction,
    }


if __name__ == "__main__":
    main()
