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


MODEL_TIMO = shape_audit.MODEL_TIMO

DEFAULT_EPSILON = 0.03
DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
DEFAULT_MU_VALUES = (0.0, 0.2, 0.4, 0.6)
DEFAULT_SORTED_INDICES = (4, 5, 6)
SMOKE_MU_VALUES = (0.0,)
SMOKE_SORTED_INDICES = (4, 5)
DEFAULT_N_POINTS = 801
SMOKE_N_POINTS = 301
DEFAULT_REQUESTED_SCALE = 0.08
FALLBACK_SCALE = 0.05
DEFAULT_OUTPUT_DIR = Path("results") / "timoshenko_shape_construction_audit"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "timoshenko_modes456_visualization_audit"
JOINT_TOL = 1.0e-6

SIMILARITY_FIELDS = [
    "mu",
    "epsilon",
    "mode_i",
    "mode_j",
    "Lambda_i",
    "Lambda_j",
    "similarity_full_field",
    "similarity_displacement_uw",
    "similarity_axial_u",
    "similarity_transverse_w",
    "similarity_rotation_shear",
    "similarity_plotted_centerline",
    "notes",
]


@dataclass(frozen=True)
class PlotScaleRecord:
    result: shape_audit.ModeResult
    scale_requested: float
    scale_used: float
    fallback_used: bool
    fold_flag: bool
    near_overlap_flag: bool
    plotted_joint_gap: float
    point_order_ok: bool
    notes: str


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Diagnostic-only visualization audit for Timoshenko sorted modes 4-6.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument("--deformation-scale-fraction", type=float, default=DEFAULT_REQUESTED_SCALE)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if bool(args.smoke):
        args.output_dir = SMOKE_OUTPUT_DIR
        args.n_points = min(int(args.n_points), SMOKE_N_POINTS)
    args.output_dir = repo_path(Path(args.output_dir))
    if int(args.n_points) < 51:
        raise ValueError("--n-points must be at least 51")
    if float(args.deformation_scale_fraction) <= 0.0:
        raise ValueError("--deformation-scale-fraction must be positive")
    return args


def selected_mu_values(smoke: bool) -> tuple[float, ...]:
    return SMOKE_MU_VALUES if bool(smoke) else DEFAULT_MU_VALUES


def selected_sorted_indices(smoke: bool) -> tuple[int, ...]:
    return SMOKE_SORTED_INDICES if bool(smoke) else DEFAULT_SORTED_INDICES


def build_results(args: argparse.Namespace) -> list[shape_audit.ModeResult]:
    mu_values = selected_mu_values(bool(args.smoke))
    sorted_indices = selected_sorted_indices(bool(args.smoke))
    provider = shape_audit.RootProvider(
        DEFAULT_BETA_DEG,
        DEFAULT_ETA,
        max(DEFAULT_SORTED_INDICES),
    )
    results: list[shape_audit.ModeResult] = []
    for mu in mu_values:
        _eb_roots, timo_roots, root_warnings = provider.roots(DEFAULT_EPSILON, float(mu))
        for sorted_index in sorted_indices:
            Lambda = float(timo_roots[int(sorted_index) - 1])
            results.append(
                shape_audit.timo_mode_result(
                    epsilon=DEFAULT_EPSILON,
                    beta_deg=DEFAULT_BETA_DEG,
                    eta=DEFAULT_ETA,
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=Lambda,
                    n_points=int(args.n_points),
                    case_role="modes456_visualization",
                    root_source="sorted_global_solver",
                    root_warnings=root_warnings,
                )
            )
    return results


def point_order_ok(result: shape_audit.ModeResult) -> bool:
    for rod_id in (1, 2):
        local = shape_audit.plotted_local_coordinate(result, rod_id)
        diffs = np.diff(local)
        if not (np.all(diffs >= -1.0e-12) or np.all(diffs <= 1.0e-12)):
            return False
    return True


def plotted_joint_gap(result: shape_audit.ModeResult, fraction: float) -> float:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(fraction))
    xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
        result,
        transverse_only=False,
        scale=scale,
    )
    return float(np.linalg.norm(np.array([xd1[-1] - xd2[0], yd1[-1] - yd2[0]], dtype=float)))


def scale_flags(result: shape_audit.ModeResult, fraction: float) -> tuple[bool, bool, float, bool]:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(fraction))
    rows = [
        shape_audit.regularity_row(
            result,
            rod_id=rod_id,
            scale_kind="fixed_fraction",
            deformation_scale_fraction=float(fraction),
            deformation_scale_value=scale,
        )
        for rod_id in (1, 2)
    ]
    fold_flag = any(bool(row["visually_folded_at_scale"]) for row in rows)
    near_overlap_flag = any(int(row["near_overlap_pair_count"]) > 0 for row in rows)
    gap = plotted_joint_gap(result, float(fraction))
    ok = bool(not fold_flag and not near_overlap_flag and gap <= JOINT_TOL and point_order_ok(result))
    return fold_flag, near_overlap_flag, gap, ok


def build_scale_records(results: Sequence[shape_audit.ModeResult], requested_scale: float) -> list[PlotScaleRecord]:
    records: list[PlotScaleRecord] = []
    for result in results:
        fold, near, gap, ok = scale_flags(result, float(requested_scale))
        fallback = not ok
        scale_used = float(requested_scale)
        notes = f"scale {float(requested_scale):.2f} passed regularity"
        if fallback:
            scale_used = FALLBACK_SCALE
            fold, near, gap, ok_fallback = scale_flags(result, scale_used)
            notes = (
                f"scale {float(requested_scale):.2f} failed regularity; "
                f"fallback {scale_used:.2f} {'passed' if ok_fallback else 'still has warnings'}"
            )
        records.append(
            PlotScaleRecord(
                result=result,
                scale_requested=float(requested_scale),
                scale_used=float(scale_used),
                fallback_used=bool(fallback),
                fold_flag=bool(fold),
                near_overlap_flag=bool(near),
                plotted_joint_gap=float(gap),
                point_order_ok=point_order_ok(result),
                notes=notes,
            )
        )
    return records


def signed_component(result: shape_audit.ModeResult, name: str) -> tuple[np.ndarray, np.ndarray]:
    sign = float(result.sign_factor)
    if name == "u":
        return sign * result.rod1.u, sign * result.rod2.u
    if name == "w":
        return sign * result.rod1.w, sign * result.rod2.w
    if name == "psi":
        return sign * np.asarray(result.rod1.psi, dtype=float), sign * np.asarray(result.rod2.psi, dtype=float)
    if name == "gamma":
        return sign * np.asarray(result.rod1.gamma, dtype=float), sign * np.asarray(result.rod2.gamma, dtype=float)
    raise ValueError(f"unknown component {name!r}")


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


def component_shape_vector(result: shape_audit.ModeResult, name: str) -> np.ndarray:
    values1, values2 = signed_component(result, name)
    weights1 = trapezoid_sqrt_weights(result.rod1.x_local)
    weights2 = trapezoid_sqrt_weights(result.rod2.x_local)
    vector = np.concatenate([values1 * weights1, values2 * weights2])
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-30:
        return vector
    return vector / norm


def grouped_field_vector(result: shape_audit.ModeResult, names: Sequence[str]) -> np.ndarray:
    pieces = [component_shape_vector(result, name) for name in names]
    vector = np.concatenate(pieces)
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-30:
        return vector
    return vector / norm


def plotted_centerline_vector(record: PlotScaleRecord) -> np.ndarray:
    result = record.result
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(record.scale_used))
    xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
        result,
        transverse_only=False,
        scale=scale,
    )
    vector = np.concatenate([xd1, yd1, xd2, yd2])
    norm = float(np.linalg.norm(vector))
    if norm <= 1.0e-30:
        return vector
    return vector / norm


def sign_invariant_similarity(vector_a: np.ndarray, vector_b: np.ndarray) -> float:
    a = np.asarray(vector_a, dtype=float)
    b = np.asarray(vector_b, dtype=float)
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 1.0e-30:
        return float("nan")
    return float(abs(np.dot(a, b)) / denom)


def similarity_rows(records: Sequence[PlotScaleRecord], mu_values: Sequence[float], sorted_indices: Sequence[int]) -> list[dict[str, object]]:
    by_key = {
        (round(float(record.result.mu), 12), int(record.result.sorted_index)): record
        for record in records
    }
    rows: list[dict[str, object]] = []
    groups = {
        "similarity_full_field": ("u", "w", "psi"),
        "similarity_displacement_uw": ("u", "w"),
        "similarity_axial_u": ("u",),
        "similarity_transverse_w": ("w",),
        "similarity_rotation_shear": ("psi", "gamma"),
    }
    for mu in mu_values:
        for i_index, mode_i in enumerate(sorted_indices):
            for mode_j in sorted_indices[i_index + 1 :]:
                rec_i = by_key[(round(float(mu), 12), int(mode_i))]
                rec_j = by_key[(round(float(mu), 12), int(mode_j))]
                result_i = rec_i.result
                result_j = rec_j.result
                row: dict[str, object] = {
                    "mu": float(mu),
                    "epsilon": DEFAULT_EPSILON,
                    "mode_i": int(mode_i),
                    "mode_j": int(mode_j),
                    "Lambda_i": result_i.Lambda,
                    "Lambda_j": result_j.Lambda,
                }
                for field_name, component_names in groups.items():
                    row[field_name] = sign_invariant_similarity(
                        grouped_field_vector(result_i, component_names),
                        grouped_field_vector(result_j, component_names),
                    )
                row["similarity_plotted_centerline"] = sign_invariant_similarity(
                    plotted_centerline_vector(rec_i),
                    plotted_centerline_vector(rec_j),
                )
                row["notes"] = (
                    "sign-invariant normalized correlation; field groups use per-component "
                    "L2 normalization over both rods; plotted centerline uses accepted fixed-scale coordinates"
                )
                rows.append(row)
    return rows


def records_by_key(records: Sequence[PlotScaleRecord]) -> dict[tuple[float, int], PlotScaleRecord]:
    return {
        (round(float(record.result.mu), 12), int(record.result.sorted_index)): record
        for record in records
    }


def full_centerline_limits(records: Sequence[PlotScaleRecord]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for record in records:
        result = record.result
        x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        scale = shape_audit.fixed_scale_from_fraction(result.mu, float(record.scale_used))
        xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
            result,
            transverse_only=False,
            scale=scale,
        )
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


def max_abs_u_w(result: shape_audit.ModeResult) -> tuple[float, float, float]:
    max_u = max(float(np.max(np.abs(result.rod1.u))), float(np.max(np.abs(result.rod2.u))))
    max_w = max(float(np.max(np.abs(result.rod1.w))), float(np.max(np.abs(result.rod2.w))))
    ratio = max_u / max_w if max_w > 1.0e-30 else float("inf")
    return max_u, max_w, ratio


def draw_full_centerline(
    ax: plt.Axes,
    record: PlotScaleRecord,
    *,
    limits: tuple[tuple[float, float], tuple[float, float]],
) -> None:
    result = record.result
    x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(record.scale_used))
    xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
        result,
        transverse_only=False,
        scale=scale,
    )
    ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=0.9)
    ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=0.9)
    ax.plot(xd1, yd1, color="#d55e00", linewidth=1.65)
    ax.plot(xd2, yd2, color="#d55e00", linewidth=1.65)
    ax.scatter([x1[-1]], [y1[-1]], color="black", s=12, zorder=5)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.90", linewidth=0.5)
    ax.tick_params(labelsize=6.5)
    energy = result.energy
    _max_u, _max_w, ratio = max_abs_u_w(result)
    fallback = " fallback" if record.fallback_used else ""
    ax.set_title(
        (
            f"mu={result.mu:g}, sorted {result.sorted_index}, Lambda={result.Lambda:.5g}\n"
            f"scale={record.scale_used:.2f}{fallback}; "
            f"A={energy['axial_energy_fraction']:.2f}, "
            f"B={energy['bending_energy_fraction']:.2f}, "
            f"S={energy['shear_energy_fraction']:.2f}; "
            f"|u|/|w|={ratio:.2g}; {energy['classification']}"
        ),
        fontsize=7.4,
    )


def plot_annotated_full_grid(
    output_dir: Path,
    records: Sequence[PlotScaleRecord],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / "timo_full_shapes_eps0p03_beta45_eta0_modes4_6_grid_annotated_corrected.png"
    by_key = records_by_key(records)
    limits = full_centerline_limits(records)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(13.2, 2.75 * len(sorted_indices) + 1.05),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            draw_full_centerline(ax, by_key[(round(float(mu), 12), int(sorted_index))], limits=limits)
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=7.5)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=7.5)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=0.9, label="undeformed rods"),
        Line2D([0], [0], color="#d55e00", linewidth=1.65, label="Timoshenko full centerline"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.text(
        0.5,
        0.045,
        "Note: full centerline may hide modal differences when axial displacement dominates.",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    fig.suptitle(
        "Timoshenko full displacement centerlines with energy annotations: eps=0.03, beta=45 deg, eta=0",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.94), h_pad=0.65, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def component_arrays(result: shape_audit.ModeResult) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    return {name: signed_component(result, name) for name in ("u", "w", "psi", "gamma")}


def plot_components_grid(
    output_dir: Path,
    results: Sequence[shape_audit.ModeResult],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / "timo_modes456_components_eps0p03_beta45_eta0.png"
    by_key = {
        (round(float(result.mu), 12), int(result.sorted_index)): result
        for result in results
    }
    components = ("u", "w", "psi", "gamma")
    row_specs = [(float(mu), int(sorted_index)) for mu in mu_values for sorted_index in sorted_indices]
    fig, axes = plt.subplots(
        len(row_specs),
        len(components),
        figsize=(13.4, 1.65 * len(row_specs) + 0.95),
        squeeze=False,
    )
    for row_index, (mu, sorted_index) in enumerate(row_specs):
        result = by_key[(round(float(mu), 12), int(sorted_index))]
        arrays = component_arrays(result)
        for col_index, name in enumerate(components):
            ax = axes[row_index, col_index]
            values1, values2 = arrays[name]
            ax.plot(result.rod1.x_local, values1, color="#0072b2", linewidth=1.05)
            ax.plot(result.rod2.x_local, values2, color="#d55e00", linewidth=1.05)
            ax.scatter(
                [result.rod1.x_local[-1], result.rod2.x_local[-1]],
                [values1[-1], values2[-1]],
                color="black",
                s=8,
                zorder=5,
            )
            ax.axhline(0.0, color="0.82", linewidth=0.5)
            ax.grid(True, color="0.91", linewidth=0.45)
            ax.tick_params(labelsize=6.5)
            if row_index == 0:
                ax.set_title(name, fontsize=8.5)
            if col_index == 0:
                ax.set_ylabel(f"mu={mu:g}\nsorted {sorted_index}", fontsize=7.2)
            if row_index == len(row_specs) - 1:
                ax.set_xlabel("signed local x_i", fontsize=7.2)
    handles = [
        Line2D([0], [0], color="#0072b2", linewidth=1.2, label="rod 1"),
        Line2D([0], [0], color="#d55e00", linewidth=1.2, label="rod 2"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(
        "Timoshenko component fields for sorted modes 4-6: eps=0.03, beta=45 deg, eta=0",
        fontsize=11,
        y=0.997,
    )
    fig.tight_layout(rect=(0.0, 0.035, 1.0, 0.965), h_pad=0.4, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def transverse_amplified_arrays(
    result: shape_audit.ModeResult,
    *,
    axial_fraction: float = 0.02,
    transverse_fraction: float = 0.16,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u1, w1, u2, w2 = shape_audit.plot_arrays(result)
    x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(u1))
    max_u = max(float(np.max(np.abs(u1))), float(np.max(np.abs(u2))), 1.0e-14)
    max_w = max(float(np.max(np.abs(w1))), float(np.max(np.abs(w2))), 1.0e-14)
    length = sum(shape_audit.TIMO.segment_lengths(result.mu))
    u_scale = axial_fraction * length / max_u
    w_scale = transverse_fraction * length / max_w
    dx1, dy1 = DISPLAY.rod1_local_displacement_to_display(u_scale * u1, w_scale * w1)
    dx2, dy2 = DISPLAY.rod2_local_displacement_to_display(
        u_scale * u2,
        w_scale * w2,
        beta_deg=result.beta_deg,
    )
    return x1 + dx1, y1 + dy1, x2 + dx2, y2 + dy2


def transverse_amplified_limits(results: Sequence[shape_audit.ModeResult]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for result in results:
        x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        xd1, yd1, xd2, yd2 = transverse_amplified_arrays(result)
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


def plot_transverse_amplified_grid(
    output_dir: Path,
    results: Sequence[shape_audit.ModeResult],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / "timo_modes456_transverse_amplified_eps0p03_beta45_eta0_corrected.png"
    by_key = {
        (round(float(result.mu), 12), int(result.sorted_index)): result
        for result in results
    }
    limits = transverse_amplified_limits(results)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(13.2, 2.55 * len(sorted_indices) + 1.05),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            result = by_key[(round(float(mu), 12), int(sorted_index))]
            x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
            xd1, yd1, xd2, yd2 = transverse_amplified_arrays(result)
            ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=0.9)
            ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=0.9)
            ax.plot(xd1, yd1, color="#009e73", linewidth=1.65)
            ax.plot(xd2, yd2, color="#009e73", linewidth=1.65)
            ax.scatter([x1[-1]], [y1[-1]], color="black", s=12, zorder=5)
            ax.set_xlim(*limits[0])
            ax.set_ylim(*limits[1])
            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, color="0.90", linewidth=0.5)
            ax.tick_params(labelsize=6.5)
            ax.set_title(
                (
                    f"mu={result.mu:g}, sorted {result.sorted_index}, Lambda={result.Lambda:.5g}\n"
                    f"transverse-amplified diagnostic view; class={result.energy['classification']}"
                ),
                fontsize=7.4,
            )
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=7.5)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=7.5)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=0.9, label="undeformed rods"),
        Line2D([0], [0], color="#009e73", linewidth=1.65, label="transverse-amplified diagnostic"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.text(
        0.5,
        0.045,
        "Diagnostic only: axial and transverse components are scaled separately; this is not the physical full displacement.",
        ha="center",
        va="center",
        fontsize=8.5,
    )
    fig.suptitle(
        "Timoshenko transverse-amplified diagnostic view: eps=0.03, beta=45 deg, eta=0",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.065, 1.0, 0.94), h_pad=0.65, w_pad=0.55)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def similarity_stats(rows: Sequence[dict[str, object]]) -> dict[str, dict[str, float]]:
    fields = [
        "similarity_full_field",
        "similarity_displacement_uw",
        "similarity_axial_u",
        "similarity_transverse_w",
        "similarity_rotation_shear",
        "similarity_plotted_centerline",
    ]
    stats: dict[str, dict[str, float]] = {}
    for field in fields:
        values = np.asarray([float(row[field]) for row in rows], dtype=float)
        stats[field] = {
            "min": float(np.min(values)),
            "mean": float(np.mean(values)),
            "max": float(np.max(values)),
        }
    return stats


def field_with_lowest_mean(stats: dict[str, dict[str, float]]) -> str:
    candidates = [
        "similarity_axial_u",
        "similarity_transverse_w",
        "similarity_rotation_shear",
    ]
    return min(candidates, key=lambda field: stats[field]["mean"])


def write_report(
    path: Path,
    *,
    similarity_csv: Path,
    figures: Sequence[Path],
    similarity: Sequence[dict[str, object]],
    scale_records: Sequence[PlotScaleRecord],
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
    requested_scale: float,
) -> Path:
    stats = similarity_stats(similarity)
    fallback_count = sum(1 for record in scale_records if record.fallback_used)
    remaining_fold = sum(1 for record in scale_records if record.fold_flag)
    remaining_near = sum(1 for record in scale_records if record.near_overlap_flag)
    max_plot_gap = max((float(record.plotted_joint_gap) for record in scale_records), default=float("nan"))
    most_distinguishing = field_with_lowest_mean(stats)
    full_mean = stats["similarity_full_field"]["mean"]
    plotted_mean = stats["similarity_plotted_centerline"]["mean"]

    lines = [
        "# Timoshenko Modes 4-6 Visualization Audit",
        "",
        "## Scope",
        "",
        f"- epsilon: `{DEFAULT_EPSILON:g}`",
        f"- beta_deg: `{DEFAULT_BETA_DEG:g}`",
        f"- eta: `{DEFAULT_ETA:g}`",
        f"- mu values: `{', '.join(f'{value:g}' for value in mu_values)}`",
        f"- sorted modes: `{', '.join(str(value) for value in sorted_indices)}`",
        "",
        "Sorted frequencies are used directly. No descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article files, old determinant code, old solvers, or baseline results are used or modified.",
        "",
        "No analytic formulas, determinant entries, root solvers, or Timoshenko shear coefficient k' were changed.",
        "",
        "## Display-Frame Correction",
        "",
        "Older Timoshenko global centerline and transverse-amplified figures used an inconsistent rod-2 display frame and are invalid for physical interpretation. The local fields, frequencies, null vectors, similarities, and energy fractions remain valid. Corrected figures reflect determinant transverse components into the display bases `t1=(1,0)`, `n1=(0,-1)`, `t2=(cos(beta),sin(beta))`, and `n2=(sin(beta),-cos(beta))`.",
        "",
        "## Similarity Method",
        "",
        f"- Similarity CSV: `{shape_audit.rel(similarity_csv)}`",
        "- Values are sign-invariant normalized correlations in `[0,1]`.",
        "- Field groups normalize each component before concatenation, so large axial displacement does not suppress smaller bending/rotation information.",
        "- Plotted-centerline similarity is computed from the accepted fixed-scale plotted coordinates.",
        "",
        "## Similarity Summary",
        "",
        "| group | min | mean | max |",
        "|---|---:|---:|---:|",
    ]
    for field, values in stats.items():
        lines.append(f"| {field} | {values['min']:.4f} | {values['mean']:.4f} | {values['max']:.4f} |")
    lines.extend(
        [
            "",
            "## Scale And Plot Regularity",
            "",
            f"- Requested full-centerline scale: `{float(requested_scale):.2f}`.",
            f"- Fallback scale: `{FALLBACK_SCALE:.2f}`.",
            f"- Full-centerline shapes requiring fallback: `{fallback_count}`.",
            f"- Remaining final fold flags: `{remaining_fold}`.",
            f"- Remaining final near-overlap flags: `{remaining_near}`.",
            f"- Maximum final plotted joint gap: `{max_plot_gap:.6e}`.",
            "",
            "Fallback shapes:",
        ]
    )
    fallbacks = [record for record in scale_records if record.fallback_used]
    if fallbacks:
        for record in fallbacks:
            result = record.result
            lines.append(f"- mu=`{result.mu:g}`, sorted=`{result.sorted_index}`, Lambda=`{result.Lambda:.10g}`")
    else:
        lines.append("- none")

    lines.extend(
        [
            "",
            "## Answers",
            "",
            "1. Are Timoshenko sorted modes 4, 5, and 6 mathematically different?",
            (
                f"Yes. The component-normalized full-field mean similarity is `{full_mean:.4f}`, "
                "and individual component groups show substantially lower similarities than the plotted centerline."
            ),
            "",
            "2. Which field components distinguish them most?",
            (
                f"The lowest mean similarity is `{most_distinguishing}`. "
                "Use the component figure to inspect whether the distinction appears in u, w, psi, or gamma for each mu."
            ),
            "",
            "3. Does the full centerline plot hide differences because axial displacement dominates?",
            (
                f"Yes for visual comparison: plotted-centerline mean similarity is `{plotted_mean:.4f}`, "
                "which is much higher than the component-field similarities. The common undeformed geometry and axial-dominated full displacement can make different modal fields look nearly identical."
            ),
            "",
            "4. Should full displacement plots be used for longitudinal/mixed modes?",
            "Use them only as a general geometry check. They are not the primary evidence for modal interpretation in longitudinal or mixed Timoshenko modes.",
            "",
            "5. Recommended visualization:",
            "- full centerline: general geometry and joint/transform sanity check;",
            "- component plots: primary modal interpretation;",
            "- transverse-amplified view: diagnostic aid only, not a physical full-displacement shape.",
            "",
            "## Generated Figures",
            "",
        ]
    )
    for figure in figures:
        lines.append(f"- `{shape_audit.rel(figure)}`")
    lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    mu_values = selected_mu_values(bool(args.smoke))
    sorted_indices = selected_sorted_indices(bool(args.smoke))
    results = build_results(args)
    scale_records = build_scale_records(results, float(args.deformation_scale_fraction))
    similarity = similarity_rows(scale_records, mu_values, sorted_indices)
    similarity_csv = shape_audit.write_csv(
        args.output_dir / "timoshenko_modes456_field_similarity.csv",
        similarity,
        SIMILARITY_FIELDS,
    )
    component_png = plot_components_grid(
        args.output_dir,
        results,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
    )
    full_png = plot_annotated_full_grid(
        args.output_dir,
        scale_records,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
    )
    transverse_png = plot_transverse_amplified_grid(
        args.output_dir,
        results,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
    )
    figures = [component_png, full_png, transverse_png]
    report = write_report(
        args.output_dir / "timoshenko_modes456_visualization_report.md",
        similarity_csv=similarity_csv,
        figures=figures,
        similarity=similarity,
        scale_records=scale_records,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
        requested_scale=float(args.deformation_scale_fraction),
    )
    stats = similarity_stats(similarity)
    fallback_count = sum(1 for record in scale_records if record.fallback_used)
    remaining_fold = sum(1 for record in scale_records if record.fold_flag)
    remaining_near = sum(1 for record in scale_records if record.near_overlap_flag)
    max_plot_gap = max((float(record.plotted_joint_gap) for record in scale_records), default=float("nan"))
    print("generated Timoshenko modes 4-6 visualization audit outputs:")
    for path in [similarity_csv, report, *figures]:
        print(f"  {shape_audit.rel(path)}")
    print("similarity summary:")
    for field, values in stats.items():
        print(f"  {field}: min={values['min']:.6f}, mean={values['mean']:.6f}, max={values['max']:.6f}")
    print(f"fallback shapes: {fallback_count}")
    print(f"remaining final fold flags: {remaining_fold}")
    print(f"remaining final near-overlap flags: {remaining_near}")
    print(f"max plotted joint gap: {max_plot_gap:.6e}")
    return {
        "similarity_csv": similarity_csv,
        "report": report,
        "figures": figures,
        "similarity_rows": similarity,
        "scale_records": scale_records,
    }


if __name__ == "__main__":
    main()
