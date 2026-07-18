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
MODEL_PREFIX = {MODEL_EB: "eb", MODEL_TIMO: "timo"}
MODEL_LABEL = {MODEL_EB: "Euler-Bernoulli", MODEL_TIMO: "Timoshenko"}
MODEL_COLOR = {MODEL_EB: "#1f77b4", MODEL_TIMO: "#d55e00"}

DEFAULT_EPSILON = 0.03
DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
DEFAULT_MU_VALUES = (0.0, 0.2, 0.4, 0.6)
DEFAULT_SORTED_INDICES = (4, 5, 6)
SMOKE_MU_VALUES = (0.0, 0.6)
SMOKE_SORTED_INDICES = (4,)
DEFAULT_N_POINTS = 801
SMOKE_N_POINTS = 301
DEFAULT_DEFORMATION_SCALE_FRACTION = 0.08
FALLBACK_DEFORMATION_SCALE_FRACTION = 0.05
DEFAULT_OUTPUT_DIR = (
    Path("results") / "eb_timo_mode_shapes_eps0p03_beta45_eta0_modes4_6"
)
SMOKE_OUTPUT_DIR = (
    Path("results") / "_smoke" / "eb_timo_mode_shapes_eps0p03_beta45_eta0_modes4_6"
)
JOINT_TOL = 1.0e-6

SUMMARY_FIELDS = [
    "model",
    "epsilon",
    "beta_deg",
    "eta",
    "mu",
    "sorted_index",
    "Lambda",
    "scale_requested",
    "scale_used",
    "fallback_used",
    "classification",
    "axial_energy_fraction",
    "bending_energy_fraction",
    "shear_energy_fraction",
    "plotted_joint_gap",
    "max_abs_kinematic_gap",
    "fold_flag",
    "near_overlap_flag",
    "point_order_ok",
    "figure_file",
    "notes",
]

REGULARITY_FIELDS = [
    "model",
    "epsilon",
    "beta_deg",
    "eta",
    "mu",
    "sorted_index",
    "Lambda",
    "rod_id",
    "scale_role",
    "scale_fraction",
    "deformation_scale",
    "min_abs_dr_plot_ds",
    "max_abs_dr_plot_ds",
    "projected_x_sign_change_count",
    "projected_x_is_monotone",
    "near_overlap_pair_count",
    "min_nonneighbor_distance",
    "fold_flag",
    "warning",
    "notes",
]


@dataclass(frozen=True)
class ShapeRecord:
    result: shape_audit.ModeResult
    scale_requested: float
    scale_used: float
    fallback_used: bool
    fold_flag: bool
    near_overlap_flag: bool
    point_order_ok: bool
    plotted_joint_gap: float
    max_abs_kinematic_gap: float
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


def mu_token(value: float) -> str:
    return value_token(float(value), min_decimals=1, max_decimals=3)


def eps_token(value: float) -> str:
    return value_token(float(value), min_decimals=2, max_decimals=4)


def beta_token(value: float) -> str:
    return value_token(float(value), min_decimals=0, max_decimals=2)


def eta_token(value: float) -> str:
    return value_token(float(value), min_decimals=0, max_decimals=3)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Diagnostic-only full EB/Timoshenko mode shapes for eps=0.03, beta=45 deg, eta=0.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument(
        "--deformation-scale-fraction",
        type=float,
        default=DEFAULT_DEFORMATION_SCALE_FRACTION,
    )
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
        eb_roots, timo_roots, root_warnings = provider.roots(DEFAULT_EPSILON, float(mu))
        for sorted_index in sorted_indices:
            root_index = int(sorted_index) - 1
            eb_lambda = float(eb_roots[root_index])
            timo_lambda = float(timo_roots[root_index])
            results.append(
                shape_audit.eb_mode_result(
                    epsilon=DEFAULT_EPSILON,
                    beta_deg=DEFAULT_BETA_DEG,
                    eta=DEFAULT_ETA,
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=eb_lambda,
                    n_points=int(args.n_points),
                    case_role="diagnostic",
                    root_source="sorted_global_solver",
                    root_warnings=root_warnings,
                )
            )
            results.append(
                shape_audit.timo_mode_result(
                    epsilon=DEFAULT_EPSILON,
                    beta_deg=DEFAULT_BETA_DEG,
                    eta=DEFAULT_ETA,
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=timo_lambda,
                    n_points=int(args.n_points),
                    case_role="diagnostic",
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
    scale_role: str,
) -> list[dict[str, object]]:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(scale_fraction))
    rows: list[dict[str, object]] = []
    for rod_id in (1, 2):
        row = shape_audit.regularity_row(
            result,
            rod_id=rod_id,
            scale_kind="fixed_fraction",
            deformation_scale_fraction=float(scale_fraction),
            deformation_scale_value=scale,
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
                "rod_id": int(rod_id),
                "scale_role": scale_role,
                "scale_fraction": float(scale_fraction),
                "deformation_scale": scale,
                "min_abs_dr_plot_ds": row["min_abs_dr_plot_ds"],
                "max_abs_dr_plot_ds": row["max_abs_dr_plot_ds"],
                "projected_x_sign_change_count": row["projected_x_sign_change_count"],
                "projected_x_is_monotone": row["projected_x_is_monotone"],
                "near_overlap_pair_count": row["near_overlap_pair_count"],
                "min_nonneighbor_distance": row["min_nonneighbor_distance"],
                "fold_flag": row["visually_folded_at_scale"],
                "warning": row["warning"],
                "notes": row["notes"],
            }
        )
    return rows


def scale_check(
    result: shape_audit.ModeResult,
    *,
    scale_fraction: float,
    scale_role: str,
) -> tuple[bool, bool, bool, float, float, list[dict[str, object]]]:
    rows = regularity_rows_for_scale(result, scale_fraction=scale_fraction, scale_role=scale_role)
    fold_flag = any(bool(row["fold_flag"]) for row in rows)
    near_overlap_flag = any(int(row["near_overlap_pair_count"]) > 0 for row in rows)
    order_ok = point_order_ok(result)
    plot_gap = plotted_joint_gap(result, scale_fraction)
    kin_gap = max_abs_kinematic_gap(result)
    passed = bool(
        not fold_flag
        and not near_overlap_flag
        and order_ok
        and plot_gap <= JOINT_TOL
        and (not isfinite(kin_gap) or kin_gap <= JOINT_TOL)
    )
    return passed, fold_flag, near_overlap_flag, plot_gap, kin_gap, rows


def individual_png_path(output_dir: Path, result: shape_audit.ModeResult) -> Path:
    prefix = MODEL_PREFIX[result.model]
    return output_dir / (
        f"{prefix}_full_shape_mu{mu_token(result.mu)}_sorted{result.sorted_index}"
        f"_eps{eps_token(result.epsilon)}_beta{beta_token(result.beta_deg)}"
        f"_eta{eta_token(result.eta)}_corrected.png"
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
        xd1, yd1, xd2, yd2 = deformed_arrays(result, scale_fraction=record.scale_used)
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


def draw_shape(
    ax: plt.Axes,
    record: ShapeRecord,
    *,
    limits: tuple[tuple[float, float], tuple[float, float]],
    title: bool,
) -> None:
    result = record.result
    x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
    xd1, yd1, xd2, yd2 = deformed_arrays(result, scale_fraction=record.scale_used)
    color = MODEL_COLOR[result.model]
    ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=0.95)
    ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=0.95)
    ax.plot(xd1, yd1, color=color, linewidth=1.75)
    ax.plot(xd2, yd2, color=color, linewidth=1.75)
    ax.scatter([x1[-1]], [y1[-1]], color="black", s=13, zorder=5)
    ax.set_xlim(*limits[0])
    ax.set_ylim(*limits[1])
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.90", linewidth=0.55)
    ax.tick_params(labelsize=7)
    if title:
        fallback = " fallback" if record.fallback_used else ""
        energy = result.energy
        classification = str(energy["classification"]).replace("_", " ")
        ax.set_title(
            (
                f"{MODEL_LABEL[result.model]} | mu={result.mu:g} | sorted {result.sorted_index}\n"
                f"Lambda={result.Lambda:.6g} | scale={record.scale_used:.2f}{fallback}\n"
                f"{classification}"
            ),
            fontsize=7.4,
        )


def write_individual_figure(record: ShapeRecord) -> None:
    record.figure_file.parent.mkdir(parents=True, exist_ok=True)
    limits = axis_limits([record])
    fig, ax = plt.subplots(figsize=(6.2, 4.1))
    draw_shape(ax, record, limits=limits, title=True)
    ax.set_xlabel("global x")
    ax.set_ylabel("global y")
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=0.95, label="undeformed rods"),
        Line2D([0], [0], color=MODEL_COLOR[record.result.model], linewidth=1.75, label="full centerline"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    ax.legend(handles=handles, loc="best", fontsize=7.5, frameon=False)
    fig.tight_layout()
    fig.savefig(record.figure_file, dpi=220, bbox_inches="tight")
    plt.close(fig)


def build_records(
    results: Sequence[shape_audit.ModeResult],
    *,
    output_dir: Path,
    scale_requested: float,
) -> tuple[list[ShapeRecord], list[dict[str, object]]]:
    records: list[ShapeRecord] = []
    regularity_rows: list[dict[str, object]] = []
    for result in results:
        passed, fold, near, plot_gap, kin_gap, rows = scale_check(
            result,
            scale_fraction=float(scale_requested),
            scale_role="requested",
        )
        regularity_rows.extend(rows)
        fallback_used = not passed
        scale_used = float(scale_requested)
        final_fold = fold
        final_near = near
        final_plot_gap = plot_gap
        final_kin_gap = kin_gap
        note_parts = [
            "full displacement centerline from u and w",
            "sorted frequency index used directly",
        ]
        if result.model == MODEL_TIMO:
            note_parts.append("corrected reflected display bases t2=(c,s), n2=(s,-c)")
        if fallback_used:
            scale_used = FALLBACK_DEFORMATION_SCALE_FRACTION
            passed_fb, fold_fb, near_fb, plot_gap_fb, kin_gap_fb, rows_fb = scale_check(
                result,
                scale_fraction=scale_used,
                scale_role="fallback",
            )
            regularity_rows.extend(rows_fb)
            final_fold = fold_fb
            final_near = near_fb
            final_plot_gap = plot_gap_fb
            final_kin_gap = kin_gap_fb
            note_parts.append(
                f"scale {scale_requested:.2f} failed regularity; "
                f"fallback {scale_used:.2f} {'passed' if passed_fb else 'still has warnings'}"
            )
        else:
            note_parts.append(f"scale {scale_requested:.2f} passed regularity")

        record = ShapeRecord(
            result=result,
            scale_requested=float(scale_requested),
            scale_used=float(scale_used),
            fallback_used=bool(fallback_used),
            fold_flag=bool(final_fold),
            near_overlap_flag=bool(final_near),
            point_order_ok=point_order_ok(result),
            plotted_joint_gap=float(final_plot_gap),
            max_abs_kinematic_gap=float(final_kin_gap),
            figure_file=individual_png_path(output_dir, result),
            notes="; ".join(note_parts),
        )
        records.append(record)
    return records, regularity_rows


def record_key(record: ShapeRecord) -> tuple[str, int, float]:
    return (record.result.model, int(record.result.sorted_index), round(float(record.result.mu), 12))


def write_model_grid(
    output_dir: Path,
    records: Sequence[ShapeRecord],
    *,
    model: str,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    model_records = [record for record in records if record.result.model == model]
    by_key = {(int(record.result.sorted_index), round(float(record.result.mu), 12)): record for record in model_records}
    output = output_dir / (
        f"{MODEL_PREFIX[model]}_full_shapes_eps{eps_token(DEFAULT_EPSILON)}"
        f"_beta{beta_token(DEFAULT_BETA_DEG)}_eta{eta_token(DEFAULT_ETA)}_modes4_6_grid_corrected.png"
    )
    limits = axis_limits(model_records)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(14.4, 3.0 * len(sorted_indices) + 0.8),
        squeeze=False,
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            record = by_key[(int(sorted_index), round(float(mu), 12))]
            draw_shape(ax, record, limits=limits, title=True)
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=8)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=0.95, label="undeformed rods"),
        Line2D([0], [0], color=MODEL_COLOR[model], linewidth=1.75, label="full centerline"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8)
    fig.suptitle(
        (
            f"{MODEL_LABEL[model]} full displacement shapes: eps={DEFAULT_EPSILON:g}, "
            f"beta=45 deg, eta=0, sorted modes 4-6"
        ),
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.04, 1.0, 0.94), h_pad=1.0, w_pad=0.8)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def write_combined_grid(
    output_dir: Path,
    records: Sequence[ShapeRecord],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    by_key = {record_key(record): record for record in records}
    output = output_dir / (
        f"eb_vs_timo_full_shapes_eps{eps_token(DEFAULT_EPSILON)}"
        f"_beta{beta_token(DEFAULT_BETA_DEG)}_eta{eta_token(DEFAULT_ETA)}_modes4_6_grid_corrected.png"
    )
    limits = axis_limits(records)
    row_specs = [(mode, model) for mode in sorted_indices for model in MODELS]
    fig, axes = plt.subplots(
        len(row_specs),
        len(mu_values),
        figsize=(14.4, 2.55 * len(row_specs) + 1.0),
        squeeze=False,
    )
    for row_index, (sorted_index, model) in enumerate(row_specs):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            record = by_key[(model, int(sorted_index), round(float(mu), 12))]
            draw_shape(ax, record, limits=limits, title=True)
            if row_index == len(row_specs) - 1:
                ax.set_xlabel("global x", fontsize=8)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=8)
    handles = [
        Line2D([0], [0], color="0.72", linestyle="--", linewidth=0.95, label="undeformed rods"),
        Line2D([0], [0], color=MODEL_COLOR[MODEL_EB], linewidth=1.75, label=MODEL_LABEL[MODEL_EB]),
        Line2D([0], [0], color=MODEL_COLOR[MODEL_TIMO], linewidth=1.75, label=MODEL_LABEL[MODEL_TIMO]),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle(
        "EB vs Timoshenko full displacement shapes: eps=0.03, beta=45 deg, eta=0",
        fontsize=11,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 0.96), h_pad=0.9, w_pad=0.8)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def summary_row(record: ShapeRecord) -> dict[str, object]:
    result = record.result
    energy = result.energy
    shear = float("nan") if result.model == MODEL_EB else energy["shear_energy_fraction"]
    return {
        "model": result.model,
        "epsilon": result.epsilon,
        "beta_deg": result.beta_deg,
        "eta": result.eta,
        "mu": result.mu,
        "sorted_index": result.sorted_index,
        "Lambda": result.Lambda,
        "scale_requested": record.scale_requested,
        "scale_used": record.scale_used,
        "fallback_used": record.fallback_used,
        "classification": energy["classification"],
        "axial_energy_fraction": energy["axial_energy_fraction"],
        "bending_energy_fraction": energy["bending_energy_fraction"],
        "shear_energy_fraction": shear,
        "plotted_joint_gap": record.plotted_joint_gap,
        "max_abs_kinematic_gap": record.max_abs_kinematic_gap,
        "fold_flag": record.fold_flag,
        "near_overlap_flag": record.near_overlap_flag,
        "point_order_ok": record.point_order_ok,
        "figure_file": shape_audit.rel(record.figure_file),
        "notes": record.notes,
    }


def write_report(
    path: Path,
    *,
    records: Sequence[ShapeRecord],
    regularity_csv: Path,
    summary_csv: Path,
    figure_paths: Sequence[Path],
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
    scale_requested: float,
) -> Path:
    passed_requested = [record for record in records if not record.fallback_used]
    fallback_records = [record for record in records if record.fallback_used]
    final_fold = [record for record in records if record.fold_flag]
    final_near = [record for record in records if record.near_overlap_flag]
    max_plot_gap = max((float(record.plotted_joint_gap) for record in records), default=float("nan"))

    lines = [
        "# EB/Timoshenko Full Mode Shape Diagnostic",
        "",
        "## Scope",
        "",
        f"- epsilon: `{DEFAULT_EPSILON:g}`",
        f"- beta_deg: `{DEFAULT_BETA_DEG:g}`",
        f"- eta: `{DEFAULT_ETA:g}`",
        f"- mu values: `{', '.join(f'{value:g}' for value in mu_values)}`",
        f"- sorted modes: `{', '.join(str(value) for value in sorted_indices)}`",
        "- Models: `Euler-Bernoulli`, `Timoshenko`",
        "",
        "The figures show full deformed centerlines from both longitudinal and transverse displacement components. They are not w-only plots and do not include component grids.",
        "",
        "The workflow uses sorted frequencies only. It does not use descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article workspaces, article figures, old determinant code, old solvers, or baseline-result edits.",
        "",
        "No analytic formulas, determinant entries, root solvers, or Timoshenko shear coefficient k' were changed.",
        "",
        "## Display-Frame Correction",
        "",
        "Older Timoshenko global centerline figures from this workflow used an inconsistent Cartesian interpretation of the determinant transverse components and must not be used for physical interpretation. Local fields, frequencies, null vectors, and energy fractions remain valid. The Timoshenko display mapping is `t1=(1,0)`, `n1=(0,-1)`, `t2=(cos(beta),sin(beta))`, and `n2=(sin(beta),-cos(beta))`, giving `dY1=-w1` and `dY2=sin(beta)*u2-cos(beta)*w2`. EB uses its own opposite determinant transverse sign, mapped through `n1=(0,1)` and `n2=(-sin(beta),cos(beta))`; both mappings produce the same physical display frame in the thin limit.",
        "",
        "## Scale Policy",
        "",
        f"- Requested fixed deformation scale fraction: `{scale_requested:.2f}`.",
        f"- Fallback scale fraction: `{FALLBACK_DEFORMATION_SCALE_FRACTION:.2f}`.",
        f"- Shapes passing at `{scale_requested:.2f}`: `{len(passed_requested)}` / `{len(records)}`.",
        f"- Shapes requiring fallback: `{len(fallback_records)}`.",
        f"- Final shapes with fold flags: `{len(final_fold)}`.",
        f"- Final shapes with near-overlap flags: `{len(final_near)}`.",
        f"- Maximum final plotted joint gap: `{max_plot_gap:.6e}`.",
        "",
    ]
    if fallback_records:
        lines.append("Fallback shapes:")
        for record in fallback_records:
            result = record.result
            lines.append(
                (
                    f"- `{MODEL_LABEL[result.model]}`, mu=`{result.mu:g}`, "
                    f"sorted=`{result.sorted_index}`, final scale=`{record.scale_used:.2f}`"
                )
            )
        lines.append("")
    else:
        lines.append("No fallback shapes were required.")
        lines.append("")

    if final_fold or final_near:
        lines.append("Remaining final regularity warnings:")
        for record in sorted(set(final_fold + final_near), key=lambda item: (item.result.model, item.result.sorted_index, item.result.mu)):
            result = record.result
            lines.append(
                (
                    f"- `{MODEL_LABEL[result.model]}`, mu=`{result.mu:g}`, sorted=`{result.sorted_index}`: "
                    f"fold={record.fold_flag}, near_overlap={record.near_overlap_flag}"
                )
            )
        lines.append("")
    else:
        lines.append("No plotted final shape has fold or near-overlap warnings.")
        lines.append("")

    lines.extend(
        [
            "Scale affects visualization only. It does not affect sorted frequencies, energy fractions, or mode classification.",
            "",
            "## Frequencies",
            "",
            "| model | mu | sorted_index | Lambda | classification | scale_used | fallback |",
            "|---|---:|---:|---:|---|---:|---|",
        ]
    )
    for record in sorted(records, key=lambda item: (item.result.model, item.result.mu, item.result.sorted_index)):
        result = record.result
        lines.append(
            (
                f"| {MODEL_LABEL[result.model]} | {result.mu:g} | {result.sorted_index} | "
                f"{result.Lambda:.10g} | {result.energy['classification']} | "
                f"{record.scale_used:.2f} | {str(record.fallback_used).lower()} |"
            )
        )

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- Summary CSV: `{shape_audit.rel(summary_csv)}`",
            f"- Regularity CSV: `{shape_audit.rel(regularity_csv)}`",
            "",
            "Generated figures:",
        ]
    )
    for figure in figure_paths:
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
    records, regularity_rows = build_records(
        results,
        output_dir=args.output_dir,
        scale_requested=float(args.deformation_scale_fraction),
    )

    for record in records:
        write_individual_figure(record)

    figure_paths: list[Path] = []
    eb_grid = write_model_grid(
        args.output_dir,
        records,
        model=MODEL_EB,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
    )
    timo_grid = write_model_grid(
        args.output_dir,
        records,
        model=MODEL_TIMO,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
    )
    figure_paths.extend([eb_grid, timo_grid])
    if not bool(args.smoke):
        figure_paths.append(
            write_combined_grid(
                args.output_dir,
                records,
                mu_values=mu_values,
                sorted_indices=sorted_indices,
            )
        )
    figure_paths.extend(record.figure_file for record in records)

    summary_rows = [summary_row(record) for record in records]
    summary_csv = shape_audit.write_csv(
        args.output_dir / "mode_shape_full_displacement_summary.csv",
        summary_rows,
        SUMMARY_FIELDS,
    )
    regularity_csv = shape_audit.write_csv(
        args.output_dir / "mode_shape_full_displacement_regularity_audit.csv",
        regularity_rows,
        REGULARITY_FIELDS,
    )
    report = write_report(
        args.output_dir / "mode_shape_full_displacement_report.md",
        records=records,
        regularity_csv=regularity_csv,
        summary_csv=summary_csv,
        figure_paths=figure_paths,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
        scale_requested=float(args.deformation_scale_fraction),
    )

    fallback_count = sum(1 for record in records if record.fallback_used)
    scale08_count = sum(
        1
        for record in records
        if abs(float(record.scale_used) - float(args.deformation_scale_fraction)) <= 1.0e-12
    )
    final_fold = sum(1 for record in records if record.fold_flag)
    final_near = sum(1 for record in records if record.near_overlap_flag)
    max_plot_gap = max((float(record.plotted_joint_gap) for record in records), default=float("nan"))

    print("generated EB/Timoshenko full mode shape diagnostic outputs:")
    for path in [report, summary_csv, regularity_csv, *figure_paths]:
        print(f"  {shape_audit.rel(path)}")
    print(f"shapes using requested scale {float(args.deformation_scale_fraction):.2f}: {scale08_count}")
    print(f"fallback shapes using {FALLBACK_DEFORMATION_SCALE_FRACTION:.2f}: {fallback_count}")
    print(f"final fold flags: {final_fold}")
    print(f"final near-overlap flags: {final_near}")
    print(f"max plotted joint gap: {max_plot_gap:.6e}")
    return {
        "report": report,
        "summary_csv": summary_csv,
        "regularity_csv": regularity_csv,
        "figure_paths": figure_paths,
        "records": records,
        "regularity_rows": regularity_rows,
    }


if __name__ == "__main__":
    main()
