from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    thickness_to_length_ratios,
)
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    analytic_null_vector,
    normalize_components,
)
from scripts.lib.thickness_mismatch_diagnostic_helpers import (  # noqa: E402
    mu_grid,
    roots_by_mu_eta,
    track_descendants_from_mu0,
    tracking_warning_rows,
)
from scripts.lib.thickness_mismatch_mac_tracking import reconstruct_components_eta  # noqa: E402


# =========================
# User-editable defaults
# =========================
DEFAULT_BRANCH_INDEX = 5
DEFAULT_BETA_DEG = 15.0
DEFAULT_EPSILON = 0.0025
DEFAULT_MU_VALUES = (0.0, 0.15)
DEFAULT_ETA_VALUES = (-0.5, 0.0, 0.5)

DEFAULT_MU_STEP = 0.005
DEFAULT_NUM_TRACKED_BRANCHES = 8
DEFAULT_NUM_SORTED_ROOTS = 12
DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
DEFAULT_NUM_SHAPE_SAMPLES = 401

DEFAULT_MAC_WARNING_THRESHOLD = 0.9
DEFAULT_MAC_MARGIN_WARNING_THRESHOLD = 0.05
DEFAULT_MAX_SORTED_POSITION_JUMP = 1
DEFAULT_THICKNESS_RATIO_LIMIT = 0.1

DEFAULT_L_TOTAL = 2.0
DEFAULT_MODE_SCALE = 0.16
DEFAULT_NORMALIZE = "max-full"
DEFAULT_DPI = 240
DEFAULT_FIGSIZE = (7.4, 4.6)
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"


COMPONENT_KEYS = ("u_left", "w_left", "u_right", "w_right")


@dataclass(frozen=True)
class TrackingCase:
    eta: float
    mu: float
    Lambda: float
    current_sorted_index: int
    tracking_step_status: str
    tracking_row: dict[str, float | int | str]


@dataclass(frozen=True)
class ShapeCase:
    eta: float
    mu: float
    Lambda: float
    current_sorted_index: int
    tracking_step_status: str
    components: dict[str, np.ndarray]
    matrix_residual_max: float
    smallest_singular_value: float
    singular_value_ratio: float
    thickness_ratio_1: float
    thickness_ratio_2: float


def number_token(value: float, *, force_one_decimal: bool = False) -> str:
    text = f"{float(value):.1f}" if force_one_decimal else f"{float(value):.10g}"
    return text.replace("-", "m").replace(".", "p")


def eta_token(value: float) -> str:
    value_f = float(value)
    if np.isclose(value_f, 0.0, rtol=0.0, atol=1e-12):
        return "0"
    prefix = "p" if value_f > 0.0 else "m"
    return prefix + number_token(abs(value_f))


def eta_list_token(values: Sequence[float]) -> str:
    return "_".join(eta_token(float(value)) for value in values)


def mu_token(value: float) -> str:
    value_f = float(value)
    return number_token(value_f, force_one_decimal=np.isclose(value_f, round(value_f), rtol=0.0, atol=1e-12))


def display_path(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot diagnostic thickness-mismatch analytic mode shapes for one "
            "descendant branch, overlaying several eta values at each target mu."
        )
    )
    parser.add_argument("--branch-index", type=int, default=DEFAULT_BRANCH_INDEX)
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mus", type=float, nargs="+", default=list(DEFAULT_MU_VALUES))
    parser.add_argument("--etas", type=float, nargs="+", default=list(DEFAULT_ETA_VALUES))
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--num-tracked-branches", type=int, default=DEFAULT_NUM_TRACKED_BRANCHES)
    parser.add_argument("--num-sorted-roots", type=int, default=DEFAULT_NUM_SORTED_ROOTS)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--root-lmax0", type=float, default=DEFAULT_ROOT_LMAX0)
    parser.add_argument("--num-shape-samples", type=int, default=DEFAULT_NUM_SHAPE_SAMPLES)
    parser.add_argument("--mac-warning-threshold", type=float, default=DEFAULT_MAC_WARNING_THRESHOLD)
    parser.add_argument("--mac-margin-warning-threshold", type=float, default=DEFAULT_MAC_MARGIN_WARNING_THRESHOLD)
    parser.add_argument("--max-sorted-position-jump", type=int, default=DEFAULT_MAX_SORTED_POSITION_JUMP)
    parser.add_argument("--thickness-ratio-limit", type=float, default=DEFAULT_THICKNESS_RATIO_LIMIT)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--mode-scale", type=float, default=DEFAULT_MODE_SCALE)
    parser.add_argument("--normalize", choices=("max-full", "max-transverse", "none"), default=DEFAULT_NORMALIZE)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def validate_args(args: argparse.Namespace) -> None:
    if int(args.branch_index) <= 0:
        raise ValueError("--branch-index must be positive.")
    if int(args.num_tracked_branches) < int(args.branch_index):
        raise ValueError("--num-tracked-branches must be at least --branch-index.")
    if int(args.num_sorted_roots) < int(args.num_tracked_branches):
        raise ValueError("--num-sorted-roots must be at least --num-tracked-branches.")
    if any(float(mu) < -1e-12 for mu in args.mus):
        raise ValueError("This diagnostic tracks from mu=0 and expects non-negative target mu values.")
    if not args.etas:
        raise ValueError("At least one eta value is required.")


def output_png_path(args: argparse.Namespace, mu: float) -> Path:
    return Path(args.output_dir) / (
        f"thickness_mismatch_shapes_branch{int(args.branch_index)}_"
        f"beta{number_token(float(args.beta_deg))}_"
        f"mu{mu_token(float(mu))}_"
        f"eta_{eta_list_token(args.etas)}.png"
    )


def output_report_path(args: argparse.Namespace, mu: float) -> Path:
    return Path(args.output_dir) / (
        f"thickness_mismatch_shapes_branch{int(args.branch_index)}_"
        f"beta{number_token(float(args.beta_deg))}_"
        f"mu{mu_token(float(mu))}_"
        f"eta_report.md"
    )


def tracking_grid_for_targets(target_mus: Sequence[float], mu_step_value: float) -> np.ndarray:
    max_mu = max(float(mu) for mu in target_mus)
    base = mu_grid(0.0, max_mu, float(mu_step_value))
    values = {round(float(mu), 12) for mu in base}
    values.update(round(float(mu), 12) for mu in target_mus)
    values.add(0.0)
    out = np.asarray(sorted(values), dtype=float)
    if not np.isclose(out[0], 0.0, rtol=0.0, atol=1e-12):
        raise RuntimeError("Tracking grid must start at mu=0.")
    return out


def find_tracking_row(
    rows: Sequence[dict[str, float | int | str]],
    *,
    mu: float,
    branch_index: int,
) -> dict[str, float | int | str]:
    matches = [
        row
        for row in rows
        if int(row["branch_index_from_mu0"]) == int(branch_index)
        and np.isclose(float(row["mu"]), float(mu), rtol=0.0, atol=1e-12)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Could not find tracking row for branch={branch_index}, mu={float(mu):g}.")
    return matches[0]


def compute_tracking_cases(
    args: argparse.Namespace,
    mu_values: np.ndarray,
) -> tuple[dict[float, dict[float, TrackingCase]], dict[float, list[dict[str, float | int | str]]]]:
    beta_rad = float(np.deg2rad(float(args.beta_deg)))
    cases: dict[float, dict[float, TrackingCase]] = {}
    rows_by_eta: dict[float, list[dict[str, float | int | str]]] = {}
    for eta in args.etas:
        roots = roots_by_mu_eta(
            beta_rad=beta_rad,
            epsilon=float(args.epsilon),
            eta=float(eta),
            mu_values=mu_values,
            n_roots=int(args.num_sorted_roots),
            root_lmax0=float(args.root_lmax0),
            root_scan_step=float(args.root_scan_step),
        )
        tracking = track_descendants_from_mu0(
            beta_rad=beta_rad,
            epsilon=float(args.epsilon),
            eta=float(eta),
            mu_values=mu_values,
            roots_by_mu=roots,
            num_descendants=int(args.num_tracked_branches),
            num_shape_samples=int(args.num_shape_samples),
            mac_warning_threshold=float(args.mac_warning_threshold),
            mac_margin_warning_threshold=float(args.mac_margin_warning_threshold),
            max_sorted_position_jump=int(args.max_sorted_position_jump),
        )
        rows_by_eta[float(eta)] = tracking.rows
        cases[float(eta)] = {}
        for target_mu in args.mus:
            row = find_tracking_row(
                tracking.rows,
                mu=float(target_mu),
                branch_index=int(args.branch_index),
            )
            cases[float(eta)][float(target_mu)] = TrackingCase(
                eta=float(eta),
                mu=float(target_mu),
                Lambda=float(row["Lambda_tracked"]),
                current_sorted_index=int(row["mac_sorted_root_index"]),
                tracking_step_status=str(row["tracking_step_status"]),
                tracking_row=row,
            )
    return cases, rows_by_eta


def build_shape_case(args: argparse.Namespace, tracking_case: TrackingCase, s_norm: np.ndarray) -> ShapeCase:
    beta_rad = float(np.deg2rad(float(args.beta_deg)))
    matrix = assemble_clamped_coupled_matrix_eta(
        tracking_case.Lambda,
        beta_rad,
        tracking_case.mu,
        float(args.epsilon),
        tracking_case.eta,
    )
    coeff, smallest, ratio = analytic_null_vector(matrix)
    components_raw = reconstruct_components_eta(
        tracking_case.Lambda,
        mu=tracking_case.mu,
        eta=tracking_case.eta,
        epsilon=float(args.epsilon),
        coeff=coeff,
        s_norm=s_norm,
    )
    components, _scale = normalize_components(
        components_raw,
        plot_kind="full",
        normalize=str(args.normalize),
    )
    matrix_residual = matrix @ coeff
    ratio1, ratio2 = thickness_to_length_ratios(float(args.epsilon), tracking_case.mu, tracking_case.eta)
    return ShapeCase(
        eta=tracking_case.eta,
        mu=tracking_case.mu,
        Lambda=tracking_case.Lambda,
        current_sorted_index=tracking_case.current_sorted_index,
        tracking_step_status=tracking_case.tracking_step_status,
        components=components,
        matrix_residual_max=float(np.max(np.abs(matrix_residual))),
        smallest_singular_value=float(smallest),
        singular_value_ratio=float(ratio),
        thickness_ratio_1=float(ratio1),
        thickness_ratio_2=float(ratio2),
    )


def component_vector(components: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate([np.asarray(components[key], dtype=float) for key in COMPONENT_KEYS])


def align_case_to_reference(case: ShapeCase, reference: ShapeCase) -> ShapeCase:
    vector = component_vector(case.components)
    reference_vector = component_vector(reference.components)
    if float(np.dot(vector, reference_vector)) >= 0.0:
        return case
    components = {key: -np.asarray(value, dtype=float) for key, value in case.components.items()}
    return ShapeCase(
        eta=case.eta,
        mu=case.mu,
        Lambda=case.Lambda,
        current_sorted_index=case.current_sorted_index,
        tracking_step_status=case.tracking_step_status,
        components=components,
        matrix_residual_max=case.matrix_residual_max,
        smallest_singular_value=case.smallest_singular_value,
        singular_value_ratio=case.singular_value_ratio,
        thickness_ratio_1=case.thickness_ratio_1,
        thickness_ratio_2=case.thickness_ratio_2,
    )


def align_cases(cases: Sequence[ShapeCase]) -> list[ShapeCase]:
    if not cases:
        return []
    reference = min(cases, key=lambda case: abs(float(case.eta)))
    return [align_case_to_reference(case, reference) for case in cases]


def arm_lengths(args: argparse.Namespace, mu: float) -> tuple[float, float]:
    ell = float(args.l_total) / 2.0
    return ell * (1.0 - float(mu)), ell * (1.0 + float(mu))


def base_coordinates(
    args: argparse.Namespace,
    s_norm: np.ndarray,
    *,
    beta_rad: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    left_length, right_length = arm_lengths(args, float(mu))
    xi = np.asarray(s_norm, dtype=float)
    x_left = left_length * xi
    y_left = np.zeros_like(x_left)
    x_right = left_length + right_length * xi * np.cos(beta_rad)
    y_right = right_length * xi * np.sin(beta_rad)
    return x_left, y_left, x_right, y_right


def deformed_coordinates(
    args: argparse.Namespace,
    components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    *,
    beta_rad: float,
    mu: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_left, y_left, x_right, y_right = base_coordinates(args, s_norm, beta_rad=beta_rad, mu=mu)
    tangent_right = np.array([np.cos(beta_rad), np.sin(beta_rad)], dtype=float)
    normal_right = np.array([-np.sin(beta_rad), np.cos(beta_rad)], dtype=float)
    u_left = np.asarray(components["u_left"], dtype=float)
    w_left = np.asarray(components["w_left"], dtype=float)
    u_right = np.asarray(components["u_right"], dtype=float)
    w_right = np.asarray(components["w_right"], dtype=float)
    x_left_def = x_left + float(args.mode_scale) * u_left
    y_left_def = y_left + float(args.mode_scale) * w_left
    x_right_def = x_right + float(args.mode_scale) * (u_right * tangent_right[0] + w_right * normal_right[0])
    y_right_def = y_right + float(args.mode_scale) * (u_right * tangent_right[1] + w_right * normal_right[1])
    return x_left_def, y_left_def, x_right_def, y_right_def


def shared_axis_limits(
    args: argparse.Namespace,
    cases: Sequence[ShapeCase],
    s_norm: np.ndarray,
    *,
    beta_rad: float,
    mu: float,
) -> tuple[tuple[float, float], tuple[float, float]]:
    x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(args, s_norm, beta_rad=beta_rad, mu=mu)
    x_values = [x_left_base, x_right_base]
    y_values = [y_left_base, y_right_base]
    for case in cases:
        x_left, y_left, x_right, y_right = deformed_coordinates(
            args,
            case.components,
            s_norm,
            beta_rad=beta_rad,
            mu=mu,
        )
        x_values.extend([x_left, x_right])
        y_values.extend([y_left, y_right])
    x_all = np.concatenate([np.asarray(value, dtype=float) for value in x_values])
    y_all = np.concatenate([np.asarray(value, dtype=float) for value in y_values])
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.08 * x_span), float(np.max(x_all) + 0.08 * x_span)),
        (float(np.min(y_all) - 0.16 * y_span), float(np.max(y_all) + 0.16 * y_span)),
    )


def plot_cases_for_mu(args: argparse.Namespace, mu: float, cases: Sequence[ShapeCase], s_norm: np.ndarray) -> Path:
    output = output_png_path(args, float(mu))
    output.parent.mkdir(parents=True, exist_ok=True)
    beta_rad = float(np.deg2rad(float(args.beta_deg)))
    x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(args, s_norm, beta_rad=beta_rad, mu=mu)
    x_limits, y_limits = shared_axis_limits(args, cases, s_norm, beta_rad=beta_rad, mu=mu)
    colors = {-0.5: "#1f77b4", 0.0: "#2f2f2f", 0.5: "#d62728"}
    fallback_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    ax.plot(x_left_base, y_left_base, color="0.72", linestyle="--", linewidth=1.1, label="undeformed")
    ax.plot(x_right_base, y_right_base, color="0.72", linestyle="--", linewidth=1.1)
    for idx, case in enumerate(cases):
        color = colors.get(round(float(case.eta), 10), fallback_colors[idx % len(fallback_colors)])
        x_left, y_left, x_right, y_right = deformed_coordinates(
            args,
            case.components,
            s_norm,
            beta_rad=beta_rad,
            mu=mu,
        )
        label = rf"$\eta={case.eta:g}$, $\Lambda={case.Lambda:.5g}$, sorted {case.current_sorted_index}"
        ax.plot(x_left, y_left, color=color, linewidth=2.0, label=label)
        ax.plot(x_right, y_right, color=color, linewidth=2.0)
    ax.scatter([x_left_base[-1]], [y_left_base[-1]], color="black", s=14, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", fontsize=8, frameon=False)
    ax.set_title(
        rf"Thickness mismatch descendant {int(args.branch_index)}: "
        rf"$\beta={float(args.beta_deg):g}^\circ$, "
        rf"$\epsilon={float(args.epsilon):g}$, $\mu={float(mu):g}$",
        fontsize=11,
    )
    fig.tight_layout()
    fig.savefig(output, dpi=int(args.dpi), bbox_inches="tight")
    plt.close(fig)
    return output


def table_lines(rows: Sequence[Sequence[str]]) -> list[str]:
    widths = [max(len(row[col]) for row in rows) for col in range(len(rows[0]))]
    lines: list[str] = []
    for idx, row in enumerate(rows):
        lines.append("| " + " | ".join(value.ljust(widths[col]) for col, value in enumerate(row)) + " |")
        if idx == 0:
            lines.append("| " + " | ".join("-" * widths[col] for col in range(len(row))) + " |")
    return lines


def write_report(
    args: argparse.Namespace,
    *,
    mu: float,
    cases: Sequence[ShapeCase],
    tracking_rows_by_eta: dict[float, list[dict[str, float | int | str]]],
    output_png: Path,
) -> Path:
    report_path = output_report_path(args, float(mu))
    report_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# Thickness-Mismatch Branch {int(args.branch_index)} Shapes vs Eta",
        "",
        "## Parameters",
        "",
        f"- beta: {float(args.beta_deg):g} deg",
        f"- epsilon: {float(args.epsilon):g}",
        f"- mu: {float(mu):g}",
        f"- eta values: {', '.join(f'{float(value):g}' for value in args.etas)}",
        f"- branch identity: descendant {int(args.branch_index)} seeded at `mu=0` for each fixed eta",
        f"- tracking: adjacent-step analytic shape MAC via `scripts/lib/thickness_mismatch_mac_tracking.py`",
        f"- PNG: `{display_path(output_png)}`",
        "",
        "Sorted position is diagnostic metadata only; it does not define the",
        "descendant branch identity.",
        "",
        "## Shape Cases",
        "",
    ]
    rows = [
        [
            "eta",
            "Lambda",
            "current sorted position",
            "tracking status",
            "2*r1/l1",
            "2*r2/l2",
            "matrix residual",
        ]
    ]
    for case in cases:
        rows.append(
            [
                f"{case.eta:g}",
                f"{case.Lambda:.10g}",
                str(case.current_sorted_index),
                case.tracking_step_status,
                f"{case.thickness_ratio_1:.6g}",
                f"{case.thickness_ratio_2:.6g}",
                f"{case.matrix_residual_max:.3e}",
            ]
        )
    lines.extend(table_lines(rows))

    lines.extend(["", "## Thin-Rod Criterion", ""])
    lines.append(f"- criterion: `2*r_i/l_i <= {float(args.thickness_ratio_limit):g}` for both rods")
    violations = [
        case
        for case in cases
        if case.thickness_ratio_1 > float(args.thickness_ratio_limit)
        or case.thickness_ratio_2 > float(args.thickness_ratio_limit)
    ]
    if violations:
        lines.append("- WARNING: criterion violations were found:")
        for case in violations:
            bad_rods = []
            if case.thickness_ratio_1 > float(args.thickness_ratio_limit):
                bad_rods.append("1")
            if case.thickness_ratio_2 > float(args.thickness_ratio_limit):
                bad_rods.append("2")
            lines.append(
                f"  - eta={case.eta:g}: rod(s) {', '.join(bad_rods)}, "
                f"2*r1/l1={case.thickness_ratio_1:.6g}, "
                f"2*r2/l2={case.thickness_ratio_2:.6g}"
            )
    else:
        lines.append("- no violations at this target mu.")

    lines.extend(["", "## Tracking Warnings Up To This Mu", ""])
    any_tracking_warning = False
    for eta in args.etas:
        rows_for_eta = [
            row
            for row in tracking_rows_by_eta[float(eta)]
            if float(row["mu"]) <= float(mu) + 1e-12
            and int(row["branch_index_from_mu0"]) == int(args.branch_index)
        ]
        warnings = tracking_warning_rows(rows_for_eta)
        if warnings:
            any_tracking_warning = True
            lines.append(f"- eta={float(eta):g}: {len(warnings)} warning row(s)")
            for row in warnings[:10]:
                lines.append(
                    f"  - mu={float(row['mu']):g}, accepted sorted position={int(row['mac_sorted_root_index'])}, "
                    f"candidate position={int(row['diagnostic_candidate_sorted_position'])}, "
                    f"candidate MAC={float(row['diagnostic_candidate_mac_to_previous']):.6g}, "
                    f"status={row['tracking_step_status']}"
                )
            if len(warnings) > 10:
                lines.append(f"  - plus {len(warnings) - 10} additional rows.")
        else:
            lines.append(f"- eta={float(eta):g}: no warnings for descendant {int(args.branch_index)} up to this mu.")
    if not any_tracking_warning:
        lines.append("")
        lines.append("No tracking warnings were found for the plotted descendant up to this target mu.")

    lines.extend(
        [
            "",
            "This diagnostic writes only results files. It does not modify article",
            "files, article figures, the old determinant, old solvers, or the FEM",
            "physical model.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def build_cases_by_mu(
    args: argparse.Namespace,
    tracking_cases: dict[float, dict[float, TrackingCase]],
) -> dict[float, list[ShapeCase]]:
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)
    out: dict[float, list[ShapeCase]] = {}
    for mu in args.mus:
        cases = [build_shape_case(args, tracking_cases[float(eta)][float(mu)], s_norm) for eta in args.etas]
        out[float(mu)] = align_cases(cases)
    return out


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    validate_args(args)
    args.output_dir = Path(args.output_dir)
    target_mus = [float(mu) for mu in args.mus]
    eta_values = [float(eta) for eta in args.etas]
    args.mus = target_mus
    args.etas = eta_values

    tracking_mu_values = tracking_grid_for_targets(target_mus, float(args.mu_step))
    tracking_cases, tracking_rows_by_eta = compute_tracking_cases(args, tracking_mu_values)
    shape_cases_by_mu = build_cases_by_mu(args, tracking_cases)
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)

    outputs: list[Path] = []
    reports: list[Path] = []
    print("Thickness-mismatch branch-shape diagnostic")
    print(
        f"branch={int(args.branch_index)}, beta={float(args.beta_deg):g} deg, "
        f"epsilon={float(args.epsilon):g}, etas={', '.join(f'{eta:g}' for eta in eta_values)}"
    )
    for mu in target_mus:
        cases = shape_cases_by_mu[float(mu)]
        output_png = plot_cases_for_mu(args, float(mu), cases, s_norm)
        report = write_report(
            args,
            mu=float(mu),
            cases=cases,
            tracking_rows_by_eta=tracking_rows_by_eta,
            output_png=output_png,
        )
        outputs.append(output_png)
        reports.append(report)
        print(f"saved shape PNG for mu={float(mu):g}: {output_png}")
        print(f"saved report for mu={float(mu):g}: {report}")
        for case in cases:
            valid = (
                case.thickness_ratio_1 <= float(args.thickness_ratio_limit)
                and case.thickness_ratio_2 <= float(args.thickness_ratio_limit)
            )
            status = "ok" if valid else "WARNING thin-rod violation"
            print(
                f"  eta={case.eta:g}: Lambda={case.Lambda:.10g}, "
                f"sorted={case.current_sorted_index}, tracking={case.tracking_step_status}, "
                f"2*r1/l1={case.thickness_ratio_1:.6g}, "
                f"2*r2/l2={case.thickness_ratio_2:.6g}, {status}"
            )

    return {"png": outputs, "reports": reports}


if __name__ == "__main__":
    main()
