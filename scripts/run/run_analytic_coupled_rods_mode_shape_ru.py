from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


# ============================================================
# USER PARAMETERS
# Change these values for ordinary runs.
# CLI arguments, if provided, override these defaults.
# ============================================================

DEFAULT_MODE_NUMBER = 5
DEFAULT_BETA = 30.0
DEFAULT_MU = 0.0
DEFAULT_EPSILON = 0.0025
DEFAULT_L_TOTAL = 2.0

DEFAULT_PLOT_KIND = "components"  # "full", "transverse", "components"
DEFAULT_MODE_SCALE = 0.22
DEFAULT_NORMALIZE = "auto"

DEFAULT_OUTPUT = None
DEFAULT_DPI = 240
DEFAULT_FIGSIZE = None
DEFAULT_SHOW = False
DEFAULT_PRINT_DIAGNOSTICS = True
DEFAULT_SAVE_SAMPLES_CSV = None
DEFAULT_SAVE_DIAGNOSTICS_CSV = None


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402
from my_project.analytic.solvers import find_first_n_roots  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    NEAR_ZERO_NORM,
    analytic_arm_energy_diagnostics,
    analytic_null_vector,
    endpoint_consistency_diagnostics,
    normalize_components,
    reconstruct_analytic_components,
    unique_sorted_roots,
)


PLOT_KINDS = ("full", "transverse", "components")
NORMALIZE_KINDS = ("auto", "max-full", "max-transverse", "none")
RESULTS_DIR = REPO_ROOT / "results"
NUM_SAMPLES = 401
RIGHT_COORDINATE_FOR_REPORTING = "joint-to-external"
DEFAULT_GEOMETRY_FIGSIZE = (8.0, 4.2)
DEFAULT_COMPONENTS_FIGSIZE = (8.0, 5.4)


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def mu_label(mu: float) -> str:
    rounded = round(float(mu))
    if abs(float(mu) - rounded) < 1e-9:
        return str(int(rounded))
    return f"{float(mu):g}"


def default_output_path(
    *,
    mode_number: int,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    plot_kind: str,
    mode_scale: float,
) -> Path:
    return RESULTS_DIR / (
        f"analytic_mode_shape_{plot_kind}_mode{int(mode_number)}_beta{filename_number_token(beta_deg)}"
        f"_mu{filename_number_token(mu_value)}_eps{filename_number_token(epsilon)}"
        f"_scale{filename_number_token(mode_scale)}_ru.png"
    )


def resolve_repo_path(value: str | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def resolve_output_path(
    value: str | None,
    *,
    mode_number: int,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    plot_kind: str,
    mode_scale: float,
) -> Path:
    resolved = resolve_repo_path(value)
    if resolved is not None:
        return resolved
    return default_output_path(
        mode_number=mode_number,
        beta_deg=beta_deg,
        mu_value=mu_value,
        epsilon=epsilon,
        plot_kind=plot_kind,
        mode_scale=mode_scale,
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Build analytic mode shapes for two rigidly coupled rods directly from "
            "the determinant nullspace. This runner does not require FEM comparison."
        ),
    )
    parser.add_argument(
        "--mode-number",
        type=int,
        default=DEFAULT_MODE_NUMBER,
        help=(
            "Analytic root/mode number in ascending Lambda order. This is the analytic "
            "root/mode number, not a FEM descendant label."
        ),
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA, help="Coupling angle in degrees.")
    parser.add_argument("--mu", type=float, default=DEFAULT_MU, help="Length-asymmetry parameter.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Slenderness/coupling epsilon.")
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL, help="Total length L1 + L2 for plotting.")
    parser.add_argument("--plot-kind", choices=PLOT_KINDS, default=DEFAULT_PLOT_KIND)
    parser.add_argument("--mode-scale", type=float, default=DEFAULT_MODE_SCALE)
    parser.add_argument("--normalize", choices=NORMALIZE_KINDS, default=DEFAULT_NORMALIZE)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("WIDTH", "HEIGHT"), default=DEFAULT_FIGSIZE)
    parser.add_argument("--show", action="store_true", default=DEFAULT_SHOW)
    parser.add_argument("--save-samples-csv", default=DEFAULT_SAVE_SAMPLES_CSV)
    parser.add_argument("--save-diagnostics-csv", default=DEFAULT_SAVE_DIAGNOSTICS_CSV)
    diagnostics_group = parser.add_mutually_exclusive_group()
    diagnostics_group.add_argument("--print-diagnostics", dest="print_diagnostics", action="store_true")
    diagnostics_group.add_argument("--no-print-diagnostics", dest="print_diagnostics", action="store_false")
    parser.set_defaults(print_diagnostics=DEFAULT_PRINT_DIAGNOSTICS)
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    if args.mode_number <= 0:
        parser.error("--mode-number must be positive.")
    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if not (-1.0 < args.mu < 1.0):
        parser.error("--mu must satisfy -1 < mu < 1 so both arm lengths are positive.")
    if args.dpi <= 0:
        parser.error("--dpi must be positive.")
    return args


def resolve_normalization(plot_kind: str, normalize: str) -> str:
    if normalize != "auto":
        return normalize
    if plot_kind == "transverse":
        return "max-transverse"
    return "max-full"


def find_mode_lambda(
    *,
    mode_number: int,
    beta_rad: float,
    mu_value: float,
    epsilon: float,
) -> tuple[float, list[float]]:
    roots = unique_sorted_roots(
        find_first_n_roots(
            beta_rad,
            mu_value,
            epsilon,
            mode_number,
        )
    )
    if len(roots) < mode_number:
        raise RuntimeError(f"Only found {len(roots)} analytic roots; cannot select mode {mode_number}.")
    return float(roots[mode_number - 1]), roots


def title_text(*, beta_deg: float, mu_value: float, epsilon: float, mode_number: int, Lambda: float) -> str:
    return (
        f"β = {beta_deg:g}°, μ = {mu_label(mu_value)}, ε = {epsilon:g}\n"
        f"аналитическая мода {int(mode_number)}, Λ = {float(Lambda):.6g}"
    )


def arm_lengths(*, mu_value: float, l_total: float) -> tuple[float, float, float]:
    ell = float(l_total) / 2.0
    return ell, ell * (1.0 - float(mu_value)), ell * (1.0 + float(mu_value))


def base_coordinates(
    *,
    s_norm: np.ndarray,
    beta_rad: float,
    mu_value: float,
    l_total: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    _, left_length, right_length = arm_lengths(mu_value=mu_value, l_total=l_total)
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
    plot_kind: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
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

    if plot_kind == "transverse":
        u_left = np.zeros_like(u_left)
        u_right = np.zeros_like(u_right)

    x_left_def = x_left + float(mode_scale) * u_left
    y_left_def = y_left + float(mode_scale) * w_left
    x_right_def = x_right + float(mode_scale) * (u_right * tangent_right[0] + w_right * normal_right[0])
    y_right_def = y_right + float(mode_scale) * (u_right * tangent_right[1] + w_right * normal_right[1])
    return x_left, y_left, x_right, y_right, x_left_def, y_left_def, x_right_def, y_right_def


def geometry_axis_limits(arrays: Sequence[np.ndarray]) -> tuple[tuple[float, float], tuple[float, float]]:
    x_values = np.concatenate([np.asarray(arrays[idx], dtype=float) for idx in range(0, len(arrays), 2)])
    y_values = np.concatenate([np.asarray(arrays[idx], dtype=float) for idx in range(1, len(arrays), 2)])
    x_span = max(float(np.max(x_values) - np.min(x_values)), 1.0)
    y_span = max(float(np.max(y_values) - np.min(y_values)), 0.5)
    return (
        (float(np.min(x_values) - 0.06 * x_span), float(np.max(x_values) + 0.06 * x_span)),
        (float(np.min(y_values) - 0.10 * y_span), float(np.max(y_values) + 0.10 * y_span)),
    )


def plot_geometry(
    output_path: Path,
    *,
    components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    beta_rad: float,
    mu_value: float,
    l_total: float,
    mode_scale: float,
    plot_kind: str,
    title: str,
    dpi: int,
    figsize: tuple[float, float] | None,
    show: bool,
) -> None:
    coords = deformed_coordinates(
        components=components,
        s_norm=s_norm,
        beta_rad=beta_rad,
        mu_value=mu_value,
        l_total=l_total,
        mode_scale=mode_scale,
        plot_kind=plot_kind,
    )
    x_left, y_left, x_right, y_right, x_left_def, y_left_def, x_right_def, y_right_def = coords
    x_limits, y_limits = geometry_axis_limits(coords)
    fig, ax = plt.subplots(figsize=figsize or DEFAULT_GEOMETRY_FIGSIZE)
    ax.plot(x_left, y_left, color="0.78", linestyle="--", linewidth=1.1, label="недеформированная геометрия")
    ax.plot(x_right, y_right, color="0.78", linestyle="--", linewidth=1.1)
    ax.plot(x_left_def, y_left_def, color="#1f77b4", linewidth=2.2, label="левое плечо")
    ax.plot(x_right_def, y_right_def, color="#ff7f0e", linewidth=2.2, label="правое плечо")
    ax.scatter([x_left[-1]], [y_left[-1]], color="black", s=14, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title, fontsize=10.6)
    ax.grid(True, alpha=0.22)
    ax.legend(fontsize=8.8, loc="best")
    fig.tight_layout(pad=0.8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def plot_components(
    output_path: Path,
    *,
    components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    title: str,
    dpi: int,
    figsize: tuple[float, float] | None,
    show: bool,
) -> None:
    fig, axes = plt.subplots(2, 1, figsize=figsize or DEFAULT_COMPONENTS_FIGSIZE, sharex=True)
    axes[0].plot(s_norm, components["u_left"], color="#1f77b4", linewidth=2.0, label="левое плечо")
    axes[0].plot(s_norm, components["u_right"], color="#ff7f0e", linewidth=2.0, label="правое плечо")
    axes[0].set_ylabel("локальная продольная")
    axes[0].grid(True, alpha=0.22)
    axes[0].legend(fontsize=9, loc="best")

    axes[1].plot(s_norm, components["w_left"], color="#1f77b4", linewidth=2.0, label="левое плечо")
    axes[1].plot(s_norm, components["w_right"], color="#ff7f0e", linewidth=2.0, label="правое плечо")
    axes[1].set_xlabel("локальная координата s / L")
    axes[1].set_ylabel("локальная поперечная")
    axes[1].grid(True, alpha=0.22)
    axes[1].legend(fontsize=9, loc="best")

    fig.suptitle(title, fontsize=10.6)
    fig.tight_layout(pad=0.8)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def write_samples_csv(path: Path, *, s_norm: np.ndarray, components: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["s_norm", "u_left", "w_left", "u_right", "w_right"]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for idx, s_value in enumerate(s_norm):
            writer.writerow(
                {
                    "s_norm": float(s_value),
                    "u_left": float(components["u_left"][idx]),
                    "w_left": float(components["w_left"][idx]),
                    "u_right": float(components["u_right"][idx]),
                    "w_right": float(components["w_right"][idx]),
                }
            )


def build_diagnostics_row(
    *,
    mode_number: int,
    beta_deg: float,
    mu_value: float,
    epsilon: float,
    l_total: float,
    Lambda: float,
    coeff: np.ndarray,
    roots_found: int,
    smallest_singular_value: float,
    singular_value_ratio: float,
    normalization: str,
    normalization_scale: float,
    output_path: Path,
    samples_csv: Path | None,
    diagnostics_csv: Path | None,
    consistency: dict[str, float],
    energy: dict[str, float],
) -> dict[str, float | int | str]:
    return {
        "mode_number": int(mode_number),
        "beta": float(beta_deg),
        "mu": float(mu_value),
        "epsilon": float(epsilon),
        "l_total": float(l_total),
        "lambda": float(Lambda),
        "roots_found": int(roots_found),
        "unknown_ordering": "A1,B1,A2,B2,P1,P2",
        "right_coordinate_convention": RIGHT_COORDINATE_FOR_REPORTING,
        "normalization": normalization,
        "normalization_scale": float(normalization_scale),
        "coeff_A1": float(coeff[0]),
        "coeff_B1": float(coeff[1]),
        "coeff_A2": float(coeff[2]),
        "coeff_B2": float(coeff[3]),
        "coeff_P1": float(coeff[4]),
        "coeff_P2": float(coeff[5]),
        "smallest_singular_value": float(smallest_singular_value),
        "singular_value_ratio": float(singular_value_ratio),
        "max_matrix_residual": float(consistency["max_matrix_residual"]),
        "max_field_residual": float(consistency["max_field_residual"]),
        "matrix_field_residual_max_abs_difference": float(consistency["matrix_field_residual_max_abs_difference"]),
        "external_clamp_residual": float(consistency["external_clamp_residual"]),
        "analytic_left_axial_energy": float(energy["left_axial_energy"]),
        "analytic_right_axial_energy": float(energy["right_axial_energy"]),
        "analytic_left_bending_energy": float(energy["left_bending_energy"]),
        "analytic_right_bending_energy": float(energy["right_bending_energy"]),
        "analytic_total_axial_energy": float(energy["total_axial_energy"]),
        "analytic_total_bending_energy": float(energy["total_bending_energy"]),
        "analytic_total_energy": float(energy["total_energy"]),
        "analytic_axial_energy_fraction": float(energy["axial_energy_fraction"]),
        "analytic_right_axial_share": float(energy["right_axial_share"]),
        "analytic_right_bending_share": float(energy["right_bending_share"]),
        "output_png": str(output_path),
        "samples_csv": "" if samples_csv is None else str(samples_csv),
        "diagnostics_csv": "" if diagnostics_csv is None else str(diagnostics_csv),
    }


def write_diagnostics_csv(path: Path, row: dict[str, float | int | str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def print_diagnostics(row: dict[str, float | int | str]) -> None:
    print("Analytic coupled-rods mode shape")
    print(f"mode_number: {int(row['mode_number'])}")
    print(f"beta: {float(row['beta']):g} deg")
    print(f"mu: {float(row['mu']):g}")
    print(f"epsilon: {float(row['epsilon']):g}")
    print(f"lambda: {float(row['lambda']):.10g}")
    print(f"unknown_ordering: {row['unknown_ordering']}")
    print(f"right_coordinate_convention: {row['right_coordinate_convention']}")
    print(f"normalization: {row['normalization']}")
    print(f"smallest_singular_value: {float(row['smallest_singular_value']):.6e}")
    print(f"singular_value_ratio: {float(row['singular_value_ratio']):.6e}")
    print(f"max_matrix_residual: {float(row['max_matrix_residual']):.6e}")
    print(f"max_field_residual: {float(row['max_field_residual']):.6e}")
    print(
        "matrix_field_residual_max_abs_difference: "
        f"{float(row['matrix_field_residual_max_abs_difference']):.6e}"
    )
    print(f"external_clamp_residual: {float(row['external_clamp_residual']):.6e}")
    print(f"analytic_axial_energy_fraction: {float(row['analytic_axial_energy_fraction']):.6e}")
    print(f"analytic_right_axial_share: {float(row['analytic_right_axial_share']):.6e}")
    print(f"analytic_right_bending_share: {float(row['analytic_right_bending_share']):.6e}")
    print(f"saved PNG: {row['output_png']}")
    if row["samples_csv"] != "":
        print(f"saved samples CSV: {row['samples_csv']}")
    if row["diagnostics_csv"] != "":
        print(f"saved diagnostics CSV: {row['diagnostics_csv']}")


def main(argv: Sequence[str] | None = None) -> dict[str, float | int | str]:
    args = parse_args(argv)
    beta_rad = float(np.deg2rad(args.beta))
    Lambda, roots = find_mode_lambda(
        mode_number=args.mode_number,
        beta_rad=beta_rad,
        mu_value=args.mu,
        epsilon=args.epsilon,
    )
    matrix = assemble_clamped_coupled_matrix(Lambda, beta_rad, args.mu, args.epsilon)
    coeff, smallest_singular_value, singular_value_ratio = analytic_null_vector(matrix)
    s_norm = np.linspace(0.0, 1.0, NUM_SAMPLES)
    components_raw = reconstruct_analytic_components(
        Lambda,
        mu_value=args.mu,
        epsilon=args.epsilon,
        coeff=coeff,
        s_norm=s_norm,
        right_coordinate=RIGHT_COORDINATE_FOR_REPORTING,
    )
    resolved_normalization = resolve_normalization(args.plot_kind, args.normalize)
    components_plot, normalization_scale = normalize_components(
        components_raw,
        plot_kind=args.plot_kind,
        normalize=resolved_normalization,
    )
    consistency = endpoint_consistency_diagnostics(
        matrix,
        coeff,
        Lambda=Lambda,
        beta_rad=beta_rad,
        mu_value=args.mu,
        epsilon=args.epsilon,
    )
    energy = analytic_arm_energy_diagnostics(
        components_plot,
        mu_value=args.mu,
        epsilon=args.epsilon,
        s_norm=s_norm,
    )

    output_path = resolve_output_path(
        args.output,
        mode_number=args.mode_number,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
        plot_kind=args.plot_kind,
        mode_scale=args.mode_scale,
    )
    samples_csv = resolve_repo_path(args.save_samples_csv)
    diagnostics_csv = resolve_repo_path(args.save_diagnostics_csv)
    figsize = None if args.figsize is None else (float(args.figsize[0]), float(args.figsize[1]))
    title = title_text(
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
        mode_number=args.mode_number,
        Lambda=Lambda,
    )

    if args.plot_kind == "components":
        plot_components(
            output_path,
            components=components_plot,
            s_norm=s_norm,
            title=title,
            dpi=args.dpi,
            figsize=figsize,
            show=bool(args.show),
        )
    else:
        plot_geometry(
            output_path,
            components=components_plot,
            s_norm=s_norm,
            beta_rad=beta_rad,
            mu_value=args.mu,
            l_total=args.l_total,
            mode_scale=args.mode_scale,
            plot_kind=args.plot_kind,
            title=title,
            dpi=args.dpi,
            figsize=figsize,
            show=bool(args.show),
        )

    if samples_csv is not None:
        write_samples_csv(samples_csv, s_norm=s_norm, components=components_plot)

    row = build_diagnostics_row(
        mode_number=args.mode_number,
        beta_deg=args.beta,
        mu_value=args.mu,
        epsilon=args.epsilon,
        l_total=args.l_total,
        Lambda=Lambda,
        coeff=coeff,
        roots_found=len(roots),
        smallest_singular_value=smallest_singular_value,
        singular_value_ratio=singular_value_ratio,
        normalization=resolved_normalization,
        normalization_scale=normalization_scale,
        output_path=output_path,
        samples_csv=samples_csv,
        diagnostics_csv=diagnostics_csv,
        consistency=consistency,
        energy=energy,
    )

    if diagnostics_csv is not None:
        write_diagnostics_csv(diagnostics_csv, row)
    if args.print_diagnostics:
        print_diagnostics(row)
    return row


if __name__ == "__main__":
    main()
