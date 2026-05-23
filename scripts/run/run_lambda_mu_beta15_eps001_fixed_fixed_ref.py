from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    DEFAULT_MU_STEPS,
    DEFAULT_N_SOLVE,
    DEFAULT_N_TRACK,
    branch_id_from_base_sorted_index,
    dense_mu_values_for_targets,
    track_mu_sweep,
)


DEFAULT_BETA = 15.0
DEFAULT_EPSILON = 0.01
DEFAULT_L_TOTAL = 2.0
DEFAULT_NUM_MODES = 6
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.9
DEFAULT_MU_STEP = 0.005
DEFAULT_Y_MAX = 11.2
DEFAULT_STYLE_SPLIT_MU = 0.6
DEFAULT_OUTPUT = None
DEFAULT_CSV_OUTPUT = None
DEFAULT_REF_CSV_OUTPUT = None
DEFAULT_SHOW = False
DEFAULT_ALLOW_LOW_MAC = False
DEFAULT_FIGURE_DPI = 600
ARTICLE_X_LABEL = r"$\mu$"
ARTICLE_Y_LABEL = r"$\Lambda$"
KNOWN_FIXED_FIXED_ALPHA_PREFIX = np.array(
    [
        4.730040744862704,
        7.853204624095838,
        10.99560783800167,
    ],
    dtype=float,
)

MAIN_FIELDNAMES = [
    "beta",
    "epsilon",
    "L_total",
    "l_base",
    "r",
    "mu",
    "mode_rank",
    "branch_id",
    "base_sorted_index",
    "current_sorted_index",
    "Lambda",
    "warning_flag",
    "mac_to_previous",
    "relative_lambda_jump",
    "relative_gap_lower",
    "relative_gap_upper",
    "region_style",
]

REFERENCE_FIELDNAMES = [
    "reference_type",
    "mode_index",
    "alpha",
    "Lambda_ref",
    "boundary_condition",
    "Lref",
    "epsilon",
    "r",
    "l_base",
    "equation",
    "line_style",
    "horizontal_in_mu",
    "plotted_y_max",
]


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def split_mu_token(value: float) -> str:
    if abs(float(value) - 0.6) <= 1e-12:
        return "mu06split"
    return f"mu{filename_number_token(value)}split"


def default_output_paths(*, beta: float, epsilon: float, L_total: float, split_mu: float) -> tuple[Path, Path, Path]:
    stem = (
        f"lambda_mu_beta{filename_number_token(beta)}"
        f"_eps{filename_number_token(epsilon)}"
        f"_fixed_fixed_L{filename_number_token(L_total)}"
        f"_compact_{split_mu_token(split_mu)}_article"
    )
    results_dir = REPO_ROOT / "results"
    return (
        results_dir / f"{stem}.png",
        results_dir / f"{stem}.csv",
        results_dir / f"{stem}_refs.csv",
    )


def resolve_repo_path(value: str | Path | None) -> Path | None:
    if value is None:
        return None
    path = Path(value)
    if path.is_absolute():
        return path
    return REPO_ROOT / path


def console_safe(text: str) -> str:
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    return str(text).encode(encoding, errors="backslashreplace").decode(encoding)


def circular_radius_from_epsilon(*, epsilon: float, L_total: float) -> float:
    l_base = 0.5 * float(L_total)
    return 2.0 * float(epsilon) * l_base


def mu_grid(*, mu_min: float, mu_max: float, mu_step: float, include_values: Sequence[float] = ()) -> np.ndarray:
    values = np.arange(float(mu_min), float(mu_max) + 0.5 * float(mu_step), float(mu_step), dtype=float)
    if values.size == 0 or abs(float(values[-1]) - float(mu_max)) > 1e-10:
        values = np.append(values, float(mu_max))
    for value in include_values:
        if float(mu_min) - 1e-12 <= float(value) <= float(mu_max) + 1e-12:
            values = np.append(values, float(value))
    values[0] = float(mu_min)
    return np.unique(np.round(values, 12))


def region_style_for_mu(mu: float, *, split_mu: float) -> str:
    return "solid" if float(mu) <= float(split_mu) + 1e-12 else "dashed"


def fixed_fixed_alpha_roots(n_roots: int) -> np.ndarray:
    mp.mp.dps = 80

    def characteristic(alpha):
        return mp.cosh(alpha) * mp.cos(alpha) - 1

    roots = []
    for mode_index in range(1, int(n_roots) + 1):
        center = (mode_index + mp.mpf("0.5")) * mp.pi
        roots.append(float(mp.findroot(characteristic, (center - mp.mpf("0.1"), center + mp.mpf("0.1")))))
    return np.asarray(roots, dtype=float)


def validate_fixed_fixed_roots(alphas: np.ndarray) -> float:
    if len(alphas) < len(KNOWN_FIXED_FIXED_ALPHA_PREFIX):
        raise RuntimeError("Need at least three fixed-fixed roots for prefix validation.")
    errors = np.abs(alphas[: len(KNOWN_FIXED_FIXED_ALPHA_PREFIX)] - KNOWN_FIXED_FIXED_ALPHA_PREFIX)
    max_error = float(np.max(errors))
    if not np.allclose(
        alphas[: len(KNOWN_FIXED_FIXED_ALPHA_PREFIX)],
        KNOWN_FIXED_FIXED_ALPHA_PREFIX,
        rtol=0.0,
        atol=5e-12,
    ):
        raise RuntimeError(
            "Computed fixed-fixed roots failed validation against the known prefix: "
            f"computed={alphas[:len(KNOWN_FIXED_FIXED_ALPHA_PREFIX)]}, "
            f"expected={KNOWN_FIXED_FIXED_ALPHA_PREFIX}, max_abs_error={max_error:.6g}"
        )
    return max_error


def fixed_fixed_reference_rows(
    *,
    epsilon: float,
    L_total: float,
    reference_count: int,
    y_max: float,
) -> tuple[list[dict[str, float | int | str]], np.ndarray, np.ndarray, float]:
    l_base = 0.5 * float(L_total)
    radius = circular_radius_from_epsilon(epsilon=float(epsilon), L_total=float(L_total))
    n_roots = max(int(reference_count), len(KNOWN_FIXED_FIXED_ALPHA_PREFIX))
    alphas = fixed_fixed_alpha_roots(n_roots)
    max_prefix_error = validate_fixed_fixed_roots(alphas)
    lambdas = alphas * l_base / float(L_total)
    keep = np.arange(int(reference_count), dtype=int)

    rows: list[dict[str, float | int | str]] = []
    for idx in keep:
        rows.append(
            {
                "reference_type": "bending_fixed_fixed_single_beam",
                "mode_index": int(idx) + 1,
                "alpha": float(alphas[idx]),
                "Lambda_ref": float(lambdas[idx]),
                "boundary_condition": "fixed-fixed",
                "Lref": float(L_total),
                "epsilon": float(epsilon),
                "r": float(radius),
                "l_base": float(l_base),
                "equation": "cosh(alpha)*cos(alpha)=1",
                "line_style": "dotted_horizontal",
                "horizontal_in_mu": 1,
                "plotted_y_max": float(y_max),
            }
        )
    return rows, alphas, lambdas, max_prefix_error


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot analytic Lambda(mu) at fixed beta, epsilon=0.01 with fixed-fixed "
            "single-beam L=2 m horizontal reference lines."
        ),
    )
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA, help="Fixed coupling angle in degrees.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON, help="Slenderness parameter epsilon.")
    parser.add_argument("--L-total", type=float, default=DEFAULT_L_TOTAL, help="Total coupled-system length in metres.")
    parser.add_argument("--num-modes", type=int, default=DEFAULT_NUM_MODES, help="Number of canonical analytic branches to draw.")
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--y-max", type=float, default=DEFAULT_Y_MAX, help="Compact upper y-axis limit; use <= 0 for automatic.")
    parser.add_argument(
        "--style-split-mu",
        type=float,
        default=DEFAULT_STYLE_SPLIT_MU,
        help="Draw coupled branches solid up to this mu and dashed after it.",
    )
    parser.add_argument("--output", default=DEFAULT_OUTPUT, help="Figure output path. Default is deterministic under results/.")
    parser.add_argument("--csv-output", default=DEFAULT_CSV_OUTPUT, help="Main curve CSV output path. Default matches the PNG stem.")
    parser.add_argument("--ref-csv-output", default=DEFAULT_REF_CSV_OUTPUT, help="Reference-line CSV output path. Default matches the PNG stem.")
    parser.add_argument("--show", action="store_true", default=DEFAULT_SHOW, help="Also display the plot window.")
    parser.add_argument(
        "--allow-low-mac",
        action="store_true",
        default=DEFAULT_ALLOW_LOW_MAC,
        help="Allow exploratory output even if canonical tracking falls below the MAC warning threshold.",
    )
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.L_total <= 0.0:
        parser.error("--L-total must be positive.")
    if args.num_modes <= 0:
        parser.error("--num-modes must be positive.")
    if args.mu_step <= 0.0:
        parser.error("--mu-step must be positive.")
    if args.mu_max < args.mu_min:
        parser.error("--mu-max must be greater than or equal to --mu-min.")
    if args.mu_min < -1e-12:
        parser.error("--mu-min must be non-negative for canonical beta-then-mu analytic tracking.")
    if not (args.mu_min - 1e-12 <= args.style_split_mu <= args.mu_max + 1e-12):
        parser.error("--style-split-mu must lie inside the plotted mu range.")
    return args


def track_and_collect_rows(
    *,
    beta: float,
    epsilon: float,
    L_total: float,
    num_modes: int,
    mu_values: np.ndarray,
    allow_low_mac: bool,
    style_split_mu: float,
) -> tuple[list[dict[str, float | int | str]], list[str], np.ndarray]:
    n_track = max(DEFAULT_N_TRACK, int(num_modes))
    tracking_mu_values = dense_mu_values_for_targets(
        mu_values,
        mu_steps=max(DEFAULT_MU_STEPS, int(len(mu_values))),
    )
    branch_ids = [branch_id_from_base_sorted_index(index) for index in range(1, int(num_modes) + 1)]
    result = track_mu_sweep(
        epsilon=float(epsilon),
        beta=float(beta),
        mu_values=tracking_mu_values,
        n_track=n_track,
        n_solve=max(DEFAULT_N_SOLVE, n_track),
        shape_metric="full",
        allow_low_mac=bool(allow_low_mac),
        required_branch_ids=branch_ids,
    )

    l_base = 0.5 * float(L_total)
    radius = circular_radius_from_epsilon(epsilon=float(epsilon), L_total=float(L_total))
    lambda_grid = np.full((int(num_modes), len(mu_values)), np.nan, dtype=float)
    rows: list[dict[str, float | int | str]] = []
    for mode_rank, branch_id in enumerate(branch_ids, start=1):
        for col, mu in enumerate(mu_values):
            point = result.point_at(branch_id, beta=float(beta), mu=float(mu))
            lambda_grid[mode_rank - 1, col] = float(point.Lambda)
            rows.append(
                {
                    "beta": float(beta),
                    "epsilon": float(epsilon),
                    "L_total": float(L_total),
                    "l_base": float(l_base),
                    "r": float(radius),
                    "mu": float(mu),
                    "mode_rank": int(mode_rank),
                    "branch_id": branch_id,
                    "base_sorted_index": int(mode_rank),
                    "current_sorted_index": int(point.current_sorted_index),
                    "Lambda": float(point.Lambda),
                    "warning_flag": str(point.warning_flag),
                    "mac_to_previous": float(point.mac_to_previous),
                    "relative_lambda_jump": float(point.relative_lambda_jump),
                    "relative_gap_lower": float(point.relative_gap_lower),
                    "relative_gap_upper": float(point.relative_gap_upper),
                    "region_style": region_style_for_mu(float(mu), split_mu=float(style_split_mu)),
                }
            )
    return rows, list(result.warnings), lambda_grid


def write_csv(path: Path, rows: list[dict[str, float | int | str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def plot_lambda_mu(
    *,
    beta: float,
    epsilon: float,
    L_total: float,
    radius: float,
    mu_values: np.ndarray,
    lambda_grid: np.ndarray,
    reference_rows: list[dict[str, float | int | str]],
    y_max: float,
    style_split_mu: float,
    output_path: Path,
    show: bool,
) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    solid_mask = mu_values <= float(style_split_mu) + 1e-12
    dashed_mask = mu_values >= float(style_split_mu) - 1e-12

    for row_idx in range(lambda_grid.shape[0]):
        color = colors[row_idx % len(colors)]
        # Style split only marks interpretation range:
        # solid for mu <= 0.6, where the thin-rod interpretation is considered acceptable;
        # dashed for mu > 0.6, shown for completeness when the shorter arm is no longer confidently thin.
        ax.plot(mu_values[solid_mask], lambda_grid[row_idx, solid_mask], lw=1.6, ls="-", color=color, zorder=3, rasterized=False)
        ax.plot(mu_values[dashed_mask], lambda_grid[row_idx, dashed_mask], lw=1.6, ls="--", color=color, zorder=3, rasterized=False)

    for row in reference_rows:
        ax.axhline(
            float(row["Lambda_ref"]),
            color="0.72",
            lw=1.25,
            ls=":",
            zorder=1,
            rasterized=False,
        )

    ax.set_xlim(float(mu_values[0]), float(mu_values[-1]))
    ax.set_ylim(0.0, float(y_max))
    ax.set_xlabel(ARTICLE_X_LABEL)
    ax.set_ylabel(ARTICLE_Y_LABEL)
    ax.grid(True, color="0.86", linewidth=0.6)
    ax.tick_params(axis="both", labelsize=10, width=0.9, length=5)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontstyle("normal")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    save_format = output_path.suffix.lstrip(".").lower() or None
    fig.savefig(output_path, format=save_format, dpi=DEFAULT_FIGURE_DPI, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)


def first_axial_lambda(*, epsilon: float, L_total: float) -> float:
    l_base = 0.5 * float(L_total)
    return float(np.sqrt(np.pi * l_base / (float(epsilon) * float(L_total))))


def main(argv: Sequence[str] | None = None) -> dict[str, Path | float | int]:
    args = parse_args(argv)
    default_figure, default_csv, default_ref_csv = default_output_paths(
        beta=float(args.beta),
        epsilon=float(args.epsilon),
        L_total=float(args.L_total),
        split_mu=float(args.style_split_mu),
    )
    output_path = resolve_repo_path(args.output) or default_figure
    csv_path = resolve_repo_path(args.csv_output) or default_csv
    ref_csv_path = resolve_repo_path(args.ref_csv_output) or default_ref_csv

    mu_values = mu_grid(
        mu_min=args.mu_min,
        mu_max=args.mu_max,
        mu_step=args.mu_step,
        include_values=[float(args.style_split_mu)],
    )
    main_rows, tracking_warnings, lambda_grid = track_and_collect_rows(
        beta=float(args.beta),
        epsilon=float(args.epsilon),
        L_total=float(args.L_total),
        num_modes=int(args.num_modes),
        mu_values=mu_values,
        allow_low_mac=bool(args.allow_low_mac),
        style_split_mu=float(args.style_split_mu),
    )

    finite_lambdas = lambda_grid[np.isfinite(lambda_grid)]
    auto_y_max = 1.08 * float(np.max(finite_lambdas)) if finite_lambdas.size else DEFAULT_Y_MAX
    y_max = float(args.y_max) if float(args.y_max) > 0.0 else max(1.0, auto_y_max)
    radius = circular_radius_from_epsilon(epsilon=float(args.epsilon), L_total=float(args.L_total))
    reference_rows, alphas, reference_lambdas, root_prefix_error = fixed_fixed_reference_rows(
        epsilon=float(args.epsilon),
        L_total=float(args.L_total),
        reference_count=int(args.num_modes),
        y_max=y_max,
    )

    write_csv(csv_path, main_rows, MAIN_FIELDNAMES)
    write_csv(ref_csv_path, reference_rows, REFERENCE_FIELDNAMES)
    plot_lambda_mu(
        beta=float(args.beta),
        epsilon=float(args.epsilon),
        L_total=float(args.L_total),
        radius=radius,
        mu_values=mu_values,
        lambda_grid=lambda_grid,
        reference_rows=reference_rows,
        y_max=y_max,
        style_split_mu=float(args.style_split_mu),
        output_path=output_path,
        show=bool(args.show),
    )

    axial_1 = first_axial_lambda(epsilon=float(args.epsilon), L_total=float(args.L_total))
    axial_inside = 0.0 <= axial_1 <= y_max

    print(f"saved figure: {output_path}")
    print(f"saved CSV: {csv_path}")
    print(f"saved reference CSV: {ref_csv_path}")
    print(f"beta: {float(args.beta):g} deg")
    print(f"epsilon: {float(args.epsilon):g}")
    print(f"L_total: {float(args.L_total):g} m")
    print(f"l_base: {0.5 * float(args.L_total):g} m")
    print(f"circular radius r = 2*epsilon*l_base: {radius:.12g} m")
    print(f"num_modes: {int(args.num_modes)}")
    print(f"mu points: {len(mu_values)}")
    print(f"y limits: [0, {y_max:g}]")
    print(f"coupled branch style split: solid for mu <= {float(args.style_split_mu):g}; dashed for mu > {float(args.style_split_mu):g}")
    print(console_safe(f'axis labels: x="{ARTICLE_X_LABEL}", y="{ARTICLE_Y_LABEL}"'))
    print("legend: absent")
    print("title: absent")
    print(f"mu split point present in grid: {bool(np.any(np.isclose(mu_values, float(args.style_split_mu), rtol=0.0, atol=1e-12)))}")
    print(f"fixed-fixed reference roots plotted: {len(reference_rows)}")
    print(f"fixed-fixed alpha prefix max abs error: {root_prefix_error:.6g}")
    print("first 5 fixed-fixed alpha roots: " + ", ".join(f"{value:.12g}" for value in alphas[:5]))
    print("first 5 fixed-fixed Lambda_ref values: " + ", ".join(f"{value:.12g}" for value in reference_lambdas[:5]))
    print("reference lines horizontal in mu: True")
    print(
        "first fixed-fixed axial single-rod Lambda "
        f"(sqrt(n*pi*l/(epsilon*Lref)), n=1): {axial_1:.12g}; inside plotted y-range: {axial_inside}"
    )
    print(f"allow_low_mac: {bool(args.allow_low_mac)}")
    if tracking_warnings:
        print(f"tracking warnings: {len(tracking_warnings)}")
        for warning in tracking_warnings[:5]:
            print(f"  {warning}")
        if len(tracking_warnings) > 5:
            print(f"  ... {len(tracking_warnings) - 5} more warnings")
    else:
        print("tracking warnings: none")
    return {
        "output_png": output_path,
        "output_figure": output_path,
        "output_csv": csv_path,
        "reference_csv": ref_csv_path,
        "beta": float(args.beta),
        "epsilon": float(args.epsilon),
        "L_total": float(args.L_total),
        "r": float(radius),
        "num_modes": int(args.num_modes),
        "y_max": float(y_max),
    }


if __name__ == "__main__":
    main()
