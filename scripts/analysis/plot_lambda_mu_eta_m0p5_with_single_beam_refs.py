from __future__ import annotations

import argparse
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

from my_project.analytic.FreqMuNet import roots_clamped_supported  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    find_first_n_roots_eta,
    thickness_mismatch_factors,
)
from my_project.analytic.solvers import fixed_fixed_lambdas  # noqa: E402
from scripts.lib.diagnostic_common import (  # noqa: E402
    inclusive_grid,
    number_text as diagnostic_number_text,
    number_token,
    write_dict_rows_csv,
)


DEFAULT_EPSILON = 0.0025
DEFAULT_ETA = -0.5
DEFAULT_BETAS = (0.0, 5.0, 10.0)
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.9
DEFAULT_MU_STEP = 0.005
DEFAULT_NUM_MODES = 6
DEFAULT_NUM_REF_MODES = 6
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
SMOKE_OUTPUT_DIR = DEFAULT_OUTPUT_DIR / "_smoke"
DEFAULT_ROOT_SCAN_STEP = 0.02
DEFAULT_ROOT_LMAX0 = 35.0
DEFAULT_Y_MARGIN = 1.15

CSV_FIELDNAMES = [
    "curve_group",
    "beta_deg",
    "epsilon",
    "eta",
    "mu",
    "mode_index",
    "Lambda",
    "boundary_condition",
    "rod_id",
    "rod_label",
    "length_factor",
    "tau",
    "notes",
]

REFERENCE_STYLE_TEMPLATES = {
    ("clamped_pinned", "thick"): {
        "linestyle": (0, (5.0, 2.2)),
        "linewidth": 1.15,
        "alpha": 0.38,
    },
    ("clamped_pinned", "thin"): {
        "linestyle": (0, (2.0, 2.4)),
        "linewidth": 1.15,
        "alpha": 0.42,
    },
    ("clamped_clamped", "thick"): {
        "linestyle": (0, (6.0, 2.0, 1.4, 2.0)),
        "linewidth": 1.25,
        "alpha": 0.42,
    },
    ("clamped_clamped", "thin"): {
        "linestyle": ":",
        "linewidth": 1.35,
        "alpha": 0.48,
    },
}


def number_text(value: float | int | str | None) -> str:
    return diagnostic_number_text(value, nonfinite_text="nan")


def resolve_output_dir(path: str | Path) -> Path:
    output_dir = Path(path)
    if output_dir.is_absolute():
        return output_dir
    return REPO_ROOT / output_dir


def output_stem(*, eta: float, epsilon: float, beta_deg: float) -> str:
    return (
        f"diagnostic_lambda_mu_eta_{number_token(eta)}"
        f"_eps{number_token(epsilon)}"
        f"_beta_{number_token(beta_deg)}deg_single_beam_refs"
    )


def validate_parameters(
    *,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
    num_modes: int,
    num_ref_modes: int,
) -> None:
    if float(epsilon) <= 0.0:
        raise ValueError("epsilon must be positive.")
    if not (-1.0 < float(eta) < 1.0):
        raise ValueError("eta must lie inside (-1, 1).")
    if int(num_modes) <= 0:
        raise ValueError("num_modes must be positive.")
    if int(num_ref_modes) <= 0:
        raise ValueError("num_ref_modes must be positive.")
    outside = mu_values[(mu_values <= -1.0) | (mu_values >= 1.0)]
    if outside.size:
        raise ValueError("mu values must stay inside (-1, 1) so both length factors are positive.")


def rod_state(mu: float, eta: float, rod_id: int) -> tuple[str, float, float]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    length_1 = 1.0 - float(mu)
    length_2 = 1.0 + float(mu)
    if rod_id == 1:
        length_label = "short" if length_1 <= length_2 else "long"
        thickness_label = "thick" if factors.tau1 >= factors.tau2 else "thin"
        return f"{length_label}_{thickness_label}", length_1, factors.tau1
    if rod_id == 2:
        length_label = "long" if length_2 >= length_1 else "short"
        thickness_label = "thick" if factors.tau2 >= factors.tau1 else "thin"
        return f"{length_label}_{thickness_label}", length_2, factors.tau2
    raise ValueError(f"unsupported rod_id {rod_id!r}")


def display_rod_label(label: str) -> str:
    return str(label).replace("_", "/")


def reference_style(boundary_condition: str, rod_label: str) -> dict[str, object]:
    thickness_label = "thick" if str(rod_label).endswith("_thick") else "thin"
    return REFERENCE_STYLE_TEMPLATES[(str(boundary_condition), thickness_label)]


def rod_labels_for_plot(mu_values: np.ndarray, eta: float) -> dict[int, str]:
    label_mu = float(mu_values[-1])
    return {rod_id: rod_state(label_mu, float(eta), rod_id)[0] for rod_id in (1, 2)}


def reference_alpha_roots(num_ref_modes: int) -> dict[str, np.ndarray]:
    return {
        "clamped_pinned": roots_clamped_supported(int(num_ref_modes)),
        "clamped_clamped": fixed_fixed_lambdas(int(num_ref_modes)),
    }


def build_reference_curves(
    *,
    mu_values: np.ndarray,
    eta: float,
    num_ref_modes: int,
) -> tuple[dict[tuple[str, int], np.ndarray], list[str]]:
    alpha_roots = reference_alpha_roots(num_ref_modes)
    references: dict[tuple[str, int], np.ndarray] = {}
    warnings: list[str] = []
    for boundary_condition, alphas in alpha_roots.items():
        for rod_id in (1, 2):
            values = np.full((int(num_ref_modes), len(mu_values)), np.nan, dtype=float)
            for col, mu in enumerate(mu_values):
                label, length_factor, tau = rod_state(float(mu), float(eta), rod_id)
                del label
                if not (np.isfinite(tau) and tau > 0.0):
                    raise ValueError(f"tau{rod_id} is not finite and positive at mu={float(mu):g}.")
                values[:, col] = alphas * np.sqrt(tau) / length_factor
            references[(boundary_condition, rod_id)] = values
            for mode_index in range(values.shape[0]):
                span = float(np.nanmax(values[mode_index]) - np.nanmin(values[mode_index]))
                scale = max(1.0, float(np.nanmax(np.abs(values[mode_index]))))
                if not np.isfinite(span) or span <= 1.0e-10 * scale:
                    warnings.append(
                        f"reference curve appears horizontal: {boundary_condition}, "
                        f"rod {rod_id}, mode {mode_index + 1}"
                    )
    return references, warnings


def solve_coupled_sorted_roots(
    *,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
    num_modes: int,
    root_scan_step: float,
    root_lmax0: float,
) -> tuple[np.ndarray, list[str], int]:
    roots_grid = np.full((len(mu_values), int(num_modes)), np.nan, dtype=float)
    warnings: list[str] = []
    missing_count = 0
    beta_rad = float(np.deg2rad(beta_deg))
    for row, mu in enumerate(mu_values):
        try:
            roots = find_first_n_roots_eta(
                beta_rad,
                float(mu),
                float(epsilon),
                float(eta),
                int(num_modes),
                Lmax0=float(root_lmax0),
                scan_step=float(root_scan_step),
            )
        except Exception as exc:  # pragma: no cover - diagnostic failure path
            roots = np.full(int(num_modes), np.nan, dtype=float)
            warnings.append(f"root solve failed at beta={float(beta_deg):g}, mu={float(mu):g}: {exc}")
        roots = np.asarray(roots, dtype=float)
        if roots.shape[0] < int(num_modes):
            padded = np.full(int(num_modes), np.nan, dtype=float)
            padded[: roots.shape[0]] = roots
            roots = padded
        roots_grid[row] = np.sort(roots[: int(num_modes)])
        missing_here = int(np.count_nonzero(~np.isfinite(roots_grid[row])))
        if missing_here:
            missing_count += missing_here
            warnings.append(
                f"missing {missing_here} sorted roots at beta={float(beta_deg):g}, mu={float(mu):g}"
            )
    return roots_grid, warnings, missing_count


def coupled_rows(
    *,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
    roots_grid: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for mu_idx, mu in enumerate(mu_values):
        for mode_idx in range(roots_grid.shape[1]):
            value = float(roots_grid[mu_idx, mode_idx])
            rows.append(
                {
                    "curve_group": "coupled_sorted",
                    "beta_deg": number_text(beta_deg),
                    "epsilon": number_text(epsilon),
                    "eta": number_text(eta),
                    "mu": number_text(float(mu)),
                    "mode_index": int(mode_idx + 1),
                    "Lambda": number_text(value),
                    "boundary_condition": "coupled",
                    "rod_id": "none",
                    "rod_label": "coupled",
                    "length_factor": "",
                    "tau": "",
                    "notes": "" if np.isfinite(value) else "missing_root",
                }
            )
    return rows


def reference_rows(
    *,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
    references: dict[tuple[str, int], np.ndarray],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for (boundary_condition, rod_id), values in references.items():
        for col, mu in enumerate(mu_values):
            label, length_factor, tau = rod_state(float(mu), float(eta), int(rod_id))
            for mode_idx in range(values.shape[0]):
                rows.append(
                    {
                        "curve_group": "single_rod_reference",
                        "beta_deg": number_text(beta_deg),
                        "epsilon": number_text(epsilon),
                        "eta": number_text(eta),
                        "mu": number_text(float(mu)),
                        "mode_index": int(mode_idx + 1),
                        "Lambda": number_text(float(values[mode_idx, col])),
                        "boundary_condition": boundary_condition,
                        "rod_id": int(rod_id),
                        "rod_label": label,
                        "length_factor": number_text(length_factor),
                        "tau": number_text(tau),
                        "notes": "Lambda_ref=alpha*sqrt(tau_i)/L_i",
                    }
                )
    return rows


def finite_max(values: np.ndarray, label: str) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        raise RuntimeError(f"No finite {label} values are available for plotting.")
    return float(np.max(finite))


def reference_max(references: dict[tuple[str, int], np.ndarray]) -> float:
    return finite_max(np.concatenate([values.ravel() for values in references.values()]), "reference")


def choose_y_max(
    *,
    coupled: np.ndarray,
    references: dict[tuple[str, int], np.ndarray],
    y_margin: float,
    ymax: float | None,
    zoom_coupled: bool,
) -> tuple[float, float, float]:
    if float(y_margin) <= 0.0:
        raise ValueError("y_margin must be positive.")
    coupled_max = finite_max(coupled, "coupled")
    refs_max = reference_max(references)
    if ymax is not None:
        y_max = float(ymax)
    elif bool(zoom_coupled):
        y_max = float(y_margin) * coupled_max
    else:
        y_max = float(y_margin) * max(coupled_max, refs_max)
    if not (np.isfinite(y_max) and y_max > 0.0):
        raise ValueError("computed y_max must be finite and positive.")
    return coupled_max, refs_max, y_max


def clip_reference_for_plot(values: np.ndarray, y_max: float, clip_ref_curves: bool) -> tuple[np.ndarray, int]:
    values = np.asarray(values, dtype=float)
    if not bool(clip_ref_curves):
        return values, 0
    mask = np.isfinite(values) & (values > float(y_max))
    return np.where(mask, np.nan, values), int(np.count_nonzero(mask))


def plot_beta_case(
    *,
    output_path: Path,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
    coupled: np.ndarray,
    references: dict[tuple[str, int], np.ndarray],
    show: bool,
    y_margin: float,
    ymax: float | None,
    zoom_coupled: bool,
    clip_ref_curves: bool,
    title_suffix: str = "",
) -> dict[str, float | int]:
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    coupled_max, refs_max, y_max = choose_y_max(
        coupled=coupled,
        references=references,
        y_margin=float(y_margin),
        ymax=ymax,
        zoom_coupled=bool(zoom_coupled),
    )
    fig, ax = plt.subplots(figsize=(10.4, 6.1))
    rod_labels = rod_labels_for_plot(mu_values, eta)

    for mode_idx in range(coupled.shape[1]):
        ax.plot(
            mu_values,
            coupled[:, mode_idx],
            color=colors[mode_idx % len(colors)],
            lw=2.1,
            ls="-",
            zorder=4,
        )

    clipped_points = 0
    for key, values in references.items():
        boundary_condition, rod_id = key
        style = reference_style(boundary_condition, rod_labels[int(rod_id)])
        plot_values, clipped_here = clip_reference_for_plot(
            values,
            y_max,
            clip_ref_curves=bool(clip_ref_curves),
        )
        clipped_points += clipped_here
        for mode_idx in range(values.shape[0]):
            ax.plot(
                mu_values,
                plot_values[mode_idx],
                color=colors[mode_idx % len(colors)],
                lw=style["linewidth"],
                ls=style["linestyle"],
                alpha=style["alpha"],
                zorder=2,
            )

    handles = [
        Line2D(
            [0],
            [0],
            color="black",
            lw=2.1,
            ls="-",
            label=f"solid: coupled sorted frequencies, first {coupled.shape[1]}",
        ),
    ]
    for key in (
        ("clamped_pinned", 1),
        ("clamped_pinned", 2),
        ("clamped_clamped", 1),
        ("clamped_clamped", 2),
    ):
        boundary_condition, rod_id = key
        rod_label = rod_labels[int(rod_id)]
        style = reference_style(boundary_condition, rod_label)
        bc_label = "CP/FP" if boundary_condition == "clamped_pinned" else "CC/FF"
        handles.append(
            Line2D(
                [0],
                [0],
                color="0.25",
                lw=style["linewidth"],
                ls=style["linestyle"],
                label=f"{bc_label} single rod {int(rod_id)}, {display_rod_label(rod_label)}",
            )
        )

    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    clipping_note = "; reference curves clipped to plotted y-range" if clip_ref_curves else ""
    rod_note = " and ".join(display_rod_label(rod_labels[rod_id]) for rod_id in (1, 2))
    ax.set_title(
        f"Diagnostic Lambda(mu), epsilon={epsilon:g}, eta={eta:g}, beta={beta_deg:g} deg\n"
        f"solid: coupled sorted frequencies; reference curves: single rods, {rod_note}"
        f"{clipping_note}{title_suffix}"
    )
    ax.set_ylim(0.0, y_max)
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(handles=handles, loc="upper left", fontsize=8.2, frameon=False, handlelength=3.8)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    return {
        "coupled_max": coupled_max,
        "reference_max": refs_max,
        "y_max": y_max,
        "clipped_reference_points": clipped_points,
    }


def run_case(
    *,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
    num_modes: int,
    num_ref_modes: int,
    output_dir: Path,
    show: bool,
    root_scan_step: float,
    root_lmax0: float,
    y_margin: float,
    ymax: float | None,
    zoom_coupled: bool,
    clip_ref_curves: bool,
    full_scale_refs: bool,
) -> dict[str, object]:
    references, reference_warnings = build_reference_curves(
        mu_values=mu_values,
        eta=eta,
        num_ref_modes=num_ref_modes,
    )
    coupled, root_warnings, missing_count = solve_coupled_sorted_roots(
        beta_deg=beta_deg,
        epsilon=epsilon,
        eta=eta,
        mu_values=mu_values,
        num_modes=num_modes,
        root_scan_step=root_scan_step,
        root_lmax0=root_lmax0,
    )
    finite_coupled = int(np.count_nonzero(np.isfinite(coupled)))
    if finite_coupled != coupled.size:
        root_warnings.append(
            f"finite sorted coupled root count {finite_coupled}/{coupled.size} at beta={float(beta_deg):g}"
        )

    stem = output_stem(eta=eta, epsilon=epsilon, beta_deg=beta_deg)
    png_path = output_dir / f"{stem}.png"
    fullscale_png_path = output_dir / f"{stem}_fullscale.png"
    csv_path = output_dir / f"{stem}.csv"
    rows = coupled_rows(
        beta_deg=beta_deg,
        epsilon=epsilon,
        eta=eta,
        mu_values=mu_values,
        roots_grid=coupled,
    )
    rows.extend(
        reference_rows(
            beta_deg=beta_deg,
            epsilon=epsilon,
            eta=eta,
            mu_values=mu_values,
            references=references,
        )
    )
    write_dict_rows_csv(csv_path, rows, CSV_FIELDNAMES)
    plot_stats = plot_beta_case(
        output_path=png_path,
        beta_deg=beta_deg,
        epsilon=epsilon,
        eta=eta,
        mu_values=mu_values,
        coupled=coupled,
        references=references,
        show=show,
        y_margin=y_margin,
        ymax=ymax,
        zoom_coupled=zoom_coupled,
        clip_ref_curves=clip_ref_curves,
    )
    fullscale_stats: dict[str, float | int] | None = None
    if bool(full_scale_refs):
        fullscale_stats = plot_beta_case(
            output_path=fullscale_png_path,
            beta_deg=beta_deg,
            epsilon=epsilon,
            eta=eta,
            mu_values=mu_values,
            coupled=coupled,
            references=references,
            show=show,
            y_margin=y_margin,
            ymax=None,
            zoom_coupled=False,
            clip_ref_curves=False,
            title_suffix="; full-scale reference view",
        )
    return {
        "beta_deg": beta_deg,
        "png": png_path,
        "fullscale_png": fullscale_png_path if bool(full_scale_refs) else None,
        "csv": csv_path,
        "missing_roots": missing_count,
        "warnings": root_warnings + reference_warnings,
        "rows": len(rows),
        "coupled_max": plot_stats["coupled_max"],
        "reference_max": plot_stats["reference_max"],
        "y_max": plot_stats["y_max"],
        "clipped_reference_points": plot_stats["clipped_reference_points"],
        "fullscale_y_max": None if fullscale_stats is None else fullscale_stats["y_max"],
    }


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic-only sorted Lambda(mu) plots for the thickness-mismatch EB "
            "model with CP/FP and CC/FF single-rod reference curves."
        )
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--betas", type=float, nargs="+", default=list(DEFAULT_BETAS))
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--num-modes", type=int, default=DEFAULT_NUM_MODES)
    parser.add_argument("--num-ref-modes", type=int, default=DEFAULT_NUM_REF_MODES)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Use a tiny grid for wiring checks and write default outputs under results/_smoke/.",
    )
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--root-lmax0", type=float, default=DEFAULT_ROOT_LMAX0)
    parser.add_argument("--ymax", type=float, default=None, help="Explicit y-axis maximum for default PNGs.")
    parser.add_argument("--y-margin", type=float, default=DEFAULT_Y_MARGIN)
    parser.add_argument(
        "--clip-ref-curves",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Clip plotted reference curves to the visible y-range; CSV data remain unclipped.",
    )
    parser.add_argument(
        "--zoom-coupled",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Use coupled sorted curves, not reference curves, to choose the default y-axis maximum.",
    )
    parser.add_argument(
        "--full-scale-refs",
        action="store_true",
        help="Also save *_fullscale.png plots with all reference curves visible.",
    )
    return parser.parse_args(argv)


def apply_smoke_defaults(args: argparse.Namespace) -> None:
    args.betas = [0.0, 5.0]
    args.mu_min = 0.0
    args.mu_max = 0.1
    args.mu_step = 0.1
    args.num_modes = min(int(args.num_modes), 3)
    args.num_ref_modes = min(int(args.num_ref_modes), 3)
    if Path(args.output_dir) == DEFAULT_OUTPUT_DIR:
        args.output_dir = SMOKE_OUTPUT_DIR


def main(argv: Sequence[str] | None = None) -> list[dict[str, object]]:
    args = parse_args(argv)
    if bool(args.smoke):
        apply_smoke_defaults(args)
    epsilon = float(args.epsilon)
    eta = float(args.eta)
    mu_values = inclusive_grid(args.mu_min, args.mu_max, args.mu_step, step_name="mu step")
    validate_parameters(
        epsilon=epsilon,
        eta=eta,
        mu_values=mu_values,
        num_modes=int(args.num_modes),
        num_ref_modes=int(args.num_ref_modes),
    )
    output_dir = resolve_output_dir(args.output_dir)

    print("diagnostic Lambda(mu) with single-beam references")
    if bool(args.smoke):
        print(f"smoke mode: outputs under {resolve_output_dir(args.output_dir)}")
    print(
        f"epsilon={epsilon:g}, eta={eta:g}, "
        f"mu={float(mu_values[0]):g}..{float(mu_values[-1]):g} step={float(args.mu_step):g}"
    )
    print("coupled curves: sorted frequencies; no descendant branch tracking")
    print("reference curves: Lambda_ref=alpha*sqrt(tau_i)/L_i, not horizontal lines")
    print(
        "default PNG y-scale: "
        f"{'coupled-spectrum zoom' if bool(args.zoom_coupled) else 'coupled+reference scale'}, "
        f"clip_ref_curves={bool(args.clip_ref_curves)}, y_margin={float(args.y_margin):g}, "
        f"ymax={'None' if args.ymax is None else f'{float(args.ymax):g}'}"
    )

    results: list[dict[str, object]] = []
    for beta_deg in args.betas:
        result = run_case(
            beta_deg=float(beta_deg),
            epsilon=epsilon,
            eta=eta,
            mu_values=mu_values,
            num_modes=int(args.num_modes),
            num_ref_modes=int(args.num_ref_modes),
            output_dir=output_dir,
            show=bool(args.show),
            root_scan_step=float(args.root_scan_step),
            root_lmax0=float(args.root_lmax0),
            y_margin=float(args.y_margin),
            ymax=None if args.ymax is None else float(args.ymax),
            zoom_coupled=bool(args.zoom_coupled),
            clip_ref_curves=bool(args.clip_ref_curves),
            full_scale_refs=bool(args.full_scale_refs),
        )
        results.append(result)
        warnings = result["warnings"]
        print(f"saved PNG: {result['png']}")
        if result["fullscale_png"] is not None:
            print(f"saved full-scale PNG: {result['fullscale_png']}")
        print(f"saved CSV: {result['csv']}")
        print(f"beta={float(beta_deg):g} deg, rows={result['rows']}, missing_roots={result['missing_roots']}")
        print(
            "visual scale: "
            f"coupled_max={float(result['coupled_max']):.12g}, "
            f"y_max={float(result['y_max']):.12g}, "
            f"max_reference={float(result['reference_max']):.12g}, "
            f"clipped_reference_points={int(result['clipped_reference_points'])}"
        )
        if result["fullscale_y_max"] is not None:
            print(f"full-scale y_max={float(result['fullscale_y_max']):.12g}")
        if warnings:
            for warning in warnings[:10]:
                print(f"WARNING: {warning}")
            if len(warnings) > 10:
                print(f"WARNING: plus {len(warnings) - 10} additional warnings")
        else:
            print(f"beta={float(beta_deg):g} deg: no root/reference warnings")

    total_missing = sum(int(result["missing_roots"]) for result in results)
    total_warnings = sum(len(result["warnings"]) for result in results)
    total_clipped = sum(int(result["clipped_reference_points"]) for result in results)
    print(f"total missing roots: {total_missing}")
    print(f"total warnings: {total_warnings}")
    print(f"total clipped reference points in default PNGs: {total_clipped}")
    print("protected areas untouched by this script: article workspace, FEM/Gmsh/CalculiX, old determinant, old solvers")
    return results


if __name__ == "__main__":
    main()
