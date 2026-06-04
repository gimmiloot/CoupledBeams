from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
from my_project.analytic.out_of_plane_fem_1d import solve_out_of_plane_fem_1d_modes  # noqa: E402
from my_project.analytic.solvers_out_of_plane import (  # noqa: E402
    find_first_n_roots_out_of_plane_with_warnings,
)


DEFAULT_EPSILON = 0.0025
DEFAULT_POISSON = 0.3
DEFAULT_MU_VALUES = (0.0, 0.3, 0.6)
DEFAULT_ETA_VALUES = (-0.5, 0.0, 0.5)
DEFAULT_BETA_MIN = 0.0
DEFAULT_BETA_MAX = 90.0
DEFAULT_BETA_STEP = 0.5
DEFAULT_N_ROOTS = 20
DEFAULT_LAMBDA_MAX = 25.0
DEFAULT_N_ELEMENTS_PER_ROD = 16
DEFAULT_ROOT_SCAN_STEP = 0.02
DEFAULT_OUTPUT_DIR = Path("results") / "out_of_plane_mode_character_beta"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "out_of_plane_mode_character_beta"

CHARACTER_COLORS = {
    "bending_dominated": "#1f77b4",
    "torsion_dominated": "#d62728",
    "mixed": "#9467bd",
    "unavailable": "0.55",
}


@dataclass(frozen=True)
class Args:
    epsilon: float
    poisson: float
    mu_values: tuple[float, ...]
    eta_values: tuple[float, ...]
    beta_min: float
    beta_max: float
    beta_step: float
    n_roots: int
    lambda_max: float
    n_elements_per_rod: int
    root_scan_step: float
    output_dir: Path
    smoke: bool


def _repo_output_dir(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


def _param_label(value: float) -> str:
    return f"{float(value):.1f}".replace("-", "m").replace(".", "p")


def _fmt(value: float) -> str:
    value_f = float(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16e}"


def _beta_grid(beta_min: float, beta_max: float, beta_step: float) -> np.ndarray:
    beta_min_f = float(beta_min)
    beta_max_f = float(beta_max)
    beta_step_f = float(beta_step)
    if not (isfinite(beta_min_f) and isfinite(beta_max_f) and isfinite(beta_step_f)):
        raise ValueError("beta grid values must be finite.")
    if beta_step_f <= 0.0:
        raise ValueError("beta_step must be positive.")
    if beta_max_f < beta_min_f:
        raise ValueError("beta_max must be greater than or equal to beta_min.")
    count = int(np.floor((beta_max_f - beta_min_f) / beta_step_f + 0.5)) + 1
    values = beta_min_f + beta_step_f * np.arange(count, dtype=float)
    if values[-1] < beta_max_f - 1e-10:
        values = np.append(values, beta_max_f)
    values[-1] = min(values[-1], beta_max_f)
    return values


def _parse_float_tuple(values: list[str], name: str) -> tuple[float, ...]:
    parsed = tuple(float(value) for value in values)
    if not parsed:
        raise ValueError(f"{name} must contain at least one value.")
    if not all(isfinite(value) for value in parsed):
        raise ValueError(f"{name} values must be finite.")
    return parsed


def classify_mode_character(bending_fraction: float, torsion_fraction: float) -> str:
    if not (isfinite(float(bending_fraction)) and isfinite(float(torsion_fraction))):
        return "unavailable"
    if float(bending_fraction) >= 0.8:
        return "bending_dominated"
    if float(torsion_fraction) >= 0.8:
        return "torsion_dominated"
    return "mixed"


def _validate_args(args: Args) -> None:
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not isfinite(args.poisson) or args.poisson <= -1.0:
        raise ValueError("poisson must be finite and greater than -1.")
    if args.n_roots <= 0:
        raise ValueError("n_roots must be positive.")
    if not isfinite(args.lambda_max) or args.lambda_max <= 0.0:
        raise ValueError("lambda_max must be positive and finite.")
    if args.n_elements_per_rod <= 0:
        raise ValueError("n_elements_per_rod must be positive.")
    if not isfinite(args.root_scan_step) or args.root_scan_step <= 0.0:
        raise ValueError("root_scan_step must be positive and finite.")
    _beta_grid(args.beta_min, args.beta_max, args.beta_step)
    for mu in args.mu_values:
        if not (-1.0 < float(mu) < 1.0):
            raise ValueError("all mu values must lie inside (-1, 1).")
    for mu in args.mu_values:
        for eta in args.eta_values:
            thickness_mismatch_factors(mu, eta)


def parse_args(argv: list[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(description="Out-of-plane mode-character Lambda(beta) diagnostics.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--mu-values", nargs="+", default=[str(value) for value in DEFAULT_MU_VALUES])
    parser.add_argument("--eta-values", nargs="+", default=[str(value) for value in DEFAULT_ETA_VALUES])
    parser.add_argument("--beta-min", type=float, default=DEFAULT_BETA_MIN)
    parser.add_argument("--beta-max", type=float, default=DEFAULT_BETA_MAX)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--n-roots", type=int, default=DEFAULT_N_ROOTS)
    parser.add_argument("--lambda-max", type=float, default=DEFAULT_LAMBDA_MAX)
    parser.add_argument("--n-elements-per-rod", type=int, default=DEFAULT_N_ELEMENTS_PER_ROD)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)

    mu_values = _parse_float_tuple(ns.mu_values, "mu-values")
    eta_values = _parse_float_tuple(ns.eta_values, "eta-values")
    beta_min = float(ns.beta_min)
    beta_max = float(ns.beta_max)
    beta_step = float(ns.beta_step)
    n_roots = int(ns.n_roots)
    output_dir = _repo_output_dir(Path(ns.output_dir))
    if ns.smoke:
        mu_values = (0.0, 0.3)
        eta_values = (0.0,)
        beta_min = 0.0
        beta_max = 90.0
        beta_step = 45.0
        n_roots = 6
        output_dir = _repo_output_dir(SMOKE_OUTPUT_DIR)

    args = Args(
        epsilon=float(ns.epsilon),
        poisson=float(ns.poisson),
        mu_values=mu_values,
        eta_values=eta_values,
        beta_min=beta_min,
        beta_max=beta_max,
        beta_step=beta_step,
        n_roots=n_roots,
        lambda_max=float(ns.lambda_max),
        n_elements_per_rod=int(ns.n_elements_per_rod),
        root_scan_step=float(ns.root_scan_step),
        output_dir=output_dir,
        smoke=bool(ns.smoke),
    )
    _validate_args(args)
    return args


def _rel_error(estimate: float, reference: float) -> float:
    if not (isfinite(float(estimate)) and isfinite(float(reference))) or abs(float(reference)) <= 0.0:
        return float("nan")
    return abs(float(estimate) - float(reference)) / abs(float(reference))


def compute_rows(args: Args, beta_values: np.ndarray) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for eta in args.eta_values:
        for mu in args.mu_values:
            for beta_deg in beta_values:
                beta_rad = float(np.deg2rad(beta_deg))
                analytic_result = find_first_n_roots_out_of_plane_with_warnings(
                    beta=beta_rad,
                    mu=float(mu),
                    epsilon=args.epsilon,
                    eta=float(eta),
                    poisson=args.poisson,
                    n_roots=args.n_roots,
                    lambda_max=args.lambda_max,
                    scan_step=args.root_scan_step,
                )
                analytic = np.asarray(analytic_result.roots, dtype=float)
                fem_warning = "none"
                try:
                    fem = solve_out_of_plane_fem_1d_modes(
                        beta=beta_rad,
                        mu=float(mu),
                        epsilon=args.epsilon,
                        eta=float(eta),
                        poisson=args.poisson,
                        n_elements_per_rod=args.n_elements_per_rod,
                        n_modes=args.n_roots,
                    )
                except Exception as exc:  # pragma: no cover - diagnostic reporting path.
                    fem = None
                    fem_warning = f"FEM mode solve failed: {exc}"

                for index in range(args.n_roots):
                    lambda_analytic = float(analytic[index]) if index < len(analytic) else float("nan")
                    if fem is None:
                        lambda_fem = float("nan")
                        bending_fraction = float("nan")
                        torsion_fraction = float("nan")
                    else:
                        lambda_fem = float(fem.lambdas[index])
                        bending_fraction = float(fem.bending_energy_fraction[index])
                        torsion_fraction = float(fem.torsion_energy_fraction[index])
                    character = classify_mode_character(bending_fraction, torsion_fraction)
                    notes = []
                    if index >= len(analytic):
                        notes.append("missing analytic root")
                    if analytic_result.warnings:
                        notes.append("analytic root warning")
                    if fem_warning != "none":
                        notes.append("FEM warning")
                    if not notes:
                        notes.append("ok")
                    rows.append(
                        {
                            "beta_deg": float(beta_deg),
                            "beta_rad": beta_rad,
                            "mu": float(mu),
                            "eta": float(eta),
                            "epsilon": args.epsilon,
                            "poisson": args.poisson,
                            "sorted_index": index + 1,
                            "Lambda_analytic": lambda_analytic,
                            "Lambda_fem_1d": lambda_fem,
                            "abs_error": abs(lambda_fem - lambda_analytic)
                            if isfinite(lambda_fem) and isfinite(lambda_analytic)
                            else float("nan"),
                            "rel_error": _rel_error(lambda_fem, lambda_analytic),
                            "bending_energy_fraction": bending_fraction,
                            "torsion_energy_fraction": torsion_fraction,
                            "mode_character": character,
                            "root_solver_warning": "; ".join(analytic_result.warnings)
                            if analytic_result.warnings
                            else "none",
                            "fem_warning": fem_warning,
                            "notes": "; ".join(notes),
                        }
                    )
    return rows


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "beta_deg",
        "beta_rad",
        "mu",
        "eta",
        "epsilon",
        "poisson",
        "sorted_index",
        "Lambda_analytic",
        "Lambda_fem_1d",
        "abs_error",
        "rel_error",
        "bending_energy_fraction",
        "torsion_energy_fraction",
        "mode_character",
        "root_solver_warning",
        "fem_warning",
        "notes",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in list(out.items()):
                if isinstance(value, float):
                    out[key] = _fmt(value)
            writer.writerow(out)


def rows_for(rows: list[dict[str, object]], *, eta: float, mu: float, sorted_index: int) -> list[dict[str, object]]:
    return [
        row
        for row in rows
        if abs(float(row["eta"]) - float(eta)) <= 1e-12
        and abs(float(row["mu"]) - float(mu)) <= 1e-12
        and int(row["sorted_index"]) == int(sorted_index)
    ]


def plot_mode_character(args: Args, rows: list[dict[str, object]], beta_values: np.ndarray, eta: float) -> Path:
    fig, axes = plt.subplots(len(args.mu_values), 1, figsize=(9.4, 3.0 * len(args.mu_values)), sharex=True)
    axes_array = np.atleast_1d(axes)
    for ax, mu in zip(axes_array, args.mu_values):
        for sorted_index in range(1, args.n_roots + 1):
            selected = sorted(rows_for(rows, eta=eta, mu=mu, sorted_index=sorted_index), key=lambda row: float(row["beta_deg"]))
            lambdas = np.asarray([float(row["Lambda_analytic"]) for row in selected], dtype=float)
            ax.plot(beta_values, lambdas, color="0.78", linewidth=0.7, zorder=1)
            for character, color in CHARACTER_COLORS.items():
                mask = np.asarray([row["mode_character"] == character for row in selected], dtype=bool)
                if np.any(mask):
                    ax.scatter(beta_values[mask], lambdas[mask], s=8, color=color, linewidths=0, zorder=2)
        ax.set_title(f"mu={mu:g}")
        ax.set_ylabel("Lambda")
        ax.grid(True, color="0.9", linewidth=0.5)
    axes_array[-1].set_xlabel("beta [deg]")
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=color, markersize=5, label=label)
        for label, color in CHARACTER_COLORS.items()
        if label != "unavailable"
    ]
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=8)
    fig.suptitle(
        "out-of-plane EB + Saint-Venant torsion mode character\n"
        f"epsilon={args.epsilon:g}, poisson={args.poisson:g}, eta={eta:g}, sorted frequencies",
        fontsize=12,
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.88, 0.93))
    path = args.output_dir / f"out_of_plane_lambda_beta_mode_character_eta_{_param_label(eta)}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_torsion_fraction(args: Args, rows: list[dict[str, object]], beta_values: np.ndarray, eta: float) -> Path:
    fig, axes = plt.subplots(len(args.mu_values), 1, figsize=(9.4, 3.0 * len(args.mu_values)), sharex=True, sharey=True)
    axes_array = np.atleast_1d(axes)
    colors = plt.cm.tab20(np.linspace(0.0, 1.0, max(args.n_roots, 1)))
    for ax, mu in zip(axes_array, args.mu_values):
        for sorted_index in range(1, args.n_roots + 1):
            selected = sorted(rows_for(rows, eta=eta, mu=mu, sorted_index=sorted_index), key=lambda row: float(row["beta_deg"]))
            values = np.asarray([float(row["torsion_energy_fraction"]) for row in selected], dtype=float)
            ax.plot(beta_values, values, color=colors[sorted_index - 1], linewidth=1.0, label=f"{sorted_index}")
        ax.axhline(0.8, color="0.35", linestyle="--", linewidth=0.8)
        ax.set_title(f"mu={mu:g}")
        ax.set_ylabel("torsion energy fraction")
        ax.grid(True, color="0.9", linewidth=0.5)
    axes_array[-1].set_xlabel("beta [deg]")
    axes_array[0].set_ylim(-0.03, 1.03)
    fig.suptitle(
        f"out-of-plane torsion fraction, eta={eta:g}, epsilon={args.epsilon:g}",
        fontsize=12,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    path = args.output_dir / f"out_of_plane_torsion_fraction_eta_{_param_label(eta)}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def character_counts(rows: list[dict[str, object]]) -> dict[str, int]:
    counts = {key: 0 for key in CHARACTER_COLORS}
    for row in rows:
        counts[str(row["mode_character"])] = counts.get(str(row["mode_character"]), 0) + 1
    return counts


def first_character_row(rows: list[dict[str, object]], character: str) -> dict[str, object] | None:
    candidates = [row for row in rows if row["mode_character"] == character and isfinite(float(row["Lambda_analytic"]))]
    if not candidates:
        return None
    return min(candidates, key=lambda row: (float(row["Lambda_analytic"]), int(row["sorted_index"]), float(row["beta_deg"])))


def beta_torsion_trend(rows: list[dict[str, object]], beta_values: np.ndarray) -> tuple[float, float, float]:
    beta0 = float(beta_values[0])
    beta1 = float(beta_values[-1])
    rows0 = [row for row in rows if abs(float(row["beta_deg"]) - beta0) <= 1e-12]
    rows1 = [row for row in rows if abs(float(row["beta_deg"]) - beta1) <= 1e-12]
    mean0 = float(np.nanmean([float(row["torsion_energy_fraction"]) for row in rows0])) if rows0 else float("nan")
    mean1 = float(np.nanmean([float(row["torsion_energy_fraction"]) for row in rows1])) if rows1 else float("nan")
    return mean0, mean1, mean1 - mean0


def write_report(args: Args, rows: list[dict[str, object]], beta_values: np.ndarray, output_paths: list[Path]) -> Path:
    report_path = args.output_dir / "out_of_plane_mode_character_beta_report.md"
    counts = character_counts(rows)
    first_eight = [row for row in rows if int(row["sorted_index"]) <= 8]
    first_eight_counts = character_counts(first_eight)
    first_torsion = first_character_row(rows, "torsion_dominated")
    first_mixed = first_character_row(rows, "mixed")
    mean0, mean1, delta = beta_torsion_trend(rows, beta_values)
    root_warning_cases = {
        (row["eta"], row["mu"], row["beta_deg"])
        for row in rows
        if row["root_solver_warning"] != "none"
    }
    fem_warning_cases = {
        (row["eta"], row["mu"], row["beta_deg"])
        for row in rows
        if row["fem_warning"] != "none"
    }

    lines = [
        "# Out-of-Plane Mode-Character Diagnostics",
        "",
        "This diagnostic classifies sorted out-of-plane EB + Saint-Venant",
        "torsion roots by matching them to the independent 1D FEM modes and",
        "using bending/torsion stiffness-energy fractions.",
        "",
        "This is mode-character diagnostics, not descendant tracking and not",
        "crossing verification.",
        "",
        "## Parameters",
        "",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        f"- beta grid: {float(beta_values[0]):g}..{float(beta_values[-1]):g} deg, "
        f"{len(beta_values)} points, nominal step {args.beta_step:g} deg",
        f"- mu values: {', '.join(f'{value:g}' for value in args.mu_values)}",
        f"- eta values: {', '.join(f'{value:g}' for value in args.eta_values)}",
        f"- sorted roots/modes requested: {args.n_roots}",
        f"- lambda_max for analytic roots: {args.lambda_max:g}",
        f"- 1D FEM elements per rod: {args.n_elements_per_rod}",
        "",
        "## Character Counts",
        "",
        f"- all rows: {counts}",
        f"- first eight sorted indices: {first_eight_counts}",
        "",
        "## Observations",
        "",
        "- For epsilon=0.0025, the first eight sorted indices remain dominated",
        f"  by bending in {first_eight_counts.get('bending_dominated', 0)}/{len(first_eight)} rows.",
    ]
    if first_torsion is None:
        lines.append("- No torsion-dominated rows were found on this grid.")
    else:
        lines.append(
            "- First torsion-dominated row by Lambda: "
            f"eta={float(first_torsion['eta']):g}, mu={float(first_torsion['mu']):g}, "
            f"beta={float(first_torsion['beta_deg']):g} deg, sorted "
            f"{int(first_torsion['sorted_index'])}, Lambda={float(first_torsion['Lambda_analytic']):.8g}, "
            f"torsion fraction={float(first_torsion['torsion_energy_fraction']):.3f}."
        )
    if first_mixed is None:
        lines.append("- No mixed rows were found on this grid.")
    else:
        lines.append(
            "- First mixed row by Lambda: "
            f"eta={float(first_mixed['eta']):g}, mu={float(first_mixed['mu']):g}, "
            f"beta={float(first_mixed['beta_deg']):g} deg, sorted "
            f"{int(first_mixed['sorted_index'])}, Lambda={float(first_mixed['Lambda_analytic']):.8g}, "
            f"torsion fraction={float(first_mixed['torsion_energy_fraction']):.3f}."
        )
    lines.extend(
        [
            f"- Mean torsion fraction at beta={float(beta_values[0]):g} deg: {mean0:.6g}",
            f"- Mean torsion fraction at beta={float(beta_values[-1]):g} deg: {mean1:.6g}",
            f"- Change over beta range: {delta:.6g}",
            "",
            "## Warnings",
            "",
            f"- analytic root warning cases: {len(root_warning_cases)}",
            f"- FEM warning cases: {len(fem_warning_cases)}",
            "",
            "## Generated Files",
            "",
        ]
    )
    for path in [*output_paths, report_path]:
        lines.append(f"- `{path.relative_to(REPO_ROOT)}`")
    lines.append("")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def run(args: Args) -> dict[str, object]:
    beta_values = _beta_grid(args.beta_min, args.beta_max, args.beta_step)
    rows = compute_rows(args, beta_values)
    output_paths: list[Path] = []
    csv_path = args.output_dir / "out_of_plane_mode_character_beta.csv"
    write_csv(csv_path, rows)
    output_paths.append(csv_path)
    for eta in args.eta_values:
        output_paths.append(plot_mode_character(args, rows, beta_values, eta))
        output_paths.append(plot_torsion_fraction(args, rows, beta_values, eta))
    report_path = write_report(args, rows, beta_values, output_paths)
    output_paths.append(report_path)

    counts = character_counts(rows)
    first_torsion = first_character_row(rows, "torsion_dominated")
    first_mixed = first_character_row(rows, "mixed")
    warning_cases = {
        (row["eta"], row["mu"], row["beta_deg"])
        for row in rows
        if row["root_solver_warning"] != "none" or row["fem_warning"] != "none"
    }
    print(f"saved mode-character CSV: {csv_path}")
    print(f"saved report: {report_path}")
    print(f"generated PNG count: {sum(1 for path in output_paths if path.suffix.lower() == '.png')}")
    print(f"beta grid: {float(beta_values[0]):g}..{float(beta_values[-1]):g} deg, {len(beta_values)} points")
    print(f"mode-character counts: {counts}")
    print(f"warning cases: {len(warning_cases)}")
    if first_torsion is not None:
        print(
            "first torsion-dominated row: "
            f"eta={float(first_torsion['eta']):g}, mu={float(first_torsion['mu']):g}, "
            f"beta={float(first_torsion['beta_deg']):g}, sorted={int(first_torsion['sorted_index'])}, "
            f"Lambda={float(first_torsion['Lambda_analytic']):.8g}"
        )
    if first_mixed is not None:
        print(
            "first mixed row: "
            f"eta={float(first_mixed['eta']):g}, mu={float(first_mixed['mu']):g}, "
            f"beta={float(first_mixed['beta_deg']):g}, sorted={int(first_mixed['sorted_index'])}, "
            f"Lambda={float(first_mixed['Lambda_analytic']):.8g}"
        )
    return {
        "rows": rows,
        "output_paths": output_paths,
        "counts": counts,
        "first_torsion": first_torsion,
        "first_mixed": first_mixed,
        "warning_cases": len(warning_cases),
    }


def main(argv: list[str] | None = None) -> dict[str, object]:
    return run(parse_args(argv))


if __name__ == "__main__":
    main()
