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

from my_project.analytic.formulas_out_of_plane import torsion_roots_uniform_eta0_beta0  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
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
DEFAULT_N_ROOTS_LOW = 8
DEFAULT_N_ROOTS_EXTENDED = 20
DEFAULT_LAMBDA_MAX_EXTENDED = 25.0
DEFAULT_ROOT_SCAN_STEP = 0.02
DEFAULT_OUTPUT_DIR = Path("results") / "out_of_plane_lambda_beta_mu_eta"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "out_of_plane_lambda_beta_mu_eta"


@dataclass(frozen=True)
class Args:
    epsilon: float
    poisson: float
    mu_values: tuple[float, ...]
    eta_values: tuple[float, ...]
    beta_min: float
    beta_max: float
    beta_step: float
    n_roots_low: int
    n_roots_extended: int
    lambda_max_extended: float
    root_scan_step: float
    output_dir: Path
    smoke: bool


@dataclass(frozen=True)
class RootCase:
    beta_deg: float
    beta_rad: float
    mu: float
    eta: float
    roots: tuple[float, ...]
    warnings: tuple[str, ...]
    notes: tuple[str, ...]


def _float_label(value: float) -> str:
    text = f"{float(value):.6g}".replace("-", "m").replace(".", "p")
    return text


def _fmt(value: float) -> str:
    value_f = float(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16e}"


def _csv_write(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            for key, value in list(out.items()):
                if isinstance(value, float):
                    out[key] = _fmt(value)
            writer.writerow(out)


def _beta_grid(beta_min: float, beta_max: float, beta_step: float) -> np.ndarray:
    beta_min_f = float(beta_min)
    beta_max_f = float(beta_max)
    beta_step_f = float(beta_step)
    if not (isfinite(beta_min_f) and isfinite(beta_max_f) and isfinite(beta_step_f)):
        raise ValueError("beta range values must be finite.")
    if beta_step_f <= 0.0:
        raise ValueError("beta_step must be positive.")
    if beta_max_f < beta_min_f:
        raise ValueError("beta_max must be greater than or equal to beta_min.")
    count = int(np.floor((beta_max_f - beta_min_f) / beta_step_f + 0.5)) + 1
    values = beta_min_f + beta_step_f * np.arange(count, dtype=float)
    if values[-1] < beta_max_f - 1e-10:
        values = np.append(values, beta_max_f)
    values[-1] = min(values[-1], beta_max_f)
    if not np.all(np.isfinite(values)):
        raise ValueError("beta grid contains non-finite values.")
    return values


def _parse_float_tuple(values: list[str], name: str) -> tuple[float, ...]:
    parsed = tuple(float(value) for value in values)
    if not parsed:
        raise ValueError(f"{name} must contain at least one value.")
    if not all(isfinite(value) for value in parsed):
        raise ValueError(f"{name} values must be finite.")
    return parsed


def _repo_output_dir(path: Path) -> Path:
    path_obj = Path(path)
    if path_obj.is_absolute():
        return path_obj
    return REPO_ROOT / path_obj


def _validate_args(args: Args) -> None:
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not isfinite(args.poisson) or args.poisson <= -1.0:
        raise ValueError("poisson must be finite and greater than -1.")
    if args.n_roots_low <= 0:
        raise ValueError("n_roots_low must be positive.")
    if args.n_roots_extended <= 0:
        raise ValueError("n_roots_extended must be positive.")
    if args.n_roots_extended < args.n_roots_low:
        raise ValueError("n_roots_extended must be at least n_roots_low.")
    if not isfinite(args.lambda_max_extended) or args.lambda_max_extended <= 0.0:
        raise ValueError("lambda_max_extended must be positive and finite.")
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
    parser = argparse.ArgumentParser(
        description="Plot sorted out-of-plane EB + Saint-Venant torsion Lambda(beta) diagnostics."
    )
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--poisson", type=float, default=DEFAULT_POISSON)
    parser.add_argument("--mu-values", nargs="+", default=[str(value) for value in DEFAULT_MU_VALUES])
    parser.add_argument("--eta-values", nargs="+", default=[str(value) for value in DEFAULT_ETA_VALUES])
    parser.add_argument("--beta-min", type=float, default=DEFAULT_BETA_MIN)
    parser.add_argument("--beta-max", type=float, default=DEFAULT_BETA_MAX)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--n-roots-low", type=int, default=DEFAULT_N_ROOTS_LOW)
    parser.add_argument("--n-roots-extended", type=int, default=DEFAULT_N_ROOTS_EXTENDED)
    parser.add_argument("--lambda-max-extended", type=float, default=DEFAULT_LAMBDA_MAX_EXTENDED)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)

    mu_values = _parse_float_tuple(ns.mu_values, "mu-values")
    eta_values = _parse_float_tuple(ns.eta_values, "eta-values")
    output_dir = _repo_output_dir(Path(ns.output_dir))
    beta_min = float(ns.beta_min)
    beta_max = float(ns.beta_max)
    beta_step = float(ns.beta_step)
    n_roots_low = int(ns.n_roots_low)
    n_roots_extended = int(ns.n_roots_extended)
    if ns.smoke:
        mu_values = (0.0, 0.3)
        eta_values = (0.0,)
        beta_min = 0.0
        beta_max = 90.0
        beta_step = 45.0
        n_roots_low = 4
        n_roots_extended = 8
        output_dir = _repo_output_dir(SMOKE_OUTPUT_DIR)

    args = Args(
        epsilon=float(ns.epsilon),
        poisson=float(ns.poisson),
        mu_values=mu_values,
        eta_values=eta_values,
        beta_min=beta_min,
        beta_max=beta_max,
        beta_step=beta_step,
        n_roots_low=n_roots_low,
        n_roots_extended=n_roots_extended,
        lambda_max_extended=float(ns.lambda_max_extended),
        root_scan_step=float(ns.root_scan_step),
        output_dir=output_dir,
        smoke=bool(ns.smoke),
    )
    _validate_args(args)
    return args


def compute_roots(args: Args, beta_values: np.ndarray) -> dict[tuple[float, float, float], RootCase]:
    cases: dict[tuple[float, float, float], RootCase] = {}
    for eta in args.eta_values:
        for mu in args.mu_values:
            for beta_deg in beta_values:
                beta_rad = float(np.deg2rad(beta_deg))
                result = find_first_n_roots_out_of_plane_with_warnings(
                    beta=beta_rad,
                    mu=float(mu),
                    epsilon=args.epsilon,
                    eta=float(eta),
                    poisson=args.poisson,
                    n_roots=args.n_roots_extended,
                    lambda_max=args.lambda_max_extended,
                    scan_step=args.root_scan_step,
                )
                notes: list[str] = []
                roots = tuple(float(root) for root in result.roots)
                if len(roots) < args.n_roots_extended:
                    notes.append("missing extended roots")
                if len(roots) < args.n_roots_low:
                    notes.append("missing low roots")
                if roots and any(not isfinite(root) for root in roots):
                    notes.append("non-finite root")
                if any(roots[i + 1] < roots[i] for i in range(len(roots) - 1)):
                    notes.append("roots not nondecreasing")
                if not notes:
                    notes.append("ok")
                cases[(float(eta), float(mu), float(beta_deg))] = RootCase(
                    beta_deg=float(beta_deg),
                    beta_rad=beta_rad,
                    mu=float(mu),
                    eta=float(eta),
                    roots=roots,
                    warnings=tuple(result.warnings),
                    notes=tuple(notes),
                )
    return cases


def root_at(case: RootCase, index_zero: int) -> float:
    if index_zero < len(case.roots):
        return float(case.roots[index_zero])
    return float("nan")


def warning_text(case: RootCase) -> str:
    return "; ".join(case.warnings) if case.warnings else "none"


def notes_text(case: RootCase) -> str:
    return "; ".join(case.notes) if case.notes else "ok"


def build_root_rows(args: Args, cases: dict[tuple[float, float, float], RootCase]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for key in sorted(cases):
        case = cases[key]
        for window, count in (("low", args.n_roots_low), ("extended", args.n_roots_extended)):
            for index_zero in range(count):
                root = root_at(case, index_zero)
                row_notes = [notes_text(case)]
                if not isfinite(root):
                    row_notes.append("missing root")
                rows.append(
                    {
                        "beta_deg": case.beta_deg,
                        "beta_rad": case.beta_rad,
                        "mu": case.mu,
                        "eta": case.eta,
                        "epsilon": args.epsilon,
                        "poisson": args.poisson,
                        "sorted_index": index_zero + 1,
                        "Lambda": root,
                        "spectrum_window": window,
                        "root_solver_warning": warning_text(case),
                        "notes": "; ".join(row_notes),
                    }
                )
    return rows


def case_series(
    cases: dict[tuple[float, float, float], RootCase],
    *,
    eta: float,
    mu: float,
    beta_values: np.ndarray,
    n_roots: int,
) -> np.ndarray:
    values = np.full((n_roots, len(beta_values)), np.nan, dtype=float)
    for beta_index, beta_deg in enumerate(beta_values):
        case = cases[(float(eta), float(mu), float(beta_deg))]
        for root_index in range(n_roots):
            values[root_index, beta_index] = root_at(case, root_index)
    return values


def set_common_y_limit(axes: np.ndarray, data: list[np.ndarray], *, upper: float | None = None) -> None:
    finite_arrays = [array[np.isfinite(array)] for array in data if np.any(np.isfinite(array))]
    if not finite_arrays:
        return
    finite_values = np.concatenate(finite_arrays)
    if len(finite_values) == 0:
        return
    y_min = max(0.0, float(np.nanmin(finite_values)) - 0.05 * float(np.nanmax(finite_values)))
    y_max = float(np.nanmax(finite_values)) * 1.05
    if upper is not None:
        y_max = max(y_max, float(upper))
    for ax in np.ravel(axes):
        ax.set_ylim(y_min, y_max)


def plot_eta_figure(
    args: Args,
    cases: dict[tuple[float, float, float], RootCase],
    beta_values: np.ndarray,
    *,
    eta: float,
    window: str,
) -> Path:
    n_roots = args.n_roots_low if window == "low" else args.n_roots_extended
    fig, axes = plt.subplots(len(args.mu_values), 1, figsize=(9.6, 3.1 * len(args.mu_values)), sharex=True)
    axes_array = np.atleast_1d(axes)
    data_arrays: list[np.ndarray] = []
    colors = plt.cm.tab20(np.linspace(0.0, 1.0, max(n_roots, 1)))
    for ax, mu in zip(axes_array, args.mu_values):
        values = case_series(cases, eta=eta, mu=mu, beta_values=beta_values, n_roots=n_roots)
        data_arrays.append(values)
        for root_index in range(n_roots):
            ax.plot(
                beta_values,
                values[root_index],
                color=colors[root_index],
                linewidth=1.2,
                label=f"sorted {root_index + 1}",
            )
        ax.set_ylabel("Lambda")
        ax.set_title(f"mu={mu:g}")
        ax.grid(True, color="0.88", linewidth=0.6)
    torsion_root = torsion_roots_uniform_eta0_beta0(1, args.epsilon, args.poisson)
    if window == "extended":
        for ax in axes_array:
            ax.axhline(torsion_root, color="0.35", linestyle="--", linewidth=0.9, label="eta=0 torsion root")
        set_common_y_limit(axes_array, data_arrays, upper=torsion_root * 1.04)
    else:
        set_common_y_limit(axes_array, data_arrays)
    axes_array[-1].set_xlabel("beta [deg]")
    handles, labels = axes_array[0].get_legend_handles_labels()
    max_legend = min(len(handles), 10 if window == "low" else 12)
    fig.legend(handles[:max_legend], labels[:max_legend], loc="upper right", frameon=False, fontsize=8)
    fig.suptitle(
        "out-of-plane EB + Saint-Venant torsion\n"
        f"epsilon={args.epsilon:g}, poisson={args.poisson:g}, eta={eta:g}, sorted frequencies",
        fontsize=12,
    )
    fig.tight_layout(rect=(0.0, 0.0, 0.88, 0.93))
    path = args.output_dir / f"out_of_plane_lambda_beta_{window}_eta_{_float_label(eta)}.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def plot_overview(
    args: Args,
    cases: dict[tuple[float, float, float], RootCase],
    beta_values: np.ndarray,
) -> Path:
    fig, axes = plt.subplots(
        len(args.eta_values),
        len(args.mu_values),
        figsize=(4.2 * len(args.mu_values), 3.0 * len(args.eta_values)),
        sharex=True,
        sharey=True,
    )
    axes_array = np.asarray(axes)
    if axes_array.ndim == 1:
        if len(args.eta_values) == 1:
            axes_array = axes_array.reshape(1, -1)
        else:
            axes_array = axes_array.reshape(-1, 1)
    colors = plt.cm.tab10(np.linspace(0.0, 1.0, max(args.n_roots_low, 1)))
    data_arrays: list[np.ndarray] = []
    for row, eta in enumerate(args.eta_values):
        for col, mu in enumerate(args.mu_values):
            ax = axes_array[row, col]
            values = case_series(cases, eta=eta, mu=mu, beta_values=beta_values, n_roots=args.n_roots_low)
            data_arrays.append(values)
            for root_index in range(args.n_roots_low):
                ax.plot(beta_values, values[root_index], color=colors[root_index], linewidth=1.0)
            ax.set_title(f"eta={eta:g}, mu={mu:g}", fontsize=9)
            ax.grid(True, color="0.9", linewidth=0.5)
            if row == len(args.eta_values) - 1:
                ax.set_xlabel("beta [deg]")
            if col == 0:
                ax.set_ylabel("Lambda")
    set_common_y_limit(axes_array, data_arrays)
    fig.suptitle("out-of-plane EB + Saint-Venant torsion, first sorted roots", fontsize=12)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
    path = args.output_dir / "out_of_plane_lambda_beta_mu_eta_overview.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def build_gap_summary(
    args: Args,
    cases: dict[tuple[float, float, float], RootCase],
    beta_values: np.ndarray,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for eta in args.eta_values:
        for mu in args.mu_values:
            best: dict[str, object] | None = None
            missing = False
            for beta_deg in beta_values:
                case = cases[(float(eta), float(mu), float(beta_deg))]
                values = np.asarray([root_at(case, index) for index in range(args.n_roots_low)], dtype=float)
                if not np.all(np.isfinite(values)):
                    missing = True
                    continue
                gaps = np.diff(values)
                if len(gaps) == 0:
                    continue
                local_index = int(np.nanargmin(gaps))
                row = {
                    "mu": float(mu),
                    "eta": float(eta),
                    "epsilon": args.epsilon,
                    "poisson": args.poisson,
                    "g_min": float(gaps[local_index]),
                    "beta_at_g_min_deg": float(beta_deg),
                    "pair_at_g_min": f"{local_index + 1}-{local_index + 2}",
                    "Lambda_lower": float(values[local_index]),
                    "Lambda_upper": float(values[local_index + 1]),
                    "notes": "close approach candidate; sorted-frequency diagnostic only",
                }
                if best is None or float(row["g_min"]) < float(best["g_min"]):
                    best = row
            if best is None:
                best = {
                    "mu": float(mu),
                    "eta": float(eta),
                    "epsilon": args.epsilon,
                    "poisson": args.poisson,
                    "g_min": float("nan"),
                    "beta_at_g_min_deg": float("nan"),
                    "pair_at_g_min": "",
                    "Lambda_lower": float("nan"),
                    "Lambda_upper": float("nan"),
                    "notes": "missing roots; no gap summary",
                }
            elif missing:
                best = dict(best)
                best["notes"] = str(best["notes"]) + "; some beta rows had missing roots"
            rows.append(best)
    return rows


def plot_gap_summary(args: Args, rows: list[dict[str, object]]) -> Path:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for eta in args.eta_values:
        selected = [row for row in rows if abs(float(row["eta"]) - float(eta)) <= 1e-12]
        selected = sorted(selected, key=lambda row: float(row["mu"]))
        ax.plot(
            [float(row["mu"]) for row in selected],
            [float(row["g_min"]) for row in selected],
            marker="o",
            linewidth=1.4,
            label=f"eta={eta:g}",
        )
    ax.set_xlabel("mu")
    ax.set_ylabel("minimum adjacent sorted gap")
    ax.set_title("Low-spectrum close-approach candidates")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(frameon=False)
    fig.tight_layout()
    path = args.output_dir / "out_of_plane_lambda_beta_min_gap_summary.png"
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return path


def write_report(
    args: Args,
    beta_values: np.ndarray,
    cases: dict[tuple[float, float, float], RootCase],
    gap_rows: list[dict[str, object]],
    output_paths: list[Path],
) -> Path:
    report_path = args.output_dir / "out_of_plane_lambda_beta_mu_eta_report.md"
    warning_cases = [case for case in cases.values() if case.warnings]
    missing_low = [case for case in cases.values() if len(case.roots) < args.n_roots_low]
    missing_extended = [case for case in cases.values() if len(case.roots) < args.n_roots_extended]
    torsion_root = torsion_roots_uniform_eta0_beta0(1, args.epsilon, args.poisson)
    extended_includes_torsion_region = args.lambda_max_extended >= torsion_root
    finite_gaps = [row for row in gap_rows if isfinite(float(row["g_min"]))]
    min_gap = min(finite_gaps, key=lambda row: float(row["g_min"])) if finite_gaps else None

    lines = [
        "# Out-of-Plane Lambda(beta) Sorted-Frequency Diagnostics",
        "",
        "These are diagnostic-only sorted-frequency plots for the out-of-plane",
        "Euler--Bernoulli bending plus Saint-Venant torsion determinant.",
        "",
        "They are not descendant-branch tracking and not crossing verification.",
        "No 3D FEM was run.",
        "",
        "## Parameters",
        "",
        f"- epsilon: {args.epsilon:g}",
        f"- poisson: {args.poisson:g}",
        f"- beta grid: {float(beta_values[0]):g}..{float(beta_values[-1]):g} deg, "
        f"{len(beta_values)} points, nominal step {args.beta_step:g} deg",
        f"- mu values: {', '.join(f'{value:g}' for value in args.mu_values)}",
        f"- eta values: {', '.join(f'{value:g}' for value in args.eta_values)}",
        f"- low spectrum roots requested: {args.n_roots_low}",
        f"- extended spectrum roots requested: {args.n_roots_extended}",
        "",
        "## Root Finder Settings",
        "",
        f"- lambda_max_extended: {args.lambda_max_extended:g}",
        f"- root scan step: {args.root_scan_step:g}",
        "- Lambda=0 is excluded by the out-of-plane root finder.",
        "",
        "## Root Availability",
        "",
        f"- parameter points solved: {len(cases)}",
        f"- root warning cases: {len(warning_cases)}",
        f"- cases missing low-spectrum roots: {len(missing_low)}",
        f"- cases missing extended-spectrum roots: {len(missing_extended)}",
        "",
        "## Observations",
        "",
        "- The first 8 roots are expected to be mostly bending-dominated for",
        "  epsilon=0.0025; the first beta=0, eta=0 torsion root is above the",
        "  low-spectrum window.",
        f"- first beta=0, eta=0 torsion root: {torsion_root:.12g}",
        f"- extended plots include the first torsion-root region: "
        f"{'yes' if extended_includes_torsion_region else 'no'}",
        "- Mode-character tags are TODO for this analytic root plotting script;",
        "  the separate 1D FEM validation CSV records bending/torsion energy",
        "  fractions.",
        "",
        "## Minimum Adjacent Gap Summary",
        "",
    ]
    if min_gap is None:
        lines.append("- no finite low-spectrum adjacent gaps were available.")
    else:
        lines.append(
            "- smallest low-spectrum adjacent sorted gap: "
            f"{float(min_gap['g_min']):.6e} at eta={float(min_gap['eta']):g}, "
            f"mu={float(min_gap['mu']):g}, beta={float(min_gap['beta_at_g_min_deg']):g} deg, "
            f"pair {min_gap['pair_at_g_min']}."
        )
        lines.append("- This is a close approach candidate only, not a crossing claim.")

    lines.extend(["", "## Generated Files", ""])
    for path in [*output_paths, report_path]:
        lines.append(f"- `{path.relative_to(REPO_ROOT)}`")
    lines.append("")

    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def run(args: Args) -> dict[str, object]:
    beta_values = _beta_grid(args.beta_min, args.beta_max, args.beta_step)
    cases = compute_roots(args, beta_values)

    output_paths: list[Path] = []
    root_rows = build_root_rows(args, cases)
    roots_csv = args.output_dir / "out_of_plane_lambda_beta_mu_eta_sorted_roots.csv"
    _csv_write(
        roots_csv,
        root_rows,
        [
            "beta_deg",
            "beta_rad",
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "sorted_index",
            "Lambda",
            "spectrum_window",
            "root_solver_warning",
            "notes",
        ],
    )
    output_paths.append(roots_csv)

    gap_rows = build_gap_summary(args, cases, beta_values)
    gap_csv = args.output_dir / "out_of_plane_lambda_beta_gap_summary.csv"
    _csv_write(
        gap_csv,
        gap_rows,
        [
            "mu",
            "eta",
            "epsilon",
            "poisson",
            "g_min",
            "beta_at_g_min_deg",
            "pair_at_g_min",
            "Lambda_lower",
            "Lambda_upper",
            "notes",
        ],
    )
    output_paths.append(gap_csv)

    for eta in args.eta_values:
        output_paths.append(plot_eta_figure(args, cases, beta_values, eta=eta, window="low"))
        output_paths.append(plot_eta_figure(args, cases, beta_values, eta=eta, window="extended"))
    output_paths.append(plot_overview(args, cases, beta_values))
    output_paths.append(plot_gap_summary(args, gap_rows))
    report_path = write_report(args, beta_values, cases, gap_rows, output_paths)
    output_paths.append(report_path)

    warning_cases = [case for case in cases.values() if case.warnings]
    missing_low = [case for case in cases.values() if len(case.roots) < args.n_roots_low]
    missing_extended = [case for case in cases.values() if len(case.roots) < args.n_roots_extended]
    finite_gaps = [row for row in gap_rows if isfinite(float(row["g_min"]))]
    min_gap = min(finite_gaps, key=lambda row: float(row["g_min"])) if finite_gaps else None

    print(f"saved sorted-root CSV: {roots_csv}")
    print(f"saved gap summary CSV: {gap_csv}")
    print(f"saved report: {report_path}")
    print(f"generated PNG count: {sum(1 for path in output_paths if path.suffix.lower() == '.png')}")
    print(f"beta grid: {float(beta_values[0]):g}..{float(beta_values[-1]):g} deg, {len(beta_values)} points")
    print(f"mu values: {', '.join(f'{value:g}' for value in args.mu_values)}")
    print(f"eta values: {', '.join(f'{value:g}' for value in args.eta_values)}")
    print(f"root warning cases: {len(warning_cases)}")
    print(f"missing low root cases: {len(missing_low)}")
    print(f"missing extended root cases: {len(missing_extended)}")
    if min_gap is not None:
        print(
            "min adjacent low-spectrum gap: "
            f"{float(min_gap['g_min']):.6e} at eta={float(min_gap['eta']):g}, "
            f"mu={float(min_gap['mu']):g}, beta={float(min_gap['beta_at_g_min_deg']):g} deg, "
            f"pair {min_gap['pair_at_g_min']}"
        )
    return {
        "output_paths": output_paths,
        "warning_cases": len(warning_cases),
        "missing_low_cases": len(missing_low),
        "missing_extended_cases": len(missing_extended),
        "min_gap": min_gap,
        "beta_count": len(beta_values),
    }


def main(argv: list[str] | None = None) -> dict[str, object]:
    return run(parse_args(argv))


if __name__ == "__main__":
    main()
