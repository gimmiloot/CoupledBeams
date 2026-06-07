from __future__ import annotations

import argparse
import csv
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
if str(SCRIPT_PATH.parent) not in sys.path:
    sys.path.insert(0, str(SCRIPT_PATH.parent))

import plot_lambda_mu_beta90_eta0p1_with_single_rod_refs as single_eta  # noqa: E402


DEFAULT_BETA_DEG = 90.0
DEFAULT_EPSILON = 0.0025
DEFAULT_ETA_VALUES = (0.1, 0.0, -0.1)
DEFAULT_MU_MIN = 0.0
DEFAULT_MU_MAX = 0.8
DEFAULT_MU_STEP = 0.01
DEFAULT_N_SYSTEM_ROOTS = 6
DEFAULT_N_REFERENCE_ROOTS = 6
DEFAULT_N_LONG_REFERENCE_ROOTS_PLOTTED = 6
DEFAULT_READABLE_YMIN = 0.0
DEFAULT_READABLE_YMAX = 13.0
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "lambda_mu_beta90_eta_scan_single_rod_refs"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "lambda_mu_beta90_eta_scan_single_rod_refs"

COMBINED_SUMMARY_NAME = "eta_comparison_summary.csv"
COMBINED_REPORT_NAME = "eta_scan_report.md"

COMBINED_SUMMARY_FIELDS = [
    "eta",
    "system_sorted_index",
    "dominant_reference_family",
    "fraction_of_mu_points",
    "number_of_family_switches",
    "notes",
]


@dataclass(frozen=True)
class ScanArgs:
    beta_deg: float
    epsilon: float
    eta_values: tuple[float, ...]
    mu_min: float
    mu_max: float
    mu_step: float
    n_system_roots: int
    n_reference_roots: int
    n_long_reference_roots_plotted: int
    readable_ymin: float
    readable_ymax: float
    auto_readable_ymax: bool
    output_dir: Path
    smoke: bool


@dataclass(frozen=True)
class EtaCaseResult:
    eta: float
    eta_token: str
    output_dir: Path
    mu_values: np.ndarray
    system_grid: np.ndarray
    system_plot_grid: np.ndarray
    warnings: Sequence[Sequence[str]]
    spike_rows: Sequence[dict[str, object]]
    branch_summary: Sequence[dict[str, object]]
    main_plot_stats: dict[str, object]
    full_plot_stats: dict[str, object]
    output_paths: Sequence[Path]


def _fmt(value: object) -> str:
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not isfinite(value_f):
        return "nan"
    return f"{value_f:.16g}"


def eta_token(value: float) -> str:
    value_f = float(value)
    if abs(value_f) < 5.0e-15:
        return "eta_0p0"
    prefix = "eta_p" if value_f > 0.0 else "eta_m"
    return prefix + f"{abs(value_f):.10g}".replace(".", "p")


def file_prefix(eta: float) -> str:
    return f"lambda_mu_beta90_eta_{eta_token(eta).removeprefix('eta_')}"


def output_dir(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def write_csv(path: Path, rows: Sequence[dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def validate_scan_args(args: ScanArgs) -> None:
    if not isfinite(args.beta_deg):
        raise ValueError("beta-deg must be finite.")
    if not isfinite(args.epsilon) or args.epsilon <= 0.0:
        raise ValueError("epsilon must be positive and finite.")
    if not args.eta_values:
        raise ValueError("eta-values must not be empty.")
    if any(not (-1.0 < eta < 1.0) for eta in args.eta_values):
        raise ValueError("all eta values must lie inside (-1, 1).")
    if args.mu_step <= 0.0:
        raise ValueError("mu-step must be positive.")
    if args.mu_max < args.mu_min:
        raise ValueError("mu-max must be greater than or equal to mu-min.")
    if args.n_system_roots <= 0:
        raise ValueError("n-system-roots must be positive.")
    if args.n_reference_roots <= 0:
        raise ValueError("n-reference-roots must be positive.")
    if args.n_long_reference_roots_plotted <= 0:
        raise ValueError("n-long-reference-roots-plotted must be positive.")
    if args.n_long_reference_roots_plotted > args.n_reference_roots:
        raise ValueError("n-long-reference-roots-plotted must not exceed n-reference-roots.")
    if not isfinite(args.readable_ymin) or not isfinite(args.readable_ymax):
        raise ValueError("readable y-limits must be finite.")
    if args.readable_ymax <= args.readable_ymin:
        raise ValueError("readable-ymax must be greater than readable-ymin.")


def single_args(scan_args: ScanArgs, eta: float, eta_dir: Path) -> single_eta.Args:
    return single_eta.Args(
        beta_deg=scan_args.beta_deg,
        epsilon=scan_args.epsilon,
        eta=float(eta),
        mu_min=scan_args.mu_min,
        mu_max=scan_args.mu_max,
        mu_step=scan_args.mu_step,
        n_system_roots=scan_args.n_system_roots,
        n_reference_roots=scan_args.n_reference_roots,
        n_reference_roots_plotted=scan_args.n_long_reference_roots_plotted,
        output_dir=eta_dir,
        smoke=scan_args.smoke,
    )


def family_style(rod_label: str, boundary_condition: str) -> dict[str, object]:
    return single_eta.family_style(rod_label, boundary_condition)


def line_label(rod_label: str, boundary_condition: str) -> str:
    bc = "CP" if boundary_condition == "clamped_pinned" else "CC"
    return f"{rod_label} {bc} refs"


def y_visible_limit(
    system_plot_grid: np.ndarray,
    reference_arrays: dict[tuple[str, str], np.ndarray],
    n_long_reference_roots_plotted: int,
) -> float:
    values: list[np.ndarray] = []
    system_finite = system_plot_grid[np.isfinite(system_plot_grid)]
    if system_finite.size:
        values.append(system_finite)
    for boundary_condition in ("clamped_pinned", "clamped_clamped"):
        arr = reference_arrays[("long", boundary_condition)][: int(n_long_reference_roots_plotted)]
        finite = arr[np.isfinite(arr)]
        if finite.size:
            values.append(finite)
    if not values:
        return 1.0
    return 1.05 * float(np.max(np.concatenate(values)))


def plot_eta_case(
    path: Path,
    *,
    scan_args: ScanArgs,
    eta: float,
    mu_values: np.ndarray,
    system_plot_grid: np.ndarray,
    reference_arrays: dict[tuple[str, str], np.ndarray],
    full_range: bool,
) -> dict[str, object]:
    fig, ax = plt.subplots(figsize=(10.8, 6.35))
    system_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for mode_index in range(scan_args.n_system_roots):
        ax.plot(
            mu_values,
            system_plot_grid[:, mode_index],
            color=system_colors[mode_index % len(system_colors)],
            linewidth=2.15,
            linestyle="-",
            zorder=5,
        )

    auto_y_limit = y_visible_limit(
        system_plot_grid,
        reference_arrays,
        scan_args.n_long_reference_roots_plotted,
    )
    y_limit = auto_y_limit if (full_range or scan_args.auto_readable_ymax) else float(scan_args.readable_ymax)
    y_min = 0.0 if full_range else float(scan_args.readable_ymin)
    short_curve_plotted = 0
    short_curve_omitted = 0
    short_curve_clipped = 0
    long_curve_plotted = 0
    all_plotted_values: list[np.ndarray] = []
    for (rod_label, boundary_condition), values in reference_arrays.items():
        style = family_style(rod_label, boundary_condition)
        if rod_label == "long":
            max_modes = min(scan_args.n_long_reference_roots_plotted, values.shape[0])
        else:
            max_modes = values.shape[0]
        for ref_index in range(max_modes):
            raw = np.asarray(values[ref_index], dtype=float)
            finite = raw[np.isfinite(raw)]
            if not finite.size:
                continue
            if full_range:
                plot_values = raw
                should_plot = True
            elif rod_label == "short":
                should_plot = bool(np.nanmin(finite) <= y_limit and np.nanmax(finite) >= y_min)
                plot_values = np.where((raw >= y_min) & (raw <= y_limit), raw, np.nan)
            else:
                should_plot = True
                plot_values = raw
            if not should_plot:
                short_curve_omitted += 1
                continue
            if rod_label == "short":
                short_curve_plotted += 1
                if np.any(np.isfinite(raw) & ((raw > y_limit) | (raw < y_min))):
                    short_curve_clipped += 1
            else:
                long_curve_plotted += 1
            all_plotted_values.append(plot_values[np.isfinite(plot_values)])
            ax.plot(
                mu_values,
                plot_values,
                color=style["color"],
                linewidth=1.22,
                linestyle=style["linestyle"],
                alpha=float(style["alpha"]),
                zorder=2,
            )

    if full_range:
        finite_arrays = [system_plot_grid[np.isfinite(system_plot_grid)]]
        finite_arrays.extend(arr[np.isfinite(arr)] for arr in reference_arrays.values())
        finite_arrays = [arr for arr in finite_arrays if arr.size]
        ymax = 1.06 * float(np.max(np.concatenate(finite_arrays))) if finite_arrays else 1.0
    else:
        ymax = y_limit
    ax.set_ylim(y_min, ymax)
    ax.set_xlabel("mu")
    ax.set_ylabel("Lambda")
    if full_range:
        scope = "full reference range"
    elif scan_args.auto_readable_ymax:
        scope = f"auto readable range: {_fmt(y_min)} <= Lambda <= {_fmt(ymax)}"
    else:
        scope = f"readable range: {_fmt(y_min)} <= Lambda <= {_fmt(ymax)}"
    ax.set_title(
        f"Lambda(mu), beta = {scan_args.beta_deg:g} deg, epsilon = {scan_args.epsilon:g}, eta = {eta:g}\n"
        f"solid: system in-plane sorted frequencies; dashed: single-rod references ({scope}; full-range plot saved separately)"
    )
    ax.grid(True, color="0.88", linewidth=0.65)
    handles = [
        Line2D([0], [0], color="black", lw=2.15, ls="-", label=f"system sorted roots, first {scan_args.n_system_roots}"),
    ]
    for rod_label, boundary_condition in (
        ("long", "clamped_pinned"),
        ("long", "clamped_clamped"),
        ("short", "clamped_pinned"),
        ("short", "clamped_clamped"),
    ):
        style = family_style(rod_label, boundary_condition)
        handles.append(
            Line2D(
                [0],
                [0],
                color=style["color"],
                lw=1.4,
                ls=style["linestyle"],
                alpha=style["alpha"],
                label=line_label(rod_label, boundary_condition),
            )
        )
    ax.legend(handles=handles, loc="upper left", fontsize=8.5, frameon=False, handlelength=3.6)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return {
        "y_visible_max": ymax,
        "y_visible_min": y_min,
        "auto_y_visible_max": auto_y_limit,
        "readable_ymax_mode": "auto" if (not full_range and scan_args.auto_readable_ymax) else ("full" if full_range else "fixed"),
        "long_reference_curves_plotted": long_curve_plotted,
        "short_reference_curves_plotted": short_curve_plotted,
        "short_reference_curves_omitted": short_curve_omitted,
        "short_reference_curves_clipped": short_curve_clipped,
        "full_range": "yes" if full_range else "no",
    }


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_eta_report(
    path: Path,
    *,
    scan_args: ScanArgs,
    eta: float,
    result: EtaCaseResult,
) -> None:
    finite_system = result.system_grid[np.isfinite(result.system_grid)]
    missing = int(np.count_nonzero(~np.isfinite(result.system_grid)))
    warnings = sum(1 for row in result.warnings for item in row if item != "none")
    stats = result.main_plot_stats
    lines = [
        f"# Lambda(mu) Single-Rod Reference Diagnostic, eta = {_fmt(eta)}",
        "",
        "This is a diagnostic-only analytic plot and nearest-frequency comparison.",
        "",
        "## Parameters",
        "",
        f"- beta_deg: {_fmt(scan_args.beta_deg)}",
        f"- epsilon: {_fmt(scan_args.epsilon)}",
        f"- eta: {_fmt(eta)}",
        f"- mu grid: {_fmt(float(result.mu_values[0]))} to {_fmt(float(result.mu_values[-1]))}, {len(result.mu_values)} points",
        f"- nominal mu step: {_fmt(0.2 if scan_args.smoke else scan_args.mu_step)}",
        f"- system sorted roots per mu: {scan_args.n_system_roots}",
        f"- reference roots saved for all four families: {scan_args.n_reference_roots}",
        f"- long-rod reference roots plotted in readable plot per CP/CC family: {scan_args.n_long_reference_roots_plotted}",
        f"- readable-ymin: {_fmt(scan_args.readable_ymin)}",
        f"- readable-ymax: {_fmt(scan_args.readable_ymax)}",
        f"- auto-readable-ymax: {'yes' if scan_args.auto_readable_ymax else 'no'}",
        "",
        "## Readable Plot Window",
        "",
        "- The readable PNG uses the configured y-range only; by default this is 0 <= Lambda <= 13.",
        "- Full-range PNGs keep the complete plotted reference range.",
        "- CSV data, nearest matches, and candidate tables are not clipped.",
        "- High short-rod reference curves are intentionally clipped or omitted from the readable view only.",
        "- This is a visualization choice, not a change in frequencies or matching.",
        f"- automatic system-plus-long-reference y_visible_max would be: {_fmt(stats['auto_y_visible_max'])}",
        f"- readable y-limit mode: {stats['readable_ymax_mode']}",
        f"- readable y_visible_min: {_fmt(stats['y_visible_min'])}",
        f"- readable y_visible_max: {_fmt(stats['y_visible_max'])}",
        f"- long-reference curves plotted in readable PNG: {stats['long_reference_curves_plotted']}",
        f"- short-reference curves plotted in readable PNG: {stats['short_reference_curves_plotted']}",
        f"- short-reference curves omitted from readable PNG: {stats['short_reference_curves_omitted']}",
        f"- short-reference curves clipped by readable y-window: {stats['short_reference_curves_clipped']}",
        "",
        "## Root And Spike Checks",
        "",
        f"- finite system roots: {int(finite_system.size)} / {result.system_grid.size}",
        f"- missing system roots: {missing}",
        f"- system root warning entries: {warnings}",
        f"- suspicious spike audit rows: {len(result.spike_rows)}",
        "",
        "## Dominant Nearest Reference Families",
        "",
    ]
    lines.extend(
        markdown_table(
            ["sorted index", "dominant family", "fraction", "switches", "switch locations"],
            [
                [
                    row["system_sorted_index"],
                    row["dominant_reference_family"],
                    row["fraction_of_mu_points"],
                    row["number_of_family_switches"],
                    row["switch_mu_locations"] or "none",
                ]
                for row in result.branch_summary
            ],
        )
    )
    lines.extend(
        [
            "",
            "## Interpretation Scope",
            "",
            "- This is frequency-proximity diagnostic bookkeeping only.",
            "- Family switches are nearest-reference changes for sorted frequencies, not descendant branch switches.",
            "- No crossing/no-crossing or strict positive-gap claim is made.",
            "",
            "## Output Files",
            "",
        ]
    )
    for output_path in result.output_paths:
        lines.append(f"- {output_path.name}")
    lines.extend(
        [
            "",
            "## Protected Scope",
            "",
            "No descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article workspace edits, old determinant edits, old solver edits, baseline-result edits, or analytic formula changes are part of this workflow.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_combined_report(path: Path, scan_args: ScanArgs, results: Sequence[EtaCaseResult]) -> None:
    lines = [
        "# Beta 90 Eta Scan Single-Rod Reference Diagnostic",
        "",
        "This report summarizes diagnostic-only sorted in-plane Lambda(mu) plots against isolated-rod references.",
        "",
        "## Parameters",
        "",
        f"- beta_deg: {_fmt(scan_args.beta_deg)}",
        f"- epsilon: {_fmt(scan_args.epsilon)}",
        f"- eta values: {', '.join(_fmt(result.eta) for result in results)}",
        f"- mu grid: {_fmt(float(results[0].mu_values[0]))} to {_fmt(float(results[0].mu_values[-1]))}, {len(results[0].mu_values)} points",
        f"- system roots per mu: {scan_args.n_system_roots}",
        f"- reference roots saved per family: {scan_args.n_reference_roots}",
        f"- long-reference roots plotted per long CP/CC family: {scan_args.n_long_reference_roots_plotted}",
        f"- readable y-range: {_fmt(scan_args.readable_ymin)} <= Lambda <= {_fmt(scan_args.readable_ymax)}",
        f"- auto-readable-ymax: {'yes' if scan_args.auto_readable_ymax else 'no'}",
        "- full-range plots keep the complete plotted reference range",
        "- CSV data and matching tables are not clipped by the readable y-range",
        "",
        "## Combined Dominant Families",
        "",
    ]
    table_rows: list[list[object]] = []
    for result in results:
        for row in result.branch_summary:
            table_rows.append(
                [
                    _fmt(result.eta),
                    row["system_sorted_index"],
                    row["dominant_reference_family"],
                    row["fraction_of_mu_points"],
                    row["number_of_family_switches"],
                ]
            )
    lines.extend(
        markdown_table(
            ["eta", "sorted index", "dominant family", "fraction", "switches"],
            table_rows,
        )
    )
    lines.extend(
        [
            "",
            "## Scope",
            "",
            "System curves are solid sorted in-plane Euler-Bernoulli frequencies. Reference curves are dashed isolated-rod CP/FP and CC/FF bending frequencies. The readable PNGs use the configured y-range only; high short-rod reference curves may be clipped or omitted from that view, while CSV data and full-range PNGs remain complete. This eta scan uses no descendant tracking, FEM, 3D FEM, Gmsh, CalculiX, article artifacts, determinant edits, old-solver edits, baseline-result edits, or strict gap/crossing verification.",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def output_paths_for(case_dir: Path, eta: float) -> dict[str, Path]:
    prefix = file_prefix(eta)
    return {
        "system": case_dir / f"{prefix}_system_in_plane_sorted.csv",
        "references": case_dir / f"{prefix}_single_rod_references.csv",
        "matches": case_dir / f"{prefix}_nearest_reference_matches.csv",
        "all_candidates": case_dir / f"{prefix}_all_candidate_matches.csv",
        "branch_summary": case_dir / f"{prefix}_branch_family_summary.csv",
        "spike_audit": case_dir / f"{prefix}_spike_audit.csv",
        "main_png": case_dir / f"{prefix}_in_plane_vs_single_rod_refs.png",
        "full_png": case_dir / f"{prefix}_in_plane_vs_single_rod_refs_full.png",
        "report": case_dir / f"{prefix}_report.md",
    }


def run_eta_case(scan_args: ScanArgs, eta: float) -> EtaCaseResult:
    token = eta_token(eta)
    case_dir = scan_args.output_dir / token
    args = single_args(scan_args, eta, case_dir)
    mu_values = single_eta.mu_grid(args)
    single_eta.validate_args(args, mu_values)
    references, references_by_mu = single_eta.build_reference_values(args, mu_values)
    system_grid, warnings, notes = single_eta.solve_system_roots(args, mu_values)
    system_plot_grid, spike_rows = single_eta.audit_and_repair_spikes(args, mu_values, system_grid, warnings, notes)
    matches, all_candidates = single_eta.compare_candidates(args, mu_values, system_grid, references_by_mu)
    branch_summary = single_eta.branch_family_summary(args, mu_values, matches)
    reference_arrays = single_eta.reference_curve_arrays(args, mu_values, references_by_mu)
    paths = output_paths_for(case_dir, eta)

    write_csv(paths["system"], single_eta.system_rows(args, mu_values, system_grid, warnings, notes), single_eta.SYSTEM_FIELDS)
    write_csv(paths["references"], single_eta.reference_rows(args, references), single_eta.REFERENCE_FIELDS)
    write_csv(paths["matches"], matches, single_eta.MATCH_FIELDS)
    write_csv(paths["all_candidates"], all_candidates, single_eta.ALL_CANDIDATE_FIELDS)
    write_csv(paths["branch_summary"], branch_summary, single_eta.BRANCH_SUMMARY_FIELDS)
    main_stats = plot_eta_case(
        paths["main_png"],
        scan_args=scan_args,
        eta=eta,
        mu_values=mu_values,
        system_plot_grid=system_plot_grid,
        reference_arrays=reference_arrays,
        full_range=False,
    )
    full_stats = plot_eta_case(
        paths["full_png"],
        scan_args=scan_args,
        eta=eta,
        mu_values=mu_values,
        system_plot_grid=system_plot_grid,
        reference_arrays=reference_arrays,
        full_range=True,
    )

    output_paths = [
        paths["system"],
        paths["references"],
        paths["matches"],
        paths["all_candidates"],
        paths["branch_summary"],
        paths["main_png"],
        paths["full_png"],
    ]
    if spike_rows:
        write_csv(paths["spike_audit"], spike_rows, single_eta.SPIKE_AUDIT_FIELDS)
        output_paths.append(paths["spike_audit"])
    pending_result = EtaCaseResult(
        eta=eta,
        eta_token=token,
        output_dir=case_dir,
        mu_values=mu_values,
        system_grid=system_grid,
        system_plot_grid=system_plot_grid,
        warnings=warnings,
        spike_rows=spike_rows,
        branch_summary=branch_summary,
        main_plot_stats=main_stats,
        full_plot_stats=full_stats,
        output_paths=[*output_paths, paths["report"]],
    )
    write_eta_report(paths["report"], scan_args=scan_args, eta=eta, result=pending_result)
    return pending_result


def combined_summary_rows(results: Sequence[EtaCaseResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        for row in result.branch_summary:
            rows.append(
                {
                    "eta": _fmt(result.eta),
                    "system_sorted_index": row["system_sorted_index"],
                    "dominant_reference_family": row["dominant_reference_family"],
                    "fraction_of_mu_points": row["fraction_of_mu_points"],
                    "number_of_family_switches": row["number_of_family_switches"],
                    "notes": row["notes"],
                }
            )
    return rows


def parse_args(argv: Sequence[str] | None = None) -> ScanArgs:
    parser = argparse.ArgumentParser(
        description="Eta scan for beta=90 sorted in-plane Lambda(mu) plots with single-rod references."
    )
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--eta-values", type=float, nargs="+", default=list(DEFAULT_ETA_VALUES))
    parser.add_argument("--mu-min", type=float, default=DEFAULT_MU_MIN)
    parser.add_argument("--mu-max", type=float, default=DEFAULT_MU_MAX)
    parser.add_argument("--mu-step", type=float, default=DEFAULT_MU_STEP)
    parser.add_argument("--n-system-roots", type=int, default=DEFAULT_N_SYSTEM_ROOTS)
    parser.add_argument("--n-reference-roots", type=int, default=DEFAULT_N_REFERENCE_ROOTS)
    parser.add_argument("--n-long-reference-roots-plotted", type=int, default=DEFAULT_N_LONG_REFERENCE_ROOTS_PLOTTED)
    parser.add_argument("--readable-ymin", type=float, default=DEFAULT_READABLE_YMIN)
    parser.add_argument("--readable-ymax", type=float, default=DEFAULT_READABLE_YMAX)
    parser.add_argument(
        "--auto-readable-ymax",
        action="store_true",
        help="Use the previous system-plus-long-reference adaptive y-limit for readable PNGs.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(argv)
    eta_values = tuple(float(value) for value in ns.eta_values)
    n_system_roots = int(ns.n_system_roots)
    n_reference_roots = int(ns.n_reference_roots)
    n_long_reference_roots_plotted = int(ns.n_long_reference_roots_plotted)
    out_dir = output_dir(Path(ns.output_dir))
    if bool(ns.smoke):
        eta_values = DEFAULT_ETA_VALUES
        n_system_roots = 4
        n_reference_roots = 4
        n_long_reference_roots_plotted = 4
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            out_dir = SMOKE_OUTPUT_DIR
    args = ScanArgs(
        beta_deg=float(ns.beta_deg),
        epsilon=float(ns.epsilon),
        eta_values=eta_values,
        mu_min=float(ns.mu_min),
        mu_max=float(ns.mu_max),
        mu_step=float(ns.mu_step),
        n_system_roots=n_system_roots,
        n_reference_roots=n_reference_roots,
        n_long_reference_roots_plotted=n_long_reference_roots_plotted,
        readable_ymin=float(ns.readable_ymin),
        readable_ymax=float(ns.readable_ymax),
        auto_readable_ymax=bool(ns.auto_readable_ymax),
        output_dir=out_dir,
        smoke=bool(ns.smoke),
    )
    validate_scan_args(args)
    return args


def run(scan_args: ScanArgs) -> dict[str, object]:
    results = [run_eta_case(scan_args, eta) for eta in scan_args.eta_values]
    combined_csv = scan_args.output_dir / COMBINED_SUMMARY_NAME
    combined_report = scan_args.output_dir / COMBINED_REPORT_NAME
    write_csv(combined_csv, combined_summary_rows(results), COMBINED_SUMMARY_FIELDS)
    write_combined_report(combined_report, scan_args, results)
    return {
        "results": results,
        "combined_csv": combined_csv,
        "combined_report": combined_report,
    }


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    scan_args = parse_args(argv)
    result = run(scan_args)
    results: Sequence[EtaCaseResult] = result["results"]
    print("diagnostic-only eta scan Lambda(mu) map with single-rod references")
    print(f"beta_deg={scan_args.beta_deg:g}, epsilon={scan_args.epsilon:g}")
    print(f"eta values: {', '.join(_fmt(case.eta) for case in results)}")
    for case in results:
        missing = int(np.count_nonzero(~np.isfinite(case.system_grid)))
        warnings = sum(1 for row in case.warnings for item in row if item != "none")
        print(
            f"eta={case.eta:g}: mu grid {len(case.mu_values)} points "
            f"{float(case.mu_values[0]):g}..{float(case.mu_values[-1]):g}; "
            f"missing roots={missing}; root warnings={warnings}; spike rows={len(case.spike_rows)}"
        )
        print(
            "  readable y-range: "
            f"mode={case.main_plot_stats['readable_ymax_mode']}; "
            f"y_visible_min={_fmt(case.main_plot_stats['y_visible_min'])}; "
            f"y_visible_max={_fmt(case.main_plot_stats['y_visible_max'])}; "
            f"auto_y_visible_max={_fmt(case.main_plot_stats['auto_y_visible_max'])}; "
            f"long curves plotted={case.main_plot_stats['long_reference_curves_plotted']}; "
            f"short curves plotted={case.main_plot_stats['short_reference_curves_plotted']}; "
            f"short curves omitted={case.main_plot_stats['short_reference_curves_omitted']}; "
            f"short curves clipped={case.main_plot_stats['short_reference_curves_clipped']}"
        )
        for row in case.branch_summary:
            print(
                f"  sorted {row['system_sorted_index']}: {row['dominant_reference_family']} "
                f"(fraction={row['fraction_of_mu_points']}, switches={row['number_of_family_switches']})"
            )
        for output_path in case.output_paths:
            print(f"  saved: {output_path}")
    print(f"combined summary: {result['combined_csv']}")
    print(f"combined report: {result['combined_report']}")
    print("solid curves: sorted in-plane system roots; dashed curves: single-rod references")
    print("no descendant tracking; no FEM/Gmsh/CalculiX; no article/formula/baseline changes")
    return result


if __name__ == "__main__":
    main()
