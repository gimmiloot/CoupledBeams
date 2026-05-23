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


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib.analytic_branch_tracking import (  # noqa: E402
    branch_id_from_base_sorted_index,
    track_path,
)


DEFAULT_EPSILON = 0.0025
DEFAULT_MUS = [0.0, 0.2, 0.5, 0.8, 0.9]
DEFAULT_BETA_MIN = 0.0
DEFAULT_BETA_MAX = 90.0
DEFAULT_BETA_STEP = 1.0
DEFAULT_NUM_ROOTS = 8
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
BRANCH_FOCUS_IDS = ["bending_desc_01", "bending_desc_02", "bending_desc_04", "bending_desc_05"]


def number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace(".", "p")


def beta_grid(beta_min: float, beta_max: float, beta_step: float) -> np.ndarray:
    if beta_step <= 0.0:
        raise ValueError("--beta-step must be positive.")
    count = int(np.floor((float(beta_max) - float(beta_min)) / float(beta_step) + 1e-12)) + 1
    values = float(beta_min) + float(beta_step) * np.arange(count, dtype=float)
    if len(values) == 0 or values[-1] < float(beta_max) - 1e-12:
        values = np.append(values, float(beta_max))
    return np.round(values, 12)


def mu_ramp_values(mu: float) -> np.ndarray:
    if abs(float(mu)) <= 1e-15:
        return np.array([0.0], dtype=float)
    steps = max(2, int(round(abs(float(mu)) / 0.01)) + 1)
    return np.linspace(0.0, float(mu), steps)


def tracking_path_for_fixed_mu(mu: float, betas: Sequence[float]) -> list[tuple[float, float, str, int]]:
    path: list[tuple[float, float, str, int]] = [(0.0, 0.0, "base", 0)]
    for step_index, mu_value in enumerate(mu_ramp_values(mu)[1:], start=1):
        path.append((0.0, float(mu_value), "mu_preset", step_index))
    beta_values = [float(beta) for beta in betas]
    for step_index, beta in enumerate(beta_values, start=1):
        if abs(beta) <= 1e-15:
            continue
        path.append((float(beta), float(mu), "beta_sweep", step_index))
    return path


def points_at(result, *, mu: float, beta: float, atol: float = 1e-10):
    return [
        point
        for point in result.points
        if abs(float(point.mu) - float(mu)) <= atol and abs(float(point.beta) - float(beta)) <= atol
    ]


def min_mac_to_beta(result, *, branch_id: str, mu: float, beta: float) -> float:
    macs = [
        float(point.mac_to_previous)
        for point in result.points_for_branch(branch_id)
        if abs(float(point.mu) - float(mu)) <= 1e-10
        and float(point.beta) <= float(beta) + 1e-10
        and point.step_type != "base"
        and np.isfinite(float(point.mac_to_previous))
    ]
    return float(min(macs)) if macs else float("nan")


def branch_ids_for_plot(num_roots: int) -> list[str]:
    return [branch_id_from_base_sorted_index(index) for index in range(1, int(num_roots) + 1)]


def run_tracking_for_mu(*, epsilon: float, mu: float, betas: Sequence[float], num_roots: int):
    n_track = max(int(num_roots), max(int(branch.rsplit("_", 1)[1]) for branch in BRANCH_FOCUS_IDS))
    n_solve = max(n_track, int(num_roots), 12)
    kwargs = {
        "epsilon": float(epsilon),
        "path": tracking_path_for_fixed_mu(float(mu), betas),
        "n_track": n_track,
        "n_solve": n_solve,
        "shape_metric": "full",
        "required_branch_ids": branch_ids_for_plot(num_roots),
        "Lmin": 0.2,
        "Lmax0": 55.0,
        "scan_step": 0.02,
        "grow_factor": 1.35,
        "max_tries": 8,
    }
    try:
        return track_path(**kwargs, allow_low_mac=False), "strict", []
    except Exception as exc:
        warnings = [f"Strict tracking failed for mu={float(mu):g}: {exc}"]
        result = track_path(**kwargs, allow_low_mac=True)
        warnings.extend(str(item) for item in result.warnings)
        return result, "allow_low_mac", warnings


def write_csv(path: Path, rows: list[dict[str, float | int | str]], fieldnames: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fieldnames), extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def collect_rows_for_mu(*, result, epsilon: float, mu: float, betas: Sequence[float], num_roots: int):
    sorted_rows: list[dict[str, float | int | str]] = []
    branch_rows: list[dict[str, float | int | str]] = []
    focus_set = set(BRANCH_FOCUS_IDS)

    for beta in betas:
        case_points = points_at(result, mu=float(mu), beta=float(beta))
        if not case_points:
            continue
        sorted_lambdas = tuple(float(value) for value in case_points[0].sorted_lambdas)
        point_by_sorted = {int(point.current_sorted_index): point for point in case_points}

        for sorted_index in range(1, int(num_roots) + 1):
            Lambda = sorted_lambdas[sorted_index - 1] if sorted_index - 1 < len(sorted_lambdas) else float("nan")
            point = point_by_sorted.get(sorted_index)
            branch_id = str(point.branch_id) if point is not None else ""
            warning_flag = str(point.warning_flag) if point is not None else ""
            min_mac = min_mac_to_beta(result, branch_id=branch_id, mu=float(mu), beta=float(beta)) if branch_id else float("nan")
            sorted_rows.append(
                {
                    "epsilon": float(epsilon),
                    "mu": float(mu),
                    "beta": float(beta),
                    "sorted_index": int(sorted_index),
                    "branch_id": branch_id,
                    "Lambda": float(Lambda),
                    "warning_flag": warning_flag,
                    "min_MAC": float(min_mac),
                }
            )

        for point in case_points:
            if str(point.branch_id) not in focus_set:
                continue
            branch_rows.append(
                {
                    "epsilon": float(epsilon),
                    "mu": float(mu),
                    "beta": float(beta),
                    "branch_id": str(point.branch_id),
                    "current_sorted_index": int(point.current_sorted_index),
                    "Lambda": float(point.Lambda),
                    "warning_flag": str(point.warning_flag),
                    "min_MAC": float(min_mac_to_beta(result, branch_id=str(point.branch_id), mu=float(mu), beta=float(beta))),
                }
            )
    return sorted_rows, branch_rows


def plot_mu(*, result, output_dir: Path, epsilon: float, mu: float, betas: Sequence[float], num_roots: int) -> Path:
    output = output_dir / (
        f"lambda_beta_eps{number_token(float(epsilon))}_mu{number_token(float(mu))}_modes{int(num_roots)}.png"
    )
    fig, ax = plt.subplots(figsize=(9.5, 5.4))
    beta_array = np.asarray(betas, dtype=float)

    for branch_id in branch_ids_for_plot(num_roots):
        lambdas = []
        for beta in beta_array:
            try:
                point = result.point_at(branch_id, beta=float(beta), mu=float(mu))
                lambdas.append(float(point.Lambda))
            except KeyError:
                lambdas.append(float("nan"))
        ax.plot(beta_array, lambdas, marker="o", linewidth=1.5, markersize=2.8, label=branch_id)

    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("Lambda")
    ax.set_title(f"Lambda(beta) at fixed mu={float(mu):g}, epsilon={float(epsilon):g}")
    ax.grid(True, alpha=0.3)
    ax.legend(ncols=2, fontsize=8)
    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return output


def adjacent_gap_metrics(sorted_rows: Sequence[dict[str, float | int | str]]) -> dict[float, dict[str, float | int]]:
    by_mu_beta: dict[tuple[float, float], list[dict[str, float | int | str]]] = {}
    for row in sorted_rows:
        by_mu_beta.setdefault((float(row["mu"]), float(row["beta"])), []).append(row)

    out: dict[float, dict[str, float | int]] = {}
    for (mu, beta), rows in by_mu_beta.items():
        values = sorted((int(row["sorted_index"]), float(row["Lambda"])) for row in rows)
        lambdas = np.asarray([value for _, value in values], dtype=float)
        if len(lambdas) < 2 or not np.all(np.isfinite(lambdas)):
            continue
        gaps = np.diff(lambdas)
        idx = int(np.argmin(gaps))
        current = out.get(mu)
        if current is None or float(gaps[idx]) < float(current["min_gap"]):
            out[mu] = {
                "min_gap": float(gaps[idx]),
                "beta_at_min_gap": float(beta),
                "lower_sorted_index": int(values[idx][0]),
                "upper_sorted_index": int(values[idx + 1][0]),
            }
    return out


def branch_exchange_metrics(branch_rows: Sequence[dict[str, float | int | str]]) -> dict[float, dict[str, float | int]]:
    by_mu_branch: dict[tuple[float, str], list[dict[str, float | int | str]]] = {}
    for row in branch_rows:
        by_mu_branch.setdefault((float(row["mu"]), str(row["branch_id"])), []).append(row)

    out: dict[float, dict[str, float | int]] = {}
    for (mu, _branch_id), rows in by_mu_branch.items():
        indices = [int(row["current_sorted_index"]) for row in sorted(rows, key=lambda item: float(item["beta"]))]
        unique_count = len(set(indices))
        max_shift = max(indices) - min(indices) if indices else 0
        metrics = out.setdefault(mu, {"branches_with_index_change": 0, "max_sorted_index_shift": 0})
        if unique_count > 1:
            metrics["branches_with_index_change"] = int(metrics["branches_with_index_change"]) + 1
        metrics["max_sorted_index_shift"] = max(int(metrics["max_sorted_index_shift"]), int(max_shift))
    return out


def report_lines(
    *,
    args: argparse.Namespace,
    sorted_rows: list[dict[str, float | int | str]],
    branch_rows: list[dict[str, float | int | str]],
    plot_paths: Sequence[Path],
    tracking_status: dict[float, str],
    warnings: Sequence[str],
) -> list[str]:
    gap_metrics = adjacent_gap_metrics(sorted_rows)
    exchange_metrics = branch_exchange_metrics(branch_rows)
    lines = [
        "# Lambda(beta) fixed-mu diagnostic",
        "",
        "Diagnostic only. No determinant, formulas.py, solvers.py, python_fem.py, FEM physical model, or article figure changes.",
        "",
        "## Configuration",
        "",
        f"- epsilon: {float(args.epsilon):g}",
        f"- mus: {', '.join(f'{float(mu):g}' for mu in args.mus)}",
        f"- beta range: {float(args.beta_min):g} .. {float(args.beta_max):g} deg, step {float(args.beta_step):g}",
        f"- num roots / tracked branches: {int(args.num_roots)}",
        "",
        "## Outputs",
        "",
        f"- sorted CSV: `{args.sorted_csv.relative_to(REPO_ROOT)}`",
        f"- branch-focused CSV: `{args.branch_csv.relative_to(REPO_ROOT)}`",
    ]
    lines.extend(f"- plot: `{path.relative_to(REPO_ROOT)}`" for path in plot_paths)
    lines.extend(["", "## Empirical Summary", ""])
    lines.append("| mu | tracking | min adjacent gap | beta at min gap | pair | focused branches with index change | max sorted-index shift |")
    lines.append("| --- | --- | ---: | ---: | --- | ---: | ---: |")
    for mu in args.mus:
        gap = gap_metrics.get(float(mu), {})
        exchange = exchange_metrics.get(float(mu), {})
        pair = ""
        if gap:
            pair = f"{int(gap['lower_sorted_index'])}-{int(gap['upper_sorted_index'])}"
        lines.append(
            f"| {float(mu):g} | {tracking_status.get(float(mu), 'unknown')} | "
            f"{float(gap.get('min_gap', np.nan)):.6g} | {float(gap.get('beta_at_min_gap', np.nan)):.6g} | "
            f"{pair} | {int(exchange.get('branches_with_index_change', 0))} | "
            f"{int(exchange.get('max_sorted_index_shift', 0))} |"
        )

    lines.extend(["", "## Interpretation", ""])
    mu0_exchange = exchange_metrics.get(0.0, {"branches_with_index_change": 0, "max_sorted_index_shift": 0})
    if int(mu0_exchange.get("branches_with_index_change", 0)) == 0:
        lines.append("- At mu=0, the tracked branches keep their sorted order in this beta sweep, so there is little/no branch-exchange evidence in Lambda(beta).")
    else:
        gap0 = gap_metrics.get(0.0, {})
        if gap0:
            lines.append(
                "- At mu=0, the tracked branches do change sorted index; the closest event occurs near "
                f"beta={float(gap0['beta_at_min_gap']):g} between sorted roots "
                f"{int(gap0['lower_sorted_index'])}-{int(gap0['upper_sorted_index'])}. "
                "So the statement 'little/no veering at mu=0' is not supported by this sorted-order diagnostic."
            )
        else:
            lines.append("- At mu=0, some tracked branches change sorted index, so the no-veering interpretation needs manual review.")

    nonzero = [float(mu) for mu in args.mus if abs(float(mu)) > 1e-12]
    visible = [
        mu
        for mu in nonzero
        if int(exchange_metrics.get(mu, {}).get("branches_with_index_change", 0)) > 0
    ]
    if visible:
        visible_text = ", ".join(f"{mu:g}" for mu in visible)
        lines.append(f"- Fixed nonzero mu cases show branch approach or sorted-index exchange for: {visible_text}.")
    else:
        lines.append("- Fixed nonzero mu cases do not show clear sorted-index exchange in the focused branches on this grid; they mainly show smooth branch separation as beta changes.")

    nonzero_with_gaps = {mu: gap_metrics[mu] for mu in nonzero if mu in gap_metrics}
    if nonzero_with_gaps:
        strongest_nonzero_mu = min(nonzero_with_gaps, key=lambda key: float(nonzero_with_gaps[key]["min_gap"]))
        lines.append(
            f"- Among mu != 0 panels, the smallest adjacent gap is at mu={strongest_nonzero_mu:g}, "
            f"beta={float(nonzero_with_gaps[strongest_nonzero_mu]['beta_at_min_gap']):g}; "
            "larger mu values are less suggestive of beta-driven branch exchange in this run."
        )

    if gap_metrics:
        strongest_mu = min(gap_metrics, key=lambda key: float(gap_metrics[key]["min_gap"]))
        lines.append(
            f"- The smallest adjacent-gap event in this run is at mu={strongest_mu:g}, "
            f"beta={float(gap_metrics[strongest_mu]['beta_at_min_gap']):g}, "
            f"between sorted roots {int(gap_metrics[strongest_mu]['lower_sorted_index'])} and "
            f"{int(gap_metrics[strongest_mu]['upper_sorted_index'])}."
        )
    lines.append(
        "- Because mu is fixed in each panel, beta mainly probes coupling-strength changes; the observed beta sweeps are useful diagnostics, "
        "but they do not replace mu sweeps as the main detuning scan."
    )
    if warnings:
        lines.extend(["", "## Tracking Warnings", ""])
        for warning in warnings[:12]:
            lines.append(f"- {warning}")
        if len(warnings) > 12:
            lines.append(f"- ... {len(warnings) - 12} more warnings")
    return lines


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plot analytic Lambda(beta) branches at fixed mu values.")
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--mus", nargs="+", type=float, default=DEFAULT_MUS)
    parser.add_argument("--beta-min", type=float, default=DEFAULT_BETA_MIN)
    parser.add_argument("--beta-max", type=float, default=DEFAULT_BETA_MAX)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--num-roots", type=int, default=DEFAULT_NUM_ROOTS)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)
    args.output_dir = Path(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.mus = [float(mu) for mu in args.mus]
    betas = beta_grid(args.beta_min, args.beta_max, args.beta_step)
    eps_token = number_token(float(args.epsilon))
    args.sorted_csv = args.output_dir / f"lambda_beta_fixed_mu_eps{eps_token}.csv"
    args.branch_csv = args.output_dir / f"lambda_beta_fixed_mu_branch_paths_eps{eps_token}.csv"
    report_path = args.output_dir / "lambda_beta_fixed_mu_report.md"

    sorted_rows: list[dict[str, float | int | str]] = []
    branch_rows: list[dict[str, float | int | str]] = []
    plot_paths: list[Path] = []
    tracking_status: dict[float, str] = {}
    warnings: list[str] = []

    for mu in args.mus:
        result, status, mu_warnings = run_tracking_for_mu(
            epsilon=float(args.epsilon),
            mu=float(mu),
            betas=betas,
            num_roots=int(args.num_roots),
        )
        tracking_status[float(mu)] = status
        warnings.extend(mu_warnings)
        mu_sorted_rows, mu_branch_rows = collect_rows_for_mu(
            result=result,
            epsilon=float(args.epsilon),
            mu=float(mu),
            betas=betas,
            num_roots=int(args.num_roots),
        )
        sorted_rows.extend(mu_sorted_rows)
        branch_rows.extend(mu_branch_rows)
        plot_paths.append(
            plot_mu(
                result=result,
                output_dir=args.output_dir,
                epsilon=float(args.epsilon),
                mu=float(mu),
                betas=betas,
                num_roots=int(args.num_roots),
            )
        )

    write_csv(
        args.sorted_csv,
        sorted_rows,
        ["epsilon", "mu", "beta", "sorted_index", "branch_id", "Lambda", "warning_flag", "min_MAC"],
    )
    write_csv(
        args.branch_csv,
        branch_rows,
        ["epsilon", "mu", "beta", "branch_id", "current_sorted_index", "Lambda", "warning_flag", "min_MAC"],
    )
    report_path.write_text(
        "\n".join(
            report_lines(
                args=args,
                sorted_rows=sorted_rows,
                branch_rows=branch_rows,
                plot_paths=plot_paths,
                tracking_status=tracking_status,
                warnings=warnings,
            )
        )
        + "\n",
        encoding="utf-8",
    )

    print(f"saved sorted CSV: {args.sorted_csv}")
    print(f"saved branch CSV: {args.branch_csv}")
    for path in plot_paths:
        print(f"saved plot: {path}")
    print(f"saved report: {report_path}")


if __name__ == "__main__":
    main()
