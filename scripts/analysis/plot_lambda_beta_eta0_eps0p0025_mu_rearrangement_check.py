from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis import audit_mu_scan_eta0_first6_rearrangement as mu_audit  # noqa: E402


EPSILON = 0.0025
ETA = 0.0
MU_VALUES = (0.001, 0.002, 0.003)
BETA_STEP = 0.1

N_DESCENDANTS_PLOT = 6
N_SORTED_GAP_ROOTS = 7

SPIKE_ABSOLUTE_RESIDUAL = 0.02
SPIKE_RESIDUAL_RATIO = 50.0
SPIKE_ABSOLUTE_STEP = 0.05
SPIKE_STEP_RATIO = 25.0
REAL_CROSSING_GAP_TOL = 1e-4

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "lambda_beta_eta0_eps0p0025_mu_rearrangement_check.csv"
OUTPUT_REPORT = OUTPUT_DIR / "lambda_beta_eta0_eps0p0025_mu_rearrangement_check_report.md"

SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

CSV_FIELDNAMES = [
    "mu",
    "eta",
    "epsilon",
    "beta_deg",
    "descendant_id",
    "Lambda",
    "sorted_position",
    "mac_to_previous",
    "tracking_warning",
    "gap_sorted1_2",
    "gap_sorted2_3",
    "gap_sorted3_4",
    "gap_sorted4_5",
    "gap_sorted5_6",
    "gap_sorted6_7",
    "min_adjacent_gap_first7_at_beta",
    "rearrangement_class",
    "spike_artifact_detected_for_descendant",
]


@dataclass(frozen=True)
class LambdaBetaCase:
    mu: float
    scan_case: mu_audit.MuScanCase
    output_png: Path
    spike_rows: list[dict[str, float | int | str]]
    smoothness_rows: list[dict[str, float | int | str]]
    pair_rows: list[dict[str, float | int | str]]


def number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace(".", "p")


def number_text(value: float) -> str:
    return f"{float(value):.12g}"


def beta_values() -> np.ndarray:
    values = np.arange(0.0, 90.0 + 0.5 * BETA_STEP, BETA_STEP, dtype=float)
    values[0] = 0.0
    values[-1] = 90.0
    return np.unique(np.round(values, 12))


def output_png_for_mu(mu: float) -> Path:
    return OUTPUT_DIR / f"lambda_beta_eta0_eps0p0025_mu{number_token(mu)}_first6.png"


def finite_values(values: Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    return array[np.isfinite(array)]


def finite_min(values: Sequence[float]) -> float:
    finite = finite_values(values)
    return float(np.min(finite)) if finite.size else float("nan")


def finite_median(values: Sequence[float]) -> float:
    finite = finite_values(values)
    return float(np.median(finite)) if finite.size else float("nan")


def rows_by_beta_descendant(result) -> dict[tuple[float, int], dict[str, float | int | str]]:
    return {
        (round(float(row["beta_deg"]), 10), int(row["descendant_id"])): row
        for row in result.rows
    }


def tracking_warning(row: Mapping[str, float | int | str]) -> str:
    text = mu_audit.tracking_warning(row)
    return text if text else "no"


def adjacent_gaps(result) -> np.ndarray:
    return np.diff(result.sorted_roots[:, :N_SORTED_GAP_ROOTS], axis=1)


def gap_at_row(gaps: np.ndarray, beta_idx: int, gap_idx: int) -> float:
    return float(gaps[beta_idx, gap_idx])


def detect_spike_artifacts(scan_case: mu_audit.MuScanCase) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    result = scan_case.result
    spike_rows: list[dict[str, float | int | str]] = []
    smoothness_rows: list[dict[str, float | int | str]] = []
    beta_grid = np.asarray(result.beta_values_deg, dtype=float)

    for desc_idx in range(N_DESCENDANTS_PLOT):
        values = np.asarray(result.tracked_lambdas[desc_idx], dtype=float)
        steps = np.abs(np.diff(values))
        step_median = finite_median(steps)
        max_step = float(np.nanmax(steps)) if steps.size else float("nan")

        residuals = np.abs(values[1:-1] - 0.5 * (values[:-2] + values[2:]))
        residual_median = finite_median(residuals)
        max_residual = float(np.nanmax(residuals)) if residuals.size else float("nan")

        step_threshold = max(SPIKE_ABSOLUTE_STEP, SPIKE_STEP_RATIO * max(step_median, 1e-12))
        residual_threshold = max(
            SPIKE_ABSOLUTE_RESIDUAL,
            SPIKE_RESIDUAL_RATIO * max(residual_median, 1e-12),
        )
        step_flag = bool(np.isfinite(max_step) and max_step > step_threshold)
        residual_flag = bool(np.isfinite(max_residual) and max_residual > residual_threshold)
        row = {
            "mu": float(scan_case.mu),
            "descendant_id": desc_idx + 1,
            "max_adjacent_step": max_step,
            "median_adjacent_step": step_median,
            "max_spike_residual": max_residual,
            "median_spike_residual": residual_median,
            "step_threshold": step_threshold,
            "residual_threshold": residual_threshold,
            "large_smooth_step_detected": "yes" if step_flag else "no",
            "spike_artifact_detected": "yes" if residual_flag else "no",
        }
        smoothness_rows.append(row)

        if residual_flag and residuals.size:
            residual_index = int(np.nanargmax(residuals)) + 1
            spike_rows.append(
                {
                    **row,
                    "spike_type": "isolated_midpoint_residual",
                    "beta_deg": float(beta_grid[residual_index]),
                }
            )
    return spike_rows, smoothness_rows


def pair_crossing_diagnostics(scan_case: mu_audit.MuScanCase) -> list[dict[str, float | int | str]]:
    result = scan_case.result
    beta_grid = np.asarray(result.beta_values_deg, dtype=float)
    rows: list[dict[str, float | int | str]] = []
    for left in range(N_DESCENDANTS_PLOT):
        for right in range(left + 1, N_DESCENDANTS_PLOT):
            diff = np.asarray(result.tracked_lambdas[right] - result.tracked_lambdas[left], dtype=float)
            abs_diff = np.abs(diff)
            min_idx = int(np.nanargmin(abs_diff))
            nonzero_signs = np.sign(diff[np.abs(diff) > REAL_CROSSING_GAP_TOL])
            sign_change = bool(np.any(nonzero_signs[:-1] * nonzero_signs[1:] < 0)) if nonzero_signs.size > 1 else False
            rows.append(
                {
                    "mu": float(scan_case.mu),
                    "descendant_left": left + 1,
                    "descendant_right": right + 1,
                    "min_abs_gap": float(abs_diff[min_idx]),
                    "beta_min_abs_gap": float(beta_grid[min_idx]),
                    "signed_gap_changes_sign": "yes" if sign_change else "no",
                    "resolved_real_crossing": "yes" if float(abs_diff[min_idx]) <= REAL_CROSSING_GAP_TOL else "no",
                }
            )
    return rows


def min_sorted_gap_summary(scan_case: mu_audit.MuScanCase) -> dict[str, float | int]:
    gaps = adjacent_gaps(scan_case.result)
    min_by_beta = np.nanmin(gaps, axis=1)
    global_idx = int(np.nanargmin(min_by_beta))
    pair_idx = int(np.nanargmin(gaps[global_idx]))
    beta_grid = np.asarray(scan_case.result.beta_values_deg, dtype=float)
    return {
        "gap": float(min_by_beta[global_idx]),
        "beta": float(beta_grid[global_idx]),
        "pair_left": pair_idx + 1,
        "pair_right": pair_idx + 2,
    }


def min_adjacent_descendant_summary(scan_case: mu_audit.MuScanCase) -> dict[str, float | int]:
    pair_rows = pair_crossing_diagnostics(scan_case)
    adjacent_rows = [
        row for row in pair_rows if int(row["descendant_right"]) == int(row["descendant_left"]) + 1
    ]
    best = min(adjacent_rows, key=lambda row: float(row["min_abs_gap"]))
    return {
        "gap": float(best["min_abs_gap"]),
        "beta": float(best["beta_min_abs_gap"]),
        "pair_left": int(best["descendant_left"]),
        "pair_right": int(best["descendant_right"]),
    }


def warning_counts(scan_case: mu_audit.MuScanCase) -> dict[str, int]:
    rows = scan_case.result.rows
    return {
        "total": len(scan_case.result.warning_rows),
        "low_mac": sum(1 for row in rows if str(row.get("low_mac", "no")) == "yes"),
        "low_margin": sum(1 for row in rows if str(row.get("low_margin", "no")) == "yes"),
        "unresolved": sum(1 for row in rows if str(row.get("unresolved_assignment", "no")) == "yes"),
        "frequency_mac_disagreement": sum(
            1 for row in rows if str(row.get("frequency_mac_disagreement", "no")) == "yes"
        ),
    }


def warning_beta_summary(scan_case: mu_audit.MuScanCase, *, limit: int = 8) -> str:
    betas = sorted({round(float(row["beta_deg"]), 6) for row in scan_case.result.warning_rows})
    if not betas:
        return "none"
    text = ", ".join(f"{beta:g}" for beta in betas[:limit])
    if len(betas) > limit:
        text += f", ... ({len(betas)} beta values)"
    return text


def plot_case(case: LambdaBetaCase) -> None:
    result = case.scan_case.result
    colors = plt.get_cmap("tab10").colors
    fig, (ax_lambda, ax_pos) = plt.subplots(
        2,
        1,
        figsize=(10.4, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.15]},
        constrained_layout=True,
    )

    row_lookup = rows_by_beta_descendant(result)
    for desc_idx in range(N_DESCENDANTS_PLOT):
        color = colors[desc_idx]
        beta = np.asarray(result.beta_values_deg, dtype=float)
        lambdas = np.asarray(result.tracked_lambdas[desc_idx], dtype=float)
        ax_lambda.plot(beta, lambdas, color=color, lw=1.55, label=f"desc {desc_idx + 1}")
        warning_beta = []
        warning_lambda = []
        for beta_idx, beta_deg in enumerate(beta):
            row = row_lookup[(round(float(beta_deg), 10), desc_idx + 1)]
            if tracking_warning(row) != "no":
                warning_beta.append(float(beta_deg))
                warning_lambda.append(float(lambdas[beta_idx]))
        if warning_beta:
            ax_lambda.scatter(
                warning_beta,
                warning_lambda,
                marker="x",
                s=16,
                color=color,
                linewidths=0.9,
                alpha=0.8,
            )
        ax_pos.step(
            beta,
            result.current_sorted_positions[desc_idx],
            where="post",
            color=color,
            lw=1.2,
            label=f"desc {desc_idx + 1}",
        )

    sorted_gap = min_sorted_gap_summary(case.scan_case)
    ax_lambda.axvline(float(sorted_gap["beta"]), color="0.25", lw=0.8, ls=":", alpha=0.8)
    ax_lambda.text(
        float(sorted_gap["beta"]),
        0.03,
        f"min sorted gap {float(sorted_gap['gap']):.4g}",
        transform=ax_lambda.get_xaxis_transform(),
        ha="left",
        va="bottom",
        fontsize=8,
        color="0.22",
        rotation=90,
    )

    ax_lambda.set_ylabel(r"$\Lambda$")
    ax_lambda.set_title(
        rf"Tracked descendant branches: $\eta=0$, $\epsilon=0.0025$, $\mu={case.mu:g}$"
    )
    ax_lambda.grid(True, color="0.86", lw=0.55)
    ax_lambda.legend(loc="best", ncol=3, fontsize=8, frameon=False)

    ax_pos.set_xlabel(r"$\beta$ (deg)")
    ax_pos.set_ylabel("sorted position")
    ax_pos.set_xlim(0.0, 90.0)
    ax_pos.set_yticks(range(1, N_DESCENDANTS_PLOT + 1))
    ax_pos.grid(True, color="0.88", lw=0.5)

    case.output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(case.output_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_csv(cases: Sequence[LambdaBetaCase]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for case in cases:
            result = case.scan_case.result
            row_lookup = rows_by_beta_descendant(result)
            gaps = adjacent_gaps(result)
            spike_descendants = {int(row["descendant_id"]) for row in case.spike_rows}
            for beta_idx, beta_deg in enumerate(result.beta_values_deg):
                beta_gaps = gaps[beta_idx]
                for desc_idx in range(N_DESCENDANTS_PLOT):
                    row = row_lookup[(round(float(beta_deg), 10), desc_idx + 1)]
                    mac = float(row["mac_to_previous"])
                    writer.writerow(
                        {
                            "mu": number_text(case.mu),
                            "eta": number_text(ETA),
                            "epsilon": number_text(EPSILON),
                            "beta_deg": number_text(float(beta_deg)),
                            "descendant_id": desc_idx + 1,
                            "Lambda": number_text(float(result.tracked_lambdas[desc_idx, beta_idx])),
                            "sorted_position": int(result.current_sorted_positions[desc_idx, beta_idx]),
                            "mac_to_previous": "" if not np.isfinite(mac) else number_text(mac),
                            "tracking_warning": tracking_warning(row),
                            "gap_sorted1_2": number_text(gap_at_row(gaps, beta_idx, 0)),
                            "gap_sorted2_3": number_text(gap_at_row(gaps, beta_idx, 1)),
                            "gap_sorted3_4": number_text(gap_at_row(gaps, beta_idx, 2)),
                            "gap_sorted4_5": number_text(gap_at_row(gaps, beta_idx, 3)),
                            "gap_sorted5_6": number_text(gap_at_row(gaps, beta_idx, 4)),
                            "gap_sorted6_7": number_text(gap_at_row(gaps, beta_idx, 5)),
                            "min_adjacent_gap_first7_at_beta": number_text(float(np.nanmin(beta_gaps))),
                            "rearrangement_class": case.scan_case.summary_row["rearrangement_class"],
                            "spike_artifact_detected_for_descendant": (
                                "yes" if desc_idx + 1 in spike_descendants else "no"
                            ),
                        }
                    )


def rearrangement_sentence(case: LambdaBetaCase) -> str:
    cls = case.scan_case.summary_row["rearrangement_class"]
    changed = case.scan_case.summary_row["descendants_with_position_changes"]
    if cls == "true_beta_rearrangement":
        return f"yes, data classify true beta rearrangement involving descendants {changed}"
    if cls == "no_rearrangement":
        return "no, first-six descendants keep sorted positions on the checked beta range"
    return f"classified as {cls}, changed descendants {changed}"


def crossing_sentence(case: LambdaBetaCase) -> str:
    sorted_gap = min_sorted_gap_summary(case.scan_case)
    pair_rows = case.pair_rows
    resolved = [row for row in pair_rows if str(row["resolved_real_crossing"]) == "yes"]
    sign_changes = [row for row in pair_rows if str(row["signed_gap_changes_sign"]) == "yes"]
    if resolved:
        pairs = ", ".join(
            f"desc {int(row['descendant_left'])}/{int(row['descendant_right'])}"
            for row in resolved
        )
        return f"resolved zero-gap crossing candidate for {pairs}"
    if sign_changes:
        pairs = ", ".join(
            f"desc {int(row['descendant_left'])}/{int(row['descendant_right'])}"
            for row in sign_changes[:4]
        )
        return (
            f"descendant ordering changes for {pairs}, but the sorted spectrum keeps a finite "
            f"minimum gap {float(sorted_gap['gap']):.6g}; treat as avoided-crossing/veering-like, not a resolved real crossing"
        )
    return (
        f"no descendant sign-change crossing and no zero sorted gap; minimum sorted gap "
        f"{float(sorted_gap['gap']):.6g}"
    )


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def write_report(cases: Sequence[LambdaBetaCase], beta_grid: np.ndarray) -> None:
    all_spikes = [row for case in cases for row in case.spike_rows]
    summary_rows: list[list[str]] = []
    for case in cases:
        sorted_gap = min_sorted_gap_summary(case.scan_case)
        adjacent_desc_gap = min_adjacent_descendant_summary(case.scan_case)
        counts = warning_counts(case.scan_case)
        summary_rows.append(
            [
                number_text(case.mu),
                case.scan_case.summary_row["rearrangement_class"],
                case.scan_case.summary_row["descendants_with_position_changes"],
                f"{float(sorted_gap['gap']):.6g} at beta={float(sorted_gap['beta']):g} (sorted {int(sorted_gap['pair_left'])}-{int(sorted_gap['pair_right'])})",
                f"{float(adjacent_desc_gap['gap']):.6g} at beta={float(adjacent_desc_gap['beta']):g} (desc {int(adjacent_desc_gap['pair_left'])}-{int(adjacent_desc_gap['pair_right'])})",
                str(counts["total"]),
                number_text(float(case.scan_case.summary_row["min_mac_continuity"])),
                "yes" if case.spike_rows else "no",
            ]
        )

    lines = [
        "# Lambda(beta) Eta=0 Mu Rearrangement Check",
        "",
        "## What Was Built",
        "",
        "Diagnostic-only analytic Euler-Bernoulli plots of tracked descendant",
        "`Lambda(beta)` branches for `mu=0.001`, `0.002`, and `0.003` at",
        "`eta=0` and `epsilon=0.0025`. The plots are not article figures.",
        "",
        "The plotted curves are descendant branches tracked by analytic shape MAC,",
        "reusing `scripts/analysis/audit_mu_scan_eta0_first6_rearrangement.py`.",
        "They are not naive sorted-root curves.",
        "",
        "## Parameters And Grid",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- mu values: {', '.join(f'{mu:g}' for mu in MU_VALUES)}",
        f"- eta: {ETA:g}",
        f"- epsilon: {EPSILON:g}",
        "- beta range: 0..90 deg",
        f"- beta step: {BETA_STEP:g} deg",
        f"- beta samples: {len(beta_grid)}",
        f"- sorted roots solved per beta: first {mu_audit.N_SORTED_ROOTS}",
        f"- roots used as MAC tracking candidates: first {mu_audit.N_TRACKING_CANDIDATE_ROOTS}",
        f"- tracked descendants: first {mu_audit.N_DESCENDANTS_TRACK}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- shape samples per arm component: {mu_audit.NUM_SHAPE_SAMPLES}",
        "",
        "## Summary",
        "",
    ]
    lines.extend(
        markdown_table(
            [
                "mu",
                "class",
                "changed descendants",
                "min sorted gap first7",
                "min adjacent descendant gap",
                "tracking warnings",
                "min MAC",
                "spike artifacts",
            ],
            summary_rows,
        )
    )
    lines.extend(["", "## Per-Mu Interpretation", ""])
    for case in cases:
        counts = warning_counts(case.scan_case)
        lines.extend(
            [
                f"### mu={case.mu:g}",
                "",
                f"- rearrangement: {rearrangement_sentence(case)}",
                f"- crossing/veering: {crossing_sentence(case)}",
                (
                    f"- tracking difficulty: {counts['total']} warning rows "
                    f"(low_mac={counts['low_mac']}, low_margin={counts['low_margin']}, "
                    f"frequency/MAC disagreement={counts['frequency_mac_disagreement']}); "
                    f"warning beta values: {warning_beta_summary(case.scan_case)}"
                ),
                (
                    "- spike/tracking artifacts: "
                    + (
                        "detected by isolated midpoint residual check"
                        if case.spike_rows
                        else "no isolated midpoint spike artifacts detected"
                    )
                ),
                f"- output PNG: `{case.output_png.relative_to(REPO_ROOT)}`",
                "",
            ]
        )
    lines.extend(
        [
            "## Spike And Tracking Artifact Handling",
            "",
            f"Spike check result: `{'spike artifacts detected' if all_spikes else 'no spike artifacts detected'}`.",
            "The script uses a finer beta step than the previous mu scan (`0.1 deg`",
            "instead of `0.25 deg`) for these three local checks. Warning points are",
            "marked directly on the plots with small `x` markers. No curve is built",
            "by switching to naive sorted roots.",
            "",
            "This refined local check does not reproduce the previous coarse-grid",
            "`0.001 <= mu <= 0.002` true-rearrangement classification from the",
            "`0.25 deg` mu scan. Treat that earlier narrow interval as",
            "grid/tracking-sensitive diagnostic evidence unless it is confirmed by",
            "a stricter localized subspace-MAC or root-coalescence analysis.",
            "",
            "## Limitations",
            "",
            "- EB analytic thickness-mismatch determinant only, evaluated at eta=0",
            "- no Timoshenko model",
            "- no FEM, Gmsh, CalculiX, or heavy 3D workflow",
            "- finite positive gaps are grid-resolved diagnostics, not exact analytic lower bounds",
            "- near multiple roots, a stricter proof of crossing vs avoided crossing may require localized subspace-MAC or symbolic/numerical root coalescence analysis",
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
        ]
    )
    for case in cases:
        lines.append(f"- PNG mu={case.mu:g}: `{case.output_png.relative_to(REPO_ROOT)}`")
    lines.append("")
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def run_case(mu: float, beta_grid: np.ndarray) -> LambdaBetaCase:
    print(f"running Lambda(beta) rearrangement check for mu={float(mu):g}")
    scan_case = mu_audit.run_mu_case(float(mu), beta_grid)
    spike_rows, smoothness_rows = detect_spike_artifacts(scan_case)
    pair_rows = pair_crossing_diagnostics(scan_case)
    case = LambdaBetaCase(
        mu=float(mu),
        scan_case=scan_case,
        output_png=output_png_for_mu(float(mu)),
        spike_rows=spike_rows,
        smoothness_rows=smoothness_rows,
        pair_rows=pair_rows,
    )
    plot_case(case)
    return case


def main() -> dict[str, object]:
    print("Lambda(beta) eta=0 mu rearrangement check")
    beta_grid = beta_values()
    cases = [run_case(float(mu), beta_grid) for mu in MU_VALUES]
    write_csv(cases)
    write_report(cases, beta_grid)
    print(f"saved CSV: {OUTPUT_CSV}")
    print(f"saved report: {OUTPUT_REPORT}")
    for case in cases:
        print(f"saved PNG: {case.output_png}")
    return {"cases": cases, "csv": OUTPUT_CSV, "report": OUTPUT_REPORT}


if __name__ == "__main__":
    main()
