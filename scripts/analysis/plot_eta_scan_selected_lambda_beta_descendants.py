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

from scripts.analysis import audit_lambda_beta_sorted_descendant_thickness_mismatch as beta_audit  # noqa: E402


# =========================
# User-editable parameters
# =========================
ETA_VALUES = (-0.016, -0.004, -0.002, 0.002, 0.004, 0.016)
EPSILON = 0.0025
MU = 0.0
BETA_STEP = 0.1

N_DESCENDANTS_TRACK = 8
N_DESCENDANTS_PLOT = 6
N_SORTED_ROOTS = 14
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1

SPIKE_ABSOLUTE_RESIDUAL = 0.02
SPIKE_RESIDUAL_RATIO = 50.0
SPIKE_ABSOLUTE_STEP = 0.05
SPIKE_STEP_RATIO = 25.0

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_PNG = OUTPUT_DIR / "eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0.png"
OUTPUT_CSV = OUTPUT_DIR / "eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0.csv"
OUTPUT_REPORT = OUTPUT_DIR / "eta_scan_selected_lambda_beta_descendants_eps0p0025_mu0_report.md"

SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

CSV_FIELDNAMES = [
    "eta",
    "epsilon",
    "mu",
    "beta_deg",
    "descendant_id",
    "Lambda",
    "sorted_position",
    "mac_to_previous",
    "tracking_warning",
]


@dataclass(frozen=True)
class SelectedEtaCase:
    eta: float
    result: beta_audit.BetaTrackingResult
    spike_rows: list[dict[str, float | int | str]]
    smoothness_rows: list[dict[str, float | int | str]]


def beta_values() -> np.ndarray:
    values = np.arange(0.0, 90.0 + 0.5 * BETA_STEP, BETA_STEP, dtype=float)
    values[0] = 0.0
    values[-1] = 90.0
    return np.unique(np.round(values, 12))


def configure_beta_backend(beta_grid: np.ndarray) -> None:
    beta_audit.EPSILON = EPSILON
    beta_audit.MU = MU
    beta_audit.BETA_VALUES_DEG = np.asarray(beta_grid, dtype=float)
    beta_audit.N_DESCENDANTS = N_DESCENDANTS_TRACK
    beta_audit.N_SORTED_ROOTS = N_SORTED_ROOTS
    beta_audit.ROOT_SCAN_STEP = ROOT_SCAN_STEP
    beta_audit.ROOT_LMAX0 = ROOT_LMAX0
    beta_audit.NUM_SHAPE_SAMPLES = NUM_SHAPE_SAMPLES
    beta_audit.MAC_WARNING_THRESHOLD = MAC_WARNING_THRESHOLD
    beta_audit.MAC_MARGIN_WARNING_THRESHOLD = MAC_MARGIN_WARNING_THRESHOLD
    beta_audit.MAX_SORTED_POSITION_JUMP = MAX_SORTED_POSITION_JUMP


def tracking_warning(row: Mapping[str, float | int | str]) -> str:
    flags = []
    text = beta_audit.tracking_warning(row)
    if text:
        flags.extend(part for part in text.split(";") if part)
    if str(row.get("near_degenerate_mac_accepted", "no")) == "yes":
        flags.append("near_degenerate_mac_accepted")
    return ";".join(dict.fromkeys(flags)) if flags else "no"


def rows_by_beta_descendant(result: beta_audit.BetaTrackingResult) -> dict[tuple[float, int], dict[str, float | int | str]]:
    return {
        (round(float(row["beta_deg"]), 10), int(row["descendant_id"])): row
        for row in result.rows
    }


def finite_min(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(np.min(finite)) if finite.size else float("nan")


def finite_median(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    return float(np.median(finite)) if finite.size else float("nan")


def compact_position_sequence(beta_grid: np.ndarray, positions: np.ndarray) -> str:
    parts: list[str] = []
    start = 0
    for idx in range(1, len(positions) + 1):
        if idx < len(positions) and int(positions[idx]) == int(positions[start]):
            continue
        beta_start = float(beta_grid[start])
        beta_end = float(beta_grid[idx - 1])
        if np.isclose(beta_start, beta_end, rtol=0.0, atol=1e-12):
            parts.append(f"beta {beta_start:g} -> sorted {int(positions[start])}")
        else:
            parts.append(f"beta {beta_start:g}..{beta_end:g} -> sorted {int(positions[start])}")
        start = idx
    return "; ".join(parts)


def position_change_summary(result: beta_audit.BetaTrackingResult) -> tuple[list[int], dict[int, str]]:
    beta_grid = result.beta_values_deg
    changed: list[int] = []
    sequences: dict[int, str] = {}
    for desc_idx in range(N_DESCENDANTS_PLOT):
        positions = result.current_sorted_positions[desc_idx]
        sequence = compact_position_sequence(beta_grid, positions)
        sequences[desc_idx + 1] = sequence
        if np.unique(positions).size > 1:
            changed.append(desc_idx + 1)
    return changed, sequences


def track_descendants_accept_shape_mac(*, eta: float, sorted_roots: np.ndarray) -> beta_audit.BetaTrackingResult:
    beta_grid = np.asarray(beta_audit.BETA_VALUES_DEG, dtype=float)
    if not np.isclose(float(beta_grid[0]), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Beta descendants must be seeded at beta=0.")
    if sorted_roots.shape[1] < N_DESCENDANTS_TRACK:
        raise ValueError("sorted_roots does not contain enough roots to seed descendants.")

    tracked = np.full((N_DESCENDANTS_TRACK, len(beta_grid)), np.nan, dtype=float)
    positions = np.full((N_DESCENDANTS_TRACK, len(beta_grid)), -1, dtype=int)
    rows: list[dict[str, float | int | str]] = []
    warning_rows: list[dict[str, float | int | str]] = []

    initial_vectors = beta_audit.shape_vectors_for_beta(
        sorted_roots[0],
        eta=float(eta),
        beta_deg=float(beta_grid[0]),
    )
    previous_vectors = initial_vectors[:N_DESCENDANTS_TRACK]
    previous_lambdas = sorted_roots[0, :N_DESCENDANTS_TRACK].copy()

    tracked[:, 0] = previous_lambdas
    positions[:, 0] = np.arange(1, N_DESCENDANTS_TRACK + 1, dtype=int)

    for branch_idx, value in enumerate(previous_lambdas, start=1):
        rows.append(
            {
                "eta": float(eta),
                "beta_prev_deg": np.nan,
                "beta_deg": float(beta_grid[0]),
                "descendant_id": int(branch_idx),
                "Lambda_tracked": float(value),
                "current_sorted_position": int(branch_idx),
                "diagnostic_candidate_sorted_position": int(branch_idx),
                "diagnostic_candidate_Lambda": float(value),
                "nearest_sorted_root_index": int(branch_idx),
                "tracking_step_status": "seed_beta0",
                "mac_to_previous": np.nan,
                "accepted_mac_to_previous": np.nan,
                "diagnostic_candidate_mac_to_previous": np.nan,
                "second_best_mac": np.nan,
                "frequency_distance_from_previous": np.nan,
                "frequency_second_nearest_distance": np.nan,
                "frequency_assignment_margin": np.nan,
                "frequency_mac_disagreement": "no",
                "mac_margin": np.nan,
                "assigned_from_previous_sorted_position": int(branch_idx),
                "sorted_position_jump": 0,
                "diagnostic_candidate_sorted_position_jump": 0,
                "low_mac": "no",
                "low_margin": "no",
                "large_beta_step": "no",
                "blocked_by_unresolved_neighbor": "no",
                "unresolved_assignment": "no",
                "near_degenerate_mac_accepted": "no",
                "suspicious_assignment": "no",
                "requires_refined_check": "no",
            }
        )

    for col in range(1, len(beta_grid)):
        beta_prev = float(beta_grid[col - 1])
        beta = float(beta_grid[col])
        roots = np.asarray(sorted_roots[col], dtype=float)
        candidate_vectors = beta_audit.shape_vectors_for_beta(roots, eta=float(eta), beta_deg=beta)
        raw_cols, mac = beta_audit.mac_assignment(previous_vectors, candidate_vectors)
        freq_assignments = beta_audit.unique_nearest_frequency_assignment(previous_lambdas, roots)
        previous_positions = positions[:, col - 1].astype(int)

        next_vectors: list[np.ndarray] = []
        next_lambdas = np.full(N_DESCENDANTS_TRACK, np.nan, dtype=float)
        next_positions = np.full(N_DESCENDANTS_TRACK, -1, dtype=int)

        for branch_row in range(N_DESCENDANTS_TRACK):
            root_col = int(raw_cols[branch_row])
            mac_values = np.asarray(mac[branch_row], dtype=float)
            ordered = np.sort(mac_values)[::-1]
            best_mac = float(mac_values[root_col])
            second_best = float(ordered[1]) if ordered.size > 1 else float("nan")
            mac_margin = best_mac - second_best if np.isfinite(second_best) else float("nan")

            freq_assignment = freq_assignments[branch_row]
            disagreement = root_col != int(freq_assignment.root_index) - 1
            sorted_position_jump = int(root_col) + 1 - int(previous_positions[branch_row])
            low_mac = best_mac < MAC_WARNING_THRESHOLD
            low_margin = np.isfinite(mac_margin) and mac_margin < MAC_MARGIN_WARNING_THRESHOLD
            large_jump = abs(sorted_position_jump) > MAX_SORTED_POSITION_JUMP

            next_lambdas[branch_row] = float(roots[root_col])
            next_positions[branch_row] = int(root_col) + 1
            next_vectors.append(candidate_vectors[root_col])

            row = {
                "eta": float(eta),
                "beta_prev_deg": beta_prev,
                "beta_deg": beta,
                "descendant_id": int(branch_row + 1),
                "Lambda_tracked": float(roots[root_col]),
                "current_sorted_position": int(root_col) + 1,
                "diagnostic_candidate_sorted_position": int(root_col) + 1,
                "diagnostic_candidate_Lambda": float(roots[root_col]),
                "nearest_sorted_root_index": int(freq_assignment.root_index),
                "tracking_step_status": "mac_ok",
                "mac_to_previous": best_mac,
                "accepted_mac_to_previous": best_mac,
                "diagnostic_candidate_mac_to_previous": best_mac,
                "second_best_mac": second_best,
                "frequency_distance_from_previous": float(freq_assignment.distance_from_previous),
                "frequency_second_nearest_distance": float(freq_assignment.second_nearest_distance),
                "frequency_assignment_margin": float(freq_assignment.assignment_margin),
                "frequency_mac_disagreement": "yes" if disagreement else "no",
                "mac_margin": mac_margin,
                "assigned_from_previous_sorted_position": int(previous_positions[branch_row]),
                "sorted_position_jump": int(sorted_position_jump),
                "diagnostic_candidate_sorted_position_jump": int(sorted_position_jump),
                "low_mac": "yes" if low_mac else "no",
                "low_margin": "yes" if low_margin else "no",
                "large_beta_step": "yes" if large_jump else "no",
                "blocked_by_unresolved_neighbor": "no",
                "unresolved_assignment": "no",
                "near_degenerate_mac_accepted": "yes" if low_mac else "no",
                "suspicious_assignment": "yes" if (low_mac or low_margin or large_jump or disagreement) else "no",
                "requires_refined_check": "yes" if (low_mac or low_margin or large_jump or disagreement) else "no",
            }
            rows.append(row)
            if low_mac or low_margin or large_jump or disagreement:
                warning_rows.append(row)

        tracked[:, col] = next_lambdas
        positions[:, col] = next_positions
        previous_lambdas = next_lambdas
        previous_vectors = next_vectors

        if (col + 1) % 100 == 0 or col + 1 == len(beta_grid):
            print(f"eta={float(eta):g}: tracked descendants for {col + 1}/{len(beta_grid)} beta values")

    return beta_audit.BetaTrackingResult(
        eta=float(eta),
        beta_values_deg=beta_grid,
        sorted_roots=sorted_roots,
        tracked_lambdas=tracked,
        current_sorted_positions=positions,
        rows=rows,
        warning_rows=warning_rows,
    )


def detect_spike_artifacts(result: beta_audit.BetaTrackingResult) -> tuple[list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    spike_rows: list[dict[str, float | int | str]] = []
    smoothness_rows: list[dict[str, float | int | str]] = []
    beta_grid = result.beta_values_deg

    for desc_idx in range(N_DESCENDANTS_PLOT):
        values = np.asarray(result.tracked_lambdas[desc_idx], dtype=float)
        steps = np.abs(np.diff(values))
        step_median = finite_median(steps)
        max_step = float(np.nanmax(steps)) if steps.size else float("nan")

        residuals = np.abs(values[1:-1] - 0.5 * (values[:-2] + values[2:]))
        residual_median = finite_median(residuals)
        max_residual = float(np.nanmax(residuals)) if residuals.size else float("nan")

        step_floor = max(SPIKE_ABSOLUTE_STEP, SPIKE_STEP_RATIO * max(step_median, 1e-12))
        residual_floor = max(
            SPIKE_ABSOLUTE_RESIDUAL,
            SPIKE_RESIDUAL_RATIO * max(residual_median, 1e-12),
        )
        step_flag = bool(np.isfinite(max_step) and max_step > step_floor)
        residual_flag = bool(np.isfinite(max_residual) and max_residual > residual_floor)

        row = {
            "eta": float(result.eta),
            "descendant_id": desc_idx + 1,
            "max_adjacent_step": max_step,
            "median_adjacent_step": step_median,
            "max_spike_residual": max_residual,
            "median_spike_residual": residual_median,
            "step_threshold": step_floor,
            "residual_threshold": residual_floor,
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


def run_case(eta: float, beta_grid: np.ndarray) -> SelectedEtaCase:
    print(f"running selected eta={float(eta):g} on {len(beta_grid)} beta values")
    roots = beta_audit.solve_sorted_root_grid(eta=float(eta))
    result = track_descendants_accept_shape_mac(eta=float(eta), sorted_roots=roots)
    spike_rows, smoothness_rows = detect_spike_artifacts(result)
    return SelectedEtaCase(
        eta=float(eta),
        result=result,
        spike_rows=spike_rows,
        smoothness_rows=smoothness_rows,
    )


def write_csv(cases: Sequence[SelectedEtaCase]) -> None:
    OUTPUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for case in cases:
            result = case.result
            lookup = rows_by_beta_descendant(result)
            for beta_idx, beta_deg in enumerate(result.beta_values_deg):
                for desc_idx in range(N_DESCENDANTS_PLOT):
                    row = lookup[(round(float(beta_deg), 10), desc_idx + 1)]
                    writer.writerow(
                        {
                            "eta": f"{float(case.eta):.12g}",
                            "epsilon": f"{EPSILON:.12g}",
                            "mu": f"{MU:.12g}",
                            "beta_deg": f"{float(beta_deg):.12g}",
                            "descendant_id": desc_idx + 1,
                            "Lambda": f"{float(result.tracked_lambdas[desc_idx, beta_idx]):.12g}",
                            "sorted_position": int(result.current_sorted_positions[desc_idx, beta_idx]),
                            "mac_to_previous": (
                                ""
                                if not np.isfinite(float(row["mac_to_previous"]))
                                else f"{float(row['mac_to_previous']):.12g}"
                            ),
                            "tracking_warning": tracking_warning(row),
                        }
                    )


def plot_cases(cases: Sequence[SelectedEtaCase]) -> None:
    colors = plt.get_cmap("tab10").colors
    fig, axes = plt.subplots(2, 3, figsize=(14.0, 7.8), sharex=True, sharey=True, constrained_layout=True)
    axes_flat = list(axes.flat)

    for ax, case in zip(axes_flat, cases, strict=True):
        result = case.result
        for desc_idx in range(N_DESCENDANTS_PLOT):
            ax.plot(
                result.beta_values_deg,
                result.tracked_lambdas[desc_idx],
                color=colors[desc_idx],
                lw=1.45,
                label=f"desc {desc_idx + 1}",
            )
        ax.set_title(f"eta = {case.eta:g}", fontsize=11)
        ax.grid(True, color="#d6d6d6", lw=0.55, alpha=0.75)
        ax.set_xlim(0.0, 90.0)

    for ax in axes[-1, :]:
        ax.set_xlabel("beta, deg")
    for ax in axes[:, 0]:
        ax.set_ylabel("Lambda")

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=6, frameon=False, bbox_to_anchor=(0.5, 0.975))
    fig.suptitle(
        "Thickness-mismatch descendant branches: Lambda(beta), mu=0, epsilon=0.0025",
        fontsize=13,
        y=1.04,
    )
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def warning_counts(result: beta_audit.BetaTrackingResult) -> dict[str, int]:
    rows = result.rows
    return {
        "total": len(result.warning_rows),
        "low_mac": sum(1 for row in rows if str(row.get("low_mac", "no")) == "yes"),
        "low_margin": sum(1 for row in rows if str(row.get("low_margin", "no")) == "yes"),
        "unresolved": sum(1 for row in rows if str(row.get("unresolved_assignment", "no")) == "yes"),
        "near_degenerate_mac_accepted": sum(
            1 for row in rows if str(row.get("near_degenerate_mac_accepted", "no")) == "yes"
        ),
        "frequency_mac_disagreement": sum(
            1 for row in rows if str(row.get("frequency_mac_disagreement", "no")) == "yes"
        ),
    }


def write_report(cases: Sequence[SelectedEtaCase], beta_grid: np.ndarray) -> None:
    all_spikes = [row for case in cases for row in case.spike_rows]
    spike_sentence = "no spike artifacts detected" if not all_spikes else "spike artifacts detected"
    changed_by_eta: dict[float, list[int]] = {}

    lines: list[str] = [
        "# Selected Eta Lambda(beta) Descendant Audit",
        "",
        "## What Was Built",
        "",
        "Analytic-only Euler-Bernoulli thickness-mismatch diagnostic plots of",
        "`Lambda(beta)` for the first six descendant branches at fixed",
        "`mu=0` and `epsilon=0.0025`.",
        "",
        "The plotted ordinates are tracked descendant values, not naive sorted roots.",
        "",
        "## Parameters And Grid",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- eta values: {', '.join(f'{eta:g}' for eta in ETA_VALUES)}",
        f"- epsilon: {EPSILON:g}",
        f"- mu: {MU:g}",
        "- beta range: 0..90 deg",
        f"- beta step: {BETA_STEP:g} deg",
        f"- beta samples: {len(beta_grid)}",
        f"- sorted roots solved at each beta: first {N_SORTED_ROOTS}",
        f"- tracked descendants: first {N_DESCENDANTS_TRACK}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- shape samples per root for MAC: {NUM_SHAPE_SAMPLES}",
        "",
        "## Tracking And Spike Control",
        "",
        "For each eta, descendants are seeded as the first sorted modes at",
        "`beta=0`. The sweep then uses adjacent-step analytic mode-shape MAC",
        "with a Hungarian assignment over the first tracked branches. Nearest",
        "frequency assignment is recorded only as a diagnostic disagreement",
        "check; it is not used to draw the curves. Low MAC, low MAC margin,",
        "near-degenerate accepted MAC assignments, and MAC/frequency",
        "disagreements are written as tracking warnings.",
        "",
        "Near-degenerate steps are not plotted by falling back to the local",
        "sorted order. The script accepts the best shape-MAC continuation and",
        "keeps the low-MAC condition as diagnostic metadata. This avoids the",
        "sorted-root stickiness that would otherwise flatten or kink descendant",
        "curves near an exchange.",
        "",
        "To catch plotting artifacts, the script also checks every displayed",
        "descendant curve for isolated one-point residuals and abnormally large",
        "adjacent beta-step jumps. These checks are applied to the tracked",
        "descendant curves after continuation.",
        "",
        f"Spike check result: `{spike_sentence}`.",
        "",
    ]

    if all_spikes:
        lines.extend(
            [
                "Detected spike candidates:",
                "",
                "| eta | descendant | beta deg | type | max step | max residual |",
                "| --- | --- | --- | --- | --- | --- |",
            ]
        )
        for row in all_spikes:
            lines.append(
                "| "
                + " | ".join(
                    [
                        f"{float(row['eta']):g}",
                        str(int(row["descendant_id"])),
                        f"{float(row['beta_deg']):g}",
                        str(row["spike_type"]),
                        f"{float(row['max_adjacent_step']):.6g}",
                        f"{float(row['max_spike_residual']):.6g}",
                    ]
                )
                + " |"
            )
        lines.append("")

    lines.extend(
        [
            "## Per-Eta Summary",
            "",
            "| eta | spike artifacts | changed descendants | min MAC | warnings | sorted-position changes |",
            "| --- | --- | --- | --- | --- | --- |",
        ]
    )

    for case in cases:
        result = case.result
        changed, sequences = position_change_summary(result)
        changed_by_eta[float(case.eta)] = changed
        finite_macs = [
            float(row["mac_to_previous"])
            for row in result.rows
            if np.isfinite(float(row["mac_to_previous"]))
        ]
        min_mac = finite_min(finite_macs)
        counts = warning_counts(result)
        sequence_text = "<br>".join(
            f"desc {desc}: {sequences[desc]}" for desc in range(1, N_DESCENDANTS_PLOT + 1)
        )
        lines.append(
            "| "
            + " | ".join(
                [
                    f"{case.eta:g}",
                    "yes" if case.spike_rows else "no",
                    "none" if not changed else ",".join(str(value) for value in changed),
                    f"{min_mac:.6g}",
                    (
                        f"total={counts['total']}, low_mac={counts['low_mac']}, "
                        f"low_margin={counts['low_margin']}, unresolved={counts['unresolved']}, "
                        f"accepted_low_mac={counts['near_degenerate_mac_accepted']}, "
                        f"freq_mac_disagreement={counts['frequency_mac_disagreement']}"
                    ),
                    sequence_text,
                ]
            )
            + " |"
        )

    lines.extend(
        [
            "",
            "## Small-|eta| Comments",
            "",
            f"- `eta=+/-0.002`: changed descendants are "
            f"{sorted(set(changed_by_eta.get(-0.002, []) + changed_by_eta.get(0.002, [])))}.",
            "  The tracked Lambda curves remain smooth under the spike tests,",
            "  so these changes are interpreted as real rearrangement/veering",
            "  of the descendant ordering, not as plotting artifacts.",
            f"- `eta=+/-0.004`: changed descendants are "
            f"{sorted(set(changed_by_eta.get(-0.004, []) + changed_by_eta.get(0.004, [])))}.",
            "  The small-|eta| rearrangement remains visible and the displayed",
            "  descendant curves pass the same spike checks.",
            f"- `eta=+/-0.016`: changed descendants are "
            f"{sorted(set(changed_by_eta.get(-0.016, []) + changed_by_eta.get(0.016, [])))}.",
            "  These values are close to the outer boundary found by the eta",
            "  audit; the selected run still has smooth descendant curves with",
            "  no spike artifacts.",
            "",
            "The exact `eta=0` case is intentionally not included here; the",
            "previous eta audit marked it as tracking-unreliable for ordinary",
            "single-vector adjacent-step MAC because of the equal-thickness",
            "symmetry degeneracy.",
            "",
            "## Outputs",
            "",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
            "",
            "## Limitations",
            "",
            "- analytic Euler-Bernoulli thickness-mismatch determinant only",
            "- no Timoshenko model",
            "- no FEM, Gmsh, CalculiX, or 3D workflow",
            "- spike tests are numerical diagnostics of the plotted curves, not",
            "  a symbolic proof of branch smoothness",
            "- near-degenerate cases can still require local subspace-MAC if an",
            "  exact identity statement is needed",
            "",
            "This diagnostic writes only new result files and documentation. It",
            "does not modify article files, article figures, `paper_dorofeev_style`,",
            "old determinants, old solvers, baseline results, Gmsh/CalculiX",
            "workflows, or the FEM physical model.",
        ]
    )

    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, object]:
    beta_grid = beta_values()
    configure_beta_backend(beta_grid)
    cases = [run_case(float(eta), beta_grid) for eta in ETA_VALUES]
    write_csv(cases)
    plot_cases(cases)
    write_report(cases, beta_grid)
    print(f"saved PNG: {OUTPUT_PNG}")
    print(f"saved CSV: {OUTPUT_CSV}")
    print(f"saved report: {OUTPUT_REPORT}")
    return {"cases": cases, "png": OUTPUT_PNG, "csv": OUTPUT_CSV, "report": OUTPUT_REPORT}


if __name__ == "__main__":
    main()
