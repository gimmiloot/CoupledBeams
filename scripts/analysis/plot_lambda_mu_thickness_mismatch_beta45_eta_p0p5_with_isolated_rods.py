from __future__ import annotations

from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib.thickness_mismatch_diagnostic_helpers import (  # noqa: E402
    isolated_clamped_supported_references,
    mu_grid,
    plot_with_validity_split,
    rods_label,
    roots_by_mu_eta,
    thickness_ratio_summary,
    track_descendants_from_mu0,
    tracking_warning_rows,
    tracking_warning_summary,
    valid_mask,
)


# =========================
# User-editable parameters
# =========================
BETA_DEG = 45.0
COMPARISON_BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.005

N_DESCENDANTS_PLOT = 6
N_DESCENDANTS_TRACK = 8
N_SORTED_SCAN = 12
N_ISOLATED_ROD_CURVES = 6
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
THICKNESS_RATIO_LIMIT = 0.1
Y_MAX = 18.0

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods.png"
OUTPUT_REPORT = OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_report.md"

MU_VALUES = mu_grid(MU_MIN, MU_MAX, MU_STEP)


def compute_tracking(beta_deg: float) -> tuple[np.ndarray, list[dict[str, float | int | str]]]:
    roots = roots_by_mu_eta(
        beta_rad=float(np.deg2rad(beta_deg)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        n_roots=N_SORTED_SCAN,
        root_lmax0=ROOT_LMAX0,
        root_scan_step=ROOT_SCAN_STEP,
    )
    result = track_descendants_from_mu0(
        beta_rad=float(np.deg2rad(beta_deg)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        roots_by_mu=roots,
        num_descendants=N_DESCENDANTS_TRACK,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
    )
    return result.tracked, result.rows


def reference_grid_distance_metrics(
    *,
    tracked: np.ndarray,
    rod1_refs: np.ndarray,
    rod2_refs: np.ndarray,
) -> dict[str, float]:
    references = np.vstack((rod1_refs, rod2_refs))
    values = np.asarray(tracked[:N_DESCENDANTS_PLOT], dtype=float)
    distances: list[float] = []
    for col in range(values.shape[1]):
        ref_column = references[:, col]
        for row in range(values.shape[0]):
            distances.append(float(np.min(np.abs(ref_column - values[row, col]))))
    arr = np.asarray(distances, dtype=float)
    return {
        "mean_nearest_distance": float(np.mean(arr)),
        "median_nearest_distance": float(np.median(arr)),
        "max_nearest_distance": float(np.max(arr)),
    }


def plot_with_isolated_rods(*, tracked: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    valid = valid_mask(epsilon=EPSILON, eta=ETA, mu_values=MU_VALUES, limit=THICKNESS_RATIO_LIMIT)
    alphas, rod1_refs, rod2_refs = isolated_clamped_supported_references(
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        n_curves=N_ISOLATED_ROD_CURVES,
    )

    fig, ax = plt.subplots(figsize=(9.8, 5.8))
    for idx in range(N_ISOLATED_ROD_CURVES):
        ax.plot(
            MU_VALUES,
            rod1_refs[idx],
            color="0.35",
            lw=1.05,
            ls=(0, (5.0, 2.2)),
            alpha=0.42,
            zorder=1,
        )
        ax.plot(
            MU_VALUES,
            rod2_refs[idx],
            color="0.55",
            lw=1.05,
            ls=(0, (2.0, 2.4)),
            alpha=0.42,
            zorder=1,
        )

    for branch in range(1, N_DESCENDANTS_PLOT + 1):
        plot_with_validity_split(
            ax,
            MU_VALUES,
            tracked[branch - 1],
            valid,
            color=colors[(branch - 1) % len(colors)],
            linewidth=2.0,
            label=f"desc {branch}",
            zorder=3,
        )

    handles = [
        Line2D([0], [0], color="black", lw=2.0, label="system descendant branches"),
        Line2D([0], [0], color="0.35", lw=1.1, ls=(0, (5.0, 2.2)), label="isolated rod 1, CS"),
        Line2D([0], [0], color="0.55", lw=1.1, ls=(0, (2.0, 2.4)), label="isolated rod 2, CS"),
        Line2D([0], [0], color="0.2", lw=1.8, ls="-", label=r"solid: $2r_i/l_i \leq 0.1$"),
        Line2D([0], [0], color="0.2", lw=1.8, ls="--", label="dashed system segment: criterion violated"),
    ]
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_title(
        rf"Thickness-mismatch descendants with isolated CS rods, "
        rf"$\eta={ETA:g}$, $\beta={BETA_DEG:g}^\circ$"
    )
    ax.set_ylim(0.0, Y_MAX)
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(handles=handles, loc="upper left", fontsize=8.8, frameon=False)
    fig.tight_layout()
    fig.savefig(OUTPUT_PNG, dpi=240, bbox_inches="tight")
    plt.close(fig)
    return alphas, rod1_refs, rod2_refs


def write_report(
    *,
    tracking_rows: list[dict[str, float | int | str]],
    alphas: np.ndarray,
    rod1_refs: np.ndarray,
    rod2_refs: np.ndarray,
    beta45_metrics: dict[str, float],
    beta15_metrics: dict[str, float],
) -> dict[str, object]:
    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    summary = thickness_ratio_summary(
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        limit=THICKNESS_RATIO_LIMIT,
    )
    warning_rows = tracking_warning_rows(tracking_rows)
    warning_summary = tracking_warning_summary(tracking_rows)
    mean_relation = "smaller" if beta45_metrics["mean_nearest_distance"] < beta15_metrics["mean_nearest_distance"] else "larger"
    median_relation = (
        "smaller" if beta45_metrics["median_nearest_distance"] < beta15_metrics["median_nearest_distance"] else "larger"
    )
    lines = [
        "# Eta=0.5 Lambda(mu) With Isolated-Rod References, beta=45",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- comparison beta for the diagnostic observation: {COMPARISON_BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu range: {MU_MIN:g}..{MU_MAX:g}, step {MU_STEP:g}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- isolated-rod curves per rod: {N_ISOLATED_ROD_CURVES}",
        f"- script: `scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods.py`",
        "",
        "## Branch Convention",
        "",
        "The plotted system curves are descendant branches seeded at `mu=0`.",
        "Sorted position is diagnostic metadata only and is not used as branch",
        "identity. Tracking warnings are retained in this report and are not",
        "drawn on the figure.",
        "",
        "## Isolated-Rod Convention",
        "",
        "The auxiliary curves use the same project convention as the beta=15",
        "diagnostic: isolated clamped-supported (CS, clamped-pinned) rods with",
        "roots `tan(alpha)=tanh(alpha)`.",
        "",
        "The thickness-mismatch Lambda normalization is included explicitly:",
        "",
        "- rod 1: `Lambda_rod_1 = alpha_n*sqrt(tau1)/(1-mu)`",
        "- rod 2: `Lambda_rod_2 = alpha_n*sqrt(tau2)/(1+mu)`",
        "",
        "Here `tau1` and `tau2` are the mass-preserving eta factors. For",
        "`eta=0`, these formulas reduce to the usual equal-radius references",
        "`alpha_n/(1-mu)` and `alpha_n/(1+mu)`.",
        "",
        f"- first CS alpha roots used: {', '.join(f'{value:.6g}' for value in alphas[:min(6, len(alphas))])}",
        "",
        "## Thin-Rod Criterion",
        "",
        f"- criterion: `2*r_i/l_i <= {THICKNESS_RATIO_LIMIT:g}` for both rods",
    ]
    if summary.has_violations:
        for segment in summary.segments:
            lines.append(
                f"- WARNING eta={ETA:g}: criterion violated on mu={segment.start_mu:g}..{segment.end_mu:g}, "
                f"rod(s) {rods_label(segment.rods)}"
            )
        lines.append(
            f"- max ratio={summary.max_ratio:.6g} at mu={summary.max_ratio_mu:g}, rod {summary.max_ratio_rod}"
        )
    else:
        lines.append(
            f"- no violations on the grid; max ratio={summary.max_ratio:.6g} "
            f"at mu={summary.max_ratio_mu:g}, rod {summary.max_ratio_rod}"
        )

    lines.extend(
        [
            "",
            "## Tracking Warnings",
            "",
            f"- low_mac={warning_summary['low_mac']}",
            f"- low_margin={warning_summary['low_margin']}",
            f"- unresolved={warning_summary['unresolved']}",
            f"- requires_refined={warning_summary['requires_refined']}",
            f"- frequency/MAC disagreements={warning_summary['frequency_mac_disagreement']}",
        ]
    )
    if warning_rows:
        lines.append("")
        lines.append("Warning rows are diagnostic only; they do not rename descendant branches:")
        for row in warning_rows[:12]:
            lines.append(
                f"- mu={float(row['mu']):g}, desc={int(row['branch_index_from_mu0'])}, "
                f"accepted pos={int(row['mac_sorted_root_index'])}, "
                f"candidate pos={int(row['diagnostic_candidate_sorted_position'])}, "
                f"candidate MAC={float(row['diagnostic_candidate_mac_to_previous']):.6g}, "
                f"status={row['tracking_step_status']}"
            )
        if len(warning_rows) > 12:
            lines.append(f"- plus {len(warning_rows) - 12} additional rows.")

    lines.extend(
        [
            "",
            "## Beta=45 vs Beta=15 Observation",
            "",
            "As a cautious visual proxy, the report compares each plotted",
            "descendant point with the nearest isolated-rod reference curve at the",
            "same `mu`. This is only a diagnostic nearest-grid-distance measure,",
            "not evidence of modal identity.",
            "",
            f"- beta=45 mean nearest distance: {beta45_metrics['mean_nearest_distance']:.6g}",
            f"- beta=15 mean nearest distance: {beta15_metrics['mean_nearest_distance']:.6g}",
            f"- beta=45 median nearest distance: {beta45_metrics['median_nearest_distance']:.6g}",
            f"- beta=15 median nearest distance: {beta15_metrics['median_nearest_distance']:.6g}",
            "",
            f"On this simple metric, beta=45 has a {mean_relation} mean nearest",
            f"distance and a {median_relation} median nearest distance than beta=15.",
            "Visually this should be read only as a diagnostic alignment tendency",
            "relative to the isolated-rod grid, not as a physical conclusion.",
            "",
            "## Output",
            "",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            "",
            "This is diagnostic-only. It does not modify article files, article",
            "figures, `paper_dorofeev_style`, the old determinant, old solvers,",
            "`formulas.py`, or the existing FEM physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    return {
        "warning_rows": len(warning_rows),
        "max_ratio": summary.max_ratio,
        "thickness_violations": summary.has_violations,
        "beta45_mean_distance": beta45_metrics["mean_nearest_distance"],
        "beta15_mean_distance": beta15_metrics["mean_nearest_distance"],
        "beta45_median_distance": beta45_metrics["median_nearest_distance"],
        "beta15_median_distance": beta15_metrics["median_nearest_distance"],
    }


def main() -> dict[str, object]:
    tracked45, tracking_rows45 = compute_tracking(BETA_DEG)
    tracked15, _tracking_rows15 = compute_tracking(COMPARISON_BETA_DEG)
    alphas, rod1_refs, rod2_refs = plot_with_isolated_rods(tracked=tracked45)
    beta45_metrics = reference_grid_distance_metrics(
        tracked=tracked45,
        rod1_refs=rod1_refs,
        rod2_refs=rod2_refs,
    )
    beta15_metrics = reference_grid_distance_metrics(
        tracked=tracked15,
        rod1_refs=rod1_refs,
        rod2_refs=rod2_refs,
    )
    report = write_report(
        tracking_rows=tracking_rows45,
        alphas=alphas,
        rod1_refs=rod1_refs,
        rod2_refs=rod2_refs,
        beta45_metrics=beta45_metrics,
        beta15_metrics=beta15_metrics,
    )
    print(f"saved beta=45 isolated-rod comparison PNG: {OUTPUT_PNG}")
    print(f"saved beta=45 isolated-rod comparison report: {OUTPUT_REPORT}")
    if report["thickness_violations"]:
        print("WARNING: thickness-ratio criterion violations detected")
    else:
        print(f"no thickness-ratio violations; max={report['max_ratio']:.6g}")
    print(f"tracking warning rows: {report['warning_rows']}")
    print(
        "nearest isolated-grid mean distance: "
        f"beta45={report['beta45_mean_distance']:.6g}, beta15={report['beta15_mean_distance']:.6g}"
    )
    print(
        "nearest isolated-grid median distance: "
        f"beta45={report['beta45_median_distance']:.6g}, beta15={report['beta15_median_distance']:.6g}"
    )
    return {"png": OUTPUT_PNG, "report": OUTPUT_REPORT, **report}


if __name__ == "__main__":
    main()
