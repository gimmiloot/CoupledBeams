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

from my_project.analytic.FreqMuNet import roots_clamped_supported  # noqa: E402
from my_project.analytic.solvers import fixed_fixed_lambdas  # noqa: E402
from my_project.analytic.formulas_thickness_mismatch import thickness_mismatch_factors  # noqa: E402
from scripts.lib.thickness_mismatch_diagnostic_helpers import (  # noqa: E402
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
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed.png"
OUTPUT_REPORT = (
    OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta45_eps0p0025_eta_p0p5_with_isolated_rods_fixed_fixed_report.md"
)

MU_VALUES = mu_grid(MU_MIN, MU_MAX, MU_STEP)


def compute_tracking() -> tuple[np.ndarray, list[dict[str, float | int | str]]]:
    roots = roots_by_mu_eta(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=MU_VALUES,
        n_roots=N_SORTED_SCAN,
        root_lmax0=ROOT_LMAX0,
        root_scan_step=ROOT_SCAN_STEP,
    )
    result = track_descendants_from_mu0(
        beta_rad=float(np.deg2rad(BETA_DEG)),
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


def isolated_fixed_fixed_references() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Clamped-clamped isolated-rod references with eta-model Lambda scaling."""
    alphas = fixed_fixed_lambdas(N_ISOLATED_ROD_CURVES)
    rod1 = np.full((N_ISOLATED_ROD_CURVES, len(MU_VALUES)), np.nan, dtype=float)
    rod2 = np.full_like(rod1, np.nan)
    for col, mu in enumerate(MU_VALUES):
        factors = thickness_mismatch_factors(float(mu), ETA)
        rod1[:, col] = alphas * np.sqrt(factors.tau1) / (1.0 - float(mu))
        rod2[:, col] = alphas * np.sqrt(factors.tau2) / (1.0 + float(mu))
    return alphas, rod1, rod2


def isolated_clamped_pinned_references() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Previous clamped-pinned references, used only for a diagnostic comparison metric."""
    alphas = roots_clamped_supported(N_ISOLATED_ROD_CURVES)
    rod1 = np.full((N_ISOLATED_ROD_CURVES, len(MU_VALUES)), np.nan, dtype=float)
    rod2 = np.full_like(rod1, np.nan)
    for col, mu in enumerate(MU_VALUES):
        factors = thickness_mismatch_factors(float(mu), ETA)
        rod1[:, col] = alphas * np.sqrt(factors.tau1) / (1.0 - float(mu))
        rod2[:, col] = alphas * np.sqrt(factors.tau2) / (1.0 + float(mu))
    return alphas, rod1, rod2


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


def plot_with_fixed_fixed_rods(*, tracked: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    valid = valid_mask(epsilon=EPSILON, eta=ETA, mu_values=MU_VALUES, limit=THICKNESS_RATIO_LIMIT)
    alphas, rod1_refs, rod2_refs = isolated_fixed_fixed_references()

    fig, ax = plt.subplots(figsize=(9.8, 5.8))
    for idx in range(N_ISOLATED_ROD_CURVES):
        ax.plot(
            MU_VALUES,
            rod1_refs[idx],
            color="0.32",
            lw=1.05,
            ls=(0, (6.0, 2.0, 1.4, 2.0)),
            alpha=0.44,
            zorder=1,
        )
        ax.plot(
            MU_VALUES,
            rod2_refs[idx],
            color="0.55",
            lw=1.05,
            ls=":",
            alpha=0.46,
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
        Line2D([0], [0], color="0.32", lw=1.1, ls=(0, (6.0, 2.0, 1.4, 2.0)), label="isolated rod 1, fixed-fixed"),
        Line2D([0], [0], color="0.55", lw=1.1, ls=":", label="isolated rod 2, fixed-fixed"),
        Line2D([0], [0], color="0.2", lw=1.8, ls="-", label=r"solid: $2r_i/l_i \leq 0.1$"),
        Line2D([0], [0], color="0.2", lw=1.8, ls="--", label="dashed system segment: criterion violated"),
    ]
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\Lambda$")
    ax.set_title(
        rf"Thickness-mismatch descendants with fixed-fixed isolated rods, "
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
    alphas_ff: np.ndarray,
    metrics_ff: dict[str, float],
    metrics_cp: dict[str, float],
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
    mean_relation = "smaller" if metrics_ff["mean_nearest_distance"] < metrics_cp["mean_nearest_distance"] else "larger"
    median_relation = "smaller" if metrics_ff["median_nearest_distance"] < metrics_cp["median_nearest_distance"] else "larger"

    lines = [
        "# Eta=0.5 Lambda(mu) With Fixed-Fixed Isolated-Rod References, beta=45",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu range: {MU_MIN:g}..{MU_MAX:g}, step {MU_STEP:g}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- isolated fixed-fixed curves per rod: {N_ISOLATED_ROD_CURVES}",
        f"- script: `scripts/analysis/plot_lambda_mu_thickness_mismatch_beta45_eta_p0p5_with_isolated_rods_fixed_fixed.py`",
        "",
        "## Branch Convention",
        "",
        "The plotted system curves are descendant branches seeded at `mu=0`.",
        "Sorted position is diagnostic metadata only and is not used as branch",
        "identity. Tracking warning markers are not drawn on the figure.",
        "",
        "## Fixed-Fixed Isolated-Rod Convention",
        "",
        "This plot changes only the isolated-rod boundary condition relative to",
        "the previous beta=45 clamped-pinned reference plot. The system branches,",
        "eta model, tracking logic, and thin-rod applicability split are kept the",
        "same.",
        "",
        "The isolated rods are clamped-clamped / fixed-fixed Euler-Bernoulli",
        "beams. Their alpha roots solve",
        "",
        "`cosh(alpha_n)*cos(alpha_n) = 1`.",
        "",
        "The thickness-mismatch Lambda normalization is kept:",
        "",
        "- rod 1: `Lambda_rod_1 = alpha_n*sqrt(tau1)/(1-mu)`",
        "- rod 2: `Lambda_rod_2 = alpha_n*sqrt(tau2)/(1+mu)`",
        "",
        "Here `tau1` and `tau2` are the mass-preserving eta factors. The only",
        "change from the clamped-pinned plot is the replacement of the alpha",
        "roots: clamped-pinned roots of `tan(alpha)=tanh(alpha)` are replaced by",
        "fixed-fixed roots of `cosh(alpha) cos(alpha)=1`.",
        "",
        f"- fixed-fixed alpha roots used: {', '.join(f'{value:.6g}' for value in alphas_ff[:min(6, len(alphas_ff))])}",
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
            "## Fixed-Fixed vs Clamped-Pinned Observation",
            "",
            "As a cautious diagnostic proxy, the report compares each plotted",
            "descendant point with the nearest isolated-rod reference curve at the",
            "same `mu`. This nearest-grid-distance measure is not modal identity",
            "evidence.",
            "",
            f"- fixed-fixed mean nearest distance: {metrics_ff['mean_nearest_distance']:.6g}",
            f"- clamped-pinned mean nearest distance: {metrics_cp['mean_nearest_distance']:.6g}",
            f"- fixed-fixed median nearest distance: {metrics_ff['median_nearest_distance']:.6g}",
            f"- clamped-pinned median nearest distance: {metrics_cp['median_nearest_distance']:.6g}",
            "",
            f"On this simple metric, the fixed-fixed grid has a {mean_relation}",
            f"mean nearest distance and a {median_relation} median nearest",
            "distance than the previous clamped-pinned grid. This is only a",
            "diagnostic visual-alignment observation.",
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
        "fixed_fixed_mean_distance": metrics_ff["mean_nearest_distance"],
        "clamped_pinned_mean_distance": metrics_cp["mean_nearest_distance"],
        "fixed_fixed_median_distance": metrics_ff["median_nearest_distance"],
        "clamped_pinned_median_distance": metrics_cp["median_nearest_distance"],
    }


def main() -> dict[str, object]:
    tracked, tracking_rows = compute_tracking()
    alphas_ff, rod1_ff, rod2_ff = plot_with_fixed_fixed_rods(tracked=tracked)
    _alphas_cp, rod1_cp, rod2_cp = isolated_clamped_pinned_references()
    metrics_ff = reference_grid_distance_metrics(tracked=tracked, rod1_refs=rod1_ff, rod2_refs=rod2_ff)
    metrics_cp = reference_grid_distance_metrics(tracked=tracked, rod1_refs=rod1_cp, rod2_refs=rod2_cp)
    report = write_report(
        tracking_rows=tracking_rows,
        alphas_ff=alphas_ff,
        metrics_ff=metrics_ff,
        metrics_cp=metrics_cp,
    )
    print(f"saved beta=45 fixed-fixed isolated-rod PNG: {OUTPUT_PNG}")
    print(f"saved beta=45 fixed-fixed isolated-rod report: {OUTPUT_REPORT}")
    if report["thickness_violations"]:
        print("WARNING: thickness-ratio criterion violations detected")
    else:
        print(f"no thickness-ratio violations; max={report['max_ratio']:.6g}")
    print(f"tracking warning rows: {report['warning_rows']}")
    print(
        "nearest isolated-grid mean distance: "
        f"fixed_fixed={report['fixed_fixed_mean_distance']:.6g}, "
        f"clamped_pinned={report['clamped_pinned_mean_distance']:.6g}"
    )
    print(
        "nearest isolated-grid median distance: "
        f"fixed_fixed={report['fixed_fixed_median_distance']:.6g}, "
        f"clamped_pinned={report['clamped_pinned_median_distance']:.6g}"
    )
    return {"png": OUTPUT_PNG, "report": OUTPUT_REPORT, **report}


if __name__ == "__main__":
    main()
