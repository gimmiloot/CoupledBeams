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
BETA_DEG = 15.0
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
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods.png"
OUTPUT_REPORT = OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_p0p5_with_isolated_rods_report.md"

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
    lines = [
        "# Eta=0.5 Lambda(mu) With Isolated-Rod References",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- mu range: {MU_MIN:g}..{MU_MAX:g}, step {MU_STEP:g}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- isolated-rod curves per rod: {N_ISOLATED_ROD_CURVES}",
        f"- y-axis cap: {Y_MAX:g}",
        "",
        "## Isolated-Rod Convention",
        "",
        "The auxiliary curves use the same convention as the existing Lambda(mu)",
        "reference curves documented in `docs/theory/assumptions.md` and used by",
        "`src/my_project/analytic/FreqMuNet.py` and",
        "`scripts/run/run_lambda_mu_fixed_beta_analytic.py`: isolated",
        "clamped-supported (CS, clamped-pinned) rods with roots",
        "`tan(alpha)=tanh(alpha)`.",
        "",
        "For the equal-radius model this gives `Lambda=alpha/(1-mu)` for rod 1",
        "and `Lambda=alpha/(1+mu)` for rod 2. In the mass-preserving",
        "thickness-mismatch Lambda scaling, the corresponding diagnostic",
        "references are",
        "",
        "- rod 1: `Lambda = alpha*sqrt(tau1)/(1-mu)`",
        "- rod 2: `Lambda = alpha*sqrt(tau2)/(1+mu)`",
        "",
        "Fixed-fixed variants were found in article diagnostics, but they are not",
        "used here because the project Lambda(mu) reference convention for this",
        "plot family is CS.",
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
            "## Output",
            "",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            "",
            "This is diagnostic-only. It does not modify article files, article",
            "figures, the old determinant, old solvers, or the existing FEM",
            "physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    return {
        "warning_rows": len(warning_rows),
        "max_ratio": summary.max_ratio,
        "thickness_violations": summary.has_violations,
        "rod1_ref_max": float(np.nanmax(rod1_refs)),
        "rod2_ref_max": float(np.nanmax(rod2_refs)),
    }


def main() -> dict[str, object]:
    tracked, tracking_rows = compute_tracking()
    alphas, rod1_refs, rod2_refs = plot_with_isolated_rods(tracked=tracked)
    report = write_report(
        tracking_rows=tracking_rows,
        alphas=alphas,
        rod1_refs=rod1_refs,
        rod2_refs=rod2_refs,
    )
    print(f"saved isolated-rod comparison PNG: {OUTPUT_PNG}")
    print(f"saved isolated-rod comparison report: {OUTPUT_REPORT}")
    if report["thickness_violations"]:
        print("WARNING: thickness-ratio criterion violations detected")
    else:
        print(f"no thickness-ratio violations; max={report['max_ratio']:.6g}")
    print(f"tracking warning rows: {report['warning_rows']}")
    print(f"max isolated rod reference Lambda: rod1={report['rod1_ref_max']:.6g}, rod2={report['rod2_ref_max']:.6g}")
    return {"png": OUTPUT_PNG, "report": OUTPUT_REPORT, **report}


if __name__ == "__main__":
    main()
