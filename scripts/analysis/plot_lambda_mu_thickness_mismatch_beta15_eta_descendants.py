from __future__ import annotations

from collections import defaultdict
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
ETA_VALUES = (-0.5, 0.0, 0.5)

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.005

N_DESCENDANTS_PLOT = 6
N_DESCENDANTS_TRACK = 8
N_SORTED_SCAN = 12
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
THICKNESS_RATIO_LIMIT = 0.1

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants.png"
OUTPUT_REPORT = (
    OUTPUT_DIR / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants_report.md"
)

MU_VALUES = mu_grid(MU_MIN, MU_MAX, MU_STEP)


def eta_title(eta: float) -> str:
    if eta > 0.0:
        meaning = "rod 2 thicker"
    elif eta < 0.0:
        meaning = "rod 1 thicker"
    else:
        meaning = "equal radii"
    return rf"$\eta={eta:g}$" + "\n" + meaning


def compute_case(eta: float) -> dict[str, object]:
    roots = roots_by_mu_eta(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=float(eta),
        mu_values=MU_VALUES,
        n_roots=N_SORTED_SCAN,
        root_lmax0=ROOT_LMAX0,
        root_scan_step=ROOT_SCAN_STEP,
    )
    tracking = track_descendants_from_mu0(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=float(eta),
        mu_values=MU_VALUES,
        roots_by_mu=roots,
        num_descendants=N_DESCENDANTS_TRACK,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
    )
    summary = thickness_ratio_summary(
        epsilon=EPSILON,
        eta=float(eta),
        mu_values=MU_VALUES,
        limit=THICKNESS_RATIO_LIMIT,
    )
    return {"roots": roots, "tracking": tracking, "summary": summary}


def warning_points(rows: list[dict[str, float | int | str]], eta: float) -> tuple[np.ndarray, np.ndarray]:
    xs: list[float] = []
    ys: list[float] = []
    for row in rows:
        if not np.isclose(float(row["eta"]), float(eta), rtol=0.0, atol=1e-12):
            continue
        if int(row["branch_index_from_mu0"]) > N_DESCENDANTS_PLOT:
            continue
        if str(row.get("requires_refined_check", "no")) != "yes":
            continue
        xs.append(float(row["mu"]))
        ys.append(float(row["Lambda_tracked"]))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def plot_cases(cases: dict[float, dict[str, object]]) -> None:
    OUTPUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(ETA_VALUES), figsize=(11.4, 4.4), sharey=True)
    axes = np.atleast_1d(axes)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, eta in zip(axes, ETA_VALUES, strict=True):
        case = cases[float(eta)]
        tracking = case["tracking"]
        tracked = tracking.tracked
        valid = valid_mask(
            epsilon=EPSILON,
            eta=float(eta),
            mu_values=MU_VALUES,
            limit=THICKNESS_RATIO_LIMIT,
        )
        for branch_idx in range(N_DESCENDANTS_PLOT):
            plot_with_validity_split(
                ax,
                MU_VALUES,
                tracked[branch_idx],
                valid,
                color=colors[branch_idx % len(colors)],
                linewidth=1.8,
                label=f"desc {branch_idx + 1}" if ax is axes[0] else None,
                zorder=3,
            )

        ax.set_title(eta_title(float(eta)), fontsize=10)
        ax.set_xlabel(r"$\mu$")
        ax.grid(True, color="0.88", linewidth=0.6)

    axes[0].set_ylabel(r"$\Lambda$")
    branch_handles = [
        Line2D([0], [0], color=colors[idx % len(colors)], lw=1.8, label=f"desc {idx + 1}")
        for idx in range(N_DESCENDANTS_PLOT)
    ]
    style_handles = [
        Line2D([0], [0], color="0.2", lw=1.8, ls="-", label=r"solid: $2r_i/l_i \leq 0.1$"),
        Line2D([0], [0], color="0.2", lw=1.8, ls="--", label="dashed: criterion violated"),
    ]
    fig.legend(
        handles=branch_handles + style_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=4,
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(
        rf"Thickness-mismatch descendant branches, $\beta={BETA_DEG:g}^\circ$, $\epsilon={EPSILON:g}$",
        fontsize=12,
    )
    fig.tight_layout(rect=(0.0, 0.17, 1.0, 0.92))
    fig.savefig(OUTPUT_PNG, dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_report(cases: dict[float, dict[str, object]]) -> dict[str, object]:
    all_rows: list[dict[str, float | int | str]] = []
    lines = [
        "# Thickness-Mismatch Lambda(mu) Descendant Diagnostic",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta values: {', '.join(f'{float(value):g}' for value in ETA_VALUES)}",
        f"- mu range: {MU_MIN:g} .. {MU_MAX:g}",
        f"- mu step: {MU_STEP:g}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- tracked descendants: first {N_DESCENDANTS_TRACK}",
        f"- sorted roots scanned at each mu: {N_SORTED_SCAN}",
        "",
        "Branches are descendant branches seeded at `mu=0`; sorted positions are",
        "diagnostic metadata only. Tracking uses adjacent-step analytic shape MAC.",
        "Low-MAC, low-margin, or large-jump candidates are warnings and do not",
        "automatically rename a branch.",
        "Tracking warnings are retained in this report but are not drawn on the",
        "presentation-style PNG; the figure shows only descendant branches and",
        "the solid/dashed thin-rod applicability split.",
        "",
        f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
        "",
        "## Thin-Rod Applicability",
        "",
        f"- criterion: `2*r_i/l_i <= {THICKNESS_RATIO_LIMIT:g}` for both rods",
    ]

    for eta in ETA_VALUES:
        case = cases[float(eta)]
        summary = case["summary"]
        if summary.has_violations:
            segment_text = "; ".join(
                f"mu={segment.start_mu:g}..{segment.end_mu:g}, rod(s)={rods_label(segment.rods)}"
                for segment in summary.segments
            )
            lines.append(
                f"- WARNING eta={float(eta):g}: violations: {segment_text}; "
                f"max ratio={summary.max_ratio:.6g} at mu={summary.max_ratio_mu:g}, "
                f"rod {summary.max_ratio_rod}"
            )
        else:
            lines.append(
                f"- eta={float(eta):g}: no violations; max ratio={summary.max_ratio:.6g} "
                f"at mu={summary.max_ratio_mu:g}, rod {summary.max_ratio_rod}"
            )
        all_rows.extend(case["tracking"].rows)

    warnings = tracking_warning_rows(all_rows)
    by_eta: dict[float, list[dict[str, float | int | str]]] = defaultdict(list)
    for row in warnings:
        by_eta[float(row["eta"])].append(row)
    lines.extend(["", "## Tracking Warnings", ""])
    if warnings:
        lines.append(f"- total warning rows: {len(warnings)}")
        for eta in ETA_VALUES:
            summary = tracking_warning_summary(cases[float(eta)]["tracking"].rows)
            eta_rows = by_eta.get(float(eta), [])
            lines.append(
                f"- eta={float(eta):g}: warning rows={len(eta_rows)}, "
                f"low_mac={summary['low_mac']}, low_margin={summary['low_margin']}, "
                f"unresolved={summary['unresolved']}, requires_refined={summary['requires_refined']}, "
                f"frequency/MAC disagreements={summary['frequency_mac_disagreement']}"
            )
            for row in eta_rows[:8]:
                lines.append(
                    f"  - mu={float(row['mu']):g}, desc={int(row['branch_index_from_mu0'])}, "
                    f"accepted position={int(row['mac_sorted_root_index'])}, "
                    f"candidate position={int(row['diagnostic_candidate_sorted_position'])}, "
                    f"candidate MAC={float(row['diagnostic_candidate_mac_to_previous']):.6g}, "
                    f"status={row['tracking_step_status']}"
                )
            if len(eta_rows) > 8:
                lines.append(f"  - plus {len(eta_rows) - 8} additional rows.")
    else:
        lines.append("- no tracking warnings on plotted/tracked descendants.")

    lines.extend(
        [
            "",
            "This diagnostic writes only results files and does not modify article",
            "files, article figures, the baseline determinant, old solvers, or the",
            "FEM physical model.",
            "",
        ]
    )
    OUTPUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")
    return {"warning_rows": len(warnings), "rows": all_rows}


def main() -> dict[str, object]:
    cases = {float(eta): compute_case(float(eta)) for eta in ETA_VALUES}
    plot_cases(cases)
    report = write_report(cases)

    print(f"saved descendant Lambda(mu) PNG: {OUTPUT_PNG}")
    print(f"saved descendant Lambda(mu) report: {OUTPUT_REPORT}")
    for eta in ETA_VALUES:
        summary = cases[float(eta)]["summary"]
        if summary.has_violations:
            for segment in summary.segments:
                print(
                    f"WARNING eta={float(eta):g}: 2*r_i/l_i violated for rod(s) "
                    f"{rods_label(segment.rods)} on mu={segment.start_mu:g}..{segment.end_mu:g}; "
                    f"max={summary.max_ratio:.6g} at mu={summary.max_ratio_mu:g}, rod={summary.max_ratio_rod}"
                )
        else:
            print(
                f"eta={float(eta):g}: no thin-rod violations; "
                f"max 2*r_i/l_i={summary.max_ratio:.6g} at mu={summary.max_ratio_mu:g}, "
                f"rod={summary.max_ratio_rod}"
            )
        warning_summary = tracking_warning_summary(cases[float(eta)]["tracking"].rows)
        print(
            f"eta={float(eta):g}: tracking warnings "
            f"low_mac={warning_summary['low_mac']}, low_margin={warning_summary['low_margin']}, "
            f"unresolved={warning_summary['unresolved']}, "
            f"requires_refined={warning_summary['requires_refined']}"
        )

    return {"png": OUTPUT_PNG, "report": OUTPUT_REPORT, **report}


if __name__ == "__main__":
    main()
