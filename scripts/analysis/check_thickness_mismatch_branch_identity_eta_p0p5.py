from __future__ import annotations

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

from my_project.analytic.formulas_thickness_mismatch import find_first_n_roots_eta  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import track_mu_branches_shape_mac  # noqa: E402


# =========================
# User-editable parameters
# =========================
BETA_DEG = 15.0
EPSILON = 0.0025
ETA = 0.5

TRACK_START_MU = 0.0
PRE_AUDIT_MU_STEP = 0.005
MU_MIN = 0.65
MU_MAX = 0.90
MU_STEP = 0.002

BRANCHES_TO_AUDIT = (5, 6, 7)
CANDIDATE_MODES = tuple(range(1, 11))
NUM_SORTED_ROOTS = 12
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
MAX_MU_STEP_FOR_CONFIDENCE = PRE_AUDIT_MU_STEP

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_CSV = OUTPUT_DIR / "thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.csv"
OUTPUT_REPORT = OUTPUT_DIR / "thickness_mismatch_branch_identity_eta_p0p5_beta15_refined_report.md"
OUTPUT_PNG = OUTPUT_DIR / "thickness_mismatch_branch_identity_eta_p0p5_beta15_refined.png"


CSV_FIELDNAMES = [
    "mu",
    "descendant_branch",
    "Lambda_descendant",
    "sorted_position",
    "best_mac_previous_step",
    "assigned_from_previous_sorted_position",
    "sorted_position_jump",
    "suspicious_assignment",
    "low_mac",
    "low_margin",
    "requires_refined_check",
    "mac_margin",
    "raw_mac_sorted_position",
    "raw_sorted_position_jump",
    "frequency_nearest_sorted_position",
    "frequency_mac_disagreement",
    "used_frequency_fallback",
    "tracking_step_status",
]


def tracking_mu_values() -> np.ndarray:
    pre_values = np.arange(TRACK_START_MU, MU_MIN, PRE_AUDIT_MU_STEP, dtype=float)
    audit_values = np.arange(MU_MIN, MU_MAX + 0.5 * MU_STEP, MU_STEP, dtype=float)
    values = sorted({round(float(value), 10) for value in pre_values} | {round(float(value), 10) for value in audit_values})
    return np.asarray(values, dtype=float)


def audit_mu_values(mu_values: np.ndarray) -> np.ndarray:
    return np.asarray(
        [
            float(mu)
            for mu in mu_values
            if float(mu) >= MU_MIN - 1e-12 and float(mu) <= MU_MAX + 1e-12
        ],
        dtype=float,
    )


def roots_for(mu: float) -> np.ndarray:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(BETA_DEG)),
        float(mu),
        EPSILON,
        ETA,
        NUM_SORTED_ROOTS,
        Lmax0=ROOT_LMAX0,
        scan_step=ROOT_SCAN_STEP,
    )
    if np.any(~np.isfinite(roots[: max(CANDIDATE_MODES)])):
        raise RuntimeError(f"Missing roots for mu={mu:g}, eta={ETA:g}.")
    return roots


def sorted_roots_by_mu(mu_values: Sequence[float]) -> dict[float, np.ndarray]:
    return {float(mu): roots_for(float(mu)) for mu in mu_values}


def row_is_audited(row: dict[str, float | int | str]) -> bool:
    return (
        int(row["branch_index_from_mu0"]) in BRANCHES_TO_AUDIT
        and float(row["mu"]) >= MU_MIN - 1e-12
        and float(row["mu"]) <= MU_MAX + 1e-12
    )


def csv_row(row: dict[str, float | int | str]) -> dict[str, float | int | str]:
    return {
        "mu": float(row["mu"]),
        "descendant_branch": int(row["branch_index_from_mu0"]),
        "Lambda_descendant": float(row["Lambda_tracked"]),
        "sorted_position": int(row["mac_sorted_root_index"]),
        "best_mac_previous_step": float(row["mac_to_previous"]),
        "assigned_from_previous_sorted_position": int(row["assigned_from_previous_sorted_position"]),
        "sorted_position_jump": int(row["sorted_position_jump"]),
        "suspicious_assignment": str(row["suspicious_assignment"]),
        "low_mac": str(row["low_mac"]),
        "low_margin": str(row["low_margin"]),
        "requires_refined_check": str(row["requires_refined_check"]),
        "mac_margin": float(row["mac_margin"]),
        "raw_mac_sorted_position": int(row["raw_mac_sorted_root_index"]),
        "raw_sorted_position_jump": int(row["raw_sorted_position_jump"]),
        "frequency_nearest_sorted_position": int(row["nearest_sorted_root_index"]),
        "frequency_mac_disagreement": str(row["frequency_mac_disagreement"]),
        "used_frequency_fallback": str(row["used_frequency_fallback"]),
        "tracking_step_status": str(row["tracking_step_status"]),
    }


def write_csv_rows(path: Path, rows: Sequence[dict[str, float | int | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def rows_by_branch(rows: Sequence[dict[str, float | int | str]]) -> dict[int, list[dict[str, float | int | str]]]:
    by_branch = {int(branch): [] for branch in BRANCHES_TO_AUDIT}
    for row in rows:
        by_branch[int(row["descendant_branch"])].append(row)
    for branch_rows in by_branch.values():
        branch_rows.sort(key=lambda item: float(item["mu"]))
    return by_branch


def format_table(rows: Sequence[Sequence[str]]) -> list[str]:
    if not rows:
        return []
    widths = [max(len(row[col]) for row in rows) for col in range(len(rows[0]))]
    lines = []
    for idx, row in enumerate(rows):
        lines.append("| " + " | ".join(value.ljust(widths[col]) for col, value in enumerate(row)) + " |")
        if idx == 0:
            lines.append("| " + " | ".join("-" * widths[col] for col in range(len(row))) + " |")
    return lines


def first_nonself_rows(
    rows: Sequence[dict[str, float | int | str]],
) -> dict[int, dict[str, float | int | str] | None]:
    first_rows: dict[int, dict[str, float | int | str] | None] = {}
    for branch in BRANCHES_TO_AUDIT:
        branch_rows = sorted(
            (row for row in rows if int(row["descendant_branch"]) == int(branch)),
            key=lambda item: float(item["mu"]),
        )
        first_rows[int(branch)] = next(
            (row for row in branch_rows if int(row["sorted_position"]) != int(branch)),
            None,
        )
    return first_rows


def audit_summary(
    rows: Sequence[dict[str, float | int | str]],
    *,
    full_rows: Sequence[dict[str, float | int | str]],
) -> dict[str, object]:
    by_branch = rows_by_branch(rows)
    positions_by_branch = {
        branch: sorted({int(row["sorted_position"]) for row in branch_rows})
        for branch, branch_rows in by_branch.items()
    }
    branch5_positions = positions_by_branch.get(5, [])
    branch5_hits_7 = [row for row in by_branch.get(5, []) if int(row["sorted_position"]) == 7]
    suspicious_rows = [row for row in rows if str(row["suspicious_assignment"]) == "yes"]
    low_mac_rows = [row for row in rows if str(row["low_mac"]) == "yes"]
    low_margin_rows = [row for row in rows if str(row["low_margin"]) == "yes"]
    refined_rows = [row for row in rows if str(row["requires_refined_check"]) == "yes"]
    jump_rows = [row for row in rows if abs(int(row["sorted_position_jump"])) > MAX_SORTED_POSITION_JUMP]
    descendants_stay_at_same_positions = all(
        positions_by_branch.get(branch, []) == [branch]
        for branch in BRANCHES_TO_AUDIT
    )
    return {
        "positions_by_branch": positions_by_branch,
        "branch5_positions": branch5_positions,
        "branch5_hits_7": branch5_hits_7,
        "suspicious_rows": suspicious_rows,
        "low_mac_rows": low_mac_rows,
        "low_margin_rows": low_margin_rows,
        "requires_refined_check_rows": refined_rows,
        "jump_rows": jump_rows,
        "descendants_stay_at_same_positions": descendants_stay_at_same_positions,
        "first_nonself_rows": first_nonself_rows(full_rows),
    }


def plot_audit(
    *,
    output: Path,
    rows: Sequence[dict[str, float | int | str]],
    roots_by_mu: dict[float, np.ndarray],
    audit_mus: np.ndarray,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    by_branch = rows_by_branch(rows)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    fig, axes = plt.subplots(2, 1, figsize=(9.8, 7.0), sharex=True, height_ratios=[2.1, 1.0])
    ax_freq, ax_pos = axes

    for offset, branch in enumerate(BRANCHES_TO_AUDIT):
        branch_rows = by_branch[int(branch)]
        mu = np.array([float(row["mu"]) for row in branch_rows], dtype=float)
        lambdas = np.array([float(row["Lambda_descendant"]) for row in branch_rows], dtype=float)
        positions = np.array([int(row["sorted_position"]) for row in branch_rows], dtype=int)
        color = colors[offset % len(colors)]
        ax_freq.plot(mu, lambdas, color=color, lw=2.0, label=f"descendant branch {branch}")
        ax_pos.step(mu, positions, where="post", color=color, lw=1.8, label=f"branch {branch}")

    for branch in BRANCHES_TO_AUDIT:
        sorted_values = [float(roots_by_mu[float(mu)][int(branch) - 1]) for mu in audit_mus]
        ax_freq.plot(
            audit_mus,
            sorted_values,
            color="0.35",
            lw=0.9,
            ls=":",
            alpha=0.65,
            label=f"sorted mode {branch}" if branch == BRANCHES_TO_AUDIT[0] else "_nolegend_",
        )
        ax_pos.axhline(float(branch), color="0.80", lw=0.8, ls=":")

    ax_freq.set_ylabel(r"$\Lambda$")
    ax_freq.set_title(
        rf"Descendant-branch identity audit, $\eta={ETA:g}$, "
        rf"$\beta={BETA_DEG:g}^\circ$, $\epsilon={EPSILON:g}$"
    )
    ax_freq.grid(True, color="0.88", linewidth=0.6)
    ax_freq.legend(loc="best", fontsize=8, frameon=False)

    ax_pos.set_xlabel(r"$\mu$")
    ax_pos.set_ylabel("current sorted position")
    ax_pos.set_yticks(range(4, 9))
    ax_pos.set_ylim(3.75, 8.25)
    ax_pos.grid(True, color="0.88", linewidth=0.6)
    ax_pos.legend(loc="best", fontsize=8, frameon=False)

    fig.text(
        0.5,
        0.01,
        "Descendant branches are seeded at mu=0; sorted position is diagnostic metadata, not branch identity.",
        ha="center",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.0, 0.03, 1.0, 1.0))
    fig.savefig(output, dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_report(
    *,
    path: Path,
    rows: Sequence[dict[str, float | int | str]],
    summary: dict[str, object],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    positions_by_branch = summary["positions_by_branch"]
    first_nonself = summary["first_nonself_rows"]
    branch5_hits_7 = summary["branch5_hits_7"]
    descendants_stay = bool(summary["descendants_stay_at_same_positions"])
    suspicious_rows = summary["suspicious_rows"]
    low_mac_rows = summary["low_mac_rows"]
    low_margin_rows = summary["low_margin_rows"]
    refined_rows = summary["requires_refined_check_rows"]
    jump_rows = summary["jump_rows"]

    position_table = [["descendant branch", "sorted positions in refined window"]]
    for branch in BRANCHES_TO_AUDIT:
        values = ", ".join(str(value) for value in positions_by_branch[branch])
        position_table.append([str(branch), values])

    warning_table = [["mu", "branch", "sorted position", "jump", "MAC", "flags"]]
    for row in list(suspicious_rows)[:20]:
        flags = ",".join(
            name
            for name in ("low_mac", "low_margin", "requires_refined_check")
            if str(row[name]) == "yes"
        )
        warning_table.append(
            [
                f"{float(row['mu']):.3f}",
                str(int(row["descendant_branch"])),
                str(int(row["sorted_position"])),
                str(int(row["sorted_position_jump"])),
                f"{float(row['best_mac_previous_step']):.6f}",
                flags or "suspicious",
            ]
        )

    context_table = [["branch", "first non-seed sorted position", "mu", "MAC", "status"]]
    for branch in BRANCHES_TO_AUDIT:
        row = first_nonself[int(branch)]
        if row is None:
            context_table.append([str(branch), "none", "-", "-", "-"])
        else:
            context_table.append(
                [
                    str(branch),
                    str(int(row["sorted_position"])),
                    f"{float(row['mu']):.3f}",
                    f"{float(row['best_mac_previous_step']):.6f}",
                    str(row["tracking_step_status"]),
                ]
            )

    lines = [
        "# Thickness-Mismatch Branch-Identity Refined Audit",
        "",
        "## Convention",
        "",
        "`branch k` means the descendant of the k-th mode shape selected at the",
        "tracking start point. In this audit the descendants are seeded at",
        "`mu=0`; the refined output window is `mu=0.65..0.90`. The sorted root",
        "index at a later mu value is diagnostic metadata only and is not used to",
        "rename the branch.",
        "",
        "A sorted-position jump larger than one in one step is a suspicious",
        "assignment. Low MAC or low MAC margin also require local review.",
        "",
        "## Parameters",
        "",
        f"- beta: {BETA_DEG:g} deg",
        f"- epsilon: {EPSILON:g}",
        f"- eta: {ETA:g}",
        f"- tracking start mu: {TRACK_START_MU:g}",
        f"- pre-audit step: {PRE_AUDIT_MU_STEP:g}",
        f"- refined window: {MU_MIN:g} .. {MU_MAX:g}",
        f"- refined step: {MU_STEP:g}",
        f"- audited descendants: {', '.join(str(value) for value in BRANCHES_TO_AUDIT)}",
        f"- tracked candidate branches: 1..{max(CANDIDATE_MODES)}",
        f"- sorted roots scanned: first {NUM_SORTED_ROOTS}",
        f"- MAC threshold: {MAC_WARNING_THRESHOLD:g}",
        f"- MAC margin threshold: {MAC_MARGIN_WARNING_THRESHOLD:g}",
        "",
        "## Refined Result",
        "",
    ]
    lines.extend(format_table(position_table))
    lines.extend(
        [
            "",
            f"- descendant branch 5 sorted positions: {', '.join(str(value) for value in positions_by_branch[5])}",
            f"- descendant branch 5 reaches sorted position 7: {'yes' if branch5_hits_7 else 'no'}",
            f"- descendants 5, 6, and 7 stay at sorted positions 5, 6, and 7 after mu about 0.7: {'yes' if descendants_stay else 'no'}",
            "",
        ]
    )

    if branch5_hits_7:
        lines.extend(
            [
                "Conclusion: the refined audit found branch 5 at sorted position 7.",
                "This would need additional shape review before a physical",
                "interpretation is made.",
                "",
            ]
        )
    elif positions_by_branch[5] != [5]:
        first_branch5 = first_nonself[5]
        lines.extend(
            [
                "Conclusion: the refined audit does not confirm a `5 -> 7` sorted",
                "position transition for descendant branch 5. In this raw-MAC run,",
                "branch 5 enters the high-mu audit window at a different sorted",
                "position, but the first such departure from its seed sorted",
                f"position occurs at `mu={float(first_branch5['mu']):g}` with",
                f"`MAC={float(first_branch5['best_mac_previous_step']):.6g}` and",
                f"status `{first_branch5['tracking_step_status']}`.",
                "That pre-window low-MAC assignment is diagnostic only and should",
                "not be interpreted physically without a dedicated refined audit",
                "around that earlier interaction.",
                "",
            ]
        )
    else:
        lines.extend(
            [
                "Conclusion: the refined audit does not confirm a `5 -> 7` sorted",
                "position transition for descendant branch 5. In this refined run,",
                "the earlier `5 -> 7` interpretation is treated as an artifact of",
                "coarser tracking and/or ambiguous shape assignment, not as branch",
                "renumbering.",
                "",
            ]
        )

    lines.extend(["## Full-Path Context", ""])
    lines.append(
        "The table below is computed on the tracking path from `mu=0` to the "
        "audit window. It is included only to show whether the high-mu result "
        "depends on an earlier suspicious assignment."
    )
    lines.append("")
    lines.extend(format_table(context_table))
    lines.extend(
        [
            "",
            "Rows with low MAC or low MAC margin are not accepted as canonical",
            "branch-identity evidence until a local refined audit resolves them.",
            "",
        ]
    )

    lines.extend(
        [
            "## Warnings",
            "",
            f"- suspicious audited rows: {len(suspicious_rows)}",
            f"- low-MAC audited rows: {len(low_mac_rows)}",
            f"- low-margin audited rows: {len(low_margin_rows)}",
            f"- sorted-position jumps larger than {MAX_SORTED_POSITION_JUMP}: {len(jump_rows)}",
            f"- rows requiring refined check: {len(refined_rows)}",
            "",
        ]
    )
    if suspicious_rows:
        lines.extend(format_table(warning_table))
        if len(suspicious_rows) > 20:
            lines.append(f"Only the first 20 of {len(suspicious_rows)} suspicious rows are shown.")
    else:
        lines.append("No suspicious assignments were found for the audited descendants on the refined grid.")

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            f"- CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
            f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
            "",
            "This diagnostic does not change the large-eta plot, article files,",
            "article figures, old determinant, old solvers, or FEM physical model.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> dict[str, object]:
    mu_values = tracking_mu_values()
    audit_mus = audit_mu_values(mu_values)
    roots = sorted_roots_by_mu(mu_values)
    result = track_mu_branches_shape_mac(
        beta_rad=float(np.deg2rad(BETA_DEG)),
        epsilon=EPSILON,
        eta=ETA,
        mu_values=mu_values,
        roots_by_mu=roots,
        num_tracked_branches=max(CANDIDATE_MODES),
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
        max_mu_step_for_confidence=MAX_MU_STEP_FOR_CONFIDENCE,
    )
    full_rows = [
        csv_row(row)
        for row in result.rows
        if int(row["branch_index_from_mu0"]) in BRANCHES_TO_AUDIT
    ]
    rows = [csv_row(row) for row in result.rows if row_is_audited(row)]
    rows.sort(key=lambda row: (float(row["mu"]), int(row["descendant_branch"])))
    full_rows.sort(key=lambda row: (float(row["mu"]), int(row["descendant_branch"])))
    summary = audit_summary(rows, full_rows=full_rows)

    write_csv_rows(OUTPUT_CSV, rows)
    plot_audit(output=OUTPUT_PNG, rows=rows, roots_by_mu=roots, audit_mus=audit_mus)
    write_report(path=OUTPUT_REPORT, rows=rows, summary=summary)

    positions_by_branch = summary["positions_by_branch"]
    descendants_stay = bool(summary["descendants_stay_at_same_positions"])
    branch5_hits_7 = summary["branch5_hits_7"]
    suspicious_rows = summary["suspicious_rows"]
    first_nonself = summary["first_nonself_rows"]

    print(f"saved refined branch-identity CSV: {OUTPUT_CSV}")
    print(f"saved refined branch-identity PNG: {OUTPUT_PNG}")
    print(f"saved refined branch-identity report: {OUTPUT_REPORT}")
    for branch in BRANCHES_TO_AUDIT:
        positions = ", ".join(str(value) for value in positions_by_branch[branch])
        print(f"descendant branch {branch}: sorted positions in refined window = {positions}")
    print(f"descendant branch 5 reaches sorted position 7: {'yes' if branch5_hits_7 else 'no'}")
    print(f"descendants 5, 6, 7 stay at sorted positions 5, 6, 7: {'yes' if descendants_stay else 'no'}")
    first_branch5 = first_nonself[5]
    if first_branch5 is not None:
        print(
            "descendant branch 5 first non-seed sorted position: "
            f"mu={float(first_branch5['mu']):g}, "
            f"sorted position={int(first_branch5['sorted_position'])}, "
            f"MAC={float(first_branch5['best_mac_previous_step']):.6g}, "
            f"status={first_branch5['tracking_step_status']}"
        )
    print(f"suspicious audited rows: {len(suspicious_rows)}")
    print("large-eta plot was not regenerated by this audit script")

    return {
        "csv": OUTPUT_CSV,
        "png": OUTPUT_PNG,
        "report": OUTPUT_REPORT,
        "positions_by_branch": positions_by_branch,
        "descendants_stay": descendants_stay,
        "branch5_hits_7": bool(branch5_hits_7),
        "suspicious_rows": len(suspicious_rows),
    }


if __name__ == "__main__":
    main()
