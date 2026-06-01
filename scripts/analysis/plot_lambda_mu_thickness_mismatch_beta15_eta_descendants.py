from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
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


# ============================
# Default diagnostic parameters
# ============================
DEFAULT_BETA_DEG_VALUES = (15.0,)
DEFAULT_EPSILON = 0.0025
DEFAULT_ETA_VALUES = (-0.5, 0.0, 0.5)

MU_MIN = 0.0
MU_MAX = 0.9
MU_STEP = 0.005

N_DESCENDANTS_PLOT = 6
N_DESCENDANTS_TRACK = 8
N_SORTED_SCAN = 12
N_SORTED_SCAN_FALLBACK = 10
ROOT_SCAN_STEP = 0.01
ROOT_LMAX0 = 35.0
NUM_SHAPE_SAMPLES = 401

MAC_WARNING_THRESHOLD = 0.9
MAC_MARGIN_WARNING_THRESHOLD = 0.05
MAX_SORTED_POSITION_JUMP = 1
THICKNESS_RATIO_LIMIT = 0.1
CLOSE_APPROACH_ABS_THRESHOLD = 0.05
CLOSE_APPROACH_REL_THRESHOLD = 0.005

DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)


def eta_title(eta: float) -> str:
    if eta > 0.0:
        meaning = "rod 2 thicker"
    elif eta < 0.0:
        meaning = "rod 1 thicker"
    else:
        meaning = "equal radii"
    return rf"$\eta = {eta:g}$" + "\n" + meaning


def number_token(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace(".", "p")


def signed_number_token(value: float) -> str:
    if float(value) > 0.0:
        return "p" + number_token(float(value))
    return number_token(float(value))


def output_paths(
    *,
    output_dir: Path,
    beta_deg: float,
    epsilon: float,
    eta_values: tuple[float, ...],
) -> tuple[Path, Path, Path, Path]:
    eta_token = "_".join(signed_number_token(float(eta)) for eta in eta_values)
    base_stem = (
        f"thickness_mismatch_lambda_mu_beta{number_token(float(beta_deg))}"
        f"_eps{number_token(float(epsilon))}_eta_{eta_token}"
    )
    descendant_stem = f"{base_stem}_descendants"
    return (
        output_dir / f"{descendant_stem}.png",
        output_dir / f"{descendant_stem}_report.md",
        output_dir / f"{descendant_stem}.csv",
        output_dir / f"{base_stem}_mac_tracking_warnings.csv",
    )


def compute_roots_for_case(
    *,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
) -> tuple[dict[float, np.ndarray], int, str]:
    last_error: RuntimeError | None = None
    for n_roots in (N_SORTED_SCAN, N_SORTED_SCAN_FALLBACK):
        if n_roots < N_DESCENDANTS_TRACK:
            continue
        try:
            roots = roots_by_mu_eta(
                beta_rad=float(np.deg2rad(beta_deg)),
                epsilon=float(epsilon),
                eta=float(eta),
                mu_values=mu_values,
                n_roots=int(n_roots),
                root_lmax0=ROOT_LMAX0,
                root_scan_step=ROOT_SCAN_STEP,
            )
        except RuntimeError as exc:
            if "Missing roots" not in str(exc):
                raise
            last_error = exc
            continue
        if int(n_roots) == N_SORTED_SCAN:
            return roots, int(n_roots), "primary scan"
        return (
            roots,
            int(n_roots),
            f"fallback after primary {N_SORTED_SCAN}-root scan hit missing high diagnostic roots: {last_error}",
        )
    if last_error is not None:
        raise last_error
    raise RuntimeError("No usable sorted-root scan size is available.")


def compute_case(
    *,
    beta_deg: float,
    epsilon: float,
    eta: float,
    mu_values: np.ndarray,
) -> dict[str, object]:
    roots, n_sorted_scan, root_scan_note = compute_roots_for_case(
        beta_deg=float(beta_deg),
        epsilon=float(epsilon),
        eta=float(eta),
        mu_values=mu_values,
    )
    tracking = track_descendants_from_mu0(
        beta_rad=float(np.deg2rad(beta_deg)),
        epsilon=float(epsilon),
        eta=float(eta),
        mu_values=mu_values,
        roots_by_mu=roots,
        num_descendants=N_DESCENDANTS_TRACK,
        num_shape_samples=NUM_SHAPE_SAMPLES,
        mac_warning_threshold=MAC_WARNING_THRESHOLD,
        mac_margin_warning_threshold=MAC_MARGIN_WARNING_THRESHOLD,
        max_sorted_position_jump=MAX_SORTED_POSITION_JUMP,
    )
    summary = thickness_ratio_summary(
        epsilon=float(epsilon),
        eta=float(eta),
        mu_values=mu_values,
        limit=THICKNESS_RATIO_LIMIT,
    )
    return {
        "roots": roots,
        "tracking": tracking,
        "summary": summary,
        "n_sorted_scan": n_sorted_scan,
        "root_scan_note": root_scan_note,
    }


def close_approach_rows(
    *,
    eta: float,
    mu_values: np.ndarray,
    tracked: np.ndarray,
) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    plotted = np.asarray(tracked[:N_DESCENDANTS_PLOT], dtype=float)
    for left in range(N_DESCENDANTS_PLOT):
        for right in range(left + 1, N_DESCENDANTS_PLOT):
            diff = plotted[left] - plotted[right]
            gaps = np.abs(diff)
            min_idx = int(np.nanargmin(gaps))
            mean_lambda = 0.5 * (abs(float(plotted[left, min_idx])) + abs(float(plotted[right, min_idx])))
            rel_gap = float(gaps[min_idx]) / mean_lambda if mean_lambda > 0.0 else np.inf
            sign_values = np.sign(diff[np.isfinite(diff)])
            sign_changes = bool(np.any(sign_values[:-1] * sign_values[1:] < 0.0)) if len(sign_values) > 1 else False
            close_flag = float(gaps[min_idx]) <= CLOSE_APPROACH_ABS_THRESHOLD or rel_gap <= CLOSE_APPROACH_REL_THRESHOLD
            rows.append(
                {
                    "eta": float(eta),
                    "desc_left": int(left + 1),
                    "desc_right": int(right + 1),
                    "mu": float(mu_values[min_idx]),
                    "gap": float(gaps[min_idx]),
                    "relative_gap": float(rel_gap),
                    "sign_change_on_grid": "yes" if sign_changes else "no",
                    "close_approach_flag": "yes" if close_flag else "no",
                }
            )
    rows.sort(key=lambda row: float(row["gap"]))
    return rows


def csv_fieldnames(rows: list[dict[str, float | int | str]]) -> list[str]:
    fieldnames: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row:
            if key not in seen:
                seen.add(key)
                fieldnames.append(key)
    return fieldnames


def write_tracking_csv(
    *,
    rows: list[dict[str, float | int | str]],
    warnings: list[dict[str, float | int | str]],
    output_csv: Path,
    output_warnings_csv: Path,
) -> None:
    fieldnames = csv_fieldnames(rows)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    with output_warnings_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(warnings)


def existing_beta15_warning_count() -> int | None:
    report_path = (
        DEFAULT_OUTPUT_DIR
        / "thickness_mismatch_lambda_mu_beta15_eps0p0025_eta_m0p5_0_p0p5_descendants_report.md"
    )
    if not report_path.exists():
        return None
    for line in report_path.read_text(encoding="utf-8").splitlines():
        prefix = "- total warning rows: "
        if line.startswith(prefix):
            return int(line.removeprefix(prefix).strip())
    return None


def plot_cases(
    *,
    cases: dict[float, dict[str, object]],
    beta_deg: float,
    epsilon: float,
    eta_values: tuple[float, ...],
    mu_values: np.ndarray,
    output_png: Path,
) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, len(eta_values), figsize=(11.4, 4.4), sharey=True)
    axes = np.atleast_1d(axes)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    for ax, eta in zip(axes, eta_values, strict=True):
        case = cases[float(eta)]
        tracking = case["tracking"]
        tracked = tracking.tracked
        valid = valid_mask(
            epsilon=float(epsilon),
            eta=float(eta),
            mu_values=mu_values,
            limit=THICKNESS_RATIO_LIMIT,
        )
        for branch_idx in range(N_DESCENDANTS_PLOT):
            plot_with_validity_split(
                ax,
                mu_values,
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
        rf"Thickness-mismatch descendant branches, $\beta = {float(beta_deg):g}^\circ$, $\epsilon = {float(epsilon):g}$",
        fontsize=12,
    )
    fig.tight_layout(rect=(0.0, 0.17, 1.0, 0.92))
    fig.savefig(output_png, dpi=240, bbox_inches="tight")
    plt.close(fig)


def write_report(
    *,
    cases: dict[float, dict[str, object]],
    beta_deg: float,
    epsilon: float,
    eta_values: tuple[float, ...],
    mu_values: np.ndarray,
    output_png: Path,
    output_report: Path,
    output_csv: Path,
    output_warnings_csv: Path,
) -> dict[str, object]:
    all_rows: list[dict[str, float | int | str]] = []
    close_rows_by_eta: dict[float, list[dict[str, float | int | str]]] = {}
    lines = [
        "# Thickness-Mismatch Lambda(mu) Descendant Diagnostic",
        "",
        "## Parameters",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- beta: {float(beta_deg):g} deg",
        f"- epsilon: {float(epsilon):g}",
        f"- eta values: {', '.join(f'{float(value):g}' for value in eta_values)}",
        f"- mu range: {MU_MIN:g} .. {MU_MAX:g}",
        f"- mu step: {MU_STEP:g}",
        f"- plotted descendants: first {N_DESCENDANTS_PLOT}",
        f"- tracked descendants: first {N_DESCENDANTS_TRACK}",
        f"- sorted roots scanned at each mu: primary {N_SORTED_SCAN}; fallback {N_SORTED_SCAN_FALLBACK} if high diagnostic roots are missing",
        "",
        "Branches are descendant branches seeded at `mu=0`; sorted positions are",
        "diagnostic metadata only. Tracking uses adjacent-step analytic shape MAC.",
        "Low-MAC, low-margin, or large-jump candidates are warnings and do not",
        "automatically rename a branch.",
        "Tracking warnings are retained in this report but are not drawn on the",
        "presentation-style PNG; the figure shows only descendant branches and",
        "the solid/dashed thin-rod applicability split.",
        "",
        f"- PNG: `{output_png.relative_to(REPO_ROOT)}`",
        f"- tracking CSV: `{output_csv.relative_to(REPO_ROOT)}`",
        f"- warning CSV: `{output_warnings_csv.relative_to(REPO_ROOT)}`",
    ]

    for eta in eta_values:
        case = cases[float(eta)]
        lines.append(
            f"- root scan eta={float(eta):g}: used {int(case['n_sorted_scan'])} sorted roots per mu "
            f"({case['root_scan_note']})"
        )

    lines.extend(
        [
            "",
            "## Thin-Rod Applicability",
            "",
            f"- criterion: `2*r_i/l_i <= {THICKNESS_RATIO_LIMIT:g}` for both rods",
        ]
    )

    for eta in eta_values:
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
        close_rows_by_eta[float(eta)] = close_approach_rows(
            eta=float(eta),
            mu_values=mu_values,
            tracked=case["tracking"].tracked,
        )

    warnings = tracking_warning_rows(all_rows)
    by_eta: dict[float, list[dict[str, float | int | str]]] = defaultdict(list)
    for row in warnings:
        by_eta[float(row["eta"])].append(row)
    lines.extend(["", "## Tracking Warnings", ""])
    if warnings:
        lines.append(f"- total warning rows: {len(warnings)}")
        for eta in eta_values:
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

    write_tracking_csv(
        rows=all_rows,
        warnings=warnings,
        output_csv=output_csv,
        output_warnings_csv=output_warnings_csv,
    )
    lines.extend(
        [
            "",
            "## CSV Outputs",
            "",
            f"- tracking rows: {len(all_rows)}",
            f"- warning rows: {len(warnings)}",
        ]
    )

    lines.extend(
        [
            "",
            "## Near-Crossing / Close-Approach Diagnostics",
            "",
            f"- close-approach flag: min gap <= {CLOSE_APPROACH_ABS_THRESHOLD:g} or relative gap <= {CLOSE_APPROACH_REL_THRESHOLD:g}",
            "- these are visual/diagnostic gap checks between the plotted descendant branches, not sorted-root identity changes.",
        ]
    )
    for eta in eta_values:
        eta_rows = close_rows_by_eta[float(eta)]
        flagged = [row for row in eta_rows if str(row["close_approach_flag"]) == "yes"]
        sign_changes = [row for row in eta_rows if str(row["sign_change_on_grid"]) == "yes"]
        lines.append(
            f"- eta={float(eta):g}: close flags={len(flagged)}, "
            f"sampled sign changes={len(sign_changes)}; three smallest gaps:"
        )
        for row in eta_rows[:3]:
            lines.append(
                f"  - desc {int(row['desc_left'])}/desc {int(row['desc_right'])}: "
                f"gap={float(row['gap']):.6g} at mu={float(row['mu']):g}, "
                f"rel={float(row['relative_gap']):.6g}, "
                f"sign_change={row['sign_change_on_grid']}, close_flag={row['close_approach_flag']}"
            )

    if not np.isclose(float(beta_deg), 15.0, rtol=0.0, atol=1e-12):
        beta15_warnings = existing_beta15_warning_count()
        lines.extend(["", "## Comparison With Existing Beta=15 Diagnostic", ""])
        lines.append(
            "- no beta=15 recomputation was run for this comparison; it reads the existing report when available."
        )
        lines.append(
            "- the thin-rod applicability split is the same for this eta/epsilon/mu grid because the project criterion "
            "`2*r_i/l_i <= 0.1` does not depend on beta."
        )
        if beta15_warnings is None:
            lines.append("- existing beta=15 warning count was not found in the local report.")
        else:
            lines.append(
                f"- existing beta=15 total tracking warning rows: {beta15_warnings}; "
                f"this beta={float(beta_deg):g} run total: {len(warnings)}."
            )

    lines.extend(
        [
            "",
            "This is an analytic-only Euler-Bernoulli thickness-mismatch diagnostic.",
            "It does not add Timoshenko curves, FEM markers, Gmsh, CalculiX, or 3D FEM runs.",
            "This diagnostic writes only results files and does not modify article",
            "files, article figures, the baseline determinant, old solvers, or the",
            "FEM physical model.",
            "",
        ]
    )
    output_report.parent.mkdir(parents=True, exist_ok=True)
    output_report.write_text("\n".join(lines), encoding="utf-8")
    return {"warning_rows": len(warnings), "rows": all_rows, "close_rows_by_eta": close_rows_by_eta}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot thickness-mismatch EB descendant Lambda(mu) branches for one or more beta angles. "
            "The default reproduces the historical beta=15 diagnostic output."
        )
    )
    parser.add_argument(
        "--betas",
        type=float,
        nargs="+",
        default=list(DEFAULT_BETA_DEG_VALUES),
        help="Beta angles in degrees. Default: 15.",
    )
    parser.add_argument(
        "--epsilon",
        type=float,
        default=DEFAULT_EPSILON,
        help="Base radius parameter epsilon. Default: 0.0025.",
    )
    parser.add_argument(
        "--eta-values",
        type=float,
        nargs="+",
        default=list(DEFAULT_ETA_VALUES),
        help="Eta panel values. Default: -0.5 0 0.5.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for PNG and Markdown outputs. Default: results.",
    )
    return parser.parse_args(argv)


def run_beta_case(
    *,
    beta_deg: float,
    epsilon: float,
    eta_values: tuple[float, ...],
    output_dir: Path,
) -> dict[str, object]:
    mu_values = mu_grid(MU_MIN, MU_MAX, MU_STEP)
    output_png, output_report, output_csv, output_warnings_csv = output_paths(
        output_dir=output_dir,
        beta_deg=float(beta_deg),
        epsilon=float(epsilon),
        eta_values=eta_values,
    )
    cases = {
        float(eta): compute_case(
            beta_deg=float(beta_deg),
            epsilon=float(epsilon),
            eta=float(eta),
            mu_values=mu_values,
        )
        for eta in eta_values
    }
    plot_cases(
        cases=cases,
        beta_deg=float(beta_deg),
        epsilon=float(epsilon),
        eta_values=eta_values,
        mu_values=mu_values,
        output_png=output_png,
    )
    report = write_report(
        cases=cases,
        beta_deg=float(beta_deg),
        epsilon=float(epsilon),
        eta_values=eta_values,
        mu_values=mu_values,
        output_png=output_png,
        output_report=output_report,
        output_csv=output_csv,
        output_warnings_csv=output_warnings_csv,
    )

    print(f"saved descendant Lambda(mu) PNG: {output_png}")
    print(f"saved descendant Lambda(mu) report: {output_report}")
    print(f"saved descendant Lambda(mu) tracking CSV: {output_csv}")
    print(f"saved descendant Lambda(mu) warning CSV: {output_warnings_csv}")
    for eta in eta_values:
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
        close_rows = report["close_rows_by_eta"][float(eta)]
        close_flags = sum(1 for row in close_rows if str(row["close_approach_flag"]) == "yes")
        min_gap = float(close_rows[0]["gap"]) if close_rows else np.nan
        min_gap_mu = float(close_rows[0]["mu"]) if close_rows else np.nan
        print(
            f"eta={float(eta):g}: close-approach flags={close_flags}, "
            f"smallest plotted descendant gap={min_gap:.6g} at mu={min_gap_mu:g}"
        )

    return {"png": output_png, "report": output_report, "csv": output_csv, "warning_csv": output_warnings_csv, **report}


def main(argv: list[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    eta_values = tuple(float(value) for value in args.eta_values)
    results = [
        run_beta_case(
            beta_deg=float(beta),
            epsilon=float(args.epsilon),
            eta_values=eta_values,
            output_dir=Path(args.output_dir),
        )
        for beta in args.betas
    ]
    return {"results": results}


if __name__ == "__main__":
    main()
