from __future__ import annotations

import csv
from dataclasses import dataclass
import importlib.util
import math
from pathlib import Path
import sys
from types import ModuleType

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt


BETA_DEG = 15.0
MU = 0.0
ETA = 0.0
EPSILON_VALUES = [0.01, 0.025, 0.05]
N_BRANCHES = 6
MAC_ELIGIBLE = {"strong", "moderate"}
THIN_ROD_LIMIT = 0.1

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_DIR = REPO_ROOT / "results"
SOURCE_CSV = RESULTS_DIR / "3d_fem_rigid_joint_thickness_trend.csv"
OUTPUT_STEM = "rigid_joint_equal_thickness_beta15_frequency_comparison_eps_trend"
OUTPUT_CSV = RESULTS_DIR / f"{OUTPUT_STEM}.csv"
OUTPUT_PNG = RESULTS_DIR / f"{OUTPUT_STEM}.png"
OUTPUT_REPORT = RESULTS_DIR / f"{OUTPUT_STEM}_report.md"


@dataclass(frozen=True)
class PlotRow:
    epsilon: float
    branch: int
    lambda_eb: float
    lambda_timoshenko: float
    lambda_fem: float
    fem_solid_mode: str
    mac: float
    mac_strength: str
    rel_error_eb: float
    rel_error_timoshenko: float
    closer_model: str
    plotted: bool
    exclusion_reason: str
    eb_thin_rod_valid: bool
    diameter_to_length: float
    omega_over_cutoff: float


def finite(value: object) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float("nan")
    return number if math.isfinite(number) else float("nan")


def bool_token(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def fmt_float(value: float, digits: int = 6) -> str:
    if not math.isfinite(float(value)):
        return ""
    return f"{float(value):.{digits}g}"


def fmt_percent(value: float) -> str:
    if not math.isfinite(float(value)):
        return "n/a"
    return f"{100.0 * float(value):.2f}%"


def eps_key(value: float) -> str:
    return f"{float(value):.12g}"


def load_csv_rows() -> list[dict[str, str]]:
    if not SOURCE_CSV.exists():
        raise FileNotFoundError(
            f"Missing source CSV {SOURCE_CSV}. Run the rigid-joint thickness trend audit first."
        )
    with SOURCE_CSV.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def load_script_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {name} from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def analytic_lookup_from_audit() -> dict[tuple[str, int], tuple[float, float, float]]:
    """Analytic-only fallback; this does not run Gmsh or CalculiX."""

    module = load_script_module(
        "rigid_joint_thickness_trend_audit_for_plot_fallback",
        REPO_ROOT / "scripts" / "analysis" / "audit_3d_fem_rigid_joint_thickness_trend.py",
    )
    module.EPSILON_VALUES = list(EPSILON_VALUES)
    module.N_ANALYTIC_MODES = max(N_BRANCHES, getattr(module, "N_ANALYTIC_MODES", N_BRANCHES))
    pj = module.load_point_joint_module()
    module.configure_point_joint_module(pj)
    eb_reference = module.load_eb_reference_module()
    modes, warnings = module.compute_analytic_modes(pj, eb_reference)
    if warnings:
        print("Analytic fallback warnings:")
        for warning in warnings:
            print(f"- {warning}")
    lookup: dict[tuple[str, int], tuple[float, float, float]] = {}
    for mode in modes:
        lookup[(eps_key(mode.epsilon), int(mode.mode))] = (
            float(mode.lambda_eb),
            float(mode.lambda_timo),
            float(mode.omega_timo_over_cutoff),
        )
    return lookup


def build_plot_rows(source_rows: list[dict[str, str]]) -> list[PlotRow]:
    candidates: dict[tuple[str, int], dict[str, str]] = {}
    for row in source_rows:
        if row.get("row_kind") != "mac_match":
            continue
        epsilon = finite(row.get("epsilon"))
        branch = int(finite(row.get("analytic_mode"))) if math.isfinite(finite(row.get("analytic_mode"))) else -1
        if eps_key(epsilon) not in {eps_key(value) for value in EPSILON_VALUES}:
            continue
        if branch < 1 or branch > N_BRANCHES:
            continue
        candidates[(eps_key(epsilon), branch)] = row

    missing_keys = [
        (eps_key(epsilon), branch)
        for epsilon in EPSILON_VALUES
        for branch in range(1, N_BRANCHES + 1)
        if (eps_key(epsilon), branch) not in candidates
    ]
    if missing_keys:
        missing_text = ", ".join(f"epsilon={epsilon}, branch={branch}" for epsilon, branch in missing_keys)
        raise RuntimeError(
            "The trend CSV is missing required MAC-matched FEM rows: "
            f"{missing_text}. Re-run the diagnostic trend audit before plotting."
        )

    need_analytic_fallback = any(
        not math.isfinite(finite(row.get("Lambda_EB")))
        or not math.isfinite(finite(row.get("Lambda_Timoshenko")))
        or not math.isfinite(finite(row.get("Omega_Timoshenko_over_cutoff")))
        for row in candidates.values()
    )
    analytic_lookup = analytic_lookup_from_audit() if need_analytic_fallback else {}

    plot_rows: list[PlotRow] = []
    for epsilon in EPSILON_VALUES:
        for branch in range(1, N_BRANCHES + 1):
            row = candidates[(eps_key(epsilon), branch)]
            lambda_eb = finite(row.get("Lambda_EB"))
            lambda_timoshenko = finite(row.get("Lambda_Timoshenko"))
            omega_over_cutoff = finite(row.get("Omega_Timoshenko_over_cutoff"))
            if not (
                math.isfinite(lambda_eb)
                and math.isfinite(lambda_timoshenko)
                and math.isfinite(omega_over_cutoff)
            ):
                fallback = analytic_lookup.get((eps_key(epsilon), branch))
                if fallback is None:
                    raise RuntimeError(f"Could not recover analytic values for epsilon={epsilon}, branch={branch}")
                lambda_eb, lambda_timoshenko, omega_over_cutoff = fallback

            mac_strength = row.get("best_MAC_strength", "").strip().lower()
            ambiguous = bool_token(row.get("different_best_solid_modes"))
            lambda_fem = finite(row.get("Lambda_FEM"))
            mac = finite(row.get("selected_MAC"))
            if ambiguous:
                plotted = False
                exclusion_reason = "ambiguous_different_best_solid_modes"
            elif mac_strength not in MAC_ELIGIBLE:
                plotted = False
                exclusion_reason = "weak_mac" if mac_strength == "weak" else "missing_mac_strength"
            elif not math.isfinite(lambda_fem):
                plotted = False
                exclusion_reason = "missing_fem_lambda"
            else:
                plotted = True
                exclusion_reason = ""

            diameter_to_length = 4.0 * float(epsilon)
            plot_rows.append(
                PlotRow(
                    epsilon=float(epsilon),
                    branch=branch,
                    lambda_eb=lambda_eb,
                    lambda_timoshenko=lambda_timoshenko,
                    lambda_fem=lambda_fem,
                    fem_solid_mode=row.get("selected_solid_mode", "").strip(),
                    mac=mac,
                    mac_strength=mac_strength,
                    rel_error_eb=finite(row.get("rel_error_EB")),
                    rel_error_timoshenko=finite(row.get("rel_error_Timoshenko")),
                    closer_model=row.get("closer_model_after_MAC", "").strip(),
                    plotted=plotted,
                    exclusion_reason=exclusion_reason,
                    eb_thin_rod_valid=diameter_to_length <= THIN_ROD_LIMIT + 1.0e-12,
                    diameter_to_length=diameter_to_length,
                    omega_over_cutoff=omega_over_cutoff,
                )
            )
    return plot_rows


def write_plot_csv(rows: list[PlotRow]) -> None:
    fieldnames = [
        "epsilon",
        "branch",
        "Lambda_EB",
        "Lambda_Timoshenko",
        "Lambda_FEM",
        "FEM_solid_mode",
        "MAC",
        "MAC_strength",
        "rel_error_EB",
        "rel_error_Timoshenko",
        "closer_model",
        "plotted",
        "exclusion_reason",
        "EB_thin_rod_valid",
        "diameter_to_length",
        "Omega_over_cutoff",
    ]
    with OUTPUT_CSV.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "epsilon": fmt_float(row.epsilon, 8),
                    "branch": row.branch,
                    "Lambda_EB": fmt_float(row.lambda_eb, 12),
                    "Lambda_Timoshenko": fmt_float(row.lambda_timoshenko, 12),
                    "Lambda_FEM": fmt_float(row.lambda_fem, 12),
                    "FEM_solid_mode": row.fem_solid_mode,
                    "MAC": fmt_float(row.mac, 12),
                    "MAC_strength": row.mac_strength,
                    "rel_error_EB": fmt_float(row.rel_error_eb, 12),
                    "rel_error_Timoshenko": fmt_float(row.rel_error_timoshenko, 12),
                    "closer_model": row.closer_model,
                    "plotted": str(row.plotted),
                    "exclusion_reason": row.exclusion_reason,
                    "EB_thin_rod_valid": str(row.eb_thin_rod_valid),
                    "diameter_to_length": fmt_float(row.diameter_to_length, 8),
                    "Omega_over_cutoff": fmt_float(row.omega_over_cutoff, 12),
                }
            )


def rows_for_epsilon(rows: list[PlotRow], epsilon: float) -> list[PlotRow]:
    return [row for row in rows if abs(row.epsilon - float(epsilon)) <= 1.0e-12]


def eligible_rows(rows: list[PlotRow], epsilon: float) -> list[PlotRow]:
    return [row for row in rows_for_epsilon(rows, epsilon) if row.plotted]


def mean(values: list[float]) -> float:
    finite_values = [float(value) for value in values if math.isfinite(float(value))]
    if not finite_values:
        return float("nan")
    return sum(finite_values) / float(len(finite_values))


def trend_summary(rows: list[PlotRow], epsilon: float) -> dict[str, object]:
    eligible = eligible_rows(rows, epsilon)
    mean_eb = mean([row.rel_error_eb for row in eligible])
    mean_timo = mean([row.rel_error_timoshenko for row in eligible])
    ratio = mean_timo / mean_eb if mean_eb > 0.0 and math.isfinite(mean_eb) else float("nan")
    return {
        "eligible": len(eligible),
        "weak": sum(1 for row in rows_for_epsilon(rows, epsilon) if row.exclusion_reason == "weak_mac"),
        "ambiguous": sum(
            1 for row in rows_for_epsilon(rows, epsilon) if row.exclusion_reason.startswith("ambiguous")
        ),
        "mean_eb": mean_eb,
        "mean_timo": mean_timo,
        "ratio": ratio,
        "closer_eb": sum(1 for row in eligible if row.closer_model == "EB"),
        "closer_timo": sum(1 for row in eligible if row.closer_model == "Timoshenko"),
        "max_omega_over_cutoff": max(
            [row.omega_over_cutoff for row in rows_for_epsilon(rows, epsilon) if math.isfinite(row.omega_over_cutoff)],
            default=float("nan"),
        ),
    }


def eb_status(row: PlotRow) -> str:
    if row.diameter_to_length < THIN_ROD_LIMIT - 1.0e-12:
        return "valid"
    if abs(row.diameter_to_length - THIN_ROD_LIMIT) <= 1.0e-12:
        return "boundary-valid"
    return "outside"


def cutoff_status(max_ratio: float) -> str:
    if not math.isfinite(float(max_ratio)):
        return "unknown"
    if max_ratio >= 1.0:
        return "violation"
    if max_ratio >= 0.8:
        return "warning"
    return "ok"


def make_plot(rows: list[PlotRow]) -> None:
    fig, axes = plt.subplots(1, len(EPSILON_VALUES), figsize=(12.0, 4.0), sharey=True, constrained_layout=True)
    if len(EPSILON_VALUES) == 1:
        axes = [axes]

    for axis, epsilon in zip(axes, EPSILON_VALUES):
        eps_rows = rows_for_epsilon(rows, epsilon)
        branches = [row.branch for row in eps_rows]
        axis.plot(
            branches,
            [row.lambda_eb for row in eps_rows],
            color="#1f77b4",
            marker="o",
            linewidth=1.8,
            markersize=5,
            label="Euler--Bernoulli",
        )
        axis.plot(
            branches,
            [row.lambda_timoshenko for row in eps_rows],
            color="#d95f02",
            marker="s",
            linewidth=1.8,
            markersize=5,
            label="Timoshenko",
        )
        plotted = [row for row in eps_rows if row.plotted]
        excluded = [row for row in eps_rows if not row.plotted and math.isfinite(row.lambda_fem)]
        axis.scatter(
            [row.branch for row in plotted],
            [row.lambda_fem for row in plotted],
            color="#111111",
            marker="D",
            s=58,
            zorder=4,
            label="3D FEM MAC-eligible",
        )
        if excluded:
            axis.scatter(
                [row.branch for row in excluded],
                [row.lambda_fem for row in excluded],
                color="#777777",
                marker="x",
                s=46,
                alpha=0.35,
                linewidths=1.4,
                zorder=3,
                label="3D FEM weak/ambiguous",
            )

        summary = trend_summary(rows, epsilon)
        axis.set_title(
            f"epsilon={epsilon:g}\n"
            f"Timo/EB mean error={fmt_float(float(summary['ratio']), 4)}"
        )
        axis.set_xlabel("branch / mode index")
        axis.set_xticks(range(1, N_BRANCHES + 1))
        axis.grid(True, linewidth=0.5, alpha=0.3)
    axes[0].set_ylabel("Lambda")

    handles, labels = axes[0].get_legend_handles_labels()
    if len(EPSILON_VALUES) > 1:
        extra_handles, extra_labels = axes[-1].get_legend_handles_labels()
        for handle, label in zip(extra_handles, extra_labels):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    fig.legend(handles, labels, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.08))
    fig.suptitle("Rigid end-face 3D FEM thickness trend at beta=15 deg", y=1.18, fontsize=12)
    fig.savefig(OUTPUT_PNG, dpi=240, bbox_inches="tight")
    plt.close(fig)


def markdown_table(headers: list[str], rows: list[list[str]]) -> list[str]:
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join("---" for _ in headers) + " |",
    ]
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")
    return lines


def write_report(rows: list[PlotRow]) -> None:
    summary_by_eps = {epsilon: trend_summary(rows, epsilon) for epsilon in EPSILON_VALUES}
    excluded_rows = [row for row in rows if not row.plotted]

    lines: list[str] = [
        "# Rigid-Joint Equal-Thickness beta=15 Frequency Trend",
        "",
        "This diagnostic plot uses the existing rigid end-face thickness-trend run. "
        "Gmsh and CalculiX are not rerun by this plotting script.",
        "",
        "## Data Source",
        "",
        f"- Source CSV: `{SOURCE_CSV.relative_to(REPO_ROOT)}`",
        f"- Plot-data CSV: `{OUTPUT_CSV.relative_to(REPO_ROOT)}`",
        f"- PNG: `{OUTPUT_PNG.relative_to(REPO_ROOT)}`",
        "",
        "The source FEM formulation is the full rigid end-face engineering point-joint: "
        "two separate solid cylinders, outer ends clamped, all inner end-face nodes "
        "coupled by CalculiX `*RIGID BODY` to one common `JOINT_REF`, with no planar "
        "constraint, no patch tuning, and no fused volume.",
        "",
        "This is not a patch-radius calibration. Exact coincidence with a 1D analytic "
        "point-joint is not expected; the diagnostic question is comparative, namely "
        "whether the Timoshenko reference becomes closer than Euler--Bernoulli as "
        "rod thickness increases.",
        "",
        "## Parameters",
        "",
        f"- `beta = {BETA_DEG:g} deg`",
        f"- `mu = {MU:g}`",
        f"- `eta = {ETA:g}`",
        f"- `epsilon = {', '.join(fmt_float(value, 6) for value in EPSILON_VALUES)}`",
        f"- branches plotted: first `{N_BRANCHES}` analytic modes",
        "",
        "## Matching Policy",
        "",
        "FEM points use the MAC-like centerline match recorded in the source trend audit. "
        "Rows with `strong` or `moderate` MAC are main plotted FEM points. Weak or "
        "ambiguous rows are not used in the main comparison; they are shown only as "
        "faint diagnostic markers in the PNG and are marked `plotted=False` in the "
        "plot-data CSV.",
        "",
        "## Excluded FEM Rows",
        "",
    ]

    if excluded_rows:
        lines.extend(
            markdown_table(
                ["epsilon", "branch", "FEM mode", "MAC", "strength", "reason"],
                [
                    [
                        fmt_float(row.epsilon, 6),
                        str(row.branch),
                        row.fem_solid_mode,
                        fmt_float(row.mac, 5),
                        row.mac_strength or "missing",
                        row.exclusion_reason,
                    ]
                    for row in excluded_rows
                ],
            )
        )
    else:
        lines.append("No branches in the plotted range were excluded.")

    lines.extend(
        [
            "",
            "Branches above the plotted range are intentionally outside this figure, not exclusions.",
            "",
            "## EB Applicability",
            "",
        ]
    )
    applicability_rows = []
    for epsilon in EPSILON_VALUES:
        first = rows_for_epsilon(rows, epsilon)[0]
        applicability_rows.append(
            [
                fmt_float(epsilon, 6),
                fmt_float(first.diameter_to_length, 4),
                eb_status(first),
            ]
        )
    lines.extend(markdown_table(["epsilon", "2*r/l = 4*epsilon", "EB thin-rod status"], applicability_rows))
    lines.extend(
        [
            "",
            "The EB diameter criterion is reported separately from the Timoshenko cut-off. "
            "The plot uses consistent line styles; it does not draw vertical applicability "
            "lines and does not dash Timoshenko curves because of the EB thin-rod criterion.",
            "",
            "## Timoshenko Cut-Off",
            "",
        ]
    )
    cutoff_rows = []
    for epsilon in EPSILON_VALUES:
        max_ratio = float(summary_by_eps[epsilon]["max_omega_over_cutoff"])
        cutoff_rows.append(
            [
                fmt_float(epsilon, 6),
                fmt_float(max_ratio, 5),
                cutoff_status(max_ratio),
            ]
        )
    lines.extend(markdown_table(["epsilon", "max Omega/Omega_c", "status"], cutoff_rows))

    lines.extend(
        [
            "",
            "## Trend Summary",
            "",
        ]
    )
    trend_rows = []
    for epsilon in EPSILON_VALUES:
        summary = summary_by_eps[epsilon]
        trend_rows.append(
            [
                fmt_float(epsilon, 6),
                str(summary["eligible"]),
                fmt_percent(float(summary["mean_eb"])),
                fmt_percent(float(summary["mean_timo"])),
                fmt_float(float(summary["ratio"]), 5),
                str(summary["closer_eb"]),
                str(summary["closer_timo"]),
            ]
        )
    lines.extend(
        markdown_table(
            [
                "epsilon",
                "eligible",
                "mean rel error vs EB",
                "mean rel error vs Timo",
                "Timo/EB",
                "closer EB",
                "closer Timo",
            ],
            trend_rows,
        )
    )

    ratios = [float(summary_by_eps[epsilon]["ratio"]) for epsilon in EPSILON_VALUES]
    decreasing = all(
        math.isfinite(left) and math.isfinite(right) and right < left
        for left, right in zip(ratios, ratios[1:])
    )
    thick_ratio = ratios[-1]
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "Within the first-six-branch plotting window, the Timoshenko/EB mean-error "
            f"ratio is {'monotonically decreasing' if decreasing else 'not strictly monotone'} "
            "over the selected epsilon values. The thickest case is the clearest: "
            f"at `epsilon={EPSILON_VALUES[-1]:g}`, the ratio is `{fmt_float(thick_ratio, 5)}`, "
            "so the Timoshenko reference is substantially closer to the full rigid "
            "end-face 3D FEM benchmark than Euler--Bernoulli for the eligible modes.",
            "",
            "This is diagnostic validation evidence, not final article proof. Mesh "
            "convergence, visual mode-shape review, and a final assignment policy remain "
            "separate prerequisites.",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    source_rows = load_csv_rows()
    rows = build_plot_rows(source_rows)
    write_plot_csv(rows)
    make_plot(rows)
    write_report(rows)
    print(f"Read {SOURCE_CSV.relative_to(REPO_ROOT)}")
    print("Did not rerun Gmsh or CalculiX")
    print(f"Wrote {OUTPUT_CSV.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_PNG.relative_to(REPO_ROOT)}")
    print(f"Wrote {OUTPUT_REPORT.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
