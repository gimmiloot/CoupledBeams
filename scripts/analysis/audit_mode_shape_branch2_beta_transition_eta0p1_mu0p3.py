from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
import sys
from typing import Sequence

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

from my_project.analytic.formulas_thickness_mismatch import assemble_clamped_coupled_matrix_eta  # noqa: E402
from scripts.analysis import plot_mode_shapes_eta_beta_scan as shape_scan  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import analytic_null_vector  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import mac_value  # noqa: E402


DEFAULT_ETA = 0.1
DEFAULT_MU = 0.3
DEFAULT_EPSILON = 0.0025
DEFAULT_OVERVIEW_BETA_START = 0.0
DEFAULT_OVERVIEW_BETA_END = 11.0
DEFAULT_OVERVIEW_BETA_STEP = 1.0
DEFAULT_REFINED_BETA_START = 4.0
DEFAULT_REFINED_BETA_END = 7.0
DEFAULT_REFINED_BETA_STEP = 0.05
DEFAULT_NUM_SORTED_ROOTS = 8
DEFAULT_NUM_TRACKED_BRANCHES = 8
DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
DEFAULT_NUM_SHAPE_SAMPLES = 401
DEFAULT_MAC_WARNING_THRESHOLD = 0.9
DEFAULT_MAC_MARGIN_WARNING_THRESHOLD = 0.05
DEFAULT_MAX_SORTED_POSITION_JUMP = 1
DEFAULT_NEAR_ZERO_GAP_TOL = 1e-7
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"

TARGET_BRANCHES = (2, 3)
CSV_SORTED_INDICES = tuple(range(1, 7))
REQUIRED_STATEMENT = (
    "A descendant-branch sorted-index change is not by itself a frequency crossing. "
    "A crossing claim requires Lambda_3^sorted - Lambda_2^sorted = 0 or unresolved near-zero gap."
)

CSV_FIELDNAMES = [
    "beta_deg",
    "sorted_index",
    "Lambda_sorted",
    "gap_to_next",
    "descendant_branch_2_sorted_index",
    "descendant_branch_2_Lambda",
    "descendant_branch_3_sorted_index",
    "descendant_branch_3_Lambda",
    "gap_23",
    "mac_branch2_to_previous",
    "mac_branch2_to_beta0",
    "mac_branch3_to_previous",
    "mac_branch3_to_beta0",
    "coeff_norm_branch2",
    "coeff_norm_branch3",
    "warning",
    "notes",
]


@dataclass(frozen=True)
class BranchDiagnostics:
    branch: int
    beta_deg: float
    sorted_index: int
    Lambda: float
    mac_to_previous: float
    mac_to_beta0: float
    coeff_norm: float
    tracking_warning: str


@dataclass(frozen=True)
class AuditResult:
    beta_values: np.ndarray
    sorted_roots: np.ndarray
    branch_rows: dict[int, list[shape_scan.TrackingRow]]
    branch_diagnostics: dict[int, dict[float, BranchDiagnostics]]
    csv_path: Path
    sorted_plot_path: Path
    gap_plot_path: Path
    tracking_map_path: Path
    mac_plot_path: Path
    report_path: Path
    min_gap_23: float
    min_gap_beta: float
    branch2_changes_sorted_index: bool
    branch2_change_beta: float | None
    gap_23_positive: bool
    tracking_warning_count: int
    classification: str


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Audit the sorted-frequency gap behind the eta=0.1, mu=0.3 "
            "descendant-branch-2 beta-scan sorted-index change."
        ),
    )
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--overview-beta-start", type=float, default=DEFAULT_OVERVIEW_BETA_START)
    parser.add_argument("--overview-beta-end", type=float, default=DEFAULT_OVERVIEW_BETA_END)
    parser.add_argument("--overview-beta-step", type=float, default=DEFAULT_OVERVIEW_BETA_STEP)
    parser.add_argument("--refined-beta-start", type=float, default=DEFAULT_REFINED_BETA_START)
    parser.add_argument("--refined-beta-end", type=float, default=DEFAULT_REFINED_BETA_END)
    parser.add_argument("--refined-beta-step", type=float, default=DEFAULT_REFINED_BETA_STEP)
    parser.add_argument("--num-sorted-roots", type=int, default=DEFAULT_NUM_SORTED_ROOTS)
    parser.add_argument("--num-tracked-branches", type=int, default=DEFAULT_NUM_TRACKED_BRANCHES)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--root-lmax0", type=float, default=DEFAULT_ROOT_LMAX0)
    parser.add_argument("--num-shape-samples", type=int, default=DEFAULT_NUM_SHAPE_SAMPLES)
    parser.add_argument("--mac-warning-threshold", type=float, default=DEFAULT_MAC_WARNING_THRESHOLD)
    parser.add_argument("--mac-margin-warning-threshold", type=float, default=DEFAULT_MAC_MARGIN_WARNING_THRESHOLD)
    parser.add_argument("--max-sorted-position-jump", type=int, default=DEFAULT_MAX_SORTED_POSITION_JUMP)
    parser.add_argument("--near-zero-gap-tol", type=float, default=DEFAULT_NEAR_ZERO_GAP_TOL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    validate_args(args)
    args.output_dir = resolve_repo_path(args.output_dir)
    return args


def validate_args(args: argparse.Namespace) -> None:
    if not (-1.0 < float(args.mu) < 1.0):
        raise ValueError("--mu must lie inside (-1, 1).")
    if not (-1.0 < float(args.eta) < 1.0):
        raise ValueError("--eta must lie inside (-1, 1).")
    if float(args.epsilon) <= 0.0:
        raise ValueError("--epsilon must be positive.")
    if float(args.overview_beta_start) < -1e-12:
        raise ValueError("--overview-beta-start must be non-negative.")
    if float(args.refined_beta_start) < -1e-12:
        raise ValueError("--refined-beta-start must be non-negative.")
    if float(args.overview_beta_step) <= 0.0 or float(args.refined_beta_step) <= 0.0:
        raise ValueError("beta steps must be positive.")
    if float(args.overview_beta_end) < float(args.overview_beta_start):
        raise ValueError("--overview-beta-end must be >= --overview-beta-start.")
    if float(args.refined_beta_end) < float(args.refined_beta_start):
        raise ValueError("--refined-beta-end must be >= --refined-beta-start.")
    if int(args.num_sorted_roots) < max(CSV_SORTED_INDICES) + 1:
        raise ValueError("--num-sorted-roots must include one extra root for gap_to_next.")
    if int(args.num_tracked_branches) < max(TARGET_BRANCHES):
        raise ValueError("--num-tracked-branches must include branches 2 and 3.")
    if int(args.num_sorted_roots) < int(args.num_tracked_branches):
        raise ValueError("--num-sorted-roots must be at least --num-tracked-branches.")
    if int(args.num_shape_samples) < 5:
        raise ValueError("--num-shape-samples must be at least 5.")
    if float(args.root_scan_step) <= 0.0 or float(args.root_lmax0) <= 0.0:
        raise ValueError("root scan parameters must be positive.")
    if float(args.near_zero_gap_tol) < 0.0:
        raise ValueError("--near-zero-gap-tol must be non-negative.")


def resolve_repo_path(path: str | Path) -> Path:
    value = Path(path)
    if value.is_absolute():
        return value
    return REPO_ROOT / value


def number_token(value: float) -> str:
    return shape_scan.number_token(float(value))


def fmt(value: float | int | str | None) -> str:
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    value_f = float(value)
    return "" if not np.isfinite(value_f) else f"{value_f:.12g}"


def path_stem(args: argparse.Namespace) -> str:
    return (
        f"mode_shape_branch2_transition_audit_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}"
    )


def output_paths(args: argparse.Namespace) -> dict[str, Path]:
    output_dir = Path(args.output_dir)
    param_suffix = (
        f"eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_eps_{number_token(float(args.epsilon))}"
    )
    return {
        "csv": output_dir / f"{path_stem(args)}.csv",
        "sorted_plot": output_dir / f"mode_shape_branch2_transition_sorted_lambda_{param_suffix}.png",
        "gap_plot": output_dir / f"mode_shape_branch2_transition_gap23_{param_suffix}.png",
        "tracking_map": output_dir / f"mode_shape_branch2_transition_tracking_map_{param_suffix}.png",
        "mac_plot": output_dir / f"mode_shape_branch2_transition_mac_{param_suffix}.png",
        "report": output_dir / f"{path_stem(args)}.md",
    }


def tracking_namespace(args: argparse.Namespace, *, branch: int, beta_step: float) -> argparse.Namespace:
    return argparse.Namespace(
        eta=float(args.eta),
        mu=float(args.mu),
        epsilon=float(args.epsilon),
        beta_step=float(beta_step),
        branch=int(branch),
        output_dir=Path(args.output_dir),
        num_tracked_branches=int(args.num_tracked_branches),
        num_sorted_roots=int(args.num_sorted_roots),
        root_scan_step=float(args.root_scan_step),
        root_lmax0=float(args.root_lmax0),
        num_shape_samples=int(args.num_shape_samples),
        mac_warning_threshold=float(args.mac_warning_threshold),
        mac_margin_warning_threshold=float(args.mac_margin_warning_threshold),
        max_sorted_position_jump=int(args.max_sorted_position_jump),
    )


def refined_tracking_grid(args: argparse.Namespace) -> np.ndarray:
    return shape_scan.inclusive_grid(0.0, float(args.refined_beta_end), float(args.refined_beta_step))


def refined_report_grid(args: argparse.Namespace, tracking_grid: np.ndarray) -> np.ndarray:
    mask = (
        tracking_grid >= float(args.refined_beta_start) - 1e-12
    ) & (
        tracking_grid <= float(args.refined_beta_end) + 1e-12
    )
    values = tracking_grid[mask]
    if values.size == 0:
        raise RuntimeError("Refined beta grid is empty.")
    return values


def lookup_tracking_rows(rows: Sequence[shape_scan.TrackingRow]) -> dict[float, shape_scan.TrackingRow]:
    return {round(float(row.beta_deg), 10): row for row in rows}


def coefficient_norm(args: argparse.Namespace, *, beta_deg: float, Lambda: float) -> float:
    matrix = assemble_clamped_coupled_matrix_eta(
        float(Lambda),
        float(np.deg2rad(float(beta_deg))),
        float(args.mu),
        float(args.epsilon),
        float(args.eta),
    )
    coeff, _smallest, _ratio = analytic_null_vector(matrix)
    return float(np.linalg.norm(coeff))


def build_branch_diagnostics(
    args: argparse.Namespace,
    *,
    beta_values: np.ndarray,
    sorted_roots: np.ndarray,
    branch_rows: dict[int, list[shape_scan.TrackingRow]],
    s_norm: np.ndarray,
) -> dict[int, dict[float, BranchDiagnostics]]:
    anchor_vectors = {}
    beta0_vectors = shape_scan.shape_vectors_for_beta(
        tracking_namespace(args, branch=2, beta_step=float(args.refined_beta_step)),
        sorted_roots[0],
        beta_deg=float(beta_values[0]),
        s_norm=s_norm,
    )
    for branch in TARGET_BRANCHES:
        anchor_vectors[int(branch)] = beta0_vectors[int(branch) - 1]

    roots_by_beta = {round(float(beta), 10): sorted_roots[idx] for idx, beta in enumerate(beta_values)}
    diagnostics: dict[int, dict[float, BranchDiagnostics]] = {int(branch): {} for branch in TARGET_BRANCHES}
    for branch in TARGET_BRANCHES:
        for row in branch_rows[int(branch)]:
            beta_key = round(float(row.beta_deg), 10)
            vectors = shape_scan.shape_vectors_for_beta(
                tracking_namespace(args, branch=int(branch), beta_step=float(args.refined_beta_step)),
                roots_by_beta[beta_key],
                beta_deg=float(row.beta_deg),
                s_norm=s_norm,
            )
            current_vector = vectors[int(row.sorted_index) - 1]
            diagnostics[int(branch)][beta_key] = BranchDiagnostics(
                branch=int(branch),
                beta_deg=float(row.beta_deg),
                sorted_index=int(row.sorted_index),
                Lambda=float(row.Lambda),
                mac_to_previous=float(row.mac_to_previous),
                mac_to_beta0=float(mac_value(anchor_vectors[int(branch)], current_vector)),
                coeff_norm=coefficient_norm(args, beta_deg=float(row.beta_deg), Lambda=float(row.Lambda)),
                tracking_warning=row.tracking_warning,
            )
    return diagnostics


def warning_for_row(
    args: argparse.Namespace,
    *,
    gap_23: float,
    branch2: BranchDiagnostics,
    branch3: BranchDiagnostics,
) -> str:
    warnings: list[str] = []
    if not np.isfinite(float(gap_23)):
        warnings.append("nonfinite_gap_23")
    elif float(gap_23) <= 0.0:
        warnings.append("nonpositive_gap_23")
    elif float(gap_23) <= float(args.near_zero_gap_tol):
        warnings.append("near_zero_gap_23")
    for label, branch in (("branch2", branch2), ("branch3", branch3)):
        if branch.tracking_warning:
            warnings.append(f"{label}:{branch.tracking_warning}")
    return ";".join(warnings) if warnings else "no"


def write_csv(
    args: argparse.Namespace,
    *,
    beta_values: np.ndarray,
    sorted_roots: np.ndarray,
    branch_diagnostics: dict[int, dict[float, BranchDiagnostics]],
    path: Path,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDNAMES)
        writer.writeheader()
        for beta_idx, beta_deg in enumerate(beta_values):
            beta_key = round(float(beta_deg), 10)
            roots = np.asarray(sorted_roots[beta_idx], dtype=float)
            gap_23 = float(roots[2] - roots[1])
            branch2 = branch_diagnostics[2][beta_key]
            branch3 = branch_diagnostics[3][beta_key]
            warning = warning_for_row(args, gap_23=gap_23, branch2=branch2, branch3=branch3)
            for sorted_index in CSV_SORTED_INDICES:
                lambda_sorted = float(roots[int(sorted_index) - 1])
                gap_to_next = (
                    float(roots[int(sorted_index)] - roots[int(sorted_index) - 1])
                    if int(sorted_index) < len(roots)
                    else float("nan")
                )
                writer.writerow(
                    {
                        "beta_deg": fmt(float(beta_deg)),
                        "sorted_index": int(sorted_index),
                        "Lambda_sorted": fmt(lambda_sorted),
                        "gap_to_next": fmt(gap_to_next),
                        "descendant_branch_2_sorted_index": int(branch2.sorted_index),
                        "descendant_branch_2_Lambda": fmt(branch2.Lambda),
                        "descendant_branch_3_sorted_index": int(branch3.sorted_index),
                        "descendant_branch_3_Lambda": fmt(branch3.Lambda),
                        "gap_23": fmt(gap_23),
                        "mac_branch2_to_previous": fmt(branch2.mac_to_previous),
                        "mac_branch2_to_beta0": fmt(branch2.mac_to_beta0),
                        "mac_branch3_to_previous": fmt(branch3.mac_to_previous),
                        "mac_branch3_to_beta0": fmt(branch3.mac_to_beta0),
                        "coeff_norm_branch2": fmt(branch2.coeff_norm),
                        "coeff_norm_branch3": fmt(branch3.coeff_norm),
                        "warning": warning,
                        "notes": (
                            "sorted-frequency row; descendant sorted index is diagnostic metadata; "
                            "crossing check uses gap_23"
                        ),
                    }
                )


def branch_positions(branch_diagnostics: dict[int, dict[float, BranchDiagnostics]], beta_values: np.ndarray, branch: int) -> np.ndarray:
    return np.array(
        [branch_diagnostics[int(branch)][round(float(beta), 10)].sorted_index for beta in beta_values],
        dtype=int,
    )


def branch_lambdas(branch_diagnostics: dict[int, dict[float, BranchDiagnostics]], beta_values: np.ndarray, branch: int) -> np.ndarray:
    return np.array(
        [branch_diagnostics[int(branch)][round(float(beta), 10)].Lambda for beta in beta_values],
        dtype=float,
    )


def branch_macs(
    branch_diagnostics: dict[int, dict[float, BranchDiagnostics]],
    beta_values: np.ndarray,
    branch: int,
    key: str,
) -> np.ndarray:
    values = []
    for beta in beta_values:
        item = branch_diagnostics[int(branch)][round(float(beta), 10)]
        values.append(getattr(item, key))
    return np.asarray(values, dtype=float)


def plot_sorted_lambda(
    args: argparse.Namespace,
    *,
    beta_values: np.ndarray,
    sorted_roots: np.ndarray,
    branch_diagnostics: dict[int, dict[float, BranchDiagnostics]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    lambda2 = sorted_roots[:, 1]
    lambda3 = sorted_roots[:, 2]
    branch2_lambda = branch_lambdas(branch_diagnostics, beta_values, 2)
    branch2_pos = branch_positions(branch_diagnostics, beta_values, 2)
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.plot(beta_values, lambda2, label="sorted 2", color="#006BA4", linewidth=1.8)
    ax.plot(beta_values, lambda3, label="sorted 3", color="#FF800E", linewidth=1.8)
    ax.scatter(
        beta_values,
        branch2_lambda,
        c=branch2_pos,
        cmap=plt.get_cmap("viridis", int(args.num_sorted_roots)),
        edgecolors="black",
        s=24,
        label="descendant 2 selected root",
        zorder=4,
    )
    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("Lambda")
    ax.set_title("Sorted Lambda_2/Lambda_3 with descendant-2 selected roots")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_gap(
    args: argparse.Namespace,
    *,
    beta_values: np.ndarray,
    gap_23: np.ndarray,
    min_idx: int,
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    ax.plot(beta_values, gap_23, color="#C85200", linewidth=1.9)
    ax.scatter([beta_values[min_idx]], [gap_23[min_idx]], color="black", s=30, zorder=5)
    ax.axhline(0.0, color="0.55", linestyle=":", linewidth=1.0)
    annotation = f"min gap={gap_23[min_idx]:.6g} at beta={beta_values[min_idx]:.3g} deg"
    if float(gap_23[min_idx]) > float(args.near_zero_gap_tol):
        annotation += "\npositive sorted gap: mode veering / avoided crossing candidate, not a crossing"
    else:
        annotation += "\nunresolved near-zero sorted gap"
    ax.annotate(
        annotation,
        xy=(beta_values[min_idx], gap_23[min_idx]),
        xytext=(0.03, 0.78),
        textcoords="axes fraction",
        arrowprops={"arrowstyle": "->", "color": "0.25", "linewidth": 0.8},
        fontsize=9,
    )
    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("Lambda_3^sorted - Lambda_2^sorted")
    ax.set_title("Sorted gap between roots 2 and 3")
    ax.grid(True, color="0.88", linewidth=0.6)
    fig.tight_layout()
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_tracking_map(
    *,
    beta_values: np.ndarray,
    branch_diagnostics: dict[int, dict[float, BranchDiagnostics]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 4.3))
    for branch, color in ((2, "#006BA4"), (3, "#FF800E")):
        ax.step(
            beta_values,
            branch_positions(branch_diagnostics, beta_values, branch),
            where="post",
            linewidth=1.9,
            marker="o",
            markersize=2.8,
            color=color,
            label=f"descendant {branch}",
        )
    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("selected sorted index")
    ax.set_yticks(range(1, 7))
    ax.set_title("Descendant selected sorted index over refined beta grid")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def plot_mac(
    *,
    beta_values: np.ndarray,
    branch_diagnostics: dict[int, dict[float, BranchDiagnostics]],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    styles = {
        (2, "mac_to_previous"): ("#006BA4", "-", "desc 2 to previous"),
        (2, "mac_to_beta0"): ("#006BA4", "--", "desc 2 to beta0"),
        (3, "mac_to_previous"): ("#FF800E", "-", "desc 3 to previous"),
        (3, "mac_to_beta0"): ("#FF800E", "--", "desc 3 to beta0"),
    }
    for (branch, key), (color, linestyle, label) in styles.items():
        values = branch_macs(branch_diagnostics, beta_values, branch, key)
        ax.plot(beta_values, values, color=color, linestyle=linestyle, linewidth=1.7, label=label)
    ax.set_xlabel("beta (deg)")
    ax.set_ylabel("MAC")
    ax.set_ylim(-0.02, 1.02)
    ax.set_title("Shape MAC diagnostics for descendants 2 and 3")
    ax.grid(True, color="0.88", linewidth=0.6)
    ax.legend(loc="best", frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output, dpi=230, bbox_inches="tight")
    plt.close(fig)


def first_change_beta(beta_values: np.ndarray, positions: np.ndarray, *, base_position: int) -> float | None:
    changed = np.flatnonzero(positions != int(base_position))
    if changed.size == 0:
        return None
    return float(beta_values[int(changed[0])])


def classify(
    *,
    branch2_changes: bool,
    gap_positive: bool,
    min_gap: float,
    near_zero_tol: float,
    tracking_warning_count: int,
) -> str:
    if not gap_positive or float(min_gap) <= float(near_zero_tol):
        return "unresolved near-zero/crossing candidate"
    if int(tracking_warning_count) > 0:
        return "positive-gap close approach with tracking ambiguity; mode-veering candidate unresolved, not a crossing"
    if branch2_changes:
        return "avoided crossing / mode-veering candidate, not a true crossing"
    return "resolved positive gap without descendant-2 sorted-index change"


def write_report(
    args: argparse.Namespace,
    *,
    result: AuditResult,
    overview_beta_values: np.ndarray,
    overview_branch2_rows: Sequence[shape_scan.TrackingRow],
) -> None:
    overview_positions = [int(row.sorted_index) for row in overview_branch2_rows]
    overview_change_beta = first_change_beta(
        overview_beta_values,
        np.asarray(overview_positions, dtype=int),
        base_position=2,
    )
    warning_lines = []
    for branch in TARGET_BRANCHES:
        for beta in result.beta_values:
            item = result.branch_diagnostics[int(branch)][round(float(beta), 10)]
            if item.tracking_warning:
                warning_lines.append(
                    f"- descendant {int(branch)}, beta={float(beta):.12g} deg: {item.tracking_warning}"
                )
    overview_lines = [
        f"- beta={float(beta):g} deg: descendant 2 sorted index {pos}"
        for beta, pos in zip(overview_beta_values, overview_positions, strict=True)
    ]
    lines = [
        "# Mode-Shape Branch 2 Beta Transition Audit",
        "",
        "## Parameters",
        "",
        f"- eta: `{float(args.eta):g}`",
        f"- mu: `{float(args.mu):g}`",
        f"- epsilon: `{float(args.epsilon):g}`",
        f"- overview beta range: `{float(args.overview_beta_start):g}..{float(args.overview_beta_end):g} deg`, step `{float(args.overview_beta_step):g} deg`",
        f"- refined beta range: `{float(args.refined_beta_start):g}..{float(args.refined_beta_end):g} deg`, step `{float(args.refined_beta_step):g} deg`",
        f"- sorted roots solved: `{int(args.num_sorted_roots)}`",
        f"- descendants tracked: `{int(args.num_tracked_branches)}`",
        "- model: diagnostic analytic Euler-Bernoulli thickness-mismatch determinant",
        "",
        "## Answers",
        "",
        f"- Does descendant branch 2 change sorted index on the refined grid? `{'yes' if result.branch2_changes_sorted_index else 'no'}`.",
        (
            f"- Change beta on the refined grid: `{result.branch2_change_beta:.12g} deg`."
            if result.branch2_change_beta is not None
            else "- Change beta on the refined grid: `none`."
        ),
        (
            f"- Coarse 1-degree overview first shows descendant branch 2 at a different sorted index at beta `{overview_change_beta:.12g} deg`."
            if overview_change_beta is not None
            else "- Coarse 1-degree overview shows no descendant-2 sorted-index change."
        ),
        f"- Minimum `gap_23 = Lambda_3^sorted - Lambda_2^sorted`: `{result.min_gap_23:.12g}` at beta `{result.min_gap_beta:.12g} deg`.",
        f"- Is `gap_23` positive over the refined interval? `{'yes' if result.gap_23_positive else 'no'}`.",
        f"- Distinct refined tracking-warning rows for descendants 2/3: `{int(result.tracking_warning_count)}`.",
        f"- Classification: `{result.classification}`.",
        (
            "- Consistent with avoided crossing / mode veering? `yes, as a diagnostic candidate`: "
            "the descendant labels exchange sorted positions while the sorted frequency gap remains positive."
            if result.branch2_changes_sorted_index and result.gap_23_positive
            else (
                "- Consistent with avoided crossing / mode veering? `unresolved diagnostic candidate`: "
                "the sorted gap is positive and small, but refined branch tracking has warnings near the minimum."
                if result.gap_23_positive and result.tracking_warning_count > 0
                else "- Consistent with avoided crossing / mode veering? `not indicated by this check`."
            )
        ),
        "- Does this contradict the previous no-crossing statement? `no`; the previous statement concerns sorted-frequency crossings, while this audit treats descendant sorted-index changes as shape-label diagnostics.",
        f"- Required caution: {REQUIRED_STATEMENT}",
        "",
        "## Overview Branch-2 Sorted Position",
        "",
        *overview_lines,
        "",
        "## Refined Tracking Warnings",
        "",
        *(warning_lines if warning_lines else ["- none"]),
        "",
        "## Outputs",
        "",
        f"- CSV: `{display_path(result.csv_path)}`",
        f"- sorted Lambda plot: `{display_path(result.sorted_plot_path)}`",
        f"- gap plot: `{display_path(result.gap_plot_path)}`",
        f"- tracking map: `{display_path(result.tracking_map_path)}`",
        f"- MAC plot: `{display_path(result.mac_plot_path)}`",
        "",
        "## Protected Areas",
        "",
        "No article workspace files, `main.tex`, article figures, FEM/Gmsh/CalculiX files, old determinants, old solvers, or baseline results were modified.",
        "",
    ]
    result.report_path.parent.mkdir(parents=True, exist_ok=True)
    result.report_path.write_text("\n".join(lines), encoding="utf-8")


def display_path(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def run_audit(args: argparse.Namespace) -> AuditResult:
    paths = output_paths(args)
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)

    tracking_grid = refined_tracking_grid(args)
    report_grid = refined_report_grid(args, tracking_grid)
    base_tracking_args = tracking_namespace(args, branch=2, beta_step=float(args.refined_beta_step))
    sorted_roots_tracking = shape_scan.solve_sorted_root_grid(base_tracking_args, tracking_grid)
    branch_rows = {
        int(branch): shape_scan.track_descendants_by_beta(
            tracking_namespace(args, branch=int(branch), beta_step=float(args.refined_beta_step)),
            tracking_grid,
            sorted_roots_tracking,
            s_norm,
        )
        for branch in TARGET_BRANCHES
    }
    branch_diagnostics_all = build_branch_diagnostics(
        args,
        beta_values=tracking_grid,
        sorted_roots=sorted_roots_tracking,
        branch_rows=branch_rows,
        s_norm=s_norm,
    )

    refined_indices = [
        int(np.flatnonzero(np.isclose(tracking_grid, beta, rtol=0.0, atol=1e-12))[0])
        for beta in report_grid
    ]
    refined_sorted_roots = sorted_roots_tracking[np.asarray(refined_indices, dtype=int)]
    refined_branch_diagnostics = {
        branch: {
            round(float(beta), 10): branch_diagnostics_all[branch][round(float(beta), 10)]
            for beta in report_grid
        }
        for branch in TARGET_BRANCHES
    }

    write_csv(
        args,
        beta_values=report_grid,
        sorted_roots=refined_sorted_roots,
        branch_diagnostics=refined_branch_diagnostics,
        path=paths["csv"],
    )

    gap_23 = refined_sorted_roots[:, 2] - refined_sorted_roots[:, 1]
    min_idx = int(np.nanargmin(gap_23))
    min_gap = float(gap_23[min_idx])
    min_beta = float(report_grid[min_idx])
    gap_positive = bool(np.all(gap_23 > 0.0))
    branch2_positions = branch_positions(refined_branch_diagnostics, report_grid, 2)
    branch2_change = first_change_beta(report_grid, branch2_positions, base_position=2)
    branch2_changes = branch2_change is not None
    tracking_warning_count = sum(
        1
        for branch in TARGET_BRANCHES
        for beta in report_grid
        if refined_branch_diagnostics[int(branch)][round(float(beta), 10)].tracking_warning
    )
    classification = classify(
        branch2_changes=branch2_changes,
        gap_positive=gap_positive,
        min_gap=min_gap,
        near_zero_tol=float(args.near_zero_gap_tol),
        tracking_warning_count=tracking_warning_count,
    )

    plot_sorted_lambda(
        args,
        beta_values=report_grid,
        sorted_roots=refined_sorted_roots,
        branch_diagnostics=refined_branch_diagnostics,
        output=paths["sorted_plot"],
    )
    plot_gap(args, beta_values=report_grid, gap_23=gap_23, min_idx=min_idx, output=paths["gap_plot"])
    plot_tracking_map(beta_values=report_grid, branch_diagnostics=refined_branch_diagnostics, output=paths["tracking_map"])
    plot_mac(beta_values=report_grid, branch_diagnostics=refined_branch_diagnostics, output=paths["mac_plot"])

    overview_grid = shape_scan.inclusive_grid(
        float(args.overview_beta_start),
        float(args.overview_beta_end),
        float(args.overview_beta_step),
    )
    overview_args = tracking_namespace(args, branch=2, beta_step=float(args.overview_beta_step))
    overview_roots = shape_scan.solve_sorted_root_grid(overview_args, overview_grid)
    overview_branch2 = shape_scan.track_descendants_by_beta(overview_args, overview_grid, overview_roots, s_norm)

    result = AuditResult(
        beta_values=report_grid,
        sorted_roots=refined_sorted_roots,
        branch_rows=branch_rows,
        branch_diagnostics=refined_branch_diagnostics,
        csv_path=paths["csv"],
        sorted_plot_path=paths["sorted_plot"],
        gap_plot_path=paths["gap_plot"],
        tracking_map_path=paths["tracking_map"],
        mac_plot_path=paths["mac_plot"],
        report_path=paths["report"],
        min_gap_23=min_gap,
        min_gap_beta=min_beta,
        branch2_changes_sorted_index=branch2_changes,
        branch2_change_beta=branch2_change,
        gap_23_positive=gap_positive,
        tracking_warning_count=tracking_warning_count,
        classification=classification,
    )
    write_report(args, result=result, overview_beta_values=overview_grid, overview_branch2_rows=overview_branch2)
    return result


def main(argv: Sequence[str] | None = None) -> AuditResult:
    args = parse_args(argv)
    result = run_audit(args)
    print("mode-shape branch 2 transition audit:")
    print(f"CSV: {display_path(result.csv_path)}")
    print(f"sorted Lambda plot: {display_path(result.sorted_plot_path)}")
    print(f"gap plot: {display_path(result.gap_plot_path)}")
    print(f"tracking map: {display_path(result.tracking_map_path)}")
    print(f"MAC plot: {display_path(result.mac_plot_path)}")
    print(f"report: {display_path(result.report_path)}")
    print(f"branch 2 changes sorted index: {'yes' if result.branch2_changes_sorted_index else 'no'}")
    if result.branch2_change_beta is not None:
        print(f"branch 2 first changed sorted index at beta={result.branch2_change_beta:.12g} deg")
    print(f"min gap_23={result.min_gap_23:.12g} at beta={result.min_gap_beta:.12g} deg")
    print(f"gap_23 positive over refined interval: {'yes' if result.gap_23_positive else 'no'}")
    print(f"tracking warning rows over refined interval: {int(result.tracking_warning_count)}")
    print(f"classification: {result.classification}")
    return result


if __name__ == "__main__":
    main()
