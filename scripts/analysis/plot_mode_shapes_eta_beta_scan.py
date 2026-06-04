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
from matplotlib.lines import Line2D
import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    find_first_n_roots_eta,
    thickness_to_length_ratios,
)
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    analytic_null_vector,
    normalize_components,
)
from scripts.lib.thickness_mismatch_mac_tracking import (  # noqa: E402
    analytic_shape_vectors_for_roots,
    mac_assignment,
    reconstruct_components_eta,
    unique_nearest_frequency_assignment,
)


DEFAULT_ETA = 0.1
DEFAULT_MU = 0.3
DEFAULT_EPSILON = 0.0025
DEFAULT_BETA_START = 0.0
DEFAULT_BETA_END = 11.0
DEFAULT_BETA_STEP = 1.0
DEFAULT_BRANCH = 2
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results"
SMOKE_OUTPUT_DIR = DEFAULT_OUTPUT_DIR / "_smoke"

DEFAULT_NUM_TRACKED_BRANCHES = 8
DEFAULT_NUM_SORTED_ROOTS = 12
DEFAULT_ROOT_SCAN_STEP = 0.01
DEFAULT_ROOT_LMAX0 = 35.0
DEFAULT_NUM_SHAPE_SAMPLES = 401
DEFAULT_MAC_WARNING_THRESHOLD = 0.9
DEFAULT_MAC_MARGIN_WARNING_THRESHOLD = 0.05
DEFAULT_MAX_SORTED_POSITION_JUMP = 1
DEFAULT_L_TOTAL = 2.0
DEFAULT_MODE_SCALE = 0.14
DEFAULT_NORMALIZE = "max-full"
DEFAULT_DPI = 240
DEFAULT_GRID_COLUMNS = 4

COEFF_NAMES = ("A1", "B1", "A2", "B2", "P1", "P2")
COMPONENT_KEYS = ("u_left", "w_left", "u_right", "w_right")
SUMMARY_FIELDNAMES = [
    "eta",
    "mu",
    "epsilon",
    "beta_deg",
    "descendant_branch",
    "branch_seed",
    "sorted_index",
    "Lambda",
    "coeff_A1",
    "coeff_B1",
    "coeff_A2",
    "coeff_B2",
    "coeff_P1",
    "coeff_P2",
    "mac_to_previous",
    "accepted_mac_to_previous",
    "second_best_mac",
    "mac_margin",
    "tracking_status",
    "tracking_warning",
    "matrix_residual_max",
    "smallest_singular_value",
    "singular_value_ratio",
    "thickness_ratio_1",
    "thickness_ratio_2",
    "output_png",
    "warning",
    "notes",
]


@dataclass(frozen=True)
class TrackingRow:
    beta_deg: float
    descendant_branch: int
    Lambda: float
    sorted_index: int
    tracking_status: str
    mac_to_previous: float
    accepted_mac_to_previous: float
    second_best_mac: float
    mac_margin: float
    frequency_assignment_margin: float
    tracking_warning: str


@dataclass(frozen=True)
class ShapeCase:
    beta_deg: float
    descendant_branch: int
    sorted_index: int
    Lambda: float
    coeff: np.ndarray
    components: dict[str, np.ndarray]
    tracking_status: str
    mac_to_previous: float
    accepted_mac_to_previous: float
    second_best_mac: float
    mac_margin: float
    tracking_warning: str
    matrix_residual_max: float
    smallest_singular_value: float
    singular_value_ratio: float
    thickness_ratio_1: float
    thickness_ratio_2: float
    output_png: Path
    warning: str
    notes: str


def number_token(value: float) -> str:
    return f"{float(value):.10g}".replace("-", "m").replace("+", "").replace(".", "p")


def number_text(value: float | int | str) -> str:
    if isinstance(value, str):
        return value
    value_f = float(value)
    return "" if not np.isfinite(value_f) else f"{value_f:.12g}"


def beta_label(value: float) -> str:
    value_f = float(value)
    rounded = round(value_f)
    if np.isclose(value_f, rounded, rtol=0.0, atol=1e-12):
        return str(int(rounded))
    return f"{value_f:g}"


def beta_token(value: float) -> str:
    return number_token(float(value))


def output_stem(args: argparse.Namespace) -> str:
    return (
        f"mode_shape_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_"
        f"eps_{number_token(float(args.epsilon))}_"
        f"branch{int(args.branch)}"
    )


def output_png_path(args: argparse.Namespace, beta_deg: float) -> Path:
    return Path(args.output_dir) / f"{output_stem(args)}_beta_{beta_token(float(beta_deg))}deg.png"


def output_grid_path(args: argparse.Namespace, beta_values: Sequence[float]) -> Path:
    return (
        Path(args.output_dir)
        / f"mode_shapes_eta_{number_token(float(args.eta))}_"
        f"mu_{number_token(float(args.mu))}_"
        f"eps_{number_token(float(args.epsilon))}_"
        f"branch{int(args.branch)}_"
        f"beta_{beta_token(float(beta_values[0]))}_to_{beta_token(float(beta_values[-1]))}deg_grid.png"
    )


def output_summary_path(args: argparse.Namespace) -> Path:
    return Path(args.output_dir) / f"{output_stem(args)}_summary.csv"


def resolve_repo_path(path: str | Path) -> Path:
    value = Path(path)
    if value.is_absolute():
        return value
    return REPO_ROOT / value


def inclusive_grid(start: float, stop: float, step: float) -> np.ndarray:
    start_f = float(start)
    stop_f = float(stop)
    step_f = float(step)
    if step_f <= 0.0:
        raise ValueError("beta step must be positive.")
    if stop_f < start_f:
        raise ValueError("beta end must be greater than or equal to beta start.")
    values = np.arange(start_f, stop_f + 0.5 * step_f, step_f, dtype=float)
    values = values[values <= stop_f + 1e-12]
    if values.size == 0 or not np.isclose(values[0], start_f, rtol=0.0, atol=1e-12):
        values = np.insert(values, 0, start_f)
    if not np.isclose(values[-1], stop_f, rtol=0.0, atol=1e-12):
        values = np.append(values, stop_f)
    values[0] = start_f
    values[-1] = stop_f
    return np.unique(np.round(values, 12))


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot diagnostic analytic Euler-Bernoulli thickness-mismatch mode shapes "
            "for one descendant branch over a beta scan."
        ),
    )
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--beta-start", type=float, default=DEFAULT_BETA_START)
    parser.add_argument("--beta-end", type=float, default=DEFAULT_BETA_END)
    parser.add_argument("--beta-step", type=float, default=DEFAULT_BETA_STEP)
    parser.add_argument("--branch", type=int, default=DEFAULT_BRANCH)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--num-tracked-branches", type=int, default=DEFAULT_NUM_TRACKED_BRANCHES)
    parser.add_argument("--num-sorted-roots", type=int, default=DEFAULT_NUM_SORTED_ROOTS)
    parser.add_argument("--root-scan-step", type=float, default=DEFAULT_ROOT_SCAN_STEP)
    parser.add_argument("--root-lmax0", type=float, default=DEFAULT_ROOT_LMAX0)
    parser.add_argument("--num-shape-samples", type=int, default=DEFAULT_NUM_SHAPE_SAMPLES)
    parser.add_argument("--mac-warning-threshold", type=float, default=DEFAULT_MAC_WARNING_THRESHOLD)
    parser.add_argument("--mac-margin-warning-threshold", type=float, default=DEFAULT_MAC_MARGIN_WARNING_THRESHOLD)
    parser.add_argument("--max-sorted-position-jump", type=int, default=DEFAULT_MAX_SORTED_POSITION_JUMP)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--mode-scale", type=float, default=DEFAULT_MODE_SCALE)
    parser.add_argument("--normalize", choices=("max-full", "max-transverse", "none"), default=DEFAULT_NORMALIZE)
    parser.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    parser.add_argument("--grid-columns", type=int, default=DEFAULT_GRID_COLUMNS)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Use a tiny beta grid for wiring checks and write default outputs under results/_smoke/.",
    )
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if bool(args.smoke):
        apply_smoke_defaults(args)
    validate_args(args)
    args.output_dir = resolve_repo_path(args.output_dir)
    return args


def apply_smoke_defaults(args: argparse.Namespace) -> None:
    args.beta_start = 0.0
    args.beta_end = 5.0
    args.beta_step = 5.0
    args.num_shape_samples = min(int(args.num_shape_samples), 81)
    args.num_tracked_branches = max(int(args.branch), min(int(args.num_tracked_branches), 4))
    args.num_sorted_roots = max(int(args.num_tracked_branches), min(int(args.num_sorted_roots), 5))
    if Path(args.output_dir) == DEFAULT_OUTPUT_DIR:
        args.output_dir = SMOKE_OUTPUT_DIR


def validate_args(args: argparse.Namespace) -> None:
    if int(args.branch) <= 0:
        raise ValueError("--branch must be positive.")
    if int(args.num_tracked_branches) < int(args.branch):
        raise ValueError("--num-tracked-branches must be at least --branch.")
    if int(args.num_sorted_roots) < int(args.num_tracked_branches):
        raise ValueError("--num-sorted-roots must be at least --num-tracked-branches.")
    if float(args.epsilon) <= 0.0:
        raise ValueError("--epsilon must be positive.")
    if not (-1.0 < float(args.mu) < 1.0):
        raise ValueError("--mu must lie inside (-1, 1).")
    if not (-1.0 < float(args.eta) < 1.0):
        raise ValueError("--eta must lie inside (-1, 1).")
    if float(args.beta_start) < -1e-12:
        raise ValueError("--beta-start must be non-negative.")
    if float(args.beta_step) <= 0.0:
        raise ValueError("--beta-step must be positive.")
    if int(args.num_shape_samples) < 5:
        raise ValueError("--num-shape-samples must be at least 5.")
    if float(args.root_scan_step) <= 0.0:
        raise ValueError("--root-scan-step must be positive.")
    if float(args.root_lmax0) <= 0.0:
        raise ValueError("--root-lmax0 must be positive.")
    if float(args.l_total) <= 0.0:
        raise ValueError("--l-total must be positive.")
    if float(args.mode_scale) <= 0.0:
        raise ValueError("--mode-scale must be positive.")
    if int(args.dpi) <= 0:
        raise ValueError("--dpi must be positive.")
    if int(args.grid_columns) <= 0:
        raise ValueError("--grid-columns must be positive.")


def yesno(value: bool) -> str:
    return "yes" if bool(value) else "no"


def tracking_warning(row: dict[str, float | int | str]) -> str:
    flags: list[str] = []
    for key, label in (
        ("low_mac", "low_mac"),
        ("low_margin", "low_margin"),
        ("large_beta_step", "large_beta_step"),
        ("blocked_by_unresolved_neighbor", "blocked_by_unresolved_neighbor"),
        ("unresolved_assignment", "unresolved_assignment"),
        ("frequency_mac_disagreement", "frequency_mac_disagreement"),
    ):
        if str(row.get(key, "no")) == "yes":
            flags.append(label)
    status = str(row.get("tracking_step_status", ""))
    if status not in {"", "seed_beta0", "mac_ok"} and status not in flags:
        flags.append(status)
    return ";".join(flags)


def solve_sorted_root_grid(args: argparse.Namespace, beta_values: np.ndarray) -> np.ndarray:
    roots = np.full((len(beta_values), int(args.num_sorted_roots)), np.nan, dtype=float)
    for idx, beta_deg in enumerate(beta_values):
        values = find_first_n_roots_eta(
            float(np.deg2rad(float(beta_deg))),
            float(args.mu),
            float(args.epsilon),
            float(args.eta),
            int(args.num_sorted_roots),
            Lmax0=float(args.root_lmax0),
            scan_step=float(args.root_scan_step),
        )
        values = np.asarray(values, dtype=float)
        if values.shape[0] < int(args.num_sorted_roots):
            padded = np.full(int(args.num_sorted_roots), np.nan, dtype=float)
            padded[: values.shape[0]] = values
            values = padded
        if np.any(~np.isfinite(values[: int(args.num_tracked_branches)])):
            raise RuntimeError(
                f"Missing roots needed for branch tracking at beta={float(beta_deg):g} deg."
            )
        roots[idx] = values[: int(args.num_sorted_roots)]
    return roots


def shape_vectors_for_beta(
    args: argparse.Namespace,
    roots: np.ndarray,
    *,
    beta_deg: float,
    s_norm: np.ndarray,
) -> list[np.ndarray]:
    finite_roots = np.asarray(roots, dtype=float)
    if np.any(~np.isfinite(finite_roots)):
        raise RuntimeError(f"Non-finite root found at beta={float(beta_deg):g} deg.")
    return analytic_shape_vectors_for_roots(
        finite_roots,
        beta_rad=float(np.deg2rad(float(beta_deg))),
        mu=float(args.mu),
        epsilon=float(args.epsilon),
        eta=float(args.eta),
        s_norm=s_norm,
    )


def build_seed_rows(args: argparse.Namespace, beta_deg: float, lambdas: np.ndarray) -> list[dict[str, float | int | str]]:
    rows: list[dict[str, float | int | str]] = []
    for branch_idx, value in enumerate(lambdas, start=1):
        rows.append(
            {
                "beta_prev_deg": np.nan,
                "beta_deg": float(beta_deg),
                "descendant_branch": int(branch_idx),
                "Lambda_tracked": float(value),
                "current_sorted_position": int(branch_idx),
                "diagnostic_candidate_sorted_position": int(branch_idx),
                "diagnostic_candidate_Lambda": float(value),
                "nearest_sorted_root_index": int(branch_idx),
                "tracking_step_status": "seed_beta0",
                "mac_to_previous": np.nan,
                "accepted_mac_to_previous": np.nan,
                "second_best_mac": np.nan,
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
                "suspicious_assignment": "no",
                "requires_refined_check": "no",
            }
        )
    return rows


def track_descendants_by_beta(
    args: argparse.Namespace,
    beta_values: np.ndarray,
    sorted_roots: np.ndarray,
    s_norm: np.ndarray,
) -> list[TrackingRow]:
    n_track = int(args.num_tracked_branches)
    if not np.isclose(float(beta_values[0]), 0.0, rtol=0.0, atol=1e-12):
        raise ValueError("Beta descendant tracking must start at beta=0.")
    if sorted_roots.shape[1] < n_track:
        raise ValueError("sorted root grid does not contain enough roots.")

    initial_vectors = shape_vectors_for_beta(args, sorted_roots[0], beta_deg=float(beta_values[0]), s_norm=s_norm)
    previous_vectors = initial_vectors[:n_track]
    previous_lambdas = np.asarray(sorted_roots[0, :n_track], dtype=float)
    previous_positions = np.arange(1, n_track + 1, dtype=int)
    max_beta_step = float(args.beta_step)

    all_rows = build_seed_rows(args, float(beta_values[0]), previous_lambdas)

    for col in range(1, len(beta_values)):
        beta_prev = float(beta_values[col - 1])
        beta = float(beta_values[col])
        roots = np.asarray(sorted_roots[col], dtype=float)
        candidate_vectors = shape_vectors_for_beta(args, roots, beta_deg=beta, s_norm=s_norm)
        raw_cols, mac = mac_assignment(previous_vectors, candidate_vectors)
        freq_assignments = unique_nearest_frequency_assignment(previous_lambdas, roots)

        second_best_values = np.full(n_track, np.nan, dtype=float)
        candidate_macs = np.full(n_track, np.nan, dtype=float)
        accepted_macs = np.full(n_track, np.nan, dtype=float)
        mac_margins = np.full(n_track, np.nan, dtype=float)
        candidate_jumps = np.zeros(n_track, dtype=int)
        low_mac_flags = np.zeros(n_track, dtype=bool)
        low_margin_flags = np.zeros(n_track, dtype=bool)
        large_step_flags = np.zeros(n_track, dtype=bool)
        unresolved = np.zeros(n_track, dtype=bool)

        for branch_row in range(n_track):
            candidate_col = int(raw_cols[branch_row])
            row_macs = np.sort(mac[branch_row])[::-1]
            second_best = float(row_macs[1]) if len(row_macs) > 1 else np.nan
            candidate_mac = float(mac[branch_row, candidate_col])
            mac_margin = candidate_mac - second_best if np.isfinite(second_best) else np.nan
            candidate_jump = int(candidate_col) + 1 - int(previous_positions[branch_row])
            low_mac = candidate_mac < float(args.mac_warning_threshold)
            low_margin = np.isfinite(mac_margin) and mac_margin < float(args.mac_margin_warning_threshold)
            large_step = abs(beta - beta_prev) > max_beta_step + 1e-12

            second_best_values[branch_row] = second_best
            candidate_macs[branch_row] = candidate_mac
            mac_margins[branch_row] = mac_margin
            candidate_jumps[branch_row] = int(candidate_jump)
            low_mac_flags[branch_row] = bool(low_mac)
            low_margin_flags[branch_row] = bool(low_margin)
            large_step_flags[branch_row] = bool(large_step)
            unresolved[branch_row] = (
                abs(candidate_jump) > int(args.max_sorted_position_jump)
                or low_mac
                or low_margin
                or large_step
            )

        blocked_by_unresolved_neighbor = np.zeros(n_track, dtype=bool)
        while True:
            retained_positions = {
                int(previous_positions[row])
                for row in range(n_track)
                if bool(unresolved[row])
            }
            updated = unresolved.copy()
            for branch_row in range(n_track):
                if bool(updated[branch_row]):
                    continue
                candidate_position = int(raw_cols[branch_row]) + 1
                if candidate_position in retained_positions and candidate_position != int(previous_positions[branch_row]):
                    updated[branch_row] = True
                    blocked_by_unresolved_neighbor[branch_row] = True
            if np.array_equal(updated, unresolved):
                break
            unresolved = updated

        accepted_cols = np.array(
            [
                int(previous_positions[row]) - 1 if bool(unresolved[row]) else int(raw_cols[row])
                for row in range(n_track)
            ],
            dtype=int,
        )

        current_lambdas = np.full(n_track, np.nan, dtype=float)
        current_positions = np.full(n_track, -1, dtype=int)
        current_vectors: list[np.ndarray] = []

        for branch_row in range(n_track):
            candidate_col = int(raw_cols[branch_row])
            root_col = int(accepted_cols[branch_row])
            if root_col < 0 or root_col >= len(roots):
                raise RuntimeError(
                    "Cannot retain previous canonical sorted position "
                    f"{int(previous_positions[branch_row])} at beta={beta:g} deg."
                )
            accepted_mac = float(mac[branch_row, root_col])
            accepted_macs[branch_row] = accepted_mac
            frequency = freq_assignments[branch_row]
            disagreement = int(frequency.root_index) != int(candidate_col) + 1
            status = "unresolved_assignment" if bool(unresolved[branch_row]) else "mac_ok"
            row = {
                "beta_prev_deg": beta_prev,
                "beta_deg": beta,
                "descendant_branch": int(branch_row) + 1,
                "Lambda_tracked": float(roots[root_col]),
                "current_sorted_position": int(root_col) + 1,
                "diagnostic_candidate_sorted_position": int(candidate_col) + 1,
                "diagnostic_candidate_Lambda": float(roots[candidate_col]),
                "nearest_sorted_root_index": int(frequency.root_index),
                "nearest_sorted_Lambda": float(roots[int(frequency.root_index) - 1]),
                "tracking_step_status": status,
                "mac_to_previous": float(candidate_macs[branch_row]),
                "accepted_mac_to_previous": accepted_mac,
                "second_best_mac": float(second_best_values[branch_row]),
                "frequency_distance_from_previous": float(frequency.distance_from_previous),
                "frequency_second_nearest_distance": float(frequency.second_nearest_distance),
                "frequency_assignment_margin": float(frequency.assignment_margin),
                "frequency_mac_disagreement": yesno(disagreement),
                "mac_margin": float(mac_margins[branch_row]),
                "assigned_from_previous_sorted_position": int(previous_positions[branch_row]),
                "sorted_position_jump": int(root_col) + 1 - int(previous_positions[branch_row]),
                "diagnostic_candidate_sorted_position_jump": int(candidate_jumps[branch_row]),
                "low_mac": yesno(bool(low_mac_flags[branch_row])),
                "low_margin": yesno(bool(low_margin_flags[branch_row])),
                "large_beta_step": yesno(bool(large_step_flags[branch_row])),
                "blocked_by_unresolved_neighbor": yesno(bool(blocked_by_unresolved_neighbor[branch_row])),
                "unresolved_assignment": yesno(bool(unresolved[branch_row])),
                "suspicious_assignment": yesno(bool(unresolved[branch_row])),
                "requires_refined_check": yesno(bool(unresolved[branch_row])),
            }
            all_rows.append(row)
            current_lambdas[branch_row] = float(roots[root_col])
            current_positions[branch_row] = int(root_col) + 1
            current_vectors.append(candidate_vectors[root_col])

        previous_lambdas = current_lambdas
        previous_positions = current_positions
        previous_vectors = current_vectors

    target_rows: list[TrackingRow] = []
    for beta in beta_values:
        matches = [
            row
            for row in all_rows
            if int(row["descendant_branch"]) == int(args.branch)
            and np.isclose(float(row["beta_deg"]), float(beta), rtol=0.0, atol=1e-12)
        ]
        if len(matches) != 1:
            raise RuntimeError(f"Could not locate tracked branch {int(args.branch)} at beta={float(beta):g} deg.")
        row = matches[0]
        target_rows.append(
            TrackingRow(
                beta_deg=float(row["beta_deg"]),
                descendant_branch=int(row["descendant_branch"]),
                Lambda=float(row["Lambda_tracked"]),
                sorted_index=int(row["current_sorted_position"]),
                tracking_status=str(row["tracking_step_status"]),
                mac_to_previous=float(row["mac_to_previous"]),
                accepted_mac_to_previous=float(row["accepted_mac_to_previous"]),
                second_best_mac=float(row["second_best_mac"]),
                mac_margin=float(row["mac_margin"]),
                frequency_assignment_margin=float(row["frequency_assignment_margin"]),
                tracking_warning=tracking_warning(row),
            )
        )
    return target_rows


def component_vector(components: dict[str, np.ndarray]) -> np.ndarray:
    return np.concatenate([np.asarray(components[key], dtype=float) for key in COMPONENT_KEYS])


def build_shape_cases(
    args: argparse.Namespace,
    tracking_rows: Sequence[TrackingRow],
    s_norm: np.ndarray,
) -> list[ShapeCase]:
    cases: list[ShapeCase] = []
    previous_components: dict[str, np.ndarray] | None = None
    for row in tracking_rows:
        matrix = assemble_clamped_coupled_matrix_eta(
            row.Lambda,
            float(np.deg2rad(float(row.beta_deg))),
            float(args.mu),
            float(args.epsilon),
            float(args.eta),
        )
        coeff, smallest, ratio = analytic_null_vector(matrix)
        if not np.all(np.isfinite(coeff)):
            raise RuntimeError(f"Non-finite null-vector coefficient at beta={row.beta_deg:g} deg.")
        components_raw = reconstruct_components_eta(
            row.Lambda,
            mu=float(args.mu),
            eta=float(args.eta),
            epsilon=float(args.epsilon),
            coeff=coeff,
            s_norm=s_norm,
        )
        components, _scale = normalize_components(
            components_raw,
            plot_kind="full",
            normalize=str(args.normalize),
        )
        if previous_components is not None:
            if float(np.dot(component_vector(components), component_vector(previous_components))) < 0.0:
                components = {key: -np.asarray(value, dtype=float) for key, value in components.items()}
                coeff = -coeff
        previous_components = components

        matrix_residual = np.asarray(matrix, dtype=float) @ np.asarray(coeff, dtype=float)
        residual_max = float(np.max(np.abs(matrix_residual)))
        ratio1, ratio2 = thickness_to_length_ratios(float(args.epsilon), float(args.mu), float(args.eta))

        warning_parts = []
        if row.tracking_warning:
            warning_parts.append(row.tracking_warning)
        if not np.isfinite(row.Lambda):
            warning_parts.append("nonfinite_lambda")
        if residual_max > 1e-7:
            warning_parts.append("matrix_residual_gt_1e-7")
        if ratio1 > 0.1 or ratio2 > 0.1:
            warning_parts.append("thin_rod_ratio_gt_0p1")

        notes = (
            "eta-aware beta shape-MAC descendant tracking seeded at beta=0 for fixed mu and eta; "
            f"components normalized by {args.normalize}; coefficient sign aligned to previous beta"
        )
        cases.append(
            ShapeCase(
                beta_deg=float(row.beta_deg),
                descendant_branch=int(row.descendant_branch),
                sorted_index=int(row.sorted_index),
                Lambda=float(row.Lambda),
                coeff=np.asarray(coeff, dtype=float),
                components={key: np.asarray(value, dtype=float) for key, value in components.items()},
                tracking_status=row.tracking_status,
                mac_to_previous=float(row.mac_to_previous),
                accepted_mac_to_previous=float(row.accepted_mac_to_previous),
                second_best_mac=float(row.second_best_mac),
                mac_margin=float(row.mac_margin),
                tracking_warning=row.tracking_warning,
                matrix_residual_max=residual_max,
                smallest_singular_value=float(smallest),
                singular_value_ratio=float(ratio),
                thickness_ratio_1=float(ratio1),
                thickness_ratio_2=float(ratio2),
                output_png=output_png_path(args, float(row.beta_deg)),
                warning=";".join(warning_parts) if warning_parts else "no",
                notes=notes,
            )
        )
    return cases


def arm_lengths(args: argparse.Namespace) -> tuple[float, float]:
    ell = float(args.l_total) / 2.0
    return ell * (1.0 - float(args.mu)), ell * (1.0 + float(args.mu))


def base_coordinates(
    args: argparse.Namespace,
    s_norm: np.ndarray,
    *,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    left_length, right_length = arm_lengths(args)
    beta_rad = float(np.deg2rad(float(beta_deg)))
    xi = np.asarray(s_norm, dtype=float)
    x_left = left_length * xi
    y_left = np.zeros_like(x_left)
    x_right = left_length + right_length * xi * np.cos(beta_rad)
    y_right = right_length * xi * np.sin(beta_rad)
    return x_left, y_left, x_right, y_right


def deformed_coordinates(
    args: argparse.Namespace,
    components: dict[str, np.ndarray],
    s_norm: np.ndarray,
    *,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x_left, y_left, x_right, y_right = base_coordinates(args, s_norm, beta_deg=beta_deg)
    beta_rad = float(np.deg2rad(float(beta_deg)))
    tangent_right = np.array([np.cos(beta_rad), np.sin(beta_rad)], dtype=float)
    normal_right = np.array([-np.sin(beta_rad), np.cos(beta_rad)], dtype=float)
    u_left = np.asarray(components["u_left"], dtype=float)
    w_left = np.asarray(components["w_left"], dtype=float)
    u_right = np.asarray(components["u_right"], dtype=float)
    w_right = np.asarray(components["w_right"], dtype=float)
    x_left_def = x_left + float(args.mode_scale) * u_left
    y_left_def = y_left + float(args.mode_scale) * w_left
    x_right_def = x_right + float(args.mode_scale) * (u_right * tangent_right[0] + w_right * normal_right[0])
    y_right_def = y_right + float(args.mode_scale) * (u_right * tangent_right[1] + w_right * normal_right[1])
    return x_left_def, y_left_def, x_right_def, y_right_def


def shared_axis_limits(
    args: argparse.Namespace,
    cases: Sequence[ShapeCase],
    s_norm: np.ndarray,
) -> tuple[tuple[float, float], tuple[float, float]]:
    x_values: list[np.ndarray] = []
    y_values: list[np.ndarray] = []
    for case in cases:
        x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(args, s_norm, beta_deg=case.beta_deg)
        x_left_def, y_left_def, x_right_def, y_right_def = deformed_coordinates(
            args,
            case.components,
            s_norm,
            beta_deg=case.beta_deg,
        )
        x_values.extend([x_left_base, x_right_base, x_left_def, x_right_def])
        y_values.extend([y_left_base, y_right_base, y_left_def, y_right_def])
    x_all = np.concatenate([np.asarray(value, dtype=float) for value in x_values])
    y_all = np.concatenate([np.asarray(value, dtype=float) for value in y_values])
    x_span = max(float(np.max(x_all) - np.min(x_all)), 1.0)
    y_span = max(float(np.max(y_all) - np.min(y_all)), 0.4)
    return (
        (float(np.min(x_all) - 0.06 * x_span), float(np.max(x_all) + 0.06 * x_span)),
        (float(np.min(y_all) - 0.18 * y_span), float(np.max(y_all) + 0.18 * y_span)),
    )


def draw_shape_case(
    ax: plt.Axes,
    args: argparse.Namespace,
    case: ShapeCase,
    s_norm: np.ndarray,
    *,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
    include_labels: bool,
) -> None:
    x_left_base, y_left_base, x_right_base, y_right_base = base_coordinates(args, s_norm, beta_deg=case.beta_deg)
    x_left_def, y_left_def, x_right_def, y_right_def = deformed_coordinates(
        args,
        case.components,
        s_norm,
        beta_deg=case.beta_deg,
    )
    undeformed_label = "undeformed" if include_labels else None
    deformed_label = "deformed descendant" if include_labels else None
    joint_label = "joint" if include_labels else None
    ax.plot(x_left_base, y_left_base, color="0.70", linestyle="--", linewidth=1.0, label=undeformed_label)
    ax.plot(x_right_base, y_right_base, color="0.70", linestyle="--", linewidth=1.0)
    ax.plot(x_left_def, y_left_def, color="#006BA4", linewidth=2.0, label=deformed_label)
    ax.plot(x_right_def, y_right_def, color="#006BA4", linewidth=2.0)
    ax.scatter([x_left_base[-1]], [y_left_base[-1]], color="black", s=18, zorder=5, label=joint_label)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(*axis_limits[0])
    ax.set_ylim(*axis_limits[1])
    ax.grid(True, color="0.88", linewidth=0.6)


def single_title(args: argparse.Namespace, case: ShapeCase) -> str:
    return (
        "Euler-Bernoulli mode shape, "
        f"eta={float(args.eta):g}, mu={float(args.mu):g}, epsilon={float(args.epsilon):g}, "
        f"beta={case.beta_deg:g} deg, descendant branch {case.descendant_branch}, "
        f"sorted index {case.sorted_index}, Lambda={case.Lambda:.8g}"
    )


def plot_single_cases(
    args: argparse.Namespace,
    cases: Sequence[ShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
) -> list[Path]:
    output_paths: list[Path] = []
    for case in cases:
        output = case.output_png
        output.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(7.8, 4.8))
        draw_shape_case(ax, args, case, s_norm, axis_limits=axis_limits, include_labels=True)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.legend(loc="best", fontsize=8, frameon=False)
        ax.set_title(single_title(args, case), fontsize=10)
        fig.tight_layout()
        fig.savefig(output, dpi=int(args.dpi), bbox_inches="tight")
        plt.close(fig)
        output_paths.append(output)
    return output_paths


def shape_legend_handles() -> list[Line2D]:
    return [
        Line2D([0], [0], color="0.70", linestyle="--", linewidth=1.0, label="undeformed"),
        Line2D([0], [0], color="#006BA4", linewidth=2.0, label="deformed descendant"),
        Line2D([0], [0], color="black", marker="o", linestyle="None", markersize=4, label="joint"),
    ]


def plot_grid(
    args: argparse.Namespace,
    cases: Sequence[ShapeCase],
    s_norm: np.ndarray,
    axis_limits: tuple[tuple[float, float], tuple[float, float]],
) -> Path:
    n_cases = len(cases)
    n_cols = min(int(args.grid_columns), n_cases)
    n_rows = int(np.ceil(n_cases / n_cols))
    output = output_grid_path(args, [case.beta_deg for case in cases])
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.8 * n_cols, 1.95 * n_rows + 0.6), squeeze=False)
    for idx, ax in enumerate(axes.ravel()):
        if idx >= n_cases:
            ax.axis("off")
            continue
        case = cases[idx]
        draw_shape_case(ax, args, case, s_norm, axis_limits=axis_limits, include_labels=False)
        ax.set_title(
            (
                f"beta={case.beta_deg:g} deg, Lambda={case.Lambda:.5g}\n"
                f"desc {case.descendant_branch}, sorted {case.sorted_index}"
            ),
            fontsize=9,
        )
        ax.set_xlabel("x", fontsize=8)
        ax.set_ylabel("y", fontsize=8)
        ax.tick_params(labelsize=7)

    fig.suptitle(
        (
            f"Diagnostic Euler-Bernoulli mode shapes: eta={float(args.eta):g}, "
            f"mu={float(args.mu):g}, epsilon={float(args.epsilon):g}, "
            f"descendant branch {int(args.branch)}"
        ),
        fontsize=13,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.93), h_pad=0.9, w_pad=0.6)
    fig.savefig(output, dpi=int(args.dpi), bbox_inches="tight")
    plt.close(fig)
    return output


def write_summary_csv(args: argparse.Namespace, cases: Sequence[ShapeCase]) -> Path:
    output = output_summary_path(args)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        for case in cases:
            coeff = {name: float(value) for name, value in zip(COEFF_NAMES, case.coeff, strict=True)}
            writer.writerow(
                {
                    "eta": number_text(float(args.eta)),
                    "mu": number_text(float(args.mu)),
                    "epsilon": number_text(float(args.epsilon)),
                    "beta_deg": number_text(case.beta_deg),
                    "descendant_branch": int(case.descendant_branch),
                    "branch_seed": "beta0_fixed_mu_eta",
                    "sorted_index": int(case.sorted_index),
                    "Lambda": number_text(case.Lambda),
                    "coeff_A1": number_text(coeff["A1"]),
                    "coeff_B1": number_text(coeff["B1"]),
                    "coeff_A2": number_text(coeff["A2"]),
                    "coeff_B2": number_text(coeff["B2"]),
                    "coeff_P1": number_text(coeff["P1"]),
                    "coeff_P2": number_text(coeff["P2"]),
                    "mac_to_previous": number_text(case.mac_to_previous),
                    "accepted_mac_to_previous": number_text(case.accepted_mac_to_previous),
                    "second_best_mac": number_text(case.second_best_mac),
                    "mac_margin": number_text(case.mac_margin),
                    "tracking_status": case.tracking_status,
                    "tracking_warning": case.tracking_warning if case.tracking_warning else "no",
                    "matrix_residual_max": number_text(case.matrix_residual_max),
                    "smallest_singular_value": number_text(case.smallest_singular_value),
                    "singular_value_ratio": number_text(case.singular_value_ratio),
                    "thickness_ratio_1": number_text(case.thickness_ratio_1),
                    "thickness_ratio_2": number_text(case.thickness_ratio_2),
                    "output_png": display_path(case.output_png),
                    "warning": case.warning,
                    "notes": case.notes,
                }
            )
    return output


def display_path(path: Path) -> str:
    try:
        return Path(path).resolve().relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return str(path)


def print_summary(cases: Sequence[ShapeCase], *, grid_path: Path, csv_path: Path) -> None:
    print("tracked descendant summary:")
    print("beta_deg, descendant_branch, sorted_index, Lambda, warning")
    for case in cases:
        print(
            f"{case.beta_deg:g}, {case.descendant_branch}, {case.sorted_index}, "
            f"{case.Lambda:.12g}, {case.warning}"
        )
    sorted_positions = [int(case.sorted_index) for case in cases]
    changed = len(set(sorted_positions)) > 1
    print(f"sorted label rearrangement over beta: {'yes' if changed else 'no'}")
    print(f"saved grid PNG: {display_path(grid_path)}")
    print(f"saved summary CSV: {display_path(csv_path)}")


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    beta_values = inclusive_grid(float(args.beta_start), float(args.beta_end), float(args.beta_step))
    s_norm = np.linspace(0.0, 1.0, int(args.num_shape_samples), dtype=float)

    sorted_roots = solve_sorted_root_grid(args, beta_values)
    tracking_rows = track_descendants_by_beta(args, beta_values, sorted_roots, s_norm)
    cases = build_shape_cases(args, tracking_rows, s_norm)
    axis_limits = shared_axis_limits(args, cases, s_norm)
    png_paths = plot_single_cases(args, cases, s_norm, axis_limits)
    grid_path = plot_grid(args, cases, s_norm, axis_limits)
    csv_path = write_summary_csv(args, cases)
    print_summary(cases, grid_path=grid_path, csv_path=csv_path)
    return {
        "cases": cases,
        "png_paths": png_paths,
        "grid_path": grid_path,
        "csv_path": csv_path,
    }


if __name__ == "__main__":
    main()
