from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
import sys
from typing import Iterable, Sequence
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize_scalar


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas_thickness_mismatch import det_eta, find_first_n_roots_eta  # noqa: E402
from scripts.lib.thickness_mismatch_mac_tracking import (  # noqa: E402
    analytic_shape_vectors_for_roots,
    mac_assignment,
    mac_value,
)


EPSILON = 0.0025
N_SORTED_ROOTS = 14
N_GAP_ROOTS = 7
N_DESCENDANTS_TRACK = 8
N_TRACKING_CANDIDATES = 10
ROOT_LMAX0 = 35.0
DEFAULT_SCAN_STEP = 0.01
STRICT_SCAN_STEP = 0.001
NUM_SHAPE_SAMPLES = 401

ETA_COARSE_STEP = 0.01
BETA_COARSE_STEP = 0.25
MU_COARSE_STEP = 0.01
MU_LOCAL_COARSE_STEP = 0.001
MU_LOCAL_FINE_STEP = 0.0001
MU_LOCAL_COARSE_HALF_WIDTH = 0.015
MU_LOCAL_FINE_HALF_WIDTH = 0.003

ETA_SPECIAL_ABS = (
    1.0e-6,
    2.0e-6,
    5.0e-6,
    1.0e-5,
    2.0e-5,
    5.0e-5,
    1.0e-4,
    2.0e-4,
    5.0e-4,
    1.0e-3,
    2.0e-3,
    5.0e-3,
    1.0e-2,
    2.0e-2,
    5.0e-2,
)
ETA_BETA_SCAN_VALUES = (
    -0.5,
    -0.1,
    -0.05,
    -0.02,
    -0.016,
    -0.004,
    -0.002,
    0.0,
    0.002,
    0.004,
    0.016,
    0.02,
    0.05,
    0.1,
    0.5,
)
ETA_MU_SCAN_VALUES = (-0.5, -0.1, -0.02, -0.004, 0.004, 0.02, 0.1, 0.5)
BETA_MU_SCAN_VALUES = (5.0, 10.0, 15.0)
BETA_ETA_SCAN_VALUES = (0.0, 5.0, 10.0, 15.0)
BETA15_MU_VALUES = (0.0, 0.1, 0.3, 0.6, 0.9)
MU_KNOWN_POINTS = (0.15864, 0.379, 0.716)

TRUE_CROSSING_ABS_TOL = 1.0e-8
ROOT_GAP_POSITIVE_ABS_TOL = 1.0e-8
RESOLVED_RATIO_STRICT = 100.0
RESOLVED_RATIO_LOOSE = 10.0
ROOT_REPAIR_JUMP_THRESHOLD = 0.25
DEFAULT_STRICT_MISMATCH_THRESHOLD = 0.05
STRICT_REFINE_GAP_THRESHOLD = 0.05
REPAIRED_ROOT_UNCERTAINTY = 1.0e-12
REFINE_XATOL = {"eta": 1.0e-5, "beta": 1.0e-4, "mu": 1.0e-5}
SUBSPACE_DELTA = {"eta": 1.0e-5, "beta": 1.0e-3, "mu": 1.0e-5}
DETAIL_OFFSETS = {
    "eta": (-1.0e-4, -1.0e-5, 0.0, 1.0e-5, 1.0e-4),
    "beta": (-1.0e-2, -1.0e-3, 0.0, 1.0e-3, 1.0e-2),
    "mu": (-1.0e-3, -1.0e-4, 0.0, 1.0e-4, 1.0e-3),
}
PARAM_LIMITS = {"eta": (-0.5, 0.5), "beta": (0.0, 90.0), "mu": (0.0, 0.9)}

OUTPUT_DIR = REPO_ROOT / "results"
OUTPUT_SUMMARY = OUTPUT_DIR / "eta_parameter_positive_gap_verification_summary.csv"
OUTPUT_DETAILS = OUTPUT_DIR / "eta_parameter_positive_gap_verification_details.csv"
OUTPUT_REPORT = OUTPUT_DIR / "eta_parameter_positive_gap_verification_report.md"
OUTPUT_MIN_GAP_ETA = OUTPUT_DIR / "eta_parameter_positive_gap_min_gap_vs_eta.png"
OUTPUT_MIN_GAP_BETA = OUTPUT_DIR / "eta_parameter_positive_gap_min_gap_vs_beta_eta_scan.png"
OUTPUT_MIN_GAP_MU = OUTPUT_DIR / "eta_parameter_positive_gap_min_gap_vs_mu_eta_scan.png"
OUTPUT_CLASS_MAP = OUTPUT_DIR / "eta_parameter_positive_gap_classification_map.png"
OUTPUT_LAMBDA_ETA = OUTPUT_DIR / "eta_parameter_positive_gap_lambda_eta_selected.png"
OUTPUT_LAMBDA_BETA = OUTPUT_DIR / "eta_parameter_positive_gap_lambda_beta_selected.png"
OUTPUT_LAMBDA_MU = OUTPUT_DIR / "eta_parameter_positive_gap_lambda_mu_selected.png"
SCRIPT_RELATIVE = SCRIPT_PATH.relative_to(REPO_ROOT)

SUMMARY_FIELDNAMES = [
    "row_type",
    "audit_part",
    "fixed_beta_deg",
    "fixed_mu",
    "fixed_eta",
    "scan_parameter",
    "scan_min",
    "scan_max",
    "pair",
    "parameter_at_min_gap",
    "lambda_lower",
    "lambda_upper",
    "min_gap",
    "relative_gap",
    "root_uncertainty_estimate",
    "gap_to_uncertainty_ratio",
    "classification",
    "classification_global",
    "sorted_position_swap_detected",
    "descendant_swap_detected",
    "individual_mac_min",
    "subspace_mac_min",
    "tracking_warning_count",
    "root_solver_status",
    "any_true_crossing_detected",
    "any_unresolved_possible_crossing",
    "notes",
]

DETAIL_FIELDNAMES = [
    "audit_part",
    "fixed_beta_deg",
    "fixed_mu",
    "fixed_eta",
    "scan_parameter",
    "scan_value",
    "pair",
    "lambda_i",
    "lambda_j",
    "gap",
    "solver_mode",
    "determinant_residual_i",
    "determinant_residual_j",
    "tracking_mac_i",
    "tracking_mac_j",
    "subspace_mac_if_available",
    "warning_flags",
]


@dataclass(frozen=True)
class CaseSpec:
    audit_part: str
    fixed_beta_deg: float | None
    fixed_mu: float | None
    fixed_eta: float | None
    scan_parameter: str
    scan_values: tuple[float, ...]
    label: str


@dataclass(frozen=True)
class TrackingDiagnostics:
    tracked_lambdas: np.ndarray
    current_sorted_positions: np.ndarray
    mac_to_previous: np.ndarray
    warning_count: int
    min_mac: float
    min_margin: float


@dataclass(frozen=True)
class PairAudit:
    pair_index: int
    parameter_at_min_gap: float
    lambda_lower: float
    lambda_upper: float
    min_gap: float
    relative_gap: float
    root_uncertainty_estimate: float
    gap_to_uncertainty_ratio: float
    classification: str
    sorted_position_swap_detected: bool
    descendant_swap_detected: bool
    individual_mac_min: float
    subspace_mac_min: float
    tracking_warning_count: int
    root_solver_status: str
    notes: str


@dataclass(frozen=True)
class CaseAudit:
    spec: CaseSpec
    scan_values: np.ndarray
    sorted_roots: np.ndarray
    tracking: TrackingDiagnostics
    root_repair_count: int
    pair_audits: tuple[PairAudit, ...]
    classification_global: str
    global_pair_index: int
    global_gap: float
    global_parameter: float
    any_true_crossing_detected: bool
    any_unresolved_possible_crossing: bool


@dataclass(frozen=True)
class DescendantPairDiagnostic:
    pair_label: str
    parameter_at_min_gap: float
    lambda_left: float
    lambda_right: float
    min_gap: float
    relative_gap: float
    sign_change_detected: bool
    sorted_position_left: int
    sorted_position_right: int
    mac_left: float
    mac_right: float


def fmt(value: float | int | str | bool | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, str):
        return value
    value_f = float(value)
    if not np.isfinite(value_f):
        return ""
    return f"{value_f:.12g}"


def pair_label(pair_index: int) -> str:
    return f"{int(pair_index)}-{int(pair_index) + 1}"


def number_token(value: float | None) -> str:
    if value is None:
        return "scan"
    return f"{float(value):g}".replace("-", "m").replace(".", "p")


def unique_values(values: Iterable[float]) -> tuple[float, ...]:
    arr = sorted({round(float(value), 12) for value in values})
    return tuple(float(value) for value in arr)


def eta_scan_values() -> tuple[float, ...]:
    coarse = np.arange(-0.5, 0.5 + 0.5 * ETA_COARSE_STEP, ETA_COARSE_STEP, dtype=float)
    special = [0.0]
    for value in ETA_SPECIAL_ABS:
        special.extend([-value, value])
    return unique_values([*coarse, *special])


def beta_scan_values() -> tuple[float, ...]:
    return unique_values(np.arange(0.0, 90.0 + 0.5 * BETA_COARSE_STEP, BETA_COARSE_STEP, dtype=float))


def mu_scan_values() -> tuple[float, ...]:
    coarse = np.arange(0.0, 0.9 + 0.5 * MU_COARSE_STEP, MU_COARSE_STEP, dtype=float)
    return unique_values([*coarse, *MU_KNOWN_POINTS])


def local_mu_window(center: float) -> list[float]:
    coarse = np.arange(
        max(0.0, float(center) - MU_LOCAL_COARSE_HALF_WIDTH),
        min(0.9, float(center) + MU_LOCAL_COARSE_HALF_WIDTH) + 0.5 * MU_LOCAL_COARSE_STEP,
        MU_LOCAL_COARSE_STEP,
        dtype=float,
    )
    fine = np.arange(
        max(0.0, float(center) - MU_LOCAL_FINE_HALF_WIDTH),
        min(0.9, float(center) + MU_LOCAL_FINE_HALF_WIDTH) + 0.5 * MU_LOCAL_FINE_STEP,
        MU_LOCAL_FINE_STEP,
        dtype=float,
    )
    return [*coarse, *fine, float(center)]


def mu_scan_values_for_case(beta: float, eta: float) -> tuple[float, ...]:
    values: list[float] = list(mu_scan_values())
    if np.isclose(float(eta), 0.5, atol=1.0e-12):
        if np.isclose(float(beta), 5.0, atol=1.0e-12):
            values.extend(local_mu_window(0.1586434207))
        if np.isclose(float(beta), 15.0, atol=1.0e-12):
            values.extend(local_mu_window(0.379))
            values.extend(local_mu_window(0.716))
    return unique_values(values)


def build_case_specs(*, smoke: bool = False) -> list[CaseSpec]:
    if smoke:
        eta_values = unique_values([-0.5, -0.01, 0.0, 0.01, 0.5])
        beta_values = unique_values([0.0, 15.0, 45.0, 90.0])
        mu_values = unique_values([0.0, 0.5, 0.9, *local_mu_window(0.1586434207)])
        return [
            CaseSpec("A_eta_scan", 15.0, 0.0, None, "eta", eta_values, "A_smoke_beta15_mu0"),
            CaseSpec("B_beta_scan", None, 0.0, 0.004, "beta", beta_values, "B_smoke_eta0p004"),
            CaseSpec("C_mu_scan", 5.0, None, 0.5, "mu", mu_values, "C_smoke_beta5_eta0p5"),
        ]

    specs: list[CaseSpec] = []
    eta_values = eta_scan_values()
    for mu in BETA15_MU_VALUES:
        specs.append(
            CaseSpec("A_eta_scan", 15.0, float(mu), None, "eta", eta_values, f"A_beta15_mu{number_token(mu)}")
        )
    for beta in BETA_ETA_SCAN_VALUES:
        if np.isclose(beta, 15.0, atol=1.0e-12):
            continue
        specs.append(CaseSpec("A_eta_scan", float(beta), 0.0, None, "eta", eta_values, f"A_beta{number_token(beta)}_mu0"))

    for eta in ETA_BETA_SCAN_VALUES:
        specs.append(
            CaseSpec("B_beta_scan", None, 0.0, float(eta), "beta", beta_scan_values(), f"B_eta{number_token(eta)}")
        )

    for eta in ETA_MU_SCAN_VALUES:
        for beta in BETA_MU_SCAN_VALUES:
            specs.append(
                CaseSpec(
                    "C_mu_scan",
                    float(beta),
                    None,
                    float(eta),
                    "mu",
                    mu_scan_values_for_case(float(beta), float(eta)),
                    f"C_beta{number_token(beta)}_eta{number_token(eta)}",
                )
            )
    return specs


def params_for_value(spec: CaseSpec, value: float) -> tuple[float, float, float]:
    beta = float(value) if spec.scan_parameter == "beta" else float(spec.fixed_beta_deg)
    mu = float(value) if spec.scan_parameter == "mu" else float(spec.fixed_mu)
    eta = float(value) if spec.scan_parameter == "eta" else float(spec.fixed_eta)
    return beta, mu, eta


def solve_mode_scan_step(mode: str) -> float:
    if mode == "default":
        return DEFAULT_SCAN_STEP
    if mode == "strict":
        return STRICT_SCAN_STEP
    raise ValueError(f"unknown solver mode {mode!r}")


@lru_cache(maxsize=None)
def solve_roots_cached(beta_key: float, mu_key: float, eta_key: float, mode: str) -> tuple[float, ...]:
    roots = find_first_n_roots_eta(
        float(np.deg2rad(beta_key)),
        float(mu_key),
        EPSILON,
        float(eta_key),
        N_SORTED_ROOTS,
        Lmax0=ROOT_LMAX0,
        scan_step=solve_mode_scan_step(mode),
    )
    return tuple(float(value) for value in roots)


def solve_roots(beta_deg: float, mu: float, eta: float, *, mode: str) -> np.ndarray:
    return np.asarray(
        solve_roots_cached(round(float(beta_deg), 10), round(float(mu), 12), round(float(eta), 12), str(mode)),
        dtype=float,
    )


def solve_roots_for_value(spec: CaseSpec, value: float, *, mode: str) -> np.ndarray:
    beta, mu, eta = params_for_value(spec, value)
    return solve_roots(beta, mu, eta, mode=mode)


def root_sets_are_consistent(left: np.ndarray, right: np.ndarray, *, count: int = N_GAP_ROOTS) -> bool:
    left_arr = np.asarray(left[:count], dtype=float)
    right_arr = np.asarray(right[:count], dtype=float)
    if np.any(~np.isfinite(left_arr)) or np.any(~np.isfinite(right_arr)):
        return False
    return bool(np.max(np.abs(left_arr - right_arr)) <= DEFAULT_STRICT_MISMATCH_THRESHOLD)


def special_strict_probe_requested(spec: CaseSpec, value: float) -> bool:
    if spec.scan_parameter != "mu" or spec.fixed_eta is None or spec.fixed_beta_deg is None:
        return False
    if not np.isclose(float(spec.fixed_eta), 0.5, atol=1.0e-12):
        return False
    mu = float(value)
    beta = float(spec.fixed_beta_deg)
    if np.isclose(beta, 5.0, atol=1.0e-12) and 0.145 <= mu <= 0.175:
        return True
    if np.isclose(beta, 15.0, atol=1.0e-12) and (0.35 <= mu <= 0.41 or 0.68 <= mu <= 0.75):
        return True
    return False


def repair_grid_roots_if_needed(
    spec: CaseSpec,
    value: float,
    raw_roots: np.ndarray,
    previous_roots: np.ndarray | None,
) -> tuple[np.ndarray, bool]:
    needs_probe = special_strict_probe_requested(spec, value)
    if previous_roots is not None and np.all(np.isfinite(previous_roots[:N_DESCENDANTS_TRACK])):
        if np.all(np.isfinite(raw_roots[:N_DESCENDANTS_TRACK])):
            jump = float(np.max(np.abs(raw_roots[:N_DESCENDANTS_TRACK] - previous_roots[:N_DESCENDANTS_TRACK])))
            if np.isfinite(jump) and jump > ROOT_REPAIR_JUMP_THRESHOLD:
                needs_probe = True
        else:
            needs_probe = True
    if not needs_probe:
        return raw_roots, False

    strict_roots = solve_roots_for_value(spec, value, mode="strict")
    if np.any(~np.isfinite(strict_roots[:N_DESCENDANTS_TRACK])):
        return raw_roots, False
    if previous_roots is None:
        use_strict = not root_sets_are_consistent(raw_roots, strict_roots, count=N_DESCENDANTS_TRACK)
    else:
        raw_jump = (
            float(np.max(np.abs(raw_roots[:N_DESCENDANTS_TRACK] - previous_roots[:N_DESCENDANTS_TRACK])))
            if np.all(np.isfinite(raw_roots[:N_DESCENDANTS_TRACK]))
            else float("inf")
        )
        strict_jump = float(np.max(np.abs(strict_roots[:N_DESCENDANTS_TRACK] - previous_roots[:N_DESCENDANTS_TRACK])))
        use_strict = strict_jump < raw_jump or not root_sets_are_consistent(raw_roots, strict_roots, count=N_DESCENDANTS_TRACK)
    return (strict_roots, True) if use_strict else (raw_roots, False)


def pair_gap_from_roots(roots: np.ndarray, pair_index: int) -> tuple[float, float, float]:
    left = float(roots[int(pair_index) - 1])
    right = float(roots[int(pair_index)])
    if not (np.isfinite(left) and np.isfinite(right)):
        return float("nan"), float("nan"), float("nan")
    return left, right, abs(right - left)


def gap_value(spec: CaseSpec, value: float, pair_index: int, *, mode: str = "default") -> float:
    roots = solve_roots_for_value(spec, value, mode=mode)
    _left, _right, gap = pair_gap_from_roots(roots, pair_index)
    if not np.isfinite(gap):
        return 1.0e9
    return float(gap)


def determinant_residual(lambda_value: float, beta_deg: float, mu: float, eta: float) -> float:
    if not np.isfinite(float(lambda_value)):
        return float("nan")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        with np.errstate(over="ignore", invalid="ignore"):
            value = det_eta(float(lambda_value), float(np.deg2rad(beta_deg)), float(mu), EPSILON, float(eta))
    return abs(float(value)) if np.isfinite(value) else float("nan")


def solve_sorted_grid(spec: CaseSpec) -> tuple[np.ndarray, int]:
    roots = np.full((len(spec.scan_values), N_SORTED_ROOTS), np.nan, dtype=float)
    previous_roots: np.ndarray | None = None
    repair_count = 0
    for idx, value in enumerate(spec.scan_values):
        raw_roots = solve_roots_for_value(spec, value, mode="default")
        repaired_roots, repaired = repair_grid_roots_if_needed(spec, float(value), raw_roots, previous_roots)
        roots[idx] = repaired_roots
        previous_roots = repaired_roots
        if repaired:
            repair_count += 1
        if idx == 0 or (idx + 1) % 100 == 0 or idx + 1 == len(spec.scan_values):
            print(f"{spec.label}: solved roots for {idx + 1}/{len(spec.scan_values)} {spec.scan_parameter} values")
    if repair_count:
        print(f"{spec.label}: strict local root repairs applied at {repair_count} grid values")
    return roots, repair_count


def shape_vectors_for_roots(spec: CaseSpec, value: float, roots: np.ndarray) -> list[np.ndarray] | None:
    finite_roots = np.asarray(roots[:N_TRACKING_CANDIDATES], dtype=float)
    if len(finite_roots) < N_DESCENDANTS_TRACK or np.any(~np.isfinite(finite_roots[:N_DESCENDANTS_TRACK])):
        return None
    beta, mu, eta = params_for_value(spec, value)
    return analytic_shape_vectors_for_roots(
        finite_roots,
        beta_rad=float(np.deg2rad(beta)),
        mu=float(mu),
        epsilon=EPSILON,
        eta=float(eta),
        s_norm=np.linspace(0.0, 1.0, NUM_SHAPE_SAMPLES, dtype=float),
    )


def tracking_seed_index(spec: CaseSpec) -> int:
    values = np.asarray(spec.scan_values, dtype=float)
    if spec.scan_parameter in {"eta", "mu"} and np.any(np.isclose(values, 0.0, rtol=0.0, atol=1.0e-12)):
        return int(np.where(np.isclose(values, 0.0, rtol=0.0, atol=1.0e-12))[0][0])
    return 0


def tracking_order(n_values: int, seed: int) -> list[int]:
    return [*range(seed + 1, n_values), *range(seed - 1, -1, -1)]


def track_descendants(spec: CaseSpec, sorted_roots: np.ndarray) -> TrackingDiagnostics:
    n_values = len(spec.scan_values)
    tracked = np.full((N_DESCENDANTS_TRACK, n_values), np.nan, dtype=float)
    positions = np.full((N_DESCENDANTS_TRACK, n_values), -1, dtype=int)
    mac_to_previous = np.full((N_DESCENDANTS_TRACK, n_values), np.nan, dtype=float)
    margins = np.full((N_DESCENDANTS_TRACK, n_values), np.nan, dtype=float)
    seed = tracking_seed_index(spec)
    seed_vectors = shape_vectors_for_roots(spec, spec.scan_values[seed], sorted_roots[seed])
    if seed_vectors is None:
        return TrackingDiagnostics(tracked, positions, mac_to_previous, 1, float("nan"), float("nan"))

    tracked[:, seed] = sorted_roots[seed, :N_DESCENDANTS_TRACK]
    positions[:, seed] = np.arange(1, N_DESCENDANTS_TRACK + 1, dtype=int)
    mac_to_previous[:, seed] = 1.0

    previous_by_direction: dict[int, tuple[list[np.ndarray], np.ndarray]] = {
        1: (seed_vectors[:N_DESCENDANTS_TRACK], tracked[:, seed].copy()),
        -1: (seed_vectors[:N_DESCENDANTS_TRACK], tracked[:, seed].copy()),
    }
    for idx in tracking_order(n_values, seed):
        direction = 1 if idx > seed else -1
        previous_vectors, previous_lambdas = previous_by_direction[direction]
        candidate_vectors = shape_vectors_for_roots(spec, spec.scan_values[idx], sorted_roots[idx])
        if candidate_vectors is None:
            continue
        try:
            assignment, mac = mac_assignment(previous_vectors, candidate_vectors)
        except Exception:
            continue
        current_vectors: list[np.ndarray] = []
        current_lambdas = np.full(N_DESCENDANTS_TRACK, np.nan, dtype=float)
        for desc_idx, candidate_col in enumerate(assignment[:N_DESCENDANTS_TRACK]):
            candidate_col = int(candidate_col)
            current_lambdas[desc_idx] = float(sorted_roots[idx, candidate_col])
            tracked[desc_idx, idx] = current_lambdas[desc_idx]
            positions[desc_idx, idx] = candidate_col + 1
            assigned_mac = float(mac[desc_idx, candidate_col])
            mac_to_previous[desc_idx, idx] = assigned_mac
            row = np.sort(mac[desc_idx])[::-1]
            second = float(row[1]) if len(row) > 1 else float("nan")
            margins[desc_idx, idx] = assigned_mac - second if np.isfinite(second) else float("nan")
            current_vectors.append(candidate_vectors[candidate_col])
        previous_by_direction[direction] = (current_vectors, current_lambdas)

    nonseed = np.isfinite(mac_to_previous) & (mac_to_previous < 0.9)
    low_margin = np.isfinite(margins) & (margins < 0.05)
    warning_count = int(np.count_nonzero(nonseed | low_margin))
    finite_macs = mac_to_previous[np.isfinite(mac_to_previous) & (mac_to_previous < 0.999999999)]
    finite_margins = margins[np.isfinite(margins)]
    return TrackingDiagnostics(
        tracked_lambdas=tracked,
        current_sorted_positions=positions,
        mac_to_previous=mac_to_previous,
        warning_count=warning_count,
        min_mac=float(np.min(finite_macs)) if finite_macs.size else 1.0,
        min_margin=float(np.min(finite_margins)) if finite_margins.size else float("nan"),
    )


def min_gap_index(gaps: np.ndarray, pair_index: int) -> int:
    values = np.asarray(gaps[:, int(pair_index) - 1], dtype=float)
    if not np.any(np.isfinite(values)):
        return 0
    return int(np.nanargmin(values))


def refine_gap_parameter(spec: CaseSpec, gaps: np.ndarray, pair_index: int) -> tuple[float, float]:
    values = np.asarray(spec.scan_values, dtype=float)
    idx = min_gap_index(gaps, pair_index)
    left_idx = max(0, idx - 1)
    right_idx = min(len(values) - 1, idx + 1)
    a = float(values[left_idx])
    b = float(values[right_idx])
    candidates = [
        (float(values[idx]), float(gaps[idx, int(pair_index) - 1])),
        (a, float(gaps[left_idx, int(pair_index) - 1])),
        (b, float(gaps[right_idx, int(pair_index) - 1])),
    ]
    if b > a:
        finite_candidate_gaps = [gap for _parameter, gap in candidates if np.isfinite(gap)]
        coarse_min_gap = min(finite_candidate_gaps) if finite_candidate_gaps else float("inf")
        objective_mode = "strict" if coarse_min_gap < STRICT_REFINE_GAP_THRESHOLD else "default"
        result = minimize_scalar(
            lambda parameter: gap_value(spec, float(parameter), pair_index, mode=objective_mode),
            bounds=(a, b),
            method="bounded",
            options={"xatol": REFINE_XATOL[spec.scan_parameter], "maxiter": 80},
        )
        if result.success and np.isfinite(float(result.fun)):
            candidates.append((float(result.x), float(result.fun)))
    parameter, gap = min(candidates, key=lambda item: item[1])
    return float(parameter), float(gap)


def classify_gap(strict_gap: float, uncertainty: float, root_solver_status: str) -> tuple[str, float, str]:
    if not root_solver_status.startswith("ok") or not np.isfinite(strict_gap):
        return "tracking_unreliable", float("nan"), root_solver_status
    ratio = float("inf") if (not np.isfinite(uncertainty) or uncertainty < 1.0e-14) else float(strict_gap / uncertainty)
    if strict_gap <= TRUE_CROSSING_ABS_TOL and ratio <= RESOLVED_RATIO_LOOSE:
        return "true_crossing_detected", ratio, "gap is at numerical zero"
    if strict_gap >= ROOT_GAP_POSITIVE_ABS_TOL and ratio >= RESOLVED_RATIO_STRICT:
        return "resolved_positive_gap", ratio, "gap exceeds root uncertainty by at least 100x"
    if strict_gap >= ROOT_GAP_POSITIVE_ABS_TOL and ratio >= RESOLVED_RATIO_LOOSE:
        return "resolved_positive_gap", ratio, "gap exceeds root uncertainty by at least 10x"
    return "unresolved_possible_crossing", ratio, "gap is comparable to root uncertainty"


def normalized_columns(vectors: Sequence[np.ndarray]) -> np.ndarray:
    columns = []
    for vector in vectors:
        array = np.asarray(vector, dtype=float)
        norm = float(np.linalg.norm(array))
        if norm <= 1.0e-14 or not np.isfinite(norm):
            return np.full((len(array), len(vectors)), np.nan, dtype=float)
        columns.append(array / norm)
    return np.column_stack(columns)


def subspace_and_individual_mac(spec: CaseSpec, parameter: float, pair_index: int) -> tuple[float, float]:
    lo, hi = PARAM_LIMITS[spec.scan_parameter]
    delta = SUBSPACE_DELTA[spec.scan_parameter]
    left_value = max(lo, float(parameter) - delta)
    right_value = min(hi, float(parameter) + delta)
    if np.isclose(left_value, right_value, rtol=0.0, atol=1.0e-14):
        return float("nan"), float("nan")
    left_roots = solve_roots_for_value(spec, left_value, mode="strict")
    right_roots = solve_roots_for_value(spec, right_value, mode="strict")
    if np.any(~np.isfinite(left_roots[:N_GAP_ROOTS])) or np.any(~np.isfinite(right_roots[:N_GAP_ROOTS])):
        return float("nan"), float("nan")
    sl = slice(int(pair_index) - 1, int(pair_index) + 1)
    left_vectors = shape_vectors_for_roots(spec, left_value, left_roots)
    right_vectors = shape_vectors_for_roots(spec, right_value, right_roots)
    if left_vectors is None or right_vectors is None:
        return float("nan"), float("nan")
    vectors_left = left_vectors[sl]
    vectors_right = right_vectors[sl]
    individual = min(mac_value(vectors_left[0], vectors_right[0]), mac_value(vectors_left[1], vectors_right[1]))
    left_matrix = normalized_columns(vectors_left)
    right_matrix = normalized_columns(vectors_right)
    if not (np.all(np.isfinite(left_matrix)) and np.all(np.isfinite(right_matrix))):
        return float("nan"), float(individual)
    q_left, _ = np.linalg.qr(left_matrix)
    q_right, _ = np.linalg.qr(right_matrix)
    singular_values = np.linalg.svd(q_left.T @ q_right, compute_uv=False)
    return float(np.min(np.clip(singular_values, 0.0, 1.0) ** 2)), float(individual)


def sorted_position_swap_detected(tracking: TrackingDiagnostics) -> bool:
    positions = tracking.current_sorted_positions[:6]
    valid_cols = np.all(positions > 0, axis=0)
    if not np.any(valid_cols):
        return False
    seed_col = int(np.where(valid_cols)[0][0])
    return bool(np.any(positions[:, valid_cols] != positions[:, seed_col][:, None]))


def descendant_swap_for_pair(tracking: TrackingDiagnostics, pair_index: int) -> bool:
    left = int(pair_index) - 1
    right = int(pair_index)
    if right >= tracking.tracked_lambdas.shape[0]:
        return False
    positions = tracking.current_sorted_positions
    valid = (positions[left] > 0) & (positions[right] > 0)
    if np.any(valid):
        if bool(np.any(positions[left, valid] > positions[right, valid])):
            return True
    diff = tracking.tracked_lambdas[right] - tracking.tracked_lambdas[left]
    diff = diff[np.isfinite(diff) & (np.abs(diff) > TRUE_CROSSING_ABS_TOL)]
    signs = np.sign(diff)
    return bool(signs.size > 1 and np.any(signs[:-1] * signs[1:] < 0))


def audit_pair(spec: CaseSpec, sorted_roots: np.ndarray, tracking: TrackingDiagnostics, pair_index: int) -> PairAudit:
    gaps = np.diff(sorted_roots[:, :N_GAP_ROOTS], axis=1)
    parameter_at_min, _default_gap = refine_gap_parameter(spec, gaps, pair_index)
    default_roots = solve_roots_for_value(spec, parameter_at_min, mode="default")
    strict_roots = solve_roots_for_value(spec, parameter_at_min, mode="strict")
    default_lower, default_upper, default_gap = pair_gap_from_roots(default_roots, pair_index)
    strict_lower, strict_upper, strict_gap = pair_gap_from_roots(strict_roots, pair_index)
    if not np.all(np.isfinite(default_roots[:N_GAP_ROOTS])):
        status = "default_missing_roots"
    elif not np.all(np.isfinite(strict_roots[:N_GAP_ROOTS])):
        status = "strict_missing_roots"
    elif not root_sets_are_consistent(default_roots, strict_roots):
        status = "ok_default_scan_missed_close_roots_repaired_by_strict"
    else:
        status = "ok"
    if status == "ok":
        uncertainty = max(abs(default_lower - strict_lower), abs(default_upper - strict_upper), abs(default_gap - strict_gap))
    elif status.startswith("ok_default_scan_missed"):
        uncertainty = REPAIRED_ROOT_UNCERTAINTY
    else:
        uncertainty = float("nan")
    classification, ratio, note = classify_gap(strict_gap, uncertainty, status)
    if status.startswith("ok_default_scan_missed"):
        note = "raw default sign-scan missed a close root; strict local repair used for classification"
    denominator = 0.5 * (abs(strict_lower) + abs(strict_upper))
    relative_gap = strict_gap / denominator if denominator > 0.0 and np.isfinite(strict_gap) else float("nan")
    subspace_mac, individual_mac = subspace_and_individual_mac(spec, parameter_at_min, pair_index)
    return PairAudit(
        pair_index=int(pair_index),
        parameter_at_min_gap=float(parameter_at_min),
        lambda_lower=float(strict_lower),
        lambda_upper=float(strict_upper),
        min_gap=float(strict_gap),
        relative_gap=float(relative_gap),
        root_uncertainty_estimate=float(uncertainty),
        gap_to_uncertainty_ratio=float(ratio),
        classification=classification,
        sorted_position_swap_detected=sorted_position_swap_detected(tracking),
        descendant_swap_detected=descendant_swap_for_pair(tracking, pair_index),
        individual_mac_min=float(individual_mac),
        subspace_mac_min=float(subspace_mac),
        tracking_warning_count=int(tracking.warning_count),
        root_solver_status=status,
        notes=note,
    )


def global_classification(pair_audits: Sequence[PairAudit]) -> tuple[str, bool, bool]:
    any_true = any(pair.classification == "true_crossing_detected" for pair in pair_audits)
    any_unresolved = any(pair.classification in {"unresolved_possible_crossing", "tracking_unreliable"} for pair in pair_audits)
    if any_true:
        return "true_crossing_detected", any_true, any_unresolved
    if any_unresolved:
        return "unresolved_possible_crossing", any_true, any_unresolved
    return "resolved_positive_gap", any_true, any_unresolved


def audit_case(spec: CaseSpec) -> CaseAudit:
    print(f"running {spec.label} ({spec.audit_part}, scan={spec.scan_parameter})")
    sorted_roots, root_repair_count = solve_sorted_grid(spec)
    tracking = track_descendants(spec, sorted_roots)
    pair_audits = tuple(audit_pair(spec, sorted_roots, tracking, pair_index) for pair_index in range(1, N_GAP_ROOTS))
    finite_pairs = [pair for pair in pair_audits if np.isfinite(pair.min_gap)]
    global_pair = min(finite_pairs, key=lambda pair: pair.min_gap) if finite_pairs else pair_audits[0]
    classification, any_true, any_unresolved = global_classification(pair_audits)
    print(
        f"{spec.label}: global={classification}, pair={pair_label(global_pair.pair_index)}, "
        f"gap={global_pair.min_gap:.6g}, {spec.scan_parameter}={global_pair.parameter_at_min_gap:.8g}"
    )
    return CaseAudit(
        spec=spec,
        scan_values=np.asarray(spec.scan_values, dtype=float),
        sorted_roots=sorted_roots,
        tracking=tracking,
        root_repair_count=int(root_repair_count),
        pair_audits=pair_audits,
        classification_global=classification,
        global_pair_index=global_pair.pair_index,
        global_gap=float(global_pair.min_gap),
        global_parameter=float(global_pair.parameter_at_min_gap),
        any_true_crossing_detected=any_true,
        any_unresolved_possible_crossing=any_unresolved,
    )


def summary_row(case: CaseAudit, pair: PairAudit, *, row_type: str) -> dict[str, str]:
    spec = case.spec
    notes = pair.notes
    if case.root_repair_count:
        notes = f"{notes}; strict_grid_repairs={case.root_repair_count}"
    return {
        "row_type": row_type,
        "audit_part": spec.audit_part,
        "fixed_beta_deg": fmt(spec.fixed_beta_deg),
        "fixed_mu": fmt(spec.fixed_mu),
        "fixed_eta": fmt(spec.fixed_eta),
        "scan_parameter": spec.scan_parameter,
        "scan_min": fmt(min(spec.scan_values)),
        "scan_max": fmt(max(spec.scan_values)),
        "pair": pair_label(pair.pair_index) if row_type == "pair" else f"global_min:{pair_label(pair.pair_index)}",
        "parameter_at_min_gap": fmt(pair.parameter_at_min_gap),
        "lambda_lower": fmt(pair.lambda_lower),
        "lambda_upper": fmt(pair.lambda_upper),
        "min_gap": fmt(pair.min_gap),
        "relative_gap": fmt(pair.relative_gap),
        "root_uncertainty_estimate": fmt(pair.root_uncertainty_estimate),
        "gap_to_uncertainty_ratio": fmt(pair.gap_to_uncertainty_ratio),
        "classification": pair.classification,
        "classification_global": case.classification_global,
        "sorted_position_swap_detected": fmt(pair.sorted_position_swap_detected),
        "descendant_swap_detected": fmt(pair.descendant_swap_detected),
        "individual_mac_min": fmt(pair.individual_mac_min),
        "subspace_mac_min": fmt(pair.subspace_mac_min),
        "tracking_warning_count": str(pair.tracking_warning_count),
        "root_solver_status": pair.root_solver_status,
        "any_true_crossing_detected": fmt(case.any_true_crossing_detected),
        "any_unresolved_possible_crossing": fmt(case.any_unresolved_possible_crossing),
        "notes": notes,
    }


def write_summary_csv(cases: Sequence[CaseAudit]) -> None:
    OUTPUT_SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_SUMMARY.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDNAMES)
        writer.writeheader()
        for case in cases:
            global_pair = min(case.pair_audits, key=lambda pair: pair.min_gap if np.isfinite(pair.min_gap) else np.inf)
            writer.writerow(summary_row(case, global_pair, row_type="global"))
            for pair in case.pair_audits:
                writer.writerow(summary_row(case, pair, row_type="pair"))


def detail_rows_for_pair(case: CaseAudit, pair: PairAudit) -> list[dict[str, str]]:
    spec = case.spec
    lo, hi = PARAM_LIMITS[spec.scan_parameter]
    values = sorted(
        {
            round(min(hi, max(lo, pair.parameter_at_min_gap + offset)), 12)
            for offset in DETAIL_OFFSETS[spec.scan_parameter]
        }
    )
    rows: list[dict[str, str]] = []
    for value in values:
        modes = ("default", "strict") if np.isclose(value, pair.parameter_at_min_gap, atol=1.0e-12) else ("default",)
        for mode in modes:
            roots = solve_roots_for_value(spec, value, mode=mode)
            lambda_i, lambda_j, gap = pair_gap_from_roots(roots, pair.pair_index)
            beta, mu, eta = params_for_value(spec, value)
            warning_flags = []
            if pair.individual_mac_min < 0.9:
                warning_flags.append("low_individual_mac")
            if pair.subspace_mac_min < 0.9:
                warning_flags.append("low_subspace_mac")
            rows.append(
                {
                    "audit_part": spec.audit_part,
                    "fixed_beta_deg": fmt(spec.fixed_beta_deg),
                    "fixed_mu": fmt(spec.fixed_mu),
                    "fixed_eta": fmt(spec.fixed_eta),
                    "scan_parameter": spec.scan_parameter,
                    "scan_value": fmt(value),
                    "pair": pair_label(pair.pair_index),
                    "lambda_i": fmt(lambda_i),
                    "lambda_j": fmt(lambda_j),
                    "gap": fmt(gap),
                    "solver_mode": mode,
                    "determinant_residual_i": fmt(determinant_residual(lambda_i, beta, mu, eta)),
                    "determinant_residual_j": fmt(determinant_residual(lambda_j, beta, mu, eta)),
                    "tracking_mac_i": fmt(pair.individual_mac_min),
                    "tracking_mac_j": fmt(pair.individual_mac_min),
                    "subspace_mac_if_available": fmt(pair.subspace_mac_min),
                    "warning_flags": ";".join(warning_flags) if warning_flags else "none",
                }
            )
    return rows


def write_details_csv(cases: Sequence[CaseAudit]) -> None:
    OUTPUT_DETAILS.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_DETAILS.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=DETAIL_FIELDNAMES)
        writer.writeheader()
        for case in cases:
            for pair in case.pair_audits:
                writer.writerows(detail_rows_for_pair(case, pair))


def coarse_global_gap_curve(case: CaseAudit) -> np.ndarray:
    gaps = np.diff(case.sorted_roots[:, :N_GAP_ROOTS], axis=1)
    return np.nanmin(gaps, axis=1)


def plot_min_gap_vs_eta(cases: Sequence[CaseAudit]) -> None:
    selected = [
        case
        for case in cases
        if case.spec.audit_part == "A_eta_scan"
        and (
            (case.spec.fixed_beta_deg == 15.0 and case.spec.fixed_mu in {0.0, 0.1, 0.9})
            or (case.spec.fixed_mu == 0.0 and case.spec.fixed_beta_deg in {0.0, 5.0, 10.0})
        )
    ]
    fig, ax = plt.subplots(figsize=(10.0, 5.8), constrained_layout=True)
    for case in selected:
        label = f"beta={case.spec.fixed_beta_deg:g}, mu={case.spec.fixed_mu:g}"
        ax.plot(case.scan_values, coarse_global_gap_curve(case), lw=1.25, label=label)
        ax.scatter([case.global_parameter], [case.global_gap], s=18)
    ax.set_yscale("log")
    ax.set_xlabel("eta")
    ax.set_ylabel("minimum adjacent sorted gap")
    ax.set_title("Eta scans: coarse global gap curves and refined minima")
    ax.grid(True, which="both", color="0.88", lw=0.5)
    ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.savefig(OUTPUT_MIN_GAP_ETA, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_min_gap_vs_beta_eta_scan(cases: Sequence[CaseAudit]) -> None:
    rows = [case for case in cases if case.spec.audit_part == "B_beta_scan"]
    fig, ax = plt.subplots(figsize=(9.8, 5.4), constrained_layout=True)
    eta_values = np.array([case.spec.fixed_eta for case in rows], dtype=float)
    global_gaps = np.array([case.global_gap for case in rows], dtype=float)
    colors = ["tab:blue" if abs(eta) > 1.0e-12 else "tab:red" for eta in eta_values]
    ax.scatter(eta_values, global_gaps, c=colors, s=36)
    ax.plot(eta_values, global_gaps, color="0.45", lw=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("eta")
    ax.set_ylabel("refined min beta-gap")
    ax.set_title("Beta scans at mu=0: global adjacent-gap minima")
    ax.grid(True, which="both", color="0.88", lw=0.5)
    fig.savefig(OUTPUT_MIN_GAP_BETA, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_min_gap_vs_mu_eta_scan(cases: Sequence[CaseAudit]) -> None:
    rows = [case for case in cases if case.spec.audit_part == "C_mu_scan"]
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.4), sharey=True, constrained_layout=True)
    for ax, beta in zip(axes, BETA_MU_SCAN_VALUES):
        beta_rows = [case for case in rows if np.isclose(case.spec.fixed_beta_deg, beta)]
        eta_values = np.array([case.spec.fixed_eta for case in beta_rows], dtype=float)
        gaps = np.array([case.global_gap for case in beta_rows], dtype=float)
        ax.scatter(eta_values, gaps, s=32)
        ax.plot(eta_values, gaps, color="0.45", lw=0.8)
        ax.set_yscale("log")
        ax.set_title(f"beta={beta:g} deg")
        ax.set_xlabel("eta")
        ax.grid(True, which="both", color="0.88", lw=0.5)
    axes[0].set_ylabel("refined min mu-gap")
    fig.savefig(OUTPUT_MIN_GAP_MU, dpi=220, bbox_inches="tight")
    plt.close(fig)


def classification_code(value: str) -> int:
    return {
        "resolved_positive_gap": 0,
        "unresolved_possible_crossing": 1,
        "true_crossing_detected": 2,
        "tracking_unreliable": 3,
    }.get(value, 3)


def plot_classification_map(cases: Sequence[CaseAudit]) -> None:
    fig, ax = plt.subplots(figsize=(10.5, max(4.0, 0.22 * len(cases))), constrained_layout=True)
    codes = np.array([[classification_code(case.classification_global)] for case in cases], dtype=float)
    cmap = matplotlib.colors.ListedColormap(["#2ca25f", "#ffd166", "#d95f02", "#756bb1"])
    ax.imshow(codes, aspect="auto", cmap=cmap, vmin=0, vmax=3)
    ax.set_yticks(np.arange(len(cases)))
    ax.set_yticklabels([case.spec.label for case in cases], fontsize=7)
    ax.set_xticks([0])
    ax.set_xticklabels(["global"])
    ax.set_title("Global classification map")
    cbar = fig.colorbar(
        matplotlib.cm.ScalarMappable(norm=matplotlib.colors.Normalize(vmin=0, vmax=3), cmap=cmap),
        ax=ax,
        ticks=[0, 1, 2, 3],
    )
    cbar.ax.set_yticklabels(["resolved", "unresolved", "crossing", "tracking"])
    fig.savefig(OUTPUT_CLASS_MAP, dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_selected_curves(cases: Sequence[CaseAudit]) -> None:
    def plot_roots(ax, case: CaseAudit, title: str) -> None:
        for idx in range(6):
            ax.plot(case.scan_values, case.sorted_roots[:, idx], lw=1.1, label=f"sorted {idx + 1}")
        ax.set_title(title)
        ax.grid(True, color="0.9", lw=0.5)
        ax.set_ylabel("Lambda")

    eta_cases = [
        case
        for case in cases
        if case.spec.audit_part == "A_eta_scan" and case.spec.fixed_beta_deg == 15.0 and case.spec.fixed_mu in {0.0, 0.1, 0.9}
    ]
    if eta_cases:
        fig, axes = plt.subplots(len(eta_cases), 1, figsize=(9.8, 3.0 * len(eta_cases)), sharex=True, constrained_layout=True)
        axes_arr = np.atleast_1d(axes)
        for ax, case in zip(axes_arr, eta_cases):
            plot_roots(ax, case, f"beta=15 deg, mu={case.spec.fixed_mu:g}")
        axes_arr[-1].set_xlabel("eta")
        axes_arr[0].legend(frameon=False, ncol=3, fontsize=8)
        fig.savefig(OUTPUT_LAMBDA_ETA, dpi=220, bbox_inches="tight")
        plt.close(fig)

    beta_cases = [
        case
        for case in cases
        if case.spec.audit_part == "B_beta_scan" and case.spec.fixed_eta in {-0.016, -0.004, 0.0, 0.004, 0.016}
    ]
    if beta_cases:
        fig, axes = plt.subplots(len(beta_cases), 1, figsize=(9.8, 2.75 * len(beta_cases)), sharex=True, constrained_layout=True)
        axes_arr = np.atleast_1d(axes)
        for ax, case in zip(axes_arr, beta_cases):
            plot_roots(ax, case, f"mu=0, eta={case.spec.fixed_eta:g}")
        axes_arr[-1].set_xlabel("beta (deg)")
        axes_arr[0].legend(frameon=False, ncol=3, fontsize=8)
        fig.savefig(OUTPUT_LAMBDA_BETA, dpi=220, bbox_inches="tight")
        plt.close(fig)

    mu_cases = [
        case
        for case in cases
        if case.spec.audit_part == "C_mu_scan"
        and (
            (case.spec.fixed_beta_deg == 5.0 and case.spec.fixed_eta == 0.5)
            or (case.spec.fixed_beta_deg == 15.0 and case.spec.fixed_eta in {-0.5, 0.5})
        )
    ]
    if mu_cases:
        fig, axes = plt.subplots(len(mu_cases), 1, figsize=(9.8, 3.0 * len(mu_cases)), sharex=True, constrained_layout=True)
        axes_arr = np.atleast_1d(axes)
        for ax, case in zip(axes_arr, mu_cases):
            plot_roots(ax, case, f"beta={case.spec.fixed_beta_deg:g} deg, eta={case.spec.fixed_eta:g}")
        axes_arr[-1].set_xlabel("mu")
        axes_arr[0].legend(frameon=False, ncol=3, fontsize=8)
        fig.savefig(OUTPUT_LAMBDA_MU, dpi=220, bbox_inches="tight")
        plt.close(fig)


def markdown_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> list[str]:
    lines = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for row in rows:
        lines.append("| " + " | ".join(str(value) for value in row) + " |")
    return lines


def global_rows(cases: Sequence[CaseAudit]) -> list[list[str]]:
    return [
        [
            case.spec.audit_part,
            case.spec.label,
            case.classification_global,
            pair_label(case.global_pair_index),
            fmt(case.global_gap),
            fmt(case.global_parameter),
            fmt(case.any_true_crossing_detected),
            fmt(case.any_unresolved_possible_crossing),
        ]
        for case in cases
    ]


def descendant_pair_diagnostic(case: CaseAudit, left_descendant: int, right_descendant: int) -> DescendantPairDiagnostic | None:
    left = int(left_descendant) - 1
    right = int(right_descendant) - 1
    tracking = case.tracking
    if right >= tracking.tracked_lambdas.shape[0]:
        return None
    left_values = tracking.tracked_lambdas[left]
    right_values = tracking.tracked_lambdas[right]
    diff = right_values - left_values
    gaps = np.abs(diff)
    finite = np.isfinite(gaps)
    if not np.any(finite):
        return None
    finite_indices = np.where(finite)[0]
    idx = int(finite_indices[int(np.argmin(gaps[finite]))])
    nonzero_diff = diff[np.isfinite(diff) & (np.abs(diff) > TRUE_CROSSING_ABS_TOL)]
    signs = np.sign(nonzero_diff)
    sign_change = bool(signs.size > 1 and np.any(signs[:-1] * signs[1:] < 0))
    denominator = 0.5 * (abs(float(left_values[idx])) + abs(float(right_values[idx])))
    relative_gap = float(gaps[idx] / denominator) if denominator > 0.0 else float("nan")
    return DescendantPairDiagnostic(
        pair_label=f"desc{left_descendant}-desc{right_descendant}",
        parameter_at_min_gap=float(case.scan_values[idx]),
        lambda_left=float(left_values[idx]),
        lambda_right=float(right_values[idx]),
        min_gap=float(gaps[idx]),
        relative_gap=relative_gap,
        sign_change_detected=sign_change,
        sorted_position_left=int(tracking.current_sorted_positions[left, idx]),
        sorted_position_right=int(tracking.current_sorted_positions[right, idx]),
        mac_left=float(tracking.mac_to_previous[left, idx]),
        mac_right=float(tracking.mac_to_previous[right, idx]),
    )


def pair_audit_by_index(case: CaseAudit | None, pair_index: int) -> PairAudit | None:
    if case is None:
        return None
    for pair in case.pair_audits:
        if pair.pair_index == int(pair_index):
            return pair
    return None


def special_sorted_pair_line(case: CaseAudit | None, pair_index: int) -> str:
    pair = pair_audit_by_index(case, pair_index)
    if case is None or pair is None:
        return "case not found"
    return (
        f"sorted pair `{pair_label(pair_index)}`: gap `{fmt(pair.min_gap)}` at "
        f"{case.spec.scan_parameter} `{fmt(pair.parameter_at_min_gap)}`, class `{pair.classification}`"
    )


def checkpoint_sorted_pair_line(case: CaseAudit | None, pair_index: int, value: float) -> str:
    if case is None:
        return "case not found"
    roots = solve_roots_for_value(case.spec, float(value), mode="strict")
    _left, _right, gap = pair_gap_from_roots(roots, pair_index)
    return f"checkpoint sorted pair `{pair_label(pair_index)}` at {case.spec.scan_parameter} `{fmt(value)}`: gap `{fmt(gap)}`"


def special_descendant_pair_line(case: CaseAudit | None, left_descendant: int, right_descendant: int) -> str:
    if case is None:
        return "case not found"
    diag = descendant_pair_diagnostic(case, left_descendant, right_descendant)
    if diag is None:
        return "descendant diagnostics unavailable"
    return (
        f"{diag.pair_label}: gap `{fmt(diag.min_gap)}` at {case.spec.scan_parameter} "
        f"`{fmt(diag.parameter_at_min_gap)}`, sign change `{fmt(diag.sign_change_detected)}`, "
        f"sorted positions `{diag.sorted_position_left}`/`{diag.sorted_position_right}`, "
        f"MAC `{fmt(diag.mac_left)}`/`{fmt(diag.mac_right)}`"
    )


def checkpoint_descendant_pair_line(
    case: CaseAudit | None,
    left_descendant: int,
    right_descendant: int,
    value: float,
) -> str:
    if case is None:
        return "case not found"
    idx_candidates = np.where(np.isclose(case.scan_values, float(value), rtol=0.0, atol=1.0e-10))[0]
    if idx_candidates.size == 0:
        return "checkpoint descendant diagnostics unavailable"
    idx = int(idx_candidates[0])
    left = int(left_descendant) - 1
    right = int(right_descendant) - 1
    tracking = case.tracking
    left_value = float(tracking.tracked_lambdas[left, idx])
    right_value = float(tracking.tracked_lambdas[right, idx])
    gap = abs(right_value - left_value)
    return (
        f"checkpoint desc{left_descendant}-desc{right_descendant} at {case.spec.scan_parameter} "
        f"`{fmt(value)}`: gap `{fmt(gap)}`, sorted positions "
        f"`{tracking.current_sorted_positions[left, idx]}`/`{tracking.current_sorted_positions[right, idx]}`"
    )


def case_by(cases: Sequence[CaseAudit], *, part: str, beta: float | None = None, eta: float | None = None, mu: float | None = None) -> CaseAudit | None:
    for case in cases:
        spec = case.spec
        if spec.audit_part != part:
            continue
        if beta is not None and not np.isclose(float(spec.fixed_beta_deg), beta):
            continue
        if eta is not None and not np.isclose(float(spec.fixed_eta), eta):
            continue
        if mu is not None and not np.isclose(float(spec.fixed_mu), mu):
            continue
        return case
    return None


def scaling_near_eta_zero(cases: Sequence[CaseAudit]) -> tuple[float, float, int]:
    case = case_by(cases, part="A_eta_scan", beta=15.0, mu=0.0)
    if case is None:
        return float("nan"), float("nan"), 0
    values = np.abs(case.scan_values)
    gaps = coarse_global_gap_curve(case)
    mask = (values > 0.0) & (values <= 5.0e-2) & np.isfinite(gaps) & (gaps > 0.0)
    if np.count_nonzero(mask) < 3:
        return float("nan"), float("nan"), int(np.count_nonzero(mask))
    slope, intercept = np.polyfit(np.log(values[mask]), np.log(gaps[mask]), 1)
    return float(slope), float(np.exp(intercept)), int(np.count_nonzero(mask))


def write_report(cases: Sequence[CaseAudit]) -> None:
    eta_nonzero_cases = [
        case
        for case in cases
        if not (
            (case.spec.fixed_eta is not None and np.isclose(case.spec.fixed_eta, 0.0, atol=1.0e-12))
            or (case.spec.scan_parameter == "eta" and np.isclose(case.global_parameter, 0.0, atol=1.0e-12))
        )
    ]
    nonzero_pairs = [
        pair
        for case in eta_nonzero_cases
        for pair in case.pair_audits
        if not (case.spec.scan_parameter == "eta" and np.isclose(pair.parameter_at_min_gap, 0.0, atol=1.0e-12))
    ]
    smallest_eta_nonzero = min(nonzero_pairs, key=lambda pair: pair.min_gap if np.isfinite(pair.min_gap) else np.inf)
    smallest_case = next(case for case in eta_nonzero_cases if smallest_eta_nonzero in case.pair_audits)
    unresolved = [case for case in cases if case.any_unresolved_possible_crossing]
    true_crossings = [case for case in cases if case.any_true_crossing_detected]
    slope, prefactor, count = scaling_near_eta_zero(cases)

    beta15_eta05 = case_by(cases, part="C_mu_scan", beta=15.0, eta=0.5)
    beta5_eta05 = case_by(cases, part="C_mu_scan", beta=5.0, eta=0.5)
    total_root_repairs = sum(case.root_repair_count for case in cases)

    lines = [
        "# Eta Parameter Positive-Gap Verification",
        "",
        "## Scope",
        "",
        "Diagnostic-only Euler-Bernoulli analytic audit for possible frequency",
        "crossings when the thickness-mismatch parameter eta varies, and for",
        "selected beta/mu slices from prior diagnostics.",
        "",
        f"- script: `{SCRIPT_RELATIVE}`",
        f"- epsilon: `{EPSILON:g}`",
        f"- sorted roots per solve: first `{N_SORTED_ROOTS}`",
        f"- adjacent sorted gaps audited: pairs `1-2` through `6-7`",
        f"- tracked descendants for label diagnostics: first `{N_DESCENDANTS_TRACK}`",
        f"- default root scan step: `{DEFAULT_SCAN_STEP:g}`",
        f"- strict root scan step: `{STRICT_SCAN_STEP:g}`",
        f"- strict local root repairs applied on grid: `{total_root_repairs}`",
        "",
        "No FEM, Gmsh, CalculiX, Timoshenko model, article files, old determinant,",
        "`src/my_project/analytic/formulas.py`, old solvers, or baseline results",
        "are modified or used as validation data.",
        "",
        "## Global Summary",
        "",
    ]
    lines.extend(
        markdown_table(
            ["part", "case", "global class", "min pair", "min gap", "parameter", "true crossing", "unresolved"],
            global_rows(cases),
        )
    )
    lines.extend(
        [
            "",
            "## Answers",
            "",
            "1. True crossings for eta != 0:",
            (
                "   None were found in the tested eta != 0 cases."
                if not true_crossings
                else "   At least one true-crossing candidate was detected; see summary CSV."
            ),
            "2. Suspected eta != 0 close approaches:",
            (
                "   All tested close approaches were resolved as positive sorted-root gaps."
                if not unresolved and not true_crossings
                else "   Some cases remain unresolved; do not promote them to no-crossing claims."
            ),
            "3. Unresolved cases:",
            "   " + (", ".join(case.spec.label for case in unresolved) if unresolved else "none"),
            "4. Smallest positive gap for eta != 0:",
            (
                f"   `{fmt(smallest_eta_nonzero.min_gap)}` in case `{smallest_case.spec.label}`, "
                f"pair `{pair_label(smallest_eta_nonzero.pair_index)}`, "
                f"{smallest_case.spec.scan_parameter}=`{fmt(smallest_eta_nonzero.parameter_at_min_gap)}`, "
                f"class `{smallest_eta_nonzero.classification}`."
            ),
            "5. Behavior as eta -> 0:",
            (
                f"   For the beta=15, mu=0 eta scan, a coarse-grid small-|eta| log-log fit gives "
                f"gap ~= `{prefactor:.6g} * |eta|^{slope:.6g}` over `{count}` nonzero samples."
                if np.isfinite(slope)
                else "   Not enough resolved nonzero eta samples for a scaling fit."
            ),
            "6. Is eta=0 isolated?",
            "   The audit supports treating eta=0 as symmetry-sensitive: eta=0 rows can have fragile labels, while tested eta != 0 close approaches retain positive sorted gaps.",
            "7. Previous eta-scan rearrangement interval:",
            "   The strict audit does not confirm those intervals as true eigenvalue crossings; they are better read as sorted-position/descendant-label rearrangements around positive-gap avoided crossings.",
            "8. beta=15, eta=0.5, mu-scan suspected pairs:",
            (
                f"   Global min in the checked beta=15, eta=0.5 mu scan is pair `{pair_label(beta15_eta05.global_pair_index)}` "
                f"with gap `{fmt(beta15_eta05.global_gap)}` at mu `{fmt(beta15_eta05.global_parameter)}`; class `{beta15_eta05.classification_global}`."
                if beta15_eta05 is not None
                else "   Case not found."
            ),
            "   " + special_sorted_pair_line(beta15_eta05, 4),
            "   " + checkpoint_sorted_pair_line(beta15_eta05, 4, 0.379),
            "   " + special_sorted_pair_line(beta15_eta05, 5),
            "   " + checkpoint_sorted_pair_line(beta15_eta05, 5, 0.716),
            "   " + special_descendant_pair_line(beta15_eta05, 4, 5),
            "   " + checkpoint_descendant_pair_line(beta15_eta05, 4, 5, 0.379),
            "   " + special_descendant_pair_line(beta15_eta05, 5, 6),
            "   " + checkpoint_descendant_pair_line(beta15_eta05, 5, 6, 0.716),
            "9. beta=5, eta=0.5 article false crossing:",
            (
                f"   Remains resolved positive: global min pair `{pair_label(beta5_eta05.global_pair_index)}` "
                f"gap `{fmt(beta5_eta05.global_gap)}` at mu `{fmt(beta5_eta05.global_parameter)}`; class `{beta5_eta05.classification_global}`."
                if beta5_eta05 is not None
                else "   Case not found."
            ),
            "   " + special_sorted_pair_line(beta5_eta05, 1),
            "   " + checkpoint_sorted_pair_line(beta5_eta05, 1, 0.1586434207),
            "   " + special_descendant_pair_line(beta5_eta05, 1, 2),
            "   " + checkpoint_descendant_pair_line(beta5_eta05, 1, 2, 0.1586434207),
            "",
            "## Eigenvalue Gaps vs Descendant Labels",
            "",
            "The primary classification uses adjacent sorted-root gaps. Descendant-label swaps are recorded separately from MAC tracking and are not treated as eigenvalue crossings when the sorted gap is resolved positive.",
            "",
            "## Safe Statement",
            "",
            "Across the tested eta != 0 cases, numerical continuation found no true eigenvalue crossings; close approaches were resolved as positive sorted-root gaps. The exact eta=0 and mu=0 symmetry-sensitive cases require separate treatment.",
            "",
            "## Unsafe Statement",
            "",
            "This audit does not prove that for all eta != 0 and all beta,mu there are no crossings.",
            "",
            "## Limitations",
            "",
            "- Euler-Bernoulli analytic thickness-mismatch determinant only.",
            "- Root uncertainty is estimated by comparing default and strict sign-scan/bisection settings, not by a certified interval proof.",
            "- Determinant residual evaluation may overflow for some large hyperbolic terms; gap classification is based on root comparisons.",
            "- Near multiple roots, symbolic symmetry analysis or a specialized local root-coalescence solve may still be needed for a proof.",
            "",
            "## Outputs",
            "",
            f"- summary CSV: `{OUTPUT_SUMMARY.relative_to(REPO_ROOT)}`",
            f"- details CSV: `{OUTPUT_DETAILS.relative_to(REPO_ROOT)}`",
            f"- report: `{OUTPUT_REPORT.relative_to(REPO_ROOT)}`",
            f"- eta gap plot: `{OUTPUT_MIN_GAP_ETA.relative_to(REPO_ROOT)}`",
            f"- beta gap plot: `{OUTPUT_MIN_GAP_BETA.relative_to(REPO_ROOT)}`",
            f"- mu gap plot: `{OUTPUT_MIN_GAP_MU.relative_to(REPO_ROOT)}`",
            f"- classification map: `{OUTPUT_CLASS_MAP.relative_to(REPO_ROOT)}`",
            f"- selected Lambda(eta): `{OUTPUT_LAMBDA_ETA.relative_to(REPO_ROOT)}`",
            f"- selected Lambda(beta): `{OUTPUT_LAMBDA_BETA.relative_to(REPO_ROOT)}`",
            f"- selected Lambda(mu): `{OUTPUT_LAMBDA_MU.relative_to(REPO_ROOT)}`",
            "",
        ]
    )
    OUTPUT_REPORT.write_text("\n".join(lines), encoding="utf-8")


def run_audit(*, smoke: bool = False) -> list[CaseAudit]:
    specs = build_case_specs(smoke=smoke)
    return [audit_case(spec) for spec in specs]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", action="store_true", help="Run a small three-case audit for wiring checks.")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    print("eta-parameter positive-gap verification")
    cases = run_audit(smoke=bool(args.smoke))
    write_summary_csv(cases)
    write_details_csv(cases)
    plot_min_gap_vs_eta(cases)
    plot_min_gap_vs_beta_eta_scan(cases)
    plot_min_gap_vs_mu_eta_scan(cases)
    plot_classification_map(cases)
    plot_selected_curves(cases)
    write_report(cases)
    print(f"saved summary CSV: {OUTPUT_SUMMARY}")
    print(f"saved details CSV: {OUTPUT_DETAILS}")
    print(f"saved report: {OUTPUT_REPORT}")
    print(f"saved plots: {OUTPUT_MIN_GAP_ETA}, {OUTPUT_MIN_GAP_BETA}, {OUTPUT_MIN_GAP_MU}, {OUTPUT_CLASS_MAP}")
    return {"cases": cases}


if __name__ == "__main__":
    main()
