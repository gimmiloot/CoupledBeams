from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
import re
from typing import Iterable, Sequence
import warnings

import numpy as np
from scipy.optimize import linear_sum_assignment, minimize_scalar

from my_project.analytic.formulas import assemble_clamped_coupled_matrix
from my_project.analytic.solvers import find_first_n_roots
from scripts.lib.analytic_coupled_rods_shapes import (
    NEAR_ZERO_NORM,
    analytic_null_vector,
    concatenate_components,
    normalize_components,
    reconstruct_analytic_components,
    unique_sorted_roots,
)


DEFAULT_N_TRACK = 12
DEFAULT_N_SOLVE = 20
DEFAULT_FREQ_WEIGHT = 0.03
DEFAULT_MAC_WARNING_THRESHOLD = 0.8
DEFAULT_BETA_STEPS = 31
DEFAULT_MU_STEPS = 91
DEFAULT_MAX_REFINEMENT_DEPTH = 8
DEFAULT_MIN_BETA_STEP = 1e-3
DEFAULT_MIN_MU_STEP = 1e-4
DEFAULT_SVD_CANDIDATE_THRESHOLD = 1e-5
DEFAULT_NUM_SAMPLES = 401
DEFAULT_LMIN = 0.2
DEFAULT_LMAX0 = 55.0
DEFAULT_SCAN_STEP = 0.02
DEFAULT_GROW_FACTOR = 1.35
DEFAULT_MAX_TRIES = 8
DEFAULT_BRANCH_PREFIX = "bending_desc"
RIGHT_COORDINATE_FOR_TRACKING = "joint-to-external"
SHAPE_METRICS = ("full", "transverse")
DEBUG_FIELDNAMES = [
    "epsilon",
    "beta",
    "mu",
    "step_type",
    "step_index",
    "branch_id",
    "base_sorted_index",
    "current_sorted_index",
    "lambda",
    "Lambda",
    "mac_to_previous",
    "relative_lambda_jump",
    "relative_gap_lower",
    "relative_gap_upper",
    "warning_flag",
    "numerical_warning",
]
FAILURE_FIELDNAMES = [
    "epsilon",
    "previous_beta",
    "previous_mu",
    "current_beta",
    "current_mu",
    "branch_id",
    "previous_current_sorted_index",
    "previous_lambda",
    "candidate_sorted_index",
    "candidate_lambda",
    "mac_to_previous",
    "relative_lambda_difference",
    "assignment_cost",
    "selected_candidate",
    "warning_flag",
    "numerical_warning",
]
DEBUG_DIR = Path(__file__).resolve().parents[2] / "results" / "debug"


@dataclass
class AnalyticModeState:
    epsilon: float
    beta: float
    mu: float
    current_sorted_index: int
    Lambda: float
    coeff: np.ndarray = field(repr=False, compare=False)
    components: dict[str, np.ndarray] = field(repr=False, compare=False)
    shape_vector: np.ndarray = field(repr=False, compare=False)
    smallest_singular_value: float
    singular_value_ratio: float
    numerical_warning: str = ""


@dataclass
class BranchPoint:
    epsilon: float
    beta: float
    mu: float
    branch_id: str
    base_sorted_index: int
    current_sorted_index: int
    Lambda: float
    mac_to_previous: float
    relative_lambda_jump: float
    relative_gap_lower: float
    relative_gap_upper: float
    warning_flag: str
    step_type: str
    step_index: int
    coeff: np.ndarray = field(repr=False, compare=False)
    components: dict[str, np.ndarray] = field(repr=False, compare=False)
    smallest_singular_value: float = field(default=np.nan)
    singular_value_ratio: float = field(default=np.nan)
    sorted_lambdas: tuple[float, ...] = field(default_factory=tuple, repr=False)
    numerical_warning: str = ""

    def __post_init__(self) -> None:
        assert_point_matches_current_sorted_lambda(self)

    def debug_row(self) -> dict[str, float | int | str]:
        return {
            "epsilon": float(self.epsilon),
            "beta": float(self.beta),
            "mu": float(self.mu),
            "step_type": self.step_type,
            "step_index": int(self.step_index),
            "branch_id": self.branch_id,
            "base_sorted_index": int(self.base_sorted_index),
            "current_sorted_index": int(self.current_sorted_index),
            "lambda": float(self.Lambda),
            "Lambda": float(self.Lambda),
            "mac_to_previous": float(self.mac_to_previous),
            "relative_lambda_jump": float(self.relative_lambda_jump),
            "relative_gap_lower": float(self.relative_gap_lower),
            "relative_gap_upper": float(self.relative_gap_upper),
            "warning_flag": self.warning_flag,
            "numerical_warning": self.numerical_warning,
        }


@dataclass
class BranchTrackingResult:
    points: list[BranchPoint]
    warnings: list[str]
    summary: dict[str, float | int | str]

    def points_for_branch(self, branch_id: str) -> list[BranchPoint]:
        return [point for point in self.points if point.branch_id == branch_id]

    def final_point(self, branch_id: str) -> BranchPoint:
        points = self.points_for_branch(branch_id)
        if not points:
            raise KeyError(f"Branch id was not tracked: {branch_id}")
        return points[-1]

    def point_at(self, branch_id: str, *, beta: float, mu: float, atol: float = 1e-10) -> BranchPoint:
        matches = [
            point
            for point in self.points_for_branch(branch_id)
            if abs(float(point.beta) - float(beta)) <= atol and abs(float(point.mu) - float(mu)) <= atol
        ]
        if not matches:
            raise KeyError(f"No point for {branch_id} at beta={beta:g}, mu={mu:g}.")
        return matches[-1]

    def lambda_grid(self, branch_ids: Sequence[str], mu_values: Sequence[float], *, beta: float) -> np.ndarray:
        out = np.full((len(branch_ids), len(mu_values)), np.nan, dtype=float)
        for row, branch_id in enumerate(branch_ids):
            for col, mu in enumerate(mu_values):
                out[row, col] = self.point_at(branch_id, beta=beta, mu=float(mu)).Lambda
        return out

    def debug_rows(self, branch_id: str | None = None) -> list[dict[str, float | int | str]]:
        points = self.points if branch_id is None else self.points_for_branch(branch_id)
        return [point.debug_row() for point in points]

    def write_debug_csv(self, path: str | Path, branch_id: str | None = None) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        rows = self.debug_rows(branch_id=branch_id)
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=DEBUG_FIELDNAMES)
            writer.writeheader()
            writer.writerows(rows)
        return output


@dataclass
class AssignmentDiagnostics:
    mac: np.ndarray = field(repr=False)
    relative_frequency_difference: np.ndarray = field(repr=False)
    cost: np.ndarray = field(repr=False)
    assignment: np.ndarray = field(repr=False)
    current_states: list[AnalyticModeState] = field(repr=False)
    sorted_lambdas: tuple[float, ...]


def assert_point_matches_current_sorted_lambda(point: BranchPoint, *, rtol: float = 1e-9, atol: float = 1e-9) -> None:
    if not point.sorted_lambdas:
        return
    idx = int(point.current_sorted_index) - 1
    if idx < 0 or idx >= len(point.sorted_lambdas):
        raise RuntimeError(
            f"{point.branch_id} has current_sorted_index={point.current_sorted_index}, "
            f"but only {len(point.sorted_lambdas)} sorted roots are available."
        )
    expected = float(point.sorted_lambdas[idx])
    if not np.isclose(float(point.Lambda), expected, rtol=rtol, atol=atol):
        raise RuntimeError(
            f"{point.branch_id} indexing inconsistency at beta={point.beta:g}, mu={point.mu:g}: "
            f"Lambda={point.Lambda:.12g} does not match sorted root "
            f"{point.current_sorted_index}={expected:.12g}."
        )


def branch_id_from_base_sorted_index(index: int, *, prefix: str = DEFAULT_BRANCH_PREFIX) -> str:
    return f"{prefix}_{int(index):02d}"


def base_sorted_index_from_branch_id(branch_id: str) -> int:
    match = re.search(r"(\d+)$", str(branch_id))
    if match is None:
        raise ValueError(f"Cannot infer base sorted index from branch id: {branch_id!r}")
    return int(match.group(1))


def validate_shape_metric(shape_metric: str) -> str:
    if shape_metric not in SHAPE_METRICS:
        raise ValueError(f"Unknown shape metric {shape_metric!r}; expected one of {SHAPE_METRICS}.")
    return shape_metric


def shape_vector_from_components(components: dict[str, np.ndarray], *, shape_metric: str) -> np.ndarray:
    validate_shape_metric(shape_metric)
    if shape_metric == "full":
        vector = concatenate_components(components)
    else:
        vector = np.concatenate(
            [
                np.asarray(components["w_left"], dtype=float),
                np.asarray(components["w_right"], dtype=float),
            ]
        )
    norm = float(np.linalg.norm(vector))
    if norm <= NEAR_ZERO_NORM:
        return np.full_like(vector, np.nan, dtype=float)
    return np.asarray(vector, dtype=float) / norm


def solve_sorted_lambdas(
    *,
    beta: float,
    mu: float,
    epsilon: float,
    n_roots: int,
    Lmin: float = DEFAULT_LMIN,
    Lmax0: float = DEFAULT_LMAX0,
    scan_step: float = DEFAULT_SCAN_STEP,
    grow_factor: float = DEFAULT_GROW_FACTOR,
    max_tries: int = DEFAULT_MAX_TRIES,
) -> list[float]:
    roots, _ = solve_sorted_lambdas_with_warning(
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        n_roots=n_roots,
        Lmin=Lmin,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=grow_factor,
        max_tries=max_tries,
    )
    return roots


def solve_sorted_lambdas_with_warning(
    *,
    beta: float,
    mu: float,
    epsilon: float,
    n_roots: int,
    Lmin: float = DEFAULT_LMIN,
    Lmax0: float = DEFAULT_LMAX0,
    scan_step: float = DEFAULT_SCAN_STEP,
    grow_factor: float = DEFAULT_GROW_FACTOR,
    max_tries: int = DEFAULT_MAX_TRIES,
) -> tuple[list[float], str]:
    with warnings.catch_warnings(record=True) as captured:
        warnings.simplefilter("always", RuntimeWarning)
        roots = find_first_n_roots(
            float(np.deg2rad(beta)),
            float(mu),
            float(epsilon),
            int(n_roots),
            Lmin=float(Lmin),
            Lmax0=float(Lmax0),
            scan_step=float(scan_step),
            grow_factor=float(grow_factor),
            max_tries=int(max_tries),
        )
    warning_messages = []
    for item in captured:
        message = str(item.message)
        if "overflow" in message or "invalid value" in message or "nan" in message.lower():
            warning_messages.append(message)
    warning_text = "; ".join(dict.fromkeys(warning_messages))
    return unique_sorted_roots(roots), warning_text


def mode_state_from_lambda(
    *,
    epsilon: float,
    beta: float,
    mu: float,
    current_sorted_index: int,
    Lambda: float,
    s_norm: np.ndarray,
    shape_metric: str,
    numerical_warning: str = "",
) -> AnalyticModeState:
    beta_rad = float(np.deg2rad(beta))
    matrix = assemble_clamped_coupled_matrix(float(Lambda), beta_rad, float(mu), float(epsilon))
    coeff, smallest_singular_value, singular_value_ratio = analytic_null_vector(matrix)
    components_raw = reconstruct_analytic_components(
        float(Lambda),
        mu_value=float(mu),
        epsilon=float(epsilon),
        coeff=coeff,
        s_norm=s_norm,
        right_coordinate=RIGHT_COORDINATE_FOR_TRACKING,
    )
    components, _ = normalize_components(components_raw, plot_kind="full", normalize="max-full")
    shape_vector = shape_vector_from_components(components, shape_metric=shape_metric)
    return AnalyticModeState(
        epsilon=float(epsilon),
        beta=float(beta),
        mu=float(mu),
        current_sorted_index=int(current_sorted_index),
        Lambda=float(Lambda),
        coeff=coeff,
        components=components,
        shape_vector=shape_vector,
        smallest_singular_value=float(smallest_singular_value),
        singular_value_ratio=float(singular_value_ratio),
        numerical_warning=str(numerical_warning),
    )


def solve_mode_states(
    *,
    epsilon: float,
    beta: float,
    mu: float,
    n_solve: int,
    s_norm: np.ndarray,
    shape_metric: str,
    Lmin: float = DEFAULT_LMIN,
    Lmax0: float = DEFAULT_LMAX0,
    scan_step: float = DEFAULT_SCAN_STEP,
    grow_factor: float = DEFAULT_GROW_FACTOR,
    max_tries: int = DEFAULT_MAX_TRIES,
) -> tuple[list[AnalyticModeState], tuple[float, ...]]:
    roots, numerical_warning = solve_sorted_lambdas_with_warning(
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        n_roots=n_solve,
        Lmin=Lmin,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=grow_factor,
        max_tries=max_tries,
    )
    states = [
        mode_state_from_lambda(
            epsilon=epsilon,
            beta=beta,
            mu=mu,
            current_sorted_index=idx + 1,
            Lambda=root,
            s_norm=s_norm,
            shape_metric=shape_metric,
            numerical_warning=numerical_warning,
        )
        for idx, root in enumerate(roots)
    ]
    return states, tuple(float(root) for root in roots)


def smallest_singular_value_at_lambda(*, Lambda: float, beta: float, mu: float, epsilon: float) -> float:
    with np.errstate(over="ignore", invalid="ignore"):
        matrix = assemble_clamped_coupled_matrix(float(Lambda), float(np.deg2rad(beta)), float(mu), float(epsilon))
    if not np.all(np.isfinite(matrix)):
        return float("inf")
    try:
        return float(np.linalg.svd(matrix, compute_uv=False)[-1])
    except np.linalg.LinAlgError:
        return float("inf")


def svd_candidate_bracket(center: float, sorted_lambdas: Sequence[float]) -> tuple[float, float]:
    roots = [float(value) for value in sorted_lambdas if np.isfinite(value)]
    lower = max([value for value in roots if value < float(center)], default=max(0.05, float(center) - 2.5))
    upper = min([value for value in roots if value > float(center)], default=float(center) + 2.5)
    margin = max(1e-4, 1e-5 * max(abs(float(center)), 1.0))
    left = max(0.05, float(lower) + margin)
    right = max(left + margin, float(upper) - margin)
    return left, right


def svd_refined_candidate(
    *,
    beta: float,
    mu: float,
    epsilon: float,
    center: float,
    sorted_lambdas: Sequence[float],
    threshold: float = DEFAULT_SVD_CANDIDATE_THRESHOLD,
) -> float | None:
    left, right = svd_candidate_bracket(float(center), sorted_lambdas)
    if not np.isfinite(left) or not np.isfinite(right) or right <= left:
        return None
    result = minimize_scalar(
        lambda value: smallest_singular_value_at_lambda(Lambda=float(value), beta=beta, mu=mu, epsilon=epsilon),
        bounds=(left, right),
        method="bounded",
        options={"xatol": 1e-11},
    )
    if not result.success or not np.isfinite(result.x) or not np.isfinite(result.fun):
        return None
    if float(result.fun) > float(threshold):
        return None
    candidate = float(result.x)
    if any(abs(candidate - float(root)) <= 1e-6 for root in sorted_lambdas):
        return None
    return candidate


def augment_mode_states_with_svd_candidates(
    *,
    previous_points: Sequence[BranchPoint],
    current_states: Sequence[AnalyticModeState],
    sorted_lambdas: Sequence[float],
    beta: float,
    mu: float,
    epsilon: float,
    s_norm: np.ndarray,
    shape_metric: str,
    required_branch_ids: Sequence[str] | None,
) -> tuple[list[AnalyticModeState], tuple[float, ...]]:
    required_set = None if required_branch_ids is None else {str(branch_id) for branch_id in required_branch_ids}
    centers = [
        point.Lambda
        for point in previous_points
        if required_set is None or point.branch_id in required_set
    ]
    augmented_roots = list(float(root) for root in sorted_lambdas)
    for center in centers:
        candidate = svd_refined_candidate(
            beta=beta,
            mu=mu,
            epsilon=epsilon,
            center=float(center),
            sorted_lambdas=augmented_roots,
        )
        if candidate is not None:
            augmented_roots.append(candidate)
            augmented_roots = unique_sorted_roots(augmented_roots)
    if len(augmented_roots) == len(sorted_lambdas):
        return list(current_states), tuple(float(root) for root in sorted_lambdas)
    numerical_warning = current_states[0].numerical_warning if current_states else ""
    states = [
        mode_state_from_lambda(
            epsilon=epsilon,
            beta=beta,
            mu=mu,
            current_sorted_index=idx + 1,
            Lambda=root,
            s_norm=s_norm,
            shape_metric=shape_metric,
            numerical_warning=numerical_warning,
        )
        for idx, root in enumerate(augmented_roots)
    ]
    return states, tuple(float(root) for root in augmented_roots)


def relative_gap(sorted_lambdas: Sequence[float], current_sorted_index: int, *, side: str) -> float:
    idx = int(current_sorted_index) - 1
    if idx < 0 or idx >= len(sorted_lambdas):
        return float("nan")
    Lambda = float(sorted_lambdas[idx])
    denominator = max(abs(Lambda), NEAR_ZERO_NORM)
    if side == "lower":
        if idx == 0:
            return float("nan")
        return float((Lambda - float(sorted_lambdas[idx - 1])) / denominator)
    if side == "upper":
        if idx + 1 >= len(sorted_lambdas):
            return float("nan")
        return float((float(sorted_lambdas[idx + 1]) - Lambda) / denominator)
    raise ValueError(f"Unknown gap side: {side}")


def warning_for_mac(mac: float, *, threshold: float) -> str:
    if not np.isfinite(mac):
        return "needs_review"
    if float(mac) < float(threshold):
        return "needs_review"
    return "ok"


def low_mac_warning_text(point: BranchPoint) -> str:
    numerical = f" numerical_warning={point.numerical_warning}" if point.numerical_warning else ""
    return (
        f"branch_id={point.branch_id} epsilon={point.epsilon:g} beta={point.beta:g} mu={point.mu:g} "
        f"current_sorted_index={point.current_sorted_index} Lambda={point.Lambda:.12g} "
        f"warning_flag={point.warning_flag} MAC={point.mac_to_previous:.6g}{numerical}"
    )


def make_branch_point(
    *,
    state: AnalyticModeState,
    branch_id: str,
    base_sorted_index: int,
    mac_to_previous: float,
    relative_lambda_jump: float,
    sorted_lambdas: Sequence[float],
    step_type: str,
    step_index: int,
    mac_warning_threshold: float,
) -> BranchPoint:
    return BranchPoint(
        epsilon=float(state.epsilon),
        beta=float(state.beta),
        mu=float(state.mu),
        branch_id=str(branch_id),
        base_sorted_index=int(base_sorted_index),
        current_sorted_index=int(state.current_sorted_index),
        Lambda=float(state.Lambda),
        mac_to_previous=float(mac_to_previous),
        relative_lambda_jump=float(relative_lambda_jump),
        relative_gap_lower=relative_gap(sorted_lambdas, state.current_sorted_index, side="lower"),
        relative_gap_upper=relative_gap(sorted_lambdas, state.current_sorted_index, side="upper"),
        warning_flag=warning_for_mac(mac_to_previous, threshold=mac_warning_threshold),
        step_type=str(step_type),
        step_index=int(step_index),
        coeff=state.coeff.copy(),
        components={key: np.asarray(value, dtype=float).copy() for key, value in state.components.items()},
        smallest_singular_value=float(state.smallest_singular_value),
        singular_value_ratio=float(state.singular_value_ratio),
        sorted_lambdas=tuple(float(value) for value in sorted_lambdas),
        numerical_warning=state.numerical_warning,
    )


def branch_ids_for_count(n_track: int, *, branch_prefix: str = DEFAULT_BRANCH_PREFIX) -> list[str]:
    return [branch_id_from_base_sorted_index(idx, prefix=branch_prefix) for idx in range(1, int(n_track) + 1)]


def initialize_base_points(
    *,
    epsilon: float,
    n_track: int,
    n_solve: int,
    s_norm: np.ndarray,
    shape_metric: str,
    branch_prefix: str,
    mac_warning_threshold: float,
    Lmin: float,
    Lmax0: float,
    scan_step: float,
    grow_factor: float,
    max_tries: int,
) -> tuple[list[BranchPoint], tuple[float, ...]]:
    states, sorted_lambdas = solve_mode_states(
        epsilon=epsilon,
        beta=0.0,
        mu=0.0,
        n_solve=n_solve,
        s_norm=s_norm,
        shape_metric=shape_metric,
        Lmin=Lmin,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=grow_factor,
        max_tries=max_tries,
    )
    if len(states) < n_track:
        raise RuntimeError(
            f"Only found {len(states)} analytic roots at beta=0, mu=0, epsilon={epsilon:g}; "
            f"cannot initialize {n_track} tracked branches."
        )
    points = [
        make_branch_point(
            state=states[idx],
            branch_id=branch_id_from_base_sorted_index(idx + 1, prefix=branch_prefix),
            base_sorted_index=idx + 1,
            mac_to_previous=1.0,
            relative_lambda_jump=0.0,
            sorted_lambdas=sorted_lambdas,
            step_type="base",
            step_index=0,
            mac_warning_threshold=mac_warning_threshold,
        )
        for idx in range(n_track)
    ]
    return points, sorted_lambdas


def mac_matrix(
    previous_points: Sequence[BranchPoint],
    current_states: Sequence[AnalyticModeState],
    *,
    shape_metric: str,
) -> np.ndarray:
    out = np.zeros((len(previous_points), len(current_states)), dtype=float)
    for row, previous in enumerate(previous_points):
        prev_vector = shape_vector_from_components(previous.components, shape_metric=shape_metric)
        for col, current in enumerate(current_states):
            denominator = float(np.linalg.norm(prev_vector) * np.linalg.norm(current.shape_vector))
            if denominator <= NEAR_ZERO_NORM:
                out[row, col] = np.nan
            else:
                out[row, col] = abs(float(np.dot(prev_vector, current.shape_vector) / denominator))
    return out


def signed_shape_dot(previous: BranchPoint, current: AnalyticModeState, *, shape_metric: str) -> float:
    prev_vector = shape_vector_from_components(previous.components, shape_metric=shape_metric)
    denominator = float(np.linalg.norm(prev_vector) * np.linalg.norm(current.shape_vector))
    if denominator <= NEAR_ZERO_NORM:
        return float("nan")
    return float(np.dot(prev_vector, current.shape_vector) / denominator)


def align_state_to_previous(state: AnalyticModeState, previous: BranchPoint, *, shape_metric: str) -> AnalyticModeState:
    dot = signed_shape_dot(previous, state, shape_metric=shape_metric)
    if not np.isfinite(dot) or dot >= 0.0:
        return state
    state.coeff = -state.coeff
    state.components = {key: -np.asarray(value, dtype=float) for key, value in state.components.items()}
    state.shape_vector = -state.shape_vector
    return state


def assignment_diagnostics(
    *,
    previous_points: Sequence[BranchPoint],
    current_states: Sequence[AnalyticModeState],
    sorted_lambdas: Sequence[float],
    freq_weight: float,
    shape_metric: str,
) -> AssignmentDiagnostics:
    if len(current_states) < len(previous_points):
        raise RuntimeError(
            f"Need at least {len(previous_points)} current roots for assignment; got {len(current_states)}."
        )
    mac = mac_matrix(previous_points, current_states, shape_metric=shape_metric)
    mac_for_cost = np.where(np.isfinite(mac), mac, 0.0)
    prev_lambdas = np.array([point.Lambda for point in previous_points], dtype=float)
    cur_lambdas = np.array([state.Lambda for state in current_states], dtype=float)
    relative_frequency_difference = np.abs(prev_lambdas[:, None] - cur_lambdas[None, :]) / np.maximum(
        np.abs(prev_lambdas[:, None]),
        1.0,
    )
    cost = (1.0 - mac_for_cost) + float(freq_weight) * relative_frequency_difference
    rows, cols = linear_sum_assignment(cost)
    assignment = np.full(len(previous_points), -1, dtype=int)
    for row, col in zip(rows, cols):
        assignment[int(row)] = int(col)
    if np.any(assignment < 0):
        raise RuntimeError("Failed to assign all analytic branches.")
    return AssignmentDiagnostics(
        mac=mac,
        relative_frequency_difference=relative_frequency_difference,
        cost=cost,
        assignment=assignment,
        current_states=list(current_states),
        sorted_lambdas=tuple(float(value) for value in sorted_lambdas),
    )


def assign_next_points(
    *,
    previous_points: Sequence[BranchPoint],
    current_states: Sequence[AnalyticModeState],
    sorted_lambdas: Sequence[float],
    step_type: str,
    step_index: int,
    freq_weight: float,
    mac_warning_threshold: float,
    shape_metric: str,
) -> list[BranchPoint]:
    diagnostics = assignment_diagnostics(
        previous_points=previous_points,
        current_states=current_states,
        sorted_lambdas=sorted_lambdas,
        freq_weight=freq_weight,
        shape_metric=shape_metric,
    )
    next_points: list[BranchPoint] = []
    for row, col in enumerate(diagnostics.assignment):
        previous = previous_points[row]
        state = align_state_to_previous(current_states[int(col)], previous, shape_metric=shape_metric)
        rel_jump = abs(float(state.Lambda) - float(previous.Lambda)) / max(abs(float(previous.Lambda)), NEAR_ZERO_NORM)
        next_points.append(
            make_branch_point(
                state=state,
                branch_id=previous.branch_id,
                base_sorted_index=previous.base_sorted_index,
                mac_to_previous=float(diagnostics.mac[row, int(col)]),
                relative_lambda_jump=float(rel_jump),
                sorted_lambdas=sorted_lambdas,
                step_type=step_type,
                step_index=step_index,
                mac_warning_threshold=mac_warning_threshold,
            )
        )
    return next_points


def build_next_points_from_assignment(
    *,
    source_points: Sequence[BranchPoint],
    current_states: Sequence[AnalyticModeState],
    sorted_lambdas: Sequence[float],
    diagnostics: AssignmentDiagnostics,
    step_type: str,
    step_index: int,
    mac_warning_threshold: float,
    shape_metric: str,
) -> list[BranchPoint]:
    next_points: list[BranchPoint] = []
    for row, col in enumerate(diagnostics.assignment):
        previous = source_points[row]
        state = align_state_to_previous(current_states[int(col)], previous, shape_metric=shape_metric)
        rel_jump = abs(float(state.Lambda) - float(previous.Lambda)) / max(abs(float(previous.Lambda)), NEAR_ZERO_NORM)
        next_points.append(
            make_branch_point(
                state=state,
                branch_id=previous.branch_id,
                base_sorted_index=previous.base_sorted_index,
                mac_to_previous=float(diagnostics.mac[row, int(col)]),
                relative_lambda_jump=float(rel_jump),
                sorted_lambdas=sorted_lambdas,
                step_type=step_type,
                step_index=step_index,
                mac_warning_threshold=mac_warning_threshold,
            )
        )
    return next_points


def path_summary(points: Sequence[BranchPoint], warnings: Sequence[str]) -> dict[str, float | int | str]:
    nonbase_macs = [point.mac_to_previous for point in points if point.step_type != "base" and np.isfinite(point.mac_to_previous)]
    return {
        "point_count": int(len(points)),
        "warning_count": int(len(warnings)),
        "min_mac_to_previous": float(min(nonbase_macs)) if nonbase_macs else np.nan,
        "mean_mac_to_previous": float(np.mean(nonbase_macs)) if nonbase_macs else np.nan,
    }


def branch_row_index(previous_points: Sequence[BranchPoint], branch_id: str) -> int:
    for idx, point in enumerate(previous_points):
        if point.branch_id == branch_id:
            return idx
    raise KeyError(f"Branch id was not found in previous points: {branch_id}")


def blocking_low_mac_points(
    points: Sequence[BranchPoint],
    *,
    required_branch_ids: Sequence[str] | None,
) -> list[BranchPoint]:
    required_set = None if required_branch_ids is None else {str(branch_id) for branch_id in required_branch_ids}
    return [
        point
        for point in points
        if point.warning_flag != "ok" and (required_set is None or point.branch_id in required_set)
    ]


def failure_diagnostic_rows(
    *,
    previous_points: Sequence[BranchPoint],
    current_beta: float,
    current_mu: float,
    diagnostics: AssignmentDiagnostics,
    branch_id: str,
    mac_warning_threshold: float,
) -> list[dict[str, float | int | str]]:
    row = branch_row_index(previous_points, branch_id)
    previous = previous_points[row]
    selected_col = int(diagnostics.assignment[row])
    rows: list[dict[str, float | int | str]] = []
    for col, state in enumerate(diagnostics.current_states):
        mac = float(diagnostics.mac[row, col])
        rel_diff = float(diagnostics.relative_frequency_difference[row, col])
        rows.append(
            {
                "epsilon": float(previous.epsilon),
                "previous_beta": float(previous.beta),
                "previous_mu": float(previous.mu),
                "current_beta": float(current_beta),
                "current_mu": float(current_mu),
                "branch_id": str(branch_id),
                "previous_current_sorted_index": int(previous.current_sorted_index),
                "previous_lambda": float(previous.Lambda),
                "candidate_sorted_index": int(state.current_sorted_index),
                "candidate_lambda": float(state.Lambda),
                "mac_to_previous": mac,
                "relative_lambda_difference": rel_diff,
                "assignment_cost": float(diagnostics.cost[row, col]),
                "selected_candidate": int(col == selected_col),
                "warning_flag": warning_for_mac(mac, threshold=mac_warning_threshold),
                "numerical_warning": str(previous.numerical_warning or state.numerical_warning),
            }
        )
    return rows


def write_failure_diagnostic_csv(
    *,
    previous_points: Sequence[BranchPoint],
    current_beta: float,
    current_mu: float,
    diagnostics: AssignmentDiagnostics,
    branch_id: str,
    mac_warning_threshold: float,
) -> Path | None:
    rows = failure_diagnostic_rows(
        previous_points=previous_points,
        current_beta=current_beta,
        current_mu=current_mu,
        diagnostics=diagnostics,
        branch_id=branch_id,
        mac_warning_threshold=mac_warning_threshold,
    )
    if not rows:
        return None
    first = rows[0]
    path = DEBUG_DIR / (
        f"analytic_tracking_failure_{branch_id}"
        f"_eps{filename_number_token(float(first['epsilon']))}"
        f"_prevbeta{filename_number_token(float(first['previous_beta']))}"
        f"_prevmu{filename_number_token(float(first['previous_mu']))}"
        f"_curbeta{filename_number_token(float(current_beta))}"
        f"_curmu{filename_number_token(float(current_mu))}.csv"
    )
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=FAILURE_FIELDNAMES)
            writer.writeheader()
            writer.writerows(rows)
    except OSError:
        return None
    return path


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace("+", "").replace(".", "p")


def low_mac_failure_message(
    *,
    previous_points: Sequence[BranchPoint],
    current_beta: float,
    current_mu: float,
    diagnostics: AssignmentDiagnostics,
    branch_id: str,
    selected_point: BranchPoint,
    debug_csv: Path | None,
) -> str:
    row = branch_row_index(previous_points, branch_id)
    previous = previous_points[row]
    mac_row = diagnostics.mac[row]
    freq_row = diagnostics.relative_frequency_difference[row]
    selected_col = int(diagnostics.assignment[row])
    finite_mac = np.where(np.isfinite(mac_row), mac_row, -np.inf)
    best_mac_col = int(np.argmax(finite_mac))
    finite_freq = np.where(np.isfinite(freq_row), freq_row, np.inf)
    nearest_freq_col = int(np.argmin(finite_freq))
    best_mac_state = diagnostics.current_states[best_mac_col]
    nearest_freq_state = diagnostics.current_states[nearest_freq_col]
    selected_state = diagnostics.current_states[selected_col]
    debug_text = str(debug_csv) if debug_csv is not None else "not saved"
    numerical_warning = previous.numerical_warning or selected_state.numerical_warning or "none"
    return (
        "Low-MAC analytic branch assignment is not canonical after adaptive refinement. "
        f"branch_id={branch_id}; previous beta={previous.beta:g}, mu={previous.mu:g}, "
        f"Lambda={previous.Lambda:.12g}, current_sorted_index={previous.current_sorted_index}; "
        f"current beta={current_beta:g}, mu={current_mu:g}; "
        f"selected candidate index={selected_state.current_sorted_index}, Lambda={selected_state.Lambda:.12g}, "
        f"MAC={selected_point.mac_to_previous:.6g}, assignment_cost={diagnostics.cost[row, selected_col]:.6g}; "
        f"best-MAC candidate index={best_mac_state.current_sorted_index}, Lambda={best_mac_state.Lambda:.12g}, "
        f"MAC={diagnostics.mac[row, best_mac_col]:.6g}; "
        f"nearest-frequency candidate index={nearest_freq_state.current_sorted_index}, "
        f"Lambda={nearest_freq_state.Lambda:.12g}, "
        f"relative_lambda_difference={diagnostics.relative_frequency_difference[row, nearest_freq_col]:.6g}; "
        f"debug_csv={debug_text}. "
        f"numerical_warning={numerical_warning}. "
        "Recommendation: increase --max-refinement-depth, reduce --min-beta-step/--min-mu-step, "
        "or use --allow-low-mac for exploratory diagnostics only."
    )


def can_refine_transition(
    previous_points: Sequence[BranchPoint],
    *,
    current_beta: float,
    current_mu: float,
    depth: int,
    max_refinement_depth: int,
    min_beta_step: float,
    min_mu_step: float,
) -> bool:
    if depth >= int(max_refinement_depth):
        return False
    previous = previous_points[0]
    beta_step = abs(float(current_beta) - float(previous.beta))
    mu_step = abs(float(current_mu) - float(previous.mu))
    beta_refinable = beta_step > float(min_beta_step)
    mu_refinable = mu_step > float(min_mu_step)
    return beta_refinable or mu_refinable


def midpoint_transition(previous_points: Sequence[BranchPoint], *, current_beta: float, current_mu: float) -> tuple[float, float]:
    previous = previous_points[0]
    return 0.5 * (float(previous.beta) + float(current_beta)), 0.5 * (float(previous.mu) + float(current_mu))


def track_path(
    *,
    epsilon: float,
    path: Sequence[tuple[float, float, str, int]],
    n_track: int = DEFAULT_N_TRACK,
    n_solve: int = DEFAULT_N_SOLVE,
    freq_weight: float = DEFAULT_FREQ_WEIGHT,
    mac_warning_threshold: float = DEFAULT_MAC_WARNING_THRESHOLD,
    shape_metric: str = "full",
    branch_prefix: str = DEFAULT_BRANCH_PREFIX,
    num_samples: int = DEFAULT_NUM_SAMPLES,
    Lmin: float = DEFAULT_LMIN,
    Lmax0: float = DEFAULT_LMAX0,
    scan_step: float = DEFAULT_SCAN_STEP,
    grow_factor: float = DEFAULT_GROW_FACTOR,
    max_tries: int = DEFAULT_MAX_TRIES,
    allow_low_mac: bool = False,
    required_branch_ids: Sequence[str] | None = None,
    max_refinement_depth: int = DEFAULT_MAX_REFINEMENT_DEPTH,
    min_beta_step: float = DEFAULT_MIN_BETA_STEP,
    min_mu_step: float = DEFAULT_MIN_MU_STEP,
) -> BranchTrackingResult:
    validate_shape_metric(shape_metric)
    if n_track < 1:
        raise ValueError("n_track must be positive.")
    if n_solve < n_track:
        raise ValueError("n_solve must be at least n_track.")
    if not path:
        raise ValueError("Tracking path must contain at least one point.")
    first_beta, first_mu, _, _ = path[0]
    if abs(float(first_beta)) > 1e-12 or abs(float(first_mu)) > 1e-12:
        raise ValueError("Analytic branch tracking must start at beta=0, mu=0.")

    s_norm = np.linspace(0.0, 1.0, int(num_samples))
    base_points, _ = initialize_base_points(
        epsilon=epsilon,
        n_track=n_track,
        n_solve=n_solve,
        s_norm=s_norm,
        shape_metric=shape_metric,
        branch_prefix=branch_prefix,
        mac_warning_threshold=mac_warning_threshold,
        Lmin=Lmin,
        Lmax0=Lmax0,
        scan_step=scan_step,
        grow_factor=grow_factor,
        max_tries=max_tries,
    )
    all_points = list(base_points)
    previous_points = base_points

    def solve_and_assign(
        *,
        source_points: Sequence[BranchPoint],
        beta: float,
        mu: float,
        step_type: str,
        step_index: int,
    ) -> tuple[list[BranchPoint], AssignmentDiagnostics]:
        current_states, sorted_lambdas = solve_mode_states(
            epsilon=epsilon,
            beta=float(beta),
            mu=float(mu),
            n_solve=n_solve,
            s_norm=s_norm,
            shape_metric=shape_metric,
            Lmin=Lmin,
            Lmax0=Lmax0,
            scan_step=scan_step,
            grow_factor=grow_factor,
            max_tries=max_tries,
        )
        diagnostics = assignment_diagnostics(
            previous_points=source_points,
            current_states=current_states,
            sorted_lambdas=sorted_lambdas,
            freq_weight=freq_weight,
            shape_metric=shape_metric,
        )
        provisional_points = build_next_points_from_assignment(
            source_points=source_points,
            current_states=current_states,
            sorted_lambdas=sorted_lambdas,
            diagnostics=diagnostics,
            step_type=step_type,
            step_index=step_index,
            mac_warning_threshold=mac_warning_threshold,
            shape_metric=shape_metric,
        )
        if blocking_low_mac_points(provisional_points, required_branch_ids=required_branch_ids) and not allow_low_mac:
            current_states, sorted_lambdas = augment_mode_states_with_svd_candidates(
                previous_points=source_points,
                current_states=current_states,
                sorted_lambdas=sorted_lambdas,
                beta=beta,
                mu=mu,
                epsilon=epsilon,
                s_norm=s_norm,
                shape_metric=shape_metric,
                required_branch_ids=required_branch_ids,
            )
            diagnostics = assignment_diagnostics(
                previous_points=source_points,
                current_states=current_states,
                sorted_lambdas=sorted_lambdas,
                freq_weight=freq_weight,
                shape_metric=shape_metric,
            )
        next_points = build_next_points_from_assignment(
            source_points=source_points,
            current_states=current_states,
            sorted_lambdas=sorted_lambdas,
            diagnostics=diagnostics,
            step_type=step_type,
            step_index=step_index,
            mac_warning_threshold=mac_warning_threshold,
            shape_metric=shape_metric,
        )
        return next_points, diagnostics

    def advance_transition(
        *,
        source_points: Sequence[BranchPoint],
        beta: float,
        mu: float,
        step_type: str,
        step_index: int,
        depth: int,
    ) -> list[list[BranchPoint]]:
        next_points, diagnostics = solve_and_assign(
            source_points=source_points,
            beta=beta,
            mu=mu,
            step_type=step_type,
            step_index=step_index,
        )
        blocking = blocking_low_mac_points(next_points, required_branch_ids=required_branch_ids)
        if not blocking or allow_low_mac:
            return [next_points]
        if can_refine_transition(
            source_points,
            current_beta=beta,
            current_mu=mu,
            depth=depth,
            max_refinement_depth=max_refinement_depth,
            min_beta_step=min_beta_step,
            min_mu_step=min_mu_step,
        ):
            mid_beta, mid_mu = midpoint_transition(source_points, current_beta=beta, current_mu=mu)
            first_segments = advance_transition(
                source_points=source_points,
                beta=mid_beta,
                mu=mid_mu,
                step_type=step_type,
                step_index=step_index,
                depth=depth + 1,
            )
            midpoint_points = first_segments[-1]
            second_segments = advance_transition(
                source_points=midpoint_points,
                beta=beta,
                mu=mu,
                step_type=step_type,
                step_index=step_index,
                depth=depth + 1,
            )
            return first_segments + second_segments

        selected = blocking[0]
        debug_csv = write_failure_diagnostic_csv(
            previous_points=source_points,
            current_beta=float(beta),
            current_mu=float(mu),
            diagnostics=diagnostics,
            branch_id=selected.branch_id,
            mac_warning_threshold=mac_warning_threshold,
        )
        raise RuntimeError(
            low_mac_failure_message(
                previous_points=source_points,
                current_beta=float(beta),
                current_mu=float(mu),
                diagnostics=diagnostics,
                branch_id=selected.branch_id,
                selected_point=selected,
                debug_csv=debug_csv,
            )
        )

    for beta, mu, step_type, step_index in path[1:]:
        segments = advance_transition(
            source_points=previous_points,
            beta=float(beta),
            mu=float(mu),
            step_type=step_type,
            step_index=int(step_index),
            depth=0,
        )
        for next_points in segments:
            all_points.extend(next_points)
            previous_points = next_points

    warnings = [low_mac_warning_text(point) for point in all_points if point.warning_flag != "ok"]
    required_set = None if required_branch_ids is None else {str(branch_id) for branch_id in required_branch_ids}
    blocking_warnings = [
        low_mac_warning_text(point)
        for point in all_points
        if point.warning_flag != "ok" and (required_set is None or point.branch_id in required_set)
    ]
    if blocking_warnings and not allow_low_mac:
        preview = "; ".join(blocking_warnings[:3])
        if len(blocking_warnings) > 3:
            preview += f"; ... ({len(blocking_warnings)} total)"
        raise RuntimeError(
            "Low-MAC analytic branch assignment is not canonical. "
            "Recommendation: try increasing --mu-steps or --beta-steps, "
            "or run with --allow-low-mac / allow_low_mac=True for exploratory diagnostics only. "
            f"Warnings: {preview}"
        )
    return BranchTrackingResult(points=all_points, warnings=warnings, summary=path_summary(all_points, warnings))


def beta_then_mu_path(*, target_beta: float, target_mu: float, beta_steps: int, mu_steps: int) -> list[tuple[float, float, str, int]]:
    if beta_steps < 2:
        raise ValueError("beta_steps must be at least 2.")
    if mu_steps < 2:
        raise ValueError("mu_steps must be at least 2.")
    path: list[tuple[float, float, str, int]] = [(0.0, 0.0, "base", 0)]
    if abs(float(target_beta)) > 1e-15:
        for step_index, beta in enumerate(np.linspace(0.0, float(target_beta), int(beta_steps))[1:], start=1):
            path.append((float(beta), 0.0, "beta", step_index))
    if abs(float(target_mu)) > 1e-15:
        for step_index, mu in enumerate(np.linspace(0.0, float(target_mu), int(mu_steps))[1:], start=1):
            path.append((float(target_beta), float(mu), "mu", step_index))
    return path


def dense_mu_values_for_targets(target_mus: Sequence[float], *, mu_steps: int) -> np.ndarray:
    targets = np.asarray(list(target_mus), dtype=float)
    if targets.size == 0:
        return np.array([0.0], dtype=float)
    if np.any(targets < -1e-12):
        chunks = [np.linspace(0.0, float(target), int(mu_steps)) for target in targets]
    else:
        max_mu = max(float(np.max(targets)), 0.0)
        chunks = [np.linspace(0.0, max_mu, int(mu_steps))]
    chunks.append(targets)
    chunks.append(np.array([0.0], dtype=float))
    return np.unique(np.round(np.concatenate(chunks), 12))


def mu_sweep_path(*, beta: float, mu_values: Sequence[float], beta_steps: int) -> list[tuple[float, float, str, int]]:
    if beta_steps < 2:
        raise ValueError("beta_steps must be at least 2.")
    mu_array = np.asarray(list(mu_values), dtype=float)
    if mu_array.size == 0:
        raise ValueError("mu_values must contain at least one value.")
    if np.any(mu_array < -1e-12):
        raise ValueError("track_mu_sweep currently expects non-negative mu values.")

    path: list[tuple[float, float, str, int]] = [(0.0, 0.0, "base", 0)]
    if abs(float(beta)) > 1e-15:
        for step_index, beta_value in enumerate(np.linspace(0.0, float(beta), int(beta_steps))[1:], start=1):
            path.append((float(beta_value), 0.0, "beta", step_index))
    for step_index, mu in enumerate(np.unique(np.round(np.concatenate(([0.0], mu_array)), 12))[1:], start=1):
        path.append((float(beta), float(mu), "mu", step_index))
    return path


def track_beta_then_mu(
    *,
    epsilon: float,
    target_beta: float,
    target_mu: float,
    n_track: int = DEFAULT_N_TRACK,
    n_solve: int = DEFAULT_N_SOLVE,
    freq_weight: float = DEFAULT_FREQ_WEIGHT,
    mac_warning_threshold: float = DEFAULT_MAC_WARNING_THRESHOLD,
    shape_metric: str = "full",
    branch_prefix: str = DEFAULT_BRANCH_PREFIX,
    beta_steps: int = DEFAULT_BETA_STEPS,
    mu_steps: int = DEFAULT_MU_STEPS,
    num_samples: int = DEFAULT_NUM_SAMPLES,
    allow_low_mac: bool = False,
    required_branch_ids: Sequence[str] | None = None,
    max_refinement_depth: int = DEFAULT_MAX_REFINEMENT_DEPTH,
    min_beta_step: float = DEFAULT_MIN_BETA_STEP,
    min_mu_step: float = DEFAULT_MIN_MU_STEP,
) -> BranchTrackingResult:
    return track_path(
        epsilon=epsilon,
        path=beta_then_mu_path(target_beta=target_beta, target_mu=target_mu, beta_steps=beta_steps, mu_steps=mu_steps),
        n_track=n_track,
        n_solve=n_solve,
        freq_weight=freq_weight,
        mac_warning_threshold=mac_warning_threshold,
        shape_metric=shape_metric,
        branch_prefix=branch_prefix,
        num_samples=num_samples,
        allow_low_mac=allow_low_mac,
        required_branch_ids=required_branch_ids,
        max_refinement_depth=max_refinement_depth,
        min_beta_step=min_beta_step,
        min_mu_step=min_mu_step,
    )


def track_mu_sweep(
    *,
    epsilon: float,
    beta: float,
    mu_values: Sequence[float],
    n_track: int = DEFAULT_N_TRACK,
    n_solve: int = DEFAULT_N_SOLVE,
    freq_weight: float = DEFAULT_FREQ_WEIGHT,
    mac_warning_threshold: float = DEFAULT_MAC_WARNING_THRESHOLD,
    shape_metric: str = "full",
    branch_prefix: str = DEFAULT_BRANCH_PREFIX,
    beta_steps: int = DEFAULT_BETA_STEPS,
    num_samples: int = DEFAULT_NUM_SAMPLES,
    allow_low_mac: bool = False,
    required_branch_ids: Sequence[str] | None = None,
    max_refinement_depth: int = DEFAULT_MAX_REFINEMENT_DEPTH,
    min_beta_step: float = DEFAULT_MIN_BETA_STEP,
    min_mu_step: float = DEFAULT_MIN_MU_STEP,
) -> BranchTrackingResult:
    return track_path(
        epsilon=epsilon,
        path=mu_sweep_path(beta=beta, mu_values=mu_values, beta_steps=beta_steps),
        n_track=n_track,
        n_solve=n_solve,
        freq_weight=freq_weight,
        mac_warning_threshold=mac_warning_threshold,
        shape_metric=shape_metric,
        branch_prefix=branch_prefix,
        num_samples=num_samples,
        allow_low_mac=allow_low_mac,
        required_branch_ids=required_branch_ids,
        max_refinement_depth=max_refinement_depth,
        min_beta_step=min_beta_step,
        min_mu_step=min_mu_step,
    )


def nearby_sorted_lambdas(point: BranchPoint, indices: Iterable[int]) -> dict[int, float]:
    out: dict[int, float] = {}
    for index in indices:
        idx = int(index) - 1
        if 0 <= idx < len(point.sorted_lambdas):
            out[int(index)] = float(point.sorted_lambdas[idx])
    return out


__all__ = [
    "BranchPoint",
    "BranchTrackingResult",
    "DEFAULT_BETA_STEPS",
    "DEFAULT_BRANCH_PREFIX",
    "DEFAULT_FREQ_WEIGHT",
    "DEFAULT_MAC_WARNING_THRESHOLD",
    "DEFAULT_MAX_REFINEMENT_DEPTH",
    "DEFAULT_MIN_BETA_STEP",
    "DEFAULT_MIN_MU_STEP",
    "DEFAULT_MU_STEPS",
    "DEFAULT_N_SOLVE",
    "DEFAULT_N_TRACK",
    "assert_point_matches_current_sorted_lambda",
    "branch_id_from_base_sorted_index",
    "base_sorted_index_from_branch_id",
    "beta_then_mu_path",
    "dense_mu_values_for_targets",
    "nearby_sorted_lambdas",
    "solve_sorted_lambdas",
    "track_beta_then_mu",
    "track_mu_sweep",
]
