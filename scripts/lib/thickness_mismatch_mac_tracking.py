from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import linear_sum_assignment

from my_project.analytic.formulas_thickness_mismatch import (
    assemble_clamped_coupled_matrix_eta,
    thickness_mismatch_factors,
)


NEAR_ZERO_NORM = 1e-14


@dataclass(frozen=True)
class ShapeMacTrackingResult:
    tracked: np.ndarray
    current_sorted_indices: np.ndarray
    rows: list[dict[str, float | int | str]]
    warning_rows: list[dict[str, float | int | str]]


@dataclass(frozen=True)
class FrequencyAssignment:
    branch_index: int
    root_index: int
    distance_from_previous: float
    second_nearest_distance: float

    @property
    def assignment_margin(self) -> float:
        return self.second_nearest_distance - self.distance_from_previous


def mac_value(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    numerator = abs(float(np.dot(a, b))) ** 2
    denominator = float(np.dot(a, a)) * float(np.dot(b, b))
    return numerator / denominator if denominator > NEAR_ZERO_NORM else 0.0


def analytic_null_vector_eta(
    Lambda: float,
    *,
    beta_rad: float,
    mu: float,
    epsilon: float,
    eta: float,
) -> np.ndarray:
    matrix = assemble_clamped_coupled_matrix_eta(float(Lambda), beta_rad, float(mu), float(epsilon), float(eta))
    _, _singular_values, vh = np.linalg.svd(matrix)
    coeff = vh[-1, :].astype(float)
    if coeff.size:
        pivot = int(np.argmax(np.abs(coeff)))
        if coeff[pivot] < 0.0:
            coeff = -coeff
    return coeff


def reconstruct_components_eta(
    Lambda: float,
    *,
    mu: float,
    eta: float,
    epsilon: float,
    coeff: np.ndarray,
    s_norm: np.ndarray,
) -> dict[str, np.ndarray]:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    xi = np.asarray(s_norm, dtype=float)

    z1 = float(Lambda) * (1.0 - float(mu)) * xi / np.sqrt(factors.tau1)
    theta1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu)) * xi
    w_left = A1 * (np.cos(z1) - np.cosh(z1)) + B1 * (np.sin(z1) - np.sinh(z1))
    u_left = P1 * np.sin(theta1)

    z2_external_to_joint = -float(Lambda) * (1.0 + float(mu)) * xi / np.sqrt(factors.tau2)
    theta2_external_to_joint = -float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu)) * xi
    w_right = A2 * (np.cos(z2_external_to_joint) - np.cosh(z2_external_to_joint)) + B2 * (
        np.sin(z2_external_to_joint) - np.sinh(z2_external_to_joint)
    )
    u_right = P2 * np.sin(theta2_external_to_joint)

    return {
        "u_left": u_left,
        "w_left": w_left,
        "u_right": u_right[::-1],
        "w_right": w_right[::-1],
    }


def vector_from_components(components: Mapping[str, np.ndarray]) -> np.ndarray:
    return np.concatenate(
        [
            np.asarray(components["u_left"], dtype=float),
            np.asarray(components["w_left"], dtype=float),
            np.asarray(components["u_right"], dtype=float),
            np.asarray(components["w_right"], dtype=float),
        ]
    )


def analytic_shape_vector_eta(
    Lambda: float,
    *,
    beta_rad: float,
    mu: float,
    epsilon: float,
    eta: float,
    s_norm: np.ndarray,
) -> np.ndarray:
    coeff = analytic_null_vector_eta(
        float(Lambda),
        beta_rad=beta_rad,
        mu=float(mu),
        epsilon=float(epsilon),
        eta=float(eta),
    )
    components = reconstruct_components_eta(
        float(Lambda),
        mu=float(mu),
        eta=float(eta),
        epsilon=float(epsilon),
        coeff=coeff,
        s_norm=s_norm,
    )
    return vector_from_components(components)


def analytic_shape_vectors_for_roots(
    roots: Sequence[float],
    *,
    beta_rad: float,
    mu: float,
    epsilon: float,
    eta: float,
    s_norm: np.ndarray,
) -> list[np.ndarray]:
    return [
        analytic_shape_vector_eta(
            float(root),
            beta_rad=beta_rad,
            mu=float(mu),
            epsilon=float(epsilon),
            eta=float(eta),
            s_norm=s_norm,
        )
        for root in roots
    ]


def unique_nearest_frequency_assignment(previous: np.ndarray, candidates: np.ndarray) -> list[FrequencyAssignment]:
    costs = np.abs(np.asarray(previous, dtype=float)[:, None] - np.asarray(candidates, dtype=float)[None, :])
    pairs: list[tuple[float, int, int]] = [
        (float(costs[row, col]), int(row), int(col))
        for row in range(costs.shape[0])
        for col in range(costs.shape[1])
    ]
    pairs.sort(key=lambda item: item[0])

    assigned_rows: set[int] = set()
    assigned_cols: set[int] = set()
    assignment: dict[int, int] = {}
    for _distance, row, col in pairs:
        if row in assigned_rows or col in assigned_cols:
            continue
        assignment[row] = col
        assigned_rows.add(row)
        assigned_cols.add(col)
        if len(assignment) == len(previous):
            break

    if len(assignment) != len(previous):
        raise RuntimeError("Could not assign all tracked branches by nearest frequency.")

    info: list[FrequencyAssignment] = []
    for row in range(len(previous)):
        col = assignment[row]
        sorted_distances = np.sort(costs[row])
        second = float(sorted_distances[1]) if len(sorted_distances) > 1 else np.inf
        info.append(
            FrequencyAssignment(
                branch_index=row + 1,
                root_index=col + 1,
                distance_from_previous=float(costs[row, col]),
                second_nearest_distance=second,
            )
        )
    return info


def mac_assignment(previous_vectors: Sequence[np.ndarray], candidate_vectors: Sequence[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mac = np.array(
        [
            [mac_value(previous_vector, candidate_vector) for candidate_vector in candidate_vectors]
            for previous_vector in previous_vectors
        ],
        dtype=float,
    )
    rows, cols = linear_sum_assignment(1.0 - mac)
    assignment = np.full(len(previous_vectors), -1, dtype=int)
    for row, col in zip(rows, cols, strict=True):
        assignment[int(row)] = int(col)
    if np.any(assignment < 0):
        raise RuntimeError("Could not assign all tracked branches by shape MAC.")
    return assignment, mac


def track_mu_branches_shape_mac(
    *,
    beta_rad: float,
    epsilon: float,
    eta: float,
    mu_values: Sequence[float],
    roots_by_mu: Mapping[float, np.ndarray],
    num_tracked_branches: int,
    num_shape_samples: int = 401,
    mac_warning_threshold: float = 0.8,
    mac_margin_warning_threshold: float = 0.05,
    max_sorted_position_jump: int = 1,
    max_mu_step_for_confidence: float | None = None,
    allow_frequency_fallback: bool = False,
) -> ShapeMacTrackingResult:
    mu_grid = np.asarray(mu_values, dtype=float)
    if mu_grid.ndim != 1 or len(mu_grid) < 1:
        raise ValueError("mu_values must be a non-empty one-dimensional sequence.")
    first_roots = np.asarray(roots_by_mu[float(mu_grid[0])], dtype=float)
    if len(first_roots) < int(num_tracked_branches):
        raise ValueError("The first root set does not contain enough roots to seed the requested branches.")

    s_norm = np.linspace(0.0, 1.0, int(num_shape_samples), dtype=float)
    shape_vectors_by_mu = {
        float(mu): analytic_shape_vectors_for_roots(
            np.asarray(roots_by_mu[float(mu)], dtype=float),
            beta_rad=float(beta_rad),
            mu=float(mu),
            epsilon=float(epsilon),
            eta=float(eta),
            s_norm=s_norm,
        )
        for mu in mu_grid
    }

    tracked = np.full((int(num_tracked_branches), len(mu_grid)), np.nan, dtype=float)
    current_sorted_indices = np.full((int(num_tracked_branches), len(mu_grid)), -1, dtype=int)
    tracked[:, 0] = first_roots[: int(num_tracked_branches)]
    current_sorted_indices[:, 0] = np.arange(1, int(num_tracked_branches) + 1, dtype=int)

    rows: list[dict[str, float | int | str]] = []
    warning_rows: list[dict[str, float | int | str]] = []
    previous_lambdas = tracked[:, 0].copy()
    previous_vectors = shape_vectors_by_mu[float(mu_grid[0])][: int(num_tracked_branches)]

    for branch_idx, value in enumerate(previous_lambdas, start=1):
        rows.append(
            {
                "eta": float(eta),
                "mu_prev": np.nan,
                "mu": float(mu_grid[0]),
                "branch_index_from_mu0": int(branch_idx),
                "Lambda_tracked": float(value),
                "mac_sorted_root_index": int(branch_idx),
                "raw_mac_sorted_root_index": int(branch_idx),
                "diagnostic_candidate_sorted_position": int(branch_idx),
                "diagnostic_candidate_Lambda": float(value),
                "nearest_sorted_root_index": int(branch_idx),
                "tracking_step_status": "seed_mu0",
                "mac_to_previous": np.nan,
                "accepted_mac_to_previous": np.nan,
                "diagnostic_candidate_mac_to_previous": np.nan,
                "raw_mac_to_previous": np.nan,
                "second_best_mac": np.nan,
                "frequency_distance_from_previous": np.nan,
                "frequency_second_nearest_distance": np.nan,
                "frequency_assignment_margin": np.nan,
                "frequency_mac_disagreement": "no",
                "used_frequency_fallback": "no",
                "mac_margin": np.nan,
                "raw_mac_margin": np.nan,
                "assigned_from_previous_sorted_position": int(branch_idx),
                "sorted_position_jump": 0,
                "raw_sorted_position_jump": 0,
                "diagnostic_candidate_sorted_position_jump": 0,
                "low_mac": "no",
                "low_margin": "no",
                "large_mu_step": "no",
                "blocked_by_unresolved_neighbor": "no",
                "unresolved_assignment": "no",
                "suspicious_assignment": "no",
                "requires_refined_check": "no",
            }
        )

    for col in range(1, len(mu_grid)):
        mu_prev = float(mu_grid[col - 1])
        mu = float(mu_grid[col])
        roots = np.asarray(roots_by_mu[mu], dtype=float)
        candidate_vectors = shape_vectors_by_mu[mu]
        if len(candidate_vectors) < int(num_tracked_branches):
            raise RuntimeError(f"Not enough shape candidates at mu={mu:g}, eta={float(eta):g}.")

        raw_mac_cols, mac = mac_assignment(previous_vectors, candidate_vectors)
        freq_assignments = unique_nearest_frequency_assignment(previous_lambdas, roots)
        fallback_to_frequency = bool(allow_frequency_fallback) and any(
            float(mac[row, int(raw_mac_cols[row])]) < float(mac_warning_threshold)
            for row in range(int(num_tracked_branches))
        )
        candidate_cols = (
            np.array([int(assignment.root_index) - 1 for assignment in freq_assignments], dtype=int)
            if fallback_to_frequency
            else raw_mac_cols
        )
        previous_positions = current_sorted_indices[:, col - 1].astype(int)
        second_best_values = np.full(int(num_tracked_branches), np.nan, dtype=float)
        candidate_macs = np.full(int(num_tracked_branches), np.nan, dtype=float)
        raw_best_macs = np.full(int(num_tracked_branches), np.nan, dtype=float)
        mac_margins = np.full(int(num_tracked_branches), np.nan, dtype=float)
        raw_mac_margins = np.full(int(num_tracked_branches), np.nan, dtype=float)
        candidate_jumps = np.zeros(int(num_tracked_branches), dtype=int)
        raw_jumps = np.zeros(int(num_tracked_branches), dtype=int)
        low_mac_flags = np.zeros(int(num_tracked_branches), dtype=bool)
        low_margin_flags = np.zeros(int(num_tracked_branches), dtype=bool)
        large_mu_step_flags = np.zeros(int(num_tracked_branches), dtype=bool)
        unresolved = np.zeros(int(num_tracked_branches), dtype=bool)

        for branch_row in range(int(num_tracked_branches)):
            candidate_col = int(candidate_cols[branch_row])
            raw_root_col = int(raw_mac_cols[branch_row])
            row_macs = np.sort(mac[branch_row])[::-1]
            second_best = float(row_macs[1]) if len(row_macs) > 1 else np.nan
            candidate_mac = float(mac[branch_row, candidate_col])
            raw_best_mac = float(mac[branch_row, raw_root_col])
            mac_margin = candidate_mac - second_best if np.isfinite(second_best) else np.nan
            raw_mac_margin = raw_best_mac - second_best if np.isfinite(second_best) else np.nan
            candidate_jump = int(candidate_col) + 1 - int(previous_positions[branch_row])
            raw_jump = int(raw_root_col) + 1 - int(previous_positions[branch_row])
            low_mac = candidate_mac < float(mac_warning_threshold)
            low_margin = np.isfinite(raw_mac_margin) and raw_mac_margin < float(mac_margin_warning_threshold)
            large_mu_step = (
                max_mu_step_for_confidence is not None
                and abs(mu - mu_prev) > float(max_mu_step_for_confidence) + 1e-12
            )

            second_best_values[branch_row] = second_best
            candidate_macs[branch_row] = candidate_mac
            raw_best_macs[branch_row] = raw_best_mac
            mac_margins[branch_row] = mac_margin
            raw_mac_margins[branch_row] = raw_mac_margin
            candidate_jumps[branch_row] = int(candidate_jump)
            raw_jumps[branch_row] = int(raw_jump)
            low_mac_flags[branch_row] = bool(low_mac)
            low_margin_flags[branch_row] = bool(low_margin)
            large_mu_step_flags[branch_row] = bool(large_mu_step)
            unresolved[branch_row] = (
                abs(candidate_jump) > int(max_sorted_position_jump)
                or abs(raw_jump) > int(max_sorted_position_jump)
                or low_mac
                or low_margin
                or large_mu_step
            )

        blocked_by_unresolved_neighbor = np.zeros(int(num_tracked_branches), dtype=bool)
        while True:
            retained_positions = {
                int(previous_positions[row])
                for row in range(int(num_tracked_branches))
                if bool(unresolved[row])
            }
            updated = unresolved.copy()
            for branch_row in range(int(num_tracked_branches)):
                if bool(updated[branch_row]):
                    continue
                candidate_position = int(candidate_cols[branch_row]) + 1
                if (
                    candidate_position in retained_positions
                    and candidate_position != int(previous_positions[branch_row])
                ):
                    updated[branch_row] = True
                    blocked_by_unresolved_neighbor[branch_row] = True
            if np.array_equal(updated, unresolved):
                break
            unresolved = updated

        accepted_cols = np.array(
            [
                int(previous_positions[row]) - 1 if bool(unresolved[row]) else int(candidate_cols[row])
                for row in range(int(num_tracked_branches))
            ],
            dtype=int,
        )
        current_lambdas = np.full(int(num_tracked_branches), np.nan, dtype=float)
        current_vectors: list[np.ndarray] = []

        for branch_row in range(int(num_tracked_branches)):
            candidate_col = int(candidate_cols[branch_row])
            raw_root_col = int(raw_mac_cols[branch_row])
            root_col = int(accepted_cols[branch_row])
            frequency = freq_assignments[branch_row]
            second_best = float(second_best_values[branch_row])
            candidate_mac = float(candidate_macs[branch_row])
            raw_best_mac = float(raw_best_macs[branch_row])
            mac_margin = float(mac_margins[branch_row])
            raw_mac_margin = float(raw_mac_margins[branch_row])
            previous_sorted_position = int(previous_positions[branch_row])
            diagnostic_candidate_jump = int(candidate_jumps[branch_row])
            raw_sorted_position_jump = int(raw_jumps[branch_row])
            low_mac = bool(low_mac_flags[branch_row])
            low_margin = bool(low_margin_flags[branch_row])
            large_mu_step = bool(large_mu_step_flags[branch_row])
            unresolved_assignment = bool(unresolved[branch_row])
            if root_col < 0 or root_col >= len(roots):
                raise RuntimeError(
                    "Cannot retain previous canonical sorted position "
                    f"{previous_sorted_position} at mu={mu:g}; only {len(roots)} roots are available."
                )
            sorted_position_jump = int(root_col) + 1 - previous_sorted_position
            accepted_mac = float(mac[branch_row, root_col])
            suspicious = unresolved_assignment
            disagreement = int(frequency.root_index) != raw_root_col + 1
            if fallback_to_frequency:
                status = "frequency_fallback_low_mac"
            elif unresolved_assignment:
                status = "unresolved_assignment"
            else:
                status = "mac_ok"
            row = {
                "eta": float(eta),
                "mu_prev": mu_prev,
                "mu": mu,
                "branch_index_from_mu0": int(branch_row) + 1,
                "Lambda_tracked": float(roots[root_col]),
                "mac_sorted_root_index": int(root_col) + 1,
                "raw_mac_sorted_root_index": int(raw_root_col) + 1,
                "diagnostic_candidate_sorted_position": int(candidate_col) + 1,
                "diagnostic_candidate_Lambda": float(roots[candidate_col]),
                "nearest_sorted_root_index": int(frequency.root_index),
                "nearest_sorted_Lambda": float(roots[int(frequency.root_index) - 1]),
                "tracking_step_status": status,
                "mac_to_previous": candidate_mac,
                "accepted_mac_to_previous": accepted_mac,
                "diagnostic_candidate_mac_to_previous": candidate_mac,
                "raw_mac_to_previous": raw_best_mac,
                "second_best_mac": second_best,
                "frequency_distance_from_previous": float(frequency.distance_from_previous),
                "frequency_second_nearest_distance": float(frequency.second_nearest_distance),
                "frequency_assignment_margin": float(frequency.assignment_margin),
                "frequency_mac_disagreement": "yes" if disagreement else "no",
                "used_frequency_fallback": "yes" if fallback_to_frequency else "no",
                "mac_margin": float(mac_margin),
                "raw_mac_margin": float(raw_mac_margin),
                "assigned_from_previous_sorted_position": previous_sorted_position,
                "sorted_position_jump": int(sorted_position_jump),
                "raw_sorted_position_jump": int(raw_sorted_position_jump),
                "diagnostic_candidate_sorted_position_jump": int(diagnostic_candidate_jump),
                "low_mac": "yes" if low_mac else "no",
                "low_margin": "yes" if low_margin else "no",
                "large_mu_step": "yes" if large_mu_step else "no",
                "blocked_by_unresolved_neighbor": "yes" if blocked_by_unresolved_neighbor[branch_row] else "no",
                "unresolved_assignment": "yes" if unresolved_assignment else "no",
                "suspicious_assignment": "yes" if suspicious else "no",
                "requires_refined_check": "yes" if suspicious else "no",
            }
            rows.append(row)
            if disagreement or suspicious:
                warning_rows.append(row)
            current_lambdas[branch_row] = float(roots[root_col])
            current_sorted_indices[branch_row, col] = int(root_col) + 1
            current_vectors.append(candidate_vectors[root_col])

        tracked[:, col] = current_lambdas
        previous_lambdas = current_lambdas
        previous_vectors = current_vectors

    return ShapeMacTrackingResult(
        tracked=tracked,
        current_sorted_indices=current_sorted_indices,
        rows=rows,
        warning_rows=warning_rows,
    )


__all__ = [
    "FrequencyAssignment",
    "ShapeMacTrackingResult",
    "analytic_null_vector_eta",
    "analytic_shape_vector_eta",
    "analytic_shape_vectors_for_roots",
    "mac_value",
    "reconstruct_components_eta",
    "track_mu_branches_shape_mac",
    "unique_nearest_frequency_assignment",
    "vector_from_components",
]
