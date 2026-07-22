from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import brentq, linear_sum_assignment, minimize_scalar

from scripts.lib import general_spectrum_completeness as complete


BRANCH_CONTINUATION_ALGORITHM_VERSION = "branch_informed_continuation_v1"


@dataclass(frozen=True)
class BlockMap:
    bending_rows: tuple[int, ...]
    bending_columns: tuple[int, ...]
    axial_rows: tuple[int, ...]
    axial_columns: tuple[int, ...]


# These maps only expose the exact beta=0 block structure of the unchanged
# production matrices.  They also document the two different coefficient
# orderings used by the EB and Timoshenko implementations.
BETA0_BLOCK_MAPS: dict[str, BlockMap] = {
    complete.MODEL_EB: BlockMap((0, 2, 3, 4), (0, 1, 2, 3), (1, 5), (4, 5)),
    complete.MODEL_TIMO: BlockMap((0, 2, 3, 4), (0, 1, 3, 4), (1, 5), (2, 5)),
}


SEED_REFINED_TO_NEW_ROOT = "seed_refined_to_new_root"
SEED_REFINED_TO_EXISTING_ROOT = "seed_refined_to_existing_root"
SEED_REJECTED_NONSTATIONARY = "seed_rejected_nonstationary"
SEED_REJECTED_BOUNDARY_MINIMUM = "seed_rejected_boundary_minimum"
SEED_REJECTED_INSUFFICIENT_SINGULARITY = "seed_rejected_insufficient_singularity"
SEED_REJECTED_INDEPENDENT_DISAGREEMENT = "seed_rejected_independent_disagreement"

_LEGACY_SEED_STATUS_ALIASES = {
    "seed_resolved_to_existing_root": SEED_REFINED_TO_EXISTING_ROOT,
    "seed_rejected_not_stationary": SEED_REJECTED_NONSTATIONARY,
    "seed_rejected_full_matrix_svd": SEED_REJECTED_INSUFFICIENT_SINGULARITY,
    "seed_rejected_unstable_refinement": SEED_REJECTED_INDEPENDENT_DISAGREEMENT,
}


@dataclass(frozen=True)
class ContinuationSettings:
    requested_roots: int = 12
    candidate_roots: int = 20
    verification_candidate_roots: int = 24
    lambda_min: float = 0.2
    lambda_max: float = 40.0
    parent_scan_step: float = 0.01
    guard_scan_step: float = 0.05
    seed_half_width: float = 0.055
    cluster_gap_absolute: float = 0.04
    cluster_gap_relative: float = 3.0e-3
    beta_initial_step_deg: float = 0.5
    beta_min_step_deg: float = 0.0125
    beta_max_step_deg: float = 5.0
    beta_growth_factor: float = 1.6
    beta_shrink_factor: float = 0.5
    sigma_accept: float = 5.0e-6
    sigma_ratio_accept: float = 5.0e-4
    stationary_tolerance: float = 2.0e-4
    independent_root_tolerance: float = 2.0e-4
    root_match_tolerance: float = 2.0e-4
    mac_accept: float = 0.25
    subspace_mac_accept: float = 0.5
    run_global_guard: bool = True
    allow_strict_fallback: bool = True

    def validate(self, *, k_max: int = 10) -> None:
        if self.requested_roots < k_max + 2:
            raise ValueError("requested_roots must be at least k_max + 2")
        if self.candidate_roots < self.requested_roots:
            raise ValueError("candidate_roots must be at least requested_roots")
        if self.verification_candidate_roots <= self.candidate_roots:
            raise ValueError("verification_candidate_roots must exceed candidate_roots")
        if not 0.0 < self.lambda_min < self.lambda_max:
            raise ValueError("lambda bounds are invalid")
        positive = (
            self.parent_scan_step,
            self.guard_scan_step,
            self.seed_half_width,
            self.cluster_gap_absolute,
            self.cluster_gap_relative,
            self.beta_initial_step_deg,
            self.beta_min_step_deg,
            self.beta_max_step_deg,
            self.beta_growth_factor,
            self.beta_shrink_factor,
            self.sigma_accept,
            self.sigma_ratio_accept,
            self.stationary_tolerance,
            self.independent_root_tolerance,
            self.root_match_tolerance,
            self.mac_accept,
            self.subspace_mac_accept,
        )
        if any(not math.isfinite(float(value)) or value <= 0.0 for value in positive):
            raise ValueError("continuation settings must be finite and positive")
        if self.beta_min_step_deg > self.beta_initial_step_deg:
            raise ValueError("beta_min_step_deg must not exceed beta_initial_step_deg")
        if self.beta_initial_step_deg > self.beta_max_step_deg:
            raise ValueError("beta_initial_step_deg must not exceed beta_max_step_deg")
        if self.beta_growth_factor <= 1.0 or not 0.0 < self.beta_shrink_factor < 1.0:
            raise ValueError("invalid adaptive beta factors")
        if not 0.0 < self.mac_accept <= 1.0 or not 0.0 < self.subspace_mac_accept <= 1.0:
            raise ValueError("MAC thresholds must lie inside (0, 1]")

    def strict_settings(self) -> complete.SearchSettings:
        return complete.SearchSettings(
            requested_roots=self.requested_roots,
            candidate_roots=self.candidate_roots,
            verification_candidate_roots=self.verification_candidate_roots,
            lambda_min=self.lambda_min,
            lambda_max=self.lambda_max,
            scan_step=self.parent_scan_step,
            sigma_accept=self.sigma_accept,
            sigma_ratio_accept=self.sigma_ratio_accept,
            root_dedup_tol=self.root_match_tolerance,
            root_match_tol=self.root_match_tolerance,
            close_gap_threshold=self.cluster_gap_absolute,
            seed_half_width=self.seed_half_width,
        )


@dataclass
class BranchOperationCounts:
    parent_matrix_evaluations: int = 0
    parent_block_SVD_calls: int = 0
    local_matrix_evaluations: int = 0
    local_full_6x6_SVD_calls: int = 0
    local_scalar_projection_evaluations: int = 0
    local_reduced_matrix_SVD_calls: int = 0
    local_bounded_minimizations: int = 0
    guard_matrix_evaluations: int = 0
    guard_full_6x6_SVD_calls: int = 0
    guard_determinant_evaluations: int = 0
    guard_recovered_roots: int = 0
    strict_fallback_runs: int = 0
    strict_fallback_matrix_evaluations: int = 0
    strict_fallback_full_6x6_SVD_calls: int = 0
    force_verification_runs: int = 0
    force_verification_matrix_evaluations: int = 0
    force_verification_full_6x6_SVD_calls: int = 0
    beta_steps_attempted: int = 0
    beta_steps_accepted: int = 0
    beta_step_reductions: int = 0
    seed_windows_attempted: int = 0
    seed_windows_accepted: int = 0
    seed_windows_merged: int = 0
    seed_windows_rejected: int = 0
    cluster_steps_attempted: int = 0
    cluster_steps_resolved: int = 0
    cache_hits: int = 0
    cache_misses: int = 0

    def add(self, other: "BranchOperationCounts") -> None:
        for name in self.__dataclass_fields__:
            setattr(self, name, int(getattr(self, name)) + int(getattr(other, name)))


@dataclass(frozen=True)
class ParentBranch:
    model: str
    branch_id: str
    family: str
    family_index: int
    Lambda: float
    right_vector: tuple[float, ...]
    left_vector: tuple[float, ...]
    sigma_1: float
    sigma_ratio: float
    block_nullity: int
    full_nullity: int
    source: str
    off_block_max_abs: float = 0.0


@dataclass(frozen=True)
class RootRefinementResult:
    branch_id: str
    seed: float
    window_left: float
    window_right: float
    projected_candidate: float
    primary_candidate: float
    narrow_candidate: float
    independent_candidate: float
    Lambda: float
    sigma_1: float
    sigma_2: float
    sigma_ratio: float
    stationary: bool
    stable_refinement: bool
    independent_agreement: bool
    full_matrix_svd_pass: bool
    acceptance_status: str
    right_vector: tuple[float, ...]
    left_vector: tuple[float, ...]
    nullity: int
    notes: str = ""


@dataclass(frozen=True)
class ContinuedBranch:
    model: str
    branch_id: str
    parent_family: str
    beta_deg: float
    predicted_Lambda: float
    Lambda: float
    right_vector: tuple[float, ...]
    left_vector: tuple[float, ...]
    sigma_1: float
    sigma_ratio: float
    self_MAC: float
    nullity: int
    cluster_id: str
    refinement_status: str
    detection_source: str


@dataclass(frozen=True)
class BranchCluster:
    cluster_id: str
    beta_from_deg: float
    beta_to_deg: float
    branch_ids: tuple[str, ...]
    expected_dimension: int
    found_dimension: int
    subspace_MAC: float
    reduced_candidates: tuple[float, ...]
    resolved: bool
    status: str


@dataclass(frozen=True)
class ContinuationStep:
    model: str
    beta_from_deg: float
    beta_to_deg: float
    attempted_step_deg: float
    accepted: bool
    status: str
    branches: tuple[ContinuedBranch, ...]
    clusters: tuple[BranchCluster, ...]
    refinements: tuple[RootRefinementResult, ...]


@dataclass(frozen=True)
class GlobalGuardResult:
    guard_limit: float
    sign_candidates: tuple[float, ...]
    sigma_candidates: tuple[float, ...]
    unmatched_candidates: tuple[float, ...]
    unresolved_intervals: tuple[str, ...]
    passed: bool
    status: str


@dataclass(frozen=True)
class BranchContinuationResult:
    algorithm_version: str
    model: str
    geometry: complete.Geometry
    settings: ContinuationSettings
    parent_branches: tuple[ParentBranch, ...]
    branches: tuple[ContinuedBranch, ...]
    steps: tuple[ContinuationStep, ...]
    guard: GlobalGuardResult
    primary_vs_verification: tuple[dict[str, object], ...]
    strict_fallback_audit: tuple[dict[str, object], ...]
    k10_guard_resolved: bool
    full12_resolved: bool
    spectrum_status: str
    exclusion_reason: str
    oracle_agreement: bool | None
    force_verification_agreement: bool | None
    cache_status: str
    operations: BranchOperationCounts

    @property
    def values(self) -> tuple[float, ...]:
        return tuple(branch.Lambda for branch in sorted(self.branches, key=lambda item: item.Lambda))


def _canonical_svd(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]:
    scaled = complete.row_normalized(np.asarray(matrix, dtype=float))
    u, singular, vh = np.linalg.svd(scaled, full_matrices=True)
    right = np.asarray(vh[-1], dtype=float)
    left = np.asarray(u[:, -1], dtype=float)
    pivot = int(np.argmax(np.abs(right)))
    if right[pivot] < 0.0:
        right = -right
        left = -left
    sigma_1 = float(singular[-1])
    sigma_2 = float(singular[-2]) if singular.size >= 2 else float("inf")
    ratio = sigma_1 / sigma_2 if sigma_2 > 0.0 else float("inf")
    nullity = int(np.count_nonzero(singular <= 1.0e-10))
    return left, right, singular, sigma_1, float(ratio), max(1, nullity)


def _canonical_block_svd(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]:
    # Row-wise normalization is appropriate for the coupled 6x6 residual, but
    # it would divide a scalar axial equation by its own near-zero value and
    # erase the beta=0 block root.  A single block scale preserves exact zeros.
    raw = np.asarray(matrix, dtype=float)
    scale = float(np.linalg.norm(raw))
    scaled = raw / scale if math.isfinite(scale) and scale > 0.0 else raw
    u, singular, vh = np.linalg.svd(scaled, full_matrices=True)
    right = np.asarray(vh[-1], dtype=float)
    left = np.asarray(u[:, -1], dtype=float)
    pivot = int(np.argmax(np.abs(right)))
    if right[pivot] < 0.0:
        right = -right
        left = -left
    sigma_1 = float(singular[-1])
    sigma_2 = float(singular[-2]) if singular.size >= 2 else float("inf")
    ratio = sigma_1 / sigma_2 if sigma_2 > 0.0 else float("inf")
    nullity = int(np.count_nonzero(singular <= 1.0e-10))
    return left, right, singular, sigma_1, float(ratio), max(1, nullity)


def _subspace_mac(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> float:
    if not left or not right or len(left) != len(right):
        return float("nan")
    a = np.column_stack([np.asarray(value, dtype=float) for value in left])
    b = np.column_stack([np.asarray(value, dtype=float) for value in right])
    qa, _ = np.linalg.qr(a)
    qb, _ = np.linalg.qr(b)
    singular = np.linalg.svd(qa.T @ qb, compute_uv=False)
    return float(np.min(np.clip(singular, 0.0, 1.0) ** 2))


class _Evaluator:
    def __init__(
        self,
        provider: complete.MatrixProvider,
        operations: BranchOperationCounts,
        category: str,
    ) -> None:
        self.provider = provider
        self.operations = operations
        self.category = category
        self._matrix_cache: dict[float, np.ndarray] = {}
        self._svd_cache: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]] = {}
        self._det_cache: dict[float, float] = {}

    def matrix(self, value: float) -> np.ndarray:
        key = round(float(value), 13)
        if key not in self._matrix_cache:
            matrix = np.asarray(self.provider(float(value)), dtype=float)
            if matrix.shape != (6, 6) or not np.all(np.isfinite(matrix)):
                raise ValueError("nonfinite or non-6x6 characteristic matrix")
            self._matrix_cache[key] = matrix
            if self.category == "local":
                self.operations.local_matrix_evaluations += 1
            elif self.category == "guard":
                self.operations.guard_matrix_evaluations += 1
            else:
                self.operations.parent_matrix_evaluations += 1
        return self._matrix_cache[key]

    def svd(self, value: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]:
        key = round(float(value), 13)
        if key not in self._svd_cache:
            self._svd_cache[key] = _canonical_svd(self.matrix(value))
            if self.category == "local":
                self.operations.local_full_6x6_SVD_calls += 1
            elif self.category == "guard":
                self.operations.guard_full_6x6_SVD_calls += 1
        return self._svd_cache[key]

    def sigma(self, value: float) -> float:
        return self.svd(value)[3]

    def determinant(self, value: float) -> float:
        key = round(float(value), 13)
        if key not in self._det_cache:
            self._det_cache[key] = float(np.linalg.det(complete.row_normalized(self.matrix(value))))
            if self.category == "guard":
                self.operations.guard_determinant_evaluations += 1
        return self._det_cache[key]


def _off_block_max_abs(matrix: np.ndarray, block_map: BlockMap) -> float:
    allowed = np.zeros((6, 6), dtype=bool)
    allowed[np.ix_(block_map.bending_rows, block_map.bending_columns)] = True
    allowed[np.ix_(block_map.axial_rows, block_map.axial_columns)] = True
    outside = np.abs(np.asarray(matrix, dtype=float)[~allowed])
    return float(np.max(outside)) if outside.size else 0.0


def _factorized_eta0_parents(
    model: str,
    geometry: complete.Geometry,
    settings: ContinuationSettings,
    operations: BranchOperationCounts,
) -> tuple[ParentBranch, ...]:
    count = settings.candidate_roots
    if model == complete.MODEL_TIMO:
        factorized = complete.straight_factorized.factorized_straight_spectrum(
            geometry.epsilon_0, geometry.mu, count
        )
        family_values = [
            (
                "axial" if item.family == "axial" else "bending",
                item.family_index,
                item.value,
                f"{factorized.algorithm_version}:{item.detection_source}",
            )
            for item in factorized.roots[:count]
        ]
    else:
        axial = np.sqrt(
            np.arange(1, count + 1, dtype=float) * np.pi / (2.0 * geometry.epsilon_0)
        )
        bending = np.asarray(complete.fixed_fixed_eb_roots(count), dtype=float) / 2.0
        family_values = [
            *(('axial', index, float(value), 'exact_axial_formula') for index, value in enumerate(axial, start=1)),
            *(('bending', index, float(value), 'fixed_fixed_EB_bending_family') for index, value in enumerate(bending, start=1)),
        ]
        family_values = sorted(family_values, key=lambda item: (item[2], item[0], item[1]))[:count]
    beta0 = complete.Geometry(geometry.epsilon_0, 0.0, geometry.mu, 0.0)
    provider = complete.model_matrix_provider(model, beta0)
    block_map = BETA0_BLOCK_MAPS[model]
    roots: list[ParentBranch] = []
    for family, family_index, value, source in family_values:
        rows = block_map.axial_rows if family == "axial" else block_map.bending_rows
        columns = block_map.axial_columns if family == "axial" else block_map.bending_columns
        matrix = np.asarray(provider(float(value)), dtype=float)
        operations.parent_matrix_evaluations += 1
        left_block, right_block, _singular, _block_sigma, _block_ratio, block_nullity = _canonical_block_svd(
            matrix[np.ix_(rows, columns)]
        )
        operations.parent_block_SVD_calls += 1
        right = np.zeros(6, dtype=float)
        left = np.zeros(6, dtype=float)
        right[np.asarray(columns, dtype=int)] = right_block
        left[np.asarray(rows, dtype=int)] = left_block
        _full_left, _full_right, _full_singular, full_sigma, full_ratio, full_nullity = _canonical_svd(matrix)
        operations.parent_block_SVD_calls += 1
        roots.append(
            ParentBranch(
                model=model,
                branch_id=f"{family[0].upper()}{int(family_index):03d}",
                family=family,
                family_index=int(family_index),
                Lambda=float(value),
                right_vector=tuple(float(item) for item in right),
                left_vector=tuple(float(item) for item in left),
                sigma_1=full_sigma,
                sigma_ratio=full_ratio,
                block_nullity=block_nullity,
                full_nullity=full_nullity,
                source=source,
                off_block_max_abs=_off_block_max_abs(matrix, block_map),
            )
        )
    roots.sort(key=lambda item: (item.Lambda, item.family, item.family_index))
    return tuple(roots[:count])


def _block_spectrum(
    model: str,
    geometry: complete.Geometry,
    settings: ContinuationSettings,
    operations: BranchOperationCounts,
) -> tuple[ParentBranch, ...]:
    if abs(geometry.eta) <= 1.0e-14:
        return _factorized_eta0_parents(model, geometry, settings, operations)
    beta0 = complete.Geometry(geometry.epsilon_0, 0.0, geometry.mu, geometry.eta)
    provider = complete.model_matrix_provider(model, beta0)
    block_map = BETA0_BLOCK_MAPS[model]
    families = (
        ("bending", block_map.bending_rows, block_map.bending_columns),
        ("axial", block_map.axial_rows, block_map.axial_columns),
    )
    roots: list[ParentBranch] = []
    for family, rows, columns in families:
        cache: dict[float, tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]] = {}

        def block(value: float) -> np.ndarray:
            operations.parent_matrix_evaluations += 1
            full = np.asarray(provider(float(value)), dtype=float)
            return full[np.ix_(rows, columns)]

        def diagnostics(value: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, int]:
            key = round(float(value), 13)
            if key not in cache:
                cache[key] = _canonical_block_svd(block(value))
                operations.parent_block_SVD_calls += 1
            return cache[key]

        def sigma(value: float) -> float:
            return diagnostics(value)[3]

        def determinant(value: float) -> float:
            matrix = block(value)
            scale = float(np.linalg.norm(matrix))
            scaled = matrix / scale if math.isfinite(scale) and scale > 0.0 else matrix
            return float(np.linalg.det(scaled))

        requested_grid = np.arange(
            settings.lambda_min,
            settings.lambda_max + 0.5 * settings.parent_scan_step,
            settings.parent_scan_step,
            dtype=float,
        )
        valid_grid: list[float] = []
        det_values: list[float] = []
        for value in requested_grid:
            try:
                det_value = determinant(float(value))
            except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
                # The Timoshenko basis has a documented finite real-wave
                # domain for a given epsilon.  Roots below its first invalid
                # point remain valid; the invalid tail is not sampled.
                break
            valid_grid.append(float(value))
            det_values.append(det_value)
        grid = np.asarray(valid_grid, dtype=float)
        det = np.asarray(det_values, dtype=float)
        candidates: list[float] = []
        for index in range(len(grid) - 1):
            left, right = float(grid[index]), float(grid[index + 1])
            f_left, f_right = float(det[index]), float(det[index + 1])
            if not (math.isfinite(f_left) and math.isfinite(f_right)):
                continue
            if f_left == 0.0:
                candidates.append(left)
            elif f_left * f_right < 0.0:
                try:
                    candidates.append(float(brentq(determinant, left, right, xtol=1.0e-12, rtol=1.0e-12)))
                except (ValueError, RuntimeError):
                    pass
        coarse_step = max(2.0 * settings.parent_scan_step, 0.02)
        coarse_upper = float(grid[-1]) if grid.size else settings.lambda_min
        coarse_count = int(math.floor((coarse_upper - settings.lambda_min) / coarse_step)) + 1
        coarse = settings.lambda_min + coarse_step * np.arange(max(coarse_count, 1), dtype=float)
        sigmas = np.asarray([sigma(float(value)) for value in coarse], dtype=float)
        for index in range(1, len(coarse) - 1):
            if sigmas[index] <= sigmas[index - 1] and sigmas[index] <= sigmas[index + 1] and sigmas[index] <= 0.08:
                fit = minimize_scalar(
                    sigma,
                    bounds=(float(coarse[index - 1]), float(coarse[index + 1])),
                    method="bounded",
                    options={"xatol": 1.0e-12},
                )
                candidates.append(float(fit.x))
        accepted: list[float] = []
        for candidate in sorted(candidates):
            fit = minimize_scalar(
                sigma,
                bounds=(max(settings.lambda_min, candidate - 1.5 * settings.parent_scan_step),
                        min(settings.lambda_max, candidate + 1.5 * settings.parent_scan_step)),
                method="bounded",
                options={"xatol": 1.0e-13},
            )
            # A sign-bracketed axial root can be much sharper than the result
            # of a generic bounded minimum on its V-shaped singular value.
            value = min((float(candidate), float(fit.x)), key=sigma)
            _left, _right, _singular, sigma_1, ratio, _nullity = diagnostics(value)
            if sigma_1 > settings.sigma_accept or ratio > settings.sigma_ratio_accept:
                continue
            duplicate = next(
                (index for index, old in enumerate(accepted) if abs(value - old) <= settings.root_match_tolerance),
                None,
            )
            if duplicate is not None:
                if sigma(value) < sigma(accepted[duplicate]):
                    accepted[duplicate] = value
                continue
            accepted.append(value)
        for family_index, value in enumerate(accepted, start=1):
            left_block, right_block, singular, sigma_1, ratio, block_nullity = diagnostics(value)
            right = np.zeros(6, dtype=float)
            left = np.zeros(6, dtype=float)
            right[np.asarray(columns, dtype=int)] = right_block
            left[np.asarray(rows, dtype=int)] = left_block
            full_left, full_right, full_singular, full_sigma, full_ratio, full_nullity = _canonical_svd(provider(value))
            operations.parent_matrix_evaluations += 1
            operations.parent_block_SVD_calls += 1
            # At a cross-family degeneracy the full SVD basis is not a stable
            # branch identity; the embedded exact block vector is.
            del full_left, full_right, singular, full_singular
            roots.append(
                ParentBranch(
                    model=model,
                    branch_id=f"{family[0].upper()}{family_index:03d}",
                    family=family,
                    family_index=family_index,
                    Lambda=value,
                    right_vector=tuple(float(item) for item in right),
                    left_vector=tuple(float(item) for item in left),
                    sigma_1=full_sigma,
                    sigma_ratio=full_ratio,
                    block_nullity=block_nullity,
                    full_nullity=full_nullity,
                    source="exact_beta0_block_nullspace",
                    off_block_max_abs=_off_block_max_abs(provider(value), block_map),
                )
            )
    roots.sort(key=lambda item: (item.Lambda, item.family, item.family_index))
    if len(roots) < settings.candidate_roots:
        raise RuntimeError(f"beta=0 block solve found only {len(roots)} roots")
    return tuple(roots[: settings.candidate_roots])


def _parent_as_branches(parent: Sequence[ParentBranch]) -> tuple[ContinuedBranch, ...]:
    return tuple(
        ContinuedBranch(
            model=item.model,
            branch_id=item.branch_id,
            parent_family=item.family,
            beta_deg=0.0,
            predicted_Lambda=item.Lambda,
            Lambda=item.Lambda,
            right_vector=item.right_vector,
            left_vector=item.left_vector,
            sigma_1=item.sigma_1,
            sigma_ratio=item.sigma_ratio,
            self_MAC=1.0,
            nullity=item.full_nullity,
            cluster_id="",
            refinement_status=SEED_REFINED_TO_NEW_ROOT,
            detection_source="exact_beta0_block_solve",
        )
        for item in parent
    )


def _predict(branch: ContinuedBranch, earlier: ContinuedBranch | None, beta_to: float) -> float:
    if earlier is None or abs(branch.beta_deg - earlier.beta_deg) <= 1.0e-14:
        return branch.Lambda
    slope = (branch.Lambda - earlier.Lambda) / (branch.beta_deg - earlier.beta_deg)
    return float(branch.Lambda + slope * (beta_to - branch.beta_deg))


def _refine_seed(
    evaluator: _Evaluator,
    branch: ContinuedBranch,
    seed: float,
    settings: ContinuationSettings,
    operations: BranchOperationCounts,
    half_width: float | None = None,
) -> RootRefinementResult:
    operations.seed_windows_attempted += 1
    half = (
        max(float(half_width), 2.0 * settings.stationary_tolerance)
        if half_width is not None
        else max(settings.seed_half_width, 2.5 * abs(seed - branch.Lambda))
    )
    left = max(settings.lambda_min, seed - half)
    right = min(settings.lambda_max, seed + half)
    u0 = np.asarray(branch.left_vector, dtype=float)
    v0 = np.asarray(branch.right_vector, dtype=float)

    def projected_signed(value: float) -> float:
        operations.local_scalar_projection_evaluations += 1
        matrix = complete.row_normalized(evaluator.matrix(value))
        return float(u0 @ matrix @ v0)

    projected_fit = minimize_scalar(
        lambda value: abs(projected_signed(value)),
        bounds=(left, right),
        method="bounded",
        options={"xatol": 1.0e-11},
    )
    operations.local_bounded_minimizations += 1
    projected_candidate = float(projected_fit.x)
    search_half = min(half, max(0.006, 0.18 * half))
    search_left = max(left, projected_candidate - search_half)
    search_right = min(right, projected_candidate + search_half)
    local_grid = np.linspace(search_left, search_right, 13)
    determinant = np.asarray([evaluator.determinant(float(value)) for value in local_grid])
    bracketed: list[float] = []
    for index in range(len(local_grid) - 1):
        if not (math.isfinite(determinant[index]) and math.isfinite(determinant[index + 1])):
            continue
        if determinant[index] * determinant[index + 1] < 0.0:
            try:
                bracketed.append(
                    float(
                        brentq(
                            evaluator.determinant,
                            float(local_grid[index]),
                            float(local_grid[index + 1]),
                            xtol=1.0e-13,
                            rtol=1.0e-13,
                        )
                    )
                )
            except (ValueError, RuntimeError):
                pass
    sigma_fit = minimize_scalar(
        evaluator.sigma,
        bounds=(search_left, search_right),
        method="bounded",
        options={"xatol": 1.0e-12},
    )
    operations.local_bounded_minimizations += 1
    primary_options = [projected_candidate, float(sigma_fit.x), *bracketed]
    primary = min(primary_options, key=evaluator.sigma)
    narrow_half = max(0.25 * search_half, 8.0 * settings.stationary_tolerance)
    narrow_left = max(left, primary - narrow_half)
    narrow_right = min(right, primary + narrow_half)
    narrow_fit = minimize_scalar(
        evaluator.sigma, bounds=(narrow_left, narrow_right), method="bounded", options={"xatol": 2.5e-13}
    )
    operations.local_bounded_minimizations += 1
    independent_left = max(left, projected_candidate - 0.93 * search_half)
    independent_right = min(right, projected_candidate + 0.89 * search_half)
    if independent_right <= independent_left:
        independent_left, independent_right = left, right
    narrow_value = float(narrow_fit.x)
    value = min((primary, narrow_value), key=evaluator.sigma)
    independent_grid = np.linspace(independent_left, independent_right, 14)
    independent_det = np.asarray([evaluator.determinant(float(item)) for item in independent_grid])
    independent_bracketed: list[float] = []
    for index in range(len(independent_grid) - 1):
        if not (math.isfinite(independent_det[index]) and math.isfinite(independent_det[index + 1])):
            continue
        if independent_det[index] * independent_det[index + 1] < 0.0:
            try:
                independent_bracketed.append(
                    float(
                        brentq(
                            evaluator.determinant,
                            float(independent_grid[index]),
                            float(independent_grid[index + 1]),
                            xtol=2.0e-13,
                            rtol=2.0e-13,
                        )
                    )
                )
            except (ValueError, RuntimeError):
                pass
    independent_value = (
        min(independent_bracketed, key=lambda item: abs(item - value))
        if independent_bracketed
        else float(sigma_fit.x)
    )
    # In highly ill-conditioned high-mode coordinates a bounded sigma search
    # can find a numerically deeper false sub-valley.  Prefer the stationary
    # determinant root when the independently phased bracket confirms it.
    if abs(primary - independent_value) <= settings.independent_root_tolerance:
        value = primary
    left_vector, right_vector, singular, sigma_1, ratio, nullity = evaluator.svd(value)
    sigma_2 = float(singular[-2]) if singular.size >= 2 else float("inf")
    margin = max(settings.stationary_tolerance, 1.0e-5 * (right - left))
    stationary = bool(left + margin < primary < right - margin and narrow_left + margin < value < narrow_right - margin)
    stable = (
        abs(narrow_value - primary) <= settings.stationary_tolerance
        or abs(primary - independent_value) <= settings.independent_root_tolerance
    )
    independent = abs(value - independent_value) <= settings.independent_root_tolerance
    svd_pass = sigma_1 <= settings.sigma_accept and ratio <= settings.sigma_ratio_accept
    boundary_minimum = bool(
        primary <= left + margin
        or primary >= right - margin
        or value <= narrow_left + margin
        or value >= narrow_right - margin
    )
    if boundary_minimum:
        status = SEED_REJECTED_BOUNDARY_MINIMUM
    elif not stationary:
        status = SEED_REJECTED_NONSTATIONARY
    elif not stable or not independent:
        status = SEED_REJECTED_INDEPENDENT_DISAGREEMENT
    elif not svd_pass:
        status = SEED_REJECTED_INSUFFICIENT_SINGULARITY
    else:
        status = SEED_REFINED_TO_NEW_ROOT
    if status == SEED_REFINED_TO_NEW_ROOT:
        operations.seed_windows_accepted += 1
    else:
        operations.seed_windows_rejected += 1
    return RootRefinementResult(
        branch_id=branch.branch_id,
        seed=float(seed),
        window_left=left,
        window_right=right,
        projected_candidate=float(projected_fit.x),
        primary_candidate=primary,
        narrow_candidate=value,
        independent_candidate=independent_value,
        Lambda=value,
        sigma_1=sigma_1,
        sigma_2=sigma_2,
        sigma_ratio=ratio,
        stationary=stationary,
        stable_refinement=stable,
        independent_agreement=independent,
        full_matrix_svd_pass=svd_pass,
        acceptance_status=status,
        right_vector=tuple(float(item) for item in right_vector),
        left_vector=tuple(float(item) for item in left_vector),
        nullity=nullity,
    )


def _cluster_groups(branches: Sequence[ContinuedBranch], settings: ContinuationSettings) -> list[list[int]]:
    order = sorted(range(len(branches)), key=lambda index: branches[index].Lambda)
    groups: list[list[int]] = []
    for index in order:
        if not groups:
            groups.append([index])
            continue
        previous = branches[groups[-1][-1]]
        current = branches[index]
        gap = current.Lambda - previous.Lambda
        threshold = max(settings.cluster_gap_absolute, settings.cluster_gap_relative * max(current.Lambda, 1.0))
        if gap <= threshold:
            groups[-1].append(index)
        else:
            groups.append([index])
    return groups


def _mark_duplicate_seed_refinements(
    refinements: Sequence[RootRefinementResult],
    expected_dimension: int,
    settings: ContinuationSettings,
    operations: BranchOperationCounts,
) -> list[RootRefinementResult]:
    """Mark extra seeds in one stationary basin without creating extra roots."""
    merged = list(refinements)
    accepted_order = sorted(
        (
            index
            for index, item in enumerate(merged)
            if item.acceptance_status == SEED_REFINED_TO_NEW_ROOT
        ),
        key=lambda index: merged[index].Lambda,
    )
    for previous_index, current_index in zip(accepted_order, accepted_order[1:]):
        previous_refinement = merged[previous_index]
        current_refinement = merged[current_index]
        if (
            current_refinement.Lambda - previous_refinement.Lambda
            <= settings.root_match_tolerance
            and max(previous_refinement.nullity, current_refinement.nullity) < expected_dimension
        ):
            merged[current_index] = replace(
                current_refinement,
                acceptance_status=SEED_REFINED_TO_EXISTING_ROOT,
            )
            operations.seed_windows_accepted -= 1
            operations.seed_windows_merged += 1
    return merged


def _reduced_candidates(
    evaluator: _Evaluator,
    branches: Sequence[ContinuedBranch],
    left: float,
    right: float,
    expected: int,
    operations: BranchOperationCounts,
) -> tuple[float, ...]:
    u, _ = np.linalg.qr(np.column_stack([np.asarray(item.left_vector) for item in branches]))
    v, _ = np.linalg.qr(np.column_stack([np.asarray(item.right_vector) for item in branches]))

    def reduced_sigma(value: float) -> float:
        operations.local_reduced_matrix_SVD_calls += 1
        reduced = u.T @ complete.row_normalized(evaluator.matrix(value)) @ v
        return float(np.linalg.svd(reduced, compute_uv=False)[-1])

    grid_count = max(65, 32 * expected + 1)
    grid = np.linspace(left, right, grid_count)
    values = np.asarray([reduced_sigma(float(item)) for item in grid])
    candidates: list[float] = []
    for index in range(1, len(grid) - 1):
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]:
            fit = minimize_scalar(
                reduced_sigma,
                bounds=(float(grid[index - 1]), float(grid[index + 1])),
                method="bounded",
                options={"xatol": 1.0e-11},
            )
            operations.local_bounded_minimizations += 1
            if all(abs(float(fit.x) - old) > 2.0e-5 for old in candidates):
                candidates.append(float(fit.x))
    selected = sorted(candidates, key=lambda value: reduced_sigma(value))[: max(expected, 1)]
    return tuple(sorted(selected))


def _attempt_step(
    model: str,
    geometry: complete.Geometry,
    beta_from: float,
    beta_to: float,
    previous: Sequence[ContinuedBranch],
    earlier_by_id: Mapping[str, ContinuedBranch],
    settings: ContinuationSettings,
    operations: BranchOperationCounts,
) -> ContinuationStep:
    target_geometry = complete.Geometry(geometry.epsilon_0, beta_to, geometry.mu, geometry.eta)
    evaluator = _Evaluator(complete.model_matrix_provider(model, target_geometry), operations, "local")
    groups = _cluster_groups(previous, settings)
    all_refinements: list[RootRefinementResult] = []
    all_clusters: list[BranchCluster] = []
    candidates: list[tuple[ContinuedBranch, float, RootRefinementResult, str]] = []
    for group_number, indices in enumerate(groups, start=1):
        members = [previous[index] for index in indices]
        predictions = [_predict(item, earlier_by_id.get(item.branch_id), beta_to) for item in members]
        cluster_id = f"beta{beta_to:.8g}_C{group_number:03d}" if len(members) > 1 else ""
        reduced: tuple[float, ...] = ()
        if len(members) > 1:
            operations.cluster_steps_attempted += 1
            combined_left = max(settings.lambda_min, min(predictions) - 1.5 * settings.seed_half_width)
            combined_right = min(settings.lambda_max, max(predictions) + 1.5 * settings.seed_half_width)
            reduced = _reduced_candidates(
                evaluator, members, combined_left, combined_right, len(members), operations
            )
        member_refinements: list[RootRefinementResult] = []
        member_seeds = [
            reduced[index] if index < len(reduced) else prediction
            for index, prediction in enumerate(predictions)
        ]
        for member_index, (member, prediction, seed) in enumerate(
            zip(members, predictions, member_seeds)
        ):
            cluster_half_width: float | None = None
            if len(members) > 1:
                neighbor_gaps = [
                    abs(seed - value)
                    for index, value in enumerate(member_seeds)
                    if index != member_index and abs(seed - value) > 1.0e-10
                ]
                if neighbor_gaps:
                    cluster_half_width = 0.45 * min(neighbor_gaps)
            refined = _refine_seed(
                evaluator,
                member,
                seed,
                settings,
                operations,
                half_width=cluster_half_width,
            )
            member_refinements.append(refined)
            all_refinements.append(refined)
        if len(member_refinements) > 1:
            member_refinements = _mark_duplicate_seed_refinements(
                member_refinements,
                len(members),
                settings,
                operations,
            )
            all_refinements[-len(member_refinements) :] = member_refinements
        cluster_mac = 1.0
        resolved = all(item.acceptance_status == SEED_REFINED_TO_NEW_ROOT for item in member_refinements)
        if len(members) > 1 and resolved:
            refined_values = sorted(item.Lambda for item in member_refinements)
            distinct = all(
                refined_values[index + 1] - refined_values[index] > settings.root_match_tolerance
                for index in range(len(refined_values) - 1)
            )
            if not distinct:
                shared_nullity = max(item.nullity for item in member_refinements)
                resolved = shared_nullity >= len(members)
            cluster_mac = _subspace_mac(
                [item.right_vector for item in members],
                [item.right_vector for item in member_refinements],
            )
            resolved = math.isfinite(cluster_mac) and cluster_mac >= settings.subspace_mac_accept
        if len(members) > 1:
            if resolved:
                operations.cluster_steps_resolved += 1
            all_clusters.append(
                BranchCluster(
                    cluster_id=cluster_id,
                    beta_from_deg=beta_from,
                    beta_to_deg=beta_to,
                    branch_ids=tuple(item.branch_id for item in members),
                    expected_dimension=len(members),
                    found_dimension=sum(item.acceptance_status == SEED_REFINED_TO_NEW_ROOT for item in member_refinements),
                    subspace_MAC=cluster_mac,
                    reduced_candidates=reduced,
                    resolved=resolved,
                    status="cluster_resolved" if resolved else "cluster_unresolved",
                )
            )
        for member, prediction, refined in zip(members, predictions, member_refinements):
            candidates.append((member, prediction, refined, cluster_id))
    required_ids = {
        item.branch_id for item in sorted(previous, key=lambda value: value.Lambda)[:11]
    }
    full_prefix_accepted = all(
        item[2].acceptance_status == SEED_REFINED_TO_NEW_ROOT for item in candidates
    ) and all(item.resolved for item in all_clusters)
    required_accepted = all(
        item[2].acceptance_status == SEED_REFINED_TO_NEW_ROOT
        for item in candidates
        if item[0].branch_id in required_ids
    )
    required_clusters_resolved = all(
        item.resolved
        for item in all_clusters
        if required_ids.intersection(item.branch_ids)
    )
    accepted = required_accepted and required_clusters_resolved
    current: list[ContinuedBranch] = []
    if accepted:
        old = list(previous) if full_prefix_accepted else [item for item in previous if item.branch_id in required_ids]
        active_candidates = (
            list(candidates)
            if full_prefix_accepted
            else [item for item in candidates if item[0].branch_id in required_ids]
        )
        new_vectors = [np.asarray(item[2].right_vector) for item in active_candidates]
        mac = np.asarray(
            [[complete.coefficient_self_mac(item.right_vector, vector) for vector in new_vectors] for item in old]
        )
        cost = 1.0 - np.nan_to_num(mac, nan=0.0)
        old_predictions = [
            _predict(item, earlier_by_id.get(item.branch_id), beta_to) for item in old
        ]
        for row, (item, prediction) in enumerate(zip(old, old_predictions)):
            frequency_scale = max(settings.seed_half_width, 0.02 * max(item.Lambda, 1.0))
            for column, candidate in enumerate(active_candidates):
                cost[row, column] += 0.75 * min(
                    abs(candidate[2].Lambda - prediction) / frequency_scale,
                    20.0,
                )
        rows, columns = linear_sum_assignment(cost)
        assignment = {int(row): int(column) for row, column in zip(rows, columns)}
        for old_index, member in enumerate(old):
            candidate = active_candidates[assignment[old_index]]
            refined = candidate[2]
            individual_mac = complete.coefficient_self_mac(member.right_vector, refined.right_vector)
            if individual_mac < settings.mac_accept and not candidate[3]:
                accepted = False
            current.append(
                ContinuedBranch(
                    model=model,
                    branch_id=member.branch_id,
                    parent_family=member.parent_family,
                    beta_deg=beta_to,
                    predicted_Lambda=old_predictions[old_index],
                    Lambda=refined.Lambda,
                    right_vector=refined.right_vector,
                    left_vector=refined.left_vector,
                    sigma_1=refined.sigma_1,
                    sigma_ratio=refined.sigma_ratio,
                    self_MAC=individual_mac,
                    nullity=refined.nullity,
                    cluster_id=candidate[3],
                    refinement_status=refined.acceptance_status,
                    detection_source="cluster_reduced_window" if candidate[3] else "isolated_projected_window",
                )
            )
    return ContinuationStep(
        model=model,
        beta_from_deg=beta_from,
        beta_to_deg=beta_to,
        attempted_step_deg=beta_to - beta_from,
        accepted=accepted,
        status=(
            "accepted_local_continuation"
            if accepted and full_prefix_accepted
            else (
                "accepted_K10_prefix_optional_tail_unresolved"
                if accepted
                else "rejected_local_continuation"
            )
        ),
        branches=tuple(current if accepted else ()),
        clusters=tuple(all_clusters),
        refinements=tuple(all_refinements),
    )


def _branches_from_strict(
    model: str,
    beta: float,
    previous: Sequence[ContinuedBranch],
    result: complete.CompleteSpectrumResult,
) -> tuple[ContinuedBranch, ...]:
    roots = list(result.primary.roots[: len(previous)])
    if len(roots) < len(previous):
        return ()
    mac = np.asarray(
        [[complete.coefficient_self_mac(item.right_vector, root.null_vector) for root in roots] for item in previous]
    )
    scale = max(root.Lambda for root in roots)
    for row, item in enumerate(previous):
        for column, root in enumerate(roots):
            mac[row, column] -= 0.02 * abs(item.Lambda - root.Lambda) / max(scale, 1.0)
    rows, columns = linear_sum_assignment(1.0 - np.nan_to_num(mac, nan=0.0))
    assignment = {int(row): int(column) for row, column in zip(rows, columns)}
    out: list[ContinuedBranch] = []
    provider = complete.model_matrix_provider(model, result.geometry)
    for index, item in enumerate(previous):
        root = roots[assignment[index]]
        left, right, _singular, sigma, ratio, nullity = _canonical_svd(provider(root.Lambda))
        out.append(
            ContinuedBranch(
                model=model,
                branch_id=item.branch_id,
                parent_family=item.parent_family,
                beta_deg=beta,
                predicted_Lambda=item.Lambda,
                Lambda=root.Lambda,
                right_vector=tuple(float(value) for value in right),
                left_vector=tuple(float(value) for value in left),
                sigma_1=sigma,
                sigma_ratio=ratio,
                self_MAC=complete.coefficient_self_mac(item.right_vector, right),
                nullity=nullity,
                cluster_id=root.root_cluster_id,
                refinement_status=SEED_REFINED_TO_NEW_ROOT,
                detection_source="triggered_general_strict_fallback",
            )
        )
    return tuple(out)


def _guard(
    model: str,
    geometry: complete.Geometry,
    branches: Sequence[ContinuedBranch],
    settings: ContinuationSettings,
    operations: BranchOperationCounts,
) -> GlobalGuardResult:
    ordered = sorted(branches, key=lambda item: item.Lambda)
    if len(ordered) < 11:
        return GlobalGuardResult(0.0, (), (), (), ("fewer_than_11_candidates",), False, "guard_unresolved")
    guard_limit = (
        0.5 * (ordered[10].Lambda + ordered[11].Lambda)
        if len(ordered) >= 12
        else ordered[10].Lambda + max(3.0 * settings.root_match_tolerance, 0.02)
    )
    evaluator = _Evaluator(complete.model_matrix_provider(model, geometry), operations, "guard")
    step = settings.guard_scan_step
    grids = (
        np.arange(settings.lambda_min, guard_limit + 0.5 * step, step),
        np.arange(settings.lambda_min + 0.5 * step, guard_limit + 0.5 * step, step),
    )
    sign_candidates: list[float] = []
    sigma_candidates: list[float] = []
    for grid in grids:
        det = np.asarray([evaluator.determinant(float(value)) for value in grid])
        for index in range(len(grid) - 1):
            if not (math.isfinite(det[index]) and math.isfinite(det[index + 1])):
                continue
            if det[index] * det[index + 1] < 0.0:
                try:
                    candidate = float(brentq(evaluator.determinant, float(grid[index]), float(grid[index + 1])))
                except (ValueError, RuntimeError):
                    continue
                if evaluator.sigma(candidate) <= settings.sigma_accept:
                    sign_candidates.append(candidate)
        sigma = np.asarray([evaluator.sigma(float(value)) for value in grid])
        for index in range(1, len(grid) - 1):
            if sigma[index] <= sigma[index - 1] and sigma[index] <= sigma[index + 1] and sigma[index] <= 0.04:
                fit = minimize_scalar(
                    evaluator.sigma,
                    bounds=(float(grid[index - 1]), float(grid[index + 1])),
                    method="bounded",
                    options={"xatol": 1.0e-11},
                )
                _left, _right, _singular, candidate_sigma, candidate_ratio, _nullity = evaluator.svd(float(fit.x))
                if candidate_sigma <= settings.sigma_accept and candidate_ratio <= settings.sigma_ratio_accept:
                    sigma_candidates.append(float(fit.x))

    def unique(values: Sequence[float]) -> tuple[float, ...]:
        out: list[float] = []
        for value in sorted(values):
            if not out or abs(value - out[-1]) > settings.root_match_tolerance:
                out.append(value)
        return tuple(out)

    sign_unique = unique(sign_candidates)
    sigma_unique = unique(sigma_candidates)
    candidates = unique((*sign_unique, *sigma_unique))
    accepted = [item.Lambda for item in ordered if item.Lambda <= guard_limit + settings.root_match_tolerance]
    unmatched = tuple(
        value for value in candidates if all(abs(value - root) > 3.0 * settings.root_match_tolerance for root in accepted)
    )
    unresolved: list[str] = []
    if len(accepted) < 11:
        unresolved.append("fewer_than_11_roots_below_guard")
    if unmatched:
        unresolved.append("unmatched_guard_candidate")
    return GlobalGuardResult(
        guard_limit=guard_limit,
        sign_candidates=sign_unique,
        sigma_candidates=sigma_unique,
        unmatched_candidates=unmatched,
        unresolved_intervals=tuple(unresolved),
        passed=not unresolved,
        status="guard_resolved" if not unresolved else "guard_unresolved",
    )


def _compare_values(
    primary: Sequence[float], verification: Sequence[float], tolerance: float
) -> tuple[tuple[dict[str, object], ...], bool]:
    rows: list[dict[str, object]] = []
    agreement = len(primary) >= 11 and len(verification) >= 11
    for index in range(min(len(primary), len(verification), 12)):
        difference = abs(float(primary[index]) - float(verification[index]))
        passed = difference <= tolerance
        if index < 11:
            agreement = agreement and passed
        rows.append(
            {
                "sorted_index": index + 1,
                "primary_Lambda": float(primary[index]),
                "verification_Lambda": float(verification[index]),
                "absolute_difference": difference,
                "within_tolerance": passed,
            }
        )
    return tuple(rows), bool(agreement)


def resolve_branch_spectrum(
    model: str,
    geometry: complete.Geometry,
    *,
    settings: ContinuationSettings | None = None,
    force_strict_verification: bool = False,
) -> BranchContinuationResult:
    active = settings or ContinuationSettings()
    active.validate()
    geometry.validate()
    if model not in complete.SUPPORTED_MODELS:
        raise ValueError(f"unsupported model: {model}")
    operations = BranchOperationCounts()
    parent = _block_spectrum(model, geometry, active, operations)
    branches = _parent_as_branches(parent)
    steps: list[ContinuationStep] = []
    fallback_rows: list[dict[str, object]] = []
    reached_target = abs(geometry.beta_deg) <= 1.0e-14
    if geometry.beta_deg > 1.0e-14:
        beta = 0.0
        step_size = min(active.beta_initial_step_deg, geometry.beta_deg)
        earlier: tuple[ContinuedBranch, ...] = ()
        while beta < geometry.beta_deg - 1.0e-12:
            target = min(geometry.beta_deg, beta + step_size)
            operations.beta_steps_attempted += 1
            earlier_map = {item.branch_id: item for item in earlier}
            attempted = _attempt_step(model, geometry, beta, target, branches, earlier_map, active, operations)
            steps.append(attempted)
            if attempted.accepted:
                earlier = branches
                branches = attempted.branches
                beta = target
                operations.beta_steps_accepted += 1
                step_size = min(active.beta_max_step_deg, step_size * active.beta_growth_factor)
                continue
            if step_size * active.beta_shrink_factor >= active.beta_min_step_deg:
                step_size *= active.beta_shrink_factor
                operations.beta_step_reductions += 1
                continue
            if not active.allow_strict_fallback:
                fallback_rows.append(
                    {
                        "beta_deg": target,
                        "trigger": attempted.status,
                        "strict_spectrum_status": "disabled",
                        "strict_exclusion_reason": "strict_fallback_disabled",
                        "recovered_branch_count": 0,
                        "direct_seed_acceptance_present": False,
                    }
                )
                break
            fallback_geometry = complete.Geometry(geometry.epsilon_0, target, geometry.mu, geometry.eta)
            strict = complete.resolve_general_spectrum(
                model,
                fallback_geometry,
                settings=active.strict_settings(),
                continuation_seeds=[item.Lambda for item in branches],
            )
            operations.strict_fallback_runs += 1
            operations.strict_fallback_matrix_evaluations += strict.operations.characteristic_matrix_evaluations
            operations.strict_fallback_full_6x6_SVD_calls += strict.operations.full_6x6_SVD_calls
            recovered = _branches_from_strict(model, target, branches, strict)
            fallback_rows.append(
                {
                    "beta_deg": target,
                    "trigger": attempted.status,
                    "strict_spectrum_status": strict.spectrum_status,
                    "strict_exclusion_reason": strict.exclusion_reason,
                    "recovered_branch_count": len(recovered),
                    "direct_seed_acceptance_present": any(
                        "direct_full_matrix_SVD" in source
                        for root in strict.primary.roots
                        for source in root.detection_sources
                    ),
                }
            )
            if len(recovered) != len(branches):
                break
            earlier = branches
            branches = recovered
            beta = target
            operations.beta_steps_accepted += 1
            step_size = active.beta_initial_step_deg
        reached_target = beta >= geometry.beta_deg - 1.0e-12
    ordered = tuple(sorted(branches, key=lambda item: item.Lambda))
    guard = (
        _guard(model, geometry, ordered, active, operations)
        if active.run_global_guard and reached_target
        else GlobalGuardResult(
            0.0,
            (),
            (),
            (),
            ("global_guard_disabled" if not active.run_global_guard else "continuation_target_not_reached",),
            False,
            "guard_disabled" if not active.run_global_guard else "guard_unresolved",
        )
    )
    if reached_target and active.run_global_guard and guard.unmatched_candidates:
        provider = complete.model_matrix_provider(model, geometry)
        recovered = list(ordered)
        for recovery_index, value in enumerate(guard.unmatched_candidates, start=1):
            if any(abs(value - item.Lambda) <= active.root_match_tolerance for item in recovered):
                continue
            left, right, _singular, sigma, ratio, nullity = _canonical_svd(provider(value))
            if sigma > active.sigma_accept or ratio > active.sigma_ratio_accept:
                continue
            recovered.append(
                ContinuedBranch(
                    model=model,
                    branch_id=f"GUARD_RECOVERED_{recovery_index:03d}",
                    parent_family="guard_recovery",
                    beta_deg=geometry.beta_deg,
                    predicted_Lambda=value,
                    Lambda=value,
                    right_vector=tuple(float(item) for item in right),
                    left_vector=tuple(float(item) for item in left),
                    sigma_1=sigma,
                    sigma_ratio=ratio,
                    self_MAC=1.0,
                    nullity=nullity,
                    cluster_id="",
                    refinement_status=SEED_REFINED_TO_NEW_ROOT,
                    detection_source="cheap_global_guard_full_matrix_recovery",
                )
            )
            operations.guard_recovered_roots += 1
        ordered = tuple(sorted(recovered, key=lambda item: item.Lambda))
        guard = _guard(model, geometry, ordered, active, operations)
    local_values = tuple(item.Lambda for item in ordered)
    force_agreement: bool | None = None
    comparison_rows: tuple[dict[str, object], ...] = ()
    if force_strict_verification:
        strict = complete.resolve_general_spectrum(
            model,
            geometry,
            settings=active.strict_settings(),
            continuation_seeds=(),
        )
        operations.force_verification_runs += 1
        operations.force_verification_matrix_evaluations += strict.operations.characteristic_matrix_evaluations
        operations.force_verification_full_6x6_SVD_calls += strict.operations.full_6x6_SVD_calls
        comparison_rows, force_agreement = _compare_values(
            local_values, [root.Lambda for root in strict.verification.roots], active.independent_root_tolerance
        )
    oracle_agreement: bool | None = None
    if abs(geometry.beta_deg) <= 1.0e-14 and abs(geometry.eta) <= 1.0e-14:
        oracle = complete.straight_oracle_values(model, geometry, active.requested_roots)
        _oracle_rows, oracle_agreement = _compare_values(local_values, oracle, active.independent_root_tolerance)
    first11 = ordered[:11]
    no_seed_only = all(item.refinement_status == SEED_REFINED_TO_NEW_ROOT for item in first11)
    final_branch_ids = {item.branch_id for item in ordered[:11]}
    clusters_resolved = all(
        cluster.resolved
        for step in steps if step.accepted
        for cluster in step.clusters
        if final_branch_ids.intersection(cluster.branch_ids)
    )
    k10 = (
        len(first11) == 11
        and all(item.sigma_1 <= active.sigma_accept and item.sigma_ratio <= active.sigma_ratio_accept for item in first11)
        and no_seed_only
        and clusters_resolved
        and guard.passed
        and oracle_agreement is not False
        and force_agreement is not False
        and reached_target
        and all(not bool(row.get("direct_seed_acceptance_present")) for row in fallback_rows)
    )
    first12 = ordered[:12]
    full12 = (
        k10
        and len(first12) == 12
        and all(item.sigma_1 <= active.sigma_accept and item.sigma_ratio <= active.sigma_ratio_accept for item in first12)
        and (not comparison_rows or all(bool(row["within_tolerance"]) for row in comparison_rows[:12]))
    )
    reasons: list[str] = []
    if len(first11) < 11:
        reasons.append("fewer_than_11_roots")
    if not guard.passed:
        reasons.append("global_guard_unresolved")
    if not clusters_resolved:
        reasons.append("cluster_continuation_unresolved")
    if oracle_agreement is False:
        reasons.append("straight_oracle_disagreement_first11")
    if force_agreement is False:
        reasons.append("force_strict_disagreement_first11")
    if not no_seed_only:
        reasons.append("seed_only_candidate_present")
    if not reached_target:
        reasons.append("continuation_target_not_reached")
    return BranchContinuationResult(
        algorithm_version=BRANCH_CONTINUATION_ALGORITHM_VERSION,
        model=model,
        geometry=geometry,
        settings=active,
        parent_branches=parent,
        branches=ordered,
        steps=tuple(steps),
        guard=guard,
        primary_vs_verification=comparison_rows,
        strict_fallback_audit=tuple(fallback_rows),
        k10_guard_resolved=k10,
        full12_resolved=full12,
        spectrum_status="K10_guard_resolved" if k10 else "unresolved",
        exclusion_reason=";".join(reasons),
        oracle_agreement=oracle_agreement,
        force_verification_agreement=force_agreement,
        cache_status="not_requested",
        operations=operations,
    )


def _result_from_dict(data: Mapping[str, object], cache_status: str) -> BranchContinuationResult:
    def parent(item: Mapping[str, object]) -> ParentBranch:
        return ParentBranch(**dict(item))  # type: ignore[arg-type]

    def branch(item: Mapping[str, object]) -> ContinuedBranch:
        return ContinuedBranch(**dict(item))  # type: ignore[arg-type]

    def refinement(item: Mapping[str, object]) -> RootRefinementResult:
        payload = dict(item)
        payload["acceptance_status"] = _LEGACY_SEED_STATUS_ALIASES.get(
            str(payload.get("acceptance_status", "")),
            str(payload.get("acceptance_status", "")),
        )
        return RootRefinementResult(**payload)  # type: ignore[arg-type]

    def cluster(item: Mapping[str, object]) -> BranchCluster:
        return BranchCluster(**dict(item))  # type: ignore[arg-type]

    def step(item: Mapping[str, object]) -> ContinuationStep:
        payload = dict(item)
        payload["branches"] = tuple(branch(value) for value in payload.get("branches", ()))
        payload["clusters"] = tuple(cluster(value) for value in payload.get("clusters", ()))
        payload["refinements"] = tuple(refinement(value) for value in payload.get("refinements", ()))
        return ContinuationStep(**payload)  # type: ignore[arg-type]

    geometry = complete.Geometry(**dict(data["geometry"]))  # type: ignore[arg-type]
    settings = ContinuationSettings(**dict(data["settings"]))  # type: ignore[arg-type]
    guard_payload = dict(data["guard"])  # type: ignore[arg-type]
    guard = GlobalGuardResult(**guard_payload)  # type: ignore[arg-type]
    operations = BranchOperationCounts(**dict(data["operations"]))  # type: ignore[arg-type]
    return BranchContinuationResult(
        algorithm_version=str(data["algorithm_version"]),
        model=str(data["model"]),
        geometry=geometry,
        settings=settings,
        parent_branches=tuple(parent(item) for item in data.get("parent_branches", ())),  # type: ignore[arg-type]
        branches=tuple(branch(item) for item in data.get("branches", ())),  # type: ignore[arg-type]
        steps=tuple(step(item) for item in data.get("steps", ())),  # type: ignore[arg-type]
        guard=guard,
        primary_vs_verification=tuple(dict(item) for item in data.get("primary_vs_verification", ())),  # type: ignore[arg-type]
        strict_fallback_audit=tuple(dict(item) for item in data.get("strict_fallback_audit", ())),  # type: ignore[arg-type]
        k10_guard_resolved=bool(data["k10_guard_resolved"]),
        full12_resolved=bool(data["full12_resolved"]),
        spectrum_status=str(data["spectrum_status"]),
        exclusion_reason=str(data.get("exclusion_reason", "")),
        oracle_agreement=data.get("oracle_agreement"),  # type: ignore[arg-type]
        force_verification_agreement=data.get("force_verification_agreement"),  # type: ignore[arg-type]
        cache_status=cache_status,
        operations=operations,
    )


class BranchContinuationCache:
    def __init__(
        self,
        cache_dir: Path,
        *,
        reuse_cache: bool = True,
        force_recompute: bool = False,
        verification_scope: str = "primary",
    ) -> None:
        if verification_scope not in {"primary", "force_strict_verification"}:
            raise ValueError("verification_scope must be primary or force_strict_verification")
        self.cache_dir = Path(cache_dir)
        self.reuse_cache = bool(reuse_cache)
        self.force_recompute = bool(force_recompute)
        self.verification_scope = verification_scope
        self.last_load_status = "not_checked"

    @staticmethod
    def identity(
        model: str,
        geometry: complete.Geometry,
        settings: ContinuationSettings,
        verification_scope: str,
    ) -> dict[str, object]:
        return {
            "algorithm_version": BRANCH_CONTINUATION_ALGORITHM_VERSION,
            "model": model,
            "geometry": {
                "epsilon_0": float(geometry.epsilon_0),
                "beta_deg": float(geometry.beta_deg),
                "mu": float(geometry.mu),
                "eta": float(geometry.eta),
            },
            "beta_path": {
                "start_beta_deg": 0.0,
                "target_beta_deg": float(geometry.beta_deg),
                "policy": "adaptive_local_continuation",
            },
            "beta_targets": [0.0, float(geometry.beta_deg)],
            "settings": asdict(settings),
            "strict_fallback_settings": asdict(settings.strict_settings()),
            "verification_scope": verification_scope,
        }

    def path(self, model: str, geometry: complete.Geometry, settings: ContinuationSettings) -> Path:
        identity = self.identity(model, geometry, settings, self.verification_scope)
        encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"))
        digest = hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:24]
        token = "eb" if model == complete.MODEL_EB else "timo"
        return self.cache_dir / self.verification_scope / f"branch_spectrum_{token}_{digest}.json"

    def load(
        self, model: str, geometry: complete.Geometry, settings: ContinuationSettings
    ) -> BranchContinuationResult | None:
        if self.force_recompute or not self.reuse_cache:
            self.last_load_status = "force_recompute" if self.force_recompute else "cache_disabled"
            return None
        path = self.path(model, geometry, settings)
        if not path.exists():
            self.last_load_status = "miss"
            return None
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            self.last_load_status = "invalid_cache_payload"
            return None
        identity = self.identity(model, geometry, settings, self.verification_scope)
        if payload.get("identity") != identity:
            self.last_load_status = "stale_cache_identity"
            return None
        self.last_load_status = "hit"
        return _result_from_dict(payload["result"], "hit")

    def save(self, result: BranchContinuationResult) -> Path:
        path = self.path(result.model, result.geometry, result.settings)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "identity": self.identity(result.model, result.geometry, result.settings, self.verification_scope),
            "result": asdict(result),
        }
        path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")), encoding="utf-8")
        return path

    def resolve(
        self,
        model: str,
        geometry: complete.Geometry,
        settings: ContinuationSettings,
    ) -> BranchContinuationResult:
        cached = self.load(model, geometry, settings)
        if cached is not None:
            cached.operations.cache_hits += 1
            return cached
        load_status = self.last_load_status
        result = resolve_branch_spectrum(
            model,
            geometry,
            settings=settings,
            force_strict_verification=self.verification_scope == "force_strict_verification",
        )
        result.operations.cache_misses += 1
        result = BranchContinuationResult(**{**result.__dict__, "cache_status": load_status})
        self.save(result)
        return result


__all__ = [
    "BETA0_BLOCK_MAPS",
    "BRANCH_CONTINUATION_ALGORITHM_VERSION",
    "BlockMap",
    "BranchCluster",
    "BranchContinuationCache",
    "BranchContinuationResult",
    "BranchOperationCounts",
    "ContinuationSettings",
    "ContinuationStep",
    "ContinuedBranch",
    "GlobalGuardResult",
    "ParentBranch",
    "RootRefinementResult",
    "SEED_REFINED_TO_EXISTING_ROOT",
    "SEED_REFINED_TO_NEW_ROOT",
    "SEED_REJECTED_BOUNDARY_MINIMUM",
    "SEED_REJECTED_INDEPENDENT_DISAGREEMENT",
    "SEED_REJECTED_INSUFFICIENT_SINGULARITY",
    "SEED_REJECTED_NONSTATIONARY",
    "resolve_branch_spectrum",
]
