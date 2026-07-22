from __future__ import annotations

from dataclasses import asdict, dataclass, field, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from scripts.analysis.compare_single_rod_eb_timoshenko import fixed_fixed_eb_roots
from scripts.lib import straight_rod_factorized_spectrum as straight_factorized
from scripts.lib import variable_length_timoshenko as TIMO
from src.my_project.analytic.formulas_thickness_mismatch import (
    assemble_clamped_coupled_matrix_eta,
)


GENERAL_SPECTRUM_ALGORITHM_VERSION = "general_complete_svd_v1"
AUTO_SPECTRUM_ALGORITHM_VERSION = "auto_complete_spectrum_v1"

MODEL_EB = "Euler-Bernoulli"
MODEL_TIMO = "Timoshenko"
SUPPORTED_MODELS = (MODEL_EB, MODEL_TIMO)

DEFAULT_SCAN_START = 0.2
DEFAULT_SCAN_STEP = 0.01
DEFAULT_LOCAL_POINTS = 17
DEFAULT_ADAPTIVE_DEPTH = 2
DEFAULT_SIGMA_PREFILTER = 5.0e-2
DEFAULT_SIGMA_ACCEPT = 5.0e-6
DEFAULT_SIGMA_RATIO_ACCEPT = 5.0e-4
DEFAULT_NULLITY_SIGMA = 1.0e-10
DEFAULT_ROOT_DEDUP_TOL = 2.0e-4
DEFAULT_ROOT_MATCH_TOL = 2.0e-4
DEFAULT_CLOSE_GAP = 1.0e-2
DEFAULT_SEED_HALF_WIDTH = 3.0e-2
DEFAULT_BRENT_XTOL = 1.0e-12
DEFAULT_BRENT_RTOL = 1.0e-12


MatrixProvider = Callable[[float], np.ndarray]


@dataclass(frozen=True)
class Geometry:
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float

    def validate(self) -> None:
        values = (self.epsilon_0, self.beta_deg, self.mu, self.eta)
        if not all(math.isfinite(float(value)) for value in values):
            raise ValueError("geometry parameters must be finite")
        if self.epsilon_0 <= 0.0:
            raise ValueError("epsilon_0 must be positive")
        if not 0.0 <= self.beta_deg <= 90.0:
            raise ValueError("beta_deg must lie in [0, 90]")
        if not -1.0 < self.mu < 1.0:
            raise ValueError("mu must lie inside (-1, 1)")
        if not -1.0 < self.eta < 1.0:
            raise ValueError("eta must lie inside (-1, 1)")


@dataclass(frozen=True)
class SearchSettings:
    requested_roots: int = 12
    candidate_roots: int = 20
    verification_candidate_roots: int = 24
    lambda_min: float = DEFAULT_SCAN_START
    lambda_max: float | None = None
    scan_step: float = DEFAULT_SCAN_STEP
    shifted_grid_phase: float = 0.5
    local_points: int = DEFAULT_LOCAL_POINTS
    adaptive_depth: int = DEFAULT_ADAPTIVE_DEPTH
    sigma_prefilter: float = DEFAULT_SIGMA_PREFILTER
    sigma_accept: float = DEFAULT_SIGMA_ACCEPT
    sigma_ratio_accept: float = DEFAULT_SIGMA_RATIO_ACCEPT
    nullity_sigma: float = DEFAULT_NULLITY_SIGMA
    root_dedup_tol: float = DEFAULT_ROOT_DEDUP_TOL
    root_match_tol: float = DEFAULT_ROOT_MATCH_TOL
    close_gap_threshold: float = DEFAULT_CLOSE_GAP
    seed_half_width: float = DEFAULT_SEED_HALF_WIDTH
    brent_xtol: float = DEFAULT_BRENT_XTOL
    brent_rtol: float = DEFAULT_BRENT_RTOL
    max_upper_growth_tries: int = 4
    upper_growth_factor: float = 1.35

    def validate(self, *, k_max: int = 10) -> None:
        if self.requested_roots < int(k_max) + 2:
            raise ValueError("requested_roots must be at least k_max + 2")
        if self.candidate_roots < self.requested_roots:
            raise ValueError("candidate_roots must be at least requested_roots")
        if self.verification_candidate_roots <= self.candidate_roots:
            raise ValueError("verification_candidate_roots must exceed candidate_roots")
        positive = {
            "lambda_min": self.lambda_min,
            "scan_step": self.scan_step,
            "sigma_prefilter": self.sigma_prefilter,
            "sigma_accept": self.sigma_accept,
            "sigma_ratio_accept": self.sigma_ratio_accept,
            "nullity_sigma": self.nullity_sigma,
            "root_dedup_tol": self.root_dedup_tol,
            "root_match_tol": self.root_match_tol,
            "close_gap_threshold": self.close_gap_threshold,
            "seed_half_width": self.seed_half_width,
            "brent_xtol": self.brent_xtol,
            "brent_rtol": self.brent_rtol,
            "upper_growth_factor": self.upper_growth_factor,
        }
        if any(not math.isfinite(float(value)) or float(value) <= 0.0 for value in positive.values()):
            raise ValueError("scan/refinement/SVD settings must be finite and positive")
        if self.lambda_max is not None and self.lambda_max <= self.lambda_min:
            raise ValueError("lambda_max must exceed lambda_min")
        if not 0.0 <= self.shifted_grid_phase < 1.0:
            raise ValueError("shifted_grid_phase must lie in [0, 1)")
        if self.local_points < 9 or self.local_points % 2 == 0:
            raise ValueError("local_points must be an odd integer >= 9")
        if self.adaptive_depth < 0 or self.max_upper_growth_tries < 1:
            raise ValueError("adaptive depth and upper-growth tries are invalid")
        if self.upper_growth_factor <= 1.0:
            raise ValueError("upper_growth_factor must exceed one")


@dataclass
class OperationCounts:
    characteristic_matrix_evaluations: int = 0
    normalized_determinant_evaluations: int = 0
    full_6x6_SVD_calls: int = 0
    shifted_grid_evaluations: int = 0
    half_step_grid_evaluations: int = 0
    adaptive_interval_subdivisions: int = 0
    bounded_sigma_minimizations: int = 0
    Brent_or_bisection_evaluations: int = 0
    projected_residual_evaluations: int = 0
    continuation_seed_windows: int = 0
    EB_seed_windows_for_Timo: int = 0
    local_recovery_attempts: int = 0
    accepted_recovered_roots: int = 0
    rejected_false_valleys: int = 0
    multiplicity_perturbation_runs: int = 0
    independent_verification_runs: int = 0
    cache_hits: int = 0
    cache_misses: int = 0

    def add(self, other: "OperationCounts") -> None:
        for name in self.__dataclass_fields__:
            setattr(self, name, int(getattr(self, name)) + int(getattr(other, name)))


@dataclass(frozen=True)
class MatrixDiagnostics:
    Lambda: float
    determinant: float
    sigma_1: float
    sigma_2: float
    sigma_3: float
    sigma_ratio: float
    null_vector: tuple[float, ...]
    null_vector_pivot: int
    finite_matrix_status: str


@dataclass(frozen=True)
class RootCandidate:
    Lambda: float
    detection_sources: tuple[str, ...]
    diagnostics: MatrixDiagnostics
    interval_left: float
    interval_right: float
    interior_minimum: bool
    acceptance_status: str
    notes: str = ""


@dataclass(frozen=True)
class RootRecord:
    Lambda: float
    sorted_index: int
    root_cluster_id: str
    cluster_member_index: int
    cluster_size: int
    detection_sources: tuple[str, ...]
    determinant: float
    sigma_1: float
    sigma_2: float
    sigma_ratio: float
    null_vector: tuple[float, ...]
    null_vector_pivot: int
    self_MAC: float
    detected_nullity: int
    track_multiplicity: int
    multiplicity_status: str
    acceptance_status: str
    notes: str = ""


@dataclass(frozen=True)
class SearchConfigurationResult:
    configuration: str
    scan_step: float
    grid_phase: float
    candidate_root_target: int
    lambda_upper: float
    roots: tuple[RootRecord, ...]
    candidates: tuple[RootCandidate, ...]
    interval_rows: tuple[dict[str, object], ...]
    unresolved_intervals: tuple[str, ...]
    operations: OperationCounts


@dataclass(frozen=True)
class CompleteSpectrumResult:
    algorithm_version: str
    model: str
    geometry: Geometry
    settings: SearchSettings
    primary: SearchConfigurationResult
    verification: SearchConfigurationResult
    roots: tuple[RootRecord, ...]
    primary_vs_verification: tuple[dict[str, object], ...]
    independent_agreement: bool
    root11_available: bool
    root12_available: bool
    root12_boundary_warning: bool
    spectrum_status: str
    exclusion_reason: str
    cache_status: str
    operations: OperationCounts

    @property
    def values(self) -> tuple[float, ...]:
        return tuple(root.Lambda for root in self.roots)


def row_normalized(matrix: np.ndarray) -> np.ndarray:
    out = np.asarray(matrix, dtype=float).copy()
    if out.ndim != 2 or out.shape[0] != out.shape[1]:
        raise ValueError("characteristic matrix must be square")
    norms = np.linalg.norm(out, axis=1)
    finite = norms[np.isfinite(norms)]
    scale = float(np.max(finite)) if finite.size else 1.0
    # A structurally zero row can contain round-off at an exact axial root.
    # Dividing that O(eps) row by its own norm would manufacture an O(1) row
    # and destroy the nullspace.  Rows above the machine-level threshold keep
    # ordinary row normalization; only numerical zero rows are left unscaled.
    zero_row_threshold = 64.0 * np.finfo(float).eps * max(scale, 1.0)
    norms[~np.isfinite(norms) | (norms <= zero_row_threshold)] = 1.0
    return out / norms[:, None]


def coefficient_self_mac(left: Sequence[float], right: Sequence[float]) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    if a.shape != b.shape or a.size == 0 or not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return float("nan")
    denom = float(np.dot(a, a) * np.dot(b, b))
    if denom <= 0.0:
        return float("nan")
    return float(np.clip(abs(float(np.dot(a, b))) ** 2 / denom, 0.0, 1.0))


def model_matrix_provider(model: str, geometry: Geometry) -> MatrixProvider:
    geometry.validate()
    if model == MODEL_EB:
        beta_rad = float(np.deg2rad(geometry.beta_deg))

        def eb_matrix(value: float) -> np.ndarray:
            return assemble_clamped_coupled_matrix_eta(
                float(value), beta_rad, geometry.mu, geometry.epsilon_0, geometry.eta
            )

        return eb_matrix
    if model == MODEL_TIMO:

        def timo_matrix(value: float) -> np.ndarray:
            matrix, _warnings = TIMO.timo_coupling_matrix(
                float(value), geometry.beta_deg, geometry.mu, geometry.epsilon_0, geometry.eta
            )
            return matrix

        return timo_matrix
    raise ValueError(f"unsupported model: {model}")


class _MatrixEvaluator:
    def __init__(self, provider: MatrixProvider, operations: OperationCounts) -> None:
        self.provider = provider
        self.operations = operations
        self._diagnostic_cache: dict[float, MatrixDiagnostics] = {}
        self._determinant_cache: dict[float, float] = {}

    def _matrix(self, value: float) -> np.ndarray | None:
        self.operations.characteristic_matrix_evaluations += 1
        try:
            matrix = np.asarray(self.provider(float(value)), dtype=float)
        except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
            return None
        if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1] or not np.all(np.isfinite(matrix)):
            return None
        return matrix

    def determinant(self, value: float) -> float:
        key = round(float(value), 14)
        if key in self._diagnostic_cache:
            return self._diagnostic_cache[key].determinant
        if key in self._determinant_cache:
            return self._determinant_cache[key]
        matrix = self._matrix(value)
        self.operations.normalized_determinant_evaluations += 1
        if matrix is None:
            return float("nan")
        try:
            result = float(np.linalg.det(row_normalized(matrix)))
        except (ValueError, np.linalg.LinAlgError):
            result = float("nan")
        self._determinant_cache[key] = result
        return result

    def diagnostics(self, value: float) -> MatrixDiagnostics:
        key = round(float(value), 14)
        cached = self._diagnostic_cache.get(key)
        if cached is not None:
            return cached
        matrix = self._matrix(value)
        self.operations.normalized_determinant_evaluations += 1
        self.operations.full_6x6_SVD_calls += 1
        if matrix is None:
            result = MatrixDiagnostics(
                float(value), float("nan"), float("inf"), float("inf"), float("inf"), float("inf"), (), -1, "nonfinite_matrix"
            )
            self._diagnostic_cache[key] = result
            return result
        try:
            scaled = row_normalized(matrix)
            determinant = float(np.linalg.det(scaled))
            _u, singular_values, vh = np.linalg.svd(scaled, full_matrices=True)
        except (ValueError, np.linalg.LinAlgError):
            result = MatrixDiagnostics(
                float(value), float("nan"), float("inf"), float("inf"), float("inf"), float("inf"), (), -1, "svd_failure"
            )
            self._diagnostic_cache[key] = result
            return result
        sigma_1 = float(singular_values[-1]) if singular_values.size else float("inf")
        sigma_2 = float(singular_values[-2]) if singular_values.size >= 2 else float("inf")
        sigma_3 = float(singular_values[-3]) if singular_values.size >= 3 else float("inf")
        ratio = sigma_1 / sigma_2 if math.isfinite(sigma_2) and sigma_2 > 0.0 else float("inf")
        vector = np.asarray(vh[-1], dtype=float) if vh.size else np.array([], dtype=float)
        norm = float(np.linalg.norm(vector))
        if norm > 0.0 and math.isfinite(norm):
            vector = vector / norm
        pivot = int(np.argmax(np.abs(vector))) if vector.size else -1
        if pivot >= 0 and vector[pivot] < 0.0:
            vector = -vector
        result = MatrixDiagnostics(
            Lambda=float(value),
            determinant=determinant,
            sigma_1=sigma_1,
            sigma_2=sigma_2,
            sigma_3=sigma_3,
            sigma_ratio=float(ratio),
            null_vector=tuple(float(item) for item in vector),
            null_vector_pivot=pivot,
            finite_matrix_status="finite",
        )
        self._diagnostic_cache[key] = result
        self._determinant_cache[key] = determinant
        return result


def _regular_grid(start: float, upper: float, step: float, phase: float) -> np.ndarray:
    first = float(start) + float(phase) * float(step)
    if first > upper:
        return np.array([], dtype=float)
    count = int(math.floor((float(upper) - first) / float(step))) + 1
    return first + float(step) * np.arange(count, dtype=float)


def _initial_upper(settings: SearchSettings, target: int) -> float:
    if settings.lambda_max is not None:
        return float(settings.lambda_max)
    return max(22.0, 1.05 * float(target + 2))


def _candidate_is_accepted(diag: MatrixDiagnostics, settings: SearchSettings) -> bool:
    if diag.finite_matrix_status != "finite" or not math.isfinite(diag.sigma_1):
        return False
    if diag.sigma_1 > settings.sigma_accept:
        return False
    nullity_two = (
        diag.sigma_2 <= settings.nullity_sigma
        and math.isfinite(diag.sigma_3)
        and diag.sigma_3 > 0.0
        and diag.sigma_2 / diag.sigma_3 <= settings.sigma_ratio_accept
    )
    return nullity_two or diag.sigma_ratio <= settings.sigma_ratio_accept


def _brent_candidate(
    evaluator: _MatrixEvaluator,
    left: float,
    right: float,
    source: str,
    settings: SearchSettings,
) -> RootCandidate | None:
    evaluations = 0

    def objective(value: float) -> float:
        nonlocal evaluations
        evaluations += 1
        return evaluator.determinant(float(value))

    try:
        root = float(
            brentq(
                objective,
                float(left),
                float(right),
                xtol=settings.brent_xtol,
                rtol=settings.brent_rtol,
                maxiter=180,
            )
        )
    except (RuntimeError, ValueError):
        return None
    finally:
        evaluator.operations.Brent_or_bisection_evaluations += int(evaluations)
    diag = evaluator.diagnostics(root)
    status = "accepted_full_matrix_svd" if _candidate_is_accepted(diag, settings) else "rejected_full_matrix_svd"
    return RootCandidate(root, (source,), diag, float(left), float(right), True, status)


def _sigma_candidate(
    evaluator: _MatrixEvaluator,
    left: float,
    right: float,
    source: str,
    settings: SearchSettings,
) -> RootCandidate | None:
    if right <= left:
        return None
    evaluator.operations.bounded_sigma_minimizations += 1
    result = minimize_scalar(
        lambda value: evaluator.diagnostics(float(value)).sigma_1,
        bounds=(float(left), float(right)),
        method="bounded",
        options={"xatol": max(settings.brent_xtol, 1.0e-13), "maxiter": 180},
    )
    if not result.success or not math.isfinite(float(result.x)):
        return None
    root = float(result.x)
    width = float(right) - float(left)
    edge_tol = max(1.0e-10, 1.0e-4 * width)
    interior = root > float(left) + edge_tol and root < float(right) - edge_tol
    # Repeat on a narrower window as an independent local refinement.
    half = max(width / 8.0, 8.0 * settings.brent_xtol)
    narrow_left = max(float(left), root - half)
    narrow_right = min(float(right), root + half)
    evaluator.operations.bounded_sigma_minimizations += 1
    repeated = minimize_scalar(
        lambda value: evaluator.diagnostics(float(value)).sigma_1,
        bounds=(narrow_left, narrow_right),
        method="bounded",
        options={"xatol": max(settings.brent_xtol, 1.0e-13), "maxiter": 180},
    )
    if repeated.success and math.isfinite(float(repeated.fun)) and float(repeated.fun) <= float(result.fun):
        root = float(repeated.x)
    diag = evaluator.diagnostics(root)
    accepted = interior and _candidate_is_accepted(diag, settings)
    if not accepted:
        evaluator.operations.rejected_false_valleys += 1
    status = "accepted_full_matrix_svd" if accepted else ("rejected_boundary_minimum" if not interior else "rejected_false_valley")
    return RootCandidate(root, (source,), diag, float(left), float(right), interior, status)


def _adaptive_interval_candidates(
    evaluator: _MatrixEvaluator,
    left: float,
    right: float,
    source: str,
    settings: SearchSettings,
) -> tuple[list[RootCandidate], dict[str, object]]:
    evaluator.operations.local_recovery_attempts += 1
    subdivisions = 2 ** int(settings.adaptive_depth)
    point_count = subdivisions * (int(settings.local_points) - 1) + 1
    evaluator.operations.adaptive_interval_subdivisions += max(0, subdivisions - 1)
    grid = np.linspace(float(left), float(right), point_count, dtype=float)
    diagnostics = [evaluator.diagnostics(float(value)) for value in grid]
    determinant = np.asarray([item.determinant for item in diagnostics], dtype=float)
    sigma = np.asarray([item.sigma_1 for item in diagnostics], dtype=float)
    ratio = np.asarray([item.sigma_ratio for item in diagnostics], dtype=float)
    candidates: list[RootCandidate] = []
    for index in range(len(grid) - 1):
        f_left = float(determinant[index])
        f_right = float(determinant[index + 1])
        if not (math.isfinite(f_left) and math.isfinite(f_right)):
            continue
        if f_left == 0.0:
            diag = diagnostics[index]
            status = "accepted_full_matrix_svd" if _candidate_is_accepted(diag, settings) else "rejected_full_matrix_svd"
            candidates.append(
                RootCandidate(float(grid[index]), (source + ":grid_zero",), diag, float(grid[index]), float(grid[index + 1]), True, status)
            )
        elif f_left * f_right < 0.0:
            candidate = _brent_candidate(
                evaluator, float(grid[index]), float(grid[index + 1]), source + ":sign_change", settings
            )
            if candidate is not None:
                candidates.append(candidate)
    for index in range(1, len(grid) - 1):
        current = float(sigma[index])
        current_ratio = float(ratio[index])
        if not (math.isfinite(current) and math.isfinite(current_ratio)):
            continue
        if current > settings.sigma_prefilter or current_ratio > 10.0 * settings.sigma_ratio_accept:
            continue
        is_sigma_min = current <= float(sigma[index - 1]) and current <= float(sigma[index + 1])
        is_ratio_min = current_ratio <= float(ratio[index - 1]) and current_ratio <= float(ratio[index + 1])
        if not (is_sigma_min or is_ratio_min):
            continue
        candidate = _sigma_candidate(
            evaluator, float(grid[index - 1]), float(grid[index + 1]), source + ":sigma_valley", settings
        )
        if candidate is not None:
            candidates.append(candidate)
    finite_sigma = sigma[np.isfinite(sigma)]
    row = {
        "Lambda_left": float(left),
        "Lambda_right": float(right),
        "determinant_sign_left": int(np.sign(determinant[0])) if math.isfinite(float(determinant[0])) else 0,
        "determinant_sign_right": int(np.sign(determinant[-1])) if math.isfinite(float(determinant[-1])) else 0,
        "sigma_left": float(sigma[0]),
        "sigma_midpoint": float(sigma[len(sigma) // 2]),
        "sigma_right": float(sigma[-1]),
        "minimum_sigma": float(np.min(finite_sigma)) if finite_sigma.size else float("inf"),
        "minimum_sigma_ratio": float(np.min(ratio[np.isfinite(ratio)])) if np.any(np.isfinite(ratio)) else float("inf"),
        "subdivision_depth": int(settings.adaptive_depth),
        "trigger_reasons": source,
        "local_minima_count": sum(
            math.isfinite(float(sigma[index]))
            and float(sigma[index]) <= float(sigma[index - 1])
            and float(sigma[index]) <= float(sigma[index + 1])
            for index in range(1, len(sigma) - 1)
        ),
        "resolution_status": "pending",
    }
    return candidates, row


def _merge_candidates(candidates: Sequence[RootCandidate], settings: SearchSettings) -> list[RootCandidate]:
    accepted = [item for item in candidates if item.acceptance_status == "accepted_full_matrix_svd"]
    accepted.sort(key=lambda item: (item.Lambda, item.diagnostics.sigma_1, item.detection_sources))
    merged: list[RootCandidate] = []
    for candidate in accepted:
        if not merged:
            merged.append(candidate)
            continue
        previous = merged[-1]
        mac = coefficient_self_mac(previous.diagnostics.null_vector, candidate.diagnostics.null_vector)
        previous_configurations = {source.split("_", 1)[0] for source in previous.detection_sources}
        candidate_configurations = {source.split("_", 1)[0] for source in candidate.detection_sources}
        intervals_overlap = max(previous.interval_left, candidate.interval_left) <= min(
            previous.interval_right, candidate.interval_right
        )
        compatible_history = intervals_overlap or bool(previous_configurations & candidate_configurations) or (
            any("sign_change" in source for source in previous.detection_sources)
            and any("sign_change" in source for source in candidate.detection_sources)
        )
        if (
            abs(candidate.Lambda - previous.Lambda) <= settings.root_dedup_tol
            and math.isfinite(mac)
            and mac >= 0.98
            and compatible_history
        ):
            best = candidate if candidate.diagnostics.sigma_1 < previous.diagnostics.sigma_1 else previous
            sources = tuple(sorted(set(previous.detection_sources + candidate.detection_sources)))
            merged[-1] = replace(best, detection_sources=sources)
        else:
            merged.append(candidate)
    return merged


def _root_records(candidates: Sequence[RootCandidate], settings: SearchSettings) -> tuple[RootRecord, ...]:
    raw: list[RootRecord] = []
    for candidate in candidates:
        nullity = 2 if (
            candidate.diagnostics.sigma_2 <= settings.nullity_sigma
            and math.isfinite(candidate.diagnostics.sigma_3)
            and candidate.diagnostics.sigma_3 > 0.0
            and candidate.diagnostics.sigma_2 / candidate.diagnostics.sigma_3 <= settings.sigma_ratio_accept
        ) else 1
        multiplicity_status = "verified_nullity_2" if nullity == 2 else "simple_root"
        copies = nullity
        for _copy in range(copies):
            raw.append(
                RootRecord(
                    Lambda=float(candidate.Lambda),
                    sorted_index=0,
                    root_cluster_id="",
                    cluster_member_index=1,
                    cluster_size=copies,
                    detection_sources=candidate.detection_sources,
                    determinant=candidate.diagnostics.determinant,
                    sigma_1=candidate.diagnostics.sigma_1,
                    sigma_2=candidate.diagnostics.sigma_2,
                    sigma_ratio=candidate.diagnostics.sigma_ratio,
                    null_vector=candidate.diagnostics.null_vector,
                    null_vector_pivot=candidate.diagnostics.null_vector_pivot,
                    self_MAC=1.0,
                    detected_nullity=nullity,
                    track_multiplicity=copies,
                    multiplicity_status=multiplicity_status,
                    acceptance_status=candidate.acceptance_status,
                    notes=candidate.notes,
                )
            )
    raw.sort(key=lambda item: (item.Lambda, item.null_vector_pivot, item.detection_sources))
    clusters: list[list[int]] = []
    for index, root in enumerate(raw):
        if not clusters or root.Lambda - raw[clusters[-1][-1]].Lambda > settings.close_gap_threshold:
            clusters.append([index])
        else:
            clusters[-1].append(index)
    result = list(raw)
    cluster_counter = 0
    for indices in clusters:
        if len(indices) <= 1:
            continue
        cluster_counter += 1
        cluster_id = f"root_cluster_{cluster_counter:03d}"
        for member, index in enumerate(indices, start=1):
            result[index] = replace(
                result[index],
                root_cluster_id=cluster_id,
                cluster_member_index=member,
                cluster_size=len(indices),
            )
    return tuple(replace(root, sorted_index=index) for index, root in enumerate(result, start=1))


def _global_candidates(
    evaluator: _MatrixEvaluator,
    settings: SearchSettings,
    *,
    configuration: str,
    scan_step: float,
    phases: Sequence[float],
    upper: float,
    seed_roots: Sequence[float],
    seed_source: str,
) -> tuple[list[RootCandidate], list[dict[str, object]], list[str]]:
    candidates: list[RootCandidate] = []
    intervals: list[tuple[float, float, str]] = []
    interval_rows: list[dict[str, object]] = []
    unresolved: list[str] = []
    for phase in phases:
        grid = _regular_grid(settings.lambda_min, upper, scan_step, phase)
        diagnostics = [evaluator.diagnostics(float(value)) for value in grid]
        if phase != 0.0:
            evaluator.operations.shifted_grid_evaluations += len(grid)
        if configuration == "verification":
            evaluator.operations.half_step_grid_evaluations += len(grid)
        determinants = np.asarray([item.determinant for item in diagnostics], dtype=float)
        sigma = np.asarray([item.sigma_1 for item in diagnostics], dtype=float)
        ratio = np.asarray([item.sigma_ratio for item in diagnostics], dtype=float)
        source_prefix = f"{configuration}_{'shifted' if phase else 'unshifted'}"
        for index in range(len(grid) - 1):
            left = float(grid[index])
            right = float(grid[index + 1])
            f_left = float(determinants[index])
            f_right = float(determinants[index + 1])
            if not (math.isfinite(f_left) and math.isfinite(f_right)):
                continue
            if f_left == 0.0:
                diag = diagnostics[index]
                on_scan_boundary = left <= settings.lambda_min + settings.root_match_tol
                status = (
                    "accepted_full_matrix_svd"
                    if not on_scan_boundary and _candidate_is_accepted(diag, settings)
                    else ("rejected_boundary_minimum" if on_scan_boundary else "rejected_full_matrix_svd")
                )
                candidates.append(
                    RootCandidate(
                        Lambda=left,
                        detection_sources=(source_prefix + ":grid_zero",),
                        diagnostics=diag,
                        interval_left=left,
                        interval_right=right,
                        interior_minimum=True,
                        acceptance_status=status,
                    )
                )
            elif f_left * f_right < 0.0:
                candidate = _brent_candidate(evaluator, left, right, source_prefix + ":sign_change", settings)
                if candidate is not None:
                    candidates.append(candidate)
        for index in range(1, len(grid) - 1):
            current = float(sigma[index])
            current_ratio = float(ratio[index])
            if not (math.isfinite(current) and math.isfinite(current_ratio)):
                continue
            sigma_min = current <= float(sigma[index - 1]) and current <= float(sigma[index + 1])
            ratio_min = current_ratio <= float(ratio[index - 1]) and current_ratio <= float(ratio[index + 1])
            deep = current <= settings.sigma_prefilter and current_ratio <= 10.0 * settings.sigma_ratio_accept
            if deep and (sigma_min or ratio_min):
                intervals.append(
                    (float(grid[index - 1]), float(grid[index + 1]), source_prefix + ":global_sigma_valley")
                )
    for seed in seed_roots:
        if not math.isfinite(float(seed)) or not settings.lambda_min < float(seed) < upper:
            continue
        half_width = max(settings.seed_half_width, 2.0 * scan_step)
        intervals.append(
            (
                max(settings.lambda_min, float(seed) - half_width),
                min(upper, float(seed) + half_width),
                seed_source,
            )
        )
        seed_diag = evaluator.diagnostics(float(seed))
        seed_status = (
            "accepted_full_matrix_svd"
            if _candidate_is_accepted(seed_diag, settings)
            else "rejected_seed_residual"
        )
        candidates.append(
            RootCandidate(
                Lambda=float(seed),
                detection_sources=(seed_source + ":direct_full_matrix_SVD",),
                diagnostics=seed_diag,
                interval_left=max(settings.lambda_min, float(seed) - half_width),
                interval_right=min(upper, float(seed) + half_width),
                interior_minimum=True,
                acceptance_status=seed_status,
                notes="seed is only accepted after its own full-matrix SVD check",
            )
        )
        if seed_source == "continuation_seed_window":
            evaluator.operations.continuation_seed_windows += 1
        elif seed_source == "EB_seed_window_for_Timo":
            evaluator.operations.EB_seed_windows_for_Timo += 1
    merged_intervals: list[tuple[float, float, set[str]]] = []
    for left, right, reason in sorted(intervals, key=lambda item: (item[0], item[1], item[2])):
        if merged_intervals and float(left) <= merged_intervals[-1][1] + settings.root_match_tol:
            old_left, old_right, reasons = merged_intervals[-1]
            reasons.add(str(reason))
            merged_intervals[-1] = (old_left, max(old_right, float(right)), reasons)
        else:
            merged_intervals.append((float(left), float(right), {str(reason)}))
    for left, right, reasons in merged_intervals:
        reason = "+".join(sorted(reasons))
        accepted_before = [
            item
            for item in candidates
            if item.acceptance_status == "accepted_full_matrix_svd"
            and float(left) - settings.root_match_tol <= item.Lambda <= float(right) + settings.root_match_tol
        ]
        seeds_inside = sum(float(left) <= float(seed) <= float(right) for seed in seed_roots)
        grid_zero_requires_pair_check = any(
            any("grid_zero" in source for source in item.detection_sources)
            for item in accepted_before
        )
        if (
            accepted_before
            and len(accepted_before) >= max(1, seeds_inside)
            and not grid_zero_requires_pair_check
        ):
            interval_rows.append(
                {
                    "Lambda_left": left,
                    "Lambda_right": right,
                    "determinant_sign_left": "",
                    "determinant_sign_right": "",
                    "sigma_left": "",
                    "sigma_midpoint": "",
                    "sigma_right": "",
                    "minimum_sigma": min(item.diagnostics.sigma_1 for item in accepted_before),
                    "minimum_sigma_ratio": min(item.diagnostics.sigma_ratio for item in accepted_before),
                    "subdivision_depth": 0,
                    "trigger_reasons": reason,
                    "local_minima_count": len(accepted_before),
                    "resolution_status": "resolved_existing_candidate",
                }
            )
            continue
        local_candidates, row = _adaptive_interval_candidates(evaluator, left, right, reason, settings)
        candidates.extend(local_candidates)
        accepted = [
            item
            for item in candidates
            if item.acceptance_status == "accepted_full_matrix_svd"
            and float(left) - settings.root_match_tol <= item.Lambda <= float(right) + settings.root_match_tol
        ]
        if accepted:
            row["resolution_status"] = "resolved_root_interval"
        elif (
            float(row["minimum_sigma"]) <= settings.sigma_accept
            and float(row["minimum_sigma_ratio"]) <= 10.0 * settings.sigma_ratio_accept
        ):
            row["resolution_status"] = "unresolved_low_sigma_interval"
            unresolved.append(f"{left:.12g}:{right:.12g}:{reason}")
        else:
            row["resolution_status"] = "closed_false_minimum"
        interval_rows.append(row)
    return candidates, interval_rows, unresolved


def _run_configuration(
    provider: MatrixProvider,
    settings: SearchSettings,
    *,
    configuration: str,
    candidate_target: int,
    scan_step: float,
    phases: Sequence[float],
    seed_roots: Sequence[float],
    seed_source: str,
) -> SearchConfigurationResult:
    operations = OperationCounts()
    if configuration == "verification":
        operations.independent_verification_runs += 1
    evaluator = _MatrixEvaluator(provider, operations)
    upper = _initial_upper(settings, candidate_target)
    selected_candidates: list[RootCandidate] = []
    interval_rows: list[dict[str, object]] = []
    unresolved: list[str] = []
    roots: tuple[RootRecord, ...] = ()
    for _attempt in range(settings.max_upper_growth_tries):
        candidates, rows, unresolved_rows = _global_candidates(
            evaluator,
            settings,
            configuration=configuration,
            scan_step=scan_step,
            phases=phases,
            upper=upper,
            seed_roots=seed_roots,
            seed_source=seed_source,
        )
        selected_candidates = _merge_candidates(candidates, settings)
        roots = _root_records(selected_candidates, settings)
        interval_rows = rows
        unresolved = unresolved_rows
        if len(roots) >= candidate_target:
            break
        upper *= settings.upper_growth_factor
    operations.accepted_recovered_roots = sum(
        not all("sign_change" in source for source in candidate.detection_sources)
        for candidate in selected_candidates
    )
    keep_count = int(candidate_target) + 4
    if len(roots) > keep_count:
        cutoff = roots[keep_count - 1].Lambda + settings.root_match_tol
        roots = tuple(roots[:keep_count])
        selected_candidates = [item for item in selected_candidates if item.Lambda <= cutoff]
    return SearchConfigurationResult(
        configuration=configuration,
        scan_step=float(scan_step),
        grid_phase=float(phases[0]) if len(phases) == 1 else -1.0,
        candidate_root_target=int(candidate_target),
        lambda_upper=float(upper),
        roots=roots,
        candidates=tuple(selected_candidates),
        interval_rows=tuple(interval_rows),
        unresolved_intervals=tuple(unresolved),
        operations=operations,
    )


def _compare_configurations(
    primary: SearchConfigurationResult,
    verification: SearchConfigurationResult,
    settings: SearchSettings,
) -> tuple[tuple[dict[str, object], ...], bool]:
    rows: list[dict[str, object]] = []
    agreement = len(primary.roots) >= settings.requested_roots and len(verification.roots) >= settings.requested_roots
    for index in range(settings.requested_roots):
        left = primary.roots[index] if index < len(primary.roots) else None
        right = verification.roots[index] if index < len(verification.roots) else None
        if left is None or right is None:
            rows.append(
                {
                    "sorted_index": index + 1,
                    "Lambda_primary": left.Lambda if left else "",
                    "Lambda_verification": right.Lambda if right else "",
                    "absolute_difference": "",
                    "self_MAC": "",
                    "multiplicity_agreement": False,
                    "status": "missing_root",
                }
            )
            agreement = False
            continue
        difference = abs(left.Lambda - right.Lambda)
        mac = coefficient_self_mac(left.null_vector, right.null_vector)
        multiplicity_agreement = left.detected_nullity == right.detected_nullity
        status_ok = difference <= settings.root_match_tol and multiplicity_agreement and (
            not math.isfinite(mac) or mac >= 0.90 or left.cluster_size > 1 or right.cluster_size > 1
        )
        if not status_ok:
            agreement = False
        rows.append(
            {
                "sorted_index": index + 1,
                "Lambda_primary": left.Lambda,
                "Lambda_verification": right.Lambda,
                "absolute_difference": difference,
                "self_MAC": mac,
                "multiplicity_agreement": multiplicity_agreement,
                "status": "pass" if status_ok else "disagreement",
            }
        )
    return tuple(rows), bool(agreement)


def resolve_matrix_spectrum(
    matrix_provider: MatrixProvider,
    *,
    settings: SearchSettings | None = None,
    primary_seeds: Sequence[float] = (),
    verification_seeds: Sequence[float] = (),
    primary_seed_source: str = "continuation_seed_window",
    verification_seed_source: str = "independent_seed_window",
    model: str = "synthetic",
    geometry: Geometry | None = None,
    cache_status: str = "not_requested",
) -> CompleteSpectrumResult:
    active = settings or SearchSettings()
    active.validate()
    synthetic_geometry = geometry or Geometry(1.0, 0.0, 0.0, 0.0)
    primary = _run_configuration(
        matrix_provider,
        active,
        configuration="primary",
        candidate_target=active.candidate_roots,
        scan_step=active.scan_step,
        phases=(0.0, active.shifted_grid_phase),
        seed_roots=primary_seeds,
        seed_source=primary_seed_source,
    )
    verification = _run_configuration(
        matrix_provider,
        active,
        configuration="verification",
        candidate_target=active.verification_candidate_roots,
        scan_step=0.5 * active.scan_step,
        phases=(active.shifted_grid_phase,),
        seed_roots=verification_seeds,
        seed_source=verification_seed_source,
    )
    comparison, agreement = _compare_configurations(primary, verification, active)
    roots = tuple(primary.roots[: active.requested_roots])
    root11 = len(roots) >= 11
    root12 = len(roots) >= 12
    boundary_warning = len(primary.roots) < active.candidate_roots or len(verification.roots) < active.verification_candidate_roots
    reasons: list[str] = []
    if len(roots) < active.requested_roots:
        reasons.append(f"found_only_{len(roots)}_of_{active.requested_roots}")
    if any(root.sigma_1 > active.sigma_accept for root in roots):
        reasons.append("full_matrix_SVD_failure")
    if not agreement:
        reasons.append("unresolved_independent_search_disagreement")
    if boundary_warning:
        reasons.append("candidate_boundary_warning")
    audit_limit = (
        roots[-1].Lambda + max(active.seed_half_width, 2.0 * active.scan_step)
        if roots
        else float("inf")
    )

    def relevant_interval(entry: str) -> bool:
        try:
            return float(str(entry).split(":", 1)[0]) <= audit_limit
        except (TypeError, ValueError):
            return True

    relevant_unresolved = [
        entry
        for entry in (*primary.unresolved_intervals, *verification.unresolved_intervals)
        if relevant_interval(entry)
    ]
    if relevant_unresolved:
        reasons.append("unresolved_low_sigma_interval")
    status = "resolved_complete" if not reasons else "unresolved"
    operations = OperationCounts()
    operations.add(primary.operations)
    operations.add(verification.operations)
    return CompleteSpectrumResult(
        algorithm_version=GENERAL_SPECTRUM_ALGORITHM_VERSION,
        model=str(model),
        geometry=synthetic_geometry,
        settings=active,
        primary=primary,
        verification=verification,
        roots=roots,
        primary_vs_verification=comparison,
        independent_agreement=agreement,
        root11_available=root11,
        root12_available=root12,
        root12_boundary_warning=boundary_warning,
        spectrum_status=status,
        exclusion_reason=";".join(dict.fromkeys(reasons)),
        cache_status=cache_status,
        operations=operations,
    )


def resolve_general_spectrum(
    model: str,
    geometry: Geometry,
    *,
    settings: SearchSettings | None = None,
    continuation_seeds: Sequence[float] = (),
    eb_seed_roots: Sequence[float] = (),
) -> CompleteSpectrumResult:
    active = settings or SearchSettings()
    active.validate()
    geometry.validate()
    if model not in SUPPORTED_MODELS:
        raise ValueError(f"unsupported model: {model}")
    seeds = list(float(value) for value in continuation_seeds if math.isfinite(float(value)))
    seed_source = "continuation_seed_window"
    if abs(geometry.beta_deg) <= 1.0e-14 and abs(geometry.eta) <= 1.0e-14:
        # This is an allowed oracle/seed use, not a hard-coded acceptance: each
        # seed is still accepted only after evaluation of the unchanged full
        # 6x6 matrix and the independent configuration does not reuse it.
        seeds.extend(straight_oracle_values(model, geometry, active.candidate_roots))
        seed_source = "straight_oracle_seed"
    if model == MODEL_TIMO and eb_seed_roots:
        seeds.extend(float(value) for value in eb_seed_roots if math.isfinite(float(value)))
        seed_source = "EB_seed_window_for_Timo"
    result = resolve_matrix_spectrum(
        model_matrix_provider(model, geometry),
        settings=active,
        primary_seeds=seeds,
        verification_seeds=(),
        primary_seed_source=seed_source,
        verification_seed_source="independent_no_primary_objects",
        model=model,
        geometry=geometry,
    )
    return annotate_coalesced_tracks(result, continuation_seeds, active)


def annotate_coalesced_tracks(
    result: CompleteSpectrumResult,
    continuation_seeds: Sequence[float],
    settings: SearchSettings | None = None,
) -> CompleteSpectrumResult:
    active = settings or result.settings
    if result.primary.roots:
        # Track annotations describe the current continuation request and are
        # therefore not persistent spectral facts.  Strip a cached annotation
        # before assigning the current seeds; exact SVD nullity metadata is
        # intentionally left untouched.
        updated_primary = [
            replace(
                root,
                track_multiplicity=1,
                multiplicity_status="simple_root",
                notes="",
            )
            if root.multiplicity_status == "coalesced_track_cluster"
            else root
            for root in result.primary.roots
        ]
        assignments: dict[int, int] = {}
        for seed in continuation_seeds:
            if not math.isfinite(float(seed)):
                continue
            index = min(
                range(len(updated_primary)),
                key=lambda item: abs(updated_primary[item].Lambda - float(seed)),
            )
            if abs(updated_primary[index].Lambda - float(seed)) <= active.seed_half_width:
                assignments[index] = assignments.get(index, 0) + 1
        for index, track_count in assignments.items():
            root = updated_primary[index]
            if track_count >= 2 and root.detected_nullity == 1:
                updated_primary[index] = replace(
                    root,
                    track_multiplicity=track_count,
                    multiplicity_status="coalesced_track_cluster",
                    notes="multiple continuation tracks converge; algebraic multiplicity is not asserted",
                )
        updated_first = tuple(updated_primary[: active.requested_roots])
        result = replace(
            result,
            primary=replace(result.primary, roots=tuple(updated_primary)),
            roots=updated_first,
        )
    return result


def straight_oracle_values(model: str, geometry: Geometry, count: int) -> tuple[float, ...]:
    if abs(geometry.beta_deg) > 1.0e-14 or abs(geometry.eta) > 1.0e-14:
        raise ValueError("straight oracle is available only at beta=0, eta=0")
    if model == MODEL_TIMO:
        return tuple(straight_factorized.factorized_straight_spectrum(
            geometry.epsilon_0, geometry.mu, int(count)
        ).values)
    if model == MODEL_EB:
        axial = np.sqrt(np.arange(1, int(count) + 1, dtype=float) * np.pi / (2.0 * geometry.epsilon_0))
        bending = np.asarray(fixed_fixed_eb_roots(int(count)), dtype=float) / 2.0
        return tuple(sorted([*(float(value) for value in axial), *(float(value) for value in bending)])[: int(count)])
    raise ValueError(f"unsupported model: {model}")


def _result_to_payload(result: CompleteSpectrumResult, identity: Mapping[str, object]) -> dict[str, object]:
    return {
        "algorithm_version": GENERAL_SPECTRUM_ALGORITHM_VERSION,
        "identity": dict(identity),
        "result": asdict(result),
    }


def _matrix_diagnostics_from_dict(data: Mapping[str, object]) -> MatrixDiagnostics:
    return MatrixDiagnostics(
        Lambda=float(data["Lambda"]),
        determinant=float(data["determinant"]),
        sigma_1=float(data["sigma_1"]),
        sigma_2=float(data["sigma_2"]),
        sigma_3=float(data.get("sigma_3", float("inf"))),
        sigma_ratio=float(data["sigma_ratio"]),
        null_vector=tuple(float(value) for value in data.get("null_vector", ())),
        null_vector_pivot=int(data["null_vector_pivot"]),
        finite_matrix_status=str(data["finite_matrix_status"]),
    )


def _candidate_from_dict(data: Mapping[str, object]) -> RootCandidate:
    return RootCandidate(
        Lambda=float(data["Lambda"]),
        detection_sources=tuple(str(value) for value in data.get("detection_sources", ())),
        diagnostics=_matrix_diagnostics_from_dict(data["diagnostics"]),  # type: ignore[arg-type]
        interval_left=float(data["interval_left"]),
        interval_right=float(data["interval_right"]),
        interior_minimum=bool(data["interior_minimum"]),
        acceptance_status=str(data["acceptance_status"]),
        notes=str(data.get("notes", "")),
    )


def _root_from_dict(data: Mapping[str, object]) -> RootRecord:
    return RootRecord(
        Lambda=float(data["Lambda"]),
        sorted_index=int(data["sorted_index"]),
        root_cluster_id=str(data.get("root_cluster_id", "")),
        cluster_member_index=int(data.get("cluster_member_index", 1)),
        cluster_size=int(data.get("cluster_size", 1)),
        detection_sources=tuple(str(value) for value in data.get("detection_sources", ())),
        determinant=float(data["determinant"]),
        sigma_1=float(data["sigma_1"]),
        sigma_2=float(data["sigma_2"]),
        sigma_ratio=float(data["sigma_ratio"]),
        null_vector=tuple(float(value) for value in data.get("null_vector", ())),
        null_vector_pivot=int(data["null_vector_pivot"]),
        self_MAC=float(data.get("self_MAC", 1.0)),
        detected_nullity=int(data.get("detected_nullity", 1)),
        track_multiplicity=int(data.get("track_multiplicity", 1)),
        multiplicity_status=str(data.get("multiplicity_status", "simple_root")),
        acceptance_status=str(data["acceptance_status"]),
        notes=str(data.get("notes", "")),
    )


def _operations_from_dict(data: Mapping[str, object]) -> OperationCounts:
    kwargs = {name: int(data.get(name, 0)) for name in OperationCounts.__dataclass_fields__}
    return OperationCounts(**kwargs)


def _configuration_from_dict(data: Mapping[str, object]) -> SearchConfigurationResult:
    return SearchConfigurationResult(
        configuration=str(data["configuration"]),
        scan_step=float(data["scan_step"]),
        grid_phase=float(data["grid_phase"]),
        candidate_root_target=int(data["candidate_root_target"]),
        lambda_upper=float(data["lambda_upper"]),
        roots=tuple(_root_from_dict(item) for item in data.get("roots", ())),  # type: ignore[arg-type]
        candidates=tuple(_candidate_from_dict(item) for item in data.get("candidates", ())),  # type: ignore[arg-type]
        interval_rows=tuple(dict(item) for item in data.get("interval_rows", ())),  # type: ignore[arg-type]
        unresolved_intervals=tuple(str(value) for value in data.get("unresolved_intervals", ())),
        operations=_operations_from_dict(data.get("operations", {})),  # type: ignore[arg-type]
    )


def _result_from_payload(payload: Mapping[str, object], cache_status: str) -> CompleteSpectrumResult:
    data = payload["result"]
    if not isinstance(data, Mapping):
        raise ValueError("invalid cache result payload")
    settings = SearchSettings(**dict(data["settings"]))  # type: ignore[arg-type]
    geometry = Geometry(**dict(data["geometry"]))  # type: ignore[arg-type]
    return CompleteSpectrumResult(
        algorithm_version=str(data["algorithm_version"]),
        model=str(data["model"]),
        geometry=geometry,
        settings=settings,
        primary=_configuration_from_dict(data["primary"]),  # type: ignore[arg-type]
        verification=_configuration_from_dict(data["verification"]),  # type: ignore[arg-type]
        roots=tuple(_root_from_dict(item) for item in data.get("roots", ())),  # type: ignore[arg-type]
        primary_vs_verification=tuple(dict(item) for item in data.get("primary_vs_verification", ())),  # type: ignore[arg-type]
        independent_agreement=bool(data["independent_agreement"]),
        root11_available=bool(data["root11_available"]),
        root12_available=bool(data["root12_available"]),
        root12_boundary_warning=bool(data["root12_boundary_warning"]),
        spectrum_status=str(data["spectrum_status"]),
        exclusion_reason=str(data.get("exclusion_reason", "")),
        cache_status=cache_status,
        operations=_operations_from_dict(data.get("operations", {})),  # type: ignore[arg-type]
    )


class GeneralSpectrumCache:
    def __init__(self, cache_dir: Path, *, reuse_cache: bool = True, force_recompute: bool = False) -> None:
        self.cache_dir = Path(cache_dir)
        self.reuse_cache = bool(reuse_cache)
        self.force_recompute = bool(force_recompute)
        self.last_load_status = "not_checked"

    @staticmethod
    def identity(model: str, geometry: Geometry, settings: SearchSettings) -> dict[str, object]:
        # The filename intentionally excludes the algorithm version so that a
        # cache from an older algorithm is found and rejected explicitly.
        return {
            "model": str(model),
            "geometry": asdict(geometry),
            "requested_root_count": settings.requested_roots,
            "candidate_root_count": settings.candidate_roots,
            "verification_candidate_root_count": settings.verification_candidate_roots,
            "lambda_scan_bounds": [settings.lambda_min, settings.lambda_max],
            "base_scan_step": settings.scan_step,
            "shifted_grid_phase": settings.shifted_grid_phase,
            "local_refinement_settings": [settings.local_points, settings.adaptive_depth],
            "SVD_acceptance_settings": [
                settings.sigma_prefilter,
                settings.sigma_accept,
                settings.sigma_ratio_accept,
                settings.nullity_sigma,
            ],
            "root_tolerances": [settings.root_dedup_tol, settings.root_match_tol],
        }

    def path(self, model: str, geometry: Geometry, settings: SearchSettings) -> Path:
        encoded = json.dumps(self.identity(model, geometry, settings), sort_keys=True, separators=(",", ":"))
        digest = hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:24]
        token = "eb" if model == MODEL_EB else "timo"
        return self.cache_dir / f"general_spectrum_{token}_{digest}.json"

    def load(self, model: str, geometry: Geometry, settings: SearchSettings) -> CompleteSpectrumResult | None:
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
        if payload.get("algorithm_version") != GENERAL_SPECTRUM_ALGORITHM_VERSION:
            self.last_load_status = "stale_cache_algorithm_version"
            return None
        if payload.get("identity") != self.identity(model, geometry, settings):
            self.last_load_status = "stale_cache_settings"
            return None
        self.last_load_status = "hit"
        return _result_from_payload(payload, "hit")

    def save(self, result: CompleteSpectrumResult) -> Path:
        path = self.path(result.model, result.geometry, result.settings)
        path.parent.mkdir(parents=True, exist_ok=True)
        identity = self.identity(result.model, result.geometry, result.settings)
        payload = _result_to_payload(result, identity)
        path.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")), encoding="utf-8")
        return path

    def resolve(
        self,
        model: str,
        geometry: Geometry,
        settings: SearchSettings,
        *,
        continuation_seeds: Sequence[float] = (),
        eb_seed_roots: Sequence[float] = (),
    ) -> CompleteSpectrumResult:
        cached = self.load(model, geometry, settings)
        if cached is not None:
            operations = replace(cached.operations)
            operations.cache_hits += 1
            cached = replace(cached, operations=operations, cache_status="hit")
            return annotate_coalesced_tracks(cached, continuation_seeds, settings)
        load_status = self.last_load_status
        result = resolve_general_spectrum(
            model,
            geometry,
            settings=settings,
            continuation_seeds=continuation_seeds,
            eb_seed_roots=eb_seed_roots,
        )
        operations = replace(result.operations)
        operations.cache_misses += 1
        result = replace(result, operations=operations, cache_status=load_status)
        self.save(result)
        return result


__all__ = [
    "AUTO_SPECTRUM_ALGORITHM_VERSION",
    "CompleteSpectrumResult",
    "GENERAL_SPECTRUM_ALGORITHM_VERSION",
    "GeneralSpectrumCache",
    "Geometry",
    "MODEL_EB",
    "MODEL_TIMO",
    "MatrixDiagnostics",
    "OperationCounts",
    "RootCandidate",
    "RootRecord",
    "SearchConfigurationResult",
    "SearchSettings",
    "coefficient_self_mac",
    "annotate_coalesced_tracks",
    "model_matrix_provider",
    "resolve_general_spectrum",
    "resolve_matrix_spectrum",
    "row_normalized",
    "straight_oracle_values",
]
