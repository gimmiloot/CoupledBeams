"""RLB-1 diagnostic pilot for one rigid joint in the straight ``beta=0`` limit.

The workflow is intentionally narrow.  It derives the joint from the frozen
RLB coordinate contract, assembles two physical Reddy-beam transfers, and
searches only the straight-limit spectrum.  Nonzero angles are used solely in
algebraic joint and virtual-work checks.  No two-arm Ritz discretization,
finite-element model, parameter sweep, or figure generation is implemented.

The model manifest is written before any joint, random, or spectral
calculation.  The primary coupled inventories are seed-free and are frozen
before the independent fixed--fixed and stepped references are constructed.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass, field, replace
import hashlib
import importlib
import json
import math
from pathlib import Path
import subprocess
import sys
from typing import Any, Callable, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import brentq, minimize_scalar


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.analysis.laminated_beams import (  # noqa: E402
    audit_reddy_table_4_3_3_discrepancies as source_audit,
)
from scripts.lib import reddy_symmetric_laminated_beam as rlb  # noqa: E402


FloatArray = NDArray[np.float64]
MatrixProvider = Callable[[float], FloatArray]

DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_coupled_beta0_pilot"
)

ALGORITHM_VERSION = "rlb_1_beta0_seed_free_inventory_v3"
RANDOM_SEED_JOINT = 2026082501
RANDOM_SEED_VIRTUAL_WORK = 2026082502
ORIGINAL_PREIMPLEMENTATION_MANIFEST_SHA256 = (
    "B77E28DB87BD721F7A92DAD57DABBC80BB870845771C9EE04AA3D26DB90EB071"
)
TASK_INITIAL_GIT_STATE = {
    "top_level": "D:/PHD/CoupledBeams/CoupledBeams",
    "branch": "main",
    "head": "10ab253456604be15676f9b8ae5a4ece97b200b2",
    "last_commit": "10ab253 Version 0.4.0.1",
    "status_short": "",
    "provenance": "captured_by_the_mandatory_pre_edit_git_audit",
}

STATUS_JOINT = "RLB-1J-JOINT-THEORY"
STATUS_HOMOGENEOUS = "RLB-1A-BETA0-HOMOGENEOUS"
STATUS_STEPPED = "RLB-1B-BETA0-STEPPED"
STATUS_INVENTORY = "RLB-1-BETA0-ROOT-INVENTORY"
STATUS_OVERALL = "OVERALL"

THRESHOLDS: dict[str, float] = {
    "coordinate_J_equality": 1.0e-14,
    "rank_singular_relative": 1.0e-12,
    "virtual_work_normalized_residual": 1.0e-12,
    "beta0_isolated_spectral_relative": 1.0e-9,
    "cluster_center_relative": 1.0e-9,
    "joint_compatibility_residual": 1.0e-10,
    "joint_equilibrium_residual": 1.0e-9,
    "boundary_residual": 1.0e-9,
    "energy_identity": 1.0e-8,
    "root_singular_ratio": 1.0e-9,
}


def _finite_positive(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _relative_difference(left: float, right: float) -> float:
    return abs(float(left) - float(right)) / max(
        abs(float(left)), abs(float(right)), np.finfo(float).tiny
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _json_value(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value).replace("\\", "/")
    if isinstance(value, tuple):
        return [_json_value(item) for item in value]
    if isinstance(value, list):
        return [_json_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, float) and not math.isfinite(value):
        return "inf" if value > 0.0 else ("-inf" if value < 0.0 else "nan")
    return value


def _write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_json_value(dict(payload)), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _write_csv(path: Path, rows: Iterable[Mapping[str, Any]], fields: Sequence[str] | None = None) -> None:
    data = [dict(row) for row in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        ordered: list[str] = []
        for row in data:
            for key in row:
                if key not in ordered:
                    ordered.append(key)
        fields = ordered or ("status",)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields), extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        for row in data:
            writer.writerow({key: _csv_value(row.get(key, "")) for key in fields})


def _csv_value(value: Any) -> Any:
    if isinstance(value, (tuple, list, dict, np.ndarray)):
        return json.dumps(_json_value(value), ensure_ascii=False, separators=(",", ":"))
    if isinstance(value, np.generic):
        return value.item()
    return value


@dataclass(frozen=True)
class SearchPolicy:
    """Frozen seed-free scan and root-quality settings."""

    requested_roots: int = 12
    guard_roots: int = 1
    omega_bar_min: float = 1.0e-8
    omega_bar_max: float = 1600.0
    primary_scan_points: int = 8001
    verification_scan_points: int = 16001
    primary_phases: tuple[float, ...] = (0.0,)
    verification_phases: tuple[float, ...] = (0.5,)
    sigma_prefilter: float = 1.0e-5
    sigma_ratio_tolerance: float = THRESHOLDS["root_singular_ratio"]
    rank_relative_tolerance: float = THRESHOLDS["rank_singular_relative"]
    boundary_residual_tolerance: float = THRESHOLDS["boundary_residual"]
    root_xtol_bar: float = 1.0e-11
    root_rtol: float = 8.0 * np.finfo(float).eps
    dedup_atol_bar: float = 5.0e-10
    dedup_rtol: float = 5.0e-12
    cluster_atol_bar: float = 1.0e-10
    cluster_rtol: float = 1.0e-10
    inventory_match_relative: float = THRESHOLDS["beta0_isolated_spectral_relative"]
    post_guard_tail_bar: float = 2.0
    reference_detector_reconciliation_atol_bar: float = 2.0e-10
    reference_detector_reconciliation_rtol: float = 5.0e-8

    def __post_init__(self) -> None:
        if self.requested_roots != 12 or self.guard_roots != 1:
            raise ValueError("RLB-1 requires exactly 12 roots and root 13 as guard")
        if not 0.0 <= self.omega_bar_min < self.omega_bar_max:
            raise ValueError("invalid omega_bar scan interval")
        if self.primary_scan_points < 33 or self.verification_scan_points <= self.primary_scan_points:
            raise ValueError("verification scan must be finer than the primary scan")
        for name in (
            "sigma_prefilter", "sigma_ratio_tolerance", "rank_relative_tolerance",
            "boundary_residual_tolerance", "root_xtol_bar", "root_rtol",
            "dedup_atol_bar", "dedup_rtol", "cluster_atol_bar", "cluster_rtol",
            "inventory_match_relative", "post_guard_tail_bar",
            "reference_detector_reconciliation_atol_bar",
            "reference_detector_reconciliation_rtol",
        ):
            _finite_positive(getattr(self, name), name)
        for phase in (*self.primary_phases, *self.verification_phases):
            if not 0.0 <= float(phase) < 1.0:
                raise ValueError("grid phases must lie in [0,1)")

    @property
    def required_slots(self) -> int:
        return self.requested_roots + self.guard_roots


@dataclass(frozen=True)
class PositiveEquilibration:
    matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray


@dataclass(frozen=True)
class MatrixDiagnostics:
    omega_bar: float
    omega: float
    raw_matrix: FloatArray
    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray
    raw_determinant: float
    raw_det_sign: float
    raw_logabsdet: float
    scaled_determinant: float
    scaled_det_sign: float
    scaled_logabsdet: float
    raw_singular_values: FloatArray
    scaled_singular_values: FloatArray
    raw_sigma_min: float
    raw_sigma_max: float
    raw_sigma_ratio: float
    scaled_sigma_min: float
    scaled_sigma_max: float
    scaled_sigma_ratio: float
    raw_condition_number: float
    scaled_condition_number: float
    detected_nullity: int
    root_gate_nullity: int
    scaled_null_vector: FloatArray
    physical_null_vector: FloatArray
    scaled_null_residual: float
    raw_boundary_null_residual: float
    finite: bool


@dataclass(frozen=True)
class RootCandidate:
    case_id: str
    builder_id: str
    scan_id: str
    omega_bar: float
    detection_sources: tuple[str, ...]
    interval_left_bar: float
    interval_right_bar: float
    interior_minimum: bool
    diagnostics: MatrixDiagnostics
    accepted: bool
    rejection_reason: str
    canonical: bool = True
    merge_group_id: str = ""


@dataclass(frozen=True)
class RootEvent:
    event_id: str
    omega_bar: float
    omega: float
    multiplicity: int
    detected_nullity: int
    candidate: RootCandidate
    cluster_id: str = ""
    cluster_semantics: str = "ISOLATED"
    cluster_multiplicity: int = 1
    cluster_total_nullity: int = 1
    cluster_center_omega_bar: float = math.nan


@dataclass(frozen=True)
class SpectrumSlot:
    sorted_slot: int
    role: str
    repeated_root_slot: int
    event: RootEvent


@dataclass(frozen=True)
class ScanResult:
    scan_id: str
    candidates: tuple[RootCandidate, ...]
    rejected_candidates: tuple[RootCandidate, ...]
    events: tuple[RootEvent, ...]
    slots: tuple[SpectrumSlot, ...]


@dataclass(frozen=True)
class RootInventory:
    case_id: str
    builder_id: str
    frequency_scale: float
    policy: SearchPolicy
    primary: ScanResult
    verification: ScanResult
    slots: tuple[SpectrumSlot, ...]
    independent_agreement: bool
    maximum_primary_verification_relative: float
    unresolved_low_sigma_count: int
    guard_available: bool
    guard_not_at_scan_boundary: bool
    status: str
    inventory_sha256: str


def positive_equilibrate(matrix: FloatArray) -> PositiveEquilibration:
    """Apply only positive finite row/column factors to a square matrix."""

    scaling = _coupled_module().positively_equilibrate_matrix(matrix)
    return PositiveEquilibration(
        np.asarray(scaling.scaled_matrix, dtype=float),
        np.asarray(scaling.row_factors, dtype=float),
        np.asarray(scaling.column_factors, dtype=float),
    )


def _singular_ratio(values: FloatArray) -> float:
    return float(values[-1] / values[0]) if values.size and values[0] > 0.0 else 0.0


def boundary_matrix_diagnostics(
    omega_bar: float,
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    *,
    rank_relative_tolerance: float = THRESHOLDS["rank_singular_relative"],
    root_ratio_tolerance: float = THRESHOLDS["root_singular_ratio"],
) -> MatrixDiagnostics:
    """Evaluate raw/scaled determinants, SVDs, nullity, and null residuals."""

    scale = _finite_positive(frequency_scale, "frequency_scale")
    omega_bar_value = float(omega_bar)
    omega = omega_bar_value / scale
    raw = np.asarray(matrix_provider(omega), dtype=float)
    equilibration = positive_equilibrate(raw)
    scaled = equilibration.matrix
    raw_sign, raw_log = np.linalg.slogdet(raw)
    scaled_sign, scaled_log = np.linalg.slogdet(scaled)
    raw_singular = np.linalg.svd(raw, compute_uv=False)
    _u, scaled_singular, vh = np.linalg.svd(scaled, full_matrices=True)
    scaled_vector = np.asarray(vh[-1], dtype=float)
    physical_vector = equilibration.column_factors * scaled_vector
    physical_norm = float(np.linalg.norm(physical_vector))
    if physical_norm > 0.0:
        physical_vector = physical_vector / physical_norm
    scaled_residual = float(
        np.linalg.norm(scaled @ scaled_vector)
        / max(np.linalg.norm(scaled) * np.linalg.norm(scaled_vector), np.finfo(float).tiny)
    )
    raw_residual = float(
        np.linalg.norm(raw @ physical_vector)
        / max(np.linalg.norm(raw) * np.linalg.norm(physical_vector), np.finfo(float).tiny)
    )
    scaled_max = float(scaled_singular[0])
    relative_singular = (
        scaled_singular / scaled_max if scaled_max > 0.0 else np.zeros_like(scaled_singular)
    )
    strict_nullity = int(np.count_nonzero(relative_singular <= rank_relative_tolerance))
    gate_nullity = int(np.count_nonzero(relative_singular <= root_ratio_tolerance))
    raw_min, raw_max = float(raw_singular[-1]), float(raw_singular[0])
    scaled_min = float(scaled_singular[-1])
    return MatrixDiagnostics(
        omega_bar=omega_bar_value,
        omega=omega,
        raw_matrix=raw,
        scaled_matrix=scaled,
        row_factors=equilibration.row_factors,
        column_factors=equilibration.column_factors,
        raw_determinant=float(np.linalg.det(raw)),
        raw_det_sign=float(raw_sign),
        raw_logabsdet=float(raw_log),
        scaled_determinant=float(np.linalg.det(scaled)),
        scaled_det_sign=float(scaled_sign),
        scaled_logabsdet=float(scaled_log),
        raw_singular_values=np.asarray(raw_singular, dtype=float),
        scaled_singular_values=np.asarray(scaled_singular, dtype=float),
        raw_sigma_min=raw_min,
        raw_sigma_max=raw_max,
        raw_sigma_ratio=(raw_min / raw_max if raw_max > 0.0 else 0.0),
        scaled_sigma_min=scaled_min,
        scaled_sigma_max=scaled_max,
        scaled_sigma_ratio=(scaled_min / scaled_max if scaled_max > 0.0 else 0.0),
        raw_condition_number=(raw_max / raw_min if raw_min > 0.0 else math.inf),
        scaled_condition_number=(scaled_max / scaled_min if scaled_min > 0.0 else math.inf),
        detected_nullity=strict_nullity,
        root_gate_nullity=gate_nullity,
        scaled_null_vector=scaled_vector,
        physical_null_vector=physical_vector,
        scaled_null_residual=scaled_residual,
        raw_boundary_null_residual=raw_residual,
        finite=bool(
            np.all(np.isfinite(raw))
            and np.all(np.isfinite(scaled))
            and np.all(np.isfinite(scaled_singular))
        ),
    )


class _DiagnosticEvaluator:
    def __init__(self, provider: MatrixProvider, frequency_scale: float, policy: SearchPolicy) -> None:
        self.provider = provider
        self.frequency_scale = frequency_scale
        self.policy = policy
        self.cache: dict[float, MatrixDiagnostics] = {}

    def diagnostics(self, omega_bar: float) -> MatrixDiagnostics:
        key = float(omega_bar)
        if key not in self.cache:
            self.cache[key] = boundary_matrix_diagnostics(
                key,
                self.provider,
                self.frequency_scale,
                rank_relative_tolerance=self.policy.rank_relative_tolerance,
                root_ratio_tolerance=self.policy.sigma_ratio_tolerance,
            )
        return self.cache[key]

    def determinant_objective(self, omega_bar: float) -> float:
        diagnostic = self.diagnostics(omega_bar)
        if diagnostic.scaled_det_sign == 0.0:
            return 0.0
        if not math.isfinite(diagnostic.scaled_logabsdet):
            return math.nan
        dimension = diagnostic.scaled_matrix.shape[0]
        return float(
            diagnostic.scaled_det_sign
            * math.exp(diagnostic.scaled_logabsdet / dimension)
        )

    def sigma_ratio(self, omega_bar: float) -> float:
        return self.diagnostics(omega_bar).scaled_sigma_ratio


def _candidate_quality(diagnostic: MatrixDiagnostics, policy: SearchPolicy) -> tuple[bool, str]:
    if not diagnostic.finite:
        return False, "NONFINITE_MATRIX"
    if diagnostic.scaled_sigma_ratio > policy.sigma_ratio_tolerance:
        return False, "ROOT_SINGULAR_RATIO_FAIL"
    if diagnostic.scaled_null_residual > policy.boundary_residual_tolerance:
        return False, "BOUNDARY_NULL_RESIDUAL_FAIL"
    if diagnostic.raw_boundary_null_residual > policy.boundary_residual_tolerance:
        return False, "PHYSICAL_RAW_BOUNDARY_NULL_RESIDUAL_FAIL"
    if diagnostic.detected_nullity < 1:
        return False, "NULLITY_UNRESOLVED_AT_1E-12"
    return True, ""


def _grid(policy: SearchPolicy, *, points: int, phase: float) -> FloatArray:
    step = (policy.omega_bar_max - policy.omega_bar_min) / (points - 1)
    start = policy.omega_bar_min + phase * step
    count = int(math.floor((policy.omega_bar_max - start) / step)) + 1
    return start + step * np.arange(count, dtype=float)


def _root_candidate(
    *,
    evaluator: _DiagnosticEvaluator,
    policy: SearchPolicy,
    case_id: str,
    builder_id: str,
    scan_id: str,
    source: str,
    left: float,
    right: float,
    omega_bar: float,
    interior: bool,
) -> RootCandidate:
    polished = float(omega_bar)
    # A determinant bracket can locate a very narrow root well while its
    # absolute/relative stopping rule still leaves sigma_min just above the
    # stricter 1e-12 nullity threshold.  Refine in a local offset coordinate;
    # using the absolute frequency directly would reintroduce a scale-dependent
    # sqrt(eps)*|omega| termination error.
    local_half_width = max(
        1.0e-4 * (right - left),
        64.0 * policy.root_xtol_bar,
        256.0 * abs(float(np.spacing(polished))),
    )
    delta_left = max(left - polished, -local_half_width)
    delta_right = min(right - polished, local_half_width)
    if delta_left < delta_right:
        try:
            local_fit = minimize_scalar(
                lambda delta: evaluator.sigma_ratio(polished + float(delta)) ** 2,
                bounds=(delta_left, delta_right),
                method="bounded",
                options={
                    "xatol": max(
                        np.finfo(float).tiny,
                        policy.root_xtol_bar * 1.0e-4,
                    ),
                    "maxiter": 240,
                },
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError):
            local_fit = None
        if local_fit is not None and local_fit.success and math.isfinite(float(local_fit.x)):
            fitted = polished + float(local_fit.x)
            if evaluator.sigma_ratio(fitted) <= evaluator.sigma_ratio(polished):
                polished = fitted
    # Absolute-frequency solvers stop at a scale-dependent tolerance.  Inspect
    # adjacent representable numbers so a root that is exactly representable
    # is not rejected merely because the returned value is one ULP away.
    for _iteration in range(8):
        neighbours = (
            polished,
            float(np.nextafter(polished, -math.inf)),
            float(np.nextafter(polished, math.inf)),
        )
        bounded = [value for value in neighbours if left <= value <= right]
        best = min(bounded, key=lambda value: evaluator.sigma_ratio(value))
        if best == polished:
            break
        polished = best
    diagnostic = evaluator.diagnostics(polished)
    accepted, reason = _candidate_quality(diagnostic, policy)
    if not interior:
        accepted, reason = False, "BOUNDARY_MINIMUM"
    return RootCandidate(
        case_id=case_id,
        builder_id=builder_id,
        scan_id=scan_id,
        omega_bar=polished,
        detection_sources=(source,),
        interval_left_bar=float(left),
        interval_right_bar=float(right),
        interior_minimum=bool(interior),
        diagnostics=diagnostic,
        accepted=accepted,
        rejection_reason=reason,
    )


def _scan_candidates(
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
    *,
    case_id: str,
    builder_id: str,
    scan_id: str,
    points: int,
    phases: Sequence[float],
) -> tuple[list[RootCandidate], _DiagnosticEvaluator]:
    evaluator = _DiagnosticEvaluator(matrix_provider, frequency_scale, policy)
    candidates: list[RootCandidate] = []
    for phase in phases:
        grid = _grid(policy, points=points, phase=float(phase))
        diagnostics = [evaluator.diagnostics(float(value)) for value in grid]
        determinants = np.asarray(
            [evaluator.determinant_objective(float(value)) for value in grid], dtype=float
        )
        sigma = np.asarray([item.scaled_sigma_ratio for item in diagnostics], dtype=float)
        source_prefix = f"{scan_id}_phase_{phase:g}"
        for index in range(grid.size - 1):
            left, right = float(grid[index]), float(grid[index + 1])
            f_left, f_right = float(determinants[index]), float(determinants[index + 1])
            if not (math.isfinite(f_left) and math.isfinite(f_right)):
                candidates.append(
                    RootCandidate(
                        case_id, builder_id, scan_id, left,
                        (source_prefix + ":determinant_interval",), left, right,
                        True, diagnostics[index], False, "NONFINITE_DETERMINANT_INTERVAL",
                    )
                )
                continue
            if f_left == 0.0:
                candidates.append(
                    _root_candidate(
                        evaluator=evaluator, policy=policy, case_id=case_id,
                        builder_id=builder_id, scan_id=scan_id,
                        source=source_prefix + ":determinant_grid_zero",
                        left=left, right=right, omega_bar=left,
                        interior=left > policy.omega_bar_min + policy.root_xtol_bar,
                    )
                )
            elif f_left * f_right < 0.0:
                try:
                    root = float(
                        brentq(
                            evaluator.determinant_objective,
                            left,
                            right,
                            xtol=policy.root_xtol_bar,
                            rtol=policy.root_rtol,
                            maxiter=180,
                        )
                    )
                except (RuntimeError, ValueError, FloatingPointError):
                    diagnostic = diagnostics[index]
                    candidates.append(
                        RootCandidate(
                            case_id, builder_id, scan_id, left,
                            (source_prefix + ":determinant_bracket",), left, right,
                            True, diagnostic, False, "BRENT_FAILURE",
                        )
                    )
                else:
                    candidates.append(
                        _root_candidate(
                            evaluator=evaluator, policy=policy, case_id=case_id,
                            builder_id=builder_id, scan_id=scan_id,
                            source=source_prefix + ":determinant_bracket",
                            left=left, right=right, omega_bar=root, interior=True,
                        )
                    )
        for index in range(1, grid.size - 1):
            current = float(sigma[index])
            if not math.isfinite(current):
                continue
            if current > float(sigma[index - 1]) or current > float(sigma[index + 1]):
                continue
            left, right = float(grid[index - 1]), float(grid[index + 1])
            try:
                first = minimize_scalar(
                    lambda value: evaluator.sigma_ratio(float(value)) ** 2,
                    bounds=(left, right),
                    method="bounded",
                    options={"xatol": policy.root_xtol_bar, "maxiter": 180},
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                candidates.append(
                    RootCandidate(
                        case_id, builder_id, scan_id, float(grid[index]),
                        (source_prefix + ":sigma_ratio_minimum",), left, right,
                        True, diagnostics[index], False, "MINIMIZER_EXCEPTION",
                    )
                )
                continue
            if not first.success or not math.isfinite(float(first.x)):
                candidates.append(
                    RootCandidate(
                        case_id, builder_id, scan_id, float(grid[index]),
                        (source_prefix + ":sigma_ratio_minimum",), left, right,
                        True, diagnostics[index], False, "MINIMIZER_FAILURE",
                    )
                )
                continue
            root = float(first.x)
            half = max((right - left) / 8.0, 8.0 * policy.root_xtol_bar)
            narrow_left, narrow_right = max(left, root - half), min(right, root + half)
            delta_left = narrow_left - root
            delta_right = narrow_right - root
            try:
                second = minimize_scalar(
                    lambda delta: evaluator.sigma_ratio(root + float(delta)) ** 2,
                    bounds=(delta_left, delta_right),
                    method="bounded",
                    options={"xatol": policy.root_xtol_bar, "maxiter": 180},
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError):
                candidates.append(
                    RootCandidate(
                        case_id, builder_id, scan_id, root,
                        (source_prefix + ":sigma_ratio_minimum",), left, right,
                        True, evaluator.diagnostics(root), False, "NESTED_MINIMIZER_EXCEPTION",
                    )
                )
                continue
            if second.success and math.isfinite(float(second.fun)) and float(second.fun) <= float(first.fun):
                root = root + float(second.x)
            edge = max(4.0 * policy.root_xtol_bar, 1.0e-5 * (right - left))
            interior = left + edge < root < right - edge
            candidate = _root_candidate(
                evaluator=evaluator, policy=policy, case_id=case_id,
                builder_id=builder_id, scan_id=scan_id,
                source=source_prefix + ":sigma_ratio_minimum",
                left=left, right=right, omega_bar=root, interior=interior,
            )
            if not candidate.accepted and candidate.rejection_reason == "ROOT_SINGULAR_RATIO_FAIL":
                candidate = replace(candidate, rejection_reason="FALSE_SIGMA_VALLEY")
            candidates.append(candidate)
    return candidates, evaluator


def _candidate_close(left: RootCandidate, right: RootCandidate, policy: SearchPolicy) -> bool:
    tolerance = policy.dedup_atol_bar + policy.dedup_rtol * max(
        abs(left.omega_bar), abs(right.omega_bar)
    )
    return abs(left.omega_bar - right.omega_bar) <= tolerance


def _merge_candidates(
    candidates: Sequence[RootCandidate], policy: SearchPolicy
) -> tuple[list[RootCandidate], list[RootCandidate]]:
    accepted = sorted(
        (candidate for candidate in candidates if candidate.accepted),
        key=lambda item: (item.omega_bar, item.diagnostics.scaled_sigma_ratio),
    )
    canonical: list[RootCandidate] = []
    duplicates: list[RootCandidate] = []
    for candidate in accepted:
        match_index = next(
            (
                index
                for index, previous in enumerate(canonical)
                if _candidate_close(previous, candidate, policy)
            ),
            None,
        )
        if match_index is None:
            canonical.append(candidate)
            continue
        previous = canonical[match_index]
        best, duplicate = (
            (candidate, previous)
            if candidate.diagnostics.scaled_sigma_ratio < previous.diagnostics.scaled_sigma_ratio
            else (previous, candidate)
        )
        sources = tuple(sorted(set(previous.detection_sources + candidate.detection_sources)))
        merge_id = f"merge_{match_index + 1:04d}"
        canonical[match_index] = replace(
            best, detection_sources=sources, merge_group_id=merge_id, canonical=True
        )
        duplicates.append(
            replace(
                duplicate,
                accepted=False,
                rejection_reason="DUPLICATE_DETECTION_SAME_ROOT",
                canonical=False,
                merge_group_id=merge_id,
            )
        )
    rejected = [candidate for candidate in candidates if not candidate.accepted]
    rejected.extend(duplicates)
    canonical.sort(key=lambda item: item.omega_bar)
    return canonical, rejected


def _events_and_slots(
    canonical: Sequence[RootCandidate], policy: SearchPolicy
) -> tuple[tuple[RootEvent, ...], tuple[SpectrumSlot, ...]]:
    raw_events = [
        RootEvent(
            event_id=f"event_{index:04d}",
            omega_bar=candidate.omega_bar,
            omega=candidate.diagnostics.omega,
            multiplicity=candidate.diagnostics.detected_nullity,
            detected_nullity=candidate.diagnostics.detected_nullity,
            candidate=candidate,
            cluster_center_omega_bar=candidate.omega_bar,
        )
        for index, candidate in enumerate(canonical, start=1)
    ]
    groups: list[list[int]] = []
    for index, event in enumerate(raw_events):
        if not groups:
            groups.append([index])
            continue
        previous = raw_events[groups[-1][-1]]
        tolerance = policy.cluster_atol_bar + policy.cluster_rtol * max(
            abs(previous.omega_bar), abs(event.omega_bar)
        )
        if event.omega_bar - previous.omega_bar <= tolerance:
            groups[-1].append(index)
        else:
            groups.append([index])
    events = list(raw_events)
    cluster_counter = 0
    for indices in groups:
        multiplicity = sum(events[index].multiplicity for index in indices)
        if multiplicity <= 1:
            continue
        cluster_counter += 1
        cluster_id = f"cluster_{cluster_counter:04d}"
        centre = math.fsum(
            events[index].omega_bar * events[index].multiplicity for index in indices
        ) / multiplicity
        total_nullity = sum(events[index].detected_nullity for index in indices)
        exact = len(indices) == 1 and events[indices[0]].multiplicity > 1
        semantics = "EXACT_DEGENERATE_SUBSPACE" if exact else "NEAR_DEGENERATE_CLUSTER"
        for index in indices:
            events[index] = replace(
                events[index],
                cluster_id=cluster_id,
                cluster_semantics=semantics,
                cluster_multiplicity=multiplicity,
                cluster_total_nullity=total_nullity,
                cluster_center_omega_bar=centre,
            )
    slots: list[SpectrumSlot] = []
    for event in events:
        for repeated in range(1, event.multiplicity + 1):
            sorted_slot = len(slots) + 1
            role = "FIRST_12" if sorted_slot <= policy.requested_roots else "ROOT_13_GUARD"
            slots.append(SpectrumSlot(sorted_slot, role, repeated, event))
    retained = slots[: policy.required_slots]
    if retained and retained[-1].event.cluster_id:
        guard_cluster = retained[-1].event.cluster_id
        retained = [
            slot
            for slot in slots
            if slot.sorted_slot <= policy.required_slots or slot.event.cluster_id == guard_cluster
        ]
        retained = [
            replace(slot, role=("GUARD_CLUSTER_COMPLETION" if slot.sorted_slot > policy.required_slots else slot.role))
            for slot in retained
        ]
    return tuple(events), tuple(retained)


def _run_scan(
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
    *,
    case_id: str,
    builder_id: str,
    scan_id: str,
    points: int,
    phases: Sequence[float],
) -> ScanResult:
    raw_candidates, _evaluator = _scan_candidates(
        matrix_provider,
        frequency_scale,
        policy,
        case_id=case_id,
        builder_id=builder_id,
        scan_id=scan_id,
        points=points,
        phases=phases,
    )
    canonical, rejected = _merge_candidates(raw_candidates, policy)
    events, slots = _events_and_slots(canonical, policy)
    return ScanResult(
        scan_id=scan_id,
        candidates=tuple(canonical),
        rejected_candidates=tuple(rejected),
        events=events,
        slots=slots,
    )


def _inventory_digest_payload(slots: Sequence[SpectrumSlot]) -> list[dict[str, Any]]:
    return [
        {
            "sorted_slot": slot.sorted_slot,
            "omega_bar": format(slot.event.omega_bar, ".17g"),
            "multiplicity": slot.event.multiplicity,
            "nullity": slot.event.detected_nullity,
            "cluster_id": slot.event.cluster_id,
        }
        for slot in slots
    ]


def seed_free_root_inventory(
    matrix_provider: MatrixProvider,
    frequency_scale: float,
    policy: SearchPolicy,
    *,
    case_id: str,
    builder_id: str,
) -> RootInventory:
    """Build primary and independent verification inventories without seeds."""

    primary = _run_scan(
        matrix_provider, frequency_scale, policy,
        case_id=case_id, builder_id=builder_id, scan_id="primary",
        points=policy.primary_scan_points, phases=policy.primary_phases,
    )
    verification = _run_scan(
        matrix_provider, frequency_scale, policy,
        case_id=case_id, builder_id=builder_id, scan_id="verification",
        points=policy.verification_scan_points, phases=policy.verification_phases,
    )
    required = policy.required_slots
    primary_slots = primary.slots[:required]
    verification_slots = verification.slots[:required]
    agreement = len(primary_slots) >= required and len(verification_slots) >= required
    maximum_relative = 0.0
    if agreement:
        for left, right in zip(primary_slots, verification_slots, strict=True):
            relative = _relative_difference(left.event.omega_bar, right.event.omega_bar)
            left_clustered = bool(left.event.cluster_id)
            right_clustered = bool(right.event.cluster_id)
            cluster_relative = _relative_difference(
                left.event.cluster_center_omega_bar,
                right.event.cluster_center_omega_bar,
            )
            clustered = left_clustered or right_clustered
            comparison_relative = cluster_relative if clustered else relative
            maximum_relative = max(maximum_relative, comparison_relative)
            frequency_pass = (
                cluster_relative <= THRESHOLDS["cluster_center_relative"]
                if clustered
                else relative <= policy.inventory_match_relative
            )
            if (
                not frequency_pass
                or left.event.multiplicity != right.event.multiplicity
                or left.event.detected_nullity != right.event.detected_nullity
                or left_clustered != right_clustered
                or cluster_relative > THRESHOLDS["cluster_center_relative"]
                or left.event.cluster_multiplicity != right.event.cluster_multiplicity
                or left.event.cluster_total_nullity != right.event.cluster_total_nullity
            ):
                agreement = False
    guard_available = len(primary_slots) >= required
    guard_not_boundary = bool(
        guard_available
        and primary_slots[required - 1].event.omega_bar
        <= policy.omega_bar_max - policy.post_guard_tail_bar
    )
    guard_frequency = (
        primary_slots[required - 1].event.omega_bar if guard_available else policy.omega_bar_max
    )
    accepted_events = (*primary.events, *verification.events)

    def genuinely_unresolved(candidate: RootCandidate) -> bool:
        if candidate.interval_left_bar > guard_frequency:
            return False
        if any(
            (
                abs(candidate.omega_bar - event.omega_bar)
                <= policy.dedup_atol_bar
                + policy.dedup_rtol * max(abs(candidate.omega_bar), abs(event.omega_bar))
            )
            or (
                candidate.interval_left_bar - policy.root_xtol_bar
                <= event.omega_bar
                <= candidate.interval_right_bar + policy.root_xtol_bar
            )
            for event in accepted_events
        ):
            return False
        if candidate.rejection_reason == "DUPLICATE_DETECTION_SAME_ROOT":
            return False
        if candidate.rejection_reason in {
            "NONFINITE_MATRIX",
            "NONFINITE_DETERMINANT_INTERVAL",
            "BRENT_FAILURE",
            "MINIMIZER_EXCEPTION",
            "MINIMIZER_FAILURE",
            "NESTED_MINIMIZER_EXCEPTION",
        }:
            return True
        return candidate.diagnostics.scaled_sigma_ratio <= policy.sigma_prefilter

    unresolved = sum(
        genuinely_unresolved(candidate)
        for candidate in (*primary.rejected_candidates, *verification.rejected_candidates)
    )
    if not guard_available or len(verification_slots) < required or not agreement:
        status = "FAIL"
    elif guard_not_boundary and unresolved == 0:
        status = "PASS"
    else:
        status = "PARTIAL_PASS"
    payload = _inventory_digest_payload(primary.slots)
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest().upper()
    return RootInventory(
        case_id=case_id,
        builder_id=builder_id,
        frequency_scale=float(frequency_scale),
        policy=policy,
        primary=primary,
        verification=verification,
        slots=primary.slots,
        independent_agreement=agreement,
        maximum_primary_verification_relative=maximum_relative,
        unresolved_low_sigma_count=int(unresolved),
        guard_available=guard_available,
        guard_not_at_scan_boundary=guard_not_boundary,
        status=status,
        inventory_sha256=digest,
    )


def _coupled_module() -> Any:
    return importlib.import_module("scripts.lib.reddy_symmetric_coupled_beams")


def _selected_benchmarks() -> tuple[dict[str, Any], dict[str, Any]]:
    """Return the frozen with-RI ``a/h=20`` 0-degree and cross-ply cases."""

    contract, cases = source_audit.load_all_direct_cases()
    selected: dict[str, Any] = {}
    for laminate_id in ("0_deg", "cross_ply_0_90_s"):
        matches = [
            case
            for case in cases
            if case.laminate_id == laminate_id
            and case.a_over_h == 20
            and case.rotary_inertia
            and case.boundary_condition == "CC"
        ]
        if len(matches) != 1:
            raise RuntimeError(
                f"Expected one frozen CC/with-RI/a-h-20 benchmark for {laminate_id}; "
                f"found {len(matches)}"
            )
        selected[laminate_id] = matches[0]
    return contract, selected


def _properties_payload(case: Any) -> dict[str, Any]:
    properties = case.properties
    return {
        "laminate_id": case.laminate_id,
        "stack_degrees": [float(ply.angle_deg) for ply in case.laminate.plies],
        "total_thickness": float(case.laminate.thickness),
        "width": float(properties.width),
        "length_total": float(case.length),
        "a_over_h": int(case.a_over_h),
        "A_axial": float(properties.A),
        "D_bending": float(properties.D),
        "S_shear": float(properties.S),
        "mass_per_length": float(properties.m),
        "rotary_inertia_per_length": float(properties.J),
        "K": float(properties.K),
        "K_provenance": str(case.K_provenance),
        "frequency_scale_omega_bar_per_omega": float(case.frequency_scale),
    }


def build_model_manifest(policy: SearchPolicy | None = None) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build the threshold/model/search contract without running any gate."""

    active = policy or SearchPolicy()
    contract, selected = _selected_benchmarks()
    core_path = ROOT / "scripts" / "lib" / "reddy_symmetric_laminated_beam.py"
    geometry_path = ROOT / "scripts" / "lib" / "reddy_inplane_geometry.py"
    coupled_path = ROOT / "scripts" / "lib" / "reddy_symmetric_coupled_beams.py"
    pilot_path = Path(__file__).resolve()
    source_path = source_audit.SOURCE_CONTRACT_PATH
    total_length = float(selected["cross_ply_0_90_s"].length)
    manifest = {
        "schema_version": 1,
        "algorithm_version": ALGORITHM_VERSION,
        "stage": "RLB-1J_RLB-1A_RLB-1B_beta0_pilot",
        "diagnostic_scope": "diagnostic beta=0 coupled-joint baseline",
        "original_preimplementation_manifest_sha256": (
            ORIGINAL_PREIMPLEMENTATION_MANIFEST_SHA256
        ),
        "precomputation_threshold_freeze_sha256": _active_contract_hash(active),
        "active_threshold_search_contract_sha256": _active_contract_hash(active),
        "written_before_joint_random_or_spectral_calculation": True,
        "thresholds": dict(THRESHOLDS),
        "search_policy": asdict(active),
        "frequency_coordinate": {
            "name": "omega_bar",
            "definition": "omega_bar = omega*L_total^2*sqrt(I0/(E2*h^3))",
            "linear_map": True,
            "reference_length": "L_total=L1+L2",
        },
        "frozen_state_order": ["u", "w", "psi", "N", "Q", "M"],
        "frozen_equations": [
            "u_prime=N/A",
            "w_prime=Q/S-psi",
            "psi_prime=M/D",
            "N_prime=-m*omega^2*u",
            "Q_prime=-m*omega^2*w",
            "M_prime=Q-J*omega^2*psi",
        ],
        "coordinate_contract": {
            "k": "-E_Z",
            "local_x": "outer_clamp_to_joint",
            "canonical_helper": "scripts/lib/reddy_inplane_geometry.py",
        },
        "joint": {
            "type": "ideal_point_massless_rigid",
            "beta_rad_spectral": 0.0,
            "state_order": ["u", "w", "psi", "N", "Q", "M"],
            "clamp_basis": [
                [0, 0, 0], [0, 0, 0], [0, 0, 0],
                [1, 0, 0], [0, 1, 0], [0, 0, 1],
            ],
            "R0_beta0_diagonal": [-1, -1, 1, 1, 1, -1],
        },
        "benchmarks": {
            key: _properties_payload(case) for key, case in selected.items()
        },
        "homogeneous_splits": [
            {"split_id": "equal", "L1": 0.5 * total_length, "L2": 0.5 * total_length},
            {"split_id": "unequal_35_65", "L1": 0.35 * total_length, "L2": 0.65 * total_length},
        ],
        "stepped_case": {
            "arm_1": "0_deg",
            "arm_2": "cross_ply_0_90_s",
            "L1": 0.5 * total_length,
            "L2": 0.5 * total_length,
            "same_width_thickness_material_K": True,
        },
        "root_policy": {
            "primary_is_seed_free": True,
            "reference_roots_used_as_seeds": False,
            "reference_builders_run_after_primary_inventory_freeze": True,
            "detectors": ["determinant_sign_brackets", "normalized_sigma_minima"],
            "candidate_polishing": "local_offset_squared_sigma_plus_adjacent_float",
            "reference_detector_reconciliation": {
                "purpose": "reconcile_RLB0_determinant_and_SVD_approximations_only",
                "atol_bar": active.reference_detector_reconciliation_atol_bar,
                "rtol": active.reference_detector_reconciliation_rtol,
                "requires_unique_cross_method_match": True,
                "determinant_candidates_are_never_merged_with_each_other": True,
            },
            "multiplicity": "strict_scaled_SVD_nullity",
            "condition_number_is_not_a_root_gate": True,
        },
        "matrix_scaling": {
            "formula": "diag(row_factors) @ raw @ diag(column_factors)",
            "multipliers": "finite_strictly_positive",
            "row_norm_denominator_relative_floor": float(
                _coupled_module().MATRIX_SCALING_RELATIVE_FLOOR
            ),
            "column_policy": "downscale_only_never_amplify",
            "floor_purpose": "do_not_inflate_a_nearly_zero_root_row_or_column",
        },
        "source_contract": {
            "path": str(source_path.relative_to(ROOT)).replace("\\", "/"),
            "sha256": _sha256(source_path),
            "contract_version": contract.get("schema_version", ""),
        },
        "scientific_modules": {
            str(core_path.relative_to(ROOT)).replace("\\", "/"): _sha256(core_path),
            str(geometry_path.relative_to(ROOT)).replace("\\", "/"): _sha256(geometry_path),
            str(coupled_path.relative_to(ROOT)).replace("\\", "/"): _sha256(coupled_path),
            str(pilot_path.relative_to(ROOT)).replace("\\", "/"): _sha256(pilot_path),
        },
        "protected_artifact_policy": {
            "required_RLB0_RLB0A_RLB1G_inputs": "explicit_SHA256_before_after",
            "legacy_single_beam_results_tree": "recursive_SHA256_before_after",
            "Reddy_source_PDF": "explicit_SHA256_before_after",
        },
        "explicit_exclusions": [
            "beta_nonzero_spectra", "two_arm_Ritz", "FEM", "torsion", "damping",
            "complex_roots", "parameter_sweep", "figures", "branch_tracking",
        ],
    }
    return manifest, selected


def write_model_manifest(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    policy: SearchPolicy | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    manifest, selected = build_model_manifest(policy)
    _write_json(Path(output_dir) / "model_manifest.json", manifest)
    return manifest, selected


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--manifest-only",
        action="store_true",
        help="write the frozen model/threshold/search manifest and stop before all gates",
    )
    return parser.parse_args(argv)


def run_pilot(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    manifest_only: bool = False,
    policy: SearchPolicy | None = None,
) -> dict[str, Any]:
    """Write the manifest first, then execute the requested beta=0 pilot."""

    active = policy or SearchPolicy()
    model_manifest, selected = write_model_manifest(Path(output_dir), active)
    if manifest_only:
        run_manifest = {
            "algorithm_version": ALGORITHM_VERSION,
            "execution_mode": "manifest-only",
            "model_manifest_sha256": _sha256(Path(output_dir) / "model_manifest.json"),
            "original_preimplementation_manifest_sha256": (
                ORIGINAL_PREIMPLEMENTATION_MANIFEST_SHA256
            ),
            "precomputation_threshold_freeze_sha256": model_manifest[
                "precomputation_threshold_freeze_sha256"
            ],
            "formal_gates_executed": False,
            "spectral_calculations_executed": False,
            "statuses": {
                STATUS_JOINT: "NOT_RUN_MANIFEST_ONLY",
                STATUS_HOMOGENEOUS: "NOT_RUN_MANIFEST_ONLY",
                STATUS_STEPPED: "NOT_RUN_MANIFEST_ONLY",
                STATUS_INVENTORY: "NOT_RUN_MANIFEST_ONLY",
                STATUS_OVERALL: "NOT_RUN_MANIFEST_ONLY",
            },
        }
        _write_json(Path(output_dir) / "run_manifest.json", run_manifest)
        return run_manifest
    return _execute_full_pilot(Path(output_dir), active, model_manifest, selected)


@dataclass(frozen=True)
class CoupledCaseSpec:
    case_id: str
    category: str
    split_id: str
    order_id: str
    length_1: float
    properties_1: rlb.BeamProperties
    length_2: float
    properties_2: rlb.BeamProperties
    total_length: float
    frequency_scale: float


def _joint_matrix_checks(coupled: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    maximum_equality = 0.0
    minimum_rank = 6
    fixed_degrees = (-90.0, -30.0, 0.0, 30.0, 90.0)
    rng = np.random.default_rng(RANDOM_SEED_JOINT)
    angles = [math.radians(value) for value in fixed_degrees]
    angles.extend(float(value) for value in rng.uniform(-0.5 * math.pi, 0.5 * math.pi, 100))
    for index, beta_rad in enumerate(angles):
        invariant = np.asarray(coupled.joint_matrix_from_physical_maps(beta_rad), dtype=float)
        closed = np.asarray(coupled.joint_matrix_closed_form(beta_rad), dtype=float)
        difference = float(np.max(np.abs(invariant - closed)))
        singular = np.linalg.svd(invariant, compute_uv=False)
        rank = int(np.count_nonzero(singular / singular[0] > THRESHOLDS["rank_singular_relative"]))
        maximum_equality = max(maximum_equality, difference)
        minimum_rank = min(minimum_rank, rank)
        rows.append(
            {
                "check": "invariant_vs_closed_form",
                "case_index": index,
                "beta_deg": math.degrees(beta_rad),
                "angle_role": "fixed" if index < len(fixed_degrees) else "fixed_seed_random",
                "max_abs_difference": difference,
                "matrix_rank": rank,
                "minimum_relative_singular_value": float(singular[-1] / singular[0]),
                "equality_tolerance": THRESHOLDS["coordinate_J_equality"],
                "rank_relative_tolerance": THRESHOLDS["rank_singular_relative"],
                "status": "PASS" if difference <= THRESHOLDS["coordinate_J_equality"] and rank == 6 else "FAIL",
            }
        )
    r0_maximum = 0.0
    matrix_zero = np.asarray(coupled.joint_matrix(0.0), dtype=float)
    for index in range(100):
        state_1 = rng.normal(size=6)
        joined = np.concatenate((state_1, np.asarray(coupled.R0_BETA0) @ state_1))
        residual = float(np.linalg.norm(matrix_zero @ joined))
        scale = max(float(np.linalg.norm(matrix_zero) * np.linalg.norm(joined)), np.finfo(float).tiny)
        normalized = residual / scale
        r0_maximum = max(r0_maximum, normalized)
        rows.append(
            {
                "check": "beta0_R0_identity",
                "case_index": index,
                "beta_deg": 0.0,
                "angle_role": "fixed_seed_random_state",
                "max_abs_difference": residual,
                "matrix_rank": 6,
                "minimum_relative_singular_value": "",
                "normalized_residual": normalized,
                "equality_tolerance": THRESHOLDS["coordinate_J_equality"],
                "rank_relative_tolerance": THRESHOLDS["rank_singular_relative"],
                "status": "PASS" if normalized <= THRESHOLDS["coordinate_J_equality"] else "FAIL",
            }
        )
    beta90_expected = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        ],
        dtype=float,
    )
    beta90_difference = float(
        np.max(np.abs(np.asarray(coupled.joint_matrix(0.5 * math.pi)) - beta90_expected))
    )
    rows.append(
        {
            "check": "beta90_scalar_limit",
            "case_index": 0,
            "beta_deg": 90.0,
            "angle_role": "fixed_limit",
            "max_abs_difference": beta90_difference,
            "matrix_rank": int(np.linalg.matrix_rank(beta90_expected)),
            "minimum_relative_singular_value": "",
            "equality_tolerance": THRESHOLDS["coordinate_J_equality"],
            "rank_relative_tolerance": THRESHOLDS["rank_singular_relative"],
            "status": "PASS" if beta90_difference <= THRESHOLDS["coordinate_J_equality"] else "FAIL",
        }
    )
    clamp_difference = float(
        np.max(
            np.abs(
                np.asarray(coupled.CLAMP_BASIS)
                - np.vstack((np.zeros((3, 3)), np.eye(3)))
            )
        )
    )
    rows.append(
        {
            "check": "outer_clamp_basis",
            "case_index": 0,
            "beta_deg": "",
            "angle_role": "exact_matrix",
            "max_abs_difference": clamp_difference,
            "matrix_rank": int(np.linalg.matrix_rank(coupled.CLAMP_BASIS)),
            "minimum_relative_singular_value": 1.0,
            "equality_tolerance": 0.0,
            "rank_relative_tolerance": THRESHOLDS["rank_singular_relative"],
            "status": "PASS" if clamp_difference == 0.0 and np.linalg.matrix_rank(coupled.CLAMP_BASIS) == 3 else "FAIL",
        }
    )
    summary = {
        "invariant_closed_form_case_count": len(angles),
        "beta0_R0_random_state_count": 100,
        "total_check_row_count": len(rows),
        "maximum_invariant_closed_form_absolute_difference": maximum_equality,
        "minimum_joint_rank": minimum_rank,
        "maximum_beta0_R0_normalized_residual": r0_maximum,
        "beta90_maximum_absolute_difference": beta90_difference,
        "clamp_basis_maximum_absolute_difference": clamp_difference,
        "pass": all(row["status"] == "PASS" for row in rows),
    }
    return rows, summary


def _virtual_work_checks(coupled: Any) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    from scripts.lib import reddy_inplane_geometry as geometry

    rng = np.random.default_rng(RANDOM_SEED_VIRTUAL_WORK)
    angles = [0.0, math.radians(30.0), 0.5 * math.pi]
    angles.extend(float(value) for value in rng.uniform(-0.5 * math.pi, 0.5 * math.pi, 1000))
    rows: list[dict[str, Any]] = []
    max_absolute = max_normalized = max_pairing = 0.0
    for index, beta_rad in enumerate(angles):
        bases = geometry.reddy_inplane_geometry(math.degrees(beta_rad))
        delta_global = rng.normal(size=2)
        delta_vector = delta_global[0] * geometry.E_X + delta_global[1] * geometry.E_Y
        delta_psi = float(rng.normal())
        virtual_1 = np.array(
            [np.dot(delta_vector, bases.t1), np.dot(delta_vector, bases.n1), delta_psi]
        )
        virtual_2 = np.array(
            [np.dot(delta_vector, bases.t2), np.dot(delta_vector, bases.n2), delta_psi]
        )
        force_1_global = rng.normal(size=2)
        force_1_vector = force_1_global[0] * geometry.E_X + force_1_global[1] * geometry.E_Y
        force_2_vector = -force_1_vector
        moment_1 = float(rng.normal())
        resultants_1 = np.array(
            [np.dot(force_1_vector, bases.t1), np.dot(force_1_vector, bases.n1), moment_1]
        )
        resultants_2 = np.array(
            [np.dot(force_2_vector, bases.t2), np.dot(force_2_vector, bases.n2), -moment_1]
        )
        check = coupled.joint_virtual_work_check(
            beta_rad, resultants_1, resultants_2, virtual_1, virtual_2
        )
        max_absolute = max(max_absolute, check.absolute_residual)
        max_normalized = max(max_normalized, check.normalized_residual)
        max_pairing = max(max_pairing, check.pairing_absolute_difference)
        rows.append(
            {
                "case_index": index,
                "beta_deg": math.degrees(beta_rad),
                "angle_role": "fixed" if index < 3 else "fixed_seed_random",
                "local_total": check.local_total,
                "global_total": check.global_total,
                "pairing_absolute_difference": check.pairing_absolute_difference,
                "absolute_residual": check.absolute_residual,
                "scale": check.scale,
                "normalized_residual": check.normalized_residual,
                "tolerance": THRESHOLDS["virtual_work_normalized_residual"],
                "status": "PASS" if check.normalized_residual <= THRESHOLDS["virtual_work_normalized_residual"] else "FAIL",
            }
        )
    return rows, {
        "case_count": len(rows),
        "random_seed": RANDOM_SEED_VIRTUAL_WORK,
        "maximum_absolute_residual": max_absolute,
        "maximum_normalized_residual": max_normalized,
        "maximum_pairing_absolute_difference": max_pairing,
        "pass": all(row["status"] == "PASS" for row in rows),
    }


def _case_specs(selected: Mapping[str, Any]) -> tuple[list[CoupledCaseSpec], list[CoupledCaseSpec]]:
    zero = selected["0_deg"]
    cross = selected["cross_ply_0_90_s"]
    total = float(cross.length)
    scale = float(cross.frequency_scale)
    homogeneous: list[CoupledCaseSpec] = []
    for laminate_id, case in (("0_deg", zero), ("cross_ply_0_90_s", cross)):
        for split_id, fraction in (("equal", 0.5), ("unequal_35_65", 0.35)):
            length_1 = fraction * total
            homogeneous.append(
                CoupledCaseSpec(
                    case_id=f"homogeneous__{laminate_id}__{split_id}",
                    category="homogeneous",
                    split_id=split_id,
                    order_id=f"{laminate_id}|{laminate_id}",
                    length_1=length_1,
                    properties_1=case.properties,
                    length_2=total - length_1,
                    properties_2=case.properties,
                    total_length=total,
                    frequency_scale=float(case.frequency_scale),
                )
            )
    stepped = [
        CoupledCaseSpec(
            case_id="stepped__0_deg__cross_ply",
            category="stepped", split_id="equal", order_id="0_deg|cross_ply_0_90_s",
            length_1=0.5 * total, properties_1=zero.properties,
            length_2=0.5 * total, properties_2=cross.properties,
            total_length=total, frequency_scale=scale,
        ),
        CoupledCaseSpec(
            case_id="stepped__cross_ply__0_deg",
            category="stepped", split_id="equal", order_id="cross_ply_0_90_s|0_deg",
            length_1=0.5 * total, properties_1=cross.properties,
            length_2=0.5 * total, properties_2=zero.properties,
            total_length=total, frequency_scale=scale,
        ),
    ]
    return homogeneous, stepped


def _coupled_provider(coupled: Any, spec: CoupledCaseSpec) -> MatrixProvider:
    joint_beta0 = np.asarray(coupled.joint_matrix(0.0), dtype=float)

    def provider(omega: float) -> FloatArray:
        clamp_to_joint = coupled.coupled_clamp_to_joint_map(
            omega,
            spec.length_1,
            spec.properties_1,
            spec.length_2,
            spec.properties_2,
        )
        return np.asarray(joint_beta0 @ clamp_to_joint, dtype=float)
    return provider


def _direct_homogeneous_provider(coupled: Any, spec: CoupledCaseSpec) -> MatrixProvider:
    def provider(omega: float) -> FloatArray:
        return np.asarray(
            coupled.direct_fixed_fixed_boundary_matrix(
                omega, spec.total_length, spec.properties_1
            ),
            dtype=float,
        )
    return provider


def _direct_stepped_provider(coupled: Any, spec: CoupledCaseSpec) -> MatrixProvider:
    def provider(omega: float) -> FloatArray:
        return np.asarray(
            coupled.direct_stepped_boundary_matrix(
                omega, spec.length_1, spec.properties_1,
                spec.length_2, spec.properties_2,
            ),
            dtype=float,
        )
    return provider


def _inventory_groups(inventory: RootInventory) -> list[dict[str, Any]]:
    groups: list[dict[str, Any]] = []
    by_key: dict[str, dict[str, Any]] = {}
    for slot in inventory.slots:
        key = slot.event.cluster_id or slot.event.event_id
        if key not in by_key:
            record = {
                "group_id": key,
                "first_slot": slot.sorted_slot,
                "last_slot": slot.sorted_slot,
                "semantics": slot.event.cluster_semantics,
                "center_omega_bar": slot.event.cluster_center_omega_bar,
                "multiplicity": slot.event.cluster_multiplicity,
                "total_nullity": slot.event.cluster_total_nullity,
            }
            by_key[key] = record
            groups.append(record)
        else:
            by_key[key]["last_slot"] = slot.sorted_slot
    return groups


def _compare_inventories(
    coupled_inventory: RootInventory,
    reference_inventory: RootInventory,
    *,
    case_id: str,
    comparison_kind: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    left_groups = _inventory_groups(coupled_inventory)
    right_groups = _inventory_groups(reference_inventory)
    rows: list[dict[str, Any]] = []
    maximum_relative = 0.0
    passed = len(left_groups) == len(right_groups)
    for index in range(max(len(left_groups), len(right_groups))):
        left = left_groups[index] if index < len(left_groups) else None
        right = right_groups[index] if index < len(right_groups) else None
        if left is None or right is None:
            rows.append(
                {
                    "case_id": case_id,
                    "comparison_kind": comparison_kind,
                    "comparison_index": index + 1,
                    "comparison_unit": "MISSING",
                    "coupled_first_slot": "" if left is None else left["first_slot"],
                    "reference_first_slot": "" if right is None else right["first_slot"],
                    "status": "FAIL",
                }
            )
            passed = False
            continue
        clustered = left["multiplicity"] > 1 or right["multiplicity"] > 1
        threshold = (
            THRESHOLDS["cluster_center_relative"]
            if clustered
            else THRESHOLDS["beta0_isolated_spectral_relative"]
        )
        relative = _relative_difference(left["center_omega_bar"], right["center_omega_bar"])
        maximum_relative = max(maximum_relative, relative)
        row_pass = bool(
            relative <= threshold
            and left["multiplicity"] == right["multiplicity"]
            and left["total_nullity"] == right["total_nullity"]
            and left["first_slot"] == right["first_slot"]
            and left["last_slot"] == right["last_slot"]
        )
        passed = passed and row_pass
        rows.append(
            {
                "case_id": case_id,
                "comparison_kind": comparison_kind,
                "comparison_index": index + 1,
                "comparison_unit": "CLUSTER" if clustered else "ISOLATED",
                "coupled_first_slot": left["first_slot"],
                "coupled_last_slot": left["last_slot"],
                "reference_first_slot": right["first_slot"],
                "reference_last_slot": right["last_slot"],
                "coupled_semantics": left["semantics"],
                "reference_semantics": right["semantics"],
                "coupled_omega_bar_or_center": left["center_omega_bar"],
                "reference_omega_bar_or_center": right["center_omega_bar"],
                "relative_difference": relative,
                "tolerance": threshold,
                "coupled_multiplicity": left["multiplicity"],
                "reference_multiplicity": right["multiplicity"],
                "coupled_total_nullity": left["total_nullity"],
                "reference_total_nullity": right["total_nullity"],
                "status": "PASS" if row_pass else "FAIL",
            }
        )
    return rows, {
        "comparison_kind": comparison_kind,
        "maximum_relative_difference": maximum_relative,
        "coupled_group_count": len(left_groups),
        "reference_group_count": len(right_groups),
        "coupled_slots_through_guard": min(13, len(coupled_inventory.slots)),
        "reference_slots_through_guard": min(13, len(reference_inventory.slots)),
        "pass": passed,
    }


def _reconcile_reference_detector_detections(
    detections: Sequence[Any],
    *,
    frequency_scale: float,
    policy: SearchPolicy,
) -> tuple[tuple[tuple[Any, ...], ...], dict[str, int]]:
    """Reconcile only unique cross-method approximations of one root.

    Determinant detections are canonical and are never merged with each other.
    A non-determinant detection is attached only when exactly one still
    unmatched determinant detection lies inside the frozen reconciliation
    window.  Ambiguous or unmatched detections remain separate.
    """

    scale = _finite_positive(frequency_scale, "frequency_scale")
    ordered = sorted(detections, key=lambda item: float(item.omega))
    groups: list[list[Any]] = [
        [item]
        for item in ordered
        if str(item.detection).startswith("determinant")
    ]
    supplemental = [
        item
        for item in ordered
        if not str(item.detection).startswith("determinant")
    ]
    reconciled = ambiguous = unmatched = 0
    absolute = policy.reference_detector_reconciliation_atol_bar / scale
    relative = policy.reference_detector_reconciliation_rtol
    for item in supplemental:
        matches: list[int] = []
        for index, group in enumerate(groups):
            # One-to-one, cross-method matching only.
            if len(group) != 1 or not str(group[0].detection).startswith("determinant"):
                continue
            determinant = group[0]
            if abs(float(item.omega) - float(determinant.omega)) <= (
                absolute
                + relative
                * max(abs(float(item.omega)), abs(float(determinant.omega)))
            ):
                matches.append(index)
        if len(matches) == 1:
            groups[matches[0]].append(item)
            reconciled += 1
        else:
            groups.append([item])
            if matches:
                ambiguous += 1
            else:
                unmatched += 1
    groups.sort(key=lambda group: min(float(item.omega) for item in group))
    return tuple(tuple(group) for group in groups), {
        "raw_detection_count": len(ordered),
        "reconciled_cross_method_count": reconciled,
        "ambiguous_cross_method_count": ambiguous,
        "unmatched_supplemental_count": unmatched,
        "reconciled_group_count": len(groups),
    }


def _union_reference(
    case: Any,
    policy: SearchPolicy,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Build the post-freeze axial-FF plus bending-CC family union."""

    scale = float(case.frequency_scale)
    axial = rlb.exact_axial_modes(case.properties, case.length, "FF", n_modes=24)
    bending_detections = rlb.find_bending_roots(
        case.properties,
        case.length,
        "CC",
        omega_max=policy.omega_bar_max / scale,
        n_roots=None,
        omega_min=0.0,
        scan_points=policy.verification_scan_points,
        sigma_ratio_tolerance=policy.sigma_ratio_tolerance,
        root_xtol=max(policy.root_xtol_bar / scale, 1.0e-15),
        dedup_rtol=1.0e-12,
    )
    # RLB-0 returns determinant-bracket and SVD-minimum detections.  Retain
    # both methods as diagnostics, but reconcile only a unique cross-method
    # pair inside the separately frozen reference-detector window.
    bending_groups, reconciliation = _reconcile_reference_detector_detections(
        bending_detections,
        frequency_scale=scale,
        policy=policy,
    )
    bending = tuple(
        min(
            group,
            key=lambda item: (
                0 if item.detection.startswith("determinant") else 1,
                item.boundary_residual,
                item.sigma_ratio,
            ),
        )
        for group in bending_groups
    )[:24]
    clusters = rlb.union_subsystem_spectra(
        [item.omega for item in axial],
        [item.omega for item in bending],
        atol=policy.cluster_atol_bar / scale,
        rtol=policy.cluster_rtol,
    )
    rows: list[dict[str, Any]] = []
    slot = 0
    for cluster_index, cluster in enumerate(clusters, start=1):
        if slot >= policy.required_slots and cluster.multiplicity == 1:
            break
        semantics = (
            "ISOLATED"
            if cluster.multiplicity == 1
            else (
                "EXACT_DEGENERATE_SUBSPACE"
                if cluster.exact_degeneracy
                else "NEAR_DEGENERATE_CLUSTER"
            )
        )
        family_content = "+".join(sorted(member.subsystem for member in cluster.members))
        for member_index, member in enumerate(cluster.members, start=1):
            slot += 1
            rows.append(
                {
                    "sorted_slot": slot,
                    "role": (
                        "FIRST_12"
                        if slot <= 12
                        else ("ROOT_13_GUARD" if slot == 13 else "GUARD_CLUSTER_COMPLETION")
                    ),
                    "cluster_index": cluster_index,
                    "cluster_semantics": semantics,
                    "cluster_multiplicity": cluster.multiplicity,
                    "cluster_total_nullity": cluster.multiplicity,
                    "cluster_center_omega_bar": cluster.representative_omega * scale,
                    "cluster_member_slot": member_index,
                    "omega": member.omega,
                    "omega_bar": member.omega * scale,
                    "family": member.subsystem,
                    "family_index": member.subsystem_index,
                    "family_content": family_content,
                }
            )
        if slot >= policy.required_slots:
            break
    return rows, {
        "axial_root_count": len(axial),
        "bending_root_count": len(bending),
        "reference_detector_reconciliation": reconciliation,
        "union_slots_retained": len(rows),
        "guard_available": len(rows) >= policy.required_slots,
    }


def _compare_inventory_to_union(
    inventory: RootInventory,
    union_rows: Sequence[Mapping[str, Any]],
    *,
    case_id: str,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    direct_groups = _inventory_groups(inventory)
    union_groups: list[dict[str, Any]] = []
    for row in union_rows:
        cluster_index = int(row["cluster_index"])
        if union_groups and union_groups[-1]["cluster_index"] == cluster_index:
            union_groups[-1]["last_slot"] = int(row["sorted_slot"])
            continue
        union_groups.append(
            {
                "cluster_index": cluster_index,
                "first_slot": int(row["sorted_slot"]),
                "last_slot": int(row["sorted_slot"]),
                "semantics": str(row["cluster_semantics"]),
                "center_omega_bar": float(row["cluster_center_omega_bar"]),
                "multiplicity": int(row["cluster_multiplicity"]),
                "total_nullity": int(row["cluster_total_nullity"]),
                "family_content": str(row["family_content"]),
            }
        )
    rows: list[dict[str, Any]] = []
    maximum = 0.0
    passed = (
        len(inventory.slots) >= 13
        and len(union_rows) >= 13
        and len(direct_groups) == len(union_groups)
    )
    for index in range(max(len(direct_groups), len(union_groups))):
        left = direct_groups[index] if index < len(direct_groups) else None
        right = union_groups[index] if index < len(union_groups) else None
        if left is None or right is None:
            rows.append(
                {
                    "case_id": case_id,
                    "comparison_kind": "direct_fixed_fixed_vs_axial_bending_union",
                    "comparison_index": index + 1,
                    "comparison_unit": "MISSING",
                    "status": "FAIL",
                }
            )
            passed = False
            continue
        clustered = left["multiplicity"] > 1 or right["multiplicity"] > 1
        threshold = (
            THRESHOLDS["cluster_center_relative"]
            if clustered
            else THRESHOLDS["beta0_isolated_spectral_relative"]
        )
        relative = _relative_difference(left["center_omega_bar"], right["center_omega_bar"])
        maximum = max(maximum, relative)
        row_pass = bool(
            relative <= threshold
            and left["multiplicity"] == right["multiplicity"]
            and left["total_nullity"] == right["total_nullity"]
            and left["first_slot"] == right["first_slot"]
            and left["last_slot"] == right["last_slot"]
        )
        passed = passed and row_pass
        rows.append(
            {
                "case_id": case_id,
                "comparison_kind": "direct_fixed_fixed_vs_axial_bending_union",
                "comparison_index": index + 1,
                "comparison_unit": "CLUSTER" if clustered else "ISOLATED",
                "direct_first_slot": left["first_slot"],
                "direct_last_slot": left["last_slot"],
                "union_first_slot": right["first_slot"],
                "union_last_slot": right["last_slot"],
                "direct_semantics": left["semantics"],
                "union_semantics": right["semantics"],
                "direct_omega_bar_or_center": left["center_omega_bar"],
                "union_omega_bar_or_center": right["center_omega_bar"],
                "relative_difference": relative,
                "tolerance": threshold,
                "direct_multiplicity": left["multiplicity"],
                "union_multiplicity": right["multiplicity"],
                "direct_total_nullity": left["total_nullity"],
                "union_total_nullity": right["total_nullity"],
                "family_content": right["family_content"],
                "status": "PASS" if row_pass else "FAIL",
            }
        )
    return rows, {
        "maximum_relative_difference": maximum,
        "root_count_through_guard": min(len(inventory.slots), len(union_rows), 13),
        "pass": passed,
    }


def _candidate_row(candidate: RootCandidate) -> dict[str, Any]:
    diagnostic = candidate.diagnostics
    return {
        "case_id": candidate.case_id,
        "builder_id": candidate.builder_id,
        "scan_id": candidate.scan_id,
        "omega_bar": candidate.omega_bar,
        "omega": diagnostic.omega,
        "detection_sources": candidate.detection_sources,
        "interval_left_bar": candidate.interval_left_bar,
        "interval_right_bar": candidate.interval_right_bar,
        "interior_minimum": candidate.interior_minimum,
        "raw_determinant": diagnostic.raw_determinant,
        "raw_det_sign": diagnostic.raw_det_sign,
        "raw_logabsdet": diagnostic.raw_logabsdet,
        "scaled_determinant": diagnostic.scaled_determinant,
        "scaled_det_sign": diagnostic.scaled_det_sign,
        "scaled_logabsdet": diagnostic.scaled_logabsdet,
        "raw_sigma_min": diagnostic.raw_sigma_min,
        "raw_sigma_max": diagnostic.raw_sigma_max,
        "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
        "scaled_sigma_min": diagnostic.scaled_sigma_min,
        "scaled_sigma_max": diagnostic.scaled_sigma_max,
        "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
        "raw_condition_number": diagnostic.raw_condition_number,
        "scaled_condition_number": diagnostic.scaled_condition_number,
        "detected_nullity": diagnostic.detected_nullity,
        "root_gate_nullity": diagnostic.root_gate_nullity,
        "scaled_null_residual": diagnostic.scaled_null_residual,
        "raw_boundary_null_residual": diagnostic.raw_boundary_null_residual,
        "row_factor_min": float(np.min(diagnostic.row_factors)),
        "row_factor_max": float(np.max(diagnostic.row_factors)),
        "column_factor_min": float(np.min(diagnostic.column_factors)),
        "column_factor_max": float(np.max(diagnostic.column_factors)),
        "merge_group_id": candidate.merge_group_id,
        "canonical": candidate.canonical,
        "accepted": candidate.accepted,
        "rejection_reason": candidate.rejection_reason,
    }


def _slot_row(spec: CoupledCaseSpec, inventory: RootInventory, slot: SpectrumSlot) -> dict[str, Any]:
    event = slot.event
    diagnostic = event.candidate.diagnostics
    return {
        "case_id": spec.case_id,
        "category": spec.category,
        "split_id": spec.split_id,
        "order_id": spec.order_id,
        "builder_id": inventory.builder_id,
        "sorted_slot": slot.sorted_slot,
        "role": slot.role,
        "root_event_id": event.event_id,
        "repeated_root_slot": slot.repeated_root_slot,
        "cluster_id": event.cluster_id,
        "cluster_semantics": event.cluster_semantics,
        "event_multiplicity": event.multiplicity,
        "event_nullity": event.detected_nullity,
        "cluster_multiplicity": event.cluster_multiplicity,
        "cluster_total_nullity": event.cluster_total_nullity,
        "cluster_center_omega_bar": event.cluster_center_omega_bar,
        "omega": event.omega,
        "omega_bar": event.omega_bar,
        "detection_sources": event.candidate.detection_sources,
        "raw_determinant": diagnostic.raw_determinant,
        "scaled_determinant": diagnostic.scaled_determinant,
        "raw_sigma_min": diagnostic.raw_sigma_min,
        "raw_sigma_max": diagnostic.raw_sigma_max,
        "raw_sigma_ratio": diagnostic.raw_sigma_ratio,
        "scaled_sigma_min": diagnostic.scaled_sigma_min,
        "scaled_sigma_max": diagnostic.scaled_sigma_max,
        "scaled_sigma_ratio": diagnostic.scaled_sigma_ratio,
        "raw_condition_number": diagnostic.raw_condition_number,
        "scaled_condition_number": diagnostic.scaled_condition_number,
        "scaled_null_residual": diagnostic.scaled_null_residual,
        "raw_boundary_null_residual": diagnostic.raw_boundary_null_residual,
        "inventory_status": inventory.status,
        "accepted": inventory.status == "PASS",
    }


def _representative_residual(
    coupled: Any,
    spec: CoupledCaseSpec,
    inventory: RootInventory,
    slot: SpectrumSlot,
    *,
    family: str,
) -> dict[str, Any]:
    diagnostic = slot.event.candidate.diagnostics
    reactions = np.asarray(diagnostic.physical_null_vector, dtype=float)
    reaction_1, reaction_2 = reactions[:3], reactions[3:]
    omega = slot.event.omega
    map_1 = coupled.arm_clamp_map(omega, spec.length_1, spec.properties_1)
    map_2 = coupled.arm_clamp_map(omega, spec.length_2, spec.properties_2)
    joint_1 = np.asarray(map_1 @ reaction_1, dtype=float)
    joint_2 = np.asarray(map_2 @ reaction_2, dtype=float)
    physical = coupled.physical_joint_residuals(0.0, joint_1, joint_2)

    displacement_absolute = float(np.linalg.norm(physical.displacement))
    rotation_absolute = float(np.linalg.norm(physical.rotation))
    force_absolute = float(np.linalg.norm(physical.force))
    moment_absolute = float(np.linalg.norm(physical.moment))
    reference_length = spec.total_length
    displacement_scale = max(
        math.hypot(joint_1[0], joint_1[1]),
        math.hypot(joint_2[0], joint_2[1]),
        reference_length * abs(joint_1[2]),
        reference_length * abs(joint_2[2]),
        np.finfo(float).tiny,
    )
    rotation_scale = displacement_scale / reference_length
    force_scale = max(
        math.hypot(joint_1[3], joint_1[4]),
        math.hypot(joint_2[3], joint_2[4]),
        math.hypot(reaction_1[0], reaction_1[1]),
        math.hypot(reaction_2[0], reaction_2[1]),
        abs(joint_1[5]) / reference_length,
        abs(joint_2[5]) / reference_length,
        abs(reaction_1[2]) / reference_length,
        abs(reaction_2[2]) / reference_length,
        np.finfo(float).tiny,
    )
    moment_scale = force_scale * reference_length
    compatibility = max(
        displacement_absolute / displacement_scale,
        rotation_absolute / rotation_scale,
    )
    equilibrium = max(force_absolute / force_scale, moment_absolute / moment_scale)

    initial_1 = np.asarray(coupled.CLAMP_BASIS @ reaction_1, dtype=float)
    initial_2 = np.asarray(coupled.CLAMP_BASIS @ reaction_2, dtype=float)
    outer_absolute = max(
        float(np.linalg.norm(initial_1[:3])), float(np.linalg.norm(initial_2[:3]))
    )
    outer_scale = max(
        float(np.linalg.norm(initial_1)), float(np.linalg.norm(initial_2)), np.finfo(float).tiny
    )
    outer_normalized = outer_absolute / outer_scale

    coordinates_1 = np.linspace(0.0, spec.length_1, 1601)
    coordinates_2 = np.linspace(0.0, spec.length_2, 1601)
    states_1 = np.vstack(
        [
            coupled.arm_transfer_matrix(omega, float(x), spec.properties_1) @ initial_1
            for x in coordinates_1
        ]
    )
    states_2 = np.vstack(
        [
            coupled.arm_transfer_matrix(omega, float(x), spec.properties_2) @ initial_2
            for x in coordinates_2
        ]
    )
    energy_1 = rlb.combined_modal_energies(
        coordinates_1, states_1, spec.properties_1, omega
    )
    energy_2 = rlb.combined_modal_energies(
        coordinates_2, states_2, spec.properties_2, omega
    )
    total_stiffness = energy_1.stiffness_integral + energy_2.stiffness_integral
    total_inertia = energy_1.inertia_integral + energy_2.inertia_integral
    energy_error = abs(total_stiffness - total_inertia) / max(
        abs(total_stiffness), abs(total_inertia), np.finfo(float).tiny
    )
    passed = bool(
        compatibility <= THRESHOLDS["joint_compatibility_residual"]
        and equilibrium <= THRESHOLDS["joint_equilibrium_residual"]
        and outer_normalized <= THRESHOLDS["boundary_residual"]
        and diagnostic.raw_boundary_null_residual <= THRESHOLDS["boundary_residual"]
        and energy_error <= THRESHOLDS["energy_identity"]
    )
    return {
        "case_id": spec.case_id,
        "family": family,
        "sorted_slot": slot.sorted_slot,
        "omega": omega,
        "omega_bar": slot.event.omega_bar,
        "cluster_semantics": slot.event.cluster_semantics,
        "reaction_coefficients": reactions,
        "joint_state_1": joint_1,
        "joint_state_2": joint_2,
        "displacement_residual_absolute": displacement_absolute,
        "rotation_residual_absolute": rotation_absolute,
        "force_residual_absolute": force_absolute,
        "moment_residual_absolute": moment_absolute,
        "joint_compatibility_normalized": compatibility,
        "joint_equilibrium_normalized": equilibrium,
        "outer_clamp_residual_absolute": outer_absolute,
        "outer_clamp_residual_normalized": outer_normalized,
        "boundary_null_residual": diagnostic.raw_boundary_null_residual,
        "total_modal_mass": energy_1.modal_mass + energy_2.modal_mass,
        "total_stiffness_integral": total_stiffness,
        "total_inertia_integral": total_inertia,
        "energy_identity_relative": energy_error,
        "compatibility_tolerance": THRESHOLDS["joint_compatibility_residual"],
        "equilibrium_tolerance": THRESHOLDS["joint_equilibrium_residual"],
        "boundary_tolerance": THRESHOLDS["boundary_residual"],
        "energy_tolerance": THRESHOLDS["energy_identity"],
        "status": "PASS" if passed else "FAIL",
    }


def _matrix_records(inventories: Iterable[RootInventory]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for inventory in inventories:
        seen: set[str] = set()
        for slot in inventory.slots:
            event = slot.event
            if event.event_id in seen:
                continue
            seen.add(event.event_id)
            diagnostic = event.candidate.diagnostics
            rows.append(
                {
                    "case_id": inventory.case_id,
                    "builder_id": inventory.builder_id,
                    "root_event_id": event.event_id,
                    "omega_bar": event.omega_bar,
                    "omega": event.omega,
                    "raw_matrix": diagnostic.raw_matrix,
                    "positively_equilibrated_matrix": diagnostic.scaled_matrix,
                    "row_factors": diagnostic.row_factors,
                    "column_factors": diagnostic.column_factors,
                    "raw_singular_values": diagnostic.raw_singular_values,
                    "scaled_singular_values": diagnostic.scaled_singular_values,
                    "physical_right_null_vector": diagnostic.physical_null_vector,
                }
            )
    return rows


def _active_contract_hash(policy: SearchPolicy) -> str:
    payload = {"thresholds": THRESHOLDS, "search_policy": asdict(policy)}
    return hashlib.sha256(
        json.dumps(_json_value(payload), sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest().upper()


def _preservation_hashes() -> dict[str, str]:
    required_paths = {
        ROOT / "scripts" / "lib" / "reddy_symmetric_laminated_beam.py",
        ROOT / "scripts" / "lib" / "reddy_inplane_geometry.py",
        ROOT / "scripts" / "analysis" / "laminated_beams" / "validate_reddy_symmetric_single_beam.py",
        ROOT / "scripts" / "analysis" / "laminated_beams" / "audit_reddy_table_4_3_3_discrepancies.py",
        ROOT / "docs" / "laminated_beams" / "reddy_inplane_coordinate_contract.md",
        ROOT / "docs" / "laminated_beams" / "reddy_ch4_source_contract.md",
        ROOT / "docs" / "laminated_beams" / "reddy_symmetric_single_beam_validation.md",
        ROOT / "docs" / "laminated_beams" / "reddy_table_4_3_3_discrepancy_audit.md",
        ROOT / "tests" / "test_reddy_inplane_geometry.py",
        ROOT / "tests" / "test_reddy_symmetric_laminated_beam.py",
        ROOT / "tests" / "test_reddy_table_4_3_3_discrepancy_audit.py",
        ROOT / "tests" / "data" / "reddy_ch4_table_4_3_3.json",
        ROOT
        / "docs"
        / "literature"
        / "pdf"
        / "EB__Mechanics_of_Laminated_Composite_Plates_and_Shells_-JN_Reddy.pdf",
    }
    legacy_results = (
        ROOT
        / "results"
        / "laminated_beams"
        / "reddy_symmetric_single_beam"
    )
    required_paths.update(
        path for path in legacy_results.rglob("*") if path.is_file()
    )
    return {
        str(path.relative_to(ROOT)).replace("\\", "/"): _sha256(path)
        for path in sorted(required_paths)
        if path.is_file()
    }


def _report_text(summary: Mapping[str, Any]) -> str:
    statuses = summary["statuses"]
    inventory = summary.get(
        "inventory",
        {
            "inventory_count": 0,
            "minimum_slot_count_through_guard": 0,
            "maximum_physical_nullity": 0,
            "unresolved_low_sigma_count": 0,
            "maximum_primary_verification_relative": math.inf,
        },
    )
    residuals = summary.get(
        "representative_joint_states",
        {
            "maximum_joint_compatibility_normalized": math.inf,
            "maximum_joint_equilibrium_normalized": math.inf,
            "maximum_outer_clamp_normalized": math.inf,
            "maximum_boundary_null_residual": math.inf,
            "maximum_energy_identity_relative": math.inf,
        },
    )
    reconciliation = inventory.get("reference_detector_reconciliation", {})
    reconciled_pairs = sum(
        int(item["reconciled_cross_method_count"])
        for item in reconciliation.values()
    )
    ambiguous_pairs = sum(
        int(item["ambiguous_cross_method_count"])
        for item in reconciliation.values()
    )
    return f"""# RLB-1: жёсткий узел и прямой предел beta=0

## Область проверки

Рассматриваются два симметрично слоистых стержня Reddy с идеальным точечным узлом. Спектральные расчёты выполнены только при `beta=0`. Ненулевые углы использованы в алгебраической проверке матрицы узла и граничной работы.

## Матрица узла и граничная работа

Инвариантный построитель использует физические отображения перемещений, поворотов, сил и моментов. Независимый построитель использует коэффициенты `cos(beta)` и `sin(beta)`. Максимальное расхождение матриц равно `{summary['joint']['maximum_invariant_closed_form_absolute_difference']:.6e}`. Максимальный нормированный остаток граничной работы равен `{summary['virtual_work']['maximum_normalized_residual']:.6e}`.

## Политика поиска корней

Первичный поиск не использует начальные приближения из опорных спектров. Он объединяет интервалы смены знака определителя и локальные минимумы отношения `sigma_min/sigma_max`. Первые 12 положительных корней дополнены корнем 13, который служит проверкой полноты. Опорные fixed--fixed и stepped inventories построены после фиксации хешей первичных coupled inventories.

Все `{inventory['inventory_count']}` coupled и reference inventories содержат по `{inventory['minimum_slot_count_through_guard']}` спектральных мест до guard включительно. В этом инвентаре физические roots изолированы, максимальная nullity равна `{inventory['maximum_physical_nullity']}`, а число неразрешённых low-sigma кандидатов до guard равно `{inventory['unresolved_low_sigma_count']}`. Максимальное расхождение primary и verification inventories равно `{inventory['maximum_primary_verification_relative']:.6e}`. Отклонённые кандидаты выше guard сохраняются отдельно и не входят в этот gate.

В post-freeze axial-plus-bending references согласованы `{reconciled_pairs}` уникальных cross-method пар RLB-0; число неоднозначных пар равно `{ambiguous_pairs}`. Два determinant-кандидата этим правилом не объединяются.

## Прямой однородный предел

Проверены равное и неравное разбиения общей длины для ламинатов `0_deg` и `[0/90/90/0]`. Coupled inventory сопоставлен с одной балкой fixed--fixed и с объединением независимых axial и bending families. Максимальное относительное расхождение равно `{summary['homogeneous']['maximum_relative_difference']:.6e}`.

## Материальный интерфейс

Для порядка `0_deg|cross_ply` coupled inventory сопоставлен с независимым stepped builder. Отдельно проверено зеркальное изменение порядка материалов. Максимальное относительное расхождение равно `{summary['stepped']['maximum_relative_difference']:.6e}`.

## Представительные состояния узла

Для шести изолированных roots восстановлены физические состояния обоих плеч. Максимальные нормированные невязки совместности, равновесия и внешних заделок равны соответственно `{residuals['maximum_joint_compatibility_normalized']:.6e}`, `{residuals['maximum_joint_equilibrium_normalized']:.6e}` и `{residuals['maximum_outer_clamp_normalized']:.6e}`. Максимальная граничная null-невязка равна `{residuals['maximum_boundary_null_residual']:.6e}`, а максимальная относительная невязка энергетического тождества — `{residuals['maximum_energy_identity_relative']:.6e}`.

## Статусы

- `{STATUS_JOINT}`: `{statuses[STATUS_JOINT]}`;
- `{STATUS_HOMOGENEOUS}`: `{statuses[STATUS_HOMOGENEOUS]}`;
- `{STATUS_STEPPED}`: `{statuses[STATUS_STEPPED]}`;
- `{STATUS_INVENTORY}`: `{statuses[STATUS_INVENTORY]}`;
- `{STATUS_OVERALL}`: `{statuses[STATUS_OVERALL]}`.

## Ограничения

Полученный результат является диагностической coupled-joint baseline при `beta=0`. Он не подтверждает спектр углового узла. На этом этапе не применялись two-arm Ritz solver, FEM, torsion, damping, complex roots, branch tracking, parameter sweep и построение фигур.
"""


def _git_state() -> dict[str, str]:
    def command(*arguments: str) -> str:
        completed = subprocess.run(
            ["git", *arguments], cwd=ROOT, check=True, capture_output=True, text=True
        )
        return completed.stdout.rstrip("\r\n")

    return {
        "top_level": command("rev-parse", "--show-toplevel").replace("\\", "/"),
        "branch": command("branch", "--show-current"),
        "head": command("rev-parse", "HEAD"),
        "last_commit": command("log", "-1", "--oneline"),
        "status_short": command("status", "--short", "--untracked-files=all"),
    }


def _all_inventory_candidates(inventories: Iterable[RootInventory]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    accepted: list[dict[str, Any]] = []
    rejected: list[dict[str, Any]] = []
    for inventory in inventories:
        for scan in (inventory.primary, inventory.verification):
            accepted.extend(_candidate_row(candidate) for candidate in scan.candidates)
            rejected.extend(_candidate_row(candidate) for candidate in scan.rejected_candidates)
    return accepted, rejected


def _status_from_comparisons(
    comparison_summaries: Sequence[Mapping[str, Any]],
    inventories: Sequence[RootInventory],
) -> str:
    if not all(bool(summary["pass"]) for summary in comparison_summaries):
        return "FAIL"
    if any(inventory.status == "FAIL" for inventory in inventories):
        return "FAIL"
    if all(inventory.status == "PASS" for inventory in inventories):
        return "PASS"
    return "PARTIAL_PASS"


def _execute_full_pilot(
    output_dir: Path,
    policy: SearchPolicy,
    model_manifest: Mapping[str, Any],
    selected: Mapping[str, Any],
) -> dict[str, Any]:
    coupled = _coupled_module()
    run_start_git = _git_state()
    preservation_before = _preservation_hashes()

    joint_rows, joint_summary = _joint_matrix_checks(coupled)
    virtual_rows, virtual_summary = _virtual_work_checks(coupled)
    joint_pass = bool(joint_summary["pass"] and virtual_summary["pass"])
    _write_csv(output_dir / "joint_matrix_checks.csv", joint_rows)
    _write_csv(output_dir / "virtual_work_checks.csv", virtual_rows)

    if not joint_pass:
        for filename in (
            "homogeneous_beta0_roots.csv", "homogeneous_beta0_comparison.csv",
            "stepped_beta0_roots.csv", "stepped_beta0_comparison.csv",
            "reflection_checks.csv", "representative_joint_residuals.csv",
            "root_candidates.csv", "rejected_candidates.csv",
        ):
            _write_csv(output_dir / filename, [], fields=("status",))
        statuses = {
            STATUS_JOINT: "FAIL",
            STATUS_HOMOGENEOUS: "FAIL",
            STATUS_STEPPED: "FAIL",
            STATUS_INVENTORY: "FAIL",
            STATUS_OVERALL: "FAIL",
        }
        summary = {
            "statuses": statuses,
            "joint": joint_summary,
            "virtual_work": virtual_summary,
            "homogeneous": {"maximum_relative_difference": math.inf},
            "stepped": {"maximum_relative_difference": math.inf},
        }
        (output_dir / "report.md").write_text(_report_text(summary), encoding="utf-8")
        run_manifest = {
            "algorithm_version": ALGORITHM_VERSION,
            "execution_mode": "full-stopped-after-joint-gate",
            "statuses": statuses,
            "initial_git": dict(TASK_INITIAL_GIT_STATE),
            "run_start_git": run_start_git,
            "final_git": _git_state(),
            "preservation_before": preservation_before,
            "preservation_after": _preservation_hashes(),
            "protected_file_count": len(preservation_before),
            "summary": summary,
        }
        _write_json(output_dir / "run_manifest.json", run_manifest)
        return run_manifest

    homogeneous_specs, stepped_specs = _case_specs(selected)
    all_specs = {spec.case_id: spec for spec in (*homogeneous_specs, *stepped_specs)}

    # Primary coupled inventories are completed and hashed before any direct
    # fixed--fixed, stepped, or separated-family reference is constructed.
    coupled_inventories: dict[str, RootInventory] = {}
    execution_order: list[dict[str, Any]] = []
    for spec in (*homogeneous_specs, *stepped_specs):
        inventory = seed_free_root_inventory(
            _coupled_provider(coupled, spec),
            spec.frequency_scale,
            policy,
            case_id=spec.case_id,
            builder_id="coupled_physical_J_block_assembly",
        )
        coupled_inventories[spec.case_id] = inventory
        execution_order.append(
            {
                "step": "primary_coupled_inventory_frozen",
                "case_id": spec.case_id,
                "inventory_sha256": inventory.inventory_sha256,
                "reference_builder_called_before_this_step": False,
            }
        )

    primary_freeze_payload = {
        case_id: inventory.inventory_sha256
        for case_id, inventory in sorted(coupled_inventories.items())
    }
    primary_freeze_sha256 = hashlib.sha256(
        json.dumps(primary_freeze_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest().upper()
    execution_order.append(
        {
            "step": "all_primary_coupled_inventories_frozen",
            "aggregate_sha256": primary_freeze_sha256,
        }
    )

    direct_inventories: dict[str, RootInventory] = {}
    union_by_laminate: dict[str, tuple[list[dict[str, Any]], dict[str, Any]]] = {}
    for spec in homogeneous_specs:
        direct = seed_free_root_inventory(
            _direct_homogeneous_provider(coupled, spec),
            spec.frequency_scale,
            policy,
            case_id=spec.case_id,
            builder_id="direct_one_transfer_fixed_fixed",
        )
        direct_inventories[spec.case_id] = direct
        laminate_id = spec.order_id.split("|", 1)[0]
        if laminate_id not in union_by_laminate:
            union_by_laminate[laminate_id] = _union_reference(selected[laminate_id], policy)
    for spec in stepped_specs:
        direct_inventories[spec.case_id] = seed_free_root_inventory(
            _direct_stepped_provider(coupled, spec),
            spec.frequency_scale,
            policy,
            case_id=spec.case_id,
            builder_id="direct_stepped_R0_no_joint_builder",
        )
    execution_order.append(
        {
            "step": "post_freeze_direct_references_completed",
            "reference_roots_used_as_seeds": False,
        }
    )

    homogeneous_comparison_rows: list[dict[str, Any]] = []
    homogeneous_summaries: list[dict[str, Any]] = []
    for spec in homogeneous_specs:
        rows, summary = _compare_inventories(
            coupled_inventories[spec.case_id],
            direct_inventories[spec.case_id],
            case_id=spec.case_id,
            comparison_kind="coupled_vs_direct_one_transfer_fixed_fixed",
        )
        homogeneous_comparison_rows.extend(rows)
        homogeneous_summaries.append({**summary, "case_id": spec.case_id})
        laminate_id = spec.order_id.split("|", 1)[0]
        union_rows, _union_summary = union_by_laminate[laminate_id]
        union_comparison, union_summary = _compare_inventory_to_union(
            direct_inventories[spec.case_id], union_rows, case_id=spec.case_id
        )
        homogeneous_comparison_rows.extend(union_comparison)
        homogeneous_summaries.append({**union_summary, "case_id": spec.case_id})

    stepped_comparison_rows: list[dict[str, Any]] = []
    stepped_summaries: list[dict[str, Any]] = []
    for spec in stepped_specs:
        rows, summary = _compare_inventories(
            coupled_inventories[spec.case_id],
            direct_inventories[spec.case_id],
            case_id=spec.case_id,
            comparison_kind="coupled_vs_direct_stepped_R0",
        )
        stepped_comparison_rows.extend(rows)
        stepped_summaries.append({**summary, "case_id": spec.case_id})

    reflection_rows, reflection_coupled_summary = _compare_inventories(
        coupled_inventories[stepped_specs[0].case_id],
        coupled_inventories[stepped_specs[1].case_id],
        case_id="reflection__coupled",
        comparison_kind="coupled_0_cross_vs_cross_0",
    )
    direct_reflection_rows, reflection_direct_summary = _compare_inventories(
        direct_inventories[stepped_specs[0].case_id],
        direct_inventories[stepped_specs[1].case_id],
        case_id="reflection__direct_stepped",
        comparison_kind="direct_0_cross_vs_cross_0",
    )
    reflection_rows.extend(direct_reflection_rows)
    stepped_summaries.extend(
        (
            {**reflection_coupled_summary, "case_id": "reflection__coupled"},
            {**reflection_direct_summary, "case_id": "reflection__direct_stepped"},
        )
    )

    representative_rows: list[dict[str, Any]] = []
    for laminate_id in ("0_deg", "cross_ply_0_90_s"):
        spec = next(
            item
            for item in homogeneous_specs
            if item.split_id == "equal" and item.order_id.startswith(laminate_id + "|")
        )
        union_rows, _union_summary = union_by_laminate[laminate_id]
        for family in ("axial", "bending"):
            reference = next(
                (
                    row for row in union_rows
                    if row["family"] == family and int(row["cluster_multiplicity"]) == 1
                ),
                None,
            )
            if reference is None:
                continue
            slot_index = int(reference["sorted_slot"]) - 1
            if slot_index < len(coupled_inventories[spec.case_id].slots):
                representative_rows.append(
                    _representative_residual(
                        coupled,
                        spec,
                        coupled_inventories[spec.case_id],
                        coupled_inventories[spec.case_id].slots[slot_index],
                        family=family,
                    )
                )
    primary_stepped = stepped_specs[0]
    isolated_stepped = [
        slot
        for slot in coupled_inventories[primary_stepped.case_id].slots
        if slot.event.cluster_multiplicity == 1
    ][:2]
    representative_rows.extend(
        _representative_residual(
            coupled,
            primary_stepped,
            coupled_inventories[primary_stepped.case_id],
            slot,
            family="stepped_unclassified_isolated",
        )
        for slot in isolated_stepped
    )

    homogeneous_inventories = [
        inventory
        for spec in homogeneous_specs
        for inventory in (
            coupled_inventories[spec.case_id], direct_inventories[spec.case_id]
        )
    ]
    stepped_inventories = [
        inventory
        for spec in stepped_specs
        for inventory in (
            coupled_inventories[spec.case_id], direct_inventories[spec.case_id]
        )
    ]
    all_inventories = [*homogeneous_inventories, *stepped_inventories]

    homogeneous_rep_pass = len([row for row in representative_rows if row["case_id"].startswith("homogeneous")]) >= 4 and all(
        row["status"] == "PASS"
        for row in representative_rows
        if row["case_id"].startswith("homogeneous")
    )
    stepped_rep_pass = len([row for row in representative_rows if row["case_id"].startswith("stepped")]) >= 2 and all(
        row["status"] == "PASS"
        for row in representative_rows
        if row["case_id"].startswith("stepped")
    )

    homogeneous_status = _status_from_comparisons(
        homogeneous_summaries, homogeneous_inventories
    )
    if not homogeneous_rep_pass:
        homogeneous_status = "FAIL"
    stepped_status = _status_from_comparisons(stepped_summaries, stepped_inventories)
    if not stepped_rep_pass:
        stepped_status = "FAIL"
    inventory_status = (
        "FAIL"
        if any(inventory.status == "FAIL" for inventory in all_inventories)
        else (
            "PASS"
            if all(inventory.status == "PASS" for inventory in all_inventories)
            else "PARTIAL_PASS"
        )
    )
    statuses = {
        STATUS_JOINT: "PASS",
        STATUS_HOMOGENEOUS: homogeneous_status,
        STATUS_STEPPED: stepped_status,
        STATUS_INVENTORY: inventory_status,
    }
    statuses[STATUS_OVERALL] = (
        "PASS"
        if all(statuses[key] == "PASS" for key in (STATUS_JOINT, STATUS_HOMOGENEOUS, STATUS_STEPPED, STATUS_INVENTORY))
        else "FAIL"
    )

    homogeneous_root_rows = [
        _slot_row(spec, inventory, slot)
        for spec in homogeneous_specs
        for inventory in (coupled_inventories[spec.case_id], direct_inventories[spec.case_id])
        for slot in inventory.slots
    ]
    stepped_root_rows = [
        _slot_row(spec, inventory, slot)
        for spec in stepped_specs
        for inventory in (coupled_inventories[spec.case_id], direct_inventories[spec.case_id])
        for slot in inventory.slots
    ]
    candidate_rows, rejected_rows = _all_inventory_candidates(all_inventories)

    _write_csv(output_dir / "homogeneous_beta0_roots.csv", homogeneous_root_rows)
    _write_csv(output_dir / "homogeneous_beta0_comparison.csv", homogeneous_comparison_rows)
    _write_csv(output_dir / "stepped_beta0_roots.csv", stepped_root_rows)
    _write_csv(output_dir / "stepped_beta0_comparison.csv", stepped_comparison_rows)
    _write_csv(output_dir / "reflection_checks.csv", reflection_rows)
    _write_csv(output_dir / "representative_joint_residuals.csv", representative_rows)
    _write_csv(output_dir / "root_candidates.csv", candidate_rows)
    _write_csv(output_dir / "rejected_candidates.csv", rejected_rows)
    _write_json(
        output_dir / "root_boundary_matrices.json",
        {"records": _matrix_records(all_inventories)},
    )

    homogeneous_max = max(
        (float(item["maximum_relative_difference"]) for item in homogeneous_summaries),
        default=math.inf,
    )
    stepped_max = max(
        (float(item["maximum_relative_difference"]) for item in stepped_summaries),
        default=math.inf,
    )
    inventory_details = {
        f"{inventory.builder_id}:{inventory.case_id}": {
            "status": inventory.status,
            "slot_count_through_guard": len(inventory.slots),
            "guard_omega_bar": (
                inventory.slots[policy.required_slots - 1].event.omega_bar
                if len(inventory.slots) >= policy.required_slots
                else math.nan
            ),
            "physical_cluster_count": len(
                {
                    slot.event.cluster_id
                    for slot in inventory.slots
                    if slot.event.cluster_semantics != "ISOLATED"
                }
            ),
            "maximum_physical_nullity": max(
                (slot.event.cluster_total_nullity for slot in inventory.slots),
                default=0,
            ),
            "maximum_primary_verification_relative": (
                inventory.maximum_primary_verification_relative
            ),
            "unresolved_low_sigma_count": inventory.unresolved_low_sigma_count,
        }
        for inventory in all_inventories
    }
    representative_summary = {
        "case_count": len(representative_rows),
        "maximum_joint_compatibility_normalized": max(
            (float(row["joint_compatibility_normalized"]) for row in representative_rows),
            default=math.inf,
        ),
        "maximum_joint_equilibrium_normalized": max(
            (float(row["joint_equilibrium_normalized"]) for row in representative_rows),
            default=math.inf,
        ),
        "maximum_outer_clamp_normalized": max(
            (float(row["outer_clamp_residual_normalized"]) for row in representative_rows),
            default=math.inf,
        ),
        "maximum_boundary_null_residual": max(
            (float(row["boundary_null_residual"]) for row in representative_rows),
            default=math.inf,
        ),
        "maximum_energy_identity_relative": max(
            (float(row["energy_identity_relative"]) for row in representative_rows),
            default=math.inf,
        ),
        "pass": bool(representative_rows) and all(
            row["status"] == "PASS" for row in representative_rows
        ),
    }
    summary = {
        "statuses": statuses,
        "joint": joint_summary,
        "virtual_work": virtual_summary,
        "homogeneous": {
            "maximum_relative_difference": homogeneous_max,
            "comparisons": homogeneous_summaries,
            "representative_states_pass": homogeneous_rep_pass,
        },
        "stepped": {
            "maximum_relative_difference": stepped_max,
            "comparisons": stepped_summaries,
            "representative_states_pass": stepped_rep_pass,
        },
        "inventory": {
            "primary_freeze_sha256": primary_freeze_sha256,
            "primary_inventory_hashes": primary_freeze_payload,
            "inventory_count": len(all_inventories),
            "minimum_slot_count_through_guard": min(
                (len(inventory.slots) for inventory in all_inventories), default=0
            ),
            "maximum_slot_count_through_guard": max(
                (len(inventory.slots) for inventory in all_inventories), default=0
            ),
            "physical_cluster_count": sum(
                int(detail["physical_cluster_count"])
                for detail in inventory_details.values()
            ),
            "maximum_physical_nullity": max(
                (int(detail["maximum_physical_nullity"]) for detail in inventory_details.values()),
                default=0,
            ),
            "unresolved_low_sigma_count": sum(
                int(detail["unresolved_low_sigma_count"])
                for detail in inventory_details.values()
            ),
            "reference_detector_reconciliation": {
                laminate_id: union_by_laminate[laminate_id][1][
                    "reference_detector_reconciliation"
                ]
                for laminate_id in sorted(union_by_laminate)
            },
            "case_details": inventory_details,
            "maximum_primary_verification_relative": max(
                inventory.maximum_primary_verification_relative for inventory in all_inventories
            ),
        },
        "representative_joint_states": representative_summary,
        "representative_joint_state_count": len(representative_rows),
    }
    (output_dir / "report.md").write_text(_report_text(summary), encoding="utf-8")

    preservation_after = _preservation_hashes()
    if preservation_after != preservation_before:
        raise RuntimeError("A protected RLB-0/source/coordinate artifact changed during the pilot")
    final_git = _git_state()
    output_hashes = {
        path.name: _sha256(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file() and path.name != "run_manifest.json"
    }
    run_manifest = {
        "algorithm_version": ALGORITHM_VERSION,
        "execution_mode": "full-beta0-pilot",
        "diagnostic_scope": "diagnostic beta=0 coupled-joint baseline",
        "original_preimplementation_manifest_sha256": (
            ORIGINAL_PREIMPLEMENTATION_MANIFEST_SHA256
        ),
        "precomputation_threshold_freeze_sha256": _active_contract_hash(policy),
        "active_threshold_search_contract_sha256": _active_contract_hash(policy),
        "model_manifest_sha256": _sha256(output_dir / "model_manifest.json"),
        "thresholds": THRESHOLDS,
        "search_policy": asdict(policy),
        "execution_order": execution_order,
        "reference_roots_used_as_seeds": False,
        "beta_nonzero_spectral_calculations": False,
        "initial_git": dict(TASK_INITIAL_GIT_STATE),
        "run_start_git": run_start_git,
        "final_git": final_git,
        "preservation_before": preservation_before,
        "preservation_after": preservation_after,
        "protected_file_count": len(preservation_before),
        "protected_files_unchanged": True,
        "statuses": statuses,
        "summary": summary,
        "outputs_sha256": output_hashes,
        "explicit_nonexecuted": [
            "beta_nonzero_spectrum", "two_arm_Ritz", "FEM", "torsion", "damping",
            "complex_roots", "parameter_study", "figures", "commit", "push",
        ],
    }
    _write_json(output_dir / "run_manifest.json", run_manifest)
    return run_manifest


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = run_pilot(args.output_dir, manifest_only=args.manifest_only)
    statuses = manifest.get("statuses", {})
    for key in (STATUS_JOINT, STATUS_HOMOGENEOUS, STATUS_STEPPED, STATUS_INVENTORY, STATUS_OVERALL):
        if key in statuses:
            print(f"{key}: {statuses[key]}")
    return 0 if statuses.get(STATUS_OVERALL) != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
