from __future__ import annotations

from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from scipy.optimize import brentq, minimize_scalar

from scripts.lib import variable_length_timoshenko as TIMO


ALGORITHM_VERSION = "factorized_straight_spectrum_v2"

# Frozen extraction from the existing Timoshenko unknown ordering
# (A1, B1, P1, A2, B2, P2).  These indices are an exact beta=0, eta=0
# block permutation of TIMO.timo_coupling_matrix, not a new determinant.
BENDING_ROWS = (0, 2, 3, 4)
BENDING_COLUMNS = (0, 1, 3, 4)
AXIAL_ROWS = (1, 5)
AXIAL_COLUMNS = (2, 5)

# Audit-local settings.  They deliberately do not modify the shared solver.
BENDING_SCAN_START = TIMO.ROOT_SCAN_START
BENDING_SCAN_STEP = TIMO.ROOT_SCAN_STEP
# The independent sign and SVD paths can converge a few micro-units apart in
# ill-conditioned high-mu coordinates.  This audit-local tolerance merges only
# roots inside the bending family; cross-family roots are never deduplicated.
BENDING_DEDUP_TOL = 1.0e-5
BENDING_BRENT_XTOL = 1.0e-12
BENDING_BRENT_RTOL = 1.0e-12
BENDING_SVD_PREFILTER = 2.0e-3
BENDING_SVD_ACCEPT = 2.0e-8
BENDING_SVD_XATOL = 1.0e-12
BENDING_GROW_FACTOR = 1.35
BENDING_MAX_GROWTH_TRIES = 5
BLOCK_OFFDIAGONAL_ABS_TOL = 1.0e-12
AXIAL_BLOCK_SVD_ACCEPT = 2.0e-8
FULL_MATRIX_SVD_ACCEPT = 5.0e-8
FULL_MATRIX_SVD_ROW_SCALE_OFFSET = 0.25 * TIMO.ROOT_SCAN_STEP
CROSS_FAMILY_CLUSTER_ABS_TOL = TIMO.ROOT_SCAN_STEP


@dataclass(frozen=True)
class FamilyRoot:
    value: float
    family: str
    family_index: int
    detection_source: str
    block_sigma_min: float
    multiplicity_group: str = ""
    within_family_duplicate: bool = False
    cross_family_cluster: bool = False
    cluster_gap: float = float("nan")


@dataclass(frozen=True)
class BendingRootSearch:
    roots: tuple[float, ...]
    sources: tuple[str, ...]
    sigma_minima: tuple[float, ...]
    scan_upper: float
    determinant_evaluations: int
    svd_refinements: int
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class FactorizedSpectrum:
    epsilon: float
    mu: float
    requested_roots: int
    roots: tuple[FamilyRoot, ...]
    bending_search: BendingRootSearch
    cache_status: str
    algorithm_version: str = ALGORITHM_VERSION

    @property
    def values(self) -> tuple[float, ...]:
        return tuple(item.value for item in self.roots)

    @property
    def families(self) -> tuple[str, ...]:
        return tuple(item.family for item in self.roots)


def exact_axial_roots(epsilon: float, count: int) -> np.ndarray:
    epsilon_f = float(epsilon)
    if epsilon_f <= 0.0:
        raise ValueError("epsilon must be positive")
    indices = np.arange(1, int(count) + 1, dtype=float)
    return np.sqrt(indices * np.pi / (2.0 * epsilon_f))


def straight_blocks(
    Lambda: float,
    mu: float,
    epsilon: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]]:
    matrix, warnings = TIMO.timo_coupling_matrix(
        float(Lambda),
        0.0,
        float(mu),
        float(epsilon),
        0.0,
    )
    bending = matrix[np.ix_(BENDING_ROWS, BENDING_COLUMNS)]
    axial = matrix[np.ix_(AXIAL_ROWS, AXIAL_COLUMNS)]
    return matrix, bending, axial, tuple(warnings)


def off_block_max_abs(matrix: np.ndarray) -> float:
    allowed = np.zeros((6, 6), dtype=bool)
    allowed[np.ix_(BENDING_ROWS, BENDING_COLUMNS)] = True
    allowed[np.ix_(AXIAL_ROWS, AXIAL_COLUMNS)] = True
    outside = np.abs(np.asarray(matrix, dtype=float)[~allowed])
    return float(np.max(outside)) if outside.size else 0.0


def singular_values(matrix: np.ndarray) -> tuple[float, float]:
    values = np.linalg.svd(TIMO.row_normalized(np.asarray(matrix, dtype=float)), compute_uv=False)
    smallest = float(values[-1]) if values.size else float("nan")
    second = float(values[-2]) if values.size >= 2 else float("nan")
    return smallest, second


def fixed_row_scale_singular_values(
    matrix: np.ndarray,
    reference_matrix: np.ndarray,
) -> tuple[float, float]:
    """Scale by nearby nonsingular row norms without amplifying a zero row."""

    scaled = np.asarray(matrix, dtype=float).copy()
    reference_norms = np.linalg.norm(np.asarray(reference_matrix, dtype=float), axis=1)
    finite_positive = reference_norms[np.isfinite(reference_norms) & (reference_norms > 0.0)]
    fallback = float(np.max(finite_positive)) if finite_positive.size else 1.0
    for index, norm in enumerate(reference_norms):
        divisor = float(norm) if math.isfinite(float(norm)) and float(norm) > 0.0 else fallback
        scaled[index, :] /= divisor
    values = np.linalg.svd(scaled, compute_uv=False)
    largest = float(values[0]) if values.size else float("nan")
    if not math.isfinite(largest) or largest <= 0.0:
        return float("nan"), float("nan")
    smallest = float(values[-1] / largest)
    second = float(values[-2] / largest) if values.size >= 2 else float("nan")
    return smallest, second


def bending_block_det(Lambda: float, mu: float, epsilon: float) -> float:
    try:
        _matrix, bending, _axial, _warnings = straight_blocks(Lambda, mu, epsilon)
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan")
    if not np.all(np.isfinite(bending)):
        return float("nan")
    return TIMO.normalized_det(bending)


def bending_block_singular_values(Lambda: float, mu: float, epsilon: float) -> tuple[float, float]:
    try:
        _matrix, bending, _axial, _warnings = straight_blocks(Lambda, mu, epsilon)
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan"), float("nan")
    if not np.all(np.isfinite(bending)):
        return float("nan"), float("nan")
    return singular_values(bending)


def axial_block_singular_values(Lambda: float, mu: float, epsilon: float) -> tuple[float, float]:
    try:
        _matrix, _bending, axial, _warnings = straight_blocks(Lambda, mu, epsilon)
        _reference_matrix, _reference_bending, reference_axial, _reference_warnings = straight_blocks(
            float(Lambda) + FULL_MATRIX_SVD_ROW_SCALE_OFFSET,
            mu,
            epsilon,
        )
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan"), float("nan")
    if not np.all(np.isfinite(axial)):
        return float("nan"), float("nan")
    return fixed_row_scale_singular_values(axial, reference_axial)


def axial_block_det(Lambda: float, mu: float, epsilon: float) -> float:
    try:
        _matrix, _bending, axial, _warnings = straight_blocks(Lambda, mu, epsilon)
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan")
    return float(np.linalg.det(axial)) if np.all(np.isfinite(axial)) else float("nan")


def full_matrix_singular_values(Lambda: float, mu: float, epsilon: float) -> tuple[float, float]:
    try:
        matrix, _bending, _axial, _warnings = straight_blocks(Lambda, mu, epsilon)
        reference, _reference_bending, _reference_axial, _reference_warnings = straight_blocks(
            float(Lambda) + FULL_MATRIX_SVD_ROW_SCALE_OFFSET,
            mu,
            epsilon,
        )
    except (FloatingPointError, OverflowError, ValueError, np.linalg.LinAlgError):
        return float("nan"), float("nan")
    if not np.all(np.isfinite(matrix)):
        return float("nan"), float("nan")
    return fixed_row_scale_singular_values(matrix, reference)


def _deduplicate_candidates(
    candidates: Sequence[tuple[float, str, float]],
    tolerance: float,
) -> list[tuple[float, str, float]]:
    result: list[tuple[float, str, float]] = []
    for value, source, sigma in sorted(candidates, key=lambda item: (item[0], item[1])):
        if result and abs(value - result[-1][0]) <= float(tolerance):
            previous = result[-1]
            merged_sources = "+".join(sorted(set(previous[1].split("+") + source.split("+"))))
            result[-1] = (
                value if sigma < previous[2] else previous[0],
                merged_sources,
                min(sigma, previous[2]),
            )
        else:
            result.append((float(value), str(source), float(sigma)))
    return result


def _refine_svd_minimum(
    func: Callable[[float], tuple[float, float]],
    left: float,
    right: float,
) -> tuple[float, float]:
    result = minimize_scalar(
        lambda value: func(float(value))[0],
        bounds=(float(left), float(right)),
        method="bounded",
        options={"xatol": BENDING_SVD_XATOL, "maxiter": 160},
    )
    if not result.success or not math.isfinite(float(result.fun)):
        return float("nan"), float("nan")
    return float(result.x), float(result.fun)


def find_bending_roots(
    epsilon: float,
    mu: float,
    count: int,
    *,
    scan_start: float = BENDING_SCAN_START,
    scan_step: float = BENDING_SCAN_STEP,
    upper: float | None = None,
) -> BendingRootSearch:
    """Find roots of the exact 4x4 block by signs plus local SVD minima."""

    count_i = int(count)
    if count_i <= 0:
        return BendingRootSearch((), (), (), float(scan_start), 0, 0, ())
    epsilon_f = float(epsilon)
    mu_f = float(mu)
    current_upper = float(
        upper
        if upper is not None
        else max(18.0, 0.5 * math.pi * (count_i + 1))
    )
    determinant_evaluations = 0
    svd_refinements = 0
    warnings: list[str] = []
    accepted: list[tuple[float, str, float]] = []

    det_func = lambda value: bending_block_det(value, mu_f, epsilon_f)
    svd_func = lambda value: bending_block_singular_values(value, mu_f, epsilon_f)

    for _attempt in range(BENDING_MAX_GROWTH_TRIES):
        grid = np.arange(float(scan_start), current_upper + float(scan_step), float(scan_step), dtype=float)
        values = np.asarray([det_func(float(value)) for value in grid], dtype=float)
        determinant_evaluations += len(grid)
        candidates: list[tuple[float, str, float]] = []

        for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if not (math.isfinite(float(f_left)) and math.isfinite(float(f_right))):
                continue
            candidate: float | None = None
            if float(f_left) == 0.0:
                candidate = float(left)
            elif float(f_left) * float(f_right) < 0.0:
                try:
                    candidate = float(
                        brentq(
                            det_func,
                            float(left),
                            float(right),
                            xtol=BENDING_BRENT_XTOL,
                            rtol=BENDING_BRENT_RTOL,
                            maxiter=160,
                        )
                    )
                except (ValueError, RuntimeError):
                    candidate = None
            if candidate is not None:
                sigma, _second = svd_func(candidate)
                if math.isfinite(sigma) and sigma <= BENDING_SVD_ACCEPT:
                    candidates.append((candidate, "sign_change", sigma))

        absolute = np.abs(values)
        for index in range(1, len(grid) - 1):
            local = float(absolute[index])
            if not math.isfinite(local) or local > BENDING_SVD_PREFILTER:
                continue
            if local > float(absolute[index - 1]) or local > float(absolute[index + 1]):
                continue
            svd_refinements += 1
            candidate, sigma = _refine_svd_minimum(
                svd_func,
                float(grid[index - 1]),
                float(grid[index + 1]),
            )
            if math.isfinite(candidate) and math.isfinite(sigma) and sigma <= BENDING_SVD_ACCEPT:
                candidates.append((candidate, "svd_minimum", sigma))

        accepted = _deduplicate_candidates(candidates, BENDING_DEDUP_TOL)
        if len(accepted) >= count_i:
            break
        current_upper *= BENDING_GROW_FACTOR

    if len(accepted) < count_i:
        warnings.append(
            f"found_only_{len(accepted)}_of_{count_i}_bending_roots_below_{current_upper:.8g}"
        )
    selected = accepted[:count_i]
    return BendingRootSearch(
        roots=tuple(item[0] for item in selected),
        sources=tuple(item[1] for item in selected),
        sigma_minima=tuple(item[2] for item in selected),
        scan_upper=float(current_upper),
        determinant_evaluations=int(determinant_evaluations),
        svd_refinements=int(svd_refinements),
        warnings=tuple(warnings),
    )


def _annotate_cross_family_clusters(roots: Sequence[FamilyRoot]) -> tuple[FamilyRoot, ...]:
    result = list(roots)
    group_index = 0
    for index in range(len(result) - 1):
        left = result[index]
        right = result[index + 1]
        gap = abs(right.value - left.value)
        if left.family == right.family or gap > CROSS_FAMILY_CLUSTER_ABS_TOL:
            continue
        group_index += 1
        group = f"cross_family_{group_index:03d}"
        result[index] = replace(
            left,
            multiplicity_group=group,
            cross_family_cluster=True,
            cluster_gap=gap,
        )
        result[index + 1] = replace(
            right,
            multiplicity_group=group,
            cross_family_cluster=True,
            cluster_gap=gap,
        )
    return tuple(result)


def factorized_straight_spectrum(
    epsilon: float,
    mu: float,
    count: int,
) -> FactorizedSpectrum:
    """Return the exact-family union for beta=0, eta=0 with multiplicity."""

    count_i = int(count)
    axial = exact_axial_roots(float(epsilon), count_i)
    # mu only moves the artificial internal joint in this straight homogeneous
    # rod.  Solve the bending family in the canonical mu=0 coordinates, then
    # verify every root in the requested-mu block below.  This removes a purely
    # coordinate-conditioning drift without importing a characteristic
    # equation or using a different physical model.
    bending = find_bending_roots(float(epsilon), 0.0, count_i)
    components: list[FamilyRoot] = []
    for family_index, value in enumerate(axial, start=1):
        sigma, _second = axial_block_singular_values(float(value), float(mu), float(epsilon))
        components.append(
            FamilyRoot(
                value=float(value),
                family="axial",
                family_index=family_index,
                detection_source="exact_axial_formula",
                block_sigma_min=sigma,
            )
        )
    for family_index, (value, source, sigma) in enumerate(
        zip(bending.roots, bending.sources, bending.sigma_minima),
        start=1,
    ):
        verified_sigma, _verified_second = bending_block_singular_values(
            float(value), float(mu), float(epsilon)
        )
        components.append(
            FamilyRoot(
                value=float(value),
                family="bending_Timoshenko",
                family_index=family_index,
                detection_source=(
                    source
                    if abs(float(mu)) <= 1.0e-14
                    else f"{source}+canonical_mu0_verified_at_mu"
                ),
                block_sigma_min=float(verified_sigma if math.isfinite(verified_sigma) else sigma),
            )
        )
    family_priority = {"axial": 0, "bending_Timoshenko": 1}
    components.sort(
        key=lambda item: (item.value, family_priority[item.family], item.family_index)
    )
    selected = _annotate_cross_family_clusters(components[:count_i])
    return FactorizedSpectrum(
        epsilon=float(epsilon),
        mu=float(mu),
        requested_roots=count_i,
        roots=selected,
        bending_search=bending,
        cache_status="miss",
    )


def _settings_payload(
    epsilon: float,
    mu: float,
    count: int,
    n_candidate_roots: int,
) -> dict[str, object]:
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "epsilon": float(epsilon),
        "mu": float(mu),
        "bending_search_mu": 0.0,
        "requested_mu_block_verification": True,
        "beta_deg": 0.0,
        "eta": 0.0,
        "n_spectrum_roots": int(count),
        "requested_roots": int(count),
        "n_candidate_roots": int(n_candidate_roots),
        "bending_rows": list(BENDING_ROWS),
        "bending_columns": list(BENDING_COLUMNS),
        "axial_rows": list(AXIAL_ROWS),
        "axial_columns": list(AXIAL_COLUMNS),
        "scan_start": BENDING_SCAN_START,
        "scan_step": BENDING_SCAN_STEP,
        "root_dedup_tol": BENDING_DEDUP_TOL,
        "brent_xtol": BENDING_BRENT_XTOL,
        "brent_rtol": BENDING_BRENT_RTOL,
        "svd_prefilter": BENDING_SVD_PREFILTER,
        "svd_accept": BENDING_SVD_ACCEPT,
        "svd_xatol": BENDING_SVD_XATOL,
        "initial_upper": max(18.0, 0.5 * math.pi * (int(count) + 1)),
        "grow_factor": BENDING_GROW_FACTOR,
        "max_growth_tries": BENDING_MAX_GROWTH_TRIES,
        "cross_family_cluster_abs_tol": CROSS_FAMILY_CLUSTER_ABS_TOL,
    }


def _cache_key(settings: dict[str, object]) -> str:
    encoded = json.dumps(settings, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _serialize(result: FactorizedSpectrum, settings: dict[str, object]) -> dict[str, object]:
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "settings": settings,
        "spectrum": {
            "epsilon": result.epsilon,
            "mu": result.mu,
            "requested_roots": result.requested_roots,
            "roots": [item.__dict__ for item in result.roots],
            "bending_search": result.bending_search.__dict__,
        },
    }


def _deserialize(payload: dict[str, object], cache_status: str) -> FactorizedSpectrum:
    spectrum = dict(payload["spectrum"])  # type: ignore[arg-type]
    roots = tuple(FamilyRoot(**dict(item)) for item in spectrum["roots"])  # type: ignore[arg-type]
    search_payload = dict(spectrum["bending_search"])  # type: ignore[arg-type]
    search = BendingRootSearch(
        roots=tuple(float(value) for value in search_payload["roots"]),
        sources=tuple(str(value) for value in search_payload["sources"]),
        sigma_minima=tuple(float(value) for value in search_payload["sigma_minima"]),
        scan_upper=float(search_payload["scan_upper"]),
        determinant_evaluations=int(search_payload["determinant_evaluations"]),
        svd_refinements=int(search_payload["svd_refinements"]),
        warnings=tuple(str(value) for value in search_payload["warnings"]),
    )
    return FactorizedSpectrum(
        epsilon=float(spectrum["epsilon"]),
        mu=float(spectrum["mu"]),
        requested_roots=int(spectrum["requested_roots"]),
        roots=roots,
        bending_search=search,
        cache_status=cache_status,
    )


@dataclass
class FactorizedSpectrumCache:
    directory: Path
    reuse: bool = True
    force_recompute: bool = False

    def path_for(
        self,
        epsilon: float,
        mu: float,
        count: int,
        n_candidate_roots: int | None = None,
    ) -> Path:
        settings = _settings_payload(
            epsilon,
            mu,
            count,
            int(n_candidate_roots if n_candidate_roots is not None else count),
        )
        return Path(self.directory) / f"{_cache_key(settings)}.json"

    def get(
        self,
        epsilon: float,
        mu: float,
        count: int,
        n_candidate_roots: int | None = None,
    ) -> FactorizedSpectrum:
        candidate_count = int(
            n_candidate_roots if n_candidate_roots is not None else count
        )
        settings = _settings_payload(epsilon, mu, count, candidate_count)
        path = self.path_for(epsilon, mu, count, candidate_count)
        stale = False
        if self.reuse and not self.force_recompute and path.exists():
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError, TypeError, ValueError):
                payload = {}
            if payload.get("algorithm_version") != ALGORITHM_VERSION:
                stale = True
            elif payload.get("settings") == settings:
                return _deserialize(payload, "hit")

        result = factorized_straight_spectrum(epsilon, mu, count)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            json.dumps(_serialize(result, settings), sort_keys=True, indent=2) + "\n",
            encoding="utf-8",
        )
        status = (
            "stale_cache_algorithm_version"
            if stale
            else ("force_recomputed" if self.force_recompute else "miss")
        )
        return replace(result, cache_status=status)


__all__ = [
    "ALGORITHM_VERSION",
    "AXIAL_BLOCK_SVD_ACCEPT",
    "AXIAL_COLUMNS",
    "AXIAL_ROWS",
    "BENDING_COLUMNS",
    "BENDING_ROWS",
    "BENDING_SVD_ACCEPT",
    "BLOCK_OFFDIAGONAL_ABS_TOL",
    "CROSS_FAMILY_CLUSTER_ABS_TOL",
    "FULL_MATRIX_SVD_ACCEPT",
    "BendingRootSearch",
    "FactorizedSpectrum",
    "FactorizedSpectrumCache",
    "FamilyRoot",
    "axial_block_singular_values",
    "axial_block_det",
    "bending_block_det",
    "bending_block_singular_values",
    "exact_axial_roots",
    "factorized_straight_spectrum",
    "fixed_row_scale_singular_values",
    "find_bending_roots",
    "full_matrix_singular_values",
    "off_block_max_abs",
    "singular_values",
    "straight_blocks",
]
