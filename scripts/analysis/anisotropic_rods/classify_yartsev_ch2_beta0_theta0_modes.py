"""Classify the fixed ``beta=theta=0`` Chapter-2 spectrum by exact blocks.

This narrow diagnostic keeps both rods, both external book-slope clamps, and
the unchanged rigid joint.  At the orthotropic straight limit the accepted
``state_corrected`` matrix separates exactly into a Timoshenko-bending block
and a generalized-torsion block.  Classification is made by membership in
those partial spectra; no energy criterion, MAC, or parameter continuation is
used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    Geometry,
    RodPoint,
    RootResult,
    _state_scales,
    cantilever_clamp_matrix,
    find_elastic_roots,
    hms_dx_209_material,
    make_rod_point,
    physical_state_trajectory_from_initial,
)


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_beta0_theta0_mode_classification"
)
REFERENCE_SPECTRUM_CSV = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_spectral_applicability_screening"
    / "theta2_small_grid"
    / "screening_spectra.csv"
)
ARTICLE_TEX = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_article_ru"
    / "drafts"
    / "article_figure_2_split"
    / "Dorofeev_references_english_revised_fig2_split.tex"
)

FORMULATION = "state_corrected"
CLAMP = "book_slope_clamp"
MATERIAL_MODE = "elastic"
BETA_DEG = 0.0
THETA_DEG = 0.0
REFERENCE_LENGTH_M = 0.4
REFERENCE_A_M = 0.005
REFERENCE_B_M = 0.020
REFERENCE_EX_PA = 191.0e9
SHEAR_FACTOR = 5.0 / 6.0
FULL_ROOT_COUNT = 7
PLOTTED_ROOT_COUNT = 6
PARTIAL_ROOT_COUNT = 7
SCAN_STEP_HZ = 0.5
INITIAL_MAX_HZ = 5_000.0
MAX_HZ = 100_000.0
ROOT_DETERMINANT_TOLERANCE = 1.0e-8
ROOT_SINGULAR_TOLERANCE = 1.0e-8
BLOCK_RELATIVE_TOLERANCE = 2.0e-14
SPECTRAL_UNION_RELATIVE_TOLERANCE = 1.0e-8
COMPONENT_SUPPRESSION_TOLERANCE = 1.0e-8
REFERENCE_REPEATABILITY_TOLERANCE = 1.0e-8
NEAR_COINCIDENCE_RELATIVE_TOLERANCE = 2.0e-8
KERNEL_RELATIVE_TOLERANCE = 1.0e-8
SHAPE_SAMPLES = 401

# Full D rows: joint conditions.  Full D columns: clamp reactions
# [Q1, M1, MT1, Q2, M2, MT2].
BENDING_ROW_INDICES = (0, 2, 3, 5)
TORSION_ROW_INDICES = (1, 4)
BENDING_COLUMN_INDICES = (0, 1, 3, 4)
TORSION_COLUMN_INDICES = (2, 5)
FULL_ROW_PERMUTATION = BENDING_ROW_INDICES + TORSION_ROW_INDICES
FULL_COLUMN_PERMUTATION = BENDING_COLUMN_INDICES + TORSION_COLUMN_INDICES
TWO_ARM_STATE_PERMUTATION = (0, 1, 3, 4, 6, 7, 9, 10, 2, 5, 8, 11)

BENDING_STATE_INDICES = (0, 1, 3, 4)
TORSION_STATE_INDICES = (2, 5)


@dataclass(frozen=True)
class CasePreset:
    case_id: str
    mu: float
    tau: float


@dataclass(frozen=True)
class SpectrumResult:
    roots: tuple[RootResult, ...]
    quality: tuple[dict[str, Any], ...]
    matrix_evaluations: int
    raw_quality_evaluations: int
    runtime_seconds: float


@dataclass(frozen=True)
class TaggedPartialRoot:
    family: str
    partial_mode_number: int
    root: RootResult
    quality: dict[str, Any]


CASES = (
    CasePreset("A", 0.0, 0.0),
    CasePreset("B", 0.5, -0.2),
)


class CountingBoundaryBuilder:
    def __init__(self, factory: Callable[[complex], NDArray[np.complex128]]) -> None:
        self.factory = factory
        self.evaluations = 0

    def __call__(
        self, omega: complex, _point: RodPoint, _formulation: str
    ) -> NDArray[np.complex128]:
        self.evaluations += 1
        return self.factory(omega)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify the two fixed beta=theta=0 article spectra by exact "
            "Chapter-2 bending/torsion block separation."
        )
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--reuse-data", action="store_true")
    group.add_argument("--force-recompute", action="store_true")
    return parser.parse_args(argv)


def article_geometry(mu: float, tau: float) -> tuple[Geometry, Geometry]:
    """Return the exact article volume-preserving rectangular geometry."""

    denominator = 1.0 + tau**2 + 2.0 * mu * tau
    if not (-1.0 < mu < 1.0) or denominator <= 0.0:
        raise ValueError("invalid article geometry parameters")
    scale = REFERENCE_A_M / math.sqrt(denominator)
    a_1 = scale * (1.0 - tau)
    a_2 = scale * (1.0 + tau)
    length_1 = REFERENCE_LENGTH_M * (1.0 - mu)
    length_2 = REFERENCE_LENGTH_M * (1.0 + mu)
    return (
        Geometry(
            a=a_1,
            b=4.0 * a_1,
            length=length_1,
            shear_factor=SHEAR_FACTOR,
        ),
        Geometry(
            a=a_2,
            b=4.0 * a_2,
            length=length_2,
            shear_factor=SHEAR_FACTOR,
        ),
    )


def case_points(case: CasePreset) -> tuple[RodPoint, RodPoint]:
    material = hms_dx_209_material()
    geometry_1, geometry_2 = article_geometry(case.mu, case.tau)
    return (
        make_rod_point(
            THETA_DEG,
            material_mode=MATERIAL_MODE,
            geometry=geometry_1,
            material=material,
        ),
        make_rod_point(
            THETA_DEG,
            material_mode=MATERIAL_MODE,
            geometry=geometry_2,
            material=material,
        ),
    )


def lambda_reference(omega_rad_s: float | NDArray[np.float64]) -> float | NDArray[np.float64]:
    """Use the fixed article normalization, never a geometry-dependent one."""

    omega = np.asarray(omega_rad_s, dtype=float)
    if np.any(omega < 0.0):
        raise ValueError("omega must be nonnegative")
    material = hms_dx_209_material()
    area_0 = REFERENCE_A_M * REFERENCE_B_M
    inertia_y_0 = REFERENCE_A_M**3 * REFERENCE_B_M / 12.0
    value = (
        material.rho
        * area_0
        * omega**2
        * REFERENCE_LENGTH_M**4
        / (REFERENCE_EX_PA * inertia_y_0)
    ) ** 0.25
    return float(value) if value.ndim == 0 else value


def split_raw_matrix(
    raw_matrix: NDArray[np.complex128],
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """Permute the accepted full matrix and return its two diagonal blocks."""

    raw = np.asarray(raw_matrix, dtype=np.complex128)
    if raw.shape != (6, 6):
        raise ValueError("the coupled Chapter-2 matrix must have shape (6, 6)")
    permuted = raw[np.ix_(FULL_ROW_PERMUTATION, FULL_COLUMN_PERMUTATION)]
    bending = permuted[:4, :4]
    torsion = permuted[4:, 4:]
    return permuted, bending, torsion


def matrix_factories(
    point_1: RodPoint, point_2: RodPoint
) -> dict[str, tuple[Callable[[complex], NDArray[np.complex128]], Callable[[complex], NDArray[np.complex128]]]]:
    def full_raw(omega: complex) -> NDArray[np.complex128]:
        return coupled_boundary_matrix_raw(omega, 0.0, point_1, point_2)

    def bending_raw(omega: complex) -> NDArray[np.complex128]:
        return split_raw_matrix(full_raw(omega))[1]

    def torsion_raw(omega: complex) -> NDArray[np.complex128]:
        return split_raw_matrix(full_raw(omega))[2]

    def scaled(factory: Callable[[complex], NDArray[np.complex128]]) -> Callable[[complex], NDArray[np.complex128]]:
        return lambda omega: equilibrate_matrix(factory(omega))[0]

    return {
        "full": (full_raw, scaled(full_raw)),
        "bending": (bending_raw, scaled(bending_raw)),
        "torsion": (torsion_raw, scaled(torsion_raw)),
    }


def matrix_quality(matrix: NDArray[np.complex128]) -> dict[str, float]:
    value = np.asarray(matrix, dtype=np.complex128)
    singular = np.linalg.svd(value, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    determinant_abs = float(abs(np.linalg.det(value)))
    determinant_scale = max(sigma_max ** value.shape[0], np.finfo(float).tiny)
    return {
        "determinant_abs": determinant_abs,
        "normalized_determinant_residual": determinant_abs / determinant_scale,
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "relative_singular_residual": sigma_min / sigma_max if sigma_max else 0.0,
    }


def root_quality(
    root: RootResult, raw_factory: Callable[[complex], NDArray[np.complex128]]
) -> dict[str, Any]:
    raw = matrix_quality(raw_factory(root.omega))
    status_ok = not root.status.startswith("rejected")
    scaled_ok = (
        root.determinant_residual <= ROOT_DETERMINANT_TOLERANCE
        and root.relative_singular_residual <= ROOT_SINGULAR_TOLERANCE
        and status_ok
    )
    raw_ok = (
        raw["normalized_determinant_residual"] <= ROOT_DETERMINANT_TOLERANCE
        and raw["relative_singular_residual"] <= ROOT_SINGULAR_TOLERANCE
        and status_ok
    )
    basis = "scaled" if scaled_ok else "physical_raw" if raw_ok else "none"
    return {
        "quality_status": "PASS" if scaled_ok or raw_ok else "FAIL",
        "quality_basis": basis,
        "root_status": root.status,
        "scaled_determinant_residual": root.determinant_residual,
        "scaled_relative_singular_residual": root.relative_singular_residual,
        "physical_raw_normalized_determinant_residual": raw[
            "normalized_determinant_residual"
        ],
        "physical_raw_relative_singular_residual": raw[
            "relative_singular_residual"
        ],
        "accepted_determinant_residual": (
            root.determinant_residual
            if scaled_ok
            else raw["normalized_determinant_residual"]
        ),
        "accepted_relative_singular_residual": (
            root.relative_singular_residual
            if scaled_ok
            else raw["relative_singular_residual"]
        ),
    }


def validate_roots(
    roots: Sequence[RootResult], quality: Sequence[Mapping[str, Any]], expected: int
) -> None:
    if len(roots) != expected or len(quality) != expected:
        raise RuntimeError(f"expected {expected} roots, got {len(roots)}")
    frequencies = np.asarray([item.frequency_hz for item in roots], dtype=float)
    if np.any(~np.isfinite(frequencies)) or np.any(frequencies <= 0.0):
        raise RuntimeError("spectrum has nonfinite or nonpositive roots")
    if np.any(np.diff(frequencies) <= 0.0):
        raise RuntimeError("spectrum is not strictly increasing")
    relative_gaps = np.diff(frequencies) / np.maximum(frequencies[:-1], frequencies[1:])
    if np.any(relative_gaps <= NEAR_COINCIDENCE_RELATIVE_TOLERANCE):
        raise RuntimeError("duplicate roots were found inside one partial family")
    failed = [index + 1 for index, item in enumerate(quality) if item["quality_status"] != "PASS"]
    if failed:
        raise RuntimeError(f"root-quality gate failed at positions {failed}")


def solve_spectrum(
    point_for_contract: RodPoint,
    raw_factory: Callable[[complex], NDArray[np.complex128]],
    scaled_factory: Callable[[complex], NDArray[np.complex128]],
    root_count: int,
) -> SpectrumResult:
    builder = CountingBoundaryBuilder(scaled_factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point_for_contract,
        FORMULATION,  # type: ignore[arg-type]
        num_roots=root_count,
        scan_step_hz=SCAN_STEP_HZ,
        initial_max_hz=INITIAL_MAX_HZ,
        max_hz=MAX_HZ,
        boundary_matrix_builder=builder,
    )
    quality = tuple(root_quality(root, raw_factory) for root in roots)
    validate_roots(roots, quality, root_count)
    return SpectrumResult(
        roots=tuple(roots),
        quality=quality,
        matrix_evaluations=builder.evaluations,
        raw_quality_evaluations=len(roots),
        runtime_seconds=time.perf_counter() - started,
    )


def merge_partial_spectra(
    bending: SpectrumResult, torsion: SpectrumResult
) -> tuple[TaggedPartialRoot, ...]:
    tagged = [
        TaggedPartialRoot("bending", index, root, quality)
        for index, (root, quality) in enumerate(zip(bending.roots, bending.quality), start=1)
    ]
    tagged.extend(
        TaggedPartialRoot("torsion", index, root, quality)
        for index, (root, quality) in enumerate(zip(torsion.roots, torsion.quality), start=1)
    )
    tagged.sort(key=lambda item: (item.root.frequency_hz, item.family))
    return tuple(tagged)


def cross_family_near_coincidences(
    merged: Sequence[TaggedPartialRoot],
) -> list[tuple[int, int, float]]:
    result: list[tuple[int, int, float]] = []
    for left_index, (left, right) in enumerate(zip(merged, merged[1:])):
        if left.family == right.family:
            continue
        relative = abs(left.root.frequency_hz - right.root.frequency_hz) / max(
            left.root.frequency_hz, right.root.frequency_hz
        )
        if relative <= NEAR_COINCIDENCE_RELATIVE_TOLERANCE:
            result.append((left_index, left_index + 1, relative))
    return result


def minimum_cross_family_gap(
    merged: Sequence[TaggedPartialRoot],
) -> tuple[float, float, str]:
    candidates: list[tuple[float, float, str]] = []
    for left, right in zip(merged, merged[1:]):
        if left.family == right.family:
            continue
        absolute = abs(left.root.frequency_hz - right.root.frequency_hz)
        relative = absolute / max(left.root.frequency_hz, right.root.frequency_hz)
        label = (
            f"{left.family}_{left.partial_mode_number}/"
            f"{right.family}_{right.partial_mode_number}"
        )
        candidates.append((absolute, relative, label))
    if not candidates:
        raise RuntimeError("partial union contains no cross-family neighbors")
    return min(candidates, key=lambda item: item[0])


def block_separation_row(
    case: CasePreset,
    frequency_hz: float,
    raw: NDArray[np.complex128],
) -> dict[str, Any]:
    permuted, _, _ = split_raw_matrix(raw)
    upper = permuted[:4, 4:]
    lower = permuted[4:, :4]
    off_values = np.concatenate((upper.ravel(), lower.ravel()))
    max_abs = float(np.max(np.abs(off_values))) if off_values.size else 0.0
    frobenius = float(np.linalg.norm(off_values))
    full_frobenius = float(np.linalg.norm(permuted))
    relative = frobenius / max(full_frobenius, np.finfo(float).tiny)
    return {
        "case_id": case.case_id,
        "mu": case.mu,
        "tau": case.tau,
        "beta_deg": BETA_DEG,
        "theta_deg": THETA_DEG,
        "frequency_hz": frequency_hz,
        "full_matrix_frobenius_norm": full_frobenius,
        "off_block_max_abs": max_abs,
        "off_block_frobenius_norm": frobenius,
        "off_block_relative_frobenius_norm": relative,
        "relative_tolerance": BLOCK_RELATIVE_TOLERANCE,
        "status": "PASS" if relative <= BLOCK_RELATIVE_TOLERANCE else "FAIL",
        "row_permutation": ",".join(map(str, FULL_ROW_PERMUTATION)),
        "column_permutation": ",".join(map(str, FULL_COLUMN_PERMUTATION)),
        "two_arm_state_permutation": ",".join(map(str, TWO_ARM_STATE_PERMUTATION)),
    }


def fixed_physical_scaling(
    point_1: RodPoint, point_2: RodPoint
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return frequency-independent row and reaction scales for full ``D``.

    Dynamic max-norm equilibration is appropriate for locating determinant
    sign changes, but it deliberately rescales a nearly zero row.  A null
    vector at an already identified root therefore uses fixed physical state
    scales instead.  This preserves the small row that identifies the partial
    block and avoids mixing quantities with different units.
    """

    scales_1 = _state_scales(point_1)
    scales_2 = _state_scales(point_2)
    row_scales = np.asarray(
        [
            max(scales_1[0], scales_2[0]),  # w continuity
            1.0,  # Phi compatibility
            1.0,  # psi compatibility
            max(scales_1[3], scales_2[3]),  # Q equilibrium
            max(scales_1[5], scales_2[5]),  # M_T equilibrium
            max(scales_1[4], scales_2[4]),  # M equilibrium
        ],
        dtype=float,
    )
    reaction_scales = np.concatenate((scales_1[3:], scales_2[3:]))
    return row_scales, reaction_scales


def fixed_scaled_full_matrix(
    raw_matrix: NDArray[np.complex128], point_1: RodPoint, point_2: RodPoint
) -> tuple[NDArray[np.complex128], NDArray[np.float64]]:
    row_scales, reaction_scales = fixed_physical_scaling(point_1, point_2)
    scaled = (
        np.asarray(raw_matrix, dtype=np.complex128)
        * reaction_scales[np.newaxis, :]
        / row_scales[:, np.newaxis]
    )
    return scaled, reaction_scales


def kernel_nullity(
    raw_matrix: NDArray[np.complex128], point_1: RodPoint, point_2: RodPoint
) -> tuple[int, float]:
    scaled = fixed_scaled_full_matrix(raw_matrix, point_1, point_2)[0]
    singular = np.linalg.svd(scaled, compute_uv=False)
    relative = singular / max(float(singular[0]), np.finfo(float).tiny)
    return int(np.count_nonzero(relative <= KERNEL_RELATIVE_TOLERANCE)), float(relative[-1])


def shape_block_norms(
    omega: complex,
    point_1: RodPoint,
    point_2: RodPoint,
    raw_matrix: NDArray[np.complex128],
) -> dict[str, float]:
    """Return dimensionless bending/torsion amplitudes of a full null vector."""

    scaled, reaction_scales = fixed_scaled_full_matrix(raw_matrix, point_1, point_2)
    _, singular, vh = np.linalg.svd(scaled)
    scaled_reactions = vh[-1, :].conj()
    physical_reactions = reaction_scales * scaled_reactions
    reaction_norm = float(np.linalg.norm(physical_reactions))
    matrix_norm = float(np.linalg.norm(raw_matrix, ord=2))
    null_residual = float(np.linalg.norm(raw_matrix @ physical_reactions)) / max(
        matrix_norm * reaction_norm, np.finfo(float).tiny
    )

    grid = np.linspace(0.0, 1.0, SHAPE_SAMPLES)
    bending_squared = 0.0
    torsion_squared = 0.0
    for offset, point in ((0, point_1), (3, point_2)):
        initial = cantilever_clamp_matrix(point, CLAMP, scaled=False) @ physical_reactions[
            offset : offset + 3
        ]
        physical = physical_state_trajectory_from_initial(omega, point, initial, grid)
        dimensionless = physical / _state_scales(point)[np.newaxis, :]
        length_weight = point.geometry.length / REFERENCE_LENGTH_M
        bending_density = np.sum(np.abs(dimensionless[:, BENDING_STATE_INDICES]) ** 2, axis=1)
        torsion_density = np.sum(np.abs(dimensionless[:, TORSION_STATE_INDICES]) ** 2, axis=1)
        bending_squared += length_weight * float(np.trapezoid(bending_density, grid))
        torsion_squared += length_weight * float(np.trapezoid(torsion_density, grid))

    bending_norm = math.sqrt(max(bending_squared, 0.0))
    torsion_norm = math.sqrt(max(torsion_squared, 0.0))
    full_norm = math.hypot(bending_norm, torsion_norm)
    if not math.isfinite(full_norm) or full_norm <= 0.0:
        raise RuntimeError("invalid dimensionless mode-shape amplitude norm")
    return {
        "bending_block_relative_norm": bending_norm / full_norm,
        "torsion_block_relative_norm": torsion_norm / full_norm,
        "full_null_vector_relative_residual": null_residual,
        "full_scaled_smallest_singular_ratio": float(singular[-1] / singular[0]),
    }


def reference_rows(case: CasePreset) -> list[dict[str, str]]:
    if not REFERENCE_SPECTRUM_CSV.exists():
        raise FileNotFoundError(f"missing reference spectrum: {REFERENCE_SPECTRUM_CSV}")
    with REFERENCE_SPECTRUM_CSV.open("r", encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    selected = [
        row
        for row in rows
        if row.get("model_id") == "T0"
        and math.isclose(float(row["mu"]), case.mu, rel_tol=0.0, abs_tol=1.0e-14)
        and math.isclose(float(row["tau"]), case.tau, rel_tol=0.0, abs_tol=1.0e-14)
        and math.isclose(float(row["beta_deg"]), 0.0, rel_tol=0.0, abs_tol=1.0e-14)
        and math.isclose(float(row["theta_deg"]), 0.0, rel_tol=0.0, abs_tol=1.0e-14)
    ]
    selected.sort(key=lambda row: int(row["mode"]))
    if len(selected) < FULL_ROOT_COUNT:
        raise RuntimeError(f"reference T0 spectrum is incomplete for case {case.case_id}")
    return selected[:FULL_ROOT_COUNT]


def compute_case(case: CasePreset) -> dict[str, Any]:
    point_1, point_2 = case_points(case)
    if point_1.properties.Sbar16 != 0.0 or point_2.properties.Sbar16 != 0.0:
        raise RuntimeError("theta=0 must give Sbar16 exactly equal to zero")
    factories = matrix_factories(point_1, point_2)
    spectra: dict[str, SpectrumResult] = {}
    for family, root_count in (
        ("full", FULL_ROOT_COUNT),
        ("bending", PARTIAL_ROOT_COUNT),
        ("torsion", PARTIAL_ROOT_COUNT),
    ):
        raw_factory, scaled_factory = factories[family]
        spectra[family] = solve_spectrum(
            point_1, raw_factory, scaled_factory, root_count
        )

    merged = merge_partial_spectra(spectra["bending"], spectra["torsion"])
    relevant_union = merged[:FULL_ROOT_COUNT]
    coincidences = cross_family_near_coincidences(relevant_union)
    minimum_gap_hz, minimum_gap_relative, minimum_gap_pair = minimum_cross_family_gap(
        relevant_union
    )
    if coincidences:
        details = ", ".join(
            f"indices {left + 1}/{right + 1}: {relative:.3e}"
            for left, right, relative in coincidences
        )
        raise RuntimeError(
            "a cross-family near-coincidence requires explicit multiplicity handling: "
            + details
        )

    full = spectra["full"]
    union = merged[:FULL_ROOT_COUNT]
    full_frequencies = np.asarray([root.frequency_hz for root in full.roots])
    union_frequencies = np.asarray([item.root.frequency_hz for item in union])
    union_relative = np.abs(full_frequencies - union_frequencies) / np.maximum(
        np.abs(full_frequencies), np.abs(union_frequencies)
    )
    max_union_relative = float(np.max(union_relative))
    if max_union_relative > SPECTRAL_UNION_RELATIVE_TOLERANCE:
        raise RuntimeError(
            f"spectral-union gate failed for case {case.case_id}: {max_union_relative:.3e}"
        )

    audit_rows = [
        block_separation_row(
            case,
            500.0,
            factories["full"][0](complex(2.0 * np.pi * 500.0)),
        )
    ]
    comparison_rows: list[dict[str, Any]] = []
    classification_rows: list[dict[str, Any]] = []
    max_forbidden = 0.0
    for position, (full_root, full_quality, partial, relative) in enumerate(
        zip(full.roots, full.quality, union, union_relative), start=1
    ):
        omega = full_root.omega
        raw_full = factories["full"][0](omega)
        audit_rows.append(block_separation_row(case, full_root.frequency_hz, raw_full))
        nullity, kernel_relative = kernel_nullity(raw_full, point_1, point_2)
        shape = shape_block_norms(omega, point_1, point_2, raw_full)
        forbidden = (
            shape["torsion_block_relative_norm"]
            if partial.family == "bending"
            else shape["bending_block_relative_norm"]
        )
        max_forbidden = max(max_forbidden, forbidden)
        if forbidden > COMPONENT_SUPPRESSION_TOLERANCE:
            raise RuntimeError(
                f"component-suppression gate failed for case {case.case_id}, "
                f"position {position}: {forbidden:.3e}"
            )
        lambda_full = float(lambda_reference(float(np.real(omega))))
        lambda_partial = float(
            lambda_reference(float(np.real(partial.root.omega)))
        )
        common = {
            "case_id": case.case_id,
            "mu": case.mu,
            "tau": case.tau,
            "beta_deg": BETA_DEG,
            "theta_deg": THETA_DEG,
            "sorted_position": position,
            "Lambda_full": lambda_full,
            "frequency_full_hz": full_root.frequency_hz,
            "partial_family": partial.family,
            "partial_mode_number": partial.partial_mode_number,
            "Lambda_partial": lambda_partial,
            "frequency_partial_hz": partial.root.frequency_hz,
            "relative_full_partial_difference": float(relative),
            "root_quality_status": "PASS",
            "full_quality_basis": full_quality["quality_basis"],
            "partial_quality_basis": partial.quality["quality_basis"],
            "full_accepted_determinant_residual": full_quality[
                "accepted_determinant_residual"
            ],
            "full_accepted_singular_residual": full_quality[
                "accepted_relative_singular_residual"
            ],
            "partial_accepted_determinant_residual": partial.quality[
                "accepted_determinant_residual"
            ],
            "partial_accepted_singular_residual": partial.quality[
                "accepted_relative_singular_residual"
            ],
            "kernel_nullity_at_relative_1e-8": nullity,
            "kernel_smallest_relative_singular_value": kernel_relative,
            **shape,
            "forbidden_component_relative_norm": forbidden,
            "component_suppression_status": "PASS",
        }
        comparison_rows.append(common)
        if position <= PLOTTED_ROOT_COUNT:
            classification_rows.append(dict(common))

    if any(row["status"] != "PASS" for row in audit_rows):
        raise RuntimeError(f"block-separation gate failed for case {case.case_id}")

    saved = reference_rows(case)
    reference_frequency = np.asarray([float(row["frequency_hz"]) for row in saved])
    reference_lambda = np.asarray([float(row["lambda_ref"]) for row in saved])
    computed_lambda = np.asarray(
        [float(lambda_reference(float(np.real(root.omega)))) for root in full.roots]
    )
    reference_frequency_relative = np.abs(full_frequencies - reference_frequency) / np.maximum(
        np.abs(full_frequencies), np.abs(reference_frequency)
    )
    reference_lambda_relative = np.abs(computed_lambda - reference_lambda) / np.maximum(
        np.abs(computed_lambda), np.abs(reference_lambda)
    )
    max_reference_relative = float(
        max(np.max(reference_frequency_relative), np.max(reference_lambda_relative))
    )
    if max_reference_relative > REFERENCE_REPEATABILITY_TOLERANCE:
        raise RuntimeError(
            f"saved T0 repeatability gate failed for case {case.case_id}: "
            f"{max_reference_relative:.3e}"
        )

    all_quality = [quality for spectrum in spectra.values() for quality in spectrum.quality]
    return {
        "case": case,
        "points": (point_1, point_2),
        "spectra": spectra,
        "classification_rows": classification_rows,
        "comparison_rows": comparison_rows,
        "audit_rows": audit_rows,
        "diagnostics": {
            "max_off_block_absolute": max(row["off_block_max_abs"] for row in audit_rows),
            "max_off_block_relative": max(
                row["off_block_relative_frobenius_norm"] for row in audit_rows
            ),
            "max_spectral_union_relative_difference": max_union_relative,
            "max_accepted_determinant_residual": max(
                item["accepted_determinant_residual"] for item in all_quality
            ),
            "max_accepted_singular_residual": max(
                item["accepted_relative_singular_residual"] for item in all_quality
            ),
            "max_forbidden_component_relative_norm": max_forbidden,
            "max_saved_T0_repeatability_relative_difference": max_reference_relative,
            "cross_family_near_coincidence_count": len(coincidences),
            "minimum_cross_family_gap_hz": minimum_gap_hz,
            "minimum_cross_family_relative_gap": minimum_gap_relative,
            "minimum_cross_family_gap_pair": minimum_gap_pair,
        },
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def git_context() -> dict[str, str]:
    def run(*args: str) -> str:
        completed = subprocess.run(
            ["git", *args],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
        )
        return completed.stdout.strip()

    return {
        "root": run("rev-parse", "--show-toplevel"),
        "branch": run("branch", "--show-current"),
        "head": run("rev-parse", "HEAD"),
        "status_short": run("status", "--short"),
    }


def geometry_manifest(point: RodPoint) -> dict[str, float]:
    geometry = point.geometry
    return {
        "length_m": geometry.length,
        "a_m": geometry.a,
        "b_m": geometry.b,
        "area_m2": geometry.area,
        "I_y_m4": geometry.I_y,
        "I_p_m4": geometry.I_p,
        "shear_factor": geometry.shear_factor,
    }


def max_gate_values(case_results: Sequence[Mapping[str, Any]]) -> dict[str, float]:
    diagnostics = [result["diagnostics"] for result in case_results]
    keys = (
        "max_off_block_absolute",
        "max_off_block_relative",
        "max_spectral_union_relative_difference",
        "max_accepted_determinant_residual",
        "max_accepted_singular_residual",
        "max_forbidden_component_relative_norm",
        "max_saved_T0_repeatability_relative_difference",
    )
    return {key: max(float(item[key]) for item in diagnostics) for key in keys}


def build_report(manifest: Mapping[str, Any], classification_rows: Sequence[Mapping[str, Any]]) -> str:
    cases = {str(item["case_id"]): [] for item in classification_rows}
    for row in classification_rows:
        cases[str(row["case_id"])].append(row)
    lines = [
        "# Классификация частот при beta = theta = 0",
        "",
        "Статус диагностического workflow: **PASS**.",
        "",
        "Рассмотрена модель T0: два жестко соединенных прямоугольных стержня HMS/DX-209, "
        "постановка `state_corrected`, защемление `book_slope_clamp` и обобщенное кручение. "
        "Слово «разделение» ниже означает разделение подсистем, а не разъединение стержней.",
        "",
        "## Точное разделение",
        "",
        "Для каждого стержня принят вектор состояния `[w, psi, Phi, Q, M, M_T]`. "
        "При theta = 0 величина Sbar16 равна нулю, поэтому после перестановки "
        "`[w, psi, Q, M, Phi, M_T]` матрица состояния распадается на изгибный "
        "блок `[w, psi, Q, M]` и крутильный блок `[Phi, M_T]`. При beta = 0 "
        "условия узла также разделяются: к изгибному блоку относятся "
        "`w1-w2=0`, `psi1+psi2=0`, `Q1+Q2=0`, `M1-M2=0`; к крутильному — "
        "`Phi1+Phi2=0`, `M_T1-M_T2=0`.",
        "",
        "Для полной характеристической матрицы строки переставлены как "
        f"`{list(FULL_ROW_PERMUTATION)}`, столбцы реакций — как "
        f"`{list(FULL_COLUMN_PERMUTATION)}`. Перестановка полного двухстержневого "
        f"состояния: `{list(TWO_ARM_STATE_PERMUTATION)}`.",
        "",
        "Изгибный блок сохраняет деформацию сдвига и инерцию вращения модели "
        "Тимошенко; крутильный блок сохраняет принятую обобщенную жесткость C_T. "
        "Оба блока содержат оба стержня, два внешних защемления и внутренний жесткий узел.",
        "",
        "## Первые шесть значений",
        "",
        "| k | Case A: Lambda_k | Case A: тип | Case B: Lambda_k | Case B: тип |",
        "|---:|---:|:---|---:|:---|",
    ]
    family_ru = {
        "bending": "изгибная",
        "torsion": "крутильная",
        "bending_torsion_near_coincidence": "совпадение изгибной и крутильной",
    }
    by_case_position = {
        (str(row["case_id"]), int(row["sorted_position"])): row
        for row in classification_rows
    }
    for position in range(1, PLOTTED_ROOT_COUNT + 1):
        a = by_case_position[("A", position)]
        b = by_case_position[("B", position)]
        lines.append(
            f"| {position} | {float(a['Lambda_full']):.9f} | "
            f"{family_ru[str(a['partial_family'])]} ({int(a['partial_mode_number'])}) | "
            f"{float(b['Lambda_full']):.9f} | "
            f"{family_ru[str(b['partial_family'])]} ({int(b['partial_mode_number'])}) |"
        )
    case_diagnostics = {
        str(item["case_id"]): item["diagnostics"] for item in manifest["cases"]
    }
    diagnostic_a = case_diagnostics["A"]
    diagnostic_b = case_diagnostics["B"]
    maxima = manifest["gate_results"]["raw_maxima"]
    lines.extend(
        [
            "",
            "Число в скобках — номер в соответствующем парциальном семействе.",
            "Практически точных совпадений двух семейств не обнаружено. "
            f"Минимальный межсемейный зазор равен "
            f"{diagnostic_a['minimum_cross_family_gap_hz']:.9f} Гц "
            f"(относительный зазор {diagnostic_a['minimum_cross_family_relative_gap']:.9e}) "
            f"для Case A и {diagnostic_b['minimum_cross_family_gap_hz']:.9f} Гц "
            f"({diagnostic_b['minimum_cross_family_relative_gap']:.9e}) для Case B.",
            "",
            "## Проверки",
            "",
            f"- Block separation: **PASS**; max raw off-block = "
            f"`{maxima['max_off_block_absolute']:.16e}`, max relative = "
            f"`{maxima['max_off_block_relative']:.16e}`.",
            f"- Spectral union: **PASS**; max relative residual = "
            f"`{maxima['max_spectral_union_relative_difference']:.16e}`.",
            f"- Root quality: **PASS**; max accepted determinant residual = "
            f"`{maxima['max_accepted_determinant_residual']:.16e}`, max accepted singular residual = "
            f"`{maxima['max_accepted_singular_residual']:.16e}`.",
            f"- Component suppression: **PASS**; max forbidden-block amplitude ratio = "
            f"`{maxima['max_forbidden_component_relative_norm']:.16e}`.",
            f"- Reproducibility: **PASS**; max relative difference from saved verified T0 rows = "
            f"`{maxima['max_saved_T0_repeatability_relative_difference']:.16e}`. "
            "Режим `--reuse-data` проверяет хеши научных CSV без запуска root solver.",
            "",
            "Амплитудная проверка форм является только подтверждением блочной классификации. "
            "Для каждого стержня физическое состояние делится на масштабы "
            "`[L, 1, 1, Ex*Iy/L^2, Ex*Iy/L, C_T/L]`. Затем отдельно интегрируются "
            "квадраты безразмерных компонент `[w, psi, Q, M]` и `[Phi, M_T]` "
            "по длине обоих стержней. Энергетические доли не используются.",
            "",
            "## Ограничения",
            "",
            "Результат относится только к beta = theta = 0 и только к двум геометриям "
            "(mu, tau) = (0, 0) и (0.5, -0.2). Продолжение по theta или beta, MAC и "
            "энергетическая классификация не выполнялись. Принадлежность форм при "
            "изменении параметров и наследственность форм не устанавливались.",
            "",
            "## Предлагаемый LaTeX-блок (не внесен в статью)",
            "",
            "```latex",
            r"В ортотропном прямом случае $\beta=\theta=0$ полная частотная задача точно распадается на независимые изгибную и крутильную подсистемы. При $(\mu,\tau)=(0,0)$ к изгибному семейству относятся значения с номерами $k=1,2,4,6$, а к крутильному семейству --- значения с номерами $k=3,5$. При $(\mu,\tau)=(0{,}5,-0{,}2)$ изгибными являются значения с номерами $k=1,2,3,5$, а крутильными --- значения с номерами $k=4,6$. Классификация выполнена только для двух указанных геометрий при $\beta=\theta=0$ и не означает сохранения принадлежности формы при изменении $\theta$ или $\beta$.",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


def run_computation(output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    results = [compute_case(case) for case in CASES]
    classification_rows = [row for result in results for row in result["classification_rows"]]
    comparison_rows = [row for result in results for row in result["comparison_rows"]]
    audit_rows = [row for result in results for row in result["audit_rows"]]

    classification_path = output_dir / "mode_classification.csv"
    comparison_path = output_dir / "full_partial_spectrum_comparison.csv"
    audit_path = output_dir / "block_separation_audit.csv"
    write_csv(classification_path, classification_rows)
    write_csv(comparison_path, comparison_rows)
    write_csv(audit_path, audit_rows)

    raw_maxima = max_gate_values(results)
    manifest: dict[str, Any] = {
        "workflow": "yartsev_ch2_beta0_theta0_mode_classification",
        "status": "PASS",
        "created_utc_note": "runtime timestamp intentionally omitted for reproducible metadata",
        "git": git_context(),
        "scientific_contract": {
            "model": "T0 Chapter-2 state_corrected Timoshenko plus generalized torsion",
            "material": "HMS/DX-209",
            "material_mode": MATERIAL_MODE,
            "theta_1_deg": THETA_DEG,
            "theta_2_deg": THETA_DEG,
            "beta_deg": BETA_DEG,
            "clamp": CLAMP,
            "joint": "unchanged J_book(beta)",
            "full_roots": FULL_ROOT_COUNT,
            "reported_roots": PLOTTED_ROOT_COUNT,
            "partial_roots_per_family": PARTIAL_ROOT_COUNT,
            "classification_basis": "exact partial-block spectral membership",
            "energy_classification": False,
            "MAC": False,
            "parameter_continuation": False,
        },
        "normalization": {
            "definition": "Lambda=(rho*A0*omega^2*l^4/(Ex0*Iy0))^(1/4)",
            "rho_kg_m3": hms_dx_209_material().rho,
            "a0_m": REFERENCE_A_M,
            "b0_m": REFERENCE_B_M,
            "A0_m2": REFERENCE_A_M * REFERENCE_B_M,
            "Iy0_m4": REFERENCE_A_M**3 * REFERENCE_B_M / 12.0,
            "l_m": REFERENCE_LENGTH_M,
            "Ex0_pa": REFERENCE_EX_PA,
        },
        "block_structure": {
            "rod_state_order": ["w", "psi", "Phi", "Q", "M", "M_T"],
            "bending_state": ["w", "psi", "Q", "M"],
            "torsion_state": ["Phi", "M_T"],
            "full_row_permutation_zero_based": list(FULL_ROW_PERMUTATION),
            "full_column_permutation_zero_based": list(FULL_COLUMN_PERMUTATION),
            "two_arm_state_permutation_zero_based": list(TWO_ARM_STATE_PERMUTATION),
            "bending_joint_conditions": [
                "w1-w2=0",
                "psi1+psi2=0",
                "Q1+Q2=0",
                "M1-M2=0",
            ],
            "torsion_joint_conditions": ["Phi1+Phi2=0", "M_T1-M_T2=0"],
        },
        "thresholds": {
            "root_normalized_determinant_residual": ROOT_DETERMINANT_TOLERANCE,
            "root_relative_singular_residual": ROOT_SINGULAR_TOLERANCE,
            "block_relative_frobenius": BLOCK_RELATIVE_TOLERANCE,
            "spectral_union_relative": SPECTRAL_UNION_RELATIVE_TOLERANCE,
            "component_suppression_relative": COMPONENT_SUPPRESSION_TOLERANCE,
            "saved_T0_repeatability_relative": REFERENCE_REPEATABILITY_TOLERANCE,
            "near_coincidence_relative": NEAR_COINCIDENCE_RELATIVE_TOLERANCE,
            "kernel_relative_singular": KERNEL_RELATIVE_TOLERANCE,
        },
        "cases": [],
        "gate_results": {
            "block_separation": "PASS",
            "spectral_union": "PASS",
            "root_quality": "PASS",
            "component_suppression": "PASS",
            "reproducibility_against_saved_T0": "PASS",
            "raw_maxima": raw_maxima,
        },
        "sources": {
            "reference_spectrum_csv": str(REFERENCE_SPECTRUM_CSV.relative_to(ROOT)),
            "reference_spectrum_sha256": sha256(REFERENCE_SPECTRUM_CSV),
            "article_tex_read_only": str(ARTICLE_TEX.relative_to(ROOT)),
            "article_tex_sha256": sha256(ARTICLE_TEX) if ARTICLE_TEX.exists() else None,
        },
        "outputs": {
            "mode_classification_csv": str(classification_path.relative_to(ROOT)),
            "full_partial_spectrum_comparison_csv": str(comparison_path.relative_to(ROOT)),
            "block_separation_audit_csv": str(audit_path.relative_to(ROOT)),
            "manifest_json": str((output_dir / "manifest.json").relative_to(ROOT)),
            "report_md": str((output_dir / "report.md").relative_to(ROOT)),
        },
        "scientific_csv_sha256": {
            classification_path.name: sha256(classification_path),
            comparison_path.name: sha256(comparison_path),
            audit_path.name: sha256(audit_path),
        },
        "total_runtime_seconds": time.perf_counter() - started,
    }
    for result in results:
        point_1, point_2 = result["points"]
        spectra: Mapping[str, SpectrumResult] = result["spectra"]
        manifest["cases"].append(
            {
                **asdict(result["case"]),
                "geometry_1": geometry_manifest(point_1),
                "geometry_2": geometry_manifest(point_2),
                "Sbar16_arm_1": float(np.real(point_1.properties.Sbar16)),
                "Sbar16_arm_2": float(np.real(point_2.properties.Sbar16)),
                "runtimes_seconds": {
                    name: spectrum.runtime_seconds for name, spectrum in spectra.items()
                },
                "matrix_evaluations": {
                    name: {
                        "scaled_solver": spectrum.matrix_evaluations,
                        "physical_raw_quality": spectrum.raw_quality_evaluations,
                    }
                    for name, spectrum in spectra.items()
                },
                "diagnostics": result["diagnostics"],
            }
        )

    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    report_path = output_dir / "report.md"
    report_path.write_text(build_report(manifest, classification_rows), encoding="utf-8")
    return manifest


def validate_reused_data(output_dir: Path) -> dict[str, Any]:
    manifest_path = output_dir / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"--reuse-data requires {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS":
        raise RuntimeError("saved manifest status is not PASS")
    for filename, expected_hash in manifest["scientific_csv_sha256"].items():
        path = output_dir / filename
        if not path.exists() or sha256(path) != expected_hash:
            raise RuntimeError(f"saved scientific CSV hash mismatch: {path}")
    with (output_dir / "mode_classification.csv").open(
        "r", encoding="utf-8", newline=""
    ) as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != len(CASES) * PLOTTED_ROOT_COUNT:
        raise RuntimeError("saved mode-classification row count is invalid")
    if any(row["root_quality_status"] != "PASS" for row in rows):
        raise RuntimeError("saved mode-classification quality is not PASS")
    return manifest


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = args.output_dir.resolve()
    if args.reuse_data:
        manifest = validate_reused_data(output_dir)
        print("Mode classification: PASS (reused data; scientific root solves = 0)")
    else:
        manifest = run_computation(output_dir)
        print("Mode classification: PASS")
    print(f"Output: {output_dir}")
    print(
        "Max spectral-union relative residual: "
        f"{manifest['gate_results']['raw_maxima']['max_spectral_union_relative_difference']:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
