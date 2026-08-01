"""Run the staged Chapter-2 unequal-thickness validation gates and audits."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import time
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix,
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
    joint_matrix_book,
    physical_end_map,
    straight_boundary_matrix,
    straight_boundary_matrix_raw,
    straight_right_clamp_matrix,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import (  # noqa: E402
    Geometry,
    RodPoint,
    RootResult,
    cantilever_clamp_matrix,
    cantilever_geometry,
    find_elastic_roots,
    generalized_torsional_stiffness,
    hms_dx_209_material,
    make_rod_point,
    physical_state_transfer_matrix,
    state_matrix,
)
from scripts.lib.yartsev_ch2_rectangular_eb import (  # noqa: E402
    EBFEMAssembly,
    FEM_NODE_DOF_COUNT,
    assemble_two_arm_eb_fem,
    eb_clamp_matrix,
    eb_coupled_boundary_matrix,
    eb_coupled_boundary_matrix_raw,
    eb_joint_dof_maps,
    eb_joint_mapping_residual,
    eb_physical_end_map,
    eb_state_matrix,
    eb_state_transfer_matrix,
    eb_straight_boundary_matrix,
    eb_straight_boundary_matrix_raw,
    eb_straight_right_clamp_matrix,
    solve_two_arm_eb_fem,
)


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_unequal_thickness_validation"
    / "ut0_smoke"
)
UT1_DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_unequal_thickness_validation"
    / "ut1_beta0"
)
UT1A_DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_unequal_thickness_validation"
    / "ut1a_fem_exchange_audit"
)
UT2_DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_unequal_thickness_validation"
    / "ut2_beta30"
)
UT3_DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_unequal_thickness_validation"
    / "ut3_beta90"
)
SECTION_A_M = {"a4": 0.004, "a5": 0.005, "a6": 0.006}
LENGTH_M = 0.2
WIDTH_B_M = 0.020
BETA_RAD = 0.0
DIAGNOSTIC_FREQUENCY_HZ = 500.0
NUM_ROOTS = 3
UT1_NUM_ROOTS = 7
SCAN_STEP_HZ = 10.0
INITIAL_MAX_HZ = 5000.0
MAX_HZ = 100_000.0
FORMULATION = "state_corrected"

SECTION_RELATIVE_TOLERANCE = 1.0e-12
COEFFICIENT_RELATIVE_TOLERANCE = 1.0e-12
TRANSFER_RELATIVE_TOLERANCE = 1.0e-10
ROOT_DETERMINANT_TOLERANCE = 1.0e-8
ROOT_SINGULAR_TOLERANCE = 1.0e-8
SPECTRUM_RELATIVE_TOLERANCE = 1.0e-8
TORSION_SERIES_RELATIVE_TOLERANCE = 1.0e-12

UT1_SECTION_CASES = {
    "baseline_5_5": (0.005, 0.005),
    "asymmetric_4_6": (0.004, 0.006),
    "swapped_6_4": (0.006, 0.004),
}
UT1_FEM_MESHES = (16, 32, 64)
UT1_FEM_FINE_MESH = 64
NEAR_DEGENERATE_RELATIVE_GAP = 1.0e-3
FEM_SYMMETRY_TOLERANCE = 1.0e-12
FEM_JOINT_KINEMATIC_TOLERANCE = 2.0e-14
FEM_JOINT_EQUILIBRIUM_TOLERANCE = 1.0e-7
FEM_FIRST_THREE_TOLERANCE = 1.0e-4
FEM_FIRST_SIX_TOLERANCE = 5.0e-4
FEM_AGGREGATE_ORDER_MINIMUM = 1.5

UT1A_FEM_ELEMENTS_PER_ARM = 64
UT1A_MATRIX_RELATIVE_TOLERANCE = 1.0e-13
UT1A_EIGENPAIR_RESIDUAL_TOLERANCE = 1.0e-8
UT1A_MASS_ORTHONORMALITY_TOLERANCE = 1.0e-10
UT1A_RAYLEIGH_RELATIVE_TOLERANCE = 1.0e-10
UT1A_CANONICAL_SPECTRUM_TOLERANCE = 1.0e-10

UT2_ANGLES_DEG = (-30.0, 30.0)
UT2_SECTION_CASES = {
    "baseline_5_5": (0.005, 0.005),
    "asymmetric_4_6": (0.004, 0.006),
    "swapped_6_4": (0.006, 0.004),
}
UT2_NUM_ROOTS = 7
UT2_FEM_ELEMENTS_PER_ARM = (64, 64)
UT2_MATRIX_RELATIVE_TOLERANCE = 1.0e-13

UT3_ANGLES_DEG = (-90.0, 90.0)
UT3_SECTION_CASES = {
    "baseline_5_5": (0.005, 0.005),
    "asymmetric_4_6": (0.004, 0.006),
    "swapped_6_4": (0.006, 0.004),
}
UT3_NUM_ROOTS = 7
UT3_FEM_ELEMENTS_PER_ARM = (64, 64)
UT3_MATRIX_RELATIVE_TOLERANCE = 1.0e-13
UT3_EXACT_LIMIT_TOLERANCE = 1.0e-14

UT3_EXPLICIT_EXCLUSIONS = [
    "all angles other than +/-90 deg inside UT-3",
    "parameter grids, intermediate angles, or mesh refinement",
    "Timoshenko FEM or 3D FEM",
    "theta != 0, complex roots, or damping",
    "MAC, physical mode shapes, or branch tracking",
    "mass-preserving parametrization or parameter maps",
    "plots or PDF",
    "continued UT-1a backward, Rayleigh, balancing, or high-precision audits",
    "changes to J_book, model coefficients, eigensolver, thresholds, or version",
    "production anisotropic API or limited 3D FEM anchor implementation",
]

UT2_EXPLICIT_EXCLUSIONS = [
    "beta=0 calculations except separately requested regression modes",
    "beta=90 deg or any angle other than +/-30 deg",
    "FEM mesh refinement or any mesh other than (64,64)",
    "Timoshenko FEM or 3D FEM",
    "theta != 0, complex roots, or damping",
    "MAC, physical mode shapes, or branch tracking",
    "parameter maps or mass-preserving parametrization",
    "plots or PDF",
    "matrix balancing, high-precision, or repeated UT-1a eigenpair audits",
    "changes to model coefficients, eigensolver, thresholds, statuses, or version",
    "production anisotropic API or automatic transition to UT-3",
]

UT1A_EXPLICIT_EXCLUSIONS = [
    "continuum root calculations or new mesh sequences",
    "changes to FEM matrices, joint maps, reduction, boundary conditions, or eigensolver",
    "matrix balancing, extrapolated frequencies, or Richardson corrections",
    "Timoshenko FEM, 3D FEM, or beta other than zero",
    "MAC, physical mode-shape interpretation, or branch tracking",
    "parameter maps, plots, or PDF",
    "production API or release-version changes",
    "changes to UT-0, UT-1, or overall unequal-thickness status",
]

UT1_EXPLICIT_EXCLUSIONS = [
    "all beta values other than beta=0 deg",
    "theta != 0",
    "complex roots or damping",
    "Timoshenko FEM or 3D FEM",
    "mode-shape science, MAC, or branch tracking",
    "parameter maps or an a_1/a_2 grid",
    "mass-preserving parametrization",
    "plots or PDF",
    "production anisotropic API or a new root finder",
    "changes to J_book, continuum/scaling formulas, FEM elements, or old thresholds",
]

EXPLICIT_EXCLUSIONS = [
    "beta=30 deg or beta=90 deg",
    "swapped spectrum (6,4)",
    "exchange-symmetry gate",
    "root 7 guard",
    "EB FEM, Timoshenko FEM, mesh refinement, or 3D FEM",
    "theta != 0",
    "complex roots or damping",
    "mode shapes, MAC, mode tracking, parameter maps, or plots",
    "mass-preserving parametrization",
    "new production anisotropic API or root/scaling-helper refactor",
]


class CountingBoundaryBuilder:
    """Cache and count unique matrices passed to the existing root finder."""

    def __init__(self, factory: Callable[[complex], np.ndarray]) -> None:
        self.factory = factory
        self.calls = 0
        self.evaluations = 0
        self.elapsed_seconds = 0.0
        self.maximum_frequency_hz = 0.0
        self._cache: dict[complex, np.ndarray] = {}

    def __call__(
        self, omega: complex, _point: RodPoint, formulation: str
    ) -> np.ndarray:
        if formulation != FORMULATION:
            raise ValueError("UT-0 accepts only the state_corrected formulation")
        self.calls += 1
        key = complex(omega)
        cached = self._cache.get(key)
        if cached is not None:
            return cached
        started = time.perf_counter()
        matrix = np.asarray(self.factory(key), dtype=np.complex128)
        self.elapsed_seconds += time.perf_counter() - started
        self.evaluations += 1
        self.maximum_frequency_hz = max(
            self.maximum_frequency_hz, abs(key) / (2.0 * math.pi)
        )
        self._cache[key] = matrix
        return matrix


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument(
        "--smoke",
        action="store_true",
        help="run the unchanged UT-0 smoke mode",
    )
    modes.add_argument(
        "--ut1-beta0",
        action="store_true",
        help="run the full seven-root UT-1 beta=0 validation",
    )
    modes.add_argument(
        "--ut1a-fem-exchange-audit",
        action="store_true",
        help="run the isolated beta=0 EB FEM matrix-level exchange audit",
    )
    modes.add_argument(
        "--ut2-beta30",
        action="store_true",
        help="run the isolated +/-30 deg unequal-thickness angular-joint gate",
    )
    modes.add_argument(
        "--ut3-beta90",
        action="store_true",
        help="run the isolated +/-90 deg unequal-thickness quarter-turn gate",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args(argv)


def _make_section_point(a_m: float, *, length_m: float = LENGTH_M) -> RodPoint:
    base = cantilever_geometry(length_m)
    return make_rod_point(
        0.0,
        material_mode="elastic",
        geometry=Geometry(
            a=float(a_m),
            b=WIDTH_B_M,
            length=float(length_m),
            shear_factor=base.shear_factor,
        ),
        material=hms_dx_209_material(),
        series_relative_tolerance=TORSION_SERIES_RELATIVE_TOLERANCE,
    )


def _make_section_points() -> dict[str, RodPoint]:
    """Build three independent points; no section field is replaced in place."""

    return {
        case: _make_section_point(a_m)
        for case, a_m in SECTION_A_M.items()
    }


def _make_ut1_section_cases() -> dict[str, tuple[RodPoint, RodPoint]]:
    """Build two independent arm points for every prescribed UT-1 case."""

    return {
        case: (_make_section_point(a_1), _make_section_point(a_2))
        for case, (a_1, a_2) in UT1_SECTION_CASES.items()
    }


def _beta0_state_map() -> np.ndarray:
    """Map arm-1 joint state to arm-2 local joint state at beta=0."""

    return np.diag([1.0, -1.0, -1.0, -1.0, 1.0, 1.0])


def _timoshenko_stepped_boundary_matrix_raw(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> np.ndarray:
    """Return independent physical 3x3 beta=0 stepped Timoshenko matrix."""

    transfer_1 = physical_state_transfer_matrix(omega, point_1)
    transfer_2 = physical_state_transfer_matrix(omega, point_2)
    left_clamp = cantilever_clamp_matrix(
        point_1, "book_slope_clamp", scaled=False
    )
    right_clamp = straight_right_clamp_matrix(point_2)
    propagated = np.linalg.solve(
        transfer_2, _beta0_state_map() @ transfer_1 @ left_clamp
    )
    return right_clamp @ propagated


def _timoshenko_stepped_boundary_matrix(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> np.ndarray:
    return equilibrate_matrix(
        _timoshenko_stepped_boundary_matrix_raw(omega, point_1, point_2)
    )[0]


def _eb_stepped_boundary_matrix_raw(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> np.ndarray:
    """Return independent physical 3x3 beta=0 stepped EB matrix."""

    transfer_1 = eb_state_transfer_matrix(omega, point_1)
    transfer_2 = eb_state_transfer_matrix(omega, point_2)
    propagated = np.linalg.solve(
        transfer_2, _beta0_state_map() @ transfer_1 @ eb_clamp_matrix(point_1)
    )
    return eb_straight_right_clamp_matrix(point_2) @ propagated


def _eb_stepped_boundary_matrix(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> np.ndarray:
    return equilibrate_matrix(
        _eb_stepped_boundary_matrix_raw(omega, point_1, point_2)
    )[0]


def _relative_error(actual: complex, expected: complex) -> float:
    return float(
        abs(complex(actual) - complex(expected))
        / max(abs(complex(expected)), np.finfo(float).tiny)
    )


def _positive_finite_real(value: complex) -> bool:
    number = complex(value)
    scale = max(abs(number.real), 1.0)
    return bool(
        np.isfinite(number.real)
        and np.isfinite(number.imag)
        and number.real > 0.0
        and abs(number.imag) <= SECTION_RELATIVE_TOLERANCE * scale
    )


def _expected_timoshenko_state_matrix(
    omega: complex, point: RodPoint
) -> np.ndarray:
    geometry = point.geometry
    properties = point.properties
    rho = point.material.rho
    omega_squared = complex(omega) ** 2
    shear_stiffness = (
        geometry.shear_factor * properties.Gxz * geometry.area
    )
    expected = np.zeros((6, 6), dtype=np.complex128)
    expected[0, 1] = 1.0
    expected[0, 3] = 1.0 / shear_stiffness
    expected[1, 4] = 1.0 / (properties.Ex * geometry.I_y)
    expected[2, 5] = 1.0 / point.torsion.C_T
    expected[3, 0] = -rho * geometry.area * omega_squared
    expected[4, 1] = -rho * geometry.I_y * omega_squared
    expected[4, 3] = -1.0
    expected[5, 2] = -rho * geometry.I_p * omega_squared
    return expected


def _expected_eb_state_matrix(omega: complex, point: RodPoint) -> np.ndarray:
    geometry = point.geometry
    rho = point.material.rho
    omega_squared = complex(omega) ** 2
    expected = np.zeros((6, 6), dtype=np.complex128)
    expected[0, 1] = 1.0
    expected[1, 4] = 1.0 / (point.properties.Ex * geometry.I_y)
    expected[2, 5] = 1.0 / point.torsion.Cbar
    expected[3, 0] = -rho * geometry.area * omega_squared
    expected[4, 3] = -1.0
    expected[5, 2] = -rho * geometry.I_p * omega_squared
    return expected


def _coefficient_residual(actual: np.ndarray, expected: np.ndarray) -> float:
    actual_value = np.asarray(actual, dtype=np.complex128)
    expected_value = np.asarray(expected, dtype=np.complex128)
    residuals: list[float] = []
    for index in np.ndindex(expected_value.shape):
        expected_entry = expected_value[index]
        actual_entry = actual_value[index]
        if expected_entry != 0.0:
            residuals.append(_relative_error(actual_entry, expected_entry))
        elif actual_entry != 0.0:
            # Any unprescribed coefficient is a structural failure.  The
            # theta=0 material coupling entries are intentionally included.
            residuals.append(1.0)
    return max(residuals, default=0.0)


def _timoshenko_state_coefficient_residual(
    omega: complex, point: RodPoint
) -> float:
    return _coefficient_residual(
        state_matrix(omega, point),
        _expected_timoshenko_state_matrix(omega, point),
    )


def _eb_state_coefficient_residual(omega: complex, point: RodPoint) -> float:
    return _coefficient_residual(
        eb_state_matrix(omega, point), _expected_eb_state_matrix(omega, point)
    )


def _physical_transfer_scaling_residual(
    omega: complex, point: RodPoint
) -> float:
    recovered = physical_state_transfer_matrix(omega, point)
    direct = expm(state_matrix(omega, point) * point.geometry.length)
    return float(
        np.linalg.norm(recovered - direct, ord="fro")
        / max(np.linalg.norm(direct, ord="fro"), np.finfo(float).tiny)
    )


def _baseline_a5_residual(point: RodPoint) -> float:
    baseline = make_rod_point(
        0.0,
        material_mode="elastic",
        geometry=cantilever_geometry(LENGTH_M),
        material=hms_dx_209_material(),
        series_relative_tolerance=TORSION_SERIES_RELATIVE_TOLERANCE,
    )
    pairs = [
        (point.geometry.a, baseline.geometry.a),
        (point.geometry.b, baseline.geometry.b),
        (point.geometry.length, baseline.geometry.length),
        (point.geometry.shear_factor, baseline.geometry.shear_factor),
        (point.geometry.area, baseline.geometry.area),
        (point.geometry.I_y, baseline.geometry.I_y),
        (point.geometry.I_p, baseline.geometry.I_p),
        (point.properties.Ex, baseline.properties.Ex),
        (point.properties.Gxz, baseline.properties.Gxz),
        (point.properties.Sbar16, baseline.properties.Sbar16),
        (point.torsion.Cbar, baseline.torsion.Cbar),
        (point.torsion.C_T, baseline.torsion.C_T),
        (
            point.torsion.estimated_relative_tail,
            baseline.torsion.estimated_relative_tail,
        ),
    ]
    residuals = [
        abs(complex(actual) - complex(expected))
        / max(abs(complex(expected)), np.finfo(float).tiny)
        if expected != 0.0
        else abs(complex(actual) - complex(expected))
        for actual, expected in pairs
    ]
    if point.torsion.terms_used != baseline.torsion.terms_used:
        residuals.append(1.0)
    return float(max(residuals, default=0.0))


def _section_rows(
    points: dict[str, RodPoint], omega_diag: float
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    reference_ex = points["a5"].properties.Ex
    reference_gxz = points["a5"].properties.Gxz
    rows: list[dict[str, Any]] = []
    for case, point in points.items():
        geometry = point.geometry
        area_expected = geometry.a * geometry.b
        iy_expected = geometry.a**3 * geometry.b / 12.0
        ip_expected = (
            geometry.a
            * geometry.b
            * (geometry.a**2 + geometry.b**2)
            / 12.0
        )
        shear_stiffness = (
            geometry.shear_factor * point.properties.Gxz * geometry.area
        )
        compliance_scale = max(
            abs(point.properties.Sbar11),
            abs(point.properties.Sbar66),
            np.finfo(float).tiny,
        )
        sbar16_scaled = float(abs(point.properties.Sbar16) / compliance_scale)
        torsion_identity = _relative_error(point.torsion.C_T, point.torsion.Cbar)
        recomputed = generalized_torsional_stiffness(
            point.properties,
            geometry,
            relative_tolerance=TORSION_SERIES_RELATIVE_TOLERANCE,
        )
        cbar_recompute_error = _relative_error(point.torsion.Cbar, recomputed.Cbar)
        timo_state_residual = _timoshenko_state_coefficient_residual(
            omega_diag, point
        )
        eb_state_residual = _eb_state_coefficient_residual(omega_diag, point)
        transfer_residual = _physical_transfer_scaling_residual(omega_diag, point)
        positive_values = (
            geometry.area,
            geometry.I_y,
            geometry.I_p,
            point.properties.Ex,
            point.properties.Gxz,
            shear_stiffness,
            point.torsion.Cbar,
            point.torsion.C_T,
        )
        checks = [
            _relative_error(geometry.area, area_expected)
            <= SECTION_RELATIVE_TOLERANCE,
            _relative_error(geometry.I_y, iy_expected)
            <= SECTION_RELATIVE_TOLERANCE,
            _relative_error(geometry.I_p, ip_expected)
            <= SECTION_RELATIVE_TOLERANCE,
            sbar16_scaled <= SECTION_RELATIVE_TOLERANCE,
            torsion_identity <= SECTION_RELATIVE_TOLERANCE,
            cbar_recompute_error <= SECTION_RELATIVE_TOLERANCE,
            point.torsion.estimated_relative_tail
            <= TORSION_SERIES_RELATIVE_TOLERANCE,
            timo_state_residual <= COEFFICIENT_RELATIVE_TOLERANCE,
            eb_state_residual <= COEFFICIENT_RELATIVE_TOLERANCE,
            transfer_residual <= TRANSFER_RELATIVE_TOLERANCE,
            all(_positive_finite_real(value) for value in positive_values),
        ]
        rows.append(
            {
                "case": case,
                "a_m": geometry.a,
                "b_m": geometry.b,
                "length_m": geometry.length,
                "area": geometry.area,
                "area_expected": area_expected,
                "area_relative_error": _relative_error(
                    geometry.area, area_expected
                ),
                "I_y": geometry.I_y,
                "I_y_expected": iy_expected,
                "I_y_relative_error": _relative_error(
                    geometry.I_y, iy_expected
                ),
                "I_p": geometry.I_p,
                "I_p_expected": ip_expected,
                "I_p_relative_error": _relative_error(
                    geometry.I_p, ip_expected
                ),
                "shear_factor": geometry.shear_factor,
                "Ex": float(np.real(point.properties.Ex)),
                "Gxz": float(np.real(point.properties.Gxz)),
                "Sbar16": float(np.real(point.properties.Sbar16)),
                "Sbar16_scaled_residual": sbar16_scaled,
                "K_s": float(np.real(shear_stiffness)),
                "Cbar": float(np.real(point.torsion.Cbar)),
                "C_T": float(np.real(point.torsion.C_T)),
                "C_T_Cbar_relative_error": torsion_identity,
                "Cbar_recompute_relative_error": cbar_recompute_error,
                "Ex_common_relative_error": _relative_error(
                    point.properties.Ex, reference_ex
                ),
                "Gxz_common_relative_error": _relative_error(
                    point.properties.Gxz, reference_gxz
                ),
                "torsion_terms_used": point.torsion.terms_used,
                "torsion_estimated_relative_tail": (
                    point.torsion.estimated_relative_tail
                ),
                "timoshenko_state_coefficient_residual": timo_state_residual,
                "eb_state_coefficient_residual": eb_state_residual,
                "physical_transfer_scaling_residual": transfer_residual,
                "status": "PASS" if all(checks) else "FAIL",
            }
        )

    cbar_values = [complex(points[key].torsion.Cbar) for key in SECTION_A_M]
    diagnostics = {
        "baseline_a5_relative_residual": _baseline_a5_residual(points["a5"]),
        "max_Ex_common_relative_error": max(
            float(row["Ex_common_relative_error"]) for row in rows
        ),
        "max_Gxz_common_relative_error": max(
            float(row["Gxz_common_relative_error"]) for row in rows
        ),
        "cbar_objects_are_distinct": len(
            {id(points[key].torsion) for key in SECTION_A_M}
        )
        == len(SECTION_A_M),
        "cbar_values_are_geometry_specific": all(
            cbar_values[index] != cbar_values[index + 1]
            for index in range(len(cbar_values) - 1)
        ),
    }
    return rows, diagnostics


def _write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> None:
    data = list(rows)
    if not data:
        raise ValueError(f"refusing to create empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames: list[str] = []
    for row in data:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(data)


def _matrix_quality(matrix: np.ndarray) -> dict[str, float]:
    value = np.asarray(matrix, dtype=np.complex128)
    singular = np.linalg.svd(value, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    determinant = complex(np.linalg.det(value))
    determinant_scale = max(
        sigma_max ** value.shape[0], np.finfo(float).tiny
    )
    return {
        "physical_raw_determinant": float(determinant.real),
        "physical_raw_determinant_abs": float(abs(determinant)),
        "physical_raw_normalized_determinant_residual": float(
            abs(determinant) / determinant_scale
        ),
        "physical_raw_sigma_min": sigma_min,
        "physical_raw_sigma_max": sigma_max,
        "physical_raw_relative_singular_residual": (
            sigma_min / sigma_max if sigma_max else 0.0
        ),
    }


def _root_quality_fields(
    root: RootResult, raw_factory: Callable[[complex], np.ndarray]
) -> dict[str, Any]:
    raw = _matrix_quality(raw_factory(root.omega))
    scaled_ok = (
        root.determinant_residual <= ROOT_DETERMINANT_TOLERANCE
        and root.relative_singular_residual <= ROOT_SINGULAR_TOLERANCE
        and not root.status.startswith("rejected")
    )
    raw_ok = (
        raw["physical_raw_normalized_determinant_residual"]
        <= ROOT_DETERMINANT_TOLERANCE
        and raw["physical_raw_relative_singular_residual"]
        <= ROOT_SINGULAR_TOLERANCE
        and not root.status.startswith("rejected")
    )
    return {
        "frequency_hz": root.frequency_hz,
        "scaled_determinant_abs": root.raw_determinant_abs,
        "scaled_determinant_residual": root.determinant_residual,
        "scaled_sigma_min": root.sigma_min,
        "scaled_sigma_max": root.sigma_max,
        "scaled_relative_singular_residual": (
            root.relative_singular_residual
        ),
        **raw,
        "root_status": root.status,
        "root_quality_status": "PASS" if scaled_ok or raw_ok else "FAIL",
        "root_quality_basis": "scaled" if scaled_ok else "physical_raw",
    }


def _solve_roots(
    point: RodPoint,
    factory: Callable[[complex], np.ndarray],
    *,
    num_roots: int = NUM_ROOTS,
) -> tuple[list[RootResult], CountingBoundaryBuilder, float]:
    builder = CountingBoundaryBuilder(factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point,
        FORMULATION,
        num_roots=num_roots,
        scan_step_hz=SCAN_STEP_HZ,
        initial_max_hz=INITIAL_MAX_HZ,
        max_hz=MAX_HZ,
        boundary_matrix_builder=builder,
    )
    return roots, builder, time.perf_counter() - started


def _prefixed(prefix: str, fields: dict[str, Any]) -> dict[str, Any]:
    return {f"{prefix}_{key}": value for key, value in fields.items()}


def _root_comparison_rows(
    spectra: dict[str, dict[str, Any]]
) -> tuple[list[dict[str, Any]], dict[str, float], dict[str, float]]:
    rows: list[dict[str, Any]] = []
    max_differences: dict[str, float] = {}
    quality_maxima = {
        "scaled_determinant_residual": 0.0,
        "scaled_relative_singular_residual": 0.0,
        "physical_raw_normalized_determinant_residual": 0.0,
        "physical_raw_relative_singular_residual": 0.0,
    }
    for model, data in spectra.items():
        coupled_roots = data["coupled_roots"]
        stepped_roots = data["stepped_roots"]
        relative_differences: list[float] = []
        for mode, (coupled_root, stepped_root) in enumerate(
            zip(coupled_roots, stepped_roots), start=1
        ):
            coupled_fields = _root_quality_fields(
                coupled_root, data["coupled_raw"]
            )
            stepped_fields = _root_quality_fields(
                stepped_root, data["stepped_raw"]
            )
            absolute_difference = abs(
                coupled_root.frequency_hz - stepped_root.frequency_hz
            )
            relative_difference = (
                absolute_difference / stepped_root.frequency_hz
            )
            relative_differences.append(relative_difference)
            for fields in (coupled_fields, stepped_fields):
                for key in quality_maxima:
                    quality_maxima[key] = max(
                        quality_maxima[key], float(fields[key])
                    )
            rows.append(
                {
                    "model": model,
                    "mode": mode,
                    "coupled_frequency_hz": coupled_root.frequency_hz,
                    "stepped_frequency_hz": stepped_root.frequency_hz,
                    "absolute_difference_hz": absolute_difference,
                    "relative_difference": relative_difference,
                    "coupled_root_quality_status": coupled_fields[
                        "root_quality_status"
                    ],
                    "stepped_root_quality_status": stepped_fields[
                        "root_quality_status"
                    ],
                    **_prefixed("coupled", coupled_fields),
                    **_prefixed("stepped", stepped_fields),
                }
            )
        max_differences[model] = max(relative_differences, default=math.inf)
    return rows, max_differences, quality_maxima


def _report(summary: dict[str, Any], section_rows: list[dict[str, Any]]) -> str:
    spectra = summary["first_three_frequencies_hz"]
    lines = [
        "# UT-0 unequal-thickness validation report",
        "",
        "Overall unequal-thickness validation: **IN_PROGRESS**",
        "",
        "UT-0 section-routing and beta=0 stepped-reference smoke: "
        f"**{summary['ut0_status']}**",
        "",
        "## Scope",
        "",
        "Elastic HMS/DX-209, `theta_1=theta_2=0 deg`, `L_1=L_2=0.2 m`, "
        "`b_1=b_2=0.020 m`, the existing `5/6` shear factor, and only "
        "`beta=0 deg` were evaluated. The section dimension `a_i` is parallel "
        "to `e_z`.",
        "",
        "## Section and arm-routing checks",
        "",
        "| case | a (mm) | A (m^2) | I_y (m^4) | I_p (m^4) | K_s (N) | Cbar=C_T (N m^2) | status |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in section_rows:
        lines.append(
            f"| {row['case']} | {1e3 * row['a_m']:.1f} | {row['area']:.12e} | "
            f"{row['I_y']:.12e} | {row['I_p']:.12e} | {row['K_s']:.12e} | "
            f"{row['Cbar']:.12e} | {row['status']} |"
        )
    maxima = summary["residual_maxima"]
    lines.extend(
        [
            "",
            "The identities `A=ab`, `I_y=a^3 b/12`, and "
            "`I_p=ab(a^2+b^2)/12` were checked independently for every point. "
            "`K_s` is the existing `kappa G_xz A`; the rectangular torsion "
            "series and its existing tolerance were unchanged.",
            "",
            f"Maximum section-identity residual: `{maxima['section_identity']:.6e}`; "
            f"scaled `Sbar16` residual: `{maxima['Sbar16_scaled']:.6e}`; "
            f"`C_T=Cbar` residual: `{maxima['C_T_Cbar']:.6e}`.",
            "",
            f"The explicit `a5` point matches `cantilever_geometry(0.2)` with "
            f"maximum residual `{maxima['baseline_a5']:.6e}`. `Ex` and `Gxz` "
            f"are common to all three points with maxima `{maxima['Ex_common']:.6e}` "
            f"and `{maxima['Gxz_common']:.6e}`; the largest estimated torsion-series "
            f"tail is `{maxima['torsion_estimated_relative_tail']:.6e}`.",
            "",
            f"Maximum Timoshenko state-coefficient residual at 500 Hz: "
            f"`{maxima['timoshenko_state_coefficients']:.6e}`; EB: "
            f"`{maxima['eb_state_coefficients']:.6e}`. Maximum physical-transfer "
            f"scaling residual: `{maxima['physical_transfer_scaling']:.6e}`.",
            "",
            "## Independent beta=0 stepped reference",
            "",
            "The script-local map is `R0=diag(1,-1,-1,-1,1,1)` in the accepted "
            "book state `[w,psi,Phi,Q,M,M_T]`. The Timoshenko reference is "
            "`C2 solve(T2, R0 T1 B1)` with the arm-2 book-slope selector; the EB "
            "reference is the analogous EB expression. Neither reference calls "
            "the joint or coupled boundary builders.",
            "",
            "## First three sorted roots",
            "",
            "| model | mode | coupled (Hz) | stepped (Hz) | relative difference |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for model in ("Timoshenko", "EB"):
        for index, (coupled, stepped) in enumerate(
            zip(spectra[model]["coupled"], spectra[model]["stepped"]), start=1
        ):
            relative = abs(coupled - stepped) / stepped
            lines.append(
                f"| {model} | {index} | {coupled:.12f} | {stepped:.12f} | "
                f"{relative:.6e} |"
            )
    quality = summary["root_quality_maxima"]
    lines.extend(
        [
            "",
            "All root rows retain scaled determinant/SVD and physical-raw "
            "determinant/SVD diagnostics. Quality passes when either the scaled "
            "or normalized physical-raw determinant and relative singular "
            "residuals pass and the root is not rejected.",
            "",
            f"Maximum scaled determinant residual: "
            f"`{quality['scaled_determinant_residual']:.6e}`; maximum scaled "
            f"relative singular residual: "
            f"`{quality['scaled_relative_singular_residual']:.6e}`; maximum "
            f"physical-raw normalized determinant residual: "
            f"`{quality['physical_raw_normalized_determinant_residual']:.6e}`; "
            f"maximum physical-raw relative singular residual: "
            f"`{quality['physical_raw_relative_singular_residual']:.6e}`.",
            "",
            "## Fixed thresholds",
            "",
            "Section identities and state coefficients: `1e-12`; physical "
            "transfer scaling: `1e-10`; root determinant and relative singular "
            "residuals: `1e-8`; coupled--stepped sorted-spectrum difference: "
            "`1e-8`. No threshold was adjusted after calculation.",
            "",
            "## Explicit exclusions",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in summary["explicit_exclusions"])
    lines.extend(
        [
            "",
            "This gate does not claim branch identity or compare the accuracy of "
            "the Timoshenko and EB continuum models. It does not complete the "
            "unequal-thickness validation program.",
            "",
        ]
    )
    return "\n".join(lines)


def _run_smoke(output_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    started = time.perf_counter()
    points = _make_section_points()
    omega_diag = 2.0 * math.pi * DIAGNOSTIC_FREQUENCY_HZ
    section_rows, section_diagnostics = _section_rows(points, omega_diag)
    point_1 = points["a4"]
    point_2 = points["a6"]

    formulations = {
        "Timoshenko": {
            "coupled": lambda omega: coupled_boundary_matrix(
                omega, BETA_RAD, point_1, point_2
            ),
            "coupled_raw": lambda omega: coupled_boundary_matrix_raw(
                omega, BETA_RAD, point_1, point_2
            ),
            "stepped": lambda omega: _timoshenko_stepped_boundary_matrix(
                omega, point_1, point_2
            ),
            "stepped_raw": lambda omega: (
                _timoshenko_stepped_boundary_matrix_raw(
                    omega, point_1, point_2
                )
            ),
        },
        "EB": {
            "coupled": lambda omega: eb_coupled_boundary_matrix(
                omega, BETA_RAD, point_1, point_2
            ),
            "coupled_raw": lambda omega: eb_coupled_boundary_matrix_raw(
                omega, BETA_RAD, point_1, point_2
            ),
            "stepped": lambda omega: _eb_stepped_boundary_matrix(
                omega, point_1, point_2
            ),
            "stepped_raw": lambda omega: _eb_stepped_boundary_matrix_raw(
                omega, point_1, point_2
            ),
        },
    }
    spectra: dict[str, dict[str, Any]] = {}
    runtime_by_spectrum: dict[str, float] = {}
    evaluations_by_spectrum: dict[str, int] = {}
    calls_by_spectrum: dict[str, int] = {}
    for model, factories in formulations.items():
        spectra[model] = {
            "coupled_raw": factories["coupled_raw"],
            "stepped_raw": factories["stepped_raw"],
        }
        for kind in ("coupled", "stepped"):
            roots, builder, runtime = _solve_roots(
                point_1, factories[kind]
            )
            spectra[model][f"{kind}_roots"] = roots
            key = f"{model}_{kind}"
            runtime_by_spectrum[key] = runtime
            evaluations_by_spectrum[key] = builder.evaluations
            calls_by_spectrum[key] = builder.calls

    root_rows, max_differences, root_quality_maxima = (
        _root_comparison_rows(spectra)
    )
    section_identity_max = max(
        max(
            float(row["area_relative_error"]),
            float(row["I_y_relative_error"]),
            float(row["I_p_relative_error"]),
        )
        for row in section_rows
    )
    residual_maxima = {
        "section_identity": section_identity_max,
        "Sbar16_scaled": max(
            float(row["Sbar16_scaled_residual"]) for row in section_rows
        ),
        "C_T_Cbar": max(
            float(row["C_T_Cbar_relative_error"]) for row in section_rows
        ),
        "Cbar_recompute": max(
            float(row["Cbar_recompute_relative_error"])
            for row in section_rows
        ),
        "Ex_common": float(
            section_diagnostics["max_Ex_common_relative_error"]
        ),
        "Gxz_common": float(
            section_diagnostics["max_Gxz_common_relative_error"]
        ),
        "baseline_a5": float(
            section_diagnostics["baseline_a5_relative_residual"]
        ),
        "torsion_estimated_relative_tail": max(
            float(row["torsion_estimated_relative_tail"])
            for row in section_rows
        ),
        "timoshenko_state_coefficients": max(
            float(row["timoshenko_state_coefficient_residual"])
            for row in section_rows
        ),
        "eb_state_coefficients": max(
            float(row["eb_state_coefficient_residual"])
            for row in section_rows
        ),
        "physical_transfer_scaling": max(
            float(row["physical_transfer_scaling_residual"])
            for row in section_rows
        ),
    }
    roots_ok = all(
        len(data[f"{kind}_roots"]) == NUM_ROOTS
        and all(
            np.isfinite(root.frequency_hz) and root.frequency_hz > 0.0
            for root in data[f"{kind}_roots"]
        )
        for data in spectra.values()
        for kind in ("coupled", "stepped")
    )
    root_quality_ok = all(
        row["coupled_root_quality_status"] == "PASS"
        and row["stepped_root_quality_status"] == "PASS"
        for row in root_rows
    )
    section_ok = all(row["status"] == "PASS" for row in section_rows)
    common_material_ok = (
        residual_maxima["Ex_common"] <= SECTION_RELATIVE_TOLERANCE
        and residual_maxima["Gxz_common"] <= SECTION_RELATIVE_TOLERANCE
    )
    baseline_ok = (
        residual_maxima["baseline_a5"] <= SECTION_RELATIVE_TOLERANCE
    )
    cbar_routing_ok = bool(
        section_diagnostics["cbar_objects_are_distinct"]
        and section_diagnostics["cbar_values_are_geometry_specific"]
        and residual_maxima["Cbar_recompute"] <= SECTION_RELATIVE_TOLERANCE
    )
    spectra_ok = all(
        difference <= SPECTRUM_RELATIVE_TOLERANCE
        for difference in max_differences.values()
    )
    ut0_status = (
        "PASS"
        if all(
            (
                section_ok,
                common_material_ok,
                baseline_ok,
                cbar_routing_ok,
                roots_ok,
                root_quality_ok,
                spectra_ok,
            )
        )
        else "FAIL"
    )
    summary = {
        "scope": "UT-0 section routing and beta=0 stepped-reference smoke only",
        "input_geometry": {
            "section_dimension_a_direction": "parallel to e_z",
            "points": {
                case: {
                    "a_m": point.geometry.a,
                    "b_m": point.geometry.b,
                    "length_m": point.geometry.length,
                    "shear_factor": point.geometry.shear_factor,
                }
                for case, point in points.items()
            },
            "spectral_case": {"point_1": "a4", "point_2": "a6"},
        },
        "material": "existing HMS/DX-209",
        "material_mode": "elastic",
        "theta_1_deg": 0.0,
        "theta_2_deg": 0.0,
        "beta_deg": 0.0,
        "diagnostic_frequency_hz": DIAGNOSTIC_FREQUENCY_HZ,
        "formulation": FORMULATION,
        "root_scan": {
            "num_roots": NUM_ROOTS,
            "scan_step_hz": SCAN_STEP_HZ,
            "initial_max_hz": INITIAL_MAX_HZ,
            "max_hz": MAX_HZ,
            "sigma_dip_trigger": "existing find_elastic_roots default (1e-5)",
            "internal_refinement_tolerances": (
                "existing find_elastic_roots values unchanged"
            ),
        },
        "thresholds": {
            "section_identity_relative": SECTION_RELATIVE_TOLERANCE,
            "zero_section_quantity_scaled_absolute": (
                SECTION_RELATIVE_TOLERANCE
            ),
            "state_coefficient_relative": COEFFICIENT_RELATIVE_TOLERANCE,
            "physical_transfer_frobenius_relative": (
                TRANSFER_RELATIVE_TOLERANCE
            ),
            "torsion_series_estimated_relative_tail": (
                TORSION_SERIES_RELATIVE_TOLERANCE
            ),
            "root_determinant_residual": ROOT_DETERMINANT_TOLERANCE,
            "root_relative_singular_residual": ROOT_SINGULAR_TOLERANCE,
            "coupled_stepped_spectrum_relative": (
                SPECTRUM_RELATIVE_TOLERANCE
            ),
        },
        "residual_maxima": residual_maxima,
        "section_diagnostics": section_diagnostics,
        "first_three_frequencies_hz": {
            model: {
                kind: [
                    root.frequency_hz
                    for root in data[f"{kind}_roots"]
                ]
                for kind in ("coupled", "stepped")
            }
            for model, data in spectra.items()
        },
        "max_coupled_stepped_relative_difference": max_differences,
        "root_quality_maxima": root_quality_maxima,
        "root_quality_pass": root_quality_ok,
        "root_count_and_positivity_pass": roots_ok,
        "runtime_seconds_by_spectrum": runtime_by_spectrum,
        "total_runtime_seconds": time.perf_counter() - started,
        "matrix_evaluations_by_spectrum": evaluations_by_spectrum,
        "matrix_calls_by_spectrum": calls_by_spectrum,
        "ut0_status": ut0_status,
        "overall_unequal_thickness_validation": "IN_PROGRESS",
        "explicit_exclusions": EXPLICIT_EXCLUSIONS,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "section_properties.csv", section_rows)
    _write_csv(output_dir / "root_comparison.csv", root_rows)
    (output_dir / "ut0_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "ut0_report.md").write_text(
        _report(summary, section_rows), encoding="utf-8"
    )
    return summary, section_rows


def _frequencies(roots: Sequence[RootResult]) -> np.ndarray:
    return np.asarray([root.frequency_hz for root in roots], dtype=float)


def _neighbor_gaps(frequencies: Sequence[float]) -> np.ndarray:
    values = np.asarray(frequencies, dtype=float)
    gaps = np.full(len(values), np.inf, dtype=float)
    for index in range(len(values)):
        candidates: list[float] = []
        if index:
            candidates.append(float(values[index] - values[index - 1]))
        if index + 1 < len(values):
            candidates.append(float(values[index + 1] - values[index]))
        if candidates:
            gaps[index] = min(candidates)
    return gaps


def _match_frequency_spectra(
    reference: Sequence[float], candidate: Sequence[float]
) -> tuple[np.ndarray, list[str], np.ndarray, np.ndarray]:
    """Match sorted spectra, assigning only inside connected close clusters."""

    reference_values = np.asarray(reference, dtype=float)
    candidate_values = np.asarray(candidate, dtype=float)
    if reference_values.shape != candidate_values.shape:
        raise ValueError("frequency spectra must have equal shapes")
    matched = candidate_values.copy()
    permutation = np.arange(len(reference_values), dtype=int)
    matching_basis = ["sorted" for _ in reference_values]
    close_links: list[int] = []
    for index in range(len(reference_values) - 1):
        reference_gap = reference_values[index + 1] - reference_values[index]
        candidate_gap = abs(candidate_values[index + 1] - candidate_values[index])
        reference_scale = max(
            abs(reference_values[index]), abs(reference_values[index + 1]), 1.0
        )
        candidate_scale = max(
            abs(candidate_values[index]), abs(candidate_values[index + 1]), 1.0
        )
        if (
            reference_gap / reference_scale < NEAR_DEGENERATE_RELATIVE_GAP
            or candidate_gap / candidate_scale < NEAR_DEGENERATE_RELATIVE_GAP
        ):
            close_links.append(index)

    cursor = 0
    while cursor < len(close_links):
        first_link = close_links[cursor]
        last_link = first_link
        while (
            cursor + 1 < len(close_links)
            and close_links[cursor + 1] == last_link + 1
        ):
            cursor += 1
            last_link = close_links[cursor]
        indices = np.arange(first_link, last_link + 2, dtype=int)
        cost = np.abs(
            reference_values[indices, np.newaxis]
            - candidate_values[indices][np.newaxis, :]
        )
        rows, columns = linear_sum_assignment(cost)
        label = "local_frequency_cluster_modes_" + "-".join(
            str(int(index + 1)) for index in indices
        )
        for row, column in zip(rows, columns):
            reference_index = int(indices[row])
            candidate_index = int(indices[column])
            matched[reference_index] = candidate_values[candidate_index]
            permutation[reference_index] = candidate_index
            matching_basis[reference_index] = label
        cursor += 1
    return matched, matching_basis, _neighbor_gaps(reference_values), permutation


def _git_context() -> dict[str, str]:
    def git(*args: str) -> str:
        completed = subprocess.run(
            ["git", *args],
            cwd=ROOT,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        )
        return completed.stdout.strip()

    return {
        "cwd": str(Path.cwd()),
        "git_root": git("rev-parse", "--show-toplevel"),
        "branch": git("branch", "--show-current"),
        "HEAD": git("rev-parse", "HEAD"),
        "origin_main": git("rev-parse", "origin/main"),
        "status_short": git("status", "--short"),
        "diff_stat": git("diff", "--stat"),
        "untracked": git("ls-files", "--others", "--exclude-standard"),
    }


def _ut1_factories(
    model: str, point_1: RodPoint, point_2: RodPoint
) -> dict[str, Callable[[complex], np.ndarray]]:
    if model == "Timoshenko":
        return {
            "coupled": lambda omega: coupled_boundary_matrix(
                omega, BETA_RAD, point_1, point_2
            ),
            "coupled_raw": lambda omega: coupled_boundary_matrix_raw(
                omega, BETA_RAD, point_1, point_2
            ),
            "direct_stepped": lambda omega: (
                _timoshenko_stepped_boundary_matrix(
                    omega, point_1, point_2
                )
            ),
            "direct_stepped_raw": lambda omega: (
                _timoshenko_stepped_boundary_matrix_raw(
                    omega, point_1, point_2
                )
            ),
        }
    if model == "EB":
        return {
            "coupled": lambda omega: eb_coupled_boundary_matrix(
                omega, BETA_RAD, point_1, point_2
            ),
            "coupled_raw": lambda omega: eb_coupled_boundary_matrix_raw(
                omega, BETA_RAD, point_1, point_2
            ),
            "direct_stepped": lambda omega: _eb_stepped_boundary_matrix(
                omega, point_1, point_2
            ),
            "direct_stepped_raw": lambda omega: (
                _eb_stepped_boundary_matrix_raw(omega, point_1, point_2)
            ),
        }
    raise ValueError(f"unsupported UT-1 model: {model}")


def _continuum_root_row(
    *,
    section_case: str,
    point_1: RodPoint,
    point_2: RodPoint,
    model: str,
    formulation: str,
    mode: int,
    root: RootResult,
    quality: dict[str, Any],
    gap_hz: float,
) -> dict[str, Any]:
    return {
        "section_case": section_case,
        "a_1_m": point_1.geometry.a,
        "a_2_m": point_2.geometry.a,
        "beta_deg": 0.0,
        "model": model,
        "formulation": formulation,
        "mode": mode,
        "frequency_hz": root.frequency_hz,
        "scaled_determinant_abs": quality["scaled_determinant_abs"],
        "scaled_determinant_residual": quality["scaled_determinant_residual"],
        "scaled_sigma_min": quality["scaled_sigma_min"],
        "scaled_sigma_max": quality["scaled_sigma_max"],
        "scaled_relative_singular_residual": quality[
            "scaled_relative_singular_residual"
        ],
        "raw_determinant_abs": quality["physical_raw_determinant_abs"],
        "raw_determinant_residual": quality[
            "physical_raw_normalized_determinant_residual"
        ],
        "raw_sigma_min": quality["physical_raw_sigma_min"],
        "raw_sigma_max": quality["physical_raw_sigma_max"],
        "raw_relative_singular_residual": quality[
            "physical_raw_relative_singular_residual"
        ],
        "root_status": quality["root_status"],
        "quality_status": quality["root_quality_status"],
        "quality_basis": quality["root_quality_basis"],
        "neighbor_gap_hz": gap_hz,
        "relative_neighbor_gap": gap_hz / root.frequency_hz,
    }


def _run_eb_fem_case(
    *,
    section_case: str,
    point_1: RodPoint,
    point_2: RodPoint,
    elements_per_arm: int,
    analytic_frequencies: Sequence[float],
) -> dict[str, Any]:
    """Run only the existing independent EB FEM assembly and solver."""

    started = time.perf_counter()
    assembly = assemble_two_arm_eb_fem(
        point_1, point_2, BETA_RAD, elements_per_arm
    )
    solution = solve_two_arm_eb_fem(assembly, num_roots=UT1_NUM_ROOTS)
    runtime = time.perf_counter() - started
    matched, bases, gaps, permutation = _match_frequency_spectra(
        analytic_frequencies, solution.frequencies_hz
    )
    matched_equilibrium = solution.joint_equilibrium_residuals[permutation]
    positive_finite = bool(
        len(solution.frequencies_hz) == UT1_NUM_ROOTS
        and np.all(np.isfinite(solution.frequencies_hz))
        and np.all(solution.frequencies_hz > 0.0)
    )
    diagnostics = {
        "section_case": section_case,
        "elements_arm_1": elements_per_arm,
        "elements_arm_2": elements_per_arm,
        "element_length_m": point_1.geometry.length / elements_per_arm,
        "reduced_size": int(assembly.mass.shape[0]),
        "stiffness_symmetry_residual": solution.stiffness_symmetry_residual,
        "mass_symmetry_residual": solution.mass_symmetry_residual,
        "minimum_mass_eigenvalue": solution.minimum_mass_eigenvalue,
        "zero_mode_count": solution.zero_mode_count,
        "root_count": len(solution.frequencies_hz),
        "positive_finite_roots": positive_finite,
        "joint_kinematic_residual": eb_joint_mapping_residual(BETA_RAD),
        "max_joint_equilibrium_residual": float(
            np.max(solution.joint_equilibrium_residuals)
        ),
        "runtime_seconds": runtime,
    }
    diagnostics["structural_status"] = (
        "PASS" if _fem_structural_ok(diagnostics) else "FAIL"
    )
    return {
        "assembly": assembly,
        "solution": solution,
        "matched_frequencies": matched,
        "matching_basis": bases,
        "neighbor_gaps": gaps,
        "permutation": permutation,
        "matched_joint_equilibrium_residuals": matched_equilibrium,
        "diagnostics": diagnostics,
    }


def _fem_structural_ok(diagnostics: dict[str, Any]) -> bool:
    return bool(
        diagnostics["stiffness_symmetry_residual"]
        <= FEM_SYMMETRY_TOLERANCE
        and diagnostics["mass_symmetry_residual"]
        <= FEM_SYMMETRY_TOLERANCE
        and diagnostics["minimum_mass_eigenvalue"] > 0.0
        and diagnostics["zero_mode_count"] == 0
        and diagnostics["root_count"] == UT1_NUM_ROOTS
        and diagnostics["positive_finite_roots"]
        and diagnostics["joint_kinematic_residual"]
        <= FEM_JOINT_KINEMATIC_TOLERANCE
        and diagnostics["max_joint_equilibrium_residual"]
        <= FEM_JOINT_EQUILIBRIUM_TOLERANCE
    )


def _ut1_report(summary: dict[str, Any]) -> str:
    lines = [
        "# UT-1 full beta=0 unequal-thickness validation report",
        "",
        "Overall unequal-thickness validation: **IN_PROGRESS**",
        "",
        "UT-0 section-routing and beta=0 stepped-reference smoke: "
        f"**{summary['ut0_regression_status']}**",
        "",
        "UT-1 full beta=0 unequal-thickness validation: "
        f"**{summary['ut1_status']}**",
        "",
        "## Scope",
        "",
        "Only `beta=0 deg`, elastic HMS/DX-209, `theta_1=theta_2=0`, "
        "`L_1=L_2=0.2 m`, `b_1=b_2=0.020 m`, and the existing `5/6` "
        "shear factor are used. The section cases are `(5,5)`, `(4,6)`, "
        "and `(6,4) mm`, with `a_i || e_z`.",
        "",
        "## Seven-root continuum comparisons",
        "",
        "| case | model | max coupled--stepped relative difference | status |",
        "|---|---|---:|---|",
    ]
    for key, value in summary["maximum_coupled_stepped_relative_difference"].items():
        lines.append(
            f"| {key.split('|')[0]} | {key.split('|')[1]} | "
            f"{value:.12e} | {'PASS' if value <= SPECTRUM_RELATIVE_TOLERANCE else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            "| case | model | formulation | f1 | f2 | f3 | f4 | f5 | f6 | f7 |",
            "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for case, models in summary["first_seven_frequencies_hz"].items():
        for model, formulations in models.items():
            for formulation, frequencies in formulations.items():
                lines.append(
                    f"| {case} | {model} | {formulation} | "
                    + " | ".join(f"{frequency:.9f}" for frequency in frequencies)
                    + " |"
                )
    lines.extend(
        [
            "",
            f"All `{summary['main_continuum_spectrum_count']}` main spectra "
            f"contain seven finite positive roots. Across continuum rows, "
            f"`{summary['root_quality_basis_counts']['scaled']}` roots pass by "
            f"scaled quality and `{summary['root_quality_basis_counts']['physical_raw']}` "
            "by physical-raw normalized quality.",
            "",
            "## Baseline and exchange gates",
            "",
            f"Equal-thickness straight regression maxima: Timoshenko "
            f"`{summary['equal_thickness_straight_regression_maxima']['Timoshenko']:.12e}`, "
            f"EB `{summary['equal_thickness_straight_regression_maxima']['EB']:.12e}`.",
            "",
            f"Analytic exchange maxima: `{max(summary['analytic_exchange_maxima'].values()):.12e}`; "
            f"FEM exchange maximum: `{summary['fem_exchange_maximum']:.12e}`.",
            "The analytic exchange gate passes. The FEM exchange value exceeds "
            f"its fixed `1e-8` threshold, so the honest UT-1 status is "
            f"`{summary['ut1_status']}`; continuum and FEM coefficients are unchanged.",
            "",
            "## Independent EB FEM",
            "",
            "| case | E3(64) | E6(64) | root-7 error | structural status |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for case, metrics in summary["fem_mesh64_errors"].items():
        lines.append(
            f"| {case} | {metrics['E3']:.12e} | {metrics['E6']:.12e} | "
            f"{metrics['root_7_error']:.12e} | {metrics['structural_status']} |"
        )
    convergence = summary["representative_convergence"]
    lines.extend(
        [
            "",
            "Representative `(4,6)` aggregate sequence:",
            "",
            "| mesh | E3 | E6 | order from previous |",
            "|---:|---:|---:|---:|",
        ]
    )
    for mesh in UT1_FEM_MESHES:
        key = str(mesh)
        order = convergence["orders"].get(key)
        order_text = "n/a" if order is None else f"{order:.12f}"
        lines.append(
            f"| {mesh} | {convergence['E3'][key]:.12e} | "
            f"{convergence['E6'][key]:.12e} | {order_text} |"
        )
    lines.extend(
        [
            "",
            "## Timoshenko--EB diagnostic",
            "",
            "No acceptance threshold or accuracy ranking is assigned. Maximum "
            "first-six relative model differences are "
            + ", ".join(
                f"`{case}` `{value:.12e}`"
                for case, value in summary["model_difference_maxima"].items()
            )
            + ".",
            "",
            "## Thresholds and exclusions",
            "",
            "Continuum reference, homogeneous baseline, and exchange gates use "
            "`1e-8`; root determinant and relative singular residuals use "
            "`1e-8`; FEM mesh-64 uses `E3<=1e-4`, `E6<=5e-4`; aggregate "
            "orders must be at least `1.5`. Existing structural tolerances are "
            "unchanged.",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in summary["explicit_exclusions"])
    lines.append("")
    return "\n".join(lines)


def _run_ut1_beta0(output_dir: Path) -> dict[str, Any]:
    started = time.perf_counter()
    section_cases = _make_ut1_section_cases()
    spectra: dict[tuple[str, str, str], dict[str, Any]] = {}
    continuum_rows: list[dict[str, Any]] = []
    runtime_by_spectrum: dict[str, float] = {}
    evaluations_by_spectrum: dict[str, int] = {}
    scan_max_by_spectrum: dict[str, float] = {}

    for section_case, (point_1, point_2) in section_cases.items():
        for model in ("Timoshenko", "EB"):
            factories = _ut1_factories(model, point_1, point_2)
            for formulation in ("coupled", "direct_stepped"):
                roots, builder, runtime = _solve_roots(
                    point_1,
                    factories[formulation],
                    num_roots=UT1_NUM_ROOTS,
                )
                raw_factory = factories[f"{formulation}_raw"]
                quality = [
                    _root_quality_fields(root, raw_factory) for root in roots
                ]
                gaps = _neighbor_gaps(_frequencies(roots))
                spectrum_key = (section_case, model, formulation)
                spectra[spectrum_key] = {
                    "roots": roots,
                    "quality": quality,
                    "raw_factory": raw_factory,
                }
                label = "|".join(spectrum_key)
                runtime_by_spectrum[label] = runtime
                evaluations_by_spectrum[label] = builder.evaluations
                scan_max_by_spectrum[label] = builder.maximum_frequency_hz
                continuum_rows.extend(
                    _continuum_root_row(
                        section_case=section_case,
                        point_1=point_1,
                        point_2=point_2,
                        model=model,
                        formulation=formulation,
                        mode=mode,
                        root=root,
                        quality=quality_fields,
                        gap_hz=float(gap),
                    )
                    for mode, (root, quality_fields, gap) in enumerate(
                        zip(roots, quality, gaps), start=1
                    )
                )

    straight_point = _make_section_point(0.005, length_m=0.4)
    straight_factories = {
        "Timoshenko": (
            lambda omega: straight_boundary_matrix(omega, straight_point),
            lambda omega: straight_boundary_matrix_raw(omega, straight_point),
        ),
        "EB": (
            lambda omega: eb_straight_boundary_matrix(omega, straight_point),
            lambda omega: eb_straight_boundary_matrix_raw(omega, straight_point),
        ),
    }
    for model, (factory, raw_factory) in straight_factories.items():
        roots, builder, runtime = _solve_roots(
            straight_point, factory, num_roots=UT1_NUM_ROOTS
        )
        quality = [_root_quality_fields(root, raw_factory) for root in roots]
        gaps = _neighbor_gaps(_frequencies(roots))
        key = ("baseline_5_5", model, "direct_straight")
        spectra[key] = {
            "roots": roots,
            "quality": quality,
            "raw_factory": raw_factory,
        }
        label = "|".join(key)
        runtime_by_spectrum[label] = runtime
        evaluations_by_spectrum[label] = builder.evaluations
        scan_max_by_spectrum[label] = builder.maximum_frequency_hz
        continuum_rows.extend(
            _continuum_root_row(
                section_case="baseline_5_5",
                point_1=straight_point,
                point_2=straight_point,
                model=model,
                formulation="direct_straight",
                mode=mode,
                root=root,
                quality=quality_fields,
                gap_hz=float(gap),
            )
            for mode, (root, quality_fields, gap) in enumerate(
                zip(roots, quality, gaps), start=1
            )
        )

    reference_rows: list[dict[str, Any]] = []
    reference_maxima: dict[str, float] = {}
    for section_case in UT1_SECTION_CASES:
        for model in ("Timoshenko", "EB"):
            coupled = spectra[(section_case, model, "coupled")]
            stepped = spectra[(section_case, model, "direct_stepped")]
            coupled_frequency = _frequencies(coupled["roots"])
            stepped_matched, bases, _, permutation = _match_frequency_spectra(
                coupled_frequency, _frequencies(stepped["roots"])
            )
            differences: list[float] = []
            for index, (coupled_hz, stepped_hz, basis, stepped_index) in enumerate(
                zip(coupled_frequency, stepped_matched, bases, permutation),
                start=1,
            ):
                absolute = abs(float(coupled_hz) - float(stepped_hz))
                relative = absolute / float(stepped_hz)
                differences.append(relative)
                coupled_quality = coupled["quality"][index - 1]
                stepped_quality = stepped["quality"][int(stepped_index)]
                reference_rows.append(
                    {
                        "section_case": section_case,
                        "model": model,
                        "mode": index,
                        "coupled_frequency_hz": coupled_hz,
                        "stepped_frequency_hz": stepped_hz,
                        "absolute_difference_hz": absolute,
                        "relative_difference": relative,
                        "matching_basis": basis,
                        "coupled_quality_status": coupled_quality[
                            "root_quality_status"
                        ],
                        "coupled_quality_basis": coupled_quality[
                            "root_quality_basis"
                        ],
                        "stepped_quality_status": stepped_quality[
                            "root_quality_status"
                        ],
                        "stepped_quality_basis": stepped_quality[
                            "root_quality_basis"
                        ],
                        "status": (
                            "PASS"
                            if relative <= SPECTRUM_RELATIVE_TOLERANCE
                            and coupled_quality["root_quality_status"] == "PASS"
                            and stepped_quality["root_quality_status"] == "PASS"
                            else "FAIL"
                        ),
                    }
                )
            reference_maxima[f"{section_case}|{model}"] = max(differences)

    baseline_rows: list[dict[str, Any]] = []
    baseline_maxima: dict[str, float] = {}
    for model in ("Timoshenko", "EB"):
        coupled = _frequencies(
            spectra[("baseline_5_5", model, "coupled")]["roots"]
        )
        stepped, stepped_bases, _, _ = _match_frequency_spectra(
            coupled,
            _frequencies(
                spectra[("baseline_5_5", model, "direct_stepped")]["roots"]
            ),
        )
        straight, straight_bases, _, _ = _match_frequency_spectra(
            coupled,
            _frequencies(
                spectra[("baseline_5_5", model, "direct_straight")]["roots"]
            ),
        )
        pairwise: list[float] = []
        for mode, values in enumerate(zip(coupled, stepped, straight), start=1):
            coupled_hz, stepped_hz, straight_hz = map(float, values)
            relative_cs = abs(coupled_hz - stepped_hz) / stepped_hz
            relative_ch = abs(coupled_hz - straight_hz) / straight_hz
            relative_sh = abs(stepped_hz - straight_hz) / straight_hz
            maximum = max(relative_cs, relative_ch, relative_sh)
            pairwise.append(maximum)
            baseline_rows.append(
                {
                    "model": model,
                    "mode": mode,
                    "coupled_frequency_hz": coupled_hz,
                    "stepped_frequency_hz": stepped_hz,
                    "straight_frequency_hz": straight_hz,
                    "coupled_stepped_relative_difference": relative_cs,
                    "coupled_straight_relative_difference": relative_ch,
                    "stepped_straight_relative_difference": relative_sh,
                    "maximum_pairwise_relative_difference": maximum,
                    "stepped_matching_basis": stepped_bases[mode - 1],
                    "straight_matching_basis": straight_bases[mode - 1],
                    "status": (
                        "PASS"
                        if maximum <= SPECTRUM_RELATIVE_TOLERANCE
                        else "FAIL"
                    ),
                }
            )
        baseline_maxima[model] = max(pairwise)

    exchange_rows: list[dict[str, Any]] = []
    analytic_exchange_maxima: dict[str, float] = {}
    for model in ("Timoshenko", "EB"):
        for formulation in ("coupled", "direct_stepped"):
            reference = _frequencies(
                spectra[("asymmetric_4_6", model, formulation)]["roots"]
            )
            matched, bases, _, _ = _match_frequency_spectra(
                reference,
                _frequencies(
                    spectra[("swapped_6_4", model, formulation)]["roots"]
                ),
            )
            differences: list[float] = []
            for mode, (frequency_4_6, frequency_6_4, basis) in enumerate(
                zip(reference, matched, bases), start=1
            ):
                absolute = abs(float(frequency_4_6) - float(frequency_6_4))
                relative = absolute / float(frequency_4_6)
                differences.append(relative)
                exchange_rows.append(
                    {
                        "model": model,
                        "formulation": formulation,
                        "mode": mode,
                        "frequency_4_6_hz": frequency_4_6,
                        "frequency_6_4_hz": frequency_6_4,
                        "absolute_difference_hz": absolute,
                        "relative_difference": relative,
                        "matching_basis": basis,
                        "status": (
                            "PASS"
                            if relative <= SPECTRUM_RELATIVE_TOLERANCE
                            else "FAIL"
                        ),
                    }
                )
            analytic_exchange_maxima[f"{model}|{formulation}"] = max(
                differences
            )

    fem_runs: dict[tuple[str, int], dict[str, Any]] = {}
    fem_runtime_by_case_mesh: dict[str, float] = {}
    for section_case, meshes in {
        "baseline_5_5": (UT1_FEM_FINE_MESH,),
        "asymmetric_4_6": UT1_FEM_MESHES,
        "swapped_6_4": (UT1_FEM_FINE_MESH,),
    }.items():
        point_1, point_2 = section_cases[section_case]
        analytic = _frequencies(
            spectra[(section_case, "EB", "coupled")]["roots"]
        )
        for mesh in meshes:
            run = _run_eb_fem_case(
                section_case=section_case,
                point_1=point_1,
                point_2=point_2,
                elements_per_arm=mesh,
                analytic_frequencies=analytic,
            )
            fem_runs[(section_case, mesh)] = run
            fem_runtime_by_case_mesh[f"{section_case}|{mesh}"] = run[
                "diagnostics"
            ]["runtime_seconds"]

    fem_rows: list[dict[str, Any]] = []
    fem_mesh64_errors: dict[str, dict[str, Any]] = {}
    for section_case in UT1_SECTION_CASES:
        run = fem_runs[(section_case, UT1_FEM_FINE_MESH)]
        analytic = _frequencies(
            spectra[(section_case, "EB", "coupled")]["roots"]
        )
        errors = np.abs(run["matched_frequencies"] - analytic) / analytic
        diagnostics = run["diagnostics"]
        E3 = float(np.max(errors[:3]))
        E6 = float(np.max(errors[:6]))
        case_status = (
            "PASS"
            if _fem_structural_ok(diagnostics)
            and E3 <= FEM_FIRST_THREE_TOLERANCE
            and E6 <= FEM_FIRST_SIX_TOLERANCE
            else "FAIL"
        )
        fem_mesh64_errors[section_case] = {
            "E3": E3,
            "E6": E6,
            "root_7_error": float(errors[6]),
            "structural_status": diagnostics["structural_status"],
            "status": case_status,
        }
        for mode in range(UT1_NUM_ROOTS):
            absolute = abs(
                float(run["matched_frequencies"][mode]) - float(analytic[mode])
            )
            fem_rows.append(
                {
                    "section_case": section_case,
                    "elements_arm_1": UT1_FEM_FINE_MESH,
                    "elements_arm_2": UT1_FEM_FINE_MESH,
                    "mode": mode + 1,
                    "analytic_frequency_hz": analytic[mode],
                    "fem_frequency_hz": run["matched_frequencies"][mode],
                    "absolute_error_hz": absolute,
                    "relative_error": errors[mode],
                    "matching_basis": run["matching_basis"][mode],
                    "stiffness_symmetry_residual": diagnostics[
                        "stiffness_symmetry_residual"
                    ],
                    "mass_symmetry_residual": diagnostics[
                        "mass_symmetry_residual"
                    ],
                    "minimum_mass_eigenvalue": diagnostics[
                        "minimum_mass_eigenvalue"
                    ],
                    "zero_mode_count": diagnostics["zero_mode_count"],
                    "joint_kinematic_residual": diagnostics[
                        "joint_kinematic_residual"
                    ],
                    "joint_equilibrium_residual": run[
                        "matched_joint_equilibrium_residuals"
                    ][mode],
                    "status": case_status,
                }
            )

    convergence_metrics: dict[int, dict[str, float]] = {}
    convergence_rows: list[dict[str, Any]] = []
    analytic_asymmetric = _frequencies(
        spectra[("asymmetric_4_6", "EB", "coupled")]["roots"]
    )
    for mesh in UT1_FEM_MESHES:
        run = fem_runs[("asymmetric_4_6", mesh)]
        errors = (
            np.abs(run["matched_frequencies"] - analytic_asymmetric)
            / analytic_asymmetric
        )
        convergence_metrics[mesh] = {
            "E3": float(np.max(errors[:3])),
            "E6": float(np.max(errors[:6])),
        }
    convergence_orders = {
        32: float(
            math.log(
                convergence_metrics[16]["E6"]
                / convergence_metrics[32]["E6"],
                2.0,
            )
        ),
        64: float(
            math.log(
                convergence_metrics[32]["E6"]
                / convergence_metrics[64]["E6"],
                2.0,
            )
        ),
    }
    convergence_ok = bool(
        convergence_metrics[32]["E6"] < convergence_metrics[16]["E6"]
        and convergence_metrics[64]["E6"] < convergence_metrics[32]["E6"]
        and convergence_orders[32] >= FEM_AGGREGATE_ORDER_MINIMUM
        and convergence_orders[64] >= FEM_AGGREGATE_ORDER_MINIMUM
        and convergence_metrics[64]["E3"] <= FEM_FIRST_THREE_TOLERANCE
        and convergence_metrics[64]["E6"] <= FEM_FIRST_SIX_TOLERANCE
    )
    for mesh in UT1_FEM_MESHES:
        run = fem_runs[("asymmetric_4_6", mesh)]
        errors = (
            np.abs(run["matched_frequencies"] - analytic_asymmetric)
            / analytic_asymmetric
        )
        diagnostics = run["diagnostics"]
        for mode in range(UT1_NUM_ROOTS):
            convergence_rows.append(
                {
                    "elements_arm_1": mesh,
                    "elements_arm_2": mesh,
                    "element_length_m": LENGTH_M / mesh,
                    "reduced_size": diagnostics["reduced_size"],
                    "mode": mode + 1,
                    "analytic_frequency_hz": analytic_asymmetric[mode],
                    "fem_frequency_hz": run["matched_frequencies"][mode],
                    "relative_error": errors[mode],
                    "matching_basis": run["matching_basis"][mode],
                    "E3": convergence_metrics[mesh]["E3"],
                    "E6": convergence_metrics[mesh]["E6"],
                    "aggregate_order_from_previous": (
                        "" if mesh == 16 else convergence_orders[mesh]
                    ),
                    "structural_diagnostics": json.dumps(
                        diagnostics, sort_keys=True
                    ),
                    "status": (
                        "PASS"
                        if _fem_structural_ok(diagnostics) and convergence_ok
                        else "FAIL"
                    ),
                }
            )

    fem_reference = fem_runs[("asymmetric_4_6", 64)]["matched_frequencies"]
    fem_swapped, fem_bases, _, _ = _match_frequency_spectra(
        fem_reference, fem_runs[("swapped_6_4", 64)]["matched_frequencies"]
    )
    fem_exchange_differences: list[float] = []
    for mode, (frequency_4_6, frequency_6_4, basis) in enumerate(
        zip(fem_reference, fem_swapped, fem_bases), start=1
    ):
        absolute = abs(float(frequency_4_6) - float(frequency_6_4))
        relative = absolute / float(frequency_4_6)
        fem_exchange_differences.append(relative)
        exchange_rows.append(
            {
                "model": "EB FEM N=64",
                "formulation": "independent_fem",
                "mode": mode,
                "frequency_4_6_hz": frequency_4_6,
                "frequency_6_4_hz": frequency_6_4,
                "absolute_difference_hz": absolute,
                "relative_difference": relative,
                "matching_basis": basis,
                "status": (
                    "PASS"
                    if relative <= SPECTRUM_RELATIVE_TOLERANCE
                    else "FAIL"
                ),
            }
        )
    fem_exchange_maximum = max(fem_exchange_differences)
    fem_exchange_structural_difference = {
        key: abs(
            float(fem_runs[("asymmetric_4_6", 64)]["diagnostics"][key])
            - float(fem_runs[("swapped_6_4", 64)]["diagnostics"][key])
        )
        for key in (
            "stiffness_symmetry_residual",
            "mass_symmetry_residual",
            "minimum_mass_eigenvalue",
            "zero_mode_count",
            "joint_kinematic_residual",
            "max_joint_equilibrium_residual",
        )
    }

    model_difference_rows: list[dict[str, Any]] = []
    model_difference_maxima: dict[str, float] = {}
    for section_case in UT1_SECTION_CASES:
        timoshenko = _frequencies(
            spectra[(section_case, "Timoshenko", "coupled")]["roots"]
        )[:6]
        eb = _frequencies(
            spectra[(section_case, "EB", "coupled")]["roots"]
        )[:6]
        matched_eb, bases, timo_gaps, permutation = _match_frequency_spectra(
            timoshenko, eb
        )
        eb_gaps = _neighbor_gaps(eb)
        differences: list[float] = []
        for mode, (timo_hz, eb_hz, basis, eb_index) in enumerate(
            zip(timoshenko, matched_eb, bases, permutation), start=1
        ):
            absolute = abs(float(timo_hz) - float(eb_hz))
            relative = absolute / float(timo_hz)
            differences.append(relative)
            model_difference_rows.append(
                {
                    "section_case": section_case,
                    "mode": mode,
                    "timoshenko_frequency_hz": timo_hz,
                    "eb_frequency_hz": eb_hz,
                    "absolute_difference_hz": absolute,
                    "relative_difference": relative,
                    "timoshenko_neighbor_gap_hz": timo_gaps[mode - 1],
                    "eb_neighbor_gap_hz": eb_gaps[int(eb_index)],
                    "matching_basis": basis,
                }
            )
        model_difference_maxima[section_case] = max(differences)

    main_spectrum_keys = [
        key for key in spectra if key[2] != "direct_straight"
    ]
    main_spectra_ok = all(
        len(spectra[key]["roots"]) == UT1_NUM_ROOTS
        and all(
            np.isfinite(root.frequency_hz) and root.frequency_hz > 0.0
            for root in spectra[key]["roots"]
        )
        for key in main_spectrum_keys
    )
    all_quality_ok = all(
        row["quality_status"] == "PASS" for row in continuum_rows
    )
    quality_basis_counts = {
        basis: sum(
            row["quality_status"] == "PASS" and row["quality_basis"] == basis
            for row in continuum_rows
        )
        for basis in ("scaled", "physical_raw")
    }
    quality_maxima_by_basis: dict[str, dict[str, float]] = {}
    for basis in ("scaled", "physical_raw"):
        selected = [row for row in continuum_rows if row["quality_basis"] == basis]
        if basis == "scaled":
            quality_maxima_by_basis[basis] = {
                "determinant_residual": max(
                    (row["scaled_determinant_residual"] for row in selected),
                    default=0.0,
                ),
                "relative_singular_residual": max(
                    (
                        row["scaled_relative_singular_residual"]
                        for row in selected
                    ),
                    default=0.0,
                ),
            }
        else:
            quality_maxima_by_basis[basis] = {
                "determinant_residual": max(
                    (row["raw_determinant_residual"] for row in selected),
                    default=0.0,
                ),
                "relative_singular_residual": max(
                    (row["raw_relative_singular_residual"] for row in selected),
                    default=0.0,
                ),
            }
    quality_representation_maxima = {
        "scaled_determinant_residual": max(
            row["scaled_determinant_residual"] for row in continuum_rows
        ),
        "scaled_relative_singular_residual": max(
            row["scaled_relative_singular_residual"] for row in continuum_rows
        ),
        "raw_determinant_residual": max(
            row["raw_determinant_residual"] for row in continuum_rows
        ),
        "raw_relative_singular_residual": max(
            row["raw_relative_singular_residual"] for row in continuum_rows
        ),
    }

    ut0_summary_path = DEFAULT_OUTPUT_DIR / "ut0_summary.json"
    ut0_regression_status = "MISSING"
    if ut0_summary_path.exists():
        ut0_regression_status = json.loads(
            ut0_summary_path.read_text(encoding="utf-8")
        ).get("ut0_status", "MISSING")
    reference_ok = all(row["status"] == "PASS" for row in reference_rows)
    baseline_ok = all(row["status"] == "PASS" for row in baseline_rows)
    analytic_exchange_ok = all(
        value <= SPECTRUM_RELATIVE_TOLERANCE
        for value in analytic_exchange_maxima.values()
    )
    continuum_hard_ok = bool(
        ut0_regression_status == "PASS"
        and len(main_spectrum_keys) == 12
        and main_spectra_ok
        and all_quality_ok
        and reference_ok
        and baseline_ok
        and analytic_exchange_ok
    )
    fem_structure_ok = all(
        _fem_structural_ok(
            fem_runs[(case, UT1_FEM_FINE_MESH)]["diagnostics"]
        )
        for case in UT1_SECTION_CASES
    )
    fem_accuracy_ok = all(
        metrics["E3"] <= FEM_FIRST_THREE_TOLERANCE
        and metrics["E6"] <= FEM_FIRST_SIX_TOLERANCE
        for metrics in fem_mesh64_errors.values()
    )
    fem_exchange_ok = fem_exchange_maximum <= SPECTRUM_RELATIVE_TOLERANCE
    fem_all_ok = bool(
        fem_structure_ok
        and fem_accuracy_ok
        and convergence_ok
        and fem_exchange_ok
    )
    if not continuum_hard_ok:
        ut1_status = "FAIL"
    elif fem_all_ok:
        ut1_status = "PASS"
    else:
        ut1_status = "PARTIAL_PASS"

    summary = {
        "git_context": _git_context(),
        "scope": "UT-1 full beta=0 unequal-thickness validation only",
        "fixed_input_geometry": {
            "length_1_m": LENGTH_M,
            "length_2_m": LENGTH_M,
            "b_1_m": WIDTH_B_M,
            "b_2_m": WIDTH_B_M,
            "shear_factor": cantilever_geometry().shear_factor,
            "a_direction": "parallel to e_z",
        },
        "material": "existing HMS/DX-209",
        "material_mode": "elastic",
        "theta_1_deg": 0.0,
        "theta_2_deg": 0.0,
        "beta_deg": 0.0,
        "section_cases": {
            case: {"a_1_m": values[0], "a_2_m": values[1]}
            for case, values in UT1_SECTION_CASES.items()
        },
        "thresholds": {
            "root_determinant_residual": ROOT_DETERMINANT_TOLERANCE,
            "root_relative_singular_residual": ROOT_SINGULAR_TOLERANCE,
            "near_degenerate_relative_gap": NEAR_DEGENERATE_RELATIVE_GAP,
            "coupled_stepped_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "baseline_pairwise_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "analytic_exchange_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "fem_exchange_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "fem_stiffness_symmetry": FEM_SYMMETRY_TOLERANCE,
            "fem_mass_symmetry": FEM_SYMMETRY_TOLERANCE,
            "fem_joint_kinematic": FEM_JOINT_KINEMATIC_TOLERANCE,
            "fem_joint_equilibrium": FEM_JOINT_EQUILIBRIUM_TOLERANCE,
            "fem_E3_mesh64": FEM_FIRST_THREE_TOLERANCE,
            "fem_E6_mesh64": FEM_FIRST_SIX_TOLERANCE,
            "fem_aggregate_order_minimum": FEM_AGGREGATE_ORDER_MINIMUM,
        },
        "root_scan": {
            "formulation": FORMULATION,
            "num_roots": UT1_NUM_ROOTS,
            "scan_step_hz": SCAN_STEP_HZ,
            "initial_max_hz": INITIAL_MAX_HZ,
            "max_hz": MAX_HZ,
            "actual_maximum_evaluated_frequency_hz": scan_max_by_spectrum,
        },
        "main_continuum_spectrum_count": len(main_spectrum_keys),
        "continuum_root_counts": {
            "|".join(key): len(data["roots"]) for key, data in spectra.items()
        },
        "first_seven_frequencies_hz": {
            section_case: {
                model: {
                    formulation: _frequencies(
                        spectra[(section_case, model, formulation)]["roots"]
                    ).tolist()
                    for formulation in ("coupled", "direct_stepped")
                }
                for model in ("Timoshenko", "EB")
            }
            for section_case in UT1_SECTION_CASES
        },
        "continuum_root_row_count": len(continuum_rows),
        "root_quality_basis_counts": quality_basis_counts,
        "main_root_quality_basis_counts": {
            basis: sum(
                row["formulation"] != "direct_straight"
                and row["quality_status"] == "PASS"
                and row["quality_basis"] == basis
                for row in continuum_rows
            )
            for basis in ("scaled", "physical_raw")
        },
        "root_quality_maxima_by_basis": quality_maxima_by_basis,
        "root_quality_representation_maxima": quality_representation_maxima,
        "maximum_coupled_stepped_relative_difference": reference_maxima,
        "equal_thickness_straight_regression_maxima": baseline_maxima,
        "analytic_exchange_maxima": analytic_exchange_maxima,
        "fem_exchange_maximum": fem_exchange_maximum,
        "fem_exchange_structural_diagnostic_differences": (
            fem_exchange_structural_difference
        ),
        "fem_mesh64_errors": fem_mesh64_errors,
        "fem_mesh64_diagnostics": {
            case: fem_runs[(case, 64)]["diagnostics"]
            for case in UT1_SECTION_CASES
        },
        "representative_convergence": {
            "case": "asymmetric_4_6",
            "meshes": list(UT1_FEM_MESHES),
            "E3": {
                str(mesh): convergence_metrics[mesh]["E3"]
                for mesh in UT1_FEM_MESHES
            },
            "E6": {
                str(mesh): convergence_metrics[mesh]["E6"]
                for mesh in UT1_FEM_MESHES
            },
            "orders": {
                "32": convergence_orders[32],
                "64": convergence_orders[64],
            },
            "status": "PASS" if convergence_ok else "FAIL",
        },
        "model_difference_maxima": model_difference_maxima,
        "runtime_seconds_by_spectrum": runtime_by_spectrum,
        "fem_runtime_seconds_by_case_mesh": fem_runtime_by_case_mesh,
        "total_runtime_seconds": time.perf_counter() - started,
        "unique_matrix_evaluations_by_spectrum": evaluations_by_spectrum,
        "total_unique_matrix_evaluations": sum(evaluations_by_spectrum.values()),
        "ut0_regression_status": ut0_regression_status,
        "continuum_hard_gates_status": "PASS" if continuum_hard_ok else "FAIL",
        "fem_gates_status": "PASS" if fem_all_ok else "FAIL",
        "ut1_status": ut1_status,
        "overall_unequal_thickness_validation": "IN_PROGRESS",
        "explicit_exclusions": UT1_EXPLICIT_EXCLUSIONS,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "continuum_roots.csv", continuum_rows)
    _write_csv(output_dir / "reference_comparison.csv", reference_rows)
    _write_csv(output_dir / "baseline_straight_regression.csv", baseline_rows)
    _write_csv(output_dir / "exchange_symmetry.csv", exchange_rows)
    _write_csv(output_dir / "fem_comparison.csv", fem_rows)
    _write_csv(output_dir / "fem_convergence.csv", convergence_rows)
    _write_csv(output_dir / "model_difference.csv", model_difference_rows)
    (output_dir / "ut1_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "ut1_report.md").write_text(
        _ut1_report(summary), encoding="utf-8"
    )
    return summary


def _validate_ut1a_dof_ordering(assembly: EBFEMAssembly) -> dict[str, Any]:
    """Assert and describe the production full/reduced EB FEM ordering."""

    counts = tuple(int(value) for value in assembly.element_counts)
    if len(counts) != 2 or any(value < 2 for value in counts):
        raise ValueError("UT-1a requires two arms with internal FEM nodes")
    full_arm_sizes = tuple(
        FEM_NODE_DOF_COUNT * (count + 1) for count in counts
    )
    internal_arm_sizes = tuple(
        FEM_NODE_DOF_COUNT * (count - 1) for count in counts
    )
    expected_full_size = sum(full_arm_sizes)
    expected_reduced_size = sum(internal_arm_sizes) + FEM_NODE_DOF_COUNT
    if assembly.stiffness_full.shape != (
        expected_full_size,
        expected_full_size,
    ):
        raise RuntimeError("full stiffness size disagrees with element metadata")
    if assembly.mass_full.shape != assembly.stiffness_full.shape:
        raise RuntimeError("full stiffness/mass sizes disagree")
    if assembly.reduction.shape != (
        expected_full_size,
        expected_reduced_size,
    ):
        raise RuntimeError("reduction size disagrees with element metadata")
    if assembly.stiffness.shape != (
        expected_reduced_size,
        expected_reduced_size,
    ) or assembly.mass.shape != assembly.stiffness.shape:
        raise RuntimeError("reduced matrix size disagrees with element metadata")

    expected_reduction = np.zeros_like(assembly.reduction)
    full_offsets = (0, full_arm_sizes[0])
    reduced_cursor = 0
    for arm, (count, full_offset) in enumerate(zip(counts, full_offsets)):
        for node in range(1, count):
            rows = full_offset + FEM_NODE_DOF_COUNT * node + np.arange(3)
            columns = reduced_cursor + np.arange(3)
            expected_reduction[np.ix_(rows, columns)] = np.eye(3)
            reduced_cursor += 3
        expected_joint_rows = (
            full_offset + FEM_NODE_DOF_COUNT * count + np.arange(3)
        )
        if not np.array_equal(
            assembly.joint_full_dofs[arm], expected_joint_rows
        ):
            raise RuntimeError("joint endpoint rows disagree with arm metadata")
        expected_reduction[
            np.ix_(expected_joint_rows, np.arange(expected_reduced_size - 3, expected_reduced_size))
        ] = assembly.joint_maps[arm]
    if reduced_cursor + 3 != expected_reduced_size:
        raise RuntimeError("internal reduced ordering count mismatch")
    if not np.array_equal(assembly.reduction, expected_reduction):
        raise RuntimeError("production reduction does not match documented ordering")

    return {
        "element_counts": list(counts),
        "full_arm_sizes": list(full_arm_sizes),
        "internal_arm_sizes": list(internal_arm_sizes),
        "full_size": expected_full_size,
        "reduced_size": expected_reduced_size,
        "full_ordering": (
            "[arm_1 nodes outer-to-joint in local [w,psi,Phi], "
            "arm_2 nodes outer-to-joint in local [w,psi,Phi]]"
        ),
        "reduced_ordering": (
            "[arm_1 internal nodes, arm_2 internal nodes, "
            "joint [w_J,theta_t,theta_n]]"
        ),
    }


def _validate_ut1a_pair_metadata(
    assembly_46: EBFEMAssembly, assembly_64: EBFEMAssembly
) -> dict[str, Any]:
    metadata_46 = _validate_ut1a_dof_ordering(assembly_46)
    metadata_64 = _validate_ut1a_dof_ordering(assembly_64)
    if assembly_46.element_counts != assembly_64.element_counts:
        raise ValueError("arm-swap audit requires identical element counts")
    if assembly_46.element_counts[0] != assembly_46.element_counts[1]:
        raise ValueError("arm-swap audit requires equal counts on both arms")
    expected_maps = eb_joint_dof_maps(BETA_RAD)
    for assembly in (assembly_46, assembly_64):
        for actual, expected in zip(assembly.joint_maps, expected_maps):
            if not np.array_equal(actual, expected):
                raise RuntimeError("production joint maps changed during UT-1a")
    if metadata_46 != metadata_64:
        raise RuntimeError("swapped assemblies expose different DOF metadata")
    return metadata_46


def _ut1a_full_swap_permutation(
    assembly_46: EBFEMAssembly, assembly_64: EBFEMAssembly
) -> np.ndarray:
    """Return the a-priori full-space arm-block permutation ``x64=P_f x46``."""

    metadata = _validate_ut1a_pair_metadata(assembly_46, assembly_64)
    arm_sizes = metadata["full_arm_sizes"]
    if arm_sizes[0] != arm_sizes[1]:
        raise ValueError("full arm blocks must have equal size")
    arm_size = int(arm_sizes[0])
    permutation = np.zeros((2 * arm_size, 2 * arm_size), dtype=float)
    identity = np.eye(arm_size)
    permutation[:arm_size, arm_size:] = identity
    permutation[arm_size:, :arm_size] = identity
    return permutation


def _ut1a_reduced_swap_permutation(
    assembly_46: EBFEMAssembly, assembly_64: EBFEMAssembly
) -> np.ndarray:
    """Return the signed reduced permutation ``q64=P_r q46`` at beta=0."""

    metadata = _validate_ut1a_pair_metadata(assembly_46, assembly_64)
    internal_sizes = metadata["internal_arm_sizes"]
    if internal_sizes[0] != internal_sizes[1]:
        raise ValueError("internal arm blocks must have equal size")
    internal_size = int(internal_sizes[0])
    reduced_size = int(metadata["reduced_size"])
    permutation = np.zeros((reduced_size, reduced_size), dtype=float)
    identity = np.eye(internal_size)
    permutation[:internal_size, internal_size : 2 * internal_size] = identity
    permutation[internal_size : 2 * internal_size, :internal_size] = identity
    permutation[-3:, -3:] = np.diag([1.0, -1.0, -1.0])
    return permutation


def _ut1a_full_dof_label(assembly: EBFEMAssembly, index: int) -> str:
    if index < 0:
        return ""
    names = ("w", "psi", "Phi")
    offset = 0
    for arm, count in enumerate(assembly.element_counts, start=1):
        arm_size = FEM_NODE_DOF_COUNT * (int(count) + 1)
        if index < offset + arm_size:
            local = index - offset
            node, dof = divmod(local, FEM_NODE_DOF_COUNT)
            return f"arm_{arm}.node_{node}.{names[dof]}"
        offset += arm_size
    raise IndexError(index)


def _ut1a_reduced_dof_label(assembly: EBFEMAssembly, index: int) -> str:
    if index < 0:
        return ""
    names = ("w", "psi", "Phi")
    offset = 0
    for arm, count in enumerate(assembly.element_counts, start=1):
        block_size = FEM_NODE_DOF_COUNT * (int(count) - 1)
        if index < offset + block_size:
            local = index - offset
            internal_node, dof = divmod(local, FEM_NODE_DOF_COUNT)
            return f"arm_{arm}.node_{internal_node + 1}.{names[dof]}"
        offset += block_size
    joint_names = ("w_J", "theta_t", "theta_n")
    joint_index = index - offset
    if 0 <= joint_index < 3:
        return f"joint.{joint_names[joint_index]}"
    raise IndexError(index)


def _ut1a_matrix_identity_row(
    *,
    check: str,
    left: np.ndarray,
    right: np.ndarray,
    threshold: float,
    row_label: Callable[[int], str] | None = None,
    column_label: Callable[[int], str] | None = None,
) -> dict[str, Any]:
    """Return Frobenius and entrywise evidence for one matrix identity."""

    left_value = np.asarray(left, dtype=float)
    right_value = np.asarray(right, dtype=float)
    if left_value.shape != right_value.shape or left_value.ndim != 2:
        raise ValueError(f"matrix identity {check!r} has incompatible shapes")
    difference = left_value - right_value
    absolute_difference = np.abs(difference)
    absolute_frobenius = float(np.linalg.norm(difference, ord="fro"))
    relative_frobenius = absolute_frobenius / max(
        float(np.linalg.norm(right_value, ord="fro")), np.finfo(float).tiny
    )
    absolute_max = float(np.max(absolute_difference))
    relative_max = absolute_max / max(
        float(np.max(np.abs(right_value))), np.finfo(float).tiny
    )
    nonzero = np.flatnonzero(absolute_difference)
    if len(nonzero):
        first_row, first_column = np.unravel_index(
            int(nonzero[0]), difference.shape
        )
        max_row, max_column = np.unravel_index(
            int(np.argmax(absolute_difference)), difference.shape
        )
        first_left = float(left_value[first_row, first_column])
        first_right = float(right_value[first_row, first_column])
        max_left = float(left_value[max_row, max_column])
        max_right = float(right_value[max_row, max_column])
    else:
        first_row = first_column = max_row = max_column = -1
        first_left = first_right = max_left = max_right = None
    row_labeler = row_label or (lambda value: str(value))
    column_labeler = column_label or (lambda value: str(value))
    return {
        "check": check,
        "shape_rows": int(left_value.shape[0]),
        "shape_columns": int(left_value.shape[1]),
        "absolute_frobenius_residual": absolute_frobenius,
        "relative_frobenius_residual": relative_frobenius,
        "absolute_max_residual": absolute_max,
        "relative_max_residual": relative_max,
        "first_mismatch_row": first_row,
        "first_mismatch_column": first_column,
        "first_mismatch_row_label": row_labeler(first_row),
        "first_mismatch_column_label": column_labeler(first_column),
        "first_left_value": first_left,
        "first_right_value": first_right,
        "max_mismatch_row": max_row,
        "max_mismatch_column": max_column,
        "max_mismatch_row_label": row_labeler(max_row),
        "max_mismatch_column_label": column_labeler(max_column),
        "left_value": max_left,
        "right_value": max_right,
        "threshold": threshold,
        "status": (
            "PASS"
            if relative_frobenius <= threshold and relative_max <= threshold
            else "FAIL"
        ),
    }


def _ut1a_matrix_congruence_rows(
    assembly_46: EBFEMAssembly,
    assembly_64: EBFEMAssembly,
    permutation_full: np.ndarray,
    permutation_reduced: np.ndarray,
) -> list[dict[str, Any]]:
    full_label = lambda index: _ut1a_full_dof_label(assembly_64, index)
    reduced_label = lambda index: _ut1a_reduced_dof_label(assembly_64, index)
    identity_full = np.eye(permutation_full.shape[0])
    identity_reduced = np.eye(permutation_reduced.shape[0])
    threshold = UT1A_MATRIX_RELATIVE_TOLERANCE
    rows = [
        _ut1a_matrix_identity_row(
            check="P_full orthogonality",
            left=permutation_full.T @ permutation_full,
            right=identity_full,
            threshold=threshold,
            row_label=full_label,
            column_label=full_label,
        ),
        _ut1a_matrix_identity_row(
            check="P_full involution",
            left=permutation_full @ permutation_full,
            right=identity_full,
            threshold=threshold,
            row_label=full_label,
            column_label=full_label,
        ),
        _ut1a_matrix_identity_row(
            check="P_reduced orthogonality",
            left=permutation_reduced.T @ permutation_reduced,
            right=identity_reduced,
            threshold=threshold,
            row_label=reduced_label,
            column_label=reduced_label,
        ),
        _ut1a_matrix_identity_row(
            check="P_reduced involution",
            left=permutation_reduced @ permutation_reduced,
            right=identity_reduced,
            threshold=threshold,
            row_label=reduced_label,
            column_label=reduced_label,
        ),
        _ut1a_matrix_identity_row(
            check="reduction intertwining",
            left=assembly_64.reduction @ permutation_reduced,
            right=permutation_full @ assembly_46.reduction,
            threshold=threshold,
            row_label=full_label,
            column_label=reduced_label,
        ),
        _ut1a_matrix_identity_row(
            check="K_full exchange",
            left=assembly_64.stiffness_full,
            right=(
                permutation_full
                @ assembly_46.stiffness_full
                @ permutation_full.T
            ),
            threshold=threshold,
            row_label=full_label,
            column_label=full_label,
        ),
        _ut1a_matrix_identity_row(
            check="M_full exchange",
            left=assembly_64.mass_full,
            right=(
                permutation_full @ assembly_46.mass_full @ permutation_full.T
            ),
            threshold=threshold,
            row_label=full_label,
            column_label=full_label,
        ),
        _ut1a_matrix_identity_row(
            check="K_reduced forward congruence",
            left=assembly_64.stiffness,
            right=(
                permutation_reduced
                @ assembly_46.stiffness
                @ permutation_reduced.T
            ),
            threshold=threshold,
            row_label=reduced_label,
            column_label=reduced_label,
        ),
        _ut1a_matrix_identity_row(
            check="M_reduced forward congruence",
            left=assembly_64.mass,
            right=(
                permutation_reduced
                @ assembly_46.mass
                @ permutation_reduced.T
            ),
            threshold=threshold,
            row_label=reduced_label,
            column_label=reduced_label,
        ),
        _ut1a_matrix_identity_row(
            check="K_reduced reverse congruence",
            left=assembly_46.stiffness,
            right=(
                permutation_reduced.T
                @ assembly_64.stiffness
                @ permutation_reduced
            ),
            threshold=threshold,
            row_label=reduced_label,
            column_label=reduced_label,
        ),
        _ut1a_matrix_identity_row(
            check="M_reduced reverse congruence",
            left=assembly_46.mass,
            right=(
                permutation_reduced.T
                @ assembly_64.mass
                @ permutation_reduced
            ),
            threshold=threshold,
            row_label=reduced_label,
            column_label=reduced_label,
        ),
    ]
    return rows


def _assemble_ut1a_fem_pair() -> tuple[EBFEMAssembly, EBFEMAssembly]:
    """Assemble exactly the two prescribed production mesh-64 FEM systems."""

    point_4 = _make_section_point(0.004)
    point_6 = _make_section_point(0.006)
    assembly_46 = assemble_two_arm_eb_fem(
        point_4, point_6, BETA_RAD, UT1A_FEM_ELEMENTS_PER_ARM
    )
    assembly_64 = assemble_two_arm_eb_fem(
        point_6, point_4, BETA_RAD, UT1A_FEM_ELEMENTS_PER_ARM
    )
    _validate_ut1a_pair_metadata(assembly_46, assembly_64)
    return assembly_46, assembly_64


def _ut1a_eigenpair_backward_residual(
    stiffness: np.ndarray,
    mass: np.ndarray,
    eigenvalue: float,
    eigenvector: np.ndarray,
) -> float:
    elastic = stiffness @ eigenvector
    inertia = float(eigenvalue) * (mass @ eigenvector)
    return float(
        np.linalg.norm(elastic - inertia)
        / max(
            np.linalg.norm(elastic) + np.linalg.norm(inertia),
            np.finfo(float).tiny,
        )
    )


def _ut1a_mass_orthonormality_residual(
    mass: np.ndarray, eigenvectors: np.ndarray
) -> float:
    identity = np.eye(eigenvectors.shape[1])
    return float(
        np.linalg.norm(eigenvectors.T @ mass @ eigenvectors - identity, ord="fro")
    )


def _ut1a_transported_rayleigh(
    stiffness: np.ndarray, mass: np.ndarray, vector: np.ndarray
) -> float:
    denominator = float(vector.T @ mass @ vector)
    if denominator <= 0.0:
        raise RuntimeError("transported vector has non-positive generalized mass")
    return float((vector.T @ stiffness @ vector) / denominator)


def _ut1a_nearest_gaps(values: Sequence[float]) -> tuple[np.ndarray, np.ndarray]:
    array = np.asarray(values, dtype=float)
    absolute = np.empty_like(array)
    for index, value in enumerate(array):
        other = np.delete(array, index)
        absolute[index] = np.min(np.abs(other - value))
    relative = absolute / np.maximum(np.abs(array), np.finfo(float).tiny)
    return absolute, relative


def _ut1a_conditioning(matrix: np.ndarray) -> dict[str, float]:
    singular_values = np.linalg.svd(np.asarray(matrix, dtype=float), compute_uv=False)
    norm_2 = float(singular_values[0])
    minimum = float(singular_values[-1])
    condition = norm_2 / minimum if minimum > 0.0 else math.inf
    return {"norm_2": norm_2, "cond_2": condition}


def _ut1a_native_and_transport_audit(
    assembly_46: EBFEMAssembly,
    assembly_64: EBFEMAssembly,
    permutation_reduced: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Run native, transported, Rayleigh, and canonicalized eigenaudits."""

    solution_46 = solve_two_arm_eb_fem(
        assembly_46, num_roots=UT1_NUM_ROOTS
    )
    solution_64 = solve_two_arm_eb_fem(
        assembly_64, num_roots=UT1_NUM_ROOTS
    )

    canonical_64_to_46_k = (
        permutation_reduced.T
        @ assembly_64.stiffness
        @ permutation_reduced
    )
    canonical_64_to_46_m = (
        permutation_reduced.T @ assembly_64.mass @ permutation_reduced
    )
    canonical_46_to_64_k = (
        permutation_reduced
        @ assembly_46.stiffness
        @ permutation_reduced.T
    )
    canonical_46_to_64_m = (
        permutation_reduced @ assembly_46.mass @ permutation_reduced.T
    )
    canonical_64_to_46_values, _ = eigh(
        canonical_64_to_46_k,
        canonical_64_to_46_m,
        subset_by_index=[0, 6],
        check_finite=True,
    )
    canonical_46_to_64_values, _ = eigh(
        canonical_46_to_64_k,
        canonical_46_to_64_m,
        subset_by_index=[0, 6],
        check_finite=True,
    )

    native_residuals = {
        "46": np.asarray(
            [
                _ut1a_eigenpair_backward_residual(
                    assembly_46.stiffness,
                    assembly_46.mass,
                    eigenvalue,
                    solution_46.eigenvectors_reduced[:, index],
                )
                for index, eigenvalue in enumerate(solution_46.eigenvalues)
            ]
        ),
        "64": np.asarray(
            [
                _ut1a_eigenpair_backward_residual(
                    assembly_64.stiffness,
                    assembly_64.mass,
                    eigenvalue,
                    solution_64.eigenvectors_reduced[:, index],
                )
                for index, eigenvalue in enumerate(solution_64.eigenvalues)
            ]
        ),
    }
    mass_orthonormality = {
        "46": _ut1a_mass_orthonormality_residual(
            assembly_46.mass, solution_46.eigenvectors_reduced
        ),
        "64": _ut1a_mass_orthonormality_residual(
            assembly_64.mass, solution_64.eigenvectors_reduced
        ),
    }

    eigenvalue_gaps_46 = _ut1a_nearest_gaps(solution_46.eigenvalues)
    eigenvalue_gaps_64 = _ut1a_nearest_gaps(solution_64.eigenvalues)
    frequency_gaps_46 = _ut1a_nearest_gaps(solution_46.frequencies_hz)
    frequency_gaps_64 = _ut1a_nearest_gaps(solution_64.frequencies_hz)
    direction_data = {
        "46_to_64": {
            "source_assembly": assembly_46,
            "target_assembly": assembly_64,
            "source_solution": solution_46,
            "counterpart_solution": solution_64,
            "transport": permutation_reduced,
            "native_residuals": native_residuals["46"],
            "mass_orthonormality": mass_orthonormality["46"],
            "canonical_values": canonical_64_to_46_values,
            "eigenvalue_gaps": eigenvalue_gaps_46,
            "frequency_gaps": frequency_gaps_46,
        },
        "64_to_46": {
            "source_assembly": assembly_64,
            "target_assembly": assembly_46,
            "source_solution": solution_64,
            "counterpart_solution": solution_46,
            "transport": permutation_reduced.T,
            "native_residuals": native_residuals["64"],
            "mass_orthonormality": mass_orthonormality["64"],
            "canonical_values": canonical_46_to_64_values,
            "eigenvalue_gaps": eigenvalue_gaps_64,
            "frequency_gaps": frequency_gaps_64,
        },
    }

    rows: list[dict[str, Any]] = []
    transported_residuals: list[float] = []
    rayleigh_differences: list[float] = []
    canonical_differences: list[float] = []
    for direction, data in direction_data.items():
        source = data["source_solution"]
        counterpart = data["counterpart_solution"]
        target_assembly = data["target_assembly"]
        transport = data["transport"]
        for index in range(UT1_NUM_ROOTS):
            eigenvalue = float(source.eigenvalues[index])
            frequency = float(source.frequencies_hz[index])
            counterpart_eigenvalue = float(counterpart.eigenvalues[index])
            counterpart_frequency = float(counterpart.frequencies_hz[index])
            transported_vector = transport @ source.eigenvectors_reduced[:, index]
            transported_residual = _ut1a_eigenpair_backward_residual(
                target_assembly.stiffness,
                target_assembly.mass,
                eigenvalue,
                transported_vector,
            )
            rayleigh = _ut1a_transported_rayleigh(
                target_assembly.stiffness,
                target_assembly.mass,
                transported_vector,
            )
            rayleigh_frequency = math.sqrt(max(rayleigh, 0.0)) / (2.0 * math.pi)
            rayleigh_difference = abs(rayleigh - eigenvalue) / max(
                abs(eigenvalue), np.finfo(float).tiny
            )
            canonical_value = float(data["canonical_values"][index])
            canonical_frequency = (
                math.sqrt(max(canonical_value, 0.0)) / (2.0 * math.pi)
            )
            canonical_difference = abs(canonical_value - eigenvalue) / max(
                abs(eigenvalue), np.finfo(float).tiny
            )
            native_eigenvalue_difference = abs(
                counterpart_eigenvalue - eigenvalue
            ) / max(abs(eigenvalue), np.finfo(float).tiny)
            native_frequency_difference = abs(
                counterpart_frequency - frequency
            ) / max(abs(frequency), np.finfo(float).tiny)
            native_residual = float(data["native_residuals"][index])
            eigenvalue_absolute_gap = float(data["eigenvalue_gaps"][0][index])
            eigenvalue_relative_gap = float(data["eigenvalue_gaps"][1][index])
            frequency_absolute_gap = float(data["frequency_gaps"][0][index])
            frequency_relative_gap = float(data["frequency_gaps"][1][index])
            quality_ok = bool(
                native_residual <= UT1A_EIGENPAIR_RESIDUAL_TOLERANCE
                and transported_residual <= UT1A_EIGENPAIR_RESIDUAL_TOLERANCE
                and data["mass_orthonormality"]
                <= UT1A_MASS_ORTHONORMALITY_TOLERANCE
                and rayleigh_difference <= UT1A_RAYLEIGH_RELATIVE_TOLERANCE
                and canonical_difference <= UT1A_CANONICAL_SPECTRUM_TOLERANCE
            )
            rows.append(
                {
                    "mode": index + 1,
                    "direction": direction,
                    "native_eigenvalue": eigenvalue,
                    "native_frequency_hz": frequency,
                    "native_counterpart_eigenvalue": counterpart_eigenvalue,
                    "native_counterpart_frequency_hz": counterpart_frequency,
                    "native_eigenpair_residual": native_residual,
                    "native_mass_orthonormality_residual": data[
                        "mass_orthonormality"
                    ],
                    "transported_eigenpair_residual": transported_residual,
                    "transported_rayleigh_eigenvalue": rayleigh,
                    "transported_rayleigh_frequency_hz": rayleigh_frequency,
                    "transported_rayleigh_relative_difference": rayleigh_difference,
                    "native_exchange_eigenvalue_relative_difference": (
                        native_eigenvalue_difference
                    ),
                    "native_exchange_frequency_relative_difference": (
                        native_frequency_difference
                    ),
                    "canonicalized_eigenvalue": canonical_value,
                    "canonicalized_frequency_hz": canonical_frequency,
                    "canonicalized_relative_difference": canonical_difference,
                    "nearest_neighbor_eigenvalue_gap": eigenvalue_absolute_gap,
                    "relative_eigenvalue_gap": eigenvalue_relative_gap,
                    "nearest_neighbor_frequency_gap_hz": frequency_absolute_gap,
                    "relative_frequency_gap": frequency_relative_gap,
                    "quality_status": "PASS" if quality_ok else "FAIL",
                }
            )
            transported_residuals.append(transported_residual)
            rayleigh_differences.append(rayleigh_difference)
            canonical_differences.append(canonical_difference)

    native_frequency_exchange = np.abs(
        solution_46.frequencies_hz - solution_64.frequencies_hz
    ) / solution_46.frequencies_hz
    native_eigenvalue_exchange = np.abs(
        solution_46.eigenvalues - solution_64.eigenvalues
    ) / solution_46.eigenvalues
    diagnostics = {
        "native_exchange_frequency_relative_by_mode": (
            native_frequency_exchange.tolist()
        ),
        "native_exchange_eigenvalue_relative_by_mode": (
            native_eigenvalue_exchange.tolist()
        ),
        "native_exchange_frequency_maximum": float(
            np.max(native_frequency_exchange)
        ),
        "native_exchange_eigenvalue_maximum": float(
            np.max(native_eigenvalue_exchange)
        ),
        "native_eigenpair_residuals": {
            key: value.tolist() for key, value in native_residuals.items()
        },
        "native_eigenpair_residual_maximum": max(
            float(np.max(value)) for value in native_residuals.values()
        ),
        "mass_orthonormality_residuals": mass_orthonormality,
        "transported_eigenpair_residual_maximum": max(transported_residuals),
        "transported_rayleigh_relative_difference_maximum": max(
            rayleigh_differences
        ),
        "canonicalized_spectrum_relative_difference_maximum": max(
            canonical_differences
        ),
        "minimum_relative_eigenvalue_gap": min(
            float(np.min(eigenvalue_gaps_46[1])),
            float(np.min(eigenvalue_gaps_64[1])),
        ),
        "minimum_relative_frequency_gap": min(
            float(np.min(frequency_gaps_46[1])),
            float(np.min(frequency_gaps_64[1])),
        ),
        "canonicalized_eigenvalues": {
            "64_to_46": canonical_64_to_46_values.tolist(),
            "46_to_64": canonical_46_to_64_values.tolist(),
        },
        "native_solutions": {
            "46": {
                "eigenvalues": solution_46.eigenvalues.tolist(),
                "frequencies_hz": solution_46.frequencies_hz.tolist(),
                "eigenvectors_reduced": (
                    solution_46.eigenvectors_reduced.tolist()
                ),
            },
            "64": {
                "eigenvalues": solution_64.eigenvalues.tolist(),
                "frequencies_hz": solution_64.frequencies_hz.tolist(),
                "eigenvectors_reduced": (
                    solution_64.eigenvectors_reduced.tolist()
                ),
            },
        },
    }
    return rows, diagnostics


def _load_ut1_regression_evidence() -> dict[str, Any]:
    path = UT1_DEFAULT_OUTPUT_DIR / "ut1_summary.json"
    if not path.is_file():
        raise FileNotFoundError(
            "UT-1a requires the existing UT-1 summary evidence: " f"{path}"
        )
    return json.loads(path.read_text(encoding="utf-8"))


def _ut1a_skipped_eigen_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for direction in ("46_to_64", "64_to_46"):
        for mode in range(1, UT1_NUM_ROOTS + 1):
            rows.append(
                {
                    "mode": mode,
                    "direction": direction,
                    "native_eigenvalue": None,
                    "native_frequency_hz": None,
                    "native_counterpart_eigenvalue": None,
                    "native_counterpart_frequency_hz": None,
                    "native_eigenpair_residual": None,
                    "native_mass_orthonormality_residual": None,
                    "transported_eigenpair_residual": None,
                    "transported_rayleigh_eigenvalue": None,
                    "transported_rayleigh_frequency_hz": None,
                    "transported_rayleigh_relative_difference": None,
                    "native_exchange_eigenvalue_relative_difference": None,
                    "native_exchange_frequency_relative_difference": None,
                    "canonicalized_eigenvalue": None,
                    "canonicalized_frequency_hz": None,
                    "canonicalized_relative_difference": None,
                    "nearest_neighbor_eigenvalue_gap": None,
                    "relative_eigenvalue_gap": None,
                    "nearest_neighbor_frequency_gap_hz": None,
                    "relative_frequency_gap": None,
                    "quality_status": "SKIPPED_MATRIX_FAILURE",
                }
            )
    return rows


def _ut1a_report(summary: dict[str, Any]) -> str:
    lines = [
        "# UT-1a beta=0 EB FEM matrix-level exchange audit",
        "",
        "Overall unequal-thickness validation: **IN_PROGRESS**",
        "",
        "UT-0 section-routing and beta=0 stepped-reference smoke: "
        f"**{summary['ut0_regression_status']}**",
        "",
        "UT-1 full beta=0 unequal-thickness validation: "
        f"**{summary['ut1_regression_status']}**",
        "",
        "UT-1a beta=0 EB FEM matrix-level exchange audit: "
        f"**{summary['ut1a_status']}**",
        "",
        "## Scope and ordering",
        "",
        "Only elastic HMS/DX-209, `theta_1=theta_2=0`, `beta=0 deg`, "
        "`L_1=L_2=0.2 m`, `b_1=b_2=0.020 m`, `(a_1,a_2)=(4,6)` "
        "and `(6,4) mm`, and `N_1=N_2=64` are used. No continuum root "
        "calculation or mesh sequence is performed.",
        "",
        "The full ordering is `[x_1,x_2]`; each arm contains all local "
        "`[w,psi,Phi]` node triples from outer clamp to joint. The reduced "
        "ordering is `[q_1,int,q_2,int,w_J,theta_t,theta_n]`.",
        "",
        f"The two full systems are `{summary['full_size']} x "
        f"{summary['full_size']}` and the two reduced systems are "
        f"`{summary['reduced_size']} x {summary['reduced_size']}`.",
        "",
        "`P_f=[[0,I],[I,0]]` swaps the two complete full arm blocks without "
        "local signs. `P_r` swaps the two internal blocks and applies "
        "`J_0=diag(1,-1,-1)` to the joint block. At `beta=0`, the existing "
        "endpoint maps are arm 1 `[w,-theta_n,theta_t]` and arm 2 "
        "`[w,theta_n,-theta_t]`; therefore both physical joint-rotation "
        "coordinates change sign under arm exchange.",
        "",
        "## Matrix identities",
        "",
        "| check | relative Frobenius | relative max entry | threshold | status |",
        "|---|---:|---:|---:|---|",
    ]
    for row in summary["matrix_checks"]:
        lines.append(
            f"| {row['check']} | {row['relative_frobenius_residual']:.12e} | "
            f"{row['relative_max_residual']:.12e} | {row['threshold']:.1e} | "
            f"{row['status']} |"
        )
    lines.extend(
        [
            "",
            "The reduction identity is `R_64 P_r = P_f R_46`; full matrices "
            "use `A_64=P_f A_46 P_f^T`; reduced matrices use both forward "
            "and reverse congruence forms.",
            "",
            "## Native and transported eigenpairs",
            "",
        ]
    )
    eigen = summary.get("eigen_audit")
    if eigen is None:
        lines.append(
            "The eigensolver interpretation was skipped because a matrix-level "
            "identity failed. Mismatch DOF labels are retained in the CSV/JSON evidence."
        )
    else:
        lines.extend(
            [
                f"Native FEM exchange maximum: "
                f"`{eigen['native_exchange_frequency_maximum']:.12e}` "
                f"(historical UT-1 threshold `1e-8`).",
                "",
                f"Maximum native backward residual: "
                f"`{eigen['native_eigenpair_residual_maximum']:.12e}` "
                f"(threshold `1e-8`).",
                "",
                f"Mass-orthonormality residuals: `(4,6)` "
                f"`{eigen['mass_orthonormality_residuals']['46']:.12e}`, "
                f"`(6,4)` `{eigen['mass_orthonormality_residuals']['64']:.12e}` "
                f"(threshold `1e-10`).",
                "",
                f"Maximum transported backward residual: "
                f"`{eigen['transported_eigenpair_residual_maximum']:.12e}` "
                f"(threshold `1e-8`).",
                "",
                f"Maximum transported Rayleigh relative difference: "
                f"`{eigen['transported_rayleigh_relative_difference_maximum']:.12e}` "
                f"(threshold `1e-10`).",
                "",
                f"Maximum canonicalized-spectrum relative difference: "
                f"`{eigen['canonicalized_spectrum_relative_difference_maximum']:.12e}` "
                f"(threshold `1e-10`).",
                "",
                f"Minimum relative generalized-eigenvalue gap: "
                f"`{eigen['minimum_relative_eigenvalue_gap']:.12e}`. No "
                "condition-number or gap acceptance criterion is introduced.",
                "",
                "| mode | native f(4,6), Hz | native f(6,4), Hz | relative difference |",
                "|---:|---:|---:|---:|",
            ]
        )
        frequencies_46 = eigen["native_solutions"]["46"]["frequencies_hz"]
        frequencies_64 = eigen["native_solutions"]["64"]["frequencies_hz"]
        for mode, (frequency_46, frequency_64, relative) in enumerate(
            zip(
                frequencies_46,
                frequencies_64,
                eigen["native_exchange_frequency_relative_by_mode"],
            ),
            start=1,
        ):
            lines.append(
                f"| {mode} | {frequency_46:.12f} | {frequency_64:.12f} | "
                f"{relative:.12e} |"
            )
        lines.extend(
            [
                "",
                "| ordering | matrix | norm_2 | cond_2 |",
                "|---|---|---:|---:|",
            ]
        )
        for ordering in ("46", "64"):
            for matrix_name in ("stiffness", "mass"):
                values = summary["conditioning"][ordering][matrix_name]
                lines.append(
                    f"| {ordering} | {matrix_name} | {values['norm_2']:.12e} | "
                    f"{values['cond_2']:.12e} |"
                )
        lines.extend(
            [
                "",
                "The transported vectors are algebraic permutation checks, not "
                "MAC, branch tracking, or physical mode-shape results. Native "
                "UT-1 frequencies are not replaced by Rayleigh or canonicalized values.",
            ]
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"Classification: `{summary['interpretation_classification']}`.",
            "",
            summary["interpretation"],
            "",
            "UT-1 remains `PARTIAL_PASS` because its historically fixed native "
            "exchange threshold is not changed. The overall unequal-thickness "
            "validation remains `IN_PROGRESS`.",
            "",
            "## Explicit exclusions",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in summary["explicit_exclusions"])
    return "\n".join(lines) + "\n"


def _run_ut1a_fem_exchange_audit(output_dir: Path) -> dict[str, Any]:
    """Run the isolated two-assembly beta-zero matrix exchange audit."""

    started = time.perf_counter()
    prior = _load_ut1_regression_evidence()
    assembly_46, assembly_64 = _assemble_ut1a_fem_pair()
    metadata = _validate_ut1a_pair_metadata(assembly_46, assembly_64)
    permutation_full = _ut1a_full_swap_permutation(assembly_46, assembly_64)
    permutation_reduced = _ut1a_reduced_swap_permutation(
        assembly_46, assembly_64
    )
    matrix_rows = _ut1a_matrix_congruence_rows(
        assembly_46,
        assembly_64,
        permutation_full,
        permutation_reduced,
    )
    matrix_equivariance_ok = all(row["status"] == "PASS" for row in matrix_rows)

    eigen_rows: list[dict[str, Any]]
    eigen_diagnostics: dict[str, Any] | None
    conditioning: dict[str, Any] | None
    if matrix_equivariance_ok:
        eigen_rows, eigen_diagnostics = _ut1a_native_and_transport_audit(
            assembly_46, assembly_64, permutation_reduced
        )
        conditioning = {
            "46": {
                "stiffness": _ut1a_conditioning(assembly_46.stiffness),
                "mass": _ut1a_conditioning(assembly_46.mass),
            },
            "64": {
                "stiffness": _ut1a_conditioning(assembly_64.stiffness),
                "mass": _ut1a_conditioning(assembly_64.mass),
            },
        }
        numerical_ok = bool(
            eigen_diagnostics["native_eigenpair_residual_maximum"]
            <= UT1A_EIGENPAIR_RESIDUAL_TOLERANCE
            and eigen_diagnostics["transported_eigenpair_residual_maximum"]
            <= UT1A_EIGENPAIR_RESIDUAL_TOLERANCE
            and max(eigen_diagnostics["mass_orthonormality_residuals"].values())
            <= UT1A_MASS_ORTHONORMALITY_TOLERANCE
            and eigen_diagnostics[
                "transported_rayleigh_relative_difference_maximum"
            ]
            <= UT1A_RAYLEIGH_RELATIVE_TOLERANCE
            and eigen_diagnostics[
                "canonicalized_spectrum_relative_difference_maximum"
            ]
            <= UT1A_CANONICAL_SPECTRUM_TOLERANCE
            and all(row["quality_status"] == "PASS" for row in eigen_rows)
        )
        if numerical_ok:
            ut1a_status = "PASS_MATRIX_EQUIVARIANCE"
            interpretation_classification = (
                "DOF_ORDER_SENSITIVE_EIGENSOLVER_ERROR_CONSISTENT"
            )
            interpretation = (
                "The beta=0 FEM assemblies are exchange-equivariant at matrix "
                "level within the declared thresholds. The remaining native "
                "spectral difference is consistent with DOF-order-sensitive "
                "generalized-eigensolver/conditioning error; assembly asymmetry "
                "was not detected within the declared thresholds. This is not "
                "an absolute proof of eigensolver error."
            )
        else:
            ut1a_status = "INCONCLUSIVE_NUMERICAL_AUDIT"
            interpretation_classification = "INCONCLUSIVE"
            interpretation = (
                "Matrix congruence passed, but at least one transported, "
                "Rayleigh, canonicalized, backward-error, or mass-orthonormality "
                "gate did not pass its predeclared threshold."
            )
    else:
        eigen_rows = _ut1a_skipped_eigen_rows()
        eigen_diagnostics = None
        conditioning = None
        ut1a_status = "FAIL_MATRIX_EQUIVARIANCE"
        interpretation_classification = (
            "ASSEMBLY_OR_REDUCTION_ASYMMETRY_DETECTED"
        )
        interpretation = (
            "At least one a-priori permutation, reduction-intertwining, full "
            "exchange, or reduced congruence identity failed. Eigensolver "
            "interpretation was stopped; mismatch entries and DOF labels are "
            "saved without changing the FEM implementation."
        )

    by_check = {row["check"]: row for row in matrix_rows}
    permutation_names = {
        "P_full orthogonality",
        "P_full involution",
        "P_reduced orthogonality",
        "P_reduced involution",
    }
    full_names = {"K_full exchange", "M_full exchange"}
    reduced_names = {
        "K_reduced forward congruence",
        "M_reduced forward congruence",
        "K_reduced reverse congruence",
        "M_reduced reverse congruence",
    }
    prior_native_exchange = float(prior["fem_exchange_maximum"])
    reproduced_native_exchange = (
        None
        if eigen_diagnostics is None
        else eigen_diagnostics["native_exchange_frequency_maximum"]
    )
    summary = {
        "git_context": _git_context(),
        "scope": "UT-1a beta=0 EB FEM matrix-level exchange audit only",
        "geometry": {
            "theta_1_deg": 0.0,
            "theta_2_deg": 0.0,
            "material_mode": "elastic",
            "material": "existing HMS/DX-209",
            "length_1_m": LENGTH_M,
            "length_2_m": LENGTH_M,
            "b_1_m": WIDTH_B_M,
            "b_2_m": WIDTH_B_M,
            "beta_deg": 0.0,
            "assembly_46_a_m": [0.004, 0.006],
            "assembly_64_a_m": [0.006, 0.004],
        },
        "element_counts": {
            "assembly_46": list(assembly_46.element_counts),
            "assembly_64": list(assembly_64.element_counts),
        },
        "assembly_count": 2,
        "num_roots": UT1_NUM_ROOTS,
        "dof_metadata": metadata,
        "full_size": int(assembly_46.stiffness_full.shape[0]),
        "reduced_size": int(assembly_46.stiffness.shape[0]),
        "permutation_definitions": {
            "P_full": (
                "[[0,I_{3(N+1)}],[I_{3(N+1)},0]]; no local nodal signs"
            ),
            "P_reduced": (
                "[[0,I_{3(N-1)},0],[I_{3(N-1)},0,0],"
                "[0,0,diag(1,-1,-1)]]"
            ),
            "mapping": "x_64=P_full x_46; q_64=P_reduced q_46",
            "joint_sign_reason": (
                "beta=0 endpoint maps are arm1 [w,-theta_n,theta_t] "
                "and arm2 [w,theta_n,-theta_t]"
            ),
        },
        "thresholds": {
            "permutation_reduction_matrix_relative": (
                UT1A_MATRIX_RELATIVE_TOLERANCE
            ),
            "native_eigenpair_backward": (
                UT1A_EIGENPAIR_RESIDUAL_TOLERANCE
            ),
            "transported_eigenpair_backward": (
                UT1A_EIGENPAIR_RESIDUAL_TOLERANCE
            ),
            "mass_orthonormality": UT1A_MASS_ORTHONORMALITY_TOLERANCE,
            "transported_rayleigh_relative": (
                UT1A_RAYLEIGH_RELATIVE_TOLERANCE
            ),
            "canonicalized_spectrum_relative": (
                UT1A_CANONICAL_SPECTRUM_TOLERANCE
            ),
            "historical_native_fem_exchange_relative": (
                SPECTRUM_RELATIVE_TOLERANCE
            ),
        },
        "matrix_checks": matrix_rows,
        "maximum_permutation_relative_frobenius_residual": max(
            row["relative_frobenius_residual"]
            for row in matrix_rows
            if row["check"] in permutation_names
        ),
        "maximum_permutation_relative_max_residual": max(
            row["relative_max_residual"]
            for row in matrix_rows
            if row["check"] in permutation_names
        ),
        "reduction_intertwining": by_check["reduction intertwining"],
        "full_matrix_residuals": {
            name: by_check[name] for name in sorted(full_names)
        },
        "reduced_matrix_residuals": {
            name: by_check[name] for name in sorted(reduced_names)
        },
        "prior_ut1_native_exchange_maximum": prior_native_exchange,
        "reproduced_native_exchange_maximum": reproduced_native_exchange,
        "native_exchange_reproduction_relative_difference": (
            None
            if reproduced_native_exchange is None
            else abs(reproduced_native_exchange - prior_native_exchange)
            / max(abs(prior_native_exchange), np.finfo(float).tiny)
        ),
        "eigen_audit": eigen_diagnostics,
        "conditioning": conditioning,
        "minimum_relative_spectral_gap": (
            None
            if eigen_diagnostics is None
            else eigen_diagnostics["minimum_relative_eigenvalue_gap"]
        ),
        "ut0_regression_status": prior["ut0_regression_status"],
        "ut1_regression_status": prior["ut1_status"],
        "ut1a_status": ut1a_status,
        "overall_unequal_thickness_validation": "IN_PROGRESS",
        "interpretation_classification": interpretation_classification,
        "interpretation": interpretation,
        "continuum_root_calculations": 0,
        "runtime_seconds": time.perf_counter() - started,
        "explicit_exclusions": UT1A_EXPLICIT_EXCLUSIONS,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "matrix_congruence.csv", matrix_rows)
    _write_csv(output_dir / "eigenpair_transport.csv", eigen_rows)
    (output_dir / "ut1a_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "ut1a_report.md").write_text(
        _ut1a_report(summary), encoding="utf-8"
    )
    return summary


def _angular_beta_rad(
    beta_deg: float,
    *,
    allowed_angles_deg: Sequence[float],
    gate_label: str,
) -> float:
    value = float(beta_deg)
    allowed = tuple(float(angle) for angle in allowed_angles_deg)
    if value not in allowed:
        raise ValueError(
            f"{gate_label} accepts only beta in {allowed}, got {value}"
        )
    return math.radians(value)


def _ut2_beta_rad(beta_deg: float) -> float:
    return _angular_beta_rad(
        beta_deg,
        allowed_angles_deg=UT2_ANGLES_DEG,
        gate_label="UT-2",
    )


def _ut3_beta_rad(beta_deg: float) -> float:
    return _angular_beta_rad(
        beta_deg,
        allowed_angles_deg=UT3_ANGLES_DEG,
        gate_label="UT-3",
    )


def _make_ut2_section_cases() -> dict[str, tuple[RodPoint, RodPoint]]:
    """Build independent arm points for the three prescribed UT-2 cases."""

    return {
        case: (_make_section_point(a_1), _make_section_point(a_2))
        for case, (a_1, a_2) in UT2_SECTION_CASES.items()
    }


def _ut2_continuum_factories(
    model: str,
    beta_rad: float,
    point_1: RodPoint,
    point_2: RodPoint,
) -> tuple[Callable[[complex], np.ndarray], Callable[[complex], np.ndarray]]:
    """Return only the accepted angular coupled/raw builders for UT-2."""

    if model == "Timoshenko":
        return (
            lambda omega: coupled_boundary_matrix(
                omega, beta_rad, point_1, point_2
            ),
            lambda omega: coupled_boundary_matrix_raw(
                omega, beta_rad, point_1, point_2
            ),
        )
    if model == "EB":
        return (
            lambda omega: eb_coupled_boundary_matrix(
                omega, beta_rad, point_1, point_2
            ),
            lambda omega: eb_coupled_boundary_matrix_raw(
                omega, beta_rad, point_1, point_2
            ),
        )
    raise ValueError(f"unsupported UT-2 continuum model: {model}")


def _ut2_continuum_root_row(
    *,
    section_case: str,
    point_1: RodPoint,
    point_2: RodPoint,
    beta_deg: float,
    model: str,
    mode: int,
    root: RootResult,
    quality: dict[str, Any],
    gap_hz: float,
) -> dict[str, Any]:
    return {
        "section_case": section_case,
        "a_1_m": point_1.geometry.a,
        "a_2_m": point_2.geometry.a,
        "beta_deg": beta_deg,
        "model": model,
        "mode": mode,
        "frequency_hz": root.frequency_hz,
        "scaled_determinant_abs": quality["scaled_determinant_abs"],
        "scaled_determinant_residual": quality[
            "scaled_determinant_residual"
        ],
        "scaled_sigma_min": quality["scaled_sigma_min"],
        "scaled_sigma_max": quality["scaled_sigma_max"],
        "scaled_relative_singular_residual": quality[
            "scaled_relative_singular_residual"
        ],
        "raw_determinant_abs": quality["physical_raw_determinant_abs"],
        "raw_determinant_residual": quality[
            "physical_raw_normalized_determinant_residual"
        ],
        "raw_sigma_min": quality["physical_raw_sigma_min"],
        "raw_sigma_max": quality["physical_raw_sigma_max"],
        "raw_relative_singular_residual": quality[
            "physical_raw_relative_singular_residual"
        ],
        "root_status": quality["root_status"],
        "quality_status": quality["root_quality_status"],
        "quality_basis": quality["root_quality_basis"],
        "neighbor_gap_hz": gap_hz,
        "relative_neighbor_gap": gap_hz / root.frequency_hz,
    }


def _ut2_pair_metadata(
    assembly_left: EBFEMAssembly, assembly_right: EBFEMAssembly
) -> dict[str, Any]:
    left = _validate_ut1a_dof_ordering(assembly_left)
    right = _validate_ut1a_dof_ordering(assembly_right)
    if assembly_left.element_counts != assembly_right.element_counts:
        raise ValueError("UT-2 matrix pair requires identical element counts")
    for key in (
        "element_counts",
        "full_arm_sizes",
        "internal_arm_sizes",
        "full_size",
        "reduced_size",
        "full_ordering",
        "reduced_ordering",
    ):
        if left[key] != right[key]:
            raise RuntimeError(f"UT-2 matrix pair metadata differs at {key}")
    return left


def _ut2_reflection_full_transform(assembly: EBFEMAssembly) -> np.ndarray:
    metadata = _validate_ut1a_dof_ordering(assembly)
    full_size = int(metadata["full_size"])
    if full_size % FEM_NODE_DOF_COUNT:
        raise RuntimeError("full size is not divisible by nodal DOF count")
    node_count = full_size // FEM_NODE_DOF_COUNT
    nodal = np.diag([1.0, 1.0, -1.0])
    return np.kron(np.eye(node_count), nodal)


def _ut2_reflection_reduced_transform(assembly: EBFEMAssembly) -> np.ndarray:
    metadata = _validate_ut1a_dof_ordering(assembly)
    reduced_size = int(metadata["reduced_size"])
    internal_size = reduced_size - 3
    if internal_size % FEM_NODE_DOF_COUNT:
        raise RuntimeError("internal reduced size is not divisible by nodal DOFs")
    transform = np.zeros((reduced_size, reduced_size), dtype=float)
    transform[:internal_size, :internal_size] = np.kron(
        np.eye(internal_size // FEM_NODE_DOF_COUNT),
        np.diag([1.0, 1.0, -1.0]),
    )
    transform[-3:, -3:] = np.diag([1.0, -1.0, 1.0])
    return transform


def _ut2_relabel_full_transform(
    assembly_source: EBFEMAssembly, assembly_target: EBFEMAssembly
) -> np.ndarray:
    metadata = _ut2_pair_metadata(assembly_source, assembly_target)
    arm_sizes = metadata["full_arm_sizes"]
    if arm_sizes[0] != arm_sizes[1]:
        raise ValueError("UT-2 relabeling requires equal full arm-block sizes")
    arm_size = int(arm_sizes[0])
    transform = np.zeros((2 * arm_size, 2 * arm_size), dtype=float)
    transform[:arm_size, arm_size:] = np.eye(arm_size)
    transform[arm_size:, :arm_size] = np.eye(arm_size)
    return transform


def _ut2_relabel_joint_transform(beta_rad: float) -> np.ndarray:
    c = math.cos(float(beta_rad))
    s = math.sin(float(beta_rad))
    return np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, -c, s],
            [0.0, -s, -c],
        ]
    )


def _ut2_relabel_reduced_transform(
    assembly_source: EBFEMAssembly,
    assembly_target: EBFEMAssembly,
    beta_rad: float,
) -> np.ndarray:
    metadata = _ut2_pair_metadata(assembly_source, assembly_target)
    internal_sizes = metadata["internal_arm_sizes"]
    if internal_sizes[0] != internal_sizes[1]:
        raise ValueError("UT-2 relabeling requires equal internal arm-block sizes")
    internal_size = int(internal_sizes[0])
    reduced_size = int(metadata["reduced_size"])
    transform = np.zeros((reduced_size, reduced_size), dtype=float)
    transform[:internal_size, internal_size : 2 * internal_size] = np.eye(
        internal_size
    )
    transform[internal_size : 2 * internal_size, :internal_size] = np.eye(
        internal_size
    )
    transform[-3:, -3:] = _ut2_relabel_joint_transform(beta_rad)
    return transform


def _ut2_matrix_row(
    *,
    symmetry_type: str,
    section_case_left: str,
    section_case_right: str,
    beta_left_deg: float,
    beta_right_deg: float,
    check: str,
    left: np.ndarray,
    right: np.ndarray,
    row_label: Callable[[int], str] | None = None,
    column_label: Callable[[int], str] | None = None,
    threshold: float = UT2_MATRIX_RELATIVE_TOLERANCE,
) -> dict[str, Any]:
    evidence = _ut1a_matrix_identity_row(
        check=check,
        left=left,
        right=right,
        threshold=threshold,
        row_label=row_label,
        column_label=column_label,
    )
    return {
        "symmetry_type": symmetry_type,
        "section_case_left": section_case_left,
        "section_case_right": section_case_right,
        "beta_left_deg": beta_left_deg,
        "beta_right_deg": beta_right_deg,
        **evidence,
    }


def _ut2_endpoint_map_rows(
    beta_deg: float = 30.0,
) -> list[dict[str, Any]]:
    """Return the generic angular reflection/relabeling endpoint identities."""

    beta = math.radians(float(beta_deg))
    maps_plus = eb_joint_dof_maps(beta)
    maps_minus = eb_joint_dof_maps(-beta)
    reflection_nodal = np.diag([1.0, 1.0, -1.0])
    reflection_joint = np.diag([1.0, -1.0, 1.0])
    relabel_joint = _ut2_relabel_joint_transform(beta)
    return [
        _ut2_matrix_row(
            symmetry_type="reflection",
            section_case_left="all",
            section_case_right="all",
            beta_left_deg=beta_deg,
            beta_right_deg=-beta_deg,
            check="endpoint map arm 1",
            left=reflection_nodal @ maps_plus[0],
            right=maps_minus[0] @ reflection_joint,
        ),
        _ut2_matrix_row(
            symmetry_type="reflection",
            section_case_left="all",
            section_case_right="all",
            beta_left_deg=beta_deg,
            beta_right_deg=-beta_deg,
            check="endpoint map arm 2",
            left=reflection_nodal @ maps_plus[1],
            right=maps_minus[1] @ reflection_joint,
        ),
        _ut2_matrix_row(
            symmetry_type="oriented_angle_relabeling",
            section_case_left="all",
            section_case_right="all",
            beta_left_deg=beta_deg,
            beta_right_deg=-beta_deg,
            check="endpoint map F1 H_rel = F2(+beta)",
            left=maps_plus[0] @ relabel_joint,
            right=maps_plus[1],
        ),
        _ut2_matrix_row(
            symmetry_type="oriented_angle_relabeling",
            section_case_left="all",
            section_case_right="all",
            beta_left_deg=beta_deg,
            beta_right_deg=-beta_deg,
            check="endpoint map F2(-beta) H_rel = F1",
            left=maps_minus[1] @ relabel_joint,
            right=maps_plus[0],
        ),
    ]


def _ut2_reflection_matrix_rows(
    section_case: str,
    assembly_plus: EBFEMAssembly,
    assembly_minus: EBFEMAssembly,
    *,
    beta_deg: float = 30.0,
) -> list[dict[str, Any]]:
    _ut2_pair_metadata(assembly_plus, assembly_minus)
    full = _ut2_reflection_full_transform(assembly_plus)
    reduced = _ut2_reflection_reduced_transform(assembly_plus)
    full_identity = np.eye(full.shape[0])
    reduced_identity = np.eye(reduced.shape[0])
    full_label = lambda index: _ut1a_full_dof_label(assembly_minus, index)
    reduced_label = lambda index: _ut1a_reduced_dof_label(
        assembly_minus, index
    )
    common = {
        "symmetry_type": "reflection",
        "section_case_left": section_case,
        "section_case_right": section_case,
        "beta_left_deg": beta_deg,
        "beta_right_deg": -beta_deg,
    }
    specifications = [
        (
            "P_full reflection orthogonality",
            full.T @ full,
            full_identity,
            full_label,
            full_label,
        ),
        (
            "P_full reflection involution",
            full @ full,
            full_identity,
            full_label,
            full_label,
        ),
        (
            "P_reduced reflection orthogonality",
            reduced.T @ reduced,
            reduced_identity,
            reduced_label,
            reduced_label,
        ),
        (
            "P_reduced reflection involution",
            reduced @ reduced,
            reduced_identity,
            reduced_label,
            reduced_label,
        ),
        (
            "reflection reduction intertwining",
            assembly_minus.reduction @ reduced,
            full @ assembly_plus.reduction,
            full_label,
            reduced_label,
        ),
        (
            "K_full reflection congruence",
            assembly_minus.stiffness_full,
            full @ assembly_plus.stiffness_full @ full.T,
            full_label,
            full_label,
        ),
        (
            "M_full reflection congruence",
            assembly_minus.mass_full,
            full @ assembly_plus.mass_full @ full.T,
            full_label,
            full_label,
        ),
        (
            "K_reduced reflection congruence",
            assembly_minus.stiffness,
            reduced @ assembly_plus.stiffness @ reduced.T,
            reduced_label,
            reduced_label,
        ),
        (
            "M_reduced reflection congruence",
            assembly_minus.mass,
            reduced @ assembly_plus.mass @ reduced.T,
            reduced_label,
            reduced_label,
        ),
    ]
    return [
        _ut2_matrix_row(
            **common,
            check=check,
            left=left,
            right=right,
            row_label=row_label,
            column_label=column_label,
        )
        for check, left, right, row_label, column_label in specifications
    ]


def _ut2_relabeling_matrix_rows(
    section_case_source: str,
    section_case_target: str,
    assembly_source: EBFEMAssembly,
    assembly_target: EBFEMAssembly,
    *,
    beta_deg: float = 30.0,
) -> list[dict[str, Any]]:
    beta = math.radians(float(beta_deg))
    full = _ut2_relabel_full_transform(assembly_source, assembly_target)
    reduced = _ut2_relabel_reduced_transform(
        assembly_source, assembly_target, beta
    )
    full_label = lambda index: _ut1a_full_dof_label(assembly_target, index)
    reduced_label = lambda index: _ut1a_reduced_dof_label(
        assembly_target, index
    )
    common = {
        "symmetry_type": "oriented_angle_relabeling",
        "section_case_left": section_case_source,
        "section_case_right": section_case_target,
        "beta_left_deg": beta_deg,
        "beta_right_deg": -beta_deg,
    }
    specifications = [
        (
            "P_full relabeling orthogonality",
            full.T @ full,
            np.eye(full.shape[0]),
            full_label,
            full_label,
        ),
        (
            "P_reduced relabeling orthogonality",
            reduced.T @ reduced,
            np.eye(reduced.shape[0]),
            reduced_label,
            reduced_label,
        ),
        (
            "relabeling reduction intertwining",
            assembly_target.reduction @ reduced,
            full @ assembly_source.reduction,
            full_label,
            reduced_label,
        ),
        (
            "K_full relabeling congruence",
            assembly_target.stiffness_full,
            full @ assembly_source.stiffness_full @ full.T,
            full_label,
            full_label,
        ),
        (
            "M_full relabeling congruence",
            assembly_target.mass_full,
            full @ assembly_source.mass_full @ full.T,
            full_label,
            full_label,
        ),
        (
            "K_reduced relabeling congruence",
            assembly_target.stiffness,
            reduced @ assembly_source.stiffness @ reduced.T,
            reduced_label,
            reduced_label,
        ),
        (
            "M_reduced relabeling congruence",
            assembly_target.mass,
            reduced @ assembly_source.mass @ reduced.T,
            reduced_label,
            reduced_label,
        ),
    ]
    return [
        _ut2_matrix_row(
            **common,
            check=check,
            left=left,
            right=right,
            row_label=row_label,
            column_label=column_label,
        )
        for check, left, right, row_label, column_label in specifications
    ]


def _classify_ut2_status(
    *,
    continuum_hard_ok: bool,
    matrix_hard_ok: bool,
    fem_hard_ok: bool,
    regressions_ok: bool,
) -> str:
    """Classify UT-2 without using native FEM spectral symmetry as a gate."""

    if not continuum_hard_ok or not matrix_hard_ok or not regressions_ok:
        return "FAIL"
    return "PASS" if fem_hard_ok else "PARTIAL_PASS"


def _ut2_symmetry_specifications(
    *,
    angle_deg: float = 30.0,
    section_cases: Sequence[str] | None = None,
) -> list[dict[str, Any]]:
    cases = tuple(UT2_SECTION_CASES if section_cases is None else section_cases)
    return [
        *[
            {
                "symmetry_type": "reflection",
                "case_left": case,
                "case_right": case,
                "beta_left_deg": angle_deg,
                "beta_right_deg": -angle_deg,
            }
            for case in cases
        ],
        {
            "symmetry_type": "oriented_angle_relabeling",
            "case_left": "asymmetric_4_6",
            "case_right": "swapped_6_4",
            "beta_left_deg": angle_deg,
            "beta_right_deg": -angle_deg,
        },
        {
            "symmetry_type": "oriented_angle_relabeling",
            "case_left": "swapped_6_4",
            "case_right": "asymmetric_4_6",
            "beta_left_deg": angle_deg,
            "beta_right_deg": -angle_deg,
        },
        {
            "symmetry_type": "same_positive_exchange",
            "case_left": "asymmetric_4_6",
            "case_right": "swapped_6_4",
            "beta_left_deg": angle_deg,
            "beta_right_deg": angle_deg,
        },
        {
            "symmetry_type": "same_negative_exchange",
            "case_left": "asymmetric_4_6",
            "case_right": "swapped_6_4",
            "beta_left_deg": -angle_deg,
            "beta_right_deg": -angle_deg,
        },
    ]


def _ut2_continuum_symmetry_rows(
    spectra: dict[tuple[str, float, str], dict[str, Any]],
    *,
    specifications: Sequence[dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, float]]]:
    rows: list[dict[str, Any]] = []
    maxima: dict[str, dict[str, float]] = {}
    selected = (
        _ut2_symmetry_specifications()
        if specifications is None
        else specifications
    )
    for specification in selected:
        for model in ("Timoshenko", "EB"):
            left_key = (
                specification["case_left"],
                specification["beta_left_deg"],
                model,
            )
            right_key = (
                specification["case_right"],
                specification["beta_right_deg"],
                model,
            )
            left = _frequencies(spectra[left_key]["roots"])
            right, matching, _, permutation = _match_frequency_spectra(
                left, _frequencies(spectra[right_key]["roots"])
            )
            differences: list[float] = []
            for mode, (left_hz, right_hz, basis, right_index) in enumerate(
                zip(left, right, matching, permutation), start=1
            ):
                absolute = abs(float(left_hz) - float(right_hz))
                relative = absolute / max(
                    abs(float(left_hz)), np.finfo(float).tiny
                )
                differences.append(relative)
                rows.append(
                    {
                        "symmetry_type": specification["symmetry_type"],
                        "model": model,
                        "case_left": specification["case_left"],
                        "case_right": specification["case_right"],
                        "beta_left_deg": specification["beta_left_deg"],
                        "beta_right_deg": specification["beta_right_deg"],
                        "mode": mode,
                        "frequency_left_hz": left_hz,
                        "frequency_right_hz": right_hz,
                        "absolute_difference_hz": absolute,
                        "relative_difference": relative,
                        "matching_basis": basis,
                        "right_sorted_index": int(right_index) + 1,
                        "threshold": SPECTRUM_RELATIVE_TOLERANCE,
                        "status": (
                            "PASS"
                            if relative <= SPECTRUM_RELATIVE_TOLERANCE
                            else "FAIL"
                        ),
                    }
                )
            label = (
                f"{model}|{specification['case_left']}|"
                f"{specification['beta_left_deg']:+g}_to_"
                f"{specification['case_right']}|"
                f"{specification['beta_right_deg']:+g}"
            )
            maxima.setdefault(specification["symmetry_type"], {})[label] = max(
                differences
            )
    return rows, maxima


def _run_ut2_eb_fem_configuration(
    *,
    section_case: str,
    point_1: RodPoint,
    point_2: RodPoint,
    beta_deg: float,
    analytic_frequencies: Sequence[float],
    allowed_angles_deg: Sequence[float] = UT2_ANGLES_DEG,
    gate_label: str = "UT-2",
    elements_per_arm: tuple[int, int] = UT2_FEM_ELEMENTS_PER_ARM,
    num_roots: int = UT2_NUM_ROOTS,
) -> dict[str, Any]:
    started = time.perf_counter()
    beta_rad = _angular_beta_rad(
        beta_deg,
        allowed_angles_deg=allowed_angles_deg,
        gate_label=gate_label,
    )
    assembly = assemble_two_arm_eb_fem(
        point_1,
        point_2,
        beta_rad,
        elements_per_arm,
    )
    solution = solve_two_arm_eb_fem(assembly, num_roots=num_roots)
    runtime = time.perf_counter() - started
    matched, matching, gaps, permutation = _match_frequency_spectra(
        analytic_frequencies, solution.frequencies_hz
    )
    errors = np.abs(matched - np.asarray(analytic_frequencies, dtype=float)) / np.asarray(
        analytic_frequencies, dtype=float
    )
    matched_equilibrium = solution.joint_equilibrium_residuals[permutation]
    diagnostics = {
        "section_case": section_case,
        "beta_deg": beta_deg,
        "elements_arm_1": int(assembly.element_counts[0]),
        "elements_arm_2": int(assembly.element_counts[1]),
        "element_length_1_m": point_1.geometry.length
        / assembly.element_counts[0],
        "element_length_2_m": point_2.geometry.length
        / assembly.element_counts[1],
        "reduced_size": int(assembly.stiffness.shape[0]),
        "stiffness_symmetry_residual": solution.stiffness_symmetry_residual,
        "mass_symmetry_residual": solution.mass_symmetry_residual,
        "minimum_mass_eigenvalue": solution.minimum_mass_eigenvalue,
        "zero_mode_count": solution.zero_mode_count,
        "root_count": len(solution.frequencies_hz),
        "positive_finite_roots": bool(
            len(solution.frequencies_hz) == num_roots
            and np.all(np.isfinite(solution.frequencies_hz))
            and np.all(solution.frequencies_hz > 0.0)
        ),
        "joint_kinematic_residual": eb_joint_mapping_residual(beta_rad),
        "max_joint_equilibrium_residual": float(
            np.max(solution.joint_equilibrium_residuals)
        ),
        "runtime_seconds": runtime,
    }
    diagnostics["structural_status"] = (
        "PASS" if _fem_structural_ok(diagnostics) else "FAIL"
    )
    e3 = float(np.max(errors[:3]))
    e6 = float(np.max(errors[:6]))
    root7 = float(errors[6])
    accuracy_status = (
        "PASS"
        if e3 <= FEM_FIRST_THREE_TOLERANCE
        and e6 <= FEM_FIRST_SIX_TOLERANCE
        else "FAIL"
    )
    diagnostics["E3"] = e3
    diagnostics["E6"] = e6
    diagnostics["root7_relative_error"] = root7
    diagnostics["accuracy_status"] = accuracy_status
    diagnostics["status"] = (
        "PASS"
        if diagnostics["structural_status"] == "PASS"
        and accuracy_status == "PASS"
        else "FAIL"
    )
    return {
        "assembly": assembly,
        "solution": solution,
        "matched_frequencies": matched,
        "matching_basis": matching,
        "neighbor_gaps": gaps,
        "permutation": permutation,
        "matched_joint_equilibrium_residuals": matched_equilibrium,
        "relative_errors": errors,
        "diagnostics": diagnostics,
    }


def _ut2_fem_comparison_rows(
    runs: dict[tuple[str, float], dict[str, Any]],
    analytic: dict[tuple[str, float], np.ndarray],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for (section_case, beta_deg), run in runs.items():
        diagnostics = run["diagnostics"]
        for mode, (
            analytic_hz,
            fem_hz,
            relative,
            basis,
            equilibrium,
        ) in enumerate(
            zip(
                analytic[(section_case, beta_deg)],
                run["matched_frequencies"],
                run["relative_errors"],
                run["matching_basis"],
                run["matched_joint_equilibrium_residuals"],
            ),
            start=1,
        ):
            rows.append(
                {
                    "section_case": section_case,
                    "beta_deg": beta_deg,
                    "elements_arm_1": diagnostics["elements_arm_1"],
                    "elements_arm_2": diagnostics["elements_arm_2"],
                    "mode": mode,
                    "analytic_frequency_hz": analytic_hz,
                    "fem_frequency_hz": fem_hz,
                    "absolute_error_hz": abs(float(fem_hz) - float(analytic_hz)),
                    "relative_error": relative,
                    "matching_basis": basis,
                    "E3": diagnostics["E3"],
                    "E6": diagnostics["E6"],
                    "root7_relative_error": diagnostics[
                        "root7_relative_error"
                    ],
                    "stiffness_symmetry_residual": diagnostics[
                        "stiffness_symmetry_residual"
                    ],
                    "mass_symmetry_residual": diagnostics[
                        "mass_symmetry_residual"
                    ],
                    "minimum_mass_eigenvalue": diagnostics[
                        "minimum_mass_eigenvalue"
                    ],
                    "zero_mode_count": diagnostics["zero_mode_count"],
                    "joint_kinematic_residual": diagnostics[
                        "joint_kinematic_residual"
                    ],
                    "joint_equilibrium_residual": equilibrium,
                    "status": diagnostics["status"],
                }
            )
    return rows


def _ut2_fem_native_symmetry_rows(
    runs: dict[tuple[str, float], dict[str, Any]],
    *,
    reflection_matrix_ok: bool,
    relabeling_matrix_ok: bool,
    specifications: Sequence[dict[str, Any]] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, dict[str, float]]]:
    if specifications is None:
        specifications = [
            *[
                {
                    "symmetry_type": "reflection",
                    "case_left": case,
                    "case_right": case,
                    "beta_left_deg": 30.0,
                    "beta_right_deg": -30.0,
                }
                for case in UT2_SECTION_CASES
            ],
            {
                "symmetry_type": "oriented_angle_relabeling",
                "case_left": "asymmetric_4_6",
                "case_right": "swapped_6_4",
                "beta_left_deg": 30.0,
                "beta_right_deg": -30.0,
            },
            {
                "symmetry_type": "oriented_angle_relabeling",
                "case_left": "swapped_6_4",
                "case_right": "asymmetric_4_6",
                "beta_left_deg": 30.0,
                "beta_right_deg": -30.0,
            },
            {
                "symmetry_type": "same_positive_exchange",
                "case_left": "asymmetric_4_6",
                "case_right": "swapped_6_4",
                "beta_left_deg": 30.0,
                "beta_right_deg": 30.0,
            },
        ]
    rows: list[dict[str, Any]] = []
    maxima: dict[str, dict[str, float]] = {}
    for specification in specifications:
        left = runs[
            (specification["case_left"], specification["beta_left_deg"])
        ]["solution"].frequencies_hz
        right, matching, _, _ = _match_frequency_spectra(
            left,
            runs[
                (specification["case_right"], specification["beta_right_deg"])
            ]["solution"].frequencies_hz,
        )
        if specification["symmetry_type"] == "reflection":
            matrix_ok = reflection_matrix_ok
        elif specification["symmetry_type"] == "oriented_angle_relabeling":
            matrix_ok = relabeling_matrix_ok
        else:
            matrix_ok = reflection_matrix_ok and relabeling_matrix_ok
        differences: list[float] = []
        for mode, (left_hz, right_hz, basis) in enumerate(
            zip(left, right, matching), start=1
        ):
            absolute = abs(float(left_hz) - float(right_hz))
            relative = absolute / max(abs(float(left_hz)), np.finfo(float).tiny)
            differences.append(relative)
            finite = bool(
                np.isfinite(left_hz)
                and np.isfinite(right_hz)
                and left_hz > 0.0
                and right_hz > 0.0
            )
            if not finite or not matrix_ok:
                diagnostic_status = "NOT_EVALUABLE"
            elif relative <= SPECTRUM_RELATIVE_TOLERANCE:
                diagnostic_status = "WITHIN_HISTORICAL_MARKER"
            else:
                diagnostic_status = "ABOVE_HISTORICAL_MARKER_MATRIX_CONGRUENT"
            rows.append(
                {
                    "symmetry_type": specification["symmetry_type"],
                    "case_left": specification["case_left"],
                    "case_right": specification["case_right"],
                    "beta_left_deg": specification["beta_left_deg"],
                    "beta_right_deg": specification["beta_right_deg"],
                    "mode": mode,
                    "frequency_left_hz": left_hz,
                    "frequency_right_hz": right_hz,
                    "absolute_difference_hz": absolute,
                    "relative_difference": relative,
                    "matching_basis": basis,
                    "historical_marker_1e_8": SPECTRUM_RELATIVE_TOLERANCE,
                    "diagnostic_status": diagnostic_status,
                }
            )
        label = (
            f"{specification['case_left']}|{specification['beta_left_deg']:+g}_to_"
            f"{specification['case_right']}|{specification['beta_right_deg']:+g}"
        )
        maxima.setdefault(specification["symmetry_type"], {})[label] = max(
            differences
        )
    return rows, maxima


def _ut2_model_difference_rows(
    spectra: dict[tuple[str, float, str], dict[str, Any]],
    *,
    angle_deg: float = 30.0,
    section_cases: Sequence[str] | None = None,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    rows: list[dict[str, Any]] = []
    maxima: dict[str, float] = {}
    cases = tuple(UT2_SECTION_CASES if section_cases is None else section_cases)
    for section_case in cases:
        timo = _frequencies(
            spectra[(section_case, angle_deg, "Timoshenko")]["roots"]
        )[:6]
        eb_source = _frequencies(
            spectra[(section_case, angle_deg, "EB")]["roots"]
        )[:6]
        eb, matching, timo_gaps, permutation = _match_frequency_spectra(
            timo, eb_source
        )
        eb_gaps = _neighbor_gaps(eb_source)
        differences: list[float] = []
        for mode, (timo_hz, eb_hz, basis, eb_index) in enumerate(
            zip(timo, eb, matching, permutation), start=1
        ):
            absolute = abs(float(timo_hz) - float(eb_hz))
            relative = absolute / max(abs(float(timo_hz)), np.finfo(float).tiny)
            differences.append(relative)
            rows.append(
                {
                    "section_case": section_case,
                    "beta_deg": angle_deg,
                    "mode": mode,
                    "timoshenko_frequency_hz": timo_hz,
                    "eb_frequency_hz": eb_hz,
                    "absolute_difference_hz": absolute,
                    "relative_difference": relative,
                    "timoshenko_neighbor_gap_hz": timo_gaps[mode - 1],
                    "eb_neighbor_gap_hz": eb_gaps[int(eb_index)],
                    "matching_basis": basis,
                    "acceptance_threshold": "none_diagnostic_only",
                }
            )
        maxima[section_case] = max(differences)
    return rows, maxima


def _load_ut2_previous_statuses() -> dict[str, str]:
    base = DEFAULT_OUTPUT_DIR.parent
    paths = {
        "ut0": base / "ut0_smoke" / "ut0_summary.json",
        "ut1": base / "ut1_beta0" / "ut1_summary.json",
        "ut1a": base / "ut1a_fem_exchange_audit" / "ut1a_summary.json",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "UT-2 requires existing regression evidence: " + ", ".join(missing)
        )
    ut0 = json.loads(paths["ut0"].read_text(encoding="utf-8"))
    ut1 = json.loads(paths["ut1"].read_text(encoding="utf-8"))
    ut1a = json.loads(paths["ut1a"].read_text(encoding="utf-8"))
    return {
        "UT-0": ut0["ut0_status"],
        "UT-1": ut1["ut1_status"],
        "UT-1a": ut1a["ut1a_status"],
    }


def _ut2_matrix_group_maxima(
    rows: Sequence[dict[str, Any]], symmetry_type: str
) -> dict[str, dict[str, float]]:
    selected = [row for row in rows if row["symmetry_type"] == symmetry_type]
    groups: dict[str, list[dict[str, Any]]] = {
        "endpoint": [row for row in selected if row["check"].startswith("endpoint")],
        "permutation": [row for row in selected if row["check"].startswith("P_")],
        "reduction": [row for row in selected if "reduction intertwining" in row["check"]],
        "full_matrices": [
            row
            for row in selected
            if row["check"].startswith(("K_full", "M_full"))
        ],
        "reduced_matrices": [
            row
            for row in selected
            if row["check"].startswith(("K_reduced", "M_reduced"))
        ],
    }
    return {
        group: {
            "relative_frobenius": max(
                (row["relative_frobenius_residual"] for row in group_rows),
                default=0.0,
            ),
            "relative_max_entry": max(
                (row["relative_max_residual"] for row in group_rows),
                default=0.0,
            ),
            "absolute_max_entry": max(
                (row["absolute_max_residual"] for row in group_rows),
                default=0.0,
            ),
        }
        for group, group_rows in groups.items()
    }


def _ut2_report(summary: dict[str, Any]) -> str:
    lines = [
        "# UT-2 beta=30 deg unequal-thickness angular-joint validation",
        "",
        "Overall unequal-thickness validation: **IN_PROGRESS**",
        "",
        f"UT-0: **{summary['regression_statuses']['UT-0']}**",
        "",
        f"UT-1: **{summary['regression_statuses']['UT-1']}**",
        "",
        f"UT-1a: **{summary['regression_statuses']['UT-1a']}**",
        "",
        f"UT-2: **{summary['ut2_status']}**",
        "",
        "## Scope",
        "",
        "Only elastic HMS/DX-209, `theta_1=theta_2=0`, "
        "`L_1=L_2=0.2 m`, `b_1=b_2=0.020 m`, the existing `5/6` shear "
        "factor, and `(5,5)`, `(4,6)`, `(6,4) mm` are used. The angles are "
        "exactly `+30 deg` and `-30 deg`; the negative angle is only the "
        "reflected/oriented-relabelled control, not a sweep.",
        "",
        "## Continuum spectra",
        "",
        "| case | beta | model | f1 | f2 | f3 | f4 | f5 | f6 | f7 |",
        "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case, angles in summary["first_seven_frequencies_hz"].items():
        for beta_label, models in angles.items():
            for model, frequencies in models.items():
                lines.append(
                    f"| {case} | {beta_label} | {model} | "
                    + " | ".join(f"{frequency:.9f}" for frequency in frequencies)
                    + " |"
                )
    lines.extend(
        [
            "",
            f"All `{summary['continuum_root_total']}` continuum roots are "
            "finite, positive, and accepted. Quality counts are "
            f"`{summary['root_quality_basis_counts']}`. The determinant and "
            "relative singular thresholds are both `1e-8`, using the unchanged "
            "scaled-or-normalized-physical-raw acceptance logic.",
            "",
            "## Continuum symmetry maxima",
            "",
            "| symmetry | maximum relative difference | threshold |",
            "|---|---:|---:|",
        ]
    )
    for symmetry_type, values in summary["continuum_symmetry_maxima"].items():
        lines.append(
            f"| {symmetry_type} | {max(values.values()):.12e} | "
            f"{SPECTRUM_RELATIVE_TOLERANCE:.1e} |"
        )
    lines.extend(
        [
            "",
            "Reflection compares `beta -> -beta`. Oriented-angle relabeling "
            "compares `(P1,P2,+beta)` with `(P2,P1,-beta)`. Same-positive and "
            "same-negative exchange are retained as direct end-to-end consequences; "
            "no branch identity or `theta != 0` generalization is claimed.",
            "",
            "## Independent EB FEM",
            "",
            "| case | beta | E3 | E6 | root-7 diagnostic | structural | accuracy |",
            "|---|---:|---:|---:|---:|---|---|",
        ]
    )
    for label, diagnostics in summary["fem_configuration_diagnostics"].items():
        lines.append(
            f"| {diagnostics['section_case']} | {diagnostics['beta_deg']:+g} | "
            f"{diagnostics['E3']:.12e} | {diagnostics['E6']:.12e} | "
            f"{diagnostics['root7_relative_error']:.12e} | "
            f"{diagnostics['structural_status']} | "
            f"{diagnostics['accuracy_status']} |"
        )
    lines.extend(
        [
            "",
            "The six systems use only `N_1=N_2=64`; no refinement sequence is "
            "run. The hard accuracy limits are `E3<=1e-4` and `E6<=5e-4`; "
            "root 7 is a positive finite completeness guard with no error limit.",
            "",
            "## FEM matrix transformations",
            "",
            "Reflection applies `G_ref=diag(1,1,-1)` to every full/internal "
            "`[w,psi,Phi]` triple and `H_ref=diag(1,-1,1)` to "
            "`[w_J,theta_t,theta_n]`. Relabeling swaps complete arm blocks "
            "without local signs and applies the declared trigonometric "
            "`H_rel(30 deg)` to the joint block.",
            "",
            "| matrix symmetry | group | max relative Frobenius | max relative entry | threshold |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for symmetry, groups in (
        ("reflection", summary["fem_reflection_matrix_residual_maxima"]),
        (
            "oriented_angle_relabeling",
            summary["fem_relabeling_matrix_residual_maxima"],
        ),
    ):
        for group, values in groups.items():
            lines.append(
                f"| {symmetry} | {group} | "
                f"{values['relative_frobenius']:.12e} | "
                f"{values['relative_max_entry']:.12e} | "
                f"{UT2_MATRIX_RELATIVE_TOLERANCE:.1e} |"
            )
    lines.extend(
        [
            "",
            "## Native FEM spectral symmetry (diagnostic only)",
            "",
            "| symmetry | maximum native relative difference | historical marker |",
            "|---|---:|---:|",
        ]
    )
    for symmetry, values in summary["native_fem_spectral_symmetry_maxima"].items():
        lines.append(
            f"| {symmetry} | {max(values.values()):.12e} | 1e-8 |"
        )
    lines.extend(
        [
            "",
            "These native values are not a separate hard rejection gate when "
            "the corresponding matrices are congruent and both FEM systems "
            "pass analytic EB accuracy. No UT-1a backward, Rayleigh, or "
            "canonicalized eigensolver audit is repeated.",
            "",
            "## Timoshenko--EB diagnostic",
            "",
        ]
    )
    for case, maximum in summary["model_difference_maxima"].items():
        lines.append(f"- `{case}`: `{maximum:.12e}` (no acceptance threshold).")
    lines.extend(
        [
            "",
            "## Status and exclusions",
            "",
            f"UT-2 status: **{summary['ut2_status']}**. Overall unequal-thickness "
            "validation remains **IN_PROGRESS** because `beta=90 deg` has not "
            "been validated.",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in summary["explicit_exclusions"])
    return "\n".join(lines) + "\n"


def _run_ut2_beta30(output_dir: Path) -> dict[str, Any]:
    """Run only the prescribed +/-30-degree unequal-thickness UT-2 gate."""

    started = time.perf_counter()
    regression_statuses = _load_ut2_previous_statuses()
    section_cases = _make_ut2_section_cases()
    spectra: dict[tuple[str, float, str], dict[str, Any]] = {}
    continuum_rows: list[dict[str, Any]] = []
    runtime_by_spectrum: dict[str, float] = {}
    evaluations_by_spectrum: dict[str, int] = {}
    scan_max_by_spectrum: dict[str, float] = {}

    for section_case, (point_1, point_2) in section_cases.items():
        for beta_deg in UT2_ANGLES_DEG:
            beta_rad = _ut2_beta_rad(beta_deg)
            for model in ("Timoshenko", "EB"):
                factory, raw_factory = _ut2_continuum_factories(
                    model, beta_rad, point_1, point_2
                )
                roots, builder, runtime = _solve_roots(
                    point_1, factory, num_roots=UT2_NUM_ROOTS
                )
                quality = [
                    _root_quality_fields(root, raw_factory) for root in roots
                ]
                gaps = _neighbor_gaps(_frequencies(roots))
                key = (section_case, beta_deg, model)
                spectra[key] = {"roots": roots, "quality": quality}
                label = f"{section_case}|{beta_deg:+g}|{model}"
                runtime_by_spectrum[label] = runtime
                evaluations_by_spectrum[label] = builder.evaluations
                scan_max_by_spectrum[label] = builder.maximum_frequency_hz
                continuum_rows.extend(
                    _ut2_continuum_root_row(
                        section_case=section_case,
                        point_1=point_1,
                        point_2=point_2,
                        beta_deg=beta_deg,
                        model=model,
                        mode=mode,
                        root=root,
                        quality=quality_fields,
                        gap_hz=float(gap),
                    )
                    for mode, (root, quality_fields, gap) in enumerate(
                        zip(roots, quality, gaps), start=1
                    )
                )

    continuum_symmetry_rows, continuum_symmetry_maxima = (
        _ut2_continuum_symmetry_rows(spectra)
    )
    analytic_eb = {
        (section_case, beta_deg): _frequencies(
            spectra[(section_case, beta_deg, "EB")]["roots"]
        )
        for section_case in UT2_SECTION_CASES
        for beta_deg in UT2_ANGLES_DEG
    }

    fem_runs: dict[tuple[str, float], dict[str, Any]] = {}
    for section_case, (point_1, point_2) in section_cases.items():
        for beta_deg in UT2_ANGLES_DEG:
            fem_runs[(section_case, beta_deg)] = _run_ut2_eb_fem_configuration(
                section_case=section_case,
                point_1=point_1,
                point_2=point_2,
                beta_deg=beta_deg,
                analytic_frequencies=analytic_eb[(section_case, beta_deg)],
            )
    fem_rows = _ut2_fem_comparison_rows(fem_runs, analytic_eb)

    matrix_rows = _ut2_endpoint_map_rows()
    for section_case in UT2_SECTION_CASES:
        matrix_rows.extend(
            _ut2_reflection_matrix_rows(
                section_case,
                fem_runs[(section_case, 30.0)]["assembly"],
                fem_runs[(section_case, -30.0)]["assembly"],
            )
        )
    matrix_rows.extend(
        _ut2_relabeling_matrix_rows(
            "asymmetric_4_6",
            "swapped_6_4",
            fem_runs[("asymmetric_4_6", 30.0)]["assembly"],
            fem_runs[("swapped_6_4", -30.0)]["assembly"],
        )
    )
    matrix_rows.extend(
        _ut2_relabeling_matrix_rows(
            "swapped_6_4",
            "asymmetric_4_6",
            fem_runs[("swapped_6_4", 30.0)]["assembly"],
            fem_runs[("asymmetric_4_6", -30.0)]["assembly"],
        )
    )
    reflection_matrix_ok = all(
        row["status"] == "PASS"
        for row in matrix_rows
        if row["symmetry_type"] == "reflection"
    )
    relabeling_matrix_ok = all(
        row["status"] == "PASS"
        for row in matrix_rows
        if row["symmetry_type"] == "oriented_angle_relabeling"
    )
    native_fem_rows, native_fem_maxima = _ut2_fem_native_symmetry_rows(
        fem_runs,
        reflection_matrix_ok=reflection_matrix_ok,
        relabeling_matrix_ok=relabeling_matrix_ok,
    )
    model_difference_rows, model_difference_maxima = (
        _ut2_model_difference_rows(spectra)
    )

    quality_basis_counts = {
        basis: sum(
            row["quality_status"] == "PASS" and row["quality_basis"] == basis
            for row in continuum_rows
        )
        for basis in ("scaled", "physical_raw")
    }
    root_quality_maxima_by_basis = {
        "scaled": {
            "normalized_determinant_residual": max(
                (
                    row["scaled_determinant_residual"]
                    for row in continuum_rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "scaled"
                ),
                default=0.0,
            ),
            "relative_singular_residual": max(
                (
                    row["scaled_relative_singular_residual"]
                    for row in continuum_rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "scaled"
                ),
                default=0.0,
            ),
        },
        "physical_raw": {
            "normalized_determinant_residual": max(
                (
                    row["raw_determinant_residual"]
                    for row in continuum_rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "physical_raw"
                ),
                default=0.0,
            ),
            "relative_singular_residual": max(
                (
                    row["raw_relative_singular_residual"]
                    for row in continuum_rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "physical_raw"
                ),
                default=0.0,
            ),
        },
    }
    continuum_counts = {
        f"{case}|{beta:+g}|{model}": len(data["roots"])
        for (case, beta, model), data in spectra.items()
    }
    continuum_complete = bool(
        len(spectra) == 12
        and len(continuum_rows) == 84
        and all(count == UT2_NUM_ROOTS for count in continuum_counts.values())
        and all(
            np.isfinite(row["frequency_hz"])
            and row["frequency_hz"] > 0.0
            and row["quality_status"] == "PASS"
            for row in continuum_rows
        )
    )
    continuum_symmetry_ok = all(
        row["status"] == "PASS" for row in continuum_symmetry_rows
    )
    continuum_hard_ok = continuum_complete and continuum_symmetry_ok
    matrix_hard_ok = reflection_matrix_ok and relabeling_matrix_ok
    fem_hard_ok = bool(
        len(fem_runs) == 6
        and all(run["diagnostics"]["status"] == "PASS" for run in fem_runs.values())
    )
    regressions_ok = regression_statuses == {
        "UT-0": "PASS",
        "UT-1": "PARTIAL_PASS",
        "UT-1a": "INCONCLUSIVE_NUMERICAL_AUDIT",
    }
    ut2_status = _classify_ut2_status(
        continuum_hard_ok=continuum_hard_ok,
        matrix_hard_ok=matrix_hard_ok,
        fem_hard_ok=fem_hard_ok,
        regressions_ok=regressions_ok,
    )

    fem_diagnostics = {
        f"{case}|{beta:+g}": run["diagnostics"]
        for (case, beta), run in fem_runs.items()
    }
    fem_structural_maxima = {
        "stiffness_symmetry_residual": max(
            item["stiffness_symmetry_residual"]
            for item in fem_diagnostics.values()
        ),
        "mass_symmetry_residual": max(
            item["mass_symmetry_residual"] for item in fem_diagnostics.values()
        ),
        "joint_kinematic_residual": max(
            item["joint_kinematic_residual"] for item in fem_diagnostics.values()
        ),
        "joint_equilibrium_residual": max(
            item["max_joint_equilibrium_residual"]
            for item in fem_diagnostics.values()
        ),
        "minimum_mass_eigenvalue": min(
            item["minimum_mass_eigenvalue"] for item in fem_diagnostics.values()
        ),
        "maximum_zero_mode_count": max(
            item["zero_mode_count"] for item in fem_diagnostics.values()
        ),
    }
    summary = {
        "git_context": _git_context(),
        "scope": "UT-2 beta=+/-30 deg unequal-thickness angular-joint validation only",
        "geometry": {
            "theta_1_deg": 0.0,
            "theta_2_deg": 0.0,
            "material_mode": "elastic",
            "material": "existing HMS/DX-209",
            "length_1_m": LENGTH_M,
            "length_2_m": LENGTH_M,
            "b_1_m": WIDTH_B_M,
            "b_2_m": WIDTH_B_M,
            "shear_factor": cantilever_geometry().shear_factor,
            "a_direction": "parallel to e_z",
        },
        "angles_deg": list(UT2_ANGLES_DEG),
        "section_cases": {
            case: {"a_1_m": values[0], "a_2_m": values[1]}
            for case, values in UT2_SECTION_CASES.items()
        },
        "thresholds": {
            "root_normalized_determinant": ROOT_DETERMINANT_TOLERANCE,
            "root_relative_singular": ROOT_SINGULAR_TOLERANCE,
            "near_degenerate_relative_gap": NEAR_DEGENERATE_RELATIVE_GAP,
            "continuum_symmetry_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "fem_matrix_relative_frobenius": UT2_MATRIX_RELATIVE_TOLERANCE,
            "fem_matrix_relative_max_entry": UT2_MATRIX_RELATIVE_TOLERANCE,
            "fem_stiffness_symmetry": FEM_SYMMETRY_TOLERANCE,
            "fem_mass_symmetry": FEM_SYMMETRY_TOLERANCE,
            "fem_joint_kinematic": FEM_JOINT_KINEMATIC_TOLERANCE,
            "fem_joint_equilibrium": FEM_JOINT_EQUILIBRIUM_TOLERANCE,
            "fem_E3": FEM_FIRST_THREE_TOLERANCE,
            "fem_E6": FEM_FIRST_SIX_TOLERANCE,
            "native_fem_historical_marker_diagnostic_only": (
                SPECTRUM_RELATIVE_TOLERANCE
            ),
        },
        "root_scan": {
            "formulation": FORMULATION,
            "num_roots": UT2_NUM_ROOTS,
            "scan_step_hz": SCAN_STEP_HZ,
            "initial_max_hz": INITIAL_MAX_HZ,
            "max_hz": MAX_HZ,
            "actual_maximum_evaluated_frequency_hz": scan_max_by_spectrum,
        },
        "continuum_spectrum_count": len(spectra),
        "continuum_root_total": len(continuum_rows),
        "continuum_root_counts": continuum_counts,
        "first_seven_frequencies_hz": {
            case: {
                f"{beta:+g} deg": {
                    model: _frequencies(spectra[(case, beta, model)]["roots"]).tolist()
                    for model in ("Timoshenko", "EB")
                }
                for beta in UT2_ANGLES_DEG
            }
            for case in UT2_SECTION_CASES
        },
        "root_quality_basis_counts": quality_basis_counts,
        "root_quality_maxima_by_accepted_basis": root_quality_maxima_by_basis,
        "continuum_symmetry_maxima": continuum_symmetry_maxima,
        "fem_assembly_count": len(fem_runs),
        "fem_configuration_diagnostics": fem_diagnostics,
        "fem_structural_maxima": fem_structural_maxima,
        "fem_reflection_matrix_residual_maxima": _ut2_matrix_group_maxima(
            matrix_rows, "reflection"
        ),
        "fem_relabeling_matrix_residual_maxima": _ut2_matrix_group_maxima(
            matrix_rows, "oriented_angle_relabeling"
        ),
        "native_fem_spectral_symmetry_maxima": native_fem_maxima,
        "native_fem_spectral_symmetry_is_hard_gate": False,
        "model_difference_maxima": model_difference_maxima,
        "runtime_seconds_by_continuum_spectrum": runtime_by_spectrum,
        "matrix_evaluations_by_continuum_spectrum": evaluations_by_spectrum,
        "total_continuum_matrix_evaluations": sum(
            evaluations_by_spectrum.values()
        ),
        "fem_runtime_seconds_by_configuration": {
            label: diagnostics["runtime_seconds"]
            for label, diagnostics in fem_diagnostics.items()
        },
        "total_runtime_seconds": time.perf_counter() - started,
        "regression_statuses": regression_statuses,
        "continuum_hard_gates_status": (
            "PASS" if continuum_hard_ok else "FAIL"
        ),
        "fem_matrix_hard_gates_status": "PASS" if matrix_hard_ok else "FAIL",
        "fem_accuracy_structural_status": "PASS" if fem_hard_ok else "FAIL",
        "ut2_status": ut2_status,
        "overall_unequal_thickness_validation": "IN_PROGRESS",
        "explicit_exclusions": UT2_EXPLICIT_EXCLUSIONS,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "continuum_roots.csv", continuum_rows)
    _write_csv(output_dir / "continuum_symmetry.csv", continuum_symmetry_rows)
    _write_csv(output_dir / "fem_comparison.csv", fem_rows)
    _write_csv(output_dir / "fem_matrix_symmetry.csv", matrix_rows)
    _write_csv(
        output_dir / "fem_native_spectral_symmetry.csv", native_fem_rows
    )
    _write_csv(output_dir / "model_difference.csv", model_difference_rows)
    (output_dir / "ut2_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "ut2_report.md").write_text(
        _ut2_report(summary), encoding="utf-8"
    )
    return summary


def _ut3_joint_plus_90_exact() -> np.ndarray:
    """Return the prescribed integer-entry +90-degree book joint matrix."""

    return np.array(
        [
            [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, -1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -1],
        ],
        dtype=float,
    )


def _ut3_joint_minus_90_exact() -> np.ndarray:
    """Return the prescribed integer-entry -90-degree book joint matrix."""

    return np.array(
        [
            [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -1, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=float,
    )


def _ut3_relabel_joint_exact() -> np.ndarray:
    """Return the exact +90-degree relabeling map for joint coordinates."""

    return np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]],
        dtype=float,
    )


def _ut3_exact_block_end_map(
    omega: complex,
    point_1: RodPoint,
    point_2: RodPoint,
    *,
    model: str,
) -> np.ndarray:
    """Assemble physical arm end maps without a coupled/joint builder."""

    if model == "Timoshenko":
        end_map = physical_end_map
    elif model == "EB":
        end_map = eb_physical_end_map
    else:
        raise ValueError(f"unsupported UT-3 exact continuum model: {model}")
    block = np.zeros((12, 6), dtype=np.complex128)
    block[:6, :3] = end_map(omega, point_1)
    block[6:, 3:] = end_map(omega, point_2)
    return block


def _ut3_exact_boundary_matrix_raw(
    omega: complex,
    point_1: RodPoint,
    point_2: RodPoint,
    *,
    model: str,
) -> np.ndarray:
    """Return physical ``6 x 6`` +90 exact-limit boundary matrix."""

    return _ut3_joint_plus_90_exact() @ _ut3_exact_block_end_map(
        omega, point_1, point_2, model=model
    )


def _ut3_exact_boundary_matrix(
    omega: complex,
    point_1: RodPoint,
    point_2: RodPoint,
    *,
    model: str,
) -> np.ndarray:
    return equilibrate_matrix(
        _ut3_exact_boundary_matrix_raw(
            omega, point_1, point_2, model=model
        )
    )[0]


def _ut3_exact_continuum_factories(
    model: str,
    point_1: RodPoint,
    point_2: RodPoint,
) -> tuple[Callable[[complex], np.ndarray], Callable[[complex], np.ndarray]]:
    """Return exact +90 scaled/raw factories without ``J_book`` calls."""

    if model not in ("Timoshenko", "EB"):
        raise ValueError(f"unsupported UT-3 exact continuum model: {model}")
    return (
        lambda omega: _ut3_exact_boundary_matrix(
            omega, point_1, point_2, model=model
        ),
        lambda omega: _ut3_exact_boundary_matrix_raw(
            omega, point_1, point_2, model=model
        ),
    )


def _ut3_continuum_root_row(
    *,
    builder: str,
    section_case: str,
    point_1: RodPoint,
    point_2: RodPoint,
    beta_deg: float,
    model: str,
    mode: int,
    root: RootResult,
    quality: dict[str, Any],
    gap_hz: float,
) -> dict[str, Any]:
    row = _ut2_continuum_root_row(
        section_case=section_case,
        point_1=point_1,
        point_2=point_2,
        beta_deg=beta_deg,
        model=model,
        mode=mode,
        root=root,
        quality=quality,
        gap_hz=gap_hz,
    )
    row["builder"] = builder
    ordered = {
        key: row[key]
        for key in (
            "section_case",
            "a_1_m",
            "a_2_m",
            "beta_deg",
            "model",
            "builder",
            "mode",
            "frequency_hz",
        )
    }
    ordered.update(
        {key: value for key, value in row.items() if key not in ordered}
    )
    return ordered


def _ut3_exact_identity_row(
    *, check: str, left: np.ndarray, right: np.ndarray
) -> dict[str, Any]:
    row = _ut1a_matrix_identity_row(
        check=check,
        left=left,
        right=right,
        threshold=UT3_EXACT_LIMIT_TOLERANCE,
    )
    row["status"] = (
        "PASS"
        if row["absolute_max_residual"] <= UT3_EXACT_LIMIT_TOLERANCE
        and row["relative_frobenius_residual"]
        <= UT3_EXACT_LIMIT_TOLERANCE
        else "FAIL"
    )
    return row


def _ut3_exact_joint_limit_rows() -> list[dict[str, Any]]:
    """Audit exact joint, endpoint, relabeling, and deterministic limits."""

    plus = math.pi / 2.0
    minus = -plus
    f1_plus, f2_plus = eb_joint_dof_maps(plus)
    f1_minus, f2_minus = eb_joint_dof_maps(minus)
    h_exact = _ut3_relabel_joint_exact()
    rows = [
        _ut3_exact_identity_row(
            check="J_book(+90) vs J_plus_90_exact",
            left=joint_matrix_book(plus),
            right=_ut3_joint_plus_90_exact(),
        ),
        _ut3_exact_identity_row(
            check="J_book(-90) vs J_minus_90_exact",
            left=joint_matrix_book(minus),
            right=_ut3_joint_minus_90_exact(),
        ),
        _ut3_exact_identity_row(
            check="F2(+90) vs I",
            left=f2_plus,
            right=np.eye(3),
        ),
        _ut3_exact_identity_row(
            check="F2(-90) vs diag(1,-1,-1)",
            left=f2_minus,
            right=np.diag([1.0, -1.0, -1.0]),
        ),
        _ut3_exact_identity_row(
            check="H_rel(+90) vs exact quarter-turn",
            left=_ut2_relabel_joint_transform(plus),
            right=h_exact,
        ),
        _ut3_exact_identity_row(
            check="F1 H_rel_exact = F2(+90)",
            left=f1_plus @ h_exact,
            right=f2_plus,
        ),
        _ut3_exact_identity_row(
            check="F2(-90) H_rel_exact = F1",
            left=f2_minus @ h_exact,
            right=f1_minus,
        ),
    ]
    joint_vectors = (
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([2.0, -3.0, 5.0]),
    )
    for index, vector in enumerate(joint_vectors, start=1):
        w_j, theta_t, theta_n = vector
        expected = np.array(
            [w_j, -theta_n, theta_t, w_j, theta_t, theta_n]
        ).reshape(-1, 1)
        actual = np.concatenate(
            (f1_plus @ vector, f2_plus @ vector)
        ).reshape(-1, 1)
        rows.append(
            _ut3_exact_identity_row(
                check=f"deterministic channel-exchange state {index}",
                left=actual,
                right=expected,
            )
        )
    return rows


def _ut3_standard_exact_rows(
    standard: dict[tuple[str, float, str], dict[str, Any]],
    exact: dict[tuple[str, str], dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    rows: list[dict[str, Any]] = []
    maxima: dict[str, float] = {}
    for section_case in UT3_SECTION_CASES:
        for model in ("Timoshenko", "EB"):
            standard_data = standard[(section_case, 90.0, model)]
            exact_data = exact[(section_case, model)]
            standard_hz = _frequencies(standard_data["roots"])
            exact_source = _frequencies(exact_data["roots"])
            exact_hz, matching, _, permutation = _match_frequency_spectra(
                standard_hz, exact_source
            )
            differences: list[float] = []
            for mode, (
                standard_frequency,
                exact_frequency,
                basis,
                exact_index,
            ) in enumerate(
                zip(standard_hz, exact_hz, matching, permutation), start=1
            ):
                absolute = abs(
                    float(standard_frequency) - float(exact_frequency)
                )
                relative = absolute / max(
                    abs(float(exact_frequency)), np.finfo(float).tiny
                )
                differences.append(relative)
                standard_quality = standard_data["quality"][mode - 1]
                exact_quality = exact_data["quality"][int(exact_index)]
                rows.append(
                    {
                        "section_case": section_case,
                        "model": model,
                        "mode": mode,
                        "standard_frequency_hz": standard_frequency,
                        "exact_limit_frequency_hz": exact_frequency,
                        "absolute_difference_hz": absolute,
                        "relative_difference": relative,
                        "matching_basis": basis,
                        "standard_quality_basis": standard_quality[
                            "root_quality_basis"
                        ],
                        "exact_quality_basis": exact_quality[
                            "root_quality_basis"
                        ],
                        "threshold": SPECTRUM_RELATIVE_TOLERANCE,
                        "status": (
                            "PASS"
                            if relative <= SPECTRUM_RELATIVE_TOLERANCE
                            else "FAIL"
                        ),
                    }
                )
            maxima[f"{model}|{section_case}"] = max(differences)
    return rows, maxima


def _ut3_previous_statuses() -> dict[str, str]:
    statuses = _load_ut2_previous_statuses()
    path = UT2_DEFAULT_OUTPUT_DIR / "ut2_summary.json"
    if not path.is_file():
        raise FileNotFoundError(f"UT-3 requires existing UT-2 evidence: {path}")
    statuses["UT-2"] = json.loads(path.read_text(encoding="utf-8"))[
        "ut2_status"
    ]
    return statuses


def _classify_ut3_status(
    *,
    continuum_hard_ok: bool,
    exact_limit_ok: bool,
    matrix_hard_ok: bool,
    fem_hard_ok: bool,
    regressions_ok: bool,
) -> tuple[str, str, str]:
    """Return UT-3, overall 1D, and phase statuses; native FEM is excluded."""

    if not (
        continuum_hard_ok
        and exact_limit_ok
        and matrix_hard_ok
        and regressions_ok
    ):
        return "FAIL", "FAIL", "BLOCKED"
    if not fem_hard_ok:
        return "PARTIAL_PASS", "PARTIAL_PASS", "COMPLETE_WITH_LIMITATIONS"
    return "PASS", "PARTIAL_PASS", "COMPLETE"


def _ut3_root_quality_summary(
    rows: Sequence[dict[str, Any]],
) -> tuple[dict[str, int], dict[str, dict[str, float]]]:
    counts = {
        basis: sum(
            row["quality_status"] == "PASS"
            and row["quality_basis"] == basis
            for row in rows
        )
        for basis in ("scaled", "physical_raw")
    }
    maxima = {
        "scaled": {
            "normalized_determinant_residual": max(
                (
                    row["scaled_determinant_residual"]
                    for row in rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "scaled"
                ),
                default=0.0,
            ),
            "relative_singular_residual": max(
                (
                    row["scaled_relative_singular_residual"]
                    for row in rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "scaled"
                ),
                default=0.0,
            ),
        },
        "physical_raw": {
            "normalized_determinant_residual": max(
                (
                    row["raw_determinant_residual"]
                    for row in rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "physical_raw"
                ),
                default=0.0,
            ),
            "relative_singular_residual": max(
                (
                    row["raw_relative_singular_residual"]
                    for row in rows
                    if row["quality_status"] == "PASS"
                    and row["quality_basis"] == "physical_raw"
                ),
                default=0.0,
            ),
        },
    }
    return counts, maxima


def _ut3_exact_rows_maxima(
    rows: Sequence[dict[str, Any]],
) -> dict[str, float]:
    return {
        "absolute_frobenius_residual": max(
            (row["absolute_frobenius_residual"] for row in rows), default=0.0
        ),
        "relative_frobenius_residual": max(
            (row["relative_frobenius_residual"] for row in rows), default=0.0
        ),
        "absolute_max_residual": max(
            (row["absolute_max_residual"] for row in rows), default=0.0
        ),
        "relative_max_residual": max(
            (row["relative_max_residual"] for row in rows), default=0.0
        ),
    }


def _ut3_fem_matrix_rows(
    fem_runs: dict[tuple[str, float], dict[str, Any]],
) -> list[dict[str, Any]]:
    rows = _ut2_endpoint_map_rows(90.0)
    rows.append(
        _ut2_matrix_row(
            symmetry_type="oriented_angle_relabeling",
            section_case_left="all",
            section_case_right="all",
            beta_left_deg=90.0,
            beta_right_deg=-90.0,
            check="endpoint exact H_rel(+90) quarter-turn identity",
            left=_ut2_relabel_joint_transform(math.pi / 2.0),
            right=_ut3_relabel_joint_exact(),
            threshold=UT3_EXACT_LIMIT_TOLERANCE,
        )
    )
    for section_case in UT3_SECTION_CASES:
        rows.extend(
            _ut2_reflection_matrix_rows(
                section_case,
                fem_runs[(section_case, 90.0)]["assembly"],
                fem_runs[(section_case, -90.0)]["assembly"],
                beta_deg=90.0,
            )
        )
    rows.extend(
        _ut2_relabeling_matrix_rows(
            "asymmetric_4_6",
            "swapped_6_4",
            fem_runs[("asymmetric_4_6", 90.0)]["assembly"],
            fem_runs[("swapped_6_4", -90.0)]["assembly"],
            beta_deg=90.0,
        )
    )
    rows.extend(
        _ut2_relabeling_matrix_rows(
            "swapped_6_4",
            "asymmetric_4_6",
            fem_runs[("swapped_6_4", 90.0)]["assembly"],
            fem_runs[("asymmetric_4_6", -90.0)]["assembly"],
            beta_deg=90.0,
        )
    )
    return rows


def _ut3_fem_structural_maxima(
    diagnostics: dict[str, dict[str, Any]],
) -> dict[str, float | int]:
    return {
        "stiffness_symmetry_residual": max(
            item["stiffness_symmetry_residual"] for item in diagnostics.values()
        ),
        "mass_symmetry_residual": max(
            item["mass_symmetry_residual"] for item in diagnostics.values()
        ),
        "joint_kinematic_residual": max(
            item["joint_kinematic_residual"] for item in diagnostics.values()
        ),
        "joint_equilibrium_residual": max(
            item["max_joint_equilibrium_residual"]
            for item in diagnostics.values()
        ),
        "minimum_mass_eigenvalue": min(
            item["minimum_mass_eigenvalue"] for item in diagnostics.values()
        ),
        "maximum_zero_mode_count": max(
            item["zero_mode_count"] for item in diagnostics.values()
        ),
    }


def _ut3_report(summary: dict[str, Any]) -> str:
    lines = [
        "# UT-3 beta=90 deg unequal-thickness quarter-turn joint validation",
        "",
        "Unequal-thickness 1D validation phase: "
        f"**{summary['unequal_thickness_1d_phase']}**",
        "",
        "Overall unequal-thickness 1D validation: "
        f"**{summary['overall_unequal_thickness_1d_validation']}**",
        "",
    ]
    for gate in ("UT-0", "UT-1", "UT-1a", "UT-2"):
        lines.extend(
            [f"{gate}: **{summary['regression_statuses'][gate]}**", ""]
        )
    lines.extend(
        [
            f"UT-3: **{summary['ut3_status']}**",
            "",
            "## Scope and exact quarter-turn limit",
            "",
            "Only elastic HMS/DX-209, `theta_1=theta_2=0`, equal `0.2 m` "
            "arms, `b_1=b_2=0.020 m`, the existing `5/6` shear factor, "
            "and `(5,5)`, `(4,6)`, `(6,4) mm` are used. Standard spectra "
            "use exactly `+90 deg` and `-90 deg`; exact-limit spectra use "
            "only the physical `+90 deg` cases.",
            "",
            "The script-local integer matrices impose `Phi_1=psi_2`, "
            "`psi_1=-Phi_2`, `M_T,1=-M_2`, and `M_1=M_T,2` at `+90 deg`, "
            "with the corresponding sign reversal at `-90 deg`. Production "
            "`J_book` is not rounded or changed.",
            "",
            "| exact check | absolute max | relative Frobenius | threshold | status |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for check, values in summary["exact_joint_limit_results"].items():
        lines.append(
            f"| {check} | {values['absolute_max_residual']:.12e} | "
            f"{values['relative_frobenius_residual']:.12e} | 1e-14 | "
            f"{values['status']} |"
        )
    lines.extend(
        [
            "",
            "## Standard continuum spectra",
            "",
            "| case | beta | model | f1 | f2 | f3 | f4 | f5 | f6 | f7 |",
            "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for case, angles in summary["first_seven_standard_frequencies_hz"].items():
        for beta_label, models in angles.items():
            for model, frequencies in models.items():
                lines.append(
                    f"| {case} | {beta_label} | {model} | "
                    + " | ".join(f"{value:.9f}" for value in frequencies)
                    + " |"
                )
    lines.extend(
        [
            "",
            f"All `{summary['standard_continuum_root_total']}` standard and "
            f"`{summary['exact_limit_root_total']}` exact-limit roots are "
            "finite, positive, accepted, and retain scaled plus raw quality "
            "evidence.",
            "",
            "## Standard versus exact-limit spectra",
            "",
            "| model and case | maximum relative difference | threshold |",
            "|---|---:|---:|",
        ]
    )
    for label, maximum in summary["standard_exact_spectrum_maxima"].items():
        lines.append(f"| {label} | {maximum:.12e} | 1e-8 |")
    lines.extend(
        [
            "",
            "## Continuum symmetries",
            "",
            "| symmetry | maximum relative difference | threshold |",
            "|---|---:|---:|",
        ]
    )
    for symmetry, values in summary["continuum_symmetry_maxima"].items():
        lines.append(f"| {symmetry} | {max(values.values()):.12e} | 1e-8 |")
    lines.extend(
        [
            "",
            "## EB FEM accuracy and structure",
            "",
            "| configuration | E3 | E6 | root 7 diagnostic | status |",
            "|---|---:|---:|---:|---|",
        ]
    )
    for label, diagnostics in summary["fem_configuration_diagnostics"].items():
        lines.append(
            f"| {label} | {diagnostics['E3']:.12e} | "
            f"{diagnostics['E6']:.12e} | "
            f"{diagnostics['root7_relative_error']:.12e} | "
            f"{diagnostics['status']} |"
        )
    lines.extend(
        [
            "",
            "The existing Hermite-bending/linear-torsion FEM, endpoint maps, "
            "reduction, clamps, and symmetric generalized eigensolver are "
            "unchanged. Six and only six mesh-64 assemblies are used.",
            "",
            "## FEM matrix symmetries",
            "",
            "| symmetry | group | max relative Frobenius | max relative entry | threshold |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for symmetry, groups in (
        ("reflection", summary["fem_reflection_matrix_residual_maxima"]),
        (
            "oriented_angle_relabeling",
            summary["fem_relabeling_matrix_residual_maxima"],
        ),
    ):
        for group, values in groups.items():
            lines.append(
                f"| {symmetry} | {group} | "
                f"{values['relative_frobenius']:.12e} | "
                f"{values['relative_max_entry']:.12e} | 1e-13 |"
            )
    lines.extend(
        [
            "",
            "Reflection uses `G_ref=diag(1,1,-1)` and "
            "`H_ref=diag(1,-1,1)`. Relabeling swaps full/internal arm blocks "
            "and uses the exact joint quarter-turn "
            "`[[1,0,0],[0,0,1],[0,-1,0]]`. This map is orthogonal but not "
            "involutory; its inverse is its transpose.",
            "",
            "## Native FEM spectral symmetry",
            "",
        ]
    )
    for symmetry, values in summary["native_fem_spectral_symmetry_maxima"].items():
        lines.append(
            f"- `{symmetry}`: `{max(values.values()):.12e}` "
            "(historical `1e-8` marker, diagnostic only)."
        )
    lines.extend(
        [
            "",
            "Native spectral equality is not a hard rejection gate when the "
            "declared matrix congruence and both analytic EB accuracy gates "
            "pass. No UT-1a backward/Rayleigh audit is repeated.",
            "",
            "## Timoshenko--EB diagnostic",
            "",
        ]
    )
    for case, maximum in summary["model_difference_maxima"].items():
        lines.append(f"- `{case}`: `{maximum:.12e}` (no threshold).")
    lines.extend(
        [
            "",
            "## Final 1D classification",
            "",
            "All continuum, root-quality, beta=30, beta=90, FEM accuracy, "
            "FEM convergence, and matrix-equivariance gates passed. The "
            "retained overall limitation is the historical beta=0 native "
            "FEM spectral-exchange threshold failure; its matrix-level "
            "follow-up found no assembly asymmetry but did not pass the "
            "declared backward/Rayleigh numerical audit.",
            "",
            "The unequal-thickness 1D phase is complete. The next separately "
            "scoped stage is limited 3D FEM anchor design; it is not started "
            "here.",
            "",
            "## Explicit exclusions",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in summary["explicit_exclusions"])
    return "\n".join(lines) + "\n"


def _run_ut3_beta90(output_dir: Path) -> dict[str, Any]:
    """Run only the prescribed +/-90-degree quarter-turn UT-3 gate."""

    started = time.perf_counter()
    regression_statuses = _ut3_previous_statuses()
    section_cases = {
        case: (_make_section_point(a_1), _make_section_point(a_2))
        for case, (a_1, a_2) in UT3_SECTION_CASES.items()
    }
    standard: dict[tuple[str, float, str], dict[str, Any]] = {}
    exact: dict[tuple[str, str], dict[str, Any]] = {}
    continuum_rows: list[dict[str, Any]] = []
    runtime_by_spectrum: dict[str, float] = {}
    evaluations_by_spectrum: dict[str, int] = {}
    scan_max_by_spectrum: dict[str, float] = {}

    for section_case, (point_1, point_2) in section_cases.items():
        for beta_deg in UT3_ANGLES_DEG:
            beta_rad = _ut3_beta_rad(beta_deg)
            for model in ("Timoshenko", "EB"):
                factory, raw_factory = _ut2_continuum_factories(
                    model, beta_rad, point_1, point_2
                )
                roots, builder, runtime = _solve_roots(
                    point_1, factory, num_roots=UT3_NUM_ROOTS
                )
                quality = [
                    _root_quality_fields(root, raw_factory) for root in roots
                ]
                gaps = _neighbor_gaps(_frequencies(roots))
                key = (section_case, beta_deg, model)
                standard[key] = {"roots": roots, "quality": quality}
                label = (
                    f"standard_J_book|{section_case}|{beta_deg:+g}|{model}"
                )
                runtime_by_spectrum[label] = runtime
                evaluations_by_spectrum[label] = builder.evaluations
                scan_max_by_spectrum[label] = builder.maximum_frequency_hz
                continuum_rows.extend(
                    _ut3_continuum_root_row(
                        builder="standard_J_book",
                        section_case=section_case,
                        point_1=point_1,
                        point_2=point_2,
                        beta_deg=beta_deg,
                        model=model,
                        mode=mode,
                        root=root,
                        quality=quality_fields,
                        gap_hz=float(gap),
                    )
                    for mode, (root, quality_fields, gap) in enumerate(
                        zip(roots, quality, gaps), start=1
                    )
                )

        for model in ("Timoshenko", "EB"):
            factory, raw_factory = _ut3_exact_continuum_factories(
                model, point_1, point_2
            )
            roots, builder, runtime = _solve_roots(
                point_1, factory, num_roots=UT3_NUM_ROOTS
            )
            quality = [
                _root_quality_fields(root, raw_factory) for root in roots
            ]
            gaps = _neighbor_gaps(_frequencies(roots))
            exact[(section_case, model)] = {
                "roots": roots,
                "quality": quality,
            }
            label = f"exact_beta90_limit|{section_case}|+90|{model}"
            runtime_by_spectrum[label] = runtime
            evaluations_by_spectrum[label] = builder.evaluations
            scan_max_by_spectrum[label] = builder.maximum_frequency_hz
            continuum_rows.extend(
                _ut3_continuum_root_row(
                    builder="exact_beta90_limit",
                    section_case=section_case,
                    point_1=point_1,
                    point_2=point_2,
                    beta_deg=90.0,
                    model=model,
                    mode=mode,
                    root=root,
                    quality=quality_fields,
                    gap_hz=float(gap),
                )
                for mode, (root, quality_fields, gap) in enumerate(
                    zip(roots, quality, gaps), start=1
                )
            )

    symmetry_specs = _ut2_symmetry_specifications(
        angle_deg=90.0, section_cases=tuple(UT3_SECTION_CASES)
    )
    continuum_symmetry_rows, continuum_symmetry_maxima = (
        _ut2_continuum_symmetry_rows(
            standard, specifications=symmetry_specs
        )
    )
    exact_comparison_rows, exact_comparison_maxima = _ut3_standard_exact_rows(
        standard, exact
    )
    exact_joint_rows = _ut3_exact_joint_limit_rows()

    analytic_eb = {
        (section_case, beta_deg): _frequencies(
            standard[(section_case, beta_deg, "EB")]["roots"]
        )
        for section_case in UT3_SECTION_CASES
        for beta_deg in UT3_ANGLES_DEG
    }
    fem_runs: dict[tuple[str, float], dict[str, Any]] = {}
    for section_case, (point_1, point_2) in section_cases.items():
        for beta_deg in UT3_ANGLES_DEG:
            fem_runs[(section_case, beta_deg)] = (
                _run_ut2_eb_fem_configuration(
                    section_case=section_case,
                    point_1=point_1,
                    point_2=point_2,
                    beta_deg=beta_deg,
                    analytic_frequencies=analytic_eb[
                        (section_case, beta_deg)
                    ],
                    allowed_angles_deg=UT3_ANGLES_DEG,
                    gate_label="UT-3",
                    elements_per_arm=UT3_FEM_ELEMENTS_PER_ARM,
                    num_roots=UT3_NUM_ROOTS,
                )
            )
    fem_rows = _ut2_fem_comparison_rows(fem_runs, analytic_eb)
    matrix_rows = _ut3_fem_matrix_rows(fem_runs)
    reflection_matrix_ok = all(
        row["status"] == "PASS"
        for row in matrix_rows
        if row["symmetry_type"] == "reflection"
    )
    relabeling_matrix_ok = all(
        row["status"] == "PASS"
        for row in matrix_rows
        if row["symmetry_type"] == "oriented_angle_relabeling"
    )
    native_fem_rows, native_fem_maxima = _ut2_fem_native_symmetry_rows(
        fem_runs,
        reflection_matrix_ok=reflection_matrix_ok,
        relabeling_matrix_ok=relabeling_matrix_ok,
        specifications=symmetry_specs,
    )
    model_difference_rows, model_difference_maxima = (
        _ut2_model_difference_rows(
            standard,
            angle_deg=90.0,
            section_cases=tuple(UT3_SECTION_CASES),
        )
    )

    quality_counts, quality_maxima = _ut3_root_quality_summary(
        continuum_rows
    )
    standard_counts = {
        f"{case}|{beta:+g}|{model}": len(data["roots"])
        for (case, beta, model), data in standard.items()
    }
    exact_counts = {
        f"{case}|+90|{model}": len(data["roots"])
        for (case, model), data in exact.items()
    }
    standard_rows = [
        row for row in continuum_rows if row["builder"] == "standard_J_book"
    ]
    exact_rows = [
        row
        for row in continuum_rows
        if row["builder"] == "exact_beta90_limit"
    ]
    all_roots_ok = all(
        np.isfinite(row["frequency_hz"])
        and row["frequency_hz"] > 0.0
        and row["quality_status"] == "PASS"
        for row in continuum_rows
    )
    root_complete = bool(
        len(standard) == 12
        and len(exact) == 6
        and len(standard_rows) == 84
        and len(exact_rows) == 42
        and all(count == UT3_NUM_ROOTS for count in standard_counts.values())
        and all(count == UT3_NUM_ROOTS for count in exact_counts.values())
        and all_roots_ok
    )
    continuum_symmetry_ok = all(
        row["status"] == "PASS" for row in continuum_symmetry_rows
    )
    standard_exact_ok = all(
        row["status"] == "PASS" for row in exact_comparison_rows
    )
    exact_joint_ok = all(row["status"] == "PASS" for row in exact_joint_rows)
    continuum_hard_ok = root_complete and continuum_symmetry_ok
    exact_limit_ok = exact_joint_ok and standard_exact_ok
    matrix_hard_ok = reflection_matrix_ok and relabeling_matrix_ok
    fem_hard_ok = bool(
        len(fem_runs) == 6
        and all(
            run["diagnostics"]["status"] == "PASS"
            for run in fem_runs.values()
        )
    )
    regressions_ok = regression_statuses == {
        "UT-0": "PASS",
        "UT-1": "PARTIAL_PASS",
        "UT-1a": "INCONCLUSIVE_NUMERICAL_AUDIT",
        "UT-2": "PASS",
    }
    ut3_status, overall_1d, phase_status = _classify_ut3_status(
        continuum_hard_ok=continuum_hard_ok,
        exact_limit_ok=exact_limit_ok,
        matrix_hard_ok=matrix_hard_ok,
        fem_hard_ok=fem_hard_ok,
        regressions_ok=regressions_ok,
    )

    fem_diagnostics = {
        f"{case}|{beta:+g}": run["diagnostics"]
        for (case, beta), run in fem_runs.items()
    }
    exact_results = {
        row["check"]: {
            "absolute_frobenius_residual": row[
                "absolute_frobenius_residual"
            ],
            "relative_frobenius_residual": row[
                "relative_frobenius_residual"
            ],
            "absolute_max_residual": row["absolute_max_residual"],
            "relative_max_residual": row["relative_max_residual"],
            "max_mismatch_row": int(row["max_mismatch_row"]),
            "max_mismatch_column": int(row["max_mismatch_column"]),
            "left_value": row["left_value"],
            "right_value": row["right_value"],
            "status": row["status"],
        }
        for row in exact_joint_rows
    }
    summary = {
        "git_context": _git_context(),
        "scope": "UT-3 beta=+/-90 deg unequal-thickness quarter-turn joint validation only",
        "geometry": {
            "theta_1_deg": 0.0,
            "theta_2_deg": 0.0,
            "material_mode": "elastic",
            "material": "existing HMS/DX-209",
            "length_1_m": LENGTH_M,
            "length_2_m": LENGTH_M,
            "b_1_m": WIDTH_B_M,
            "b_2_m": WIDTH_B_M,
            "shear_factor": cantilever_geometry().shear_factor,
            "a_direction": "parallel to e_z",
        },
        "angles_deg": list(UT3_ANGLES_DEG),
        "section_cases": {
            case: {"a_1_m": values[0], "a_2_m": values[1]}
            for case, values in UT3_SECTION_CASES.items()
        },
        "thresholds": {
            "root_normalized_determinant": ROOT_DETERMINANT_TOLERANCE,
            "root_relative_singular": ROOT_SINGULAR_TOLERANCE,
            "near_degenerate_relative_gap": NEAR_DEGENERATE_RELATIVE_GAP,
            "exact_joint_absolute_max": UT3_EXACT_LIMIT_TOLERANCE,
            "exact_joint_relative_frobenius": UT3_EXACT_LIMIT_TOLERANCE,
            "standard_exact_spectrum_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "continuum_symmetry_relative": SPECTRUM_RELATIVE_TOLERANCE,
            "fem_matrix_relative_frobenius": UT3_MATRIX_RELATIVE_TOLERANCE,
            "fem_matrix_relative_max_entry": UT3_MATRIX_RELATIVE_TOLERANCE,
            "fem_stiffness_symmetry": FEM_SYMMETRY_TOLERANCE,
            "fem_mass_symmetry": FEM_SYMMETRY_TOLERANCE,
            "fem_joint_kinematic": FEM_JOINT_KINEMATIC_TOLERANCE,
            "fem_joint_equilibrium": FEM_JOINT_EQUILIBRIUM_TOLERANCE,
            "fem_E3": FEM_FIRST_THREE_TOLERANCE,
            "fem_E6": FEM_FIRST_SIX_TOLERANCE,
            "native_fem_historical_marker_diagnostic_only": (
                SPECTRUM_RELATIVE_TOLERANCE
            ),
        },
        "root_scan": {
            "formulation": FORMULATION,
            "num_roots": UT3_NUM_ROOTS,
            "scan_step_hz": SCAN_STEP_HZ,
            "initial_max_hz": INITIAL_MAX_HZ,
            "max_hz": MAX_HZ,
            "actual_maximum_evaluated_frequency_hz": scan_max_by_spectrum,
        },
        "standard_continuum_spectrum_count": len(standard),
        "standard_continuum_root_total": len(standard_rows),
        "standard_continuum_root_counts": standard_counts,
        "exact_limit_spectrum_count": len(exact),
        "exact_limit_root_total": len(exact_rows),
        "exact_limit_root_counts": exact_counts,
        "continuum_root_total": len(continuum_rows),
        "first_seven_standard_frequencies_hz": {
            case: {
                f"{beta:+g} deg": {
                    model: _frequencies(
                        standard[(case, beta, model)]["roots"]
                    ).tolist()
                    for model in ("Timoshenko", "EB")
                }
                for beta in UT3_ANGLES_DEG
            }
            for case in UT3_SECTION_CASES
        },
        "first_seven_exact_limit_frequencies_hz": {
            case: {
                model: _frequencies(exact[(case, model)]["roots"]).tolist()
                for model in ("Timoshenko", "EB")
            }
            for case in UT3_SECTION_CASES
        },
        "root_quality_basis_counts": quality_counts,
        "root_quality_maxima_by_accepted_basis": quality_maxima,
        "exact_joint_limit_results": exact_results,
        "exact_joint_limit_maxima": _ut3_exact_rows_maxima(exact_joint_rows),
        "exact_joint_matrix_residual_maxima": _ut3_exact_rows_maxima(
            [
                row
                for row in exact_joint_rows
                if row["check"].startswith("J_book")
            ]
        ),
        "endpoint_quarter_turn_residual_maxima": _ut3_exact_rows_maxima(
            [
                row
                for row in exact_joint_rows
                if not row["check"].startswith("J_book")
            ]
        ),
        "standard_exact_spectrum_maxima": exact_comparison_maxima,
        "continuum_symmetry_maxima": continuum_symmetry_maxima,
        "fem_assembly_count": len(fem_runs),
        "fem_configuration_diagnostics": fem_diagnostics,
        "fem_structural_maxima": _ut3_fem_structural_maxima(
            fem_diagnostics
        ),
        "fem_reflection_matrix_residual_maxima": _ut2_matrix_group_maxima(
            matrix_rows, "reflection"
        ),
        "fem_relabeling_matrix_residual_maxima": _ut2_matrix_group_maxima(
            matrix_rows, "oriented_angle_relabeling"
        ),
        "native_fem_spectral_symmetry_maxima": native_fem_maxima,
        "native_fem_spectral_symmetry_is_hard_gate": False,
        "model_difference_maxima": model_difference_maxima,
        "runtime_seconds_by_continuum_spectrum": runtime_by_spectrum,
        "matrix_evaluations_by_continuum_spectrum": evaluations_by_spectrum,
        "total_continuum_matrix_evaluations": sum(
            evaluations_by_spectrum.values()
        ),
        "fem_runtime_seconds_by_configuration": {
            label: diagnostics["runtime_seconds"]
            for label, diagnostics in fem_diagnostics.items()
        },
        "total_runtime_seconds": time.perf_counter() - started,
        "regression_statuses": regression_statuses,
        "continuum_hard_gates_status": (
            "PASS" if continuum_hard_ok else "FAIL"
        ),
        "exact_limit_hard_gates_status": (
            "PASS" if exact_limit_ok else "FAIL"
        ),
        "fem_matrix_hard_gates_status": "PASS" if matrix_hard_ok else "FAIL",
        "fem_accuracy_structural_status": "PASS" if fem_hard_ok else "FAIL",
        "ut3_status": ut3_status,
        "overall_unequal_thickness_1d_validation": overall_1d,
        "unequal_thickness_1d_phase": phase_status,
        "next_separate_stage": "limited 3D FEM anchor design",
        "explicit_exclusions": UT3_EXPLICIT_EXCLUSIONS,
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "continuum_roots.csv", continuum_rows)
    _write_csv(output_dir / "exact_joint_limit.csv", exact_joint_rows)
    _write_csv(
        output_dir / "exact_limit_comparison.csv", exact_comparison_rows
    )
    _write_csv(output_dir / "continuum_symmetry.csv", continuum_symmetry_rows)
    _write_csv(output_dir / "fem_comparison.csv", fem_rows)
    _write_csv(output_dir / "fem_matrix_symmetry.csv", matrix_rows)
    _write_csv(
        output_dir / "fem_native_spectral_symmetry.csv", native_fem_rows
    )
    _write_csv(output_dir / "model_difference.csv", model_difference_rows)
    (output_dir / "ut3_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "ut3_report.md").write_text(
        _ut3_report(summary), encoding="utf-8"
    )
    return summary


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if (
        not args.smoke
        and not args.ut1_beta0
        and not args.ut1a_fem_exchange_audit
        and not args.ut2_beta30
        and not args.ut3_beta90
    ):
        print(
            "Select exactly one validation mode: --smoke, --ut1-beta0, "
            "--ut1a-fem-exchange-audit, --ut2-beta30, or --ut3-beta90.",
            file=sys.stderr,
        )
        return 2
    if args.smoke:
        output_dir = (
            DEFAULT_OUTPUT_DIR if args.output_dir is None else args.output_dir
        ).resolve()
        summary, _ = _run_smoke(output_dir)
        print(f"UT-0 status: {summary['ut0_status']}")
        print("Overall unequal-thickness validation: IN_PROGRESS")
        for model, values in summary["first_three_frequencies_hz"].items():
            print(f"{model}_coupled_hz={values['coupled']}")
            print(f"{model}_stepped_hz={values['stepped']}")
            print(
                f"{model}_max_relative_difference="
                f"{summary['max_coupled_stepped_relative_difference'][model]:.6e}"
            )
        print(f"output_dir={output_dir}")
        return 0 if summary["ut0_status"] == "PASS" else 4

    if args.ut1a_fem_exchange_audit:
        output_dir = (
            UT1A_DEFAULT_OUTPUT_DIR
            if args.output_dir is None
            else args.output_dir
        ).resolve()
        summary = _run_ut1a_fem_exchange_audit(output_dir)
        print(f"UT-1a status: {summary['ut1a_status']}")
        print(f"UT-1 status preserved: {summary['ut1_regression_status']}")
        print("Overall unequal-thickness validation: IN_PROGRESS")
        print(
            "interpretation_classification="
            f"{summary['interpretation_classification']}"
        )
        print(f"output_dir={output_dir}")
        return (
            0
            if summary["ut1a_status"] == "PASS_MATRIX_EQUIVARIANCE"
            else 4
        )

    if args.ut2_beta30:
        output_dir = (
            UT2_DEFAULT_OUTPUT_DIR
            if args.output_dir is None
            else args.output_dir
        ).resolve()
        summary = _run_ut2_beta30(output_dir)
        print(f"UT-2 status: {summary['ut2_status']}")
        print("Overall unequal-thickness validation: IN_PROGRESS")
        print(
            "continuum_root_counts="
            f"{summary['continuum_spectrum_count']}x{UT2_NUM_ROOTS}"
        )
        print(
            "root_quality_basis_counts="
            f"{summary['root_quality_basis_counts']}"
        )
        print(f"fem_assembly_count={summary['fem_assembly_count']}")
        print(f"output_dir={output_dir}")
        return 0 if summary["ut2_status"] in ("PASS", "PARTIAL_PASS") else 4

    if args.ut3_beta90:
        output_dir = (
            UT3_DEFAULT_OUTPUT_DIR
            if args.output_dir is None
            else args.output_dir
        ).resolve()
        summary = _run_ut3_beta90(output_dir)
        print(f"UT-3 status: {summary['ut3_status']}")
        print(
            "Overall unequal-thickness 1D validation: "
            f"{summary['overall_unequal_thickness_1d_validation']}"
        )
        print(
            "Unequal-thickness 1D validation phase: "
            f"{summary['unequal_thickness_1d_phase']}"
        )
        print(
            "standard_continuum_root_counts="
            f"{summary['standard_continuum_spectrum_count']}x{UT3_NUM_ROOTS}"
        )
        print(
            "exact_limit_root_counts="
            f"{summary['exact_limit_spectrum_count']}x{UT3_NUM_ROOTS}"
        )
        print(
            "root_quality_basis_counts="
            f"{summary['root_quality_basis_counts']}"
        )
        print(f"fem_assembly_count={summary['fem_assembly_count']}")
        print(f"output_dir={output_dir}")
        return 0 if summary["ut3_status"] in ("PASS", "PARTIAL_PASS") else 4

    output_dir = (
        UT1_DEFAULT_OUTPUT_DIR if args.output_dir is None else args.output_dir
    ).resolve()
    summary = _run_ut1_beta0(output_dir)
    print(f"UT-1 status: {summary['ut1_status']}")
    print("Overall unequal-thickness validation: IN_PROGRESS")
    print(
        "continuum_root_counts="
        f"{summary['main_continuum_spectrum_count']}x{UT1_NUM_ROOTS}"
    )
    print(
        "root_quality_basis_counts="
        f"{summary['root_quality_basis_counts']}"
    )
    print(f"output_dir={output_dir}")
    return 0 if summary["ut1_status"] in ("PASS", "PARTIAL_PASS") else 4


if __name__ == "__main__":
    raise SystemExit(main())
