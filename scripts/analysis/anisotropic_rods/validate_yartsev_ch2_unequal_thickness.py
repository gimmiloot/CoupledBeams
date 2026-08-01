"""Run the UT-0 unequal-thickness section-routing and beta=0 smoke gate."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
from collections.abc import Callable, Iterable, Sequence
from pathlib import Path
from typing import Any

import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib.yartsev_ch2_coupled_rods import (  # noqa: E402
    coupled_boundary_matrix,
    coupled_boundary_matrix_raw,
    equilibrate_matrix,
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
    eb_clamp_matrix,
    eb_coupled_boundary_matrix,
    eb_coupled_boundary_matrix_raw,
    eb_state_matrix,
    eb_state_transfer_matrix,
    eb_straight_right_clamp_matrix,
)


DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "anisotropic_rods"
    / "yartsev_ch2_unequal_thickness_validation"
    / "ut0_smoke"
)
SECTION_A_M = {"a4": 0.004, "a5": 0.005, "a6": 0.006}
LENGTH_M = 0.2
WIDTH_B_M = 0.020
BETA_RAD = 0.0
DIAGNOSTIC_FREQUENCY_HZ = 500.0
NUM_ROOTS = 3
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
        self._cache[key] = matrix
        return matrix


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="run the only implemented UT-0 smoke mode",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def _make_section_points() -> dict[str, RodPoint]:
    """Build three independent points; no section field is replaced in place."""

    base = cantilever_geometry(LENGTH_M)
    material = hms_dx_209_material()
    points: dict[str, RodPoint] = {}
    for case, a_m in SECTION_A_M.items():
        geometry = Geometry(
            a=a_m,
            b=WIDTH_B_M,
            length=LENGTH_M,
            shear_factor=base.shear_factor,
        )
        points[case] = make_rod_point(
            0.0,
            material_mode="elastic",
            geometry=geometry,
            material=material,
            series_relative_tolerance=TORSION_SERIES_RELATIVE_TOLERANCE,
        )
    return points


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
    point: RodPoint, factory: Callable[[complex], np.ndarray]
) -> tuple[list[RootResult], CountingBoundaryBuilder, float]:
    builder = CountingBoundaryBuilder(factory)
    started = time.perf_counter()
    roots = find_elastic_roots(
        point,
        FORMULATION,
        num_roots=NUM_ROOTS,
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


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.smoke:
        print(
            "UT-0 implements only smoke mode; rerun with --smoke.",
            file=sys.stderr,
        )
        return 2
    output_dir = args.output_dir.resolve()
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


if __name__ == "__main__":
    raise SystemExit(main())
