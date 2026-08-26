from __future__ import annotations

import ast
import csv
import json
import math
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal
from scipy.optimize import brentq

from scripts.analysis.laminated_beams import (
    reddy_four_ply_isotropic_postprocess as iso_postprocess,
)
from scripts.analysis.laminated_beams import (
    validate_reddy_four_ply_isotropic_limit as iso_audit,
)
from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as legacy_rect
from scripts.lib import reddy_symmetric_coupled_beams as coupled_rlb
from scripts.lib import reddy_symmetric_laminated_beam as rlb
from scripts.lib import variable_length_timoshenko as legacy_circular


ROOT = Path(__file__).resolve().parents[1]
CONTRACT_PATH = ROOT / "tests" / "data" / "reddy_four_ply_isotropic_limit_cases.json"
COMPARATOR_PATH = ROOT / "scripts" / "lib" / "isotropic_rectangular_timoshenko_coupled_beams.py"
RESULT_DIRECTORY = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_four_ply_isotropic_limit_validation"
)

ANGLE_TOLERANCE = 1.0e-12
SECTION_TOLERANCE = 1.0e-11
REGULAR_LEGACY_TOLERANCE = 1.0e-11
CUTOFF_LEGACY_TOLERANCE = 1.0e-8
REGULAR_LOCAL_TOLERANCE = 1.0e-10
CUTOFF_LOCAL_TOLERANCE = 1.0e-8
CIRCULAR_ROOT_SCAN_START = 0.2
CIRCULAR_ROOT_SCAN_STOP = 45.0
CIRCULAR_ROOT_SCAN_STEP = 0.005
CIRCULAR_ROOT_COUNT = 6
CIRCULAR_ROOT_RELATIVE_TOLERANCE = 1.0e-11

PRIMARY_STACK = (0.0, 90.0, 90.0, 0.0)
CONTROL_STACKS = (
    (0.0, 0.0, 0.0, 0.0),
    (17.0, -38.0, -38.0, 17.0),
)
ALL_STACKS = (PRIMARY_STACK, *CONTROL_STACKS)
ISOTROPIC_ANGLES_DEG = (0.0, 17.0, 45.0, 90.0, -38.0)


def _contract() -> dict[str, object]:
    return json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))


def _material(contract: dict[str, object]) -> rlb.LaminaMaterial:
    values = contract["material"]
    assert isinstance(values, dict)
    E = float(values["E"])
    nu = float(values["nu"])
    shear_modulus = E / (2.0 * (1.0 + nu))
    return rlb.LaminaMaterial(
        E1=E,
        E2=E,
        nu12=nu,
        G12=shear_modulus,
        G13=shear_modulus,
        G23=shear_modulus,
        rho=float(values["rho"]),
        name="four-ply isotropic-limit material",
    )


def _geometry(contract: dict[str, object], geometry_id: str) -> dict[str, float]:
    geometries = contract["geometries"]
    assert isinstance(geometries, dict)
    geometry = geometries[geometry_id]
    assert isinstance(geometry, dict)
    return {name: float(value) for name, value in geometry.items()}


def _four_ply_laminate(
    material: rlb.LaminaMaterial,
    stack: tuple[float, float, float, float],
    thickness: float,
) -> rlb.LaminateSection:
    ply_thickness = thickness / 4.0
    return rlb.integrate_laminate(
        tuple(rlb.Ply(material, angle, ply_thickness) for angle in stack)
    )


def _matrix_relative_residual(actual: np.ndarray, expected: np.ndarray) -> float:
    numerator = float(np.linalg.norm(np.asarray(actual) - np.asarray(expected), ord="fro"))
    denominator = max(float(np.linalg.norm(np.asarray(expected), ord="fro")), np.finfo(float).tiny)
    return numerator / denominator


def _relative_residual(actual: float, expected: float) -> float:
    return abs(float(actual) - float(expected)) / max(abs(float(expected)), np.finfo(float).tiny)


def _endpoint_state(columns: dict[str, np.ndarray]) -> np.ndarray:
    """Return a 6-by-3 endpoint matrix in legacy/RLB common state-name order."""

    return np.vstack(
        [
            np.asarray(columns[name], dtype=float)
            for name in ("u", "w", "psi", "N", "Q", "M")
        ]
    )


def _column_space_projector(matrix: np.ndarray) -> np.ndarray:
    values = np.asarray(matrix, dtype=float)
    column_norms = np.linalg.norm(values, axis=0)
    assert np.all(np.isfinite(column_norms)) and np.all(column_norms > 0.0)
    basis, triangular = np.linalg.qr(
        values / column_norms[None, :], mode="reduced"
    )
    assert np.linalg.matrix_rank(triangular) == matrix.shape[1]
    return basis @ basis.T


def _row_normalized_determinant(matrix: np.ndarray) -> float:
    values = np.asarray(matrix, dtype=float)
    norms = np.linalg.norm(values, axis=1)
    assert np.all(np.isfinite(norms)) and np.all(norms > 0.0)
    return float(np.linalg.det(values / norms[:, None]))


def test_frozen_json_contract_contains_exactly_the_eight_declared_cases() -> None:
    contract = _contract()
    assert contract["contract_version"] == "rlb_1c_iso_v1"
    assert contract["spectral_stack_deg"] == [0.0, 90.0, 90.0, 0.0]
    assert contract["algebraic_control_stacks_deg"] == [
        [0.0, 0.0, 0.0, 0.0],
        [17.0, -38.0, -38.0, 17.0],
    ]
    assert contract["material"] == {
        "E": 1.0,
        "rho": 1.0,
        "nu": 0.3,
        "K": 5.0 / 6.0,
    }
    assert contract["geometries"] == {
        "G20": {
            "width": 0.2,
            "thickness": 0.05,
            "width_to_thickness": 4.0,
            "L_ref_to_thickness": 20.0,
        },
        "G10": {
            "width": 0.4,
            "thickness": 0.1,
            "width_to_thickness": 4.0,
            "L_ref_to_thickness": 10.0,
        },
    }
    expected_cases = [
        {"case_id": "ISO-01", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": 0.0},
        {"case_id": "ISO-02", "geometry": "G20", "L1": 0.7, "L2": 1.3, "beta_deg": 0.0},
        {"case_id": "ISO-03", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": 30.0},
        {"case_id": "ISO-04", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": 90.0},
        {"case_id": "ISO-05", "geometry": "G20", "L1": 1.0, "L2": 1.0, "beta_deg": -30.0},
        {"case_id": "ISO-06", "geometry": "G20", "L1": 0.7, "L2": 1.3, "beta_deg": 30.0},
        {"case_id": "ISO-07", "geometry": "G20", "L1": 1.3, "L2": 0.7, "beta_deg": 30.0},
        {"case_id": "ISO-08", "geometry": "G10", "L1": 1.0, "L2": 1.0, "beta_deg": 30.0},
    ]
    cases = contract["cases"]
    assert isinstance(cases, list)
    assert len(cases) == 8
    assert cases == expected_cases
    assert len({case["case_id"] for case in cases}) == 8


def test_isotropic_lamina_reciprocity_q_positive_definiteness_and_angle_invariance() -> None:
    contract = _contract()
    material = _material(contract)
    values = contract["material"]
    assert isinstance(values, dict)
    E = float(values["E"])
    nu = float(values["nu"])
    shear_modulus = E / (2.0 * (1.0 + nu))
    expected_q = np.array(
        [
            [E / (1.0 - nu**2), nu * E / (1.0 - nu**2), 0.0],
            [nu * E / (1.0 - nu**2), E / (1.0 - nu**2), 0.0],
            [0.0, 0.0, shear_modulus],
        ]
    )
    expected_shear = shear_modulus * np.eye(2)

    assert material.nu21 == pytest.approx(nu, rel=0.0, abs=2.0e-16)
    assert material.nu12 / material.E1 == pytest.approx(
        material.nu21 / material.E2, rel=0.0, abs=2.0e-16
    )
    assert_allclose(material.Q, expected_q, rtol=2.0e-15, atol=0.0)
    assert_allclose(material.Q, material.Q.T, rtol=0.0, atol=0.0)
    assert np.min(np.linalg.eigvalsh(material.Q)) > 0.0
    assert np.min(np.linalg.eigvalsh(material.Q_shear)) > 0.0

    for angle_deg in ISOTROPIC_ANGLES_DEG:
        qbar = rlb.transformed_reduced_stiffness(material, angle_deg)
        shear_qbar = rlb.transformed_transverse_shear_stiffness(material, angle_deg)
        assert _matrix_relative_residual(qbar, expected_q) <= ANGLE_TOLERANCE
        assert _matrix_relative_residual(shear_qbar, expected_shear) <= ANGLE_TOLERANCE
        assert_allclose(qbar, qbar.T, rtol=0.0, atol=2.0e-16)
        assert_allclose(shear_qbar, shear_qbar.T, rtol=0.0, atol=2.0e-16)


@pytest.mark.parametrize("geometry_id", ["G20", "G10"])
def test_four_equal_plies_all_control_stacks_match_analytic_abd_shear_and_mass(
    geometry_id: str,
) -> None:
    contract = _contract()
    material = _material(contract)
    geometry = _geometry(contract, geometry_id)
    thickness = geometry["thickness"]
    density = material.rho
    q = material.Q
    shear_q = material.Q_shear
    expected_interfaces = np.asarray(
        [-thickness / 2.0, -thickness / 4.0, 0.0, thickness / 4.0, thickness / 2.0]
    )
    expected_a = thickness * q
    expected_d = thickness**3 * q / 12.0
    expected_shear = thickness * shear_q
    expected_i0 = density * thickness
    expected_i2 = density * thickness**3 / 12.0

    laminates = [
        _four_ply_laminate(material, stack, thickness)
        for stack in ALL_STACKS
    ]
    for laminate in laminates:
        assert len(laminate.plies) == 4
        assert_array_equal(
            np.asarray([ply.thickness for ply in laminate.plies]),
            np.full(4, thickness / 4.0),
        )
        assert_array_equal(laminate.z_interfaces, expected_interfaces)
        assert _matrix_relative_residual(laminate.A, expected_a) <= ANGLE_TOLERANCE
        assert _matrix_relative_residual(laminate.D, expected_d) <= ANGLE_TOLERANCE
        assert _matrix_relative_residual(laminate.shear, expected_shear) <= ANGLE_TOLERANCE
        assert np.linalg.norm(laminate.B, ord="fro") / (
            np.linalg.norm(q, ord="fro") * thickness**2
        ) <= ANGLE_TOLERANCE
        assert _relative_residual(laminate.I0, expected_i0) <= ANGLE_TOLERANCE
        assert _relative_residual(laminate.I2, expected_i2) <= ANGLE_TOLERANCE
        assert abs(laminate.I1) / (density * thickness**2) <= ANGLE_TOLERANCE
        assert rlb.check_laminate_symmetry(laminate).is_symmetric

    reference = laminates[0]
    for control in laminates[1:]:
        for name in ("A", "B", "D", "shear"):
            difference = np.asarray(getattr(control, name)) - np.asarray(
                getattr(reference, name)
            )
            if name == "B":
                residual = float(np.linalg.norm(difference, ord="fro")) / (
                    float(np.linalg.norm(q, ord="fro")) * thickness**2
                )
            else:
                residual = _matrix_relative_residual(
                    getattr(control, name), getattr(reference, name)
                )
            assert residual <= ANGLE_TOLERANCE
        mass_scales = {
            "I0": density * thickness,
            "I1": density * thickness**2,
            "I2": density * thickness**3,
        }
        for name, scale in mass_scales.items():
            assert (
                abs(float(getattr(control, name)) - float(getattr(reference, name)))
                / scale
                <= ANGLE_TOLERANCE
            )


@pytest.mark.parametrize("geometry_id", ["G20", "G10"])
def test_four_ply_reduction_equals_rectangular_ea_ei_kga_rhoa_and_rhoi(
    geometry_id: str,
) -> None:
    contract = _contract()
    material = _material(contract)
    material_contract = contract["material"]
    assert isinstance(material_contract, dict)
    geometry = _geometry(contract, geometry_id)
    E = float(material_contract["E"])
    nu = float(material_contract["nu"])
    density = float(material_contract["rho"])
    correction = float(material_contract["K"])
    width = geometry["width"]
    thickness = geometry["thickness"]
    area = width * thickness
    second_moment = width * thickness**3 / 12.0
    shear_modulus = E / (2.0 * (1.0 + nu))
    expected = np.asarray(
        [
            E * area,
            E * second_moment,
            correction * shear_modulus * area,
            density * area,
            density * second_moment,
        ]
    )
    direct = legacy_rect.rectangular_section(
        E=E,
        nu=nu,
        rho=density,
        width=width,
        thickness=thickness,
        K=correction,
    )
    assert direct.area == pytest.approx(area, rel=0.0, abs=0.0)
    assert direct.I_y == pytest.approx(second_moment, rel=0.0, abs=0.0)
    assert direct.I_y != pytest.approx(thickness * width**3 / 12.0)
    assert direct.width == width and direct.thickness == thickness
    assert direct.section_kind == "rectangle_width_ey_thickness_ez"
    assert_allclose(
        [direct.EA, direct.EI, direct.KGA, direct.rhoA, direct.rhoI],
        expected,
        rtol=2.0e-15,
        atol=0.0,
    )

    reduced_sets: list[np.ndarray] = []
    for stack in ALL_STACKS:
        laminate = _four_ply_laminate(material, stack, thickness)
        properties = rlb.reduce_to_beam_properties(
            laminate,
            width=width,
            K=correction,
        )
        reduced = np.asarray(
            [properties.A, properties.D, properties.S, properties.m, properties.J]
        )
        reduced_sets.append(reduced)
        assert np.max(np.abs((reduced - expected) / expected)) <= SECTION_TOLERANCE
        for reduction in (
            properties.axial_reduction,
            properties.bending_reduction,
            properties.shear_reduction_before_K,
        ):
            assert reduction is not None
            assert reduction.relative_difference <= SECTION_TOLERANCE
    for reduced in reduced_sets[1:]:
        assert_allclose(reduced, reduced_sets[0], rtol=SECTION_TOLERANCE, atol=0.0)


def test_independent_comparator_has_no_forbidden_physics_imports() -> None:
    source = COMPARATOR_PATH.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(COMPARATOR_PATH))
    imported_modules: set[str] = set()
    dynamic_import_calls: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_modules.add(node.module or "")
        elif isinstance(node, ast.Call):
            if isinstance(node.func, ast.Name) and node.func.id == "__import__":
                dynamic_import_calls.append("__import__")
            if isinstance(node.func, ast.Attribute) and node.func.attr == "import_module":
                dynamic_import_calls.append("import_module")

    forbidden_fragments = (
        "reddy_symmetric_laminated_beam",
        "reddy_inplane_geometry",
        "reddy_symmetric_coupled_beams",
        "reddy_symmetric_coupled_beams_ritz",
        "variable_length_timoshenko",
        "analytic_coupled_rods_shapes",
        "yartsev",
        "shellbuckling",
        "scipy",
    )
    assert not dynamic_import_calls
    assert not {
        module
        for module in imported_modules
        if any(fragment.casefold() in module.casefold() for fragment in forbidden_fragments)
    }
    assert imported_modules <= {"__future__", "dataclasses", "math", "numpy", "numpy.typing"}


@pytest.mark.parametrize(
    ("cutoff_factor", "expected_regime"),
    [
        (0.8, legacy_rect.TIMO_REGIME_MIXED),
        (1.0, legacy_rect.TIMO_REGIME_CUTOFF),
        (1.2, legacy_rect.TIMO_REGIME_TWO_TRIG),
    ],
)
@pytest.mark.parametrize("beta_deg", [0.0, 30.0])
def test_generic_engine_matches_frozen_circular_basis_endpoints_and_coupled_matrix(
    cutoff_factor: float,
    expected_regime: str,
    beta_deg: float,
) -> None:
    epsilon = 0.05
    circular_reference = legacy_circular.section_from_epsilon(epsilon)
    circular_generic = legacy_rect.circular_section(
        E=legacy_circular.E,
        nu=legacy_circular.NU,
        rho=legacy_circular.RHO,
        radius=circular_reference.radius,
        K=circular_reference.kappa,
    )
    lambda_cutoff = legacy_circular.lambda_cutoff(epsilon, circular_reference)
    lambda_value = cutoff_factor * lambda_cutoff
    omega = legacy_circular.project_omega(lambda_value, epsilon)
    old_basis = legacy_circular.timo_basis(lambda_value, epsilon, circular_reference)
    new_basis = legacy_rect.timoshenko_spatial_basis(omega, circular_generic)
    tolerance = CUTOFF_LEGACY_TOLERANCE if cutoff_factor == 1.0 else REGULAR_LEGACY_TOLERANCE

    assert old_basis.regime == expected_regime
    assert new_basis.regime == expected_regime
    assert_allclose(
        [new_basis.a, new_basis.b, new_basis.h, new_basis.q, new_basis.z_a, new_basis.z_b, new_basis.alpha],
        [old_basis.a, old_basis.b, old_basis.h, old_basis.q, old_basis.z_a, old_basis.z_b, old_basis.alpha],
        rtol=tolerance,
        atol=tolerance * max(1.0, abs(old_basis.z_a), abs(old_basis.z_b)),
    )
    assert_allclose(
        legacy_rect.bending_column_scales(new_basis),
        legacy_circular.bending_column_scales(old_basis),
        rtol=tolerance,
        atol=tolerance,
    )

    for x in (-1.0, -0.31, 0.0, 0.47, 1.0):
        old_bending = legacy_circular.bending_endpoint_columns(x, old_basis)
        new_bending = legacy_rect.clamped_bending_columns(x, new_basis)
        for name in ("w", "psi", "w_prime", "psi_prime"):
            assert_allclose(
                new_bending[name], old_bending[name], rtol=tolerance, atol=tolerance
            )
        old_endpoint = legacy_circular.endpoint_columns(
            x, omega, old_basis, circular_reference
        )
        new_endpoint = legacy_rect.clamped_endpoint_columns(
            x, omega, circular_generic, basis=new_basis
        )
        for name in ("u", "w", "psi", "N", "Q", "M"):
            assert_allclose(
                new_endpoint[name], old_endpoint[name], rtol=tolerance, atol=tolerance
            )

    old_matrix, old_warnings = legacy_circular.timo_coupling_matrix(
        lambda_value,
        beta_deg,
        0.0,
        epsilon,
        0.0,
    )
    new_assembly = legacy_rect.assemble_legacy_coupled_boundary(
        omega,
        circular_generic,
        1.0,
        circular_generic,
        1.0,
        beta_deg=beta_deg,
    )
    assert new_assembly.warnings == old_warnings
    assert_allclose(new_assembly.matrix, old_matrix, rtol=tolerance, atol=tolerance)
    assert _row_normalized_determinant(new_assembly.matrix) == pytest.approx(
        _row_normalized_determinant(old_matrix), rel=tolerance, abs=tolerance
    )


def test_circular_adapter_first_six_roots_match_frozen_legacy_case() -> None:
    """Bounded adapter gate; this is not an isotropic RLB spectrum search."""

    epsilon = 0.05
    beta_deg = 30.0
    old_section = legacy_circular.section_from_epsilon(epsilon)
    generic_section = legacy_rect.circular_section(
        E=legacy_circular.E,
        nu=legacy_circular.NU,
        rho=legacy_circular.RHO,
        radius=old_section.radius,
        K=old_section.kappa,
    )

    def generic_determinant(lambda_value: float) -> float:
        omega = legacy_circular.project_omega(lambda_value, epsilon)
        return legacy_rect.legacy_coupled_characteristic_determinant(
            omega,
            generic_section,
            legacy_circular.L_SEGMENT,
            generic_section,
            legacy_circular.L_SEGMENT,
            beta_deg=beta_deg,
            scaled=True,
        )

    # This fixed grid is independent of the frozen root inventory.  Only sign
    # brackets of the generic comparator determinant enter the refinements.
    grid = np.arange(
        CIRCULAR_ROOT_SCAN_START,
        CIRCULAR_ROOT_SCAN_STOP + 0.5 * CIRCULAR_ROOT_SCAN_STEP,
        CIRCULAR_ROOT_SCAN_STEP,
    )
    determinant_values = np.array(
        [generic_determinant(float(value)) for value in grid], dtype=float
    )
    generic_roots: list[float] = []
    for left, right, value_left, value_right in zip(
        grid[:-1], grid[1:], determinant_values[:-1], determinant_values[1:]
    ):
        if not (math.isfinite(float(value_left)) and math.isfinite(float(value_right))):
            continue
        candidate: float | None = None
        if value_left == 0.0:
            candidate = float(left)
        elif np.signbit(value_left) != np.signbit(value_right):
            candidate = float(
                brentq(
                    generic_determinant,
                    float(left),
                    float(right),
                    xtol=1.0e-13,
                    rtol=1.0e-13,
                    maxiter=160,
                )
            )
        if candidate is None:
            continue
        if generic_roots and abs(candidate - generic_roots[-1]) <= 1.0e-7:
            continue
        generic_roots.append(candidate)
        if len(generic_roots) == CIRCULAR_ROOT_COUNT:
            break

    generic_root_array = np.asarray(generic_roots, dtype=float)
    assert generic_root_array.shape == (CIRCULAR_ROOT_COUNT,)
    assert np.all(np.isfinite(generic_root_array))
    assert np.all(generic_root_array > 0.0)

    frozen_roots, _warnings = legacy_circular.timo_sorted_roots(
        beta_deg,
        0.0,
        epsilon,
        CIRCULAR_ROOT_COUNT,
        eta=0.0,
    )
    frozen_roots = np.asarray(frozen_roots, dtype=float)
    assert frozen_roots.shape == (CIRCULAR_ROOT_COUNT,)
    assert np.all(np.isfinite(frozen_roots))
    relative_differences = np.abs(generic_root_array - frozen_roots) / np.abs(
        frozen_roots
    )
    assert float(np.max(relative_differences)) <= CIRCULAR_ROOT_RELATIVE_TOLERANCE


def _primary_rlb_and_rectangular_sections() -> tuple[
    rlb.BeamProperties, legacy_rect.SectionProperties
]:
    contract = _contract()
    material_contract = contract["material"]
    assert isinstance(material_contract, dict)
    material = _material(contract)
    geometry = _geometry(contract, "G20")
    laminate = _four_ply_laminate(material, PRIMARY_STACK, geometry["thickness"])
    properties = rlb.reduce_to_beam_properties(
        laminate,
        width=geometry["width"],
        K=float(material_contract["K"]),
    )
    section = legacy_rect.rectangular_section(
        E=float(material_contract["E"]),
        nu=float(material_contract["nu"]),
        rho=float(material_contract["rho"]),
        width=geometry["width"],
        thickness=geometry["thickness"],
        K=float(material_contract["K"]),
    )
    return properties, section


@pytest.mark.parametrize(
    "cutoff_factor", [0.23, 0.5, 0.73, 0.99, 1.0, 1.01, 1.5, 2.2]
)
def test_local_clamped_arm_endpoint_spaces_match_under_analytic_sign_and_coordinate_maps(
    cutoff_factor: float,
) -> None:
    properties, section = _primary_rlb_and_rectangular_sections()
    omega = cutoff_factor * legacy_rect.omega_cutoff(section)
    basis = legacy_rect.timoshenko_spatial_basis(omega, section)
    # This preliminary local-space gate is deliberately evaluated on a short,
    # fixed audit arm.  It spans all regimes without allowing exponentially
    # imbalanced hyperbolic columns to consume the digits needed by the strict
    # principal-angle check; it is not a coupled spectral calculation.
    length = 0.25
    at_zero = _endpoint_state(
        legacy_rect.clamped_endpoint_columns(0.0, omega, section, basis=basis)
    )
    at_positive = _endpoint_state(
        legacy_rect.clamped_endpoint_columns(length, omega, section, basis=basis)
    )
    at_negative = _endpoint_state(
        legacy_rect.clamped_endpoint_columns(-length, omega, section, basis=basis)
    )
    arm1_map = np.diag([1.0, 1.0, -1.0, 1.0, 1.0, -1.0])
    arm2_map = np.diag([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])
    transfer = rlb.combined_transfer_matrix(omega, length, properties)
    tolerance = (
        CUTOFF_LOCAL_TOLERANCE
        if 0.99 <= cutoff_factor <= 1.01
        else REGULAR_LOCAL_TOLERANCE
    )

    arm1_initial = arm1_map @ at_zero
    arm1_endpoint = arm1_map @ at_positive
    arm2_initial = arm2_map @ at_zero
    arm2_endpoint = arm2_map @ at_negative
    state_scale = rlb.combined_state_scale(properties, length)
    inverse_state_scale = np.diag(1.0 / np.diag(state_scale))
    clamp_basis = np.zeros((6, 3))
    clamp_basis[3:, :] = np.eye(3)
    rlb_endpoint_space = inverse_state_scale @ transfer @ clamp_basis
    assert_allclose(
        transfer @ arm1_initial,
        arm1_endpoint,
        rtol=tolerance,
        atol=tolerance * max(1.0, float(np.max(np.abs(arm1_endpoint)))),
    )
    assert_allclose(
        transfer @ arm2_initial,
        arm2_endpoint,
        rtol=tolerance,
        atol=tolerance * max(1.0, float(np.max(np.abs(arm2_endpoint)))),
    )
    for legacy_endpoint_space, legacy_initial_state in (
        (inverse_state_scale @ arm1_endpoint, arm1_initial),
        (inverse_state_scale @ arm2_endpoint, arm2_initial),
    ):
        assert np.linalg.matrix_rank(rlb_endpoint_space) == 3
        assert np.linalg.matrix_rank(legacy_endpoint_space) == 3
        # At the common outer clamp the kinematic rows vanish.  The lower
        # 3-by-3 block is therefore the exact, analytically oriented map from
        # legacy coefficients to the canonical RLB reaction basis.  Applying
        # its inverse before QR avoids losing the mixed-regime subspace when
        # the two hyperbolic endpoint columns are severely imbalanced.
        column_map = np.asarray(legacy_initial_state[3:, :], dtype=float)
        assert np.linalg.matrix_rank(column_map) == 3
        assert math.isfinite(float(np.linalg.cond(column_map)))
        mapping_residual = np.linalg.norm(
            rlb_endpoint_space @ column_map - legacy_endpoint_space, ord="fro"
        ) / np.linalg.norm(legacy_endpoint_space, ord="fro")
        assert mapping_residual <= tolerance
        best_column_map, _residuals, fit_rank, _singular_values = np.linalg.lstsq(
            rlb_endpoint_space, legacy_endpoint_space, rcond=None
        )
        assert fit_rank == 3
        assert np.linalg.matrix_rank(best_column_map) == 3
        assert math.isfinite(float(np.linalg.cond(best_column_map)))
        best_mapping_residual = np.linalg.norm(
            rlb_endpoint_space @ best_column_map - legacy_endpoint_space,
            ord="fro",
        ) / np.linalg.norm(legacy_endpoint_space, ord="fro")
        assert best_mapping_residual <= tolerance
        assert _matrix_relative_residual(best_column_map, column_map) <= tolerance
        assert math.isfinite(float(np.linalg.cond(rlb_endpoint_space)))
        assert math.isfinite(float(np.linalg.cond(legacy_endpoint_space)))
        reconditioned_legacy_space = legacy_endpoint_space @ np.linalg.inv(
            best_column_map
        )
        rlb_projector = _column_space_projector(rlb_endpoint_space)
        legacy_projector = _column_space_projector(reconditioned_legacy_space)
        assert _matrix_relative_residual(rlb_projector, legacy_projector) <= tolerance
        rlb_basis, _ = np.linalg.qr(
            rlb_endpoint_space
            / np.linalg.norm(rlb_endpoint_space, axis=0)[None, :],
            mode="reduced",
        )
        legacy_basis, _ = np.linalg.qr(
            reconditioned_legacy_space
            / np.linalg.norm(reconditioned_legacy_space, axis=0)[None, :],
            mode="reduced",
        )
        overlap_singular_values = np.linalg.svd(
            rlb_basis.T @ legacy_basis, compute_uv=False
        )
        assert 1.0 - float(np.min(overlap_singular_values)) <= tolerance


def test_legacy_to_rlb_state_maps_satisfy_state_equations_and_boundary_virtual_work() -> None:
    properties, section = _primary_rlb_and_rectangular_sections()
    omega = 0.73 * legacy_rect.omega_cutoff(section)
    arm1_map = np.diag([1.0, 1.0, -1.0, 1.0, 1.0, -1.0])
    arm2_map = np.diag([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])
    old_state_matrix = np.array(
        [
            [0.0, 0.0, 0.0, 1.0 / section.EA, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 1.0 / section.KGA, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / section.EI],
            [-section.rhoA * omega**2, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -section.rhoA * omega**2, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -section.rhoI * omega**2, 0.0, -1.0, 0.0],
        ]
    )
    rlb_state_matrix = rlb.combined_state_matrix(omega, properties)
    assert_allclose(rlb_state_matrix @ arm1_map, arm1_map @ old_state_matrix, rtol=2.0e-14, atol=2.0e-14)
    assert_allclose(rlb_state_matrix @ arm2_map, -arm2_map @ old_state_matrix, rtol=2.0e-14, atol=2.0e-14)

    rng = np.random.default_rng(20260826)
    for _ in range(128):
        old_motion = rng.normal(size=3)
        old_resultants = rng.normal(size=3)
        arm1_motion = np.diag([1.0, 1.0, -1.0]) @ old_motion
        arm1_resultants = np.diag([1.0, 1.0, -1.0]) @ old_resultants
        arm2_motion = -old_motion
        arm2_resultants = old_resultants
        old_pairing = float(old_resultants @ old_motion)
        assert float(arm1_resultants @ arm1_motion) == pytest.approx(old_pairing, abs=2.0e-14)
        # The legacy arm-2 joint is the negative-coordinate endpoint.
        assert float(arm2_resultants @ arm2_motion) == pytest.approx(-old_pairing, abs=2.0e-14)


@pytest.mark.parametrize("beta_deg", [0.0, 30.0, 90.0, -30.0])
def test_analytic_state_maps_transform_legacy_joint_rows_to_frozen_rlb_rows(
    beta_deg: float,
) -> None:
    _properties, section = _primary_rlb_and_rectangular_sections()
    omega = 0.77 * legacy_rect.omega_cutoff(section)
    length_1, length_2 = 0.7, 1.3
    assembly = legacy_rect.assemble_legacy_coupled_boundary(
        omega,
        section,
        length_1,
        section,
        length_2,
        beta_deg=beta_deg,
    )
    endpoint_1 = _endpoint_state(
        legacy_rect.clamped_endpoint_columns(
            length_1, omega, section, basis=assembly.basis_1
        )
    )
    endpoint_2 = _endpoint_state(
        legacy_rect.clamped_endpoint_columns(
            -length_2, omega, section, basis=assembly.basis_2
        )
    )
    arm1_map = np.diag([1.0, 1.0, -1.0, 1.0, 1.0, -1.0])
    arm2_map = np.diag([-1.0, -1.0, -1.0, 1.0, 1.0, 1.0])
    coefficient_to_rlb_joint_state = np.zeros((12, 6))
    coefficient_to_rlb_joint_state[:6, :3] = arm1_map @ endpoint_1
    coefficient_to_rlb_joint_state[6:, 3:] = arm2_map @ endpoint_2
    rlb_rows = coupled_rlb.joint_matrix_closed_form(math.radians(beta_deg)) @ coefficient_to_rlb_joint_state
    old_rows = np.asarray(assembly.matrix)
    analytically_reordered_old_rows = np.vstack(
        [
            old_rows[1],
            old_rows[0],
            -old_rows[2],
            old_rows[5],
            -old_rows[4],
            -old_rows[3],
        ]
    )
    assert_allclose(rlb_rows, analytically_reordered_old_rows, rtol=2.0e-13, atol=2.0e-13)


def _synthetic_root_policy(
    *,
    requested_roots: int,
    guard_roots: int,
    Omega_min: float = 0.2,
    Omega_max: float = 2.8,
    primary_scan_points: int = 53,
    verification_scan_points: int = 105,
    cluster_atol_Omega: float = 1.0e-8,
    local_close_pair_guard_subintervals: int = 0,
) -> iso_audit.SearchPolicy:
    """Small deterministic policy for model-neutral root-engine regressions."""

    return iso_audit.SearchPolicy(
        requested_roots=requested_roots,
        guard_roots=guard_roots,
        Omega_min=Omega_min,
        Omega_max=Omega_max,
        primary_scan_points=primary_scan_points,
        verification_scan_points=verification_scan_points,
        primary_phases=(0.0,),
        verification_phases=(0.5,),
        sigma_prefilter=1.0e-5,
        root_singular_ratio=1.0e-9,
        nullity_relative_threshold=1.0e-12,
        boundary_null_residual=1.0e-9,
        root_xtol_Omega=1.0e-12,
        root_rtol=8.0 * np.finfo(float).eps,
        dedup_atol_Omega=1.0e-8,
        dedup_rtol=1.0e-10,
        cluster_atol_Omega=cluster_atol_Omega,
        cluster_rtol=0.0,
        post_guard_tail_Omega=0.2,
        local_close_pair_guard_subintervals=local_close_pair_guard_subintervals,
    )


def _synthetic_diagonal_provider(*entries: object):
    def provider(omega: float) -> np.ndarray:
        diagonal = [
            float(entry(omega)) if callable(entry) else float(entry)
            for entry in entries
        ]
        return np.diag(diagonal)

    return provider


def _synthetic_spectral_provider(*eigenvalues: object):
    """Dense symmetric matrix with prescribed frequency-dependent eigenvalues."""

    random = np.random.default_rng(20260827)
    orthogonal, _triangular = np.linalg.qr(random.normal(size=(6, 6)))

    def provider(omega: float) -> np.ndarray:
        diagonal = [
            float(entry(omega)) if callable(entry) else float(entry)
            for entry in eigenvalues
        ]
        return orthogonal @ np.diag(diagonal) @ orthogonal.T

    return provider


def test_synthetic_root_engine_finds_simple_sign_roots_and_serializes_guard(
    tmp_path: Path,
) -> None:
    policy = _synthetic_root_policy(requested_roots=1, guard_roots=1)
    provider = _synthetic_diagonal_provider(
        lambda value: (value - 1.0) * (value - 2.0),
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    inventory = iso_audit.seed_free_root_inventory(
        provider,
        1.0,
        policy,
        case_id="SYN-SIMPLE",
        builder_id="SYNTHETIC",
    )

    assert inventory.status == "PASS"
    assert inventory.independent_agreement
    assert inventory.guard_available
    assert inventory.guard_not_at_scan_boundary
    assert inventory.unresolved_low_sigma_count == 0
    assert len(inventory.slots) == 2
    assert_allclose(
        [slot.event.Omega for slot in inventory.slots],
        [1.0, 2.0],
        rtol=0.0,
        atol=2.0e-10,
    )
    assert [slot.role for slot in inventory.slots] == [
        "FIRST_12",
        "ROOT_13_GUARD",
    ]
    assert all(slot.event.detected_nullity == 1 for slot in inventory.slots)
    assert any(
        "determinant_bracket" in source
        for candidate in inventory.primary.candidates
        for source in candidate.detection_sources
    )

    root_rows = [
        iso_audit._root_row(inventory, slot) for slot in inventory.slots
    ]
    verification_rows = [
        iso_audit._candidate_row(candidate)
        for candidate in inventory.verification.candidates
    ]
    iso_audit._write_csv(tmp_path / "roots.csv", root_rows)
    iso_audit._write_csv(tmp_path / "verification_candidates.csv", verification_rows)
    with (tmp_path / "roots.csv").open(encoding="utf-8", newline="") as stream:
        serialized_roots = list(csv.DictReader(stream))
    with (tmp_path / "verification_candidates.csv").open(
        encoding="utf-8", newline=""
    ) as stream:
        serialized_verification = list(csv.DictReader(stream))
    assert [row["role"] for row in serialized_roots] == [
        "FIRST_12",
        "ROOT_13_GUARD",
    ]
    assert all(row["case_inventory_sha256"] for row in serialized_roots)
    assert serialized_verification
    assert {row["scan_id"] for row in serialized_verification} == {"verification"}


def test_synthetic_root_engine_detects_even_double_root_without_sign_change() -> None:
    policy = _synthetic_root_policy(requested_roots=2, guard_roots=1)
    provider = _synthetic_spectral_provider(
        lambda value: (value - 1.0) * (value - 2.0),
        lambda value: value - 1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    inventory = iso_audit.seed_free_root_inventory(
        provider,
        1.0,
        policy,
        case_id="SYN-DOUBLE",
        builder_id="SYNTHETIC",
    )

    assert inventory.status == "PASS"
    assert len(inventory.slots) == 3
    double_slots = inventory.slots[:2]
    assert_allclose(
        [slot.event.Omega for slot in double_slots],
        [1.0, 1.0],
        rtol=0.0,
        atol=2.0e-10,
    )
    assert all(slot.event.multiplicity == 2 for slot in double_slots)
    assert all(slot.event.detected_nullity == 2 for slot in double_slots)
    assert all(
        slot.event.cluster_semantics == "EXACT_DEGENERATE_SUBSPACE"
        for slot in double_slots
    )
    assert inventory.slots[2].role == "ROOT_13_GUARD"
    assert inventory.slots[2].event.Omega == pytest.approx(2.0, abs=2.0e-10)
    assert any(
        "sigma_ratio_minimum" in source
        for candidate in inventory.primary.candidates
        if candidate.diagnostics.detected_nullity == 2
        for source in candidate.detection_sources
    )


def test_synthetic_root_engine_groups_two_close_distinct_roots() -> None:
    policy = _synthetic_root_policy(
        requested_roots=2,
        guard_roots=1,
        primary_scan_points=261,
        verification_scan_points=521,
        cluster_atol_Omega=3.0e-2,
    )
    provider = _synthetic_spectral_provider(
        lambda value: value - 1.0,
        lambda value: value - 1.02,
        lambda value: value - 2.0,
        1.0,
        1.0,
        1.0,
    )
    inventory = iso_audit.seed_free_root_inventory(
        provider,
        1.0,
        policy,
        case_id="SYN-CLOSE",
        builder_id="SYNTHETIC",
    )

    assert inventory.status == "PASS"
    assert len(inventory.slots) == 3
    close_slots = inventory.slots[:2]
    assert_allclose(
        [slot.event.Omega for slot in close_slots],
        [1.0, 1.02],
        rtol=0.0,
        atol=2.0e-10,
    )
    assert close_slots[0].event.cluster_id == close_slots[1].event.cluster_id
    assert close_slots[0].event.cluster_id
    assert all(slot.event.multiplicity == 1 for slot in close_slots)
    assert all(slot.event.detected_nullity == 1 for slot in close_slots)
    assert all(slot.event.cluster_multiplicity == 2 for slot in close_slots)
    assert all(slot.event.cluster_total_nullity == 2 for slot in close_slots)
    assert all(
        slot.event.cluster_semantics == "NEAR_DEGENERATE_CLUSTER"
        for slot in close_slots
    )
    assert close_slots[0].event.cluster_center_Omega == pytest.approx(
        1.01, abs=2.0e-10
    )


def test_local_close_pair_guard_resolves_two_roots_inside_one_coarse_cell() -> None:
    policy = _synthetic_root_policy(
        requested_roots=2,
        guard_roots=1,
        primary_scan_points=27,
        verification_scan_points=53,
        cluster_atol_Omega=5.0e-2,
        local_close_pair_guard_subintervals=32,
    )
    provider = _synthetic_spectral_provider(
        lambda value: value - 1.03,
        lambda value: value - 1.06,
        lambda value: value - 2.0,
        1.0,
        1.0,
        1.0,
    )
    inventory = iso_audit.seed_free_root_inventory(
        provider,
        1.0,
        policy,
        case_id="SYN-CLOSE-PAIR-GUARD",
        builder_id="SYNTHETIC",
    )

    assert inventory.status == "PASS"
    assert_allclose(
        [slot.event.Omega for slot in inventory.slots],
        [1.03, 1.06, 2.0],
        rtol=0.0,
        atol=2.0e-10,
    )
    assert any(
        "determinant_close_pair_guard" in source
        for candidate in inventory.primary.candidates
        for source in candidate.detection_sources
    )


def test_cross_detector_duplicate_uses_parallel_nullspace_not_wider_frequency_tolerance() -> None:
    policy = _synthetic_root_policy(requested_roots=1, guard_roots=0)
    provider = _synthetic_spectral_provider(
        lambda value: value - 1.0,
        lambda value: value - 1.02,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    at_first = iso_audit._boundary_matrix_diagnostics(1.0, provider, 1.0, policy)
    at_second = iso_audit._boundary_matrix_diagnostics(1.02, provider, 1.0, policy)

    determinant_candidate = iso_audit.RootCandidate(
        case_id="SYN-DETECTION",
        builder_id="SYNTHETIC",
        scan_id="primary",
        Omega=1.0,
        detection_sources=("primary:determinant_bracket",),
        interval_left_Omega=0.9,
        interval_right_Omega=1.1,
        interior_minimum=True,
        diagnostics=at_first,
        accepted=True,
        rejection_reason="",
    )
    same_root_sigma_candidate = iso_audit.RootCandidate(
        case_id="SYN-DETECTION",
        builder_id="SYNTHETIC",
        scan_id="primary",
        Omega=1.0 + 1.0e-5,
        detection_sources=("primary:sigma_ratio_minimum",),
        interval_left_Omega=0.8,
        interval_right_Omega=1.2,
        interior_minimum=True,
        diagnostics=at_first,
        accepted=True,
        rejection_reason="",
    )
    distinct_root_sigma_candidate = iso_audit.RootCandidate(
        case_id="SYN-DETECTION",
        builder_id="SYNTHETIC",
        scan_id="primary",
        Omega=1.02,
        detection_sources=("primary:sigma_ratio_minimum",),
        interval_left_Omega=0.8,
        interval_right_Omega=1.2,
        interior_minimum=True,
        diagnostics=at_second,
        accepted=True,
        rejection_reason="",
    )

    assert not iso_audit._candidate_close(
        determinant_candidate, same_root_sigma_candidate, policy
    )
    assert iso_audit._same_root_detection(
        determinant_candidate, same_root_sigma_candidate, policy
    )
    assert not iso_audit._same_root_detection(
        determinant_candidate, distinct_root_sigma_candidate, policy
    )


def test_synthetic_false_sigma_valley_is_rejected_and_remains_auditable() -> None:
    policy = _synthetic_root_policy(
        requested_roots=1,
        guard_roots=0,
        Omega_min=0.5,
        Omega_max=1.5,
        primary_scan_points=101,
        verification_scan_points=201,
    )
    provider = _synthetic_spectral_provider(
        lambda value: 1.0e-6 + (value - 1.0) ** 2,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    scan = iso_audit._run_scan(
        provider,
        1.0,
        policy,
        case_id="SYN-FALSE-VALLEY",
        builder_id="SYNTHETIC",
        scan_id="primary",
        points=policy.primary_scan_points,
        phases=policy.primary_phases,
    )

    assert not scan.candidates
    false_valleys = [
        candidate
        for candidate in scan.rejected_candidates
        if candidate.rejection_reason == "FALSE_SIGMA_VALLEY"
    ]
    assert false_valleys
    assert min(
        candidate.diagnostics.scaled_sigma_ratio for candidate in false_valleys
    ) <= policy.sigma_prefilter
    assert all(not candidate.accepted for candidate in false_valleys)


def test_synthetic_strict_nullity_zero_candidate_cannot_silently_disappear() -> None:
    policy = _synthetic_root_policy(
        requested_roots=1,
        guard_roots=0,
        Omega_min=0.5,
        Omega_max=1.5,
        primary_scan_points=101,
        verification_scan_points=201,
    )
    provider = _synthetic_spectral_provider(
        lambda value: 1.0e-10 + (value - 1.0) ** 2,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    scan = iso_audit._run_scan(
        provider,
        1.0,
        policy,
        case_id="SYN-STRICT-NULLITY",
        builder_id="SYNTHETIC",
        scan_id="primary",
        points=policy.primary_scan_points,
        phases=policy.primary_phases,
    )

    assert not scan.candidates
    unresolved = [
        candidate
        for candidate in scan.rejected_candidates
        if candidate.rejection_reason == "STRICT_NULLITY_UNRESOLVED"
    ]
    assert unresolved
    assert all(candidate.diagnostics.root_gate_nullity == 1 for candidate in unresolved)
    assert all(candidate.diagnostics.detected_nullity == 0 for candidate in unresolved)
    assert all(not candidate.accepted for candidate in unresolved)


def test_positive_scaling_preserves_rank_determinant_sign_and_zero_set() -> None:
    raw = np.diag([1.0e6, 1.0e-4, -3.0, 4.0, 5.0, 6.0])
    scaling = iso_audit._positive_equilibrate_model_neutral(raw)
    scaled = scaling.scaled_matrix
    assert np.all(np.isfinite(scaling.row_factors))
    assert np.all(np.isfinite(scaling.column_factors))
    assert np.all(scaling.row_factors > 0.0)
    assert np.all(scaling.column_factors > 0.0)
    assert np.linalg.matrix_rank(raw) == np.linalg.matrix_rank(scaled) == 6
    assert np.linalg.slogdet(raw)[0] == np.linalg.slogdet(scaled)[0]
    expected_scaled_determinant = (
        np.prod(scaling.row_factors)
        * np.linalg.det(raw)
        * np.prod(scaling.column_factors)
    )
    assert np.linalg.det(scaled) == pytest.approx(
        expected_scaled_determinant, rel=2.0e-14, abs=0.0
    )

    singular = np.array(raw, copy=True)
    singular[-1, -1] = 0.0
    singular_scaling = iso_audit._positive_equilibrate_model_neutral(singular)
    assert np.linalg.det(singular) == 0.0
    assert np.linalg.det(singular_scaling.scaled_matrix) == 0.0
    assert np.linalg.matrix_rank(singular) == np.linalg.matrix_rank(
        singular_scaling.scaled_matrix
    )


def _postprocess_synthetic_contract(
    *,
    Omega_min: float = 0.2,
    Omega_max: float = 2.8,
    primary_scan_points: int = 53,
    verification_scan_points: int = 105,
    close_pair_subintervals: int = 32,
) -> dict[str, object]:
    return {
        "thresholds": {"root_singular_ratio": 1.0e-9},
        "root_policy": {
            "Omega_min": Omega_min,
            "Omega_max": Omega_max,
            "primary_scan_points": primary_scan_points,
            "verification_scan_points": verification_scan_points,
            "primary_phase": 0.0,
            "verification_phase": 0.5,
            "sigma_prefilter": 1.0e-5,
            "nullity_relative_threshold": 1.0e-12,
            "root_xtol_Omega": 1.0e-12,
            "dedup_atol_Omega": 1.0e-8,
            "dedup_rtol": 1.0e-10,
            "cluster_atol_Omega": 1.0e-8,
            "cluster_rtol": 0.0,
            "post_guard_tail_Omega": 0.2,
            "local_close_pair_guard_subintervals": close_pair_subintervals,
        },
    }


@pytest.mark.parametrize(
    ("sources", "sigma_ratio", "nullity", "expected"),
    [
        (("determinant_sign_bracket",), 1.0e-12, 1, True),
        (("determinant_zero_grid_certified",), 1.0e-12, 1, True),
        (("sigma_ratio_minimum",), 1.0e-12, 1, False),
        (("sigma_ratio_minimum",), 1.0e-12, 2, True),
        (("determinant_sign_bracket",), 1.0e-12, 0, False),
        (("determinant_sign_bracket",), 1.0e-6, 1, False),
        (("determinant_zero_grid",), 1.0e-12, 1, False),
    ],
)
def test_direct_control_candidate_acceptance_requires_certified_simple_crossing(
    sources: tuple[str, ...], sigma_ratio: float, nullity: int, expected: bool
) -> None:
    assert (
        iso_postprocess._direct_candidate_is_canonical_root(
            sources, sigma_ratio, nullity, 1.0e-9
        )
        is expected
    )


def test_direct_control_grid_zero_certification_rejects_numerical_plateau() -> None:
    isolated, isolated_rejections = iso_postprocess._certified_zero_grid_candidates(
        np.array([0.0, 1.0, 2.0]),
        np.array([-1.0, 0.0, 1.0]),
        source_prefix="determinant",
    )
    plateau, plateau_rejections = iso_postprocess._certified_zero_grid_candidates(
        np.array([0.0, 1.0, 2.0, 3.0]),
        np.array([-1.0, 0.0, 0.0, 1.0]),
        source_prefix="determinant",
    )
    assert isolated == [(1.0, "determinant_zero_grid_certified")]
    assert isolated_rejections == 0
    assert plateau == []
    assert plateau_rejections == 1


def test_direct_control_scanner_keeps_even_double_and_same_cell_simple_roots() -> None:
    contract = _postprocess_synthetic_contract(
        primary_scan_points=27,
        verification_scan_points=53,
    )
    provider = _synthetic_spectral_provider(
        lambda value: 1.0e-4 * (value - 1.03),
        lambda value: 1.0e-4 * (value - 1.03),
        lambda value: 1.0e-4 * (value - 1.06),
        lambda value: 1.0e-4 * (value - 2.0),
        1.0,
        1.0,
    )
    audit: dict[str, object] = {}
    events = iso_postprocess._scan_small_matrix(
        provider,
        1.0,
        contract,
        points=27,
        phase=0.0,
        audit=audit,
    )
    assert_allclose(
        [event.Omega for event in events[:3]],
        [1.03, 1.06, 2.0],
        rtol=0.0,
        atol=5.0e-9,
    )
    assert [event.multiplicity for event in events[:3]] == [2, 1, 1]
    assert any(
        "determinant_close_pair_sign_bracket" in event.sources
        for event in events
    )
    assert int(audit["canonical_event_count"]) == len(events)


def test_direct_control_distinct_axial_like_roots_are_not_merged_by_nullvector() -> None:
    contract = _postprocess_synthetic_contract(close_pair_subintervals=0)
    provider = _synthetic_spectral_provider(
        lambda value: (value - 1.0) * (value - 2.0),
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    events = iso_postprocess._scan_small_matrix(
        provider,
        1.0,
        contract,
        points=53,
        phase=0.0,
    )
    assert_allclose(
        [event.Omega for event in events],
        [1.0, 2.0],
        rtol=0.0,
        atol=2.0e-9,
    )


def test_direct_control_post_guard_zero_plateau_does_not_inflate_inventory() -> None:
    roots = tuple(0.4 + 0.1 * index for index in range(13))

    def first_eigenvalue(value: float) -> float:
        if value >= 2.2:
            return 0.0
        return math.prod(value - root for root in roots)

    contract = _postprocess_synthetic_contract(
        Omega_min=0.2,
        Omega_max=2.8,
        primary_scan_points=261,
        verification_scan_points=521,
        close_pair_subintervals=0,
    )
    provider = _synthetic_spectral_provider(
        first_eigenvalue,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
    )
    audit: dict[str, object] = {}
    events = iso_postprocess._scan_small_matrix(
        provider,
        1.0,
        contract,
        points=261,
        phase=0.0,
        audit=audit,
    )
    blocks = iso_postprocess._small_events_to_blocks("SYN-PLATEAU", events, contract)
    assert len(events) == 13
    assert sum(block.root_count for block in blocks) == 13
    assert int(audit["canonical_event_count"]) == 13
    assert max(event.Omega for event in events) < 2.2


def test_closing_classification_preserves_all_auxiliary_numeric_fields() -> None:
    rows: list[dict[str, str]] = []
    for index in range(104):
        rows.append(
            {
                "comparison_kind": iso_postprocess.SCIENTIFIC_COMPARISON_KIND,
                "case_id": f"ISO-{index // 13 + 1:02d}",
                "status": "PASS",
                "reason": "",
                "relative_difference": "1e-12",
                "alignment_slot": str(index % 13 + 1),
            }
        )
    auxiliary_counts = (28, 26, 28, 28, 28, 26)
    for kind, count in zip(
        iso_postprocess.DIRECT_BETA0_AUXILIARY_KINDS,
        auxiliary_counts,
        strict=True,
    ):
        for index in range(count):
            row = {
                "comparison_kind": kind,
                "case_id": "ISO-01",
                "status": "PASS" if index % 2 == 0 else "FAIL",
                "reason": f"old-reason-{index}",
            }
            row.update(
                {
                    field: f"{index + field_index / 1000.0:.17g}"
                    for field_index, field in enumerate(
                        iso_postprocess.DIRECT_NUMERIC_FIELDS
                    )
                }
            )
            rows.append(row)

    before = [dict(row) for row in rows]
    numerical_hash_before = iso_postprocess._direct_numeric_projection_sha256(rows)
    classified, summary = iso_postprocess.classify_existing_spectrum_rows(rows)
    numerical_hash_after = iso_postprocess._direct_numeric_projection_sha256(classified)

    assert summary["scientific_row_count"] == 104
    assert summary["auxiliary_row_count"] == 164
    assert numerical_hash_after == numerical_hash_before
    for original, current in zip(before, classified, strict=True):
        if original["comparison_kind"] == iso_postprocess.SCIENTIFIC_COMPARISON_KIND:
            assert current["status"] == "PASS"
            assert current["used_for_scientific_status"] == "true"
            continue
        assert current["status"] == "SUPERSEDED_DIAGNOSTIC_ARTIFACT"
        assert current["scientific_role"] == "AUXILIARY_DIAGNOSTIC"
        assert current["used_for_scientific_status"] == "false"
        assert current["original_status"] == original["status"]
        assert current["original_reason"] == original["reason"]
        assert current["reason"] == iso_postprocess.DIRECT_BETA0_QUALIFICATION_REASON
        for field in iso_postprocess.DIRECT_NUMERIC_FIELDS:
            assert current[field] == original[field]


def test_bounded_refinement_uses_only_the_supplied_sign_changing_bracket() -> None:
    determinant_calls: list[float] = []
    diagnostics_calls: list[float] = []

    def determinant(value: float) -> float:
        determinant_calls.append(float(value))
        return float(value) - 2.0

    def diagnostics(value: float) -> dict[str, float]:
        diagnostics_calls.append(float(value))
        return {"scaled_sigma_ratio": 0.0, "raw_boundary_null_residual": 0.0}

    result = iso_postprocess.bounded_sign_change_refinement(
        old_Omega=2.1,
        bracket_left_Omega=1.0,
        bracket_right_Omega=3.0,
        frequency_scale=4.0,
        determinant_Omega=determinant,
        diagnostics_Omega=diagnostics,
    )

    assert result["refined_Omega"] == pytest.approx(2.0, rel=0.0, abs=1.0e-15)
    assert result["refined_omega"] == pytest.approx(0.5, rel=0.0, abs=1.0e-15)
    assert determinant_calls
    assert all(1.0 <= value <= 3.0 for value in determinant_calls)
    assert diagnostics_calls == pytest.approx([2.0], rel=0.0, abs=1.0e-15)


def test_closing_entry_point_cannot_call_superseded_or_global_helpers() -> None:
    module_tree = ast.parse(
        (
            ROOT
            / "scripts"
            / "analysis"
            / "laminated_beams"
            / "reddy_four_ply_isotropic_postprocess.py"
        ).read_text(encoding="utf-8")
    )
    function = next(
        node
        for node in module_tree.body
        if isinstance(node, ast.FunctionDef) and node.name == "run_postprocessing"
    )
    called_names = {
        node.func.id
        for node in ast.walk(function)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    assert called_names.isdisjoint(
        {
            "_beta0_direct_rows",
            "_verified_small_spectrum",
            "_scan_small_matrix",
            "_process_selected_modes",
            "seed_free_root_inventory",
        }
    )


def test_saved_symmetry_gate_allows_only_the_declared_legacy_slot12_exception(
    tmp_path: Path,
) -> None:
    rows: list[dict[str, object]] = []
    for method in ("RLB", "LEGACY_RECTANGULAR"):
        for kind in (
            "ANGLE_REFLECTION_PLUS30_MINUS30",
            "ARM_EXCHANGE_0P7_1P3",
        ):
            for slot in range(1, 14):
                declared_exception = (
                    method == "LEGACY_RECTANGULAR"
                    and kind == "ARM_EXCHANGE_0P7_1P3"
                    and slot == 12
                )
                rows.append(
                    {
                        "method": method,
                        "comparison_kind": kind,
                        "alignment_slot": slot,
                        "relative_difference": (
                            1.7901527995729936e-9 if declared_exception else 1.0e-12
                        ),
                        "tolerance": 1.0e-9,
                        "status": "FAIL" if declared_exception else "PASS",
                    }
                )
    iso_postprocess._write_csv(tmp_path / "symmetry_checks.csv", rows)
    summary = iso_postprocess._validate_saved_symmetry_evidence(
        tmp_path, _contract()
    )
    assert summary["row_count"] == 52
    assert summary["pre_refinement_pass_rows"] == 51
    assert summary["sole_pre_refinement_exception"]["alignment_slot"] == 12

    rows[0]["status"] = "FAIL"
    iso_postprocess._write_csv(tmp_path / "symmetry_checks.csv", rows)
    with pytest.raises(RuntimeError, match="unexpected failure"):
        iso_postprocess._validate_saved_symmetry_evidence(tmp_path, _contract())


@pytest.mark.skipif(
    not (RESULT_DIRECTORY / "run_manifest.json").is_file(),
    reason="ignored closing evidence is not present in a clean checkout",
)
def test_existing_closing_evidence_has_exact_frozen_inventories_and_final_status() -> None:
    rlb_manifest = json.loads(
        (RESULT_DIRECTORY / "rlb_inventory_manifest.json").read_text(encoding="utf-8")
    )
    legacy_manifest = json.loads(
        (RESULT_DIRECTORY / "legacy_rectangular_inventory_manifest.json").read_text(
            encoding="utf-8"
        )
    )
    rlb_rows = iso_postprocess._read_csv(RESULT_DIRECTORY / "rlb_roots.csv")
    legacy_rows = iso_postprocess._read_csv(
        RESULT_DIRECTORY / "legacy_rectangular_roots.csv"
    )
    assert len(rlb_rows) == len(legacy_rows) == 104
    assert (
        iso_postprocess._semantic_inventory_hash(rlb_rows, rlb_manifest)
        == iso_postprocess.EXPECTED_RLB_INVENTORY_SHA256
    )
    assert (
        iso_postprocess._semantic_inventory_hash(legacy_rows, legacy_manifest)
        == iso_postprocess.EXPECTED_LEGACY_INVENTORY_SHA256
    )
    for manifest in (rlb_manifest, legacy_manifest):
        assert manifest["seed_free"] is True
        assert manifest["cross_model_roots_read"] is False
        for name, expected in manifest["output_sha256"].items():
            assert iso_postprocess._sha256_file(RESULT_DIRECTORY / name) == expected

    run_manifest = json.loads(
        (RESULT_DIRECTORY / "run_manifest.json").read_text(encoding="utf-8")
    )
    final_status = json.loads(
        (RESULT_DIRECTORY / "final_status.json").read_text(encoding="utf-8")
    )
    assert run_manifest["global_root_searches_run_in_closing_stage"] == 0
    assert run_manifest["ritz_solves_run_in_closing_stage"] == 0
    assert run_manifest["direct_fixed_fixed_solves_run_in_closing_stage"] == 0
    assert run_manifest["scientific_evidence"]["spectrum"]["comparison_count"] == 104
    assert run_manifest["scientific_evidence"]["spectrum"]["status"] == "PASS"
    assert run_manifest["scientific_evidence"]["modes"]["status"] == "PASS"
    assert final_status["scientific_statuses"] == iso_postprocess.SCIENTIFIC_STATUSES
    assert final_status["SCIENTIFIC_OVERALL"] == (
        iso_postprocess.SCIENTIFIC_OVERALL_STATUS
    )
    assert final_status["auxiliary_statuses"]["AUX-LEGACY-ARM-EXCHANGE"] == (
        iso_postprocess.ARM_EXCHANGE_PASS_STATUS
    )


@pytest.mark.skipif(
    not (RESULT_DIRECTORY / "rlb_root_candidates.csv").is_file(),
    reason="ignored frozen candidate evidence is not present in a clean checkout",
)
def test_existing_iso04_close_pair_and_mode_evidence_remain_distinct_and_passing() -> None:
    close_pair = iso_postprocess._verify_iso04_close_pair(RESULT_DIRECTORY)
    assert close_pair["status"] == "PASS"
    for method in ("RLB", "LEGACY_RECTANGULAR"):
        for scan in ("primary", "verification"):
            values = close_pair[method][scan]["Omega"]
            assert len(values) == 2
            assert 0.02 < values[1] - values[0] < 0.04

    mode_evidence = iso_postprocess._validate_mode_evidence(
        RESULT_DIRECTORY, _contract()
    )
    assert mode_evidence["pointwise_comparison_rows"] == 14
    assert mode_evidence["physical_residual_rows"] == 28
    assert mode_evidence["maximum_one_minus_MAC"] <= 1.0e-6
    assert mode_evidence["maximum_joint_compatibility_normalized"] <= 1.0e-10
    assert mode_evidence["maximum_joint_equilibrium_normalized"] <= 1.0e-9
    assert mode_evidence["maximum_energy_identity_relative"] <= 1.0e-8
