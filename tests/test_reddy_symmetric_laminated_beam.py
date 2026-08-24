from __future__ import annotations

import ast
from dataclasses import replace
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.integrate import simpson
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib import reddy_symmetric_laminated_beam as beam  # noqa: E402


E2 = 1.0
E1 = 25.0 * E2
G12 = G13 = 0.5 * E2
G23 = 0.2 * E2
NU12 = 0.25
RHO = 1.0
SOURCE_K = 5.0 / 6.0
SOURCE_JSON = ROOT / "tests" / "data" / "reddy_ch4_table_4_3_3.json"

STACKS = {
    "0°": (0.0,),
    "90°": (90.0,),
    "(0/90)_s": (0.0, 90.0, 90.0, 0.0),
    "(45/-45)_s": (45.0, -45.0, -45.0, 45.0),
}


def _material(**overrides: float) -> beam.OrthotropicLamina:
    values = {
        "E1": E1,
        "E2": E2,
        "nu12": NU12,
        "G12": G12,
        "G13": G13,
        "G23": G23,
        "rho": RHO,
    }
    values.update(overrides)
    return beam.OrthotropicLamina(**values)


def _laminate(
    angles: tuple[float, ...] | list[float],
    *,
    thickness: float = 1.0,
    material: beam.OrthotropicLamina | None = None,
) -> beam.LaminateSection:
    mat = _material() if material is None else material
    ply_thickness = thickness / len(angles)
    return beam.integrate_laminate(
        [beam.Ply(mat, angle, ply_thickness) for angle in angles]
    )


def _properties(
    angles: tuple[float, ...] | list[float] = (0.0, 90.0, 90.0, 0.0),
    *,
    thickness: float = 0.2,
    width: float = 0.7,
    K: float = SOURCE_K,
) -> beam.BeamProperties:
    return beam.reduce_to_beam_properties(
        _laminate(angles, thickness=thickness),
        width=width,
        K=K,
    )


def _mass_weighted_bending_mac(
    x: np.ndarray,
    first: np.ndarray,
    second: np.ndarray,
    properties: beam.BeamProperties,
) -> float:
    inner = simpson(
        properties.m * first[:, 0] * second[:, 0]
        + properties.J * first[:, 1] * second[:, 1],
        x=x,
    )
    first_mass = beam.bending_modal_mass(x, first, properties)
    second_mass = beam.bending_modal_mass(x, second, properties)
    return float(abs(inner) ** 2 / (first_mass * second_mass))


def _amplitude_ratio(x: np.ndarray, states: np.ndarray, k: float) -> float:
    sine = np.sin(k * x)
    cosine = np.cos(k * x)
    w_amplitude = float(sine @ states[:, 0] / (sine @ sine))
    psi_amplitude = float(cosine @ states[:, 1] / (cosine @ cosine))
    return psi_amplitude / w_amplitude


def _source_frequency_scale(laminate: beam.LaminateSection, length: float) -> float:
    h = laminate.thickness
    return length**2 * math.sqrt(laminate.I0 / (E2 * h**3))


def _direct_transfer_root_near(
    source_omega_bar: float,
    scale: float,
    length: float,
    properties: beam.BeamProperties,
    boundary_condition: str,
) -> float:
    center = source_omega_bar / scale

    def determinant(omega: float) -> float:
        return beam.bending_characteristic_determinant(
            omega, length, properties, boundary_condition
        )

    def verify(root: float) -> float:
        diagnostic = beam.bending_root_diagnostic(
            root,
            length,
            properties,
            boundary_condition,
            sigma_ratio_tolerance=1.0e-10,
            detection="source_bracket",
        )
        assert diagnostic.accepted
        assert diagnostic.sigma_ratio <= 1.0e-10
        return root

    for relative_half_width in (0.005, 0.02, 0.08, 0.20):
        lower = center * (1.0 - relative_half_width)
        upper = center * (1.0 + relative_half_width)
        f_lower = determinant(lower)
        f_upper = determinant(upper)
        if f_lower == 0.0:
            return verify(lower)
        if f_upper == 0.0:
            return verify(upper)
        if math.isfinite(f_lower) and math.isfinite(f_upper) and f_lower * f_upper < 0.0:
            root = float(
                brentq(
                    determinant,
                    lower,
                    upper,
                    xtol=1.0e-14,
                    rtol=4.0 * np.finfo(float).eps,
                )
            )
            return verify(root)
    raise AssertionError(
        f"No direct-transfer sign-changing root near omega_bar={source_omega_bar}, "
        f"BC={boundary_condition}."
    )


def _source_classical_characteristic_frequency(
    properties: beam.BeamProperties,
    length: float,
    boundary_condition: str,
    source_roots: dict[str, object],
) -> float:
    # These are the source's printed classical characteristic constants.
    alpha = float(source_roots[boundary_condition]["value"])
    k = alpha / length
    S, D, mass, rotary = properties.S, properties.D, properties.m, properties.J
    if rotary == 0.0:
        omega_squared = S * D * k**4 / (mass * (S + D * k**2))
    else:
        coefficient_a = mass * rotary
        coefficient_b = rotary * S * k**2 + mass * (S + D * k**2)
        coefficient_c = S * D * k**4
        discriminant = coefficient_b**2 - 4.0 * coefficient_a * coefficient_c
        large = 0.5 * (coefficient_b + math.sqrt(discriminant)) / coefficient_a
        omega_squared = coefficient_c / (coefficient_a * large)
    return math.sqrt(omega_squared)


def test_lamina_q_reciprocity_and_input_validation() -> None:
    material = _material()
    denominator = 1.0 - material.nu12 * material.nu21
    expected = np.array(
        [
            [material.E1 / denominator, material.nu12 * material.E2 / denominator, 0.0],
            [material.nu12 * material.E2 / denominator, material.E2 / denominator, 0.0],
            [0.0, 0.0, material.G12],
        ]
    )
    assert material.nu21 == pytest.approx(material.nu12 * material.E2 / material.E1)
    assert material.nu12 / material.E1 == pytest.approx(material.nu21 / material.E2)
    assert_allclose(beam.lamina_Q(material), expected, rtol=2.0e-15, atol=0.0)
    assert_allclose(beam.lamina_Q(material), beam.lamina_Q(material).T, atol=0.0)
    assert np.all(np.linalg.eigvalsh(beam.lamina_Q(material)) > 0.0)

    with pytest.raises(ValueError):
        _material(E1=0.0)
    with pytest.raises(ValueError):
        _material(rho=-1.0)
    with pytest.raises(ValueError):
        _material(nu12=6.0)
    with pytest.raises(ValueError):
        beam.Ply(material, 0.0, 0.0)
    with pytest.raises(ValueError):
        beam.Ply(material, math.nan, 1.0)
    with pytest.raises(ValueError):
        beam.integrate_laminate([])
    with pytest.raises(ValueError):
        beam.transformed_reduced_stiffness(material, math.nan)
    with pytest.raises(ValueError):
        beam.transformed_transverse_shear_stiffness(material, math.inf)


def test_qbar_zero_ninety_and_angle_reversal_symmetry() -> None:
    material = _material()
    Q = beam.lamina_Q(material)
    at_zero = beam.lamina_Qbar(material, 0.0)
    at_ninety = beam.lamina_Qbar(material, 90.0)
    assert_allclose(at_zero, Q, rtol=0.0, atol=2.0e-15)
    assert_allclose(
        at_ninety,
        np.array(
            [
                [Q[1, 1], Q[0, 1], 0.0],
                [Q[0, 1], Q[0, 0], 0.0],
                [0.0, 0.0, Q[2, 2]],
            ]
        ),
        rtol=2.0e-15,
        atol=2.0e-14,
    )
    assert_allclose(
        beam.lamina_shear_Qbar(material, 0.0),
        np.diag([material.G23, material.G13]),
        atol=0.0,
    )
    assert_allclose(
        beam.lamina_shear_Qbar(material, 90.0),
        np.diag([material.G13, material.G23]),
        rtol=0.0,
        atol=3.0e-17,
    )

    positive = beam.lamina_Qbar(material, 37.0)
    negative = beam.lamina_Qbar(material, -37.0)
    assert_allclose(positive[np.ix_([0, 1, 2], [0, 1, 2])], positive.T)
    assert_allclose(positive[[0, 1, 2], [0, 1, 2]], negative[[0, 1, 2], [0, 1, 2]])
    assert positive[0, 1] == pytest.approx(negative[0, 1])
    assert positive[0, 2] == pytest.approx(-negative[0, 2])
    assert positive[1, 2] == pytest.approx(-negative[1, 2])
    shear_positive = beam.lamina_shear_Qbar(material, 37.0)
    shear_negative = beam.lamina_shear_Qbar(material, -37.0)
    assert_allclose(np.diag(shear_positive), np.diag(shear_negative))
    assert shear_positive[0, 1] == pytest.approx(-shear_negative[0, 1])
    for angle in (-90.0, -45.0, 0.0, 27.0, 90.0):
        assert np.all(np.linalg.eigvalsh(beam.lamina_Qbar(material, angle)) > 0.0)
        assert np.all(
            np.linalg.eigvalsh(beam.lamina_shear_Qbar(material, angle)) > 0.0
        )


@pytest.mark.parametrize("used_label", tuple(STACKS))
def test_laminate_thickness_symmetry_mass_and_positive_matrices(
    used_label: str,
) -> None:
    thickness = 0.36
    laminate = _laminate(STACKS[used_label], thickness=thickness)
    assert laminate.thickness == pytest.approx(thickness, abs=2.0e-16)
    assert laminate.z_interfaces[0] == pytest.approx(-thickness / 2.0)
    assert laminate.z_interfaces[-1] == pytest.approx(thickness / 2.0)
    assert np.all(np.diff(laminate.z_interfaces) > 0.0)
    symmetry = beam.check_laminate_symmetry(laminate, tolerance=1.0e-12)
    assert symmetry.is_symmetric
    assert symmetry.B_relative <= 1.0e-14
    assert symmetry.I1_relative <= 1.0e-14
    assert laminate.I0 == pytest.approx(RHO * thickness)
    assert laminate.I1 == pytest.approx(0.0, abs=2.0e-18)
    assert laminate.I2 == pytest.approx(RHO * thickness**3 / 12.0, rel=2.0e-15)
    for matrix in (laminate.A, laminate.D, laminate.shear):
        assert_allclose(matrix, matrix.T, rtol=0.0, atol=2.0e-15)
        assert np.min(np.linalg.eigvalsh(matrix)) > 0.0


def test_compliance_schur_single_ply_and_invalid_reductions() -> None:
    matrix = np.array([[8.0, 1.5, -0.2], [1.5, 5.0, 0.4], [-0.2, 0.4, 3.0]])
    for index in range(3):
        reduction = beam.effective_stiffness_reduction(matrix, index)
        assert reduction.compliance_value == pytest.approx(reduction.schur_value, rel=2.0e-15)
        assert reduction.relative_difference <= 2.0e-15

    material = _material()
    h, width, K = 0.4, 0.73, 0.79
    laminate = _laminate((0.0,), thickness=h, material=material)
    properties = beam.reduce_to_beam_properties(
        laminate, width=width, K=K, reduction_tolerance=1.0e-13
    )
    assert properties.A == pytest.approx(width * material.E1 * h, rel=2.0e-15)
    assert properties.D == pytest.approx(
        width * material.E1 * h**3 / 12.0, rel=2.0e-15
    )
    assert properties.S == pytest.approx(K * width * material.G13 * h, rel=2.0e-15)
    assert properties.m == pytest.approx(width * material.rho * h, rel=2.0e-15)
    assert properties.J == pytest.approx(
        width * material.rho * h**3 / 12.0, rel=2.0e-15
    )
    for reduction in (
        properties.axial_reduction,
        properties.bending_reduction,
        properties.shear_reduction_before_K,
    ):
        assert reduction is not None
        assert reduction.relative_difference <= 2.0e-15

    with pytest.raises(ValueError):
        beam.effective_stiffness_reduction(np.ones((2, 3)), 0)
    with pytest.raises(ValueError):
        beam.effective_stiffness_reduction(np.array([[1.0, 2.0], [0.0, 1.0]]), 0)
    with pytest.raises(IndexError):
        beam.effective_stiffness_reduction(np.eye(2), 2)
    with pytest.raises(ValueError):
        beam.reduce_to_beam_properties(laminate, width=0.0, K=K)
    with pytest.raises(ValueError):
        beam.reduce_to_beam_properties(laminate, width=width, K=0.0)
    nonsymmetric = _laminate((0.0, 90.0), thickness=h)
    with pytest.raises(ValueError, match="Only symmetric laminates"):
        beam.reduce_to_beam_properties(nonsymmetric, width=width, K=K)


@pytest.mark.parametrize("angle", (0.0, 90.0), ids=("0deg", "90deg"))
def test_width_scaling_partition_and_frequency_invariance(angle: float) -> None:
    h = 0.4
    one_ply = _laminate((angle,), thickness=h)
    eight_plies = _laminate((angle,) * 8, thickness=h)
    for name in ("A", "B", "D", "shear"):
        assert_allclose(
            getattr(one_ply, name),
            getattr(eight_plies, name),
            rtol=3.0e-15,
            atol=3.0e-16,
        )
    for name in ("I0", "I1", "I2"):
        assert getattr(one_ply, name) == pytest.approx(
            getattr(eight_plies, name), rel=3.0e-15, abs=3.0e-18
        )

    narrow = beam.reduce_to_beam_properties(one_ply, width=0.4, K=SOURCE_K)
    wide = beam.reduce_to_beam_properties(one_ply, width=1.8, K=SOURCE_K)
    partitioned = beam.reduce_to_beam_properties(eight_plies, width=0.4, K=SOURCE_K)
    factor = wide.width / narrow.width
    for name in ("A", "D", "S", "m", "J"):
        assert getattr(wide, name) == pytest.approx(factor * getattr(narrow, name))
        assert getattr(partitioned, name) == pytest.approx(
            getattr(narrow, name), rel=4.0e-15
        )
    narrow_modes = beam.hinged_hinged_modes(narrow, 1.2, 3)
    wide_modes = beam.hinged_hinged_modes(wide, 1.2, 3)
    partitioned_modes = beam.hinged_hinged_modes(partitioned, 1.2, 3)
    assert_allclose(
        [mode.omega for mode in narrow_modes],
        [mode.omega for mode in wide_modes],
        rtol=3.0e-14,
        atol=0.0,
    )
    assert_allclose(
        [mode.omega for mode in narrow_modes],
        [mode.omega for mode in partitioned_modes],
        rtol=3.0e-14,
        atol=0.0,
    )


def test_bending_state_coefficients_static_transfer_semigroup_and_rank() -> None:
    properties = _properties()
    omega = 1.37
    expected = np.array(
        [
            [0.0, -1.0, 1.0 / properties.S, 0.0],
            [0.0, 0.0, 0.0, 1.0 / properties.D],
            [-properties.m * omega**2, 0.0, 0.0, 0.0],
            [0.0, -properties.J * omega**2, 1.0, 0.0],
        ]
    )
    assert_allclose(beam.bending_state_matrix(omega, properties), expected, atol=0.0)
    assert_allclose(beam.bending_transfer_matrix(omega, 0.0, properties), np.eye(4))

    x = 0.83
    static_expected = np.array(
        [
            [
                1.0,
                -x,
                x / properties.S - x**3 / (6.0 * properties.D),
                -x**2 / (2.0 * properties.D),
            ],
            [0.0, 1.0, x**2 / (2.0 * properties.D), x / properties.D],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, x, 1.0],
        ]
    )
    assert_allclose(
        beam.bending_transfer_matrix(0.0, x, properties),
        static_expected,
        rtol=3.0e-14,
        atol=3.0e-14,
    )
    left_length, right_length = 0.37, 0.52
    assert_allclose(
        beam.bending_transfer_matrix(omega, left_length + right_length, properties),
        beam.bending_transfer_matrix(omega, right_length, properties)
        @ beam.bending_transfer_matrix(omega, left_length, properties),
        rtol=4.0e-14,
        atol=4.0e-14,
    )

    mode = beam.hinged_hinged_algebraic_modes(properties, 1.0, 1)[0]
    reduced = beam.bending_boundary_matrix(mode.omega, 1.0, properties, "HH")
    full = beam.bending_full_boundary_matrix(mode.omega, 1.0, properties, "HH")
    assert np.linalg.matrix_rank(reduced, tol=1.0e-10) == 1
    assert np.linalg.matrix_rank(full, tol=1.0e-10) == 3
    generic = beam.bending_full_boundary_matrix(
        0.73 * mode.omega, 1.0, properties, "HH"
    )
    assert np.linalg.matrix_rank(generic, tol=1.0e-10) == 4


def test_hh_n1_to_n3_both_branches_transfer_ratio_mac_mass_and_energy() -> None:
    properties = _properties()
    length = 1.0
    x = np.linspace(0.0, length, 801)
    modes = beam.hinged_hinged_modes(properties, length, 3)
    assert [(mode.n, mode.branch) for mode in modes] == [
        (1, "lower"),
        (1, "upper"),
        (2, "lower"),
        (2, "upper"),
        (3, "lower"),
        (3, "upper"),
    ]

    for mode in modes:
        k = mode.wavenumber
        expected_ratio = (
            properties.m * mode.omega_squared - properties.S * k**2
        ) / (properties.S * k)
        assert mode.psi_over_w == pytest.approx(expected_ratio, rel=2.0e-14)
        second_equation_residual = (
            -properties.D * k**2 * mode.psi_over_w
            - properties.S * (k + mode.psi_over_w)
            + properties.J * mode.omega_squared * mode.psi_over_w
        )
        scale = max(
            abs(properties.D * k**2 * mode.psi_over_w),
            abs(properties.S * (k + mode.psi_over_w)),
            abs(properties.J * mode.omega_squared * mode.psi_over_w),
            1.0,
        )
        assert abs(second_equation_residual) / scale <= 2.0e-13

        diagnostic = beam.bending_root_diagnostic(
            mode.omega, length, properties, "HH", sigma_ratio_tolerance=1.0e-10
        )
        assert diagnostic.accepted
        assert diagnostic.sigma_ratio <= 1.0e-10

        exact = beam.hinged_hinged_exact_shape(mode, properties, length, x)
        transferred = beam.bending_mode_shape(
            mode.omega, properties, length, "HH", x
        )
        for shape in (exact, transferred):
            threshold = 64.0 * np.finfo(float).eps * np.max(
                np.abs(shape.states[:, 0])
            )
            first_significant = np.flatnonzero(
                np.abs(shape.states[:, 0]) > threshold
            )[0]
            assert shape.states[first_significant, 0] > 0.0
        assert exact.boundary_residual <= 1.0e-11
        assert transferred.boundary_residual <= 1.0e-10
        assert beam.bending_modal_mass(x, exact.states, properties) == pytest.approx(
            1.0, rel=3.0e-13
        )
        assert beam.bending_modal_mass(
            x, transferred.states, properties
        ) == pytest.approx(1.0, rel=3.0e-13)
        assert _amplitude_ratio(x, exact.states, k) == pytest.approx(
            mode.psi_over_w, rel=3.0e-13, abs=3.0e-13
        )
        assert _amplitude_ratio(x, transferred.states, k) == pytest.approx(
            mode.psi_over_w, rel=3.0e-11, abs=3.0e-11
        )
        assert _mass_weighted_bending_mac(
            x, exact.states, transferred.states, properties
        ) >= 1.0 - 1.0e-10
        assert exact.energies.identity_relative_error <= 1.0e-10
        assert transferred.energies.identity_relative_error <= 1.0e-10


def test_source_contract_identity_roles_and_missing_rows() -> None:
    data = json.loads(SOURCE_JSON.read_text(encoding="utf-8"))
    assert data["source_status"] == "PASS_WITH_DOCUMENTED_SOURCE_RECONSTRUCTION"
    policy = data["approved_source_policy"]
    assert policy["cross_ply"] == {
        "printed_label": "(90/0)_s",
        "used_label": "(0/90)_s",
        "used_stack_degrees": [0, 90, 90, 0],
        "correction_status": "CORRECTED_BY_INTERNAL_SOURCE_CROSSCHECK",
        "evidence": ["Table 4.2.4", "Figures 4.3.3–4.3.4"],
    }
    assert policy["shear_correction_factor"]["expression"] == "5/6"
    assert policy["shear_correction_factor"]["value"] == pytest.approx(SOURCE_K)
    assert (
        policy["shear_correction_factor"]["provenance"]
        == "INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY"
    )
    assert (
        policy["equation_4_3_50a"]["status"]
        == "PRINTED_INTERMEDIATE_FORMULA_INCONSISTENCY"
    )
    material = data["material_definition"]
    assert material["source_formula"] == "Eq. (4.2.25)"
    assert material["E1_over_E2"] == E1 / E2
    assert material["G12_over_E2"] == G12 / E2
    assert material["G13_over_E2"] == G13 / E2
    assert material["G23_over_E2"] == G23 / E2
    assert material["nu12"] == NU12
    assert material["dimensionless_realization"] == {
        "E2": E2,
        "rho": RHO,
        "thickness": 1.0,
        "width": 1.0,
    }
    normalization = data["normalization"]
    assert normalization["source_formula"] == "Eq. (4.3.52)"
    assert (
        normalization["angular_frequency"]
        == "omega_bar = omega_1*length^2*sqrt(I0/(E2*thickness^3))"
    )
    assert normalization["I0_definition"] == "I0 = integral(rho dz)"
    ordering = data["coordinate_and_ordering"]
    assert ordering["membrane_and_bending_components"] == ["xx", "yy", "xy"]
    assert ordering["transverse_shear_components"] == ["yz", "xz"]
    assert ordering["coordinates"]["z_reddy"] == (
        "thickness coordinate, positive downward"
    )
    assert ordering["reddy_thickness_coordinate"]["top_surface"] == (
        "z_reddy = -thickness/2"
    )
    assert ordering["reddy_thickness_coordinate"]["bottom_surface"] == (
        "z_reddy = +thickness/2"
    )
    assert ordering["project_storage_convention"]["coordinate"] == (
        "z_project = -z_reddy"
    )
    assert ordering["ply_stack_direction"] == (
        "project storage: bottom-to-top in increasing z_project"
    )
    assert ordering["bending_state"] == ["w", "psi_b", "Q", "M"]
    assert ordering["axial_state"] == ["u", "N"]
    assert ordering["combined_state"] == ["u", "w", "psi_b", "N", "Q", "M"]
    source_roots = data["source_classical_characteristic_roots"]
    assert source_roots["source_item"] == "Table 4.2.3"
    assert source_roots["printed_page"] == 184
    assert source_roots["pdf_index_zero_based"] == 206
    assert source_roots["HH"] == {"expression": "pi", "value": math.pi}
    assert source_roots["CC"] == {"expression": "printed 4.730", "value": 4.73}
    assert source_roots["CF"] == {"expression": "printed 1.875", "value": 1.875}
    assert data["laminate_definitions"]["0_deg"]["printed_label"] == "0"
    assert data["laminate_definitions"]["0_deg"]["used_label"] == "0°"
    assert data["laminate_definitions"]["90_deg"]["printed_label"] == "90"
    assert data["laminate_definitions"]["90_deg"]["used_label"] == "90°"

    records = data["records"]
    keys = {
        (
            row["laminate_id"],
            row["row_role"],
            row["boundary_condition"],
            row["a_over_h"],
        )
        for row in records
    }
    assert len(records) == len(keys) == 180
    published = [row for row in records if row["source_status"] == "PUBLISHED"]
    missing = [
        row for row in records if row["source_status"] == "NOT_PUBLISHED_IN_TABLE"
    ]
    assert len(published) == 144
    assert all(
        row["printed_label"] == "0"
        for row in records
        if row["laminate_id"] == "0_deg"
    )
    assert all(
        row["printed_label"] == "90"
        for row in records
        if row["laminate_id"] == "90_deg"
    )
    assert len(missing) == 36
    assert all(row["source_value"] is None for row in missing)
    assert all(not row["include_in_source_pass_fail"] for row in missing)
    assert all(row["printed_tolerance"] is None for row in missing)
    assert {
        row["used_label"] for row in missing
    } == {"(0/90)_s", "(45/-45)_s"}
    assert {
        row["row_role"] for row in missing
    } == {
        "fsdt_frequency_without_RI",
        "source_classical_characteristic_without_RI",
    }
    included = [row for row in records if row["include_in_source_pass_fail"]]
    assert len(included) == 108
    assert all(row["quantity"] == "omega_bar" for row in included)
    assert all(row["printed_tolerance"] == 0.0005 for row in included)
    assert all(
        row["row_role"]
        in {
            "fsdt_frequency_with_RI",
            "source_classical_characteristic_with_RI",
            "fsdt_frequency_without_RI",
            "source_classical_characteristic_without_RI",
        }
        for row in included
    )


def test_table_4_3_3_partial_pass_mismatch_set_is_explicit_and_stable() -> None:
    data = json.loads(SOURCE_JSON.read_text(encoding="utf-8"))
    included = [row for row in data["records"] if row["include_in_source_pass_fail"]]
    laminate_cache: dict[str, beam.LaminateSection] = {}
    properties_cache: dict[tuple[str, bool], beam.BeamProperties] = {}
    direct_root_cache: dict[tuple[str, bool, str, int], float] = {}
    comparisons: list[dict[str, object]] = []

    for row in included:
        used_label = row["used_label"]
        with_rotary_inertia = not row["row_role"].endswith("without_RI")
        if used_label not in laminate_cache:
            laminate_cache[used_label] = _laminate(STACKS[used_label], thickness=1.0)
        laminate = laminate_cache[used_label]
        property_key = (used_label, with_rotary_inertia)
        if property_key not in properties_cache:
            properties = beam.reduce_to_beam_properties(
                laminate, width=1.0, K=SOURCE_K
            )
            properties_cache[property_key] = (
                properties
                if with_rotary_inertia
                else beam.without_rotary_inertia(properties)
            )
        properties = properties_cache[property_key]
        length = float(row["a_over_h"])
        scale = _source_frequency_scale(laminate, length)
        role = row["row_role"]
        if role.startswith("fsdt_frequency"):
            root_key = (
                used_label,
                with_rotary_inertia,
                row["boundary_condition"],
                row["a_over_h"],
            )
            if root_key not in direct_root_cache:
                direct_root_cache[root_key] = _direct_transfer_root_near(
                    float(row["source_value"]),
                    scale,
                    length,
                    properties,
                    row["boundary_condition"],
                )
            omega = direct_root_cache[root_key]
        elif role in {
            "source_classical_characteristic_with_RI",
            "source_classical_characteristic_without_RI",
        }:
            omega = _source_classical_characteristic_frequency(
                properties,
                length,
                row["boundary_condition"],
                data["source_classical_characteristic_roots"],
            )
        else:
            raise AssertionError(f"Unexpected included source role: {role}")
        computed = omega * scale
        error = abs(computed - float(row["source_value"]))
        comparisons.append(
            {
                "key": (
                    used_label,
                    role,
                    row["boundary_condition"],
                    row["a_over_h"],
                ),
                "tier": row["benchmark_tier"],
                "source": float(row["source_value"]),
                "computed": computed,
                "error": error,
                "passes": error <= float(row["printed_tolerance"]) + 2.0e-12,
            }
        )

    expected_mismatches = {
        ("0°", "fsdt_frequency_with_RI", "CF", 100),
        ("0°", "fsdt_frequency_with_RI", "CF", 20),
        ("0°", "fsdt_frequency_with_RI", "CF", 10),
        ("0°", "fsdt_frequency_without_RI", "CC", 100),
        ("0°", "fsdt_frequency_without_RI", "CC", 20),
        ("0°", "fsdt_frequency_without_RI", "CC", 10),
        ("0°", "fsdt_frequency_without_RI", "CF", 100),
        ("0°", "fsdt_frequency_without_RI", "CF", 20),
        ("0°", "fsdt_frequency_without_RI", "CF", 10),
        ("90°", "fsdt_frequency_with_RI", "CF", 10),
        ("90°", "fsdt_frequency_without_RI", "CC", 100),
        ("90°", "fsdt_frequency_without_RI", "CC", 20),
        ("90°", "fsdt_frequency_without_RI", "CC", 10),
        ("90°", "fsdt_frequency_without_RI", "CF", 20),
        ("90°", "fsdt_frequency_without_RI", "CF", 10),
        ("(0/90)_s", "fsdt_frequency_with_RI", "CF", 20),
        ("(0/90)_s", "fsdt_frequency_with_RI", "CF", 10),
        ("(45/-45)_s", "fsdt_frequency_with_RI", "CC", 20),
        ("(45/-45)_s", "fsdt_frequency_with_RI", "CC", 10),
        ("(45/-45)_s", "fsdt_frequency_with_RI", "CF", 20),
        ("(45/-45)_s", "fsdt_frequency_with_RI", "CF", 10),
        (
            "(45/-45)_s",
            "source_classical_characteristic_with_RI",
            "CC",
            100,
        ),
    }
    observed_mismatches = {
        comparison["key"]
        for comparison in comparisons
        if not comparison["passes"]
    }
    assert len(comparisons) == 108
    assert observed_mismatches == expected_mismatches
    mismatches = [item for item in comparisons if not item["passes"]]
    assert len(mismatches) == 22
    assert sum(item["tier"] == "A" for item in mismatches) == 15
    assert sum(item["tier"] == "B" for item in mismatches) == 7
    assert sum(item["passes"] for item in comparisons) == 86

    approved_hh_ri = [
        item
        for item in comparisons
        if item["key"][1] == "fsdt_frequency_with_RI"
        and item["key"][2] == "HH"
    ]
    assert len(approved_hh_ri) == 12
    assert all(item["passes"] for item in approved_hh_ri)

    by_key = {item["key"]: item for item in comparisons}
    assert by_key[
        ("0°", "fsdt_frequency_without_RI", "CC", 20)
    ]["computed"] == pytest.approx(25.3394906128, rel=2.0e-11)
    assert by_key[
        ("0°", "fsdt_frequency_with_RI", "CF", 10)
    ]["computed"] == pytest.approx(4.55981341042, rel=2.0e-11)
    assert by_key[
        ("(0/90)_s", "fsdt_frequency_with_RI", "CF", 10)
    ]["computed"] == pytest.approx(4.17837401438, rel=2.0e-11)
    assert by_key[
        ("(45/-45)_s", "source_classical_characteristic_with_RI", "CC", 100)
    ]["computed"] == pytest.approx(8.53162774515, rel=2.0e-11)
    assert len(direct_root_cache) == 54


def test_eq_4_3_51_is_secondary_check_of_independent_cc_transfer_root() -> None:
    laminate = _laminate((0.0,), thickness=1.0)
    properties = beam.reduce_to_beam_properties(
        laminate, width=1.0, K=SOURCE_K
    )
    length = 20.0
    scale = _source_frequency_scale(laminate, length)
    root = _direct_transfer_root_near(
        25.327, scale, length, properties, "CC"
    )
    transfer_diagnostic = beam.bending_root_diagnostic(
        root, length, properties, "CC", sigma_ratio_tolerance=1.0e-10
    )
    source_characteristic = beam.reddy_eq_4_3_51_characteristic(
        root, length, properties
    )
    assert transfer_diagnostic.sigma_ratio <= 1.0e-10
    assert abs(source_characteristic.imag) <= 1.0e-12
    assert abs(source_characteristic.real) <= 1.0e-9


@pytest.mark.parametrize("boundary_condition", ["FF", "FC"])
def test_first_six_axial_modes_transfer_shapes_mass_energy_and_independence(
    boundary_condition: str,
) -> None:
    properties = _properties()
    length = 1.3
    x = np.linspace(0.0, length, 801)
    modes = beam.exact_axial_modes(properties, length, boundary_condition, 6)
    altered_bending = replace(
        properties,
        D=3.7 * properties.D,
        S=0.41 * properties.S,
        J=2.3 * properties.J,
    )
    altered_modes = beam.exact_axial_modes(
        altered_bending, length, boundary_condition, 6
    )
    trial_omega = 2.17
    assert_allclose(
        beam.axial_state_matrix(trial_omega, properties),
        np.array(
            [
                [0.0, 1.0 / properties.A],
                [-properties.m * trial_omega**2, 0.0],
            ]
        ),
        rtol=0.0,
        atol=0.0,
    )
    assert_allclose(
        beam.axial_transfer_matrix(trial_omega, 0.0, properties),
        np.eye(2),
        rtol=0.0,
        atol=0.0,
    )
    assert_allclose(
        [mode.omega for mode in modes],
        [mode.omega for mode in altered_modes],
        rtol=0.0,
        atol=0.0,
    )

    for mode in modes:
        reduced_boundary = beam.axial_boundary_matrix(
            mode.omega, length, properties, boundary_condition
        )
        assert abs(float(reduced_boundary[0, 0])) <= 2.0e-13
        assert_allclose(
            beam.axial_transfer_matrix(
                mode.omega, length, altered_bending
            ),
            beam.axial_transfer_matrix(mode.omega, length, properties),
            rtol=0.0,
            atol=0.0,
        )
        shape = beam.axial_exact_shape(mode, properties, length, x)
        assert shape.boundary_residual <= 2.0e-12
        assert beam.axial_modal_mass(x, shape.states, properties) == pytest.approx(
            1.0, rel=4.0e-13
        )
        assert shape.energies.identity_relative_error <= 1.0e-11

        initial = shape.states[0]
        sample_indices = (0, 137, 400, 800)
        transferred = np.vstack(
            [
                beam.axial_transfer_matrix(
                    mode.omega, float(x[index]), properties
                )
                @ initial
                for index in sample_indices
            ]
        )
        assert_allclose(
            transferred,
            shape.states[list(sample_indices)],
            rtol=3.0e-12,
            atol=3.0e-12,
        )


def test_combined_block_identity_union_and_nondegenerate_purity() -> None:
    properties = _properties()
    length = 1.0
    omega = 1.23
    permutation = [0, 3, 1, 2, 4, 5]
    expected_state = np.zeros((6, 6))
    expected_state[:2, :2] = beam.axial_state_matrix(omega, properties)
    expected_state[2:, 2:] = beam.bending_state_matrix(omega, properties)
    combined_state = beam.combined_state_matrix(omega, properties)
    assert_allclose(
        combined_state[np.ix_(permutation, permutation)],
        expected_state,
        rtol=0.0,
        atol=0.0,
    )
    expected_transfer = np.zeros((6, 6))
    expected_transfer[:2, :2] = beam.axial_transfer_matrix(
        omega, length, properties
    )
    expected_transfer[2:, 2:] = beam.bending_transfer_matrix(
        omega, length, properties
    )
    combined_transfer = beam.combined_transfer_matrix(omega, length, properties)
    assert_allclose(
        combined_transfer[np.ix_(permutation, permutation)],
        expected_transfer,
        rtol=5.0e-14,
        atol=5.0e-14,
    )

    axial = beam.exact_axial_modes(properties, length, "FF", 3)
    bending = tuple(
        mode for mode in beam.hinged_hinged_modes(properties, length, 3)
        if mode.branch == "lower"
    )
    clusters = beam.combined_spectrum_union(
        [mode.omega for mode in axial],
        [mode.omega for mode in bending],
        atol=1.0e-13,
        rtol=1.0e-12,
    )
    flattened = [
        member.omega for cluster in clusters for member in cluster.members
    ]
    assert_allclose(
        flattened,
        sorted([mode.omega for mode in axial] + [mode.omega for mode in bending]),
        rtol=0.0,
        atol=0.0,
    )
    assert sum(cluster.multiplicity for cluster in clusters) == 6

    axial_root_matrix = beam.combined_boundary_matrix(
        axial[0].omega, length, properties, "FF", "HH"
    )
    _, singular, vh = np.linalg.svd(axial_root_matrix)
    assert singular[-1] / singular[0] <= 1.0e-12
    axial_null = vh[-1]
    assert np.linalg.norm(axial_null[1:]) <= 1.0e-10 * abs(axial_null[0])

    bending_root_matrix = beam.combined_boundary_matrix(
        bending[0].omega, length, properties, "FF", "HH"
    )
    _, singular, vh = np.linalg.svd(bending_root_matrix)
    assert singular[-1] / singular[0] <= 1.0e-12
    bending_null = vh[-1]
    assert abs(bending_null[0]) <= 1.0e-10 * np.linalg.norm(bending_null[1:])

    near = beam.combined_spectrum_union(
        [1.0], [1.0 + 1.0e-9], atol=0.0, rtol=2.0e-9
    )
    assert len(near) == 1
    assert near[0].multiplicity == 2
    assert not near[0].exact_degeneracy
    assert {member.subsystem for member in near[0].members} == {"axial", "bending"}


@pytest.fixture(
    scope="module",
    params=(
        pytest.param(
            ((0.0, 90.0, 90.0, 0.0), "FF", "CC"),
            id="cross-ply-FF-CC",
        ),
        pytest.param(
            ((45.0, -45.0, -45.0, 45.0), "FC", "CF"),
            id="angle-ply-FC-CF",
        ),
    ),
)
def required_combined_case(request: pytest.FixtureRequest) -> tuple[
    beam.BeamProperties, str, str, tuple[beam.RootDiagnostic, ...]
]:
    angles, axial_boundary_condition, bending_boundary_condition = request.param
    properties = _properties(angles)
    roots = beam.find_bending_roots(
        properties,
        1.0,
        bending_boundary_condition,
        omega_max=100.0,
        n_roots=6,
        scan_points=2001,
        sigma_ratio_tolerance=1.0e-9,
    )
    return properties, axial_boundary_condition, bending_boundary_condition, roots


def test_required_physical_combined_cases_are_exact_subsystem_unions(
    required_combined_case: tuple[
        beam.BeamProperties, str, str, tuple[beam.RootDiagnostic, ...]
    ],
) -> None:
    properties, axial_bc, bending_bc, bending_roots = required_combined_case
    length = 1.0
    left, right, admissible_basis = beam.combined_boundary_operators(
        axial_bc, bending_bc
    )
    expected_left = np.eye(6)[[0, 1, 2]]
    expected_right = (
        np.eye(6)[[0, 1, 2]] if (axial_bc, bending_bc) == ("FF", "CC")
        else np.eye(6)[[3, 4, 5]]
    )
    assert_allclose(left, expected_left, rtol=0.0, atol=0.0)
    assert_allclose(right, expected_right, rtol=0.0, atol=0.0)
    assert_allclose(left @ admissible_basis, np.zeros((3, 3)), atol=0.0)
    assert np.linalg.matrix_rank(admissible_basis) == 3

    arbitrary_omega = 1.234
    combined = beam.combined_boundary_matrix(
        arbitrary_omega, length, properties, axial_bc, bending_bc
    )
    expected = np.zeros((3, 3))
    expected[0, 0] = beam.axial_boundary_matrix(
        arbitrary_omega, length, properties, axial_bc
    )[0, 0]
    expected[1:, 1:] = beam.bending_boundary_matrix(
        arbitrary_omega, length, properties, bending_bc
    )
    assert_allclose(combined, expected, rtol=2.0e-14, atol=2.0e-14)

    axial_modes = beam.exact_axial_modes(properties, length, axial_bc, 6)
    assert len(bending_roots) == 6
    assert all(root.accepted and root.sigma_ratio <= 1.0e-9 for root in bending_roots)
    axial_omegas = [mode.omega for mode in axial_modes]
    bending_omegas = [root.omega for root in bending_roots]
    clusters = beam.combined_spectrum_union(
        axial_omegas, bending_omegas, atol=1.0e-13, rtol=1.0e-12
    )
    members = [member for cluster in clusters for member in cluster.members]
    assert len(members) == 12
    assert_allclose(
        [member.omega for member in members],
        sorted(axial_omegas + bending_omegas),
        rtol=0.0,
        atol=0.0,
    )
    assert sum(member.subsystem == "axial" for member in members) == 6
    assert sum(member.subsystem == "bending" for member in members) == 6

    for family, omega in (
        [("axial", value) for value in axial_omegas]
        + [("bending", value) for value in bending_omegas]
    ):
        characteristic = beam.combined_boundary_matrix(
            omega, length, properties, axial_bc, bending_bc
        )
        _u, singular, vh = np.linalg.svd(characteristic)
        assert singular[-1] / singular[0] <= 2.0e-9
        assert np.count_nonzero(singular <= 2.0e-8 * singular[0]) == 1
        null_vector = vh[-1]
        if family == "axial":
            assert np.linalg.norm(null_vector[1:]) <= 1.0e-10
        else:
            assert abs(float(null_vector[0])) <= 1.0e-10


def test_combined_exact_degeneracy_nullity_and_subspace_semantics() -> None:
    base = _properties()
    length = 1.0
    bending_mode = beam.hinged_hinged_algebraic_modes(base, length, 1)[0]
    tuned_A = base.m * (bending_mode.omega * length / math.pi) ** 2
    properties = replace(base, A=tuned_A)
    axial_mode = beam.exact_axial_modes(properties, length, "FF", 1)[0]
    assert axial_mode.omega == pytest.approx(bending_mode.omega, rel=3.0e-15)

    combined = beam.combined_boundary_matrix(
        bending_mode.omega, length, properties, "FF", "HH"
    )
    singular = np.linalg.svd(combined, compute_uv=False)
    assert singular[-1] <= 1.0e-12 * singular[0]
    assert singular[-2] <= 1.0e-12 * singular[0]
    assert np.linalg.matrix_rank(combined, tol=1.0e-10) == 1
    full = beam.combined_full_boundary_matrix(
        bending_mode.omega, length, properties, "FF", "HH"
    )
    assert np.linalg.matrix_rank(full, tol=1.0e-10) == 4

    exact_clusters = beam.combined_spectrum_union(
        [axial_mode.omega],
        [bending_mode.omega],
        atol=1.0e-13,
        rtol=1.0e-13,
    )
    assert len(exact_clusters) == 1
    assert exact_clusters[0].multiplicity == 2
    assert exact_clusters[0].exact_degeneracy

    bending_initial = beam.hinged_hinged_exact_shape(
        bending_mode,
        properties,
        length,
        np.array([0.0, length]),
        mass_normalize=False,
    ).states[0]
    subsystem_basis = beam.combined_degeneracy_subspace(
        axial_vectors=np.array([0.0, 1.0]),
        bending_vectors=bending_initial,
    )
    assert subsystem_basis.shape == (6, 2)
    assert_allclose(subsystem_basis.T @ subsystem_basis, np.eye(2), atol=2.0e-15)
    assert np.count_nonzero(np.abs(subsystem_basis[[0, 3], 1]) > 1.0e-14) == 0
    assert np.count_nonzero(
        np.abs(subsystem_basis[[1, 2, 4, 5], 0]) > 1.0e-14
    ) == 0

    rotation = np.array(
        [[math.cos(0.37), -math.sin(0.37)], [math.sin(0.37), math.cos(0.37)]]
    )
    mixed_basis = subsystem_basis @ rotation
    assert_allclose(
        mixed_basis @ mixed_basis.T,
        subsystem_basis @ subsystem_basis.T,
        rtol=2.0e-15,
        atol=2.0e-15,
    )


def test_required_combined_exact_degeneracy_preserves_nullspace_and_mass_metric(
    required_combined_case: tuple[
        beam.BeamProperties, str, str, tuple[beam.RootDiagnostic, ...]
    ],
) -> None:
    base, axial_bc, bending_bc, bending_roots = required_combined_case
    length = 1.0
    target = bending_roots[0].omega
    axial_wavenumber = math.pi if axial_bc == "FF" else 0.5 * math.pi
    properties = replace(
        base,
        A=base.m * (target / axial_wavenumber) ** 2,
    )
    axial_mode = beam.exact_axial_modes(properties, length, axial_bc, 1)[0]
    assert axial_mode.omega == pytest.approx(target, rel=3.0e-15)

    reduced = beam.combined_boundary_matrix(
        target, length, properties, axial_bc, bending_bc
    )
    singular = np.linalg.svd(reduced, compute_uv=False)
    assert np.count_nonzero(singular <= 2.0e-8 * singular[0]) == 2
    full = beam.combined_full_boundary_matrix(
        target, length, properties, axial_bc, bending_bc
    )
    full_singular = np.linalg.svd(full, compute_uv=False)
    assert np.count_nonzero(full_singular <= 2.0e-8 * full_singular[0]) == 2

    exact = beam.combined_spectrum_union(
        [axial_mode.omega], [target], atol=1.0e-13, rtol=1.0e-13
    )
    assert len(exact) == 1
    assert exact[0].multiplicity == 2
    assert exact[0].exact_degeneracy

    x = np.linspace(0.0, length, 801)
    axial_shape = beam.axial_exact_shape(axial_mode, properties, length, x)
    bending_shape = beam.bending_mode_shape(
        target, properties, length, bending_bc, x
    )
    embedded_axial = beam.embed_axial_shape(axial_shape.states)
    embedded_bending = beam.embed_bending_shape(bending_shape.states)
    subsystem_basis = beam.combined_degeneracy_subspace(
        axial_vectors=axial_shape.states[0],
        bending_vectors=bending_shape.states[0],
    )
    assert subsystem_basis.shape == (6, 2)
    assert_allclose(subsystem_basis.T @ subsystem_basis, np.eye(2), atol=3.0e-15)
    assert_allclose(subsystem_basis[[1, 2, 4, 5], 0], 0.0, atol=0.0)
    assert_allclose(subsystem_basis[[0, 3], 1], 0.0, atol=0.0)

    inverse_scale = np.diag(
        1.0 / np.diag(beam.combined_state_scale(properties, length))
    )
    for column in subsystem_basis.T:
        scaled_initial = inverse_scale @ column
        relative_residual = np.linalg.norm(full @ scaled_initial) / (
            np.linalg.norm(full, ord=2) * np.linalg.norm(scaled_initial)
        )
        assert relative_residual <= 2.0e-9

    def mass_inner(first: np.ndarray, second: np.ndarray) -> float:
        density = properties.m * (
            first[:, 0] * second[:, 0] + first[:, 1] * second[:, 1]
        ) + properties.J * first[:, 2] * second[:, 2]
        return float(simpson(density, x=x))

    shape_basis = np.stack([embedded_axial, embedded_bending], axis=2)
    mass_gram = np.array(
        [
            [mass_inner(shape_basis[:, :, i], shape_basis[:, :, j]) for j in range(2)]
            for i in range(2)
        ]
    )
    assert_allclose(mass_gram, np.eye(2), rtol=2.0e-11, atol=2.0e-11)

    rotation = np.array(
        [[math.cos(0.41), -math.sin(0.41)], [math.sin(0.41), math.cos(0.41)]]
    )
    rotated_shapes = np.einsum("nsi,ij->nsj", shape_basis, rotation)
    rotated_mass_gram = np.array(
        [
            [
                mass_inner(rotated_shapes[:, :, i], rotated_shapes[:, :, j])
                for j in range(2)
            ]
            for i in range(2)
        ]
    )
    assert_allclose(rotated_mass_gram, np.eye(2), rtol=3.0e-11, atol=3.0e-11)
    for index in range(2):
        energies = beam.combined_modal_energies(
            x, rotated_shapes[:, :, index], properties, target
        )
        assert energies.modal_mass == pytest.approx(1.0, rel=3.0e-11)
        assert energies.identity_relative_error <= 2.0e-8


def test_smoke_runner_cannot_pass_unexecuted_combined_gate(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from scripts.analysis.laminated_beams import (  # noqa: PLC0415
        validate_reddy_symmetric_single_beam as runner,
    )

    tier_summary = {
        "comparison_count": 0,
        "pass_count": 0,
        "fail_count": 0,
        "maximum_absolute_error": 0.0,
        "failures": [],
    }
    monkeypatch.setattr(
        runner, "validate_source_contract", lambda _output: {"contract": {}}
    )
    monkeypatch.setattr(runner, "_scientific_module", lambda: object())
    monkeypatch.setattr(
        runner,
        "_laminate_validation",
        lambda *_args: (
            {"cross_ply_0_90_s": (object(), object())},
            [{"status": "PASS"}],
            {
                "pass": True,
                "maximum_symmetry_relative": 0.0,
                "maximum_compliance_schur_relative": 0.0,
                "maximum_partition_difference": 0.0,
            },
        ),
    )
    monkeypatch.setattr(
        runner,
        "_source_comparisons",
        lambda *_args: (
            [{"status": "PASS"}],
            [
                {
                    "sigma_ratio": 0.0,
                    "condition_number": 1.0,
                    "accepted": True,
                }
            ],
            {"A": dict(tier_summary), "B": dict(tier_summary)},
        ),
    )
    monkeypatch.setattr(
        runner,
        "_hinged_hinged_gate",
        lambda *_args: (
            [
                {
                    "sigma_ratio": 0.0,
                    "condition_number": 1.0,
                    "accepted": True,
                }
            ],
            [{"status": "PASS"}],
            {
                "pass": True,
                "maximum_analytic_transfer_relative": 0.0,
                "maximum_one_minus_MAC": 0.0,
                "maximum_energy_identity_relative": 0.0,
            },
        ),
    )
    monkeypatch.setattr(
        runner,
        "_eq_4_3_51_check",
        lambda *_args: {"status": "PASS", "relative_difference": 0.0},
    )
    monkeypatch.setattr(
        runner,
        "_axial_gate",
        lambda *_args: (
            [{"status": "PASS"}],
            [{"status": "PASS"}],
            {"pass": True, "maximum_analytic_transfer_relative": 0.0},
        ),
    )
    monkeypatch.setattr(
        runner,
        "_limit_checks",
        lambda *_args: ([{"status": "PASS"}], {"pass": True}),
    )
    monkeypatch.setattr(
        runner, "_scale_invariance", lambda *_args: {"status": "PASS"}
    )
    written_csv: dict[str, list[dict[str, object]]] = {}
    monkeypatch.setattr(
        runner,
        "_write_csv",
        lambda path, rows: written_csv.__setitem__(path.name, list(rows)),
    )
    monkeypatch.setattr(runner, "_write_json", lambda *_args: None)

    summary = runner.run_full_validation(tmp_path, smoke=True)
    assert summary["combined"].get("smoke_subset") is True
    assert summary["combined"]["pass"] is False
    assert summary["statuses"]["RLB-0C-COMBINED-UNION"] != "PASS"
    assert summary["statuses"]["OVERALL"] != "PASS_WITH_SOURCE_QUALIFICATIONS"
    for filename in (
        "combined_roots.csv",
        "combined_union_comparison.csv",
        "mode_shapes_combined.csv",
    ):
        assert written_csv[filename][0]["status"] == "SMOKE_NOT_RUN"


def test_import_isolation_from_yartsev_and_circular_workflows() -> None:
    module_path = ROOT / "scripts" / "lib" / "reddy_symmetric_laminated_beam.py"
    tree = ast.parse(module_path.read_text(encoding="utf-8"))
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    forbidden_fragments = (
        "yartsev",
        "circular",
        "variable_length_timoshenko",
        "analytic_branch_tracking",
    )
    assert not {
        name
        for name in imported
        if any(fragment in name.lower() for fragment in forbidden_fragments)
    }
    source_lower = module_path.read_text(encoding="utf-8").lower()
    assert "scripts.lib.yartsev" not in source_lower
    assert "scripts.lib.circular" not in source_lower


def test_results_directory_is_gitignored_when_git_metadata_is_available() -> None:
    git = shutil.which("git")
    if git is None or not (ROOT / ".git").exists():
        pytest.skip("Git metadata is unavailable.")
    completed = subprocess.run(
        [
            git,
            "check-ignore",
            "-q",
            "results/laminated_beams/reddy_symmetric_single_beam/",
        ],
        cwd=ROOT,
        check=False,
    )
    assert completed.returncode == 0
