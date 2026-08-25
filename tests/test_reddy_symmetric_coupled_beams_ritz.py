"""Focused gates for the independent two-arm RLB Rayleigh--Ritz model.

The frozen RLB-1 beta=0 transfer roots are deliberately used only after each
Ritz eigenproblem has been assembled and solved.  The N=16 guard does not pass
the required first-13 bridge tolerance, so the nonzero-beta spectral, MAC, and
symmetry tests below are explicitly skipped rather than bypassing that gate.
"""

from __future__ import annotations

import ast
import inspect
import math
from pathlib import Path

import numpy as np
from numpy.testing import assert_allclose
import pytest

from scripts.analysis.laminated_beams import (
    pilot_reddy_symmetric_coupled_beams_beta0 as beta0_pilot,
)
from scripts.lib import reddy_symmetric_coupled_beams_ritz as ritz
from scripts.lib import reddy_symmetric_laminated_beam as single


ROOT = Path(__file__).resolve().parents[1]
BRIDGE_TOLERANCE = 1.0e-8
NONZERO_BLOCK_REASON = (
    "RLB-1C0 beta=0 hard gate fails for the first 13 roots at the allowed "
    "N=16 guard; the specification forbids proceeding to nonzero-beta spectra"
)

# Existing committed RLB-1 beta=0 transfer inventories through root 13.
FROZEN_CROSS_PLY_ROOTS_THROUGH_GUARD = np.array(
    [
        22.67206118744846,
        50.61286710541817,
        83.29900985509377,
        117.87342669122222,
        153.2219871010395,
        188.6755267183602,
        224.03164905362817,
        226.78523261343506,
        259.2125392941688,
        294.2274722772023,
        329.0885957451078,
        363.8240553255954,
        398.4504559370047,
    ]
)
FROZEN_STEPPED_ROOTS_THROUGH_GUARD = np.array(
    [
        23.9158372842782,
        53.91664489847483,
        89.3108965084892,
        126.2839945093899,
        165.51250350976832,
        202.95573042836187,
        242.94976674041186,
        270.15742581052865,
        279.687109257467,
        319.49507406026686,
        356.2164998331029,
        395.0881961705233,
        432.35104074013026,
    ]
)


def _synthetic_properties(*, stepped: bool = False) -> single.BeamProperties:
    if stepped:
        return single.BeamProperties(
            A=8.2, D=0.22, S=1.1, m=1.3, J=0.002, K=5.0 / 6.0
        )
    return single.BeamProperties(
        A=24.7, D=0.81, S=2.6, m=1.3, J=0.002, K=5.0 / 6.0
    )


@pytest.fixture(scope="module")
def beta0_cases():
    """Return the two frozen bridge cases without reading transfer roots."""

    _manifest, selected = beta0_pilot._selected_benchmarks()
    homogeneous, stepped = beta0_pilot._case_specs(selected)
    cross_ply = next(
        case
        for case in homogeneous
        if case.case_id == "homogeneous__cross_ply_0_90_s__equal"
    )
    material_step = next(
        case
        for case in stepped
        if case.case_id == "stepped__0_deg__cross_ply"
    )
    return {"cross_ply": cross_ply, "stepped": material_step}


@pytest.fixture(scope="module")
def beta0_spectra(beta0_cases):
    """Solve Ritz orders independently before any frozen-root comparison."""

    result = {}
    for case_name, case in beta0_cases.items():
        result[case_name] = {
            order: ritz.solve_coupled_ritz_spectrum(
                case.properties_1,
                case.length_1,
                case.properties_2,
                case.length_2,
                0.0,
                order,
            )
            for order in (8, 10, 12, 16)
        }
    return result


def _relative_frequency_errors(values: np.ndarray, reference: np.ndarray) -> np.ndarray:
    return np.abs(values - reference) / np.maximum(
        np.abs(reference), np.finfo(float).tiny
    )


def test_ritz_module_imports_only_frozen_rlb_inputs() -> None:
    path = ROOT / "scripts" / "lib" / "reddy_symmetric_coupled_beams_ritz.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    project_modules: set[str] = set()
    imported_names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            project_modules.update(
                alias.name for alias in node.names if alias.name.startswith("scripts.")
            )
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith("scripts."):
                project_modules.add(node.module)
                imported_names.update(alias.name for alias in node.names)

    assert project_modules == {
        "scripts.lib",
        "scripts.lib.reddy_symmetric_laminated_beam",
    }
    assert imported_names == {"BeamProperties", "reddy_inplane_geometry"}
    forbidden = ("yartsev", "circular", "shellbuckling", "coupled_rods")
    import_text = " ".join(sorted(project_modules | imported_names)).lower()
    assert not any(fragment in import_text for fragment in forbidden)

    call_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }
    call_attributes = {
        node.func.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
    }
    assert "expm" not in call_names | call_attributes


def test_shifted_legendre_basis_imposes_only_outer_clamp_value() -> None:
    points = np.array([0.0, 0.25, 1.0])
    basis = ritz.shifted_legendre_clamped_basis(16, points)

    assert basis.values.shape == (16, 3)
    assert basis.derivatives.shape == (16, 3)
    assert_allclose(basis.values[:, 0], 0.0, atol=0.0, rtol=0.0)
    assert_allclose(basis.values[:, -1], 1.0, atol=2.0e-14, rtol=0.0)
    assert np.linalg.matrix_rank(
        np.column_stack((basis.values[:, -1], basis.derivatives[:, -1]))
    ) == 2


def test_fixed_quadrature_exactly_integrates_all_required_polynomial_degrees() -> None:
    nodes, weights = ritz.gauss_legendre_rule()
    assert len(nodes) == len(weights) == ritz.GAUSS_LEGENDRE_ORDER == 64
    # For N=16 the largest energy-product degree is 2N=32.  The fixed rule is
    # exact through degree 2*64-1=127, so this checks the entire required range.
    for degree in range(2 * ritz.RITZ_GUARD_ORDER + 1):
        computed = float(np.dot(weights, nodes**degree))
        assert computed == pytest.approx(1.0 / (degree + 1), abs=5.0e-14)
    with pytest.raises(ValueError, match="fixes Gauss--Legendre order"):
        ritz.gauss_legendre_rule(48)


@pytest.mark.parametrize("beta_deg", [0.0, 30.0, 90.0, -30.0])
def test_global_endpoint_maps_and_constraint_have_exact_three_row_contract(
    beta_deg: float,
) -> None:
    beta_rad = math.radians(beta_deg)
    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)
    expected_1 = np.array([[1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0]])
    expected_2 = np.array(
        [[-cosine, -sine, 0.0], [-sine, cosine, 0.0], [0.0, 0.0, 1.0]]
    )
    assert_allclose(ritz.global_kinematic_map(beta_rad, 1), expected_1, atol=2e-15)
    assert_allclose(ritz.global_kinematic_map(beta_rad, 2), expected_2, atol=2e-15)

    constraint = ritz.kinematic_constraint_matrix(beta_rad, 6)
    assert constraint.shape == (3, 36)
    assert np.linalg.matrix_rank(constraint) == 3
    assert list(inspect.signature(ritz.kinematic_constraint_matrix).parameters) == [
        "beta_rad",
        "order",
    ]
    source = inspect.getsource(ritz.kinematic_constraint_matrix)
    assert "force_vector" not in source
    assert "joint_matrix" not in source
    assert "properties" not in source


def test_constraint_rank_and_svd_nullspace_for_fixed_and_random_angles() -> None:
    rng = np.random.default_rng(20260825)
    angles = [0.0, math.radians(30.0), math.radians(90.0), math.radians(-30.0)]
    angles.extend(rng.uniform(-math.pi / 2.0, math.pi / 2.0, size=25))
    for beta_rad in angles:
        reduction = ritz.constraint_nullspace(float(beta_rad), 8)
        assert reduction.rank == 3
        assert reduction.nullspace.shape == (48, 45)
        assert reduction.orthonormality_residual <= 1.0e-12
        assert reduction.constraint_nullspace_residual <= 1.0e-12
        assert_allclose(
            reduction.nullspace.T @ reduction.nullspace,
            np.eye(45),
            atol=2.0e-14,
            rtol=0.0,
        )
        assert_allclose(
            reduction.constraint @ reduction.nullspace,
            np.zeros((3, 45)),
            atol=2.0e-14,
            rtol=0.0,
        )


def test_ritz_energy_matrices_and_reduced_eigenproblem_pass_algebraic_gates() -> None:
    properties_1 = _synthetic_properties()
    properties_2 = _synthetic_properties(stepped=True)
    spectrum = ritz.solve_coupled_ritz_spectrum(
        properties_1, 1.7, properties_2, 2.1, 0.0, 8
    )

    matrices = spectrum.matrices
    assert matrices.stiffness.shape == matrices.mass.shape == (48, 48)
    assert matrices.stiffness_symmetry_residual <= 1.0e-12
    assert matrices.mass_symmetry_residual <= 1.0e-12
    assert_allclose(matrices.stiffness, matrices.stiffness.T, atol=0.0, rtol=0.0)
    assert_allclose(matrices.mass, matrices.mass.T, atol=0.0, rtol=0.0)
    assert matrices.full_stiffness_minimum_eigenvalue > 0.0
    assert matrices.full_mass_minimum_eigenvalue > 0.0
    assert spectrum.reduced_stiffness_spd
    assert spectrum.reduced_mass_spd
    assert spectrum.reduced_stiffness_minimum_eigenvalue > 0.0
    assert spectrum.reduced_mass_minimum_eigenvalue > 0.0
    assert spectrum.zero_or_negative_mode_count == 0
    assert np.all(np.diff(spectrum.omegas) >= 0.0)
    assert spectrum.maximum_eigenpair_backward_residual <= 1.0e-9
    assert spectrum.maximum_rayleigh_relative_residual <= 1.0e-9
    assert spectrum.mass_orthonormality_residual <= 1.0e-12
    assert spectrum.maximum_constraint_residual <= 1.0e-12
    assert_allclose(
        spectrum.coefficients.T @ matrices.mass @ spectrum.coefficients,
        np.eye(45),
        atol=2.0e-11,
        rtol=0.0,
    )


def test_mode_reconstruction_outer_clamps_constitutive_fields_and_energy() -> None:
    properties_1 = _synthetic_properties()
    properties_2 = _synthetic_properties(stepped=True)
    spectrum = ritz.solve_coupled_ritz_spectrum(
        properties_1, 1.7, properties_2, 2.1, 0.0, 8
    )
    mode = spectrum.coefficients[:, 0]
    points = np.array([0.0, 0.25, 1.0])
    arm1 = ritz.evaluate_arm_mode(mode, 8, 1, points, 1.7, properties_1)
    arm2 = ritz.evaluate_arm_mode(mode, 8, 2, points, 2.1, properties_2)

    for fields, properties in ((arm1, properties_1), (arm2, properties_2)):
        assert fields.u[0] == pytest.approx(0.0, abs=0.0)
        assert fields.w[0] == pytest.approx(0.0, abs=0.0)
        assert fields.psi[0] == pytest.approx(0.0, abs=0.0)
        assert_allclose(fields.N, properties.A * fields.u_prime, atol=0.0, rtol=0.0)
        assert_allclose(
            fields.Q,
            properties.S * (fields.w_prime + fields.psi),
            atol=0.0,
            rtol=0.0,
        )
        assert_allclose(fields.M, properties.D * fields.psi_prime, atol=0.0, rtol=0.0)

    residual = ritz.joint_residuals_from_ritz_mode(
        mode, 8, 0.0, 1.7, properties_1, 2.1, properties_2
    )
    assert residual.displacement_normalized <= 1.0e-12
    assert residual.rotation_normalized <= 1.0e-12
    assert float(mode @ spectrum.matrices.mass @ mode) == pytest.approx(1.0, abs=2e-12)
    energy = ritz.modal_energy_diagnostics(
        mode, spectrum.omegas[0], 8, 1.7, properties_1, 2.1, properties_2
    )
    assert energy.total_mass_norm == pytest.approx(1.0, abs=2.0e-12)
    assert energy.energy_identity_relative <= 1.0e-9
    assert sum(
        (energy.T_axial_share, energy.T_transverse_share, energy.T_rotation_share)
    ) == pytest.approx(1.0, abs=2.0e-14)
    assert sum(
        (energy.U_axial_share, energy.U_bending_share, energy.U_shear_share)
    ) == pytest.approx(1.0, abs=2.0e-14)


def test_low_mode_natural_joint_equilibrium_converges_without_constraint_rows(
    beta0_cases, beta0_spectra
) -> None:
    case = beta0_cases["cross_ply"]
    residuals = {}
    for order in (8, 12):
        spectrum = beta0_spectra["cross_ply"][order]
        residuals[order] = []
        for mode_index in (0, 1):
            check = ritz.joint_residuals_from_ritz_mode(
                spectrum.coefficients[:, mode_index],
                order,
                0.0,
                case.length_1,
                case.properties_1,
                case.length_2,
                case.properties_2,
            )
            residuals[order].append(max(check.force_normalized, check.moment_normalized))
            assert check.displacement_normalized <= 1.0e-12
            assert check.rotation_normalized <= 1.0e-12

    assert max(residuals[12]) <= 1.0e-8
    assert max(residuals[12]) < 1.0e-3 * max(residuals[8])


def test_mass_mac_sign_invariance_and_cluster_subspace_overlap() -> None:
    mass = np.diag([1.0, 2.0, 3.0, 4.0])
    first = np.array([1.0, 0.0, 1.0, 0.0])
    second = np.array([0.0, 1.0, 0.0, 1.0])
    assert ritz.mass_mac(first, first, mass) == pytest.approx(1.0)
    assert ritz.mass_mac(first, -first, mass) == pytest.approx(1.0)
    assert ritz.mass_mac(first, second, mass) == pytest.approx(0.0)

    left = np.column_stack((first, second))
    right = left @ np.array([[2.0, -1.0], [0.5, 3.0]])
    orthonormal = ritz.mass_orthonormalize(left, mass)
    assert_allclose(orthonormal.T @ mass @ orthonormal, np.eye(2), atol=2.0e-14)
    overlap = ritz.compare_mass_subspaces(left, right, mass)
    assert_allclose(overlap.singular_values, np.ones(2), atol=2.0e-14)
    assert overlap.minimum_singular_value == pytest.approx(1.0, abs=2.0e-14)
    assert overlap.maximum_principal_angle_rad <= 2.0e-7


@pytest.mark.parametrize(
    ("case_name", "reference"),
    [
        ("cross_ply", FROZEN_CROSS_PLY_ROOTS_THROUGH_GUARD),
        ("stepped", FROZEN_STEPPED_ROOTS_THROUGH_GUARD),
    ],
)
def test_beta0_n16_bridge_gate_failure_is_detected_not_hidden(
    case_name: str, reference: np.ndarray, beta0_cases, beta0_spectra
) -> None:
    case = beta0_cases[case_name]
    spectrum = beta0_spectra[case_name][16]
    ritz_roots = spectrum.omegas[:13] * case.frequency_scale
    relative = _relative_frequency_errors(ritz_roots, reference)
    violating_positions = tuple(np.flatnonzero(relative > BRIDGE_TOLERANCE) + 1)

    assert len(ritz_roots) == len(reference) == 13
    assert violating_positions == (12, 13)
    assert float(np.max(relative)) > BRIDGE_TOLERANCE
    assert not bool(np.all(relative <= BRIDGE_TOLERANCE))


@pytest.mark.parametrize("case_name", ["cross_ply", "stepped"])
def test_n16_guard_was_legitimately_triggered_by_n10_n12_nonconvergence(
    case_name: str, beta0_cases, beta0_spectra
) -> None:
    case = beta0_cases[case_name]
    roots_10 = beta0_spectra[case_name][10].omegas[:13] * case.frequency_scale
    roots_12 = beta0_spectra[case_name][12].omegas[:13] * case.frequency_scale
    relative = _relative_frequency_errors(roots_12, roots_10)
    assert float(np.max(relative)) > BRIDGE_TOLERANCE


def test_full_first13_natural_equilibrium_gate_is_also_reported_as_unpassed(
    beta0_cases, beta0_spectra
) -> None:
    maxima = []
    for case_name, case in beta0_cases.items():
        spectrum = beta0_spectra[case_name][16]
        for mode_index in range(13):
            check = ritz.joint_residuals_from_ritz_mode(
                spectrum.coefficients[:, mode_index],
                16,
                0.0,
                case.length_1,
                case.properties_1,
                case.length_2,
                case.properties_2,
            )
            maxima.append(max(check.force_normalized, check.moment_normalized))
    assert max(maxima) > 1.0e-8


@pytest.mark.skip(reason=NONZERO_BLOCK_REASON)
def test_beta30_transfer_ritz_spectrum_contract() -> None:
    """Blocked until RLB-1C0 passes; no nonzero transfer/Ritz solve is allowed."""


@pytest.mark.skip(reason=NONZERO_BLOCK_REASON)
def test_beta90_transfer_ritz_spectrum_contract() -> None:
    """Blocked until RLB-1C0 passes; no nonzero transfer/Ritz solve is allowed."""


@pytest.mark.skip(reason=NONZERO_BLOCK_REASON)
def test_nonzero_beta_isolated_mac_and_cluster_subspace_contract() -> None:
    """Blocked until RLB-1C0 passes; no nonzero mode comparison is allowed."""


@pytest.mark.skip(reason=NONZERO_BLOCK_REASON)
def test_beta_reflection_and_arm_exchange_symmetry_contract() -> None:
    """Blocked until RLB-1C0 passes; no nonzero symmetry spectrum is allowed."""
