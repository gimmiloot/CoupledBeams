"""Targeted gates for the narrow RLB rigid-joint and beta=0 assembly.

The tests in this file deliberately avoid a physical spectral sweep.  They
check the invariant joint construction, virtual work, transfer assembly,
algebraic beta=0 reference identities, and the frozen physical roots through
the guard.  Root-search completeness and cluster handling are exercised with
deterministic synthetic matrices whose complete spectra are prescribed.
"""

from __future__ import annotations

import ast
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from numpy.testing import assert_allclose
import pytest
from scipy.linalg import expm

from scripts.analysis.laminated_beams import (
    pilot_reddy_symmetric_coupled_beams_beta0 as pilot,
)
from scripts.lib import reddy_symmetric_coupled_beams as coupled
from scripts.lib import reddy_symmetric_laminated_beam as single


ROOT = Path(__file__).resolve().parents[1]
JOINT_EQUALITY_TOLERANCE = 1.0e-14
VIRTUAL_WORK_TOLERANCE = 1.0e-12

FROZEN_0_DEG_ROOTS_THROUGH_GUARD = (
    25.32726452561564,
    57.75894427342482,
    96.10202658801227,
    137.05445711814514,
    179.20916099762718,
    221.6837223945868,
    264.12944283329995,
    306.4007476212575,
    314.1592653589793,
    348.4742600341009,
    390.35388421759313,
    432.06558308684197,
    473.62996426712886,
)
FROZEN_CROSS_PLY_ROOTS_THROUGH_GUARD = (
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
)
FROZEN_STEPPED_ROOTS_THROUGH_GUARD = (
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
)


def _properties(*, stepped: bool = False) -> single.BeamProperties:
    """Return moderate positive synthetic RLB properties for matrix tests."""

    if stepped:
        return single.BeamProperties(
            A=8.2,
            D=0.22,
            S=1.1,
            m=1.3,
            J=0.002,
            K=5.0 / 6.0,
        )
    return single.BeamProperties(
        A=24.7,
        D=0.81,
        S=2.6,
        m=1.3,
        J=0.002,
        K=5.0 / 6.0,
    )


def _expected_joint_matrix(beta_rad: float) -> np.ndarray:
    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)
    return np.array(
        [
            [1, 0, 0, 0, 0, 0, cosine, sine, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, -sine, cosine, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, -cosine, -sine, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, sine, -cosine, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        ],
        dtype=float,
    )


def _compatible_arm2_state(beta_rad: float, arm1_state: np.ndarray) -> np.ndarray:
    """Construct state 2 from compatibility and physical equilibrium."""

    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)
    u1, w1, psi1, N1, Q1, M1 = np.asarray(arm1_state, dtype=float)
    return np.array(
        [
            -cosine * u1 + sine * w1,
            -sine * u1 - cosine * w1,
            psi1,
            cosine * N1 - sine * Q1,
            sine * N1 + cosine * Q1,
            -M1,
        ]
    )


def _relative_difference(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), np.finfo(float).tiny)


def test_coupled_module_imports_only_the_two_canonical_rlb_helpers() -> None:
    module_path = ROOT / "scripts" / "lib" / "reddy_symmetric_coupled_beams.py"
    source = module_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    project_imports: set[str] = set()
    imported_aliases: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("scripts."):
                    project_imports.add(alias.name)
        elif isinstance(node, ast.ImportFrom) and node.module:
            if node.module.startswith("scripts."):
                project_imports.add(node.module)
                imported_aliases.update(alias.name for alias in node.names)

    assert project_imports == {"scripts.lib"}
    assert imported_aliases == {
        "reddy_inplane_geometry",
        "reddy_symmetric_laminated_beam",
    }
    forbidden = (
        "yartsev",
        "circular",
        "shellbuckling",
        "variable_length_timoshenko",
    )
    imported_text = " ".join(
        [*project_imports, *imported_aliases]
    ).lower()
    assert not any(fragment in imported_text for fragment in forbidden)


@pytest.mark.parametrize("beta_deg", [-90.0, -30.0, 0.0, 30.0, 90.0])
def test_invariant_closed_form_and_production_joint_matrices_are_identical(
    beta_deg: float,
) -> None:
    beta_rad = math.radians(beta_deg)
    invariant = coupled.joint_matrix_from_physical_maps(beta_rad)
    closed = coupled.joint_matrix_closed_form(beta_rad)

    assert_allclose(
        invariant,
        closed,
        rtol=0.0,
        atol=JOINT_EQUALITY_TOLERANCE,
    )
    assert_allclose(closed, _expected_joint_matrix(beta_rad), rtol=0.0, atol=0.0)
    assert_allclose(coupled.joint_matrix(beta_rad), invariant, rtol=0.0, atol=0.0)
    assert np.linalg.matrix_rank(invariant, tol=1.0e-12) == 6


def test_joint_builder_equality_for_one_hundred_fixed_seed_random_angles() -> None:
    rng = np.random.default_rng(2026082501)
    maximum = 0.0
    for beta_rad in rng.uniform(-0.5 * math.pi, 0.5 * math.pi, size=100):
        invariant = coupled.joint_matrix_from_physical_maps(float(beta_rad))
        closed = coupled.joint_matrix_closed_form(float(beta_rad))
        maximum = max(maximum, float(np.max(np.abs(invariant - closed))))
        assert np.linalg.matrix_rank(invariant, tol=1.0e-12) == 6
    assert maximum <= JOINT_EQUALITY_TOLERANCE


def test_joint_builders_require_a_finite_radian_angle() -> None:
    for builder in (
        coupled.joint_matrix_from_physical_maps,
        coupled.joint_matrix_closed_form,
        coupled.joint_matrix,
    ):
        with pytest.raises(ValueError, match="beta_rad"):
            builder(math.nan)
        with pytest.raises(ValueError, match="beta_rad"):
            builder(math.inf)


def test_beta_zero_R0_identity_for_random_joint_states() -> None:
    expected_R0 = np.diag([-1.0, -1.0, 1.0, 1.0, 1.0, -1.0])
    assert_allclose(coupled.R0_BETA0, expected_R0, rtol=0.0, atol=0.0)

    joint = coupled.joint_matrix(0.0)
    rng = np.random.default_rng(2026082502)
    for _ in range(100):
        state_1 = rng.normal(size=6)
        state_2 = coupled.R0_BETA0 @ state_1
        assert_allclose(
            joint @ np.concatenate((state_1, state_2)),
            np.zeros(6),
            rtol=0.0,
            atol=2.0e-15,
        )


def test_beta_zero_and_ninety_scalar_limits_are_canonical() -> None:
    beta_zero = coupled.joint_matrix_closed_form(0.0)
    assert_allclose(beta_zero, _expected_joint_matrix(0.0), rtol=0.0, atol=0.0)

    beta_ninety = coupled.joint_matrix_closed_form(0.5 * math.pi)
    expected_ninety = np.array(
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
    assert_allclose(beta_ninety, expected_ninety, rtol=0.0, atol=7.0e-17)


def test_physical_and_scalar_joint_residuals_vanish_for_compatible_states() -> None:
    rng = np.random.default_rng(2026082503)
    beta_values = np.concatenate(
        (
            np.radians([0.0, 30.0, 90.0]),
            rng.uniform(-0.5 * math.pi, 0.5 * math.pi, size=197),
        )
    )
    for beta_rad in beta_values:
        state_1 = rng.normal(size=6)
        state_2 = _compatible_arm2_state(float(beta_rad), state_1)
        residual = coupled.physical_joint_residuals(beta_rad, state_1, state_2)
        for vector in (
            residual.displacement,
            residual.rotation,
            residual.force,
            residual.moment,
        ):
            assert np.linalg.norm(vector) <= 4.0e-15
        scalar = coupled.scalar_joint_residuals_from_physical_maps(
            beta_rad, state_1, state_2
        )
        matrix_scalar = coupled.joint_matrix(beta_rad) @ np.concatenate(
            (state_1, state_2)
        )
        assert_allclose(scalar, matrix_scalar, rtol=0.0, atol=3.0e-15)
        assert np.linalg.norm(scalar) <= 6.0e-15


def test_virtual_work_gate_for_one_thousand_fixed_seed_cases() -> None:
    rng = np.random.default_rng(2026082504)
    beta_values = np.concatenate(
        (
            np.radians([0.0, 30.0, 90.0]),
            rng.uniform(-0.5 * math.pi, 0.5 * math.pi, size=997),
        )
    )
    maximum_normalized = 0.0
    maximum_pairing = 0.0
    for beta_rad in beta_values:
        state_1 = rng.normal(size=6)
        state_2 = _compatible_arm2_state(float(beta_rad), state_1)
        check = coupled.joint_virtual_work_check(
            beta_rad,
            state_1[3:],
            state_2[3:],
            state_1[:3],
            state_2[:3],
        )
        maximum_normalized = max(maximum_normalized, check.normalized_residual)
        maximum_pairing = max(maximum_pairing, check.pairing_absolute_difference)
        assert check.pairing_absolute_difference <= 2.0e-14
        assert check.normalized_residual <= VIRTUAL_WORK_TOLERANCE
        assert abs(check.local_total) == pytest.approx(check.absolute_residual)
        assert check.global_total == pytest.approx(
            check.local_total, rel=0.0, abs=2.0e-14
        )
    assert maximum_normalized <= VIRTUAL_WORK_TOLERANCE
    assert maximum_pairing <= 2.0e-14


def test_local_global_pairing_is_preserved_without_joint_constraints() -> None:
    rng = np.random.default_rng(2026082505)
    for beta_rad in (0.0, math.pi / 6.0, math.pi / 2.0):
        resultants_1 = rng.normal(size=3)
        resultants_2 = rng.normal(size=3)
        motion_1 = rng.normal(size=3)
        motion_2 = rng.normal(size=3)
        check = coupled.joint_virtual_work_check(
            beta_rad, resultants_1, resultants_2, motion_1, motion_2
        )
        assert check.pairing_absolute_difference <= 4.0e-15


def test_clamp_basis_selector_rank_and_outer_kinematic_residual() -> None:
    expected_basis = np.vstack((np.zeros((3, 3)), np.eye(3)))
    expected_selector = np.hstack((np.eye(3), np.zeros((3, 3))))
    assert_allclose(coupled.CLAMP_BASIS, expected_basis, rtol=0.0, atol=0.0)
    assert_allclose(coupled.CLAMP_SELECTOR, expected_selector, rtol=0.0, atol=0.0)
    assert np.linalg.matrix_rank(coupled.CLAMP_BASIS) == 3
    assert np.linalg.matrix_rank(coupled.CLAMP_SELECTOR) == 3
    assert_allclose(
        coupled.CLAMP_SELECTOR @ coupled.CLAMP_BASIS,
        np.zeros((3, 3)),
        rtol=0.0,
        atol=0.0,
    )

    rng = np.random.default_rng(2026082506)
    for _ in range(20):
        initial_state = coupled.CLAMP_BASIS @ rng.normal(size=3)
        assert_allclose(initial_state[:3], np.zeros(3), rtol=0.0, atol=0.0)


@pytest.mark.parametrize("stepped", [False, True])
def test_arm_state_transfer_and_clamp_map_reuse_single_beam_physics(
    stepped: bool,
) -> None:
    properties = _properties(stepped=stepped)
    omega = 0.137
    length = 3.2
    state_matrix = coupled.arm_state_matrix(omega, properties)

    assert_allclose(
        state_matrix,
        single.combined_state_matrix(omega, properties),
        rtol=0.0,
        atol=0.0,
    )
    assert_allclose(
        coupled.arm_transfer_matrix(omega, length, properties),
        single.combined_transfer_matrix(omega, length, properties),
        rtol=0.0,
        atol=0.0,
    )
    assert_allclose(
        coupled.arm_transfer_matrix(omega, length, properties),
        expm(state_matrix * length),
        rtol=3.0e-14,
        atol=3.0e-14,
    )
    assert_allclose(
        coupled.arm_clamp_map(omega, length, properties),
        coupled.arm_transfer_matrix(omega, length, properties)
        @ coupled.CLAMP_BASIS,
        rtol=0.0,
        atol=0.0,
    )


def test_coupled_clamp_map_and_boundary_matrix_have_exact_block_assembly() -> None:
    properties_1 = _properties()
    properties_2 = _properties(stepped=True)
    omega = 0.071
    beta_rad = 0.0
    length_1, length_2 = 7.0, 13.0

    block = coupled.coupled_clamp_to_joint_map(
        omega, length_1, properties_1, length_2, properties_2
    )
    expected = np.zeros((12, 6))
    expected[:6, :3] = coupled.arm_clamp_map(omega, length_1, properties_1)
    expected[6:, 3:] = coupled.arm_clamp_map(omega, length_2, properties_2)
    assert_allclose(block, expected, rtol=0.0, atol=0.0)
    assert_allclose(
        coupled.coupled_boundary_matrix(
            omega,
            beta_rad,
            length_1,
            properties_1,
            length_2,
            properties_2,
        ),
        coupled.joint_matrix(beta_rad) @ expected,
        rtol=0.0,
        atol=0.0,
    )


def test_positive_equilibration_records_finite_positive_exact_factors() -> None:
    raw = np.array(
        [
            [1.0e-9, 2.0, 0.0],
            [3.0e7, -4.0e-4, 5.0],
            [0.0, 6.0, 7.0e3],
        ]
    )
    scaled = coupled.positively_equilibrate_matrix(raw)
    assert np.all(np.isfinite(scaled.row_factors))
    assert np.all(np.isfinite(scaled.column_factors))
    assert np.all(scaled.row_factors > 0.0)
    assert np.all(scaled.column_factors > 0.0)
    assert_allclose(
        scaled.scaled_matrix,
        np.diag(scaled.row_factors) @ raw @ np.diag(scaled.column_factors),
        rtol=2.0e-15,
        atol=2.0e-15,
    )
    assert np.linalg.det(scaled.scaled_matrix) == pytest.approx(
        np.linalg.det(raw) * scaled.determinant_factor,
        rel=3.0e-14,
    )


def test_positive_scaling_preserves_singular_zeros_and_recovers_raw_null_vector() -> None:
    raw = np.array(
        [
            [1.0e12, 0.0, 2.0],
            [0.0, 3.0e-8, 0.0],
            [0.0, 0.0, 0.0],
        ]
    )
    scaling = coupled.positively_equilibrate_matrix(raw)
    diagnostics = coupled.boundary_matrix_diagnostics(raw)

    # Use one explicit absolute algebraic threshold.  The default relative
    # threshold is intentionally changed by equilibration and is therefore not
    # an invariant rank definition for this deliberately ill-scaled example.
    assert np.linalg.matrix_rank(raw, tol=1.0e-20) == 2
    assert np.linalg.matrix_rank(scaling.scaled_matrix, tol=1.0e-20) == 2
    assert np.linalg.det(raw) == 0.0
    assert np.linalg.det(scaling.scaled_matrix) == 0.0
    assert diagnostics.raw_determinant == 0.0
    assert diagnostics.scaled_determinant == 0.0
    assert diagnostics.raw_sigma_ratio == 0.0
    assert diagnostics.scaled_sigma_ratio == 0.0
    assert diagnostics.relative_singular_residual == 0.0
    assert diagnostics.boundary_null_residual <= 2.0e-16
    assert np.linalg.norm(raw @ diagnostics.raw_right_null_vector) <= 2.0e-16


def test_scaling_floor_does_not_hide_a_nearly_zero_characteristic_row() -> None:
    raw = np.diag([2.0e-20, 1.0, 1.0e3])
    scaling = coupled.positively_equilibrate_matrix(raw)
    singular = np.linalg.svd(scaling.scaled_matrix, compute_uv=False)

    assert coupled.MATRIX_SCALING_RELATIVE_FLOOR == pytest.approx(
        float(np.finfo(float).eps ** 0.25), rel=0.0, abs=0.0
    )
    assert singular[-1] / singular[0] <= 2.0e-15


@pytest.mark.parametrize("stepped", [False, True])
def test_direct_fixed_fixed_scaling_recovers_the_exact_first_axial_root(
    stepped: bool,
) -> None:
    properties = _properties(stepped=stepped)
    length = 20.0
    omega = math.pi * math.sqrt(properties.A / properties.m) / length
    raw = coupled.direct_fixed_fixed_boundary_matrix(omega, length, properties)
    diagnostic = coupled.boundary_matrix_diagnostics(raw)

    assert abs(raw[0, 0]) <= 1.0e-14
    assert diagnostic.scaled_sigma_ratio <= 1.0e-12
    assert diagnostic.boundary_null_residual <= 1.0e-12


@pytest.mark.parametrize("split", [(10.0, 10.0), (7.0, 13.0)])
@pytest.mark.parametrize("stepped", [False, True])
def test_homogeneous_beta0_assembly_matches_one_fixed_fixed_beam_without_roots(
    split: tuple[float, float],
    stepped: bool,
) -> None:
    properties = _properties(stepped=stepped)
    length_1, length_2 = split
    total_length = length_1 + length_2
    right_local_kinematic_sign = np.diag([-1.0, -1.0, 1.0])

    for omega in (0.0, 0.01, 0.05, 0.1):
        assembled = coupled.coupled_boundary_matrix(
            omega, 0.0, length_1, properties, length_2, properties
        )
        stepped_form = coupled.direct_stepped_boundary_matrix(
            omega, length_1, properties, length_2, properties
        )
        one_beam = coupled.direct_fixed_fixed_boundary_matrix(
            omega, total_length, properties
        )

        assert_allclose(
            stepped_form,
            right_local_kinematic_sign @ one_beam,
            rtol=2.0e-11,
            atol=2.0e-9,
        )
        assert _relative_difference(
            float(np.linalg.det(assembled)), float(np.linalg.det(stepped_form))
        ) <= 2.0e-10


def test_stepped_beta0_coupled_and_independent_direct_determinants_agree() -> None:
    properties_1 = _properties()
    properties_2 = _properties(stepped=True)
    for omega in (0.0, 0.01, 0.05, 0.1):
        assembled = coupled.coupled_boundary_matrix(
            omega, 0.0, 10.0, properties_1, 10.0, properties_2
        )
        direct = coupled.direct_stepped_boundary_matrix(
            omega, 10.0, properties_1, 10.0, properties_2
        )
        assert _relative_difference(
            float(np.linalg.det(assembled)), float(np.linalg.det(direct))
        ) <= 1.0e-10


def test_equal_length_stepped_reflection_preserves_the_beta0_spectrum() -> None:
    properties_0 = _properties()
    properties_cross_ply = _properties(stepped=True)
    for omega in (0.0, 0.01, 0.05, 0.1):
        forward = coupled.direct_stepped_boundary_matrix(
            omega, 10.0, properties_0, 10.0, properties_cross_ply
        )
        reflected = coupled.direct_stepped_boundary_matrix(
            omega, 10.0, properties_cross_ply, 10.0, properties_0
        )
        assert _relative_difference(
            float(np.linalg.det(forward)), float(np.linalg.det(reflected))
        ) <= 1.0e-10


@pytest.mark.parametrize(
    ("laminate_id", "roots"),
    [
        ("0_deg", FROZEN_0_DEG_ROOTS_THROUGH_GUARD),
        ("cross_ply_0_90_s", FROZEN_CROSS_PLY_ROOTS_THROUGH_GUARD),
    ],
)
@pytest.mark.parametrize("split_id", ["equal", "unequal_35_65"])
def test_frozen_homogeneous_physical_roots_pass_coupled_and_direct_gates(
    laminate_id: str,
    roots: tuple[float, ...],
    split_id: str,
) -> None:
    _, selected = pilot._selected_benchmarks()
    homogeneous, _ = pilot._case_specs(selected)
    spec = next(
        item
        for item in homogeneous
        if item.order_id == f"{laminate_id}|{laminate_id}"
        and item.split_id == split_id
    )
    providers = (
        pilot._coupled_provider(coupled, spec),
        pilot._direct_homogeneous_provider(coupled, spec),
    )

    assert len(roots) == 13
    for omega_bar in roots:
        for provider in providers:
            diagnostic = pilot.boundary_matrix_diagnostics(
                omega_bar, provider, spec.frequency_scale
            )
            assert diagnostic.detected_nullity == 1
            assert diagnostic.scaled_sigma_ratio <= 1.0e-9
            assert diagnostic.raw_boundary_null_residual <= 1.0e-9


@pytest.mark.parametrize(
    "case_id",
    ["stepped__0_deg__cross_ply", "stepped__cross_ply__0_deg"],
)
def test_frozen_stepped_physical_roots_pass_coupled_direct_and_reflection_gates(
    case_id: str,
) -> None:
    _, selected = pilot._selected_benchmarks()
    _, stepped = pilot._case_specs(selected)
    spec = next(item for item in stepped if item.case_id == case_id)
    providers = (
        pilot._coupled_provider(coupled, spec),
        pilot._direct_stepped_provider(coupled, spec),
    )

    assert len(FROZEN_STEPPED_ROOTS_THROUGH_GUARD) == 13
    for omega_bar in FROZEN_STEPPED_ROOTS_THROUGH_GUARD:
        for provider in providers:
            diagnostic = pilot.boundary_matrix_diagnostics(
                omega_bar, provider, spec.frequency_scale
            )
            assert diagnostic.detected_nullity == 1
            assert diagnostic.scaled_sigma_ratio <= 1.0e-9
            assert diagnostic.raw_boundary_null_residual <= 1.0e-9


def test_direct_reference_builders_do_not_call_joint_or_coupled_assembly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def forbidden(*_args: object, **_kwargs: object) -> np.ndarray:
        raise AssertionError("direct reference called a forbidden joint/coupled helper")

    monkeypatch.setattr(coupled, "joint_matrix_from_physical_maps", forbidden)
    monkeypatch.setattr(coupled, "joint_matrix_closed_form", forbidden)
    monkeypatch.setattr(coupled, "joint_matrix", forbidden)
    monkeypatch.setattr(coupled, "coupled_boundary_matrix", forbidden)
    monkeypatch.setattr(coupled, "coupled_clamp_to_joint_map", forbidden)

    properties_1 = _properties()
    properties_2 = _properties(stepped=True)
    fixed = coupled.direct_fixed_fixed_boundary_matrix(0.03, 20.0, properties_1)
    stepped = coupled.direct_stepped_boundary_matrix(
        0.03, 10.0, properties_1, 10.0, properties_2
    )
    assert fixed.shape == (3, 3)
    assert stepped.shape == (3, 3)


def test_reference_detector_reconciliation_is_cross_method_and_unambiguous() -> None:
    policy = pilot.SearchPolicy()
    determinant = SimpleNamespace(
        omega=1.0,
        detection="determinant_bracket",
        boundary_residual=1.0e-16,
        sigma_ratio=1.0e-16,
    )
    svd = SimpleNamespace(
        omega=1.0 + 1.0e-9,
        detection="svd_minimum",
        boundary_residual=1.0e-10,
        sigma_ratio=1.0e-10,
    )
    groups, summary = pilot._reconcile_reference_detector_detections(
        [svd, determinant],
        frequency_scale=1.0,
        policy=policy,
    )
    assert len(groups) == 1
    assert len(groups[0]) == 2
    assert summary["reconciled_cross_method_count"] == 1

    second_determinant = SimpleNamespace(
        omega=1.0 + 2.0e-9,
        detection="determinant_bracket",
        boundary_residual=1.0e-16,
        sigma_ratio=1.0e-16,
    )
    ambiguous_svd = SimpleNamespace(
        omega=1.0 + 1.0e-9,
        detection="svd_minimum",
        boundary_residual=1.0e-10,
        sigma_ratio=1.0e-10,
    )
    ambiguous_groups, ambiguous_summary = (
        pilot._reconcile_reference_detector_detections(
            [determinant, second_determinant, ambiguous_svd],
            frequency_scale=1.0,
            policy=policy,
        )
    )
    assert len(ambiguous_groups) == 3
    assert sum(
        str(item.detection).startswith("determinant")
        for group in ambiguous_groups
        for item in group
    ) == 2
    assert ambiguous_summary["reconciled_cross_method_count"] == 0
    assert ambiguous_summary["ambiguous_cross_method_count"] == 1


def _synthetic_inventory_matrix(omega: float) -> np.ndarray:
    """Return a bounded 6x6 matrix with 13 slots in 12 root events.

    The roots are the integers 1 through 12.  Both independent 2x2 blocks are
    singular at five, giving exact nullity and multiplicity two without a
    determinant sign change.  All other roots have nullity one.
    """

    first_roots = np.array([1.0, 3.0, 5.0, 7.0, 9.0, 11.0])
    second_roots = np.array([2.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0])
    first = float(np.prod(np.tanh((float(omega) - first_roots) / 0.4)))
    second = float(np.prod(np.tanh((float(omega) - second_roots) / 0.4)))
    matrix = np.eye(6)
    matrix[:2, :2] = np.array([[1.0, 1.0], [1.0, 1.0 + first]])
    matrix[2:4, 2:4] = np.array([[1.0, 1.0], [1.0, 1.0 + second]])
    return matrix


def test_synthetic_inventory_detects_sign_roots_exact_multiplicity_and_guard() -> None:
    policy = pilot.SearchPolicy(
        omega_bar_min=0.0,
        omega_bar_max=14.0,
        primary_scan_points=1401,
        verification_scan_points=2801,
        root_xtol_bar=1.0e-13,
        post_guard_tail_bar=1.0,
    )
    inventory = pilot.seed_free_root_inventory(
        _synthetic_inventory_matrix,
        1.0,
        policy,
        case_id="synthetic_nullity_inventory",
        builder_id="prescribed_6x6_matrix",
    )

    assert inventory.status == "PASS"
    assert inventory.independent_agreement is True
    assert inventory.guard_available is True
    assert inventory.guard_not_at_scan_boundary is True
    assert inventory.unresolved_low_sigma_count == 0
    assert inventory.maximum_primary_verification_relative <= 1.0e-9
    assert len(inventory.slots) == 13
    assert [slot.sorted_slot for slot in inventory.slots] == list(range(1, 14))
    assert [slot.role for slot in inventory.slots[:12]] == ["FIRST_12"] * 12
    assert inventory.slots[12].role == "ROOT_13_GUARD"

    expected_slots = [1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 10, 11, 12]
    assert_allclose(
        [slot.event.omega_bar for slot in inventory.slots],
        expected_slots,
        rtol=0.0,
        atol=2.0e-9,
    )
    exact = [event for event in inventory.primary.events if abs(event.omega_bar - 5.0) < 2.0e-9]
    assert len(exact) == 1
    assert exact[0].multiplicity == 2
    assert exact[0].detected_nullity == 2
    assert exact[0].cluster_semantics == "EXACT_DEGENERATE_SUBSPACE"
    assert exact[0].cluster_multiplicity == 2
    assert exact[0].cluster_total_nullity == 2
    assert [slot.repeated_root_slot for slot in inventory.slots[4:6]] == [1, 2]


def test_synthetic_exact_double_root_has_no_sign_change_but_nullity_two() -> None:
    left = np.linalg.det(_synthetic_inventory_matrix(5.0 - 1.0e-4))
    right = np.linalg.det(_synthetic_inventory_matrix(5.0 + 1.0e-4))
    assert left * right > 0.0

    diagnostic = pilot.boundary_matrix_diagnostics(
        5.0,
        _synthetic_inventory_matrix,
        1.0,
    )
    assert diagnostic.detected_nullity == 2
    assert diagnostic.root_gate_nullity == 2
    assert diagnostic.scaled_sigma_ratio <= 1.0e-15
    assert diagnostic.raw_boundary_null_residual <= 1.0e-15


def test_local_offset_polish_resolves_a_narrow_nullity_gate() -> None:
    root = math.pi

    def narrow_row_matrix(omega: float) -> np.ndarray:
        return np.diag([math.sin(float(omega)), 1.0, 1.0e3])

    policy = pilot.SearchPolicy(
        omega_bar_min=3.0,
        omega_bar_max=3.3,
        primary_scan_points=101,
        verification_scan_points=201,
        root_xtol_bar=1.0e-11,
    )
    evaluator = pilot._DiagnosticEvaluator(narrow_row_matrix, 1.0, policy)
    candidate = pilot._root_candidate(
        evaluator=evaluator,
        policy=policy,
        case_id="narrow_polish",
        builder_id="prescribed_diagonal_matrix",
        scan_id="unit",
        source="determinant_bracket",
        left=3.0,
        right=3.3,
        omega_bar=root + 2.0e-12,
        interior=True,
    )

    assert candidate.accepted is True
    assert candidate.diagnostics.detected_nullity == 1
    assert candidate.diagnostics.scaled_sigma_ratio <= 1.0e-12
    assert abs(candidate.omega_bar - root) <= 2.0e-13


_CLOSE_ROOT_VALUES = np.array(
    [1.0, 2.0, 3.0, 4.0, 5.0, 5.005, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
)


def _synthetic_close_root_matrix(omega: float) -> np.ndarray:
    """Return 13 simple roots, including a resolved near-degenerate pair."""

    first_roots = _CLOSE_ROOT_VALUES[::2]
    second_roots = _CLOSE_ROOT_VALUES[1::2]
    width = 0.4
    first = float(np.prod(np.tanh((float(omega) - first_roots) / width)))
    second = float(np.prod(np.tanh((float(omega) - second_roots) / width)))
    matrix = np.eye(6)
    matrix[:2, :2] = np.array([[1.0, 1.0], [1.0, 1.0 + first]])
    matrix[2:4, 2:4] = np.array([[1.0, 1.0], [1.0, 1.0 + second]])
    return matrix


def test_synthetic_inventory_resolves_close_roots_as_a_near_cluster() -> None:
    policy = pilot.SearchPolicy(
        omega_bar_min=0.0,
        omega_bar_max=14.0,
        primary_scan_points=3501,
        verification_scan_points=7001,
        primary_phases=(0.0, 0.5),
        verification_phases=(0.25,),
        root_xtol_bar=1.0e-13,
        cluster_atol_bar=6.0e-3,
        cluster_rtol=1.0e-12,
        post_guard_tail_bar=1.0,
    )
    inventory = pilot.seed_free_root_inventory(
        _synthetic_close_root_matrix,
        1.0,
        policy,
        case_id="synthetic_close_roots",
        builder_id="prescribed_close_root_6x6_matrix",
    )

    assert inventory.status == "PASS"
    assert inventory.independent_agreement is True
    assert inventory.unresolved_low_sigma_count == 0
    assert len(inventory.slots) == 13
    assert_allclose(
        [slot.event.omega_bar for slot in inventory.slots],
        _CLOSE_ROOT_VALUES,
        rtol=0.0,
        atol=2.0e-9,
    )
    clustered_events = [
        event
        for event in inventory.primary.events
        if event.cluster_semantics == "NEAR_DEGENERATE_CLUSTER"
    ]
    assert len(clustered_events) == 2
    assert {event.multiplicity for event in clustered_events} == {1}
    assert {event.detected_nullity for event in clustered_events} == {1}
    assert {event.cluster_multiplicity for event in clustered_events} == {2}
    assert {event.cluster_total_nullity for event in clustered_events} == {2}
    assert len({event.cluster_id for event in clustered_events}) == 1
    expected_center = float(math.fsum(_CLOSE_ROOT_VALUES[4:6]) / 2.0)
    assert_allclose(
        [event.cluster_center_omega_bar for event in clustered_events],
        [expected_center, expected_center],
        rtol=0.0,
        atol=2.0e-9,
    )
