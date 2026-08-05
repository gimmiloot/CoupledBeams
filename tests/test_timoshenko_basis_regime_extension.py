from __future__ import annotations

import math

import numpy as np
import pytest
from scipy.linalg import expm

from scripts.lib import variable_length_timoshenko as timo


def _state_matrix(omega: float, section: timo.Section) -> np.ndarray:
    """Independent first-order form; it does not call the closed-form basis."""
    k_stiff = section.shear_stiffness
    b_stiff = section.bending_stiffness
    m = section.mass_per_length
    j = section.rotary_inertia_per_length
    return np.array(
        [
            [0.0, 1.0, 1.0 / k_stiff, 0.0],
            [0.0, 0.0, 0.0, 1.0 / b_stiff],
            [-m * omega**2, 0.0, 0.0, 0.0],
            [0.0, -j * omega**2, -1.0, 0.0],
        ],
        dtype=float,
    )


def _state_space_bending_columns(
    x: float,
    omega: float,
    section: timo.Section,
    scales: np.ndarray,
) -> dict[str, np.ndarray]:
    initial = np.array(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, section.shear_stiffness * scales[1]],
            [section.bending_stiffness * scales[0], 0.0],
        ],
        dtype=float,
    )
    state = expm(_state_matrix(omega, section) * float(x)) @ initial
    w, psi, shear, moment = state
    return {
        "w": w,
        "psi": psi,
        "w_prime": psi + shear / section.shear_stiffness,
        "psi_prime": moment / section.bending_stiffness,
    }


def _old_mixed_columns(x: float, basis: timo.TimoshenkoBasis) -> dict[str, np.ndarray]:
    assert basis.regime == timo.TIMO_REGIME_MIXED
    a, b, h, q = basis.a, basis.b, basis.h, basis.q
    ax, bx = a * float(x), b * float(x)
    return {
        "w": np.array([np.cosh(ax) - np.cos(bx), np.sinh(ax) - (h / q) * np.sin(bx)]),
        "psi": np.array([h * np.sinh(ax) + q * np.sin(bx), h * (np.cosh(ax) - np.cos(bx))]),
        "w_prime": np.array(
            [a * np.sinh(ax) + b * np.sin(bx), a * np.cosh(ax) - (h / q) * b * np.cos(bx)]
        ),
        "psi_prime": np.array(
            [h * a * np.cosh(ax) + q * b * np.cos(bx), h * (a * np.sinh(ax) + b * np.sin(bx))]
        ),
    }


@pytest.mark.parametrize(
    "cutoff_factor,expected_regime",
    [
        (0.8, timo.TIMO_REGIME_MIXED),
        (1.0, timo.TIMO_REGIME_CUTOFF),
        (1.2, timo.TIMO_REGIME_TWO_TRIG),
    ],
)
def test_closed_form_columns_match_independent_state_space_oracle(
    cutoff_factor: float,
    expected_regime: str,
) -> None:
    epsilon = 0.05
    section = timo.section_from_epsilon(epsilon)
    Lambda = cutoff_factor * timo.lambda_cutoff(epsilon, section)
    basis = timo.timo_basis(Lambda, epsilon, section)
    assert basis.regime == expected_regime
    omega = timo.project_omega(Lambda, epsilon)
    for x in (-0.75, -0.1, 0.0, 0.4, 1.0):
        closed = timo.bending_endpoint_columns(x, basis)
        oracle = _state_space_bending_columns(x, omega, section, timo.bending_column_scales(basis))
        for field in closed:
            assert np.allclose(closed[field], oracle[field], rtol=2.0e-11, atol=2.0e-11)


def test_clamped_reduced_columns_have_canonical_initial_derivatives_and_rank() -> None:
    epsilon = 0.05
    section = timo.section_from_epsilon(epsilon)
    cutoff = timo.lambda_cutoff(epsilon, section)
    for Lambda in (0.7 * cutoff, cutoff, 1.3 * cutoff):
        columns = timo.bending_endpoint_columns(0.0, timo.timo_basis(Lambda, epsilon, section))
        scales = timo.bending_column_scales(timo.timo_basis(Lambda, epsilon, section))
        assert np.array_equal(columns["w"], np.zeros(2))
        assert np.array_equal(columns["psi"], np.zeros(2))
        assert np.allclose(columns["w_prime"], (0.0, scales[1]), rtol=0.0, atol=2.0e-11)
        assert np.allclose(columns["psi_prime"], (scales[0], 0.0), rtol=0.0, atol=2.0e-11)
        initial_derivatives = np.vstack((columns["w_prime"], columns["psi_prime"]))
        assert np.linalg.matrix_rank(initial_derivatives) == 2


def test_closed_form_basis_is_continuous_across_cutoff() -> None:
    epsilon = 0.05
    section = timo.section_from_epsilon(epsilon)
    cutoff = timo.lambda_cutoff(epsilon, section)
    exact = timo.bending_endpoint_columns(0.73, timo.timo_basis(cutoff, epsilon, section))
    for relative_offset in (1.0e-5, 1.0e-7):
        below = timo.bending_endpoint_columns(
            0.73, timo.timo_basis(cutoff * (1.0 - relative_offset), epsilon, section)
        )
        above = timo.bending_endpoint_columns(
            0.73, timo.timo_basis(cutoff * (1.0 + relative_offset), epsilon, section)
        )
        scale = relative_offset * 200.0
        for field in exact:
            assert np.allclose(below[field], exact[field], rtol=scale, atol=scale)
            assert np.allclose(above[field], exact[field], rtol=scale, atol=scale)


def test_characteristic_roots_change_sign_continuously_at_cutoff() -> None:
    epsilon = 0.05
    section = timo.section_from_epsilon(epsilon)
    cutoff = timo.lambda_cutoff(epsilon, section)
    below = timo.timo_basis(cutoff * (1.0 - 1.0e-7), epsilon, section)
    exact = timo.timo_basis(cutoff, epsilon, section)
    above = timo.timo_basis(cutoff * (1.0 + 1.0e-7), epsilon, section)
    assert below.z_a > 0.0
    assert exact.z_a == pytest.approx(0.0, abs=2.0e-12)
    assert above.z_a < 0.0
    assert below.z_b < 0.0 and exact.z_b < 0.0 and above.z_b < 0.0
    assert below.z_b == pytest.approx(exact.z_b, rel=5.0e-7)
    assert above.z_b == pytest.approx(exact.z_b, rel=5.0e-7)


@pytest.mark.parametrize("cutoff_factor", [0.8, 1.0, 1.2])
def test_w_psi_derivatives_and_physical_moment_shear_relations(cutoff_factor: float) -> None:
    epsilon = 0.05
    section = timo.section_from_epsilon(epsilon)
    Lambda = cutoff_factor * timo.lambda_cutoff(epsilon, section)
    basis = timo.timo_basis(Lambda, epsilon, section)
    x = 0.43
    step = 2.0e-6
    values = timo.bending_endpoint_columns(x, basis)
    left = timo.bending_endpoint_columns(x - step, basis)
    right = timo.bending_endpoint_columns(x + step, basis)
    numerical_w_prime = (right["w"] - left["w"]) / (2.0 * step)
    numerical_psi_prime = (right["psi"] - left["psi"]) / (2.0 * step)
    assert np.allclose(numerical_w_prime, values["w_prime"], rtol=2.0e-8, atol=2.0e-8)
    assert np.allclose(numerical_psi_prime, values["psi_prime"], rtol=2.0e-8, atol=2.0e-8)
    endpoint = timo.endpoint_columns(x, 0.7, basis, section)
    assert np.allclose(
        endpoint["Q"][:2],
        section.shear_stiffness * (values["w_prime"] - values["psi"]),
        rtol=2.0e-13,
        atol=2.0e-13,
    )
    assert np.allclose(
        endpoint["M"][:2],
        section.bending_stiffness * values["psi_prime"],
        rtol=2.0e-13,
        atol=2.0e-13,
    )


def test_coupled_matrix_supports_mixed_arm_regimes_and_both_arms_above_cutoff() -> None:
    epsilon, mu, eta = 0.05, 0.7, -0.5
    factors = timo.tau_factors(mu, eta)
    sections = (
        timo.section_from_epsilon_tau(epsilon, factors.tau1),
        timo.section_from_epsilon_tau(epsilon, factors.tau2),
    )
    cutoffs = tuple(timo.lambda_cutoff(epsilon, section) for section in sections)
    between = 0.5 * (min(cutoffs) + max(cutoffs))
    regimes = {timo.timo_basis(between, epsilon, section).regime for section in sections}
    assert regimes == {timo.TIMO_REGIME_MIXED, timo.TIMO_REGIME_TWO_TRIG}
    mixed_matrix, _ = timo.timo_coupling_matrix(between, 45.0, mu, epsilon, eta)
    assert np.all(np.isfinite(mixed_matrix))
    above = 1.1 * max(cutoffs)
    assert all(
        timo.timo_basis(above, epsilon, section).regime == timo.TIMO_REGIME_TWO_TRIG
        for section in sections
    )
    above_matrix, _ = timo.timo_coupling_matrix(above, 45.0, mu, epsilon, eta)
    assert np.all(np.isfinite(above_matrix))


def test_system_matrix_is_continuous_and_not_forced_singular_at_arm_cutoff() -> None:
    epsilon, mu, eta = 0.05, 0.3, 0.2
    factors = timo.tau_factors(mu, eta)
    section1 = timo.section_from_epsilon_tau(epsilon, factors.tau1)
    cutoff = timo.lambda_cutoff(epsilon, section1)
    matrices = [
        timo.timo_coupling_matrix(cutoff * factor, 37.0, mu, epsilon, eta)[0]
        for factor in (1.0 - 1.0e-7, 1.0, 1.0 + 1.0e-7)
    ]
    scale = max(1.0, float(np.max(np.abs(matrices[1]))))
    assert np.max(np.abs(matrices[0] - matrices[1])) / scale < 2.0e-5
    assert np.max(np.abs(matrices[2] - matrices[1])) / scale < 2.0e-5
    singular_values = np.linalg.svd(timo.row_normalized(matrices[1]), compute_uv=False)
    assert singular_values[-1] > 1.0e-8


def test_mixed_regime_is_only_an_invertible_column_rescaling_of_legacy_basis() -> None:
    epsilon = 0.02
    section = timo.section_from_epsilon(epsilon)
    basis = timo.timo_basis(8.0, epsilon, section)
    assert basis.regime == timo.TIMO_REGIME_MIXED
    scales = timo.bending_column_scales(basis)
    assert np.all(np.isfinite(scales))
    assert np.all(np.abs(scales) > 0.0)
    current = timo.bending_endpoint_columns(0.61, basis)
    legacy = _old_mixed_columns(0.61, basis)
    for field in current:
        assert np.allclose(current[field], legacy[field], rtol=3.0e-13, atol=3.0e-13)


def test_spatial_roots_satisfy_frozen_characteristic_polynomial() -> None:
    epsilon = 0.05
    section = timo.section_from_epsilon(epsilon)
    cutoff = timo.lambda_cutoff(epsilon, section)
    for Lambda in (0.8 * cutoff, cutoff, 1.2 * cutoff):
        basis = timo.timo_basis(Lambda, epsilon, section)
        omega = timo.project_omega(Lambda, epsilon)
        B = section.bending_stiffness
        K = section.shear_stiffness
        m = section.mass_per_length
        j = section.rotary_inertia_per_length
        c2 = omega**2 * (B * m / K + j)
        c0 = j * m * omega**4 / K - m * omega**2
        scale = max(abs(B * basis.z_b**2), abs(c2 * basis.z_b), abs(c0), 1.0)
        for z in (basis.z_a, basis.z_b):
            assert abs(B * z**2 + c2 * z + c0) <= 3.0e-15 * scale


def test_basis_evaluator_version_is_explicit() -> None:
    assert timo.TIMOSHENKO_BASIS_EVALUATOR_VERSION == "closed_form_regime_complete_v2"
    assert timo.TIMO_REGIME_CUTOFF == "cutoff_limit"
