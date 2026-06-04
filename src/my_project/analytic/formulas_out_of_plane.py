from __future__ import annotations

from dataclasses import dataclass
from math import cos, isfinite, pi, sin, sqrt

import numpy as np

from .formulas_thickness_mismatch import thickness_mismatch_factors


@dataclass(frozen=True)
class OutOfPlaneFactors:
    """Dimensionless factors for the out-of-plane EB + torsion determinant."""

    Lambda: float
    beta: float
    mu: float
    epsilon: float
    eta: float
    poisson: float
    L1: float
    L2: float
    tau1: float
    tau2: float
    alpha1: float
    alpha2: float
    gamma1: float
    gamma2: float
    a1: float
    a2: float
    b1: float
    b2: float
    c1: float
    c2: float
    e1: float
    e2: float
    chi_T: float
    sqrt_E_over_G: float


def _finite_float(name: str, value: float) -> float:
    value_f = float(value)
    if not isfinite(value_f):
        raise ValueError(f"{name} must be finite.")
    return value_f


def _positive_float(name: str, value: float) -> float:
    value_f = _finite_float(name, value)
    if value_f <= 0.0:
        raise ValueError(f"{name} must be positive.")
    return value_f


def out_of_plane_factors(
    Lambda: float,
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
) -> OutOfPlaneFactors:
    """Return factors for the out-of-plane determinant.

    ``beta`` is in radians. ``Lambda`` must be positive; ``Lambda=0`` is not
    treated as a physical root of this determinant.
    """

    Lambda_f = _positive_float("Lambda", Lambda)
    beta_f = _finite_float("beta", beta)
    epsilon_f = _positive_float("epsilon", epsilon)
    poisson_f = _finite_float("poisson", poisson)
    if poisson_f <= -1.0:
        raise ValueError("poisson must be greater than -1.")

    thickness = thickness_mismatch_factors(mu, eta)
    for name, value in (("tau1", thickness.tau1), ("tau2", thickness.tau2)):
        if not isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be positive and finite.")

    L1 = 1.0 - thickness.mu
    L2 = 1.0 + thickness.mu
    sqrt_E_over_G = sqrt(2.0 * (1.0 + poisson_f))
    sqrt_G_over_E = 1.0 / sqrt_E_over_G

    alpha1 = Lambda_f * L1 / sqrt(thickness.tau1)
    alpha2 = Lambda_f * L2 / sqrt(thickness.tau2)
    gamma1 = sqrt_E_over_G * epsilon_f * Lambda_f**2 * L1
    gamma2 = sqrt_E_over_G * epsilon_f * Lambda_f**2 * L2

    return OutOfPlaneFactors(
        Lambda=Lambda_f,
        beta=beta_f,
        mu=thickness.mu,
        epsilon=epsilon_f,
        eta=thickness.eta,
        poisson=poisson_f,
        L1=float(L1),
        L2=float(L2),
        tau1=thickness.tau1,
        tau2=thickness.tau2,
        alpha1=float(alpha1),
        alpha2=float(alpha2),
        gamma1=float(gamma1),
        gamma2=float(gamma2),
        a1=thickness.tau1 ** (-0.5),
        a2=thickness.tau2 ** (-0.5),
        b1=thickness.tau1**3,
        b2=thickness.tau2**3,
        c1=thickness.tau1 ** 2.5,
        c2=thickness.tau2 ** 2.5,
        e1=thickness.tau1**4,
        e2=thickness.tau2**4,
        chi_T=2.0 * epsilon_f * sqrt_G_over_E,
        sqrt_E_over_G=float(sqrt_E_over_G),
    )


def _bending_values(alpha: float) -> tuple[float, float, float, float]:
    c = np.cos(alpha)
    s = np.sin(alpha)
    ch = np.cosh(alpha)
    sh = np.sinh(alpha)
    c_minus = c - ch
    s_minus = s - sh
    c_plus = c + ch
    s_plus = s + sh
    return float(c_minus), float(s_minus), float(c_plus), float(s_plus)


def assemble_out_of_plane_matrix(
    Lambda: float,
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
) -> np.ndarray:
    """Assemble the 6x6 out-of-plane determinant matrix.

    The unknown order is ``X=(A1, B1, A2, B2, P1, P2)^T``. ``beta`` is in
    radians. The row order is the compact determinant order documented in
    ``docs/theory/out_of_plane_eb_torsion.md``:
    z, psi compatibility, phi compatibility, shear, torsion/moment, bending
    moment.
    """

    f = out_of_plane_factors(Lambda, beta, mu, epsilon, eta, poisson)

    cd1, sd1, cp1, sp1 = _bending_values(f.alpha1)
    cd2, sd2, cp2, sp2 = _bending_values(f.alpha2)
    sg1, cg1 = sin(f.gamma1), cos(f.gamma1)
    sg2, cg2 = sin(f.gamma2), cos(f.gamma2)
    cb, sb = cos(f.beta), sin(f.beta)

    matrix = np.zeros((6, 6), dtype=float)

    # Z1 - Z2 = 0
    matrix[0, 0] = cd1
    matrix[0, 1] = sd1
    matrix[0, 2] = -cd2
    matrix[0, 3] = -sd2

    # Psi1 + Psi2 cos(beta) - Lambda a2 R2 sin(beta) = 0
    matrix[1, 2] = f.Lambda * f.a2 * sp2 * sb
    matrix[1, 3] = -f.Lambda * f.a2 * cd2 * sb
    matrix[1, 4] = sg1
    matrix[1, 5] = sg2 * cb

    # -Lambda a1 R1 - Psi2 sin(beta) - Lambda a2 R2 cos(beta) = 0
    matrix[2, 0] = f.Lambda * f.a1 * sp1
    matrix[2, 1] = -f.Lambda * f.a1 * cd1
    matrix[2, 2] = f.Lambda * f.a2 * sp2 * cb
    matrix[2, 3] = -f.Lambda * f.a2 * cd2 * cb
    matrix[2, 5] = -sg2 * sb

    # c1 V1 + c2 V2 = 0
    matrix[3, 0] = f.c1 * sd1
    matrix[3, 1] = -f.c1 * cp1
    matrix[3, 2] = f.c2 * sd2
    matrix[3, 3] = -f.c2 * cp2

    # chi_T e1 D1 - chi_T e2 D2 cos(beta) + b2 K2 sin(beta) = 0
    matrix[4, 2] = -f.b2 * cp2 * sb
    matrix[4, 3] = -f.b2 * sp2 * sb
    matrix[4, 4] = f.chi_T * f.e1 * cg1
    matrix[4, 5] = -f.chi_T * f.e2 * cg2 * cb

    # -b1 K1 + b2 K2 cos(beta) + chi_T e2 D2 sin(beta) = 0
    matrix[5, 0] = f.b1 * cp1
    matrix[5, 1] = f.b1 * sp1
    matrix[5, 2] = -f.b2 * cp2 * cb
    matrix[5, 3] = -f.b2 * sp2 * cb
    matrix[5, 5] = f.chi_T * f.e2 * cg2 * sb

    return matrix


def det_out_of_plane(
    Lambda: float,
    beta: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
) -> float:
    """Return ``det M_perp`` for the out-of-plane determinant."""

    matrix = assemble_out_of_plane_matrix(Lambda, beta, mu, epsilon, eta, poisson)
    return float(np.linalg.det(matrix))


def assemble_out_of_plane_beta0_blocks(
    Lambda: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    poisson: float = 0.3,
) -> tuple[np.ndarray, np.ndarray]:
    """Return the beta=0 bending and torsion blocks.

    The returned blocks use bending rows ``[0, 2, 3, 5]`` and torsion rows
    ``[1, 4]`` from the full determinant matrix.
    """

    matrix = assemble_out_of_plane_matrix(
        Lambda=Lambda,
        beta=0.0,
        mu=mu,
        epsilon=epsilon,
        eta=eta,
        poisson=poisson,
    )
    bending = matrix[np.ix_([0, 2, 3, 5], [0, 1, 2, 3])]
    torsion = matrix[np.ix_([1, 4], [4, 5])]
    return bending, torsion


def torsion_roots_uniform_eta0_beta0(n: int, epsilon: float, poisson: float = 0.3) -> float:
    """Return the eta=0, beta=0 torsion root for mode number ``n``."""

    n_i = int(n)
    if n_i <= 0:
        raise ValueError("n must be a positive integer.")
    epsilon_f = _positive_float("epsilon", epsilon)
    poisson_f = _finite_float("poisson", poisson)
    if poisson_f <= -1.0:
        raise ValueError("poisson must be greater than -1.")
    sqrt_E_over_G = sqrt(2.0 * (1.0 + poisson_f))
    return float(sqrt(n_i * pi / (2.0 * epsilon_f * sqrt_E_over_G)))


def fixed_fixed_bending_roots_total_length_two() -> np.ndarray:
    """First three eta=0, beta=0 bending roots for total length two."""

    clamped_clamped_alphas = np.array(
        [
            4.730040744862704,
            7.853204624095838,
            10.99560783800167,
        ],
        dtype=float,
    )
    return clamped_clamped_alphas / 2.0


__all__ = [
    "OutOfPlaneFactors",
    "assemble_out_of_plane_beta0_blocks",
    "assemble_out_of_plane_matrix",
    "det_out_of_plane",
    "fixed_fixed_bending_roots_total_length_two",
    "out_of_plane_factors",
    "torsion_roots_uniform_eta0_beta0",
]
