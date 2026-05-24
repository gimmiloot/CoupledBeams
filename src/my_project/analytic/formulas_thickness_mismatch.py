from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ThicknessMismatchFactors:
    mu: float
    eta: float
    denom: float
    tau1: float
    tau2: float

    @property
    def mass_factor(self) -> float:
        return (1.0 - self.mu) * self.tau1**2 + (1.0 + self.mu) * self.tau2**2


def thickness_mismatch_factors(mu: float, eta: float) -> ThicknessMismatchFactors:
    """Mass-preserving radius factors for the diagnostic thickness-mismatch model."""
    mu_f = float(mu)
    eta_f = float(eta)
    if not (-1.0 < mu_f < 1.0):
        raise ValueError("mu must lie inside (-1, 1) for positive segment lengths.")
    if not (-1.0 < eta_f < 1.0):
        raise ValueError("eta must lie inside (-1, 1) for positive radii.")

    denom_sq = 1.0 + 2.0 * mu_f * eta_f + eta_f**2
    if denom_sq <= 0.0:
        raise ValueError("mass-preserving radius denominator is not positive.")
    denom = float(np.sqrt(denom_sq))
    tau1 = (1.0 - eta_f) / denom
    tau2 = (1.0 + eta_f) / denom
    return ThicknessMismatchFactors(mu=mu_f, eta=eta_f, denom=denom, tau1=float(tau1), tau2=float(tau2))


def local_epsilons(epsilon: float, mu: float, eta: float) -> tuple[float, float]:
    factors = thickness_mismatch_factors(mu, eta)
    return float(epsilon) * factors.tau1, float(epsilon) * factors.tau2


def thickness_to_length_ratios(epsilon: float, mu: float, eta: float) -> tuple[float, float]:
    """Diameter-to-length ratios ``2*r_i/l_i`` for the eta diagnostic model."""
    factors = thickness_mismatch_factors(mu, eta)
    eps_f = float(epsilon)
    ratio1 = 4.0 * eps_f * factors.tau1 / (1.0 - factors.mu)
    ratio2 = 4.0 * eps_f * factors.tau2 / (1.0 + factors.mu)
    return float(ratio1), float(ratio2)


def thin_rod_validity(epsilon: float, mu: float, eta: float, limit: float = 0.1) -> bool:
    """Return True when both rods satisfy the diameter criterion."""
    ratio1, ratio2 = thickness_to_length_ratios(epsilon, mu, eta)
    limit_f = float(limit)
    return ratio1 <= limit_f and ratio2 <= limit_f


def assemble_clamped_coupled_matrix_eta(
    Lambda: float,
    beta: float,
    mu: float,
    eps: float,
    eta: float,
) -> np.ndarray:
    """Diagnostic 6x6 determinant matrix for mass-preserving radius mismatch.

    The unknown ordering is the same as the base analytic model:
    (A1, B1, A2, B2, P1, P2).  At eta=0 this matrix reduces directly to
    `assemble_clamped_coupled_matrix` from `formulas.py`.
    """
    Lambda_f = float(Lambda)
    beta_f = float(beta)
    eps_f = float(eps)
    factors = thickness_mismatch_factors(mu, eta)

    L1f, L2f = 1.0 - factors.mu, 1.0 + factors.mu
    x1 = Lambda_f * L1f / np.sqrt(factors.tau1)
    x2 = Lambda_f * L2f / np.sqrt(factors.tau2)

    Cd1 = np.cos(x1) - np.cosh(x1)
    Sd1 = np.sin(x1) - np.sinh(x1)
    Cs1 = np.cos(x1) + np.cosh(x1)
    Ss1 = np.sin(x1) + np.sinh(x1)

    Cd2 = np.cos(x2) - np.cosh(x2)
    Sd2 = np.sin(x2) - np.sinh(x2)
    Cs2 = np.cos(x2) + np.cosh(x2)
    Ss2 = np.sin(x2) + np.sinh(x2)

    th1 = eps_f * (Lambda_f**2) * L1f
    th2 = eps_f * (Lambda_f**2) * L2f

    cb, sb = np.cos(beta_f), np.sin(beta_f)

    rot1 = factors.tau1 ** (-0.5)
    rot2 = factors.tau2 ** (-0.5)
    mom1 = factors.tau1**3
    mom2 = factors.tau2**3
    shear1 = factors.tau1 ** 2.5
    shear2 = factors.tau2 ** 2.5
    axial1 = factors.tau1**2
    axial2 = factors.tau2**2

    return np.array(
        [
            [Cd1, Sd1, -Cd2 * cb, Sd2 * cb, 0.0, np.sin(th2) * sb],
            [0.0, 0.0, Cd2 * sb, -Sd2 * sb, np.sin(th1), np.sin(th2) * cb],
            [-rot1 * Ss1, rot1 * Cd1, -rot2 * Ss2, -rot2 * Cd2, 0.0, 0.0],
            [-mom1 * Cs1, -mom1 * Ss1, mom2 * Cs2, -mom2 * Ss2, 0.0, 0.0],
            [
                -eps_f * Lambda_f * shear1 * Sd1,
                eps_f * Lambda_f * shear1 * Cs1,
                -eps_f * Lambda_f * shear2 * Sd2 * cb,
                -eps_f * Lambda_f * shear2 * Cs2 * cb,
                0.0,
                -axial2 * np.cos(th2) * sb,
            ],
            [
                0.0,
                0.0,
                eps_f * Lambda_f * shear2 * Sd2 * sb,
                eps_f * Lambda_f * shear2 * Cs2 * sb,
                axial1 * np.cos(th1),
                -axial2 * np.cos(th2) * cb,
            ],
        ],
        dtype=float,
    )


def det_eta(Lambda: float, beta: float, mu: float, epsilon: float, eta: float) -> float:
    return float(np.linalg.det(assemble_clamped_coupled_matrix_eta(Lambda, beta, mu, epsilon, eta)))


def find_roots_scan_bisect_eta(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float,
    n_roots: int,
    Lmin: float,
    Lmax: float,
    scan_step: float,
    bisect_iters: int = 70,
) -> list[float]:
    grid = np.arange(float(Lmin), float(Lmax) + float(scan_step), float(scan_step), dtype=float)
    vals = np.array([det_eta(L, beta, mu, epsilon, eta) for L in grid], dtype=float)

    roots: list[float] = []

    def sgn(x: float) -> int:
        return 1 if x > 0 else (-1 if x < 0 else 0)

    for i in range(len(grid) - 1):
        a, b = float(grid[i]), float(grid[i + 1])
        fa, fb = float(vals[i]), float(vals[i + 1])
        sa, sb = sgn(fa), sgn(fb)

        if sa == 0:
            roots.append(a)
        elif sa * sb < 0:
            left, right = a, b
            fl = fa
            for _ in range(int(bisect_iters)):
                mid = 0.5 * (left + right)
                fm = det_eta(mid, beta, mu, epsilon, eta)
                sm = sgn(fm)
                if sm == 0:
                    left = right = mid
                    break
                if sgn(fl) * sm < 0:
                    right = mid
                else:
                    left = mid
                    fl = fm
            roots.append(0.5 * (left + right))

        if len(roots) >= int(n_roots):
            break

    roots = sorted(set(float(root) for root in roots))
    return roots[: int(n_roots)]


def find_first_n_roots_eta(
    beta: float,
    mu: float,
    epsilon: float,
    eta: float,
    n_roots: int,
    Lmin: float = 0.2,
    Lmax0: float = 35.0,
    scan_step: float = 0.02,
    grow_factor: float = 1.35,
    max_tries: int = 7,
) -> np.ndarray:
    Lmax = float(Lmax0)
    roots: list[float] = []
    for _ in range(int(max_tries)):
        roots = find_roots_scan_bisect_eta(
            beta,
            mu,
            epsilon,
            eta,
            int(n_roots),
            float(Lmin),
            Lmax,
            float(scan_step),
        )
        if len(roots) >= int(n_roots):
            break
        Lmax *= float(grow_factor)

    out = np.full(int(n_roots), np.nan, dtype=float)
    out[: len(roots)] = roots
    return out


__all__ = [
    "ThicknessMismatchFactors",
    "assemble_clamped_coupled_matrix_eta",
    "det_eta",
    "find_first_n_roots_eta",
    "find_roots_scan_bisect_eta",
    "local_epsilons",
    "thickness_mismatch_factors",
    "thickness_to_length_ratios",
    "thin_rod_validity",
]
