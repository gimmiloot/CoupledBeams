from dataclasses import dataclass

import numpy as np


@dataclass
class BeamParams:
    E: float
    rho: float
    r: float
    L_total: float  # fixed: L1+L2

    @property
    def S(self) -> float:
        return np.pi * self.r**2

    @property
    def I(self) -> float:
        return np.pi * self.r**4 / 4.0

    @property
    def L_base(self) -> float:
        return 0.5 * self.L_total

    @property
    def eps(self) -> float:
        return np.sqrt(self.I / self.S) / self.L_base

    @property
    def f_scale(self) -> float:
        return np.sqrt(self.E * self.I / (self.rho * self.S)) / (2.0 * np.pi)


def assemble_clamped_coupled_matrix(
    Lambda: float,
    beta: float,
    mu: float,
    eps: float,
) -> np.ndarray:
    L1f, L2f = 1.0 - mu, 1.0 + mu
    x1, x2 = Lambda * L1f, Lambda * L2f

    Cd1 = np.cos(x1) - np.cosh(x1)
    Sd1 = np.sin(x1) - np.sinh(x1)
    Cs1 = np.cos(x1) + np.cosh(x1)
    Ss1 = np.sin(x1) + np.sinh(x1)

    Cd2 = np.cos(x2) - np.cosh(x2)
    Sd2 = np.sin(x2) - np.sinh(x2)
    Cs2 = np.cos(x2) + np.cosh(x2)
    Ss2 = np.sin(x2) + np.sinh(x2)

    th1 = eps * (Lambda**2) * L1f
    th2 = eps * (Lambda**2) * L2f

    cb, sb = np.cos(beta), np.sin(beta)

    return np.array(
        [
            [Cd1, Sd1, -Cd2 * cb, Sd2 * cb, 0.0, np.sin(th2) * sb],
            [0.0, 0.0, Cd2 * sb, -Sd2 * sb, np.sin(th1), np.sin(th2) * cb],
            [-Ss1, Cd1, -Ss2, -Cd2, 0.0, 0.0],
            [-Cs1, -Ss1, Cs2, -Ss2, 0.0, 0.0],
            [
                -eps * Lambda * Sd1,
                eps * Lambda * Cs1,
                -eps * Lambda * Sd2 * cb,
                -eps * Lambda * Cs2 * cb,
                0.0,
                -np.cos(th2) * sb,
            ],
            [0.0, 0.0, eps * Lambda * Sd2 * sb, eps * Lambda * Cs2 * sb, np.cos(th1), -np.cos(th2) * cb],
        ],
        dtype=float,
    )


def det_clamped_coupled(Lambda: float, beta: float, mu: float, eps: float) -> float:
    return float(np.linalg.det(assemble_clamped_coupled_matrix(Lambda, beta, mu, eps)))


def frequency_scale(params: BeamParams) -> float:
    return params.f_scale / (params.L_base**2)


def lambdas_to_frequencies(lambdas: np.ndarray, params: BeamParams) -> np.ndarray:
    return np.asarray(lambdas, dtype=float) ** 2 * frequency_scale(params)


def segment_lengths(params: BeamParams, mu: float) -> tuple[float, float]:
    Lbase = params.L_base
    return Lbase * (1.0 - mu), Lbase * (1.0 + mu)


__all__ = [
    "BeamParams",
    "assemble_clamped_coupled_matrix",
    "det_clamped_coupled",
    "frequency_scale",
    "lambdas_to_frequencies",
    "segment_lengths",
]
