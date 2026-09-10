"""Physical EB adapter for one massless in-plane rotational spring.

Axes point from the clamps to the joint; beta_rad is in radians.
Theory: docs/laminated_beams/inplane_rotational_spring_joint.md.
This module does not implement a Reddy/FSDT spring or change rigid defaults.
"""
from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.linalg import block_diag, expm

from scripts.lib.reddy_inplane_geometry import STATE_ORDER
from scripts.lib.reddy_symmetric_coupled_beams import (
    joint_matrix as rigid_joint_matrix,
    positively_equilibrate_matrix,
)


def _finite(value: float, name: str, *, positive: bool = False) -> float:
    value = float(value)
    if not np.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


@dataclass(frozen=True)
class EBArm:
    """A = axial rigidity, D = flexural rigidity, m = mass/length."""

    A: float
    D: float
    m: float
    L: float

    def __post_init__(self) -> None:
        for name in ("A", "D", "m", "L"):
            object.__setattr__(self, name, _finite(getattr(self, name), name, positive=True))


@dataclass(frozen=True)
class Joint:
    """RIGID is an exact constraint; SPRING requires a finite stiffness."""

    mode: Literal["SPRING", "RIGID"]
    k_theta: float | None = None

    def __post_init__(self) -> None:
        if self.mode == "RIGID":
            if self.k_theta is not None:
                raise ValueError("RIGID requires k_theta=None")
        elif self.mode == "SPRING":
            if self.k_theta is None:
                raise ValueError("SPRING requires explicit k_theta, including zero")
            stiffness = _finite(self.k_theta, "k_theta")
            if stiffness < 0:
                raise ValueError("k_theta must be nonnegative")
            object.__setattr__(self, "k_theta", stiffness)
        else:
            raise ValueError("mode must be SPRING or RIGID")


def joint_matrix(beta_rad: float, joint: Joint) -> np.ndarray:
    """Physical 6x12 matrix; only row 3 differs from the old rigid joint."""
    beta_rad = _finite(beta_rad, "beta_rad")
    matrix = rigid_joint_matrix(beta_rad).copy()
    if joint.mode == "SPRING":
        matrix[2] = 0.0
        matrix[2, [2, 5, 8]] = [joint.k_theta, 1.0, -joint.k_theta]
    return matrix


def scalar_joint_residuals(states: np.ndarray, beta_rad: float, joint: Joint) -> np.ndarray:
    """Six physical residuals, written independently of the matrix assembly."""
    beta_rad = _finite(beta_rad, "beta_rad")
    c, s = np.cos(beta_rad), np.sin(beta_rad)
    u1, w1, p1, n1, q1, m1, u2, w2, p2, n2, q2, m2 = np.asarray(states)
    law = p1 - p2 if joint.mode == "RIGID" else m1 + joint.k_theta * (p1 - p2)
    return np.array([u1 + c*u2 + s*w2, w1 - s*u2 + c*w2, law,
                     n1 - c*n2 - s*q2, q1 + s*n2 - c*q2, m1 + m2])


def state_matrix(omega: float, arm: EBArm) -> np.ndarray:
    omega = _finite(omega, "omega")
    if omega < 0:
        raise ValueError("omega must be nonnegative")
    matrix = np.zeros((6, 6))
    matrix[0, 3] = 1 / arm.A
    matrix[1, 2] = -1
    matrix[2, 5] = 1 / arm.D
    matrix[3, 0] = matrix[4, 1] = -arm.m * omega**2
    matrix[5, 4] = 1
    return matrix


def state_scale(arm: EBArm) -> np.ndarray:
    """Same dimensional scales as the project's physical-state beam adapter."""
    return np.array([arm.L, arm.L, 1, arm.A, arm.D/arm.L**2, arm.D/arm.L])


def _scaled_transfer(omega: float, arm: EBArm) -> np.ndarray:
    scale = state_scale(arm)
    matrix = state_matrix(omega, arm) * scale[None, :] / scale[:, None]
    return expm(matrix * arm.L)


def transfer_matrix(omega: float, arm: EBArm) -> np.ndarray:
    scale = state_scale(arm)
    return scale[:, None] * _scaled_transfer(omega, arm) / scale[None, :]


def clamp_to_joint_map(omega: float, arm: EBArm) -> np.ndarray:
    """Maps the three physical clamp reactions (N,Q,M) to the endpoint state."""
    return transfer_matrix(omega, arm)[:, 3:]


@dataclass(frozen=True)
class BoundaryAssembly:
    physical: np.ndarray
    dimensionless: np.ndarray
    endpoint_map: np.ndarray
    reaction_scales: np.ndarray
    row_units: np.ndarray
    row_factors: np.ndarray


def boundary_assembly(omega: float, arm1: EBArm, arm2: EBArm,
                      beta_rad: float, joint: Joint, reference: EBArm) -> BoundaryAssembly:
    """B_dimless = diag(row_factors) B_physical diag(reaction_scales).

    reference.L and reference.D set l and D_ref. SPRING row 3 first uses
    the moment unit D_ref/l, then divides by max(1,kappa_theta).
    There is no division by k_theta and no arithmetic with infinity.
    """
    endpoint_blocks = [state_scale(arm)[:, None] * _scaled_transfer(omega, arm)[:, 3:]
                       for arm in (arm1, arm2)]
    endpoint_map = block_diag(*endpoint_blocks)
    reactions = np.concatenate([state_scale(arm)[3:] for arm in (arm1, arm2)])
    moment = reference.D / reference.L
    force = reference.D / reference.L**2
    units = np.array([reference.L, reference.L,
                      1.0 if joint.mode == "RIGID" else moment, force, force, moment])
    factors = 1 / units
    if joint.mode == "SPRING":
        kappa = joint.k_theta * reference.L / reference.D
        factors[2] /= max(1.0, kappa)
    reacted = joint_matrix(beta_rad, joint) @ endpoint_map
    return BoundaryAssembly(reacted / reactions[None, :], factors[:, None] * reacted,
                            endpoint_map, reactions, units, factors)


def endpoint_diagnostics(assembly: BoundaryAssembly, beta_rad: float, joint: Joint,
                         reference: EBArm, rank_rtol: float = 1e-12) -> dict:
    """Recover every detected null vector, not a mode shape or branch label.

    The helper's column factors act on dimensionless clamp reactions. They
    must be applied before conversion to physical reactions and endpoint states.
    Physical residuals use fixed moment/force units, without the kappa divisor.
    """
    eq = positively_equilibrate_matrix(assembly.dimensionless)
    _, singular, vh = np.linalg.svd(eq.scaled_matrix)
    ratios = singular / singular[0]
    nullity = int(np.count_nonzero(ratios <= rank_rtol))
    unit_state = np.tile([reference.L, reference.L, 1, reference.D/reference.L**2,
                          reference.D/reference.L**2, reference.D/reference.L], 2)
    vectors = []
    for vector in vh[-max(1, nullity):]:
        reactions_hat = eq.column_factors * vector
        states = assembly.endpoint_map @ reactions_hat
        amplitude = np.max(np.abs(states / unit_state))
        if not np.isfinite(amplitude) or amplitude <= 0:
            raise ValueError("invalid physical endpoint amplitude")
        states = states / amplitude
        residuals = scalar_joint_residuals(states, beta_rad, joint)
        normalized = residuals / assembly.row_units
        norm = np.linalg.norm(assembly.dimensionless, ord="fro") * np.linalg.norm(reactions_hat)
        vectors.append({
            "endpoint_states": states.tolist(),
            "physical_clamp_reactions": (assembly.reaction_scales*reactions_hat/amplitude).tolist(),
            "physical_residuals": residuals.tolist(),
            "normalized_physical_residuals": normalized.tolist(),
            "boundary_residual": float(np.linalg.norm(assembly.dimensionless@reactions_hat)/norm),
            "scaled_residual": float(np.linalg.norm(eq.scaled_matrix@vector) /
                                     (np.linalg.norm(eq.scaled_matrix, ord="fro")*np.linalg.norm(vector))),
            "delta_psi": float(states[2] - states[8]),
            "moments_in_reference_units": (states[[5, 11]]/(reference.D/reference.L)).tolist(),
        })
    return {"sigma_ratio": float(ratios[-1]), "nullity": nullity, "vectors": vectors}
