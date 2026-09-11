"""Symmetric RLB spring assembly and a coefficient-only RLB -> EB limit.

The reduced section is immutable. epsilon_limit is a diagnostic coefficient
parameter, not a ply thickness or a claim of laminate realizability.
Physical state: [u,w,psi,N,Q,M]; local axes run from clamps to joint.
"""
from dataclasses import dataclass

import numpy as np
from scipy.linalg import block_diag, expm

from scripts.lib import reddy_symmetric_laminated_beam as native
from scripts.lib.inplane_rotational_spring_eb import BoundaryAssembly, EBArm, Joint, joint_matrix


def benchmark_section():
    """Four actual equal 0-degree plies; K=5/6 is this benchmark's input."""
    shear = 1 / (2 * (1 + .3))
    material = native.OrthotropicLamina(1., 1., .3, shear, shear, shear, 1.)
    section = native.integrate_laminate(tuple(native.Ply(material, 0., .05/4) for _ in range(4)))
    properties = native.reduce_to_beam_properties(section, width=.20, K=5/6)
    return section, properties


@dataclass(frozen=True)
class LimitArm:
    properties: native.BeamProperties
    L: float
    epsilon_limit: float

    def __post_init__(self):
        if not isinstance(self.properties, native.BeamProperties):
            raise TypeError("properties must be a reduced BeamProperties section")
        if not np.isfinite(self.L) or self.L <= 0:
            raise ValueError("L must be finite and positive")
        if not np.isfinite(self.epsilon_limit) or not 0 <= self.epsilon_limit <= 1:
            raise ValueError("epsilon_limit must be finite in [0,1]")

    @property
    def invS(self):
        return self.epsilon_limit / self.properties.S

    @property
    def J(self):
        return self.epsilon_limit * self.properties.J


def state_matrix(omega: float, arm: LimitArm) -> np.ndarray:
    """Use native RLB physics, modifying exactly the two limit coefficients.

epsilon_limit=0 follows this same construction, without an EB dispatch,
infinite S, or modification of the constitutive API.
"""
    matrix = native.combined_state_matrix(omega, arm.properties).copy()
    matrix[1, 4] = arm.invS
    matrix[5, 2] = -arm.J * omega**2
    return matrix


def state_scale(arm: LimitArm) -> np.ndarray:
    return np.diag(native.combined_state_scale(arm.properties, arm.L))


def scaled_transfer(omega: float, arm: LimitArm) -> np.ndarray:
    scale = state_scale(arm)
    return expm(state_matrix(omega, arm) * scale[None, :] / scale[:, None] * arm.L)


def transfer_matrix(omega: float, arm: LimitArm) -> np.ndarray:
    scale = state_scale(arm)
    return scale[:, None] * scaled_transfer(omega, arm) / scale[None, :]


def clamp_to_joint_map(omega: float, arm: LimitArm) -> np.ndarray:
    return transfer_matrix(omega, arm)[:, 3:]


def boundary_assembly(omega: float, arm1: LimitArm, arm2: LimitArm,
                      beta_rad: float, joint: Joint, reference: EBArm) -> BoundaryAssembly:
    """RLB endpoint maps, common physical spring and fixed reference units.

The returned maps can feed EB's *generic* endpoint_diagnostics: that helper
uses only the supplied matrices and never constructs an EB transfer.
"""
    block1 = state_scale(arm1)[:, None] * scaled_transfer(omega, arm1)[:, 3:]
    # Identity is intentional: no array-valued dataclass equality is needed.
    block2 = block1 if arm2 is arm1 else state_scale(arm2)[:, None] * scaled_transfer(omega, arm2)[:, 3:]
    endpoints = block_diag(block1, block2)
    reactions = np.concatenate([state_scale(arm)[3:] for arm in (arm1, arm2)])
    moment, force = reference.D/reference.L, reference.D/reference.L**2
    units = np.array([reference.L, reference.L,
                      1. if joint.mode == "RIGID" else moment, force, force, moment])
    factors = 1 / units
    if joint.mode == "SPRING":
        factors[2] /= max(1., joint.k_theta * reference.L/reference.D)
    reacted = joint_matrix(beta_rad, joint) @ endpoints
    return BoundaryAssembly(reacted/reactions[None, :], factors[:, None]*reacted,
                            endpoints, reactions, units, factors)
