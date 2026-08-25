"""Physical three-dimensional coordinate contract for the RLB in-plane model.

The local coordinate of each arm points from its outer clamp to the future
connection point.  This module defines geometry and physical vector mappings
only; it contains no connection equations, characteristic determinants, or
frequency calculations.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable

import numpy as np
from numpy.typing import NDArray


Vector3 = NDArray[np.float64]
Matrix3 = NDArray[np.float64]

_GEOMETRY_ATOL = 5.0e-14


def _readonly_vector3(values: Iterable[float], *, name: str) -> Vector3:
    vector = np.asarray(tuple(values), dtype=float)
    if vector.shape != (3,):
        raise ValueError(f"{name} must contain exactly three components")
    if not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must contain only finite components")
    result = np.array(vector, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _readonly_matrix3(columns: tuple[Vector3, Vector3, Vector3]) -> Matrix3:
    matrix = np.column_stack(columns)
    matrix.setflags(write=False)
    return matrix


def _finite_scalar(value: float, *, name: str) -> float:
    scalar = float(value)
    if not math.isfinite(scalar):
        raise ValueError(f"{name} must be finite")
    return scalar


E_X = _readonly_vector3((1.0, 0.0, 0.0), name="E_X")
E_Y = _readonly_vector3((0.0, 1.0, 0.0), name="E_Y")
E_Z = _readonly_vector3((0.0, 0.0, 1.0), name="E_Z")
k = _readonly_vector3(-E_Z, name="k")

STATE_ORDER = ("u", "w", "psi", "N", "Q", "M")


@dataclass(frozen=True, eq=False)
class ArmBasis:
    """Physical RLB basis of one arm in the fixed global view basis."""

    t: Vector3
    n: Vector3

    def __post_init__(self) -> None:
        t = _readonly_vector3(self.t, name="t")
        n = _readonly_vector3(self.n, name="n")
        object.__setattr__(self, "t", t)
        object.__setattr__(self, "n", n)

        if not math.isclose(
            float(np.linalg.norm(t)), 1.0, rel_tol=0.0, abs_tol=_GEOMETRY_ATOL
        ):
            raise ValueError("t must be a unit vector")
        if not math.isclose(
            float(np.linalg.norm(n)), 1.0, rel_tol=0.0, abs_tol=_GEOMETRY_ATOL
        ):
            raise ValueError("n must be a unit vector")
        if abs(float(np.dot(t, n))) > _GEOMETRY_ATOL:
            raise ValueError("t and n must be orthogonal")
        if not np.allclose(n, np.cross(k, t), rtol=0.0, atol=_GEOMETRY_ATOL):
            raise ValueError("n must equal k cross t")

    @property
    def e_x(self) -> Vector3:
        return self.t

    @property
    def e_y(self) -> Vector3:
        return E_Z

    @property
    def e_z(self) -> Vector3:
        return self.n

    @property
    def inplane_matrix(self) -> NDArray[np.float64]:
        """Return the column matrix ``[t, n]``."""

        matrix = np.column_stack((self.t, self.n))
        matrix.setflags(write=False)
        return matrix

    @property
    def reddy_matrix(self) -> Matrix3:
        """Return the right-handed column matrix ``[e_x, e_y, e_z]``."""

        return _readonly_matrix3((self.e_x, self.e_y, self.e_z))


@dataclass(frozen=True, eq=False)
class ReddyInPlaneGeometry:
    """The two physical arm bases and the second outer-clamp ray."""

    beta_deg: float
    g2: Vector3
    arm1: ArmBasis
    arm2: ArmBasis

    def __post_init__(self) -> None:
        beta_deg = _finite_scalar(self.beta_deg, name="beta_deg")
        g2 = _readonly_vector3(self.g2, name="g2")
        object.__setattr__(self, "beta_deg", beta_deg)
        object.__setattr__(self, "g2", g2)

        beta_rad = math.radians(beta_deg)
        cosine = math.cos(beta_rad)
        sine = math.sin(beta_rad)
        expected_t1 = E_X
        expected_n1 = -E_Y
        expected_g2 = cosine * E_X + sine * E_Y
        expected_t2 = -cosine * expected_t1 + sine * expected_n1
        expected_n2 = -sine * expected_t1 - cosine * expected_n1
        vectors = (
            (self.arm1.t, expected_t1, "arm1.t"),
            (self.arm1.n, expected_n1, "arm1.n"),
            (g2, expected_g2, "g2"),
            (self.arm2.t, expected_t2, "arm2.t"),
            (self.arm2.n, expected_n2, "arm2.n"),
        )
        for actual, expected, name in vectors:
            if not np.allclose(actual, expected, rtol=0.0, atol=_GEOMETRY_ATOL):
                raise ValueError(f"{name} is inconsistent with beta_deg")

    @property
    def t1(self) -> Vector3:
        return self.arm1.t

    @property
    def n1(self) -> Vector3:
        return self.arm1.n

    @property
    def t2(self) -> Vector3:
        return self.arm2.t

    @property
    def n2(self) -> Vector3:
        return self.arm2.n


_ARM1 = ArmBasis(t=E_X, n=-E_Y)


def arm1_basis() -> ArmBasis:
    """Return the first-arm basis, with local x from outer clamp to connection."""

    return _ARM1


def reddy_inplane_geometry(beta_deg: float) -> ReddyInPlaneGeometry:
    """Construct the requested physical RLB geometry for an angle in degrees."""

    beta_value = _finite_scalar(beta_deg, name="beta_deg")
    beta_rad = math.radians(beta_value)
    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)

    arm1 = arm1_basis()
    g2 = cosine * E_X + sine * E_Y
    t2 = -cosine * arm1.t + sine * arm1.n
    n2 = -sine * arm1.t - cosine * arm1.n
    arm2 = ArmBasis(t=t2, n=n2)
    return ReddyInPlaneGeometry(
        beta_deg=beta_value,
        g2=g2,
        arm1=arm1,
        arm2=arm2,
    )


def arm2_basis(beta_deg: float) -> ArmBasis:
    """Return the second-arm physical basis, directed outer clamp to connection."""

    return reddy_inplane_geometry(beta_deg).arm2


def inplane_vector(basis: ArmBasis, axial: float, transverse: float) -> Vector3:
    """Map two local components to a physical global vector."""

    axial_value = _finite_scalar(axial, name="axial")
    transverse_value = _finite_scalar(transverse, name="transverse")
    return _readonly_vector3(
        axial_value * basis.t + transverse_value * basis.n,
        name="inplane_vector",
    )


def reddy_width_vector(basis: ArmBasis, component: float) -> Vector3:
    """Map a component along the physical Reddy ``e_y=-k=E_Z`` axis."""

    component_value = _finite_scalar(component, name="component")
    return _readonly_vector3(component_value * basis.e_y, name="reddy_width_vector")


def displacement_vector(basis: ArmBasis, u: float, w: float) -> Vector3:
    """Return ``d = u t + w n``."""

    return inplane_vector(basis, u, w)


def rotation_vector(basis: ArmBasis, psi: float) -> Vector3:
    """Return ``vartheta = psi e_y = -psi k``."""

    return reddy_width_vector(basis, psi)


def force_vector(basis: ArmBasis, axial_force: float, shear_force: float) -> Vector3:
    """Return ``F = N t + Q n``."""

    return inplane_vector(basis, axial_force, shear_force)


def moment_vector(basis: ArmBasis, bending_moment: float) -> Vector3:
    """Return ``m = M e_y = -M k``."""

    return reddy_width_vector(basis, bending_moment)
