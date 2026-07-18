from __future__ import annotations

"""Shared local-field to display-geometry mapping for in-plane mode shapes.

The Euler--Bernoulli and Timoshenko determinants use opposite transverse-field
sign conventions.  Both mappings below describe the same displayed rods; they
only convert already reconstructed local fields into Cartesian display values.
"""

from dataclasses import dataclass
import math

import numpy as np


@dataclass(frozen=True)
class DisplayBasis:
    tangent: np.ndarray
    normal: np.ndarray


@dataclass(frozen=True)
class DisplayCurve:
    x_base: np.ndarray
    y_base: np.ndarray
    displacement_x: np.ndarray
    displacement_y: np.ndarray
    x_deformed: np.ndarray
    y_deformed: np.ndarray


def rod1_display_basis() -> DisplayBasis:
    return DisplayBasis(
        tangent=np.array([1.0, 0.0], dtype=float),
        normal=np.array([0.0, -1.0], dtype=float),
    )


def rod2_display_basis(beta_deg: float) -> DisplayBasis:
    beta_rad = math.radians(float(beta_deg))
    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)
    return DisplayBasis(
        tangent=np.array([cosine, sine], dtype=float),
        normal=np.array([sine, -cosine], dtype=float),
    )


def eb_rod1_display_basis() -> DisplayBasis:
    return DisplayBasis(
        tangent=np.array([1.0, 0.0], dtype=float),
        normal=np.array([0.0, 1.0], dtype=float),
    )


def eb_rod2_display_basis(beta_deg: float) -> DisplayBasis:
    beta_rad = math.radians(float(beta_deg))
    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)
    return DisplayBasis(
        tangent=np.array([cosine, sine], dtype=float),
        normal=np.array([-sine, cosine], dtype=float),
    )


def _validated_fields(u: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    u_array = np.asarray(u, dtype=float)
    w_array = np.asarray(w, dtype=float)
    if u_array.shape != w_array.shape:
        raise ValueError("u and w must have the same shape")
    return u_array, w_array


def _display_displacement(
    u: np.ndarray,
    w: np.ndarray,
    basis: DisplayBasis,
) -> tuple[np.ndarray, np.ndarray]:
    u_array, w_array = _validated_fields(u, w)
    displacement_x = u_array * basis.tangent[0] + w_array * basis.normal[0]
    displacement_y = u_array * basis.tangent[1] + w_array * basis.normal[1]
    return displacement_x, displacement_y


def rod1_local_displacement_to_display(u: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return _display_displacement(u, w, rod1_display_basis())


def rod2_local_displacement_to_display(
    u: np.ndarray,
    w: np.ndarray,
    *,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    return _display_displacement(u, w, rod2_display_basis(float(beta_deg)))


def eb_rod1_local_displacement_to_display(
    u: np.ndarray,
    w: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    return _display_displacement(u, w, eb_rod1_display_basis())


def eb_rod2_local_displacement_to_display(
    u: np.ndarray,
    w: np.ndarray,
    *,
    beta_deg: float,
) -> tuple[np.ndarray, np.ndarray]:
    return _display_displacement(u, w, eb_rod2_display_basis(float(beta_deg)))


def _rod1_fields_to_display(
    x1: np.ndarray,
    u1: np.ndarray,
    w1: np.ndarray,
    *,
    scale: float,
    basis: DisplayBasis,
) -> DisplayCurve:
    x_base = np.asarray(x1, dtype=float)
    u_array, w_array = _validated_fields(u1, w1)
    if x_base.shape != u_array.shape:
        raise ValueError("x1, u1, and w1 must have the same shape")
    y_base = np.zeros_like(x_base)
    displacement_x, displacement_y = _display_displacement(u_array, w_array, basis)
    return DisplayCurve(
        x_base=x_base,
        y_base=y_base,
        displacement_x=displacement_x,
        displacement_y=displacement_y,
        x_deformed=x_base + float(scale) * displacement_x,
        y_deformed=y_base + float(scale) * displacement_y,
    )


def _rod2_fields_to_display(
    x2_grid: np.ndarray,
    u2: np.ndarray,
    w2: np.ndarray,
    *,
    l2: float,
    x_joint: float,
    scale: float,
    basis: DisplayBasis,
) -> DisplayCurve:
    x2_array = np.asarray(x2_grid, dtype=float)
    u_array, w_array = _validated_fields(u2, w2)
    if x2_array.shape != u_array.shape:
        raise ValueError("x2_grid, u2, and w2 must have the same shape")
    xi2 = x2_array + float(l2)
    x_base = float(x_joint) + xi2 * basis.tangent[0]
    y_base = xi2 * basis.tangent[1]
    displacement_x, displacement_y = _display_displacement(u_array, w_array, basis)
    return DisplayCurve(
        x_base=x_base,
        y_base=y_base,
        displacement_x=displacement_x,
        displacement_y=displacement_y,
        x_deformed=x_base + float(scale) * displacement_x,
        y_deformed=y_base + float(scale) * displacement_y,
    )


def rod1_local_fields_to_display(
    x1: np.ndarray,
    u1: np.ndarray,
    w1: np.ndarray,
    *,
    scale: float,
) -> DisplayCurve:
    return _rod1_fields_to_display(
        x1,
        u1,
        w1,
        scale=scale,
        basis=rod1_display_basis(),
    )


def rod2_local_fields_to_display(
    x2_grid: np.ndarray,
    u2: np.ndarray,
    w2: np.ndarray,
    *,
    l2: float,
    x_joint: float,
    beta_deg: float,
    scale: float,
) -> DisplayCurve:
    return _rod2_fields_to_display(
        x2_grid,
        u2,
        w2,
        l2=l2,
        x_joint=x_joint,
        scale=scale,
        basis=rod2_display_basis(float(beta_deg)),
    )


def eb_rod1_local_fields_to_display(
    x1: np.ndarray,
    u1: np.ndarray,
    w1: np.ndarray,
    *,
    scale: float,
) -> DisplayCurve:
    return _rod1_fields_to_display(
        x1,
        u1,
        w1,
        scale=scale,
        basis=eb_rod1_display_basis(),
    )


def eb_rod2_local_fields_to_display(
    x2_grid: np.ndarray,
    u2: np.ndarray,
    w2: np.ndarray,
    *,
    l2: float,
    x_joint: float,
    beta_deg: float,
    scale: float,
) -> DisplayCurve:
    return _rod2_fields_to_display(
        x2_grid,
        u2,
        w2,
        l2=l2,
        x_joint=x_joint,
        scale=scale,
        basis=eb_rod2_display_basis(float(beta_deg)),
    )
