"""Narrow Chapter-2 helper for two rods with an ideal rigid angular joint.

The state of each arm is the accepted book-facing corrected state
``[w, psi_b, Phi_t, Q, M, M_T]``.  This module adds only joint geometry,
notation conversion, physical end maps, and the coupled/reference boundary
matrices needed by the small elastic pilot.  It is not a production API.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from scripts.lib.yartsev_ch2_monoclinic_rod import (
    RodPoint,
    cantilever_clamp_matrix,
    physical_state_transfer_matrix,
)


BOOK_STATE_SIZE = 6


@dataclass(frozen=True)
class JointBasis:
    """Project bases expressed in one fixed Cartesian frame."""

    e_z: NDArray[np.float64]
    t_1: NDArray[np.float64]
    n_1: NDArray[np.float64]
    t_2: NDArray[np.float64]
    n_2: NDArray[np.float64]


@dataclass(frozen=True)
class MatrixScaling:
    """Positive left/right diagonal factors used for equilibration."""

    row_factors: NDArray[np.float64]
    column_factors: NDArray[np.float64]


def notation_transform_matrix() -> NDArray[np.float64]:
    """Return ``P`` in ``y_old = P @ y_book``."""

    return np.diag([1.0, -1.0, 1.0, 1.0, -1.0, 1.0])


def two_arm_notation_transform_matrix() -> NDArray[np.float64]:
    """Return ``block_diag(P, P)`` without a SciPy dependency."""

    matrix = np.zeros((12, 12), dtype=float)
    matrix[:6, :6] = notation_transform_matrix()
    matrix[6:, 6:] = notation_transform_matrix()
    return matrix


def joint_basis(beta: float) -> JointBasis:
    """Return the unchanged project basis convention for angle ``beta``."""

    c = float(np.cos(beta))
    s = float(np.sin(beta))
    e_z = np.array([0.0, 0.0, 1.0])
    t_1 = np.array([1.0, 0.0, 0.0])
    n_1 = np.cross(e_z, t_1)
    t_2 = -c * t_1 + s * n_1
    n_2 = -s * t_1 - c * n_1
    return JointBasis(e_z=e_z, t_1=t_1, n_1=n_1, t_2=t_2, n_2=n_2)


def physical_rotation_vector(
    Phi_t: float, psi_b: float, tangent: NDArray[np.float64], normal: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Return the physical book rotation ``Phi t - psi n``."""

    return float(Phi_t) * tangent - float(psi_b) * normal


def physical_moment_vector(
    M_T: float, M: float, tangent: NDArray[np.float64], normal: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Return the physical book moment ``M_T t - M n``."""

    return float(M_T) * tangent - float(M) * normal


def rotation_components(
    vector: NDArray[np.float64], tangent: NDArray[np.float64], normal: NDArray[np.float64]
) -> tuple[float, float]:
    """Return book components ``(Phi, psi)`` of an in-plane vector."""

    value = np.asarray(vector, dtype=float)
    return float(value @ tangent), float(-(value @ normal))


def moment_components(
    vector: NDArray[np.float64], tangent: NDArray[np.float64], normal: NDArray[np.float64]
) -> tuple[float, float]:
    """Return book resultants ``(M_T, M)`` of an in-plane moment."""

    value = np.asarray(vector, dtype=float)
    return float(value @ tangent), float(-(value @ normal))


def joint_matrix_book(beta: float) -> NDArray[np.float64]:
    """Return the exact ``6 x 12`` rigid-joint matrix in book notation."""

    c = float(np.cos(beta))
    s = float(np.sin(beta))
    return np.array(
        [
            [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, -s, c, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, c, s, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, s, -c],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, -c, -s],
        ],
        dtype=float,
    )


def joint_matrix_old(beta: float) -> NDArray[np.float64]:
    """Return the old-project joint matrix with rows oriented for the sign gate.

    The third and sixth scalar conditions in the old theory note are multiplied
    by ``-1``.  They are the same homogeneous conditions, and this explicit row
    orientation makes ``J_book = J_old @ block_diag(P, P)`` exact.
    """

    c = float(np.cos(beta))
    s = float(np.sin(beta))
    return np.array(
        [
            [1, 0, 0, 0, 0, 0, -1, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, s, c, 0, 0, 0],
            [0, -1, 0, 0, 0, 0, 0, -c, s, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, -s, -c],
            [0, 0, 0, 0, -1, 0, 0, 0, 0, 0, c, -s],
        ],
        dtype=float,
    )


def physical_end_map(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
    """Map physical book-clamp reactions to the physical state at the joint."""

    return physical_state_transfer_matrix(omega, point) @ cantilever_clamp_matrix(
        point, "book_slope_clamp", scaled=False
    )


def _block_end_map(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    matrix = np.zeros((12, 6), dtype=np.complex128)
    matrix[:6, :3] = physical_end_map(omega, point_1)
    matrix[6:, 3:] = physical_end_map(omega, point_2)
    return matrix


def coupled_boundary_matrix_raw(
    omega: complex, beta: float, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    """Return physical ``D_coupled = J_book @ block_diag(H_1,H_2)``."""

    return joint_matrix_book(beta) @ _block_end_map(omega, point_1, point_2)


def equilibrate_matrix(
    matrix: NDArray[np.complex128],
) -> tuple[NDArray[np.complex128], MatrixScaling]:
    """Equilibrate by positive row/column max norms without moving determinant zeros."""

    value = np.asarray(matrix, dtype=np.complex128)
    row_norms = np.max(np.abs(value), axis=1)
    row_factors = np.ones(value.shape[0], dtype=float)
    nonzero_rows = row_norms > np.finfo(float).tiny
    row_factors[nonzero_rows] = 1.0 / row_norms[nonzero_rows]
    row_scaled = row_factors[:, np.newaxis] * value
    column_norms = np.max(np.abs(row_scaled), axis=0)
    column_factors = np.ones(value.shape[1], dtype=float)
    nonzero_columns = column_norms > np.finfo(float).tiny
    column_factors[nonzero_columns] = 1.0 / column_norms[nonzero_columns]
    scaled = row_scaled * column_factors[np.newaxis, :]
    return scaled, MatrixScaling(row_factors, column_factors)


def coupled_boundary_matrix(
    omega: complex,
    beta: float,
    point_1: RodPoint,
    point_2: RodPoint,
    *,
    scaled: bool = True,
) -> NDArray[np.complex128]:
    """Return the raw or positively equilibrated coupled ``6 x 6`` matrix."""

    raw = coupled_boundary_matrix_raw(omega, beta, point_1, point_2)
    return equilibrate_matrix(raw)[0] if scaled else raw


def straight_right_clamp_matrix(point: RodPoint) -> NDArray[np.complex128]:
    """Select ``w``, ``w'=psi+Q/K_s``, and ``Phi`` at the right end."""

    shear_stiffness = (
        point.geometry.shear_factor
        * point.properties.Gxz
        * point.geometry.area
    )
    matrix = np.zeros((3, 6), dtype=np.complex128)
    matrix[0, 0] = 1.0
    matrix[1, 1] = 1.0
    matrix[1, 3] = 1.0 / shear_stiffness
    matrix[2, 2] = 1.0
    return matrix


def straight_boundary_matrix_raw(
    omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Return the independent both-end book-slope-clamped ``3 x 3`` matrix."""

    return (
        straight_right_clamp_matrix(point)
        @ physical_state_transfer_matrix(omega, point)
        @ cantilever_clamp_matrix(point, "book_slope_clamp", scaled=False)
    )


def straight_boundary_matrix(
    omega: complex, point: RodPoint, *, scaled: bool = True
) -> NDArray[np.complex128]:
    """Return the raw or positively equilibrated straight reference matrix."""

    raw = straight_boundary_matrix_raw(omega, point)
    return equilibrate_matrix(raw)[0] if scaled else raw


def kinematic_joint_residual(
    beta: float, joint_rotation: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Return the two rotation rows for one compatible global rotation."""

    basis = joint_basis(beta)
    Phi_1, psi_1 = rotation_components(joint_rotation, basis.t_1, basis.n_1)
    Phi_2, psi_2 = rotation_components(joint_rotation, basis.t_2, basis.n_2)
    state = np.zeros(12, dtype=float)
    state[[1, 2, 7, 8]] = [psi_1, Phi_1, psi_2, Phi_2]
    return joint_matrix_book(beta)[1:3] @ state


def moment_joint_residual(
    beta: float, moment_1: NDArray[np.float64]
) -> NDArray[np.float64]:
    """Return the two moment rows for ``m_2=-m_1``."""

    basis = joint_basis(beta)
    moment_2 = -np.asarray(moment_1, dtype=float)
    M_T_1, M_1 = moment_components(moment_1, basis.t_1, basis.n_1)
    M_T_2, M_2 = moment_components(moment_2, basis.t_2, basis.n_2)
    state = np.zeros(12, dtype=float)
    state[[4, 5, 10, 11]] = [M_1, M_T_1, M_2, M_T_2]
    return joint_matrix_book(beta)[4:6] @ state


def joint_virtual_work(
    Q_1: float,
    Q_2: float,
    delta_w_1: float,
    delta_w_2: float,
    moment_1: NDArray[np.float64],
    moment_2: NDArray[np.float64],
    delta_rotation_1: NDArray[np.float64],
    delta_rotation_2: NDArray[np.float64],
) -> float:
    """Evaluate the invariant joint virtual-work expression."""

    return float(
        Q_1 * delta_w_1
        + Q_2 * delta_w_2
        + np.asarray(moment_1) @ np.asarray(delta_rotation_1)
        + np.asarray(moment_2) @ np.asarray(delta_rotation_2)
    )


__all__ = [
    "BOOK_STATE_SIZE",
    "JointBasis",
    "MatrixScaling",
    "coupled_boundary_matrix",
    "coupled_boundary_matrix_raw",
    "equilibrate_matrix",
    "joint_basis",
    "joint_matrix_book",
    "joint_matrix_old",
    "joint_virtual_work",
    "kinematic_joint_residual",
    "moment_components",
    "moment_joint_residual",
    "notation_transform_matrix",
    "physical_end_map",
    "physical_moment_vector",
    "physical_rotation_vector",
    "rotation_components",
    "straight_boundary_matrix",
    "straight_boundary_matrix_raw",
    "straight_right_clamp_matrix",
    "two_arm_notation_transform_matrix",
]
