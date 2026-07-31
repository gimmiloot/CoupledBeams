"""Rectangular orthotropic Euler--Bernoulli comparator for the Chapter-2 gate.

This script-owned helper deliberately stays narrower than a production API.  It
uses the accepted book state ``[w, psi, Phi, Q, M, M_T]``, the existing rigid
joint matrix, and the existing rectangular Saint-Venant stiffness.  It also
contains an independent displacement-based 1D FEM used only as a validation
comparator.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigh, expm
from scipy.optimize import brentq

from scripts.lib.yartsev_ch2_coupled_rods import (
    equilibrate_matrix,
    joint_matrix_book,
)
from scripts.lib.yartsev_ch2_monoclinic_rod import RodPoint


BOOK_STATE_SIZE = 6
FEM_NODE_DOF_COUNT = 3


@dataclass(frozen=True)
class EBFEMAssembly:
    """Full and constrained matrices for the independent two-arm FEM."""

    stiffness_full: NDArray[np.float64]
    mass_full: NDArray[np.float64]
    reduction: NDArray[np.float64]
    stiffness: NDArray[np.float64]
    mass: NDArray[np.float64]
    joint_full_dofs: tuple[NDArray[np.int64], NDArray[np.int64]]
    joint_maps: tuple[NDArray[np.float64], NDArray[np.float64]]
    elements_per_arm: int | tuple[int, int]
    element_counts: tuple[int, int]


@dataclass(frozen=True)
class EBFEMResult:
    """Lowest FEM eigenpairs and numerical diagnostics."""

    eigenvalues: NDArray[np.float64]
    frequencies_hz: NDArray[np.float64]
    eigenvectors_reduced: NDArray[np.float64]
    eigenvectors_full: NDArray[np.float64]
    stiffness_symmetry_residual: float
    mass_symmetry_residual: float
    minimum_mass_eigenvalue: float
    zero_mode_count: int
    joint_equilibrium_residuals: NDArray[np.float64]


def require_rectangular_orthotropic_eb(point: RodPoint) -> None:
    """Reject inputs outside the deliberately narrow elastic ``theta=0`` gate."""

    if point.material_mode != "elastic":
        raise ValueError("the rectangular EB comparator is elastic only")
    if abs(point.properties.theta_deg) > 1.0e-14:
        raise ValueError("the rectangular EB comparator is restricted to theta=0")
    compliance_scale = max(abs(point.properties.Sbar11), np.finfo(float).tiny)
    if abs(point.properties.Sbar16) > 1.0e-12 * compliance_scale:
        raise ValueError("theta=0 orthotropy requires Sbar16=0")
    stiffness_scale = max(abs(point.torsion.Cbar), np.finfo(float).tiny)
    if abs(point.torsion.C_T - point.torsion.Cbar) > 1.0e-12 * stiffness_scale:
        raise ValueError("theta=0 requires C_T=Cbar=C_SV")


def saint_venant_stiffness(point: RodPoint) -> float:
    """Return the real generalized Saint-Venant stiffness ``C_SV=Cbar``."""

    require_rectangular_orthotropic_eb(point)
    value = complex(point.torsion.Cbar)
    if abs(value.imag) > 1.0e-12 * max(abs(value.real), 1.0):
        raise ValueError("elastic Cbar must be real")
    return float(value.real)


def eb_state_matrix(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
    """Return the physical EB state matrix in accepted book ordering."""

    require_rectangular_orthotropic_eb(point)
    geometry = point.geometry
    rho = point.material.rho
    bending_stiffness = complex(point.properties.Ex) * geometry.I_y
    torsional_stiffness = saint_venant_stiffness(point)
    omega2 = complex(omega) ** 2

    matrix = np.zeros((BOOK_STATE_SIZE, BOOK_STATE_SIZE), dtype=np.complex128)
    matrix[0, 1] = 1.0
    matrix[1, 4] = 1.0 / bending_stiffness
    matrix[2, 5] = 1.0 / torsional_stiffness
    matrix[3, 0] = -rho * geometry.area * omega2
    matrix[4, 3] = -1.0
    matrix[5, 2] = -rho * geometry.I_p * omega2
    return matrix


def eb_state_transfer_matrix(
    omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Return the full physical-state transfer over one arm."""

    return expm(eb_state_matrix(omega, point) * point.geometry.length)


def eb_clamp_matrix(point: RodPoint) -> NDArray[np.complex128]:
    """Map clamp reactions ``[Q,M,M_T]`` to the EB book state."""

    require_rectangular_orthotropic_eb(point)
    matrix = np.zeros((BOOK_STATE_SIZE, 3), dtype=np.complex128)
    matrix[3:, :] = np.eye(3, dtype=np.complex128)
    return matrix


def eb_physical_end_map(
    omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Map physical clamp reactions to the physical joint-end EB state."""

    return eb_state_transfer_matrix(omega, point) @ eb_clamp_matrix(point)


def _block_end_map(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    matrix = np.zeros((12, 6), dtype=np.complex128)
    matrix[:6, :3] = eb_physical_end_map(omega, point_1)
    matrix[6:, 3:] = eb_physical_end_map(omega, point_2)
    return matrix


def eb_coupled_boundary_matrix_raw(
    omega: complex, beta: float, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    """Return physical ``D_EB=J_book @ blockdiag(H_1,H_2)``."""

    return joint_matrix_book(beta) @ _block_end_map(omega, point_1, point_2)


def eb_coupled_boundary_matrix(
    omega: complex,
    beta: float,
    point_1: RodPoint,
    point_2: RodPoint,
    *,
    scaled: bool = True,
) -> NDArray[np.complex128]:
    """Return the raw or positively equilibrated coupled EB matrix."""

    raw = eb_coupled_boundary_matrix_raw(omega, beta, point_1, point_2)
    return equilibrate_matrix(raw)[0] if scaled else raw


def eb_straight_right_clamp_matrix(point: RodPoint) -> NDArray[np.complex128]:
    """Select fixed right-end ``w``, ``psi``, and ``Phi``."""

    require_rectangular_orthotropic_eb(point)
    selector = np.zeros((3, BOOK_STATE_SIZE), dtype=np.complex128)
    selector[:, :3] = np.eye(3, dtype=np.complex128)
    return selector


def eb_straight_boundary_matrix_raw(
    omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Return the independent fixed--fixed straight EB boundary matrix."""

    return (
        eb_straight_right_clamp_matrix(point)
        @ eb_state_transfer_matrix(omega, point)
        @ eb_clamp_matrix(point)
    )


def eb_straight_boundary_matrix(
    omega: complex, point: RodPoint, *, scaled: bool = True
) -> NDArray[np.complex128]:
    """Return raw or positively equilibrated fixed--fixed straight matrix."""

    raw = eb_straight_boundary_matrix_raw(omega, point)
    return equilibrate_matrix(raw)[0] if scaled else raw


def eb_bending_coupled_boundary_matrix_raw(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    """Return the ``beta=0`` bending block for independent family checks."""

    full = eb_coupled_boundary_matrix_raw(omega, 0.0, point_1, point_2)
    return full[np.ix_([0, 2, 3, 5], [0, 1, 3, 4])]


def eb_bending_coupled_boundary_matrix(
    omega: complex, point_1: RodPoint, point_2: RodPoint
) -> NDArray[np.complex128]:
    """Return the positively equilibrated ``beta=0`` bending block."""

    return equilibrate_matrix(
        eb_bending_coupled_boundary_matrix_raw(omega, point_1, point_2)
    )[0]


def fixed_fixed_bending_dimensionless_roots(
    num_roots: int,
) -> NDArray[np.float64]:
    """Solve ``cosh(x) cos(x)-1=0`` for its first positive elastic roots."""

    if num_roots < 1:
        raise ValueError("num_roots must be positive")
    roots: list[float] = []
    for mode in range(1, num_roots + 1):
        center = (mode + 0.5) * np.pi
        root = brentq(
            lambda value: np.cosh(value) * np.cos(value) - 1.0,
            center - 0.24 * np.pi,
            center + 0.24 * np.pi,
            xtol=1.0e-14,
            rtol=4.0 * np.finfo(float).eps,
        )
        roots.append(float(root))
    return np.asarray(roots)


def fixed_fixed_bending_frequencies_hz(
    point: RodPoint, num_roots: int
) -> NDArray[np.float64]:
    """Return exact fixed--fixed rectangular EB bending frequencies."""

    require_rectangular_orthotropic_eb(point)
    geometry = point.geometry
    bending_stiffness = float(np.real(point.properties.Ex)) * geometry.I_y
    factor = np.sqrt(bending_stiffness / (point.material.rho * geometry.area))
    roots = fixed_fixed_bending_dimensionless_roots(num_roots)
    omega = roots**2 * factor / geometry.length**2
    return omega / (2.0 * np.pi)


def fixed_fixed_torsion_frequencies_hz(
    point: RodPoint, num_roots: int
) -> NDArray[np.float64]:
    """Return exact fixed--fixed Saint-Venant torsion frequencies."""

    if num_roots < 1:
        raise ValueError("num_roots must be positive")
    geometry = point.geometry
    wave_speed = np.sqrt(
        saint_venant_stiffness(point) / (point.material.rho * geometry.I_p)
    )
    modes = np.arange(1, num_roots + 1, dtype=float)
    omega = modes * np.pi * wave_speed / geometry.length
    return omega / (2.0 * np.pi)


def eb_element_matrices(
    point: RodPoint, element_length: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return one ``[w,psi,Phi]`` EB plus torsion element in SI units."""

    require_rectangular_orthotropic_eb(point)
    if element_length <= 0.0:
        raise ValueError("element_length must be positive")
    length = float(element_length)
    geometry = point.geometry
    rho = point.material.rho
    bending_stiffness = float(np.real(point.properties.Ex)) * geometry.I_y

    bending_k = bending_stiffness / length**3 * np.array(
        [
            [12.0, 6.0 * length, -12.0, 6.0 * length],
            [6.0 * length, 4.0 * length**2, -6.0 * length, 2.0 * length**2],
            [-12.0, -6.0 * length, 12.0, -6.0 * length],
            [6.0 * length, 2.0 * length**2, -6.0 * length, 4.0 * length**2],
        ]
    )
    bending_m = rho * geometry.area * length / 420.0 * np.array(
        [
            [156.0, 22.0 * length, 54.0, -13.0 * length],
            [22.0 * length, 4.0 * length**2, 13.0 * length, -3.0 * length**2],
            [54.0, 13.0 * length, 156.0, -22.0 * length],
            [-13.0 * length, -3.0 * length**2, -22.0 * length, 4.0 * length**2],
        ]
    )
    torsion_k = saint_venant_stiffness(point) / length * np.array(
        [[1.0, -1.0], [-1.0, 1.0]]
    )
    torsion_m = rho * geometry.I_p * length / 6.0 * np.array(
        [[2.0, 1.0], [1.0, 2.0]]
    )

    stiffness = np.zeros((6, 6), dtype=float)
    mass = np.zeros((6, 6), dtype=float)
    bending_dofs = [0, 1, 3, 4]
    torsion_dofs = [2, 5]
    stiffness[np.ix_(bending_dofs, bending_dofs)] = bending_k
    mass[np.ix_(bending_dofs, bending_dofs)] = bending_m
    stiffness[np.ix_(torsion_dofs, torsion_dofs)] = torsion_k
    mass[np.ix_(torsion_dofs, torsion_dofs)] = torsion_m
    return stiffness, mass


def eb_joint_dof_maps(beta: float) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Map common ``[w_J,theta_t,theta_n]`` to each arm's endpoint DOFs."""

    c = float(np.cos(beta))
    s = float(np.sin(beta))
    arm_1 = np.array(
        [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]]
    )
    arm_2 = np.array(
        [[1.0, 0.0, 0.0], [0.0, s, c], [0.0, -c, s]]
    )
    return arm_1, arm_2


def eb_joint_mapping_residual(beta: float) -> float:
    """Return the maximum kinematic ``J_book`` residual of the FEM maps."""

    maps = eb_joint_dof_maps(beta)
    transform = np.zeros((12, 3), dtype=float)
    for arm, mapping in enumerate(maps):
        offset = 6 * arm
        transform[offset + 0] = mapping[0]
        transform[offset + 1] = mapping[1]
        transform[offset + 2] = mapping[2]
    residual = joint_matrix_book(beta)[:3] @ transform
    return float(np.max(np.abs(residual)))


def assemble_two_arm_eb_fem(
    point_1: RodPoint,
    point_2: RodPoint,
    beta: float,
    elements_per_arm: int | tuple[int, int],
) -> EBFEMAssembly:
    """Assemble an independent two-arm EB plus torsion FEM by congruence.

    Passing one integer retains the original equal-count mesh.  A two-integer
    tuple gives the arm-specific counts without changing element matrices,
    joint maps, clamps, or the eigensolver contract.
    """

    require_rectangular_orthotropic_eb(point_1)
    require_rectangular_orthotropic_eb(point_2)
    if isinstance(elements_per_arm, (int, np.integer)) and not isinstance(
        elements_per_arm, (bool, np.bool_)
    ):
        element_counts = (int(elements_per_arm), int(elements_per_arm))
    elif isinstance(elements_per_arm, tuple) and len(elements_per_arm) == 2:
        if any(
            not isinstance(value, (int, np.integer))
            or isinstance(value, (bool, np.bool_))
            for value in elements_per_arm
        ):
            raise TypeError("arm-specific element counts must be integers")
        element_counts = (int(elements_per_arm[0]), int(elements_per_arm[1]))
    else:
        raise TypeError("elements_per_arm must be an integer or a two-integer tuple")
    if any(count < 1 for count in element_counts):
        raise ValueError("each arm must have at least one element")

    nodes_per_arm = tuple(count + 1 for count in element_counts)
    full_dofs_per_arm = tuple(
        FEM_NODE_DOF_COUNT * node_count for node_count in nodes_per_arm
    )
    arm_offsets = (0, full_dofs_per_arm[0])
    total_full_dofs = sum(full_dofs_per_arm)
    stiffness_full = np.zeros((total_full_dofs, total_full_dofs), dtype=float)
    mass_full = np.zeros_like(stiffness_full)

    for point, element_count, arm_offset in zip(
        (point_1, point_2), element_counts, arm_offsets
    ):
        element_k, element_m = eb_element_matrices(
            point, point.geometry.length / element_count
        )
        for element in range(element_count):
            start = arm_offset + FEM_NODE_DOF_COUNT * element
            dofs = np.arange(start, start + 2 * FEM_NODE_DOF_COUNT)
            stiffness_full[np.ix_(dofs, dofs)] += element_k
            mass_full[np.ix_(dofs, dofs)] += element_m

    reduced_size = (
        FEM_NODE_DOF_COUNT * sum(count - 1 for count in element_counts) + 3
    )
    reduction = np.zeros((total_full_dofs, reduced_size), dtype=float)
    cursor = 0
    for element_count, arm_offset in zip(element_counts, arm_offsets):
        for node in range(1, element_count):
            rows = arm_offset + FEM_NODE_DOF_COUNT * node + np.arange(3)
            columns = cursor + np.arange(3)
            reduction[np.ix_(rows, columns)] = np.eye(3)
            cursor += 3

    joint_maps = eb_joint_dof_maps(beta)
    joint_full_dofs: list[NDArray[np.int64]] = []
    joint_columns = np.arange(cursor, cursor + 3)
    for mapping, element_count, arm_offset in zip(
        joint_maps, element_counts, arm_offsets
    ):
        rows = arm_offset + FEM_NODE_DOF_COUNT * element_count + np.arange(3)
        reduction[np.ix_(rows, joint_columns)] = mapping
        joint_full_dofs.append(rows.astype(np.int64))
    if cursor + 3 != reduced_size:
        raise RuntimeError("internal reduced-DOF count mismatch")

    stiffness = reduction.T @ stiffness_full @ reduction
    mass = reduction.T @ mass_full @ reduction
    return EBFEMAssembly(
        stiffness_full=stiffness_full,
        mass_full=mass_full,
        reduction=reduction,
        stiffness=stiffness,
        mass=mass,
        joint_full_dofs=(joint_full_dofs[0], joint_full_dofs[1]),
        joint_maps=joint_maps,
        elements_per_arm=elements_per_arm,
        element_counts=element_counts,
    )


def _relative_symmetry_residual(matrix: NDArray[np.float64]) -> float:
    scale = max(float(np.linalg.norm(matrix)), np.finfo(float).tiny)
    return float(np.linalg.norm(matrix - matrix.T) / scale)


def solve_two_arm_eb_fem(
    assembly: EBFEMAssembly, *, num_roots: int = 7
) -> EBFEMResult:
    """Solve the lowest constrained FEM modes and audit joint equilibrium."""

    if num_roots < 1 or num_roots >= assembly.stiffness.shape[0]:
        raise ValueError("num_roots must be positive and below reduced size")
    minimum_mass_eigenvalue = float(np.linalg.eigvalsh(assembly.mass)[0])
    eigenvalues, eigenvectors = eigh(
        assembly.stiffness,
        assembly.mass,
        subset_by_index=[0, num_roots - 1],
        check_finite=True,
    )
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    zero_mode_count = int(np.count_nonzero(eigenvalues <= 1.0e-12 * scale))
    clipped = np.maximum(eigenvalues, 0.0)
    frequencies = np.sqrt(clipped) / (2.0 * np.pi)
    full_vectors = assembly.reduction @ eigenvectors

    equilibrium: list[float] = []
    for mode, eigenvalue in enumerate(eigenvalues):
        vector = full_vectors[:, mode]
        elastic_force = assembly.stiffness_full @ vector
        inertia_force = float(eigenvalue) * assembly.mass_full @ vector
        residual = (
            elastic_force
            - inertia_force
        )
        endpoint_forces = [
            residual[dofs] for dofs in assembly.joint_full_dofs
        ]
        joint_residual = sum(
            mapping.T @ force
            for mapping, force in zip(assembly.joint_maps, endpoint_forces)
        )
        reaction_scale = max(
            sum(float(np.linalg.norm(force)) for force in endpoint_forces),
            float(np.linalg.norm(elastic_force)) + float(np.linalg.norm(inertia_force)),
            np.finfo(float).tiny,
        )
        equilibrium.append(float(np.linalg.norm(joint_residual) / reaction_scale))

    return EBFEMResult(
        eigenvalues=np.asarray(eigenvalues, dtype=float),
        frequencies_hz=np.asarray(frequencies, dtype=float),
        eigenvectors_reduced=np.asarray(eigenvectors, dtype=float),
        eigenvectors_full=np.asarray(full_vectors, dtype=float),
        stiffness_symmetry_residual=_relative_symmetry_residual(assembly.stiffness),
        mass_symmetry_residual=_relative_symmetry_residual(assembly.mass),
        minimum_mass_eigenvalue=minimum_mass_eigenvalue,
        zero_mode_count=zero_mode_count,
        joint_equilibrium_residuals=np.asarray(equilibrium, dtype=float),
    )


__all__ = [
    "BOOK_STATE_SIZE",
    "EBFEMAssembly",
    "EBFEMResult",
    "FEM_NODE_DOF_COUNT",
    "assemble_two_arm_eb_fem",
    "eb_bending_coupled_boundary_matrix",
    "eb_bending_coupled_boundary_matrix_raw",
    "eb_clamp_matrix",
    "eb_coupled_boundary_matrix",
    "eb_coupled_boundary_matrix_raw",
    "eb_element_matrices",
    "eb_joint_dof_maps",
    "eb_joint_mapping_residual",
    "eb_physical_end_map",
    "eb_state_matrix",
    "eb_state_transfer_matrix",
    "eb_straight_boundary_matrix",
    "eb_straight_boundary_matrix_raw",
    "eb_straight_right_clamp_matrix",
    "fixed_fixed_bending_dimensionless_roots",
    "fixed_fixed_bending_frequencies_hz",
    "fixed_fixed_torsion_frequencies_hz",
    "require_rectangular_orthotropic_eb",
    "saint_venant_stiffness",
    "solve_two_arm_eb_fem",
]
