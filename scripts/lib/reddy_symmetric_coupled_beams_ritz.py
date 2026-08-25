"""Independent two-arm Rayleigh--Ritz model for symmetric Reddy beams.

This narrowly scoped module discretizes the physical energies of two arms and
imposes only the three rigid-joint *kinematic* constraints.  It deliberately
contains no transfer matrix, matrix exponential, determinant, root finder, or
joint force-equilibrium row.  Reduced single-arm properties and the frozen
physical coordinate maps are the only project-level scientific inputs.

Coefficient blocks are ordered as
``[u1, w1, psi1, u2, w2, psi2]`` with ``order`` coefficients per block.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

import numpy as np
from numpy.polynomial import Legendre
from numpy.typing import ArrayLike, NDArray
from scipy.linalg import eigh

from scripts.lib import reddy_inplane_geometry as inplane_geometry
from scripts.lib.reddy_symmetric_laminated_beam import BeamProperties


FloatArray = NDArray[np.float64]

GAUSS_LEGENDRE_ORDER = 64
RITZ_ORDERS = (4, 6, 8, 10, 12)
RITZ_GUARD_ORDER = 16
COEFFICIENT_BLOCK_ORDER = ("u1", "w1", "psi1", "u2", "w2", "psi2")
CONSTRAINT_RANK_RELATIVE_THRESHOLD = 1.0e-12


def _readonly(values: ArrayLike) -> FloatArray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _positive(value: float, *, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def _finite_angle(beta_rad: float) -> float:
    result = float(beta_rad)
    if not math.isfinite(result):
        raise ValueError("beta_rad must be finite")
    return result


def _order(value: int) -> int:
    if not isinstance(value, (int, np.integer)) or int(value) < 1:
        raise ValueError("order must be a positive integer")
    return int(value)


def _symmetry_residual(matrix: FloatArray) -> float:
    denominator = max(float(np.linalg.norm(matrix)), np.finfo(float).tiny)
    return float(np.linalg.norm(matrix - matrix.T) / denominator)


def _relative_residual(numerator: float, denominator: float) -> float:
    return float(numerator / max(denominator, np.finfo(float).tiny))


def field_slice(order: int, arm: int, field: str) -> slice:
    """Return one coefficient-block slice in the frozen six-block ordering."""

    size = _order(order)
    if arm not in (1, 2):
        raise ValueError("arm must be 1 or 2")
    try:
        field_offset = {"u": 0, "w": 1, "psi": 2}[str(field)]
    except KeyError as error:
        raise ValueError("field must be 'u', 'w', or 'psi'") from error
    block = (arm - 1) * 3 + field_offset
    return slice(block * size, (block + 1) * size)


@dataclass(frozen=True, eq=False)
class BasisEvaluation:
    """Values and derivatives with respect to ``xi``; shape ``(N,n)``."""

    values: FloatArray
    derivatives: FloatArray


@dataclass(frozen=True, eq=False)
class CoupledRitzMatrices:
    """Unconstrained energy Hessians and their assembly diagnostics."""

    order: int
    gauss_order: int
    length_1: float
    length_2: float
    properties_1: BeamProperties
    properties_2: BeamProperties
    stiffness: FloatArray
    mass: FloatArray
    stiffness_symmetry_residual: float
    mass_symmetry_residual: float
    full_stiffness_minimum_eigenvalue: float
    full_mass_minimum_eigenvalue: float


@dataclass(frozen=True, eq=False)
class ConstraintReduction:
    """Three-row joint constraint and an orthonormal basis of its kernel."""

    beta_rad: float
    constraint: FloatArray
    singular_values: FloatArray
    rank: int
    threshold: float
    nullspace: FloatArray
    orthonormality_residual: float
    constraint_nullspace_residual: float


@dataclass(frozen=True, eq=False)
class RitzSpectrum:
    """Mass-orthonormal constrained Ritz eigenspectrum."""

    matrices: CoupledRitzMatrices
    reduction: ConstraintReduction
    reduced_stiffness: FloatArray
    reduced_mass: FloatArray
    eigenvalues: FloatArray
    omegas: FloatArray
    reduced_eigenvectors: FloatArray
    coefficients: FloatArray
    reduced_stiffness_minimum_eigenvalue: float
    reduced_mass_minimum_eigenvalue: float
    reduced_stiffness_spd: bool
    reduced_mass_spd: bool
    zero_or_negative_mode_count: int
    maximum_eigenpair_backward_residual: float
    maximum_rayleigh_relative_residual: float
    mass_orthonormality_residual: float
    maximum_constraint_residual: float


@dataclass(frozen=True, eq=False)
class ArmModeFields:
    """One arm's fields, physical derivatives, and constitutive resultants."""

    xi: FloatArray
    x: FloatArray
    u: FloatArray
    w: FloatArray
    psi: FloatArray
    u_prime: FloatArray
    w_prime: FloatArray
    psi_prime: FloatArray
    N: FloatArray
    Q: FloatArray
    M: FloatArray


@dataclass(frozen=True, eq=False)
class JointResiduals:
    """Kinematic and naturally recovered physical joint residuals."""

    displacement_components: FloatArray
    rotation_component: float
    force_components: FloatArray
    moment_component: float
    displacement_absolute: float
    rotation_absolute: float
    force_absolute: float
    moment_absolute: float
    kinematic_scale: float
    force_scale: float
    moment_scale: float
    displacement_normalized: float
    rotation_normalized: float
    force_normalized: float
    moment_normalized: float


@dataclass(frozen=True)
class EnergyDiagnostics:
    """Six additive energy diagnostics; no modal-family classification."""

    T_axial: float
    T_transverse: float
    T_rotation: float
    U_axial: float
    U_bending: float
    U_shear: float
    total_mass_norm: float
    total_stiffness_norm: float
    inertia_norm: float
    energy_identity_relative: float
    T_axial_share: float
    T_transverse_share: float
    T_rotation_share: float
    U_axial_share: float
    U_bending_share: float
    U_shear_share: float


@dataclass(frozen=True, eq=False)
class SubspaceComparison:
    """Mass-inner-product overlap of two modal subspaces."""

    singular_values: FloatArray
    principal_angles_rad: FloatArray
    minimum_singular_value: float
    maximum_principal_angle_rad: float


def shifted_legendre_clamped_basis(
    order: int, xi: Sequence[float] | FloatArray
) -> BasisEvaluation:
    """Evaluate ``phi_n=xi*P_n(2*xi-1)`` and ``d phi_n/d xi``.

    Only the outer-clamp value at ``xi=0`` is imposed.  Endpoint values and
    slopes at ``xi=1`` remain unrestricted.
    """

    size = _order(order)
    points = np.asarray(xi, dtype=float)
    if points.ndim != 1 or not np.all(np.isfinite(points)):
        raise ValueError("xi must be a finite one-dimensional array")
    if np.any(points < 0.0) or np.any(points > 1.0):
        raise ValueError("xi points must lie in [0,1]")
    coordinate = 2.0 * points - 1.0
    values = np.empty((size, points.size), dtype=float)
    derivatives = np.empty_like(values)
    for degree in range(size):
        polynomial = Legendre.basis(degree)
        raw = polynomial(coordinate)
        raw_derivative = 2.0 * polynomial.deriv()(coordinate)
        values[degree] = points * raw
        derivatives[degree] = raw + points * raw_derivative
    return BasisEvaluation(_readonly(values), _readonly(derivatives))


def gauss_legendre_rule(order: int = GAUSS_LEGENDRE_ORDER) -> tuple[FloatArray, FloatArray]:
    """Return the fixed rule on ``[0,1]`` used by the Ritz implementation."""

    if int(order) != GAUSS_LEGENDRE_ORDER:
        raise ValueError(
            f"RLB-1C fixes Gauss--Legendre order at {GAUSS_LEGENDRE_ORDER}"
        )
    raw_nodes, raw_weights = np.polynomial.legendre.leggauss(GAUSS_LEGENDRE_ORDER)
    return _readonly(0.5 * (raw_nodes + 1.0)), _readonly(0.5 * raw_weights)


def global_kinematic_map(beta_rad: float, arm: int) -> FloatArray:
    """Return ``[d_X,d_Y,psi]^T = G_i [u,w,psi]^T`` from geometry maps."""

    angle = _finite_angle(beta_rad)
    geometry = inplane_geometry.reddy_inplane_geometry(math.degrees(angle))
    if arm == 1:
        basis = geometry.arm1
    elif arm == 2:
        basis = geometry.arm2
    else:
        raise ValueError("arm must be 1 or 2")
    result = np.zeros((3, 3), dtype=float)
    result[:2, :2] = np.asarray(basis.inplane_matrix[:2, :], dtype=float)
    # Both physical rotations are components along the common axis -k=E_Z.
    result[2, 2] = 1.0
    return _readonly(result)


def kinematic_constraint_matrix(beta_rad: float, order: int) -> FloatArray:
    """Build the three-row endpoint constraint ``G1*q1J-G2*q2J=0``."""

    size = _order(order)
    endpoint = shifted_legendre_clamped_basis(size, np.array([1.0])).values[:, 0]
    selector_1 = np.zeros((3, 6 * size), dtype=float)
    selector_2 = np.zeros_like(selector_1)
    for index, field in enumerate(("u", "w", "psi")):
        selector_1[index, field_slice(size, 1, field)] = endpoint
        selector_2[index, field_slice(size, 2, field)] = endpoint
    constraint = (
        np.asarray(global_kinematic_map(beta_rad, 1)) @ selector_1
        - np.asarray(global_kinematic_map(beta_rad, 2)) @ selector_2
    )
    return _readonly(constraint)


def constraint_nullspace(
    beta_rad: float,
    order: int,
    *,
    relative_threshold: float = CONSTRAINT_RANK_RELATIVE_THRESHOLD,
) -> ConstraintReduction:
    """Compute the canonical SVD nullspace of the three kinematic rows."""

    threshold_relative = _positive(relative_threshold, name="relative_threshold")
    constraint = np.asarray(kinematic_constraint_matrix(beta_rad, order), dtype=float)
    _left, singular_values, right = np.linalg.svd(constraint, full_matrices=True)
    reference = float(singular_values[0]) if singular_values.size else 0.0
    threshold = threshold_relative * reference
    rank = int(np.count_nonzero(singular_values > threshold))
    nullspace = np.asarray(right[rank:, :].T, dtype=float)
    identity = np.eye(nullspace.shape[1])
    orthonormality = _relative_residual(
        float(np.linalg.norm(nullspace.T @ nullspace - identity)),
        max(float(np.linalg.norm(identity)), 1.0),
    )
    cz = _relative_residual(
        float(np.linalg.norm(constraint @ nullspace)),
        float(np.linalg.norm(constraint) * np.linalg.norm(nullspace)),
    )
    return ConstraintReduction(
        beta_rad=float(beta_rad),
        constraint=_readonly(constraint),
        singular_values=_readonly(singular_values),
        rank=rank,
        threshold=threshold,
        nullspace=_readonly(nullspace),
        orthonormality_residual=orthonormality,
        constraint_nullspace_residual=cz,
    )


def _arm_energy_matrices(
    properties: BeamProperties,
    length: float,
    values: FloatArray,
    derivatives: FloatArray,
    weights: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    L = _positive(length, name="length")

    def inner(left: FloatArray, right: FloatArray) -> FloatArray:
        return (left * weights[None, :]) @ right.T

    gram_0 = inner(values, values)
    gram_1 = inner(derivatives, derivatives)
    derivative_value = inner(derivatives, values)
    size = values.shape[0]
    stiffness = np.zeros((3 * size, 3 * size), dtype=float)
    mass = np.zeros_like(stiffness)
    u = slice(0, size)
    w = slice(size, 2 * size)
    psi = slice(2 * size, 3 * size)
    stiffness[u, u] = properties.A / L * gram_1
    stiffness[w, w] = properties.S / L * gram_1
    stiffness[w, psi] = properties.S * derivative_value
    stiffness[psi, w] = properties.S * derivative_value.T
    stiffness[psi, psi] = (
        properties.D / L * gram_1 + properties.S * L * gram_0
    )
    mass[u, u] = properties.m * L * gram_0
    mass[w, w] = properties.m * L * gram_0
    mass[psi, psi] = properties.J * L * gram_0
    return stiffness, mass


def assemble_coupled_ritz_matrices(
    properties_1: BeamProperties,
    length_1: float,
    properties_2: BeamProperties,
    length_2: float,
    order: int,
    *,
    gauss_order: int = GAUSS_LEGENDRE_ORDER,
) -> CoupledRitzMatrices:
    """Assemble physical two-arm energy matrices before joint constraints."""

    size = _order(order)
    L1 = _positive(length_1, name="length_1")
    L2 = _positive(length_2, name="length_2")
    nodes, weights = gauss_legendre_rule(gauss_order)
    basis = shifted_legendre_clamped_basis(size, nodes)
    K1, M1 = _arm_energy_matrices(
        properties_1, L1, basis.values, basis.derivatives, weights
    )
    K2, M2 = _arm_energy_matrices(
        properties_2, L2, basis.values, basis.derivatives, weights
    )
    stiffness = np.zeros((6 * size, 6 * size), dtype=float)
    mass = np.zeros_like(stiffness)
    stiffness[: 3 * size, : 3 * size] = K1
    stiffness[3 * size :, 3 * size :] = K2
    mass[: 3 * size, : 3 * size] = M1
    mass[3 * size :, 3 * size :] = M2
    stiffness_symmetry = _symmetry_residual(stiffness)
    mass_symmetry = _symmetry_residual(mass)
    # Retain the exact energy-Hessian representation while removing only
    # floating-point antisymmetry from independently evaluated products.
    stiffness = 0.5 * (stiffness + stiffness.T)
    mass = 0.5 * (mass + mass.T)
    return CoupledRitzMatrices(
        order=size,
        gauss_order=gauss_order,
        length_1=L1,
        length_2=L2,
        properties_1=properties_1,
        properties_2=properties_2,
        stiffness=_readonly(stiffness),
        mass=_readonly(mass),
        stiffness_symmetry_residual=stiffness_symmetry,
        mass_symmetry_residual=mass_symmetry,
        full_stiffness_minimum_eigenvalue=float(np.linalg.eigvalsh(stiffness)[0]),
        full_mass_minimum_eigenvalue=float(np.linalg.eigvalsh(mass)[0]),
    )


def solve_coupled_ritz_spectrum(
    properties_1: BeamProperties,
    length_1: float,
    properties_2: BeamProperties,
    length_2: float,
    beta_rad: float,
    order: int,
    *,
    gauss_order: int = GAUSS_LEGENDRE_ORDER,
    rank_relative_threshold: float = CONSTRAINT_RANK_RELATIVE_THRESHOLD,
) -> RitzSpectrum:
    """Solve ``Z.T K Z v = omega**2 Z.T M Z v`` without transfer data."""

    matrices = assemble_coupled_ritz_matrices(
        properties_1,
        length_1,
        properties_2,
        length_2,
        order,
        gauss_order=gauss_order,
    )
    reduction = constraint_nullspace(
        beta_rad, order, relative_threshold=rank_relative_threshold
    )
    if reduction.rank != 3:
        raise ArithmeticError(f"kinematic constraint rank is {reduction.rank}, expected 3")
    Z = np.asarray(reduction.nullspace)
    stiffness = np.asarray(matrices.stiffness)
    mass = np.asarray(matrices.mass)
    reduced_stiffness = 0.5 * (Z.T @ stiffness @ Z + Z.T @ stiffness.T @ Z)
    reduced_mass = 0.5 * (Z.T @ mass @ Z + Z.T @ mass.T @ Z)
    stiffness_eigenvalues = np.linalg.eigvalsh(reduced_stiffness)
    mass_eigenvalues = np.linalg.eigvalsh(reduced_mass)
    stiffness_spd = bool(stiffness_eigenvalues[0] > 0.0)
    mass_spd = bool(mass_eigenvalues[0] > 0.0)
    if not stiffness_spd or not mass_spd:
        raise ArithmeticError("constrained Ritz stiffness and mass must be positive definite")
    eigenvalues, reduced_vectors = eigh(
        reduced_stiffness, reduced_mass, check_finite=True
    )
    tolerance = 100.0 * np.finfo(float).eps * max(
        float(np.max(np.abs(eigenvalues))), 1.0
    )
    nonpositive = int(np.count_nonzero(eigenvalues <= tolerance))
    if nonpositive:
        raise ArithmeticError("constrained Ritz spectrum contains zero/negative modes")
    coefficients = Z @ reduced_vectors
    # Fix each isolated coordinate vector's sign deterministically.  This does
    # not assign a physical meaning to native vectors in a repeated eigenspace.
    for column in range(coefficients.shape[1]):
        pivot = int(np.argmax(np.abs(coefficients[:, column])))
        if coefficients[pivot, column] < 0.0:
            coefficients[:, column] *= -1.0
            reduced_vectors[:, column] *= -1.0

    backward: list[float] = []
    rayleigh: list[float] = []
    constraints: list[float] = []
    for index, eigenvalue in enumerate(eigenvalues):
        vector = reduced_vectors[:, index]
        residual = reduced_stiffness @ vector - eigenvalue * reduced_mass @ vector
        denominator = (
            np.linalg.norm(reduced_stiffness) + abs(float(eigenvalue)) * np.linalg.norm(reduced_mass)
        ) * np.linalg.norm(vector)
        backward.append(_relative_residual(float(np.linalg.norm(residual)), float(denominator)))
        coefficient = coefficients[:, index]
        numerator = float(coefficient @ stiffness @ coefficient)
        denominator_mass = float(coefficient @ mass @ coefficient)
        quotient = numerator / denominator_mass
        rayleigh.append(
            abs(quotient - float(eigenvalue))
            / max(abs(quotient), abs(float(eigenvalue)), np.finfo(float).tiny)
        )
        constraints.append(
            _relative_residual(
                float(np.linalg.norm(reduction.constraint @ coefficient)),
                float(np.linalg.norm(reduction.constraint) * np.linalg.norm(coefficient)),
            )
        )
    mass_orthogonality = _relative_residual(
        float(np.linalg.norm(coefficients.T @ mass @ coefficients - np.eye(coefficients.shape[1]))),
        math.sqrt(float(coefficients.shape[1])),
    )
    return RitzSpectrum(
        matrices=matrices,
        reduction=reduction,
        reduced_stiffness=_readonly(reduced_stiffness),
        reduced_mass=_readonly(reduced_mass),
        eigenvalues=_readonly(eigenvalues),
        omegas=_readonly(np.sqrt(eigenvalues)),
        reduced_eigenvectors=_readonly(reduced_vectors),
        coefficients=_readonly(coefficients),
        reduced_stiffness_minimum_eigenvalue=float(stiffness_eigenvalues[0]),
        reduced_mass_minimum_eigenvalue=float(mass_eigenvalues[0]),
        reduced_stiffness_spd=stiffness_spd,
        reduced_mass_spd=mass_spd,
        zero_or_negative_mode_count=nonpositive,
        maximum_eigenpair_backward_residual=max(backward, default=0.0),
        maximum_rayleigh_relative_residual=max(rayleigh, default=0.0),
        mass_orthonormality_residual=mass_orthogonality,
        maximum_constraint_residual=max(constraints, default=0.0),
    )


def evaluate_arm_mode(
    coefficients: ArrayLike,
    order: int,
    arm: int,
    xi: Sequence[float] | FloatArray,
    length: float,
    properties: BeamProperties,
) -> ArmModeFields:
    """Evaluate one Ritz arm and recover physical constitutive resultants."""

    size = _order(order)
    vector = np.asarray(coefficients, dtype=float)
    if vector.shape != (6 * size,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"coefficients must have shape ({6 * size},)")
    L = _positive(length, name="length")
    points = np.asarray(xi, dtype=float)
    basis = shifted_legendre_clamped_basis(size, points)

    def field(name: str) -> tuple[FloatArray, FloatArray]:
        block = vector[field_slice(size, arm, name)]
        values = block @ basis.values
        derivatives = block @ basis.derivatives / L
        return np.asarray(values, dtype=float), np.asarray(derivatives, dtype=float)

    u, u_prime = field("u")
    w, w_prime = field("w")
    psi, psi_prime = field("psi")
    axial = properties.A * u_prime
    shear = properties.S * (w_prime + psi)
    moment = properties.D * psi_prime
    return ArmModeFields(
        xi=_readonly(points),
        x=_readonly(L * points),
        u=_readonly(u),
        w=_readonly(w),
        psi=_readonly(psi),
        u_prime=_readonly(u_prime),
        w_prime=_readonly(w_prime),
        psi_prime=_readonly(psi_prime),
        N=_readonly(axial),
        Q=_readonly(shear),
        M=_readonly(moment),
    )


def joint_residuals_from_ritz_mode(
    coefficients: ArrayLike,
    order: int,
    beta_rad: float,
    length_1: float,
    properties_1: BeamProperties,
    length_2: float,
    properties_2: BeamProperties,
) -> JointResiduals:
    """Recover compatibility and natural force/moment equilibrium at the joint."""

    angle = _finite_angle(beta_rad)
    L1 = _positive(length_1, name="length_1")
    L2 = _positive(length_2, name="length_2")
    geometry = inplane_geometry.reddy_inplane_geometry(math.degrees(angle))
    endpoints = np.array([0.0, 1.0])
    arm1 = evaluate_arm_mode(coefficients, order, 1, endpoints, L1, properties_1)
    arm2 = evaluate_arm_mode(coefficients, order, 2, endpoints, L2, properties_2)
    q1 = np.array([arm1.u[-1], arm1.w[-1], arm1.psi[-1]])
    q2 = np.array([arm2.u[-1], arm2.w[-1], arm2.psi[-1]])
    global1 = np.asarray(global_kinematic_map(angle, 1)) @ q1
    global2 = np.asarray(global_kinematic_map(angle, 2)) @ q2
    displacement = global1[:2] - global2[:2]
    rotation = float(global1[2] - global2[2])
    force1 = inplane_geometry.force_vector(geometry.arm1, arm1.N[-1], arm1.Q[-1])
    force2 = inplane_geometry.force_vector(geometry.arm2, arm2.N[-1], arm2.Q[-1])
    force = np.asarray(force1 + force2, dtype=float)
    moment = float(arm1.M[-1] + arm2.M[-1])
    total_length = L1 + L2
    kinematic_scale = max(
        float(np.linalg.norm(global1[:2])),
        float(np.linalg.norm(global2[:2])),
        total_length * abs(float(global1[2])),
        total_length * abs(float(global2[2])),
        np.finfo(float).tiny,
    )
    rotation_scale = kinematic_scale / total_length

    # A joint can be a nodal point.  The scale therefore also contains the
    # outer-clamp reactions; it is frozen independently of the residual.
    force_candidates: list[float] = []
    for basis, fields in ((geometry.arm1, arm1), (geometry.arm2, arm2)):
        for endpoint_index in (0, 1):
            endpoint_force = inplane_geometry.force_vector(
                basis, fields.N[endpoint_index], fields.Q[endpoint_index]
            )
            force_candidates.extend(
                [
                    float(np.linalg.norm(endpoint_force)),
                    abs(float(fields.M[endpoint_index])) / total_length,
                ]
            )
    force_scale = max(*force_candidates, np.finfo(float).tiny)
    moment_scale = force_scale * total_length
    displacement_absolute = float(np.linalg.norm(displacement))
    rotation_absolute = abs(rotation)
    force_absolute = float(np.linalg.norm(force))
    moment_absolute = abs(moment)
    return JointResiduals(
        displacement_components=_readonly(displacement),
        rotation_component=rotation,
        force_components=_readonly(force),
        moment_component=moment,
        displacement_absolute=displacement_absolute,
        rotation_absolute=rotation_absolute,
        force_absolute=force_absolute,
        moment_absolute=moment_absolute,
        kinematic_scale=kinematic_scale,
        force_scale=force_scale,
        moment_scale=moment_scale,
        displacement_normalized=displacement_absolute / kinematic_scale,
        rotation_normalized=rotation_absolute / rotation_scale,
        force_normalized=force_absolute / force_scale,
        moment_normalized=moment_absolute / moment_scale,
    )


def modal_energy_diagnostics(
    coefficients: ArrayLike,
    omega: float,
    order: int,
    length_1: float,
    properties_1: BeamProperties,
    length_2: float,
    properties_2: BeamProperties,
) -> EnergyDiagnostics:
    """Evaluate exact-quadrature component shares and the energy identity."""

    size = _order(order)
    vector = np.asarray(coefficients, dtype=float)
    if vector.shape != (6 * size,):
        raise ValueError(f"coefficients must have shape ({6 * size},)")
    frequency = _positive(omega, name="omega")
    nodes, weights = gauss_legendre_rule()
    kinetic = np.zeros(3, dtype=float)
    potential = np.zeros(3, dtype=float)
    for arm, length, properties in (
        (1, float(length_1), properties_1),
        (2, float(length_2), properties_2),
    ):
        fields = evaluate_arm_mode(vector, size, arm, nodes, length, properties)
        kinetic[0] += length * float(np.dot(weights, properties.m * fields.u**2))
        kinetic[1] += length * float(np.dot(weights, properties.m * fields.w**2))
        kinetic[2] += length * float(np.dot(weights, properties.J * fields.psi**2))
        potential[0] += length * float(np.dot(weights, properties.A * fields.u_prime**2))
        potential[1] += length * float(np.dot(weights, properties.D * fields.psi_prime**2))
        potential[2] += length * float(
            np.dot(weights, properties.S * (fields.w_prime + fields.psi) ** 2)
        )
    mass_total = float(np.sum(kinetic))
    stiffness_total = float(np.sum(potential))
    inertia = frequency**2 * mass_total
    identity = abs(stiffness_total - inertia) / max(
        abs(stiffness_total), abs(inertia), np.finfo(float).tiny
    )
    return EnergyDiagnostics(
        T_axial=float(kinetic[0]),
        T_transverse=float(kinetic[1]),
        T_rotation=float(kinetic[2]),
        U_axial=float(potential[0]),
        U_bending=float(potential[1]),
        U_shear=float(potential[2]),
        total_mass_norm=mass_total,
        total_stiffness_norm=stiffness_total,
        inertia_norm=inertia,
        energy_identity_relative=identity,
        T_axial_share=float(kinetic[0] / mass_total),
        T_transverse_share=float(kinetic[1] / mass_total),
        T_rotation_share=float(kinetic[2] / mass_total),
        U_axial_share=float(potential[0] / stiffness_total),
        U_bending_share=float(potential[1] / stiffness_total),
        U_shear_share=float(potential[2] / stiffness_total),
    )


def mass_inner_product(left: ArrayLike, right: ArrayLike, mass: ArrayLike) -> float:
    """Return the real coefficient-space physical mass inner product."""

    left_vector = np.asarray(left, dtype=float)
    right_vector = np.asarray(right, dtype=float)
    matrix = np.asarray(mass, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("mass must be square")
    if left_vector.shape != (matrix.shape[0],) or right_vector.shape != left_vector.shape:
        raise ValueError("vectors and mass dimensions are inconsistent")
    return float(left_vector @ matrix @ right_vector)


def mass_mac(left: ArrayLike, right: ArrayLike, mass: ArrayLike) -> float:
    """Return mass MAC for two nonzero vectors."""

    cross = mass_inner_product(left, right, mass)
    left_norm = mass_inner_product(left, left, mass)
    right_norm = mass_inner_product(right, right, mass)
    if left_norm <= 0.0 or right_norm <= 0.0:
        raise ValueError("mass norms must be positive")
    return float(cross * cross / (left_norm * right_norm))


def mass_orthonormalize(vectors: ArrayLike, mass: ArrayLike) -> FloatArray:
    """Return a mass-orthonormal basis spanning the supplied columns."""

    basis = np.asarray(vectors, dtype=float)
    matrix = np.asarray(mass, dtype=float)
    if basis.ndim == 1:
        basis = basis[:, None]
    if basis.ndim != 2 or matrix.shape != (basis.shape[0], basis.shape[0]):
        raise ValueError("vectors and mass dimensions are inconsistent")
    gram = 0.5 * (basis.T @ matrix @ basis + basis.T @ matrix.T @ basis)
    values, vectors_gram = np.linalg.eigh(gram)
    threshold = 100.0 * np.finfo(float).eps * max(float(np.max(values)), 1.0)
    if np.any(values <= threshold):
        raise ValueError("modal subspace columns must be linearly independent")
    transform = vectors_gram @ np.diag(1.0 / np.sqrt(values)) @ vectors_gram.T
    return _readonly(basis @ transform)


def compare_mass_subspaces(
    left_vectors: ArrayLike, right_vectors: ArrayLike, mass: ArrayLike
) -> SubspaceComparison:
    """Return singular overlaps and principal angles of equal-dimensional spaces."""

    matrix = np.asarray(mass, dtype=float)
    left = mass_orthonormalize(left_vectors, matrix)
    right = mass_orthonormalize(right_vectors, matrix)
    if left.shape[1] != right.shape[1]:
        raise ValueError("subspaces must have equal dimension")
    overlap = left.T @ matrix @ right
    singular_values = np.linalg.svd(overlap, compute_uv=False)
    singular_values = np.clip(singular_values, 0.0, 1.0)
    angles = np.arccos(singular_values)
    return SubspaceComparison(
        singular_values=_readonly(singular_values),
        principal_angles_rad=_readonly(angles),
        minimum_singular_value=float(np.min(singular_values)),
        maximum_principal_angle_rad=float(np.max(angles)),
    )


__all__ = [
    "BasisEvaluation",
    "COEFFICIENT_BLOCK_ORDER",
    "CONSTRAINT_RANK_RELATIVE_THRESHOLD",
    "ConstraintReduction",
    "CoupledRitzMatrices",
    "EnergyDiagnostics",
    "GAUSS_LEGENDRE_ORDER",
    "JointResiduals",
    "RITZ_GUARD_ORDER",
    "RITZ_ORDERS",
    "RitzSpectrum",
    "SubspaceComparison",
    "ArmModeFields",
    "assemble_coupled_ritz_matrices",
    "compare_mass_subspaces",
    "constraint_nullspace",
    "evaluate_arm_mode",
    "field_slice",
    "gauss_legendre_rule",
    "global_kinematic_map",
    "joint_residuals_from_ritz_mode",
    "kinematic_constraint_matrix",
    "mass_inner_product",
    "mass_mac",
    "mass_orthonormalize",
    "modal_energy_diagnostics",
    "shifted_legendre_clamped_basis",
    "solve_coupled_ritz_spectrum",
]
