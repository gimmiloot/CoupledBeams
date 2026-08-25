"""Narrow RLB helpers for two symmetric laminated beams with a rigid joint.

The physical state of each arm is ``[u, w, psi, N, Q, M]`` and its local
coordinate points from the outer clamp to the joint.  Geometry and physical
vector maps come exclusively from :mod:`reddy_inplane_geometry`; single-arm
state and transfer matrices come exclusively from
:mod:`reddy_symmetric_laminated_beam`.

This module contains no material reduction, laminate integration, root finder,
angle sweep, mode-shape plotting, or legacy coupled-rod implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import ArrayLike, NDArray

from scripts.lib import reddy_inplane_geometry as inplane_geometry
from scripts.lib import reddy_symmetric_laminated_beam as single_beam


FloatArray = NDArray[np.float64]
STATE_ORDER = inplane_geometry.STATE_ORDER
MATRIX_SCALING_RELATIVE_FLOOR = float(np.finfo(float).eps ** 0.25)


def _readonly(array: ArrayLike) -> FloatArray:
    result = np.array(array, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _finite_scalar(value: float, *, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _finite_vector(value: ArrayLike, *, size: int, name: str) -> FloatArray:
    result = np.asarray(value, dtype=float)
    if result.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},)")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return np.array(result, dtype=float, copy=True)


def _finite_matrix(value: ArrayLike, *, name: str) -> FloatArray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 2:
        raise ValueError(f"{name} must be two-dimensional")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return np.array(result, dtype=float, copy=True)


_clamp_basis = np.zeros((6, 3), dtype=float)
_clamp_basis[3:, :] = np.eye(3)
CLAMP_BASIS = _readonly(_clamp_basis)

_clamp_selector = np.zeros((3, 6), dtype=float)
_clamp_selector[:, :3] = np.eye(3)
CLAMP_SELECTOR = _readonly(_clamp_selector)

R0_BETA0 = _readonly(np.diag([-1.0, -1.0, 1.0, 1.0, 1.0, -1.0]))


@dataclass(frozen=True, eq=False)
class PhysicalJointResiduals:
    """Physical global residual vectors of the four invariant conditions."""

    displacement: FloatArray
    rotation: FloatArray
    force: FloatArray
    moment: FloatArray


@dataclass(frozen=True)
class VirtualWorkCheck:
    """Local/global virtual-work comparison for one prescribed variation."""

    local_total: float
    global_total: float
    pairing_absolute_difference: float
    absolute_residual: float
    scale: float
    normalized_residual: float


@dataclass(frozen=True, eq=False)
class PositiveMatrixScaling:
    """Positive diagonal row/column equilibration and its recorded factors.

    ``scaled_matrix = diag(row_factors) @ raw @ diag(column_factors)``.
    Every multiplier is finite and strictly positive.  Row-norm denominators
    are floored at ``eps**(1/4)`` times the largest row norm.  Column factors
    only downscale norms greater than one and never amplify a small column.
    This prevents a nearly zero characteristic row or column at a root from
    being inflated to unit norm.  If all rows are exactly zero, all
    multipliers are one.
    """

    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray

    @property
    def determinant_factor(self) -> float:
        return float(np.prod(self.row_factors) * np.prod(self.column_factors))


@dataclass(frozen=True, eq=False)
class BoundaryMatrixDiagnostics:
    """Raw/scaled characteristic-matrix diagnostics and null candidate."""

    raw_matrix: FloatArray
    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray
    raw_determinant: float
    scaled_determinant: float
    raw_sigma_min: float
    raw_sigma_max: float
    raw_sigma_ratio: float
    raw_singular_values: FloatArray
    scaled_sigma_min: float
    scaled_sigma_max: float
    scaled_sigma_ratio: float
    scaled_singular_values: FloatArray
    relative_singular_residual: float
    boundary_null_residual: float
    raw_right_null_vector: FloatArray
    scaled_right_null_vector: FloatArray

    @property
    def raw_condition_number(self) -> float:
        if self.raw_sigma_min == 0.0:
            return math.inf
        return self.raw_sigma_max / self.raw_sigma_min

    @property
    def scaled_condition_number(self) -> float:
        if self.scaled_sigma_min == 0.0:
            return math.inf
        return self.scaled_sigma_max / self.scaled_sigma_min


def _geometry_from_beta_rad(
    beta_rad: float,
) -> inplane_geometry.ReddyInPlaneGeometry:
    angle_rad = _finite_scalar(beta_rad, name="beta_rad")
    return inplane_geometry.reddy_inplane_geometry(math.degrees(angle_rad))


def physical_joint_residuals(
    beta_rad: float,
    arm1_state: ArrayLike,
    arm2_state: ArrayLike,
) -> PhysicalJointResiduals:
    """Return ``d1-d2``, ``vartheta1-vartheta2``, ``F1+F2``, ``m1+m2``."""

    geometry = _geometry_from_beta_rad(beta_rad)
    state_1 = _finite_vector(arm1_state, size=6, name="arm1_state")
    state_2 = _finite_vector(arm2_state, size=6, name="arm2_state")

    displacement_1 = inplane_geometry.displacement_vector(
        geometry.arm1, state_1[0], state_1[1]
    )
    displacement_2 = inplane_geometry.displacement_vector(
        geometry.arm2, state_2[0], state_2[1]
    )
    rotation_1 = inplane_geometry.rotation_vector(geometry.arm1, state_1[2])
    rotation_2 = inplane_geometry.rotation_vector(geometry.arm2, state_2[2])
    force_1 = inplane_geometry.force_vector(
        geometry.arm1, state_1[3], state_1[4]
    )
    force_2 = inplane_geometry.force_vector(
        geometry.arm2, state_2[3], state_2[4]
    )
    moment_1 = inplane_geometry.moment_vector(geometry.arm1, state_1[5])
    moment_2 = inplane_geometry.moment_vector(geometry.arm2, state_2[5])

    return PhysicalJointResiduals(
        displacement=_readonly(displacement_1 - displacement_2),
        rotation=_readonly(rotation_1 - rotation_2),
        force=_readonly(force_1 + force_2),
        moment=_readonly(moment_1 + moment_2),
    )


def scalar_joint_residuals_from_physical_maps(
    beta_rad: float,
    arm1_state: ArrayLike,
    arm2_state: ArrayLike,
) -> FloatArray:
    """Project physical residuals in the canonical ``t1,n1,e_y1`` rows."""

    geometry = _geometry_from_beta_rad(beta_rad)
    residuals = physical_joint_residuals(beta_rad, arm1_state, arm2_state)
    return _readonly(
        [
            np.dot(geometry.t1, residuals.displacement),
            np.dot(geometry.n1, residuals.displacement),
            np.dot(geometry.arm1.e_y, residuals.rotation),
            np.dot(geometry.t1, residuals.force),
            np.dot(geometry.n1, residuals.force),
            np.dot(geometry.arm1.e_y, residuals.moment),
        ]
    )


def joint_matrix_from_physical_maps(beta_rad: float) -> FloatArray:
    """Build the canonical 6x12 joint matrix from physical vector maps.

    No trigonometric coefficient is written in this builder.  Each column is
    obtained by applying one unit state to the invariant physical residuals.
    """

    angle_rad = _finite_scalar(beta_rad, name="beta_rad")
    matrix = np.empty((6, 12), dtype=float)
    for column in range(12):
        joint_state = np.zeros(12, dtype=float)
        joint_state[column] = 1.0
        matrix[:, column] = scalar_joint_residuals_from_physical_maps(
            angle_rad, joint_state[:6], joint_state[6:]
        )
    return _readonly(matrix)


def joint_matrix_closed_form(beta_rad: float) -> FloatArray:
    """Build the independently written canonical 6x12 ``c/s`` matrix."""

    angle_rad = _finite_scalar(beta_rad, name="beta_rad")
    cosine = math.cos(angle_rad)
    sine = math.sin(angle_rad)
    return _readonly(
        [
            [1, 0, 0, 0, 0, 0, cosine, sine, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, -sine, cosine, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, -cosine, -sine, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, sine, -cosine, 0],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1],
        ]
    )


def joint_matrix(beta_rad: float) -> FloatArray:
    """Return the production joint matrix using the invariant construction."""

    return joint_matrix_from_physical_maps(beta_rad)


def joint_virtual_work_check(
    beta_rad: float,
    resultants_1: ArrayLike,
    resultants_2: ArrayLike,
    virtual_motion_1: ArrayLike,
    virtual_motion_2: ArrayLike,
) -> VirtualWorkCheck:
    """Compare local and physical-global virtual work for both arms.

    Resultants are ordered ``[N,Q,M]`` and virtual motions are ordered
    ``[delta_u,delta_w,delta_psi]``.  Compatibility and equilibrium are not
    imposed inside this function, so callers can audit their construction.
    """

    geometry = _geometry_from_beta_rad(beta_rad)
    forces_1 = _finite_vector(resultants_1, size=3, name="resultants_1")
    forces_2 = _finite_vector(resultants_2, size=3, name="resultants_2")
    motion_1 = _finite_vector(virtual_motion_1, size=3, name="virtual_motion_1")
    motion_2 = _finite_vector(virtual_motion_2, size=3, name="virtual_motion_2")

    local_terms = np.concatenate((forces_1 * motion_1, forces_2 * motion_2))
    local_total = float(math.fsum(float(value) for value in local_terms))

    global_total = 0.0
    for basis, resultants, motion in (
        (geometry.arm1, forces_1, motion_1),
        (geometry.arm2, forces_2, motion_2),
    ):
        force = inplane_geometry.force_vector(basis, resultants[0], resultants[1])
        delta_displacement = inplane_geometry.displacement_vector(
            basis, motion[0], motion[1]
        )
        moment = inplane_geometry.moment_vector(basis, resultants[2])
        delta_rotation = inplane_geometry.rotation_vector(basis, motion[2])
        global_total += float(
            np.dot(force, delta_displacement) + np.dot(moment, delta_rotation)
        )

    scale = float(math.fsum(abs(float(value)) for value in local_terms))
    absolute_residual = abs(local_total)
    normalized_residual = 0.0 if scale == 0.0 else absolute_residual / scale
    return VirtualWorkCheck(
        local_total=local_total,
        global_total=global_total,
        pairing_absolute_difference=abs(local_total - global_total),
        absolute_residual=absolute_residual,
        scale=scale,
        normalized_residual=normalized_residual,
    )


def arm_state_matrix(
    omega: float, properties: single_beam.BeamProperties
) -> FloatArray:
    """Return the verified physical RLB state matrix of one arm."""

    return np.asarray(single_beam.combined_state_matrix(omega, properties), dtype=float)


def arm_transfer_matrix(
    omega: float,
    length: float,
    properties: single_beam.BeamProperties,
) -> FloatArray:
    """Return the verified physical transfer matrix ``expm(A*length)``."""

    return np.asarray(
        single_beam.combined_transfer_matrix(omega, length, properties), dtype=float
    )


def arm_clamp_map(
    omega: float,
    length: float,
    properties: single_beam.BeamProperties,
) -> FloatArray:
    """Map outer-clamp reactions ``[N(0),Q(0),M(0)]`` to the joint state."""

    return arm_transfer_matrix(omega, length, properties) @ CLAMP_BASIS


def coupled_clamp_to_joint_map(
    omega: float,
    length_1: float,
    properties_1: single_beam.BeamProperties,
    length_2: float,
    properties_2: single_beam.BeamProperties,
) -> FloatArray:
    """Return ``block_diag(H1,H2)`` in the canonical 12-state ordering."""

    map_1 = arm_clamp_map(omega, length_1, properties_1)
    map_2 = arm_clamp_map(omega, length_2, properties_2)
    combined = np.zeros((12, 6), dtype=float)
    combined[:6, :3] = map_1
    combined[6:, 3:] = map_2
    return combined


def coupled_boundary_matrix(
    omega: float,
    beta_rad: float,
    length_1: float,
    properties_1: single_beam.BeamProperties,
    length_2: float,
    properties_2: single_beam.BeamProperties,
) -> FloatArray:
    """Return the physical raw 6x6 coupled boundary matrix."""

    return joint_matrix(beta_rad) @ coupled_clamp_to_joint_map(
        omega, length_1, properties_1, length_2, properties_2
    )


def direct_fixed_fixed_boundary_matrix(
    omega: float,
    total_length: float,
    properties: single_beam.BeamProperties,
) -> FloatArray:
    """Independent one-transfer fixed--fixed reference matrix."""

    return (
        CLAMP_SELECTOR
        @ arm_transfer_matrix(omega, total_length, properties)
        @ CLAMP_BASIS
    )


def direct_stepped_boundary_matrix(
    omega: float,
    length_1: float,
    properties_1: single_beam.BeamProperties,
    length_2: float,
    properties_2: single_beam.BeamProperties,
) -> FloatArray:
    """Independent beta=0 direct reference for a material step.

    This builder deliberately does not call either joint-matrix builder, the
    coupled boundary matrix, or the coupled block assembly.
    """

    transfer_1 = arm_transfer_matrix(omega, length_1, properties_1)
    transfer_2 = arm_transfer_matrix(omega, length_2, properties_2)
    joint_state_from_reactions_1 = R0_BETA0 @ transfer_1 @ CLAMP_BASIS
    outer_state_2 = np.linalg.solve(transfer_2, joint_state_from_reactions_1)
    return CLAMP_SELECTOR @ outer_state_2


def positively_equilibrate_matrix(matrix: ArrayLike) -> PositiveMatrixScaling:
    """Apply continuous positive row/column scaling with a relative floor.

    A per-row factor ``1 / ||row||`` is unsuitable at a characteristic root
    whose complete row tends to zero: roundoff would inflate that row to unit
    norm and conceal the singularity.  The shared relative norm floor keeps
    the scaling finite and continuous while preserving determinant zeros.
    """

    raw = _finite_matrix(matrix, name="matrix")
    row_norms = np.linalg.norm(raw, axis=1)
    if not np.all(np.isfinite(row_norms)):
        raise ArithmeticError("row norms are not finite")
    row_factors = np.ones(raw.shape[0], dtype=float)
    row_reference = float(np.max(row_norms)) if row_norms.size else 0.0
    if row_reference > 0.0:
        row_floor = MATRIX_SCALING_RELATIVE_FLOOR * row_reference
        row_factors = 1.0 / np.maximum(row_norms, row_floor)
    row_scaled = row_factors[:, None] * raw

    column_norms = np.linalg.norm(row_scaled, axis=0)
    if not np.all(np.isfinite(column_norms)):
        raise ArithmeticError("column norms are not finite")
    # After row scaling, column multipliers are deliberately capped at one.
    # A second upscaling pass would otherwise re-inflate the near-zero row via
    # its near-zero column and conceal an axial or repeated root again.
    column_factors = 1.0 / np.maximum(column_norms, 1.0)
    scaled = row_scaled * column_factors[None, :]

    if not (
        np.all(np.isfinite(row_factors))
        and np.all(row_factors > 0.0)
        and np.all(np.isfinite(column_factors))
        and np.all(column_factors > 0.0)
    ):
        raise ArithmeticError("equilibration factors must be finite and positive")
    return PositiveMatrixScaling(
        scaled_matrix=_readonly(scaled),
        row_factors=_readonly(row_factors),
        column_factors=_readonly(column_factors),
    )


def boundary_matrix_diagnostics(matrix: ArrayLike) -> BoundaryMatrixDiagnostics:
    """Return determinant, singular and recovered-null-vector diagnostics."""

    raw = _finite_matrix(matrix, name="matrix")
    if raw.shape[0] != raw.shape[1]:
        raise ValueError("boundary matrix must be square")
    scaling = positively_equilibrate_matrix(raw)
    scaled = np.asarray(scaling.scaled_matrix)

    raw_singular_values = np.linalg.svd(raw, compute_uv=False)
    _left, scaled_singular_values, scaled_vh = np.linalg.svd(
        scaled, full_matrices=False
    )
    raw_sigma_max = float(raw_singular_values[0])
    raw_sigma_min = float(raw_singular_values[-1])
    scaled_sigma_max = float(scaled_singular_values[0])
    scaled_sigma_min = float(scaled_singular_values[-1])
    raw_sigma_ratio = (
        raw_sigma_min / raw_sigma_max if raw_sigma_max > 0.0 else 0.0
    )
    scaled_sigma_ratio = (
        scaled_sigma_min / scaled_sigma_max if scaled_sigma_max > 0.0 else 0.0
    )

    scaled_null = np.array(scaled_vh[-1, :], dtype=float, copy=True)
    raw_null = scaling.column_factors * scaled_null
    raw_null_norm = float(np.linalg.norm(raw_null))
    if raw_null_norm == 0.0 or not math.isfinite(raw_null_norm):
        raise ArithmeticError("recovered raw null vector is not finite and nonzero")
    raw_null /= raw_null_norm
    boundary_numerator = float(np.linalg.norm(raw @ raw_null))
    boundary_denominator = raw_sigma_max * float(np.linalg.norm(raw_null))
    boundary_null_residual = (
        boundary_numerator / boundary_denominator
        if boundary_denominator > 0.0
        else 0.0
    )

    return BoundaryMatrixDiagnostics(
        raw_matrix=_readonly(raw),
        scaled_matrix=_readonly(scaled),
        row_factors=scaling.row_factors,
        column_factors=scaling.column_factors,
        raw_determinant=float(np.linalg.det(raw)),
        scaled_determinant=float(np.linalg.det(scaled)),
        raw_sigma_min=raw_sigma_min,
        raw_sigma_max=raw_sigma_max,
        raw_sigma_ratio=raw_sigma_ratio,
        raw_singular_values=_readonly(raw_singular_values),
        scaled_sigma_min=scaled_sigma_min,
        scaled_sigma_max=scaled_sigma_max,
        scaled_sigma_ratio=scaled_sigma_ratio,
        scaled_singular_values=_readonly(scaled_singular_values),
        relative_singular_residual=scaled_sigma_ratio,
        boundary_null_residual=boundary_null_residual,
        raw_right_null_vector=_readonly(raw_null),
        scaled_right_null_vector=_readonly(scaled_null),
    )


__all__ = [
    "BoundaryMatrixDiagnostics",
    "CLAMP_BASIS",
    "CLAMP_SELECTOR",
    "MATRIX_SCALING_RELATIVE_FLOOR",
    "PhysicalJointResiduals",
    "PositiveMatrixScaling",
    "R0_BETA0",
    "STATE_ORDER",
    "VirtualWorkCheck",
    "arm_clamp_map",
    "arm_state_matrix",
    "arm_transfer_matrix",
    "boundary_matrix_diagnostics",
    "coupled_boundary_matrix",
    "coupled_clamp_to_joint_map",
    "direct_fixed_fixed_boundary_matrix",
    "direct_stepped_boundary_matrix",
    "joint_matrix",
    "joint_matrix_closed_form",
    "joint_matrix_from_physical_maps",
    "joint_virtual_work_check",
    "physical_joint_residuals",
    "positively_equilibrate_matrix",
    "scalar_joint_residuals_from_physical_maps",
]
