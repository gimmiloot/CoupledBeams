"""Independent closed-form isotropic Timoshenko two-arm comparator.

This narrowly scoped module adapts the project's verified legacy circular-beam
solution space to generic dimensional isotropic sections.  Its primary factory
is :func:`rectangular_section`, with ``I_y = b h**3 / 12`` for in-plane
bending.  :func:`circular_section` exists only to support numerical
backward-compatibility checks against the frozen circular implementation.

The implementation is intentionally independent of every Reddy, Yartsev,
ShellBuckling, and historical circular module.  It uses no transfer matrix or
matrix exponential.  Each clamped arm is represented by two exact Timoshenko
bending columns and one exact axial-sine column.

Legacy local-coordinate contract
--------------------------------

The first local coordinate runs from its outer clamp to the joint,
``x_1 in [0, L_1]``.  The second legacy coordinate has the opposite direction,
``x_2 in [0, -L_2]``.  The joint endpoint is consequently evaluated at
``x_2 = -L_2``.  The local shear convention is

``Q = KGA * (w' - psi)``.

The six coupled coefficients are ordered as
``[bend_1a, bend_1b, axial_1, bend_2a, bend_2b, axial_2]``.  The six boundary
rows retain the canonical legacy orientation exactly; no signed row fitting is
performed.  Diagnostic equilibration uses only finite positive diagonal row
and column factors, so it cannot move determinant zeros.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from numpy.typing import ArrayLike, NDArray


FloatArray = NDArray[np.float64]

TIMO_REGIME_MIXED = "mixed_hyperbolic_trigonometric"
TIMO_REGIME_CUTOFF = "cutoff_limit"
TIMO_REGIME_TWO_TRIG = "two_trigonometric_pairs"
TIMO_BASIS_COLUMN_SCALE_CAP = 1.0e3
MATRIX_SCALING_RELATIVE_FLOOR = float(np.finfo(float).eps ** 0.25)
TIMOSHENKO_BASIS_EVALUATOR_VERSION = (
    "generic_dimensional_closed_form_regime_complete_v1"
)
STATE_ORDER = ("u", "w", "psi", "N", "Q", "M")
LEGACY_SHEAR_CONVENTION = "Q=KGA*(dw_dx-psi)"
LEGACY_LOCAL_COORDINATE_CONTRACT = (
    "x1: outer_clamp_to_joint [0,+L1]; "
    "x2: outer_clamp_to_joint_in_parameter_order [0,-L2]"
)

COEFFICIENT_ORDER = (
    "arm1_bending_0",
    "arm1_bending_1",
    "arm1_axial",
    "arm2_bending_0",
    "arm2_bending_1",
    "arm2_axial",
)
BOUNDARY_ROW_ORDER = (
    "legacy_transverse_displacement_compatibility",
    "legacy_axial_displacement_compatibility",
    "section_rotation_compatibility",
    "bending_moment_equilibrium",
    "global_transverse_force_equilibrium",
    "global_axial_force_equilibrium",
)


def _finite(value: float, name: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite; got {value!r}.")
    return result


def _positive(value: float, name: str) -> float:
    result = _finite(value, name)
    if result <= 0.0:
        raise ValueError(f"{name} must be positive; got {value!r}.")
    return result


def _poisson_ratio(value: float) -> float:
    result = _finite(value, "nu")
    if not -1.0 < result < 0.5:
        raise ValueError("nu must lie in (-1, 0.5) for a stable isotropic solid.")
    return result


def _readonly(values: ArrayLike) -> FloatArray:
    result = np.array(values, dtype=float, copy=True)
    result.setflags(write=False)
    return result


def _finite_matrix(values: ArrayLike, *, name: str) -> FloatArray:
    result = np.asarray(values, dtype=float)
    if result.ndim != 2 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite two-dimensional matrix.")
    return np.array(result, dtype=float, copy=True)


@dataclass(frozen=True)
class SectionProperties:
    """Dimensional isotropic section data used directly by the comparator.

    ``second_moment`` is the bending second moment about the physical section
    rotation axis.  For the primary rectangular case this is ``I_y``.
    """

    E: float
    nu: float
    rho: float
    K: float
    area: float
    second_moment: float
    section_kind: str
    width: float | None = None
    thickness: float | None = None
    radius: float | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "E", _positive(self.E, "E"))
        object.__setattr__(self, "nu", _poisson_ratio(self.nu))
        object.__setattr__(self, "rho", _positive(self.rho, "rho"))
        object.__setattr__(self, "K", _positive(self.K, "K"))
        object.__setattr__(self, "area", _positive(self.area, "area"))
        object.__setattr__(
            self,
            "second_moment",
            _positive(self.second_moment, "second_moment"),
        )
        kind = str(self.section_kind)
        if not kind:
            raise ValueError("section_kind must not be empty.")
        object.__setattr__(self, "section_kind", kind)
        for field_name in ("width", "thickness", "radius"):
            value = getattr(self, field_name)
            if value is not None:
                object.__setattr__(self, field_name, _positive(value, field_name))

    @property
    def shear_modulus(self) -> float:
        return self.E / (2.0 * (1.0 + self.nu))

    @property
    def axial_stiffness(self) -> float:
        return self.E * self.area

    @property
    def bending_stiffness(self) -> float:
        return self.E * self.second_moment

    @property
    def shear_stiffness(self) -> float:
        return self.K * self.shear_modulus * self.area

    @property
    def mass_per_length(self) -> float:
        return self.rho * self.area

    @property
    def rotary_inertia_per_length(self) -> float:
        return self.rho * self.second_moment

    # Compact aliases used in manifests and analytic comparisons.
    @property
    def G(self) -> float:
        return self.shear_modulus

    @property
    def EA(self) -> float:
        return self.axial_stiffness

    @property
    def EI(self) -> float:
        return self.bending_stiffness

    @property
    def KGA(self) -> float:
        return self.shear_stiffness

    @property
    def rhoA(self) -> float:
        return self.mass_per_length

    @property
    def rhoI(self) -> float:
        return self.rotary_inertia_per_length

    @property
    def inertia(self) -> float:
        return self.second_moment

    @property
    def kappa(self) -> float:
        return self.K

    @property
    def I_y(self) -> float:
        return self.second_moment

    @property
    def b(self) -> float | None:
        return self.width

    @property
    def h(self) -> float | None:
        return self.thickness


def circular_shear_coefficient(nu: float) -> float:
    """Return the shear coefficient used by the frozen circular comparator."""

    value = _poisson_ratio(nu)
    return (6.0 + 12.0 * value + 6.0 * value**2) / (
        7.0 + 12.0 * value + 4.0 * value**2
    )


def rectangular_section(
    *,
    E: float,
    nu: float,
    rho: float,
    width: float,
    thickness: float,
    K: float,
) -> SectionProperties:
    """Build the primary rectangular section with ``area=b*h`` and ``I_y=b*h^3/12``."""

    b = _positive(width, "width")
    h = _positive(thickness, "thickness")
    return SectionProperties(
        E=E,
        nu=nu,
        rho=rho,
        K=K,
        area=b * h,
        second_moment=b * h**3 / 12.0,
        section_kind="rectangle_width_ey_thickness_ez",
        width=b,
        thickness=h,
    )


def circular_section(
    *,
    E: float,
    nu: float,
    rho: float,
    radius: float,
    K: float | None = None,
) -> SectionProperties:
    """Build a circular section solely for frozen-legacy backcompat controls."""

    r = _positive(radius, "radius")
    correction = circular_shear_coefficient(nu) if K is None else _positive(K, "K")
    return SectionProperties(
        E=E,
        nu=nu,
        rho=rho,
        K=correction,
        area=math.pi * r**2,
        second_moment=math.pi * r**4 / 4.0,
        section_kind="circle_backcompat_only",
        radius=r,
    )


@dataclass(frozen=True)
class TimoshenkoSpatialBasis:
    """Regime-complete spatial-root data for one clamped bending solution."""

    regime: str
    a: float
    b: float
    h: float
    q: float
    z_a: float
    z_b: float
    omega: float
    alpha: float
    warnings: tuple[str, ...]


@dataclass(frozen=True, eq=False)
class PositiveMatrixScaling:
    """Positive diagonal scaling with ``scaled=diag(r) raw diag(c)``."""

    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray

    @property
    def determinant_factor(self) -> float:
        return float(np.prod(self.row_factors) * np.prod(self.column_factors))


@dataclass(frozen=True, eq=False)
class BoundaryMatrixDiagnostics:
    """Raw/scaled determinant, SVD, and recovered-null-vector diagnostics."""

    raw_matrix: FloatArray
    scaled_matrix: FloatArray
    row_factors: FloatArray
    column_factors: FloatArray
    raw_determinant: float
    raw_det_sign: float
    raw_logabsdet: float
    scaled_determinant: float
    scaled_det_sign: float
    scaled_logabsdet: float
    raw_singular_values: FloatArray
    scaled_singular_values: FloatArray
    raw_sigma_min: float
    raw_sigma_max: float
    raw_sigma_ratio: float
    scaled_sigma_min: float
    scaled_sigma_max: float
    scaled_sigma_ratio: float
    raw_condition_number: float
    scaled_condition_number: float
    scaled_null_residual: float
    raw_boundary_null_residual: float
    raw_right_null_vector: FloatArray
    scaled_right_null_vector: FloatArray


@dataclass(frozen=True, eq=False)
class CoupledBoundaryAssembly:
    """Canonical legacy boundary matrix and both independent arm bases."""

    matrix: FloatArray
    basis_1: TimoshenkoSpatialBasis
    basis_2: TimoshenkoSpatialBasis
    warnings: tuple[str, ...]


@dataclass(frozen=True, eq=False)
class CoupledBoundaryEvaluation:
    assembly: CoupledBoundaryAssembly
    diagnostics: BoundaryMatrixDiagnostics


@dataclass(frozen=True, eq=False)
class ModeCoefficients:
    coefficients: FloatArray
    smallest_scaled_singular_value: float
    scaled_singular_ratio: float
    raw_boundary_residual: float
    warnings: tuple[str, ...]


@dataclass(frozen=True, eq=False)
class ArmFields:
    """Legacy local fields and exact derivative/resultant values on one arm."""

    x: FloatArray
    u: FloatArray
    w: FloatArray
    psi: FloatArray
    u_prime: FloatArray
    w_prime: FloatArray
    psi_prime: FloatArray
    u_second: FloatArray
    w_second: FloatArray
    psi_second: FloatArray
    N: FloatArray
    Q: FloatArray
    M: FloatArray
    N_prime: FloatArray
    Q_prime: FloatArray
    M_prime: FloatArray


@dataclass(frozen=True, eq=False)
class CoupledModeFields:
    """Two legacy-coordinate clamped-arm fields for one supplied frequency."""

    omega: float
    beta_deg: float
    coefficients: FloatArray
    section_1: SectionProperties
    section_2: SectionProperties
    basis_1: TimoshenkoSpatialBasis
    basis_2: TimoshenkoSpatialBasis
    arm_1: ArmFields
    arm_2: ArmFields
    boundary_diagnostics: BoundaryMatrixDiagnostics
    warnings: tuple[str, ...]


def omega_cutoff(section: SectionProperties) -> float:
    """Return the positive cutoff ``sqrt(KGA/rhoI)``."""

    return math.sqrt(section.shear_stiffness / section.rotary_inertia_per_length)


def axial_wavenumber(omega: float, section: SectionProperties) -> float:
    """Return ``omega*sqrt(rhoA/EA)`` for the independent axial equation."""

    frequency = _positive(omega, "omega")
    return frequency * math.sqrt(section.mass_per_length / section.axial_stiffness)


def timoshenko_spatial_basis(
    omega: float, section: SectionProperties
) -> TimoshenkoSpatialBasis:
    """Construct stable mixed, cutoff, or two-trigonometric bending roots.

    The squared spatial roots satisfy

    ``EI*z**2 + omega**2*(EI*rhoA/KGA + rhoI)*z``
    ``+ rhoI*rhoA*omega**4/KGA - rhoA*omega**2 = 0``.
    """

    frequency = _positive(omega, "omega")
    shear = section.shear_stiffness
    bending = section.bending_stiffness
    mass = section.mass_per_length
    rotary = section.rotary_inertia_per_length

    c2 = frequency**2 * (bending * mass / shear + rotary)
    c0 = rotary * mass * frequency**4 / shear - mass * frequency**2
    discriminant = c2**2 - 4.0 * bending * c0
    discriminant_scale = max(c2**2, abs(4.0 * bending * c0), np.finfo(float).tiny)
    if discriminant < -64.0 * np.finfo(float).eps * discriminant_scale:
        raise ArithmeticError(
            "Timoshenko spatial-root discriminant is materially negative at "
            f"omega={frequency:.17g}."
        )
    square_root = math.sqrt(max(0.0, discriminant))
    stable_product_root = -0.5 * (c2 + square_root)
    if stable_product_root == 0.0:
        raise ArithmeticError("The positive-frequency spatial-root product vanished.")
    z_b = stable_product_root / bending
    z_a = c0 / stable_product_root
    if z_a < z_b:
        z_a, z_b = z_b, z_a

    root_scale = max(1.0, abs(z_a), abs(z_b))
    zero_tolerance = 64.0 * np.finfo(float).eps * root_scale
    alpha = mass * frequency**2 / shear
    warnings: list[str] = []
    if abs(z_a) <= zero_tolerance:
        regime = TIMO_REGIME_CUTOFF
        a = 0.0
        b = math.sqrt(max(0.0, -z_b))
        h = alpha
        q = b - alpha / b
    elif z_a > 0.0 and z_b < 0.0:
        regime = TIMO_REGIME_MIXED
        a = math.sqrt(z_a)
        b = math.sqrt(-z_b)
        h = a + alpha / a
        q = b - alpha / b
    elif z_a < 0.0 and z_b < 0.0:
        regime = TIMO_REGIME_TWO_TRIG
        a = math.sqrt(-z_a)
        b = math.sqrt(-z_b)
        h = a - alpha / a
        q = b - alpha / b
    else:
        raise ArithmeticError(
            "Unsupported Timoshenko spatial-root regime: "
            f"z_a={z_a:.17g}, z_b={z_b:.17g}, omega={frequency:.17g}."
        )
    q_scale = max(1.0, abs(b), abs(alpha / b))
    if abs(q) <= 64.0 * np.finfo(float).eps * q_scale:
        warnings.append(
            "The second trigonometric coupling q is close to zero at "
            f"omega={frequency:.17g}."
        )
    return TimoshenkoSpatialBasis(
        regime=regime,
        a=float(a),
        b=float(b),
        h=float(h),
        q=float(q),
        z_a=float(z_a),
        z_b=float(z_b),
        omega=frequency,
        alpha=float(alpha),
        warnings=tuple(warnings),
    )


def bending_column_scales(basis: TimoshenkoSpatialBasis) -> FloatArray:
    """Return the frozen cutoff-continuous column orientation/scales."""

    if basis.regime == TIMO_REGIME_MIXED:
        scale_0 = basis.a**2 + basis.b**2
        raw_scale_1 = basis.a - (basis.h / basis.q) * basis.b
        scale_1 = math.copysign(
            min(abs(raw_scale_1), TIMO_BASIS_COLUMN_SCALE_CAP), raw_scale_1
        )
    elif basis.regime == TIMO_REGIME_CUTOFF:
        scale_0 = basis.b**2
        scale_1 = -TIMO_BASIS_COLUMN_SCALE_CAP
    elif basis.regime == TIMO_REGIME_TWO_TRIG:
        scale_0 = basis.b**2 - basis.a**2
        raw_scale_1 = basis.a - (basis.h / basis.q) * basis.b
        scale_1 = -min(abs(raw_scale_1), TIMO_BASIS_COLUMN_SCALE_CAP)
    else:
        raise ValueError(f"Unknown Timoshenko basis regime: {basis.regime!r}.")
    if not (
        math.isfinite(scale_0)
        and math.isfinite(scale_1)
        and scale_0 != 0.0
        and scale_1 != 0.0
    ):
        raise ArithmeticError(
            f"Singular bending-column scaling in regime {basis.regime!r}."
        )
    return np.array([scale_0, scale_1], dtype=float)


def clamped_bending_columns(
    x: float, basis: TimoshenkoSpatialBasis
) -> dict[str, FloatArray]:
    """Evaluate the two exact bending columns that satisfy ``w(0)=psi(0)=0``."""

    coordinate = _finite(x, "x")
    a, b, h, q = basis.a, basis.b, basis.h, basis.q
    bx = b * coordinate
    cos_bx = float(np.cos(bx))
    sin_bx = float(np.sin(bx))

    if basis.regime == TIMO_REGIME_CUTOFF:
        b2 = b**2
        one_minus_cos = 2.0 * float(np.sin(0.5 * bx)) ** 2
        canonical = {
            "w": np.array([one_minus_cos / b2, sin_bx / b]),
            "psi": np.array(
                [
                    (basis.alpha * coordinate + q * sin_bx) / b2,
                    -q * one_minus_cos / b,
                ]
            ),
            "w_prime": np.array([sin_bx / b, cos_bx]),
            "psi_prime": np.array(
                [
                    (basis.alpha + q * b * cos_bx) / b2,
                    -q * sin_bx,
                ]
            ),
        }
    elif basis.regime == TIMO_REGIME_MIXED:
        ax = a * coordinate
        cosh_ax = float(np.cosh(ax))
        sinh_ax = float(np.sinh(ax))
        cosh_minus_cos = (
            2.0 * float(np.sinh(0.5 * ax)) ** 2
            + 2.0 * float(np.sin(0.5 * bx)) ** 2
        )
        d0 = a**2 + b**2
        ratio = h / q
        d1 = a - ratio * b
        canonical = {
            "w": np.array(
                [
                    cosh_minus_cos / d0,
                    (sinh_ax - ratio * sin_bx) / d1,
                ]
            ),
            "psi": np.array(
                [
                    (h * sinh_ax + q * sin_bx) / d0,
                    h * cosh_minus_cos / d1,
                ]
            ),
            "w_prime": np.array(
                [
                    (a * sinh_ax + b * sin_bx) / d0,
                    (a * cosh_ax - ratio * b * cos_bx) / d1,
                ]
            ),
            "psi_prime": np.array(
                [
                    (h * a * cosh_ax + q * b * cos_bx) / d0,
                    h * (a * sinh_ax + b * sin_bx) / d1,
                ]
            ),
        }
    elif basis.regime == TIMO_REGIME_TWO_TRIG:
        ax = a * coordinate
        cos_ax = float(np.cos(ax))
        sin_ax = float(np.sin(ax))
        cos_difference = -2.0 * float(np.sin(0.5 * (ax + bx))) * float(
            np.sin(0.5 * (ax - bx))
        )
        d0 = b**2 - a**2
        ratio = h / q
        d1 = a - ratio * b
        canonical = {
            "w": np.array(
                [
                    cos_difference / d0,
                    (sin_ax - ratio * sin_bx) / d1,
                ]
            ),
            "psi": np.array(
                [
                    (-h * sin_ax + q * sin_bx) / d0,
                    h * cos_difference / d1,
                ]
            ),
            "w_prime": np.array(
                [
                    (-a * sin_ax + b * sin_bx) / d0,
                    (a * cos_ax - ratio * b * cos_bx) / d1,
                ]
            ),
            "psi_prime": np.array(
                [
                    (-h * a * cos_ax + q * b * cos_bx) / d0,
                    h * (-a * sin_ax + b * sin_bx) / d1,
                ]
            ),
        }
    else:
        raise ValueError(f"Unknown Timoshenko basis regime: {basis.regime!r}.")

    scales = bending_column_scales(basis)
    return {
        field: np.asarray(values * scales, dtype=float)
        for field, values in canonical.items()
    }


def clamped_endpoint_columns(
    x: float,
    omega: float,
    section: SectionProperties,
    *,
    basis: TimoshenkoSpatialBasis | None = None,
) -> dict[str, FloatArray]:
    """Return the three clamped endpoint columns in bending/bending/axial order."""

    coordinate = _finite(x, "x")
    frequency = _positive(omega, "omega")
    active_basis = (
        timoshenko_spatial_basis(frequency, section) if basis is None else basis
    )
    if not math.isclose(active_basis.omega, frequency, rel_tol=2.0e-14, abs_tol=0.0):
        raise ValueError("basis.omega is inconsistent with omega.")
    bending = clamped_bending_columns(coordinate, active_basis)
    fields = {
        name: np.zeros(3, dtype=float)
        for name in ("u", "w", "psi", "u_prime", "w_prime", "psi_prime")
    }
    for name in ("w", "psi", "w_prime", "psi_prime"):
        fields[name][:2] = bending[name]
    axial_k = axial_wavenumber(frequency, section)
    fields["u"][2] = math.sin(axial_k * coordinate)
    fields["u_prime"][2] = axial_k * math.cos(axial_k * coordinate)
    fields["N"] = section.axial_stiffness * fields["u_prime"]
    fields["Q"] = section.shear_stiffness * (
        fields["w_prime"] - fields["psi"]
    )
    fields["M"] = section.bending_stiffness * fields["psi_prime"]
    return fields


def _basis_warnings(*bases: TimoshenkoSpatialBasis) -> tuple[str, ...]:
    return tuple(dict.fromkeys(item for basis in bases for item in basis.warnings))


def assemble_legacy_coupled_boundary(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
) -> CoupledBoundaryAssembly:
    """Assemble the independent canonical legacy ``6 x 6`` boundary matrix.

    Arm 1 is evaluated at ``x_1=L_1`` and arm 2 at ``x_2=-L_2``.  ``beta_deg``
    is explicit at the public boundary; no ambiguous angle parameter is used.
    """

    frequency = _positive(omega, "omega")
    l1 = _positive(length_1, "length_1")
    l2 = _positive(length_2, "length_2")
    angle_deg = _finite(beta_deg, "beta_deg")
    basis_1 = timoshenko_spatial_basis(frequency, section_1)
    basis_2 = timoshenko_spatial_basis(frequency, section_2)
    arm_1 = clamped_endpoint_columns(l1, frequency, section_1, basis=basis_1)
    arm_2 = clamped_endpoint_columns(-l2, frequency, section_2, basis=basis_2)
    beta_rad = math.radians(angle_deg)
    cosine = math.cos(beta_rad)
    sine = math.sin(beta_rad)

    matrix = np.zeros((6, 6), dtype=float)
    matrix[0, 0:3] = arm_1["w"]
    matrix[0, 3:6] = -cosine * arm_2["w"] + sine * arm_2["u"]
    matrix[1, 0:3] = arm_1["u"]
    matrix[1, 3:6] = -sine * arm_2["w"] - cosine * arm_2["u"]
    matrix[2, 0:3] = arm_1["psi"]
    matrix[2, 3:6] = -arm_2["psi"]
    matrix[3, 0:3] = arm_1["M"]
    matrix[3, 3:6] = -arm_2["M"]
    matrix[4, 0:3] = -arm_1["Q"]
    matrix[4, 3:6] = cosine * arm_2["Q"] - sine * arm_2["N"]
    matrix[5, 0:3] = arm_1["N"]
    matrix[5, 3:6] = -sine * arm_2["Q"] - cosine * arm_2["N"]
    if not np.all(np.isfinite(matrix)):
        raise ArithmeticError("Legacy coupled boundary matrix is not finite.")
    return CoupledBoundaryAssembly(
        matrix=_readonly(matrix),
        basis_1=basis_1,
        basis_2=basis_2,
        warnings=_basis_warnings(basis_1, basis_2),
    )


def assemble_legacy_coupled_boundary_rad(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_rad: float,
) -> CoupledBoundaryAssembly:
    """Radian-explicit wrapper for :func:`assemble_legacy_coupled_boundary`."""

    angle_rad = _finite(beta_rad, "beta_rad")
    return assemble_legacy_coupled_boundary(
        omega,
        section_1,
        length_1,
        section_2,
        length_2,
        beta_deg=math.degrees(angle_rad),
    )


def legacy_coupled_boundary_matrix_raw(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
) -> FloatArray:
    """Return a copy of the canonical unscaled legacy boundary matrix."""

    assembly = assemble_legacy_coupled_boundary(
        omega,
        section_1,
        length_1,
        section_2,
        length_2,
        beta_deg=beta_deg,
    )
    return np.array(assembly.matrix, dtype=float, copy=True)


def positively_equilibrate_matrix(matrix: ArrayLike) -> PositiveMatrixScaling:
    """Apply only finite positive row/column factors to a matrix."""

    raw = _finite_matrix(matrix, name="matrix")
    row_norms = np.linalg.norm(raw, axis=1)
    row_factors = np.ones(raw.shape[0], dtype=float)
    row_reference = float(np.max(row_norms)) if row_norms.size else 0.0
    if row_reference > 0.0:
        row_floor = MATRIX_SCALING_RELATIVE_FLOOR * row_reference
        row_factors = 1.0 / np.maximum(row_norms, row_floor)
    row_scaled = row_factors[:, None] * raw

    column_norms = np.linalg.norm(row_scaled, axis=0)
    column_factors = 1.0 / np.maximum(column_norms, 1.0)
    scaled = row_scaled * column_factors[None, :]
    if not (
        np.all(np.isfinite(row_factors))
        and np.all(row_factors > 0.0)
        and np.all(np.isfinite(column_factors))
        and np.all(column_factors > 0.0)
        and np.all(np.isfinite(scaled))
    ):
        raise ArithmeticError("Equilibration factors and matrix must be finite and positive.")
    return PositiveMatrixScaling(
        scaled_matrix=_readonly(scaled),
        row_factors=_readonly(row_factors),
        column_factors=_readonly(column_factors),
    )


def boundary_matrix_diagnostics(matrix: ArrayLike) -> BoundaryMatrixDiagnostics:
    """Compute raw/scaled determinants, singular values, and null residuals."""

    raw = _finite_matrix(matrix, name="boundary matrix")
    if raw.shape[0] != raw.shape[1] or raw.shape[0] == 0:
        raise ValueError("boundary matrix must be nonempty and square.")
    scaling = positively_equilibrate_matrix(raw)
    scaled = np.asarray(scaling.scaled_matrix, dtype=float)
    raw_sign, raw_log = np.linalg.slogdet(raw)
    scaled_sign, scaled_log = np.linalg.slogdet(scaled)
    raw_singular = np.linalg.svd(raw, compute_uv=False)
    _left, scaled_singular, scaled_vh = np.linalg.svd(scaled, full_matrices=False)

    scaled_null = np.array(scaled_vh[-1, :], dtype=float, copy=True)
    raw_null = np.asarray(scaling.column_factors) * scaled_null
    raw_norm = float(np.linalg.norm(raw_null))
    if raw_norm == 0.0 or not math.isfinite(raw_norm):
        raise ArithmeticError("Recovered raw null vector is not finite and nonzero.")
    raw_null /= raw_norm
    scaled_norm = float(np.linalg.norm(scaled_null))
    scaled_null /= scaled_norm
    pivot = int(np.argmax(np.abs(raw_null)))
    if raw_null[pivot] < 0.0:
        raw_null = -raw_null
        scaled_null = -scaled_null

    raw_max = float(raw_singular[0])
    raw_min = float(raw_singular[-1])
    scaled_max = float(scaled_singular[0])
    scaled_min = float(scaled_singular[-1])
    scaled_residual = float(
        np.linalg.norm(scaled @ scaled_null)
        / max(scaled_max * np.linalg.norm(scaled_null), np.finfo(float).tiny)
    )
    raw_residual = float(
        np.linalg.norm(raw @ raw_null)
        / max(raw_max * np.linalg.norm(raw_null), np.finfo(float).tiny)
    )
    return BoundaryMatrixDiagnostics(
        raw_matrix=_readonly(raw),
        scaled_matrix=_readonly(scaled),
        row_factors=scaling.row_factors,
        column_factors=scaling.column_factors,
        raw_determinant=float(np.linalg.det(raw)),
        raw_det_sign=float(raw_sign),
        raw_logabsdet=float(raw_log),
        scaled_determinant=float(np.linalg.det(scaled)),
        scaled_det_sign=float(scaled_sign),
        scaled_logabsdet=float(scaled_log),
        raw_singular_values=_readonly(raw_singular),
        scaled_singular_values=_readonly(scaled_singular),
        raw_sigma_min=raw_min,
        raw_sigma_max=raw_max,
        raw_sigma_ratio=raw_min / raw_max if raw_max > 0.0 else 0.0,
        scaled_sigma_min=scaled_min,
        scaled_sigma_max=scaled_max,
        scaled_sigma_ratio=scaled_min / scaled_max if scaled_max > 0.0 else 0.0,
        raw_condition_number=raw_max / raw_min if raw_min > 0.0 else math.inf,
        scaled_condition_number=(
            scaled_max / scaled_min if scaled_min > 0.0 else math.inf
        ),
        scaled_null_residual=scaled_residual,
        raw_boundary_null_residual=raw_residual,
        raw_right_null_vector=_readonly(raw_null),
        scaled_right_null_vector=_readonly(scaled_null),
    )


def evaluate_legacy_coupled_boundary(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
) -> CoupledBoundaryEvaluation:
    """Return the canonical assembly together with full raw/scaled diagnostics."""

    assembly = assemble_legacy_coupled_boundary(
        omega,
        section_1,
        length_1,
        section_2,
        length_2,
        beta_deg=beta_deg,
    )
    return CoupledBoundaryEvaluation(
        assembly=assembly,
        diagnostics=boundary_matrix_diagnostics(assembly.matrix),
    )


def legacy_coupled_boundary_matrix_scaled(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
) -> FloatArray:
    """Return the positively equilibrated characteristic matrix."""

    evaluation = evaluate_legacy_coupled_boundary(
        omega,
        section_1,
        length_1,
        section_2,
        length_2,
        beta_deg=beta_deg,
    )
    return np.array(evaluation.diagnostics.scaled_matrix, dtype=float, copy=True)


def legacy_coupled_characteristic_determinant(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
    scaled: bool = True,
) -> float:
    """Return a raw or positively scaled determinant without searching roots."""

    raw = legacy_coupled_boundary_matrix_raw(
        omega,
        section_1,
        length_1,
        section_2,
        length_2,
        beta_deg=beta_deg,
    )
    if not scaled:
        return float(np.linalg.det(raw))
    scaled_matrix = positively_equilibrate_matrix(raw).scaled_matrix
    return float(np.linalg.det(scaled_matrix))


def coupled_mode_coefficients(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
) -> ModeCoefficients:
    """Return a deterministic SVD null candidate in canonical coefficient order."""

    evaluation = evaluate_legacy_coupled_boundary(
        omega,
        section_1,
        length_1,
        section_2,
        length_2,
        beta_deg=beta_deg,
    )
    diagnostics = evaluation.diagnostics
    return ModeCoefficients(
        coefficients=diagnostics.raw_right_null_vector,
        smallest_scaled_singular_value=diagnostics.scaled_sigma_min,
        scaled_singular_ratio=diagnostics.scaled_sigma_ratio,
        raw_boundary_residual=diagnostics.raw_boundary_null_residual,
        warnings=evaluation.assembly.warnings,
    )


def evaluate_clamped_arm_fields(
    omega: float,
    section: SectionProperties,
    basis: TimoshenkoSpatialBasis,
    coefficients: ArrayLike,
    x_values: ArrayLike,
) -> ArmFields:
    """Evaluate exact legacy local fields for one three-column clamped arm."""

    frequency = _positive(omega, "omega")
    if not math.isclose(basis.omega, frequency, rel_tol=2.0e-14, abs_tol=0.0):
        raise ValueError("basis.omega is inconsistent with omega.")
    coeff = np.asarray(coefficients, dtype=float)
    if coeff.shape != (3,) or not np.all(np.isfinite(coeff)):
        raise ValueError("coefficients must be a finite vector of length three.")
    coordinates = np.asarray(x_values, dtype=float)
    if coordinates.ndim != 1 or coordinates.size == 0 or not np.all(np.isfinite(coordinates)):
        raise ValueError("x_values must be a nonempty finite one-dimensional array.")

    bending_coeff = coeff[:2]
    w = np.empty_like(coordinates)
    psi = np.empty_like(coordinates)
    w_prime = np.empty_like(coordinates)
    psi_prime = np.empty_like(coordinates)
    for index, coordinate in enumerate(coordinates):
        columns = clamped_bending_columns(float(coordinate), basis)
        w[index] = float(columns["w"] @ bending_coeff)
        psi[index] = float(columns["psi"] @ bending_coeff)
        w_prime[index] = float(columns["w_prime"] @ bending_coeff)
        psi_prime[index] = float(columns["psi_prime"] @ bending_coeff)

    axial_k = axial_wavenumber(frequency, section)
    axial_amplitude = float(coeff[2])
    u = axial_amplitude * np.sin(axial_k * coordinates)
    u_prime = axial_amplitude * axial_k * np.cos(axial_k * coordinates)
    u_second = -(axial_k**2) * u
    normal = section.axial_stiffness * u_prime
    shear = section.shear_stiffness * (w_prime - psi)
    moment = section.bending_stiffness * psi_prime

    normal_prime = -section.mass_per_length * frequency**2 * u
    shear_prime = -section.mass_per_length * frequency**2 * w
    moment_prime = -shear - section.rotary_inertia_per_length * frequency**2 * psi
    w_second = psi_prime + shear_prime / section.shear_stiffness
    psi_second = moment_prime / section.bending_stiffness
    return ArmFields(
        x=_readonly(coordinates),
        u=_readonly(u),
        w=_readonly(w),
        psi=_readonly(psi),
        u_prime=_readonly(u_prime),
        w_prime=_readonly(w_prime),
        psi_prime=_readonly(psi_prime),
        u_second=_readonly(u_second),
        w_second=_readonly(w_second),
        psi_second=_readonly(psi_second),
        N=_readonly(normal),
        Q=_readonly(shear),
        M=_readonly(moment),
        N_prime=_readonly(normal_prime),
        Q_prime=_readonly(shear_prime),
        M_prime=_readonly(moment_prime),
    )


def legacy_coupled_mode_fields(
    omega: float,
    section_1: SectionProperties,
    length_1: float,
    section_2: SectionProperties,
    length_2: float,
    *,
    beta_deg: float,
    coefficients: ArrayLike | None = None,
    n_points: int = 801,
) -> CoupledModeFields:
    """Reconstruct both determinant-mode fields without a transfer or Ritz model."""

    frequency = _positive(omega, "omega")
    l1 = _positive(length_1, "length_1")
    l2 = _positive(length_2, "length_2")
    angle_deg = _finite(beta_deg, "beta_deg")
    points = int(n_points)
    if points < 2:
        raise ValueError("n_points must be at least two.")
    evaluation = evaluate_legacy_coupled_boundary(
        frequency,
        section_1,
        l1,
        section_2,
        l2,
        beta_deg=angle_deg,
    )
    if coefficients is None:
        coeff = np.asarray(evaluation.diagnostics.raw_right_null_vector, dtype=float)
    else:
        coeff = np.asarray(coefficients, dtype=float)
        if coeff.shape != (6,) or not np.all(np.isfinite(coeff)):
            raise ValueError("coefficients must be a finite vector of length six.")
        norm = float(np.linalg.norm(coeff))
        if norm == 0.0:
            raise ValueError("coefficients must be nonzero.")
        coeff = coeff / norm
        pivot = int(np.argmax(np.abs(coeff)))
        if coeff[pivot] < 0.0:
            coeff = -coeff

    x_1 = np.linspace(0.0, l1, points, dtype=float)
    x_2 = np.linspace(0.0, -l2, points, dtype=float)
    assembly = evaluation.assembly
    arm_1 = evaluate_clamped_arm_fields(
        frequency, section_1, assembly.basis_1, coeff[:3], x_1
    )
    arm_2 = evaluate_clamped_arm_fields(
        frequency, section_2, assembly.basis_2, coeff[3:], x_2
    )
    return CoupledModeFields(
        omega=frequency,
        beta_deg=angle_deg,
        coefficients=_readonly(coeff),
        section_1=section_1,
        section_2=section_2,
        basis_1=assembly.basis_1,
        basis_2=assembly.basis_2,
        arm_1=arm_1,
        arm_2=arm_2,
        boundary_diagnostics=evaluation.diagnostics,
        warnings=assembly.warnings,
    )


# Clear generic aliases for the validation CLI; every public angle remains explicit.
spatial_basis = timoshenko_spatial_basis
coupled_boundary_matrix_raw = legacy_coupled_boundary_matrix_raw
coupled_boundary_matrix_scaled = legacy_coupled_boundary_matrix_scaled
coupled_characteristic_determinant = legacy_coupled_characteristic_determinant
coupled_mode_fields = legacy_coupled_mode_fields


__all__ = [
    "ArmFields",
    "BOUNDARY_ROW_ORDER",
    "BoundaryMatrixDiagnostics",
    "COEFFICIENT_ORDER",
    "CoupledBoundaryAssembly",
    "CoupledBoundaryEvaluation",
    "CoupledModeFields",
    "LEGACY_LOCAL_COORDINATE_CONTRACT",
    "LEGACY_SHEAR_CONVENTION",
    "MATRIX_SCALING_RELATIVE_FLOOR",
    "ModeCoefficients",
    "PositiveMatrixScaling",
    "SectionProperties",
    "STATE_ORDER",
    "TIMO_BASIS_COLUMN_SCALE_CAP",
    "TIMO_REGIME_CUTOFF",
    "TIMO_REGIME_MIXED",
    "TIMO_REGIME_TWO_TRIG",
    "TIMOSHENKO_BASIS_EVALUATOR_VERSION",
    "TimoshenkoSpatialBasis",
    "assemble_legacy_coupled_boundary",
    "assemble_legacy_coupled_boundary_rad",
    "axial_wavenumber",
    "bending_column_scales",
    "boundary_matrix_diagnostics",
    "circular_section",
    "circular_shear_coefficient",
    "clamped_bending_columns",
    "clamped_endpoint_columns",
    "coupled_boundary_matrix_raw",
    "coupled_boundary_matrix_scaled",
    "coupled_characteristic_determinant",
    "coupled_mode_coefficients",
    "coupled_mode_fields",
    "evaluate_clamped_arm_fields",
    "evaluate_legacy_coupled_boundary",
    "legacy_coupled_boundary_matrix_raw",
    "legacy_coupled_boundary_matrix_scaled",
    "legacy_coupled_characteristic_determinant",
    "legacy_coupled_mode_fields",
    "omega_cutoff",
    "positively_equilibrate_matrix",
    "rectangular_section",
    "spatial_basis",
    "timoshenko_spatial_basis",
]
