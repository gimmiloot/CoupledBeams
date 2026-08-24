"""Reddy FSDT model of one straight symmetric laminated beam.

The module is deliberately independent of the historical coupled-rod models in
this repository.  It implements only a straight, undamped, symmetric laminate:

* classical-lamination integration of ``A``, ``B``, ``D`` and transverse-shear
  matrices;
* static elimination of unused laminate resultants to obtain the one-dimensional
  axial, bending and shear stiffnesses;
* the bending state ``[w, psi_b, Q, M]`` and the independent axial state
  ``[u, N]``;
* their block-diagonal union in the state ``[u, w, psi_b, N, Q, M]``.

Angles are measured from the global beam ``x`` axis to the lamina material-1
axis.  Engineering shear strains are used.  The transverse-shear ordering is
always ``(yz, xz)``; consequently the beam shear direction is index 1.  Ply
storage uses the project coordinate ``z_project`` increasing bottom-to-top,
whereas the cited Reddy beam pages use a downward-positive coordinate; thus
``z_project = -z_Reddy`` when comparing interface coordinates literally.

The shear-correction factor is never selected inside this module.  Callers must
pass ``K`` explicitly when reducing a laminate to beam properties.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Literal, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import simpson
from scipy.linalg import expm
from scipy.optimize import brentq, minimize_scalar


FloatArray = NDArray[np.float64]
ComplexArray = NDArray[np.complex128]
BendingBoundaryCondition = Literal[
    "hinged_hinged", "clamped_clamped", "clamped_free", "HH", "CC", "CF"
]
AxialBoundaryCondition = Literal["fixed_fixed", "fixed_free", "FF", "FC"]


_BENDING_BC_ALIASES = {
    "hinged_hinged": "hinged_hinged",
    "hh": "hinged_hinged",
    "clamped_clamped": "clamped_clamped",
    "cc": "clamped_clamped",
    "clamped_free": "clamped_free",
    "cf": "clamped_free",
}
_AXIAL_BC_ALIASES = {
    "fixed_fixed": "fixed_fixed",
    "ff": "fixed_fixed",
    "fixed_free": "fixed_free",
    "fc": "fixed_free",
}


def _positive(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and positive; got {value!r}.")
    return out


def _nonnegative(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out < 0.0:
        raise ValueError(f"{name} must be finite and nonnegative; got {value!r}.")
    return out


def _finite(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out):
        raise ValueError(f"{name} must be finite; got {value!r}.")
    return out


def _canonical_bending_bc(value: str) -> str:
    try:
        return _BENDING_BC_ALIASES[str(value).lower()]
    except KeyError as exc:
        raise ValueError(
            "bending boundary condition must be HH/hinged_hinged, "
            "CC/clamped_clamped, or CF/clamped_free."
        ) from exc


def _canonical_axial_bc(value: str) -> str:
    try:
        return _AXIAL_BC_ALIASES[str(value).lower()]
    except KeyError as exc:
        raise ValueError(
            "axial boundary condition must be FF/fixed_fixed or FC/fixed_free."
        ) from exc


@dataclass(frozen=True)
class OrthotropicLamina:
    """Plane-stress orthotropic lamina data in material axes.

    ``G23`` and ``G13`` correspond respectively to the local shear components
    ``(yz, xz)`` at zero ply angle.  Density is volumetric.
    """

    E1: float
    E2: float
    nu12: float
    G12: float
    G13: float
    G23: float
    rho: float
    name: str = ""

    def __post_init__(self) -> None:
        for field in ("E1", "E2", "G12", "G13", "G23", "rho"):
            object.__setattr__(self, field, _positive(getattr(self, field), field))
        nu12 = float(self.nu12)
        if not math.isfinite(nu12):
            raise ValueError("nu12 must be finite.")
        object.__setattr__(self, "nu12", nu12)
        if 1.0 - nu12 * self.nu21 <= 0.0:
            raise ValueError("The plane-stress denominator 1 - nu12*nu21 must be positive.")

    @property
    def nu21(self) -> float:
        return self.nu12 * self.E2 / self.E1

    @property
    def Q(self) -> FloatArray:
        return lamina_reduced_stiffness(self)

    @property
    def Q_shear(self) -> FloatArray:
        return lamina_transverse_shear_stiffness(self)


# A descriptive alias retained for callers that prefer ``Material`` terminology.
LaminaMaterial = OrthotropicLamina


@dataclass(frozen=True)
class Ply:
    """One constant-property ply; stacks are supplied bottom-to-top."""

    material: OrthotropicLamina
    angle_deg: float
    thickness: float
    label: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.material, OrthotropicLamina):
            raise TypeError("material must be an OrthotropicLamina.")
        object.__setattr__(self, "angle_deg", float(self.angle_deg))
        if not math.isfinite(self.angle_deg):
            raise ValueError("angle_deg must be finite.")
        object.__setattr__(self, "thickness", _positive(self.thickness, "thickness"))


@dataclass(frozen=True)
class LaminateSection:
    """Integrated laminate quantities per unit width.

    ``A``, ``B`` and ``D`` use the membrane ordering ``(xx, yy, xy)``.
    ``shear`` uses ``(yz, xz)``.  ``z_interfaces`` has one more entry than
    ``plies`` and runs from ``-h/2`` to ``+h/2``.
    """

    plies: tuple[Ply, ...]
    z_interfaces: FloatArray
    A: FloatArray
    B: FloatArray
    D: FloatArray
    shear: FloatArray
    I0: float
    I1: float
    I2: float

    @property
    def thickness(self) -> float:
        return float(self.z_interfaces[-1] - self.z_interfaces[0])

    @property
    def shear_matrix(self) -> FloatArray:
        return self.shear


@dataclass(frozen=True)
class SymmetryCheck:
    B_norm: float
    B_relative: float
    I1_absolute: float
    I1_relative: float
    tolerance: float
    is_symmetric: bool


@dataclass(frozen=True)
class EffectiveReduction:
    """Two algebraically equivalent static-condensation evaluations."""

    primary_index: int
    compliance_value: float
    schur_value: float
    absolute_difference: float
    relative_difference: float

    @property
    def value(self) -> float:
        return self.compliance_value


@dataclass(frozen=True)
class BeamProperties:
    """One-dimensional force/moment, mass and shear properties.

    ``A``, ``D`` and ``S`` are respectively axial, bending and corrected shear
    stiffnesses.  ``m`` and ``J`` are translational mass and rotary inertia per
    unit beam length.  ``K`` records the explicit shear-correction input.
    """

    A: float
    D: float
    S: float
    m: float
    J: float
    K: float
    width: float = 1.0
    axial_reduction: EffectiveReduction | None = None
    bending_reduction: EffectiveReduction | None = None
    shear_reduction_before_K: EffectiveReduction | None = None

    def __post_init__(self) -> None:
        for field in ("A", "D", "S", "m", "K", "width"):
            object.__setattr__(self, field, _positive(getattr(self, field), field))
        object.__setattr__(self, "J", _nonnegative(self.J, "J"))

    @property
    def axial_stiffness(self) -> float:
        return self.A

    @property
    def bending_stiffness(self) -> float:
        return self.D

    @property
    def shear_stiffness(self) -> float:
        return self.S

    @property
    def mass_per_length(self) -> float:
        return self.m

    @property
    def rotary_inertia_per_length(self) -> float:
        return self.J


@dataclass(frozen=True)
class HingedHingedMode:
    n: int
    branch: str
    wavenumber: float
    omega_squared: float
    omega: float
    psi_over_w: float


@dataclass(frozen=True)
class BendingDispersionMode:
    """One algebraic FSDT frequency branch at an arbitrary positive wavenumber."""

    branch: str
    wavenumber: float
    omega_squared: float
    omega: float
    psi_over_w: float


@dataclass(frozen=True)
class AxialMode:
    n: int
    boundary_condition: str
    wavenumber: float
    omega: float


@dataclass(frozen=True)
class RootDiagnostic:
    omega: float
    determinant: float
    sigma_min: float
    sigma_max: float
    sigma_ratio: float
    condition_number: float
    boundary_residual: float
    accepted: bool
    detection: str


@dataclass(frozen=True)
class ModalEnergies:
    modal_mass: float
    stiffness_integral: float
    inertia_integral: float
    strain_energy: float
    kinetic_energy: float
    identity_absolute_error: float
    identity_relative_error: float


@dataclass(frozen=True)
class ModeShape:
    x: FloatArray
    states: FloatArray
    normalization_factor: float
    boundary_residual: float
    energies: ModalEnergies


@dataclass(frozen=True)
class SpectrumMember:
    omega: float
    subsystem: str
    subsystem_index: int


@dataclass(frozen=True)
class SpectrumCluster:
    representative_omega: float
    members: tuple[SpectrumMember, ...]
    exact_degeneracy: bool

    @property
    def multiplicity(self) -> int:
        return len(self.members)


@dataclass(frozen=True)
class ReddyEq4351Parameters:
    """Auxiliary quantities for the source-secondary Eq. (4.3.51) check."""

    p: float
    q: float
    r: float
    lambda_value: complex
    mu_value: complex
    S11: complex
    S22: complex


def lamina_reduced_stiffness(material: OrthotropicLamina) -> FloatArray:
    """Return the plane-stress reduced stiffness ``Q`` in (11, 22, 12) order."""

    denominator = 1.0 - material.nu12 * material.nu21
    Q11 = material.E1 / denominator
    Q22 = material.E2 / denominator
    Q12 = material.nu12 * material.E2 / denominator
    return np.array(
        [[Q11, Q12, 0.0], [Q12, Q22, 0.0], [0.0, 0.0, material.G12]],
        dtype=float,
    )


def lamina_transverse_shear_stiffness(material: OrthotropicLamina) -> FloatArray:
    """Return local transverse shear stiffness in fixed ``(yz, xz)`` order."""

    return np.diag([material.G23, material.G13]).astype(float)


def transformed_reduced_stiffness(
    material: OrthotropicLamina, angle_deg: float
) -> FloatArray:
    """Return transformed membrane stiffness ``Qbar`` in (xx, yy, xy) order."""

    theta = math.radians(_finite(angle_deg, "angle_deg"))
    m = math.cos(theta)
    n = math.sin(theta)
    m2, n2 = m * m, n * n
    m4, n4 = m2 * m2, n2 * n2
    Q = lamina_reduced_stiffness(material)
    Q11, Q22, Q12, Q66 = Q[0, 0], Q[1, 1], Q[0, 1], Q[2, 2]
    mn2 = m2 * n2
    Qbar11 = Q11 * m4 + 2.0 * (Q12 + 2.0 * Q66) * mn2 + Q22 * n4
    Qbar22 = Q11 * n4 + 2.0 * (Q12 + 2.0 * Q66) * mn2 + Q22 * m4
    Qbar12 = (Q11 + Q22 - 4.0 * Q66) * mn2 + Q12 * (m4 + n4)
    Qbar16 = (Q11 - Q12 - 2.0 * Q66) * m2 * m * n - (
        Q22 - Q12 - 2.0 * Q66
    ) * m * n2 * n
    Qbar26 = (Q11 - Q12 - 2.0 * Q66) * m * n2 * n - (
        Q22 - Q12 - 2.0 * Q66
    ) * m2 * m * n
    Qbar66 = (Q11 + Q22 - 2.0 * Q12 - 2.0 * Q66) * mn2 + Q66 * (
        m4 + n4
    )
    return np.array(
        [
            [Qbar11, Qbar12, Qbar16],
            [Qbar12, Qbar22, Qbar26],
            [Qbar16, Qbar26, Qbar66],
        ],
        dtype=float,
    )


def transformed_transverse_shear_stiffness(
    material: OrthotropicLamina, angle_deg: float
) -> FloatArray:
    """Return transformed transverse shear stiffness in ``(yz, xz)`` order."""

    theta = math.radians(_finite(angle_deg, "angle_deg"))
    m = math.cos(theta)
    n = math.sin(theta)
    Q44 = material.G23
    Q55 = material.G13
    Qbar44 = Q44 * m * m + Q55 * n * n
    Qbar55 = Q44 * n * n + Q55 * m * m
    Qbar45 = (Q55 - Q44) * m * n
    return np.array([[Qbar44, Qbar45], [Qbar45, Qbar55]], dtype=float)


# Short source-notation aliases.
lamina_Q = lamina_reduced_stiffness
lamina_Qbar = transformed_reduced_stiffness
lamina_shear_Qbar = transformed_transverse_shear_stiffness


def integrate_laminate(plies: Sequence[Ply]) -> LaminateSection:
    """Integrate a bottom-to-top stack in ``z_project=-z_Reddy`` coordinates."""

    ply_tuple = tuple(plies)
    if not ply_tuple:
        raise ValueError("At least one ply is required.")
    if any(not isinstance(ply, Ply) for ply in ply_tuple):
        raise TypeError("Every stack entry must be a Ply.")
    thickness = math.fsum(ply.thickness for ply in ply_tuple)
    interfaces = [-0.5 * thickness]
    for ply in ply_tuple:
        interfaces.append(interfaces[-1] + ply.thickness)
    # Remove accumulated endpoint roundoff without changing internal interfaces.
    interfaces[-1] = 0.5 * thickness
    z = np.asarray(interfaces, dtype=float)
    A = np.zeros((3, 3), dtype=float)
    B = np.zeros((3, 3), dtype=float)
    D = np.zeros((3, 3), dtype=float)
    shear = np.zeros((2, 2), dtype=float)
    I0 = I1 = I2 = 0.0
    for index, ply in enumerate(ply_tuple):
        z0, z1 = float(z[index]), float(z[index + 1])
        Qbar = transformed_reduced_stiffness(ply.material, ply.angle_deg)
        Qbar_s = transformed_transverse_shear_stiffness(ply.material, ply.angle_deg)
        dz1 = z1 - z0
        dz2 = z1 * z1 - z0 * z0
        dz3 = z1**3 - z0**3
        A += Qbar * dz1
        B += 0.5 * Qbar * dz2
        D += (1.0 / 3.0) * Qbar * dz3
        shear += Qbar_s * dz1
        I0 += ply.material.rho * dz1
        I1 += 0.5 * ply.material.rho * dz2
        I2 += (1.0 / 3.0) * ply.material.rho * dz3
    return LaminateSection(
        plies=ply_tuple,
        z_interfaces=z,
        A=A,
        B=B,
        D=D,
        shear=shear,
        I0=float(I0),
        I1=float(I1),
        I2=float(I2),
    )


laminate_stiffness = integrate_laminate


def check_laminate_symmetry(
    laminate: LaminateSection, *, tolerance: float = 1.0e-10
) -> SymmetryCheck:
    """Check the symmetric-beam requirements ``B ~= 0`` and ``I1 ~= 0``."""

    tol = _positive(tolerance, "tolerance")
    h = laminate.thickness
    B_norm = float(np.linalg.norm(laminate.B, ord="fro"))
    B_scale = max(float(np.linalg.norm(laminate.A, ord="fro")) * h, np.finfo(float).tiny)
    I1_abs = abs(float(laminate.I1))
    I1_scale = max(abs(float(laminate.I0)) * h, np.finfo(float).tiny)
    B_relative = B_norm / B_scale
    I1_relative = I1_abs / I1_scale
    return SymmetryCheck(
        B_norm=B_norm,
        B_relative=B_relative,
        I1_absolute=I1_abs,
        I1_relative=I1_relative,
        tolerance=tol,
        is_symmetric=B_relative <= tol and I1_relative <= tol,
    )


def effective_stiffness_reduction(
    matrix: ArrayLike, primary_index: int
) -> EffectiveReduction:
    """Reduce a symmetric constitutive matrix by compliance and Schur routes."""

    array = np.asarray(matrix, dtype=float)
    if array.ndim != 2 or array.shape[0] != array.shape[1]:
        raise ValueError("matrix must be square.")
    if not np.all(np.isfinite(array)):
        raise ValueError("matrix must contain only finite values.")
    if not np.allclose(array, array.T, rtol=1.0e-11, atol=1.0e-13):
        raise ValueError("matrix must be symmetric.")
    n = array.shape[0]
    p = int(primary_index)
    if p < 0 or p >= n:
        raise IndexError("primary_index is outside the matrix.")
    unit = np.zeros(n, dtype=float)
    unit[p] = 1.0
    compliance_column = np.linalg.solve(array, unit)
    compliance = float(compliance_column[p])
    if compliance <= 0.0:
        raise ValueError("The selected compliance must be positive.")
    compliance_value = 1.0 / compliance
    remainder = [index for index in range(n) if index != p]
    if remainder:
        coupling = array[p, remainder]
        secondary = array[np.ix_(remainder, remainder)]
        schur_value = float(array[p, p] - coupling @ np.linalg.solve(secondary, coupling))
    else:
        schur_value = float(array[p, p])
    absolute = abs(compliance_value - schur_value)
    relative = absolute / max(abs(compliance_value), abs(schur_value), np.finfo(float).tiny)
    return EffectiveReduction(
        primary_index=p,
        compliance_value=float(compliance_value),
        schur_value=schur_value,
        absolute_difference=absolute,
        relative_difference=relative,
    )


def reduce_to_beam_properties(
    laminate: LaminateSection,
    *,
    width: float,
    K: float,
    symmetry_tolerance: float = 1.0e-10,
    reduction_tolerance: float = 1.0e-10,
) -> BeamProperties:
    """Reduce a symmetric laminate to the independent 1-D beam properties.

    Membrane and bending use component ``xx`` (index 0).  Shear uses component
    ``xz`` (index 1 of the fixed ``(yz, xz)`` ordering).  Secondary generalized
    resultants are statically eliminated.  ``K`` is applied only after the
    uncorrected laminate shear reduction.
    """

    b = _positive(width, "width")
    correction = _positive(K, "K")
    symmetry = check_laminate_symmetry(laminate, tolerance=symmetry_tolerance)
    if not symmetry.is_symmetric:
        raise ValueError(
            "Only symmetric laminates are implemented: "
            f"B_relative={symmetry.B_relative:.3e}, "
            f"I1_relative={symmetry.I1_relative:.3e}."
        )
    axial = effective_stiffness_reduction(laminate.A, 0)
    bending = effective_stiffness_reduction(laminate.D, 0)
    shear = effective_stiffness_reduction(laminate.shear, 1)
    for name, reduction in (("A", axial), ("D", bending), ("shear", shear)):
        if reduction.relative_difference > reduction_tolerance:
            raise ArithmeticError(
                f"{name} compliance/Schur reduction mismatch: "
                f"{reduction.relative_difference:.3e}."
            )
    return BeamProperties(
        A=b * axial.value,
        D=b * bending.value,
        S=correction * b * shear.value,
        m=b * laminate.I0,
        J=b * laminate.I2,
        K=correction,
        width=b,
        axial_reduction=axial,
        bending_reduction=bending,
        shear_reduction_before_K=shear,
    )


reduce_laminate_to_beam = reduce_to_beam_properties


def without_rotary_inertia(properties: BeamProperties) -> BeamProperties:
    """Return a diagnostic copy with ``J=0`` and all stiffnesses unchanged."""

    return BeamProperties(
        A=properties.A,
        D=properties.D,
        S=properties.S,
        m=properties.m,
        J=0.0,
        K=properties.K,
        width=properties.width,
        axial_reduction=properties.axial_reduction,
        bending_reduction=properties.bending_reduction,
        shear_reduction_before_K=properties.shear_reduction_before_K,
    )


def bending_state_matrix(omega: float, properties: BeamProperties) -> FloatArray:
    """Return the physical FSDT state matrix for ``[w, psi_b, Q, M]``.

    The frozen equations are

    ``w' = Q/S - psi_b``, ``psi_b' = M/D``,
    ``Q' = -m*omega**2*w``, and ``M' = Q - J*omega**2*psi_b``.
    """

    w = _nonnegative(omega, "omega")
    w2 = w * w
    return np.array(
        [
            [0.0, -1.0, 1.0 / properties.S, 0.0],
            [0.0, 0.0, 0.0, 1.0 / properties.D],
            [-properties.m * w2, 0.0, 0.0, 0.0],
            [0.0, -properties.J * w2, 1.0, 0.0],
        ],
        dtype=float,
    )


def bending_state_scale(properties: BeamProperties, length: float) -> FloatArray:
    """Return a dimensional scale for ``[w, psi_b, Q, M]``."""

    L = _positive(length, "length")
    return np.diag([L, 1.0, properties.D / L**2, properties.D / L]).astype(float)


def bending_transfer_matrix(
    omega: float, length: float, properties: BeamProperties
) -> FloatArray:
    """Propagate the physical bending state through ``length``."""

    L = _nonnegative(length, "length")
    return np.asarray(expm(bending_state_matrix(omega, properties) * L), dtype=float)


def bending_scaled_transfer_matrix(
    omega: float, length: float, properties: BeamProperties
) -> FloatArray:
    """Propagate the dimensionless state used by boundary diagnostics."""

    L = _positive(length, "length")
    scale = bending_state_scale(properties, L)
    inverse_scale = np.diag(1.0 / np.diag(scale))
    scaled_matrix = inverse_scale @ bending_state_matrix(omega, properties) @ scale
    return np.asarray(expm(scaled_matrix * L), dtype=float)


def bending_boundary_operators(
    boundary_condition: BendingBoundaryCondition | str,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return left selector, right selector, and left-admissible basis.

    The basis columns parameterize an initial dimensionless state satisfying the
    left conditions, so ``right @ transfer @ basis`` is the reduced 2x2
    characteristic matrix.
    """

    bc = _canonical_bending_bc(boundary_condition)
    if bc == "hinged_hinged":
        selector = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
        basis = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]])
        return selector.copy(), selector.copy(), basis
    clamp = np.array([[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]])
    basis = np.array([[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    if bc == "clamped_clamped":
        return clamp.copy(), clamp.copy(), basis
    free = np.array([[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]])
    return clamp.copy(), free, basis


def bending_boundary_matrix(
    omega: float,
    length: float,
    properties: BeamProperties,
    boundary_condition: BendingBoundaryCondition | str,
) -> FloatArray:
    """Return the independent 2x2 transfer-matrix characteristic matrix."""

    _left, right, basis = bending_boundary_operators(boundary_condition)
    return right @ bending_scaled_transfer_matrix(omega, length, properties) @ basis


def bending_full_boundary_matrix(
    omega: float,
    length: float,
    properties: BeamProperties,
    boundary_condition: BendingBoundaryCondition | str,
) -> FloatArray:
    """Return all four endpoint constraints acting on the initial scaled state."""

    left, right, _basis = bending_boundary_operators(boundary_condition)
    transfer = bending_scaled_transfer_matrix(omega, length, properties)
    return np.vstack([left, right @ transfer])


def equilibrated_matrix(matrix: ArrayLike) -> FloatArray:
    """Apply finite nonzero row/column scalings without changing singular roots."""

    out = np.asarray(matrix, dtype=float).copy()
    if out.ndim != 2:
        raise ValueError("matrix must be two-dimensional.")
    row_norms = np.linalg.norm(out, axis=1)
    nonzero_rows = row_norms > np.finfo(float).tiny
    out[nonzero_rows] /= row_norms[nonzero_rows, None]
    column_norms = np.linalg.norm(out, axis=0)
    nonzero_columns = column_norms > np.finfo(float).tiny
    out[:, nonzero_columns] /= column_norms[None, nonzero_columns]
    return out


def _row_normalized_matrix(matrix: ArrayLike) -> FloatArray:
    """Scale rows only, preserving right-nullspace coordinates."""

    out = np.asarray(matrix, dtype=float).copy()
    row_norms = np.linalg.norm(out, axis=1)
    nonzero_rows = row_norms > np.finfo(float).tiny
    out[nonzero_rows] /= row_norms[nonzero_rows, None]
    return out


def bending_characteristic_determinant(
    omega: float,
    length: float,
    properties: BeamProperties,
    boundary_condition: BendingBoundaryCondition | str,
) -> float:
    matrix = equilibrated_matrix(
        bending_boundary_matrix(omega, length, properties, boundary_condition)
    )
    return float(np.linalg.det(matrix))


def bending_root_diagnostic(
    omega: float,
    length: float,
    properties: BeamProperties,
    boundary_condition: BendingBoundaryCondition | str,
    *,
    sigma_ratio_tolerance: float = 1.0e-7,
    detection: str = "evaluation",
) -> RootDiagnostic:
    """Evaluate scaled determinant, singularity ratio, and endpoint residual."""

    matrix = equilibrated_matrix(
        bending_boundary_matrix(omega, length, properties, boundary_condition)
    )
    singular = np.linalg.svd(matrix, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    ratio = sigma_min / sigma_max if sigma_max > 0.0 else 0.0
    condition = sigma_max / sigma_min if sigma_min > 0.0 else math.inf
    _u, _s, vh = np.linalg.svd(matrix)
    vector = vh[-1]
    residual = float(np.linalg.norm(matrix @ vector))
    tolerance = _positive(sigma_ratio_tolerance, "sigma_ratio_tolerance")
    return RootDiagnostic(
        omega=float(omega),
        determinant=float(np.linalg.det(matrix)),
        sigma_min=sigma_min,
        sigma_max=sigma_max,
        sigma_ratio=ratio,
        condition_number=condition,
        boundary_residual=residual,
        accepted=ratio <= tolerance,
        detection=str(detection),
    )


def find_bending_roots(
    properties: BeamProperties,
    length: float,
    boundary_condition: BendingBoundaryCondition | str,
    *,
    omega_max: float,
    n_roots: int | None = None,
    omega_min: float = 0.0,
    scan_points: int = 4001,
    sigma_ratio_tolerance: float = 1.0e-7,
    root_xtol: float = 1.0e-11,
    dedup_rtol: float = 1.0e-8,
) -> tuple[RootDiagnostic, ...]:
    """Find positive transfer roots using determinant brackets plus SVD minima.

    The determinant supplies ordinary sign-changing roots.  Interior local
    minima of the smallest-singular-value ratio provide an independent route for
    even-multiplicity or nearly coincident candidates.  Every returned root must
    pass the requested full reduced-matrix SVD gate.
    """

    L = _positive(length, "length")
    lower = _nonnegative(omega_min, "omega_min")
    upper = _positive(omega_max, "omega_max")
    if upper <= lower:
        raise ValueError("omega_max must exceed omega_min.")
    count = None if n_roots is None else int(n_roots)
    if count is not None and count <= 0:
        raise ValueError("n_roots must be positive when supplied.")
    points = int(scan_points)
    if points < 33:
        raise ValueError("scan_points must be at least 33.")
    start = max(lower, np.finfo(float).eps * max(1.0, upper))
    grid = np.linspace(start, upper, points)

    def determinant(value: float) -> float:
        return bending_characteristic_determinant(
            value, L, properties, boundary_condition
        )

    def sigma_ratio(value: float) -> float:
        matrix = equilibrated_matrix(
            bending_boundary_matrix(value, L, properties, boundary_condition)
        )
        singular = np.linalg.svd(matrix, compute_uv=False)
        return float(singular[-1] / singular[0]) if singular[0] > 0.0 else 0.0

    determinants = np.asarray([determinant(float(value)) for value in grid])
    ratios = np.asarray([sigma_ratio(float(value)) for value in grid])
    candidates: list[tuple[float, str]] = []
    for index in range(points - 1):
        left, right = float(grid[index]), float(grid[index + 1])
        f_left, f_right = determinants[index], determinants[index + 1]
        if not (math.isfinite(float(f_left)) and math.isfinite(float(f_right))):
            continue
        if f_left == 0.0:
            candidates.append((left, "determinant_grid_zero"))
        elif f_left * f_right < 0.0:
            root = brentq(determinant, left, right, xtol=root_xtol, rtol=4.0 * np.finfo(float).eps)
            candidates.append((float(root), "determinant_bracket"))
    for index in range(1, points - 1):
        if ratios[index] <= ratios[index - 1] and ratios[index] <= ratios[index + 1]:
            fit = minimize_scalar(
                sigma_ratio,
                bounds=(float(grid[index - 1]), float(grid[index + 1])),
                method="bounded",
                options={"xatol": root_xtol},
            )
            if fit.success and float(fit.fun) <= sigma_ratio_tolerance:
                candidates.append((float(fit.x), "svd_minimum"))
    accepted: list[RootDiagnostic] = []
    for value, detection in sorted(candidates):
        diagnostic = bending_root_diagnostic(
            value,
            L,
            properties,
            boundary_condition,
            sigma_ratio_tolerance=sigma_ratio_tolerance,
            detection=detection,
        )
        if not diagnostic.accepted:
            continue
        duplicate = any(
            abs(diagnostic.omega - previous.omega)
            <= root_xtol + dedup_rtol * max(diagnostic.omega, previous.omega)
            for previous in accepted
        )
        if not duplicate:
            accepted.append(diagnostic)
    accepted.sort(key=lambda item: item.omega)
    if count is not None:
        accepted = accepted[:count]
    return tuple(accepted)


def bending_dispersion_branches(
    properties: BeamProperties, wavenumber: float
) -> tuple[BendingDispersionMode, ...]:
    """Return the FSDT algebraic branch(es) for any prescribed ``k > 0``.

    The roots are obtained from the frozen 2x2 harmonic system, equivalently

    ``m*J*x**2 - (J*S*k**2 + m*(S + D*k**2))*x + S*D*k**4 = 0``,

    where ``x = omega**2``.  This helper is also suitable for evaluating source
    classical-characteristic diagnostics at a separately supplied ``k=alpha/L``;
    it does not assign or validate a boundary condition by itself.
    """

    k = _positive(wavenumber, "wavenumber")
    S, D, mass, rotary = properties.S, properties.D, properties.m, properties.J
    if rotary == 0.0:
        omega_squared = S * D * k**4 / (mass * (S + D * k**2))
        ratio = (mass * omega_squared - S * k**2) / (S * k)
        return (
            BendingDispersionMode(
                branch="finite_no_RI",
                wavenumber=k,
                omega_squared=omega_squared,
                omega=math.sqrt(omega_squared),
                psi_over_w=ratio,
            ),
        )
    coefficient_a = mass * rotary
    middle = rotary * S * k**2 + mass * (S + D * k**2)
    coefficient_c = S * D * k**4
    discriminant = middle * middle - 4.0 * coefficient_a * coefficient_c
    scale = middle * middle + 4.0 * coefficient_a * coefficient_c
    if discriminant < -64.0 * np.finfo(float).eps * scale:
        raise ArithmeticError("The HH frequency quadratic has a negative discriminant.")
    radical = math.sqrt(max(0.0, discriminant))
    large = 0.5 * (middle + radical) / coefficient_a
    small = coefficient_c / (coefficient_a * large)
    values = sorted((small, large))
    modes = []
    for branch, omega_squared in zip(("lower", "upper"), values, strict=True):
        ratio = (mass * omega_squared - S * k**2) / (S * k)
        modes.append(
            BendingDispersionMode(
                branch=branch,
                wavenumber=k,
                omega_squared=omega_squared,
                omega=math.sqrt(omega_squared),
                psi_over_w=ratio,
            )
        )
    return tuple(modes)


fsdt_dispersion_branches = bending_dispersion_branches


def hinged_hinged_algebraic_modes(
    properties: BeamProperties, length: float, n: int
) -> tuple[HingedHingedMode, ...]:
    """Return the exact HH algebraic branch(es) for spatial index ``n``."""

    L = _positive(length, "length")
    mode_index = int(n)
    if mode_index <= 0 or mode_index != n:
        raise ValueError("n must be a positive integer.")
    branches = bending_dispersion_branches(properties, mode_index * math.pi / L)
    return tuple(
        HingedHingedMode(
            n=mode_index,
            branch=mode.branch,
            wavenumber=mode.wavenumber,
            omega_squared=mode.omega_squared,
            omega=mode.omega,
            psi_over_w=mode.psi_over_w,
        )
        for mode in branches
    )


def hinged_hinged_modes(
    properties: BeamProperties, length: float, n_modes: int = 3
) -> tuple[HingedHingedMode, ...]:
    """Return both exact algebraic branches for ``n=1..n_modes`` (RI case)."""

    count = int(n_modes)
    if count <= 0:
        raise ValueError("n_modes must be positive.")
    return tuple(
        mode
        for n in range(1, count + 1)
        for mode in hinged_hinged_algebraic_modes(properties, length, n)
    )


def hinged_hinged_exact_shape(
    mode: HingedHingedMode,
    properties: BeamProperties,
    length: float,
    x: ArrayLike,
    *,
    mass_normalize: bool = True,
) -> ModeShape:
    """Evaluate ``w=W sin(kx)``, ``psi_b=Psi cos(kx)`` and resultants."""

    L = _positive(length, "length")
    coordinates = _coordinate_array(x, L)
    k = mode.wavenumber
    ratio = mode.psi_over_w
    amplitude = 1.0
    if mass_normalize:
        modal_mass_per_amplitude = 0.5 * L * (properties.m + properties.J * ratio**2)
        amplitude = 1.0 / math.sqrt(modal_mass_per_amplitude)
    w = amplitude * np.sin(k * coordinates)
    psi = amplitude * ratio * np.cos(k * coordinates)
    shear_force = properties.S * amplitude * (k + ratio) * np.cos(k * coordinates)
    moment = -properties.D * k * amplitude * ratio * np.sin(k * coordinates)
    states = np.column_stack([w, psi, shear_force, moment])
    energies = bending_modal_energies(coordinates, states, properties, mode.omega)
    residual = max(
        abs(float(w[0])),
        abs(float(w[-1])),
        abs(float(moment[0])),
        abs(float(moment[-1])),
    )
    return ModeShape(
        x=coordinates,
        states=states,
        normalization_factor=amplitude,
        boundary_residual=residual,
        energies=energies,
    )


def reddy_eq_4_3_51_parameters(
    omega: float, length: float, properties: BeamProperties
) -> ReddyEq4351Parameters:
    """Build the printed-source auxiliary quantities for Eq. (4.3.51).

    This helper is a *secondary diagnostic only*.  Boundary matrices and roots
    must be built from the physical first-order state, not from the printed
    intermediate clamped-boundary Eq. (4.3.50a).
    """

    del length  # Kept in the signature to mirror the characteristic helper.
    frequency = _nonnegative(omega, "omega")
    frequency_squared = frequency * frequency
    p = properties.D
    q = (properties.D * properties.m / properties.S + properties.J) * frequency_squared
    r = (1.0 - properties.J * frequency_squared / properties.S) * properties.m * frequency_squared
    radical = np.lib.scimath.sqrt(q * q + 4.0 * p * r)
    lambda_value = complex(np.lib.scimath.sqrt((q + radical) / (2.0 * p)))
    mu_value = complex(np.lib.scimath.sqrt((-q + radical) / (2.0 * p)))
    S11 = mu_value * (properties.m * frequency_squared - lambda_value**2 * properties.S)
    S22 = lambda_value * (properties.m * frequency_squared + mu_value**2 * properties.S)
    return ReddyEq4351Parameters(
        p=p,
        q=q,
        r=r,
        lambda_value=lambda_value,
        mu_value=mu_value,
        S11=S11,
        S22=S22,
    )


def reddy_eq_4_3_51_characteristic(
    omega: float, length: float, properties: BeamProperties
) -> complex:
    """Evaluate Eq. (4.3.51) solely as a source-secondary CC diagnostic."""

    L = _positive(length, "length")
    data = reddy_eq_4_3_51_parameters(omega, L, properties)
    if abs(data.S11) == 0.0 or abs(data.S22) == 0.0:
        return complex(np.nan, np.nan)
    lamL = data.lambda_value * L
    muL = data.mu_value * L
    return complex(
        -2.0
        + 2.0 * np.cos(lamL) * np.cosh(muL)
        + np.sin(lamL)
        * np.sinh(muL)
        * (data.S22 / data.S11 - data.S11 / data.S22)
    )


def _coordinate_array(x: ArrayLike, length: float) -> FloatArray:
    coordinates = np.asarray(x, dtype=float)
    if coordinates.ndim != 1 or coordinates.size < 2:
        raise ValueError("x must be a one-dimensional array with at least two points.")
    if not np.all(np.isfinite(coordinates)) or np.any(np.diff(coordinates) <= 0.0):
        raise ValueError("x must be finite and strictly increasing.")
    tolerance = 64.0 * np.finfo(float).eps * max(1.0, length)
    if coordinates[0] < -tolerance or coordinates[-1] > length + tolerance:
        raise ValueError("x must lie inside [0, length].")
    return coordinates


def _integral(values: ArrayLike, x: FloatArray) -> float:
    values_array = np.asarray(values, dtype=float)
    return float(simpson(values_array, x=x))


def bending_modal_mass(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties
) -> float:
    coordinates = np.asarray(x, dtype=float)
    state_array = np.asarray(states)
    if state_array.shape != (coordinates.size, 4):
        raise ValueError("states must have shape (len(x), 4).")
    density = properties.m * np.abs(state_array[:, 0]) ** 2 + properties.J * np.abs(
        state_array[:, 1]
    ) ** 2
    return _integral(np.real(density), coordinates)


def normalize_bending_mode(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties
) -> tuple[FloatArray, float]:
    state_array = np.asarray(states, dtype=float)
    mass = bending_modal_mass(x, state_array, properties)
    if mass <= 0.0:
        raise ValueError("Bending modal mass must be positive.")
    factor = 1.0 / math.sqrt(mass)
    return state_array * factor, factor


def bending_modal_energies(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties, omega: float
) -> ModalEnergies:
    coordinates = np.asarray(x, dtype=float)
    state_array = np.asarray(states)
    if state_array.shape != (coordinates.size, 4):
        raise ValueError("states must have shape (len(x), 4).")
    modal_mass = bending_modal_mass(coordinates, state_array, properties)
    stiffness_density = np.abs(state_array[:, 2]) ** 2 / properties.S + np.abs(
        state_array[:, 3]
    ) ** 2 / properties.D
    stiffness = _integral(np.real(stiffness_density), coordinates)
    inertia = float(omega) ** 2 * modal_mass
    error = abs(stiffness - inertia)
    relative = error / max(abs(stiffness), abs(inertia), np.finfo(float).tiny)
    return ModalEnergies(
        modal_mass=modal_mass,
        stiffness_integral=stiffness,
        inertia_integral=inertia,
        strain_energy=0.5 * stiffness,
        kinetic_energy=0.5 * inertia,
        identity_absolute_error=error,
        identity_relative_error=relative,
    )


def bending_mode_shape(
    omega: float,
    properties: BeamProperties,
    length: float,
    boundary_condition: BendingBoundaryCondition | str,
    x: ArrayLike,
    *,
    mass_normalize: bool = True,
) -> ModeShape:
    """Recover a transfer-matrix bending mode and endpoint residual."""

    L = _positive(length, "length")
    coordinates = _coordinate_array(x, L)
    _left, right, basis = bending_boundary_operators(boundary_condition)
    characteristic = right @ bending_scaled_transfer_matrix(omega, L, properties) @ basis
    # Row scaling preserves the right null vector; column equilibration would
    # require an explicit back-transformation of the modal coefficients.
    _u, _s, vh = np.linalg.svd(_row_normalized_matrix(characteristic))
    coefficient = vh[-1]
    initial_scaled = basis @ coefficient
    significant = np.flatnonzero(np.abs(initial_scaled) > 64.0 * np.finfo(float).eps)
    if significant.size and initial_scaled[significant[0]] < 0.0:
        initial_scaled *= -1.0
    scale = bending_state_scale(properties, L)
    inverse_scale = np.diag(1.0 / np.diag(scale))
    scaled_system = inverse_scale @ bending_state_matrix(omega, properties) @ scale
    state_rows = np.vstack(
        [scale @ expm(scaled_system * float(position)) @ initial_scaled for position in coordinates]
    )
    displacement_scale = float(np.max(np.abs(state_rows[:, 0])))
    displacement_threshold = 64.0 * np.finfo(float).eps * displacement_scale
    significant_displacement = np.flatnonzero(
        np.abs(state_rows[:, 0]) > displacement_threshold
    )
    if (
        significant_displacement.size
        and state_rows[significant_displacement[0], 0] < 0.0
    ):
        state_rows *= -1.0
    factor = 1.0
    if mass_normalize:
        state_rows, factor = normalize_bending_mode(coordinates, state_rows, properties)
    left, right, _basis = bending_boundary_operators(boundary_condition)
    initial_hat = inverse_scale @ state_rows[0]
    final_hat = inverse_scale @ state_rows[-1]
    residual = float(
        max(np.linalg.norm(left @ initial_hat), np.linalg.norm(right @ final_hat))
    )
    energies = bending_modal_energies(coordinates, state_rows, properties, omega)
    return ModeShape(
        x=coordinates,
        states=np.asarray(state_rows, dtype=float),
        normalization_factor=factor,
        boundary_residual=residual,
        energies=energies,
    )


def axial_state_matrix(omega: float, properties: BeamProperties) -> FloatArray:
    """Return the independent symmetric-CLT axial state matrix for ``[u, N]``."""

    frequency = _nonnegative(omega, "omega")
    return np.array(
        [[0.0, 1.0 / properties.A], [-properties.m * frequency**2, 0.0]], dtype=float
    )


def axial_state_scale(properties: BeamProperties, length: float) -> FloatArray:
    _positive(length, "length")
    return np.diag([float(length), properties.A]).astype(float)


def axial_transfer_matrix(
    omega: float, length: float, properties: BeamProperties
) -> FloatArray:
    L = _nonnegative(length, "length")
    return np.asarray(expm(axial_state_matrix(omega, properties) * L), dtype=float)


def axial_scaled_transfer_matrix(
    omega: float, length: float, properties: BeamProperties
) -> FloatArray:
    L = _positive(length, "length")
    scale = axial_state_scale(properties, L)
    inverse_scale = np.diag(1.0 / np.diag(scale))
    return np.asarray(
        expm(inverse_scale @ axial_state_matrix(omega, properties) @ scale * L), dtype=float
    )


def axial_boundary_operators(
    boundary_condition: AxialBoundaryCondition | str,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    bc = _canonical_axial_bc(boundary_condition)
    fixed = np.array([[1.0, 0.0]])
    right = fixed.copy() if bc == "fixed_fixed" else np.array([[0.0, 1.0]])
    basis = np.array([[0.0], [1.0]])
    return fixed, right, basis


def axial_boundary_matrix(
    omega: float,
    length: float,
    properties: BeamProperties,
    boundary_condition: AxialBoundaryCondition | str,
) -> FloatArray:
    _left, right, basis = axial_boundary_operators(boundary_condition)
    return right @ axial_scaled_transfer_matrix(omega, length, properties) @ basis


def exact_axial_modes(
    properties: BeamProperties,
    length: float,
    boundary_condition: AxialBoundaryCondition | str,
    n_modes: int = 3,
) -> tuple[AxialMode, ...]:
    """Return exact roots of the fixed-fixed or fixed-free axial subsystem."""

    L = _positive(length, "length")
    bc = _canonical_axial_bc(boundary_condition)
    count = int(n_modes)
    if count <= 0:
        raise ValueError("n_modes must be positive.")
    wave_speed = math.sqrt(properties.A / properties.m)
    modes = []
    for n in range(1, count + 1):
        k = n * math.pi / L if bc == "fixed_fixed" else (n - 0.5) * math.pi / L
        modes.append(AxialMode(n=n, boundary_condition=bc, wavenumber=k, omega=k * wave_speed))
    return tuple(modes)


def axial_exact_shape(
    mode: AxialMode,
    properties: BeamProperties,
    length: float,
    x: ArrayLike,
    *,
    mass_normalize: bool = True,
) -> ModeShape:
    """Return the exact axial ``[u,N]`` shape in the generic ModeShape container."""

    L = _positive(length, "length")
    coordinates = _coordinate_array(x, L)
    amplitude = math.sqrt(2.0 / (properties.m * L)) if mass_normalize else 1.0
    displacement = amplitude * np.sin(mode.wavenumber * coordinates)
    force = properties.A * amplitude * mode.wavenumber * np.cos(mode.wavenumber * coordinates)
    states = np.column_stack([displacement, force])
    energies = axial_modal_energies(coordinates, states, properties, mode.omega)
    if mode.boundary_condition == "fixed_fixed":
        residual = max(abs(float(displacement[0])), abs(float(displacement[-1])))
    else:
        residual = max(abs(float(displacement[0])), abs(float(force[-1] / properties.A)))
    return ModeShape(
        x=coordinates,
        states=states,
        normalization_factor=amplitude,
        boundary_residual=residual,
        energies=energies,
    )


def axial_modal_mass(x: ArrayLike, states: ArrayLike, properties: BeamProperties) -> float:
    coordinates = np.asarray(x, dtype=float)
    state_array = np.asarray(states)
    if state_array.shape != (coordinates.size, 2):
        raise ValueError("states must have shape (len(x), 2).")
    return _integral(properties.m * np.abs(state_array[:, 0]) ** 2, coordinates)


def normalize_axial_mode(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties
) -> tuple[FloatArray, float]:
    state_array = np.asarray(states, dtype=float)
    mass = axial_modal_mass(x, state_array, properties)
    if mass <= 0.0:
        raise ValueError("Axial modal mass must be positive.")
    factor = 1.0 / math.sqrt(mass)
    return state_array * factor, factor


def axial_modal_energies(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties, omega: float
) -> ModalEnergies:
    coordinates = np.asarray(x, dtype=float)
    state_array = np.asarray(states)
    if state_array.shape != (coordinates.size, 2):
        raise ValueError("states must have shape (len(x), 2).")
    modal_mass = axial_modal_mass(coordinates, state_array, properties)
    stiffness = _integral(np.abs(state_array[:, 1]) ** 2 / properties.A, coordinates)
    inertia = float(omega) ** 2 * modal_mass
    error = abs(stiffness - inertia)
    relative = error / max(abs(stiffness), abs(inertia), np.finfo(float).tiny)
    return ModalEnergies(
        modal_mass=modal_mass,
        stiffness_integral=stiffness,
        inertia_integral=inertia,
        strain_energy=0.5 * stiffness,
        kinetic_energy=0.5 * inertia,
        identity_absolute_error=error,
        identity_relative_error=relative,
    )


def combined_state_matrix(omega: float, properties: BeamProperties) -> FloatArray:
    """Return the exact block-diagonal physics in interleaved combined ordering.

    The ordering is ``[u, w, psi_b, N, Q, M]``.  There are no axial--bending
    coupling entries.
    """

    frequency = _nonnegative(omega, "omega")
    w2 = frequency * frequency
    return np.array(
        [
            [0.0, 0.0, 0.0, 1.0 / properties.A, 0.0, 0.0],
            [0.0, 0.0, -1.0, 0.0, 1.0 / properties.S, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0 / properties.D],
            [-properties.m * w2, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, -properties.m * w2, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, -properties.J * w2, 0.0, 1.0, 0.0],
        ],
        dtype=float,
    )


def combined_state_scale(properties: BeamProperties, length: float) -> FloatArray:
    L = _positive(length, "length")
    return np.diag(
        [L, L, 1.0, properties.A, properties.D / L**2, properties.D / L]
    ).astype(float)


def combined_transfer_matrix(
    omega: float, length: float, properties: BeamProperties
) -> FloatArray:
    L = _nonnegative(length, "length")
    return np.asarray(expm(combined_state_matrix(omega, properties) * L), dtype=float)


def combined_scaled_transfer_matrix(
    omega: float, length: float, properties: BeamProperties
) -> FloatArray:
    L = _positive(length, "length")
    scale = combined_state_scale(properties, L)
    inverse_scale = np.diag(1.0 / np.diag(scale))
    return np.asarray(
        expm(inverse_scale @ combined_state_matrix(omega, properties) @ scale * L), dtype=float
    )


def embed_axial_state(state: ArrayLike) -> FloatArray:
    """Embed one axial ``[u,N]`` vector in combined state coordinates."""

    vector = np.asarray(state, dtype=float)
    if vector.shape != (2,):
        raise ValueError("axial state must have shape (2,).")
    out = np.zeros(6, dtype=float)
    out[[0, 3]] = vector
    return out


def embed_bending_state(state: ArrayLike) -> FloatArray:
    """Embed one bending ``[w,psi_b,Q,M]`` vector in combined coordinates."""

    vector = np.asarray(state, dtype=float)
    if vector.shape != (4,):
        raise ValueError("bending state must have shape (4,).")
    out = np.zeros(6, dtype=float)
    out[[1, 2, 4, 5]] = vector
    return out


def combined_boundary_operators(
    axial_boundary_condition: AxialBoundaryCondition | str,
    bending_boundary_condition: BendingBoundaryCondition | str,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return combined left/right selectors and the 3-D admissible basis."""

    axial_left, axial_right, axial_basis = axial_boundary_operators(axial_boundary_condition)
    bending_left, bending_right, bending_basis = bending_boundary_operators(
        bending_boundary_condition
    )
    left = np.zeros((3, 6), dtype=float)
    right = np.zeros((3, 6), dtype=float)
    left[0, [0, 3]] = axial_left[0]
    right[0, [0, 3]] = axial_right[0]
    left[np.ix_([1, 2], [1, 2, 4, 5])] = bending_left
    right[np.ix_([1, 2], [1, 2, 4, 5])] = bending_right
    basis = np.zeros((6, 3), dtype=float)
    basis[[0, 3], 0] = axial_basis[:, 0]
    basis[np.ix_([1, 2, 4, 5], [1, 2])] = bending_basis
    return left, right, basis


def combined_boundary_matrix(
    omega: float,
    length: float,
    properties: BeamProperties,
    axial_boundary_condition: AxialBoundaryCondition | str,
    bending_boundary_condition: BendingBoundaryCondition | str,
) -> FloatArray:
    """Return the reduced 3x3 block-diagonal characteristic matrix."""

    _left, right, basis = combined_boundary_operators(
        axial_boundary_condition, bending_boundary_condition
    )
    return right @ combined_scaled_transfer_matrix(omega, length, properties) @ basis


def combined_full_boundary_matrix(
    omega: float,
    length: float,
    properties: BeamProperties,
    axial_boundary_condition: AxialBoundaryCondition | str,
    bending_boundary_condition: BendingBoundaryCondition | str,
) -> FloatArray:
    """Return all six endpoint constraints on the initial combined state."""

    left, right, _basis = combined_boundary_operators(
        axial_boundary_condition, bending_boundary_condition
    )
    transfer = combined_scaled_transfer_matrix(omega, length, properties)
    return np.vstack([left, right @ transfer])


def embed_axial_shape(states: ArrayLike) -> FloatArray:
    """Embed sampled axial states of shape ``(n,2)`` in combined coordinates."""

    state_array = np.asarray(states, dtype=float)
    if state_array.ndim != 2 or state_array.shape[1] != 2:
        raise ValueError("axial states must have shape (n, 2).")
    out = np.zeros((state_array.shape[0], 6), dtype=float)
    out[:, [0, 3]] = state_array
    return out


def embed_bending_shape(states: ArrayLike) -> FloatArray:
    """Embed sampled bending states of shape ``(n,4)`` in combined coordinates."""

    state_array = np.asarray(states, dtype=float)
    if state_array.ndim != 2 or state_array.shape[1] != 4:
        raise ValueError("bending states must have shape (n, 4).")
    out = np.zeros((state_array.shape[0], 6), dtype=float)
    out[:, [1, 2, 4, 5]] = state_array
    return out


def combined_modal_mass(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties
) -> float:
    """Return mass norm for sampled ``[u,w,psi_b,N,Q,M]`` states."""

    coordinates = np.asarray(x, dtype=float)
    state_array = np.asarray(states)
    if state_array.shape != (coordinates.size, 6):
        raise ValueError("states must have shape (len(x), 6).")
    density = properties.m * (
        np.abs(state_array[:, 0]) ** 2 + np.abs(state_array[:, 1]) ** 2
    ) + properties.J * np.abs(state_array[:, 2]) ** 2
    return _integral(np.real(density), coordinates)


def normalize_combined_mode(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties
) -> tuple[FloatArray, float]:
    state_array = np.asarray(states, dtype=float)
    mass = combined_modal_mass(x, state_array, properties)
    if mass <= 0.0:
        raise ValueError("Combined modal mass must be positive.")
    factor = 1.0 / math.sqrt(mass)
    return state_array * factor, factor


def combined_modal_energies(
    x: ArrayLike, states: ArrayLike, properties: BeamProperties, omega: float
) -> ModalEnergies:
    """Evaluate the additive axial+bending energy identity."""

    coordinates = np.asarray(x, dtype=float)
    state_array = np.asarray(states)
    if state_array.shape != (coordinates.size, 6):
        raise ValueError("states must have shape (len(x), 6).")
    modal_mass = combined_modal_mass(coordinates, state_array, properties)
    stiffness_density = (
        np.abs(state_array[:, 3]) ** 2 / properties.A
        + np.abs(state_array[:, 4]) ** 2 / properties.S
        + np.abs(state_array[:, 5]) ** 2 / properties.D
    )
    stiffness = _integral(np.real(stiffness_density), coordinates)
    inertia = float(omega) ** 2 * modal_mass
    error = abs(stiffness - inertia)
    relative = error / max(abs(stiffness), abs(inertia), np.finfo(float).tiny)
    return ModalEnergies(
        modal_mass=modal_mass,
        stiffness_integral=stiffness,
        inertia_integral=inertia,
        strain_energy=0.5 * stiffness,
        kinetic_energy=0.5 * inertia,
        identity_absolute_error=error,
        identity_relative_error=relative,
    )


def union_subsystem_spectra(
    axial_omegas: Iterable[float],
    bending_omegas: Iterable[float],
    *,
    atol: float = 1.0e-10,
    rtol: float = 1.0e-8,
) -> tuple[SpectrumCluster, ...]:
    """Return sorted union clusters while retaining every multiplicity slot."""

    absolute_tolerance = _nonnegative(atol, "atol")
    relative_tolerance = _nonnegative(rtol, "rtol")
    members = [
        SpectrumMember(_nonnegative(value, "axial omega"), "axial", index)
        for index, value in enumerate(axial_omegas)
    ] + [
        SpectrumMember(_nonnegative(value, "bending omega"), "bending", index)
        for index, value in enumerate(bending_omegas)
    ]
    members.sort(key=lambda item: (item.omega, item.subsystem, item.subsystem_index))
    if not members:
        return ()
    grouped: list[list[SpectrumMember]] = [[members[0]]]
    for member in members[1:]:
        previous = grouped[-1][-1]
        threshold = absolute_tolerance + relative_tolerance * max(
            abs(member.omega), abs(previous.omega)
        )
        if abs(member.omega - previous.omega) <= threshold:
            grouped[-1].append(member)
        else:
            grouped.append([member])
    clusters = []
    exact_scale = 8.0 * np.finfo(float).eps
    for group in grouped:
        frequencies = [member.omega for member in group]
        representative = float(math.fsum(frequencies) / len(frequencies))
        exact = bool(
            max(frequencies) - min(frequencies)
            <= exact_scale * max(1.0, max(map(abs, frequencies)))
        )
        clusters.append(
            SpectrumCluster(
                representative_omega=representative,
                members=tuple(group),
                exact_degeneracy=exact and len(group) > 1,
            )
        )
    return tuple(clusters)


def combined_degeneracy_subspace(
    *,
    axial_vectors: ArrayLike | None = None,
    bending_vectors: ArrayLike | None = None,
    orthonormalize: bool = True,
) -> FloatArray:
    """Embed subsystem bases without assigning coupling meaning to SVD mixtures.

    Input columns are independent vectors in axial dimension 2 or bending
    dimension 4.  The returned columns span the same combined-state subspace.
    Any arbitrary vector selected from an exact/near-degenerate SVD nullspace is
    therefore treated as a coordinate choice, not as physical axial--bending
    coupling.
    """

    columns: list[FloatArray] = []
    if axial_vectors is not None:
        axial = np.asarray(axial_vectors, dtype=float)
        if axial.ndim == 1:
            axial = axial[:, None]
        if axial.ndim != 2 or axial.shape[0] != 2:
            raise ValueError("axial_vectors must have shape (2, n).")
        columns.extend(embed_axial_state(axial[:, index]) for index in range(axial.shape[1]))
    if bending_vectors is not None:
        bending = np.asarray(bending_vectors, dtype=float)
        if bending.ndim == 1:
            bending = bending[:, None]
        if bending.ndim != 2 or bending.shape[0] != 4:
            raise ValueError("bending_vectors must have shape (4, n).")
        columns.extend(
            embed_bending_state(bending[:, index]) for index in range(bending.shape[1])
        )
    if not columns:
        return np.empty((6, 0), dtype=float)
    basis = np.column_stack(columns)
    rank = int(np.linalg.matrix_rank(basis))
    if rank != basis.shape[1]:
        raise ValueError("Input subsystem vectors are linearly dependent.")
    if orthonormalize:
        basis, _r = np.linalg.qr(basis, mode="reduced")
    return np.asarray(basis, dtype=float)


# Intent-revealing alias used by combined-spectrum validation code.
combined_spectrum_union = union_subsystem_spectra


__all__ = [
    "AxialMode",
    "BendingDispersionMode",
    "BeamProperties",
    "EffectiveReduction",
    "HingedHingedMode",
    "LaminaMaterial",
    "LaminateSection",
    "ModalEnergies",
    "ModeShape",
    "OrthotropicLamina",
    "Ply",
    "ReddyEq4351Parameters",
    "RootDiagnostic",
    "SpectrumCluster",
    "SpectrumMember",
    "SymmetryCheck",
    "axial_boundary_matrix",
    "axial_boundary_operators",
    "axial_exact_shape",
    "axial_modal_energies",
    "axial_modal_mass",
    "axial_scaled_transfer_matrix",
    "axial_state_matrix",
    "axial_state_scale",
    "axial_transfer_matrix",
    "bending_boundary_matrix",
    "bending_boundary_operators",
    "bending_characteristic_determinant",
    "bending_dispersion_branches",
    "bending_full_boundary_matrix",
    "bending_modal_energies",
    "bending_modal_mass",
    "bending_mode_shape",
    "bending_root_diagnostic",
    "bending_scaled_transfer_matrix",
    "bending_state_matrix",
    "bending_state_scale",
    "bending_transfer_matrix",
    "check_laminate_symmetry",
    "combined_boundary_matrix",
    "combined_boundary_operators",
    "combined_full_boundary_matrix",
    "combined_modal_energies",
    "combined_modal_mass",
    "combined_spectrum_union",
    "combined_degeneracy_subspace",
    "combined_scaled_transfer_matrix",
    "combined_state_matrix",
    "combined_state_scale",
    "combined_transfer_matrix",
    "effective_stiffness_reduction",
    "embed_axial_state",
    "embed_axial_shape",
    "embed_bending_state",
    "embed_bending_shape",
    "equilibrated_matrix",
    "exact_axial_modes",
    "find_bending_roots",
    "fsdt_dispersion_branches",
    "hinged_hinged_algebraic_modes",
    "hinged_hinged_exact_shape",
    "hinged_hinged_modes",
    "integrate_laminate",
    "lamina_Q",
    "lamina_Qbar",
    "lamina_reduced_stiffness",
    "lamina_shear_Qbar",
    "lamina_transverse_shear_stiffness",
    "laminate_stiffness",
    "normalize_axial_mode",
    "normalize_bending_mode",
    "normalize_combined_mode",
    "reddy_eq_4_3_51_characteristic",
    "reddy_eq_4_3_51_parameters",
    "reduce_laminate_to_beam",
    "reduce_to_beam_properties",
    "transformed_reduced_stiffness",
    "transformed_transverse_shear_stiffness",
    "union_subsystem_spectra",
    "without_rotary_inertia",
]
