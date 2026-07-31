"""Diagnostic reproduction model for Yartsev Chapter 2, one monoclinic rod.

This module is intentionally script-owned.  It implements only the monoclinic
Timoshenko bending--torsion model used for the Chapter-2 source-reproduction
gate; it is not a production anisotropic API and contains no coupled-rod,
Euler--Bernoulli, Saint-Venant, or FEM model.

Source conventions
------------------
* Complex moduli follow Chapter 1, equations (1.32) and (1.34):
  ``M* = Re(M*) * (1 + 1j * eta)``.
* Material rotation follows equation (1.52), and ``Gxy_bar`` follows (1.56).
* Generalized torsional stiffness follows equations (2.8)--(2.10).
* ``state_corrected`` uses equations (2.1)--(2.10) with the dimensionally
  necessary ``I_y`` restored in (2.1).
* ``eliminated_printed`` preserves the printed signs of ``d0`` and ``f0``
  after (2.16).  ``eliminated_corrected`` uses the internally consistent
  positive signs documented by the source audit.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Literal, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import expm
from scipy.optimize import brentq, linear_sum_assignment, minimize_scalar, root


MaterialMode = Literal["elastic", "book_complex"]
Formulation = Literal[
    "state_corrected", "eliminated_corrected", "eliminated_printed"
]
ClampVariant = Literal["book_slope_clamp", "timoshenko_section_clamp"]
BoundaryMatrixBuilder = Callable[
    [complex, "RodPoint", Formulation], NDArray[np.complex128]
]


@dataclass(frozen=True)
class Geometry:
    """Book specimen geometry, in SI units."""

    a: float = 9.76e-3
    b: float = 2.524e-2
    length: float = 0.295
    shear_factor: float = 5.0 / 6.0

    @property
    def area(self) -> float:
        return self.a * self.b

    @property
    def I_y(self) -> float:
        return self.a**3 * self.b / 12.0

    @property
    def I_p(self) -> float:
        return self.a * self.b * (self.a**2 + self.b**2) / 12.0


@dataclass(frozen=True)
class BookMaterial:
    """T-53(VM)-78/PN-609-21M data from Chapter-1 Table 1.2."""

    E1_real: float = 28.4e9
    E2_real: float = 11.6e9
    G12_real: float = 4.4e9
    G13_real: float = 3.7e9
    G23_real: float = 2.6e9
    nu12: float = 0.22
    rho: float = 1650.0
    eta1: float = 4.3e-3
    eta2: float = 13.9e-3
    eta12: float = 2.0e-2
    eta13: float = 2.74e-2
    eta23: float = 3.1e-2
    ply_thickness: float = 3.5e-4

    def moduli(
        self, mode: MaterialMode, *, loss_scale: float = 1.0
    ) -> tuple[complex, complex, complex, complex, complex]:
        if mode not in ("elastic", "book_complex"):
            raise ValueError(f"unsupported material mode: {mode}")
        if not 0.0 <= loss_scale <= 1.0:
            raise ValueError("loss_scale must lie in [0, 1]")
        scale = 0.0 if mode == "elastic" else float(loss_scale)
        return (
            self.E1_real * (1.0 + 1j * scale * self.eta1),
            self.E2_real * (1.0 + 1j * scale * self.eta2),
            self.G12_real * (1.0 + 1j * scale * self.eta12),
            self.G13_real * (1.0 + 1j * scale * self.eta13),
            self.G23_real * (1.0 + 1j * scale * self.eta23),
        )


def hms_dx_209_material() -> BookMaterial:
    """HMS/DX-209 data from Chapter-1 Table 1.1 (printed p.45)."""

    return BookMaterial(
        E1_real=191.0e9,
        E2_real=5.0e9,
        G12_real=3.0e9,
        G13_real=3.0e9,
        G23_real=2.5e9,
        nu12=0.279,
        rho=1580.0,
        eta1=7.8e-4,
        eta2=6.7e-3,
        eta12=1.16e-2,
        eta13=1.16e-2,
        eta23=1.15e-2,
        ply_thickness=2.0e-4,
    )


def cantilever_geometry(length: float = 0.2) -> Geometry:
    """Rectangular specimen geometry used in Chapter 2 section 2.2."""

    if length <= 0.0:
        raise ValueError("cantilever length must be positive")
    return Geometry(a=5.0e-3, b=2.0e-2, length=float(length))


@dataclass(frozen=True)
class RotatedProperties:
    theta_deg: float
    S11: complex
    S22: complex
    S12: complex
    S66: complex
    S44: complex
    S55: complex
    Sbar11: complex
    Sbar16: complex
    Sbar66: complex
    Sbar55: complex
    Ex: complex
    Gxy: complex
    Gxz: complex
    Gxy_bar: complex
    mu_x_xy: complex
    mu_xy_x: complex


@dataclass(frozen=True)
class TorsionalStiffness:
    Cbar: complex
    C_T: complex
    series_sum: complex
    terms_used: int
    estimated_relative_tail: float


@dataclass(frozen=True)
class RodPoint:
    geometry: Geometry
    material: BookMaterial
    material_mode: MaterialMode
    loss_scale: float
    properties: RotatedProperties
    torsion: TorsionalStiffness


@dataclass(frozen=True)
class BoundaryQuality:
    determinant: complex
    sigma_min: float
    sigma_max: float
    relative_singular_residual: float


@dataclass(frozen=True)
class RootResult:
    omega: complex
    frequency_hz: float
    determinant_residual: float
    raw_determinant_abs: float
    sigma_min: float
    sigma_max: float
    relative_singular_residual: float
    min_neighbor_distance_hz: float = float("nan")
    refinements: int = 0
    status: str = "accepted"


def rotate_material(
    theta_deg: float,
    *,
    material: BookMaterial | None = None,
    material_mode: MaterialMode = "elastic",
    loss_scale: float = 1.0,
) -> RotatedProperties:
    """Rotate the Chapter-1 compliance entries exactly as in equation (1.52)."""

    mat = material or BookMaterial()
    E1, E2, G12, G13, G23 = mat.moduli(material_mode, loss_scale=loss_scale)
    theta = np.deg2rad(float(theta_deg))
    m = float(np.cos(theta))
    n = float(np.sin(theta))

    S11 = 1.0 / E1
    S22 = 1.0 / E2
    S12 = -mat.nu12 / E1
    S66 = 1.0 / G12
    S44 = 1.0 / G23
    S55 = 1.0 / G13

    Sbar11 = (
        m**4 * S11
        + 2.0 * m**2 * n**2 * S12
        + n**4 * S22
        + m**2 * n**2 * S66
    )
    Sbar16 = (
        2.0 * m**3 * n * S11
        + 2.0 * m * n * (n**2 - m**2) * S12
        - 2.0 * m * n**3 * S22
        + m * n * (n**2 - m**2) * S66
    )
    Sbar66 = (
        4.0 * m**2 * n**2 * S11
        - 8.0 * m**2 * n**2 * S12
        + 4.0 * m**2 * n**2 * S22
        + (m**2 - n**2) ** 2 * S66
    )
    Sbar55 = n**2 * S44 + m**2 * S55

    Ex = 1.0 / Sbar11
    Gxy = 1.0 / Sbar66
    Gxz = 1.0 / Sbar55
    Gxy_bar = 1.0 / (Sbar66 - Sbar16**2 / Sbar11)
    mu_x_xy = Gxy * Sbar16
    mu_xy_x = Ex * Sbar16

    return RotatedProperties(
        theta_deg=float(theta_deg),
        S11=S11,
        S22=S22,
        S12=S12,
        S66=S66,
        S44=S44,
        S55=S55,
        Sbar11=Sbar11,
        Sbar16=Sbar16,
        Sbar66=Sbar66,
        Sbar55=Sbar55,
        Ex=Ex,
        Gxy=Gxy,
        Gxz=Gxz,
        Gxy_bar=Gxy_bar,
        mu_x_xy=mu_x_xy,
        mu_xy_x=mu_xy_x,
    )


def _odd_fifth_series(
    base_argument: complex,
    *,
    relative_tolerance: float,
    max_terms: int,
) -> tuple[complex, int, float]:
    if relative_tolerance <= 0.0:
        raise ValueError("relative_tolerance must be positive")
    if max_terms < 1:
        raise ValueError("max_terms must be positive")
    if np.real(base_argument) <= 0.0:
        raise ValueError("torsional-series argument must have positive real part")

    total = 0.0j
    # For Re(z)>0, |tanh(z)| is bounded by coth(Re(z)).  Combined with
    # the integral bound for the odd inverse-fifth-power tail, this gives a
    # conservative stopping test for the weakly complex book material.
    tanh_bound = 1.0 / np.tanh(float(np.real(base_argument)))
    relative_tail = float("inf")
    for q in range(max_terms):
        odd = 2 * q + 1
        argument = odd * base_argument
        # Direct complex tanh overflows internally for very large positive
        # real parts although its value is already unity to machine precision.
        tanh_value = 1.0 + 0.0j if np.real(argument) > 20.0 else np.tanh(argument)
        total += tanh_value / odd**5
        next_odd = 2 * (q + 1) + 1
        tail_bound = tanh_bound * (
            next_odd**-5 + 1.0 / (8.0 * next_odd**4)
        )
        scale = max(abs(total), np.finfo(float).tiny)
        relative_tail = float(tail_bound / scale)
        if relative_tail <= relative_tolerance:
            return total, q + 1, relative_tail
    raise RuntimeError(
        "torsional series did not reach the requested tolerance "
        f"within {max_terms} terms (estimated relative tail={relative_tail:.3e})"
    )


def generalized_torsional_stiffness(
    properties: RotatedProperties,
    geometry: Geometry | None = None,
    *,
    relative_tolerance: float = 1.0e-12,
    max_terms: int = 100_000,
) -> TorsionalStiffness:
    """Evaluate equations (2.8)--(2.10), retaining the printed series form."""

    geom = geometry or Geometry()
    sqrt_gxz_over_gbar = np.sqrt(properties.Gxz / properties.Gxy_bar)
    base_argument = (
        np.pi * geom.b / (2.0 * geom.a) * sqrt_gxz_over_gbar
    )
    series_sum, terms_used, relative_tail = _odd_fifth_series(
        base_argument,
        relative_tolerance=relative_tolerance,
        max_terms=max_terms,
    )
    correction = (
        192.0
        / np.pi**5
        * geom.a
        / geom.b
        * np.sqrt(properties.Gxy_bar / properties.Gxz)
        * series_sum
    )
    Cbar = geom.a**3 * geom.b * properties.Gxy_bar / 3.0 * (1.0 - correction)
    C_T = Cbar / (
        1.0
        + Cbar
        * properties.Sbar16**2
        / (4.0 * properties.Sbar11 * geom.I_y)
    )
    return TorsionalStiffness(
        Cbar=Cbar,
        C_T=C_T,
        series_sum=series_sum,
        terms_used=terms_used,
        estimated_relative_tail=relative_tail,
    )


def make_rod_point(
    theta_deg: float,
    *,
    material_mode: MaterialMode = "elastic",
    loss_scale: float = 1.0,
    geometry: Geometry | None = None,
    material: BookMaterial | None = None,
    series_relative_tolerance: float = 1.0e-12,
    series_max_terms: int = 100_000,
) -> RodPoint:
    geom = geometry or Geometry()
    mat = material or BookMaterial()
    props = rotate_material(
        theta_deg,
        material=mat,
        material_mode=material_mode,
        loss_scale=loss_scale,
    )
    torsion = generalized_torsional_stiffness(
        props,
        geom,
        relative_tolerance=series_relative_tolerance,
        max_terms=series_max_terms,
    )
    return RodPoint(
        geometry=geom,
        material=mat,
        material_mode=material_mode,
        loss_scale=float(loss_scale),
        properties=props,
        torsion=torsion,
    )


def state_matrix(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
    """Physical state matrix for y=[w, psi_b, Phi_t, Q, M, M_T]."""

    g = point.geometry
    p = point.properties
    rho = point.material.rho
    omega2 = complex(omega) ** 2
    matrix = np.zeros((6, 6), dtype=np.complex128)
    matrix[0, 1] = 1.0
    matrix[0, 3] = 1.0 / (g.shear_factor * p.Gxz * g.area)
    matrix[1, 4] = 1.0 / (p.Ex * g.I_y)
    matrix[1, 5] = -p.Sbar16 / (2.0 * g.I_y)
    matrix[2, 4] = -p.Sbar16 / (2.0 * g.I_y)
    matrix[2, 5] = 1.0 / point.torsion.C_T
    matrix[3, 0] = -rho * g.area * omega2
    matrix[4, 1] = -rho * g.I_y * omega2
    matrix[4, 3] = -1.0
    matrix[5, 2] = -rho * g.I_p * omega2
    return matrix


def _state_scales(point: RodPoint) -> NDArray[np.float64]:
    g = point.geometry
    bending = abs(point.properties.Ex) * g.I_y
    torsion = abs(point.torsion.C_T)
    return np.array(
        [
            g.length,
            1.0,
            1.0,
            bending / g.length**2,
            bending / g.length,
            torsion / g.length,
        ],
        dtype=float,
    )


def scaled_state_matrix(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
    scales = _state_scales(point)
    physical = state_matrix(omega, point)
    return point.geometry.length * (physical * scales[np.newaxis, :]) / scales[:, np.newaxis]


def physical_state_transfer_matrix(
    omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Return the full-length corrected transfer matrix in physical states.

    This is the narrow public bridge needed by the coupled-rod diagnostic.  It
    keeps the private state scales centralized here and does not redefine the
    accepted ``state_corrected`` equations.
    """

    scales = _state_scales(point)
    scaled_transfer = expm(scaled_state_matrix(omega, point))
    return (
        scales[:, np.newaxis]
        * scaled_transfer
        / scales[np.newaxis, :]
    )


def state_boundary_matrix(omega: complex, point: RodPoint) -> NDArray[np.complex128]:
    transfer = expm(scaled_state_matrix(omega, point))
    return transfer[3:6, 0:3]


def cantilever_clamp_matrix(
    point: RodPoint,
    clamp_variant: ClampVariant,
    *,
    scaled: bool = False,
) -> NDArray[np.complex128]:
    """Return the left-end reaction-to-state matrix for a cantilever.

    The physical matrix maps ``[Q(0), M(0), M_T(0)]`` to
    ``[w, psi_b, Phi_t, Q, M, M_T]``.  ``scaled=True`` returns the
    equivalent well-conditioned matrix used with :func:`scaled_state_matrix`.
    """

    matrix = np.zeros((6, 3), dtype=np.complex128)
    matrix[3:, :] = np.eye(3, dtype=np.complex128)
    if clamp_variant == "book_slope_clamp":
        shear_stiffness = (
            point.geometry.shear_factor
            * point.properties.Gxz
            * point.geometry.area
        )
        matrix[1, 0] = -1.0 / shear_stiffness
    elif clamp_variant != "timoshenko_section_clamp":
        raise ValueError(f"unsupported clamp variant: {clamp_variant}")
    if not scaled:
        return matrix
    state_scales = _state_scales(point)
    return (
        matrix
        * state_scales[3:][np.newaxis, :]
        / state_scales[:, np.newaxis]
    )


def state_cantilever_boundary_matrix(
    omega: complex, point: RodPoint, clamp_variant: ClampVariant
) -> NDArray[np.complex128]:
    """Return ``S_f exp(A L) B_clamp`` in scaled state coordinates."""

    transfer = expm(scaled_state_matrix(omega, point))
    return transfer[3:6, :] @ cantilever_clamp_matrix(
        point, clamp_variant, scaled=True
    )


def partial_bending_scaled_system(
    omega: complex, point: RodPoint
) -> tuple[NDArray[np.complex128], NDArray[np.float64]]:
    """Narrow equation-(2.14) state system using free ``Ex`` and ``Gxz``."""

    g = point.geometry
    p = point.properties
    rho = point.material.rho
    omega2 = complex(omega) ** 2
    physical = np.zeros((4, 4), dtype=np.complex128)
    physical[0, 1] = 1.0
    physical[0, 2] = 1.0 / (g.shear_factor * p.Gxz * g.area)
    physical[1, 3] = 1.0 / (p.Ex * g.I_y)
    physical[2, 0] = -rho * g.area * omega2
    physical[3, 1] = -rho * g.I_y * omega2
    physical[3, 2] = -1.0
    bending = abs(p.Ex) * g.I_y
    scales = np.array(
        [g.length, 1.0, bending / g.length**2, bending / g.length],
        dtype=float,
    )
    scaled = g.length * (physical * scales[np.newaxis, :]) / scales[:, np.newaxis]
    return scaled, scales


def partial_bending_boundary_matrix(
    omega: complex,
    point: RodPoint,
    formulation: Formulation = "state_corrected",
) -> NDArray[np.complex128]:
    """Equation-(2.14) characteristic matrix with the printed slope clamp."""

    if formulation != "state_corrected":
        raise ValueError("partial bending is implemented only in the state form")
    system, scales = partial_bending_scaled_system(omega, point)
    physical_clamp = np.zeros((4, 2), dtype=np.complex128)
    physical_clamp[2:, :] = np.eye(2, dtype=np.complex128)
    shear_stiffness = (
        point.geometry.shear_factor
        * point.properties.Gxz
        * point.geometry.area
    )
    physical_clamp[1, 0] = -1.0 / shear_stiffness
    scaled_clamp = (
        physical_clamp * scales[2:][np.newaxis, :] / scales[:, np.newaxis]
    )
    return expm(system)[2:4, :] @ scaled_clamp


def partial_bending_mode_shape(
    omega: complex,
    point: RodPoint,
    x_over_length: NDArray[np.float64] | Sequence[float],
) -> NDArray[np.complex128]:
    """Normalized ``[w/L, psi_b]`` partial-bending fields."""

    grid = np.asarray(x_over_length, dtype=float)
    if np.any(grid < 0.0) or np.any(grid > 1.0):
        raise ValueError("x_over_length must lie in [0, 1]")
    system, scales = partial_bending_scaled_system(omega, point)
    physical_clamp = np.zeros((4, 2), dtype=np.complex128)
    physical_clamp[2:, :] = np.eye(2, dtype=np.complex128)
    shear_stiffness = (
        point.geometry.shear_factor
        * point.properties.Gxz
        * point.geometry.area
    )
    physical_clamp[1, 0] = -1.0 / shear_stiffness
    scaled_clamp = (
        physical_clamp * scales[2:][np.newaxis, :] / scales[:, np.newaxis]
    )
    initial = scaled_clamp @ _null_vector(
        partial_bending_boundary_matrix(omega, point)
    )
    fields = np.vstack([(expm(system * x) @ initial)[:2] for x in grid])
    norm = np.linalg.norm(fields)
    return fields / norm if norm else fields


def eliminated_coefficients(
    omega: complex, point: RodPoint, *, printed: bool
) -> dict[str, complex]:
    g = point.geometry
    p = point.properties
    rho = point.material.rho
    omega2 = complex(omega) ** 2
    common = rho * omega2
    sign = -1.0 if printed else 1.0
    return {
        "a0": common * (1.0 / p.Ex + 1.0 / (g.shear_factor * p.Gxz)),
        "b0": rho**2
        * omega2**2
        / (p.Ex * g.shear_factor * p.Gxz)
        - common * g.area / (g.I_y * p.Ex),
        "c0": -common * p.mu_x_xy * g.I_p / (2.0 * p.Gxy * g.I_y),
        "d0": sign * common * g.I_p / point.torsion.Cbar,
        "e0": p.mu_xy_x / 2.0,
        "f0": sign * common * p.mu_xy_x / (2.0 * g.shear_factor * p.Gxz),
        "ag": common / (g.shear_factor * p.Gxz),
        "bg": common * (1.0 / p.Ex + 1.0 / (g.shear_factor * p.Gxz)),
        "cg": -common * p.mu_x_xy * g.I_p / (2.0 * p.Gxy * g.I_y),
    }


def scaled_eliminated_system(
    omega: complex, point: RodPoint, *, printed: bool
) -> tuple[NDArray[np.complex128], NDArray[np.complex128], NDArray[np.complex128]]:
    """Dimensionless (2.16) system and free-end matrices from (2.17)."""

    c = eliminated_coefficients(omega, point, printed=printed)
    length = point.geometry.length
    physical = np.zeros((6, 6), dtype=np.complex128)
    physical[0, 1] = 1.0
    physical[1, 2] = 1.0
    physical[2, 3] = 1.0
    physical[3, 0] = -c["b0"]
    physical[3, 2] = -c["a0"]
    physical[3, 5] = -c["c0"]
    physical[4, 5] = 1.0
    physical[5, 1] = -c["f0"]
    physical[5, 3] = -c["e0"]
    physical[5, 4] = -c["d0"]

    scales = np.array(
        [length, 1.0, 1.0 / length, 1.0 / length**2, 1.0, 1.0 / length],
        dtype=float,
    )
    scaled = length * (physical * scales[np.newaxis, :]) / scales[:, np.newaxis]

    boundary = np.zeros((3, 6), dtype=np.complex128)
    boundary[0, 0] = c["ag"] * length**2
    boundary[0, 2] = 1.0
    boundary[1, 1] = c["bg"] * length**2
    boundary[1, 3] = 1.0
    boundary[1, 4] = c["cg"] * length**2
    boundary[2, 5] = 1.0

    unknowns = np.zeros((6, 3), dtype=np.complex128)
    unknowns[0, 0] = 1.0
    unknowns[2, 0] = -c["ag"] * length**2
    unknowns[1, 1] = 1.0
    unknowns[3, 1] = -c["bg"] * length**2
    unknowns[4, 2] = 1.0
    unknowns[3, 2] = -c["cg"] * length**2
    return scaled, boundary, unknowns


def eliminated_boundary_matrix(
    omega: complex, point: RodPoint, *, printed: bool
) -> NDArray[np.complex128]:
    matrix, boundary, unknowns = scaled_eliminated_system(
        omega, point, printed=printed
    )
    return boundary @ expm(matrix) @ unknowns


def eliminated_cantilever_boundary_matrix(
    omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Corrected (2.16)--(2.18) matrix for the printed slope clamp."""

    matrix, free_boundary, _ = scaled_eliminated_system(
        omega, point, printed=False
    )
    clamp_unknowns = np.zeros((6, 3), dtype=np.complex128)
    clamp_unknowns[2, 0] = 1.0
    clamp_unknowns[3, 1] = 1.0
    clamp_unknowns[5, 2] = 1.0
    return free_boundary @ expm(matrix) @ clamp_unknowns


def cantilever_boundary_matrix(
    omega: complex,
    point: RodPoint,
    clamp_variant: ClampVariant,
    formulation: Formulation = "state_corrected",
) -> NDArray[np.complex128]:
    """Cantilever characteristic matrix for one supported formulation."""

    if formulation == "state_corrected":
        return state_cantilever_boundary_matrix(omega, point, clamp_variant)
    if formulation == "eliminated_corrected" and clamp_variant == "book_slope_clamp":
        return eliminated_cantilever_boundary_matrix(omega, point)
    if formulation == "eliminated_printed":
        raise ValueError("the printed-sign variant is not used for the cantilever gate")
    raise ValueError(
        f"formulation {formulation!r} is not implemented for {clamp_variant!r}"
    )


def boundary_matrix(
    omega: complex, point: RodPoint, formulation: Formulation
) -> NDArray[np.complex128]:
    if formulation == "state_corrected":
        return state_boundary_matrix(omega, point)
    if formulation == "eliminated_corrected":
        return eliminated_boundary_matrix(omega, point, printed=False)
    if formulation == "eliminated_printed":
        return eliminated_boundary_matrix(omega, point, printed=True)
    raise ValueError(f"unsupported formulation: {formulation}")


def boundary_quality(
    omega: complex,
    point: RodPoint,
    formulation: Formulation,
    *,
    boundary_matrix_builder: BoundaryMatrixBuilder | None = None,
) -> BoundaryQuality:
    builder = boundary_matrix if boundary_matrix_builder is None else boundary_matrix_builder
    matrix = builder(omega, point, formulation)
    singular = np.linalg.svd(matrix, compute_uv=False)
    sigma_max = float(singular[0])
    sigma_min = float(singular[-1])
    relative = sigma_min / sigma_max if sigma_max > 0.0 else 0.0
    return BoundaryQuality(
        determinant=complex(np.linalg.det(matrix)),
        sigma_min=sigma_min,
        sigma_max=sigma_max,
        relative_singular_residual=float(relative),
    )


def _real_determinant(
    omega: float,
    point: RodPoint,
    formulation: Formulation,
    boundary_matrix_builder: BoundaryMatrixBuilder | None = None,
) -> float:
    value = boundary_quality(
        float(omega),
        point,
        formulation,
        boundary_matrix_builder=boundary_matrix_builder,
    ).determinant
    tolerance = 2.0e-9 * max(abs(value.real), 1.0e-30)
    if abs(value.imag) > tolerance:
        raise RuntimeError(
            "elastic determinant acquired a non-negligible imaginary part: "
            f"det={value!r}, omega={omega:.16g}"
        )
    return float(value.real)


def _root_result(
    omega: complex,
    point: RodPoint,
    formulation: Formulation,
    *,
    refinements: int,
    status: str,
    boundary_matrix_builder: BoundaryMatrixBuilder | None = None,
) -> RootResult:
    quality = boundary_quality(
        omega,
        point,
        formulation,
        boundary_matrix_builder=boundary_matrix_builder,
    )
    builder = boundary_matrix if boundary_matrix_builder is None else boundary_matrix_builder
    matrix_size = builder(omega, point, formulation).shape[0]
    determinant_scale = max(quality.sigma_max**matrix_size, np.finfo(float).tiny)
    return RootResult(
        omega=complex(omega),
        frequency_hz=float(np.real(omega) / (2.0 * np.pi)),
        determinant_residual=float(abs(quality.determinant) / determinant_scale),
        raw_determinant_abs=float(abs(quality.determinant)),
        sigma_min=quality.sigma_min,
        sigma_max=quality.sigma_max,
        relative_singular_residual=quality.relative_singular_residual,
        refinements=int(refinements),
        status=status,
    )


def _append_unique_root(roots: list[RootResult], candidate: RootResult) -> None:
    for existing in roots:
        scale = max(abs(existing.omega), abs(candidate.omega), 1.0)
        if abs(existing.omega - candidate.omega) <= 2.0e-8 * scale:
            return
    roots.append(candidate)
    roots.sort(key=lambda item: item.frequency_hz)


def find_elastic_roots(
    point: RodPoint,
    formulation: Formulation,
    *,
    num_roots: int,
    scan_step_hz: float = 10.0,
    initial_max_hz: float = 1000.0,
    max_hz: float = 100_000.0,
    sigma_dip_trigger: float = 1.0e-5,
    boundary_matrix_builder: BoundaryMatrixBuilder | None = None,
) -> list[RootResult]:
    """Find positive roots by determinant sign changes with an SVD dip guard."""

    if point.material_mode != "elastic":
        raise ValueError("find_elastic_roots requires an elastic RodPoint")
    if num_roots < 1:
        raise ValueError("num_roots must be positive")
    if scan_step_hz <= 0.0 or initial_max_hz <= 0.0:
        raise ValueError("scan ranges must be positive")

    roots: list[RootResult] = []
    first_hz = max(1.0e-3, 0.05 * scan_step_hz)
    last_hz = first_hz
    last_omega = 2.0 * np.pi * last_hz
    last_quality = boundary_quality(
        last_omega,
        point,
        formulation,
        boundary_matrix_builder=boundary_matrix_builder,
    )
    last_det = _real_determinant(
        last_omega, point, formulation, boundary_matrix_builder
    )
    samples: list[tuple[float, float, float]] = [
        (last_omega, last_det, last_quality.relative_singular_residual)
    ]
    current_max_hz = max(initial_max_hz, first_hz + scan_step_hz)

    while len(roots) < num_roots:
        if current_max_hz > max_hz * (1.0 + 1.0e-12):
            raise RuntimeError(
                f"only {len(roots)} positive roots found below {max_hz:g} Hz"
            )
        grid_hz = np.arange(last_hz + scan_step_hz, current_max_hz, scan_step_hz)
        grid_hz = np.append(grid_hz, current_max_hz)
        for frequency_hz in grid_hz:
            omega = float(2.0 * np.pi * frequency_hz)
            quality = boundary_quality(
                omega,
                point,
                formulation,
                boundary_matrix_builder=boundary_matrix_builder,
            )
            determinant = _real_determinant(
                omega, point, formulation, boundary_matrix_builder
            )
            samples.append((omega, determinant, quality.relative_singular_residual))
            if last_det == 0.0 or determinant == 0.0 or last_det * determinant < 0.0:
                left, right = last_omega, omega
                if last_det == 0.0:
                    left = max(first_hz * 2.0 * np.pi, last_omega - 0.25 * (omega - last_omega))
                if determinant == 0.0:
                    right = omega + 0.25 * (omega - last_omega)
                try:
                    solved, info = brentq(
                        lambda value: _real_determinant(
                            value, point, formulation, boundary_matrix_builder
                        ),
                        left,
                        right,
                        xtol=1.0e-8,
                        rtol=2.0e-13,
                        maxiter=150,
                        full_output=True,
                        disp=False,
                    )
                except ValueError:
                    solved = float(omega if determinant == 0.0 else last_omega)
                    info = None
                candidate = _root_result(
                    solved,
                    point,
                    formulation,
                    refinements=0 if info is None else info.iterations,
                    status="accepted_brent",
                    boundary_matrix_builder=boundary_matrix_builder,
                )
                _append_unique_root(roots, candidate)
            last_hz = float(frequency_hz)
            last_omega = omega
            last_det = determinant

        # A sharp singular-value dip without a sampled sign change triggers a
        # local determinant refinement.  Sigma_min is never accepted alone.
        for index in range(1, len(samples) - 1):
            omega_mid, _, sigma_mid = samples[index]
            left_sample = samples[index - 1]
            right_sample = samples[index + 1]
            if not (
                sigma_mid < sigma_dip_trigger
                and sigma_mid < 0.05 * min(left_sample[2], right_sample[2])
                and left_sample[1] * right_sample[1] > 0.0
            ):
                continue
            optimized = minimize_scalar(
                lambda value: boundary_quality(
                    value,
                    point,
                    formulation,
                    boundary_matrix_builder=boundary_matrix_builder,
                ).relative_singular_residual,
                bounds=(left_sample[0], right_sample[0]),
                method="bounded",
                options={"xatol": 1.0e-8, "maxiter": 120},
            )
            if not optimized.success:
                continue
            candidate = _root_result(
                float(optimized.x),
                point,
                formulation,
                refinements=int(optimized.nfev),
                status="sigma_dip_det_verified",
                boundary_matrix_builder=boundary_matrix_builder,
            )
            if (
                candidate.relative_singular_residual <= 1.0e-9
                and candidate.determinant_residual <= 1.0e-8
            ):
                _append_unique_root(roots, candidate)

        if len(roots) >= num_roots:
            break
        current_max_hz = min(max_hz, 2.0 * current_max_hz)
        if last_hz >= max_hz:
            current_max_hz = max_hz * (1.0 + 1.0e-9)

    roots = sorted(roots, key=lambda item: item.frequency_hz)[:num_roots]
    for index, item in enumerate(roots):
        distances: list[float] = []
        if index:
            distances.append(item.frequency_hz - roots[index - 1].frequency_hz)
        if index + 1 < len(roots):
            distances.append(roots[index + 1].frequency_hz - item.frequency_hz)
        distance = min(distances) if distances else float("inf")
        roots[index] = replace(item, min_neighbor_distance_hz=float(distance))
    return roots


def solve_complex_root(
    point: RodPoint,
    formulation: Formulation,
    initial_omega: complex,
    *,
    tolerance: float = 1.0e-10,
    max_evaluations: int = 400,
    boundary_matrix_builder: BoundaryMatrixBuilder | None = None,
) -> RootResult:
    if point.material_mode != "book_complex":
        raise ValueError("solve_complex_root requires a book_complex RodPoint")

    initial = complex(initial_omega)

    def residual(values: Sequence[float]) -> NDArray[np.float64]:
        omega = complex(float(values[0]), float(values[1]))
        determinant = boundary_quality(
            omega,
            point,
            formulation,
            boundary_matrix_builder=boundary_matrix_builder,
        ).determinant
        return np.array([determinant.real, determinant.imag], dtype=float)

    solved = root(
        residual,
        np.array([initial.real, max(initial.imag, 0.0)], dtype=float),
        method="hybr",
        options={"xtol": tolerance, "maxfev": max_evaluations},
    )
    omega = complex(float(solved.x[0]), float(solved.x[1]))
    candidate = _root_result(
        omega,
        point,
        formulation,
        refinements=int(solved.nfev),
        status=(
            "accepted_complex"
            if solved.success
            else "accepted_complex_quality_after_solver_warning"
        ),
        boundary_matrix_builder=boundary_matrix_builder,
    )
    physical = candidate.omega.real > 0.0 and candidate.omega.imag >= -1.0e-9
    quality_ok = candidate.relative_singular_residual <= 2.0e-7
    if not physical or not quality_ok:
        return replace(candidate, status="rejected_complex_quality")
    if candidate.omega.imag < 0.0:
        return replace(
            candidate,
            omega=complex(candidate.omega.real, 0.0),
            status="accepted_complex_imag_clipped",
        )
    return candidate


def continue_loss_root(
    point_factory: Callable[[float], RodPoint],
    formulation: Formulation,
    elastic_omega: float,
    *,
    loss_steps: Sequence[float] = (0.2, 0.4, 0.6, 0.8, 1.0),
    initial_predictor: complex | None = None,
    boundary_matrix_builder: BoundaryMatrixBuilder | None = None,
) -> RootResult:
    predictor = complex(elastic_omega if initial_predictor is None else initial_predictor)
    result: RootResult | None = None
    for scale in loss_steps:
        result = solve_complex_root(
            point_factory(float(scale)),
            formulation,
            predictor,
            boundary_matrix_builder=boundary_matrix_builder,
        )
        if result.status == "rejected_complex_quality":
            raise RuntimeError(
                f"complex continuation failed at loss_scale={scale:g}, "
                f"omega0={predictor!r}"
            )
        predictor = result.omega
    assert result is not None
    return result


def modal_loss_factors(omega: complex) -> tuple[float, float, float]:
    real = float(np.real(omega))
    imag = float(np.imag(omega))
    if real <= 0.0:
        return float("nan"), float("nan"), float("nan")
    ratio = imag / real
    exact = 2.0 * ratio / (1.0 - ratio**2)
    approximate = 2.0 * ratio
    relative_difference = (
        abs(exact - approximate) / abs(exact) if exact != 0.0 else 0.0
    )
    return float(exact), float(approximate), float(relative_difference)


def rigid_body_nullity(point: RodPoint, formulation: Formulation = "state_corrected") -> int:
    matrix = boundary_matrix(0.0, point, formulation)
    return int(3 - np.linalg.matrix_rank(matrix, tol=1.0e-12))


def decoupled_boundary_factors(
    omega: complex, point: RodPoint
) -> tuple[complex, complex, complex]:
    """Return full, bending, and torsion factors when Sbar16 is zero."""

    matrix = state_boundary_matrix(omega, point)
    bending = complex(np.linalg.det(matrix[np.ix_([0, 1], [0, 1])]))
    torsion = complex(matrix[2, 2])
    return complex(np.linalg.det(matrix)), bending, torsion


def decoupled_cantilever_boundary_factors(
    omega: complex,
    point: RodPoint,
    clamp_variant: ClampVariant = "book_slope_clamp",
) -> tuple[complex, complex, complex]:
    """Full, bending, and torsion factors when cantilever coupling vanishes."""

    matrix = state_cantilever_boundary_matrix(omega, point, clamp_variant)
    bending = complex(np.linalg.det(matrix[np.ix_([0, 1], [0, 1])]))
    torsion = complex(matrix[2, 2])
    return complex(np.linalg.det(matrix)), bending, torsion


def cantilever_zero_frequency_nullity(
    point: RodPoint, clamp_variant: ClampVariant
) -> int:
    """Nullity of the 3x3 cantilever boundary matrix at zero frequency."""

    matrix = state_cantilever_boundary_matrix(0.0, point, clamp_variant)
    return int(3 - np.linalg.matrix_rank(matrix, tol=1.0e-12))


def free_free_torsion_omega(point: RodPoint, mode_number: int) -> complex:
    if mode_number < 1:
        raise ValueError("torsion mode_number starts at 1")
    g = point.geometry
    return (
        mode_number
        * np.pi
        / g.length
        * np.sqrt(point.torsion.C_T / (point.material.rho * g.I_p))
    )


def fixed_free_torsion_omega(
    point: RodPoint, mode_number: int, *, partial: bool = False
) -> complex:
    """Analytical fixed-free torsion root; partial uses equation (2.15)."""

    if mode_number < 1:
        raise ValueError("torsion mode_number starts at 1")
    stiffness = point.torsion.Cbar if partial else point.torsion.C_T
    g = point.geometry
    return (
        (2 * mode_number - 1)
        * np.pi
        / (2.0 * g.length)
        * np.sqrt(stiffness / (point.material.rho * g.I_p))
    )


def partial_torsion_mode_shape(
    mode_number: int,
    x_over_length: NDArray[np.float64] | Sequence[float],
) -> NDArray[np.float64]:
    """Normalized fixed-free equation-(2.15) twist shape."""

    if mode_number < 1:
        raise ValueError("torsion mode_number starts at 1")
    grid = np.asarray(x_over_length, dtype=float)
    if np.any(grid < 0.0) or np.any(grid > 1.0):
        raise ValueError("x_over_length must lie in [0, 1]")
    shape = np.sin((2 * mode_number - 1) * np.pi * grid / 2.0)
    norm = np.linalg.norm(shape)
    return shape / norm if norm else shape


def with_gxz_scale(point: RodPoint, scale: float) -> RodPoint:
    """Return the artificial shear-rigid diagnostic point ``Gxz *= scale``."""

    if scale <= 0.0:
        raise ValueError("Gxz scale must be positive")
    properties = replace(point.properties, Gxz=point.properties.Gxz * scale)
    torsion = generalized_torsional_stiffness(properties, point.geometry)
    return replace(point, properties=properties, torsion=torsion)


def _null_vector(matrix: NDArray[np.complex128]) -> NDArray[np.complex128]:
    _, _, vh = np.linalg.svd(matrix)
    vector = vh.conj().T[:, -1]
    pivot = int(np.argmax(np.abs(vector)))
    if abs(vector[pivot]) > 0.0:
        vector *= np.exp(-1j * np.angle(vector[pivot]))
    norm = np.linalg.norm(vector)
    return vector / norm if norm else vector


def mode_shape(
    omega: complex,
    point: RodPoint,
    formulation: Formulation,
    x_over_length: NDArray[np.float64] | Sequence[float],
) -> NDArray[np.complex128]:
    grid = np.asarray(x_over_length, dtype=float)
    if np.any(grid < 0.0) or np.any(grid > 1.0):
        raise ValueError("x_over_length must lie in [0, 1]")
    if formulation == "state_corrected":
        system = scaled_state_matrix(omega, point)
        initial = np.zeros(6, dtype=np.complex128)
        initial[:3] = _null_vector(state_boundary_matrix(omega, point))
        fields = np.vstack([(expm(system * x) @ initial)[:3] for x in grid])
    else:
        printed = formulation == "eliminated_printed"
        system, _, unknowns = scaled_eliminated_system(
            omega, point, printed=printed
        )
        initial = unknowns @ _null_vector(
            eliminated_boundary_matrix(omega, point, printed=printed)
        )
        raw = np.vstack([expm(system * x) @ initial for x in grid])
        fields = raw[:, [0, 1, 4]]
    norm = np.linalg.norm(fields)
    return fields / norm if norm else fields


def cantilever_state_trajectory(
    omega: complex,
    point: RodPoint,
    clamp_variant: ClampVariant,
    x_over_length: NDArray[np.float64] | Sequence[float],
) -> NDArray[np.complex128]:
    """Reconstruct the raw scaled corrected state of a cantilever mode."""

    grid = np.asarray(x_over_length, dtype=float)
    if np.any(grid < 0.0) or np.any(grid > 1.0):
        raise ValueError("x_over_length must lie in [0, 1]")
    characteristic = state_cantilever_boundary_matrix(
        omega, point, clamp_variant
    )
    initial = cantilever_clamp_matrix(
        point, clamp_variant, scaled=True
    ) @ _null_vector(characteristic)
    system = scaled_state_matrix(omega, point)
    return np.vstack([expm(system * x) @ initial for x in grid])


def cantilever_mode_shape(
    omega: complex,
    point: RodPoint,
    clamp_variant: ClampVariant,
    x_over_length: NDArray[np.float64] | Sequence[float],
) -> NDArray[np.complex128]:
    """Return normalized ``[w/L, psi_b, Phi_t]`` cantilever fields."""

    fields = cantilever_state_trajectory(
        omega, point, clamp_variant, x_over_length
    )[:, :3]
    norm = np.linalg.norm(fields)
    return fields / norm if norm else fields


def cantilever_boundary_residuals(
    omega: complex, point: RodPoint, clamp_variant: ClampVariant
) -> dict[str, complex]:
    """Return the six physical boundary quantities and clamp diagnostics."""

    scaled = cantilever_state_trajectory(
        omega, point, clamp_variant, np.array([0.0, 1.0])
    )
    physical = scaled * _state_scales(point)[np.newaxis, :]
    left, right = physical
    shear_stiffness = (
        point.geometry.shear_factor
        * point.properties.Gxz
        * point.geometry.area
    )
    slope_left = left[1] + left[3] / shear_stiffness
    return {
        "w_0": left[0],
        "psi_b_0": left[1],
        "Phi_t_0": left[2],
        "Q_0": left[3],
        "Q_0_over_Ks": left[3] / shear_stiffness,
        "w_prime_0": slope_left,
        "slope_compensation": left[1] + left[3] / shear_stiffness,
        "Q_L": right[3],
        "M_L": right[4],
        "M_T_L": right[5],
    }


def cantilever_energy_fractions(
    omega: complex,
    point: RodPoint,
    clamp_variant: ClampVariant,
    *,
    samples: int = 201,
) -> tuple[float, float, float]:
    """Diagnostic bending/shear/torsion strain-energy fractions."""

    grid = np.linspace(0.0, 1.0, samples)
    scaled = cantilever_state_trajectory(omega, point, clamp_variant, grid)
    physical = scaled * _state_scales(point)[np.newaxis, :]
    psi = physical[:, 1]
    shear = physical[:, 3]
    moment = physical[:, 4]
    torque = physical[:, 5]
    g = point.geometry
    p = point.properties
    curvature = moment / (p.Ex * g.I_y) - p.Sbar16 * torque / (2.0 * g.I_y)
    twist = -p.Sbar16 * moment / (2.0 * g.I_y) + torque / point.torsion.C_T
    shear_angle = shear / (g.shear_factor * p.Gxz * g.area)
    bending_density = abs(p.Ex) * g.I_y * np.abs(curvature) ** 2
    shear_density = (
        g.shear_factor * abs(p.Gxz) * g.area * np.abs(shear_angle) ** 2
    )
    torsion_density = abs(point.torsion.C_T) * np.abs(twist) ** 2
    bending = float(np.trapezoid(bending_density, grid))
    shear_energy = float(np.trapezoid(shear_density, grid))
    torsion = float(np.trapezoid(torsion_density, grid))
    total = bending + shear_energy + torsion
    if total <= 0.0:
        return 0.0, 0.0, 0.0
    return bending / total, shear_energy / total, torsion / total


def modal_assurance(
    left: NDArray[np.complex128], right: NDArray[np.complex128]
) -> float:
    a = np.asarray(left, dtype=np.complex128).ravel()
    b = np.asarray(right, dtype=np.complex128).ravel()
    denominator = float(np.vdot(a, a).real * np.vdot(b, b).real)
    if denominator <= 0.0:
        return 0.0
    return float(abs(np.vdot(a, b)) ** 2 / denominator)


def assign_modes_by_mac(
    previous_shapes: Sequence[NDArray[np.complex128]],
    current_shapes: Sequence[NDArray[np.complex128]],
) -> tuple[list[int], NDArray[np.float64]]:
    """Assign current sorted shapes to previous continuity labels by MAC."""

    if len(previous_shapes) != len(current_shapes):
        raise ValueError("previous and current shape counts must agree")
    mac_matrix = np.array(
        [
            [modal_assurance(previous, current) for current in current_shapes]
            for previous in previous_shapes
        ],
        dtype=float,
    )
    tracked_indices, sorted_indices = linear_sum_assignment(-mac_matrix)
    assignment = [-1] * len(previous_shapes)
    for tracked_index, sorted_index in zip(tracked_indices, sorted_indices):
        assignment[int(tracked_index)] = int(sorted_index)
    return assignment, mac_matrix


def state_to_eliminated_state(
    scaled_state: NDArray[np.complex128], omega: complex, point: RodPoint
) -> NDArray[np.complex128]:
    """Map a scaled corrected state to scaled [w,w',w'',w''',Phi,Phi']."""

    g = point.geometry
    p = point.properties
    rho = point.material.rho
    state_scales = _state_scales(point)
    y = state_scales * np.asarray(scaled_state, dtype=np.complex128)
    w, psi, phi, shear, moment, torque = y
    omega2 = complex(omega) ** 2
    w1 = psi + shear / (g.shear_factor * p.Gxz * g.area)
    w2 = (
        moment / (p.Ex * g.I_y)
        - p.Sbar16 * torque / (2.0 * g.I_y)
        - rho * omega2 * w / (g.shear_factor * p.Gxz)
    )
    w3 = (
        (-shear - rho * g.I_y * omega2 * psi) / (p.Ex * g.I_y)
        + p.Sbar16 * rho * g.I_p * omega2 * phi / (2.0 * g.I_y)
        - rho * omega2 * w1 / (g.shear_factor * p.Gxz)
    )
    phi1 = -p.Sbar16 * moment / (2.0 * g.I_y) + torque / point.torsion.C_T
    physical = np.array([w, w1, w2, w3, phi, phi1], dtype=np.complex128)
    eliminated_scales = np.array(
        [g.length, 1.0, 1.0 / g.length, 1.0 / g.length**2, 1.0, 1.0 / g.length],
        dtype=float,
    )
    return physical / eliminated_scales


def formulation_mapping_residual(
    omega: complex, point: RodPoint, *, samples: int = 11
) -> float:
    """Compare corrected state propagation with corrected (2.16) propagation."""

    state_system = scaled_state_matrix(omega, point)
    state_initial = np.zeros(6, dtype=np.complex128)
    state_initial[:3] = _null_vector(state_boundary_matrix(omega, point))
    eliminated_system, _, _ = scaled_eliminated_system(
        omega, point, printed=False
    )
    eliminated_initial = state_to_eliminated_state(state_initial, omega, point)
    residual = 0.0
    for x in np.linspace(0.0, 1.0, samples):
        mapped = state_to_eliminated_state(
            expm(state_system * x) @ state_initial, omega, point
        )
        propagated = expm(eliminated_system * x) @ eliminated_initial
        scale = max(np.linalg.norm(mapped), np.linalg.norm(propagated), 1.0)
        residual = max(residual, float(np.linalg.norm(mapped - propagated) / scale))
    return residual


def cantilever_formulation_mapping_residual(
    omega: complex, point: RodPoint, *, samples: int = 11
) -> float:
    """Compare corrected state and eliminated fields for the slope clamp."""

    state_characteristic = state_cantilever_boundary_matrix(
        omega, point, "book_slope_clamp"
    )
    state_initial = cantilever_clamp_matrix(
        point, "book_slope_clamp", scaled=True
    ) @ _null_vector(state_characteristic)
    eliminated_system, _, _ = scaled_eliminated_system(
        omega, point, printed=False
    )
    eliminated_unknowns = np.zeros((6, 3), dtype=np.complex128)
    eliminated_unknowns[2, 0] = 1.0
    eliminated_unknowns[3, 1] = 1.0
    eliminated_unknowns[5, 2] = 1.0
    eliminated_initial = eliminated_unknowns @ _null_vector(
        eliminated_cantilever_boundary_matrix(omega, point)
    )
    mapped_initial = state_to_eliminated_state(state_initial, omega, point)
    denominator = np.vdot(eliminated_initial, eliminated_initial)
    factor = (
        np.vdot(eliminated_initial, mapped_initial) / denominator
        if abs(denominator) > 0.0
        else 1.0 + 0.0j
    )
    residual = 0.0
    for x in np.linspace(0.0, 1.0, samples):
        mapped = state_to_eliminated_state(
            expm(scaled_state_matrix(omega, point) * x) @ state_initial,
            omega,
            point,
        )
        propagated = factor * expm(eliminated_system * x) @ eliminated_initial
        scale = max(np.linalg.norm(mapped), np.linalg.norm(propagated), 1.0)
        residual = max(residual, float(np.linalg.norm(mapped - propagated) / scale))
    return residual


__all__ = [
    "BookMaterial",
    "BoundaryMatrixBuilder",
    "BoundaryQuality",
    "ClampVariant",
    "Formulation",
    "Geometry",
    "MaterialMode",
    "RodPoint",
    "RootResult",
    "RotatedProperties",
    "TorsionalStiffness",
    "boundary_matrix",
    "boundary_quality",
    "assign_modes_by_mac",
    "cantilever_boundary_matrix",
    "cantilever_boundary_residuals",
    "cantilever_clamp_matrix",
    "cantilever_energy_fractions",
    "cantilever_formulation_mapping_residual",
    "cantilever_geometry",
    "cantilever_mode_shape",
    "cantilever_state_trajectory",
    "cantilever_zero_frequency_nullity",
    "continue_loss_root",
    "decoupled_boundary_factors",
    "decoupled_cantilever_boundary_factors",
    "eliminated_boundary_matrix",
    "eliminated_coefficients",
    "find_elastic_roots",
    "fixed_free_torsion_omega",
    "formulation_mapping_residual",
    "free_free_torsion_omega",
    "generalized_torsional_stiffness",
    "hms_dx_209_material",
    "make_rod_point",
    "modal_assurance",
    "modal_loss_factors",
    "mode_shape",
    "partial_bending_boundary_matrix",
    "partial_bending_mode_shape",
    "partial_bending_scaled_system",
    "partial_torsion_mode_shape",
    "physical_state_transfer_matrix",
    "rigid_body_nullity",
    "rotate_material",
    "scaled_eliminated_system",
    "scaled_state_matrix",
    "solve_complex_root",
    "state_boundary_matrix",
    "state_cantilever_boundary_matrix",
    "state_matrix",
    "state_to_eliminated_state",
    "with_gxz_scale",
]
