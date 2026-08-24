"""Independent Rayleigh--Ritz audit of Reddy Table 4.3.3 discrepancies.

This diagnostic deliberately does not use a transfer matrix or a root finder
to construct its Rayleigh--Ritz spectrum.  It discretizes the two-field FSDT
energy with polynomial trial functions built from shifted Legendre
polynomials.  The existing RLB-0 library is used only to construct the reduced
beam properties and, *after* each Ritz solve, to obtain an independent
transfer-matrix comparison.

For ``W = w / length`` and ``xi = x / length`` the dimensionless functionals
are

``U_bar = integral(psi_,xi**2 + gamma*(W_,xi + psi)**2) dxi``

and

``T_bar = integral(W**2 + rho_r*psi**2) dxi``,

where ``gamma = S*length**2/D`` and ``rho_r = J/(m*length**2)``.  Their
generalized eigenvalue is ``lambda = omega**2*m*length**4/D``.  When ``J=0``
the rotational algebraic block is statically condensed exactly; it is never
regularized with an artificial rotary mass.

The script is a stable, diagnostic-only entry point because it introduces an
independent discretization and a new output contract.  It does not modify the
production RLB-0 validator or any historical result directory.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass, replace
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any, Iterable, Literal, Sequence

import numpy as np
from numpy.polynomial import Legendre
from numpy.typing import NDArray
from scipy.linalg import cholesky, eigh, expm, solve_triangular
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.lib import reddy_symmetric_laminated_beam as rlb  # noqa: E402


FloatArray = NDArray[np.float64]
BoundaryCondition = Literal["HH", "CC", "CF"]

SOURCE_CONTRACT_PATH = ROOT / "tests" / "data" / "reddy_ch4_table_4_3_3.json"
SOURCE_CONTRACT_SHA256 = (
    "3F479D04DBD601C4BBFEB9D32463956A5201D72ECC41A067B8234979D885C388"
)
SOURCE_PDF_PATH = (
    ROOT
    / "docs"
    / "literature"
    / "pdf"
    / "EB__Mechanics_of_Laminated_Composite_Plates_and_Shells_-JN_Reddy.pdf"
)
SOURCE_PDF_SHA256 = (
    "A8D4D0FA67C7073D6EB48903B868BFA157875DC60F54484376A5B510547B37EA"
)
LEGACY_TRANSFER_INVENTORY_PATH = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_single_beam"
    / "bending_source_comparison.csv"
)
LEGACY_TRANSFER_INVENTORY_SHA256 = (
    "73D4F0D966A6F502F48D3900375B6FC1F0BACE67AC636FC1AA7D9893063B7D92"
)
DEFAULT_OUTPUT_DIR = (
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_single_beam"
    / "source_discrepancy_audit_v1"
)
PROTECTED_LEGACY_OUTPUT_SHA256 = {
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_single_beam"
    / "report.md": (
        "BA293BD5EB88EF0F130D28873F14BC48864BF4C674B43DB8857977B80F7ACE99"
    ),
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_single_beam"
    / "bending_modes.png": (
        "90FC11CE1806CB4F9646964CEC4F00110C8D41BF4A3E17B3B907ADC3F9D971BD"
    ),
    ROOT
    / "results"
    / "laminated_beams"
    / "reddy_symmetric_single_beam"
    / "combined_spectrum.png": (
        "E229EF0EBDCCC571B52DF95E565C54D0C99E2D610B2DD89A85FEFB6525653E35"
    ),
}

REQUIRED_OUTPUTS = (
    "ritz_convergence.csv",
    "transfer_ritz_comparison.csv",
    "rotary_inertia_monotonicity.csv",
    "source_inconsistency_inventory.csv",
    "cc_characteristic_check.csv",
    "cf_characteristic_check.csv",
    "report.md",
    "run_manifest.json",
)

RITZ_ORDERS = (4, 6, 8, 10, 12)
RITZ_GUARD_ORDER = 16
GAUSS_LEGENDRE_ORDER = 48
K_SOURCE = 5.0 / 6.0
K_PROVENANCE = "INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY"
NUMERICAL_ALLOWANCE = 5.0e-9

TRANSFER_RITZ_RELATIVE_TOLERANCE = 1.0e-8
RITZ_CONVERGENCE_RELATIVE_TOLERANCE = 1.0e-8
ALGEBRAIC_RESIDUAL_TOLERANCE = 1.0e-10
ENERGY_IDENTITY_TOLERANCE = 1.0e-9
NATURAL_BOUNDARY_RESIDUAL_TOLERANCE = 1.0e-8
ESSENTIAL_BOUNDARY_RESIDUAL_TOLERANCE = 1.0e-12
TRANSFER_SINGULAR_RATIO_TOLERANCE = 1.0e-9
RI_MONOTONICITY_RELATIVE_TOLERANCE = 1.0e-12

STATUS_MODEL = "RLB-0A-MATHEMATICAL-MODEL"
STATUS_RITZ = "RLB-0A-INDEPENDENT-RITZ"
STATUS_RI = "RLB-0A-ROTARY-INERTIA-MONOTONICITY"
STATUS_TABLE = "RLB-0A-TABLE-4.3.3-REPRODUCTION"
STATUS_OVERALL = "OVERALL-SCIENTIFIC"

DISCREPANCY_CLASSES = frozenset(
    {
        "MATCHES_SOURCE_WITHIN_PRINTED_TOLERANCE",
        "SOURCE_ROUNDING_OR_HIDDEN_INTERMEDIATE_ROUNDING",
        "SOURCE_TABLE_VIOLATES_ROTARY_INERTIA_MONOTONICITY",
        "SOURCE_TABLE_DISAGREES_WITH_TRANSFER_AND_RITZ",
        "PRINTED_CHARACTERISTIC_DIAGNOSTIC_ROUNDING",
        "UNRESOLVED_MODEL_DISCREPANCY",
    }
)

FULL_PRECISION_CLASSICAL_ROOTS = {
    "CC": 4.730040744862704,
    "CF": 1.875104068711961,
}

PROBLEM_4_14_DIAGNOSTIC = {
    "status": "UNDEFINED_SYMBOL_IN_PRINTED_PROBLEM_4_14",
    "printed_page": 237,
    "pdf_index_zero_based": 259,
    "pdf_index_one_based": 260,
    "rotary_inertia_statement": "rotary inertia neglected",
    "undefined_symbol": "beta",
    "policy": (
        "Non-executable secondary diagnostic; beta is not fitted or evaluated."
    ),
    "R1": (
        "beta*I0_tilde*omega^2 + K*G_xz_b*b*h "
        "- lambda^2*E_xx_b*I_yy"
    ),
    "R2": (
        "beta*I0_tilde*omega^2 + K*G_xz_b*b*h "
        "+ mu^2*E_xx_b*I_yy"
    ),
    "S1": "I0_tilde*omega^2 - lambda^2*K*G_xz_b*b*h",
    "S2": "I0_tilde*omega^2 + mu^2*K*G_xz_b*b*h",
}


@dataclass(frozen=True)
class BenchmarkCase:
    """One published direct-FSDT Table 4.3.3 row."""

    case_id: str
    laminate_id: str
    benchmark_tier: str
    printed_label: str
    used_label: str
    boundary_condition: BoundaryCondition
    a_over_h: int
    rotary_inertia: bool
    row_role: str
    source_value: float | None
    source_status: str
    printed_tolerance: float | None
    K: float
    K_provenance: str
    length: float
    laminate: Any
    properties: rlb.BeamProperties

    @property
    def frequency_scale(self) -> float:
        """Multiplier from dimensional omega to the source omega-bar."""

        return float(
            self.length**2
            * math.sqrt(
                self.laminate.I0
                / self.laminate.thickness**3
            )
        )


@dataclass(frozen=True)
class RitzMatrices:
    """Dimensionless Galerkin matrices and their exact two-field blocks."""

    order: int
    boundary_condition: BoundaryCondition
    gauss_order: int
    gamma: float
    rotary_mass_ratio: float
    stiffness: FloatArray
    mass: FloatArray
    K_ww: FloatArray
    K_wpsi: FloatArray
    K_psiw: FloatArray
    K_psipsi: FloatArray
    M_ww: FloatArray
    M_psipsi: FloatArray
    stiffness_symmetry_residual: float
    mass_symmetry_residual: float
    stiffness_minimum_eigenvalue: float
    stiffness_spd: bool
    full_mass_minimum_eigenvalue: float | None
    full_mass_spd: bool | None
    translational_mass_minimum_eigenvalue: float
    translational_mass_spd: bool
    psi_mass_zero_residual: float
    j0_mass_structure_pass: bool


@dataclass(frozen=True)
class RitzMode:
    """Fundamental Ritz eigenpair with audit-quality residuals."""

    order: int
    boundary_condition: BoundaryCondition
    gauss_order: int
    eigenvalue: float
    rayleigh_eigenvalue: float
    rayleigh_eigenvalue_relative_difference: float
    omega_squared: float
    omega: float
    coefficients: FloatArray
    static_condensation_used: bool
    static_condensation_condition_number: float | None
    static_condensation_residual: float | None
    stiffness_condition_number: float
    mass_condition_number: float | None
    condensed_stiffness_condition_number: float | None
    generalized_eigen_residual: float
    matrix_norm_backward_residual: float
    recovered_full_mode_residual: float
    energy_identity_relative: float
    essential_boundary_residual: float
    natural_boundary_residual: float
    modal_mass: float
    strain_energy_twice: float
    kinetic_energy_twice_times_lambda: float
    gamma: float
    stiffness_symmetry_residual: float
    mass_symmetry_residual: float
    stiffness_minimum_eigenvalue: float
    stiffness_spd: bool
    full_mass_minimum_eigenvalue: float | None
    full_mass_spd: bool | None
    translational_mass_minimum_eigenvalue: float
    translational_mass_spd: bool
    psi_mass_zero_residual: float
    j0_mass_structure_pass: bool
    zero_mode_count: int
    lowest_positive_eigenvalue: float


@dataclass(frozen=True)
class TransferResult:
    """Independent physical transfer-matrix fundamental root."""

    omega: float
    determinant: float
    sigma_min: float
    sigma_max: float
    sigma_ratio: float
    condition_number: float
    boundary_residual: float | None
    accepted: bool
    detection: str
    bracket_left: float | None = None
    bracket_right: float | None = None
    historical_status: str | None = None


@dataclass(frozen=True)
class CaseAudit:
    """All independent data used to classify one published direct row."""

    case: BenchmarkCase
    modes: dict[int, RitzMode]
    legacy_transfer: TransferResult
    transfer: TransferResult
    ritz_convergence_relative: float
    guard_relative: float | None
    converged_order: int
    used_convergence_relative: float
    ritz_status: str
    model_status: str

    @property
    def converged_mode(self) -> RitzMode:
        return self.modes[self.converged_order]


def _positive(value: float, name: str) -> float:
    out = float(value)
    if not math.isfinite(out) or out <= 0.0:
        raise ValueError(f"{name} must be finite and positive.")
    return out


def _canonical_boundary_condition(value: str) -> BoundaryCondition:
    key = str(value).strip().upper().replace("-", "_")
    aliases = {
        "HH": "HH",
        "HINGED_HINGED": "HH",
        "CC": "CC",
        "CLAMPED_CLAMPED": "CC",
        "CF": "CF",
        "CLAMPED_FREE": "CF",
    }
    try:
        return aliases[key]  # type: ignore[return-value]
    except KeyError as error:
        raise ValueError(f"Unsupported bending boundary condition: {value!r}.") from error


def _relative_difference(value: float, reference: float) -> float:
    return abs(float(value) - float(reference)) / max(
        abs(float(reference)), sys.float_info.min
    )


def _matrix_relative_residual(matrix: FloatArray, vector: FloatArray) -> float:
    numerator = float(np.linalg.norm(matrix @ vector))
    denominator = max(
        float(np.linalg.norm(matrix) * np.linalg.norm(vector)),
        sys.float_info.min,
    )
    return numerator / denominator


def _symmetry_residual(matrix: FloatArray) -> float:
    return float(
        np.linalg.norm(matrix - matrix.T)
        / max(np.linalg.norm(matrix), sys.float_info.min)
    )


def _field_multiplier(
    boundary_condition: BoundaryCondition,
    field: Literal["w", "psi"],
    xi: FloatArray,
) -> tuple[FloatArray, FloatArray]:
    """Return the essential-BC multiplier and its xi derivative."""

    if boundary_condition == "HH":
        if field == "w":
            return xi * (1.0 - xi), 1.0 - 2.0 * xi
        return np.ones_like(xi), np.zeros_like(xi)
    if boundary_condition == "CC":
        return xi * (1.0 - xi), 1.0 - 2.0 * xi
    return xi, np.ones_like(xi)


def shifted_legendre_essential_basis(
    order: int,
    xi: Sequence[float] | FloatArray,
    boundary_condition: str,
    field: Literal["w", "psi"],
) -> tuple[FloatArray, FloatArray]:
    """Evaluate a shifted-Legendre trial basis and its xi derivative.

    The returned arrays have shape ``(order, n_points)``.  Essential boundary
    conditions are imposed by a polynomial multiplier; natural conditions are
    left to the variational statement.
    """

    if not isinstance(order, int) or order < 1:
        raise ValueError("order must be a positive integer.")
    bc = _canonical_boundary_condition(boundary_condition)
    points = np.asarray(xi, dtype=float)
    if points.ndim != 1 or not np.all(np.isfinite(points)):
        raise ValueError("xi must be a finite one-dimensional array.")
    if np.any(points < 0.0) or np.any(points > 1.0):
        raise ValueError("xi points must lie in [0, 1].")
    multiplier, multiplier_derivative = _field_multiplier(bc, field, points)
    coordinate = 2.0 * points - 1.0
    values = np.empty((order, points.size), dtype=float)
    derivatives = np.empty_like(values)
    for degree in range(order):
        polynomial = Legendre.basis(degree)
        raw = polynomial(coordinate)
        raw_derivative = 2.0 * polynomial.deriv()(coordinate)
        values[degree] = multiplier * raw
        derivatives[degree] = (
            multiplier_derivative * raw + multiplier * raw_derivative
        )
    return values, derivatives


def assemble_ritz_matrices(
    properties: rlb.BeamProperties,
    length: float,
    boundary_condition: str,
    order: int,
    *,
    gauss_order: int = GAUSS_LEGENDRE_ORDER,
) -> RitzMatrices:
    """Assemble the independent polynomial FSDT generalized eigenproblem."""

    L = _positive(length, "length")
    if not isinstance(order, int) or order < 2:
        raise ValueError("order must be an integer of at least two.")
    if gauss_order != GAUSS_LEGENDRE_ORDER:
        raise ValueError(
            f"This audit fixes Gauss--Legendre order at {GAUSS_LEGENDRE_ORDER}."
        )
    bc = _canonical_boundary_condition(boundary_condition)
    raw_nodes, raw_weights = np.polynomial.legendre.leggauss(gauss_order)
    xi = 0.5 * (raw_nodes + 1.0)
    weights = 0.5 * raw_weights
    w_values, w_derivatives = shifted_legendre_essential_basis(
        order, xi, bc, "w"
    )
    psi_values, psi_derivatives = shifted_legendre_essential_basis(
        order, xi, bc, "psi"
    )

    def inner(left: FloatArray, right: FloatArray) -> FloatArray:
        return (left * weights[None, :]) @ right.T

    gamma = float(properties.S * L**2 / properties.D)
    rotary_mass_ratio = float(properties.J / (properties.m * L**2))
    raw_K_ww = gamma * inner(w_derivatives, w_derivatives)
    raw_K_wpsi = gamma * inner(w_derivatives, psi_values)
    raw_K_psiw = gamma * inner(psi_values, w_derivatives)
    raw_K_psipsi = inner(psi_derivatives, psi_derivatives) + gamma * inner(
        psi_values, psi_values
    )
    raw_M_ww = inner(w_values, w_values)
    raw_M_psipsi = rotary_mass_ratio * inner(psi_values, psi_values)
    zeros = np.zeros_like(raw_M_ww)
    raw_stiffness = np.block(
        [[raw_K_ww, raw_K_wpsi], [raw_K_psiw, raw_K_psipsi]]
    )
    raw_mass = np.block([[raw_M_ww, zeros], [zeros, raw_M_psipsi]])
    stiffness_symmetry_residual = _symmetry_residual(raw_stiffness)
    mass_symmetry_residual = _symmetry_residual(raw_mass)
    # After auditing both independently integrated coupling blocks, retain the
    # energy-Hessian representation used by the symmetric variational form.
    energy_stiffness = np.block(
        [
            [raw_K_ww, raw_K_wpsi],
            [raw_K_wpsi.T.copy(), raw_K_psipsi],
        ]
    )
    stiffness = 0.5 * (energy_stiffness + energy_stiffness.T)
    mass = 0.5 * (raw_mass + raw_mass.T)
    K_ww = stiffness[:order, :order]
    K_wpsi = stiffness[:order, order:]
    K_psiw = stiffness[order:, :order]
    K_psipsi = stiffness[order:, order:]
    M_ww = mass[:order, :order]
    M_psipsi = mass[order:, order:]
    stiffness_eigenvalues = np.linalg.eigvalsh(stiffness)
    translational_mass_eigenvalues = np.linalg.eigvalsh(M_ww)
    full_mass_eigenvalues = (
        None if properties.J == 0.0 else np.linalg.eigvalsh(mass)
    )
    psi_mass_zero_residual = float(
        np.linalg.norm(M_psipsi)
        / max(np.linalg.norm(M_ww), sys.float_info.min)
    )
    return RitzMatrices(
        order=order,
        boundary_condition=bc,
        gauss_order=gauss_order,
        gamma=gamma,
        rotary_mass_ratio=rotary_mass_ratio,
        stiffness=stiffness,
        mass=mass,
        K_ww=K_ww,
        K_wpsi=K_wpsi,
        K_psiw=K_psiw,
        K_psipsi=K_psipsi,
        M_ww=M_ww,
        M_psipsi=M_psipsi,
        stiffness_symmetry_residual=stiffness_symmetry_residual,
        mass_symmetry_residual=mass_symmetry_residual,
        stiffness_minimum_eigenvalue=float(stiffness_eigenvalues[0]),
        stiffness_spd=bool(stiffness_eigenvalues[0] > 0.0),
        full_mass_minimum_eigenvalue=(
            None
            if full_mass_eigenvalues is None
            else float(full_mass_eigenvalues[0])
        ),
        full_mass_spd=(
            None
            if full_mass_eigenvalues is None
            else bool(full_mass_eigenvalues[0] > 0.0)
        ),
        translational_mass_minimum_eigenvalue=float(
            translational_mass_eigenvalues[0]
        ),
        translational_mass_spd=bool(translational_mass_eigenvalues[0] > 0.0),
        psi_mass_zero_residual=psi_mass_zero_residual,
        j0_mass_structure_pass=bool(
            properties.J != 0.0 or psi_mass_zero_residual == 0.0
        ),
    )


def evaluate_ritz_mode(
    mode: RitzMode,
    xi: Sequence[float] | FloatArray,
) -> FloatArray:
    """Return dimensionless physical state ``[W, psi, Qbar, Mbar]``.

    ``Qbar = Q*length**2/D = gamma*(W_,xi + psi)`` and
    ``Mbar = M*length/D = psi_,xi``.
    """

    points = np.asarray(xi, dtype=float)
    w_values, w_derivatives = shifted_legendre_essential_basis(
        mode.order, points, mode.boundary_condition, "w"
    )
    psi_values, psi_derivatives = shifted_legendre_essential_basis(
        mode.order, points, mode.boundary_condition, "psi"
    )
    w_coefficients = mode.coefficients[: mode.order]
    psi_coefficients = mode.coefficients[mode.order :]
    W = w_coefficients @ w_values
    W_derivative = w_coefficients @ w_derivatives
    psi = psi_coefficients @ psi_values
    psi_derivative = psi_coefficients @ psi_derivatives
    Qbar = mode.gamma * (W_derivative + psi)
    return np.column_stack([W, psi, Qbar, psi_derivative])


def _boundary_residuals(mode: RitzMode) -> tuple[float, float]:
    states = evaluate_ritz_mode(mode, np.array([0.0, 1.0]))
    denominator = max(float(np.linalg.norm(states)), sys.float_info.min)
    if mode.boundary_condition == "HH":
        essential = states[:, 0]
        natural = states[:, 3]
    elif mode.boundary_condition == "CC":
        essential = states[:, :2].reshape(-1)
        natural = np.zeros(1)
    else:
        essential = states[0, :2]
        natural = states[1, 2:]
    return (
        float(np.linalg.norm(essential) / denominator),
        float(np.linalg.norm(natural) / denominator),
    )


def solve_ritz_fundamental(
    properties: rlb.BeamProperties,
    length: float,
    boundary_condition: str,
    order: int,
    *,
    gauss_order: int = GAUSS_LEGENDRE_ORDER,
) -> RitzMode:
    """Solve the fundamental Ritz mode without any transfer/root machinery."""

    matrices = assemble_ritz_matrices(
        properties,
        length,
        boundary_condition,
        order,
        gauss_order=gauss_order,
    )
    static_condensation_used = properties.J == 0.0
    condensation_condition: float | None = None
    condensation_residual: float | None = None
    condensed_condition: float | None = None
    mass_condition: float | None = None

    def solve_spd_generalized(
        stiffness: FloatArray, mass: FloatArray
    ) -> tuple[FloatArray, FloatArray]:
        """Solve by an explicit mass-orthonormal congruence."""

        factor = cholesky(mass, lower=True, check_finite=True)
        inverse_factor = solve_triangular(
            factor,
            np.eye(factor.shape[0]),
            lower=True,
            check_finite=True,
        )
        orthonormal_stiffness = inverse_factor @ stiffness @ inverse_factor.T
        orthonormal_stiffness = 0.5 * (
            orthonormal_stiffness + orthonormal_stiffness.T
        )
        values, orthonormal_vectors = eigh(
            orthonormal_stiffness,
            check_finite=True,
        )
        physical_vectors = solve_triangular(
            factor.T,
            orthonormal_vectors,
            lower=False,
            check_finite=True,
        )
        return values, physical_vectors

    if static_condensation_used:
        condensation_condition = float(np.linalg.cond(matrices.K_psipsi))
        rotation_from_w = np.linalg.solve(
            matrices.K_psipsi, matrices.K_psiw
        )
        condensed = matrices.K_ww - matrices.K_wpsi @ rotation_from_w
        condensed = 0.5 * (condensed + condensed.T)
        condensed_condition = float(np.linalg.cond(condensed))
        eigenvalues, eigenvectors = solve_spd_generalized(
            condensed, matrices.M_ww
        )
        zero_tolerance = 100.0 * np.finfo(float).eps * max(
            float(np.max(np.abs(eigenvalues))), 1.0
        )
        positive = np.flatnonzero(eigenvalues > zero_tolerance)
        if positive.size == 0:
            raise RuntimeError("The condensed Ritz problem has no positive root.")
        index = int(positive[0])
        w_coefficients = eigenvectors[:, index]
        psi_coefficients = -rotation_from_w @ w_coefficients
        coefficients = np.concatenate([w_coefficients, psi_coefficients])
        eigenvalue = float(eigenvalues[index])
        static_vector = (
            matrices.K_psiw @ w_coefficients
            + matrices.K_psipsi @ psi_coefficients
        )
        static_denominator = max(
            float(
                np.linalg.norm(matrices.K_psiw) * np.linalg.norm(w_coefficients)
                + np.linalg.norm(matrices.K_psipsi)
                * np.linalg.norm(psi_coefficients)
            ),
            sys.float_info.min,
        )
        condensation_residual = float(
            np.linalg.norm(static_vector) / static_denominator
        )
    else:
        mass_condition = float(np.linalg.cond(matrices.mass))
        eigenvalues, eigenvectors = solve_spd_generalized(
            matrices.stiffness, matrices.mass
        )
        zero_tolerance = 100.0 * np.finfo(float).eps * max(
            float(np.max(np.abs(eigenvalues))), 1.0
        )
        positive = np.flatnonzero(eigenvalues > zero_tolerance)
        if positive.size == 0:
            raise RuntimeError("The full Ritz problem has no positive root.")
        index = int(positive[0])
        coefficients = eigenvectors[:, index]
        eigenvalue = float(eigenvalues[index])

    pivot = int(np.argmax(np.abs(coefficients)))
    if coefficients[pivot] < 0.0:
        coefficients = -coefficients
    modal_mass = float(coefficients @ matrices.mass @ coefficients)
    if not math.isfinite(modal_mass) or modal_mass <= 0.0:
        raise RuntimeError("The recovered Ritz mode has non-positive modal mass.")
    coefficients = coefficients / math.sqrt(modal_mass)
    modal_mass = float(coefficients @ matrices.mass @ coefficients)
    strain = float(coefficients @ matrices.stiffness @ coefficients)
    rayleigh_eigenvalue = float(strain / modal_mass)
    kinetic_lambda = float(eigenvalue * modal_mass)
    full_residual_vector = (
        matrices.stiffness @ coefficients
        - eigenvalue * (matrices.mass @ coefficients)
    )
    residual_denominator = max(
        float(
            np.linalg.norm(matrices.stiffness @ coefficients)
            + abs(eigenvalue) * np.linalg.norm(matrices.mass @ coefficients)
        ),
        sys.float_info.min,
    )
    full_residual = float(np.linalg.norm(full_residual_vector) / residual_denominator)
    matrix_norm_denominator = max(
        float(
            (
                np.linalg.norm(matrices.stiffness)
                + abs(eigenvalue) * np.linalg.norm(matrices.mass)
            )
            * np.linalg.norm(coefficients)
        ),
        sys.float_info.min,
    )
    matrix_norm_backward_residual = float(
        np.linalg.norm(full_residual_vector) / matrix_norm_denominator
    )
    energy_residual = abs(strain - kinetic_lambda) / max(
        abs(strain), abs(kinetic_lambda), sys.float_info.min
    )
    omega_squared = float(
        eigenvalue * properties.D / (properties.m * float(length) ** 4)
    )
    if not math.isfinite(omega_squared) or omega_squared <= 0.0:
        raise RuntimeError("The Ritz eigenfrequency is not finite and positive.")

    provisional = RitzMode(
        order=order,
        boundary_condition=matrices.boundary_condition,
        gauss_order=gauss_order,
        eigenvalue=eigenvalue,
        rayleigh_eigenvalue=rayleigh_eigenvalue,
        rayleigh_eigenvalue_relative_difference=_relative_difference(
            rayleigh_eigenvalue, eigenvalue
        ),
        omega_squared=omega_squared,
        omega=math.sqrt(omega_squared),
        coefficients=coefficients,
        static_condensation_used=static_condensation_used,
        static_condensation_condition_number=condensation_condition,
        static_condensation_residual=condensation_residual,
        stiffness_condition_number=float(np.linalg.cond(matrices.stiffness)),
        mass_condition_number=mass_condition,
        condensed_stiffness_condition_number=condensed_condition,
        generalized_eigen_residual=full_residual,
        matrix_norm_backward_residual=matrix_norm_backward_residual,
        recovered_full_mode_residual=full_residual,
        energy_identity_relative=float(energy_residual),
        essential_boundary_residual=0.0,
        natural_boundary_residual=0.0,
        modal_mass=modal_mass,
        strain_energy_twice=strain,
        kinetic_energy_twice_times_lambda=kinetic_lambda,
        gamma=matrices.gamma,
        stiffness_symmetry_residual=matrices.stiffness_symmetry_residual,
        mass_symmetry_residual=matrices.mass_symmetry_residual,
        stiffness_minimum_eigenvalue=matrices.stiffness_minimum_eigenvalue,
        stiffness_spd=matrices.stiffness_spd,
        full_mass_minimum_eigenvalue=matrices.full_mass_minimum_eigenvalue,
        full_mass_spd=matrices.full_mass_spd,
        translational_mass_minimum_eigenvalue=matrices.translational_mass_minimum_eigenvalue,
        translational_mass_spd=matrices.translational_mass_spd,
        psi_mass_zero_residual=matrices.psi_mass_zero_residual,
        j0_mass_structure_pass=matrices.j0_mass_structure_pass,
        zero_mode_count=int(np.count_nonzero(eigenvalues <= zero_tolerance)),
        lowest_positive_eigenvalue=float(eigenvalues[int(positive[0])]),
    )
    essential, natural = _boundary_residuals(provisional)
    return replace(
        provisional,
        essential_boundary_residual=essential,
        natural_boundary_residual=natural,
    )


def solve_ritz_sequence(
    properties: rlb.BeamProperties,
    length: float,
    boundary_condition: str,
    orders: Iterable[int] = RITZ_ORDERS,
) -> dict[int, RitzMode]:
    """Solve a deterministic sequence of independent Ritz spaces."""

    requested = tuple(int(order) for order in orders)
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("orders must be a non-empty sequence without duplicates.")
    return {
        order: solve_ritz_fundamental(
            properties,
            length,
            boundary_condition,
            order,
        )
        for order in requested
    }


def _load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as stream:
        return json.load(
            stream,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ValueError(f"Non-finite JSON constant is forbidden: {value}")
            ),
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def _require_sha256(path: Path, expected: str, label: str) -> str:
    """Reject a changed frozen input before consuming any of its values."""

    if not path.is_file():
        raise FileNotFoundError(f"Frozen {label} is missing: {path}")
    actual = _sha256(path)
    if actual != expected:
        raise RuntimeError(
            f"Frozen {label} SHA-256 mismatch: {actual} != {expected}."
        )
    return actual


def _require_protected_legacy_outputs() -> dict[str, str]:
    verified: dict[str, str] = {}
    for path, expected in PROTECTED_LEGACY_OUTPUT_SHA256.items():
        actual = _require_sha256(path, expected, "protected legacy output")
        verified[str(path.relative_to(ROOT)).replace("\\", "/")] = actual
    return verified


def _source_key(record: dict[str, Any]) -> tuple[str, str, int, str]:
    return (
        str(record["laminate_id"]),
        str(record["boundary_condition"]),
        int(record["a_over_h"]),
        str(record["row_role"]),
    )


def load_legacy_transfer_inventory(
    contract: dict[str, Any],
    path: Path = LEGACY_TRANSFER_INVENTORY_PATH,
) -> tuple[dict[tuple[str, str, int, str], dict[str, str]], dict[str, Any]]:
    """Read, validate, and hash the immutable RLB-0 transfer inventory.

    The source keys, source values, printed tolerances, and explicit ``K`` are
    cross-checked against the frozen JSON before any transfer value is used.
    The CSV is never written by this audit.
    """

    _require_sha256(
        path,
        LEGACY_TRANSFER_INVENTORY_SHA256,
        "RLB-0 transfer inventory",
    )
    with path.open("r", encoding="utf-8", newline="") as stream:
        csv_rows = list(csv.DictReader(stream))
    direct_rows = [
        row
        for row in csv_rows
        if row.get("row_role")
        in {"fsdt_frequency_with_RI", "fsdt_frequency_without_RI"}
    ]
    inventory: dict[tuple[str, str, int, str], dict[str, str]] = {}
    for row in direct_rows:
        key = _source_key(row)
        if key in inventory:
            raise RuntimeError(f"Duplicate direct-FSDT transfer key: {key!r}.")
        inventory[key] = row

    expected_records = [
        record
        for record in contract["records"]
        if record["row_role"]
        in {"fsdt_frequency_with_RI", "fsdt_frequency_without_RI"}
    ]
    mismatches: list[str] = []
    for record in expected_records:
        key = _source_key(record)
        row = inventory.get(key)
        if row is None:
            mismatches.append(f"missing key {key!r}")
            continue
        source_value = record["source_value"]
        csv_source = row.get("source_value", "").strip()
        if source_value is None:
            if csv_source:
                mismatches.append(f"source-null mismatch at {key!r}")
        elif not math.isclose(
            float(csv_source), float(source_value), rel_tol=0.0, abs_tol=0.0
        ):
            mismatches.append(f"source-value mismatch at {key!r}")
        exact_checks = {
            "source_status": str(record["source_status"]),
            "K_provenance": str(record["K_provenance"]),
        }
        for field, expected in exact_checks.items():
            if row.get(field, "") != expected:
                mismatches.append(f"{field} mismatch at {key!r}")
        numeric_checks = {
            "K": record["K"],
            "printed_tolerance": record["printed_tolerance"],
        }
        for field, expected_value in numeric_checks.items():
            csv_value = row.get(field, "").strip()
            if expected_value is None:
                if csv_value:
                    mismatches.append(f"{field} null mismatch at {key!r}")
                continue
            expected = float(expected_value)
            if not math.isclose(
                float(csv_value or "nan"), expected, rel_tol=0.0, abs_tol=0.0
            ):
                mismatches.append(f"{field} mismatch at {key!r}")
    expected_keys = {_source_key(record) for record in expected_records}
    extra_keys = sorted(set(inventory) - expected_keys)
    if extra_keys:
        mismatches.append(f"unexpected transfer keys: {extra_keys!r}")
    if mismatches:
        raise RuntimeError(
            "The authoritative transfer inventory does not match the frozen "
            "source contract: " + "; ".join(mismatches)
        )
    metadata = {
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "sha256": _sha256(path),
        "direct_row_count": len(direct_rows),
        "published_direct_row_count": sum(
            bool(row.get("source_value", "").strip()) for row in direct_rows
        ),
        "source_null_project_diagnostic_count": sum(
            not bool(row.get("source_value", "").strip()) for row in direct_rows
        ),
        "source_key_validation": "PASS",
    }
    return inventory, metadata


def _validate_direct_case_inventory(cases: Sequence[BenchmarkCase]) -> None:
    """Enforce the frozen 72/54/36/24 Table 4.3.3 audit contract."""

    failures: list[str] = []
    expected_laminates = {
        "0_deg",
        "90_deg",
        "cross_ply_0_90_s",
        "angle_ply_45_minus45_s",
    }
    if {case.laminate_id for case in cases} != expected_laminates:
        failures.append("laminate IDs")
    if {case.boundary_condition for case in cases} != {"HH", "CC", "CF"}:
        failures.append("boundary-condition IDs")
    if {case.a_over_h for case in cases} != {10, 20, 100}:
        failures.append("a/h values")
    if len({case.case_id for case in cases}) != 72:
        failures.append("unique case IDs")
    if sum(case.source_value is not None for case in cases) != 54:
        failures.append("54 published direct rows")
    if sum(case.source_value is None for case in cases) != 18:
        failures.append("18 source-null diagnostics")
    if sum(case.rotary_inertia for case in cases) != 36:
        failures.append("36 with-RI cases")
    if sum(not case.rotary_inertia for case in cases) != 36:
        failures.append("36 no-RI cases")
    for boundary_condition in ("HH", "CC", "CF"):
        if sum(
            case.boundary_condition == boundary_condition for case in cases
        ) != 24:
            failures.append(f"24 {boundary_condition} model cases")
        if sum(
            case.boundary_condition == boundary_condition
            and case.source_value is not None
            for case in cases
        ) != 18:
            failures.append(f"18 published {boundary_condition} rows")
    ri_pairs: dict[tuple[str, str, int], set[bool]] = {}
    for case in cases:
        key = (case.laminate_id, case.boundary_condition, case.a_over_h)
        ri_pairs.setdefault(key, set()).add(case.rotary_inertia)
    if len(ri_pairs) != 36 or any(
        states != {False, True} for states in ri_pairs.values()
    ):
        failures.append("36 complete RI/no-RI pairs")
    if failures:
        raise RuntimeError(
            "Frozen direct-case inventory violates the audit contract: "
            + "; ".join(failures)
        )


def load_all_direct_cases(
    source_contract_path: Path = SOURCE_CONTRACT_PATH,
) -> tuple[dict[str, Any], list[BenchmarkCase]]:
    """Load all direct rows, retaining source-null project diagnostics."""

    _require_sha256(
        source_contract_path,
        SOURCE_CONTRACT_SHA256,
        "Table 4.3.3 source contract",
    )
    contract = _load_json(source_contract_path)
    material = rlb.LaminaMaterial(
        E1=25.0,
        E2=1.0,
        G12=0.5,
        G13=0.5,
        G23=0.2,
        nu12=0.25,
        rho=1.0,
        name="Reddy Eq. (4.2.25) dimensionless realization",
    )
    laminates: dict[str, tuple[Any, rlb.BeamProperties]] = {}
    for laminate_id, definition in contract["laminate_definitions"].items():
        stack = tuple(float(angle) for angle in definition["stack_degrees"])
        ply_thickness = 1.0 / len(stack)
        laminate = rlb.integrate_laminate(
            [
                rlb.Ply(material, theta_deg, ply_thickness)
                for theta_deg in stack
            ]
        )
        properties = rlb.reduce_to_beam_properties(
            laminate, width=1.0, K=K_SOURCE
        )
        laminates[laminate_id] = laminate, properties

    cases: list[BenchmarkCase] = []
    for record in contract["records"]:
        row_role = str(record["row_role"])
        if row_role not in {
            "fsdt_frequency_with_RI",
            "fsdt_frequency_without_RI",
        }:
            continue
        laminate, base_properties = laminates[record["laminate_id"]]
        with_ri = row_role == "fsdt_frequency_with_RI"
        properties = base_properties if with_ri else rlb.without_rotary_inertia(
            base_properties
        )
        ratio = int(record["a_over_h"])
        bc = _canonical_boundary_condition(record["boundary_condition"])
        ri_label = "with_RI" if with_ri else "without_RI"
        cases.append(
            BenchmarkCase(
                case_id=f"{record['laminate_id']}__{bc}__a_h_{ratio}__{ri_label}",
                laminate_id=str(record["laminate_id"]),
                benchmark_tier=str(record["benchmark_tier"]),
                printed_label=str(record["printed_label"]),
                used_label=str(record["used_label"]),
                boundary_condition=bc,
                a_over_h=ratio,
                rotary_inertia=with_ri,
                row_role=row_role,
                source_value=(
                    None
                    if record["source_value"] is None
                    else float(record["source_value"])
                ),
                source_status=str(record["source_status"]),
                printed_tolerance=(
                    None
                    if record["printed_tolerance"] is None
                    else float(record["printed_tolerance"])
                ),
                K=float(record["K"]),
                K_provenance=str(record["K_provenance"]),
                length=float(ratio * laminate.thickness),
                laminate=laminate,
                properties=properties,
            )
        )
    expected = 72
    if len(cases) != expected:
        raise RuntimeError(
            f"Frozen source contract yielded {len(cases)} direct rows; expected {expected}."
        )
    _validate_direct_case_inventory(cases)
    return contract, cases


def load_published_direct_cases(
    source_contract_path: Path = SOURCE_CONTRACT_PATH,
) -> tuple[dict[str, Any], list[BenchmarkCase]]:
    """Load the 54 published direct-FSDT Table 4.3.3 rows."""

    contract, all_cases = load_all_direct_cases(source_contract_path)
    cases = [case for case in all_cases if case.source_value is not None]
    if len(cases) != 54:
        raise RuntimeError(
            f"Frozen source contract yielded {len(cases)} published direct rows; "
            "expected 54."
        )
    return contract, cases


def _transfer_result_from_row(row: dict[str, str]) -> TransferResult:
    required = (
        "root_omega",
        "scaled_determinant",
        "sigma_min",
        "sigma_max",
        "relative_singular_residual",
        "condition_number",
    )
    missing = [field for field in required if not row.get(field, "").strip()]
    if missing:
        raise RuntimeError(
            f"Transfer inventory row lacks diagnostics {missing!r}: "
            f"{_source_key(row)!r}."
        )
    root_status = row.get("root_status", "")
    sigma_ratio = float(row["relative_singular_residual"])
    return TransferResult(
        omega=float(row["root_omega"]),
        determinant=float(row["scaled_determinant"]),
        sigma_min=float(row["sigma_min"]),
        sigma_max=float(row["sigma_max"]),
        sigma_ratio=sigma_ratio,
        condition_number=float(row["condition_number"]),
        # The normalized smallest singular value is the minimum endpoint
        # residual over unit initial free-state vectors.
        boundary_residual=sigma_ratio,
        accepted=(
            root_status == "ACCEPTED"
            and math.isfinite(sigma_ratio)
            and sigma_ratio <= TRANSFER_SINGULAR_RATIO_TOLERANCE
        ),
        detection=f"immutable_RLB0_inventory:{root_status}",
        historical_status=row.get("status", ""),
    )


def _mode_quality_pass(mode: RitzMode) -> bool:
    return (
        mode.stiffness_symmetry_residual <= ALGEBRAIC_RESIDUAL_TOLERANCE
        and mode.mass_symmetry_residual <= ALGEBRAIC_RESIDUAL_TOLERANCE
        and mode.stiffness_spd
        and mode.translational_mass_spd
        and (mode.full_mass_spd is True or mode.static_condensation_used)
        and mode.j0_mass_structure_pass
        and mode.zero_mode_count == 0
        and mode.matrix_norm_backward_residual <= ALGEBRAIC_RESIDUAL_TOLERANCE
        and (
            not mode.static_condensation_used
            or mode.recovered_full_mode_residual
            <= ALGEBRAIC_RESIDUAL_TOLERANCE
        )
        and mode.rayleigh_eigenvalue_relative_difference
        <= ENERGY_IDENTITY_TOLERANCE
        and mode.energy_identity_relative <= ENERGY_IDENTITY_TOLERANCE
        and mode.essential_boundary_residual
        <= ESSENTIAL_BOUNDARY_RESIDUAL_TOLERANCE
        and mode.natural_boundary_residual
        <= NATURAL_BOUNDARY_RESIDUAL_TOLERANCE
        and (
            mode.static_condensation_residual is None
            or mode.static_condensation_residual
            <= ALGEBRAIC_RESIDUAL_TOLERANCE
        )
    )


def audit_published_case(
    case: BenchmarkCase,
    transfer_row: dict[str, str],
    *,
    include_guard: bool = False,
) -> CaseAudit:
    """Solve Ritz first, then compare with the immutable transfer inventory."""

    modes = solve_ritz_sequence(
        case.properties,
        case.length,
        case.boundary_condition,
        RITZ_ORDERS,
    )
    convergence = _relative_difference(modes[12].omega, modes[10].omega)
    guard_relative: float | None = None
    if (
        include_guard
        or convergence > RITZ_CONVERGENCE_RELATIVE_TOLERANCE
        or not _mode_quality_pass(modes[12])
    ):
        modes[RITZ_GUARD_ORDER] = solve_ritz_fundamental(
            case.properties,
            case.length,
            case.boundary_condition,
            RITZ_GUARD_ORDER,
        )
        guard_relative = _relative_difference(
            modes[RITZ_GUARD_ORDER].omega, modes[12].omega
        )

    use_guard = (
        convergence > RITZ_CONVERGENCE_RELATIVE_TOLERANCE
        or not _mode_quality_pass(modes[12])
    )
    converged_order = RITZ_GUARD_ORDER if use_guard else 12
    if use_guard and guard_relative is None:
        raise RuntimeError(f"N=16 guard was required but absent for {case.case_id}.")
    used_convergence = float(guard_relative) if use_guard else convergence
    converged_mode = modes[converged_order]
    if (
        used_convergence <= RITZ_CONVERGENCE_RELATIVE_TOLERANCE
        and _mode_quality_pass(converged_mode)
    ):
        ritz_status = "PASS"
    elif (
        guard_relative is not None
        and guard_relative <= RITZ_CONVERGENCE_RELATIVE_TOLERANCE
        and _mode_quality_pass(modes[RITZ_GUARD_ORDER])
    ):
        ritz_status = "PARTIAL_PASS"
    else:
        ritz_status = "FAIL"

    # Transfer data are consumed only after the Ritz order and eigenpair have
    # been selected from Ritz convergence/quality information alone.
    legacy_transfer = _transfer_result_from_row(transfer_row)
    transfer = refine_transfer_root_after_ritz(case, legacy_transfer)
    transfer_relative = _relative_difference(converged_mode.omega, transfer.omega)
    model_status = (
        "PASS"
        if (
            ritz_status == "PASS"
            and legacy_transfer.accepted
            and transfer.accepted
            and transfer_relative <= TRANSFER_RITZ_RELATIVE_TOLERANCE
        )
        else (
            "PARTIAL_PASS"
            if (
                ritz_status == "PARTIAL_PASS"
                and legacy_transfer.accepted
                and transfer.accepted
                and transfer_relative <= TRANSFER_RITZ_RELATIVE_TOLERANCE
            )
            else "FAIL"
        )
    )
    return CaseAudit(
        case=case,
        modes=modes,
        legacy_transfer=legacy_transfer,
        transfer=transfer,
        ritz_convergence_relative=convergence,
        guard_relative=guard_relative,
        converged_order=converged_order,
        used_convergence_relative=used_convergence,
        ritz_status=ritz_status,
        model_status=model_status,
    )


def audit_direct_cases(
    cases: Sequence[BenchmarkCase],
    transfer_inventory: dict[tuple[str, str, int, str], dict[str, str]],
    *,
    include_guard: bool = False,
) -> list[CaseAudit]:
    """Run the independent audit in frozen source-record order."""

    audits: list[CaseAudit] = []
    for case in cases:
        key = (
            case.laminate_id,
            case.boundary_condition,
            case.a_over_h,
            case.row_role,
        )
        try:
            transfer_row = transfer_inventory[key]
        except KeyError as error:
            raise RuntimeError(f"No transfer inventory row for {case.case_id}.") from error
        audits.append(
            audit_published_case(
                case,
                transfer_row,
                include_guard=include_guard,
            )
        )
    return audits


def audit_published_cases(
    cases: Sequence[BenchmarkCase],
    transfer_inventory: dict[tuple[str, str, int, str], dict[str, str]],
    *,
    include_guard: bool = False,
) -> list[CaseAudit]:
    """Compatibility name for callers passing only published direct rows."""

    if any(case.source_value is None for case in cases):
        raise ValueError("audit_published_cases received a source-null case.")
    return audit_direct_cases(
        cases,
        transfer_inventory,
        include_guard=include_guard,
    )


def _mode_frequency_bar(case: BenchmarkCase, mode: RitzMode) -> float:
    return float(mode.omega * case.frequency_scale)


def _transfer_frequency_bar(audit: CaseAudit) -> float:
    return float(audit.transfer.omega * audit.case.frequency_scale)


def _legacy_transfer_frequency_bar(audit: CaseAudit) -> float:
    return float(audit.legacy_transfer.omega * audit.case.frequency_scale)


def build_ritz_convergence_rows(audits: Sequence[CaseAudit]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for audit in audits:
        reference = audit.modes[12].omega
        ordered = sorted(audit.modes)
        for index, order in enumerate(ordered):
            mode = audit.modes[order]
            row_status = "PASS" if _mode_quality_pass(mode) else "FAIL"
            previous_order: int | str = "" if index == 0 else ordered[index - 1]
            relative_to_previous: float | str = ""
            convergence_status = "NOT_EVALUATED"
            if previous_order != "":
                relative_to_previous = _relative_difference(
                    mode.omega, audit.modes[int(previous_order)].omega
                )
                convergence_status = (
                    "PASS"
                    if float(relative_to_previous)
                    <= RITZ_CONVERGENCE_RELATIVE_TOLERANCE
                    else "FAIL"
                )
            rows.append(
                {
                    "case_id": audit.case.case_id,
                    "laminate_id": audit.case.laminate_id,
                    "boundary_condition": audit.case.boundary_condition,
                    "a_over_h": audit.case.a_over_h,
                    "rotary_inertia": audit.case.rotary_inertia,
                    "row_role": audit.case.row_role,
                    "order": order,
                    "guard_only": order == RITZ_GUARD_ORDER,
                    "gauss_order": mode.gauss_order,
                    "omega": mode.omega,
                    "omega_bar": _mode_frequency_bar(audit.case, mode),
                    "solver_eigenvalue": mode.eigenvalue,
                    "rayleigh_eigenvalue": mode.rayleigh_eigenvalue,
                    "rayleigh_eigenvalue_relative_difference": mode.rayleigh_eigenvalue_relative_difference,
                    "relative_to_N12": _relative_difference(mode.omega, reference),
                    "previous_order": previous_order,
                    "relative_to_previous_order": relative_to_previous,
                    "convergence_status": convergence_status,
                    "selected_for_comparison": order == audit.converged_order,
                    "generalized_eigen_residual": mode.generalized_eigen_residual,
                    "matrix_norm_backward_residual": mode.matrix_norm_backward_residual,
                    "recovered_full_mode_residual": mode.recovered_full_mode_residual,
                    "static_condensation_used": mode.static_condensation_used,
                    "static_condensation_condition_number": mode.static_condensation_condition_number,
                    "static_condensation_residual": mode.static_condensation_residual,
                    "energy_identity_relative": mode.energy_identity_relative,
                    "essential_boundary_residual": mode.essential_boundary_residual,
                    "natural_boundary_residual": mode.natural_boundary_residual,
                    "modal_mass": mode.modal_mass,
                    "stiffness_condition_number": mode.stiffness_condition_number,
                    "mass_condition_number": mode.mass_condition_number,
                    "condensed_stiffness_condition_number": mode.condensed_stiffness_condition_number,
                    "stiffness_symmetry_residual": mode.stiffness_symmetry_residual,
                    "mass_symmetry_residual": mode.mass_symmetry_residual,
                    "stiffness_minimum_eigenvalue": mode.stiffness_minimum_eigenvalue,
                    "stiffness_spd": mode.stiffness_spd,
                    "full_mass_minimum_eigenvalue": mode.full_mass_minimum_eigenvalue,
                    "full_mass_spd": mode.full_mass_spd,
                    "translational_mass_minimum_eigenvalue": mode.translational_mass_minimum_eigenvalue,
                    "translational_mass_spd": mode.translational_mass_spd,
                    "psi_mass_zero_residual": mode.psi_mass_zero_residual,
                    "j0_mass_structure_pass": mode.j0_mass_structure_pass,
                    "zero_mode_count": mode.zero_mode_count,
                    "lowest_positive_eigenvalue": mode.lowest_positive_eigenvalue,
                    "mode_quality_status": row_status,
                    "status": row_status,
                }
            )
    return rows


def _base_comparison_row(audit: CaseAudit) -> dict[str, Any]:
    case = audit.case
    if case.source_value is None or case.printed_tolerance is None:
        raise ValueError(
            "transfer_ritz_comparison rows require a published source value."
        )
    transfer_value = _transfer_frequency_bar(audit)
    legacy_transfer_value = _legacy_transfer_frequency_bar(audit)
    ritz_values = {
        order: _mode_frequency_bar(case, audit.modes[order])
        for order in RITZ_ORDERS
    }
    used_mode = audit.converged_mode
    ritz_used = _mode_frequency_bar(case, used_mode)
    source_vs_ritz_absolute = abs(case.source_value - ritz_used)
    source_vs_transfer_absolute = abs(case.source_value - transfer_value)
    source_vs_legacy_transfer_absolute = abs(
        case.source_value - legacy_transfer_value
    )
    original_status = (
        "PASS"
        if source_vs_legacy_transfer_absolute
        <= case.printed_tolerance + NUMERICAL_ALLOWANCE
        else "FAIL"
    )
    if original_status != audit.legacy_transfer.historical_status:
        raise RuntimeError(
            "Recomputed legacy source comparison status does not match the "
            f"immutable inventory for {case.case_id}: {original_status!r} != "
            f"{audit.legacy_transfer.historical_status!r}."
        )
    original_status = str(audit.legacy_transfer.historical_status)
    row = {
        "case_id": case.case_id,
        "laminate_id": case.laminate_id,
        "benchmark_tier": case.benchmark_tier,
        "printed_label": case.printed_label,
        "used_label": case.used_label,
        "boundary_condition": case.boundary_condition,
        "a_over_h": case.a_over_h,
        "rotary_inertia": case.rotary_inertia,
        "row_role": case.row_role,
        "K": case.K,
        "K_provenance": case.K_provenance,
        "source_value": case.source_value,
        "source_status": case.source_status,
        "printed_tolerance": case.printed_tolerance,
        "transfer_value": transfer_value,
        "transfer_value_normalization": "omega*L^2*sqrt(I0/h^3)",
        "transfer_root_omega_dimensional": audit.transfer.omega,
        "legacy_transfer_value": legacy_transfer_value,
        "legacy_transfer_root_omega_dimensional": audit.legacy_transfer.omega,
        "refined_transfer_value": transfer_value,
        "legacy_vs_refined_transfer_relative": _relative_difference(
            audit.legacy_transfer.omega, audit.transfer.omega
        ),
        "ritz_value_N4": ritz_values[4],
        "ritz_value_N6": ritz_values[6],
        "ritz_value_N8": ritz_values[8],
        "ritz_value_N10": ritz_values[10],
        "ritz_value_N12": ritz_values[12],
        "ritz_value_N16_guard": (
            ""
            if RITZ_GUARD_ORDER not in audit.modes
            else _mode_frequency_bar(case, audit.modes[RITZ_GUARD_ORDER])
        ),
        "ritz_used_order": audit.converged_order,
        "ritz_converged_value": ritz_used,
        "transfer_vs_ritz_relative": _relative_difference(
            ritz_used, transfer_value
        ),
        "ritz_convergence_relative": audit.used_convergence_relative,
        "N10_vs_N12_relative": audit.ritz_convergence_relative,
        "N16_guard_relative": (
            "" if audit.guard_relative is None else audit.guard_relative
        ),
        "source_vs_ritz_absolute": source_vs_ritz_absolute,
        "source_vs_ritz_relative": _relative_difference(ritz_used, case.source_value),
        "source_vs_transfer_absolute": source_vs_transfer_absolute,
        "source_vs_transfer_relative": _relative_difference(
            transfer_value, case.source_value
        ),
        "source_vs_legacy_transfer_absolute": source_vs_legacy_transfer_absolute,
        "source_vs_legacy_transfer_relative": _relative_difference(
            legacy_transfer_value, case.source_value
        ),
        "transfer_sigma_ratio": audit.transfer.sigma_ratio,
        "transfer_determinant": audit.transfer.determinant,
        "transfer_sigma_min": audit.transfer.sigma_min,
        "transfer_sigma_max": audit.transfer.sigma_max,
        "transfer_condition_number": audit.transfer.condition_number,
        "transfer_boundary_residual": audit.transfer.boundary_residual,
        "transfer_accepted": audit.transfer.accepted,
        "transfer_root_status": audit.transfer.detection,
        "transfer_refinement_bracket_left": audit.transfer.bracket_left,
        "transfer_refinement_bracket_right": audit.transfer.bracket_right,
        "transfer_refinement_bracket_left_omega_dimensional": (
            audit.transfer.bracket_left
        ),
        "transfer_refinement_bracket_right_omega_dimensional": (
            audit.transfer.bracket_right
        ),
        "legacy_transfer_sigma_ratio": audit.legacy_transfer.sigma_ratio,
        "legacy_transfer_accepted": audit.legacy_transfer.accepted,
        "legacy_transfer_boundary_residual": (
            audit.legacy_transfer.boundary_residual
        ),
        "legacy_transfer_root_status": audit.legacy_transfer.detection,
        "ritz_generalized_eigen_residual": used_mode.generalized_eigen_residual,
        "ritz_recovered_full_mode_residual": used_mode.recovered_full_mode_residual,
        "ritz_static_condensation_residual": used_mode.static_condensation_residual,
        "ritz_energy_identity_relative": used_mode.energy_identity_relative,
        "ritz_essential_boundary_residual": used_mode.essential_boundary_residual,
        "ritz_natural_boundary_residual": used_mode.natural_boundary_residual,
        "ritz_status": audit.ritz_status,
        "model_status": audit.model_status,
        "ritz_mode_quality_pass": _mode_quality_pass(used_mode),
        # Preserve the original Table-vs-transfer meaning of status.  New audit
        # conclusions live in separate columns below.
        "status": original_status,
        "independent_audit_status": (
            "PASS" if audit.model_status == "PASS" else audit.model_status
        ),
        "discrepancy_class": "",
    }
    return row


def build_rotary_inertia_rows(
    audits: Sequence[CaseAudit],
    transfer_inventory: dict[tuple[str, str, int, str], dict[str, str]],
) -> tuple[list[dict[str, Any]], set[tuple[str, str, int]]]:
    """Check physical RI monotonicity, including source-null project rows."""

    audits_by_key = {
        (
            audit.case.laminate_id,
            audit.case.boundary_condition,
            audit.case.a_over_h,
            audit.case.rotary_inertia,
        ): audit
        for audit in audits
    }
    pair_keys = sorted(
        {
            (
                audit.case.laminate_id,
                audit.case.boundary_condition,
                audit.case.a_over_h,
            )
            for audit in audits
        }
    )
    rows: list[dict[str, Any]] = []
    source_violations: set[tuple[str, str, int]] = set()
    for laminate_id, bc, ratio in pair_keys:
        with_key = (laminate_id, bc, ratio, "fsdt_frequency_with_RI")
        without_key = (laminate_id, bc, ratio, "fsdt_frequency_without_RI")
        if with_key not in transfer_inventory or without_key not in transfer_inventory:
            continue
        with_transfer_row = transfer_inventory[with_key]
        without_transfer_row = transfer_inventory[without_key]
        legacy_transfer_with = float(with_transfer_row["computed_value"])
        legacy_transfer_without = float(without_transfer_row["computed_value"])
        with_audit = audits_by_key.get((laminate_id, bc, ratio, True))
        without_audit = audits_by_key.get((laminate_id, bc, ratio, False))
        if with_audit is None or without_audit is None:
            raise RuntimeError(
                "The all-model audit is missing an RI pair for "
                f"{(laminate_id, bc, ratio)!r}."
            )
        transfer_with = _transfer_frequency_bar(with_audit)
        transfer_without = _transfer_frequency_bar(without_audit)
        model_transfer_monotone = transfer_with <= transfer_without * (
            1.0 + RI_MONOTONICITY_RELATIVE_TOLERANCE
        )
        ritz_with = _mode_frequency_bar(
            with_audit.case, with_audit.converged_mode
        )
        ritz_without = _mode_frequency_bar(
            without_audit.case, without_audit.converged_mode
        )
        model_ritz_monotone = ritz_with <= ritz_without * (
            1.0 + RI_MONOTONICITY_RELATIVE_TOLERANCE
        )
        source_with_text = with_transfer_row.get("source_value", "").strip()
        source_without_text = without_transfer_row.get("source_value", "").strip()
        source_pair_published = bool(source_with_text and source_without_text)
        source_with: float | str = ""
        source_without: float | str = ""
        source_monotone: bool | str = ""
        source_violation = False
        if source_pair_published:
            source_with = float(source_with_text)
            source_without = float(source_without_text)
            source_monotone = source_with <= source_without * (
                1.0 + RI_MONOTONICITY_RELATIVE_TOLERANCE
            )
            source_violation = not bool(source_monotone)
            if source_violation:
                source_violations.add((laminate_id, bc, ratio))
        model_pass = model_transfer_monotone and (
            model_ritz_monotone in (True, "")
        )
        rows.append(
            {
                "laminate_id": laminate_id,
                "boundary_condition": bc,
                "a_over_h": ratio,
                "source_pair_published": source_pair_published,
                "source_with_RI": source_with,
                "source_without_RI": source_without,
                "source_without_minus_with": (
                    ""
                    if not source_pair_published
                    else float(source_without) - float(source_with)
                ),
                "source_monotonicity_numerical_tolerance": RI_MONOTONICITY_RELATIVE_TOLERANCE,
                "source_literal_values_monotonicity": source_monotone,
                "source_violation": source_violation,
                "transfer_with_RI": transfer_with,
                "transfer_without_RI": transfer_without,
                "transfer_without_minus_with": transfer_without - transfer_with,
                "transfer_monotonicity": model_transfer_monotone,
                "legacy_transfer_with_RI": legacy_transfer_with,
                "legacy_transfer_without_RI": legacy_transfer_without,
                "legacy_transfer_without_minus_with": (
                    legacy_transfer_without - legacy_transfer_with
                ),
                "ritz_with_RI": ritz_with,
                "ritz_without_RI": ritz_without,
                "ritz_without_minus_with": ritz_without - ritz_with,
                "ritz_monotonicity": model_ritz_monotone,
                "status": "PASS" if model_pass else "FAIL",
            }
        )
    return rows, source_violations


def classify_comparison_rows(
    rows: list[dict[str, Any]],
    source_monotonicity_violations: set[tuple[str, str, int]],
) -> None:
    """Assign the approved evidence class without altering ``status``."""

    for row in rows:
        if (
            str(row["laminate_id"]),
            str(row["boundary_condition"]),
            int(row["a_over_h"]),
        ) in source_monotonicity_violations:
            classification = "SOURCE_TABLE_VIOLATES_ROTARY_INERTIA_MONOTONICITY"
        elif float(row["source_vs_ritz_absolute"]) <= (
            float(row["printed_tolerance"]) + NUMERICAL_ALLOWANCE
        ):
            classification = "MATCHES_SOURCE_WITHIN_PRINTED_TOLERANCE"
        elif (
            float(row["transfer_vs_ritz_relative"])
            <= TRANSFER_RITZ_RELATIVE_TOLERANCE
            and float(row["ritz_convergence_relative"])
            <= RITZ_CONVERGENCE_RELATIVE_TOLERANCE
            and float(row["ritz_energy_identity_relative"])
            <= ENERGY_IDENTITY_TOLERANCE
            and bool(row["ritz_mode_quality_pass"])
            and float(row["ritz_natural_boundary_residual"])
            <= NATURAL_BOUNDARY_RESIDUAL_TOLERANCE
            and row["model_status"] == "PASS"
        ):
            classification = "SOURCE_TABLE_DISAGREES_WITH_TRANSFER_AND_RITZ"
        else:
            classification = "UNRESOLVED_MODEL_DISCREPANCY"
        if classification not in DISCREPANCY_CLASSES:
            raise AssertionError(f"Unapproved discrepancy class: {classification}")
        row["discrepancy_class"] = classification


def _characteristic_record(
    contract: dict[str, Any],
    case: BenchmarkCase,
) -> dict[str, Any] | None:
    role = (
        "source_classical_characteristic_with_RI"
        if case.rotary_inertia
        else "source_classical_characteristic_without_RI"
    )
    for record in contract["records"]:
        if (
            record["laminate_id"] == case.laminate_id
            and record["boundary_condition"] == case.boundary_condition
            and int(record["a_over_h"]) == case.a_over_h
            and record["row_role"] == role
        ):
            return record
    return None


def _classical_dispersion_value(
    case: BenchmarkCase,
    characteristic_root: float,
) -> float:
    mode = rlb.bending_dispersion_branches(
        case.properties,
        float(characteristic_root) / case.length,
    )[0]
    return float(mode.omega * case.frequency_scale)


def _real_characteristic_value(value: complex | float) -> float:
    result = complex(value)
    if not (math.isfinite(result.real) and math.isfinite(result.imag)):
        return math.nan
    if abs(result.imag) > 1.0e-8 * max(abs(result.real), 1.0):
        return math.nan
    return float(result.real)


def _bracketed_root_near(
    target: float,
    function: Any,
    *,
    relative_window: float = 0.04,
) -> tuple[float, float, float]:
    lower = max(target * (1.0 - relative_window), target * 0.5)
    upper = target * (1.0 + relative_window)
    points = np.linspace(lower, upper, 1201)
    values = [_real_characteristic_value(function(float(point))) for point in points]
    brackets: list[tuple[float, float]] = []
    for left, right, left_value, right_value in zip(
        points[:-1], points[1:], values[:-1], values[1:], strict=True
    ):
        if not (math.isfinite(left_value) and math.isfinite(right_value)):
            continue
        if left_value == 0.0:
            point = float(left)
            return point, point, point
        if left_value * right_value < 0.0:
            brackets.append((float(left), float(right)))
    if not brackets:
        raise RuntimeError(f"No characteristic sign bracket near omega={target:.17g}.")
    bracket = min(
        brackets,
        key=lambda pair: abs(0.5 * (pair[0] + pair[1]) - target),
    )
    root = float(
        brentq(
            lambda omega: _real_characteristic_value(function(omega)),
            bracket[0],
            bracket[1],
            xtol=max(1.0e-15, abs(target) * 2.0e-14),
            rtol=4.0 * sys.float_info.epsilon,
        )
    )
    return root, bracket[0], bracket[1]


def reddy_cc_characteristic_local(
    omega: float,
    length: float,
    properties: rlb.BeamProperties,
) -> tuple[float, str]:
    """Evaluate Eq. (4.3.51), including its exact analytic ``J=0`` limit.

    The zero-rotary-inertia branch uses the identities
    ``S11=-S*mu**3`` and ``S22=S*lambda**3``.  No tiny surrogate value of
    ``J`` is introduced.
    """

    frequency = _positive(omega, "omega")
    L = _positive(length, "length")
    frequency_squared = frequency**2
    p = properties.D
    q = (properties.D * properties.m / properties.S + properties.J) * (
        frequency_squared
    )
    r = (
        1.0 - properties.J * frequency_squared / properties.S
    ) * properties.m * frequency_squared
    radicand = q * q + 4.0 * p * r
    if not math.isfinite(radicand) or radicand <= 0.0:
        return math.nan, "NON_REAL_CHARACTERISTIC_PARAMETERS"
    radical = math.sqrt(radicand)
    lambda_squared = (q + radical) / (2.0 * p)
    # Product identity avoids cancellation in ``(-q + radical)/(2D)``.
    mu_squared = r / (p * lambda_squared)
    if lambda_squared <= 0.0 or mu_squared <= 0.0:
        return math.nan, "NON_POSITIVE_CHARACTERISTIC_PARAMETERS"
    lambda_value = math.sqrt(lambda_squared)
    mu_value = math.sqrt(mu_squared)
    x = lambda_value * L
    y = mu_value * L
    sech_y = 0.0 if y > 700.0 else 1.0 / math.cosh(y)
    if properties.J == 0.0:
        lambda_cubed = lambda_value**3
        mu_cubed = mu_value**3
        scale = max(lambda_cubed, mu_cubed)
        value = (
            2.0
            * (lambda_cubed / scale)
            * (mu_cubed / scale)
            * (math.cos(x) - sech_y)
            + (
                (mu_cubed / scale) ** 2
                - (lambda_cubed / scale) ** 2
            )
            * math.sin(x)
            * math.tanh(y)
        )
        branch = "SCALED_ANALYTIC_J0_LIMIT_NO_TINY_J"
    else:
        S11 = mu_value * (
            properties.m * frequency_squared
            - lambda_squared * properties.S
        )
        S22 = lambda_value * (
            properties.m * frequency_squared
            + mu_squared * properties.S
        )
        if S11 == 0.0 or S22 == 0.0:
            return math.nan, "SINGULAR_FINITE_J_COEFFICIENT"
        scale = max(abs(S11), abs(S22))
        value = (
            2.0
            * (S11 / scale)
            * (S22 / scale)
            * (math.cos(x) - sech_y)
            + ((S22 / scale) ** 2 - (S11 / scale) ** 2)
            * math.sin(x)
            * math.tanh(y)
        )
        branch = "SCALED_FINITE_J_EQ_4_3_51"
    return float(value), branch


def _local_bending_state_matrix(
    omega: float,
    properties: rlb.BeamProperties,
) -> FloatArray:
    """Literal physical matrix for ``[w, psi_b, Q, M]``."""

    frequency = _positive(omega, "omega")
    frequency_squared = frequency**2
    return np.array(
        [
            [0.0, -1.0, 1.0 / properties.S, 0.0],
            [0.0, 0.0, 0.0, 1.0 / properties.D],
            [-properties.m * frequency_squared, 0.0, 0.0, 0.0],
            [0.0, -properties.J * frequency_squared, 1.0, 0.0],
        ],
        dtype=float,
    )


def _equilibrated_local_matrix(matrix: FloatArray) -> FloatArray:
    out = np.asarray(matrix, dtype=float).copy()
    row_norms = np.linalg.norm(out, axis=1)
    rows = row_norms > sys.float_info.min
    out[rows] /= row_norms[rows, None]
    column_norms = np.linalg.norm(out, axis=0)
    columns = column_norms > sys.float_info.min
    out[:, columns] /= column_norms[None, columns]
    return out


def cf_physical_boundary_matrix_local(
    omega: float,
    length: float,
    properties: rlb.BeamProperties,
) -> FloatArray:
    """Derive ``C_f exp(A L) B_c`` directly from the physical CF state."""

    L = _positive(length, "length")
    scale = np.diag(
        [L, 1.0, properties.D / L**2, properties.D / L]
    )
    inverse_scale = np.diag(1.0 / np.diag(scale))
    scaled_state_matrix = (
        inverse_scale @ _local_bending_state_matrix(omega, properties) @ scale
    )
    transfer = expm(scaled_state_matrix * L)
    B_c = np.array(
        [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        dtype=float,
    )
    C_f = np.array(
        [[0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
        dtype=float,
    )
    return np.asarray(C_f @ transfer @ B_c, dtype=float)


def cc_physical_boundary_matrix_local(
    omega: float,
    length: float,
    properties: rlb.BeamProperties,
) -> FloatArray:
    """Derive the physical CC matrix with clamps at both endpoints."""

    L = _positive(length, "length")
    scale = np.diag(
        [L, 1.0, properties.D / L**2, properties.D / L]
    )
    inverse_scale = np.diag(1.0 / np.diag(scale))
    scaled_state_matrix = (
        inverse_scale @ _local_bending_state_matrix(omega, properties) @ scale
    )
    transfer = expm(scaled_state_matrix * L)
    B_c = np.array(
        [[0.0, 0.0], [0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        dtype=float,
    )
    C_c = np.array(
        [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0]],
        dtype=float,
    )
    return np.asarray(C_c @ transfer @ B_c, dtype=float)


def hh_physical_boundary_matrix_local(
    omega: float,
    length: float,
    properties: rlb.BeamProperties,
) -> FloatArray:
    """Derive the physical HH matrix with ``w=M=0`` at both ends."""

    L = _positive(length, "length")
    scale = np.diag(
        [L, 1.0, properties.D / L**2, properties.D / L]
    )
    inverse_scale = np.diag(1.0 / np.diag(scale))
    scaled_state_matrix = (
        inverse_scale @ _local_bending_state_matrix(omega, properties) @ scale
    )
    transfer = expm(scaled_state_matrix * L)
    B_h = np.array(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.0, 0.0]],
        dtype=float,
    )
    C_h = np.array(
        [[1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 1.0]],
        dtype=float,
    )
    return np.asarray(C_h @ transfer @ B_h, dtype=float)


def physical_boundary_matrix_local(
    omega: float,
    length: float,
    properties: rlb.BeamProperties,
    boundary_condition: str,
) -> FloatArray:
    bc = _canonical_boundary_condition(boundary_condition)
    if bc == "HH":
        return hh_physical_boundary_matrix_local(omega, length, properties)
    if bc == "CC":
        return cc_physical_boundary_matrix_local(omega, length, properties)
    return cf_physical_boundary_matrix_local(omega, length, properties)


def refine_transfer_root_after_ritz(
    case: BenchmarkCase,
    legacy: TransferResult,
) -> TransferResult:
    """Refine the frozen physical transfer root without any source-value input."""

    if case.boundary_condition == "HH":
        omega = rlb.bending_dispersion_branches(
            case.properties, math.pi / case.length
        )[0].omega
        bracket_left = bracket_right = None
        detection = "exact_HH_dispersion_post_Ritz"
    else:
        determinant = lambda frequency: float(
            np.linalg.det(
                _equilibrated_local_matrix(
                    physical_boundary_matrix_local(
                        frequency,
                        case.length,
                        case.properties,
                        case.boundary_condition,
                    )
                )
            )
        )
        omega, bracket_left, bracket_right = _bracketed_root_near(
            legacy.omega,
            determinant,
            relative_window=0.01,
        )
        detection = "audit_local_physical_determinant_Brent_post_Ritz"
    matrix = _equilibrated_local_matrix(
        physical_boundary_matrix_local(
            omega, case.length, case.properties, case.boundary_condition
        )
    )
    singular = np.linalg.svd(matrix, compute_uv=False)
    ratio = float(singular[-1] / singular[0])
    return TransferResult(
        omega=float(omega),
        determinant=float(np.linalg.det(matrix)),
        sigma_min=float(singular[-1]),
        sigma_max=float(singular[0]),
        sigma_ratio=ratio,
        condition_number=(
            float(singular[0] / singular[-1])
            if singular[-1] > 0.0
            else math.inf
        ),
        # Minimum normalized free-end/clamp endpoint residual over the two
        # admissible initial-state amplitudes.
        boundary_residual=ratio,
        accepted=ratio <= TRANSFER_SINGULAR_RATIO_TOLERANCE,
        detection=detection,
        bracket_left=bracket_left,
        bracket_right=bracket_right,
    )


def cf_physical_characteristic_local(
    omega: float,
    length: float,
    properties: rlb.BeamProperties,
) -> float:
    matrix = _equilibrated_local_matrix(
        cf_physical_boundary_matrix_local(omega, length, properties)
    )
    return float(np.linalg.det(matrix))


def _first_root_sign_scan(
    function: Any,
    upper: float,
    *,
    physical_matrix: Any | None = None,
    scan_points: int = 4001,
) -> dict[str, Any]:
    """Inventory all sign roots below a known post-Ritz upper bound."""

    maximum = _positive(upper, "upper")
    minimum = max(maximum * 1.0e-6, 1.0e-14)
    points = np.linspace(minimum, maximum, scan_points)
    values = [_real_characteristic_value(function(float(point))) for point in points]
    raw_roots: list[float] = []
    for left, right, left_value, right_value in zip(
        points[:-1], points[1:], values[:-1], values[1:], strict=True
    ):
        if not (math.isfinite(left_value) and math.isfinite(right_value)):
            continue
        if left_value == 0.0:
            raw_roots.append(float(left))
        elif left_value * right_value < 0.0:
            raw_roots.append(
                float(
                    brentq(
                        lambda value: _real_characteristic_value(function(value)),
                        float(left),
                        float(right),
                        xtol=max(1.0e-15, maximum * 1.0e-14),
                        rtol=4.0 * sys.float_info.epsilon,
                    )
                )
            )
    unique: list[float] = []
    for root in sorted(raw_roots):
        if not unique or _relative_difference(root, unique[-1]) > 1.0e-8:
            unique.append(root)
    accepted: list[float] = []
    rejected: list[float] = []
    minimum_sigma_ratios: list[float] = []
    if physical_matrix is None:
        accepted = unique.copy()
    else:
        for root in unique:
            matrix = _equilibrated_local_matrix(physical_matrix(root))
            singular = np.linalg.svd(matrix, compute_uv=False)
            ratio = float(singular[-1] / singular[0])
            minimum_sigma_ratios.append(ratio)
            if ratio <= TRANSFER_SINGULAR_RATIO_TOLERANCE:
                accepted.append(root)
            else:
                rejected.append(root)
    return {
        "scan_lower": minimum,
        "scan_upper": maximum,
        "scan_points": scan_points,
        "raw_sign_root_count": len(unique),
        "accepted_physical_root_count": len(accepted),
        "rejected_spurious_root_count": len(rejected),
        "first_accepted_root": "" if not accepted else accepted[0],
        "accepted_roots": accepted,
        "rejected_roots": rejected,
        "candidate_sigma_ratios": minimum_sigma_ratios,
    }


def _classical_source_fields(
    contract: dict[str, Any],
    case: BenchmarkCase,
) -> dict[str, Any]:
    record = _characteristic_record(contract, case)
    printed_root = float(
        contract["source_classical_characteristic_roots"][
            case.boundary_condition
        ]["value"]
    )
    exact_root = FULL_PRECISION_CLASSICAL_ROOTS[case.boundary_condition]
    printed_value = _classical_dispersion_value(case, printed_root)
    exact_value = _classical_dispersion_value(case, exact_root)
    source_value = None if record is None else record["source_value"]
    tolerance = None if record is None else record["printed_tolerance"]
    source_status = "NOT_PUBLISHED_IN_TABLE"
    source_absolute_error: float | str = ""
    source_exact_absolute_error: float | str = ""
    discrepancy_class = ""
    if source_value is not None and tolerance is not None:
        source_absolute_error = abs(float(source_value) - printed_value)
        source_exact_absolute_error = abs(float(source_value) - exact_value)
        source_status = (
            "PASS"
            if source_absolute_error
            <= float(tolerance) + NUMERICAL_ALLOWANCE
            else "FAIL"
        )
        if source_status == "PASS":
            discrepancy_class = "MATCHES_SOURCE_WITHIN_PRINTED_TOLERANCE"
        else:
            # Qualified source-procedure label only: the source row is based
            # on a printed approximate characteristic input and may include
            # hidden intermediate rounding.  It is not a claim that replacing
            # 4.730/1.875 by a full-precision root corrects the printed value.
            discrepancy_class = "PRINTED_CHARACTERISTIC_DIAGNOSTIC_ROUNDING"
    return {
        "source_classical_row_role": (
            "" if record is None else record["row_role"]
        ),
        "source_classical_value": (
            "" if source_value is None else float(source_value)
        ),
        "source_classical_printed_tolerance": (
            "" if tolerance is None else float(tolerance)
        ),
        "printed_characteristic_root": printed_root,
        "full_precision_characteristic_root": exact_root,
        "printed_root_dispersion_value": printed_value,
        "full_precision_root_dispersion_value": exact_value,
        "printed_vs_full_precision_relative": _relative_difference(
            printed_value, exact_value
        ),
        "source_vs_printed_root_absolute": source_absolute_error,
        "source_vs_full_precision_root_absolute": source_exact_absolute_error,
        "source_classical_status": source_status,
        "source_classical_discrepancy_class": discrepancy_class,
        "source_classical_classification_qualification": (
            ""
            if source_status != "FAIL"
            else (
                "Qualified source-procedure class; both printed-root and "
                "full-precision-root errors are retained, and no corrected "
                "source reproduction is claimed."
            )
        ),
    }


def build_cc_characteristic_rows(
    audits: Sequence[CaseAudit],
    contract: dict[str, Any],
) -> list[dict[str, Any]]:
    """Compare Eq. (4.3.51), transfer, Ritz, and printed CC diagnostics."""

    rows: list[dict[str, Any]] = []
    for audit in audits:
        if audit.case.boundary_condition != "CC":
            continue
        case = audit.case
        ritz_mode = audit.converged_mode
        error = ""
        eq_root: float | str = ""
        bracket_left: float | str = ""
        bracket_right: float | str = ""
        eq_relative: float | str = ""
        characteristic_at_transfer: float | str = ""
        coefficient_branch = ""
        first_scan: dict[str, Any] = {
            "scan_lower": "",
            "scan_upper": "",
            "scan_points": 0,
            "raw_sign_root_count": 0,
            "accepted_physical_root_count": 0,
            "rejected_spurious_root_count": 0,
            "first_accepted_root": "",
        }
        try:
            characteristic = lambda omega: reddy_cc_characteristic_local(
                omega, case.length, case.properties
            )[0]
            eq_root_value, bracket_left, bracket_right = _bracketed_root_near(
                ritz_mode.omega,
                characteristic,
            )
            eq_root = eq_root_value
            _value, coefficient_branch = reddy_cc_characteristic_local(
                eq_root_value, case.length, case.properties
            )
            eq_relative = _relative_difference(eq_root_value, audit.transfer.omega)
            characteristic_at_transfer = _real_characteristic_value(
                characteristic(audit.transfer.omega)
            )
            first_scan = _first_root_sign_scan(
                characteristic,
                ritz_mode.omega * 1.05,
                physical_matrix=lambda omega: cc_physical_boundary_matrix_local(
                    omega, case.length, case.properties
                ),
            )
            first_root = first_scan["first_accepted_root"]
            no_earlier_physical_root = (
                first_root != ""
                and _relative_difference(float(first_root), eq_root_value)
                <= TRANSFER_RITZ_RELATIVE_TOLERANCE
            )
            eq_pass = (
                eq_relative <= TRANSFER_RITZ_RELATIVE_TOLERANCE
                and no_earlier_physical_root
                and first_scan["rejected_spurious_root_count"] == 0
            )
        except (RuntimeError, ValueError, FloatingPointError) as exception:
            eq_pass = False
            error = f"{type(exception).__name__}: {exception}"
        ritz_matrix = _equilibrated_local_matrix(
            cc_physical_boundary_matrix_local(
                ritz_mode.omega, case.length, case.properties
            )
        )
        ritz_singular = np.linalg.svd(ritz_matrix, compute_uv=False)
        ritz_sigma_ratio = float(ritz_singular[-1] / ritz_singular[0])
        rows.append(
            {
                "case_id": case.case_id,
                "laminate_id": case.laminate_id,
                "a_over_h": case.a_over_h,
                "rotary_inertia": case.rotary_inertia,
                "physical_boundary_conditions": (
                    "w=psi_b=0 at x=0 and x=length"
                ),
                "transfer_value": _transfer_frequency_bar(audit),
                "ritz_value_N12": _mode_frequency_bar(case, audit.modes[12]),
                "ritz_used_order": audit.converged_order,
                "ritz_converged_value": _mode_frequency_bar(case, ritz_mode),
                "eq_4_3_51_value": (
                    "" if eq_root == "" else float(eq_root) * case.frequency_scale
                ),
                "eq_4_3_51_vs_transfer_relative": eq_relative,
                "eq_4_3_51_characteristic_at_transfer": characteristic_at_transfer,
                "eq_4_3_51_coefficient_branch": coefficient_branch,
                "J_zero_exact_limit_used": case.properties.J == 0.0,
                "tiny_J_surrogate_used": False,
                "eq_4_3_51_bracket_left": bracket_left,
                "eq_4_3_51_bracket_right": bracket_right,
                "transfer_vs_ritz_relative": _relative_difference(
                    ritz_mode.omega, audit.transfer.omega
                ),
                "physical_CC_determinant_at_ritz": float(np.linalg.det(ritz_matrix)),
                "physical_CC_sigma_min_at_ritz": float(ritz_singular[-1]),
                "physical_CC_sigma_max_at_ritz": float(ritz_singular[0]),
                "physical_CC_sigma_ratio_at_ritz": ritz_sigma_ratio,
                "ritz_natural_boundary_residual": ritz_mode.natural_boundary_residual,
                "first_root_scan_lower": first_scan["scan_lower"],
                "first_root_scan_upper": first_scan["scan_upper"],
                "first_root_scan_points": first_scan["scan_points"],
                "first_root_scan_raw_count": first_scan["raw_sign_root_count"],
                "first_root_scan_accepted_physical_count": first_scan[
                    "accepted_physical_root_count"
                ],
                "first_root_scan_rejected_spurious_count": first_scan[
                    "rejected_spurious_root_count"
                ],
                "first_root_scan_first_accepted": first_scan[
                    "first_accepted_root"
                ],
                "first_root_scan_accepted_roots_omega": json.dumps(
                    first_scan.get("accepted_roots", []), separators=(",", ":")
                ),
                "first_root_scan_rejected_roots_omega": json.dumps(
                    first_scan.get("rejected_roots", []), separators=(",", ":")
                ),
                "first_root_scan_candidate_sigma_ratios": json.dumps(
                    first_scan.get("candidate_sigma_ratios", []),
                    separators=(",", ":"),
                ),
                **_classical_source_fields(contract, case),
                "error": error,
                "status": (
                    "PASS"
                    if (
                        eq_pass
                        and ritz_sigma_ratio
                        <= TRANSFER_SINGULAR_RATIO_TOLERANCE
                        and audit.model_status == "PASS"
                    )
                    else "FAIL"
                ),
            }
        )
    return rows


def build_cf_characteristic_rows(
    audits: Sequence[CaseAudit],
    contract: dict[str, Any],
) -> list[dict[str, Any]]:
    """Expose the physical CF endpoint check and source classical diagnostic.

    Reddy Eq. (4.3.51) is a CC formula and is not repurposed for CF.  The CF
    hook therefore evaluates the independent transfer boundary singularity at
    the Ritz value and records both the printed and full-precision classical
    diagnostic roots without using either as the physical solver source.
    """

    rows: list[dict[str, Any]] = []
    for audit in audits:
        if audit.case.boundary_condition != "CF":
            continue
        case = audit.case
        ritz_mode = audit.converged_mode
        characteristic = lambda omega: cf_physical_characteristic_local(
            omega, case.length, case.properties
        )
        first_scan = _first_root_sign_scan(
            characteristic,
            ritz_mode.omega * 1.05,
            physical_matrix=lambda omega: cf_physical_boundary_matrix_local(
                omega, case.length, case.properties
            ),
        )
        first_root = first_scan["first_accepted_root"]
        first_root_relative: float | str = ""
        if first_root != "":
            first_root_relative = _relative_difference(
                float(first_root), audit.transfer.omega
            )
        ritz_matrix = _equilibrated_local_matrix(
            cf_physical_boundary_matrix_local(
                ritz_mode.omega, case.length, case.properties
            )
        )
        ritz_singular = np.linalg.svd(ritz_matrix, compute_uv=False)
        ritz_sigma_ratio = float(ritz_singular[-1] / ritz_singular[0])
        rows.append(
            {
                "case_id": case.case_id,
                "laminate_id": case.laminate_id,
                "a_over_h": case.a_over_h,
                "rotary_inertia": case.rotary_inertia,
                "physical_left_boundary_conditions": "w=psi_b=0 at x=0",
                "physical_right_boundary_conditions": "Q=M=0 at x=length",
                "eq_4_3_51_applicable": False,
                "transfer_value": _transfer_frequency_bar(audit),
                "ritz_value_N12": _mode_frequency_bar(case, audit.modes[12]),
                "ritz_used_order": audit.converged_order,
                "ritz_converged_value": _mode_frequency_bar(case, ritz_mode),
                "local_physical_first_root_value": (
                    ""
                    if first_root == ""
                    else float(first_root) * case.frequency_scale
                ),
                "local_physical_first_root_vs_transfer_relative": first_root_relative,
                "transfer_vs_ritz_relative": _relative_difference(
                    ritz_mode.omega, audit.transfer.omega
                ),
                "physical_CF_determinant_at_ritz": float(np.linalg.det(ritz_matrix)),
                "physical_CF_sigma_min_at_ritz": float(ritz_singular[-1]),
                "physical_CF_sigma_max_at_ritz": float(ritz_singular[0]),
                "physical_CF_sigma_ratio_at_ritz": ritz_sigma_ratio,
                "ritz_natural_boundary_residual": ritz_mode.natural_boundary_residual,
                "first_root_scan_lower": first_scan["scan_lower"],
                "first_root_scan_upper": first_scan["scan_upper"],
                "first_root_scan_points": first_scan["scan_points"],
                "first_root_scan_raw_count": first_scan["raw_sign_root_count"],
                "first_root_scan_accepted_physical_count": first_scan[
                    "accepted_physical_root_count"
                ],
                "first_root_scan_rejected_spurious_count": first_scan[
                    "rejected_spurious_root_count"
                ],
                "first_root_scan_first_accepted": first_scan[
                    "first_accepted_root"
                ],
                "first_root_scan_accepted_roots_omega": json.dumps(
                    first_scan.get("accepted_roots", []), separators=(",", ":")
                ),
                "first_root_scan_rejected_roots_omega": json.dumps(
                    first_scan.get("rejected_roots", []), separators=(",", ":")
                ),
                "first_root_scan_candidate_sigma_ratios": json.dumps(
                    first_scan.get("candidate_sigma_ratios", []),
                    separators=(",", ":"),
                ),
                **_classical_source_fields(contract, case),
                "status": (
                    "PASS"
                    if (
                        first_root_relative != ""
                        and float(first_root_relative)
                        <= TRANSFER_RITZ_RELATIVE_TOLERANCE
                        and first_scan["rejected_spurious_root_count"] == 0
                        and ritz_sigma_ratio
                        <= TRANSFER_SINGULAR_RATIO_TOLERANCE
                        and audit.model_status == "PASS"
                    )
                    else "FAIL"
                ),
            }
        )
    return rows


def build_source_inconsistency_inventory(
    comparison_rows: Sequence[dict[str, Any]],
    cc_rows: Sequence[dict[str, Any]],
    cf_rows: Sequence[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Collect direct and source-classical inconsistencies without relabeling."""

    rows: list[dict[str, Any]] = []
    for comparison in comparison_rows:
        if comparison["discrepancy_class"] == (
            "MATCHES_SOURCE_WITHIN_PRINTED_TOLERANCE"
        ):
            continue
        rows.append(
            {
                "inventory_kind": "published_direct_FSDT",
                "case_id": comparison["case_id"],
                "laminate_id": comparison["laminate_id"],
                "boundary_condition": comparison["boundary_condition"],
                "a_over_h": comparison["a_over_h"],
                "rotary_inertia": comparison["rotary_inertia"],
                "row_role": comparison["row_role"],
                "source_value": comparison["source_value"],
                "comparison_value": comparison["ritz_converged_value"],
                "secondary_value": comparison["transfer_value"],
                "source_vs_comparison_absolute": comparison[
                    "source_vs_ritz_absolute"
                ],
                "source_vs_comparison_relative": comparison[
                    "source_vs_ritz_relative"
                ],
                "printed_tolerance": comparison["printed_tolerance"],
                "original_status": comparison["status"],
                "independent_audit_status": comparison[
                    "independent_audit_status"
                ],
                "discrepancy_class": comparison["discrepancy_class"],
                "qualification": "",
            }
        )
    for characteristic in [*cc_rows, *cf_rows]:
        if characteristic["source_classical_status"] != "FAIL":
            continue
        rows.append(
            {
                "inventory_kind": "source_classical_characteristic_diagnostic",
                "case_id": characteristic["case_id"],
                "laminate_id": characteristic["laminate_id"],
                "boundary_condition": (
                    "CC" if "eq_4_3_51_value" in characteristic else "CF"
                ),
                "a_over_h": characteristic["a_over_h"],
                "rotary_inertia": characteristic["rotary_inertia"],
                "row_role": characteristic["source_classical_row_role"],
                "source_value": characteristic["source_classical_value"],
                "comparison_value": characteristic[
                    "printed_root_dispersion_value"
                ],
                "secondary_value": characteristic[
                    "full_precision_root_dispersion_value"
                ],
                "source_vs_comparison_absolute": characteristic[
                    "source_vs_printed_root_absolute"
                ],
                "source_vs_comparison_relative": _relative_difference(
                    float(characteristic["printed_root_dispersion_value"]),
                    float(characteristic["source_classical_value"]),
                ),
                "printed_tolerance": characteristic[
                    "source_classical_printed_tolerance"
                ],
                "original_status": characteristic["source_classical_status"],
                "independent_audit_status": characteristic["status"],
                "discrepancy_class": characteristic[
                    "source_classical_discrepancy_class"
                ],
                "qualification": characteristic[
                    "source_classical_classification_qualification"
                ],
                "source_vs_full_precision_root_absolute": characteristic[
                    "source_vs_full_precision_root_absolute"
                ],
            }
        )
    return rows


def _aggregate_status(values: Sequence[str]) -> str:
    if values and all(value == "PASS" for value in values):
        return "PASS"
    if any(value == "FAIL" for value in values):
        return "FAIL"
    return "PARTIAL_PASS"


def _finite_max(values: Iterable[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    return max(finite, default=0.0)


def _finite_min(values: Iterable[float]) -> float:
    finite = [float(value) for value in values if math.isfinite(float(value))]
    return min(finite, default=0.0)


def build_summary(
    audits: Sequence[CaseAudit],
    comparison_rows: Sequence[dict[str, Any]],
    rotary_rows: Sequence[dict[str, Any]],
    inconsistency_rows: Sequence[dict[str, Any]],
    cc_rows: Sequence[dict[str, Any]],
    cf_rows: Sequence[dict[str, Any]],
    *,
    smoke: bool,
) -> dict[str, Any]:
    converged_modes = [audit.converged_mode for audit in audits]
    ritz_status = _aggregate_status([audit.ritz_status for audit in audits])
    model_status = _aggregate_status([audit.model_status for audit in audits])
    ri_status = (
        "PASS" if rotary_rows and all(row["status"] == "PASS" for row in rotary_rows)
        else "FAIL"
    )
    unresolved = sum(
        row["discrepancy_class"] == "UNRESOLVED_MODEL_DISCREPANCY"
        for row in inconsistency_rows
    )
    source_issue_count = sum(
        row["discrepancy_class"]
        != "MATCHES_SOURCE_WITHIN_PRINTED_TOLERANCE"
        for row in comparison_rows
    ) + sum(
        row["source_classical_status"] == "FAIL" for row in [*cc_rows, *cf_rows]
    )
    if unresolved:
        table_status = "FAIL"
    elif source_issue_count:
        table_status = "PARTIAL_PASS_SOURCE_TABLE_INCONSISTENCY"
    else:
        table_status = "PASS"
    model_gates_pass = (
        model_status == "PASS"
        and ritz_status == "PASS"
        and ri_status == "PASS"
        and all(row["status"] == "PASS" for row in cc_rows)
        and all(row["status"] == "PASS" for row in cf_rows)
        and unresolved == 0
    )
    if smoke:
        overall = "PARTIAL_PASS" if model_gates_pass else "FAIL"
    elif model_gates_pass:
        overall = "PASS_WITH_SOURCE_TABLE_LIMITATIONS"
    elif model_status == "FAIL" or ritz_status == "FAIL" or ri_status == "FAIL":
        overall = "FAIL"
    else:
        overall = "PARTIAL_PASS"

    discrepancy_counts = {
        classification: sum(
            row["discrepancy_class"] == classification
            for row in comparison_rows
        )
        + sum(
            row["source_classical_status"] == "FAIL"
            and row["source_classical_discrepancy_class"] == classification
            for row in [*cc_rows, *cf_rows]
        )
        for classification in sorted(DISCREPANCY_CLASSES)
    }
    matrix_gates = {
        "maximum_stiffness_symmetry_residual": _finite_max(
            mode.stiffness_symmetry_residual for mode in converged_modes
        ),
        "maximum_mass_symmetry_residual": _finite_max(
            mode.mass_symmetry_residual for mode in converged_modes
        ),
        "minimum_stiffness_eigenvalue": _finite_min(
            mode.stiffness_minimum_eigenvalue for mode in converged_modes
        ),
        "all_stiffness_SPD": all(mode.stiffness_spd for mode in converged_modes),
        "minimum_full_mass_eigenvalue_J_positive": _finite_min(
            mode.full_mass_minimum_eigenvalue
            for mode in converged_modes
            if mode.full_mass_minimum_eigenvalue is not None
        ),
        "all_full_mass_SPD_J_positive": all(
            mode.full_mass_spd is True
            for mode in converged_modes
            if not mode.static_condensation_used
        ),
        "minimum_translational_mass_eigenvalue": _finite_min(
            mode.translational_mass_minimum_eigenvalue for mode in converged_modes
        ),
        "all_translational_mass_SPD": all(
            mode.translational_mass_spd for mode in converged_modes
        ),
        "maximum_J0_psi_mass_zero_residual": _finite_max(
            mode.psi_mass_zero_residual
            for mode in converged_modes
            if mode.static_condensation_used
        ),
        "all_J0_mass_structures_exact": all(
            mode.j0_mass_structure_pass
            for mode in converged_modes
            if mode.static_condensation_used
        ),
        "maximum_zero_mode_count": max(
            (mode.zero_mode_count for mode in converged_modes), default=0
        ),
        "minimum_lowest_positive_eigenvalue": _finite_min(
            mode.lowest_positive_eigenvalue for mode in converged_modes
        ),
    }
    numerical = {
        "all_model_case_count": len(audits),
        "published_direct_case_count": len(comparison_rows),
        "rotary_pair_count": len(rotary_rows),
        "CC_characteristic_case_count": len(cc_rows),
        "CF_characteristic_case_count": len(cf_rows),
        "maximum_transfer_vs_Ritz_relative": _finite_max(
            row["transfer_vs_ritz_relative"] for row in comparison_rows
        ),
        "maximum_N10_vs_N12_relative": _finite_max(
            audit.ritz_convergence_relative for audit in audits
        ),
        "maximum_used_convergence_relative": _finite_max(
            audit.used_convergence_relative for audit in audits
        ),
        "N16_guard_case_count": sum(
            RITZ_GUARD_ORDER in audit.modes for audit in audits
        ),
        "N16_used_case_count": sum(
            audit.converged_order == RITZ_GUARD_ORDER for audit in audits
        ),
        "maximum_legacy_vs_refined_transfer_relative": _finite_max(
            _relative_difference(
                audit.legacy_transfer.omega, audit.transfer.omega
            )
            for audit in audits
        ),
        "maximum_solver_vs_Rayleigh_relative": _finite_max(
            mode.rayleigh_eigenvalue_relative_difference for mode in converged_modes
        ),
        "maximum_generalized_eigen_residual": _finite_max(
            mode.generalized_eigen_residual for mode in converged_modes
        ),
        "maximum_matrix_norm_backward_residual": _finite_max(
            mode.matrix_norm_backward_residual for mode in converged_modes
        ),
        "maximum_recovered_full_mode_residual": _finite_max(
            mode.recovered_full_mode_residual for mode in converged_modes
        ),
        "maximum_recovered_full_mode_residual_J0": _finite_max(
            mode.recovered_full_mode_residual
            for mode in converged_modes
            if mode.static_condensation_used
        ),
        "maximum_raw_coordinate_residual_J_positive": _finite_max(
            mode.generalized_eigen_residual
            for mode in converged_modes
            if not mode.static_condensation_used
        ),
        "maximum_energy_identity_relative": _finite_max(
            mode.energy_identity_relative for mode in converged_modes
        ),
        "maximum_natural_boundary_residual": _finite_max(
            mode.natural_boundary_residual for mode in converged_modes
        ),
        "maximum_static_condensation_residual": _finite_max(
            mode.static_condensation_residual
            for mode in converged_modes
            if mode.static_condensation_residual is not None
        ),
        "source_literal_RI_violation_pair_count": sum(
            bool(row["source_violation"]) for row in rotary_rows
        ),
        "source_literal_RI_violation_member_count": sum(
            row["discrepancy_class"]
            == "SOURCE_TABLE_VIOLATES_ROTARY_INERTIA_MONOTONICITY"
            for row in comparison_rows
        ),
        "original_direct_source_fail_count": sum(
            row["status"] == "FAIL" for row in comparison_rows
        ),
        "independent_unresolved_count": unresolved,
    }
    tier_counts: dict[str, dict[str, Any]] = {}
    for row in comparison_rows:
        tier = str(row["benchmark_tier"])
        entry = tier_counts.setdefault(
            tier,
            {
                "published_direct_count": 0,
                "historical_PASS": 0,
                "historical_FAIL": 0,
                "discrepancy_classes": {
                    classification: 0
                    for classification in sorted(DISCREPANCY_CLASSES)
                },
            },
        )
        entry["published_direct_count"] += 1
        entry[f"historical_{row['status']}"] += 1
        entry["discrepancy_classes"][row["discrepancy_class"]] += 1
    return {
        "statuses": {
            STATUS_MODEL: model_status,
            STATUS_RITZ: ritz_status,
            STATUS_RI: ri_status,
            STATUS_TABLE: table_status,
            STATUS_OVERALL: overall,
        },
        "matrix_gates": matrix_gates,
        "numerical": numerical,
        "discrepancy_counts": discrepancy_counts,
        "tier_counts": tier_counts,
        "smoke": smoke,
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        result = float(value)
        return result if math.isfinite(result) else None
    if isinstance(value, Path):
        return str(value)
    return value


def _write_csv(path: Path, rows: Sequence[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty required CSV: {path.name}.")
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        ""
                        if row.get(key) is None
                        else row.get(key, "")
                    )
                    for key in fields
                }
            )


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            _json_safe(payload),
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
        + "\n",
        encoding="utf-8",
    )


def _report_text(
    summary: dict[str, Any],
    transfer_metadata: dict[str, Any],
) -> str:
    statuses = summary["statuses"]
    numerical = summary["numerical"]
    matrices = summary["matrix_gates"]
    tier_lines = "\n".join(
        (
            f"- {tier}: {values['historical_PASS']} historical PASS, "
            f"{values['historical_FAIL']} historical FAIL, "
            f"{values['published_direct_count']} direct rows"
        )
        for tier, values in sorted(summary["tier_counts"].items())
    )
    return f"""# RLB-0A-S1: независимая проверка Table 4.3.3 Reddy

## Результат

```text
{STATUS_MODEL}: {statuses[STATUS_MODEL]}
{STATUS_RITZ}: {statuses[STATUS_RITZ]}
{STATUS_RI}: {statuses[STATUS_RI]}
{STATUS_TABLE}: {statuses[STATUS_TABLE]}
{STATUS_OVERALL}: {statuses[STATUS_OVERALL]}
```

Проверка ограничена одним прямым симметрично слоистым стержнем. Физическая
модель RLB-0 и опубликованные значения источника не изменяются и не
подгоняются на этом этапе. Все выводы условны относительно ранее принятого
`K=5/6` с provenance `INFERRED_FROM_INTERNAL_NUMERICAL_CONSISTENCY`.

## Независимая постановка Рэлея—Ритца

При `xi=x/length` и `W=w/length` собраны формы

```text
K = integral[psi_,xi^2 + gamma*(W_,xi+psi)^2] dxi
M = integral[W^2 + rho_r*psi^2] dxi
gamma = S*length^2/D
rho_r = J/(m*length^2)
lambda = omega^2*m*length^4/D.
```

Shifted-Legendre polynomials умножаются на `xi*(1-xi)` или `xi`, чтобы
наложить только существенные условия HH, CC и CF. Естественные условия
следуют из вариационной постановки. Для порядков `{RITZ_ORDERS}` используется
фиксированная {GAUSS_LEGENDRE_ORDER}-точечная квадратура Гаусса—Лежандра. При
`J=0` вращательный блок статически конденсируется точно, после чего
восстанавливается полный двухполевой вектор. Искусственная малая вращательная
масса не вводится. Обобщённая собственная задача решается посредством явного
mass-orthonormal congruence и обратного преобразования. Собственное значение
решателя и Rayleigh quotient сохраняются раздельно.

## Независимое сравнение

Неизменяемый transfer inventory `{transfer_metadata['path']}` загружен и
проверен по хешу до цикла. Его численные transfer values не используются до
того, как для соответствующего case Ritz convergence/quality независимо
определили выбранный порядок и eigenpair. SHA-256 inventory равен
`{transfer_metadata['sha256']}`. Все {transfer_metadata['direct_row_count']}
direct source keys прошли сверку с замороженным JSON.

- все модельные direct cases: {numerical['all_model_case_count']}
- опубликованные direct rows: {numerical['published_direct_case_count']}
- максимум transfer/Ritz: {numerical['maximum_transfer_vs_Ritz_relative']:.6e}
- максимум immutable-legacy/refined-transfer: {numerical['maximum_legacy_vs_refined_transfer_relative']:.6e}
- максимум N10/N12: {numerical['maximum_N10_vs_N12_relative']:.6e}
- максимум выбранной сходимости N12/N16: {numerical['maximum_used_convergence_relative']:.6e}
- число рассчитанных N16 guards: {numerical['N16_guard_case_count']}
- число случаев с выбранным N16: {numerical['N16_used_case_count']}
- максимум residual естественных условий: {numerical['maximum_natural_boundary_residual']:.6e}
- максимум generalized-eigen residual: {numerical['maximum_generalized_eigen_residual']:.6e}
- максимум полного residual при точной конденсации `J=0`: {numerical['maximum_recovered_full_mode_residual_J0']:.6e}
- максимум raw-coordinate residual при `J>0`: {numerical['maximum_raw_coordinate_residual_J_positive']:.6e}
- максимум solver/Rayleigh: {numerical['maximum_solver_vs_Rayleigh_relative']:.6e}
- максимум energy-identity residual: {numerical['maximum_energy_identity_relative']:.6e}

Жёсткий порог transfer/Ritz равен
`{TRANSFER_RITZ_RELATIVE_TOLERANCE:.1e}`. Порядок N12 не заменяется скрыто:
N16 используется только как явно маркированный guard.

Локальный energy tolerance `{ENERGY_IDENTITY_TOLERANCE:.1e}` является порогом
численной обусловленности независимой Ritz-проверки и не изменяет printed
tolerance источника. Algebraic gate использует matrix-norm backward residual
для всех мод и дополнительно полный восстановленный residual для точной
конденсации при `J=0`. Raw-coordinate residual при `J>0` сохраняется отдельно
из-за обусловленности координат матрицы массы.

## Матричные gates

- максимум stiffness symmetry residual: {matrices['maximum_stiffness_symmetry_residual']:.6e}
- все матрицы жёсткости после essential BC положительно определены: {matrices['all_stiffness_SPD']}
- все полные матрицы массы при `J>0` положительно определены: {matrices['all_full_mass_SPD_J_positive']}
- все translational mass blocks положительно определены: {matrices['all_translational_mass_SPD']}
- максимум residual нулевого psi-mass block при `J=0`: {matrices['maximum_J0_psi_mass_zero_residual']:.6e}
- максимальное число нулевых мод: {matrices['maximum_zero_mode_count']}

## Вращательная инерция и классификация источника

С вращательной инерцией и без неё проверены все
{numerical['rotary_pair_count']} комбинаций ламината, граничных условий и
гибкости, включая 18 source-null multilayer diagnostics. Transfer и Ritz
spectra удовлетворяют физическому условию монотонности. Буквальные значения
источника нарушают его в
{numerical['source_literal_RI_violation_pair_count']} парах, то есть для
{numerical['source_literal_RI_violation_member_count']} опубликованных строк.
Этот класс назначается обоим членам пары даже при историческом `status=PASS`;
исходный статус сохраняется.

Раздельные исторические итоги published direct rows:

{tier_lines}

Direct row вне printed tolerance получает класс несогласованности таблицы
только после прохождения transfer/Ritz, выбранного N12/N16 convergence guard,
algebraic, energy и endpoint gates. Класс округления без независимого
основания не назначается.

## Диагностика CC и CF

Для CC локально проверяется Eq. (4.3.51). Значение
`mu^2=r/(D*lambda^2)` вычисляется из произведения корней. При конечном `J`
используется устойчивая форма с `sech` и `tanh`, масштабированная на
`max(|S11|,|S22|)^2`. При `J=0` точно применяются
`S11=-S*mu^3`, `S22=S*lambda^3` и соответствующая масштабированная форма
`G0`; малое `J` не вводится. Для CF независимо строится
`det(C_f exp(A*length) B_c)` из physical state `[w,psi_b,Q,M]` и условий
`w=psi_b=0` слева, `Q=M=0` справа. Оба CSV содержат inventory sign roots ниже
physical fundamental с пятипроцентным запасом и отвергнутые spurious
candidates, включая их root locations и singular ratios.

Класс `PRINTED_CHARACTERISTIC_DIAGNOSTIC_ROUNDING` является только
квалифицированным описанием source procedure. Сохранены ошибки для
напечатанного и полноточного корней; замена напечатанного корня не объявляется
исправлением source value.

## Problem 4.14 Reddy

Problem 4.14 имеет статус `{PROBLEM_4_14_DIAGNOSTIC['status']}` (напечатанная
с. {PROBLEM_4_14_DIAGNOSTIC['printed_page']}, zero-based PDF index
{PROBLEM_4_14_DIAGNOSTIC['pdf_index_zero_based']}). Символ `beta` в указанном
контексте не определён. Он не подбирается и не вычисляется, поэтому эта
вторичная диагностика не является исполняемым уравнением.

## Ограничения

- Совпадение Рэлея—Ритца и transfer roots проверяет реализованную одномерную
  модель, но не доказывает правильность каждой строки печатной таблицы.
- Проверены только fundamental direct-FSDT rows Table 4.3.3.
- Classical-characteristic rows остаются вторичными diagnostics и не служат
  источником roots.
- Coupled rods, angular joint, torsion, damping, FEM и несимметричные
  ламинаты не вводятся.
- Исторические HH/RI controls участвовали в прежнем обосновании принятого
  `K=5/6`; их совпадение не рассматривается как независимая повторная
  калибровка или validation значения `K` в RLB-0A-S1.
"""


def _validate_output_directory(output_dir: Path, *, smoke: bool) -> None:
    resolved = output_dir.resolve()
    protected = LEGACY_TRANSFER_INVENTORY_PATH.parent.resolve()
    if not smoke and resolved != DEFAULT_OUTPUT_DIR.resolve():
        raise ValueError(
            "A full RLB-0A-S1 run may write only to the required "
            f"audit directory: {DEFAULT_OUTPUT_DIR}."
        )
    if resolved == protected:
        raise ValueError(
            "The new audit output directory must not equal the existing RLB-0 "
            "result directory."
        )
    if output_dir.exists():
        extra = [
            item.name
            for item in output_dir.iterdir()
            if item.name not in REQUIRED_OUTPUTS
        ]
        if extra:
            raise RuntimeError(
                "The audit output directory contains non-contract entries; "
                f"refusing to alter them: {sorted(extra)!r}."
            )


def run_audit(
    output_dir: Path = DEFAULT_OUTPUT_DIR,
    *,
    smoke: bool = False,
    include_n16_guard: bool = False,
) -> dict[str, Any]:
    """Run the independent audit and write exactly the eight contract files."""

    _validate_output_directory(output_dir, smoke=smoke)
    protected_outputs_before = _require_protected_legacy_outputs()
    source_pdf_sha_before = _require_sha256(
        SOURCE_PDF_PATH,
        SOURCE_PDF_SHA256,
        "Reddy source PDF",
    )
    contract, all_cases = load_all_direct_cases()
    transfer_inventory, transfer_metadata = load_legacy_transfer_inventory(
        contract
    )
    transfer_sha_before = transfer_metadata["sha256"]
    if smoke:
        cases = [
            case
            for case in all_cases
            if case.a_over_h == 20
        ]
    else:
        cases = all_cases
    audits = audit_direct_cases(
        cases,
        transfer_inventory,
        include_guard=include_n16_guard,
    )
    published_audits = [
        audit for audit in audits if audit.case.source_value is not None
    ]
    convergence_rows = build_ritz_convergence_rows(audits)
    comparison_rows = [_base_comparison_row(audit) for audit in published_audits]
    rotary_rows, source_monotonicity_violations = build_rotary_inertia_rows(
        audits, transfer_inventory
    )
    classify_comparison_rows(comparison_rows, source_monotonicity_violations)
    cc_rows = build_cc_characteristic_rows(audits, contract)
    cf_rows = build_cf_characteristic_rows(audits, contract)
    inconsistency_rows = build_source_inconsistency_inventory(
        comparison_rows, cc_rows, cf_rows
    )
    summary = build_summary(
        audits,
        comparison_rows,
        rotary_rows,
        inconsistency_rows,
        cc_rows,
        cf_rows,
        smoke=smoke,
    )
    transfer_sha_after = _sha256(LEGACY_TRANSFER_INVENTORY_PATH)
    if transfer_sha_after != transfer_sha_before:
        raise RuntimeError("The immutable legacy transfer inventory changed during audit.")
    source_contract_sha_after = _require_sha256(
        SOURCE_CONTRACT_PATH,
        SOURCE_CONTRACT_SHA256,
        "Table 4.3.3 source contract",
    )
    source_pdf_sha_after = _require_sha256(
        SOURCE_PDF_PATH,
        SOURCE_PDF_SHA256,
        "Reddy source PDF",
    )
    protected_outputs_after = _require_protected_legacy_outputs()
    if protected_outputs_after != protected_outputs_before:
        raise RuntimeError("A protected legacy RLB-0 output changed during audit.")

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(output_dir / "ritz_convergence.csv", convergence_rows)
    _write_csv(output_dir / "transfer_ritz_comparison.csv", comparison_rows)
    _write_csv(output_dir / "rotary_inertia_monotonicity.csv", rotary_rows)
    _write_csv(
        output_dir / "source_inconsistency_inventory.csv", inconsistency_rows
    )
    _write_csv(output_dir / "cc_characteristic_check.csv", cc_rows)
    _write_csv(output_dir / "cf_characteristic_check.csv", cf_rows)
    (output_dir / "report.md").write_text(
        _report_text(summary, transfer_metadata),
        encoding="utf-8",
    )
    generated_hashes = {
        name: _sha256(output_dir / name)
        for name in REQUIRED_OUTPUTS
        if name != "run_manifest.json"
    }
    manifest = {
        "audit_id": "RLB-0A-S1",
        "command_mode": "smoke" if smoke else "full",
        "statuses": summary["statuses"],
        "summary": summary,
        "source_contract": {
            "path": str(SOURCE_CONTRACT_PATH.relative_to(ROOT)).replace("\\", "/"),
            "sha256": source_contract_sha_after,
            "expected_sha256": SOURCE_CONTRACT_SHA256,
            "preserved": source_contract_sha_after == SOURCE_CONTRACT_SHA256,
        },
        "source_pdf": {
            "path": str(SOURCE_PDF_PATH.relative_to(ROOT)).replace("\\", "/"),
            "sha256_before": source_pdf_sha_before,
            "sha256_after": source_pdf_sha_after,
            "expected_sha256": SOURCE_PDF_SHA256,
            "preserved": source_pdf_sha_after == source_pdf_sha_before,
        },
        "immutable_transfer_inventory": {
            **transfer_metadata,
            "sha256_after": transfer_sha_after,
            "preserved": transfer_sha_after == transfer_sha_before,
        },
        "protected_legacy_outputs": {
            "sha256_before": protected_outputs_before,
            "sha256_after": protected_outputs_after,
            "preserved": protected_outputs_after == protected_outputs_before,
        },
        "Ritz": {
            "orders": list(RITZ_ORDERS),
            "guard_order": RITZ_GUARD_ORDER,
            "guard_policy": (
                "automatic only after N10/N12 convergence or N12-quality failure; optional explicit all-case guard"
            ),
            "gauss_legendre_order": GAUSS_LEGENDRE_ORDER,
            "coordinates": "xi=x/length, W=w/length",
            "gamma": "S*length^2/D",
            "rotary_mass_ratio": "J/(m*length^2)",
            "eigenvalue": "omega^2*m*length^4/D",
            "J0_policy": "exact static condensation; no tiny J",
            "solver": "mass-orthonormal congruence plus exact back-transform",
        },
        "thresholds": {
            "transfer_vs_Ritz_relative": TRANSFER_RITZ_RELATIVE_TOLERANCE,
            "N10_vs_N12_relative": RITZ_CONVERGENCE_RELATIVE_TOLERANCE,
            "algebraic_residual": ALGEBRAIC_RESIDUAL_TOLERANCE,
            "energy_identity": ENERGY_IDENTITY_TOLERANCE,
            "natural_boundary_residual": NATURAL_BOUNDARY_RESIDUAL_TOLERANCE,
            "essential_boundary_residual": ESSENTIAL_BOUNDARY_RESIDUAL_TOLERANCE,
            "transfer_singular_ratio": TRANSFER_SINGULAR_RATIO_TOLERANCE,
            "RI_monotonicity_relative": RI_MONOTONICITY_RELATIVE_TOLERANCE,
            "source_printed_tolerance_policy": (
                "unchanged per frozen JSON; numerical allowance 5e-9"
            ),
        },
        "problem_4_14": PROBLEM_4_14_DIAGNOSTIC,
        "allowed_discrepancy_classes": sorted(DISCREPANCY_CLASSES),
        "outputs": list(REQUIRED_OUTPUTS),
        "output_sha256": generated_hashes,
        "run_manifest_self_hash_embedded": False,
        "old_RLB0_outputs_modified": False,
    }
    _write_json(output_dir / "run_manifest.json", manifest)
    return manifest


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Independent shifted-Legendre Rayleigh--Ritz audit of Reddy "
            "Table 4.3.3 direct-FSDT discrepancies."
        ),
        allow_abbrev=False,
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Run only the frozen a/h=20 cross-section of the audit.",
    )
    parser.add_argument(
        "--include-n16-guard",
        action="store_true",
        help="Evaluate the labeled N=16 guard for every selected case.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = run_audit(
        args.output_dir,
        smoke=args.smoke,
        include_n16_guard=args.include_n16_guard,
    )
    for key, value in manifest["statuses"].items():
        print(f"{key}: {value}")
    return 0 if manifest["statuses"][STATUS_OVERALL] != "FAIL" else 1


if __name__ == "__main__":
    raise SystemExit(main())
