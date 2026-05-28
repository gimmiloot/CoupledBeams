from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np
from scipy.optimize import brentq


E = 1.0
RHO = 1.0
NU = 0.3
L_SEGMENT = 1.0
ROOT_SCAN_START = 0.2
ROOT_SCAN_STEP = 0.01
ROOT_DEDUP_TOL = 1.0e-7


@dataclass(frozen=True)
class TauFactors:
    mu: float
    eta: float
    denom: float
    tau1: float
    tau2: float

    @property
    def mass_factor(self) -> float:
        return (1.0 - self.mu) * self.tau1**2 + (1.0 + self.mu) * self.tau2**2


@dataclass(frozen=True)
class Section:
    radius: float
    area: float
    inertia: float
    shear_modulus: float
    kappa: float

    @property
    def shear_stiffness(self) -> float:
        return self.kappa * self.shear_modulus * self.area

    @property
    def bending_stiffness(self) -> float:
        return E * self.inertia

    @property
    def mass_per_length(self) -> float:
        return RHO * self.area

    @property
    def rotary_inertia_per_length(self) -> float:
        return RHO * self.inertia


@dataclass(frozen=True)
class TimoshenkoBasis:
    a: float
    b: float
    h: float
    q: float
    warnings: tuple[str, ...]


@dataclass(frozen=True)
class TimoshenkoModeCoefficients:
    coeff: np.ndarray
    smallest_singular_value: float
    singular_value_ratio: float
    warnings: tuple[str, ...]


def circular_shear_coefficient(nu: float) -> float:
    return (6.0 + 12.0 * nu + 6.0 * nu**2) / (7.0 + 12.0 * nu + 4.0 * nu**2)


def tau_factors(mu: float, eta: float) -> TauFactors:
    mu_f = float(mu)
    eta_f = float(eta)
    if not (-1.0 < mu_f < 1.0):
        raise ValueError("mu must lie inside (-1, 1) for positive segment lengths.")
    if not (-1.0 < eta_f < 1.0):
        raise ValueError("eta must lie inside (-1, 1) for positive radii.")
    denom_sq = 1.0 + 2.0 * mu_f * eta_f + eta_f**2
    if denom_sq <= 0.0:
        raise ValueError("tau denominator is not positive.")
    denom = math.sqrt(denom_sq)
    return TauFactors(
        mu=mu_f,
        eta=eta_f,
        denom=denom,
        tau1=(1.0 - eta_f) / denom,
        tau2=(1.0 + eta_f) / denom,
    )


def section_from_epsilon_tau(epsilon: float, tau: float, kappa: float | None = None) -> Section:
    tau_f = float(tau)
    if tau_f <= 0.0:
        raise ValueError("tau must be positive.")
    radius = 2.0 * float(epsilon) * L_SEGMENT * tau_f
    area = math.pi * radius**2
    inertia = math.pi * radius**4 / 4.0
    shear_modulus = E / (2.0 * (1.0 + NU))
    return Section(
        radius=radius,
        area=area,
        inertia=inertia,
        shear_modulus=shear_modulus,
        kappa=circular_shear_coefficient(NU) if kappa is None else float(kappa),
    )


def section_from_epsilon(epsilon: float, kappa: float | None = None) -> Section:
    return section_from_epsilon_tau(float(epsilon), 1.0, kappa=kappa)


def project_omega(lambda_value: float, epsilon: float) -> float:
    return float(epsilon) * float(lambda_value) ** 2


def omega_cutoff(section: Section) -> float:
    return math.sqrt(section.shear_stiffness / section.rotary_inertia_per_length)


def lambda_cutoff(epsilon: float, section: Section) -> float:
    return math.sqrt(omega_cutoff(section) / float(epsilon))


def segment_lengths(mu: float) -> tuple[float, float]:
    return 1.0 - float(mu), 1.0 + float(mu)


def row_normalized(matrix: np.ndarray) -> np.ndarray:
    out = np.array(matrix, dtype=float, copy=True)
    norms = np.linalg.norm(out, axis=1)
    for index, row_norm in enumerate(norms):
        if row_norm > 0.0 and math.isfinite(float(row_norm)):
            out[index, :] /= row_norm
    return out


def normalized_det(matrix: np.ndarray) -> float:
    return float(np.linalg.det(row_normalized(matrix)))


def analytic_null_vector(matrix: np.ndarray) -> TimoshenkoModeCoefficients:
    scaled = row_normalized(np.asarray(matrix, dtype=float))
    _u, singular_values, vh = np.linalg.svd(scaled)
    coeff = vh[-1, :].astype(float)
    if coeff.size:
        pivot = int(np.argmax(np.abs(coeff)))
        if coeff[pivot] < 0.0:
            coeff = -coeff
    smallest = float(singular_values[-1])
    ratio = (
        float(singular_values[-1] / singular_values[-2])
        if len(singular_values) >= 2 and abs(float(singular_values[-2])) > 1.0e-14
        else float("nan")
    )
    return TimoshenkoModeCoefficients(
        coeff=coeff,
        smallest_singular_value=smallest,
        singular_value_ratio=ratio,
        warnings=(),
    )


def timo_basis(Lambda: float, epsilon: float, section: Section) -> TimoshenkoBasis:
    omega = project_omega(float(Lambda), float(epsilon))
    k_stiff = section.shear_stiffness
    b_stiff = section.bending_stiffness
    m = section.mass_per_length
    j = section.rotary_inertia_per_length
    c2 = omega**2 * (b_stiff * m / k_stiff + j)
    c0 = j * m * omega**4 / k_stiff - m * omega**2
    roots = np.roots([b_stiff, c2, c0])
    warnings: list[str] = []
    scale = max(1.0, float(np.max(np.abs(roots))))
    if np.max(np.abs(np.imag(roots))) > 1.0e-8 * scale:
        warnings.append(f"Timoshenko lambda^2 roots have non-negligible imaginary parts at Lambda={Lambda:.8g}.")
    z_values = np.real(roots)
    positives = [float(z) for z in z_values if z > 1.0e-12]
    negatives = [float(z) for z in z_values if z < -1.0e-12]
    if len(positives) != 1 or len(negatives) != 1:
        warnings.append(f"Expected one positive and one negative lambda^2 root at Lambda={Lambda:.8g}.")
        z_positive = float(np.max(z_values))
        z_negative = float(np.min(z_values))
        if z_positive <= 0.0 or z_negative >= 0.0:
            raise ValueError(warnings[-1])
    else:
        z_positive = positives[0]
        z_negative = negatives[0]
    a = math.sqrt(z_positive)
    b = math.sqrt(-z_negative)
    h = a + m * omega**2 / (k_stiff * a)
    q = b - m * omega**2 / (k_stiff * b)
    if abs(q) <= 1.0e-12 * max(1.0, abs(b)):
        warnings.append(f"Timoshenko trigonometric coupling q is close to zero at Lambda={Lambda:.8g}.")
    return TimoshenkoBasis(a=float(a), b=float(b), h=float(h), q=float(q), warnings=tuple(warnings))


def bending_endpoint_columns(x: float, basis: TimoshenkoBasis) -> dict[str, np.ndarray]:
    a = basis.a
    b = basis.b
    h = basis.h
    q = basis.q
    ax = a * float(x)
    bx = b * float(x)
    cosh_ax = np.cosh(ax)
    sinh_ax = np.sinh(ax)
    cos_bx = np.cos(bx)
    sin_bx = np.sin(bx)
    return {
        "w": np.array([cosh_ax - cos_bx, sinh_ax - (h / q) * sin_bx], dtype=float),
        "psi": np.array([h * sinh_ax + q * sin_bx, h * (cosh_ax - cos_bx)], dtype=float),
        "w_prime": np.array([a * sinh_ax + b * sin_bx, a * cosh_ax - (h / q) * b * cos_bx], dtype=float),
        "psi_prime": np.array([h * a * cosh_ax + q * b * cos_bx, h * (a * sinh_ax + b * sin_bx)], dtype=float),
    }


def endpoint_columns(x: float, theta: float, basis: TimoshenkoBasis, section: Section) -> dict[str, np.ndarray]:
    bend = bending_endpoint_columns(float(x), basis)
    w = np.zeros(3, dtype=float)
    psi = np.zeros(3, dtype=float)
    w_prime = np.zeros(3, dtype=float)
    psi_prime = np.zeros(3, dtype=float)
    u = np.zeros(3, dtype=float)
    u_prime = np.zeros(3, dtype=float)
    w[:2] = bend["w"]
    psi[:2] = bend["psi"]
    w_prime[:2] = bend["w_prime"]
    psi_prime[:2] = bend["psi_prime"]
    u[2] = np.sin(theta * float(x))
    u_prime[2] = theta * np.cos(theta * float(x))
    q_force = section.shear_stiffness * (w_prime - psi)
    moment = section.bending_stiffness * psi_prime
    normal = E * section.area * u_prime
    return {"w": w, "u": u, "psi": psi, "Q": q_force, "M": moment, "N": normal}


def timo_coupling_matrix(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    *,
    kappa: float | None = None,
) -> tuple[np.ndarray, tuple[str, ...]]:
    factors = tau_factors(float(mu), float(eta))
    section1 = section_from_epsilon_tau(float(epsilon), factors.tau1, kappa=kappa)
    section2 = section_from_epsilon_tau(float(epsilon), factors.tau2, kappa=kappa)
    basis1 = timo_basis(float(Lambda), float(epsilon), section1)
    basis2 = timo_basis(float(Lambda), float(epsilon), section2)
    theta = project_omega(float(Lambda), float(epsilon))
    l1, l2 = segment_lengths(factors.mu)
    rod1 = endpoint_columns(l1, theta, basis1, section1)
    rod2 = endpoint_columns(-l2, theta, basis2, section2)
    cb = math.cos(math.radians(float(beta_deg)))
    sb = math.sin(math.radians(float(beta_deg)))
    matrix = np.zeros((6, 6), dtype=float)
    matrix[0, 0:3] = rod1["w"]
    matrix[0, 3:6] = -cb * rod2["w"] + sb * rod2["u"]
    matrix[1, 0:3] = rod1["u"]
    matrix[1, 3:6] = -sb * rod2["w"] - cb * rod2["u"]
    matrix[2, 0:3] = rod1["psi"]
    matrix[2, 3:6] = -rod2["psi"]
    matrix[3, 0:3] = rod1["M"]
    matrix[3, 3:6] = -rod2["M"]
    matrix[4, 0:3] = -rod1["Q"]
    matrix[4, 3:6] = cb * rod2["Q"] - sb * rod2["N"]
    matrix[5, 0:3] = rod1["N"]
    matrix[5, 3:6] = -sb * rod2["Q"] - cb * rod2["N"]
    return matrix, basis1.warnings + basis2.warnings


def timo_mode_coefficients(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    *,
    kappa: float | None = None,
) -> TimoshenkoModeCoefficients:
    matrix, matrix_warnings = timo_coupling_matrix(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        kappa=kappa,
    )
    mode = analytic_null_vector(matrix)
    return TimoshenkoModeCoefficients(
        coeff=mode.coeff,
        smallest_singular_value=mode.smallest_singular_value,
        singular_value_ratio=mode.singular_value_ratio,
        warnings=tuple(matrix_warnings),
    )


def trapezoid_integral(values: np.ndarray, coordinates: np.ndarray) -> float:
    y_values = np.asarray(values, dtype=float)
    x_values = np.asarray(coordinates, dtype=float)
    integrate = getattr(np, "trapezoid", None)
    if integrate is not None:
        return float(integrate(y_values, x_values))
    return float(np.sum(0.5 * (y_values[1:] + y_values[:-1]) * np.diff(x_values)))


def evaluate_timo_rod_fields(
    Lambda: float,
    epsilon: float,
    section: Section,
    basis: TimoshenkoBasis,
    coeff_triplet: np.ndarray,
    x_values: np.ndarray,
) -> dict[str, np.ndarray]:
    c0, c1, axial = [float(value) for value in coeff_triplet]
    x_array = np.asarray(x_values, dtype=float)
    bend_coeff = np.array([c0, c1], dtype=float)
    theta = project_omega(float(Lambda), float(epsilon))

    w = np.empty_like(x_array, dtype=float)
    psi = np.empty_like(x_array, dtype=float)
    w_prime = np.empty_like(x_array, dtype=float)
    psi_prime = np.empty_like(x_array, dtype=float)
    for index, x_value in enumerate(x_array):
        columns = bending_endpoint_columns(float(x_value), basis)
        w[index] = float(columns["w"] @ bend_coeff)
        psi[index] = float(columns["psi"] @ bend_coeff)
        w_prime[index] = float(columns["w_prime"] @ bend_coeff)
        psi_prime[index] = float(columns["psi_prime"] @ bend_coeff)

    u = axial * np.sin(theta * x_array)
    u_prime = axial * theta * np.cos(theta * x_array)
    return {
        "x": x_array,
        "u": u,
        "w": w,
        "psi": psi,
        "u_prime": u_prime,
        "w_prime": w_prime,
        "psi_prime": psi_prime,
        "Q": section.shear_stiffness * (w_prime - psi),
        "M": section.bending_stiffness * psi_prime,
        "N": E * section.area * u_prime,
    }


def timo_mode_fields(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    *,
    coeff: np.ndarray | None = None,
    n_points: int = 801,
    kappa: float | None = None,
) -> dict[str, object]:
    factors = tau_factors(float(mu), float(eta))
    section1 = section_from_epsilon_tau(float(epsilon), factors.tau1, kappa=kappa)
    section2 = section_from_epsilon_tau(float(epsilon), factors.tau2, kappa=kappa)
    basis1 = timo_basis(float(Lambda), float(epsilon), section1)
    basis2 = timo_basis(float(Lambda), float(epsilon), section2)
    coeff_array = (
        np.asarray(coeff, dtype=float)
        if coeff is not None
        else timo_mode_coefficients(
            float(Lambda),
            float(beta_deg),
            factors.mu,
            float(epsilon),
            factors.eta,
            kappa=kappa,
        ).coeff
    )
    if coeff_array.shape != (6,):
        raise ValueError("Timoshenko mode coefficient vector must have length 6.")
    l1, l2 = segment_lengths(factors.mu)
    x1 = np.linspace(0.0, l1, int(n_points), dtype=float)
    x2 = np.linspace(0.0, -l2, int(n_points), dtype=float)
    return {
        "coeff": coeff_array,
        "factors": factors,
        "section1": section1,
        "section2": section2,
        "basis1": basis1,
        "basis2": basis2,
        "rod1": evaluate_timo_rod_fields(float(Lambda), float(epsilon), section1, basis1, coeff_array[0:3], x1),
        "rod2": evaluate_timo_rod_fields(float(Lambda), float(epsilon), section2, basis2, coeff_array[3:6], x2),
        "warnings": tuple(dict.fromkeys([*basis1.warnings, *basis2.warnings])),
    }


def _safe_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > 1.0e-30 else float("nan")


def _integral_over_rod(values: np.ndarray, x_values: np.ndarray) -> float:
    return trapezoid_integral(np.asarray(values, dtype=float), np.abs(np.asarray(x_values, dtype=float)))


def _rod_energy_terms(fields: dict[str, np.ndarray], section: Section) -> dict[str, float]:
    x_values = np.asarray(fields["x"], dtype=float)
    bending = 0.5 * section.bending_stiffness * _integral_over_rod(np.asarray(fields["psi_prime"]) ** 2, x_values)
    shear_strain = np.asarray(fields["w_prime"], dtype=float) - np.asarray(fields["psi"], dtype=float)
    shear = 0.5 * section.shear_stiffness * _integral_over_rod(shear_strain**2, x_values)
    axial = 0.5 * E * section.area * _integral_over_rod(np.asarray(fields["u_prime"]) ** 2, x_values)
    kinetic_trans = 0.5 * section.mass_per_length * _integral_over_rod(np.asarray(fields["w"]) ** 2, x_values)
    kinetic_axial = 0.5 * section.mass_per_length * _integral_over_rod(np.asarray(fields["u"]) ** 2, x_values)
    kinetic_rot = 0.5 * section.rotary_inertia_per_length * _integral_over_rod(
        np.asarray(fields["psi"]) ** 2,
        x_values,
    )
    return {
        "U_b": max(float(bending), 0.0),
        "U_s": max(float(shear), 0.0),
        "U_a": max(float(axial), 0.0),
        "T_trans": max(float(kinetic_trans), 0.0),
        "T_axial": max(float(kinetic_axial), 0.0),
        "T_rot": max(float(kinetic_rot), 0.0),
    }


def timo_energy_partition(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
    *,
    coeff: np.ndarray | None = None,
    n_points: int = 801,
    kappa: float | None = None,
) -> dict[str, float]:
    fields = timo_mode_fields(
        float(Lambda),
        float(beta_deg),
        float(mu),
        float(epsilon),
        float(eta),
        coeff=coeff,
        n_points=int(n_points),
        kappa=kappa,
    )
    rod1 = _rod_energy_terms(fields["rod1"], fields["section1"])  # type: ignore[arg-type]
    rod2 = _rod_energy_terms(fields["rod2"], fields["section2"])  # type: ignore[arg-type]
    U_b_total = rod1["U_b"] + rod2["U_b"]
    U_s_total = rod1["U_s"] + rod2["U_s"]
    U_a_total = rod1["U_a"] + rod2["U_a"]
    U_total = U_b_total + U_s_total + U_a_total
    T_trans_total = rod1["T_trans"] + rod2["T_trans"]
    T_axial_total = rod1["T_axial"] + rod2["T_axial"]
    T_rot_total = rod1["T_rot"] + rod2["T_rot"]
    T_total = T_trans_total + T_axial_total + T_rot_total
    return {
        "U_b1": rod1["U_b"],
        "U_s1": rod1["U_s"],
        "U_a1": rod1["U_a"],
        "U_b2": rod2["U_b"],
        "U_s2": rod2["U_s"],
        "U_a2": rod2["U_a"],
        "U_total": U_total,
        "U_b_total": U_b_total,
        "U_s_total": U_s_total,
        "U_a_total": U_a_total,
        "rod1_energy_fraction": _safe_ratio(rod1["U_b"] + rod1["U_s"] + rod1["U_a"], U_total),
        "rod2_energy_fraction": _safe_ratio(rod2["U_b"] + rod2["U_s"] + rod2["U_a"], U_total),
        "shear_fraction": _safe_ratio(U_s_total, U_total),
        "bending_fraction": _safe_ratio(U_b_total, U_total),
        "axial_fraction": _safe_ratio(U_a_total, U_total),
        "rod1_shear_fraction": _safe_ratio(rod1["U_s"], U_total),
        "rod2_shear_fraction": _safe_ratio(rod2["U_s"], U_total),
        "T_trans1": rod1["T_trans"],
        "T_axial1": rod1["T_axial"],
        "T_rot1": rod1["T_rot"],
        "T_trans2": rod2["T_trans"],
        "T_axial2": rod2["T_axial"],
        "T_rot2": rod2["T_rot"],
        "T_total": T_total,
        "T_trans_total": T_trans_total,
        "T_axial_total": T_axial_total,
        "T_rot_total": T_rot_total,
        "kinetic_trans_fraction": _safe_ratio(T_trans_total, T_total),
        "kinetic_axial_fraction": _safe_ratio(T_axial_total, T_total),
        "kinetic_rot_fraction": _safe_ratio(T_rot_total, T_total),
    }


def timo_det_for_scan(
    Lambda: float,
    beta_deg: float,
    mu: float,
    epsilon: float,
    eta: float = 0.0,
) -> float:
    try:
        matrix, _warnings = timo_coupling_matrix(
            float(Lambda),
            float(beta_deg),
            float(mu),
            float(epsilon),
            float(eta),
        )
    except (FloatingPointError, ValueError, OverflowError):
        return float("nan")
    if not np.all(np.isfinite(matrix)):
        return float("nan")
    return normalized_det(matrix)


def find_roots_by_sign_scan(
    func: Callable[[float], float],
    n_roots: int,
    *,
    start: float,
    upper: float,
    scan_step: float,
    grow_factor: float = 1.35,
    max_tries: int = 8,
) -> tuple[np.ndarray, list[str]]:
    warnings: list[str] = []
    roots: list[float] = []
    current_upper = float(upper)
    for _ in range(max_tries):
        roots.clear()
        grid = np.arange(float(start), current_upper + scan_step, scan_step)
        values = np.array([func(float(x)) for x in grid], dtype=float)
        for left, right, f_left, f_right in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if not (np.isfinite(f_left) and np.isfinite(f_right)):
                continue
            candidate: float | None = None
            if f_left == 0.0:
                candidate = float(left)
            elif f_left * f_right < 0.0:
                try:
                    candidate = float(brentq(func, float(left), float(right), xtol=1.0e-12, rtol=1.0e-12, maxiter=120))
                except ValueError as exc:
                    warnings.append(f"Skipped bracket [{left:.8g}, {right:.8g}]: {exc}")
                    continue
            if candidate is None:
                continue
            if roots and abs(candidate - roots[-1]) <= ROOT_DEDUP_TOL:
                continue
            roots.append(candidate)
            if len(roots) >= n_roots:
                return np.array(roots[:n_roots], dtype=float), warnings
        current_upper *= grow_factor
    warnings.append(f"Found only {len(roots)} Timoshenko roots below Lambda={current_upper:.8g}.")
    out = np.full(int(n_roots), np.nan, dtype=float)
    out[: len(roots)] = roots[:n_roots]
    return out, warnings


def timo_sorted_roots(
    beta_deg: float,
    mu: float,
    epsilon: float,
    n_roots: int,
    *,
    eta: float = 0.0,
) -> tuple[np.ndarray, list[str]]:
    factors = tau_factors(float(mu), float(eta))
    l1, l2 = segment_lengths(factors.mu)
    rod_scales = [
        math.sqrt(factors.tau1) / max(l1, 0.08),
        math.sqrt(factors.tau2) / max(l2, 0.08),
    ]
    upper = max(18.0, 1.35 * math.pi * (float(n_roots) + 4) * max(rod_scales))
    return find_roots_by_sign_scan(
        lambda value: timo_det_for_scan(
            value,
            float(beta_deg),
            factors.mu,
            float(epsilon),
            factors.eta,
        ),
        n_roots,
        start=ROOT_SCAN_START,
        upper=upper,
        scan_step=ROOT_SCAN_STEP,
    )
