"""Diagnostic planar M-H/Timoshenko source blocks for one isotropic rectangle.

No production coefficient defaults. c is independent and dimensionless;
q=(u,c,w,theta). Source I is I_y=b*h**3/12, never the polar I_p.
Rucka (2)--(13), Jang (1),(5),(6),(A1)--(A9); see the canonical theory note.
The existing rectangular helper supplies all Timoshenko section coefficients.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import (
    SectionProperties, rectangular_section,
)

EQUATIONS_VERSION = "planar-mh-tim-source-energy-v1"
DOF_ORDER = ("u", "c", "w", "theta")
STRAIN_ORDER = ("u_x", "c", "c_x", "w_x-theta", "theta_x")
COEFFICIENT_UNITS = {
    "C": (1, 1, -2), "H": (1, 3, -2), "m": (1, -1, 0),
    "j": (1, 1, 0), "B": (1, 3, -2), "S": (1, 1, -2),
    "r": (1, 1, 0),  # exponents of kg,m,s
}


def _positive(value, name):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return number


def rectangle_moments(b, h):
    """Centered moments, also supports exact Fraction inputs."""
    if b <= 0 or h <= 0:
        raise ValueError("Rectangle dimensions must be positive")
    area = b*h
    iy, iz = b*h**3/12, h*b**3/12
    return {"A": area, "Qy": area*0, "Qz": area*0, "Iyz": area*0,
            "Iy": iy, "Iz": iz, "Ip": iy+iz}


@dataclass(frozen=True)
class SourceModel:
    section: SectionProperties
    mh_shear_factor: float
    mh_inertia_factor: float
    tim_rotary_factor: float
    variant: str

    def __post_init__(self):
        if self.section.section_kind != "rectangle_width_ey_thickness_ez":
            raise ValueError("This audit accepts rectangular sections only")
        for name in ("mh_shear_factor", "mh_inertia_factor", "tim_rotary_factor"):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        if not self.variant:
            raise ValueError("An explicit diagnostic variant identifier is required")

    @property
    def coefficients(self):
        s = self.section
        return {"C": s.EA/(1-s.nu**2), "H": self.mh_shear_factor*s.G*s.inertia,
                "m": s.rhoA, "j": self.mh_inertia_factor*s.rhoI,
                "B": s.EI, "S": s.KGA, "r": self.tim_rotary_factor*s.rhoI}


def source_model(parameters, *, mh_shear_factor, mh_inertia_factor,
                 tim_shear_factor, tim_rotary_factor, variant):
    section = rectangular_section(E=parameters["E"], nu=parameters["nu"],
        rho=parameters["rho"], width=parameters["b"], thickness=parameters["h"],
        K=tim_shear_factor)
    return SourceModel(section, mh_shear_factor, mh_inertia_factor,
                       tim_rotary_factor, variant)


def strain_operator(k):
    """Rucka D, with exp(i(k*x-omega*t)); signs local to this audit."""
    return np.array([[1j*k, 0, 0, 0], [0, 1, 0, 0], [0, 1j*k, 0, 0],
                     [0, 0, 1j*k, -1], [0, 0, 0, 1j*k]], dtype=complex)


def energy_matrices(model):
    p, nu = model.coefficients, model.section.nu
    elastic = np.diag([p["C"], p["C"], p["H"], p["S"], p["B"]])
    elastic[0, 1] = elastic[1, 0] = nu*p["C"]
    mass = np.diag([p["m"], p["j"], p["m"], p["r"]])
    return elastic, mass


def fourier_matrices(model, k):
    elastic, mass = energy_matrices(model)
    d = strain_operator(k)
    return d.conj().T @ elastic @ d, mass


def boundary_quantities(model, *, ux, c, cx, wx, theta, thetax):
    p = model.coefficients
    return {"N": p["C"]*(ux+model.section.nu*c), "R": p["H"]*cx,
            "Q": p["S"]*(wx-theta), "M": p["B"]*thetax}


@dataclass(frozen=True)
class DispersionBlock:
    """Two source branches only; polynomial in lambda=omega**2, s=k**2."""
    a1: float
    d0: float
    d1: float
    b1: float
    p1: float
    p2: float
    labels: tuple[str, str]

    @property
    def cutoff_hz(self):
        return math.sqrt(self.d0)/(2*math.pi)

    def temporal(self, k):
        """Stable roots and analytic d omega/d k. No finite differences."""
        k = float(k)
        if not math.isfinite(k) or k < 0:
            raise ValueError("k must be finite and nonnegative")
        s = k*k
        difference = (self.a1-self.d1)*s-self.d0
        gap = math.hypot(difference, 2*self.b1*k)
        trace = self.d0+(self.a1+self.d1)*s
        high = (trace+gap)/2
        product = s*(self.p1+self.p2*s)
        low = product/high
        if gap == 0:
            raise ArithmeticError("Degenerate branch derivative is not defined")
        high_prime = ((self.a1+self.d1)*2*k +
            (difference*(self.a1-self.d1)*2*k+4*self.b1**2*k)/gap)/2
        product_prime = 2*k*(self.p1+2*self.p2*s)
        low_prime = (product_prime-low*high_prime)/high
        lambdas = (low, high)
        velocities = (math.sqrt(self.p1/self.d0) if k == 0 else
                      low_prime/(2*math.sqrt(low)), high_prime/(2*math.sqrt(high)))
        return [{"branch": label, "omega_squared": value,
                 "frequency_hz": math.sqrt(value)/(2*math.pi),
                 "group_velocity_m_s": velocity}
                for label, value, velocity in zip(self.labels, lambdas, velocities)]

    def spatial(self, frequency_hz):
        """Both spatial roots, including evanescence, via stable quadratic."""
        f = float(frequency_hz)
        if not math.isfinite(f) or f < 0:
            raise ValueError("Frequency must be finite and nonnegative")
        lam = (2*math.pi*f)**2
        if abs(lam-self.d0) <= 32*np.finfo(float).eps*max(lam, self.d0):
            lam = self.d0  # analytic cutoff limit, rather than a tiny fake k
        a, b, c = self.p2, self.p1-(self.a1+self.d1)*lam, lam*(lam-self.d0)
        discriminant = b*b-4*a*c
        if discriminant < -64*np.finfo(float).eps*(b*b+abs(4*a*c)):
            raise ArithmeticError("Negative spatial-root discriminant")
        q = -.5*(b+math.copysign(math.sqrt(max(0., discriminant)), b))
        roots = (0., 0.) if q == 0 else (q/a, c/q)
        roots = sorted(roots, reverse=True)  # acoustic has larger k^2 here
        rows = []
        for index, (label, s) in enumerate(zip(self.labels, roots)):
            scale = max(abs(roots[0]), abs(roots[1]), 1.)
            if abs(s) <= 64*np.finfo(float).eps*scale:
                s = 0.
            k, attenuation = math.sqrt(max(s, 0.)), math.sqrt(max(-s, 0.))
            propagating = s > 0
            group = self.temporal(k)[index]["group_velocity_m_s"] if propagating else None
            if s == 0:
                group = self.temporal(0.)[index]["group_velocity_m_s"]
            terms = (a*s*s, b*s, c)
            residual = abs(sum(terms))/max(sum(abs(v) for v in terms), 1.)
            rows.append({"branch": label, "frequency_hz": f, "k_squared_per_m2": s,
                "wavenumber_per_m": k, "attenuation_per_m": attenuation,
                "group_velocity_m_s": group,
                "phase_velocity_m_s": (2*math.pi*f/k if propagating else None),
                "state": "PROPAGATING" if propagating else
                         "EVANESCENT" if s < 0 else "ZERO_OR_CUTOFF",
                "polynomial_scaled_residual": residual})
        return rows


def blocks(model):
    p, nu = model.coefficients, model.section.nu
    a, d0, d1 = p["C"]/p["m"], p["C"]/p["j"], p["H"]/p["j"]
    axial = DispersionBlock(a, d0, d1, nu*p["C"]/math.sqrt(p["m"]*p["j"]),
                           a*d0*(1-nu**2), a*d1, ("axial", "contraction"))
    a, d0, d1 = p["S"]/p["m"], p["S"]/p["r"], p["B"]/p["r"]
    bending = DispersionBlock(a, d0, d1, p["S"]/math.sqrt(p["m"]*p["r"]),
                             0., a*d1, ("bending", "shear"))
    return axial, bending


def limits(model):
    p = model.coefficients
    speed = math.sqrt(model.section.E/model.section.rho)
    return {"axial_acoustic_speed_m_s": speed,
        "contraction_cutoff_hz": blocks(model)[0].cutoff_hz,
        "shear_cutoff_hz": blocks(model)[1].cutoff_hz,
        "axial_omega2_k4_coefficient": model.section.nu**2*(p["H"]-p["j"]*speed**2)/p["m"],
        "bending_omega_over_k2_m2_s": math.sqrt(p["B"]/p["m"]),
        "mh_high_k_speeds_m_s": sorted([math.sqrt(p["C"]/p["m"]), math.sqrt(p["H"]/p["j"])]),
        "tim_high_k_speeds_m_s": sorted([math.sqrt(p["S"]/p["m"]), math.sqrt(p["B"]/p["r"])])}
