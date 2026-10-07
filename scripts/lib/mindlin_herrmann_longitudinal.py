"""Diagnostic planar M-H/Timoshenko source blocks for one isotropic rectangle.

Source factors stay explicit; the separate project preset uses the frozen
rectangular Timoshenko K=5/6 contract. c is independent and dimensionless;
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
FINITE_ROD_VERSION = "jang-project-cc-bounded-basis-qr-v1"
PROJECT_RECTANGULAR_KAPPA = 5/6
PROJECT_VARIANT = "project_jang_reduced_rectangular"
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


def project_jang_reduced_rectangular(section):
    """Selected Jang closure, not a recovered Jang source numeric input.

    K is the accepted RLB-1C-ISO rectangular contract (5/6), independently
    recorded in rectangular_isotropic_models_vs_beta_note.md. Keep the
    existing section and hence bending coefficients unchanged. No Ng factors.
    """
    if section.K != PROJECT_RECTANGULAR_KAPPA:
        raise ValueError("Project preset requires the accepted rectangular K=5/6")
    return SourceModel(section, section.K, 1., 1., PROJECT_VARIANT)


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


def harmonic_state_matrix(model, omega, block="mh"):
    """exp(i*omega*t), states (u,c,N,R) or (w,theta,Q,M).

    From Hamilton variation and resultants; see the finite-rod note.
    SourceModel supplies coefficients, not an angular-joint contract.
    """
    p, nu = model.coefficients, model.section.nu
    lam = _positive(omega, "omega")**2
    if block == "mh":
        return np.array([[0, -nu, 1/p["C"], 0], [0, 0, 0, 1/p["H"]],
            [-p["m"]*lam, 0, 0, 0], [0, p["C"]*(1-nu**2)-p["j"]*lam, nu, 0]])
    if block == "timoshenko":
        return np.array([[0, 1, 1/p["S"], 0], [0, 0, 0, 1/p["B"]],
            [-p["m"]*lam, 0, 0, 0], [0, -p["r"]*lam, -1, 0]])
    raise ValueError("Unknown single-rod block")


def full_harmonic_state_matrix(model, omega):
    """Ordering (u,c,w,theta,N,R,Q,M), derived source energy blocks."""
    matrix = np.zeros((8, 8))
    matrix[np.ix_((0, 1, 4, 5), (0, 1, 4, 5))] = harmonic_state_matrix(model, omega)
    matrix[np.ix_((2, 3, 6, 7), (2, 3, 6, 7))] = harmonic_state_matrix(model, omega, "timoshenko")
    return matrix


def finite_state_basis(model, length, omega, x, block="mh", derivative=0):
    """Bounded analytic basis below the optical cutoff, no sinh/cosh.

    Spatial roots come from the PDE dispersion polynomial. Mode amplitudes
    come from the first second-order PDE, independently of state expm.
    Columns: acoustic cos/sin, left/right anchored evanescent exponentials.
    No claim to a general above-cutoff finite-spectrum solver.
    """
    length = _positive(length, "length")
    if derivative not in (0, 1):
        raise ValueError("State derivative must be 0 or 1")
    if block not in ("mh", "timoshenko"):
        raise ValueError("Unknown single-rod block")
    omega = _positive(omega, "omega")
    dispersion = blocks(model)[block == "timoshenko"]
    if omega >= 2*math.pi*dispersion.cutoff_hz:
        raise ValueError("Finite bounded basis is restricted below optical cutoff")
    roots = dispersion.spatial(omega/(2*math.pi))
    k = roots[0]["wavenumber_per_m"]
    alpha = roots[1]["attenuation_per_m"]
    if k <= 0 or alpha <= 0:
        raise ArithmeticError("Expected one acoustic and one evanescent root")
    points = np.atleast_1d(np.asarray(x, dtype=float))
    if not np.all(np.isfinite(points)) or np.any(points < 0) or np.any(points > length):
        raise ValueError("Basis coordinates must lie on the finite rod")
    p, nu = model.coefficients, model.section.nu
    decoupled = block == "mh" and nu == 0
    a = p["C"] if block == "mh" else p["S"]
    sign = 1 if block == "mh" else -1
    coupling = nu*p["C"] if block == "mh" else p["S"]
    if not decoupled:
        trig_ratio = (a*k*k-p["m"]*omega**2)/(sign*coupling*k)
        exp_ratio = (a*alpha*alpha+p["m"]*omega**2)/(sign*coupling*alpha)
        scales = np.array([1/math.hypot(1, length*trig_ratio)]*2 +
                          [1/math.hypot(1, length*exp_ratio)]*2)
    else:
        scales = np.array([1., 1., 1/length, 1/length])
    def displacements(order):
        cosine, sine = np.cos(k*points), np.sin(k*points)
        cycle = ((cosine, sine), (-sine, cosine), (-cosine, -sine))
        co, si = cycle[order]
        left = (-alpha)**order*np.exp(-alpha*points)
        right = alpha**order*np.exp(-alpha*(length-points))
        first = np.column_stack((k**order*co, k**order*si, left, right))
        if decoupled:
            first[:, 2:] = 0
            second = np.column_stack((points*0, points*0, left, right))
        else:
            second = np.column_stack((trig_ratio*k**order*si,
                -trig_ratio*k**order*co, exp_ratio*left, -exp_ratio*right))
        return first*scales, second*scales
    first, second = displacements(derivative)
    first_prime, second_prime = displacements(derivative+1)
    force = p["C"]*(first_prime+nu*second) if block == "mh" else p["S"]*(first_prime-second)
    moment = (p["H"] if block == "mh" else p["B"])*second_prime
    values = np.stack((first, second, force, moment), axis=1)
    return values[0] if np.ndim(x) == 0 else values


def finite_boundary_matrix(model, length, omega, block="mh"):
    endpoints = finite_state_basis(model, length, omega, [0., length], block)
    matrix = np.concatenate((endpoints[0, :2], endpoints[1, :2]))
    matrix[[1, 3]] *= length  # c/theta dimensionless, compare L*c with u/w
    row_norm = np.linalg.norm(matrix, axis=1)
    if np.any(row_norm == 0):
        raise ArithmeticError("Zero essential-boundary row")
    return matrix/row_norm[:, None]


def transfer_boundary_matrix(model, length, omega, block="mh", exponent_cap=1., max_steps=512):
    """Independent state shooting with exact short-step expm and QR.

    No product of exponentially ill-conditioned full transfer matrices.
    Positive QR diagonal factors preserve zeros and determinant signs.
    Impedance scaling balances displacement and force units.
    """
    from scipy.linalg import expm
    length = _positive(length, "length")
    matrix = harmonic_state_matrix(model, omega, block)
    rate = float(np.max(np.abs(np.linalg.eigvals(matrix))))
    p = model.coefficients
    elastic = p["C"] if block == "mh" else p["S"]
    gradient = p["H"] if block == "mh" else p["B"]
    scales = np.array([1., length, 1/(elastic*rate), length/(gradient*rate)])
    balanced = matrix*scales[:, None]/scales[None, :]
    steps = max(1, math.ceil(rate*length/_positive(exponent_cap, "exponent_cap")))
    if steps > max_steps:
        raise ArithmeticError("Independent shooting step budget exceeded")
    step = expm(balanced*(length/steps))
    frame = np.vstack((np.zeros((2, 2)), np.eye(2)))
    for _ in range(steps):
        frame, triangular = np.linalg.qr(step@frame, mode="reduced")
        signs = np.where(np.diag(triangular) >= 0, 1., -1.)
        frame *= signs[None, :]
    return frame[:2], {"steps": steps, "max_step_exponent": rate*length/steps,
        "step_condition": float(np.linalg.cond(step)), "state_scales": scales.tolist(),
        "method": "independent harmonic-state expm/positive-diagonal QR"}


def finite_roots(model, length, block, omega_min, omega_max, policy):
    """Bounded determinant sign search; completeness is certified separately."""
    from scipy.optimize import brentq
    evaluations = 0
    def determinant(omega):
        nonlocal evaluations
        evaluations += 1
        return float(np.linalg.det(finite_boundary_matrix(model, length, omega, block)))
    intervals = policy["scan_intervals"]
    nodes = np.linspace(omega_min, omega_max, intervals+1)
    samples = [determinant(w) for w in nodes]
    records = []
    for left, right, fl, fr in zip(nodes[:-1], nodes[1:], samples[:-1], samples[1:]):
        if fl*fr > 0:
            continue
        before = evaluations
        root, info = brentq(determinant, left, right, xtol=policy["root_xtol"],
            rtol=policy["root_rtol"], full_output=True)
        if records and abs(root-records[-1]["omega"]) <= 10*policy["root_xtol"]:
            continue
        matrix = finite_boundary_matrix(model, length, root, block)
        singular = np.linalg.svd(matrix, compute_uv=False)
        records.append({"omega": root, "frequency_hz": root/(2*math.pi),
            "bracket_omega": [float(left), float(right)], "bracket_determinants": [fl, fr],
            "evaluations": evaluations-before, "iterations": info.iterations,
            "converged": info.converged, "determinant": float(np.linalg.det(matrix)),
            "singular_ratio": float(singular[-1]/singular[0]),
            "nonzero_singular_condition": float(singular[0]/singular[-2])})
    return records, {"range_omega": [omega_min, omega_max], "scan_intervals": intervals,
        "evaluations": evaluations, "brackets_found": len(records), "failed_intervals": []}


def finite_count_upper_bound(model, length, omega, block, young_eta=.18):
    """Min-max certificate: at most this many CC eigenvalues <= omega.

    MH: a Young-inequality lower quadratic form, exact scalar Dirichlet
    spectra. Timoshenko: relax theta end constraints; exact simply-supported
    spectrum (including the uniform-rotation optical mode) is a lower bound.
    Found independent modes saturating this bound establish completeness;
    a sign scan alone is never used as a completeness claim.
    """
    p, nu = model.coefficients, model.section.nu
    if block == "mh":
        if not nu**2 < young_eta < 1:
            raise ValueError("Young certificate requires nu^2 < eta < 1")
        axial = p["C"]*(1-young_eta)/p["m"]
        normal = p["C"]*(1-nu**2/young_eta)
        axial_count = math.floor(omega*length/(math.pi*math.sqrt(axial)))
        excess = p["j"]*omega**2-normal
        contraction_count = 0 if excess < 0 else math.floor(length/math.pi*math.sqrt(excess/p["H"]))
        return {"upper_count": axial_count+contraction_count,
            "method": "Young lower form + exact two scalar Dirichlet spectra",
            "young_eta": young_eta, "axial_count": axial_count,
            "contraction_count": contraction_count, "lower_contraction_cutoff_hz": math.sqrt(normal/p["j"])/(2*math.pi)}
    if block != "timoshenko":
        raise ValueError("Unknown single-rod block")
    dispersion = blocks(model)[1]
    count = int(omega >= 2*math.pi*dispersion.cutoff_hz)  # uniform theta mode
    rows = []
    for n in range(1, 10001):
        branches = dispersion.temporal(n*math.pi/length)
        local_count = sum(r["omega_squared"] <= omega**2 for r in branches)
        if local_count == 0:
            break
        count += local_count
        rows.append({"n": n, "frequencies_hz": [r["frequency_hz"] for r in branches]})
    else:
        raise ArithmeticError("Simply-supported count budget exceeded")
    return {"upper_count": count, "method": "relaxed rotation BC / exact simply-supported spectrum",
        "uniform_rotation_cutoff_hz": dispersion.cutoff_hz, "lower_modes": rows}


def finite_mode(model, length, omega, block, order=200):
    """Mass-normalized analytical shape and scaled ODE/energy/BC diagnostics."""
    matrix = finite_boundary_matrix(model, length, omega, block)
    _, _, right = np.linalg.svd(matrix)
    coefficients = right[-1]
    nodes, weights = np.polynomial.legendre.leggauss(order)
    points, weights = (nodes+1)*length/2, weights*length/2
    values = finite_state_basis(model, length, omega, points, block)@coefficients
    gradients = finite_state_basis(model, length, omega, points, block, 1)@coefficients
    p = model.coefficients
    second_mass = p["j"] if block == "mh" else p["r"]
    mass = float(weights@(p["m"]*values[:, 0]**2+second_mass*values[:, 1]**2))
    if block == "mh":
        strain = p["C"]*(gradients[:, 0]**2+2*model.section.nu*gradients[:, 0]*values[:, 1]+values[:, 1]**2)+p["H"]*gradients[:, 1]**2
    else:
        strain = p["B"]*gradients[:, 1]**2+p["S"]*(gradients[:, 0]-values[:, 1])**2
    energy = float(weights@strain)
    endpoint = finite_state_basis(model, length, omega, [0., length], block)@coefficients
    qscale = np.array([1., length])
    boundary = float(np.max(np.abs(endpoint[:, :2]*qscale)))/float(np.max(np.abs(values[:, :2]*qscale)))
    state_matrix = harmonic_state_matrix(model, omega, block)
    rhs = values@state_matrix.T
    scale = np.maximum(np.max(np.abs(rhs), axis=0)+np.max(np.abs(gradients), axis=0), 1e-30)
    residual = float(np.max(np.abs(gradients-rhs)/scale))
    return {"coefficients": coefficients/math.sqrt(mass), "points": points,
        "weights": weights, "values": values/math.sqrt(mass),
        "diagnostics": {"boundary_scaled_residual": boundary,
            "equation_scaled_residual": residual,
            "energy_omega_squared": energy/mass,
            "energy_relative_error": abs(energy/mass/omega**2-1), "mass_before_normalization": mass}}
