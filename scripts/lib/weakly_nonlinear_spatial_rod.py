"""Isolated seven-field spatial reduced rod; exact and cubic mass-form audits.

Field order is (u,w,v,Phi,psi,theta,c), with a=(Phi,-psi,theta).
This is the explicitly adopted nonlinear reduced energy, not a new source
attribution or a modification of the verified linear production solvers.
The exact evaluator uses Rodrigues and the SO(3) right Jacobian.  The lazy
symbolic derivation uses sparse rational polynomials and needs no optional
CAS installation.  Path A varies the quartic action; Path B independently
expands matrix-exponential balances, including their coordinate covector.
No mass inverse, time integration, root finding or modal reduction is here.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
import math
from numbers import Integral, Rational
from typing import Mapping

import numpy as np


MODEL_VERSION = "seven-field-reduced-v0-quartic-action-v1"
FIELD_ORDER = ("u", "w", "v", "Phi", "psi", "theta", "c")
JET_ORDER = ("q", "qs", "qt", "qss", "qst", "qtt")
COEFFICIENT_ORDER = ("m", "jp", "jb", "C", "H", "S", "Bp", "Bb", "CT", "nu")
COEFFICIENT_UNITS = {
    "m": (1, -1, 0), "jp": (1, 1, 0), "jb": (1, 1, 0),
    "C": (1, 1, -2), "S": (1, 1, -2),
    "H": (1, 3, -2), "Bp": (1, 3, -2), "Bb": (1, 3, -2),
    "CT": (1, 3, -2), "nu": (0, 0, 0),  # exponents of kg,m,s
}
_SUFFIXES = ("", "_s", "_t", "_ss", "_st", "_tt")
SYMBOL_ORDER = tuple(f"{name}{suffix}" for suffix in _SUFFIXES for name in FIELD_ORDER) + COEFFICIENT_ORDER
_INDEX = {name: i for i, name in enumerate(SYMBOL_ORDER)}
_N_JETS = len(FIELD_ORDER) * len(JET_ORDER)
_CAP = 4


def _degree(monomial):
    return sum(index < _N_JETS for index in monomial)


class Polynomial:
    """Sparse exact ring with jet degree <=4 and symbolic constant coefficients.

    A monomial is a sorted tuple of repeated SYMBOL_ORDER indices; values are
    Fraction coefficients.  Material coefficients have amplitude degree zero.
    Every multiplication discards only degrees above four, the action's cap.
    Differentiating this capped action gives exactly the cubic residual.
    """

    __slots__ = ("terms",)

    def __init__(self, value=0):
        if isinstance(value, Polynomial):
            self.terms = dict(value.terms)
        elif isinstance(value, Mapping):
            self.terms = {tuple(key): Fraction(coef) for key, coef in value.items() if coef}
        elif isinstance(value, Rational):
            self.terms = {(): Fraction(value)} if value else {}
        else:
            raise TypeError("Polynomial constants must be exact rational numbers")

    @classmethod
    def symbol(cls, name):
        return cls({(_INDEX[name],): Fraction(1)})

    @classmethod
    def deserialize(cls, data):
        if data["symbol_order"] != list(SYMBOL_ORDER):
            raise ValueError("Polynomial symbol order differs from the model contract")
        return cls({tuple(row[0]): Fraction(row[1], row[2]) for row in data["terms"]})

    def serialize(self):
        return {"symbol_order": list(SYMBOL_ORDER), "terms": [[list(key), value.numerator, value.denominator]
                for key, value in sorted(self.terms.items(), key=lambda item: (_degree(item[0]), item[0]))]}

    def __bool__(self):
        return bool(self.terms)

    def __eq__(self, other):
        try:
            return self.terms == Polynomial(other).terms
        except TypeError:
            return False

    def __add__(self, other):
        other = Polynomial(other)
        terms = dict(self.terms)
        for key, value in other.terms.items():
            terms[key] = terms.get(key, Fraction(0)) + value
            if not terms[key]:
                del terms[key]
        return Polynomial(terms)

    __radd__ = __add__

    def __neg__(self):
        return Polynomial({key: -value for key, value in self.terms.items()})

    def __sub__(self, other):
        return self + -Polynomial(other)

    def __rsub__(self, other):
        return Polynomial(other) + -self

    def __mul__(self, other):
        other = Polynomial(other)
        terms = {}
        for left, lv in self.terms.items():
            for right, rv in other.terms.items():
                if _degree(left) + _degree(right) > _CAP:
                    continue
                key = tuple(sorted(left + right))
                terms[key] = terms.get(key, Fraction(0)) + lv * rv
        return Polynomial(terms)

    __rmul__ = __mul__

    def __truediv__(self, other):
        if not isinstance(other, Rational) or not other:
            raise TypeError("Only division by a nonzero exact rational is supported")
        return self * (Fraction(1) / other)

    def __pow__(self, power):
        if not isinstance(power, Integral) or power < 0:
            raise ValueError("Polynomial power must be a nonnegative integer")
        result = Polynomial(1)
        for _ in range(power):
            result = result * self
        return result

    def truncate(self, degree):
        return Polynomial({key: value for key, value in self.terms.items() if _degree(key) <= degree})

    def homogeneous(self, degree):
        return Polynomial({key: value for key, value in self.terms.items() if _degree(key) == degree})

    def derivative(self, symbol):
        index = _INDEX[symbol] if isinstance(symbol, str) else int(symbol)
        terms = {}
        for key, value in self.terms.items():
            count = key.count(index)
            if count:
                remaining = list(key)
                remaining.remove(index)
                terms[tuple(remaining)] = count * value
        return Polynomial(terms)

    def total_derivative(self, axis):
        if axis not in ("s", "t"):
            raise ValueError("Total derivative axis must be s or t")
        offsets = {"s": {0: 7, 7: 21, 14: 28}, "t": {0: 14, 7: 28, 14: 35}}[axis]
        result = Polynomial()
        active = set(index for key in self.terms for index in key if index < _N_JETS)
        for index in active:
            block = index // 7 * 7
            if block not in offsets:
                raise ValueError("A third derivative jet would be required")
            target = offsets[block] + index % 7
            result = result + self.derivative(index) * Polynomial({(target,): Fraction(1)})
        return result

    def substitute(self, mapping):
        replacements = {_INDEX[name]: Polynomial(value) for name, value in mapping.items()}
        result = Polynomial()
        for key, value in self.terms.items():
            term = Polynomial(value)
            for index in key:
                term = term * replacements.get(index, Polynomial({(index,): Fraction(1)}))
            result = result + term
        return result

    def evaluate(self, values):
        values = [values[name] if name in values else None for name in SYMBOL_ORDER]
        result = 0.0
        for key, coefficient in self.terms.items():
            term = float(coefficient)
            for index in key:
                if values[index] is None:
                    raise ValueError(f"Missing polynomial value: {SYMBOL_ORDER[index]}")
                term *= values[index]
            result += term
        return result

    def __str__(self):
        if not self.terms:
            return "0"
        chunks = []
        for key, coefficient in sorted(self.terms.items(), key=lambda item: (_degree(item[0]), item[0])):
            factors = []
            for index in sorted(set(key)):
                count = key.count(index)
                factors.append(SYMBOL_ORDER[index] + (f"**{count}" if count > 1 else ""))
            chunks.append(str(coefficient) + ("*" + "*".join(factors) if factors else ""))
        return " + ".join(chunks).replace("+ -", "- ")


def _cross(a, b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])


def _add(a, b):
    return tuple(x+y for x, y in zip(a, b))


def _scale(scalar, vector):
    return tuple(scalar*x for x in vector)


def _dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Polynomial())


def _hat(a):
    zero = Polynomial()
    return ((zero, -a[2], a[1]), (a[2], zero, -a[0]), (-a[1], a[0], zero))


def _identity():
    return tuple(tuple(Polynomial(int(i == j)) for j in range(3)) for i in range(3))


def _transpose(matrix):
    return tuple(zip(*matrix))


def _matmul(a, b):
    return tuple(tuple(_dot(row, col) for col in _transpose(b)) for row in a)


def _matvec(a, b):
    return tuple(_dot(row, b) for row in a)


def _matrix_add(a, b):
    return tuple(_add(row, other) for row, other in zip(a, b))


def _matrix_scale(coef, matrix):
    return tuple(_scale(coef, row) for row in matrix)


def _truncate_vector(vector, degree):
    return tuple(value.truncate(degree) for value in vector)


def _vee(matrix):
    # Taking the antisymmetric part retains independently obtained series
    # cancellation evidence rather than assuming the product skew in code.
    return ((matrix[2][1]-matrix[1][2])/2, (matrix[0][2]-matrix[2][0])/2,
            (matrix[1][0]-matrix[0][1])/2)


@dataclass(frozen=True)
class SymbolicModel:
    symbols: Mapping[str, Polynomial]
    T4: Polynomial
    V4: Polynomial
    residual_a: tuple[Polynomial, ...]
    residual_b: tuple[Polynomial, ...]
    flux_a: tuple[Polynomial, ...]
    flux_b: tuple[Polynomial, ...]
    body_b: tuple[Polynomial, ...]
    gamma_a: tuple[Polynomial, ...]
    chi_a: tuple[Polynomial, ...]
    omega_a: tuple[Polynomial, ...]
    gamma_b: tuple[Polynomial, ...]
    chi_b: tuple[Polynomial, ...]
    omega_b: tuple[Polynomial, ...]
    right_jacobian_b: tuple[tuple[Polynomial, ...], ...]
    rotation_b: tuple[tuple[Polynomial, ...], ...]

    @property
    def L4(self):
        return self.T4 - self.V4

    def comparisons(self):
        return [{"field": name, "degree": degree,
                 "difference_terms": (a-b).homogeneous(degree).serialize()["terms"],
                 "status": "PASS" if not (a-b).homogeneous(degree) else "FAIL"}
                for name, a, b in zip(FIELD_ORDER, self.residual_a, self.residual_b)
                for degree in (1, 2, 3)]


@lru_cache(maxsize=1)
def derive_polynomials():
    """Generate quartic action and two independent cubic derivations once.

    A uses vector-series kinematic measures before variation.  B obtains
    chi/Omega/Jr from derivatives of the matrix exponential itself, then
    applies the body balances and coordinate virtual-work transformation.
    No A residual or A flux is an input to B.
    """
    p = {name: Polynomial.symbol(name) for name in SYMBOL_ORDER}
    q = tuple(p[name] for name in FIELD_ORDER)
    qs = tuple(p[name+"_s"] for name in FIELD_ORDER)
    qt = tuple(p[name+"_t"] for name in FIELD_ORDER)
    qtt = tuple(p[name+"_tt"] for name in FIELD_ORDER)
    a = (q[3], -q[4], q[5])
    ass = (qs[3], -qs[4], qs[5])
    at = (qt[3], -qt[4], qt[5])
    e1 = (Polynomial(1), Polynomial(), Polynomial())
    # A: explicit vector kinematic expansion, through third order.
    base = _add(e1, qs[:3])
    gamma_a = _truncate_vector(_add(_add(_add(qs[:3], _scale(-1, _cross(a, base))),
                       _scale(Fraction(1, 2), _cross(a, _cross(a, base)))),
                       _scale(Fraction(-1, 6), _cross(a, _cross(a, _cross(a, base))))), 3)
    def axial_series(direction):
        return _add(_add(direction, _scale(Fraction(-1, 2), _cross(a, direction))),
                    _scale(Fraction(1, 6), _cross(a, _cross(a, direction))))
    chi_a, omega_a = axial_series(ass), axial_series(at)
    c = q[6]
    inertia = (p["jb"]+p["jp"]*(1+c)**2, p["jb"], p["jp"]*(1+c)**2)
    stiffness = (p["C"], p["S"], p["S"])
    curvature = (p["CT"], p["Bb"], p["Bp"])
    T4 = (p["m"]*_dot(qt[:3], qt[:3])+p["jp"]*qt[6]**2 +
          sum((j*x*x for j, x in zip(inertia, omega_a)), Polynomial()))/2
    V4 = (sum((d*x*x for d, x in zip(stiffness, gamma_a)), Polynomial())+
          p["C"]*c*c+p["H"]*qs[6]**2+
          sum((d*x*x for d, x in zip(curvature, chi_a)), Polynomial()))/2 + p["nu"]*p["C"]*c*gamma_a[0]
    L4 = T4-V4
    residual_a = tuple(L4.derivative(name+"_t").total_derivative("t") +
                       L4.derivative(name+"_s").total_derivative("s")-L4.derivative(name)
                       for name in FIELD_ORDER)
    flux_a = tuple(V4.derivative(name+"_s") for name in FIELD_ORDER)
    # B: R=exp(K) polynomial, matrix differentiation R^T R_{s,t}.
    K = _hat(a)
    R = _identity()
    power = _identity()
    for n in range(1, 4):
        power = _matmul(power, K)
        R = _matrix_add(R, _matrix_scale(Fraction(1, math.factorial(n)), power))
    Rt = _transpose(R)
    R_s = tuple(tuple(x.total_derivative("s") for x in row) for row in R)
    R_t = tuple(tuple(x.total_derivative("t") for x in row) for row in R)
    chi_b = _truncate_vector(_vee(_matmul(Rt, R_s)), 3)
    omega_b = _truncate_vector(_vee(_matmul(Rt, R_t)), 3)
    gamma_b = _truncate_vector(_add(_matvec(Rt, base), _scale(-1, e1)), 3)
    # Differentiation with respect to a uses the sign a2=-psi explicitly.
    da_names = ("Phi", "psi", "theta")
    columns = []
    for name, sign in zip(da_names, (1, -1, 1)):
        R_a = tuple(tuple(sign*x.derivative(name) for x in row) for row in R)
        columns.append(_truncate_vector(_vee(_matmul(Rt, R_a)), 2))
    Jr = _transpose(tuple(columns))
    N = _add(tuple(d*x for d, x in zip(stiffness, gamma_b)),
             (p["nu"]*p["C"]*c, Polynomial(), Polynomial()))
    M = tuple(d*x for d, x in zip(curvature, chi_b))
    ell = _truncate_vector(tuple(j*x for j, x in zip(inertia, omega_b)), 3)
    g = _add(e1, gamma_b)
    body_b = _truncate_vector(_add(_add(_add(tuple(x.total_derivative("t") for x in ell),
                 _cross(omega_b, ell)), _scale(-1, tuple(x.total_derivative("s") for x in M))),
                 _scale(-1, _add(_cross(chi_b, M), _cross(g, N)))), 3)
    F = _truncate_vector(_matvec(R, N), 3)
    ez = _truncate_vector(_matvec(_transpose(Jr), body_b), 3)
    ez = (ez[0], -ez[1], ez[2])
    ec = (p["jp"]*qtt[6]-p["H"]*p["c_ss"]+p["C"]*(c+p["nu"]*gamma_b[0])-
          p["jp"]*(1+c)*(omega_b[0]**2+omega_b[2]**2)).truncate(3)
    residual_b = tuple(p["m"]*qtt[i]-F[i].total_derivative("s") for i in range(3)) + ez + (ec,)
    pz = _truncate_vector(_matvec(_transpose(Jr), M), 3)
    flux_b = F + (pz[0], -pz[1], pz[2], p["H"]*qs[6])
    return SymbolicModel(p, T4, V4, residual_a, residual_b, flux_a, flux_b, body_b,
                         gamma_a, chi_a, omega_a, gamma_b, chi_b, omega_b, Jr, R)


@dataclass(frozen=True)
class RodCoefficients:
    m: float
    jp: float
    jb: float
    C: float
    H: float
    S: float
    Bp: float
    Bb: float
    CT: float
    nu: float

    def __post_init__(self):
        for name in COEFFICIENT_ORDER[:-1]:
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not math.isfinite(self.nu) or not -1 < self.nu < 0.5:
            raise ValueError("nu must lie in the physical isotropic interval (-1,0.5)")

    @classmethod
    def rectangular(cls, E, rho, nu, b, h, torsional_stiffness, kappa=5/6):
        """Original-section coefficients; CT must be explicitly supplied.

        CT is not silently replaced by G*Ip. kappa is the project value;
        accepting an explicit positive value supports algebraic diagnostics,
        not fitting or a new source prescription.
        """
        for name, value in (("E", E), ("rho", rho), ("b", b), ("h", h), ("kappa", kappa)):
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not -1 < nu < 0.5:
            raise ValueError("nu must lie in (-1,0.5)")
        area, ip, ib = b*h, b*h**3/12, h*b**3/12
        G = E/(2*(1+nu))
        return cls(rho*area, rho*ip, rho*ib, E*area/(1-nu**2), kappa*G*ip,
                   kappa*G*area, E*ip, E*ib, torsional_stiffness, nu)

    def values(self):
        return {name: getattr(self, name) for name in COEFFICIENT_ORDER}


@dataclass(frozen=True)
class FieldJet:
    q: np.ndarray
    qs: np.ndarray
    qt: np.ndarray
    qss: np.ndarray
    qst: np.ndarray
    qtt: np.ndarray

    def __post_init__(self):
        for name in JET_ORDER:
            value = np.asarray(getattr(self, name), dtype=float)
            if value.shape != (7,) or not np.all(np.isfinite(value)):
                raise ValueError(f"{name} must contain seven finite real values")
            object.__setattr__(self, name, value)

    def scaled(self, amplitude):
        return FieldJet(*(getattr(self, name)*amplitude for name in JET_ORDER))

    def values(self):
        return {field+suffix: float(value) for name, suffix in zip(JET_ORDER, _SUFFIXES)
                for field, value in zip(FIELD_ORDER, getattr(self, name))}


def skew(vector):
    a = np.asarray(vector)
    return np.array(((0, -a[2], a[1]), (a[2], 0, -a[0]), (-a[1], a[0], 0)), dtype=a.dtype)


def _so3_coefficients(t):
    """A,B,D and their analytic derivatives with respect to |a|^2."""
    if t < 1e-3:
        # Analytic entire series: avoids both cancellation and norm division.
        values, derivatives = [], []
        for offset in (1, 2, 3):
            values.append(sum((-1)**n*t**n/math.factorial(2*n+offset) for n in range(10)))
            derivatives.append(sum(n*(-1)**n*t**(n-1)/math.factorial(2*n+offset) for n in range(1, 10)))
        return tuple(values + derivatives)
    r = math.sqrt(t)
    sine, cosine = math.sin(r), math.cos(r)
    return (sine/r, (1-cosine)/t, (r-sine)/(r*t),
            (r*cosine-sine)/(2*r**3), (r*sine-2*(1-cosine))/(2*r**4),
            (r*(1-cosine)-3*(r-sine))/(2*r**5))


def rotation_and_right_jacobian(a, direction=None):
    """Exact Rodrigues/Jr and optional analytic directional derivatives.

    The ten-term entire-series branch is only a stable evaluation near zero;
    it is not the cubic model.  No finite-difference derivative is used.
    """
    a = np.asarray(a, dtype=float)
    if a.shape != (3,) or not np.all(np.isfinite(a)):
        raise ValueError("Rotation vector must contain three finite values")
    A, B, D, Ap, Bp, Dp = _so3_coefficients(float(a@a))
    K = skew(a)
    K2 = K@K
    R = np.eye(3)+A*K+B*K2
    Jr = np.eye(3)-B*K+D*K2
    if direction is None:
        return R, Jr
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (3,) or not np.all(np.isfinite(direction)):
        raise ValueError("Rotation direction must contain three finite values")
    dK = skew(direction)
    dt = 2*float(a@direction)
    dK2 = dK@K+K@dK
    dR = Ap*dt*K+A*dK+Bp*dt*K2+B*dK2
    dJr = -Bp*dt*K-B*dK+Dp*dt*K2+D*dK2
    return R, Jr, dR, dJr


def full_evaluate(jet: FieldJet, coefficients: RodCoefficients):
    """Full untruncated reduced energies, balances and coordinate fluxes.

    Residual order follows FIELD_ORDER; flux order is conjugate to q_s.
    The body moment residual is separately exposed so that dropping Jr^T
    can be detected by negative controls.  Boundary moment flux p_z is a
    coordinate covector, not the physical moment before transformation.
    """
    p = coefficients
    P = np.diag((1., -1., 1.))
    a, ass, at = P@jet.q[3:6], P@jet.qs[3:6], P@jet.qt[3:6]
    R, Jr, R_s, J_s = rotation_and_right_jacobian(a, ass)
    _, _, _, J_t = rotation_and_right_jacobian(a, at)
    e1 = np.array((1., 0., 0.))
    g = R.T@(e1+jet.qs[:3])
    gamma = g-e1
    gamma_s = R_s.T@(e1+jet.qs[:3])+R.T@jet.qss[:3]
    chi, omega = Jr@ass, Jr@at
    chi_s = J_s@ass+Jr@(P@jet.qss[3:6])
    omega_t = J_t@at+Jr@(P@jet.qtt[3:6])
    c, cs, ct = jet.q[6], jet.qs[6], jet.qt[6]
    diagonal = np.array((p.jb+p.jp*(1+c)**2, p.jb, p.jp*(1+c)**2))
    diagonal_c = np.array((2*p.jp*(1+c), 0., 2*p.jp*(1+c)))
    DG, DK = np.array((p.C, p.S, p.S)), np.array((p.CT, p.Bb, p.Bp))
    N, N_s = DG*gamma, DG*gamma_s
    N[0] += p.nu*p.C*c
    N_s[0] += p.nu*p.C*cs
    M, M_s = DK*chi, DK*chi_s
    ell = diagonal*omega
    ell_t = diagonal*omega_t+diagonal_c*ct*omega
    body = ell_t+np.cross(omega, ell)-M_s-np.cross(chi, M)-np.cross(g, N)
    Eu = p.m*jet.qtt[:3]-R_s@N-R@N_s
    Ez = P@Jr.T@body
    Ec = p.jp*jet.qtt[6]-p.H*jet.qss[6]+p.C*(c+p.nu*gamma[0])-p.jp*(1+c)*(omega[0]**2+omega[2]**2)
    F, pz, Rc = R@N, P@Jr.T@M, p.H*cs
    T = (p.m*float(jet.qt[:3]@jet.qt[:3])+p.jp*ct**2+float(omega@(diagonal*omega)))/2
    V = (float(gamma@(DG*gamma))+p.C*c*c+p.H*cs*cs+float(chi@(DK*chi)))/2+p.nu*p.C*c*gamma[0]
    mass = np.zeros((7, 7))
    mass[:3, :3] = p.m*np.eye(3)
    mass[3:6, 3:6] = P@Jr.T@np.diag(diagonal)@Jr@P
    mass[6, 6] = p.jp
    return {"residual": np.concatenate((Eu, Ez, (Ec,))), "flux": np.concatenate((F, pz, (Rc,))),
            "T": T, "V": V, "R": R, "Jr": Jr, "Gamma": gamma, "chi": chi, "Omega": omega,
            "body_residual": body, "N": N, "M_body": M, "Rc": Rc, "mass_matrix": mass}


def polynomial_evaluate(jet: FieldJet, coefficients: RodCoefficients, model=None, path="a"):
    """Quartic-action/cubic-coordinate evaluator; no hidden full-model label."""
    model = derive_polynomials() if model is None else model
    if path not in ("a", "b"):
        raise ValueError("Polynomial path must be a or b")
    values = jet.values() | coefficients.values()
    residuals = model.residual_a if path == "a" else model.residual_b
    fluxes = model.flux_a if path == "a" else model.flux_b
    return {"residual": np.array([x.evaluate(values) for x in residuals]),
            "flux": np.array([x.evaluate(values) for x in fluxes]),
            "T": model.T4.evaluate(values), "V": model.V4.evaluate(values)}


def quartic_mass_matrix(jet: FieldJet, coefficients: RodCoefficients, model=None):
    model = derive_polynomials() if model is None else model
    values = jet.values() | coefficients.values()
    return np.array([[model.T4.derivative(a+"_t").derivative(b+"_t").evaluate(values)
                      for b in FIELD_ORDER] for a in FIELD_ORDER])
