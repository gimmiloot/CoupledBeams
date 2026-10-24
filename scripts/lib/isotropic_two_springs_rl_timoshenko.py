"""Dimensional isotropic RL/Tim diagnostic with TWO exact/finite springs.

Canonical note: docs/theory/isotropic_two_springs_rl_timoshenko.md.
State (u,w,psi,N,Q,M), both x from outer support to joint. Existing arm
bases/geometry remain unchanged. No damping, laminates or production dispatch.
"""
from dataclasses import dataclass
from functools import cached_property
import math
import numpy as np
from scipy.linalg import block_diag, cholesky, eigh, null_space, solve_triangular
from scipy.optimize import brentq
from scripts.lib import bishop_longitudinal as axial
from scripts.lib import mindlin_herrmann_longitudinal as tim
from scripts.lib import mindlin_herrmann_timoshenko_joint as geometry
from scripts.lib import contraction_joint_sensitivity as old_joint
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section
from scripts.lib.reddy_symmetric_coupled_beams_ritz import shifted_legendre_clamped_basis

RIGID = "RIGID"
VARIANTS = ("RL_FULL", "RL_PLANAR_LEGACY", "ELEMENTARY_LEGACY")


@dataclass(frozen=True)
class Model:
    epsilon: float = .005
    nu: float = .3
    aspect: float = 4.
    E: float = 1.
    rho: float = 1.
    l: float = 1.
    ks: float = 5/6
    variant: str = "RL_FULL"

    def __post_init__(self):
        if self.variant not in VARIANTS or min(self.epsilon, self.aspect, self.E,
                self.rho, self.l, self.ks) <= 0 or not -1 < self.nu < .5:
            raise ValueError("Invalid explicit isotropic model")

    @cached_property
    def section(self):
        h = math.sqrt(12)*self.epsilon*self.l
        return rectangular_section(E=self.E, rho=self.rho, nu=self.nu,
                                   width=self.aspect*h, thickness=h, K=self.ks)

    @cached_property
    def carrier(self):
        # Only the unchanged Timoshenko block is called; no MH axial physics.
        return tim.SourceModel(self.section, self.ks, 1., 1., "two_springs_tim_carrier")

    @property
    def J(self):
        return self.nu**2*self.section.rhoI * {
            "RL_FULL": 1+self.aspect**2, "RL_PLANAR_LEGACY": 1., "ELEMENTARY_LEGACY": 0.}[self.variant]

    @property
    def time(self):
        return self.l**2*math.sqrt(self.section.rhoA/self.section.EI)

    def omega(self, lam):
        return np.asarray(lam)**2/self.time

    def Lambda(self, omega):
        return np.sqrt(np.asarray(omega)*self.time)


@dataclass(frozen=True)
class Case:
    beta: float
    mu: float = .01
    k_joint: float | str = RIGID  # dimensionless k*l/EI
    k_end: float | str = RIGID

    def __post_init__(self):
        if not math.isfinite(self.beta) or not abs(self.mu) < 1:
            raise ValueError("Finite angle and positive lengths required")
        for k in (self.k_joint, self.k_end):
            if k != RIGID and (not math.isfinite(float(k)) or float(k) < 0):
                raise ValueError("Use nonnegative finite stiffness or exact RIGID")

    def lengths(self, m):
        return (m.l*(1-self.mu), m.l*(1+self.mu))


def basis(m, length, omega, x, derivative=0):
    s = m.section
    segment = axial.Segment(length, s.EA, s.rhoA, H=0., J=m.J)
    points = np.atleast_1d(x)
    b = np.zeros((len(points), 6, 6))
    b[:, 0, :2] = axial.basis(segment, omega, points, derivative)
    b[:, 3, :2] = (s.EA-m.J*omega**2)*axial.basis(segment, omega, points, derivative+1)
    b[:, (1, 2, 4, 5), 2:] = tim.finite_state_basis(m.carrier, length, omega,
                                                   points, "timoshenko", derivative)
    return b[0] if np.ndim(x) == 0 else b


def state_matrix(m, omega):
    h = np.zeros((6, 6))
    h[0, 3], h[3, 0] = 1/(m.section.EA-m.J*omega**2), -m.section.rhoA*omega**2
    h[np.ix_((1, 2, 4, 5), (1, 2, 4, 5))] = tim.harmonic_state_matrix(m.carrier, omega, "timoshenko")
    return h


def joint_matrix(m, case):
    c = old_joint.Case(case.beta, variant="elementary",
                       rotation="RIGID" if case.k_joint == RIGID else "SPRING")
    k = 0. if case.k_joint == RIGID else case.k_joint*m.section.EI/m.l
    return old_joint.joint_matrix(c, k)


def external_matrix(m, case, arm):
    e = np.eye(6)[:3].copy()
    if arm == 1 and case.k_end != RIGID:
        e[2] = 0.
        e[2, 2], e[2, 5] = case.k_end*m.section.EI/m.l, -1.
    return e


def boundary(m, case, omega):
    ends = [basis(m, L, omega, [0., L]) for L in case.lengths(m)]
    b = np.zeros((12, 12))
    for i, values in enumerate(ends):
        b[3*i:3*i+3, 6*i:6*i+6] = external_matrix(m, case, i)@values[0]
    b[6:] = joint_matrix(m, case)@block_diag(*(v[1] for v in ends))
    return b


def balance(b):
    """Positive row/column equilibration, with explicit coefficient pullback."""
    row = 1/np.linalg.norm(b, axis=1)
    first = row[:, None]*b
    col = 1/np.linalg.norm(first, axis=0)
    return first*col, row, col


def trace_map(m, case):
    # All endpoints q=(q10,q1L,q20,q2L); external translations fixed.
    nr = 1 if case.k_joint == RIGID else 2
    ne = int(case.k_end != RIGID)
    p = np.zeros((12, 2+nr+ne))
    for i, f in enumerate(geometry.frames(case.beta)):
        p[6*i+3, :2], p[6*i+4, :2] = f.t, f.n
        p[6*i+5, 2+(i if nr == 2 else 0)] = 1.
    if ne:
        p[8, -1] = 1.
    jump, end = p[5]-p[11], p[8]
    spring = np.zeros((p.shape[1],)*2)
    for k, v in ((case.k_joint, jump), (case.k_end, end)):
        if k != RIGID:
            spring += k*m.section.EI/m.l*np.outer(v, v)
    return p, spring


def arm_dtn(m, length, omega):
    b = basis(m, length, omega, [0., length])
    q = np.vstack((b[0, :3], b[1, :3]))
    force = np.vstack((-b[0, 3:], b[1, 3:]))
    return np.linalg.solve(q.T, force.T).T, float(np.linalg.cond(balance(q)[0]))


def count_at(m, case, omega, catalogs, policy):
    poles = []
    for L in case.lengths(m):
        cat = catalogs[L]
        if omega >= cat["upper"]:
            raise ArithmeticError("Query outside certified arm-pole interval")
        poles += cat["poles"]
        d = m.section.EA-m.J*cat["upper"]**2
        if d <= 0:
            raise ArithmeticError("Outside RL acoustic interval")
        count = math.floor(cat["upper"]*math.sqrt(m.section.rhoA/d)*L/math.pi)
        poles += [math.sqrt(m.section.EA*(n*math.pi/L)**2 /
                  (m.section.rhoA+m.J*(n*math.pi/L)**2)) for n in range(1, count+1)]
    w = float(omega)
    exclusion = policy["pole_exclusion"]*max(1., w)
    for pole in sorted(poles):
        if abs(w-pole) <= exclusion:
            w = pole+2*exclusion
    pairs = [arm_dtn(m, L, w) for L in case.lengths(m)]
    p, spring = trace_map(m, case)
    a = p.T@block_diag(*(r[0] for r in pairs))@p+spring
    # A positive congruence only; never change the physical characteristic.
    scale = np.sqrt(m.l/m.section.EI)*np.array([m.l, m.l]+[1.]*(len(a)-2))
    a = a*scale[:, None]*scale
    skew = np.linalg.norm(a-a.T)/np.linalg.norm(a)
    vals = np.linalg.eigvalsh((a+a.T)/2)
    margin = min(abs(vals))/max(abs(vals))
    if skew > policy["schur_tol"] or margin < policy["inertia_margin"] or max(r[1] for r in pairs) > policy["condition_max"]:
        raise ArithmeticError(f"Unresolved DtN inertia: skew={skew}, margin={margin}")
    j0 = sum(pole < w for pole in poles)
    return w, j0+int(sum(vals < 0)), {"requested": omega, "effective": w,
        "J0": j0, "eigenvalues": vals.tolist(), "skew": float(skew),
        "margin": float(margin), "conditions": [r[1] for r in pairs]}


def prefix(m, case, catalogs, policy, neighbour=None):
    """Scoped adaptation of LONG count isolation to the five/four/three traces."""
    samples, records = [], []
    def query(w):
        if len(samples) >= policy["max_count_queries"]:
            raise ArithmeticError("Count budget exhausted")
        w, n, d = count_at(m, case, w, catalogs, policy)
        samples.append({"count": n, **d})
        return w, n
    def det(w):
        return float(np.linalg.det(balance(boundary(m, case, w))[0]))
    def isolate(a, b, ca, cb, depth=0):
        if ca >= 3 or cb == ca:
            return
        if cb < ca or depth > policy["max_depth"]:
            raise ArithmeticError("Unresolved count isolation")
        if cb-ca == 1 and det(a)*det(b) < 0:
            w, info = brentq(det, a, b, xtol=policy["root_xtol"],
                            rtol=policy["root_rtol"], full_output=True)
            records.append({"omega": w, "bracket": [a, b], "counts": [ca, cb],
                            "iterations": info.iterations})
            return
        mid, cm = query((a+b)/2)
        if not a < mid < b:
            raise ArithmeticError("Pole-shifted node outside bracket")
        isolate(a, mid, ca, cm, depth+1)
        isolate(mid, b, cm, cb, depth+1)
    a, ca = query(m.omega(policy["Lambda_lower"]))
    upper = m.omega(policy["Lambda_upper"])
    b, cb = query(upper)
    if ca != 0 or cb < 3:
        raise ArithmeticError("Declared bounds do not cover the three-root prefix")
    nodes = [(a, ca)]
    if neighbour is not None:
        for x in (np.mean(neighbour[:2]), np.mean(neighbour[1:3])):
            if a < x < b:
                nodes.append(query(float(x)))
    nodes.append((b, cb))
    for (lo, cl), (hi, ch) in zip(nodes[:-1], nodes[1:]):
        isolate(lo, hi, cl, ch)
    if len(records) != 3 or np.any(np.diff([r["omega"] for r in records]) <= policy["root_xtol"]):
        raise ArithmeticError("Incomplete, duplicate or multiple prefix")
    return records, {"status": "PASS", "samples": samples, "computed_roots": 3}


def mode(m, case, omega, order=120):
    raw = boundary(m, case, omega)
    b, _, col = balance(raw)
    _, singular, vh = np.linalg.svd(b)
    coef = (col*vh[-1]).reshape(2, 6)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    mass = energy = 0.
    arms = []
    s = m.section
    for L, c in zip(case.lengths(m), coef):
        x, wt = (nodes+1)*L/2, weights*L/2
        y, dy = (basis(m, L, omega, x, d)@c for d in (0, 1))
        ends = basis(m, L, omega, [0., L])@c
        mass += wt@(s.rhoA*(y[:, 0]**2+y[:, 1]**2)+s.rhoI*y[:, 2]**2+m.J*dy[:, 0]**2)
        energy += wt@(s.EA*dy[:, 0]**2+s.EI*dy[:, 2]**2+s.KGA*(dy[:, 1]-y[:, 2])**2)
        arms.append({"x": x, "weights": wt, "y": y, "dy": dy, "ends": ends})
    delta = arms[0]["ends"][1, 2]-arms[1]["ends"][1, 2]
    for k, v in ((case.k_joint, delta), (case.k_end, arms[1]["ends"][0, 2])):
        if k != RIGID:
            energy += k*s.EI/m.l*v*v
    for arm in arms:
        for key in ("y", "dy", "ends"):
            arm[key] /= math.sqrt(mass)
    coef /= math.sqrt(mass)
    qscale = max(np.max(abs(a["y"][:, :3]*[1, 1, m.l])) for a in arms)
    fscale = max(np.max(abs(a["y"][:, 3:]*[1, 1, 1/m.l])) for a in arms)
    jr = joint_matrix(m, case)@np.concatenate([a["ends"][1] for a in arms])
    js = abs(jr)/np.array([qscale, qscale, fscale, fscale,
                  qscale/m.l if case.k_joint == RIGID else fscale*m.l, fscale*m.l])
    er = [external_matrix(m, case, i)@a["ends"][0] for i, a in enumerate(arms)]
    es = [abs(r)/[qscale, qscale, qscale/m.l if i == 0 or case.k_end == RIGID else fscale*m.l]
          for i, r in enumerate(er)]
    scales = np.array([1, 1, m.l, qscale/fscale, qscale/fscale, qscale/(fscale*m.l)])
    ode = max(np.max(abs((a["dy"]-a["y"]@state_matrix(m, omega).T)*scales)) /
              max(np.max(abs(a["dy"]*scales)), 1e-30) for a in arms)
    return {"coefficients": coef.tolist(), "arms": arms, "diagnostics": {
        "mass": float(mass/math.sqrt(mass)**2), "raw_mass": float(mass),
        "energy_error": float(abs(energy/mass/omega**2-1)),
        "singular_ratio": float(singular[-1]/singular[0]),
        "nonzero_condition": float(singular[0]/singular[-2]),
        "joint_scaled": js.tolist(), "external_scaled": np.array(es).tolist(),
        "joint_absolute": jr.tolist(), "external_absolute": np.array(er).tolist(),
        "ode_error": float(ode), "delta_psi": float(delta/math.sqrt(mass))}}


def inner(m, first, second, beta1=None, beta2=None):
    value = 0.
    for i, (a, b) in enumerate(zip(first["arms"], second["arms"])):
        v, w = a["y"], b["y"]
        trans = v[:, :2]*w[:, :2]
        if beta1 is not None:
            f1, f2 = geometry.frames(beta1)[i], geometry.frames(beta2)[i]
            trans = (v[:, :2]@f1.translation.T)*(w[:, :2]@f2.translation.T)
        value += a["weights"]@(m.section.rhoA*trans.sum(axis=1)+m.section.rhoI*v[:, 2]*w[:, 2]+m.J*a["dy"][:, 0]*b["dy"][:, 0])
    return float(value)


def ritz(m, case, order, count=3):
    """Independent direct energies; reuse the existing clamped polynomial basis.

    The extra constant psi2 column is present ONLY for finite end springs.
    Force/moment equations are not imposed as Ritz constraints.
    """
    size = 6*order+int(case.k_end != RIGID)
    k, mass = np.zeros((size, size)), np.zeros((size, size))
    traces = []
    nodes, weights = np.polynomial.legendre.leggauss(max(64, order+2))
    b = shifted_legendre_clamped_basis(order, (nodes+1)/2)
    for i, L in enumerate(case.lengths(m)):
        v, d = np.zeros((len(nodes), 3, size)), np.zeros((len(nodes), 3, size))
        end, outer = np.zeros((3, size)), np.zeros((3, size))
        for f in range(3):
            sl = slice((3*i+f)*order, (3*i+f+1)*order)
            v[:, f, sl], d[:, f, sl] = b.values.T, b.derivatives.T/L
            end[f, sl] = 1.
        if i == 1 and case.k_end != RIGID:
            v[:, 2, -1] = end[2, -1] = outer[2, -1] = 1.
        wt = weights*L/2
        def gram(a, b):
            return a.T@(wt[:, None]*b)
        s = m.section
        k += s.EA*gram(d[:, 0], d[:, 0])+s.EI*gram(d[:, 2], d[:, 2])+s.KGA*gram(d[:, 1]-v[:, 2], d[:, 1]-v[:, 2])
        mass += s.rhoA*(gram(v[:, 0], v[:, 0])+gram(v[:, 1], v[:, 1]))+s.rhoI*gram(v[:, 2], v[:, 2])+m.J*gram(d[:, 0], d[:, 0])
        traces.append(end)
        if i == 1 and case.k_end != RIGID:
            k += case.k_end*s.EI/m.l*np.outer(outer[2], outer[2])
    g = [np.column_stack((f.t, f.n)) for f in geometry.frames(case.beta)]
    constraint = g[0]@traces[0][:2]-g[1]@traces[1][:2]
    jump = traces[0][2]-traces[1][2]
    if case.k_joint == RIGID:
        constraint = np.vstack((constraint, jump))
    else:
        k += case.k_joint*m.section.EI/m.l*np.outer(jump, jump)
    z = null_space(constraint)
    kr, mr = z.T@k@z, z.T@mass@z
    lk = cholesky((kr+kr.T)/2, lower=True)
    cholesky((mr+mr.T)/2, lower=True)
    a = solve_triangular(lk, mr, lower=True)
    a = solve_triangular(lk, a.T, lower=True).T
    val, v = eigh((a+a.T)/2, subset_by_index=[len(a)-count, len(a)-1])
    val, v = val[::-1], v[:, ::-1]
    coef = z@solve_triangular(lk.T, v, lower=False)/np.sqrt(val)
    return np.sqrt(1/val), coef, {"mass_gram_error": float(np.max(abs(coef.T@mass@coef-np.eye(count)))),
        "energy_error": float(np.max(abs(np.diag(coef.T@k@coef)*val-1))),
        "constraint_residual": float(np.max(abs(constraint@coef))), "positive_K_M": True}
