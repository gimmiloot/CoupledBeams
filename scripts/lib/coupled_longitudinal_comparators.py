"""Bounded elementary/planar Love frame comparators; shared verified Tim basis.

No c/R coordinates, no production MH changes. Local state (u,w,theta,N,Q,M).
Love boundary variation gives N=EA*u_x+J*u_xtt=(EA-J*omega**2)*U_x.
The independent count is Dirichlet poles plus the energy Schur index.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.optimize import brentq

from scripts.lib import bishop_longitudinal as bishop
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint

VERSION = "planar-reduced-axial-shared-tim-schur-v1"
NAMES = ("elementary", "rayleigh_love_planar")
STATE_ORDER = ("u", "w", "theta", "N", "Q", "M")
MIRROR = np.array([1., -1., -1., 1., -1., -1.])


def segment(model, length, variant):
    if variant not in NAMES:
        raise ValueError("Explicit elementary or planar Love variant required")
    s = model.section
    return bishop.Segment(length, s.EA, s.rhoA, H=0.,
        J=0. if variant == "elementary" else s.nu**2*s.rhoI)


def transform(frame):
    return frame.nodal_transform[np.ix_((0, 2, 3), (0, 2, 3))]


def joint_matrix(frames):
    g1, g2 = map(transform, frames)
    z = np.zeros((3, 3))
    return np.block([[g1, z, -g2, z], [z, g1, z, g2]])


def basis(model, length, omega, x, variant, derivative=0):
    """Six-by-six analytic basis; no copied Timoshenko formulas."""
    s = segment(model, length, variant)
    points = np.atleast_1d(x)
    result = np.zeros((len(points), 6, 6))
    result[:, 0, :2] = bishop.basis(s, omega, points, derivative)
    result[:, 3, :2] = (s.EA-s.J*omega**2)*bishop.basis(s, omega, points, derivative+1)
    result[:, (1, 2, 4, 5), 2:] = mh.finite_state_basis(model, length, omega,
        points, "timoshenko", derivative)
    return result[0] if np.ndim(x) == 0 else result


def state(model, length, omega, coefficients, x, variant, derivative=0):
    return basis(model, length, omega, x, variant, derivative)@np.asarray(coefficients)


def boundary(model, lengths, omega, variant, frames):
    matrix = np.zeros((12, 12))
    ends = np.zeros((12, 12))
    for arm, length in enumerate(lengths):
        cols = slice(6*arm, 6*(arm+1))
        matrix[3*arm:3*(arm+1), cols] = basis(model, length, omega, 0., variant)[:3]
        ends[6*arm:6*(arm+1), cols] = basis(model, length, omega, length, variant)
    matrix[6:] = joint_matrix(frames)@ends
    return matrix


def scaled_boundary(model, lengths, omega, variant, frames):
    raw = boundary(model, lengths, omega, variant, frames)
    return raw/np.linalg.norm(raw, axis=1)[:, None]


def arm_stiffness(model, length, omega, variant):
    values = basis(model, length, omega, [0., length], variant)
    matrix = np.concatenate((values[0, :3], values[1, :3]))
    norms = np.linalg.norm(matrix, axis=1)
    coefficients = np.linalg.solve(matrix/norms[:, None],
        np.vstack((np.zeros((3, 3)), np.eye(3)))/norms[:, None])
    return values[1, 3:]@coefficients, float(np.linalg.cond(matrix/norms[:, None]))


def axial_poles(model, length, variant, upper):
    s = segment(model, length, variant)
    if s.EA-s.J*upper**2 <= 0:
        raise ArithmeticError("Planar Love accumulation frequency exceeded")
    number = int(math.floor(upper*math.sqrt(s.m/(s.EA-s.J*upper**2))*length/math.pi))
    return [math.sqrt(s.EA*(n*math.pi/length)**2/(s.m+s.J*(n*math.pi/length)**2))
            for n in range(1, number+1)]


def count_at(model, lengths, omega, variant, frames, tim_catalogs, policy):
    """Off-pole count; query nodes shift, physical determinant roots do not.

    Coincident global root/arm pole may change count. Never assert equality
    of counts on the two sides. Return actual effective query frequency.
    """
    poles = []
    for length in lengths:
        catalog = tim_catalogs[length]
        if omega >= catalog["upper"]:
            raise ArithmeticError("Outside certified Timoshenko pole coverage")
        poles.extend(catalog["poles"])
        poles.extend(axial_poles(model, length, variant, catalog["upper"]))
    w = float(omega)
    exclusion = policy["pole_exclusion_relative"]*max(1., w)
    shifts = []
    for pole in sorted(poles):
        if abs(w-pole) <= exclusion:
            shifts.append(pole)
            w = pole+2*exclusion
    schur = np.zeros((3, 3))
    conditions = []
    for length, frame in zip(lengths, frames):
        local, condition = arm_stiffness(model, length, w, variant)
        g = transform(frame)
        schur += g@local@g.T
        conditions.append(condition)
    s, length = model.section, sum(lengths)
    scales = np.sqrt([length/s.EA, length**3/s.EI, length/s.EI])
    balanced = schur*scales[:, None]*scales[None, :]
    skew = float(np.max(np.abs(balanced-balanced.T))/max(np.linalg.norm(balanced, 2), 1e-30))
    eigenvalues = np.linalg.eigvalsh((balanced+balanced.T)/2)
    margin = float(min(abs(eigenvalues))/max(abs(eigenvalues)))
    if skew > policy["schur_symmetry_tol"] or margin < policy["count_inertia_margin_min"] or max(conditions) > policy["pole_condition_max"]:
        raise ArithmeticError(f"Comparator energy count conditioning unresolved: {skew}, {margin}, {conditions}")
    j0 = sum(pole < w for pole in poles)
    negative = int(np.count_nonzero(eigenvalues < 0))
    return w, j0+negative, {"requested_omega": omega, "effective_omega": w,
        "excluded_poles": shifts, "J0": j0, "negative_inertia": negative,
        "eigenvalues": eigenvalues.tolist(), "scaled_skew_residual": skew,
        "inertia_relative_margin": margin, "arm_boundary_conditions": conditions}


def roots(model, lengths, variant, frames, lower, upper, catalogs, policy):
    evaluations, subdivisions = 0, 0
    samples, brackets, failures = [], [], []
    def determinant(w):
        nonlocal evaluations
        evaluations += 1
        return float(np.linalg.det(scaled_boundary(model, lengths, w, variant, frames)))
    def query(w):
        effective, count, diagnostics = count_at(model, lengths, w, variant, frames, catalogs, policy)
        samples.append({"count": count, **diagnostics})
        return effective, count
    def interval(a, b, ca, cb, depth):
        nonlocal subdivisions
        if cb < ca:
            raise ArithmeticError("Nonmonotone comparator count")
        if cb == ca:
            return
        if cb-ca == 1:
            fa, fb = determinant(a), determinant(b)
            if fa*fb <= 0:
                brackets.append((a, b, fa, fb, ca, cb))
                return
        if depth >= policy["max_subdivision_depth"] or subdivisions >= policy["max_subdivisions"]:
            failures.append({"bracket": [a, b], "count_jump": cb-ca, "reason": "bounded subdivision exhausted"})
            return
        subdivisions += 1
        middle, cm = query((a+b)/2)
        if not a < middle < b:
            raise ArithmeticError("Pole exclusion leaves subdivision interval")
        interval(a, middle, ca, cm, depth+1)
        interval(middle, b, cm, cb, depth+1)
    nodes = [query(float(w)) for w in np.linspace(lower, upper, policy["scan_intervals"]+1)]
    if any(b[0] <= a[0] for a, b in zip(nodes[:-1], nodes[1:])):
        raise ArithmeticError("Pole-shifted scan nodes crossed")
    for (a, ca), (b, cb) in zip(nodes[:-1], nodes[1:]):
        interval(a, b, ca, cb, 0)
    records = []
    for a, b, fa, fb, ca, cb in brackets:
        omega, info = brentq(determinant, a, b, xtol=policy["root_xtol"], rtol=policy["root_rtol"], full_output=True)
        singular = np.linalg.svd(scaled_boundary(model, lengths, omega, variant, frames), compute_uv=False)
        records.append({"omega": omega, "frequency_hz": omega/(2*math.pi),
            "bracket_omega": [a, b], "bracket_determinants": [fa, fb], "bracket_counts": [ca, cb],
            "iterations": info.iterations, "converged": info.converged,
            "singular_ratio": float(singular[-1]/singular[0]),
            "nonzero_singular_condition": float(singular[0]/singular[-2])})
    search = {"range_omega": [nodes[0][0], nodes[-1][0]], "lower_count": nodes[0][1],
        "upper_count": nodes[-1][1], "scan_intervals": policy["scan_intervals"],
        "determinant_evaluations": evaluations, "subdivisions": subdivisions,
        "count_samples": samples, "failed_intervals": failures}
    if failures or len(records) != nodes[-1][1]-nodes[0][1] or nodes[0][1] != 0:
        raise ArithmeticError(f"Comparator root inventory incomplete: {search}")
    if any(b["omega"]-a["omega"] <= policy["root_xtol"] for a, b in zip(records[:-1], records[1:])):
        raise ArithmeticError("Duplicate or unresolved multiple root")
    search["status"] = "PASS"
    return records, search


def mode(model, lengths, omega, variant, frames, order):
    boundary_matrix = scaled_boundary(model, lengths, omega, variant, frames)
    _, singular, right = np.linalg.svd(boundary_matrix)
    coefficients = right[-1].reshape(2, 6)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    mass, energy, qpeak, ppeak, clamp, equation = 0., 0., 0., 0., 0., 0.
    endpoints = []
    p, total = model.coefficients, sum(lengths)
    for length, a in zip(lengths, coefficients):
        x, weight = (nodes+1)*length/2, weights*length/2
        value = state(model, length, omega, a, x, variant)
        gradient = state(model, length, omega, a, x, variant, 1)
        s = segment(model, length, variant)
        mass += float(weight@(s.m*(value[:, 0]**2+value[:, 1]**2)+p["r"]*value[:, 2]**2+s.J*gradient[:, 0]**2))
        energy += float(weight@(s.EA*gradient[:, 0]**2+p["B"]*gradient[:, 2]**2+p["S"]*(gradient[:, 1]-value[:, 2])**2))
        end = state(model, length, omega, a, [0., length], variant)
        qpeak = max(qpeak, float(np.max(abs(value[:, :3]*[1., 1., total]))))
        ppeak = max(ppeak, float(np.max(abs(value[:, 3:]*[1., 1., 1/total]))))
        clamp = max(clamp, float(np.max(abs(end[0, :3]*[1., 1., total]))))
        endpoints.append(end[1])
        operator = np.zeros((6, 6))
        operator[0, 3] = 1/(s.EA-s.J*omega**2)
        operator[3, 0] = -s.m*omega**2
        operator[np.ix_((1, 2, 4, 5), (1, 2, 4, 5))] = mh.harmonic_state_matrix(model, omega, "timoshenko")
        rhs = value@operator.T
        denominator = np.maximum(np.max(abs(rhs), axis=0)+np.max(abs(gradient), axis=0), 1e-30)
        equation = max(equation, float(np.max(abs(gradient-rhs)/denominator)))
    residual = joint_matrix(frames)@np.concatenate(endpoints)
    scaled = (
        abs(residual)*np.array([1., 1., total, 1., 1., 1/total])/np.array([qpeak]*3+[ppeak]*3))
    return {"coefficients": (coefficients/math.sqrt(mass)).tolist(), "diagnostics": {
        "mass_norm": mass/math.sqrt(mass)**2, "energy_relative_error": abs(energy/mass/omega**2-1),
        "clamp_scaled_residual": clamp/qpeak, "equation_scaled_residual": equation,
        "joint_residual_scaled": dict(zip(("d_X", "d_Y", "theta", "F_X", "F_Y", "M_node"), scaled.tolist())),
        "singular_ratio": float(singular[-1]/singular[0]),
        "nonzero_singular_condition": float(singular[0]/singular[-2])}}


def validate_mode(record, policy):
    d = record["diagnostics"]
    if (max(d["joint_residual_scaled"].values()) > policy["boundary_joint_scaled_tol"] or
        d["clamp_scaled_residual"] > policy["boundary_joint_scaled_tol"] or
        d["equation_scaled_residual"] > policy["equation_scaled_tol"] or
        d["energy_relative_error"] > policy["energy_relative_tol"] or
        d["nonzero_singular_condition"] > policy["nonzero_singular_condition_max"] or
        d["singular_ratio"] > policy["boundary_joint_scaled_tol"]):
        raise ArithmeticError(f"Comparator residual gate failed: {d}")


def overlap(first, second, weights, first_zero=None, second_zero=None):
    """Separate geometric metric; null pairs are None, never fake overlaps."""
    first, second = np.asarray(first), np.asarray(second)
    weights = np.asarray(weights)
    norms1 = np.sqrt(np.sum(first**2*weights, axis=-1))
    norms2 = np.sqrt(np.sum(second**2*weights, axis=-1))
    products = (first*weights)@second.T
    zeros1 = np.zeros(len(first), bool) if first_zero is None else first_zero
    zeros2 = np.zeros(len(second), bool) if second_zero is None else second_zero
    result = []
    for i in range(len(first)):
        result.append([None if zeros1[i] or zeros2[j] or norms1[i]*norms2[j] == 0 else
            min(1., float((products[i, j]/(norms1[i]*norms2[j]))**2)) for j in range(len(second))])
    return result, norms1.tolist(), norms2.tolist()


def contraction_diagnostic(c, nu_ux, weights, zero_scale=0.):
    norms = [math.sqrt(float(np.asarray(weights)@np.asarray(field)**2)) for field in (c, nu_ux, np.asarray(c)+nu_ux)]
    denominator = norms[0]+norms[1]
    return {"c_norm": norms[0], "nu_u_x_norm": norms[1], "sum_norm": norms[2],
        "numerical_zero_scale": zero_scale, "D_c": norms[2]/denominator if denominator > zero_scale else None,
        "status": "DEFINED" if denominator > zero_scale else "NOT_DEFINED_SMALL_FIELD"}
