"""Published reduced frame closure and common general-angle finite assembly.

Rucka (14),(15),(26)--(28): contraction/rotation are scalar common nodal
DOFs. Their compatibility is adopted reduced 1D closure, not a 3D joint
elasticity theorem. Arm physics comes only from the unchanged single-rod
module. Both positive local coordinates point from outer clamp to joint.
"""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

from scripts.lib import mindlin_herrmann_longitudinal as mh

VERSION = "jang-reduced-frame-beta0-transparency-v1"
GENERAL_VERSION = "jang-reduced-frame-general-angle-schur-count-v1"
STATE_ORDER = ("u", "c", "w", "theta", "N", "R", "Q", "M")
JOINT_ROWS = ("d_X", "d_Y", "c", "theta", "F_X", "F_Y", "R_node", "M_node")
BLOCK_INDICES = {"mh": (0, 1, 4, 5), "timoshenko": (2, 3, 6, 7)}
# Reflection maps local arm 2 fields to the single global axial coordinate.
REFLECTION = np.array([-1., 1., -1., 1., 1., -1., 1., -1.])


@dataclass(frozen=True)
class Frame:
    """Physical t,n in drawing (EX right,EY up); positive theta about -EZ.

    t x n = -EZ. Ux=u-z*theta makes theta a clockwise section rotation.
    c and theta are invariant under proper in-plane element rotations.
    This is a physical mapping, separate from legacy negative-x storage.
    """
    t: tuple[float, float]
    n: tuple[float, float]

    def __post_init__(self):
        a = self.translation
        if (not np.all(np.isfinite(a)) or
            not np.allclose(a.T@a, np.eye(2), rtol=0, atol=1e-14) or
            not math.isclose(float(np.linalg.det(a)), -1., rel_tol=0, abs_tol=1e-14)):
            raise ValueError("Frame must be orthonormal with t cross n = -EZ")

    @property
    def translation(self):
        return np.column_stack((self.t, self.n))

    @property
    def nodal_transform(self):
        """q_global=(dX,c,dY,theta)=T q_local; conjugate forces use T too."""
        matrix = np.eye(4)
        matrix[np.ix_((0, 2), (0, 2))] = self.translation
        return matrix


def beta0_frames():
    return frames(0.)


def frames(beta_deg):
    """Project geometry, not a new angle convention or Reddy arm physics.

    Positive beta turns joint->right-clamp ray upward. Local x points from
    each clamp to joint. Only the existing geometry helper's t,n are used;
    Jang theta's physical sign remains about k=-EZ, as in the beta0 audit.
    """
    from scripts.lib.reddy_inplane_geometry import reddy_inplane_geometry
    geometry = reddy_inplane_geometry(beta_deg)
    return tuple(Frame(tuple(arm.t[:2]), tuple(arm.n[:2])) for arm in (geometry.arm1, geometry.arm2))


MIRROR_STATE = np.array([1., 1., -1., -1., 1., 1., -1., -1.])


def reflected_frames(original):
    """Reflect EX/EY geometry in EX; keep t cross n=-EZ via n*=-S n.

    Local state transforms by MIRROR_STATE: theta/M change sign, c/R do
    not. c is a scalar, theta is a signed planar rotation (pseudoscalar).
    """
    mirror = np.diag([1., -1.])
    return tuple(Frame(tuple(mirror@f.t), tuple(-mirror@f.n)) for f in original)


def endpoint_sign(end):
    """[p delta q]_0^L: outward nodal efforts are -p at 0, +p at L."""
    if end not in ("left", "right"):
        raise ValueError("Endpoint must be left or right")
    return -1. if end == "left" else 1.


def nodal_efforts(state, frame, end):
    return endpoint_sign(end)*frame.nodal_transform@np.asarray(state)[4:]


def joint_residual(left, right, frames=None, ends=("right", "right")):
    """Eight invariant rows: d1-d2,c1-c2,theta1-theta2; sum outward efforts."""
    f1, f2 = beta0_frames() if frames is None else frames
    q1, q2 = f1.nodal_transform@np.asarray(left)[:4], f2.nodal_transform@np.asarray(right)[:4]
    difference = (q1-q2)[[0, 2, 1, 3]]
    balance = (nodal_efforts(left, f1, ends[0])+nodal_efforts(right, f2, ends[1]))[[0, 2, 1, 3]]
    return np.concatenate((difference, balance))


def joint_matrix(frames=None, ends=("right", "right")):
    identity = np.eye(16)
    return np.column_stack([joint_residual(v[:8], v[8:], frames, ends) for v in identity])


def virtual_work(states, variation, frames=None, ends=("right", "right")):
    """Arbitrary compatible variation, without first imposing equilibrium."""
    frames = beta0_frames() if frames is None else frames
    variation = np.asarray(variation)  # global order dX,c,dY,theta
    local_work = sum(endpoint_sign(end)*np.dot(np.asarray(state)[4:], frame.nodal_transform.T@variation)
                     for state, frame, end in zip(states, frames, ends))
    global_work = np.dot(sum(nodal_efforts(s, f, e) for s, f, e in zip(states, frames, ends)), variation)
    return float(local_work), float(global_work)


def _lengths(length, split, beta_deg=0.):
    if beta_deg != 0:
        raise ValueError("Only beta=0 finite assembly is authorized and implemented")
    if not math.isfinite(length) or length <= 0 or not math.isfinite(split) or not 0 < split < 1:
        raise ValueError("Positive total length and internal split required")
    return length*split, length*(1-split)


def arm_basis(model, length, omega, x):
    """q-first state/coefficients (MH4,Tim4), reusing bounded arm basis."""
    result = np.zeros((8, 8))
    for block, cols in (("mh", range(4)), ("timoshenko", range(4, 8))):
        result[np.ix_(BLOCK_INDICES[block], cols)] = mh.finite_state_basis(model, length, omega, x, block)
    return result


def boundary_matrix(model, length, split, omega, *, beta_deg=0.):
    """Frozen beta0 API guard; delegates to the ONE general assembly."""
    l1, l2 = _lengths(length, split, beta_deg)
    return frame_boundary_matrix(model, (l1, l2), omega, beta_deg=beta_deg)


def frame_boundary_matrix(model, lengths, omega, *, beta_deg=0., arm_frames=None):
    """16x16: two four-field clamps and eight invariant joint rows.

    Explicit frames are for rigid-body remapping/swap/reflection audits;
    they use exactly the same joint operator and arm coefficients.
    """
    if len(lengths) != 2 or any(not math.isfinite(l) or l <= 0 for l in lengths):
        raise ValueError("Two positive arm lengths required")
    l1, l2 = lengths
    active_frames = frames(beta_deg) if arm_frames is None else arm_frames
    left, right = arm_basis(model, l1, omega, 0.), arm_basis(model, l2, omega, 0.)
    matrix = np.zeros((16, 16))
    matrix[:4, :8], matrix[4:8, 8:] = left[:4], right[:4]
    joined = np.zeros((16, 16))
    joined[:8, :8] = arm_basis(model, l1, omega, l1)
    joined[8:, 8:] = arm_basis(model, l2, omega, l2)
    matrix[8:] = joint_matrix(active_frames)@joined
    return matrix


def scaled_frame_matrix(model, lengths, omega, *, beta_deg=0., arm_frames=None):
    raw = frame_boundary_matrix(model, lengths, omega, beta_deg=beta_deg, arm_frames=arm_frames)
    row_norms = np.linalg.norm(raw, axis=1)
    if np.any(row_norms == 0) or not np.all(np.isfinite(raw)):
        raise ArithmeticError("Invalid frame boundary matrix")
    return raw/row_norms[:, None]


def arm_dynamic_stiffness(model, length, omega):
    """Exact q0=0 end Dirichlet-to-Neumann map, away from CC poles.

    Independent of the 16x16 joint determinant; energy gives symmetric
    4x4 K on local end q=(u,c,w,theta), nodal p=(N,R,Q,M).
    """
    stiffness = np.zeros((4, 4))
    condition = 0.
    for block, indices in (("mh", (0, 1)), ("timoshenko", (2, 3))):
        basis = mh.finite_state_basis(model, length, omega, [0., length], block)
        boundary = np.concatenate((basis[0, :2], basis[1, :2]))
        scales = np.array([1., length, 1., length])
        balanced = boundary*scales[:, None]
        norms = np.linalg.norm(balanced, axis=1)
        balanced /= norms[:, None]
        rhs = np.vstack((np.zeros((2, 2)), np.eye(2)))*scales[:, None]/norms[:, None]
        coefficients = np.linalg.solve(balanced, rhs)
        stiffness[np.ix_(indices, indices)] = basis[1, 2:]@coefficients
        condition = max(condition, float(np.linalg.cond(balanced)))
    return stiffness, condition


def nodal_schur_matrix(model, lengths, omega, *, beta_deg=0., arm_frames=None):
    active_frames = frames(beta_deg) if arm_frames is None else arm_frames
    matrix, conditions = np.zeros((4, 4)), []
    for length, frame in zip(lengths, active_frames):
        local, condition = arm_dynamic_stiffness(model, length, omega)
        transform = frame.nodal_transform
        matrix += transform@local@transform.T
        conditions.append(condition)
    # Fixed positive congruence: balances energy, preserves inertia exactly.
    p, length = model.coefficients, sum(lengths)
    scales = np.array([math.sqrt(length/p["C"]), (p["C"]*p["H"])**(-.25),
                       math.sqrt(length**3/p["B"]), math.sqrt(length/p["B"])])
    balanced = matrix*scales[:, None]*scales[None, :]
    skew = float(np.max(np.abs(balanced-balanced.T)))/max(float(np.linalg.norm(balanced, 2)), 1e-30)
    # Roundoff-only symmetric part for a theoretically symmetric energy map.
    symmetric = (balanced+balanced.T)/2
    eigenvalues = np.linalg.eigvalsh(symmetric)
    return symmetric, {"eigenvalues": eigenvalues.tolist(), "scaled_skew_residual": skew,
        "negative_inertia": int(np.count_nonzero(eigenvalues < 0)),
        "arm_boundary_conditions": conditions,
        "inertia_relative_margin": float(np.min(np.abs(eigenvalues))/np.max(np.abs(eigenvalues)))}


def frame_roots(model, lengths, beta_deg, lower, upper, policy, count, arm_frames=None):
    """Count-certified bounded intervals + primary 16x16 determinant roots.

    Count callback is the separately derived energy/Schur index, not this
    determinant. One bounded subdivision mechanism resolves missed/close
    simple roots. Unresolved multiplicity/conditioning stops the pilot.
    """
    from scipy.optimize import brentq
    evaluations, subdivisions = 0, 0
    failures, brackets, count_samples = [], [], []
    def determinant(w):
        nonlocal evaluations
        evaluations += 1
        return float(np.linalg.det(scaled_frame_matrix(model, lengths, w,
            beta_deg=beta_deg, arm_frames=arm_frames)))
    def certified(w):
        value, diagnostic = count(w, beta_deg, arm_frames)
        count_samples.append({"omega": w, "count": value, **diagnostic})
        return value
    def interval(a, b, ca, cb, depth):
        nonlocal subdivisions
        number = cb-ca
        if number < 0:
            raise ArithmeticError("Energy count is not monotone")
        if number == 0:
            return
        if number == 1:
            fa, fb = determinant(a), determinant(b)
            if fa*fb <= 0:
                brackets.append((a, b, fa, fb, ca, cb))
                return
        if depth >= policy["max_subdivision_depth"] or subdivisions >= policy["max_subdivisions"]:
            failures.append({"range_omega": [a, b], "count_jump": number,
                "reason": "unresolved cluster/multiplicity or sign/conditioning within fixed budget"})
            return
        subdivisions += 1
        middle = (a+b)/2
        cm = certified(middle)
        interval(a, middle, ca, cm, depth+1)
        interval(middle, b, cm, cb, depth+1)
    nodes = np.linspace(lower, upper, policy["scan_intervals"]+1)
    counts = [certified(float(w)) for w in nodes]
    for a, b, ca, cb in zip(nodes[:-1], nodes[1:], counts[:-1], counts[1:]):
        interval(float(a), float(b), ca, cb, 0)
    records = []
    for a, b, fa, fb, ca, cb in brackets:
        omega, info = brentq(determinant, a, b, xtol=policy["root_xtol"],
                            rtol=policy["root_rtol"], full_output=True)
        singular = np.linalg.svd(scaled_frame_matrix(model, lengths, omega,
            beta_deg=beta_deg, arm_frames=arm_frames), compute_uv=False)
        records.append({"omega": omega, "frequency_hz": omega/(2*math.pi),
            "bracket_omega": [a, b], "bracket_determinants": [fa, fb],
            "bracket_counts": [ca, cb], "iterations": info.iterations, "converged": info.converged,
            "singular_ratio": float(singular[-1]/singular[0]),
            "nonzero_singular_condition": float(singular[0]/singular[-2])})
    if failures or len(records) != counts[-1]-counts[0]:
        raise ArithmeticError(f"Frame root count unresolved: found={len(records)}, counts={counts[0],counts[-1]}, failures={failures}")
    if any(b["omega"]-a["omega"] <= policy["root_xtol"] for a, b in zip(records[:-1], records[1:])):
        raise ArithmeticError("Duplicated or unresolved multiple roots")
    return records, {"range_omega": [lower, upper], "scan_intervals": policy["scan_intervals"],
        "determinant_evaluations": evaluations, "subdivisions": subdivisions,
        "lower_count": counts[0], "upper_count": counts[-1], "count_samples": count_samples,
        "failed_intervals": failures, "status": "PASS"}


def arm_state(model, length, omega, coefficients, x, derivative=0):
    """Full local state/gradient, unchanged independent local arm blocks."""
    points = np.atleast_1d(np.asarray(x, dtype=float))
    result = np.zeros((len(points), 8))
    for block, cols in (("mh", slice(0, 4)), ("timoshenko", slice(4, 8))):
        result[:, BLOCK_INDICES[block]] = mh.finite_state_basis(model, length, omega,
            points, block, derivative)@np.asarray(coefficients)[cols]
    return result[0] if np.ndim(x) == 0 else result


def frame_mode(model, lengths, beta_deg, omega, order=200, arm_frames=None):
    """Mass-normalized full mode and dimensionally scaled residual gates."""
    active_frames = frames(beta_deg) if arm_frames is None else arm_frames
    boundary = scaled_frame_matrix(model, lengths, omega, beta_deg=beta_deg, arm_frames=arm_frames)
    _, singular, right = np.linalg.svd(boundary)
    coefficients = right[-1].reshape(2, 8)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    p, mass, energy = model.coefficients, 0., 0.
    mass_density = np.array([p["m"], p["j"], p["m"], p["r"]])
    qpeak, fpeak, clamp, ode = 0., 0., 0., 0.
    cpeak, rpeak = 0., 0.
    endpoints, local_work = [], []
    for length, a in zip(lengths, coefficients):
        x, weight = (nodes+1)*length/2, weights*length/2
        value = arm_state(model, length, omega, a, x)
        gradient = arm_state(model, length, omega, a, x, 1)
        mass += float(weight@(value[:, :4]**2@mass_density))
        energy += float(weight@(p["C"]*(gradient[:, 0]**2+2*model.section.nu*gradient[:, 0]*value[:, 1]+value[:, 1]**2)
            +p["H"]*gradient[:, 1]**2+p["B"]*gradient[:, 3]**2+p["S"]*(gradient[:, 2]-value[:, 3])**2))
        end = arm_state(model, length, omega, a, [0., length])
        qpeak = max(qpeak, float(np.max(np.abs(value[:, :4]*[1., sum(lengths), 1., sum(lengths)]))))
        fpeak = max(fpeak, float(np.max(np.abs(value[:, 4:]*[1., 1/sum(lengths), 1., 1/sum(lengths)]))))
        cpeak = max(cpeak, float(np.max(np.abs(value[:, 1]))))
        rpeak = max(rpeak, float(np.max(np.abs(value[:, 5]))))
        clamp = max(clamp, float(np.max(np.abs(end[0, :4]*[1., sum(lengths), 1., sum(lengths)]))))
        endpoints.append(end[1])
        rhs = value@mh.full_harmonic_state_matrix(model, omega).T
        denominator = np.maximum(np.max(np.abs(rhs), axis=0)+np.max(np.abs(gradient), axis=0), 1e-30)
        ode = max(ode, float(np.max(np.abs(gradient-rhs)/denominator)))
        local_work.append(float(end[1, :4]@end[1, 4:]))
    residual = joint_residual(*endpoints, active_frames)
    length = sum(lengths)
    scaled = (
        np.abs(residual)*np.array([1., 1., length, length, 1., 1., 1/length, 1/length]) /
        np.array([qpeak]*4+[fpeak]*4))
    return {"coefficients": coefficients/math.sqrt(mass), "diagnostics": {
        "mass_norm": mass/math.sqrt(mass)**2, "mass_before_normalization": mass,
        "energy_relative_error": abs(energy/mass/omega**2-1),
        "clamp_scaled_residual": clamp/qpeak, "equation_scaled_residual": ode,
        "joint_residual_raw": dict(zip(JOINT_ROWS, (residual/math.sqrt(mass)).tolist())),
        "joint_residual_scaled": dict(zip(JOINT_ROWS, scaled.tolist())),
        "c_compatibility_own_amplitude": abs(float(residual[2]))/cpeak if cpeak else 0.,
        "R_balance_own_amplitude": abs(float(residual[6]))/rpeak if rpeak else 0.,
        "joint_work_relative": abs(sum(local_work))/max(energy, 1e-30),
        "singular_ratio": float(singular[-1]/singular[0]),
        "nonzero_singular_condition": float(singular[0]/singular[-2])}}


def block_matrix(model, length, split, omega, block):
    full = boundary_matrix(model, length, split, omega)
    if block == "mh":
        rows, cols = (0, 1, 4, 5, 8, 10, 12, 14), (*range(4), *range(8, 12))
    elif block == "timoshenko":
        rows, cols = (2, 3, 6, 7, 9, 11, 13, 15), (*range(4, 8), *range(12, 16))
    else:
        raise ValueError("Unknown block")
    return full[np.ix_(rows, cols)]


def scaled_block_matrix(model, length, split, omega, block):
    raw = block_matrix(model, length, split, omega, block)
    norms = np.linalg.norm(raw, axis=1)
    if np.any(norms == 0) or not np.all(np.isfinite(norms)):
        raise ArithmeticError("Invalid joint boundary scaling")
    return raw/norms[:, None]  # strictly positive, no signed fitting


def roots(model, length, split, block, lower, upper, policy):
    """Bounded independent two-arm sign scan; min-max count checked by CLI."""
    from scipy.optimize import brentq
    evaluations = 0
    def determinant(omega):
        nonlocal evaluations
        evaluations += 1
        return float(np.linalg.det(scaled_block_matrix(model, length, split, omega, block)))
    nodes = np.linspace(lower, upper, policy["scan_intervals"]+1)
    samples = [determinant(w) for w in nodes]
    records, failures = [], []
    for a, b, fa, fb in zip(nodes[:-1], nodes[1:], samples[:-1], samples[1:]):
        if fa*fb > 0:
            continue
        before = evaluations
        try:
            omega, info = brentq(determinant, a, b, xtol=policy["root_xtol"],
                                rtol=policy["root_rtol"], full_output=True)
        except (ValueError, RuntimeError) as exc:
            failures.append({"bracket_omega": [float(a), float(b)], "reason": str(exc)})
            continue
        if records and abs(omega-records[-1]["omega"]) <= 10*policy["root_xtol"]:
            continue
        singular = np.linalg.svd(scaled_block_matrix(model, length, split, omega, block), compute_uv=False)
        records.append({"omega": omega, "frequency_hz": omega/(2*math.pi),
            "bracket_omega": [float(a), float(b)], "bracket_determinants": [fa, fb],
            "evaluations": evaluations-before, "iterations": info.iterations,
            "converged": info.converged, "singular_ratio": float(singular[-1]/singular[0]),
            "nonzero_singular_condition": float(singular[0]/singular[-2])})
    return records, {"range_omega": [lower, upper], "scan_intervals": policy["scan_intervals"],
        "evaluations": evaluations, "brackets_found": len(records), "failed_intervals": failures}


def mode_coefficients(model, length, split, omega, block, order=200):
    matrix = scaled_block_matrix(model, length, split, omega, block)
    _, _, right = np.linalg.svd(matrix)
    coefficients = right[-1].reshape(2, 4)
    points, weights = quadrature(length, split, order)
    values = mode_state(model, length, split, omega, block, coefficients, points)
    p = model.coefficients
    mass = float(weights@(p["m"]*values[:, 0]**2+(p["j"] if block == "mh" else p["r"])*values[:, 1]**2))
    return coefficients/math.sqrt(mass)


def quadrature(length, split, order):
    l1, l2 = _lengths(length, split)
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return (np.concatenate(((nodes+1)*l1/2, l1+(nodes+1)*l2/2)),
            np.concatenate((weights*l1/2, weights*l2/2)))


def mode_state(model, length, split, omega, block, coefficients, points):
    """Single-coordinate state on X in [0,L]. Right local x2=L-X.

    Each 4-state reflection is (-first,+second,+force,-second effort).
    Inactive block is absent, not numerical noise from a full nullvector.
    """
    l1, l2 = _lengths(length, split)
    points = np.atleast_1d(np.asarray(points, dtype=float))
    if np.any(points < 0) or np.any(points > length):
        raise ValueError("Points outside straight rod")
    values = np.empty((len(points), 4))
    left = points <= l1
    for mask, size, x, a in ((left, l1, points[left], coefficients[0]),
                             (~left, l2, length-points[~left], coefficients[1])):
        if len(x):
            values[mask] = mh.finite_state_basis(model, size, omega, np.clip(x, 0., size), block)@a
    values[~left] *= np.array([-1., 1., 1., -1.])
    return values


def full_state(values, block):
    result = np.zeros((len(values), 8))
    result[:, BLOCK_INDICES[block]] = values
    return result


def mode_diagnostics(model, length, split, omega, block, coefficients, reference_omega, reference_coefficients, order):
    points, weights = quadrature(length, split, order)
    values = mode_state(model, length, split, omega, block, coefficients, points)
    direct = mh.finite_state_basis(model, length, reference_omega, points, block)@reference_coefficients
    p = model.coefficients
    mass_weights = np.array([p["m"], p["j"] if block == "mh" else p["r"]])
    norm = float(weights@(values[:, :2]**2@mass_weights))
    refnorm = float(weights@(direct[:, :2]**2@mass_weights))
    overlap = float(weights@((values[:, :2]*direct[:, :2])@mass_weights))/math.sqrt(norm*refnorm)
    sign = 1. if overlap >= 0 else -1.
    values *= sign
    l2_errors = [math.sqrt(float(weights@(values[:, i]-direct[:, i])**2)/float(weights@direct[:, i]**2)) for i in range(2)]
    l1, l2 = _lengths(length, split)
    endpoints = [mh.finite_state_basis(model, size, omega, size, block)@a*sign
                 for size, a in zip((l1, l2), coefficients)]
    states = [full_state(v[None, :], block)[0] for v in endpoints]
    residual = joint_residual(*states)
    # Component amplitude scaling in global state, same dimensions throughout.
    peak = np.max(np.abs(values), axis=0)
    boundary = np.concatenate([mh.finite_state_basis(model, size, omega, 0., block)@a
                               for size, a in zip((l1, l2), coefficients)]).reshape(2, 4)
    clamp = float(np.max(np.abs(boundary[:, :2])/peak[:2]))
    if block == "mh":
        scales = np.array([peak[0], 1., peak[1], 1., peak[2], 1., peak[3], 1.])
    else:
        scales = np.array([1., peak[0], 1., peak[1], 1., peak[2], 1., peak[3]])
    scaled_residual = np.abs(residual)/scales
    # Compare each interface state's ordinary/resultant components independently.
    ref = mh.finite_state_basis(model, length, reference_omega, l1, block)@reference_coefficients
    interface_error = max(float(np.max(np.abs(endpoints[0]-ref)/peak)),
                          float(np.max(np.abs(endpoints[1]*[-1, 1, 1, -1]-ref)/peak)))
    return {"mass_overlap_abs": abs(overlap), "mass_MAC": min(1., overlap**2),
        "component_L2_relative": dict(zip(("u", "c") if block == "mh" else ("w", "theta"), l2_errors)),
        "mass_norm": norm, "clamp_scaled_residual": clamp,
        "joint_residual_raw": dict(zip(JOINT_ROWS, residual.tolist())),
        "joint_residual_scaled": dict(zip(JOINT_ROWS, scaled_residual.tolist())),
        "interface_reference_scaled_error": interface_error, "sign_alignment": sign}


def segmented_transfer(model, length, split, omega, block, exponent_cap=1., max_steps=512):
    """Independent propagation through positive global segments, short expm/QR.

    Tests the semigroup on an essential initial subspace; no unsafe product
    of full-length transfer matrices (MH alpha*L is large).
    """
    from scipy.linalg import expm
    lengths = _lengths(length, split)
    matrix = mh.harmonic_state_matrix(model, omega, block)
    rate = float(np.max(np.abs(np.linalg.eigvals(matrix))))
    p = model.coefficients
    elastic, gradient = (p["C"], p["H"]) if block == "mh" else (p["S"], p["B"])
    scales = np.array([1., length, 1/(elastic*rate), length/(gradient*rate)])
    balanced = matrix*scales[:, None]/scales[None, :]
    counts = [max(1, math.ceil(rate*d/exponent_cap)) for d in lengths]
    if sum(counts) > max_steps:
        raise ArithmeticError("Segmented QR propagation step budget exceeded")
    frame = np.vstack((np.zeros((2, 2)), np.eye(2)))
    conditions = []
    for d, count in zip(lengths, counts):
        step = expm(balanced*d/count)
        conditions.append(float(np.linalg.cond(step)))
        for _ in range(count):
            frame, triangular = np.linalg.qr(step@frame, mode="reduced")
            frame *= np.where(np.diag(triangular) >= 0, 1., -1.)[None, :]
    direct, _ = mh.transfer_boundary_matrix(model, length, omega, block, exponent_cap, max_steps)
    return {"steps_by_segment": counts, "max_step_condition": max(conditions),
        "direct_projected_boundary_difference": float(np.max(np.abs(frame[:2]-direct))),
        "projected_boundary_singular_ratio": float(np.linalg.svd(frame[:2], compute_uv=False)[-1]/np.linalg.norm(frame[:2], 2))}
