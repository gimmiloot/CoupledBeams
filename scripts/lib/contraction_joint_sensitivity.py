"""Isolated c-closure diagnostic, reusing unchanged production arm physics.

Theory/count proof: docs/theory/contraction_joint_sensitivity.md.
No production dispatch, no FEM, no parameter-map or branch-tracking API.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import block_diag
from scipy.optimize import brentq

from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as rigid
from scripts.lib import coupled_longitudinal_comparators as reduced


@dataclass(frozen=True)
class Case:
    beta: float
    closure: str = "CONTINUOUS"
    variant: str = "mindlin_herrmann"
    rotation: str = "SPRING"

    def __post_init__(self):
        if self.closure not in ("FREE", "CONTINUOUS", "BLOCKED"):
            raise ValueError("Explicit c closure required")
        if self.variant not in ("mindlin_herrmann", *reduced.NAMES):
            raise ValueError("Unknown axial theory")
        if self.rotation not in ("SPRING", "RIGID"):
            raise ValueError("Unknown rotation closure")

    @property
    def is_mh(self):
        return self.variant == "mindlin_herrmann"

    @property
    def nq(self):
        return 4 if self.is_mh else 3

    @property
    def key(self):
        return f"b{self.beta:g}_{self.variant}_{self.closure if self.is_mh else 'NA'}_{self.rotation}"


def basis(model, length, omega, x, case, derivative=0):
    if not case.is_mh:
        return reduced.basis(model, length, omega, x, case.variant, derivative)
    points = np.atleast_1d(x)
    values = np.zeros((len(points), 8, 8))
    for block, cols in (("mh", range(4)), ("timoshenko", range(4, 8))):
        values[:, rigid.BLOCK_INDICES[block], cols.start:cols.stop] = mh.finite_state_basis(
            model, length, omega, points, block, derivative)
    return values[0] if np.ndim(x) == 0 else values


def trace_map(case):
    """Local endpoint q=P z; z=(dX,dY,theta1[,theta2],c trace...)."""
    nr = 1 if case.rotation == "RIGID" else 2
    nc = {"FREE": 2, "CONTINUOUS": 1, "BLOCKED": 0}[case.closure] if case.is_mh else 0
    nq, nz = case.nq, 2+nr+nc
    matrix = np.zeros((2*nq, nz))
    for i, frame in enumerate(rigid.frames(case.beta)):
        u, w, theta = (0, 2, 3) if case.is_mh else (0, 1, 2)
        matrix[i*nq+u, :2] = frame.t
        matrix[i*nq+w, :2] = frame.n
        matrix[i*nq+theta, 2+(i if nr == 2 else 0)] = 1
        if nc:
            matrix[i*nq+1, 2+nr+(i if nc == 2 else 0)] = 1
    jump = matrix[nq-1]-matrix[2*nq-1]
    return matrix, jump


def joint_matrix(case, stiffness):
    """Rows dX,dY,FX,FY,rotation1,rotation2[,c row1,c row2].

    Both arm endpoints are local right ends, so both outward signs are +.
    Rigid mode uses delta theta=0 and M1+M2=0 exactly, never a large k.
    """
    n = case.nq
    operator = np.zeros((2*n, 4*n), dtype=np.result_type(stiffness, float))
    frames = rigid.frames(case.beta)
    u, w, theta = (0, 2, 3) if case.is_mh else (0, 1, 2)
    for i, frame in enumerate(frames):
        offset = 2*n*i
        sign = 1 if i == 0 else -1
        operator[:2, offset+u] = sign*np.array(frame.t)
        operator[:2, offset+w] = sign*np.array(frame.n)
        operator[2:4, offset+n+u] = frame.t
        operator[2:4, offset+n+w] = frame.n
    if case.rotation == "RIGID":
        operator[4, [theta, 2*n+theta]] = [1, -1]
        operator[5, [n+theta, 3*n+theta]] = 1
    else:
        operator[4, [theta, 2*n+theta, n+theta]] = [stiffness, -stiffness, 1]
        operator[5, [theta, 2*n+theta, 3*n+theta]] = [-stiffness, stiffness, 1]
    if case.is_mh:
        if case.closure == "FREE":
            operator[6, 5] = operator[7, 13] = 1
        elif case.closure == "CONTINUOUS":
            operator[6, [1, 9]] = [1, -1]
            operator[7, [5, 13]] = 1
        else:
            operator[6, 1] = operator[7, 9] = 1
    return operator


def assemble_boundary(endpoint_bases, case, stiffness):
    n = case.nq
    matrix = np.zeros((4*n, 4*n), dtype=np.result_type(*endpoint_bases, stiffness))
    ends = np.zeros_like(matrix)
    for i, values in enumerate(endpoint_bases):
        cols = slice(2*n*i, 2*n*(i+1))
        matrix[n*i:n*(i+1), cols] = values[0, :n]
        ends[cols, cols] = values[1]
    matrix[2*n:] = joint_matrix(case, stiffness)@ends
    return matrix


def boundary(model, lengths, omega, case, stiffness, scaled=True):
    raw = assemble_boundary([basis(model, l, omega, [0., l], case) for l in lengths], case, stiffness)
    return raw/np.linalg.norm(raw, axis=1)[:, None] if scaled else raw


def schur(model, lengths, omega, case, stiffness):
    """New admissible trace space: P.T diag(K_arm) P + k v v.T.

    The old common-theta/common-c 4x4 count is not used.
    """
    matrices, conditions = [], []
    for length in lengths:
        k, cond = (rigid.arm_dynamic_stiffness(model, length, omega) if case.is_mh else
                   reduced.arm_stiffness(model, length, omega, case.variant))
        matrices.append(k)
        conditions.append(cond)
    pmap, jump = trace_map(case)
    raw = pmap.T@block_diag(*matrices)@pmap+stiffness*np.outer(jump, jump)
    p, length = model.coefficients, sum(lengths)
    nr = 1 if case.rotation == "RIGID" else 2
    nc = pmap.shape[1]-2-nr
    scales = np.array([math.sqrt(length/p["C"]), math.sqrt(length**3/p["B"])]
        +[math.sqrt(length/p["B"])]*nr+[(p["C"]*p["H"])**(-.25)]*nc)
    balanced = raw*scales[:, None]*scales[None, :]
    skew = float(np.linalg.norm(balanced-balanced.T, 2)/max(np.linalg.norm(balanced, 2), 1e-30))
    eigenvalues = np.linalg.eigvalsh((balanced+balanced.T)/2)
    return raw, {"dimension": len(raw), "eigenvalues": eigenvalues.tolist(),
        "skew": skew, "arm_conditions": conditions,
        "margin": float(min(abs(eigenvalues))/max(abs(eigenvalues))),
        "negative_inertia": int(np.count_nonzero(eigenvalues < 0))}


def count(model, lengths, omega, case, stiffness, catalog, policy):
    poles = []
    for length in lengths:
        tim = catalog["timoshenko"]
        ceiling = tim["search"]["range_omega"][1]
        poles.extend(r["omega"] for r in tim["roots"])
        if case.is_mh:
            axial = catalog["mh"]
            ceiling = min(ceiling, axial["search"]["range_omega"][1])
            poles.extend(r["omega"] for r in axial["roots"])
        else:
            poles.extend(reduced.axial_poles(model, length, case.variant, ceiling))
        if omega >= ceiling or length != .5:
            raise ArithmeticError("Outside reused G20 half-arm pole catalog")
    effective = float(omega)
    exclusion = policy["pole_exclusion_relative"]*max(1., omega)
    for pole in sorted(poles):
        if abs(effective-pole) <= exclusion:
            effective = pole+2*exclusion
    _, diag = schur(model, lengths, effective, case, stiffness)
    if (diag["skew"] > policy["schur_symmetry_tol"] or
        diag["margin"] < policy["count_inertia_margin_min"] or
        max(diag["arm_conditions"]) > policy["pole_condition_max"]):
        raise ArithmeticError(f"Unresolved count conditioning: {diag}")
    j0 = sum(p < effective for p in poles)
    return effective, j0+diag["negative_inertia"], {"requested_omega": omega,
        "effective_omega": effective, "J0": j0, **diag}


def prefix_roots(model, lengths, case, stiffness, catalog, policy):
    """Count-guided prefix isolation plus scipy Brent on the full matrix.

    This small adapter is necessary for the changed trace spaces. It reuses
    arm DtN/bases and certified poles, does no global scan or root deletion.
    """
    samples, brackets, roots = [], [], []
    def query(w):
        if len(samples) >= policy["max_count_queries"]:
            raise ArithmeticError("Count-query budget exhausted")
        effective, number, diag = count(model, lengths, w, case, stiffness, catalog, policy)
        samples.append({"count": number, **diag})
        return effective, number
    def determinant(w):
        return float(np.linalg.det(boundary(model, lengths, w, case, stiffness)))
    def isolate(a, b, ca, cb, depth):
        if ca >= policy["guard_position"] or cb == ca:
            return
        if cb < ca or depth > policy["max_isolation_depth"]:
            raise ArithmeticError("Nonmonotone or unresolved count isolation")
        if cb-ca == 1:
            fa, fb = determinant(a), determinant(b)
            if fa*fb > 0:
                raise ArithmeticError("Count-one interval lacks primary determinant sign bracket")
            w, info = brentq(determinant, a, b, xtol=policy["root_xtol"],
                             rtol=policy["root_rtol"], full_output=True)
            brackets.append({"bracket": [a, b], "counts": [ca, cb], "signs": [fa, fb],
                "omega": w, "iterations": info.iterations, "converged": info.converged})
            roots.append(w)
            return
        middle, cm = query((a+b)/2)
        if not a < middle < b:
            raise ArithmeticError("Pole-shifted count endpoint left interval")
        isolate(a, middle, ca, cm, depth+1)
        isolate(middle, b, cm, cb, depth+1)
    lower, c0 = query(policy["omega_lower"])
    upper, c1 = query(2*math.pi*policy["frequency_ceiling_fstar"])
    if c0 != 0 or c1 < policy["guard_position"]:
        raise ArithmeticError("Declared interval does not cover first seven roots")
    isolate(lower, upper, c0, c1, 0)
    if len(roots) != policy["guard_position"] or np.any(np.diff(roots) <= policy["root_xtol"]):
        raise ArithmeticError("Incomplete/non-simple prefix; no reordering permitted")
    return roots, {"status": "PASS", "count_samples": samples, "brackets": brackets,
        "lower_count": c0, "upper_count": c1, "computed_roots": len(roots)}


def field_diagnostics(model, lengths, case, stiffness, evaluator, order, elastic_stiffness=None):
    """Physical Hermitian mass/elastic forms, endpoint residuals and fields."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    p, n = model.coefficients, case.nq
    total = sum(lengths)
    mass = energy = cmass = axial_energy = 0.
    qpeak = ppeak = clamp = ode_error = 0.
    ends, values_all, gradients_all, weights_all = [], [], [], []
    for arm, length in enumerate(lengths):
        x, weight = (nodes+1)*length/2, weights*length/2
        value, gradient, operator = evaluator(arm, x)
        endpoint, _, _ = evaluator(arm, np.array([0., length]))
        u, w, theta = (0, 2, 3) if case.is_mh else (0, 1, 2)
        squares = abs(value)**2
        if case.is_mh:
            local_cmass = weight@(p["j"]*squares[:, 1])
            axial_density = p["C"]*(abs(gradient[:, 0])**2+2*model.section.nu*
                np.real(gradient[:, 0].conj()*value[:, 1])+squares[:, 1])+p["H"]*abs(gradient[:, 1])**2
        else:
            seg = reduced.segment(model, length, case.variant)
            local_cmass = weight@(seg.J*abs(gradient[:, 0])**2)
            axial_density = seg.EA*abs(gradient[:, 0])**2
        mass += float(weight@(p["m"]*(squares[:, u]+squares[:, w])+p["r"]*squares[:, theta])+local_cmass)
        cmass += float(local_cmass)
        axial_energy += float(weight@axial_density)
        energy += float(weight@(axial_density+p["B"]*abs(gradient[:, theta])**2+
                                p["S"]*abs(gradient[:, w]-value[:, theta])**2))
        qscales = np.array([1, total, 1, total] if case.is_mh else [1, 1, total])
        pscales = 1/qscales
        qpeak = max(qpeak, float(np.max(abs(value[:, :n])*qscales)))
        ppeak = max(ppeak, float(np.max(abs(value[:, n:])*pscales)))
        clamp = max(clamp, float(np.max(abs(endpoint[0, :n])*qscales)))
        rhs = value@operator.T
        scale = np.maximum(np.max(abs(rhs), axis=0)+np.max(abs(gradient), axis=0), 1e-30)
        ode_error = max(ode_error, float(np.max(abs(gradient-rhs)/scale)))
        ends.append(endpoint[1])
        values_all.append(value)
        gradients_all.append(gradient)
        weights_all.append(weight)
    jump = ends[0][n-1]-ends[1][n-1]
    spring_k = float(np.real(stiffness)) if elastic_stiffness is None else elastic_stiffness
    spring = spring_k*abs(jump)**2 if case.rotation == "SPRING" else 0.
    residual = joint_matrix(case, stiffness)@np.concatenate(ends)
    row_scale = [1/qpeak]*2+[1/ppeak]*2
    row_scale += [total/qpeak, 1/(total*ppeak)] if case.rotation == "RIGID" else [1/(total*ppeak)]*2
    if case.is_mh:
        row_scale += {"FREE": [1/(total*ppeak)]*2,
                      "CONTINUOUS": [total/qpeak, 1/(total*ppeak)],
                      "BLOCKED": [total/qpeak]*2}[case.closure]
    return {"mass": mass, "contraction_or_gradient_mass": cmass, "bulk_energy": energy,
        "axial_energy": axial_energy, "spring_energy": float(spring),
        "elastic_energy": float(energy+spring), "delta_theta": jump,
        "joint_raw": residual, "joint_scaled": abs(residual)*row_scale,
        "clamp_scaled": clamp/qpeak, "ode_scaled": ode_error,
        "end_states": np.array(ends), "values": np.array(values_all),
        "gradients": np.array(gradients_all), "weights": np.array(weights_all)}


def state_operator(model, omega, case):
    if case.is_mh:
        return mh.full_harmonic_state_matrix(model, omega)
    s = reduced.segment(model, .5, case.variant)
    operator = np.zeros((6, 6))
    operator[0, 3] = 1/(s.EA-s.J*omega**2)
    operator[3, 0] = -s.m*omega**2
    operator[np.ix_((1, 2, 4, 5), (1, 2, 4, 5))] = mh.harmonic_state_matrix(model, omega, "timoshenko")
    return operator


def elastic_mode(model, lengths, omega, case, stiffness, order):
    matrix = boundary(model, lengths, omega, case, stiffness)
    _, singular, right = np.linalg.svd(matrix)
    coefficients = right[-1].reshape(2, 2*case.nq)
    operator = state_operator(model, omega, case)
    def evaluator(arm, x):
        return (basis(model, lengths[arm], omega, x, case)@coefficients[arm],
                basis(model, lengths[arm], omega, x, case, 1)@coefficients[arm], operator)
    initial = field_diagnostics(model, lengths, case, stiffness, evaluator, order)
    coefficients /= math.sqrt(initial["mass"])
    diag = field_diagnostics(model, lengths, case, stiffness, evaluator, order)
    diag["energy_relative"] = abs(diag["elastic_energy"]/omega**2-1)
    diag["singular_ratio"] = float(singular[-1]/singular[0])
    diag["nonzero_condition"] = float(singular[0]/singular[-2])
    diag["normalization_multiplier"] = 1/math.sqrt(initial["mass"])
    return coefficients, diag


def mass_vectors(model, case, diag, common=False):
    p, y, weights = model.coefficients, diag["values"], diag["weights"]
    if common or not case.is_mh:
        indices = [0, 2, 3] if case.is_mh else [0, 1, 2]
        vector = (y[:, :, indices]*np.sqrt(weights[:, :, None]*[p["m"], p["m"], p["r"]])).ravel()
        if not common and case.variant == "rayleigh_love_planar":
            vector = np.r_[vector, (diag["gradients"][:, :, 0]*np.sqrt(weights*model.section.nu**2*p["j"])).ravel()]
    else:
        vector = (y[:, :, :4]*np.sqrt(weights[:, :, None]*[p["m"], p["j"], p["m"], p["r"]])).ravel()
    return vector


def validate_mode(diag, policy):
    if (diag["clamp_scaled"] > policy["boundary_tol"] or
        max(diag["joint_scaled"]) > policy["boundary_tol"] or
        diag["ode_scaled"] > policy["ode_tol"] or
        diag["energy_relative"] > policy["energy_tol"] or
        diag["singular_ratio"] > policy["boundary_tol"] or
        diag["nonzero_condition"] > policy["nonzero_condition_max"]):
        raise ArithmeticError("Elastic mode quality gate failed: "+str({k:v for k,v in diag.items() if np.ndim(v) == 0}))


def structural_checks(model, lengths, stiffness):
    """Independent random virtual-work check on admissible variations."""
    from scipy.linalg import null_space
    rng = np.random.default_rng(19261)
    records = []
    for beta in (0., 45.):
        for variant, closures in (("mindlin_herrmann", ("FREE", "CONTINUOUS", "BLOCKED")),
                                  ("elementary", ("CONTINUOUS",)), ("rayleigh_love_planar", ("CONTINUOUS",))):
            for closure in closures:
                case = Case(beta, closure, variant)
                operator = joint_matrix(case, stiffness)
                scaled = operator/np.linalg.norm(operator, axis=1)[:, None]
                pmap, jump = trace_map(case)
                allowed = null_space(scaled)
                error = 0.
                for _ in range(32):
                    state = allowed@rng.normal(size=allowed.shape[1])
                    delta = pmap@rng.normal(size=pmap.shape[1])
                    n = case.nq
                    q = np.r_[state[:n], state[2*n:3*n]]
                    forces = np.r_[state[n:2*n], state[3*n:]]
                    spring_work = stiffness*(q[n-1]-q[2*n-1])*(delta[n-1]-delta[2*n-1])
                    work = forces@delta
                    error = max(error, float(abs(work+spring_work)/max(1., abs(work), abs(spring_work))))
                rank = int(np.linalg.matrix_rank(scaled))
                if rank != 2*case.nq or np.linalg.matrix_rank(pmap) != pmap.shape[1] or error > 2e-13:
                    raise ArithmeticError("Rank/virtual-work structural gate failed")
                raw, diag = schur(model, lengths, .001, case, stiffness)
                # Below first mode: positive index check supplements the analytic coercivity proof.
                if min(diag["eigenvalues"]) <= 0 or diag["skew"] > 1e-10:
                    raise ArithmeticError("Low-frequency positive-energy gate failed")
                # Independent integral identity for the new trace reduction,
                # before any eigenvalue calculation (not just matrix symmetry).
                check_omega = .7
                z = rng.normal(size=pmap.shape[1])
                endpoint_q = (pmap@z).reshape(2,case.nq)
                amplitudes = []
                for i,length in enumerate(lengths):
                    ends = basis(model,length,check_omega,[0.,length],case)
                    system = np.concatenate((ends[0,:case.nq],ends[1,:case.nq]))
                    amplitudes.append(np.linalg.solve(system,np.r_[np.zeros(case.nq),endpoint_q[i]]))
                op = state_operator(model,check_omega,case)
                def extension(arm,x):
                    return (basis(model,lengths[arm],check_omega,x,case)@amplitudes[arm],
                        basis(model,lengths[arm],check_omega,x,case,1)@amplitudes[arm],op)
                energy = field_diagnostics(model,lengths,case,stiffness,extension,200)
                dynamic, _ = schur(model,lengths,check_omega,case,stiffness)
                integral = energy["elastic_energy"]-check_omega**2*energy["mass"]
                contraction = float(z@dynamic@z)
                extension_error = abs(integral-contraction)/max(abs(integral),abs(contraction),1e-30)
                if extension_error>5e-8:
                    raise ArithmeticError("New trace-space harmonic-extension energy identity failed")
                records.append({"case": case.key, "rank": rank, "trace_dimension": pmap.shape[1],
                    "operator": operator.tolist(), "trace_map": pmap.tolist(),
                    "virtual_work_error": error, "low_frequency_schur": diag,
                    "harmonic_extension_check":{"omega":check_omega,"integral_a_minus_lambda_m":integral,
                        "trace_schur_form":contraction,"relative_error":extension_error}})
    elastic, mass = mh.energy_matrices(model)
    return {"status": "PASS", "records": records,
        "local_energy_eigenvalues": np.linalg.eigvalsh(elastic).tolist(),
        "local_mass_eigenvalues": np.linalg.eigvalsh(mass).tolist(),
        "endpoint_signs": [rigid.endpoint_sign("left"), rigid.endpoint_sign("right")]}


class ComplexBasis:
    """Bounded eigen-exponentials of the SAME MH first-order state operator.

    Fixed elastic eigenvalue labels, eigenvector pivots, endpoint anchors and
    positive row scales define a local smooth basis. No long transfer product.
    Only the inertial parameter -omega**2 becomes p**2 (exp(p*t)).
    Existing production real-only functions are never altered or monkeypatched.
    """
    def __init__(self, model, lengths, omega, case, stiffness, tref):
        from scipy.linalg import eig
        if not case.is_mh:
            raise ValueError("Complex confirmation is scoped to MH closures")
        self.model, self.lengths, self.omega, self.case = model, lengths, omega, case
        self.stiffness, self.tref = stiffness, tref
        self.elastic_operator = mh.full_harmonic_state_matrix(model, omega)
        c = model.coefficients
        rate = max(abs(np.linalg.eigvals(self.elastic_operator)))
        length = lengths[0]
        self.scales = np.array([1.,length,1.,length,1/(c["C"]*rate),
            length/(c["H"]*rate),1/(c["S"]*rate),length/(c["B"]*rate)])
        balanced = self.elastic_operator*self.scales[:,None]/self.scales[None,:]
        self.rates, vectors = eig(balanced)
        self.pivots = np.argmax(abs(vectors), axis=0)
        self.anchor_right = self.rates.real > 1e-6
        raw = self.raw_boundary(1j*omega, 0.)
        self.row_scales = 1/np.linalg.norm(raw, axis=1)

    def operator(self, p):
        result = self.elastic_operator.astype(complex)
        coefficients = self.model.coefficients
        change = p*p+self.omega**2
        for row, col, name in ((4,0,"m"),(5,1,"j"),(6,2,"m"),(7,3,"r")):
            result[row,col] += coefficients[name]*change
        return result

    def modes(self, p):
        from scipy.linalg import eig
        from scipy.optimize import linear_sum_assignment
        operator = self.operator(p)
        balanced = operator*self.scales[:,None]/self.scales[None,:]
        rates, vectors = eig(balanced)
        _, assignment = linear_sum_assignment(abs(self.rates[:,None]-rates[None,:]))
        rates, vectors = rates[assignment], vectors[:,assignment]
        pivots = vectors[self.pivots, np.arange(8)]
        if min(abs(pivots)) < 1e-4:
            raise ArithmeticError("Complex eigenbasis left its fixed-pivot neighbourhood")
        vectors = vectors/pivots[None,:]/self.scales[:,None]
        return rates, vectors

    def values(self, p, arm, x, derivative=0):
        rates, vectors = self.modes(p)
        points = np.atleast_1d(x)
        anchors = np.where(self.anchor_right, self.lengths[arm], 0.)
        factors = np.exp((points[:,None]-anchors)*rates)*rates**derivative
        return vectors[None,:,:]*factors[:,None,:]

    def raw_boundary(self, p, ctheta):
        return assemble_boundary([self.values(p, arm, [0.,length]) for arm,length in enumerate(self.lengths)],
                                 self.case, self.stiffness+ctheta*p)

    def matrix(self, z, ctheta, derivative=False):
        p = z/self.tref
        raw = self.raw_boundary(p, ctheta)*self.row_scales[:,None]
        if not derivative:
            return raw
        # Fixed real step in nondimensional z; bases/scales are not re-anchored.
        h = 1e-5
        derivative_matrix = (self.raw_boundary((z+h)/self.tref, ctheta)-
            self.raw_boundary((z-h)/self.tref, ctheta))/(2*h)*self.row_scales[:,None]
        return raw, derivative_matrix


def complex_mode(model, lengths, case, omega, alpha_predictor, stiffness, ctheta,
                 tref, elastic_diag, gap, policy):
    from scripts.lib import inplane_kelvin_voigt as kv
    provider = ComplexBasis(model, lengths, omega, case, stiffness, tref)
    criteria = {"max_steps":policy["complex_newton_max_steps"],
        "newton_residual":policy["complex_newton_residual"], "newton_step":policy["complex_newton_step"]}
    def correct(seed, damping):
        matrix = lambda z, derivative=False: provider.matrix(z, damping, derivative)
        return kv.correct(matrix, seed, kv.right_null(matrix(seed)), criteria)
    elastic = correct(1j*omega*tref, 0.)
    p_elastic = elastic["z"]/tref
    if elastic["status"] != "CONVERGED" or abs(p_elastic-1j*omega)/omega > policy["complex_elastic_relative_tol"]:
        raise ArithmeticError("Independent complex-basis elastic recovery failed")
    corrected = correct((-alpha_predictor+1j*omega)*tref, ctheta)
    if corrected["status"] != "CONVERGED":
        raise ArithmeticError("Local KV corrector failed")
    p = corrected["z"]/tref
    coefficients = corrected["a"].reshape(2,8)
    def evaluator(arm, x):
        return (provider.values(p,arm,x)@coefficients[arm],
                provider.values(p,arm,x,1)@coefficients[arm], provider.operator(p))
    diag = field_diagnostics(model,lengths,case,stiffness+ctheta*p,evaluator,
                             policy["quadrature_order"], elastic_stiffness=stiffness)
    coefficients /= math.sqrt(diag["mass"])
    diag = field_diagnostics(model,lengths,case,stiffness+ctheta*p,evaluator,
                             policy["quadrature_order"], elastic_stiffness=stiffness)
    refined = field_diagnostics(model,lengths,case,stiffness+ctheta*p,evaluator,
                             policy["quadrature_check_order"], elastic_stiffness=stiffness)
    damping = ctheta*abs(diag["delta_theta"])**2
    pencil = p*p*diag["mass"]+p*damping+diag["elastic_energy"]
    energy_error = abs(pencil)/(abs(p*p)*diag["mass"]+abs(p)*damping+diag["elastic_energy"])
    alpha_exact = damping/(2*diag["mass"])
    alpha_error = abs(-p.real/alpha_exact-1)
    a = mass_vectors(model,case,elastic_diag)
    b = mass_vectors(model,case,diag)
    overlap = abs(np.vdot(a,b))**2/(np.vdot(a,a).real*np.vdot(b,b).real)
    predictor_error = abs(-p.real/alpha_predictor-1)
    singular = np.linalg.svd(provider.matrix(corrected["z"],ctheta),compute_uv=False)
    quadrature_error = max(abs(refined[name]/diag[name]-1) for name in ("mass","elastic_energy"))
    if (p.real >= 0 or p.imag <= 0 or abs(p-1j*omega) > policy["complex_neighbour_radius_fraction"]*gap
        or diag["clamp_scaled"] > policy["boundary_tol"] or max(diag["joint_scaled"]) > policy["boundary_tol"]
        or diag["ode_scaled"] > policy["ode_tol"] or energy_error > policy["energy_tol"]
        or alpha_error > policy["energy_tol"] or quadrature_error > policy["energy_tol"]
        or predictor_error > policy["complex_predictor_relative_tol"] or overlap < policy["shape_mac_min"]
        or singular[-1]/singular[0] > policy["boundary_tol"]
        or singular[0]/singular[-2] > policy["nonzero_condition_max"]):
        raise ArithmeticError(f"Complex physical gate failed: p={p}, energy={energy_error}, alpha={alpha_error}, "
            f"joint={max(diag['joint_scaled'])}, clamp={diag['clamp_scaled']}, ode={diag['ode_scaled']}, predictor={predictor_error}")
    return {"status":"PASS", "p":p,"z":p*tref,"coefficients":coefficients,
        "alpha":-p.real,"omega_d":p.imag,"alpha_predictor":alpha_predictor,
        "predictor_relative_error":predictor_error,"elastic_basis_p":p_elastic,
        "elastic_frequency_relative_error":abs(p_elastic-1j*omega)/omega,
        "mass_MAC_with_elastic":float(overlap), "energy_identity_relative":float(energy_error),
        "alpha_identity_relative":float(alpha_error),"quadrature_relative":float(quadrature_error),
        "singular_ratio":float(singular[-1]/singular[0]),
        "nonzero_condition":float(singular[0]/singular[-2]),
        "newton_history":corrected["history"],"newton_steps":corrected["steps"],
        "diagnostics":diag}, provider
