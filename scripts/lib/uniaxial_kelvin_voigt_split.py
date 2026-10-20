"""Linear-only spatial extension of the accepted planar FREE/KV joint.

Canonical spatial q=(u,w,v,Phi,psi,theta,c), y=(q,f). No historical
rigid/c-continuous or nonlinear model is modified. The out-of-plane arm
basis reuses rotated-section Tim and scalar elementary torsion functions.
"""
from __future__ import annotations

from fractions import Fraction
from types import SimpleNamespace
import math
import numpy as np
from scipy.linalg import block_diag, expm, null_space
from scipy.optimize import brentq

from scripts.analysis import verify_weakly_nonlinear_spatial_rod as spatial
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as frames
from scripts.lib import contraction_joint_sensitivity as planar
from scripts.lib import yartsev_ch2_coupled_rods as book_joint
from scripts.lib import bishop_longitudinal as scalar
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section

PQ = [0, 6, 1, 5]  # Existing MH order u,c,w,theta, explicit diagnostic permutation.
OQ = [2, 4, 3]     # Existing out-of-plane order v,psi,Phi.
PY, OY = PQ+[j+7 for j in PQ], OQ+[j+7 for j in OQ]
PC = PY+[j+14 for j in PY]
OC = OY+[j+14 for j in OY]
PR = [0, 1, 3, 4, 10, 11, 12, 13]
OR = [2, 6, 7, 5, 8, 9]
OSIGN = np.array([-1., 1., 1., -1., 1., 1.])
NORMAL = np.array([0., 0., -1.])


def quadratic_action(saved):
    """Exact degree-two extraction of the hash-verified serialized action."""
    t = rod.Polynomial.deserialize(saved["polynomials"]["T4"]).homogeneous(2)
    v = rod.Polynomial.deserialize(saved["polynomials"]["V4"]).homogeneous(2)
    s = {n: rod.Polynomial.symbol(n) for n in rod.SYMBOL_ORDER}
    half = Fraction(1, 2)
    expected_t = half*(s["m"]*(s["u_t"]**2+s["w_t"]**2+s["v_t"]**2)
        +s["jp"]*s["c_t"]**2+(s["jp"]+s["jb"])*s["Phi_t"]**2
        +s["jb"]*s["psi_t"]**2+s["jp"]*s["theta_t"]**2)
    expected_v = half*(s["C"]*(s["u_s"]**2+2*s["nu"]*s["u_s"]*s["c"]+s["c"]**2)
        +s["H"]*s["c_s"]**2+s["S"]*(s["w_s"]-s["theta"])**2+s["Bp"]*s["theta_s"]**2
        +s["S"]*(s["v_s"]-s["psi"])**2+s["Bb"]*s["psi_s"]**2+s["CT"]*s["Phi_s"]**2)
    if t != expected_t or v != expected_v:
        raise ArithmeticError("Saved quadratic action differs from the declared linear limit")
    coefficient = rod.RodCoefficients(**saved["coefficients"])
    harmonic, matrices = spatial.linear_state_from_action(coefficient, SimpleNamespace(T4=t, V4=v))
    if matrices != saved["linear_action_matrices"]:
        raise ArithmeticError("Quadratic Hessians differ from historical linear-action matrices")
    mixed = {name: float(max(np.max(abs(np.asarray(a)[np.ix_(PQ, OQ)])),
                             np.max(abs(np.asarray(a)[np.ix_(OQ, PQ)])))) for name, a in matrices.items()}
    if any(mixed.values()):
        raise ArithmeticError("Mixed quadratic coefficients")
    return coefficient, lambda p: harmonic(-1j*p), {"status": "PASS", "exact_T2": t.serialize(),
        "exact_V2": v.serialize(), "matrices": matrices, "mixed_blocks": mixed}


def joint_operator(beta, stiffness, transforms=None):
    """14x28 physical-vector conditions, with finite normal rotational impedance.

    Rows: dXYZ, FXYZ, rotationXY, momentXY, normal moment1/2, R1/R2.
    Transform rows contain physical moments dual to the named rotations.
    """
    transforms = spatial.frame_transforms(beta) if transforms is None else transforms
    result = np.zeros((14, 28), dtype=np.result_type(stiffness, float))
    for arm, transform in enumerate(transforms):
        q, f = slice(14*arm, 14*arm+7), slice(14*arm+7, 14*arm+14)
        sign = 1 if arm == 0 else -1
        result[:3, q] = sign*transform[:3]
        result[3:6, f] = transform[:3]
        result[6:8, q] = sign*transform[3:5]
        result[8:10, f] = transform[3:5]
        result[10+arm, f] = NORMAL@transform[3:6]
        result[12+arm, 14*arm+13] = 1.
    jump = np.zeros(28)
    jump[5], jump[19] = 1., -1.
    result[10] += stiffness*jump
    result[11] -= stiffness*jump
    return result


def trace_map(beta, transforms=None):
    """14 local endpoint coordinates from (dXYZ,alphaXY,theta1,theta2,c1,c2)."""
    transforms = spatial.frame_transforms(beta) if transforms is None else transforms
    result = np.zeros((14, 9))
    for i, transform in enumerate(transforms):
        global_map = np.zeros((7, 9))
        global_map[:3, :3] = np.eye(3)
        global_map[3:5, 3:5] = np.eye(2)
        global_map[5, 5+i] = -1.  # global alphaZ=-theta
        global_map[6, 7+i] = 1.
        result[7*i:7*i+7] = transform.T@global_map
    return result


def out_joint(beta, stiffness=0., rigid=False):
    if rigid:
        raw = spatial.joint_operator(beta)[np.ix_([2, 3, 4, 9, 10, 11], OC)]
    else:
        raw = joint_operator(beta, stiffness)[np.ix_(OR, OC)]
    return raw*OSIGN[:, None]


def assemble(endpoint_bases, joint):
    n = endpoint_bases[0].shape[1]//2
    result = np.zeros((4*n, 4*n), dtype=np.result_type(*endpoint_bases, joint))
    for i, values in enumerate(endpoint_bases):
        result[i*n:(i+1)*n, i*2*n:(i+1)*2*n] = values[0, :n]
    result[2*n:] = joint@block_diag(*[v[1] for v in endpoint_bases])
    return result


def _rank_fraction(matrix):
    """Exact rational row reduction for the 0/90 degree endpoint checks."""
    a = [[Fraction(v) for v in row] for row in matrix]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(rank, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        value = a[rank][col]
        a[rank] = [v/value for v in a[rank]]
        for i in range(rank+1, len(a)):
            value = a[i][col]
            a[i] = [v-value*w for v, w in zip(a[i], a[rank])]
        rank += 1
        if rank == len(a):
            break
    return rank


def structural_checks(operator, coefficients, mh_model, kv_result, policy):
    rng = np.random.default_rng(20261019)
    k, c = kv_result["k_theta"], kv_result["c_theta"]
    source_mh = [r for r in kv_result["complex_modes"] if r["variant"] == "mindlin_herrmann"]
    checks = []
    for beta in (0., 45., 90.):
        transforms = spatial.frame_transforms(beta)
        frame_errors = [max(np.max(abs(t.T@t-np.eye(7))),
                            abs(np.linalg.det(t[:3, :3])-1),
                            np.max(abs(np.cross(t[:3, 0], t[:3, 1])-NORMAL))) for t in transforms]
        trace = trace_map(beta)
        power_error = duality_error = dissipation_error = 0.
        J = joint_operator(beta, k)
        null = null_space(J/np.linalg.norm(J, axis=1)[:, None])
        for _ in range(32):
            y = null@rng.normal(size=null.shape[1])
            z = rng.normal(size=9)
            variation = (trace@z).reshape(2, 7)
            local = sum(y[14*i+7:14*i+14]@variation[i] for i in range(2))
            jump, delta_jump = y[5]-y[19], variation[0, 5]-variation[1, 5]
            power_error = max(power_error, abs(local+k*jump*delta_jump))
            viscous_efforts = y.copy()
            viscous_efforts[12] -= c*delta_jump
            viscous_efforts[26] += c*delta_jump
            actual_power = sum(viscous_efforts[14*i+7:14*i+14]@variation[i] for i in range(2))
            dissipation_error = max(dissipation_error, abs(actual_power+k*jump*delta_jump+c*delta_jump**2))
            for t in transforms:
                f, dq = rng.normal(size=(2, 7))
                duality_error = max(duality_error, abs(f@dq-(t@f)@(t@dq)))
        # Swap the arm data and their physical frames, retaining global axes.
        swapped = joint_operator(beta, k, transforms[::-1])[:, list(range(14, 28))+list(range(14))]
        rows = list(range(10))+[11, 10, 13, 12]
        signs = np.ones(14); signs[:3] = -1.; signs[6:8] = -1.
        swap_error = float(np.max(abs(swapped-J[rows]*signs[:, None])))
        exact_rank = None
        exact_out = None
        if beta in (0., 90.):
            exact_transforms = [np.rint(t) for t in transforms]
            exact_rank = _rank_fraction(joint_operator(beta, 2., exact_transforms).astype(int))
            exact_out = float(np.max(abs(out_joint(beta)-np.rint(book_joint.joint_matrix_book(math.radians(beta))))))
        record = {"beta": beta, "frames": transforms, "frame_errors": frame_errors,
            "joint_rank": int(np.linalg.matrix_rank(J/np.linalg.norm(J, axis=1)[:, None])),
            "essential_trace_dimension": int(np.linalg.matrix_rank(trace)),
            "duality_absolute": duality_error, "joint_virtual_work_absolute": power_error,
            "power_balance_absolute": dissipation_error,
            "arm_swap_max_error": swap_error, "exact_integer_rank": exact_rank,
            "endpoint_exact_out_error": exact_out, "spectral_matrices": []}
        for source_mode in source_mh:
            p = complex(source_mode["p"]["real"], source_mode["p"]["imag"])
            omega = source_mode["omega_elastic"]
            reference = planar.ComplexBasis(mh_model, [.5, .5], omega, planar.Case(beta, "FREE"), k, kv_result["t_ref"])
            a = operator(p)
            units = np.array([1., 1., 1., coefficients.S, coefficients.Bb, coefficients.CT])
            ao = a[np.ix_(OY, OY)]
            balanced = ao*units[None, :]/units[:, None]
            out_end = np.array([units[:, None]*expm(balanced*x) for x in (0., .5)])
            ends = []
            for arm in range(2):
                end = np.zeros((2, 14, 14), dtype=complex)
                for ix in range(2):
                    end[ix][np.ix_(PY, PY)] = reference.values(p, arm, [0., .5])[ix]
                    end[ix][np.ix_(OY, OY)] = out_end[ix]
                ends.append(end)
            full = assemble(ends, joint_operator(beta, k+p*c))
            prows = PQ+[j+7 for j in PQ]+[j+14 for j in PR]
            orows = OQ+[j+7 for j in OQ]+[j+14 for j in OR]
            bp, bo = full[np.ix_(prows, PC)], full[np.ix_(orows, OC)].copy()
            bo[6:] *= OSIGN[:, None]
            rigid_bo = assemble([out_end, out_end], out_joint(beta, rigid=True))
            stiffness_derivative = joint_operator(beta, 1.)-joint_operator(beta, 0.)
            record["spectral_matrices"].append({"p": p,
                "local_mixed_max": float(max(np.max(abs(a[np.ix_(PY, OY)])), np.max(abs(a[np.ix_(OY, PY)])))),
                "boundary_mixed_max": float(max(np.max(abs(full[np.ix_(prows, OC)])), np.max(abs(full[np.ix_(orows, PC)])))),
                "parallel_old_matrix_max": float(np.max(abs(bp-reference.raw_boundary(p, c)))),
                "out_rigid_matrix_max": float(np.max(abs(bo-rigid_bo))),
                "out_parameter_derivative_max": float(np.max(abs(stiffness_derivative[np.ix_(OR, OC)]))),
                "joint_rank": int(np.linalg.matrix_rank(joint_operator(beta, k+p*c))),
                "out_joint_book_max": float(np.max(abs(out_joint(beta)-book_joint.joint_matrix_book(math.radians(beta))))),
                "parallel_operator_max": float(np.max(abs(a[np.ix_(PY, PY)]-reference.operator(p))))})
        # Mutants intentionally violate the stipulated one-axis law. No roots.
        direction = np.zeros(14); direction[5], direction[12], direction[3] = 1., -1., .25
        mutant = np.outer(direction, direction)
        qp, qo = PQ+[j+7 for j in PQ], OQ+[j+7 for j in OQ]
        mixed_mutant = float(np.max(abs(mutant[np.ix_(qp, qo)])))
        bad_joint = out_joint(beta).copy(); bad_joint[4, 2] += .25
        out_mutant = float(np.max(abs(bad_joint-out_joint(beta))))
        record["negative_controls"] = {"rotation_coupling_mixed": mixed_mutant, "out_moment_change": out_mutant}
        checks.append(record)
        metrics = [*frame_errors, power_error, duality_error, dissipation_error, swap_error]
        for row in record["spectral_matrices"]:
            metrics += [v for name, v in row.items() if name.endswith("_max")]
        if (max(metrics) > policy["structure_tol"] or record["joint_rank"] != 14
            or record["essential_trace_dimension"] != 9 or exact_rank not in (None, 14)
            or mixed_mutant == 0 or out_mutant == 0):
            raise ArithmeticError(f"Structural gate failed at beta={beta}: {record}")
    return {"status": "PASS", "geometries": checks, "canonical_state": "u,w,v,Phi,psi,theta,c,N,Qw,Qv,MT,Mpsi,Mtheta,R",
        "parallel_indices": PY, "out_indices": OY, "parallel_joint_rows": PR, "out_joint_rows": OR,
        "out_row_signs": OSIGN, "outer_clamp": "all seven coordinates zero; no slope constraints"}


class OutOfPlane:
    """One fixed G20 section; no angular sweep or repeated KV spectrum API."""
    def __init__(self, coefficients, operator):
        self.c, self.operator = coefficients, operator
        self.tim = mh.project_jang_reduced_rectangular(rectangular_section(
            E=1., rho=1., nu=.3, width=.05, thickness=.20, K=5/6))
        self.torsion = scalar.Segment(.5, coefficients.CT, coefficients.jp+coefficients.jb, H=0., J=0.)

    def basis(self, omega, x, derivative=0):
        x = np.atleast_1d(x)
        result = np.zeros((len(x), 6, 6))
        result[:, [0, 1, 3, 4], :4] = mh.finite_state_basis(self.tim, .5, omega, x, "timoshenko", derivative)
        result[:, 2, 4:] = scalar.basis(self.torsion, omega, x, derivative)
        result[:, 5, 4:] = self.c.CT*scalar.basis(self.torsion, omega, x, derivative+1)
        return result

    def boundary(self, omega, rigid=True, stiffness=0.):
        end = self.basis(omega, [0., .5])
        raw = assemble([end, end], out_joint(45., stiffness, rigid))
        return raw/np.linalg.norm(raw, axis=1)[:, None]

    def arm_stiffness(self, omega):
        end = self.basis(omega, [0., .5])
        boundary = np.concatenate((end[0, :3], end[1, :3]))
        norms = np.linalg.norm(boundary, axis=1)
        coeff = np.linalg.solve(boundary/norms[:, None], np.vstack((np.zeros((3, 3)), np.eye(3)))/norms[:, None])
        return end[1, 3:]@coeff, float(np.linalg.cond(boundary/norms[:, None])), coeff

    def schur(self, omega):
        arm, condition, coefficients = self.arm_stiffness(omega)
        maps = []
        for frame in frames.frames(45.):
            maps.append(np.array([[1., 0., 0.], [0., -frame.n[0], frame.t[0]], [0., -frame.n[1], frame.t[1]]]).T)
        trace = np.vstack(maps)
        raw = trace.T@block_diag(arm, arm)@trace
        scales = np.sqrt([1/self.c.S, 1/self.c.CT, 1/self.c.Bb])
        balanced = raw*scales[:, None]*scales[None, :]
        skew = np.linalg.norm(balanced-balanced.T, 2)/np.linalg.norm(balanced, 2)
        eig = np.linalg.eigvalsh((balanced+balanced.T)/2)
        return raw, {"skew": float(skew), "condition": condition, "margin": float(min(abs(eig))/max(abs(eig))),
                     "eigenvalues": eig, "negative_inertia": int(sum(eig < 0))}, coefficients, maps

    def pole_catalog(self, policy):
        roots, search = mh.finite_roots(self.tim, .5, "timoshenko", policy["omega_lower"], policy["omega_upper"],
            {**policy, "scan_intervals": policy["arm_scan_intervals"]})
        upper = mh.finite_count_upper_bound(self.tim, .5, policy["omega_upper"], "timoshenko")
        lower = mh.finite_count_upper_bound(self.tim, .5, policy["omega_lower"], "timoshenko")
        speed = math.sqrt(self.c.CT/(self.c.jp+self.c.jb))
        torsion = [n*math.pi*speed/.5 for n in range(1, int(policy["omega_upper"]*.5/(math.pi*speed))+1)]
        status = "PASS" if len(roots) == upper["upper_count"] and lower["upper_count"] == 0 else "SPECTRAL_UNRESOLVED"
        return {"status": status, "tim_roots": roots, "tim_search": search, "tim_count_bound": upper,
                "tim_lower_count": lower, "torsion_exact_poles": torsion,
                "seventh_trial_upper_bound": 4*math.pi*speed/.5,
                "all_poles_with_arm_multiplicity": sorted(2*([r["omega"] for r in roots]+torsion))}

    def prefix(self, catalog, policy, checkpoint):
        """N(w)=J0(w)+negative inertia(S(w)); only the first7 are isolated.

        A narrow adapter of the existing count-guided Brent path: the trace
        now contains (vJ,alphaX,alphaY), not the MH contraction coordinates.
        """
        samples, roots = checkpoint.setdefault("count_samples", []), checkpoint.setdefault("roots", [])
        poles = catalog["all_poles_with_arm_multiplicity"]
        def query(w):
            if len(samples) >= policy["max_count_queries"]:
                raise ArithmeticError("Count query budget exhausted")
            effective = w
            exclusion = policy["pole_exclusion_relative"]*max(1., w)
            for pole in poles:
                if abs(effective-pole) <= exclusion:
                    effective = pole+2*exclusion
            _, d, _, _ = self.schur(effective)
            if d["skew"] > policy["schur_symmetry_tol"] or d["margin"] < policy["count_inertia_margin_min"] or d["condition"] > policy["pole_condition_max"]:
                raise ArithmeticError(f"Schur count unresolved: {d}")
            j0 = sum(pole < effective for pole in poles)
            count = j0+d["negative_inertia"]
            samples.append({"requested_omega": w, "effective_omega": effective, "J0": j0, "count": count, **d})
            return effective, count
        def isolate(a, b, ca, cb, depth=0):
            if ca >= 7 or ca == cb:
                return
            if cb < ca or depth > policy["max_isolation_depth"]:
                raise ArithmeticError("Count isolation unresolved")
            if cb-ca == 1:
                f = lambda w: float(np.linalg.det(self.boundary(w)))
                fa, fb = f(a), f(b)
                if fa*fb > 0:
                    raise ArithmeticError("Count-one interval lacks determinant sign bracket")
                root, info = brentq(f, a, b, xtol=policy["root_xtol"], rtol=policy["root_rtol"], full_output=True)
                roots.append({"omega": root, "bracket": [a, b], "counts": [ca, cb], "iterations": info.iterations,
                              "converged": info.converged})
                return
            middle, cm = query((a+b)/2)
            if not a < middle < b:
                raise ArithmeticError("Pole avoidance left isolation interval")
            isolate(a, middle, ca, cm, depth+1); isolate(middle, b, cm, cb, depth+1)
        lower, c0 = query(policy["omega_lower"]); upper, c1 = query(policy["omega_upper"])
        if c0 != 0 or c1 < 7:
            raise ArithmeticError("Declared range does not certify the first7")
        isolate(lower, upper, c0, c1)
        if len(roots) != 7 or np.any(np.diff([r["omega"] for r in roots]) <= policy["root_xtol"]):
            raise ArithmeticError("Incomplete/non-simple prefix")
        checkpoint.update(status="PASS", lower_count=c0, upper_count=c1)
        return roots

    def fields(self, omega, coefficients, order):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        x, weights = (nodes+1)*.25, weights*.25
        values = np.array([self.basis(omega, x)@a for a in coefficients])
        gradients = np.array([self.basis(omega, x, 1)@a for a in coefficients])
        end = np.array([self.basis(omega, [0., .5])@a for a in coefficients])
        c = self.c
        mass = np.sum(weights*(c.m*values[:, :, 0]**2+c.jb*values[:, :, 1]**2+(c.jp+c.jb)*values[:, :, 2]**2))
        energy = np.sum(weights*(c.S*(gradients[:, :, 0]-values[:, :, 1])**2+c.Bb*gradients[:, :, 1]**2+c.CT*gradients[:, :, 2]**2))
        a = self.operator(1j*omega)[np.ix_(OY, OY)].real
        rhs = values@a.T
        ode = np.max(abs(gradients-rhs)/np.maximum(np.max(abs(gradients)+abs(rhs), axis=(0, 1)), 1e-30))
        qpeak, fpeak = np.max(abs(values[:, :, :3])), np.max(abs(values[:, :, 3:]))
        joint = out_joint(45., rigid=True)@end[:, 1].ravel()
        return {"x": x, "weights": weights, "values": values, "gradients": gradients,
            "mass": float(mass), "energy": float(energy), "energy_relative": float(abs(energy/(omega**2*mass)-1)),
            "ode_scaled": float(ode), "clamp_scaled": float(np.max(abs(end[:, 0, :3]))/qpeak),
            "joint_raw": joint, "joint_scaled": abs(joint)/np.r_[[qpeak]*3, [fpeak]*3], "end_states": end}

    def mode(self, omega, stiffness, policy):
        boundary = self.boundary(omega)
        _, singular, vh = np.linalg.svd(boundary)
        coefficients = vh[-1].reshape(2, 6)
        first = self.fields(omega, coefficients, policy["quadrature_order"])
        coefficients /= math.sqrt(first["mass"])
        d = self.fields(omega, coefficients, policy["quadrature_order"])
        fine = self.fields(omega, coefficients, policy["quadrature_check_order"])
        qerror = max(abs(fine[n]/d[n]-1) for n in ("mass", "energy"))
        kv_joint = out_joint(45., stiffness)@d["end_states"][:, 1].ravel()
        # Compare in exactly the SAME coordinates, arithmetic and row scales.
        # Separate real/complex norm reductions can differ by one ulp even
        # for identical physical rows. This does not relax the exact-zero gate.
        end_basis = self.basis(omega, [0., .5]).astype(complex)
        rigid_raw = assemble([end_basis, end_basis], out_joint(45., rigid=True).astype(complex))
        kv_raw = assemble([end_basis, end_basis], out_joint(45., stiffness).astype(complex))
        equality = float(np.max(abs((kv_raw-rigid_raw)/np.linalg.norm(rigid_raw, axis=1)[:, None])))
        d.update(quadrature_relative=qerror, kv_joint_raw=kv_joint,
            rigid_kv_boundary_max=equality, alpha_exact=0., delta_theta_exact=0.,
            complex_energy_relative=abs(-omega**2*d["mass"]+d["energy"])/(omega**2*d["mass"]+d["energy"]),
            singular_ratio=float(singular[-1]/singular[0]), nonzero_condition=float(singular[0]/singular[-2]))
        if (max(d["joint_scaled"]) > policy["boundary_tol"] or d["clamp_scaled"] > policy["boundary_tol"]
            or d["ode_scaled"] > policy["ode_tol"] or d["energy_relative"] > policy["energy_tol"]
            or qerror > policy["energy_tol"] or equality != 0 or d["singular_ratio"] > policy["boundary_tol"]
            or d["nonzero_condition"] > policy["nonzero_condition_max"]):
            raise ArithmeticError("Physical spectral gate failed")
        return coefficients, d

    def green_check(self, omega, order=200):
        raw, _, coefficients, maps = self.schur(omega)
        z = np.array([.7, -.3, .2])
        amplitudes = np.array([coefficients@mapping@z for mapping in maps])
        d = self.fields(omega, amplitudes, order)
        lhs, rhs = z@raw@z, d["energy"]-omega**2*d["mass"]
        return float(abs(lhs-rhs)/max(abs(lhs), abs(rhs), 1e-30))
