"""Local elementary/planar RL complex continuation; unchanged Tim/Jang signs.

Only the four requested G20 targets are exposed by the diagnostic CLI.
Assembly, mass/energy quadrature and augmented Newton reuse existing helpers.
The bounded basis extends the previous MH approach to six states, with an
analytic derivative obtained from differentiated, fixed-pivot eigenpairs.
"""
from __future__ import annotations

import math
import numpy as np
from scipy.linalg import eig, block_diag
from scipy.optimize import linear_sum_assignment

from scripts.lib import contraction_joint_sensitivity as old
from scripts.lib import coupled_longitudinal_comparators as reduced
from scripts.lib import inplane_kelvin_voigt as kv


def complex_operator(model, case, p, derivative=False):
    """y=(u,w,theta,N,Q,M); N_RL=(EA+J*p**2)*u', N'=m*p**2*u."""
    if case.is_mh:
        raise ValueError("MH roots are read from the accepted bundle, never solved here")
    seg, c = reduced.segment(model, .5, case.variant), model.coefficients
    result = np.zeros((6, 6), dtype=complex)
    axial = seg.EA+seg.J*p*p
    if derivative:
        result[0, 3] = -2*seg.J*p/axial**2
        result[3, 0] = 2*seg.m*p
        result[4, 1], result[5, 2] = 2*c["m"]*p, 2*c["r"]*p
    else:
        # Copy the verified undamped operator, replacing its inertial entries.
        # The real API requires omega>0; all four frequency-dependent entries
        # are overwritten below, so this reference value introduces no scale.
        result[:] = old.state_operator(model, 1., case)
        result[0, 3], result[3, 0] = 1/axial, seg.m*p*p
        result[4, 1], result[5, 2] = c["m"]*p*p, c["r"]*p*p
    return result


class ReducedComplexBasis:
    def __init__(self, model, lengths, omega, case, stiffness, tref):
        self.model, self.lengths, self.omega, self.case = model, lengths, omega, case
        self.stiffness, self.tref = stiffness, tref
        a = self.operator(1j*omega)
        c, length = model.coefficients, lengths[0]
        rate = max(abs(np.linalg.eigvals(a)))
        self.scales = np.array([1., 1., length, 1/(model.section.EA*rate),
                               1/(c["S"]*rate), length/(c["B"]*rate)])
        self.rates, vectors = eig(self.balance(a))
        self.pivots = np.argmax(abs(vectors), axis=0)
        self.anchor_right = self.rates.real > 1e-6
        self.row_scales = np.ones(12)
        self.row_scales = 1/np.linalg.norm(self.raw_boundary(1j*omega, 0.), axis=1)

    def operator(self, p, derivative=False):
        return complex_operator(self.model, self.case, p, derivative)

    def balance(self, a):
        return a*self.scales[:, None]/self.scales[None, :]

    def modes(self, p, derivative=False):
        a = self.balance(self.operator(p))
        rates, vectors = eig(a)
        _, order = linear_sum_assignment(abs(self.rates[:, None]-rates[None, :]))
        rates, vectors = rates[order], vectors[:, order]
        pivot = vectors[self.pivots, np.arange(6)]
        if min(abs(pivot)) < 1e-4:
            raise ArithmeticError("Fixed eigenvector pivot neighbourhood exceeded")
        vectors /= pivot[None, :]
        if not derivative:
            return rates, vectors/self.scales[:, None]
        ap = self.balance(self.operator(p, True))
        vp, rp = np.zeros_like(vectors), np.zeros(6, dtype=complex)
        for j in range(6):
            v = vectors[:, j]
            gauge = np.eye(6)[self.pivots[j]]
            bordered = np.block([[a-rates[j]*np.eye(6), -v[:, None]],
                                 [gauge[None, :], np.zeros((1, 1))]])
            solution = np.linalg.solve(bordered, np.r_[-ap@v, 0.])
            vp[:, j], rp[j] = solution[:6], solution[6]
        return rates, vectors/self.scales[:, None], rp, vp/self.scales[:, None]

    def values(self, p, arm, x, derivative=0, spectral_derivative=False):
        modes = self.modes(p, spectral_derivative)
        rates, vectors = modes[:2]
        distance = np.atleast_1d(x)[:, None]-np.where(self.anchor_right, self.lengths[arm], 0.)
        exponential = np.exp(distance*rates)
        value = vectors[None, :, :]*(exponential*rates**derivative)[:, None, :]
        if not spectral_derivative:
            return value
        if derivative:
            raise ValueError("Mixed derivatives are not needed by this local provider")
        rp, vp = modes[2:]
        value_p = (vp[None, :, :]+vectors[None, :, :]*(distance*rp)[:, None, :])*exponential[:, None, :]
        return value, value_p

    def raw_boundary(self, p, ctheta, derivative=False):
        ends = [self.values(p, i, [0., length], spectral_derivative=derivative)
                for i, length in enumerate(self.lengths)]
        stiffness = self.stiffness+p*ctheta
        if not derivative:
            return old.assemble_boundary(ends, self.case, stiffness)
        values, slopes = zip(*ends)
        value = old.assemble_boundary(values, self.case, stiffness)
        slope = old.assemble_boundary(slopes, self.case, stiffness)
        joint_p = ctheta*(old.joint_matrix(self.case, 1.)-old.joint_matrix(self.case, 0.))
        slope[6:] += joint_p@block_diag(*[v[1] for v in values])
        return value, slope

    def matrix(self, z, ctheta, derivative=False):
        result = self.raw_boundary(z/self.tref, ctheta, derivative)
        if not derivative:
            return result*self.row_scales[:, None]
        value, slope = result
        return value*self.row_scales[:, None], slope*self.row_scales[:, None]/self.tref


def saved_elastic_diagnostics(model, lengths, case, root, stiffness, order):
    """Evaluate saved coefficients; no SVD, root search or shape reconstruction."""
    coefficients, omega = np.array(root["coefficients"]), root["omega"]
    def evaluator(arm, x):
        return (old.basis(model, lengths[arm], omega, x, case)@coefficients[arm],
                old.basis(model, lengths[arm], omega, x, case, 1)@coefficients[arm],
                old.state_operator(model, omega, case))
    return old.field_diagnostics(model, lengths, case, stiffness, evaluator, order)


def overlap(a, b):
    return float(abs(np.vdot(a, b))**2/(np.vdot(a, a).real*np.vdot(b, b).real))


def projection_correspondence(model, first_case, second_case, first, second, policy):
    """Common u,w,theta projection only, NOT the full mass of either theory."""
    va = [old.mass_vectors(model, first_case, d, common=True) for d in first]
    vb = [old.mass_vectors(model, second_case, d, common=True) for d in second]
    matrix = np.array([[overlap(a, b) for b in vb] for a in va])
    margins = [matrix[i, i]-max(np.max(np.delete(matrix[i], i)),
                               np.max(np.delete(matrix[:, i], i))) for i in range(6)]
    accepted = bool(min(np.diag(matrix)) >= policy["shape_mac_min"] and
                    min(margins) >= policy["shape_margin_min"])
    return {"metric": "common physical u,w,theta projection weighted by m,m,r; excludes c and u' inertia",
            "matrix": matrix, "diagonal": np.diag(matrix), "margins": margins,
            "status": "PASS" if accepted else "NUMERICAL_UNRESOLVED",
            "semantics": "fixed-case correspondence; no branch continuation or root permutation"}


def provider_preflight(provider, ctheta, policy):
    z = 1j*provider.omega*provider.tref
    h = policy["derivative_check_step_z"]
    errors = []
    # Check the real and imaginary differentiation directions near both seeds.
    for damping in (0., ctheta):
        zz = z-(.001 if damping else 0.)
        _, analytic = provider.matrix(zz, damping, True)
        for direction in (1., 1j):
            numeric = (provider.matrix(zz+h*direction, damping)-
                       provider.matrix(zz-h*direction, damping))/(2*h*direction)
            errors.append(float(np.linalg.norm(numeric-analytic)/np.linalg.norm(analytic)))
    real_error = float(np.max(abs(provider.operator(1j*provider.omega)-
                                 old.state_operator(provider.model, provider.omega, provider.case))))
    return {"derivative_relative_errors": errors, "elastic_operator_absolute_error": real_error,
            "status": "PASS" if max(errors) <= policy["derivative_relative_tol"] and real_error < 1e-12
            else "NUMERICAL_UNRESOLVED"}


def complex_target(provider, alpha_pred, elastic_diags, index, gap, ctheta, policy, ledger):
    """One c=0 provider recovery and one damped corrector; no retries."""
    criteria = {"max_steps": policy["complex_newton_max_steps"],
                "newton_residual": policy["complex_newton_residual"],
                "newton_step": policy["complex_newton_step"]}
    def correct(seed, damping, kind):
        ledger[kind] += 1
        limit = policy["max_new_damped_correctors" if damping else "max_elastic_provider_correctors"]
        if ledger[kind] > limit:
            raise RuntimeError("Declared root-corrector budget exceeded")
        matrix = lambda z, derivative=False: provider.matrix(z, damping, derivative)
        return kv.correct(matrix, seed, kv.right_null(matrix(seed)), criteria)
    omega, tref = provider.omega, provider.tref
    recovered = correct(1j*omega*tref, 0., "elastic_provider_correctors")
    recovery_error = abs(recovered["z"]/tref-1j*omega)/omega
    record = {"elastic_recovery": recovered, "elastic_frequency_relative_error": recovery_error,
              "status": "NUMERICAL_UNRESOLVED"}
    if recovered["status"] != "CONVERGED" or recovery_error > policy["complex_elastic_relative_tol"]:
        return record, None
    solved = correct((-alpha_pred+1j*omega)*tref, ctheta, "new_damped_correctors")
    p, coefficients = solved["z"]/tref, solved["a"].reshape(2, 6)
    record.update(p=p, newton=solved)
    def evaluator(arm, x):
        return (provider.values(p, arm, x)@coefficients[arm],
                provider.values(p, arm, x, 1)@coefficients[arm], provider.operator(p))
    def diagnostics(order):
        return old.field_diagnostics(provider.model, provider.lengths, provider.case,
            provider.stiffness+ctheta*p, evaluator, order, elastic_stiffness=provider.stiffness)
    initial_mass = diagnostics(policy["quadrature_order"])["mass"]
    coefficients /= math.sqrt(initial_mass)
    d, check = diagnostics(policy["quadrature_order"]), diagnostics(policy["quadrature_check_order"])
    damping = ctheta*abs(d["delta_theta"])**2
    energy_residual = p*p*d["mass"]+p*damping+d["elastic_energy"]
    energy = abs(energy_residual)/(abs(p*p)*d["mass"]+abs(p)*damping+d["elastic_energy"])
    alpha_identity = abs(-p.real/(damping/(2*d["mass"]))-1)
    prediction = abs(-p.real/alpha_pred-1)
    vector = old.mass_vectors(provider.model, provider.case, d)
    macs = [overlap(vector, old.mass_vectors(provider.model, provider.case, e)) for e in elastic_diags]
    margin = macs[index]-max(np.delete(macs, index))
    quad = max(abs(check[n]/d[n]-1) for n in ("mass", "elastic_energy"))
    singular = np.linalg.svd(provider.matrix(solved["z"], ctheta), compute_uv=False)
    gates = {
        "newton": solved["status"] == "CONVERGED", "sign": p.real <= 0 and p.imag > 0,
        "local_radius": abs(p-1j*omega) <= policy["complex_neighbour_radius_fraction"]*gap,
        "joint": max(d["joint_scaled"]) <= policy["boundary_tol"],
        "clamp": d["clamp_scaled"] <= policy["boundary_tol"], "ode": d["ode_scaled"] <= policy["ode_tol"],
        "mass": abs(d["mass"]-1) <= policy["mass_normalization_tol"],
        "energy": energy <= policy["energy_tol"], "alpha_identity": alpha_identity <= policy["energy_tol"],
        "quadrature": quad <= policy["energy_tol"], "prediction": prediction <= policy["complex_predictor_relative_tol"],
        "identity": macs[index] >= policy["shape_mac_min"] and margin >= policy["shape_margin_min"],
        "singularity": singular[-1]/singular[0] <= policy["boundary_tol"],
        "simple": singular[0]/singular[-2] <= policy["nonzero_condition_max"]}
    record.update(status="PASS" if all(gates.values()) else "NUMERICAL_UNRESOLVED", gates=gates,
        coefficients=coefficients, diagnostics=d, alpha=-p.real, omega_d=p.imag,
        zeta=-p.real/abs(p), alpha_predictor=alpha_pred, predictor_relative_error=prediction,
        energy_identity_residual=energy_residual, energy_identity_relative=energy,
        alpha_identity_relative=alpha_identity, quadrature_relative=quad,
        mass_MAC_with_elastic=macs[index], all_elastic_MAC=macs, identity_margin=margin,
        singular_ratio=singular[-1]/singular[0], nonzero_condition=singular[0]/singular[-2],
        basis_reference_rates=provider.rates, basis_pivots=provider.pivots, basis_scales=provider.scales,
        initial_mass=initial_mass)
    return record, d
