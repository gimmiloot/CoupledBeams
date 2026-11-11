"""CS boundary adapter and immutable CC verification for square partial spectra.

Reuses the accepted finite Timoshenko basis, root locator, modal recovery,
harmonic matrix and SS count from mindlin_herrmann_longitudinal. The private
module copy changes ONLY the boundary selector; the original module and all
coupled solvers remain untouched. No Euler--Bernoulli frequencies are used.
"""
from __future__ import annotations

import importlib.util
import math
import sys
import numpy as np
from scipy.linalg import expm
from scripts.lib import mindlin_herrmann_longitudinal as mh

VERSION = "tim-partial-cs-boundary-adapter-v1"


def cs_boundary_matrix(model, length, omega, block="timoshenko"):
    if block != "timoshenko":
        raise ValueError("CS adapter supports the Timoshenko bending block only")
    ends = mh.finite_state_basis(model, length, omega, [0., length], block)
    matrix = np.vstack((ends[0, 0], ends[0, 1], ends[1, 0], ends[1, 3]))
    norms = np.linalg.norm(matrix, axis=1)
    if np.any(norms == 0):
        raise ArithmeticError("Zero CS boundary row")
    return matrix/norms[:, None]


_spec = importlib.util.spec_from_file_location(__name__+"._cs", mh.__file__)
_cs = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = _cs
_spec.loader.exec_module(_cs)
_cs.finite_boundary_matrix = cs_boundary_matrix


def boundary_matrix(model, length, omega, bc):
    if bc == "CS":
        return cs_boundary_matrix(model, length, omega)
    if bc == "CC":
        return mh.finite_boundary_matrix(model, length, omega, "timoshenko")
    raise ValueError("Expected CS or CC")


def qr_boundary(model, length, omega, bc, exponent_cap=1.):
    """Existing short-step expm/positive-QR algorithm, with CS/CC end rows.

Unlike the closed-form finite basis, this evaluates the canonical first-order
state matrix. It checks a supplied frequency and never searches for a root.
"""
    if bc not in ("CS", "CC"):
        raise ValueError("Expected CS or CC")
    matrix = mh.harmonic_state_matrix(model, omega, "timoshenko")
    rate = float(np.max(np.abs(np.linalg.eigvals(matrix))))
    p = model.coefficients
    scales = np.array([1., length, 1/(p["S"]*rate), length/(p["B"]*rate)])
    balanced = matrix*scales[:, None]/scales[None, :]
    steps = max(1, math.ceil(rate*length/exponent_cap))
    if steps > 512:
        raise ArithmeticError("Independent shooting step budget exceeded")
    step = expm(balanced*length/steps)
    frame = np.vstack((np.zeros((2, 2)), np.eye(2)))
    for _ in range(steps):
        frame, triangular = np.linalg.qr(step@frame, mode="reduced")
        frame *= np.where(np.diag(triangular) >= 0, 1., -1.)[None, :]
    right = frame[[0, 1] if bc == "CC" else [0, 3]]
    singular = np.linalg.svd(right, compute_uv=False)
    return dict(singular_ratio=float(singular[-1]/singular[0]), steps=steps,
                exponent_cap=exponent_cap, method="canonical state expm / positive QR")


def recover(model, length, omega, bc, order):
    engine = _cs if bc == "CS" else mh
    mode = engine.finite_mode(model, length, omega, "timoshenko", order)
    coef = mode["coefficients"]
    mode["gradients"] = mh.finite_state_basis(model, length, omega, mode["points"], "timoshenko", 1)@coef
    mode["ends"] = mh.finite_state_basis(model, length, omega, [0., length], "timoshenko")@coef
    if bc == "CS":
        # The reused modal routine's CC boundary diagnostic is replaced by the
        # actual CS residual; all its ODE/energy/mass calculations are retained.
        y, ends = mode["values"], mode["ends"]
        qpeak = max(float(np.max(abs(y[:, 0]))), length*float(np.max(abs(y[:, 1]))))
        mpeak = float(np.max(abs(y[:, 3])))
        mode["diagnostics"]["boundary_scaled_residual"] = max(
            abs(ends[0, 0])/qpeak, length*abs(ends[0, 1])/qpeak,
            abs(ends[1, 0])/qpeak, abs(ends[1, 3])/mpeak)
    return mode


def verify_modes(model, length, roots, bc, policy):
    """Verify supplied roots; no CC root search or coupled-spectrum evaluation."""
    modes, rows = [], []
    for rank, root in enumerate(roots, 1):
        omega = root["omega"]
        low = recover(model, length, omega, bc, policy["quadrature_order"])
        high = recover(model, length, omega, bc, policy["quadrature_check_order"])
        s = np.linalg.svd(boundary_matrix(model, length, omega, bc), compute_uv=False)
        nullity = int(np.sum(s <= policy["nullity_tol"]*s[0]))
        qr = qr_boundary(model, length, omega, bc)
        qr_fine = qr_boundary(model, length, omega, bc, .5)
        d = dict(low["diagnostics"])
        d.update(rank=rank, omega=omega, singular_ratio=float(s[-1]/s[0]),
                 nonzero_condition=float(s[0]/s[-2]), nullity=nullity,
                 qr_singular_ratio=qr["singular_ratio"], qr_fine_singular_ratio=qr_fine["singular_ratio"],
                 quadrature_mass_relative_change=abs(high["diagnostics"]["mass_before_normalization"]/d["mass_before_normalization"]-1),
                 quadrature_energy_relative_change=abs(high["diagnostics"]["energy_omega_squared"]/d["energy_omega_squared"]-1))
        gates = dict(boundary=d["boundary_scaled_residual"]<=policy["boundary_tol"],
            ode=d["equation_scaled_residual"]<=policy["ode_tol"],
            energy=max(d["energy_relative_error"],high["diagnostics"]["energy_relative_error"])<=policy["energy_tol"],
            singular=d["singular_ratio"]<=policy["singular_tol"], simple=nullity==1,
            condition=d["nonzero_condition"]<=policy["nonzero_condition_max"],
            independent=max(qr["singular_ratio"],qr_fine["singular_ratio"])<=policy["singular_tol"],
            quadrature=d["quadrature_mass_relative_change"]<=policy["mass_gram_tol"] and d["quadrature_energy_relative_change"]<=policy["energy_tol"])
        d["gates"] = gates
        if not all(gates.values()):
            raise ArithmeticError(f"{bc} L={length}, rank={rank}: {d}")
        modes.append(low); rows.append(d)
    if modes:
        p=model.coefficients
        vectors=np.array([np.r_[np.sqrt(m["weights"]*p["m"])*m["values"][:,0],
                                np.sqrt(m["weights"]*p["r"])*m["values"][:,1]] for m in modes])
        error=float(np.max(abs(vectors@vectors.T-np.eye(len(modes)))))
        if error>policy["mass_gram_tol"]:
            raise ArithmeticError(f"{bc} L={length}: mass Gram {error}")
    else:
        error=0.
    return modes, dict(status="PASS",gram_max_error=error,modes=rows)


def solve_cs(model, length, cc, omega_top, policy):
    """Saturate the existing SS min-max upper count, with at most one band guard.

    V_CC subset V_CS subset V_SS. Thus lambda_SS,n <= lambda_CS,n <=
    lambda_CC,n for this ONE isolated beam. No coupled-system bound is claimed.
    The copied locator is unchanged; only its boundary matrix is CS.
    """
    lower=.001*math.pi
    lower_bound=mh.finite_count_upper_bound(model,length,lower,"timoshenko")
    if lower_bound["upper_count"]:
        raise ArithmeticError("Uncertified lower search cutoff")
    band_bound=mh.finite_count_upper_bound(model,length,omega_top,"timoshenko")
    count=band_bound["upper_count"]
    if count==0:
        return [],dict(status="PASS",lower_count=0,upper_count=0,band_upper_count=0,
                       upper=omega_top,method="SS relaxation excludes all roots in band",search_calls=0)
    if len(cc["poles"])<count:
        raise ArithmeticError("CC catalogue insufficient to bound CS search")
    upper=cc["poles"][count-1]*(1+policy["cluster_probe_relative"])
    bound=mh.finite_count_upper_bound(model,length,upper,"timoshenko")
    if bound["upper_count"]!=count:
        raise ArithmeticError("Chosen CS ceiling does not saturate the SS count")
    search_policy=dict(scan_intervals=policy["pole_scan_intervals"],root_xtol=policy["root_xtol"],root_rtol=policy["root_rtol"])
    roots,search=_cs.finite_roots(model,length,"timoshenko",lower,upper,search_policy)
    if len(roots)!=count or not all(r["converged"] for r in roots):
        raise ArithmeticError(f"CS count not saturated: L={length}, found {len(roots)}, bound {count}")
    return roots,dict(status="PASS",lower_count=0,upper_count=count,band_upper_count=count,
                     upper=upper,method="found modes saturate SS relaxation upper bound",search_calls=1,search=search,bound=bound)
