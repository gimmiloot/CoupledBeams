"""NSL-1: straight cross-ply FSDT beam, with extension/bending coupling.

Independent of symmetric RLB production physics. Layer integration is imported,
never reimplemented. State (u,w,psi,N,Q,M), Ux=u+z*psi, gamma=w'+psi.
No angular joint, spring, dashpot or damped-spectrum implementation.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
import numpy as np
from scipy.linalg import eigh, expm
from scipy.optimize import brentq

from scripts.lib import reddy_symmetric_laminated_beam as rlb
from scripts.analysis.laminated_beams.audit_reddy_table_4_3_3_discrepancies import (
    shifted_legendre_essential_basis,
)

SOURCE_BACHOO = "SOURCE_BACHOO"
PROJECT_REDDY_COUPLED = "PROJECT_REDDY_COUPLED"
REVERSE = np.diag([1., -1., -1., 1., -1., -1.])


def positive_definite(matrix):
    """Check symmetry and SPD after dimensionally consistent diagonal scaling."""
    a = np.asarray(matrix, dtype=float)
    if not np.isfinite(a).all() or np.any(np.diag(a) <= 0):
        raise ValueError("Nonfinite or nonpositive constitutive/mass diagonal")
    scale = np.sqrt(np.diag(a))
    balanced = a/scale[:, None]/scale[None, :]
    symmetry = float(np.max(abs(balanced-balanced.T)))
    if symmetry > 1e-11:
        raise ValueError("Constitutive/mass matrix is not symmetric")
    eigenvalues = np.linalg.eigvalsh(balanced)
    if eigenvalues[0] <= 0:
        raise ValueError("Constitutive/mass matrix is not positive definite")
    return {"symmetry": symmetry, "scaled_min_eigenvalue": float(eigenvalues[0])}


@dataclass(frozen=True)
class Properties:
    A: float
    B: float
    D: float
    S: float
    m: float
    j1: float
    j2: float
    variant: str
    K: float
    width: float
    provenance: str

    def __post_init__(self):
        if self.variant not in (SOURCE_BACHOO, PROJECT_REDDY_COUPLED):
            raise ValueError("Unknown NSL constitutive variant")
        if not all(math.isfinite(v) for v in (self.A, self.B, self.D, self.S,
                                              self.m, self.j1, self.j2, self.K, self.width)):
            raise ValueError("Nonfinite NSL properties")
        if min(self.S, self.K, self.width) <= 0:
            raise ValueError("S, K and width must be positive")
        positive_definite([[self.A, self.B], [self.B, self.D]])
        positive_definite(self.mass_matrix())

    def mass_matrix(self):
        return np.array([[self.m, 0., self.j1], [0., self.m, 0.], [self.j1, 0., self.j2]])

    def energy_matrices(self):
        gradient = np.array([[self.A, 0., self.B], [0., self.S, 0.], [self.B, 0., self.D]])
        mixed = np.zeros((3, 3)); mixed[1, 2] = self.S
        value = np.diag([0., 0., self.S])
        return gradient, mixed, value, self.mass_matrix()


def _supported(section):
    # Even a balanced off-axis stack may hide layer-level spatial couplings.
    # NSL-1 accepts only principal-axis 0/90 plies; it does not discard them.
    for ply in section.plies:
        if abs(math.sin(2*math.radians(ply.angle_deg))) > 1e-12:
            raise ValueError("NSL-1 supports cross-ply/coaxial 0/90 layers only; off-axis spatial coupling is unsupported")
    for name in ("A", "B", "D"):
        a = getattr(section, name)
        scale = np.linalg.norm(section.A)*section.thickness**{"A": 0, "B": 1, "D": 2}[name]
        if max(abs(a[0, 2]), abs(a[1, 2]))/scale > 1e-12:
            raise ValueError("Unsupported extension/bending-twist coupling in "+name)
    if abs(section.shear[0, 1])/np.linalg.norm(section.shear) > 1e-12:
        raise ValueError("Unsupported coupled transverse shear")


def reduce_laminate(section, *, width, K, variant, provenance):
    """Two distinct reductions of the same existing integrated section.

    PROJECT: retained ABD coordinates (0,3); zero secondary resultants
    Ny,Nxy,My,Mxy, with no additional inertia assigned to eliminated strains.
    SOURCE: direct b*A11,b*B11,b*D11 and K*b*As55 (cross-ply only).
    """
    _supported(section)
    if not math.isfinite(width) or width <= 0 or not math.isfinite(K) or K <= 0:
        raise ValueError("width and K must be positive")
    abd = width*np.block([[section.A, section.B], [section.B, section.D]])
    full_check = positive_definite(abd)
    aa, bb = [0, 3], [1, 2, 4, 5]
    primary = abd[np.ix_(aa, aa)]
    coupling = abd[np.ix_(aa, bb)]
    secondary = abd[np.ix_(bb, bb)]
    eliminated = -np.linalg.solve(secondary, coupling.T)
    schur = primary+coupling@eliminated
    columns = np.linalg.solve(abd, np.eye(6)[:, aa])
    compliance_route = np.linalg.solve(columns[aa], np.eye(2))
    scale = np.sqrt(np.diag(schur))
    route_error = float(np.max(abs(schur-compliance_route)/scale[:, None]/scale[None, :]))
    if route_error > 1e-11:
        raise ArithmeticError("Coupled compliance/Schur mismatch")
    if variant == SOURCE_BACHOO:
        effective = primary
        shear = K*width*section.shear[1, 1]
    elif variant == PROJECT_REDDY_COUPLED:
        effective = schur
        shear = K*width*rlb.effective_stiffness_reduction(section.shear, 1).value
    else:
        raise ValueError("Explicit SOURCE_BACHOO or PROJECT_REDDY_COUPLED required")
    properties = Properties(float(effective[0, 0]), float(effective[0, 1]),
        float(effective[1, 1]), float(shear), width*section.I0, width*section.I1,
        width*section.I2, variant, K, width, provenance)
    audit = {"variant": variant, "ABD": abd, "primary_indices": aa, "secondary_indices": bb,
        "z_interfaces": section.z_interfaces,
        "plies_bottom_to_top": [asdict(ply) for ply in section.plies],
        "full_SPD": full_check, "effective_SPD": positive_definite(effective),
        "mass_SPD": positive_definite(properties.mass_matrix()),
        "coupled_schur": schur, "compliance_route": compliance_route,
        "reduction_relative_error": route_error, "eliminated_strain_map": eliminated,
        "effective": effective, "properties": asdict(properties),
        "units": {"A": "force", "B": "force*length", "D": "force*length^2",
                  "S": "force", "m": "mass/length", "j1": "mass", "j2": "mass*length"},
        "coordinate": "z_project increasing bottom-to-top, -h/2 to h/2"}
    return properties, audit


def source_bachoo_properties(section, *, width, K, provenance):
    return reduce_laminate(section, width=width, K=K, variant=SOURCE_BACHOO, provenance=provenance)


def project_reddy_coupled_properties(section, *, width, K, provenance):
    return reduce_laminate(section, width=width, K=K, variant=PROJECT_REDDY_COUPLED, provenance=provenance)


def state_matrix(p, properties):
    """Hamilton equations derived from the quadratic energy Hessians.

    f=G*q'+C*q; f'=C.T*q'+(D+p^2*M)*q. No inverse of the full ABD is formed.
    """
    g, c, d, mass = properties.energy_matrices()
    gi = np.linalg.solve(g, np.eye(3))
    gc = np.linalg.solve(g, c)
    return np.block([[-gc, gi], [d+p*p*mass-c.T@gc, c.T@gi]])


def control_state_matrix(p, a):
    """Independent displayed NSL state formula for regression, same state order."""
    delta = a.A*a.D-a.B*a.B
    return np.array([
        [0, 0, 0, a.D/delta, 0, -a.B/delta],
        [0, 0, -1, 0, 1/a.S, 0],
        [0, 0, 0, -a.B/delta, 0, a.A/delta],
        [a.m*p*p, 0, a.j1*p*p, 0, 0, 0],
        [0, a.m*p*p, 0, 0, 0, 0],
        [a.j1*p*p, 0, a.j2*p*p, 0, 1, 0]])


def cc_basis(order, xi):
    # This prior RLB single-beam basis imposes both endpoint VALUES exactly.
    return shifted_legendre_essential_basis(order, np.asarray(xi), "CC", "w")


def ritz_matrices(properties, length, order, quadrature=64):
    """Independent direct assembly of V,T; does not call a state/transfer API."""
    nodes, weights = np.polynomial.legendre.leggauss(quadrature)
    xi, weights = (nodes+1)/2, weights*length/2
    values, derivative = cc_basis(order, xi)
    derivative /= length
    v = (values*weights)@values.T
    d = (derivative*weights)@derivative.T
    h = (derivative*weights)@values.T
    a = properties
    zero = np.zeros_like(v)
    stiffness = np.block([[a.A*d, zero, a.B*d],
                          [zero, a.S*d, a.S*h],
                          [a.B*d, a.S*h.T, a.D*d+a.S*v]])
    mass = np.block([[a.m*v, zero, a.j1*v],
                     [zero, a.m*v, zero], [a.j1*v, zero, a.j2*v]])
    return stiffness, mass


def ritz_spectrum(properties, length, order, quadrature=64):
    stiffness, mass = ritz_matrices(properties, length, order, quadrature)
    positive_definite(stiffness); positive_definite(mass)
    # Constant nonsingular coordinate scaling; energies and eigenvalues unchanged.
    scales = 1/np.sqrt(np.diag(mass))
    kb, mb = stiffness*scales[:, None]*scales[None, :], mass*scales[:, None]*scales[None, :]
    lam, vectors = eigh(kb, mb, subset_by_index=[0, 6])
    if min(lam) <= 0:
        raise ArithmeticError("Nonpositive Ritz eigenvalue")
    coefficients = scales[:, None]*vectors
    residual = stiffness@coefficients-(mass@coefficients)*lam
    relative = np.linalg.norm(residual, axis=0)/(np.linalg.norm(stiffness@coefficients, axis=0)+np.linalg.norm((mass@coefficients)*lam, axis=0))
    return np.sqrt(lam), coefficients, {"eigen_residual": relative,
        "mass_gram_error": float(np.max(abs(coefficients.T@mass@coefficients-np.eye(7))))}


def ritz_fields(coefficients, length, x):
    order = len(coefficients)//3
    v, d = cc_basis(order, np.asarray(x)/length)
    return v.T@coefficients.reshape(3, order).T, d.T@coefficients.reshape(3, order).T/length


class SingleBeam:
    """Bounded NSL-1 propagation/count adapter for the one clamped beam."""
    def __init__(self, properties, length):
        self.a, self.length = properties, float(length)
        self.scale = np.array([length, length, 1., properties.D/length**2,
                               properties.D/length**2, properties.D/length])

    def balanced(self, omega):
        return state_matrix(1j*omega, self.a).real*self.scale[None, :]/self.scale[:, None]

    def basis(self, omega, x, derivative=False):
        h = self.balanced(omega)
        values = np.array([expm(h*(xx-self.length/2)) for xx in np.atleast_1d(x)])
        if derivative:
            values = h@values
        return self.scale[None, :, None]*values

    def boundary(self, omega):
        end = self.basis(omega, [0., self.length])[:, :3, :].reshape(6, 6)
        return end/np.linalg.norm(end, axis=1)[:, None]

    def segment_bound(self, segments):
        """Young + vector Poincare lower bound for every auxiliary DD segment."""
        a, ell = self.a, self.length/segments
        gradient = np.array([[a.A, 0., a.B], [0., a.S/2, 0.], [a.B, 0., a.D]])
        mass = a.mass_matrix()
        derivative_min = float(eigh(gradient, mass, eigvals_only=True)[0])
        negative_max = a.S*np.linalg.solve(mass, np.eye(3)[:, 2])[2]
        lower_squared = derivative_min*(np.pi/ell)**2-negative_max
        return {"segments": segments, "length": ell, "derivative_min": derivative_min,
                "negative_max": negative_max, "omega_squared_lower": float(lower_squared)}

    def count_matrix(self, omega, segments):
        """Exact DtN condensation, NOT a finite-element discretization."""
        t = expm(self.balanced(omega)*(self.length/segments))
        a, b, c, d = t[:3, :3], t[:3, 3:], t[3:, :3], t[3:, 3:]
        ba, bi = np.linalg.solve(b, a), np.linalg.solve(b, np.eye(3))
        # The balancing is work-dual up to one common positive energy factor.
        local = np.block([[ba, -bi], [c-d@ba, d@bi]])
        size = 3*(segments-1)
        matrix = np.zeros((size, size))
        for i in range(segments-1):
            block = slice(3*i, 3*i+3)
            matrix[block, block] = local[:3, :3]+local[3:, 3:]
            if i < segments-2:
                next_block = slice(3*i+3, 3*i+6)
                matrix[block, next_block] = local[:3, 3:]
                matrix[next_block, block] = local[3:, :3]
        # Congruence makes the inertia computation insensitive to trace units.
        scaling = 1/np.sqrt(abs(np.diag(matrix)))
        balanced = matrix*scaling[:, None]*scaling[None, :]
        skew = np.max(abs(balanced-balanced.T))/max(np.max(abs(balanced)), 1e-30)
        eigenvalues = np.linalg.eigvalsh((balanced+balanced.T)/2)
        return {"count": int(sum(eigenvalues < 0)), "symmetry": float(skew),
                "margin": float(min(abs(eigenvalues))/max(abs(eigenvalues))),
                "eigenvalues": eigenvalues, "segment_DtN": local}

    def green_identity(self, omega, segments):
        """Independent physical integral checks the signs/work scale of DtN."""
        ell=self.length/segments
        h=self.balanced(omega)
        t=expm(h*ell)
        z=np.array([.7,.2,-.3,.1,-.4,.5])
        f0=np.linalg.solve(t[:3,3:],z[3:]-t[:3,:3]@z[:3])
        initial=np.r_[z[:3],f0]
        nodes,weights=np.polynomial.legendre.leggauss(32)
        x,weights=(nodes+1)*ell/2,weights*ell/2
        y=np.array([self.scale*(expm(h*xx)@initial) for xx in x])
        dy=(y/self.scale)@h.T*self.scale
        q,dq=y[:,:3],dy[:,:3]
        a=self.a
        stiffness=weights@(a.A*dq[:,0]**2+2*a.B*dq[:,0]*dq[:,2]+a.D*dq[:,2]**2+a.S*(dq[:,1]+q[:,2])**2)
        mass=np.sum(weights*np.einsum("ni,ij,nj->n",q,a.mass_matrix(),q))
        end=t@initial
        work=(z[3:]@end[3:]-z[:3]@f0)*(a.D/self.length)
        return float(abs(work-(stiffness-omega**2*mass))/max(abs(work),abs(stiffness-omega**2*mass)))

    def prefix(self, policy, checkpoint):
        bound = self.segment_bound(policy["count_segments"])
        checkpoint["segment_coercivity"] = bound
        trial = 7*np.pi*np.sqrt(self.a.S/self.a.m)/self.length
        checkpoint["seventh_trial_bound"] = trial
        if not trial < policy["omega_upper"] or not bound["omega_squared_lower"] > policy["omega_upper"]**2:
            raise ArithmeticError("Preset spectral range/count coercivity not certified")
        samples, roots = checkpoint.setdefault("count_samples", []), checkpoint.setdefault("roots", [])
        def count(w):
            if len(samples) >= policy["count_query_budget"]:
                raise ArithmeticError("Count budget exhausted")
            d = self.count_matrix(w, policy["count_segments"])
            samples.append({"omega": w, **d})
            if d["symmetry"] > policy["count_symmetry_tol"] or d["margin"] < policy["count_margin_min"]:
                raise ArithmeticError("Unresolved exact-segment inertia")
            return d["count"]
        def isolate(lo, hi, nlo, nhi, depth=0):
            if nlo >= 7 or nhi == nlo:
                return
            if nhi < nlo or depth > policy["isolation_depth"]:
                raise ArithmeticError("Prefix isolation unresolved")
            if nhi-nlo == 1:
                f = lambda w: float(np.linalg.det(self.boundary(w)))
                if f(lo)*f(hi) >= 0:
                    raise ArithmeticError("Count-one interval lacks a simple determinant bracket")
                w, info = brentq(f, lo, hi, xtol=policy["root_xtol"], rtol=policy["root_rtol"], full_output=True)
                if not info.converged:
                    raise ArithmeticError("Brent did not converge")
                roots.append({"omega": w, "bracket": [lo, hi], "counts": [nlo, nhi], "iterations": info.iterations})
            else:
                middle = (lo+hi)/2
                nm = count(middle)
                isolate(lo, middle, nlo, nm, depth+1)
                isolate(middle, hi, nm, nhi, depth+1)
        lo, hi = policy["omega_lower"], policy["omega_upper"]
        nlo, nhi = count(lo), count(hi)
        if nlo != 0 or nhi < 7:
            raise ArithmeticError("Fixed range does not certify the requested prefix")
        isolate(lo, hi, nlo, nhi)
        if len(roots) != 7 or np.any(np.diff([x["omega"] for x in roots]) <= policy["root_xtol"]):
            raise ArithmeticError("Incomplete or duplicated prefix")
        checkpoint.update(status="PASS", lower_count=nlo, upper_count=nhi)
        return roots

    def evaluate(self, omega, coefficients, order):
        a = self.a
        x, weights = np.polynomial.legendre.leggauss(order)
        x, weights = (x+1)*self.length/2, weights*self.length/2
        y = self.basis(omega, x)@coefficients
        derivative = self.basis(omega, x, True)@coefficients
        q, dq = y[:, :3], derivative[:, :3]
        mass = float(np.sum(weights*np.einsum("ni,ij,nj->n", q, a.mass_matrix(), q)))
        density = a.A*dq[:, 0]**2+2*a.B*dq[:, 0]*dq[:, 2]+a.D*dq[:, 2]**2+a.S*(dq[:, 1]+q[:, 2])**2
        stiffness = float(weights@density)
        physical_rhs = y@control_state_matrix(1j*omega, a).real.T
        ode = float(np.max(abs(derivative-physical_rhs)/np.maximum(np.max(abs(derivative)+abs(physical_rhs), axis=0), 1e-30)))
        ends = self.basis(omega, [0., self.length])@coefficients
        bc = float(np.max(abs(ends[:, :3])/np.maximum(np.max(abs(q), axis=0), 1e-30)))
        return {"x": x, "weights": weights, "states": y, "derivatives": derivative,
                "mass": mass, "stiffness": stiffness, "energy_relative": abs(stiffness/(omega**2*mass)-1),
                "ode_scaled": ode, "boundary_scaled": bc, "end_states": ends}

    def mode(self, omega, policy):
        _, singular, vh = np.linalg.svd(self.boundary(omega))
        coefficient = vh[-1].copy()
        d = self.evaluate(omega, coefficient, policy["mode_quadrature"])
        coefficient /= math.sqrt(d["mass"])
        d = self.evaluate(omega, coefficient, policy["mode_quadrature"])
        fine = self.evaluate(omega, coefficient, policy["mode_quadrature_check"])
        d["quadrature_relative"] = max(abs(fine[n]/d[n]-1) for n in ("mass", "stiffness"))
        d["singular_ratio"] = float(singular[-1]/singular[0])
        d["nonzero_condition"] = float(singular[0]/singular[-2])
        return coefficient, d
