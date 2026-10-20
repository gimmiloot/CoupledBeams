"""NSL-3A complex adapter: unchanged NSL H, centered basis and physical joint.

Only kv.correct/right_null/joint_matrix are shared with the historical KV
implementation. No EB/RLB internal operator is used for an NSL arm.
Row and column scales are frozen at the elastic seed, as required by NSL-2R.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import block_diag, expm, expm_frechet
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as elastic
from scripts.lib import inplane_kelvin_voigt as kv


def state_derivative(p, properties):
    """dH/dp from the NSL energy mass Hessian, including both j1 entries."""
    h = np.zeros((6, 6), complex)
    h[3:, :3] = 2*p*properties.mass_matrix()
    return h


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.linalg.norm(b), 1e-300))


class Provider:
    """Full 12x12 analytic boundary pencil; centered coefficients are retained."""
    def __init__(self, properties, lengths, beta, k, c, t_ref, omega0):
        if k < 0 or c < 0:
            raise ValueError("Passive k>=0,c>=0 required")
        if any(a.variant != "PROJECT_REDDY_COUPLED" for a in properties):
            raise ValueError("NSL-3A requires PROJECT_REDDY_COUPLED")
        self.properties, self.lengths = tuple(properties), tuple(lengths)
        self.arms = [single.SingleBeam(a, length) for a, length in zip(properties, lengths)]
        self.beta, self.k, self.c, self.t_ref = beta, k, c, t_ref
        self.calls = dict(B=0, B_z=0, expm=0, frechet=0, recoveries=0)
        raw, _ = self.raw(1j*omega0, c_override=0.)
        self.rows = np.linalg.norm(raw, axis=1)
        self.columns = np.linalg.norm(raw/self.rows[:, None], axis=0)
        if min(self.rows.min(), self.columns.min()) <= 0:
            raise ArithmeticError("Singular scaling")

    def basis(self, p, arm, x, derivative=False, offset=None):
        a, length = arm.a, arm.length
        scale = np.diag(arm.scale)
        if offset is not None:
            a = elastic.shift_properties(a, offset)
            scale = elastic.axis_transform(offset)@scale
        h = np.linalg.solve(scale, single.state_matrix(p, a)@scale)
        hp = np.linalg.solve(scale, state_derivative(p, a)@scale)
        values, derivatives = [], []
        for xx in np.atleast_1d(x):
            distance = xx-length/2
            self.calls["expm"] += 1
            values.append(scale@expm(h*distance))
            if derivative:
                self.calls["frechet"] += 1
                derivatives.append(scale@expm_frechet(
                    h*distance, hp*distance, compute_expm=False))
        return np.asarray(values), np.asarray(derivatives) if derivative else None

    def raw(self, p, derivative=False, offsets=None, c_override=None):
        c = self.c if c_override is None else c_override
        transforms = [np.eye(6)]*2 if offsets is None else [
            elastic.axis_transform(z) for z in offsets]
        inv = [np.linalg.solve(t, np.eye(6)) for t in transforms]
        values = [self.basis(p, arm, [0., arm.length], derivative,
                            None if offsets is None else offsets[i])
                  for i, arm in enumerate(self.arms)]
        j = kv.joint_matrix(p, self.beta, self.k, c)@block_diag(*inv)
        end = block_diag(*(v[0][1] for v in values))
        clamps = block_diag(*(ti[:3]@v[0][0] for ti, v in zip(inv, values)))
        b = np.vstack([clamps, j@end])
        if not derivative:
            return b, None
        jp = kv.joint_matrix(p, self.beta, self.k, c, derivative=True)@block_diag(*inv)
        endp = block_diag(*(v[1][1] for v in values))
        clamp_p = block_diag(*(ti[:3]@v[1][0] for ti, v in zip(inv, values)))
        return b, np.vstack([clamp_p, jp@end+j@endp])

    def matrices(self, z, *, derivative=False):
        self.calls["B"] += 1
        self.calls["B_z"] += int(derivative)
        b, bp = self.raw(z/self.t_ref, derivative)
        scales = self.rows[:, None]*self.columns[None, :]
        return b/scales, bp/scales/self.t_ref if derivative else None

    def null(self, z):
        """Hermitian right SVD null; undo fixed column scaling explicitly."""
        b, _ = self.matrices(z)
        v = kv.right_null(b)
        return v/self.columns, v

    def evaluate(self, p, coefficients, order):
        nodes, weights = np.polynomial.legendre.leggauss(order)
        fields = []
        components = dict(axial=0., mixed=0., bending=0., shear=0.)
        mass_components = dict(translation=0., mixed=0., rotation=0.)
        for arm, a, coef in zip(self.arms, self.properties, np.asarray(coefficients).reshape(2, 6)):
            x, weight = (nodes+1)*arm.length/2, weights*arm.length/2
            basis, _ = self.basis(p, arm, x)
            y = basis@coef
            # Independently displayed balance/constitutive matrix audits H.
            dy = (single.state_matrix(p, a)@basis)@coef
            control = y@single.control_state_matrix(p, a).T
            end_basis, _ = self.basis(p, arm, [0., arm.length])
            end = end_basis@coef
            q, dq = y[:, :3], dy[:, :3]
            for key, val in (
                ("axial", a.A*abs(dq[:, 0])**2),
                ("mixed", 2*a.B*np.real(dq[:, 0].conj()*dq[:, 2])),
                ("bending", a.D*abs(dq[:, 2])**2),
                ("shear", a.S*abs(dq[:, 1]+q[:, 2])**2)):
                components[key] += float(weight@val)
            for key, val in (
                ("translation", a.m*(abs(q[:, 0])**2+abs(q[:, 1])**2)),
                ("mixed", 2*a.j1*np.real(q[:, 0].conj()*q[:, 2])),
                ("rotation", a.j2*abs(q[:, 2])**2)):
                mass_components[key] += float(weight@val)
            ode_rows = np.max(abs(dy-control), axis=0)/np.maximum(
                np.max(abs(dy), axis=0)+np.max(abs(control), axis=0), 1e-30)
            fields.append(dict(x=x, weights=weight, states=y, derivatives=dy,
                               end_states=end, ode_rows=ode_rows))
        end = np.concatenate([f["end_states"][1] for f in fields])
        peaks = np.concatenate([np.maximum(np.max(abs(f["states"]), axis=0), 1e-30) for f in fields])
        j = kv.joint_matrix(p, self.beta, self.k, self.c)
        raw_joint = kv.scalar_conditions(end, p, self.beta, self.k, self.c)
        joint_rows = abs(raw_joint)/np.maximum(abs(j)@peaks, 1e-30)
        # Also verify the second independent moment law, not just moment sum.
        law2 = end[11]-(self.k+self.c*p)*(end[2]-end[8])
        law2_scale = peaks[11]+abs(self.k+self.c*p)*(peaks[2]+peaks[8])
        clamp_rows = np.concatenate([abs(f["end_states"][0, :3])/peaks[6*i:6*i+3]
                                    for i, f in enumerate(fields)])
        delta = end[2]-end[8]
        components["spring"] = float(self.k*abs(delta)**2)
        m, stiffness = sum(mass_components.values()), sum(components.values())
        dissipation = float(self.c*abs(delta)**2)
        identity = p*p*m+p*dissipation+stiffness
        energy = float(abs(identity)/(abs(p*p*m)+abs(p*dissipation)+abs(stiffness)))
        return dict(arms=fields, mass=m, mass_components=mass_components,
                    stiffness=stiffness, stiffness_components=components,
                    delta_psi=delta, C_phi=dissipation, energy_identity=identity,
                    energy_residual=energy, alpha_energy=dissipation/(2*m),
                    clamp_rows=clamp_rows, joint_rows=joint_rows,
                    joint_physical=raw_joint, second_moment_physical=law2,
                    second_moment_scaled=float(abs(law2)/max(law2_scale, 1e-30)),
                    ode_residual=float(max(max(f["ode_rows"]) for f in fields)))

    def recover(self, z, order, check_order):
        self.calls["recoveries"] += 1
        coefficients, _ = self.null(z)
        first = self.evaluate(z/self.t_ref, coefficients, order)
        coefficients /= np.sqrt(first["mass"])
        data = self.evaluate(z/self.t_ref, coefficients, order)
        fine = self.evaluate(z/self.t_ref, coefficients, check_order)
        data["quadrature_relative"] = max(abs(fine[k]/data[k]-1) for k in ("mass", "stiffness"))
        return coefficients, data


def pairing(first, second, properties):
    return sum(np.sum(a["weights"]*np.einsum("ni,ij,nj->n",
        a["states"][:, :3].conj(), p.mass_matrix(), b["states"][:, :3]))
        for a, b, p in zip(first["arms"], second["arms"], properties))


def mac(first, second, properties):
    return float(abs(pairing(first, second, properties))**2/(first["mass"]*second["mass"]))


def reflection_rows(beta):
    """B S = W B in physical rows; derived from the six scalar NSL conditions."""
    c, s = np.cos(beta), np.sin(beta)
    w = np.zeros((6, 6))
    w[:2, :2] = [[c, -s], [-s, -c]]
    w[2, [2, 5]] = [1., -1.]
    w[3:5, 3:5] = [[-c, s], [s, c]]
    w[5, 5] = -1.
    rq = single.REVERSE[:3, :3]
    wc = np.block([[np.zeros((3, 3)), rq], [rq, np.zeros((3, 3))]])
    return block_diag(wc, w)


def physical_gates(data, policy):
    return dict(clamps=float(max(data["clamp_rows"])) <= policy["physical_tol"],
        joint=float(max(data["joint_rows"])) <= policy["physical_tol"],
        second_moment=data["second_moment_scaled"] <= policy["physical_tol"],
        compatibility=float(max(data["joint_rows"][:2])) <= policy["compatibility_tol"],
        ode=data["ode_residual"] <= policy["ode_tol"],
        energy=data["energy_residual"] <= policy["energy_tol"],
        mass=abs(data["mass"]-1) <= policy["mass_tol"],
        quadrature=data.get("quadrature_relative", 0.) <= policy["quadrature_tol"])
