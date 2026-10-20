"""NSL-2 elastic angular joint; NSL-1 physics and RLB coordinates unchanged.

No constitutive reduction is duplicated here. The short-segment inertia
counter is the NSL-1 energy construction assembled on the coupled trace space.
Only real undamped frequencies are searched.
"""
from dataclasses import replace
import numpy as np
from scipy.linalg import block_diag, eigh, expm, null_space
from scripts.lib import nonsymmetric_laminated_beam as nsl
from scripts.lib.inplane_rotational_spring_eb import Joint, joint_matrix
from scripts.lib.reddy_symmetric_coupled_beams_ritz import (
    global_kinematic_map, shifted_legendre_clamped_basis,
)

REFLECTION = np.block([[np.zeros((6, 6)), nsl.REVERSE],
                       [nsl.REVERSE, np.zeros((6, 6))]])
CONTINUOUS = np.diag([-1., -1., 1., 1., 1., -1.])


def axis_transform(z0):
    """y*=T y; forces are work-dual, not independently rotated."""
    a = np.eye(3); a[0, 2] = z0
    return block_diag(a, np.linalg.solve(a.T, np.eye(3)))


def shift_properties(a, z0):
    return replace(a, B=a.B-z0*a.A, D=a.D-2*z0*a.B+z0*z0*a.A,
                   j1=a.j1-z0*a.m, j2=a.j2-2*z0*a.j1+z0*z0*a.m,
                   provenance=a.provenance+f"; reference shift {z0}")


def covariant_joint(beta, joint, offsets):
    """Pull back to the ORIGINAL physical center and moment about that point."""
    t = block_diag(*(axis_transform(z) for z in offsets))
    return np.linalg.solve(t.T, joint_matrix(beta, joint).T).T


def trace_maps(beta, joint):
    """q_i(L)=E_i z: common global translations, independent/common rotations."""
    size = 3 if joint.mode == "RIGID" else 4
    maps = []
    for arm in (1, 2):
        e = np.zeros((3, size))
        e[:2, :2] = global_kinematic_map(beta, arm)[:2, :2].T
        e[2, 2 if size == 3 else arm+1] = 1
        maps.append(e)
    jump = maps[0][2]-maps[1][2]
    return maps, jump


def arm_ritz(a, length, order, quadrature):
    """Direct energy assembly; independent of H, propagation and force rows."""
    x, weights = np.polynomial.legendre.leggauss(quadrature)
    basis = shifted_legendre_clamped_basis(order, (x+1)/2)
    v, d, weights = basis.values, basis.derivatives/length, weights*length/2
    vv, dd, dv = (v*weights)@v.T, (d*weights)@d.T, (d*weights)@v.T
    zero = np.zeros_like(vv)
    k = np.block([[a.A*dd, zero, a.B*dd], [zero, a.S*dd, a.S*dv],
                  [a.B*dd, a.S*dv.T, a.D*dd+a.S*vv]])
    m = np.block([[a.m*vv, zero, a.j1*vv], [zero, a.m*vv, zero],
                  [a.j1*vv, zero, a.j2*vv]])
    return k, m


def ritz_matrices(properties, lengths, beta, joint, order, quadrature=64):
    matrices = [arm_ritz(a, length, order, quadrature) for a, length in zip(properties, lengths)]
    k, m = (block_diag(*(v[j] for v in matrices)) for j in (0, 1))
    end = np.kron(np.eye(3), np.ones((1, order)))
    maps = [global_kinematic_map(beta, i) for i in (1, 2)]
    constraint = np.hstack([maps[0]@end, -maps[1]@end])
    jump = constraint[2].copy()
    if joint.mode != "RIGID":
        constraint = constraint[:2]
        k += joint.k_theta*np.outer(jump, jump)
    z = null_space(constraint)
    return k, m, z, constraint, jump


def ritz_spectrum(properties, lengths, beta, joint, order, quadrature=64):
    k, m, z, c, jump = ritz_matrices(properties, lengths, beta, joint, order, quadrature)
    kr, mr = z.T@k@z, z.T@m@z
    spd = {"K": nsl.positive_definite(kr), "M": nsl.positive_definite(mr)}
    scales = 1/np.sqrt(np.diag(mr))
    lam, vectors = eigh(kr*scales[:, None]*scales, mr*scales[:, None]*scales,
                        subset_by_index=[0, 6])
    reduced = scales[:, None]*vectors
    coefficients = z@reduced
    residual = kr@reduced-(mr@reduced)*lam
    return np.sqrt(lam), coefficients, {
        "SPD": spd, "constraint_rank": int(np.linalg.matrix_rank(c)),
        "constraint_residual": float(np.max(abs(c@coefficients))),
        "mass_gram_error": float(np.max(abs(coefficients.T@m@coefficients-np.eye(7)))),
        "eigen_residual": np.linalg.norm(residual, axis=0)/
        (np.linalg.norm(kr@reduced, axis=0)+np.linalg.norm((mr@reduced)*lam, axis=0)),
        "spring_derivatives": (jump@coefficients)**2}


def ritz_fields(coefficients, length, x):
    order = len(coefficients)//3
    b = shifted_legendre_clamped_basis(order, np.asarray(x)/length)
    return b.values.T@coefficients.reshape(3, order).T, b.derivatives.T@coefficients.reshape(3, order).T/length


class CoupledBeams(nsl.SingleBeam):
    """Reuse NSL-1 isolated-prefix algorithm with the *new* trace-space count."""
    def __init__(self, properties, lengths, beta, joint):
        self.properties, self.lengths = tuple(properties), tuple(lengths)
        self.arms = [nsl.SingleBeam(a, length) for a, length in zip(properties, lengths)]
        self.beta, self.joint = float(beta), joint
        # Inherited prefix trial uses 7*pi*sqrt(S/m)/length. We supply the
        # conservative value for 4 trial functions on each arm (>=7 total).
        self.a = properties[0]
        trial = max(4*np.pi*np.sqrt(a.S/a.m)/l for a, l in zip(properties, lengths))
        self.length = 7*np.pi*np.sqrt(self.a.S/self.a.m)/trial

    def boundary(self, omega, normalized=True, offsets=None):
        blocks = [a.basis(omega, [0., a.length]) for a in self.arms]
        end = block_diag(*(b[1] for b in blocks))
        clamps = block_diag(*(b[0, :3] for b in blocks))
        j = joint_matrix(self.beta, self.joint)
        if offsets is not None:
            # Independently exponentiate shifted physical H with the original
            # column coordinates transported by T. Do not just return J@basis.
            transformed, clamp_blocks = [], []
            for arm, z0 in zip(self.arms, offsets):
                t = axis_transform(z0); ai = shift_properties(arm.a, z0)
                h = nsl.state_matrix(1j*omega, ai).real
                scale = t@np.diag(arm.scale)
                ti = np.linalg.solve(t, np.eye(6))
                # Stable, fixed work-dual state scaling (not spectral-dependent).
                balanced = np.linalg.solve(scale, h@scale)
                b = [scale@expm(balanced*(x-arm.length/2)) for x in (0., arm.length)]
                transformed.append(b[1]); clamp_blocks.append(ti[:3]@b[0])
            end, clamps = block_diag(*transformed), block_diag(*clamp_blocks)
            j = covariant_joint(self.beta, self.joint, offsets)
        matrix = np.vstack([clamps, j@end])
        if normalized:
            matrix = matrix/np.linalg.norm(matrix, axis=1)[:, None]
        return matrix

    def segment_bound(self, segments):
        arms = [a.segment_bound(segments) for a in self.arms]
        return {"omega_squared_lower": min(a["omega_squared_lower"] for a in arms), "arms": arms}

    def count_matrix(self, omega, segments, schur_audit=False):
        maps, jump = trace_maps(self.beta, self.joint)
        nj = len(jump); ni = 6*(segments-1)
        matrix = np.zeros((ni+nj, ni+nj))
        # Common displacement units (L_ref=1); force units are work-dual.
        # All main/control arms have L=1 and the same D. Use physical DtN
        # here so the assembly also permits unequal arm scales.
        locals_ = []
        for i, arm in enumerate(self.arms):
            transfer = expm(arm.balanced(omega)*arm.length/segments)
            a, b, c, d = transfer[:3, :3], transfer[:3, 3:], transfer[3:, :3], transfer[3:, 3:]
            ba, bi = np.linalg.solve(b, a), np.linalg.solve(b, np.eye(3))
            local = np.block([[ba, -bi], [c-d@ba, d@bi]])
            # Convert balanced displacement/force into physical work.
            qs = np.tile(arm.scale[:3], 2)
            local = (arm.a.D/arm.length)*local/qs[:, None]/qs
            locals_.append(local)
            for segment in range(segments):
                indices, blocks = [], []
                for node in (segment, segment+1):
                    if node == 0:
                        indices.append([]); blocks.append(np.zeros((3, 0)))
                    elif node == segments:
                        indices.append(list(range(ni, ni+nj))); blocks.append(maps[i])
                    else:
                        start = i*3*(segments-1)+3*(node-1)
                        indices.append(list(range(start, start+3))); blocks.append(np.eye(3))
                for aidx in (0, 1):
                    for bidx in (0, 1):
                        matrix[np.ix_(indices[aidx], indices[bidx])] += (
                            blocks[aidx].T@local[3*aidx:3*aidx+3, 3*bidx:3*bidx+3]@blocks[bidx])
        if self.joint.mode != "RIGID":
            matrix[ni:, ni:] += self.joint.k_theta*np.outer(jump, jump)
        scales = 1/np.sqrt(abs(np.diag(matrix)))
        balanced = matrix*scales[:, None]*scales
        skew = np.max(abs(balanced-balanced.T))/max(np.max(abs(balanced)), 1e-30)
        ev = np.linalg.eigvalsh((balanced+balanced.T)/2)
        result = {"count": int(sum(ev < 0)), "symmetry": float(skew),
                  "margin": float(min(abs(ev))/max(abs(ev))),
                  "joint_trace_dimension": nj, "assembled_dimension": len(ev)}
        if schur_audit:
            ii, ij, jj = balanced[:ni, :ni], balanced[:ni, ni:], balanced[ni:, ni:]
            schur = jj-ij.T@np.linalg.solve(ii, ij)
            internal = np.linalg.eigvalsh((ii+ii.T)/2)
            joint = np.linalg.eigvalsh((schur+schur.T)/2)
            result.update(joint_DtN=schur, internal_inertia=int(sum(internal < 0)),
                          joint_inertia=int(sum(joint < 0)),
                          schur_inertia_sum=int(sum(internal < 0)+sum(joint < 0)),
                          joint_DtN_symmetry=float(np.max(abs(schur-schur.T))/np.max(abs(schur))))
        return result

    def evaluate(self, omega, coefficients, order):
        arms = [a.evaluate(omega, coefficient, order) for a, coefficient in zip(self.arms, coefficients.reshape(2, 6))]
        end = np.concatenate([d["end_states"][1] for d in arms])
        peaks = np.concatenate([np.maximum(np.max(abs(d["states"]), axis=0), 1e-30) for d in arms])
        j = joint_matrix(self.beta, self.joint)
        jr = j@end
        # Physical row residual scaled by the largest participating state
        # over the arms, not only a nearly-zero endpoint value.
        joint_residual = float(np.max(abs(jr)/np.maximum(abs(j)@peaks, 1e-30)))
        bc = max(float(np.max(abs(d["end_states"][0, :3])/peaks[6*i:6*i+3])) for i, d in enumerate(arms))
        jump = end[2]-end[8]
        spring = 0 if self.joint.mode == "RIGID" else self.joint.k_theta*jump*jump
        components = dict(axial=0., mixed=0., bending=0., shear=0., spring=float(spring))
        for a, d in zip(self.properties, arms):
            q, dq, weights = d["states"][:, :3], d["derivatives"][:, :3], d["weights"]
            for key, value in (
                ("axial", a.A*dq[:, 0]**2), ("mixed", 2*a.B*dq[:, 0]*dq[:, 2]),
                ("bending", a.D*dq[:, 2]**2), ("shear", a.S*(dq[:, 1]+q[:, 2])**2)):
                components[key] += float(weights@value)
        mass = sum(d["mass"] for d in arms); stiffness = sum(components.values())
        return {"arms": arms, "mass": mass, "stiffness": stiffness, "stiffness_components": components,
                "delta_psi": float(jump), "spring_energy": float(spring/2),
                "energy_relative": float(abs(stiffness/(omega*omega*mass)-1)),
                "boundary_scaled": bc, "joint_scaled": joint_residual,
                "joint_physical_residual": jr, "ode_scaled": max(d["ode_scaled"] for d in arms)}

    # SingleBeam.mode uses only boundary/evaluate, therefore also mass-normalizes
    # both arms with their complete NSL mass, including j1 cross inertia.


def inner(left, right, properties):
    """Physical mass pairing on common arm grids, including cross inertia."""
    return float(sum(np.sum(a["weights"]*np.einsum("ni,ij,nj->n",
        a["states"][:, :3], p.mass_matrix(), b["states"][:, :3]))
        for a, b, p in zip(left["arms"], right["arms"], properties)))


def reflection_metrics(mode, properties):
    flipped = {"arms": []}
    for d in reversed(mode["arms"]):
        flipped["arms"].append({**d, "states": d["states"]@nsl.REVERSE})
    value = inner(mode, flipped, properties)/mode["mass"]
    parity = 1 if value >= 0 else -1
    return {"reflection_pairing": value, "candidate_parity": parity,
            "relative_error": float(np.sqrt(max(0., 2-2*abs(value))))}
