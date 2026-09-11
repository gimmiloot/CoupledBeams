"""Mass-normalized EB modes in the existing physical inward-coordinate state.

No circular-rod physics or frequency relabeling. Reflection exchanges arms
and reverses w, psi, Q and M. Its two exact classes can be resolved separately.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import block_diag, expm
from scipy.optimize import linear_sum_assignment

from scripts.lib import inplane_rotational_spring_eb as eb

REFLECTION = np.diag([1., -1., -1., 1., -1., -1.])


def quadrature(nodes=129):
    if nodes < 5 or nodes % 2 != 1:
        raise ValueError("Composite Simpson quadrature needs an odd node count >=5")
    xi = np.linspace(0., 1., nodes)
    weights = np.ones(nodes)
    weights[1:-1:2] = 4
    weights[2:-1:2] = 2
    return xi, weights / (3*(nodes-1))


def arm_states(omega, arm, reactions, xi):
    """Exact constant-coefficient EB functions, reactions in physical units.

    This is exp(H*x) C_clamp r, evaluated together at all material points.
    The four bending functions retain psi=-w', M=-D*w'', Q=M'.
    """
    if not np.isfinite(omega) or omega < 0:
        raise ValueError("omega must be finite and nonnegative")
    xi = np.asarray(xi, dtype=float)
    if not np.all(np.isfinite(xi)) or np.any((xi < 0) | (xi > 1)):
        raise ValueError("xi must lie in [0,1]")
    n, q, moment = np.asarray(reactions)
    x = arm.L*xi
    y = np.zeros((len(x), 6), dtype=np.result_type(reactions, float))
    if omega == 0:
        y[:, 0] = n*x/arm.A
        y[:, 1] = -moment*x*x/(2*arm.D)-q*x**3/(6*arm.D)
        y[:, 2] = moment*x/arm.D+q*x*x/(2*arm.D)
        y[:, 3], y[:, 4], y[:, 5] = n, q, moment+q*x
        return y
    a = omega*np.sqrt(arm.m/arm.A)
    z = (arm.m*omega**2/arm.D)**.25
    sh, ch, sn, cs = np.sinh(z*x), np.cosh(z*x), np.sin(z*x), np.cos(z*x)
    y[:, 0] = n*np.sin(a*x)/(arm.A*a)
    y[:, 3] = n*np.cos(a*x)
    y[:, 1] = -moment*(ch-cs)/(2*arm.D*z*z)-q*(sh-sn)/(2*arm.D*z**3)
    y[:, 2] = moment*(sh+sn)/(2*arm.D*z)+q*(ch-cs)/(2*arm.D*z*z)
    y[:, 5] = moment*(ch+cs)/2+q*(sh+sn)/(2*z)
    y[:, 4] = moment*z*(sh-sn)/2+q*(ch+cs)/2
    return y


def mass_vector(states, arm, weights):
    """Fixed material coordinates, local u,w; no rotations/forces in EB mass."""
    displacements = np.asarray(states)[..., :2]
    return (displacements*np.sqrt(arm.m*arm.L*np.asarray(weights))[None, :, None]).ravel()


def normalize(states, reactions, arm, weights):
    vector = mass_vector(states, arm, weights)
    mass = float(np.vdot(vector, vector).real)
    if not np.isfinite(mass) or mass <= 0:
        raise ValueError("nonpositive modal mass")
    factor = np.sqrt(mass)
    return states/factor, reactions/factor, vector/factor, mass


def reflect(states):
    return np.asarray(states)[::-1] @ REFLECTION


def mac_matrix(left, right):
    left, right = np.asarray(left), np.asarray(right)
    norms = np.sum(abs(left)**2, axis=1)[:, None]*np.sum(abs(right)**2, axis=1)[None, :]
    if np.any(norms <= 0):
        raise ValueError("zero mass norm")
    return np.abs(left.conj() @ right.T)**2/norms


def assign(left, right, left_classes=None, right_classes=None):
    """Global bijective MAC assignment, frequency absent from the cost."""
    matrix = mac_matrix(left, right)
    cost = 1-matrix
    if left_classes is not None:
        cost = cost.copy()
        cost[np.asarray(left_classes)[:, None] != np.asarray(right_classes)[None, :]] = 1e6
    rows, columns = linear_sum_assignment(cost)
    if len(rows) != len(left):
        raise ValueError("insufficient candidates")
    margins = []
    for i, j in zip(rows, columns):
        eligible = cost[i] < 1e5
        eligible[j] = False
        competing = max(matrix[i, eligible], default=0.)
        margins.append(float(matrix[i, j]-competing))
    return columns, matrix, np.asarray(margins)


def orthonormal_subspace(vectors, frequencies):
    """Only a genuinely repeated eigenvalue permits rotating its mode basis.

    The caller must establish multiplicity independently. Distinct resolved
    frequencies are rejected here rather than silently mixed.
    """
    frequencies = np.asarray(frequencies)
    if not np.all(frequencies == frequencies[0]):
        raise ValueError("cannot mix distinct resolved eigenfrequencies")
    vectors = np.asarray(vectors)
    gram = vectors.conj() @ vectors.T
    values, basis = np.linalg.eigh(gram)
    if min(values) <= 1e-12*max(values):
        raise ValueError("linearly dependent subspace")
    transform = (basis/np.sqrt(values)) @ basis.conj().T
    # Rows store modes, while G_ij=<v_i,v_j>; transpose the Hermitian
    # inverse square root so complex phase/basis changes also whiten G.
    return transform.T @ vectors


def principal_correlations(left, right):
    """Input rows already mass-orthonormal; basis-independent diagnostic."""
    return np.linalg.svd(np.asarray(left).conj() @ np.asarray(right).T, compute_uv=False)


def class_conditions(beta_rad, k_theta, parity):
    if parity not in (-1, 1):
        raise ValueError("reflection class is +1 or -1")
    c, s = np.cos(beta_rad/2), np.sin(beta_rad/2)
    if parity == 1:
        return np.array([[c,-s,0,0,0,0],[0,0,0,s,c,0],[0,0,2*k_theta,0,0,1.]])
    return np.array([[s,c,0,0,0,0],[0,0,0,c,-s,0],[0,0,0,0,0,1.]])


def class_matrix(endpoint_map, beta_rad, k_theta, arm, parity):
    scales = np.array([arm.L, arm.D/arm.L**2, arm.D/arm.L])
    if parity == 1:
        scales[2] *= max(1.,2*k_theta*arm.L/arm.D)
    return (class_conditions(beta_rad,k_theta,parity) @ endpoint_map[:6,:3])/scales[:,None]


def symmetry_equivalence(beta_rad, k_theta, arm):
    """Row-space identity of full six constraints and the two exact classes."""
    eye = np.eye(6)
    transform = np.block([[eye,eye],[REFLECTION,-REFLECTION]])/np.sqrt(2)
    separated = block_diag(class_conditions(beta_rad,k_theta,1),
                           class_conditions(beta_rad,k_theta,-1)) @ transform.T
    joint = eb.joint_matrix(beta_rad,eb.Joint("SPRING",k_theta))
    # Solve for the row transformation after fixed dimensional normalization.
    units = np.array([arm.L,arm.L,arm.D/arm.L,arm.D/arm.L**2,arm.D/arm.L**2,arm.D/arm.L])
    separated_units = np.tile([arm.L,arm.D/arm.L**2,arm.D/arm.L],2)
    a,b = separated/separated_units[:,None],joint/units[:,None]
    row_map = np.linalg.solve(b@b.T,b@a.T).T
    return dict(row_error=float(np.linalg.norm(a-row_map@b)/np.linalg.norm(a)),
                row_rank=int(np.linalg.matrix_rank(row_map)),
                reflection_commutator=float(np.linalg.norm(REFLECTION@eb.state_matrix(.5,arm)-eb.state_matrix(.5,arm)@REFLECTION)))


def recover(assembly, omega, beta_rad, joint, arm, nodes=129, parity=None):
    """Undo equilibration and reaction scales before reconstructing both arms.

    endpoint_diagnostics performs exactly that inverse transformation. A
    reflection projection selects an invariant class only at this frequency;
    it never combines solutions at different frequencies.
    """
    endpoint = eb.endpoint_diagnostics(assembly,beta_rad,joint,arm)
    if endpoint['nullity'] < 1 or endpoint['sigma_ratio'] > 1e-9:
        raise ValueError("FULL_MATRIX_ROOT_GATE")
    xi, weights = quadrature(nodes)
    possibilities = []
    for record in endpoint['vectors']:
        reactions = np.asarray(record['physical_clamp_reactions']).reshape(2,3)
        if parity is not None:
            mirrored = reactions[::-1]*np.array([1,-1,-1])
            projected = (reactions+parity*mirrored)/2
            if np.linalg.norm(projected) < 1e-6*np.linalg.norm(reactions):
                continue
            reactions = projected
        states = np.array([arm_states(omega,arm,r,xi) for r in reactions])
        states,reactions,vector,mass = normalize(states,reactions,arm,weights)
        mirrored = mass_vector(reflect(states),arm,weights)
        symmetry = float(np.vdot(vector,mirrored).real)
        eta = 1 if symmetry >= 0 else -1
        defect = float(np.linalg.norm(mirrored-eta*vector))
        state_units = np.tile([arm.L,arm.L,1,arm.D/arm.L**2,arm.D/arm.L**2,arm.D/arm.L],2)
        ends = states[:,-1,:].ravel()
        amplitude = max(abs(ends/state_units))
        residual = abs(eb.scalar_joint_residuals(ends/amplitude,beta_rad,joint)/assembly.row_units)
        hat = reactions.ravel()/assembly.reaction_scales
        boundary_residual = np.linalg.norm(assembly.dimensionless@hat)/(np.linalg.norm(assembly.dimensionless)*np.linalg.norm(hat))
        failures = []
        if max(residual)>1e-9 or max(residual[:2])>1e-10:failures.append('RECONSTRUCTED_PHYSICAL_GATE')
        if boundary_residual>1e-9:failures.append('RECONSTRUCTED_NULL_GATE')
        if defect>1e-6:failures.append('SYMMETRY_UNRESOLVED')
        possibilities.append(dict(states=states,reactions=reactions,vector=vector,mass_before_normalization=mass,
            symmetry_class=eta,symmetry_defect=defect,physical_residuals=residual.tolist(),
            null_residual=float(boundary_residual),sigma_ratio=endpoint['sigma_ratio'],
            detected_nullity=endpoint['nullity'],failures=failures))
    if not possibilities:raise ValueError('NO_CLASS_IN_NULLSPACE')
    return min(possibilities,key=lambda d:max(d['physical_residuals']))
