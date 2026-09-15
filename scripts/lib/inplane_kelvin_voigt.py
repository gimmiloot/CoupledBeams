"""Narrow complex EB/RLB adapter; exp(p*t), fixed scales, full joint matrix.

Source: inplane_kelvin_voigt_joint_theory.tex, eqs. Hp, Jmatrix, scaledB,
newton, derivativeB, modalidentity. Elastic APIs remain real-only.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import block_diag, expm, expm_frechet

from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib.inplane_rotational_spring_eb_modes import quadrature, mac_matrix
from scripts.lib.inplane_spring_modes import REFERENCE, section

T_REF = REFERENCE.L**2 * np.sqrt(REFERENCE.m / REFERENCE.D)
F_REF = REFERENCE.D / REFERENCE.L**2
M_REF = REFERENCE.D / REFERENCE.L
CRITERIA = dict(null_residual=1e-9, sigma_ratio=1e-9, physical_residual=1e-9,
    compatibility=1e-10, energy_residual=1e-6, frequency_rtol=1e-6, MAC=.95,
    a_atol=1e-8, a_rtol=1e-3, matrix_atol=1e-12, H_rtol=1e-12,
    transfer_rtol=1e-9, derivative_rtol=1e-7, quadrature_rtol=1e-6,
    symmetry_defect=1e-6, inactive_delta=1e-8, newton_residual=1e-13,
    newton_step=1e-11, simple_sigma_separation=1e-8, max_steps=20,
    max_retries_per_branch=2, max_intermediate=4, max_complex_evaluations=2000)


def real(value, name, *, positive=False):
    if np.iscomplexobj(value):
        raise ValueError(f'{name} must be real')
    value = float(value)
    if not np.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(f'invalid {name}')
    return value


def spectral(value):
    value = complex(value)
    if not np.isfinite(value):
        raise ValueError('finite complex spectral argument required')
    return value


@dataclass(frozen=True)
class Arm:
    model: str
    A: float
    D: float
    m: float
    L: float
    invS: float
    J: float

    def __post_init__(self):
        if self.model not in ('EB', 'RLB'):
            raise ValueError('model must be EB or RLB')
        for name in ('A', 'D', 'm', 'L', 'invS', 'J'):
            object.__setattr__(self, name, real(getattr(self, name), name,
                               positive=name in ('A', 'D', 'm', 'L')))
        if self.model == 'EB' and (self.invS or self.J):
            raise ValueError('classical EB requires invS=J=0')

    @classmethod
    def reduced(cls, model, properties, L=1.):
        return cls(model, properties.A, properties.D, properties.m, L,
                   1/properties.S if model == 'RLB' else 0.,
                   properties.J if model == 'RLB' else 0.)

    def scale(self):
        return np.array([self.L, self.L, 1., self.A, self.D/self.L**2, self.D/self.L])


def state_matrix(p, arm, *, derivative=False):
    p = spectral(p)
    H = np.zeros((6, 6), dtype=complex)
    if not derivative:
        H[0, 3], H[1, 2], H[1, 4] = 1/arm.A, -1, arm.invS
        H[2, 5], H[5, 4] = 1/arm.D, 1
    factor = 2*p if derivative else p*p
    H[3, 0] = H[4, 1] = arm.m*factor
    H[5, 2] = arm.J*factor
    return H


def joint_matrix(p, beta_rad, k_theta, c_theta, *, derivative=False):
    p = spectral(p)
    if np.iscomplexobj(beta_rad) or not np.isfinite(beta_rad):
        raise ValueError('beta_rad must be finite real radians')
    k, c = real(k_theta, 'k_theta'), real(c_theta, 'c_theta')
    if derivative:
        result = np.zeros((6, 12), dtype=complex)
        result[2, [2, 8]] = [c, -c]
    else:
        result = eb.joint_matrix(beta_rad, eb.Joint('SPRING', k)).astype(complex)
        result[2, 2] += c*p
        result[2, 8] -= c*p
    return result


def scalar_conditions(states, p, beta_rad, k, c):
    u1,w1,psi1,n1,q1,m1,u2,w2,psi2,n2,q2,m2 = np.asarray(states).ravel()
    cb, sb = np.cos(beta_rad), np.sin(beta_rad)
    return np.array([u1+cb*u2+sb*w2, w1-sb*u2+cb*w2,
        m1+(k+c*p)*(psi1-psi2), n1-cb*n2-sb*q2,
        q1+sb*n2-cb*q2, m1+m2], dtype=complex)


@dataclass
class Calls:
    B: int = 0
    B_z: int = 0
    expm: int = 0
    frechet: int = 0
    shape_expm: int = 0
    direct_shape_expm: int = 0
    recoveries: int = 0
    corrections: int = 0
    limit: int | None = None

    def matrix(self, derivative):
        required = 2 if derivative else 1
        if self.limit is not None and self.B+self.B_z+required > self.limit:
            raise RuntimeError('COST_LIMIT')
        self.B += 1
        self.B_z += int(derivative)

    def snapshot(self):
        return asdict(self)


class Provider:
    """Analytic B_hat(z); all row/column/state scales are independent of z."""
    def __init__(self, arms, beta_rad, kappa_theta, d_theta, calls=None):
        self.arms = tuple(arms)
        if len(self.arms) != 2 or self.arms[0].model != self.arms[1].model:
            raise ValueError('two arms of one theory required')
        self.beta = float(beta_rad)
        self.kappa = real(kappa_theta, 'kappa_theta')
        self.d = real(d_theta, 'd_theta')
        self.k, self.c = self.kappa*M_REF, self.d*M_REF*T_REF
        joint_matrix(0j, beta_rad, self.k, self.c)  # also validates radians
        self.calls = calls if calls is not None else Calls()
        self.reaction_scales = np.tile([F_REF, F_REF, M_REF], 2)
        self.row_units = np.array([REFERENCE.L, REFERENCE.L, M_REF, F_REF, F_REF, M_REF])
        self.state_units = np.tile([REFERENCE.L, REFERENCE.L, 1., F_REF, F_REF, M_REF], 2)

    def transfer(self, z, arm, *, derivative=False, x=None):
        p = spectral(z)/T_REF
        x = arm.L if x is None else real(x, 'x')
        scale = arm.scale()
        X = state_matrix(p, arm)*scale[None, :]/scale[:, None]*x
        if derivative:
            X_z = state_matrix(p, arm, derivative=True)*scale[None, :]/scale[:, None]*x/T_REF
            self.calls.frechet += 1
            T, T_z = expm_frechet(X, X_z)
            return (scale[:, None]*T/scale[None, :],
                    scale[:, None]*T_z/scale[None, :])
        self.calls.expm += 1
        return scale[:, None]*expm(X)/scale[None, :]

    def matrices(self, z, *, derivative=False):
        self.calls.matrix(derivative)
        values = [self.transfer(z, a, derivative=derivative) for a in self.arms]
        p = spectral(z)/T_REF
        joint = joint_matrix(p, self.beta, self.k, self.c)
        if derivative:
            maps = block_diag(*(v[0][:, 3:] for v in values))
            maps_z = block_diag(*(v[1][:, 3:] for v in values))
        else:
            maps = block_diag(*(v[:, 3:] for v in values))
        physical = joint@maps
        scale = self.reaction_scales[None, :]/self.row_units[:, None]
        B = physical*scale
        if not derivative:
            return B, None
        B_z = (joint_matrix(p, self.beta, self.k, self.c, derivative=True)@maps/T_REF
               +joint@maps_z)*scale
        return B, B_z


def right_null(B):
    return np.linalg.svd(B)[2].conj().T[:, -1]


def correct(matrix, z, initial_a, criteria=CRITERIA):
    """One augmented complex Newton corrector, also testable on a scalar pencil."""
    z = spectral(z)
    a = np.asarray(initial_a, dtype=complex).copy()
    norm = np.linalg.norm(a)
    if not np.isfinite(norm) or norm == 0:
        raise ValueError('finite nonzero initial reactions required')
    a /= norm
    h = a.copy()
    history, last = [], None
    for iteration in range(criteria['max_steps']+1):
        B, B_z = matrix(z, derivative=True)
        f = B@a
        # Scalar matrices have |B*a|/(|B|*|a|)=1 away from zero.
        if B.shape == (1, 1):
            residual = float(abs(f[0])/max(1., np.linalg.norm(B_z)*np.linalg.norm(a)))
        else:
            residual = float(np.linalg.norm(f)/(np.linalg.norm(B)*np.linalg.norm(a)))
        gauge = np.vdot(h, a)-1
        if residual <= criteria['newton_residual'] and abs(gauge) <= criteria['newton_residual']:
            if last is None or abs(last) <= criteria['newton_step']*max(1., abs(z)):
                return dict(z=z, a=a, steps=iteration, history=history,
                            last_delta_z=last, status='CONVERGED')
        if iteration == criteria['max_steps']:
            break
        augmented = np.block([[B, (B_z@a)[:, None]], [h.conj()[None, :], np.zeros((1, 1))]])
        change = np.linalg.solve(augmented, -np.r_[f, gauge])
        a += change[:-1]
        z += change[-1]
        last = complex(change[-1])
        history.append(dict(z_real=z.real, z_imag=z.imag, delta_real=last.real,
                            delta_imag=last.imag, residual=residual))
    return dict(z=z, a=a, steps=criteria['max_steps'], history=history,
                last_delta_z=last, status='NEWTON_LIMIT')


def mass_vector(states, arms, weights):
    return np.concatenate([(y[:, :3]*np.sqrt(a.L*weights[:, None]
        *np.array([a.m, a.m, a.J]))).ravel() for y, a in zip(states, arms)])


def recover(provider, z, reactions_hat, *, nodes=129, direct_check=False):
    """Complex full states from clamp reactions, with one total mass factor."""
    provider.calls.recoveries += 1
    xi, weights = quadrature(nodes)
    reactions = np.asarray(reactions_hat, complex)*provider.reaction_scales
    states, errors = [], []
    p = spectral(z)/T_REF
    for arm, r in zip(provider.arms, reactions.reshape(2, 3)):
        scale = arm.scale()
        H = state_matrix(p, arm)*scale[None, :]/scale[:, None]
        provider.calls.shape_expm += 1
        step = expm(H*arm.L/(nodes-1))
        initial = np.r_[np.zeros(3), r]/scale
        work = initial.copy()
        values = np.empty((nodes, 6), complex)
        for i in range(nodes):
            values[i] = scale*work
            if i+1 < nodes:
                work = step@work
        if direct_check:
            for i in (nodes//4, nodes//2, nodes-1):
                provider.calls.direct_shape_expm += 1
                exact = expm(H*arm.L*xi[i])@initial
                errors.append(float(np.linalg.norm(values[i]/scale-exact)/np.linalg.norm(exact)))
        states.append(values)
    states = np.asarray(states)
    vector = mass_vector(states, provider.arms, weights)
    mass = float(np.vdot(vector, vector).real)
    if not np.isfinite(mass) or mass <= 0:
        raise ValueError('invalid modal mass')
    factor = np.sqrt(mass)
    return dict(states=states/factor, reactions=reactions/factor,
        a=np.asarray(reactions_hat)/factor, vector=vector/factor,
        mass_before_normalization=mass, direct_errors=errors)


def diagnose(provider, z, shape):
    B, B_z = provider.matrices(z, derivative=True)
    U, singular, Vh = np.linalg.svd(B)
    a, y = shape['a'], shape['states']
    p = spectral(z)/T_REF
    ends = y[:, -1, :].ravel()
    amplitude = np.max(abs(ends/provider.state_units))
    physical = abs(scalar_conditions(ends/amplitude, p, provider.beta,
                                   provider.k, provider.c)/provider.row_units)
    _, weights = quadrature(y.shape[1])
    v = mass_vector(y, provider.arms, weights)
    M = float(np.vdot(v, v).real)
    delta = complex(y[0, -1, 2]-y[1, -1, 2])
    K = sum(float(np.dot(weights, abs(yi[:, 3])**2/arm.A
        +abs(yi[:, 5])**2/arm.D+arm.invS*abs(yi[:, 4])**2))*arm.L
        for yi, arm in zip(y, provider.arms))+provider.k*abs(delta)**2
    C = provider.c*abs(delta)**2
    r_E = abs(p*p*M+p*C+K)/(abs(p)**2*M+abs(p)*C+K)
    alpha_energy = C/(2*M)
    right = Vh.conj().T[:, -1]
    coupling = float(abs(np.vdot(U[:, -1], B_z@right)))
    conjugate, _ = provider.matrices(np.conj(z))
    mirror = y[::-1]*np.array([1, -1, -1, 1, -1, -1])
    mv = mass_vector(mirror, provider.arms, weights)
    reflection = float(np.vdot(v, mv).real/M)
    eta = 1 if reflection >= 0 else -1
    return dict(null_residual=float(np.linalg.norm(B@a)/(np.linalg.norm(B)*np.linalg.norm(a))),
        sigma_ratio=float(singular[-1]/singular[0]),
        next_sigma_ratio=float(singular[-2]/singular[0]),
        physical_residuals=physical.tolist(), clamp_residual=float(np.max(abs(y[:, 0, :3]))),
        conjugate_residual=float(np.linalg.norm(conjugate@a.conj())/(np.linalg.norm(conjugate)*np.linalg.norm(a))),
        M_phi=M, K_phi=K, C_phi=C, Delta_psi=delta, r_E=float(r_E),
        alpha_energy=float(alpha_energy), a_energy=float(alpha_energy*T_REF),
        omega_identity_residual=float(abs(p.imag**2-(K/M-p.real**2))/max(p.imag**2, K/M)),
        symmetry_class=eta, symmetry_defect=float(np.linalg.norm(mv-eta*v)/np.linalg.norm(v)),
        left_Bz_right=coupling, inverse_coupling=None if coupling == 0 else 1/coupling)


def failures(diag, z, role, Omega0, MAC):
    gates = [('null_residual','null_residual'), ('sigma_ratio','sigma_ratio'),
             ('r_E','energy_residual'), ('symmetry_defect','symmetry_defect'),
             ('conjugate_residual','null_residual'), ('clamp_residual','compatibility'),
             ('omega_identity_residual','energy_residual')]
    out = [name.upper() for name, gate in gates if diag[name] > CRITERIA[gate]]
    if max(diag['physical_residuals']) > CRITERIA['physical_residual']:
        out.append('PHYSICAL_GATE')
    if max(diag['physical_residuals'][:2]) > CRITERIA['compatibility']:
        out.append('COMPATIBILITY_GATE')
    if MAC < CRITERIA['MAC']:
        out.append('TRACKING_AMBIGUOUS')
    if diag['next_sigma_ratio'] < CRITERIA['simple_sigma_separation']:
        out.append('POSSIBLE_MULTIPLICITY')
    a = -z.real
    if z.imag <= 0 or abs(z.imag) <= 1e-6*abs(z):
        out.append('NONOSCILLATORY_OR_AXIS_APPROACH')
    if abs(a-diag['a_energy']) > CRITERIA['a_atol']+CRITERIA['a_rtol']*max(abs(a), abs(diag['a_energy'])):
        out.append('DECAY_ENERGY_GATE')
    if role == 'INACTIVE':
        if abs(a)>CRITERIA['a_atol'] or abs(diag['Delta_psi'])>CRITERIA['inactive_delta']:
            out.append('INACTIVE_GATE')
        if abs(z.imag-Omega0)>CRITERIA['frequency_rtol']*Omega0:
            out.append('INACTIVE_FREQUENCY_GATE')
    return out
