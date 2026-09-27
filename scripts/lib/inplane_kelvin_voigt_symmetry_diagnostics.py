"""D16 diagnostic only: exact EB reflection blocks of the unchanged K12 pencil.

No production solver replacement. Fixed K12 reference units are used in Newton;
row/column equilibration below is a separately labelled rank diagnostic only.
"""
from dataclasses import dataclass, asdict
import numpy as np
from scipy.linalg import expm, block_diag
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_rotational_spring_eb_modes as modes

F = modes.REFLECTION.copy()
R = F[3:, 3:].copy()
HALF_UNITS = np.array([kv.REFERENCE.L, kv.F_REF, kv.M_REF])
SELECTION = (('A_STRONG', 0., 'sorted_05'), ('C_WEAK_ACTIVE', 75., 'sorted_05'))
LOCAL_BETAS = (70., 72.5, 75., 77.5, 80.)


@dataclass
class Calls(kv.Calls):
    full_B: int = 0
    full_B_z: int = 0
    half_B: int = 0
    half_B_z: int = 0
    analytic_transfer: int = 0
    half_recoveries: int = 0
    budget: int = 1000

    def cost(self):
        # A Frechet call returns both T and T_z: charge two expm equivalents.
        return (self.B+self.B_z+self.expm+2*self.frechet+self.shape_expm
                +self.direct_shape_expm+self.analytic_transfer)

    def reserve(self, maximum):
        if self.cost()+maximum > self.budget:
            raise RuntimeError('COST_LIMIT')

    def snapshot(self):
        return dict(asdict(self), total_build_equivalents=self.cost())


class FullProvider(kv.Provider):
    def __init__(self, arms, beta_rad, kappa_theta, d_theta, calls):
        if arms[0] != arms[1] or any(a.model != 'EB' for a in arms):
            raise ValueError('diagnostic requires identical EB arms')
        super().__init__(arms, beta_rad, kappa_theta, d_theta, calls)

    def matrices(self, z, *, derivative=False):
        self.calls.reserve(6 if derivative else 3)
        self.calls.full_B += 1
        self.calls.full_B_z += int(derivative)
        return super().matrices(z, derivative=derivative)


def conditions(p, beta, k, c, eta, *, derivative=False):
    """Reuse the verified K09 half conditions; complexify only K_J=k+c*p."""
    if eta not in (-1, 1):
        raise ValueError('eta must be +1 or -1')
    if derivative:
        value = np.zeros((3, 6), complex)
        if eta == 1:
            value[2, 2] = 2*c
        return value
    return modes.class_conditions(beta, k+c*kv.spectral(p), eta).astype(complex)


def lift(a, eta):
    if eta not in (-1, 1):
        raise ValueError('eta must be +1 or -1')
    a = np.asarray(a, complex)
    return np.r_[a, eta*R@a]


def project(a):
    a = np.asarray(a, complex).reshape(2, 3)
    plus = (a[0]+R@a[1])/2
    minus = (a[0]-R@a[1])/2
    return lift(plus, 1), lift(minus, -1)


def row_transform(beta):
    """On K12 dimensionless full residuals: plus rows then minus rows."""
    c, s = np.cos(beta/2), np.sin(beta/2)
    return np.array([[c/2,-s/2,0,0,0,0], [0,0,0,s/2,c/2,0],
        [0,0,1,0,0,-.5], [s/2,c/2,0,0,0,0],
        [0,0,0,c/2,-s/2,0], [0,0,0,0,0,.5]])


class HalfProvider:
    def __init__(self, full, eta):
        if eta not in (-1, 1):
            raise ValueError('eta must be +1 or -1')
        self.full, self.eta = full, eta

    def matrices(self, z, *, derivative=False):
        f = self.full
        f.calls.reserve(4 if derivative else 2)
        f.calls.matrix(derivative)
        f.calls.half_B += 1
        f.calls.half_B_z += int(derivative)
        value = f.transfer(z, f.arms[0], derivative=derivative)
        C = conditions(z/kv.T_REF, f.beta, f.k, f.c, self.eta)
        scales = f.reaction_scales[None, :3]/HALF_UNITS[:, None]
        if not derivative:
            return (C@value[:, 3:])*scales, None
        T, Tz = value
        Cp = conditions(z/kv.T_REF, f.beta, f.k, f.c, self.eta, derivative=True)
        return (C@T[:, 3:])*scales, (Cp@T[:, 3:]/kv.T_REF+C@Tz[:, 3:])*scales


def spectrum(B, Bz=None, a=None):
    U, s, Vh = np.linalg.svd(B)
    if a is None:
        a = Vh.conj().T[:, -1]
    residual = B@a
    rows = np.linalg.norm(B, axis=1)
    row_scaled = B/np.where(rows > 0, rows, 1)[:, None]
    cols = np.linalg.norm(row_scaled, axis=0)
    balanced = row_scaled/np.where(cols > 0, cols, 1)[None, :]
    bs = np.linalg.svd(balanced, compute_uv=False)
    result = dict(singular_values=s.tolist(), ratios=(s/s[0]).tolist(),
        balanced_ratios=(bs/bs[0]).tolist(), row_norms=rows.tolist(),
        column_norms_after_row_scaling=cols.tolist(), matrix_norm=float(np.linalg.norm(B)),
        null_raw_norm=float(np.linalg.norm(residual)),
        null_residual=float(np.linalg.norm(residual)/(np.linalg.norm(B)*np.linalg.norm(a))))
    if Bz is not None:
        coupling = float(abs(np.vdot(U[:, -1], Bz@Vh.conj().T[:, -1])))
        result.update(left_Bz_right=coupling, inverse_coupling=1/coupling if coupling else None)
    return result


def algebra_check(full, z):
    B, Bz = full.matrices(z, derivative=True)
    halves = [HalfProvider(full, eta).matrices(z, derivative=True) for eta in (1, -1)]
    V = np.block([[np.eye(3), np.eye(3)], [R, -R]])
    W = row_transform(full.beta)
    expected = block_diag(*(v[0] for v in halves))
    expected_z = block_diag(*(v[1] for v in halves))
    got = W@B@V
    off = got.copy(); off[:3, :3] = 0; off[3:, 3:] = 0
    H = kv.state_matrix(z/kv.T_REF, full.arms[0])
    result = dict(z=z, block_relative=float(np.linalg.norm(got-expected)/np.linalg.norm(expected)),
        block_absolute=float(np.linalg.norm(got-expected)),
        derivative_relative=float(np.linalg.norm(W@Bz@V-expected_z)/np.linalg.norm(expected_z)),
        off_block_relative=float(np.linalg.norm(off)/np.linalg.norm(got)),
        FH_HF_norm=float(np.linalg.norm(F@H-H@F)), row_transform_rank=int(np.linalg.matrix_rank(W)))
    result['accepted'] = bool(max(result['block_relative'], result['derivative_relative'],
        result['off_block_relative']) < kv.CRITERIA['H_rtol'] and result['FH_HF_norm'] == 0
        and result['row_transform_rank'] == 6)
    return result


def normalize_mirrored(full, first, a, eta):
    states = np.array([first, eta*first@F])
    _, weights = kv.quadrature(len(first))
    v = kv.mass_vector(states, full.arms, weights)
    mass = float(np.vdot(v, v).real)
    hat = lift(a, eta)
    return dict(states=states/np.sqrt(mass), a=hat/np.sqrt(mass),
        reactions=hat*full.reaction_scales/np.sqrt(mass), vector=v/np.sqrt(mass),
        mass_before_normalization=mass)


def recover_half(half, z, a, nodes=129):
    full, eta = half.full, half.eta
    full.calls.reserve(1)
    full.calls.recoveries += 1
    full.calls.half_recoveries += 1
    full.calls.shape_expm += 1
    arm = full.arms[0]; scale = arm.scale()
    H = kv.state_matrix(z/kv.T_REF, arm)*scale[None, :]/scale[:, None]
    step = expm(H*arm.L/(nodes-1))
    work = np.r_[np.zeros(3), np.asarray(a)*full.reaction_scales[:3]]/scale
    first = np.empty((nodes, 6), complex)
    for i in range(nodes):
        first[i] = scale*work
        if i+1 < nodes:
            work = step@work
    return normalize_mirrored(full, first, a, eta)


def physical_details(half, z, shape):
    f = half.full; ends = shape['states'][:, -1, :]
    amplitude = float(np.max(abs(ends.ravel()/f.state_units)))
    full_raw = kv.scalar_conditions(ends.ravel(), z/kv.T_REF, f.beta, f.k, f.c)
    half_raw = conditions(z/kv.T_REF, f.beta, f.k, f.c, half.eta)@ends[0]
    return dict(full_raw=full_raw, half_raw=half_raw, endpoint_amplitude=amplitude,
        full_normalized=(abs(full_raw/f.row_units)/amplitude).tolist(),
        half_normalized=(abs(half_raw/HALF_UNITS)/amplitude).tolist())


def beta_trigger(*, opposite_block_near, same_block_suspect, second_sigma_unexplained):
    return bool(opposite_block_near or same_block_suspect or second_sigma_unexplained)


def local_beta_plan(trigger):
    # This helper cannot expand the allowed grid or localize a crossing.
    return LOCAL_BETAS if trigger else ()


def closed_transfer(p, arm, x=None, *, derivative=False):
    """Conditional D16 control only, independently evaluated EB matrix functions.

    Axial H_a^2=q^2 I, q=p*sqrt(m/A); bending H_b^4=lambda^4 I,
    lambda^4=-m*p^2/D. Principal fourth root is used; lambda->i*lambda
    leaves the polynomial in H_b unchanged. Derivative is with respect to p.
    No production transfer or matrix exponential is called here.
    """
    if arm.model != 'EB' or arm.invS or arm.J:
        raise ValueError('closed control is EB only')
    p = kv.spectral(p); x = arm.L if x is None else float(x)
    H = kv.state_matrix(p, arm); Hp = kv.state_matrix(p, arm, derivative=True)
    if p == 0:
        T = np.eye(6, dtype=complex)
        power = np.eye(6, dtype=complex)
        for j in range(1, 4):
            power = power@H*x/j
            T += power
        return (T, np.zeros_like(T)) if derivative else T
    q = p*np.sqrt(arm.m/arm.A); qp = np.sqrt(arm.m/arm.A)
    axial = [0, 3]; bend = [1, 2, 4, 5]
    Ha = H[np.ix_(axial, axial)]; Hap = Hp[np.ix_(axial, axial)]
    Ta = np.cosh(q*x)*np.eye(2)+np.sinh(q*x)/q*Ha
    Tap = (qp*x*np.sinh(q*x)*np.eye(2)
           +qp*(x*np.cosh(q*x)/q-np.sinh(q*x)/q**2)*Ha+np.sinh(q*x)/q*Hap)
    lam = complex(-arm.m*p*p/arm.D)**.25
    lp = lam/(2*p); t = lam*x
    f = np.array([(np.cosh(t)+np.cos(t))/2, (np.sinh(t)+np.sin(t))/2,
                  (np.cosh(t)-np.cos(t))/2, (np.sinh(t)-np.sin(t))/2])
    # f'_0=f_3, f'_1=f_0, f'_2=f_1, f'_3=f_2.
    Hb = H[np.ix_(bend, bend)]; Hbp = Hp[np.ix_(bend, bend)]
    power = np.eye(4, dtype=complex); power_p = np.zeros((4, 4), complex)
    Tb = np.zeros((4, 4), complex); Tbp = np.zeros_like(Tb)
    for j in range(4):
        coef = f[j]/lam**j
        coef_p = lp/lam**j*(x*f[(j-1) % 4]-j*f[j]/lam)
        Tb += coef*power
        Tbp += coef_p*power+coef*power_p
        power_p = power_p@Hb+power@Hbp
        power = power@Hb
    T = np.zeros((6, 6), complex); Tp = np.zeros_like(T)
    T[np.ix_(axial, axial)], T[np.ix_(bend, bend)] = Ta, Tb
    Tp[np.ix_(axial, axial)], Tp[np.ix_(bend, bend)] = Tap, Tbp
    return (T, Tp) if derivative else T


class AnalyticHalfProvider(HalfProvider):
    """Shared verified K17 assembly, also used by the D17 production dispatcher.

    Historical ClosedHalfProvider below retains its conditional audit contract.
    """
    def matrices(self, z, *, derivative=False):
        f = self.full
        f.calls.reserve(3 if derivative else 2)
        f.calls.matrix(derivative)
        f.calls.half_B += 1; f.calls.half_B_z += int(derivative)
        f.calls.analytic_transfer += 1
        value = closed_transfer(z/kv.T_REF, f.arms[0], derivative=derivative)
        C = conditions(z/kv.T_REF, f.beta, f.k, f.c, self.eta)
        scales = f.reaction_scales[None, :3]/HALF_UNITS[:, None]
        if not derivative:
            return (C@value[:, 3:])*scales, None
        T, Tp = value
        Cp = conditions(z/kv.T_REF, f.beta, f.k, f.c, self.eta, derivative=True)
        return (C@T[:, 3:])*scales, ((Cp@T[:, 3:]+C@Tp[:, 3:])/kv.T_REF)*scales


class ClosedHalfProvider(AnalyticHalfProvider):
    """Historical diagnostic wrapper: requires the K17 transfer trigger."""
    def __init__(self, full, eta, *, triggered):
        if not triggered:
            raise ValueError('conditional analytic diagnostic not triggered')
        super().__init__(full, eta)


class FrozenBalanced:
    """Positive equilibration at the input predictor, fixed for the whole solve.

    Uses already dimensionless rows and columns. No z-dependent scaling enters
    Newton or its derivative; physical residual definitions remain K12's.
    """
    def __init__(self, half, z):
        self.half = half
        B, _ = half.matrices(z)
        rows = np.linalg.norm(B, axis=1)
        self.rows = np.where(rows > 0, rows, 1.)
        cols = np.linalg.norm(B/self.rows[:, None], axis=0)
        self.cols = np.where(cols > 0, cols, 1.)

    def matrices(self, z, *, derivative=False):
        B, Bz = self.half.matrices(z, derivative=derivative)
        scale = self.rows[:, None]*self.cols[None, :]
        return B/scale, None if Bz is None else Bz/scale

    def reactions(self, b):
        return np.asarray(b)/self.cols


def recover_closed(half, z, a, nodes=129):
    f = half.full; arm = f.arms[0]
    f.calls.reserve(nodes)
    f.calls.recoveries += 1; f.calls.half_recoveries += 1
    initial = np.r_[np.zeros(3), np.asarray(a)*f.reaction_scales[:3]]
    xi, _ = kv.quadrature(nodes)
    first = []
    for x in xi*arm.L:
        f.calls.analytic_transfer += 1
        first.append(closed_transfer(z/kv.T_REF, arm, x)@initial)
    return normalize_mirrored(f, np.asarray(first), a, half.eta)
