"""One global horizontal spring, no outer guide or rotation spring.

Canonical note: docs/theory/isotropic_two_springs_horizontal_spring.md.
Reuse RL_FULL/Tim arm equations and the canonical node; zero-inclusive
bounded verification is deliberately separate from the old positive maps.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import block_diag, cholesky, eigh, expm, null_space, qr, solve_triangular
from scipy.optimize import brentq
from scripts.lib import isotropic_two_springs_rl_timoshenko as base

Model, RIGID = base.Model, base.RIGID
LOW_FREQUENCY_BASIS_LAMBDA = .25


@dataclass(frozen=True)
class Case:
    beta: float
    mu: float = .01
    k_joint: float | str = .1
    kappa_horizontal: float | str = 1.

    def __post_init__(self):
        base.Case(self.beta, self.mu, self.k_joint, self.kappa_horizontal)

    def lengths(self, m):
        return m.l*(1-self.mu), m.l*(1+self.mu)


def cs(case):
    # Exact cardinal angles are mathematical limits, not small-angle cutoffs.
    b = math.radians(case.beta)
    c, s = math.cos(b), math.sin(b)
    if case.beta % 180 == 0:
        s = 0.
    if (case.beta-90) % 180 == 0:
        c = 0.
    return c, s


def expected_nullity(case):
    return int(case.k_joint == 0 and (case.kappa_horizontal == 0 or cs(case)[1] == 0))


def state_matrix(m, omega):
    if omega:
        return base.state_matrix(m, omega)
    h = np.zeros((6, 6))
    h[0, 3] = 1/m.section.EA
    h[1, 2], h[1, 4] = 1., 1/m.section.KGA
    h[2, 5], h[5, 4] = 1/m.section.EI, -1.
    return h


def basis(m, length, omega, x, derivative=0):
    if omega and m.Lambda(omega) >= LOW_FREQUENCY_BASIS_LAMBDA:
        return base.basis(m, length, omega, x, derivative)
    a = state_matrix(m, omega)
    scale = np.diag([m.l, m.l, 1., m.section.EI/m.l**2,
                     m.section.EI/m.l**2, m.section.EI/m.l])
    if omega:
        # Same canonical ODE, evaluated in a static-regular fundamental basis.
        # At small frequency the anchored-wave basis loses digits through
        # cancellation; positive state scaling keeps expm well conditioned.
        ah = np.linalg.solve(scale,a@scale)
        ad = np.linalg.matrix_power(ah,derivative)
        out = np.array([scale@ad@expm(ah*t) for t in np.atleast_1d(x)])
        return out[0] if np.ndim(x)==0 else out
    powers = [np.linalg.matrix_power(a, j) for j in range(4)]
    out = np.array([sum((t**j/math.factorial(j))*powers[j]
                       for j in range(4))@np.linalg.matrix_power(a, derivative)@scale
                    for t in np.atleast_1d(x)])
    return out[0] if np.ndim(x) == 0 else out


def external_matrix(m, case, arm):
    if arm == 0:
        return np.eye(6)[:3]
    c, s = cs(case)
    e = np.zeros((3, 6))
    if case.kappa_horizontal == RIGID:
        e[0, :2], e[1, 3:5] = [c, s], [-s, c]
    else:
        k = case.kappa_horizontal*m.section.EI/m.l**3
        e[:2, :2] = k*np.outer([c, s], [c, s])
        e[0, 3], e[1, 4] = -1., -1.
    e[2, 5] = 1.
    return e


def boundary(m, case, omega):
    ends = [basis(m, L, omega, [0., L]) for L in case.lengths(m)]
    a = np.zeros((12, 12))
    for i, v in enumerate(ends):
        a[3*i:3*i+3, 6*i:6*i+6] = external_matrix(m, case, i)@v[0]
    a[6:] = base.joint_matrix(m, case)@block_diag(*(v[1] for v in ends))
    return a


def trace_map(m, case):
    nr = 1 if case.k_joint == RIGID else 2
    nt = 1 if case.kappa_horizontal == RIGID else 2
    p = np.zeros((12, 2+nr+nt+1))
    for i, f in enumerate(base.geometry.frames(case.beta)):
        p[6*i+3, :2], p[6*i+4, :2] = f.t, f.n
        p[6*i+5, 2+(i if nr == 2 else 0)] = 1.
    c, s = cs(case)
    if nt == 1:
        p[6:8, 2+nr] = [-s, c]
    else:
        p[6:8, 2+nr:2+nr+2] = np.eye(2)
    p[8, -1] = 1.
    jump, horizontal = p[5]-p[11], c*p[6]+s*p[7]
    spring = np.zeros((p.shape[1],)*2)
    if case.k_joint != RIGID:
        spring += case.k_joint*m.section.EI/m.l*np.outer(jump, jump)
    if case.kappa_horizontal != RIGID:
        spring += case.kappa_horizontal*m.section.EI/m.l**3*np.outer(horizontal, horizontal)
    return p, spring


def arm_dtn(m, L, omega):
    b = basis(m, L, omega, [0., L])
    q = np.vstack((b[0, :3], b[1, :3]))
    f = np.vstack((-b[0, 3:], b[1, 3:]))
    return np.linalg.solve(q.T, f.T).T, float(np.linalg.cond(base.balance(q)[0]))


def assembled(m, case, omega):
    pairs = [arm_dtn(m, L, omega) for L in case.lengths(m)]
    p, spring = trace_map(m, case)
    return p.T@block_diag(*(a[0] for a in pairs))@p+spring, [a[1] for a in pairs]


def congruence(m, case):
    nr = 1 if case.k_joint == RIGID else 2
    nt = 1 if case.kappa_horizontal == RIGID else 2
    return math.sqrt(m.l/m.section.EI)*np.array([m.l]*2+[1.]*nr+[m.l]*nt+[1.])


def static_kernel(m, case, tol=1e-11):
    raw = boundary(m, case, 0.)
    sv = np.linalg.svd(base.balance(raw)[0], compute_uv=False)
    a, _ = assembled(m, case, 0.)
    # Inertia is unchanged by positive congruence. Equilibrate static diagonal
    # energies: otherwise EA/EI can make a small POSITIVE bending stiffness
    # look like a null eigenvalue under a global relative rank threshold.
    if np.any(np.diag(a) <= 0):
        raise ArithmeticError('Nonpositive diagonal in static energy matrix')
    d = 1/np.sqrt(np.diag(a))
    vals = np.linalg.eigvalsh((a+a.T)/2*d[:, None]*d)
    n_boundary = int(sum(sv < tol*sv[0]))
    n_energy = int(sum(abs(vals) < tol*max(abs(vals))))
    expected = expected_nullity(case)
    ok = n_boundary == n_energy == expected and min(vals) >= -tol*max(abs(vals))
    return {"status": "PASS" if ok else "FAIL", "expected": expected,
            "boundary_nullity": n_boundary, "energy_nullity": n_energy,
            "boundary_singular_values": sv.tolist(), "energy_eigenvalues": vals.tolist(),
            "energy_congruence":d.tolist(),
            "static_skew": float(np.linalg.norm(a-a.T)/np.linalg.norm(a))}


def poles(m, case, catalogs):
    out = []
    for L in case.lengths(m):
        cat = catalogs[L]
        out.extend(cat['poles'])
        p = m.section.EA-m.J*cat['upper']**2
        if p <= 0:
            raise ArithmeticError('RL acoustic interval exceeded')
        nmax = math.floor(cat['upper']*math.sqrt(m.section.rhoA/p)*L/math.pi)
        out.extend(math.sqrt(m.section.EA*(n*math.pi/L)**2 /
                    (m.section.rhoA+m.J*(n*math.pi/L)**2)) for n in range(1, nmax+1))
    return sorted(out)


def count_at(m, case, omega, catalogs, policy):
    if omega <= 0:
        raise ValueError('Use static_kernel at omega=0; count is right-sided for omega>0')
    pp = poles(m, case, catalogs)
    w = float(omega)
    exclusion = policy['pole_exclusion']*max(1., w)
    for pole in pp:
        if abs(w-pole) <= exclusion:
            w = pole+2*exclusion
    if w >= min(cat['upper'] for cat in catalogs.values()):
        raise ArithmeticError('Outside certified pole catalog')
    a, conditions = assembled(m, case, w)
    d = congruence(m, case)
    a = a*d[:, None]*d
    skew = np.linalg.norm(a-a.T)/np.linalg.norm(a)
    vals = np.linalg.eigvalsh((a+a.T)/2)
    margin = min(abs(vals))/max(abs(vals))
    if skew > policy['schur_tol'] or margin < policy['inertia_margin'] or max(conditions) > policy['condition_max']:
        raise ArithmeticError(f'Unresolved inertia: omega={w}, skew={skew}, margin={margin}')
    j0 = sum(p < w for p in pp)
    n = j0+int(sum(vals < 0))
    return w, n, {'requested': omega, 'effective': w, 'J0': j0, 'count': n,
                   'eigenvalues': vals.tolist(), 'skew': float(skew),
                   'margin': float(margin), 'conditions': conditions}


def prefix(m, case, catalogs, policy, count=3):
    gate = static_kernel(m, case, policy['kernel_rank_relative_tol'])
    if gate['status'] != 'PASS':
        raise ArithmeticError(f'Zero gate failed: {gate}')
    n0 = gate['expected']
    samples, roots = [], []
    def query(w):
        if len(samples) >= policy['max_count_queries']:
            raise ArithmeticError('Count query budget exceeded')
        w, n, diag = count_at(m, case, w, catalogs, policy)
        samples.append(diag)
        return w, n
    def det(w):
        return float(np.linalg.det(base.balance(boundary(m, case, w))[0]))
    def isolate(a, b, ca, cb, depth=0):
        if ca >= count or ca == cb:
            return
        if cb < ca or depth > policy['max_depth']:
            raise ArithmeticError('Count isolation failed')
        if cb-ca == 1 and det(a)*det(b) < 0:
            w, info = brentq(det, a, b, xtol=policy['root_xtol'], rtol=policy['root_rtol'], full_output=True)
            roots.append({'omega': w, 'bracket': [a,b], 'counts': [ca,cb],
                          'physical_position': ca+1, 'positive_position': ca-n0+1,
                          'iterations': info.iterations})
            return
        mid, cm = query((a+b)/2)
        if not a < mid < b:
            raise ArithmeticError('Pole shift outside bracket')
        isolate(a,mid,ca,cm,depth+1)
        isolate(mid,b,cm,cb,depth+1)
    # The lower endpoint is established by inertia, not a positive Lambda floor.
    a, ca = query(float(m.omega(.25)))
    while ca > n0:
        a, ca = query(a/4)
    if ca != n0:
        raise ArithmeticError('Near-zero count contradicts the static kernel')
    b, cb = query(float(m.omega(policy['Lambda_upper'])))
    if cb < count:
        raise ArithmeticError('Upper bound misses required physical prefix')
    isolate(a,b,ca,cb)
    if len(roots) != count-n0:
        raise ArithmeticError('Incomplete physical prefix')
    records = ([{'omega':0., 'physical_position':1, 'positive_position':None,
                 'source':'exact_static_kernel'}] if n0 else [])+roots
    return records, {'status':'PASS', 'zero_gate':gate, 'samples':samples,
                      'physical_count':count, 'positive_count':len(roots)}


def mode(m, case, omega, order=120):
    raw = boundary(m, case, omega)
    balanced, _, col = base.balance(raw)
    _, singular, vh = np.linalg.svd(balanced)
    coef = (col*vh[-1]).reshape(2,6)
    if omega == 0:
        coef[:] = 0.
        coef[1,1],coef[1,2] = -case.lengths(m)[1]/m.l,1.
    nodes, weights = np.polynomial.legendre.leggauss(order)
    mass = energy = 0.
    arms = []
    sec = m.section
    for i,(L, co) in enumerate(zip(case.lengths(m),coef)):
        x, wt = (nodes+1)*L/2, weights*L/2
        y, dy = (basis(m,L,omega,x,d)@co for d in (0,1))
        ends = basis(m,L,omega,[0.,L])@co
        if omega == 0:
            if not expected_nullity(case):
                raise ValueError('No zero mode in this case')
            y[:], dy[:], ends[:] = 0., 0., 0.
            if i == 1:
                y[:,1], y[:,2], dy[:,1] = x-L, 1., 1.
                ends[:,1], ends[:,2] = [-L,0.], 1.
        mass += wt@(sec.rhoA*(y[:,0]**2+y[:,1]**2)+sec.rhoI*y[:,2]**2+m.J*dy[:,0]**2)
        energy += wt@(sec.EA*dy[:,0]**2+sec.EI*dy[:,2]**2+sec.KGA*(dy[:,1]-y[:,2])**2)
        arms.append({'x':x,'weights':wt,'y':y,'dy':dy,'ends':ends})
    delta = arms[0]['ends'][1,2]-arms[1]['ends'][1,2]
    c,s = cs(case)
    horizontal = c*arms[1]['ends'][0,0]+s*arms[1]['ends'][0,1]
    if case.k_joint != RIGID:
        energy += case.k_joint*sec.EI/m.l*delta**2
    if case.kappa_horizontal != RIGID:
        energy += case.kappa_horizontal*sec.EI/m.l**3*horizontal**2
    for arm in arms:
        for key in ('y','dy','ends'):
            arm[key] /= math.sqrt(mass)
    coef /= math.sqrt(mass)
    qscale = max(np.max(abs(a['y'][:,:3]*[1,1,m.l])) for a in arms)
    fscale = max(np.max(abs(a['y'][:,3:]*[1,1,1/m.l])) for a in arms)
    if omega == 0:
        fscale = sec.EI*qscale/m.l**3
    jr = base.joint_matrix(m,case)@np.concatenate([a['ends'][1] for a in arms])
    js = abs(jr)/np.array([qscale,qscale,fscale,fscale,
                  qscale/m.l if case.k_joint == RIGID else fscale*m.l,fscale*m.l])
    er = [external_matrix(m,case,i)@a['ends'][0] for i,a in enumerate(arms)]
    es = [abs(er[0])/[qscale,qscale,qscale/m.l],
          abs(er[1])/[qscale if case.kappa_horizontal==RIGID else fscale,fscale,fscale*m.l]]
    scales = np.array([1,1,m.l,qscale/fscale,qscale/fscale,qscale/(fscale*m.l)])
    ode = max(np.max(abs((a['dy']-a['y']@state_matrix(m,omega).T)*scales)) /
              max(np.max(abs(a['dy']*scales)),qscale/m.l) for a in arms)
    diag = {'mass':1.,'raw_mass':float(mass),'energy':float(energy/mass),
            'energy_error':float(abs(energy/mass/omega**2-1)) if omega else float(energy/mass),
            'singular_ratio':float(singular[-1]/singular[0]),
            'nonzero_condition':float(singular[0]/singular[-2]),
            'joint_scaled':js.tolist(),'external_scaled':np.array(es).tolist(),
            'joint_absolute':jr.tolist(),'external_absolute':np.array(er).tolist(),
            'ode_error':float(ode),'delta_psi':float(delta/math.sqrt(mass)),
            'zero_absolute_checks':omega==0}
    return {'coefficients':coef.tolist(),'arms':arms,'diagnostics':diag}


inner = base.inner


def ritz_assembly(m, case, order):
    """Independent full energies; constant u2,w2,psi2 are all admissible.

    Keep the weighted strain operator B (K=B.T B), so the exact zero
    mechanism can be tested without assuming that K is positive definite.
    """
    size = 6*order+3
    nodes, weights = np.polynomial.legendre.leggauss(max(64,order+2))
    bv = base.shifted_legendre_clamped_basis(order,(nodes+1)/2)
    mass = np.zeros((size,size))
    strains, traces, fields = [], [], []
    for i,L in enumerate(case.lengths(m)):
        v,d = np.zeros((len(nodes),3,size)),np.zeros((len(nodes),3,size))
        end,outer = np.zeros((3,size)),np.zeros((3,size))
        for f in range(3):
            sl = slice((3*i+f)*order,(3*i+f+1)*order)
            v[:,f,sl],d[:,f,sl],end[f,sl] = bv.values.T,bv.derivatives.T/L,1.
            if i == 1:
                v[:,f,6*order+f] = end[f,6*order+f] = outer[f,6*order+f] = 1.
        wt = weights*L/2
        sec = m.section
        for factor, op in ((sec.EA,d[:,0]),(sec.EI,d[:,2]),(sec.KGA,d[:,1]-v[:,2])):
            strains.append(np.sqrt(factor*wt)[:,None]*op)
        for factor,op in ((sec.rhoA,v[:,0]),(sec.rhoA,v[:,1]),(sec.rhoI,v[:,2]),(m.J,d[:,0])):
            mass += op.T@((factor*wt)[:,None]*op)
        traces.append(end)
        fields.append((v,d,wt))
    frames = [np.column_stack((f.t,f.n)) for f in base.geometry.frames(case.beta)]
    constraint = frames[0]@traces[0][:2]-frames[1]@traces[1][:2]
    jump = traces[0][2]-traces[1][2]
    if case.k_joint == RIGID:
        constraint = np.vstack((constraint,jump))
    elif case.k_joint:
        strains.append(math.sqrt(case.k_joint*sec.EI/m.l)*jump[None,:])
    c,s = cs(case)
    h = c*outer[0]+s*outer[1]
    if case.kappa_horizontal == RIGID:
        constraint = np.vstack((constraint,h))
    elif case.kappa_horizontal:
        strains.append(math.sqrt(case.kappa_horizontal*sec.EI/m.l**3)*h[None,:])
    z = null_space(constraint)
    return np.vstack(strains),mass,constraint,z,fields


def ritz(m, case, order, count=3, rank_tol=1e-11, kernel_only=False):
    b,mass,constraint,z,fields = ritz_assembly(m,case,order)
    br,mr = b@z,z.T@mass@z
    cholesky((mr+mr.T)/2,lower=True)
    scale = 1/np.linalg.norm(br,axis=0)
    _,sv,vh = np.linalg.svd(br*scale,full_matrices=False)
    n0 = int(sum(sv < rank_tol*sv[0]))
    if n0 != expected_nullity(case):
        raise ArithmeticError(f'Ritz kernel mismatch: {n0} != {expected_nullity(case)}')
    zero = np.empty((len(mr),0))
    if n0:
        # Exact affine mode, represented in the independently assembled space.
        exact = np.zeros(len(mass))
        L2 = case.lengths(m)[1]
        exact[4*order],exact[6*order+1],exact[6*order+2] = L2,-L2,1.
        zero = (z.T@exact)[:,None]
        zero /= math.sqrt(float(zero.T@mr@zero))
    residual = float(np.linalg.norm(br@zero)/np.linalg.norm(br)) if n0 else 0.
    diag = {'nullity':n0,'strain_singular_values':sv.tolist(),
            'kernel_residual':residual,'mass_condition':float(np.linalg.cond(mr)),
            'mass_positive':True,'admissible_dimension':len(mr)}
    if kernel_only:
        return diag
    # Deflate a VERIFIED nullspace in the M metric, with no penalty or epsilon.
    q = null_space(zero.T@mr) if n0 else np.eye(len(mr))
    bp,mp = br@q,q.T@mr@q
    _,r = qr(bp,mode='economic')
    # Inverse stiffness pencil in its positive subspace, using the strain QR.
    a = solve_triangular(r.T,mp,lower=True)
    a = solve_triangular(r.T,a.T,lower=True).T
    npos = count-n0
    val,v = eigh((a+a.T)/2,subset_by_index=[len(a)-npos,len(a)-1])
    val,v = val[::-1],v[:,::-1]
    cp = q@solve_triangular(r,v,lower=False)/np.sqrt(val)
    co = z@np.column_stack((zero,cp))
    omega = np.r_[np.zeros(n0),np.sqrt(1/val)]
    gram = co.T@mass@co
    diag.update(mass_gram_error=float(np.max(abs(gram-np.eye(count)))),
                energy_error=float(np.max(abs(np.sum((b@co[:,n0:])**2,axis=0)/omega[n0:]**2-1))),
                constraint_residual=float(np.max(abs(constraint@co))),
                matrix_residual=float(np.linalg.norm(z.T@(b.T@(b@co)-(mass@co)*omega**2)) /
                    ((np.linalg.norm(br)**2+max(omega**2)*np.linalg.norm(mr))*np.linalg.norm(co))),
                stiffness_condition=float(np.linalg.cond(r)**2))
    return omega,co,diag


def ritz_fields(m,case,order,co,quadrature=120):
    nodes,weights = np.polynomial.legendre.leggauss(quadrature)
    bv = base.shifted_legendre_clamped_basis(order,(nodes+1)/2)
    arms=[]
    for i,L in enumerate(case.lengths(m)):
        y,dy = np.zeros((quadrature,6)),np.zeros((quadrature,6))
        for f in range(3):
            cf = co[(3*i+f)*order:(3*i+f+1)*order]
            y[:,f],dy[:,f] = bv.values.T@cf,bv.derivatives.T@cf/L
            if i==1:
                y[:,f] += co[6*order+f]
        arms.append({'x':(nodes+1)*L/2,'weights':weights*L/2,'y':y,'dy':dy})
    return {'arms':arms}


def angle_pencil_coefficients(m,case,omega):
    """Exact beta-radian Taylor coefficients A0,A1,A2 at beta=0.

    Includes BOTH node rotation and the explicit horizontal spring rotation.
    A2 denotes the beta**2 coefficient (half the second derivative).
    """
    if case.beta != 0:
        raise ValueError('Local straight-configuration calculation only')
    p0,spring = trace_map(m,case)
    p1,p2 = np.zeros_like(p0),np.zeros_like(p0)
    p1[9,1],p1[10,0] = -1.,-1.
    p2[9,0],p2[10,1] = .5,-.5
    nr = 1 if case.k_joint==RIGID else 2
    if case.kappa_horizontal==RIGID:
        p1[6,2+nr],p2[7,2+nr] = -1.,-.5
    d = block_diag(*(arm_dtn(m,L,omega)[0] for L in case.lengths(m)))
    a0 = p0.T@d@p0+spring
    a1 = p1.T@d@p0+p0.T@d@p1
    a2 = p2.T@d@p0+p0.T@d@p2+p1.T@d@p1
    if case.kappa_horizontal != RIGID:
        k = case.kappa_horizontal*m.section.EI/m.l**3
        h0,h1,h2 = p0[6],p0[7],-.5*p0[6]
        a1 += k*(np.outer(h0,h1)+np.outer(h1,h0))
        a2 += k*(np.outer(h0,h2)+np.outer(h2,h0)+np.outer(h1,h1))
    return a0,a1,a2


def curvature_prediction(m,case,omega,modal):
    """Full-model perturbation, not a frequency fit or leading xi model.

    Requires the simple bending mode at beta=0 and no axial resonance.
    Input modal fields have the FULL RL/Tim mass norm equal to one.
    """
    if case.beta != 0 or omega <= 0:
        raise ValueError('Positive straight bending eigenpair required')
    p,_ = trace_map(m,case)
    endpoint = np.concatenate([a['ends'][:,:3].ravel() for a in modal['arms']])
    q0 = np.linalg.lstsq(p,endpoint,rcond=None)[0]
    a0,a1,a2 = angle_pencil_coefficients(m,case,omega)
    nr = 1 if case.k_joint == RIGID else 2
    axial_indices = [0] if case.kappa_horizontal==RIGID else [0,2+nr]
    aa = a0[np.ix_(axial_indices,axial_indices)]
    q1 = np.zeros(len(q0))
    q1[axial_indices] = np.linalg.solve(aa,-(a1@q0)[axial_indices])
    c2 = float(q0@(a2@q0+a1@q1))
    compliance = m.section.EA-m.J*omega**2
    qa = omega*math.sqrt(m.section.rhoA/compliance)
    ca = math.tan(2*m.l*qa)/(compliance*qa)
    out = {'C2_direct':c2,'mass_norm':inner(m,modal,modal),
           'axial_dynamic_compliance':ca,'axial_phase':2*m.l*qa,
           'axial_block_eigenvalues':np.linalg.eigvalsh(aa).tolist(),
           'first_order_residual':float(np.linalg.norm(a0@q1+a1@q0)),
           'q0':q0.tolist(),'q1':q1.tolist()}
    if case.kappa_horizontal == 0:
        gamma = q1[2+nr]+q0[2+nr+1]
        out['gamma0'] = float(gamma)
        out['axial_compliance_from_schur'] = float(np.linalg.inv(aa)[-1,-1])
    return out
