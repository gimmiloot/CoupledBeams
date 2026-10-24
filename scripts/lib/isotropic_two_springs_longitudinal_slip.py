"""Directed longitudinal compliance: a separate physical branch of LONG.

No executable entry point. The immutable parent engine is loaded into a PRIVATE
module namespace, so its existing count isolation, basis, and mode recovery can
use branch-specific physical hooks without mutating the canonical module.
The only physical changes are c1 -> k_parallel*c1 and its energy/trace space.
"""
from dataclasses import dataclass
import importlib.util
import math
import sys
import numpy as np
from scipy.linalg import null_space, cholesky, solve_triangular, eigh
from scripts.lib import isotropic_two_springs_rl_timoshenko as base

Model = base.Model
RIGID = base.RIGID


@dataclass(frozen=True)
class Case(base.Case):
    chi: float | str = RIGID

    def __post_init__(self):
        super().__post_init__()
        if self.chi != RIGID and (not math.isfinite(self.chi) or self.chi < 0):
            raise ValueError('chi must be nonnegative finite or exact RIGID')


def joint_matrix(m, case):
    a = base.joint_matrix(m, case).copy()
    if case.chi != RIGID:
        # N1 + k_parallel*c1 = 0; retain c2 and both force balances.
        a[0] *= case.chi*m.section.EA/m.l
        a[0, 3] += 1.
    return a


def trace_map(m, case):
    if case.chi == RIGID:
        return base.trace_map(m, case)
    if case.k_joint == RIGID or case.k_end == RIGID:
        raise ValueError('Bounded slip branch requires finite rotational springs')
    # q=(q10,q1L,q20,q2L), independent (u1,u2,w2,psi1,psi2,psi20).
    c,s = math.cos(math.radians(case.beta)),math.sin(math.radians(case.beta))
    p = np.zeros((12,6));p[3,0]=1;p[4,1:3]=[s,-c]
    p[9,1]=p[10,2]=p[5,3]=p[11,4]=p[8,5]=1.
    jump,end = p[5]-p[11],p[8];slip=p[3]+c*p[9]+s*p[10]
    spring = (case.k_joint*m.section.EI/m.l*np.outer(jump,jump)
        +case.k_end*m.section.EI/m.l*np.outer(end,end)
        +case.chi*m.section.EA/m.l*np.outer(slip,slip))
    return p,spring


_spec=importlib.util.spec_from_file_location(__name__+'._engine',base.__file__)
_engine=importlib.util.module_from_spec(_spec)
sys.modules[_spec.name]=_engine
_spec.loader.exec_module(_engine)
_engine.joint_matrix=joint_matrix
_engine.trace_map=trace_map
prefix=_engine.prefix
count_at=_engine.count_at
boundary=_engine.boundary
inner=base.inner


def mode(m, case, omega, order=120):
    result=_engine.mode(m,case,omega,order)
    if case.chi != RIGID:
        a,b=[arm['ends'][1] for arm in result['arms']]
        c,s=math.cos(math.radians(case.beta)),math.sin(math.radians(case.beta))
        slip=a[0]+c*b[0]+s*b[1]
        sec=m.section
        rod=sum(arm['weights']@(sec.EA*arm['dy'][:,0]**2+sec.EI*arm['dy'][:,2]**2
             +sec.KGA*(arm['dy'][:,1]-arm['y'][:,2])**2) for arm in result['arms'])
        energy=rod+case.k_joint*sec.EI/m.l*(a[2]-b[2])**2+case.k_end*sec.EI/m.l*result['arms'][1]['ends'][0,2]**2+case.chi*sec.EA/m.l*slip**2
        d=result['diagnostics'];d['energy_error']=float(abs(energy/omega**2-1))
        fscale=max(np.max(abs(arm['y'][:,3:]*[1,1,1/m.l])) for arm in result['arms'])
        d['joint_scaled'][0]=abs(d['joint_absolute'][0])/fscale
        d['longitudinal_slip']=float(slip)
    return result


def ritz(m, case, order, count=3):
    """Independent direct energies; reuse the existing clamped polynomial basis.

    The extra constant psi2 column is present ONLY for finite end springs.
    Force/moment equations are not imposed as Ritz constraints.
    """
    from scripts.lib.reddy_symmetric_coupled_beams_ritz import shifted_legendre_clamped_basis
    geometry=base.geometry
    size = 6*order+int(case.k_end != RIGID)
    k, mass = np.zeros((size, size)), np.zeros((size, size))
    traces = []
    nodes, weights = np.polynomial.legendre.leggauss(max(64, order+2))
    b = shifted_legendre_clamped_basis(order, (nodes+1)/2)
    for i, L in enumerate(case.lengths(m)):
        v, d = np.zeros((len(nodes), 3, size)), np.zeros((len(nodes), 3, size))
        end, outer = np.zeros((3, size)), np.zeros((3, size))
        for f in range(3):
            sl = slice((3*i+f)*order, (3*i+f+1)*order)
            v[:, f, sl], d[:, f, sl] = b.values.T, b.derivatives.T/L
            end[f, sl] = 1.
        if i == 1 and case.k_end != RIGID:
            v[:, 2, -1] = end[2, -1] = outer[2, -1] = 1.
        wt = weights*L/2
        def gram(a, b):
            return a.T@(wt[:, None]*b)
        s = m.section
        k += s.EA*gram(d[:, 0], d[:, 0])+s.EI*gram(d[:, 2], d[:, 2])+s.KGA*gram(d[:, 1]-v[:, 2], d[:, 1]-v[:, 2])
        mass += s.rhoA*(gram(v[:, 0], v[:, 0])+gram(v[:, 1], v[:, 1]))+s.rhoI*gram(v[:, 2], v[:, 2])+m.J*gram(d[:, 0], d[:, 0])
        traces.append(end)
        if i == 1 and case.k_end != RIGID:
            k += case.k_end*s.EI/m.l*np.outer(outer[2], outer[2])
    g = [np.column_stack((f.t, f.n)) for f in geometry.frames(case.beta)]
    constraint = g[0]@traces[0][:2]-g[1]@traces[1][:2]
    if case.chi != RIGID:
        slip=constraint[0].copy()
        k += case.chi*m.section.EA/m.l*np.outer(slip,slip)
        constraint=constraint[1:2]
    jump = traces[0][2]-traces[1][2]
    if case.k_joint == RIGID:
        constraint = np.vstack((constraint, jump))
    else:
        k += case.k_joint*m.section.EI/m.l*np.outer(jump, jump)
    z = null_space(constraint)
    kr, mr = z.T@k@z, z.T@mass@z
    lk = cholesky((kr+kr.T)/2, lower=True)
    cholesky((mr+mr.T)/2, lower=True)
    a = solve_triangular(lk, mr, lower=True)
    a = solve_triangular(lk, a.T, lower=True).T
    val, v = eigh((a+a.T)/2, subset_by_index=[len(a)-count, len(a)-1])
    val, v = val[::-1], v[:, ::-1]
    coef = z@solve_triangular(lk.T, v, lower=False)/np.sqrt(val)
    return np.sqrt(1/val), coef, {"mass_gram_error": float(np.max(abs(coef.T@mass@coef-np.eye(count)))),
        "energy_error": float(np.max(abs(np.diag(coef.T@k@coef)*val-1))),
        "constraint_residual": float(np.max(abs(constraint@coef))), "positive_K_M": True}
