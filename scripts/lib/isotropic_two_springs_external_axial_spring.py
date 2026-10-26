"""Local axial spring at arm 2 outer end; both joint translations stay exact.

Canonical extension: docs/theory/isotropic_two_springs_external_axial_spring.md.
Only the outer axial row, admissible traces and spring energy change.
Existing arm fields, pole count and root isolation are reused privately.
"""
from dataclasses import dataclass
import importlib.util
import math
import sys
import numpy as np
from scipy.linalg import cholesky,eigh,null_space,solve_triangular
from scripts.lib import isotropic_two_springs_rl_timoshenko as base

Model=base.Model
RIGID=base.RIGID
joint_matrix=base.joint_matrix
inner=base.inner


@dataclass(frozen=True)
class Case(base.Case):
    chi_u: float | str=RIGID

    def __post_init__(self):
        super().__post_init__()
        if self.chi_u!=RIGID and (not math.isfinite(self.chi_u) or self.chi_u<0):
            raise ValueError('chi_u must be finite nonnegative or exact RIGID')


def external_matrix(m,case,arm):
    row=base.external_matrix(m,case,arm).copy()
    if arm==1 and case.chi_u!=RIGID:
        row[0]=0.;row[0,0]=case.chi_u*m.section.EA/m.l;row[0,3]=-1.
    return row


def trace_map(m,case):
    p,spring=base.trace_map(m,case)
    if case.chi_u==RIGID:return p,spring
    # q=(q10,q1L,q20,q2L); append the LOCAL u2(0), not a global displacement.
    p=np.column_stack((p,np.eye(12)[:,6]));s=np.zeros((p.shape[1],p.shape[1]))
    s[:-1,:-1]=spring;s[-1,-1]=case.chi_u*m.section.EA/m.l
    return p,s


_spec=importlib.util.spec_from_file_location(__name__+'._engine',base.__file__)
_engine=importlib.util.module_from_spec(_spec);sys.modules[_spec.name]=_engine
_spec.loader.exec_module(_engine)
_engine.external_matrix=external_matrix
_engine.trace_map=trace_map
prefix=_engine.prefix
count_at=_engine.count_at
boundary=_engine.boundary


def energy_components(m,case,arms):
    s=m.section
    rod=sum(np.array([a['weights']@(s.EA*a['dy'][:,0]**2),
        a['weights']@(s.EI*a['dy'][:,2]**2),a['weights']@(s.KGA*(a['dy'][:,1]-a['y'][:,2])**2)]) for a in arms)
    jump=arms[0]['ends'][1,2]-arms[1]['ends'][1,2]
    joint=0. if case.k_joint==RIGID else case.k_joint*s.EI/m.l*jump**2
    end=0. if case.k_end==RIGID else case.k_end*s.EI/m.l*arms[1]['ends'][0,2]**2
    axial_end=0. if case.chi_u==RIGID else case.chi_u*s.EA/m.l*arms[1]['ends'][0,0]**2
    return np.r_[rod,joint,end,axial_end]  # quadratic energies, twice V


def mode(m,case,omega,order=120):
    result=_engine.mode(m,case,omega,order)
    if case.chi_u!=RIGID:
        a=result['arms'];d=result['diagnostics'];energy=energy_components(m,case,a)
        d['energy_error']=float(abs(energy.sum()/omega**2-1))
        fscale=max(np.max(abs(v['y'][:,3:]*[1,1,1/m.l])) for v in a)
        d['external_scaled'][1][0]=float(abs(d['external_absolute'][1][0])/fscale)
        d['u20']=float(a[1]['ends'][0,0]);d['N20']=float(a[1]['ends'][0,3])
        d['outer_axial_spring_quadratic_energy']=float(energy[-1])
    return result


def ritz(m,case,order,count=3):
    """Existing direct-energy space plus one constant u2 column when free."""
    rotation=int(case.k_end!=RIGID);axial=int(case.chi_u!=RIGID)
    size=6*order+rotation+axial;k=np.zeros((size,size));mass=np.zeros_like(k)
    traces=[];nodes,weights=np.polynomial.legendre.leggauss(max(64,order+2))
    b=base.shifted_legendre_clamped_basis(order,(nodes+1)/2)
    for i,L in enumerate(case.lengths(m)):
        v=np.zeros((len(nodes),3,size));d=np.zeros_like(v)
        end=np.zeros((3,size));outer=np.zeros_like(end)
        for field in range(3):
            sl=slice((3*i+field)*order,(3*i+field+1)*order)
            v[:,field,sl],d[:,field,sl]=b.values.T,b.derivatives.T/L;end[field,sl]=1.
        if i==1:
            if rotation:v[:,2,6*order]=end[2,6*order]=outer[2,6*order]=1.
            if axial:v[:,0,-1]=end[0,-1]=outer[0,-1]=1.
        wt=weights*L/2;s=m.section
        def gram(a,b):return a.T@(wt[:,None]*b)
        k+=s.EA*gram(d[:,0],d[:,0])+s.EI*gram(d[:,2],d[:,2])+s.KGA*gram(d[:,1]-v[:,2],d[:,1]-v[:,2])
        mass+=s.rhoA*(gram(v[:,0],v[:,0])+gram(v[:,1],v[:,1]))+s.rhoI*gram(v[:,2],v[:,2])+m.J*gram(d[:,0],d[:,0])
        traces.append(end)
        if i==1:
            if rotation:k+=case.k_end*s.EI/m.l*np.outer(outer[2],outer[2])
            if axial:k+=case.chi_u*s.EA/m.l*np.outer(outer[0],outer[0])
    g=[np.column_stack((f.t,f.n)) for f in base.geometry.frames(case.beta)]
    constraint=g[0]@traces[0][:2]-g[1]@traces[1][:2]
    jump=traces[0][2]-traces[1][2]
    if case.k_joint==RIGID:constraint=np.vstack((constraint,jump))
    else:k+=case.k_joint*m.section.EI/m.l*np.outer(jump,jump)
    Z=null_space(constraint);Kr,Mr=Z.T@k@Z,Z.T@mass@Z
    L=cholesky((Kr+Kr.T)/2,lower=True);cholesky((Mr+Mr.T)/2,lower=True)
    a=solve_triangular(L,Mr,lower=True);a=solve_triangular(L,a.T,lower=True).T
    val,v=eigh((a+a.T)/2,subset_by_index=[len(a)-count,len(a)-1]);val,v=val[::-1],v[:,::-1]
    coef=Z@solve_triangular(L.T,v,lower=False)/np.sqrt(val)
    return np.sqrt(1/val),coef,dict(mass_gram_error=float(np.max(abs(coef.T@mass@coef-np.eye(count)))),
        energy_error=float(np.max(abs(np.diag(coef.T@k@coef)*val-1))),
        constraint_residual=float(np.max(abs(constraint@coef))),positive_K_M=True)


def leading_scale(chi_u):
    """xi_eff/xi for the unchanged leading Ritz: xi_eff²/2=xi²/(2+1/chi)."""
    if chi_u==RIGID:return 1.
    if not math.isfinite(chi_u) or chi_u<0:raise ValueError('Invalid axial stiffness')
    return math.sqrt(2*chi_u/(1+2*chi_u))
