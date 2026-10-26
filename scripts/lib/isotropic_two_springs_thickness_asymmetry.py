"""Individual rectangular sections, fixed volume and reference spring stiffness.

Canonical extension: docs/theory/isotropic_two_springs_thickness_asymmetry.md.
Unchanged state (u,w,psi,N,Q,M), inward coordinates, exact translations.
Reuses the original arm functions and its count/Brent isolation in a private
namespace, as in the existing slip adapter; never mutates the parent module.
"""
from dataclasses import dataclass, replace
from functools import cached_property
import importlib.util
import math
import sys
import numpy as np
from scipy.linalg import block_diag, cholesky, eigh, null_space, solve_triangular
from scripts.lib import isotropic_two_springs_rl_timoshenko as base
from scripts.lib import isotropic_two_springs_leading as leading


def ratios(mu,tau):
    if not abs(mu)<1 or not abs(tau)<1:
        raise ValueError('Positive lengths and thicknesses require |mu|,|tau|<1')
    d=1+tau*tau+2*mu*tau
    h=np.array([1-tau,1+tau])/math.sqrt(d)
    area=h*h
    return h,area,area*area


@dataclass(frozen=True)
class Model(base.Model):
    mu: float = 0.
    tau: float = 0.

    def __post_init__(self):
        super().__post_init__()
        ratios(self.mu,self.tau)
        if self.variant!='RL_FULL' or self.aspect!=4:
            raise ValueError('This study is RL_FULL with rectangular aspect 4 only')

    @cached_property
    def arms(self):
        return tuple(base.Model(epsilon=self.epsilon*t,nu=self.nu,aspect=self.aspect,
                    E=self.E,rho=self.rho,l=self.l,ks=self.ks,variant=self.variant)
                     for t in ratios(self.mu,self.tau)[0])

    @property
    def S_ax(self):
        area=ratios(self.mu,self.tau)[1]
        return (1-self.mu)/area[0]+(1+self.mu)/area[1]


@dataclass(frozen=True)
class Case(base.Case):
    tau: float = 0.

    def __post_init__(self):
        super().__post_init__()
        ratios(self.mu,self.tau)


def compatible(m,case):
    if m.mu!=case.mu or m.tau!=case.tau:
        raise ValueError('Model/case geometry mismatch')


# These operators act on physical traces, and spring coefficients use EI0.
joint_matrix=base.joint_matrix
external_matrix=base.external_matrix
trace_map=base.trace_map
balance=base.balance


def boundary(m,case,omega):
    compatible(m,case)
    ends=[base.basis(a,L,omega,[0.,L]) for a,L in zip(m.arms,case.lengths(m))]
    b=np.zeros((12,12))
    for i,v in enumerate(ends):
        b[3*i:3*i+3,6*i:6*i+6]=external_matrix(m,case,i)@v[0]
    b[6:]=joint_matrix(m,case)@block_diag(*(v[1] for v in ends))
    return b


def pole_key(arm,length):
    # E,rho,nu,aspect,ks are fixed by the enclosing immutable study contract.
    return f'{length:.17g}|{arm.epsilon:.17g}'


def count_at(m,case,omega,catalogs,policy):
    compatible(m,case)
    poles=[]
    for arm,L in zip(m.arms,case.lengths(m)):
        cat=catalogs[pole_key(arm,L)]
        if omega>=cat['upper']:
            raise ArithmeticError('Query outside verified arm-pole interval')
        poles+=cat['poles']
        d=arm.section.EA-arm.J*cat['upper']**2
        if d<=0:raise ArithmeticError('Outside individual RL acoustic interval')
        nmax=math.floor(cat['upper']*math.sqrt(arm.section.rhoA/d)*L/math.pi)
        poles += [math.sqrt(arm.section.EA*(n*math.pi/L)**2 /
                    (arm.section.rhoA+arm.J*(n*math.pi/L)**2)) for n in range(1,nmax+1)]
    w=float(omega)
    exclusion=policy['pole_exclusion']*max(1.,w)
    for pole in sorted(poles):
        if abs(w-pole)<=exclusion:w=pole+2*exclusion
    pairs=[base.arm_dtn(a,L,w) for a,L in zip(m.arms,case.lengths(m))]
    p,spring=trace_map(m,case)
    a=p.T@block_diag(*(r[0] for r in pairs))@p+spring
    scale=np.sqrt(m.l/m.section.EI)*np.array([m.l,m.l]+[1.]*(len(a)-2))
    a=a*scale[:,None]*scale
    skew=np.linalg.norm(a-a.T)/np.linalg.norm(a)
    vals=np.linalg.eigvalsh((a+a.T)/2)
    margin=min(abs(vals))/max(abs(vals))
    if skew>policy['schur_tol'] or margin<policy['inertia_margin'] or max(r[1] for r in pairs)>policy['condition_max']:
        raise ArithmeticError(f'Unresolved individual-section DtN inertia: {skew}, {margin}')
    j0=sum(p<w for p in poles)
    return w,j0+int(sum(vals<0)),dict(requested=omega,effective=w,J0=j0,
            eigenvalues=vals.tolist(),skew=float(skew),margin=float(margin),conditions=[r[1] for r in pairs])


_spec=importlib.util.spec_from_file_location(__name__+'._engine',base.__file__)
_engine=importlib.util.module_from_spec(_spec)
sys.modules[_spec.name]=_engine
_spec.loader.exec_module(_engine)
_engine.boundary=boundary
_engine.count_at=count_at
prefix=_engine.prefix


def energies(m,case,arms):
    mass_components=[];elastic=[]
    for a,arm in zip(arms,m.arms):
        y,dy,wt=a['y'],a['dy'],a['weights'];s=arm.section
        mass_components.append([wt@(s.rhoA*y[:,0]**2),wt@(s.rhoA*y[:,1]**2),
                                wt@(s.rhoI*y[:,2]**2),wt@(arm.J*dy[:,0]**2)])
        elastic.append([wt@(s.EA*dy[:,0]**2),wt@(s.EI*dy[:,2]**2),wt@(s.KGA*(dy[:,1]-y[:,2])**2)])
    springs=[]
    delta=arms[0]['ends'][1,2]-arms[1]['ends'][1,2]
    for k,v in ((case.k_joint,delta),(case.k_end,arms[1]['ends'][0,2])):
        springs.append(0. if k==base.RIGID else k*m.section.EI/m.l*v*v)
    return np.array(mass_components),np.array(elastic),np.array(springs)


def mode(m,case,omega,order=120):
    compatible(m,case)
    b,_,col=balance(boundary(m,case,omega))
    _,singular,vh=np.linalg.svd(b)
    coef=(col*vh[-1]).reshape(2,6)
    nodes,weights=np.polynomial.legendre.leggauss(order)
    arms=[]
    for arm,L,c in zip(m.arms,case.lengths(m),coef):
        x,wt=(nodes+1)*L/2,weights*L/2
        y,dy=(base.basis(arm,L,omega,x,d)@c for d in (0,1))
        ends=base.basis(arm,L,omega,[0.,L])@c
        arms.append(dict(x=x,weights=wt,y=y,dy=dy,ends=ends))
    mc,ec,sp=energies(m,case,arms)
    mass=float(mc.sum());energy=float(ec.sum()+sp.sum())
    if mass<=0 or energy<=0:raise ArithmeticError('Nonpositive energy/mass')
    for a in arms:
        for key in ('y','dy','ends'):a[key]/=math.sqrt(mass)
    coef/=math.sqrt(mass)
    qscale=max(np.max(abs(a['y'][:,:3]*[1,1,m.l])) for a in arms)
    fscale=max(np.max(abs(a['y'][:,3:]*[1,1,1/m.l])) for a in arms)
    jr=joint_matrix(m,case)@np.concatenate([a['ends'][1] for a in arms])
    js=abs(jr)/np.array([qscale,qscale,fscale,fscale,
                         qscale/m.l if case.k_joint==base.RIGID else fscale*m.l,fscale*m.l])
    er=[external_matrix(m,case,i)@a['ends'][0] for i,a in enumerate(arms)]
    es=[abs(r)/[qscale,qscale,qscale/m.l if i==0 or case.k_end==base.RIGID else fscale*m.l] for i,r in enumerate(er)]
    scales=np.array([1,1,m.l,qscale/fscale,qscale/fscale,qscale/(fscale*m.l)])
    ode=max(np.max(abs((a['dy']-a['y']@base.state_matrix(arm,omega).T)*scales)) /
            max(np.max(abs(a['dy']*scales)),1e-30) for a,arm in zip(arms,m.arms))
    return dict(coefficients=coef.tolist(),arms=arms,diagnostics=dict(
        mass=mass/math.sqrt(mass)**2,raw_mass=mass,energy_error=abs(energy/mass/omega**2-1),
        singular_ratio=float(singular[-1]/singular[0]),nonzero_condition=float(singular[0]/singular[-2]),
        joint_scaled=js.tolist(),external_scaled=np.array(es).tolist(),joint_absolute=jr.tolist(),
        external_absolute=np.array(er).tolist(),ode_error=float(ode),
        mass_components=(mc/mass).tolist(),elastic_components=(ec/mass).tolist(),spring_components=(sp/mass).tolist()))


def inner(m,first,second):
    value=0.
    for a,b,arm in zip(first['arms'],second['arms'],m.arms):
        v,w=a['y'],b['y'];s=arm.section
        value+=a['weights']@(s.rhoA*(v[:,:2]*w[:,:2]).sum(axis=1)+s.rhoI*v[:,2]*w[:,2]+arm.J*a['dy'][:,0]*b['dy'][:,0])
    return float(value)


def ritz(m,case,order,count=3):
    """Direct energies, independent of transfer bases/count; original trial space."""
    compatible(m,case)
    size=6*order+int(case.k_end!=base.RIGID)
    K,M=np.zeros((size,size)),np.zeros((size,size));traces=[]
    nodes,weights=np.polynomial.legendre.leggauss(max(64,order+2))
    b=base.shifted_legendre_clamped_basis(order,(nodes+1)/2)
    for i,(L,arm) in enumerate(zip(case.lengths(m),m.arms)):
        v,d=np.zeros((len(nodes),3,size)),np.zeros((len(nodes),3,size))
        end,outer=np.zeros((3,size)),np.zeros((3,size))
        for f in range(3):
            sl=slice((3*i+f)*order,(3*i+f+1)*order)
            v[:,f,sl],d[:,f,sl]=b.values.T,b.derivatives.T/L
            end[f,sl]=1.
        if i==1 and case.k_end!=base.RIGID:v[:,2,-1]=end[2,-1]=outer[2,-1]=1.
        wt=weights*L/2;s=arm.section
        def gram(a,b):return a.T@(wt[:,None]*b)
        K+=s.EA*gram(d[:,0],d[:,0])+s.EI*gram(d[:,2],d[:,2])+s.KGA*gram(d[:,1]-v[:,2],d[:,1]-v[:,2])
        M+=s.rhoA*(gram(v[:,0],v[:,0])+gram(v[:,1],v[:,1]))+s.rhoI*gram(v[:,2],v[:,2])+arm.J*gram(d[:,0],d[:,0])
        traces.append(end)
        if i==1 and case.k_end!=base.RIGID:K+=case.k_end*m.section.EI/m.l*np.outer(outer[2],outer[2])
    g=[np.column_stack((f.t,f.n)) for f in base.geometry.frames(case.beta)]
    constraint=g[0]@traces[0][:2]-g[1]@traces[1][:2]
    jump=traces[0][2]-traces[1][2]
    if case.k_joint==base.RIGID:constraint=np.vstack((constraint,jump))
    else:K+=case.k_joint*m.section.EI/m.l*np.outer(jump,jump)
    Z=null_space(constraint);Kr,Mr=Z.T@K@Z,Z.T@M@Z
    L=cholesky((Kr+Kr.T)/2,lower=True);cholesky((Mr+Mr.T)/2,lower=True)
    A=solve_triangular(L,Mr,lower=True);A=solve_triangular(L,A.T,lower=True).T
    val,v=eigh((A+A.T)/2,subset_by_index=[len(A)-count,len(A)-1]);val,v=val[::-1],v[:,::-1]
    c=Z@solve_triangular(L.T,v,lower=False)/np.sqrt(val)
    return np.sqrt(1/val),c,dict(mass_gram_error=float(np.max(abs(c.T@M@c-np.eye(count)))),
        energy_error=float(np.max(abs(np.diag(c.T@K@c)*val-1))),constraint_residual=float(np.max(abs(constraint@c))),positive_K_M=True)


class Leading(leading.Ritz):
    """Weighted version of the same integrated-Legendre weak problem."""
    def __init__(self,k_joint,k_end,mu=0.,tau=0.,order=20,quadrature=64):
        self.tau=tau
        super().__init__(k_joint,k_end,mu,order,quadrature)
        _,self.area,self.inertia=ratios(mu,tau)
        self.S_ax=sum(L/a for L,a in zip(self.lengths,self.area))
        self.arm_K=[];self.arm_M=[]
        for i in range(2):
            sl=slice(i*order,(i+1)*order)
            K,M=np.zeros_like(self.Kb),np.zeros_like(self.M)
            K[sl,sl]=self.Kb[sl,sl];M[sl,sl]=self.M[sl,sl]
            self.arm_K.append(K);self.arm_M.append(M)
        self.Kb=sum(j*K for j,K in zip(self.inertia,self.arm_K))
        self.M=sum(a*M for a,M in zip(self.area,self.arm_M))
        self.K0=self.Kb+self.Kj+self.Ke
        self.matrix_checks['symmetry']=max(float(np.max(abs(A-A.T))) for A in (self.K0,self.M,self.P))
        cholesky(self.Z.T@self.M@self.Z,lower=True)

    def solve(self,xi,count=3,node_fixed=False):
        # Parent solves the identical matrix K0 + t P after rescaling its xi.
        # Its unweighted strong residuals are inapplicable and are REPLACED below.
        result=super().solve(xi*math.sqrt(sum(self.lengths)/self.S_ax),count,node_fixed)
        c=result['coefficients'];zeta=result['zeta'];t=0. if node_fixed else xi*xi/self.S_ax
        natural=[];ode=[];components=[]
        for j in range(count):
            v=c[:,j];end=[[e@v for e in a] for a in self.ends];outer=[[e@v for e in a] for a in self.outer]
            i1,i2=self.inertia;a1,a2=self.area
            residual=np.array([-i2*outer[1][2]+self.k_end*outer[1][1],
                i1*end[0][2]+self.k_joint*(self.jump@v),i2*end[1][2]-self.k_joint*(self.jump@v),
                i1*end[0][3]-i2*end[1][3]+t*(self.y@v)])
            if node_fixed:residual=residual[:-1]
            norm=max(1.,np.max(abs(np.array(end))),np.max(abs(np.array(outer))),abs(t*(self.y@v)))
            natural.append(float(max(abs(residual))/norm));components.append((residual/norm).tolist())
            num=den=0.
            for i,L in enumerate(self.lengths):
                f=leading.fields(self.order,self.nodes,L,i);ci=v[i*self.order:(i+1)*self.order]
                A=self.inertia[i]*(f[4]@ci);B=zeta[j]*self.area[i]*(f[0]@ci)
                num+=self.weights@((A-B)**2)*L;den+=self.weights@(A*A+B*B)*L
            ode.append(float(math.sqrt(num/den)))
        result['diagnostics'].update(natural_boundary=max(natural),strong_ode=max(ode),natural_components=components)
        y=result['node_value']
        result.update(s_xi=xi*y*y/(self.S_ax*zeta),C_at_zero=y[0]**2/(2*self.S_ax*zeta[0]) if xi==0 and not node_fixed else None)
        return result

    def tau_derivative_at_zero(self,result):
        if self.tau!=0:raise ValueError('The derivative formula is local to tau=0')
        da=np.array([-2*(1+self.mu),2*(1-self.mu)]);di=2*da
        c=result['coefficients'][:,0];zeta=result['zeta'][0];mass=c@self.M@c
        bending=sum(v*(c@K@c) for v,K in zip(di,self.arm_K))/mass
        redistribution=-zeta*sum(v*(c@M@c) for v,M in zip(da,self.arm_M))/mass
        return dict(bending=float(bending),mass=float(redistribution),axial_direct=0.,
                    total=float(bending+redistribution),relative_log_frequency=float((bending+redistribution)/(2*zeta)))
