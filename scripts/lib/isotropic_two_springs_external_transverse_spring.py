"""Local transverse spring at arm 2 outer end; exact axial support and joint.

Canonical extension: docs/theory/isotropic_two_springs_external_transverse_spring.md.
No change to arm equations, coordinate frames or rotational springs.
"""
from dataclasses import dataclass
import importlib.util
import math
import sys
import numpy as np
from scipy.linalg import cholesky,eigh,null_space,solve_triangular
from scripts.lib import isotropic_two_springs_rl_timoshenko as base
from scripts.lib import isotropic_two_springs_external_axial_spring as axial_support
from scripts.lib import isotropic_two_springs_leading as leading

Model=base.Model
RIGID=base.RIGID
joint_matrix=base.joint_matrix
inner=base.inner


@dataclass(frozen=True)
class Case(base.Case):
    k_w: float | str=RIGID  # k_w dimensional * l^3/(EI)

    def __post_init__(self):
        super().__post_init__()
        if self.k_w!=RIGID and (not math.isfinite(self.k_w) or self.k_w<0):
            raise ValueError('k_w must be finite nonnegative or exact RIGID')


def external_matrix(m,case,arm):
    row=base.external_matrix(m,case,arm).copy()
    if arm==1 and case.k_w!=RIGID:
        row[1]=0.;row[1,1]=case.k_w*m.section.EI/m.l**3;row[1,4]=-1.
    return row


def trace_map(m,case):
    p,spring=base.trace_map(m,case)
    if case.k_w==RIGID:return p,spring
    p=np.column_stack((p,np.eye(12)[:,7]))  # local w2(0), u2(0) stays zero
    s=np.zeros((p.shape[1],p.shape[1]));s[:-1,:-1]=spring
    s[-1,-1]=case.k_w*m.section.EI/m.l**3
    return p,s


_spec=importlib.util.spec_from_file_location(__name__+'._engine',base.__file__)
_engine=importlib.util.module_from_spec(_spec);sys.modules[_spec.name]=_engine
_spec.loader.exec_module(_engine)
_engine.external_matrix=external_matrix;_engine.trace_map=trace_map
prefix=_engine.prefix
count_at=_engine.count_at
boundary=_engine.boundary


def energy_components(m,case,arms):
    # Reuse the checked rod and rotational energies with NO axial spring.
    rigid_axial=axial_support.Case(case.beta,case.mu,case.k_joint,case.k_end,RIGID)
    e=axial_support.energy_components(m,rigid_axial,arms)
    e[-1]=0. if case.k_w==RIGID else case.k_w*m.section.EI/m.l**3*arms[1]['ends'][0,1]**2
    return e  # twice V, final entry is transverse spring energy


def mode(m,case,omega,order=120):
    result=_engine.mode(m,case,omega,order)
    if case.k_w!=RIGID:
        arms=result['arms'];d=result['diagnostics'];energy=energy_components(m,case,arms)
        d['energy_error']=float(abs(energy.sum()/omega**2-1))
        fscale=max(np.max(abs(a['y'][:,3:]*[1,1,1/m.l])) for a in arms)
        d['external_scaled'][1][1]=float(abs(d['external_absolute'][1][1])/fscale)
        d['w20']=float(arms[1]['ends'][0,1]);d['Q20']=float(arms[1]['ends'][0,4])
        d['outer_transverse_quadratic_energy']=float(energy[-1])
    return result


def ritz(m,case,order,count=3):
    """Canonical direct-energy polynomial space plus constant w2 when free."""
    rotation=int(case.k_end!=RIGID);transverse=int(case.k_w!=RIGID)
    size=6*order+rotation+transverse;k=np.zeros((size,size));mass=np.zeros_like(k)
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
            if transverse:v[:,1,-1]=end[1,-1]=outer[1,-1]=1.
        wt=weights*L/2;s=m.section
        def gram(a,b):return a.T@(wt[:,None]*b)
        k+=s.EA*gram(d[:,0],d[:,0])+s.EI*gram(d[:,2],d[:,2])+s.KGA*gram(d[:,1]-v[:,2],d[:,1]-v[:,2])
        mass+=s.rhoA*(gram(v[:,0],v[:,0])+gram(v[:,1],v[:,1]))+s.rhoI*gram(v[:,2],v[:,2])+m.J*gram(d[:,0],d[:,0])
        traces.append(end)
        if i==1:
            if rotation:k+=case.k_end*s.EI/m.l*np.outer(outer[2],outer[2])
            if transverse:k+=case.k_w*s.EI/m.l**3*np.outer(outer[1],outer[1])
    g=[np.column_stack((f.t,f.n)) for f in base.geometry.frames(case.beta)]
    constraint=g[0]@traces[0][:2]-g[1]@traces[1][:2];jump=traces[0][2]-traces[1][2]
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


@dataclass
class Leading(leading.Ritz):
    k_w: float | str=RIGID
    # Opt-in numerical assembly verified by the single-case strong-residual
    # diagnostic. The old default remains available for archive reproduction.
    analytic_curvature: bool=False

    def __post_init__(self):
        if self.k_w!=RIGID and (not math.isfinite(self.k_w) or self.k_w<0):raise ValueError('Invalid transverse spring')
        super().__post_init__()
        if self.analytic_curvature:
            # Integral_0^L B_n'' B_m'' ds = delta_nm/((2n+1)*L**3).
            # This is exactly the existing integrated-Legendre space; avoid
            # quadrature roundoff exciting tiny high-degree coefficients.
            self.Kb[:]=0.
            for arm,L in enumerate(self.lengths):
                offset=int(arm==1)
                for n in range(self.order-offset):
                    index=arm*self.order+offset+n
                    self.Kb[index,index]=1./((2*n+1)*L**3)
            self.K0=self.Kb+self.Kj+self.Ke
            self.matrix_checks['symmetry']=max(float(np.max(abs(a-a.T))) for a in (self.K0,self.M,self.P))
        if self.k_w==RIGID:return
        # Preserve the old integrated basis and add one constant on arm 2.
        # Its curvature vanishes, but all consistent mass cross terms remain.
        self.M=np.pad(self.M,((0,1),(0,1)));self.Kb=np.pad(self.Kb,((0,1),(0,1)))
        L=self.lengths[1];f=leading.fields(self.order,self.nodes,L,1)[0]
        cross=f.T@(self.weights*L);sl=slice(self.order,2*self.order)
        self.M[sl,-1]=cross;self.M[-1,sl]=cross;self.M[-1,-1]=L
        for collection in (self.ends,self.outer):
            for i in range(2):
                for derivative in range(4):collection[i][derivative]=np.r_[collection[i][derivative],float(i==1 and derivative==0)]
        self.jump=self.ends[0][1]-self.ends[1][1];self.end_rotation=self.outer[1][1]
        self.y=self.ends[1][0];self.P=np.outer(self.y,self.y)
        self.Kj=self.k_joint*np.outer(self.jump,self.jump);self.Ke=self.k_end*np.outer(self.end_rotation,self.end_rotation)
        self.Kw=self.k_w*np.outer(self.outer[1][0],self.outer[1][0]);self.K0=self.Kb+self.Kj+self.Ke+self.Kw
        self.constraint=(self.ends[0][0]+self.ends[1][0])[None,:];self.Z=null_space(self.constraint)
        cholesky(self.Z.T@self.M@self.Z,lower=True)
        self.matrix_checks.update(symmetry=max(float(np.max(abs(a-a.T))) for a in (self.K0,self.M,self.P)),
            constraint_rank=int(np.linalg.matrix_rank(self.constraint)),P_rank=int(np.linalg.matrix_rank(self.P)),dimension=self.Z.shape[1])

    def solve(self,xi,count=3,node_fixed=False):
        if node_fixed:raise ValueError('Node-fixed exploration is outside this study')
        result=super().solve(xi,count=count)
        if self.k_w==RIGID:return result
        # The inherited eigensolve uses the augmented matrices. Its original
        # field diagnostics assume w2(0)=0; replace only those diagnostics by
        # the complete augmented fields and the newly derived natural row.
        natural=[];strong=[];essential=[];outer_force=[];t=xi*xi/2
        for j,zeta in enumerate(result['zeta']):
            c=result['coefficients'][:,j];end=np.array([[e@c for e in arm] for arm in self.ends]);outer=np.array([[e@c for e in arm] for arm in self.outer])
            residual=np.array([outer[1,3]+self.k_w*outer[1,0],-outer[1,2]+self.k_end*outer[1,1],
                end[0,2]+self.k_joint*(self.jump@c),end[1,2]-self.k_joint*(self.jump@c),end[0,3]-end[1,3]+t*(self.y@c)])
            scale=max(1.,np.max(abs(end)),np.max(abs(outer)),abs(t*(self.y@c)),abs(self.k_w*outer[1,0]))
            natural.append(float(np.max(abs(residual))/scale));outer_force.append(float(abs(residual[0])/scale))
            essential.append(max(abs(outer[0,0]),abs(outer[0,1]),abs(end[0,0]+end[1,0])))
            num=den=0.
            for i,L in enumerate(self.lengths):
                f=leading.fields(self.order,self.nodes,L,i);ci=c[i*self.order:(i+1)*self.order]
                W=f[0]@ci+(c[-1] if i==1 else 0.);D4=f[4]@ci
                num+=self.weights@((D4-zeta*W)**2)*L;den+=self.weights@(D4**2+(zeta*W)**2)*L
            strong.append(float(np.sqrt(num/den)))
        result['diagnostics'].update(natural_boundary=max(natural),strong_ode=max(strong),essential=max(essential),outer_transverse_force=max(outer_force))
        return result
