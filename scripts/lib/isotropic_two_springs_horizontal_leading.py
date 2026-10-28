"""HS-L1 leading horizontal-spring problem; no full-model modifications.

K=Kb+kj*jump'*jump+tau*D'*D, D=W2(0)-W2(L2).
Only clamp 1 and W1(L1)+W2(L2)=0 are essential; saturation adds D=0.
Integrated Legendre fields are reused, but the trace space and natural rows
are built here independently of the old guided transverse-spring problem.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.linalg import cholesky, eigh, null_space, solve_triangular, expm
from scipy.optimize import brentq
from scripts.lib.isotropic_two_springs_leading import fields

RIGID = 'RIGID'


@dataclass
class Leading:
    k_joint: float = .1
    mu: float = .01
    order: int = 24
    quadrature: int = 80

    def __post_init__(self):
        if self.k_joint <= 0 or abs(self.mu)>=1 or self.order<4:
            raise ValueError('HS-L1 finite positive internal spring contract')
        if self.quadrature < self.order+3:
            raise ValueError('Insufficient polynomial mass quadrature')
        self.lengths=(1-self.mu,1+self.mu); n=2*self.order+1
        t,wt=np.polynomial.legendre.leggauss(self.quadrature)
        self.nodes=(t+1)/2;self.weights=wt/2
        self.Kb=np.zeros((n,n));self.M=np.zeros((n,n))
        self.ends=[];self.outer=[];self.values=[]
        for arm,L in enumerate(self.lengths):
            values=self.field_matrices(arm,self.nodes)
            self.values.append(values)
            self.M+=values[0].T@((self.weights*L)[:,None]*values[0])
            offset=int(arm==1)
            for degree in range(self.order-offset):
                index=arm*self.order+offset+degree
                self.Kb[index,index]=1/((2*degree+1)*L**3)
            ends=self.field_matrices(arm,[0.,1.])
            self.outer.append(np.array([a[0] for a in ends[:4]]))
            self.ends.append(np.array([a[1] for a in ends[:4]]))
        self.jump=self.ends[0][1]-self.ends[1][1]
        self.D=self.outer[1][0]-self.ends[1][0]
        self.P=np.outer(self.D,self.D)
        self.K0=self.Kb+self.k_joint*np.outer(self.jump,self.jump)
        self.constraint=(self.ends[0][0]+self.ends[1][0])[None,:]
        self.Z=null_space(self.constraint)

    def field_matrices(self,arm,points):
        arrays=fields(self.order,points,self.lengths[arm],arm)
        out=[]
        for k,a in enumerate(arrays):
            v=np.zeros((len(points),2*self.order+1))
            v[:,arm*self.order:(arm+1)*self.order]=a
            if arm==1 and k==0:v[:,-1]=1
            out.append(v)
        return out

    def solve(self,tau,count=3):
        saturated=tau==RIGID
        if not saturated and (not np.isfinite(tau) or tau<0):raise ValueError('Invalid tau')
        C=np.vstack((self.constraint,self.D)) if saturated else self.constraint
        Z=null_space(C) if saturated else self.Z
        K=self.K0 if saturated else self.K0+tau*self.P
        Kr,Mr=Z.T@K@Z,Z.T@self.M@Z
        mass_chol=cholesky((Mr+Mr.T)/2,lower=True)
        chol=cholesky((Kr+Kr.T)/2,lower=True)
        a=solve_triangular(chol,Mr,lower=True)
        a=solve_triangular(chol,a.T,lower=True).T
        inv,v=eigh((a+a.T)/2,subset_by_index=[len(a)-count,len(a)-1])
        inv,v=inv[::-1],v[:,::-1];zeta=1/inv
        coef=Z@solve_triangular(chol.T,v,lower=False)/np.sqrt(inv)
        residual=Z.T@(K@coef-(self.M@coef)*zeta)
        scale=np.linalg.norm(Z.T@K@coef,axis=0)+np.linalg.norm(Z.T@self.M@coef,axis=0)*zeta
        diagnostics=dict(symmetry=float(max(np.max(abs(K-K.T)),np.max(abs(self.M-self.M.T)))),
            mass_positive=True,stiffness_positive=True,mass_cholesky_min=float(np.min(np.diag(mass_chol))),
            condition_mass=float(np.linalg.cond(Mr)),condition_stiffness=float(np.linalg.cond(Kr)),
            dimension=Z.shape[1],constraint_rank=int(np.linalg.matrix_rank(C)),
            matrix_residual=float(np.max(np.linalg.norm(residual,axis=0)/scale)),
            energy=float(np.max(abs(np.diag(coef.T@K@coef)/zeta-1))),
            mass_gram=float(np.max(abs(coef.T@self.M@coef-np.eye(count)))),
            essential=float(np.max(abs(C@coef))),gap12_relative=float(np.sqrt(zeta[1]/zeta[0])-1))
        per_mode=[]
        for j,z in enumerate(zeta):
            c=coef[:,j];e=np.array([a@c for a in self.ends]);o=np.array([a@c for a in self.outer])
            jump=self.jump@c;d=self.D@c
            moment=[o[1,2],e[0,2]+self.k_joint*jump,e[1,2]-self.k_joint*jump]
            # Exact saturation eliminates the multiplier using the exterior force.
            forces=[e[0,3]-e[1,3]+o[1,3]] if saturated else [o[1,3]+tau*d,e[0,3]-e[1,3]-tau*d]
            natural=np.array(moment+forces)
            bc_scale=max(1.,np.max(abs(e)),np.max(abs(o)),abs(tau*d) if not saturated else 0.)
            numerator=denominator=0.;arm_res=[]
            for arm,L in enumerate(self.lengths):
                W=self.values[arm][0]@c;D4=self.values[arm][4]@c
                ni=float(self.weights@((D4-z*W)**2)*L)
                di=float(self.weights@(D4**2+(z*W)**2)*L)
                numerator+=ni;denominator+=di;arm_res.append(math.sqrt(ni/di))
            per_mode.append(dict(natural_rows=natural.tolist(),boundary_scale=float(bc_scale),
                natural=float(np.max(abs(natural))/bc_scale),strong=math.sqrt(numerator/denominator),
                strong_arms=arm_res,strong_numerator=math.sqrt(numerator),strong_denominator=math.sqrt(denominator),
                D=float(d),reaction=float(-o[1,3]) if saturated else float(tau*d)))
        diagnostics.update(natural_boundary=max(r['natural'] for r in per_mode),strong_ode=max(r['strong'] for r in per_mode))
        return dict(zeta=zeta,coefficients=coef,diagnostics=diagnostics,per_mode=per_mode)


def boundary_matrix(zeta,tau,k_joint=.1,mu=.01):
    """Independent exact constant-coefficient EB ODE, eight boundary rows.

    Unknowns are the four initial state values on each arm; no polynomial
    Ritz coefficients or assembled energies enter this check.
    """
    A=np.zeros((4,4));A[0,1]=A[1,2]=A[2,3]=1;A[3,0]=zeta
    o=[np.eye(8)[:4],np.eye(8)[4:]]
    e=[expm(A*L)@a for L,a in zip((1-mu,1+mu),o)]
    jump=e[0][1]-e[1][1];D=o[1][0]-e[1][0]
    rows=[o[0][0],o[0][1],e[0][0]+e[1][0],o[1][2],e[0][2]+k_joint*jump,e[1][2]-k_joint*jump]
    if tau==RIGID:rows += [D,e[0][3]-e[1][3]+o[1][3]]
    else:rows += [o[1][3]+tau*D,e[0][3]-e[1][3]-tau*D]
    B=np.array(rows);return B/np.linalg.norm(B,axis=1)[:,None]


def independent_roots(approximate,tau,k_joint=.1,mu=.01):
    roots=[]
    for z in approximate:
        fun=lambda x:np.linalg.det(boundary_matrix(x,tau,k_joint,mu))
        lo,hi=z*(1-1e-4),z*(1+1e-4)
        if fun(lo)*fun(hi)>=0:raise ArithmeticError('Independent EB root not bracketed')
        root=brentq(fun,lo,hi,xtol=1e-13,rtol=2e-14)
        singular=np.linalg.svd(boundary_matrix(root,tau,k_joint,mu),compute_uv=False)
        roots.append(dict(zeta=float(root),relative_frequency_error=float(math.sqrt(root/z)-1),singular_ratio=float(singular[-1]/singular[0]),bracket=[lo,hi]))
    return roots
