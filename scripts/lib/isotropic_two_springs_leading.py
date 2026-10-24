"""Leading beta/epsilon problem; not the full RL_FULL/Timoshenko solver.

Canonical weak form: docs/theory/isotropic_two_springs_angle_asymptotics.md
(14)--(21). Eigenvalue zeta=Lambda**4; omega_lim=sqrt(zeta).
P is y*y, with no factor 1/2. K(xi)=K0+xi**2/(L1+L2)*P.
Use the integrated-Legendre conditioning pattern of the existing Ritz checks.
The polynomial spaces equal t**2*P_n on arm 1 and t*P_n on arm 2,
but orthogonal curvatures avoid high-derivative cancellation. Natural
force/moment conditions are not constrained.
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import Legendre
from scipy.linalg import cholesky, eigh, null_space, solve_triangular


def fields(order, points, length, arm):
    """Values and first four s derivatives; endpoint positions remain free."""
    t = np.asarray(points, dtype=float)
    derivatives = [[] for _ in range(5)]
    if arm == 1:
        for k in range(5):derivatives[k].append(t if k==0 else np.ones_like(t) if k==1 else np.zeros_like(t))
    for n in range(order-int(arm==1)):
        poly=Legendre.basis(n);p1=poly.integ();p2=p1.integ();x=2*t-1
        values=[(p2(x)-p2(-1))/4-t*p1(-1)/2,
                (p1(x)-p1(-1))/2,poly(x),2*poly.deriv()(x),4*poly.deriv(2)(x)]
        for k in range(5):derivatives[k].append(values[k])
    return [np.array(a).T/length**k for k,a in enumerate(derivatives)]


@dataclass
class Ritz:
    k_joint: float
    k_end: float
    mu: float = .01
    order: int = 12
    quadrature: int = 64

    def __post_init__(self):
        if min(self.k_joint,self.k_end) <= 0 or not abs(self.mu)<1 or self.order<4:
            raise ValueError('Leading finite-positive-spring contract')
        if self.quadrature < self.order+3:
            raise ValueError('Quadrature must integrate the polynomial mass exactly')
        self.lengths = (1-self.mu,1+self.mu)
        n = 2*self.order
        self.M = np.zeros((n,n)); self.Kb = np.zeros((n,n))
        self.ends = []; self.outer = []
        nodes,weights = np.polynomial.legendre.leggauss(self.quadrature)
        self.nodes = (nodes+1)/2; self.weights = weights/2
        for i,L in enumerate(self.lengths):
            sl = slice(i*self.order,(i+1)*self.order)
            f = fields(self.order,self.nodes,L,i)
            self.M[sl,sl] = f[0].T@((self.weights*L)[:,None]*f[0])
            self.Kb[sl,sl] = f[2].T@((self.weights*L)[:,None]*f[2])
            e = [np.zeros(n) for _ in range(4)]
            o = [np.zeros(n) for _ in range(4)]
            fe = fields(self.order,[0.,1.],L,i)
            for k in range(4):e[k][sl],o[k][sl] = fe[k][1],fe[k][0]
            self.ends.append(e);self.outer.append(o)
        self.jump = self.ends[0][1]-self.ends[1][1]
        self.end_rotation = self.outer[1][1]
        self.y = self.ends[1][0]
        self.P = np.outer(self.y,self.y)
        self.Kj = self.k_joint*np.outer(self.jump,self.jump)
        self.Ke = self.k_end*np.outer(self.end_rotation,self.end_rotation)
        self.K0 = self.Kb+self.Kj+self.Ke
        self.constraint = (self.ends[0][0]+self.ends[1][0])[None,:]
        self.Z = null_space(self.constraint)
        self.matrix_checks = dict(
            symmetry=max(float(np.max(abs(a-a.T))) for a in (self.K0,self.M,self.P)),
            constraint_rank=int(np.linalg.matrix_rank(self.constraint)),
            P_rank=int(np.linalg.matrix_rank(self.P)),
            dimension=n-1, mass_cholesky=True)
        cholesky(self.Z.T@self.M@self.Z,lower=True)

    def solve(self, xi, count=3, node_fixed=False):
        """Lowest count Ritz positions via stiffness-whitened symmetric solve.

        The exact node-fixed limit uses another constraint, never a penalty.
        Low roots are the largest eigenvalues of K**(-1/2) M K**(-1/2).
        """
        if not np.isfinite(xi) or xi<0:raise ValueError('xi must be finite and nonnegative')
        constraint = np.vstack((self.constraint,self.y)) if node_fixed else self.constraint
        Z = null_space(constraint) if node_fixed else self.Z
        t = 0. if node_fixed else xi*xi/sum(self.lengths)
        K = self.K0+t*self.P
        Kr,Mr = Z.T@K@Z,Z.T@self.M@Z
        L = cholesky((Kr+Kr.T)/2,lower=True)
        a = solve_triangular(L,Mr,lower=True)
        a = solve_triangular(L,a.T,lower=True).T
        inverse,v = eigh((a+a.T)/2,subset_by_index=[len(a)-count,len(a)-1])
        inverse,v = inverse[::-1],v[:,::-1]
        zeta = 1/inverse
        coef = Z@solve_triangular(L.T,v,lower=False)/np.sqrt(inverse)
        q = Z.T@(K@coef-(self.M@coef)*zeta)
        scale = np.linalg.norm(Z.T@K@coef,axis=0)+np.linalg.norm(Z.T@self.M@coef,axis=0)*zeta
        y = self.y@coef
        diagnostic = dict(**self.matrix_checks,
            eigen_residual=float(np.max(np.linalg.norm(q,axis=0)/scale)),
            mass_gram=float(np.max(abs(coef.T@self.M@coef-np.eye(count)))),
            energy=float(np.max(abs(np.diag(coef.T@K@coef)/zeta-1))),
            constraint=float(np.max(abs(constraint@coef))),
            gap12_relative=float(np.sqrt(zeta[1]/zeta[0])-1))
        diagnostic.update(constraint_rank=int(np.linalg.matrix_rank(constraint)), dimension=Z.shape[1])
        natural=[];ode=[];essential=[]
        for j in range(count):
            c=coef[:,j];end=[[e@c for e in arm] for arm in self.ends]
            outer=[[e@c for e in arm] for arm in self.outer]
            residual=np.array([-outer[1][2]+self.k_end*outer[1][1],
                end[0][2]+self.k_joint*(self.jump@c),
                end[1][2]-self.k_joint*(self.jump@c),
                end[0][3]-end[1][3]+t*(self.y@c)])
            if node_fixed:residual=residual[:-1]
            norm=max(1.,np.max(abs(np.array(end))),np.max(abs(np.array(outer))),abs(t*(self.y@c)))
            natural.append(float(np.max(abs(residual))/norm))
            essential.append(max(abs(outer[0][0]),abs(outer[0][1]),abs(outer[1][0]),abs(end[0][0]+end[1][0])))
            num=den=0.
            for i,Li in enumerate(self.lengths):
                f=fields(self.order,self.nodes,Li,i);ci=c[i*self.order:(i+1)*self.order]
                W,D4=f[0]@ci,f[4]@ci
                num+=self.weights@((D4-zeta[j]*W)**2)*Li
                den+=self.weights@(D4**2+(zeta[j]*W)**2)*Li
            ode.append(float(np.sqrt(num/den)))
        diagnostic.update(natural_boundary=max(natural),strong_ode=max(ode),essential=max(essential))
        return dict(zeta=zeta,omega_lim=np.sqrt(zeta),Lambda=zeta**.25,coefficients=coef,
            node_value=y,s_xi=xi*y*y/(sum(self.lengths)*zeta),
            C_at_zero=y[0]**2/(2*sum(self.lengths)*zeta[0]) if xi==0 and not node_fixed else None,
            diagnostics=diagnostic)
