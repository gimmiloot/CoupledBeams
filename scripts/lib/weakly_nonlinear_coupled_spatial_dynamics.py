"""Two endpoint-lifted arms, assembled from the unchanged seven-field action.

The seven shared coordinates are (d_X,d_Y,d_Z,alpha_X,alpha_Y,alpha_Z,c).
All transformations are constant variational pullbacks. No Cartesian moment
is substituted for the covector conjugate to alpha. No IVP is run here.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular, eigh
from scipy.spatial.transform import Rotation

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_spatial_dynamics import SpatialGalerkin, FIELDS
from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin
from scripts.lib.reddy_inplane_geometry import reddy_inplane_geometry, k
from scripts.lib import yartsev_ch2_monoclinic_rod as book

VERSION = "coupled-seven-field-endpoint-lifted-quartic-v1"


def rectangular_coefficients(b, h, length=1.):
    """Unchanged isotropic material and existing generalized torsion series."""
    shear = 1/(2*(1+.3))
    material = book.BookMaterial(E1_real=1., E2_real=1., G12_real=shear,
        G13_real=shear, G23_real=shear, nu12=.3, rho=1.,
        eta1=0., eta2=0., eta12=0., eta13=0., eta23=0.)
    point = book.make_rod_point(0., geometry=book.Geometry(a=b,b=h,length=length),
        material=material, material_mode="elastic")
    coefficient = rod.RodCoefficients.rectangular(1.,1.,.3,b,h,float(point.torsion.C_T.real))
    return coefficient, {"CT":coefficient.CT,"series_terms":point.torsion.terms_used,
                         "series_relative_tail":point.torsion.estimated_relative_tail}


class EndpointArm(SpatialGalerkin):
    """p-1 unchanged Shen functions plus s/L, with outer values zero only."""
    def __init__(self, coefficients, p, length, model):
        self.interior_degree = int(p)
        # Parent n=p is required for the extra lift. Its two extra Gauss points
        # retain its existing lower bound, and exactly integrate degree 4p.
        super().__init__(coefficients,p+1,length,nq=2*p+3,model=model,whiten=False)
        self.p = int(p)

    def _raw_basis(self, points, derivative):
        template = SimpleNamespace(p=self.interior_degree,n=self.interior_degree-1,length=self.length)
        interior = PlanarGalerkin._raw_basis(template,points,derivative)
        lift = (np.asarray(points)/self.length if derivative==0 else
                np.full(len(points),1/self.length) if derivative==1 else np.zeros(len(points)))
        return np.column_stack((interior,lift))


def force_tangents(arm, q, velocity, acceleration):
    """Local H(V)+g_q+M_q*a and g_v, using the assembled acceleration.

    An isolated-arm Jacobian would use the wrong endpoint acceleration. These
    are the same compiled action derivatives as SpatialGalerkin.jacobian.
    """
    values = arm._kinetic_variables(q,velocity)
    gq = arm._inertia_q.evaluate(values).reshape(7,4,arm.nq)
    gv = arm._inertia_v.evaluate(values).reshape(7,7,arm.nq)
    mq = arm._mass_q.evaluate(values[:4]).reshape(3,3,4,arm.nq)
    actual = arm.reconstruct(acceleration)
    for i,f in enumerate((3,4,5)):
        gq[f] += np.einsum("jkx,xj->kx",mq[i],actual[:,(3,4,5)],optimize=True)
    first = arm.potential(q,hessian=True)["hessian"].copy()
    second = np.zeros_like(first)
    for i,f in enumerate(FIELDS):
        left = arm.B[i]; ls=arm.slices[f]
        for j,fj in enumerate((3,4,5,6)):
            if np.any(gq[i,j]):first[ls,arm.slices[FIELDS[fj]]] += left.T@((arm.weights*gq[i,j])[:,None]*arm.B[fj])
        for j,fj in enumerate(FIELDS):
            if np.any(gv[i,j]):second[ls,arm.slices[fj]] += left.T@((arm.weights*gv[i,j])[:,None]*arm.B[j])
    return first,second


class CoupledSpatialGalerkin:
    """Full seven-field variational assembly; endpoint rotations use T7.T."""
    def __init__(self, coefficients, p, lengths=(.8,1.2), beta_deg=30., model=None):
        if model is None:raise ValueError("A verified frozen action is required")
        self.p=int(p);self.lengths=tuple(lengths);self.beta_deg=float(beta_deg)
        if self.p<2 or len(coefficients)!=2 or min(lengths)<=0:raise ValueError("Two positive arms and p>=2 required")
        self.arms=tuple(EndpointArm(c,p,L,model) for c,L in zip(coefficients,lengths))
        geometry=reddy_inplane_geometry(beta_deg)
        self.bases=tuple(np.column_stack((a.t,a.n,k)) for a in (geometry.arm1,geometry.arm2))
        self.T7=[]
        for B in self.bases:
            T=np.zeros((7,7));T[:3,:3]=B;T[3:6,3:6]=B@np.diag([1.,-1.,1.]);T[6,6]=1
            self.T7.append(T)
        n=p-1;self.ndof=14*n+7;self.joint_ids=np.arange(14*n,self.ndof)
        self.raw_maps=[];groups={"plane":[],"out":[],"rot":[]}
        for i in range(2):
            P=np.zeros((7*p,self.ndof))
            for f in range(7):
                ids=np.arange(i*7*n+f*n,i*7*n+(f+1)*n)
                P[f*p:f*p+n,ids]=np.eye(n)
                P[(f+1)*p-1,self.joint_ids]=self.T7[i].T[f]
                groups["out" if f in (2,3,4) else "plane"].extend(ids)
                if f in (3,4,5):groups["rot"].extend(ids)
            self.raw_maps.append(P)
        groups["plane"].extend(self.joint_ids[[0,1,5,6]])
        groups["out"].extend(self.joint_ids[[2,3,4]])
        groups["rot"].extend(self.joint_ids[[3,4,5]])
        self.plane=np.array(sorted(groups["plane"]));self.out=np.array(sorted(groups["out"]));self.rot=np.array(sorted(groups["rot"]))
        self.constant=np.setdiff1d(np.arange(self.ndof),self.rot)
        raw_mass=sum(P.T@a.M0@P for a,P in zip(self.arms,self.raw_maps))
        # Blockwise reversible whitening preserves the reflection eigenspaces.
        self.W=np.zeros_like(raw_mass)
        for ids in (self.plane,self.out):
            lower=np.linalg.cholesky(raw_mass[np.ix_(ids,ids)])
            self.W[np.ix_(ids,ids)]=solve_triangular(lower.T,np.eye(len(ids)),lower=False)
        self.raw_mass=raw_mass;self.maps=tuple(P@self.W for P in self.raw_maps)
        self.M0=self.W.T@raw_mass@self.W
        self.rot_maps=tuple(P[a._rot_ids][:,self.rot] for a,P in zip(self.arms,self.maps))
        self.constant_factor=cho_factor(self.M0[np.ix_(self.constant,self.constant)],lower=True)
        self._mass_key=self._mass_value=self._mass_factor=None
        self.reset_counters();self.K=self.potential(np.zeros(self.ndof),hessian=True)["hessian"]
        self.state_plane=np.r_[self.plane,self.plane+self.ndof]
        self.state_out=np.r_[self.out,self.out+self.ndof]

    def reset_counters(self):self.rhs_calls=self.jacobian_calls=self.mass_factorizations=0
    def counters(self):return {k:getattr(self,k) for k in ("rhs_calls","jacobian_calls","mass_factorizations")}
    def local(self,q):return tuple(P@q for P in self.maps)
    def raw(self,q):return self.W@q
    def from_raw(self,q):return np.linalg.solve(self.W,q)

    def potential(self,q,gradient=True,hessian=False):
        result={"V":0.}
        if gradient or hessian:result["gradient"]=np.zeros(self.ndof)
        if hessian:result["hessian"]=np.zeros((self.ndof,self.ndof))
        for a,P,qi in zip(self.arms,self.maps,self.local(q)):
            part=a.potential(qi,gradient=gradient,hessian=hessian);result["V"]+=part["V"]
            if gradient or hessian:result["gradient"]+=P.T@part["gradient"]
            if hessian:result["hessian"]+=P.T@part["hessian"]@P
        return result

    def rotational_mass(self,q):
        if self._mass_key is not None and np.array_equal(q,self._mass_key):return self._mass_value,self._mass_factor
        matrix=sum(P.T@a._rotational_mass(qi)[0]@P for a,P,qi in zip(self.arms,self.rot_maps,self.local(q)))
        factor=cho_factor(matrix,lower=True,check_finite=False)
        self._mass_key=q.copy();self._mass_value=matrix;self._mass_factor=factor;self.mass_factorizations+=1
        return matrix,factor

    def mass_matrix(self,q):
        M=self.M0.copy();M[np.ix_(self.rot,self.rot)]=self.rotational_mass(q)[0];return M

    def solve_mass(self,q,rhs):
        result=np.empty_like(rhs)
        result[self.constant]=cho_solve(self.constant_factor,rhs[self.constant],check_finite=False)
        result[self.rot]=cho_solve(self.rotational_mass(q)[1],rhs[self.rot],check_finite=False)
        return result

    def inertial_terms(self,q,v):return sum(P.T@a.inertial_terms(qi,vi) for a,P,qi,vi in zip(self.arms,self.maps,self.local(q),self.local(v)))
    def kinetic(self,q,v):return sum(a.kinetic(qi,vi) for a,qi,vi in zip(self.arms,self.local(q),self.local(v)))
    def energy(self,q,v):return self.kinetic(q,v)+self.potential(q,gradient=False)["V"]
    def acceleration(self,q,v):return self.solve_mass(q,-self.potential(q)["gradient"]-self.inertial_terms(q,v))

    def rhs(self,t,state):
        del t;self.rhs_calls+=1
        if state.shape!=(2*self.ndof,) or not np.isfinite(state).all():raise FloatingPointError("Finite full coupled state required")
        return np.r_[state[self.ndof:],self.acceleration(state[:self.ndof],state[self.ndof:])]

    def jacobian(self,t,state):
        del t;self.jacobian_calls+=1
        q,v=state[:self.ndof],state[self.ndof:];acc=self.acceleration(q,v)
        first=np.zeros((self.ndof,self.ndof));second=np.zeros_like(first)
        for a,P,qi,vi,ai in zip(self.arms,self.maps,self.local(q),self.local(v),self.local(acc)):
            A,B=force_tangents(a,qi,vi,ai);first+=P.T@A@P;second+=P.T@B@P
        J=np.zeros((2*self.ndof,2*self.ndof));J[:self.ndof,self.ndof:]=np.eye(self.ndof)
        J[self.ndof:,:self.ndof]=self.solve_mass(q,-first)
        J[self.ndof:,self.ndof:]=self.solve_mass(q,-second)
        return J

    def embed_plane(self,x):
        result=np.zeros(2*self.ndof);result[self.state_plane]=x;return result
    def plane_rhs(self,t,x):return self.rhs(t,self.embed_plane(x))[self.state_plane]
    def plane_jacobian(self,t,x):return self.jacobian(t,self.embed_plane(x))[np.ix_(self.state_plane,self.state_plane)]

    def joint_diagnostics(self,q,v):
        acceleration=self.acceleration(q,v);forces=[];powers=[];points=[];orientations=[];values=[]
        for a,P,T,B,qi,vi,ai in zip(self.arms,self.maps,self.T7,self.bases,self.local(q),self.local(v),self.local(acceleration)):
            residual=a.mass_matrix(qi)@ai+a.inertial_terms(qi,vi)+a.potential(qi)["gradient"]
            endpoint=residual[np.arange(a.p-1,a.ndof,a.p)]
            forces.append(T@endpoint);node=a.reconstruct(qi,np.array([a.length]))[0];speed=a.reconstruct(vi,np.array([a.length]))[0]
            points.append(B@node[:3]);orientations.append(B@Rotation.from_rotvec(node[3:6]*[1,-1,1]).as_matrix()@B.T)
            values.append(node[6]);powers.append(float(endpoint@speed))
        return {"translation_difference":float(np.max(abs(points[0]-points[1]))),
            "orientation_difference":float(np.max(abs(orientations[0]-orientations[1]))),"contraction_difference":float(abs(values[0]-values[1])),
            "generalized_nodal_balance":np.sum(forces,axis=0),"joint_power":sum(powers),
            "effort_qualification":"discrete action covectors, not raw Cartesian physical moments"}

    def energy_rate(self,q,v):
        acc=self.acceleration(q,v);rate=0.
        for a,qi,vi,ai in zip(self.arms,self.local(q),self.local(v),self.local(acc)):
            rate+=a.energy_rate(qi,vi,ai)
        return rate

    def body_force(self,global_acceleration):
        local=[]
        for a,B in zip(self.arms,self.bases):
            components=a.coefficients.m*(B.T@global_acceleration);force=np.zeros(a.ndof)
            for j,f in enumerate(FIELDS[:3]):force[a.slices[f]]=components[j]*(a.B[j].T@a.weights)
            local.append(force)
        return sum(P.T@f for P,f in zip(self.maps,local))

    def safety(self,q):
        from scripts.lib.nlsp_spatial_1d_program import safety_check
        return [safety_check(a,qi) for a,qi in zip(self.arms,self.local(q))]

    def eigenpairs(self,plane=True):
        ids=self.plane if plane else self.out
        values,vectors=eigh(self.K[np.ix_(ids,ids)],self.M0[np.ix_(ids,ids)],check_finite=False)
        if values.min()<=0:raise ArithmeticError("Nonpositive constrained reference stiffness")
        return np.sqrt(values),vectors

    def reference_energy_parts(self,position,velocity):
        result=np.zeros((2,2))
        for i,(a,q,v) in enumerate(zip(self.arms,self.local(position),self.local(velocity))):
            f=a.reconstruct(q);s=a.reconstruct(q,derivative=1);speed=a.reconstruct(v);c=a.coefficients
            result[i,0]=float(a.weights@(c.S*(s[:,2]-f[:,4])**2+c.Bb*s[:,4]**2+c.m*speed[:,2]**2+c.jb*speed[:,4]**2))
            result[i,1]=float(a.weights@(c.CT*s[:,3]**2+(c.jp+c.jb)*speed[:,3]**2))
        return result
