"""Variational spatial discretization of the audited planar quartic action.

Four independent fields (u,w,theta,c) use the same Shen Legendre space with
only essential endpoint values fixed.  The potential and its derivatives are
compiled from the accepted seven-field Polynomial artifact, restricted to the
invariant plane.  No trigonometric full-model PDE is substituted for the cubic
action.  All polynomial/CAS work finishes before any RHS evaluation.

The coefficient coordinates may be whitened by the constant resting mass.
This is a reversible basis change, not modal reduction.  The nonlinear theta
mass is assembled and solved at each state.  Its derivative is retained in the
analytic RHS Jacobian; neither a constant-mass substitution nor an inverse
Taylor series is used.  Time integration and cache orchestration belong to
the single diagnostic CLI, not this helper.
"""
from __future__ import annotations

import math
from typing import Callable

import numpy as np
from numpy.polynomial.legendre import Legendre, leggauss, legvander
from scipy.linalg import cho_factor, cho_solve, eigh, solve_triangular

from scripts.lib import weakly_nonlinear_spatial_rod as rod


VERSION = "audited-planar-quartic-shen-variable-mass-v2-lazy-potential"
FIELDS = ("u", "w", "theta", "c")
PLANAR_SPATIAL_INDICES = (0, 1, 5, 6)
_LOCAL_POTENTIAL_NAMES = ("u_s", "w_s", "theta", "c", "c_s", "theta_s")
_LOCAL_FIELDS = (0, 1, 2, 3, 3, 2)
_JET_SUFFIXES = ("", "_s", "_t", "_ss", "_st", "_tt")


def _restrict_to_plane(polynomial):
    return polynomial.substitute({field+suffix: 0 for field in ("v", "Phi", "psi")
                                  for suffix in _JET_SUFFIXES})


class _CompiledPolynomials:
    """Small numeric coefficient matrix over shared monomial features.

    Only symbols explicitly listed in names and constant model coefficients
    are admitted.  Unknown active variables raise rather than disappear.
    Numeric derivatives are produced by exact polynomial differentiation at
    construction, not by finite differences in a dynamics call.
    """

    def __init__(self, polynomials, names, coefficients):
        self.names = tuple(names)
        positions = {name: index for index, name in enumerate(self.names)}
        values = coefficients.values()
        rows, all_exponents = [], set()
        for polynomial in polynomials:
            terms = {}
            for monomial, rational in polynomial.terms.items():
                exponent = [0]*len(self.names)
                coefficient = float(rational)
                for index in monomial:
                    name = rod.SYMBOL_ORDER[index]
                    if name in positions:
                        exponent[positions[name]] += 1
                    elif name in values:
                        coefficient *= values[name]
                    else:
                        raise ValueError(f"Unexpected active polynomial variable {name}")
                key = tuple(exponent)
                terms[key] = terms.get(key, 0.)+coefficient
            rows.append(terms)
            all_exponents.update(terms)
        keys = sorted(all_exponents)
        self.exponents = np.array(keys, dtype=np.int8)
        self.coefficients = np.array([[row.get(key, 0.) for key in keys] for row in rows])
        self.max_power = int(self.exponents.max(initial=0))
        # The exponent pattern is fixed by the audited polynomials.  Finding
        # these selectors on every RHS call does no state-dependent work.
        self._multipliers = tuple(
            (index, power, selected)
            for index in range(len(self.names))
            for power in range(1, self.max_power+1)
            if (selected := np.flatnonzero(self.exponents[:, index] == power)).size)

    def evaluate(self, variables):
        variables = np.asarray(variables)
        if variables.ndim != 2 or variables.shape[0] != len(self.names):
            raise ValueError("Local variable array has the wrong shape")
        features = np.ones((len(self.exponents), variables.shape[1]), dtype=variables.dtype)
        powers = [np.ones_like(variables)]
        for _ in range(self.max_power):
            powers.append(powers[-1]*variables)
        for index, power, selected in self._multipliers:
            features[selected] *= powers[power][index]
        return self.coefficients@features


class PlanarGalerkin:
    """Four-field Shen discretization of the accepted cubic coordinate model.

    p is the maximum polynomial degree; each field has p-1 independent
    coefficients and ndof=4(p-1).  Gauss quadrature integrates through degree
    2*nq-1.  The quartic action's maximum spatial degree is <=4*p, so nq>=2*p+1
    is sufficient.  No endpoint slope constraint or nonlinear product filter
    is imposed.
    """

    fields = FIELDS

    def __init__(self, coefficients: rod.RodCoefficients, p: int, length=1.,
                 nq=None, model=None, whiten=True):
        if not isinstance(p, (int, np.integer)) or p < 2:
            raise ValueError("Maximum polynomial degree p must be an integer >=2")
        if not math.isfinite(length) or length <= 0:
            raise ValueError("Rod length must be finite and positive")
        nq = 2*p+1 if nq is None else int(nq)
        if nq < 2*p+1:
            raise ValueError("Quartic action requires nq>=2*p+1; reduced integration is not allowed")
        self.coefficients, self.p, self.length, self.nq = coefficients, int(p), float(length), nq
        self.n, self.ndof = self.p-1, 4*(self.p-1)
        self.whiten = bool(whiten)
        self.slices = {field: slice(i*self.n, (i+1)*self.n) for i, field in enumerate(FIELDS)}
        self._indices = {field: np.arange(s.start, s.stop) for field, s in self.slices.items()}
        xi, weights = leggauss(self.nq)
        self.x, self.weights = (xi+1)*self.length/2, weights*self.length/2
        raw = self._raw_basis(self.x, 0)
        gram = raw.T@(self.weights[:, None]*raw)
        base_masses = (coefficients.m, coefficients.m, coefficients.jp, coefficients.jp)
        self._transforms, self._constant_masses, self._constant_factors = [], [], []
        for mass in base_masses:
            resting = mass*gram
            lower = np.linalg.cholesky(resting)
            transform = solve_triangular(lower.T, np.eye(self.n), lower=False) if self.whiten else np.eye(self.n)
            constant = transform.T@resting@transform
            self._transforms.append(transform)
            self._constant_masses.append(constant)
            self._constant_factors.append(cho_factor(constant, lower=True, check_finite=False))
        self.B = tuple(raw@transform for transform in self._transforms)
        self.D = tuple(self._raw_basis(self.x, 1)@transform for transform in self._transforms)
        self.D2 = tuple(self._raw_basis(self.x, 2)@transform for transform in self._transforms)
        self.M0 = np.zeros((self.ndof, self.ndof))
        for field, mass in zip(FIELDS, self._constant_masses):
            self.M0[self.slices[field], self.slices[field]] = mass
        self.model = rod.derive_polynomials() if model is None else model
        potential = _restrict_to_plane(self.model.V4)
        kinetic = _restrict_to_plane(self.model.T4)
        symbols = self.model.symbols
        expected = (symbols["m"]*(symbols["u_t"]**2+symbols["w_t"]**2)+
                    symbols["jp"]*symbols["c_t"]**2+
                    symbols["jp"]*(1+symbols["c"])**2*symbols["theta_t"]**2)/2
        if kinetic != expected:
            raise ValueError("Audited planar kinetic restriction differs from the expected quartic action")
        gradients = [potential.derivative(name) for name in _LOCAL_POTENTIAL_NAMES]
        hessians = [gradient.derivative(name) for gradient in gradients for name in _LOCAL_POTENTIAL_NAMES]
        # Ordinary RHS requests need V/gradient, not the 36 local Hessian
        # entries.  Compile separate paths from the same exact derivatives.
        self._potential_energy = _CompiledPolynomials([potential],
                                                      _LOCAL_POTENTIAL_NAMES, coefficients)
        self._potential_gradient = _CompiledPolynomials([potential]+gradients,
                                                        _LOCAL_POTENTIAL_NAMES, coefficients)
        self._potential_hessian = _CompiledPolynomials(hessians,
                                                       _LOCAL_POTENTIAL_NAMES, coefficients)
        names = tuple(field+suffix for suffix in _JET_SUFFIXES for field in FIELDS)
        residuals = [_restrict_to_plane(self.model.residual_a[index]) for index in PLANAR_SPATIAL_INDICES]
        self._residual = _CompiledPolynomials(residuals, names, coefficients)
        self._potential_matrices = (self.D[0], self.D[1], self.B[2], self.B[3], self.D[3], self.D[2])
        self._potential_cache_q = None
        self._potential_cache = None
        self._theta_mass_cache_c = None
        self._theta_mass_cache_matrix = None
        self._theta_mass_cache_factor = None
        self._linear_modes = {}
        self.reset_counters()
        self.K = self.potential(np.zeros(self.ndof), hessian=True)["hessian"]

    @property
    def linear_stiffness(self):
        return self.K

    @property
    def resting_mass(self):
        return self.M0

    def reset_counters(self):
        self.rhs_calls = 0
        self.jacobian_calls = 0
        self.mass_factorizations = 0
        self.force_evaluations = 0
        self.linear_eigendecompositions = 0

    def counters(self):
        return {name: getattr(self, name) for name in
                ("rhs_calls", "jacobian_calls", "mass_factorizations", "force_evaluations", "linear_eigendecompositions")}

    def _raw_basis(self, points, derivative):
        points = np.asarray(points, dtype=float)
        if points.ndim != 1 or np.any(points < 0) or np.any(points > self.length):
            raise ValueError("Physical basis points must lie in [0,L]")
        xi = 2*points/self.length-1
        if derivative == 0:
            values = legvander(xi, self.p)
            return values[:, :self.n]-values[:, 2:self.n+2]
        if derivative not in (1, 2):
            raise ValueError("Supported basis derivative orders are 0,1,2")
        return np.column_stack([(Legendre.basis(n)-Legendre.basis(n+2)).deriv(derivative)(xi)
                                for n in range(self.n)])*(2/self.length)**derivative

    def basis_at(self, points, derivative=0):
        raw = self._raw_basis(points, derivative)
        return {field: raw@transform for field, transform in zip(FIELDS, self._transforms)}

    def _check_coordinate(self, value):
        value = np.asarray(value)
        if value.shape != (self.ndof,):
            raise ValueError(f"Expected {self.ndof} coefficient coordinates")
        return value

    def raw_coefficients(self,coordinate):
        """Physical Shen coefficients before the reversible resting-mass scaling.

        Flat ordering is the same four independent field blocks as coordinate.
        This representation is suitable for saved cross-resolution snapshots;
        it is not a change to the physical fields or a projection onto modes.
        """
        coordinate=self._check_coordinate(coordinate)
        return np.concatenate([transform@coordinate[self.slices[field]]
                               for field,transform in zip(FIELDS,self._transforms)])

    def from_raw_coefficients(self,raw):
        """Inverse of raw_coefficients using solves, with no mass approximation."""
        raw=self._check_coordinate(raw)
        return np.concatenate([solve_triangular(transform,raw[self.slices[field]],lower=False,
                                               check_finite=False)
                               for field,transform in zip(FIELDS,self._transforms)])

    def reconstruct(self, coordinate, points=None, derivative=0):
        coordinate = self._check_coordinate(coordinate)
        if points is None:
            matrices = self.B if derivative == 0 else self.D if derivative == 1 else self.D2 if derivative == 2 else None
            if matrices is None:
                raise ValueError("Supported derivative orders are 0,1,2")
        else:
            matrices = tuple(self.basis_at(points, derivative).values())
        return np.column_stack([matrix@coordinate[self.slices[field]] for field, matrix in zip(FIELDS, matrices)])

    def reconstruct_series(self, coordinate_rows, points=None, derivative=0):
        coordinate_rows = np.asarray(coordinate_rows)
        if coordinate_rows.ndim != 2 or coordinate_rows.shape[1] != self.ndof:
            raise ValueError("Coefficient series must have shape (nt,ndof)")
        matrices = self.basis_at(self.x if points is None else points, derivative)
        return np.stack([coordinate_rows[:, self.slices[field]]@matrices[field].T for field in FIELDS], axis=2)

    def project(self, fields: np.ndarray | Callable):
        values = np.asarray(fields(self.x) if callable(fields) else fields)
        if values.shape != (self.nq, 4):
            raise ValueError("Physical projection fields must have shape (nq,4)")
        coordinate = np.empty(self.ndof)
        for i, field in enumerate(FIELDS):
            matrix = self.B[i]
            gram = matrix.T@(self.weights[:, None]*matrix)
            coordinate[self.slices[field]] = cho_solve(cho_factor(gram, lower=True, check_finite=False),
                    matrix.T@(self.weights*values[:, i]), check_finite=False)
        return coordinate

    def _local_variables(self, coordinate):
        coordinate = self._check_coordinate(coordinate)
        return np.vstack([matrix@coordinate[self.slices[FIELDS[field]]]
                          for field, matrix in zip(_LOCAL_FIELDS, self._potential_matrices)])

    def _local_potential(self, coordinate, gradient=True, hessian=False):
        coordinate = self._check_coordinate(coordinate)
        if self._potential_cache_q is None or not np.array_equal(coordinate, self._potential_cache_q):
            self.force_evaluations += 1
            self._potential_cache_q = coordinate.copy()
            self._potential_cache = {"variables": self._local_variables(coordinate),
                                     "energy": None, "gradient": None, "hessian": None}
        cache = self._potential_cache
        if gradient and cache["gradient"] is None:
            values = self._potential_gradient.evaluate(cache["variables"])
            cache["energy"], cache["gradient"] = float(self.weights@values[0]), values[1:]
        elif cache["energy"] is None:
            values = self._potential_energy.evaluate(cache["variables"])
            cache["energy"] = float(self.weights@values[0])
        if hessian and cache["hessian"] is None:
            values = self._potential_hessian.evaluate(cache["variables"])
            cache["hessian"] = values.reshape(6, 6, self.nq)
        return cache["energy"], cache["gradient"], cache["hessian"]

    def potential(self, coordinate, gradient=True, hessian=False):
        energy, local_gradient, local_hessian = self._local_potential(coordinate, gradient, hessian)
        result = {"V": energy}
        if gradient:
            force = np.zeros(self.ndof)
            for field, matrix, values in zip(_LOCAL_FIELDS, self._potential_matrices, local_gradient):
                force[self.slices[FIELDS[field]]] += matrix.T@(self.weights*values)
            result["gradient"] = force
        if hessian:
            tangent = np.zeros((self.ndof,self.ndof))
            for k in range(6):
                left = self._potential_matrices[k]
                ls = self.slices[FIELDS[_LOCAL_FIELDS[k]]]
                for ell in range(k, 6):
                    right = self._potential_matrices[ell]
                    rs = self.slices[FIELDS[_LOCAL_FIELDS[ell]]]
                    block = left.T@((self.weights*local_hessian[k,ell])[:,None]*right)
                    tangent[ls, rs] += block
                    if k != ell:
                        tangent[rs, ls] += block.T
            result["hessian"] = tangent
        return result

    def _theta_mass(self, coordinate):
        cc = coordinate[self.slices["c"]]
        if self._theta_mass_cache_c is not None and np.array_equal(cc,self._theta_mass_cache_c):
            return self._theta_mass_cache_matrix, self._theta_mass_cache_factor
        c = self.B[3]@cc
        if not np.all(np.isfinite(c)) or np.min(1+c) <= 0:
            raise FloatingPointError("Contraction scale 1+c must remain finite and positive")
        matrix = self.coefficients.jp*self.B[2].T@((self.weights*(1+c)**2)[:,None]*self.B[2])
        factor = cho_factor(matrix, lower=True, check_finite=False)
        self.mass_factorizations += 1
        self._theta_mass_cache_c = cc.copy()
        self._theta_mass_cache_matrix, self._theta_mass_cache_factor = matrix, factor
        return matrix, factor

    def mass_matrix(self, coordinate):
        coordinate = self._check_coordinate(coordinate)
        matrix = self.M0.copy()
        matrix[self.slices["theta"],self.slices["theta"]] = self._theta_mass(coordinate)[0]
        return matrix

    def _solve_mass(self, coordinate, rhs):
        result = np.empty_like(rhs)
        for i, field in enumerate(FIELDS):
            factor = self._theta_mass(coordinate)[1] if field == "theta" else self._constant_factors[i]
            result[self.slices[field]] = cho_solve(factor,rhs[self.slices[field]],check_finite=False)
        return result

    def inertial_terms(self, coordinate, velocity):
        coordinate, velocity = self._check_coordinate(coordinate),self._check_coordinate(velocity)
        c = self.B[3]@coordinate[self.slices["c"]]
        ct = self.B[3]@velocity[self.slices["c"]]
        thetat = self.B[2]@velocity[self.slices["theta"]]
        result = np.zeros(self.ndof)
        jp = self.coefficients.jp
        result[self.slices["theta"]] = 2*jp*self.B[2].T@(self.weights*(1+c)*ct*thetat)
        result[self.slices["c"]] = -jp*self.B[3].T@(self.weights*(1+c)*thetat**2)
        return result

    def acceleration(self, coordinate, velocity, linear=False):
        coordinate,velocity = self._check_coordinate(coordinate),self._check_coordinate(velocity)
        if linear:
            result = np.empty(self.ndof)
            force = -self.K@coordinate
            for i,field in enumerate(FIELDS):
                result[self.slices[field]]=cho_solve(self._constant_factors[i],force[self.slices[field]],check_finite=False)
            return result
        force = self.potential(coordinate)["gradient"]+self.inertial_terms(coordinate,velocity)
        return self._solve_mass(coordinate,-force)

    def rhs(self, time, state):
        del time
        self.rhs_calls += 1
        state=np.asarray(state)
        if state.shape != (2*self.ndof,) or not np.all(np.isfinite(state)):
            raise FloatingPointError("Dynamics state must contain finite coordinates and velocities")
        coordinate,velocity=state[:self.ndof],state[self.ndof:]
        return np.concatenate((velocity,self.acceleration(coordinate,velocity)))

    def linear_rhs(self,time,state):
        """Constant resting-action control, independent of nonlinear RHS work."""
        del time
        state=np.asarray(state)
        if state.shape != (2*self.ndof,):
            raise ValueError("Linear state has the wrong size")
        return np.concatenate((state[self.ndof:],self.acceleration(state[:self.ndof],state[self.ndof:],linear=True)))

    def linear_jacobian(self,time=None,state=None):
        del time,state
        result=np.zeros((2*self.ndof,2*self.ndof))
        result[:self.ndof,self.ndof:]=np.eye(self.ndof)
        for i,field in enumerate(FIELDS):
            rows=self.slices[field]
            result[self.ndof+rows.start:self.ndof+rows.stop,:self.ndof]=cho_solve(
                self._constant_factors[i],-self.K[rows],check_finite=False)
        return result

    def jacobian(self, time, state):
        """Analytic RHS Jacobian, including dM/dc multiplied by acceleration."""
        del time
        self.jacobian_calls += 1
        state=np.asarray(state)
        coordinate,velocity=state[:self.ndof],state[self.ndof:]
        potential=self.potential(coordinate,hessian=True)
        acceleration=self._solve_mass(coordinate,-potential["gradient"]-self.inertial_terms(coordinate,velocity))
        c=self.B[3]@coordinate[self.slices["c"]]
        ct=self.B[3]@velocity[self.slices["c"]]
        thetat=self.B[2]@velocity[self.slices["theta"]]
        thetaacc=self.B[2]@acceleration[self.slices["theta"]]
        t,csl=self.slices["theta"],self.slices["c"]
        Bt,Bc,W,jp=self.B[2],self.B[3],self.weights,self.coefficients.jp
        tangent=potential["hessian"].copy()
        tangent[t,csl] += 2*jp*Bt.T@((W*(ct*thetat+(1+c)*thetaacc))[:,None]*Bc)
        tangent[csl,csl] += -jp*Bc.T@((W*thetat**2)[:,None]*Bc)
        damping=np.zeros((self.ndof,self.ndof))
        damping[t,t] = 2*jp*Bt.T@((W*(1+c)*ct)[:,None]*Bt)
        damping[t,csl] = 2*jp*Bt.T@((W*(1+c)*thetat)[:,None]*Bc)
        damping[csl,t] = -2*jp*Bc.T@((W*(1+c)*thetat)[:,None]*Bt)
        result=np.zeros((2*self.ndof,2*self.ndof))
        result[:self.ndof,self.ndof:]=np.eye(self.ndof)
        result[self.ndof:,:self.ndof]=self._solve_mass(coordinate,-tangent)
        result[self.ndof:,self.ndof:]=self._solve_mass(coordinate,-damping)
        return result

    def energy(self, coordinate, velocity):
        fields=self.reconstruct(coordinate)
        speeds=self.reconstruct(velocity)
        p=self.coefficients
        kinetic=(p.m*(speeds[:,0]**2+speeds[:,1]**2)+p.jp*speeds[:,3]**2+
                 p.jp*(1+fields[:,3])**2*speeds[:,2]**2)/2
        return float(self.weights@kinetic+self.potential(coordinate,gradient=False)["V"])

    def energy_rate(self,coordinate,velocity,acceleration=None):
        velocity=self._check_coordinate(velocity)
        acceleration=self.acceleration(coordinate,velocity) if acceleration is None else self._check_coordinate(acceleration)
        fields=self.reconstruct(coordinate)
        speeds=self.reconstruct(velocity)
        mass_rate_power=self.coefficients.jp*float(self.weights@((1+fields[:,3])*speeds[:,3]*speeds[:,2]**2))
        return float(velocity@(self.mass_matrix(coordinate)@acceleration+self.potential(coordinate)["gradient"])+mass_rate_power)

    def weak_residual(self,coordinate,velocity,acceleration):
        """Independent test-function projection of audited continuum residuals."""
        local=[]
        for values,derivative in ((coordinate,0),(coordinate,1),(velocity,0),
                                  (coordinate,2),(velocity,1),(acceleration,0)):
            local.extend(self.reconstruct(values,derivative=derivative).T)
        residual=self._residual.evaluate(np.asarray(local))
        return np.concatenate([matrix.T@(self.weights*values) for matrix,values in zip(self.B,residual)])

    def linear_eigenpairs(self,block=None):
        if block not in (None,"mh","timoshenko"):
            raise ValueError("Linear block must be mh,timoshenko or None")
        if block not in self._linear_modes:
            names=FIELDS if block is None else ("u","c") if block=="mh" else ("w","theta")
            ids=np.concatenate([self._indices[name] for name in names])
            values,vectors=eigh(self.K[np.ix_(ids,ids)],self.M0[np.ix_(ids,ids)],check_finite=False)
            if np.any(values<=0) or not np.all(np.isfinite(values)):
                raise ArithmeticError("Fixed-fixed linear action has a nonpositive eigenvalue")
            expanded=np.zeros((self.ndof,len(values)))
            expanded[ids]=vectors
            omega=np.sqrt(values)
            self._linear_modes[block]={"eigenvalues":values,"omega":omega,"frequency_hz":omega/(2*np.pi),"vectors":expanded}
            self.linear_eigendecompositions+=1
        return self._linear_modes[block]

    def linear_reference(self,coordinate0,velocity0,times):
        """Exact-in-time full semidiscrete linear solution of the same projection."""
        coordinate0,velocity0=self._check_coordinate(coordinate0),self._check_coordinate(velocity0)
        times=np.asarray(times,dtype=float)
        if times.ndim!=1 or not np.all(np.isfinite(times)):
            raise ValueError("Linear reference times must be finite and one-dimensional")
        modes=self.linear_eigenpairs()
        omega,vectors=modes["omega"],modes["vectors"]
        initial=vectors.T@(self.M0@coordinate0)
        speed=vectors.T@(self.M0@velocity0)
        angles=times[:,None]*omega[None,:]
        cosine,sine=np.cos(angles),np.sin(angles)
        coordinate=(cosine*initial+sine*(speed/omega))@vectors.T
        velocity=(-sine*(initial*omega)+cosine*speed)@vectors.T
        return {"times":times,"q":coordinate,"velocity":velocity}

    def diagnostics(self,coordinate,velocity=None):
        fields=self.reconstruct(coordinate)
        gradients=self.reconstruct(coordinate,derivative=1)
        us,ws,theta,c=gradients[:,0],gradients[:,1],fields[:,2],fields[:,3]
        gamma1=us+theta*ws-theta**2/2-us*theta**2/2
        gamma2=ws-theta-theta*us-ws*theta**2/2+theta**3/6
        theta_mass=self._theta_mass(coordinate)[0]
        mass0=self._constant_masses[2]
        eigenvalues=eigh(theta_mass,mass0,eigvals_only=True,check_finite=False)
        relative_min=min(1.,float(eigenvalues.min()))
        relative_max=max(1.,float(eigenvalues.max()))
        result={"min_one_plus_c":float(np.min(1+c)),"max_abs_c":float(np.max(np.abs(c))),
                "max_abs_theta":float(np.max(np.abs(theta))),
                "max_abs_retained_axial_strain":float(np.max(np.abs(gamma1))),
                "max_abs_retained_shear_strain":float(np.max(np.abs(gamma2))),
                "max_normalized_curvature":float(self.length*np.max(np.abs(gradients[:,2]))),
                "max_normalized_contraction_gradient":float(self.length*np.max(np.abs(gradients[:,3]))),
                "relative_mass_min_eigenvalue":relative_min,"relative_mass_max_eigenvalue":relative_max,
                "relative_mass_condition":relative_max/relative_min,
                "mass_positive":relative_min>0}
        if velocity is not None:
            result["energy"]=self.energy(coordinate,velocity)
        return result
