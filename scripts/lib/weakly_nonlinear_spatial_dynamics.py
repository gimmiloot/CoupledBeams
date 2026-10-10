"""Seven-field Shen discretization of the frozen quartic spatial action.

This helper compiles exact derivatives of a supplied, serialized action. It
does not derive a new action, integrate an IVP, solve a static problem, or
select an excitation. The resting-mass whitening is reversible and retains
all seven independent endpoint-constrained Shen spaces. The complete coupled
rotational mass and all kinetic coordinate forces are retained.
"""
from __future__ import annotations

import math
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import cho_factor, cho_solve, eigh, solve_triangular

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin, _CompiledPolynomials

VERSION = "spatial-quartic-shen-exact-kinetic-derivatives-v1"
FIELDS = rod.FIELD_ORDER
_ROT = (3, 4, 5)
_KIN_Q = (3, 4, 5, 6)
_POT_NAMES = ("u_s", "w_s", "v_s", "Phi", "psi", "theta", "c",
              "Phi_s", "psi_s", "theta_s", "c_s")
_POT_FIELDS = (0, 1, 2, 3, 4, 5, 6, 3, 4, 5, 6)
_POT_DERIV = (1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1)
_KIN_NAMES = tuple(FIELDS[i] for i in _KIN_Q) + tuple(f + "_t" for f in FIELDS)
_MASS_NAMES = tuple(FIELDS[i] for i in _KIN_Q)
_JET_NAMES = tuple(f + suffix for suffix in ("", "_s", "_t", "_ss", "_st", "_tt") for f in FIELDS)


class _SparseCompiled:
    """Use the established monomial compiler only for nonzero outputs."""

    def __init__(self, polynomials, names, coefficients):
        self.size = len(polynomials)
        self.indices = np.array([i for i, p in enumerate(polynomials) if p], dtype=int)
        self.compiled = (_CompiledPolynomials([polynomials[i] for i in self.indices], names, coefficients)
                         if len(self.indices) else None)
        self.names = tuple(names)

    def evaluate(self, variables):
        result = np.zeros((self.size, variables.shape[1]), dtype=np.asarray(variables).dtype)
        if self.compiled is not None:
            result[self.indices] = self.compiled.evaluate(variables)
        return result


class SpatialGalerkin:
    """Seven independent fields, canonical ordering, and exact local Jacobian.

    A supplied frozen model must contain T4, V4, and residual_a. No fallback
    symbolic derivation is permitted. p is the maximum polynomial degree,
    ndof=7*(p-1), and nq>=2*p+1 integrates the retained quartic action.
    """

    fields = FIELDS
    _raw_basis = PlanarGalerkin._raw_basis

    def __init__(self, coefficients: rod.RodCoefficients, p: int, length=1.,
                 nq=None, model=None, whiten=True):
        if model is None:
            raise ValueError("Supply the hash-verified frozen serialized quartic action")
        if not isinstance(p, (int, np.integer)) or p < 2:
            raise ValueError("Maximum polynomial degree p must be an integer >=2")
        if not math.isfinite(length) or length <= 0:
            raise ValueError("Rod length must be finite and positive")
        nq = 2*p+1 if nq is None else int(nq)
        if nq < 2*p+1:
            raise ValueError("Quartic action requires nq>=2*p+1; reduced integration is not allowed")
        self.coefficients, self.p, self.length, self.nq = coefficients, int(p), float(length), nq
        self.n, self.ndof, self.whiten = self.p-1, 7*(self.p-1), bool(whiten)
        self.slices = {f: slice(i*self.n, (i+1)*self.n) for i, f in enumerate(FIELDS)}
        self._indices = {f: np.arange(s.start, s.stop) for f, s in self.slices.items()}
        self._rot_ids = np.concatenate([self._indices[FIELDS[i]] for i in _ROT])
        xi, weights = leggauss(nq)
        self.x, self.weights = (xi+1)*self.length/2, weights*self.length/2
        raw = self._raw_basis(self.x, 0)
        gram = raw.T@(self.weights[:, None]*raw)
        masses = (coefficients.m, coefficients.m, coefficients.m,
                  coefficients.jp+coefficients.jb, coefficients.jb, coefficients.jp, coefficients.jp)
        self._transforms, self._constant_masses, self._constant_factors = [], [], []
        for mass in masses:
            resting = mass*gram
            lower = np.linalg.cholesky(resting)
            transform = solve_triangular(lower.T, np.eye(self.n), lower=False) if whiten else np.eye(self.n)
            constant = transform.T@resting@transform
            self._transforms.append(transform)
            self._constant_masses.append(constant)
            self._constant_factors.append(cho_factor(constant, lower=True, check_finite=False))
        self.B = tuple(raw@a for a in self._transforms)
        self.D = tuple(self._raw_basis(self.x, 1)@a for a in self._transforms)
        self.D2 = tuple(self._raw_basis(self.x, 2)@a for a in self._transforms)
        self.M0 = np.zeros((self.ndof, self.ndof))
        for f, m in zip(FIELDS, self._constant_masses):
            self.M0[self.slices[f], self.slices[f]] = m
        self.model = model
        potential, kinetic = model.V4, model.T4
        # Reject a kinetic artifact with gradients or nonquadratic velocities.
        velocity_indices = {rod.SYMBOL_ORDER.index(f+"_t") for f in FIELDS}
        allowed_indices = {rod.SYMBOL_ORDER.index(f) for f in _MASS_NAMES} | velocity_indices
        for monomial in kinetic.terms:
            if sum(i in velocity_indices for i in monomial) != 2:
                raise ValueError("Frozen T4 must be homogeneous of degree two in velocities")
            if any(i < 42 and i not in allowed_indices for i in monomial):
                raise ValueError("Unsupported frozen kinetic coordinate/gradient")
        gradient = [potential.derivative(f) for f in _POT_NAMES]
        hessian = [g.derivative(f) for g in gradient for f in _POT_NAMES]
        self._potential_energy = _SparseCompiled([potential], _POT_NAMES, coefficients)
        self._potential_gradient = _SparseCompiled([potential]+gradient, _POT_NAMES, coefficients)
        self._potential_hessian = _SparseCompiled(hessian, _POT_NAMES, coefficients)
        self._potential_matrices = tuple((self.D if d else self.B)[i] for i, d in zip(_POT_FIELDS, _POT_DERIV))
        momenta = [kinetic.derivative(f+"_t") for f in FIELDS]
        mass = [[a.derivative(f+"_t") for f in FIELDS] for a in momenta]
        inertia = [sum((a.derivative(f)*rod.Polynomial.symbol(f+"_t") for f in FIELDS), rod.Polynomial())
                   - kinetic.derivative(FIELDS[i]) for i, a in enumerate(momenta)]
        # Translation and c retain their constant independent mass blocks.
        for i in (0, 1, 2, 6):
            for j in range(7):
                expected = rod.Polynomial.symbol("m" if i < 3 else "jp") if i == j else rod.Polynomial()
                if mass[i][j] != expected or mass[j][i] != expected:
                    raise ValueError("Frozen kinetic mass differs from its constant translation/contraction blocks")
        self._kinetic = _SparseCompiled([kinetic]+[kinetic.derivative(f) for f in _MASS_NAMES]+momenta,
                                       _KIN_NAMES, coefficients)
        self._mass = _SparseCompiled([mass[i][j] for i in _ROT for j in _ROT], _MASS_NAMES, coefficients)
        self._inertia = _SparseCompiled(inertia, _KIN_NAMES, coefficients)
        self._inertia_q = _SparseCompiled([g.derivative(f) for g in inertia for f in _MASS_NAMES],
                                         _KIN_NAMES, coefficients)
        self._inertia_v = _SparseCompiled([g.derivative(f+"_t") for g in inertia for f in FIELDS],
                                         _KIN_NAMES, coefficients)
        self._mass_q = _SparseCompiled([mass[i][j].derivative(f) for i in _ROT for j in _ROT for f in _MASS_NAMES],
                                      _MASS_NAMES, coefficients)
        self._residual = _SparseCompiled(list(model.residual_a), _JET_NAMES, coefficients)
        self._potential_cache_q = self._potential_cache = None
        self._mass_cache_q = self._mass_cache_matrix = self._mass_cache_factor = None
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
        self.rhs_calls = self.jacobian_calls = self.mass_factorizations = 0
        self.force_evaluations = self.linear_eigendecompositions = 0

    def counters(self):
        return {name: getattr(self, name) for name in
                ("rhs_calls", "jacobian_calls", "mass_factorizations", "force_evaluations", "linear_eigendecompositions")}

    def _check_coordinate(self, value):
        value = np.asarray(value)
        if value.shape != (self.ndof,):
            raise ValueError(f"Expected {self.ndof} coefficient coordinates")
        return value

    def basis_at(self, points, derivative=0):
        raw = self._raw_basis(points, derivative)
        return {f: raw@a for f, a in zip(FIELDS, self._transforms)}

    def raw_coefficients(self, coordinate):
        coordinate = self._check_coordinate(coordinate)
        return np.concatenate([a@coordinate[self.slices[f]] for f, a in zip(FIELDS, self._transforms)])

    def from_raw_coefficients(self, raw):
        raw = self._check_coordinate(raw)
        return np.concatenate([solve_triangular(a, raw[self.slices[f]], lower=False, check_finite=False)
                               for f, a in zip(FIELDS, self._transforms)])

    def reconstruct(self, coordinate, points=None, derivative=0):
        coordinate = self._check_coordinate(coordinate)
        if points is None:
            matrices = self.B if derivative == 0 else self.D if derivative == 1 else self.D2 if derivative == 2 else None
            if matrices is None:
                raise ValueError("Supported derivative orders are 0,1,2")
        else:
            matrices = tuple(self.basis_at(points, derivative).values())
        return np.column_stack([a@coordinate[self.slices[f]] for f, a in zip(FIELDS, matrices)])

    def reconstruct_series(self, rows, points=None, derivative=0):
        rows = np.asarray(rows)
        if rows.ndim != 2 or rows.shape[1] != self.ndof:
            raise ValueError("Coefficient series must have shape (nt,ndof)")
        matrices = self.basis_at(self.x if points is None else points, derivative)
        return np.stack([rows[:, self.slices[f]]@matrices[f].T for f in FIELDS], axis=2)

    def project(self, values):
        values = np.asarray(values(self.x) if callable(values) else values)
        if values.shape != (self.nq, 7):
            raise ValueError("Physical projection fields must have shape (nq,7)")
        coordinate = np.empty(self.ndof)
        for i, f in enumerate(FIELDS):
            a = self.B[i]
            gram = a.T@(self.weights[:, None]*a)
            coordinate[self.slices[f]] = cho_solve(cho_factor(gram, lower=True, check_finite=False),
                                                   a.T@(self.weights*values[:, i]), check_finite=False)
        return coordinate

    def _local_potential(self, coordinate, gradient=True, hessian=False):
        coordinate = self._check_coordinate(coordinate)
        if self._potential_cache_q is None or not np.array_equal(coordinate, self._potential_cache_q):
            self.force_evaluations += 1
            variables = np.vstack([a@coordinate[self.slices[FIELDS[i]]] for i, a in zip(_POT_FIELDS, self._potential_matrices)])
            self._potential_cache_q = coordinate.copy()
            self._potential_cache = {"variables": variables, "energy": None, "gradient": None, "hessian": None}
        cache = self._potential_cache
        if gradient and cache["gradient"] is None:
            values = self._potential_gradient.evaluate(cache["variables"])
            cache["energy"], cache["gradient"] = float(self.weights@values[0]), values[1:]
        elif cache["energy"] is None:
            cache["energy"] = float(self.weights@self._potential_energy.evaluate(cache["variables"])[0])
        if hessian and cache["hessian"] is None:
            cache["hessian"] = self._potential_hessian.evaluate(cache["variables"]).reshape(11, 11, self.nq)
        return cache["energy"], cache["gradient"], cache["hessian"]

    def potential(self, coordinate, gradient=True, hessian=False):
        energy, local_gradient, local_hessian = self._local_potential(coordinate, gradient, hessian)
        result = {"V": energy}
        if gradient:
            force = np.zeros(self.ndof)
            for i, a, value in zip(_POT_FIELDS, self._potential_matrices, local_gradient):
                force[self.slices[FIELDS[i]]] += a.T@(self.weights*value)
            result["gradient"] = force
        if hessian:
            tangent = np.zeros((self.ndof, self.ndof))
            for k in range(11):
                ls, left = self.slices[FIELDS[_POT_FIELDS[k]]], self._potential_matrices[k]
                for j in range(k, 11):
                    if not np.any(local_hessian[k, j]):
                        continue
                    rs, right = self.slices[FIELDS[_POT_FIELDS[j]]], self._potential_matrices[j]
                    block = left.T@((self.weights*local_hessian[k, j])[:, None]*right)
                    tangent[ls, rs] += block
                    if k != j:
                        tangent[rs, ls] += block.T
            result["hessian"] = tangent
        return result

    def _mass_variables(self, coordinate):
        return np.vstack([self.B[i]@coordinate[self.slices[FIELDS[i]]] for i in _KIN_Q])

    def _kinetic_variables(self, coordinate, velocity):
        return np.vstack((self._mass_variables(coordinate), self.reconstruct(velocity).T))

    def _rotational_mass(self, coordinate):
        key = np.concatenate([coordinate[self.slices[FIELDS[i]]] for i in _KIN_Q])
        if self._mass_cache_q is not None and np.array_equal(key, self._mass_cache_q):
            return self._mass_cache_matrix, self._mass_cache_factor
        variables = self._mass_variables(coordinate)
        if not np.isfinite(variables).all() or np.min(1+variables[3]) <= 0:
            raise FloatingPointError("Contraction scale 1+c must remain finite and positive")
        values = self._mass.evaluate(variables).reshape(3, 3, self.nq)
        matrix = np.empty((3*self.n, 3*self.n))
        for i, field_i in enumerate(_ROT):
            for j, field_j in enumerate(_ROT):
                matrix[i*self.n:(i+1)*self.n, j*self.n:(j+1)*self.n] = (
                    self.B[field_i].T@((self.weights*values[i, j])[:, None]*self.B[field_j]))
        factor = cho_factor(matrix, lower=True, check_finite=False)
        self.mass_factorizations += 1
        self._mass_cache_q, self._mass_cache_matrix, self._mass_cache_factor = key.copy(), matrix, factor
        return matrix, factor

    def mass_matrix(self, coordinate):
        coordinate = self._check_coordinate(coordinate)
        result = self.M0.copy()
        result[np.ix_(self._rot_ids, self._rot_ids)] = self._rotational_mass(coordinate)[0]
        return result

    def _solve_mass(self, coordinate, rhs):
        result = np.empty_like(rhs)
        for i in (0, 1, 2, 6):
            s = self.slices[FIELDS[i]]
            result[s] = cho_solve(self._constant_factors[i], rhs[s], check_finite=False)
        result[self._rot_ids] = cho_solve(self._rotational_mass(coordinate)[1], rhs[self._rot_ids], check_finite=False)
        return result

    def inertial_terms(self, coordinate, velocity):
        coordinate, velocity = self._check_coordinate(coordinate), self._check_coordinate(velocity)
        local = self._inertia.evaluate(self._kinetic_variables(coordinate, velocity))
        return np.concatenate([a.T@(self.weights*value) for a, value in zip(self.B, local)])

    def kinetic(self, coordinate, velocity):
        coordinate, velocity = self._check_coordinate(coordinate), self._check_coordinate(velocity)
        local = self._kinetic.evaluate(self._kinetic_variables(coordinate, velocity))
        return float(self.weights@local[0])

    def acceleration(self, coordinate, velocity, linear=False):
        coordinate, velocity = self._check_coordinate(coordinate), self._check_coordinate(velocity)
        if linear:
            force = -self.K@coordinate
            return np.concatenate([cho_solve(a, force[self.slices[f]], check_finite=False)
                                   for f, a in zip(FIELDS, self._constant_factors)])
        return self._solve_mass(coordinate, -self.potential(coordinate)["gradient"]-self.inertial_terms(coordinate, velocity))

    def rhs(self, time, state):
        del time
        self.rhs_calls += 1
        state = np.asarray(state)
        if state.shape != (2*self.ndof,) or not np.isfinite(state).all():
            raise FloatingPointError("Dynamics state must contain finite coordinates and velocities")
        return np.concatenate((state[self.ndof:], self.acceleration(state[:self.ndof], state[self.ndof:])))

    def jacobian(self, time, state):
        """Exact action derivatives, including coordinate-dependent mass times a."""
        del time
        self.jacobian_calls += 1
        state = np.asarray(state)
        if state.shape != (2*self.ndof,) or not np.isfinite(state).all():
            raise FloatingPointError("Dynamics state must contain finite coordinates and velocities")
        coordinate, velocity = state[:self.ndof], state[self.ndof:]
        potential = self.potential(coordinate, hessian=True)
        acceleration = self._solve_mass(coordinate, -potential["gradient"]-self.inertial_terms(coordinate, velocity))
        variables = self._kinetic_variables(coordinate, velocity)
        gq = self._inertia_q.evaluate(variables).reshape(7, 4, self.nq)
        gv = self._inertia_v.evaluate(variables).reshape(7, 7, self.nq)
        mq = self._mass_q.evaluate(variables[:4]).reshape(3, 3, 4, self.nq)
        physical_acceleration = self.reconstruct(acceleration)
        for i, field_i in enumerate(_ROT):
            gq[field_i] += np.einsum("jkx,xj->kx", mq[i], physical_acceleration[:, _ROT], optimize=True)
        tangent = potential["hessian"].copy()
        velocity_tangent = np.zeros((self.ndof, self.ndof))
        for i in range(7):
            ls, left = self.slices[FIELDS[i]], self.B[i]
            for k, field_k in enumerate(_KIN_Q):
                if np.any(gq[i, k]):
                    tangent[ls, self.slices[FIELDS[field_k]]] += left.T@((self.weights*gq[i, k])[:, None]*self.B[field_k])
            for j in range(7):
                if np.any(gv[i, j]):
                    velocity_tangent[ls, self.slices[FIELDS[j]]] += left.T@((self.weights*gv[i, j])[:, None]*self.B[j])
        result = np.zeros((2*self.ndof, 2*self.ndof))
        result[:self.ndof, self.ndof:] = np.eye(self.ndof)
        result[self.ndof:, :self.ndof] = self._solve_mass(coordinate, -tangent)
        result[self.ndof:, self.ndof:] = self._solve_mass(coordinate, -velocity_tangent)
        return result

    def linear_rhs(self, time, state):
        del time
        state = np.asarray(state)
        if state.shape != (2*self.ndof,):
            raise ValueError("Linear state has the wrong size")
        return np.concatenate((state[self.ndof:], self.acceleration(state[:self.ndof], state[self.ndof:], linear=True)))

    def linear_jacobian(self, time=None, state=None):
        del time, state
        result = np.zeros((2*self.ndof, 2*self.ndof))
        result[:self.ndof, self.ndof:] = np.eye(self.ndof)
        for f, factor in zip(FIELDS, self._constant_factors):
            s = self.slices[f]
            result[self.ndof+s.start:self.ndof+s.stop, :self.ndof] = cho_solve(factor, -self.K[s], check_finite=False)
        return result

    def energy(self, coordinate, velocity):
        return self.kinetic(coordinate, velocity)+self.potential(coordinate, gradient=False)["V"]

    def energy_rate(self, coordinate, velocity, acceleration=None):
        coordinate, velocity = self._check_coordinate(coordinate), self._check_coordinate(velocity)
        acceleration = self.acceleration(coordinate, velocity) if acceleration is None else self._check_coordinate(acceleration)
        values = self._kinetic.evaluate(self._kinetic_variables(coordinate, velocity))
        physical_v, physical_a = self.reconstruct(velocity), self.reconstruct(acceleration)
        power = sum(values[1+k]*physical_v[:, i] for k, i in enumerate(_KIN_Q))
        power += np.sum(values[5:]*physical_a.T, axis=0)
        return float(self.weights@power+velocity@self.potential(coordinate)["gradient"])

    def weak_residual(self, coordinate, velocity, acceleration):
        local = []
        for values, derivative in ((coordinate, 0), (coordinate, 1), (velocity, 0),
                                  (coordinate, 2), (velocity, 1), (acceleration, 0)):
            local.extend(self.reconstruct(values, derivative=derivative).T)
        residual = self._residual.evaluate(np.asarray(local))
        return np.concatenate([a.T@(self.weights*value) for a, value in zip(self.B, residual)])

    def linear_eigenpairs(self, block=None):
        blocks = {"mh": ("u", "c"), "timoshenko": ("w", "theta"),
                  "outplane": ("v", "psi"), "torsion": ("Phi",)}
        if block is not None and block not in blocks:
            raise ValueError("Linear block must be mh,timoshenko,outplane,torsion or None")
        if block not in self._linear_modes:
            names = FIELDS if block is None else blocks[block]
            ids = np.concatenate([self._indices[f] for f in names])
            values, vectors = eigh(self.K[np.ix_(ids, ids)], self.M0[np.ix_(ids, ids)], check_finite=False)
            if np.any(values <= 0) or not np.isfinite(values).all():
                raise ArithmeticError("Fixed-fixed linear action has a nonpositive eigenvalue")
            expanded = np.zeros((self.ndof, len(values)))
            expanded[ids] = vectors
            omega = np.sqrt(values)
            self._linear_modes[block] = {"eigenvalues": values, "omega": omega,
                                         "frequency_hz": omega/(2*np.pi), "vectors": expanded}
            self.linear_eigendecompositions += 1
        return self._linear_modes[block]

    def linear_reference(self, coordinate0, velocity0, times):
        coordinate0, velocity0 = self._check_coordinate(coordinate0), self._check_coordinate(velocity0)
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not np.isfinite(times).all():
            raise ValueError("Linear reference times must be finite and one-dimensional")
        modes = self.linear_eigenpairs()
        omega, vectors = modes["omega"], modes["vectors"]
        initial, speed = vectors.T@(self.M0@coordinate0), vectors.T@(self.M0@velocity0)
        angles = times[:, None]*omega
        cosine, sine = np.cos(angles), np.sin(angles)
        return {"times": times,
                "q": (cosine*initial+sine*(speed/omega))@vectors.T,
                "velocity": (-sine*(initial*omega)+cosine*speed)@vectors.T}

    def mass_spectral_bounds(self, coordinate):
        """Pointwise Loewner bounds relative to the resting rotational metric.

        The same positive spatial quadrature preserves these bounds for the
        assembled Galerkin mass. They are conservative bounds, not an extra
        continuum eigenvalue solve or a change to the solved variable mass.
        """
        variables = self._mass_variables(self._check_coordinate(coordinate))
        local = self._mass.evaluate(variables).reshape(3, 3, self.nq).transpose(2, 0, 1)
        resting = np.sqrt((self.coefficients.jp+self.coefficients.jb,
                           self.coefficients.jb, self.coefficients.jp))
        relative = local/resting[None, :, None]/resting[None, None, :]
        values = np.linalg.eigvalsh(relative)
        low, high = min(1., float(values.min())), max(1., float(values.max()))
        return {"relative_mass_lower_bound": low, "relative_mass_upper_bound": high,
                "relative_mass_condition_upper_bound": high/low if low > 0 else float("inf"),
                "mass_positive": low > 0,
                "method": "positive-quadrature Loewner bounds of local rotational quartic mass"}

    def diagnostics(self, coordinate, velocity=None):
        values, gradients = self.reconstruct(coordinate), self.reconstruct(coordinate, derivative=1)
        a = values[:, (3, 4, 5)]*np.array((1., -1., 1.))
        ass = gradients[:, (3, 4, 5)]*np.array((1., -1., 1.))
        qs = gradients[:, :3]
        e1 = np.zeros_like(qs); e1[:, 0] = 1.
        cross = np.cross
        gamma = qs-cross(a, e1+qs)+cross(a, cross(a, e1+qs))/2-cross(a, cross(a, cross(a, e1)))/6
        chi = ass-cross(a, ass)/2+cross(a, cross(a, ass))/6
        relative = eigh(self._rotational_mass(coordinate)[0], self.M0[np.ix_(self._rot_ids, self._rot_ids)],
                        eigvals_only=True, check_finite=False)
        low, high = min(1., float(relative.min())), max(1., float(relative.max()))
        result = {"min_one_plus_c": float(np.min(1+values[:, 6])), "max_abs_c": float(np.max(abs(values[:, 6]))),
                  "max_abs_theta": float(np.max(abs(values[:, 5]))), "max_abs_rotation": float(np.max(np.linalg.norm(a, axis=1))),
                  "max_abs_axial_gradient": float(np.max(abs(qs[:, 0]))),
                  "max_abs_transverse_gradient": float(np.max(abs(qs[:, 1:]))),
                  "max_abs_retained_axial_strain": float(np.max(abs(gamma[:, 0]))),
                  "max_abs_retained_shear_strain": float(np.max(abs(gamma[:, 1:]))),
                  "max_normalized_curvature": float(self.length*np.max(abs(chi))),
                  "max_L_abs_curvature": float(self.length*np.max(abs(ass))),
                  "max_normalized_contraction_gradient": float(self.length*np.max(abs(gradients[:, 6]))),
                  "relative_mass_min_eigenvalue": low, "relative_mass_max_eigenvalue": high,
                  "relative_mass_condition": high/low, "mass_positive": low > 0}
        if velocity is not None:
            result["energy"] = self.energy(coordinate, velocity)
        return result
