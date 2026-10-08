"""Exact-time diagnostic of the leading axial response of the audited action.

This companion does not advance the nonlinear model. Its prescribed bending
background is the same continuous analytical Timoshenko eigenpair for every
Shen resolution. All coordinates of the two-field M-H Galerkin space remain
in the generalized eigendecomposition used to evaluate its matrix function.
No ODE integrator, filtering, modal reduction or initial correction is used.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import weakly_nonlinear_planar_dynamics as planar
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section

VERSION = "audited-planar-second-order-axial-exact-time-v1"
FIELDS = ("u", "c")
_SUFFIXES = ("", "_s", "_t", "_ss", "_st", "_tt")
_DERIVATION_MODEL = None
_DERIVATION_VALUE = None


def _weighted_order(polynomial, order):
    """Extract epsilon_a order, with u,c order two and w,theta order one."""
    weights = {name+suffix: weight for name, weight in
               (("u", 2), ("c", 2), ("w", 1), ("theta", 1))
               for suffix in _SUFFIXES}
    return rod.Polynomial({monomial: coefficient
                           for monomial, coefficient in polynomial.terms.items()
                           if sum(weights.get(rod.SYMBOL_ORDER[index], 0)
                                  for index in monomial) == order})


def derive_second_order(model=None):
    """Exact polynomial extraction and sign checks, cached once per action."""
    global _DERIVATION_MODEL, _DERIVATION_VALUE
    model = rod.derive_polynomials() if model is None else model
    if model is _DERIVATION_MODEL:
        return _DERIVATION_VALUE
    inactive = {name+suffix: 0 for name in ("v", "Phi", "psi") for suffix in _SUFFIXES}
    residuals = tuple(expression.substitute(inactive) for expression in model.residual_a)
    u_equation = _weighted_order(residuals[0], 2)
    c_equation = _weighted_order(residuals[6], 2)
    bending_second_order = tuple(_weighted_order(residuals[index], 2) for index in (1, 5))
    bending_zero = {name+suffix: 0 for name in ("w", "theta") for suffix in _SUFFIXES}
    axial_zero = {name+suffix: 0 for name in FIELDS for suffix in _SUFFIXES}
    linear_u, linear_c = u_equation.substitute(bending_zero), c_equation.substitute(bending_zero)
    forcing_u, forcing_c = -u_equation.substitute(axial_zero), -c_equation.substitute(axial_zero)
    s = model.symbols
    Z = ((s["C"]-s["S"])*s["theta"]*s["w_s"]
         +(s["S"]-s["C"]/2)*s["theta"]**2)
    D = s["theta"]*s["w_s"]-s["theta"]**2/2
    checks = {
        "linear_u": linear_u == s["m"]*s["u_tt"]-s["C"]*s["u_ss"]-s["nu"]*s["C"]*s["c_s"],
        "linear_c": linear_c == s["jp"]*s["c_tt"]-s["H"]*s["c_ss"]+s["C"]*(s["c"]+s["nu"]*s["u_s"]),
        "forcing_u": forcing_u == Z.total_derivative("s"),
        "forcing_c": forcing_c == -s["nu"]*s["C"]*D+s["jp"]*s["theta_t"]**2,
        "bending_second_order_zero": all(expression == 0 for expression in bending_second_order),
        "unknowns_absent_from_forcing": all(
            not any(rod.SYMBOL_ORDER[index] in axial_zero for index in monomial)
            for expression in (forcing_u, forcing_c) for monomial in expression.terms),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"Protected-action second-order extraction failed: {checks}")
    value = {"u_equation": u_equation, "c_equation": c_equation,
             "linear_u": linear_u, "linear_c": linear_c,
             "forcing_u": forcing_u, "forcing_c": forcing_c,
             "Z": Z, "D": D, "bending_second_order": bending_second_order,
             "checks": checks, "model_version": rod.MODEL_VERSION}
    _DERIVATION_MODEL, _DERIVATION_VALUE = model, value
    return value


def sinc_unscaled(argument):
    """sin(x)/x with its removable value at zero, no pi normalization."""
    argument = np.asarray(argument, dtype=float)
    result = np.ones_like(argument)
    np.divide(np.sin(argument), argument, out=result, where=argument != 0)
    return result


def response_kernel(omega, driving_omega, times, derivative=0):
    """Stable zero-IC response to cos(driving_omega*t), and time derivatives.

    Broadcasting is explicit NumPy broadcasting: use omega[None,:] and
    times[:,None] for a time-by-coordinate matrix. No approximate detuning
    is replaced by zero, including when the arguments are near resonance.
    """
    if derivative not in (0, 1, 2):
        raise ValueError("Analytical time derivative must be 0,1 or 2")
    omega, driving_omega, times = np.broadcast_arrays(
        np.asarray(omega, dtype=float), np.asarray(driving_omega, dtype=float),
        np.asarray(times, dtype=float))
    if not np.all(np.isfinite(omega)) or not np.all(np.isfinite(driving_omega)) or not np.all(np.isfinite(times)):
        raise ValueError("Response arguments must be finite")
    if np.any(omega < 0) or np.any(driving_omega < 0):
        raise ValueError("Response frequencies must be nonnegative")
    a, b = (omega+driving_omega)/2, (omega-driving_omega)/2
    at, bt = a*times, b*times
    sa, sb = sinc_unscaled(at), sinc_unscaled(bt)
    if derivative == 1:
        return times/2*(np.cos(at)*sb+sa*np.cos(bt))
    value = times**2/2*sa*sb
    if derivative == 2:
        return np.cos(at)*np.cos(bt)-(a*a+b*b)*value
    return value


@dataclass(frozen=True)
class AnalyticBackground:
    """Common continuous first bending mode, scaled by h0, not by A."""
    source_model: object
    length: float
    omega: float
    h0: float
    normalized_coefficients: np.ndarray
    source_bundle: str
    normalization: str = "one common signed w(L/2) normalization, max|w_hat|=1"

    @property
    def T1(self):
        return 2*math.pi/self.omega

    def evaluate(self, points, derivative=0):
        points = np.asarray(points, dtype=float)
        if derivative not in (0, 1, 2):
            raise ValueError("Background derivative must be 0,1 or 2")
        if points.ndim != 1 or np.any(points < 0) or np.any(points > self.length):
            raise ValueError("Background points must lie in [0,L]")
        order = min(derivative, 1)
        states = (mh.finite_state_basis(self.source_model, self.length, self.omega,
                                      points, "timoshenko", order)
                  @ self.normalized_coefficients)
        if derivative == 2:
            states = states @ mh.harmonic_state_matrix(self.source_model, self.omega, "timoshenko").T
        return self.h0*states[:, :2]

    def as_dict(self):
        return {"omega1": self.omega, "T1": self.T1, "h0": self.h0, "L": self.length,
                "source_bundle": self.source_bundle, "normalization": self.normalization,
                "normalized_analytic_coefficients": self.normalized_coefficients.tolist(),
                "amplitude_convention": "epsilon_a=A/h0; W=h0*w_hat, Theta=h0*theta_hat; physical(u,c)=epsilon_a^2*(u2,c2)"}


def background_from_pilot(config, root=None):
    """Recover the original analytic eigenpair, without new frequency roots."""
    root = Path(__file__).resolve().parents[2] if root is None else Path(root)
    geometry = config["material_geometry"]
    reference_path = root/config["linear_reference_bundle"]
    reference = json.loads((reference_path/"result.json").read_text(encoding="utf8"))
    omega = reference["timoshenko"]["roots"][0]["omega"]
    section = rectangular_section(E=geometry["E"], rho=geometry["rho"], nu=geometry["nu"],
                                  width=geometry["b"], thickness=geometry["h"], K=5/6)
    model = mh.project_jang_reduced_rectangular(section)
    mode = mh.finite_mode(model, geometry["L"], omega, "timoshenko")
    peak = (mh.finite_state_basis(model, geometry["L"], omega, [geometry["L"]/2], "timoshenko")
            @ mode["coefficients"])[0, 0]
    background = AnalyticBackground(model, geometry["L"], omega, geometry["h"],
                                   np.asarray(mode["coefficients"])/peak,
                                   config["linear_reference_bundle"])
    grid = np.linspace(0, background.length, 501)
    if abs(np.max(abs(background.evaluate(grid)[:, 0]))/background.h0-1) > 2e-12:
        raise ArithmeticError("Common continuous first-mode normalization failed")
    return background


class SecondOrderAxial:
    """All 2(p-1) coordinates of the linear forced M-H diagnostic."""
    fields = FIELDS

    def __init__(self, coefficients: rod.RodCoefficients, p: int,
                 background: AnalyticBackground, model=None, nq=None):
        if not isinstance(p, (int, np.integer)) or p < 2:
            raise ValueError("Maximum degree p must be an integer >=2")
        self.coefficients, self.p, self.background = coefficients, int(p), background
        self.length, self.n, self.ndof = background.length, self.p-1, 2*(self.p-1)
        self.matrix_nq = 2*self.p+1
        self.nq = self.matrix_nq if nq is None else int(nq)
        self.model = rod.derive_polynomials() if model is None else model
        self.derivation = derive_second_order(self.model)
        self.disc = planar.PlanarGalerkin(coefficients, self.p, length=self.length,
                                        nq=self.matrix_nq, model=self.model)
        self.slices = {"u": slice(0, self.n), "c": slice(self.n, self.ndof)}
        self.transforms = tuple(self.disc._transforms[index] for index in (0, 3))
        self.B = tuple(self.disc.B[index] for index in (0, 3))
        self.D = tuple(self.disc.D[index] for index in (0, 3))
        self.points, self.weights = self.disc.x, self.disc.weights
        Bu, Bc = self.B
        Du, Dc = self.D
        weights, pcoef = self.weights, self.coefficients
        self.M = np.zeros((self.ndof, self.ndof))
        self.K = np.zeros_like(self.M)
        self.M[self.slices["u"], self.slices["u"]] = pcoef.m*Bu.T@(weights[:, None]*Bu)
        self.M[self.slices["c"], self.slices["c"]] = pcoef.jp*Bc.T@(weights[:, None]*Bc)
        self.K[self.slices["u"], self.slices["u"]] = pcoef.C*Du.T@(weights[:, None]*Du)
        self.K[self.slices["u"], self.slices["c"]] = pcoef.nu*pcoef.C*Du.T@(weights[:, None]*Bc)
        self.K[self.slices["c"], self.slices["u"]] = self.K[self.slices["u"], self.slices["c"]].T
        self.K[self.slices["c"], self.slices["c"]] = pcoef.C*Bc.T@(weights[:, None]*Bc)+pcoef.H*Dc.T@(weights[:, None]*Dc)
        indices = np.r_[self.disc._indices["u"], self.disc._indices["c"]]
        old_mass, old_stiffness = self.disc.M0[np.ix_(indices, indices)], self.disc.K[np.ix_(indices, indices)]
        self.matrix_checks = {
            "mass_max_absolute_difference": float(np.max(abs(self.M-old_mass))),
            "mass_relative_difference": float(np.linalg.norm(self.M-old_mass)/np.linalg.norm(old_mass)),
            "stiffness_max_absolute_difference": float(np.max(abs(self.K-old_stiffness))),
            "stiffness_relative_difference": float(np.linalg.norm(self.K-old_stiffness)/np.linalg.norm(old_stiffness)),
            "mass_relative_symmetry": float(np.linalg.norm(self.M-self.M.T)/np.linalg.norm(self.M)),
            "stiffness_relative_symmetry": float(np.linalg.norm(self.K-self.K.T)/np.linalg.norm(self.K)),
        }
        np.linalg.cholesky(self.M)
        np.linalg.cholesky(self.K)
        self.f0, self.f2 = self.assemble_forcing(self.nq)
        self.eigenvalues, self.vectors = eigh(self.K, self.M, check_finite=False)
        self.eigen_decompositions = 1
        if np.any(self.eigenvalues <= 0) or not np.all(np.isfinite(self.eigenvalues)):
            raise ArithmeticError("Fixed-fixed M-H diagnostic has nonpositive eigenvalues")
        self.omega = np.sqrt(self.eigenvalues)
        self.driving_omega = 2*self.background.omega
        self.b0, self.b2 = self.vectors.T@self.f0, self.vectors.T@self.f2
        residual = self.K@self.vectors-self.M@self.vectors*self.eigenvalues[None, :]
        self.eigen_checks = {
            "coordinates_retained": self.ndof, "eigenvectors_retained": self.vectors.shape[1],
            "mass_orthogonality_max_absolute": float(np.max(abs(self.vectors.T@self.M@self.vectors-np.eye(self.ndof)))),
            "eigenpair_relative_residual": float(np.linalg.norm(residual)/(np.linalg.norm(self.K@self.vectors)+np.linalg.norm(self.M@self.vectors*self.eigenvalues[None, :]))),
            "omega_min": float(self.omega[0]), "omega_max": float(self.omega[-1]),
            "minimum_absolute_detuning": float(np.min(abs(self.omega-self.driving_omega))),
        }
        self.exact_time_evaluations = 0


    def basis_at(self, points, derivative=0):
        values = self.disc.basis_at(points, derivative)
        return {name: values[name] for name in FIELDS}

    def assemble_forcing(self, nq, strong=False):
        """Continuous analytic forcing, in the same physical coefficient basis."""
        if int(nq) < self.p+1:
            raise ValueError("Forcing quadrature must resolve the test polynomial")
        nodes, weights = np.polynomial.legendre.leggauss(int(nq))
        points, weights = (nodes+1)*self.length/2, weights*self.length/2
        basis, gradients = self.basis_at(points), self.basis_at(points, 1)
        W, Theta = self.background.evaluate(points).T
        Ws, Thetas = self.background.evaluate(points, 1).T
        Wss = self.background.evaluate(points, 2)[:, 0]
        p = self.coefficients
        Z = (p.C-p.S)*Theta*Ws+(p.S-p.C/2)*Theta**2
        D = Theta*Ws-Theta**2/2
        Fc0 = (-p.nu*p.C*D+p.jp*self.background.omega**2*Theta**2)/2
        Fc2 = (-p.nu*p.C*D-p.jp*self.background.omega**2*Theta**2)/2
        if strong:
            Zs = (p.C-p.S)*(Thetas*Ws+Theta*Wss)+(2*p.S-p.C)*Theta*Thetas
            axial = basis["u"].T@(weights*Zs/2)
        else:
            axial = -gradients["u"].T@(weights*Z/2)
        return np.r_[axial, basis["c"].T@(weights*Fc0)], np.r_[axial, basis["c"].T@(weights*Fc2)]

    def forcing_checks(self, quadratures=None):
        quadratures = (self.nq, 3*self.p+7, 4*self.p+9) if quadratures is None else tuple(quadratures)
        rows = []
        previous = None
        for nq in quadratures:
            weak = self.assemble_forcing(nq)
            strong = self.assemble_forcing(nq, strong=True)
            row = {"nq": nq, "relative_changes": {}, "strong_weak_relative": {}, "strong_weak_absolute": {}}
            for name, source, strong_source in zip(("f0", "f2"), weak, strong):
                scale = max(np.linalg.norm(source), 1e-30)
                row["strong_weak_relative"][name] = float(np.linalg.norm(source-strong_source)/scale)
                row["strong_weak_absolute"][name] = float(np.max(abs(source-strong_source)))
                if previous is not None:
                    row["relative_changes"][name] = float(np.linalg.norm(source-previous[name])/scale)
            rows.append(row)
            previous = dict(zip(("f0", "f2"), weak))
        return {"quadrature_is_exact_for_matrices_not_analytic_background": True, "rows": rows}

    def modal_response(self, times, derivative=0):
        times = np.asarray(times, dtype=float)
        if times.ndim != 1 or not np.all(np.isfinite(times)):
            raise ValueError("Exact-time evaluation requires a finite one-dimensional time grid")
        return (response_kernel(self.omega[None, :], 0., times[:, None], derivative)*self.b0[None, :]
                +response_kernel(self.omega[None, :], self.driving_omega, times[:, None], derivative)*self.b2[None, :])

    def evaluate(self, times, derivative=0):
        self.exact_time_evaluations += 1
        return self.modal_response(times, derivative)@self.vectors.T

    def reconstruct(self, coordinate, points=None, derivative=0):
        coordinate = np.asarray(coordinate)
        if coordinate.shape != (self.ndof,):
            raise ValueError("Expected all two-field M-H coordinates")
        basis = self.B if points is None and derivative == 0 else tuple(self.basis_at(self.points if points is None else points, derivative).values())
        return np.column_stack([matrix@coordinate[self.slices[name]] for name, matrix in zip(FIELDS, basis)])

    def reconstruct_series(self, coordinates, points=None, derivative=0):
        coordinates = np.asarray(coordinates)
        if coordinates.ndim != 2 or coordinates.shape[1] != self.ndof:
            raise ValueError("Coefficient histories must retain every coordinate")
        basis = self.B if points is None and derivative == 0 else tuple(self.basis_at(self.points if points is None else points, derivative).values())
        return np.stack([coordinates[:, self.slices[name]]@matrix.T for name, matrix in zip(FIELDS, basis)], axis=2)

    def forcing(self, times):
        times = np.asarray(times, dtype=float)
        return self.f0[None, :]+np.cos(self.driving_omega*times[:, None])*self.f2[None, :]

    def equation_and_power_checks(self, times):
        coordinate, velocity, acceleration = (self.evaluate(times, order) for order in (0, 1, 2))
        mass_term, stiffness_term, source = acceleration@self.M.T, coordinate@self.K.T, self.forcing(times)
        residual = mass_term+stiffness_term-source
        scale = np.linalg.norm(mass_term, axis=1)+np.linalg.norm(stiffness_term, axis=1)+np.linalg.norm(source, axis=1)
        energy_rate = np.einsum("ti,ti->t", velocity, mass_term+stiffness_term)
        forcing_power = np.einsum("ti,ti->t", velocity, source)
        power_scale = np.linalg.norm(velocity, axis=1)*scale
        return {"equation_max_absolute": float(np.max(abs(residual))),
                "equation_scaled_max": float(np.max(np.linalg.norm(residual, axis=1)/np.maximum(scale, 1e-30))),
                "power_identity_max_absolute": float(np.max(abs(energy_rate-forcing_power))),
                "power_identity_scaled_max": float(np.max(abs(energy_rate-forcing_power)/np.maximum(power_scale, 1e-30))),
                "initial_coordinate_max_absolute": float(np.max(abs(self.evaluate([0.], 0)))),
                "initial_velocity_max_absolute": float(np.max(abs(self.evaluate([0.], 1)))),
                "initial_acceleration_relative_residual": float(np.linalg.norm(self.M@self.evaluate([0.], 2)[0]-self.f0-self.f2)/max(np.linalg.norm(self.f0+self.f2), 1e-30))}

    def counters(self):
        return {"mh_generalized_eigendecompositions": self.eigen_decompositions,
                "exact_time_evaluations": self.exact_time_evaluations,
                "coordinates_retained": self.ndof, "time_integrations": 0,
                "modal_reduction": False, "root_solves": 0}


    @classmethod
    def from_saved(cls, coefficients, p, background, path, model=None):
        """Restore a verified full eigensystem without a second eigensolve.

        The caller validates the artifact hash/provenance. This method also
        checks every restored matrix/vector dimension, the original physical
        basis, the M-H operator, continuous forcing and spectral identities.
        It restores every coordinate; no selection or modal filtering occurs.
        """
        if not isinstance(p, (int, np.integer)) or p < 2:
            raise ValueError("Maximum degree p must be an integer >=2")
        result = cls.__new__(cls)
        result.coefficients, result.p, result.background = coefficients, int(p), background
        result.length, result.n, result.ndof = background.length, int(p)-1, 2*(int(p)-1)
        result.matrix_nq = result.nq = 2*int(p)+1
        result.model = rod.derive_polynomials() if model is None else model
        result.derivation = derive_second_order(result.model)
        result.disc = planar.PlanarGalerkin(coefficients, int(p), length=result.length,
                                            nq=result.matrix_nq, model=result.model)
        result.slices = {"u": slice(0, result.n), "c": slice(result.n, result.ndof)}
        result.transforms = tuple(result.disc._transforms[index] for index in (0, 3))
        result.B = tuple(result.disc.B[index] for index in (0, 3))
        result.D = tuple(result.disc.D[index] for index in (0, 3))
        result.points, result.weights = result.disc.x, result.disc.weights
        required = ("M", "K", "f0", "f2", "omega", "vectors", "b0", "b2",
                    "transform_u", "transform_c")
        with np.load(Path(path), allow_pickle=False) as saved:
            if any(name not in saved for name in required):
                raise ValueError("Saved response is missing a full spectral artifact")
            arrays = {name: np.asarray(saved[name]).copy() for name in required}
        shapes = {"M": (result.ndof, result.ndof), "K": (result.ndof, result.ndof),
                  "vectors": (result.ndof, result.ndof),
                  "transform_u": (result.n, result.n), "transform_c": (result.n, result.n)}
        for name, array in arrays.items():
            if array.shape != shapes.get(name, (result.ndof,)) or not np.all(np.isfinite(array)):
                raise ValueError(f"Saved response shape/finiteness mismatch: {name}")
        for name in ("M", "K", "f0", "f2", "omega", "vectors", "b0", "b2"):
            setattr(result, name, arrays[name])
        if np.any(result.omega <= 0):
            raise ValueError("Saved fixed-fixed response has nonpositive frequencies")
        Bu, Bc = result.B
        Du, Dc = result.D
        weights, c = result.weights, coefficients
        expected_M, expected_K = np.zeros_like(result.M), np.zeros_like(result.K)
        us, cs = result.slices["u"], result.slices["c"]
        expected_M[us, us] = c.m*Bu.T@(weights[:, None]*Bu)
        expected_M[cs, cs] = c.jp*Bc.T@(weights[:, None]*Bc)
        expected_K[us, us] = c.C*Du.T@(weights[:, None]*Du)
        expected_K[us, cs] = c.nu*c.C*Du.T@(weights[:, None]*Bc)
        expected_K[cs, us] = expected_K[us, cs].T
        expected_K[cs, cs] = c.C*Bc.T@(weights[:, None]*Bc)+c.H*Dc.T@(weights[:, None]*Dc)
        expected_f0, expected_f2 = result.assemble_forcing(result.nq)
        checks = {}
        # The same 2e-11 operator/eigen audit tolerance as the diagnostic
        # contract; not a changed trajectory-convergence acceptance gate.
        tolerance = 2e-11
        for name, actual, expected in (
                ("mass", result.M, expected_M), ("stiffness", result.K, expected_K),
                ("transform_u", arrays["transform_u"], result.transforms[0]),
                ("transform_c", arrays["transform_c"], result.transforms[1]),
                ("forcing0", result.f0, expected_f0), ("forcing2", result.f2, expected_f2),
                ("modal_force0", result.b0, result.vectors.T@result.f0),
                ("modal_force2", result.b2, result.vectors.T@result.f2)):
            relative = float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected), 1e-30))
            checks[name+"_relative_difference"] = relative
            if relative > tolerance:
                raise ValueError(f"Saved response differs from declared model/background: {name}")
        result.eigenvalues = result.omega**2
        result.driving_omega = 2*background.omega
        result.eigen_decompositions, result.exact_time_evaluations = 0, 0
        residual = result.K@result.vectors-result.M@result.vectors*result.eigenvalues[None, :]
        orthogonality = float(np.max(abs(result.vectors.T@result.M@result.vectors-np.eye(result.ndof))))
        relative_residual = float(np.linalg.norm(residual)/(np.linalg.norm(result.K@result.vectors)
                                  +np.linalg.norm(result.M@result.vectors*result.eigenvalues[None, :])))
        if orthogonality > tolerance or relative_residual > tolerance:
            raise ValueError("Saved response spectral orthogonality/residual check failed")
        np.linalg.cholesky(result.M)
        np.linalg.cholesky(result.K)
        indices = np.r_[result.disc._indices["u"], result.disc._indices["c"]]
        old_M = result.disc.M0[np.ix_(indices, indices)]
        old_K = result.disc.K[np.ix_(indices, indices)]
        result.matrix_checks = {
            "mass_max_absolute_difference": float(np.max(abs(result.M-old_M))),
            "mass_relative_difference": float(np.linalg.norm(result.M-old_M)/np.linalg.norm(old_M)),
            "stiffness_max_absolute_difference": float(np.max(abs(result.K-old_K))),
            "stiffness_relative_difference": float(np.linalg.norm(result.K-old_K)/np.linalg.norm(old_K)),
            "mass_relative_symmetry": float(np.linalg.norm(result.M-result.M.T)/np.linalg.norm(result.M)),
            "stiffness_relative_symmetry": float(np.linalg.norm(result.K-result.K.T)/np.linalg.norm(result.K)),
        }
        result.eigen_checks = {
            "coordinates_retained": result.ndof, "eigenvectors_retained": result.vectors.shape[1],
            "mass_orthogonality_max_absolute": orthogonality,
            "eigenpair_relative_residual": relative_residual,
            "omega_min": float(result.omega[0]), "omega_max": float(result.omega[-1]),
            "minimum_absolute_detuning": float(np.min(abs(result.omega-result.driving_omega))),
        }
        result.restore_checks = {"source": str(path), "eigendecompositions": 0,
                                 "tolerance": tolerance, **checks}
        return result

