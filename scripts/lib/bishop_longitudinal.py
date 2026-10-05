"""Diagnostic-only local longitudinal rod; see docs/theory/bishop_literature_reproduction.md.

No import of, or changes to, the project's bending/angled-joint solvers.
State order is (U, U', N, P); Gamma in Marais is -N. H=0 is second order.
"""
from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class Segment:
    L: float
    EA: float
    m: float
    H: float = 0.0
    J: float = 0.0

    def __post_init__(self):
        if not all(np.isfinite(x) for x in (self.L, self.EA, self.m, self.H, self.J)):
            raise ValueError("Nonfinite segment")
        if min(self.L, self.EA, self.m) <= 0 or min(self.H, self.J) < 0:
            raise ValueError("Require L, EA, m > 0 and H, J >= 0")


def circular_segment(L, E, rho, nu, radius, model="bishop"):
    if model not in ("wave", "rayleigh", "bishop") or not (0 <= nu < .5) or radius <= 0:
        raise ValueError("Unsupported circular isotropic parameters/model")
    area, polar = np.pi * radius**2, np.pi * radius**4 / 2
    j = nu**2 * rho * polar if model != "wave" else 0.0
    h = nu**2 * E / (2 * (1 + nu)) * polar if model == "bishop" else 0.0
    return Segment(L, E * area, rho * area, h, j)


def speed_segment(L, c, nu, diameter, model="bishop"):
    """Equation divided by rho*A: m=1 is a normalization, not a material density."""
    if model not in ("wave", "rayleigh", "bishop") or not (0 <= nu < .5):
        raise ValueError("Unsupported model/Poisson ratio")
    if c <= 0 or diameter <= 0:
        raise ValueError("Positive c and diameter required")
    j = nu**2 * diameter**2 / 8 if model != "wave" else 0.0
    h = c**2 * j / (2 * (1 + nu)) if model == "bishop" else 0.0
    return Segment(L, c**2, 1.0, h, j)


def wave_numbers(s, omega):
    """a>0 and b>0: roots +a,-a,+ib,-ib; rationalization avoids cancellation."""
    if s.H <= 0:
        raise ValueError("Fourth-order wave numbers require H>0")
    q = (s.EA - s.J * omega**2) / (2 * s.H)
    t = s.m * omega**2 / s.H
    radical = np.sqrt(q * q + t + 0j)
    if np.real(q) >= 0:
        a2 = radical + q
        b2 = t / a2
    else:
        b2 = radical - q
        a2 = t / b2
    a, b = np.sqrt(a2), np.sqrt(b2)
    if not np.iscomplexobj(omega):
        return float(a.real), float(b.real)
    return a, b


def basis(s, omega, x, derivative=0):
    """Bounded basis derivatives, columns are coefficients; local x in [0,L]."""
    x = np.atleast_1d(x)
    if derivative < 0 or derivative > 4:
        raise ValueError("Derivative must be 0..4")
    if s.H == 0:
        d = s.EA - s.J * omega**2
        if not np.iscomplexobj(omega) and d <= 0:
            raise ValueError("Rayleigh-Love cutoff: no oscillatory basis at/above cutoff")
        b = np.sqrt(s.m * omega**2 / d)
    else:
        a, b = wave_numbers(s, omega)
    if omega == 0:
        pair = (np.ones_like(x), x) if derivative == 0 else (
            np.zeros_like(x), np.ones_like(x) if derivative == 1 else np.zeros_like(x))
    else:
        # Exact sign cycle avoids a spurious cos(pi/2) at a free endpoint.
        cosine, sine = np.cos(b*x), np.sin(b*x)
        cycle = [(cosine, sine), (-sine, cosine), (-cosine, -sine), (sine, -cosine)]
        pair = tuple(b**derivative * value for value in cycle[derivative % 4])
    columns = list(pair)
    if s.H > 0:
        columns += [(-a)**derivative * np.exp(-a*x),
                    a**derivative * np.exp(-a*(s.L-x))]
    return np.stack(columns, axis=-1)


def state_basis(s, omega, x):
    u, du, ddu, dddu = [basis(s, omega, [x], d)[0] for d in range(4)]
    return np.array([u, du, (s.EA-s.J*omega**2)*du-s.H*dddu, s.H*ddu])


def boundary_rows(kind, fourth_order):
    if fourth_order:
        return {"C": (0, 1), "F": (2, 3), "UP": (0, 3)}[kind]
    if kind == "C":
        raise ValueError("C=U=U'=0 belongs to fourth order; use U for H=0")
    return {"U": (0,), "F": (2,)}[kind]


def boundary_matrix(segments, f_hz, ends, balanced=True):
    fourth = segments[0].H > 0
    if any((s.H > 0) != fourth for s in segments):
        raise ValueError("Mixed-order interfaces are outside this diagnostic contract")
    width = 4 if fourth else 2
    omega = 2*np.pi*f_hz
    length, ea = sum(s.L for s in segments), max(s.EA for s in segments)
    state_scale = np.array([1., 1/length, ea/length, ea])
    dtype = complex if np.iscomplexobj(f_hz) else float
    matrix = np.zeros((width*len(segments), width*len(segments)), dtype=dtype)
    row = 0
    for end, index, x in ((ends[0], 0, 0), (ends[1], len(segments)-1, segments[-1].L)):
        rows = boundary_rows(end, fourth)
        b = state_basis(segments[index], omega, x)/state_scale[:, None]
        matrix[row:row+len(rows), index*width:(index+1)*width] = b[list(rows)]
        row += len(rows)
    continuity = [0, 1, 2, 3] if fourth else [0, 2]
    for i in range(len(segments)-1):
        left = state_basis(segments[i], omega, segments[i].L)/state_scale[:, None]
        right = state_basis(segments[i+1], omega, 0)/state_scale[:, None]
        matrix[row:row+width, i*width:(i+1)*width] = left[continuity]
        matrix[row:row+width, (i+1)*width:(i+2)*width] = -right[continuity]
        row += width
    columns = np.ones(len(matrix))
    if balanced:
        matrix = matrix/np.maximum(np.linalg.norm(matrix, axis=1), 1e-300)[:, None]
        # A column can physically vanish at a free-free second-order eigenvalue.
        # Never amplify that null column (and its roundoff) to unit magnitude.
        columns = 1/np.maximum(np.linalg.norm(matrix, axis=0), 1.)
        matrix = matrix*columns
    return matrix, columns


def characteristic(segments, f_hz, ends):
    return np.linalg.det(boundary_matrix(segments, f_hz, ends)[0])


def bounded_roots(function, bounds, intervals, expected, contract):
    """One scan and at most one doubled scan. Every bracket/error is retained."""
    attempts = []
    for attempt in range(2):
        calls = 0
        def counted(f):
            nonlocal calls
            calls += 1
            value = float(function(f))
            if not np.isfinite(value):
                raise ValueError("Nonfinite characteristic")
            return value
        xs = np.linspace(*bounds, intervals*2**attempt+1)
        record = {"intervals": len(xs)-1, "bounds_hz": list(bounds), "brackets": [],
                  "errors": [], "reason": "initial" if attempt == 0 else "count/error gate: one doubled scan"}
        ys = []
        for x in xs:
            try:
                ys.append(counted(x))
            except (ValueError, FloatingPointError) as error:
                ys.append(np.nan)
                record["errors"].append({"frequency_hz": float(x), "error": str(error)})
        roots = []
        for a, b, fa, fb in zip(xs[:-1], xs[1:], ys[:-1], ys[1:]):
            if not np.isfinite(fa+fb) or fa*fb > 0:
                continue
            item = {"a_hz": float(a), "b_hz": float(b)}
            before = calls
            try:
                root = brentq(counted, a, b, xtol=contract["root_xtol_hz"],
                              rtol=contract["root_rtol"], maxiter=contract["root_maxiter"])
                if not roots or abs(root-roots[-1]) > 1e-5:
                    roots.append(root)
                item.update(root_hz=root, status="CONVERGED")
            except (ValueError, RuntimeError) as error:
                item.update(status="UNRESOLVED", error=str(error))
                record["errors"].append(item.copy())
            item["evaluations"] = calls-before
            record["brackets"].append(item)
        record.update(evaluations=calls, roots_hz=roots)
        attempts.append(record)
        if len(roots) == expected and not record["errors"]:
            return np.array(roots), {"status": "PASS", "attempts": attempts}
    return np.array(roots), {"status": "UNRESOLVED", "attempts": attempts}


def explicit_frequencies(s, count):
    k = np.arange(1, count+1)*np.pi/s.L
    return np.sqrt((s.EA*k*k+s.H*k**4)/(s.m+s.J*k*k))/(2*np.pi)


def popov_characteristic(s, f_hz, end):
    """Direct printed (13), scaled by positive a^(4j+2), no large cosh."""
    if end not in ("C", "F"):
        raise ValueError("Equation (13) only covers C-C and F-F")
    a, b = wave_numbers(s, 2*np.pi*f_hz)
    j = 0 if end == "C" else 1
    r = b/a
    sech = 2*np.exp(-a*s.L)/(1+np.exp(-2*a*s.L))
    return (2*r**(2*j+1)*(np.cos(b*s.L)-sech)
            + (-1)**j*(r**(4*j+2)-1)*np.sin(b*s.L)*np.tanh(a*s.L))


@lru_cache(maxsize=8)
def quadrature(order):
    return np.polynomial.legendre.leggauss(order)


def mode_coefficients(segments, frequency, ends):
    b, columns = boundary_matrix(segments, frequency, ends)
    _, singular, vh = np.linalg.svd(b)
    coeff = columns*vh[-1]
    width = 4 if segments[0].H > 0 else 2
    coeff = coeff.reshape(len(segments), width)
    # Orient by the first nonzero interior displacement, as in Marais Fig. 2.
    near_left = basis(segments[0], 2*np.pi*frequency, [segments[0].L*.01])[0]@coeff[0]
    coeff *= 1 if near_left >= 0 else -1
    return coeff, {"singular_values": singular.tolist(),
                   "singular_ratio": float(singular[-1]/singular[0]),
                   "nonnull_condition": float(singular[0]/singular[-2]),
                   "balanced_bc_residual": float(np.linalg.norm(b@vh[-1], ord=np.inf))}


def verify_modes(segments, frequencies, ends, quadrature_order=256):
    """Modal profiles, all physical BCs, ODE, energy, mass Gram, both norms."""
    coefficients, diagnostics = zip(*(mode_coefficients(segments, f, ends) for f in frequencies))
    coefficients = np.array(coefficients)
    count = len(frequencies)
    mass, energy = np.zeros((count, count)), np.zeros(count)
    ode = np.zeros(count)
    maxima = np.zeros(count)
    z, weights = quadrature(quadrature_order)
    for i, s in enumerate(segments):
        # Split off physical boundary layers; no missed quadrature mass at long aL.
        if s.H > 0:
            amin = min(wave_numbers(s, 2*np.pi*f)[0] for f in frequencies)
            layer = min(s.L/2, 20/amin)
            cuts = sorted(set([0., layer, s.L-layer, s.L]))
        else:
            cuts = [0., s.L]
        x = np.concatenate([(left+right)/2+(right-left)/2*z for left, right in zip(cuts[:-1], cuts[1:])])
        w = np.concatenate([(right-left)/2*weights for left, right in zip(cuts[:-1], cuts[1:])])
        fields = np.array([[basis(s, 2*np.pi*f, x, d)@coefficients[n, i]
                            for d in range(5)] for n, f in enumerate(frequencies)])
        mass += (fields[:, 0]*w)@fields[:, 0].T*s.m + (fields[:, 1]*w)@fields[:, 1].T*s.J
        energy += np.sum(w*(s.EA*fields[:, 1]**2+s.H*fields[:, 2]**2), axis=1)
        for n, f in enumerate(frequencies):
            omega = 2*np.pi*f
            terms = np.array([s.H*fields[n, 4], (s.J*omega**2-s.EA)*fields[n, 2],
                              -s.m*omega**2*fields[n, 0]])
            ode[n] = max(ode[n], np.max(np.abs(terms.sum(axis=0)))/np.max(np.abs(terms).sum(axis=0)))
            grid = np.linspace(0, s.L, 401)
            slopes = basis(s, omega, grid, 1)@coefficients[n, i]
            extrema = [0., s.L]
            for left, right, vl, vr in zip(grid[:-1], grid[1:], slopes[:-1], slopes[1:]):
                if vl*vr < 0:
                    extrema.append(brentq(lambda y: float(basis(s, omega, [y], 1)[0]@coefficients[n, i]), left, right))
            maxima[n] = max(maxima[n], np.max(np.abs(basis(s, omega, extrema)@coefficients[n, i])))
    norms = np.sqrt(np.diag(mass))
    gram = mass/np.outer(norms, norms)
    length, ea = sum(s.L for s in segments), max(s.EA for s in segments)
    profiles = []
    for n, (f, diag) in enumerate(zip(frequencies, diagnostics)):
        physical_b, _ = boundary_matrix(segments, f, ends, balanced=False)
        # Each physical row is in U, L*chi, L*N/EAref, P/EAref units.
        residuals = physical_b@coefficients[n].ravel()/maxima[n]
        diag.update(frequency_hz=float(f), ode_relative=float(ode[n]),
                    physical_bc_rows=residuals.tolist(), physical_bc_max=float(np.max(np.abs(residuals))),
                    energy_omega2=float(energy[n]/mass[n, n]),
                    energy_relative=float(abs(energy[n]/mass[n, n]/(2*np.pi*f)**2-1)),
                    mass_norm_before=float(norms[n]), max_abs_before=float(maxima[n]))
        offset = 0.
        for i, s in enumerate(segments):
            x = np.linspace(0, s.L, 601)
            u = basis(s, 2*np.pi*f, x)@coefficients[n, i]
            du = basis(s, 2*np.pi*f, x, 1)@coefficients[n, i]
            for xx, uu, slope in zip(x, u, du):
                profiles.append({"mode": n+1, "segment": i+1, "x_m": float(offset+xx),
                                 "Y": float(uu/maxima[n]), "U_mass": float(uu/norms[n]),
                                 "dU_mass_dx": float(slope/norms[n])})
            offset += s.L
    return {"modes": list(diagnostics), "mass_gram": gram.tolist(),
            "mass_orthogonality_max": float(np.max(np.abs(gram-np.eye(count)))),
            "mass_norm_squared_after": np.diag(gram).tolist(),
            "max_abs_after": (maxima/maxima).tolist(),
            "quadrature_order_per_subinterval": quadrature_order,
            "coefficients_mass_normalized": (coefficients/norms[:, None, None]).tolist(),
            "bc_definition": "physical matrix rows / continuous max|U|; states scaled by (1,1/L,EAref/L,EAref)",
            "ode_definition": "max|HU4+(Jw2-EA)U2-mw2U| / max(sum of absolute terms)",
            "orthogonality_definition": "max|G-I|, Gmn=integral(mUmUn+JUm'Un')/(mass_norm_m*mass_norm_n)"}, profiles


def marais_argument_count(segments, bounds=(1., 31000.), height=25., horizontal_samples=2048):
    """Finite contour count, not an interval-arithmetic theorem or general WW code."""
    lo, hi = bounds
    corners = [lo-1j*height, hi-1j*height, hi+1j*height, lo+1j*height, lo-1j*height]
    vertical_samples = max(64, horizontal_samples//32)
    points = np.concatenate([np.linspace(a, b, horizontal_samples if i%2 == 0 else vertical_samples, endpoint=False)
                             for i, (a, b) in enumerate(zip(corners[:-1], corners[1:]))])
    signs = np.array([np.linalg.slogdet(boundary_matrix(segments, f, ("C", "F"))[0])[0] for f in points])
    steps = np.angle(np.roll(signs, -1)/signs)
    winding = float(steps.sum()/(2*np.pi))
    # Branch points of sqrt((EA-J*w²)^2+4*H*m*w²) must lie outside the box.
    branch_points = []
    for s in segments:
        w2 = np.roots([s.J**2, -2*s.EA*s.J+4*s.H*s.m, s.EA**2])
        for value in w2:
            f = np.sqrt(complex(value))/(2*np.pi)
            branch_points.extend([f, -f])
    branch_clear = all(not (lo <= f.real <= hi and abs(f.imag) <= height) for f in branch_points)
    return {"count": int(round(winding)), "winding": winding,
            "max_phase_step_rad": float(max(abs(steps))), "evaluations": len(points),
            "bounds_hz": list(bounds), "height_hz": height,
            "horizontal_samples": horizontal_samples, "vertical_samples": vertical_samples,
            "branch_points_outside_contour": branch_clear,
            "branch_points_f_hz": [[f.real, f.imag] for f in branch_points]}


def independent_marais_transfer(parameters, seeds, dps):
    """Independent first-order expm shooting at five targets only, exact decimal inputs."""
    import mpmath as mp
    with mp.workdps(dps):
        E, rho, nu = [mp.mpf(str(parameters[k])) for k in ("E", "rho", "nu")]
        G = E/(2*(1+nu))
        arms = []
        for length, radius in zip(parameters["lengths"], parameters["radii"]):
            r, length = mp.mpf(str(radius)), mp.mpf(str(length))
            area, polar = mp.pi*r*r, mp.pi*r**4/2
            arms.append((length, E*area, rho*area, nu*nu*G*polar, nu*nu*rho*polar))
        total_L = sum(a[0] for a in arms)
        ea_ref = max(a[1] for a in arms)
        scales = [1, 1/total_L, ea_ref/total_L, ea_ref]
        calls = 0
        def determinant(f):
            nonlocal calls
            calls += 1
            omega = 2*mp.pi*f
            transfer = mp.eye(4)
            for length, ea, mass, h, j in arms:
                a = mp.matrix([[0, 1, 0, 0], [0, 0, 0, 1/h],
                               [-mass*omega**2, 0, 0, 0], [0, ea-j*omega**2, -1, 0]])
                for i in range(4):
                    for k in range(4):
                        a[i, k] *= length*scales[k]/scales[i]
                transfer = mp.expm(a)*transfer
            return transfer[2, 2]*transfer[3, 3]-transfer[2, 3]*transfer[3, 2]
        rows = []
        for seed in seeds:
            before = calls
            f = mp.mpf(str(seed))
            def addressed(candidate):
                if abs(candidate-f) > mp.mpf('.1'):
                    raise ValueError("Addressed evaluation escaped fixed +/-0.1 Hz window")
                return determinant(candidate)
            try:
                root = mp.findroot(addressed, (f-mp.mpf('.05'), f+mp.mpf('.05')),
                                   maxsteps=30, tol=mp.mpf(10)**(-(dps-15)))
                if abs(root-f) > mp.mpf('.1'):
                    raise ValueError("Addressed refinement escaped target window")
                rows.append({"seed_hz": seed, "frequency_hz": mp.nstr(root, dps-5),
                             "evaluation_bounds_hz": [mp.nstr(f-mp.mpf('.1'), 20), mp.nstr(f+mp.mpf('.1'), 20)],
                             "relative_to_double": float(abs(root-f)/root), "status": "PASS",
                             "determinant_absolute": mp.nstr(abs(determinant(root)), 8),
                             "evaluations": calls-before})
            except (ValueError, ZeroDivisionError) as error:
                rows.append({"seed_hz": seed, "status": "UNRESOLVED", "error": str(error),
                             "evaluation_bounds_hz": [mp.nstr(f-mp.mpf('.1'), 20), mp.nstr(f+mp.mpf('.1'), 20)],
                             "evaluations": calls-before})
        return {"dps": dps, "method": "independent scaled (U,chi,N,P) matrix exponential",
                "evaluations": calls, "roots": rows}
