"""D17 production entry for targeted EB KV modes, with explicit routing.

The historical K12 Provider remains intact. New full solves use direct expm
for T and Frechet only for its derivative. Exact identical arms use the
single verified K17 analytic half assembly; no duplicated beam/joint law.
"""
from dataclasses import dataclass
import numpy as np
from scipy.linalg import expm, expm_frechet
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as symmetry

ROUTING_VERSION = 'eb-kv-routing-v1'
ROOT_AGREEMENT = dict(relative=1e-9, absolute_near_zero=1e-10)


@dataclass(frozen=True)
class Config:
    """Clamped EB arms and the one massless rotational KV joint only.

    mu records an explicit geometric construction, when supplied. A nonzero
    mu can never select a reduced path, even if tiny length changes round off.
    Other supports/arm laws are outside this API, not silently approximated.
    """
    arms: tuple
    beta_rad: float
    kappa_theta: float
    d_theta: float
    mu: float | None = None

    def __post_init__(self):
        object.__setattr__(self, 'arms', tuple(self.arms))
        if len(self.arms) != 2 or any(a.model != 'EB' or a.invS or a.J for a in self.arms):
            raise ValueError('routing v1 supports two clamped classical EB arms only')
        kv.joint_matrix(0j, self.beta_rad, self.kappa_theta*kv.M_REF,
                        self.d_theta*kv.M_REF*kv.T_REF)
        if self.mu is not None and (not np.isfinite(self.mu) or abs(self.mu) >= 1):
            raise ValueError('finite |mu| < 1 required')

    @property
    def identical(self):
        return self.arms[0] == self.arms[1] and (self.mu is None or self.mu == 0)


def route(config, solver_path='auto'):
    if solver_path not in ('auto', 'reduced', 'full'):
        raise ValueError('solver_path must be auto, reduced or full')
    if solver_path == 'reduced' and not config.identical:
        raise ValueError('forced reduced requires structurally identical arms and mu=0')
    return ('SYMMETRY_REDUCED' if solver_path == 'reduced' or
            (solver_path == 'auto' and config.identical) else 'FULL_TWO_ARM')


class FullProvider(kv.Provider):
    """Bounded consistency change: authoritative direct T, same analytic dH."""
    def __init__(self, config, calls):
        super().__init__(config.arms, config.beta_rad, config.kappa_theta, config.d_theta, calls)

    def transfer(self, z, arm, *, derivative=False, x=None):
        p = kv.spectral(z)/kv.T_REF
        x = arm.L if x is None else kv.real(x, 'x')
        scale = arm.scale()
        X = kv.state_matrix(p, arm)*scale[None, :]/scale[:, None]*x
        self.calls.reserve(3 if derivative else 1)
        self.calls.expm += 1
        T = expm(X)
        physical = scale[:, None]*T/scale[None, :]
        if not derivative:
            return physical
        Xz = kv.state_matrix(p, arm, derivative=True)*scale[None, :]/scale[:, None]*x/kv.T_REF
        self.calls.frechet += 1
        Tz = expm_frechet(X, Xz, compute_expm=False)
        return physical, scale[:, None]*Tz/scale[None, :]

    def matrices(self, z, *, derivative=False):
        self.calls.reserve(8 if derivative else 3)
        self.calls.full_B += 1
        self.calls.full_B_z += int(derivative)
        return super().matrices(z, derivative=derivative)


def agreement(value, reference):
    absolute = float(abs(value-reference))
    scale = float(abs(reference))
    relative = absolute/scale if scale else None
    # Near zero uses the declared absolute criterion; elsewhere relative only.
    accepted = (absolute <= ROOT_AGREEMENT['absolute_near_zero'] if scale < 1e-6
                else relative <= ROOT_AGREEMENT['relative'])
    return dict(absolute=absolute, relative=relative, accepted=bool(accepted))


def solve_mode(config, predictor, *, eta=None, solver_path='auto', seed_states=None,
               elastic_z=None, calls=None):
    """One targeted correction, no retry/search/continuation orchestration.

    Both routes return normalized *two-arm physical* states and reactions.
    An eta=-1 identical-arm solution is validated at an explicit elastic_z,
    without any complex Newton call. sorted labels are external metadata.
    """
    path = route(config, solver_path)
    if path == 'SYMMETRY_REDUCED' and eta not in (-1, 1):
        raise ValueError('reduced mode requires explicit eta=+1 or -1')
    if eta is not None and eta not in (-1, 1):
        raise ValueError('eta must be +1, -1 or None')
    calls = calls if calls is not None else symmetry.Calls(budget=2000)
    full = FullProvider(config, calls)
    inactive = bool(config.identical and eta == -1)
    z = kv.spectral(predictor)
    if inactive:
        if elastic_z is None or complex(elastic_z).real != 0 or complex(elastic_z).imag <= 0:
            raise ValueError('exact inactive reuse requires a positive pure-imaginary elastic eigenvalue')
        z = complex(elastic_z)  # reuse the elastic eigenvalue, never clip a solved complex root
    matrix = symmetry.AnalyticHalfProvider(full, eta) if path == 'SYMMETRY_REDUCED' else full
    frozen = symmetry.FrozenBalanced(matrix, z)
    B0, _ = frozen.matrices(z)
    initial = kv.right_null(B0)
    if inactive:
        correction = dict(z=z, a=initial, steps=0, history=[], last_delta_z=None,
                          status='EXACT_INACTIVE_BY_SYMMETRY')
    else:
        correction = kv.correct(frozen.matrices, z, initial)
        calls.corrections += correction['steps']
    z = correction['z']
    reactions_hat = frozen.reactions(correction['a'])
    if path == 'SYMMETRY_REDUCED':
        shape = symmetry.recover_closed(matrix, z, reactions_hat)
    else:
        calls.reserve(2)
        shape = kv.recover(full, z, reactions_hat)
    MAC = None
    if seed_states is not None:
        seed_states = np.asarray(seed_states)
        if seed_states.shape != shape['states'].shape:
            raise ValueError('seed and result require the same material grid and two-arm state ordering')
        _, weights = kv.quadrature(seed_states.shape[1])
        seed_vector = kv.mass_vector(seed_states, config.arms, weights)
        overlap = np.vdot(seed_vector, shape['vector'])
        if abs(overlap):
            phase = overlap.conjugate()/abs(overlap)
            for name in ('states','a','reactions','vector'):
                shape[name] *= phase
        MAC = float(kv.mac_matrix([seed_vector], [shape['vector']])[0,0])
    diag = kv.diagnose(full, z, shape)
    full_diag = diag.copy()
    primary_B, primary_Bz = matrix.matrices(z, derivative=True)
    primary_a = shape['a'][:3] if path == 'SYMMETRY_REDUCED' else shape['a']
    primary = symmetry.spectrum(primary_B, primary_Bz, primary_a)
    balanced_B, balanced_Bz = frozen.matrices(z, derivative=True)
    balanced = symmetry.spectrum(balanced_B, balanced_Bz)
    reduced_physical = None
    if path == 'SYMMETRY_REDUCED':
        diag = dict(diag, null_residual=primary['null_residual'],
            sigma_ratio=primary['ratios'][-1], next_sigma_ratio=primary['ratios'][-2])
        reduced_physical = symmetry.physical_details(matrix, z, shape)
    gates = kv.failures(diag, z, 'INACTIVE' if inactive else 'ACTIVE',
                        complex(elastic_z).imag if elastic_z is not None else z.imag,
                        MAC if MAC is not None else 1.)
    if config.identical and eta is not None and diag['symmetry_class'] != eta:
        gates.append('SYMMETRY_CLASS_MISMATCH')
    if not config.identical:
        # Reflection is not a symmetry of unequal arms: its defect is diagnostic,
        # not an admissibility condition. No projection is performed.
        gates = [g for g in gates if g != 'SYMMETRY_DEFECT']
    if reduced_physical and max(reduced_physical['half_normalized']) > kv.CRITERIA['physical_residual']:
        gates.append('REDUCED_PHYSICAL_GATE')
    if not inactive and correction['status'] != 'CONVERGED':
        gates.append(correction['status'])
    root_failures = [g for g in gates if g in ('NULL_RESIDUAL','SIGMA_RATIO','CONJUGATE_RESIDUAL',
        'NONOSCILLATORY_OR_AXIS_APPROACH','NEWTON_LIMIT')]
    form_failures = [g for g in gates if g not in root_failures and g != 'POSSIBLE_MULTIPLICITY']
    return dict(z=z, p=z/kv.T_REF, shape=shape, solver_path=path,
        eta=eta if config.identical else None, symmetry_reduced=path=='SYMMETRY_REDUCED',
        symmetry_status='EXACT_IDENTICAL' if config.identical else 'NON_IDENTICAL',
        activity_status='EXACT_INACTIVE_BY_SYMMETRY' if inactive else 'ACTIVE_OR_UNCLASSIFIED',
        complex_newton_calls=0 if inactive else 1, correction=correction,
        diagnostics=diag, full_diagnostics=full_diag, primary_spectrum=primary,
        frozen_spectrum=balanced, reduced_physical=reduced_physical, MAC=MAC,
        accepted=not gates, failures=gates, root_failures=root_failures, form_failures=form_failures,
        root_equation_status='PASS' if not root_failures else 'FAIL',
        form_recovery_status='PASS' if not form_failures else 'FAIL',
        rank_status='QUALIFIED' if 'POSSIBLE_MULTIPLICITY' in gates else 'PASS',
        frozen_row_scales=frozen.rows, frozen_column_scales=frozen.cols,
        symmetry_gate_applicable=config.identical, calls=calls.snapshot())
