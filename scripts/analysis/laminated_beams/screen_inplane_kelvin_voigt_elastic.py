"""D14: four independent EB elastic spectra, no damped solver or tracking.

The new contract is a participation table, not complex continuation or a map.
Reuse K11 roots/forms, its real detector and its two-arm reconstruction.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict
from dataclasses import asdict
import hashlib
import io
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
import numpy as np
import scipy
from scripts.lib import inplane_spring_modes as mechanics
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import inplane_rotational_spring_eb_modes as modes
from scripts.analysis.laminated_beams import pilot_inplane_rotational_spring_eb as pilot
from scripts.analysis.laminated_beams.check_inplane_spring_robustness import intervals

SOURCE = ROOT / 'results/laminated_beams/inplane_spring_robustness'
K12 = ROOT / 'results/laminated_beams/inplane_kelvin_voigt_pilot'
OUTPUT = ROOT / 'results/laminated_beams/inplane_kelvin_voigt_elastic_screening'
BETAS = (0., 5., 45., 75.)
REF, T_REF = mechanics.REFERENCE, float(mechanics.FREQUENCY_SCALE)
CRITERIA = dict(null_residual=1e-9, sigma_ratio=1e-9, physical_residual=1e-9,
                compatibility=1e-10, rank_rtol=1e-12, simple_sigma_separation=1e-8,
                symmetry_residual=1e-6, mass_rtol=1e-6, frequency_interpretation=1e-6,
                K12_slope_relative=1e-3, inactive_delta=1e-8,
                max_point_matrices=6000, max_root_retries=2,
                small_participation_descriptive=1e-8)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def atomic(path, text):
    path = Path(path)
    tmp = path.with_name(path.name + '.tmp')
    tmp.write_text(text, encoding='utf-8')
    os.replace(tmp, path)


def write_json(path, value):
    atomic(path, json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False))


def protected_hashes():
    return {p.relative_to(ROOT).as_posix(): sha(p) for directory in (SOURCE, K12)
            for p in sorted(directory.iterdir()) if p.suffix in ('.json', '.csv', '.npz')}


def configuration():
    integrated, p = mechanics.section()
    expected = dict(A=.011, D=2.979166666666667e-6, S=.0032051282051282055,
                    m=.01, J=2.083333333333334e-6)
    for key, value in expected.items():
        np.testing.assert_allclose(getattr(p, key), value, rtol=1e-12, atol=0.)
    assert np.linalg.norm(integrated.B) < 1e-12 and abs(integrated.I1) < 1e-12
    config = dict(model='EB', mu=0., L1=1., L2=1., b=.2, h=.05,
                  layup='H/L/L/H', contrast=.4, ply_angles_deg=[0.]*4,
                  ply_thickness=.0125, material_base=dict(E1=1.1, E2=.9, nu12=.3,
                  G12=1/2.6, G13=1/2.6, G23=1/2.6, rho=1.),
                  material_factors=[1.4, .6, .6, 1.4], K=5/6, properties=asdict(p),
                  B=integrated.B.tolist(), I1=float(integrated.I1),
                  kappa_theta=1., d_theta=0., k_theta=REF.D/REF.L,
                  D_ref=REF.D, m_ref=REF.m, l_ref=REF.L, t_ref=T_REF,
                  beta_deg=list(BETAS), nodes_per_arm=129)
    return p, config


def participation(delta, omega, mass):
    if not np.isfinite([abs(delta), omega, mass]).all() or omega <= 0 or mass <= 0:
        raise ValueError('finite positive frequency and modal mass required')
    s = REF.D / REF.L * abs(delta)**2 / (omega**2 * mass)
    Omega = omega * T_REF
    return dict(s_joint=float(s), a_slope_pred=float(.5 * Omega**2 * s),
                zeta_slope_pred=float(.5 * Omega * s))


def classify(eta, symmetry_residual, *, identical_arms, root_confirmed):
    if not root_confirmed or symmetry_residual > CRITERIA['symmetry_residual']:
        raise ValueError('unconfirmed root or symmetry')
    if identical_arms and eta == -1:
        return 'EXACT_INACTIVE_BY_SYMMETRY'
    return 'ACTIVE'


class ElasticMatrices:
    """Real EB only; counts both transfers actually used by boundary_assembly."""
    def __init__(self, beta, properties):
        if beta not in BETAS:
            raise ValueError('no automatic angular refinement')
        self.pair = mechanics.arms('EB', 0., properties)
        self.beta = math.radians(beta)
        self.joint = eb.Joint('SPRING', REF.D / REF.L)
        self.cache = OrderedDict()
        self.full = self.blocks = self.transfer_expm = 0

    def tick(self):
        if self.full + self.blocks >= CRITERIA['max_point_matrices']:
            raise RuntimeError('COST_LIMIT')

    def assembly(self, omega):
        if np.iscomplexobj(omega) or not np.isfinite(omega) or omega <= 0:
            raise ValueError('positive real elastic frequency required')
        if omega not in self.cache:
            self.tick()
            self.full += 1
            self.transfer_expm += 2
            self.cache[omega] = mechanics.assembly(omega, self.pair, self.beta, self.joint)
            if len(self.cache) > 128:
                self.cache.popitem(last=False)
        return self.cache[omega]

    def matrix(self, omega, eta):
        a = self.assembly(omega)
        self.tick()
        self.blocks += 1
        return modes.class_matrix(a.endpoint_map, self.beta, self.joint.k_theta, REF, eta)


def diagnose_shape(context, omega, states, reactions):
    """Check saved/recovered real forms without projecting or changing them."""
    if np.iscomplexobj(states) or np.iscomplexobj(reactions):
        raise ValueError('real-only source required; imaginary components may not be discarded')
    _, weights = modes.quadrature(states.shape[1])
    vector = mechanics.physical_vector(states, context.pair, weights)
    mass = float(np.vdot(vector, vector).real)
    units = np.array([REF.L, REF.L, 1., REF.D/REF.L**2, REF.D/REF.L**2, REF.D/REF.L])
    y, reflected = states/units, modes.reflect(states)/units
    eta = 1 if np.vdot(y, reflected).real >= 0 else -1
    defect = float(np.linalg.norm(reflected-eta*y)/np.linalg.norm(y))
    assembly = context.assembly(omega)
    hat = reactions.ravel()/assembly.reaction_scales
    root_residual = float(np.linalg.norm(assembly.dimensionless @ hat) /
                          (np.linalg.norm(assembly.dimensionless)*np.linalg.norm(hat)))
    eq = eb.positively_equilibrate_matrix(assembly.dimensionless)
    singular = np.linalg.svd(eq.scaled_matrix, compute_uv=False)
    sigma, next_sigma = float(singular[-1]/singular[0]), float(singular[-2]/singular[0])
    nullity = int(np.count_nonzero(singular/singular[0] <= CRITERIA['rank_rtol']))
    ends = states[:, -1, :].ravel()
    amplitude = np.max(abs(ends/np.tile(units, 2)))
    residuals = abs(eb.scalar_joint_residuals(ends/amplitude, context.beta, context.joint)/assembly.row_units)
    failures = []
    for name, value in [('null_residual', root_residual), ('sigma_ratio', sigma),
                        ('symmetry_residual', defect), ('physical_residual', max(residuals)),
                        ('compatibility', max(residuals[:2]))]:
        if value > CRITERIA[name]:
            failures.append(name.upper())
    if nullity != 1 or next_sigma < CRITERIA['simple_sigma_separation']:
        failures.append('NUMERICAL_UNRESOLVED_MULTIPLICITY')
    if abs(mass-1) > CRITERIA['mass_rtol']:
        failures.append('MASS_NORMALIZATION')
    return dict(mass_M=mass, symmetry_eta=eta, symmetry_residual=defect,
                root_residual=root_residual, sigma_ratio=sigma, next_sigma_ratio=next_sigma,
                detected_nullity=nullity, physical_residual=float(max(residuals)),
                physical_residuals=residuals.tolist(), failures=failures)


def validate_sources(config):
    old, damped = read_json(SOURCE/'run_manifest.json'), read_json(K12/'run_manifest.json')
    for saved in (old['preflight'], damped['configuration']):
        assert saved['layup'] == config['layup'] and saved['contrast'] == config['contrast']
        for key in ('A', 'D', 'S', 'm', 'J', 'K', 'width'):
            np.testing.assert_allclose(saved['properties'][key], config['properties'][key], rtol=1e-12, atol=0.)
    assert old['preflight']['ply_angles_deg'] == [0]*4
    assert old['preflight']['ply_thickness'] == .0125
    for key in ('b', 'h', 'D_ref', 'm_ref', 'l_ref'):
        assert old['parameters'][key] == config[key]
    assert old['parameters']['frequency_scale'] == T_REF
    for key in ('L1', 'L2', 'b', 'h', 'kappa_theta', 'D_ref', 'm_ref', 't_ref'):
        assert damped['configuration'][key] == config[key]
    return old['created_HEAD'], damped['initial_HEAD']


def saved_pool(beta, source):
    point = source['points'].get(f'EB_m0_k1_b{beta:g}')
    if point is None:
        return [], None
    assert point['model'] == 'EB' and point['mu'] == 0 and point['kappa'] == 1
    assert set(point['errors']) <= {'MISSING_TARGET_OR_GUARD', 'GUARD_QUALIFIED'}
    rows = point['roots']
    assert all(r['current_sorted_position'] == i for i, r in enumerate(rows, 1))
    assert all(a['Omega'] < b['Omega'] for a, b in zip(rows, rows[1:]))
    for row in rows:
        assert row['root_status'] == 'CONFIRMED' and not row['failures']
        assert row['multiplicity'] == row['detected_nullity'] == 1
        assert abs(row['Lambda']**2-row['Omega']) < 1e-12*row['Omega']
        assert abs(row['omega']*T_REF-row['Omega']) < 1e-12*row['Omega']
    # This partial K11 point has only two roots, with no rejected lower event.
    if point['status'] != 'CONFIRMED':
        assert beta == 45. and len(rows) == 2
        assert all(c['accepted'] for scan in point['search'] for c in scan['candidates'])
    csv_rows = [r for r in read_csv(SOURCE/'verified_roots.csv') if r['point_id'] == point['point_id']]
    assert [float(r['Omega']) for r in csv_rows] == [r['Omega'] for r in rows]
    return [dict(r, reuse_status='REUSED_ROOT_AND_FORM') for r in rows], point


def search_missing(context, beta, pool, previous, source, audit):
    """One bounded extension per missing group, existing detector/refiner only."""
    near = min((p for p in source['points'].values() if p['model'] == 'EB' and p['mu'] == 0
                and p['kappa'] == 1 and p['status'] == 'CONFIRMED' and len(p['roots']) >= 7),
               key=lambda p: abs(p['beta_deg']-beta))
    guesses = [r['Omega'] for r in near['roots'][:7]]
    upper = max(guesses) * 1.02 + .1
    lower = previous['search_upper']-.1 if previous else .01
    search, suspects, additions = [], [], []
    audit.update(lower=lower, upper=upper, predictor_beta=near['beta_deg'],
                 source_lower_coverage=previous['search_upper'] if previous else None,
                 scans=search, newly_calculated_roots=0, retries=0)
    for eta in (1, -1):
        provider = lambda omega: context.matrix(omega, eta)
        candidates = []
        for lo, hi, count in intervals(guesses, upper):
            if hi <= lower:
                continue
            candidates.extend(pilot.roots._scan_candidates(provider, T_REF, pilot.policy(max(lo, lower), hi),
                case_id=f'EB_screen_b{beta:g}', builder_id=f'EB_class_{eta}',
                scan_id='MISSING_ONLY', points=count, phases=(0.,))[0])
        candidates, proof = pilot.reconcile_local_detections(candidates, provider)
        accepted, ambiguous = pilot.consolidate(candidates)
        suspects.extend(c for c in candidates if pilot.suspicious(c))
        if ambiguous:
            raise ValueError('NUMERICAL_UNRESOLVED_DETECTION_CLUSTER')
        search.append(dict(eta=eta, candidates=[pilot.candidate_record(c) for c in candidates], reconciliations=proof))
        for c in accepted:
            if c.diagnostics.detected_nullity != 1:
                raise ValueError('NUMERICAL_UNRESOLVED_MULTIPLICITY')
            additions.append(dict(Omega=c.omega_bar, omega=c.omega_bar/T_REF,
                Lambda=math.sqrt(c.omega_bar), symmetry_class=eta,
                reuse_status='NEW_ELASTIC_ROOT', source='CURRENT_REAL_EB_CHARACTERISTIC_MATRIX'))
        audit['newly_calculated_roots'] = len(additions)
    combined = sorted(pool+additions, key=lambda r: r['Omega'])
    if len(combined) < 7:
        raise ValueError('NUMERICAL_UNRESOLVED_MISSING_TARGET_OR_GUARD')
    target = combined[5]['Omega']
    unresolved = [pilot.candidate_record(c) for c in suspects if c.interval_left_bar <= target]
    if unresolved:
        raise ValueError('NUMERICAL_UNRESOLVED_BELOW_TARGET')
    audit.update(unresolved_below_target=unresolved,
                 guard_suspects=[pilot.candidate_record(c) for c in suspects if c.interval_left_bar > target])
    return combined


def k12_control(rows):
    source = read_csv(K12/'modal_results.csv')
    result = dict(relative_criterion=CRITERIA['K12_slope_relative'], source=str((K12/'modal_results.csv').relative_to(ROOT)))
    for role, position in [('ACTIVE', 1), ('INACTIVE', 2)]:
        current = next(r for r in rows if r['beta_deg'] == 5 and r['sorted_mode'] == f'sorted_{position:02d}')
        subset = [r for r in source if r['model'] == 'EB' and r['role'] == role]
        zero = next(r for r in subset if float(r['d_theta']) == 0)
        small = min((r for r in subset if float(r['d_theta']) > 0), key=lambda r: float(r['d_theta']))
        assert zero['status'] == small['status'] == 'CONFIRMED'
        assert float(zero['beta0_deg']) == float(small['beta0_deg']) == 5.
        assert float(zero['Omega0']) == current['Omega']
        assert float(small['kappa_theta']) == 1.
        observed = float(small['a_decay'])/float(small['d_theta'])
        entry = dict(source_key=small['key'], d_theta=float(small['d_theta']),
                     observed_a_over_d=observed, predicted=current['a_slope_pred'],
                     elastic_Delta_source=float(zero['Delta_psi_real']),
                     source_Delta_imag=float(zero['Delta_psi_imag']))
        if role == 'ACTIVE':
            error = abs(observed-current['a_slope_pred'])/current['a_slope_pred']
            entry.update(relative_difference=error, accepted=error <= CRITERIA['K12_slope_relative'])
        else:
            entry.update(max_abs_a=max(abs(float(r['a_decay'])) for r in subset),
                max_abs_Delta=max(abs(complex(float(r['Delta_psi_real']), float(r['Delta_psi_imag']))) for r in subset),
                accepted=current['activity_status'] == 'EXACT_INACTIVE_BY_SYMMETRY'
                and max(abs(float(r['a_decay'])) for r in subset) <= 1e-8)
        result[role] = entry
    return result


def compute(output=OUTPUT):
    output = Path(output).resolve()
    if any(output.is_relative_to(p.resolve()) for p in (SOURCE, K12)) or not output.is_relative_to(ROOT/'results'):
        raise ValueError('source data read-only; output must be a new results directory')
    checkpoint = output/'diagnostics.json'
    if checkpoint.exists():
        data = read_json(checkpoint)
        assert data['protected_sources'] == protected_hashes(), 'source data changed'
        assert data['criteria'] == CRITERIA
        if set(data['groups']) == {str(b) for b in BETAS} and data.get('output_hashes') and (output/'run_manifest.json').exists():
            assert all(sha(output/name) == digest for name, digest in data['output_hashes'].items())
            return dict(missing_only=True, new_root_groups=0, new_elastic_solver_calls=0,
                        form_recoveries=0, matrix_builds=0, new_positive_d_complex_roots=0)
    else:
        data = dict(groups={}, criteria=CRITERIA, protected_sources=protected_hashes(),
                    initial_HEAD=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
                    git_status_at_invocation=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
                    new_positive_d_complex_roots=0)
    started = time.perf_counter()
    properties, config = configuration()
    data['source_heads'] = validate_sources(config)
    data['configuration'] = config
    data['configuration_sha256'] = hashlib.sha256(json.dumps(config, sort_keys=True).encode()).hexdigest()
    output.mkdir(parents=True, exist_ok=True)
    source = read_json(SOURCE/'diagnostics.json')
    shapes = dict(np.load(output/'screening_shapes.npz')) if (output/'screening_shapes.npz').exists() else {}
    with np.load(SOURCE/'shapes.npz', allow_pickle=False) as saved_shapes:
        for beta in BETAS:
            if str(beta) in data['groups']:
                continue
            tick = time.perf_counter()
            context = ElasticMatrices(beta, properties)
            group = dict(beta_deg=beta, rows=[], status='NUMERICAL_UNRESOLVED', retries=0,
                         form_recoveries=0, form_analytic_calls=0, quadrature_checks=[])
            try:
                pool, previous = saved_pool(beta, source)
                group['source_point_status'] = previous['status'] if previous else None
                group['new_root_groups'] = int(len(pool) < 7)
                if len(pool) < 7:
                    group['search'] = {}
                    pool = search_missing(context, beta, pool, previous, source, group['search'])
                group['guard'] = {k: pool[6][k] for k in ('Omega', 'omega', 'reuse_status')}
                guard = eb.endpoint_diagnostics(context.assembly(pool[6]['omega']), context.beta, context.joint, REF)
                group['guard'].update(sigma_ratio=guard['sigma_ratio'], nullity=guard['nullity'])
                group['guard']['status'] = 'CONFIRMED' if guard['nullity'] == 1 and guard['sigma_ratio'] <= 1e-9 else 'GUARD_QUALIFIED'
                for position, record in enumerate(pool[:6], 1):
                    omega = record['omega']
                    key = f'b{beta:g}_sorted_{position:02d}'
                    if record['reuse_status'] == 'REUSED_ROOT_AND_FORM':
                        states = saved_shapes[record['shape_key']+'__states']
                        reactions = saved_shapes[record['shape_key']+'__reactions']
                        provenance = f"{SOURCE.relative_to(ROOT).as_posix()}/verified_roots.csv#{record['shape_key']}"
                    else:
                        recovered = mechanics.recover(context.assembly(omega), omega, context.pair,
                            context.beta, context.joint, parity=record['symmetry_class'])
                        states, reactions = recovered['states'], recovered['reactions']
                        group['form_recoveries'] += 1
                        group['form_analytic_calls'] += recovered['reconstruction_calls']['analytic']
                        provenance = 'NEW_REAL_EB_ROOT_AND_FORM_D14'
                    diagnostic = diagnose_shape(context, omega, states, reactions)
                    if diagnostic['failures']:
                        group.setdefault('rejected', []).append(dict(position=position, Omega=record['Omega'], diagnostics=diagnostic))
                        raise ValueError('NUMERICAL_UNRESOLVED_ROOT_OR_SHAPE')
                    assert diagnostic['symmetry_eta'] == record['symmetry_class']
                    delta = float(states[0, -1, 2]-states[1, -1, 2])
                    row = dict(model='EB', beta_deg=beta, kappa_theta=1., d_theta=0.,
                        sorted_mode=f'sorted_{position:02d}', source_provenance=provenance,
                        omega=omega, Omega=record['Omega'], Lambda=record['Lambda'],
                        psi1_joint=float(states[0, -1, 2]), psi2_joint=float(states[1, -1, 2]),
                        Delta_psi=delta, abs_Delta_psi=abs(delta),
                        activity_status=classify(diagnostic['symmetry_eta'], diagnostic['symmetry_residual'],
                                                identical_arms=context.pair[0] == context.pair[1], root_confirmed=True),
                        **participation(delta, omega, diagnostic['mass_M']),
                        **{k: v for k, v in diagnostic.items() if k not in ('failures', 'physical_residuals')},
                        root_status='CONFIRMED', reuse_status=record['reuse_status'], shape_key=key, notes='')
                    if row['activity_status'] == 'ACTIVE' and row['s_joint'] < CRITERIA['small_participation_descriptive']:
                        row['notes'] = 'SMALL_PARTICIPATION; not proof of exact inactivity'
                    shapes[key+'__states'], shapes[key+'__reactions'] = states, reactions
                    group['rows'].append(row)
                    group.setdefault('physical_conditions', {})[key] = diagnostic['physical_residuals']
                    if beta == 5 and position <= 2:
                        values = [mechanics.states_along_arm(omega, arm, r, nodes=257) for arm, r in zip(context.pair, reactions)]
                        fine = np.array([v[0] for v in values])
                        _, weights = modes.quadrature(257)
                        v = mechanics.physical_vector(fine, context.pair, weights)
                        fine_mass = float(np.vdot(v, v).real)
                        relative = abs(fine_mass-diagnostic['mass_M'])/diagnostic['mass_M']
                        group['quadrature_checks'].append(dict(sorted_mode=row['sorted_mode'], mass_257=fine_mass,
                            relative_mass_difference=relative, accepted=relative <= CRITERIA['mass_rtol']))
                        group['form_recoveries'] += 1
                        group['form_analytic_calls'] += 2
                group['status'] = 'CONFIRMED' if group['guard']['status'] == 'CONFIRMED' else 'TARGET_CONFIRMED_GUARD_QUALIFIED'
            except (ValueError, RuntimeError) as error:
                group['failure'] = str(error)
            group.update(full_B=context.full, symmetry_B=context.blocks, transfer_expm=context.transfer_expm,
                         seconds=time.perf_counter()-tick)
            data['groups'][str(beta)] = group
            with (output/'screening_shapes.npz.tmp').open('wb') as stream:
                np.savez_compressed(stream, **shapes)
            os.replace(output/'screening_shapes.npz.tmp', output/'screening_shapes.npz')
            write_json(checkpoint, data)
            print(beta, group['status'], 'states', len(group['rows']), 'B', context.full+context.blocks, flush=True)
    rows = [r for group in data['groups'].values() for r in group['rows']]
    if sum(r['beta_deg'] == 5 for r in rows) >= 2:
        data['K12_control'] = k12_control(rows)
    data['protected_sources_unchanged'] = data['protected_sources'] == protected_hashes()
    assert data['protected_sources_unchanged']
    stream = io.StringIO(newline='')
    if rows:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    atomic(output/'elastic_screening.csv', stream.getvalue())
    data['output_hashes'] = {name: sha(output/name) for name in ('elastic_screening.csv', 'screening_shapes.npz')}
    data['total_seconds'] = data.get('total_seconds', 0.) + time.perf_counter()-started
    write_json(checkpoint, data)
    counters = {key: sum(g.get(key, 0) for g in data['groups'].values()) for key in
                ('full_B', 'symmetry_B', 'transfer_expm', 'form_recoveries', 'form_analytic_calls', 'new_root_groups', 'retries')}
    counters.update(new_elastic_solver_calls=counters['new_root_groups'],
        new_elastic_modal_roots=sum(r['reuse_status'] == 'NEW_ELASTIC_ROOT' for r in rows),
        reused_modal_roots=sum(r['reuse_status'] == 'REUSED_ROOT_AND_FORM' for r in rows),
        all_new_roots_including_guard=sum(g.get('search', {}).get('newly_calculated_roots', 0) for g in data['groups'].values()),
        new_positive_d_complex_roots=0)
    manifest = dict(initial_HEAD=data['initial_HEAD'], git_status_at_invocation=data['git_status_at_invocation'],
        configuration=config, configuration_sha256=data['configuration_sha256'], criteria=CRITERIA,
        environment=dict(executable=sys.executable, python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        source_heads=data['source_heads'], source_hashes=data['protected_sources'],
        protected_sources_unchanged=True, spectrum_semantics='sorted_positions',
        counters=counters, runtime_seconds=data['total_seconds'],
        K12_control=data.get('K12_control'), output_hashes=data['output_hashes'],
        working_tree_source_sha256={Path(__file__).relative_to(ROOT).as_posix(): sha(__file__)},
        scope='EB elastic screening only; no cross-beta identity, RLB, or new positive-d roots', tests='NOT_RUN_YET')
    write_json(output/'run_manifest.json', manifest)
    return counters


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute', action='store_true', required=True, help='missing-only; completed groups are immutable')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    print(json.dumps(compute(args.output), indent=2))


if __name__ == '__main__':
    main()
