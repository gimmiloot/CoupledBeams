"""D15: six prescribed EB KV roots from three K15 ACTIVE seeds.

Diagnostic-only selection/continuation contract. All complex physics, Newton,
mass normalization and energy gates are the unchanged K12 implementation.
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
from pathlib import Path
import os
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.analysis.laminated_beams import screen_inplane_kelvin_voigt_elastic as screen
from scripts.analysis.laminated_beams.pilot_inplane_kelvin_voigt import clean

SOURCE = screen.OUTPUT
OUTPUT = ROOT/'results/laminated_beams/inplane_kelvin_voigt_targeted_weak_damping'
SELECTION = (('A_STRONG', 0., 'sorted_05'), ('B_INTERMEDIATE', 45., 'sorted_02'),
             ('C_WEAK_ACTIVE', 75., 'sorted_05'))
D_VALUES = (.001, .005)
BUDGET = 1000
ATTEMPTS = 2
WEAK_COMPARISON_RTOL = .01  # descriptive only, and only for d=.001
ARRAYS = ('states', 'reactions', 'a', 'vector')


def protected_hashes():
    paths = [p for folder in (SOURCE, screen.K12) for p in sorted(folder.iterdir())
             if p.suffix in ('.csv', '.json', '.npz')]
    paths += [ROOT/'scripts/lib/inplane_kelvin_voigt.py',
              ROOT/'docs/laminated_beams/inplane_kelvin_voigt_joint_theory.tex']
    return {p.relative_to(ROOT).as_posix(): screen.sha(p) for p in paths}


def select_seeds(rows):
    selected = []
    for state_id, beta, mode in SELECTION:
        matching = [r for r in rows if float(r['beta_deg']) == beta and r['sorted_mode'] == mode]
        if len(matching) != 1:
            raise ValueError('missing/duplicate predetermined seed')
        row = matching[0]
        if (row['model'] != 'EB' or row['activity_status'] != 'ACTIVE' or
            row['root_status'] != 'CONFIRMED' or float(row['d_theta']) != 0 or
            float(row['kappa_theta']) != 1):
            raise ValueError('predetermined seed must be confirmed ACTIVE EB at d=0')
        selected.append(dict(state_id=state_id, source_row=row, beta_deg=beta,
            elastic_sorted_mode=mode, Omega0=float(row['Omega']), omega0=float(row['omega']),
            Lambda0=float(row['Lambda']), a_slope_pred=float(row['a_slope_pred']),
            zeta_slope_pred=float(row['zeta_slope_pred']), shape_key=row['shape_key']))
    active = [float(r['zeta_slope_pred']) for r in rows if r['activity_status'] == 'ACTIVE']
    assert selected[0]['zeta_slope_pred'] == max(active)
    assert selected[2]['zeta_slope_pred'] == min(x for x in active if x > 0)
    return selected


def provider(properties, seed, d, calls):
    if d not in (0., *D_VALUES) or (seed['state_id'], seed['beta_deg'], seed['elastic_sorted_mode']) not in SELECTION:
        raise ValueError('outside predetermined states/d values')
    arm = kv.Arm.reduced('EB', properties)
    return kv.Provider((arm, arm), math.radians(seed['beta_deg']), 1., d, calls)


def predictor(seed, d, previous_z=None):
    if d not in D_VALUES:
        raise ValueError('no additional d values')
    if d == D_VALUES[0]:
        return complex(-seed['a_slope_pred']*d, seed['Omega0'])
    if previous_z is None:
        raise ValueError('second target requires the first converged root')
    return complex(previous_z)


def observables(seed, d, z):
    if d not in D_VALUES or z.imag <= 0:
        raise ValueError('positive-imaginary principal root at a prescribed d required')
    p = z/kv.T_REF
    a, Omega = -z.real, z.imag
    alpha, omega = -p.real, p.imag
    zeta = alpha/math.hypot(alpha, omega)
    nondim = a/math.hypot(a, Omega)
    np.testing.assert_allclose(zeta, nondim, atol=0., rtol=1e-14)
    shift = (Omega-seed['Omega0'])/seed['Omega0']
    dimensional_shift = (omega-seed['omega0'])/seed['omega0']
    np.testing.assert_allclose(shift, dimensional_shift, atol=1e-14, rtol=1e-10)
    a_error = (a/d-seed['a_slope_pred'])/seed['a_slope_pred']
    zeta_error = (zeta/d-seed['zeta_slope_pred'])/seed['zeta_slope_pred']
    status = ('WEAK_DAMPING_CLOSE' if max(abs(a_error), abs(zeta_error)) <= WEAK_COMPARISON_RTOL
              and a > 0 and zeta > 0 else 'WEAK_DAMPING_DEVIATION') if d == .001 else 'DESCRIPTIVE_ONLY'
    return dict(p_re=p.real, p_im=p.imag, z_re=z.real, z_im=z.imag,
        alpha=alpha, omega_d=omega, a=a, Omega_d=Omega, Lambda_d=math.sqrt(Omega),
        zeta=zeta, zeta_dimensionless=nondim,
        a_over_d=a/d, zeta_over_d=zeta/d,
        relative_a_slope_error=a_error, relative_zeta_slope_error=zeta_error,
        frequency_shift=Omega-seed['Omega0'], relative_frequency_shift=shift,
        relative_dimensional_frequency_shift=dimensional_shift,
        frequency_shift_over_d2=shift/d**2, weak_damping_status=status)


def json_write(path, value):
    screen.write_json(path, clean(value))


def save(output, data, shapes, calls, started, previous_seconds):
    data['calls'] = calls.snapshot()
    data['runtime_seconds'] = previous_seconds + time.perf_counter()-started
    with (output/'complex_shapes.npz.tmp').open('wb') as stream:
        np.savez_compressed(stream, **shapes)
    os.replace(output/'complex_shapes.npz.tmp', output/'complex_shapes.npz')
    keys = [f'{sid}_d{i}' for sid, _, _ in SELECTION for i in (1, 2)]
    rows = [data['points'][key]['row'] for key in keys if 'row' in data['points'].get(key, {})]
    stream = io.StringIO(newline='')
    if rows:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)
    screen.atomic(output/'targeted_weak_damping.csv', stream.getvalue())
    json_write(output/'diagnostics.json', data)


def compute(output=OUTPUT, *, retry_unresolved=False):
    output = Path(output).resolve()
    if output != OUTPUT.resolve():
        raise ValueError('this bounded workflow writes only its own output directory')
    checkpoint = output/'diagnostics.json'
    hashes = protected_hashes()
    if checkpoint.exists():
        data = screen.read_json(checkpoint)
        assert data['protected_sources'] == hashes
        assert data['K12_criteria'] == kv.CRITERIA and data['D_VALUES'] == list(D_VALUES)
        eligible = [key for key, point in data['points'].items()
                    if not point['accepted'] and len(point['attempts']) < ATTEMPTS]
        if data.get('finished') and not (retry_unresolved and eligible):
            assert all(screen.sha(output/name) == h for name, h in data['output_hashes'].items())
            return dict(new_principal_roots=0, B=0, B_z=0, form_recoveries=0,
                        auxiliary_positive_d_roots=0, missing_only=True)
    else:
        data = dict(initial_HEAD=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            git_status_at_invocation=subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True),
            protected_sources=hashes, K12_criteria=kv.CRITERIA, D_VALUES=list(D_VALUES),
            descriptive_weak_rtol=WEAK_COMPARISON_RTOL, points={}, seeds={}, calls={},
            principal_positive_d_roots_requested=6, auxiliary_positive_d_roots=0,
            inactive_roots=0, RLB_roots=0, additional_beta=0, additional_d_values=0)
    started = time.perf_counter()
    initial_computed = {key for key, point in data['points'].items() if 'row' in point}
    initial_accepted = sum(p['accepted'] for p in data['points'].values())
    old_seconds = data.get('runtime_seconds', 0.)
    calls = kv.Calls(**data['calls'])
    calls.limit = BUDGET
    properties, config = screen.configuration()  # same native ply reduction
    old = screen.read_json(SOURCE/'run_manifest.json')
    assert config == old['configuration']
    data['configuration'] = dict(config, beta_deg=[s[1] for s in SELECTION], d_theta=list(D_VALUES))
    data['c_values'] = [d*kv.M_REF*kv.T_REF for d in D_VALUES]
    selected = select_seeds(screen.read_csv(SOURCE/'elastic_screening.csv'))
    shapes = dict(np.load(output/'complex_shapes.npz')) if (output/'complex_shapes.npz').exists() else {}
    retry_records = {}
    if retry_unresolved:
        for key in list(data['points']):
            point = data['points'][key]
            if not point['accepted'] and len(point['attempts']) < ATTEMPTS:
                retry_records[key] = data['points'].pop(key)
                for field in ARRAYS:
                    if key+'__'+field in shapes:
                        shapes[key+'__attempt1__'+field] = shapes[key+'__'+field].copy()
    data['finished'] = False
    output.mkdir(parents=True, exist_ok=True)
    with np.load(SOURCE/'screening_shapes.npz', allow_pickle=False) as archive:
        for seed in selected:
            sid = seed['state_id']
            elastic_states = archive[seed['shape_key']+'__states']
            reactions = archive[seed['shape_key']+'__reactions'].ravel()
            zero = provider(properties, seed, 0., calls)
            elastic_a = reactions/zero.reaction_scales
            _, weights = kv.quadrature(elastic_states.shape[1])
            elastic_vector = kv.mass_vector(elastic_states, zero.arms, weights)
            if sid not in data['seeds']:
                B, _ = zero.matrices(1j*seed['Omega0'])
                residual = float(np.linalg.norm(B@elastic_a)/(np.linalg.norm(B)*np.linalg.norm(elastic_a)))
                if residual > kv.CRITERIA['null_residual']:
                    raise ValueError('ELASTIC_SEED_MAPPING_GATE')
                data['seeds'][sid] = dict(seed, mapped_B_residual=residual,
                    mass=float(np.vdot(elastic_vector, elastic_vector).real),
                    source_HEAD=old['initial_HEAD'])
            previous_z, previous_a = None, elastic_a
            previous_vector = elastic_vector
            for index, d in enumerate(D_VALUES, 1):
                key = f'{sid}_d{index}'
                if key in data['points']:
                    point = data['points'][key]
                    if not point.get('accepted'):
                        break
                    row = point['row']
                    previous_z = complex(row['z_re'], row['z_im'])
                    previous_a, previous_vector = shapes[key+'__a'], shapes[key+'__vector']
                    continue
                target = provider(properties, seed, d, calls)
                before = calls.snapshot()
                point_start = time.perf_counter()
                prior = retry_records.get(key)
                point = dict(state_id=sid, d_theta=d, attempts=prior['attempts'][:] if prior else [], accepted=False)
                if prior:
                    point['prior_record'] = prior
                data['points'][key] = point
                try:
                    zguess = predictor(seed, d, previous_z)
                    guess_a = previous_a
                    for attempt in range(len(point['attempts']), ATTEMPTS):
                        if attempt:
                            zguess = (complex(prior['row']['z_re'], prior['row']['z_im']) if prior and 'row' in prior
                                      else complex(-seed['a_slope_pred']*d, seed['Omega0']))
                            guess_a = kv.right_null(target.matrices(zguess)[0])
                        try:
                            correction = kv.correct(target.matrices, zguess, guess_a)
                        except np.linalg.LinAlgError as error:
                            point['attempts'].append(dict(attempt=attempt+1, error=str(error)))
                            continue
                        calls.corrections += correction['steps']
                        point['attempts'].append(dict(attempt=attempt+1, predictor=zguess, correction=correction))
                        if correction['status'] == 'CONVERGED':
                            break
                    if not point['attempts'] or 'correction' not in point['attempts'][-1]:
                        raise RuntimeError('NUMERICAL_UNRESOLVED')
                    z = correction['z']
                    shape = kv.recover(target, z, correction['a'])
                    overlap = np.vdot(elastic_vector, shape['vector'])
                    if abs(overlap):
                        phase = overlap.conjugate()/abs(overlap)
                        for field in ARRAYS:
                            shape[field] *= phase
                    MAC = float(kv.mac_matrix([elastic_vector], [shape['vector']])[0, 0])
                    sequential_MAC = float(kv.mac_matrix([previous_vector], [shape['vector']])[0, 0])
                    diag = kv.diagnose(target, z, shape)
                    full, _ = target.matrices(z)
                    conjugate, _ = target.matrices(z.conjugate())
                    absolute = float(np.linalg.norm(conjugate-full.conj()))
                    criterion = kv.CRITERIA['matrix_atol'] + kv.CRITERIA['transfer_rtol']*np.linalg.norm(full)
                    point['conjugate_matrix_check'] = dict(absolute=absolute,
                        relative=absolute/float(np.linalg.norm(full)), accepted=bool(absolute <= criterion))
                    failures = kv.failures(diag, z, 'ACTIVE', seed['Omega0'], MAC)
                    if correction['status'] != 'CONVERGED':
                        failures.append(correction['status'])
                    if -z.real <= kv.CRITERIA['a_atol']:
                        failures.append('ACTIVE_DECAY_UNRESOLVED')
                    if absolute > criterion:
                        failures.append('CONJUGATE_MATRIX_GATE')
                    root_failures = [f for f in failures if f != 'TRACKING_AMBIGUOUS']
                    status = 'ROOT_ACCEPTED' if not root_failures else 'NUMERICAL_UNRESOLVED'
                    delta = diag['Delta_psi']
                    last = correction['last_delta_z'] or 0j
                    row = dict(state_id=sid, model='EB', beta_deg=seed['beta_deg'],
                        elastic_sorted_mode=seed['elastic_sorted_mode'], d_theta=d, c_theta=target.c,
                        elastic_omega=seed['omega0'], elastic_Omega=seed['Omega0'], elastic_Lambda=seed['Lambda0'],
                        a_slope_pred=seed['a_slope_pred'], zeta_slope_pred=seed['zeta_slope_pred'],
                        **observables(seed, d, z), MAC=MAC, sequential_MAC=sequential_MAC,
                        Delta_psi_re=delta.real, Delta_psi_im=delta.imag, abs_Delta_psi=abs(delta),
                        M_phi=diag['M_phi'], K_phi=diag['K_phi'], C_phi=diag['C_phi'],
                        alpha_energy=diag['alpha_energy'], alpha_energy_discrepancy=(-z.real/kv.T_REF)-diag['alpha_energy'],
                        root_residual=diag['null_residual'], sigma_ratio=diag['sigma_ratio'],
                        physical_residual=max(diag['physical_residuals']), energy_residual=diag['r_E'],
                        conjugate_residual=diag['conjugate_residual'],
                        Newton_iterations=sum(a.get('correction', {}).get('steps', 0) for a in point['attempts']),
                        last_correction_re=last.real, last_correction_im=last.imag,
                        B_calls=calls.B-before['B']+(prior.get('row', {}).get('B_calls', 0) if prior else 0),
                        B_z_calls=calls.B_z-before['B_z']+(prior.get('row', {}).get('B_z_calls', 0) if prior else 0),
                        solver_status=status, mode_status='MODE_CONFIRMED' if MAC >= kv.CRITERIA['MAC'] else 'MODE_IDENTITY_WARNING',
                        K12_gate_status='CONFIRMED' if not failures else 'QUALIFIED', failures=';'.join(failures),
                        source_provenance=seed['source_row']['source_provenance'], shape_key=key)
                    point.update(row=row, diagnostics=diag, accepted=not failures)
                    for field in ARRAYS:
                        shapes[key+'__'+field] = shape[field]
                    previous_z, previous_a, previous_vector = z, shape['a'], shape['vector']
                except (RuntimeError, ValueError) as error:
                    point['error'] = str(error)
                point['calls'] = {k: getattr(calls, k)-before[k] for k in before if k != 'limit'}
                point['seconds'] = time.perf_counter()-point_start
                save(output, data, shapes, calls, started, old_seconds)
                print(key, 'ACCEPTED' if point['accepted'] else 'NUMERICAL_UNRESOLVED',
                      point.get('row', {}).get('weak_damping_status', ''), flush=True)
                if not point['accepted']:
                    break
    data['finished'] = True  # unresolved branches remain stopped, never silently retried
    data['principal_positive_d_roots_computed'] = sum('row' in p for p in data['points'].values())
    data['principal_positive_d_roots_accepted'] = sum(p['accepted'] for p in data['points'].values())
    data['principal_positive_d_roots_unresolved'] = sum(not p['accepted'] for p in data['points'].values())
    data['principal_targets_not_run'] = [f'{sid}_d{i}' for sid, _, _ in SELECTION for i in (1, 2)
                                         if f'{sid}_d{i}' not in data['points']]
    data['stage_status'] = ('COMPLETED_SIX_ROOTS' if data['principal_positive_d_roots_accepted'] == 6
                            else 'PARTIAL_NUMERICAL_QUALIFICATIONS')
    data['protected_sources_unchanged'] = protected_hashes() == hashes
    assert data['protected_sources_unchanged']
    save(output, data, shapes, calls, started, old_seconds)
    data['output_hashes'] = {name: screen.sha(output/name) for name in ('targeted_weak_damping.csv', 'complex_shapes.npz')}
    json_write(checkpoint, data)
    manifest = {key: data[key] for key in ('initial_HEAD', 'git_status_at_invocation', 'protected_sources',
        'protected_sources_unchanged', 'configuration', 'K12_criteria', 'D_VALUES', 'c_values', 'seeds',
        'calls', 'runtime_seconds', 'principal_positive_d_roots_requested', 'principal_positive_d_roots_computed',
        'principal_positive_d_roots_accepted', 'principal_positive_d_roots_unresolved', 'principal_targets_not_run', 'stage_status',
        'auxiliary_positive_d_roots', 'inactive_roots', 'RLB_roots', 'additional_beta', 'additional_d_values', 'output_hashes')}
    manifest.update(environment=dict(executable=sys.executable, python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
        task_scope='three independent ACTIVE EB seeds; continuation only in prescribed d, no cross-beta identity',
        descriptive_weak_rtol=WEAK_COMPARISON_RTOL, max_attempts=ATTEMPTS, matrix_budget=BUDGET,
        working_tree_source_sha256={Path(__file__).relative_to(ROOT).as_posix(): screen.sha(__file__)},
        tests='NOT_RUN_YET', seed_reuse=3, new_elastic_roots=0)
    json_write(output/'run_manifest.json', manifest)
    return dict(new_principal_roots=len({k for k, p in data['points'].items() if 'row' in p}-initial_computed),
                newly_accepted=data['principal_positive_d_roots_accepted']-initial_accepted,
                **calls.snapshot(), auxiliary_positive_d_roots=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute', action='store_true', required=True, help='six prescribed targets, missing-only')
    parser.add_argument('--retry-unresolved', action='store_true', help='one remaining attempt at the same target; prior failures retained')
    args = parser.parse_args()
    print(json.dumps(compute(retry_unresolved=args.retry_unresolved), indent=2))


if __name__ == '__main__':
    main()
