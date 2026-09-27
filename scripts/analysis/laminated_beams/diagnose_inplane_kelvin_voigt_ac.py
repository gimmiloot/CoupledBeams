"""D16: diagnose only saved A/C at d=.001; never continue the K16 pilot."""
from __future__ import annotations
import argparse
import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.analysis.laminated_beams import screen_inplane_kelvin_voigt_elastic as screen
from scripts.analysis.laminated_beams.pilot_inplane_kelvin_voigt import clean

K15 = screen.OUTPUT
K16 = ROOT/'results/laminated_beams/inplane_kelvin_voigt_targeted_weak_damping'
OUTPUT = ROOT/'results/laminated_beams/inplane_kelvin_voigt_ac_diagnostics'


def protected_hashes():
    paths = [p for folder in (K15, K16) for p in sorted(folder.iterdir()) if p.is_file()]
    paths += [ROOT/'scripts/lib/inplane_kelvin_voigt.py',
              ROOT/'scripts/lib/inplane_rotational_spring_eb_modes.py',
              ROOT/'docs/laminated_beams/inplane_kelvin_voigt_joint_theory.tex']
    return {p.relative_to(ROOT).as_posix(): screen.sha(p) for p in paths}


def write_json(path, value):
    screen.write_json(path, clean(value))


def select_inputs(data):
    result = []
    for sid, beta, mode in sd.SELECTION:
        point = data['points'][sid+'_d1']
        row = point['row']
        if (row['beta_deg'], row['elastic_sorted_mode'], row['d_theta'], row['model']) != (beta, mode, .001, 'EB'):
            raise ValueError('A/C selection differs from the prescribed K16 targets')
        if point['accepted'] or 'PHYSICAL_GATE' not in row['failures']:
            raise ValueError('expected historical rejected K16 candidate')
        result.append((sid, point, data['seeds'][sid]))
    return result


def provider(properties, sid, beta, d, calls):
    if d != .001 or (sid, beta, 'sorted_05') not in sd.SELECTION:
        raise ValueError('only A/C at d=.001; no B, .005 or new target')
    arm = kv.Arm.reduced('EB', properties)
    return sd.FullProvider((arm, arm), np.deg2rad(beta), 1., d, calls)


def phase_align(shape, reference):
    overlap = np.vdot(reference, shape['vector'])
    phase = overlap.conjugate()/abs(overlap) if abs(overlap) else 1.
    for name in ('states', 'a', 'reactions', 'vector'):
        shape[name] *= phase


def audit_point(full, original, seed, archive, elastic, calls):
    row = original['row']; key = row['shape_key']
    zold = complex(row['z_re'], row['z_im'])
    ahat = archive[key+'__a']
    plus, minus = sd.project(ahat)
    leakage = float(np.linalg.norm(minus)/np.linalg.norm(plus))
    half = sd.HalfProvider(full, 1)
    B, Bz = half.matrices(zold, derivative=True)
    old_block = sd.spectrum(B, Bz)
    calls.reserve(2)
    projected_shape = kv.recover(full, zold, plus)
    projected = kv.diagnose(full, zold, projected_shape)
    old_half = sd.recover_half(half, zold, kv.right_null(B))
    old_half_diagnostics = kv.diagnose(full, zold, old_half)
    correction = kv.correct(half.matrices, zold, kv.right_null(B))
    calls.corrections += correction['steps']
    z = correction['z']
    shape = sd.recover_half(half, z, correction['a'])
    _, weights = kv.quadrature()
    elastic_vector = kv.mass_vector(elastic[seed['shape_key']+'__states'], full.arms, weights)
    phase_align(shape, elastic_vector)
    diag = kv.diagnose(full, z, shape)
    mac = float(kv.mac_matrix([elastic_vector], [shape['vector']])[0, 0])
    mac_old = float(kv.mac_matrix([archive[key+'__vector']], [shape['vector']])[0, 0])
    Bplus, Bplus_z = half.matrices(z, derivative=True)
    Bminus, Bminus_z = sd.HalfProvider(full, -1).matrices(z, derivative=True)
    fullB, fullBz = full.matrices(z, derivative=True)
    blocks = dict(full=sd.spectrum(fullB, fullBz, shape['a']),
        plus=sd.spectrum(Bplus, Bplus_z, shape['a'][:3]),
        minus=sd.spectrum(Bminus, Bminus_z))
    conjugate, _ = half.matrices(z.conjugate())
    same_algorithm, _ = half.matrices(z)
    conjugate_error = float(np.linalg.norm(conjugate-same_algorithm.conj())/np.linalg.norm(same_algorithm))
    physical = max(diag['physical_residuals'])
    # Triggers declared before the optional stage; no automatic extra root solve.
    needs_transfer = (abs(z-zold)/abs(zold) > kv.CRITERIA['newton_step'] or
        physical > kv.CRITERIA['physical_residual'] or
        (leakage < kv.CRITERIA['symmetry_defect'] and row['physical_residual'] > kv.CRITERIA['physical_residual']))
    result = dict(original_row=row, original_diagnostics=original['diagnostics'],
        original_attempts=original['attempts'], seed=seed, same_p_plus=old_block,
        same_p_half_recovery=old_half_diagnostics, projection=dict(norm_plus=float(np.linalg.norm(plus)),
        norm_minus=float(np.linalg.norm(minus)), symmetry_leakage=leakage,
        reconstruction_error=float(np.linalg.norm(plus+minus-ahat)/np.linalg.norm(ahat))),
        projected_diagnostics=projected, correction=correction, half_diagnostics=diag,
        half_physical=sd.physical_details(half, z, shape), spectra=blocks,
        MAC_elastic_half=mac, MAC_full_half=mac_old, half_conjugacy_relative=conjugate_error,
        full_gate_failures_on_half_shape=kv.failures(diag, z, 'ACTIVE', seed['Omega0'], mac),
        analytic_transfer_trigger=bool(needs_transfer), analytic_transfer_ran=False,
        attempts=1)
    result['row'] = dict(state=row['state_id'], beta_deg=row['beta_deg'], d_theta=.001,
        p_full_re=zold.real/kv.T_REF, p_full_im=zold.imag/kv.T_REF,
        p_half_re=z.real/kv.T_REF, p_half_im=z.imag/kv.T_REF,
        z_full_re=zold.real, z_full_im=zold.imag, z_half_re=z.real, z_half_im=z.imag,
        root_difference=abs(z-zold)/kv.T_REF, relative_root_difference=abs(z-zold)/abs(zold),
        full_rB=row['root_residual'], half_rB=blocks['plus']['null_residual'],
        full_sigma1_ratio=diag['sigma_ratio'], full_sigma2_ratio=diag['next_sigma_ratio'],
        plus_sigma1_ratio=blocks['plus']['ratios'][-1], plus_sigma2_ratio=blocks['plus']['ratios'][-2],
        minus_sigma1_ratio=blocks['minus']['ratios'][-1],
        full_physical_residual_original=row['physical_residual'],
        full_physical_residual_projected=max(projected['physical_residuals']),
        full_physical_residual_lifted_half=physical,
        symmetry_defect_original=original['diagnostics']['symmetry_defect'],
        symmetry_leakage=leakage, MAC_full_half=mac_old, MAC_elastic_half=mac,
        energy_residual_half=diag['r_E'], alpha_energy_error_half=-z.real/kv.T_REF-diag['alpha_energy'],
        diagnostic_status='REDUCED_PROBLEM_NUMERICALLY_UNRESOLVED' if physical > kv.CRITERIA['physical_residual'] else 'HALF_AND_FULL_CONFIRMED')
    return result, shape


def save(data, shapes, calls, started, previous_seconds):
    data['calls'] = calls.snapshot()
    data['runtime_seconds'] = previous_seconds+time.perf_counter()-started
    with (OUTPUT/'ac_reduced_shapes.npz.tmp').open('wb') as stream:
        np.savez_compressed(stream, **shapes)
    os.replace(OUTPUT/'ac_reduced_shapes.npz.tmp', OUTPUT/'ac_reduced_shapes.npz')
    for filename, rows in [('ac_symmetry_diagnostics.csv', [p['row'] for p in data['points'].values()]),
        ('singular_values.csv', [dict(state=sid, block=name, sorted_singular_index=i+1,
            singular_value=value, ratio=s['ratios'][i], balanced_ratio=s['balanced_ratios'][i])
            for sid, p in data['points'].items() for name, s in p['spectra'].items()
            for i, value in enumerate(s['singular_values'])])]:
        stream = io.StringIO(newline='')
        if rows:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
            writer.writeheader(); writer.writerows(rows)
        screen.atomic(OUTPUT/filename, stream.getvalue())
    write_json(OUTPUT/'diagnostics.json', data)


def transfer_control(full, result, shape):
    """Conditional comparison at the unchanged candidate; no new root here."""
    z = complex(result['row']['z_half_re'], result['row']['z_half_im'])
    arm = full.arms[0]; p = z/kv.T_REF
    full.calls.reserve(6)
    direct = full.transfer(z, arm)
    frechet, frechet_z = full.transfer(z, arm, derivative=True)
    full.calls.analytic_transfer += 1
    closed, closed_p = sd.closed_transfer(p, arm, derivative=True)
    scale = arm.scale()
    norm = np.linalg.norm(closed*scale[None, :]/scale[:, None])
    comparisons = {}
    half = sd.HalfProvider(full, 1)
    initial = np.r_[np.zeros(3), shape['reactions'][:3]]
    for name, T in [('direct_expm', direct), ('frechet_expm', frechet), ('closed_form', closed)]:
        end = T@initial
        lifted = dict(states=np.array([[end], [end@sd.F]]))
        error = (T-closed)*scale[None, :]/scale[:, None]
        comparisons[name] = dict(transfer_relative=float(np.linalg.norm(error)/norm),
            transfer_absolute=float(np.linalg.norm(error)),
            physical=sd.physical_details(half, z, lifted))
    dz = (frechet_z-closed_p/kv.T_REF)*scale[None, :]/scale[:, None]
    comparisons['frechet_derivative_relative'] = float(np.linalg.norm(dz)/np.linalg.norm(closed_p/kv.T_REF*scale[None, :]/scale[:, None]))
    comparisons['step_expm_physical'] = result['half_physical']
    return comparisons


def closed_attempt(full, result, original_shape, elastic_vector):
    """Last allowed attempt at the same root, after an observed transfer trigger."""
    zold = complex(result['row']['z_full_re'], result['row']['z_full_im'])
    half = sd.ClosedHalfProvider(full, 1, triggered=result['analytic_transfer_trigger'])
    frozen = sd.FrozenBalanced(half, zold)
    initial = kv.right_null(frozen.matrices(zold)[0])
    correction = kv.correct(frozen.matrices, zold, initial)
    full.calls.corrections += correction['steps']
    z, a = correction['z'], frozen.reactions(correction['a'])
    shape = sd.recover_closed(half, z, a)
    phase_align(shape, elastic_vector)
    diag = kv.diagnose(full, z, shape)
    B, Bz = half.matrices(z, derivative=True)
    minus = sd.ClosedHalfProvider(full, -1, triggered=True)
    Bm, Bmz = minus.matrices(z, derivative=True)
    conjugate = half.matrices(z.conjugate())[0]
    blocks = dict(plus_closed=sd.spectrum(B, Bz, shape['a'][:3]),
                  minus_closed=sd.spectrum(Bm, Bmz))
    # Full-B gate is retained. Rank assessment is recorded separately by class.
    mac = float(kv.mac_matrix([elastic_vector], [shape['vector']])[0, 0])
    old_mac = float(kv.mac_matrix([original_shape['vector']], [shape['vector']])[0, 0])
    gates = kv.failures(diag, z, 'ACTIVE', result['seed']['Omega0'], mac)
    closed = dict(correction=correction, frozen_row_scales=frozen.rows.tolist(),
        frozen_column_scales=frozen.cols.tolist(), half_diagnostics=diag,
        physical=sd.physical_details(half, z, shape), spectra=blocks,
        half_conjugate_residual=float(np.linalg.norm(conjugate@shape['a'][:3].conj())/(np.linalg.norm(conjugate)*np.linalg.norm(shape['a'][:3]))),
        conjugacy_relative=float(np.linalg.norm(conjugate-B.conj())/np.linalg.norm(B)),
        full_K12_gate_failures=gates, MAC_elastic_half=mac, MAC_full_half=old_mac)
    result.update(closed_attempt=closed, attempts=2)
    result['spectra'].update(blocks)
    # Preserve the first reduction result in full; do not overwrite its failure.
    result['frechet_half_row'] = result['row'].copy()
    row = result['row']
    row.update(p_half_re=z.real/kv.T_REF, p_half_im=z.imag/kv.T_REF,
        z_half_re=z.real, z_half_im=z.imag, root_difference=abs(z-zold)/kv.T_REF,
        relative_root_difference=abs(z-zold)/abs(zold), half_rB=blocks['plus_closed']['null_residual'],
        plus_sigma1_ratio=blocks['plus_closed']['ratios'][-1],
        plus_sigma2_ratio=blocks['plus_closed']['ratios'][-2],
        minus_sigma1_ratio=blocks['minus_closed']['ratios'][-1],
        full_physical_residual_lifted_half=max(diag['physical_residuals']),
        MAC_full_half=old_mac, MAC_elastic_half=mac, energy_residual_half=diag['r_E'],
        alpha_energy_error_half=-z.real/kv.T_REF-diag['alpha_energy'])
    return shape


def finalize_point(full, result, shape, elastic_rows):
    z = complex(result['row']['z_half_re'], result['row']['z_half_im'])
    B, _ = full.matrices(z)
    _, singular, Vh = np.linalg.svd(B)
    v = Vh.conj().T[:, -2]
    plus, minus = sd.project(v)
    rank = dict(second_full_singular_value=float(singular[-2]),
        second_full_right_plus_norm=float(np.linalg.norm(plus)),
        second_full_right_minus_norm=float(np.linalg.norm(minus)))
    spectra = result['closed_attempt']['spectra']
    ps, ms = spectra['plus_closed'], spectra['minus_closed']
    rank['minus_linearized_distance_z'] = ms['singular_values'][-1]/ms['left_Bz_right']
    rank['plus_linearized_correction_scale_z'] = ps['singular_values'][-1]/ps['left_Bz_right']
    # This is a rank/conditioning interpretation, not a new physical gate.
    rank['same_block_suspect'] = (ps['balanced_ratios'][-2] < kv.CRITERIA['simple_sigma_separation'] or
        rank['plus_linearized_correction_scale_z']/abs(z) > kv.CRITERIA['frequency_rtol'])
    rank['opposite_block_near'] = ms['balanced_ratios'][-1] < kv.CRITERIA['simple_sigma_separation']
    rank['second_sigma_unexplained'] = not (np.linalg.norm(minus) > .99 and
        not rank['same_block_suspect'] and not rank['opposite_block_near'])
    if result['row']['state'] == 'A_STRONG':
        # No A angular investigation is authorized or needed.
        rank['second_sigma_unexplained'] = False
    rank['local_beta_trigger'] = (result['row']['state'] == 'C_WEAK_ACTIVE' and
        sd.beta_trigger(**{k:rank[k] for k in ('same_block_suspect','opposite_block_near','second_sigma_unexplained')}))
    rank['allowed_beta_plan'] = sd.local_beta_plan(rank['local_beta_trigger'])
    rows = [r for r in elastic_rows if float(r['beta_deg']) == result['row']['beta_deg'] and int(r['symmetry_eta']) == -1]
    nearest = min(rows, key=lambda r:abs(float(r['Omega'])-z.imag))
    rank['read_only_K15_nearest_minus'] = dict(row=nearest,
        signed_gap_to_complex_Omega=float(nearest['Omega'])-z.imag,
        note='existing elastic eta=-1 root is independent of d for identical arms; not a new solve')
    # Is the recurrence itself the cause? Propagate the corrected reactions with
    # the unchanged full recovery, without any symmetry projection in that call.
    full.calls.reserve(2)
    independent = kv.recover(full, z, shape['a'])
    independent_diag = kv.diagnose(full, z, independent)
    result['independent_full_recovery_with_closed_reactions'] = independent_diag
    result['rank_interpretation'] = rank
    gates = result['closed_attempt']['full_K12_gate_failures']
    physical_failures = [g for g in gates if g != 'POSSIBLE_MULTIPLICITY']
    if result['closed_attempt']['correction']['status'] != 'CONVERGED':
        physical_failures.append('NEWTON_LIMIT')
    if ps['null_residual'] > kv.CRITERIA['null_residual'] or ps['ratios'][-1] > kv.CRITERIA['sigma_ratio']:
        physical_failures.append('HALF_ROOT_GATE')
    if result['closed_attempt']['half_conjugate_residual'] > kv.CRITERIA['null_residual']:
        physical_failures.append('HALF_CONJUGATE_GATE')
    result['diagnostic_physical_failures'] = physical_failures
    if physical_failures:
        status = 'REDUCED_PROBLEM_NUMERICALLY_UNRESOLVED'
    elif rank['same_block_suspect']:
        status = 'POSSIBLE_SAME_CLASS_MULTIPLICITY'
    elif rank['opposite_block_near']:
        status = 'POSSIBLE_OPPOSITE_CLASS_NEARBY_ROOT'
    elif result['row']['symmetry_leakage'] >= kv.CRITERIA['symmetry_defect']:
        status = 'FULL_NULLVECTOR_SYMMETRY_LEAKAGE'
    else:
        status = 'FULL_TRANSFER_RECOVERY_CONDITIONING'
    result['row']['diagnostic_status'] = status
    result['diagnostic_finalized'] = True


def compute():
    checkpoint = OUTPUT/'diagnostics.json'; hashes = protected_hashes()
    if checkpoint.exists():
        data = screen.read_json(checkpoint)
        assert data['source_hashes'] == hashes and data['criteria'] == kv.CRITERIA
        if data.get('finished'):
            assert all(screen.sha(OUTPUT/p) == h for p, h in data['output_hashes'].items())
            return dict(new_targets=0, half_resolves=0, full_B=0, half_B=0, B_z=0, expm=0, recoveries=0, missing_only=True)
        counters = {k: v for k, v in data['calls'].items() if k != 'total_build_equivalents'}
        calls = sd.Calls(**counters)
        shapes = dict(np.load(OUTPUT/'ac_reduced_shapes.npz', allow_pickle=False))
    else:
        data = dict(initial_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            initial_git_status='', git_status_at_invocation=subprocess.check_output(['git','status','--short'],text=True),
            source_hashes=hashes, criteria=kv.CRITERIA, points={}, algebra=[],
            environment=dict(executable=sys.executable, python=sys.version, numpy=np.__version__, scipy=scipy.__version__),
            new_physical_targets=0, new_d_values=0, B_resolves=0, RLB_roots=0,
            local_beta_roots=0, local_beta_ran=False, analytic_transfer_ran=False,
            fixed_half_units=sd.HALF_UNITS.tolist(), reaction_reflection=sd.R.tolist(),
            rank_diagnostic='one row-2norm then column-2norm equilibration, diagnostic only; fixed units in Newton',
            tests='NOT_RUN_YET')
        calls, shapes = sd.Calls(), {}
    started = time.perf_counter(); previous_seconds = data.get('runtime_seconds', 0.)
    old = screen.read_json(K16/'diagnostics.json')
    properties, config = screen.configuration()
    assert config == screen.read_json(K15/'run_manifest.json')['configuration']
    assert dict(old['configuration'], beta_deg=config['beta_deg'], d_theta=config['d_theta']) == config
    data['configuration'] = dict(config, beta_deg=[0., 75.], d_theta=.001)
    inputs = select_inputs(old)
    data['input_rows'] = {sid:point['row'] for sid, point, seed in inputs}
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if not data['algebra']:
        for sid, point, seed in inputs:
            full = provider(properties, sid, seed['beta_deg'], .001, calls)
            for z in (-.12+7j, -.2+80j):
                data['algebra'].append(sd.algebra_check(full, z))
    if len(data['algebra']) != 4 or not all(a['accepted'] for a in data['algebra']):
        save(data, shapes, calls, started, previous_seconds)
        raise RuntimeError('REDUCTION_ALGEBRA_GATE; no root solve allowed')
    with np.load(K16/'complex_shapes.npz', allow_pickle=False) as archive, np.load(K15/'screening_shapes.npz', allow_pickle=False) as elastic:
        for sid, point, seed in inputs:
            if sid in data['points']:
                continue
            full = provider(properties, sid, seed['beta_deg'], .001, calls)
            result, shape = audit_point(full, point, seed, archive, elastic, calls)
            data['points'][sid] = result
            for name in ('states', 'a', 'reactions', 'vector'):
                shapes[sid+'__'+name] = shape[name]
            save(data, shapes, calls, started, previous_seconds)
            print(sid, result['row'], flush=True)
    data['half_audit_ran'] = True
    for sid, point, seed in inputs:
        result = data['points'][sid]
        if result['analytic_transfer_trigger'] and not result['analytic_transfer_ran']:
            full = provider(properties, sid, seed['beta_deg'], .001, calls)
            shape = {name:shapes[sid+'__'+name] for name in ('states', 'a', 'reactions', 'vector')}
            result['transfer_control'] = transfer_control(full, result, shape)
            result['analytic_transfer_ran'] = True
            data['analytic_transfer_ran'] = True
            save(data, shapes, calls, started, previous_seconds)
            print(sid, 'TRANSFER', clean(result['transfer_control']), flush=True)
    with np.load(K16/'complex_shapes.npz', allow_pickle=False) as archive, np.load(K15/'screening_shapes.npz', allow_pickle=False) as elastic:
        for sid, point, seed in inputs:
            result = data['points'][sid]
            if not result['analytic_transfer_trigger'] or 'closed_attempt' in result:
                continue
            full = provider(properties, sid, seed['beta_deg'], .001, calls)
            _, weights = kv.quadrature()
            elastic_vector = kv.mass_vector(elastic[seed['shape_key']+'__states'], full.arms, weights)
            original_shape = {'vector':archive[point['row']['shape_key']+'__vector']}
            shape = closed_attempt(full, result, original_shape, elastic_vector)
            for name in ('states', 'a', 'reactions', 'vector'):
                shapes[sid+'__frechet__'+name] = shapes[sid+'__'+name].copy()
                shapes[sid+'__'+name] = shape[name]
            save(data, shapes, calls, started, previous_seconds)
            print(sid, 'CLOSED_ATTEMPT', result['row'], 'GATES', result['closed_attempt']['full_K12_gate_failures'],
                  'BALANCED', {k:v['balanced_ratios'] for k,v in result['closed_attempt']['spectra'].items()}, flush=True)
    for sid, point, seed in inputs:
        result = data['points'][sid]
        if not result.get('diagnostic_finalized'):
            full = provider(properties, sid, seed['beta_deg'], .001, calls)
            shape = {name:shapes[sid+'__'+name] for name in ('states','a','reactions','vector')}
            finalize_point(full, result, shape, screen.read_csv(K15/'elastic_screening.csv'))
            save(data, shapes, calls, started, previous_seconds)
            print(sid, 'FINAL', result['row']['diagnostic_status'], result['rank_interpretation'], flush=True)
    data['pending_conditional_transfer'] = [sid for sid, p in data['points'].items() if p['analytic_transfer_trigger'] and not p['analytic_transfer_ran']]
    data['local_beta_trigger'] = data['points']['C_WEAK_ACTIVE']['rank_interpretation']['local_beta_trigger']
    if data['local_beta_trigger']:
        # A trigger is a stop for review of the bounded optional plan, never an
        # implicit global scan. This run has not calculated any optional root.
        data['stage_status'] = 'CONDITIONAL_ELASTIC_DIAGNOSTIC_PENDING'
    else:
        data['stage_status'] = 'AC_DIAGNOSTIC_QUALIFICATION_COMPLETE'
    data['finished'] = not data['pending_conditional_transfer'] and not data['local_beta_trigger']
    assert protected_hashes() == hashes
    save(data, shapes, calls, started, previous_seconds)
    data['output_hashes'] = {name:screen.sha(OUTPUT/name) for name in
        ('ac_symmetry_diagnostics.csv','singular_values.csv','ac_reduced_shapes.npz')}
    data['source_hashes_unchanged'] = protected_hashes() == hashes
    write_json(checkpoint, data)
    write_json(OUTPUT/'run_manifest.json', {k:v for k,v in data.items() if k != 'points'})
    return dict(calls=calls.snapshot(), pending_conditional_transfer=data['pending_conditional_transfer'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute', action='store_true', required=True, help='missing-only A/C diagnostic')
    parser.parse_args()
    print(json.dumps(compute(), indent=2))


if __name__ == '__main__':
    main()
