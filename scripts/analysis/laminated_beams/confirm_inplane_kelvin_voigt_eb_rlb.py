"""D22: six fixed reduced complex roots and four read-only K12/K19 states.

This bounded comparison has a different input/output contract from D18:
matched K15/K22 seeds, two theories, one d, five cross-theory comparisons.
It adds no root algorithm or physical equations. --compute is missing-only.
"""
from __future__ import annotations
import argparse
from dataclasses import fields
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import scipy
from scripts.analysis.laminated_beams import complete_inplane_kelvin_voigt_weak_damping as completion
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as production
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd

screen = completion.screen
write_json = completion.old.json_write
BASE = ROOT/'results/laminated_beams'
K12 = BASE/'inplane_kelvin_voigt_pilot'
K15 = BASE/'inplane_kelvin_voigt_elastic_screening'
K19 = BASE/'inplane_kelvin_voigt_targeted_weak_damping_completion'
K22 = BASE/'inplane_kelvin_voigt_rlb_elastic_screening'
OUTPUT = BASE/'inplane_kelvin_voigt_eb_rlb_complex_confirmation'
CASES = (('R0', 5., '01'), ('R1', 0., '05'), ('R2', 45., '05'),
         ('R3', 75., '05'), ('R4', 45., '03'))
NEW_TARGETS = ('R1_RLB', 'R2_EB', 'R2_RLB', 'R3_RLB', 'R4_EB', 'R4_RLB')
D = .001
RATIO_DESCRIPTIVE_RTOL = .01  # declared before calculation; never a root gate
MAX_ATTEMPTS = 2
BUDGET = 5000  # total counted build-equivalents, no new solver tolerance


def initial_state():
    path = OUTPUT/'initial_state.json'
    if path.exists():
        return screen.read_json(path)
    paths = [p for folder in (K12, K15, K19, K22) for p in folder.iterdir() if p.is_file()]
    paths += [ROOT/p for p in ('scripts/lib/inplane_kelvin_voigt.py',
        'scripts/lib/inplane_kelvin_voigt_solver.py',
        'scripts/lib/inplane_kelvin_voigt_symmetry_diagnostics.py',
        'scripts/analysis/laminated_beams/complete_inplane_kelvin_voigt_weak_damping.py')]
    value = dict(initial_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        initial_git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        protected_sources={p.relative_to(ROOT).as_posix():screen.sha(p) for p in paths})
    OUTPUT.mkdir(parents=True, exist_ok=True); write_json(path, value)
    return value


def inputs():
    """Exact saved mode mapping, never nearest-frequency selection."""
    matches = screen.read_csv(K22/'eb_rlb_matched_comparison.csv')
    tables = dict(EB=screen.read_csv(K15/'elastic_screening.csv'),
                  RLB=screen.read_csv(K22/'rlb_elastic_screening.csv'))
    seeds, selected = {}, {}
    for cid, beta, index in CASES:
        match, = [r for r in matches if float(r['beta_deg']) == beta and r['eb_sorted_mode'] == 'sorted_'+index]
        assert match['rlb_sorted_mode'] == 'rlb_sorted_'+index and match['match_status'] == 'CONFIRMED'
        assert int(match['eta']) == 1
        selected[cid] = match
        for theory in ('EB','RLB'):
            field = 'sorted_mode' if theory == 'EB' else 'rlb_sorted_mode'
            mode = match['eb_sorted_mode' if theory == 'EB' else 'rlb_sorted_mode']
            row, = [r for r in tables[theory] if float(r['beta_deg']) == beta and r[field] == mode]
            assert row['root_status'] == 'CONFIRMED' and row['activity_status'] == 'ACTIVE'
            assert float(row['kappa_theta']) == 1 and float(row['d_theta']) == 0
            assert int(row['symmetry_eta' if theory == 'EB' else 'eta']) == 1
            assert float(row['Omega']) == float(match['Omega_'+theory])
            assert float(row['zeta_slope_pred']) == float(match['G_'+theory])
            key = cid+'_'+theory
            seeds[key] = dict(key=key, case_id=cid, theory=theory, beta_deg=beta,
                matched_mode=mode, eta=1, Omega0=float(row['Omega']), G=float(row['zeta_slope_pred']),
                a_slope_pred=float(row['a_slope_pred']), source_row=row)
    for theory, path in (('EB', K15/'screening_shapes.npz'), ('RLB', K22/'rlb_elastic_shapes.npz')):
        with np.load(path, allow_pickle=False) as shapes:
            for key in NEW_TARGETS:
                seed = seeds[key]
                if seed['theory'] == theory:
                    seed['states'] = shapes[seed['source_row']['shape_key']+'__states']
    return seeds, selected


def config_for(properties, seed, d=D):
    if seed['key'] not in NEW_TARGETS or d != D or seed['eta'] != 1:
        raise ValueError('only the six predetermined reduced eta=+1 targets at d=.001')
    cid, beta, index = next(c for c in CASES if c[0] == seed['case_id'])
    expected_mode = ('rlb_' if seed['theory'] == 'RLB' else '')+'sorted_'+index
    if (seed['beta_deg'], seed['matched_mode'], seed['key']) != (beta, expected_mode, cid+'_'+seed['theory']):
        raise ValueError('no new beta/mode/theory or asymmetric construction')
    arm = kv.Arm.reduced(seed['theory'], properties)
    return production.Config((arm,arm), np.deg2rad(beta), 1., d, mu=0.)


def observables(seed, d, z, p=None):
    p = z/kv.T_REF if p is None else p
    assert z.imag > 0 and p.imag > 0 and d > 0
    a, omega = -z.real, z.imag
    zeta = -p.real/math.hypot(p.real, p.imag)
    assert math.isclose(zeta, a/math.hypot(a, omega), rel_tol=1e-14)
    return dict(case_id=seed['case_id'], theory=seed['theory'], beta_deg=seed['beta_deg'],
        matched_mode=seed['matched_mode'], eta=1, d_theta=d, c_theta=d*kv.M_REF*kv.T_REF,
        Omega_0=seed['Omega0'], G=seed['G'], G_prediction=seed['G'],
        p_re=p.real, p_im=p.imag, z_re=z.real, z_im=z.imag, alpha=-p.real, omega_d=p.imag,
        a=a, Omega_d=omega, zeta=zeta, zeta_over_d=zeta/d,
        relative_predictor_error=(zeta/d-seed['G'])/seed['G'],
        relative_frequency_shift=(omega-seed['Omega0'])/seed['Omega0'])


def reused_rows(seeds):
    rows, sources = [], {}
    pilot = screen.read_csv(K12/'modal_results.csv')
    common = set.intersection(*[set(float(r['d_theta']) for r in pilot if
        r['model'] == t and r['role'] == 'ACTIVE' and r['status'] == 'CONFIRMED' and float(r['d_theta']) > 0)
        for t in ('EB','RLB')])
    d0 = min(common)
    for theory in ('EB','RLB'):
        r, = [r for r in pilot if r['model'] == theory and r['role'] == 'ACTIVE' and float(r['d_theta']) == d0]
        assert r['status'] == 'CONFIRMED' and float(r['beta0_deg']) == 5 and float(r['kappa_theta']) == 1
        seed = seeds['R0_'+theory]; assert float(r['Omega0']) == seed['Omega0']
        row = observables(seed, d0, complex(float(r['z_real']),float(r['z_imag'])),
                          complex(float(r['p_real']),float(r['p_imag'])))
        row.update(source='K12', source_row_key=r['key'], solver_path='HISTORICAL_FULL_K12',
            MAC_to_elastic=float(r['MAC']), root_residual=float(r['r_B']), sigma_ratio=float(r['sigma_ratio']),
            physical_residual=float(r['physical_residual']), reduced_physical_residual=None,
            lifted_full_physical_residual=None, energy_residual=float(r['r_E']),
            full_control_qualification=r['failures'], root_status='ROOT_ACCEPTED',
            M_phi=float(r['M_phi']), K_phi=float(r['K_phi']), C_phi=float(r['C_phi']),
            alpha_energy=float(r['alpha_energy']), Delta_psi_re=float(r['Delta_psi_real']),
            Delta_psi_im=float(r['Delta_psi_imag']), notes='Historical K12 minimum common d; no new evaluation')
        rows.append(row); sources[seed['key']] = r
    prior = screen.read_csv(K19/'combined_six_state_summary.csv')
    for cid, state in (('R1','A'),('R3','C')):
        r, = [r for r in prior if r['state'] == state and float(r['d_theta']) == D]
        assert r['root_status'] == 'ROOT_ACCEPTED'
        seed = seeds[cid+'_EB']; assert float(r['Omega_0']) == seed['Omega0']
        row = observables(seed,D,complex(float(r['z_re']),float(r['z_im'])),complex(float(r['p_re']),float(r['p_im'])))
        for name in ('root_residual','sigma_ratio','conjugate_residual','energy_residual',
                     'reduced_physical_residual','lifted_full_physical_residual','full_control_physical_residual',
                     'M_phi','K_phi','C_phi','alpha_energy','Delta_psi_re','Delta_psi_im'):
            row[name] = float(r[name])
        row.update(source='K19', source_row_key=state+'/.001', solver_path=r['solver_path'],
            MAC_to_elastic=float(r['MAC']), physical_residual=float(r['lifted_full_physical_residual']),
            full_control_qualification=r['full_control_qualification'], root_status=r['root_status'], notes=r['notes'])
        rows.append(row); sources[seed['key']] = r
    return rows, sources


def new_row(seed, result, control):
    row = observables(seed,D,complex(result['z']))
    diag = result['diagnostics']; delta = diag['Delta_psi']
    row.update(source='NEW', solver_path=result['solver_path'], MAC_to_elastic=result['MAC'],
        root_status='ROOT_ACCEPTED' if result['accepted'] else 'NUMERICAL_UNRESOLVED',
        root_residual=diag['null_residual'], sigma_ratio=diag['sigma_ratio'],
        next_sigma_ratio=diag['next_sigma_ratio'], conjugate_residual=diag['conjugate_residual'],
        physical_residual=max(diag['physical_residuals']),
        reduced_physical_residual=max(result['reduced_physical']['half_normalized']),
        lifted_full_physical_residual=max(diag['physical_residuals']), energy_residual=diag['r_E'],
        M_phi=diag['M_phi'], K_phi=diag['K_phi'], C_phi=diag['C_phi'], alpha_energy=diag['alpha_energy'],
        Delta_psi_re=delta.real, Delta_psi_im=delta.imag,
        full_control_physical_residual=max(control['diagnostics']['physical_residuals']) if 'diagnostics' in control else None,
        full_control_qualification=';'.join(control['failures']),
        notes='Matched elastic seed; direct d=0 to .001; full control evaluation only')
    return row


def comparisons(rows, matches):
    lookup = {r['case_id']+'_'+r['theory']:r for r in rows}
    result = []
    for cid,beta,index in CASES:
        m = matches[cid]; e = lookup.get(cid+'_EB',{}); r = lookup.get(cid+'_RLB',{})
        row = dict(case_id=cid,beta_deg=beta,eb_mode=m['eb_sorted_mode'],rlb_mode=m['rlb_sorted_mode'],
            d_theta=e.get('d_theta',D),Omega0_EB=float(m['Omega_EB']),Omega0_RLB=float(m['Omega_RLB']),
            elastic_delta_Omega=float(m['delta_Omega']),G_EB=float(m['G_EB']),G_RLB=float(m['G_RLB']),
            G_ratio=float(m['G_ratio']),elastic_delta_G=float(m['delta_G']),
            eb_source=e.get('source'),rlb_source=r.get('source'),
            eb_root_status=e.get('root_status','NOT_RUN'),rlb_root_status=r.get('root_status','NOT_RUN'))
        if e.get('root_status') == r.get('root_status') == 'ROOT_ACCEPTED':
            assert e['d_theta'] == r['d_theta']
            ratio = r['zeta']/e['zeta']; error = (ratio-row['G_ratio'])/row['G_ratio']
            row.update(zeta_EB=e['zeta'],zeta_RLB=r['zeta'],zeta_ratio=ratio,complex_delta_zeta=ratio-1,
                ratio_error=error,Omega_d_EB=e['Omega_d'],Omega_d_RLB=r['Omega_d'],
                damped_delta_Omega=(r['Omega_d']-e['Omega_d'])/e['Omega_d'],
                comparison_status='CLOSE_TO_ELASTIC_PREDICTION' if abs(error)<=RATIO_DESCRIPTIVE_RTOL else 'DESCRIPTIVE_DEVIATION',
                notes='Historical common K12 d' if cid=='R0' else 'd=.001; fixed K22 matching')
        else:
            row.update(comparison_status='NUMERICAL_UNRESOLVED',notes='No ratio substituted for an unaccepted state')
        result.append(row)
    return result


def save(data, shapes, calls):
    data['calls'] = calls.snapshot()
    new = [p['row'] for p in data['points'].values() if 'row' in p]
    data['all_states'] = data['reused_rows']+new
    completion.write_csv(OUTPUT/'new_complex_roots.csv',new)
    completion.write_csv(OUTPUT/'representative_complex_comparison.csv',comparisons(data['all_states'],data['matching_source_rows']))
    with (OUTPUT/'new_complex_shapes.npz.tmp').open('wb') as stream:
        np.savez_compressed(stream,**shapes)
    os.replace(OUTPUT/'new_complex_shapes.npz.tmp',OUTPUT/'new_complex_shapes.npz')
    write_json(OUTPUT/'diagnostics.json',data)


def compute():
    initial = initial_state()
    assert all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    checkpoint = OUTPUT/'diagnostics.json'
    if checkpoint.exists():
        data = screen.read_json(checkpoint)
        assert data['criteria'] == kv.CRITERIA and data['targets'] == list(NEW_TARGETS)
        if data.get('finished'):
            assert all(screen.sha(OUTPUT/p)==h for p,h in data['output_hashes'].items())
            return dict(missing_only=True,new_roots=0,matrix_calls=0,form_recoveries=0)
        shapes = dict(np.load(OUTPUT/'new_complex_shapes.npz',allow_pickle=False))
    else:
        data = dict(**{k:v for k,v in initial.items() if k!='initial_memory'},criteria=kv.CRITERIA,
            targets=list(NEW_TARGETS),points={},max_attempts=MAX_ATTEMPTS,budget=BUDGET,
            ratio_descriptive_rtol=RATIO_DESCRIPTIVE_RTOL,full_control_evaluations=0,
            full_newton_roots=0,auxiliary_roots=0,new_beta=0,new_d_outside_authorized=0,
            asymmetric_positive_d=0,elastic_recomputations=0,solver_changes=0,high_precision=0,
            runtime_seconds=0.,tests='NOT_RUN_YET',environment=dict(executable=sys.executable,
                python=sys.version,numpy=np.__version__,scipy=scipy.__version__))
        shapes = {}
    calls = sd.Calls(**{f.name:data.get('calls',{}).get(f.name,f.default) for f in fields(sd.Calls)})
    calls.budget = BUDGET
    properties, configuration = screen.configuration()
    seeds, matches = inputs()
    data['configuration'] = dict(configuration,model=['EB','RLB'],beta_deg=[0.,45.,75.],d_theta=D,
                                 reuse_R0_beta_deg=5.,new_targets=list(NEW_TARGETS))
    data['matching_source_rows'] = matches
    data['elastic_source_rows'] = {k:s['source_row'] for k,s in seeds.items()}
    data['reused_rows'], data['reused_source_rows'] = reused_rows(seeds)
    for key in NEW_TARGETS:
        point = data['points'].setdefault(key,dict(attempts=[],complete=False))
        if point['complete']: continue
        seed = seeds[key]; cfg = config_for(properties,seed)
        for attempt in range(len(point['attempts'])+1,MAX_ATTEMPTS+1):
            started = time.perf_counter()
            z0 = complex(-seed['a_slope_pred']*D if attempt==1 else 0.,seed['Omega0'])
            record = dict(attempt=attempt,predictor=z0,d_theta=D,solver_path='reduced',eta=1)
            try:
                result = production.solve_mode(cfg,z0,eta=1,solver_path='reduced',
                    seed_states=seed['states'],elastic_z=1j*seed['Omega0'],calls=calls)
                record['result'] = {k:v for k,v in result.items() if k!='shape'}
                for name in completion.old.ARRAYS:
                    shapes[f'{key}_attempt{attempt}__{name}'] = result['shape'][name]
                control = dict(status='NOT_RUN',failures=['NOT_RUN'],newton_calls=0)
                if result['accepted']:
                    data['full_control_evaluations'] += 1
                    try:
                        control = completion.full_control(cfg,result['z'],seed,calls)
                        for name in completion.old.ARRAYS:
                            shapes[f'{key}_full__{name}'] = control['shape'][name]
                    except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:
                        control.update(status='CONTROL_UNRESOLVED',failures=[str(error)])
                    record['full_control'] = {k:v for k,v in control.items() if k!='shape'}
                    point['complete'] = True
                point['row'] = new_row(seed,result,control)
                point['accepted'] = result['accepted']
            except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:
                record['error'] = str(error)
                point['accepted'] = False
                point['row'] = dict(case_id=seed['case_id'],theory=seed['theory'],beta_deg=seed['beta_deg'],
                    matched_mode=seed['matched_mode'],eta=1,d_theta=D,source='NEW',
                    root_status='NUMERICAL_UNRESOLVED',notes=str(error))
            record['runtime_seconds'] = time.perf_counter()-started
            data['runtime_seconds'] += record['runtime_seconds']
            point['attempts'].append(record)
            save(data,shapes,calls)
            print(key,attempt,point['row']['root_status'],record.get('error',''),flush=True)
            if point['complete']: break
        point['complete'] = True  # exhausted failures stay qualified; no automatic extra retry
        save(data,shapes,calls)
    data['finished'] = True
    data['new_principal_roots'] = len(data['points'])
    data['new_roots_accepted'] = sum(p.get('accepted',False) for p in data['points'].values())
    data['reused_complex_roots'] = len(data['reused_rows'])
    data['root_attempts'] = {k:len(p['attempts']) for k,p in data['points'].items()}
    data['Newton_iterations'] = {k:[a['result']['correction']['steps'] if 'result' in a else None for a in p['attempts']] for k,p in data['points'].items()}
    data['stage_status'] = 'COMPLEX_CONFIRMATION_COMPLETED' if data['new_roots_accepted']==6 else 'PARTIAL_NUMERICAL_UNRESOLVED'
    data['protected_sources_unchanged'] = all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    assert data['protected_sources_unchanged']
    save(data,shapes,calls)
    data['output_hashes'] = {p:screen.sha(OUTPUT/p) for p in ('new_complex_roots.csv','representative_complex_comparison.csv','new_complex_shapes.npz')}
    write_json(checkpoint,data)
    write_json(OUTPUT/'run_manifest.json',{k:v for k,v in data.items() if k not in ('points','all_states','reused_rows')})
    return dict(stage_status=data['stage_status'],new_roots=data['new_roots_accepted'],calls=calls.snapshot())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute',action='store_true',required=True)
    parser.parse_args()
    print(json.dumps(compute(),indent=2))
