"""D17 fixed architecture regression set; no parameter sweep or new targets."""
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
sys.path.insert(0,str(ROOT)); sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as solver
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as symmetry
from scripts.analysis.laminated_beams import screen_inplane_kelvin_voigt_elastic as screen
from scripts.analysis.laminated_beams.pilot_inplane_kelvin_voigt import clean

BASE = ROOT/'results/laminated_beams'
OUTPUT = BASE/'inplane_kelvin_voigt_solver_architecture'
K11 = BASE/'inplane_spring_robustness'
K12 = BASE/'inplane_kelvin_voigt_pilot'
K15 = BASE/'inplane_kelvin_voigt_elastic_screening'
K16 = BASE/'inplane_kelvin_voigt_targeted_weak_damping'
K17 = BASE/'inplane_kelvin_voigt_ac_diagnostics'
CASE_IDS = ('K12_ACTIVE','K12_INACTIVE','B_d001','B_d005','A_d001','C_d001','K11_ASYMMETRIC')


def initial_state():
    path=OUTPUT/'initial_state.json'
    if path.exists():return screen.read_json(path)
    # Capture provenance before the first output; never reconstruct missing old data.
    names=subprocess.check_output(['git','diff','--name-only'],text=True).splitlines()
    names+=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],text=True).splitlines()
    protected=[p for folder in (K11,K12,K15,K16,K17) for p in folder.iterdir()
               if p.suffix in ('.csv','.json','.npz')]
    protected += [ROOT/'scripts/lib/inplane_kelvin_voigt.py',
        ROOT/'docs/laminated_beams/inplane_kelvin_voigt_joint_theory.tex',
        ROOT/'docs/laminated_beams/inplane_kelvin_voigt_ac_diagnostics.md']
    value=dict(initial_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        initial_git_status=subprocess.check_output(['git','status','--short'],text=True),
        initial_worktree_sha256={p:screen.sha(ROOT/p) for p in names},
        protected_sources={p.relative_to(ROOT).as_posix():screen.sha(p) for p in protected})
    OUTPUT.mkdir(parents=True,exist_ok=True);json_write(path,value)
    return value


def read_shape(path,key):
    with np.load(path,allow_pickle=False) as archive:
        return archive[key+'__states']


def inputs():
    """Read exactly seven controls. No data-based choice of easier replacements."""
    old12=screen.read_csv(K12/'modal_results.csv')
    old16=screen.read_json(K16/'diagnostics.json')
    old17=screen.read_json(K17/'diagnostics.json')
    cases=[]
    for role,eta in [('ACTIVE',1),('INACTIVE',-1)]:
        row=next(r for r in old12 if r['key']==f'EB_{role}_d1')
        seed=next(r for r in old12 if r['key']==f'EB_{role}_d0')
        assert row['status']==seed['status']=='CONFIRMED'
        cases.append(dict(case_id='K12_'+role,beta_deg=5.,mu=0.,d=float(row['d_theta']),eta=eta,
            reference=complex(float(row['z_real']),float(row['z_imag'])),
            elastic_z=1j*float(seed['Omega0']),
            seed_states=read_shape(K12/'shapes.npz',seed['key']),
            source_reference=f'{(K12/"modal_results.csv").relative_to(ROOT).as_posix()}#{row["key"]}',
            source_row=row, historical_reference=None,
            paths=['auto'] if eta==-1 else ['auto','full']))
    for case_id,key in [('B_d001','B_INTERMEDIATE_d1'),('B_d005','B_INTERMEDIATE_d2'),
                         ('A_d001','A_STRONG_d1'),('C_d001','C_WEAK_ACTIVE_d1')]:
        row=old16['points'][key]['row'];sid=row['state_id'];seed=old16['seeds'][sid]
        historical=complex(row['z_re'],row['z_im'])
        reference=historical
        source=f'{(K16/"targeted_weak_damping.csv").relative_to(ROOT).as_posix()}#{key}'
        if sid in old17['points']:
            r=old17['points'][sid]['row'];reference=complex(r['z_half_re'],r['z_half_im'])
            source=f'{(K17/"ac_symmetry_diagnostics.csv").relative_to(ROOT).as_posix()}#{sid}'
        cases.append(dict(case_id=case_id,beta_deg=row['beta_deg'],mu=0.,d=row['d_theta'],eta=1,
            reference=reference,elastic_z=1j*seed['Omega0'],
            seed_states=read_shape(K15/'screening_shapes.npz',seed['shape_key']),
            source_reference=source,source_row=row,historical_reference=historical,paths=['auto','full']))
    rows=screen.read_csv(K11/'verified_roots.csv')
    row=next(r for r in rows if r['shape_key']=='EB_m0.01_k1_b5_p01')
    assert row['root_status']=='CONFIRMED' and row['model']=='EB' and float(row['mu'])==.01
    cases.append(dict(case_id='K11_ASYMMETRIC',beta_deg=5.,mu=.01,d=0.,eta=None,
        reference=1j*float(row['Omega']),elastic_z=1j*float(row['Omega']),
        seed_states=read_shape(K11/'shapes.npz',row['shape_key']),
        source_reference=f'{(K11/"verified_roots.csv").relative_to(ROOT).as_posix()}#{row["shape_key"]}',
        source_row=row,historical_reference=None,paths=['auto']))
    assert tuple(c['case_id'] for c in cases)==CASE_IDS
    return cases


def configuration(properties,case):
    if case['case_id'] not in CASE_IDS:
        raise ValueError('outside the predetermined regression set')
    if case['case_id'].startswith(('A_','C_')) and case['d']!=.001:
        raise ValueError('A/C .005 and new d are not authorized')
    fixed_d=dict(B_d001=.001,B_d005=.005,K11_ASYMMETRIC=0.)
    if case['case_id'] in fixed_d and case['d']!=fixed_d[case['case_id']]:
        raise ValueError('new d is outside the regression set')
    layout=dict(K12_ACTIVE=(5.,0.),K12_INACTIVE=(5.,0.),B_d001=(45.,0.),B_d005=(45.,0.),
                A_d001=(0.,0.),C_d001=(75.,0.),K11_ASYMMETRIC=(5.,.01))
    if (case['beta_deg'],case['mu'])!=layout[case['case_id']]:
        raise ValueError('new beta/mu is outside the regression set')
    arms=tuple(kv.Arm.reduced('EB',properties,L) for L in (1-case['mu'],1+case['mu']))
    return solver.Config(arms,np.deg2rad(case['beta_deg']),1.,case['d'],mu=case['mu'])


def json_write(path,value):
    screen.write_json(path,clean(value))


def save(data,shapes):
    stream=io.StringIO(newline='')
    rows=[p['row'] for p in data['points'].values()]
    if rows:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
    screen.atomic(OUTPUT/'solver_regression.csv',stream.getvalue())
    with (OUTPUT/'control_shapes.npz.tmp').open('wb') as f:np.savez_compressed(f,**shapes)
    os.replace(OUTPUT/'control_shapes.npz.tmp',OUTPUT/'control_shapes.npz')
    json_write(OUTPUT/'diagnostics.json',data)


def classify_stage(points):
    reduced=[p for p in points.values() if p['row']['solver_path']=='SYMMETRY_REDUCED']
    full=[p for p in points.values() if p['row']['solver_path']=='FULL_TWO_ARM']
    reduced_status=('SYMMETRY_REDUCED_PRODUCTION_PASS' if len(reduced)==6 and all(p['row']['accepted'] for p in reduced)
                    else 'SYMMETRY_REDUCED_REGRESSION_STOP')
    qualified=[]
    for p in full:
        r=p['row']
        if not r['accepted']:
            if r['case_id'] not in ('A_d001','C_d001') or not p['root_agreement']['accepted'] or p['root_failures']:
                return reduced_status,'FULL_TWO_ARM_REGRESSION_STOP'
            if set(p['failures'])-{'PHYSICAL_GATE','COMPATIBILITY_GATE','POSSIBLE_MULTIPLICITY'}:
                return reduced_status,'FULL_TWO_ARM_REGRESSION_STOP'
            qualified.append(r['case_id'])
    if len(full)!=6:return reduced_status,'FULL_TWO_ARM_INCOMPLETE'
    status='FULL_TWO_ARM_PRODUCTION_PASS_WITH_HIGH_MODE_QUALIFICATION' if qualified else 'FULL_TWO_ARM_PRODUCTION_PASS'
    if reduced_status!='SYMMETRY_REDUCED_PRODUCTION_PASS':status='FULL_TWO_ARM_REGRESSION_STOP'
    return reduced_status,status


def compute():
    initial=initial_state()
    assert all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    checkpoint=OUTPUT/'diagnostics.json'
    if checkpoint.exists():
        data=screen.read_json(checkpoint)
        assert data['routing_version']==solver.ROUTING_VERSION and data['criteria']==kv.CRITERIA
        if data.get('finished'):
            assert all(screen.sha(OUTPUT/p)==h for p,h in data['output_hashes'].items())
            return dict(missing_only=True,new_control_rows=0,matrix_calls=0,root_calls=0,form_recoveries=0)
        shapes=dict(np.load(OUTPUT/'control_shapes.npz',allow_pickle=False))
    else:
        data=dict(**initial,routing_version=solver.ROUTING_VERSION,criteria=kv.CRITERIA,
            root_agreement_criterion=solver.ROOT_AGREEMENT,points={},attempts_per_point=1,
            full_T_policy='direct expm; expm_frechet(compute_expm=False) only for derivative',
            reaction_mapping='same K17 frozen positive row/column equilibration; a=b/column_scales, r=a*reaction_units',
            full_recovery_policy='unchanged K12 two independent step-expm propagations; no symmetry projection',
            environment=dict(executable=sys.executable,python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
            A_d005='NOT_RUN',C_d005='NOT_RUN',new_d=0,new_beta=0,RLB_roots=0,new_physical_study=0,
            tests='NOT_RUN_YET')
        shapes={}
    properties,config=screen.configuration();data['section_configuration']=config
    cases=inputs();completed=0
    for case in cases:
        for path in case['paths']:
            key=case['case_id']+'_'+path
            if key in data['points']:continue
            start=time.perf_counter();cfg=configuration(properties,case)
            calls=symmetry.Calls(budget=2000)
            result=solver.solve_mode(cfg,case['reference'],eta=case['eta'],solver_path=path,
                seed_states=case['seed_states'],elastic_z=case['elastic_z'],calls=calls)
            compare=solver.agreement(result['p'],case['reference']/kv.T_REF)
            diag=result['diagnostics'];full=result['full_diagnostics'];p=result['p'];z=result['z']
            accepted=result['accepted'] and compare['accepted']
            row=dict(case_id=case['case_id'],symmetry_status=result['symmetry_status'],solver_path=result['solver_path'],
                eta=result['eta'],beta_deg=case['beta_deg'],mu=case['mu'],d_theta=case['d'],source_reference=case['source_reference'],
                p_re=p.real,p_im=p.imag,z_re=z.real,z_im=z.imag,reference_p_re=case['reference'].real/kv.T_REF,
                reference_p_im=case['reference'].imag/kv.T_REF,root_abs_diff=compare['absolute'],root_rel_diff=compare['relative'],
                root_residual=diag['null_residual'],sigma_ratio=diag['sigma_ratio'],next_sigma_ratio=diag['next_sigma_ratio'],
                physical_residual=max(diag['physical_residuals']),energy_residual=diag['r_E'],MAC=result['MAC'],
                lifted_full_physical_residual=max(full['physical_residuals']) if result['symmetry_reduced'] else None,
                reduced_physical_residual=max(result['reduced_physical']['half_normalized']) if result['reduced_physical'] else None,
                accepted=accepted,qualification=';'.join(result['failures']),activity_status=result['activity_status'],
                root_equation_status=result['root_equation_status'],form_recovery_status=result['form_recovery_status'],
                rank_status=result['rank_status'],symmetry_reduced=result['symmetry_reduced'],
                notes='Root and form gates are separate; raw full rank flag is retained where present')
            point={k:v for k,v in result.items() if k!='shape'}
            point.update(row=row,root_agreement=compare,source_row=case['source_row'],
                historical_K16_z=case['historical_reference'],runtime_seconds=time.perf_counter()-start)
            data['points'][key]=point
            for name in ('states','reactions','a','vector'):shapes[key+'__'+name]=result['shape'][name]
            save(data,shapes);completed+=1
            print(key,accepted,row['qualification'],'physical',row['physical_residual'],'energy',row['energy_residual'],flush=True)
            # One bounded pass; stop if standard/reduced controls regress. No retry.
            local_qualification=(path=='full' and case['case_id'] in ('A_d001','C_d001')
                and compare['accepted'] and not result['root_failures']
                and not set(result['failures'])-{'PHYSICAL_GATE','COMPATIBILITY_GATE','POSSIBLE_MULTIPLICITY'})
            if not accepted and not local_qualification:
                data['stop_reason']=key;break
        if data.get('stop_reason'):break
    data['reduced_status'],data['full_status']=classify_stage(data['points'])
    data['finished']=True
    data['runtime_seconds']=sum(p['runtime_seconds'] for p in data['points'].values())
    keys=['B','B_z','full_B','full_B_z','half_B','half_B_z','expm','frechet','shape_expm','analytic_transfer','recoveries','corrections','total_build_equivalents']
    data['calls']={k:sum(p['calls'][k] for p in data['points'].values()) for k in keys}
    data['complex_newton_calls']=sum(p['complex_newton_calls'] for p in data['points'].values())
    data['regression_controls']=[dict(case_id=c['case_id'],beta_deg=c['beta_deg'],mu=c['mu'],
        d_theta=c['d'],eta=c['eta'],paths=c['paths'],source_reference=c['source_reference']) for c in cases]
    data['source_hashes_unchanged']=all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    assert data['source_hashes_unchanged']
    save(data,shapes)
    data['output_hashes']={p:screen.sha(OUTPUT/p) for p in ('solver_regression.csv','control_shapes.npz')}
    json_write(checkpoint,data)
    json_write(OUTPUT/'run_manifest.json',{k:v for k,v in data.items() if k!='points'})
    return dict(new_control_rows=completed,reduced_status=data['reduced_status'],full_status=data['full_status'],calls=data['calls'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute',action='store_true',required=True,help='fixed regression controls, missing-only')
    parser.parse_args()
    print(json.dumps(compute(),indent=2))
