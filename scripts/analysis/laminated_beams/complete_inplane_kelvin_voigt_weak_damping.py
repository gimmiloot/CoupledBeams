"""D18: only A/C .001 -> .005 using frozen K18 production; four rows read-only.

Separate result/provenance contract from historical K16: full control is
evaluation-only. This orchestration adds no equations or root algorithm.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import fields
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as production
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_targeted_weak_damping as old
from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_solver_architecture as architecture

screen=old.screen
OUTPUT=ROOT/'results/laminated_beams/inplane_kelvin_voigt_targeted_weak_damping_completion'
TARGETS=(('A_STRONG',0.,'sorted_05',.005),('C_WEAK_ACTIVE',75.,'sorted_05',.005))
MAX_ATTEMPTS=2
REDUCED_BUDGET=300


class CompletionCalls(sd.Calls):
    """Count only; preserve solver numerics and cap total reduced B+B_z."""
    def matrix(self,derivative):
        # FullProvider increments full_B before this hook; HalfProvider after it.
        is_half=self.B==self.full_B+self.half_B
        if is_half and self.half_B+self.half_B_z+1+int(derivative)>REDUCED_BUDGET:
            raise RuntimeError('REDUCED_COST_LIMIT')
        super().matrix(derivative)


def initial_state():
    path=OUTPUT/'initial_state.json'
    if path.exists():return screen.read_json(path)
    folders=(old.SOURCE,old.OUTPUT,architecture.K17,architecture.OUTPUT)
    paths=[p for folder in folders for p in folder.iterdir() if p.is_file()]
    paths += [ROOT/p for p in (
        'scripts/lib/inplane_kelvin_voigt.py','scripts/lib/inplane_kelvin_voigt_solver.py',
        'scripts/lib/inplane_kelvin_voigt_symmetry_diagnostics.py',
        'scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_targeted_weak_damping.py',
        'scripts/analysis/laminated_beams/check_inplane_kelvin_voigt_solver_architecture.py',
        'docs/laminated_beams/inplane_kelvin_voigt_joint_theory.tex',
        'docs/laminated_beams/inplane_kelvin_voigt_targeted_weak_damping.md',
        'docs/laminated_beams/inplane_kelvin_voigt_ac_diagnostics.md',
        'docs/laminated_beams/inplane_kelvin_voigt_solver_architecture.md')]
    value=dict(initial_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        initial_git_status=subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
        protected_sources={p.relative_to(ROOT).as_posix():screen.sha(p) for p in paths})
    OUTPUT.mkdir(parents=True,exist_ok=True);old.json_write(path,value)
    return value


def seeds_and_sources():
    seeds=old.select_seeds(screen.read_csv(old.SOURCE/'elastic_screening.csv'))
    prior=screen.read_json(architecture.OUTPUT/'diagnostics.json')
    assert prior['reduced_status']=='SYMMETRY_REDUCED_PRODUCTION_PASS'
    assert prior['criteria']==kv.CRITERIA and prior['routing_version']==production.ROUTING_VERSION
    with np.load(old.SOURCE/'screening_shapes.npz',allow_pickle=False) as shapes:
        for seed in seeds:seed['states']=shapes[seed['shape_key']+'__states']
    return seeds,prior


def config_for(properties,seed,d):
    if (seed['state_id'],seed['beta_deg'],seed['elastic_sorted_mode'],d) not in TARGETS:
        raise ValueError('only A/C at .005; no other beta, d or mode')
    arm=kv.Arm.reduced('EB',properties)
    return production.Config((arm,arm),np.deg2rad(seed['beta_deg']),1.,d,mu=0.)


def predictor(seed,prior,attempt):
    key=seed['state_id'][0]+'_d001_auto'
    point=prior['points'][key];r=point['row']
    assert point['accepted'] and r['eta']==1 and r['d_theta']==.001
    assert r['solver_path']=='SYMMETRY_REDUCED'
    z=complex(r['z_re'],r['z_im'])
    if attempt==1:return z
    if attempt==2:return z-seed['a_slope_pred']*(.005-.001)
    raise ValueError('at most two attempts, no auxiliary d')


def full_control(cfg,z,seed,calls):
    """One full form at the given root, no corrector or solve_mode call."""
    before=calls.snapshot()
    provider=production.FullProvider(cfg,calls)
    balanced=sd.FrozenBalanced(provider,z)
    B,_=balanced.matrices(z)
    shape=kv.recover(provider,z,balanced.reactions(kv.right_null(B)))
    _,weights=kv.quadrature()
    seed_vector=kv.mass_vector(seed['states'],cfg.arms,weights)
    overlap=np.vdot(seed_vector,shape['vector'])
    if abs(overlap):
        phase=overlap.conjugate()/abs(overlap)
        for name in old.ARRAYS:shape[name]*=phase
    MAC=float(kv.mac_matrix([seed_vector],[shape['vector']])[0,0])
    diag=kv.diagnose(provider,z,shape)
    failures=kv.failures(diag,z,'ACTIVE',seed['Omega0'],MAC)
    return dict(z=z,diagnostics=diag,MAC=MAC,failures=failures,accepted=not failures,
                status='QUALIFIED' if failures else 'PASS',newton_calls=0,shape=shape,
                calls={k:calls.snapshot()[k]-before[k] for k in before if k not in ('limit','budget')})


def row_for(seed,d,result,control,origin,notes):
    r=result['row'] if 'row' in result else None
    z=complex(r['z_re'],r['z_im']) if r else complex(result['z'])
    diag=result['diagnostics'];full=control['diagnostics'];delta=diag['Delta_psi']
    if isinstance(delta,dict):delta=complex(delta['real'],delta['imag'])
    obs=old.observables(seed,d,z)
    return dict(state=seed['state_id'][0],state_id=seed['state_id'],
        sensitivity_class=dict(A='STRONG',B='INTERMEDIATE',C='WEAK ACTIVE')[seed['state_id'][0]],
        beta_deg=seed['beta_deg'],eta=1,d_theta=d,c_theta=d*kv.M_REF*kv.T_REF,
        data_origin=origin,solver_path=result['solver_path'],elastic_sorted_mode=seed['elastic_sorted_mode'],
        Omega_0=seed['Omega0'],elastic_Omega_0=seed['Omega0'],elastic_Lambda_0=seed['Lambda0'],
        **obs,a_slope_pred=seed['a_slope_pred'],zeta_slope_pred=seed['zeta_slope_pred'],
        relative_frequency_shift_over_d2=obs['frequency_shift_over_d2'],
        MAC=result['MAC'],M_phi=diag['M_phi'],K_phi=diag['K_phi'],C_phi=diag['C_phi'],
        Delta_psi_re=delta.real,Delta_psi_im=delta.imag,alpha_energy=diag['alpha_energy'],
        root_residual=diag['null_residual'],sigma_ratio=diag['sigma_ratio'],
        conjugate_residual=diag['conjugate_residual'],energy_residual=diag['r_E'],
        reduced_physical_residual=max(result['reduced_physical']['half_normalized']),
        lifted_full_physical_residual=max(diag['physical_residuals']),
        full_control_physical_residual=max(full['physical_residuals']),
        full_control_qualification=';'.join(control['failures']),
        root_status='ROOT_ACCEPTED' if result['accepted'] else 'NUMERICAL_UNRESOLVED',notes=notes)


def reused_rows(seeds,prior):
    rows=[]
    for seed in seeds:
        for d in ((.001,.005) if seed['state_id'].startswith('B_') else (.001,)):
            key=seed['state_id'][0]+('_d001' if d==.001 else '_d005')
            point=prior['points'][key+'_auto'];full=prior['points'][key+'_full']
            assert point['accepted']
            note=('K16 accepted; K18 reduced/full regression; read-only' if key.startswith('B') else
                  'K16 historical full recovery rejected; same root confirmed by K17/K18; read-only')
            rows.append(row_for(seed,d,point,full,'REUSED_K18:'+key+'_auto',note))
    return rows


def write_csv(path,rows):
    stream=io.StringIO(newline='');names=list(dict.fromkeys(k for r in rows for k in r))
    writer=csv.DictWriter(stream,fieldnames=names,lineterminator='\n')
    writer.writeheader();writer.writerows(rows);screen.atomic(path,stream.getvalue())


def save(data,shapes,calls):
    data['calls']=calls.snapshot()
    rows=data['reused_rows']+[p['row'] for p in data['points'].values() if 'row' in p]
    known={(r['state_id'],r['d_theta']) for r in rows}
    for sid,beta,mode,d in TARGETS:
        if (sid,d) not in known:
            rows.append(dict(state=sid[0],state_id=sid,beta_deg=beta,eta=1,d_theta=d,
                root_status=data['points'].get(sid,{}).get('status','NOT_RUN'),
                notes=data['points'].get(sid,{}).get('error','No accepted/computed value available')))
    rows.sort(key=lambda r:(r['state'],r['d_theta']))
    write_csv(OUTPUT/'combined_six_state_summary.csv',rows)
    write_csv(OUTPUT/'new_roots.csv',[r for r in rows if r['state_id'] in data['points'] and r['d_theta']==.005])
    with (OUTPUT/'new_complex_shapes.npz.tmp').open('wb') as f:np.savez_compressed(f,**shapes)
    os.replace(OUTPUT/'new_complex_shapes.npz.tmp',OUTPUT/'new_complex_shapes.npz')
    old.json_write(OUTPUT/'diagnostics.json',data)


def compute():
    initial=initial_state()
    assert all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    checkpoint=OUTPUT/'diagnostics.json'
    if checkpoint.exists():
        data=screen.read_json(checkpoint)
        assert data['criteria']==kv.CRITERIA
        if data.get('finished'):
            assert all(screen.sha(OUTPUT/p)==h for p,h in data['output_hashes'].items())
            return dict(missing_only=True,new_roots=0,B=0,B_z=0,form_recoveries=0,full_control_evaluations=0)
        shapes=dict(np.load(OUTPUT/'new_complex_shapes.npz',allow_pickle=False))
    else:
        data=dict(**initial,criteria=kv.CRITERIA,routing_version=production.ROUTING_VERSION,
            targets=TARGETS,points={},auxiliary_roots=0,new_beta=0,RLB=0,solver_architecture_changes=0,
            high_precision=0,full_newton_calls=0,full_control_evaluations=0,runtime_seconds=0.,
            environment=dict(executable=sys.executable,python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
            max_attempts=MAX_ATTEMPTS,reduced_matrix_budget=REDUCED_BUDGET,tests='NOT_RUN_YET')
        shapes={}
    calls=CompletionCalls(**{f.name:data.get('calls',{}).get(f.name,f.default) for f in fields(sd.Calls)})
    calls.budget=2000
    properties,config=screen.configuration()
    assert config==screen.read_json(old.SOURCE/'run_manifest.json')['configuration']
    data['elastic_source_configuration']=config
    data['configuration']=dict(config,beta_deg=[0.,75.],d_theta=.005)
    seeds,prior=seeds_and_sources();data['reused_rows']=reused_rows(seeds,prior)
    for sid,beta,mode,d in TARGETS:
        if sid in data['points'] and data['points'][sid].get('complete'):continue
        seed=next(s for s in seeds if s['state_id']==sid);cfg=config_for(properties,seed,d)
        point=data['points'].setdefault(sid,dict(attempts=[],seed_source=f'K18:{sid[0]}_d001_auto',
            source_row=prior['points'][sid[0]+'_d001_auto']['row']))
        for attempt in range(len(point['attempts'])+1,MAX_ATTEMPTS+1):
            start=time.perf_counter();z0=predictor(seed,prior,attempt)
            record=dict(attempt=attempt,predictor=z0,d_theta=d,solver_path='reduced',eta=1)
            try:
                result=production.solve_mode(cfg,z0,eta=1,solver_path='reduced',
                    seed_states=seed['states'],elastic_z=1j*seed['Omega0'],calls=calls)
                record['result']={k:v for k,v in result.items() if k!='shape'}
                for name in old.ARRAYS:shapes[f'{sid}_attempt{attempt}__{name}']=result['shape'][name]
                point.update(accepted=result['accepted'],status='ROOT_ACCEPTED' if result['accepted'] else 'NUMERICAL_UNRESOLVED')
                control=dict(diagnostics={'physical_residuals':[None]},failures=['NOT_RUN'],accepted=False,status='NOT_RUN',newton_calls=0)
                if result['accepted']:
                    data['full_control_evaluations']+=1
                    try:
                        control=full_control(cfg,result['z'],seed,calls)
                        for name in old.ARRAYS:shapes[f'{sid}_full__{name}']=control['shape'][name]
                    except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:
                        control.update(failures=[str(error)],status='CONTROL_UNRESOLVED')
                    record['full_control']={k:v for k,v in control.items() if k!='shape'}
                    point['complete']=True
                point['row']=row_for(seed,d,result,control,'NEW_REDUCED_K18',
                    'Direct .001 to .005; source K16 rejection and K17/K18 confirmation preserved')
            except (RuntimeError,ValueError,np.linalg.LinAlgError) as error:
                record['error']=str(error);point.update(accepted=False,status='NUMERICAL_UNRESOLVED',error=str(error))
            record['runtime_seconds']=time.perf_counter()-start
            data['runtime_seconds']+=record['runtime_seconds'];point['attempts'].append(record)
            save(data,shapes,calls)
            print(sid,'attempt',attempt,point['status'],record.get('error',''),flush=True)
            if point.get('complete') or 'COST_LIMIT' in record.get('error',''):break
        if not point.get('complete'):
            point['complete']=True;data['stop_reason']=sid;break
    data['finished']=True
    data['new_principal_roots']=len(data['points'])
    data['new_roots_accepted']=sum(p.get('accepted',False) for p in data['points'].values())
    data['continuation_seeds']={sid:dict(source=p['seed_source'],row=p['source_row']) for sid,p in data['points'].items()}
    data['root_attempts']={sid:len(p['attempts']) for sid,p in data['points'].items()}
    data['Newton_iterations']={sid:[a['result']['correction']['steps'] if 'result' in a else None for a in p['attempts']] for sid,p in data['points'].items()}
    data['stage_status']='SIX_STATE_COMPARISON_COMPLETED' if data['new_roots_accepted']==2 else 'PARTIAL_NUMERICAL_UNRESOLVED'
    data['protected_sources_unchanged']=all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    assert data['protected_sources_unchanged']
    save(data,shapes,calls)
    data['output_hashes']={p:screen.sha(OUTPUT/p) for p in ('new_roots.csv','combined_six_state_summary.csv','new_complex_shapes.npz')}
    old.json_write(checkpoint,data)
    old.json_write(OUTPUT/'run_manifest.json',{k:v for k,v in data.items() if k not in ('points','reused_rows')})
    return dict(stage_status=data['stage_status'],new_roots=data['new_roots_accepted'],calls=calls.snapshot())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute',required=True,action='store_true',help='only A/C .005, missing-only')
    parser.parse_args();print(json.dumps(compute(),indent=2))
