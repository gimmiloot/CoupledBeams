"""D20 bounded RLB production validation; saved K12/K11 controls only.

Unlike D17's fixed EB routing audit, this contract gates roots on RLB block
and derivative identities and records an exact coefficient-limit comparison.
It reuses the production dispatcher, K17 mathematics and K19 evaluation-only
full control; it is not a second solver or a physical parameter sweep.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import asdict, fields, replace
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT/'src'))
import numpy as np
import scipy
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as solver
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.lib import inplane_spring_modes as elastic
from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_solver_architecture as old
from scripts.analysis.laminated_beams.complete_inplane_kelvin_voigt_weak_damping import full_control

screen = old.screen
OUTPUT = old.BASE/'inplane_kelvin_voigt_rlb_solver_architecture'
K07 = old.BASE/'inplane_rotational_spring_rlb_eb_limit'
K18 = old.OUTPUT
VERSION = 'rlb-production-validation-v1'
BUDGET = 1000  # B, B_z, expm and two equivalents per Frechet, including forms
CASE_IDS = tuple(f'K12_RLB_{role}_d{i}' for role in ('ACTIVE','INACTIVE') for i in (1,2,3)) + (
    'EXACT_EB_LIMIT', 'K11_RLB_ASYMMETRIC')


def initial_state():
    path = OUTPUT/'initial_state.json'
    if path.exists():
        return screen.read_json(path)
    protected = [p for folder in (K07,old.K11,old.K12,K18,
        old.BASE/'inplane_kelvin_voigt_targeted_weak_damping_completion')
        for p in folder.iterdir() if p.suffix in ('.csv','.json','.npz')]
    protected += [ROOT/'scripts/lib/inplane_kelvin_voigt.py',
        ROOT/'scripts/lib/inplane_kelvin_voigt_symmetry_diagnostics.py']
    protected += [ROOT/'docs/laminated_beams'/name for name in (
        'inplane_kelvin_voigt_joint_theory.tex','inplane_kelvin_voigt_pilot.md',
        'inplane_spring_robustness.md','inplane_rotational_spring_rlb_eb_limit.md',
        'inplane_kelvin_voigt_solver_architecture.md','inplane_kelvin_voigt_weak_damping_parity.md')]
    value = dict(initial_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        initial_git_status=subprocess.check_output(['git','status','--short'],text=True),
        protected_sources={p.relative_to(ROOT).as_posix():screen.sha(p) for p in protected})
    OUTPUT.mkdir(parents=True,exist_ok=True)
    old.json_write(path,value)
    return value


def inputs():
    rows = {r['key']:r for r in screen.read_csv(old.K12/'modal_results.csv')}
    seed_info = screen.read_json(old.K12/'seed_manifest.json')
    cases = []
    for role, eta in (('ACTIVE',1),('INACTIVE',-1)):
        seed = rows[f'RLB_{role}_d0']
        assert seed_info['seeds']['RLB_'+role]['symmetry_class'] == eta
        for i in (1,2,3):
            row = rows[f'RLB_{role}_d{i}']
            assert row['status'] == seed['status'] == 'CONFIRMED'
            assert float(row['beta0_deg']) == 5. and float(row['kappa_theta']) == 1.
            cases.append(dict(case_id='K12_'+row['key'],beta_deg=5.,mu=0.,eta=eta,
                d=float(row['d_theta']),reference=complex(float(row['z_real']),float(row['z_imag'])),
                reference_p=complex(float(row['p_real']),float(row['p_imag'])),
                elastic_z=1j*float(seed['Omega0']),
                seed_states=old.read_shape(old.K12/'shapes.npz',seed['key']),
                reference_states=old.read_shape(old.K12/'shapes.npz',row['key']),
                source_reference=f'{(old.K12/"modal_results.csv").relative_to(ROOT).as_posix()}#{row["key"]}',
                source_root_status=row['status'],source_row=row))
    # Structural limit uses one already validated EB production root, not a new d.
    saved = screen.read_json(K18/'diagnostics.json')['points']['K12_ACTIVE_auto']
    row = saved['row']; assert row['accepted'] and row['beta_deg'] == 5.
    seed = rows['EB_ACTIVE_d0']
    cases.append(dict(case_id='EXACT_EB_LIMIT',beta_deg=5.,mu=0.,eta=1,d=row['d_theta'],
        reference=complex(row['z_re'],row['z_im']),reference_p=complex(row['p_re'],row['p_im']),
        elastic_z=1j*float(seed['Omega0']),seed_states=old.read_shape(old.K12/'shapes.npz',seed['key']),
        reference_states=old.read_shape(K18/'control_shapes.npz','K12_ACTIVE_auto'),
        source_reference=f'{(K18/"solver_regression.csv").relative_to(ROOT).as_posix()}#K12_ACTIVE_auto',
        source_root_status='ACCEPTED',source_row=row))
    # First simple, confirmed low mode at the already used 5 degrees, away from interaction.
    key = 'RLB_m0.01_k1_b5_p01'
    row = next(r for r in screen.read_csv(old.K11/'verified_roots.csv') if r['shape_key'] == key)
    mapping = next(r for r in screen.read_csv(old.K11/'control_modes.csv') if r['shape_key'] == key)
    assert row['root_status'] == mapping['local_assignment_status'] == mapping['mapping_status'] == 'CONFIRMED'
    assert int(row['multiplicity']) == int(row['detected_nullity']) == 1
    assert float(mapping['MAC']) >= kv.CRITERIA['MAC'] and float(mapping['margin']) >= .20
    shape = old.read_shape(old.K11/'shapes.npz',key)
    cases.append(dict(case_id='K11_RLB_ASYMMETRIC',beta_deg=5.,mu=.01,eta=None,d=0.,
        reference=1j*float(row['Omega']),reference_p=1j*float(row['omega']),elastic_z=1j*float(row['Omega']),
        seed_states=shape,reference_states=shape,
        source_reference=f'{(old.K11/"verified_roots.csv").relative_to(ROOT).as_posix()}#{key}',
        source_root_status=row['root_status'],source_row=row,source_mapping=mapping))
    assert tuple(c['case_id'] for c in cases) == CASE_IDS
    return cases


def configuration(properties, case):
    if case['case_id'] not in CASE_IDS or case['beta_deg'] != 5.:
        raise ValueError('outside fixed saved controls; no beta screening')
    if case['case_id']=='K11_RLB_ASYMMETRIC':
        expected=(0.,.01,None)
    else:
        key=case['case_id'].removeprefix('K12_') if case['case_id']!='EXACT_EB_LIMIT' else 'EB_ACTIVE_d1'
        source=next(r for r in screen.read_csv(old.K12/'modal_results.csv') if r['key']==key)
        expected=(float(source['d_theta']),0.,-1 if source['role']=='INACTIVE' else 1)
    if (case['d'],case['mu'],case['eta']) != expected:
        raise ValueError('new d, mu or mode is outside the regression contract')
    arms = tuple(kv.Arm.reduced('RLB',properties,L) for L in (1-case['mu'],1+case['mu']))
    if case['case_id'] == 'EXACT_EB_LIMIT':
        arms = tuple(replace(a,invS=0.,J=0.) for a in arms)
    return solver.Config(arms,np.deg2rad(case['beta_deg']),1.,case['d'],mu=case['mu'])


def configuration_record(section_config, cases):
    """Keep shared laminate/reference data, not the EB screening's scope labels."""
    shared = {k:v for k,v in section_config.items()
              if k not in ('model','mu','L1','L2','d_theta','beta_deg')}
    return dict(shared, model='RLB', beta_deg=sorted({c['beta_deg'] for c in cases}),
        cases=[dict(case_id=c['case_id'],mu=c['mu'],L1=1-c['mu'],L2=1+c['mu'],
                    d_theta=c['d'],exact_EB_coefficient_limit=c['case_id']=='EXACT_EB_LIMIT')
               for c in cases])


def difference(name, actual, reference, criterion, *, exact=False):
    absolute = float(np.linalg.norm(np.asarray(actual)-np.asarray(reference)))
    norm = float(np.linalg.norm(reference))
    relative = absolute/norm if norm else None
    accepted = bool(np.array_equal(actual,reference) if exact else
                    absolute <= kv.CRITERIA['matrix_atol'] + criterion*norm)
    return dict(quantity=name,absolute_difference=absolute,relative_difference=relative,
        criterion='exact equality' if exact else f'atol=1e-12; rtol={criterion:g}',
        status='PASS' if accepted else 'FAIL',accepted=accepted)


def matrix_checks(properties, cases, calls):
    """Stage A and coefficient-limit matrices; must finish before any corrector."""
    checks, limits = [], []
    for case in cases[:3]:
        cfg = configuration(properties,case); f = solver.FullProvider(cfg,calls)
        z = case['reference']; p=z/kv.T_REF; arm=cfg.arms[0]
        H = kv.state_matrix(p,arm); Hp=kv.state_matrix(p,arm,derivative=True)
        check = sd.algebra_check(f,z)
        expected_Hp = np.zeros((6,6),complex)
        expected_Hp[3,0]=expected_Hp[4,1]=2*arm.m*p; expected_Hp[5,2]=2*arm.J*p
        check['H_p_exact'] = bool(np.array_equal(Hp,expected_Hp))
        old_H = elastic.Arm('RLB',properties,arm.L).matrix(p.imag)
        check['elastic_H'] = difference('H at i*omega',kv.state_matrix(1j*p.imag,arm),old_H,kv.CRITERIA['H_rtol'])
        T,Tz=f.transfer(z,arm,derivative=True)
        check['T_reflection'] = difference('F*T=T*F',sd.F@T,T@sd.F,kv.CRITERIA['H_rtol'])
        # Joint nullspace invariance in the existing dimensionless state/row units.
        joint=kv.joint_matrix(p,f.beta,f.k,f.c)*f.state_units[None,:]/f.row_units[:,None]
        swap=np.block([[np.zeros((6,6)),sd.F],[sd.F,np.zeros((6,6))]])
        null=np.linalg.svd(joint)[2].conj().T[:,6:]
        check['joint_reflection_relative']=float(np.linalg.norm(joint@swap@null)/np.linalg.norm(joint))
        minus=sd.conditions(p,f.beta,f.k,f.c,-1)
        check['inactive_independent']=bool(np.array_equal(minus,sd.conditions(p,f.beta,0.,0.,-1)))
        B,Bz=f.matrices(z,derivative=True)
        h=1e-4  # fixed K18 diagnostic step; no tuning or production finite difference
        fd=(f.matrices(z+h)[0]-f.matrices(z-h)[0])/(2*h)
        check['B_z_finite_difference']=difference('B_z',Bz,fd,kv.CRITERIA['derivative_rtol'])
        check['conjugacy']=difference('B conjugacy',f.matrices(z.conjugate())[0],B.conj(),kv.CRITERIA['H_rtol'])
        scale=arm.scale()[None,:]/arm.scale()[:,None]
        check['T_authoritative']=difference('direct T request independence',f.transfer(z,arm)*scale,T*scale,0.,exact=True)
        check['accepted'] = bool(check['accepted'] and check['H_p_exact'] and check['inactive_independent']
            and check['joint_reflection_relative'] <= kv.CRITERIA['H_rtol']
            and all(check[k]['accepted'] for k in ('elastic_H','T_reflection','B_z_finite_difference','conjugacy','T_authoritative')))
        check['case_id']=case['case_id']; checks.append(check)
        # RLB remains model='RLB'; only invS and J are set to zero, exactly.
        lim_arm=replace(arm,invS=0.,J=0.); eb_arm=kv.Arm.reduced('EB',properties)
        r=solver.FullProvider(replace(cfg,arms=(lim_arm,lim_arm)),calls)
        e=solver.FullProvider(replace(cfg,arms=(eb_arm,eb_arm)),calls)
        pairs=[('H',kv.state_matrix(p,lim_arm),kv.state_matrix(p,eb_arm),True),
               ('H_p',kv.state_matrix(p,lim_arm,derivative=True),kv.state_matrix(p,eb_arm,derivative=True),True)]
        rt,rtz=r.transfer(z,lim_arm,derivative=True); et,etz=e.transfer(z,eb_arm,derivative=True)
        pairs += [('T',rt*scale,et*scale,True),('T_p',rtz*kv.T_REF*scale,etz*kv.T_REF*scale,True)]
        for eta,label in ((1,'B_plus'),(-1,'B_minus')):
            rb,rbz=solver.reduced_provider(r,eta).matrices(z,derivative=True)
            eb,ebz=solver.reduced_provider(e,eta).matrices(z,derivative=True)
            pairs += [(label,rb,eb,False),(label+'_p',rbz*kv.T_REF,ebz*kv.T_REF,False)]
        rb,rbz=r.matrices(z,derivative=True); eb,ebz=e.matrices(z,derivative=True)
        pairs += [('B_full',rb,eb,True),('B_full_p',rbz*kv.T_REF,ebz*kv.T_REF,True)]
        limits.extend(dict(difference(name,a,b,kv.CRITERIA['transfer_rtol'],exact=exact),
                           case_id=case['case_id']) for name,a,b,exact in pairs)
    return dict(checks=checks,exact_limit=limits,
        accepted=all(c['accepted'] for c in checks+limits),calls=calls.snapshot())


def shape_comparison(shape, reference, arms):
    _,w=kv.quadrature(shape['states'].shape[1])
    ref=kv.mass_vector(reference,arms,w); v=shape['vector']
    overlap=np.vdot(ref,v)
    phase=overlap.conjugate()/abs(overlap) if abs(overlap) else 1.
    aligned=shape['states']*phase
    return dict(MAC=float(kv.mac_matrix([ref],[v])[0,0]),
        mass_norm_difference=float(np.linalg.norm(v*phase-ref)),
        delta_aligned=complex(aligned[0,-1,2]-aligned[1,-1,2]),
        reference_delta=complex(reference[0,-1,2]-reference[1,-1,2]))


def record(case, cfg, result):
    diag=result['diagnostics']; z=result['z']; p=z/kv.T_REF
    compare=solver.agreement(p,case['reference_p'])
    shapes=shape_comparison(result['shape'],case['reference_states'],cfg.arms)
    failures=list(result['failures'])
    if not compare['accepted']: failures.append('REFERENCE_ROOT_MISMATCH')
    if shapes['MAC']<kv.CRITERIA['MAC']: failures.append('REFERENCE_FORM_MISMATCH')
    reduced=result.get('reduced_physical')
    row=dict(case_id=case['case_id'],model='RLB',symmetry_status='EXACT_IDENTICAL' if cfg.identical else 'NON_IDENTICAL',
        solver_path=result['solver_path'],eta=case['eta'],beta_deg=case['beta_deg'],mu=case['mu'],d_theta=case['d'],
        source_reference=case['source_reference'],source_root_status=case['source_root_status'],
        p_re=p.real,p_im=p.imag,z_re=z.real,z_im=z.imag,
        reference_p_re=case['reference_p'].real,reference_p_im=case['reference_p'].imag,
        root_abs_diff=compare['absolute'],root_rel_diff=compare['relative'],
        root_residual=diag['null_residual'],sigma_ratio=diag['sigma_ratio'],next_sigma_ratio=diag['next_sigma_ratio'],
        physical_residual=max(diag['physical_residuals']),
        reduced_physical_residual=max(reduced['half_normalized']) if reduced else None,
        lifted_full_physical_residual=max(diag['physical_residuals']) if result['symmetry_reduced'] else None,
        energy_residual=diag['r_E'],conjugate_residual=diag['conjugate_residual'],
        MAC=result['MAC'],reference_MAC=shapes['MAC'],symmetry_defect=diag['symmetry_defect'],
        inactive_by_symmetry=result['inactive_by_symmetry'],accepted=not failures,qualification=';'.join(failures),
        alpha=-p.real,a=-z.real,Omega_d=z.imag,zeta=-p.real/abs(p),
        M_phi=diag['M_phi'],K_phi=diag['K_phi'],C_phi=diag['C_phi'],alpha_energy=diag['alpha_energy'],
        Delta_psi_re=diag['Delta_psi'].real,Delta_psi_im=diag['Delta_psi'].imag,
        complex_newton_calls=result['complex_newton_calls'],newton_updates=result.get('correction',{}).get('steps',0),
        root_origin=result['root_origin'],
        notes='saved-state regression; exact coefficient limit' if case['case_id']=='EXACT_EB_LIMIT' else
              'saved-state regression; no new physical parameters')
    point={k:v for k,v in result.items() if k!='shape'}
    point.update(row=row,root_agreement=compare,shape_comparison=shapes,source_row=case['source_row'])
    if 'alpha' in case['source_row']:
        ref=case['source_row']
        point['source_observable_differences']={k:row[k]-float(ref[source]) for k,source in (
            ('alpha','alpha'),('a','a_decay'),('Omega_d','Omega_d'),('M_phi','M_phi'),('K_phi','K_phi'),('C_phi','C_phi'))}
        point['source_observable_differences']['zeta']=row['zeta']-float(ref['alpha'])/abs(case['reference_p'])
        point['source_observable_differences']['Delta_psi_aligned']=shapes['delta_aligned']-shapes['reference_delta']
    return point


def csv_write(path, rows):
    names=list(dict.fromkeys(k for r in rows for k in r)); stream=io.StringIO(newline='')
    if rows:
        writer=csv.DictWriter(stream,fieldnames=names,lineterminator='\n');writer.writeheader();writer.writerows(rows)
    screen.atomic(path,stream.getvalue())


def save(data, shapes, calls):
    data['calls']=calls.snapshot()
    csv_write(OUTPUT/'rlb_solver_regression.csv',[p['row'] for p in data['points'].values()])
    csv_write(OUTPUT/'rlb_eb_exact_limit.csv',data.get('matrix_validation',{}).get('exact_limit',[])+data.get('limit_root_form',[]))
    if shapes:
        with (OUTPUT/'control_shapes.npz.tmp').open('wb') as f: np.savez_compressed(f,**shapes)
        os.replace(OUTPUT/'control_shapes.npz.tmp',OUTPUT/'control_shapes.npz')
    data['output_hashes']={p.name:screen.sha(p) for p in OUTPUT.iterdir() if p.name in (
        'rlb_solver_regression.csv','rlb_eb_exact_limit.csv','control_shapes.npz')}
    old.json_write(OUTPUT/'diagnostics.json',data)
    old.json_write(OUTPUT/'run_manifest.json',{k:v for k,v in data.items() if k not in ('points','matrix_validation')})


def compute(matrix_only=False):
    initial=initial_state()
    assert all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    path=OUTPUT/'diagnostics.json'
    if path.exists():
        data=screen.read_json(path)
        assert data['version']==VERSION and data['criteria']==kv.CRITERIA
        assert data['production_sha256']==screen.sha(ROOT/'scripts/lib/inplane_kelvin_voigt_solver.py')
        if data.get('finished') or (matrix_only and 'matrix_validation' in data):
            assert all(screen.sha(OUTPUT/p)==h for p,h in data['output_hashes'].items())
            return dict(missing_only=True,new_rows=0,root_calls=0,matrix_calls=0,form_recoveries=0,status=data['status'])
        calls=sd.Calls(**{k:v for k,v in data['calls'].items() if k in {f.name for f in fields(sd.Calls)}})
        shapes={}
        if (OUTPUT/'control_shapes.npz').exists():
            with np.load(OUTPUT/'control_shapes.npz',allow_pickle=False) as archive: shapes=dict(archive)
    else:
        data=dict(**initial,version=VERSION,routing_version=solver.RLB_ROUTING_VERSION,criteria=kv.CRITERIA,
            root_agreement_criterion=solver.ROOT_AGREEMENT,
            production_sha256=screen.sha(ROOT/'scripts/lib/inplane_kelvin_voigt_solver.py'),
            environment=dict(executable=sys.executable,python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
            points={},runtime_seconds=0.,rlb_active_corrector_calls=0,exact_limit_corrector_calls=0,
            asymmetric_elastic_corrector_calls=0,inactive_complex_newton=0,full_evaluation_controls=0,
            asymmetric_positive_d_roots=0,new_beta=0,new_d_outside_K12=0,new_physical_parameter_study=0,
            FEM=0,high_precision=0,tests='NOT_RUN_YET',budget=BUDGET)
        calls=sd.Calls(budget=BUDGET);shapes={}
    start=time.perf_counter()
    properties,section_config=screen.configuration()
    cases=inputs();data['configuration']=configuration_record(section_config,cases)
    data['controls']=[{k:c[k] for k in ('case_id','beta_deg','mu','d','eta','source_reference','source_root_status')} for c in cases]
    data['arm_configuration']=[asdict(a) for a in configuration(properties,cases[0]).arms]
    data['d_levels']=[c['d'] for c in cases[:3]]
    try:
        if 'matrix_validation' not in data:
            data['matrix_validation']=matrix_checks(properties,cases,calls)
            data['status']='MATRIX_PASS' if data['matrix_validation']['accepted'] else 'MATRIX_FAIL_STOP'
            data['finished']=not data['matrix_validation']['accepted']
            save(data,shapes,calls)
        if data.get('finished') or matrix_only:
            data['runtime_seconds']+=time.perf_counter()-start;save(data,shapes,calls)
            return dict(status=data['status'],calls=calls.snapshot())
        for case in cases:
            key=case['case_id'];cfg=configuration(properties,case)
            if key not in data['points']:
                if key.startswith('K12_RLB_ACTIVE'): data['rlb_active_corrector_calls']+=1
                elif key=='EXACT_EB_LIMIT': data['exact_limit_corrector_calls']+=1
                elif key=='K11_RLB_ASYMMETRIC': data['asymmetric_elastic_corrector_calls']+=1
                assert data['rlb_active_corrector_calls']<=3
                result=solver.solve_mode(cfg,case['reference'],eta=case['eta'],seed_states=case['seed_states'],
                    elastic_z=case['elastic_z'],calls=calls)
                data['points'][key]=record(case,cfg,result)
                for name in ('states','reactions','a','vector'): shapes[key+'__'+name]=result['shape'][name]
                if key=='EXACT_EB_LIMIT':
                    comp=data['points'][key]['shape_comparison']
                    data['limit_root_form']=[difference('one root',result['p'],case['reference_p'],solver.ROOT_AGREEMENT['relative']),
                        dict(quantity='one form',absolute_difference=comp['mass_norm_difference'],relative_difference=comp['mass_norm_difference'],
                             criterion='phase-aligned physical mass norm; rtol=1e-9',status='PASS' if comp['mass_norm_difference']<=kv.CRITERIA['transfer_rtol'] else 'FAIL',
                             accepted=comp['mass_norm_difference']<=kv.CRITERIA['transfer_rtol'])]
                save(data,shapes,calls)
                print(key,data['points'][key]['row']['accepted'],data['points'][key]['row']['qualification'],flush=True)
                if not data['points'][key]['row']['accepted']:
                    data['stop_reason']=key;break
            if key.startswith('K12_RLB_ACTIVE') and key+'_full_control' not in data['points']:
                z=complex(data['points'][key]['row']['z_re'],data['points'][key]['row']['z_im'])
                control=full_control(cfg,z,dict(states=case['seed_states'],Omega0=case['elastic_z'].imag),calls)
                control.update(solver_path='FULL_TWO_ARM',symmetry_reduced=False,inactive_by_symmetry=False,
                    complex_newton_calls=0,root_origin='EVALUATION_ONLY_AT_REDUCED_ROOT')
                full_case=dict(case,case_id=key+'_full_control')
                point=record(full_case,cfg,control)
                reduced_states=shapes[key+'__states']
                point['full_reduced_form']=shape_comparison(control['shape'],reduced_states,cfg.arms)
                point['full_reduced_energy']={n:control['diagnostics'][n]-data['points'][key]['diagnostics'][n] for n in ('M_phi','K_phi','C_phi','alpha_energy')}
                data['points'][full_case['case_id']]=point;data['full_evaluation_controls']+=1
                for name in ('states','reactions','a','vector'): shapes[full_case['case_id']+'__'+name]=control['shape'][name]
                save(data,shapes,calls)
                if not point['row']['accepted']:
                    data['stop_reason']=full_case['case_id'];break
        complete=len(data['points'])==11 and all(p['row']['accepted'] for p in data['points'].values())
        limit_pass=all(r['accepted'] for r in data.get('limit_root_form',[])) and len(data.get('limit_root_form',[]))==2
        data['status']='RLB_KV_PRODUCTION_PASS' if complete and limit_pass else 'RLB_KV_PRODUCTION_STOP'
        data['finished']=True
    except (RuntimeError,ValueError,np.linalg.LinAlgError) as exc:
        data.update(status='RLB_KV_PRODUCTION_STOP',finished=True,stop_reason=str(exc))
    data['runtime_seconds']+=time.perf_counter()-start
    data['source_hashes_unchanged']=all(screen.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    assert data['source_hashes_unchanged']
    data['final_git_status']=subprocess.check_output(['git','status','--short'],text=True)
    save(data,shapes,calls)
    return dict(status=data['status'],rows=len(data['points']),calls=calls.snapshot(),runtime_seconds=data['runtime_seconds'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute',action='store_true',required=True,help='fixed missing-only validation')
    parser.add_argument('--matrix-only',action='store_true',help='stop after matrices, before any root call')
    args=parser.parse_args()
    print(json.dumps(compute(matrix_only=args.matrix_only),indent=2))
