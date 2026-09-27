"""D21: sparse real RLB screening and same-angle, read-only K15 EB comparison.

Different contract from K15: K21 reduced matrices/recovery and cross-theory
assignment, not its full EB assembly. Reuses the real detector, participation
formula, K11 common-coordinate matching and K12 physical diagnostics. No
complex corrector, positive-d provider or cross-beta tracking is called.
"""
from __future__ import annotations
import argparse
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
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as production
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.lib import inplane_spring_modes as mechanics
from scripts.analysis.laminated_beams import screen_inplane_kelvin_voigt_elastic as eb
from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_rlb_solver_architecture as arch

OUTPUT = eb.OUTPUT.parent/'inplane_kelvin_voigt_rlb_elastic_screening'
BETAS = eb.BETAS
VERSION = 'sparse-rlb-screening-v1'
CRITERIA = dict(eb.CRITERIA, comparison_MAC=.90, comparison_margin=.20,
                energy_residual=kv.CRITERIA['energy_residual'], max_seconds_per_group=120.)
METHOD = 'K11_COMMON_UW_MAC_GLOBAL_ASSIGNMENT_WITH_ETA_FILTER'


def initial_state():
    path=OUTPUT/'initial_state.json'
    if path.exists(): return eb.read_json(path)
    folders=(eb.SOURCE,eb.K12,eb.OUTPUT,arch.OUTPUT)
    files=[p for folder in folders for p in folder.iterdir() if p.suffix in ('.csv','.json','.npz')]
    files += [ROOT/p for p in ('scripts/lib/inplane_kelvin_voigt_solver.py',
        'scripts/lib/inplane_kelvin_voigt.py','scripts/lib/inplane_kelvin_voigt_symmetry_diagnostics.py',
        'scripts/lib/inplane_spring_modes.py','scripts/analysis/laminated_beams/screen_inplane_kelvin_voigt_elastic.py')]
    value=dict(initial_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        initial_git_status=subprocess.check_output(['git','status','--short'],text=True),
        protected_sources={p.relative_to(ROOT).as_posix():eb.sha(p) for p in files})
    OUTPUT.mkdir(parents=True,exist_ok=True);eb.write_json(path,value)
    return value


def real_array(value):
    value=np.asarray(value)
    if np.iscomplexobj(value):
        if np.any(value.imag != 0):
            raise ValueError('elastic real reduction requires identically zero imaginary part')
        value=value.real
    return value


class ElasticContext:
    def __init__(self, beta, properties):
        if beta not in BETAS: raise ValueError('no additional beta')
        arm=kv.Arm.reduced('RLB',properties)
        self.config=production.Config((arm,arm),math.radians(beta),1.,0.,mu=0.)
        assert production.route(self.config)=='SYMMETRY_REDUCED'
        self.calls=sd.Calls(budget=25000)
        self.full=production.FullProvider(self.config,self.calls)
        self.halves={eta:production.reduced_provider(self.full,eta) for eta in (1,-1)}
        self.started=time.perf_counter()

    def matrix(self, omega, eta):
        if np.iscomplexobj(omega) or not np.isfinite(omega) or omega<=0:
            raise ValueError('positive real elastic frequency only')
        if self.calls.B+self.calls.B_z>=CRITERIA['max_point_matrices']:
            raise RuntimeError('MATRIX_BUDGET')
        if time.perf_counter()-self.started>CRITERIA['max_seconds_per_group']:
            raise RuntimeError('GROUP_TIME_BUDGET')
        return real_array(self.halves[eta].matrices(1j*omega*kv.T_REF)[0])


class ReducedDiagnosticView:
    """K12 physical/energy diagnostics, using only the primary half matrix.

Full scalar joint residuals still act on both lifted arms. This view does
not assemble a full B and never invokes a root corrector.
"""
    def __init__(self, half): self.half=half
    def __getattr__(self, name): return getattr(self.half.full,name)
    def matrices(self, z, *, derivative=False):
        return self.half.matrices(z,derivative=derivative)


def saved_pool(beta, source):
    point=source['points'].get(f'RLB_m0_k1_b{beta:g}')
    if point is None: return [],None
    assert (point['model'],point['mu'],point['kappa'])==('RLB',0.,1.)
    assert set(point['errors']) <= {'MISSING_TARGET_OR_GUARD','GUARD_QUALIFIED'}
    rows=point['roots']
    assert all(a['Omega']<b['Omega'] for a,b in zip(rows,rows[1:]))
    for i,row in enumerate(rows,1):
        assert row['current_sorted_position']==i and row['root_status']=='CONFIRMED'
        assert not row['failures'] and row['multiplicity']==row['detected_nullity']==1
        assert abs(row['Lambda']**2-row['Omega'])<1e-12*row['Omega']
        assert abs(row['omega']*kv.T_REF-row['Omega'])<1e-12*row['Omega']
    if point['status']!='CONFIRMED':
        assert beta==45. and len(rows)==2
        assert all(c['accepted'] for scan in point['search'] for c in scan['candidates'])
    csv_rows=[r for r in eb.read_csv(eb.SOURCE/'verified_roots.csv') if r['point_id']==point['point_id']]
    assert [float(r['Omega']) for r in csv_rows]==[r['Omega'] for r in rows]
    return [dict(r,reuse_status='REUSED_ROOT_AND_FORM',source_provenance=
        f'{(eb.SOURCE/"verified_roots.csv").relative_to(ROOT).as_posix()}#{r["shape_key"]}') for r in rows[:7]],point


def search_missing(context, beta, pool, previous, source, audit):
    near=min((p for p in source['points'].values() if p['model']=='RLB' and p['mu']==0
              and p['kappa']==1 and p['status']=='CONFIRMED' and len(p['roots'])>=7),
             key=lambda p:abs(p['beta_deg']-beta))
    guesses=[r['Omega'] for r in near['roots'][:7]]
    upper=max(guesses)*1.02+.1
    lower=previous['search_upper']-.1 if previous else .01
    audit.update(lower=lower,upper=upper,predictor_beta=near['beta_deg'],retries=0,scans=[])
    additions=[]; suspects=[]
    for eta in (1,-1):
        provider=lambda omega:context.matrix(omega,eta)
        candidates=[]
        for lo,hi,count in eb.intervals(guesses,upper):
            if hi<=lower: continue
            candidates.extend(eb.pilot.roots._scan_candidates(provider,kv.T_REF,
                eb.pilot.policy(max(lo,lower),hi),case_id=f'RLB_screen_b{beta:g}',
                builder_id=f'K21_RLB_eta{eta}',scan_id='MISSING_ONLY',points=count,phases=(0.,))[0])
        candidates,proof=eb.pilot.reconcile_local_detections(candidates,provider)
        accepted,ambiguous=eb.pilot.consolidate(candidates)
        audit['scans'].append(dict(eta=eta,candidates=[eb.pilot.candidate_record(c) for c in candidates],reconciliations=proof))
        if ambiguous: raise ValueError('NUMERICAL_UNRESOLVED_DETECTION_CLUSTER')
        suspects.extend(c for c in candidates if eb.pilot.suspicious(c))
        for c in accepted:
            if c.diagnostics.detected_nullity!=1: raise ValueError('NUMERICAL_UNRESOLVED_MULTIPLICITY')
            additions.append(dict(Omega=c.omega_bar,omega=c.omega_bar/kv.T_REF,
                symmetry_class=eta,reuse_status='NEW_ELASTIC_ROOT',
                source_provenance='K21_REDUCED_REAL_CHARACTERISTIC_ROOT'))
    combined=sorted(pool+additions,key=lambda r:r['Omega'])
    audit['new_roots_including_guard']=len(additions)
    if len(combined)<7: raise ValueError('NUMERICAL_UNRESOLVED_MISSING_TARGET_OR_GUARD')
    unresolved=[eb.pilot.candidate_record(c) for c in suspects if c.interval_left_bar<=combined[5]['Omega']]
    audit['unresolved_below_target']=unresolved
    if unresolved: raise ValueError('NUMERICAL_UNRESOLVED_BELOW_TARGET')
    return combined[:7]


def shape_at(context, root, source_shapes):
    eta=root['symmetry_class']; z=1j*root['Omega']; half=context.halves[eta]
    if root['reuse_status']=='REUSED_ROOT_AND_FORM':
        key=root['shape_key']
        states=real_array(source_shapes[key+'__states'])
        reactions=real_array(source_shapes[key+'__reactions']).ravel()
        _,w=kv.quadrature(states.shape[1])
        shape=dict(states=states,reactions=reactions,a=reactions/context.full.reaction_scales,
                   vector=kv.mass_vector(states,context.config.arms,w))
    else:
        frozen=sd.FrozenBalanced(half,z)
        B,_=frozen.matrices(z)
        shape=sd.recover_half(half,z,frozen.reactions(kv.right_null(real_array(B))))
        for key in ('states','reactions','a','vector'): shape[key]=real_array(shape[key])
    diag=kv.diagnose(ReducedDiagnosticView(half),z,dict(shape,a=shape['a'][:3]))
    physical=sd.physical_details(half,z,shape)
    gates=kv.failures(diag,z,'INACTIVE' if eta==-1 else 'ACTIVE',root['Omega'],1.)
    if diag['symmetry_class']!=eta: gates.append('SYMMETRY_CLASS_MISMATCH')
    if max(physical['half_normalized'])>kv.CRITERIA['physical_residual']: gates.append('REDUCED_PHYSICAL_GATE')
    if abs(diag['M_phi']-1)>CRITERIA['mass_rtol']: gates.append('MASS_NORMALIZATION')
    return shape,diag,physical,gates


def modal_row(beta, index, root, shape, diag, physical, gates):
    eta=root['symmetry_class']; y=real_array(shape['states']); mass=diag['M_phi']
    delta=float(y[0,-1,2]-y[1,-1,2])
    # Structural zero is supported by the full-field reflection, not a tolerance.
    if eta==-1:
        assert np.array_equal(y[1],-y[0]@sd.F) and delta==0.
    status=eb.classify(eta,diag['symmetry_defect'],identical_arms=True,root_confirmed=not gates)
    if eta==1 and delta==0.: status='ZERO_PARTICIPATION_UNRESOLVED'
    indicators=eb.participation(delta,root['omega'],mass)
    notes=''
    if eta==1 and indicators['s_joint']<CRITERIA['small_participation_descriptive']:
        notes='SMALL_PARTICIPATION; descriptive K15 threshold only'
    return dict(model='RLB',beta_deg=beta,kappa_theta=1.,d_theta=0.,rlb_sorted_mode=f'rlb_sorted_{index:02d}',
        eta=eta,omega=root['omega'],Omega=root['Omega'],Lambda=math.sqrt(root['Omega']),mass_M=mass,
        psi1_joint=float(y[0,-1,2]),psi2_joint=float(y[1,-1,2]),Delta_psi=delta,abs_Delta_psi=abs(delta),
        P_joint=abs(delta)**2/mass,**indicators,G=indicators['zeta_slope_pred'],activity_status=status,
        symmetry_residual=diag['symmetry_defect'],root_residual=diag['null_residual'],sigma_ratio=diag['sigma_ratio'],
        next_sigma_ratio=diag['next_sigma_ratio'],physical_residual=max(diag['physical_residuals']),
        reduced_physical_residual=max(physical['half_normalized']),energy_residual=diag['r_E'],
        source_provenance=root['source_provenance'],reuse_status=root['reuse_status'],
        shape_key=f'b{beta:g}_rlb_sorted_{index:02d}',root_status='CONFIRMED',notes=notes)


def match_pairs(left, right):
    indices, mac, margin=mechanics.match(left,right,restrict_symmetry=True)
    return [(int(j),float(m),float(g),'CONFIRMED' if m>=CRITERIA['comparison_MAC']
             and g>=CRITERIA['comparison_margin'] else 'MATCH_AMBIGUOUS')
            for j,m,g in zip(indices,mac,margin)]


def compare(rows, shapes):
    reference=eb.read_csv(eb.OUTPUT/'elastic_screening.csv'); out=[]
    _,w=kv.quadrature()
    with np.load(eb.OUTPUT/'screening_shapes.npz',allow_pickle=False) as archive:
        for beta in BETAS:
            rlb=[r for r in rows if r['beta_deg']==beta]
            if len(rlb)!=6: continue
            refs=[r for r in reference if float(r['beta_deg'])==beta]
            assert len(refs)==6 and all(r['model']=='EB' and float(r['d_theta'])==0 and float(r['kappa_theta'])==1 for r in refs)
            left=[dict(symmetry_class=int(r['symmetry_eta']),common=mechanics.common_vector(archive[r['shape_key']+'__states'],w)) for r in refs]
            right=[dict(symmetry_class=r['eta'],common=mechanics.common_vector(shapes[r['shape_key']+'__states'],w)) for r in rlb]
            for ref,(j,mac,margin,status) in zip(refs,match_pairs(left,right)):
                row=rlb[j]; eta=row['eta']; assert int(ref['symmetry_eta'])==eta
                values=dict(beta_deg=beta,eta=eta,eb_sorted_mode=ref['sorted_mode'],rlb_sorted_mode=row['rlb_sorted_mode'],
                    match_method=METHOD,match_status=status,MAC_common=mac,competing_margin=margin,
                    Omega_EB=float(ref['Omega']),Omega_RLB=row['Omega'],P_joint_EB=float(ref['abs_Delta_psi'])**2/float(ref['mass_M']),
                    P_joint_RLB=row['P_joint'],s_joint_EB=float(ref['s_joint']),s_joint_RLB=row['s_joint'],
                    G_EB=float(ref['zeta_slope_pred']),G_RLB=row['G'],
                    a_slope_pred_EB=float(ref['a_slope_pred']),a_slope_pred_RLB=row['a_slope_pred'],
                    activity_EB=ref['activity_status'],activity_RLB=row['activity_status'],
                    delta_Omega=None,abs_delta_Omega=None,Omega_ratio=None,delta_G=None,abs_delta_G=None,G_ratio=None,
                    notes='same-angle comparison only; no cross-beta identity')
                if status=='CONFIRMED':
                    dw=(row['Omega']-values['Omega_EB'])/values['Omega_EB']
                    values.update(delta_Omega=dw,abs_delta_Omega=abs(dw),Omega_ratio=row['Omega']/values['Omega_EB'])
                    if eta==-1:
                        assert values['G_EB']==values['G_RLB']==0.
                        values['comparison_status']='INACTIVE_BOTH_THEORIES'
                    elif values['G_EB']!=0 and row['activity_status']=='ACTIVE':
                        dg=(row['G']-values['G_EB'])/values['G_EB']
                        values.update(delta_G=dg,abs_delta_G=abs(dg),G_ratio=row['G']/values['G_EB'],comparison_status='ACTIVE_MATCHED')
                    else: values['comparison_status']='ZERO_PARTICIPATION_UNRESOLVED'
                out.append(values)
    return out


def k12_control(rows):
    sources=eb.read_csv(eb.K12/'modal_results.csv'); result={}
    for role in ('ACTIVE','INACTIVE'):
        zero=next(r for r in sources if r['key']==f'RLB_{role}_d0')
        small=min((r for r in sources if r['model']=='RLB' and r['role']==role and float(r['d_theta'])>0),key=lambda r:float(r['d_theta']))
        current=[r for r in rows if r['beta_deg']==5 and r['Omega']==float(zero['Omega0'])]
        if len(current)!=1: result[role]=dict(status='UNRESOLVED_ELASTIC_MATCH'); continue
        row=current[0]; d=float(small['d_theta']); a=float(small['a_decay'])/d
        p=complex(float(small['p_real']),float(small['p_imag'])); zeta=-p.real/abs(p)/d
        item=dict(source_key=small['key'],d=d,rlb_sorted_mode=row['rlb_sorted_mode'],
            observed_a_over_d=a,observed_zeta_over_d=zeta,predicted_a=row['a_slope_pred'],predicted_G=row['G'],criterion=CRITERIA['K12_slope_relative'])
        if role=='ACTIVE':
            ea=abs(a/row['a_slope_pred']-1); eg=abs(zeta/row['G']-1)
            item.update(relative_a_difference=ea,relative_G_difference=eg,status='PASS' if max(ea,eg)<=item['criterion'] else 'FINITE_D_DISAGREEMENT')
        else:
            item.update(status='PASS' if row['activity_status']=='EXACT_INACTIVE_BY_SYMMETRY' and abs(float(small['a_decay']))<kv.CRITERIA['a_atol'] else 'FAIL',
                        historical_Delta_psi=[float(small['Delta_psi_real']),float(small['Delta_psi_imag'])])
        result[role]=item
    return result


def save(data, shapes):
    rows=[r for group in data['groups'].values() for r in group['rows']]
    arch.csv_write(OUTPUT/'rlb_elastic_screening.csv',rows)
    if shapes:
        with (OUTPUT/'rlb_elastic_shapes.npz.tmp').open('wb') as stream: np.savez_compressed(stream,**shapes)
        os.replace(OUTPUT/'rlb_elastic_shapes.npz.tmp',OUTPUT/'rlb_elastic_shapes.npz')
    data['output_hashes']={p.name:eb.sha(p) for p in OUTPUT.iterdir() if p.suffix in ('.csv','.npz')}
    arch.old.json_write(OUTPUT/'diagnostics.json',data)
    arch.old.json_write(OUTPUT/'run_manifest.json',{k:v for k,v in data.items() if k not in ('groups',)})


def compute():
    initial=initial_state()
    assert all(eb.sha(ROOT/p)==h for p,h in initial['protected_sources'].items())
    checkpoint=OUTPUT/'diagnostics.json'
    if checkpoint.exists():
        data=eb.read_json(checkpoint)
        assert data['version']==VERSION and data['criteria']==CRITERIA
        if data.get('finished'):
            assert all(eb.sha(OUTPUT/name)==h for name,h in data['output_hashes'].items())
            return dict(missing_only=True,new_states=0,new_root_groups=0,matrix_builds=0,form_recoveries=0,new_positive_d_complex_roots=0)
    else:
        data=dict(**initial,version=VERSION,criteria=CRITERIA,groups={},runtime_seconds=0.,
            environment=dict(executable=sys.executable,python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
            new_positive_d_complex_roots=0,complex_newton_calls=0,full_root_searches=0,
            cross_beta_tracking=0,new_d=0,beta_outside_selection=0,asymmetric_damping=0,FEM=0,
            matching_method=METHOD,tests='NOT_RUN_YET')
    start=time.perf_counter(); properties,config=eb.configuration(); config['model']='RLB'
    data['configuration']=config; data['source_heads']=eb.validate_sources(config)
    eb_config=eb.read_json(eb.OUTPUT/'run_manifest.json')['configuration']
    assert {k:v for k,v in config.items() if k!='model'}=={k:v for k,v in eb_config.items() if k!='model'}
    source=eb.read_json(eb.SOURCE/'diagnostics.json')
    shapes={}
    if (OUTPUT/'rlb_elastic_shapes.npz').exists():
        with np.load(OUTPUT/'rlb_elastic_shapes.npz') as archive: shapes=dict(archive)
    with np.load(eb.SOURCE/'shapes.npz',allow_pickle=False) as old_shapes:
        for beta in BETAS:
            key=str(beta)
            if key in data['groups']: continue
            context=ElasticContext(beta,properties); group=dict(beta_deg=beta,rows=[],diagnostics=[],search={},status='RUNNING')
            data['groups'][key]=group
            try:
                pool,previous=saved_pool(beta,source)
                group['reused_pool_size']=len(pool)
                if len(pool)<7: pool=search_missing(context,beta,pool,previous,source,group['search'])
                for index,root in enumerate(pool[:6],1):
                    shape,diag,physical,gates=shape_at(context,root,old_shapes)
                    group['diagnostics'].append(dict(Omega=root['Omega'],eta=root['symmetry_class'],diag=diag,physical=physical,failures=gates))
                    if gates: raise ValueError('NUMERICAL_UNRESOLVED:'+','.join(gates))
                    row=modal_row(beta,index,root,shape,diag,physical,gates);group['rows'].append(row)
                    for name in ('states','reactions'): shapes[row['shape_key']+'__'+name]=shape[name]
                    if beta==5 and index in (1,2):
                        B,_=context.full.matrices(1j*root['Omega'])
                        residual=float(np.linalg.norm(B@shape['a'])/(np.linalg.norm(B)*np.linalg.norm(shape['a'])))
                        assert residual<=kv.CRITERIA['null_residual']
                        group.setdefault('full_evaluation_controls',[]).append(dict(index=index,null_residual=residual))
                guard=pool[6];B=context.matrix(guard['omega'],guard['symmetry_class'])
                sigma=sd.spectrum(B)['ratios'][-1]
                group.update(guard=dict(Omega=guard['Omega'],eta=guard['symmetry_class'],sigma_ratio=sigma,
                    source_provenance=guard['source_provenance'],status='CONFIRMED' if sigma<=kv.CRITERIA['sigma_ratio'] else 'GUARD_QUALIFIED'),
                    status='CONFIRMED' if sigma<=kv.CRITERIA['sigma_ratio'] else 'TARGETS_CONFIRMED_GUARD_QUALIFIED')
            except (ValueError,RuntimeError,AssertionError) as error:
                group.update(status='NUMERICAL_UNRESOLVED',error=str(error))
            group['calls']=context.calls.snapshot();group['seconds']=time.perf_counter()-context.started
            save(data,shapes);print(beta,group['status'],len(group['rows']),group['calls']['B'],flush=True)
    rows=[r for g in data['groups'].values() for r in g['rows']]
    comparison=compare(rows,shapes);arch.csv_write(OUTPUT/'eb_rlb_matched_comparison.csv',comparison)
    data['K12_control']=k12_control(rows)
    data['reused_states']=sum(r['reuse_status']=='REUSED_ROOT_AND_FORM' for r in rows)
    data['new_states']=len(rows)-data['reused_states']
    data['new_root_groups']=sum(bool(g['search']) for g in data['groups'].values())
    data['new_roots_including_guard']=sum(g['search'].get('new_roots_including_guard',0) for g in data['groups'].values())
    data['matching_confirmed']=sum(r['match_status']=='CONFIRMED' for r in comparison)
    data['matching_ambiguous']=sum(r['match_status']!='CONFIRMED' for r in comparison)
    data['calls']={k:sum(g['calls'][k] for g in data['groups'].values()) for k in ('B','B_z','full_B','full_B_z','half_B','half_B_z','expm','frechet','shape_expm','recoveries','half_recoveries','corrections','total_build_equivalents')}
    data['status']='SCREENING_AND_COMPARISON_COMPLETED' if len(rows)==data['matching_confirmed']==24 else 'PARTIAL'
    data['finished']=True;data['runtime_seconds']+=time.perf_counter()-start
    data['final_git_status']=subprocess.check_output(['git','status','--short'],text=True)
    save(data,shapes)
    return {k:data[k] for k in ('status','reused_states','new_states','new_root_groups','matching_confirmed','matching_ambiguous','calls','runtime_seconds')}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute',action='store_true',required=True,help='bounded missing-only real screening')
    parser.parse_args();print(json.dumps(compute(),indent=2))
