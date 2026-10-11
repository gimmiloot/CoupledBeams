"""Scoped physical-state continuation of the saved C modal experiment.

The only entry point is the existing coupled pilot CLI. No action, spatial
operator, initial subspace, or solver tolerance is defined anew here.
"""
from __future__ import annotations
import argparse, csv, cProfile, hashlib, io, json, os, pstats, shutil, sys, time
from pathlib import Path
from contextlib import contextmanager
import numpy as np
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_coupled_sensitivity_analysis as a
from scripts.lib import nlsp_construction_c_sensitivity as c
from scripts.lib import nlsp_coupled_sensitivity_resume as r

CONFIG=p.ROOT/'data/input/nlsp_construction_C_longer_sensitivity.json'
OUTPUT=p.ROOT/'results/nlsp_construction_C_longer_sensitivity'
SMOKE=p.ROOT/'results/_smoke/nlsp_construction_C_longer_sensitivity'
AUTH='explicit_user_construction_C_out_period_continuation_2026_10_11'
REPEAT_AUTH='explicit_user_construction_C_output_repeat_1_2026_10_11'
REPEAT_CONFIG=p.ROOT/'data/input/nlsp_construction_C_longer_sensitivity_output_retry.json'
PARENT=p.ROOT/'results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd'
PARENT_SHA='b54b9271260b55841b260692d88f8802f430d780fe21a5d91fe049f91705f3dd'

def verify_parent():
    if p.sha(PARENT/'manifest.json')!=PARENT_SHA:raise ValueError('C parent manifest changed')
    if c.cache(PARENT)['overall']!='C_MODAL_SENSITIVITY_RESOLVED':raise ValueError('Unexpected parent result')
    return p.read(PARENT/'config.json')

def load_level(degree,model):
    pc=p.read(PARENT/'config.json');d=p.make_disc(pc,'C',degree,model)
    base=np.load(PARENT/f'C_p{degree}_base/trajectory.npz')
    tan=np.load(PARENT/f'C_p{degree}_tangent/tangent.npz')
    times=base['times'].copy();Z=tan['vectors'].copy();B=tan['basis'].copy();E=tan['metric'].copy()
    if not np.array_equal(times,tan['times']):raise ValueError('Saved base/tangent timestamps differ')
    if not np.isfinite(base['state']).all() or not np.isfinite(Z).all():raise ValueError('Nonfinite source')
    n=len(d.out)
    if B.shape!=(2*n,3) or np.max(abs(B[n:]))!=0:raise ValueError('Wrong original displacement seed space')
    if np.max(abs(B.T@E@B-np.eye(3)))>2e-10:raise ValueError('Original basis normalization changed')
    if not np.allclose(E[:n,:n],d.K[np.ix_(d.out,d.out)],rtol=2e-13,atol=2e-13):raise ValueError('Reference K mismatch')
    if not np.allclose(E[n:,n:],d.M0[np.ix_(d.out,d.out)],rtol=2e-13,atol=2e-13):raise ValueError('Reference M mismatch')
    states=base['state'].copy();energy=base['energy'].copy()
    dense=a.SavedPlane(PARENT/f'C_p{degree}_base/dense.npz',d)
    last=dense(times[-1]);difference=float(np.max(abs(last-states[-1])))
    if difference>2e-14:raise ValueError('Source dense endpoint mismatch')
    if p.read(PARENT/f'C_p{degree}_base/integration.json')['actual_accepted_end']!=times[-1]:raise ValueError('Base final time mismatch')
    if p.read(PARENT/f'C_p{degree}_tangent/integration.json')['actual_accepted_end']!=times[-1]:raise ValueError('Tangent final time mismatch')
    return d,times,states,Z,B,E,energy,{'restart_source':'saved complete physical coordinate/velocity state and sampled tangent dense endpoint',
        'source_dense_endpoint_max_difference':difference,'adaptive_solver_internal_state_saved':False,
        'restart_is_bitwise_native_path_continuation':False,'original_basis_retained_without_renormalization':True}

@contextmanager
def time_methods(d):
    originals={};durations={};counts={}
    def wrap(owner,key,label):
        fn=getattr(owner,key);originals[(owner,key)]=fn;durations[label]=0.;counts[label]=0
        def timed(*args,**kw):
            start=time.perf_counter()
            try:return fn(*args,**kw)
            finally:durations[label]+=time.perf_counter()-start;counts[label]+=1
        setattr(owner,key,timed)
    for key in ('potential','inertial_terms','rotational_mass','solve_mass','rhs','jacobian','safety','energy'):wrap(d,key,key)
    try:yield durations,counts
    finally:
        for (owner,key),fn in originals.items():setattr(owner,key,fn)

@contextmanager
def time_radau_linear_algebra():
    """Observe the installed Radau LU calls without changing their arguments."""
    from scipy.integrate._ivp import radau
    originals={key:getattr(radau,key) for key in ('lu_factor','lu_solve')};timing={};counts={}
    def wrapper(fn,key):
        timing[key]=0.;counts[key]=0
        def measured(*args,**kwargs):
            start=time.perf_counter()
            try:return fn(*args,**kwargs)
            finally:timing[key]+=time.perf_counter()-start;counts[key]+=1
        return measured
    try:
        for key,fn in originals.items():setattr(radau,key,wrapper(fn,key))
        yield timing,counts
    finally:
        for key,fn in originals.items():setattr(radau,key,fn)

def profile_saved():
    verify_parent();pc=p.read(PARENT/'provenance.json')['physical_source_config'];model,_,_=p.sources(pc)
    SMOKE.mkdir(parents=True,exist_ok=True);report={'new_ODE_calls':0,'profile_kind':'pointwise saved-state evaluations; inclusive nested timings, not production wall fractions','levels':{}}
    for degree in (12,16):
        d,times,states,Z,B,E,energy,source=load_level(degree,model)
        rows={}
        for name in ('rhs','jacobian','safety'):
            profiler=cProfile.Profile();start=time.perf_counter()
            with time_methods(d) as (timings,counts):
                profiler.enable()
                for idx in (0,1000,2000):
                    x=states[idx]
                    for j in range(12):
                        # Rebuild the assembled mass. Local arm caches may remain
                        # warm in this pointwise benchmark; production timings
                        # below measure the actual changing-stage mass costs.
                        d._mass_key=None
                        if name=='safety':d.safety(x[:d.ndof])
                        else:getattr(d,name)(float(times[idx]),x)
                profiler.disable()
            out=io.StringIO();pstats.Stats(profiler,stream=out).strip_dirs().sort_stats('cumulative').print_stats(30)
            (SMOKE/f'profile_p{degree}_{name}.txt').write_text(out.getvalue(),encoding='utf8')
            rows[name]={'evaluations':36,'seconds':time.perf_counter()-start,'inclusive_seconds':dict(timings),'calls':dict(counts)}
        old={name:p.read(PARENT/f'C_p{degree}_{name}/integration.json') for name in ('base','tangent')}
        report['levels'][str(degree)]={'configuration_DOF':d.ndof,'planar_state_DOF':len(d.state_plane),'OUT_state_DOF':len(d.state_out),'matrix_tangent_state_DOF':len(d.state_out)*3,'pointwise':rows,'saved_integrations':old,'restart':source}
    p.write(SMOKE/'performance_profile.json',report);return report

def out_operator(d,t,state):
    """Same full analytic derivative, assembling only required OUT columns.

    Local force_tangents and global acceleration are unchanged. All global
    rows remain until the coupled mass solve; no arm is solved in isolation.
    This avoids computing unused PL input columns, rather than truncating the
    output to the three initial seeds.
    """
    from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import force_tangents
    q,v=state[:d.ndof],state[d.ndof:];acc=d.acceleration(q,v);n=len(d.out)
    first=np.zeros((d.ndof,n));second=np.zeros_like(first)
    for arm,P,qi,vi,ai in zip(d.arms,d.maps,d.local(q),d.local(v),d.local(acc)):
        A,B=force_tangents(arm,qi,vi,ai)
        first+=P.T@A@P[:,d.out];second+=P.T@B@P[:,d.out]
    J=np.zeros((2*n,2*n));J[:n,n:]=np.eye(n)
    J[n:,:n]=d.solve_mass(q,-first)[d.out]
    J[n:,n:]=d.solve_mass(q,-second)[d.out]
    return J

def optimization_check():
    pc=p.read(PARENT/'provenance.json')['physical_source_config'];model,_,_=p.sources(pc);report={}
    for degree in (12,16):
        d,times,states,Z,B,E,energy,source=load_level(degree,model);rows=[]
        for idx in (0,375,745,1000,1500,2000):
            x=states[idx];old=d.jacobian(times[idx],x)[np.ix_(d.state_out,d.state_out)];new=out_operator(d,times[idx],x)
            rows.append({'time':float(times[idx]),'relative_frobenius':float(np.linalg.norm(new-old)/np.linalg.norm(old)),
                'maximum_absolute_difference':float(np.max(abs(new-old)))})
        timing={}
        for name,fn in [('full',lambda:d.jacobian(times[-1],states[-1])[np.ix_(d.state_out,d.state_out)]),
                        ('OUT_columns',lambda:out_operator(d,times[-1],states[-1]))]:
            samples=[]
            for _ in range(40):
                d._mass_key=None;start=time.perf_counter();fn();samples.append(time.perf_counter()-start)
            timing[name]=float(np.median(samples[5:]))
        if max(row['relative_frobenius'] for row in rows)>2e-13:raise ArithmeticError('OUT optimization changes operator')
        zero=np.zeros(2*d.ndof)
        J=out_operator(d,0,zero)
        zero_residual=float(np.linalg.norm(J.T@E+E@J)/np.linalg.norm(E@J))
        if zero_residual>2e-12:raise ArithmeticError('Optimized zero-base energy identity failed')
        report[str(degree)]={'saved_state_comparisons':rows,'median_unprofiled_seconds':timing,
            'speedup':timing['full']/timing['OUT_columns'],'zero_base_energy_skew_relative':zero_residual,
            'new_ODE_calls':0,'same_operator_preserved':True}
    p.write(SMOKE/'optimization_check.json',report);return report

def validate_scope(config):
    repeated=config.get('authorization')==REPEAT_AUTH
    expected=p.read(REPEAT_CONFIG if repeated else CONFIG)
    if config!=expected or config['authorization'] not in (AUTH,REPEAT_AUTH):raise ValueError('Only frozen authorized C continuation')
    if repeated:
        old=p.ROOT/config['archive_failure_bundle']
        if (config['additional_p12_base_repeat_authorized']!=1 or config['previous_numerical_seconds_charged']!=450.
            or config['user_repeat_permission']!='Разрешить один технический повтор p12'
            or p.sha(old/'manifest.json')!=config['archive_failure_manifest_sha256']):raise ValueError('Unverified output-repeat authorization')
    if (config['degree']!=12 or config['conditional_degree']!=16 or config['numerical_budget_seconds']>2700
        or config['maximum_nonlinear_trajectories']>2 or config['maximum_tangent_trajectories']>2
        or config['new_FEM_jobs'] or config['new_full_finite_perturbation'] or config['automatic_retry']
        or config['p16_new_maximum_trigger']!=.005):raise ValueError('Continuation scope expansion')
    return config

def prepare(path=CONFIG):
    config=validate_scope(p.read(path));verify_parent()
    unfinished=None
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        item=p.read(old)
        if item['config']['authorization']==config['authorization']:
            if item['config']!=config:raise ValueError('Authorized continuation already exists with different settings')
            if (old.parent/'manifest.json').exists():return old.parent,item,p.cache(old.parent)
            if (old.parent/'summary.json').exists() or any(old.parent.glob('C_p*_continuation')):raise ValueError('Unmanifested scientific attempt preserved')
            unfinished=old.parent
    if not (SMOKE/'performance_profile.json').exists() or not (SMOKE/'optimization_check.json').exists():raise ValueError('Read-only performance gate missing')
    parent_item=p.read(PARENT/'provenance.json');modal=p.read(PARENT/'frozen_modal_excitation.json')
    spectrum=p.read(PARENT/'restricted_subspace.json')
    # The stored spectrum key is read from the actual source, never rounded.
    omega=float(spectrum['omega'][0])
    target=2*np.pi/omega;t0=2*modal['Tref']
    item={'schema':'construction-C-out-period-continuation-v1','config':config,
        'physical_config':p.read(PARENT/'config.json'),'physical_source_config':parent_item['physical_source_config'],
        'source_artifacts':dict(parent_item['source_artifacts']),
        'source_manifests':{**parent_item['source_manifests'],'C_modal':{'path':config['parent_bundle'],'sha256':PARENT_SHA}},
        'code_hashes':{**p.code_hashes(),Path(__file__).relative_to(p.ROOT).as_posix():p.sha(__file__),CONFIG.relative_to(p.ROOT).as_posix():p.sha(CONFIG)},
        'interpreter':sys.executable,'threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')},
        'Tref':modal['Tref'],'omega_OUT1':omega,'physical_start':t0,'target_end':float(target),
        'initial_snapshot':p.read(SMOKE/'initial_state.json')}
    if config.get('authorization')==REPEAT_AUTH:
        item['source_manifests']['archive_failure']={'path':config['archive_failure_bundle'],'sha256':config['archive_failure_manifest_sha256']}
    if any(item['threads'][k]!='1' for k in item['threads']):raise ValueError('One numerical thread required')
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp if unfinished is None else unfinished
    if unfinished is None and b.exists():raise ValueError('Interrupted unmanifested attempt preserved')
    if unfinished is not None:
        # Preparation failed before any solver attempt. Preserve that evidence;
        # finishing its input archive is not a scientific retry.
        shutil.copyfile(b/'provenance.json',b/'preparation_before_output_copy_fix.json')
        item['preparation_completion']='output archive copy repaired before first scientific call'
    b.mkdir(parents=True,exist_ok=unfinished is not None);p.write(b/'config.json',config);p.write(b/'provenance.json',item)
    archive=b/'execution_code';archive.mkdir(exist_ok=unfinished is not None)
    for file in item['code_hashes']:shutil.copyfile(p.ROOT/file,archive/Path(file).name)
    for source in SMOKE.glob('profile*'):shutil.copyfile(source,b/source.name)
    for name in ('performance_profile.json','optimization_check.json','initial_state.json'):shutil.copyfile(SMOKE/name,b/name)
    for file in r.KERNELS:
        if p.sha(p.ROOT/file)!=parent_item['code_hashes'][file]:raise ValueError('Frozen physical kernel changed')
    profile=p.read(b/'performance_profile.json');opt=p.read(b/'optimization_check.json')
    forecasts={}
    for degree in (12,16):
        old=profile['levels'][str(degree)]['saved_integrations'];ratio=(target-t0)/t0
        forecasts[str(degree)]={'base_seconds':float(old['base']['seconds']*ratio),
            'tangent_seconds':float(old['tangent']['seconds']*ratio/opt[str(degree)]['speedup'])}
    plan={'frozen_before_new_ODE':True,'source_horizon':t0,'target_end':target,'target_Tref':target/modal['Tref'],
        'omega_OUT1_exact_saved':omega,'forecasts':forecasts,'reserve_factor':config['cost_reserve_factor'],
        'four_run_total_with_reserve':sum(sum(v.values()) for v in forecasts.values())*config['cost_reserve_factor'],
        'estimate_is_not_guarantee':True,'initial_subspace_unchanged':True,'restart_has_no_integrator_internal_history':True,
        'p16_decision_rule':{k:config[k] for k in ('p16_new_maximum_trigger','historical_spatial_curve_difference','near_end_window_Tref','near_end_trigger')}}
    p.write(b/'frozen_continuation_decision.json',plan)
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'static_Newton':0,'CCX':0,'Gmsh':0},
        'cases':{},'attempts':[],'numerical_seconds':config.get('previous_numerical_seconds_charged',0.),'stages':{'source':'PASS','performance':'PASS','restart':'PASS','p12':'NOT_RUN','p16':'NOT_RUN'},
        'deferred_theory':config['deferred_theory'],'energy_qualification':'same parent componentwise solver tolerances; no independent full-period temporal certification'}
    if config.get('authorization')==REPEAT_AUTH:s['inherited_attempts']=[{'bundle':config['archive_failure_bundle'],'nonlinear_ODE':1,'charged_seconds':450.,'status':'OUTPUT_ARCHIVE_FAILED','separate_explicit_repeat_permission':True}]
    p.finalize(b,item,s);return b,item,s

def continuation_times(old_times,Tref,target):
    dt=.001*Tref;start=old_times[-1]
    extra=start+np.arange(1,int(np.floor((target-start)/dt))+1)*dt
    extra=extra[extra<target]
    return np.r_[start,extra,target]

def merge_prefix(old_times,old_states,new_times,new_states):
    if new_times[0]!=old_times[-1] or not np.array_equal(new_states[0],old_states[-1]):raise ValueError('Restart discontinuity')
    times=np.r_[old_times,new_times[1:]];states=np.concatenate((old_states,new_states[1:]))
    if np.any(np.diff(times)<=0):raise ValueError('Nonincreasing physical time')
    return times,states

def growth_gate(times,gains,old_max,Tref,config):
    ids=times>=times[-1]-config['near_end_window_Tref']*Tref
    late=gains[ids];tau=times[ids]/Tref
    slope=float(np.polyfit(tau,late,1)[0]);rise=float(late[-1]-late[0]);tol=config['historical_spatial_curve_difference']
    growing=bool(slope>0 and rise>tol and late[-1]>=late.max()-tol)
    increase=float(gains.max()-old_max)
    return {'new_global_maximum_increase':increase,'near_end_slope_per_Tref':slope,'near_end_rise':rise,
        'near_end_window_Tref':config['near_end_window_Tref'],'clearly_rising_near_end':growing,
        'p16_scientific_trigger':bool(increase>config['p16_new_maximum_trigger'] or growing),
        'trigger_is_administrative_not_physical_threshold':True}

def post_level(b,item,d,degree,old_times,old_state,old_Z,B,E,old_energy,new_times,new_state,new_Z):
    started=time.perf_counter();times,state=merge_prefix(old_times,old_state,new_times,new_state)
    _,Z=merge_prefix(old_times,old_Z,new_times,new_Z)
    result,gains=a.gain_result(d,Z,times,B,E);_,oldgain=a.gain_result(d,old_Z,old_times,B,E)
    original=p.read(PARENT/f'C_p{degree}_gain.json')
    if abs(float(oldgain.max())-original['gain'])>3e-14:raise ArithmeticError('Historical gain is not reproduced')
    combo=np.asarray(result['initial_combination'])
    oldcombo=np.asarray(original['initial_combination'])
    if combo@oldcombo<0:combo=-combo
    response=np.einsum('tij,j->ti',Z,combo);initial=B@combo
    delta,scaling=c.perturbation_scale(d,initial,item['physical_config'])
    physical_state=np.zeros_like(state);physical_state[:,d.state_out]=response*delta
    physical,_=a.physical_observations(d,physical_state)
    q0=np.zeros(d.ndof);q0[d.out]=initial[:len(d.out)]
    initial_fields,_=a.physical_observations(d,np.r_[q0,np.zeros(d.ndof)][None,:])
    oldq=np.zeros(d.ndof);oldq[d.out]=(B@oldcombo)[:len(d.out)]
    oldfields,_=a.physical_observations(d,np.r_[oldq,np.zeros(d.ndof)][None,:])
    shape_cosine=float(np.sum(initial_fields[...,2:5]*oldfields[...,2:5])/np.sqrt(np.sum(initial_fields[...,2:5]**2)*np.sum(oldfields[...,2:5]**2)))
    gate=growth_gate(times,gains,original['gain'],item['Tref'],item['config'])
    energy=np.r_[old_energy,np.array([d.energy(x[:d.ndof],x[d.ndof:]) for x in new_state[1:]])]
    end_checks=[]
    for idx in np.linspace(len(old_times)-1,len(times)-1,9,dtype=int):
        x=state[idx];J=d.jacobian(times[idx],x)
        off=max(np.max(abs(J[np.ix_(d.state_plane,d.state_out)])),np.max(abs(J[np.ix_(d.state_out,d.state_plane)])))
        joint=d.joint_diagnostics(x[:d.ndof],x[d.ndof:]);safety=d.safety(x[:d.ndof])
        if off>2e-12 or np.max(abs(x[d.state_out]))!=0:raise ArithmeticError('Planar invariance failed')
        if max(joint[k] for k in ('translation_difference','orientation_difference','contraction_difference'))>2e-12:raise ArithmeticError('Joint identity failed')
        end_checks.append({'time':times[idx],'safety':safety,'joint':joint,'off_block_J_max':off})
    result.update(gate,degree=degree,initial_combination=combo,original_maximum=original['gain'],
        added_interval_maximum=float(gains[len(old_times):].max()),gain_at_end=float(gains[-1]),
        historical_prefix_maximum_difference=float(abs(oldgain.max()-original['gain'])),
        optimal_original_seed_cosine=float(combo@oldcombo),optimal_physical_shape_cosine=shape_cosine,
        physical_scaling=scaling,delta_E0=float(delta),physical_max_v_Phi_psi=np.max(abs(physical[...,2:5]),axis=(0,1,2)),
        Rv=float(np.max(abs(physical[...,2]))/.005),initial_basis_renormalized_at_restart=False,
        mechanical_energy_drift_from_original_t0=float(np.max(abs((energy-energy[0])/energy[0]))),
        restart_energy_relative_jump=float((energy[len(old_times)-1]-old_energy[-1])/old_energy[0]),
        temporal_new_interval='NOT_RUN: no additional IVP authorized; parent tight fixed-optimal control covers old prefix only',
        postprocessing_seconds=time.perf_counter()-started)
    np.savez_compressed(b/f'p{degree}_merged.npz',times=times,state=state,vectors=Z,basis=B,metric=E,gain=gains,energy=energy,
        optimal_normalized_out=response,physical_OUT_profiles=physical[...,2:5],initial_physical_OUT_profiles=initial_fields[...,2:5],old_initial_physical_OUT_profiles=oldfields[...,2:5])
    p.write(b/f'p{degree}_result.json',result);p.write(b/f'p{degree}_safety_checks.json',end_checks)
    with (b/f'p{degree}_gain.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.writer(stream);writer.writerow(['t','t_over_Tref','G_disp3','source','max_v','max_Phi','max_psi','energy_relative_drift'])
        for idx,t in enumerate(times):writer.writerow([t,t/item['Tref'],gains[idx],'immutable_parent_sample' if idx<len(old_times) else 'new_continuation_dense_sample',*np.max(abs(physical[idx,...,2:5]),axis=(0,1)),(energy[idx]-energy[0])/energy[0]])
    return result

def run_level(b,item,s,degree,model):
    if s['hard_stop']:return None
    d,oldtimes,oldstate,oldZ,B,E,oldenergy,restart=load_level(degree,model)
    times=continuation_times(oldtimes,item['Tref'],item['target_end']);config=item['config']
    p.write(b/f'p{degree}_restart.json',restart)
    np.savez_compressed(b/f'p{degree}_restart.npz',physical_state=oldstate[-1],tangent_matrix=oldZ[-1],original_basis=B,metric=E,physical_time=oldtimes[-1],original_initial_state=oldstate[0])
    forecast=p.read(b/'frozen_continuation_decision.json')['forecasts'][str(degree)]
    needed=sum(forecast.values())*config['cost_reserve_factor']
    if needed>config['numerical_budget_seconds']-s['numerical_seconds']:
        s['stop_reason']='Remaining stage forecast exceeds remaining numerical budget';s['overall']='CONSTRUCTION_C_LONGER_NUMERICAL_PARTIAL';return None
    timing=p.read(PARENT/f'C_p{degree}_base/time_settings.json');label=f'C_p{degree}_base_continuation'
    with c.operator_progress(d,b,label,times[-1]),time_methods(d) as (method_times,method_counts),time_radau_linear_algebra() as (LU_times,LU_counts):
        newstate,record=p.integrate_attempt(b,item,s,label,d,oldstate[-1],times,timing,time.perf_counter()+config['numerical_budget_seconds']-s['numerical_seconds'],plane=True,t_start=times[0])
    p.write(b/f'p{degree}_base_production_profile.json',{'inclusive_method_seconds':method_times,'method_calls':method_counts,'Radau_linear_algebra_seconds':LU_times,'Radau_linear_algebra_calls':LU_counts,'inclusive_timings_are_not_additive':True,'total_integration_archive_seconds':record['seconds']})
    if s['hard_stop']:return None
    # A slower-than-planned base prefix updates the remaining-stage forecast.
    multiplier=max(1.,record['seconds']/forecast['base_seconds'])
    needed=forecast['tangent_seconds']*multiplier*config['cost_reserve_factor'];remaining=config['numerical_budget_seconds']-s['numerical_seconds']
    p.write(b/f'p{degree}_tangent_resource_gate.json',{'forecast_with_reserve':needed,'remaining_seconds':remaining,'passed':needed<=remaining,'base_actual_seconds':record['seconds']})
    if needed>remaining:s.update(overall='CONSTRUCTION_C_LONGER_NUMERICAL_PARTIAL',stop_reason='Tangent continuation forecast exceeds remaining budget');return None
    ti={**item,'config':{**config,'rtol':config['tangent_rtol'],'atol_relative':config['tangent_atol_relative']}}
    oldrun=p.read(PARENT/f'C_p{degree}_tangent/integration.json')
    label=f'C_p{degree}_tangent_continuation';last=[0.];start=time.perf_counter();operator_seconds=[0.];operator_calls=[0]
    def observed_operator(t,x):
        tic=time.perf_counter();J=out_operator(d,t,x);now=time.perf_counter();operator_seconds[0]+=now-tic;operator_calls[0]+=1
        if now-last[0]>60:
            p.write(b/(label+'_progress.json'),{'last_operator_evaluation_time':t,'target_time':times[-1],'wall_seconds':now-start,'qualification':'trial operator evaluation, not accepted state'});last[0]=now
        return J
    newZ,tanrecord=a.tangent_attempt(b,ti,s,label,d,b/f'C_p{degree}_base_continuation',times,B,E,
        time.perf_counter()+remaining,max_step=oldrun['max_step'],method=oldrun['method'],sampled_dense_only=True,
        t_start=times[0],initial_matrix=oldZ[-1],operator_override=observed_operator)
    p.write(b/f'p{degree}_tangent_production_profile.json',{'OUT_operator_seconds':operator_seconds[0],'OUT_operator_calls':operator_calls[0],'total_integration_archive_seconds':tanrecord['seconds'],'linear_algebra_method':'DOP853; no implicit Radau LU for tangent'})
    if s['hard_stop']:return None
    result=post_level(b,item,d,degree,oldtimes,oldstate,oldZ,B,E,oldenergy,times,newstate,newZ)
    s['stages'][f'p{degree}']='PASS';p.finalize(b,item,s);return result

def run(b,item,s,through='all'):
    if s['completed'] or s['hard_stop']:return s
    if through=='preflight':return s
    if s['calls']['nonlinear_ODE']==0 and s['calls']['tangent_ODE']==0:
        item['execution_code_hashes']={**p.code_hashes(),Path(__file__).relative_to(p.ROOT).as_posix():p.sha(__file__)}
        archive=b/'production_execution_code';archive.mkdir(exist_ok=True)
        for file in item['execution_code_hashes']:shutil.copyfile(p.ROOT/file,archive/Path(file).name)
        p.write(b/'provenance.json',item);p.finalize(b,item,s)
    model,_,_=p.sources(item['physical_source_config'])
    result=run_level(b,item,s,12,model)
    if result is None:s['completed']=True;p.finalize(b,item,s);return s
    plan=p.read(b/'frozen_continuation_decision.json');remaining=item['config']['numerical_budget_seconds']-s['numerical_seconds']
    resource=sum(plan['forecasts']['16'].values())*item['config']['cost_reserve_factor']
    gate={'p12_trigger':result['p16_scientific_trigger'],'remaining_seconds':remaining,'p16_forecast_with_reserve':resource,
        'p16_authorized_and_affordable':bool(result['p16_scientific_trigger'] and resource<=remaining)}
    p.write(b/'p16_decision.json',gate)
    if gate['p16_authorized_and_affordable']:run_level(b,item,s,16,model)
    if not s['hard_stop']:
        s['overall']='CONSTRUCTION_C_LONGER_SENSITIVITY_COMPLETE_WITH_QUALIFICATIONS' if s['stages']['p16']=='PASS' else 'CONSTRUCTION_C_LONGER_EXPLORATORY_P12'
        if not result['p16_scientific_trigger']:s['stop_reason']='No >0.005 new global maximum and no declared near-end growth trigger; p16 not run'
        elif not gate['p16_authorized_and_affordable']:s['stop_reason']='p16 trigger met but resource forecast exceeds remaining budget'
    s['completed']=True;p.finalize(b,item,s);return s

def cache(b):
    verify_parent();return p.cache(b)

def export_saved(b):
    """Pure saved-array/record analysis: never constructs a spatial model."""
    item=p.read(b/'provenance.json');summary=p.read(b/'summary.json');rows=[]
    for degree in (12,16):
        for kind in ('base','tangent'):
            old=p.read(PARENT/f'C_p{degree}_{kind}/integration.json')
            for source,record in [('parent',old),('continuation',summary['cases'].get(f'C_p{degree}_{kind}_continuation'))]:
                if record is None:continue
                duration=record['target_end']-record.get('physical_start',0.)
                rows.append({'p':degree,'kind':kind,'source':source,'physical_duration':duration,'seconds':record['seconds'],
                    'seconds_per_physical_time':record['seconds']/duration,'seconds_per_Tref':record['seconds']/duration*item['Tref'],
                    'accepted_steps':record['accepted_steps'],'nfev':record['nfev'],'njev':record['njev'],'nlu':record['nlu'],
                    'min_step':record.get('actual_min_step'),'max_step':record.get('actual_max_step')})
    with (b/'computational_cost.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    p.write(b/'computational_cost.json',rows)
    result={'new_interval_temporal_control':'NOT_RUN','parent_tight_optimal_E0_difference':1.9335025137769784e-9,
        'parent_sampled_spatial_gain_difference':item['config']['historical_spatial_curve_difference'],
        'all_field_spatial_certification':'PARTIAL (inherited, particularly c and velocities)',
        'restart_error_independently_integrated':'NOT_RUN','restart_error_evidence':'complete join states identical; source PL dense endpoint equals saved state; same tolerances and max_step; no adaptive internal history saved',
        'old_zero_reference_reused_without_new_ODE':True,'new_full_finite_perturbation':'NOT_RUN (forbidden)',
        'new_global_stability_or_parameter_maps':'NOT_RUN'}
    if (b/'p16_merged.npz').exists():
        with np.load(b/'p12_merged.npz') as z:times,g12,physical12,initial12=(z[k].copy() for k in ('times','gain','physical_OUT_profiles','initial_physical_OUT_profiles'))
        with np.load(b/'p16_merged.npz') as z:
            if not np.array_equal(times,z['times']):raise ValueError('New physical comparison times differ')
            g16=z['gain'];physical16=z['physical_OUT_profiles'];initial16=z['initial_physical_OUT_profiles']
            result.update(spatial_gain_curve_max_difference=float(np.max(abs(g16-g12))),
                new_interval_spatial_gain_curve_max_difference=float(np.max(abs(g16[times>item['physical_start']]-g12[times>item['physical_start']]))),
                max_gain_difference=float(abs(g16.max()-g12.max())),
                physical_profile_difference_separately_optimal=np.max(abs(physical16-physical12),axis=(0,1,2)),
                initial_physical_profile_difference_separately_optimal=np.max(abs(initial16-initial12),axis=(0,1,2)),
                physical_comparison_qualification='same original projected seed-space rule; each degree uses its own optimal original input, so physical curve differences combine discretization and optimizer change')
        np.savez_compressed(b/'spatial_sensitivity.npz',times=times,gain_difference=g16-g12)
    else:result['new_interval_spatial_control']='NOT_RUN: conditional p16 gate not selected or not affordable; p12 exploratory'
    p.write(b/'numerical_qualifications.json',result)
    p.write(b/'export_provenance.json',{'sources':['summary.json','p12_merged.npz','p16_merged.npz if available'],'new_scientific_calls':0,'phase_or_amplitude_fitting':False})
    return result

def plot_saved(b):
    """At most three figures, all from actual saved numerical arrays."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    item=p.read(b/'provenance.json');figures=[];Tref=item['Tref'];available=[]
    for degree in (12,16):
        if (b/f'p{degree}_merged.npz').exists():available.append(degree)
    if not available:return []
    def finish(fig,name):
        fig.tight_layout()
        for ext in ('png','pdf'):fig.savefig(b/f'{name}.{ext}',dpi=180,bbox_inches='tight');figures.append(f'{name}.{ext}')
        plt.close(fig)
    fig,ax=plt.subplots(figsize=(8.8,4.1))
    for degree in available:
        with np.load(b/f'p{degree}_merged.npz') as z:ax.plot(z['times']/Tref,z['gain'],label=f'p{degree}, full OUT response; 3 initial displacement seeds')
    ax.axhline(1,color='grey',ls=':',label='conservative zero-base reference (theoretical)')
    ax.axvline(2,color='black',lw=.9,ls='--',label='restart at 2 Tref');ax.set(xlabel='t/Tref',ylabel='restricted energy gain Gdisp(3)');ax.legend(fontsize=8);ax.grid(alpha=.2)
    finish(fig,'C_longer_restricted_gain')
    degree=available[-1];result=p.read(b/f'p{degree}_result.json')
    with np.load(b/f'p{degree}_merged.npz') as z:times=z['times'];profiles=z['physical_OUT_profiles']
    fig,axes=plt.subplots(3,1,figsize=(9,7),sharex=True)
    for j,field in enumerate(('v','Phi','psi')):
        for arm in range(2):
            # Signed trajectories at predeclared informative arm sections.
            section=100 if j==0 else 50
            axes[j].plot(times/Tref,profiles[:,arm,section,j],label=f'arm {arm+1}, s/L={section/200:g}')
        axes[j].axvline(2,color='grey',ls='--');axes[j].set(ylabel=field+(' [Lref]' if j==0 else ' [rad]'));axes[j].grid(alpha=.2);axes[j].legend(fontsize=8)
    axes[-1].set_xlabel('t/Tref');fig.suptitle(f'p{degree}: fixed original optimal input; initial max(v,Lref Phi,Lref psi)=5e-5')
    finish(fig,'C_longer_physical_OUT_response')
    fig,axes=plt.subplots(1,2,figsize=(10,3.8))
    if len(available)==2:
        with np.load(b/'spatial_sensitivity.npz') as z:axes[0].plot(z['times']/Tref,z['gain_difference'])
        axes[0].set(xlabel='t/Tref',ylabel='Gdisp(3), p16 minus p12');axes[0].axvline(2,color='grey',ls='--')
    else:
        with np.load(b/'p12_merged.npz') as z:axes[0].plot(z['times']/Tref,z['gain']-1)
        axes[0].axhline(item['config']['historical_spatial_curve_difference'],color='grey',ls=':',label='old p12/p16 difference; not a new error bound')
        axes[0].legend(fontsize=7)
        axes[0].set(xlabel='t/Tref',ylabel='p12 gain minus 1',title='New p16 comparison NOT_RUN')
    costs=p.read(b/'computational_cost.json');labels=[];values=[]
    for row in costs:
        if row['source']=='continuation':labels.append(f"p{row['p']} {row['kind']}");values.append(row['seconds_per_Tref'])
    axes[1].bar(labels,values);axes[1].set(ylabel='wall seconds per Tref',title='Actual continuation costs');axes[1].tick_params(axis='x',rotation=25)
    for ax in axes:ax.grid(axis='y',alpha=.2)
    finish(fig,'C_longer_numerical_cost');p.write(b/'figure_provenance.json',{'figures':figures,'new_scientific_calls':0,'zero_control_label':'exact conservative reference plus previously saved numerical check; no new zero-base trajectory','smoothing':False})
    return figures

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--compute',action='store_true');group.add_argument('--report-only',type=Path);group.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--through-stage',choices=('preflight','all'),default='all')
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:b=(args.report_only or args.plot_only).resolve();s=cache(b)
    else:
        b,item,s=prepare(args.config)
        try:s=run(b,item,s,args.through_stage)
        except (ArithmeticError,ValueError,TimeoutError,RuntimeError,np.linalg.LinAlgError) as error:
            s.update(overall='BLOCKED',hard_stop=True,stop_reason=type(error).__name__+': '+str(error));p.finalize(b,item,s)
    if s['completed'] and (b/'p12_result.json').exists() and not (b/'export_provenance.json').exists():
        export_saved(b);plot_saved(b);p.finalize(b,p.read(b/'provenance.json'),s)
    print(json.dumps(p.plain({'bundle':str(b),'overall':s['overall'],'stages':s['stages'],'calls':s['calls'],'numerical_seconds':s['numerical_seconds'],'stop_reason':s.get('stop_reason')}),indent=2))
