"""Authorized, bounded continuation; reuses the frozen coupled kernel and IVPs.

No mathematical model, static initial state A, or historical result is changed.
The four-dimensional input space retains the complete output state.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_coupled_sensitivity_analysis as a
from scripts.lib import nlsp_coupled_torsion_diagnostic as g

CONFIG=p.ROOT/'data/input/nlsp_coupled_planar_sensitivity_resume.json'
OUTPUT=p.ROOT/'results/nlsp_coupled_planar_sensitivity_resume'
AUTH='explicit_user_coupled_sensitivity_resume_2026_10_10'
VERSION='coupled-planar-sensitivity-resume-v1'
KERNELS=['scripts/lib/weakly_nonlinear_spatial_rod.py','scripts/lib/weakly_nonlinear_spatial_dynamics.py',
    'scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_coupled_spatial_dynamics.py',
    'scripts/lib/reddy_inplane_geometry.py','scripts/lib/yartsev_ch2_monoclinic_rod.py']

def validate_scope(c):
    if c!=p.read(CONFIG) or c['authorization']!=AUTH:raise ValueError('Only the explicit frozen continuation is authorized')
    fixed={'horizon_Tref':2.,'mixed_horizon_Tref':.1,'saved_seed_columns':[0,2,3,5],
        'full_rtol':1e-12,'full_atol_factor':.01,'full_max_step_factor':.5,'finite_delta_factor':1.,
        'tangent_full_relative_gate':1e-5,'numerical_floor_relative':1e-8}
    if any(c.get(key)!=value for key,value in fixed.items()):raise ValueError('Frozen protocol cannot be changed')
    if c['new_FEM_jobs'] or c['automatic_retry'] or c['degree']!=16 or c['lower_degree']!=12:raise ValueError('Scope expansion')
    if c['maximum_nonlinear_trajectories']>8 or c['maximum_tangent_trajectories']>16 or c['numerical_budget_seconds']>7200:raise ValueError('Integration ceiling')
    return c

def verify_sources(c):
    paths={}
    for key,expected in [('parent','BLOCKED_BY_TANGENT_CHECK'),('diagnostic','NUMERICAL_FLOOR_IDENTIFIED')]:
        b=p.ROOT/c[key+'_bundle']
        if p.sha(b/'manifest.json')!=c[key+'_manifest_sha256']:raise ValueError(key+' manifest mismatch')
        s=p.cache(b) if key=='parent' else g.cache(b)
        if s['overall']!=expected:raise ValueError('Unexpected historical status')
        paths[key]=b
    source=p.read(paths['parent']/'A_p16_base/executed_source_hashes.json')
    for file in KERNELS:
        if p.sha(p.ROOT/file)!=source[file]:raise ValueError('Frozen kernel changed: '+file)
    return paths

def cache(b):
    b=Path(b);verify_sources(p.read(b/'config.json'));return p.cache(b)

def prepare(config):
    c=validate_scope(p.read(config));parents=verify_sources(c)
    # Authorization blocks cross-hash retries, while completed replay is read-only.
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        olditem=p.read(old)
        if olditem['config']['authorization']==AUTH:
            if olditem['config']!=c:raise ValueError('An authorization already has an attempt')
            return old.parent,olditem,cache(old.parent),parents
    pc=p.read(parents['parent']/'config.json');_,files,manifests=p.sources(pc)
    for name,b in parents.items():
        manifests[name]={'path':b.relative_to(p.ROOT).as_posix(),'sha256':p.sha(b/'manifest.json')}
    hashes=p.code_hashes();hashes[Path(__file__).relative_to(p.ROOT).as_posix()]=p.sha(__file__)
    item={'schema':VERSION,'config':c,'physical_config':pc,'source_artifacts':files,'source_manifests':manifests,
        'initial_code_hashes':hashes,'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=p.ROOT,text=True).strip(),
        'python':sys.version,'interpreter':sys.executable,'threads':1,
        'initial_snapshot':'results/_smoke/nlsp_coupled_planar_sensitivity_resume/initial_state.json'}
    fingerprint=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fingerprint
    if b.exists():raise ValueError('Preserve interrupted unmanifested attempt')
    b.mkdir(parents=True);p.write(b/'config.json',c);p.write(b/'provenance.json',item)
    folder=b/'initial_code';folder.mkdir()
    for file in hashes:shutil.copyfile(p.ROOT/file,folder/Path(file).name)
    shutil.copyfile(p.ROOT/item['initial_snapshot'],b/'initial_state.json')
    p.write(b/'frozen_excitation.json',p.read(parents['parent']/'frozen_excitation.json'))
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'stages':{key:'NOT_RUN' for key in ('source','mixed','B_base','gain','spatial','temporal','finite')},
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'static_Newton':0,'CCX':0,'Gmsh':0},'attempts':[],'cases':{},'numerical_seconds':0.,
        'execution_mode':'EXPLORATORY_NOT_CERTIFIED','admitted':False,'deferred_theory':'DEFERRED_TWO_PARAMETER_ACTION_EXPANSION',
        'historical_statuses':['BLOCKED_BY_TANGENT_CHECK','NUMERICAL_FLOOR_IDENTIFIED']}
    s['stages']['source']='PASS';p.finalize(b,item,s);return b,item,s,parents

def fixed_metric(d):
    n=len(d.out);metric=np.zeros((2*n,2*n));metric[:n,:n]=d.K[np.ix_(d.out,d.out)];metric[n:,n:]=d.M0[np.ix_(d.out,d.out)]
    np.linalg.cholesky(metric);return metric

def projected_coordinates(old,new,coordinates):
    """Nested Shen truncation with exactly the same seven joint coordinates."""
    raw=old.raw(coordinates);target=np.zeros(new.ndof)
    for arm in range(2):
        for field in range(7):
            i=arm*7*(old.p-1)+field*(old.p-1);j=arm*7*(new.p-1)+field*(new.p-1)
            target[j:j+new.p-1]=raw[i:i+new.p-1]
    target[new.joint_ids]=raw[old.joint_ids];return new.from_raw(target)

def project_basis(old,new,basis,metric):
    result=np.zeros((2*len(new.out),basis.shape[1]))
    for j in range(basis.shape[1]):
        for part in range(2):
            q=np.zeros(old.ndof);q[old.out]=basis[part*len(old.out):(part+1)*len(old.out),j]
            result[part*len(new.out):(part+1)*len(new.out),j]=projected_coordinates(old,new,q)[new.out]
    gram=result.T@metric@result;v,U=eigh(gram)
    if v.min()<=0:raise ArithmeticError('Projected input subspace is rank deficient')
    orth=result@(U@np.diag(v**-.5)@U.T)
    return orth,{'physical_projection_gram':gram,'orthonormality_error':float(np.max(abs(orth.T@metric@orth-np.eye(basis.shape[1])))),
        'no_new_modes_or_static_states':True,'projection':'nested Shen truncation; seven common endpoint coordinates preserved'}

def mixed_direction():return np.array([1.,1.,0.,0.])/np.sqrt(2.)

def raw_error_details(d,error,metric):
    n=len(d.out);W=d.W[np.ix_(d.out,d.out)];inverse=np.linalg.inv(W);raw_metric=np.zeros_like(metric)
    raw_metric[:n,:n]=inverse.T@metric[:n,:n]@inverse;raw_metric[n:,n:]=inverse.T@metric[n:,n:]@inverse
    joint=np.flatnonzero(np.isin(d.out,d.joint_ids));internal=np.setdiff1d(np.arange(n),joint)
    raw=np.r_[W@error[:n],W@error[n:]]
    return g.bilinear_partition(raw,raw_metric,{'internal':np.r_[internal,internal+n],'joint':np.r_[joint,joint+n]})

def parity_metrics(plus,minus,d,metric,delta):
    out=(plus[:,d.state_out]+minus[:,d.state_out])/delta
    n=len(d.plane);M=np.zeros((2*n,2*n));M[:n,:n]=d.K[np.ix_(d.plane,d.plane)];M[n:,n:]=d.M0[np.ix_(d.plane,d.plane)]
    return {'out_normalized_energy_error':float(a.energy_norm(out,metric).max()),
        'plane_absolute_energy_error':float(a.energy_norm(plus[:,d.state_plane]-minus[:,d.state_plane],M).max()),
        'raw_state_max':float(np.max(abs(plus[:,d.state_plane]-minus[:,d.state_plane])))}

def limits(s,c):
    if s['hard_stop']:raise RuntimeError('Stopped continuation cannot launch another integration')
    if s['calls']['nonlinear_ODE']>c['maximum_nonlinear_trajectories'] or s['calls']['tangent_ODE']>c['maximum_tangent_trajectories']:raise RuntimeError('Call ceiling exceeded')
    if s['numerical_seconds']>=c['numerical_budget_seconds']:raise TimeoutError('Total integration budget exhausted')

def load_states(case):
    with np.load(case/'trajectory.npz') as z:return z['state'].copy()

def tangent(b,item,s,label,d,case,times,basis,metric,deadline,zero=False,tight=False):
    if label in s['cases']:
        with np.load(b/label/'tangent.npz') as z:return z['vectors'].copy()
    limits(s,item['config'])
    if s['calls']['tangent_ODE']>=item['config']['maximum_tangent_trajectories']:raise TimeoutError('Tangent call ceiling')
    c=item['config'];local=dict(item);local['config']=dict(c,rtol=c['tangent_rtol']*(.1 if tight else 1),atol_relative=c['tangent_atol_relative']*(.1 if tight else 1))
    step=p.settings(d,.0025,{'rtol':1e-10,'atol_relative':1e-10})['max_step']*c['tangent_max_step_factor']*(.5 if tight else 1)
    Z,record=a.tangent_attempt(b,local,s,label,d,case,times,basis,metric,deadline,zero,max_step=step,method=c['tangent_method'],sampled_dense_only=True)
    return Z

def full(b,item,s,label,d,initial,times,deadline,plane=False):
    if label in s['cases']:return load_states(b/label)
    limits(s,item['config']);c=item['config'];timing=p.settings(d,.0025,{'rtol':c['base_rtol'],'atol_relative':c['base_atol_relative']})
    if not plane:timing.update(rtol=c['full_rtol'],atol=timing['atol']*c['full_atol_factor'],max_step=timing['max_step']*c['full_max_step_factor'])
    history,record=p.integrate_attempt(b,item,s,label,d,initial,times,timing,deadline,plane)
    return history

def halt(b,item,s,reason,overall='BLOCKED'):
    s.update(hard_stop=True,completed=True,overall=overall,stop_reason=reason);p.finalize(b,item,s);return s

def stage_mixed(b,item,s,parent,d,basis,metric,q0,deadline):
    c=item['config'];exc=p.read(b/'frozen_excitation.json')['A'];times=p.make_times(exc['Tref'],.1,c)
    with np.load(parent/'A_short_tangent/tangent.npz') as z:old=z['vectors'][:,:,c['saved_seed_columns']]
    Z=tangent(b,item,s,'A_short_tangent',d,parent/'A_p16_base',times,basis,metric,deadline)
    if s['hard_stop']:return
    sensitivity=float(max(a.energy_norm(Z[:,:,j]-old[:,:,j],metric).max() for j in range(4)))
    p.write(b/'short_tangent_preflight.json',{'method':c['tangent_method'],'relative_tangent_reference_difference':sensitivity,
        'gate':c['short_tangent_reference_relative_gate'],'all_response_components':Z.shape[1],'same_derivative_of_full_RHS':True})
    if sensitivity>c['short_tangent_reference_relative_gate']:
        halt(b,item,s,'Short independent tangent integrator check failed');return
    delta=.03*np.sqrt(2*p.read(parent/'A_p16_static.json')['energy_at_release']);direction=mixed_direction();seed=basis@direction
    base=a.SavedPlane(parent/'A_p16_base/dense.npz',d);X0=np.array([base(t) for t in times]);rows=[];histories=[];errors=[]
    for j,factor in enumerate(c['delta_factors']):
        initial=q0.copy();initial[d.state_out]+=delta*factor*seed
        p.write(b/f'A_mixed_initial_{j}.json',{'delta':delta*factor,'energy_norm_seed':float(a.energy_norm(seed,metric)),
            'planar_addition_max':float(np.max(abs(initial[d.state_plane]-q0[d.state_plane]))),
            'initial_velocity_max':float(np.max(abs(initial[d.ndof:]))),'joint':d.joint_diagnostics(initial[:d.ndof],initial[d.ndof:])})
        X=full(b,item,s,f'A_mixed_{j}',d,initial,times,deadline)
        if s['hard_stop']:return
        report=a.full_tangent_error(d,X,X0,Z,direction,delta*factor,metric)
        row={key:value for key,value in report.items() if key not in ('out_error_curve','plane_change_curve','predicted','observed')}
        row['fixed_relative_denominator']=max(float(a.energy_norm(report['predicted'],metric).max()),1.)
        row['absolute_error_before_delta_division']=abs(delta*factor)*report['normalized_out_error_max']
        err=report['observed']-report['predicted'];row.update(a.response_error_details(d,err,delta*factor,metric));row['error_groups']=g.bilinear_partition(err[np.argmax(report['out_error_curve'])],metric,
            {'coordinates':np.arange(len(d.out)),'velocities':np.arange(len(d.out),2*len(d.out))})
        rows.append(row);histories.append(X);errors.append(err)
        row['raw_internal_joint_error']=raw_error_details(d,err[np.argmax(report['out_error_curve'])],metric)
    signs=parity_metrics(histories[0],histories[2],d,metric,delta);order=a.measured_orders(rows[:2])[0]
    good=rows[1]['normalized_out_error_relative']<=c['tangent_full_relative_gate'] and rows[1]['normalized_out_error_relative']<=1.05*rows[0]['normalized_out_error_relative']
    resolved=rows[0]['normalized_out_error_relative']>10*c['numerical_floor_relative']
    if resolved:good=good and order['normalized_out_order']>=1.5
    good=good and signs['out_normalized_energy_error']<=c['sign_symmetry_relative_gate'] and signs['plane_absolute_energy_error']<=delta*c['sign_symmetry_relative_gate']
    result={'status':'PASS' if good else 'FAIL','rows':rows,'order':order,'order_resolved_above_historical_floor':resolved,
        'sign_symmetry':signs,'delta0':delta,'direction':direction,'tangent_temporal_preflight':sensitivity,'unchanged_relative_gate':c['tangent_full_relative_gate']}
    p.write(b/'mixed_check.json',result);np.savez_compressed(b/'mixed_check.npz',times=times,tau=times/exc['Tref'],errors=errors,metric=metric,delta=np.array(c['delta_factors'])*delta)
    s['stages']['mixed']=result['status'];p.finalize(b,item,s)
    if not good:halt(b,item,s,'Mixed/sign tangent-full gate failed')

def freeze_long_plan(b,item,s):
    path=b/'frozen_long_plan.json'
    if path.exists():return p.read(path)
    c=item['config'];exc=p.read(b/'frozen_excitation.json');ratio=exc['B']['Tref']/exc['A']['Tref']
    measured=s['cases']['A_short_tangent']['seconds'];parent=p.ROOT/c['parent_bundle'];diagnostic=p.ROOT/c['diagnostic_bundle']
    # The short run predates this output-only optimization. Installed DOP853
    # requires three extra RHS calls for dense output, without changing its
    # accepted state/steps. Long cases build it only at sampled increments.
    measured*=.90  # conservative planning reduction; not a runtime guarantee
    oldbase=p.read(parent/'A_p16_base/integration.json')['seconds'];oldfull=p.read(diagnostic/'tight_full_0/integration.json')['seconds']
    estimates={'B_base':oldbase*ratio,'A_tangent':measured*20,'B_tangent':measured*20*ratio,
        'A_finite':oldfull*20,'B_finite':oldfull*20*ratio,'zero_controls':15.,
        'B_p12_base':oldbase*ratio*.65,'B_p12_tangent':measured*20*ratio*.65}
    required=sum(v for key,v in estimates.items() if not key.startswith('B_p12'))
    remaining=c['numerical_budget_seconds']-s['numerical_seconds'];spatial=required+estimates['B_p12_base']+estimates['B_p12_tangent']<=remaining
    plan={'frozen_before_long_results':not any(key=='B_p16_base' for key in s['cases']),'method':c['tangent_method'],'input_dimension':4,
        'estimates_seconds':estimates,'completed_numerical_seconds':s['numerical_seconds'],'remaining_seconds':remaining,
        'primary_forecast_seconds':required,'total_forecast_seconds':required+(estimates['B_p12_base']+estimates['B_p12_tangent'] if spatial else 0),
        'spatial_B_full_horizon_planned':spatial,'A_p12_dynamic_control':'NOT_RUN: omitted before long production to respect budget',
        'planning_only_not_guarantee':True,'physical_horizon_not_shortened':True,'same_four_input_directions_each_geometry':True,
        'scope_reduction':'six to four initial states; one B p12 spatial control if forecast fits; no A p12 integrations',
        'temporal_long_control':'NOT_RUN: only short independent tangent control fits the frozen programme',
        'output_optimization':{'dense_output_only_at_requested_samples':True,'physics_or_solver_step_changes':False,
            'installed_DOP853_extra_dense_rhs_per_increment':3,'forecast_factor':.90,'actual_short_rhs':s['cases']['A_short_tangent']['nfev'],
            'actual_short_accepted_steps':s['cases']['A_short_tangent']['accepted_steps']}}
    folder=b/'long_execution_code';folder.mkdir(exist_ok=True)
    for file in (Path(__file__),p.ROOT/'scripts/lib/nlsp_coupled_sensitivity_analysis.py'):shutil.copyfile(file,folder/file.name)
    p.write(folder/'hashes.json',{file.relative_to(p.ROOT).as_posix():p.sha(file) for file in (Path(__file__),p.ROOT/'scripts/lib/nlsp_coupled_sensitivity_analysis.py')})
    p.write(path,plan);return plan

def base_checks(case,d,Tref):
    a.save_base_observations(case,d,Tref);X=load_states(case);checks=[]
    for j in np.linspace(0,len(X)-1,11,dtype=int):checks.append({'sample_index':int(j),'safety':d.safety(X[j,:d.ndof]),'joint':d.joint_diagnostics(X[j,:d.ndof],X[j,d.ndof:])})
    p.write(case/'sampled_safety_joint.json',checks)

def stage_gain(b,item,s,d,name,case,basis,metric,deadline):
    c=item['config'];T=p.read(b/'frozen_excitation.json')[name]['Tref'];times=p.make_times(T,2,c)
    Z=tangent(b,item,s,name+'_p16_tangent',d,case,times,basis,metric,deadline)
    if s['hard_stop']:return
    zero=tangent(b,item,s,name+'_zero_tangent',d,case,times,basis,metric,deadline,zero=True)
    if s['hard_stop']:return
    gram=np.einsum('tij,ik,tkl->tjl',zero,metric,zero);error=float(np.max(abs(gram-np.eye(basis.shape[1]))))
    if error>c['zero_base_gram_gate']:halt(b,item,s,'Zero-base energy gain control failed');return
    result,gain=a.gain_result(d,Z,times,basis,metric);result.update(time_at_max_Tref=result['time_at_max']/T,zero_gram_error=error,
        restriction='G_perp^(4), two spatial directions and their velocity states; all output OUT coordinates retained')
    p.write(b/(name+'_gain.json'),result);np.savez_compressed(b/(name+'_gain.npz'),times=times,tau=times/T,gain=gain,zero_gain=np.sqrt(np.linalg.eigvalsh(gram)[:,-1]))
    p.finalize(b,item,s)

def stage_finite(b,item,s,d,name,case,basis,metric,initial,deadline):
    c=item['config'];T=p.read(b/'frozen_excitation.json')[name]['Tref'];times=p.make_times(T,2,c);result=p.read(b/(name+'_gain.json'))
    direction=np.array(result['initial_combination']);E=d.energy(initial[:d.ndof],initial[d.ndof:]);delta=.03*np.sqrt(2*E)*c['finite_delta_factor']
    perturbed=initial.copy();perturbed[d.state_out]+=delta*(basis@direction)
    X=full(b,item,s,name+'_optimal_finite',d,perturbed,times,deadline)
    if s['hard_stop']:return
    with np.load(b/(name+'_p16_tangent')/'tangent.npz') as z:Z=z['vectors']
    saved=a.SavedPlane(case/'dense.npz',d);base=np.array([saved(t) for t in times]);error=a.full_tangent_error(d,X,base,Z,direction,delta,metric)
    report={key:value for key,value in error.items() if key not in ('out_error_curve','plane_change_curve','predicted','observed')}
    report.update(a.response_error_details(d,error['observed']-error['predicted'],delta,metric));report['status']='PASS' if report['normalized_out_error_relative']<=c['tangent_full_relative_gate'] else 'PARTIAL'
    report['finite_delta_not_fitted']=True;report['out_response_max_energy_norm']=float(a.energy_norm(X[:,d.state_out],metric).max())
    p.write(b/(name+'_finite_check.json'),report);np.savez_compressed(b/(name+'_finite_check.npz'),times=times,tau=times/T,error=error['out_error_curve'],
        predicted=error['predicted'],observed=error['observed'],tangent_norm=a.energy_norm(error['predicted'],metric),full_norm=a.energy_norm(error['observed'],metric),planar_change=error['plane_change_curve'])
    states=np.zeros((2,2*d.ndof));states[0,d.state_out]=basis@direction;states[1,d.state_out]=np.array(result['response_at_max']);profiles,_=a.physical_observations(d,states)
    physical,obs=a.physical_observations(d,X);report['maximum_physical_OUT_v_Phi_psi']=np.max(abs(physical[:,:,:,[2,3,4]]),axis=(0,1,2));p.write(b/(name+'_finite_check.json'),report)
    np.savez_compressed(b/(name+'_optimal_profiles.npz'),x=np.linspace(0,1,201),initial=profiles[0],response=profiles[1],delta=delta,
        finite_observations=obs,time_at_max=result['time_at_max'])
    p.finalize(b,item,s)

def stage_spatial(b,item,s,model,d16,basis16,metric16,initial16,deadline):
    c=item['config'];d12=p.make_disc(item['physical_config'],'B',12,model);metric=fixed_metric(d12)
    basis,projection=project_basis(d16,d12,basis16,metric);q0=projected_coordinates(d16,d12,initial16[:d16.ndof]);initial=np.r_[q0,np.zeros(d12.ndof)]
    # Same physical initial state, projected; no new p12 Newton equilibrium.
    x=np.linspace(0,1,501);differences=[]
    for old,new,qi,qj in zip(d16.arms,d12.arms,d16.local(initial16[:d16.ndof]),d12.local(q0)):
        differences.append(np.max(abs(old.reconstruct(qi,x*old.length)-new.reconstruct(qj,x*new.length)),axis=0))
    projection['initial_physical_field_max_change']=differences;p.write(b/'B_p12_projection.json',projection)
    T=p.read(b/'frozen_excitation.json')['B']['Tref'];times=p.make_times(T,2,c)
    X=full(b,item,s,'B_p12_base',d12,initial,times,deadline,plane=True)
    if s['hard_stop']:return
    base_checks(b/'B_p12_base',d12,T)
    Z=tangent(b,item,s,'B_p12_tangent',d12,b/'B_p12_base',times,basis,metric,deadline)
    if s['hard_stop']:return
    result,curve=a.gain_result(d12,Z,times,basis,metric);high=p.read(b/'B_gain.json');change=abs(result['gain']-high['gain'])/high['gain']
    result.update(relative_gain_change=change,time_at_max_Tref=result['time_at_max']/T,reference_p16_time_at_max_Tref=high['time_at_max_Tref'],
        status='PASS' if change<=c['gain_spatial_relative_gate'] else 'PARTIAL',qualification='B projected identical physical initial subspace/state only; A has no p12 dynamic control')
    p.write(b/'spatial_control.json',result);np.savez_compressed(b/'B_p12_gain.npz',times=times,tau=times/T,gain=curve)
    full16=load_states(b/'B_p16_base');f16,_=a.physical_observations(d16,full16);f12,_=a.physical_observations(d12,X)
    p.write(b/'B_base_spatial_difference.json',{'field_max':np.max(abs(f12-f16),axis=(0,1,2)),
        'common_physical_coordinate':True,'initial_projection_included':True,'no_coefficient_vector_subtraction':True})
    s['stages']['spatial']=result['status'];p.finalize(b,item,s)

def run(b,item,s,parents,through='all'):
    if s['completed'] or s['hard_stop']:return s
    c=item['config'];model,_,_=p.sources(item['physical_config']);dA=p.make_disc(item['physical_config'],'A',16,model)
    parent=parents['parent'];deadline=time.perf_counter()+c['numerical_budget_seconds']-s['numerical_seconds']
    with np.load(parent/'A_short_tangent/tangent.npz') as z:basisA=z['basis'][:,c['saved_seed_columns']];metricA=z['metric']
    if np.max(abs(fixed_metric(dA)-metricA))!=0:raise ValueError('Saved reference norm changed')
    with np.load(parent/'A_p16_static.npz') as z:initialA=np.r_[z['nonlinear'],np.zeros(dA.ndof)]
    if s['stages']['mixed']=='NOT_RUN':stage_mixed(b,item,s,parent,dA,basisA,metricA,initialA,deadline)
    if s['hard_stop'] or through=='mixed':return s
    plan=freeze_long_plan(b,item,s)
    if plan['primary_forecast_seconds']>plan['remaining_seconds']:return halt(b,item,s,'Read-only long-program forecast exceeds remaining integration budget','NUMERICAL_PARTIAL')
    if through=='planning':return s
    dB=p.make_disc(item['physical_config'],'B',16,model);exc=p.read(b/'frozen_excitation.json')['B'];T=exc['Tref'];times=p.make_times(T,2,c)
    if not (b/'B_p16_static.npz').exists():
        start=time.perf_counter();static=p.static_state(dB,exc);s['calls']['static_Newton']+=1
        np.savez_compressed(b/'B_p16_static.npz',**{k:static[k] for k in ('linear','nonlinear','force')});p.write(b/'B_p16_static.json',{k:v for k,v in static.items() if k not in ('linear','nonlinear','force')})
        p.write(b/'B_static_execution.json',{'seconds':time.perf_counter()-start,'not_in_ODE_budget':True,'load_frozen_in_parent':True});p.finalize(b,item,s)
    with np.load(b/'B_p16_static.npz') as z:initialB=np.r_[z['nonlinear'],np.zeros(dB.ndof)]
    full(b,item,s,'B_p16_base',dB,initialB,times,deadline,plane=True)
    if s['hard_stop']:return s
    base_checks(b/'B_p16_base',dB,T);s['stages']['B_base']='PASS';p.finalize(b,item,s)
    if (b/'B_subspace.npz').exists():
        with np.load(b/'B_subspace.npz') as z:basisB=z['basis'];metricB=z['metric']
    else:
        allbasis,metricB,subspace=a.reference_subspace(dB);basisB=allbasis[:,[0,2,3,5]]
        subspace.update(initial_dimension=4,selection='first bending and first torsional direction; displacement and velocity states; before gain results')
        np.savez_compressed(b/'B_subspace.npz',basis=basisB,metric=metricB);p.write(b/'B_subspace.json',subspace)
    np.savez_compressed(b/'A_subspace.npz',basis=basisA,metric=metricA)
    for name,d,case,basis,metric in [('A',dA,parent/'A_p16_base',basisA,metricA),('B',dB,b/'B_p16_base',basisB,metricB)]:
        stage_gain(b,item,s,d,name,case,basis,metric,deadline)
        if s['hard_stop']:return s
    s['stages']['gain']='PASS';p.finalize(b,item,s)
    # Full finite-amplitude confirmation has priority over optional controls.
    for name,d,case,basis,metric,initial in [('A',dA,parent/'A_p16_base',basisA,metricA,initialA),('B',dB,b/'B_p16_base',basisB,metricB,initialB)]:
        stage_finite(b,item,s,d,name,case,basis,metric,initial,deadline)
        if s['hard_stop']:return s
    s['stages']['finite']='PASS' if all(p.read(b/(name+'_finite_check.json'))['status']=='PASS' for name in ('A','B')) else 'PARTIAL'
    remaining=c['numerical_budget_seconds']-s['numerical_seconds'];forecast=sum(plan['estimates_seconds'][key] for key in ('B_p12_base','B_p12_tangent'))
    p.write(b/'spatial_launch_decision.json',{'originally_planned':plan['spatial_B_full_horizon_planned'],'remaining_seconds':remaining,'forecast_seconds':forecast,'no_unbounded_launch':True})
    if plan['spatial_B_full_horizon_planned'] and remaining>1.1*forecast:stage_spatial(b,item,s,model,dB,basisB,metricB,initialB,deadline)
    if s['hard_stop']:return s
    s['stages']['temporal']='PARTIAL';s['overall']='NUMERICAL_PARTIAL';s['completed']=True
    A=p.read(b/'A_gain.json');B=p.read(b/'B_gain.json');s['observed_gain_contrast']=B['gain']-A['gain']
    s['geometry_contrast_qualification']='restricted finite-time contrast observed; A spatial and long-time tangent controls not certified'
    s['scientific_claim']='method checked in four-dimensional initial space, all output DOFs; neither Floquet stability nor physical 3D joint validation'
    p.finalize(b,item,s);return s

def plot_saved(b):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b=Path(b);item=p.read(b/'provenance.json');figures=[]
    def finish(fig,name):
        fig.tight_layout()
        for ext in ('pdf','png'):fig.savefig(b/(name+'.'+ext),dpi=180,bbox_inches='tight');figures.append(name+'.'+ext)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for name,ax in zip(('A','B'),axes):
        case=(p.ROOT/item['config']['parent_bundle']/'A_p16_base') if name=='A' else b/'B_p16_base'
        if not (case/'physical_observations.npz').exists():continue
        with np.load(case/'physical_observations.npz') as z:
            for arm in range(2):ax.plot(z['tau'],z['observations'][:,arm,1,1],label=f'arm {arm+1}')
        ax.set(title=name,xlabel='t / own Tref',ylabel='local w(Li/2)');ax.legend()
    finish(fig,'planar_A_B')
    if (b/'mixed_check.npz').exists():
        fig,axes=plt.subplots(1,2,figsize=(10,3.7))
        with np.load(b/'mixed_check.npz') as z:
            n=z['errors'].shape[-1]//2
            for j in range(2):
                axes[0].semilogy(z['tau'][1:],a.energy_norm(z['errors'][j],z['metric'])[1:],label=f'delta={z["delta"][j]:.3g}')
                axes[1].semilogy(z['tau'][1:],a.energy_norm(z['errors'][j,:,:n],z['metric'][:n,:n])[1:],label=f'q, level {j}')
                axes[1].semilogy(z['tau'][1:],a.energy_norm(z['errors'][j,:,n:],z['metric'][n:,n:])[1:],ls='--',label=f'velocity, level {j}')
        for ax in axes:ax.set(xlabel='t/Tref',ylabel='fixed E0 error norm');ax.legend(fontsize=8)
        finish(fig,'mixed_tangent_full')
    if (b/'A_gain.npz').exists():
        fig,ax=plt.subplots(figsize=(8,4))
        for name in ('A','B'):
            if not (b/(name+'_gain.npz')).exists():continue
            with np.load(b/(name+'_gain.npz')) as z:ax.plot(z['tau'],z['gain'],label=name+' p16');ax.plot(z['tau'],z['zero_gain'],ls=':',label=name+' zero base')
        if (b/'B_p12_gain.npz').exists():
            with np.load(b/'B_p12_gain.npz') as z:ax.plot(z['tau'],z['gain'],ls='--',label='B p12 projected physical subspace')
        ax.set(xlabel='t / own Tref',ylabel='restricted G_perp^(4)(t)');ax.legend();finish(fig,'restricted_gains')
    if all((b/(name+'_finite_check.npz')).exists() for name in ('A','B')):
        fig,axes=plt.subplots(4,2,figsize=(10,11))
        for j,name in enumerate(('A','B')):
            with np.load(b/(name+'_finite_check.npz')) as z:
                axes[0,j].plot(z['tau'],z['tangent_norm'],label='tangent');axes[0,j].plot(z['tau'],z['full_norm'],ls='--',label='full / delta');axes[0,j].set(title=name,xlabel='t/Tref',ylabel='E0 response / initial norm');axes[0,j].legend()
            with np.load(b/(name+'_optimal_profiles.npz')) as z:
                for row,field,label in [(1,2,'v'),(2,3,'Phi'),(3,4,'psi')]:
                    for arm,color in enumerate(('C0','C1')):
                        axes[row,j].plot(z['x'],z['initial'][arm,:,field]*z['delta'],color=color,label=f'initial, arm {arm+1}')
                        axes[row,j].plot(z['x'],z['response'][arm,:,field]*z['delta'],color=color,ls='--',label=f'at gain max, arm {arm+1}')
                    axes[row,j].set(xlabel='s/Li',ylabel='local '+label+(' [rad]' if field>2 else ' [length]'));axes[row,j].legend(fontsize=7)
        finish(fig,'optimal_full_response')
    p.write(b/'figure_provenance.json',{'figures':figures,'saved_data_only':True,'smoothing_or_phase_fitting':False,
        'plot_helper_sha256':p.sha(__file__),'data_hashes':{path.relative_to(b).as_posix():p.sha(path) for path in b.glob('*.npz')}});return figures

def physical_field_max(values,fields=(2,3,4)):
    selected=np.take(values,fields,axis=-1)
    return np.max(abs(selected),axis=tuple(range(selected.ndim-1)))

def saved_safety_diagnostics(b):
    b=Path(b);target=b/'saved_state_safety_joint.json'
    if target.exists():return p.read(target)
    item=p.read(b/'provenance.json');pc=item['physical_config'];model,_,_=p.sources(pc);result={}
    for name in ('A','B'):
        d=p.make_disc(pc,name,16,model)
        for label in ((name+'_optimal_finite',) if name=='A' else ('B_p16_base','B_optimal_finite')):
            if not (b/label/'trajectory.npz').exists():continue
            with np.load(b/label/'trajectory.npz') as z:times=z['times'];states=z['state']
            rows=[d.safety(x[:d.ndof]) for x in states];extrema=[]
            for arm in range(2):
                values=[row[arm] for row in rows];summary={}
                for key in values[0]:
                    if key=='method':continue
                    if key=='mass_positive':summary[key]=all(v[key] for v in values)
                    elif key in ('min_one_plus_c','relative_mass_lower_bound'):summary[key]=min(v[key] for v in values)
                    else:summary[key]=max(v[key] for v in values)
                extrema.append(summary)
            indices=[0,len(times)//2,len(times)-1];joint=[dict(time=float(times[j]),**d.joint_diagnostics(states[j,:d.ndof],states[j,d.ndof:])) for j in indices]
            result[label]={'samples':len(times),'safety_on_saved_physical_grid_not_every_native_step':True,'per_arm':extrema,'joint_at_three_saved_states':joint,'source_sha256':p.sha(b/label/'trajectory.npz')}
    p.write(target,result);return result

def export_saved(b):
    """First export reconstructs saved physical fields; replay reads the cache.

    There are no integrations, equilibria or eigenanalyses here. The unchanged
    spatial maps are assembled only for physical reconstruction of saved
    coefficient vectors. All scientific trajectories remain byte-identical.
    """
    b=Path(b)
    if (b/'export_provenance.json').exists():return p.read(b/'export_provenance.json')
    s=cache(b)
    if not s['completed'] and not s['hard_stop']:raise ValueError('Export only a completed/stopped prefix')
    def table(name,rows):
        if not rows:return
        with (b/name).open('w',newline='',encoding='utf8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    rows=[]
    for label,record in s['cases'].items():
        rows.append({'case':label,'kind':next((x['kind'] for x in s['attempts'] if x['label']==label),'NOT_RUN'),
            **{key:record.get(key) for key in ('status','seconds','accepted_steps','nfev','njev','nlu','actual_accepted_end','target_end','energy_relative_drift_max')}})
    table('actual_execution.csv',rows)
    rows=[]
    for name in ('A','B'):
        if not (b/(name+'_gain.npz')).exists():continue
        with np.load(b/(name+'_gain.npz')) as z:
            for t,tau,gain,zero in zip(z['times'],z['tau'],z['gain'],z['zero_gain']):rows.append({'geometry':name,'physical_time':t,'t_over_own_Tref':tau,'restricted_G4':gain,'zero_base_G4':zero})
    table('gain_histories.csv',rows)
    with np.load(b/'mixed_check.npz') as z:
        rows=[];n=z['errors'].shape[-1]//2
        for j,delta in enumerate(z['delta']):
            norms=g.norm_curves(z['errors'][j],z['metric'])
            for index,(t,tau) in enumerate(zip(z['times'],z['tau'])):
                rows.append({'delta':delta,'physical_time':t,'t_over_Tref':tau,'coordinate_K_error':norms['coordinate'][index],
                    'velocity_M_error':norms['velocity'][index],'total_E0_error':norms['total'][index],'absolute_before_delta_division':abs(delta)*norms['total'][index]})
    table('mixed_error_histories.csv',rows)
    pc=p.read(b/'provenance.json')['physical_config'];model,_,_=p.sources(pc);physical={};composition={};combos=[];finite_rows=[]
    for name in ('A','B'):
        if not (b/(name+'_finite_check.npz')).exists():continue
        d=p.make_disc(pc,name,16,model);result=p.read(b/(name+'_gain.json'));control=p.read(b/(name+'_finite_check.json'));n=len(d.out)
        combination=np.array(result['initial_combination']);delta=control['delta'];index=result['sample_index_at_max'];profile_rows=[]
        for label,value in zip(('bend_coordinate','torsional_coordinate','bend_velocity','torsional_velocity'),combination):combos.append({'geometry':name,'component':label,'coefficient':value,'squared_coefficient':value*value})
        with np.load(b/(name+'_p16_tangent')/'tangent.npz') as z:initial=z['basis']@combination;response=z['vectors'][index]@combination
        with np.load(b/(name+'_zero_tangent')/'tangent.npz') as z:zero=z['vectors'][index]@combination
        q_profiles=[];v_profiles=[];parts=[]
        for vector in (initial,response,zero):
            q=np.zeros(d.ndof);v=q.copy();q[d.out]=vector[:n];v[d.out]=vector[n:]
            q_profiles.append(np.stack([arm.reconstruct(qi,np.linspace(0,arm.length,201)) for arm,qi in zip(d.arms,d.local(q))]))
            v_profiles.append(np.stack([arm.reconstruct(vi,np.linspace(0,arm.length,201)) for arm,vi in zip(d.arms,d.local(v))]))
            energy=d.reference_energy_parts(q,v);parts.append({'parts_by_arm_bending_torsion':energy,'fractions':energy/energy.sum()})
        q_profiles=np.array(q_profiles)*delta;v_profiles=np.array(v_profiles)*delta
        np.savez_compressed(b/(name+'_optimal_physical_state.npz'),x=np.linspace(0,1,201),positions=q_profiles,velocities=v_profiles,
            labels=np.array(['initial','nonlinear_base_tangent_at_gain_max','zero_base_tangent_at_same_time']),delta=delta,time_at_max=result['time_at_max'])
        physical[name]={'max_initial_OUT_v_Phi_psi':physical_field_max(q_profiles[0]),
            'max_initial_OUT_velocity_v_Phi_psi':physical_field_max(v_profiles[0]),
            'max_peak_OUT_v_Phi_psi':physical_field_max(q_profiles[1]),
            'full_out_v_over_declared_linear_w_scale':control['maximum_physical_OUT_v_Phi_psi'][0]/.0025,
            'common_displacement_scale':.0025,'no_individual_curve_normalization':True}
        composition[name]={'initial':parts[0],'nonlinear_base_at_gain_max':parts[1],'zero_base_same_input_time':parts[2]}
        with np.load(b/(name+'_finite_check.npz')) as z:
            for t,tau,pred,fullnorm,error,pl in zip(z['times'],z['tau'],z['tangent_norm'],z['full_norm'],z['error'],z['planar_change']):
                finite_rows.append({'geometry':name,'physical_time':t,'t_over_own_Tref':tau,'delta':delta,'tangent_norm':pred,'full_over_delta_norm':fullnorm,'E0_mismatch':error,'planar_change_E0':pl})
    table('optimal_initial_combinations.csv',combos);table('finite_control_histories.csv',finite_rows)
    p.write(b/'physical_response_metrics.json',physical);p.write(b/'optimal_energy_composition.json',composition)
    folder=b/'postprocessing_code';folder.mkdir(exist_ok=True)
    shutil.copyfile(__file__,folder/Path(__file__).name)
    prov={'scientific_calls':0,'source_helper_sha256':p.sha(__file__),'raw_trajectories_not_modified':True,
        'csv_files':['actual_execution.csv','gain_histories.csv','mixed_error_histories.csv','optimal_initial_combinations.csv','finite_control_histories.csv'],
        'physical_reconstruction_only':True,'all_output_DOFs_retained':True,'no_phase_amplitude_fitting':True}
    p.write(b/'export_provenance.json',prov);p.finalize(b,p.read(b/'provenance.json'),s);return prov

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--through-stage',choices=('mixed','planning','all'),default='all');args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=cache(b)
        if s['completed'] or s['hard_stop']:export_saved(b)
        if args.plot_only and not (b/'figure_provenance.json').exists():s['figures']=plot_saved(b);p.finalize(b,p.read(b/'provenance.json'),s)
    else:
        b,item,s,parents=prepare(args.config);s=run(b,item,s,parents,args.through_stage)
        if s['completed'] or s['hard_stop']:export_saved(b)
        if s['completed'] and not (b/'figure_provenance.json').exists():s['figures']=plot_saved(b);p.finalize(b,item,s)
    print(json.dumps(p.plain({'bundle':str(b),'overall':s['overall'],'stages':s['stages'],'calls':s['calls'],'numerical_seconds':s['numerical_seconds'],'stop_reason':s.get('stop_reason')}),ensure_ascii=False,indent=2))
