"""Scoped diagnostic of one saved torsional tangent check, with no gain scan.

The parent is immutable. All initial vectors, the norm and the base dense
polynomial are read directly from its checked files. No equilibrium/eigen
problem is solved and no mathematical kernel is modified.
"""
from __future__ import annotations
import argparse,csv,json,shutil,sys,time,subprocess
from pathlib import Path
import numpy as np
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_coupled_sensitivity_analysis as a

CONFIG=p.ROOT/'data/input/nlsp_coupled_torsion_tangent_diagnostic.json'
OUTPUT=p.ROOT/'results/nlsp_coupled_torsion_tangent_diagnostic'
VERSION='coupled-torsion-tangent-diagnostic-v1'
AUTH='explicit_user_torsion_tangent_diagnostic_2026_10_10'

def validate_scope(c):
    frozen={'authorization':AUTH,'geometry':'A','degree':16,'horizon_Tref':.1,'saved_torsional_seed_column':2,
        'maximum_new_ODE_integrations':6,'numerical_budget_seconds':1200,'new_static_equilibria':0,
        'new_base_trajectories_without_evidence':0,'new_eigenanalyses':0,'new_FEM_jobs':0,'automatic_retry':False}
    if any(c.get(key)!=value for key,value in frozen.items()):raise ValueError('Diagnostic scope cannot be expanded')
    return c

def finalize(b,item,s):
    p.write(b/'summary.json',s)
    p.write(b/'manifest.json',{'schema':VERSION,'parent_manifest_sha256':item['config']['parent_manifest_sha256'],
        'artifact_hashes':{f.relative_to(b).as_posix():p.sha(f) for f in sorted(b.rglob('*')) if f.is_file() and f!=b/'manifest.json'}})

def verify_parent(c):
    parent=p.ROOT/c['parent_bundle']
    if p.sha(parent/'manifest.json')!=c['parent_manifest_sha256']:raise ValueError('Parent manifest mismatch')
    summary=p.cache(parent)
    if summary['overall']!=c['historical_status_unchanged']:raise ValueError('Historical failure status mismatch')
    return parent

def prepare(config):
    c=validate_scope(p.read(config))
    if c!=p.read(CONFIG):raise ValueError('Only the frozen one-case diagnostic is authorized')
    parent=verify_parent(c)
    item={'schema':VERSION,'config':c,'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],cwd=p.ROOT,text=True).strip(),
        'python':sys.version,'initial_snapshot':'results/_smoke/nlsp_coupled_torsion_tangent_diagnostic/initial_state.json'}
    import hashlib
    fingerprint=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fingerprint
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        if p.read(old)['config']['authorization']==c['authorization'] and old.parent!=b:raise ValueError('Existing diagnostic authorization; no hash-based retry')
    if (b/'manifest.json').exists():return b,p.read(b/'provenance.json'),cache(b),parent
    if b.exists():raise ValueError('Preserve an interrupted unmanifested diagnostic')
    b.mkdir(parents=True);item['initial_code_hashes']=p.code_hashes();item['diagnostic_helper_sha256']=p.sha(Path(__file__))
    folder=b/'initial_code';folder.mkdir()
    shutil.copyfile(__file__,folder/Path(__file__).name)
    p.write(b/'provenance.json',item);p.write(b/'config.json',c)
    s={'overall':'PREPARED','stages':{'read_only':'NOT_RUN','primary_refinement':'NOT_RUN','conditional_control':'NOT_RUN'},
        'hard_stop':False,'completed':False,'attempts':[],'cases':{},'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'CCX':0,'Gmsh':0},'numerical_seconds':0.,
        'parent_status':'BLOCKED_BY_TANGENT_CHECK','maximum_new_integrations':6,
        'execution_mode':'EXPLORATORY_NOT_CERTIFIED','admitted':False}
    finalize(b,item,s);return b,item,s,parent

def cache(b):
    verify_parent(p.read(b/'config.json'))
    for path,h in p.read(b/'manifest.json')['artifact_hashes'].items():
        if p.sha(b/path)!=h:raise ValueError('Diagnostic artifact changed: '+path)
    return p.read(b/'summary.json')

def load_case(parent):
    c=p.read(parent/'config.json');model,_,_=p.sources(c);d=p.make_disc(c,'A',16,model)
    with np.load(parent/'A_short_tangent/tangent.npz') as z:times,history,basis,metric=(z[k].copy() for k in ('times','vectors','basis','metric'))
    with np.load(parent/'A_p16_static.npz') as z:base0=np.r_[z['nonlinear'],np.zeros(d.ndof)]
    delta0=c['perturbation_energy_norm_over_base']*np.sqrt(2*p.read(parent/'A_p16_static.json')['energy_at_release'])
    full=[];initial=[]
    for j in range(2):
        case=parent/f'A_check_torsional_{j}'
        with np.load(case/'trajectory.npz') as z:
            if not np.array_equal(times,z['times']):raise ValueError('Saved comparison timestamps differ')
            full.append(z['state'].copy())
        with np.load(case/'initial.npz') as z:initial.append(z['state'].copy())
    return d,times,history[:,:,2],basis,metric,base0,delta0,np.array(full),np.array(initial)

def norm_curves(error,metric):
    n=error.shape[-1]//2
    coordinate=np.einsum('...i,ij,...j->...',error[...,:n],metric[:n,:n],error[...,:n])
    velocity=np.einsum('...i,ij,...j->...',error[...,n:],metric[n:,n:],error[...,n:])
    return {'total':np.sqrt(np.maximum(coordinate+velocity,0)), 'coordinate':np.sqrt(np.maximum(coordinate,0)),
        'velocity':np.sqrt(np.maximum(velocity,0)), 'velocity_fraction':velocity/np.maximum(coordinate+velocity,1e-300)}

def bilinear_partition(x,metric,groups):
    """Exact signed row allocation; diagonal norms plus explicit cross terms.

    Cross-coupled elastic/internal-joint groups are not independent energies.
    Signed shares sum to one but can be negative. No component is dropped.
    """
    total=float(x@metric@x);allocated={};diagonal={}
    for name,ids in groups.items():
        ids=np.asarray(ids);allocated[name]=float(x[ids]@(metric@x)[ids]/max(total,1e-300))
        diagonal[name]=float(x[ids]@metric[np.ix_(ids,ids)]@x[ids])
    return {'squared_norm':total,'signed_fractions':allocated,'diagonal_projected_squared_norms':diagonal,
        'cross_squared_norm':total-sum(diagonal.values()),'qualification':'signed bilinear allocation, not independent positive group energies'}

def physical_kinetic_parts(d,velocity):
    result=[]
    for arm,local in zip(d.arms,d.local(velocity)):
        fields=arm.reconstruct(local);c=arm.coefficients
        result.append([float(arm.weights@(c.m*fields[:,2]**2)),
            float(arm.weights@((c.jp+c.jb)*fields[:,3]**2)),float(arm.weights@(c.jb*fields[:,4]**2))])
    return np.array(result)

def stage_read_only(b,parent):
    d,t,Z,basis,metric,base0,delta0,full,initial=load_case(parent);n=len(d.out)
    deltas=np.array([delta0,delta0/2]);errors=full[:,:,d.state_out]/deltas[:,None,None]-Z[None]
    scale=max(float(a.energy_norm(Z,metric).max()),1.);curves=norm_curves(errors,metric);values=curves['total'].max(axis=1)/scale
    saved=p.read(parent/'stopped_prefix_diagnostics.json')['shapes']['torsional']['rows']
    if not np.allclose(values,[r['normalized_out_error_relative'] for r in saved],rtol=2e-14,atol=0.):raise ValueError('Historical metric does not reproduce')
    if not np.allclose(metric[:n,:n],d.K[np.ix_(d.out,d.out)],rtol=1e-14,atol=0.) or not np.allclose(metric[n:,n:],d.M0[np.ix_(d.out,d.out)],rtol=1e-14,atol=0.):raise ValueError('Saved reference norm differs')
    np.linalg.cholesky(metric)
    report={'historical_relative_errors':values,'deltas':deltas,'fixed_relative_denominator':scale,
        'saved_timestamps_exact':True,'physical_end':float(t[-1]),'OUT_fields':['v','Phi','psi'],'coordinate_and_velocity_state_dimension':2*n,
        'velocity_source':'direct solver state components, never differentiated reconstructed coordinates',
        'new_ODE_integrations':0,'modal_basis_coverage':'only the saved three modes/six initial states; no new eigenanalysis',
        'rows':[],'initial_checks':[],'dense_checks':[],'pointwise_checks':[]}
    W=d.W[np.ix_(d.out,d.out)];inverse=np.linalg.inv(W);raw_metric=np.zeros_like(metric)
    raw_metric[:n,:n]=inverse.T@metric[:n,:n]@inverse;raw_metric[n:,n:]=inverse.T@metric[n:,n:]@inverse
    joint=np.flatnonzero(np.isin(d.out,d.joint_ids));internal=np.setdiff1d(np.arange(n),joint)
    groups={'internal':np.r_[internal,internal+n],'joint':np.r_[joint,joint+n]}
    displacement=np.r_[np.arange(2*(d.p-1),3*(d.p-1)),7*(d.p-1)+np.arange(2*(d.p-1),3*(d.p-1)),d.joint_ids[2]]
    translation=np.flatnonzero(np.isin(d.out,displacement));rotation=np.setdiff1d(np.arange(n),translation)
    field_groups={'displacement_coordinates':translation,'rotation_coordinates':rotation,'translational_velocities':translation+n,'rotational_velocities':rotation+n}
    gram=basis.T@metric@basis
    for j,delta in enumerate(deltas):
        error=errors[j];index=int(np.argmax(curves['total'][j]));q=np.zeros(d.ndof);v=q.copy();q[d.out]=error[index,:n];v[d.out]=error[index,n:]
        parts=d.reference_energy_parts(q,v);kinetic=physical_kinetic_parts(d,v);total=curves['total'][j,index]**2
        projection=(basis@np.linalg.solve(gram,basis.T@metric@error[index]));low=float(projection@metric@projection/total)
        velocity_rayleigh=float(np.sqrt(v@d.K@v/max(v@d.M0@v,1e-300)))
        raw_error=np.r_[W@error[index,:n],W@error[index,n:]]
        row={'delta':float(delta),'time_at_max':float(t[index]),'relative_error_max':float(values[j]),
            'absolute_error_before_dividing_by_delta':float(delta*curves['total'][j,index]),
            'coordinate_norm_at_max':float(curves['coordinate'][j,index]),'velocity_norm_at_max':float(curves['velocity'][j,index]),
            'coordinate_norm_max':float(curves['coordinate'][j].max()),'velocity_norm_max':float(curves['velocity'][j].max()),
            'velocity_fraction_at_max':float(curves['velocity_fraction'][j,index]),'arm_bending_torsion_fractions':parts/total,
            'kinetic_v_Phi_psi_fractions_by_arm_of_total':kinetic/total,'saved_low_subspace_fraction':low,
            'velocity_error_Rayleigh_omega':velocity_rayleigh,'raw_internal_joint':bilinear_partition(raw_error,raw_metric,groups),
            'raw_field_partition':bilinear_partition(raw_error,raw_metric,field_groups)}
        kinetic_contributions=error[index,n:]*(metric[n:,n:]@error[index,n:]);order=np.argsort(abs(kinetic_contributions))[::-1]
        row['top_white_velocity_DOFs']=[{'coupled_index':int(d.out[k]),'fraction_of_total':float(kinetic_contributions[k]/total)} for k in order[:10]]
        report['rows'].append(row)
        expected=base0.copy();expected[d.state_out]+=delta*basis[:,2]
        report['initial_checks'].append({'delta':float(delta),'exact_global_state_match':bool(np.array_equal(initial[j],expected)),
            'planar_addition_max':float(np.max(abs(initial[j,d.state_plane]-base0[d.state_plane]))),
            'initial_velocity_max':float(np.max(abs(initial[j,d.ndof:]))),
            'normalized_OUT_difference':float(a.energy_norm(initial[j,d.state_out]/delta-basis[:,2],metric)),
            'joint':d.joint_diagnostics(initial[j,:d.ndof],initial[j,d.ndof:])})
        case=parent/f'A_check_torsional_{j}'
        evaluated=p.evaluate_dense_records(case/'dense.npz',t)
        with np.load(case/'dense.npz') as z:steps=z['t']-z['t_old']
        norm_diff=a.energy_norm((evaluated[:,d.state_out]-full[j][:,d.state_out])/delta,metric)
        settings=p.read(case/'time_settings.json');atol=np.array(settings['atol'])[d.state_out]
        report['dense_checks'].append({'full_dense_vs_saved_normalized_max':float(norm_diff.max()),'accepted_steps':len(steps),
            'actual_min_step':float(steps.min()),'actual_max_step':float(steps.max()),'rtol':settings['rtol'],'configured_max_step':settings['max_step'],
            'atol_E0_norm_over_delta':float(a.energy_norm(atol,metric)/delta),'sampled_dense_calls_inferred':len(t)-1,
            'full_native_timestamps_available':True,'old_tangent_native_timestamps_and_dense':'NOT_AVAILABLE'})
    if not all(row['exact_global_state_match'] for row in report['initial_checks']):raise ValueError('Saved initial state mismatch')
    report['old_tangent_settings']=p.read(parent/'A_short_tangent/integration.json')
    report['old_tangent_settings']['actual_min_step']='NOT_AVAILABLE';report['old_tangent_settings']['actual_max_step']='NOT_AVAILABLE'
    for index in (0,25,50,100):
        base=p.evaluate_dense_records(parent/'A_p16_base/dense.npz',t[index:index+1])[0];x=np.zeros(2*d.ndof);x[d.state_plane]=base
        direction=np.zeros(2*d.ndof);direction[d.state_out]=Z[index];exact=d.jacobian(t[index],x)@direction
        rows=[]
        for factor in (16.,8.,4.,2.,1.,.5):
            delta=delta0*factor;central=(d.rhs(t[index],x+delta*direction)-d.rhs(t[index],x-delta*direction))/(2*delta)
            difference=central-exact
            rows.append({'delta_factor':factor,'delta':float(delta),'relative_OUT':float(np.linalg.norm(difference[d.state_out])/max(np.linalg.norm(exact[d.state_out]),1e-300)),
                'absolute_velocity_rows':float(np.linalg.norm(difference[d.out])),
                'relative_acceleration_rows':float(np.linalg.norm(difference[d.out+d.ndof])/max(np.linalg.norm(exact[d.out+d.ndof]),1e-300)),
                'raw_joint_acceleration_max':float(np.max(abs((d.W@difference[d.ndof:])[d.joint_ids[[2,3,4]]])))})
        h=delta0;central1=(d.rhs(t[index],x+h*direction)-d.rhs(t[index],x-h*direction))/(2*h)
        central2=(d.rhs(t[index],x+h/2*direction)-d.rhs(t[index],x-h/2*direction))/h
        limit=(4*central2-central1)/3
        limit_relative=float(np.linalg.norm((limit-exact)[d.state_out])/max(np.linalg.norm(exact[d.state_out]),1e-300))
        if limit_relative>2e-7:raise ValueError('Pointwise tangent operator consistency failed before ODE')
        orders=[float(np.log(first['relative_OUT']/second['relative_OUT'])/np.log(2)) for first,second in zip(rows,rows[1:])]
        report['pointwise_checks'].append({'time':float(t[index]),'rows':rows,'observed_orders':orders,'two_size_limit_relative':limit_relative})
    report['old_full_solver_counters']=[p.read(parent/f'A_check_torsional_{j}/integration.json') for j in range(2)]
    np.savez_compressed(b/'historical_error_decomposition.npz',times=t,errors=errors,deltas=deltas,metric=metric,
        total=curves['total'],coordinate=curves['coordinate'],velocity=curves['velocity'],velocity_fraction=curves['velocity_fraction'],relative_denominator=scale)
    p.write(b/'read_only_diagnostics.json',report)
    return report

def freeze_refinement(b,item,parent):
    """Freeze before ODE results; use measured old cost and unchanged IC/norm."""
    if (b/'frozen_refinement.json').exists():return p.read(b/'frozen_refinement.json')
    s=p.read(b/'summary.json')
    if s['attempts'] or sum(s['calls'].values()):raise ValueError('Cannot choose refinement after integration results')
    old=p.read(parent/'A_check_torsional_0/time_settings.json');tan=p.read(parent/'A_short_tangent/integration.json')
    # Installed Radau uses exponent 1/4; existing frozen forecasts stay intact.
    r=item['config']['primary_refinement'];factor=max(r['rtol_factor']**(-1/4),1/r['max_step_factor'])
    measured=tan['seconds']+sum(p.read(parent/f'A_check_torsional_{j}/integration.json')['seconds'] for j in range(2))
    estimate=1.5*factor*measured
    if estimate>item['config']['numerical_budget_seconds']:raise TimeoutError('Forecast exceeds diagnostic budget before integration')
    result={'rtol':old['rtol']*r['rtol_factor'],'full_atol':np.array(old['atol'])*r['atol_factor'],
        'tangent_atol':np.array(tan['atol_componentwise'])[:tan['full_response_dimension']]*r['atol_factor'],
        'atol_relative':1e-10*r['atol_factor'],'max_step':old['max_step']*r['max_step_factor'],
        'per_nonlinear_run_seconds':item['config']['per_job_seconds'],'per_tangent_run_seconds':item['config']['per_job_seconds'],
        'maximum_nonlinear_trajectories':item['config']['maximum_new_ODE_integrations'],
        'one_torsion_column':'exact saved column 2; independent linear column equation, every OUT DOF retained',
        'base':'same saved A_p16_base native dense polynomials; no base rerun',
        'reason':'100-fold rtol/atol reduction resolves observed divided-by-delta absolute-error scale; max_step halved',
        'estimated_primary_seconds':estimate,'cost_estimate_method':'old measured complete six-column tangent plus two full runs, factor max(tolerance^(-1/4), max-step ratio) times 1.5; planning only',
        'frozen_before_results':True,'no_physics_changes':True}
    p.write(b/'frozen_refinement.json',result);return result

def mismatch_matrix(old_full,new_full,old_tangent,new_tangent,deltas,metric,scale):
    """All four pairings share the historical fixed denominator."""
    result=[];curves={}
    for j,delta in enumerate(deltas):
        row={'delta':float(delta)}
        for fname,full in [('old_full',old_full),('tight_full',new_full)]:
            for tname,tangent in [('old_tangent',old_tangent),('tight_tangent',new_tangent)]:
                key=fname+'/'+tname;curve=a.energy_norm(full[j]/delta-tangent,metric)
                row[key]=float(curve.max()/scale);curves[f'{j}_{key}']=curve
        row['D_full']=float(a.energy_norm((new_full[j]-old_full[j])/delta,metric).max()/scale)
        row['D_tangent']=float(a.energy_norm(new_tangent-old_tangent,metric).max()/scale)
        result.append(row)
    return result,curves

def run_primary(b,item,s,parent):
    if s['hard_stop'] or s['stages']['primary_refinement']!='NOT_RUN':return s
    plan=freeze_refinement(b,item,parent);d,t,Z,basis,metric,base0,delta0,full,initial=load_case(parent)
    deadline=time.perf_counter()+max(0,item['config']['numerical_budget_seconds']-s['numerical_seconds'])
    source=item.copy();source['config']={**item['config'],**plan};source['source_manifests']={'parent':{'path':item['config']['parent_bundle'],'sha256':item['config']['parent_manifest_sha256']}}
    p.write(b/'execution_code_hashes.json',{**p.code_hashes(),'scripts/lib/nlsp_coupled_torsion_diagnostic.py':p.sha(Path(__file__))})
    target=b/'execution_code';target.mkdir(exist_ok=True);shutil.copyfile(__file__,target/Path(__file__).name)
    if sum(s['calls'].values())>=item['config']['maximum_new_ODE_integrations']:raise TimeoutError('ODE ceiling')
    tangent,record=a.tangent_attempt(b,source,s,'tight_tangent',d,parent/'A_p16_base',t,basis[:,[2]],metric,deadline,
        max_step=plan['max_step'],dense_archive=True)
    finalize(b,item,s)
    if s['hard_stop']:return s
    new=[]
    for j in range(2):
        if sum(s['calls'].values())>=item['config']['maximum_new_ODE_integrations']:raise TimeoutError('ODE ceiling')
        timing={'rtol':plan['rtol'],'atol':plan['full_atol'],'max_step':plan['max_step'],'normalization':'saved full componentwise atol times 0.01; unchanged field scales'}
        history,record=p.integrate_attempt(b,source,s,f'tight_full_{j}',d,initial[j].copy(),t,timing,deadline)
        finalize(b,item,s)
        if s['hard_stop']:return s
        new.append(history[:,d.state_out])
    scale=p.read(b/'read_only_diagnostics.json')['fixed_relative_denominator'];deltas=np.array([delta0,delta0/2])
    table,curves=mismatch_matrix(full[:,:,d.state_out],np.array(new),Z,tangent[:,:,0],deltas,metric,scale)
    np.savez_compressed(b/'primary_temporal_comparison.npz',times=t,deltas=deltas,metric=metric,old_tangent=Z,tight_tangent=tangent[:,:,0],
        old_full=full[:,:,d.state_out],tight_full=np.array(new),**{key.replace('/','_'):value for key,value in curves.items()})
    p.write(b/'primary_temporal_comparison.json',{'rows':table,'new_ODE_calls':3,'same_saved_base':True,'same_initial_states':True,'historical_denominator_unchanged':True})
    s['stages']['primary_refinement']='COMPLETE';s['overall']='PRIMARY_REFINEMENT_COMPLETE';finalize(b,item,s);return s

def choose_conditional_level(table):
    """Variant B: exactly the trajectory with greatest observed sensitivity."""
    index=int(np.argmax([row['D_full'] for row in table]))
    return ('tangent',None) if max(row['D_tangent'] for row in table)>=table[index]['D_full'] else ('full',index)

def run_conditional(b,item,s,parent):
    if s['hard_stop'] or s['stages']['conditional_control']!='NOT_RUN':return s
    if s['stages']['primary_refinement']!='COMPLETE':raise ValueError('Primary refinement must finish first')
    table=p.read(b/'primary_temporal_comparison.json')['rows'];which,index=choose_conditional_level(table)
    label='extra_tangent' if which=='tangent' else f'extra_full_{index}'
    prior='tight_tangent' if which=='tangent' else f'tight_full_{index}'
    plan=p.read(b/'frozen_refinement.json');level=item['config']['conditional_extra_level']
    forecast=s['cases'][prior]['seconds']*max(level['rtol_factor']**(-1/4),1/level['max_step_factor'])*1.5
    if s['numerical_seconds']+forecast>item['config']['numerical_budget_seconds']:
        s['stages']['conditional_control']='NOT_RUN_RESOURCE_LIMIT';finalize(b,item,s);return s
    if sum(s['calls'].values())>=item['config']['maximum_new_ODE_integrations']:raise TimeoutError('Total ODE ceiling')
    decision={'variant':'B','trajectory_kind':which,'delta_index':index,'case_label':label,
        'reason':'Dominant old temporal sensitivity is established; independently check the remaining small primary mismatch using only the most sensitive trajectory',
        'rtol':plan['rtol']*level['rtol_factor'],'max_step':plan['max_step']*level['max_step_factor'],
        'full_atol':np.array(plan['full_atol'])*level['atol_factor'],'atol_relative':plan['atol_relative']*level['atol_factor'],
        'estimated_seconds':forecast,'elapsed_plus_forecast_seconds':s['numerical_seconds']+forecast,
        'frozen_before_extra_results':True,'negative_delta_trajectories':0}
    if (b/'conditional_decision.json').exists():
        if p.read(b/'conditional_decision.json')!=p.plain(decision):raise ValueError('Frozen conditional decision changed')
    else:p.write(b/'conditional_decision.json',decision)
    if label in s['cases']:raise ValueError('No automatic repeat of conditional trajectory')
    phase=b/'conditional_execution_code';phase.mkdir();shutil.copyfile(__file__,phase/Path(__file__).name)
    p.write(phase/'source_hashes.json',{**p.code_hashes(),'scripts/lib/nlsp_coupled_torsion_diagnostic.py':p.sha(Path(__file__))})
    d,t,Z,basis,metric,base0,delta0,oldfull,initial=load_case(parent)
    source=item.copy();source['config']={**item['config'],**plan,**decision};source['source_manifests']={'parent':{'path':item['config']['parent_bundle'],'sha256':item['config']['parent_manifest_sha256']}}
    deadline=time.perf_counter()+max(0,item['config']['numerical_budget_seconds']-s['numerical_seconds'])
    if which=='full':
        full,record=p.integrate_attempt(b,source,s,label,d,initial[index].copy(),t,
            {'rtol':decision['rtol'],'atol':decision['full_atol'],'max_step':decision['max_step'],'normalization':'same physical scales; primary component atol times 0.1'},deadline)
        result=full[:,d.state_out]
        with np.load(b/prior/'trajectory.npz') as z:previous=z['state'][:,d.state_out]
        delta=delta0/2**index
        difference=a.energy_norm((result-previous)/delta,metric)
        with np.load(b/'tight_tangent/tangent.npz') as z:tan=z['vectors'][:,:,0]
        mismatch=a.energy_norm(result/delta-tan,metric)
    else:
        full,record=a.tangent_attempt(b,source,s,label,d,parent/'A_p16_base',t,basis[:,[2]],metric,deadline,
            max_step=decision['max_step'],dense_archive=True)
        result=full[:,:,0]
        with np.load(b/prior/'tangent.npz') as z:previous=z['vectors'][:,:,0]
        difference=a.energy_norm(result-previous,metric);mismatch=None
    finalize(b,item,s)
    if s['hard_stop']:return s
    scale=p.read(b/'read_only_diagnostics.json')['fixed_relative_denominator']
    np.savez_compressed(b/'conditional_comparison.npz',times=t,extra_OUT=result,change=difference,
        mismatch=[] if mismatch is None else mismatch,relative_denominator=scale)
    p.write(b/'conditional_comparison.json',{'variant':'B','kind':which,'delta_index':index,
        'observed_change_relative':float(difference.max()/scale),'extra_full_tight_tangent_relative':None if mismatch is None else float(mismatch.max()/scale),
        'new_ODE_calls':1,'does_not_establish_a_uniform_extra_level_pair':True})
    s['stages']['conditional_control']='COMPLETE';s['overall']='TEMPORAL_DIAGNOSTIC_DATA_COMPLETE';s['completed']=True;finalize(b,item,s);return s

def temporal_order(values):
    return float(np.log(values[0]/values[1])/np.log(2))

def postprocess_saved(b):
    """Read-only numerical diagnostics; no new dynamics or modal solve."""
    if (b/'final_metrics.json').exists():return p.read(b/'final_metrics.json')
    parent=verify_parent(p.read(b/'config.json'));d,t,Z,basis,metric,base0,delta0,old,initial=load_case(parent)
    r=p.read(b/'read_only_diagnostics.json');scale=r['fixed_relative_denominator'];deltas=np.array(r['deltas'])
    with np.load(b/'primary_temporal_comparison.npz') as z:new=z['tight_full'];newtan=z['tight_tangent']
    errors=new/deltas[:,None,None]-newtan[None];curves=norm_curves(errors,metric)
    native=[];tolerance=[];execution=[]
    for j in range(2):
        oldcase=parent/f'A_check_torsional_{j}';newcase=b/f'tight_full_{j}'
        with np.load(oldcase/'dense.npz') as z:
            index=min(int(np.searchsorted(z['t'],r['rows'][j]['time_at_max'],side='left')),len(z['t'])-1)
            start=z['y_old'][index];end=z['y_old'][index+1] if index+1<len(z['t']) else start+z['Q'][index].sum(axis=1)
        settings=p.read(oldcase/'time_settings.json');atol=np.array(settings['atol'])[d.out+d.ndof]
        relative_scale=settings['rtol']*np.maximum(abs(start[d.out+d.ndof]),abs(end[d.out+d.ndof]))
        normalized_error=old[j][:,d.state_out]/deltas[j]-Z
        error=normalized_error[int(np.argmax(norm_curves(normalized_error,metric)['total']))]
        ev=error[len(d.out):];contributions=ev*(metric[len(d.out):,len(d.out):]@ev)
        dominance=float(np.sum(contributions*atol/(atol+relative_scale))/contributions.sum())
        tolerance.append({'delta':float(deltas[j]),'absolute_part_of_velocity_error_weighted_solver_scale':dominance,
            'scope':'Radau atol+rtol*max(abs(start),abs(end)) on the saved accepted interval containing the old sampled maximum; not a measured error bound'})
        for name,case in [('old_full',oldcase),('tight_full',newcase)]:
            with np.load(case/'dense.npz') as z:native_times=np.r_[0.,z['t']]
            with np.load(b/'tight_tangent/dense.npz') as z:tan_times=z['t']
            grid=np.unique(np.r_[native_times,tan_times])
            full=p.evaluate_dense_records(case/'dense.npz',grid)[:,d.state_out]
            tan=p.evaluate_dense_records(b/'tight_tangent/dense.npz',grid)
            mismatch=a.energy_norm(full/deltas[j]-tan,metric)/scale
            native.append({'delta':float(deltas[j]),'full_case':name,'tangent_case':'tight_tangent',
                'native_union_count':len(grid),'mismatch_max_on_union':float(mismatch.max()),'time_at_max_on_union':float(grid[np.argmax(mismatch)]),
                'origin':'union of two actual accepted-end timestamp lists; counterpart values from native dense polynomials, not native output at every union time',
                'old_tangent_union_comparison':'NOT_AVAILABLE: old tangent accepted times/dense not saved'})
        evaluated=p.evaluate_dense_records(newcase/'dense.npz',t)[:,d.state_out]
        native.append({'delta':float(deltas[j]),'tight_full_dense_vs_saved_grid_normalized_max':float(a.energy_norm((evaluated-new[j])/deltas[j],metric).max())})
    result={'historical_relative_errors':r['historical_relative_errors'],
        'primary_relative_errors':curves['total'].max(axis=1)/scale,
        'primary_observed_order':temporal_order(curves['total'].max(axis=1)),
        'fixed_relative_denominator':scale,'native_timestamp_checks':native,'absolute_tolerance_diagnostics':tolerance,
        'initial_states_unchanged':True,'new_base_trajectories':0,'new_static_solves':0,'new_eigenanalyses':0,
        'roundoff_or_base_error_not_independently_bounded':True,
        'rtol_and_atol_not_isolated_by_separate_jobs':True,
        'historical_floor_and_thresholds_unchanged':True,'parent_status_unchanged':'BLOCKED_BY_TANGENT_CHECK'}
    for label in ('tight_tangent','tight_full_0','tight_full_1','extra_full_1'):
        case=b/label
        if not case.exists():continue
        stats=p.read(case/'integration.json')
        with np.load(case/'dense.npz') as z:steps=z['t']-z['t_old']
        timing=stats if label=='tight_tangent' else p.read(case/'time_settings.json')
        execution.append({'case':label,**{key:stats.get(key) for key in ('status','actual_accepted_end','accepted_steps','nfev','njev','nlu','seconds','energy_relative_drift_max')},
            'rtol':timing['rtol'],'configured_max_step':timing['max_step'],'actual_min_step':float(steps.min()),'actual_max_step':float(steps.max())})
    result['execution']=execution
    # Conditional cancellation is a diagnostic only; it does not replace data.
    result['conditional_zero_amplitude_residual']=a.two_amplitude_diagnostic(new[0]/deltas[0],new[1]/deltas[1],newtan,metric)
    np.savez_compressed(b/'tight_error_decomposition.npz',times=t,errors=errors,deltas=deltas,metric=metric,
        total=curves['total'],coordinate=curves['coordinate'],velocity=curves['velocity'],velocity_fraction=curves['velocity_fraction'],relative_denominator=scale)
    with (b/'error_time_histories.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.writer(stream);writer.writerow(['level','delta','time','relative_total','coordinate_E0_norm','velocity_E0_norm','velocity_fraction','absolute_E0_error_before_delta_division'])
        with np.load(b/'historical_error_decomposition.npz') as z:oldcurves={key:z[key] for key in ('total','coordinate','velocity','velocity_fraction')}
        for level,data in [('old',oldcurves),('tight',curves)]:
            for j,delta in enumerate(deltas):
                for k in range(len(t)):writer.writerow([level,delta,t[k],data['total'][j,k]/scale,data['coordinate'][j,k],data['velocity'][j,k],data['velocity_fraction'][j,k],delta*data['total'][j,k]])
    for filename,rows in [('temporal_pairings.csv',p.read(b/'primary_temporal_comparison.json')['rows']),('actual_execution.csv',execution)]:
        with (b/filename).open('w',newline='',encoding='utf8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    p.write(b/'final_metrics.json',result);return result

def plot_saved(b):
    """Only saved arrays. Two figures, no smoothing or scientific calls."""
    if (b/'figure_provenance.json').exists():return p.read(b/'figure_provenance.json')['figures']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with np.load(b/'historical_error_decomposition.npz') as z:old={key:z[key] for key in ('times','deltas','total','coordinate','velocity')}
    with np.load(b/'tight_error_decomposition.npz') as z:tight={key:z[key] for key in ('total','coordinate','velocity')};scale=float(z['relative_denominator'])
    Tref=p.read(verify_parent(p.read(b/'config.json'))/'frozen_excitation.json')['A']['Tref']
    names=[]
    def finish(fig,name):
        fig.tight_layout()
        for suffix in ('pdf','png'):fig.savefig(b/(name+'.'+suffix),dpi=180,bbox_inches='tight');names.append(name+'.'+suffix)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
    for j,ax in enumerate(axes):
        for group,color in [('total','black'),('coordinate','tab:blue'),('velocity','tab:orange')]:
            ax.semilogy(old['times'][1:]/Tref,old[group][j,1:]/scale,color=color,ls='--',label='old '+group)
            ax.semilogy(old['times'][1:]/Tref,tight[group][j,1:]/scale,color=color,label='tight '+group)
        ax.set(title=f'delta = {old["deltas"][j]:.6g}',xlabel='t/Tref',ylabel='Error / fixed historical denominator');ax.legend(fontsize=7)
    finish(fig,'torsion_coordinate_velocity_error')
    fig,ax=plt.subplots(figsize=(6.5,4.3));deltas=old['deltas'];before=old['total'].max(axis=1)/scale;after=tight['total'].max(axis=1)/scale
    ax.loglog(deltas,before,'o--',label='old');ax.loglog(deltas,after,'o-',label='tight, same policy at both delta')
    ax.loglog(deltas,after[0]*(deltas/deltas[0])**2,':',color='gray',label='delta squared guide, anchored to first tight point')
    if (b/'conditional_comparison.json').exists():
        extra=p.read(b/'conditional_comparison.json')
        if extra['kind']=='full':ax.loglog([deltas[extra['delta_index']]],[extra['extra_full_tight_tangent_relative']],'*',ms=11,label='extra full level, one delta only')
    ax.axhline(1e-8,color='gray',ls='--',lw=.7,label='historical floor (unchanged)')
    ax.set(xlabel='Initial perturbation E0 norm delta',ylabel='Sampled max error / fixed denominator');ax.legend(fontsize=7)
    finish(fig,'torsion_error_amplitude_refinement')
    p.write(b/'figure_provenance.json',{'figures':names,'new_ODE_calls':0,'no_smoothing_or_time_alignment':True});return names

def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--through-stage',choices=('read-only','primary','all'),default='read-only')
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=cache(b)
        if args.plot_only and not (b/'figure_provenance.json').exists():
            postprocess_saved(b);s['figures']=plot_saved(b);finalize(b,p.read(b/'provenance.json'),s)
    else:
        b,item,s,parent=prepare(args.config)
        if s['stages']['read_only']=='NOT_RUN':
            try:stage_read_only(b,parent);s['stages']['read_only']='PASS'
            except (ValueError,ArithmeticError) as error:s.update(overall='BLOCKED',hard_stop=True,reason=str(error))
            finalize(b,item,s)
        if args.through_stage!='read-only' and not s['hard_stop']:
            try:s=run_primary(b,item,s,parent)
            except TimeoutError as error:s.update(overall='BLOCKED',hard_stop=True,reason=str(error));finalize(b,item,s)
        if args.through_stage=='all' and not s['hard_stop']:s=run_conditional(b,item,s,parent)
    print(json.dumps(p.plain({'bundle':str(b),'overall':s['overall'],'stages':s['stages'],'calls':s['calls']}),ensure_ascii=False,indent=2))
