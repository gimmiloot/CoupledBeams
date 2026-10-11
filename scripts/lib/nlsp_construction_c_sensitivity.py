"""One modal-amplitude construction C preset of the existing coupled workflow.

Only orchestration, modal normalization and saved-data diagnostics are new.
The action, assembly, full/tangent operators and IVP implementations are reused.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, os, shutil, subprocess, sys, time
from contextlib import contextmanager
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_coupled_sensitivity_analysis as a
from scripts.lib import nlsp_coupled_sensitivity_resume as r

CONFIG=p.ROOT/'data/input/nlsp_construction_C_modal_sensitivity.json'
OUTPUT=p.ROOT/'results/nlsp_construction_C_modal_sensitivity'
AUTH='explicit_user_construction_C_modal_sensitivity_2026_10_11'
VERSION='construction-C-modal-displacement-sensitivity-v1'


def validate_scope(c):
    if c!=p.read(CONFIG) or c['authorization']!=AUTH:
        raise ValueError('Only the explicitly authorized frozen C preset')
    if (c['sections']!={'C':[[.05,.1],[.05,.1]]} or c['lengths']!=[.8,1.2]
        or c['beta_deg']!=30 or c['amplitude']!=.005 or c['Lref']!=1
        or c['material']!={'E':1.,'rho':1.,'nu':.3,'kappa':5/6}
        or c['physical_perturbation_fraction']!=.01 or c['degree']!=16
        or c['lower_degree']!=12 or c['horizon_Tref']!=2
        or c['maximum_nonlinear_trajectories']>3 or c['maximum_planar_trajectories']>2
        or c['maximum_tangent_trajectories']>8 or c['numerical_budget_seconds']>7200
        or c['initial_subspace_dimension']>4 or c['new_FEM_jobs'] or c['automatic_retry']):
        raise ValueError('Unauthorized expansion of C experiment')
    return c


def verify_sources(c):
    parent=p.ROOT/c['parent_bundle']
    if p.sha(parent/'manifest.json')!=c['parent_manifest_sha256']:
        raise ValueError('A/B parent manifest mismatch')
    previous=r.cache(parent)
    if previous['overall']!='NUMERICAL_PARTIAL' or not previous['completed']:
        raise ValueError('Unexpected A/B source status')
    return parent


def cache(b):
    b=Path(b);verify_sources(p.read(b/'config.json'));return p.cache(b)


def prepare(config):
    c=validate_scope(p.read(config));parent=verify_sources(c)
    # Same authorization never creates a new attempt after a code/HEAD change.
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        item=p.read(old)
        if item['config']['authorization']==AUTH:
            if item['config']!=c:raise ValueError('An authorized attempt already exists')
            return old.parent,item,cache(old.parent)
    pc=p.read(parent/'provenance.json')['physical_config']
    _,files,manifests=p.sources(pc)
    manifests['AB_resume']={'path':c['parent_bundle'],'sha256':c['parent_manifest_sha256']}
    hashes=p.code_hashes()
    for file in (Path(__file__),Path(r.__file__),p.ROOT/'scripts/lib/nlsp_coupled_torsion_diagnostic.py'):
        hashes[file.relative_to(p.ROOT).as_posix()]=p.sha(file)
    item={'schema':VERSION,'config':c,'physical_source_config':pc,'source_artifacts':files,
        'source_manifests':manifests,'code_hashes':hashes,'interpreter':sys.executable,
        'python':sys.version,'threads':1,'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        'initial_snapshot':'results/_smoke/nlsp_construction_C_modal_sensitivity/initial_state.json'}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp
    if b.exists():raise ValueError('Preserve interrupted unmanifested attempt')
    b.mkdir(parents=True);p.write(b/'config.json',c);p.write(b/'provenance.json',item)
    archive=b/'execution_code';archive.mkdir()
    for file in hashes:shutil.copyfile(p.ROOT/file,archive/Path(file).name)
    shutil.copyfile(p.ROOT/item['initial_snapshot'],b/'initial_state.json')
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'stages':{k:'NOT_RUN' for k in
        ('source','geometry','modal_initial','base','gain','zero','spatial','temporal','finite','extension')},
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'static_Newton':0,'CCX':0,'Gmsh':0},
        'attempts':[],'cases':{},'numerical_seconds':0.,'deferred_theory':c['deferred_theory']}
    s['stages']['source']='PASS';p.finalize(b,item,s);return b,item,s


def limited_modes(d,plane,count):
    ids=d.plane if plane else d.out
    lam,vec=eigh(d.K[np.ix_(ids,ids)],d.M0[np.ix_(ids,ids)],
        subset_by_index=[0,min(count,len(ids))-1],check_finite=False)
    if lam.min()<=0:raise ArithmeticError('Reference K/M are not positive definite')
    return np.sqrt(lam),vec


def planar_bending_parts(d,q):
    """Twice V2 bending/shear by arm, from the accepted linear limit.

    The complementary axial/M-H part is q.K.q minus this sum. No new law.
    """
    parts=[]
    for arm,qi in zip(d.arms,d.local(q)):
        f=arm.reconstruct(qi);ds=arm.reconstruct(qi,derivative=1);c=arm.coefficients
        parts.append(float(arm.weights@(c.S*(ds[:,1]-f[:,5])**2+c.Bp*ds[:,5]**2)))
    return np.asarray(parts)


def select_modal_state(d,c):
    omega,vec=limited_modes(d,True,c['planar_modes_inspected']);rows=[];selected=None
    for j,w in enumerate(omega):
        q=np.zeros(d.ndof);q[d.plane]=vec[:,j];total=float(q@d.K@q);parts=planar_bending_parts(d,q)
        fraction=float(parts.sum()/total);maximum,points=p.polynomial_w_max(d,q)
        rows.append({'sorted_index_one_based':j+1,'omega':float(w),'bending_shear_fraction':fraction,
            'bending_shear_energy_by_arm':parts,'linear_potential_twice':total,'unit_mode_max_w':maximum})
        if selected is None and fraction>.5:selected=j
    if selected is None:raise ArithmeticError('No low planar bending-dominant mode')
    q=np.zeros(d.ndof);q[d.plane]=vec[:,selected]
    maximum,points=p.polynomial_w_max(d,q);q*=c['amplitude']/maximum
    # Reproducible sign: the largest initial local-w extremum is positive.
    _,points=p.polynomial_w_max(d,q)
    if max(points,key=lambda row:abs(row['value']))['value']<0:q=-q
    maximum,points=p.polynomial_w_max(d,q)
    T=float(2*np.pi/omega[selected]);report={'modes':rows,'selected_sorted_index_one_based':selected+1,
        'omega':float(omega[selected]),'Tref':T,'actual_max_w':maximum,'w_extrema':points,
        'normalization':'entire mass-normalized eigenvector scaled by exact represented local-w extrema',
        'zero_initial_velocity':True,'no_static_load_or_equilibrium':True,'frozen_before_nonlinear':True,
        'q_all_components':q,'eigen_residual_relative':float(np.linalg.norm(d.K@q-omega[selected]**2*d.M0@q)/np.linalg.norm(d.K@q))}
    return q,report


def displacement_subspace(d,c,omega_pl):
    omega,vec=limited_modes(d,False,c['out_modes_inspected']);fractions=[]
    for j,w in enumerate(omega):
        speed=np.zeros(d.ndof);speed[d.out]=vec[:,j]
        parts=d.reference_energy_parts(np.zeros(d.ndof),speed)
        fractions.append(float(parts[:,1].sum()/parts.sum()))
    bending=[j for j,f in enumerate(fractions) if f<.5][:2]
    torsion=next((j for j,f in enumerate(fractions) if f>=.5),None)
    if len(bending)!=2 or torsion is None:raise ArithmeticError('Declared bending/torsional space unavailable')
    selected=bending+[torsion];n=len(d.out);basis=np.zeros((2*n,len(selected)))
    for k,j in enumerate(selected):basis[:n,k]=vec[:,j]/omega[j]
    metric=r.fixed_metric(d);error=float(np.max(abs(basis.T@metric@basis-np.eye(len(selected)))))
    if error>2e-10:raise ArithmeticError('Initial displacement basis is not E0-orthonormal')
    return basis,metric,{'selected_sorted_indices_one_based':[j+1 for j in selected],
        'omega':omega,'omega_over_planar':omega/omega_pl,'torsional_mass_fractions':fractions,
        'initial_dimension':len(selected),'reference_gram_error':error,'zero_initial_velocities':True,
        'full_output_dimension':2*n,'restriction':'S_disp^(3), no initial velocity seeds'}


def physical_coordinate_max(d,q,fields=(2,3,4),Lref=1.):
    """Exact extrema of the represented Shen-plus-lift polynomials."""
    maxima=np.zeros(len(fields))
    for arm,qi in zip(d.arms,d.local(q)):
        for k,field in enumerate(fields):
            coeff=qi[field*arm.p:(field+1)*arm.p];leg=np.zeros(arm.p+1)
            for j,value in enumerate(coeff[:-1]):leg[j]+=value;leg[j+2]-=value
            leg[:2]+=coeff[-1]/2
            f=np.polynomial.legendre.Legendre(leg);roots=f.deriv().roots()
            xx=np.r_[-1.,1.,roots.real[(abs(roots.imag)<1e-9)&(abs(roots.real)<=1)]]
            maxima[k]=max(maxima[k],float(np.max(abs(f(xx))))*(Lref if field in (3,4) else 1))
    return maxima


def perturbation_scale(d,initial_out,c):
    q=np.zeros(d.ndof);q[d.out]=initial_out[:len(d.out)]
    maxima=physical_coordinate_max(d,q,Lref=c['Lref'])
    if maxima.max()<=0 or np.any(initial_out[len(d.out):]!=0):raise ValueError('Nonzero displacement-only perturbation required')
    delta=c['physical_perturbation_fraction']*c['amplitude']/float(maxima.max())
    return delta,{'delta_E0':delta,'target_max_v_LPhi_Lpsi':c['physical_perturbation_fraction']*c['amplitude'],
        'actual_initial_v_LPhi_Lpsi':delta*maxima,'Lref':c['Lref'],'all_initial_velocities_zero':True,
        'same_scale_for_all_components':True,'no_amplitude_fitting':True}


def limits(s,c,kind):
    if s['hard_stop']:raise RuntimeError('Historical attempt has hard-stopped; no retry')
    if s['numerical_seconds']>=c['numerical_budget_seconds']:raise TimeoutError('Numerical budget exhausted')
    if kind=='tangent' and s['calls']['tangent_ODE']>=c['maximum_tangent_trajectories']:raise TimeoutError('Tangent ceiling')
    if kind=='planar' and sum(x.get('plane',False) for x in s['cases'].values())>=c['maximum_planar_trajectories']:raise TimeoutError('Planar ceiling')
    if kind in ('planar','full') and s['calls']['nonlinear_ODE']>=c['maximum_nonlinear_trajectories']:raise TimeoutError('Nonlinear ceiling')


def full_attempt(b,item,s,label,d,q,times,deadline,plane):
    limits(s,item['config'],'planar' if plane else 'full');c=item['config']
    archive_orchestration(b,label)
    timing=p.settings(d,c['amplitude'],c)
    if not plane:
        timing['rtol']=c['full_rtol'];timing['atol']*=c['full_atol_factor'];timing['max_step']*=c['full_max_step_factor']
    initial=q if q.shape==(2*d.ndof,) else np.r_[q,np.zeros(d.ndof)]
    if not plane:
        base=p.read(b/'C_p16_base/integration.json');forecast_steps=base['accepted_steps']*c['full_atol_factor']**(-.25)
        p.write(b/'frozen_full_storage_policy.json',{'frozen_before_full_result':True,
            'estimated_accepted_steps':forecast_steps,'dense_archive_peak_estimate_GiB':3.632984*forecast_steps/122659,
            'native_state_file_estimate_MiB':(forecast_steps+1)*2*d.ndof*8/2**20,
            'storage':'all native accepted states streamed as float64; all native timestamps; 2001 full comparison states; native final state',
            'dense_coefficients_archived':False,'no_solver_method_tolerance_or_sampling_change':True,
            'reason':'avoid duplicated dense arrays exceeding the unchanged 4 GiB ceiling',
            'forecast_is_not_a_bound':True})
    with operator_progress(d,b,label,times[-1]):
        state,record=p.integrate_attempt(b,item,s,label,d,initial,times,timing,deadline,plane=plane,dense_archive=plane)
    if record['status']=='PASS':
        if plane:record['actual_min_step'],record['actual_max_step']=dense_step_extrema(b/label/'dense.npz')
        else:
            with np.load(b/label/'accepted_times.npz') as z:steps=np.diff(z['times'])
            record['actual_min_step'],record['actual_max_step']=float(steps.min()),float(steps.max())
        p.write(b/label/'integration.json',record);s['cases'][label]=record;p.finalize(b,item,s)
    return state


def dense_step_extrema(path):
    with np.load(path) as z:steps=z['t']-z['t_old']
    return float(steps.min()),float(steps.max())


def tangent_attempt(b,item,s,label,d,case,times,basis,metric,deadline,zero=False,tight=False):
    limits(s,item['config'],'tangent');c=item['config'];local=dict(item);tc=dict(c)
    archive_orchestration(b,label)
    factor=c['temporal_tolerance_factor'] if tight else 1.
    tc['rtol']=c['tangent_rtol']*factor;tc['atol_relative']=c['tangent_atol_relative']*factor;local['config']=tc
    max_step=p.settings(d,c['amplitude'],c)['max_step']*(c['temporal_max_step_factor'] if tight else 1)
    with operator_progress(d,b,label,times[-1]):
        Z,record=a.tangent_attempt(b,local,s,label,d,case,times,basis,metric,deadline,zero=zero,
            max_step=max_step,method=c['tangent_method'],sampled_dense_only=True)
    p.finalize(b,item,s);return Z


def archive_orchestration(b,label):
    folder=b/(label+'_orchestration');folder.mkdir()
    files=[Path(__file__),CONFIG]
    for file in files:shutil.copyfile(file,folder/file.name)
    p.write(folder/'hashes.json',{file.relative_to(p.ROOT).as_posix():p.sha(file) for file in files})
    p.write(folder/'actual_environment.json',{'interpreter':sys.executable,'python':sys.version,
        'numpy':np.__version__,'thread_environment':{key:os.environ.get(key) for key in
        ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')},'argv':sys.argv})


@contextmanager
def operator_progress(d,b,label,end):
    """Observe existing evaluations once a minute; never evaluate an extra RHS."""
    started=time.perf_counter();next_record=[started];original={key:getattr(d,key) for key in ('rhs','jacobian')}
    def wrap(fn):
        def observed(t,state):
            result=fn(t,state);now=time.perf_counter()
            if now>=next_record[0]:
                p.write(b/(label+'_progress.json'),{'last_operator_evaluation_time':float(t),'target_time':float(end),
                    'wall_seconds':now-started,'qualification':'operator evaluation, possibly trial stage; not an accepted native state'})
                next_record[0]=now+60
            return result
        return observed
    try:
        for key,fn in original.items():setattr(d,key,wrap(fn))
        yield
    finally:
        for key,fn in original.items():setattr(d,key,fn)


def preflight(b,item,s,model):
    c=item['config'];d=p.make_disc(c,'C',c['degree'],model)
    q,modal=select_modal_state(d,c);basis,metric,subspace=displacement_subspace(d,c,modal['omega'])
    geometry={'beta_deg':d.beta_deg,'lengths':d.lengths,'coefficients':[vars(arm.coefficients) for arm in d.arms],
        'Bp_over_Bb':[arm.coefficients.Bp/arm.coefficients.Bb for arm in d.arms],
        'bases_local_to_global':d.bases,'T7_local_to_global':d.T7,'same_joint_c':True,
        'rotation_contract':'a_local=(Phi,-psi,theta); k=-EZ; n along h, width along k'}
    if not np.allclose(geometry['Bp_over_Bb'],4,rtol=1e-14):raise ArithmeticError('C orientation incorrect')
    safety=d.safety(q);joint=d.joint_diagnostics(q,np.zeros(d.ndof));mass=np.linalg.eigvalsh(d.mass_matrix(q))
    for key in ('translation_difference','orientation_difference','contraction_difference'):
        if joint[key]>2e-12:raise ArithmeticError('Modal joint compatibility failed')
    p.write(b/'geometry.json',geometry);p.write(b/'frozen_modal_excitation.json',modal);p.write(b/'restricted_subspace.json',subspace)
    p.write(b/'initial_safety.json',{'arms':safety,'joint':joint,'mass_eigenvalue_min':float(mass.min()),
        'amplitude_over_h':.05,'potential':d.potential(q,gradient=False)['V'],'energy_rate_initial':d.energy_rate(q,np.zeros(d.ndof))})
    np.savez_compressed(b/'frozen_initial.npz',q=q,velocity=np.zeros(d.ndof),basis=basis,metric=metric)
    # Pointwise cost estimate only: no extra ODE, Newton, or eigen scan.
    x=np.r_[q,np.zeros(d.ndof)];bench={}
    for name,fn in [('rhs',lambda:d.rhs(0,x)),('jacobian',lambda:d.jacobian(0,x))]:
        samples=[]
        for _ in range(4):start=time.perf_counter();fn();samples.append(time.perf_counter()-start)
        bench[name]=float(np.median(samples[1:]))
    horizon=2*modal['Tref'];old=p.read(p.ROOT/c['parent_bundle']/'summary.json')['cases']
    base_rate=max(old['B_p16_base']['seconds']/old['B_p16_base']['target_end'],old['A_optimal_finite']['seconds']/old['A_optimal_finite']['target_end']/3)
    tangent_rate=max(old['A_p16_tangent']['seconds']/old['A_p16_tangent']['target_end'],old['B_p16_tangent']['seconds']/old['B_p16_tangent']['target_end'])
    plan={'frozen_before_new_ODE':True,'pointwise_seconds':bench,'physical_horizon':horizon,
        'base_p16_estimate':1.3*base_rate*horizon,'base_p12_estimate':1.3*base_rate*horizon,
        'main_tangent_each_estimate':1.2*tangent_rate*horizon,
        'optimal_tight_tangent_estimate':1.5*tangent_rate*horizon,
        'finite_full_estimate':1.2*max(old['A_optimal_finite']['seconds']/old['A_optimal_finite']['target_end'],old['B_optimal_finite']['seconds']/old['B_optimal_finite']['target_end'])*horizon,
        'zero_each_estimate':15.,'cost_source':'actual A/B counters; geometry/excitation dependent, estimates are not guarantees',
        'chosen_input_dimension':3,'controls_priority':'p12 and one temporal response before optional full/extension',
        'maximum_numerical_seconds':7200,'extension_is_conditional_not_reserved':True}
    plan['all_initial_program_estimate']=2*plan['base_p16_estimate']+2*plan['main_tangent_each_estimate']+plan['optimal_tight_tangent_estimate']+plan['finite_full_estimate']+30
    plan['mandatory_with_controls_estimate']=plan['all_initial_program_estimate']-plan['finite_full_estimate']
    if plan['mandatory_with_controls_estimate']>7200:raise TimeoutError('Required base/gain controls forecast exceeds budget before production')
    p.write(b/'frozen_resource_plan.json',plan)
    s['stages'].update(geometry='PASS',modal_initial='PASS');s['overall']='PREFLIGHT_COMPLETE';p.finalize(b,item,s)
    return d,q,basis,metric,modal


def restore_initial(b,item,model):
    d=p.make_disc(item['config'],'C',item['config']['degree'],model)
    with np.load(b/'frozen_initial.npz') as z:q,basis,metric=(z[k].copy() for k in ('q','basis','metric'))
    return d,q,basis,metric,p.read(b/'frozen_modal_excitation.json')


def save_base_checks(case,d,T):
    a.save_base_observations(case,d,T)
    with np.load(case/'trajectory.npz') as z:state=z['state'];energy=z['energy_relative_drift']
    rows=[]
    for index in np.linspace(0,len(state)-1,21,dtype=int):
        x=state[index];rows.append({'index':int(index),'safety':d.safety(x[:d.ndof]),
            'joint':d.joint_diagnostics(x[:d.ndof],x[d.ndof:]),'energy_rate':d.energy_rate(x[:d.ndof],x[d.ndof:])})
    p.write(case/'sampled_safety_joint.json',rows)
    invariance=p.read(case/'actual_planar_invariance.json')
    return {'energy_drift_max':float(np.max(abs(energy))),'sampled_checks':rows,
        'planar_invariance':invariance}


def save_gain(b,item,s,d,name,case,times,basis,metric,deadline):
    Z=tangent_attempt(b,item,s,name+'_tangent',d,case,times,basis,metric,deadline)
    if s['hard_stop']:return None
    zero=tangent_attempt(b,item,s,name+'_zero',d,case,times,basis,metric,deadline,zero=True)
    if s['hard_stop']:return None
    gram=np.einsum('tij,ik,tkl->tjl',zero,metric,zero);zero_error=float(np.max(abs(gram-np.eye(basis.shape[1]))))
    if zero_error>item['config']['zero_base_gram_gate']:raise ArithmeticError('Zero-base norm failed')
    result,gain=a.gain_result(d,Z,times,basis,metric);result.update(zero_gram_error=zero_error,time_at_max_Tref=result['time_at_max']/p.read(b/'frozen_modal_excitation.json')['Tref'])
    p.write(b/(name+'_gain.json'),result)
    delta,scale=perturbation_scale(d,np.asarray(result['initial_state_out']),item['config'])
    predicted=Z@np.asarray(result['initial_combination']);linear=zero@np.asarray(result['initial_combination'])
    histories=[]
    for value in (predicted,linear):
        x=np.zeros((len(times),2*d.ndof));x[:,d.state_out]=delta*value
        fields,_=a.physical_observations(d,x);histories.append(fields)
    physical={'scaling':scale,'max_v_Phi_psi':r.physical_field_max(histories[0]),
        'zero_base_same_input_max_v_Phi_psi':r.physical_field_max(histories[1]),
        'Rv':float(r.physical_field_max(histories[0])[0]/item['config']['amplitude']),
        'same_input_zero_base_Rv':float(r.physical_field_max(histories[1])[0]/item['config']['amplitude'])}
    p.write(b/(name+'_physical_response.json'),physical)
    np.savez_compressed(b/(name+'_gain_histories.npz'),times=times,gain=gain,
        zero_gain=np.sqrt(np.linalg.eigvalsh(gram)[:,-1]),optimal_normalized_out=predicted,
        zero_optimal_normalized_out=linear,optimal_physical_profiles=histories[0],zero_physical_profiles=histories[1])
    return result


def run(b,item,s,through='all'):
    if s['completed'] or s['hard_stop']:return s
    if any(row.get('status')!='PASS' for row in s.get('cases',{}).values()):
        raise RuntimeError('An incomplete actual attempt exists; no automatic retry or new direction')
    c=item['config'];model,_,_=p.sources(item['physical_source_config'])
    if not (b/'frozen_initial.npz').exists():d,q,basis,metric,modal=preflight(b,item,s,model)
    else:d,q,basis,metric,modal=restore_initial(b,item,model)
    if through=='preflight':return s
    T=modal['Tref'];times=p.make_times(T,2,c);deadline=time.perf_counter()+7200-s['numerical_seconds']
    if 'C_p16_base' not in s['cases']:
        full_attempt(b,item,s,'C_p16_base',d,q,times,deadline,plane=True)
        if s['hard_stop']:return s
        checks=save_base_checks(b/'C_p16_base',d,T);p.write(b/'C_p16_base_checks.json',checks)
        s['stages']['base']='PASS';p.finalize(b,item,s)
    if through=='base':return s
    if not (b/'C_p16_gain.json').exists():
        save_gain(b,item,s,d,'C_p16',b/'C_p16_base',times,basis,metric,deadline)
        if s['hard_stop']:return s
        s['stages'].update(gain='PASS',zero='PASS');p.finalize(b,item,s)
    if through=='gain':return s
    result=p.read(b/'C_p16_gain.json');direction=np.asarray(result['initial_combination'])
    if 'C_p16_optimal_tight' not in s['cases']:
        tight=tangent_attempt(b,item,s,'C_p16_optimal_tight',d,b/'C_p16_base',times,
            (basis@direction)[:,None],metric,deadline,tight=True)
        if s['hard_stop']:return s
        with np.load(b/'C_p16_gain_histories.npz') as z:old=z['optimal_normalized_out']
        curves=r.g.norm_curves(tight[:,:,0]-old,metric);error=curves['total']
        temporal={'E0_response_change_max':float(error.max()),'norm_max_old':float(a.energy_norm(old,metric).max()),
            'norm_max_tight':float(a.energy_norm(tight[:,:,0],metric).max()),'change_time_at_max':float(times[np.argmax(error)]),
            'coordinate_change_max':float(curves['coordinate'].max()),'velocity_change_max':float(curves['velocity'].max()),
            'velocity_fraction_at_max':float(curves['velocity_fraction'][np.argmax(error)]),
            'scope':'one fixed optimal input, not a re-optimization of the whole tight propagator',
            'baseline_trajectory_temporal_refinement':'NOT_RUN'}
        np.savez_compressed(b/'temporal_sensitivity_arrays.npz',times=times,**curves)
        p.write(b/'temporal_sensitivity.json',temporal);s['stages']['temporal']='PASS';p.finalize(b,item,s)
    if not (b/'spatial_sensitivity.json').exists():
        d12=p.make_disc(c,'C',12,model);q12=r.projected_coordinates(d,d12,q);E12=r.fixed_metric(d12)
        B12,projection=r.project_basis(d,d12,basis,E12)
        p.write(b/'spatial_projection.json',projection)
        np.savez_compressed(b/'p12_projected_initial.npz',q=q12,basis=B12,metric=E12)
        if 'C_p12_base' not in s['cases']:
            full_attempt(b,item,s,'C_p12_base',d12,q12,times,deadline,plane=True)
            if s['hard_stop']:return s
            checks=save_base_checks(b/'C_p12_base',d12,T);p.write(b/'C_p12_base_checks.json',checks)
        save_gain(b,item,s,d12,'C_p12',b/'C_p12_base',times,B12,E12,deadline)
        if s['hard_stop']:return s
        spatial=compare_spatial_saved(b,d,d12,times);p.write(b/'spatial_sensitivity.json',spatial)
        s['stages']['spatial']='PASS';p.finalize(b,item,s)
    qualification=qualify_saved(b,c);p.write(b/'signal_qualification.json',qualification)
    if through=='controls':return s
    # Extension evaluated before any new long trajectory; controls/full have priority.
    if not (b/'extension_decision.json').exists():
        with np.load(b/'C_p16_gain_histories.npz') as z:gain=z['gain']
        continuing=bool(np.argmax(gain)>=len(times)-201 and gain[-1]>np.max(gain[:-201]))
        forecast=1.5*(s['cases']['C_p16_base']['seconds']+s['cases']['C_p16_tangent']['seconds'])
        # A longer sensitive direction would require the single full control
        # on the longer horizon as well. Reserve 5T, not just its 2T prefix.
        full_forecast=p.read(b/'after_gain_resource_decision.json')['full_tight_estimate'] if (b/'after_gain_resource_decision.json').exists() else p.read(b/'frozen_resource_plan.json')['finite_full_estimate']
        reserve=2.5*full_forecast if qualification['gain_resolved_in_observed_controls'] else 0
        affordable=forecast+reserve<7200-s['numerical_seconds']
        decision={'resolved':qualification['gain_resolved_in_observed_controls'],'new_maxima_in_last_point2T':continuing,
            'extension_estimate_seconds':forecast,'finite_reservation_seconds':reserve,'remaining_seconds':7200-s['numerical_seconds'],
            'affordable':affordable,'extended':False,'status':'NOT_RUN',
            'reason':'no resolved continuing growth or insufficient reserved budget; two-period experiment retained'}
        # A positive gate is not silently ignored: a continuation would need the
        # same actual final states and a dedicated t_start contract, never restart.
        if qualification['gain_resolved_in_observed_controls'] and continuing and affordable:
            decision['reason']='resolved continuing growth, but current reused IVP interface begins at zero; no unsafe restart; scoped result limited to 2T'
        p.write(b/'extension_decision.json',decision)
    if qualification['gain_resolved_in_observed_controls'] and 'C_optimal_finite' not in s['cases']:
        forecast=p.read(b/'after_gain_resource_decision.json')['full_tight_estimate'] if (b/'after_gain_resource_decision.json').exists() else p.read(b/'frozen_resource_plan.json')['finite_full_estimate']
        if forecast<7200-s['numerical_seconds']:
            delta,scale=perturbation_scale(d,np.asarray(result['initial_state_out']),c);initial=np.r_[q,np.zeros(d.ndof)]
            initial[d.state_out]+=delta*np.asarray(result['initial_state_out']);p.write(b/'finite_initial_scale.json',scale)
            full=full_attempt(b,item,s,'C_optimal_finite',d,initial,times,deadline,plane=False)
            if s['hard_stop']:return s
            with np.load(b/'C_p16_base/trajectory.npz') as z:base=z['state']
            with np.load(b/'C_p16_tangent/tangent.npz') as z:Z=z['vectors']
            finite=a.full_tangent_error(d,full,base,Z,direction,delta,metric)
            p.write(b/'finite_comparison.json',{k:v for k,v in finite.items() if k not in ('predicted','observed')})
            p.write(b/'finite_error_details.json',a.response_error_details(d,finite['observed']-finite['predicted'],delta,metric))
            np.savez_compressed(b/'finite_comparison_arrays.npz',times=times,observed=finite['observed'],predicted=finite['predicted'])
            s['stages']['finite']='PASS' if finite['normalized_out_error_relative']<=c['tangent_full_relative_gate'] else 'PARTIAL'
        else:p.write(b/'finite_NOT_RUN.json',{'reason':'forecast exceeds remaining numerical budget','seconds_estimate':forecast})
    s['overall']='C_MODAL_SENSITIVITY_RESOLVED' if qualification['gain_resolved_in_observed_controls'] else 'C_MODAL_SENSITIVITY_NUMERICAL_PARTIAL'
    if s['stages']['finite']=='PARTIAL':s['overall']='C_MODAL_SENSITIVITY_NUMERICAL_PARTIAL'
    s['completed']=True;p.finalize(b,item,s);return s


def compare_spatial_saved(b,d16,d12,times):
    with np.load(b/'C_p16_base/trajectory.npz') as z:X16=z['state']
    with np.load(b/'C_p12_base/trajectory.npz') as z:X12=z['state']
    f16,_=a.physical_observations(d16,X16);f12,_=a.physical_observations(d12,X12)
    gain16=p.read(b/'C_p16_gain.json');gain12=p.read(b/'C_p12_gain.json')
    with np.load(b/'C_p16_gain_histories.npz') as z:g16=z['gain'];p16=z['optimal_physical_profiles']
    with np.load(b/'C_p12_gain_histories.npz') as z:g12=z['gain'];p12=z['optimal_physical_profiles']
    # Compare the very same p16 optimal physical initial state after projection.
    # Re-optimized eigenvectors have arbitrary signs and are not a same-input check.
    initial=np.asarray(gain16['initial_state_out']);q0=np.zeros(d16.ndof);q0[d16.out]=initial[:len(d16.out)]
    projected=np.zeros(len(d12.state_out));projected[:len(d12.out)]=r.projected_coordinates(d16,d12,q0)[d12.out]
    with np.load(b/'C_p12_tangent/tangent.npz') as z:Z12,B12,E12=(z[k] for k in ('vectors','basis','metric'))
    weights=B12.T@E12@projected;projection_residual=a.energy_norm(projected-B12@weights,E12)
    response12=Z12@weights
    delta=p.read(b/'C_p16_physical_response.json')['scaling']['delta_E0']
    embedded=np.zeros((len(times),2*d12.ndof));embedded[:,d12.state_out]=delta*response12
    same12,_=a.physical_observations(d12,embedded)
    velocity16,_=a.physical_observations(d16,np.c_[X16[:,d16.ndof:],np.zeros_like(X16[:,d16.ndof:])])
    velocity12,_=a.physical_observations(d12,np.c_[X12[:,d12.ndof:],np.zeros_like(X12[:,d12.ndof:])])
    np.savez_compressed(b/'spatial_comparison_arrays.npz',times=times,base_p16=f16,base_p12=f12,
        optimal_p16=p16,optimal_p12=p12,same_initial_optimal_p12=same12,
        base_velocity_p16=velocity16,base_velocity_p12=velocity12)
    return {'gain_max_p16':gain16['gain'],'gain_max_p12':gain12['gain'],
        'gain_max_change':abs(gain16['gain']-gain12['gain']),'gain_curve_change_max':float(np.max(abs(g16-g12))),
        'time_max_Tref_p16':gain16['time_at_max_Tref'],'time_max_Tref_p12':gain12['time_at_max_Tref'],
        'base_field_difference_max':np.max(abs(f16-f12),axis=(0,1,2)),
        'base_velocity_difference_max':np.max(abs(velocity16-velocity12),axis=(0,1,2)),
        'optimal_physical_difference_max':np.max(abs(p16-p12),axis=(0,1,2)),
        'same_initial_physical_response_difference_max':np.max(abs(p16-same12),axis=(0,1,2)),
        'same_initial_p12_projection_residual_E0':float(projection_residual),
        'optimal_response_comparison':'separately optimized directions; no phase or amplitude fitting',
        'same_physical_input_span':'nested projection; each level uses its own fixed E0 norm',
        'p12_excitation':'physical projection of the same p16 state, no new modal excitation amplitude'}


def qualify_saved(b,c):
    gain=p.read(b/'C_p16_gain.json')['gain'];spatial=p.read(b/'spatial_sensitivity.json');temporal=p.read(b/'temporal_sensitivity.json')
    signal=gain-1;dh=spatial['gain_curve_change_max'];dt=temporal['E0_response_change_max'];factor=c['signal_over_observed_controls_planning_factor']
    return {'gain_minus_one':signal,'spatial_gain_change':dh,'temporal_optimal_response_change':dt,
        'signal_over_spatial_change':signal/max(dh,1e-300),'signal_over_temporal_change':signal/max(dt,1e-300),
        'planning_factor_fixed_before_results':factor,'gain_resolved_in_observed_controls':bool(signal>factor*max(dh,dt,c['numerical_floor_relative'])),
        'qualification':'observed p12/p16 and one tangent temporal control; not strict error bounds or infinite-dimensional stability',
        'base_temporal_control':'NOT_RUN','all_seven_field_spatial_certification':'NOT_ESTABLISHED'}


def export_saved(b):
    """Physical reconstruction from saved states only, including real velocities."""
    b=Path(b)
    if (b/'export_provenance.json').exists():return p.read(b/'export_provenance.json')
    item=p.read(b/'provenance.json');s=p.read(b/'summary.json');model,_,_=p.sources(item['physical_source_config'])
    d=p.make_disc(item['config'],'C',16,model);T=p.read(b/'frozen_modal_excitation.json')['Tref']
    rows=[]
    for label,record in s['cases'].items():
        rows.append({'case':label,**{key:record.get(key) for key in ('status','seconds','accepted_steps','nfev','njev','nlu','actual_min_step','actual_max_step','actual_accepted_end','energy_relative_drift_max')}})
    if rows:
        with (b/'actual_execution.csv').open('w',newline='',encoding='utf8') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for label in ('C_p16_base','C_optimal_finite'):
        case=b/label
        if not (case/'trajectory.npz').exists():continue
        with np.load(case/'trajectory.npz') as z:times,state,drift=(z[k] for k in ('times','state','energy_relative_drift'))
        f,obs=a.physical_observations(d,state)
        velocity,velobs=a.physical_observations(d,np.c_[state[:,d.ndof:],np.zeros_like(state[:,d.ndof:])])
        np.savez_compressed(b/(label+'_physical.npz'),times=times,tau=times/T,x=np.linspace(0,1,201),
            profiles=f,observations=obs,velocity_profiles=velocity,velocity_observations=velobs,energy_drift=drift)
        metrics={'max_abs_fields_u_w_v_Phi_psi_theta_c':np.max(abs(f),axis=(0,1,2)),
            'max_abs_velocities_u_w_v_Phi_psi_theta_c':np.max(abs(velocity),axis=(0,1,2)),
            'end_max_abs_fields':np.max(abs(f[-1]),axis=(0,1)),
            'sampled_grid':'2001 actual requested times; 201 physical coordinates per arm',
            'sampled_safety':[d.safety(x[:d.ndof]) for x in state[::100]],
            'joint_checks':[d.joint_diagnostics(state[j,:d.ndof],state[j,d.ndof:]) for j in (0,len(times)//2,len(times)-1)]}
        p.write(b/(label+'_physical_metrics.json'),metrics)
    for name in ('C_p16','C_p12'):
        if not (b/(name+'_gain_histories.npz')).exists():continue
        with np.load(b/(name+'_gain_histories.npz')) as z:times,gain,zero,profile=(z[k] for k in ('times','gain','zero_gain','optimal_physical_profiles'))
        amplitude=item['config']['amplitude']
        with (b/(name+'_gain_observations.csv')).open('w',newline='',encoding='utf8') as stream:
            writer=csv.writer(stream);writer.writerow(['time','time_over_Tref','restricted_Gdisp','zero_Gdisp','max_abs_v','max_abs_Phi','max_abs_psi','Rv'])
            for j,t in enumerate(times):
                values=r.physical_field_max(profile[j]);writer.writerow([t,t/T,gain[j],zero[j],*values,values[0]/amplitude])
        with np.load(b/(name+'_gain_histories.npz')) as z:optimal=z['optimal_normalized_out']
        with np.load(b/(name+'_tangent/tangent.npz')) as z:metric=z['metric']
        n=optimal.shape[1]//2;coordinate=a.energy_norm(optimal[:,:n],metric[:n,:n]);velocity=a.energy_norm(optimal[:,n:],metric[n:,n:])
        result=p.read(b/(name+'_gain.json'));index=result['sample_index_at_max']
        p.write(b/(name+'_fixed_optimal_energy.json'),{'coordinate_fraction_at_global_gain_max':float(coordinate[index]**2/result['gain']**2),
            'velocity_fraction_at_global_gain_max':float(velocity[index]**2/result['gain']**2),
            'reference_energy_parts_by_arm_bending_torsion':result['reference_energy_distribution_at_max_by_arm_bending_torsion'],
            'qualification':'fixed E0 perturbation norm; not a separate conserved nonlinear energy'})
        np.savez_compressed(b/(name+'_fixed_optimal_energy_histories.npz'),times=times,coordinate_norm=coordinate,velocity_norm=velocity)
    evidence={'scientific_calls':0,'physical_reconstruction_only':True,'native_velocity_state_used':True,
        'source_summary_sha256':p.sha(b/'summary.json'),'no_phase_amplitude_fitting':True,'all_OUT_components_retained':True}
    p.write(b/'export_provenance.json',evidence);return evidence


def plot_saved(b):
    """Four figures, saved-data only; no IVP, eigenanalysis or smoothing."""
    b=Path(b)
    if (b/'figure_provenance.json').exists():return p.read(b/'figure_provenance.json')['figures']
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    figures=[];T=p.read(b/'frozen_modal_excitation.json')['Tref'];c=p.read(b/'config.json')
    def finish(fig,name):
        fig.tight_layout()
        for suffix in ('pdf','png'):fig.savefig(b/(name+'.'+suffix),dpi=180,bbox_inches='tight');figures.append(name+'.'+suffix)
        plt.close(fig)
    with np.load(b/'C_p16_base_physical.npz') as z:f=z['profiles'];tau=z['tau'];x=z['x'];drift=z['energy_drift']
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    for field,ax,label in zip((0,1,5,6),axes.flat,('u / reference length','w / reference length','theta [rad]','c [1]')):
        for arm in range(2):ax.plot(x,f[0,arm,:,field],label=f'arm {arm+1}')
        ax.set(xlabel='s / Li',ylabel=label);ax.grid(alpha=.3);ax.legend()
    fig.suptitle('C: complete planar mode, max |w| = 0.005; all initial velocities zero');finish(fig,'C_modal_initial_fields')
    with np.load(b/'C_p16_gain_histories.npz') as z:out=z['optimal_physical_profiles'];zero=z['zero_physical_profiles'];gain=z['gain'];gz=z['zero_gain']
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    for arm in range(2):axes[0,0].plot(tau,f[:,arm,100,1],label=f'w, arm {arm+1} midspan')
    axes[0,0].set(ylabel='w / reference length');axes[0,0].legend()
    for k,field in enumerate((2,3,4)):
        ax=axes.flat[k+1]
        for arm in range(2):ax.plot(tau,out[:,arm,100,field],label=f'arm {arm+1}, s/Li=0.5')
        ax.set(ylabel=('v / reference length','Phi [rad]','psi [rad]')[k]);ax.legend()
    for ax in axes.flat:ax.set_xlabel('t / Tref');ax.grid(alpha=.3)
    fig.suptitle('C: nonlinear planar base and tangent OUT prediction; 1% initial deformation convention');finish(fig,'C_planar_and_OUT_motion')
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True)
    axes[0].plot(tau,gain,label='p16, S_disp^(3)');axes[0].plot(tau,gz,'--',label='p16 zero base')
    if (b/'C_p12_gain_histories.npz').exists():
        with np.load(b/'C_p12_gain_histories.npz') as z:axes[0].plot(z['times']/T,z['gain'],':',label='p12, projected same initial span')
    axes[0].set(ylabel='restricted energy gain');axes[0].legend()
    axes[1].plot(tau,np.max(abs(out[:,:,:,2]),axis=(1,2))/c['amplitude'],label='nonlinear base, fixed optimal input')
    axes[1].plot(tau,np.max(abs(zero[:,:,:,2]),axis=(1,2))/c['amplitude'],'--',label='zero base, SAME input')
    axes[1].set(xlabel='t / Tref',ylabel='max |v| / Apl');axes[1].legend()
    for ax in axes:ax.grid(alpha=.3)
    finish(fig,'C_restricted_gain_and_zero_control')
    fig,axes=plt.subplots(2,2,figsize=(10,6))
    for arm in range(2):
        axes[0,0].plot(x,out[0,arm,:,2],label=f'initial v, arm {arm+1}')
        axes[0,1].plot(x,out[0,arm,:,3],label=f'initial Phi, arm {arm+1}')
        axes[0,1].plot(x,out[0,arm,:,4],'--',label=f'initial psi, arm {arm+1}')
    axes[0,0].set(xlabel='s / Li',ylabel='v / reference length');axes[0,1].set(xlabel='s / Li',ylabel='rotation [rad]')
    axes[1,0].plot(tau,drift,label='p16 base raw energy drift')
    if (b/'C_optimal_finite_physical.npz').exists():
        with np.load(b/'C_optimal_finite_physical.npz') as z:finite=z['profiles'];fd=z['energy_drift'];finite_tau=z['tau']
        axes[1,0].plot(finite_tau,fd,label='full seven-field raw drift')
        for arm in range(2):
            axes[1,1].plot(tau,out[:,arm,100,2],label=f'tangent arm {arm+1}')
            axes[1,1].plot(finite_tau,finite[:,arm,100,2],'--',label=f'full arm {arm+1}')
        axes[1,1].set(ylabel='v(s/Li=0.5) / reference length')
    elif (b/'C_p16_optimal_tight/tangent.npz').exists():
        with np.load(b/'C_p16_optimal_tight/tangent.npz') as z:tight=z['vectors'][:,:,0]
        with np.load(b/'C_p16_gain_histories.npz') as z:old=z['optimal_normalized_out']
        with np.load(b/'frozen_initial.npz') as z:E=z['metric']
        axes[1,1].plot(tau,a.energy_norm(tight-old,E),label='tight/old fixed optimal response')
        axes[1,1].set(ylabel='E0 norm difference')
    axes[1,0].set(ylabel='(E-E0) / E0')
    for ax in axes[1]:ax.set_xlabel('t / Tref')
    for ax in axes.flat:ax.legend(fontsize=8);ax.grid(alpha=.3)
    finish(fig,'C_optimal_form_and_controls')
    p.write(b/'figure_provenance.json',{'figures':figures,'scientific_calls':0,'smoothing':False,'fitting':False,
        'OUT_curves_are_tangent_prediction_except_explicit_full_label':True})
    return figures


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--compute',action='store_true');group.add_argument('--report-only',type=Path);group.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--through-stage',choices=('preflight','base','gain','controls','all'),default='all')
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=cache(b)
    else:
        b,item,s=prepare(args.config)
        try:s=run(b,item,s,args.through_stage)
        except (ArithmeticError,ValueError,TimeoutError,RuntimeError,np.linalg.LinAlgError) as error:
            s.update(hard_stop=True,overall='BLOCKED',stop_reason=type(error).__name__+': '+str(error));p.finalize(b,item,s)
    if s.get('completed') and not (b/'export_provenance.json').exists():
        export_saved(b);plot_saved(b);p.finalize(b,p.read(b/'provenance.json'),s)
    print(json.dumps(p.plain({'bundle':str(b),'overall':s['overall'],'stages':s['stages'],'calls':s['calls'],'numerical_seconds':s['numerical_seconds'],'stop_reason':s.get('stop_reason')}),ensure_ascii=False,indent=2))
