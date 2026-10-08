"""Prepared-IC audit from immutable spectral data; no historical reruns.

The 1e-6 profile/jet admission gate precedes all nonlinear trajectories.
This workflow preserves the old zero-axial initial-value problem separately.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

if __name__ == '__main__':
    for key in ('OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'OMP_NUM_THREADS'):
        os.environ[key] = '1'
import numpy as np
from numpy.polynomial import legendre as leg

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.analysis import verify_planar_second_order_axial_response as previous
sha, read_json, write_json, save_npz = previous.sha, previous.read_json, previous.write_json, previous.save_npz
CONFIG = ROOT/'data/input/planar_prepared_initial_state.json'
OUTPUT = ROOT/'results/planar_prepared_initial_state'
VERSION = 'nlsp-prepared-ic-profiles-and-endpoint-audit-v1'


def identity(config_path=CONFIG):
    config = read_json(config_path)
    paths = ('scripts/analysis/prepare_planar_initial_state.py',
             'scripts/lib/planar_prepared_initial_state.py',
             'scripts/lib/planar_second_order_axial_response.py',
             'scripts/analysis/verify_planar_second_order_axial_response.py',
             'scripts/lib/weakly_nonlinear_planar_dynamics.py',
             'scripts/lib/weakly_nonlinear_spatial_rod.py',
             'scripts/lib/mindlin_herrmann_longitudinal.py',
             'scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py',
             'scripts/analysis/simulate_weakly_nonlinear_planar_rod.py',
             config['pilot_config'])
    pilot = read_json(ROOT/config['pilot_config'])
    source_paths = {**config['historical_bundles'], 'spectral':config['spectral_bundle'],
                    'audit':pilot['audit_bundle'], 'linear':pilot['linear_reference_bundle']}
    source = {}
    for name, path in source_paths.items():
        manifest = ROOT/path/'manifest.json'
        source[name] = {'bundle':path, 'manifest_sha256':sha(manifest) if manifest.exists() else None}
        if name in ('audit','linear'):
            target = ROOT/path/'result.json'
            source[name]['result_sha256'] = sha(target) if target.exists() else None
    value = {'version':VERSION, 'config':config, 'config_sha256':sha(config_path),
             'code_hashes':{p:sha(ROOT/p) for p in paths}, 'sources':source,
             'python':sys.version, 'dependencies':{n:importlib.metadata.version(n) for n in ('numpy','scipy','matplotlib')},
             'blas_threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')}}
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()[:16], value


def validate_cache(bundle, expected=None):
    bundle = Path(bundle)
    manifest = read_json(bundle/'manifest.json')
    if expected is not None and manifest['identity'] != expected:
        raise ValueError('Prepared-IC cache identity mismatch')
    for path, digest in manifest['artifact_hashes'].items():
        if sha(bundle/path) != digest:
            raise ValueError('Prepared-IC artifact hash mismatch: '+path)
    return read_json(bundle/'summary.json')


def profile_comparison(low, high, length, h0, policy):
    """Exact physical Legendre Gram norms and independent derivative jets."""
    degree = max(low.shape[-1], high.shape[-1])
    a, b = np.zeros((2,degree)), np.zeros((2,degree))
    a[:,:low.shape[-1]], b[:,:high.shape[-1]] = low, high
    xi, _ = leg.leggauss(100)
    rows = []
    for d in policy['spatial_derivatives']:
        for f, field in enumerate(('u','c')):
            aa, bb = leg.legder(a[f],m=d)*(2/length)**d, leg.legder(b[f],m=d)*(2/length)**d
            dif = aa-bb
            gram = length/(2*np.arange(len(dif))+1)
            abs_l2 = float(np.sqrt(np.sum(dif*dif*gram)))
            abs_max = float(np.max(abs(leg.legval(xi,dif))))
            ref_l2 = float(np.sqrt(np.sum(bb*bb*gram)))
            ref_max = float(np.max(abs(leg.legval(xi,bb))))
            scale = (h0 if field=='u' else 1.)/length**d
            floor = 1e-10*scale
            jets_a, jets_b = leg.legval([-1.,1.],aa), leg.legval([-1.,1.],bb)
            row = {'field':field,'derivative':d, 'absolute_L2':abs_l2, 'absolute_max_Gauss100':abs_max,
                   'reference_L2':ref_l2, 'reference_max_Gauss100':ref_max,
                   'relative_L2':abs_l2/max(ref_l2,floor*np.sqrt(length)),
                   'relative_max':abs_max/max(ref_max,floor),
                   'fixed_scaled_L2':abs_l2/(scale*np.sqrt(length)), 'fixed_scaled_max':abs_max/scale,
                   'fixed_scale':scale, 'endpoint_low':jets_a.tolist(),'endpoint_high':jets_b.tolist(),
                   'endpoint_absolute_difference':abs(jets_a-jets_b).tolist(),
                   'endpoint_fixed_scaled_difference':float(np.max(abs(jets_a-jets_b))/scale)}
            tol, etol = policy['relative_tolerance'], policy['endpoint_tolerance']
            row['pass'] = all(row[k]<=tol for k in ('relative_L2','relative_max','fixed_scaled_L2','fixed_scaled_max')) and row['endpoint_fixed_scaled_difference']<=etol
            rows.append(row)
    return {'rows':rows, 'pass':all(r['pass'] for r in rows), 'max_is_sampled':True,
            'derivatives_from_Legendre_not_PDE':True}


def plot_only(bundle):
    """Render saved arrays only; no preparation or time/eigen evaluations."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bundle = Path(bundle)
    validate_cache(bundle)
    if not (bundle/'profiles.npz').exists():
        return {'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0}
    data = np.load(bundle/'profiles.npz')
    x = data['s']
    fig, axs = plt.subplots(1,2,figsize=(9,3.4),layout='constrained')
    for f, field in enumerate(('u','c')):
        for part, style in (('stat','-'),('harm','--')):
            axs[f].plot(x,data['p96_'+part][:,f],style,label=part)
        axs[f].set(xlabel='s/L',ylabel=field+' (second-order coefficient)')
        axs[f].legend(frameon=False)
        axs[f].grid(alpha=.2)
    figs=bundle/'figures';figs.mkdir(exist_ok=True)
    for ext in ('pdf','png'):
        fig.savefig(figs/('periodic_profiles.'+ext),dpi=220,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
    plt.close(fig)
    if not (bundle/'profile_convergence.json').exists():
        return {'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0}
    convergence = read_json(bundle/'profile_convergence.json')
    fig,axs = plt.subplots(1,2,figsize=(9,3.4),layout='constrained')
    for f,field in enumerate(('u','c')):
        for part,style in (('stat','o-'),('harm','s--')):
            for d in (0,1,2):
                y=[];xx=[]
                for pair in convergence['pairs']:
                    row=next(r for r in pair[part]['rows'] if r['field']==field and r['derivative']==d)
                    y.append(max(row['relative_L2'],1e-18));xx.append(pair['high_p'])
                axs[f].semilogy(xx,y,style,label=f'{part}, d={d}')
        axs[f].axhline(1e-6,color='black',ls=':',lw=.8)
        axs[f].set(xlabel='higher p in neighboring comparison',ylabel=field+' relative L2 difference')
        axs[f].legend(frameon=False,fontsize=7)
        axs[f].grid(alpha=.2)
    for ext in ('pdf','png'):
        fig.savefig(figs/('profile_derivative_convergence.'+ext),dpi=220,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
    plt.close(fig)
    return {'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0}


def numeric_endpoint_audit(state, epsilon, coefficients):
    """Independently reconstructed jets, with exact essential-value contract."""
    end = np.array([0.,state.background.length])
    first, second = state.evaluate(end,epsilon,1,require_admitted=False), state.evaluate(end,epsilon,2,require_admitted=False)
    u_s,w_s,theta_s,c_s = first.T
    u_ss,w_ss,theta_ss,c_ss = second.T
    p=coefficients
    numerators = np.column_stack((p.C*u_ss+p.nu*p.C*c_s+(p.C-p.S)*theta_s*w_s,
                                  p.S*(w_ss-theta_s)+(p.C-p.S)*u_s*theta_s,
                                  p.Bp*theta_ss+p.S*w_s-(p.C-p.S)*u_s*w_s,
                                  p.H*c_ss-p.nu*p.C*u_s))
    mass=np.array([p.m,p.m,p.jp,p.jp])
    return {'epsilon_a':epsilon,'first_derivatives':first.tolist(),'second_derivatives':second.tolist(),
            'acceleration_numerators':numerators.tolist(),'accelerations':(numerators/mass).tolist(),
            'essential_values_actual':state.evaluate(end,epsilon,require_admitted=False).tolist(),
            'velocities':[0.,0.,0.,0.], 'finite_amplitude_residual_not_truncated':True}


def project_initial(state, disc, epsilon, policy):
    initial=state.evaluate(disc.x,epsilon)
    q=disc.project(initial)
    xi,ww=leg.leggauss(max(100,2*disc.p+1));x=(xi+1)*disc.length/2;ww*=disc.length/2
    rows=[]
    fixed=[epsilon*state.background.h0,epsilon*state.background.h0,
           epsilon*state.background.h0/disc.length,epsilon*state.background.h0/disc.length]
    for d in policy['spatial_derivatives']:
        ref=state.evaluate(x,epsilon,d)
        val=disc.reconstruct(q,x,d)
        jets=disc.reconstruct(q,[0.,disc.length],d)-state.evaluate(np.array([0.,disc.length]),epsilon,d)
        for f,field in enumerate(('u','w','theta','c')):
            dif=val[:,f]-ref[:,f]; scale=fixed[f]/disc.length**d
            norm=float(np.sqrt(np.dot(ww,dif*dif))); peak=float(np.max(abs(dif)))
            rnorm=float(np.sqrt(np.dot(ww,ref[:,f]**2))); rmax=float(np.max(abs(ref[:,f])))
            row={'field':field,'derivative':d,'absolute_L2':norm,'absolute_max':peak,
                 'relative_L2':norm/max(rnorm,1e-10*scale*np.sqrt(disc.length)),
                 'relative_max':peak/max(rmax,1e-10*scale),
                 'endpoint_fixed_scaled_error':float(np.max(abs(jets[:,f]))/scale)}
            row['pass']=row['relative_L2']<=policy['relative_tolerance'] and row['relative_max']<=policy['relative_tolerance'] and row['endpoint_fixed_scaled_error']<=policy['endpoint_tolerance']
            rows.append(row)
    # Formal through-cubic coefficient traces after projection require projecting
    # the common unscaled components separately. No derivative is inferred by PDE.
    bg=state.background
    def full(fields):
        return disc.project(fields)
    o1=np.zeros((disc.nq,4));o2=o1.copy();o3=o1.copy()
    pair=bg.evaluate(disc.x); uc=state.profiles.evaluate(disc.x)
    o1[:,1:3]=pair; o2[:,0]=uc[:,0];o2[:,3]=uc[:,1]
    o3[:,2]=state.correction.evaluate(disc.x)
    q1,q2,q3=map(full,(o1,o2,o3))
    d1,d2,d3=[disc.reconstruct(v,[0.,disc.length],1) for v in (q1,q2,q3)]
    dd1,dd2,dd3=[disc.reconstruct(v,[0.,disc.length],2) for v in (q1,q2,q3)]
    c=disc.coefficients
    e1=np.column_stack((np.zeros(2),c.S*(dd1[:,1]-d1[:,2]),c.Bp*dd1[:,2]+c.S*d1[:,1],np.zeros(2)))
    e2=np.column_stack((c.C*dd2[:,0]+c.nu*c.C*d2[:,3]+(c.C-c.S)*d1[:,2]*d1[:,1],np.zeros(2),np.zeros(2),c.H*dd2[:,3]-c.nu*c.C*d2[:,0]))
    e3=np.column_stack((np.zeros(2),-c.S*d3[:,2]+(c.C-c.S)*d2[:,0]*d1[:,2],c.Bp*dd3[:,2]-(c.C-c.S)*d2[:,0]*d1[:,1],np.zeros(2)))
    scales=np.array([c.C*bg.h0/disc.length**2,c.S*bg.h0/disc.length**2,c.S*bg.h0/disc.length,c.C])
    scaled=float(max(np.max(abs(e)/scales) for e in (e1,e2,e3)))
    diag=mass_safety_bounds(disc,q)
    energy=disc.energy(q,np.zeros(disc.ndof))
    result={'p':disc.p,'rows':rows,'through_cubic_coefficients':{'epsilon1':e1.tolist(),'epsilon2':e2.tolist(),'epsilon3':e3.tolist()},
            'compatibility_fixed_scaled_max':scaled,'initial_energy':float(energy), 'safety':diag,
            'pass':all(r['pass'] for r in rows) and scaled<=policy['endpoint_tolerance'],
            'common_reference_unchanged':True}
    return q,result


def mass_safety_bounds(disc,q):
    """Exact Loewner bounds for the variable theta Gram, without eigensolves."""
    f=disc.reconstruct(q);g=disc.reconstruct(q,derivative=1)
    lo=min(1.,float(np.min((1+f[:,3])**2)));hi=max(1.,float(np.max((1+f[:,3])**2)))
    return {'min_one_plus_c':float(np.min(1+f[:,3])), 'max_abs_c':float(np.max(abs(f[:,3]))),
            'max_abs_theta':float(np.max(abs(f[:,2]))), 'max_abs_axial_gradient':float(np.max(abs(g[:,0]))),
            'max_abs_transverse_gradient':float(np.max(abs(g[:,1]))),
            'max_L_abs_curvature':float(disc.length*np.max(abs(g[:,2]))),
            'relative_mass_eigenvalue_lower_bound':lo,'relative_mass_eigenvalue_upper_bound':hi,
            'relative_mass_condition_upper_bound':hi/lo,'mass_positive':lo>0,
            'method':'weighted Gram Loewner bounds; no eigenvalue decomposition'}



def stopped_preparation(config,bundle,started,reason,details=None):
    statuses={name:'NOT_RUN' for name in ('NLSP_PERIODIC_PROFILE_CONVERGENCE','NLSP_PREPARED_INITIAL_STATE','NLSP_INITIAL_COMPATIBILITY_THROUGH_CUBIC_ORDER','NLSP_COMMON_INITIAL_PROJECTION','NLSP_PREPARED_SHORT_TEMPORAL_CHECK','NLSP_PREPARED_SHORT_SPATIAL_CHECK')}
    statuses.update(NLSP_PERIODIC_AXIAL_PROFILES='PARTIAL',NLSP_PREPARED_INITIAL_STATE_PILOT='PARTIAL')
    summary={'statuses':statuses,'config':config,'stop_reason':reason,'details':details,'new_ODE_integrations':0,'new_eigendecompositions':0,'old_task_status':'PARTIAL','runtime':{'numerical_wall_seconds':time.perf_counter()-started,'limit_seconds':config['budget']['numerical_wall_seconds'],'ODE_integrations':0,'eigendecompositions':0}}
    write_json(bundle/'summary.json',summary)
    return summary

def run_compute(config,bundle):
    from scripts.lib import planar_prepared_initial_state as prep
    from scripts.lib import planar_second_order_axial_response as axial
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    started=time.perf_counter()
    runner.load_runtime()
    pilot=read_json(ROOT/config['pilot_config'])
    eps=config['amplitude_over_h']
    mandatory=[ROOT/config['spectral_bundle']/'manifest.json',ROOT/pilot['audit_bundle']/'manifest.json',ROOT/pilot['audit_bundle']/'result.json',ROOT/pilot['linear_reference_bundle']/'manifest.json',ROOT/pilot['linear_reference_bundle']/'result.json']
    mandatory += [ROOT/config['spectral_bundle']/f'models/p{p}.npz' for p in config['degrees']]
    missing=[str(p) for p in mandatory if not p.is_file()]
    if missing:
        return stopped_preparation(config,bundle,started,'DATA_UNAVAILABLE: mandatory initial-preparation source',missing)
    provenance=previous.historical_provenance(config)
    previous.validate_cache(ROOT/config['spectral_bundle'])
    provenance['spectral']={'bundle':config['spectral_bundle'],'manifest_sha256':sha(ROOT/config['spectral_bundle']/'manifest.json'),'own_manifest_verified':True}
    write_json(bundle/'source_provenance.json',provenance)
    accepted=read_json(ROOT/pilot['audit_bundle']/'result.json')
    coeff=rod.RodCoefficients(**accepted['coefficients'])
    background=axial.background_from_pilot(pilot,ROOT)
    model=rod.derive_polynomials()
    models={};parts={};physical={};checks={};payload={'s':np.linspace(0,background.length,1001)}
    for p in config['degrees']:
        if time.perf_counter()-started>=config['budget']['numerical_wall_seconds']:
            write_json(bundle/'periodic_checks.json',checks)
            return stopped_preparation(config,bundle,started,'PREPARATION_BUDGET_EXHAUSTED',{'completed_degrees':list(models)})
        saved=axial.SecondOrderAxial.from_saved(coeff,p,background,ROOT/config['spectral_bundle']/f'models/p{p}.npz',model=model)
        split=prep.periodic_parts(saved)
        if split['status']!='PASS':
            return stopped_preparation(config,bundle,started,'PREPARATION_UNRESOLVED',split)
        models[p]=saved;parts[p]=split
        physical[p]={part:prep.physical_legendre_coefficients(saved,split[part]) for part in ('stat','harm')}
        check=dict(split['checks'])
        # Independent direct system solutions retain the frequency-dependent inertia.
        if p in (config['degrees'][0],config['degrees'][-1]):
            for part,operator,rhs in (('stat',saved.K,saved.f0),('harm',saved.K-saved.driving_omega**2*saved.M,saved.f2)):
                direct=np.linalg.solve(operator,rhs)
                check[part+'_direct_relative_difference']=float(np.linalg.norm(direct-split[part])/np.linalg.norm(direct))
        t=background.T1*np.array([0.,.001,.017,.1,.7])
        periodic=split['stat'][None,:]+np.cos(saved.driving_omega*t[:,None])*split['harm'][None,:]
        free=-(np.cos(t[:,None]*saved.omega)*split['free_modal_amplitudes'])@saved.vectors.T
        old=saved.evaluate(t)
        check['per_plus_free_relative_difference']=float(np.linalg.norm(periodic+free-old)/max(np.linalg.norm(old),1e-30))
        check['coordinates_retained']=saved.ndof
        check['restore_checks']=saved.restore_checks
        checks[str(p)]=check
        save_npz(bundle/'coordinates'/f'p{p}.npz',stat=split['stat'],harm=split['harm'],free_modal_amplitudes=split['free_modal_amplitudes'],legendre_stat=physical[p]['stat'],legendre_harm=physical[p]['harm'])
        for part in ('stat','harm'):
            payload[f'p{p}_{part}']=prep.LegendreProfiles(physical[p][part],background.length).evaluate(payload['s'])
    convergence={'policy':config['profile_policy'],'pairs':[]}
    for lo,hi in zip(config['degrees'][:-1],config['degrees'][1:]):
        row={'low_p':lo,'high_p':hi}
        for part in ('stat','harm'):
            row[part]=profile_comparison(physical[lo][part],physical[hi][part],background.length,background.h0,config['profile_policy'])
        row['pass']=row['stat']['pass'] and row['harm']['pass']
        # Triangle bounds are upper estimates, not the exact periodic supremum.
        row['periodic_triangle_bounds']={}
        for field in ('u','c'):
            a=next(r for r in row['stat']['rows'] if r['field']==field and r['derivative']==0)
            b=next(r for r in row['harm']['rows'] if r['field']==field and r['derivative']==0)
            row['periodic_triangle_bounds'][field]={'L2':a['absolute_L2']+b['absolute_L2'],'max_Gauss100':a['absolute_max_Gauss100']+b['absolute_max_Gauss100'],'velocity_L2':2*background.omega*b['absolute_L2'],'velocity_max_Gauss100':2*background.omega*b['absolute_max_Gauss100']}
        convergence['pairs'].append(row)
    save_npz(bundle/'profiles.npz',**payload)
    write_json(bundle/'periodic_checks.json',checks)
    write_json(bundle/'profile_convergence.json',convergence)
    trace=prep.endpoint_trace_audit(model)
    write_json(bundle/'endpoint_derivation.json',trace)
    statuses={'NLSP_PERIODIC_AXIAL_PROFILES':'PASS','NLSP_PERIODIC_PROFILE_CONVERGENCE':'PASS' if convergence['pairs'][-1]['pass'] else 'PARTIAL',
              'NLSP_PREPARED_INITIAL_STATE':'NOT_RUN','NLSP_INITIAL_COMPATIBILITY_THROUGH_CUBIC_ORDER':'NOT_RUN',
              'NLSP_COMMON_INITIAL_PROJECTION':'NOT_RUN','NLSP_PREPARED_SHORT_TEMPORAL_CHECK':'NOT_RUN',
              'NLSP_PREPARED_SHORT_SPATIAL_CHECK':'NOT_RUN','NLSP_PREPARED_INITIAL_STATE_PILOT':'PARTIAL'}
    summary={'statuses':statuses,'config':config,'background':background.as_dict(),'coefficients':coeff.values(),'periodic_checks':checks,
             'common_reference':None,'new_ODE_integrations':0,'new_eigendecompositions':0,
             'old_task_status':'PARTIAL','old_recovery_status':'PARTIAL','new_case':config['case'],
             'actual_short_runs':[],'all_eight_new_convergence_metrics':'NOT_RUN',
             'temporal_uncertainty':'NOT_RUN; old task temporal evidence is not transferred',
             'limitations':['No exact nonlinear periodic orbit','No solution of historical zero-axial IC case','No full5T1 evidence','No higher compatibility orders','LONG CLOSED; EB/RLB-KV PAUSED; angular same-clamp reference UNAVAILABLE']}
    if not convergence['pairs'][-1]['pass']:
        summary['stop_reason']='PROFILE_PREPARATION_PARTIAL'
    else:
        ref=config['degrees'][-1]
        profiles=prep.LegendreProfiles(physical[ref]['stat']+physical[ref]['harm'],background.length)
        e=np.array([0.,background.length]);uc1=profiles.evaluate(e,1);uc2=profiles.evaluate(e,2)
        b1=background.evaluate(e,1);b2=background.evaluate(e,2)
        correction=prep.quintic_theta3(uc1[:,0],b1[:,1],b1[:,0],coeff,background.length)
        state=prep.PreparedInitialState(background,profiles,correction,admitted=True)
        write_json(bundle/'full_cubic_initial_audit.json',state.endpoint_audit(coeff,config['algebraic_amplitudes'],model=model))
        summary['common_reference']={'p':ref,'source':config['spectral_bundle'],'exact_continuum_truth':False,'same_evaluator_for_every_p':True}
        eps=config['amplitude_over_h']
        # O2 traces use actual independently differentiated numerical profiles.
        r2=np.column_stack((coeff.C*uc2[:,0]+coeff.nu*coeff.C*uc1[:,1]+(coeff.C-coeff.S)*b1[:,1]*b1[:,0],
                            coeff.H*uc2[:,1]-coeff.nu*coeff.C*uc1[:,0]))
        r1=np.column_stack((coeff.S*(b2[:,0]-b1[:,1]),coeff.Bp*b2[:,1]+coeff.S*b1[:,0]))
        c1=correction.evaluate(e,1);c2=correction.evaluate(e,2)
        r3=np.column_stack((-coeff.S*c1+(coeff.C-coeff.S)*uc1[:,0]*b1[:,1],
                            coeff.Bp*c2-(coeff.C-coeff.S)*uc1[:,0]*b1[:,0]))
        scales2=np.array([coeff.C*background.h0/background.length**2,coeff.C]);scales13=np.array([coeff.S*background.h0/background.length**2,coeff.S*background.h0/background.length])
        compatibility={'epsilon1_bending':r1.tolist(),'epsilon2_axial_contraction':r2.tolist(),'epsilon3_bending_before_correction':np.column_stack(((coeff.C-coeff.S)*uc1[:,0]*b1[:,1],-(coeff.C-coeff.S)*uc1[:,0]*b1[:,0])).tolist(),
                       'epsilon3_bending_after_correction':r3.tolist(),'scales_O2':scales2.tolist(),'scales_O1_O3':scales13.tolist(),
                       'fixed_scaled_max':float(max(np.max(abs(r2)/scales2),np.max(abs(r1)/scales13),np.max(abs(r3)/scales13))),
                       'numeric_derivatives_not_replaced_by_BVP':True,'formal_evidence':'endpoint_derivation.json',
                       'finite_amplitudes':[numeric_endpoint_audit(state,a,coeff) for a in config['algebraic_amplitudes']]}
        compatible=compatibility['fixed_scaled_max']<=config['profile_policy']['endpoint_tolerance']
        correction_info=correction.as_dict()
        vals=state.evaluate(payload['s'],eps)
        correction_size=float(np.max(abs(eps**3*correction.evaluate(payload['s'])))/np.max(abs(eps*background.evaluate(payload['s'])[:,1])))
        correction_info['physical_correction_to_main_theta_max_ratio']=correction_size
        write_json(bundle/'theta3.json',correction_info)
        write_json(bundle/'compatibility.json',compatibility)
        save_npz(bundle/'common_initial_state.npz',s=payload['s'],physical_fields=vals,U_C_legendre=profiles.coefficients,theta3_eta_power=correction.eta_power_coefficients,epsilon_a=np.array(eps))
        statuses['NLSP_PREPARED_INITIAL_STATE']='PASS' if compatible else 'PARTIAL'
        statuses['NLSP_INITIAL_COMPATIBILITY_THROUGH_CUBIC_ORDER']='PASS' if compatible else 'PARTIAL'
        summary['compatibility_label']='ENDPOINT_ACCELERATION_COMPATIBLE_THROUGH_CUBIC_ORDER' if compatible else 'UNRESOLVED'
        summary['compatibility']=compatibility;summary['theta3']=correction_info
        if compatible:
            projection={};states={}
            for p in (32,48,64):
                disc=models[p].disc
                q,result=project_initial(state,disc,eps,config['profile_policy'])
                runner.safety_check(disc,q,pilot['safety'])
                projection[str(p)]=result;states[p]=q
                save_npz(bundle/'initial_projection'/f'p{p}.npz',q=q,velocity=np.zeros(disc.ndof))
            write_json(bundle/'initial_projection.json',projection)
            primary=config['nonlinear_policy']['primary_pair'];alternate=config['nonlinear_policy']['allowed_pre_run_replacement']
            chosen=primary if all(projection[str(p)]['pass'] for p in primary) else alternate if all(projection[str(p)]['pass'] for p in alternate) else None
            decision={'primary_pair':primary,'allowed_replacement':alternate,'selected_pair':chosen,'made_before_trajectories':True,'no_independent_per_p_preparation':True}
            write_json(bundle/'pre_run_decision.json',decision)
            summary['projection']=projection;summary['pre_run_decision']=decision
            statuses['NLSP_COMMON_INITIAL_PROJECTION']='PASS' if chosen else 'PARTIAL'
            if chosen is None:
                summary['stop_reason']='COMMON_INITIAL_PROJECTION_PARTIAL: neither permitted pair resolves common initial profiles and endpoint jets'
            else:
                # This bounded gate is reached only after admission; the existing
                # autonomous four-field runner remains the sole time integrator.
                run_short_controls(state,chosen,models,pilot,config,bundle,summary,started)
        else:
            summary['stop_reason']='INITIAL_COMPATIBILITY_PARTIAL'
    recovery_path=ROOT/config['historical_bundles']['recovery']
    old_required=[recovery_path/'summary.json']+[recovery_path/'controls'/name/'trajectory.npz' for name in ('new_p32','p48_strict_short')]
    old_available=all(path.is_file() for path in old_required) and provenance.get('recovery',{}).get('status')=='PASS'
    if not old_available:
        summary['old_short_comparison']={'status':'DATA_UNAVAILABLE'}
        summary['old_short_qualification']='Unavailable old raw controls; completed preparation retained; no old rerun'
    else:
        recovery=read_json(ROOT/config['historical_bundles']['recovery']/'summary.json')
        old_controls={}
        for name in ('new_p32','p48_strict_short'):
            folder=ROOT/config['historical_bundles']['recovery']/'controls'/name
            with np.load(folder/'trajectory.npz',allow_pickle=False) as saved:
                old_controls[name]={key:saved[key].copy() for key in ('q','velocity','time')}
        if not np.array_equal(old_controls['new_p32']['time'],old_controls['p48_strict_short']['time']):
            raise ValueError('Old short controls do not share actual timestamps')
        summary['old_short_comparison']=runner.compare_histories(models[32].disc,old_controls['new_p32'],models[48].disc,old_controls['p48_strict_short'],pilot)
        summary['old_short_actual_times']={'samples':len(old_controls['new_p32']['time']),'start':float(old_controls['new_p32']['time'][0]),'end':float(old_controls['new_p32']['time'][-1]),'periods':float(old_controls['new_p32']['time'][-1]/background.T1),'exact_common_grid':True}
        for key,row in summary['old_short_comparison']['fields'].items():
            original=recovery['short_p32_p48_diagnostic']['fields'][key]
            if any(abs(row[metric]-original[metric])>2e-11*max(abs(original[metric]),1e-30) for metric in ('relative_L2','relative_max')):
                raise ArithmeticError('Old short norms changed: '+key)
            field=key.split('_')[-1];speed=key.startswith('velocity')
            fixed_scale=eps*background.h0/(background.length if field in ('theta','c') else 1.)
            if speed:
                fixed_scale*=np.sqrt(coeff.C/coeff.jp)
            row['fixed_common_physical_scale']=fixed_scale
            row['absolute_L2_over_fixed_common_scale']=row['max_time_L2_difference']/(fixed_scale*np.sqrt(background.length))
            row['absolute_max_over_fixed_common_scale']=row['max_space_time_difference']/fixed_scale
        summary['old_short_qualification']='p32 tight versus p48 allowed_extra on identical actual0..0.1T1 saved grid; historical PARTIAL, no rerun'
    summary['old_new_improvement']='NOT_ESTABLISHED: no new nonlinear trajectories' if summary['new_ODE_integrations']==0 else 'Compare absolute differences and fixed common scales; different initial tasks' 
    write_json(bundle/'old_short_comparison.json',{'comparison':summary['old_short_comparison'],'qualification':summary['old_short_qualification'],'new_metrics':'NOT_RUN'})
    summary['runtime']={'numerical_wall_seconds':time.perf_counter()-started,'limit_seconds':config['budget']['numerical_wall_seconds'],'ODE_integrations':summary['new_ODE_integrations'],'BVP_direct_validation_solves':4,'eigendecompositions':0,'saved_spectral_restores':len(models),'old_evaluator_spot_calls':len(models),'new_full_analytic_histories':0}
    write_json(bundle/'summary.json',summary)
    return summary



def run_short_controls(state, chosen, models, pilot, config, bundle, summary, started):
    """Reuse the old autonomous runner only after both initial projections pass.

    The current recorded case fails admission; this path is never entered in
    that audit. No test is permitted to advance an extra nonlinear trajectory.
    """
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    if not all(summary['projection'][str(p)]['pass'] for p in chosen):
        raise ValueError('Short integration requires both projection gates')
    runner.load_runtime()
    eps=config['amplitude_over_h']; amplitude=eps*state.background.h0
    shape=lambda x: state.evaluate(x,eps)/amplitude
    horizon=config['nonlinear_policy']['short_periods']*state.background.T1
    omega_bound=0.
    for p in chosen:
        disc=models[p].disc
        diagonal=np.diag(disc.M0)
        lower=float(np.min(diagonal-(np.sum(abs(disc.M0),axis=1)-abs(diagonal))))
        if lower<=0:
            raise ArithmeticError('Nonpositive resting-mass Gershgorin bound')
        omega_bound=max(omega_bound,float(np.sqrt(np.max(np.sum(abs(disc.K),axis=1))/(lower*pilot['safety']['min_relative_mass_eigenvalue']))))
    count=int(np.ceil(horizon*omega_bound/(2*np.pi)*16))+1
    historical=ROOT/config['historical_bundles']['recovery']/'controls/new_p32/trajectory.npz'
    with np.load(historical,allow_pickle=False) as saved:
        old_times=saved['time'].copy()
    times=np.unique(np.r_[np.linspace(0,horizon,count),old_times[old_times<=horizon]])
    write_json(bundle/'short_sampling.json',{'horizon':horizon,'samples':len(times),'omega_upper_bound':omega_bound,'samples_per_upper_bound_period':16,'bound':'Gershgorin K norm / positive M0 lower bound and accepted variable-mass lower bound','not_internal_time_error_control':True})
    histories={};deadline=started+config['budget']['numerical_wall_seconds']
    cases=[(chosen[0],'tight'),(chosen[1],'tight'),(chosen[1],'allowed_extra')]
    for p,level in cases:
        if time.perf_counter()>=deadline:
            break
        disc=models[p].disc
        history,stats=runner.integrate_case(disc,shape,state.background.as_dict(),pilot,eps,level,times,deadline)
        summary['new_ODE_integrations']+=1
        summary['actual_short_runs'].append(stats)
        actual=times[:len(history)];q,v=history[:,:disc.ndof],history[:,disc.ndof:]
        energy=np.array([disc.energy(a,b) for a,b in zip(q,v)])
        bounds=[mass_safety_bounds(disc,a) for a in q]
        stats['relative_energy_drift_max']=float(np.max(abs((energy-energy[0])/energy[0])))
        stats['mass_lower_bound_min']=min(a['relative_mass_eigenvalue_lower_bound'] for a in bounds)
        stats['mass_condition_bound_max']=max(a['relative_mass_condition_upper_bound'] for a in bounds)
        name=f'p{p}_{level}'
        save_npz(bundle/'short_controls'/name/'trajectory.npz',time=actual,q=q,velocity=v,energy=energy)
        write_json(bundle/'short_controls'/name/'case.json',stats)
        histories[name]={'q':q,'velocity':v,'time':actual}
        if stats['status']!='PASS':
            break
    def comparison(a,b):
        if a not in histories or b not in histories:
            return {'status':'NOT_RUN'}
        aa,bb=histories[a],histories[b]
        if not np.array_equal(aa['time'],bb['time']) or aa['time'][-1]!=horizon:
            return {'status':'PARTIAL','reason':'actual prefix or unequal timestamps'}
        pa,pb=int(a.split('_')[0][1:]),int(b.split('_')[0][1:])
        return runner.compare_histories(models[pa].disc,aa,models[pb].disc,bb,pilot)
    spatial=comparison(f'p{chosen[0]}_tight',f'p{chosen[1]}_tight')
    temporal=comparison(f'p{chosen[1]}_tight',f'p{chosen[1]}_allowed_extra')
    summary['new_spatial_comparison']=spatial;summary['new_temporal_comparison']=temporal
    summary['all_eight_new_convergence_metrics']=spatial
    summary['statuses']['NLSP_PREPARED_SHORT_SPATIAL_CHECK']=spatial['status']
    summary['statuses']['NLSP_PREPARED_SHORT_TEMPORAL_CHECK']=temporal['status']
    energy_pass=len(summary['actual_short_runs'])==3 and all(r['relative_energy_drift_max']<=pilot['gates']['energy_relative_drift'] for r in summary['actual_short_runs'])
    summary['statuses']['NLSP_PREPARED_INITIAL_STATE_PILOT']='PASS' if spatial['status']==temporal['status']=='PASS' and energy_pass else 'PARTIAL'
    summary['temporal_uncertainty']=temporal
    summary['stop_reason']='Bounded0.1T1 short controls complete or actual prefix saved; no automatic extension'


def manifest_for(bundle,item):
    return {'identity':item,'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'git_status_before_outputs':subprocess.check_output(['git','status','--short'],cwd=ROOT,text=True),
            'command':sys.argv,'artifact_hashes':{str(p.relative_to(bundle)):sha(p) for p in bundle.rglob('*') if p.is_file() and p.name!='manifest.json' and not p.name.endswith('.tmp')}}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true')
    mode.add_argument('--report-only',type=Path)
    mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG)
    parser.add_argument('--output-dir',type=Path,default=OUTPUT)
    args=parser.parse_args(argv)
    zeros={'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0}
    if args.report_only or args.plot_only:
        bundle=args.report_only or args.plot_only
        summary=validate_cache(bundle)
        if args.plot_only:
            plot_only(bundle)
        print(json.dumps({'bundle':str(bundle),'statuses':summary['statuses'],'this_run_counters':zeros},indent=2))
        return summary
    key,item=identity(args.config);bundle=args.output_dir/key
    if (bundle/'manifest.json').exists():
        summary=validate_cache(bundle,item)
        print(json.dumps({'bundle':str(bundle),'cache_hit':True,'statuses':summary['statuses'],'this_run_counters':zeros},indent=2))
        return summary
    bundle.mkdir(parents=True,exist_ok=True)
    summary=run_compute(read_json(args.config),bundle)
    write_json(bundle/'manifest.json',manifest_for(bundle,item))
    plot_only(bundle)
    write_json(bundle/'manifest.json',manifest_for(bundle,item))
    print(json.dumps({'bundle':str(bundle),'cache_hit':False,'statuses':summary['statuses'],'runtime':summary['runtime'],'stop_reason':summary.get('stop_reason')},indent=2))
    return summary


if __name__ == '__main__':
    main()
