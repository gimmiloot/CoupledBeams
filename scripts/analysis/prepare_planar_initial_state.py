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
    cached=validate_cache(bundle)
    if cached.get('schema')=='nlsp-prepared-one-T1-v1':
        return one_T1_plot_only(bundle,cached)
    if cached.get('schema')=='nlsp-prepared-feasibility-v1':
        return feasibility_plot_only(bundle,cached)
    if not (bundle/'profiles.npz').exists():
        return {'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0,'symbolic_derivations':0}
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
        return {'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0,'symbolic_derivations':0}
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
    return {'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0,'symbolic_derivations':0}


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



FEASIBILITY_CONFIG = ROOT/'data/input/planar_prepared_feasibility.json'
FEASIBILITY_OUTPUT = ROOT/'results/planar_prepared_feasibility'
FEASIBILITY_VERSION = 'frozen-prepared-initial-only-endpoint-L2-feasibility-v1'


def feasibility_identity(config_path):
    key,item=identity(config_path)
    config=item['config']
    item['version']=FEASIBILITY_VERSION
    for name in ('prepared_bundle','precision_evidence'):
        path=ROOT/config[name]
        manifest=read_json(path/'manifest.json')
        item['sources'][name]={'bundle':config[name],'manifest_sha256':sha(path/'manifest.json'),
                              'artifact_hashes':manifest['artifact_hashes']}
    item['dependencies']['mpmath']=importlib.metadata.version('mpmath')
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item


def authorize_feasibility(strict_table,basic,evidence):
    """Explicit authorization; never changes the physical state's admitted flag."""
    if not all(basic.values()) or not evidence:
        raise ArithmeticError('BLOCKED_BY_UNEXPLAINED_INCONSISTENCY: basic correctness or independent precision evidence missing')
    return 'STRICT_ADMITTED' if all(row['pass'] for row in strict_table) else 'EXPLORATORY_NOT_CERTIFIED'


def initial_coordinate_metrics(state,disc,q,epsilon,policy):
    """Historical norms, floors, endpoint scales, without reprojecting q."""
    xi,ww=leg.leggauss(max(100,2*disc.p+1));x=(xi+1)*disc.length/2;ww*=disc.length/2
    fixed=np.array([epsilon*state.background.h0]*2+[epsilon*state.background.h0/disc.length]*2)
    rows=[]
    for d in policy['spatial_derivatives']:
        truth=state.evaluate(x,epsilon,d,require_admitted=False)
        actual=disc.reconstruct(q,x,d)
        jets=disc.reconstruct(q,[0.,disc.length],d)-state.evaluate([0.,disc.length],epsilon,d,require_admitted=False)
        for f,field in enumerate(('u','w','theta','c')):
            dif=actual[:,f]-truth[:,f];scale=fixed[f]/disc.length**d
            norm=float(np.sqrt(ww@dif**2));peak=float(np.max(abs(dif)))
            own=float(np.sqrt(ww@truth[:,f]**2));ownmax=float(np.max(abs(truth[:,f])))
            row={'field':field,'derivative':d,'absolute_L2':norm,'absolute_max':peak,
                 'relative_L2':norm/max(own,1e-10*scale*np.sqrt(disc.length)),
                 'relative_max':peak/max(ownmax,1e-10*scale),
                 'endpoint_fixed_scaled_error':float(np.max(abs(jets[:,f]))/scale)}
            row['pass']=row['relative_L2']<=policy['relative_tolerance'] and row['relative_max']<=policy['relative_tolerance'] and row['endpoint_fixed_scaled_error']<=policy['endpoint_tolerance']
            rows.append(row)
    return {'rows':rows,'pass':all(r['pass'] for r in rows),'unchanged_historical_norms':True,
            'sampled_maximum':True,'comparison_quadrature':len(x),'essential_endpoint_values':disc.reconstruct(q,[0.,disc.length]).tolist()}


def initial_formal_compatibility(state,disc,projection,epsilon):
    """Unchanged through-cubic traces; split a linear initial-only projection."""
    from scripts.lib import planar_prepared_initial_state as prep
    eta=np.polynomial.Polynomial(state.correction.coefficients)
    power=eta(np.polynomial.Polynomial([.5,.5])).coef
    lc=leg.poly2leg(power)
    exact=prep.project_saved_legendre(lc,disc.p,disc.length,prep.UNCONSTRAINED_PROJECTION,70)
    r3=np.zeros(disc.ndof);r3[disc.slices['theta']]=exact['raw']
    q3=disc.from_raw_coefficients(r3)
    q=projection['q'];q1=np.zeros(disc.ndof);q2=q1.copy()
    for field in ('u','c'):q2[disc.slices[field]]=q[disc.slices[field]]/epsilon**2
    q1[disc.slices['w']]=q[disc.slices['w']]/epsilon
    q1[disc.slices['theta']]=(q[disc.slices['theta']]-epsilon**3*q3[disc.slices['theta']])/epsilon
    d1,d2,d3=[disc.reconstruct(a,[0.,disc.length],1) for a in (q1,q2,q3)]
    dd1,dd2,dd3=[disc.reconstruct(a,[0.,disc.length],2) for a in (q1,q2,q3)]
    c=disc.coefficients
    e1=np.column_stack((np.zeros(2),c.S*(dd1[:,1]-d1[:,2]),c.Bp*dd1[:,2]+c.S*d1[:,1],np.zeros(2)))
    e2=np.column_stack((c.C*dd2[:,0]+c.nu*c.C*d2[:,3]+(c.C-c.S)*d1[:,2]*d1[:,1],np.zeros(2),np.zeros(2),c.H*dd2[:,3]-c.nu*c.C*d2[:,0]))
    e3=np.column_stack((np.zeros(2),-c.S*d3[:,2]+(c.C-c.S)*d2[:,0]*d1[:,2],c.Bp*dd3[:,2]-(c.C-c.S)*d2[:,0]*d1[:,1],np.zeros(2)))
    scales=np.array([c.C*state.background.h0/disc.length**2,c.S*state.background.h0/disc.length**2,c.S*state.background.h0/disc.length,c.C])
    peak=float(max(np.max(abs(e)/scales) for e in (e1,e2,e3)))
    return {'through_cubic_coefficients':{'epsilon1':e1.tolist(),'epsilon2':e2.tolist(),'epsilon3':e3.tolist()},
            'compatibility_fixed_scaled_max':peak,'pass':peak<=1e-6,'Theta3_not_refitted':True}


def strong_weak_metrics(disc,q,velocity,label):
    acceleration=disc.acceleration(q,velocity)
    gradient=disc.potential(q)['gradient'];mass=disc.mass_matrix(q)@acceleration
    inertia=disc.inertial_terms(q,velocity)
    action=mass+inertia+gradient;strong=disc.weak_residual(q,velocity,acceleration)
    delta=strong-action
    local=disc._local_potential(q)[1]
    work=sum(np.linalg.norm(matrix.T@(disc.weights*value)) for matrix,value in zip(disc._potential_matrices,local))+np.linalg.norm(mass)+np.linalg.norm(inertia)
    correction=disc._solve_mass(q,delta)
    fields=disc.reconstruct(correction)
    absolute=float(np.max(abs(delta)));relative=float(np.linalg.norm(delta)/work)
    return {'state':label,'absolute_max':absolute,'difference_L2':float(np.linalg.norm(delta)),
            'uncancelled_work_scale':float(work),'relative_residual':relative,'absolute_gate':2e-12,'relative_gate':2e-12,
            'absolute_pass':absolute<=2e-12,'relative_pass':relative<=2e-12,'pass':absolute<=2e-12 and relative<=2e-12,
            'difference_vector':delta.tolist(),'point_acceleration_difference_L2':np.sqrt(disc.weights@fields**2).tolist(),
            'point_acceleration_difference_max':np.max(abs(fields),axis=0).tolist(),
            'qualification':'Point-state assembly discrepancy, not a trajectory error bound; RHS unmodified'}


def run_feasibility(config,bundle):
    """Three bounded runs of one immutable state, strict and exploratory separate."""
    import shutil
    from types import SimpleNamespace
    from scripts.lib import planar_prepared_initial_state as prep
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    started=time.perf_counter();previous_charge=config['budget']['previous_precision_seconds']
    deadline=started+config['budget']['numerical_wall_seconds']-previous_charge
    local_started=started
    runner.load_runtime();pilot=read_json(ROOT/config['pilot_config'])
    source=ROOT/config['prepared_bundle'];old=validate_cache(source)
    for path in (config['spectral_bundle'],*config['historical_bundles'].values()):validate_cache(ROOT/path)
    precision=ROOT/config['precision_evidence'];pm=read_json(precision/'manifest.json')
    for name,digest in pm['artifact_hashes'].items():
        if sha(precision/name)!=digest:raise ValueError('Precision artifact hash mismatch: '+name)
        target=bundle/'precision_evidence'/name;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(precision/name,target)
    evidence=read_json(precision/'nlsp_strong_weak_precision_20261008.json')
    if evidence['input_hashes']['source_manifest']!=sha(source/'manifest.json'):raise ValueError('Precision common-source mismatch')
    source_state,coeff,provenance=prep.load_frozen_prepared_state(source)
    audit=read_json(ROOT/pilot['audit_bundle']/'result.json')
    if evidence['input_hashes']['audit_result']!=sha(ROOT/pilot['audit_bundle']/'result.json'):raise ValueError('Precision frozen action mismatch')
    pol=audit['polynomials']
    model=SimpleNamespace(T4=rod.Polynomial.deserialize(pol['T4']),V4=rod.Polynomial.deserialize(pol['V4']),
        residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol['residuals_A']),symbols={a:rod.Polynomial.symbol(a) for a in rod.SYMBOL_ORDER})
    precise_rows=[a for a in evidence['rows'] if a.get('stage')=='refined_gauss' and a.get('precision') in (45,70)]
    independent_ok=len(precise_rows)==8 and all(float(a['relative_on_original_float64_scale'])<=2e-12 for a in precise_rows)
    eps=config['amplitude_over_h'];discs={};qs={};projection={};strict=[];basic={};identities={};original={}
    for p in config['degrees']:
        d=PlanarGalerkin(coeff,p,model=model);discs[p]=d
        copies=[prep.stable_initial_projection(source_state,d,eps,config['projection_policy'],digits) for digits in config['projection_dps']]
        rep=copies[-1];q=rep['q'];qs[p]=q
        before=initial_coordinate_metrics(source_state,d,d.project(source_state.evaluate(d.x,eps,require_admitted=False)),eps,config['profile_policy'])
        current=initial_coordinate_metrics(source_state,d,q,eps,config['profile_policy'])
        compatibility=initial_formal_compatibility(source_state,d,rep,eps)
        current.update(compatibility);current['pass']=all(r['pass'] for r in current['rows']) and compatibility['pass']
        current['policy']=config['projection_policy'];current['dps']=config['projection_dps'];current['same_coefficients_after_float64_rounding']=bool(np.array_equal(copies[0]['raw'],rep['raw']))
        current['roundtrip_relative']=rep['raw_roundtrip_relative'];current['source_jets']=rep['source_jets'].tolist();current['source_jets_decimal']=rep['source_jets_decimal']
        current['mp_endpoint_error']=rep['mp_endpoint_error'].tolist();current['float64_endpoint_error']=rep['float64_endpoint_error'].tolist()
        projection[str(p)]=current;original[str(p)]=before
        save_npz(bundle/'initial_projection'/f'p{p}.npz',q=q,velocity=np.zeros(d.ndof),raw=rep['raw'],raw_roundtrip=rep['raw_roundtrip'],source_jets=rep['source_jets'])
        write_json(bundle/'initial_projection'/f'p{p}.json',current)
        write_json(bundle/'initial_projection'/f'p{p}_decimal.json',{'raw_decimal':rep['raw_decimal'],'source_jets_decimal':rep['source_jets_decimal']})
        synthetic=d.project(np.column_stack((d.x*(1-d.x)*1e-8,d.x*(1-d.x)*1e-7,d.x*(1-d.x)*1e-8,d.x*(1-d.x)*1e-9)))
        checks=[strong_weak_metrics(d,q,np.zeros(d.ndof),'initial_zero_velocity'),strong_weak_metrics(d,q,synthetic,'previous_synthetic_velocity')]
        identities[str(p)]=checks
        strict.append({'p':p,'check':'initial_projection_and_formal_jets','pass':current['pass'],'tolerance':1e-6})
        strict.extend({'p':p,'check':a['state']+'_strong_weak','pass':a['pass'],'absolute':a['absolute_max'],'relative':a['relative_residual'],'tolerance':2e-12} for a in checks)
        runner.safety_check(d,q,pilot['safety']);diag=mass_safety_bounds(d,q)
        basic[str(p)+'_finite_q_RHS']=bool(np.all(np.isfinite(q)) and np.all(np.isfinite(d.rhs(0,np.r_[q,np.zeros(d.ndof)]))))
        basic[str(p)+'_positive_mass']=diag['mass_positive'] and diag['relative_mass_eigenvalue_lower_bound']>=pilot['safety']['min_relative_mass_eigenvalue']
        basic[str(p)+'_essential_BC']=bool(np.max(abs(d.reconstruct(q,[0.,d.length])))==0)
        basic[str(p)+'_zero_velocities']=True
    local_elapsed=time.perf_counter()-local_started+previous_charge
    if local_elapsed>config['budget']['local_precision_seconds']:raise TimeoutError('Local precision budget exhausted before ODE')
    execution=authorize_feasibility(strict,basic,independent_ok)
    statuses={'NLSP_PROJECTION_ARITHMETIC_AUDIT':'COMPLETED','NLSP_NUMERICAL_REPRESENTATION_FIX':'PASS' if all(a['pass'] for a in projection.values()) else 'PARTIAL',
        'NLSP_STRICT_INITIAL_VERIFICATION':'PASS' if execution=='STRICT_ADMITTED' else 'PARTIAL',
        'NLSP_PREPARED_FEASIBILITY_RUN':'NOT_RUN','NLSP_PREPARED_SHORT_SPATIAL_CHECK':'NOT_RUN','NLSP_PREPARED_SHORT_TEMPORAL_CHECK':'NOT_RUN'}
    summary={'schema':config['schema'],'config':config,'background':source_state.background.as_dict(),'coefficients':coeff.values(),'statuses':statuses,
        'historical_statuses':old['statuses'],'source_provenance':provenance,'projection':projection,'original_projection':original,
        'strong_weak':identities,'strict_table':strict,'execution_mode':execution,'execution_modes':{},'actual_short_runs':[],
        'new_ODE_integrations':0,'new_eigendecompositions':0,'new_BVP_solves':0,'state_admitted_flag':source_state.admitted,
        'numerical_policy_chosen_before_ODE':True,'independent_precision_evidence':independent_ok,
        'qualification':'Short finite-dimensional feasibility, not continuous-PDE convergence, periodic-orbit or out-of-plane certification'}
    write_json(bundle/'pre_run_decision.json',{'execution_mode':execution,'state_admitted_flag':source_state.admitted,'strict_table':strict,'basic_checks':basic,'independent_precision_evidence':independent_ok,'policy':config['projection_policy'],'new_dynamic_constraints':False})
    write_json(bundle/'projection_comparison.json',{'before':original,'after':projection})
    save_npz(bundle/'common_initial_state.npz',**{k:v for k,v in np.load(source/'common_initial_state.npz',allow_pickle=False).items()})
    horizon=config['short_periods']*source_state.background.T1;omega_bound=0.
    for d in discs.values():
        diagonal=np.diag(d.M0);lower=float(np.min(diagonal-(np.sum(abs(d.M0),axis=1)-abs(diagonal))))
        if lower<=0:raise ArithmeticError('Nonpositive resting-mass bound')
        omega_bound=max(omega_bound,float(np.sqrt(np.max(np.sum(abs(d.K),axis=1))/(lower*pilot['safety']['min_relative_mass_eigenvalue']))))
    count=int(np.ceil(horizon*omega_bound/(2*np.pi)*16))+1
    with np.load(ROOT/config['historical_bundles']['recovery']/'controls/new_p32/trajectory.npz') as saved:old_times=saved['time'].copy()
    times=np.unique(np.r_[np.linspace(0,horizon,count),old_times[old_times<=horizon]])
    summary['sampling']={'horizon':horizon,'samples':len(times),'omega_upper_bound':omega_bound,'samples_per_bound_period':16,'old_times_preserved':True,'output_grid_not_time_accuracy_control':True}
    summary['runtime']={'local_precision_seconds':local_elapsed,'limit_seconds':config['budget']['numerical_wall_seconds']}
    write_json(bundle/'summary.json',summary)
    print(json.dumps({'stage':'PRE_RUN','bundle':str(bundle),'execution':execution,'strict_table':strict,'samples':len(times),'local_precision_seconds':local_elapsed}),flush=True)
    histories={}
    for p,level in config['cases']:
        if time.perf_counter()>=deadline:break
        d=discs[p];name=f'p{p}_{level}'
        shape=lambda x:source_state.evaluate(x,eps,require_admitted=False)/(eps*source_state.background.h0)
        print(json.dumps({'stage':'START_ODE','case':name,'target_time':horizon,'execution_mode':execution}),flush=True)
        summary['new_ODE_integrations']+=1
        history,stats=runner.integrate_case(d,shape,source_state.background.as_dict(),pilot,eps,level,times,deadline,initial_coordinates=qs[p])
        actual=times[:len(history)];q,v=history[:,:d.ndof],history[:,d.ndof:]
        energy=np.array([d.energy(a,b) for a,b in zip(q,v)])
        bounds=[mass_safety_bounds(d,a) for a in q]
        stats.update(execution_mode=execution,projection_policy=config['projection_policy'],initial_energy=float(energy[0]),
            relative_energy_drift_max=float(np.max(abs((energy-energy[0])/energy[0]))),
            mass_lower_bound_min=min(a['relative_mass_eigenvalue_lower_bound'] for a in bounds),
            mass_condition_bound_max=max(a['relative_mass_condition_upper_bound'] for a in bounds))
        stats['safety_extrema']={key:(min(a[key] for a in bounds) if key=='min_one_plus_c' else max(a[key] for a in bounds)) for key in pilot['safety'] if key!='min_relative_mass_eigenvalue'}
        summary['actual_short_runs'].append(stats);summary['execution_modes'][name]=execution
        save_npz(bundle/'short_controls'/name/'trajectory.npz',time=actual,q=q,velocity=v,energy=energy)
        save_npz(bundle/'short_controls'/name/'internal_steps.npz',time_step=np.asarray(stats['internal_time_steps']))
        write_json(bundle/'short_controls'/name/'case.json',stats)
        histories[name]={'q':q,'velocity':v,'time':actual}
        write_json(bundle/'summary.json',summary)
        print(json.dumps({'stage':'END_ODE','case':name,'stats':stats}),flush=True)
        if stats['status']!='PASS':break
    def comparison(a,b):
        if a not in histories or b not in histories:return {'status':'NOT_RUN'}
        aa,bb=histories[a],histories[b]
        if not np.array_equal(aa['time'],bb['time']) or aa['time'][-1]!=horizon:return {'status':'PARTIAL','reason':'actual prefix; no repeated final snapshot'}
        da,db=discs[int(a.split('_')[0][1:])],discs[int(b.split('_')[0][1:])]
        result=runner.compare_histories(da,aa,db,bb,pilot)
        fixed=np.array([eps*source_state.background.h0]*2+[eps*source_state.background.h0/source_state.length]*2)
        for key,row in result.get('fields',{}).items():
            idx=('u','w','theta','c').index(key.split('_')[-1]);scale=fixed[idx]
            if key.startswith('velocity'):scale*=source_state.background.omega
            row['fixed_physical_scale']=float(scale);row['fixed_scaled_L2']=row['max_time_L2_difference']/(scale*np.sqrt(source_state.length));row['fixed_scaled_max']=row['max_space_time_difference']/scale
        result['sampling_qualification']='Shared actual dense-output samples; maxima sampled, no phase alignment. Output bound 16 points/period, not independent accuracy proof.'
        return result
    spatial=comparison('p48_tight','p64_tight');temporal=comparison('p64_tight','p64_allowed_extra')
    plot_data={'T1':np.asarray(source_state.background.T1)}
    for name,hist in histories.items():
        d=discs[int(name.split('_')[0][1:])]
        plot_data[name+'_time']=hist['time']
        plot_data[name+'_quarter_fields']=d.reconstruct_series(hist['q'],[source_state.length/4])[:,0,:]
        with np.load(bundle/'short_controls'/name/'trajectory.npz') as saved:plot_data[name+'_energy']=saved['energy'].copy()
    save_npz(bundle/'plot_data.npz',**plot_data)
    summary['new_spatial_comparison']=spatial;summary['new_temporal_comparison']=temporal
    statuses['NLSP_PREPARED_SHORT_SPATIAL_CHECK']=spatial['status'];statuses['NLSP_PREPARED_SHORT_TEMPORAL_CHECK']=temporal['status']
    complete=len(summary['actual_short_runs'])==3 and all(r['status']=='PASS' for r in summary['actual_short_runs'])
    statuses['NLSP_PREPARED_FEASIBILITY_RUN']=('COMPLETED_'+execution) if complete else 'PARTIAL'
    summary['runtime'].update(numerical_wall_seconds=time.perf_counter()-started+previous_charge,integration_seconds=sum(r['integration_seconds'] for r in summary['actual_short_runs']),ODE_integrations=summary['new_ODE_integrations'],eigendecompositions=0,BVP_solves=0)
    summary['stop_reason']='Bounded0.1T1 controls completed or actual prefix saved; strict thresholds unchanged; no full-period extension'
    write_json(bundle/'summary.json',summary)
    return summary



def feasibility_plot_only(bundle,summary):
    """Only saved physical arrays, zero preparation/ODE/eigen evaluations."""
    import matplotlib.pyplot as plt
    bundle=Path(bundle);figs=bundle/'figures';figs.mkdir(exist_ok=True)
    def finish(fig,name):
        for ext in ('pdf','png'):fig.savefig(figs/(name+'.'+ext),dpi=220,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.4),layout='constrained')
    names=[f'{field} d{d}' for d in (0,1,2) for field in ('u','w','theta','c')]
    for ax,p in zip(axes,(48,64)):
        for label,key in (('original float64','original_projection'),('initial endpoint L2','projection')):
            rows=summary[key][str(p)]['rows'];ax.semilogy(range(12),[max(r['relative_max'],r['relative_L2'],r['endpoint_fixed_scaled_error'],1e-18) for r in rows],'.-',label=label)
        ax.axhline(1e-6,color='k',ls=':',lw=.8);ax.set_xticks(range(12),names,rotation=60,fontsize=7)
        ax.set(title=f'p={p}',ylabel='initial representation errors');ax.grid(alpha=.2);ax.legend(fontsize=7,frameon=False)
    finish(fig,'initial_projection_before_after')
    if not (bundle/'plot_data.npz').exists():return {'ODE_integrations':0,'eigendecompositions':0}
    with np.load(bundle/'plot_data.npz',allow_pickle=False) as saved:data={k:saved[k] for k in saved.files}
    fig,axes=plt.subplots(2,2,figsize=(9,5.6),layout='constrained')
    for f,(field,ax) in enumerate(zip(('u','w','theta','c'),axes.flat)):
        for name,style in (('p48_tight','-'),('p64_tight','--')):
            if name+'_time' in data:ax.plot(data[name+'_time']/data['T1'],data[name+'_quarter_fields'][:,f],style,label=name.replace('_',' '))
        ax.set(xlabel='t/T1',ylabel=field+' at s=L/4');ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(fig,'prepared_short_trajectories')
    fig,ax=plt.subplots(figsize=(7,3.4),layout='constrained')
    for name in ('p48_tight','p64_tight','p64_allowed_extra'):
        if name+'_time' not in data:continue
        e=data[name+'_energy'];ax.plot(data[name+'_time']/data['T1'],(e-e[0])/e[0],label=name.replace('_',' '))
    ax.set(xlabel='t/T1',ylabel='relative semidiscrete energy drift');ax.grid(alpha=.2);ax.legend(frameon=False,fontsize=8)
    finish(fig,'prepared_short_energy')
    return {'ODE_integrations':0,'eigendecompositions':0,'BVP_solves':0}



ONE_T1_CONFIG=ROOT/'data/input/planar_prepared_one_T1.json'
ONE_T1_OUTPUT=ROOT/'results/planar_prepared_one_T1'
ONE_T1_VERSION='saved-q0-one-linear-period-memmap-v1'


def one_T1_identity(config_path):
    config=read_json(config_path);source=ROOT/config['source_bundle']
    paths=('scripts/analysis/prepare_planar_initial_state.py','scripts/analysis/simulate_weakly_nonlinear_planar_rod.py',
           'scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_spatial_rod.py',
           'scripts/lib/planar_prepared_initial_state.py',config['pilot_config'])
    manifest=read_json(source/'manifest.json')
    item={'version':ONE_T1_VERSION,'config':config,'config_sha256':sha(config_path),
          'source':{'bundle':config['source_bundle'],'manifest_sha256':sha(source/'manifest.json'),'artifact_hashes':manifest['artifact_hashes']},
          'code_hashes':{name:sha(ROOT/name) for name in paths},'python':sys.version,
          'dependencies':{name:importlib.metadata.version(name) for name in ('numpy','scipy','matplotlib')},
          'blas_threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')}}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item


def load_one_T1_source(config):
    """Own historical manifests and saved initial coordinates, no preparation."""
    if config['degrees']!=[48,64] or config['cases']!=[[48,'tight'],[64,'tight'],[64,'allowed_extra']] or config['periods']!=1 or config['amplitude_over_h']!=.05:
        raise ValueError('One-T1 authorization is restricted to the declared three cases')
    if config['execution_mode']!='EXPLORATORY_NOT_CERTIFIED' or config['state_admitted'] is not False or config['projection_policy']!='common_endpoint_constrained_L2':
        raise ValueError('Exploratory qualification and fixed initial policy required')
    if config['budget']['maximum_integrations']!=3 or config['budget']['numerical_wall_seconds']>1200:
        raise ValueError('Declared one-T1 computation budget exceeded')
    source=ROOT/config['source_bundle'];summary=validate_cache(source);manifest=read_json(source/'manifest.json')
    if summary['execution_mode']!=config['execution_mode'] or summary['state_admitted_flag'] is not False:
        raise ValueError('Source execution qualification mismatch')
    frozen=('scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_spatial_rod.py',config['pilot_config'])
    for name in frozen:
        if sha(ROOT/name)!=manifest['identity']['code_hashes'][name]:raise ValueError('Frozen physics/time input mismatch: '+name)
    q0={};v0={};settings={};times=None
    for degree in config['degrees']:
        with np.load(source/'initial_projection'/f'p{degree}.npz',allow_pickle=False) as saved:
            q0[degree]=saved['q'].copy();v0[degree]=saved['velocity'].copy()
        if q0[degree].shape!=(4*(degree-1),) or not np.all(np.isfinite(q0[degree])) or not np.array_equal(v0[degree],np.zeros_like(q0[degree])):
            raise ValueError('Invalid saved initial coordinates/velocities')
    for degree,level in config['cases']:
        name=f'p{degree}_{level}';case=read_json(source/'short_controls'/name/'case.json');settings[name]=case
        if case['p']!=degree or case['ndof']!=4*(degree-1) or case['nq']!=2*degree+1 or case['projection_policy']!=config['projection_policy'] or case['status']!='PASS':
            raise ValueError('Historical case contract mismatch: '+name)
        with np.load(source/'short_controls'/name/'trajectory.npz',allow_pickle=False) as saved:
            tt=saved['time'].copy()
            if not np.array_equal(saved['q'][0],q0[degree]) or not np.array_equal(saved['velocity'][0],v0[degree]):raise ValueError('Initial/history mismatch: '+name)
        if times is None:times=tt
        if not np.array_equal(tt,times) or np.any(np.diff(tt)<=0) or tt[0]!=0 or tt[-1]!=case['time_end']:raise ValueError('Historical actual timestamps mismatch')
    return {'summary':summary,'manifest':manifest,'source':source,'q0':q0,'v0':v0,'settings':settings,'short_time':times,
            'pilot':read_json(ROOT/config['pilot_config']),
            'provenance':{'source_bundle':config['source_bundle'],'manifest_sha256':sha(source/'manifest.json'),'own_manifest_verified':True,
                'initial_npz_sha256':{str(p):sha(source/'initial_projection'/f'p{p}.npz') for p in config['degrees']},
                'new_projection_BVP_MP_eigen_symbolic_calls':0,'field_order':['u','w','theta','c']}}


def one_T1_time_grid(source,t_end):
    sampling=source['summary']['sampling'];short=source['short_time'];bound=sampling['omega_upper_bound']
    per=sampling['samples_per_bound_period'];n=int(np.ceil(short[-1]*bound/(2*np.pi)*per))
    step=short[-1]/n
    base=np.arange(int(np.floor(t_end/step))+1,dtype=float)*step
    checkpoints=t_end*np.array([0.,.1,.25,.5,.75,1.])
    times=np.unique(np.r_[base[base<=t_end],short,checkpoints,t_end])
    if np.any(np.diff(times)<=0) or times[-1]!=t_end:raise ValueError('Invalid one-T1 output grid')
    return times,{'samples':len(times),'omega_upper_bound':bound,'samples_per_bound_period':per,'bound_step_reused':step,
                  'source_timestamps_preserved':bool(np.array_equal(times[np.searchsorted(times,short)],short)),
                  'period_fractions':[0.,.1,.25,.5,.75,1.],'output_grid_not_time_error_control':True}


def one_T1_history(bundle,name):
    path=Path(bundle)/'cases'/name;case=read_json(path/'case.json')
    times=np.load(path/'time.npy',mmap_mode='r');state=np.load(path/'state.npy',mmap_mode='r')
    n=case['actual_valid_rows'];d=case['ndof']
    if len(times)!=n or state.shape[0]<n or state.shape[1]!=2*d or (n and times[-1]!=case['time_end']):raise ValueError('Actual prefix metadata mismatch')
    return {'time':times,'q':state[:n,:d],'velocity':state[:n,d:],'case':case}


def one_T1_compare(a,b,da,db,pilot,T1,fixed,outpath,window_fractions=(.1,.25,.5,.75,1.),*,required_end=None):
    """Historical all-eight norms plus sampled time curves, in bounded blocks."""
    n=min(len(a['time']),len(b['time']));times=np.asarray(a['time'][:n])
    if n==0:return {'status':'NOT_RUN'}
    if not np.array_equal(times,b['time'][:n]):raise ValueError('Comparison requires shared actual timestamps')
    xi,weights=leg.leggauss(pilot['spatial']['comparison_quadrature']);x=(xi+1)*da.length/2;weights*=da.length/2
    l2=np.empty((n,8));peak=np.empty((n,8));location=np.empty((n,8));ref_l2=np.zeros(8);ref_peak=np.zeros(8)
    for partidx,part in enumerate(('q','velocity')):
        sl=slice(4*partidx,4*partidx+4)
        for start in range(0,n,256):
            stop=min(start+256,n);aa=da.reconstruct_series(a[part][start:stop],x);bb=db.reconstruct_series(b[part][start:stop],x)
            dif=aa-bb;ab=abs(dif)
            l2[start:stop,sl]=np.sqrt(np.einsum('tif,i,tif->tf',dif,weights,dif))
            peak[start:stop,sl]=ab.max(axis=1);location[start:stop,sl]=x[ab.argmax(axis=1)]
            ref_l2[sl]=np.maximum(ref_l2[sl],np.sqrt(np.einsum('tif,i,tif->tf',bb,weights,bb)).max(axis=0))
            ref_peak[sl]=np.maximum(ref_peak[sl],abs(bb).max(axis=(0,1)))
    floors=np.r_[np.repeat(pilot['gates']['relative_numerical_floor']*max(ref_l2[:4].max(),1e-30),4),np.repeat(pilot['gates']['relative_numerical_floor']*max(ref_l2[4:].max(),1e-30),4)]
    scale_l2=np.maximum(ref_l2,floors);scale_max=np.maximum(ref_peak,floors)
    cum_l2=np.maximum.accumulate(l2,axis=0);cum_peak=np.maximum.accumulate(peak,axis=0)
    names=[part+'_'+f for part in ('q','velocity') for f in ('u','w','theta','c')];rows={};windows=[]
    for k,name in enumerate(names):
        field=name.split('_')[-1];tol=pilot['gates']['w_theta_relative' if field in ('w','theta') else 'u_c_relative']
        i=int(l2[:,k].argmax());j=int(peak[:,k].argmax());rl=float(cum_l2[-1,k]/scale_l2[k]);rm=float(cum_peak[-1,k]/scale_max[k])
        rows[name]={'max_time_L2_difference':float(cum_l2[-1,k]),'max_space_time_difference':float(cum_peak[-1,k]),
            'reference_max_time_L2':float(ref_l2[k]),'reference_max_space_time':float(ref_peak[k]),'relative_L2':rl,'relative_max':rm,
            'numerical_floor':float(floors[k]),'floor_limited':bool(ref_l2[k]<=floors[k]),'tolerance':tol,'pass':rl<=tol and rm<=tol,
            'fixed_physical_scale':float(fixed[k]),'fixed_scaled_L2':float(cum_l2[-1,k]/(fixed[k]*np.sqrt(da.length))),
            'fixed_scaled_max':float(cum_peak[-1,k]/fixed[k]),'L2_peak_time_tau':float(times[i]/T1),
            'max_peak_time_tau':float(times[j]/T1),'max_peak_s_over_L':float(location[j,k]/da.length)}
        for tau in window_fractions:
            stop=int(np.searchsorted(times,tau*T1,side='right'))
            if not stop:continue
            index=stop-1
            if times[index]<tau*T1:continue
            windows.append({'component':name,'end_tau':tau,'absolute_cumulative_L2':float(cum_l2[index,k]),
                            'absolute_cumulative_max':float(cum_peak[index,k]),'relative_L2_full_scale':float(cum_l2[index,k]/scale_l2[k]),
                            'relative_max_full_scale':float(cum_peak[index,k]/scale_max[k]),'fixed_scaled_L2':float(cum_l2[index,k]/(fixed[k]*np.sqrt(da.length))),
                            'fixed_scaled_max':float(cum_peak[index,k]/fixed[k])})
    save_npz(outpath,time=times,d_L2=l2,d_max=peak,max_location_s=location,cumulative_L2=cum_l2,cumulative_max=cum_peak,
             full_reference_L2=ref_l2,full_reference_max=ref_peak,normalization_L2=scale_l2,normalization_max=scale_max,fixed_physical_scales=fixed)
    component_pass=all(r['pass'] for r in rows.values())
    target_end=T1 if required_end is None else required_end
    coverage=bool(times[-1]==target_end and n==len(a['time'])==len(b['time']))
    return {'status':'PASS' if component_pass and coverage else 'PARTIAL','horizon_complete':coverage,
            'required_end':float(target_end),'component_gates_pass':component_pass,'fields':rows,'windows':windows,
            'time_end':float(times[-1]),'samples':n,'normalization':'one full-comparison-horizon own scale; old per-field floor, fixed physical scale; no phase alignment',
            'maximum_semantics':'sampled Gauss100 spatial/time maxima; not continuum supremum'}


def one_T1_case_diagnostics(d,hist,pilot,bg,path,fractions):
    n=len(hist['time']);energy=np.lib.format.open_memmap(path/'energy.npy',mode='w+',dtype=float,shape=(n,))
    observations=np.empty((n,6));norms=np.empty((n,4));velocity_norms=np.empty((n,4));safety={};min_mass=1.;max_condition=1.
    for start in range(0,n,256):
        stop=min(n,start+256);q=hist['q'][start:stop];v=hist['velocity'][start:stop]
        if not(np.all(np.isfinite(q)) and np.all(np.isfinite(v))):raise ArithmeticError('NONFINITE_SAVED_OUTPUT')
        f=d.reconstruct_series(q);vv=d.reconstruct_series(v);g=d.reconstruct_series(q,derivative=1)
        obs=d.reconstruct_series(q,[d.length/4,d.length/2])
        observations[start:stop]=np.column_stack((obs[:,1,1],obs[:,0,2],obs[:,0,0],obs[:,0,3],obs[:,1,3],obs[:,1,2]))
        norms[start:stop]=np.sqrt(np.einsum('tif,i,tif->tf',f,d.weights,f));velocity_norms[start:stop]=np.sqrt(np.einsum('tif,i,tif->tf',vv,d.weights,vv))
        c,theta=f[:,:,3],f[:,:,2];us,ws=g[:,:,0],g[:,:,1]
        values={'min_one_plus_c':float((1+c).min()),'max_abs_c':float(abs(c).max()),'max_abs_theta':float(abs(theta).max()),
            'max_abs_axial_gradient':float(abs(us).max()),'max_abs_transverse_gradient':float(abs(ws).max()),'max_L_abs_curvature':float(d.length*abs(g[:,:,2]).max())}
        for key,val in values.items():safety[key]=(min(safety.get(key,val),val) if key=='min_one_plus_c' else max(safety.get(key,val),val))
        gamma1=us+ws*theta-theta**2/2-us*theta**2/2-ws*theta**3/6+theta**4/24
        gamma2=ws-theta-us*theta-ws*theta**2/2+theta**3/6+us*theta**3/6
        safety['max_quartic_Gamma1']=max(safety.get('max_quartic_Gamma1',0),float(abs(gamma1).max()))
        safety['max_quartic_Gamma2']=max(safety.get('max_quartic_Gamma2',0),float(abs(gamma2).max()))
        lo=min(1.,float(((1+c)**2).min()));hi=max(1.,float(((1+c)**2).max()))
        min_mass=min(min_mass,lo);max_condition=max(max_condition,hi/lo)
        for i,(a,b) in enumerate(zip(q,v),start):energy[i]=d.energy(a,b)
    energy.flush();drift=float(np.max(abs((energy-energy[0])/energy[0])))
    safe=safety['min_one_plus_c']>pilot['safety']['min_one_plus_c'] and all(safety[k]<=v for k,v in pilot['safety'].items() if k not in ('min_one_plus_c','min_relative_mass_eigenvalue'))
    stats={'initial_energy':float(energy[0]),'relative_energy_drift_max':drift,'mass_lower_bound_min':min_mass,'mass_condition_bound_max':max_condition,
           'safety_extrema':safety,'safety_pass':safe,'energy_and_mass_pass':safe and min_mass>=pilot['safety']['min_relative_mass_eigenvalue'] and drift<=pilot['gates']['energy_relative_drift'],
           'mass_bound_method':'same weighted-Gram Loewner bounds, no eigensolves'}
    save_npz(path/'observations.npz',time=hist['time'],observations=observations,field_L2_norms=norms,velocity_L2_norms=velocity_norms)
    ix=[int(np.searchsorted(hist['time'],tau*bg['T1'])) for tau in fractions if tau*bg['T1']<=hist['time'][-1]]
    x=np.linspace(0,d.length,501)
    save_npz(path/'snapshots.npz',s=x,time=hist['time'][ix],fields=d.reconstruct_series(hist['q'][ix],x),velocities=d.reconstruct_series(hist['velocity'][ix],x))
    return stats


def run_one_T1(config,bundle):
    from types import SimpleNamespace
    import shutil,csv
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    started=time.perf_counter();deadline=started+config['budget']['numerical_wall_seconds'];source=load_one_T1_source(config)
    runner.load_runtime();pilot=source['pilot'];old=source['summary'];bg=old['background'];coeff=rod.RodCoefficients(**old['coefficients'])
    audit_path=ROOT/pilot['audit_bundle']/'result.json'
    if sha(audit_path)!=source['manifest']['identity']['sources']['audit']['result_sha256']:raise ValueError('Frozen action archive hash mismatch')
    pol=read_json(audit_path)['polynomials']
    model=SimpleNamespace(T4=rod.Polynomial.deserialize(pol['T4']),V4=rod.Polynomial.deserialize(pol['V4']),
        residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol['residuals_A']),symbols={a:rod.Polynomial.symbol(a) for a in rod.SYMBOL_ORDER})
    discs={p:PlanarGalerkin(coeff,p,model=model) for p in config['degrees']}
    T1=bg['T1'];times,sampling=one_T1_time_grid(source,T1)
    for p,d in discs.items():
        q=source['q0'][p]
        with np.load(source['source']/'initial_projection'/f'p{p}.npz') as saved:
            raw=d.raw_coefficients(q);scale=max(float(np.max(abs(saved['raw_roundtrip']))),1e-30)
            if not np.allclose(raw,saved['raw_roundtrip'],rtol=2e-11,atol=2e-11*scale):raise ValueError('Saved whitening/basis contract mismatch')
        runner.safety_check(d,q,pilot['safety'])
        if not np.all(np.isfinite(d.rhs(0,np.r_[q,source['v0'][p]]))) or np.any(d.reconstruct(q,[0.,d.length])!=0):raise ValueError('Initial RHS or essential BC inconsistency')
    for p,level in config['cases']:
        name=f'p{p}_{level}';settings=runner.time_settings(discs[p],config['amplitude_over_h']*pilot['material_geometry']['h'],level,pilot);prev=source['settings'][name]
        if settings['rtol']!=prev['rtol'] or settings['max_step']!=prev['max_step'] or not np.array_equal(settings['atol'],prev['atol']):raise ValueError('Actual time settings mismatch: '+name)
    statuses={key:'NOT_RUN' for key in ('NLSP_PREPARED_ONE_T1_EXECUTION','NLSP_PREPARED_ONE_T1_PREFIX_REGRESSION','NLSP_PREPARED_ONE_T1_SPATIAL_CHECK','NLSP_PREPARED_ONE_T1_TEMPORAL_CHECK','NLSP_PREPARED_ONE_T1_ENERGY_AND_MASS','NLSP_PREPARED_ONE_T1_FEASIBILITY')}
    summary={'schema':config['schema'],'config':config,'background':bg,'coefficients':old['coefficients'],'source_provenance':source['provenance'],
        'source_strict_table':old['strict_table'],'source_strict_status':old['statuses']['NLSP_STRICT_INITIAL_VERIFICATION'],
        'execution_mode':config['execution_mode'],'state_admitted_flag':False,'projection_policy':config['projection_policy'],
        'sampling':sampling,'statuses':statuses,'cases':{},'prefix_regression':{},'new_ODE_integrations':0,
        'forecast_integration_seconds':{name:case['integration_seconds']*10 for name,case in source['settings'].items()},
        'new_projection_MP_BVP_eigen_symbolic_calls':0,'no_new_dynamic_constraints':True,'qualification':'One LINEAR period horizon; no periodic orbit or continuum truth claim'}
    save_npz(bundle/'time_grid.npz',time=times);write_json(bundle/'summary.json',summary)
    shutil.copy2(source['source']/'common_initial_state.npz',bundle/'common_initial_state.npz')
    for p in config['degrees']:save_npz(bundle/'initial_projection'/f'p{p}.npz',q=source['q0'][p],velocity=source['v0'][p])
    print(json.dumps({'stage':'ONE_T1_PRE_RUN','bundle':str(bundle),'T1':T1,'samples':len(times),'execution':config['execution_mode'],'forecast':summary['forecast_integration_seconds']}),flush=True)
    fixed=np.array([config['amplitude_over_h']*bg['h0']]*2+[config['amplitude_over_h']*bg['h0']/bg['L']]*2);fixed=np.r_[fixed,fixed*bg['omega1']]
    for p,level in config['cases']:
        if time.perf_counter()>=deadline:break
        name=f'p{p}_{level}';path=bundle/'cases'/name;path.mkdir(parents=True,exist_ok=True);d=discs[p]
        storage=np.lib.format.open_memmap(path/'state.npy',mode='w+',dtype=float,shape=(len(times),2*d.ndof))
        summary['new_ODE_integrations']+=1;print(json.dumps({'stage':'START_ONE_T1','case':name}),flush=True)
        history,stats=runner.integrate_case(d,None,bg,pilot,config['amplitude_over_h'],level,times,deadline,initial_coordinates=source['q0'][p],history_buffer=storage)
        storage.flush();np.save(path/'time.npy',times[:len(history)]);stats.update(actual_valid_rows=len(history),storage_allocated_rows=len(times),
            execution_mode=config['execution_mode'],state_admitted=False,projection_policy=config['projection_policy'],periods=1.,T1=T1,omega1=bg['omega1'])
        write_json(path/'case.json',stats);save_npz(path/'internal_steps.npz',dt=np.asarray(stats['internal_time_steps']))
        del history,storage
        hist=one_T1_history(bundle,name)
        stats.update(one_T1_case_diagnostics(d,hist,pilot,bg,path,config['sampling']['period_fractions']));write_json(path/'case.json',stats)
        sourcepath=source['source']/'short_controls'/name/'trajectory.npz'
        with np.load(sourcepath) as saved:
            prefix_time=saved['time'].copy();oldhist={part:saved[part].copy() for part in ('q','velocity')};oldhist['time']=prefix_time;oldenergy=saved['energy'].copy()
        indexes=np.searchsorted(hist['time'],prefix_time);available=indexes<len(hist['time']);indexes=indexes[available];prefix_time=prefix_time[available]
        if len(indexes) and not np.array_equal(hist['time'][indexes],prefix_time):raise ValueError('Source prefix timestamps not found exactly')
        newprefix={'time':prefix_time,'q':hist['q'][indexes],'velocity':hist['velocity'][indexes]};oldhist={k:v[:len(indexes)] for k,v in oldhist.items()}
        prefix=one_T1_compare(newprefix,oldhist,d,d,pilot,T1,fixed,path/'prefix_differences.npz',(.1,),required_end=source['short_time'][-1])
        en=np.load(path/'energy.npy',mmap_mode='r')[indexes]
        prefix['relative_energy_difference_max']=float(np.max(abs(en-oldenergy[:len(indexes)]))/oldenergy[0])
        prefix['initial_q_v_exact']=bool(np.array_equal(newprefix['q'][0],source['q0'][p]) and np.array_equal(newprefix['velocity'][0],source['v0'][p]))
        prefix['source_interval_complete']=len(indexes)==len(source['short_time']);prefix['time_settings_identical']=True
        prefix['pass']=prefix['status']=='PASS' and prefix['relative_energy_difference_max']<=pilot['gates']['energy_relative_drift'] and prefix['initial_q_v_exact'] and prefix['source_interval_complete']
        summary['prefix_regression'][name]=prefix;summary['cases'][name]=stats
        write_json(path/'prefix_regression.json',prefix);write_json(bundle/'summary.json',summary)
        print(json.dumps({'stage':'END_ONE_T1','case':name,'runtime':stats['integration_seconds'],'nfev':stats['nfev'],'end_tau':stats['time_end']/T1,'energy_drift':stats['relative_energy_drift_max'],'prefix_pass':prefix['pass']}),flush=True)
        del hist,oldhist,newprefix
        if stats['status']!='PASS' or not prefix['pass']:break
    pairs=(('spatial','p48_tight','p64_tight'),('temporal','p64_tight','p64_allowed_extra'))
    for label,a,b in pairs:
        if a not in summary['cases'] or b not in summary['cases']:summary[label+'_comparison']={'status':'NOT_RUN'};continue
        aa,bb=one_T1_history(bundle,a),one_T1_history(bundle,b)
        result=one_T1_compare(aa,bb,discs[int(a.split('_')[0][1:])],discs[int(b.split('_')[0][1:])],pilot,T1,fixed,bundle/(label+'_differences.npz'))
        summary[label+'_comparison']=result;write_json(bundle/(label+'_comparison.json'),result)
        for suffix,rows in (('all8',[{'component':key,**val} for key,val in result['fields'].items()]),('windows',result['windows'])):
            with (bundle/f'{label}_{suffix}.csv').open('w',newline='',encoding='utf8') as f:
                writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        del aa,bb
    complete=len(summary['cases'])==3 and all(c['status']=='PASS' for c in summary['cases'].values())
    statuses['NLSP_PREPARED_ONE_T1_EXECUTION']='COMPLETED_EXPLORATORY_NOT_CERTIFIED' if complete else 'PARTIAL'
    statuses['NLSP_PREPARED_ONE_T1_PREFIX_REGRESSION']='PASS' if len(summary['prefix_regression'])==3 and all(c['pass'] for c in summary['prefix_regression'].values()) else 'PARTIAL'
    for label in ('spatial','temporal'):statuses['NLSP_PREPARED_ONE_T1_'+label.upper()+'_CHECK']=summary[label+'_comparison']['status']
    statuses['NLSP_PREPARED_ONE_T1_ENERGY_AND_MASS']='PASS' if complete and all(c['energy_and_mass_pass'] for c in summary['cases'].values()) else 'PARTIAL'
    statuses['NLSP_PREPARED_ONE_T1_FEASIBILITY']='COMPLETED_EXPLORATORY_NOT_CERTIFIED' if complete and statuses['NLSP_PREPARED_ONE_T1_PREFIX_REGRESSION']==statuses['NLSP_PREPARED_ONE_T1_ENERGY_AND_MASS']=='PASS' else 'PARTIAL'
    summary['runtime']={'numerical_wall_seconds':time.perf_counter()-started,'integration_seconds':sum(c['integration_seconds'] for c in summary['cases'].values()),'limit_seconds':1200,'ODE_integrations':summary['new_ODE_integrations'],'projection_MP_BVP_eigen_symbolic_calls':0}
    summary['stop_reason']='Authorized one-linear-period horizon complete or accepted prefix saved; no automatic extension'
    write_json(bundle/'summary.json',summary);return summary


def one_T1_plot_only(bundle,summary):
    import matplotlib.pyplot as plt
    figs=Path(bundle)/'figures';figs.mkdir(exist_ok=True);T1=summary['background']['T1']
    def finish(fig,name):
        for ext in ('pdf','png'):fig.savefig(figs/(name+'.'+ext),dpi=200,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(11,5.6),layout='constrained')
    labels=('w(L/2)','theta(L/4)','u(L/4)','c(L/4)','c(L/2)','||c|| L2')
    for name,style in (('p48_tight','-'),('p64_tight','--')):
        if name not in summary['cases']:continue
        with np.load(Path(bundle)/'cases'/name/'observations.npz') as data:
            for k,ax in enumerate(axes.flat):
                values=data['observations'][:,k] if k<5 else data['field_L2_norms'][:,3]
                ax.plot(data['time']/T1,values,style,label=name.replace('_',' '))
                ax.set(xlabel='t/T1',ylabel=labels[k]);ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8,frameon=False);finish(fig,'prepared_motion_one_T1')
    fig,axes=plt.subplots(2,4,figsize=(12,5.5),layout='constrained')
    for label,style in (('spatial','-'),('temporal','--')):
        path=Path(bundle)/(label+'_differences.npz')
        if not path.exists():continue
        with np.load(path) as data:
            scale=summary['spatial_comparison']['fields'] if summary.get('spatial_comparison',{}).get('fields') else summary[label+'_comparison']['fields']
            keys=[part+'_'+f for part in ('q','velocity') for f in ('u','w','theta','c')]
            for k,ax in enumerate(axes.flat):
                ref=max(scale[keys[k]]['reference_max_time_L2'],scale[keys[k]]['numerical_floor'])
                ax.semilogy(data['time']/T1,np.maximum(data['d_L2'][:,k]/ref,1e-18),style,lw=.8,label=label)
                ax.set(xlabel='t/T1',ylabel=keys[k].replace('q_','')+' L2 difference');ax.grid(alpha=.2)
    fig.suptitle('L2 differences / common full-horizon p64 tight scales',fontsize=10)
    axes.flat[0].legend(fontsize=8,frameon=False);finish(fig,'prepared_differences_one_T1')
    fig,axes=plt.subplots(1,2,figsize=(9,3.4),layout='constrained')
    for name,case in summary['cases'].items():
        tt=np.load(Path(bundle)/'cases'/name/'time.npy',mmap_mode='r');energy=np.load(Path(bundle)/'cases'/name/'energy.npy',mmap_mode='r')
        axes[0].plot(tt/T1,(energy-energy[0])/energy[0],lw=.8,label=name.replace('_',' '))
    if 'spatial_comparison' in summary and (Path(bundle)/'spatial_differences.npz').exists():
        with np.load(Path(bundle)/'spatial_differences.npz') as data:
            for k,f in enumerate(('u','w','theta','c','u_t','w_t','theta_t','c_t')):
                axes[1].plot(data['time']/T1,data['cumulative_max'][:,k]/data['normalization_max'][k],lw=.8,label=f)
    axes[0].set(xlabel='t/T1',ylabel='relative energy drift');axes[1].set(xlabel='t/T1',ylabel='cumulative spatial max / full-horizon scale')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=7,frameon=False)
    finish(fig,'prepared_quality_one_T1');return {'ODE_integrations':0,'BVP_solves':0,'eigensolves':0,'symbolic_derivations':0}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true')
    mode.add_argument('--report-only',type=Path)
    mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path)
    horizon=parser.add_mutually_exclusive_group()
    horizon.add_argument('--feasibility',action='store_true',help='Explicit bounded exploratory authorization, strict tests retained')
    horizon.add_argument('--one-T1',dest='one_T1',action='store_true',help='Reuse saved q0 for exactly three exploratory runs to one linear period')
    parser.add_argument('--output-dir',type=Path)
    args=parser.parse_args(argv)
    if (args.feasibility or args.one_T1) and not args.compute:
        parser.error('Horizon selection requires --compute; report/plot read the saved mode')
    args.config=args.config or (ONE_T1_CONFIG if args.one_T1 else FEASIBILITY_CONFIG if args.feasibility else CONFIG)
    args.output_dir=args.output_dir or (ONE_T1_OUTPUT if args.one_T1 else FEASIBILITY_OUTPUT if args.feasibility else OUTPUT)
    zeros={'BVP_solves':0,'eigendecompositions':0,'analytic_history_evaluations':0,'ODE_integrations':0,'symbolic_derivations':0}
    if args.report_only or args.plot_only:
        bundle=args.report_only or args.plot_only
        summary=validate_cache(bundle)
        if args.plot_only:
            plot_only(bundle)
        print(json.dumps({'bundle':str(bundle),'statuses':summary['statuses'],'this_run_counters':zeros},indent=2))
        return summary
    key,item=(one_T1_identity(args.config) if args.one_T1 else feasibility_identity(args.config) if args.feasibility else identity(args.config));bundle=args.output_dir/key
    if (bundle/'manifest.json').exists():
        summary=validate_cache(bundle,item)
        print(json.dumps({'bundle':str(bundle),'cache_hit':True,'statuses':summary['statuses'],'this_run_counters':zeros},indent=2))
        return summary
    bundle.mkdir(parents=True,exist_ok=True)
    summary=(run_one_T1 if args.one_T1 else run_feasibility if args.feasibility else run_compute)(read_json(args.config),bundle)
    write_json(bundle/'manifest.json',manifest_for(bundle,item))
    plot_only(bundle)
    write_json(bundle/'manifest.json',manifest_for(bundle,item))
    print(json.dumps({'bundle':str(bundle),'cache_hit':False,'statuses':summary['statuses'],'runtime':summary['runtime'],'stop_reason':summary.get('stop_reason')},indent=2))
    return summary


if __name__ == '__main__':
    main()
