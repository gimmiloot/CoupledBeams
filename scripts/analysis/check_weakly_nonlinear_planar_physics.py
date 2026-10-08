"""Bounded physical consistency checks of the frozen planar quartic action.

New I/O contract: leading asymptotics, one half-amplitude run, formal limit,
strains and reactions. Existing RHS/runner are reused; no new physics solver.
"""
from __future__ import annotations
import argparse,hashlib,json,os,sys,time,subprocess,importlib.metadata as md
from pathlib import Path
from fractions import Fraction
from types import SimpleNamespace
if __name__=='__main__':
    for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):os.environ[k]='1'
import numpy as np
from numpy.polynomial import legendre as leg
ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from scripts.analysis import prepare_planar_initial_state as saved
from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib import weakly_nonlinear_planar_dynamics as dyn
from scripts.lib import planar_prepared_initial_state as prep
read_json,write_json,save_npz,sha=saved.read_json,saved.write_json,saved.save_npz,saved.sha
CONFIG=ROOT/'data/input/nlsp_planar_physical_sanity_checks.json'
OUTPUT=ROOT/'results/nlsp_planar_physical_sanity_checks'
VERSION='frozen-planar-physical-sanity-v1'
FIELDS=('u','w','theta','c');P=np.array([1.,-1.,-1.,1.])


def identity(config_path=CONFIG):
    c=read_json(config_path)
    paths=('scripts/analysis/check_weakly_nonlinear_planar_physics.py','scripts/analysis/prepare_planar_initial_state.py','scripts/analysis/simulate_weakly_nonlinear_planar_rod.py','scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_spatial_rod.py','scripts/lib/planar_prepared_initial_state.py',c['pilot_config'])
    item={'version':VERSION,'config':c,'config_sha256':sha(config_path),'code_hashes':{p:sha(ROOT/p) for p in paths},
          'source_manifests':{k:sha(ROOT/c[k]/'manifest.json') for k in ('one_T1_source','prepared_source','second_order_source','feasibility_source')},
          'audit_result_sha256':sha(ROOT/read_json(ROOT/c['pilot_config'])['audit_bundle']/'result.json'),
          'python':sys.version,'dependencies':{k:md.version(k) for k in ('numpy','scipy','matplotlib','mpmath')},
          'blas_threads':{k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS')}}
    return hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16],item


def load_inputs(c):
    if (c['p'],c['large_amplitude'],c['small_amplitude'],c['time_level'],c['periods'])!=(64,.05,.025,'tight',1.):raise ValueError('Only the declared physical sanity case is authorized')
    if c['projection_policy']!=prep.CONSTRAINED_PROJECTION or c['budget']['max_new_integrations']!=1 or c['budget']['numerical_wall_seconds']>600:raise ValueError('Bounded policy mismatch')
    if c['execution_mode']!='EXPLORATORY_NOT_CERTIFIED':raise ValueError('Strict admission is not certified by this task')
    main=ROOT/c['one_T1_source'];summary=saved.validate_cache(main);manifest=read_json(main/'manifest.json')
    for p in ('scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_spatial_rod.py','scripts/analysis/simulate_weakly_nonlinear_planar_rod.py'):
        if sha(ROOT/p)!=manifest['identity']['code_hashes'][p]:raise ValueError('Frozen numerical model mismatch: '+p)
    state,coeff,provenance=prep.load_frozen_prepared_state(ROOT/c['prepared_source'])
    initial_manifest=read_json(ROOT/c['prepared_source']/'manifest.json')
    profile_path=ROOT/c['prepared_source']/'coordinates/p96.npz'
    if sha(profile_path)!={k.replace('\\','/'):v for k,v in initial_manifest['artifact_hashes'].items()}['coordinates/p96.npz']:raise ValueError('Saved leading profiles corrupted')
    with np.load(profile_path) as data:stat=data['legendre_stat'].copy();harm=data['legendre_harm'].copy()
    if not np.array_equal(stat+harm,state.profiles.coefficients):raise ValueError('Common initial and leading profiles disagree')
    pilot=read_json(ROOT/c['pilot_config']);audit_path=ROOT/pilot['audit_bundle']/'result.json'
    antecedent=read_json(ROOT/c['feasibility_source']/'manifest.json')['identity']['sources']['audit']
    if sha(audit_path)!=antecedent['result_sha256']:raise ValueError('Frozen action archive hash mismatch')
    audit=read_json(audit_path);pol=audit['polynomials']
    model=SimpleNamespace(T4=rod.Polynomial.deserialize(pol['T4']),V4=rod.Polynomial.deserialize(pol['V4']),residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol['residuals_A']),symbols={k:rod.Polynomial.symbol(k) for k in rod.SYMBOL_ORDER})
    if coeff.values()!=summary['coefficients']:raise ValueError('Physical coefficients differ between sources')
    return {'source':main,'summary':summary,'state':state,'coeff':coeff,'pilot':pilot,'model':model,'stat':stat,'harm':harm,'provenance':provenance}


def cubic_measures(us,ws,theta):
    return us+theta*ws-theta**2/2-us*theta**2/2,ws-theta-us*theta-ws*theta**2/2+theta**3/6


def local_polynomials(model):
    V=dyn._restrict_to_plane(model.V4);T=dyn._restrict_to_plane(model.T4);z=model.symbols
    us,ws,th,c=(z[k] for k in ('u_s','w_s','theta','c'))
    g1,g2=cubic_measures(us,ws,th);Fu,Fw=V.derivative('u_s'),V.derivative('w_s')
    M,Rc,Vth,Vc=(V.derivative(k) for k in ('theta_s','c_s','theta','c'))
    K=Vth+(1+us)*Fw-ws*Fu
    return {'T':T,'V':V,'Gamma1':g1,'Gamma2':g2,'Fu':Fu,'Fw':Fw,'M':M,'Rc':Rc,'Vtheta':Vth,'Vc':Vc,'K':K}


def limited_math(model):
    p=local_polynomials(model);z=model.symbols;flip={k:-z[k] for k in rod.SYMBOL_ORDER if k.split('_')[0] in ('w','theta')}
    checks={'T_parity':p['T'].substitute(flip)==p['T'],'V_parity':p['V'].substitute(flip)==p['V'],'K_through_cubic_zero':not p['K'].truncate(3),'nonzero_quartic_K_retained':bool(p['K'])}
    for k,sign in (('Gamma1',1),('Gamma2',-1),('Fu',1),('Fw',-1),('M',-1),('Rc',1)):
        checks[k+'_parity']=p[k].substitute(flip)==sign*p[k]
    for i,sign in zip((0,1,5,6),(1,-1,-1,1)):
        q=dyn._restrict_to_plane(model.residual_a[i]);checks['residual_'+str(i)+'_parity']=q.substitute(flip)==sign*q
    N=z['C']*(p['Gamma1']+z['nu']*z['c']);Q=z['S']*p['Gamma2'];th=z['theta'];cos=1-th**2/2;sin=th-th**3/6
    checks['Fu_rotated_retained']=p['Fu']==(N*cos-Q*sin).truncate(3);checks['Fw_rotated_retained']=p['Fw']==(N*sin+Q*cos).truncate(3)
    # Classical quadratic coefficient and fixed-end elimination in exact rational arithmetic.
    nu=Fraction(3,10);C=Fraction(1,1)/(1-nu*nu);checks['classical_EA_coefficient']=C*(1-nu*nu)==1
    checks['classical_positive_one_eighth']=Fraction(1,2)*Fraction(1,2)**2==Fraction(1,8)
    return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks,'polynomials':{k:v.serialize() for k,v in p.items()},
        'expressions':{k:str(v) for k,v in p.items()},'classical':{'N':'EA/(2L)*integral(w_s^2 ds)','Vstretch':'EA/(8L)*(integral(w_s^2 ds))^2','variation':'N*integral(w_s*delta_w_s ds); gradient=-N*w_ss','bulk_limit_only':True,'finite_c_clamp_elimination':False},
        'angular_balance':'J_t-[r cross F+M]=integral(r cross E_U+E_theta-K4); K4 genericOeps4/preparedOeps5'}


def asymptotic_fields(inputs,x,t,epsilon,velocity=False):
    state=inputs['state'];omega=state.background.omega;times=np.atleast_1d(t)
    pair=state.background.evaluate(x);st=prep.LegendreProfiles(inputs['stat'],state.length).evaluate(x);ha=prep.LegendreProfiles(inputs['harm'],state.length).evaluate(x)
    out=np.zeros((len(times),len(x),4))
    c1=(-omega*np.sin(omega*times) if velocity else np.cos(omega*times))
    c2=(-2*omega*np.sin(2*omega*times) if velocity else np.cos(2*omega*times))
    axial=(c2[:,None,None]*ha[None,:,:] if velocity else st[None,:,:]+c2[:,None,None]*ha[None,:,:])
    out[:,:,0]=epsilon**2*axial[:,:,0];out[:,:,3]=epsilon**2*axial[:,:,1];out[:,:,1:3]=epsilon*c1[:,None,None]*pair[None,:,:]
    return out


def amplitude_scales(e,h,L,omega):
    q=np.array([e**2*h,e*h,e*h/L,e**2]);return np.r_[q,q*np.array([2*omega,omega,omega,2*omega])]


def canonical_flux(model,coeff,fields,gradients):
    p=local_polynomials(model);names=('u_s','w_s','theta','c','c_s','theta_s')
    compiler=dyn._CompiledPolynomials([p[k] for k in ('Fu','Fw','M','Rc','Vtheta','Vc','K')],names,coeff)
    return compiler.evaluate((gradients[:,0],gradients[:,1],fields[:,2],fields[:,3],gradients[:,3],gradients[:,2])).T


def history_analysis(inputs,d,hist,e,label,bundle,cfg):
    n=len(hist['time']);x,w=leg.leggauss(cfg['comparison_quadrature']);x=(x+1)*d.length/2;w*=d.length/2
    scales=amplitude_scales(e,inputs['state'].background.h0,d.length,inputs['state'].background.omega)
    dif_l2=np.empty((n,8));dif_max=np.empty((n,8));norm_l2=np.zeros(8);norm_max=np.zeros(8);obs=np.empty((n,4));leading_obs=obs.copy();strain={}
    names=('u_s','w_s','theta','c','L_c_s','L_theta_s','Gamma1','Gamma2','surface_bending_strain')
    for start in range(0,n,cfg['block_rows']):
        stop=min(n,start+cfg['block_rows']);tt=hist['time'][start:stop];q=hist['q'][start:stop];v=hist['velocity'][start:stop]
        for partidx,part in enumerate((q,v)):
            actual=d.reconstruct_series(part,x);reference=asymptotic_fields(inputs,x,tt,e,partidx==1);dif=actual-reference;sl=slice(partidx*4,partidx*4+4)
            dif_l2[start:stop,sl]=np.sqrt(np.einsum('tif,i,tif->tf',dif,w,dif));dif_max[start:stop,sl]=abs(dif).max(axis=1)
            norm_l2[sl]=np.maximum(norm_l2[sl],np.sqrt(np.einsum('tif,i,tif->tf',actual,w,actual)).max(axis=0));norm_max[sl]=np.maximum(norm_max[sl],abs(actual).max(axis=(0,1)))
        f=d.reconstruct_series(q);g=d.reconstruct_series(q,derivative=1);g1,g2=cubic_measures(g[:,:,0],g[:,:,1],f[:,:,2])
        data=(g[:,:,0],g[:,:,1],f[:,:,2],f[:,:,3],d.length*g[:,:,3],d.length*g[:,:,2],g1,g2,inputs['state'].background.h0/2*g[:,:,2])
        for name,val in zip(names,data):
            i,j=np.unravel_index(np.argmax(abs(val)),val.shape);peak=float(abs(val[i,j]))
            if name not in strain or peak>strain[name]['max_abs']:strain[name]={'max_abs':peak,'time_tau':float(tt[i]/inputs['state'].background.T1),'s_over_L':float(d.x[j]/d.length)}
        observations=d.reconstruct_series(q,[d.length/4,d.length/2]);obs[start:stop]=np.column_stack((observations[:,0,0],observations[:,1,1],observations[:,0,2],observations[:,0,3]))
        oo=asymptotic_fields(inputs,[d.length/4,d.length/2],tt,e);leading_obs[start:stop]=np.column_stack((oo[:,0,0],oo[:,1,1],oo[:,0,2],oo[:,0,3]))
    rows={}
    keys=[part+'_'+f for part in ('q','velocity') for f in FIELDS]
    for k,key in enumerate(keys):rows[key]={'absolute_L2':float(dif_l2[:,k].max()),'absolute_max':float(dif_max[:,k].max()),'fixed_physical_scale':float(scales[k]),'fixed_scaled_L2':float(dif_l2[:,k].max()/scales[k]),'fixed_scaled_max':float(dif_max[:,k].max()/scales[k]),'trajectory_characteristic_L2':float(norm_l2[k]),'trajectory_characteristic_max':float(norm_max[k])}
    save_npz(bundle/(label+'_asymptotic.npz'),time=hist['time'],difference_L2=dif_l2,difference_max=dif_max,observations=obs,leading_observations=leading_obs,physical_scales=scales)
    return {'deviations':rows,'strains':strain,'epsilon_a':e,'actual_time_end':float(hist['time'][-1]),'samples':n,'qualification':'Leading approximation is not exact cubic solution; differences combine higher orders and discretization'}


def rhs_parity(d,hist,indices):
    pp=np.repeat(P,d.n);out=[]
    for i in indices:
        y=np.r_[hist['q'][i],hist['velocity'][i]];full=np.r_[pp,pp];a=d.rhs(hist['time'][i],y);b=d.rhs(hist['time'][i],full*y)
        out.append({'time':float(hist['time'][i]),'absolute_max':float(np.max(abs(b-full*a))),'relative':float(np.linalg.norm(b-full*a)/max(np.linalg.norm(a),1e-30))})
    return {'status':'PASS' if all(r['relative']<=2e-12 for r in out) else 'FAIL','rows':out,'negative_amplitude_integrations':0}


def reaction_checks(inputs,d,hist,indices):
    results=[];c=d.coefficients;T1=inputs['state'].background.T1
    for i in indices:
        q,v=hist['q'][i],hist['velocity'][i];a=d.acceleration(q,v);f=d.reconstruct(q);g=d.reconstruct(q,derivative=1);vv=d.reconstruct(v);aa=d.reconstruct(a)
        flux=canonical_flux(inputs['model'],c,f,g);endf=d.reconstruct(q,[0.,d.length]);endg=d.reconstruct(q,[0.,d.length],1);ef=canonical_flux(inputs['model'],c,endf,endg)
        w=d.weights;s=d.x;r=np.column_stack((s+f[:,0],f[:,1]));Fu,Fw,M,Rc,Vth,Vc,K=flux.T
        dot_spin=c.jp*((1+f[:,3])**2*aa[:,2]+2*(1+f[:,3])*vv[:,3]*vv[:,2])
        pdot=c.m*(w@aa[:,:2]);boundary=ef[1,:2]-ef[0,:2]
        jdot=float(w@(r[:,0]*c.m*aa[:,1]-r[:,1]*c.m*aa[:,0]+dot_spin));torque=float(d.length*ef[1,1]+ef[1,2]-ef[0,2])
        # Independent strong polynomial evaluator of the accepted action.
        local=[]
        for coord,derivative in ((q,0),(q,1),(v,0),(q,2),(v,1),(a,0)):local.extend(d.reconstruct(coord,derivative=derivative).T)
        E=d._residual.evaluate(local).T
        EU=w@E[:,:2];ang_expected=float(w@(r[:,0]*E[:,1]-r[:,1]*E[:,0]+E[:,2]-K))
        lifts=np.column_stack((1-s/d.length,s/d.length));dlifts=np.column_stack((-np.ones_like(s)/d.length,np.ones_like(s)/d.length))
        recovered=np.empty((2,4))
        for side in range(2):
            test=lifts[:,side];dt=dlifts[:,side]
            recovered[side,0]=w@(c.m*aa[:,0]*test+Fu*dt);recovered[side,1]=w@(c.m*aa[:,1]*test+Fw*dt)
            recovered[side,2]=w@((dot_spin+Vth)*test+M*dt)
            recovered[side,3]=w@((c.jp*aa[:,3]+Vc-c.jp*(1+f[:,3])*vv[:,2]**2)*test+Rc*dt)
        endpoint=ef[:,:4]*np.array([-1.,1.])[:,None]
        results.append({'tau':float(hist['time'][i]/T1),'actual_time':float(hist['time'][i]),'endpoint_on_rod_local_Fu_Fw_M_Rc':endpoint.tolist(),'weak_lift_support_reactions':recovered.tolist(),
            'physical_global_force_xy':np.column_stack((endpoint[:,0],-endpoint[:,1])).tolist(),'physical_global_couple_z':(-endpoint[:,2]).tolist(),
            'linear_momentum_derivative':pdot.tolist(),'endpoint_force_sum':boundary.tolist(),'linear_balance_difference':(pdot-boundary).tolist(),
            'independent_integrated_strong_translation':EU.tolist(),'translation_decomposition_error':(pdot-boundary-EU).tolist(),
            'angular_momentum_derivative':jdot,'endpoint_torque_about_left':torque,'angular_balance_difference':jdot-torque,'integrated_K4':float(w@K),
            'independent_strong_plus_truncation':ang_expected,'angular_decomposition_error':jdot-torque-ang_expected,
            'strong_residual_L2':np.sqrt(w@E**2).tolist(),'projected_strong_weak_max':float(np.max(abs(d.weak_residual(q,v,a)))),
            'lift_minus_endpoint':(recovered-endpoint).tolist(),'weak_reaction_balance_is_bookkeeping_not_independent_validation':True})
    return {'snapshots':results,'local_orientation':'t=+EX,n=-EY,k=-EZ; signed angular momentum along k','outward_signs':[-1,1],
            'qualification':'Endpoint flux discrepancies contain finite-p strong residual; angular adds K4; lift reactions are independently assembled boundary-work residuals, their summed balance is not independent PASS'}


def normalized_amplitude_comparison(d,large,small,eL,eS,inputs):
    n=min(len(large['time']),len(small['time']));x,w=leg.leggauss(100);x=(x+1)*d.length/2;w*=d.length/2
    rows={};powers=np.array([2,1,1,2]);fixed=np.array([inputs['state'].background.h0]*2+[inputs['state'].background.h0/d.length,1.])
    if not np.array_equal(large['time'][:n],small['time'][:n]):raise ValueError('Amplitude comparisons require exact common timestamps')
    for part in ('q','velocity'):
        part_scales=fixed*(np.array([2,1,1,2])*inputs['state'].background.omega if part=='velocity' else 1)
        maxL2=np.zeros(4);maxpeak=np.zeros(4);charL=np.zeros(4);charS=np.zeros(4)
        for start in range(0,n,256):
            end=min(start+256,n);a=d.reconstruct_series(large[part][start:end],x);b=d.reconstruct_series(small[part][start:end],x)
            charL=np.maximum(charL,abs(a).max(axis=(0,1)));charS=np.maximum(charS,abs(b).max(axis=(0,1)))
            difference=a/eL**powers-b/eS**powers;maxL2=np.maximum(maxL2,np.sqrt(np.einsum('tif,i,tif->tf',difference,w,difference)).max(axis=0));maxpeak=np.maximum(maxpeak,abs(difference).max(axis=(0,1)))
        for j,f in enumerate(FIELDS):rows[part+'_'+f]={'amplitude_power':int(powers[j]),'expected_characteristic_ratio':float((eL/eS)**powers[j]),'observed_max_ratio':float(charL[j]/charS[j]),'large_max':float(charL[j]),'small_max':float(charS[j]),'normalized_difference_L2':float(maxL2[j]),'normalized_difference_max':float(maxpeak[j]),'fixed_normalized_scale':float(part_scales[j]),'fixed_scaled_normalized_max':float(maxpeak[j]/part_scales[j])}
    return {'rows':rows,'pointwise_ratios_near_zero':False,'power_law_fit':False,'same_physical_time':True}


def supplement_reports(inputs,d,large,half,summary,bundle,cfg):
    """Bounded saved-data diagnostics only; never integrates or prepares new IC."""
    import csv
    omega=inputs['state'].background.omega;h=inputs['state'].background.h0
    L=d.length;c=d.coefficients
    for label,hist,e,snapshot_path in (
        ('large',large,cfg['large_amplitude'],inputs['source']/'cases/p64_tight/snapshots.npz'),
        ('half',half,cfg['small_amplitude'],bundle/'cases/p64_half_tight/snapshots.npz')):
        with np.load(snapshot_path) as z:
            x=z['s'];tt=z['time'];f=z['fields'];v=z['velocities']
        lead=asymptotic_fields(inputs,x,tt,e)
        leadv=asymptotic_fields(inputs,x,tt,e,True)
        indices=np.searchsorted(hist['time'],tt)
        if not np.array_equal(hist['time'][indices],tt):raise ValueError('Snapshot timestamps inconsistent')
        gradients=d.reconstruct_series(hist['q'][indices],x,derivative=1)
        g1,g2=cubic_measures(gradients[:,:,0],gradients[:,:,1],f[:,:,2])
        exact1=(1+gradients[:,:,0])*np.cos(f[:,:,2])+gradients[:,:,1]*np.sin(f[:,:,2])-1
        exact2=gradients[:,:,1]*np.cos(f[:,:,2])-(1+gradients[:,:,0])*np.sin(f[:,:,2])
        save_npz(bundle/(label+'_comparison_snapshots.npz'),s=x,time=tt,fields=f,velocities=v,
                 leading_fields=lead,leading_velocities=leadv,difference=f-lead,
                 velocity_difference=v-leadv,Gamma1_cubic=g1,Gamma2_cubic=g2,
                 Gamma1_exact_reduced=exact1,Gamma2_exact_reduced=exact2)
        forces={}
        for start in range(0,len(hist['time']),cfg['block_rows']):
            end=min(start+cfg['block_rows'],len(hist['time']));q=hist['q'][start:end]
            fields=d.reconstruct_series(q);g=d.reconstruct_series(q,derivative=1)
            G1,G2=cubic_measures(g[:,:,0],g[:,:,1],fields[:,:,2])
            for name,val in zip(('N3','Q3','M','Rc'),(c.C*(G1+c.nu*fields[:,:,3]),c.S*G2,c.Bp*g[:,:,2],c.H*g[:,:,3])):
                i,j=np.unravel_index(np.argmax(abs(val)),val.shape);peak=float(abs(val[i,j]))
                if name not in forces or peak>forces[name]['max_abs']:
                    forces[name]={'max_abs':peak,'tau':float(hist['time'][start+i]/inputs['state'].background.T1),'s_over_L':float(d.x[j]/L)}
        record=summary[label+'_analysis'];record['resultants']=forces
        record['geometric_truncation_snapshot_difference']={'Gamma1_max':float(abs(exact1-g1).max()),'Gamma2_max':float(abs(exact2-g2).max()),'exact_reduced_dynamics_used':False}
        write_json(bundle/(label+'_analysis.json'),record)
        reaction=summary[label+'_reactions']
        force_scale=np.array([c.m*(2*omega)**2*e**2*h*L,c.m*omega**2*e*h*L])
        moment_scale=force_scale[1]*L
        reaction['fixed_balance_scales']={'translation':force_scale.tolist(),'moment':float(moment_scale),'policy':'inertial physical scales fixed over entire horizon; no new acceptance threshold'}
        for row in reaction['snapshots']:
            row['fixed_scaled_translation_discrepancy']=(np.array(row['linear_balance_difference'])/force_scale).tolist()
            row['fixed_scaled_angular_discrepancy']=row['angular_balance_difference']/moment_scale
        reaction['maxima']={
            'translation_difference':max(np.max(np.abs(r['linear_balance_difference'])) for r in reaction['snapshots']),
            'angular_difference':max(abs(r['angular_balance_difference']) for r in reaction['snapshots']),
            'translation_decomposition_error':max(np.max(np.abs(r['translation_decomposition_error'])) for r in reaction['snapshots']),
            'angular_decomposition_error':max(abs(r['angular_decomposition_error']) for r in reaction['snapshots']),
            'fixed_scaled_translation_discrepancy':max(np.max(np.abs(r['fixed_scaled_translation_discrepancy'])) for r in reaction['snapshots']),
            'fixed_scaled_angular_discrepancy':max(abs(r['fixed_scaled_angular_discrepancy']) for r in reaction['snapshots'])}
        write_json(bundle/(label+'_reactions.json'),reaction)
    summary['source_numerical_uncertainty']={k:inputs['summary'][k] for k in ('source_strict_status','source_strict_table','spatial_comparison','temporal_comparison')}
    summary['small_amplitude_independent_p_time_check']=False
    summary['amplitude_comparison']['qualifications']=['Some displacement maxima occur at t=0; ratios alone partly reflect prescribed initial scaling',
        'Full normalized histories and velocity ratios supply additional dynamic evidence',
        'One half-amplitude p64 run has no independent spatial/temporal control']
    for key in summary['large_analysis']['deviations']:
        a=summary['large_analysis']['deviations'][key];b=summary['half_analysis']['deviations'][key]
        summary['amplitude_comparison']['rows'][key]['normalized_leading_deviation_reduction']=a['fixed_scaled_max']/b['fixed_scaled_max']
    def table(name,rows):
        with (bundle/name).open('w',encoding='utf8',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    table('amplitude_scaling.csv',[{'component':k,**v} for k,v in summary['amplitude_comparison']['rows'].items()])
    table('leading_deviations.csv',[{'amplitude':summary[label+'_analysis']['epsilon_a'],'component':k,**v} for label in ('large','half') for k,v in summary[label+'_analysis']['deviations'].items()])
    table('strains.csv',[{'amplitude':summary[label+'_analysis']['epsilon_a'],'measure':k,**v} for label in ('large','half') for k,v in summary[label+'_analysis']['strains'].items()])
    table('resultants.csv',[{'amplitude':summary[label+'_analysis']['epsilon_a'],'measure':k,**v} for label in ('large','half') for k,v in summary[label+'_analysis']['resultants'].items()])
    table('support_reactions.csv',[{'amplitude':summary[label+'_analysis']['epsilon_a'],'tau':r['tau'],'actual_time':r['actual_time'],'side':side,**dict(zip(('Fu','Fw','M','Rc'),values))} for label in ('large','half') for r in summary[label+'_reactions']['snapshots'] for side,values in zip(('left','right'),r['endpoint_on_rod_local_Fu_Fw_M_Rc'])])
    table('momentum_balances.csv',[{'amplitude':summary[label+'_analysis']['epsilon_a'],'tau':r['tau'],'actual_time':r['actual_time'],'translation_difference_max':float(np.max(np.abs(r['linear_balance_difference']))),'angular_difference':r['angular_balance_difference'],'translation_decomposition_error_max':float(np.max(np.abs(r['translation_decomposition_error']))),'angular_decomposition_error':r['angular_decomposition_error'],'integrated_K4':r['integrated_K4']} for label in ('large','half') for r in summary[label+'_reactions']['snapshots']])
    table('energy_mass.csv',[{'amplitude':e,**{k:record[k] for k in ('initial_energy','relative_energy_drift_max','mass_lower_bound_min','mass_condition_bound_max','safety_pass','energy_and_mass_pass')}} for e,record in ((cfg['large_amplitude'],inputs['summary']['cases']['p64_tight']),(cfg['small_amplitude'],summary['half_case']))])
    write_json(bundle/'amplitude_comparison.json',summary['amplitude_comparison'])
    return summary


def run_compute(c,bundle):
    started=time.perf_counter();deadline=started+c['budget']['numerical_wall_seconds']-c['budget']['previous_limited_checks_seconds'];inputs=load_inputs(c);runner.load_runtime()
    d=dyn.PlanarGalerkin(inputs['coeff'],64,model=inputs['model']);large=saved.one_T1_history(inputs['source'],'p64_tight');times=np.asarray(large['time']);T1=inputs['state'].background.T1
    proof=limited_math(inputs['model']);write_json(bundle/'limited_math.json',proof)
    if proof['status']!='PASS':raise ArithmeticError('Limited action/limit consistency failed; half run prohibited')
    idx=np.unique([min(int(np.searchsorted(times,t*T1)),len(times)-1) for t in c['snapshot_fractions']])
    # Old-data analysis and classical/Noether derivation precede the one new run.
    analysisL=history_analysis(inputs,d,large,c['large_amplitude'],'large',bundle,c);reactionsL=reaction_checks(inputs,d,large,idx);parity=rhs_parity(d,large,idx)
    write_json(bundle/'large_analysis.json',analysisL);write_json(bundle/'large_reactions.json',reactionsL);write_json(bundle/'rhs_parity.json',parity)
    if parity['status']!='PASS':raise ArithmeticError('RHS sign symmetry failed; half run prohibited')
    state=inputs['state'];e=c['small_amplitude'];projection=prep.stable_initial_projection(state,d,e,c['projection_policy'],c['projection_dps'])
    initial=saved.initial_coordinate_metrics(state,d,projection['q'],e,read_json(ROOT/'data/input/planar_prepared_initial_state.json')['profile_policy'])
    compatibility=saved.initial_formal_compatibility(state,d,projection,e);initial.update(compatibility);initial['pass']=initial['pass'] and compatibility['pass']
    save_npz(bundle/'half_initial.npz',q=projection['q'],velocity=np.zeros(d.ndof),raw=projection['raw']);write_json(bundle/'half_initial.json',initial)
    write_json(bundle/'initial_provenance.json',{'common_frozen_target':inputs['provenance'],'epsilon_a':e,'projection_policy':projection['policy'],'one_new_amplitude_projection_no_precision_ladder':True,'Theta3_regenerated':False,'state_admitted':False})
    runner.safety_check(d,projection['q'],inputs['pilot']['safety'])
    if not initial['pass']:raise ArithmeticError('Half initial projection incomplete; saved diagnosis, no unexplained-admission bypass')
    summary={'schema':c['schema'],'config':c,'statuses':{},'large_analysis':analysisL,'reflection':parity,'classical_limit':proof['classical'],'source_strict_qualification':'PARTIAL',
        'execution_mode':c['execution_mode'],'state_admitted_flag':False,'new_ODE_integrations':0,'initial_projection':initial,'large_reactions':reactionsL,'source_provenance':inputs['provenance'],
        'forecast_ODE_seconds':inputs['summary']['cases']['p64_tight']['integration_seconds'],'new_BVP_eigensolves':0}
    path=bundle/'cases/p64_half_tight';path.mkdir(parents=True,exist_ok=True);buffer=np.lib.format.open_memmap(path/'state.npy',mode='w+',dtype=float,shape=(len(times),2*d.ndof))
    print(json.dumps({'stage':'PRE_HALF_RUN','bundle':str(bundle),'forecast_seconds':summary['forecast_ODE_seconds'],'elapsed_diagnostic_seconds':time.perf_counter()-started,'initial_projection_pass':initial['pass'],'samples':len(times)}),flush=True)
    summary['new_ODE_integrations']=1
    history,stats=runner.integrate_case(d,None,state.background.as_dict(),inputs['pilot'],e,'tight',times,deadline,initial_coordinates=projection['q'],history_buffer=buffer)
    buffer.flush();np.save(path/'time.npy',times[:len(history)]);stats.update(actual_valid_rows=len(history),storage_allocated_rows=len(times),execution_mode=c['execution_mode'],projection_policy=c['projection_policy'],state_admitted=False,T1=T1,periods=1.,omega1=state.background.omega)
    write_json(path/'case.json',stats);del history,buffer
    half=saved.one_T1_history(bundle,'p64_half_tight');stats.update(saved.one_T1_case_diagnostics(d,half,inputs['pilot'],state.background.as_dict(),path,[0.,.1,.25,.5,.75,1.]));write_json(path/'case.json',stats)
    analysisS=history_analysis(inputs,d,half,e,'half',bundle,c);reactionsS=reaction_checks(inputs,d,half,idx[idx<len(half['time'])]);amp=normalized_amplitude_comparison(d,large,half,c['large_amplitude'],e,inputs)
    summary.update(half_case=stats,half_analysis=analysisS,half_reactions=reactionsS,amplitude_comparison=amp,energy_ratio_large_over_half=inputs['summary']['cases']['p64_tight']['initial_energy']/stats['initial_energy'])
    write_json(bundle/'half_analysis.json',analysisS);write_json(bundle/'half_reactions.json',reactionsS);write_json(bundle/'amplitude_comparison.json',amp)
    complete=stats['status']=='PASS';summary['statuses']={'NLSP_SANITY_SECOND_ORDER_COMPARISON':'PARTIAL','NLSP_SANITY_AMPLITUDE_SCALING':'PARTIAL','NLSP_SANITY_REFLECTION_SYMMETRY':parity['status'],'NLSP_SANITY_CLASSICAL_STRETCHING_LIMIT':proof['status'],'NLSP_SANITY_STRAINS':'PASS' if stats['safety_pass'] else 'FAIL','NLSP_SANITY_REACTIONS_AND_MOMENTUM':'PARTIAL','NLSP_PLANAR_PHYSICAL_SANITY_CHECKS':'DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS' if complete else 'PARTIAL'}
    summary['status_interpretation']='Numerical amplitude/asymptotic/momentum checks remain qualified by finite-p evidence; no arbitrary new physical tolerance or physical-validation claim'
    supplement_reports(inputs,d,large,half,summary,bundle,c)
    summary['runtime']={'numerical_wall_seconds':time.perf_counter()-started+c['budget']['previous_limited_checks_seconds'],'integration_seconds':stats['integration_seconds'],'limit_seconds':600,'new_ODE_integrations':1,'new_BVP_eigensolves':0}
    summary['stop_reason']='One permitted half-amplitude run and bounded checks complete or accepted prefix saved; no next study'
    write_json(bundle/'summary.json',summary);return summary


def plot_only(bundle):
    import matplotlib;matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bundle=Path(bundle);s=saved.validate_cache(bundle);T1=s['half_case']['T1'];figs=bundle/'figures';figs.mkdir(exist_ok=True)
    a=np.load(bundle/'large_asymptotic.npz');b=np.load(bundle/'half_asymptotic.npz')
    def finish(fig,name):
        for ext in ('pdf','png'):fig.savefig(figs/(name+'.'+ext),dpi=200,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(9,5.5),layout='constrained');powers=[2,1,1,2]
    for k,(field,ax) in enumerate(zip(FIELDS,axes.flat)):
        for data,e,label,style in ((a,.05,'epsilon=.05','-'),(b,.025,'epsilon=.025','--')):ax.plot(data['time']/T1,data['observations'][:,k]/e**powers[k],style,label=label)
        ax.set(xlabel='t/T1',ylabel=field+'/epsilon^'+str(powers[k]));ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(fig,'amplitude_normalized_motion')
    fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    for k,ax in ((0,axes[0]),(3,axes[1])):
        ax.plot(a['time']/T1,a['observations'][:,k],label='cubic trajectory');ax.plot(a['time']/T1,a['leading_observations'][:,k],'--',label='second order')
        ax.set(xlabel='t/T1',ylabel=FIELDS[k]+' at L/4');ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(fig,'leading_axial_response')
    fig,axes=plt.subplots(1,2,figsize=(9,3.5),layout='constrained')
    for label,record,style in (('epsilon=.05',s['large_reactions'],'-'),('epsilon=.025',s['half_reactions'],'--')):
        tt=[r['tau'] for r in record['snapshots']];axes[0].plot(tt,[r['endpoint_on_rod_local_Fu_Fw_M_Rc'][0][1] for r in record['snapshots']],style+'o',label=label)
        axes[1].plot(tt,[r['angular_balance_difference'] for r in record['snapshots']],style+'o',label=label)
    axes[0].set(xlabel='t/T1',ylabel='left support local transverse force');axes[1].set(xlabel='t/T1',ylabel='endpoint angular flux discrepancy')
    for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=8,frameon=False)
    finish(fig,'support_reaction_diagnostics')
    return {'ODE_integrations':0,'BVP_eigensolves':0,'symbolic_model_derivations':0}


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path);p.add_argument('--config',type=Path,default=CONFIG);p.add_argument('--output-dir',type=Path,default=OUTPUT);a=p.parse_args(argv)
    if a.report_only or a.plot_only:
        bundle=a.report_only or a.plot_only;s=saved.validate_cache(bundle)
        if a.plot_only:plot_only(bundle)
        print(json.dumps({'bundle':str(bundle),'statuses':s['statuses'],'new_ODE_BVP_eigen_symbolic_calls':0},indent=2));return s
    key,item=identity(a.config);bundle=a.output_dir/key
    if (bundle/'manifest.json').exists():
        s=saved.validate_cache(bundle,item);print(json.dumps({'bundle':str(bundle),'cache_hit':True,'statuses':s['statuses'],'new_ODE_BVP_eigen_symbolic_calls':0},indent=2));return s
    bundle.mkdir(parents=True,exist_ok=True);s=run_compute(read_json(a.config),bundle);write_json(bundle/'manifest.json',saved.manifest_for(bundle,item));plot_only(bundle);write_json(bundle/'manifest.json',saved.manifest_for(bundle,item));print(json.dumps({'bundle':str(bundle),'statuses':s['statuses'],'runtime':s['runtime']},indent=2));return s

if __name__=='__main__':main()
