"""Pointwise modal/harmonic diagnostics; this module has no IVP entry path.

Reuses the full frozen coupled Jacobian. Frequencies, metric and both signs
of each state-space mode belong to the same p12 Galerkin problem.
"""
from __future__ import annotations
import argparse,csv,hashlib,json,shutil,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import eigh,cho_factor
from scipy.optimize import linear_sum_assignment
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_construction_c_sensitivity as c
from scripts.lib import nlsp_construction_c_longer_sensitivity as l
from scripts.lib import nlsp_coupled_sensitivity_resume as r

CONFIG=p.ROOT/'data/input/nlsp_resonant_spatial_coupling_screening.json'
OUTPUT=p.ROOT/'results/nlsp_resonant_spatial_coupling_screening'
SMOKE=p.ROOT/'results/_smoke/nlsp_resonant_spatial_coupling_screening'
AUTH='explicit_user_resonant_spatial_coupling_screening_2026_10_11'
BETA_CONFIG=p.ROOT/'data/input/nlsp_beta_resonant_coupling_screening.json'
BETA_OUTPUT=p.ROOT/'results/nlsp_beta_resonant_coupling_screening'
BETA_SMOKE=p.ROOT/'results/_smoke/nlsp_beta_resonant_coupling_screening'
BETA_AUTH='explicit_user_beta_resonant_coupling_screening_2026_10_11'

def validate_scope(config):
    if config!=p.read(CONFIG) or config['authorization']!=AUTH:raise ValueError('Only the frozen screening is authorized')
    if (config['degree']!=12 or config['amplitude']!=.005 or config['mu_points']!=[0.,.1,.2,.3,.4]
        or config['numerical_budget_seconds']>900 or config['new_ODE_integrations'] or config['new_FEM_jobs']
        or config['automatic_retry'] or config['out_modes_retained']!=5 or config['phase_samples']>64):raise ValueError('Scope expansion')
    return config

def verify_sources(config):
    for name in ('modal_parent','longer_parent'):
        path=p.ROOT/config[name]
        if p.sha(path/'manifest.json')!=config[name+'_manifest_sha256']:raise ValueError('Historical source manifest mismatch')
    l.cache(p.ROOT/config['longer_parent'])
    parent=p.read(p.ROOT/config['modal_parent']/'provenance.json')
    if any(p.sha(p.ROOT/f)!=parent['code_hashes'][f] for f in r.KERNELS):raise ValueError('Frozen physical kernel changed')

def prepare(path=CONFIG):
    config=validate_scope(p.read(path));verify_sources(config)
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        item=p.read(old)
        if item['config']['authorization']==AUTH:
            if item['config']!=config:raise ValueError('Authorization already used with different settings')
            return old.parent,item,p.cache(old.parent)
    parent=p.read(p.ROOT/config['modal_parent']/'provenance.json')
    hashes={**p.code_hashes(),Path(__file__).relative_to(p.ROOT).as_posix():p.sha(__file__),CONFIG.relative_to(p.ROOT).as_posix():p.sha(CONFIG)}
    theories=['docs/theory/weakly_nonlinear_spatial_rod.md','docs/theory/weakly_nonlinear_spatial_rod_expansion_generated.md']
    item={'schema':'nlsp-resonant-spatial-coupling-screening-v1','config':config,'physical_config':p.read(c.CONFIG),
        'physical_source_config':parent['physical_source_config'],'source_artifacts':{**parent['source_artifacts'],**{f:p.sha(p.ROOT/f) for f in theories}},
        'source_manifests':{**parent['source_manifests'],**{key:{'path':config[key],'sha256':config[key+'_manifest_sha256']} for key in ('modal_parent','longer_parent')}},
        'code_hashes':hashes,'initial_state':p.read(SMOKE/'initial_state.json'),'interpreter':sys.executable,
        'no_new_ODE_or_Floquet':True,'modal_projection':'right E0-unit complex frequency modes; left rows R.H E0; both signs retained'}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp
    if b.exists():raise ValueError('Preserve interrupted attempt')
    b.mkdir(parents=True);p.write(b/'provenance.json',item);p.write(b/'config.json',config)
    folder=b/'execution_code';folder.mkdir()
    for f in hashes:shutil.copyfile(p.ROOT/f,folder/Path(f).name)
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':5.,'preflight_charge_qualification':'5-second conservative reservation for initial p12 low-eigenpair probe; includes no ODE',
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'CCX':0,'Gmsh':0,'linear_eigenproblems':2,'pointwise_Jacobian':0},'stages':{'source':'PASS','A':'NOT_RUN','B':'NOT_RUN','C':'NOT_RUN'}}
    p.finalize(b,item,s);return b,item,s

def state_modes(K,M,omega,phi):
    """R+ = (phi/omega,i phi)/sqrt(2), R-=conj(R+); R.H E R=I."""
    n,k=phi.shape;E=np.zeros((2*n,2*n));E[:n,:n]=K;E[n:,n:]=M
    Q=phi/(np.sqrt(2)*omega[None,:]);V=1j*phi/np.sqrt(2)
    R=np.block([[Q,Q],[V,-V]]);L=R.conj().T@E
    return R,L,E

def modal_problem(model,mu,config,counter,*,beta_deg=None):
    pc={**p.read(c.CONFIG),'lengths':[1-mu,1+mu],'planar_modes_inspected':config['planar_modes_inspected']}
    if beta_deg is not None:pc['beta_deg']=beta_deg
    d=p.make_disc(pc,'C',config['degree'],model)
    q,modal=c.select_modal_state(d,pc);counter['linear_eigenproblems']+=1
    omega,phi=c.limited_modes(d,False,config['out_modes_retained']+config['out_modes_guard']);counter['linear_eigenproblems']+=1
    # Rayleigh-Ritz inside these already computed low vectors improves their
    # K/M duality without searching any additional modes or changing the space.
    K=d.K[np.ix_(d.out,d.out)];M=d.M0[np.ix_(d.out,d.out)]
    val,rotation=eigh(phi.T@K@phi,phi.T@M@phi,check_finite=False)
    phi=phi@rotation;omega=np.sqrt(val)
    profiles=[];torsion=[]
    for j in range(len(omega)):
        x=np.zeros(2*d.ndof);x[d.out]=phi[:,j]
        fields,_=l.a.physical_observations(d,x[None,:]);field=fields[0]
        sign=np.sign(field[...,2:5].reshape(-1)[np.argmax(abs(field[...,2:5]))])
        if sign<0:phi[:,j]*=-1;field=-field
        profiles.append(field)
        v=np.zeros(d.ndof);v[d.out]=phi[:,j];parts=d.reference_energy_parts(np.zeros(d.ndof),v)
        torsion.append(float(parts[:,1].sum()/parts.sum()))
    k=config['out_modes_retained'];guard_omega=omega[k:].copy();phi=phi[:,:k];omega=omega[:k]
    R,L,E=state_modes(K,M,omega,phi);duality=float(np.max(abs(L@R-np.eye(2*k))))
    if duality>config['modal_duality_gate']:raise ArithmeticError('Modal energy duality failed')
    cho_factor(K,check_finite=False);cho_factor(M,check_finite=False)
    A0=d.jacobian(0,np.zeros(2*d.ndof))[np.ix_(d.state_out,d.state_out)];counter['pointwise_Jacobian']+=1
    eigen=np.r_[1j*omega,-1j*omega];residual=float(np.linalg.norm(A0@R-R*eigen)/np.linalg.norm(R*eigen))
    if residual>2e-8:raise ArithmeticError('Zero-base state eigenpair residual')
    return {'d':d,'q':q,'Omega':modal['omega'],'modal':modal,'omega':omega,'phi':phi,'R':R,'L':L,'E':E,'A0':A0,
        'profiles':np.array(profiles),'torsional_mass_fraction':torsion,'guard_omega':guard_omega,'mu':mu,'duality':duality,'eigen_residual':residual,'beta_deg':d.beta_deg}

def harmonics(values,phases,orders):
    return np.array([np.einsum('t,tij->ij',np.exp(-1j*m*phases),values)/len(phases) for m in orders])

def projected_at(problem,state,counter):
    d=problem['d'];J=d.jacobian(0,state);counter['pointwise_Jacobian']+=1
    A=J[np.ix_(d.state_out,d.state_out)]
    return problem['L']@(A-problem['A0'])@problem['R']

def harmonic_sample(problem,factor,count,orders,counter,deadline):
    phases=np.arange(count)*2*np.pi/count;values=[];d=problem['d'];q=problem['q']*factor
    for phase in phases:
        if time.perf_counter()>deadline:raise TimeoutError('Screening numerical ceiling')
        state=np.r_[q*np.cos(phase),-problem['Omega']*q*np.sin(phase)]
        values.append(projected_at(problem,state,counter))
    return harmonics(np.array(values),phases,orders),np.array(values)

def amplitude_orders(H,amplitudes,orders):
    """Leading linear odd and quadratic even Fourier terms, finite-amplitude diagnostic."""
    powers=np.array([2 if m%2==0 else 1 for m in orders])
    scaled=H/(amplitudes[:,None,None,None]**powers[None,:,None,None])
    leading=(4*scaled[1]-scaled[0])/3;next_leading=(4*scaled[2]-scaled[1])/3
    return leading,next_leading,powers

def complex_row(value):return {'real':float(value.real),'imag':float(value.imag),'abs':float(abs(value))}

def pair_rows(problem,H,uncertainty,config):
    omega=problem['omega'];Omega=problem['Omega'];k=len(omega);orders=config['harmonics'];rows=[]
    for kind in ('self','sum','difference'):
        pairs=[(j,j) for j in range(k)] if kind=='self' else [(i,j) for i in range(k) for j in range(i+1,k)]
        for i,j in pairs:
            for m in (1,2):
                if kind=='difference':ridx,sidx=j,i;target=omega[j]-omega[i]
                else:ridx,sidx=i,k+j;target=omega[i]+omega[j]
                hm=H[orders.index(m),ridx,sidx];back=H[orders.index(-m),sidx,ridx]
                product=hm*back;coupling=float(np.sqrt(abs(product)));delta=float(m*Omega-target)
                hr=H[orders.index(0),ridx,ridx];hs=H[orders.index(0),sidx,sidx]
                shifted=delta-float((hr-hs).imag)
                small=np.array([[hr-1j*delta/2,hm],[back,hs+1j*delta/2]])
                growth=float(np.max(np.linalg.eigvals(small).real))
                error=max(float(uncertainty[orders.index(m),ridx,sidx]),float(uncertainty[orders.index(-m),sidx,ridx]),1e-13)
                resolved=bool(min(abs(hm),abs(back))>config['coupling_over_numerical_difference_gate']*error)
                growth_candidate=bool(kind!='difference' and resolved and product.real>0 and abs(product.imag)<=.05*abs(product)
                    and 2*coupling>abs(shifted) and growth>10*error)
                rows.append({'kind':kind,'i_sorted':i+1,'j_sorted':j+1,'harmonic':m,'mu':problem['mu'],
                    'Omega':float(Omega),'omega_i':float(omega[i]),'omega_j':float(omega[j]),'target_frequency':float(target),
                    'signed_detuning':delta,'absolute_detuning':abs(delta),'detuning_over_pump':delta/(m*Omega),
                    'coupling_forward':complex_row(hm),'coupling_reverse':complex_row(back),'product':complex_row(product),
                    'sqrt_abs_product':coupling,'relative_coefficient':float(abs(hm)/np.sqrt(omega[i]*omega[j])),
                    'coupling_numerical_difference':error,'coupling_resolved':resolved,
                    'mean_shifted_detuning':shifted,'two_coupling_over_bare_detuning':None if delta==0 else 2*coupling/abs(delta),
                    'two_coupling_over_mean_shifted_detuning':None if shifted==0 else 2*coupling/abs(shifted),
                    'isolated_pair_averaging_growth_proxy':growth,'direct_term_candidate':growth_candidate,
                    'classification':'difference_exchange' if kind=='difference' and resolved and product.real<0 else 'direct_parametric_candidate' if growth_candidate else 'coupling_below_numerical_indicator' if not resolved else 'off_resonance_direct_term',
                    'qualification':'isolated resonant-term algebraic proxy, not Floquet, instability threshold, or complete second-order normal form'})
    return rows

def save_problem(b,label,problem):
    np.savez_compressed(b/(label+'_modes.npz'),omega=problem['omega'],phi=problem['phi'],right=problem['R'],left=problem['L'],metric=problem['E'],
        planar_initial_shape=problem['q'],OUT_physical_profiles=problem['profiles'],x=np.linspace(0,1,201),M0=problem['d'].M0,K0=problem['d'].K)
    p.write(b/(label+'_mode_metadata.json'),{'mu':problem['mu'],'lengths':problem['d'].lengths,'degree':problem['d'].p,'Omega':problem['Omega'],
        'beta_deg':problem['d'].beta_deg,
        'planar':problem['modal'],'omega':problem['omega'],'guard_omega':problem['guard_omega'],'torsional_mass_fraction':problem['torsional_mass_fraction'],
        'state_duality_max_error':problem['duality'],'zero_operator_eigen_residual':problem['eigen_residual'],
        'named_historical_modes':'historical omega1,omega2,omega3 correspond to OUT sorted1,2,5; p16 references are not projected here'})

def stage_b(b,item,s,problem,deadline):
    config=item['config'];orders=config['harmonics'];count=config['phase_samples'];start=time.perf_counter()
    H=[];coarse=[]
    for factor in config['amplitude_factors']:
        fine,values=harmonic_sample(problem,factor,count,orders,s['calls'],deadline)
        H.append(fine);coarse.append(harmonics(values[::2],np.arange(count//2)*2*np.pi/(count//2),orders))
    H=np.array(H);coarse=np.array(coarse);amplitudes=config['amplitude']*np.array(config['amplitude_factors'])
    leading,next_leading,powers=amplitude_orders(H,amplitudes,orders)
    actual_leading=leading*amplitudes[0]**powers[:,None,None]
    amplitude_error=abs((leading-next_leading)*amplitudes[0]**powers[:,None,None]);phase_error=abs(H[0]-coarse[0])
    floor=max(problem['duality']*problem['Omega']*1e-3,1e-13)
    uncertainty=np.maximum(np.maximum(amplitude_error,phase_error),floor)
    rows=pair_rows(problem,actual_leading,uncertainty,config)
    symmetry=[]
    for phase in (0.,.3,1.4,2.9):
        q=problem['q'];Omega=problem['Omega']
        minus=projected_at(problem,np.r_[-q*np.cos(phase),Omega*q*np.sin(phase)],s['calls'])
        plus=projected_at(problem,np.r_[q*np.cos(phase+np.pi),-Omega*q*np.sin(phase+np.pi)],s['calls'])
        symmetry.append(float(np.max(abs(minus-plus))))
    k=len(problem['omega']);swap=np.r_[np.arange(k,2*k),np.arange(k)]
    conjugacy=max(float(np.max(abs(H[0,orders.index(-m)]-H[0,orders.index(m)][np.ix_(swap,swap)].conj()))) for m in (0,1,2))
    norms=np.linalg.norm(H,axis=(-2,-1));ratio_order=np.log2(norms[0]/np.maximum(norms[1],1e-300))
    diagnostic={'amplitudes':amplitudes,'harmonics':orders,'matrix_norms':norms,'observed_full_half_orders':ratio_order,
        'phase_grid_absolute_difference':float(np.max(phase_error)),'Richardson_coefficient_difference_at_full_amplitude':float(np.max(amplitude_error)),
        'negative_amplitude_phase_shift_max_difference':max(symmetry),'conjugate_branch_Fourier_difference':conjugacy,
        'leading_order_for_each_harmonic':powers,'exact_mass_inverse_higher_terms_not_new_action_coefficients':True,
        'quartic_order_control':'finite amplitude plus half/quarter and Richardson; no two-parameter expansion of original action',
        'direct_coefficient_not_complete_second_order_normal_form':True}
    np.savez_compressed(b/'C_harmonic_coefficients.npz',amplitudes=amplitudes,orders=orders,full=H,coarse=coarse,leading_per_amplitude=leading,
        next_leading_per_amplitude=next_leading,leading_at_full_amplitude=actual_leading,numerical_difference=uncertainty)
    p.write(b/'C_amplitude_order_diagnostics.json',diagnostic);p.write(b/'C_interaction_rows.json',rows)
    resolved=[row for row in rows if row['coupling_resolved']]
    if not resolved or max(symmetry)>1e-9 or conjugacy>1e-9:raise ArithmeticError('Hard gate B: unresolved coefficient or symmetry')
    s['stages']['B']='PASS';s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)
    return H[0],rows

def actual_trace(b,item,s,problem,harmonic_reference,deadline):
    started=time.perf_counter();config=item['config']
    with np.load(p.ROOT/config['longer_parent']/'p12_merged.npz') as z:times=z['times'].copy();states=z['state'].copy()
    results=[];arrays={}
    for count in (config['actual_trajectory_samples'],config['actual_trajectory_fine_samples']):
        indices=np.unique(np.rint(np.linspace(0,len(times)-1,count)).astype(int));t=times[indices];values=[]
        for index in indices:
            if time.perf_counter()>deadline:raise TimeoutError('Saved-trajectory evaluation budget')
            values.append(projected_at(problem,states[index],s['calls']))
        values=np.array(values);design=np.exp(1j*t[:,None]*problem['Omega']*np.array(config['harmonics'])[None,:])
        coeff=np.linalg.lstsq(design,values.reshape(len(t),-1),rcond=None)[0].reshape(len(config['harmonics']),values.shape[1],values.shape[2])
        fitted=np.einsum('tm,mij->tij',design,coeff);norm=float(np.linalg.norm(values))
        results.append({'samples':len(t),'actual_start':float(t[0]),'actual_end':float(t[-1]),'time_grid_max_spacing':float(np.diff(t).max()),
            'fit_residual_over_total_operator_norm':float(np.linalg.norm(values-fitted)/norm),
            'harmonic_matrix_norms':np.linalg.norm(coeff,axis=(-2,-1)),
            'coefficient_difference_over_reference_Omega':np.linalg.norm(coeff-harmonic_reference,axis=(-2,-1))/problem['Omega']})
        arrays[str(count)+'_times']=t;arrays[str(count)+'_values']=values;arrays[str(count)+'_fit']=coeff
    fine=arrays[str(config['actual_trajectory_fine_samples'])+'_fit'];coarse=arrays[str(config['actual_trajectory_samples'])+'_fit']
    np.savez_compressed(b/'C_actual_trace_projection.npz',**arrays,harmonic_ansatz_reference=harmonic_reference,orders=config['harmonics'])
    p.write(b/'C_actual_trace_diagnostics.json',{'finite_window_not_Floquet':True,'not_a_periodic_orbit_Fourier_series':True,
        'Omega_same_p12_linear_value':problem['Omega'],'sampled_velocity_states_not_derivatives':True,'fits':results,
        'coarse_fine_coefficient_absolute_difference':float(np.max(abs(fine-coarse))),
        'frequency_peaks_not_called_nonlinear_eigenfrequencies':True})
    s['numerical_seconds']+=time.perf_counter()-started;p.finalize(b,item,s)

def physical_MAC(first,second):
    """Common material coordinates and reference kinetic weights, not coefficients."""
    a=first['profiles'][:5];b=second['profiles'][:5];length=(np.array(first['d'].lengths)+np.array(second['d'].lengths))/2
    co=first['d'].arms[0].coefficients;weight=np.array([co.m,co.m,co.m,co.jp+co.jb,co.jb,co.jp,co.jp])
    quad=np.ones(201);quad[[0,-1]]=.5;W=length[:,None,None]*quad[None,:,None]*weight[None,None,:]/200
    A=(a*np.sqrt(W)).reshape(len(a),-1);B=(b*np.sqrt(W)).reshape(len(b),-1)
    return (A@B.T)**2/(np.sum(A*A,axis=1)[:,None]*np.sum(B*B,axis=1)[None,:])

def mirror_selection(problem,H,config,counter):
    """Actual equal-arm reflection: polar displacement and axial rotation."""
    if problem['mu']!=0:return {'status':'NOT_APPLICABLE: unequal arm lengths; no exact arm-exchange symmetry asserted'}
    d=problem['d'];n=d.p-1;raw=np.zeros((d.ndof,d.ndof));signs=[1,-1,1,-1,1,-1,1]
    for f,sign in enumerate(signs):
        left=np.arange(f*n,(f+1)*n);right=left+7*n
        raw[np.ix_(left,right)]=sign*np.eye(n);raw[np.ix_(right,left)]=sign*np.eye(n)
    bisector=d.bases[0][:,0]+d.bases[1][:,0];bisector/=np.linalg.norm(bisector)
    reflection=2*np.outer(bisector,bisector)-np.eye(3);reflection[2,2]=1
    joint=np.zeros((7,7));joint[:3,:3]=reflection;joint[3:6,3:6]=-reflection;joint[6,6]=1
    raw[np.ix_(d.joint_ids,d.joint_ids)]=joint;S=np.linalg.solve(d.W,raw@d.W)
    q=problem['q'];plpar=float(q@d.M0@S@q/(q@d.M0@q));par=[]
    for j in range(len(problem['omega'])):
        full=np.zeros(d.ndof);full[d.out]=problem['phi'][:,j];par.append(float(full@d.M0@S@full/(full@d.M0@full)))
    Sout=S[np.ix_(d.out,d.out)];Sstate=np.zeros_like(problem['A0']);N=len(Sout);Sstate[:N,:N]=Sout;Sstate[N:,N:]=Sout
    phase=.7;x=np.r_[q*np.cos(phase),-problem['Omega']*q*np.sin(phase)];fullS=np.zeros((2*d.ndof,2*d.ndof));fullS[:d.ndof,:d.ndof]=S;fullS[d.ndof:,d.ndof:]=S
    J=d.jacobian(0,x);transformed=d.jacobian(0,fullS@x);counter['pointwise_Jacobian']+=2
    equiv=float(np.linalg.norm(transformed-fullS@J@fullS)/np.linalg.norm(J))
    if equiv>1e-10 or min(abs(plpar),min(abs(np.array(par))))<.999: return {'status':'UNRESOLVED','operator_equivariance_relative':equiv,'planar_parity':plpar,'OUT_parities':par}
    forbidden={}
    for m in (1,2):
        expected=np.sign(plpar)**m;mask=np.outer(np.sign(par),np.sign(par))!=expected
        block=H[config['harmonics'].index(m),:len(par),:len(par)]
        forbidden[str(m)]={'forbidden_pairs_max_coefficient':float(np.max(abs(block[mask]))) if np.any(mask) else 0.,
            'allowed_pairs_max_coefficient':float(np.max(abs(block[~mask]))) if np.any(~mask) else 0.}
    return {'status':'PASS','operator_equivariance_relative':equiv,'planar_parity':plpar,'OUT_parities':par,'selection_rule':'parity_i*parity_j = parity_planar**m','coefficients':forbidden}

def stage_c(b,item,s,model,reference,deadline):
    started=time.perf_counter();config=item['config'];problems={.2:reference};results={};allrows=[];previous=reference
    for mu in (.1,0.,.3,.4):
        if mu==.3:previous=reference
        if time.perf_counter()>deadline:raise TimeoutError('Geometry screening budget')
        prob=modal_problem(model,mu,config,s['calls']);mac=physical_MAC(previous,prob);old_ids,new_ids=linear_sum_assignment(-mac)
        quality=float(min(mac[old_ids,new_ids]));correspondence={'reference_mu':previous['mu'],'MAC_matrix':mac,'matching_sorted_indices':new_ids+1,'minimum_assigned_MAC':quality,'status':'PASS' if quality>=config['correspondence_MAC_gate'] else 'UNRESOLVED',
            'identity':'diagnostic physical-shape correspondence; not canonical analytic branch_id'}
        save_problem(b,f'mu_{mu:.1f}',prob);H=[]
        for factor in (1.,.5):
            h,values=harmonic_sample(prob,factor,config['geometry_phase_samples'],config['harmonics'],s['calls'],deadline);H.append(h)
        H=np.array(H);powers=np.array([2 if m%2==0 else 1 for m in config['harmonics']]);a=config['amplitude']
        leading=(4*H[1]/(.5*a)**powers[:,None,None]-H[0]/a**powers[:,None,None])/3*a**powers[:,None,None]
        err=abs(H[0]-leading);rows=pair_rows(prob,leading,np.maximum(err,1e-13),config)
        for row in rows:row['modal_correspondence_status']=correspondence['status'];row['direct_term_candidate'] &= correspondence['status']=='PASS'
        selection=mirror_selection(prob,leading,config,s['calls'])
        if selection.get('status')=='UNRESOLVED':
            for row in rows:row['direct_term_candidate']=False
        np.savez_compressed(b/f'mu_{mu:.1f}_coefficients.npz',full=H,leading_at_full_amplitude=leading,numerical_difference=err,orders=config['harmonics'])
        p.write(b/f'mu_{mu:.1f}_correspondence.json',correspondence);p.write(b/f'mu_{mu:.1f}_selection.json',selection)
        results[str(mu)]={'Omega':prob['Omega'],'omega':prob['omega'],'correspondence':correspondence['status'],'minimum_MAC':quality,'selection':selection,'coupling_resolved_count':sum(row['coupling_resolved'] for row in rows)}
        allrows.extend(rows);problems[mu]=prob;previous=prob
    p.write(b/'mu_screening.json',results);p.write(b/'mu_interaction_rows.json',allrows)
    s['stages']['C']='PASS' if all(r['correspondence']=='PASS' for r in results.values()) else 'PARTIAL';s['numerical_seconds']+=time.perf_counter()-started;p.finalize(b,item,s)
    return allrows

def run(b,item,s,through='all'):
    if s['completed'] or s['hard_stop']:return s
    config=item['config'];start=time.perf_counter();deadline=start+config['numerical_budget_seconds']-s['numerical_seconds']
    model,_,_=p.sources(item['physical_source_config']);problem=modal_problem(model,.2,config,s['calls']);save_problem(b,'C',problem)
    with np.load(p.ROOT/config['modal_parent']/'p12_projected_initial.npz') as z:sourceq=z['q']
    q=problem['q'];d=problem['d'];mac=float((sourceq@d.M0@q)**2/(sourceq@d.M0@sourceq)/(q@d.M0@q))
    if mac<config['source_mode_MAC_gate']:raise ArithmeticError('Source planar mode mismatch')
    np.testing.assert_allclose(d.K[np.ix_(d.out,d.out)],np.load(p.ROOT/config['modal_parent']/'p12_projected_initial.npz')['metric'][:len(d.out),:len(d.out)],rtol=2e-13,atol=2e-13)
    p.write(b/'C_frequency_source_verification.json',{'source_projected_planar_MAC':mac,'Omega_p12':problem['Omega'],'omega_p12':problem['omega'],
        'known_p16_Omega':p.read(p.ROOT/config['modal_parent']/'frozen_modal_excitation.json')['omega'],
        'known_p16_selected_OUT':p.read(p.ROOT/config['modal_parent']/'restricted_subspace.json')['omega'],
        'no_mixing_p12_p16_in_operator':True,'Omega_minus_4omega1':float(problem['Omega']-4*problem['omega'][0]),
        'Omega_minus_4omega1_over_Omega':float((problem['Omega']-4*problem['omega'][0])/problem['Omega']),
        '2Omega_minus_sorted5_minus_sorted2':float(2*problem['Omega']-(problem['omega'][4]-problem['omega'][1]))})
    s['stages']['A']='PASS';s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)
    H,rows=stage_b(b,item,s,problem,deadline);actual_trace(b,item,s,problem,H,deadline)
    if through=='C-only':s.update(overall='COUPLING_PRESENT_NO_GROWTH_MECHANISM_IDENTIFIED',completed=True);p.finalize(b,item,s);return s
    if time.perf_counter()+20<deadline:
        geometry_rows=stage_c(b,item,s,model,problem,deadline);candidates=[row for row in rows+geometry_rows if row['direct_term_candidate']]
    else:geometry_rows=[];candidates=[];s['stages']['C']='NOT_RUN'
    status='RESONANT_COUPLING_CANDIDATE_IDENTIFIED' if candidates else 'NUMERICAL_PARTIAL' if s['stages']['C']=='PARTIAL' else 'COUPLING_PRESENT_NO_GROWTH_MECHANISM_IDENTIFIED'
    s.update(overall=status,completed=True,candidate_count=len(candidates));p.write(b/'candidate_rows.json',candidates)
    export_saved(b);plot_saved(b);p.finalize(b,item,s);return s

def cache(b):
    config=p.read(Path(b)/'config.json');verify_sources(config);return p.cache(Path(b))

def export_saved(b):
    rows=p.read(b/'C_interaction_rows.json')
    if (b/'mu_interaction_rows.json').exists():rows+=p.read(b/'mu_interaction_rows.json')
    fieldnames=['mu','kind','i_sorted','j_sorted','harmonic','Omega','omega_i','omega_j','signed_detuning','detuning_over_pump',
        'forward_real','forward_imag','forward_abs','reverse_real','reverse_imag','product_real','product_imag','sqrt_abs_product',
        'relative_coefficient','mean_shifted_detuning','two_coupling_over_bare_detuning','two_coupling_over_mean_shifted_detuning',
        'coupling_resolved','direct_term_candidate','classification']
    with (b/'interaction_coefficients.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=fieldnames);writer.writeheader()
        for row in sorted(rows,key=lambda r:(r['mu'],r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])):
            values={key:row.get(key) for key in fieldnames}
            for name in ('real','imag','abs'):values['forward_'+name]=row['coupling_forward'][name]
            for name in ('real','imag'):values['reverse_'+name]=row['coupling_reverse'][name];values['product_'+name]=row['product'][name]
            writer.writerow(values)
    p.write(b/'export_provenance.json',{'new_ODE_calls':0,'data_sources':['C_interaction_rows.json','mu_interaction_rows.json if available'],
        'complex_values_preserved_in_NPZ_and_real_imag_CSV':True})

def plot_saved(b):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    config=p.read(b/'config.json');plots=[]
    def finish(fig,name):
        fig.tight_layout()
        for ext in ('pdf','png'):fig.savefig(b/(name+'.'+ext),dpi=180,bbox_inches='tight');plots.append(name+'.'+ext)
        plt.close(fig)
    with np.load(b/'C_harmonic_coefficients.npz') as z:H=z['leading_at_full_amplitude'];orders=z['orders'].tolist()
    omega=np.array(p.read(b/'C_mode_metadata.json')['omega']);k=len(omega)
    fig,axes=plt.subplots(2,2,figsize=(10,8));maximum=0.
    for m in (1,2):
        maximum=max(maximum,float(np.max(abs(H[orders.index(m),:k]))))
    for row,m in enumerate((1,2)):
        for col,title in enumerate(('same-sign branches (+ to +)','opposite-sign branches (- to +)')):
            matrix=abs(H[orders.index(m),:k,col*k:(col+1)*k]);ax=axes[row,col]
            im=ax.imshow(matrix,vmin=0,vmax=maximum,cmap='viridis',origin='lower')
            ax.set(title=f'{m} Omega: '+('+ / +' if col==0 else '+ / -'),xlabel='input OUT sorted mode',ylabel='output OUT sorted mode')
            ax.set_xticks(np.arange(k),np.arange(1,k+1));ax.set_yticks(np.arange(k),np.arange(1,k+1))
    fig.subplots_adjust(left=.08,right=.83,bottom=.08,top=.91,hspace=.35,wspace=.35)
    cax=fig.add_axes([.87,.15,.025,.68]);fig.colorbar(im,cax=cax,label='absolute direct Fourier coefficient [1/time]')
    fig.suptitle('C, p12: direct coefficients; both frequency signs, E0 unit modes',fontsize=11)
    for ext in ('pdf','png'):fig.savefig(b/('C_modal_coupling_matrix.'+ext),dpi=180,bbox_inches='tight');plots.append('C_modal_coupling_matrix.'+ext)
    plt.close(fig)
    if (b/'mu_interaction_rows.json').exists():
        rows=p.read(b/'C_interaction_rows.json')+p.read(b/'mu_interaction_rows.json')
        fig,axes=plt.subplots(2,2,figsize=(10,6.7),sharex=True)
        groups=[[('self',2,2,1),('sum',1,2,1)],[('difference',2,5,2),('sum',2,3,2),('sum',1,5,2)]]
        for panel,mechanisms in enumerate(groups):
            for kind,i,j,m in mechanisms:
                data=sorted([r for r in rows if r['kind']==kind and r['i_sorted']==i and r['j_sorted']==j and r['harmonic']==m],key=lambda r:r['mu'])
                label=f'{kind} {i},{j}, m={m}'
                valid=[r for r in data if r.get('modal_correspondence_status','PASS')=='PASS'];unresolved=[r for r in data if r.get('modal_correspondence_status','PASS')!='PASS']
                for col,metric in enumerate(('signed_detuning','sqrt_abs_product')):
                    line,=axes[panel,col].plot([r['mu'] for r in valid],[r[metric] for r in valid],'o-',label=label)
                    # Local sorted values at mu=0 remain visible, but no line
                    # claims correspondence through its unresolved modes4/5.
                    axes[panel,col].plot([r['mu'] for r in unresolved],[r[metric] for r in unresolved],'o',mfc='none',mec=line.get_color(),ms=8)
            axes[panel,0].axhline(0,color='grey',ls=':');axes[panel,0].set_ylabel('bare signed detuning [1/time]')
            axes[panel,1].set_ylabel('sqrt(abs(h_forward h_reverse)) [1/time]');axes[panel,0].legend(fontsize=8)
            for ax in axes[panel]:ax.grid(alpha=.2)
        for ax in axes[-1]:ax.set_xlabel('mu; five points only')
        fig.suptitle('Open mu=0 markers: local sorted modes; correspondence UNRESOLVED',fontsize=10)
        finish(fig,'mu_resonant_coupling_screening')
    p.write(b/'figure_provenance.json',{'figures':plots,'new_scientific_calls':0,'no_interpolation_of_parameter_map':True,
        'lines_between_five_points':'visual guides only, no adaptive search','modal_indices_are_sorted_positions':True});return plots

def beta_validate_scope(config):
    if config!=p.read(BETA_CONFIG) or config['authorization']!=BETA_AUTH:
        raise ValueError('Only the frozen five-angle screening is authorized')
    if (config['degree']!=12 or config['mu']!=.1 or config['lengths']!=[.9,1.1]
        or config['beta_points']!=[15.,30.,45.,60.,75.] or config['amplitude']!=.005
        or config['numerical_budget_seconds']>300 or config['new_ODE_integrations']
        or config['new_FEM_jobs'] or config['automatic_retry']):raise ValueError('Beta scope expansion')
    return config


def beta_verify_sources(config):
    parent=p.ROOT/config['parent']
    if p.sha(parent/'manifest.json')!=config['parent_manifest_sha256']:
        raise ValueError('Immutable resonant-screening parent changed')
    return cache(parent)


def beta_prepare(path=BETA_CONFIG):
    config=beta_validate_scope(p.read(path));beta_verify_sources(config)
    for old in BETA_OUTPUT.glob('*/provenance.json') if BETA_OUTPUT.exists() else []:
        item=p.read(old)
        if item['config']['authorization']==BETA_AUTH:
            if item['config']!=config:raise ValueError('Authorization already used with different settings')
            return old.parent,item,p.cache(old.parent)
    parent=p.read(p.ROOT/config['parent']/'provenance.json')
    files={**p.code_hashes(),Path(__file__).relative_to(p.ROOT).as_posix():p.sha(__file__),
        BETA_CONFIG.relative_to(p.ROOT).as_posix():p.sha(BETA_CONFIG)}
    item={'schema':'nlsp-beta-resonant-coupling-screening-v1','config':config,
        'physical_source_config':parent['physical_source_config'],'source_artifacts':parent['source_artifacts'],
        'source_manifests':{**parent['source_manifests'],'resonant_parent':{'path':config['parent'],'sha256':config['parent_manifest_sha256']}},
        'code_hashes':files,'initial_state':p.read(BETA_SMOKE/'initial_state.json'),'interpreter':sys.executable,
        'mathematical_operator':'unchanged full analytic Jacobian, E0 dual complex modes, direct resonant-term pair algebra',
        'no_new_ODE_or_Floquet':True}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=BETA_OUTPUT/fp
    if b.exists():raise ValueError('Preserve interrupted beta attempt')
    b.mkdir(parents=True);p.write(b/'provenance.json',item);p.write(b/'config.json',config)
    folder=b/'execution_code';folder.mkdir()
    for f in files:shutil.copyfile(p.ROOT/f,folder/Path(f).name)
    summary={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':0.,
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'CCX':0,'Gmsh':0,'linear_eigenproblems':0,'pointwise_Jacobian':0,'low_subspace_Ritz':0,'isolated_pair_static_algebra':0},
        'stages':{'source':'PASS','B':'NOT_RUN','C':'NOT_RUN','D':'NOT_RUN'},'points':{}}
    p.finalize(b,item,summary);return b,item,summary


def beta_cache(b):
    beta_verify_sources(p.read(Path(b)/'config.json'));return p.cache(b)


def beta_phase_invariance(problem,H,uncertainty,config):
    """Change time origin and complex modal gauges; products/eigenvalues stay invariant."""
    base=pair_rows(problem,H,uncertainty,config);k=len(problem['omega']);orders=config['harmonics']
    alpha=np.arange(1,k+1)*.173;D=np.exp(1j*np.r_[alpha,-alpha]);transformed=[]
    for m,h in zip(orders,H):transformed.append(np.exp(1j*m*.731)*(h/D[:,None])*D[None,:])
    other=pair_rows(problem,np.array(transformed),uncertainty,config)
    error=max(max(abs(complex(x['product']['real'],x['product']['imag'])-complex(y['product']['real'],y['product']['imag'])),
        abs(x['isolated_pair_averaging_growth_proxy']-y['isolated_pair_averaging_growth_proxy'])) for x,y in zip(base,other))
    if error>config['phase_invariance_absolute_gate']:raise ArithmeticError('Pair phase invariance failed')
    return {'maximum_product_or_growth_difference':error,'time_origin_shift_radians':.731,'modal_gauges_radians':alpha,
        'products_and_pair_eigenvalues_invariant':True,'new_Jacobian_or_ODE_calls':0}


def beta_correspondence(previous,problem,config):
    mac=physical_MAC(previous,problem);old,new=linear_sum_assignment(-mac)
    assigned=mac[old,new];margin=np.array([mac[i,j]-max(np.delete(mac[i],j)) for i,j in zip(old,new)])
    x1=np.r_[previous['q'],np.zeros(previous['d'].ndof)];x2=np.r_[problem['q'],np.zeros(problem['d'].ndof)]
    f1=l.a.physical_observations(previous['d'],x1[None,:])[0][0]
    f2=l.a.physical_observations(problem['d'],x2[None,:])[0][0]
    co=previous['d'].arms[0].coefficients;weights=np.array([co.m,co.m,co.m,co.jp+co.jb,co.jb,co.jp,co.jp])
    quad=np.ones(201);quad[[0,-1]]=.5;weight=np.array(problem['d'].lengths)[:,None,None]*quad[None,:,None]*weights/200
    plmac=float(np.sum(weight*f1*f2)**2/(np.sum(weight*f1*f1)*np.sum(weight*f2*f2)))
    valid=bool(min(assigned)>=config['correspondence_MAC_gate'] and min(margin)>=config['correspondence_margin_gate'] and plmac>=config['correspondence_MAC_gate'])
    return {'status':'PASS' if valid else 'UNRESOLVED','from_beta_deg':previous['beta_deg'],'to_beta_deg':problem['beta_deg'],
        'MAC_matrix':mac,'matching_sorted_indices':new+1,'assigned_MAC':assigned,'minimum_assigned_MAC':float(min(assigned)),
        'row_assignment_margins':margin,'planar_MAC':plmac,
        'comparison_contract':'common material coordinate and local physical fields with reference kinetic weights; no forced analytic branch_id'}


def beta_enrich_rows(problem,H,full,uncertainty,config,correspondence):
    rows=pair_rows(problem,H,uncertainty,config);k=len(problem['omega']);orders=config['harmonics']
    for row in rows:
        i=row['i_sorted']-1;j=row['j_sorted']-1;m=row['harmonic']
        rindex,sindex=(j,i) if row['kind']=='difference' else (i,k+j)
        observed=[]
        for index,ridx,sidx in [(orders.index(m),rindex,sindex),(orders.index(-m),sindex,rindex)]:
            h0=abs(full[0,index,ridx,sidx]);hhalf=abs(full[1,index,ridx,sidx])
            observed.append(None if min(h0,hhalf)<=1e-13 else float(np.log2(h0/hhalf)))
        expected=1 if m==1 else 2;order_pass=all(x is not None and abs(x-expected)<=config['amplitude_order_gate'] for x in observed)
        hr=H[orders.index(0),rindex,rindex];hs=H[orders.index(0),sindex,sindex]
        h=complex(row['coupling_forward']['real'],row['coupling_forward']['imag']);back=complex(row['coupling_reverse']['real'],row['coupling_reverse']['imag'])
        eigen=np.linalg.eigvals([[hr-1j*row['signed_detuning']/2,h],[back,hs+1j*row['signed_detuning']/2]])
        row.update(beta_deg=problem['beta_deg'],amplitude_order_expected=expected,amplitude_orders_observed=observed,
            amplitude_order_gate='PASS' if order_pass else 'UNRESOLVED',modal_correspondence_status=correspondence['status'],
            torsional_mass_fractions=[problem['torsional_mass_fraction'][i],problem['torsional_mass_fraction'][j]],
            pair_eigenvalues=[complex_row(v) for v in eigen],mean_diagonal=[complex_row(hr),complex_row(hs)])
        row['direct_term_candidate'] &= order_pass and correspondence['status']=='PASS'
        growth=row['isolated_pair_averaging_growth_proxy']
        row['characteristic_growth_time']=1/growth if row['direct_term_candidate'] else None
        row['growth_time_in_reference_periods']=problem['Omega']/(2*np.pi*growth) if row['direct_term_candidate'] else None
        row['pair_classification']='candidate_exponential_growth' if row['direct_term_candidate'] else 'unresolved' if not row['coupling_resolved'] or not order_pass else 'oscillatory_exchange' if row['kind']=='difference' and row['product']['real']<0 else 'oscillatory_pair_off_resonance'
    return rows


def beta_run(b,item,summary):
    if summary['completed'] or summary['hard_stop']:return summary
    config=item['config'];started=time.perf_counter();deadline=started+config['numerical_budget_seconds']-summary['numerical_seconds']
    previous=None;anchor=None;allrows=[]
    try:
        model,_,_=p.sources(item['physical_source_config'])
        for beta in (30.,15.,45.,60.,75.):
            if time.perf_counter()>deadline:raise TimeoutError('300-second beta ceiling')
            if beta==45.:previous=anchor
            point_start=time.perf_counter();prob=modal_problem(model,.1,config,summary['calls'],beta_deg=beta)
            summary['calls']['low_subspace_Ritz']+=1
            case=b/f'beta_{beta:02.0f}';case.mkdir();save_problem(case,'C',prob)
            planar_fields=l.a.physical_observations(prob['d'],np.r_[prob['q'],np.zeros(prob['d'].ndof)][None,:])[0][0]
            np.savez_compressed(case/'planar_physical_profile.npz',fields=planar_fields,x=np.linspace(0,1,201))
            physics={'beta_deg':beta,'lengths':prob['d'].lengths,'coefficients':[vars(a.coefficients) for a in prob['d'].arms],
                'local_to_global_bases':prob['d'].bases,'T7':prob['d'].T7,'joint_checks':prob['d'].joint_diagnostics(prob['q'],np.zeros(prob['d'].ndof)),
                'safety':prob['d'].safety(prob['q']),'actual_max_w':p.polynomial_w_max(prob['d'],prob['q'])[0]}
            if not all(a['mass_positive'] for a in physics['safety']):raise ArithmeticError('Probe mass positivity failed')
            if any(physics['joint_checks'][key]>2e-12 for key in ('translation_difference','orientation_difference','contraction_difference')):
                raise ArithmeticError('Probe rigid-joint compatibility failed')
            if anchor is not None and physics['coefficients']!=[vars(a.coefficients) for a in anchor['d'].arms]:
                raise ArithmeticError('Physical coefficients changed across beta')
            p.write(case/'geometry_and_safety.json',physics)
            if anchor is None:
                with np.load(p.ROOT/config['parent']/(config['reference_mu_label']+'_modes.npz')) as z:
                    checks={'frequency_difference':float(max(abs(z['omega']-prob['omega']))),
                        'M0_difference':float(max(abs(z['M0']-prob['d'].M0).flat)),
                        'K0_difference':float(max(abs(z['K0']-prob['d'].K).flat)),
                        'initial_shape_difference':float(max(abs(z['planar_initial_shape']-prob['q'])))}
                if max(checks.values())>config['coefficient_reference_absolute_gate']:raise ArithmeticError('Saved beta30 source mismatch')
                anchor=prob;correspondence={'status':'PASS','reference':'saved beta30 mu0.1','minimum_assigned_MAC':1.,'planar_MAC':1.}
                p.write(case/'reference_checks.json',checks)
            else:correspondence=beta_correspondence(previous,prob,config)
            p.write(case/'mode_correspondence.json',correspondence)
            summary['numerical_seconds']+=time.perf_counter()-point_start
            stage_b(case,item,summary,prob,deadline)
            post_start=time.perf_counter()
            with np.load(case/'C_harmonic_coefficients.npz') as z:H=z['leading_at_full_amplitude'];full=z['full'];uncertainty=z['numerical_difference']
            if beta==30.:
                with np.load(p.ROOT/config['parent']/(config['reference_mu_label']+'_coefficients.npz')) as z:err=float(np.max(abs(H-z['leading_at_full_amplitude'])))
                checks['harmonic_coefficient_difference']=err;p.write(case/'reference_checks.json',checks)
                if err>config['coefficient_reference_absolute_gate']:raise ArithmeticError('Saved harmonic coefficient mismatch')
            p.write(case/'pair_phase_invariance.json',beta_phase_invariance(prob,H,uncertainty,config))
            rows=beta_enrich_rows(prob,H,full,uncertainty,config,correspondence)
            reflection_errors=[];d=prob['d']
            for phase in (0.,.7):
                x=np.r_[prob['q']*np.cos(phase),-prob['Omega']*prob['q']*np.sin(phase)]
                J=d.jacobian(0,x);summary['calls']['pointwise_Jacobian']+=1
                reflection_errors.append(float(max(np.max(abs(J[np.ix_(d.state_out,d.state_plane)])),np.max(abs(J[np.ix_(d.state_plane,d.state_out)])))))
            if max(reflection_errors)>1e-10:raise ArithmeticError('Planar reflection block selection failed')
            p.write(case/'planar_reflection_selection.json',{'cross_blocks_max':max(reflection_errors),'no_equal_arm_exchange_symmetry_assumed':True})
            # stage_b once, phase/gauge control twice, pair enrichment twice
            # (base classification and explicitly retained eigenvalues).
            summary['calls']['isolated_pair_static_algebra']+=5*len(rows)
            p.write(case/'interaction_rows.json',rows);allrows+=rows
            summary['points'][str(beta)]={'Omega':prob['Omega'],'omega':prob['omega'],'correspondence':correspondence,
                'resolved_couplings':sum(r['coupling_resolved'] for r in rows),'candidates':sum(r['direct_term_candidate'] for r in rows)}
            previous=prob;summary['numerical_seconds']+=time.perf_counter()-post_start
            p.write(b/'interaction_rows.json',allrows);p.finalize(b,item,summary)
        candidates=[row for row in allrows if row['direct_term_candidate']]
        partial=any(v['correspondence']['status']!='PASS' for v in summary['points'].values())
        summary.update(completed=True,overall='DIRECT_GROWTH_CANDIDATE_IDENTIFIED' if candidates else 'MODE_IDENTIFICATION_PARTIAL' if partial else 'COUPLING_PRESENT_NO_DIRECT_GROWTH_CANDIDATE',
            no_direct_growth_status=None if candidates else 'NO_DIRECT_GROWTH_CANDIDATE_IN_TESTED_BETA_SET',candidate_count=len(candidates))
        summary['stages'].update(B='PASS',C='PASS',D='PARTIAL' if partial else 'PASS')
        summary['numerical_seconds']=time.perf_counter()-started
        p.write(b/'candidate_rows.json',candidates);beta_export(b);beta_plot(b)
    except Exception as error:
        summary.update(hard_stop=True,overall='NUMERICAL_PARTIAL' if isinstance(error,TimeoutError) else 'BLOCKED',stop_reason=type(error).__name__+': '+str(error))
        summary['numerical_seconds']=time.perf_counter()-started
    p.finalize(b,item,summary);return summary


def beta_export(b):
    rows=p.read(b/'interaction_rows.json')
    names=['beta_deg','Omega','kind','i_sorted','j_sorted','harmonic','signed_detuning','detuning_over_pump','sqrt_abs_product',
        'mean_shifted_detuning','two_coupling_over_mean_shifted_detuning','coupling_resolved','amplitude_order_gate','modal_correspondence_status','pair_classification','isolated_pair_averaging_growth_proxy']
    with (b/'interaction_coefficients.csv').open('w',encoding='utf8',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=names+['h_real','h_imag','k_real','k_imag','product_real','product_imag']);writer.writeheader()
        for r0 in sorted(rows,key=lambda r:(r['beta_deg'],r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])):
            r1={key:r0[key] for key in names}
            for key,source in [('h','coupling_forward'),('k','coupling_reverse'),('product','product')]:
                for part in ('real','imag'):r1[key+'_'+part]=r0[source][part]
            writer.writerow(r1)


def beta_plot(b):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=p.read(b/'interaction_rows.json');config=p.read(b/'config.json');plots=[]
    # Include the two closest self channels and retain the historical difference
    # pair. Lines assert correspondence only when every point passes its guide.
    channels=[('self',2,2,1),('self',3,3,2),('sum',1,4,2),('difference',2,5,2)]
    for metric,name,ylabel in [('mean_shifted_detuning','beta_resonant_detunings','direct mean-corrected detuning [1/time]'),
        ('two_coupling_over_mean_shifted_detuning','beta_resonant_interactions','2 sqrt(abs(h k)) / abs(effective detuning)')]:
        fig,axes=plt.subplots(2,2,figsize=(10,6.5),sharex=True)
        for ax,(kind,i,j,m) in zip(axes.flat,channels):
            data=sorted([r for r in rows if (r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])==(kind,i,j,m)],key=lambda r:r['beta_deg'])
            for r0 in data:
                color='red' if r0['direct_term_candidate'] else 'tab:blue'
                ax.plot(r0['beta_deg'],r0[metric],'o',color=color,mfc=color if r0['modal_correspondence_status']=='PASS' else 'none')
            if all(r['modal_correspondence_status']=='PASS' for r in data):ax.plot([r['beta_deg'] for r in data],[r[metric] for r in data],color='tab:blue',alpha=.5)
            if metric=='mean_shifted_detuning':ax.axhline(0,color='grey',ls=':')
            else:
                maximum=max(r[metric] for r in data)
                ax.set_ylim(0,max(.01,maximum*1.12))
                if kind!='difference':
                    if maximum>=.9:ax.axhline(1,color='grey',ls=':')
                    else:ax.text(.97,.92,'all ratios < 1',ha='right',transform=ax.transAxes,fontsize=8)
                else:ax.text(.97,.92,'exchange; ratio is not a growth test',ha='right',transform=ax.transAxes,fontsize=8)
            ax.grid(alpha=.2)
            ax.set_title(f'{kind} {i},{j}, m={m}; local sorted modes');ax.set_ylabel(ylabel,fontsize=8);ax.set_xticks(config['beta_points'])
        for ax in axes[-1]:ax.set_xlabel('beta [degree]')
        fig.suptitle('Five angles only; open markers: identification unresolved; red: direct pair candidate',fontsize=10);fig.tight_layout()
        for ext in ('pdf','png'):fig.savefig(b/(name+'.'+ext),dpi=180,bbox_inches='tight');plots.append(name+'.'+ext)
        plt.close(fig)
    p.write(b/'figure_provenance.json',{'figures':plots,'new_scientific_calls':0,'illustrated_channels':channels,
        'sorted_positions_not_forced_branches':True,'difference_channel_ratio_not_a_growth_criterion':True});return plots


def beta_main(argv=None):
    parser=argparse.ArgumentParser(description='Five-angle pointwise diagnostic; no IVP')
    mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--compute',action='store_true')
    mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path);parser.add_argument('--config',type=Path,default=BETA_CONFIG)
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();summary=beta_cache(b)
        if args.plot_only:beta_plot(b);p.finalize(b,p.read(b/'provenance.json'),summary)
    else:b,item,summary=beta_prepare(args.config);summary=beta_run(b,item,summary)
    print(json.dumps(p.plain({'bundle':str(b),'overall':summary['overall'],'stages':summary['stages'],'calls':summary['calls'],
        'seconds':summary['numerical_seconds'],'stop_reason':summary.get('stop_reason')}),indent=2))


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--through-stage',choices=('C-only','all'),default='all');args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=cache(b)
        if args.plot_only:
            plot_saved(b);p.finalize(b,p.read(b/'provenance.json'),s)
    else:
        b,item,s=prepare(args.config)
        try:s=run(b,item,s,args.through_stage)
        except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as error:
            s.update(overall='BLOCKED',hard_stop=True,stop_reason=type(error).__name__+': '+str(error));p.finalize(b,item,s)
    print(json.dumps(p.plain({'bundle':str(b),'overall':s['overall'],'stages':s['stages'],'calls':s['calls'],'seconds':s['numerical_seconds'],'stop_reason':s.get('stop_reason')}),indent=2))
