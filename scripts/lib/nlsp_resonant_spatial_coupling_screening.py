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

def modal_problem(model,mu,config,counter):
    pc={**p.read(c.CONFIG),'lengths':[1-mu,1+mu],'planar_modes_inspected':config['planar_modes_inspected']}
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
        'profiles':np.array(profiles),'torsional_mass_fraction':torsion,'guard_omega':guard_omega,'mu':mu,'duality':duality,'eigen_residual':residual}

def harmonics(values,phases,orders):
    return np.array([np.einsum('t,tij->ij',np.exp(-1j*m*phases),values)/len(phases) for m in orders])

def projected_at(problem,state,counter):
    d=problem['d'];J=d.jacobian(0,state);counter['pointwise_Jacobian']+=1
    A=J[np.ix_(d.state_out,d.state_out)]
    return problem['L']@(A-problem['A0'])@problem['R']

def harmonic_sample(problem,factor,count,orders,counter,deadline):
    phases=np.arange(count)*2*np.pi/count;values=[];d=problem['d'];q=problem['q']*factor
    for phase in phases:
        if time.perf_counter()>deadline:raise TimeoutError('900-second screening ceiling')
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
