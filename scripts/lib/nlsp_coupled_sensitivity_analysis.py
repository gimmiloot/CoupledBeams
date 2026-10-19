"""Full-DOF tangent propagation and saved-data sensitivity diagnostics.

Only the initial seeds are restricted. The propagated response retains every
out-of-plane Galerkin DOF, including the common joint coordinates.
"""
from functools import lru_cache
from pathlib import Path
import time
import numpy as np
from scipy.integrate import Radau,DOP853
from scipy.linalg import eigh
from scipy.sparse import csr_matrix,kron,eye


class SavedPlane:
    """Existing accepted Radau cubic output, loaded once and never refitted."""
    def __init__(self,path,disc):
        with np.load(path) as z:self.old,self.end,self.initial,self.Q=(z[k].copy() for k in ('t_old','t','y_old','Q'))
        self.disc=disc
    def __call__(self,t):
        if t<self.old[0] or t>self.end[-1]:raise ValueError('Outside the actual accepted base prefix')
        i=min(int(np.searchsorted(self.end,t,side='left')),len(self.end)-1)
        x=(t-self.old[i])/(self.end[i]-self.old[i])
        return self.disc.embed_plane(self.initial[i]+self.Q[i]@np.array([x,x*x,x*x*x]))


def reference_subspace(disc):
    ids=disc.out;omega,vectors=disc.eigenpairs(plane=False);fractions=[]
    for j in range(min(12,len(omega))):
        q=np.zeros(disc.ndof);q[ids]=vectors[:,j]
        parts=disc.reference_energy_parts(np.zeros(disc.ndof),q)
        fractions.append(float(parts[:,1].sum()/parts.sum()))
    bending=[j for j,f in enumerate(fractions) if f<.5][:2]
    torsion=next((j for j,f in enumerate(fractions) if f>=.5),None)
    if len(bending)<2 or torsion is None:raise ArithmeticError('Low-mode bending/torsion initial subspace not identified')
    chosen=bending+[torsion];n=len(ids);basis=np.zeros((2*n,6))
    for k,j in enumerate(chosen):basis[:n,k]=vectors[:,j]/omega[j];basis[n:,k+3]=vectors[:,j]
    metric=np.zeros((2*n,2*n));metric[:n,:n]=disc.K[np.ix_(ids,ids)];metric[n:,n:]=disc.M0[np.ix_(ids,ids)]
    gram=basis.T@metric@basis
    if np.max(abs(gram-np.eye(6)))>2e-10:raise ArithmeticError('Reference-energy initial seeds not orthonormal')
    return basis,metric,{'mode_indices_zero_based':chosen,'omega':omega[chosen],'torsional_mass_fractions':np.array(fractions)[chosen],
        'initial_dimension':6,'full_outplane_state_dimension':2*n,'reference_gram_error':float(np.max(abs(gram-np.eye(6)))),
        'restriction':'six initial states only; every outplane response DOF is propagated','zero_state_spectrum':omega}


def tangent_attempt(b,item,s,label,d,base_case,times,basis,metric,deadline,zero=False,*,max_step=None,dense_archive=False,method='Radau',sampled_dense_only=False,t_start=0.,initial_matrix=None,operator_override=None,base_evaluator=None):
    from scripts.lib import nlsp_coupled_sensitivity_pilot as p
    if method not in ('Radau','DOP853') or (dense_archive and method!='Radau'):raise ValueError('Unsupported tangent/dense method')
    if label in s['cases']:raise ValueError('A tangent attempt already exists')
    if times[0]!=t_start or times[-1]<=t_start:raise ValueError('Tangent initial time mismatch')
    if s['calls']['tangent_ODE']>=item['config'].get('maximum_tangent_trajectories',10**9):raise TimeoutError('Declared tangent call ceiling')
    if initial_matrix is not None and (initial_matrix.shape!=basis.shape or not np.isfinite(initial_matrix).all()):raise ValueError('Wrong continued matrix shape')
    case=b/label;case.mkdir(parents=True);p.snapshot_case(case)
    row={'label':label,'kind':'tangent_ODE','status':'STARTED','automatic_retry':False,'zero_base':zero}
    s['attempts'].append(row);s['calls']['tangent_ODE']+=1;s['cases'][label]={'status':'STARTED'};p.finalize(b,item,s)
    base=None if zero else (base_evaluator if base_evaluator is not None else SavedPlane(base_case/'dense.npz',d))
    J0=d.jacobian(0,np.zeros(2*d.ndof))[np.ix_(d.state_out,d.state_out)]
    @lru_cache(maxsize=16)
    def operator(t):
        if zero:return J0
        if operator_override is not None:return operator_override(t,base(t))
        J=d.jacobian(t,base(t))
        return J[np.ix_(d.state_out,d.state_out)]
    n,k=basis.shape;start_matrix=basis if initial_matrix is None else initial_matrix
    initial=start_matrix.flatten(order='F');history=np.empty((len(times),n,k));history[0]=start_matrix;cursor=1;steps=[];records=[];actual_times=[];solver=None;failure=None;start=time.perf_counter()
    limit=min(deadline,start+item['config']['per_tangent_run_seconds'])
    atol=item['config']['atol_relative']/np.sqrt(n)/np.sqrt(np.diag(metric));atol=np.tile(atol,k)
    cutoff=max(np.sqrt(a.coefficients.C/a.coefficients.jp) for a in d.arms)
    max_step=2*np.pi/cutoff/24 if max_step is None else float(max_step)
    try:
        def rhs(t,y):
            if time.perf_counter()>limit:raise TimeoutError('Declared tangent budget')
            return (operator(float(t))@y.reshape(n,k,order='F')).flatten(order='F')
        def jac(t,y):return kron(eye(k,format='csc'),csr_matrix(operator(float(t))),format='csc')
        kw={'jac':jac} if method=='Radau' else {}
        solver=(Radau if method=='Radau' else DOP853)(rhs,float(t_start),initial,times[-1],rtol=item['config']['rtol'],atol=atol,max_step=max_step,**kw)
        while solver.status=='running':
            if time.perf_counter()>limit:raise TimeoutError('Declared tangent budget')
            old=solver.t;message=solver.step()
            if solver.status=='failed':raise ArithmeticError(message or 'Tangent integration failed')
            steps.append(solver.t-old);end=int(np.searchsorted(times,solver.t,side='right'))
            actual_times.append(float(solver.t))
            if not sampled_dense_only or end>cursor or dense_archive:
                dense=solver.dense_output()
                if dense_archive:records.append({'t_old':dense.t_old,'t':dense.t,'y_old':dense.y_old.copy(),'Q':dense.Q.copy()})
                if end>cursor:history[cursor:end]=dense(times[cursor:end]).T.reshape(end-cursor,n,k,order='F');cursor=end
    except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as error:failure=type(error).__name__+': '+str(error)
    record={'status':'PASS' if failure is None and cursor==len(times) else 'PARTIAL','failure':failure,'samples':cursor,
        'target_end':float(times[-1]),'actual_accepted_end':float(solver.t if solver else 0.),'sampled_end':float(times[cursor-1]),
        'seconds':time.perf_counter()-start,'accepted_steps':len(steps),'nfev':0 if solver is None else solver.nfev,
        'njev':0 if solver is None else solver.njev,'nlu':0 if solver is None else solver.nlu,'initial_dimension':k,
        'full_response_dimension':n,'automatic_retry':False,'operator_source':'full coupled analytic Jacobian restricted after verified invariant-block separation',
        'max_step':max_step,'rtol':item['config']['rtol'],'atol_componentwise':atol,
        'operator_cache':str(operator.cache_info()),'method':method,'sampled_dense_only':sampled_dense_only,
        'actual_min_step':float(min(steps)) if steps else None,'actual_max_step':float(max(steps)) if steps else None}
    if operator_override is not None or base_evaluator is not None:
        record['operator_source']='caller-supplied verified tangent operator; see scoped operator provenance'
        record['caller_base_evaluator']=base_evaluator is not None
        record['operator_override']=operator_override is not None
    if dense_archive and records:
        if t_start:p.save_dense_records(case/'dense.npz',records,t_start=float(t_start))
        else:p.save_dense_records(case/'dense.npz',records)
        record['actual_min_step']=float(min(steps));record['actual_max_step']=float(max(steps))
    np.savez_compressed(case/'tangent.npz',times=times[:cursor],vectors=history[:cursor],basis=basis,metric=metric)
    np.savez_compressed(case/'accepted_times.npz',times=np.r_[t_start,actual_times] if t_start else np.r_[0.,np.cumsum(steps)])
    if t_start:
        record['physical_start']=float(t_start);record['original_initial_basis_preserved']=True
        record['operator_override']=operator_override is not None
        np.savez_compressed(case/'restart_state.npz',initial_matrix=start_matrix,actual_final_matrix=solver.y.reshape(n,k,order='F') if solver is not None else start_matrix,original_basis=basis,metric=metric)
    else:
        # Preserve the native endpoint for a later physical-state restart;
        # DOP853 has no saved Radau dense-coefficient archive.
        np.savez_compressed(case/'restart_state.npz',initial_matrix=start_matrix,actual_final_matrix=solver.y.reshape(n,k,order='F') if solver is not None else start_matrix,original_basis=basis,metric=metric)
    p.write(case/'integration.json',record);row['status']=record['status'];s['cases'][label]=record;s['numerical_seconds']+=record['seconds']
    if record['status']!='PASS':s.update(overall='NUMERICAL_FEASIBILITY_PARTIAL',hard_stop=True,stop_reason=failure or 'Incomplete tangent prefix')
    p.finalize(b,item,s);return history[:cursor],record


def energy_norm(values,metric):return np.sqrt(np.maximum(0,np.einsum('...i,ij,...j->...',values,metric,values)))

def full_tangent_error(d,full,base,tangent,direction,delta,metric):
    predicted=np.einsum('tij,j->ti',tangent,direction);observed=full[:,d.state_out]/delta
    error=energy_norm(observed-predicted,metric);scale=max(float(energy_norm(predicted,metric).max()),1.)
    planar=full[:,d.state_plane]-base[:,d.state_plane]
    plmetric=np.zeros((len(d.state_plane),len(d.state_plane)));n=len(d.plane)
    plmetric[:n,:n]=d.K[np.ix_(d.plane,d.plane)];plmetric[n:,n:]=d.M0[np.ix_(d.plane,d.plane)]
    return {'normalized_out_error_max':float(error.max()),'normalized_out_error_relative':float(error.max()/scale),
        'planar_change_energy_norm_max':float(energy_norm(planar,plmetric).max()),'delta':float(delta),
        'out_error_curve':error,'plane_change_curve':energy_norm(planar,plmetric),'predicted':predicted,'observed':observed}

def measured_orders(rows):
    orders=[]
    for first,second in zip(rows,rows[1:]):
        ratio=first['delta']/second['delta']
        orders.append({'normalized_out_order':float(np.log(max(first['normalized_out_error_relative'],1e-300)/max(second['normalized_out_error_relative'],1e-300))/np.log(ratio)),
            'planar_change_order':float(np.log(max(first['planar_change_energy_norm_max'],1e-300)/max(second['planar_change_energy_norm_max'],1e-300))/np.log(ratio))})
    return orders

def two_amplitude_diagnostic(first,second,tangent,metric):
    """Read-only cancellation of the leading even normalized-response term.

    Inputs are full responses already divided by delta and delta/2. This
    conditional diagnostic neither certifies a delta->0 limit nor changes
    any frozen pass/fail gate or the original integration outputs.
    """
    estimate=(4*second-first)/3
    return {'conditional_limit_discrepancy':float(energy_norm(estimate-tangent,metric).max()),
        'observed_amplitude_change':float(energy_norm(first-second,metric).max()),
        'first_error':float(energy_norm(first-tangent,metric).max()),
        'second_error':float(energy_norm(second-tangent,metric).max()),
        'assumption':'leading normalized nonlinear term proportional to delta squared; diagnostic only',
        'does_not_override_frozen_gate':True}

def response_error_details(d,error,delta,metric):
    """Physical coordinate errors and separate fixed K/M error norms."""
    n=len(d.out);curve=energy_norm(error,metric);index=int(np.argmax(curve))
    q=np.zeros(d.ndof);v=q.copy();q[d.out]=error[index,:n];v[d.out]=error[index,n:]
    parts=d.reference_energy_parts(q,v)
    state=np.zeros((len(error),2*d.ndof));state[:,d.state_out]=delta*error
    physical,_=physical_observations(d,state)
    potential=energy_norm(error[:,:n],metric[:n,:n]);kinetic=energy_norm(error[:,n:],metric[n:,n:])
    return {'initial_relative_error':float(curve[0]),
        'error_energy_fractions_arm_bending_torsion_at_max':parts/max(parts.sum(),1e-300),
        'physical_coordinate_error_max_u_w_v_Phi_psi_theta_c':np.max(abs(physical),axis=(0,1,2)),
        'normalized_coordinate_K_error_max':float(potential.max()),'normalized_velocity_M_error_max':float(kinetic.max()),
        'velocity_error_energy_fraction_at_max_total_error':float(kinetic[index]**2/max(curve[index]**2,1e-300))}

def diagnose_stopped_prefix(b):
    """Postprocess the saved failed tangent gate, with zero new IVP/eigensolves."""
    from scripts.lib import nlsp_coupled_sensitivity_pilot as p
    import csv
    b=Path(b);s=p.read(b/'summary.json');c=p.read(b/'config.json')
    if not s['hard_stop'] or s['overall']!='BLOCKED_BY_TANGENT_CHECK':raise ValueError('Only the saved stopped tangent-check prefix')
    if (b/'stopped_prefix_diagnostics.json').exists():return p.read(b/'stopped_prefix_diagnostics.json')
    model,_,_=p.sources(c);d=p.make_disc(c,'A',16,model)
    with np.load(b/'A_short_tangent/tangent.npz') as z:times,Z,basis,metric=(z[k].copy() for k in ('times','vectors','basis','metric'))
    base=p.evaluate_dense_records(b/'A_p16_base/dense.npz',times);embedded=np.zeros((len(times),2*d.ndof));embedded[:,d.state_plane]=base
    static=p.read(b/'A_p16_static.json');delta0=c['perturbation_energy_norm_over_base']*np.sqrt(2*static['energy_at_release'])
    report={'overall':s['overall'],'actual_stop_reason':s['stop_reason'],'source_execution_manifest_sha256':p.sha(b/'stopped_execution_evidence/manifest.json'),
        'new_nonlinear_ODE':0,'new_tangent_ODE':0,'new_static_equilibria':0,'new_eigenanalyses':0,
        'amplitude_floor_not_modified':True,'gains_not_computed':True,'geometric_contrast_not_established':True,
        'shapes':{},'pure_matrix_directional_checks':{},'reference_subspace':{'initial_dimension':basis.shape[1],
        'full_out_response_dimension':basis.shape[0],'initial_basis_gram_error':float(np.max(abs(basis.T@metric@basis-np.eye(basis.shape[1]))))}}
    seed_fractions=[]
    for index in range(3):
        q=np.zeros(d.ndof);v=q.copy();v[d.out]=basis[len(d.out):,index+3]
        parts=d.reference_energy_parts(q,v);seed_fractions.append(float(parts[:,1].sum()/parts.sum()))
    report['reference_subspace']['torsional_mass_fractions']=seed_fractions
    vectors=basis[len(d.out):,3:];n=len(d.out)
    report['reference_subspace']['selected_omega']=np.sqrt(np.einsum('ik,ij,jk->k',vectors,metric[:n,:n],vectors)/np.einsum('ik,ij,jk->k',vectors,metric[n:,n:],vectors))
    arrays=b/'check_arrays';arrays.mkdir(exist_ok=True)
    table=[]
    for shape,index,field in [('bending',0,2),('torsional',2,3)]:
        direction=np.eye(6)[index];predicted=np.einsum('tij,j->ti',Z,direction);rows=[];responses=[];curves=[];plane_curves=[];details=[]
        for j in range(2):
            label=f'A_check_{shape}_{j}'
            with np.load(b/label/'trajectory.npz') as z:full=z['state'];assert np.array_equal(times,z['times'])
            delta=delta0/2**j;error=full_tangent_error(d,full,embedded,Z,direction,delta,metric)
            row={k:v for k,v in error.items() if k not in ('out_error_curve','plane_change_curve','predicted','observed')}
            row['time_at_max_error']=float(times[np.argmax(error['out_error_curve'])]);row['gate_relative_scale']=max(float(energy_norm(predicted,metric).max()),1.)
            rows.append(row);table.append({'shape':shape,'level':j,**row});responses.append(error['observed']);curves.append(error['out_error_curve']);plane_curves.append(error['plane_change_curve'])
            details.append(response_error_details(d,error['observed']-predicted,delta,metric))
        orders=measured_orders(rows);first,last=(r['normalized_out_error_relative'] for r in rows)
        trend=last<=first*1.05;resolved=first>10*c['numerical_floor_relative']
        order_ok=all(row['normalized_out_order']>=1.5 for row in orders) if resolved and last>c['numerical_floor_relative'] else True
        report['shapes'][shape]={'rows':rows,'orders':orders,'absolute_relative_gate_pass':last<=c['tangent_full_relative_gate'],
            'declared_monotone_gate_pass':trend,'order_resolved_above_10_floor':resolved,
            'observed_out_error_status':'PASS' if last<=c['tangent_full_relative_gate'] and (last<=c['numerical_floor_relative'] or trend) and order_ok else 'FAIL',
            'physical_error_diagnostics':details,
            'floor_qualification':'order is diagnostic only below 10 times the predeclared floor',
            'two_amplitude_postprocessing':two_amplitude_diagnostic(responses[0],responses[1],predicted,metric)}
        predicted_full=np.zeros_like(embedded);predicted_full[:,d.state_out]=predicted
        observed_full=np.zeros_like(embedded);observed_full[:,d.state_out]=responses[-1]
        physical_pred=d.arms[0].reconstruct_series(predicted_full[:,:d.ndof]@d.maps[0].T,np.array([d.lengths[0]/2]))[:,0,field]
        physical_obs=d.arms[0].reconstruct_series(observed_full[:,:d.ndof]@d.maps[0].T,np.array([d.lengths[0]/2]))[:,0,field]
        np.savez_compressed(arrays/(shape+'_stopped_diagnostic.npz'),times=times,tau=times/p.read(b/'frozen_excitation.json')['A']['Tref'],
            tangent=physical_pred,full_over_delta=physical_obs,error=np.array(curves),plane_change=np.array(plane_curves),delta=[delta0,delta0/2],field=field)
        # Independent central directional differences at actual base states.
        jac_rows=[];seed=np.zeros(2*d.ndof);seed[d.state_out]=basis@direction
        for j in (0,25,50,100):
            x=embedded[j];exact=d.jacobian(times[j],x)@seed;central=[]
            for h in (delta0,delta0/2):central.append((d.rhs(times[j],x+h*seed)-d.rhs(times[j],x-h*seed))/(2*h))
            extrap=(4*central[1]-central[0])/3
            jac_rows.append({'time':float(times[j]),'directional_relative':float(np.linalg.norm(extrap-exact)/max(np.linalg.norm(exact),1e-30))})
        report['pure_matrix_directional_checks'][shape]=jac_rows
    with (b/'tangent_full_amplitude_diagnostics.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(table[0]));writer.writeheader();writer.writerows(table)
    atol=np.array(p.read(b/'A_check_torsional_0/time_settings.json')['atol'])[d.state_out]
    report['component_atol_reference_norm_over_delta']=[float(energy_norm(atol,metric)/(delta0/2**j)) for j in range(2)]
    report['atol_qualification']='norm of positive component atol vector, diagnostic scale, neither achieved error nor rigorous bound; divided by the actual perturbation norm'
    with np.load(b/'A_zero_tangent/tangent.npz') as z:
        gram=np.einsum('tij,ik,tkl->tjl',z['vectors'],z['metric'],z['vectors']);report['autonomous_gram_error']=float(np.max(abs(gram-np.eye(basis.shape[1]))))
    with np.load(b/'A_p16_base/trajectory.npz') as z:t,state,E,drift=(z[k].copy() for k in ('times','state','energy','energy_relative_drift'))
    sampled_safety=[d.safety(x[:d.ndof]) for x in state];safety=[]
    for arm in range(2):
        safety.append({key:(min(r[arm][key] for r in sampled_safety) if key in ('min_one_plus_c','relative_mass_lower_bound') else max(r[arm][key] for r in sampled_safety))
            for key in ('min_one_plus_c','max_abs_c','max_rotation_vector_norm','max_abs_axial_gradient','max_abs_transverse_gradient','max_L_abs_curvature','relative_mass_lower_bound','relative_mass_upper_bound','relative_mass_condition_upper_bound')})
    report['sampled_base_safety']=safety
    with np.load(b/'A_p16_base/physical_observations.npz') as z:obs=z['observations'];tau=z['tau']
    with (b/'A_planar_motion_observations.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.writer(stream);writer.writerow(['time','tau','energy','relative_energy_drift']+[f'arm{i+1}_{field}_{position}' for i in range(2) for field in ('u','w','theta','c') for position in ('quarter','middle','three_quarters')])
        for j in range(len(t)):writer.writerow([t[j],tau[j],E[j],drift[j]]+[obs[j,i,k,field] for i in range(2) for field in (0,1,5,6) for k in range(3)])
    execution=[]
    for label,r in s['cases'].items():execution.append({'label':label,**{key:r.get(key) for key in ('status','target_end','actual_accepted_end','samples','accepted_steps','nfev','njev','nlu','seconds','energy_relative_drift_max')}})
    report['executions']=execution
    with (b/'actual_execution_summary.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(execution[0]));writer.writeheader();writer.writerows(execution)
    p.write(b/'stopped_prefix_diagnostics.json',report);return report

def gain_result(d,history,times,basis,metric):
    values=[];combinations=[]
    for Z in history:
        gram=Z.T@metric@Z;lam,vec=eigh(gram,subset_by_index=[gram.shape[0]-1,gram.shape[0]-1])
        values.append(float(np.sqrt(max(lam[0],0))));combinations.append(vec[:,0])
    gains=np.asarray(values);index=int(np.argmax(gains));combination=combinations[index]
    initial=basis@combination;response=history[index]@combination;n=len(d.out)
    q=np.zeros(d.ndof);v=q.copy();q[d.out]=response[:n];v[d.out]=response[n:]
    parts=d.reference_energy_parts(q,v)
    return {'gain':float(gains[index]),'time_at_max':float(times[index]),'sample_index_at_max':index,
        'initial_combination':combination,'initial_state_out':initial,'response_at_max':response,
        'reference_energy_distribution_at_max_by_arm_bending_torsion':parts,
        'restricted_initial_dimension':basis.shape[1],'full_outplane_response_dimension':len(initial),
        'sampled_time_maximum':True,'lower_bound_for_full_discrete_initial_space':True,
        'infinite_dimensional_supremum_established':False},gains


def physical_observations(d,states):
    x=np.linspace(0,1,201);profiles=[];observations=[]
    for i,a in enumerate(d.arms):
        raw=states[:,:d.ndof]@d.maps[i].T;velocity=states[:,d.ndof:]@d.maps[i].T
        profiles.append(a.reconstruct_series(raw,x*a.length));observations.append(a.reconstruct_series(raw,np.array([.25,.5,.75])*a.length))
    return np.stack(profiles,axis=1),np.stack(observations,axis=1)

def save_base_observations(case,d,Tref):
    from scripts.lib.nlsp_coupled_sensitivity_pilot import write
    with np.load(case/'trajectory.npz') as z:times,states=z['times'],z['state']
    ids=np.linspace(0,len(times)-1,min(11,len(times)),dtype=int);checks=[]
    for index in ids:
        J=d.jacobian(times[index],states[index]);off=max(np.max(abs(J[np.ix_(d.state_out,d.state_plane)])),np.max(abs(J[np.ix_(d.state_plane,d.state_out)])))
        rhs=float(np.max(abs(d.rhs(times[index],states[index])[d.state_out])))
        if off>2e-12 or rhs>2e-12:raise ArithmeticError('Actual planar base tangent blocks do not separate')
        checks.append({'physical_time':float(times[index]),'off_diagonal_J_max':float(off),'planar_rhs_out_max':rhs})
    write(case/'actual_planar_invariance.json',{'status':'PASS','samples':checks,'blocks_not_forcibly_zeroed':True})
    profiles,observations=physical_observations(d,states)
    snapshots=np.array([np.argmin(abs(times-t)) for t in np.linspace(0,times[-1],5)])
    np.savez_compressed(case/'physical_observations.npz',times=times,tau=times/Tref,x=np.linspace(0,1,201),
        observations=observations,profiles=profiles[snapshots],snapshot_times=times[snapshots],
        joint_coordinates=states[:,:d.ndof]@d.W[d.joint_ids].T)

def plot_saved(b):
    from scripts.lib.nlsp_coupled_sensitivity_pilot import read,write
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    b=Path(b);s=read(b/'summary.json');plots=[]
    def finish(fig,name):
        fig.tight_layout()
        for suffix in ('pdf','png'):fig.savefig(b/(name+'.'+suffix),dpi=180,bbox_inches='tight');plots.append(name+'.'+suffix)
        plt.close(fig)
    available=[]
    for name in ('A','B'):
        path=b/(name+'_p24_base')/'physical_observations.npz'
        if not path.exists():path=b/(name+'_p16_base')/'physical_observations.npz'
        if path.exists():available.append((name,path))
    if available:
        fig,axes=plt.subplots(len(available),2,figsize=(10,3.3*len(available)),squeeze=False)
        for row,(name,path) in enumerate(available):
            with np.load(path) as z:
                for arm in range(2):axes[row,arm].plot(z['tau'],z['observations'][:,arm,1,1]);axes[row,arm].set(title=f'Geometry {name}, arm {arm+1}',xlabel='t/Tref',ylabel='local w(Li/2)')
        finish(fig,'coupled_planar_motion')
    if (b/'stopped_prefix_diagnostics.json').exists():
        fig,axes=plt.subplots(2,2,figsize=(10,6.5))
        for column,shape in enumerate(('bending','torsional')):
            path=b/'check_arrays'/(shape+'_stopped_diagnostic.npz')
            with np.load(path) as z:
                axes[0,column].plot(z['tau'],z['tangent'],label='tangent')
                axes[0,column].plot(z['tau'],z['full_over_delta'],ls='--',label='full / delta, smaller delta')
                axes[0,column].set(title='A '+shape,xlabel='t/Tref',ylabel=('v' if shape=='bending' else 'Phi')+' at L1/2 per unit initial E0 norm');axes[0,column].legend(fontsize=8)
                for j in range(2):axes[1,column].semilogy(z['tau'][1:],z['error'][j,1:],label=f'delta = {z["delta"][j]:.3g}')
                axes[1,column].axhline(1e-8,color='gray',ls=':',label='predeclared floor')
                axes[1,column].set(xlabel='t/Tref',ylabel='E0 norm of full/delta - tangent',title='Saved amplitude check; no new integration');axes[1,column].legend(fontsize=8)
        finish(fig,'tangent_full_perturbations')
    elif (b/'tangent_checks.json').exists():
        checks=read(b/'tangent_checks.json');fig,axes=plt.subplots(1,2,figsize=(10,3.8))
        for ax,name in zip(axes,('A','B')):
            for shape in ('bending','torsional','mixed'):
                path=b/'check_arrays'/(name+'_'+shape+'.npz')
                if not path.exists():continue
                with np.load(path) as z:
                    ax.plot(z['tau'],z['tangent_physical_v'],label=shape+' tangent')
                    ax.plot(z['tau'],z['full_smallest_delta_physical_v'],ls='--',label=shape+' full / delta')
            ax.set(title='Geometry '+name,xlabel='t/Tref',ylabel='v1(L1/2) per E0-normalized delta');ax.legend(fontsize=7)
        finish(fig,'tangent_full_perturbations')
    if (b/'gain_results.json').exists():
        g=read(b/'gain_results.json');fig,axes=plt.subplots(1,2,figsize=(10,3.8))
        for name in ('A','B'):
            for degree in (16,24):
                path=b/(f'{name}_p{degree}_gain.npz')
                if not path.exists():continue
                with np.load(path) as z:axes[0].plot(z['tau'],z['gain'],label=f'{name}, p{degree}')
        axes[0].axhline(1,color='k',lw=.7);axes[0].set(xlabel='t/Tref',ylabel='G_perp^(6)(t)',title='Restricted initial subspace; full response DOFs');axes[0].legend()
        for name in ('A','B'):
            path=b/(name+'_finite_check.npz')
            if path.exists():
                with np.load(path) as z:axes[1].plot(z['tau'],z['tangent_norm'],label=name+' tangent');axes[1].plot(z['tau'],z['full_norm'],ls='--',label=name+' full')
        axes[1].set(xlabel='t/Tref',ylabel='norm / initial perturbation norm',title='Finite perturbation of maximizing seed');axes[1].legend(fontsize=8)
        finish(fig,'restricted_finite_time_gain')
        fig,axes=plt.subplots(2,2,figsize=(10,6))
        for row,name in enumerate(('A','B')):
            path=b/(name+'_optimal_profiles.npz')
            if not path.exists():continue
            with np.load(path) as z:
                for arm in range(2):
                    ax=axes[row,arm];ax.plot(z['x'],z['initial'][arm,:,2],label='initial v');ax.plot(z['x'],z['response'][arm,:,2],label='response v')
                    ax.set(title=f'{name}, arm {arm+1}; E0-normalized seed',xlabel='s/Li',ylabel='v per unit initial E0 norm');ax.legend()
        finish(fig,'most_amplified_perturbation_shapes')
    write(b/'figure_provenance.json',{'figures':plots,'new_scientific_calls':0,'native_3D_data':False,'no_phase_amplitude_time_fitting':True})
    return plots
