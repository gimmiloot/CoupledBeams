"""Second variation of saved quartic action and selected resonant-pair algebra.

No action rederivation, eigenproblem, ODE, Floquet or parameter search.
The independent action assembly retains the existing arm maps/quadrature.
"""
from __future__ import annotations
import argparse,hashlib,json,shutil,sys,time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scripts.lib import nlsp_resonant_spatial_coupling_screening as s
from scripts.lib.weakly_nonlinear_spatial_dynamics import _SparseCompiled,_POT_NAMES,_POT_FIELDS,_POT_DERIV,_JET_NAMES

p=s.p;rod=p.rod
CONFIG=p.ROOT/'data/input/nlsp_out_of_plane_instability_feasibility.json'
OUTPUT=p.ROOT/'results/nlsp_out_of_plane_instability_feasibility'
SMOKE=p.ROOT/'results/_smoke/nlsp_out_of_plane_instability_feasibility'
AUTH='explicit_user_out_of_plane_instability_feasibility_2026_10_11'
OUT_FIELDS=('v','Phi','psi')


def out_degree(monomial):return sum(i<42 and i%7 in (2,3,4) for i in monomial)


def select_out(polynomial,degree):
    return rod.Polynomial({m:c for m,c in polynomial.terms.items() if out_degree(m)==degree})


def restrict_plane(polynomial):
    return rod.Polynomial({m:c for m,c in polynomial.terms.items() if out_degree(m)==0})


def second_variation(model):
    T=select_out(model.T4,2);V=select_out(model.V4,2);L=T-V;checks={}
    for field in OUT_FIELDS:
        el=(L.derivative(field+'_t').total_derivative('t')+
            L.derivative(field+'_s').total_derivative('s')-L.derivative(field))
        saved=select_out(model.residual_a[rod.FIELD_ORDER.index(field)],1)
        checks[field]={'exact_residual_difference_terms':len((el-saved).terms),'linear_OUT_residual':el.serialize()}
    odd={'T':sum(out_degree(m)%2 for m in model.T4.terms),'V':sum(out_degree(m)%2 for m in model.V4.terms)}
    if any(v['exact_residual_difference_terms'] for v in checks.values()) or any(odd.values()):
        raise ArithmeticError('Hard gate A: action/residual or reflection mismatch')
    return T,V,{'exact_action_residual_checks':checks,'odd_OUT_terms':odd,
        'term_counts':{name:{str(n):len(poly.homogeneous(n).terms) for n in (2,3,4)} for name,poly in [('T',T),('V',V)]},
        'T_perp2':T.serialize(),'V_perp2':V.serialize(),'S_perp2_density':L.serialize(),
        'scope':'coefficient selection from serialized S<=4; not expansion of original untruncated action'}


class QuadraticAction:
    """L2=.5 v.T M v+v.T C eta+.5 eta.T Q eta, independently from polynomials."""
    def __init__(self,d,T,V):
        self.d=d;fields=rod.FIELD_ORDER
        M=[T.derivative(f+'_t').derivative(g+'_t') for f in fields for g in fields]
        C=[T.derivative(f+'_t').derivative(g) for f in fields for g in fields]
        M=[restrict_plane(x) for x in M];C=[restrict_plane(x) for x in C]
        Tqq=[restrict_plane(T.derivative(f).derivative(g)) for f in fields for g in fields]
        Vqq=[restrict_plane(V.derivative(f).derivative(g)) for f in _POT_NAMES for g in _POT_NAMES]
        kinetic=M+C+Tqq+[x.total_derivative('t') for x in M]+[x.total_derivative('t') for x in C]
        self.compiled=[];self.B=[];self.P=[]
        for arm,map0 in zip(d.arms,d.maps):
            self.compiled.append((_SparseCompiled(kinetic,_JET_NAMES,arm.coefficients),_SparseCompiled(Vqq,_JET_NAMES,arm.coefficients)))
            self.B.append([arm.B[i]@map0[arm.slices[f],:][:,d.out] for i,f in enumerate(fields)])
            self.P.append([(arm.D if deriv else arm.B)[i]@map0[arm.slices[fields[i]],:][:,d.out] for i,deriv in zip(_POT_FIELDS,_POT_DERIV)])

    @staticmethod
    def pullback(values,basis,weights):
        n=basis[0].shape[1];result=np.zeros((n,n))
        for i,left in enumerate(basis):
            for j,right in enumerate(basis):
                if np.any(values[i,j]):result+=left.T@((weights*values[i,j])[:,None]*right)
        return result

    def blocks(self,q,v,acceleration=None):
        if acceleration is None:acceleration=np.zeros(self.d.ndof)
        totals=[np.zeros((len(self.d.out),len(self.d.out))) for _ in range(5)]
        for arm,map0,(kinetic,potential),B,P in zip(self.d.arms,self.d.maps,self.compiled,self.B,self.P):
            qi=map0@q;vi=map0@v;ai=map0@acceleration
            variables=np.vstack([arm.reconstruct(qi).T,arm.reconstruct(qi,derivative=1).T,arm.reconstruct(vi).T,
                arm.reconstruct(qi,derivative=2).T,arm.reconstruct(vi,derivative=1).T,arm.reconstruct(ai).T])
            values=kinetic.evaluate(variables).reshape(5,7,7,arm.nq)
            for i in range(5):totals[i]+=self.pullback(values[i],B,arm.weights)
            totals[2]-=self.pullback(potential.evaluate(variables).reshape(11,11,arm.nq),P,arm.weights)
        return dict(zip(('M','C','Q','Mdot','Cdot'),totals))


def action_operators(blocks):
    M,C,Q,Mdot,Cdot=(blocks[k] for k in ('M','C','Q','Mdot','Cdot'));n=len(M)
    inv=cho_solve(cho_factor(M,lower=True,check_finite=False),np.eye(n),check_finite=False)
    A=np.block([[np.zeros_like(M),np.eye(n)],[-inv@(Cdot-Q),-inv@(Mdot+C-C.T)]])
    transform=np.block([[np.eye(n),np.zeros_like(M)],[C,M]])
    derivative=np.block([[np.zeros_like(M),np.zeros_like(M)],[Cdot,Mdot]])
    Hessian=np.block([[C.T@inv@C-Q,-C.T@inv],[-inv@C,inv]])
    J=np.block([[np.zeros_like(M),np.eye(n)],[-np.eye(n),np.zeros_like(M)]])
    return A,J@Hessian,Hessian,transform,derivative,J


def pair_matrix(h,k,delta,hr=0j,hs=0j):
    return np.array([[hr-1j*delta/2,h],[k,hs+1j*delta/2]],complex)


def selected_pair(row):
    h=complex(row['coupling_forward']['real'],row['coupling_forward']['imag'])
    k=complex(row['coupling_reverse']['real'],row['coupling_reverse']['imag'])
    hr,hs=[complex(v['real'],v['imag']) for v in row['mean_diagonal']]
    B=pair_matrix(h,k,row['signed_detuning'],hr,hs);eigen=np.linalg.eigvals(B)
    trace=np.trace(B)/2;discriminant=h*k+((B[0,0]-B[1,1])/2)**2
    W=np.diag([1/row['omega_i'],-1/row['omega_j']])
    ham=float(np.linalg.norm(B.conj().T@W+W@B)/max(np.linalg.norm(W@B),1e-300))
    shifted=row['signed_detuning']-(hr-hs).imag
    g=np.sqrt(abs(h*k));formal=pair_matrix(h,k,(hr-hs).imag,hr,hs)
    return {'source_row':row,'h':s.complex_row(h),'k':s.complex_row(k),'product':s.complex_row(h*k),
        'g':float(g),'effective_detuning':float(shifted),'eigenvalues':[s.complex_row(x) for x in eigen],
        'trace_half':s.complex_row(trace),'discriminant':s.complex_row(discriminant),
        'pseudo_Hermitian_Hamiltonian_relative_defect':ham,
        'formal_effective_detuning_zero_eigenvalues':[s.complex_row(x) for x in np.linalg.eigvals(formal)],
        'formal_resonance_not_realized_geometry':True}


def canonical_saved_sum(b,model,T,V,config,summary,deadline):
    """One historical C sum pair in momenta; raw velocity rows are not canonical."""
    source=p.ROOT/config['modal_parent']
    with np.load(source/'C_modes.npz') as z:
        q0=z['planar_initial_shape'].copy();R=z['right'][:,[0,1,5,6]].copy();L=z['left'][[0,1,5,6]].copy()
        savedK=z['K0'].copy();savedM=z['M0'].copy();omega=z['omega'][:2].copy()
    meta=p.read(source/'C_mode_metadata.json');Omega=meta['Omega']
    pc={**p.read(s.c.CONFIG),'lengths':[.8,1.2],'beta_deg':30.};d=p.make_disc(pc,'C',12,model)
    np.testing.assert_allclose(d.K,savedK,rtol=0,atol=1e-13);np.testing.assert_allclose(d.M0,savedM,rtol=0,atol=1e-13)
    action=QuadraticAction(d,T,V);zero=action_operators(action.blocks(np.zeros(d.ndof),np.zeros(d.ndof)))
    summary['calls']['canonical_action_evaluations']+=1
    S0=zero[3];Rc=S0@R;J=zero[5];Sigma=Rc.conj().T@J@Rc
    # The symplectic left dual avoids treating velocity coordinates as momenta.
    Lc=np.linalg.solve(Sigma,Rc.conj().T@J)
    dual_difference=float(np.max(abs(Lc-L@np.linalg.inv(S0))))
    if dual_difference>config['modal_duality_gate']:raise ArithmeticError('Historical canonical duality failed')
    phases=np.arange(config['canonical_phase_samples'])*2*np.pi/config['canonical_phase_samples'];orders=config['harmonics'];H=[];errors=[]
    for factor in config['amplitude_factors']:
        values=[]
        for phase in phases:
            if time.perf_counter()>deadline:raise TimeoutError('Selected canonical sum budget')
            blocks=action.blocks(q0*factor*np.cos(phase),-Omega*q0*factor*np.sin(phase))
            Hessian=action_operators(blocks)[2]
            B=-np.linalg.solve(Sigma,Rc.conj().T@(Hessian-zero[2])@Rc);values.append(B)
            errors.append(float(np.linalg.norm(B.conj().T@Sigma+Sigma@B)/max(np.linalg.norm(Sigma@B),1e-300)))
            summary['calls']['canonical_action_evaluations']+=1
        H.append(s.harmonics(np.array(values),phases,orders))
    H=np.array(H);amplitudes=config['amplitude']*np.array(config['amplitude_factors'])
    leading,next_leading,powers=s.amplitude_orders(H,amplitudes,orders);result=leading*amplitudes[0]**powers[:,None,None]
    h=result[orders.index(1),0,3];k=result[orders.index(-1),3,0]
    weight=np.diag([1/omega[0],-1/omega[1]]);small=pair_matrix(h,k,0)
    ham=float(np.linalg.norm(small.conj().T@weight+weight@small)/np.linalg.norm(weight@small))
    if max(errors)>config['Hamiltonian_relative_gate'] or ham>config['Hamiltonian_relative_gate'] or (h*k).real<=0:
        raise ArithmeticError('Canonical sum-pair Hamiltonian signs failed')
    diagnostic={'geometry':'historical C beta30 mu0.2, no search or new eigenpairs','source':config['modal_parent'],
        'canonical_h':s.complex_row(h),'canonical_k':s.complex_row(k),'canonical_product':s.complex_row(h*k),
        'canonical_g':float(np.sqrt(abs(h*k))),'projected_Hamiltonian_relative_defect':max(errors),
        'pair_Hamiltonian_relative_defect':ham,'symplectic_vs_energy_left_dual_difference':dual_difference,
        'amplitude_order_h':float(np.log2(abs(H[0,orders.index(1),0,3])/abs(H[1,orders.index(1),0,3]))),
        'quartet_required_for_real_solution':True,'canonical_state_signs':[1,2,-1,-2],
        'formal_zero_detuning_eigenvalues':[s.complex_row(v) for v in np.linalg.eigvals(small)],
        'not_an_actual_geometry_instability_result':True}
    p.write(b/'canonical_sum_channel_check.json',diagnostic)
    np.savez_compressed(b/'canonical_sum_channel_coefficients.npz',right=Rc,left=Lc,symscale=Sigma,
        full=H,leading_at_full_amplitude=result,orders=orders,amplitudes=amplitudes)
    return diagnostic


def validate_scope(config):
    if config!=p.read(CONFIG) or config['authorization']!=AUTH:raise ValueError('Only the frozen feasibility audit is authorized')
    if (config['degree']!=12 or config['beta_deg']!=15. or config['mu']!=.1 or config['amplitude']!=.005
        or config['numerical_budget_seconds']>180 or config['new_ODE_integrations'] or config['new_eigenproblems']
        or config['new_FEM_jobs'] or config['automatic_retry']):raise ValueError('Feasibility scope expansion')
    return config


def coordinate_comparison(b):
    """Saved direct coefficients in two representations, not a complete normal form."""
    saved=p.read(b/'beta15_margin_diagnostics.json')
    with np.load(b/'canonical_channel_coefficients.npz') as z:H=z['leading_at_full_amplitude'];orders=z['orders'].tolist()
    h=H[orders.index(1),0,1];k=H[orders.index(-1),1,0];hr=H[orders.index(0),0,0];hs=H[orders.index(0),1,1]
    delta=saved['Omega']-2*saved['omega_OUT2'];effective=delta-(hr-hs).imag;g=float(np.sqrt(abs(h*k)))
    pair=pair_matrix(h,k,delta,hr,hs)
    result={'canonical_h':s.complex_row(h),'canonical_k':s.complex_row(k),'canonical_mean_diagonal':[s.complex_row(hr),s.complex_row(hs)],
        'canonical_direct_g':g,'canonical_direct_effective_detuning':float(effective),
        'canonical_direct_ratio':float(2*g/abs(effective)),'canonical_direct_eigenvalues':[s.complex_row(v) for v in np.linalg.eigvals(pair)],
        'canonical_vs_saved_g_difference':g-saved['g'],'canonical_vs_saved_mean_detuning_difference':float(effective-saved['effective_detuning']),
        'canonical_direct_margin':float(abs(effective)-2*g),
        'detuning_change_as_main_frequency_fraction_if_other_terms_fixed':saved['required_detuning_change_toward_zero']/saved['Omega'],
        'interpretation':'direct O(a^2) mean coefficients are not a coordinate-invariant second-order normal form; the momentum map and off-orbit acceleration matter',
        'old_beta15_values_not_replaced':True,'complete_second_order_averaging':'NOT_RUN'}
    p.write(b/'canonical_coordinate_effects.json',result);return result


def verify_sources(config):
    for key in ('beta_parent','modal_parent'):
        if p.sha(p.ROOT/config[key]/'manifest.json')!=config[key+'_manifest_sha256']:raise ValueError('Immutable parent changed')
    s.beta_cache(p.ROOT/config['beta_parent'])


def prepare(path=CONFIG):
    config=validate_scope(p.read(path));verify_sources(config)
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        item=p.read(old)
        if item['config']['authorization']==AUTH:
            if item['config']!=config:raise ValueError('Authorization already used with different settings')
            return old.parent,item,p.cache(old.parent)
    parent=p.read(p.ROOT/config['beta_parent']/'provenance.json')
    dependencies=['scripts/lib/nlsp_instability_feasibility.py','scripts/lib/nlsp_resonant_spatial_coupling_screening.py',
        'scripts/lib/nlsp_construction_c_sensitivity.py','data/input/nlsp_construction_C_modal_sensitivity.json',CONFIG.relative_to(p.ROOT).as_posix()]
    hashes={**p.code_hashes(),**{f:p.sha(p.ROOT/f) for f in dependencies}}
    item={'schema':'nlsp-instability-feasibility-v1','config':config,'physical_source_config':parent['physical_source_config'],
        'source_artifacts':parent['source_artifacts'],'source_manifests':{**parent['source_manifests'],
            'beta_parent':{'path':config['beta_parent'],'sha256':config['beta_parent_manifest_sha256']}},
        'code_hashes':hashes,'initial_state':p.read(SMOKE/'initial_state.json'),'interpreter':sys.executable,
        'no_ODE_eigenproblem_Floquet_or_action_rederivation':True}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp
    if b.exists():raise ValueError('Preserve interrupted feasibility attempt')
    b.mkdir(parents=True);p.write(b/'config.json',config);p.write(b/'provenance.json',item);folder=b/'execution_code';folder.mkdir()
    for f in hashes:shutil.copyfile(p.ROOT/f,folder/Path(f).name)
    summary={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':0.,
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'eigenproblems':0,'CCX':0,'Gmsh':0,'pointwise_Jacobian':0,'canonical_action_evaluations':0},
        'stages':{'source':'PASS','A':'NOT_RUN','B':'NOT_RUN','C':'NOT_RUN'}}
    p.finalize(b,item,summary);return b,item,summary


def cache(b):verify_sources(p.read(Path(b)/'config.json'));return p.cache(b)


def run(b,item,summary):
    if summary['completed'] or summary['hard_stop']:return summary
    started=time.perf_counter();config=item['config'];deadline=started+config['numerical_budget_seconds']
    try:
        model=p.sources(item['physical_source_config'])[0];T,V,identity=second_variation(model);p.write(b/'quadratic_action_polynomials.json',identity)
        pc={**p.read(s.c.CONFIG),'lengths':config['lengths'],'beta_deg':config['beta_deg']};d=p.make_disc(pc,'C',12,model)
        source=p.ROOT/config['beta_parent'];folder=source/'beta_15'
        with np.load(folder/'C_modes.npz') as z:
            q0=z['planar_initial_shape'].copy();R=z['right'][:,[1,6]].copy();L=z['left'][[1,6]].copy();metric=z['metric'].copy()
            phi=z['phi'][:,1].copy();omega=float(z['omega'][1]);K0=z['K0'].copy();M0=z['M0'].copy()
        meta=p.read(folder/'C_mode_metadata.json');Omega=meta['Omega'];operator=QuadraticAction(d,T,V)
        np.testing.assert_allclose(d.K,K0,rtol=0,atol=1e-13);np.testing.assert_allclose(d.M0,M0,rtol=0,atol=1e-13)
        checks=[];arrays={};zero_transform=None
        for j,phase in enumerate([None]+config['point_phases']):
            if time.perf_counter()>deadline:raise TimeoutError('180-second feasibility ceiling')
            q=np.zeros(d.ndof) if phase is None else q0*np.cos(phase)
            velocity=np.zeros(d.ndof) if phase is None else -Omega*q0*np.sin(phase)
            acc=d.acceleration(q,velocity);blocks=operator.blocks(q,velocity,acc)
            A,canonical,H,S,Sdot,J=action_operators(blocks);jac=d.jacobian(0,np.r_[q,velocity])[np.ix_(d.state_out,d.state_out)]
            summary['calls']['pointwise_Jacobian']+=1;summary['calls']['canonical_action_evaluations']+=1
            Sinv=np.linalg.inv(S);converted=S@jac@Sinv+Sdot@Sinv
            action_error=float(np.linalg.norm(A-jac)/np.linalg.norm(jac));can_error=float(np.linalg.norm(converted-canonical)/np.linalg.norm(canonical))
            ham=float(np.linalg.norm(canonical.T@J+J@canonical)/np.linalg.norm(canonical))
            if max(action_error,can_error)>config['quadratic_action_Jacobian_relative_gate'] or ham>config['Hamiltonian_relative_gate']:
                raise ArithmeticError('Hard gate A: independent quadratic action/Jacobian mismatch')
            row={'phase':phase,'action_Jacobian_relative_difference':action_error,'canonical_transform_relative_difference':can_error,
                'Hamiltonian_relative_defect':ham,'Hessian_symmetry_relative':float(np.linalg.norm(H-H.T)/np.linalg.norm(H)),
                'minimum_transverse_mass_eigenvalue':float(np.linalg.eigvalsh(blocks['M']).min()),
                'background_acceleration':'full coupled RHS acceleration; not imposed harmonic ddq'}
            if phase is not None:
                harmonic=-Omega**2*q
                Aharm=action_operators(operator.blocks(q,velocity,harmonic))[0]
                summary['calls']['canonical_action_evaluations']+=1
                row.update(harmonic_acceleration_residual_relative=float(np.linalg.norm(acc-harmonic)/np.linalg.norm(acc)),
                    prescribed_harmonic_jet_action_vs_Jacobian_difference=float(np.linalg.norm(Aharm-jac)/np.linalg.norm(jac)))
            else:zero_transform=S
            checks.append(row)
            for key,value in {**blocks,'action_state_operator':A,'kernel_Jacobian':jac,'canonical_operator':canonical,'canonical_Hessian':H,'momentum_map':S}.items():arrays[f'point_{j}_{key}']=value
        summary['stages']['A']='PASS';p.write(b/'quadratic_action_Jacobian_checks.json',checks);np.savez_compressed(b/'quadratic_action_matrices.npz',**arrays)
        Sinv0=np.linalg.inv(zero_transform);Rc=zero_transform@R;Lc=L@Sinv0
        J=np.block([[np.zeros((len(d.out),len(d.out))),np.eye(len(d.out))],[-np.eye(len(d.out)),np.zeros((len(d.out),len(d.out)))]])
        Sigma=Rc.conj().T@J@Rc
        dual=float(np.max(abs(Lc@Rc-np.eye(2))));symdual=float(np.max(abs(Lc-np.linalg.solve(Sigma,Rc.conj().T@J))))
        if dual>config['modal_duality_gate'] or symdual>config['modal_duality_gate']:raise ArithmeticError('Canonical modal duality failed')
        base=action_operators(operator.blocks(np.zeros(d.ndof),np.zeros(d.ndof)))[1];amplitudes=config['amplitude']*np.array(config['amplitude_factors'])
        summary['calls']['canonical_action_evaluations']+=1
        phase=np.arange(config['canonical_phase_samples'])*2*np.pi/config['canonical_phase_samples'];harmonics=[];coarse=[];Hamiltonian_errors=[]
        for factor in config['amplitude_factors']:
            values=[]
            for angle in phase:
                if time.perf_counter()>deadline:raise TimeoutError('Canonical coefficient budget')
                stateq=q0*factor*np.cos(angle);statev=-Omega*q0*factor*np.sin(angle)
                Ac=action_operators(operator.blocks(stateq,statev))[1];B=Lc@(Ac-base)@Rc;values.append(B)
                Hamiltonian_errors.append(float(np.linalg.norm(B.conj().T@Sigma+Sigma@B)/max(np.linalg.norm(Sigma@B),1e-300)))
                summary['calls']['canonical_action_evaluations']+=1
            values=np.array(values);harmonics.append(s.harmonics(values,phase,config['harmonics']))
            coarse.append(s.harmonics(values[::2],phase[::2],config['harmonics']))
        harmonics=np.array(harmonics);leading,next_leading,powers=s.amplitude_orders(harmonics,amplitudes,config['harmonics'])
        Hlead=leading*amplitudes[0]**powers[:,None,None];order=config['harmonics'];h=Hlead[order.index(1),0,1];k=Hlead[order.index(-1),1,0]
        can_g=float(np.sqrt(abs(h*k)));order_h=float(np.log2(abs(harmonics[0,order.index(1),0,1])/abs(harmonics[1,order.index(1),0,1])))
        if max(Hamiltonian_errors)>config['Hamiltonian_relative_gate'] or abs(h-k.conjugate())>1e-10 or (h*k).real<=0 or abs(order_h-1)>.1:
            raise ArithmeticError('Canonical leading growth channel not confirmed')
        canonical_diag={'canonical_duality_max_difference':dual,'symplectic_dual_max_difference':symdual,
            'projected_Hamiltonian_relative_defect':max(Hamiltonian_errors),'canonical_leading_h':s.complex_row(h),
            'canonical_leading_k':s.complex_row(k),'canonical_leading_g':can_g,'amplitude_order':order_h,
            'phase_16_32_difference':float(np.max(abs(harmonics-np.array(coarse)))),
            'Richardson_difference_at_full_amplitude':float(np.max(abs((leading-next_leading)*amplitudes[0]**powers[:,None,None]))),
            'frequency':omega,'Omega':Omega,'normal_form_second_order':'NOT_RUN; canonical and velocity direct coefficients may differ by transformation terms'}
        p.write(b/'canonical_channel_check.json',canonical_diag);np.savez_compressed(b/'canonical_channel_coefficients.npz',right=Rc,left=Lc,symscale=Sigma,
            full=harmonics,leading_at_full_amplitude=Hlead,orders=order,amplitudes=amplitudes)
        rows=p.read(folder/'interaction_rows.json');row=next(r for r in rows if (r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])==('self',2,2,1))
        selected=selected_pair(row)
        if selected['pseudo_Hermitian_Hamiltonian_relative_defect']>1e-10:raise ArithmeticError('Saved self-pair Hamiltonian signs failed')
        g=selected['g'];delta=selected['effective_detuning'];margin=abs(delta)-2*g
        with np.load(folder/'C_harmonic_coefficients.npz') as z:
            uncertainty=z['numerical_difference'];full=z['full'];errors=uncertainty[order.index(1),1,6],uncertainty[order.index(-1),6,1]
            detuning_difference=float(uncertainty[order.index(0),1,1]+uncertainty[order.index(0),6,6])
        g_difference=float(max(errors));margin_difference=detuning_difference+2*g_difference
        physical={'Omega':Omega,'omega_OUT2':omega,'selected_PL_index':meta['planar']['selected_sorted_index_one_based'],
            'selected_PL_bending_shear_fraction':meta['planar']['modes'][meta['planar']['selected_sorted_index_one_based']-1]['bending_shear_fraction'],
            'max_w':p.polynomial_w_max(d,q0)[0],'planar_eigenpair_relative_residual':float(np.linalg.norm(d.K@q0-Omega**2*d.M0@q0)/np.linalg.norm(d.K@q0)),
            'joint':d.joint_diagnostics(q0,np.zeros(d.ndof)),'no_beta30_branch_identity_asserted':True,
            'g':g,'effective_detuning':delta,'two_g_over_abs_effective_detuning':2*g/abs(delta),
            'sigma_squared':g*g-delta*delta/4,'margin_abs_delta_minus_two_g':margin,
            'required_detuning_change_toward_zero':margin,'required_relative_detuning_change':margin/abs(delta),
            'required_g_increase':abs(delta)/2-g,'required_relative_g_increase':abs(delta)/(2*g)-1,
            'ideal_tuning_e_folding_time':1/g,'ideal_tuning_e_folding_reference_periods':Omega/(2*np.pi*g),
            'observed_g_numerical_difference':g_difference,'observed_effective_detuning_difference':detuning_difference,
            'observed_margin_difference':margin_difference,'margin_over_observed_difference':margin/margin_difference,
            'numerical_difference_is_not_rigorous_error_bound':True,'nonlinear_orbit_frequency_and_shape_corrections':'NOT_AVAILABLE at beta15'}
        p.write(b/'saved_beta15_pair.json',selected);p.write(b/'beta15_margin_diagnostics.json',physical)
        coordinate_comparison(b)
        # One historical sum pair, read only; do not re-screen the 250 channels.
        modal_source=p.ROOT/config['modal_parent'];r0=next(r for r in p.read(modal_source/'C_interaction_rows.json') if (r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])==('sum',1,2,1))
        with np.load(modal_source/'C_harmonic_coefficients.npz') as z:H=z['leading_at_full_amplitude']
        r0['mean_diagonal']=[s.complex_row(H[order.index(0),0,0]),s.complex_row(H[order.index(0),6,6])]
        sum_pair=selected_pair(r0);p.write(b/'saved_sum_pair.json',sum_pair)
        # Historical rows were projected in (eta,velocity); Hamiltonian
        # structure must be checked after the time-dependent momentum map.
        canonical_saved_sum(b,model,T,V,config,summary,deadline)
        summary['stages'].update(B='PASS',C='PASS');summary.update(completed=True,overall='STRUCTURAL_PARAMETRIC_GROWTH_CHANNEL_CONFIRMED',
            physical_geometry_with_instability_proved=False,true_periodic_orbit_instability_proved=False,new_geometry_candidate_selected=False)
    except Exception as error:
        summary.update(hard_stop=True,overall='NUMERICAL_PARTIAL' if isinstance(error,TimeoutError) else 'BLOCKED',stop_reason=type(error).__name__+': '+str(error))
    summary['numerical_seconds']=time.perf_counter()-started;p.finalize(b,item,summary)
    if summary['completed']:plot_saved(b);p.finalize(b,item,summary)
    return summary


def plot_saved(b):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data=p.read(b/'beta15_margin_diagnostics.json');g=data['g'];delta=data['effective_detuning']
    x=np.linspace(-.012,.012,401);y=np.linspace(0,.006,301);X,Y=np.meshgrid(x,y)
    fig,ax=plt.subplots(figsize=(8,5));ax.contourf(X,Y,(Y*Y-X*X/4)>0,levels=[-.5,.5,1.5],colors=['#edf2fa','#fff0df'])
    ax.plot(x,abs(x)/2,color='black',label='pair boundary: 2g = |effective detuning|')
    ax.plot(delta,g,'o',color='tab:blue',label='saved beta15 self2/Ω (no growth)');ax.plot(0,g,'x',ms=9,color='tab:red',label='formal tuning only; not a realized geometry')
    ax.set(xlabel='effective detuning [1/time]',ylabel='g = sqrt(|h k|) [1/time]',ylim=(0,.006),title='Leading self-pair model; no geometry/amplitude map')
    ax.text(0,.0055,'formal pair growth region',ha='center')
    ax.annotate('saved ratio = 0.9662',xy=(delta,g),xytext=(-.011,.003),arrowprops={'arrowstyle':'->'},fontsize=8)
    ax.grid(alpha=.2);ax.legend(fontsize=8,loc='lower left');fig.tight_layout()
    names=[]
    for ext in ('pdf','png'):name='pair_growth_feasibility.'+ext;fig.savefig(b/name,dpi=180,bbox_inches='tight');names.append(name)
    plt.close(fig);p.write(b/'figure_provenance.json',{'figures':names,'no_scientific_calls':True,'diagram':'algebraic pair-parameter plane, not physical geometry or amplitude map'})


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();summary=cache(b)
        if args.plot_only:plot_saved(b);p.finalize(b,p.read(b/'provenance.json'),summary)
    else:b,item,summary=prepare(args.config);summary=run(b,item,summary)
    print(json.dumps(p.plain({'bundle':str(b),'overall':summary['overall'],'stages':summary['stages'],'calls':summary['calls'],
        'seconds':summary['numerical_seconds'],'stop_reason':summary.get('stop_reason')}),indent=2))
