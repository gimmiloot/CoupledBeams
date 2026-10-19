"""Fixed beta15 PL expansion and (only after its gate) canonical pair reduction.

Reads the existing quartic action. No time integrator or model modification.
The physical maximum initial w and the fundamental modal coefficient differ.
"""
from __future__ import annotations
import argparse, hashlib, json, shutil, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import solve, eigh
from scipy.optimize import brentq
from scripts.lib import nlsp_instability_feasibility as f
from scripts.lib.weakly_nonlinear_spatial_dynamics import _SparseCompiled, _JET_NAMES, _POT_NAMES, _POT_FIELDS, _POT_DERIV

p=f.p; s=f.s; rod=f.rod
CONFIG=p.ROOT/'data/input/nlsp_beta15_nonlinear_planar_resonance.json'
OUTPUT=p.ROOT/'results/nlsp_beta15_nonlinear_planar_resonance'
SMOKE=p.ROOT/'results/_smoke/nlsp_beta15_nonlinear_planar_resonance'


class HomogeneousMassForm:
    """Weak EL of T2..4/V2..4, with original arm maps and quadrature.

    No inversion of M and no Taylor subtraction of large linear residuals.
    This is coefficient selection from the accepted polynomial action.
    """
    def __init__(self,d,model):
        self.d=d; self.parts={}
        for degree in (2,3,4):
            T=model.T4.homogeneous(degree); V=model.V4.homogeneous(degree)
            kinetic=[]
            for field in rod.FIELD_ORDER:
                mom=T.derivative(field+'_t')
                kinetic.append(mom.total_derivative('t')-T.derivative(field))
            potential=[V.derivative(x) for x in _POT_NAMES]
            self.parts[degree-1]=[_SparseCompiled(kinetic+potential,_JET_NAMES,a.coefficients) for a in d.arms]
        self.calls=0

    def residual(self,degree,q,v,acc):
        self.calls+=1; result=np.zeros(self.d.ndof)
        for arm,P,compiled in zip(self.d.arms,self.d.maps,self.parts[degree]):
            qi,vi,ai=P@q,P@v,P@acc
            variables=np.vstack([arm.reconstruct(qi).T,arm.reconstruct(qi,derivative=1).T,
                arm.reconstruct(vi).T,arm.reconstruct(qi,derivative=2).T,
                arm.reconstruct(vi,derivative=1).T,arm.reconstruct(ai).T])
            values=compiled.evaluate(variables); local=np.zeros(arm.ndof)
            for i,field in enumerate(rod.FIELD_ORDER):
                local[arm.slices[field]]+=arm.B[i].T@(arm.weights*values[i])
            for j,(i,deriv) in enumerate(zip(_POT_FIELDS,_POT_DERIV)):
                local[arm.slices[rod.FIELD_ORDER[i]]]+=(arm.D if deriv else arm.B)[i].T@(arm.weights*values[7+j])
            result+=P.T@local
        return result


def cosine_coefficients(values,phases,n):
    return np.array([np.mean(values,axis=0) if j==0 else
        2*np.mean(values*np.cos(j*phases)[:,None],axis=0) for j in range(n+1)])


def orbit(coeff,a,phase):
    """Mass-coordinate q,v,ddq of the third-order formal PL orbit."""
    om=coeff['Omega0']+coeff['Omega2']*a*a
    q=a*coeff['phi']*np.cos(phase)+a*a*(coeff['q20']+coeff['q22']*np.cos(2*phase))
    q+=a**3*(coeff['q31']*np.cos(phase)+coeff['q33']*np.cos(3*phase))
    qt=-a*coeff['phi']*np.sin(phase)-2*a*a*coeff['q22']*np.sin(2*phase)
    qt-=a**3*(coeff['q31']*np.sin(phase)+3*coeff['q33']*np.sin(3*phase))
    qtt=-a*coeff['phi']*np.cos(phase)-4*a*a*coeff['q22']*np.cos(2*phase)
    qtt-=a**3*(coeff['q31']*np.cos(phase)+9*coeff['q33']*np.cos(3*phase))
    return q,om*qt,om*om*qtt


def pl_coefficients(d,massform,q0,Omega,cfg):
    ids=d.plane; M=d.M0[np.ix_(ids,ids)]; K=d.K[np.ix_(ids,ids)]
    phi=q0/cfg['reference_amplitude']; ph=phi[ids]; phases=np.arange(cfg['PL_phases'])*2*np.pi/cfg['PL_phases']
    triplets=[(phi*np.cos(t),-Omega*phi*np.sin(t),-Omega**2*phi*np.cos(t)) for t in phases]
    R2=np.array([massform.residual(2,*x)[ids] for x in triplets]); F2=cosine_coefficients(R2,phases,3)
    q20=solve(K,-F2[0],assume_a='sym'); q22=solve(K-4*Omega**2*M,-F2[2],assume_a='sym')
    Omega1=float(ph@F2[1]/(2*Omega*(ph@M@ph)))
    def full(x):
        y=np.zeros(d.ndof);y[ids]=x;return y
    R3=[]
    for t,x in zip(phases,triplets):
        y=(full(q20+q22*np.cos(2*t)),full(-2*Omega*q22*np.sin(2*t)),full(-4*Omega**2*q22*np.cos(2*t)))
        # Exact polarization of the quadratic mass-form residual.
        cross=massform.residual(2,*[a+b for a,b in zip(x,y)])-massform.residual(2,*x)-massform.residual(2,*y)
        R3.append((cross+massform.residual(3,*x))[ids])
    F3=cosine_coefficients(np.array(R3),phases,3)
    Omega2=float(ph@F3[1]/(2*Omega*(ph@M@ph)))
    bordered=np.block([[K-Omega**2*M,(M@ph)[:,None]],[(ph@M)[None,:],np.zeros((1,1))]])
    rhs=-F3[1]+2*Omega*Omega2*M@ph
    q31=solve(bordered,np.r_[rhs,0.],assume_a='sym')[:-1]
    q33=solve(K-9*Omega**2*M,-F3[3],assume_a='sym')
    coeff={'Omega0':Omega,'Omega1':Omega1,'Omega2':Omega2,'phi':phi,
        'q20':full(q20),'q22':full(q22),'q31':full(q31),'q33':full(q33)}
    # Full fixed-p12 diagonalization diagnoses algebraic harmonic denominators;
    # it is not a geometry or frequency/root scan and computes no trajectory.
    val,vec=eigh(K,M,check_finite=False); omega=np.sqrt(val)
    denominators=[]
    for harmonic,force,response in [(0,F2[0],q20),(2,F2[2],q22),(3,F3[3],q33)]:
        den=val-harmonic**2*Omega**2; j=int(np.argmin(abs(den)))
        overlap=abs(vec[:,j]@force)/max(np.linalg.norm(vec.T@force),1e-300)
        energy_ratio=cfg['amplitudes'][-1]**(1 if harmonic<3 else 2)*abs(vec[:,j]@M@response)*omega[j]/np.sqrt(ph@K@ph)
        denominators.append({'harmonic':harmonic,'sorted_mode':j+1,'mode_omega':omega[j],
            'squared_frequency_denominator':den[j],
            'relative_frequency_gap':abs(omega[j]-harmonic*Omega)/max(omega[j],Omega),
            'forcing_modal_fraction':overlap,'correction_energy_norm_ratio_at_Amax':energy_ratio,
            'spectral_condition':float(max(abs(den))/min(abs(den)))})
    checks={'Omega1':Omega1,'Omega2':Omega2,'frequency_solvability':float(abs(ph@rhs)/max(np.linalg.norm(ph)*np.linalg.norm(rhs),1e-300)),
        'fundamental_gauge':float(abs(ph@M@q31)),
        'forcing_harmonic_leak_R2':float(np.linalg.norm(F2[[1,3]])/np.linalg.norm(F2)),
        'forcing_harmonic_leak_R3':float(np.linalg.norm(F3[[0,2]])/np.linalg.norm(F3)),
        'denominators':denominators,'full_planar_omega':omega.tolist(),
        'harmonic_solve_relative_residuals':[float(np.linalg.norm(A@x+F)/max(np.linalg.norm(F),1e-300))
            for A,x,F in [(K,q20,F2[0]),(K-4*Omega**2*M,q22,F2[2]),(K-9*Omega**2*M,q33,F3[3])]]}
    return coeff,checks,{'phases':phases,'F2':F2,'F3':F3,'planar_omega':omega,'planar_phi':vec}


def amplitude_parameter(d,coeff,A):
    def size(a):return p.polynomial_w_max(d,orbit(coeff,a,0.)[0])[0]
    a=brentq(lambda x:size(x)-A,0.,1.5*A,xtol=1e-15)
    return a,size(a)


def validate_scope(cfg):
    fixed={'degree':12,'beta_deg':15.,'mu':.1,'lengths':[.9,1.1],
        'sections':[[.05,.1],[.05,.1]],'amplitudes':[.005,.006,.008,.010],
        'reference_amplitude':.005,'new_ODE_integrations':0,'new_FEM_jobs':0,
        'automatic_retry':False,'PL_phases':16,'canonical_phases':32,
        'action_massform_relative_gate':2e-10,'harmonic_solve_relative_gate':1e-8,
        'internal_resonance_relative_gap_trigger':.01,'internal_resonance_forcing_fraction_gate':1e-8,
        'internal_resonance_correction_energy_ratio_trigger':.05,
        'expansion_correction_energy_ratio_gate':.25,'residual_order_minimum':3.5,
        'Hamiltonian_relative_gate':2e-10,'numerical_budget_seconds':600.}
    if any(cfg.get(k)!=v for k,v in fixed.items()):raise ValueError('Frozen diagnostic scope or gates changed')
    return cfg


def prepare(config=CONFIG):
    cfg=validate_scope(p.read(config))
    if cfg!=p.read(CONFIG) or cfg['amplitudes']!=[.005,.006,.008,.010] or cfg['new_ODE_integrations']!=0:
        raise ValueError('Only the fixed four-amplitude algebraic continuation is authorized')
    for key in ('beta_parent','feasibility_parent'):
        source=p.ROOT/cfg[key]
        if p.sha(source/'manifest.json')!=cfg[key+'_manifest_sha256']:raise ValueError('Parent manifest mismatch')
        p.cache(source)
    for old in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        item=p.read(old)
        if item['config']==cfg:return old.parent,item,p.cache(old.parent)
    parent=p.read(p.ROOT/cfg['feasibility_parent']/'provenance.json')
    for file in f.s.r.KERNELS:
        if p.sha(p.ROOT/file)!=parent['code_hashes'][file]:raise ValueError('Frozen physical kernel changed')
    hashes={**p.code_hashes(),__file__:p.sha(__file__)}
    hashes={Path(k).relative_to(p.ROOT).as_posix() if Path(k).is_absolute() else k:v for k,v in hashes.items()}
    item={'schema':'nlsp-beta15-nonlinear-planar-resonance-v1','config':cfg,
        'source_artifacts':parent['source_artifacts'],
        'source_manifests':{**parent['source_manifests'],**{k:{'path':cfg[k],'sha256':cfg[k+'_manifest_sha256']} for k in ('beta_parent','feasibility_parent')}},
        'code_hashes':hashes,'physical_source_config':parent['physical_source_config'],
        'initial_state':p.read(SMOKE/'initial_state.json'),'interpreter':sys.executable}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp
    if b.exists():raise ValueError('Interrupted bundle must be preserved')
    b.mkdir(parents=True);p.write(b/'provenance.json',item);p.write(b/'config.json',cfg)
    (b/'execution_code').mkdir()
    for file in hashes:shutil.copyfile(p.ROOT/file,b/'execution_code'/Path(file).name)
    summary={'overall':'PREPARED','completed':False,'numerical_seconds':0.,'calls':{'ODE':0,'FEM':0,'pointwise_Jacobian':0,'algebraic_planar_eigendecompositions':0},
        'stages':{'source':'PASS','B':'NOT_RUN','C':'NOT_RUN'}}
    p.finalize(b,item,summary);return b,item,summary


def compute(b,item,summary):
    cfg=item['config'];start=time.perf_counter()
    model,_,_=p.sources(item['physical_source_config']);pc={**p.read(s.c.CONFIG),'lengths':[.9,1.1],'beta_deg':15.}
    d=p.make_disc(pc,'C',12,model);source=p.ROOT/cfg['beta_parent']/'beta_15'
    with np.load(source/'C_modes.npz') as z:q0=z['planar_initial_shape'].copy();M0=z['M0'].copy();K0=z['K0'].copy()
    np.testing.assert_allclose(d.K,K0,atol=1e-13,rtol=0);np.testing.assert_allclose(d.M0,M0,atol=1e-13,rtol=0)
    Omega=p.read(source/'C_mode_metadata.json')['Omega'];massform=HomogeneousMassForm(d,model)
    checks=[]
    for phase in [0.,.7,1.9]:
        q=q0*np.cos(phase);v=-Omega*q0*np.sin(phase);acc=-Omega**2*q
        assembled=sum(massform.residual(n,q,v,acc) for n in (1,2,3))
        current=d.mass_matrix(q)@acc+d.inertial_terms(q,v)+d.potential(q)['gradient']
        checks.append(float(np.linalg.norm(assembled-current)/max(np.linalg.norm(current),1e-300)))
    if max(checks)>cfg['action_massform_relative_gate']:raise ArithmeticError('Independent action mass-form mismatch')
    coeff,diagnostic,raw=pl_coefficients(d,massform,q0,Omega,cfg)
    summary['calls']['algebraic_planar_eigendecompositions']+=1
    diagnostic['action_massform_relative_differences']=checks
    if max(diagnostic['harmonic_solve_relative_residuals'])>cfg['harmonic_solve_relative_gate']:
        raise ArithmeticError('Harmonic algebraic solve residual failed')
    np.savez_compressed(b/'PL_coefficients.npz',**coeff,**raw)
    p.write(b/'PL_diagnostics.json',diagnostic)
    summary['numerical_seconds']=time.perf_counter()-start
    summary['calls']['homogeneous_massform_evaluations']=massform.calls
    summary['overall']='PL_COEFFICIENTS_SAVED';p.finalize(b,item,summary)
    return d,massform,coeff,diagnostic


def baseline_rows(cfg):
    """Read-only M0 extrapolation: old velocity and canonical direct kept apart."""
    parent=p.ROOT/cfg['feasibility_parent'];old=p.read(parent/'saved_beta15_pair.json')
    can=p.read(parent/'canonical_coordinate_effects.json');delta=old['source_row']['signed_detuning'];rows=[]
    for A in cfg['amplitudes']:
        row={'A':A,'approximation':'historical harmonic/direct amplitude extrapolation only','definitions':{}}
        for name,g0,de0 in [('velocity_direct',old['g'],old['effective_detuning']),
                            ('canonical_direct',can['canonical_direct_g'],can['canonical_direct_effective_detuning'])]:
            g=g0*A/cfg['reference_amplitude'];de=delta+(de0-delta)*(A/cfg['reference_amplitude'])**2
            discriminant=g*g-de*de/4;sigma=float(np.sqrt(max(discriminant,0.)))
            row['definitions'][name]={'g':g,'effective_detuning':de,'two_g_over_absolute_detuning':2*g/abs(de),
                'discriminant':discriminant,'sigma':sigma,'margin':2*g-abs(de)}
        rows.append(row)
    return rows


def cache_saved(b):
    summary=p.cache(b);selection=Path(b)/'source_selection_provenance.json'
    if selection.exists():
        for source,hash0 in p.read(selection)['selected_parent_artifacts'].items():
            if p.sha(p.ROOT/source)!=hash0:raise ValueError('Selected parent array changed: '+source)
    return summary


def finish_preflight(b,item,summary,d,massform,coeff,diagnostic):
    """Only the residual/normalization audit of Stage B, never OUT averaging."""
    cfg=item['config'];start=time.perf_counter();deadline=start+cfg['numerical_budget_seconds']-summary['numerical_seconds']
    phases=np.arange(cfg['PL_phases'])*2*np.pi/cfg['PL_phases'];rows=[];arrays={}
    M=d.M0;K=d.K;ph=coeff['phi'];energy=lambda x:float(np.sqrt(x@K@x))
    spectral=diagnostic['denominators'][-1];small_denominator=(spectral['relative_frequency_gap']<cfg['internal_resonance_relative_gap_trigger']
        and spectral['forcing_modal_fraction']>cfg['internal_resonance_forcing_fraction_gate']
        and spectral['correction_energy_norm_ratio_at_Amax']>=cfg['internal_resonance_correction_energy_ratio_trigger'])
    if not small_denominator:
        raise ArithmeticError('Stage C implementation requires a passed, fully audited PL gate; no automatic scope expansion')
    for A in cfg['amplitudes']:
        if time.perf_counter()>deadline:raise TimeoutError('Stage B numerical budget')
        a,actual=amplitude_parameter(d,coeff,A);q,v,acc=orbit(coeff,a,0.)
        norms=[];energy_hist=[];field_res=[];states=[];velocities=[];accelerations=[]
        for fraction in (1.,.5):
            aa=a*fraction;residual=[];energies=[]
            for t in phases:
                x,y,z=orbit(coeff,aa,t);res=sum(massform.residual(n,x,y,z) for n in (1,2,3))
                residual.append(res);energies.append(d.energy(x,y))
                if fraction==1.:states.append(x);velocities.append(y);accelerations.append(z)
            norms.append(float(np.max(np.linalg.norm(residual,axis=1))))
            energy_hist.append(float((max(energies)-min(energies))/abs(np.mean(energies))))
            field_res.append(np.array(residual))
        ratio20=a*energy(coeff['q20'])/energy(ph);ratio22=a*energy(coeff['q22'])/energy(ph)
        ratio31=a*a*energy(coeff['q31'])/energy(ph);ratio33=a*a*energy(coeff['q33'])/energy(ph)
        harmonic_w=[p.polynomial_w_max(d,a*a*coeff[k])[0] for k in ('q20','q22')]
        harmonic_w += [p.polynomial_w_max(d,a**3*coeff[k])[0] for k in ('q31','q33')]
        mean_gap=spectral['mode_omega']-3*(coeff['Omega0']+coeff['Omega2']*a*a)
        row={'A':A,'fundamental_coefficient_a':a,'actual_initial_max_w':actual,'Omega_formal':coeff['Omega0']+coeff['Omega2']*a*a,
            'frequency_correction_formal':coeff['Omega2']*a*a,'frequency_accepted':False,
            'planar_mode4_minus_3Omega_formal':mean_gap,
            'initial_velocity_max':float(max(abs(v))),
            'initial_state_correction_energy_ratio':energy(q-a*ph)/(a*energy(ph)),
            'mean_energy_ratio':ratio20,'second_harmonic_energy_ratio':ratio22,
            'fundamental_cubic_energy_ratio':ratio31,'third_harmonic_energy_ratio':ratio33,
            'mean_max_w':harmonic_w[0],'second_harmonic_max_w':harmonic_w[1],
            'fundamental_cubic_max_w':harmonic_w[2],'third_harmonic_max_w':harmonic_w[3],
            'massform_residual_max':norms[0],'massform_residual_at_half_a':norms[1],
            'massform_residual_observed_order':float(np.log2(norms[0]/norms[1])),
            'relative_massform_residual_to_linear_inertial_scale':norms[0]/(a*coeff['Omega0']**2*np.linalg.norm(M@ph)),
            'energy_relative_variation':energy_hist[0],'energy_relative_variation_at_half_a':energy_hist[1],
            'energy_relative_variation_order':float(np.log2(energy_hist[0]/energy_hist[1])),
            'M1':'NOT_AVAILABLE','M2':'NOT_AVAILABLE','sigma_pair':'NOT_AVAILABLE',
            'classification':'PLANAR_SINGLE_MODE_EXPANSION_UNRESOLVED'}
        rows.append(row);label=str(A).replace('.','p')
        arrays['q_'+label]=np.array(states);arrays['v_'+label]=np.array(velocities);arrays['acc_'+label]=np.array(accelerations)
        arrays['residual_'+label]=field_res[0];arrays['residual_half_a_'+label]=field_res[1]
    np.savez_compressed(b/'PL_formal_residual_checks.npz',phases=phases,**arrays)
    p.write(b/'four_amplitude_diagnostics.json',rows);p.write(b/'M0_historical_extrapolation.json',baseline_rows(cfg))
    import csv
    with (b/'four_amplitude_diagnostics.csv').open('w',newline='',encoding='utf8') as stream:
        writer=csv.DictWriter(stream,fieldnames=rows[0].keys());writer.writeheader();writer.writerows(rows)
    crossing=np.sqrt((spectral['mode_omega']/3-coeff['Omega0'])/coeff['Omega2'])
    crossing_physical=p.polynomial_w_max(d,orbit(coeff,crossing,0.)[0])[0]
    linear_res=float(np.linalg.norm(K@ph-coeff['Omega0']**2*M@ph)/np.linalg.norm(K@ph))
    diagnostic.update({'linear_eigenpair_relative_residual':linear_res,
        'accepted_periodic_orbit':False,'near_internal_resonance':'3*Omega0 ~ planar sorted4, resolved cubic forcing',
        'formal_3to1_crossing_fundamental_a':crossing,'formal_3to1_crossing_initial_A':crossing_physical,
        'crossing_qualification':'algebraic warning from unaccepted single-mode coefficients, not an internal-resonance boundary',
        'gate_B':'FAIL','no_M1_or_M2_computed':True,
        'numerical_error':'action/weak massform differences and harmonic solve residuals are observed indicators, not rigorous bounds',
        'modal_truncation_error':'all 92 planar p12 DOFs included; no within-p12 planar modal truncation; OUT reduction NOT_RUN',
        'asymptotic_error':'near 1:3 denominator and large third-harmonic correction prevent controlled uniform single-mode approximation over requested amplitudes; no certified bound'})
    p.write(b/'PL_diagnostics.json',diagnostic)
    p.write(b/'diagnostic_completion_provenance.json',{'code_sha256':p.sha(__file__),
        'reason':'finish normalization/residual audit from saved PL coefficients; no repeat of primary coefficient solve',
        'execution_code_preserved':True,'config_changed':False})
    (b/'diagnostic_code').mkdir(exist_ok=True);shutil.copyfile(__file__,b/'diagnostic_code'/Path(__file__).name)
    summary['numerical_seconds']+=time.perf_counter()-start
    summary['calls']['homogeneous_massform_evaluations']=massform.calls+summary.get('prior_massform_evaluations',0)
    summary.update(completed=True,hard_stop=True,overall='PLANAR_SINGLE_MODE_EXPANSION_UNRESOLVED',
        stop_reason='Resolved near 1:3 planar internal resonance violates frozen isolated-PL gate; M1/M2 prohibited',
        candidate_amplitude_selected=False,accepted_nonlinear_frequency=False)
    summary['stages'].update(B='FAIL',C='NOT_RUN',D='PARTIAL_FORMAL_PL_AND_HISTORICAL_M0_ONLY')
    p.finalize(b,item,summary)
    return summary


def plot_saved(b):
    """Formal Stage-B evidence only. Never fills missing M1/M2 values."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=p.read(b/'four_amplitude_diagnostics.json');hist=p.read(b/'M0_historical_extrapolation.json');diag=p.read(b/'PL_diagnostics.json')
    A=np.array([r['A'] for r in rows]);freq=np.array([r['Omega_formal'] for r in rows])
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    ax[0].plot(A,freq,'o--',label='formal single-mode PL; gate FAILED')
    ax[0].axhline(diag['full_planar_omega'][3]/3,color='tab:red',linestyle=':',label='planar mode 4 / 3')
    ax[0].set(xlabel='initial maximum |w| [model length]',ylabel='frequency [1/model time]',title='Unaccepted frequency correction; near 1:3 denominator')
    for key,label in [('mean_energy_ratio','mean'),('second_harmonic_energy_ratio','second harmonic'),('third_harmonic_energy_ratio','third harmonic')]:
        ax[1].plot(A,[r[key] for r in rows],'o-',label=label)
    ax[1].set(xlabel='initial maximum |w| [model length]',ylabel='quadratic stiffness norm ratio to leading term',title='Formal correction sizes (no orbit certification)')
    for a in ax:a.grid(alpha=.2);a.legend(fontsize=7)
    fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(b/('planar_internal_resonance_diagnostic.'+ext),dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,4))
    for name,label in [('velocity_direct','M0 old velocity/direct'),('canonical_direct','M0 canonical/direct')]:
        ax.plot(A,[r['definitions'][name]['margin'] for r in hist],'o--',label=label)
    ax.axhline(0,color='black',linewidth=.8)
    ax.text(.03,.94,'M1/M2: NOT_AVAILABLE (planar gate failed)',fontsize=8,
        transform=ax.transAxes,va='top',bbox={'facecolor':'white','alpha':.8,'edgecolor':'none'})
    ax.set(xlabel='initial harmonic maximum |w| [model length]',ylabel='2g - |effective detuning| [1/model time]',
        title='Historical M0 extrapolation only; not an instability region')
    ax.grid(alpha=.2);ax.legend(fontsize=8);fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(b/('historical_M0_only.'+ext),dpi=180)
    plt.close(fig)
    p.write(b/'figure_provenance.json',{'no_scientific_calls':True,'not_a_confirmed_instability_region':True,
        'M1_M2_not_drawn':'NOT_AVAILABLE after hard gate B','files':['planar_internal_resonance_diagnostic.pdf','historical_M0_only.pdf']})


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',default=str(CONFIG));args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();summary=cache_saved(b);item=p.read(b/'provenance.json')
        if args.plot_only:plot_saved(b);p.finalize(b,item,summary)
    else:
        b,item,summary=prepare(Path(args.config))
        if not summary['completed']:
            if (b/'PL_coefficients.npz').exists():
                cfg=item['config'];model,_,_=p.sources(item['physical_source_config'])
                d=p.make_disc({**p.read(s.c.CONFIG),'lengths':[.9,1.1],'beta_deg':15.},'C',12,model)
                massform=HomogeneousMassForm(d,model)
                with np.load(b/'PL_coefficients.npz') as z:coeff={k:z[k].copy() for k in ('Omega0','Omega1','Omega2','phi','q20','q22','q31','q33')}
                diagnostic=p.read(b/'PL_diagnostics.json');summary['resumed_from_PL_coefficients']=True
                summary['prior_massform_evaluations']=summary['calls']['homogeneous_massform_evaluations']
            else:d,massform,coeff,diagnostic=compute(b,item,summary)
            finish_preflight(b,item,summary,d,massform,coeff,diagnostic);plot_saved(b);p.finalize(b,item,summary)
    print(json.dumps({'bundle':str(b),'summary':summary},ensure_ascii=False))
    return 0
