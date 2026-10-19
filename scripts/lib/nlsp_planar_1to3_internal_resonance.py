"""Fixed beta15, two active planar modes and slaved nonresonant corrections.

Canonical quartic resonant coefficients from the accepted action, followed
by bounded algebraic periodic-candidate checks. The separately authorized
stability mode refines the saved Psi=0 branch and uses spectral Hill algebra.
The separately authorized dynamic mode validates that saved orbit with the
existing physical ODE helpers. No mode runs FEM or changes the physical model.
"""
from __future__ import annotations
import argparse, hashlib, itertools, json, math, shutil, sys, time
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import root
from scripts.lib import nlsp_beta15_nonlinear_planar_resonance as old
from scripts.lib.weakly_nonlinear_spatial_dynamics import _SparseCompiled,_JET_NAMES,_POT_NAMES,_POT_FIELDS,_POT_DERIV

p=old.p;rod=old.rod
CONFIG=p.ROOT/'data/input/nlsp_beta15_planar_1to3_internal_resonance.json'
OUTPUT=p.ROOT/'results/nlsp_beta15_planar_1to3_internal_resonance'
SMOKE=p.ROOT/'results/_smoke/nlsp_beta15_planar_1to3_internal_resonance'
ACTIVE=(1,3)


def add(poly,key,value):
    key=tuple(sorted(key));poly[key]=poly.get(key,0j)+value


def conjugate_key(key):return tuple(sorted(i^1 for i in key))


def reality_defect(poly):
    return max((abs(v-poly.get(conjugate_key(k),0j).conjugate()) for k,v in poly.items()),default=0.)/max(max(map(abs,poly.values()),default=0.),1e-300)


def values(poly,z):return sum(v*np.prod(z[list(k)]) for k,v in poly.items())


def derivatives(poly,z):
    grad=np.zeros(len(z),complex)
    for key,c in poly.items():
        for j in set(key):
            reduced=list(key);reduced.remove(j);grad[j]+=key.count(j)*c*np.prod(z[reduced])
    return grad


class CanonicalAction:
    """Legendre expansion of the saved action, retaining M1 squared in H4.

    H3=V3-T3(q,M0^-1 Pi),
    H4=V4-T4+.5 (T3_v).T M0^-1 T3_v. No physical RHS is changed.
    """
    def __init__(self,d,model,Phi,omega):
        self.d=d;self.Phi=Phi;self.omega=omega;self.calls=0
        T3=model.T4.homogeneous(3);T4=model.T4.homogeneous(4)
        V3=model.V4.homogeneous(3);V4=model.V4.homogeneous(4)
        kinetic=[T3,T4]+[T3.derivative(x) for x in rod.FIELD_ORDER]+[T3.derivative(x+'_t') for x in rod.FIELD_ORDER]
        potential=[V3,V4]+[V3.derivative(x) for x in _POT_NAMES]
        self.compilers=[(_SparseCompiled(kinetic,_JET_NAMES,a.coefficients),_SparseCompiled(potential,_JET_NAMES,a.coefficients)) for a in d.arms]
        self.mass_factor=cho_factor(d.M0,lower=True,check_finite=False)

    def state(self,z):
        Q=(z[0::2]+z[1::2])/np.sqrt(2*self.omega)
        P=-1j*np.sqrt(self.omega/2)*(z[0::2]-z[1::2])
        return self.Phi@Q,self.Phi@P

    def parts(self,q,v,*,component='total'):
        self.calls+=1;energies=np.zeros(4,complex);gV=np.zeros(self.d.ndof,complex);gT=np.zeros_like(gV);mom=np.zeros_like(gV)
        for arm,map0,(kin,pot) in zip(self.d.arms,self.d.maps,self.compilers):
            qi,vi=map0@q,map0@v
            variables=np.vstack([arm.reconstruct(qi).T,arm.reconstruct(qi,derivative=1).T,arm.reconstruct(vi).T,
                np.zeros((21,arm.nq),dtype=np.asarray(q).dtype)])
            kt=kin.evaluate(variables);vp=pot.evaluate(variables)
            energies+=np.array([arm.weights@kt[0],arm.weights@kt[1],arm.weights@vp[0],arm.weights@vp[1]])
            localT=np.zeros(arm.ndof,complex);localP=np.zeros_like(localT);localV=np.zeros_like(localT)
            for i,field in enumerate(rod.FIELD_ORDER):
                localT[arm.slices[field]]+=arm.B[i].T@(arm.weights*kt[2+i])
                localP[arm.slices[field]]+=arm.B[i].T@(arm.weights*kt[9+i])
            for j,(i,deriv) in enumerate(zip(_POT_FIELDS,_POT_DERIV)):
                localV[arm.slices[rod.FIELD_ORDER[i]]]+=(arm.D if deriv else arm.B)[i].T@(arm.weights*vp[2+j])
            gV+=map0.T@localV;gT+=map0.T@localT;mom+=map0.T@localP
        if component=='potential':return energies,gV,np.zeros_like(mom)
        if component=='kinetic':return energies,-gT,mom
        return energies,gV-gT,mom

    def gradient_H3(self,z,component='total'):
        q,v=self.state(z);_,dq,mom=self.parts(q,v,component=component)
        dQ=self.Phi.T@dq;dP=-self.Phi.T@mom
        result=np.zeros(len(z),complex)
        result[0::2]=dQ/np.sqrt(2*self.omega)-1j*np.sqrt(self.omega/2)*dP
        result[1::2]=dQ/np.sqrt(2*self.omega)+1j*np.sqrt(self.omega/2)*dP
        return result

    def H4(self,z):
        q,v=self.state(z);(T3,T4,V3,V4),_,mom=self.parts(q,v)
        inv_mom=cho_solve(self.mass_factor,mom,check_finite=False)
        return V4-T4+.5*(mom@inv_mom)


def cubic_coefficients(action,component='total'):
    """Only H3 terms with zero/one passive variable can affect active H4."""
    n=2*len(action.omega);active=tuple(i for j in ACTIVE for i in (2*j,2*j+1));unit={i:np.eye(n,dtype=complex)[i] for i in active}
    diagonal={i:action.gradient_H3(unit[i],component) for i in active};H3={}
    for i,j in itertools.combinations_with_replacement(active,2):
        g=diagonal[i] if i==j else action.gradient_H3(unit[i]+unit[j],component)-diagonal[i]-diagonal[j]
        for k,c in enumerate(g):
            if k in active:add(H3,(i,j,k),c/3)
            else:add(H3,(i,j,k),c)
    return H3


def quartic_coefficients(action,evaluate=None):
    """Exact fourth-degree polarization, not a fit to dynamics."""
    n=2*len(action.omega);active=tuple(i for j in ACTIVE for i in (2*j,2*j+1));cache={};poly={}
    def eval_counts(counts):
        key=tuple(counts)
        if key not in cache:
            z=np.zeros(n,complex)
            for j,c in zip(active,key):z[j]=c
            cache[key]=action.H4(z) if evaluate is None else evaluate(z)
        return cache[key]
    for key in itertools.combinations_with_replacement(active,4):
        summation=0j
        for mask in range(16):
            counts=[sum(bool(mask&(1<<k)) and j==key[k] for k in range(4)) for j in active]
            summation+=(-1)**(4-mask.bit_count())*eval_counts(counts)
        denominator=math.prod(math.factorial(key.count(j)) for j in set(key));poly[key]=summation/denominator
    return poly


def cubic_elimination(H3,omega):
    """H4_eff=H4+.5{H3,W3}; {H2,W3}=-H3, {a,bar a}=-i."""
    active=set(i for j in ACTIVE for i in (2*j,2*j+1));W3={};gaps=[]
    for key,c in H3.items():
        detuning=sum((1 if i%2==0 else -1)*omega[i//2] for i in key)
        gaps.append(abs(detuning))
        if abs(detuning)<1e-10:raise ArithmeticError('Unresolved cubic homological denominator')
        W3[key]=1j*c/detuning
    pure=[(k,v) for k,v in H3.items() if set(k)<=active];wpure=[(k,W3[k]) for k,_ in pure]
    by_mode=[{} for _ in omega]
    for k,h in pure:
        for m,w in wpure:
            for j in ACTIVE:
                a,abar=2*j,2*j+1;factor=-.5j*(k.count(a)*m.count(abar)-k.count(abar)*m.count(a))
                if factor:
                    key=list(k+m);key.remove(a);key.remove(abar);add(by_mode[j],key,factor*h*w)
    for j in range(len(omega)):
        if j in ACTIVE:continue
        for sign in (0,1):
            a,b=2*j+sign,2*j+1-sign;factor=-.5j if sign==0 else .5j
            first=[(k,h) for k,h in H3.items() if a in k and len(set(k)-active)==1]
            second=[(m,W3[m]) for m,h in H3.items() if b in m and len(set(m)-active)==1]
            for k,h in first:
                for m,w in second:
                    key=list(k+m);key.remove(a);key.remove(b);add(by_mode[j],key,factor*h*w)
    indirect={}
    for poly in by_mode:
        for key,c in poly.items():add(indirect,key,c)
    return W3,indirect,by_mode,min(gaps)


def resonant_constants(poly):
    a,b=2*ACTIVE[0],2*ACTIVE[1]
    keys={'alpha22':(a,a,a+1,a+1),'alpha24':(a,a+1,b,b+1),'alpha44':(b,b,b+1,b+1),'gamma':(a,a,a,b+1)}
    return {name:poly.get(tuple(sorted(key)),0j) for name,key in keys.items()}


def load_problem(item):
    cfg=item['config'];parent=p.ROOT/cfg['planar_parent'];prior=p.read(parent/'provenance.json')
    model,_,_=p.sources(prior['physical_source_config']);pc={**p.read(old.s.c.CONFIG),'lengths':[.9,1.1],'beta_deg':15.}
    d=p.make_disc(pc,'C',12,model)
    with np.load(parent/'PL_coefficients.npz') as z:omega=z['planar_omega'].copy();vectors=z['planar_phi'].copy();previous_phi=z['phi'].copy()
    Phi=np.zeros((d.ndof,len(omega)));Phi[d.plane]=vectors
    if Phi[:,1]@d.M0@previous_phi<0:Phi[:,1]*=-1
    field=np.concatenate([a.reconstruct(P@Phi[:,3])[:,1] for a,P in zip(d.arms,d.maps)])
    if field[np.argmax(abs(field))]<0:Phi[:,3]*=-1
    return d,model,Phi,omega


def validate_scope(cfg):
    expected={'degree':12,'beta_deg':15.,'mu':.1,'lengths':[.9,1.1],'sections':[[.05,.1],[.05,.1]],
        'material':{'E':1.,'rho':1.,'nu':.3,'kappa':5/6},
        'amplitudes':[.005,.006,.008,.010],'active_planar_sorted_modes':[2,4],
        'primary_harmonics':9,'control_harmonics':7,'fixed_point_iterations':8,'coefficient_phase_samples':64,
        'independent_residual_phases':53,'independent_phase_offset':.173,
        'mass_duality_gate':2e-10,'action_relative_gate':2e-10,'canonical_reality_relative_gate':2e-9,
        'full_planar_residual_relative_gate':1e-3,'energy_relative_variation_gate':1e-4,
        'frequency_control_absolute_gate':2e-5,'nonresonant_energy_ratio_gate':.2,
        'numerical_budget_seconds':600.,'new_ODE_integrations':0,'new_FEM_jobs':0,'OUT_reduction':False,
        'automatic_geometry_or_amplitude_extension':False}
    if any(cfg.get(k)!=v for k,v in expected.items()):raise ValueError('Frozen geometry, algebraic scope or gates changed')
    return cfg


def prepare():
    cfg=validate_scope(p.read(CONFIG))
    if cfg['amplitudes']!=[.005,.006,.008,.010] or cfg['active_planar_sorted_modes']!=[2,4] or cfg['new_ODE_integrations'] or cfg['OUT_reduction']:
        raise ValueError('Fixed planar-only scope required')
    parent=p.ROOT/cfg['planar_parent'];old.cache_saved(parent)
    if p.sha(parent/'manifest.json')!=cfg['planar_parent_manifest_sha256']:raise ValueError('Historical parent changed')
    for path in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        item=p.read(path)
        if item['config']==cfg:return path.parent,item,p.cache(path.parent)
    prior=p.read(parent/'provenance.json')
    for file in old.f.s.r.KERNELS:
        if p.sha(p.ROOT/file)!=prior['code_hashes'][file]:raise ValueError('Frozen physical kernel changed')
    hashes={**p.code_hashes(),Path(__file__).relative_to(p.ROOT).as_posix():p.sha(__file__),CONFIG.relative_to(p.ROOT).as_posix():p.sha(CONFIG),
        'scripts/lib/nlsp_beta15_nonlinear_planar_resonance.py':p.sha(old.__file__)}
    item={'schema':'nlsp-planar-1to3-v1','config':cfg,'code_hashes':hashes,'source_artifacts':prior['source_artifacts'],
        'source_manifests':{**prior['source_manifests'],'planar_parent':{'path':cfg['planar_parent'],'sha256':cfg['planar_parent_manifest_sha256']}},
        'initial_state':p.read(SMOKE/'initial_state.json'),'interpreter':sys.executable}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fp
    if b.exists():raise ValueError('Interrupted attempt must be preserved')
    b.mkdir(parents=True);p.write(b/'provenance.json',item);p.write(b/'config.json',cfg);(b/'execution_code').mkdir()
    for file in hashes:shutil.copyfile(p.ROOT/file,b/'execution_code'/Path(file).name)
    summary={'overall':'PREPARED','completed':False,'numerical_seconds':0.,'calls':{'ODE':0,'FEM':0,'OUT':0,'new_eigenproblems':0},'stages':{'source':'PASS','coefficients':'NOT_RUN','periodic':'NOT_RUN','validation':'NOT_RUN'}}
    p.finalize(b,item,summary);return b,item,summary


def compute_coefficients(b,item,summary):
    start=time.perf_counter();cfg=item['config'];d,model,Phi,omega=load_problem(item)
    dual=np.max(abs(Phi.T@d.M0@Phi-np.eye(len(omega))))
    residual=np.linalg.norm(d.K@Phi-(d.M0@Phi)*omega[None,:]**2)/np.linalg.norm(d.K@Phi)
    if dual>cfg['mass_duality_gate'] or residual>cfg['linear_eigenpair_relative_gate']:raise ArithmeticError('Saved p12 modal gate')
    action=CanonicalAction(d,model,Phi,omega);H3=cubic_coefficients(action);H4=quartic_coefficients(action)
    W3,indirect,by_mode,gap=cubic_elimination(H3,omega);effective=H4.copy()
    for key,c in indirect.items():add(effective,key,c)
    defects={name:reality_defect(poly) for name,poly in [('H3',H3),('W3',W3),('H4_direct',H4),('H4_indirect',indirect),('H4_effective',effective)]}
    if max(defects.values())>cfg['canonical_reality_relative_gate']:raise ArithmeticError('Canonical conjugacy failed')
    constants={name:resonant_constants(poly) for name,poly in [('direct',H4),('indirect',indirect),('effective',effective)]}
    encode=lambda poly:[{'monomial':list(k),'real':v.real,'imag':v.imag} for k,v in sorted(poly.items())]
    result={'frequencies':omega[ list(ACTIVE)].tolist(),'delta_13':float(omega[3]-3*omega[1]),'mass_orthogonality':float(dual),'linear_eigenpair_residual':float(residual),
        'canonical_conjugacy_defects':defects,'minimum_cubic_elimination_frequency_denominator':gap,
        'constants':{label:{k:old.s.complex_row(v) for k,v in coeff.items()} for label,coeff in constants.items()},
        'canonical_convention':'a_j=(sqrt(omega_j) Q_j+i P_j/sqrt(omega_j))/sqrt(2); {a_j,bar a_k}=-i delta_jk',
        'two_active_modes_only':True,'passive_planar_modes_included':90,'H3':encode(H3),'W3':encode(W3),'H4_direct':encode(H4),'H4_indirect':encode(indirect),'H4_effective':encode(effective)}
    p.write(b/'canonical_resonant_coefficients.json',result)
    keys=list(H4);np.savez_compressed(b/'canonical_resonant_coefficients.npz',omega=omega,Phi=Phi,monomials=np.array(keys),
        direct=np.array([H4[k] for k in keys]),indirect=np.array([indirect.get(k,0j) for k in keys]),
        indirect_by_mode=np.array([[x.get(k,0j) for k in keys] for x in by_mode]))
    summary['numerical_seconds']+=time.perf_counter()-start;summary['calls']['canonical_action_evaluations']=action.calls
    summary['stages']['coefficients']='PASS';summary['overall']='COEFFICIENTS_SAVED';p.finalize(b,item,summary)
    return result


def stationary_normal_form(constants,omega,W2,W4,A,seeds):
    """Stationarity of H-Ω(I2+3I4); both signs of the real mode4 amplitude."""
    aa,ab,bb,g=(constants[k] for k in ('alpha22','alpha24','alpha44','gamma'))
    w2,w4=omega[list(ACTIVE)];scale=A/W2*np.sqrt(w2/2);candidates=[]
    def fun(u):
        r,t=u;rr=r*scale;ss=t*scale
        om=w2+2*aa*rr*rr+ab*ss*ss+3*g*rr*ss
        second=(w4-3*om)*ss+ab*rr*rr*ss+2*bb*ss**3+g*rr**3
        x=np.sqrt(2/w2)*rr;y=np.sqrt(2/w4)*ss
        # A leading-only scalar norm; full physical normalization is imposed
        # later on the reconstructed HB state, not in the normal form.
        physical=np.sqrt((W2*x)**2+(W4*y)**2)
        return [(physical-A)/A,second/(w4*scale)]
    for seed in seeds:
        result=root(fun,[1.,seed],tol=1e-11)
        if not result.success or np.linalg.norm(fun(result.x))>1e-8:continue
        r,t=result.x;rr=r*scale;ss=t*scale
        if rr<=scale*1e-5:continue
        om=w2+2*aa*rr*rr+ab*ss*ss+3*g*rr*ss
        x=np.sqrt(2/w2)*rr;y=np.sqrt(2/w4)*ss
        if any(np.linalg.norm(np.array([x,y])-np.array([z['x'],z['y']]))<1e-8*A/W2 for z in candidates):continue
        candidates.append({'x':x,'y':y,'Omega':om,'normal_form_stationary_residual':float(np.linalg.norm(fun(result.x))),
            'phase_Psi':0. if ss>=0 else float(np.pi),'normal_form_actions':[float(rr*rr),float(ss*ss)],
            'amplitude_qualification':'root seed only; preliminary root-sum-square leading profile scale, not final physical A'})
    return sorted(candidates,key=lambda z:z['y'])


class PeriodicBalance:
    """Two active slots; all other p12 Fourier slots are algebraically slaved.

    This is not a global many-active-mode normal form or a time integration.
    The two resonant slots are never divided by a linear denominator.
    """
    def __init__(self,d,model,Phi,omega,cfg,deadline):
        self.d=d;self.Phi=Phi;self.omega=omega;self.cfg=cfg;self.deadline=deadline
        self.massform=old.HomogeneousMassForm(d,model);self.evaluations=0
        self.phases=np.arange(cfg['coefficient_phase_samples'])*2*np.pi/cfg['coefficient_phase_samples']
        self.cos=np.cos(np.arange(cfg['primary_harmonics']+1)[:,None]*self.phases)
        self.sin=np.sin(np.arange(cfg['primary_harmonics']+1)[:,None]*self.phases)

    def fields(self,U,Omega,phases=None):
        if phases is None:C=self.cos[:len(U)];S=self.sin[:len(U)]
        else:C=np.cos(np.arange(len(U))[:,None]*phases);S=np.sin(np.arange(len(U))[:,None]*phases)
        n=np.arange(len(U))[:,None]
        Q=C.T@U;V=-Omega*S.T@(n*U);acc=-Omega**2*C.T@(n*n*U)
        return Q@self.Phi.T,V@self.Phi.T,acc@self.Phi.T

    def forcing(self,U,Omega):
        if time.perf_counter()>self.deadline:raise TimeoutError('600 s numerical budget')
        self.evaluations+=1;q,v,a=self.fields(U,Omega);values=[]
        for qi,vi,ai in zip(q,v,a):
            values.append(self.Phi.T@(self.massform.residual(2,qi,vi,ai)+self.massform.residual(3,qi,vi,ai)))
        return old.cosine_coefficients(np.array(values),self.phases,len(U)-1)

    def evaluate(self,x,y,Omega,harmonics,iterations=None):
        U=np.zeros((harmonics+1,len(self.omega)));U[1,1]=x;U[3,3]=y
        n=np.arange(harmonics+1)[:,None];den=self.omega[None,:]**2-(n*Omega)**2
        mask=np.ones_like(U,dtype=bool);mask[1,1]=False;mask[3,3]=False
        history=[]
        for _ in range(self.cfg['fixed_point_iterations'] if iterations is None else iterations):
            F=self.forcing(U,Omega);new=U.copy();new[mask]=-F[mask]/den[mask]
            history.append(float(np.linalg.norm((new-U)*self.omega[None,:])))
            U=new
        F=self.forcing(U,Omega);equations=den*U+F
        amplitude=old.p.polynomial_w_max(self.d,self.Phi@U.sum(axis=0))[0]
        return U,equations,amplitude,history

    def solve(self,A,seed,harmonics):
        W2=p.polynomial_w_max(self.d,self.Phi[:,1])[0];W4=p.polynomial_w_max(self.d,self.Phi[:,3])[0]
        def fun(z):
            x=z[0]*A/W2;y=z[1]*A/W4;Omega=z[2]
            U,eq,amp,_=self.evaluate(x,y,Omega,harmonics)
            return [eq[1,1]/(self.omega[1]**2*A/W2),eq[3,3]/(self.omega[3]**2*A/W4),(amp-A)/A]
        initial=[seed['x']*W2/A,seed['y']*W4/A,seed['Omega']]
        result=root(fun,initial,tol=1e-10,options={'maxfev':self.cfg['max_root_evaluations']})
        x=result.x[0]*A/W2;y=result.x[1]*A/W4;Omega=result.x[2]
        U,eq,amp,history=self.evaluate(x,y,Omega,harmonics)
        diagnostics={'success':bool(result.success),'message':str(result.message),'nfev':int(result.nfev),
            'scaled_root_residual':float(np.linalg.norm(fun(result.x))),'fixed_point_history':history,
            'A':A,'actual_initial_max_w':amp,'x':x,'y':y,'Omega':Omega,'phase_Psi':0. if y>=0 else float(np.pi),
            'harmonics':harmonics,'resonant_slots_not_inverted':[[1,2],[3,4]]}
        return U,diagnostics

    def validate(self,U,diagnostics):
        count=self.cfg['independent_residual_phases'];phases=np.arange(count)*2*np.pi/count+self.cfg['independent_phase_offset']
        Omega=diagnostics['Omega'];q,v,a=self.fields(U,Omega,phases);R=[];energy=[]
        for qi,vi,ai in zip(q,v,a):
            R.append(sum(self.massform.residual(n,qi,vi,ai) for n in (1,2,3)))
            energy.append(self.d.energy(qi,vi))
        residual=np.array(R);scale=Omega**2*np.sqrt(diagnostics['x']**2+81*diagnostics['y']**2)
        weights=U*self.omega[None,:];leading=np.zeros_like(U);leading[1,1]=U[1,1];leading[3,3]=U[3,3]
        corrections=(U-leading)*self.omega[None,:]
        nr=float(np.linalg.norm(corrections)/max(np.linalg.norm(leading*self.omega[None,:]),1e-300))
        n=np.arange(len(U))[:,None];den=abs(self.omega[None,:]**2-(n*Omega)**2);den[1,1]=np.inf;den[3,3]=np.inf
        loc=np.unravel_index(np.argmin(den),den.shape)
        eigE=.5*self.omega[list(ACTIVE)]**2*np.array([diagnostics['x']**2,diagnostics['y']**2])
        diagnostic={**diagnostics,'independent_phase_count':count,'phase_offset':self.cfg['independent_phase_offset'],
            'full_massform_residual_relative':float(np.max(np.linalg.norm(residual,axis=1))/max(scale,1e-300)),
            'energy_relative_variation':float((max(energy)-min(energy))/abs(np.mean(energy))),
            'energy_mean':float(np.mean(energy)),'nonresonant_correction_energy_norm_ratio':nr,
            'linear_active_energy_parts':eigE.tolist(),'linear_mode4_energy_fraction':float(eigE[1]/sum(eigE)),
            'energy_partition_qualification':'leading quadratic reference energies, not separately conserved nonlinear modal energies',
            'minimum_remaining_squared_frequency_denominator':float(den[loc]),'remaining_denominator_harmonic':int(loc[0]),'remaining_denominator_sorted_mode':int(loc[1]+1),
            'mean_stiffness_norm':float(np.linalg.norm(weights[0])),
            'second_harmonic_stiffness_norm':float(np.linalg.norm(weights[2])),
            'third_harmonic_stiffness_norm':float(np.linalg.norm(weights[3])),
            'reconstructed_initial_velocity_max':0.,'cosine_phase_closure_exact':True}
        return diagnostic,{'phases':phases,'q':q,'v':v,'acceleration':a,'massform_residual':residual,'energy':np.array(energy),'U':U}


def periodic_candidates(b,item,summary):
    start=time.perf_counter();cfg=item['config'];deadline=start+cfg['numerical_budget_seconds']-summary['numerical_seconds']
    d,model,Phi,omega=load_problem(item);coeff=p.read(b/'canonical_resonant_coefficients.json')
    constants={k:v['real'] for k,v in coeff['constants']['effective'].items()}
    W2=p.polynomial_w_max(d,Phi[:,1])[0];W4=p.polynomial_w_max(d,Phi[:,3])[0]
    balance=PeriodicBalance(d,model,Phi,omega,cfg,deadline);rows=[];allseeds={}
    try:
        for A in cfg['amplitudes']:
            seeds=stationary_normal_form(constants,omega,W2,W4,A,cfg['algebraic_ratio_seeds']);allseeds[str(A)]=seeds
            for j,seed in enumerate(seeds):
                label=f'A{A:.3f}_branch{j+1}';U,rootinfo=balance.solve(A,seed,cfg['primary_harmonics'])
                diagnostic,raw=balance.validate(U,rootinfo)
                p.write(b/(label+'_diagnostics.json'),diagnostic);np.savez_compressed(b/(label+'_raw.npz'),**raw)
                rows.append(diagnostic);p.write(b/'periodic_candidates.json',rows);p.write(b/'stationary_normal_form_seeds.json',allseeds)
                summary['numerical_seconds']+=time.perf_counter()-start;start=time.perf_counter();p.finalize(b,item,summary)
    except TimeoutError:
        summary.update(overall='NUMERICAL_PARTIAL',hard_stop=True,stop_reason='600 s algebraic budget; completed candidates preserved')
    summary['numerical_seconds']+=time.perf_counter()-start;summary['calls'].update(HB_forcing_evaluations=balance.evaluations,homogeneous_massform_evaluations=balance.massform.calls)
    summary['stages']['periodic']='CANDIDATES_SAVED';p.write(b/'periodic_candidates.json',rows);p.write(b/'stationary_normal_form_seeds.json',allseeds);p.finalize(b,item,summary)
    return rows


def normal_form_hamiltonian(z,omega,constants):
    r,t=z[2*ACTIVE[0]],z[2*ACTIVE[1]];I2=abs(r)**2;I4=abs(t)**2
    aa,ab,bb,g=(constants[k] for k in ('alpha22','alpha24','alpha44','gamma'))
    return float(omega[1]*I2+omega[3]*I4+aa*I2**2+ab*I2*I4+bb*I4**2+2*np.real(g*r**3*t.conjugate()))


def verify_and_control(b,item,summary):
    start=time.perf_counter();cfg=item['config'];deadline=start+cfg['numerical_budget_seconds']-summary['numerical_seconds']
    d,model,Phi,omega=load_problem(item);data=p.read(b/'canonical_resonant_coefficients.json')
    constants={k:v['real'] for k,v in data['constants']['effective'].items()}
    H3={tuple(row['monomial']):complex(row['real'],row['imag']) for row in data['H3']}
    H4={tuple(row['monomial']):complex(row['real'],row['imag']) for row in data['H4_direct']}
    W3={tuple(row['monomial']):complex(row['real'],row['imag']) for row in data['W3']}
    action=CanonicalAction(d,model,Phi,omega);balance=PeriodicBalance(d,model,Phi,omega,cfg,deadline)
    rows=p.read(b/'periodic_candidates.json');controlled=[];rawchecks=[]
    try:
        for j,row in enumerate(rows):
            label=f'A{row["A"]:.3f}_branch{sum(r["A"]==row["A"] for r in rows[:j])+1}'
            control,info=balance.solve(row['A'],row,cfg['control_harmonics']);info,raw=balance.validate(control,info)
            diff=abs(info['Omega']-row['Omega'])
            final={**row,'harmonic_control_Omega':info['Omega'],'harmonic_control_frequency_difference':diff,
                'harmonic_control_relative_mode4_change':abs(info['y']-row['y'])/max(abs(row['y']),1e-300),
                'harmonic_control_residual_relative':info['full_massform_residual_relative'],
                'criterion_passed':bool(row['success'] and info['success'] and row['full_massform_residual_relative']<=cfg['full_planar_residual_relative_gate']
                    and row['energy_relative_variation']<=cfg['energy_relative_variation_gate'] and diff<=cfg['frequency_control_absolute_gate']
                    and row['nonresonant_correction_energy_norm_ratio']<=cfg['nonresonant_energy_ratio_gate'])}
            with np.load(b/(label+'_raw.npz')) as z:U=z['U'].copy();q=z['q'].copy();v=z['v'].copy();acc=z['acceleration'].copy()
            agreement=[];joint=[];mass=[]
            for k in [0,17,35]:
                r=sum(balance.massform.residual(n,q[k],v[k],acc[k]) for n in (1,2,3))
                direct=d.mass_matrix(q[k])@acc[k]+d.inertial_terms(q[k],v[k])+d.potential(q[k])['gradient']
                agreement.append(float(np.linalg.norm(r-direct)/max(np.linalg.norm(direct),1e-300)))
                joint.append(d.joint_diagnostics(q[k],v[k]))
                # Bound with the already accepted pointwise mass matrix, no
                # frequency scan: resting mass is whitened and positive.
                mass.append(float(np.linalg.eigvalsh(d.mass_matrix(q[k])).min()))
            final['massform_direct_equation_relative_differences']=agreement
            final['joint_compatibility_max']=max(max(x['translation_difference'],x['orientation_difference'],x['contraction_difference']) for x in joint)
            final['pointwise_mass_min_eigenvalue']=min(mass)
            # Near-zero residual needs an absolute comparison as well.
            final['criterion_passed'] &= final['joint_compatibility_max']<cfg['mass_duality_gate'] and min(mass)>0
            p.write(b/(label+'_harmonic_control.json'),info);np.savez_compressed(b/(label+'_control_raw.npz'),**raw)
            controlled.append(final);p.write(b/'validated_periodic_candidates.json',controlled)
            summary['numerical_seconds']+=time.perf_counter()-start;start=time.perf_counter();p.finalize(b,item,summary)
        W2=p.polynomial_w_max(d,Phi[:,1])[0]
        prior=p.read(p.ROOT/cfg['planar_parent']/'PL_diagnostics.json')
        shift=constants['alpha22']*omega[1]/W2**2
        # Check the cubic polynomial and its full passive gradient on physical
        # active states; omitted H3 terms contain >=2 passive variables.
        for phase in [.2,.9,2.3]:
            z=np.zeros(2*len(omega),complex);z[2]=.00013*np.exp(-1j*phase);z[3]=z[2].conjugate()
            z[6]=.000025*np.exp(-3j*phase);z[7]=z[6].conjugate()
            q,v=action.state(z);parts,_,_=action.parts(q,v)
            h3=parts[2]-parts[0];h4=action.H4(z)
            rawchecks.append({'phase':phase,'H3_absolute_polynomial_difference':float(abs(values(H3,z)-h3)),
                'H4_absolute_polynomial_difference':float(abs(values(H4,z)-h4)),
                'H3_full_gradient_relative_difference':float(np.linalg.norm(derivatives(H3,z)-action.gradient_H3(z))/np.linalg.norm(action.gradient_H3(z)))})
        typeI=[];typeI_states={}
        for A in cfg['amplitudes']:
            z=np.zeros(2*len(omega),complex);z[2]=z[3]=np.sqrt(omega[1]/2)*A/W2
            gW=derivatives(W3,z);correction=np.empty_like(z);correction[0::2]=-1j*gW[1::2];correction[1::2]=1j*gW[0::2]
            mapped=z-correction;I2=abs(mapped[2])**2;I4=abs(mapped[6])**2
            typeI_states[f'A{A:.3f}']=mapped.copy()
            passive=[j for j in range(len(omega)) if j not in ACTIVE]
            canonical_passive_energy=sum(omega[j]*abs(mapped[2*j])**2 for j in passive)
            source=abs(-1j*constants['gamma']*mapped[2]**3)
            typeI.append({'A':A,'physical_initial_Q4':0.,'canonical_I2_through_quadratic_map':I2,
                'canonical_I4_through_quadratic_map':I4,'canonical_J':I2+3*I4,
                'canonical_mode4_initial_real':mapped[6].real,'canonical_mode4_initial_imag':mapped[6].imag,
                'normal_form_H':normal_form_hamiltonian(mapped,omega,constants),
                'pure_mode2_resonant_mode4_source_magnitude':source,
                'passive_initial_reference_energy_from_inverse_map':canonical_passive_energy,
                'qualification':'inverse cubic-generator map to quadratic coordinate order; free passive transients retained as a limitation; no type-I trajectory or modulation integrated'})
        verification={'self_shift_per_physical_amplitude_squared':shift,'previous_PL_Omega2':prior['Omega2'],
            'self_shift_difference':shift-prior['Omega2'],'canonical_polynomial_projection_checks':rawchecks,
            'J_invariant':'I2+3 I4, coefficient-exact commutation of resonant monomials with J',
            'phase_dependency':'2 gamma I2^(3/2) sqrt(I4) cos(3 theta2-theta4)',
            'actual_growth_or_OUT_analysis':'NOT_RUN'}
        p.write(b/'canonical_variational_checks.json',verification);p.write(b/'type_I_initial_modulation_diagnostics.json',typeI)
        np.savez_compressed(b/'type_I_canonical_initial_states.npz',**typeI_states)
        with np.load(b/'canonical_resonant_coefficients.npz') as z:
            monomials=[tuple(k) for k in z['monomials']];direct=z['direct'];by_mode=z['indirect_by_mode'];key=monomials.index((2,2,2,7))
            modal=[{'retained_sorted_modes':k,'gamma_direct_plus_cubic_elimination':float((direct[key]+by_mode[:k,key].sum()).real)} for k in [8,16,32,64,92]]
        p.write(b/'passive_mode_coefficient_convergence.json',modal)
    except TimeoutError:
        summary.update(overall='NUMERICAL_PARTIAL',completed=True,hard_stop=True,stop_reason='600 s budget during controls')
    summary['numerical_seconds']+=time.perf_counter()-start
    summary['calls']['canonical_action_evaluations']+=action.calls
    summary['calls']['HB_forcing_evaluations']+=balance.evaluations
    summary['calls']['homogeneous_massform_evaluations']+=balance.massform.calls
    if not summary.get('hard_stop'):
        main_rows=[r for r in controlled if r['phase_Psi']==0.]
        all_A=sorted(r['A'] for r in main_rows)==cfg['amplitudes']
        passed=all_A and all(r['criterion_passed'] for r in main_rows)
        summary.update(completed=True,overall='PLANAR_1TO3_RESONANT_REDUCTION_VALIDATED' if passed else 'PLANAR_1TO3_APPROXIMATION_UNRESOLVED',
            periodic_branch_identified=True,type_I_modulation_indicated=True,
            suitable_for_future_OUT_formulation=passed,qualification='validated finite-harmonic p12 type-II approximation; not orbit stability, not a prediction that type-I initial conditions lie on this branch')
        summary['stages'].update(periodic='PASS',validation='PASS' if passed else 'FAIL')
    p.write(b/'validated_periodic_candidates.json',controlled);p.finalize(b,item,summary)
    return controlled


def complete_diagnostics(b,item,summary):
    """Small coefficient decomposition and cached-candidate physical observables."""
    start=time.perf_counter();cfg=item['config'];d,model,Phi,w=load_problem(item)
    data=p.read(b/'canonical_resonant_coefficients.json');action=CanonicalAction(d,model,Phi,w)
    decode=lambda key:{tuple(r['monomial']):complex(r['real'],r['imag']) for r in data[key]}
    H3V=cubic_coefficients(action,'potential');H3T=cubic_coefficients(action,'kinetic')
    _,indV,_,_=cubic_elimination(H3V,w);_,indT,_,_=cubic_elimination(H3T,w)
    indAll=decode('H4_indirect');indCross={k:indAll.get(k,0)-indV.get(k,0)-indT.get(k,0) for k in indAll}
    V4=quartic_coefficients(action,lambda z:action.parts(*action.state(z))[0][3]);H4=decode('H4_direct')
    T4={k:H4[k]-V4[k] for k in H4}
    parts={'direct_V4':V4,'direct_rotational_Legendre_terms':T4,'indirect_V3_V3':indV,'indirect_T3_T3':indT,'indirect_cross_V3_T3':indCross}
    p.write(b/'coefficient_energy_decomposition.json',{name:{k:old.s.complex_row(v) for k,v in resonant_constants(poly).items()} for name,poly in parts.items()})
    constants={k:v['real'] for k,v in data['constants']['effective'].items()};W3=decode('W3')
    rows=p.read(b/'validated_periodic_candidates.json');checks=[]
    for j,r in enumerate(rows):
        label=f'A{r["A"]:.3f}_branch{sum(t["A"]==r["A"] for t in rows[:j])+1}'
        with np.load(b/(label+'_raw.npz')) as z:U=z['U'].copy();q=z['q'].copy();v=z['v'].copy();acc=z['acceleration'].copy();R=z['massform_residual'].copy()
        diffs=[]
        for k in [0,17,35]:
            inertia=d.inertial_terms(q[k],v[k]);grad=d.potential(q[k])['gradient'];ma=d.mass_matrix(q[k])@acc[k]
            diffs.append(float(np.linalg.norm(R[k]-ma-inertia-grad)/(np.linalg.norm(ma)+np.linalg.norm(inertia)+np.linalg.norm(grad))))
        if max(diffs)>cfg['action_relative_gate']:raise ArithmeticError('Action/full planar operator disagreement')
        r['action_massform_difference_relative_to_equation_terms']=diffs
        r['mode2_leading_max_w']=p.polynomial_w_max(d,Phi[:,1]*r['x'])[0];r['mode4_leading_max_w']=p.polynomial_w_max(d,Phi[:,3]*r['y'])[0]
        r['mean_max_w']=p.polynomial_w_max(d,Phi@U[0])[0];r['second_harmonic_max_w']=p.polynomial_w_max(d,Phi@U[2])[0]
        third=U[3].copy();third[3]=0.;r['third_harmonic_nonresonant_max_w']=p.polynomial_w_max(d,Phi@third)[0]
        r['third_harmonic_nonresonant_stiffness_norm']=float(np.linalg.norm(third*w))
        r['mode4_to_mode2_stiffness_norm_ratio']=abs(w[3]*r['y']/(w[1]*r['x']))
        r['linear_modal_detuning_remaining_from_common_frequency']=float(w[3]-3*r['Omega'])
        prior=next(t for t in p.read(p.ROOT/cfg['planar_parent']/'four_amplitude_diagnostics.json') if t['A']==r['A'])
        r['comparison_to_prior_formal_single_mode_Omega']=r['Omega']-prior['Omega_formal']
        # Distinguish an NF diagnostic of mapped HB initial data from an exact
        # canonical amplitude parameterization at all asymptotic orders.
        leadQ=np.zeros(len(w));leadQ[1]=r['x'];leadQ[3]=r['y'];zlead=np.repeat(np.sqrt(w/2)*leadQ,2).astype(complex)
        g=derivatives(W3,zlead);correction=np.empty_like(g);correction[0::2]=-1j*g[1::2];correction[1::2]=1j*g[0::2]
        physicalQ=U.sum(axis=0);mapped=np.repeat(np.sqrt(w/2)*physicalQ,2)-correction
        I2,I4=abs(mapped[2])**2,abs(mapped[6])**2
        r['NF_diagonal_detuning_at_quadratically_mapped_HB_initial_state']=float(w[3]-3*w[1]+(constants['alpha24']-6*constants['alpha22'])*I2+(2*constants['alpha44']-3*constants['alpha24'])*I4)
        r['NF_diagnostic_qualification']='quadratic inverse coordinate map only; not a complete finite-amplitude canonical reconstruction'
        checks.append({'label':label,'operator_scaled_agreement':diffs})
    profiles=[]
    for mode in ACTIVE:
        profiles.append([arm.reconstruct(P@Phi[:,mode],np.linspace(0,arm.length,121)) for arm,P in zip(d.arms,d.maps)])
    np.savez_compressed(b/'linear_mode_physical_profiles.npz',profiles=np.array(profiles),lengths=d.lengths,field_order=rod.FIELD_ORDER,omega=w[list(ACTIVE)])
    p.write(b/'validated_periodic_candidates.json',rows);p.write(b/'operator_scaled_action_checks.json',checks)
    summary['numerical_seconds']+=time.perf_counter()-start;summary['calls']['canonical_action_evaluations']+=action.calls
    p.finalize(b,item,summary)


def frequency_sensitivity_check(b,item,summary):
    """Observed stationary-root conditioning and one passive-update error.

    No new branch solve, extra amplitude or time trajectory. The prescribed
    eight-iteration production solution remains unchanged.
    """
    start=time.perf_counter();cfg=item['config'];d,model,Phi,w=load_problem(item)
    balance=PeriodicBalance(d,model,Phi,w,cfg,start+cfg['numerical_budget_seconds']-summary['numerical_seconds'])
    rows=p.read(b/'validated_periodic_candidates.json');records=[]
    W2=p.polynomial_w_max(d,Phi[:,1])[0];W4=p.polynomial_w_max(d,Phi[:,3])[0]
    for j,row in enumerate(rows):
        A=row['A'];label=f'A{A:.3f}_branch{sum(r["A"]==A for r in rows[:j])+1}'
        with np.load(b/(label+'_raw.npz')) as z:U=z['U'].copy()
        center=np.array([row['x']*W2/A,row['y']*W4/A,(row['Omega']-w[1])/w[1]])
        def fun(z):
            _,eq,amp,_=balance.evaluate(z[0]*A/W2,z[1]*A/W4,w[1]*(1+z[2]),cfg['primary_harmonics'])
            return np.array([eq[1,1]/(w[1]**2*A/W2),eq[3,3]/(w[3]**2*A/W4),(amp-A)/A])
        J=np.zeros((3,3));step=1e-6
        for k in range(3):
            delta=np.eye(3)[k]*step;J[:,k]=(fun(center+delta)-fun(center-delta))/(2*step)
        forcing=balance.forcing(U,row['Omega']);n=np.arange(len(U))[:,None];den=w[None,:]**2-(n*row['Omega'])**2
        mask=np.ones_like(U,dtype=bool);mask[1,1]=mask[3,3]=False;updated=U.copy();updated[mask]=-forcing[mask]/den[mask]
        equations=den*updated+balance.forcing(updated,row['Omega'])
        amp=p.polynomial_w_max(d,Phi@updated.sum(axis=0))[0]
        defect=np.array([equations[1,1]/(w[1]**2*A/W2),equations[3,3]/(w[3]**2*A/W4),(amp-A)/A])
        correction=np.linalg.solve(J,-defect);frequency_indicator=abs(w[1]*correction[2])
        row['scaled_stationary_root_Jacobian_condition']=float(np.linalg.cond(J))
        row['one_passive_update_frequency_correction_indicator']=float(frequency_indicator)
        row['frequency_uncertainty_qualification']='observed harmonic difference and linearized one-update indicator; not a certified error bound'
        row['criterion_passed'] &= frequency_indicator<=cfg['frequency_control_absolute_gate']
        records.append({'label':label,'scaled_Jacobian':J.tolist(),'normalized_defect_after_one_passive_update':defect.tolist(),
            'linearized_correction':correction.tolist(),'frequency_indicator':frequency_indicator})
    p.write(b/'stationary_frequency_sensitivity.json',records);p.write(b/'validated_periodic_candidates.json',rows)
    main=[r for r in rows if r['phase_Psi']==0.]
    if not all(r['criterion_passed'] for r in main):summary.update(overall='PLANAR_1TO3_APPROXIMATION_UNRESOLVED',suitable_for_future_OUT_formulation=False)
    summary['numerical_seconds']+=time.perf_counter()-start
    summary['calls']['HB_forcing_evaluations']+=balance.evaluations;summary['calls']['homogeneous_massform_evaluations']+=balance.massform.calls
    p.finalize(b,item,summary)


def plot_saved(b):
    if not (b/'validated_periodic_candidates.json').exists():
        p.write(b/'figure_availability.json',{'status':'NOT_AVAILABLE','reason':'No controlled periodic candidates; no misleading figure created'})
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=p.read(b/'validated_periodic_candidates.json');main=[r for r in rows if r['phase_Psi']==0.]
    A=[r['A'] for r in main];freq=[r['Omega'] for r in main];meta=p.read(b/'canonical_resonant_coefficients.json')
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    ax[0].plot(A,freq,'o-',label='Psi=0 HB candidate, type II')
    other=[r for r in rows if r['phase_Psi']>0]
    ax[0].plot([r['A'] for r in other],[r['Omega'] for r in other],'s',label='Psi=pi candidates; no stability claim')
    ax[0].axhline(meta['frequencies'][1]/3,color='gray',linestyle=':',label='linear omega4 / 3')
    ax[0].set(ylabel='common frequency [1/model time]',title='Fixed p12, algebraic periodic candidates')
    ax[1].plot(A,[r['linear_modal_detuning_remaining_from_common_frequency'] for r in main],'o-',label='linear omega4 - 3*Omega_HB')
    ax[1].axhline(0,color='gray',linewidth=.8);ax[1].set(ylabel='frequency difference [1/model time]',title='Linear reference detuning (not full NF detuning)')
    for a in ax:a.set_xlabel('initial max |w| [model length]');a.grid(alpha=.2);a.legend(fontsize=7)
    fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(b/('two_mode_frequencies.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
    plt.close(fig)
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    ax[0].plot(A,[r['mode2_leading_max_w'] for r in main],'o-',label='mode 2 leading |w| maximum')
    ax[0].plot(A,[r['mode4_leading_max_w'] for r in main],'o-',label='mode 4 leading |w| maximum')
    ax[0].plot(A,A,'k:',label='total initial |w| maximum A')
    ax[0].set(ylabel='maximum displacement [model length]',title='Initial type-II components, Psi=0')
    ax[1].plot(A,[100*r['linear_mode4_energy_fraction'] for r in main],'o-',label='mode 4 quadratic active-energy share')
    ax[1].set(ylabel='share [%]',ylim=(0,60),title='Initial quadratic reference energy share')
    for a in ax:a.set_xlabel('initial max |w| [model length]');a.grid(alpha=.2);a.legend(fontsize=7)
    fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(b/('two_mode_participation.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
    plt.close(fig)
    p.write(b/'figure_provenance.json',{'figures':['two_mode_frequencies.pdf','two_mode_participation.pdf'],
        'no_scientific_calls':True,'not_stability_diagram':True,'no_independent_curve_normalization':True})


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();summary=p.cache(b)
        if args.plot_only:plot_saved(b);p.finalize(b,p.read(b/'provenance.json'),summary)
    else:
        b,item,summary=prepare();was_completed=summary['completed']
        started=time.perf_counter();initial_seconds=summary['numerical_seconds']
        try:
            if not summary['completed'] and not summary.get('hard_stop') and not (b/'canonical_resonant_coefficients.json').exists():compute_coefficients(b,item,summary)
            if not summary['completed'] and not summary.get('hard_stop') and not (b/'periodic_candidates.json').exists():periodic_candidates(b,item,summary)
            if not summary['completed'] and not summary.get('hard_stop'):verify_and_control(b,item,summary)
            if not was_completed and summary['completed'] and not summary.get('hard_stop'):
                complete_diagnostics(b,item,summary);frequency_sensitivity_check(b,item,summary);plot_saved(b);p.finalize(b,item,summary)
        except (TimeoutError,ArithmeticError,np.linalg.LinAlgError) as error:
            summary.update(completed=True,hard_stop=True,overall='NUMERICAL_PARTIAL' if isinstance(error,TimeoutError) else 'BLOCKED',
                stop_reason=type(error).__name__+': '+str(error))
            summary['numerical_seconds']=max(summary['numerical_seconds'],initial_seconds+time.perf_counter()-started)
            p.finalize(b,item,summary)
    print(json.dumps(p.plain({'bundle':str(b),'summary':summary}),ensure_ascii=False))
    return 0


# Full-space HB and Hill diagnostics of the SAME fixed Type-II branch.
# Separate authorization/cache prevents this mode from changing old results.
STABILITY_CONFIG=p.ROOT/'data/input/nlsp_beta15_periodic_orbit_spatial_stability.json'
STABILITY_OUTPUT=p.ROOT/'results/nlsp_beta15_periodic_orbit_spatial_stability'
STABILITY_SMOKE=p.ROOT/'results/_smoke/nlsp_beta15_periodic_orbit_spatial_stability'


def stability_scope(cfg):
    if cfg!=p.read(STABILITY_CONFIG):raise ValueError('Frozen stability configuration required')
    if (cfg['degree']!=12 or cfg['beta_deg']!=15 or cfg['mu']!=.1 or cfg['phase_Psi']!=0
        or cfg['amplitudes']!=[.005,.006,.008,.010] or cfg['new_ODE_integrations'] or cfg['new_FEM_jobs']
        or cfg['automatic_extension'] or cfg['numerical_budget_seconds']>900):raise ValueError('Stability scope expansion')
    return cfg


def stability_prepare():
    cfg=stability_scope(p.read(STABILITY_CONFIG));parent=p.ROOT/cfg['parent']
    if p.sha(parent/'manifest.json')!=cfg['parent_manifest_sha256']:raise ValueError('Historical parent manifest changed')
    p.cache(parent);prior=p.read(parent/'provenance.json')
    for name in old.f.s.r.KERNELS:
        if p.sha(p.ROOT/name)!=prior['code_hashes'][name]:raise ValueError('Physical kernel changed: '+name)
    for existing in STABILITY_OUTPUT.glob('*/provenance.json') if STABILITY_OUTPUT.exists() else []:
        item=p.read(existing)
        if item['config']!=cfg:raise ValueError('Existing attempt has different configuration')
        return existing.parent,item,p.cache(existing.parent)
    files=[STABILITY_CONFIG.relative_to(p.ROOT).as_posix(),Path(__file__).resolve().relative_to(p.ROOT).as_posix(),
        'scripts/lib/nlsp_instability_feasibility.py','scripts/lib/nlsp_beta15_nonlinear_planar_resonance.py']
    item={'schema':'nlsp-beta15-HB-Hill-v1','config':cfg,'code_hashes':{**p.code_hashes(),**{f:p.sha(p.ROOT/f) for f in files}},
        'source_artifacts':{**prior['source_artifacts'],**{f.relative_to(p.ROOT).as_posix():p.sha(f) for f in parent.iterdir() if f.is_file()}},
        'source_manifests':{**prior['source_manifests'],'two_mode_parent':{'path':cfg['parent'],'sha256':cfg['parent_manifest_sha256']}},
        'initial_state':p.read(STABILITY_SMOKE/'initial_state.json'),'interpreter':sys.executable,
        'physical_model_unchanged':True,'type_II_only':True,'no_new_ODE_FEM':True}
    fingerprint=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=STABILITY_OUTPUT/fingerprint
    if b.exists():raise ValueError('Preserve interrupted attempt')
    b.mkdir(parents=True);p.write(b/'config.json',cfg);p.write(b/'provenance.json',item)
    folder=b/'execution_code';folder.mkdir()
    for f in item['code_hashes']:shutil.copyfile(p.ROOT/f,folder/Path(f).name)
    summary={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':0.,'cases':[],
        'stages':{'A':'PASS','B':'NOT_RUN','C':'NOT_RUN','D':'NOT_RUN','E':'NOT_RUN'},
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'FEM':0,'HB_phase_evaluations':0,'canonical_phase_evaluations':0,'pointwise_Jacobian':0,'Hill_eigenproblems':0}}
    p.finalize(b,item,summary);return b,item,summary


def hb_phases(U,Omega,Phi,phases):
    n=np.arange(len(U));C=np.cos(phases[:,None]*n);S=np.sin(phases[:,None]*n)
    return C@U@Phi.T,-Omega*S@(n[:,None]*U)@Phi.T,-Omega**2*C@(n[:,None]**2*U)@Phi.T


def physical_amplitude_gradient(d,Phi,U):
    q=Phi@U.sum(axis=0);amplitude,points=p.polynomial_w_max(d,q)
    peak=max(points,key=lambda r:abs(r['value']));i=peak['arm']-1;arm=d.arms[i]
    basis=arm._raw_basis(np.array([peak['s']]),0)[0]
    gradient=np.sign(peak['value'])*(basis@d.maps[i][arm.slices['w']]@Phi)
    return amplitude,np.tile(gradient,len(U)),peak


def full_massform_tangents(d,q,v,acc):
    from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import force_tangents
    Q=np.zeros((d.ndof,d.ndof));V=np.zeros_like(Q)
    for arm,P,qi,vi,ai in zip(d.arms,d.maps,d.local(q),d.local(v),d.local(acc)):
        first,second=force_tangents(arm,qi,vi,ai);Q+=P.T@first@P;V+=P.T@second@P
    return Q,V,d.mass_matrix(q)


class FullPlanarHB:
    """Newton AFT on EVERY planar Fourier coefficient and the frequency.

    Reversibility fixes phase with cosine coordinates. Velocity-dependent
    kinetic terms are even under simultaneous time reversal, so a turning
    point with all velocities zero is compatible with this fixed branch.
    The equation is the frozen mass-form residual, without mass inversion.
    """
    def __init__(self,d,Phi,omega,cfg,deadline):
        self.d=d;self.Phi=Phi;self.omega=omega;self.cfg=cfg;self.deadline=deadline;self.calls=0

    def check_time(self):
        if time.perf_counter()>self.deadline:raise TimeoutError('900 s scientific budget')

    def equations(self,U,Omega,A,jacobian=False):
        self.check_time();count=self.cfg['HB_phase_samples'];phases=np.arange(count)*2*np.pi/count
        n=np.arange(len(U));C=np.cos(phases[:,None]*n);S=np.sin(phases[:,None]*n)
        W=C.T/count;W[1:]*=2;der=-Omega*S*n;sec=-Omega**2*C*n*n
        q,v,acc=hb_phases(U,Omega,self.Phi,phases);R=[];Q=[];V=[];M=[]
        for qi,vi,ai in zip(q,v,acc):
            self.check_time();mass=self.d.mass_matrix(qi)
            residual=mass@ai+self.d.inertial_terms(qi,vi)+self.d.potential(qi)['gradient'];R.append(self.Phi.T@residual)
            self.calls+=1
            if jacobian:
                first,second,mass=full_massform_tangents(self.d,qi,vi,ai)
                Q.append(self.Phi.T@first@self.Phi);V.append(self.Phi.T@second@self.Phi);M.append(self.Phi.T@mass@self.Phi)
        eq=W@np.array(R);amplitude,gradient,peak=physical_amplitude_gradient(self.d,self.Phi,U)
        vector=np.r_[eq.ravel(),amplitude-A]
        if not jacobian:return vector
        Q,V,M=map(np.array,(Q,V,M));size=U.size
        J=sum(np.einsum('nt,tm,tij->nimj',W,basis,mat,optimize=True) for basis,mat in ((C,Q),(der,V),(sec,M))).reshape(size,size)
        Vm=v@self.d.M0@self.Phi;Am=acc@self.d.M0@self.Phi
        # Phi is a COMPLETE mass-orthonormal planar basis.
        frequency=W@(np.einsum('tij,tj->ti',V,Vm/Omega)+np.einsum('tij,tj->ti',M,2*Am/Omega))
        J=np.block([[J,frequency.reshape(-1,1)],[gradient[None,:],np.zeros((1,1))]])
        return vector,J,peak

    def solve(self,U,Omega,A,H):
        self.check_time();U=np.pad(U,((0,max(0,H+1-len(U))),(0,0)))[:H+1].copy();history=[]
        modal_scale=(A/p.polynomial_w_max(self.d,self.Phi[:,1])[0])*self.omega[1]/self.omega
        columns=np.r_[np.tile(modal_scale,H+1),self.omega[1]]
        rows=np.r_[((self.omega[None,:]**2+(np.arange(H+1)[:,None]*Omega)**2)*modal_scale).ravel(),A]
        condition=None
        for iteration in range(self.cfg['max_Newton_steps']):
            r,J,peak=self.equations(U,Omega,A,True);scaled=r/rows;relative=float(np.max(abs(scaled)))
            history.append({'iteration':iteration,'scaled_residual_max':relative,'Omega':float(Omega),'peak':peak})
            if relative<self.cfg['HB_projected_relative_target']:break
            Js=J*columns[None,:]/rows[:,None];step=np.linalg.solve(Js,-scaled)*columns
            # Bounded line search preserves the saved phase branch.
            for factor in (1.,.5,.25,.125):
                newU=U+factor*step[:-1].reshape(U.shape);newO=Omega+factor*step[-1]
                trial=self.equations(newU,newO,A)
                if np.max(abs(trial/rows))<relative and newU[1,1]>0 and newU[3,3]>0:
                    U,Omega=newU,newO;break
            else:raise ArithmeticError('Full HB Newton line search failed; no branch switch allowed')
        r,J,peak=self.equations(U,Omega,A,True);Js=J*columns[None,:]/rows[:,None]
        condition=float(np.linalg.cond(Js))
        return U,float(Omega),{'history':history,'scaled_projected_residual':float(np.max(abs(r/rows))),
            'scaled_Jacobian_condition':condition,'all_planar_modes':len(self.omega),'all_Fourier_unknowns':int(U.size+1),'harmonics':H}

    def validate(self,U,Omega,A):
        phases=np.arange(self.cfg['independent_phase_count'])*2*np.pi/self.cfg['independent_phase_count']+self.cfg['independent_phase_offset']
        q,v,acc=hb_phases(U,Omega,self.Phi,phases);res=[];energies=[];massmin=[];joints=[]
        for qi,vi,ai in zip(q,v,acc):
            self.check_time();M=self.d.mass_matrix(qi);res.append(M@ai+self.d.inertial_terms(qi,vi)+self.d.potential(qi)['gradient'])
            energies.append(self.d.energy(qi,vi));massmin.append(np.linalg.eigvalsh(M).min());self.calls+=1
        for qi,vi in zip(q[::13],v[::13]):
            r=self.d.joint_diagnostics(qi,vi);joints.append(max(r[k] for k in ('translation_difference','orientation_difference','contraction_difference')))
        res=np.array(res);scale=Omega**2*np.sqrt(U[1,1]**2+81*U[3,3]**2)
        amplitude,_,peak=physical_amplitude_gradient(self.d,self.Phi,U)
        diagnostics={'A':A,'Omega':Omega,'x':float(U[1,1]),'y':float(U[3,3]),'phase_Psi':0.,'actual_A':amplitude,
            'full_residual_absolute':float(np.max(np.linalg.norm(res,axis=1))),'full_residual_relative':float(np.max(np.linalg.norm(res,axis=1))/scale),
            'energy_relative_variation':float(np.ptp(energies)/abs(np.mean(energies))),
            'minimum_mass_eigenvalue':float(min(massmin)),'joint_compatibility':float(max(joints)),
            'independent_phase_count':len(phases),'independent_phase_offset':self.cfg['independent_phase_offset'],
            'initial_velocity_max':0.,'phase_closure_max':float(np.max(abs(hb_phases(U,Omega,self.Phi,np.array([0.,2*np.pi]))[0][0]-hb_phases(U,Omega,self.Phi,np.array([0.,2*np.pi]))[0][1]))),
            'last_two_harmonic_energy_norm_ratio':float(np.linalg.norm(U[-2:]*self.omega)/np.linalg.norm(U*self.omega)),
            'initial_mode4_quadratic_active_energy_fraction':float(self.omega[3]**2*U[3,3]**2/(self.omega[1]**2*U[1,1]**2+self.omega[3]**2*U[3,3]**2))}
        return diagnostics,{'U':U,'phases':phases,'q':q,'v':v,'acceleration':acc,'massform_residual':res,'energy':np.array(energies)}


def energy_modal_maps(d,ids,Phi,omega):
    """Canonical Pi, scaled to reference energy coordinates (omega Q,P)."""
    mass=d.M0[np.ix_(ids,ids)];zero=np.zeros_like(Phi)
    right=np.block([[Phi/omega,zero],[zero,mass@Phi]])
    left=np.block([[omega[:,None]*Phi.T@mass,np.zeros_like(Phi.T)],[np.zeros_like(Phi.T),Phi.T]])
    return right,left


def canonical_fourier(d,model,Phi,omega,U,Omega,cfg,deadline,ids,planar_Phi):
    import copy
    from scripts.lib import nlsp_instability_feasibility as f
    # This holder changes ONLY the projected index list, never the model/maps.
    holder=copy.copy(d);holder.out=ids
    action=f.QuadraticAction(holder,model.T4,model.V4)
    R,L=energy_modal_maps(d,ids,Phi,omega);phases=np.arange(cfg['tangent_phase_samples'])*2*np.pi/cfg['tangent_phase_samples']
    q,v,acc=hb_phases(U,Omega,planar_Phi,phases)
    values=[];ham=[];transform_errors=[];jac_errors=[];off_orbit_errors=[];cross=[]
    for index,(qi,vi,ai) in enumerate(zip(q,v,acc)):
        if time.perf_counter()>deadline:raise TimeoutError('Canonical Fourier budget')
        blocks=action.blocks(qi,vi,ai);Av,Ac,H,S,Sdot,J=f.action_operators(blocks)
        values.append(L@Ac@R)
        ham.append(np.linalg.norm(Ac.T@J+J@Ac)/max(np.linalg.norm(Ac),1e-300))
        if index in (7,31):
            inverse=np.linalg.inv(S);transform_errors.append(np.linalg.norm((Sdot+S@Av)@inverse-Ac)/max(np.linalg.norm(Ac),1e-300))
            actual=d.acceleration(qi,vi);exact=f.action_operators(action.blocks(qi,vi,actual));full=d.jacobian(0,np.r_[qi,vi])
            state_ids=np.r_[ids,ids+d.ndof];block=full[np.ix_(state_ids,state_ids)]
            jac_errors.append(np.linalg.norm(exact[0]-block)/max(np.linalg.norm(block),1e-300))
            off_orbit_errors.append(np.linalg.norm((Sdot+S@block)@inverse-Ac)/max(np.linalg.norm(Ac),1e-300))
            cross.append(float(max(np.max(abs(full[np.ix_(d.state_plane,d.state_out)])),np.max(abs(full[np.ix_(d.state_out,d.state_plane)])))))
    orders=np.arange(-cfg['tangent_Fourier_orders'],cfg['tangent_Fourier_orders']+1)
    values=np.array(values);F=np.einsum('mt,tij->mij',np.exp(-1j*orders[:,None]*phases)/len(phases),values,optimize=True)
    conjugacy=float(np.max(abs(F-F[::-1].conj()))/max(np.max(abs(F)),1e-300))
    diag={'modal_duality':float(np.max(abs(L@R-np.eye(len(L))))),'Hamiltonian_relative_defect':float(max(ham)),
        'time_derivative_transform_relative_defect':float(max(transform_errors)),'action_vs_full_Jacobian_relative_defect':float(max(jac_errors)),
        'orbit_Sdot_plus_full_Jacobian_relative_defect':float(max(off_orbit_errors)),
        'forbidden_cross_block_max':float(max(cross)),'Fourier_conjugacy_relative_defect':conjugacy,
        'tail_Fourier_relative_norm':float(np.linalg.norm(F[:3])+np.linalg.norm(F[-3:]))/max(np.linalg.norm(F),1e-300),
        'phase_count':len(phases),'Fourier_max_order':int(max(orders)),'all_linear_modes':len(omega)}
    diag['gate_C']=bool(max(diag[k] for k in ('modal_duality','Hamiltonian_relative_defect','time_derivative_transform_relative_defect',
        'action_vs_full_Jacobian_relative_defect','orbit_Sdot_plus_full_Jacobian_relative_defect',
        'forbidden_cross_block_max','Fourier_conjugacy_relative_defect'))<=cfg['operator_relative_gate'])
    return F,orders,diag,R,L


def build_hill(F,orders,Omega,N,k):
    """Substitution gives sum A_(n-m) z_m - i*n*Omega*z_n = lambda*z_n."""
    full=F.shape[-1]//2;ids=np.r_[np.arange(k),np.arange(k)+full];coeff=F[:,ids][:,:,ids]
    n=np.arange(-N,N+1);size=2*k;matrix=np.zeros((len(n)*size,len(n)*size),complex)
    for row,nr in enumerate(n):
        for col,nc in enumerate(n):
            hit=np.flatnonzero(orders==nr-nc)
            if len(hit):matrix[row*size:(row+1)*size,col*size:(col+1)*size]=coeff[hit[0]]
        block=slice(row*size,(row+1)*size);matrix[block,block]-=1j*nr*Omega*np.eye(size)
    return matrix,coeff


def hill_action(coeff,orders,Omega,N,z):
    # Independent convolution, with no use of the assembled dense matrix.
    n=np.arange(-N,N+1);x=z.reshape(len(n),-1);result=-1j*n[:,None]*Omega*x
    for index,nr in enumerate(n):
        for m,nc in enumerate(n):
            hit=np.flatnonzero(orders==nr-nc)
            if len(hit):result[index]+=coeff[hit[0]]@x[m]
    return result.ravel()


def fold_exponents(eigenvalues,Omega):
    return eigenvalues.real+1j*((eigenvalues.imag+Omega/2)%Omega-Omega/2)


def alias_edge_weights(eigen,folded,states,Omega):
    """lambda'=lambda-i*k*Omega requires z'_n=z_(n-k).

    Selection uses the shifted Fourier support, including discarded weight.
    Thus a tiny imaginary truncation error at +/-Omega/2 cannot erase a
    physical pair, and shifted replicas are not counted as new growth modes.
    """
    count,_,number=states.shape;N=(count-1)//2;norms=np.sum(abs(states)**2,axis=(0,1))
    shifts=np.rint((eigen.imag-folded.imag)/Omega).astype(int);weights=np.empty(number)
    for j,shift in enumerate(shifts):
        original=np.arange(-N,N+1);indices=original+shift
        # Both computed and alias-shifted support must be interior.
        bad=(abs(indices)>=N)|(abs(original)>=N)
        weights[j]=np.sum(abs(states[bad,:,j])**2)/norms[j]
    return weights,shifts


def hill_spectrum(F,orders,Omega,N,k,cfg,*,retain_reciprocal=False):
    from scipy.linalg import eig
    start=time.perf_counter();matrix,coeff=build_hill(F,orders,Omega,N,k)
    eigen,vectors=eig(matrix,check_finite=False,overwrite_a=False);folded=fold_exponents(eigen,Omega)
    # Central replicas with negligible outermost Fourier weight only.
    states=vectors.reshape(2*N+1,2*k,-1);norms=np.sum(abs(states)**2,axis=(0,1))
    edge,aliases=alias_edge_weights(eigen,folded,states,Omega)
    accepted=edge<=cfg['Hill_edge_weight_gate']
    found=np.flatnonzero(accepted)
    if not len(found):
        return {'N':N,'modal_count':k,'matrix_dimension':len(matrix),'matrix_memory_MiB':matrix.nbytes/2**20,
            'central_accepted_count':0,'status':'NO_INTERIOR_EIGENPAIR_AT_THIS_TRUNCATION','sigma':None,
            'seconds':time.perf_counter()-start}, {'eigenvalues':eigen,'folded_eigenvalues':folded,'aliases':aliases,'edge_weights':edge,'accepted':accepted}
    chosen=found[np.argmax(eigen[found].real)];lam=eigen[chosen];z=vectors[:,chosen]
    independent=hill_action(coeff,orders,Omega,N,z)
    residual=float(np.linalg.norm(independent-lam*z)/(np.linalg.norm(matrix)*np.linalg.norm(z)))
    local_residual=float(np.linalg.norm(independent-lam*z)/max(np.linalg.norm(independent)+abs(lam)*np.linalg.norm(z),1e-300))
    independent_difference=float(np.linalg.norm(independent-matrix@z)/max(np.linalg.norm(matrix@z),1e-300))
    quartet=max(float(np.min(abs(eigen-target))) for target in (lam.conjugate(),-lam,-lam.conjugate()))
    modal_weights=np.sum(abs(states[:,:,chosen])**2,axis=0);modal_weights=modal_weights[:k]+modal_weights[k:]
    modal_weights/=sum(modal_weights)
    record={'N':N,'modal_count':k,'matrix_dimension':len(matrix),'matrix_memory_MiB':matrix.nbytes/2**20,
        'central_accepted_count':len(found),'sigma':float(max(0.,lam.real)),'dominant_lambda':old.s.complex_row(folded[chosen]),
        'dominant_raw_lambda':old.s.complex_row(lam),'dominant_alias':int(aliases[chosen]),
        'independent_eigen_residual_relative':residual,'independent_eigen_residual_local_relative':local_residual,
        'independent_vs_dense_action_relative':independent_difference,
        'quartet_absolute_defect':quartet,'edge_weight':float(edge[chosen]),'modal_energy_coordinate_weights':modal_weights.tolist(),
        'seconds':time.perf_counter()-start,'resolved_spectrum_scope':'folded strip, retained modes, interior original AND shifted Fourier support; accepted counts are replicas'}
    arrays={'eigenvalues':eigen,'folded_eigenvalues':folded,'aliases':aliases,'edge_weights':edge,'accepted':accepted,
        'dominant_eigenvector':z,'dominant_independent_action':independent}
    if retain_reciprocal:
        # A negative-real unstable multiplier pairs with its inverse. Choosing
        # -conj(lambda) retains the SAME imaginary alias and Fourier convention.
        partner=found[np.argmin(abs(eigen[found]+lam.conjugate()))]
        stable=vectors[:,partner];sl=eigen[partner];sa=hill_action(coeff,orders,Omega,N,stable)
        record['reciprocal_raw_lambda']=old.s.complex_row(sl)
        record['reciprocal_pair_absolute_defect']=float(abs(sl+lam.conjugate()))
        record['reciprocal_independent_local_residual']=float(np.linalg.norm(sa-sl*stable)/max(np.linalg.norm(sa)+abs(sl)*np.linalg.norm(stable),1e-300))
        arrays.update(reciprocal_eigenvector=stable,reciprocal_lambda=sl,reciprocal_independent_action=sa)
    return record,arrays


def stability_compute(b,item,summary):
    cfg=item['config'];start=time.perf_counter();initial_seconds=summary['numerical_seconds']
    deadline=start+cfg['numerical_budget_seconds']-initial_seconds;parent=p.ROOT/cfg['parent']
    d,model,Phi,w=load_problem(p.read(parent/'provenance.json'))
    if (b/'linear_reference.npz').exists():
        with np.load(b/'linear_reference.npz') as z:outw=z['out_omega'].copy();outphi=z['out_Phi'].copy()
    else:
        outw,outphi=d.eigenpairs(plane=False)
        summary['calls']['linear_OUT_eigendecompositions']=summary['calls'].get('linear_OUT_eigendecompositions',0)+1
        np.savez_compressed(b/'linear_reference.npz',planar_omega=w,planar_Phi=Phi,out_omega=outw,out_Phi=outphi,K=d.K,M0=d.M0,plane=d.plane,out=d.out)
    summary['setup_seconds']=time.perf_counter()-start
    balance=FullPlanarHB(d,Phi,w,cfg,deadline);rows=p.read(parent/'validated_periodic_candidates.json')
    for A in ([] if summary.get('HB_refinement_completed') else cfg['amplitudes']):
        saved=next(r for r in rows if r['A']==A and r['phase_Psi']==0)
        index=sum(r['A']==A for r in rows[:rows.index(saved)])+1;label=f'A{A:.3f}'
        with np.load(parent/f'{label}_branch{index}_raw.npz') as z:sourceU=z['U'].copy()
        np.savez_compressed(b/(label+'_source.npz'),U=sourceU,Omega=saved['Omega'])
        local_start=time.perf_counter();U,Omega,newton=balance.solve(sourceU,saved['Omega'],A,cfg['HB_harmonics'][0])
        primary,raw=balance.validate(U,Omega,A);np.savez_compressed(b/(label+'_HB9.npz'),**raw);p.write(b/(label+'_HB9.json'),{**primary,'Newton':newton})
        Uc,Oc,control_newton=balance.solve(U,Omega,A,cfg['HB_harmonics'][1]);control,raw=balance.validate(Uc,Oc,A)
        np.savez_compressed(b/(label+'_HB11.npz'),**raw)
        control['Newton']=control_newton;control['HB_seconds']=time.perf_counter()-local_start
        control['source_Omega']=saved['Omega'];control['frequency_change_from_source']=Oc-saved['Omega']
        control['harmonic_frequency_difference']=abs(Oc-Omega)
        control['branch_active_relative_change']=float(max(abs(Uc[1,1]/sourceU[1,1]-1),abs(Uc[3,3]/sourceU[3,3]-1)))
        # The accepted orbit is HB11. HB9 measures convergence, not the
        # absolute residual of this finer representation.
        control['gate_B']=bool(control['full_residual_relative']<=cfg['HB_full_relative_gate']
            and control['energy_relative_variation']<=cfg['energy_relative_gate']
            and abs(Oc-Omega)<=cfg['frequency_control_gate'] and abs(control['actual_A']-A)<=cfg['amplitude_absolute_gate']
            and control['branch_active_relative_change']<=cfg['branch_relative_change_gate'] and control['minimum_mass_eigenvalue']>0
            and control['joint_compatibility']<=cfg['operator_relative_gate'] and Uc[1,1]>0 and Uc[3,3]>0)
        control['status']='HB_PASS' if control['gate_B'] else 'NOT_EVALUATED_HB_GATE'
        p.write(b/(label+'_HB11.json'),control);summary['cases'].append(control)
        summary['numerical_seconds']=initial_seconds+time.perf_counter()-start;summary['calls']['HB_phase_evaluations']=balance.calls;p.finalize(b,item,summary)
        print(json.dumps({'A':A,'HB':control},ensure_ascii=False),flush=True)
    summary['HB_refinement_completed']=True
    summary['stages']['B']='PASS' if all(r['gate_B'] for r in summary['cases']) else 'PARTIAL';p.finalize(b,item,summary)
    # All four HB costs are known before spectrum stages.
    for row in summary['cases']:
        if not row['gate_B']:continue
        balance.check_time();A=row['A'];label=f'A{A:.3f}';Omega=row['Omega']
        with np.load(b/(label+'_HB11.npz')) as z:U=z['U'].copy()
        stage=time.perf_counter()
        if (b/(label+'_OUT_Fourier.npz')).exists():
            with np.load(b/(label+'_OUT_Fourier.npz')) as ar:F=ar['F'].copy();orders=ar['orders'].copy()
            diag=row['OUT_canonical']
        else:
            F,orders,diag,R,L=canonical_fourier(d,model,outphi,outw,U,Omega,cfg,deadline,d.out,Phi)
            summary['calls']['canonical_phase_evaluations']+=cfg['tangent_phase_samples'];summary['calls']['pointwise_Jacobian']+=2
            np.savez_compressed(b/(label+'_OUT_Fourier.npz'),F=F,orders=orders,R=R,L=L,omega=outw,Phi=outphi)
        row['OUT_canonical']=diag
        if not diag['gate_C']:row['status']='OUT_INDETERMINATE';continue
        checks=[]
        for k in cfg['OUT_modal_counts']:
            for N in cfg['OUT_Hill_orders']:
                balance.check_time();check,arrays=hill_spectrum(F,orders,Omega,N,k,cfg);checks.append(check)
                summary['calls']['Hill_eigenproblems']+=1;np.savez_compressed(b/f'{label}_OUT_k{k}_N{N}.npz',**arrays)
                p.write(b/(label+'_OUT_Hill_checks.json'),checks)
        row['OUT_spectral_checks']=checks;last=checks[-1]
        if any(c.get('sigma') is None for c in (last,checks[-2],checks[-4])):
            row['status']='OUT_INDETERMINATE';row['OUT_spectral_convergence']=False;continue
        temporal=abs(last['sigma']-checks[-2]['sigma']);spatial=abs(last['sigma']-checks[-4]['sigma'])
        row['OUT_sigma']=last['sigma'];row['OUT_temporal_change']=temporal;row['OUT_modal_change']=spatial
        row['OUT_uncertainty_indicator']=max(temporal,spatial,cfg['spectral_numerical_floor'],last['quartet_absolute_defect'])
        converged=bool(max(temporal,spatial)<=cfg['spectral_control_absolute_gate'] and last['independent_eigen_residual_relative']<=cfg['Hill_residual_gate']
            and last['independent_eigen_residual_local_relative']<=cfg['Hill_residual_gate']
            and last['independent_vs_dense_action_relative']<=cfg['operator_relative_gate'] and last['quartet_absolute_defect']<=cfg['spectral_control_absolute_gate'])
        row['OUT_spectral_convergence']=converged
        row['status']=('OUT_GROWTH_CANDIDATE' if converged and row['OUT_sigma']>cfg['growth_uncertainty_factor']*row['OUT_uncertainty_indicator']
            else 'OUT_NO_GROWTH_IN_RESOLVED_SPECTRUM' if converged and row['OUT_sigma']<=cfg['spectral_numerical_floor'] else 'OUT_INDETERMINATE')
        row['OUT_modal_scope']='lowest 24 of 69 p12 OUT modes; no all-space stability claim'
        if row['status']=='OUT_GROWTH_CANDIDATE':
            row['e_folding_time']=1/row['OUT_sigma'];row['e_folding_periods']=Omega/(2*np.pi*row['OUT_sigma'])
        row['OUT_seconds']=time.perf_counter()-stage
        # A retained-growth sensitivity control uses the already saved HB9 orbit.
        with np.load(b/(label+'_HB9.npz')) as z:Up=z['U'].copy()
        if row['OUT_sigma']>cfg['spectral_numerical_floor'] and row['OUT_spectral_convergence']:
            primary_Omega=p.read(b/(label+'_HB9.json'))['Omega']
            Fp,_,dp,_,_=canonical_fourier(d,model,outphi,outw,Up,primary_Omega,cfg,deadline,d.out,Phi)
            cp,ap=hill_spectrum(Fp,orders,primary_Omega,7,24,cfg)
            summary['calls']['canonical_phase_evaluations']+=cfg['tangent_phase_samples'];summary['calls']['pointwise_Jacobian']+=2;summary['calls']['Hill_eigenproblems']+=1
            np.savez_compressed(b/(label+'_OUT_HB9_control.npz'),F=Fp,**ap)
            row['OUT_orbit_refinement_change']=abs(cp['sigma']-row['OUT_sigma']);row['OUT_HB9_control']=cp
            row['OUT_uncertainty_indicator']=max(row['OUT_uncertainty_indicator'],row['OUT_orbit_refinement_change'])
            if row['OUT_sigma']<=cfg['growth_uncertainty_factor']*row['OUT_uncertainty_indicator']:row['status']='OUT_INDETERMINATE'
        summary['numerical_seconds']=initial_seconds+time.perf_counter()-start;p.finalize(b,item,summary)
        print(json.dumps({'A':A,'OUT_status':row['status'],'sigma':row['OUT_sigma'],'uncertainty':row['OUT_uncertainty_indicator']},ensure_ascii=False),flush=True)
    # Limited PL control is deliberately last in the budget priority.
    for row in summary['cases']:
        if not row['gate_B']:continue
        balance.check_time();label=f'A{row["A"]:.3f}';Omega=row['Omega']
        with np.load(b/(label+'_HB11.npz')) as z:U=z['U'].copy()
        stage=time.perf_counter();F,orders,diag,R,L=canonical_fourier(d,model,Phi[d.plane],w,U,Omega,cfg,deadline,d.plane,Phi)
        summary['calls']['canonical_phase_evaluations']+=cfg['tangent_phase_samples'];summary['calls']['pointwise_Jacobian']+=2
        np.savez_compressed(b/(label+'_PL_Fourier.npz'),F=F,orders=orders,R=R,L=L,omega=w)
        row['PL_canonical']=diag;checks=[]
        if diag['gate_C']:
            for k in cfg['PL_modal_counts']:
                for N in cfg['PL_Hill_orders']:
                    balance.check_time();check,arrays=hill_spectrum(F,orders,Omega,N,k,cfg);checks.append(check)
                    summary['calls']['Hill_eigenproblems']+=1;np.savez_compressed(b/f'{label}_PL_k{k}_N{N}.npz',**arrays)
            row['PL_sigma_raw_including_neutral']=checks[-1]['sigma'];row['PL_checks']=checks
            # Time shift is a known neutral solution: canonical S*(qdot,qddot).
            from scripts.lib import nlsp_instability_feasibility as f
            import copy
            holder=copy.copy(d);holder.out=d.plane;action=f.QuadraticAction(holder,model.T4,model.V4)
            phases=np.arange(cfg['tangent_phase_samples'])*2*np.pi/cfg['tangent_phase_samples'];q,v,a=hb_phases(U,Omega,Phi,phases)
            neutral=np.array([L@f.action_operators(action.blocks(qi,vi,ai))[3]@np.r_[vi[d.plane],ai[d.plane]] for qi,vi,ai in zip(q,v,a)])
            N=7;k=16;nh=np.arange(-N,N+1);coef=np.einsum('nt,ti->ni',np.exp(-1j*nh[:,None]*phases)/len(phases),neutral,optimize=True)
            selected=np.r_[np.arange(k),np.arange(k)+len(w)];z=coef[:,selected].ravel();matrix,coeff=build_hill(F,orders,Omega,N,k)
            row['PL_time_shift_neutral_residual_relative']=float(np.linalg.norm(hill_action(coeff,orders,Omega,N,z))/max(Omega*np.linalg.norm(z),1e-300))
            with np.load(b/f'{label}_PL_k{k}_N{N}.npz') as ar:ev=ar['eigenvalues'];mask=ar['accepted']
            # Near-zero exponents are recorded separately, never declared growth.
            neutral_band=max(cfg['spectral_control_absolute_gate'],2*Omega*row['PL_time_shift_neutral_residual_relative'])
            ev=ev[mask];nonnull=ev[abs(ev)>neutral_band]
            row['PL_neutral_band']=neutral_band;row['PL_non_neutral_sigma']=float(max(0.,max(nonnull.real,default=0.)))
            row['PL_modal_scope']='lowest 16 of 92 p12 PL modes, neutral direction explicitly tested; incomplete planar spectrum'
            p.write(b/(label+'_PL_Hill_checks.json'),checks)
        row['PL_seconds']=time.perf_counter()-stage;summary['numerical_seconds']=initial_seconds+time.perf_counter()-start;p.finalize(b,item,summary)
    passed=[r for r in summary['cases'] if r['gate_B']]
    summary['stages'].update(C='NOT_RUN' if not passed else 'PASS' if all(r.get('OUT_canonical',{}).get('gate_C',False) for r in passed) else 'PARTIAL',
        D='NOT_RUN' if not passed else 'PASS' if all(r.get('OUT_spectral_convergence',False) for r in passed) else 'PARTIAL',
        E='NOT_RUN' if not passed else 'PASS' if all(r['status']!='OUT_INDETERMINATE' for r in passed) else 'PARTIAL')
    summary['overall']=('OUT_GROWTH_CANDIDATE_RESOLVED' if any(r['status']=='OUT_GROWTH_CANDIDATE' for r in passed)
        else 'HB_PERIODIC_BRANCH_PARTIAL' if len(passed)<4 else 'NO_OUT_GROWTH_IN_RESOLVED_SPECTRUM' if all(r['status']=='OUT_NO_GROWTH_IN_RESOLVED_SPECTRUM' for r in passed)
        else 'OUT_STABILITY_INDETERMINATE')
    summary.update(completed=True,numerical_seconds=initial_seconds+time.perf_counter()-start)
    summary['calls']['HB_phase_evaluations']+=balance.calls if summary.get('qualification_writer_correction') else 0
    p.write(b/'numerical_costs.json',{'seconds':summary['numerical_seconds'],'budget_seconds':900,'calls':summary['calls']})
    p.finalize(b,item,summary)


def stability_full_controls(b,item,summary):
    """Triggered, bounded full-p12 controls; no geometry/amplitude extension."""
    if (b/'full_space_controls.json').exists():return
    cfg=item['config'];initial=summary['numerical_seconds'];start=time.perf_counter();deadline=start+900-initial
    candidates=[r for r in summary['cases'] if r.get('status')=='OUT_GROWTH_CANDIDATE']
    if not candidates:return
    row=max(candidates,key=lambda r:r['OUT_sigma']);label=f'A{row["A"]:.3f}';Omega=row['Omega']
    policy={'trigger':'resolved growth in lowest 24 modes; assess remaining p12 modes and PL neutral directions',
        'OUT':[[48,7],[69,7],[69,9]],'PL':[[92,5],[92,7]],'largest_OUT_dimension':2622,'largest_OUT_dense_MiB':2622**2*16/2**20,
        'largest_PL_dimension':2760,'largest_PL_dense_MiB':2760**2*16/2**20,'forecast_seconds':150.,
        'remaining_budget_seconds':900-initial,'no_new_ODE_FEM':True,'all_thresholds_unchanged':True}
    p.write(b/'full_space_control_policy.json',policy)
    if 900-initial<policy['forecast_seconds']:return
    results={}
    for block in ('OUT','PL'):
        with np.load(b/(label+'_'+block+'_Fourier.npz')) as z:F=z['F'].copy();orders=z['orders'].copy()
        records=[]
        for k,N in policy[block]:
            if time.perf_counter()>deadline:raise TimeoutError('Full-space control budget')
            record,arrays=hill_spectrum(F,orders,Omega,N,k,cfg)
            np.savez_compressed(b/f'{label}_{block}_k{k}_N{N}.npz',**arrays);records.append(record)
            summary['calls']['Hill_eigenproblems']+=1;results[block]=records
            p.write(b/'full_space_controls.json',results);summary['numerical_seconds']=initial+time.perf_counter()-start;p.finalize(b,item,summary)
            print(json.dumps({'block':block,'k':k,'N':N,'sigma':record['sigma'],'seconds':record['seconds']},ensure_ascii=False),flush=True)
    last=results['OUT'][-1];modal=abs(last['sigma']-results['OUT'][-3]['sigma']);temporal=abs(last['sigma']-results['OUT'][-2]['sigma'])
    row['OUT_full_p12_controls']=results['OUT'];row['OUT_sigma']=last['sigma'];row['OUT_modal_change']=modal;row['OUT_temporal_change']=temporal
    row['OUT_uncertainty_indicator']=max(modal,temporal,row.get('OUT_orbit_refinement_change',0),cfg['spectral_numerical_floor'])
    row['OUT_modal_scope']='all 69 p12 OUT modes; temporal Hill N7/9, no p-refinement'
    row['OUT_spectral_convergence']=bool(max(modal,temporal)<=cfg['spectral_control_absolute_gate']
        and last['independent_eigen_residual_relative']<=cfg['Hill_residual_gate']
        and last['independent_eigen_residual_local_relative']<=cfg['Hill_residual_gate'] and last['quartet_absolute_defect']<=cfg['spectral_control_absolute_gate'])
    row['status']='OUT_GROWTH_CANDIDATE' if row['OUT_spectral_convergence'] and row['OUT_sigma']>cfg['growth_uncertainty_factor']*row['OUT_uncertainty_indicator'] else 'OUT_INDETERMINATE'
    row['e_folding_time']=1/row['OUT_sigma'];row['e_folding_periods']=Omega/(2*np.pi*row['OUT_sigma'])
    row['PL_full_p12_controls']=results['PL'];row['PL_modal_scope']='all 92 p12 PL modes, interior replicas only; N5/7'
    with np.load(b/f'{label}_PL_k92_N7.npz') as z:
        ev=z['folded_eigenvalues'][z['accepted']];row['PL_full_near_zero_eigenvalues']=[old.s.complex_row(e) for e in sorted(ev,key=abs)[:8]]
        nonzero=ev[abs(ev)>cfg['spectral_control_absolute_gate']];row['PL_non_neutral_sigma']=float(max(0.,max(nonzero.real,default=0.)))
    # Full time-shift solution, no modal projection error in the neutral check.
    from scripts.lib import nlsp_instability_feasibility as f
    import copy
    d,model,Phi,w=load_problem(p.read(p.ROOT/cfg['parent']/'provenance.json'))
    holder=copy.copy(d);holder.out=d.plane;action=f.QuadraticAction(holder,model.T4,model.V4)
    with np.load(b/(label+'_HB11.npz')) as z:U=z['U'].copy()
    with np.load(b/(label+'_PL_Fourier.npz')) as z:PLF=z['F'].copy();PLL=z['L'].copy()
    phases=np.arange(cfg['tangent_phase_samples'])*2*np.pi/cfg['tangent_phase_samples'];q,v,a=hb_phases(U,Omega,Phi,phases)
    neutral=np.array([PLL@f.action_operators(action.blocks(qi,vi,ai))[3]@np.r_[vi[d.plane],ai[d.plane]] for qi,vi,ai in zip(q,v,a)])
    n=np.arange(-11,12);coef=np.einsum('nt,ti->ni',np.exp(-1j*n[:,None]*phases)/len(phases),neutral,optimize=True)
    # Convolution only: do not allocate the full N11 PL matrix for this check.
    nz=coef.ravel();neutral_res=hill_action(PLF,orders,Omega,11,nz)
    row['PL_full_time_shift_neutral_residual_relative']=float(np.linalg.norm(neutral_res)/(Omega*np.linalg.norm(nz)))
    np.savez_compressed(b/(label+'_PL_neutral_direction.npz'),coefficients=coef,orders=n,residual=neutral_res)
    summary['calls']['canonical_phase_evaluations']+=len(phases)
    with np.load(b/'linear_reference.npz') as z:outphi=z['out_Phi'].copy();outw=z['out_omega'].copy()
    with np.load(b/f'{label}_OUT_k69_N9.npz') as z:zvec=z['dominant_eigenvector'].copy()
    state=zvec.reshape(19,138).sum(axis=0);position=np.zeros(d.ndof,complex);momentum=np.zeros_like(position)
    position[d.out]=outphi@(state[:69]/outw);momentum[d.out]=d.M0[np.ix_(d.out,d.out)]@outphi@state[69:]
    outaction=f.QuadraticAction(d,model.T4,model.V4);q0=Phi@U.sum(axis=0)
    blocks=outaction.blocks(q0,np.zeros(d.ndof));velocity=np.zeros_like(position)
    velocity[d.out]=np.linalg.solve(blocks['M'],momentum[d.out]-blocks['C']@position[d.out])
    parts=sum(d.reference_energy_parts(x.real,y.real)+d.reference_energy_parts(x.imag,y.imag) for x,y in [(position,velocity)])
    row['dangerous_OUT_reference_bending_fraction']=float(parts[:,0].sum()/parts.sum())
    row['dangerous_OUT_reference_torsion_fraction']=float(parts[:,1].sum()/parts.sum())
    row['dangerous_OUT_modal_weights']=last['modal_energy_coordinate_weights']
    profiles=[arm.reconstruct(P@position,np.linspace(0,arm.length,201)) for arm,P in zip(d.arms,d.maps)]
    np.savez_compressed(b/(label+'_dangerous_OUT_shape.npz'),position=position,velocity=velocity,modal_Hill_coefficients=zvec.reshape(19,138),
        phase_zero_profiles=np.array(profiles),reference_energy_parts=parts,omega=outw,Phi=outphi)
    # Exact saved-action selection and zero-base spectrum; no ODE control.
    _,_,reflection=f.second_variation(model);zero=f.action_operators(outaction.blocks(np.zeros(d.ndof),np.zeros(d.ndof)))
    R,L=energy_modal_maps(d,d.out,outphi,outw);A0=L@zero[1]@R
    expected=np.block([[np.zeros((69,69)),np.diag(outw)],[-np.diag(outw),np.zeros((69,69))]])
    zero_F=np.zeros((len(orders),138,138),complex);zero_F[list(orders).index(0)]=A0
    zr,za=hill_spectrum(zero_F,orders,Omega,7,8,cfg);summary['calls']['Hill_eigenproblems']+=1
    summary['calls']['canonical_phase_evaluations']+=2
    zero_diag={'action_vs_linear_relative':float(np.linalg.norm(A0-expected)/np.linalg.norm(expected)),'zero_base_Hill':zr,
        'odd_OUT_terms':reflection['odd_OUT_terms'],'exact_second_variation_residual_differences':{key:value['exact_residual_difference_terms'] for key,value in reflection['exact_action_residual_checks'].items()}}
    p.write(b/'zero_base_and_reflection_checks.json',zero_diag);np.savez_compressed(b/'zero_base_Hill.npz',**za)
    # Independent Fourier replay at phases absent from coefficient extraction.
    phases=np.arange(12)*2*np.pi/12+.271;qb,vb,ab=hb_phases(U,Omega,Phi,phases)
    with np.load(b/(label+'_OUT_Fourier.npz')) as z:F=z['F'].copy()
    errors=[]
    for ph,qi,vi,ai in zip(phases,qb,vb,ab):
        actual=L@f.action_operators(outaction.blocks(qi,vi,ai))[1]@R
        reconstructed=np.einsum('m,mij->ij',np.exp(1j*orders*ph),F)
        errors.append(float(np.linalg.norm(actual-reconstructed)/np.linalg.norm(actual)))
    summary['calls']['canonical_phase_evaluations']+=len(phases);row['OUT_independent_Fourier_reconstruction_relative']=max(errors)
    row['PL_neutral_band']=cfg['spectral_control_absolute_gate']
    row['full_control_seconds']=time.perf_counter()-start;summary['numerical_seconds']=initial+time.perf_counter()-start
    if row['status']!='OUT_GROWTH_CANDIDATE':summary.update(overall='OUT_STABILITY_INDETERMINATE')
    p.write(b/'full_space_controls.json',results);p.write(b/'numerical_costs.json',{'seconds':summary['numerical_seconds'],'budget_seconds':900,'calls':summary['calls']})
    p.finalize(b,item,summary)


def stability_plot(b):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summary=p.read(b/'summary.json');rows=summary['cases']
    if not rows:return
    fig,ax=plt.subplots(1,2,figsize=(10,4));A=[r['A'] for r in rows]
    ax[0].plot(A,[r['source_Omega'] for r in rows],'s:',label='saved two-active-slot approximation')
    ax[0].plot(A,[r['Omega'] for r in rows],'o-',label='full planar HB11')
    ax[0].set_ylabel('Omega [1/model time]');ax[0].legend(fontsize=7)
    ax[1].semilogy(A,[r['full_residual_relative'] for r in rows],'o-',label='independent phase residual')
    ax[1].axhline(summary.get('HB_gate',1e-7),color='gray',ls=':',label='frozen HB gate');ax[1].legend(fontsize=7);ax[1].set_ylabel('relative mass-form residual')
    for a in ax:a.set_xlabel('initial max |w| [model length]');a.grid(alpha=.2)
    fig.tight_layout()
    for ext in ('pdf','png'):fig.savefig(b/('HB_orbit_controls.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
    plt.close(fig)
    spectral=[r for r in rows if 'OUT_sigma' in r]
    if spectral:
        fig,ax=plt.subplots(figsize=(7,4))
        for r in spectral:
            ax.errorbar(r['A'],r['OUT_sigma'],yerr=r['OUT_uncertainty_indicator'],fmt='o' if r['status']=='OUT_GROWTH_CANDIDATE' else 'x',color='tab:blue')
        ax.plot([r['A'] for r in spectral],[r.get('PL_non_neutral_sigma',np.nan) for r in spectral],'s:',label='PL: resolved spectrum, neutral band excluded')
        ax.plot([],[],'o',color='tab:blue',label='OUT: resolved candidate; bars = observed uncertainty')
        if any(r['status']!='OUT_GROWTH_CANDIDATE' for r in spectral):ax.plot([],[],'x',color='tab:blue',label='OUT: unresolved / no growth')
        for r in rows:
            if 'OUT_sigma' not in r:
                ax.axvline(r['A'],color='gray',alpha=.2,lw=.8)
                ax.text(r['A'],.97,'HB gate:\nnot evaluated',transform=ax.get_xaxis_transform(),ha='center',va='top',fontsize=7,color='gray')
        ax.set_xlim(min(A)-.0005,max(A)+.0005);ax.axhline(0,color='gray',lw=.8)
        ax.set(xlabel='initial max |w| [model length]',ylabel='max Re lambda [1/model time]',title='p12, Type II; finite Hill/modal spectrum')
        ax.legend(fontsize=7);ax.grid(alpha=.2);fig.tight_layout()
        for ext in ('pdf','png'):fig.savefig(b/('Hill_growth_controls.'+ext),dpi=180,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(fig)


def stability_main(argv=None):
    parser=argparse.ArgumentParser(description='Fixed beta15 full-space HB / bounded spectral Hill audit')
    mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();summary=p.cache(b)
        if args.plot_only:stability_plot(b);p.finalize(b,p.read(b/'provenance.json'),summary)
    else:
        b,item,summary=stability_prepare()
        if not summary['completed'] and not summary['hard_stop']:
            start=time.perf_counter()
            try:
                stability_compute(b,item,summary)
                stability_full_controls(b,item,summary)
            except (TimeoutError,ArithmeticError,np.linalg.LinAlgError) as error:
                summary.update(completed=True,hard_stop=True,overall='NUMERICAL_PARTIAL' if isinstance(error,TimeoutError) else 'BLOCKED',
                    stop_reason=type(error).__name__+': '+str(error),numerical_seconds=time.perf_counter()-start)
                p.finalize(b,item,summary)
            stability_plot(b);p.finalize(b,item,summary)
    print(json.dumps(p.plain({'bundle':str(b),'summary':summary}),ensure_ascii=False))
    return 0


# One explicitly authorized dynamic validation of the saved A=.005 orbit.
# Uses the established physical Radau/DOP853 helpers, never a reduced ODE.
DYNAMIC_CONFIG=p.ROOT/'data/input/nlsp_beta15_dynamic_out_growth_validation.json'
DYNAMIC_OUTPUT=p.ROOT/'results/nlsp_beta15_dynamic_out_growth_validation'
DYNAMIC_SMOKE=p.ROOT/'results/_smoke/nlsp_beta15_dynamic_out_growth_validation'


def dynamic_prepare():
    cfg=p.read(DYNAMIC_CONFIG)
    if (cfg['degree']!=12 or cfg['beta_deg']!=15 or cfg['mu']!=.1 or cfg['A']!=.005 or cfg['phase_Psi']!=0
        or cfg['numerical_budget_seconds']>3600 or cfg['conditional_horizon_periods']>5 or cfg['new_FEM_jobs']
        or cfg['automatic_retry'] or cfg['automatic_geometry_or_amplitude_extension']):raise ValueError('Dynamic scope expansion')
    parent=p.ROOT/cfg['parent']
    if p.sha(parent/'manifest.json')!=cfg['parent_manifest_sha256']:raise ValueError('Historical Hill source changed')
    p.cache(parent);prior=p.read(parent/'provenance.json')
    for name in old.f.s.r.KERNELS:
        if p.sha(p.ROOT/name)!=prior['code_hashes'][name]:raise ValueError('Frozen physical kernel changed')
    for provenance in DYNAMIC_OUTPUT.glob('*/provenance.json') if DYNAMIC_OUTPUT.exists() else []:
        item=p.read(provenance)
        if item['config']!=cfg:raise ValueError('A different dynamic authorization already exists')
        return provenance.parent,item,p.cache(provenance.parent)
    files=[DYNAMIC_CONFIG.relative_to(p.ROOT).as_posix(),Path(__file__).resolve().relative_to(p.ROOT).as_posix(),
        'scripts/lib/nlsp_instability_feasibility.py','scripts/lib/nlsp_beta15_nonlinear_planar_resonance.py']
    item={'schema':'nlsp-beta15-dynamic-validation-v1','config':cfg,
        'source_manifests':{**prior['source_manifests'],'Hill_parent':{'path':cfg['parent'],'sha256':cfg['parent_manifest_sha256']}},
        'source_artifacts':{**prior['source_artifacts'],**{f.relative_to(p.ROOT).as_posix():p.sha(f) for f in parent.iterdir() if f.is_file()}},
        'code_hashes':{**p.code_hashes(),**{name:p.sha(p.ROOT/name) for name in files}},
        'initial_state':p.read(DYNAMIC_SMOKE/'initial_state.json'),'interpreter':sys.executable,
        'scientific_contract':'one fixed p12 Type-II orbit, stages conditional, no forcing/FEM/geometry extension'}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=DYNAMIC_OUTPUT/fp
    if b.exists():raise ValueError('Preserve interrupted attempt')
    b.mkdir(parents=True);p.write(b/'provenance.json',item);p.write(b/'config.json',cfg);folder=b/'execution_code';folder.mkdir()
    for name in item['code_hashes']:shutil.copyfile(p.ROOT/name,folder/Path(name).name)
    summary={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':p.read(DYNAMIC_SMOKE/'preflight.json')['scientific_preflight_seconds'],
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'FEM':0},'cases':{},'attempts':[],
        'stages':{key:'NOT_RUN' for key in ('A','B','C','D','E')},'deferred':cfg['deferred']}
    p.finalize(b,item,summary);return b,item,summary


def dynamic_problem(item):
    parent=p.ROOT/item['config']['parent'];hill=p.read(parent/'provenance.json')
    d,model,Phi,w=load_problem(p.read(p.ROOT/hill['config']['parent']/'provenance.json'))
    with np.load(parent/'A0.005_HB11.npz') as z:U=z['U'].copy()
    row=p.read(parent/'summary.json')['cases'][0]
    with np.load(parent/'A0.005_OUT_k69_N9.npz') as z:coeff=z['dominant_eigenvector'].reshape(19,138).copy()
    with np.load(parent/'A0.005_OUT_Fourier.npz') as z:R=z['R'].copy()
    lam=row['OUT_full_p12_controls'][-1]['dominant_raw_lambda'];lam=complex(lam['real'],lam['imag'])
    return d,model,Phi,U,row['Omega'],lam,coeff,R


def dynamic_metric(d,ids):
    n=len(ids);metric=np.zeros((2*n,2*n));metric[:n,:n]=d.K[np.ix_(ids,ids)];metric[n:,n:]=d.M0[np.ix_(ids,ids)]
    return metric


def out_polynomial_extrema(d,coordinate,Lref=1.):
    """Exact extrema of the represented local v/Phi/psi polynomials."""
    maxima=np.zeros(3);locations=[]
    for arm,qi in zip(d.arms,d.local(coordinate)):
        for k,field in enumerate(('v','Phi','psi')):
            coeff=qi[arm.slices[field]];leg=np.zeros(arm.p+1)
            for j,val in enumerate(coeff[:-1]):leg[j]+=val;leg[j+2]-=val
            leg[:2]+=coeff[-1]/2;poly=np.polynomial.legendre.Legendre(leg);roots=poly.deriv().roots()
            points=np.r_[-1.,1.,roots.real[(abs(roots.imag)<1e-9)&(abs(roots.real)<=1)]]
            value=float(np.max(abs(poly(points))));maxima[k]=max(maxima[k],value)
            locations.append({'field':field,'maximum':value})
    return maxima,float(np.max(maxima*np.array([1.,Lref,Lref]))),locations


def dynamic_initial(b,item,s,problem):
    if (b/'initial_states.npz').exists():return
    from scripts.lib import nlsp_instability_feasibility as f
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];start=time.perf_counter()
    q,v,a=hb_phases(U,Omega,Phi,np.array([0.]));initial=np.r_[q[0],v[0]];T=2*np.pi/Omega
    if abs(p.polynomial_w_max(d,q[0])[0]-.005)>2e-12 or np.max(abs(initial[d.state_out]))!=0:raise ArithmeticError('HB11 physical identity failed')
    action=f.QuadraticAction(d,model.T4,model.V4);S=f.action_operators(action.blocks(q[0],v[0],a[0]))[3]
    z0=np.linalg.solve(S,R@coeff.sum(axis=0));metric=dynamic_metric(d,d.out)
    B=np.c_[z0.real,z0.imag];Gram=B.T@metric@B;values,vectors=np.linalg.eigh(Gram);phase=vectors[0,-1]-1j*vectors[1,-1]
    rotated=phase*z0;real=rotated.real;normal=np.sqrt(real@metric@real);real/=normal;phase/=normal
    multiplier=np.exp(lam*T);ratio=float(values[0]/values[-1])
    if ratio>cfg['real_mode_Gram_ratio_gate'] or abs(multiplier.imag)/abs(multiplier)>cfg['real_multiplier_imaginary_relative_gate'] or multiplier.real>=0:
        raise ArithmeticError('Saved mode is not a resolved one-dimensional negative-real multiplier direction')
    full=np.zeros(2*d.ndof);full[d.state_out]=real
    physical,maxscale,_=out_polynomial_extrema(d,full[:d.ndof],cfg['Lref'])
    if maxscale<=1e-10:raise ArithmeticError('Displacement normalization is ill-conditioned')
    # Fix the arbitrary sign before any time-domain result is inspected.
    arm=d.arms[0];mid=arm.reconstruct(d.maps[0]@full[:d.ndof],np.array([arm.length/2]))[0]
    if mid[2]<0:full*=-1;real*=-1;phase*=-1
    delta=cfg['physical_OUT_fraction']*cfg['A']/maxscale;perturbed=initial+delta*full
    scaled,maxactual,_=out_polynomial_extrema(d,perturbed[:d.ndof],cfg['Lref'])
    phases=np.arange(16)*2*np.pi/16+.137;qv,vv,av=hb_phases(U,Omega,Phi,phases);errors=[]
    for ph,qi,vi,ai in zip(phases,qv,vv,av):
        modal=np.exp(1j*np.arange(-9,10)*ph)@coeff;canonical=R@modal
        physical_state=np.linalg.solve(f.action_operators(action.blocks(qi,vi,ai))[3],canonical)
        derivative=R@(np.exp(1j*np.arange(-9,10)*ph)@((lam+1j*np.arange(-9,10)*Omega)[:,None]*coeff))
        difference=physical_state[69:]-derivative[:69]
        errors.append(float(np.sqrt(np.vdot(difference,d.M0[np.ix_(d.out,d.out)]@difference).real)/np.sqrt(np.vdot(physical_state,metric@physical_state).real)/Omega))
    record={'T':T,'Omega':Omega,'A':cfg['A'],'actual_A':p.polynomial_w_max(d,q[0])[0],
        'lambda_raw':old.s.complex_row(lam),'sigma_Hill':float(lam.real),'multiplier':old.s.complex_row(multiplier),
        'Hill_source':'A0.005_OUT_k69_N9.npz, saved dominant eigenvector; no new eigensolve',
        'Fourier_indices':list(range(-9,10)),'canonical_mapping':'physical (eta,eta_dot)=S(t)^-1 R sum_n z_n exp(in Omega t)',
        'real_Gram':Gram.tolist(),'real_Gram_eigenvalues':values.tolist(),'real_Gram_small_large_ratio':ratio,
        'real_subspace_rank':1,'rank_qualification':'negative-real multiplier; tiny second direction is reconstruction uncertainty, not an independent physical unstable branch',
        'phase_combination':old.s.complex_row(phase),'phase_selected_before_ODE':True,'reference_energy_norm_unit_direction':float(np.sqrt(real@metric@real)),
        'physical_delta_energy_coordinates':float(delta),'initial_OUT_max_v_Phi_psi':scaled.tolist(),
        'actual_scaled_displacement_rotation_max':maxactual,'initial_OUT_reference_energy_norm':float(delta),
        'initial_OUT_velocity_M_norm':float(delta*np.sqrt(real[69:]@metric[69:,69:]@real[69:])),
        'kinematic_reconstruction_relative_error_max':max(errors),'source_eigenpair_local_residual':p.read(p.ROOT/cfg['parent']/'summary.json')['cases'][0]['OUT_full_p12_controls'][-1]['independent_eigen_residual_local_relative'],
        'external_force':0.,'all_initial_OUT_velocities_preserved':True,'no_geometric_imperfection':True}
    if max(errors)>cfg['Floquet_reconstruction_relative_gate']:raise ArithmeticError('Floquet physical reconstruction residual unresolved')
    np.savez_compressed(b/'initial_states.npz',HB_U=U,planar_Phi=Phi,base_state=initial,OUT_unit_state=real,
        spatial_state=perturbed,metric_OUT=metric,metric_PL=dynamic_metric(d,d.plane),Floquet_coefficients=coeff,Floquet_right=R,complex_initial=z0,phase_combination=phase)
    p.write(b/'Floquet_reconstruction.json',record);s['stages'].update(A='PASS',C='PASS');s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)


def dynamic_remaining(s,cfg):return cfg['numerical_budget_seconds']-s['numerical_seconds']


def dynamic_timing(d,cfg,control=False,full=False):
    setting=p.settings(d,cfg['A'],cfg)
    if control:
        setting['rtol']*=cfg['temporal_tolerance_factor'];setting['atol']*=cfg['temporal_tolerance_factor']
    if full:
        setting=p.settings(d,cfg['A'],{**cfg,'rtol':cfg['full_rtol'],'atol_relative':cfg['full_atol_relative']})
        setting['max_step']*=cfg['full_max_step_factor']
    return setting


def dynamic_planar_check(b,item,s,problem):
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];T=2*np.pi/Omega
    with np.load(b/'initial_states.npz') as z:initial=z['base_state'].copy();metric=z['metric_PL'].copy()
    times=np.linspace(0,T,cfg['samples_per_period']+1)
    for label,control in [('planar_primary_1T',False),('planar_refined_1T',True)]:
        if label not in s['cases']:
            p.integrate_attempt(b,item,s,label,d,initial,times,dynamic_timing(d,cfg,control=control),
                time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=True)
            if s['hard_stop']:return
    if (b/'planar_one_period_checks.json').exists():return
    start=time.perf_counter();tr=[]
    for label in ('planar_primary_1T','planar_refined_1T'):
        with np.load(b/label/'trajectory.npz') as z:tr.append(z['state'].copy())
    scale=analysis.energy_norm(initial[d.state_plane],metric);hbq,hbv,_=hb_phases(U,Omega,Phi,Omega*times);hb=np.c_[hbq,hbv]
    records=[]
    for label,state in zip(('planar_primary_1T','planar_refined_1T'),tr):
        error=state[:,d.state_plane]-hb[:,d.state_plane];closure=state[-1,d.state_plane]-initial[d.state_plane]
        record={'case':label,'HB_difference_energy_relative_max':float(analysis.energy_norm(error,metric).max()/scale),
            'closure_energy_relative':float(analysis.energy_norm(closure,metric)/scale),
            'closure_coefficient_absolute_max':float(np.max(abs(closure))),
            'closure_position_K_norm':float(analysis.energy_norm(closure[:92],metric[:92,:92])),
            'closure_velocity_M_norm':float(analysis.energy_norm(closure[92:],metric[92:,92:])),
            'energy_drift':s['cases'][label]['energy_relative_drift_max']}
        records.append(record)
    temporal=float(analysis.energy_norm(tr[0][:,d.state_plane]-tr[1][:,d.state_plane],metric).max()/scale)
    indices=np.unique(np.r_[0,np.linspace(0,len(times)-1,33,dtype=int),len(times)-1]);masses=[];joints=[];leaks=[]
    for x in tr[-1][indices]:
        masses.append(float(np.linalg.eigvalsh(d.mass_matrix(x[:d.ndof])).min()));jj=d.joint_diagnostics(x[:d.ndof],x[d.ndof:])
        joints.append(max(jj[k] for k in ('translation_difference','orientation_difference','contraction_difference')));leaks.append(float(np.max(abs(d.rhs(0,x)[d.state_out]))))
    check={'records':records,'temporal_difference_energy_relative_max':temporal,'reference_PL_energy_norm':float(scale),
        'mass_min_eigenvalue':min(masses),'joint_compatibility_max':max(joints),'planar_rhs_OUT_max':max(leaks),
        'fixed_period':T,'period_not_fitted':True,'authoritative_base':'planar_refined_1T','adaptive_solver_path_restart_not_identical':True,
        'gate_B':bool(max(r['HB_difference_energy_relative_max'] for r in records)<=cfg['planar_orbit_energy_relative_gate']
            and max(r['closure_energy_relative'] for r in records)<=cfg['planar_orbit_energy_relative_gate']
            and max(r['energy_drift'] for r in records)<=cfg['energy_drift_relative_gate'] and temporal<=cfg['temporal_energy_relative_gate']
            and min(masses)>0 and max(joints)<=cfg['joint_absolute_gate'] and max(leaks)<=cfg['reflection_absolute_gate'])}
    np.savez_compressed(b/'planar_one_period_comparison.npz',times=times,primary=tr[0],refined=tr[1],HB=hb,metric=metric)
    p.write(b/'planar_one_period_checks.json',check);s['numerical_seconds']+=time.perf_counter()-start
    s['stages']['B']='PASS' if check['gate_B'] else 'FAIL'
    if not check['gate_B']:s.update(completed=True,hard_stop=True,overall='PLANAR_ORBIT_DYNAMIC_VALIDATION_FAILED',stop_reason='One-period orbit/time/safety gate failed')
    p.finalize(b,item,s)


def join_dynamic_base(b,labels,name):
    folder=b/name
    if (folder/'dense.npz').exists():return folder
    folder.mkdir();arrays={key:[] for key in ('t_old','t','y_old','Q')};times=[];states=[]
    for j,label in enumerate(labels):
        with np.load(b/label/'dense.npz') as z:
            for key in arrays:arrays[key].append(z[key].copy())
        with np.load(b/label/'trajectory.npz') as z:
            times.append(z['times'][int(j>0):].copy());states.append(z['state'][int(j>0):].copy())
    np.savez_compressed(folder/'dense.npz',**{key:np.concatenate(values) for key,values in arrays.items()})
    np.savez_compressed(folder/'trajectory.npz',times=np.concatenate(times),state=np.concatenate(states))
    p.write(folder/'join_provenance.json',{'source_cases':labels,'physical_states_not_renormalized':True,'prefix_not_recomputed':True,
        'join_t':float(arrays['t'][0][-1]),'join_y_difference':float(np.max(abs(arrays['y_old'][1][0]-(arrays['y_old'][0][-1]+arrays['Q'][0][-1].sum(axis=1))))) if len(labels)>1 else 0.})
    return folder


def dynamic_tangent_check(b,item,s,problem):
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];T=2*np.pi/Omega
    with np.load(b/'initial_states.npz') as z:basis=z['OUT_unit_state'].copy()[:,None];metric=z['metric_OUT'].copy();initial=z['base_state'].copy();plmetric=z['metric_PL'].copy();complex_initial=z['complex_initial'].copy();phase=z['phase_combination'].item()
    if 'planar_continue_T_to_3T' not in s['cases']:
        with np.load(b/'planar_refined_1T/accepted_times.npz') as z:restart=z['actual_final_state'].copy()
        times=np.linspace(T,3*T,2*cfg['samples_per_period']+1)
        p.integrate_attempt(b,item,s,'planar_continue_T_to_3T',d,restart,times,dynamic_timing(d,cfg,control=True),
            time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=True,t_start=T)
        if s['hard_stop']:return
    base=join_dynamic_base(b,['planar_refined_1T','planar_continue_T_to_3T'],'planar_combined_3T')
    tc={**cfg,'rtol':cfg['tangent_rtol'],'atol_relative':cfg['tangent_atol_relative']};local_item={**item,'config':tc}
    for label,periods,control in [('tangent_primary_3T',3,False),('tangent_refined_1T',1,True)]:
        if label in s['cases']:continue
        config={**tc}
        if control:config['rtol']*=cfg['temporal_tolerance_factor'];config['atol_relative']*=cfg['temporal_tolerance_factor']
        times=np.linspace(0,periods*T,periods*cfg['samples_per_period']+1)
        analysis.tangent_attempt(b,{**item,'config':config},s,label,d,base,times,basis,metric,
            time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],method=cfg['tangent_method'],sampled_dense_only=True)
        if s['hard_stop']:return
    if (b/'tangent_growth_checks.json').exists():return
    start=time.perf_counter()
    with np.load(b/'tangent_primary_3T/tangent.npz') as z:times=z['times'].copy();vectors=z['vectors'][:,:,0].copy()
    with np.load(b/'tangent_refined_1T/tangent.npz') as z:refined=z['vectors'][:,:,0].copy()
    temporal=float(analysis.energy_norm(vectors[:len(refined)]-refined,metric).max())
    strobe=np.arange(4)*cfg['samples_per_period'];X=vectors[strobe];ts=times[strobe]
    expected=np.real(np.exp(lam*ts[:,None])*phase*complex_initial[None,:])
    errors=analysis.energy_norm(X-expected,metric)/analysis.energy_norm(expected,metric)
    coefficients=X@metric@basis[:,0];ampl=np.abs(coefficients/coefficients[0]);slope=float(np.polyfit(ts,np.log(ampl),1)[0])
    refined_coeff=float(refined[-1]@metric@basis[:,0]);primary_sigma=float(np.log(abs(coefficients[1]/coefficients[0]))/T)
    refined_sigma=float(np.log(abs(refined_coeff/coefficients[0]))/T);sigma_unc=abs(primary_sigma-refined_sigma)
    with np.load(base/'trajectory.npz') as z:pl=z['state'].copy()
    hq,hv,_=hb_phases(U,Omega,Phi,Omega*times);HB=np.c_[hq,hv];plscale=analysis.energy_norm(initial[d.state_plane],plmetric)
    drift=float(analysis.energy_norm(pl[:,d.state_plane]-HB[:,d.state_plane],plmetric).max()/plscale)
    record={'periods':3,'stroboscopic_times':ts.tolist(),'signed_Floquet_coefficients':coefficients.tolist(),
        'stroboscopic_amplification':ampl.tolist(),'Hill_stroboscopic_amplification':np.exp(lam.real*ts).tolist(),
        'stroboscopic_Floquet_vector_relative_errors':errors.tolist(),'sigma_dynamic':slope,'sigma_Hill':float(lam.real),
        'sigma_difference':abs(slope-lam.real),'sigma_temporal_uncertainty':sigma_unc,'tangent_temporal_energy_error':temporal,
        'planar_HB_difference_3T_energy_relative_max':drift,'negative_multiplier_signs_correct':bool(np.all(coefficients*(-1.)**np.arange(4)>0)),
        'gate_D':bool(np.max(errors)<=cfg['Floquet_stroboscopic_relative_gate'] and abs(slope-lam.real)<=cfg['growth_exponent_absolute_gate']
            and sigma_unc<=cfg['growth_temporal_absolute_gate'] and temporal<=cfg['temporal_energy_relative_gate'] and drift<=cfg['planar_orbit_energy_relative_gate']
            and np.all(coefficients*(-1.)**np.arange(4)>0)),
        'base_source':'nonlinear planar ODE accepted cubic dense output, not HB arrays'}
    p.write(b/'tangent_growth_checks.json',record);np.savez_compressed(b/'tangent_growth_comparison.npz',times=times,tangent=vectors,
        stroboscopic_times=ts,stroboscopic_actual=X,stroboscopic_Hill=expected,metric=metric,base_state=pl,HB_state=HB)
    s['numerical_seconds']+=time.perf_counter()-start;s['stages']['D']='PASS' if record['gate_D'] else 'FAIL'
    if not record['gate_D']:s.update(completed=True,hard_stop=True,overall='TANGENT_HILL_MISMATCH',stop_reason='Dynamic Floquet growth/time/background gate failed')
    else:s['overall']='TANGENT_GROWTH_CONFIRMED_FULL_PENDING'
    p.finalize(b,item,s)


def dynamic_full_forecast(b,item,s,problem):
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];start=time.perf_counter()
    with np.load(b/'initial_states.npz') as z:x=z['spatial_state'].copy()
    with np.load(b/'tangent_growth_comparison.npz') as z:
        indices=np.linspace(0,len(z['times'])-1,10,dtype=int);probe=z['base_state'][indices].copy()
        delta=p.read(b/'Floquet_reconstruction.json')['physical_delta_energy_coordinates']
        probe[:,d.state_out]=delta*z['tangent'][indices]
    timing={}
    for name in ('rhs','jacobian','safety'):
        fn=getattr(d,name);began=time.perf_counter()
        for state in probe:
            args=(state[:d.ndof],) if name=='safety' else (0,state)
            fn(*args)
        timing[name]=(time.perf_counter()-began)/10
    historical=p.ROOT/'results/nlsp_construction_C_modal_sensitivity/e7f729a57a8cdfcd/C_optimal_finite/integration.json';oldcost=p.read(historical)
    T=2*np.pi/Omega;rate=oldcost['nfev']/oldcost['target_end'];estimate=3*T*rate*(timing['rhs']+timing['safety'])
    estimate+=oldcost['nlu']/oldcost['target_end']*3*T*.01+10.
    forecast=cfg['cost_forecast_factor']*estimate;remaining=dynamic_remaining(s,cfg)
    record={'pointwise_seconds':timing,'historical_cost_source':str(historical.relative_to(p.ROOT)),
        'historical_nfev_per_model_time':rate,'forecast_3T_seconds':forecast,'remaining_seconds':remaining,
        'reserve_seconds':cfg['cost_reserve_seconds'],'affordable':forecast+cfg['cost_reserve_seconds']<=remaining,
        'no_new_ODE_for_profile':True,'forecast_not_guaranteed':True,
        'profile_samples':'ten DIFFERENT saved planar states plus their tangent perturbations; avoid identical-state cache timing'}
    p.write(b/'full_spatial_cost_forecast.json',record);s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s);return record


def dynamic_temporal_step_control(b,item,s,problem):
    """One stronger step check after equal-max-step tolerance comparisons.

    Frozen error gates are unchanged. This is a different numerical level,
    not a retry, a new amplitude or a recomputation of an accepted prefix.
    """
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    if (b/'temporal_step_control.json').exists() or s['stages']['D']!='PASS':return
    cfg=item['config'];d,model,Phi,U,Omega,lam,coeff,R=problem;T=2*np.pi/Omega
    with np.load(b/'initial_states.npz') as z:initial=z['base_state'].copy();basis=z['OUT_unit_state'].copy()[:,None];metric=z['metric_OUT'].copy();plmetric=z['metric_PL'].copy()
    policy={'trigger':'tighter tolerances used the same limiting native maximum step; strengthen time resolution',
        'additional_planar_horizon_T':1,'additional_tangent_horizon_T':1,'max_step_factor':.5,'error_gates_unchanged':True,
        'remaining_budget_seconds':dynamic_remaining(s,cfg),'no_failed_scientific_run_retried':True}
    p.write(b/'temporal_step_control_policy.json',policy);times=np.linspace(0,T,cfg['samples_per_period']+1)
    label='planar_half_step_1T'
    if label not in s['cases']:
        timing=dynamic_timing(d,cfg,control=True);timing['max_step']*=.5
        p.integrate_attempt(b,item,s,label,d,initial,times,timing,time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=True)
        if s['hard_stop']:return
    label='tangent_half_step_1T'
    if label not in s['cases']:
        tc={**cfg,'rtol':cfg['tangent_rtol']*cfg['temporal_tolerance_factor'],'atol_relative':cfg['tangent_atol_relative']*cfg['temporal_tolerance_factor']}
        maxstep=.5*s['cases']['tangent_primary_3T']['max_step']
        analysis.tangent_attempt(b,{**item,'config':tc},s,label,d,b/'planar_combined_3T',times,basis,metric,
            time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],method=cfg['tangent_method'],sampled_dense_only=True,max_step=maxstep)
        if s['hard_stop']:return
    start=time.perf_counter()
    with np.load(b/'planar_half_step_1T/trajectory.npz') as z:pl=z['state'].copy()
    with np.load(b/'planar_refined_1T/trajectory.npz') as z:pref=z['state'].copy()
    with np.load(b/'tangent_half_step_1T/tangent.npz') as z:tangent=z['vectors'][:,:,0].copy()
    with np.load(b/'tangent_primary_3T/tangent.npz') as z:tref=z['vectors'][:len(times),:,0].copy()
    plerr=float(analysis.energy_norm(pl[:,d.state_plane]-pref[:,d.state_plane],plmetric).max()/analysis.energy_norm(initial[d.state_plane],plmetric))
    terr=float(analysis.energy_norm(tangent-tref,metric).max());normal=basis[:,0]
    slope=lambda Z:float(np.log(abs(Z[-1]@metric@normal)/abs(Z[0]@metric@normal))/T)
    sigerr=abs(slope(tangent)-slope(tref));record={'planar_energy_relative_difference':plerr,'tangent_energy_difference':terr,
        'sigma_one_period_half_step':slope(tangent),'sigma_one_period_primary':slope(tref),'sigma_step_difference':sigerr,
        'gate':bool(plerr<=cfg['temporal_energy_relative_gate'] and terr<=cfg['temporal_energy_relative_gate'] and sigerr<=cfg['growth_temporal_absolute_gate'])}
    p.write(b/'temporal_step_control.json',record);s['numerical_seconds']+=time.perf_counter()-start
    if not record['gate']:s.update(hard_stop=True,completed=True,overall='DYNAMIC_NUMERICAL_PARTIAL',stop_reason='Independent half-step temporal gate failed')
    p.finalize(b,item,s)


def dynamic_full_check(b,item,s,problem):
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];T=2*np.pi/Omega
    with np.load(b/'initial_states.npz') as z:initial=z['spatial_state'].copy();basis=z['OUT_unit_state'].copy();metric=z['metric_OUT'].copy();plmetric=z['metric_PL'].copy()
    delta=p.read(b/'Floquet_reconstruction.json')['physical_delta_energy_coordinates'];label='spatial_full_1percent_3T'
    if label not in s['cases']:
        forecast=dynamic_full_forecast(b,item,s,problem)
        if not forecast['affordable']:
            s.update(completed=True,overall='TANGENT_GROWTH_CONFIRMED_FULL_PENDING',stop_reason='Full spatial 3T forecast exceeds remaining hard budget');p.finalize(b,item,s);return
        times=np.linspace(0,3*T,3*cfg['samples_per_period']+1)
        p.integrate_attempt(b,item,s,label,d,initial,times,dynamic_timing(d,cfg,full=True),
            time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=False,dense_archive=False)
        if s['hard_stop']:return
    start=time.perf_counter()
    with np.load(b/label/'trajectory.npz') as z:times=z['times'].copy();full=z['state'].copy()
    with np.load(b/'tangent_growth_comparison.npz') as z:tan=z['tangent'].copy();base=z['base_state'].copy()
    observed=full[:,d.state_out]/delta;err=analysis.energy_norm(observed-tan,metric);norm=analysis.energy_norm(tan,metric)
    plchange=analysis.energy_norm(full[:,d.state_plane]-base[:,d.state_plane],plmetric)
    strobe=np.arange(4)*cfg['samples_per_period'];coefficients=observed[strobe]@metric@basis;ampl=np.abs(coefficients/coefficients[0])
    sigma=float(np.polyfit(times[strobe],np.log(ampl),1)[0]);phys=[];tphys=[]
    for state,z in zip(full,tan):
        phys.append(out_polynomial_extrema(d,state[:d.ndof],cfg['Lref'])[0]);qt=np.zeros(d.ndof);qt[d.out]=delta*z[:69]
        tphys.append(out_polynomial_extrema(d,qt,cfg['Lref'])[0])
    physical=np.array(phys);pred=np.array(tphys);relative=float(np.max(err)/np.max(norm))
    record={'periods':3,'initial_physical_OUT_scale':cfg['A']*cfg['physical_OUT_fraction'],'delta_energy_coordinates':delta,
        'OUT_tangent_relative_error_max':relative,'OUT_tangent_absolute_energy_error_max':float(delta*err.max()),
        'stroboscopic_amplification':ampl.tolist(),'sigma_NL_projected_finite_interval':sigma,
        'growth_exponent_qualification':'finite-interval log projection, not nonlinear Floquet exponent or asymptotic rate',
        'maximum_v_Phi_psi':physical.max(axis=0).tolist(),'maximum_v_over_A':float(physical[:,0].max()/cfg['A']),
        'planar_backreaction_energy_norm_max':float(plchange.max()),'energy_drift':s['cases'][label]['energy_relative_drift_max'],
        'finite_delta_order':'NOT_TESTED_ONE_AMPLITUDE','numerical_floor_declared_relative':cfg['numerical_floor_relative'],
        'gate_E':bool(relative<=cfg['tangent_full_relative_gate'] and s['cases'][label]['energy_relative_drift_max']<=cfg['energy_drift_relative_gate']
            and np.max(abs(ampl-np.exp(lam.real*times[strobe])))<=cfg['Floquet_stroboscopic_relative_gate'])}
    np.savez_compressed(b/'full_tangent_comparison.npz',times=times,OUT_nonlinear=observed,OUT_tangent=tan,OUT_error=err,
        planar_backreaction=plchange,physical_OUT_extrema=physical,tangent_physical_OUT_extrema=pred,stroboscopic_coefficients=coefficients)
    p.write(b/'full_nonlinear_checks.json',record);s['numerical_seconds']+=time.perf_counter()-start;s['stages']['E']='PASS' if record['gate_E'] else 'FAIL'
    s.update(overall='DYNAMIC_OUT_GROWTH_CONFIRMED_P12' if record['gate_E'] else 'FULL_NONLINEAR_TANGENT_MISMATCH',hard_stop=not record['gate_E'])
    p.finalize(b,item,s)


def dynamic_plot(b):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    def save(fig,name):
        fig.tight_layout()
        for ext in ('pdf','png'):fig.savefig(b/(name+'.'+ext),dpi=170,metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else None)
        plt.close(fig)
    if (b/'planar_one_period_comparison.npz').exists():
        with np.load(b/'planar_one_period_comparison.npz') as z:t=z['times'];a=z['refined'];hb=z['HB'];metric=z['metric']
        fig,ax=plt.subplots(1,2,figsize=(10,4));rec=p.read(b/'Floquet_reconstruction.json');T=rec['T']
        if (b/'physical_observation_maps.npz').exists():
            with np.load(b/'physical_observation_maps.npz') as z:maps=z['midpoint_maps'];ids=z['state_plane']
            ndof=a.shape[1]//2
            for j in range(2):
                row=maps[j,1];ax[0].plot(t/T,a[:,:ndof]@row,label=f'ODE w{j+1}(L{j+1}/2)')
                ax[0].plot(t/T,hb[:,:ndof]@row,'--',label=f'HB11, arm {j+1}')
            from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
            error=analysis.energy_norm((a-hb)[:,ids],metric)/analysis.energy_norm(a[0,ids],metric)
            ax[0].set(ylabel='midpoint w [model length]',xlabel='t/T');ax[0].legend(fontsize=7)
            ax[1].semilogy(t/T,np.maximum(error,1e-18));ax[1].set(xlabel='t/T',ylabel='full PL relative reference-energy difference')
        else:
            idx=int(np.argmax(np.max(abs(hb[:,:hb.shape[1]//2]),axis=0)))
            ax[0].plot(t/T,a[:,idx],label='nonlinear planar ODE');ax[0].plot(t/T,hb[:,idx],'--',label='HB11')
            ax[0].set(ylabel=f'global resting-mass coordinate {idx}',xlabel='t/T');ax[0].legend(fontsize=8)
            error=np.linalg.norm(a-hb,axis=1);ax[1].semilogy(t/T,np.maximum(error,1e-18));ax[1].set(xlabel='t/T',ylabel='full coefficient difference [resting-mass state units]')
        save(fig,'planar_orbit_validation')
    if (b/'tangent_growth_checks.json').exists():
        r=p.read(b/'tangent_growth_checks.json');T=p.read(b/'Floquet_reconstruction.json')['T'];fig,ax=plt.subplots(figsize=(7,4))
        if (b/'extension_5T_checks.json').exists():
            ext=p.read(b/'extension_5T_checks.json');r={**r,'stroboscopic_times':np.arange(6)*T,
                'stroboscopic_amplification':ext['stroboscopic_tangent_amplification']}
        tau=np.array(r['stroboscopic_times'])/T;ax.plot(tau,r['stroboscopic_amplification'],'o',label='direct tangent, projected magnitude')
        fine=np.linspace(0,max(tau),150);ax.plot(fine,np.exp(r['sigma_Hill']*T*fine),'--',label='Hill exp(sigma t)')
        ax.set(xlabel='stroboscopic t/T',ylabel='Floquet component amplification',title='Signed component changes sign after every period');ax.legend();save(fig,'tangent_Floquet_validation')
    if (b/'full_tangent_comparison.npz').exists():
        file=b/'extension_5T_comparison.npz' if (b/'extension_5T_comparison.npz').exists() else b/'full_tangent_comparison.npz'
        with np.load(file) as z:t=z['times'];ph=z['physical_OUT_extrema'];th=z['tangent_physical_OUT_extrema'];obs=z['OUT_nonlinear'];tan=z['tangent'] if 'tangent' in z.files else z['OUT_tangent']
        from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
        with np.load(b/'initial_states.npz') as z:metric=z['metric_OUT']
        T=p.read(b/'Floquet_reconstruction.json')['T'];fig,ax=plt.subplots(1,2,figsize=(10,4))
        ax[0].plot(t/T,ph[:,0]/.005,label='full max |v| / A');ax[0].plot(t/T,th[:,0]/.005,'--',label='tangent');ax[0].set(xlabel='t/T',ylabel='spatial displacement / A');ax[0].legend(fontsize=8)
        n0=analysis.energy_norm(tan[0],metric);ax[1].plot(t/T,analysis.energy_norm(obs,metric)/n0,label='full OUT / initial norm');ax[1].plot(t/T,analysis.energy_norm(tan,metric)/n0,'--',label='tangent')
        ax[1].set(xlabel='t/T',ylabel='reference-energy norm ratio');ax[1].legend(fontsize=8);save(fig,'nonlinear_spatial_validation')


def dynamic_optional_controls(b,item,s,problem):
    """Only affordable 5T continuation and one authorized half-scale 1T check."""
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    if s['overall']!='DYNAMIC_OUT_GROWTH_CONFIRMED_P12':return
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];T=2*np.pi/Omega
    with np.load(b/'initial_states.npz') as z:
        initial=z['base_state'].copy();unit=z['OUT_unit_state'].copy();metric=z['metric_OUT'].copy();plmetric=z['metric_PL'].copy();z0=z['complex_initial'].copy();phase=z['phase_combination'].item()
    delta=p.read(b/'Floquet_reconstruction.json')['physical_delta_energy_coordinates']
    fullcost=s['cases']['spatial_full_1percent_3T']['seconds']
    forecast=cfg['cost_forecast_factor']*(s['cases']['planar_continue_T_to_3T']['seconds']+s['cases']['tangent_primary_3T']['seconds']*2/3+fullcost*2/3+10)
    decision={'forecast_all_5T_continuations_seconds':forecast,'remaining_seconds':dynamic_remaining(s,cfg),
        'authorized_maximum_periods':5,'gate_5T':forecast+cfg['cost_reserve_seconds']<=dynamic_remaining(s,cfg),
        'prefixes_preserved':True,'initial_Floquet_direction_not_renormalized':True}
    p.write(b/'optional_control_decision.json',decision);p.finalize(b,item,s)
    if decision['gate_5T']:
        label='planar_continue_3T_to_5T'
        if label not in s['cases']:
            with np.load(b/'planar_continue_T_to_3T/accepted_times.npz') as z:restart=z['actual_final_state'].copy()
            p.integrate_attempt(b,item,s,label,d,restart,np.linspace(3*T,5*T,2*cfg['samples_per_period']+1),dynamic_timing(d,cfg,control=True),
                time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=True,t_start=3*T)
            if s['hard_stop']:return
        base=join_dynamic_base(b,['planar_refined_1T','planar_continue_T_to_3T',label],'planar_combined_5T')
        label='tangent_continue_3T_to_5T'
        if label not in s['cases']:
            with np.load(b/'tangent_primary_3T/restart_state.npz') as z:restart=z['actual_final_matrix'].copy()
            tc={**cfg,'rtol':cfg['tangent_rtol'],'atol_relative':cfg['tangent_atol_relative']}
            analysis.tangent_attempt(b,{**item,'config':tc},s,label,d,base,np.linspace(3*T,5*T,2*cfg['samples_per_period']+1),unit[:,None],metric,
                time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],method=cfg['tangent_method'],sampled_dense_only=True,t_start=3*T,initial_matrix=restart)
            if s['hard_stop']:return
        with np.load(b/'tangent_primary_3T/tangent.npz') as z:t0=z['times'].copy();tan0=z['vectors'][:,:,0].copy()
        with np.load(b/label/'tangent.npz') as z:times=np.r_[t0,z['times'][1:]];tan=np.r_[tan0,z['vectors'][1:,:,0]]
        strobe=np.arange(6)*cfg['samples_per_period'];expected=np.real(np.exp(lam*times[strobe,None])*phase*z0[None,:])
        error=analysis.energy_norm(tan[strobe]-expected,metric)/analysis.energy_norm(expected,metric)
        with np.load(base/'trajectory.npz') as z:pl=z['state'].copy()
        hq,hv,_=hb_phases(U,Omega,Phi,Omega*times);HB=np.c_[hq,hv]
        drift=float(analysis.energy_norm(pl[:,d.state_plane]-HB[:,d.state_plane],plmetric).max()/analysis.energy_norm(initial[d.state_plane],plmetric))
        gate=bool(max(error)<=cfg['Floquet_stroboscopic_relative_gate'] and drift<=cfg['planar_orbit_energy_relative_gate'])
        p.write(b/'extension_5T_tangent_gate.json',{'errors':error.tolist(),'PL_HB_relative_max':drift,'gate':gate})
        if not gate:
            s.update(hard_stop=True,overall='DYNAMIC_NUMERICAL_PARTIAL',stop_reason='5T continuation gate failed; 3T confirmation retained');p.finalize(b,item,s);return
        label='spatial_full_continue_3T_to_5T'
        if label not in s['cases']:
            with np.load(b/'spatial_full_1percent_3T/accepted_times.npz') as z:restart=z['actual_final_state'].copy()
            p.integrate_attempt(b,item,s,label,d,restart,np.linspace(3*T,5*T,2*cfg['samples_per_period']+1),dynamic_timing(d,cfg,full=True),
                time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=False,dense_archive=False,t_start=3*T)
            if s['hard_stop']:return
        start=time.perf_counter()
        with np.load(b/'spatial_full_1percent_3T/trajectory.npz') as z:full0=z['state'].copy()
        with np.load(b/label/'trajectory.npz') as z:full=np.r_[full0,z['state'][1:]]
        observed=full[:,d.state_out]/delta;relative=float(analysis.energy_norm(observed-tan,metric).max()/analysis.energy_norm(tan,metric).max())
        signed=observed[strobe]@metric@unit;signed_tan=tan[strobe]@metric@unit
        physical=np.array([out_polynomial_extrema(d,x[:d.ndof],cfg['Lref'])[0] for x in full])
        predicted=[]
        for zz in tan:
            coordinate=np.zeros(d.ndof);coordinate[d.out]=delta*zz[:69];predicted.append(out_polynomial_extrema(d,coordinate,cfg['Lref'])[0])
        energy=np.array([d.energy(x[:d.ndof],x[d.ndof:]) for x in full]);edrift=float(np.max(abs(energy/energy[0]-1)))
        record={'periods':5,'stroboscopic_nonlinear_amplification':np.abs(signed/signed[0]).tolist(),
            'stroboscopic_tangent_amplification':np.abs(signed_tan/signed_tan[0]).tolist(),
            'sigma_tangent_dynamic':float(np.polyfit(times[strobe],np.log(abs(signed_tan/signed_tan[0])),1)[0]),
            'sigma_NL_projected_finite_interval':float(np.polyfit(times[strobe],np.log(abs(signed/signed[0])),1)[0]),
            'tangent_full_error_relative':relative,'full_energy_relative_drift':edrift,'maximum_v_Phi_psi':physical.max(axis=0).tolist(),
            'gate':bool(relative<=cfg['tangent_full_relative_gate'] and edrift<=cfg['energy_drift_relative_gate'])}
        np.savez_compressed(b/'extension_5T_comparison.npz',times=times,tangent=tan,OUT_nonlinear=observed,base_state=pl,
            physical_OUT_extrema=physical,tangent_physical_OUT_extrema=np.array(predicted),stroboscopic_times=times[strobe])
        p.write(b/'extension_5T_checks.json',record);s['numerical_seconds']+=time.perf_counter()-start
        if not record['gate']:s.update(hard_stop=True,overall='DYNAMIC_NUMERICAL_PARTIAL',stop_reason='5T nonlinear comparison failed; 3T confirmation retained');p.finalize(b,item,s);return
        p.finalize(b,item,s)
    # A single SHORT second perturbation level; never a new planar amplitude.
    remaining=dynamic_remaining(s,cfg);half_forecast=cfg['cost_forecast_factor']*(fullcost/3+10)
    decision.update(half_1T_forecast_seconds=half_forecast,remaining_after_extension_seconds=remaining,
        half_1T_affordable=half_forecast+cfg['cost_reserve_seconds']<=remaining)
    p.write(b/'optional_control_decision.json',decision)
    if not decision['half_1T_affordable']:p.finalize(b,item,s);return
    label='spatial_half_percent_1T'
    if label not in s['cases']:
        half=initial.copy();half[d.state_out]=delta*.5*unit
        p.integrate_attempt(b,item,s,label,d,half,np.linspace(0,T,cfg['samples_per_period']+1),dynamic_timing(d,cfg,full=True),
            time.perf_counter()+dynamic_remaining(s,cfg)-cfg['cost_reserve_seconds'],plane=False,dense_archive=False)
        if s['hard_stop']:return
    start=time.perf_counter()
    with np.load(b/label/'trajectory.npz') as z:half=z['state'].copy();tt=z['times'].copy()
    with np.load(b/'spatial_full_1percent_3T/trajectory.npz') as z:first=z['state'][:len(tt)].copy()
    with np.load(b/'tangent_growth_comparison.npz') as z:tan=z['tangent'][:len(tt)].copy();base=z['base_state'][:len(tt)].copy()
    e1=analysis.energy_norm(first[:,d.state_out]/delta-tan,metric);e2=analysis.energy_norm(half[:,d.state_out]/(.5*delta)-tan,metric)
    pl1=analysis.energy_norm(first[:,d.state_plane]-base[:,d.state_plane],plmetric);pl2=analysis.energy_norm(half[:,d.state_plane]-base[:,d.state_plane],plmetric)
    one=float(e1.max());two=float(e2.max());floor=cfg['numerical_floor_relative'];resolved=min(one,two)>cfg['observed_error_over_floor_for_order']*floor
    record={'periods':1,'physical_OUT_scales':[cfg['A']*.01,cfg['A']*.005],
        'normalized_OUT_errors':[one,two],'observed_normalized_OUT_order':float(np.log2(one/two)) if one>0 and two>0 else None,
        'PL_backreaction_norms':[float(pl1.max()),float(pl2.max())],
        'observed_PL_backreaction_order':float(np.log2(pl1.max()/pl2.max())) if pl1.max()>0 and pl2.max()>0 else None,
        'order_numerically_resolved':resolved,'order_status':'RESOLVED' if resolved else 'NUMERICAL_FLOOR_LIMITED',
        'does_not_override_growth_confirmation':True,'same_base_direction_phase_and_solver_policy':True}
    np.savez_compressed(b/'half_amplitude_comparison.npz',times=tt,OUT_error_first=e1,OUT_error_half=e2,PL_change_first=pl1,PL_change_half=pl2)
    p.write(b/'half_amplitude_checks.json',record);s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)


def dynamic_main(argv=None):
    parser=argparse.ArgumentParser(description='One fixed beta15 dynamic growth validation using existing ODE helpers')
    mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--through',choices=['initial','planar','tangent','full','all'],default='all');args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=p.cache(b)
        if args.plot_only:dynamic_plot(b);p.finalize(b,p.read(b/'provenance.json'),s)
    else:
        b,item,s=dynamic_prepare()
        if not s['completed'] and not s['hard_stop']:
            start=time.perf_counter();old_seconds=s['numerical_seconds']
            try:
                problem=dynamic_problem(item);s['numerical_seconds']+=time.perf_counter()-start
                dynamic_initial(b,item,s,problem)
                if args.through!='initial':dynamic_planar_check(b,item,s,problem)
                if not s['hard_stop'] and args.through in ('tangent','full','all'):dynamic_tangent_check(b,item,s,problem)
                if not s['hard_stop'] and args.through in ('tangent','full','all'):dynamic_temporal_step_control(b,item,s,problem)
                if not s['hard_stop'] and args.through in ('full','all'):dynamic_full_check(b,item,s,problem)
                if not s['hard_stop'] and args.through=='all':dynamic_optional_controls(b,item,s,problem)
                if s['hard_stop'] and s['overall']=='NUMERICAL_FEASIBILITY_PARTIAL':
                    s['helper_failure_status']=s['overall'];s['overall']='DYNAMIC_NUMERICAL_PARTIAL'
                if args.through in ('full','all') and not s['hard_stop']:
                    dynamic_postprocess(b,item,s,problem);s['completed']=True
            except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as error:
                s.update(completed=True,hard_stop=True,overall='FLOQUET_MODE_RECONSTRUCTION_UNRESOLVED' if s['stages']['C']!='PASS' else 'DYNAMIC_NUMERICAL_PARTIAL',
                    stop_reason=type(error).__name__+': '+str(error))
                s['numerical_seconds']=max(s['numerical_seconds'],old_seconds+time.perf_counter()-start)
            dynamic_plot(b);p.finalize(b,item,s)
    print(json.dumps(p.plain({'bundle':str(b),'summary':s}),ensure_ascii=False));return 0


def dynamic_postprocess(b,item,s,problem):
    """Pointwise safety, physical maps and saved Floquet-wave reconstruction.

    No IVP/eigenproblem; it is never called by completed cache replay.
    """
    if (b/'dynamic_physical_diagnostics.json').exists():return
    from scripts.lib import nlsp_instability_feasibility as f
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    d,model,Phi,U,Omega,lam,coeff,R=problem;cfg=item['config'];start=time.perf_counter();T=2*np.pi/Omega
    maps=[]
    for arm,P in zip(d.arms,d.maps):
        B=arm._raw_basis(np.array([arm.length/2]),0)
        maps.append(np.array([(B@arm._transforms[i]@P[arm.slices[field]])[0] for i,field in enumerate(rod.FIELD_ORDER)]))
    np.savez_compressed(b/'physical_observation_maps.npz',midpoint_maps=np.array(maps),state_plane=d.state_plane,state_out=d.state_out)
    p.write(b/'actual_geometry_and_coefficients.json',{'degree':d.p,'beta_deg':d.beta_deg,'lengths':d.lengths,
        'sections':cfg['sections'],'material':{'E':1.,'rho':1.,'nu':.3,'kappa':5/6},'coefficients':[arm.coefficients.values() for arm in d.arms],
        'local_bases':d.bases,'local_to_global_T7':d.T7,'all_fields':rod.FIELD_ORDER,'ndof':d.ndof,
        'planar_dof_count':len(d.plane),'OUT_dof_count':len(d.out),'physical_units':'normalized model length/time, rotations radians; Lref=1',
        'new_geometric_imperfections':False,'physical_kernels_unchanged':True})
    with np.load(b/'initial_states.npz') as z:
        delta=p.read(b/'Floquet_reconstruction.json')['physical_delta_energy_coordinates'];phase=z['phase_combination'].item();metric=z['metric_OUT'].copy();plmetric=z['metric_PL'].copy();initial=z['base_state'].copy()
    with np.load(b/'spatial_full_1percent_3T/trajectory.npz') as z:times=z['times'].copy();full=z['state'].copy()
    if 'spatial_full_continue_3T_to_5T' in s['cases'] and s['cases']['spatial_full_continue_3T_to_5T']['status']=='PASS':
        with np.load(b/'spatial_full_continue_3T_to_5T/trajectory.npz') as z:times=np.r_[times,z['times'][1:]];full=np.r_[full,z['state'][1:]]
    count=cfg['samples_per_period'];action=f.QuadraticAction(d,model.T4,model.V4);phases=np.arange(count)*2*np.pi/count
    qb,vb,ab=hb_phases(U,Omega,Phi,phases);one=[]
    for ph,q,v,acc in zip(phases,qb,vb,ab):
        canonical=R@(np.exp(1j*np.arange(-9,10)*ph)@coeff)
        S=f.action_operators(action.blocks(q,v,acc))[3]
        one.append(phase*np.exp(lam*ph/Omega)*np.linalg.solve(S,canonical))
    one=np.array(one);Hill=np.array([(np.exp(lam*T*(j//count))*one[j%count]).real for j in range(len(times))])
    hill_extrema=[];velocity_extrema=[];masses=[];joint=[];balances=[];leaks=[]
    for zz in Hill:
        q=np.zeros(d.ndof);q[d.out]=delta*zz[:69];hill_extrema.append(out_polynomial_extrema(d,q,cfg['Lref'])[0])
    for state in full:velocity_extrema.append(out_polynomial_extrema(d,state[d.ndof:],cfg['Lref'])[0])
    indices=np.unique(np.r_[np.linspace(0,len(full)-1,65,dtype=int),np.arange(int(round(times[-1]/T))+1)*count])
    for state in full[indices]:
        masses.append(float(np.linalg.eigvalsh(d.mass_matrix(state[:d.ndof])).min()));jd=d.joint_diagnostics(state[:d.ndof],state[d.ndof:])
        joint.append(max(jd[k] for k in ('translation_difference','orientation_difference','contraction_difference')))
        balances.append(float(np.max(abs(jd['generalized_nodal_balance']))))
        plane=state.copy();plane[d.state_out]=0;leaks.append(float(np.max(abs(d.rhs(0,plane)[d.state_out]))))
    basecase=b/'planar_combined_5T' if times[-1]>3*T+1e-10 else b/'planar_combined_3T'
    with np.load(basecase/'trajectory.npz') as z:base=z['state'].copy()
    backreaction=analysis.energy_norm(full[:,d.state_plane]-base[:,d.state_plane],plmetric)
    fullsettings=p.read(b/'spatial_full_1percent_3T/time_settings.json');atol=np.array(fullsettings['atol'])[d.state_out]
    floor_proxy=float(analysis.energy_norm(atol,metric)/delta)
    diag={'periods':int(round(times[-1]/T)),'maximum_Hill_v_Phi_psi':np.max(hill_extrema,axis=0).tolist(),
        'maximum_full_field_velocity_vt_Phit_psit':np.max(velocity_extrema,axis=0).tolist(),
        'full_mass_min_eigenvalue_sampled':min(masses),'full_joint_compatibility_max_sampled':max(joint),
        'full_generalized_joint_balance_max_sampled':max(balances),'exact_zero_OUT_rhs_max_on_full_PL_components':max(leaks),
        'full_PL_backreaction_energy_norm_max':float(backreaction.max()),
        'full_PL_backreaction_over_reference_PL_norm':float(backreaction.max()/analysis.energy_norm(initial[d.state_plane],plmetric)),
        'full_component_atol_reference_norm_over_delta':floor_proxy,
        'atol_proxy_qualification':'norm of positive component tolerance vector; neither achieved integration error nor rigorous bound',
        'spatial_extrema':'exact roots of represented polynomials; time maxima sampled at T/256, no smoothing or phase fitting',
        'Floquet_reconstruction':'saved N9 coefficients and original raw lambda on HB11; periodic S inverse, no new Hill or ODE',
        'numerical_floor_declared_relative':cfg['numerical_floor_relative'],'scientific_postprocessing_seconds':time.perf_counter()-start}
    np.savez_compressed(b/'Hill_physical_prediction.npz',times=times,OUT_unit_prediction=Hill,physical_extrema=np.array(hill_extrema),
        full_field_velocity_extrema=np.array(velocity_extrema),PL_backreaction_energy_curve=backreaction)
    p.write(b/'dynamic_physical_diagnostics.json',diag);s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)


# One fixed p12 -> p14 spatial check, using the preceding HB/Hill implementation.
P14_CONFIG=p.ROOT/'data/input/nlsp_beta15_out_growth_p14_convergence.json'
P14_OUTPUT=p.ROOT/'results/nlsp_beta15_out_growth_p14_convergence'
P14_SMOKE=p.ROOT/'results/_smoke/nlsp_beta15_out_growth_p14_convergence'


def physical_mass_transfer(source,target):
    """Resting-mass Galerkin projection on common material quadrature.

    This evaluates physical fields, including endpoint lifts and whitening;
    no identification of coefficient indices across discretizations is used.
    """
    cross=np.zeros((target.ndof,source.ndof));gram=np.zeros((target.ndof,target.ndof))
    for a,b,Ps,Pt in zip(source.arms,target.arms,source.maps,target.maps):
        nodes,weights=np.polynomial.legendre.leggauss(source.p+target.p+3)
        points=(nodes+1)*a.length/2;weights=weights*a.length/2;c=b.coefficients
        masses=(c.m,c.m,c.m,c.jp+c.jb,c.jb,c.jp,c.jp)
        Bs=a.basis_at(points);Bt=b.basis_at(points)
        for f,mass in zip(rod.FIELD_ORDER,masses):
            oldB=Bs[f]@Ps[a.slices[f]];newB=Bt[f]@Pt[b.slices[f]]
            cross+=mass*newB.T@(weights[:,None]*oldB)
            gram+=mass*newB.T@(weights[:,None]*newB)
    transfer=np.linalg.solve(gram,cross)
    return transfer,{'target_mass_quadrature_relative_error':float(np.linalg.norm(gram-target.M0)/np.linalg.norm(target.M0)),
        'source_mass_preservation_relative_error':float(np.linalg.norm(transfer.T@target.M0@transfer-source.M0)/np.linalg.norm(source.M0)),
        'source_stiffness_preservation_relative_error':float(np.linalg.norm(transfer.T@target.K@transfer-source.K)/np.linalg.norm(source.K)),
        'method':'physical seven-field resting-mass Galerkin projection on common material quadrature'}


def cross_resolution_MAC(source,target,transfer,old_vectors,new_vectors):
    transported=transfer@old_vectors;overlap=new_vectors.T@target.M0@transported
    oldnorm=np.diag(transported.T@target.M0@transported);newnorm=np.diag(new_vectors.T@target.M0@new_vectors)
    return abs(overlap)**2/(newnorm[:,None]*oldnorm[None,:]),overlap


def p14_prepare():
    cfg=p.read(P14_CONFIG)
    if (cfg['degree']!=14 or cfg['reference_degree']!=12 or cfg['A']!=.005 or cfg['phase_Psi']!=0
        or cfg['beta_deg']!=15 or cfg['mu']!=.1 or cfg['numerical_budget_seconds']>1800
        or cfg['new_ODE_integrations'] or cfg['new_FEM_jobs'] or cfg['automatic_extension']):raise ValueError('Fixed p14 scope required')
    for previous in P14_OUTPUT.glob('*/provenance.json') if P14_OUTPUT.exists() else []:
        item=p.read(previous)
        if item['config']!=cfg:raise ValueError('Different p14 authorization already exists')
        return previous.parent,item,p.cache(previous.parent)
    parent=p.ROOT/cfg['parent'];prior=p.read(parent/'provenance.json')
    if p.sha(parent/'manifest.json')!=cfg['parent_manifest_sha256']:raise ValueError('Historical reference changed')
    references={}
    for name in ('parent','dynamic_parent','two_mode_parent'):
        source=p.ROOT/cfg[name];p.cache(source);references[name]={'path':cfg[name],'sha256':p.sha(source/'manifest.json')}
    if any(p.sha(p.ROOT/name)!=prior['code_hashes'][name] for name in old.f.s.r.KERNELS):raise ValueError('Frozen physical kernel changed')
    files=[P14_CONFIG.relative_to(p.ROOT).as_posix(),Path(__file__).resolve().relative_to(p.ROOT).as_posix()]
    artifacts={**prior['source_artifacts']}
    for folder in (parent,p.ROOT/cfg['dynamic_parent']):
        for name in ('linear_reference.npz','A0.005_HB11.npz','A0.005_OUT_k69_N9.npz','A0.005_OUT_Fourier.npz','initial_states.npz','summary.json'):
            file=folder/name
            if file.exists():artifacts[file.relative_to(p.ROOT).as_posix()]=p.sha(file)
    item={'schema':'nlsp-beta15-p14-spatial-check-v1','config':cfg,'source_manifests':{**prior['source_manifests'],**references},
        'source_artifacts':artifacts,'code_hashes':{**p.code_hashes(),**{name:p.sha(p.ROOT/name) for name in files}},
        'initial_state':p.read(P14_SMOKE/'initial_state.json'),'interpreter':sys.executable,'physical_model_unchanged':True,
        'scope':'one fixed Type-II orbit, p14 HB11/13 and OUT Hill, displacement-only adjoint projection, no ODE/FEM'}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=P14_OUTPUT/fp;b.mkdir(parents=True)
    p.write(b/'config.json',cfg);p.write(b/'provenance.json',item);folder=b/'execution_code';folder.mkdir()
    for name in item['code_hashes']:shutil.copyfile(p.ROOT/name,folder/Path(name).name)
    shutil.copyfile(P14_SMOKE/'preflight.json',b/'cost_forecast.json')
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':p.read(P14_SMOKE/'preflight.json')['scientific_preflight_seconds'],
        'stages':{key:'NOT_RUN' for key in ('A','B','C','D','E','F')},
        'calls':{'linear_eigenproblems':0,'HB_solves':0,'Hill_eigenproblems':0,'ODE_integrations':0,'FEM':0},'deferred':cfg['deferred']}
    p.finalize(b,item,s);return b,item,s


def p14_problem(item):
    cfg=item['config'];parent=p.ROOT/cfg['parent'];prior=p.read(parent/'provenance.json')
    d12,model,Phi12,w12=load_problem(p.read(p.ROOT/prior['config']['parent']/'provenance.json'))
    pc={**p.read(old.s.c.CONFIG),'lengths':cfg['lengths'],'beta_deg':cfg['beta_deg']}
    d=p.make_disc(pc,'C',14,model)
    with np.load(parent/'A0.005_HB11.npz') as z:U12=z['U'].copy()
    with np.load(parent/'linear_reference.npz') as z:out12=np.zeros((d12.ndof,len(d12.out)));out12[d12.out]=z['out_Phi'];ow12=z['out_omega'].copy()
    return d12,d,model,Phi12,w12,U12,out12,ow12


def p14_linear(b,item,s,problem,refresh):
    if (b/'linear_reference.npz').exists():return
    d12,d,model,Phi12,w12,U12,out12,ow12=problem;cfg=item['config'];transfer,diag=physical_mass_transfer(d12,d)
    omega,v=d.eigenpairs(True);ow,ov=d.eigenpairs(False);s['calls']['linear_eigenproblems']+=2
    Phi=np.zeros((d.ndof,len(omega)));Phi[d.plane]=v;out=np.zeros((d.ndof,len(ow)));out[d.out]=ov
    mac,overlap=cross_resolution_MAC(d12,d,transfer,Phi12,Phi);outmac,outoverlap=cross_resolution_MAC(d12,d,transfer,out12,out)
    identity={};gate=True
    for name,oldindex,table,dot,shapes in [('PL2',1,mac,overlap,Phi),('PL4',3,mac,overlap,Phi),('OUT2',1,outmac,outoverlap,out)]:
        column=table[:,oldindex];chosen=int(np.argmax(column));runner=np.partition(column,-2)[-2]
        identity[name]={'p12_sorted_index':oldindex+1,'p14_sorted_index':chosen+1,'mass_MAC':float(column[chosen]),'MAC_margin':float(column[chosen]-runner)}
        gate &= column[chosen]>=cfg['mode_identity_MAC_gate'] and column[chosen]-runner>=cfg['mode_identity_margin_gate']
        if dot[chosen,oldindex]<0:shapes[:,chosen]*=-1
    profile_errors=[]
    for a,c,Pa,Pc in zip(d12.arms,d.arms,d12.maps,d.maps):
        points=np.linspace(0,a.length,103)
        for deriv in (0,1,2):
            ba=a.basis_at(points,deriv);bc=c.basis_at(points,deriv)
            for field in rod.FIELD_ORDER:
                original=ba[field]@Pa[a.slices[field]];lifted=bc[field]@Pc[c.slices[field]]@transfer
                profile_errors.append(np.linalg.norm(lifted-original)/max(np.linalg.norm(original),1e-300))
    diag.update(identity=identity,physical_profile_relative_error_max=float(max(profile_errors)),
        linear_PL_omega2=float(omega[1]),linear_PL_omega4=float(omega[3]),linear_OUT_omega2=float(ow[1]),
        planar_1to3_detuning=float(omega[3]-3*omega[1]),linear_self_OUT2_detuning=float(omega[1]-2*ow[1]),
        PL_mass_orthogonality=float(np.max(abs(Phi.T@d.M0@Phi-np.eye(len(omega))))),
        OUT_mass_orthogonality=float(np.max(abs(out.T@d.M0@out-np.eye(len(ow))))),
        PL_individual_eigen_residuals=[float(np.linalg.norm(d.K@Phi[:,j]-omega[j]**2*d.M0@Phi[:,j])/(np.linalg.norm(d.K@Phi[:,j])+omega[j]**2*np.linalg.norm(d.M0@Phi[:,j]))) for j in (1,3)],
        OUT2_individual_eigen_residual=float(np.linalg.norm(d.K@out[:,1]-ow[1]**2*d.M0@out[:,1])/(np.linalg.norm(d.K@out[:,1])+ow[1]**2*np.linalg.norm(d.M0@out[:,1]))),
        gate_A=bool(gate and identity['PL2']['p14_sorted_index']==2 and identity['PL4']['p14_sorted_index']==4 and identity['OUT2']['p14_sorted_index']==2
            and max(profile_errors)<cfg['lifting_relative_gate'] and diag['source_mass_preservation_relative_error']<cfg['lifting_relative_gate']))
    U=U12@Phi12.T@transfer.T@d.M0@Phi
    np.savez_compressed(b/'linear_reference.npz',Phi=Phi,omega=omega,out_Phi=out[d.out],out_omega=ow,M0=d.M0,K=d.K,
        plane=d.plane,out=d.out,transfer=transfer,PL_MAC=mac,OUT_MAC=outmac,initial_HB_U=U,reference_Phi=Phi12,reference_out_Phi=out12)
    p.write(b/'linear_matching.json',diag);s['stages']['A']='PASS' if diag['gate_A'] else 'FAIL'
    if not diag['gate_A']:s.update(hard_stop=True,overall='P14_BRANCH_IDENTITY_UNRESOLVED')
    refresh();print(json.dumps({'stage':'linear','identity':identity,'omega':[omega[1],omega[3],ow[1]],'lifting_error':max(profile_errors)},ensure_ascii=False),flush=True)


def p14_HB(b,item,s,problem,deadline,refresh):
    if (b/'HB_checks.json').exists():return
    d12,d,model,Phi12,w12,U12,out12,ow12=problem;cfg=item['config']
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy();omega=z['omega'].copy();U=z['initial_HB_U'].copy();transfer=z['transfer'].copy()
    reference=p.read(p.ROOT/cfg['parent']/'summary.json')['cases'][0];Omega=reference['Omega'];rows=[]
    balance=FullPlanarHB(d,Phi,omega,cfg,deadline);initial_active=np.array([U[1,1],U[3,3]])
    for H in cfg['HB_harmonics']:
        name=f'HB{H}'
        if (b/(name+'.json')).exists():
            row=p.read(b/(name+'.json'))
            with np.load(b/(name+'.npz')) as z:U=z['U'].copy()
            Omega=row['Omega']
        else:
            start=time.perf_counter();U,Omega,solver=balance.solve(U,Omega,cfg['A'],H);row,arrays=balance.validate(U,Omega,cfg['A'])
            row.update(Newton=solver,seconds=time.perf_counter()-start);s['calls']['HB_solves']+=1
            np.savez_compressed(b/(name+'.npz'),**arrays);p.write(b/(name+'.json'),row)
        rows.append(row);refresh();print(json.dumps({'stage':name,'Omega':Omega,'residual':row['full_residual_relative'],'seconds':row['seconds']}),flush=True)
    phases=np.arange(65)*2*np.pi/65+.291;oldq,oldv,_=hb_phases(U12,reference['Omega'],Phi12,phases);q,v,_=hb_phases(U,Omega,Phi,phases)
    oldlift=oldq@transfer.T;shape_difference=float(np.max(np.linalg.norm((q-oldlift)@d.M0@Phi,axis=1))/np.max(np.linalg.norm(oldlift@d.M0@Phi,axis=1)))
    branch_change=float(np.max(abs(np.array([U[1,1],U[3,3]])-initial_active)/np.maximum(abs(initial_active),1e-300)))
    gap=abs(rows[-1]['Omega']-rows[0]['Omega'])
    gate=all(r['full_residual_relative']<=cfg['HB_full_relative_gate'] and r['energy_relative_variation']<=cfg['energy_relative_gate']
        and abs(r['actual_A']-cfg['A'])<=cfg['amplitude_absolute_gate'] and r['minimum_mass_eigenvalue']>0
        and r['joint_compatibility']<=cfg['joint_absolute_gate'] and r['last_two_harmonic_energy_norm_ratio']<=cfg['HB_tail_relative_gate']
        and r['Newton']['scaled_projected_residual']<=cfg['HB_projected_relative_target'] for r in rows)
    gate &= gap<=cfg['frequency_control_gate'] and branch_change<=cfg['branch_relative_change_gate'] and shape_difference<=cfg['branch_relative_change_gate'] and U[1,1]>0 and U[3,3]>0
    harmonic_norms=np.linalg.norm(U*omega,axis=1)
    check={'rows':rows,'frequency_harmonic_difference':gap,'branch_active_relative_change':branch_change,
        'physical_waveform_mass_relative_difference':shape_difference,'harmonic_reference_energy_coefficient_norms':harmonic_norms.tolist(),
        'frequency_difference_from_p12':Omega-reference['Omega'],'gate_B':bool(gate),'type':'II','phase_Psi':0.,'phase_or_amplitude_fitting':False}
    p.write(b/'HB_checks.json',check);np.savez_compressed(b/'PL_waveform_comparison.npz',phases=phases,p12_lifted=oldlift,p14=q,p14_velocity=v,Omega12=reference['Omega'],Omega14=Omega)
    s['stages']['B']='PASS' if gate else 'FAIL'
    if not gate:s.update(hard_stop=True,overall='P14_HB_GATE_BLOCKED')
    refresh()


def p14_OUT(b,item,s,problem,deadline,refresh):
    if (b/'OUT_checks.json').exists():return
    d12,d,model,Phi12,w12,U12,out12,ow12=problem;cfg=item['config']
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy();out=z['out_Phi'].copy();ow=z['out_omega'].copy()
    with np.load(b/'HB13.npz') as z:U=z['U'].copy()
    Omega=p.read(b/'HB13.json')['Omega'];records=[]
    if (b/'OUT_Fourier.npz').exists():
        with np.load(b/'OUT_Fourier.npz') as z:F=z['F'].copy();orders=z['orders'].copy();R=z['R'].copy();L=z['L'].copy()
        diag=p.read(b/'OUT_Fourier.json')
    else:
        F,orders,diag,R,L=canonical_fourier(d,model,out,ow,U,Omega,cfg,deadline,d.out,Phi)
        np.savez_compressed(b/'OUT_Fourier.npz',F=F,orders=orders,R=R,L=L);p.write(b/'OUT_Fourier.json',diag)
    if not diag['gate_C']:s.update(hard_stop=True,overall='P14_HILL_SPECTRUM_UNRESOLVED');refresh();return
    for k,N in cfg['OUT_controls']:
        if time.perf_counter()>deadline:raise TimeoutError('p14 OUT budget')
        name=f'OUT_k{k}_N{N}'
        if (b/(name+'.json')).exists():record=p.read(b/(name+'.json'))
        else:
            record,arrays=hill_spectrum(F,orders,Omega,N,k,cfg,retain_reciprocal=k==len(ow))
            np.savez_compressed(b/(name+'.npz'),**arrays);p.write(b/(name+'.json'),record);s['calls']['Hill_eigenproblems']+=1
        records.append(record);refresh();print(json.dumps({'stage':'OUT','k':k,'N':N,'sigma':record['sigma'],'seconds':record['seconds']}),flush=True)
    # Orbit-refinement sensitivity on the same restricted space, with all
    # state-dependent mass/inertia terms, not a stiffness-only surrogate.
    with np.load(b/'HB11.npz') as z:Uc=z['U'].copy()
    Oc=p.read(b/'HB11.json')['Omega']
    holder=__import__('copy').copy(d);holder.out=d.out
    action=old.f.QuadraticAction(holder,model.T4,model.V4);phases=np.arange(cfg['tangent_phase_samples'])*2*np.pi/cfg['tangent_phase_samples']
    qc,vc,ac=hb_phases(Uc,Oc,Phi,phases);ids=np.r_[np.arange(24),np.arange(24)+len(ow)]
    values0=np.array([L[ids]@old.f.action_operators(action.blocks(q,v,a))[1]@R[:,ids] for q,v,a in zip(qc,vc,ac)])
    Fc=np.einsum('mt,tij->mij',np.exp(-1j*orders[:,None]*phases)/len(phases),values0,optimize=True)
    comparator,arrays=hill_spectrum(Fc,orders,Oc,7,24,cfg);s['calls']['Hill_eigenproblems']+=1
    p.write(b/'OUT_HB11_control.json',comparator);np.savez_compressed(b/'OUT_HB11_control.npz',**arrays)
    phases=np.arange(12)*2*np.pi/12+.271;q,v,a=hb_phases(U,Omega,Phi,phases)
    fourier_error=max(np.linalg.norm(L@old.f.action_operators(action.blocks(qi,vi,ai))[1]@R-np.einsum('m,mij->ij',np.exp(1j*orders*ph),F))/np.linalg.norm(L@old.f.action_operators(action.blocks(qi,vi,ai))[1]@R) for ph,qi,vi,ai in zip(phases,q,v,a))
    modal=abs(records[-1]['sigma']-records[1]['sigma']);temporal=abs(records[-1]['sigma']-records[-2]['sigma']);orbit=abs(records[0]['sigma']-comparator['sigma'])
    uncertainty=max(modal,temporal,orbit,cfg['spectral_numerical_floor']);last=records[-1];lam=complex(last['dominant_raw_lambda']['real'],last['dominant_raw_lambda']['imag']);rho=np.exp(lam*2*np.pi/Omega)
    gate=bool(max(modal,temporal,orbit)<=cfg['spectral_control_absolute_gate'] and fourier_error<=cfg['operator_relative_gate']
        and max(last['independent_eigen_residual_local_relative'],last['reciprocal_independent_local_residual'])<=cfg['Hill_residual_gate']
        and last['quartet_absolute_defect']<=cfg['spectral_control_absolute_gate'] and last['reciprocal_pair_absolute_defect']<=cfg['spectral_control_absolute_gate'])
    check={'records':records,'modal_change_48_full':modal,'temporal_change_N7_N9':temporal,'orbit_HB11_HB13_change_k24':orbit,
        'uncertainty_indicator':uncertainty,'independent_Fourier_relative_error':float(fourier_error),'Omega':Omega,'sigma':last['sigma'],
        'sigma_T':float(last['sigma']*2*np.pi/Omega),'multiplier':old.s.complex_row(rho),'reciprocal_multiplier':old.s.complex_row(np.exp(complex(last['reciprocal_raw_lambda']['real'],last['reciprocal_raw_lambda']['imag'])*2*np.pi/Omega)),
        'all_OUT_modes':len(ow),'gate_C':gate,'resolved_positive_growth':bool(gate and last['sigma']>cfg['growth_uncertainty_factor']*uncertainty),
        'negative_real_multiplier':bool(rho.real<0 and abs(rho.imag)/abs(rho)<cfg['real_multiplier_relative_gate'])}
    p.write(b/'OUT_checks.json',check);s['stages']['C']='PASS' if gate else 'FAIL'
    if not gate:s.update(hard_stop=True,overall='P14_HILL_SPECTRUM_UNRESOLVED')
    elif not check['resolved_positive_growth']:s.update(hard_stop=True,overall='P14_OUT_GROWTH_NOT_REPRODUCED')
    refresh()
    # Reproducible handling of the documented p14 resolution trigger. This
    # adds exactly one near-full modal control, retaining the first failed
    # qualification and ALL previous solves, with no threshold change.
    if not gate and modal>cfg['spectral_control_absolute_gate']:
        policy=b/'supplementary_modal_policy.json'
        if not policy.exists():
            p.write(policy,{'trigger':'48/full81 modal difference exceeds unchanged gate',
                'one_additional_control':[72,7],'reference_full_controls_already_computed':[[81,7],[81,9]],
                'spectral_control_absolute_gate_unchanged':cfg['spectral_control_absolute_gate'],
                'forecast_seconds':20,'budget_seconds':cfg['numerical_budget_seconds'],
                'no_repeated_HB_or_Hill_prefix':'one new 72-mode solve; original qualification and costs retained'})
            refresh()
        p14_modal_control(b,item,s,deadline,refresh)


def real_Floquet_direction(canonical,S,metric,reference=None):
    physical=np.linalg.solve(S,canonical);B=np.c_[physical.real,physical.imag];Gram=B.T@metric@B
    values,vectors=np.linalg.eigh(Gram)
    if values[-1]<=0:raise ArithmeticError('Zero reconstructed physical Floquet mode')
    phase=(vectors[0,-1]-1j*vectors[1,-1])/np.sqrt(values[-1]);unit=(phase*physical).real
    if reference is not None and unit@metric@reference<0:phase=-phase;unit=-unit
    return unit,phase,float(max(0.,values[0])/values[-1])


def p14_modal_control(b,item,s,deadline,refresh):
    """One triggered 72/full81 check; preserve the insufficient 48-mode result."""
    cfg=item['config'];check=p.read(b/'OUT_checks.json');policy=p.read(b/'supplementary_modal_policy.json')
    if policy['one_additional_control']!=[72,7] or check.get('last_modal_change_72_full') is not None:return
    last=check['records'][-1]
    other_gates=(max(check['temporal_change_N7_N9'],check['orbit_HB11_HB13_change_k24'])<=cfg['spectral_control_absolute_gate']
        and check['independent_Fourier_relative_error']<=cfg['operator_relative_gate']
        and max(last['independent_eigen_residual_local_relative'],last['reciprocal_independent_local_residual'])<=cfg['Hill_residual_gate']
        and last['quartet_absolute_defect']<=cfg['spectral_control_absolute_gate']
        and last['reciprocal_pair_absolute_defect']<=cfg['spectral_control_absolute_gate'])
    if not other_gates:return
    if deadline-time.perf_counter()<policy['forecast_seconds']:return
    for original,archive in [('OUT_checks.json','OUT_initial_48_full_qualification.json'),('summary.json','initial_modal_qualification_summary.json')]:
        if not (b/archive).exists():shutil.copyfile(b/original,b/archive)
    with np.load(b/'OUT_Fourier.npz') as z:F=z['F'].copy();orders=z['orders'].copy()
    record,arrays=hill_spectrum(F,orders,check['Omega'],7,72,cfg);s['calls']['Hill_eigenproblems']+=1
    p.write(b/'OUT_k72_N7.json',record);np.savez_compressed(b/'OUT_k72_N7.npz',**arrays)
    change=abs(record['sigma']-check['sigma']);check['additional_modal_control']=record;check['last_modal_change_72_full']=change
    uncertainty=max(change,check['temporal_change_N7_N9'],check['orbit_HB11_HB13_change_k24'],cfg['spectral_numerical_floor'])
    check['initial_gate_C_preserved_in']='OUT_initial_48_full_qualification.json';check['uncertainty_indicator']=uncertainty
    check['gate_C']=bool(other_gates and change<=cfg['spectral_control_absolute_gate']
        and record['independent_eigen_residual_local_relative']<=cfg['Hill_residual_gate'])
    check['resolved_positive_growth']=bool(check['gate_C'] and check['sigma']>cfg['growth_uncertainty_factor']*uncertainty)
    p.write(b/'OUT_checks.json',check);s['stages']['C']='PASS' if check['gate_C'] else 'FAIL'
    if check['resolved_positive_growth']:s.update(hard_stop=False,overall='PREPARED')
    refresh();print(json.dumps({'stage':'supplementary_modal_control','k':72,'sigma':record['sigma'],'change_full81':change,'gate':check['gate_C'],'seconds':record['seconds']}),flush=True)


def p14_projection(b,item,s,problem,deadline,refresh):
    """Symplectic dual from the reciprocal Floquet solution, no monodromy IVP."""
    if (b/'projection_checks.json').exists():return
    d12,d,model,Phi12,w12,U12,out12,ow12=problem;cfg=item['config']
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy();ow=z['out_omega'].copy();out=z['out_Phi'].copy();transfer=z['transfer'].copy()
    with np.load(b/'HB13.npz') as z:U=z['U'].copy()
    with np.load(b/'OUT_Fourier.npz') as z:R=z['R'].copy()
    Omega=p.read(b/'HB13.json')['Omega'];action=old.f.QuadraticAction(d,model.T4,model.V4)
    q0,v0,a0=hb_phases(U,Omega,Phi,np.array([0.]));_,Ac,_,S,_,J=old.f.action_operators(action.blocks(q0[0],v0[0],a0[0]))
    E=dynamic_metric(d,d.out);T=2*np.pi/Omega;records=[];states=[]
    with np.load(p.ROOT/cfg['dynamic_parent']/'initial_states.npz') as z:reference=z['OUT_unit_state'].copy()
    lift_out=transfer[np.ix_(d.out,d12.out)];ref=np.r_[lift_out@reference[:len(d12.out)],lift_out@reference[len(d12.out):]]
    disp=np.r_[out[:,1],np.zeros(len(d.out))];position=np.zeros(d.ndof);position[d.out]=disp[:len(d.out)]
    scale=.01*cfg['A']/out_polynomial_extrema(d,position,cfg['Lref'])[1];disp*=scale;disp_c=S@disp
    for N in (7,9):
        row=p.read(b/f'OUT_k81_N{N}.json')
        with np.load(b/f'OUT_k81_N{N}.npz') as z:cu=z['dominant_eigenvector'].reshape(2*N+1,162).copy();cs=z['reciprocal_eigenvector'].reshape(2*N+1,162).copy()
        lu=complex(row['dominant_raw_lambda']['real'],row['dominant_raw_lambda']['imag']);ls=complex(row['reciprocal_raw_lambda']['real'],row['reciprocal_raw_lambda']['imag'])
        ru,pu,ranku=real_Floquet_direction(R@cu.sum(axis=0),S,E,ref);rs,ps,ranks=real_Floquet_direction(R@cs.sum(axis=0),S,E)
        uc=S@ru;sc=S@rs;den=float(sc@J@uc)
        if abs(den)<np.finfo(float).eps*np.linalg.norm(sc)*np.linalg.norm(uc):raise ArithmeticError('Reciprocal pair cannot define a resolved symplectic dual')
        left=sc@J/den;left_phys=left@S;alpha=float(left@disp_c);dualnorm=float(np.sqrt(left_phys@np.linalg.solve(E,left_phys)))
        phases=np.arange(33)*2*np.pi/33+.137;q,v,a=hb_phases(U,Omega,Phi,phases);resu=[];ress=[];wronskians=[]
        for ph,qi,vi,ai in zip(phases,q,v,a):
            if time.perf_counter()>deadline:raise TimeoutError('Floquet projection diagnostics budget')
            _,operator,_,_,_,_=old.f.action_operators(action.blocks(qi,vi,ai))
            n=np.arange(-N,N+1);eu=np.exp((lu/Omega+1j*n)*ph);es=np.exp((ls/Omega+1j*n)*ph)
            u=(pu*R@(eu@cu)).real;stable=(ps*R@(es@cs)).real
            udot=(pu*R@((eu*(lu+1j*n*Omega))@cu)).real;sdot=(ps*R@((es*(ls+1j*n*Omega))@cs)).real
            resu.append(np.linalg.norm(udot-operator@u)/max(np.linalg.norm(udot)+np.linalg.norm(operator@u),1e-300))
            ress.append(np.linalg.norm(sdot-operator@stable)/max(np.linalg.norm(sdot)+np.linalg.norm(operator@stable),1e-300))
            wronskians.append(stable@J@u/den)
        record={'N':N,'unstable_real_Gram_ratio':ranku,'stable_real_Gram_ratio':ranks,'symplectic_denominator':den,
            'left_right_duality_error':abs(float(left@uc)-1.),'left_stable_zero_error':abs(float(left@sc)),
            'physical_dual_E0_norm':dualnorm,'unstable_coefficient':alpha,'initial_OUT_E0_norm':float(np.sqrt(disp@E@disp)),
            'coefficient_over_initial_E0_norm':alpha/float(np.sqrt(disp@E@disp)),
            'unstable_pointwise_residual':float(max(resu)),'stable_pointwise_residual':float(max(ress)),
            'symplectic_Wronskian_relative_error':float(max(abs(np.array(wronskians)-1))),
            'p12_p14_unstable_state_E0_MAC':float((ru@E@ref)**2/((ru@E@ru)*(ref@E@ref)))}
        records.append(record);states.append((ru,rs,left,left_phys,pu,ps));refresh()
    change=abs(records[-1]['unstable_coefficient']-records[0]['unstable_coefficient'])/max(abs(records[-1]['unstable_coefficient']),1e-300)
    last=records[-1];valid=all(max(r['unstable_real_Gram_ratio'],r['stable_real_Gram_ratio'])<=cfg['real_Gram_ratio_gate']
        and max(r['unstable_pointwise_residual'],r['stable_pointwise_residual'],r['symplectic_Wronskian_relative_error'])<=cfg['Floquet_pointwise_relative_gate']
        and max(r['left_right_duality_error'],r['left_stable_zero_error'])<=cfg['operator_relative_gate'] for r in records)
    valid &= change<=cfg['projection_relative_control_gate']
    uncertainty=max(abs(records[-1]['unstable_coefficient']-records[0]['unstable_coefficient']),
        cfg['spectral_numerical_floor']*last['physical_dual_E0_norm']*last['initial_OUT_E0_norm'])
    nonzero=valid and abs(last['unstable_coefficient'])>cfg['growth_uncertainty_factor']*uncertainty
    ru,rs,left,left_phys,pu,ps=states[-1];fullru=np.zeros(2*d.ndof);fullru[d.state_out]=ru;fulldisp=np.zeros(2*d.ndof);fulldisp[d.state_out]=disp
    parts=d.reference_energy_parts(fullru[:d.ndof],fullru[d.ndof:]);profiles=np.array([arm.reconstruct(P@fullru[:d.ndof],np.linspace(0,arm.length,201)) for arm,P in zip(d.arms,d.maps)])
    oldfull=np.zeros(2*d12.ndof);oldfull[d12.state_out]=reference
    oldprofiles=np.array([arm.reconstruct(P@oldfull[:d12.ndof],np.linspace(0,arm.length,201)) for arm,P in zip(d12.arms,d12.maps)])
    diag={'records':records,'projection_relative_N7_N9_change':change,'projection_uncertainty_indicator':uncertainty,
        'projection_valid':bool(valid),'displacement_only_excites_growing_subspace':bool(nonzero),
        'projection_method':'left_u = canonical_stable.T J / (canonical_stable.T J canonical_unstable); symplectic reciprocal pair, not Euclidean right-vector overlap',
        'projection_is_not_an_energy_fraction':True,'short_tangent_control':'NOT_RUN_ANALYTIC_DUAL_SUFFICIENT',
        'initial_physical_OUT_scale':.01*cfg['A'],'all_initial_OUT_velocities_zero':True,
        'dangerous_OUT_bending_fraction':float(parts[:,0].sum()/parts.sum()),'dangerous_OUT_torsion_fraction':float(parts[:,1].sum()/parts.sum()),
        'unstable_state_MAC':last['p12_p14_unstable_state_E0_MAC'],'gate_D':bool(last['p12_p14_unstable_state_E0_MAC']>=cfg['mode_identity_MAC_gate'])}
    p.write(b/'projection_checks.json',diag);np.savez_compressed(b/'physical_Floquet_and_displacement.npz',OUT_unit=ru,stable_OUT_unit=rs,
        canonical_left=left,physical_left=left_phys,displacement_only_OUT=disp,displacement_only_full_state=fulldisp,
        metric=E,S=S,J=J,p12_lifted_unit=ref,p14_profiles=profiles,p12_profiles=oldprofiles,reference_energy_parts=parts,
        unstable_phase=pu,stable_phase=ps)
    s['stages']['D']='PASS' if diag['gate_D'] else 'FAIL';s['stages']['E']='PASS' if valid else 'PARTIAL'
    if not diag['gate_D']:s.update(hard_stop=True,overall='P14_BRANCH_IDENTITY_UNRESOLVED')
    refresh();print(json.dumps({'stage':'physical_matching_projection','MAC':diag['unstable_state_MAC'],'coefficient':last['unstable_coefficient'],'valid':bool(valid),'excites':bool(nonzero)}),flush=True)


def p14_PL(b,item,s,problem,deadline,refresh):
    if (b/'PL_checks.json').exists():return
    d12,d,model,Phi12,w12,U12,out12,ow12=problem;cfg=item['config']
    if deadline-time.perf_counter()<120+cfg['cost_reserve_seconds']:
        s['stages']['F']='NOT_RUN_RESOURCE_GATE';refresh();return
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy();omega=z['omega'].copy()
    with np.load(b/'HB13.npz') as z:U=z['U'].copy()
    Omega=p.read(b/'HB13.json')['Omega']
    F,orders,diag,R,L=canonical_fourier(d,model,Phi[d.plane],omega,U,Omega,cfg,deadline,d.plane,Phi)
    np.savez_compressed(b/'PL_Fourier.npz',F=F,orders=orders,R=R,L=L);p.write(b/'PL_Fourier.json',diag)
    records=[]
    for k,N in cfg['PL_controls']:
        if time.perf_counter()>deadline:raise TimeoutError('p14 secondary PL control budget')
        record,arrays=hill_spectrum(F,orders,Omega,N,k,cfg);records.append(record);s['calls']['Hill_eigenproblems']+=1
        np.savez_compressed(b/f'PL_k{k}_N{N}.npz',**arrays);p.write(b/f'PL_k{k}_N{N}.json',record)
        ev=arrays['folded_eigenvalues'][arrays['accepted']];band=cfg['spectral_control_absolute_gate'];nonzero=ev[abs(ev)>band]
        record['non_neutral_sigma']=float(max(0.,max(nonzero.real,default=0.)))
        record['near_zero_eigenvalues']=[old.s.complex_row(value) for value in sorted(ev,key=abs)[:8]]
        refresh();print(json.dumps({'stage':'PL','k':k,'N':N,'non_neutral_sigma':record['non_neutral_sigma'],'seconds':record['seconds']}),flush=True)
    holder=__import__('copy').copy(d);holder.out=d.plane;action=old.f.QuadraticAction(holder,model.T4,model.V4)
    phases=np.arange(cfg['tangent_phase_samples'])*2*np.pi/cfg['tangent_phase_samples'];q,v,a=hb_phases(U,Omega,Phi,phases)
    neutral=np.array([L@old.f.action_operators(action.blocks(qi,vi,ai))[3]@np.r_[vi[d.plane],ai[d.plane]] for qi,vi,ai in zip(q,v,a)])
    n=np.arange(-13,14);coef=np.einsum('nt,ti->ni',np.exp(-1j*n[:,None]*phases)/len(phases),neutral,optimize=True)
    residual=hill_action(F,orders,Omega,13,coef.ravel());neutral_res=float(np.linalg.norm(residual)/(Omega*np.linalg.norm(coef)))
    p.write(b/'PL_checks.json',{'records':records,'canonical':diag,'time_shift_neutral_relative_residual':neutral_res,
        'non_neutral_sigma':records[-1]['non_neutral_sigma'],'neutral_band':cfg['spectral_control_absolute_gate'],
        'gate_F':bool(diag['gate_C'] and neutral_res<=cfg['operator_relative_gate'] and max(r['independent_eigen_residual_local_relative'] for r in records)<=cfg['Hill_residual_gate']),
        'scope':'all 108 p14 PL modes, interior replicas N5/7; neutral band separated, not a global stability proof'})
    np.savez_compressed(b/'PL_time_shift_neutral.npz',coefficients=coef,orders=n,residual=residual)
    s['stages']['F']='PASS' if p.read(b/'PL_checks.json')['gate_F'] else 'PARTIAL';refresh()


def p14_comparison(b,item,s):
    cfg=item['config'];source=p.read(p.ROOT/cfg['parent']/'summary.json')['cases'][0]
    hb=p.read(b/'HB13.json');out=p.read(b/'OUT_checks.json');proj=p.read(b/'projection_checks.json')
    T12=2*np.pi/source['Omega'];T14=2*np.pi/hb['Omega'];sigma12=source['OUT_sigma'];sigma14=out['sigma']
    comparison={'p12':{'Omega':source['Omega'],'sigma':sigma12,'sigma_T':sigma12*T12,
        'multiplier_real':-float(np.exp(sigma12*T12)),'active_PL4_fraction':source['initial_mode4_quadratic_active_energy_fraction'],
        'HB_residual':source['full_residual_relative'],'e_folding_periods':1/(sigma12*T12),
        'OUT_bending_fraction':source['dangerous_OUT_reference_bending_fraction'],'OUT_torsion_fraction':source['dangerous_OUT_reference_torsion_fraction']},
        'p14':{'Omega':hb['Omega'],'sigma':sigma14,'sigma_T':sigma14*T14,'multiplier_real':out['multiplier']['real'],
        'active_PL4_fraction':hb['initial_mode4_quadratic_active_energy_fraction'],'HB_residual':hb['full_residual_relative'],
        'e_folding_periods':1/(sigma14*T14),'OUT_bending_fraction':proj['dangerous_OUT_bending_fraction'],'OUT_torsion_fraction':proj['dangerous_OUT_torsion_fraction']},
        'sigma_relative_difference':abs(sigma14-sigma12)/sigma12,'sigma_T_relative_difference':abs(sigma14*T14-sigma12*T12)/(sigma12*T12),
        'unstable_OUT_state_MAC':proj['unstable_state_MAC'],'p12_results_reused_not_recomputed':True,
        'two_p_values_are_not_a_proof_of_p_convergence':True}
    strong=out['resolved_positive_growth'] and out['negative_real_multiplier'] and proj['gate_D'] and comparison['sigma_relative_difference']<=cfg['strong_p_agreement_relative_gate']
    s['overall']='P14_OUT_GROWTH_REPRODUCED' if strong else 'P14_OUT_GROWTH_REPRODUCED_CONVERGENCE_PARTIAL'
    s['completed']=True;p.write(b/'comparison.json',comparison)


def p14_physical_diagnostics(b,item,s,problem):
    if (b/'physical_observations.json').exists():return
    d12,d,model,Phi12,w12,U12,out12,ow12=problem;start=time.perf_counter();cfg=item['config']
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy();omega=z['omega'].copy()
    with np.load(b/'HB13.npz') as z:U=z['U'].copy()
    observations={};profiles=[];midpoints=[]
    for name,disc,shapes,frequencies,coeff in [('p12',d12,Phi12,w12,U12),('p14',d,Phi,omega,U)]:
        harmonic_q=coeff@shapes.T;init=harmonic_q.sum(axis=0);energy_by_mode=.5*frequencies**2*coeff.sum(axis=0)**2
        norms=np.linalg.norm(coeff*frequencies,axis=1)
        observations[name]={'harmonic_E0_position_norms':norms.tolist(),'harmonic_norm_ratios_to_fundamental':(norms/norms[1]).tolist(),
            'harmonic_exact_max_w':[p.polynomial_w_max(disc,row)[0] for row in harmonic_q],
            'initial_quadratic_modal_energy':energy_by_mode.tolist(),'initial_quadratic_total_energy':float(energy_by_mode.sum()),
            'initial_PL2_fraction_of_total_reference_energy':float(energy_by_mode[1]/energy_by_mode.sum()),
            'initial_PL4_fraction_of_total_reference_energy':float(energy_by_mode[3]/energy_by_mode.sum()),
            'initial_other_modes_fraction_of_total_reference_energy':float(1-(energy_by_mode[1]+energy_by_mode[3])/energy_by_mode.sum())}
        profiles.append(np.array([arm.reconstruct(P@init,np.linspace(0,arm.length,201)) for arm,P in zip(disc.arms,disc.maps)]))
    with np.load(b/'PL_waveform_comparison.npz') as z:oldq=z['p12_lifted'].copy();q=z['p14'].copy();phases=z['phases'].copy()
    midpoint_map=np.array([(arm.basis_at(np.array([arm.length/2]))['w']@P[arm.slices['w']])[0] for arm,P in zip(d.arms,d.maps)])
    projection=p.read(b/'projection_checks.json');last=projection['records'][-1]
    # Add a conservative pointwise-residual indicator to the N7/N9/floor
    # indicators. This is not a rigorous coefficient-error bound.
    projection['initial_N7_N9_floor_uncertainty_indicator']=projection['projection_uncertainty_indicator']
    extra=last['initial_OUT_E0_norm']*last['physical_dual_E0_norm']*max(last['unstable_pointwise_residual'],last['stable_pointwise_residual'],last['symplectic_Wronskian_relative_error'])
    projection['pointwise_residual_projection_uncertainty_indicator']=extra
    projection['projection_uncertainty_indicator']=max(extra,projection['projection_uncertainty_indicator'])
    projection['displacement_only_excites_growing_subspace']=bool(projection['projection_valid'] and abs(last['unstable_coefficient'])>cfg['growth_uncertainty_factor']*projection['projection_uncertainty_indicator'])
    projection['uncertainty_qualification']='observed Fourier and pointwise residual indicators, not a rigorous adjoint error bound'
    with np.load(b/'physical_Floquet_and_displacement.npz') as z:
        ru=z['OUT_unit'].copy();disp=z['displacement_only_full_state'].copy()
    q0,v0,_=hb_phases(U,p.read(b/'HB13.json')['Omega'],Phi,np.array([0.]));base=np.r_[q0[0],v0[0]];initial=base+disp
    phys=out_polynomial_extrema(d,initial[:d.ndof],cfg['Lref'])[0]
    jj=d.joint_diagnostics(initial[:d.ndof],initial[d.ndof:]);unit=np.zeros(d.ndof);unit[d.out]=ru[:len(d.out)]
    unit_scale=.01*cfg['A']/out_polynomial_extrema(d,unit,cfg['Lref'])[1]
    source_delta=p.read(p.ROOT/cfg['dynamic_parent']/'Floquet_reconstruction.json')['physical_delta_energy_coordinates']
    observations.update(displacement_only_initial_OUT_max_v_Phi_psi=phys.tolist(),
        displacement_only_initial_velocity_max=float(np.max(abs(initial[d.state_out[len(d.out):]]))),
        displacement_only_joint_compatibility=max(jj[key] for key in ('translation_difference','orientation_difference','contraction_difference')),
        reference_initial_Floquet_physical_scale=source_delta,p14_initial_Floquet_physical_scale=unit_scale,
        angular_units='radians',length_time_units='normalized model',no_new_ODE=True)
    np.savez_compressed(b/'physical_observations.npz',p12_initial_profiles=profiles[0],p14_initial_profiles=profiles[1],phases=phases,
        p12_midpoint_w=oldq@midpoint_map.T,p14_midpoint_w=q@midpoint_map.T,
        displacement_only_initial_base_state=base,displacement_only_initial_full_state=initial)
    p.write(b/'projection_checks.json',projection);p.write(b/'physical_observations.json',observations)
    s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)


def p14_mode_resolution_check(b,item,s):
    """Physical-state MAC for saved temporal and final modal comparisons."""
    if (b/'OUT_mode_resolution.json').exists():return
    start=time.perf_counter()
    with np.load(b/'OUT_Fourier.npz') as z:R=z['R'].copy()
    with np.load(b/'physical_Floquet_and_displacement.npz') as z:S=z['S'].copy();E=z['metric'].copy();reference=z['OUT_unit'].copy()
    checks={}
    for k,N in ((81,7),(72,7)):
        with np.load(b/f'OUT_k{k}_N{N}.npz') as z:coeff=z['dominant_eigenvector'].reshape(2*N+1,2*k).sum(axis=0)
        ids=np.r_[np.arange(k),np.arange(k)+81];state,_,_=real_Floquet_direction(R[:,ids]@coeff,S,E,reference)
        checks[f'k{k}_N{N}_vs_full81_N9_MAC']=float((state@E@reference)**2/((state@E@state)*(reference@E@reference)))
    checks['gate']=bool(min(checks.values())>=item['config']['mode_identity_MAC_gate'])
    p.write(b/'OUT_mode_resolution.json',checks)
    cmp=p.read(b/'comparison.json');cmp['uncertainty_indicator_p12']=p.read(p.ROOT/item['config']['parent']/'summary.json')['cases'][0]['OUT_uncertainty_indicator']
    cmp['uncertainty_indicator_p14']=p.read(b/'OUT_checks.json')['uncertainty_indicator'];p.write(b/'comparison.json',cmp)
    if not checks['gate']:s['overall']='P14_OUT_GROWTH_REPRODUCED_CONVERGENCE_PARTIAL'
    s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)


def p14_plot(b):
    if not (b/'comparison.json').exists():return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with np.load(b/'PL_waveform_comparison.npz') as z:oldq=z['p12_lifted'].copy();q=z['p14'].copy();phases=z['phases'].copy()
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy();M=z['M0'].copy()
    oldm=oldq@M@Phi;newm=q@M@Phi;cmp=p.read(b/'comparison.json')
    fig,ax=plt.subplots(1,3,figsize=(13,3.7))
    if (b/'physical_observations.npz').exists():
        with np.load(b/'physical_observations.npz') as z:wm12=z['p12_midpoint_w'].copy();wm14=z['p14_midpoint_w'].copy()
        for arm in range(2):
            ax[0].plot(phases/(2*np.pi),wm12[:,arm],label=f'p12 arm {arm+1}')
            ax[0].plot(phases/(2*np.pi),wm14[:,arm],'--',label=f'p14 arm {arm+1}')
    ax[0].set(xlabel='phase / (2 pi)',ylabel='midpoint w [model length]');ax[0].legend(fontsize=8)
    for j,label in ((1,'PL2'),(3,'PL4')):
        ax[1].plot(phases/(2*np.pi),oldm[:,j],label='p12 lifted '+label)
        ax[1].plot(phases/(2*np.pi),newm[:,j],'--',label='p14 '+label)
    ax[1].set(xlabel='phase / (2 pi)',ylabel='mass-normalized modal displacement');ax[1].legend(fontsize=8)
    x=np.arange(2);shares=[cmp['p12']['active_PL4_fraction'],cmp['p14']['active_PL4_fraction']]
    ax[2].bar(x,100*np.array(shares));ax[2].set(xticks=x,xticklabels=['p12','p14'],ylabel='PL4 active reference energy share [%]')
    fig.tight_layout();fig.savefig(b/'p12_p14_planar_branch.png',dpi=170);fig.savefig(b/'p12_p14_planar_branch.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    with np.load(b/'physical_Floquet_and_displacement.npz') as z:p12=z['p12_profiles'];p14=z['p14_profiles']
    if (b/'physical_observations.json').exists():
        obs=p.read(b/'physical_observations.json');p12=p12*obs['reference_initial_Floquet_physical_scale'];p14=p14*obs['p14_initial_Floquet_physical_scale']
    fig,ax=plt.subplots(1,2,figsize=(10,3.7))
    sigma=[cmp['p12']['sigma'],cmp['p14']['sigma']]
    errors=[cmp.get('uncertainty_indicator_p12',0),cmp.get('uncertainty_indicator_p14',0)]
    ax[0].errorbar([12,14],sigma,yerr=errors,fmt='o',capsize=4)
    ax[0].set_ylim(0,max(sigma)*1.25);ax[0].ticklabel_format(axis='y',style='plain',useOffset=False)
    ax[0].text(.04,.92,f'observed difference: {100*cmp["sigma_relative_difference"]:.5f}%',transform=ax[0].transAxes)
    ax[0].set(xticks=[12,14],xlabel='Galerkin degree p',ylabel='OUT growth exponent [1/model time]')
    for arm in range(2):
        x=np.linspace(0,1,201);ax[1].plot(x,p12[arm,:,2],label=f'p12 arm {arm+1}');ax[1].plot(x,p14[arm,:,2],'--',label=f'p14 arm {arm+1}')
    ax[1].set(xlabel='material s/L',ylabel='v [model length], OUT scale 1% A');ax[1].legend(fontsize=8)
    fig.tight_layout();fig.savefig(b/'p12_p14_OUT_growth_mode.png',dpi=170);fig.savefig(b/'p12_p14_OUT_growth_mode.pdf',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)


def p14_main(argv=None):
    parser=argparse.ArgumentParser(description='One fixed p14 HB/Hill spatial convergence check, no ODE/FEM')
    modes=parser.add_mutually_exclusive_group(required=True);modes.add_argument('--compute',action='store_true');modes.add_argument('--report-only',type=Path);modes.add_argument('--plot-only',type=Path)
    parser.add_argument('--through',choices=['linear','hb','out','all'],default='all');args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=p.cache(b)
        if args.plot_only:p14_plot(b);p.finalize(b,p.read(b/'provenance.json'),s)
    else:
        b,item,s=p14_prepare()
        triggered_modal=(not s['completed'] and s['hard_stop'] and s['overall']=='P14_HILL_SPECTRUM_UNRESOLVED'
            and (b/'supplementary_modal_policy.json').exists())
        if not s['completed'] and (not s['hard_stop'] or triggered_modal):
            start=time.perf_counter();offset=s['numerical_seconds'];deadline=start+item['config']['numerical_budget_seconds']-offset-item['config']['cost_reserve_seconds']
            def refresh():s['numerical_seconds']=offset+time.perf_counter()-start;p.finalize(b,item,s)
            try:
                if triggered_modal:p14_modal_control(b,item,s,deadline,refresh)
                problem=p14_problem(item);p14_linear(b,item,s,problem,refresh)
                if not s['hard_stop'] and args.through!='linear':p14_HB(b,item,s,problem,deadline,refresh)
                if not s['hard_stop'] and args.through in ('out','all'):p14_OUT(b,item,s,problem,deadline,refresh)
                if not s['hard_stop'] and args.through=='all':
                    p14_projection(b,item,s,problem,deadline,refresh)
                    if not s['hard_stop']:
                        p14_PL(b,item,s,problem,deadline,refresh);p14_comparison(b,item,s);refresh();p14_physical_diagnostics(b,item,s,problem);p14_mode_resolution_check(b,item,s)
            except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as error:
                s.update(hard_stop=True,completed=True,overall='NUMERICAL_PARTIAL',stop_reason=type(error).__name__+': '+str(error));refresh()
            p14_plot(b);p.finalize(b,item,s)
    print(json.dumps({'bundle':str(b),'status':s['overall'],'stages':s['stages'],'seconds':s['numerical_seconds'],'calls':s['calls']},ensure_ascii=False));return 0


# Final fixed-regime diagnostics and supervisor-report assets; no new model.
SUPERVISOR_CONFIG=p.ROOT/'data/input/nlsp_beta15_pre_supervisor_checks.json'
SUPERVISOR_OUTPUT=p.ROOT/'results/nlsp_beta15_pre_supervisor_checks'
SUPERVISOR_SMOKE=p.ROOT/'results/_smoke/nlsp_beta15_pre_supervisor_checks'


def supervisor_prepare():
    cfg=p.read(SUPERVISOR_CONFIG)
    if cfg['degree']!=14 or cfg['A']!=.005 or cfg['phase_Psi']!=0 or cfg['numerical_budget_seconds']>600 or cfg['new_full_nonlinear_calls'] or cfg['new_FEM_jobs'] or cfg['new_HB_Hill_eigensolves']:raise ValueError('Fixed final-check scope required')
    for path in SUPERVISOR_OUTPUT.glob('*/provenance.json') if SUPERVISOR_OUTPUT.exists() else []:
        item=p.read(path)
        if item['config']!=cfg:raise ValueError('Different final-check scope')
        return path.parent,item,p.cache(path.parent)
    sources={};artifacts={}
    for key in ('p14_parent','p12_Hill_parent','p12_dynamic_parent'):
        root=p.ROOT/cfg[key];p.cache(root);sources[key]={'path':cfg[key],'sha256':p.sha(root/'manifest.json')}
        for file in root.iterdir():
            if file.is_file():artifacts[file.relative_to(p.ROOT).as_posix()]=p.sha(file)
    prior=p.read(p.ROOT/cfg['p14_parent']/'provenance.json')
    sources={**prior['source_manifests'],**sources};artifacts={**prior['source_artifacts'],**artifacts}
    if any(p.sha(p.ROOT/name)!=prior['code_hashes'][name] for name in old.f.s.r.KERNELS):raise ValueError('Physical kernel changed')
    names=[SUPERVISOR_CONFIG.relative_to(p.ROOT).as_posix(),Path(__file__).resolve().relative_to(p.ROOT).as_posix()]
    item={'schema':'nlsp-beta15-pre-supervisor-checks-v1','config':cfg,'source_manifests':sources,'source_artifacts':artifacts,
        'code_hashes':{**p.code_hashes(),**{name:p.sha(p.ROOT/name) for name in names}},'initial_state':p.read(SUPERVISOR_SMOKE/'initial_state.json'),
        'scientific_contract':'fixed p14 neutral residual audit and one conditional displacement-only tangent IVP, no full nonlinear/HB/Hill/FEM solves'}
    fp=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=SUPERVISOR_OUTPUT/fp;b.mkdir(parents=True)
    p.write(b/'config.json',cfg);p.write(b/'provenance.json',item);folder=b/'execution_code';folder.mkdir()
    for name in item['code_hashes']:shutil.copyfile(p.ROOT/name,folder/Path(name).name)
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'numerical_seconds':p.read(SUPERVISOR_SMOKE/'preflight.json')['scientific_seconds'],
        'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'FEM':0,'HB_Hill_eigensolves':0},'cases':{},'attempts':[],
        'stages':{'neutral':'NOT_RUN','projection':'NOT_RUN','tangent':'NOT_RUN','report':'NOT_RUN'},'deferred':cfg['deferred']}
    p.finalize(b,item,s);return b,item,s


def supervisor_problem(item):
    parent=p.ROOT/item['config']['p14_parent'];prior=p.read(parent/'provenance.json');problem=p14_problem(prior)
    d12,d,model,Phi12,w12,U12,out12,ow12=problem
    with np.load(parent/'linear_reference.npz') as z:Phi=z['Phi'].copy();omega=z['omega'].copy();out=z['out_Phi'].copy();ow=z['out_omega'].copy()
    with np.load(parent/'HB13.npz') as z:U=z['U'].copy()
    return d,model,Phi,U,p.read(parent/'HB13.json')['Omega'],omega,out,ow


def orbit_with_jerk(U,Omega,Phi,phases):
    q,v,a=hb_phases(U,Omega,Phi,phases);n=np.arange(len(U))
    jerk=Omega**3*np.sin(phases[:,None]*n)@(n[:,None]**3*U)@Phi.T
    return q,v,a,jerk


def supervisor_neutral_audit(b,item,s,problem,deadline):
    if (b/'neutral_diagnostics.json').exists():return
    import copy
    d,model,Phi,U,Omega,omega,out,ow=problem;cfg=item['config'];parent=p.ROOT/cfg['p14_parent'];start=time.perf_counter()
    with np.load(parent/'PL_Fourier.npz') as z:F=z['F'].copy();orders=z['orders'].copy();R=z['R'].copy();L=z['L'].copy()
    holder=copy.copy(d);holder.out=d.plane;action=old.f.QuadraticAction(holder,model.T4,model.V4)
    count=cfg['neutral_phase_samples'];phases=np.arange(count)*2*np.pi/count;q,v,a,jerk=orbit_with_jerk(U,Omega,Phi,phases)
    neutral=[];canonical_residual=[];mapped_rdot=[];identity_errors=[];velocity_residual=[];mass_residual=[];off_orbit=[]
    for qi,vi,ai,ji in zip(q,v,a,jerk):
        if time.perf_counter()>deadline:raise TimeoutError('Neutral audit budget')
        M=d.mass_matrix(qi);actual=d.acceleration(qi,vi);r=ai-actual
        Q,V,_=full_massform_tangents(d,qi,vi,ai);Rdot=Q@vi+V@ai+M@ji
        blocks=action.blocks(qi,vi,ai);Av,Ac,H,S,Sdot,J=old.f.action_operators(blocks)
        rdot=np.linalg.solve(blocks['M'],Rdot[d.plane]-blocks['Mdot']@r[d.plane])
        full=d.jacobian(0,np.r_[qi,vi]);Jv=full[np.ix_(d.state_plane,d.state_plane)]
        z=np.r_[vi[d.plane],ai[d.plane]];zdot=np.r_[ai[d.plane],ji[d.plane]]
        Z=L@S@z;cres=L@(Sdot@z+S@zdot-Ac@S@z)
        mapder=L@S@np.r_[np.zeros(len(d.plane)),rdot]
        delta=L@((Sdot+S@Jv)-Ac@S)@z
        identity_errors.append(L@S@(zdot-Jv@z)-mapder)
        neutral.append(Z);canonical_residual.append(cres);mapped_rdot.append(mapder);off_orbit.append(delta)
        velocity_residual.append(np.r_[np.zeros(len(d.plane)),r[d.plane]])
        mass_residual.append(M@ai+d.inertial_terms(qi,vi)+d.potential(qi)['gradient'])
    neutral,canonical_residual,mapped_rdot,off_orbit,identity_errors=map(np.array,(neutral,canonical_residual,mapped_rdot,off_orbit,identity_errors))
    scale=Omega*np.sqrt(np.mean(np.sum(abs(neutral)**2,axis=1)));n=np.arange(-33,34)
    C=np.exp(-1j*n[:,None]*phases)/count
    coefficients=C@neutral;rescoeff=C@canonical_residual;records=[]
    for N in cfg['neutral_Fourier_controls']:
        indices=np.flatnonzero(abs(n)<=N);z=coefficients[indices];pred=hill_action(F,orders,Omega,N,z.ravel()).reshape(len(indices),-1)
        error=pred+rescoeff[indices]
        records.append({'N':N,'raw_Hill_neutral_relative_residual':float(np.linalg.norm(pred)/(Omega*np.linalg.norm(z))),
            'orbit_residual_corrected_Hill_relative_error':float(np.linalg.norm(error)/(Omega*np.linalg.norm(z))),
            'neutral_discarded_coefficient_relative_norm':float(np.linalg.norm(coefficients[abs(n)>N])/np.linalg.norm(coefficients))})
    # Independent central differences of the ORIGINAL first-order residual.
    ph=np.arange(cfg['independent_phase_count'])*2*np.pi/cfg['independent_phase_count']+cfg['independent_phase_offset'];qq,vv,aa,jj=orbit_with_jerk(U,Omega,Phi,ph)
    fd_records=[]
    for step in cfg['residual_derivative_phase_steps']:
        qp,vp,ap,_=orbit_with_jerk(U,Omega,Phi,ph+step);qm,vm,am,_=orbit_with_jerk(U,Omega,Phi,ph-step);errors=[]
        for qi,vi,ai,ji,qpi,vpi,api,qmi,vmi,ami in zip(qq,vv,aa,jj,qp,vp,ap,qm,vm,am):
            Jfull=d.jacobian(0,np.r_[qi,vi]);expected=np.r_[ai,ji]-Jfull@np.r_[vi,ai]
            residual_plus=np.r_[np.zeros(d.ndof),api-d.acceleration(qpi,vpi)]
            residual_minus=np.r_[np.zeros(d.ndof),ami-d.acceleration(qmi,vmi)]
            observed=(residual_plus-residual_minus)*Omega/(2*step)
            S=old.f.action_operators(action.blocks(qi,vi,ai))[3]
            errors.append(np.linalg.norm(L@S@(expected-observed)[d.state_plane])/scale)
        fd_records.append({'phase_step':step,'first_order_residual_derivative_identity_relative_error_max':float(max(errors))})
    # Identify the existing near-zero dominant PL eigenvector without a new solve.
    with np.load(parent/'PL_k108_N7.npz') as z:vector=z['dominant_eigenvector'].reshape(15,216).copy()
    row=p.read(parent/'PL_k108_N7.json');raw=row['dominant_raw_lambda'];fold=row['dominant_lambda'];shift=int(round((raw['imag']-fold['imag'])/Omega))
    shifted=np.zeros_like(coefficients);shifted[np.arange(-7,8)+shift+33]=vector
    mac=float(abs(np.vdot(coefficients,shifted))**2/(np.linalg.norm(coefficients)**2*np.linalg.norm(shifted)**2))
    normalized_harmonics=np.linalg.norm(rescoeff,axis=1)/scale;dominant=int(abs(n[np.argmax(normalized_harmonics)]))
    diag={'raw_source_N13_residual':p.read(parent/'PL_checks.json')['time_shift_neutral_relative_residual'],
        'canonical_pointwise_neutral_relative_RMS':float(np.linalg.norm(canonical_residual)/np.sqrt(count)/scale),
        'mapped_orbit_residual_derivative_relative_RMS':float(np.linalg.norm(mapped_rdot)/np.sqrt(count)/scale),
        'off_orbit_canonical_difference_relative_RMS':float(np.linalg.norm(off_orbit)/np.sqrt(count)/scale),
        'analytic_first_order_identity_relative_error_max':float(np.linalg.norm(identity_errors,axis=1).max()/scale),
        'canonical_identity_decomposition_relative_error_max':float(np.linalg.norm(canonical_residual-mapped_rdot-off_orbit,axis=1).max()/scale),
        'Fourier_controls':records,'independent_FD_controls':fd_records,'neutral_scale_Omega_RMS_Z':float(scale),
        'dominant_residual_derivative_harmonic':dominant,'saved_near_zero_PL_mode_time_shift_MAC':mac,
        'saved_non_neutral_sigma':p.read(parent/'PL_checks.json')['non_neutral_sigma'],
        'original_gate_unchanged':cfg['neutral_relative_gate'],
        'variational_identity_gate':bool(max(np.linalg.norm(identity_errors,axis=1).max()/scale,np.linalg.norm(canonical_residual-mapped_rdot-off_orbit,axis=1).max()/scale,
            records[-1]['orbit_residual_corrected_Hill_relative_error'])<=cfg['neutral_relative_gate']),
        'raw_exact_periodic_neutral_gate':bool(records[-1]['raw_Hill_neutral_relative_residual']<=cfg['neutral_relative_gate']),
        'PL_status':'PL_STABILITY_REMAINS_PARTIAL','no_new_Hill_eigensolve':True}
    if diag['raw_exact_periodic_neutral_gate'] and diag['variational_identity_gate']:diag['PL_status']='PL_NEUTRAL_CONTROL_RESOLVED'
    p.write(b/'neutral_diagnostics.json',diag)
    np.savez_compressed(b/'neutral_audit.npz',phases=phases,canonical_time_shift=neutral,canonical_time_shift_residual=canonical_residual,
        mapped_first_order_residual_derivative=mapped_rdot,canonical_off_orbit_difference=off_orbit,velocity_residual=np.array(velocity_residual),
        massform_residual=np.array(mass_residual),Fourier_orders=n,neutral_coefficients=coefficients,residual_derivative_coefficients=rescoeff,
        shifted_saved_PL_mode=shifted,identity_errors=identity_errors)
    s['stages']['neutral']=diag['PL_status'];s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)
    print(json.dumps({'neutral_status':diag['PL_status'],'identity_gate':diag['variational_identity_gate'],'dominant_harmonic':dominant,'MAC':mac,'controls':records},ensure_ascii=False),flush=True)


def supervisor_projection(b,item,s,problem):
    if (b/'OUT2_prediction.json').exists():return
    d,model,Phi,U,Omega,omega,out,ow=problem;cfg=item['config'];parent=p.ROOT/cfg['p14_parent'];start=time.perf_counter()
    with np.load(parent/'physical_Floquet_and_displacement.npz') as z:
        ru=z['OUT_unit'].copy();rs=z['stable_OUT_unit'].copy();x=z['displacement_only_OUT'].copy();S=z['S'].copy();J=z['J'].copy();E=z['metric'].copy()
    den=(S@rs)@J@(S@ru);cu=float((S@rs)@J@(S@x)/den);cs=float((S@ru)@J@(S@x)/(-den));rest=x-cu*ru-cs*rs
    components=[x,cu*ru,cs*rs,rest];scales=[]
    for comp in components:
        position=np.zeros(d.ndof);position[d.out]=comp[:len(d.out)];scales.append(out_polynomial_extrema(d,position,1)[0])
    scales=np.array(scales);refnorm=float(np.sqrt(x@E@x));sigma=p.read(parent/'OUT_checks.json')['sigma'];T=2*np.pi/Omega
    source=p.read(parent/'projection_checks.json');gain=np.exp(sigma*T*np.array(cfg['forecast_periods']))
    prediction={'initial_OUT2_E0_norm':refnorm,'growing_coefficient_unit_E0':cu,'stable_coefficient_unit_E0':cs,
        'projection_reference_difference':cu-source['records'][-1]['unstable_coefficient'],'projection_uncertainty_indicator':source['projection_uncertainty_indicator'],
        'component_order':['total initial','unstable initial','stable initial','remaining initial'],
        'initial_physical_maxima_v_Phi_psi':scales.tolist(),'initial_growing_combined_scale_over_total':float(scales[1].max()/scales[0].max()),
        'initial_growing_v_scale_over_total_v':float(scales[1,0]/scales[0,0]),
        'initial_component_E0_norms':[float(np.sqrt(comp@E@comp)) for comp in components],
        'remaining_E0_norm_over_initial':float(np.sqrt(rest@E@rest)/refnorm),'not_energy_fractions':True,
        'forecast_periods':cfg['forecast_periods'],'unstable_stroboscopic_gains':gain.tolist(),
        'unstable_stroboscopic_physical_maxima':(gain[:,None]*scales[1]).tolist(),
        'forecast_not_a_total_OUT_response':True,'forecast_10T_20T_not_integrated':True,
        'spectral_total_response_unavailable_reason':'only critical unstable/reciprocal eigenvectors were saved; center-mode eigenvectors and a complete Floquet decomposition are unavailable',
        'short_tangent_is_needed':True,'initial_OUT_velocities_zero':bool(np.max(abs(x[len(d.out):]))==0)}
    p.write(b/'OUT2_prediction.json',prediction);np.savez_compressed(b/'OUT2_initial_components.npz',total=x,unstable=cu*ru,stable=cs*rs,remaining=rest,
        ru=ru,rs=rs,metric=E,S=S,J=J,canonical_displacement=S@x)
    s['stages']['projection']='PASS';s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)
    print(json.dumps({'initial_growing_scale_ratio':prediction['initial_growing_combined_scale_over_total'],
        'initial_growing_v_ratio':prediction['initial_growing_v_scale_over_total_v'],'remaining_norm_ratio':prediction['remaining_E0_norm_over_initial'],'gain_forecast':gain.tolist()}),flush=True)


class SavedCanonicalFourier:
    """Real evaluation of the previously verified canonical Fourier operator."""
    def __init__(self,F,orders,Omega):
        zero=np.flatnonzero(orders==0)[0];positive=np.flatnonzero(orders>0)
        self.harmonics=orders[positive];self.Omega=Omega
        self.matrices=np.concatenate((F[zero:zero+1].real,2*F[positive].real,-2*F[positive].imag),axis=0)
    def __call__(self,t):
        phase=self.harmonics*self.Omega*t
        return np.einsum('m,mij->ij',np.r_[1.,np.cos(phase),np.sin(phase)],self.matrices,optimize=True)


def supervisor_tangent_control(b,item,s,problem,deadline):
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    if (b/'OUT2_tangent_checks.json').exists():return
    d,model,Phi,U,Omega,omega,out,ow=problem;cfg=item['config'];parent=p.ROOT/cfg['p14_parent'];start=time.perf_counter();T=2*np.pi/Omega
    with np.load(parent/'OUT_Fourier.npz') as z:F=z['F'].copy();orders=z['orders'].copy();R=z['R'].copy();L=z['L'].copy()
    with np.load(b/'OUT2_initial_components.npz') as z:x0=z['total'].copy();S0=z['S'].copy();E=z['metric'].copy();ru=z['ru'].copy();rs=z['rs'].copy()
    action=old.f.QuadraticAction(d,model.T4,model.V4);operator=SavedCanonicalFourier(F,orders,Omega)
    mapped=L@S0@x0;scale=float(np.linalg.norm(mapped));basis=(mapped/scale)[:,None];metric=np.eye(len(mapped))
    # Audit the saved representation against the FULL analytic Jacobian at
    # independent phases before using it in a direct time-domain IVP.
    operator_errors=[];conjugate_errors=[]
    for phase in (.219,1.719,4.319):
        q,v,a=hb_phases(U,Omega,Phi,np.array([phase]));blocks=action.blocks(q[0],v[0],a[0]);_,_,_,S,Sdot,_=old.f.action_operators(blocks)
        full=d.jacobian(0,np.r_[q[0],v[0]])[np.ix_(d.state_out,d.state_out)]
        actual=L@(Sdot+S@full)@np.linalg.inv(S)@R;value=operator(phase/Omega)
        operator_errors.append(np.linalg.norm(actual-value)/np.linalg.norm(actual))
        direct=np.einsum('m,mij->ij',np.exp(1j*orders*phase),F);conjugate_errors.append(np.linalg.norm(direct-value)/np.linalg.norm(direct))
    probe=np.linspace(.173,2.9*T,30);before=time.perf_counter()
    for t in probe:operator(t)@basis[:,0]
    rhs_seconds=(time.perf_counter()-before)/len(probe)
    oldrun=p.read(p.ROOT/cfg['p12_dynamic_parent']/'tangent_primary_3T/integration.json')
    forecast=1.5*rhs_seconds*oldrun['nfev']+20
    policy={'reason':'saved critical pair cannot reconstruct the remaining 35% E0-norm center component',
        'source':'saved full81-mode canonical Fourier operator, 96 phases,m=-16..16,from full analytic Jacobian/action',
        'pointwise_relative_Jacobian_comparison':float(max(operator_errors)),'real_Fourier_evaluation_relative_error':float(max(conjugate_errors)),
        'rhs_seconds':rhs_seconds,'forecast_3T_seconds':forecast,'remaining_seconds':600-s['numerical_seconds'],
        'physical_initial_scale':5e-5,'OUT_velocities_zero':True,'canonical_initial_norm':scale,
        'no_new_HB_Hill_or_planar_ODE':True,'same_frozen_model':True}
    p.write(b/'OUT2_tangent_policy.json',policy)
    s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)
    if max(operator_errors)>cfg['canonical_identity_relative_gate'] or forecast+cfg['cost_reserve_seconds']>600-s['numerical_seconds']:
        s['stages']['tangent']='NOT_RUN_OPERATOR_OR_RESOURCE_GATE';p.finalize(b,item,s);return
    def base(t):
        q,v,_=hb_phases(U,Omega,Phi,np.array([Omega*t]));return np.r_[q[0],v[0]]
    def override(t,state):return operator(t)
    tc={**cfg,'rtol':cfg['tangent_rtol'],'atol_relative':cfg['tangent_atol_relative']}
    label='OUT2_displacement_only_3T'
    if label not in s['cases']:
        analysis.tangent_attempt(b,{**item,'config':tc},s,label,d,parent,np.linspace(0,3*T,3*cfg['samples_per_period']+1),basis,metric,deadline,
            method=cfg['tangent_method'],sampled_dense_only=True,operator_override=override,base_evaluator=base)
        if s['hard_stop']:return
    with np.load(b/label/'tangent.npz') as z:times=z['times'].copy();Y=scale*z['vectors'][:,:,0].copy()
    # A true physical-state restart in the canonical variables, original seed
    # unchanged; extension does not create another initial condition.
    extension_forecast=1.5*s['cases'][label]['seconds']*2/3+15
    den=(S0@rs)@np.block([[np.zeros((81,81)),np.eye(81)],[-np.eye(81),np.zeros((81,81))]])@(S0@ru)
    J=np.block([[np.zeros((81,81)),np.eye(81)],[-np.eye(81),np.zeros((81,81))]])
    left=(S0@rs)@J/den;prediction=p.read(b/'OUT2_prediction.json');cu=prediction['growing_coefficient_unit_E0'];sigma=p.read(parent/'OUT_checks.json')['sigma']
    strobe=np.arange(4)*cfg['samples_per_period'];actual=Y[strobe]@R.T@left
    expected=cu*(-1.)**np.arange(4)*np.exp(sigma*times[strobe]);gate=np.max(abs(actual-expected))/max(abs(expected))<cfg['projection_relative_gate']
    p.write(b/'OUT2_extension_policy.json',{'forecast_2T_seconds':extension_forecast,'remaining_seconds':600-s['numerical_seconds'],
        'projected_growth_3T_gate':bool(gate),'stroboscopic_relative_error':float(np.max(abs(actual-expected))/max(abs(expected))),
        'conditional_extension':bool(gate and extension_forecast+cfg['cost_reserve_seconds']<600-s['numerical_seconds'])})
    if gate and extension_forecast+cfg['cost_reserve_seconds']<600-s['numerical_seconds']:
        continued='OUT2_continue_3T_to_5T'
        if continued not in s['cases']:
            with np.load(b/label/'restart_state.npz') as z:restart=z['actual_final_matrix'].copy()
            analysis.tangent_attempt(b,{**item,'config':tc},s,continued,d,parent,np.linspace(3*T,5*T,2*cfg['samples_per_period']+1),basis,metric,deadline,
                method=cfg['tangent_method'],sampled_dense_only=True,t_start=3*T,initial_matrix=restart,operator_override=override,base_evaluator=base)
            if s['hard_stop']:return
        with np.load(b/continued/'tangent.npz') as z:times=np.r_[times,z['times'][1:]];Y=np.r_[Y,scale*z['vectors'][1:,:,0]]
    start=time.perf_counter();physical=[];extrema=[];total=[];grow=[];rest=[]
    with np.load(parent/'OUT_k81_N9.npz') as z:coeff=z['dominant_eigenvector'].reshape(19,162).copy()
    raw=p.read(parent/'OUT_k81_N9.json')['dominant_raw_lambda'];lam=complex(raw['real'],raw['imag'])
    with np.load(parent/'physical_Floquet_and_displacement.npz') as z:phase=z['unstable_phase'].item()
    for t,value in zip(times,Y):
        q,v,a=hb_phases(U,Omega,Phi,np.array([Omega*t]));S=old.f.action_operators(action.blocks(q[0],v[0],a[0]))[3]
        observed=np.linalg.solve(S,R@value);unstable=cu*np.real(phase*np.exp(lam*t)*np.linalg.solve(S,R@(np.exp(1j*np.arange(-9,10)*Omega*t)@coeff)))
        physical.append(observed);grow.append(unstable);rest.append(observed-unstable)
        qout=np.zeros(d.ndof);qout[d.out]=observed[:81];extrema.append(out_polynomial_extrema(d,qout,1)[0])
    physical,grow,rest,extrema=map(np.array,(physical,grow,rest,extrema));periods=int(round(times[-1]/T));strobe=np.arange(periods+1)*cfg['samples_per_period']
    coeffs=Y[strobe]@R.T@left;ampl=np.abs(coeffs/coeffs[0]);expected_gain=np.exp(sigma*times[strobe])
    norm=analysis.energy_norm(physical,E);grow_norm=analysis.energy_norm(grow,E);rest_norm=analysis.energy_norm(rest,E)
    diag={'periods':periods,'stroboscopic_growth_projection_amplification':ampl.tolist(),'Hill_component_prediction':expected_gain.tolist(),
        'finite_interval_projected_lograte':float(np.polyfit(times[strobe],np.log(ampl),1)[0]),'not_a_full_nonlinear_Floquet_exponent':True,
        'maximum_v_Phi_psi':extrema.max(axis=0).tolist(),'maximum_v_over_A':float(extrema[:,0].max()/cfg['A']),
        'initial_E0_norm':float(norm[0]),'maximum_total_E0_amplification':float(norm.max()/norm[0]),
        'stroboscopic_total_E0_amplification':(norm[strobe]/norm[0]).tolist(),
        'growing_and_remaining_E0_norm_ratio_at_end':float(grow_norm[-1]/rest_norm[-1]),
        'growth_projection_relative_error':float(max(abs(ampl-expected_gain))/max(expected_gain)),
        'time_maxima_sampled_T_over_256':True,'background':'saved p14 HB13 periodic orbit; canonical Fourier operator including mass/inertia/Sdot',
        'one_selected_trajectory_with_conditional_restart':True}
    np.savez_compressed(b/'OUT2_tangent_physical.npz',times=times,physical_OUT=physical,canonical_energy_state=Y,unstable_component=grow,remaining_component=rest,
        physical_extrema=extrema,total_E0_norm=norm,unstable_E0_norm=grow_norm,remaining_E0_norm=rest_norm,stroboscopic_times=times[strobe],
        projected_growth_coefficient=coeffs,metric=E)
    p.write(b/'OUT2_tangent_checks.json',diag);s['stages']['tangent']='PASS' if diag['growth_projection_relative_error']<cfg['projection_relative_gate'] else 'PARTIAL'
    s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)
    print(json.dumps({'tangent':diag,'seconds':s['numerical_seconds']},ensure_ascii=False),flush=True)


def supervisor_report_data(b,item,s):
    """Saved-data scientific assets, no physical solver or new integration."""
    if (b/'report_data.npz').exists():return
    start=time.perf_counter();cfg=item['config'];parent=p.ROOT/cfg['p14_parent'];dynamic=p.ROOT/cfg['p12_dynamic_parent']
    geometry=p.read(dynamic/'actual_geometry_and_coefficients.json');bases=np.array(geometry['local_bases']);lengths=np.array(geometry['lengths'])
    with np.load(parent/'physical_observations.npz') as z:pl=z['p14_initial_profiles'].copy()
    with np.load(parent/'physical_Floquet_and_displacement.npz') as z:out=z['p14_profiles'].copy()
    out*=p.read(parent/'physical_observations.json')['p14_initial_Floquet_physical_scale']
    rest=np.array([(-length*B[:,0])[None,:]+np.linspace(0,length,len(pl[0]))[:,None]*B[:,0] for length,B in zip(lengths,bases)])
    pl_displacement=np.array([fields[:,:3]@B.T for fields,B in zip(pl,bases)])
    out_displacement=np.array([fields[:,:3]@B.T for fields,B in zip(out,bases)])
    extension=p.read(dynamic/'extension_5T_checks.json');hill=p.read(dynamic/'Floquet_reconstruction.json')
    periods=np.arange(6);old_Hill=np.exp(hill['sigma_Hill']*hill['T']*periods)
    arrays={'bases':bases,'lengths':lengths,'rest_centerlines':rest,'PL_displacement':pl_displacement,'OUT_displacement':out_displacement,
        'PL_fields':pl,'OUT_fields_at_scale_1percent_A':out,'old_periods':periods,'old_Hill_gain':old_Hill,
        'old_tangent_gain':np.array(extension['stroboscopic_tangent_amplification']),
        'old_full_gain':np.array(extension['stroboscopic_nonlinear_amplification'])}
    if (b/'OUT2_tangent_physical.npz').exists():
        with np.load(b/'OUT2_tangent_physical.npz') as z:
            arrays.update(new_times=z['times'].copy(),new_physical_extrema=z['physical_extrema'].copy(),new_total_E0=z['total_E0_norm'].copy(),
                new_unstable_E0=z['unstable_E0_norm'].copy(),new_remaining_E0=z['remaining_E0_norm'].copy())
        T=2*np.pi/p.read(parent/'HB13.json')['Omega'];arrays['new_tau']=arrays['new_times']/T
        extras={}
        for periods0 in (3,5):
            ids=arrays['new_tau']<=periods0+1e-12
            extras[f'maximum_v_over_A_through_{periods0}T']=float(arrays['new_physical_extrema'][ids,0].max()/cfg['A'])
            extras[f'maximum_total_E0_gain_through_{periods0}T']=float(arrays['new_total_E0'][ids].max()/arrays['new_total_E0'][0])
        check=p.read(b/'OUT2_tangent_checks.json');check.update(extras);p.write(b/'OUT2_tangent_checks.json',check)
    np.savez_compressed(b/'report_data.npz',**arrays)
    p.write(b/'report_figure_conventions.json',{'geometry':geometry,'PL_visual_magnification':15,'OUT_visual_magnification':2000,
        'OUT_initial_physical_scale':5e-5,'no_phase_or_amplitude_fit':True,
        'centroid_shape_uses_saved_full_field_profiles':True,'rotations_retained_in_data_but_not_all_displayed_on_centerlines':True})
    s['numerical_seconds']+=time.perf_counter()-start;p.finalize(b,item,s)


def supervisor_figures(b):
    """Render only completed scientific data, safe for temporary-copy replay."""
    if not (b/'report_data.npz').exists():return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Arc
    with np.load(b/'report_data.npz') as z:data={key:z[key].copy() for key in z.files}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10})
    def save(fig,name):
        fig.savefig(b/(name+'.png'),dpi=200,bbox_inches='tight')
        fig.savefig(b/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None});plt.close(fig)
    B=data['bases'];L=data['lengths'];rest=data['rest_centerlines'];fig,ax=plt.subplots(figsize=(8.5,2.5))
    for i,line in enumerate(rest):
        ax.plot(line[:,0],line[:,1],color='black',lw=2)
        clamp=line[0,:2];normal=B[i,:2,1];tangent=B[i,:2,0]
        bar=np.array([clamp-.07*normal,clamp+.07*normal]);ax.plot(bar[:,0],bar[:,1],color='black',lw=2)
        for offset in np.linspace(-.06,.06,5):
            a=clamp+offset*normal;end=a-.035*tangent+.02*normal;ax.plot([a[0],end[0]],[a[1],end[1]],color='black',lw=.8)
        midpoint=line[len(line)//2,:2]
        for vector,label in [(tangent,rf'$t_{i+1},s_{i+1}$'),(normal,rf'$n_{i+1}$')]:
            end=midpoint+.14*vector;ax.annotate('',xy=end,xytext=midpoint,arrowprops={'arrowstyle':'->','lw':1})
            ax.text(*(end+.025*normal),label,fontsize=10)
        label_point=midpoint-.13*normal if i==0 else midpoint-.16*normal-.24*tangent
        ax.text(*label_point,rf'$L_{i+1}={L[i]:.1f}$',ha='center',fontsize=11)
    ax.plot(0,0,'ko',ms=5);ax.plot([0,.38],[0,0],':',color='gray');ax.add_patch(Arc((0,0),.55,.55,theta1=0,theta2=15,lw=1.2))
    ax.text(.30,.032,r'$\beta=15^\circ$',fontsize=11);ax.text(-.08,-.15,r'$\otimes\ k=-E_Z,\quad v\parallel k$',ha='center',fontsize=10)
    ax.annotate('',xy=(-.70,-.27),xytext=(-.85,-.27),arrowprops={'arrowstyle':'->'});ax.text(-.67,-.27,r'$E_X$',va='center')
    ax.annotate('',xy=(-.85,-.13),xytext=(-.85,-.27),arrowprops={'arrowstyle':'->'});ax.text(-.85,-.10,r'$E_Y$',ha='center')
    ax.set_aspect('equal');ax.set_xlim(-1.05,1.20);ax.set_ylim(-.33,.48);ax.axis('off');save(fig,'supervisor_geometry')
    fig,(ax,ax3)=plt.subplots(1,2,figsize=(9,2.8))
    for i,base in enumerate(rest):
        shape=base+15*data['PL_displacement'][i];ax.plot(base[:,0],base[:,1],':',color='gray');ax.plot(shape[:,0],shape[:,1],lw=1.8,label=f'Плечо {i+1}')
        out=base+2000*data['OUT_displacement'][i];ax3.plot(base[:,0],base[:,2],':',color='gray');ax3.plot(out[:,0],out[:,2],lw=1.8)
    ax.set(xlabel='X',ylabel='Y',title='а) Плоская орбита, масштаб ×15');ax.set_aspect('equal');ax.legend(fontsize=8)
    ax3.set(xlabel='X',ylabel='Z',title='б) OUT-форма: проекция XZ, ×2000')
    ax3.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax3.grid(alpha=.2)
    fig.tight_layout();save(fig,'supervisor_physical_modes')
    fig,ax=plt.subplots(1,2,figsize=(9,3.1));tau=data['old_periods']
    ax[0].plot(tau,data['old_Hill_gain'],'-',color='black',label='Хилл, p12')
    ax[0].plot(tau,data['old_tangent_gain'],'o',mfc='none',ms=7,label='Касательная динамика, p12')
    ax[0].plot(tau,data['old_full_gain'],'x',ms=6,label='Полная нелинейная, p12')
    ax[0].set(xlabel='t/T',ylabel='Усиление Floquet-компоненты');ax[0].legend(fontsize=8)
    if 'new_tau' in data:ax[1].plot(data['new_tau'],100*data['new_physical_extrema'][:,0]/.005,color='tab:blue')
    ax[1].set(xlabel='t/T',ylabel='max |v| / A, %',title='OUT2: начальные скорости нулевые');ax[1].set_ylim(bottom=0)
    fig.tight_layout();save(fig,'supervisor_growth_confirmation')


def supervisor_main(argv=None):
    parser=argparse.ArgumentParser(description='Fixed neutral/OUT2 final checks and saved report assets')
    mode=parser.add_mutually_exclusive_group(required=True);mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--through',choices=['audit','tangent','all'],default='all');args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=p.cache(b)
        if args.plot_only:supervisor_figures(b);p.finalize(b,p.read(b/'provenance.json'),s)
    else:
        b,item,s=supervisor_prepare()
        if not s['completed'] and not s['hard_stop']:
            start=time.perf_counter();initial=s['numerical_seconds'];deadline=start+600-initial-item['config']['cost_reserve_seconds']
            try:
                science_done=all((b/name).exists() for name in ('neutral_diagnostics.json','OUT2_prediction.json','OUT2_tangent_checks.json','report_data.npz'))
                if not science_done:
                    problem=supervisor_problem(item);s['numerical_seconds']+=time.perf_counter()-start
                    supervisor_neutral_audit(b,item,s,problem,deadline);supervisor_projection(b,item,s,problem)
                    if args.through in ('tangent','all'):supervisor_tangent_control(b,item,s,problem,deadline)
                    if args.through=='all':supervisor_report_data(b,item,s,problem)
                if args.through=='all':
                    supervisor_figures(b);s['stages']['report']='ASSETS_COMPLETE'
                    s.update(completed=True,overall='PRE_SUPERVISOR_CHECKS_COMPLETE_WITH_PL_PARTIAL')
            except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as error:
                s.update(hard_stop=True,overall='NUMERICAL_PARTIAL',stop_reason=type(error).__name__+': '+str(error));s['numerical_seconds']=max(s['numerical_seconds'],initial+time.perf_counter()-start)
            p.finalize(b,item,s)
    print(json.dumps({'bundle':str(b),'stages':s['stages'],'seconds':s['numerical_seconds'],'status':s['overall']},ensure_ascii=False));return 0
