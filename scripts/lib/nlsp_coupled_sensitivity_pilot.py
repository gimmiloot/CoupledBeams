"""Bounded coupled-rod initial-perturbation pilot; no native FEM capability.

One new workflow/contract: two lifted arms, Jacobian tangent dynamics, and a
finite-dimensional initial-state subspace. This reuses frozen action kernels,
the existing Newton/Radau/dense-record machinery, and introduces no model law.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib,json,os,shutil,subprocess,sys,time
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from scipy.linalg import cho_factor,cho_solve,eigh
from scipy.integrate import Radau
from scipy.sparse import csr_matrix,kron,eye
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_spatial_dynamics import SpatialGalerkin,FIELDS
from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import CoupledSpatialGalerkin,rectangular_coefficients
from scripts.analysis.verify_nlsp_nonlinear_static_3d_fem import fem2_static_newton
from scripts.lib.nlsp_spatial_1d_program import STATIC_POLICY
from scripts.lib.nlsp_fem3b_continuation import save_dense_records,evaluate_dense_records

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'data/input/nlsp_coupled_planar_sensitivity_pilot.json'
OUTPUT=ROOT/'results/nlsp_coupled_planar_sensitivity_pilot'
AUTH='explicit_user_coupled_planar_sensitivity_2026_10_10'
VERSION='coupled-planar-initial-perturbation-pilot-v1'

def plain(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,dict):return {str(k):plain(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [plain(v) for v in x]
    return x
def read(path):return json.loads(Path(path).read_text(encoding='utf8'))
def write(path,data):Path(path).write_text(json.dumps(plain(data),ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')
def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def validate_config(c):
    expected=read(CONFIG)
    if c!=expected or c['authorization']!=AUTH:raise ValueError('Only the frozen two-geometry pilot is authorized')
    if c['new_FEM_jobs'] or c['parameter_maps'] or c['Floquet'] or c['automatic_retry']:raise ValueError('No expanded or native programme permitted')
    return c

def sources(c):
    selected={};manifests={}
    for name,path,files in [('action',c['action_bundle'],['result.json']),
        ('previous_spatial',c['previous_spatial_bundle'],['summary.json','stage_a_checks.json','config.json'])]:
        source=ROOT/path;manifest=read(source/'manifest.json');manifests[name]={'path':path,'sha256':sha(source/'manifest.json')}
        hashes=manifest.get('artifact_hashes',manifest.get('artifacts',{}))
        for file in files:
            p=source/file
            if hashes.get(file)!=sha(p):raise ValueError('Missing/corrupt saved source '+str(p))
            selected[p.relative_to(ROOT).as_posix()]=sha(p)
    if read(ROOT/c['previous_spatial_bundle']/'summary.json')['overall']!='NUMERICAL_PARTIAL':raise ValueError('Unexpected previous spatial status')
    raw=read(ROOT/c['action_bundle']/'result.json')['polynomials']
    model=SimpleNamespace(T4=rod.Polynomial.deserialize(raw['T4']),V4=rod.Polynomial.deserialize(raw['V4']),
        residual_a=tuple(rod.Polynomial.deserialize(p) for p in raw['residuals_A']))
    return model,selected,manifests

def make_disc(c,name,p,model,beta=None,lengths=None):
    L=tuple(c['lengths'] if lengths is None else lengths)
    coeff=[rectangular_coefficients(b,h,l)[0] for (b,h),l in zip(c['sections'][name],L)]
    return CoupledSpatialGalerkin(coeff,p,L,c['beta_deg'] if beta is None else beta,model)

def code_hashes():
    files=['scripts/lib/weakly_nonlinear_spatial_rod.py','scripts/lib/weakly_nonlinear_spatial_dynamics.py',
        'scripts/lib/weakly_nonlinear_planar_dynamics.py','scripts/lib/weakly_nonlinear_coupled_spatial_dynamics.py',
        'scripts/lib/nlsp_coupled_sensitivity_pilot.py','scripts/lib/yartsev_ch2_monoclinic_rod.py',
        'scripts/lib/nlsp_coupled_sensitivity_analysis.py','scripts/analysis/pilot_nlsp_coupled_planar_sensitivity.py',
        'scripts/lib/reddy_inplane_geometry.py','scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py',
        'scripts/lib/nlsp_fem3b_continuation.py']
    return {p:sha(ROOT/p) for p in files}

def finalize(b,item,summary):
    write(b/'summary.json',summary)
    write(b/'manifest.json',{'schema':item.get('schema',VERSION),'source_manifests':item['source_manifests'],
        'artifact_hashes':{p.relative_to(b).as_posix():sha(p) for p in sorted(b.rglob('*')) if p.is_file() and p!=b/'manifest.json'}})

def cache(b):
    b=Path(b)
    for p,h in read(b/'manifest.json')['artifact_hashes'].items():
        if sha(b/p)!=h:raise ValueError('Immutable pilot artifact changed: '+p)
    item=read(b/'provenance.json')
    for source in item['source_manifests'].values():
        if sha(ROOT/source['path']/'manifest.json')!=source['sha256']:raise ValueError('Historical source manifest changed')
    for p,h in item['source_artifacts'].items():
        if sha(ROOT/p)!=h:raise ValueError('Historical source artifact changed')
    return read(b/'summary.json')

def prepare(config):
    c=validate_config(read(config));model,selected,manifests=sources(c)
    item={'version':VERSION,'config':c,'source_artifacts':selected,'source_manifests':manifests,
        'HEAD':subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(),
        'python':sys.version,'initial_snapshot':'results/_smoke/nlsp_coupled_planar_sensitivity/initial_state.json'}
    fingerprint=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16];b=OUTPUT/fingerprint
    for existing in OUTPUT.glob('*/provenance.json') if OUTPUT.exists() else []:
        if read(existing)['config']['authorization']==AUTH and existing.parent!=b:raise ValueError('An attempt exists; code/HEAD changes do not authorize a retry')
    if (b/'manifest.json').exists():return b,read(b/'provenance.json'),cache(b),model
    if b.exists():raise ValueError('Interrupted unmanifested attempt; preserve evidence')
    b.mkdir(parents=True);item['initial_code_sha256']=code_hashes();write(b/'provenance.json',item);write(b/'config.json',c)
    s={'overall':'PREPARED','completed':False,'hard_stop':False,'stages':{k:'NOT_RUN' for k in 'ABCDEF'},
       'attempts':[],'calls':{'nonlinear_ODE':0,'tangent_ODE':0,'CCX':0,'Gmsh':0},'cases':{},'numerical_seconds':0.,
       'deferred_theory':c['deferred_theory'],'previous_spatial_qualification':'NUMERICAL_PARTIAL','execution_mode':c['execution_mode'],'admitted':False}
    finalize(b,item,s);return b,item,s,model

def physical_polynomial_map(disc,single,coordinates):
    """Exact split of one contained polynomial, with the existing reversal."""
    locals=[];sign=np.array([-1.,-1.,1.,-1.,-1.,1.,1.])
    for i,a in enumerate(disc.arms):
        x=a.x if i==0 else sum(disc.lengths)-a.x
        fields=single.reconstruct(coordinates,x)*(1 if i==0 else sign)
        locals.append(a.project(fields))
    raw=np.zeros(disc.ndof)
    for i,(a,q) in enumerate(zip(disc.arms,locals)):
        for f in range(7):raw[i*7*(disc.p-1)+f*(disc.p-1):i*7*(disc.p-1)+(f+1)*(disc.p-1)]=q[a.slices[FIELDS[f]]][:-1]
    raw[disc.joint_ids]=disc.T7[0]@disc.arms[0].reconstruct(locals[0],np.array([disc.lengths[0]]))[0]
    return disc.from_raw(raw)

def assembly_gates(c,model):
    result={'status':'PASS','geometry':{},'per_geometry':{},'straight_control':{}}
    first,_=rectangular_coefficients(.1,.05);swapped,_=rectangular_coefficients(.05,.1)
    if not (first.m==swapped.m and np.isclose(first.Bp,swapped.Bb,rtol=2e-14) and np.isclose(first.CT,swapped.CT,rtol=2e-10)):
        raise ArithmeticError('GEOMETRY_SWAP_NOT_EQUIVALENT')
    result['geometry']={'equal_area_mass':True,'bending_stiffnesses_exchanged':True,'weak_directions_perpendicular':True,
        'CT_swap_relative':abs(first.CT-swapped.CT)/first.CT,'coefficients_original':first.values(),'coefficients_swapped':swapped.values(),
        'reflection_preserved':True,'c_joint_closure_3D_validated':False}
    rng=np.random.default_rng(83)
    for name in ('A','B'):
        d=make_disc(c,name,8,model);q=rng.normal(size=d.ndof)*2e-7;v=rng.normal(size=d.ndof)*1e-7;x=np.r_[q,v]
        h=rng.normal(size=len(x));h/=np.linalg.norm(h);eps=1e-8
        J=d.jacobian(0,x);central=(d.rhs(0,x+eps*h)-d.rhs(0,x-eps*h))/(2*eps)
        jac=float(np.linalg.norm(J@h-central)/np.linalg.norm(central));plane=x.copy();plane[d.state_out]=0;PJ=d.jacobian(0,plane)
        off=max(np.max(abs(PJ[np.ix_(d.state_out,d.state_plane)])),np.max(abs(PJ[np.ix_(d.state_plane,d.state_out)])))
        joint=d.joint_diagnostics(q,v);rate=abs(d.energy_rate(q,v));parity=np.ones(len(x));parity[d.state_out]=-1
        mirror=float(np.linalg.norm(d.rhs(0,parity*x)-parity*d.rhs(0,x))/np.linalg.norm(d.rhs(0,x)))
        plrhs=float(np.max(abs(d.rhs(0,plane)[d.state_out])))
        lower=float(np.linalg.eigvalsh(d.M0).min());stiff=float(np.linalg.eigvalsh(d.K).min())
        if jac>2e-7 or off>2e-12 or plrhs>2e-12 or mirror>2e-12 or lower<=0 or stiff<=0 or rate>2e-12 or np.max(abs(joint['generalized_nodal_balance']))>2e-12:
            raise ArithmeticError('Required variational/Jacobian/reflection gate failed')
        # Independent linear density, including all four linear families.
        V=K=0.
        for a,qi,vi in zip(d.arms,d.local(q),d.local(v)):
            f=a.reconstruct(qi);s=a.reconstruct(qi,derivative=1);speed=a.reconstruct(vi);z=a.coefficients
            V+=float(a.weights@(z.C*s[:,0]**2+2*z.nu*z.C*s[:,0]*f[:,6]+z.C*f[:,6]**2+z.H*s[:,6]**2+
                z.S*((s[:,1]-f[:,5])**2+(s[:,2]-f[:,4])**2)+z.Bp*s[:,5]**2+z.Bb*s[:,4]**2+z.CT*s[:,3]**2))
            K+=float(a.weights@(z.m*np.sum(speed[:,:3]**2,axis=1)+(z.jp+z.jb)*speed[:,3]**2+z.jb*speed[:,4]**2+z.jp*(speed[:,5]**2+speed[:,6]**2)))
        linear_error=max(abs(V-q@d.K@q)/V,abs(K-v@d.M0@v)/K)
        if linear_error>2e-12:raise ArithmeticError('Independent linear densities disagree')
        result['per_geometry'][name]={'jacobian_directional_relative':jac,'off_diagonal_J_max':float(off),'planar_rhs_out_max':plrhs,
            'reflection_rhs_relative':mirror,'energy_rate_absolute':rate,'joint':joint,'M0_min_eigenvalue':lower,'K0_min_eigenvalue':stiff,'linear_density_relative':linear_error}
    d=CoupledSpatialGalerkin((first,first),16,(1.,1.),0.,model);single=SpatialGalerkin(first,32,length=2.,model=model)
    # Degree four is contained in both spaces; no new symbolic derivation.
    fields=lambda s: np.column_stack([s*(2-s)*(s-1)**(j%3)*(.0001 if j<3 else .001) for j in range(7)])
    qsingle=single.project(fields);vsingle=single.project(lambda s:fields(s)*.27)
    q=physical_polynomial_map(d,single,qsingle);v=physical_polynomial_map(d,single,vsingle)
    controls={'V':(d.potential(q,gradient=False)['V'],single.potential(qsingle,gradient=False)['V']),
              'T':(d.kinetic(q,v),single.kinetic(qsingle,vsingle)),
              'linear_weak_stiffness':(q@d.K@q,qsingle@single.K@qsingle),
              'linear_weak_mass':(v@d.M0@v,vsingle@single.M0@vsingle)}
    errors={key:abs(a-b)/max(abs(a),abs(b),1e-30) for key,(a,b) in controls.items()}
    if max(errors.values())>2e-12:raise ArithmeticError('Artificial straight joint action mismatch')
    frequencies={}
    for family,fs,nodes,block in [('mh',(0,6),(0,6),'mh'),('inplane',(1,5),(1,5),'timoshenko'),
        ('outplane',(2,4),(2,4),'outplane'),('torsion',(3,),(3,),'torsion')]:
        ids=np.r_[*[np.arange(i*7*15+f*15,i*7*15+(f+1)*15) for i in range(2) for f in fs],d.joint_ids[list(nodes)]]
        val=eigh(d.K[np.ix_(ids,ids)],d.M0[np.ix_(ids,ids)],subset_by_index=[0,1],check_finite=False)[0]
        old=single.linear_eigenpairs(block)['omega'][:2];new=np.sqrt(val);relative=np.max(abs(new-old)/old)
        frequencies[family]={'coupled':new,'single':old,'relative_max':float(relative)}
        if relative>1e-3:raise ArithmeticError('Artificial straight joint low modes unresolved at predeclared 0.1% control')
    result['straight_control']={'length':2.,'coupled_degree_per_arm':16,'single_degree':32,'action_relative_errors':errors,
        'frequencies':frequencies,'low_mode_control':1e-3,'same_external_value_BC':True,'slope_constraints':False}
    return result

def polynomial_w_max(d,q):
    candidates=[]
    for i,(a,qi) in enumerate(zip(d.arms,d.local(q))):
        coeff=qi[a.slices['w']];leg=np.zeros(a.p+1)
        for j,val in enumerate(coeff[:-1]):leg[j]+=val;leg[j+2]-=val
        leg[:2]+=coeff[-1]/2
        f=np.polynomial.legendre.Legendre(leg);roots=f.deriv().roots();x=np.r_[-1.,1.,roots.real[(abs(roots.imag)<1e-9)&(abs(roots.real)<=1)]]
        j=int(np.argmax(abs(f(x))));candidates.append({'arm':i+1,'s':float((x[j]+1)*a.length/2),'value':float(f(x[j]))})
    return max(abs(row['value']) for row in candidates),candidates

def frozen_excitation(c,name,model):
    p=max(c['degrees']);d=make_disc(c,name,p,model);unit=d.body_force(np.asarray(c['global_body_direction'],float))
    response=cho_solve(cho_factor(d.K,lower=True),unit);maximum,points=polynomial_w_max(d,response)
    target=c['linear_local_w_over_min_h']*min(h for b,h in c['sections'][name]);g=target/maximum
    omega,_=d.eigenpairs(plane=True);T=2*np.pi/omega[0]
    return {'reference_degree':p,'body_acceleration_magnitude':float(g),'global_direction':c['global_body_direction'],
        'target_maximum_local_w':target,'normalization':'maximum of represented linear local-w polynomial on either arm / min_i h_i; stationary points and endpoints',
        'unit_compliance_max':maximum,'maximum_candidates':points,'omega_planar_1':float(omega[0]),'Tref':float(T),
        'target_time_end':float(c['horizon_Tref']*T),'frozen_before_sensitivity':True,'no_load_fitting':True}

def settings(d,amplitude,c):
    raw=np.empty(d.ndof)
    for i,L in enumerate(d.lengths):
        for f in range(7):raw[i*7*(d.p-1)+f*(d.p-1):i*7*(d.p-1)+(f+1)*(d.p-1)]=amplitude if f<3 else amplitude/L
    raw[d.joint_ids]=[amplitude]*3+[amplitude/min(d.lengths)]*4
    inverse=np.linalg.inv(d.W);scale=np.linalg.norm(inverse*raw[None,:],axis=1)/np.sqrt(d.p)
    cutoff=max(np.sqrt(a.coefficients.C/a.coefficients.jp) for a in d.arms)
    return {'rtol':c['rtol'],'atol':c['atol_relative']*np.r_[scale,scale*cutoff],
        'max_step':2*np.pi/cutoff/24,'cutoff':float(cutoff),'normalization':'fixed raw field scales mapped by inverse resting-mass whitening, divided by sqrt(p); componentwise vector'}

def static_state(d,excitation):
    force=d.body_force(np.asarray(excitation['global_direction'])*excitation['body_acceleration_magnitude'])
    linear=cho_solve(cho_factor(d.K,lower=True),force)
    nl=fem2_static_newton(d,force,STATIC_POLICY)
    if nl['status']!='PASS':raise ArithmeticError('Nonlinear static equilibrium failed: '+str(nl.get('reason')))
    nonlinear=nl['coordinate'];d.safety(nonlinear)
    if np.max(abs(nonlinear[d.out]))>2e-12:raise ArithmeticError('Static plane leaked into outplane variables')
    return {'linear':linear,'nonlinear':nonlinear,'force':force,'newton':{k:v for k,v in nl.items() if k!='coordinate'},
        'energy_at_release':d.energy(nonlinear,np.zeros(d.ndof)),
        'relative_static_correction':float(np.linalg.norm(nonlinear-linear)/np.linalg.norm(linear)),
        'linear_maximum_local_w':polynomial_w_max(d,linear)[0],'nonlinear_maximum_local_w':polynomial_w_max(d,nonlinear)[0]}

def snapshot_case(case):
    folder=case/'execution_code';folder.mkdir(parents=True)
    hashes=code_hashes()
    for path,digest in hashes.items():shutil.copyfile(ROOT/path,folder/Path(path).name)
    write(case/'executed_source_hashes.json',hashes)

def integrate_attempt(b,item,s,label,d,initial,times,timing,deadline,plane=False,*,dense_archive=True,t_start=0.):
    """Save accepted prefix, attempt before call, and never retry a real run."""
    if label in s['cases']:raise ValueError('Attempt exists; no automatic nonlinear retry')
    if (b/label).exists():raise FileExistsError('Reserved attempt directory already exists')
    if times[0]!=t_start or times[-1]<=t_start:raise ValueError('Initial physical time mismatch')
    if s['calls']['nonlinear_ODE']>=item['config']['maximum_nonlinear_trajectories']:raise TimeoutError('Declared nonlinear call ceiling')
    case=b/label;case.mkdir(parents=True);snapshot_case(case)
    attempt={'label':label,'status':'STARTED','kind':'nonlinear_ODE','automatic_retry':False,'plane_invariant_subsystem':plane}
    s['attempts'].append(attempt);s['calls']['nonlinear_ODE']+=1;s['cases'][label]={'status':'STARTED'};finalize(b,item,s)
    write(case/'time_settings.json',timing);np.savez_compressed(case/'initial.npz',state=initial,times=times)
    ids=d.state_plane if plane else np.arange(2*d.ndof);history=np.empty((len(times),len(ids)));history[0]=initial[ids];cursor=1;records=[];steps=[];started=time.perf_counter();solver=None;failure=None
    d.reset_counters();limit=min(deadline,started+item['config']['per_nonlinear_run_seconds'])
    # Optional C finite-control storage: every native accepted state is streamed
    # without retaining/duplicating all cubic dense coefficients in RAM.
    # Historical callers retain the identical default dense archive and solve.
    native_times=[];native_stream=None
    if not dense_archive:
        native_stream=(case/'accepted_states.bin').open('wb')
        native_stream.write(np.asarray(history[0],dtype='<f8').tobytes())
    try:
        fun=d.plane_rhs if plane else d.rhs;jac=d.plane_jacobian if plane else d.jacobian
        def rhs(t,y):
            if time.perf_counter()>limit:raise TimeoutError('Declared integration budget')
            d.safety(d.embed_plane(y)[:d.ndof] if plane else y[:d.ndof]);return fun(t,y)
        solver=Radau(rhs,float(t_start),history[0],times[-1],jac=jac,rtol=timing['rtol'],atol=np.asarray(timing['atol'])[ids],max_step=timing['max_step'])
        while solver.status=='running':
            if time.perf_counter()>limit:raise TimeoutError('Declared integration budget')
            before=solver.t;message=solver.step()
            if solver.status=='failed':raise ArithmeticError(message or 'Radau failure')
            dense=solver.dense_output();steps.append(solver.t-before)
            native_times.append(float(solver.t))
            if dense_archive:
                records.append({'t_old':dense.t_old,'t':dense.t,'y_old':dense.y_old.copy(),'Q':dense.Q.copy()})
            else:
                native_stream.write(np.asarray(solver.y,dtype='<f8').tobytes())
            end=int(np.searchsorted(times,solver.t,side='right'))
            if end>cursor:history[cursor:end]=dense(times[cursor:end]).T;cursor=end
    except (ArithmeticError,ValueError,TimeoutError,np.linalg.LinAlgError) as error:failure=type(error).__name__+': '+str(error)
    finally:
        if native_stream is not None:native_stream.close()
    full=np.zeros((cursor,2*d.ndof));full[:,ids]=history[:cursor]
    energies=np.array([d.energy(x[:d.ndof],x[d.ndof:]) for x in full]);E0=energies[0]
    np.savez_compressed(case/'trajectory.npz',times=times[:cursor],state=full,energy=energies,energy_relative_drift=(energies-E0)/max(abs(E0),1e-30))
    if records:
        if t_start:save_dense_records(case/'dense.npz',records,t_start=float(t_start))
        else:save_dense_records(case/'dense.npz',records)
    if not dense_archive:
        actual_final=initial.copy()
        if solver is not None:actual_final[ids]=solver.y
        np.savez_compressed(case/'accepted_times.npz',times=np.r_[t_start,native_times],actual_final_state=actual_final)
        write(case/'accepted_state_layout.json',{'dtype':'<f8','shape':[len(native_times)+1,len(ids)],
            'state_ids':ids,'native_times':'accepted_times.npz','states':'accepted_states.bin',
            'all_accepted_coordinates_and_velocities_saved':True,'dense_coefficients_saved':False,
            'comparison_samples':'trajectory.npz','method_and_tolerances_unchanged':True})
    record={'status':'PASS' if failure is None and cursor==len(times) else 'PARTIAL','failure':failure,
        'target_end':float(times[-1]),'actual_accepted_end':float(solver.t if solver else 0),'sampled_end':float(times[cursor-1]),
        'samples':cursor,'accepted_steps':len(steps),'nfev':0 if solver is None else solver.nfev,'njev':0 if solver is None else solver.njev,'nlu':0 if solver is None else solver.nlu,
        'seconds':time.perf_counter()-started,'counters':d.counters(),'energy_relative_drift_max':float(np.max(abs((energies-E0)/max(abs(E0),1e-30)))),
        'plane':plane,'state_ids':ids,'zero_external_force_after_release':True,'automatic_retry':False,'safety_final':d.safety(full[-1,:d.ndof])}
    if t_start:
        record['physical_start']=float(t_start)
        final=initial.copy() if solver is None else d.embed_plane(solver.y) if plane else solver.y
        np.savez_compressed(case/'accepted_times.npz',times=np.r_[t_start,native_times],actual_final_state=final)
        record['actual_min_step']=float(min(steps)) if steps else None;record['actual_max_step']=float(max(steps)) if steps else None
    elif dense_archive:
        # Initial runs need the same exact physical restart checkpoint as
        # continued runs. This writer does not change the solver path.
        final=initial.copy() if solver is None else d.embed_plane(solver.y) if plane else solver.y.copy()
        np.savez_compressed(case/'accepted_times.npz',times=np.r_[t_start,native_times],actual_final_state=final)
        record['actual_min_step']=float(min(steps)) if steps else None;record['actual_max_step']=float(max(steps)) if steps else None
    write(case/'integration.json',record);attempt.update(status=record['status']);s['cases'][label]=record;s['numerical_seconds']+=record['seconds']
    if record['status']!='PASS':s.update(hard_stop=True,overall='NUMERICAL_FEASIBILITY_PARTIAL',stop_reason=failure or 'Incomplete prefix')
    finalize(b,item,s);return full,record


def make_times(Tref,end_fraction,c):
    count=int(round(end_fraction/c['sampling_Tref_fraction']))
    return np.linspace(0,end_fraction*Tref,count+1)

def verify_source_snapshot(case):
    for path,h in read(case/'executed_source_hashes.json').items():
        if sha(ROOT/path)!=h:raise ValueError('A scientific source changed after the saved attempt: '+path)

def tangent_validation(b,item,s,name,d,base_case,static,excitation,deadline):
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    c=item['config'];times=make_times(excitation['Tref'],c['tangent_check_Tref'],c)
    basis,metric,subspace=analysis.reference_subspace(d)
    tangent,record=analysis.tangent_attempt(b,item,s,name+'_short_tangent',d,base_case,times,basis,metric,deadline)
    if s['hard_stop']:return None
    zero,zero_record=analysis.tangent_attempt(b,item,s,name+'_zero_tangent',d,base_case,times,basis,metric,deadline,zero=True)
    if s['hard_stop']:return None
    gram=np.einsum('tij,ik,tkl->tjl',zero,metric,zero)
    autonomous_error=float(np.max(abs(gram-np.eye(basis.shape[1]))))
    if autonomous_error>1e-8:raise ArithmeticError('Autonomous reference energy norm not preserved by tangent integrator')
    state0=np.r_[static['nonlinear'],np.zeros(d.ndof)]
    base=evaluate_dense_records(base_case/'dense.npz',times)
    embedded=np.zeros((len(times),2*d.ndof));embedded[:,d.state_plane]=base
    norm=np.sqrt(2*static['energy_at_release']);delta0=c['perturbation_energy_norm_over_base']*norm
    timing=settings(d,excitation['target_maximum_local_w'],c);rows={};positive_mixed=None
    directions={'bending':np.eye(6)[0],'torsional':np.eye(6)[2],
                'mixed':np.array([1.,0.,1.,1.,0.,1.])/2}
    arrays=b/'check_arrays';arrays.mkdir(exist_ok=True)
    for shape,direction in directions.items():
        factors=c['decreasing_delta_factors'] if shape=='mixed' else c['decreasing_delta_factors'][:2]
        shape_rows=[];smallest=None
        for j,factor in enumerate(factors):
            delta=delta0*factor;initial=state0.copy();initial[d.state_out]+=delta*(basis@direction)
            full,run=integrate_attempt(b,item,s,f'{name}_check_{shape}_{j}',d,initial,times,timing,deadline)
            if s['hard_stop']:return None
            errors=analysis.full_tangent_error(d,full,embedded,tangent,direction,delta,metric)
            shape_rows.append({k:v for k,v in errors.items() if k not in ('out_error_curve','plane_change_curve','predicted','observed')})
            smallest=errors
            if shape=='mixed' and j==0:positive_mixed=full
        orders=analysis.measured_orders(shape_rows)
        # Numerical-floor cases are qualified, not falsely fitted to a power law.
        last=shape_rows[-1]['normalized_out_error_relative'];first=shape_rows[0]['normalized_out_error_relative']
        resolved=first>10*c['numerical_floor_relative']
        order_ok=all(row['normalized_out_order']>=1.5 for row in orders) if resolved and last>c['numerical_floor_relative'] else True
        if last>c['tangent_full_relative_gate'] or (last>c['numerical_floor_relative'] and last>first*1.05) or not order_ok:
            raise ArithmeticError('Tangent/full amplitude limit failed: '+name+' '+shape)
        rows[shape]={'rows':shape_rows,'orders':orders,'status':'PASS' if resolved and order_ok else 'PASS_WITH_NUMERICAL_FLOOR_QUALIFICATION',
            'asymptotic_order_measured':resolved and order_ok,'frozen_relative_gate':c['tangent_full_relative_gate']}
        pred=np.zeros((len(times),2*d.ndof));obs=pred.copy();pred[:,d.state_out]=smallest['predicted'];obs[:,d.state_out]=smallest['observed']
        pred_value=d.arms[0].reconstruct_series(pred[:,:d.ndof]@d.maps[0].T,np.array([d.lengths[0]/2]))[:,0,2]
        obs_value=d.arms[0].reconstruct_series(obs[:,:d.ndof]@d.maps[0].T,np.array([d.lengths[0]/2]))[:,0,2]
        np.savez_compressed(arrays/(name+'_'+shape+'.npz'),times=times,tau=times/excitation['Tref'],
            tangent_physical_v=pred_value,full_smallest_delta_physical_v=obs_value,
            out_error_curve=smallest['out_error_curve'],plane_change_curve=smallest['plane_change_curve'])
    negative=state0.copy();negative[d.state_out]-=delta0*(basis@directions['mixed'])
    minus,run=integrate_attempt(b,item,s,name+'_check_mixed_negative',d,negative,times,timing,deadline)
    if s['hard_stop']:return None
    out_parity=analysis.energy_norm((positive_mixed[:,d.state_out]+minus[:,d.state_out])/delta0,metric).max()
    pl_parity=np.max(abs(positive_mixed[:,d.state_plane]-minus[:,d.state_plane]))
    if out_parity>1e-7 or pl_parity>1e-10:raise ArithmeticError('Full nonlinear reflection parity failed')
    return {'status':'PASS_WITH_QUALIFICATIONS','subspace':subspace,'delta0':float(delta0),
        'delta0_over_base_reference_energy_norm':c['perturbation_energy_norm_over_base'],
        'shapes':rows,'autonomous_norm_gram_error':autonomous_error,'full_odd_parity_relative':float(out_parity),
        'full_even_planar_parity_absolute':float(pl_parity),'coverage_Tref':c['tangent_check_Tref']}


PLANNING_STOP='Preflight forecast exceeds declared budget; 2Tref production not authorized to shorten automatically'

def reassess_planning_stop(b,item,s):
    """Correct only a pre-production cost forecast, never a scientific failure.

    The successful cost probe is reused. Its actual attempt, first forecast,
    first stop and original manifest remain evidence of the planning history.
    No horizon, numerical gate, tolerance, budget or attempt cap is changed.
    """
    correction=b/'cost_preflight_reassessment.json'
    if not s['hard_stop'] or s.get('stop_reason')!=PLANNING_STOP or not correction.exists():return False
    if s['completed'] or s['calls']!={'nonlinear_ODE':1,'tangent_ODE':0,'CCX':0,'Gmsh':0}:return False
    if set(s['cases'])!={'cost_probe'} or len(s['attempts'])!=1:return False
    if s['attempts'][0].get('label')!='cost_probe' or s['attempts'][0].get('status')!='PASS':return False
    if s['cases']['cost_probe']['status']!='PASS' or s['cases']['cost_probe'].get('failure'):return False
    if any(s['stages'][k]!='NOT_RUN' for k in 'CDEF'):return False
    c=validate_config(item['config']);r=read(correction)
    if r['original_cost_estimate_sha256']!=sha(b/'cost_preflight.json') or r['saved_probe_sha256']!=sha(b/'cost_probe/integration.json'):return False
    if r['new_ODE_calls']!=0 or r['budget_seconds']!=c['numerical_budget_seconds']:return False
    if not r['physical_horizons_preserved'] or r['scientific_gates_or_tolerances_changed']:return False
    if r['decision']!='PRODUCTION_WITHIN_ORIGINAL_BUDGET' or not 0<r['total_planning_seconds']<=r['budget_seconds']:return False
    archive=b/'first_cost_decision';archive.mkdir()
    for file in ('summary.json','manifest.json','config.json','cost_preflight.json'):
        shutil.copyfile(b/file,archive/file)
    write(b/'planning_decision_ledger.json',[
        {'decision':'STOP_BEFORE_PRODUCTION','forecast_seconds':read(b/'cost_preflight.json')['total_planning_seconds'],
         'reason':PLANNING_STOP,'preserved_in':'first_cost_decision'},
        {'decision':r['decision'],'reassessment_sha256':sha(correction),'new_scientific_calls':0,
         'reason':'Use each geometry physical horizon and separate measured RHS/Jacobian/LU costs with saved probe counters; remove guessed combined multipliers',
         'budget_horizon_tolerances_caps_unchanged':True}])
    s.update(hard_stop=False,overall='PREPARED',cost_preflight='PASS',cost_reassessment_sha256=sha(correction))
    s.pop('stop_reason',None);finalize(b,item,s);return True

def run_program(b,item,s,model,through='all'):
    reassess_planning_stop(b,item,s)
    if s['completed'] or s['hard_stop']:return s
    c=item['config'];start=time.perf_counter();deadline=start+max(0,c['numerical_budget_seconds']-s['numerical_seconds'])
    if s['stages']['A']=='NOT_RUN':
        try:
            gates=assembly_gates(c,model);write(b/'assembly_gates.json',gates);s['stages']['A']='PASS';s['stages']['B']='PASS'
        except (ArithmeticError,ValueError,np.linalg.LinAlgError) as error:
            s.update(overall='BLOCKED_BY_ASSEMBLY',hard_stop=True,stop_reason=str(error));write(b/'assembly_gates.json',{'status':'FAIL','reason':str(error)});finalize(b,item,s);return s
        excitation={name:frozen_excitation(c,name,model) for name in ('A','B')};write(b/'frozen_excitation.json',excitation)
        finalize(b,item,s)
    if through=='preflight':return s
    excitation=read(b/'frozen_excitation.json')
    if 'cost_probe' not in s['cases']:
        d=make_disc(c,'A',c['cost_probe_degree'],model);static=static_state(d,excitation['A']);write(b/'cost_probe_static.json',{k:v for k,v in static.items() if k not in ('linear','nonlinear','force')})
        timing=settings(d,excitation['A']['target_maximum_local_w'],c);times=make_times(excitation['A']['Tref'],c['cost_probe_Tref'],c)
        state=np.r_[static['nonlinear'],np.zeros(d.ndof)]
        full,probe=integrate_attempt(b,item,s,'cost_probe',d,state,times,timing,deadline,plane=True)
        if s['hard_stop']:return s
        # Measure both target matrix kernels; do not use a universal p^3 guess.
        kernels={}
        for name in ('A','B'):
            for degree in c['degrees']:
                disc=make_disc(c,name,degree,model);q=cho_solve(cho_factor(disc.K,lower=True),disc.body_force(np.array([0.,-1.,0.])))*excitation[name]['body_acceleration_magnitude'];X=np.r_[q,np.zeros(disc.ndof)]
                t0=time.perf_counter()
                for j in range(8):disc.rhs(0,X);disc.jacobian(0,X)
                kernels[f'{name}_p{degree}']=(time.perf_counter()-t0)/8
        q0=np.r_[static['nonlinear'],np.zeros(d.ndof)];t0=time.perf_counter()
        for j in range(8):d.rhs(0,q0);d.jacobian(0,q0)
        reference=(time.perf_counter()-t0)/8
        # Production: four bases, two full finite checks, sixteen short checks,
        # six nonzero-base and two autonomous tangent integrations. This is a
        # planning estimate; step counts and LU reuse may change on longer runs.
        base_est={key:probe['seconds']/c['cost_probe_Tref']*c['horizon_Tref']*value/reference for key,value in kernels.items()}
        estimate=sum(base_est.values())+2.0*(base_est['A_p24']+base_est['B_p24'])+16*c['tangent_check_Tref']/c['horizon_Tref']*max(base_est['A_p16'],base_est['B_p16'])+1.5*sum(base_est.values())
        plan={'probe_seconds':probe['seconds'],'probe_accepted_steps':probe['accepted_steps'],'target_kernel_seconds':kernels,
            'probe_kernel_seconds':reference,'per_base_planning_seconds':base_est,'total_planning_seconds':estimate,
            'budget_seconds':c['numerical_budget_seconds'],'target_horizon_Tref':2.,'horizon_not_changed':True,
            'estimate_is_not_guarantee':True,'nonlinear_attempt_ceiling':c['maximum_nonlinear_trajectories'],
            'staged_gate':'no production if the estimate exceeds the finite budget'}
        write(b/'cost_preflight.json',plan)
        if estimate>c['numerical_budget_seconds']:
            s.update(overall='NUMERICAL_FEASIBILITY_PARTIAL',hard_stop=True,stop_reason=PLANNING_STOP)
            finalize(b,item,s);return s
        s['cost_preflight']='PASS';finalize(b,item,s)
    if through=='cost':return s
    from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis
    checks=read(b/'tangent_checks.json') if (b/'tangent_checks.json').exists() else {}
    gains=read(b/'gain_results.json') if (b/'gain_results.json').exists() else {}
    for name in ('A','B'):
        d=make_disc(c,name,16,model);label=name+'_p16_base';case=b/label
        if label not in s['cases']:
            static=static_state(d,excitation[name]);np.savez_compressed(b/(name+'_p16_static.npz'),**{k:static[k] for k in ('linear','nonlinear','force')});write(b/(name+'_p16_static.json'),{k:v for k,v in static.items() if k not in ('linear','nonlinear','force')})
            state=np.r_[static['nonlinear'],np.zeros(d.ndof)];timing=settings(d,excitation[name]['target_maximum_local_w'],c)
            times=make_times(excitation[name]['Tref'],2.,c);full,record=integrate_attempt(b,item,s,label,d,state,times,timing,deadline,plane=True)
            if s['hard_stop']:return s
            analysis.save_base_observations(case,d,excitation[name]['Tref'])
        else:
            verify_source_snapshot(case)
            with np.load(b/(name+'_p16_static.npz')) as z:static={key:z[key] for key in z.files}
            static.update(read(b/(name+'_p16_static.json')))
        if name not in checks:
            try:checks[name]=tangent_validation(b,item,s,name,d,case,static,excitation[name],deadline)
            except (ArithmeticError,ValueError,np.linalg.LinAlgError) as error:
                s['stages']['C']='PARTIAL';s['stages']['D']='FAIL'
                s.update(overall='BLOCKED_BY_TANGENT_CHECK',hard_stop=True,stop_reason=str(error));write(b/'tangent_check_failure.json',{'geometry':name,'reason':str(error)});finalize(b,item,s);return s
            if s['hard_stop']:return s
            write(b/'tangent_checks.json',checks);finalize(b,item,s)
    s['stages']['C']='PASS';s['stages']['D']='PASS_WITH_QUALIFICATIONS';finalize(b,item,s)
    for name in ('A','B'):
        gains.setdefault(name,{})
        for degree in c['degrees']:
            d=make_disc(c,name,degree,model);label=f'{name}_p{degree}_base';case=b/label
            if label not in s['cases']:
                static=static_state(d,excitation[name]);np.savez_compressed(b/(name+f'_p{degree}_static.npz'),**{k:static[k] for k in ('linear','nonlinear','force')});write(b/(name+f'_p{degree}_static.json'),{k:v for k,v in static.items() if k not in ('linear','nonlinear','force')})
                timing=settings(d,excitation[name]['target_maximum_local_w'],c);state=np.r_[static['nonlinear'],np.zeros(d.ndof)];times=make_times(excitation[name]['Tref'],2.,c)
                full,run=integrate_attempt(b,item,s,label,d,state,times,timing,deadline,plane=True)
                if s['hard_stop']:return s
                analysis.save_base_observations(case,d,excitation[name]['Tref'])
            basis,metric,subspace=analysis.reference_subspace(d);times=make_times(excitation[name]['Tref'],2.,c)
            tangent_label=f'{name}_p{degree}_tangent'
            if tangent_label not in s['cases']:
                history,run=analysis.tangent_attempt(b,item,s,tangent_label,d,case,times,basis,metric,deadline)
                if s['hard_stop']:return s
            else:
                with np.load(b/tangent_label/'tangent.npz') as z:history=z['vectors']
            result,curve=analysis.gain_result(d,history,times,basis,metric);result['subspace']=subspace
            gains[name][str(degree)]=result;np.savez_compressed(b/f'{name}_p{degree}_gain.npz',times=times,tau=times/excitation[name]['Tref'],gain=curve)
            write(b/'gain_results.json',gains);finalize(b,item,s)
        low,high=gains[name]['16'],gains[name]['24'];change=abs(high['gain']-low['gain'])/high['gain']
        gains[name]['gain_spatial_relative_change']=change;gains[name]['gain_spatial_status']='PASS' if change<=c['gain_spatial_relative_gate'] else 'PARTIAL'
        d=make_disc(c,name,24,model);basis,metric,subspace=analysis.reference_subspace(d)
        with np.load(b/(name+'_p24_static.npz')) as z:initial=np.r_[z['nonlinear'],np.zeros(d.ndof)]
        static=read(b/(name+'_p24_static.json'));delta=checks[name]['delta0']/4
        direction=np.asarray(high['initial_combination']);initial[d.state_out]+=delta*(basis@direction)
        times=make_times(excitation[name]['Tref'],2.,c)
        label=name+'_optimal_finite'
        if label not in s['cases']:
            full,run=integrate_attempt(b,item,s,label,d,initial,times,settings(d,excitation[name]['target_maximum_local_w'],c),deadline)
            if s['hard_stop']:return s
        else:
            with np.load(b/label/'trajectory.npz') as z:full=z['state']
        with np.load(b/(name+'_p24_tangent')/'tangent.npz') as z:tangent=z['vectors']
        base=evaluate_dense_records(b/(name+'_p24_base')/'dense.npz',times);embedded=np.zeros_like(full);embedded[:,d.state_plane]=base
        error=analysis.full_tangent_error(d,full,embedded,tangent,direction,delta,metric)
        gains[name]['finite_check']={k:v for k,v in error.items() if k not in ('out_error_curve','plane_change_curve','predicted','observed')}
        gains[name]['finite_check']['status']='PASS' if error['normalized_out_error_relative']<=c['tangent_full_relative_gate'] else 'PARTIAL'
        predicted=analysis.energy_norm(error['predicted'],metric);actual=analysis.energy_norm(error['observed'],metric)
        np.savez_compressed(b/(name+'_finite_check.npz'),times=times,tau=times/excitation[name]['Tref'],tangent_norm=predicted,full_norm=actual,error=error['out_error_curve'],planar_change=error['plane_change_curve'])
        seeds=np.zeros((2,2*d.ndof));seeds[0,d.state_out]=basis@direction;seeds[1,d.state_out]=np.asarray(high['response_at_max'])
        profiles,obs=analysis.physical_observations(d,seeds)
        np.savez_compressed(b/(name+'_optimal_profiles.npz'),x=np.linspace(0,1,201),initial=profiles[0],response=profiles[1],delta=delta,time_at_max=high['time_at_max'])
        write(b/'gain_results.json',gains);finalize(b,item,s)
    s['stages']['E']='COMPLETE_WITH_QUALIFICATIONS';s['stages']['F']='PASS' if all(gains[n]['finite_check']['status']=='PASS' for n in ('A','B')) else 'PARTIAL'
    s['overall']='PILOT_METHOD_CONFIRMED' if all(gains[n]['gain_spatial_status']=='PASS' and gains[n]['finite_check']['status']=='PASS' for n in ('A','B')) else 'NUMERICAL_PARTIAL'
    s['completed']=True;s['scientific_claim']='restricted initial-state sensitivity, not Floquet stability, imperfection sensitivity, or 3D validation'
    s['figures']=analysis.plot_saved(b);finalize(b,item,s);return s


def main(argv=None):
    argv=list(sys.argv[1:] if argv is None else argv)
    if '--pre-supervisor-checks' in argv:
        from scripts.lib.nlsp_planar_1to3_internal_resonance import supervisor_main
        return supervisor_main([arg for arg in argv if arg!='--pre-supervisor-checks'])
    if '--beta15-p14-convergence' in argv:
        from scripts.lib.nlsp_planar_1to3_internal_resonance import p14_main
        return p14_main([arg for arg in argv if arg!='--beta15-p14-convergence'])
    if '--beta15-dynamic-out-growth' in argv:
        from scripts.lib.nlsp_planar_1to3_internal_resonance import dynamic_main
        return dynamic_main([arg for arg in argv if arg!='--beta15-dynamic-out-growth'])
    if '--beta15-periodic-spatial-stability' in argv:
        from scripts.lib.nlsp_planar_1to3_internal_resonance import stability_main
        return stability_main([arg for arg in argv if arg!='--beta15-periodic-spatial-stability'])
    if '--planar-1to3-internal-resonance' in argv:
        from scripts.lib.nlsp_planar_1to3_internal_resonance import main as internal_main
        return internal_main([arg for arg in argv if arg!='--planar-1to3-internal-resonance'])
    if '--beta15-nonlinear-planar-resonance' in argv:
        from scripts.lib.nlsp_beta15_nonlinear_planar_resonance import main as resonance_main
        return resonance_main([arg for arg in argv if arg!='--beta15-nonlinear-planar-resonance'])
    if '--instability-feasibility' in argv:
        from scripts.lib.nlsp_instability_feasibility import main as feasibility_main
        return feasibility_main([arg for arg in argv if arg!='--instability-feasibility'])
    if '--beta-resonant-screening' in argv:
        from scripts.lib.nlsp_resonant_spatial_coupling_screening import beta_main
        return beta_main([arg for arg in argv if arg!='--beta-resonant-screening'])
    if '--resonant-screening' in argv:
        from scripts.lib.nlsp_resonant_spatial_coupling_screening import main as screening_main
        return screening_main([arg for arg in argv if arg!='--resonant-screening'])
    if '--construction-c-longer' in argv:
        from scripts.lib.nlsp_construction_c_longer_sensitivity import main as longer_main
        return longer_main([arg for arg in argv if arg!='--construction-c-longer'])
    if '--construction-c' in argv:
        from scripts.lib.nlsp_construction_c_sensitivity import main as construction_main
        return construction_main([arg for arg in argv if arg!='--construction-c'])
    if '--resume-sensitivity' in argv:
        from scripts.lib.nlsp_coupled_sensitivity_resume import main as resume_main
        return resume_main([arg for arg in argv if arg!='--resume-sensitivity'])
    if '--torsion-diagnostic' in argv:
        from scripts.lib.nlsp_coupled_torsion_diagnostic import main as diagnostic_main
        return diagnostic_main([arg for arg in argv if arg!='--torsion-diagnostic'])
    parser=argparse.ArgumentParser(description=__doc__);mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--compute',action='store_true');mode.add_argument('--report-only',type=Path);mode.add_argument('--plot-only',type=Path)
    parser.add_argument('--config',type=Path,default=CONFIG);parser.add_argument('--through-stage',choices=('preflight','cost','all'),default='all')
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only:
        b=(args.report_only or args.plot_only).resolve();s=cache(b)
        if args.plot_only and not (b/'figure_provenance.json').exists():
            from scripts.lib.nlsp_coupled_sensitivity_analysis import plot_saved
            s['figures']=plot_saved(b);finalize(b,read(b/'provenance.json'),s)
    else:
        b,item,s,model=prepare(args.config);s=run_program(b,item,s,model,args.through_stage)
    print(json.dumps(plain({'bundle':str(b),'overall':s['overall'],'stages':s['stages'],'calls':s['calls'],'stop_reason':s.get('stop_reason')}),ensure_ascii=False,indent=2))
