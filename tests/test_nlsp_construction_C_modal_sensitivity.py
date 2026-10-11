"""C preset contracts and saved-state algebra; no real IVP or FEM in tests."""
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.lib import nlsp_construction_c_sensitivity as c
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_coupled_sensitivity_analysis as a


@pytest.fixture(scope='module')
def disc():
    config=p.read(c.CONFIG);model,_,_=p.sources(p.read(p.CONFIG))
    return p.make_disc(config,'C',16,model)


def test_frozen_authorization_and_call_budget():
    config=c.validate_scope(p.read(c.CONFIG))
    assert config['authorization'] not in (p.AUTH,'explicit_user_coupled_sensitivity_resume_2026_10_10')
    assert config['maximum_nonlinear_trajectories']==3 and config['maximum_planar_trajectories']==2
    assert config['maximum_tangent_trajectories']==8 and config['numerical_budget_seconds']==7200


@pytest.mark.parametrize('key,value',[('amplitude',.01),('degree',24),('new_FEM_jobs',1),('automatic_retry',True)])
def test_cannot_change_scope(key,value):
    config=p.read(c.CONFIG);config[key]=value
    with pytest.raises(ValueError):c.validate_scope(config)


def test_geometry_axes_and_stiffness(disc):
    assert disc.lengths==(.8,1.2) and disc.beta_deg==30
    for arm,B in zip(disc.arms,disc.bases):
        assert np.isclose(arm.coefficients.Bp/arm.coefficients.Bb,4,rtol=1e-14)
        assert np.allclose(B.T@B,np.eye(3),atol=1e-15) and np.linalg.det(B)>0
        assert np.array_equal(B[:,2],[0,0,-1])


def test_linear_density_decomposition(disc):
    q=np.zeros(disc.ndof);q[disc.plane]=np.sin(np.arange(len(disc.plane)))*1e-6
    bending=c.planar_bending_parts(disc,q).sum();axial=0
    for arm,qi in zip(disc.arms,disc.local(q)):
        f=arm.reconstruct(qi);ds=arm.reconstruct(qi,derivative=1);co=arm.coefficients
        axial+=arm.weights@(co.C*(ds[:,0]**2+f[:,6]**2+2*co.nu*f[:,6]*ds[:,0])+co.H*ds[:,6]**2)
    assert np.isclose(bending+axial,q@disc.K@q,rtol=2e-12)


def test_mode_selection_scaling_and_joint(disc):
    q,report=c.select_modal_state(disc,p.read(c.CONFIG))
    assert report['selected_sorted_index_one_based']==1
    assert report['modes'][0]['bending_shear_fraction']>.5
    assert np.isclose(p.polynomial_w_max(disc,q)[0],.005,rtol=3e-15)
    assert np.all(q[disc.out]==0) and report['zero_initial_velocity']
    assert report['eigen_residual_relative']<1e-8
    joint=disc.joint_diagnostics(q,np.zeros(disc.ndof))
    assert joint['translation_difference']<2e-12 and joint['orientation_difference']<2e-12
    for arm,qi in zip(disc.arms,disc.local(q)):
        assert np.max(abs(arm.reconstruct(qi,np.array([0.]))))<1e-17


def test_displacement_only_orthonormal_seeds(disc):
    B,E,report=c.displacement_subspace(disc,p.read(c.CONFIG),.31338933639565253)
    assert B.shape==(186,3) and np.all(B[len(disc.out):]==0)
    assert np.max(abs(B.T@E@B-np.eye(3)))<2e-10
    assert np.linalg.eigvalsh(E).min()>0
    assert report['selected_sorted_indices_one_based']==[1,2,5]


def test_physical_perturbation_scale_and_zero_velocities(disc):
    B,E,_=c.displacement_subspace(disc,p.read(c.CONFIG),.31338933639565253)
    initial=B@np.array([.6,0,.8]);delta,report=c.perturbation_scale(disc,initial,p.read(c.CONFIG))
    assert np.isclose(np.max(report['actual_initial_v_LPhi_Lpsi']),5e-5,rtol=2e-15)
    assert report['all_initial_velocities_zero'] and delta>0
    bad=initial.copy();bad[-1]=1
    with pytest.raises(ValueError):c.perturbation_scale(disc,bad,p.read(c.CONFIG))


def test_zero_base_exact_modal_evolution_energy(disc):
    B,E,_=c.displacement_subspace(disc,p.read(c.CONFIG),.31338933639565253)
    n=len(disc.out);times=np.linspace(0,40,11);omega=np.sqrt(np.einsum('ik,ij,jk->k',B[:n],E[:n,:n],B[:n])/np.einsum('ik,ij,jk->k',B[:n],E[n:,n:],B[:n]))
    Z=np.zeros((len(times),2*n,3));Z[:,:n]=B[:n][None]*np.cos(times[:,None]*omega)[:,None]
    Z[:,n:]=-B[:n][None]*(omega*np.sin(times[:,None]*omega))[:,None]
    gram=np.einsum('tij,ik,tkl->tjl',Z,E,Z)
    assert np.max(abs(gram-np.eye(3)))<2e-10


def test_same_tangent_operator_from_full_jacobian(disc):
    q,_=c.select_modal_state(disc,p.read(c.CONFIG));X=np.r_[q,np.zeros(disc.ndof)]
    J=disc.jacobian(0,X)
    assert np.max(abs(J[np.ix_(disc.state_out,disc.state_plane)]))==0
    assert np.max(abs(J[np.ix_(disc.state_plane,disc.state_out)]))==0


def test_C_directional_rhs_agrees_with_tangent_without_ODE(disc):
    q,_=c.select_modal_state(disc,p.read(c.CONFIG));X=np.r_[q,np.zeros(disc.ndof)]
    B,_,_=c.displacement_subspace(disc,p.read(c.CONFIG),.31338933639565253)
    seed=np.zeros(2*disc.ndof);seed[disc.state_out]=B[:,2]
    eps=1e-9;central=(disc.rhs(0,X+eps*seed)-disc.rhs(0,X-eps*seed))/(2*eps)
    exact=disc.jacobian(0,X)@seed
    assert np.linalg.norm(central-exact)/np.linalg.norm(exact)<2e-7


def test_gain_keeps_all_output_components():
    d=SimpleNamespace(out=np.arange(2),ndof=2,reference_energy_parts=lambda q,v:np.ones((2,2)))
    B=np.eye(4)[:,:1];Z=np.stack([B,B.copy()]);Z[1,3,0]=2
    report,G=a.gain_result(d,Z,np.array([0,1]),B,np.eye(4))
    assert np.isclose(G[-1],np.sqrt(5)) and report['full_outplane_response_dimension']==4


def test_no_retry_and_ceilings():
    config=p.read(c.CONFIG);s={'hard_stop':True,'numerical_seconds':0,'calls':{'nonlinear_ODE':0,'tangent_ODE':0},'cases':{}}
    with pytest.raises(RuntimeError):c.limits(s,config,'tangent')
    s['hard_stop']=False;s['calls']['tangent_ODE']=8
    with pytest.raises(TimeoutError):c.limits(s,config,'tangent')


def test_incomplete_attempt_never_launches_next_case(monkeypatch):
    summary={'completed':False,'hard_stop':False,'cases':{'C_p16_base':{'status':'STARTED'}}}
    monkeypatch.setattr(p,'sources',lambda *x:pytest.fail('new action/model reconstruction'))
    with pytest.raises(RuntimeError):c.run(None,None,summary)


def test_progress_observation_preserves_evaluation(tmp_path):
    calls=[]
    fn=lambda t,state:calls.append(t) or state*2
    d=SimpleNamespace(rhs=fn,jacobian=fn)
    with c.operator_progress(d,tmp_path,'test',1):
        assert np.array_equal(d.rhs(.25,np.array([3.])),[6.])
    assert d.rhs is fn and calls==[.25]
    assert 'not an accepted' in p.read(tmp_path/'test_progress.json')['qualification']


def test_streamed_native_states_preserve_samples_and_legacy_default(tmp_path,monkeypatch):
    import time
    class MockRadau:
        def __init__(self,fun,t0,initial,end,**kw):
            self.t=t0;self.y=initial.copy();self.status='running';self.nfev=self.njev=self.nlu=0
        def step(self):
            self.old=self.t;self.t+=.25
            if self.t==1:self.status='finished'
        def dense_output(self):
            class Dense:
                t_old=self.old;t=self.t;y_old=self.y.copy();Q=np.zeros((4,3))
                def __call__(this,times):return np.tile(this.y_old[:,None],(1,len(times)))
            return Dense()
    d=SimpleNamespace(ndof=2,rhs=lambda *x:None,jacobian=lambda *x:None,reset_counters=lambda:None,
        counters=lambda:{},energy=lambda q,v:float(q@q+v@v),safety=lambda q:{})
    monkeypatch.setattr(p,'Radau',MockRadau);monkeypatch.setattr(p,'snapshot_case',lambda case:None)
    monkeypatch.setattr(p,'finalize',lambda *x:None)
    outputs=[]
    for label,dense in [('legacy',True),('streamed',False)]:
        s={'cases':{},'calls':{'nonlinear_ODE':0},'attempts':[],'numerical_seconds':0.}
        item={'config':{'maximum_nonlinear_trajectories':1,'per_nonlinear_run_seconds':10}}
        state,record=p.integrate_attempt(tmp_path,item,s,label,d,np.arange(4,dtype=float),np.array([0.,.5,1.]),
            {'rtol':1e-10,'atol':np.ones(4),'max_step':.1},time.perf_counter()+10,dense_archive=dense)
        outputs.append(state);assert record['status']=='PASS'
    assert np.array_equal(*outputs) and (tmp_path/'legacy/dense.npz').exists()
    assert not (tmp_path/'streamed/dense.npz').exists()
    layout=p.read(tmp_path/'streamed/accepted_state_layout.json')
    native=np.memmap(tmp_path/'streamed/accepted_states.bin',dtype=layout['dtype'],shape=tuple(layout['shape']),mode='r')
    assert np.array_equal(native,np.tile(np.arange(4),(5,1)))
    with np.load(tmp_path/'streamed/accepted_times.npz') as z:
        assert np.array_equal(z['times'],[0,.25,.5,.75,1])
        assert np.array_equal(z['actual_final_state'],np.arange(4))


def test_source_parent_and_kernels_preserved():
    parent=c.verify_sources(p.read(c.CONFIG))
    assert p.sha(parent/'manifest.json')==p.read(c.CONFIG)['parent_manifest_sha256']


def test_cached_run_avoids_all_numerical_calls(monkeypatch):
    summary={'completed':True,'hard_stop':False}
    monkeypatch.setattr(p,'sources',lambda *x:pytest.fail('source reconstruction on cache'))
    monkeypatch.setattr(p,'integrate_attempt',lambda *x:pytest.fail('new full IVP'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x,**kw:pytest.fail('new tangent'))
    assert c.run(None,None,summary)==summary


def test_dispatch_keeps_old_cli_routes(monkeypatch):
    seen=[];monkeypatch.setattr(c,'main',lambda args:seen.append(args))
    p.main(['--construction-c','--compute','--through-stage','preflight'])
    assert seen==[['--compute','--through-stage','preflight']]


def test_fixed_signal_control_arithmetic(tmp_path,monkeypatch):
    p.write(tmp_path/'C_p16_gain.json',{'gain':1.02})
    p.write(tmp_path/'spatial_sensitivity.json',{'gain_curve_change_max':.003})
    p.write(tmp_path/'temporal_sensitivity.json',{'E0_response_change_max':1e-5})
    result=c.qualify_saved(tmp_path,p.read(c.CONFIG))
    assert not result['gain_resolved_in_observed_controls']
    assert np.isclose(result['signal_over_spatial_change'],.02/.003)
    p.write(tmp_path/'spatial_sensitivity.json',{'gain_curve_change_max':.0001})
    assert c.qualify_saved(tmp_path,p.read(c.CONFIG))['gain_resolved_in_observed_controls']


@pytest.mark.parametrize('mode',['--report-only','--plot-only'])
def test_cache_cli_uses_no_solvers(tmp_path,monkeypatch,mode):
    p.write(tmp_path/'export_provenance.json',{'scientific_calls':0})
    summary={'completed':True,'overall':'C_MODAL_SENSITIVITY_NUMERICAL_PARTIAL','stages':{},'calls':{},'numerical_seconds':123.}
    monkeypatch.setattr(c,'cache',lambda b:summary)
    monkeypatch.setattr(p,'make_disc',lambda *x:pytest.fail('model on cache'))
    monkeypatch.setattr(p,'sources',lambda *x:pytest.fail('action on cache'))
    monkeypatch.setattr(p,'static_state',lambda *x:pytest.fail('static on cache'))
    monkeypatch.setattr(p,'integrate_attempt',lambda *x:pytest.fail('full on cache'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x,**kw:pytest.fail('tangent on cache'))
    c.main([mode,str(tmp_path)])
