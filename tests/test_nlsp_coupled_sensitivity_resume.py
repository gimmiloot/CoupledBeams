"""Continuation contracts and algebra; no ODE or native FEM in tests."""
import inspect
import numpy as np
import pytest
from scripts.lib import nlsp_coupled_sensitivity_resume as r
from scripts.lib import nlsp_coupled_sensitivity_pilot as p
from scripts.lib import nlsp_coupled_sensitivity_analysis as a

def test_scope_and_independent_authorization():
    c=r.validate_scope(p.read(r.CONFIG))
    assert c['authorization'] not in (p.AUTH,'explicit_user_torsion_tangent_diagnostic_2026_10_10')
    assert c['maximum_nonlinear_trajectories']==8 and c['maximum_tangent_trajectories']==16
    assert c['numerical_budget_seconds']==7200 and c['new_FEM_jobs']==0
    assert c['horizon_Tref']==2 and c['degree']==16 and c['lower_degree']==12

@pytest.mark.parametrize('key,value',[('degree',24),('authorization','code_hash_changed'),('new_FEM_jobs',1),('maximum_nonlinear_trajectories',9)])
def test_scope_cannot_expand(key,value):
    c=p.read(r.CONFIG);c[key]=value
    with pytest.raises(ValueError):r.validate_scope(c)

def test_mixed_input_has_bending_and_torsion_and_zero_initial_speed():
    seed=r.mixed_direction()
    assert np.isclose(seed@seed,1) and np.all(seed[:2]!=0) and np.all(seed[2:]==0)
    assert p.read(r.CONFIG)['saved_seed_columns']==[0,2,3,5]

def test_parent_source_preservation_and_frozen_kernel():
    paths=r.verify_sources(p.read(r.CONFIG))
    assert p.read(paths['parent']/'summary.json')['overall']=='BLOCKED_BY_TANGENT_CHECK'
    assert p.read(paths['diagnostic']/'summary.json')['overall']=='NUMERICAL_FLOOR_IDENTIFIED'

@pytest.fixture(scope='module')
def discs():
    c=p.read(p.CONFIG);model,_,_=p.sources(c)
    return p.make_disc(c,'B',16,model),p.make_disc(c,'B',12,model)

def test_fixed_positive_energy_metric(discs):
    d,_=discs;metric=r.fixed_metric(d);n=len(d.out)
    assert np.linalg.eigvalsh(metric).min()>0
    assert np.array_equal(metric[:n,:n],d.K[np.ix_(d.out,d.out)])
    assert np.array_equal(metric[n:,n:],d.M0[np.ix_(d.out,d.out)])

def test_geometry_mapping_and_joint_projection(discs):
    d16,d12=discs
    assert d16.lengths==(.8,1.2) and d16.beta_deg==30
    assert d16.arms[0].coefficients.Bp<d16.arms[0].coefficients.Bb
    assert d16.arms[1].coefficients.Bp>d16.arms[1].coefficients.Bb
    raw=np.zeros(d16.ndof);raw[d16.joint_ids]=np.arange(7)*1e-5
    q16=d16.from_raw(raw);q12=r.projected_coordinates(d16,d12,q16)
    assert np.max(abs(d12.raw(q12)[d12.joint_ids]-raw[d16.joint_ids]))<1e-18
    for T,q in zip(d12.T7,d12.local(q12)):
        local=np.array([q[(f+1)*d12.p-1] for f in range(7)])
        assert np.max(abs(T@local-raw[d16.joint_ids]))<1e-18

def test_projected_subspace_energy_orthonormal_and_complete_output(discs):
    d16,d12=discs;basis,_,_=a.reference_subspace(d16);metric=r.fixed_metric(d12)
    projected,report=r.project_basis(d16,d12,basis[:,[0,2,3,5]],metric)
    assert projected.shape==(len(d12.state_out),4)
    assert report['orthonormality_error']<2e-10  # existing reference_subspace contract
    assert np.max(abs(projected.T@metric@projected-np.eye(4)))<2e-10

def test_planar_parity_algebra(discs):
    d,_=discs;plus=np.zeros((3,2*d.ndof));minus=plus.copy()
    plus[:,d.state_plane]=1e-6;minus[:,d.state_plane]=1e-6
    plus[:,d.state_out]=2e-8;minus[:,d.state_out]=-2e-8
    report=r.parity_metrics(plus,minus,d,r.fixed_metric(d),1e-5)
    assert report['out_normalized_energy_error']==0 and report['plane_absolute_energy_error']==0

def test_restricted_gain_uses_full_output(discs):
    d,_=discs;n=len(d.state_out);basis=np.eye(n)[:,:2];history=np.stack([basis,basis.copy()]);history[1,-1,0]=3
    report,gain=a.gain_result(d,history,np.array([0.,1.]),basis,np.eye(n))
    assert np.isclose(gain[-1],np.sqrt(10)) and report['full_outplane_response_dimension']==n
    assert report['restricted_initial_dimension']==2 and not report['infinite_dimensional_supremum_established']

def test_zero_base_energy_control_algebra():
    w=np.array([.4,1.7]);t=np.linspace(0,100,300);Z=np.zeros((len(t),4,4))
    for j in range(2):
        Z[:,j,j]=np.cos(w[j]*t);Z[:,j+2,j]=-np.sin(w[j]*t)
        Z[:,j,j+2]=np.sin(w[j]*t);Z[:,j+2,j+2]=np.cos(w[j]*t)
    gram=np.einsum('tij,tik->tjk',Z,Z)
    assert np.max(abs(gram-np.eye(4)))<1e-15

def test_frozen_full_accuracy_not_blanket_tangent_policy():
    c=p.read(r.CONFIG)
    assert c['full_rtol']==1e-12 and c['full_atol_factor']==.01 and c['full_max_step_factor']==.5
    assert c['tangent_rtol']==1e-11 and c['tangent_method']=='DOP853'
    assert inspect.signature(a.tangent_attempt).parameters['method'].default=='Radau'

def test_no_retry_after_stop_or_ceiling():
    c=p.read(r.CONFIG);s={'hard_stop':True,'calls':{'nonlinear_ODE':0,'tangent_ODE':0},'numerical_seconds':0}
    with pytest.raises(RuntimeError):r.limits(s,c)
    s['hard_stop']=False;s['calls']['nonlinear_ODE']=9
    with pytest.raises(RuntimeError):r.limits(s,c)

def test_cached_completed_run_has_zero_scientific_calls(monkeypatch):
    s={'completed':True,'hard_stop':False}
    monkeypatch.setattr(p,'sources',lambda *x:pytest.fail('source assembly on completed replay'))
    monkeypatch.setattr(p,'integrate_attempt',lambda *x:pytest.fail('new IVP on cache'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x:pytest.fail('new tangent on cache'))
    assert r.run(None,None,s,None)==s

def test_resume_dispatch_preserves_historical_modes(monkeypatch):
    seen=[];monkeypatch.setattr(r,'main',lambda args:seen.append(args))
    p.main(['--resume-sensitivity','--report-only','saved_bundle'])
    assert seen==[['--report-only','saved_bundle']]

def test_output_sampling_skips_dense_work_without_changing_samples(tmp_path,monkeypatch):
    """Mock accepted steps, not a scientific or synthetic ODE solve."""
    from types import SimpleNamespace
    import time
    instances=[]
    class AcceptedSteps:
        def __init__(self,fun,t0,initial,end,**kw):
            self.t=0.;self.status='running';self.initial=initial.copy();self.nfev=self.njev=self.nlu=0;self.dense_calls=0;instances.append(self)
        def step(self):
            self.t+=.25
            if self.t==1.:self.status='finished'
        def dense_output(self):
            self.dense_calls+=1
            return lambda t:np.tile(self.initial[:,None],(1,len(t)))
    d=SimpleNamespace(ndof=2,state_out=np.arange(4),arms=[SimpleNamespace(coefficients=SimpleNamespace(C=1.,jp=1.))],jacobian=lambda t,x:np.zeros((4,4)))
    monkeypatch.setattr(a,'DOP853',AcceptedSteps);monkeypatch.setattr(p,'snapshot_case',lambda case:None);monkeypatch.setattr(p,'finalize',lambda *x:None)
    item={'config':{'rtol':1e-11,'atol_relative':1e-11,'per_tangent_run_seconds':1}}
    outputs=[]
    for j,sampled in enumerate((False,True)):
        s={'cases':{},'attempts':[],'calls':{'tangent_ODE':0},'numerical_seconds':0.}
        output,_=a.tangent_attempt(tmp_path,item,s,str(j),d,None,np.array([0.,1.]),np.eye(4)[:,:2],np.eye(4),time.perf_counter()+2,
            zero=True,method='DOP853',sampled_dense_only=sampled)
        outputs.append(output)
    assert instances[0].dense_calls==4 and instances[1].dense_calls==1
    assert np.array_equal(outputs[0],outputs[1])

def test_physical_field_selection_keeps_section_axis():
    values=np.tile(np.arange(7),(2,201,1));values[1,170,3]=-40
    assert np.array_equal(r.physical_field_max(values),[2,40,4])

@pytest.mark.parametrize('mode',['--report-only','--plot-only'])
def test_completed_cached_cli_modes_call_no_scientific_functions(tmp_path,monkeypatch,mode):
    p.write(tmp_path/'figure_provenance.json',{'figures':[]});p.write(tmp_path/'export_provenance.json',{'scientific_calls':0})
    s={'completed':True,'hard_stop':False,'overall':'NUMERICAL_PARTIAL','stages':{},'calls':{'nonlinear_ODE':6,'tangent_ODE':5},'numerical_seconds':6520.}
    monkeypatch.setattr(r,'cache',lambda b:s)
    monkeypatch.setattr(p,'integrate_attempt',lambda *x:pytest.fail('cached nonlinear IVP'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x:pytest.fail('cached tangent IVP'))
    monkeypatch.setattr(p,'static_state',lambda *x:pytest.fail('cached equilibrium'))
    monkeypatch.setattr(p,'make_disc',lambda *x:pytest.fail('cached model reconstruction'))
    r.main([mode,str(tmp_path)])

def test_updated_resource_preflight_blocks_before_physics(tmp_path,monkeypatch):
    case=tmp_path/'B_p12_base';case.mkdir();p.write(case/'NOT_RUN_resource_preflight.json',{'new_integrations':0})
    item={'config':{'maximum_nonlinear_trajectories':8}};s={'cases':{},'calls':{'nonlinear_ODE':6},'attempts':[]}
    # Atomic case creation is before recording/starting an IVP. The reserved
    # administrative directory therefore cannot masquerade as a real attempt.
    with pytest.raises(FileExistsError):p.integrate_attempt(tmp_path,item,s,'B_p12_base',None,None,None,None,None)
    assert s['calls']['nonlinear_ODE']==6 and s['attempts']==[]
