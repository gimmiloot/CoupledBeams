"""No physical IVPs: saved states and synthetic continuation contracts only."""
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.lib import nlsp_construction_c_longer_sensitivity as l
from scripts.lib import nlsp_coupled_sensitivity_analysis as a
from scripts.lib import nlsp_coupled_sensitivity_pilot as p

def test_frozen_scope():
    config=p.read(l.CONFIG);assert l.validate_scope(config)==config
    for key,value in [('numerical_budget_seconds',2701),('degree',16),('automatic_retry',True),('p16_new_maximum_trigger',.004)]:
        with pytest.raises(ValueError):l.validate_scope({**config,key:value})

def test_physical_time_grid():
    old=np.linspace(0,40,2001);times=l.continuation_times(old,20,80.83)
    assert times[0]==40 and times[-1]==80.83 and np.all(np.diff(times)>0)
    assert np.max(np.diff(times))<=.020000000001

def test_join_preserves_all_prefix_values():
    t,x=l.merge_prefix(np.array([0.,1.]),np.array([[1.,2.],[3.,4.]]),np.array([1.,2.]),np.array([[3.,4.],[5.,6.]]))
    assert np.array_equal(t,[0,1,2]) and np.array_equal(x,[[1,2],[3,4],[5,6]])
    with pytest.raises(ValueError):l.merge_prefix(np.array([0.,1.]),np.ones((2,2)),np.array([1.,2.]),np.zeros((2,2)))

def test_gain_not_renormalized_at_restart():
    d=SimpleNamespace(out=np.array([0]),ndof=1,reference_energy_parts=lambda q,v:np.array([[q[0]**2,v[0]**2]]))
    B=np.array([[1.],[0.]]);Z=np.array([[[2.],[0.]],[[3.],[0.]]])
    result,gain=a.gain_result(d,Z,np.array([4.,5.]),B,np.eye(2))
    assert np.array_equal(gain,[2.,3.]) and np.array_equal(result['initial_state_out'],B[:,0])

def test_predeclared_p16_gate():
    config=p.read(l.CONFIG);times=np.linspace(0,4,4001)
    small=np.ones(len(times))*1.053
    assert not l.growth_gate(times,small,1.053,1.,config)['p16_scientific_trigger']
    big=small.copy();big[100]=1.06
    assert l.growth_gate(times,big,1.053,1.,config)['p16_scientific_trigger']
    rising=small.copy();rising[times>3.8]+=np.linspace(0,.002,(times>3.8).sum())
    assert l.growth_gate(times,rising,1.053,1.,config)['clearly_rising_near_end']

@pytest.fixture(scope='module')
def saved_model():
    pc=p.read(l.PARENT/'provenance.json')['physical_source_config'];return p.sources(pc)[0]

@pytest.mark.parametrize('degree',[12,16])
def test_saved_restart_identity_and_same_operator(saved_model,degree):
    d,t,state,Z,B,E,energy,info=l.load_level(degree,saved_model)
    assert info['source_dense_endpoint_max_difference']<2e-14
    assert np.array_equal(Z[0],B) and np.array_equal(state[0,d.state_out],np.zeros(len(d.state_out)))
    for idx in (0,745,2000):
        old=d.jacobian(t[idx],state[idx])[np.ix_(d.state_out,d.state_out)]
        assert np.allclose(old,l.out_operator(d,t[idx],state[idx]),rtol=2e-13,atol=2e-13)
    zero=l.out_operator(d,0,np.zeros(2*d.ndof));assert np.linalg.norm(zero.T@E+E@zero)/np.linalg.norm(E@zero)<2e-12

def test_historical_manifest_unchanged():
    assert p.sha(l.PARENT/'manifest.json')==l.PARENT_SHA

def test_shared_time_validation_has_no_ODE(tmp_path):
    with pytest.raises(ValueError):p.integrate_attempt(tmp_path,{}, {'cases':{}},'bad',None,None,np.array([1.,2.]),{},0,t_start=0)
    with pytest.raises(ValueError):a.tangent_attempt(tmp_path,{}, {'cases':{}},'bad',None,None,np.array([1.,2.]),None,None,0,t_start=0)

def test_tangent_restart_uses_endpoint_not_basis(tmp_path,monkeypatch):
    calls=[]
    class FakeDOP:
        def __init__(self,fn,t,y,end,**kw):
            calls.append((t,y.copy()));self.t=t;self.y=y.copy();self.end=end;self.status='running';self.nfev=self.njev=self.nlu=0
        def step(self):self.t=self.end;self.status='finished'
        def dense_output(self):return lambda t:np.repeat(self.y[:,None],len(np.atleast_1d(t)),axis=1)
    d=SimpleNamespace(ndof=1,state_out=np.array([0,1]),arms=[SimpleNamespace(coefficients=SimpleNamespace(C=1,jp=1))],jacobian=lambda *x:np.array([[0,1],[-1,0]]))
    monkeypatch.setattr(a,'DOP853',FakeDOP);monkeypatch.setattr(a,'SavedPlane',lambda *x:lambda t:np.zeros(2))
    monkeypatch.setattr(p,'snapshot_case',lambda *x:None);monkeypatch.setattr(p,'finalize',lambda *x:None)
    s={'calls':{'tangent_ODE':0},'cases':{},'attempts':[],'numerical_seconds':0.}
    item={'config':{'per_tangent_run_seconds':20,'rtol':1e-11,'atol_relative':1e-11}}
    basis=np.array([[1.],[0.]]);restart=np.array([[.3],[.9]])
    history,record=a.tangent_attempt(tmp_path,item,s,'restart',d,tmp_path,np.array([4.,5.]),basis,np.eye(2),1e30,
        method='DOP853',t_start=4.,initial_matrix=restart)
    assert calls[0][0]==4 and np.array_equal(calls[0][1],restart[:,0])
    assert np.array_equal(history[0],restart)
    with np.load(tmp_path/'restart/tangent.npz') as z:assert np.array_equal(z['basis'],basis)

def test_restart_preserves_solver_time_and_physical_state(tmp_path,monkeypatch):
    seen=[]
    class FakeRadau:
        def __init__(self,fun,t,y,end,**kw):
            seen.append((t,y.copy(),kw));self.t=t;self.y=y.copy();self.end=end;self.status='running';self.nfev=self.njev=self.nlu=0
        def step(self):self.t=self.end;self.status='finished'
        def dense_output(self):
            return SimpleNamespace(t_old=4.,t=5.,y_old=self.y.copy(),Q=np.zeros((4,3)),__call__=None)
    class Dense:
        t_old=4.;t=5.;Q=np.zeros((4,3))
        def __init__(self,y):self.y_old=y
        def __call__(self,t):return np.repeat(self.y_old[:,None],len(np.atleast_1d(t)),axis=1)
    FakeRadau.dense_output=lambda self:Dense(self.y)
    monkeypatch.setattr(p,'Radau',FakeRadau);monkeypatch.setattr(p,'snapshot_case',lambda *x:None);monkeypatch.setattr(p,'finalize',lambda *x:None)
    d=SimpleNamespace(ndof=2,state_plane=np.arange(4),reset_counters=lambda:None,safety=lambda *x:[],plane_rhs=lambda *x:np.zeros(4),plane_jacobian=lambda *x:np.eye(4),energy=lambda *x:1.,counters=lambda:{},embed_plane=lambda x:x)
    item={'config':{'maximum_nonlinear_trajectories':2,'per_nonlinear_run_seconds':20}}
    s={'calls':{'nonlinear_ODE':0},'cases':{},'attempts':[],'numerical_seconds':0.}
    initial=np.array([.1,.2,.3,.4]);timing={'rtol':1e-10,'atol':np.ones(4)*1e-10,'max_step':.01}
    states,record=p.integrate_attempt(tmp_path,item,s,'restart',d,initial,np.array([4.,5.]),timing,1e30,plane=True,t_start=4.)
    assert seen[0][0]==4 and np.array_equal(seen[0][1],initial)
    assert np.array_equal(states[0],initial) and record['physical_start']==4
    with np.load(tmp_path/'restart/accepted_times.npz') as z:assert np.array_equal(z['times'],[4,5]) and np.array_equal(z['actual_final_state'],initial)

@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_cached_modes_zero_scientific_calls(monkeypatch,mode,tmp_path):
    summary={'completed':True,'hard_stop':False,'overall':'DONE','stages':{},'calls':{},'numerical_seconds':0}
    monkeypatch.setattr(l,'cache',lambda *x:summary);monkeypatch.setattr(l,'prepare',lambda *x:(tmp_path,{},summary))
    monkeypatch.setattr(p,'sources',lambda *x:pytest.fail('scientific reconstruction on replay'))
    monkeypatch.setattr(p,'make_disc',lambda *x:pytest.fail('model construction on replay'))
    monkeypatch.setattr(p,'integrate_attempt',lambda *x,**kw:pytest.fail('new nonlinear integration'))
    monkeypatch.setattr(a,'tangent_attempt',lambda *x,**kw:pytest.fail('new tangent integration'))
    l.main([mode] if mode=='--compute' else [mode,str(tmp_path)])

def test_no_retry_after_hard_stop(monkeypatch,tmp_path):
    summary={'completed':False,'hard_stop':True}
    monkeypatch.setattr(p,'sources',lambda *x:pytest.fail('new work after solver failure'))
    assert l.run(tmp_path,{},summary) is summary

def test_dense_absolute_start_contract(tmp_path):
    from scripts.lib.nlsp_fem3b_continuation import save_dense_records,evaluate_dense_records
    records=[{'t_old':4.,'t':5.,'y_old':np.array([2.]),'Q':np.array([[1.,0.,0.]])}]
    with pytest.raises(ValueError):save_dense_records(tmp_path/'wrong.npz',records)
    save_dense_records(tmp_path/'right.npz',records,t_start=4.)
    assert np.array_equal(evaluate_dense_records(tmp_path/'right.npz',np.array([4.,4.5,5.])),[[2.],[2.5],[3.]])

def test_single_explicit_output_repeat_keeps_budget_and_old_failure():
    config=p.read(l.REPEAT_CONFIG);assert l.validate_scope(config)==config
    assert config['numerical_budget_seconds']==2700 and config['previous_numerical_seconds_charged']==450
    assert config['additional_p12_base_repeat_authorized']==1
    source=p.ROOT/config['archive_failure_bundle']
    assert p.sha(source/'manifest.json')==config['archive_failure_manifest_sha256']
    assert p.read(source/'summary.json')['calls']['nonlinear_ODE']==1
    assert p.read(source/'summary.json')['hard_stop']
