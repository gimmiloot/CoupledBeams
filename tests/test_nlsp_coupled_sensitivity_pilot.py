"""Scoped source/cache/metric tests, with no hidden nonlinear trajectories."""
from pathlib import Path
import copy,json
import time
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.lib import nlsp_coupled_sensitivity_pilot as pilot
from scripts.lib import nlsp_coupled_sensitivity_analysis as analysis

def test_config_has_only_two_frozen_geometries_and_no_FEM():
    c=pilot.read(pilot.CONFIG);assert pilot.validate_config(c) is c
    assert c['degrees']==[16,24] and c['horizon_Tref']==2
    assert c['new_FEM_jobs']==0 and c['maximum_nonlinear_trajectories']==23
    assert c['deferred_theory']=='DEFERRED_TWO_PARAMETER_ACTION_EXPANSION'

@pytest.mark.parametrize('key,value',[('beta_deg',0),('horizon_Tref',5),('rtol',1e-6),('Floquet',True),('automatic_retry',True),('new_FEM_jobs',1),('initial_subspace_dimension',200)])
def test_unannounced_config_changes_rejected(key,value):
    c=pilot.read(pilot.CONFIG);c[key]=value
    with pytest.raises(ValueError):pilot.validate_config(c)

def test_common_times_do_not_shift_phase_or_shorten_horizon():
    c=pilot.read(pilot.CONFIG);t=pilot.make_times(7.3,2,c)
    assert t[0]==0 and t[-1]==14.6 and len(t)==2001
    assert np.all(np.diff(t)>0)

def test_reference_energy_norm_and_gain_preserve_physical_weighting():
    metric=np.diag([4.,.25]);basis=np.diag([.5,2.]);history=np.array([basis,2*basis])
    class Fake:
        out=np.array([0]);ndof=1
        def reference_energy_parts(self,q,v):return np.array([[q[0]**2*4,v[0]**2*.25]])
    result,curve=analysis.gain_result(Fake(),history,np.array([0.,1.]),basis,metric)
    assert np.allclose(curve,[1.,2.]) and result['restricted_initial_dimension']==2
    assert result['lower_bound_for_full_discrete_initial_space'] and not result['infinite_dimensional_supremum_established']

def test_dense_base_uses_actual_cubic_and_preserves_prefix(tmp_path):
    p=tmp_path/'dense.npz';np.savez(p,t_old=[0.],t=[2.],y_old=[[1.,2.]],Q=[[[3.,4.,5.],[6.,7.,8.]]])
    class Fake:
        def embed_plane(self,y):return y
    dense=analysis.SavedPlane(p,Fake());assert np.allclose(dense(1.),[1+1.5+1+.625,2+3+1.75+1])
    with pytest.raises(ValueError):dense(2.01)

def test_orders_are_of_normalized_odd_response_and_even_planar_change():
    rows=[{'delta':d,'normalized_out_error_relative':d*d,'planar_change_energy_norm_max':3*d*d} for d in [1.,.5,.25]]
    for row in analysis.measured_orders(rows):assert np.isclose(row['normalized_out_order'],2) and np.isclose(row['planar_change_order'],2)

def test_two_amplitude_postprocessing_is_conditional_and_does_not_override_fail():
    metric=np.eye(2);tangent=np.array([[1.,0.],[0.,1.]])
    floor=np.array([[0.,.1],[.1,0.]]);cubic=np.array([[.02,0.],[0.,.02]])
    first=tangent+floor+cubic;second=tangent+floor+cubic/4
    result=analysis.two_amplitude_diagnostic(first,second,tangent,metric)
    assert np.isclose(result['conditional_limit_discrepancy'],.1)
    assert result['does_not_override_frozen_gate'] and 'diagnostic only' in result['assumption']

def test_sparse_batch_layout_is_columnwise_and_retains_full_response():
    from scipy.sparse import kron,eye,csr_matrix
    A=np.array([[1.,2.],[3.,4.]]);Z=np.arange(6.).reshape(2,3)
    assert np.allclose(kron(eye(3),csr_matrix(A))@Z.flatten(order='F'),(A@Z).flatten(order='F'))
    dense=np.stack([Z.flatten(order='F'),(2*Z).flatten(order='F')])
    assert np.array_equal(dense.reshape(2,2,3,order='F'),np.stack([Z,2*Z]))

def test_cached_report_and_plot_do_not_prepare_or_compute(tmp_path,monkeypatch):
    summary={'overall':'NUMERICAL_PARTIAL','stages':{},'calls':{},'stop_reason':None}
    monkeypatch.setattr(pilot,'cache',lambda b:summary)
    monkeypatch.setattr(pilot,'prepare',lambda *a:pytest.fail('unexpected prepare/scientific computation'))
    monkeypatch.setattr(pilot,'run_program',lambda *a:pytest.fail('unexpected new trajectories'))
    (tmp_path/'figure_provenance.json').write_text('{}')
    pilot.main(['--report-only',str(tmp_path)]);pilot.main(['--plot-only',str(tmp_path)])

def test_attempt_guard_independent_of_new_code_hash():
    source=Path(pilot.__file__).read_text(encoding='utf8')
    assert 'code/HEAD changes do not authorize a retry' in source
    assert "if label in s['cases']:raise ValueError" in source
    assert 'no production if the estimate exceeds' in source

def planning_fixture(tmp_path):
    c=pilot.read(pilot.CONFIG);item={'config':c,'source_manifests':{}}
    s={'overall':'NUMERICAL_FEASIBILITY_PARTIAL','completed':False,'hard_stop':True,'stop_reason':pilot.PLANNING_STOP,
       'calls':{'nonlinear_ODE':1,'tangent_ODE':0,'CCX':0,'Gmsh':0},'cases':{'cost_probe':{'status':'PASS','failure':None}},
       'attempts':[{'label':'cost_probe','status':'PASS'}],'stages':dict.fromkeys('CDEF','NOT_RUN')}
    (tmp_path/'cost_probe').mkdir();pilot.write(tmp_path/'cost_probe/integration.json',{'status':'PASS'})
    pilot.write(tmp_path/'cost_preflight.json',{'total_planning_seconds':10000})
    for name,data in [('summary',s),('config',c),('manifest',{'first_decision':True})]:pilot.write(tmp_path/(name+'.json'),data)
    r={'original_cost_estimate_sha256':pilot.sha(tmp_path/'cost_preflight.json'),'saved_probe_sha256':pilot.sha(tmp_path/'cost_probe/integration.json'),
       'new_ODE_calls':0,'budget_seconds':7200,'physical_horizons_preserved':True,'scientific_gates_or_tolerances_changed':False,
       'decision':'PRODUCTION_WITHIN_ORIGINAL_BUDGET','total_planning_seconds':5590}
    pilot.write(tmp_path/'cost_preflight_reassessment.json',r);return item,s,r

def test_planning_correction_preserves_first_decision_without_new_scientific_calls(tmp_path):
    item,s,r=planning_fixture(tmp_path);first=(tmp_path/'summary.json').read_bytes();forecast=(tmp_path/'cost_preflight.json').read_bytes()
    assert pilot.reassess_planning_stop(tmp_path,item,s)
    assert not s['hard_stop'] and s['calls']['nonlinear_ODE']==1
    assert (tmp_path/'first_cost_decision/summary.json').read_bytes()==first
    assert (tmp_path/'cost_preflight.json').read_bytes()==forecast
    assert not pilot.reassess_planning_stop(tmp_path,item,s)

@pytest.mark.parametrize('change',['solver_failure','failed_probe','production_exists','corrupt_hash','budget_changed','missing_reassessment'])
def test_planning_correction_cannot_bypass_actual_failure_or_change_budget(tmp_path,change):
    item,s,r=planning_fixture(tmp_path)
    if change=='solver_failure':s['stop_reason']='Actual Radau failure'
    elif change=='failed_probe':s['attempts'][0]['status']='FAIL'
    elif change=='production_exists':s['calls']['nonlinear_ODE']=2
    elif change=='corrupt_hash':r['saved_probe_sha256']='corrupt'
    elif change=='budget_changed':r['budget_seconds']=10000
    elif change=='missing_reassessment':(tmp_path/'cost_preflight_reassessment.json').unlink()
    if change!='missing_reassessment':pilot.write(tmp_path/'cost_preflight_reassessment.json',r)
    assert not pilot.reassess_planning_stop(tmp_path,item,s)
    assert s['hard_stop'] and not (tmp_path/'first_cost_decision').exists()

@pytest.mark.parametrize('fail',[False,True])
def test_attempt_records_actual_prefix_without_retry_or_tail_padding(tmp_path,monkeypatch,fail):
    class Disc:
        ndof=1;state_plane=np.array([0,1])
        def reset_counters(self):pass
        def counters(self):return {}
        def safety(self,q):return []
        def embed_plane(self,x):return x
        def plane_rhs(self,t,x):pytest.fail('Synthetic fixture must not evaluate physical RHS')
        def plane_jacobian(self,t,x):pytest.fail('Synthetic fixture must not evaluate physical Jacobian')
        def energy(self,q,v):return float(q@q+v@v)
    class Fake:
        def __init__(self,fun,t,y,bound,**kwargs):
            self.t=t;self.y=y;self.status='running';self.nfev=self.njev=self.nlu=0
        def step(self):
            if fail and self.t==.1:self.status='failed';return 'Synthetic native failure'
            self.old=self.t;self.t+=.1
            if self.t>=.2:self.status='finished'
        def dense_output(self):
            class Dense:
                def __call__(self,times):return np.tile(self.y_old[:,None],(1,len(times)))
            dense=Dense();dense.t_old=self.old;dense.t=self.t;dense.y_old=self.y.copy();dense.Q=np.zeros((2,3));return dense
    monkeypatch.setattr(pilot,'Radau',Fake);monkeypatch.setattr(pilot,'snapshot_case',lambda case:None)
    c=pilot.read(pilot.CONFIG);item={'config':c,'source_manifests':{}}
    s={'cases':{},'attempts':[],'calls':{'nonlinear_ODE':0},'numerical_seconds':0.,'hard_stop':False}
    states,record=pilot.integrate_attempt(tmp_path,item,s,'synthetic',Disc(),np.array([1.,0.]),np.array([0.,.1,.2]),
        {'rtol':1e-10,'atol':[1e-10,1e-10],'max_step':.1},time.perf_counter()+10,plane=True)
    count=2 if fail else 3
    assert len(states)==count and record['sampled_end']==(.1 if fail else .2)
    assert s['calls']['nonlinear_ODE']==1 and len(s['attempts'])==1
    assert s['hard_stop']==fail
    with np.load(tmp_path/'synthetic/trajectory.npz') as z:assert len(z['times'])==count
    with pytest.raises(ValueError,match='no automatic nonlinear retry'):
        pilot.integrate_attempt(tmp_path,item,s,'synthetic',Disc(),np.array([1.,0.]),np.array([0.,.1,.2]),{},time.perf_counter()+10)

def test_saved_actual_torsional_failure_is_not_promoted_or_retried(monkeypatch):
    b=pilot.OUTPUT/'d024ee84989ce6b5'
    if not (b/'stopped_prefix_diagnostics.json').exists():pytest.skip('Local generated stopped-prefix evidence unavailable')
    s=pilot.cache(b);r=pilot.read(b/'stopped_prefix_diagnostics.json')
    assert s['overall']=='BLOCKED_BY_TANGENT_CHECK' and s['hard_stop']
    assert s['calls']=={'nonlinear_ODE':6,'tangent_ODE':2,'CCX':0,'Gmsh':0}
    assert s['stages']['C']=='PARTIAL' and s['stages']['D']=='FAIL'
    assert r['shapes']['torsional']['absolute_relative_gate_pass']
    assert not r['shapes']['torsional']['declared_monotone_gate_pass']
    assert r['shapes']['torsional']['observed_out_error_status']=='FAIL'
    monkeypatch.setattr(pilot,'make_disc',lambda *a,**k:pytest.fail('Read-only saved diagnostics must not assemble or solve'))
    assert analysis.diagnose_stopped_prefix(b)==r
    for case in s['cases'].values():
        assert case['status']=='PASS' and case['actual_accepted_end']==case['target_end']
    monkeypatch.setattr(pilot,'Radau',lambda *a,**k:pytest.fail('No retry'))
    assert pilot.run_program(b,pilot.read(b/'provenance.json'),s,None) is s
