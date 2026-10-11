"""Diagnostic algebra, scope and cached dispatch; no physical IVP in tests."""
import numpy as np
import pytest
from scripts.lib import nlsp_coupled_torsion_diagnostic as g
from scripts.lib import nlsp_coupled_sensitivity_pilot as p

def test_scope_preserves_parent_failure_and_has_only_one_case():
    c=p.read(g.CONFIG)
    assert c['geometry']=='A' and c['degree']==16 and c['horizon_Tref']==.1
    assert c['maximum_new_ODE_integrations']==6 and c['numerical_budget_seconds']==1200
    assert c['new_static_equilibria']==c['new_eigenanalyses']==c['new_FEM_jobs']==0
    assert not c['automatic_retry'] and c['historical_status_unchanged']=='BLOCKED_BY_TANGENT_CHECK'

@pytest.mark.parametrize('key,value',[('geometry','B'),('degree',24),('maximum_new_ODE_integrations',7),('new_static_equilibria',1),('authorization','new-code-hash')])
def test_new_hash_or_config_does_not_expand_authorized_scope(key,value):
    c=p.read(g.CONFIG);c[key]=value
    with pytest.raises(ValueError,match='scope cannot be expanded'):g.validate_scope(c)

def test_coordinate_velocity_squared_norm_partition():
    metric=np.diag([4.,9.,1.,16.]);error=np.array([[1.,2.,3.,4.]])
    r=g.norm_curves(error,metric)
    assert np.isclose(r['total'][0]**2,305)
    assert np.isclose(r['coordinate'][0]**2+r['velocity'][0]**2,r['total'][0]**2)
    assert np.isclose(r['velocity_fraction'][0],265/305)

def test_joint_internal_partition_keeps_signed_cross_terms():
    x=np.array([1.,-.2]);metric=np.array([[1.,.9],[.9,1.]])
    r=g.bilinear_partition(x,metric,{'internal':[0],'joint':[1]})
    assert np.isclose(sum(r['signed_fractions'].values()),1.)
    assert r['signed_fractions']['joint']<0
    assert np.isclose(r['squared_norm'],sum(r['diagonal_projected_squared_norms'].values())+r['cross_squared_norm'])

def test_comparison_has_all_four_pairings_on_one_fixed_denominator():
    metric=np.eye(2);deltas=np.array([2.,1.]);old=np.array([[[2.,0.]],[[1.,0.]]]);new=old.copy();new[:,:,1]=deltas[:,None]*.1
    table,_=g.mismatch_matrix(old,new,np.array([[1.,0.]]),np.array([[1.,.04]]),deltas,metric,2.)
    for row in table:
        assert row['old_full/old_tangent']==0 and np.isclose(row['tight_full/tight_tangent'],.03)
        assert np.isclose(row['D_full'],.05) and np.isclose(row['D_tangent'],.02)

def test_common_grid_vector_extraction_preserves_time_axis():
    full=np.arange(2*5*8).reshape(2,5,8);ids=np.array([2,4,6])
    assert full[0][:,ids].shape==(5,3) and full[:,:,ids].shape==(2,5,3)

def test_conditional_variant_B_selects_only_largest_observed_sensitivity():
    assert g.choose_conditional_level([{'D_full':5.,'D_tangent':1.},{'D_full':8.,'D_tangent':1.}])==('full',1)
    assert g.choose_conditional_level([{'D_full':5.,'D_tangent':9.},{'D_full':8.,'D_tangent':9.}])==('tangent',None)

def test_refinement_frozen_before_new_integrations(tmp_path):
    parent=p.ROOT/p.read(g.CONFIG)['parent_bundle'];item={'config':p.read(g.CONFIG)}
    p.write(tmp_path/'summary.json',{'attempts':[],'calls':{'nonlinear_ODE':0,'tangent_ODE':0}})
    plan=g.freeze_refinement(tmp_path,item,parent)
    assert plan['rtol']==1e-12 and np.isclose(plan['max_step'],.0018023482469251596)
    assert plan['estimated_primary_seconds']<1200 and plan['frozen_before_results']
    assert np.array_equal(np.array(plan['full_atol']),np.array(p.read(parent/'A_check_torsional_0/time_settings.json')['atol'])*.01)

def test_new_freeze_cannot_use_results(tmp_path):
    p.write(tmp_path/'summary.json',{'attempts':[{'label':'already_run'}],'calls':{'nonlinear_ODE':1}})
    with pytest.raises(ValueError,match='after integration results'):g.freeze_refinement(tmp_path,{'config':p.read(g.CONFIG)},p.ROOT/p.read(g.CONFIG)['parent_bundle'])

def test_diagnostic_dispatch_preserves_old_default(monkeypatch):
    called=[];monkeypatch.setattr(g,'main',lambda args:called.append(args))
    p.main(['--torsion-diagnostic','--compute','--through-stage','read-only'])
    assert called==[['--compute','--through-stage','read-only']]

def test_cached_report_does_not_prepare_or_integrate(tmp_path,monkeypatch):
    s={'overall':'UNRESOLVED','stages':{},'calls':{}}
    monkeypatch.setattr(g,'cache',lambda b:s)
    monkeypatch.setattr(g,'prepare',lambda *args:pytest.fail('Cached report prepared computation'))
    monkeypatch.setattr(g,'run_primary',lambda *args:pytest.fail('Cached report integrated'))
    g.main(['--report-only',str(tmp_path)])

def test_actual_saved_refinement_preserves_initial_states_and_parent():
    b=g.OUTPUT/'a10f1110ffbce7ca'
    if not (b/'final_metrics.json').exists():pytest.skip('Local generated diagnostic unavailable')
    s=g.cache(b);metrics=p.read(b/'final_metrics.json');parent=g.verify_parent(p.read(g.CONFIG))
    assert s['calls']=={'nonlinear_ODE':3,'tangent_ODE':1,'CCX':0,'Gmsh':0}
    assert p.read(parent/'summary.json')['overall']=='BLOCKED_BY_TANGENT_CHECK'
    assert np.isclose(metrics['primary_observed_order'],2,atol=.001)
    table=p.read(b/'primary_temporal_comparison.json')['rows']
    assert all(row['D_full']>100*row['D_tangent'] for row in table)
    for j in range(2):
        with np.load(parent/f'A_check_torsional_{j}/initial.npz') as old,np.load(b/f'tight_full_{j}/initial.npz') as new:
            assert np.array_equal(old['state'],new['state']) and np.array_equal(old['times'],new['times'])
    with np.load(b/'tight_full_1/initial.npz') as old,np.load(b/'extra_full_1/initial.npz') as new:assert np.array_equal(old['state'],new['state'])

def test_cached_plot_does_not_reassemble_or_integrate(tmp_path,monkeypatch):
    s={'overall':'NUMERICAL_FLOOR_IDENTIFIED','stages':{},'calls':{}}
    (tmp_path/'figure_provenance.json').write_text('{}')
    monkeypatch.setattr(g,'cache',lambda b:s)
    monkeypatch.setattr(g,'postprocess_saved',lambda *args:pytest.fail('Cached plot recalculated diagnostics'))
    monkeypatch.setattr(g,'plot_saved',lambda *args:pytest.fail('Cached plot rendered again'))
    g.main(['--plot-only',str(tmp_path)])
