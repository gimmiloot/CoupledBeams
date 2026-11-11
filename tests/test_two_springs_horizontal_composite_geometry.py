"""HS-D2 geometry identity, archive-only comparisons and root evidence."""
import math
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_composite_geometry as g

FOLDER=g.output_folder()


def test_exact_scope_and_unchanged_numerical_contract():
    cfg=g.read(g.CONFIG);old=g.read(g.d1.CONFIG)
    assert [c['id'] for c in cfg['cases']]==g.expected_ids()
    assert len(cfg['cases'])==34
    assert sum(c['role']=='control_I' for c in cfg['cases'])==9
    assert cfg['policy']==old['policy'] and cfg['orders']==[20,24]
    assert cfg['new_full_solves_allowed']==0 and cfg['fit_parameters']==0
    assert sum(map(len,cfg['angles_deg'].values()))==2880
    assert g.d1.output_folder().name=='60fc0dd404dcf167'


def test_support_scaling_and_cache_identity_are_model_specific():
    cfg=g.read(g.CONFIG);cases={c['id']:c for c in cfg['cases']}
    for name,value in [('B1_A2_e0.0025',16),('B0_A2',4),('B1_A2_e0.006',25/9),('B4_A2',25/9)]:
        assert cases[name]['kappa_s']==pytest.approx(value,rel=1e-14)
    a,b=cases['B0_A1'],cases['B3_A1_a2']
    assert a['aspect']!=b['aspect'] and g.mechanical_identity(a)==g.mechanical_identity(b)
    assert g.mechanical_identity(a)!=g.mechanical_identity(cases['B1_A1_e0.0025'])
    assert g.point_identity('C',.1,.01,30,4)!=g.point_identity('C',.1,.01,30,g.RIGID)
    assert g.point_identity('C',.1,.01,30,4)!=g.point_identity('C',.1,-.01,30,4)
    assert g.point_identity('L1',.1,.01,4)!=g.point_identity('L1',.1,.01,16)


def test_all_actual_full_nodes_and_frequency_units():
    rows=g.rows_csv(FOLDER/'comparison.csv');cfg=g.read(g.CONFIG)
    assert len(rows)==2880
    for case in cfg['cases']:
        rr=[r for r in rows if r['id']==case['id']]
        assert [float(r['beta_deg']) for r in rr]==cfg['angles_deg'][case['id']]
        for r in rr:
            assert float(r['omega_primary'])==pytest.approx(case['epsilon']*math.sqrt(float(r['zeta_primary'])),rel=1e-14)
            assert float(r['Lambda_primary'])**4==pytest.approx(float(r['zeta_primary']),rel=1e-14)
            assert float(r['omega2_full'])>float(r['omega_full'])>0
            if case['role']=='control_I':assert r['zeta_comp']=='' and r['omega_comp']=='' and r['primary']=='C_I'
            else:assert float(r['zeta_comp'])>0


def test_fifteen_distinct_common_limits_and_normalized_forms():
    coefficients=g.read(FOLDER/'coefficients.json');assert len(coefficients)==15
    assert len({(r['k_joint'],r['mu']) for r in coefficients})==15
    for r in coefficients:
        for order in r['orders']:
            assert order['common_limit_frequency_error']<=1e-7
            assert order['global_composite_lower_bound']>0
            assert order['C2']==pytest.approx(-order['zeta_C_zero_plus']*order['A']*order['integral_V'],rel=2e-14)
        assert abs(r['C2_Richardson_relative_error'])<=5e-5
        assert (FOLDER/'points'/(r['inner_limit_key']+'.npz')).is_file()
        assert (FOLDER/'points'/(r['outer_limit_key']+'.npz')).is_file()


def test_direct_threshold_evidence_and_reference_semantics():
    predictions=g.read(FOLDER/'predictions.json');assert len(predictions)==68
    assert len(g.rows_csv(FOLDER/'threshold_comparison.csv'))==68
    for r in predictions:
        assert r['high_at_zero'] is False
        assert r['terminal_status']=='TERMINAL_HIGH_RANGE'
        for c in r['crossings']:
            assert c['signs'][0]*c['signs'][1]<0
            assert c['bracket'][1]-c['bracket'][0]==pytest.approx(2e-6,rel=1e-8)
            assert abs(c['residual'])<=1e-9
            assert abs(c['differences'][-1]['relative_error'])<=5e-5
    m=g.read(FOLDER/'manifest.json');assert m['new_full_solves']==0 and m['full_thresholds']==68


def test_D1_reproduction_and_aspect_caveat():
    r=g.read(FOLDER/'HS_D1_regression.json')
    assert r['status']=='PASS' and len(r['curves'])==524 and len(r['thresholds'])==12
    assert max(v for row in r['curves'] for v in row['errors'].values())<=1e-7
    assert max(row['absolute_error_deg'] for row in r['thresholds'])<=2e-6
    aspect=g.rows_csv(FOLDER/'aspect_comparison.csv');assert len(aspect)==4
    assert 'UNRESOLVED' in next(r for r in aspect if r['regime']=='A2')['qualification']
    assert all(float(r['primary_relative_difference'])==0 for r in aspect)


def test_error_statistics_keep_sign_and_absolute_errors_separate():
    s=g.error_statistics([-3,1,2])
    assert s['max_abs']==3 and s['mean_abs']==2 and s['median_abs']==2
    assert s['p95_abs']==pytest.approx(2.9) and s['signed_bias']==0


def test_tiny_decrease_is_diagnosed_without_replacing_spectra():
    d=g.read(FOLDER/'local_neighbour_diagnostics.json')
    assert d['status']=='PASS' and d['original_spectra_preserved']
    assert d['new_full_solves']==d['new_ritz_solves']==0
    for event in d['events']:
        assert event['event']['decrease']
        assert event['boundary_relative_change']>0
        assert all(v['analytic_slope_log_per_degree']>0 for v in event['values'])
        rr=g.rows_csv(FOLDER/'comparison.csv')
        for v in event['values']:
            row=next(r for r in rr if r['id']==event['event']['id'] and float(r['beta_deg'])==v['beta_deg'])
            assert float(row['omega_primary'])==v['omega_original']


def test_plot_only_forbids_all_spectral_solutions(monkeypatch,tmp_path):
    for name in ['contract.json','comparison.csv','frequency_statistics.csv','threshold_comparison.csv','aspect_common_points.csv']:
        (tmp_path/name).write_bytes((FOLDER/name).read_bytes())
    def forbidden(*args,**kwargs):raise AssertionError('Solver in plot_only')
    monkeypatch.setattr(g.d1.inner.Leading,'solve',forbidden)
    monkeypatch.setattr(g.d1.outer.FiniteAngle,'solve',forbidden)
    monkeypatch.setattr(g.GeometryStudy,'__init__',forbidden)
    g.render(tmp_path)
    assert g.read(tmp_path/'figures_manifest.json')['new_spectral_solves']==0
    assert len(list((tmp_path/'figures').glob('*.pdf')))==7
