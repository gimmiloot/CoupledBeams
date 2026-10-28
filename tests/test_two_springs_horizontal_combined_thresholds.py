"""HS-LA6 regression checks: saved forecasts, no full-model computations."""
import json
import math
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_combined_thresholds as study

CFG=json.loads(study.CONFIG.read_text())
FOLDER=study.BASE/study.identity(CFG)[0]
PRED=json.loads((FOLDER/'apriori_predictions.json').read_text())


def test_targets_below_saturation_and_isolated_roots():
    for r in PRED['thresholds']:
        assert r['target_zeta']==pytest.approx(r['q']**2*r['zeta_outer'],rel=1e-14)
        assert r['delta']>0 and r['target_zeta']<r['zeta_infty']
        assert r['tau_lower']<r['tau']<r['tau_upper']
        assert r['sign_lower']<-CFG['policy']['bracket_sign_min']
        assert r['sign_upper']>CFG['policy']['bracket_sign_min']
        assert abs(r['equation_residual'])<CFG['policy']['threshold_equation_tol']


def test_large_two_term_root_and_scaling():
    for r in PRED['thresholds']:
        assert r['discriminant']>0 and 0<r['tau_other']<r['tau_two']
        assert r['a']/r['tau_two']+r['b']/r['tau_two']**2==pytest.approx(r['delta'],rel=1e-13)
        assert r['tau_two']==pytest.approx(3.385989401 if r['q']==.9 else 8.463701164,abs=5e-8)
        assert r['tau_two_relative_error']==pytest.approx(r['tau_two']/r['tau']-1,abs=1e-14)
    for q in CFG['levels']:
        a,b=[next(r for r in PRED['angles'] if r['id']==name and r['q']==q) for name in ['H5','H6']]
        for key in ['beta_exact','beta_two','beta_one']:assert a[key]/b[key]==pytest.approx(math.sqrt(2),rel=1e-14)


def test_independent_prefix_and_physical_convergence():
    for row in json.loads((FOLDER/'leading_points.json').read_text()):
        d=row['diagnostics']
        for key in ['strong_ode','natural_boundary','essential','matrix_residual']:assert d[key]<1e-9
        assert d['mass_gram']<5e-7 and d['energy']<5e-8
        assert row['zeta'][0]>0 and np.diff(row['zeta']).min()>0
    checks=json.loads((FOLDER/'independent_checks.json').read_text())
    for row in checks:
        for i,r in enumerate(row['checks']):
            assert r['count_lower']['count']==i and r['count_upper']['count']==i+1
            if 'frequency_error' in r:assert r['frequency_error']<1e-7
    for r in PRED['convergence']:
        assert max(r['frequency_relative_error'],r['tau_relative_error'])<1e-7
        assert r['shape_mass_error']<1e-6


def test_predictions_sealed_before_full_reference_access(monkeypatch):
    events=json.loads((FOLDER/'access_order.json').read_text())
    assert [r['phase'] for r in events]==['leading_sources_only','predictions_sealed','full_reference_values_parsed']
    assert events[1]['sha256']==study.sha(FOLDER/'apriori_predictions.json')
    assert not PRED['full_data_used']
    s=study.Study(CFG);s.events=[dict(phase='unsealed')]
    def forbidden(*args,**kwargs):raise AssertionError('premature full-reference access')
    monkeypatch.setattr(study,'read_csv',forbidden)
    with pytest.raises(ArithmeticError,match='frozen first'):s.compare_archives([],[])


def test_additive_signed_error_accounting_and_reference_precision():
    rows=study.read_csv(FOLDER/'comparison.csv');parts=study.read_csv(FOLDER/'error_components.csv')
    for r in rows:
        p=next(x for x in parts if x['id']==r['id'] and x['q']==r['q'])
        beta=float(r['beta_full'])
        assert float(r['full_bracket_lower'])<beta<float(r['full_bracket_upper'])
        assert float(r['full_bracket_upper'])-float(r['full_bracket_lower'])==pytest.approx(2e-6,abs=1e-14)
        for method in ['exact','two','one','old']:
            delta=float(r['beta_'+method])-beta
            assert float(r[method+'_absolute_error_deg'])==pytest.approx(abs(delta),abs=1e-14)
            assert float(r[method+'_relative_error_percent'])==pytest.approx(100*delta/beta,abs=1e-12)
        curve,level,asym=[float(p[key+'_percent']) for key in ['curve','level','asymptotic']]
        assert curve<0<level and asym>0
        assert curve+level==pytest.approx(float(r['exact_relative_error_percent']),abs=1e-12)
        assert curve+level+asym==pytest.approx(float(r['two_relative_error_percent']),abs=1e-12)


def test_four_physical_thresholds_and_no_full_or_outer_solves():
    manifest=json.loads((FOLDER/'manifest.json').read_text())
    assert manifest['new_full_model_solves']==manifest['new_outer_solves']==0
    assert manifest['full_thresholds_reused']==4 and len(PRED['angles'])==4
    assert CFG['k_joint']==.1 and CFG['levels']==[.9,.95]
    assert [c['id'] for c in CFG['constructions']]==['H5','H6']
    assert manifest['new_leading_prefixes']==50 and manifest['reused_leading_prefixes']==4


def test_cached_run_cannot_recompute(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('cached result must not solve')
    monkeypatch.setattr(study.Leading,'solve',forbidden)
    monkeypatch.setattr(study,'brentq',forbidden)
    assert study.run()==FOLDER
