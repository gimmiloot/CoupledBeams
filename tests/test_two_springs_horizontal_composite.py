"""HS-D1 checks: mass derivative, positivity, archived predictions and scope."""
import math
from pathlib import Path
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_composite as c

FOLDER=c.output_folder()


@pytest.mark.parametrize('k,mu',[(.1,.01),(1,.01),(1,.1),(10,.01),(10,-.1)])
def test_C2_mass_derivative_and_common_limit(k,mu):
    model=c.outer.FiniteAngle(k,'RIGID',0,mu,20)
    sol=model.solve();x=sol['coefficients'][:,0];z=sol['zeta'][0]
    coefficient=-z*float(model.A@x)*float(model.H@x)
    assert coefficient==pytest.approx(z*float(x@model.cross@x)/2,rel=2e-13)
    h=math.radians(.1)
    changed=c.outer.FiniteAngle(k,'RIGID',.1,mu,20)
    assert np.array_equal(model.K,changed.K)
    assert np.allclose(changed.M-model.M,(math.cos(h)-1)*model.cross,atol=2e-16,rtol=0)
    assert z==pytest.approx(c.inner.Leading(k,mu,20).solve('RIGID')['zeta'][0],rel=2e-12)
    expected=-z*float(x@(changed.M-model.M)@x)/h**2
    assert expected==pytest.approx(coefficient,rel=1e-6)


def test_positive_composite_and_no_control_composition():
    rows=c.outer.rows_csv(FOLDER/'model_curves.csv')
    coefficients={r['id']:r for r in c.read(FOLDER/'coefficients.json')}
    for r in rows:
        if r['id']=='B0_A2':
            assert r['zeta_comp']=='' and r['primary']=='C_I'
            assert float(r['zeta_primary'])==float(r['zeta_C'])
            continue
        z=float(r['zeta_comp']);bound=coefficients[r['id']]['orders'][-1]['global_composite_lower_bound']
        assert z>=bound*(1-1e-11)>0
        assert float(r['zeta_L1'])<=z+1e-11 and z<=float(r['zeta_C'])+1e-11
        if float(r['beta_deg'])==0:
            assert r['zero_C_semantics']=='right_limit_0plus'
            assert z==pytest.approx(float(r['zeta_L1']),abs=2e-12)


def test_coefficients_unfitted_and_two_orders():
    rows=c.read(FOLDER/'coefficients.json')
    for r in rows:
        assert abs(r['C2_Richardson_relative_error'])<5e-5
        for p in r['orders']:
            assert p['C2']>0 and abs(p['linear']['G'])<1e-9
            assert p['a']==pytest.approx(p['linear']['a'],rel=1e-8)
            assert p['reaction']['kkt_residual']<1e-9
        for key in ['a','b','C2']:assert r['orders'][0][key]==pytest.approx(r['orders'][1][key],rel=1e-7)
    base=next(r for r in rows if r['id']=='B0_A1')['orders'][-1]
    assert base['b']==pytest.approx(-.639617945455777,rel=1e-12)


def test_finite_tau_not_replaced_and_prefix_gates():
    maximum=0.
    for p in (FOLDER/'points').glob('*.json'):
        r=c.read(p)
        assert r['status']=='PASS' and r['frequency_convergence']<=1e-7
        for row in r['orders']:
            assert all(z>0 for z in row['zeta']) and np.all(np.diff(row['zeta'])>0)
        if r['kind']=='L1' and r['parameter']!='RIGID':maximum=max(maximum,r['parameter'])
    assert maximum==pytest.approx(49348.02200544679)
    rows=c.read(FOLDER/'preflight.json')['values']
    for r in rows:
        if r['beta_deg']==90:assert r['zeta_L1']<r['zeta_infty']-1e-5


def test_roots_and_slopes_are_direct_and_bracketed():
    rows=c.read(FOLDER/'predictions.json');assert len(rows)==12
    for r in rows:
        assert r['status']=='PASS' and r['count_up']==1 and r['count_down']==0
        for root in r['crossings']:
            assert root['signs'][0]*root['signs'][1]<0
            assert root['bracket'][1]-root['bracket'][0]==pytest.approx(2e-6,rel=1e-8)
            assert abs(root['residual'])<=1e-9
            assert abs(root['differences'][-1]['relative_error'])<=5e-5
    provenance=c.read(FOLDER/'prediction_provenance.json')
    assert provenance['fit_parameters']==0 and provenance['prior_HSB_knowledge']
    assert provenance['predictions_sha256']==c.sha(FOLDER/'predictions.json')


def test_scope_and_frozen_sources():
    cfg=c.read(c.CONFIG);assert len(cfg['cases'])==6 and cfg['accuracy_requirement_percent'] is None
    assert cfg['new_full_solves_allowed']==0
    for path,digest in cfg['sources'].items():assert c.sha(c.ROOT/path)==digest
    m=c.read(FOLDER/'manifest.json');assert m['new_full_solves']==0 and m['full_archived_grid_points']==524


def test_plot_only_has_no_solver_calls(monkeypatch,tmp_path):
    for name in ['comparison.csv','interval_errors.csv','threshold_comparison.csv']:
        (tmp_path/name).write_bytes((FOLDER/name).read_bytes())
    def forbidden(*args,**kwargs):raise AssertionError('Spectral calculation in plot_only')
    monkeypatch.setattr(c.inner.Leading,'solve',forbidden)
    monkeypatch.setattr(c.outer.FiniteAngle,'solve',forbidden)
    c.render(tmp_path)
    assert c.read(tmp_path/'figures_manifest.json')['new_spectral_solves']==0
    assert len(list((tmp_path/'figures').glob('*.pdf')))==4
