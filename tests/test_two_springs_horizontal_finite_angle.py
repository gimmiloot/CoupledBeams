"""HS-C analytical identities, exact limits, independent count and provenance."""
import math
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_finite_angle as f

CFG=f.read(f.CONFIG)
FOLDER=f.OUTPUT/'6b15c780c138c0ed'


def test_kinematics_and_mass_bilinear_identity():
    rng=np.random.default_rng(2857)
    for beta in [0,2,30,60,90]:
        c,s=f.trig(beta)
        for A,B,V,Y in rng.normal(size=(10,4)):
            U=s*A;Wend=-c*A
            assert c*U+s*Wend==pytest.approx(0,abs=1e-15)
            assert A-s*U+c*Wend==pytest.approx(0,abs=2e-15)
            assert c*U+s*(V-c*A)==pytest.approx(s*V,abs=2e-15)
            assert s*s*A*B+(V-c*A)*(Y-c*B)==pytest.approx(V*Y-c*(A*Y+B*V)+A*B,abs=3e-15)


@pytest.mark.parametrize('mu',[-.1,0,.01,.1])
def test_exact_zero_and_right_angle_limits(mu):
    old=f.Leading(.1,mu,20);outer=f.RightAngle(.1,mu,20)
    zero=f.FiniteAngle(.1,4,0,mu,20).solve()
    zeroplus=f.FiniteAngle(.1,'RIGID',0,mu,20).solve()
    right=f.FiniteAngle(.1,'RIGID',90,mu,20)
    assert np.allclose(zero['zeta'],old.solve(0)['zeta'],rtol=1e-11,atol=0)
    assert np.allclose(zeroplus['zeta'],old.solve('RIGID')['zeta'],rtol=1e-11,atol=0)
    assert zeroplus['zeta'][0]>zero['zeta'][0]
    assert np.array_equal(right.K,outer.K) and np.array_equal(right.M,outer.M) and np.array_equal(right.C,outer.C)
    finite=f.FiniteAngle(.1,4,90,mu,20)
    assert finite.Z.shape[1]==right.Z.shape[1]+1
    assert finite.solve()['zeta'][0]<right.solve()['zeta'][0]


def test_independent_boundary_and_count_without_full_solver(monkeypatch):
    from scripts.lib import isotropic_two_springs_horizontal_spring as full
    def forbidden(*args,**kwargs):raise AssertionError('Full model calls prohibited')
    monkeypatch.setattr(full,'prefix',forbidden)
    model=f.FiniteAngle(1,400,30,.1,24);sol=model.solve()
    checks=f.independent(sol['zeta'],1,400,30,.1,CFG['policy'])
    for i,row in enumerate(checks):
        assert row['lower']['count']==i and row['upper']['count']==i+1
        assert row['frequency_error']<CFG['policy']['independent_frequency_tol']
    d=sol['diagnostics']
    assert d['natural']<1e-9 and d['strong']<1e-9 and d['mass_square_identity_error']<1e-12


def test_fixed_chi_family_is_not_fixed_kappa():
    cases=[r for r in CFG['cases'] if r['id'] in ['B0_A2','B1_A2_e0.0025','B1_A2_e0.006']]
    assert len(cases)==3 and len({c['kappa_s'] for c in cases})==3
    for c in cases:assert c['kappa_s']*c['epsilon']**2==pytest.approx(.0001)


def test_accepted_low_prefixes_and_source_preservation():
    records=f.read(FOLDER/'points.json');assert records
    for p in records.values():
        assert p['status']=='PASS' and np.all(np.diff(p['zeta'])>0)
        assert p['frequency_convergence']<=1e-7 and p['shape_convergence']<=1e-6
        for i,row in enumerate(p['orders'][-1]['independent']):
            assert row['lower']['count']==i and row['upper']['count']==i+1
    for path,digest in CFG['sources'].items():assert f.sha(f.ROOT/path)==digest


def test_conditional_thresholds_and_zero_plus_semantics():
    rows=f.read(FOLDER/'apriori_predictions.json');assert len(rows)==60
    zero=[r for r in rows if r['status']=='ALREADY_HIGH_AT_ZERO_PLUS_NOT_FULL_THRESHOLD']
    assert zero and all(r['beta_deg'] is None for r in zero)
    skipped=[r for r in rows if r['status']=='SKIPPED_LATE_ACCURACY_GATE']
    assert skipped and all(r['beta_deg'] is None for r in skipped)
    for r in rows:
        if r['status']=='PASS':
            assert len(r['crossings'])==1
            c=r['crossings'][0];assert c['signs'][0]*c['signs'][1]<0 and abs(c['residual'])<=1e-9


def test_predictions_not_calibrated_by_full_thresholds():
    p=f.read(FOLDER/'prediction_provenance.json')
    assert p['predictions_written_before_threshold_table_read'] and p['fit_parameters']==0
    assert p['prior_HSB_threshold_knowledge']
    assert p['predictions_sha256']==f.sha(FOLDER/'apriori_predictions.json')
    m=f.read(FOLDER/'manifest.json');assert m['new_full_solves']==0 and m['status']=='PASS'
    assert m['reused_full_common_points']==90 and m['reused_full_thresholds']==18
