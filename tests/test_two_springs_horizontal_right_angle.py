"""Bounded HS-LA5 regressions; consume saved spectra, never solve full model."""
import csv
import json
import math
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_right_angle as study

CFG=json.loads(study.CONFIG.read_text())
FOLDER=study.BASE/study.identity(CFG)[0]
POINTS=json.loads((FOLDER/'point_records.json').read_text())


def csvrows(name):
    with open(FOLDER/name,encoding='utf8') as stream:return list(csv.DictReader(stream))


def test_symbolic_characteristic_and_mass_identity():
    assert study.symbolic_checks()['status']=='PASS'
    b=json.loads((FOLDER/'benchmark.json').read_text())
    assert b['omega']==pytest.approx(.0078805648576462,rel=1e-12)
    assert abs(b['characteristic_residual'])<1e-12
    assert b['zeta']<b['V_first_zeta']
    p=next(r for r in POINTS if r['k_joint']==0 and r['order']==24)
    assert p['zeta'][1]==pytest.approx((math.pi/1.01)**4,rel=1e-12)
    # The first and third positions belong to the cantilever with tip mass.
    for ev in [p['zeta'][0],p['zeta'][2]]:
        z=ev**.25*.99
        residual=1+math.cosh(z)*math.cos(z)+(1.01/.99)*z*(math.sinh(z)*math.cos(z)-math.cosh(z)*math.sin(z))
        assert abs(residual)<1e-10
    assert p['per_mode'][0]['mass_V']<1e-20


def test_exact_space_and_stiffness_mass_congruences():
    rng=np.random.default_rng(527)
    for row in POINTS:
        model=study.RightAngle(row['k_joint'],order=row['order'])
        ids=model.identities()
        assert ids['rank_constraints']==2
        assert max(ids['errors'].values())<CFG['policy']['form_identity_tol']
        assert ids['mass_difference_extreme_eigenvalues'][0]<0<ids['mass_difference_extreme_eigenvalues'][1]
        x=model.Z@rng.normal(size=model.Z.shape[1]);y=model.Z@rng.normal(size=model.Z.shape[1])
        expected=-(model.A@x)*(model.H@y)-(model.A@y)*(model.H@x)
        assert x@(model.Minner-model.M)@y==pytest.approx(expected,abs=1e-13)
        assert row['transformed_inner_residual']<1e-9


def test_first_three_independent_positions_and_no_zero():
    for row in POINTS:
        assert row['diagnostics']['zero_nullity']==0
        assert row['zeta'][0]>0 and min(np.diff(row['zeta']))>0
        if row['order']!=24:continue
        checks=row['independent']
        assert checks['zero_nullity']==checks['count_near_zero']['count']==0
        for i,r in enumerate(checks['roots']):
            assert r['lower']['count']==i and r['upper']['count']==i+1
            assert r['frequency_error']<CFG['policy']['independent_frequency_tol']
            assert r['singular_ratio']<1e-9


def test_mass_accounting_and_variational_bounds():
    for row in csvrows('energies.csv'):
        r={k:float(v) for k,v in row.items() if v and k!='natural_rows'}
        assert r['mass_F']+r['mass_V']+r['mass_rigid_translation']==pytest.approx(1,abs=5e-7)
        assert r['mass_inner_on_outer']==pytest.approx(1+r['mass_cross_inner'],abs=5e-7)
        assert r['outer_mass_on_inner']+r['cross_on_inner']==pytest.approx(r['inner_norm'],abs=5e-7)
        point=next(p for p in POINTS if p['k_joint']==r['k_joint'] and p['order']==r['order'])
        assert point['zeta'][0]<=r['rayleigh_outer_on_inner']*(1+1e-10)
        assert r['inner_zeta']<=r['rayleigh_inner_on_outer']*(1+1e-10)


def test_physical_residuals_and_order_convergence():
    for row in POINTS:
        d=row['diagnostics']
        assert d['stiffness_positive'] and d['mass_positive']
        for key in ['strong','natural','essential','matrix_residual']:assert d[key]<1e-9
        assert d['mass_gram']<5e-7 and d['energy']<5e-8
    for row in csvrows('convergence.csv'):
        assert float(row['frequency_error'])<1e-7
        assert float(row['shape_mass_error'])<1e-6


def test_comparison_normalization_scope_and_preserved_references():
    assert CFG['k_joint']==[0,.1,10] and CFG['orders']==[20,24]
    manifest=json.loads((FOLDER/'manifest.json').read_text())
    assert manifest['new_outer_ritz_prefixes']==6 and manifest['physical_cases']==3
    assert manifest['new_full_model_solves']==0 and manifest['independent_ode_positive_roots']==9
    for r in csvrows('comparison.csv'):
        inner,out,full=[float(r[k]) for k in ['omega_inner','omega_outer','omega_full']]
        assert float(r['inner_error_percent'])==pytest.approx(100*(inner/full-1),abs=1e-12)
        assert float(r['outer_error_percent'])==pytest.approx(100*(out/full-1),abs=1e-12)
        assert abs(out-full)/abs(inner-full)<.02
        assert out==pytest.approx(CFG['epsilon']*math.sqrt(next(p['zeta'][0] for p in POINTS if p['k_joint']==float(r['k_joint']) and p['order']==24)),rel=1e-14)
    exact=json.loads((FOLDER/'exact_support_context.json').read_text())
    assert exact['omega_full_H6']>exact['omega_full_H5']
    assert exact['support_difference_percent']<.00002


def test_archive_verification_never_recomputes(monkeypatch):
    def forbidden(*args,**kwargs):raise AssertionError('cached archive must not solve')
    monkeypatch.setattr(study.RightAngle,'solve',forbidden)
    assert study.run()==FOLDER
