"""HS-B geometry identity, nondimensional experiment and accepted prefix tests."""
import csv
import math
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_geometry_map as b

CFG=b.read(b.CONFIG)
FOLDER=b.OUTPUT/'7b2b5ba8c33abe48'


def table(name):
    with (FOLDER/name).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))


def test_exact_sparse_variant_contract():
    vs=b.variants(CFG)
    assert len(vs)==34
    assert {g:sum(v['group']==g for v in vs) for g in ['B0','B1','B2','B3','B4']}==dict(B0=4,B1=8,B2=16,B3=4,B4=2)
    for v in vs:
        m,g=b.geometry(v,CFG)
        assert g['status']=='PASS' and g['limiter']<=.1
        assert g['L1']+g['L2']==pytest.approx(2)
        assert g['I_y']/g['A']==pytest.approx(v['epsilon']**2)
        assert m.J==pytest.approx(m.rho*m.nu**2*g['I_y']*(1+v['aspect']**2),rel=2e-15)
        assert g['k_joint_dimensional']*m.l/g['EI']==pytest.approx(v['k_joint'])
        if v['chi_s']!='RIGID':
            assert v['kappa_horizontal']*v['epsilon']**2==pytest.approx(v['chi_s'])
            assert g['k_horizontal_dimensional']*m.l/g['EA']==pytest.approx(v['chi_s'])
    assert {v['mu'] for v in vs}=={-.1,-.01,0,.01,.1}


def test_aspect_only_changes_gradient_inertia_in_normalized_operator():
    v=next(v for v in b.variants(CFG) if v['id']=='B0_A1')
    m4,g4=b.geometry(v,CFG);m2,g2=b.geometry({**v,'aspect':2},CFG)
    for key in ['EA','KGA','rhoA','rhoI']:
        assert g4[key]/g4['EI']==pytest.approx(g2[key]/g2['EI'],rel=3e-15)
    assert (m4.J/g4['EI'])/(m2.J/g2['EI'])==pytest.approx(17/5)
    # Identical bending ratios do not authorize sharing a different arm catalog.
    assert b.arm_signature(m4,.99)!=b.arm_signature(m2,.99)
    assert b.arm_signature(m4,.99)!=b.arm_signature(m4,1.01)


def test_grid_scale_and_origin_for_every_geometry():
    grids=b.read(FOLDER/'grid.json')
    for v in b.variants(CFG):
        nodes=grids[v['id']];angles={n['beta_deg'] for n in nodes}
        assert set(range(0,91,2))<=angles
        assert set(np.arange(80,90.25,.5))<=angles
        scale2=2*v['epsilon']**2+(0 if v['chi_s']=='RIGID' else 1/v['kappa_horizontal'])
        for n in nodes:
            for origin in n['origins']:
                if origin.startswith('scale_'):
                    assert math.radians(n['beta_deg'])**2/scale2==pytest.approx(float(origin[6:]),rel=4e-15)


def test_complete_physical_prefixes_and_units():
    rows=table('spectra.csv');grid=b.read(FOLDER/'grid.json');done={(r['id'],float(r['beta_deg'])) for r in rows}
    for name,nodes in grid.items():assert all((name,n['beta_deg']) in done for n in nodes)
    for r in rows:
        rec=b.read(b.ROOT/r['point_file'])
        b.a.pilot.validate_record(rec,b.a.h.Case(**rec['case']),CFG['policy'])
        assert int(r['zero_count'])==0
        for i in [1,2,3]:
            assert float(r[f'omega{i}'])==pytest.approx(float(r['epsilon'])*float(r[f'Lambda{i}'])**2,rel=3e-15)
    assert len({r['id'] for r in rows})==34


def test_catalog_physical_signatures_and_reuse():
    signatures={}
    for p in (FOLDER/'variants').glob('*/catalog_identity.json'):
        for r in b.read(p).values():
            assert b.sha(b.ROOT/r['file'])==r['sha256']
            cat=b.read(b.ROOT/r['file'])
            assert cat['signature']==r['signature'] and cat['status']=='PASS'
            assert cat['catalog']['count_bound']['upper_count']==len(cat['catalog']['poles'])
            if r['file'] in signatures:assert signatures[r['file']]==r['signature']
            signatures[r['file']]=r['signature']


def test_roots_both_targets_direct_brackets():
    rows=b.read(FOLDER/'thresholds.json');assert len(rows)==136
    for r in rows:
        assert r['status']=='TERMINAL_HIGH_RANGE' and len(r['crossings'])==1
        root=r['crossings'][0]
        assert root['direction']=='up' and root['status']=='PASS'
        assert root['signs'][0]*root['signs'][1]<0
        assert min(map(abs,root['signs']))>CFG['analysis']['sign_resolution_relative']
        assert root['bracket'][1]-root['bracket'][0]<=2.000001e-6
        assert abs(root['equation_residual'])<=CFG['analysis']['threshold_residual_tol']


def test_combined_unreachable_target_is_a_status_not_an_angle():
    rows=table('predictions.csv');assert len(rows)==136
    absent=[r for r in rows if r['combined_status']=='TARGET_AT_OR_ABOVE_LEADING_SATURATION']
    assert absent
    for r in absent:assert r['combined_beta_deg']=='' and r['combined_error_percent']==''
    for r in rows:
        assert float(r['full_xi'])==pytest.approx(math.radians(float(r['full_beta_deg']))/float(r['epsilon']))
        assert r['qualification']=='FORMAL_EXTENSION_POINTWISE_CHECKED_NO_UNIFORM_ERROR_BOUND'


def test_missing_only_point_does_not_resolve(monkeypatch):
    s=b.Study();child=s.children['B1_A1_e0.0025']
    def forbidden(*args,**kwargs):raise AssertionError('Accepted spectrum was recomputed')
    monkeypatch.setattr(child.engine,'solve',forbidden)
    monkeypatch.setattr(b.a.h,'prefix',forbidden)
    assert child.point(child.v,0.,'frozen_grid')['status']=='PASS'
    assert child.new_this_run==0


def test_no_threshold_work_after_detected_interior_maximum(monkeypatch):
    s=b.GeometryStudy.__new__(b.GeometryStudy);s.v={'id':'synthetic'};s.cfg=CFG
    pairs=[({'beta_deg':x},{'omega':[w,2,3]}) for x,w in [(0,1),(45,1.5),(90,1.4)]]
    monkeypatch.setattr(s,'series',lambda c:pairs);calls=[]
    monkeypatch.setattr(s,'extrema',lambda c:calls.append('diagnostic'))
    monkeypatch.setattr(s,'stop',lambda *args:calls.append(args[-1]))
    monkeypatch.setattr(s,'thresholds',lambda *args:pytest.fail('Threshold work after stop'))
    s.investigate();assert calls==['diagnostic','STOP_INTERIOR_MAXIMUM']
