"""HS-A grid provenance, all-crossings semantics, full spectra and archive reuse."""
import csv
import json
import math
import shutil
from pathlib import Path
import numpy as np
import pytest
from scripts.lib import two_springs_horizontal_angle_map as a

CFG=a.read(a.CONFIG)
FOLDER=a.OUTPUT/'2857f883fb745757'


def table(name):
    with (FOLDER/name).open(encoding='utf8',newline='') as f:return list(csv.DictReader(f))


def test_scaled_grid_radians_and_node_provenance():
    assert CFG['joint_stiffnesses']==[.1,1,10]
    assert CFG['horizontal_kappa']==[4,40,400,4000,20000,'RIGID']
    for k,chi in zip(CFG['horizontal_kappa'],CFG['horizontal_chi']):
        nodes=a.base_nodes(k,CFG)
        assert all(0<=n['beta_deg']<=90 for n in nodes)
        assert len({n['beta_deg'] for n in nodes})==len(nodes)
        assert set(range(0,91,2)).issubset(n['beta_deg'] for n in nodes)
        assert set(np.arange(80,90.25,.5)).issubset(n['beta_deg'] for n in nodes)
        if k!='RIGID':assert chi==pytest.approx(CFG['material']['epsilon']**2*k)
        for n in nodes:
            for origin in n['origins']:
                if origin.startswith('scale_'):
                    tau=float(origin[6:]);b=math.radians(n['beta_deg'])
                    recovered=b*b/(2*.005**2+(0 if k=='RIGID' else 1/k))
                    assert recovered==pytest.approx(tau,rel=2e-15)
    merged=a.merge_nodes([(1.,'one'),(1.+1e-13,'two'),(1.+1e-9,'three')],1e-12)
    assert len(merged)==2 and merged[0]['origins']==['one','two']


def test_all_crossings_and_terminal_tail_not_first_hit():
    angles=[0,10,20,30,40,50,90]
    values=[1,3,1,3,1,3,3]
    brackets=a.crossing_brackets(angles,values,2)
    assert [d for l,r,d in brackets]==['up','down','up','down','up']
    runs,start,status=a.high_segments(angles,values,2)
    assert runs==[[10,10],[30,30],[50,90]] and start==50
    assert status=='TERMINAL_HIGH_RANGE'
    assert a.high_segments(angles,[3,3,1,3,1,3,1],2)[1:]==(None,'NO_TERMINAL_HIGH_RANGE')
    assert a.high_segments(angles,[1]*6+[3],2)[1:]==(90,'TERMINAL_POINT_ONLY')
    assert a.high_segments(angles,[3]*7,2)[1:]==(0,'ALL_CHECKED_NODES_HIGH')


def test_preflight_policy_and_all_frozen_nodes_present():
    p=a.read(FOLDER/'preflight.json');assert p['status']=='PASS' and len(p['checks'])==6
    assert len(p['regressions'])==2
    old=a.read(a.ROOT/'data/input/isotropic_two_springs_horizontal_pilot.json')
    assert CFG['policy']==old['policy']
    assert CFG['frequency_map']['K_plot']==2 and CFG['frequency_map']['K_guard']==3
    assert CFG['frequency_map']['spectrum_semantics']=='sorted_positions'
    assert CFG['frequency_map']['strict_audit_default'] is False
    rows=table('spectra.csv');done={(r['id'],float(r['beta_deg'])) for r in rows}
    grid=a.read(FOLDER/'grid.json');assert len(grid)==18 and sum(map(len,grid.values()))==1335
    for name,nodes in grid.items():
        for n in nodes:assert (name,n['beta_deg']) in done


def test_physical_prefixes_units_and_preservation():
    rows=table('spectra.csv')
    for r in rows:
        rec=a.read(a.ROOT/r['point_file'])
        a.pilot.validate_record(rec,a.h.Case(**rec['case']),CFG['policy'])
        assert int(r['zero_count'])==0
        for i in [1,2,3]:
            assert float(r[f'omega{i}'])==pytest.approx(.005*float(r[f'Lambda{i}'])**2,rel=2e-15)
            assert float(r[f'Omega{i}'])==pytest.approx(float(r[f'Lambda{i}'])**2,rel=2e-15)
    initial=a.read(a.ROOT/CFG['initial_workspace'])
    assert all(a.sha(a.ROOT/p)==s for p,s in initial['archives'].items())
    for p,s in CFG['sources'].items():assert a.sha(a.ROOT/p)==s


def test_thresholds_all_targets_brackets_and_h5_h6():
    thresholds=a.read(FOLDER/'thresholds.json')
    assert len(thresholds)==72
    for r in thresholds:
        assert r['status']=='TERMINAL_HIGH_RANGE'
        assert len(r['crossings'])==1 and not r['high_at_zero']
        root=r['crossings'][0];lo,hi=root['bracket']
        assert root['status']=='PASS' and root['direction']=='up'
        assert hi-lo<=2.000001e-6
        assert root['signs'][0]*root['signs'][1]<0
        assert min(map(abs,root['signs']))>CFG['analysis']['sign_resolution_relative']
        assert abs(root['equation_residual'])<=CFG['analysis']['threshold_residual_tol']
        assert lo<=root['beta_deg']<=hi
        assert r['subsequent_checked_nodes']>1
    for name in {r['id'] for r in thresholds}:
        for q in [.9,.95]:
            t=[r for r in thresholds if r['id']==name and r['q']==q]
            assert t[0]['omega_target']==t[1]['omega_target']
            assert abs(t[0]['terminal_beta_root_deg']-t[1]['terminal_beta_root_deg'])<=2e-6
    assert all(r['status']=='PASS' for r in a.read(FOLDER/'threshold_regression.json'))


def test_final_monotonicity_and_independent_checks():
    m=a.read(FOLDER/'manifest.json');assert m['status']=='PASS'
    assert m['accepted_unique_points']==m['new_unique_points']+m['reused_unique_points']
    assert all(r['status']=='PASS' for r in a.read(FOLDER/'stiffness_monotonicity.json'))
    for r in a.read(FOLDER/'independent.json').values():
        for c in r['checks']:assert c['frequency_error']<=CFG['policy']['ritz_frequency_tol']
    for r in table('summary.csv'):
        assert float(r['beta_max_grid'])==90 and float(r['drop_to_90_percent'])==0
    for name,digest in m['artifacts'].items():assert a.sha(FOLDER/name)==digest
    assert all(r['status']=='PASS' for r in a.read(FOLDER/'limit_checks.json')['cardinal_cases'])


def test_missing_only_reads_accepted_point_without_solver(monkeypatch):
    s=a.Study();c=s.constructions[0]
    before=len(s.entries)
    def forbidden(*args,**kwargs):raise AssertionError('Repeated accepted spectral solve')
    monkeypatch.setattr(s.engine,'solve',forbidden)
    monkeypatch.setattr(a.h,'prefix',forbidden)
    r=s.point(c,0.,'frozen_grid')
    monkeypatch.setattr(a,'brentq',forbidden)
    s.analyze()
    assert r['status']=='PASS' and len(s.entries)==before and s.new_this_run==0


def test_atomic_publication_and_windows_retry(tmp_path,monkeypatch):
    target=tmp_path/'point.json';a.atomic(target,{'status':'old'})
    original=a.os.replace;attempts=[]
    def flaky(source,destination):
        attempts.append(1)
        if len(attempts)==1:raise PermissionError('temporary sharing violation')
        return original(source,destination)
    monkeypatch.setattr(a.os,'replace',flaky)
    a.atomic(target,{'status':'PASS'})
    assert len(attempts)==2 and a.read(target)['status']=='PASS'
    assert not target.with_suffix('.json.tmp').exists()


def test_plot_only_has_no_spectral_execution(tmp_path,monkeypatch):
    from scripts.lib.two_springs_horizontal_angle_figures import render
    from PIL import Image
    def forbidden(*args,**kwargs):raise AssertionError('Spectral solve in plot_only')
    for name in ['prefix','mode','ritz','count_at']:monkeypatch.setattr(a.h,name,forbidden)
    for name in ['spectra.csv','summary.csv','thresholds.csv']:shutil.copy2(FOLDER/name,tmp_path/name)
    render(tmp_path)
    assert a.read(tmp_path/'figures_manifest.json')['new_eigenvalue_solves']==0
    assert len(list((tmp_path/'figures').glob('*.pdf')))==3
    for p in (tmp_path/'figures').glob('*.png'):
        with Image.open(p) as im:assert min(im.info['dpi'])>=300
