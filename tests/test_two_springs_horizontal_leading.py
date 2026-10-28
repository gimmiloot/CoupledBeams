"""HS-L1: independent variational rows, exact saturation and accepted data."""
import csv
import json
import math
from pathlib import Path
import numpy as np
import pytest
from scripts.lib.isotropic_two_springs_horizontal_leading import Leading,independent_roots
from scripts.lib import two_springs_horizontal_leading_validation as workflow


@pytest.mark.parametrize('tau',[0.,3.,'RIGID'])
def test_variational_boundary_signs(tau):
    m=Leading(order=6)
    C=np.vstack((m.constraint,m.D)) if tau=='RIGID' else m.constraint
    from scipy.linalg import null_space
    Z=null_space(C);rng=np.random.default_rng(42);c=Z@rng.normal(size=Z.shape[1]);v=Z@rng.normal(size=Z.shape[1])
    t=0 if tau=='RIGID' else tau
    lhs=float(c@(m.K0+t*m.P)@v)
    volume=sum(float(m.weights@((f[4]@c)*(f[0]@v))*L) for f,L in zip(m.values,m.lengths))
    e=np.array([a@c for a in m.ends]);o=np.array([a@c for a in m.outer]);ev=np.array([a@v for a in m.ends]);ov=np.array([a@v for a in m.outer])
    jump=m.jump@c
    boundary=-o[1,2]*ov[1,1]+(o[1,3]+t*(m.D@c))*ov[1,0]
    boundary+=(e[0,2]+m.k_joint*jump)*ev[0,1]+(e[1,2]-m.k_joint*jump)*ev[1,1]
    boundary-=(e[0,3]-e[1,3]-t*(m.D@c))*ev[0,0]
    assert lhs==pytest.approx(volume+boundary,abs=3e-12)


def test_free_trace_rank_one_and_exact_limit():
    m=Leading();assert np.linalg.matrix_rank(m.P)==1
    assert np.linalg.norm(m.outer[1][0]@m.Z)>0
    assert np.linalg.norm(m.outer[1][1]@m.Z)>0
    s=m.solve('RIGID');c=s['coefficients'][:,0]
    assert abs(m.D@c)<1e-12
    # The exact connection does not fix the two transverse displacements to zero.
    assert abs(m.outer[1][0]@c)>.01
    assert abs(s['per_mode'][0]['reaction'])>.01
    assert s['diagnostics']['natural_boundary']<1e-9


@pytest.mark.parametrize('tau',[0,2.9124184246057174,6.135887525019995,'RIGID'])
def test_independent_ode_and_order_convergence(tau):
    a=Leading(order=20).solve(tau);b=Leading(order=24).solve(tau)
    assert max(abs(np.sqrt(a['zeta']/b['zeta'])-1))<1e-7
    assert b['diagnostics']['strong_ode']<1e-9
    assert b['diagnostics']['energy']<5e-8
    assert max(abs(r['relative_frequency_error']) for r in independent_roots(b['zeta'],tau))<1e-7


def test_monotone_bounded_spectrum_and_coefficients():
    m=Leading();values=np.array([m.solve(t)['zeta'] for t in [0,1,10,100,10000,'RIGID']])
    assert np.min(np.diff(values,axis=0))>0
    c=m.solve(0)['coefficients'][:,0];a=.25*(m.D@c)**2;b=.5*(m.D@c)**2
    assert a/b==.5
    assert a==pytest.approx(.5872162446895,rel=1e-10)


def archive():
    cfg=json.loads(workflow.CONFIG.read_text());identifier,_=workflow.identity(cfg)
    return workflow.BASE/identifier


def test_threshold_brackets_levels_and_no_global_map():
    folder=archive();m=json.loads((folder/'manifest.json').read_text())
    assert m['status']=='PASS' and m['full_new_points']==44 and m['full_new_this_invocation']==0
    rows=list(csv.DictReader(open(folder/'thresholds.csv')))
    assert len(rows)==4
    for r in rows:
        assert float(r['bracket_F_lower'])<0<float(r['bracket_F_upper'])
        assert float(r['beta_full_bracket_upper'])-float(r['beta_full_bracket_lower'])<2.01e-6
        assert abs(float(r['equation_residual']))<2e-8
        assert float(r['beta_lead_deg'])<float(r['beta_full_deg'])
    for q in ['0.9','0.95']:
        a,b=[next(r for r in rows if r['id']==id and r['q']==q) for id in ['H5','H6']]
        assert float(a['beta_lead_deg'])/float(b['beta_lead_deg'])==pytest.approx(math.sqrt(2),rel=1e-14)
    new=[r for r in csv.DictReader(open(folder/'full_points.csv')) if r['origin']!='HS-P1']
    assert max(float(r['beta_deg']) for r in new)<2


def test_plot_only_has_no_spectral_calculations(tmp_path,monkeypatch):
    folder=archive()
    for name in ['comparison.csv','thresholds.csv']:(tmp_path/name).write_bytes((folder/name).read_bytes())
    def forbidden(*a,**k):raise AssertionError('plot_only spectral call')
    monkeypatch.setattr(Leading,'__post_init__',forbidden)
    monkeypatch.setattr(workflow.full.Run,'solve',forbidden)
    workflow.plot_only(tmp_path)
    assert len(list((tmp_path/'figures').glob('*.pdf')))==2
    from PIL import Image
    for p in (tmp_path/'figures').glob('*.png'):
        with Image.open(p) as image:assert min(image.info['dpi'])>=300
