"""Preservation and scope of the small-axial-stiffness preset."""
import json
import numpy as np
import pytest
from scripts.lib import two_springs_external_axial_transition as w


@pytest.mark.parametrize('k',[0.,.01,'RIGID',1e-6,.1])
def test_main_map_rejects_unarchived_endpoint_or_unrequested_stiffness(k):
    e=object.__new__(w.Experiment);e.archive={};e.cfg={'new_chi_values':[1e-5,2.5e-5,1e-4,1e-3]}
    with pytest.raises(ArithmeticError,match='cannot be recomputed'):e.point('J0E4',k,45.)


def test_saved_map_retains_endpoint_spectra_and_valid_suffixes():
    p=w.BASE/'2b5e8e8520c386b7'
    if not (p/'calculation_complete.json').exists():pytest.skip('Map incomplete')
    c=json.loads((p/'contract.json').read_text())['config'];rows=w.read_csv(p/'spectrum.csv');s=w.read_csv(p/'summary.csv')
    by={(r['pair'],w.chi(r['chi_u']),float(r['beta_deg'])):r for r in rows}
    old=w.read_csv(w.ROOT/c['sources']['ranges']['folder']/'spectrum.csv')
    for r in old:
        if w.chi(r['chi_u']) not in c['archived_chi_values']:continue
        v=by[r['pair'],w.chi(r['chi_u']),float(r['beta_deg'])]
        assert all(v[n]==r[n] for n in ['omega1','omega2','omega3','Lambda1','Lambda2','Lambda3','point_sha256','shape_sha256'])
    assert len(s)==12
    for r in s:
        rr=[v for v in rows if (v['pair'],v['chi_u'])==(r['pair'],r['chi_u'])]
        b=np.array([float(v['beta_deg']) for v in rr]);om=np.array([float(v['omega1']) for v in rr])
        for d in (10,5):
            q=w.previous.ranges.high_ranges(b,om,d)
            assert q['status']==r[f'status_{d}']
            assert q['beta_high']==(float(r[f'beta_high_{d}']) if r[f'beta_high_{d}'] else None)
            for run in q['intervals']:
                assert run['entry_resolution_deg'] is None or run['entry_resolution_deg']<=.25
                assert run['exit_resolution_deg'] is None or run['exit_resolution_deg']<=.25
    for r in rows:
        point=json.loads((w.ROOT/r['point_file']).read_text())
        assert point['count']['samples'][0]['count']==0
        assert [a['counts'] for a in point['roots']]==[[0,1],[1,2],[2,3]]
        assert all(all(a['gates'].values()) for a in point['roots'])


def test_axial_and_bending_stiffness_normalizations():
    cfg=json.loads(w.CONFIG.read_text());m=w.physics.Model(**cfg['material'])
    for k,expected in zip(cfg['chi_values'],[0.,.4,1.,4.,40.,400.]):
        ku=k*m.section.EA/m.l
        assert ku*m.l**3/m.section.EI==pytest.approx(expected,rel=1e-14)


def test_addressed_free_limit_and_fixed_angle_brackets():
    p=w.BASE/'2b5e8e8520c386b7'
    if not (p/'addressed_complete.json').exists():pytest.skip('Addressed checks incomplete')
    f=w.read_csv(p/'free_limit_addressed.csv');b=w.read_csv(p/'transition_refined_brackets.csv')
    assert len(f)==len(b)==6
    for r in f:
        assert float(r['omega_free'])<=float(r['omega_1e_6'])<=float(r['free_form_Rayleigh_upper'])*(1+1e-8)
        assert 0<float(r['relative_1e_6'])<float(r['relative_1e_5'])
    for r in b:
        assert float(r['D0_lower'])<0<float(r['D0_upper'])
        assert float(r['log10_width'])<=.125+1e-12
