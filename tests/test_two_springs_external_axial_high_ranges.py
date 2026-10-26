"""Bounded axial extension: support identity, archive preservation and suffixes."""
import json
import numpy as np
import pytest
from scripts.lib import two_springs_external_axial_high_ranges as w


def test_transverse_compliance_is_not_an_axial_reference():
    c=dict(beta=0.,mu=.01,k_joint=.1,k_end=10.,variant='RL_FULL',k_w=1.)
    with pytest.raises(ArithmeticError,match='not axial'):w.archive_case(c)
    c['k_w']='RIGID'
    assert w.archive_case(c).k_w=='RIGID'


def test_exact_axial_support_cannot_be_recomputed():
    e=object.__new__(w.Experiment);e.archive={}
    with pytest.raises(ArithmeticError,match='rigid map'):e.point('J0E4','RIGID',22.)


def test_unplanned_finite_point_cannot_be_solved():
    e=object.__new__(w.Experiment);e.archive={};e.allowed={};e.cfg={'map':{'parameter_grid':[22.]}}
    with pytest.raises(ArithmeticError,match='Unplanned'):e.point('J0E4',.1,21.1)


def test_integer_zero_and_float_zero_have_one_cache_contract():
    from dataclasses import asdict
    e=object.__new__(w.Experiment);e.cfg={'mu':.01,'cases':{'J0E4':[.1,10.]}}
    assert json.dumps(asdict(e.case('J0E4',0,40)),sort_keys=True)==json.dumps(asdict(e.case('J0E4',0.,40.)),sort_keys=True)


def test_completed_extension_preserves_sources_and_entire_terminal_suffix():
    folder=w.BASE/'242dfa8318f411c5'
    if not (folder/'calculation_complete.json').exists():pytest.skip('Bounded extension not complete')
    rows=w.read_csv(folder/'spectrum.csv');summary=w.read_csv(folder/'summary.csv')
    cfg=json.loads((folder/'contract.json').read_text())['config'];plan=json.loads((folder/'refinement_plan.json').read_text())
    by={(r['pair'],w.chi(r['chi_u']),float(r['beta_deg'])):r for r in rows}
    prior=w.read_csv(w.ROOT/cfg['sources']['axial']['folder']/'spectrum.csv')
    for r in prior:
        v=by[r['pair'],w.chi(r['chi_u']),float(r['beta_deg'])]
        assert all(v[n]==r[n] for n in ('omega1','omega2','omega3','Lambda1','Lambda2','Lambda3','point_sha256','shape_sha256'))
    assert len(summary)==10
    for r in summary:
        rr=[v for v in rows if (v['pair'],v['chi_u'])==(r['pair'],r['chi_u'])]
        b=np.array([float(v['beta_deg']) for v in rr]);om=np.array([float(v['omega1']) for v in rr])
        for d in (10,5):
            q=w.ranges.high_ranges(b,om,d)
            assert q['status']==r[f'status_{d}']
            assert q['beta_high']==(float(r[f'beta_high_{d}']) if r[f'beta_high_{d}'] else None)
            if q['beta_high'] is not None:
                assert np.all(om[b>=q['beta_high']]>=om.max()*(1-d/100))
                assert q['resolution'] is None or q['resolution']<=.25
    for r in rows:
        rec=json.loads((w.ROOT/r['point_file']).read_text())
        assert [v['counts'] for v in rec['roots']]==[[0,1],[1,2],[2,3]]
        assert all(all(v['gates'].values()) for v in rec['roots'])
        if r['origin']=='new_study':
            assert r['chi_u']!='RIGID'
            b=float(r['beta_deg'])
            if b not in cfg['map']['parameter_grid']:
                assert any(p['pair']==r['pair'] and p['chi_u']==w.chi(r['chi_u']) and b in p['new_beta_deg'] for p in plan)
