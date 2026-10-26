"""Engineering range semantics: terminal suffix, re-entry and omega scaling."""
import numpy as np
import pytest
from scripts.lib.two_springs_high_frequency_angle_ranges import high_ranges


def test_terminal_uses_last_reentry_not_first_crossing():
    r=high_ranges([0,10,20,30,40,90],[1,2,1.8,1.96,1.99,2],2)
    assert r['first_high']==10 and r['beta_high']==30
    assert r['terminal_nodes']==3 and len(r['intervals'])==2
    assert r['intervals'][0]['singleton'] and not r['intervals'][0]['terminal']


def test_no_terminal_range_retains_separate_intervals():
    r=high_ranges([0,10,20,30,90],[1,2,2,1.9,1.7],2)
    assert r['status']=='NO_TERMINAL_HIGH_RANGE' and r['beta_high'] is None
    assert r['intervals'][0]['start_deg']==10 and r['intervals'][0]['end_deg']==20
    assert r['beta_max']==10 and r['R90']==pytest.approx(15)


def test_only_endpoint_is_not_a_wide_range():
    r=high_ranges([0,20,80,90],[1,1.5,1.8,2],2)
    assert r['status']=='TERMINAL_POINT_ONLY' and r['beta_high']==90
    assert r['terminal_nodes']==1 and r['resolution']==10


def test_small_total_gain_can_make_every_node_high():
    r=high_ranges([0,13,57,90],[1,1.004,1.009,1.003],2)
    assert r['beta_high']==0 and r['G_max_percent']==pytest.approx(.9)
    assert r['resolution'] is None and r['beta_max']==57


def test_percent_is_circular_frequency_not_lambda():
    lam=np.array([.8,.98,1.]);r=high_ranges([0,10,90],lam**2,2)
    assert r['R'][1]==pytest.approx(3.96) and r['beta_high']==90


def test_inclusive_exact_threshold_and_nonuniform_resolution():
    r=high_ranges([0,8.75,9,90],[1,1.8,1.9,2],5)
    assert r['beta_high']==9 and r['resolution']==.25


def test_initial_90_percent_does_not_skip_a_later_dip():
    q=high_ranges([0,1,10,90],[.91,.85,.98,1],10)
    assert q['first_high']==0 and q['beta_high']==10
    assert len(q['intervals'])==2


def test_90_percent_can_start_at_zero_while_95_percent_does_not():
    b=[0,.25,.5,90];w=[.92,.94,.96,1]
    q90=high_ranges(b,w,10);q95=high_ranges(b,w,5)
    assert q90['beta_high']==0 and q90['resolution'] is None
    assert q95['beta_high']==.5 and q95['resolution']==.25


def test_revised_preset_cannot_solve_an_unapproved_point():
    from scripts.lib import two_springs_high_frequency_angle_ranges as w
    exp=object.__new__(w.Experiment)
    exp.revised=True;exp.cfg={'map':{'parameter_grid':[40]}};exp.archive={};exp.allowed={}
    with pytest.raises(ArithmeticError,match='approved'):
        exp.point('I_J01',1.,40)


def test_revised_archive_preserves_old_points_and_checks_whole_tail():
    import json
    from scripts.lib import two_springs_high_frequency_angle_ranges as w
    folder=w.REVISED_BASE/'b7c2787461cb3cda';prior=w.BASE/'18570e576270e0b4'
    if not (folder/'calculation_complete.json').exists():pytest.skip('Revised processing not yet complete')
    rows=w.read_csv(folder/'spectrum.csv');old=w.read_csv(prior/'spectrum.csv');summary=w.read_csv(folder/'summary.csv')
    key=lambda r:(r['pair'],w.stiffness(r['k_w']),float(r['beta_deg']))
    by={key(r):r for r in rows};before={key(r):r for r in old}
    for k,r in before.items():
        assert all(by[k][n]==r[n] for n in ['omega1','omega2','omega3','Lambda1','Lambda2','Lambda3','point_sha256','shape_sha256'])
    plan=json.loads((folder/'refinement_plan.json').read_text())
    for k in by.keys()-before.keys():
        assert any((r['pair'],r['k_w'])==k[:2] and k[2] in r['new_beta_deg'] and r['reason']=='crossing_10_percent' for r in plan)
    assert len(summary)==20
    for r in summary:
        group=sorted([v for v in rows if (v['pair'],v['k_w'])==(r['pair'],r['k_w'])],key=lambda v:float(v['beta_deg']))
        assert float(r['beta_high_10'])<=float(r['beta_high_5'])
        for d in (10,5):
            start=float(r[f'beta_high_{d}']);threshold=float(r['omega_max'])*(1-d/100)
            assert all(float(v['omega1'])>=threshold for v in group if float(v['beta_deg'])>=start)
            if start>0:
                outside=[v for v in group if float(v['beta_deg'])<start][-1]
                assert float(outside['omega1'])<threshold and start-float(outside['beta_deg'])<=.25+1e-12
            else:assert r[f'resolution_{d}']==''


@pytest.mark.parametrize('b',[[0,10,20],[0,10,10,90],[10,20,90]])
def test_invalid_coverage_is_rejected(b):
    with pytest.raises(ArithmeticError):high_ranges(b,np.ones(len(b)),2)


def test_saved_grid_ranges_provenance_and_refinement_authorization():
    import json
    from scripts.lib import two_springs_high_frequency_angle_ranges as w
    folder=w.BASE/'18570e576270e0b4'
    if not (folder/'calculation_complete.json').exists():pytest.skip('The bounded study is not complete')
    rows=w.read_csv(folder/'spectrum.csv');summary=w.read_csv(folder/'summary.csv');plan=json.loads((folder/'refinement_plan.json').read_text())
    assert len(summary)==20
    for r in summary:
        group=[v for v in rows if (v['pair'],v['k_w'])==(r['pair'],r['k_w'])]
        for level in (2,5):
            q=high_ranges([float(v['beta_deg']) for v in group],[float(v['omega1']) for v in group],level)
            assert q['status']==r[f'status_{level}']
            assert q['beta_high']==(float(r[f'beta_high_{level}']) if r[f'beta_high_{level}'] else None)
            assert q['resolution'] is None or q['resolution']<=.25
    for r in rows:
        rec=json.loads((w.ROOT/r['point_file']).read_text());beta=float(r['beta_deg'])
        assert rec['case']['beta']==beta and rec['status']=='PASS'
        assert [x['counts'] for x in rec['roots']]==[[0,1],[1,2],[2,3]]
        assert float(r['omega1'])==rec['roots'][0]['omega']
        if r['origin']=='new_study' and beta<=20:
            assert any(p['pair']==r['pair'] and p['k_w']==w.stiffness(r['k_w']) and beta in p['new_beta_deg'] for p in plan)
