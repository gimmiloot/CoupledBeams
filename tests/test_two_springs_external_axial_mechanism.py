"""Bounded transition preset and its scientific acceptance distinctions."""
import hashlib
import json
import numpy as np
import pytest
from scripts.lib import two_springs_external_axial_mechanism as w

FOLDER=w.BASE/'45c7729f346df3c2'


@pytest.mark.parametrize('k',[0.,1e-4,1e-3,.000283,'RIGID'])
def test_unplanned_new_map_is_rejected(k):
    e=object.__new__(w.Experiment);e.archive={};e.cfg={'new_chi_values':[.00024,.00027,.0003]}
    with pytest.raises(ArithmeticError,match='three-stiffness scope'):e.point('J0E4',k,45.)


def test_frozen_theory_and_physics():
    c=json.loads(w.CONFIG.read_text())
    for path,digest in c['frozen_files'].items():assert hashlib.sha256((w.ROOT/path).read_bytes()).hexdigest()==digest


def test_source_points_and_previous_sign_brackets_are_exactly_reused(tmp_path):
    if not FOLDER.exists():pytest.skip('Local research archive absent')
    rows=w.read_csv(FOLDER/'reference_spectrum.csv')
    assert len(rows)==638
    for r in rows:
        p=w.ROOT/r['point_file'];rec=json.loads(p.read_text())
        assert hashlib.sha256(p.read_bytes()).hexdigest()==r['point_sha256']
        assert all(float(r[f'omega{i}'])==v['omega'] for i,v in enumerate(rec['roots'],1))
    assert len(w.read_csv(FOLDER/'archive_bracket_regression.csv'))==6
    # A clean continuation reproduces the bounded scientific trigger without
    # depending on an unrecorded manual plan or requesting spectral work.
    for name in ('summary.csv','curvature_summary.csv'):
        (tmp_path/name).write_bytes((FOLDER/name).read_bytes())
    e=object.__new__(w.Experiment);e.cfg=json.loads(w.CONFIG.read_text());e.folder=tmp_path;e.a=w.api()
    summary=w.read_csv(FOLDER/'summary.csv')
    for r in summary:
        r['chi_u']=float(r['chi_u']);r['beta_max']=float(r['beta_max'])
    e.register_diagnostic_plan(summary)
    plan=json.loads((tmp_path/'diagnostic_stiffness_plan.json').read_text())
    assert plan['defined_before_extra_maps'] and plan['chi_values']==[.000276,.000283] and plan['no_further_expansion']


def test_reflection_and_full_model_curvature_are_distinct_from_leading():
    if not (FOLDER/'curvature_summary.csv').exists():pytest.skip('Local analysis absent')
    ref=w.read_csv(FOLDER/'reflection_checks.csv')
    assert len(ref)==30
    assert max(float(r[k]) for r in ref for k in ('joint_error','external_error','ode_error'))<1e-13
    curv=w.read_csv(FOLDER/'curvature_summary.csv')
    assert len(curv)==6 and all(r['status']!='UNRESOLVED' for r in curv)
    for r in curv:
        if float(r['chi_u'])<.0003:assert float(r['C_upper'])<0
        else:assert float(r['C_lower'])>0
    optional=w.read_csv(FOLDER/'diagnostic_curvature_summary.csv')
    assert sum(r['status']=='UNRESOLVED' for r in optional)==2
    for r in optional:
        if r['status']=='UNRESOLVED':assert float(r['C_lower'])<=0<=float(r['C_upper'])


def test_threshold_brackets_have_resolved_signs_and_bounded_width():
    if not (FOLDER/'threshold_brackets.csv').exists():pytest.skip('Local analysis absent')
    main=w.read_csv(FOLDER/'threshold_brackets.csv');extra=w.read_csv(FOLDER/'additional_threshold_brackets.csv')
    assert len(main)==6 and len(extra)==2
    for r in main+extra:
        assert float(r['D_lower']) < -1e-8 < 1e-8 < float(r['D_upper'])
        assert float(r['chi_upper'])/float(r['chi_lower'])-1<=.015
    r30=next(r for r in main if r['pair']=='J4E0' and float(r['beta_deg'])==30)
    r90=next(r for r in main if r['pair']=='J4E0' and float(r['beta_deg'])==90)
    assert float(r30['chi_upper'])<float(r90['chi_lower'])


def test_all_suffixes_including_decreasing_curves():
    if not (FOLDER/'summary.csv').exists():pytest.skip('Local maps absent')
    for prefix in ('','diagnostic_'):
        rows=w.read_csv(FOLDER/(prefix+'spectrum.csv'));summary=w.read_csv(FOLDER/(prefix+'summary.csv'))
        for s in summary:
            rr=[r for r in rows if (r['pair'],r['chi_u'])==(s['pair'],s['chi_u'])]
            b=np.array([float(r['beta_deg']) for r in rr]);om=np.array([float(r['omega1']) for r in rr])
            assert max(abs(om-.005*np.array([float(r['Lambda1']) for r in rr])**2)/om)<2e-14
            for d in (10,5):
                result=w.previous.ranges.high_ranges(b,om,d)
                assert result['status']==s[f'status_{d}']
                assert result['beta_high']==(float(s[f'beta_high_{d}']) if s[f'beta_high_{d}'] else None)
            assert float(s['beta_high_10'])==0 # Does not imply positive angular growth.


def test_new_prefixes_have_no_duplicates_and_pass_all_gates():
    if not (FOLDER/'analysis_complete.json').exists():pytest.skip('Local analysis absent')
    rows=[]
    for name in ('spectrum','diagnostic_spectrum','addressed_spectrum'):rows+=w.read_csv(FOLDER/(name+'.csv'))
    assert len(rows)==len({(r['pair'],float(r['chi_u']),float(r['beta_deg'])) for r in rows})==711
    for r in rows:
        d=json.loads((w.ROOT/r['point_file']).read_text())
        assert d['count']['samples'][0]['count']==0
        assert [x['counts'] for x in d['roots']]==[[0,1],[1,2],[2,3]]
        assert all(all(x['gates'].values()) for x in d['roots'])


def test_no_interpolated_extrema_and_no_competing_maxima_in_sample():
    if not (FOLDER/'extremum_evidence.csv').exists():pytest.skip('Local analysis absent')
    rows=w.read_csv(FOLDER/'spectrum.csv')+w.read_csv(FOLDER/'diagnostic_spectrum.csv')
    for r in w.read_csv(FOLDER/'extremum_evidence.csv'):
        assert (r['pair'],float(r['chi_u']),float(r['beta_max_grid'])) in {(v['pair'],float(v['chi_u']),float(v['beta_deg'])) for v in rows}
        assert int(r['interior_maxima'])<=1 and int(r['interior_minima'])==0
