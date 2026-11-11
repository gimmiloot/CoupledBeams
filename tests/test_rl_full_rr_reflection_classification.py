"""Projection algebra and conservative classification, without spectral solves.

Synthetic arbitrary fields test algebra, not the eigen-equations. Classification
decision tests explicitly supply validated-physics flags. A saved physical mode
separately tests the unchanged production ODE/energy/trace diagnostics.
"""
import copy
import json
from pathlib import Path

import numpy as np
import pytest
from scripts.analysis import classify_rl_full_rr_reflection_modes as c


@pytest.fixture
def model():
    return c.make_model(**c.EXPECTED_MODEL)


@pytest.fixture
def arbitrary():
    rng = np.random.default_rng(71)
    x = (np.arange(12)+.5)/12
    return dict(arms=[dict(x=x.copy(), weights=np.full(12, 1/12),
        y=rng.normal(size=(12, 6)), dy=rng.normal(size=(12, 6)),
        ends=rng.normal(size=(2, 6))) for _ in range(2)])


def ideal_checks(d):
    """Decision-only fixture; does not assert arbitrary fields are eigenmodes."""
    d = copy.deepcopy(d)
    for name in ('original', 'plus', 'minus'):
        if d[name]['mass'] > 0:
            d[name].update(ode=0., joint=0., external=0., energy_error=0.)
    return d


POLICY = dict(ode_tol=1e-9, boundary_tol=1e-9, energy_tol=5e-8)


@pytest.mark.parametrize('complex_fields', [False, True])
def test_projectors_full_mass_and_reconstruction(model, arbitrary, complex_fields):
    if complex_fields:
        arbitrary = c.combine(arbitrary, c.sections.reflect_mode(arbitrary), 1+2j, .2-.5j)
    with c.no_spectral_calls() as calls:
        plus, minus = c.project(arbitrary, 1), c.project(arbitrary, -1)
        reflection = c.sections.reflect_mode(arbitrary)
        assert c.relative_fields(c.sections.reflect_mode(reflection), arbitrary, arbitrary) == 0
        assert c.relative_fields(c.project(plus, 1), plus, arbitrary) < 1e-15
        assert c.relative_fields(c.project(minus, -1), minus, arbitrary) < 1e-15
        assert c.relative_fields(c.combine(plus, minus), arbitrary, arbitrary) < 1e-15
        assert c.mass(model, reflection) == pytest.approx(c.mass(model, arbitrary), rel=2e-15)
        assert abs(c.mass_inner(model, plus, minus)) < c.mass(model, arbitrary)*1e-15
        assert c.mass(model, plus)+c.mass(model, minus) == pytest.approx(c.mass(model, arbitrary), rel=2e-15)
    assert not calls


@pytest.mark.parametrize('field,column,coefficient', [('y',2,'rhoI'),('dy',0,'J')])
def test_metric_retains_rotation_and_longitudinal_gradient(model, arbitrary, field, column, coefficient):
    for arm in arbitrary['arms']:
        for key in c.FIELDS: arm[key] *= 0
        arm[field][:, column] = 1
    expected = 2*(model.arms[0].J if coefficient == 'J' else model.arms[0].section.rhoI)
    assert c.mass(model, arbitrary) == pytest.approx(expected, rel=1e-15)
    assert expected > 0


@pytest.mark.parametrize('sign,status', [(1,'PURE_PLUS'),(-1,'PURE_MINUS')])
def test_synthetic_class_and_global_sign_phase(model, arbitrary, sign, status):
    case = c.PHYSICS.Case(30.,0.,'RIGID','RIGID',0.)
    mode = c.project(arbitrary, sign)
    for factor in (1, -1, 1j):
        v = c.combine(mode, mode, factor, 0)
        d = ideal_checks(c.analyze_mode(model, case, {'omega': .4}, v))
        assert c.classify(d, POLICY, 1e-9)[1] == status
        assert d['identities_pass']


def test_mixed_requires_independent_multiplicity_and_physics(model, arbitrary):
    case = c.PHYSICS.Case(30.,0.,'RIGID','RIGID',0.)
    plus, minus = c.project(arbitrary,1), c.project(arbitrary,-1)
    mixed = c.combine(plus,minus,1/np.sqrt(c.mass(model,plus)),1/np.sqrt(c.mass(model,minus)))
    d = ideal_checks(c.analyze_mode(model,case,{'omega':.4},mixed))
    assert d['xi_plus'] == pytest.approx(.5)
    assert c.classify(d,POLICY,1e-9)[1] == 'UNRESOLVED_MIXED'
    assert c.classify(d,POLICY,1e-9,multiplicity=2)[1] == 'UNRESOLVED_MIXED'
    assert c.classify(d,POLICY,1e-9,multiplicity=2,multiplicity_confirmed=True)[1] == 'MIXED_DEGENERATE'
    d['minus']['joint'] = 1e-3
    assert c.classify(d,POLICY,1e-9,multiplicity=2,multiplicity_confirmed=True)[1] == 'INVALID_PROJECTION'


def test_small_complement_not_promoted_to_independent_eigenmode(model, arbitrary):
    case = c.PHYSICS.Case(30.,0.,'RIGID','RIGID',0.)
    plus,minus=c.project(arbitrary,1),c.project(arbitrary,-1)
    mixed=c.combine(plus,minus,1,1e-13)
    d=ideal_checks(c.analyze_mode(model,case,{'omega':.4},mixed))
    d['minus']['energy_error']=.5
    assert c.classify(d,POLICY,1e-9)[1]=='PURE_PLUS'
    d['plus']['ode']=1e-3
    assert c.classify(d,POLICY,1e-9)[1]=='INVALID_PROJECTION'


def test_close_distinct_frequencies_never_establish_multiplicity():
    roots=[dict(omega=.4,multiplicity=2,event_id=7),dict(omega=.4+1e-12,multiplicity=2,event_id=7)]
    record=dict(roots=roots,search=dict(events=[dict(event_id=7,nullity=2,multiplicity=2)]))
    assert not c.confirmed_multiple_event(record,roots[0])
    record['roots'][1]['omega']=.4
    assert c.confirmed_multiple_event(record,roots[0])


def test_threshold_requires_observed_separation():
    selected=c.threshold_interval([dict(delta_plus=1e-12,delta_minus=2),dict(delta_plus=2,delta_minus=1e-11)])
    assert selected['lower'] > 1e-11 and selected['upper'] < 2
    assert c.threshold_interval([dict(delta_plus=1,delta_minus=1.7)])['chosen'] is None


def row(beta,rank,value,sign,status=None):
    return dict(beta_deg=beta,sorted_position=rank,Lambda=value,symmetry_class=sign,
                classification_status=status or ('PURE_PLUS' if sign=='+' else 'PURE_MINUS'),
                source_shape_key=f'mode{rank}',source_npz='synthetic.npz')


def test_class_rank_and_crossing_from_actual_nodes_only():
    rows=[row(6.,1,3.8,'-'),row(6.,2,3.9,'+'),row(6.1,1,3.9,'+'),row(6.1,2,4.,'-')]
    spec,coverage,gaps,brackets=c.spectra_and_crossings(rows)
    assert all(s['class_rank']==1 for s in spec)
    assert len(brackets)==1
    assert (brackets[0]['beta_left'],brackets[0]['beta_right'])==(6.,6.1)
    assert brackets[0]['order_exchanged']
    assert all(x['first_available'] for x in coverage)
    assert {s['Lambda'] for s in spec}=={r['Lambda'] for r in rows}


def test_unresolved_low_rank_blocks_false_class_minimum():
    rows=[row(0.,1,1.,'', 'INVALID_PROJECTION'),row(0.,2,2.,'+'),row(0.,3,3.,'-')]
    _,coverage,gaps,brackets=c.spectra_and_crossings(rows)
    assert not any(x['first_available'] for x in coverage)
    assert gaps[0]['gap'] is None and not brackets


def test_guard_forbids_basis_and_count():
    with c.no_spectral_calls() as calls:
        with pytest.raises(AssertionError,match='saved fields'):
            c.PHYSICS.base.basis(None,None,None,None)
        with pytest.raises(AssertionError,match='saved fields'):
            c.PHYSICS.count_at(None,None,None,None,None)
    assert len(calls)==2


def test_saved_physical_mode_has_valid_projection_without_solver():
    if not c.DEFAULT_INVENTORY.exists(): pytest.skip('Local saved archive unavailable')
    item=json.loads(c.DEFAULT_INVENTORY.read_text())['points'][0]
    record=json.loads((c.ROOT/item['record']).read_text())
    model=c.make_model(**record['contract']['model']);case=c.PHYSICS.Case(**record['contract']['case'])
    with np.load(c.ROOT/item['npz'],allow_pickle=False) as arrays:
        mode=dict(arms=[{k:arrays[f'mode1_arm{i}_{k}'] for k in ('x','weights',*c.FIELDS)} for i in (1,2)])
    with c.no_spectral_calls() as calls:
        d=c.analyze_mode(model,case,record['roots'][0],mode)
        assert c.classify(d,record['contract']['policy'],1e-9)[1] in ('PURE_PLUS','PURE_MINUS')
        assert d['original']['energy_error'] < record['contract']['policy']['energy_tol']
    assert not calls
