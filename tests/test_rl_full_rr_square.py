"""Section extension contracts; historical physical kernels are frozen."""
from dataclasses import asdict
import copy
import numpy as np
import pytest
from scripts.lib import rl_full_rr_aspect as a
from scripts.lib import rl_full_rr_sections_prefix as e
from scripts.analysis import map_rl_full_rr as m
from scripts.lib.rl_full_rr_square_checks import check_coefficients


def test_square_full_inertia_and_reference_width_cancellation():
    square=a.make_model(aspect=1.,mu=.1,tau=0.)
    prior=a.make_model(aspect=4.,mu=.1,tau=0.)
    assert square.time==prior.time
    assert square.section.rhoA/prior.section.rhoA==pytest.approx(.25)
    assert square.section.EI/prior.section.EI==pytest.approx(.25)
    assert square.J/prior.J==pytest.approx(1/34)
    assert square.J==pytest.approx(2*square.nu**2*square.section.rhoI)
    assert square.J!=pytest.approx(square.nu**2*square.section.rhoI,abs=1e-15)
    assert square.arms[0].epsilon==prior.arms[0].epsilon
    with pytest.raises(ValueError):a.Model(aspect=1,variant='RL_PLANAR_LEGACY')
    with pytest.raises(ValueError):a.Model(aspect=1,variant='ELEMENTARY_LEGACY')
    with pytest.raises(ValueError):a.frozen.Model(aspect=1)


def test_aspect_four_retains_exact_original_model_and_boundary():
    original=a.frozen.Model(mu=.1,tau=-.1)
    current=a.make_model(aspect=4.,mu=.1,tau=-.1)
    assert type(current) is a.frozen.Model and asdict(current)==asdict(original)
    case=e.p.Case(30.,.1,'RIGID','RIGID',-.1)
    for w in [.02,.3,1.2]:
        assert np.array_equal(e.p.boundary(current,case,w),e.p.boundary(original,case,w))


@pytest.mark.parametrize('mu,tau',[(0,0),(.8,0),(.1,-.8),(.1,.8),(-.1,.1)])
def test_coefficients_volume_and_slenderness(mu,tau):
    model=a.make_model(aspect=1.,mu=mu,tau=tau);case=e.p.Case(5.,mu,'RIGID','RIGID',tau)
    check=check_coefficients(model,case)
    assert check['status']=='PASS'
    for arm in check['geometry']['arms']:
        assert arm['L_over_h']==arm['L_over_b']
        assert arm['J']==pytest.approx(model.rho*model.nu**2*arm['I_P'])


def test_not_slender_geometry_is_marked_not_rejected():
    model=a.make_model(aspect=1,mu=.9,tau=0)
    g=m.geometry(model,e.p.Case(5,.9,'RIGID','RIGID',0))
    assert not g['slenderness_geometric_pass']
    assert g['arms'][0]['L_over_h']<10


def test_aspect_identity_and_catalog_contract_mismatch():
    _,policy=m.configuration(m.ROOT/'data/input/rl_full_rr_square_maps.json')
    sq=a.make_model(aspect=1,mu=.1,tau=0)
    old=a.make_model(aspect=4,mu=.1,tau=0)
    case=e.p.Case(5,.1,'RIGID','RIGID',0)
    assert m.physical_key(e.point_contract(sq,case,8,9,policy))!=m.physical_key(e.point_contract(old,case,8,9,policy))
    catalogs={e.p.pole_key(arm,L):dict(upper=2.,poles=[],arm_contract=e.arm_contract(arm,L)) for arm,L in zip(old.arms,case.lengths(old))}
    with pytest.raises(e.PrefixFailure,match='CATALOG_ARM_CONTRACT_MISMATCH'):
        e.PhysicalCounter(sq,case,catalogs,policy)


def test_control_gate_blocks_map_on_failure(monkeypatch):
    from scripts.lib import rl_full_rr_square_checks as controls
    c,policy=m.configuration(m.ROOT/'data/input/rl_full_rr_square_maps.json')
    monkeypatch.setattr(controls,'control_gate',lambda *x:dict(status='STOP_CONTROL',folder='diagnostic'))
    def no_solve(*args,**kwargs):raise AssertionError('Map solved before gate')
    monkeypatch.setattr(m,'solve',no_solve)
    with pytest.raises(RuntimeError,match='maps not started'):m.run(c,policy,c['maps'])


def test_numpy_geometry_metadata_roundtrips_without_solver(tmp_path):
    model=a.make_model(aspect=1.,mu=.1,tau=.8)
    value=m.geometry(model,e.p.Case(30.,.1,'RIGID','RIGID',.8))
    path=tmp_path/'geometry.json';m.publish(path,value)
    assert m.read(path)['slenderness_geometric_pass'] is True
