"""Map orchestration tests: no spectral sweep or root search."""
import csv
import json
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.analysis.laminated_beams import plot_inplane_rotational_spring_eb_beta as m


def test_grid_and_explicit_joint_states():
    ticks=m.grid_tenths()
    assert len(ticks)==len(set(ticks))==201
    assert all(t in ticks for t in (0,100,300,900))
    assert ticks[:101]==list(range(101))
    assert ticks[101:141]==list(range(105,301,5))
    assert ticks[141:]==list(range(310,901,10))
    assert len(m.STATES)==5 and dict(m.STATES)["RIGID"] is None
    assert 10000 not in dict(m.STATES).values()


@pytest.mark.parametrize("label",[s for s,_ in m.STATES])
def test_constant_dimensional_spring_and_normalization(label):
    left,right=m.case(label,0),m.case(label,90)
    assert left["k_theta"]==right["k_theta"]
    if label=="RIGID":assert left["mode"]=="RIGID" and left["k_theta"] is None
    else:assert left["k_theta"]==dict(m.STATES)[label]*(.2*.05**3/12)
    assert m.FS==pytest.approx(np.sqrt((1*.2*.05)/(1*.2*.05**3/12)))
    assert right["beta_rad"]==np.pi/2


@pytest.mark.parametrize("label",[s for s,_ in m.STATES])
def test_cached_physical_assembly_matches_existing_eb(label):
    p=m.case(label,17.)
    cache=m.Transfers();provider=m.Provider(p,cache)
    actual=provider.assembly(.3)
    expected=m.eb.boundary_assembly(.3,m.ARM,m.ARM,p["beta_rad"],m.eb.Joint(p["mode"],p["k_theta"]),m.ARM)
    np.testing.assert_array_equal(actual.dimensionless,expected.dimensionless)
    np.testing.assert_array_equal(actual.endpoint_map,expected.endpoint_map)
    second=m.Provider(m.case(label,18.),cache);second(.3)
    assert cache.expm_calls==1 and cache.hits==1


def test_window_merging_preserves_spectral_multiplicity():
    events=[SimpleNamespace(diagnostics=SimpleNamespace(detected_nullity=n)) for n in (1,2,1)]
    slots=m.slots_from_events(events)
    assert slots==[events[0],events[1],events[1],events[2]]
    assert m.merge_windows([(1,3),(2,4),(7,8)])==[[1,4],[7,8]]


def test_existing_qualified_point_is_not_recomputed(monkeypatch):
    point=m.case("k1",.1)
    state={"points":{point["point_id"]:{"status":"TARGET_CONFIRMED_GUARD_QUALIFIED"}}}
    def forbidden(*args,**kwargs):raise AssertionError("recomputation")
    monkeypatch.setattr(m,"solve_point",forbidden)
    assert not m.compute_point(state,point,None)


def test_recovery_budgets_and_repeated_error_limit():
    p=m.case("k1",.1);key="k1:TARGET_UNCONFIRMED"
    state=dict(recovery_points=[],failure_recovery_counts={})
    assert m.recovery_allowed(state,p,"TARGET_UNCONFIRMED")
    state["failure_recovery_counts"][key]=3
    assert not m.recovery_allowed(state,p,"TARGET_UNCONFIRMED")
    state["failure_recovery_counts"]={};state["recovery_points"]=list(range(30))
    assert not m.recovery_allowed(state,p,"TARGET_UNCONFIRMED")
    state["recovery_points"]=[p["point_id"]]
    assert not m.recovery_allowed(state,p,"TARGET_UNCONFIRMED")


def test_gaps_and_separate_guard_quality():
    rows=[dict(state="k1",sorted_position=1,beta_deg=b,Lambda=y,status=status)
          for b,y,status in [(0,2.,"COMPLETED"),(.1,999.,"TARGET_UNCONFIRMED"),
                             (.2,2.1,"TARGET_CONFIRMED_GUARD_QUALIFIED")]]
    rows.append(dict(state="k1",sorted_position=7,beta_deg=.2,Lambda=88.,status="TARGET_CONFIRMED_GUARD_QUALIFIED"))
    x,y=m.curve_arrays(rows,"k1",1)
    assert x[:3].tolist()==[0.,.1,.2]
    assert y[0]==2 and np.isnan(y[1]) and y[2]==2.1


def test_provider_cost_limit(monkeypatch):
    p=m.Provider(m.case("k1",.1),m.Transfers());p.builds=6000
    with pytest.raises(m.CostLimit):p(.3)


def test_failed_added_angle_remains_an_explicit_gap():
    point=m.case("k1",.15,"ADDED")
    state={"points":{point["point_id"]:dict(case=point,status="TARGET_UNCONFIRMED",
        rows=[dict(sorted_position=1,Lambda=999.)],origin="NEW_COMPUTATION")}}
    rows=m.table_rows(state)
    assert len(rows)==6 and all(r["Lambda"] is None for r in rows)
    x,y=m.curve_arrays(rows,"k1",1)
    assert .15 in x and np.isnan(y[list(x).index(.15)])


def test_audit_uses_real_angular_steps_and_detects_stiffness_violation():
    state={"points":{}}
    for label in ("k0","k0.1"):
        for beta in (0.,.1,1.):
            p=m.case(label,beta)
            state["points"][p["point_id"]]=dict(case=p,status="COMPLETED",
                rows=[dict(Omega=j*np.exp(.01*beta)) for j in range(1,7)])
    result=m.audit(state)
    assert result==dict(neighbour_flags=[],stiffness_violations=[])
    state["points"]["k0.1_b1"]["rows"][5]["Omega"]*=.95
    violations=m.audit(state)["stiffness_violations"]
    assert len(violations)==1 and violations[0]["beta_deg"]==1.
    assert violations[0]["relative_increment"][5]==pytest.approx(-.05)


def test_plot_only_does_not_call_solver_or_mutate_spectra(tmp_path,monkeypatch):
    import matplotlib.figure
    def forbidden(*args,**kwargs):raise AssertionError("numerical call in plot-only")
    for obj,name in [(m,"Provider"),(m.pilot.roots,"_scan_candidates"),(m.eb,"_scaled_transfer"),
                     (np.linalg,"svd"),(np.linalg,"det")]:
        monkeypatch.setattr(obj,name,forbidden)
    monkeypatch.setattr(matplotlib.figure.Figure,"savefig",lambda *a,**kw:None)
    rows=[dict(state=label,sorted_position=j,beta_deg=0,Lambda=2+j,status="COMPLETED") for label,_ in m.STATES for j in range(1,7)]
    path=tmp_path/"spectrum_roots.csv"
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    (tmp_path/"run_manifest.json").write_text('{}')
    (tmp_path/"diagnostics.json").write_text('{"untouched":true}')
    before=path.read_bytes();checkpoint=(tmp_path/"diagnostics.json").read_bytes()
    result=m.render(tmp_path)
    assert result["matrix_calls"]==result["root_calls"]==0
    assert path.read_bytes()==before and (tmp_path/"diagnostics.json").read_bytes()==checkpoint
