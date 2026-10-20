"""Two-target preset and three saved-data figures; tests never solve roots."""
import csv
import json
import numpy as np
import pytest
from scripts.analysis.nonsymmetric_laminated_beams import confirm_angle_complex as run


@pytest.fixture(scope="module")
def context():
    return run.load_context()


def forbid(*args,**kwargs):
    raise AssertionError("No eigenvalue solve or model recovery is allowed")


def test_preset_keeps_nsl3a_physics_tolerances_and_two_seeds(context):
    config,old,source,props,ref,man,kv,screen,states=context
    original=run.read(run.pilot.CONFIG)["policy"]
    for name,value in original.items():
        if name not in config["policy_overrides"]:
            assert config["policy"][name]==value
    assert len(states)==config["max_complex_solves"]==2
    assert config["max_elastic_solves"]==0
    assert config["policy"]["max_auxiliary_per_failed_mode"]==0
    assert states[0]["accepted_cross_layout_pair"]["MAC"]>.99
    assert all(s["count"]["counts"]==[4,5] for s in states)
    assert all(s["saved"]["coefficients"] for s in states)


def test_parameterized_existing_preflight_without_corrector(context,monkeypatch):
    config,old,_,props,ref,_,_,_,states=context
    monkeypatch.setattr(run.pilot.kv,"correct",forbid)
    checks=run.pilot.structural(config,old,props,ref,states)
    assert checks["status"]=="PASS"
    assert run.beta15_power(ref)["status"]=="PASS"
    for state in states:
        pr=run.pilot.provider(config,old,props,ref,state)
        assert pr.beta==np.deg2rad(15)
        assert pr.properties[0].B!=0 and pr.c==ref["c_theta"]


def test_saved_rotation_tables_preserve_mass_signs_and_endpoints(context,tmp_path,monkeypatch):
    config,old,source,_,_,man,kv,screen,_=context
    monkeypatch.setattr(run.pilot.model,"Provider",forbid)
    monkeypatch.setattr(run.elastic.CoupledBeams,"prefix",forbid)
    monkeypatch.setattr(run.elastic.CoupledBeams,"mode",forbid)
    data=run.prepare_saved_tables(config,old,source,screen,kv,man,tmp_path)
    assert data["figure1_rows"]==10 and data["figure2_rows"]==4*2*162
    for pos in (4,6):
        pair=[r for r in data["pairs"] if r["position"]==pos]
        a,b=pair
        assert a["display_factor"]*b["display_factor"]*a["signed_mass_overlap"]>0
        assert min(r["pair_MAC"] for r in pair)>.98
        assert abs(a["delta_plot"])<1e-8 and abs(b["delta_plot"])>7.
        assert a["actual_zeta"]==0 and b["actual_zeta"]>0
    assert data["heatmap"]["reused_without_rendering"]


def test_envelopes_use_damped_frequency_and_distinct_times(tmp_path):
    states=[dict(id="C1",case="MIRROR",status="PASS",alpha=.002,omega_d=.9,zeta=.002/np.hypot(.002,.9)),
            dict(id="C2",case="SAME_LOCAL",status="PASS",alpha=.003,omega_d=.8,zeta=.003/np.hypot(.003,.8))]
    config=run.read(run.CONFIG)
    modes=run.envelope_tables(states,config,tmp_path)
    assert modes[0]["N_1e"]==.9/(2*np.pi*.002)
    with (tmp_path/"figure3_envelopes.csv").open() as f:rows=list(csv.DictReader(f))
    assert len(rows)==802
    for case in ("MIRROR","SAME_LOCAL"):
        values=[r for r in rows if r["case"]==case]
        assert float(values[0]["amplitude"])==1.
        assert float(values[-1]["N"])==200
        assert np.all(np.diff([float(r["amplitude"]) for r in values])<0)
    assert rows[200]["time"]!=rows[601]["time"]
    # Ratios of alpha and zeta must not be conflated.
    states[0].update(alpha_pred=.002,zeta_pred=.002/.9)
    states[1].update(alpha_pred=.003,zeta_pred=.003/.8)
    ratio=run.comparison(states)
    assert ratio["alpha"]["ratio"]!=ratio["zeta"]["ratio"]


def completed():
    if not (run.OUTPUT/"current.json").exists():
        pytest.skip("Two-target calculation not yet performed; no test solve")
    return run.OUTPUT/run.read(run.OUTPUT/"current.json")["fingerprint"]


def test_two_saved_roots_pass_without_other_solves():
    out=completed();result=run.validate_cache(out)
    assert result["status"]=="PASS"
    assert result["solve_ledger"]["primary_complex_calls"]==2
    assert result["solve_ledger"]["auxiliary_complex_calls"]==0
    assert result["solve_ledger"]["elastic_root_solves"]==0
    for r in result["states"]:
        assert all(r["gates"].values())
        assert r["alpha"]>0 and r["omega_d"]>0 and r["MAC"]>.95
        assert r["alpha_energy_relative_error"]<1e-6


def test_cache_and_plot_only_never_call_solver_or_provider(monkeypatch,capsys,tmp_path):
    out=completed()
    monkeypatch.setattr(run.pilot.kv,"correct",forbid)
    monkeypatch.setattr(run.pilot.model,"Provider",forbid)
    monkeypatch.setattr(run.elastic.CoupledBeams,"prefix",forbid)
    monkeypatch.setattr(run.sys,"argv",["confirm_angle_complex.py","--compute"])
    run.main()
    assert '"new_complex_solves": 0' in capsys.readouterr().out
    # Render the same three figures, never a new physical calculation.
    manifest=run.figures.render(out,tmp_path,dpi=80)
    assert manifest["figure_count"]==3 and manifest["new_eigenvalue_solves"]==0
    assert len(list(tmp_path.glob("*.png")))==len(list(tmp_path.glob("*.pdf")))==3
