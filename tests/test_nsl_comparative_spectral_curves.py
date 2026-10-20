"""NSL-VIS-1 contracts, physical postprocessing and zero-solve replay."""
import json
from pathlib import Path
import sys

import numpy as np
import pytest

from scripts.analysis.nonsymmetric_laminated_beams import comparative_spectral_curves as run
from scripts.analysis.nonsymmetric_laminated_beams import verify_coupled_joint as old
from scripts.lib import nonsymmetric_laminated_joint as joint


def test_fixed_grids_and_material_control():
    config=run.read(run.CONFIG)
    beta,kappa=run.grids(config)
    assert beta==list(range(91)) and len(kappa)==61
    assert kappa[0]==0 and kappa[-1]==100 and 1 in kappa
    assert config["production_order"]==44 and config["ritz_orders"]==[36,44]
    assert config["K_guard"]==7 and config["c_theta"]==0
    c=old.read(old.CONFIG);props,_=old.materials(c)
    assert run.material_contract(props)["status"]=="PASS"
    assert props["BASE"].B==0


def test_vectorized_physical_energy_matches_existing_rayleigh_evaluator():
    c=old.read(old.CONFIG);props,_=old.materials(c)
    pair=[props["HHLL"],props["LLHH"]]
    law=joint.Joint("SPRING",c["reference"]["D"])
    # Arbitrary real fields also exercise nonzero mixed inertia after axis shift.
    pair=[joint.shift_properties(a,.003) for a in pair]
    coeff=np.random.default_rng(15).normal(size=(6*8,7))
    d=run.physical(coeff,pair,[1,1],law,64)
    previous=old.energy_ritz_frequencies(pair,[1,1],law,coeff,64)
    assert np.allclose(d["omega"],previous,rtol=2e-14,atol=0)
    k,m,_,_,_=joint.ritz_matrices(pair,[1,1],np.pi/6,law,8)
    assert np.allclose(np.diag(coeff.T@k@coeff),d["energy"],rtol=2e-14,atol=0)
    assert np.allclose(coeff.T@m@coeff,d["gram"],rtol=2e-13,atol=1e-14)


@pytest.mark.parametrize("case",["BASE","MIRROR","SAME_LOCAL"])
def test_reflection_and_independent_spring_rotations(case):
    c=old.read(old.CONFIG);props,_=old.materials(c)
    pair=[props[n] for n in c["cases"][case]]
    for law in (joint.Joint("SPRING",0),joint.Joint("SPRING",c["reference"]["D"]),joint.Joint("RIGID")):
        d=run.structural(pair,[1,1],np.pi/6,law,8,case)
        assert d["pass_structure"]
        assert d["constraint_rows"]==(3 if law.mode=="RIGID" else 2)
        assert d["reflection_claim"]==(case!="SAME_LOCAL")


def test_cluster_check_uses_subspace_not_arbitrary_individual_modes():
    c=old.read(old.CONFIG);props,_=old.materials(c)
    pair=[props["HHLL"],props["HHLL"]]
    law=joint.Joint("SPRING",c["reference"]["D"])
    coeff=np.random.default_rng(2).normal(size=(48,2))
    d=run.physical(coeff,pair,[1,1],law,64)
    scale=np.linalg.cholesky(d["gram"])
    coeff=coeff@np.linalg.inv(scale.T)
    a=run.physical(coeff,pair,[1,1],law,64)
    b=run.physical(coeff@np.array([[.6,-.8],[.8,.6]]),pair,[1,1],law,64)
    a["omega"]=b["omega"]=np.array([1.,1.])
    result=run.compare_spaces(a,b,pair,[1,1],run.read(run.CONFIG)["policy"])
    assert result["maximum_one_minus_MAC"]<1e-12
    assert result["groups"][0]["dimension"]==2
    assert not result["groups"][0]["individual_identity"]


def test_saved_maps_and_replay_do_not_call_solvers(monkeypatch,tmp_path,capsys):
    if not (run.OUTPUT/"current.json").exists():
        pytest.skip("Map production not yet complete")
    out=run.OUTPUT/run.read(run.OUTPUT/"current.json")["fingerprint"]
    result=run.validate(out)
    assert result["status"]=="PASS"
    assert result["point_counts"]==dict(angle=273,spring=183,rigid=3,unique=456,shared_angle_spring=3,benchmark_points=20)
    assert result["solve_ledger"]==dict(ritz_eigenproblems=912,boundary_root_solves=0,complex_solves=0)
    def forbidden(*args,**kwargs):
        pytest.fail("Plot/replay must not call an eigenvalue solver")
    monkeypatch.setattr(joint,"ritz_spectrum",forbidden)
    monkeypatch.setattr(joint.CoupledBeams,"prefix",forbidden)
    monkeypatch.setattr(joint,"ritz_matrices",forbidden)
    monkeypatch.setattr(sys,"argv",["curves","--compute"])
    run.main()
    assert '"new_eigenvalue_solves": 0' in capsys.readouterr().out
    rendered=run.render(out,tmp_path)
    assert rendered["new_eigenvalue_solves"]==0
    assert len(list(tmp_path.glob("*.png")))==len(list(tmp_path.glob("*.pdf")))==2
