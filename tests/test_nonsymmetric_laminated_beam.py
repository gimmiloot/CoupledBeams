"""NSL-1 structural tests; cached physical modes are tested without new solves."""
from dataclasses import replace
import json
import numpy as np
import pytest
from scripts.lib import nonsymmetric_laminated_beam as nsl
from scripts.lib import reddy_symmetric_laminated_beam as old
from scripts.analysis.nonsymmetric_laminated_beams import verify_single_beam as run


@pytest.fixture(scope="module")
def config():
    return run.read(run.CONFIG)


@pytest.fixture(scope="module")
def case(config):
    section,_=run.hhll(config)
    prop,audit=nsl.project_reddy_coupled_properties(section,width=.2,K=5/6,provenance="test")
    return prop,section,audit


def test_old_reducer_still_rejects_nonsymmetric(case):
    with pytest.raises(ValueError,match="Only symmetric"):
        old.reduce_to_beam_properties(case[1],width=.2,K=5/6)


def test_schur_zero_secondary_resultants_and_energy(case):
    prop,_,audit=case
    q=np.array([.13,-.7])
    allq=np.zeros(6)
    allq[[0,3]]=q
    allq[[1,2,4,5]]=audit["eliminated_strain_map"]@q
    f=audit["ABD"]@allq
    np.testing.assert_allclose(f[[1,2,4,5]],0.,atol=1e-17)
    assert abs((allq@f)/(q@audit["effective"]@q)-1)<1e-12
    np.testing.assert_allclose(audit["coupled_schur"],audit["compliance_route"],rtol=1e-12)
    assert prop.B<0 and prop.j1==0


def test_two_variants_remain_distinct(case):
    source,_=nsl.source_bachoo_properties(case[1],width=.2,K=5/6,provenance="test source")
    project=case[0]
    assert source.A>project.A
    assert source.D>project.D
    assert source.variant!=project.variant
    assert abs(source.S/project.S-1)<1e-14


@pytest.mark.parametrize("angle",[15.,45.,-45.])
def test_offaxis_is_explicitly_rejected(case,angle):
    sec=old.integrate_laminate([replace(p,angle_deg=angle) for p in case[1].plies])
    with pytest.raises(ValueError,match="off-axis"):
        nsl.project_reddy_coupled_properties(sec,width=.2,K=5/6,provenance="negative")


@pytest.mark.parametrize("p",[0.,1.3j,-.2+1.7j])
def test_energy_derived_state_control_and_work_duality(case,p):
    prop=case[0]
    # Nonzero j1 must be represented as well, even though the main case has j1=0.
    prop=replace(prop,j1=.2*np.sqrt(prop.m*prop.j2))
    h=nsl.state_matrix(p,prop)
    np.testing.assert_allclose(h,nsl.control_state_matrix(p,prop),rtol=1e-13,atol=1e-12)
    canonical=np.block([[np.zeros((3,3)),np.eye(3)],[-np.eye(3),np.zeros((3,3))]])
    assert np.linalg.norm(h.T@canonical+canonical@h)<1e-8
    rev=replace(prop,B=-prop.B,j1=-prop.j1)
    np.testing.assert_allclose(nsl.state_matrix(p,rev),nsl.REVERSE@h@nsl.REVERSE,rtol=1e-14,atol=1e-14)


def test_density_first_moment_is_independent_of_elastic_B(case):
    prop,section,_=case
    assert prop.B!=0 and prop.j1==0
    plies=list(section.plies)
    plies[0]=replace(plies[0],material=replace(plies[0].material,rho=2.))
    sec=old.integrate_laminate(plies)
    a,_=nsl.project_reddy_coupled_properties(sec,width=.2,K=5/6,provenance="density control")
    reverse=old.integrate_laminate(tuple(reversed(plies)))
    b,_=nsl.project_reddy_coupled_properties(reverse,width=.2,K=5/6,provenance="density reverse")
    assert a.j1<0
    assert abs(a.j1/b.j1+1)<1e-14
    assert abs(a.B/b.B+1)<1e-14


def test_material_analytic_and_symmetric_limits(config):
    _,_,audit=run.coefficient_audit(config)
    assert audit["status"]=="PASS"
    assert audit["maximum_error"]<1e-11
    rows=audit["source_eq6"]
    np.testing.assert_allclose([r["B_times_n"] for r in rows],rows[0]["B_times_n"],rtol=1e-13)


def test_symbolic_source_discrepancy_is_not_hidden():
    audit=run.symbolic_audit()
    assert audit["source_eq11_abc"]=="EXACT_MATCH"
    assert audit["source_eq11d"]=="PRINTED_IDENTITY_FAIL"
    assert audit["source_no_silent_repair"]


def test_ritz_basis_exact_clamps_and_independence(case,monkeypatch):
    v,_=nsl.cc_basis(12,[0.,1.])
    np.testing.assert_array_equal(v,np.zeros_like(v))
    def forbidden(*args,**kwargs):
        pytest.fail("Ritz called the propagation model")
    monkeypatch.setattr(nsl,"state_matrix",forbidden)
    k,m=nsl.ritz_matrices(case[0],1.,12)
    nsl.positive_definite(k);nsl.positive_definite(m)
    assert np.max(abs(k-k.T))<1e-14
    assert np.max(abs(m-m.T))<1e-14


def test_exact_segment_coercivity_count_and_green(case,config):
    beam=nsl.SingleBeam(case[0],1.)
    p=config["policy"]
    assert beam.segment_bound(p["count_segments"])["omega_squared_lower"]>p["omega_upper"]**2
    assert beam.count_matrix(0.,p["count_segments"])["count"]==0
    assert beam.green_identity(.7,p["count_segments"])<1e-10


@pytest.fixture(scope="module")
def accepted(config):
    fp,identity=run.identity(config)
    directory=run.OUTPUT/fp
    if not (directory/"result.json").exists():
        pytest.skip("NSL spectrum not yet computed")
    result=run.read(directory/"result.json")
    if result.get("spectral",{}).get("status")!="PASS":
        pytest.skip("No accepted NSL numerical spectrum")
    manifest=run.read(directory/"manifest.json")
    assert manifest["identity"]==identity
    assert all(run.sha(directory/n)==h for n,h in manifest["artifact_hashes"].items())
    return result


@pytest.mark.parametrize("position",range(7))
def test_saved_mode_quality_and_reversal(accepted,position):
    mode=accepted["spectral"]["modes"][position]
    assert all(mode["gates"].values())
    assert abs(mode["diagnostics"]["mass"]-1)<1e-12
    assert mode["guard"]==(position==6)
    assert mode["diagnostics"]["reversal_energy"]<1e-8
    properties=nsl.Properties(**accepted["properties"])
    nodes,weights=np.polynomial.legendre.leggauss(160)
    x,weights=(nodes+1)/2,weights/2
    q,dq=nsl.ritz_fields(np.array(mode["ritz_coefficients"]),1.,x)
    mass=np.sum(weights*np.einsum("ni,ij,nj->n",q,properties.mass_matrix(),q))
    energy=weights@(properties.A*dq[:,0]**2+2*properties.B*dq[:,0]*dq[:,2]
        +properties.D*dq[:,2]**2+properties.S*(dq[:,1]+q[:,2])**2)
    assert abs(energy/(mode["ritz_omega"]**2*mass)-1)<1e-8
    assert abs(mass-1)<1e-12
    endpoints,_=nsl.ritz_fields(np.array(mode["ritz_coefficients"]),1.,[0.,1.])
    np.testing.assert_array_equal(endpoints,np.zeros((2,3)))


def test_saved_prefix_complete_and_no_repeated_search(accepted):
    s=accepted["spectral"]
    assert s["search"]["lower_count"]==0
    assert s["search"]["upper_count"]>=7
    assert len(s["search"]["roots"])==7
    assert all(r["counts"][1]-r["counts"][0]==1 for r in s["search"]["roots"])
    assert s["mass_gram_error"]<1e-8
    assert accepted["solve_ledger"]=={"ritz_eigenproblems":4,"propagation_prefix_searches":1,
        "propagation_roots":7,"reverse_root_searches":0,"complex_root_solves":0}


def test_final_cache_never_calls_eigensolvers(accepted,monkeypatch,capsys):
    def forbidden(*args,**kwargs):
        pytest.fail("Accepted cache repeated a numerical solve")
    monkeypatch.setattr(nsl,"ritz_spectrum",forbidden)
    monkeypatch.setattr(nsl.SingleBeam,"prefix",forbidden)
    monkeypatch.setattr(nsl.SingleBeam,"mode",forbidden)
    monkeypatch.setattr("sys.argv",["verify_single_beam.py","--compute"])
    run.main()
    assert '"new_eigenvalue_solves": 0' in capsys.readouterr().out
