"""C0 algorithm/geometry gates, before any pilot. Historical tests are unchanged."""
import json
from types import SimpleNamespace
import numpy as np
import pytest
from scipy.linalg import expm
from scripts.lib import rl_full_rr_sections_prefix as e
from scripts.analysis import verify_rl_full_rr_geometry as run


@pytest.fixture
def policy():
    return run.read_config()["policy"]


def test_individual_coefficients_reference_normalization_and_cache_identity(policy):
    keys=[]
    for mu,tau in [(0.,0.),(.01,0.),(.1,-.1),(-.1,.1),(0.,.1)]:
        m=e.p.Model(mu=mu,tau=tau);c=e.p.Case(30.,mu,"RIGID","RIGID",tau)
        contract=e.point_contract(m,c,8,9,policy);keys.append(e.digest(contract))
        h,a,i=e.p.ratios(mu,tau)
        assert m.time==e.p.base.Model().time
        for j,arm in enumerate(contract["arms"]):
            assert arm["EA"]==pytest.approx(m.section.EA*a[j])
            assert arm["EI"]==pytest.approx(m.section.EI*i[j])
            assert arm["KGA"]==pytest.approx(m.section.KGA*a[j])
            assert arm["rhoA"]==pytest.approx(m.section.rhoA*a[j])
            assert arm["rhoI"]==pytest.approx(m.section.rhoI*i[j])
            assert arm["J"]==pytest.approx(m.J*i[j])
    assert len(set(keys))==len(keys)
    with pytest.raises(ValueError,match="mismatch"):
        e.point_contract(e.p.Model(mu=.1),e.p.Case(0.,0.),8,9,policy)
    # State and numerical meaning of the frozen homogeneous matrices preserved.
    m=e.p.Model(mu=.01);c=e.p.Case(30.,.01)
    np.testing.assert_array_equal(e.p.boundary(m,c,.08),e.p.base.boundary(m.arms[0],c,.08))


def test_catalog_full_arm_contract_and_effective_domain_preflight(policy,monkeypatch):
    m=e.p.Model(mu=0,tau=.1);c=e.p.Case(30.,0.,"RIGID","RIGID",.1)
    cats={e.p.pole_key(a,L):dict(upper=1.,poles=[.99999995],arm_contract=e.arm_contract(a,L))
          for a,L in zip(m.arms,c.lengths(m))}
    assert len(cats)==2
    counter=e.PhysicalCounter(m,c,cats,policy)
    monkeypatch.setattr(e.p,"count_at",lambda *a:pytest.fail("Out-of-domain request reached physics"))
    with pytest.raises(e.PrefixFailure,match="COUNT_OUTSIDE_CATALOG"):counter(.99999996)
    list(cats.values())[0]["arm_contract"]["J"]*=2
    with pytest.raises(e.PrefixFailure,match="CATALOG_ARM_CONTRACT_MISMATCH"):
        e.PhysicalCounter(m,c,cats,policy)


@pytest.mark.parametrize("beta",[0.,30.,90.])
def test_reflection_maps_global_vectors_pseudorotation_and_energy(beta):
    m=e.p.Model(mu=.1,tau=-.1);c=e.p.Case(beta,.1,"RIGID","RIGID",-.1)
    mr=e.p.Model(mu=-.1,tau=.1);cr=e.p.Case(beta,-.1,"RIGID","RIGID",.1)
    b=np.deg2rad(beta);S=np.array([[-np.cos(b),-np.sin(b)],[-np.sin(b),np.cos(b)]])
    rng=np.random.default_rng(19)
    arms=[dict(x=np.linspace(0,L,8),weights=np.ones(8)*L/8,y=rng.normal(size=(8,6)),
               dy=rng.normal(size=(8,6)),ends=rng.normal(size=(2,6))) for L in c.lengths(m)]
    mode=dict(arms=arms);ref=e.reflect_mode(mode);frames=e.p.base.geometry.frames(beta)
    for j in range(2):
        for cols in (slice(0,2),slice(3,5)):
            before=arms[1-j]["y"][:,cols]@frames[1-j].translation.T
            after=ref["arms"][j]["y"][:,cols]@frames[j].translation.T
            np.testing.assert_allclose(after,before@S.T,atol=1e-15)
        np.testing.assert_allclose(ref["arms"][j]["y"][:,2],np.linalg.det(S)*arms[1-j]["y"][:,2],atol=1e-15)
    mc,kc,_=e.p.energies(m,c,arms);mm,kk,_=e.p.energies(mr,cr,ref["arms"])
    np.testing.assert_allclose(mc,mm[::-1],rtol=3e-15)
    np.testing.assert_allclose(kc,kk[::-1],rtol=3e-15)
    assert e.p.inner(m,mode,mode)==pytest.approx(e.p.inner(mr,ref,ref),rel=3e-15)


def test_uniform_positive_s_transfer_semigroup_independent_of_cut():
    m=e.p.base.Model();A=e.p.base.state_matrix(m,.08)
    direct=expm(2*A)
    for mu in (0.,.01,.1):
        composed=expm((1+mu)*A)@expm((1-mu)*A)
        assert np.linalg.norm(composed-direct)/np.linalg.norm(direct)<1e-12


def test_full_repeated_kernel_uses_individual_masses_and_rejects_fake_physics(policy,monkeypatch):
    m=e.p.Model(mu=.1,tau=-.1);case=e.p.Case(0.,.1,"RIGID","RIGID",-.1)
    matrix=np.diag([0.,0.,*([1.]*10)])
    monkeypatch.setattr(e.p,"boundary",lambda *a:matrix)
    monkeypatch.setattr(e.p,"balance",lambda b:(b,np.ones(12),np.ones(12)))
    modes=e.recovered_modes(m,case,dict(omega=.1,multiplicity=2),120,policy)
    gram=np.array([[e.p.inner(m,a,b) for b in modes] for a in modes])
    np.testing.assert_allclose(gram,np.eye(2),atol=3e-14)
    assert e.mode_gates(m,modes,policy)["status"]=="MODE_QUALITY_FAILED"


def test_C0_failure_prevents_C1(tmp_path,monkeypatch):
    called=[]
    monkeypatch.setattr(run,"OUTPUT",tmp_path)
    monkeypatch.setattr(run,"provenance",lambda c:dict(synthetic_test=True))
    monkeypatch.setattr(run,"reference",lambda c:(tmp_path,{},[]))
    def failed(c,stage,mu,tau,beta,folder,key,*args,**kwargs):
        called.append(stage)
        return folder/"failure.json",dict(status="SYNTHETIC_PHYSICAL_FAILURE",elapsed_seconds=0.)
    monkeypatch.setattr(run,"solve_point",failed)
    monkeypatch.setattr(run,"entry",lambda *a:dict(origin="synthetic"))
    run.run(run.read_config())
    result=json.loads(next(tmp_path.glob("*/manifest.json")).read_text())
    assert called==["C0"] and result["C0"]=="STOP_C0" and result["C1"]=="NOT_RUN"

