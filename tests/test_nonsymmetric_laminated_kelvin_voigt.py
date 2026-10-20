"""NSL-3A structural and saved-result regressions; no eigenvalue solves."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.linalg import block_diag
from scripts.analysis.nonsymmetric_laminated_beams import verify_kelvin_voigt_pilot as run
from scripts.lib import nonsymmetric_laminated_kelvin_voigt as model
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as elastic
from scripts.lib import inplane_kelvin_voigt as kv


@pytest.fixture(scope="module")
def saved():
    config=run.read(run.CONFIG)
    old,data,properties,ref,manifests=run.context(config)
    return config,old,properties,ref,run.seeds(config,data),manifests


def test_fixed_scope_and_normalization(saved):
    config,old,props,ref,states,_=saved
    assert [(s["case"],s["position"],s["role"]) for s in states] == [
        ("MIRROR",4,"PROTECTED"),("MIRROR",6,"PROTECTED"),("MIRROR",5,"ACTIVE"),
        ("SAME_LOCAL",4,"ACTIVE"),("SAME_LOCAL",6,"ACTIVE")]
    assert config["d_theta"]==.001 and config["beta_deg"]==30.
    assert ref["k_theta"]==old["reference"]["D"]/old["reference"]["length"]
    assert np.isclose(ref["c_theta"],.001*kv.M_REF*kv.T_REF,rtol=1e-14,atol=0)


def test_preflight_without_root_search(saved,monkeypatch):
    config,old,props,ref,states,_=saved
    def forbidden(*args,**kwargs):raise AssertionError("Unexpected root search")
    monkeypatch.setattr(kv,"correct",forbidden)
    assert run.structural(config,old,props,ref,states)["status"]=="PASS"


def test_cross_inertia_derivative_and_energy_are_preserved(saved):
    _,_,props,_,_,_=saved
    a=elastic.shift_properties(props["HHLL"],props["HHLL"].B/props["HHLL"].A)
    assert abs(a.B)<1e-18 and a.j1!=0.
    p=-.002+.73j; h=1e-6
    fd=(single.state_matrix(p+h,a)-single.state_matrix(p-h,a))/(2*h)
    assert model.relative(fd,model.state_derivative(p,a))<1e-9
    assert model.state_derivative(p,a)[3,2]==2*p*a.j1
    assert model.state_derivative(p,a)[5,0]==2*p*a.j1
    q=np.array([.2+.3j,.1-.4j,.7+.2j])
    mass=np.vdot(q,a.mass_matrix()@q).real
    no_cross=np.vdot(q,replace(a,j1=0.).mass_matrix()@q).real
    assert abs(mass-no_cross)>1e-6


def test_complex_axis_pullback_contains_original_moment(saved):
    config,old,props,ref,states,_=saved
    pr=run.provider(config,old,props,ref,states[3]);p=-.001+.69j
    offsets=[a.B/a.A for a in pr.properties]
    t=block_diag(*(elastic.axis_transform(z) for z in offsets))
    j=kv.joint_matrix(p,pr.beta,pr.k,pr.c)
    star=np.linalg.solve(t.T,j.T).T
    assert star[2,3]==offsets[0]  # M_original = M_star+z0*N_star.
    assert model.relative(pr.raw(p)[0],pr.raw(p,offsets=offsets)[0])<1e-11
    assert model.relative(j,star)>1e-5


def test_fixed_column_scales_and_complex_right_null(saved):
    config,old,props,ref,states,_=saved
    pr=run.provider(config,old,props,ref,states[4])
    rows,columns=pr.rows.copy(),pr.columns.copy()
    z=(-1e-5+1j*states[4]["omega"])*pr.t_ref
    coefficient,v=pr.null(z)
    b,_=pr.matrices(z)
    sigma=np.linalg.svd(b,compute_uv=False)[-1]
    assert np.linalg.norm(b.conj().T@b@v-sigma**2*v)<1e-12
    assert np.array_equal(coefficient,v/columns)
    pr.matrices(z+.001,derivative=True)
    assert np.array_equal(rows,pr.rows) and np.array_equal(columns,pr.columns)


def test_full_cp_and_passivity(saved):
    config,old,props,ref,states,_=saved
    pr=run.provider(config,old,props,ref,states[3]);p=-.02+.69j
    j=kv.joint_matrix(p,pr.beta,pr.k,pr.c)
    fake=kv.joint_matrix(1j*p.imag,pr.beta,pr.k,pr.c)
    assert abs(j[2,2]-fake[2,2]-pr.c*p.real)<1e-20
    assert j[2,2]!=fake[2,2]
    with pytest.raises(ValueError,match="Passive"):
        model.Provider(pr.properties,pr.lengths,pr.beta,pr.k,-pr.c,pr.t_ref,states[3]["omega"])


@pytest.mark.parametrize("position",[0,1])
def test_protected_sector_removes_both_k_and_c(saved,position):
    config,old,props,ref,states,_=saved
    pr=run.provider(config,old,props,ref,states[position]);p=-.001+.7j
    j=kv.joint_matrix(p,pr.beta,pr.k,pr.c)
    zero=kv.joint_matrix(p,pr.beta,0.,0.)
    lift=np.vstack([np.eye(6),-single.REVERSE])
    assert np.max(abs((j-zero)@lift))==0.
    full=block_diag(*(single.state_matrix(p,a) for a in pr.properties))
    assert model.relative(full@elastic.REFLECTION,elastic.REFLECTION@full)<1e-12


def test_saved_pilot_has_exactly_three_new_roots_and_preserves_sources(saved):
    config,_,_,_,_,manifests=saved
    fp,_=run.identity(config,manifests)
    path=run.OUTPUT/fp/"result.json"
    if not path.exists():pytest.skip("Pilot not yet computed; no automatic root solve")
    result=run.read(path)
    assert result["status"]=="PASS" and result["protected_sources_unchanged"]
    ledger=result["solve_ledger"]
    assert ledger["primary_complex_calls"]==3
    assert ledger["auxiliary_complex_calls"]==ledger["protected_newton_calls"]==ledger["elastic_root_solves"]==0
    for row in result["states"]:
        assert all(row["gates"].values())
        if row["role"]=="PROTECTED":
            assert row["alpha"]==0 and row["predictor_relative_difference"] is None
            assert row["correction"]["steps"]==0
        else:
            assert row["alpha"]>0 and row["omega_d"]>0


def test_validated_cache_never_calls_complex_provider_or_corrector(saved,monkeypatch,capsys):
    config,_,_,_,_,manifests=saved
    fp,_=run.identity(config,manifests)
    path=run.OUTPUT/fp/"manifest.json"
    if not path.exists():pytest.skip("Pilot not yet computed")
    before=run.sha(path)
    def forbidden(*args,**kwargs):raise AssertionError("Cache triggered numerical work")
    monkeypatch.setattr(model,"Provider",forbidden)
    monkeypatch.setattr(kv,"correct",forbidden)
    monkeypatch.setattr(run.sys,"argv",["verify_kelvin_voigt_pilot.py","--compute"])
    run.main()
    text=capsys.readouterr().out
    assert '"cache": "VALIDATED"' in text and '"new_complex_solves": 0' in text
    assert run.sha(path)==before
