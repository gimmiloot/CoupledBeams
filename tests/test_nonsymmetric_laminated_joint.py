"""NSL-2 new physics and numerical path; no historical spectrum recomputation."""
from dataclasses import replace
import numpy as np
import pytest
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as joint
from scripts.analysis.nonsymmetric_laminated_beams import verify_coupled_joint as run


@pytest.fixture(scope="module")
def setup():
    config=run.read(run.CONFIG)
    props,_=run.materials(config)
    return config,props


@pytest.mark.parametrize("beta",[0,30,90])
@pytest.mark.parametrize("rigid",[False,True])
def test_covariant_physical_node_and_rank(setup,beta,rigid):
    config,p=setup
    law=joint.Joint("RIGID") if rigid else joint.Joint("SPRING",config["reference"]["D"])
    model=joint.CoupledBeams([p["HHLL"],p["LLHH"]],[1.,1.],np.deg2rad(beta),law)
    j=joint.joint_matrix(model.beta,law)
    assert np.linalg.matrix_rank(j/np.linalg.norm(j,axis=1)[:,None])==6
    offsets=[a.B/a.A for a in model.properties]
    assert run.relative(model.boundary(.37,False),model.boundary(.37,False,offsets))<1e-11
    # Negative control: moving the reference line without moving the physical
    # center/moment back is a different node, not an invariant representation.
    assert run.relative(j,joint.covariant_joint(model.beta,law,offsets))>1e-5


def test_general_axis_shift_preserves_cross_inertia_and_work(setup):
    _,p=setup
    a=replace(p["HHLL"],j1=2e-5)
    z=.013;t=joint.axis_transform(z)
    star=joint.shift_properties(a,z)
    h=single.state_matrix(.2+.7j,a)
    assert run.relative(single.state_matrix(.2+.7j,star)@t,t@h)<1e-12
    for seed in range(3):
        y=np.random.default_rng(seed).normal(size=6);ys=t@y
        assert abs(y[:3]@y[3:]-ys[:3]@ys[3:])<1e-12
        assert abs(y[:3]@a.mass_matrix()@y[:3]-ys[:3]@star.mass_matrix()@ys[:3])<1e-15


def test_material_layout_not_mass_or_D_change(setup):
    _,p=setup
    for field in ("A","D","S","m","j2"):
        assert np.allclose([getattr(a,field) for a in p.values()],getattr(p["BASE"],field),rtol=1e-13,atol=0)
    assert p["BASE"].B==0 and np.isclose(p["HHLL"].B,-p["LLHH"].B,rtol=1e-13,atol=0)
    assert all(a.j1==0 for a in p.values())


def test_full_reflection_and_straight_orientation(setup):
    _,p=setup
    from scipy.linalg import block_diag
    h=single.state_matrix(.7j,p["HHLL"])
    hr=single.state_matrix(.7j,p["LLHH"])
    mirror=block_diag(h,hr);same=block_diag(h,h)
    s=joint.REFLECTION
    assert np.array_equal(s@s,np.eye(12))
    assert run.relative(mirror@s,s@mirror)<1e-14
    assert run.relative(same@s,s@same)>1e-5
    assert run.relative(hr,-joint.CONTINUOUS@h@joint.CONTINUOUS)<1e-14


def test_variational_spring_has_only_two_essential_node_constraints(setup):
    config,p=setup
    a=[p["HHLL"],p["LLHH"]];k=config["reference"]["D"]
    hinge=joint.ritz_matrices(a,[1,1],np.pi/6,joint.Joint("SPRING",0.),5)
    spring=joint.ritz_matrices(a,[1,1],np.pi/6,joint.Joint("SPRING",k),5)
    rigid=joint.ritz_matrices(a,[1,1],np.pi/6,joint.Joint("RIGID"),5)
    assert spring[3].shape==(2,30) and rigid[3].shape==(3,30)
    assert np.allclose(spring[0]-hinge[0],k*np.outer(spring[4],spring[4]),atol=1e-17)
    assert np.array_equal(spring[1],hinge[1])


@pytest.mark.parametrize("rigid",[False,True])
def test_new_trace_space_inertia_and_green(setup,rigid):
    config,p=setup
    law=joint.Joint("RIGID") if rigid else joint.Joint("SPRING",config["reference"]["D"])
    model=joint.CoupledBeams([p["HHLL"],p["LLHH"]],[1,1],np.pi/6,law)
    d=model.count_matrix(.47,config["policy"]["count_segments"],True)
    assert d["joint_trace_dimension"]==(3 if rigid else 4)
    assert d["count"]==d["schur_inertia_sum"]
    assert max(d["symmetry"],d["joint_DtN_symmetry"])<1e-10
    assert model.segment_bound(64)["omega_squared_lower"]>config["policy"]["omega_upper"]**2
    assert model.arms[0].green_identity(.47,64)<1e-8


def test_saved_nsl1_bundle_and_protected_sources(setup):
    config,_=setup
    _,identity=run.identity(config)
    assert identity["nsl1_verified_bundle"]=="8004648ae38d120e"


def test_single_axis_boundary_at_saved_roots_independent_exponentials(setup):
    """Changing only H is insufficient: compare both physical clamp rows."""
    from scipy.linalg import expm
    config,p=setup
    a=p["HHLL"];beam=single.SingleBeam(a,1.)
    z=a.B/a.A;t=joint.axis_transform(z);ti=np.linalg.solve(t,np.eye(6))
    scale=t@np.diag(beam.scale);star=joint.shift_properties(a,z)
    saved=run.read(run.previous.OUTPUT/config["nsl1_bundle"]/"result.json")
    for mode in saved["spectral"]["modes"]:
        w=mode["omega"]
        hs=np.linalg.solve(scale,single.state_matrix(1j*w,star).real@scale)
        actual=np.vstack([ti[:3]@scale@expm(hs*(x-.5)) for x in (0.,1.)])
        original=beam.basis(w,[0,1])[:,:3].reshape(6,6)
        assert run.relative(actual,original)<1e-11
        coeff=np.asarray(mode["coefficients"])
        y=beam.basis(w,[.13,.48,.87])@coeff
        ys=y@t.T
        assert np.allclose(np.einsum("ni,ij,nj->n",y[:,:3],a.mass_matrix(),y[:,:3]),
                           np.einsum("ni,ij,nj->n",ys[:,:3],star.mass_matrix(),ys[:,:3]),rtol=1e-13,atol=1e-15)
