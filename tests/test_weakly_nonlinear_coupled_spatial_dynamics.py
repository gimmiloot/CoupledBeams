"""Action/mapping checks only; no IVP or native FEM job is run by tests."""
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import (
    CoupledSpatialGalerkin,EndpointArm,rectangular_coefficients)

ROOT=Path(__file__).resolve().parents[1]

@pytest.fixture(scope="module")
def model():
    raw=json.loads((ROOT/'results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json').read_text())["polynomials"]
    return SimpleNamespace(T4=rod.Polynomial.deserialize(raw["T4"]),V4=rod.Polynomial.deserialize(raw["V4"]),
        residual_a=tuple(rod.Polynomial.deserialize(p) for p in raw["residuals_A"]))

@pytest.fixture(scope="module")
def coupled(model):
    c,_=rectangular_coefficients(.1,.05)
    return CoupledSpatialGalerkin((c,c),4,model=model)

def test_rectangular_swap_preserves_mass_torsion_and_exchanges_bending():
    a,_=rectangular_coefficients(.1,.05);b,_=rectangular_coefficients(.05,.1)
    assert a.m==b.m
    assert np.isclose(a.Bp,b.Bb,rtol=1e-14)
    assert np.isclose(a.Bb,b.Bp,rtol=1e-14)
    assert np.isclose(a.CT,b.CT,rtol=2e-10)
    assert a.Bp<a.Bb and b.Bp>b.Bb

def test_endpoint_lift_and_unconstrained_derivatives(coupled):
    a=coupled.arms[0];B=a._raw_basis(np.array([0.,a.length]),0)
    assert np.array_equal(B[0],np.zeros(a.p))
    assert np.array_equal(B[1,:-1],np.zeros(a.p-1)) and B[1,-1]==1
    assert a._raw_basis(np.array([0.,a.length]),1)[0,-1]==1/a.length
    assert np.any(abs(a._raw_basis(np.array([0.]),1)[0,:-1])>0)
    assert coupled.ndof==14*(a.p-1)+7

def test_T7_covector_virtual_work(coupled):
    rng=np.random.default_rng(21);q=rng.normal(size=7);f=rng.normal(size=7)
    for T,B in zip(coupled.T7,coupled.bases):
        assert np.allclose(T.T@T,np.eye(7),atol=2e-15)
        assert np.linalg.det(B)>0
        assert np.isclose(f@(T.T@q),(T@f)@q,atol=2e-15)
    assert not np.allclose(coupled.T7[0],coupled.T7[1])

def test_joint_compatibility_balance_and_zero_work(coupled):
    rng=np.random.default_rng(5);q=rng.normal(size=coupled.ndof)*2e-7;v=rng.normal(size=coupled.ndof)*1e-7
    d=coupled.joint_diagnostics(q,v)
    assert d['translation_difference']<1e-15 and d['orientation_difference']<1e-15
    assert d['contraction_difference']<1e-15
    assert np.max(abs(d['generalized_nodal_balance']))<1e-12
    assert abs(d['joint_power'])<1e-15
    assert abs(coupled.energy_rate(q,v))<1e-15

def test_mass_is_full_action_congruence(coupled):
    q=np.sin(np.arange(coupled.ndof))*3e-7
    exact=sum(P.T@a.mass_matrix(qi)@P for a,P,qi in zip(coupled.arms,coupled.maps,coupled.local(q)))
    assert np.allclose(coupled.mass_matrix(q),exact,rtol=2e-13,atol=2e-13)
    assert np.linalg.eigvalsh(exact).min()>0
    assert np.linalg.norm(exact-coupled.M0)>1e-7

def test_analytic_full_jacobian_includes_shared_node(coupled):
    rng=np.random.default_rng(7);x=rng.normal(size=2*coupled.ndof)*2e-7;direction=rng.normal(size=len(x));direction/=np.linalg.norm(direction)
    eps=1e-8;central=(coupled.rhs(0,x+eps*direction)-coupled.rhs(0,x-eps*direction))/(2*eps)
    actual=coupled.jacobian(0,x)@direction
    assert np.linalg.norm(actual-central)/np.linalg.norm(central)<2e-7

def test_plane_invariance_and_jacobian_separation(coupled):
    q=np.zeros(coupled.ndof);v=q.copy();q[coupled.plane]=np.sin(coupled.plane)*1e-7;v[coupled.plane]=np.cos(coupled.plane)*2e-7;x=np.r_[q,v]
    assert np.max(abs(coupled.rhs(0,x)[coupled.state_out]))<1e-14
    J=coupled.jacobian(0,x)
    assert np.max(abs(J[np.ix_(coupled.state_out,coupled.state_plane)]))<1e-12
    assert np.max(abs(J[np.ix_(coupled.state_plane,coupled.state_out)]))<1e-12

def test_mirror_parity_and_odd_response(coupled):
    rng=np.random.default_rng(30);x=rng.normal(size=2*coupled.ndof)*1e-7;s=np.ones(len(x));s[coupled.state_out]=-1
    assert np.allclose(coupled.rhs(0,s*x),s*coupled.rhs(0,x),rtol=1e-12,atol=1e-13)
    assert np.isclose(coupled.energy(x[:coupled.ndof],x[coupled.ndof:]),coupled.energy((s*x)[:coupled.ndof],(s*x)[coupled.ndof:]),rtol=1e-13)

def test_global_dead_body_load_projects_to_both_axial_and_transverse(coupled):
    f=coupled.body_force(np.array([0.,-1.,0.]));assert np.max(abs(f[coupled.out]))<1e-14
    assert np.linalg.norm(f)>0
    components=coupled.bases[1].T@np.array([0.,-1.,0.])
    assert np.isclose(components[0],.5) and np.isclose(components[1],-np.sqrt(3)/2)

def test_reference_energy_parts_match_full_outplane_operators(coupled):
    q=np.zeros(coupled.ndof);v=q.copy();q[coupled.out]=np.sin(coupled.out)*1e-7;v[coupled.out]=np.cos(coupled.out)*1e-7
    part=coupled.reference_energy_parts(q,v).sum()
    full=q@coupled.K@q+v@coupled.M0@v
    assert np.isclose(part,full,rtol=1e-12)
