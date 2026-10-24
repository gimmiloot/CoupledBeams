"""Bounded action/orbit checks. Tests never solve an orbit or eigenproblem."""
from pathlib import Path
import shutil

import numpy as np
import pytest

from scripts.lib import nlsp_six_vs_seven_planar_orbit as n


@pytest.fixture(scope='module')
def problem():
    cfg=n.io.read(n.CONFIG)
    if not (n.ROOT/cfg['linear_parent']/'manifest.json').exists():
        pytest.skip('Requires the optional local p14 linear bundle')
    n.verified_sources(cfg)
    return n.build_problem(cfg),cfg


@pytest.fixture(scope='module')
def saved():
    folders=list(n.OUTPUT.glob('*/manifest.json'))
    if not folders:
        pytest.skip('Requires the completed local six-field orbit bundle')
    b=folders[-1].parent
    return b,n.io.cache(b)


def test_exact_action_planar_formula_and_nu_zero(problem):
    p,cfg=problem; a=p['algebra']
    assert a['exact_action_identity'] and a['exact_planar_identity']
    assert a['nu0_restricted_action_identity'] and a['residual_max_degree']==3
    assert all(a['correction_terms_by_degree'][str(j)]>0 for j in (2,3,4))
    # No c jet occurs anywhere in the explicit physical six-field action/residual.
    for poly in (p['model6'].T4,p['model6'].V4,*p['model6'].residual_a):
        assert all(not (i<42 and i%7==6) for key in poly.terms for i in key)


def test_only_six_physical_fields_and_saved_linear_system(problem):
    p,cfg=problem; d=p['six']
    assert (d.ndof,len(d.plane),len(d.out))==(162,81,81)
    assert not set(d.plane).intersection(d.out)
    assert sorted(np.r_[d.plane,d.out])==list(range(162))
    for a,P in zip(d.arms,d.maps):
        assert a.nq==31 and 2*a.nq-1>=4*14
        np.testing.assert_array_equal(P[a.slices['c']],0.)
    assert n.linear.relative(d.K,d.saved['K'])<2e-10
    assert n.linear.relative(d.M0,d.saved['M0'])<2e-10
    assert d.joint_diagnostics(np.zeros(162))['contraction_difference'] is None


def test_independent_action_derivatives_spatial_inertia_and_planar_invariance(problem):
    p,cfg=problem
    check=n.operator_checks(p['six'],p['model6'],p['old_action'],cfg)
    plane,*space=check['directional_checks']
    assert plane['planar_mass_constant']==plane['planar_inertia']==plane['OUT_residual']==plane['cross_tangents']==0
    for row in space:
        assert row['spatial_inertia_norm']>0 and row['spatial_mass_change']>0
        assert row['residual_direction']<2e-7
        assert row['tangent_direction_second_variation']<2e-7
        assert row['energy_power_identity']<2e-9


def test_seed_transfer_discards_c_only_in_physical_field_space(problem):
    p,_=problem; d7=p['seven']; raw=np.zeros(d7.ndof)
    for i,label in enumerate(d7.saved['labels']):
        if ':c:' in label or label=='joint:c':raw[i]=.003
    contraction_only=np.linalg.solve(d7.W,raw)
    np.testing.assert_allclose(p['transfer']@contraction_only,0,atol=1e-14)
    assert p['transfer_profile_relative']<2e-10
    assert p['Ucommon'][1,1]>0 and p['Ucommon'][3,3]>0


def test_full_HB_jacobian_at_arbitrary_saved_coefficients(problem,saved):
    p,cfg=problem; b,_=saved
    z=n.load_npz(b/'HB11.npz'); U=z['U']; O=float(z['Omega'])
    solver=n.SixPlanarHB(p['six'],p['Phi6'],p['omega6'],cfg,n.time.perf_counter()+60)
    r,J,_=solver.equations(U,O,cfg['A'],True)
    rng=np.random.default_rng(19)
    direction=np.r_[rng.normal(size=U.size)*1e-7,1e-4]
    h=1e-3; du=direction[:-1].reshape(U.shape); do=direction[-1]
    fd=(solver.equations(U+h*du,O+h*do,cfg['A'])-solver.equations(U-h*du,O-h*do,cfg['A']))/(2*h)
    assert n.linear.relative(fd,J@direction)<2e-7


def test_saved_orbit_gates_and_temporal_resolution(saved):
    b,s=saved
    assert s['operator_status']==s['orbit_status']=='PASS'
    assert s['spatial_stability']=='NOT_EVALUATED'
    assert s['calls']==dict(linear_eigensolves=0,HB_solves=2,HB13_retries=0,ODE=0,Hill=0,Floquet=0,FEM=0)
    seed=n.io.read(b/'transferred_seed.json')
    assert seed['full_residual_relative']>.1
    for H in (11,13):
        row=n.io.read(b/f'HB{H}.json')
        assert row['acceptance']['passed']
        assert row['Newton']['all_planar_modes']==81
        assert row['Newton']['all_Fourier_unknowns']==81*(H+1)+1
        assert row['Newton']['scaled_projected_residual']<=1e-12
        assert row['phase_Psi']==0 and row['Q21']>0 and row['Q43']>0
        assert abs(row['actual_A']-.005)<=2e-12
        assert row['independent_phase_count']==128
        assert row['residual_harmonic_limit']==3*H
    assert s['comparison']['harmonic_convergence']['same_phase_branch']
    assert s['comparison']['harmonic_convergence']['frequency_absolute_difference']<2e-7


def test_residual_time_derivative_is_not_identity_defect(saved):
    b,_=saved; z=n.load_npz(b/'HB13.npz'); row=n.io.read(b/'HB13.json')
    nmode=z['residual_orders']; phase=z['phases']; coeff=z['residual_coefficients']
    expected=np.exp(1j*phase[:,None]*nmode)@(1j*nmode[:,None]*float(z['Omega'])*coeff)
    np.testing.assert_allclose(expected.real,z['residual_time_derivative'],rtol=2e-12,atol=1e-25)
    assert row['residual_time_derivative_max']>0
    assert row['time_shift_identity_relative_defect']<2e-9
    # High discarded harmonics are included, not silently projected away.
    assert np.linalg.norm(coeff[abs(nmode)>13])>0


def test_future_reconstruction_and_full_spatial_tangents(saved):
    b,_=saved; d=n.load_operator(b)
    q,v,a=n.reconstruct_orbit(b,np.array([.27,3.1]))
    assert q.shape==v.shape==a.shape==(2,162)
    for qi,vi,ai in zip(q,v,a):
        Rq,Rv,Ra=d.tangents(qi,vi,ai)
        assert Rq.shape==Rv.shape==Ra.shape==(162,162)
        np.testing.assert_allclose(Ra,d.mass_matrix(qi),atol=1e-14)
    period=float(n.load_npz(b/'HB13.npz')['Omega'])
    q0,v0,a0=n.reconstruct_orbit(b,[0.,2*np.pi/period])
    np.testing.assert_allclose(q0[0],q0[1],atol=1e-16)
    np.testing.assert_allclose(v0,0,atol=1e-16)


def test_cache_and_reports_do_not_build_or_solve(saved,monkeypatch,tmp_path):
    b,_=saved
    def forbidden(*args,**kwargs):raise AssertionError('No scientific recomputation allowed')
    monkeypatch.setattr(n,'build_problem',forbidden)
    monkeypatch.setattr(n.SixPlanarHB,'solve',forbidden)
    monkeypatch.setattr(n.linear.LinearCoupled,'eigenpairs',forbidden)
    found,_,cached=n.run(); assert found==b and cached
    copy=tmp_path/'report'; shutil.copytree(b,copy)
    n.render_saved(copy); n.io.cache(copy)


def test_corrupt_saved_result_is_rejected(saved,tmp_path):
    b,_=saved; copy=tmp_path/'corrupt'; shutil.copytree(b,copy)
    with (copy/'HB13.json').open('a',encoding='utf8') as stream:stream.write(' ')
    with pytest.raises(ValueError,match='artifact changed'):n.io.cache(copy)
