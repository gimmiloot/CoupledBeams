"""Targeted OUT adapters and saved-data checks; no spectral or orbit solves."""
import shutil

import numpy as np
import pytest

from scripts.lib import nlsp_six_vs_seven_out_stability as n


@pytest.fixture(scope='module')
def saved():
    folders=list(n.OUTPUT.glob('*/manifest.json'))
    if not folders:pytest.skip('Requires optional local completed OUT bundle')
    b=folders[-1].parent
    return b,n.io.cache(b)


@pytest.fixture(scope='module')
def problem(saved):
    cfg=n.io.read(n.CONFIG)
    return n.load_problem(cfg),cfg


def test_saved_full_OUT_basis_and_physical_transfer(problem):
    p,cfg=problem;d=p['d']
    assert (d.ndof,len(d.plane),len(d.out))==(162,81,81)
    assert max(p['checks'].values())<2e-9
    assert p['Phi'].shape==(81,81)
    np.testing.assert_allclose(p['L']@p['R'],np.eye(162),atol=2e-12)
    assert p['second_variation']['odd_OUT_terms']==dict(T=0,V=0)


def test_variational_tangent_and_Sdot_are_consistent(problem):
    p,cfg=problem;d=p['d'];q,v,a=n.orbit_states(p['orbits'][13],[.371])
    o=n.canonical_blocks(p['action'].blocks(q[0],v[0],a[0]));Q,V,M=d.tangents(q[0],v[0],a[0]);ix=np.ix_(d.out,d.out)
    assert n.linear.relative(o['K'],Q[ix])<2e-9
    assert n.linear.relative(o['C'],V[ix])<2e-9
    assert n.linear.relative(o['M'],M[ix])<2e-9
    actual=np.linalg.solve(o['S'].T,(o['Sdot']+o['S']@o['Av']).T).T
    assert n.linear.relative(actual,o['Ac'])<2e-9
    assert np.linalg.norm(o['Sdot'])>0
    assert n.linear.relative(np.linalg.solve(o['S'].T,(o['S']@o['Av']).T).T,o['Ac'])>2e-9
    assert not np.allclose(o['B'],o['C'],rtol=1e-5,atol=1e-14)


def test_independent_hill_convolution_and_alias_sign():
    rng=np.random.default_rng(84);orders=np.arange(-2,3);F=rng.normal(size=(5,4,4))+1j*rng.normal(size=(5,4,4))
    H,coeff=n.hb.build_hill(F,orders,.43,2,2);z=rng.normal(size=20)+1j*rng.normal(size=20)
    np.testing.assert_allclose(n.hb.hill_action(coeff,orders,.43,2,z),H@z,rtol=2e-14,atol=2e-14)
    eigen=np.array([-.02+1.1j]);folded=n.hb.fold_exponents(eigen,.43)
    assert folded[0].real==-.02
    states=np.zeros((7,2,1),complex);states[3,0,0]=1
    edge,alias=n.hb.alias_edge_weights(eigen,folded,states,.43)
    assert alias[0]==3 and edge[0]==1  # shifted central support lands at the edge


def test_saved_fourier_and_base_residual_are_separate(saved):
    b,s=saved
    for H in (11,13):
        d=s['operator_checks'][str(H)]
        assert d['count']==128 and d['order']==20
        assert d['operator_pass'] and d['Fourier_pass']
        assert d['independent_reconstruction']<2e-9
        assert d['minimum_mass_Cholesky_pivot']>0 and d['mass_change_from_rest']>0
        assert d['operator_checks']['forbidden_PL_OUT']==0
        assert d['operator_checks']['base_residual_norm']>0
        assert d['operator_checks']['base_acceleration_identity']<2e-9
        z=n.orbit.load_npz(b/f'HB{H}_Fourier_128_20.npz')
        assert z['orders'][0]<=-18 and z['orders'][-1]>=18
        np.testing.assert_allclose(z['F'],z['F'][::-1].conj(),atol=1e-13)


def test_four_saved_decompositions_keep_signed_spectra_and_partners(saved):
    b,s=saved
    assert s['calls']==dict(Hill=4,linear_eigensolves=0,HB=0,PL_Hill=0,ODE=0,FEM=0,adjoint_initial_projection=0)
    assert [(r['HB'],r['k'],r['N']) for r in s['tasks']]==[(13,72,7),(13,81,7),(13,81,9),(11,81,9)]
    for r in s['tasks']:
        z=n.orbit.load_npz(b/f"HB{r['HB']}_OUT{r['k']}_N{r['N']}.npz")
        for name in ('target','reciprocal'):
            mode=r[name];lam=complex(z[name+'_lambda'])
            assert mode['sigma']==lam.real
            assert abs(n.uncomplex(mode['rho'])-np.exp(lam*r['period']))<1e-13
            assert mode['local_residual']<1e-9 and mode['edge_weight']<1e-6
            assert mode['pointwise_residual']<.002
        assert r['reciprocal']['sigma']<0
        assert r['identity_established']
        assert r['reciprocal_pair_defect']<5e-6


def test_complex_subspace_is_not_forced_into_real_line():
    E=np.eye(4);y=np.array([1,1j,0,0]);ref=np.array([0,1,0,0])
    mac,ratio,dimension,aligned=n.subspace_match(y,ref,E,1e-8)
    assert dimension==2 and ratio==1 and mac==1
    np.testing.assert_allclose(aligned,ref)


def test_physical_recovery_and_observed_indicators(saved,problem):
    b,s=saved;p,cfg=problem;z=n.orbit.load_npz(b/'HB13_OUT81_N9.npz')
    cu=z['target_eigenvector'].reshape(19,162)
    y=np.linalg.solve(z['S0'],p['R']@cu.sum(axis=0))
    np.testing.assert_allclose(y,z['target_physical_complex'],rtol=1e-12,atol=1e-12)
    c=s['comparison'];r=s['tasks']
    expected=max(abs(r[2]['target']['sigma']-r[1]['target']['sigma']),abs(r[1]['target']['sigma']-r[0]['target']['sigma']),abs(r[2]['target']['sigma']-r[3]['target']['sigma']),1e-8)
    assert expected==c['eta6']
    assert c['eta_diff']==c['eta6']+c['eta7']
    assert c['negative_real_multiplier'] and c['growth_removed_antiperiodic_relative_defect']<1e-7


def test_cache_report_and_plot_never_solve(saved,monkeypatch,tmp_path):
    b,s=saved
    def forbidden(*args,**kwargs):raise AssertionError('Forbidden numerical call')
    monkeypatch.setattr(n,'load_problem',forbidden)
    monkeypatch.setattr(n,'spectral_task',forbidden)
    monkeypatch.setattr(n,'eig',forbidden)
    found,_,cached=n.run();assert found==b and cached
    copy=tmp_path/'report';shutil.copytree(b,copy)
    n.render_saved(copy);n.io.cache(copy)


def test_tampered_artifact_is_rejected(saved,tmp_path):
    b,s=saved;copy=tmp_path/'corrupt';shutil.copytree(b,copy)
    with (copy/'comparison.json').open('a') as f:f.write(' ')
    with pytest.raises(ValueError,match='artifact changed'):n.io.cache(copy)
