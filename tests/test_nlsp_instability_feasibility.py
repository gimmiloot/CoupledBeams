"""Exact action selection, canonical structure and pair signs, without ODE/FEM."""
import numpy as np
import pytest
from scripts.lib import nlsp_instability_feasibility as f


@pytest.fixture(scope='module')
def source_model():
    config=f.p.read(f.CONFIG);parent=f.p.read(f.p.ROOT/config['beta_parent']/'provenance.json')
    return f.p.sources(parent['physical_source_config'])[0]


def test_exact_second_variation_matches_saved_linear_OUT_residual(source_model):
    T,V,result=f.second_variation(source_model)
    assert all(v['exact_residual_difference_terms']==0 for v in result['exact_action_residual_checks'].values())
    assert not any(result['odd_OUT_terms'].values())
    assert T.homogeneous(3) and V.homogeneous(3) and T.homogeneous(4) and V.homogeneous(4)
    assert all(f.out_degree(m)==2 for m in (T-V).terms)


def test_reflection_selection_is_coefficient_exact(source_model):
    for poly in (source_model.T4,source_model.V4):
        reflected=f.rod.Polynomial({m:c*(-1)**f.out_degree(m) for m,c in poly.terms.items()})
        assert reflected==poly and not f.select_out(poly,1)


def test_quadratic_action_and_existing_Jacobian_at_saved_beta15(source_model,monkeypatch):
    def forbidden(*a,**k):pytest.fail('Scientific trajectory or eigensolve forbidden')
    for owner,name in [(f.p,'integrate_attempt'),(f.s.l.a,'tangent_attempt'),(f.s.c,'limited_modes')]:monkeypatch.setattr(owner,name,forbidden)
    config=f.p.read(f.CONFIG);T,V,_=f.second_variation(source_model)
    pc={**f.p.read(f.s.c.CONFIG),'lengths':config['lengths'],'beta_deg':15.}
    d=f.p.make_disc(pc,'C',12,source_model)
    with np.load(f.p.ROOT/config['beta_parent']/'beta_15/C_modes.npz') as z:q=z['planar_initial_shape']*.7
    velocity=-q*.2;acc=d.acceleration(q,velocity);operator=f.QuadraticAction(d,T,V)
    A,Ac,H,S,Sdot,J=f.action_operators(operator.blocks(q,velocity,acc))
    expected=d.jacobian(0,np.r_[q,velocity])[np.ix_(d.state_out,d.state_out)]
    np.testing.assert_allclose(A,expected,atol=3e-9,rtol=2e-10)
    Sinv=np.linalg.inv(S)
    np.testing.assert_allclose(S@A@Sinv+Sdot@Sinv,Ac,atol=3e-9,rtol=2e-10)
    assert np.linalg.norm(Ac.T@J+J@Ac)/np.linalg.norm(Ac)<2e-10
    assert np.linalg.norm(H-H.T)/np.linalg.norm(H)<2e-10


def test_self_pair_discriminant_mean_shift_and_real_growing_solution():
    g=.0038;shift=.0012;delta=-.0054
    B=f.pair_matrix(1j*g,-1j*g,delta,1j*shift,-1j*shift)
    effective=delta-2*shift
    np.testing.assert_allclose(np.linalg.eigvals(B)**2,g*g-effective**2/4,atol=1e-18)
    tuned=f.pair_matrix(1j*g,-1j*g,2*shift,1j*shift,-1j*shift)
    physical=np.array([1+1j,1-1j])
    np.testing.assert_allclose(tuned@physical,g*physical,atol=1e-17)
    assert physical[1]==physical[0].conjugate()
    assert np.max(np.linalg.eigvals(tuned).real)==pytest.approx(g)


def test_pair_phase_gauge_preserves_product_and_eigenvalues():
    B=f.pair_matrix(.0038j,-.0038j,-.0079)
    D=np.diag(np.exp(1j*np.array([.7,-.7])))
    transformed=D.conj().T@B@D
    assert transformed[0,1]*transformed[1,0]==pytest.approx(B[0,1]*B[1,0],abs=1e-18)
    np.testing.assert_allclose(np.sort_complex(np.linalg.eigvals(transformed)**2),np.sort_complex(np.linalg.eigvals(B)**2),atol=1e-18)


def test_weighted_sum_pair_Hamiltonian_signs_and_difference_exchange():
    h=.002j;w1=.08;w2=.22;k=h.conjugate()*w2/w1
    B=f.pair_matrix(h,k,0);W=np.diag([1/w1,-1/w2])
    np.testing.assert_allclose(B.conj().T@W+W@B,0,atol=1e-16)
    assert np.max(np.linalg.eigvals(B).real)>0
    exchange=f.pair_matrix(h,-k,0)
    np.testing.assert_allclose(np.linalg.eigvals(exchange).real,0,atol=1e-17)


def test_saved_beta15_coefficients_and_safety(source_model):
    config=f.p.read(f.CONFIG);folder=f.p.ROOT/config['beta_parent']/'beta_15'
    row=next(r for r in f.p.read(folder/'interaction_rows.json') if (r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])==('self',2,2,1))
    pair=f.selected_pair(row)
    assert pair['g']==pytest.approx(.0038214554367286607,rel=1e-12)
    assert pair['effective_detuning']==pytest.approx(-.007910428073002419,rel=1e-12)
    assert pair['pseudo_Hermitian_Hamiltonian_relative_defect']<1e-12
    assert max(abs(v['real']) for v in pair['eigenvalues'])<1e-14
    assert max(v['real'] for v in pair['formal_effective_detuning_zero_eigenvalues'])==pytest.approx(pair['g'])
    with np.load(folder/'C_harmonic_coefficients.npz') as z:
        H=z['leading_at_full_amplitude'];orders=z['orders'].tolist()
        assert complex(row['coupling_forward']['real'],row['coupling_forward']['imag'])==H[orders.index(1),1,6]


def test_frozen_physical_kernels_and_feasibility_scope():
    config=f.p.read(f.CONFIG);assert f.validate_scope(config)==config
    for key,value in [('degree',16),('beta_deg',30.),('new_ODE_integrations',1),('new_eigenproblems',1),('numerical_budget_seconds',181)]:
        with pytest.raises(ValueError):f.validate_scope({**config,key:value})
    for key in ('beta_parent','modal_parent'):
        parent=f.p.ROOT/config[key]
        assert f.p.sha(parent/'manifest.json')==config[key+'_manifest_sha256']
        provenance=f.p.read(parent/'provenance.json')
        assert all(f.p.sha(f.p.ROOT/name)==provenance['code_hashes'][name] for name in f.s.r.KERNELS)


def test_canonical_saved_pairs_and_reality_without_new_evaluations():
    bundles=list(f.OUTPUT.glob('*/canonical_sum_channel_check.json'))
    if not bundles:pytest.skip('Requires optional completed local feasibility bundle')
    b=bundles[0].parent;config=f.p.read(f.CONFIG)
    for name in ('canonical_channel_check.json','canonical_sum_channel_check.json'):
        data=f.p.read(b/name)
        assert data['projected_Hamiltonian_relative_defect']<config['Hamiltonian_relative_gate']
    data=f.p.read(b/'canonical_sum_channel_check.json')
    assert data['pair_Hamiltonian_relative_defect']<config['Hamiltonian_relative_gate']
    assert data['canonical_product']['real']>0
    with np.load(b/'canonical_sum_channel_coefficients.npz') as z:
        R=z['right'];assert np.array_equal(R[:,2:],R[:,:2].conj())
        physical=R@np.array([1+.4j,.2-.3j,1-.4j,.2+.3j])
        np.testing.assert_allclose(physical.imag,0,atol=1e-14)


def test_raw_velocity_projection_is_not_assumed_canonical():
    config=f.p.read(f.CONFIG);source=f.p.ROOT/config['modal_parent'];orders=config['harmonics']
    row=next(r for r in f.p.read(source/'C_interaction_rows.json') if (r['kind'],r['i_sorted'],r['j_sorted'],r['harmonic'])==('sum',1,2,1))
    with np.load(source/'C_harmonic_coefficients.npz') as z:H=z['leading_at_full_amplitude']
    row['mean_diagonal']=[f.s.complex_row(H[orders.index(0),0,0]),f.s.complex_row(H[orders.index(0),6,6])]
    assert f.selected_pair(row)['pseudo_Hermitian_Hamiltonian_relative_defect']>config['Hamiltonian_relative_gate']


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_cached_modes_without_action_rebuild_or_scientific_calls(monkeypatch,tmp_path,mode):
    summary={'overall':'DONE','completed':True,'hard_stop':False,'stages':{},'calls':{},'numerical_seconds':0.}
    monkeypatch.setattr(f,'prepare',lambda *x:(tmp_path,{},summary));monkeypatch.setattr(f,'cache',lambda *x:summary)
    def forbidden(*a,**k):pytest.fail('New scientific call in cache')
    for owner,name in [(f.p,'sources'),(f.p,'make_disc'),(f.p,'integrate_attempt'),(f.s.l.a,'tangent_attempt'),
                       (f,'second_variation'),(f,'QuadraticAction'),(f.s.c,'limited_modes')]:monkeypatch.setattr(owner,name,forbidden)
    plots=[];monkeypatch.setattr(f,'plot_saved',lambda b:plots.append(b))
    monkeypatch.setattr(f.p,'read',lambda *x:{});monkeypatch.setattr(f.p,'finalize',lambda *x:None)
    f.main([mode] if mode=='--compute' else [mode,str(tmp_path)])
    assert len(plots)==(mode=='--plot-only')
