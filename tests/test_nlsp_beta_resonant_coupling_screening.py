"""Five-angle scope and modal pair algebra; no ODE/FEM or new eigenproblems."""
from pathlib import Path
import numpy as np
import pytest
from scripts.lib import nlsp_resonant_spatial_coupling_screening as s


def test_frozen_five_angle_scope_and_authorization():
    config=s.p.read(s.BETA_CONFIG)
    assert s.beta_validate_scope(config)==config
    for key,value in [('mu',.2),('beta_points',[15.,30.,45.,60.,75.,80.]),('degree',16),
                      ('amplitude',.006),('new_ODE_integrations',1),('numerical_budget_seconds',301)]:
        with pytest.raises(ValueError):s.beta_validate_scope({**config,key:value})


def test_mean_detuning_can_remove_bare_pair_growth():
    config=s.p.read(s.BETA_CONFIG);orders=config['harmonics'];H=np.zeros((5,2,2),complex)
    H[orders.index(1),0,1]=.01j;H[orders.index(-1),1,0]=-.01j
    H[orders.index(0)]=np.diag([.01j,-.01j])
    problem={'omega':np.array([.5]),'Omega':.99,'mu':.1}
    row=s.pair_rows(problem,H,np.zeros_like(H.real),config)[0]
    assert row['product']['real']>0
    assert row['two_coupling_over_bare_detuning']>1
    assert row['two_coupling_over_mean_shifted_detuning']<1
    assert abs(row['isolated_pair_averaging_growth_proxy'])<1e-14
    assert not row['direct_term_candidate']


def test_pair_products_and_growth_are_phase_and_modal_gauge_invariant():
    config=s.p.read(s.BETA_CONFIG);H=np.zeros((5,4,4),complex);orders=config['harmonics']
    H[orders.index(1),0,2]=.02j;H[orders.index(-1),2,0]=-.02j
    H[orders.index(2),1,0]=.001j;H[orders.index(-2),0,1]=.001j
    problem={'omega':np.array([.5,1.5]),'Omega':1.,'mu':.1}
    result=s.beta_phase_invariance(problem,H,np.zeros_like(H.real),config)
    assert result['maximum_product_or_growth_difference']<1e-15


def test_unresolved_amplitude_order_excludes_pair_candidate():
    config=s.p.read(s.BETA_CONFIG);H=np.zeros((5,2,2),complex);orders=config['harmonics']
    H[orders.index(1),0,1]=.02j;H[orders.index(-1),1,0]=-.02j
    problem={'omega':np.array([.5]),'Omega':1.,'mu':.1,'beta_deg':30.,'torsional_mass_fraction':[.1]}
    full=np.array([H,H/8,H/64])  # Cubic scaling fails the leading O(a) gate.
    rows=s.beta_enrich_rows(problem,H,full,np.zeros_like(H.real),config,{'status':'PASS'})
    row=next(r for r in rows if r['harmonic']==1)
    assert row['amplitude_order_gate']=='UNRESOLVED'
    assert row['pair_classification']=='unresolved' and not row['direct_term_candidate']


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_beta_cached_modes_make_no_scientific_calls(monkeypatch,tmp_path,mode):
    summary={'completed':True,'hard_stop':False,'overall':'DONE','stages':{},'calls':{},'numerical_seconds':1.}
    monkeypatch.setattr(s,'beta_cache',lambda *x:summary)
    monkeypatch.setattr(s,'beta_prepare',lambda *x:(tmp_path,{},summary))
    def forbidden(*args,**kwargs):pytest.fail('scientific call during beta cache')
    for owner,name in [(s.p,'sources'),(s.p,'make_disc'),(s.p,'integrate_attempt'),(s.l.a,'tangent_attempt'),
                       (s,'modal_problem'),(s,'projected_at'),(s,'harmonic_sample')]:monkeypatch.setattr(owner,name,forbidden)
    plots=[];monkeypatch.setattr(s,'beta_plot',lambda b:plots.append(b))
    monkeypatch.setattr(s.p,'read',lambda *x:{})
    monkeypatch.setattr(s.p,'finalize',lambda *x:None)
    s.beta_main([mode] if mode=='--compute' else [mode,str(tmp_path)])
    assert len(plots)==(mode=='--plot-only')


@pytest.fixture(scope='module')
def saved_bundle():
    bundles=list(s.BETA_OUTPUT.glob('*/summary.json'))
    if not bundles:pytest.skip('Optional real-data checks require completed local beta bundle')
    assert len(bundles)==1
    return bundles[0].parent


def test_saved_geometry_axes_coefficients_and_joint_compatibility(saved_bundle):
    coefficients=None
    for beta in (15,30,45,60,75):
        data=s.p.read(saved_bundle/f'beta_{beta:02d}'/'geometry_and_safety.json')
        assert data['lengths']==[.9,1.1] and data['beta_deg']==beta
        B=np.array(data['local_to_global_bases']);angle=np.deg2rad(beta)
        np.testing.assert_allclose(B[1,:,0],[-np.cos(angle),-np.sin(angle),0],atol=2e-15)
        np.testing.assert_allclose(B[1,:,1],[-np.sin(angle),np.cos(angle),0],atol=2e-15)
        np.testing.assert_allclose(B[0,:,2],[0,0,-1],atol=2e-15)
        if coefficients is None:coefficients=data['coefficients']
        assert data['coefficients']==coefficients
        for key in ('translation_difference','orientation_difference','contraction_difference'):
            assert data['joint_checks'][key]<2e-12
        assert data['actual_max_w']==pytest.approx(.005,rel=2e-14)


def test_saved_modal_mass_orthogonality_duality_and_both_frequency_signs(saved_bundle):
    for beta in (15,30,45,60,75):
        with np.load(saved_bundle/f'beta_{beta:02d}'/'C_modes.npz') as z:
            R=z['right'];L=z['left'];E=z['metric'];phi=z['phi'];n=E.shape[0]//2;k=len(z['omega'])
            np.testing.assert_allclose(phi.T@E[n:,n:]@phi,np.eye(k),atol=2e-10)
            np.testing.assert_allclose(L@R,np.eye(2*k),atol=2e-10)
            np.testing.assert_allclose(L,R.conj().T@E,atol=2e-14)
            assert np.array_equal(R[:,k:],R[:,:k].conj())


def test_saved_control_point_and_symmetry_gates(saved_bundle):
    check=s.p.read(saved_bundle/'beta_30/reference_checks.json')
    assert all(value<1e-10 for value in check.values())
    for beta in (15,30,45,60,75):
        folder=saved_bundle/f'beta_{beta:02d}'
        assert s.p.read(folder/'planar_reflection_selection.json')['cross_blocks_max']<1e-10
        diagnostic=s.p.read(folder/'C_amplitude_order_diagnostics.json')
        assert diagnostic['conjugate_branch_Fourier_difference']<1e-10
        assert diagnostic['phase_grid_absolute_difference']<1e-10


def test_unresolved_planar_identity_is_not_forced(saved_bundle):
    metadata=s.p.read(saved_bundle/'beta_15/C_mode_metadata.json')
    correspondence=s.p.read(saved_bundle/'beta_15/mode_correspondence.json')
    assert metadata['planar']['selected_sorted_index_one_based']==2
    assert correspondence['planar_MAC']<s.p.read(s.BETA_CONFIG)['correspondence_MAC_gate']
    assert correspondence['status']=='UNRESOLVED'


def test_saved_source_manifest_and_physical_kernels_unchanged():
    config=s.p.read(s.BETA_CONFIG);parent=s.p.ROOT/config['parent']
    assert s.p.sha(parent/'manifest.json')==config['parent_manifest_sha256']
    provenance=s.p.read(parent/'provenance.json')
    assert all(s.p.sha(s.p.ROOT/f)==provenance['code_hashes'][f] for f in s.r.KERNELS)
