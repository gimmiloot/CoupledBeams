"""Fixed-geometry PL diagnostics and hard-stop/cache contract; no ODE/FEM."""
from types import SimpleNamespace
import shutil
import numpy as np
import pytest
from scripts.lib import nlsp_beta15_nonlinear_planar_resonance as n


@pytest.fixture(scope='module')
def saved():
    folders=list(n.OUTPUT.glob('*/four_amplitude_diagnostics.json'))
    if not folders:pytest.skip('Requires the optional local beta15 algebraic diagnostic bundle')
    b=folders[0].parent
    with np.load(b/'PL_coefficients.npz') as z:coeff={k:z[k].copy() for k in ('Omega0','Omega1','Omega2','phi','q20','q22','q31','q33')}
    return b,coeff,n.p.read(b/'PL_diagnostics.json'),n.p.read(b/'four_amplitude_diagnostics.json')


def test_saved_linear_identity_and_full_p12_space(saved):
    _,coeff,diag,_=saved;cfg=n.p.read(n.CONFIG)
    with np.load(n.p.ROOT/cfg['beta_parent']/'beta_15/C_modes.npz') as z:
        # Compare the same division; divide-then-multiply is not bit-exact.
        np.testing.assert_array_equal(coeff['phi'],z['planar_initial_shape']/.005)
    assert coeff['Omega0']==pytest.approx(.434820730491217,abs=1e-15)
    assert diag['linear_eigenpair_relative_residual']<2e-8
    assert len(diag['full_planar_omega'])==92
    assert max(diag['action_massform_relative_differences'])<2e-10


def test_PL_solvability_against_scalar_variational_benchmark():
    # Verifies factors/signs of the general algorithm, not a rod replacement.
    omega=1.7;alpha=.4;beta=.9
    class ScalarMassForm:
        def residual(self,n,q,v,acc):return alpha*q*q if n==2 else beta*q*q*q
    d=SimpleNamespace(plane=np.array([0]),M0=np.eye(1),K=np.array([[omega**2]]),ndof=1)
    cfg=n.p.read(n.CONFIG);coeff,diag,_=n.pl_coefficients(d,ScalarMassForm(),np.array([.005]),omega,cfg)
    expected=3*beta/(8*omega)-5*alpha**2/(12*omega**3)
    assert coeff['Omega2']==pytest.approx(expected,rel=1e-13)
    assert abs(coeff['Omega1'])<1e-15
    np.testing.assert_allclose(coeff['q20'],-alpha/(2*omega**2),atol=1e-14)
    np.testing.assert_allclose(coeff['q22'],alpha/(6*omega**2),atol=1e-14)
    assert diag['frequency_solvability']<1e-12


def test_phase_reversibility_and_zero_initial_velocities(saved):
    _,coeff,_,rows=saved
    for r in rows:
        a=r['fundamental_coefficient_a'];q,v,acc=n.orbit(coeff,a,.7);qm,vm,am=n.orbit(coeff,a,-.7)
        np.testing.assert_array_equal(q,qm);np.testing.assert_array_equal(v,-vm);np.testing.assert_array_equal(acc,am)
        assert max(abs(n.orbit(coeff,a,0.)[1]))==0


def test_physical_amplitude_differs_from_fundamental_gauge(saved):
    _,coeff,diag,rows=saved
    assert abs(diag['fundamental_gauge'])<1e-12
    assert abs(coeff['Omega1'])<1e-12
    for r in rows:
        assert r['actual_initial_max_w']==pytest.approx(r['A'],abs=1e-14)
        assert r['fundamental_coefficient_a']!=r['A']
        assert r['initial_state_correction_energy_ratio']>0


def test_formal_residual_and_energy_orders_do_not_override_internal_resonance(saved):
    _,_,diag,rows=saved
    assert max(diag['harmonic_solve_relative_residuals'])<1e-8
    for r in rows:
        assert 3.5<r['massform_residual_observed_order']<4.5
        assert 2.5<r['energy_relative_variation_order']<3.5
        assert not r['frequency_accepted']
    near=diag['denominators'][-1]
    assert near['sorted_mode']==4 and near['harmonic']==3
    assert near['relative_frequency_gap']<.01
    assert near['forcing_modal_fraction']>1e-8
    assert near['correction_energy_norm_ratio_at_Amax']>.05
    assert rows[0]['third_harmonic_energy_ratio']>.05
    assert not diag['accepted_periodic_orbit'] and diag['gate_B']=='FAIL'


def test_hard_gate_blocks_M1_M2_and_amplitude_selection(saved):
    b,_,_,rows=saved;summary=n.p.read(b/'summary.json')
    assert summary['overall']=='PLANAR_SINGLE_MODE_EXPANSION_UNRESOLVED'
    assert summary['stages']['C']=='NOT_RUN' and not summary['candidate_amplitude_selected']
    assert all(r['M1']==r['M2']==r['sigma_pair']=='NOT_AVAILABLE' for r in rows)
    assert not list(b.glob('*effective_pair*'))


def test_two_historical_M0_definitions_are_separate_and_preserved(saved):
    b,_,_,_=saved;cfg=n.p.read(n.CONFIG);rows=n.baseline_rows(cfg)
    assert rows==n.p.read(b/'M0_historical_extrapolation.json')
    first=rows[0]['definitions']
    assert first['velocity_direct']['effective_detuning']==pytest.approx(-.007910428073002419,abs=1e-15)
    assert first['canonical_direct']['effective_detuning']==pytest.approx(-.007867007108532394,abs=1e-15)
    assert rows[-1]['definitions']['velocity_direct']['sigma']==0
    assert rows[-1]['definitions']['canonical_direct']['sigma']>0


def test_sources_physical_kernels_and_historical_bundles(saved):
    b,_,_,_=saved;provenance=n.p.read(b/'provenance.json');n.p.cache(b)
    for file in n.f.s.r.KERNELS:
        assert n.p.sha(n.p.ROOT/file)==provenance['code_hashes'][file]
    for ref in provenance['source_manifests'].values():
        assert n.p.sha(n.p.ROOT/ref['path']/'manifest.json')==ref['sha256']


def test_scope_and_gates_frozen():
    cfg=n.p.read(n.CONFIG);assert n.validate_scope(cfg)==cfg
    for key,value in [('degree',16),('beta_deg',30.),('amplitudes',[.005,.006]),('new_ODE_integrations',1),('numerical_budget_seconds',601),('internal_resonance_correction_energy_ratio_trigger',.3)]:
        with pytest.raises(ValueError):n.validate_scope({**cfg,key:value})


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_real_cache_replay_without_scientific_calls(saved,monkeypatch,mode,tmp_path):
    original,_,_,_=saved
    original_manifest=n.p.sha(original/'manifest.json')
    b=tmp_path/'cached_bundle';shutil.copytree(original,b)
    # Plot replay writes PDF metadata and its manifest. Exercise it on a real
    # copied cache, preserving historical source bundles byte-for-byte.
    monkeypatch.setattr(n,'prepare',lambda *a:(b,n.p.read(b/'provenance.json'),n.p.cache(b)))
    def forbidden(*a,**k):pytest.fail('Scientific computation in cached replay')
    for owner,name in [(n.p,'sources'),(n.p,'make_disc'),(n,'compute'),(n,'finish_preflight'),(n,'pl_coefficients'),
                       (n.p,'integrate_attempt'),(n.s.l.a,'tangent_attempt'),(n.s.c,'limited_modes')]:
        monkeypatch.setattr(owner,name,forbidden)
    n.main([mode] if mode=='--compute' else [mode,str(b)])
    assert n.p.read(b/'summary.json')['calls']['ODE']==n.p.read(b/'summary.json')['calls']['FEM']==0
    assert n.p.sha(original/'manifest.json')==original_manifest
