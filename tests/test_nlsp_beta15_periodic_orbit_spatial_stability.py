"""Fixed p12 Type-II HB / Hill checks; cache replay uses temporary copies."""
import json
import shutil
import time
import numpy as np
import pytest
from scripts.lib import nlsp_planar_1to3_internal_resonance as n
from scripts.lib import nlsp_instability_feasibility as f


@pytest.fixture(scope='module')
def saved():
    b=n.STABILITY_OUTPUT/'6b93bc782316296d'
    if not (b/'manifest.json').exists():pytest.skip('Scoped scientific bundle not available')
    return b,n.p.read(b/'summary.json'),n.p.read(b/'provenance.json')


@pytest.fixture(scope='module')
def problem(saved):
    _,_,item=saved
    return n.load_problem(n.p.read(n.p.ROOT/item['config']['parent']/'provenance.json'))


def test_saved_linear_eigenpairs_and_mass_normalization(saved,problem):
    b,_,_=saved;d,_,Phi,w=problem
    with np.load(b/'linear_reference.npz') as z:
        np.testing.assert_array_equal(w,z['planar_omega']);np.testing.assert_array_equal(Phi,z['planar_Phi'])
        assert w[1]==pytest.approx(.4348207304902722,abs=1e-14)
        assert w[3]==pytest.approx(1.3063293949106531,abs=1e-14)
        np.testing.assert_allclose(Phi.T@d.M0@Phi,np.eye(92),atol=2e-10)
        np.testing.assert_allclose(z['out_Phi'].T@d.M0[np.ix_(d.out,d.out)]@z['out_Phi'],np.eye(69),atol=2e-10)


def test_four_amplitudes_full_space_branch_and_harmonic_gates(saved):
    b,s,item=saved;cfg=item['config'];assert [r['A'] for r in s['cases']]==cfg['amplitudes']
    for r in s['cases']:
        assert r['phase_Psi']==0 and r['x']>0 and r['y']>0
        assert r['Newton']['all_planar_modes']==92 and r['Newton']['all_Fourier_unknowns']==1105
        assert r['branch_active_relative_change']<cfg['branch_relative_change_gate']
        assert r['harmonic_frequency_difference']<cfg['frequency_control_gate']
        expected=(r['full_residual_relative']<=cfg['HB_full_relative_gate'] and r['energy_relative_variation']<=cfg['energy_relative_gate'])
        assert r['gate_B']==expected
        if not r['gate_B']:assert r['status']=='NOT_EVALUATED_HB_GATE' and 'OUT_sigma' not in r
    assert s['cases'][0]['gate_B']


def test_physical_amplitude_phase_reversal_boundary_joint_mass(saved,problem):
    b,s,_=saved;d,_,Phi,_=problem
    for r in s['cases']:
        with np.load(b/f'A{r["A"]:.3f}_HB11.npz') as z:U=z['U'].copy()
        q,v,a=n.hb_phases(U,r['Omega'],Phi,np.array([0.,.71,-.71,2*np.pi]))
        assert n.p.polynomial_w_max(d,q[0])[0]==pytest.approx(r['A'],abs=2e-12)
        np.testing.assert_array_equal(v[0],np.zeros(d.ndof));np.testing.assert_array_equal(q[1],q[2]);np.testing.assert_array_equal(v[1],-v[2])
        np.testing.assert_allclose(q[0],q[3],atol=1e-15)
        assert r['joint_compatibility']<2e-9 and r['minimum_mass_eigenvalue']>0
        for arm,P in zip(d.arms,d.maps):np.testing.assert_allclose(arm.reconstruct(P@q[0],np.array([0.])),0,atol=1e-14)


def test_independent_full_planar_residual_and_energy(saved):
    b,s,_=saved
    for r in s['cases']:
        with np.load(b/f'A{r["A"]:.3f}_HB11.npz') as z:
            assert len(z['phases'])==53
            assert z['phases'][0]==pytest.approx(.173)
            scale=r['Omega']**2*np.sqrt(r['x']**2+81*r['y']**2)
            assert np.max(np.linalg.norm(z['massform_residual'],axis=1))/scale==pytest.approx(r['full_residual_relative'])
            assert np.ptp(z['energy'])/abs(np.mean(z['energy']))==pytest.approx(r['energy_relative_variation'])


def test_full_HB_analytic_directional_derivative(saved,problem):
    b,s,item=saved;d,_,Phi,w=problem;r=s['cases'][0]
    with np.load(b/'A0.005_HB11.npz') as z:U=z['U'].copy()
    hb=n.FullPlanarHB(d,Phi,w,item['config'],time.perf_counter()+30)
    _,J,_=hb.equations(U,r['Omega'],.005,True)
    rng=np.random.default_rng(2);du=rng.normal(size=U.shape)*1e-5/w;do=.1;step=1e-6
    fd=(hb.equations(U+step*du,r['Omega']+step*do,.005)-hb.equations(U-step*du,r['Omega']-step*do,.005))/(2*step)
    assert np.linalg.norm(fd-J@np.r_[du.ravel(),do])/np.linalg.norm(fd)<2e-8


def test_canonical_action_Jacobian_and_Sdot(saved,problem):
    b,s,_=saved;d,model,Phi,_=problem;r=s['cases'][0]
    with np.load(b/'A0.005_HB11.npz') as z:U=z['U'].copy()
    q,v,a=n.hb_phases(U,r['Omega'],Phi,np.array([.73]));q,v=q[0],v[0]
    action=f.QuadraticAction(d,model.T4,model.V4);actual=d.acceleration(q,v)
    Av,Ac,H,S,Sdot,J=f.action_operators(action.blocks(q,v,actual))
    block=d.jacobian(0,np.r_[q,v])[np.ix_(d.state_out,d.state_out)]
    assert np.linalg.norm(Av-block)/np.linalg.norm(block)<2e-9
    inv=np.linalg.inv(S)
    assert np.linalg.norm((Sdot+S@Av)@inv-Ac)/np.linalg.norm(Ac)<2e-9
    assert np.linalg.norm(S@Av@inv-Ac)>100*np.linalg.norm((Sdot+S@Av)@inv-Ac)
    assert np.linalg.norm(Ac.T@J+J@Ac)/np.linalg.norm(Ac)<2e-9


def test_reflection_exact_second_variation_and_decoupling(saved,problem):
    b,s,_=saved;d,model,_,_=problem;_,_,checks=f.second_variation(model)
    assert checks['odd_OUT_terms']=={'T':0,'V':0}
    assert all(r['exact_residual_difference_terms']==0 for r in checks['exact_action_residual_checks'].values())
    assert s['cases'][0]['OUT_canonical']['forbidden_cross_block_max']==0
    assert s['cases'][0]['PL_canonical']['forbidden_cross_block_max']==0


def test_canonical_duality_Fourier_reality_and_Hamiltonian(saved):
    b,s,item=saved
    for block in ('OUT','PL'):
        with np.load(b/f'A0.005_{block}_Fourier.npz') as z:
            np.testing.assert_allclose(z['L']@z['R'],np.eye(len(z['L'])),atol=2e-9)
            np.testing.assert_array_equal(z['F'],z['F'][::-1].conj())
        assert s['cases'][0][block+'_canonical']['gate_C']


def test_Hill_sign_by_constant_oscillator_and_phase_invariance():
    Omega=.4;orders=np.arange(-2,3);F=np.zeros((5,2,2),complex);F[2]=[[0,.17],[-.17,0]]
    H,c=n.build_hill(F,orders,Omega,2,1);e=np.linalg.eigvals(H)
    target=np.array([sign*.17-idx*Omega for idx in range(-2,3) for sign in (-1,1)])
    np.testing.assert_allclose(np.sort(e.imag),np.sort(target),atol=1e-14)
    # Time translation multiplies A_m by exp(i*m*phase); Hill similarity.
    F[1]=[[0,.002],[-.003,0]];F[3]=F[1].conj();ph=.37
    H0,_=n.build_hill(F,orders,Omega,2,1);H1,_=n.build_hill(F*np.exp(1j*orders*ph)[:,None,None],orders,Omega,2,1)
    D=np.diag(np.repeat(np.exp(1j*np.arange(-2,3)*ph),2))
    np.testing.assert_allclose(H1,D@H0@D.conj().T,atol=2e-16)


def test_alias_handling_near_fundamental_strip_boundary():
    Omega=.4;ev=np.array([.003+1j*(Omega/2+2e-10),.003-1j*(Omega/2+2e-10)])
    folded=n.fold_exponents(ev,Omega);assert np.all(abs(folded.imag)<=Omega/2)
    states=np.zeros((11,2,2),complex);states[5]=1
    weights,aliases=n.alias_edge_weights(ev,folded,states,Omega)
    np.testing.assert_array_equal(aliases,[1,-1]);assert np.max(weights)==0
    # An original boundary vector cannot become admissible by alias shifting.
    states[:]=0;states[0]=1;weights,_=n.alias_edge_weights(ev,folded,states,Omega)
    np.testing.assert_array_equal(weights,[1.,1.])


def test_saved_Hill_eigenpair_independent_action_and_quartet(saved):
    b,s,item=saved;r=s['cases'][0];check=r['OUT_full_p12_controls'][-1]
    assert check['independent_eigen_residual_relative']<item['config']['Hill_residual_gate']
    assert check['independent_eigen_residual_local_relative']<item['config']['Hill_residual_gate']
    assert check['independent_vs_dense_action_relative']<item['config']['operator_relative_gate']
    assert check['quartet_absolute_defect']<item['config']['spectral_control_absolute_gate']
    assert check['edge_weight']<item['config']['Hill_edge_weight_gate']
    assert abs(check['dominant_lambda']['imag'])<=r['Omega']/2
    with np.load(b/'A0.005_OUT_k69_N9.npz') as z:
        vec=z['dominant_eigenvector'];raw=complex(check['dominant_raw_lambda']['real'],check['dominant_raw_lambda']['imag'])
        assert np.linalg.norm(z['dominant_independent_action']-raw*vec)<1e-9


def test_Hill_harmonic_modal_and_orbit_controls(saved):
    _,s,item=saved;r=s['cases'][0];cfg=item['config']
    assert r['OUT_temporal_change']<cfg['spectral_control_absolute_gate']
    assert r['OUT_modal_change']<cfg['spectral_control_absolute_gate']
    assert r['OUT_orbit_refinement_change']<cfg['spectral_control_absolute_gate']
    assert r['OUT_independent_Fourier_reconstruction_relative']<cfg['operator_relative_gate']
    assert r['OUT_sigma']>cfg['growth_uncertainty_factor']*r['OUT_uncertainty_indicator']
    assert r['OUT_full_p12_controls'][-1]['modal_count']==69


def test_PL_neutral_direction_and_non_neutral_classification(saved):
    _,s,item=saved;r=s['cases'][0]
    assert r['PL_full_time_shift_neutral_residual_relative']<1e-6
    assert r['PL_neutral_band']==item['config']['spectral_control_absolute_gate']
    assert r['PL_non_neutral_sigma']<item['config']['spectral_numerical_floor']
    assert min(e['abs'] for e in r['PL_full_near_zero_eigenvalues'])<1e-6


def test_zero_base_and_no_scientific_scope_expansion(saved):
    b,s,item=saved;z=n.p.read(b/'zero_base_and_reflection_checks.json')
    assert z['action_vs_linear_relative']<item['config']['operator_relative_gate']
    assert z['zero_base_Hill']['sigma']<item['config']['spectral_numerical_floor']
    assert s['calls']['nonlinear_ODE']==s['calls']['tangent_ODE']==s['calls']['FEM']==0
    assert s['numerical_seconds']<900 and n.stability_scope(item['config'])==item['config']


def test_physical_kernels_and_every_historical_parent_hash(saved):
    b,_,item=saved;source=n.p.ROOT/item['config']['parent']
    assert n.p.sha(source/'manifest.json')==item['config']['parent_manifest_sha256']
    initial=item['initial_state']
    assert all(n.p.sha(source/name)==sha for name,sha in initial['parent_files'].items())
    assert all(n.p.sha(n.p.ROOT/name)==item['code_hashes'][name] for name in n.old.f.s.r.KERNELS)


def test_no_ODE_or_FEM_call_in_authorized_scientific_functions():
    import ast
    import inspect
    forbidden={'solve_ivp','integrate_attempt','tangent_attempt','run_gmsh','run_calculix','subprocess'}
    for function in (n.stability_compute,n.stability_full_controls,n.canonical_fourier,n.hill_spectrum):
        tree=ast.parse(inspect.getsource(function))
        calls={node.func.id if isinstance(node.func,ast.Name) else node.func.attr for node in ast.walk(tree)
            if isinstance(node,ast.Call) and isinstance(node.func,(ast.Name,ast.Attribute))}
        assert not calls&forbidden


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_cached_compute_report_plot_without_scientific_calls(saved,tmp_path,monkeypatch,mode):
    b,_,_=saved;original=n.p.sha(b/'manifest.json');copy=tmp_path/'cache';shutil.copytree(b,copy)
    item=n.p.read(copy/'provenance.json');cached=n.p.cache(copy)
    monkeypatch.setattr(n,'stability_prepare',lambda:(copy,item,cached))
    def forbidden(*args,**kwargs):pytest.fail('Scientific call during cache replay')
    for owner,name in [(n,'stability_compute'),(n,'stability_full_controls'),(n,'load_problem'),(n,'canonical_fourier'),(n,'hill_spectrum'),
        (n.p,'sources'),(n.p,'make_disc'),(n.p,'integrate_attempt')]:monkeypatch.setattr(owner,name,forbidden)
    n.stability_main([mode] if mode=='--compute' else [mode,str(copy)])
    assert n.p.sha(b/'manifest.json')==original
    assert n.p.cache(copy)['calls']['nonlinear_ODE']==0
