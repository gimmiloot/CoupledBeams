"""Scoped p14 saved-data/algebraic checks; no new HB, Hill or ODE production."""
import ast
import inspect
import shutil
import numpy as np
import pytest
from scripts.lib import nlsp_planar_1to3_internal_resonance as n


@pytest.fixture(scope='module')
def saved():
    b=n.P14_OUTPUT/'2f77d2952671122b'
    if not (b/'manifest.json').exists():pytest.skip('Scoped p14 bundle unavailable')
    return b,n.p.read(b/'summary.json'),n.p.read(b/'provenance.json')


@pytest.fixture(scope='module')
def problem(saved):return n.p14_problem(saved[2])


def test_scope_and_frozen_physical_coefficients(saved,problem):
    b,s,item=saved;d12,d,_,_,_,_,_,_=problem;cfg=item['config']
    assert d12.p==12 and d.p==14 and d.lengths==(.9,1.1) and d.beta_deg==15
    assert cfg['A']==.005 and cfg['phase_Psi']==0 and cfg['mu']==.1
    assert cfg['HB_harmonics']==[11,13] and cfg['numerical_budget_seconds']==1800
    assert len(d.plane)==108 and len(d.out)==81 and d.ndof==189
    for a,c in zip(d12.arms,d.arms):assert a.coefficients.values()==c.coefficients.values()
    assert all(n.p.sha(n.p.ROOT/name)==sha for name,sha in item['initial_state']['physical_kernel_hashes'].items())


def test_exact_saved_p12_reference_outputs(saved):
    b,_,item=saved;row=n.p.read(n.p.ROOT/item['config']['parent']/'summary.json')['cases'][0]
    assert row['Omega']==.4351155867232466 and row['OUT_sigma']==.0034114100665946817
    dynamic=n.p.read(n.p.ROOT/item['config']['dynamic_parent']/'extension_5T_checks.json')
    assert dynamic['stroboscopic_nonlinear_amplification'][-1]==1.2792923254340933
    assert dynamic['sigma_tangent_dynamic']==.0034114054948789644


def test_physical_Galerkin_transfer_on_independent_points(saved,problem):
    b,_,item=saved;d12,d,_,_,_,_,_,_=problem
    with np.load(b/'linear_reference.npz') as z:T=z['transfer'].copy()
    q=np.random.default_rng(41).normal(size=d12.ndof)*1e-7;lift=T@q
    for old,new,Ps,Pt in zip(d12.arms,d.arms,d12.maps,d.maps):
        points=np.linspace(0,1,37)**1.3*old.length
        for derivative in (0,1,2):
            original=old.reconstruct(Ps@q,points,derivative);mapped=new.reconstruct(Pt@lift,points,derivative)
            assert np.linalg.norm(original-mapped)/np.linalg.norm(original)<item['config']['lifting_relative_gate']
    np.testing.assert_allclose(T.T@d.M0@T,d12.M0,atol=2e-10)
    assert abs(q@d12.K@q-lift@d.K@lift)/abs(q@d12.K@q)<item['config']['lifting_relative_gate']


def test_linear_mode_identity_mass_orthogonality_and_eigen_residuals(saved):
    b,_,item=saved;diag=n.p.read(b/'linear_matching.json');cfg=item['config']
    assert diag['gate_A']
    for name in ('PL2','PL4','OUT2'):
        row=diag['identity'][name]
        assert row['mass_MAC']>=cfg['mode_identity_MAC_gate'] and row['MAC_margin']>=cfg['mode_identity_margin_gate']
    assert max(diag['PL_mass_orthogonality'],diag['OUT_mass_orthogonality'])<cfg['operator_relative_gate']
    assert max(diag['PL_individual_eigen_residuals']+[diag['OUT2_individual_eigen_residual']])<cfg['operator_relative_gate']


def test_HB_full_p14_branch_gauge_and_independent_residual(saved,problem):
    b,_,item=saved;d12,d,_,_,_,_,_,_=problem;cfg=item['config'];check=n.p.read(b/'HB_checks.json');assert check['gate_B']
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy()
    for H in (11,13):
        row=n.p.read(b/f'HB{H}.json')
        with np.load(b/f'HB{H}.npz') as z:
            U=z['U'].copy();phases=z['phases'];res=z['massform_residual'];energy=z['energy']
            assert len(phases)==61 and phases[0]==pytest.approx(.173)
            assert not np.any(np.isclose(phases[0],np.arange(80)*2*np.pi/80))
            scale=row['Omega']**2*np.sqrt(U[1,1]**2+81*U[3,3]**2)
            assert np.linalg.norm(res,axis=1).max()/scale==pytest.approx(row['full_residual_relative'])
            assert np.ptp(energy)/abs(energy.mean())==pytest.approx(row['energy_relative_variation'])
        q,v,_=n.hb_phases(U,row['Omega'],Phi,np.array([0.]))
        assert n.p.polynomial_w_max(d,q[0])[0]==pytest.approx(.005,abs=cfg['amplitude_absolute_gate'])
        np.testing.assert_array_equal(v[0],np.zeros(d.ndof));assert U[1,1]>0 and U[3,3]>0
        assert row['Newton']['all_planar_modes']==108 and row['full_residual_relative']<=cfg['HB_full_relative_gate']
        assert row['minimum_mass_eigenvalue']>0 and row['joint_compatibility']<cfg['joint_absolute_gate']
    assert check['frequency_harmonic_difference']<cfg['frequency_control_gate']
    assert check['branch_active_relative_change']<cfg['branch_relative_change_gate']


def test_p14_residual_at_further_independent_phases(saved,problem):
    b,_,item=saved;_,d,_,_,_,_,_,_=problem
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy()
    with np.load(b/'HB13.npz') as z:U=z['U'].copy()
    row=n.p.read(b/'HB13.json');phases=np.array([.037,.521,1.913,3.287,5.173])
    q,v,a=n.hb_phases(U,row['Omega'],Phi,phases)
    R=np.array([d.mass_matrix(qi)@ai+d.inertial_terms(qi,vi)+d.potential(qi)['gradient'] for qi,vi,ai in zip(q,v,a)])
    scale=row['Omega']**2*np.sqrt(U[1,1]**2+81*U[3,3]**2)
    assert np.max(np.linalg.norm(R,axis=1))/scale<item['config']['HB_full_relative_gate']


def test_p14_boundary_joint_reversal_and_mass(saved,problem):
    b,_,item=saved;_,d,_,_,_,_,_,_=problem
    with np.load(b/'physical_observations.npz') as z:initial=z['displacement_only_initial_full_state'].copy()
    for arm,P in zip(d.arms,d.maps):np.testing.assert_allclose(arm.reconstruct(P@initial[:d.ndof],np.array([0.])),0,atol=1e-14)
    joint=d.joint_diagnostics(initial[:d.ndof],initial[d.ndof:])
    assert max(joint[key] for key in ('translation_difference','orientation_difference','contraction_difference'))<item['config']['joint_absolute_gate']
    assert np.linalg.eigvalsh(d.mass_matrix(initial[:d.ndof])).min()>0
    rhs=d.rhs(0,initial);reverse=initial.copy();reverse[d.ndof:]*=-1
    reversed_rhs=d.rhs(0,reverse)
    np.testing.assert_allclose(reversed_rhs[:d.ndof],-rhs[:d.ndof],atol=2e-12)
    np.testing.assert_allclose(reversed_rhs[d.ndof:],rhs[d.ndof:],atol=2e-12)


def test_action_full_Jacobian_Sdot_and_reflection(saved,problem):
    b,_,item=saved;_,d,model,_,_,_,_,_=problem
    with np.load(b/'linear_reference.npz') as z:Phi=z['Phi'].copy()
    with np.load(b/'HB13.npz') as z:U=z['U'].copy()
    q,v,a=n.hb_phases(U,n.p.read(b/'HB13.json')['Omega'],Phi,np.array([.623]))
    actual=d.acceleration(q[0],v[0]);A,Ac,_,S,Sdot,J=n.old.f.action_operators(n.old.f.QuadraticAction(d,model.T4,model.V4).blocks(q[0],v[0],actual))
    full=d.jacobian(0,np.r_[q[0],v[0]]);block=full[np.ix_(d.state_out,d.state_out)];gate=item['config']['operator_relative_gate']
    assert np.linalg.norm(A-block)/np.linalg.norm(block)<gate
    assert np.linalg.norm((Sdot+S@A)@np.linalg.inv(S)-Ac)/np.linalg.norm(Ac)<gate
    assert np.linalg.norm(Ac.T@J+J@Ac)/np.linalg.norm(Ac)<gate
    assert np.max(abs(full[np.ix_(d.state_plane,d.state_out)]))==0
    assert np.max(abs(full[np.ix_(d.state_out,d.state_plane)]))==0


def test_Hill_indexing_signs_and_independent_action():
    omega=.22;Omega=.435;F=np.zeros((3,2,2),complex);orders=np.arange(-1,2)
    F[1]=np.array([[0,omega],[-omega,0]])
    matrix,coeff=n.build_hill(F,orders,Omega,2,1)
    v=np.zeros((5,2),complex);v[3]=[1,1j];lam=1j*omega-1j*Omega
    np.testing.assert_allclose(matrix@v.ravel(),lam*v.ravel(),atol=1e-15)
    np.testing.assert_allclose(n.hill_action(coeff,orders,Omega,2,v.ravel()),matrix@v.ravel(),atol=1e-15)
    np.testing.assert_allclose(n.fold_exponents(np.array([lam,lam+2j*Omega]),Omega),[lam,lam],atol=1e-15)


def test_saved_Hill_eigenpair_and_reciprocal_residuals(saved):
    b,_,item=saved;cfg=item['config'];out=n.p.read(b/'OUT_checks.json');assert out['gate_C']
    with np.load(b/'OUT_Fourier.npz') as z:F=z['F'].copy();orders=z['orders'].copy()
    row=n.p.read(b/'OUT_k81_N9.json')
    with np.load(b/'OUT_k81_N9.npz') as z:
        for key,lvalue in [('dominant_eigenvector',row['dominant_raw_lambda']),('reciprocal_eigenvector',row['reciprocal_raw_lambda'])]:
            vector=z[key];lam=complex(lvalue['real'],lvalue['imag']);action=n.hill_action(F,orders,out['Omega'],9,vector)
            assert np.linalg.norm(action-lam*vector)/(np.linalg.norm(action)+abs(lam)*np.linalg.norm(vector))<cfg['Hill_residual_gate']
    rho=complex(out['multiplier']['real'],out['multiplier']['imag']);stable=complex(out['reciprocal_multiplier']['real'],out['reciprocal_multiplier']['imag'])
    assert rho.real<0 and abs(rho.imag/rho)<cfg['real_multiplier_relative_gate']
    assert abs(rho.conjugate()*stable-1)<cfg['operator_relative_gate']


def test_resolution_failure_preserved_and_triggered_control_passes_same_gate(saved):
    b,s,item=saved;cfg=item['config'];old=n.p.read(b/'OUT_initial_48_full_qualification.json');new=n.p.read(b/'OUT_checks.json')
    assert not old['gate_C'] and old['modal_change_48_full']>cfg['spectral_control_absolute_gate']
    assert new['modal_change_48_full']==old['modal_change_48_full']
    assert new['last_modal_change_72_full']<cfg['spectral_control_absolute_gate']
    assert new['temporal_change_N7_N9']<cfg['spectral_control_absolute_gate']
    assert new['resolved_positive_growth'] and new['all_OUT_modes']==81
    assert n.p.read(b/'supplementary_modal_policy.json')['spectral_control_absolute_gate_unchanged']==cfg['spectral_control_absolute_gate']


def test_unstable_mode_matching_and_sigma_per_period(saved):
    b,_,item=saved;cfg=item['config'];cmp=n.p.read(b/'comparison.json')
    assert cmp['unstable_OUT_state_MAC']>cfg['mode_identity_MAC_gate']
    for degree in ('p12','p14'):
        row=cmp[degree];assert row['sigma_T']==pytest.approx(row['sigma']*2*np.pi/row['Omega'],rel=1e-14)
        assert row['e_folding_periods']==pytest.approx(1/row['sigma_T'],rel=1e-14)
    assert cmp['sigma_relative_difference']<cfg['strong_p_agreement_relative_gate']
    controls=n.p.read(b/'OUT_mode_resolution.json');assert controls['gate']
    assert controls['k81_N7_vs_full81_N9_MAC']>=cfg['mode_identity_MAC_gate']
    assert controls['k72_N7_vs_full81_N9_MAC']>=cfg['mode_identity_MAC_gate']


def test_displacement_only_initialization_and_actual_base_identity(saved,problem):
    b,_,item=saved;_,d,_,_,_,_,_,_=problem
    with np.load(b/'physical_observations.npz') as z:base=z['displacement_only_initial_base_state'];initial=z['displacement_only_initial_full_state']
    np.testing.assert_array_equal(initial[d.state_plane],base[d.state_plane])
    np.testing.assert_array_equal(initial[d.state_out[len(d.out):]],np.zeros(len(d.out)))
    assert n.out_polynomial_extrema(d,initial[:d.ndof],1)[1]==pytest.approx(5e-5,abs=2e-15)


def test_correct_symplectic_left_projection_and_duality(saved):
    b,_,item=saved;diag=n.p.read(b/'projection_checks.json');assert diag['projection_valid']
    with np.load(b/'physical_Floquet_and_displacement.npz') as z:
        ru=z['S']@z['OUT_unit'];rs=z['S']@z['stable_OUT_unit'];left=z['canonical_left'];disp=z['S']@z['displacement_only_OUT']
        np.testing.assert_allclose(left,rs@z['J']/(rs@z['J']@ru),atol=2e-12)
        assert left@ru==pytest.approx(1,abs=2e-12) and abs(left@rs)<2e-12
        assert left@disp==pytest.approx(diag['records'][-1]['unstable_coefficient'],rel=2e-12)
    assert diag['displacement_only_excites_growing_subspace']
    assert abs(diag['records'][-1]['unstable_coefficient'])>item['config']['growth_uncertainty_factor']*diag['projection_uncertainty_indicator']
    assert diag['projection_is_not_an_energy_fraction']


def test_Floquet_phase_invariance_of_real_projector(saved):
    b,_,item=saved
    with np.load(b/'OUT_Fourier.npz') as z:R=z['R'].copy()
    with np.load(b/'OUT_k81_N9.npz') as z:cu=R@z['dominant_eigenvector'].reshape(19,162).sum(axis=0);cs=R@z['reciprocal_eigenvector'].reshape(19,162).sum(axis=0)
    with np.load(b/'physical_Floquet_and_displacement.npz') as z:S=z['S'];E=z['metric'];J=z['J'];reference=z['OUT_unit'];disp=z['displacement_only_OUT']
    projections=[]
    for phase in (0.,.71,-1.23):
        u,_,_=n.real_Floquet_direction(cu*np.exp(1j*phase),S,E,reference)
        stable,_,_=n.real_Floquet_direction(cs*np.exp(-.3j*phase),S,E)
        uc=S@u;sc=S@stable;left=sc@J/(sc@J@uc);projections.append(u*(left@S@disp))
    np.testing.assert_allclose(projections[0],projections[1],rtol=2e-10,atol=2e-12)
    np.testing.assert_allclose(projections[0],projections[2],rtol=2e-10,atol=2e-12)


def test_secondary_PL_control_keeps_unpassed_neutral_gate(saved):
    b,s,item=saved;row=n.p.read(b/'PL_checks.json');cfg=item['config']
    assert row['non_neutral_sigma']<cfg['spectral_numerical_floor']
    expected=(row['canonical']['gate_C'] and row['time_shift_neutral_relative_residual']<=cfg['operator_relative_gate']
        and max(r['independent_eigen_residual_local_relative'] for r in row['records'])<=cfg['Hill_residual_gate'])
    assert row['gate_F']==expected and s['stages']['F']==('PASS' if expected else 'PARTIAL')


def test_history_preservation_budget_and_no_ODE_FEM(saved):
    b,s,item=saved
    assert s['numerical_seconds']<item['config']['numerical_budget_seconds']
    assert s['calls']['ODE_integrations']==s['calls']['FEM']==0
    assert all(n.p.sha(n.p.ROOT/name)==sha for name,sha in item['initial_state']['historical_file_hashes'].items())
    forbidden={'solve_ivp','integrate_attempt','tangent_attempt','Radau','DOP853','run_gmsh','run_calculix'}
    for fun in (n.p14_problem,n.p14_linear,n.p14_HB,n.p14_OUT,n.p14_modal_control,n.p14_projection,n.p14_PL,n.p14_physical_diagnostics):
        tree=ast.parse(inspect.getsource(fun));calls={node.func.id if isinstance(node.func,ast.Name) else node.func.attr for node in ast.walk(tree) if isinstance(node,ast.Call) and isinstance(node.func,(ast.Name,ast.Attribute))}
        assert not calls&forbidden


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_p14_cached_replay_without_scientific_calls(saved,tmp_path,monkeypatch,mode):
    b,s,item=saved;before=n.p.sha(b/'manifest.json');copy=tmp_path/'bundle';shutil.copytree(b,copy)
    monkeypatch.setattr(n,'p14_prepare',lambda:(copy,n.p.read(copy/'provenance.json'),n.p.cache(copy)))
    def forbidden(*args,**kwargs):pytest.fail('Scientific call during completed p14 replay')
    for name in ('p14_problem','p14_linear','p14_HB','p14_OUT','p14_modal_control','p14_projection','p14_PL','p14_comparison','p14_physical_diagnostics','p14_mode_resolution_check','hill_spectrum','canonical_fourier'):
        monkeypatch.setattr(n,name,forbidden)
    n.p14_main([mode] if mode=='--compute' else [mode,str(copy)])
    assert n.p.sha(b/'manifest.json')==before and n.p.cache(copy)['calls']==s['calls']
