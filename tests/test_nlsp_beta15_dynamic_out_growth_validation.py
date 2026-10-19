"""Saved-data regression of one fixed dynamic experiment; no new IVPs."""
import ast
import inspect
import shutil
import numpy as np
import pytest
from scripts.lib import nlsp_planar_1to3_internal_resonance as n
from scripts.lib import nlsp_coupled_sensitivity_analysis as a
from scripts.lib import nlsp_instability_feasibility as f


@pytest.fixture(scope='module')
def saved():
    b=n.DYNAMIC_OUTPUT/'49dc29161a5dda78'
    if not (b/'manifest.json').exists():pytest.skip('Scoped dynamic bundle unavailable')
    return b,n.p.read(b/'summary.json'),n.p.read(b/'provenance.json')


@pytest.fixture(scope='module')
def problem(saved):return n.dynamic_problem(saved[2])


def test_HB11_source_and_initial_identity(saved,problem):
    b,_,item=saved;d,model,Phi,U,Omega,lam,coeff,R=problem
    with np.load(b/'initial_states.npz') as z:
        np.testing.assert_array_equal(z['HB_U'],U);np.testing.assert_array_equal(z['planar_Phi'],Phi)
        q,v,_=n.hb_phases(U,Omega,Phi,np.array([0.]));np.testing.assert_array_equal(z['base_state'],np.r_[q[0],v[0]])
        assert np.max(abs(z['base_state'][d.state_out]))==0
    assert n.p.polynomial_w_max(d,q[0])[0]==pytest.approx(.005,abs=2e-12)
    assert Omega==pytest.approx(.4351155867232466,abs=1e-15)
    assert n.p.read(b/'Floquet_reconstruction.json')['T']==2*np.pi/Omega


def test_fixed_geometry_coefficients_and_scope(saved,problem):
    _,_,item=saved;d,_,_,_,_,_,_,_=problem;cfg=item['config']
    assert d.p==12 and d.beta_deg==15 and d.lengths==(.9,1.1)
    for arm in d.arms:assert arm.coefficients.Bp/arm.coefficients.Bb==pytest.approx(4)
    assert cfg['A']==.005 and cfg['phase_Psi']==0 and cfg['numerical_budget_seconds']==3600
    assert cfg['conditional_horizon_periods']==5 and not cfg['automatic_retry'] and cfg['new_FEM_jobs']==0


def test_phase_and_zero_PL_velocities_not_arbitrarily_reset(saved,problem):
    b,_,_=saved;d,_,Phi,U,Omega,_,_,_=problem
    q,v,_=n.hb_phases(U,Omega,Phi,np.array([0.,.71,-.71]))
    # Exact parity is checked before BLAS multiplication. Physical reconstruction
    # has a dimension-dependent floating-point summation bound, not bitwise parity.
    h=np.arange(len(U));angle=.71*h
    np.testing.assert_array_equal(np.cos(angle),np.cos(-angle))
    np.testing.assert_array_equal(np.sin(angle),-np.sin(-angle))
    operations=len(U)+Phi.shape[1];eps=np.finfo(float).eps
    gamma=operations*eps/(1-operations*eps)
    position_bound=2*gamma*(abs(U).sum(axis=0)@abs(Phi.T))
    velocity_bound=2*gamma*Omega*((h[:,None]*abs(U)).sum(axis=0)@abs(Phi.T))
    assert np.all(abs(q[1]-q[2])<=position_bound)
    assert np.all(abs(v[1]+v[2])<=velocity_bound)
    assert np.max(abs(v[0]))==0
    assert n.p.read(b/'Floquet_reconstruction.json')['all_initial_OUT_velocities_preserved']


def test_canonical_to_physical_Floquet_initial_mapping(saved,problem):
    b,_,_=saved;d,model,Phi,U,Omega,lam,coeff,R=problem
    q,v,acc=n.hb_phases(U,Omega,Phi,np.array([0.]));S=f.action_operators(f.QuadraticAction(d,model.T4,model.V4).blocks(q[0],v[0],acc[0]))[3]
    with np.load(b/'initial_states.npz') as z:
        phase=z['phase_combination'].item();unit=z['OUT_unit_state'];canonical=R@coeff.sum(axis=0)
        np.testing.assert_allclose(S@unit,np.real(phase*canonical),atol=2e-12)
        np.testing.assert_allclose(z['complex_initial'],np.linalg.solve(S,canonical),atol=2e-12)
        assert np.max(abs(unit[len(d.out):]))>0
    rec=n.p.read(b/'Floquet_reconstruction.json')
    assert rec['kinematic_reconstruction_relative_error_max']<.002
    assert rec['source_eigenpair_local_residual']<1e-9


def test_negative_real_multiplier_and_Gram_condition(saved,problem):
    b,_,item=saved;_,_,_,_,Omega,lam,_,_=problem;rec=n.p.read(b/'Floquet_reconstruction.json');mul=np.exp(lam*2*np.pi/Omega)
    assert mul.real<0 and abs(mul.imag)/abs(mul)<item['config']['real_multiplier_imaginary_relative_gate']
    assert rec['real_Gram_small_large_ratio']<item['config']['real_mode_Gram_ratio_gate']
    assert rec['real_subspace_rank']==1 and rec['phase_selected_before_ODE']
    assert rec['Fourier_indices']==list(range(-9,10))


def test_physical_OUT_scaling_and_preserved_speeds(saved,problem):
    b,_,item=saved;d,_,_,_,_,_,_,_=problem;rec=n.p.read(b/'Floquet_reconstruction.json')
    with np.load(b/'initial_states.npz') as z:
        state=z['spatial_state'];base=z['base_state'];unit=z['OUT_unit_state'];E=z['metric_OUT'];delta=rec['physical_delta_energy_coordinates']
        np.testing.assert_array_equal(state[d.state_plane],base[d.state_plane])
        np.testing.assert_allclose(state[d.state_out],delta*unit,atol=1e-20)
        assert np.sqrt(unit@E@unit)==pytest.approx(1,abs=2e-12)
        assert np.sqrt(state[d.state_out]@E@state[d.state_out])==pytest.approx(delta,rel=1e-12)
        assert np.max(abs(state[d.state_out[len(d.out):]]))>0
        _,scale,_=n.out_polynomial_extrema(d,state[:d.ndof],item['config']['Lref']);assert scale==pytest.approx(5e-5,abs=2e-15)


def test_mass_norm_fixed_and_reference_matrices_positive(saved,problem):
    b,_,_=saved;d,_,_,_,_,_,_,_=problem
    with np.load(b/'initial_states.npz') as z:
        np.testing.assert_array_equal(z['metric_OUT'],n.dynamic_metric(d,d.out))
        np.testing.assert_array_equal(z['metric_PL'],n.dynamic_metric(d,d.plane))
        assert np.linalg.eigvalsh(z['metric_OUT']).min()>0
        assert np.linalg.eigvalsh(z['metric_PL']).min()>0


def test_joint_clamps_and_reflection_invariance(saved,problem):
    b,_,item=saved;d,_,_,_,_,_,_,_=problem
    with np.load(b/'initial_states.npz') as z:x=z['spatial_state'].copy();base=z['base_state'].copy()
    reflection=np.ones(2*d.ndof);reflection[d.state_out]=-1
    assert np.max(abs(d.rhs(0,base)[d.state_out]))==0
    np.testing.assert_allclose(d.rhs(0,reflection*x),reflection*d.rhs(0,x),rtol=2e-12,atol=2e-12)
    joint=d.joint_diagnostics(x[:d.ndof],x[d.ndof:])
    assert max(joint[k] for k in ('translation_difference','orientation_difference','contraction_difference'))<item['config']['joint_absolute_gate']
    for arm,P in zip(d.arms,d.maps):np.testing.assert_allclose(arm.reconstruct(P@x[:d.ndof],np.array([0.])),0,atol=1e-14)


def test_independent_planar_period_energy_mass_time_gates(saved):
    b,_,item=saved;cfg=item['config'];r=n.p.read(b/'planar_one_period_checks.json');assert r['gate_B']
    assert r['period_not_fitted'] and r['planar_rhs_OUT_max']==0
    assert r['mass_min_eigenvalue']>0 and r['joint_compatibility_max']<cfg['joint_absolute_gate']
    for row in r['records']:
        assert row['HB_difference_energy_relative_max']<cfg['planar_orbit_energy_relative_gate']
        assert row['closure_energy_relative']<cfg['planar_orbit_energy_relative_gate']
        assert row['energy_drift']<cfg['energy_drift_relative_gate']
    assert r['temporal_difference_energy_relative_max']<cfg['temporal_energy_relative_gate']


def test_exact_native_restart_and_dense_join(saved,problem):
    b,s,_=saved;d,_,_,_,_,_,_,_=problem
    with np.load(b/'planar_refined_1T/accepted_times.npz') as z:native=z['actual_final_state'].copy();native_times=z['times'].copy()
    assert np.all(np.diff(native_times)>0) and native_times[0]==0
    with np.load(b/'planar_continue_T_to_3T/initial.npz') as z:np.testing.assert_array_equal(z['state'],native)
    join=n.p.read(b/'planar_combined_3T/join_provenance.json')
    assert join['physical_states_not_renormalized'] and join['prefix_not_recomputed'] and join['join_y_difference']<2e-12
    dense=a.SavedPlane(b/'planar_combined_3T/dense.npz',d)
    np.testing.assert_allclose(dense(join['join_t']),native,atol=2e-14)


def test_dynamic_tangent_Floquet_structure_and_time_control(saved):
    b,s,item=saved
    if not (b/'tangent_growth_checks.json').exists():pytest.skip('Tangent stage not completed')
    r=n.p.read(b/'tangent_growth_checks.json');cfg=item['config']
    expected=(max(r['stroboscopic_Floquet_vector_relative_errors'])<=cfg['Floquet_stroboscopic_relative_gate'] and
        r['sigma_difference']<=cfg['growth_exponent_absolute_gate'] and r['sigma_temporal_uncertainty']<=cfg['growth_temporal_absolute_gate'] and
        r['tangent_temporal_energy_error']<=cfg['temporal_energy_relative_gate'] and r['planar_HB_difference_3T_energy_relative_max']<=cfg['planar_orbit_energy_relative_gate'] and r['negative_multiplier_signs_correct'])
    assert r['gate_D']==expected
    if expected:
        assert r['sigma_dynamic']>0 and r['negative_multiplier_signs_correct']
        np.testing.assert_allclose(r['stroboscopic_amplification'],r['Hill_stroboscopic_amplification'],atol=.002)
    assert 'not HB arrays' in r['base_source']


def test_temporal_policies_declared_and_preserved(saved):
    b,s,item=saved;cfg=item['config']
    primary=n.p.read(b/'planar_primary_1T/time_settings.json');refined=n.p.read(b/'planar_refined_1T/time_settings.json')
    assert refined['rtol']==primary['rtol']*cfg['temporal_tolerance_factor']
    np.testing.assert_allclose(refined['atol'],np.array(primary['atol'])*cfg['temporal_tolerance_factor'],rtol=1e-14)
    assert primary['max_step']==refined['max_step']
    assert n.p.read(b/'pre_ODE_policy.json')['gates']==cfg


def test_independent_half_step_controls_keep_error_gates(saved):
    b,s,item=saved;cfg=item['config'];record=n.p.read(b/'temporal_step_control.json');assert record['gate']
    assert record['planar_energy_relative_difference']<cfg['temporal_energy_relative_gate']
    assert record['tangent_energy_difference']<cfg['temporal_energy_relative_gate']
    assert record['sigma_step_difference']<cfg['growth_temporal_absolute_gate']
    assert s['cases']['tangent_half_step_1T']['max_step']==.5*s['cases']['tangent_primary_3T']['max_step']


def test_full_nonlinear_tangent_consistency_and_physical_scale(saved):
    b,s,item=saved
    if not (b/'full_nonlinear_checks.json').exists():pytest.skip('Full stage not run by gate or budget')
    cfg=item['config'];r=n.p.read(b/'full_nonlinear_checks.json');rec=n.p.read(b/'Floquet_reconstruction.json')
    assert r['initial_physical_OUT_scale']==5e-5 and r['delta_energy_coordinates']==rec['physical_delta_energy_coordinates']
    with np.load(b/'full_tangent_comparison.npz') as z:
        E=np.load(b/'initial_states.npz')['metric_OUT'];err=a.energy_norm(z['OUT_nonlinear']-z['OUT_tangent'],E)
        np.testing.assert_allclose(z['OUT_error'],err,atol=1e-14)
        relative=float(err.max()/a.energy_norm(z['OUT_tangent'],E).max())
        assert relative==pytest.approx(r['OUT_tangent_relative_error_max'],rel=1e-12)
        assert z['physical_OUT_extrema'][:,0].max()==pytest.approx(r['maximum_v_Phi_psi'][0])
    assert r['maximum_v_over_A']==pytest.approx(r['maximum_v_Phi_psi'][0]/.005)


def test_conditional_5T_continuation_preserves_original_normalization(saved):
    b,s,_=saved
    if 'tangent_continue_3T_to_5T' not in s['cases']:pytest.skip('5T resource gate not passed')
    with np.load(b/'tangent_primary_3T/restart_state.npz') as z:last=z['actual_final_matrix'].copy();original=z['original_basis'].copy()
    with np.load(b/'tangent_continue_3T_to_5T/restart_state.npz') as z:
        np.testing.assert_array_equal(z['initial_matrix'],last);np.testing.assert_array_equal(z['original_basis'],original)
    r=n.p.read(b/'extension_5T_tangent_gate.json');assert r['gate']


def test_half_control_same_direction_and_frozen_phase(saved):
    b,s,_=saved
    if not (b/'half_amplitude_checks.json').exists():pytest.skip('Short half-scale resource gate not passed')
    r=n.p.read(b/'half_amplitude_checks.json');assert r['same_base_direction_phase_and_solver_policy'] and r['periods']==1
    with np.load(b/'initial_states.npz') as z:x=z['spatial_state'].copy();base=z['base_state'].copy()
    with np.load(b/'spatial_half_percent_1T/initial.npz') as z:half=z['state'].copy()
    np.testing.assert_allclose(half-base,.5*(x-base),atol=2e-20)
    assert r['order_status'] in ('RESOLVED','NUMERICAL_FLOOR_LIMITED')


def test_budget_actual_integrations_and_no_FEM(saved):
    b,s,item=saved
    assert s['numerical_seconds']<=item['config']['numerical_budget_seconds']
    assert s['calls']['FEM']==0 and s['calls']['nonlinear_ODE']<=8 and s['calls']['tangent_ODE']<=4
    for row in s['cases'].values():
        assert row['status'] in ('PASS','PARTIAL','STARTED')
        if row['status']=='PASS':assert row['actual_accepted_end']==row['target_end']
    if (b/'full_spatial_cost_forecast.json').exists():assert n.p.read(b/'full_spatial_cost_forecast.json')['no_new_ODE_for_profile']


def test_every_historical_artifact_and_physical_kernel_preserved(saved):
    _,_,item=saved;initial=item['initial_state']
    assert all(n.p.sha(n.p.ROOT/name)==sha for name,sha in initial['parent_file_hashes'].items())
    assert all(n.p.sha(n.p.ROOT/name)==item['code_hashes'][name] for name in n.old.f.s.r.KERNELS)
    for name in ('scripts/lib/nlsp_resonant_spatial_coupling_screening.py','data/input/nlsp_beta15_periodic_orbit_spatial_stability.json','docs/numerics/nlsp_beta15_periodic_orbit_spatial_stability.md'):
        assert n.p.sha(n.p.ROOT/name)==initial['initial_file_hashes'][name]


def test_no_HB_Hill_FEM_or_external_force_calls_in_dynamic_runner():
    forbidden={'solve_ivp','hill_spectrum','stability_compute','periodic_candidates','full_massform_tangents','body_force','run_gmsh','run_calculix'}
    for function in (n.dynamic_initial,n.dynamic_planar_check,n.dynamic_tangent_check,n.dynamic_full_check,n.dynamic_optional_controls):
        tree=ast.parse(inspect.getsource(function));calls={node.func.id if isinstance(node.func,ast.Name) else node.func.attr for node in ast.walk(tree)
            if isinstance(node,ast.Call) and isinstance(node.func,(ast.Name,ast.Attribute))}
        assert not calls&forbidden


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_cached_replay_on_temporary_copy_without_scientific_calls(saved,tmp_path,monkeypatch,mode):
    b,s,item=saved
    if not s['completed'] and not s['hard_stop']:pytest.skip('Scientific task still active')
    original=n.p.sha(b/'manifest.json');copied=tmp_path/'cache';shutil.copytree(b,copied)
    monkeypatch.setattr(n,'dynamic_prepare',lambda:(copied,n.p.read(copied/'provenance.json'),n.p.cache(copied)))
    def forbidden(*args,**kwargs):pytest.fail('Scientific call during cached replay')
    for owner,name in [(n,'dynamic_problem'),(n,'dynamic_initial'),(n,'dynamic_planar_check'),(n,'dynamic_tangent_check'),(n,'dynamic_temporal_step_control'),(n,'dynamic_full_check'),(n,'dynamic_optional_controls'),(n,'dynamic_postprocess'),
        (n.p,'sources'),(n.p,'make_disc'),(n.p,'integrate_attempt'),(a,'tangent_attempt')]:monkeypatch.setattr(owner,name,forbidden)
    n.dynamic_main([mode] if mode=='--compute' else [mode,str(copied)])
    assert n.p.sha(b/'manifest.json')==original
    assert n.p.cache(copied)['calls']==s['calls']
