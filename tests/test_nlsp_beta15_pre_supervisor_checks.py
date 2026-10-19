"""Fixed saved-orbit checks; no production integration during regression."""
import ast
import inspect
import shutil
import numpy as np
import pytest
from scripts.lib import nlsp_planar_1to3_internal_resonance as n


@pytest.fixture(scope='module')
def saved():
    b=n.SUPERVISOR_OUTPUT/'d3f6375a27dba883'
    if not (b/'manifest.json').exists():pytest.skip('Final-check bundle unavailable')
    return b,n.p.read(b/'summary.json'),n.p.read(b/'provenance.json')


@pytest.fixture(scope='module')
def problem(saved):return n.supervisor_problem(saved[2])


def test_frozen_sources_and_all_historical_files(saved):
    b,s,item=saved
    for name,sha in item['initial_state']['physical_kernel_hashes'].items():assert n.p.sha(n.p.ROOT/name)==sha
    for name,sha in item['initial_state']['historical_file_hashes'].items():assert n.p.sha(n.p.ROOT/name)==sha
    for source in item['source_manifests'].values():assert n.p.sha(n.p.ROOT/source['path']/'manifest.json')==source['sha256']
    assert item['config']['degree']==14 and item['config']['A']==.005 and item['config']['phase_Psi']==0
    parent=n.p.ROOT/item['config']['p14_parent'];row=n.p.read(parent/'OUT_checks.json')
    assert row['sigma']==.003411315507705602
    assert n.p.read(parent/'HB13.json')['Omega']==.43511545830378884
    old=n.p.read(n.p.ROOT/item['config']['p12_Hill_parent']/'summary.json')['cases'][0]
    assert old['OUT_sigma']==.0034114100665946817 and old['Omega']==.4351155867232466


def test_physical_OUT2_normalization_and_zero_velocities(saved,problem):
    b,s,item=saved;d,*_=problem
    with np.load(b/'OUT2_initial_components.npz') as z:x=z['total'].copy()
    np.testing.assert_array_equal(x[len(d.out):],np.zeros(len(d.out)))
    q=np.zeros(d.ndof);q[d.out]=x[:len(d.out)]
    assert n.out_polynomial_extrema(d,q,1)[1]==pytest.approx(5e-5,abs=2e-15)
    assert max(d.joint_diagnostics(q,np.zeros(d.ndof))[key] for key in ('translation_difference','orientation_difference','contraction_difference'))<2e-12


def test_symplectic_duality_and_saved_projection(saved):
    b,_,item=saved;row=n.p.read(b/'OUT2_prediction.json')
    with np.load(b/'OUT2_initial_components.npz') as z:
        ru,rs,S,J,E,x=(z[key] for key in ('ru','rs','S','J','metric','total'))
        left=(S@rs)@J/((S@rs)@J@(S@ru))
        assert left@S@ru==pytest.approx(1,abs=2e-12) and abs(left@S@rs)<2e-12
        assert ru@E@ru==pytest.approx(1,abs=2e-12)
        assert left@S@x==pytest.approx(-1.8779480638240486e-7,rel=2e-12)
        np.testing.assert_allclose(z['unstable'],ru*(left@S@x),rtol=2e-12,atol=1e-16)
        np.testing.assert_allclose(x,z['unstable']+z['stable']+z['remaining'],atol=1e-20)
    assert abs(row['projection_reference_difference'])<1e-18 and row['not_energy_fractions']


@pytest.mark.parametrize('factor',[2.,-.7,3.+1.2j])
def test_dual_physical_projector_is_invariant_to_vector_normalization(saved,factor):
    b,_,_=saved
    with np.load(b/'OUT2_initial_components.npz') as z:
        S,J,x=z['S'],z['J'],z['total'];ru=factor*z['ru'];rs=z['rs']/factor
        left=(S@rs)@J/((S@rs)@J@(S@ru))
        np.testing.assert_allclose(ru*(left@S@x),z['unstable'],rtol=2e-12,atol=1e-18)


def test_neutral_cause_and_unchanged_failed_raw_gate(saved):
    b,_,item=saved;row=n.p.read(b/'neutral_diagnostics.json');gate=item['config']['neutral_relative_gate']
    assert gate==row['original_gate_unchanged']==2e-9
    assert row['dominant_residual_derivative_harmonic']==12
    assert row['raw_source_N13_residual']==7.585992727452457e-8
    assert row['variational_identity_gate'] and not row['raw_exact_periodic_neutral_gate']
    assert row['PL_status']=='PL_STABILITY_REMAINS_PARTIAL'
    assert row['saved_near_zero_PL_mode_time_shift_MAC']>.99999999
    assert row['saved_non_neutral_sigma']<2e-9
    for value in row['Fourier_controls']:
        if value['N']>=13:
            assert value['raw_Hill_neutral_relative_residual']>gate
            assert value['orbit_residual_corrected_Hill_relative_error']<gate
    # N=11 misses harmonic 12 and cannot close the unchanged physical gate.
    assert row['Fourier_controls'][1]['raw_Hill_neutral_relative_residual']<gate


def test_saved_variational_identity_including_canonical_time_derivative(saved):
    b,_,_=saved;diag=n.p.read(b/'neutral_diagnostics.json')
    with np.load(b/'neutral_audit.npz') as z:
        scale=diag['neutral_scale_Omega_RMS_Z']
        error=z['canonical_time_shift_residual']-z['mapped_first_order_residual_derivative']-z['canonical_off_orbit_difference']
        assert np.linalg.norm(error,axis=1).max()/scale<2e-9
        assert np.linalg.norm(z['identity_errors'],axis=1).max()/scale<2e-9
        ids=abs(z['Fourier_orders'])>=13
        assert np.linalg.norm(z['residual_derivative_coefficients'][ids])>0
    assert max(r['first_order_residual_derivative_identity_relative_error_max'] for r in diag['independent_FD_controls'])<2e-9


def test_original_first_order_residual_identity_at_new_independent_phase(saved,problem):
    import copy
    b,_,item=saved;d,model,Phi,U,Omega,*_=problem;phase=.437;h=.001
    parent=n.p.ROOT/item['config']['p14_parent']
    with np.load(parent/'PL_Fourier.npz') as z:L=z['L'].copy()
    holder=copy.copy(d);holder.out=d.plane;action=n.old.f.QuadraticAction(holder,model.T4,model.V4)
    q,v,a,j=n.orbit_with_jerk(U,Omega,Phi,np.array([phase]))
    qp,vp,ap,_=n.orbit_with_jerk(U,Omega,Phi,np.array([phase+h]))
    qm,vm,am,_=n.orbit_with_jerk(U,Omega,Phi,np.array([phase-h]))
    direct=np.r_[a[0],j[0]]-d.jacobian(0,np.r_[q[0],v[0]])@np.r_[v[0],a[0]]
    fd=np.r_[np.zeros(d.ndof),(ap[0]-d.acceleration(qp[0],vp[0])-am[0]+d.acceleration(qm[0],vm[0]))*Omega/(2*h)]
    S=n.old.f.action_operators(action.blocks(q[0],v[0],a[0]))[3]
    scale=n.p.read(b/'neutral_diagnostics.json')['neutral_scale_Omega_RMS_Z']
    assert np.linalg.norm(L@S@(direct-fd)[d.state_plane])/scale<2e-9


def test_Fourier_conventions_real_operator_and_multiplier_forecast(saved):
    b,_,item=saved;parent=n.p.ROOT/item['config']['p14_parent'];Omega=n.p.read(parent/'HB13.json')['Omega']
    with np.load(parent/'OUT_Fourier.npz') as z:F=z['F'].copy();orders=z['orders'].copy()
    operator=n.SavedCanonicalFourier(F,orders,Omega)
    for phase in (.127,1.413,4.731):
        exact=np.einsum('m,mij->ij',np.exp(1j*orders*phase),F)
        assert np.linalg.norm(exact-operator(phase/Omega))/np.linalg.norm(exact)<2e-12
    diag=n.p.read(b/'OUT2_prediction.json');rho=n.p.read(parent/'OUT_checks.json')['multiplier']
    multiplier=complex(rho['real'],rho['imag']);assert multiplier.real<0
    np.testing.assert_allclose(diag['unstable_stroboscopic_gains'],abs(multiplier)**np.array([3,5,10,20]),rtol=2e-12)
    assert diag['forecast_10T_20T_not_integrated'] and diag['forecast_not_a_total_OUT_response']


def test_one_short_tangent_control_checkpoint_and_observables(saved):
    b,s,item=saved;diag=n.p.read(b/'OUT2_tangent_checks.json');policy=n.p.read(b/'OUT2_tangent_policy.json')
    assert policy['OUT_velocities_zero'] and policy['pointwise_relative_Jacobian_comparison']<2e-9
    assert diag['one_selected_trajectory_with_conditional_restart']
    assert s['calls']=={'nonlinear_ODE':0,'tangent_ODE':2,'FEM':0,'HB_Hill_eigensolves':0}
    assert s['cases']['OUT2_continue_3T_to_5T']['original_initial_basis_preserved']
    assert diag['growth_projection_relative_error']<item['config']['projection_relative_gate']
    assert diag['maximum_v_over_A_through_5T']==pytest.approx(.002863587272859495)
    assert diag['stroboscopic_growth_projection_amplification'][-1]==pytest.approx(1.2792854285835593)
    assert diag['maximum_total_E0_amplification']==pytest.approx(1.2179856168644745)
    with np.load(b/'OUT2_tangent_physical.npz') as z:
        assert len(z['times'])==1281 and np.all(np.diff(z['times'])>0)
        assert z['total_E0_norm'].max()/z['total_E0_norm'][0]==pytest.approx(diag['maximum_total_E0_amplification'])
    assert s['numerical_seconds']<600


def test_no_new_full_nonlinear_HB_Hill_FEM_or_expanded_scope(saved):
    forbidden={'solve_ivp','integrate_attempt','Radau','DOP853','run_gmsh','run_calculix','hb_newton','hill_spectrum'}
    for fun in (n.supervisor_problem,n.supervisor_neutral_audit,n.supervisor_projection,n.supervisor_tangent_control,n.supervisor_report_data):
        tree=ast.parse(inspect.getsource(fun))
        calls={node.func.id if isinstance(node.func,ast.Name) else node.func.attr for node in ast.walk(tree)
               if isinstance(node,ast.Call) and isinstance(node.func,(ast.Name,ast.Attribute))}
        assert not calls&forbidden
    cfg=saved[2]['config'];assert cfg['new_full_nonlinear_calls']==cfg['new_FEM_jobs']==cfg['new_HB_Hill_eigensolves']==0


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_completed_cache_replay_on_temporary_copy_without_science(saved,tmp_path,monkeypatch,mode):
    b,s,item=saved;before=n.p.sha(b/'manifest.json');copy=tmp_path/'bundle';shutil.copytree(b,copy)
    monkeypatch.setattr(n,'supervisor_prepare',lambda:(copy,n.p.read(copy/'provenance.json'),n.p.cache(copy)))
    def forbidden(*args,**kwargs):pytest.fail('Scientific call during completed replay')
    for name in ('supervisor_problem','supervisor_neutral_audit','supervisor_projection','supervisor_tangent_control','supervisor_report_data'):
        monkeypatch.setattr(n,name,forbidden)
    n.supervisor_main([mode] if mode=='--compute' else [mode,str(copy)])
    assert n.p.sha(b/'manifest.json')==before and n.p.cache(copy)['calls']==s['calls']


def test_report_profiles_are_physical_and_geometry_is_not_an_acute_V(saved):
    b,_,_=saved
    with np.load(b/'report_data.npz') as z:
        assert z['rest_centerlines'][0,0,0]<0 and z['rest_centerlines'][1,0,0]>0
        assert z['rest_centerlines'][1,0,1]>0
        assert z['bases'][0,2,2]==z['bases'][1,2,2]==-1
        np.testing.assert_allclose(z['PL_displacement'][:,:,2],0,atol=1e-15)
        assert np.max(abs(z['OUT_displacement'][:,:,2]))>0
    report=n.p.ROOT/'docs/reports/nlsp_beta15_parametric_instability/report.tex'
    assert report.exists() and len(report.read_text(encoding='utf-8'))>7000
    assert 'Тип' in report.read_text(encoding='utf-8') or 'тип' in report.read_text(encoding='utf-8')
