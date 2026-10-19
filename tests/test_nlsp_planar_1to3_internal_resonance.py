"""Two-active-mode variational reduction, periodic residual and cache checks."""
import itertools, math
import numpy as np
import pytest
from scripts.lib import nlsp_planar_1to3_internal_resonance as n


@pytest.fixture(scope='module')
def saved():
    files=list(n.OUTPUT.glob('*/validated_periodic_candidates.json'))
    if not files:pytest.skip('Optional local fixed-beta15 diagnostic bundle is absent')
    b=files[0].parent
    return b,n.p.read(b/'canonical_resonant_coefficients.json'),n.p.read(b/'validated_periodic_candidates.json')


def decode(data,key):return {tuple(r['monomial']):complex(r['real'],r['imag']) for r in data[key]}


def test_saved_p12_frequencies_and_mass_orthogonality(saved):
    b,data,_=saved
    assert data['frequencies']==pytest.approx([.4348207304902722,1.3063293949106531],abs=1e-14)
    assert data['delta_13']==pytest.approx(.0018672034398365955,abs=1e-14)
    assert data['mass_orthogonality']<2e-10 and data['linear_eigenpair_residual']<2e-8
    assert data['passive_planar_modes_included']==90
    with np.load(b/'canonical_resonant_coefficients.npz') as z:assert z['Phi'].shape==(161,92)


def test_homological_sign_and_known_scalar_variational_limit():
    omega=1.7;alpha=.4;beta=.9;w=np.array([.9,omega,2.4,5.13]);H3={};H4={}
    for k in range(4):
        key=(2,)*(3-k)+(3,)*k;H3[tuple(sorted(key))]=alpha*math.comb(3,k)/(3*(2*omega)**1.5)
    for k in range(5):
        key=(2,)*(4-k)+(3,)*k;H4[tuple(sorted(key))]=beta*math.comb(4,k)/(4*(2*omega)**2)
    W,correction,_,gap=n.cubic_elimination(H3,w)
    coefficient=(H4[(2,2,3,3)]+correction[(2,2,3,3)]).real
    assert coefficient*omega==pytest.approx(3*beta/(8*omega)-5*alpha**2/(12*omega**3),rel=1e-13)
    for key,c in H3.items():
        D=sum((1 if i%2==0 else -1)*w[i//2] for i in key)
        assert 1j*D*W[key]+c==pytest.approx(0,abs=1e-14)
    assert gap>0


def test_projection_and_Legendre_terms_from_existing_action(saved):
    b,data,_=saved;checks=n.p.read(b/'canonical_variational_checks.json')
    assert abs(checks['self_shift_difference'])<1e-8
    assert all(r['H3_full_gradient_relative_difference']<2e-10 for r in checks['canonical_polynomial_projection_checks'])
    assert all(r['H4_absolute_polynomial_difference']<1e-20 for r in checks['canonical_polynomial_projection_checks'])
    decomposition=n.p.read(b/'coefficient_energy_decomposition.json')
    total=sum(v['gamma']['real'] for v in decomposition.values())
    assert total==pytest.approx(data['constants']['effective']['gamma']['real'],abs=1e-8)


def test_indirect_terms_and_conjugate_Hamiltonian_coefficients(saved):
    _,data,_=saved
    assert max(data['canonical_conjugacy_defects'].values())<2e-9
    assert data['constants']['direct']['gamma']['real']<0
    assert data['constants']['indirect']['gamma']['real']>0
    assert data['constants']['effective']['gamma']['real']==pytest.approx(-2858.5841701043296,rel=1e-12)
    assert data['minimum_cubic_elimination_frequency_denominator']>.04


def test_resonant_phase_invariance_and_J_identity(saved):
    _,data,_=saved;constants={k:v['real'] for k,v in data['constants']['effective'].items()}
    w=np.zeros(92);w[1],w[3]=data['frequencies'];z=np.zeros(184,complex)
    z[2]=.001*np.exp(-.7j);z[3]=z[2].conjugate();z[6]=.0004*np.exp(.2j);z[7]=z[6].conjugate()
    initial=n.normal_form_hamiltonian(z,w,constants);z[2]*=np.exp(-.6j);z[3]=z[2].conjugate();z[6]*=np.exp(-1.8j);z[7]=z[6].conjugate()
    assert n.normal_form_hamiltonian(z,w,constants)==pytest.approx(initial,abs=1e-18)
    # Every retained monomial commutes with J=I2+3I4 coefficient-exactly.
    for key in [(2,2,3,3),(2,3,6,7),(6,6,7,7),(2,2,2,7),(3,3,3,6)]:
        assert sum((1 if j%2==0 else -1)*(1 if j//2==1 else 3) for j in key)==0


def test_physical_amplitude_boundary_joint_and_time_reversal(saved):
    b,_,rows=saved;d,model,Phi,w=n.load_problem(n.p.read(b/'provenance.json'))
    for j,r in enumerate(rows):
        label=f'A{r["A"]:.3f}_branch{sum(t["A"]==r["A"] for t in rows[:j])+1}'
        with np.load(b/(label+'_raw.npz')) as z:U=z['U'].copy()
        q0=Phi@U.sum(axis=0)
        assert n.p.polynomial_w_max(d,q0)[0]==pytest.approx(r['A'],abs=1e-12)
        nidx=np.arange(len(U));phase=.73
        qplus=Phi@(np.cos(nidx*phase)@U);qminus=Phi@(np.cos(-nidx*phase)@U)
        vplus=-r['Omega']*Phi@((nidx*np.sin(nidx*phase))@U)
        vminus=-r['Omega']*Phi@((nidx*np.sin(-nidx*phase))@U)
        np.testing.assert_array_equal(qplus,qminus);np.testing.assert_array_equal(vplus,-vminus)
        assert r['joint_compatibility_max']<2e-10 and r['pointwise_mass_min_eigenvalue']>0
        for arm,P in zip(d.arms,d.maps):np.testing.assert_allclose(arm.reconstruct(P@q0,np.array([0.])),0,atol=1e-14)


def test_full_independent_phase_residual_and_harmonic_control(saved):
    _,_,rows=saved;cfg=n.p.read(n.CONFIG)
    for r in rows:
        assert r['independent_phase_count']==53 and r['phase_offset']==.173
        assert r['full_massform_residual_relative']<=cfg['full_planar_residual_relative_gate']
        assert r['energy_relative_variation']<=cfg['energy_relative_variation_gate']
        assert max(r['action_massform_difference_relative_to_equation_terms'])<=cfg['action_relative_gate']
        assert r['harmonic_control_frequency_difference']<=cfg['frequency_control_absolute_gate']
        assert r['one_passive_update_frequency_correction_indicator']<=cfg['frequency_control_absolute_gate']
        assert r['criterion_passed']


def test_1to3_resonant_slot_is_active_and_not_inverted(saved):
    _,_,rows=saved
    for r in rows:assert r['resonant_slots_not_inverted']==[[1,2],[3,4]]
    assert all(r['remaining_denominator_sorted_mode']!=4 or r['remaining_denominator_harmonic']!=3 for r in rows)
    main=[r for r in rows if r['phase_Psi']==0.]
    assert [r['A'] for r in main]==[.005,.006,.008,.010]
    assert main[-1]['linear_mode4_energy_fraction']>.5


def test_type_I_inverse_map_and_nonzero_resonant_source(saved):
    b,_,_=saved;rows=n.p.read(b/'type_I_initial_modulation_diagnostics.json')
    with np.load(b/'type_I_canonical_initial_states.npz') as z:
        for r in rows:
            state=z[f'A{r["A"]:.3f}'];assert state.shape==(184,)
            assert abs(state[6])**2==pytest.approx(r['canonical_I4_through_quadratic_map'],rel=1e-14)
            # Independent polynomial accumulation orders need the already
            # frozen canonical conjugacy gate, not bitwise float equality.
            defect=np.max(abs(state[1::2]-state[0::2].conj()))/max(np.max(abs(state)),1e-300)
            assert defect<n.p.read(n.CONFIG)['canonical_reality_relative_gate']
            assert r['physical_initial_Q4']==0 and r['pure_mode2_resonant_mode4_source_magnitude']>0
            assert r['passive_initial_reference_energy_from_inverse_map']>0


def test_historical_kernels_and_single_mode_results_unchanged(saved):
    b,_,_=saved;cfg=n.p.read(n.CONFIG);assert n.p.sha(n.p.ROOT/cfg['planar_parent']/'manifest.json')==cfg['planar_parent_manifest_sha256']
    initial=n.p.read(b/'provenance.json')
    assert all(n.p.sha(n.p.ROOT/name)==initial['code_hashes'][name] for name in n.old.f.s.r.KERNELS)
    assert n.old.cache_saved(n.p.ROOT/cfg['planar_parent'])['overall']=='PLANAR_SINGLE_MODE_EXPANSION_UNRESOLVED'


def test_scope_and_numerical_gates_are_frozen():
    cfg=n.p.read(n.CONFIG);assert n.validate_scope(cfg)==cfg
    for key,value in [('degree',16),('beta_deg',30),('amplitudes',[.005]),('OUT_reduction',True),('numerical_budget_seconds',601),('full_planar_residual_relative_gate',.1)]:
        with pytest.raises(ValueError):n.validate_scope({**cfg,key:value})


@pytest.mark.parametrize('mode',['--compute','--report-only','--plot-only'])
def test_real_cached_replay_without_scientific_calls(saved,monkeypatch,mode,tmp_path):
    b,_,_=saved
    import shutil
    original=n.p.sha(b/'manifest.json');copied=tmp_path/'historical_cache_copy'
    shutil.copytree(b,copied);item=n.p.read(copied/'provenance.json');cached=n.p.cache(copied)
    monkeypatch.setattr(n,'prepare',lambda:(copied,item,cached))
    def forbidden(*a,**k):pytest.fail('Scientific call in cached report/plot')
    for owner,name in [(n,'compute_coefficients'),(n,'periodic_candidates'),(n,'verify_and_control'),(n,'complete_diagnostics'),
        (n,'frequency_sensitivity_check'),(n,'load_problem'),(n.p,'sources'),(n.p,'make_disc'),(n.p,'integrate_attempt'),
        (n.old.s.l.a,'tangent_attempt'),(n.old.s.c,'limited_modes')]:monkeypatch.setattr(owner,name,forbidden)
    n.main([mode] if mode=='--compute' else [mode,str(copied)])
    summary=n.p.read(copied/'summary.json');assert summary['calls']['ODE']==summary['calls']['FEM']==summary['calls']['OUT']==0
    assert n.p.sha(b/'manifest.json')==original
