"""D15 selection, normalization and saved-result checks; no repeated root solve."""
import numpy as np
import pytest

from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_targeted_weak_damping as run


@pytest.fixture(scope='module')
def saved():
    if not (run.OUTPUT/'diagnostics.json').exists():
        pytest.skip('local D15 artifacts unavailable')
    return run.screen.read_json(run.OUTPUT/'diagnostics.json')


@pytest.fixture(scope='module')
def seeds():
    return run.select_seeds(run.screen.read_csv(run.SOURCE/'elastic_screening.csv'))


def test_three_predetermined_active_seeds_and_six_requested_targets(seeds):
    assert [(s['state_id'], s['beta_deg'], s['elastic_sorted_mode']) for s in seeds] == list(run.SELECTION)
    assert [(s['beta_deg'], s['elastic_sorted_mode']) for s in seeds] == [(0., 'sorted_05'), (45., 'sorted_02'), (75., 'sorted_05')]
    assert all(s['source_row']['activity_status'] == 'ACTIVE' for s in seeds)
    assert run.D_VALUES == (.001, .005)
    assert len({(s['state_id'], d) for s in seeds for d in run.D_VALUES}) == 6
    assert run.ATTEMPTS == 2 and run.BUDGET == 1000


def test_inactive_or_other_seed_is_not_substituted():
    rows = run.screen.read_csv(run.SOURCE/'elastic_screening.csv')
    target = next(r for r in rows if r['beta_deg'] == '0.0' and r['sorted_mode'] == 'sorted_05')
    target['activity_status'] = 'EXACT_INACTIVE_BY_SYMMETRY'
    with pytest.raises(ValueError, match='ACTIVE'):
        run.select_seeds(rows)


@pytest.mark.parametrize('d,expected',[(.001,1.443375672974065e-7),(.005,7.216878364870325e-7)])
def test_dimensional_conversion_and_eb_only(seeds, d, expected, monkeypatch):
    p = run.kv.section()[1]
    original = run.kv.Arm.reduced
    def eb_only(model, properties, L=1.):
        assert model == 'EB'
        return original(model, properties, L)
    monkeypatch.setattr(run.kv.Arm, 'reduced', eb_only)
    provider = run.provider(p, seeds[0], d, run.kv.Calls())
    assert provider.c == pytest.approx(expected, rel=1e-14)
    assert provider.k == run.kv.M_REF
    assert all(a.J == a.invS == 0 for a in provider.arms)
    with pytest.raises(ValueError, match='predetermined'):
        run.provider(p, seeds[0], .002, run.kv.Calls())


def test_predictors_use_elastic_slope_then_previous_complex_root(seeds):
    s = seeds[0]
    expected = -s['a_slope_pred']*.001 + 1j*s['Omega0']
    assert run.predictor(s, .001) == expected
    previous = -.2 + 76.5j
    assert run.predictor(s, .005, previous) == previous
    with pytest.raises(ValueError):
        run.predictor(s, .005)
    with pytest.raises(ValueError):
        run.predictor(s, .003)


def test_observables_against_independent_arithmetic(seeds):
    s = seeds[1]; d = .001
    z = -.021 + 1j*(s['Omega0']+.002)
    r = run.observables(s, d, z)
    expected_zeta = .021/np.sqrt(.021**2+(s['Omega0']+.002)**2)
    assert r['zeta'] == pytest.approx(expected_zeta)
    assert r['a_over_d'] == pytest.approx(21.)
    assert r['zeta_over_d'] == pytest.approx(expected_zeta/d)
    assert r['relative_frequency_shift'] == pytest.approx(.002/s['Omega0'])
    assert r['frequency_shift_over_d2'] == pytest.approx(.002/s['Omega0']/d**2)
    assert r['omega_d'] == z.imag/run.kv.T_REF
    assert r['omega_d'] != abs(z/run.kv.T_REF)
    assert run.observables(s, .005, z)['weak_damping_status'] == 'DESCRIPTIVE_ONLY'


def test_full_matrix_and_conjugate_partner_of_saved_candidates(saved, seeds):
    properties = run.kv.section()[1]
    with np.load(run.OUTPUT/'complex_shapes.npz') as archive:
        for key, point in saved['points'].items():
            row = point['row']; seed = next(s for s in seeds if s['state_id'] == row['state_id'])
            provider = run.provider(properties, seed, row['d_theta'], run.kv.Calls())
            z = complex(row['z_re'], row['z_im'])
            B = provider.matrices(z)[0]; conjugate = provider.matrices(z.conjugate())[0]
            a = archive[key+'__a']
            assert B.shape == (6, 6)
            assert np.linalg.norm(B@a)/(np.linalg.norm(B)*np.linalg.norm(a)) < run.kv.CRITERIA['null_residual']
            np.testing.assert_allclose(conjugate, B.conj(), rtol=run.kv.CRITERIA['transfer_rtol'], atol=run.kv.CRITERIA['matrix_atol'])
            assert z.imag > 0
            # Full-B residual alone cannot erase the recorded physical gate.
            if not point['accepted']:
                assert 'PHYSICAL_GATE' in row['failures']


def test_energy_identity_and_alpha_from_found_complex_form(saved):
    for point in saved['points'].values():
        r = point['row']; p = complex(r['p_re'], r['p_im'])
        M, K, C = r['M_phi'], r['K_phi'], r['C_phi']
        residual = abs(p*p*M+p*C+K)/(abs(p)**2*M+abs(p)*C+K)
        assert residual == pytest.approx(r['energy_residual'], rel=1e-6, abs=1e-15)
        assert residual <= run.kv.CRITERIA['energy_residual']
        assert r['alpha_energy'] == pytest.approx(C/(2*M), rel=1e-14)
        assert abs(r['alpha']-C/(2*M))*run.kv.T_REF <= run.kv.CRITERIA['a_atol'] + run.kv.CRITERIA['a_rtol']*max(r['a'], C/(2*M)*run.kv.T_REF)


def test_mac_and_mass_use_corresponding_seed_not_frequency_nearest(saved, seeds):
    properties = run.kv.section()[1]
    with np.load(run.SOURCE/'screening_shapes.npz') as elastic, np.load(run.OUTPUT/'complex_shapes.npz') as complex_shapes:
        for key, point in saved['points'].items():
            r = point['row']; s = next(s for s in seeds if s['state_id'] == r['state_id'])
            provider = run.provider(properties, s, r['d_theta'], run.kv.Calls())
            _, weights = run.kv.quadrature()
            left = run.kv.mass_vector(elastic[s['shape_key']+'__states'], provider.arms, weights)
            y = complex_shapes[key+'__states']
            right = run.kv.mass_vector(y, provider.arms, weights)
            MAC = abs(np.vdot(left, right))**2/(np.vdot(left,left).real*np.vdot(right,right).real)
            assert MAC == pytest.approx(r['MAC'], abs=1e-14)
            assert np.vdot(right, right).real == pytest.approx(1., abs=1e-12)
            assert r['mode_status'] == 'MODE_CONFIRMED' and MAC >= .95
            assert r['abs_Delta_psi'] == pytest.approx(abs(y[0,-1,2]-y[1,-1,2]))


def test_failed_gates_retained_after_bounded_retry_and_branch_stops(saved):
    assert saved['principal_positive_d_roots_requested'] == 6
    assert saved['principal_positive_d_roots_computed'] == 4
    assert saved['principal_positive_d_roots_accepted'] == 2
    assert set(saved['principal_targets_not_run']) == {'A_STRONG_d2', 'C_WEAK_ACTIVE_d2'}
    for sid in ('A_STRONG', 'C_WEAK_ACTIVE'):
        point = saved['points'][sid+'_d1']
        assert len(point['attempts']) == 2 and not point['accepted']
        assert point['prior_record']['row']['solver_status'] == 'NUMERICAL_UNRESOLVED'
        assert point['row']['weak_damping_status'] == 'WEAK_DAMPING_CLOSE'
        assert point['row']['solver_status'] == 'NUMERICAL_UNRESOLVED'
        assert point['row']['K12_gate_status'] == 'QUALIFIED'
    assert saved['calls']['B']+saved['calls']['B_z'] <= 1000
    assert all(a['correction']['steps'] <= run.kv.CRITERIA['max_steps'] for p in saved['points'].values() for a in p['attempts'])


def test_no_inactive_rlb_or_auxiliary_roots(saved):
    for key in ('inactive_roots','RLB_roots','additional_beta','additional_d_values','auxiliary_positive_d_roots'):
        assert saved[key] == 0
    assert all(p['row']['model'] == 'EB' and p['d_theta'] in run.D_VALUES for p in saved['points'].values())
    assert saved['K12_criteria'] == run.kv.CRITERIA


@pytest.mark.parametrize('retry', [False, True])
def test_missing_only_and_exhausted_retry_make_no_computational_calls(saved, monkeypatch, retry):
    def forbidden(*args, **kwargs):
        pytest.fail('completed/exhausted targets must not be recomputed')
    for name in ('correct','recover','diagnose'):
        monkeypatch.setattr(run.kv, name, forbidden)
    monkeypatch.setattr(run.kv.Provider, 'matrices', forbidden)
    before = {p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    result = run.compute(retry_unresolved=retry)
    assert result['new_principal_roots'] == result['B'] == result['B_z'] == result['form_recoveries'] == 0
    assert before == {p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    assert saved['protected_sources'] == run.protected_hashes()
