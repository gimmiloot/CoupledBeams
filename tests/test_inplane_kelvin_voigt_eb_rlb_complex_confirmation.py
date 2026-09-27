"""D22 scope, saved physical checks and read-only reproducibility; no roots."""
import ast
import inspect
import math
import numpy as np
import pytest
from scripts.analysis.laminated_beams import confirm_inplane_kelvin_voigt_eb_rlb as run


@pytest.fixture(scope='module')
def saved():
    if not (run.OUTPUT/'diagnostics.json').exists():
        pytest.skip('local D22 results unavailable; do not recreate old data')
    return run.screen.read_json(run.OUTPUT/'diagnostics.json')


@pytest.fixture(scope='module')
def inputs():
    return run.inputs()


def test_exact_scope_and_cost(saved):
    assert run.CASES == (('R0',5.,'01'),('R1',0.,'05'),('R2',45.,'05'),('R3',75.,'05'),('R4',45.,'03'))
    assert run.NEW_TARGETS == ('R1_RLB','R2_EB','R2_RLB','R3_RLB','R4_EB','R4_RLB')
    assert tuple(saved['points']) == run.NEW_TARGETS
    assert saved['new_principal_roots'] == saved['new_roots_accepted'] == 6
    assert saved['reused_complex_roots'] == 4
    assert saved['criteria'] == run.kv.CRITERIA
    assert saved['ratio_descriptive_rtol'] == .01
    for key in ('auxiliary_roots','full_newton_roots','new_beta','new_d_outside_authorized',
                'asymmetric_positive_d','elastic_recomputations','solver_changes','high_precision'):
        assert saved[key] == 0
    assert saved['calls']['total_build_equivalents'] <= run.BUDGET
    assert saved['full_control_evaluations'] == 6


@pytest.mark.parametrize('key',run.NEW_TARGETS)
def test_matched_seed_routing_and_actual_attempts(saved,inputs,key):
    seeds,matches = inputs; s = seeds[key]; m = matches[s['case_id']]
    props = run.kv.section()[1]
    cfg = run.config_for(props,s)
    assert cfg.identical and cfg.mu == 0 and cfg.d_theta == .001
    assert cfg.arms[0].L == cfg.arms[1].L == 1
    assert run.production.route(cfg,'reduced') == 'SYMMETRY_REDUCED'
    assert s['eta'] == 1 and s['Omega0'] == float(m['Omega_'+s['theory']])
    assert s['G'] == float(m['G_'+s['theory']])
    assert m['match_status'] == 'CONFIRMED'
    for bad_d in (.005,.002,.0005,0.):
        with pytest.raises(ValueError): run.config_for(props,s,bad_d)
    for bad in (dict(s,beta_deg=46.),dict(s,eta=-1),dict(s,matched_mode='sorted_02')):
        with pytest.raises(ValueError): run.config_for(props,bad)
    attempts = saved['points'][key]['attempts']
    assert 1 <= len(attempts) <= 2
    for attempt in attempts:
        assert attempt['d_theta'] == .001 and attempt['eta'] == 1 and attempt['solver_path'] == 'reduced'
        r = attempt['result']
        assert r['solver_path'] == 'SYMMETRY_REDUCED' and r['eta'] == 1
        assert r['complex_newton_calls'] == 1 and r['model'] == s['theory']


def test_four_reused_states_are_exact_read_only(saved,inputs):
    seeds,_ = inputs
    rows,source = run.reused_rows(seeds)
    assert rows == saved['reused_rows'] and source == saved['reused_source_rows']
    assert [(r['case_id'],r['theory'],r['source']) for r in rows] == [
        ('R0','EB','K12'),('R0','RLB','K12'),('R1','EB','K19'),('R3','EB','K19')]
    assert rows[0]['d_theta'] == rows[1]['d_theta'] == .0005423772776686932
    for r in rows:
        s = source[r['case_id']+'_'+r['theory']]
        suffix = ('real','imag') if r['source']=='K12' else ('re','im')
        assert r['z_re'] == float(s['z_'+suffix[0]]) and r['z_im'] == float(s['z_'+suffix[1]])
        assert r['p_re'] == float(s['p_'+suffix[0]]) and r['p_im'] == float(s['p_'+suffix[1]])
        if r['source']=='K19':
            assert 'K16 historical full recovery rejected' in r['notes']
            assert 'K17/K18' in r['notes']
        with pytest.raises(ValueError): run.config_for(run.kv.section()[1],seeds[r['case_id']+'_'+r['theory']])


@pytest.mark.parametrize('key',run.NEW_TARGETS)
def test_saved_new_roots_forms_and_energy(saved,inputs,key):
    point = saved['points'][key]; r = point['row']; attempt = point['attempts'][-1]
    result = attempt['result']; diag = result['diagnostics']; control = attempt['full_control']
    assert result['accepted'] and not result['failures'] and r['root_status'] == 'ROOT_ACCEPTED'
    gates = run.kv.CRITERIA
    for field,gate in (('root_residual','null_residual'),('sigma_ratio','sigma_ratio'),
        ('reduced_physical_residual','physical_residual'),('lifted_full_physical_residual','physical_residual'),
        ('energy_residual','energy_residual'),('conjugate_residual','null_residual')):
        assert r[field] <= gates[gate]
    assert r['MAC_to_elastic'] >= gates['MAC']
    assert control['newton_calls'] == 0 and control['calls']['corrections'] == 0
    assert control['calls']['recoveries'] == 1
    assert control['z'] == result['z']
    assert control['status'] == 'PASS'  # actual six controls, separate from primary acceptance
    seed = inputs[0][key]; cfg = run.config_for(run.kv.section()[1],seed)
    with np.load(run.OUTPUT/'new_complex_shapes.npz',allow_pickle=False) as archive:
        y = archive[key+f"_attempt{attempt['attempt']}__states"]
    np.testing.assert_array_equal(y[1],y[0]@run.sd.F)
    np.testing.assert_array_equal(y[:,0,:3],np.zeros((2,3)))
    _,w = run.kv.quadrature()
    mass = sum(a.L*np.dot(w,a.m*(abs(v[:,0])**2+abs(v[:,1])**2)+a.J*abs(v[:,2])**2) for a,v in zip(cfg.arms,y))
    assert mass == pytest.approx(1.,abs=1e-12) and mass == pytest.approx(r['M_phi'],abs=1e-14)
    delta = y[0,-1,2]-y[1,-1,2]
    assert delta == pytest.approx(complex(r['Delta_psi_re'],r['Delta_psi_im']),rel=1e-14)
    stiffness = sum(a.L*np.dot(w,abs(v[:,3])**2/a.A+abs(v[:,5])**2/a.D+a.invS*abs(v[:,4])**2) for a,v in zip(cfg.arms,y))+run.kv.M_REF*abs(delta)**2
    damping = r['c_theta']*abs(delta)**2
    assert stiffness == pytest.approx(r['K_phi'],rel=1e-13)
    assert damping == pytest.approx(r['C_phi'],rel=1e-13)
    p = complex(r['p_re'],r['p_im'])
    energy = abs(p*p*mass+p*damping+stiffness)/(abs(p)**2*mass+abs(p)*damping+stiffness)
    assert energy <= gates['energy_residual']
    assert abs(-p.real-damping/(2*mass)) <= gates['a_atol']/run.kv.T_REF+gates['a_rtol']*abs(p.real)
    v = run.kv.mass_vector(y,cfg.arms,w); v0 = run.kv.mass_vector(seed['states'],cfg.arms,w)
    assert float(run.kv.mac_matrix([v0],[v])[0,0]) == pytest.approx(r['MAC_to_elastic'],abs=1e-14)
    if r['theory']=='RLB':
        assert cfg.arms[0].J > 0 and cfg.arms[0].invS > 0
    else:
        assert cfg.arms[0].J == cfg.arms[0].invS == 0


@pytest.mark.parametrize('cid',[c[0] for c in run.CASES])
def test_comparison_ratios_and_frequency(saved,cid):
    rows = run.comparisons(saved['all_states'],saved['matching_source_rows'])
    assert len(rows) == 5
    r = next(r for r in rows if r['case_id']==cid)
    e,b = [next(s for s in saved['all_states'] if s['case_id']==cid and s['theory']==t) for t in ('EB','RLB')]
    assert r['G_ratio'] == pytest.approx(r['G_RLB']/r['G_EB'],rel=1e-15)
    assert r['zeta_ratio'] == b['zeta']/e['zeta']
    assert r['ratio_error'] == (r['zeta_ratio']/r['G_ratio']-1) or math.isclose(r['ratio_error'],r['zeta_ratio']/r['G_ratio']-1,abs_tol=2e-16)
    assert r['elastic_delta_G'] == pytest.approx(r['G_ratio']-1,abs=2e-16)
    assert r['complex_delta_zeta'] == r['zeta_ratio']-1
    assert r['elastic_delta_Omega'] == (r['Omega0_RLB']-r['Omega0_EB'])/r['Omega0_EB']
    assert r['damped_delta_Omega'] == (b['Omega_d']-e['Omega_d'])/e['Omega_d']
    for s in (e,b):
        assert s['Omega_d'] == s['z_im'] and s['Omega_d'] != abs(complex(s['z_re'],s['z_im']))
        assert s['zeta_over_d'] == s['zeta']/s['d_theta']
        assert s['relative_predictor_error'] == (s['zeta_over_d']-s['G'])/s['G']
    # Numerical validity is stored independently of this descriptive label.
    changed = [dict(s,zeta=s['zeta']*2) if s['case_id']==cid and s['theory']=='RLB' else s for s in saved['all_states']]
    altered = next(x for x in run.comparisons(changed,saved['matching_source_rows']) if x['case_id']==cid)
    assert altered['comparison_status'] == 'DESCRIPTIVE_DEVIATION'
    assert altered['eb_root_status'] == altered['rlb_root_status'] == 'ROOT_ACCEPTED'


def test_no_full_newton_in_control_or_orchestration():
    source = inspect.getsource(run.completion.full_control)
    assert 'solve_mode(' not in source and 'kv.correct(' not in source
    tree = ast.parse(inspect.getsource(run.compute))
    calls = [n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='solve_mode']
    assert len(calls)==1
    args = {v.arg:ast.literal_eval(v.value) for v in calls[0].keywords if v.arg in ('eta','solver_path')}
    assert args == dict(eta=1,solver_path='reduced')


def test_source_bytes_and_no_solver_changes(saved):
    assert saved['protected_sources_unchanged']
    assert all(run.screen.sha(run.ROOT/p)==h for p,h in saved['protected_sources'].items())
    assert any(p.endswith('inplane_kelvin_voigt_solver.py') for p in saved['protected_sources'])


def test_missing_only_zero_scientific_calls_and_unchanged_files(saved,monkeypatch):
    def forbidden(*a,**k): pytest.fail('missing-only must not perform scientific work')
    for obj,name in ((run.production,'solve_mode'),(run.production.FullProvider,'matrices'),
        (run.kv,'correct'),(run.kv,'recover'),(run.sd,'recover_half'),(run.sd,'recover_closed'),
        (run,'inputs'),(run.screen,'configuration'),(run.completion,'full_control')):
        monkeypatch.setattr(obj,name,forbidden)
    before = {p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    assert run.compute() == dict(missing_only=True,new_roots=0,matrix_calls=0,form_recoveries=0)
    assert before == {p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
