"""D18 scope and saved evidence; no repeat of any physical root or full control."""
import ast
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from scripts.analysis.laminated_beams import complete_inplane_kelvin_voigt_weak_damping as run


@pytest.fixture(scope='module')
def saved():
    if not (run.OUTPUT/'diagnostics.json').exists():pytest.skip('local completion data unavailable')
    return run.screen.read_json(run.OUTPUT/'diagnostics.json')


@pytest.fixture(scope='module')
def sources():return run.seeds_and_sources()


def test_exact_two_targets_reduced_eb_only(sources):
    seeds,_=sources;p=run.kv.section()[1]
    assert run.TARGETS==(('A_STRONG',0.,'sorted_05',.005),('C_WEAK_ACTIVE',75.,'sorted_05',.005))
    for seed in seeds:
        if seed['state_id'].startswith('B_'):
            with pytest.raises(ValueError):run.config_for(p,seed,.005)
        else:
            cfg=run.config_for(p,seed,.005)
            assert run.production.route(cfg,'reduced')=='SYMMETRY_REDUCED'
            assert cfg.mu==0 and all(a.model=='EB' and a.J==a.invS==0 for a in cfg.arms)
        for d in (.001,.002,.003,.004,.01):
            with pytest.raises(ValueError):run.config_for(p,seed,d)
        with pytest.raises(ValueError):run.config_for(p,dict(seed,beta_deg=4.),.005)


def test_confirmed_continuation_and_single_fallback(sources,saved):
    seeds,prior=sources
    for seed in (seeds[0],seeds[2]):
        sid=seed['state_id'];r=prior['points'][sid[0]+'_d001_auto']['row']
        z=complex(r['z_re'],r['z_im'])
        assert run.predictor(seed,prior,1)==z
        assert run.predictor(seed,prior,2)==z-seed['a_slope_pred']*.004
        with pytest.raises(ValueError):run.predictor(seed,prior,3)
        q=saved['points'][sid];assert q['source_row']==r and len(q['attempts'])==1
        a=q['attempts'][0]
        assert complex(a['predictor']['real'],a['predictor']['imag'])==z
        assert a['solver_path']=='reduced' and a['eta']==1 and a['d_theta']==.005


def test_total_reduced_budget_stops_without_new_algorithm():
    c=run.CompletionCalls(B=150,B_z=150,half_B=150,half_B_z=150)
    with pytest.raises(RuntimeError,match='REDUCED_COST_LIMIT'):c.matrix(False)
    assert c.B==150
    # A full evaluation does not spend the reduced budget.
    c.full_B+=1;c.matrix(False);assert c.B==151


def test_reuse_is_read_only_exact_four_rows(sources,saved,monkeypatch):
    def forbidden(*a,**k):pytest.fail('reuse must not evaluate matrices, forms or roots')
    for obj,name in ((run.production,'solve_mode'),(run.production.FullProvider,'matrices'),
                     (run.kv,'recover'),(run.kv,'correct')):monkeypatch.setattr(obj,name,forbidden)
    rows=run.reused_rows(*sources)
    assert rows==saved['reused_rows'] and len(rows)==4
    assert [(r['state'],r['d_theta']) for r in rows]==[('A',.001),('B',.001),('B',.005),('C',.001)]
    for r in rows:
        key=r['data_origin'].split(':')[1];old=sources[1]['points'][key]['row']
        assert (r['z_re'],r['z_im'])==(old['z_re'],old['z_im'])
        if r['state']!='B':assert 'K16 historical full recovery rejected' in r['notes'] and 'K17/K18' in r['notes']


@pytest.mark.parametrize('sid',['A_STRONG','C_WEAK_ACTIVE'])
def test_new_root_production_gates_energy_mass_and_lift(saved,sid):
    point=saved['points'][sid];r=point['row'];a=point['attempts'][-1];q=a['result']
    assert point['accepted'] and q['accepted'] and q['failures']==[]
    assert q['eta']==1 and q['solver_path']=='SYMMETRY_REDUCED'
    assert q['correction']['status']=='CONVERGED' and 0<q['correction']['steps']<=run.kv.CRITERIA['max_steps']
    assert q['diagnostics']['sigma_ratio']<=run.kv.CRITERIA['sigma_ratio']
    assert max(q['reduced_physical']['half_normalized'])<=run.kv.CRITERIA['physical_residual']
    assert max(q['diagnostics']['physical_residuals'])<=run.kv.CRITERIA['physical_residual']
    p=complex(r['p_re'],r['p_im']);M,K,C=r['M_phi'],r['K_phi'],r['C_phi']
    residual=abs(p*p*M+p*C+K)/(abs(p)**2*M+abs(p)*C+K)
    assert residual==pytest.approx(r['energy_residual'],rel=1e-6,abs=1e-15)
    assert residual<=run.kv.CRITERIA['energy_residual']
    assert r['alpha']==-p.real and r['alpha_energy']==pytest.approx(C/(2*M))
    assert abs((r['alpha']-r['alpha_energy'])*run.kv.T_REF)<=run.kv.CRITERIA['a_atol']+run.kv.CRITERIA['a_rtol']*r['a']
    with np.load(run.OUTPUT/'new_complex_shapes.npz') as archive:
        y=archive[sid+'_attempt1__states'];v=archive[sid+'_attempt1__vector']
        np.testing.assert_array_equal(y[1],y[0]@run.sd.F)
        assert np.vdot(v,v).real==pytest.approx(1.,abs=1e-12)
    assert r['MAC']>=run.kv.CRITERIA['MAC']
    full=a['full_control'];assert full['newton_calls']==0
    assert full['z']==q['z']
    assert full['failures']==([] if sid=='A_STRONG' else ['POSSIBLE_MULTIPLICITY'])


@pytest.mark.parametrize('state,d',[(s,d) for s in 'ABC' for d in (.001,.005)])
def test_six_state_observables_from_actual_complex_components(state,d):
    rows=run.screen.read_csv(run.OUTPUT/'combined_six_state_summary.csv');assert len(rows)==6
    row=next(r for r in rows if r['state']==state and float(r['d_theta'])==d)
    n=lambda key:float(row[key])
    assert row['root_status']=='ROOT_ACCEPTED'
    assert n('a_over_d')==pytest.approx(-n('z_re')/d)
    zeta=-n('z_re')/np.hypot(n('z_re'),n('z_im'))
    assert n('zeta_over_d')==pytest.approx(zeta/d)
    assert n('zeta')==pytest.approx(-n('p_re')/np.hypot(n('p_re'),n('p_im')))
    assert n('Omega_d')==n('z_im')  # not absolute eigenvalue
    shift=(n('z_im')-n('Omega_0'))/n('Omega_0')
    assert n('relative_frequency_shift')==pytest.approx(shift)
    assert n('relative_frequency_shift_over_d2')==pytest.approx(shift/d**2)
    for prefix in ('a','zeta'):
        assert n('relative_'+prefix+'_slope_error')==pytest.approx(n(prefix+'_over_d')/n(prefix+'_slope_pred')-1,abs=1e-15)


def test_full_control_has_no_newton_evaluation_only(monkeypatch):
    # Synthetic plumbing test: no second full evaluation at a scientific point.
    events=[]
    def forbidden(*a,**k):pytest.fail('no full Newton allowed')
    monkeypatch.setattr(run.production,'solve_mode',forbidden);monkeypatch.setattr(run.kv,'correct',forbidden)
    provider=SimpleNamespace(arms=('fake','fake'))
    monkeypatch.setattr(run.production,'FullProvider',lambda *a:provider)
    class Balanced:
        def __init__(self,*a):events.append('B_scaling')
        def matrices(self,z):events.append('B');return np.eye(2),None
        def reactions(self,b):return b
    monkeypatch.setattr(run.sd,'FrozenBalanced',Balanced)
    def recover(*args):
        events.append('recover');return {k:np.ones(2,complex) for k in run.old.ARRAYS}
    monkeypatch.setattr(run.kv,'recover',recover)
    monkeypatch.setattr(run.kv,'mass_vector',lambda *a:np.ones(2,complex))
    monkeypatch.setattr(run.kv,'diagnose',lambda *a:events.append('diagnose') or {})
    monkeypatch.setattr(run.kv,'failures',lambda *a:[])
    q=run.full_control(provider,1j,dict(states=np.ones(2),Omega0=1),run.CompletionCalls())
    assert events==['B_scaling','B','recover','diagnose'] and q['newton_calls']==0


def test_scope_counts_and_architecture_unchanged(saved):
    assert saved['new_principal_roots']==saved['new_roots_accepted']==2
    assert saved['full_control_evaluations']==2
    assert saved['calls']['half_B']+saved['calls']['half_B_z']<=300
    for key in ('auxiliary_roots','new_beta','RLB','solver_architecture_changes','high_precision','full_newton_calls'):assert saved[key]==0
    assert saved['criteria']==run.kv.CRITERIA
    assert all(run.screen.sha(run.ROOT/p)==h for p,h in saved['protected_sources'].items())
    # Dispatcher can only be called from the restricted reduced orchestration.
    tree=ast.parse(Path(run.__file__).read_text())
    calls=[n for n in ast.walk(tree) if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute) and n.func.attr=='solve_mode']
    assert len(calls)==1
    assert next(k.value.value for k in calls[0].keywords if k.arg=='solver_path')=='reduced'


def test_missing_only_zero_computation_and_unchanged_bytes(saved,monkeypatch):
    def forbidden(*a,**k):pytest.fail('finished work must not be repeated')
    for obj,name in ((run.production,'solve_mode'),(run,'full_control'),(run.kv,'correct'),(run.kv,'recover'),(run.production.FullProvider,'matrices')):
        monkeypatch.setattr(obj,name,forbidden)
    before={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    result=run.compute()
    assert result['missing_only'] and all(result[k]==0 for k in ('new_roots','B','B_z','form_recoveries','full_control_evaluations'))
    assert before=={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
