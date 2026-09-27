"""D17: routing and fixed saved controls, no repeated spectral search."""
from dataclasses import replace
import json
import numpy as np
import pytest
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as solver
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_solver_architecture as run


@pytest.fixture(scope='module')
def properties():
    return kv.section()[1]


@pytest.fixture(scope='module')
def saved():
    if not (run.OUTPUT/'diagnostics.json').exists():
        pytest.skip('local D17 control artifacts unavailable')
    return run.screen.read_json(run.OUTPUT/'diagnostics.json')


@pytest.fixture(scope='module', autouse=True)
def counters():
    from unittest.mock import patch
    original=sd.Calls; ledger=[]
    def factory(**kwargs):
        obj=original(**kwargs); ledger.append(obj); return obj
    with patch.object(sd,'Calls',factory):
        yield
    keys=('B','B_z','full_B','half_B','expm','frechet','analytic_transfer','recoveries','corrections')
    print('\nD17_TEST_CALLS',json.dumps({k:sum(getattr(v,k) for v in ledger) for k in keys}),
          'BUILD_EQUIVALENTS',sum(v.cost() for v in ledger))


@pytest.mark.parametrize('field',['L','A','D','m'])
def test_exact_routing_never_uses_tolerance(properties,field):
    arm=kv.Arm.reduced('EB',properties)
    cfg=solver.Config((arm,arm),0.,1.,.001,mu=0.)
    assert solver.route(cfg)=='SYMMETRY_REDUCED'
    assert solver.route(cfg,'full')=='FULL_TWO_ARM'
    other=replace(arm,**{field:np.nextafter(getattr(arm,field),np.inf)})
    unequal=replace(cfg,arms=(arm,other))
    assert solver.route(unequal)=='FULL_TWO_ARM'
    with pytest.raises(ValueError,match='structurally identical'):
        solver.route(unequal,'reduced')
    assert solver.route(replace(cfg,mu=1e-30))=='FULL_TWO_ARM'
    assert solver.route(replace(cfg,mu=.01))=='FULL_TWO_ARM'
    # D20 extends the same API to RLB, but a mixed-theory pair is still invalid.
    with pytest.raises(ValueError,match='same EB or RLB theory'):
        solver.Config((arm,replace(arm,model='RLB')),0.,1.,.001)


@pytest.mark.parametrize('case_id',['A_d001','C_d001'])
def test_authoritative_transfer_and_derivative(properties,case_id,monkeypatch):
    case=next(c for c in run.inputs() if c['case_id']==case_id)
    cfg=run.configuration(properties,case); f=solver.FullProvider(cfg,sd.Calls())
    z=case['reference']; arm=cfg.arms[0]; units=arm.scale()
    original=solver.expm_frechet; flags=[]
    def frechet(*args,**kwargs):
        flags.append(kwargs.get('compute_expm')); return original(*args,**kwargs)
    monkeypatch.setattr(solver,'expm_frechet',frechet)
    T,Tz=f.transfer(z,arm,derivative=True)
    np.testing.assert_array_equal(T,f.transfer(z,arm))
    analytic=sd.closed_transfer(z/kv.T_REF,arm)
    scale=units[None,:]/units[:,None]
    assert np.linalg.norm((T-analytic)*scale)/np.linalg.norm(analytic*scale)<kv.CRITERIA['transfer_rtol']
    # Finite difference is only an independent derivative test, never the solver derivative.
    h=1e-4
    fd=(f.transfer(z+h,arm)-f.transfer(z-h,arm))/(2*h)
    assert np.linalg.norm((Tz-fd)*scale)/np.linalg.norm(Tz*scale)<kv.CRITERIA['derivative_rtol']
    assert flags==[False]


@pytest.mark.parametrize('case_id',run.CASE_IDS)
def test_saved_control_gates_and_physical_reactions(saved,properties,case_id):
    case=next(c for c in run.inputs() if c['case_id']==case_id)
    cfg=run.configuration(properties,case)
    with np.load(run.OUTPUT/'control_shapes.npz',allow_pickle=False) as archive:
        for path in case['paths']:
            key=case_id+'_'+path; point=saved['points'][key]; row=point['row']
            y=archive[key+'__states']; a=archive[key+'__a']; r=archive[key+'__reactions']
            assert row['solver_path']==solver.route(cfg,path)
            assert point['root_agreement']['accepted'] and not point['root_failures']
            assert not point['form_failures']
            assert point['root_equation_status']==point['form_recovery_status']=='PASS'
            if cfg.identical:
                assert point['diagnostics']['symmetry_class']==case['eta']
            assert row['physical_residual']<=kv.CRITERIA['physical_residual']
            assert row['energy_residual']<=kv.CRITERIA['energy_residual']
            assert row['MAC']>=kv.CRITERIA['MAC']
            assert point['correction']['steps']==0  # saved roots passed substitution
            np.testing.assert_allclose(r,a*np.tile([kv.F_REF,kv.F_REF,kv.M_REF],2),rtol=1e-14)
            np.testing.assert_allclose(y[:,0,3:].ravel(),r,rtol=1e-14,atol=1e-15)
            np.testing.assert_array_equal(y[:,0,:3],np.zeros((2,3)))
            _,weights=kv.quadrature()
            v=kv.mass_vector(y,cfg.arms,weights)
            assert np.vdot(v,v).real==pytest.approx(1.,abs=1e-12)
            if row['solver_path']=='SYMMETRY_REDUCED':
                np.testing.assert_array_equal(y[1],case['eta']*y[0]@sd.F)
                assert row['reduced_physical_residual']<=kv.CRITERIA['physical_residual']
            if key=='C_d001_full':
                assert not row['accepted'] and point['failures']==['POSSIBLE_MULTIPLICITY']
                assert row['next_sigma_ratio']<kv.CRITERIA['simple_sigma_separation']
                assert row['rank_status']=='QUALIFIED'
            else:
                assert row['accepted']


def test_inactive_reuse_has_no_complex_newton(saved,properties,monkeypatch):
    case=next(c for c in run.inputs() if c['case_id']=='K12_INACTIVE')
    # Exercise the routing branch without repeating reconstruction/diagnostics.
    def forbidden(*args,**kwargs):pytest.fail('inactive branch must not call Newton')
    monkeypatch.setattr(kv,'correct',forbidden)
    class ReachedRecovery(Exception):pass
    def stop(half,z,a):
        assert z==case['elastic_z'] and half.eta==-1
        raise ReachedRecovery
    monkeypatch.setattr(sd,'recover_closed',stop)
    with pytest.raises(ReachedRecovery):
        solver.solve_mode(run.configuration(properties,case),case['reference'],eta=-1,elastic_z=case['elastic_z'])
    point=saved['points']['K12_INACTIVE_auto']
    assert point['complex_newton_calls']==0 and point['row']['z_re']==0
    assert point['activity_status']=='EXACT_INACTIVE_BY_SYMMETRY'


def test_asymmetric_not_projected_or_symmetry_rejected(saved):
    p=saved['points']['K11_ASYMMETRIC_auto']
    assert not p['symmetry_reduced'] and not p['symmetry_gate_applicable'] and p['eta'] is None
    with np.load(run.OUTPUT/'control_shapes.npz') as archive:
        y=archive['K11_ASYMMETRIC_auto__states']
        assert np.linalg.norm(y[1]-y[0]@sd.F)/np.linalg.norm(y)>1e-4
    assert p['accepted']


@pytest.mark.parametrize('path',['auto','full'])
def test_fixed_equilibration_unscales_reaction_coordinates(properties,path):
    case=next(c for c in run.inputs() if c['case_id']=='A_d001')
    f=solver.FullProvider(run.configuration(properties,case),sd.Calls())
    matrix=sd.AnalyticHalfProvider(f,1) if path=='auto' else f
    frozen=sd.FrozenBalanced(matrix,case['reference'])
    B,_=matrix.matrices(case['reference']); balanced,_=frozen.matrices(case['reference'])
    b=np.arange(1,B.shape[1]+1)*(1+2j)
    np.testing.assert_allclose(balanced@b,B@frozen.reactions(b)/frozen.rows,rtol=1e-13)


def test_fixed_scope_and_preserved_sources(saved,properties):
    cases=run.inputs()
    assert len(cases)==7 and len(saved['points'])==12
    assert saved['criteria']==kv.CRITERIA
    assert saved['A_d005']==saved['C_d005']=='NOT_RUN'
    assert saved['new_d']==saved['new_beta']==saved['RLB_roots']==saved['new_physical_study']==0
    for case in cases:
        if case['case_id'].startswith(('A_','C_')):
            with pytest.raises(ValueError):run.configuration(properties,dict(case,d=.005))
        with pytest.raises(ValueError):run.configuration(properties,dict(case,beta_deg=89.))
        if case['case_id'].startswith('B_'):
            with pytest.raises(ValueError):run.configuration(properties,dict(case,d=.002))
    assert all(run.screen.sha(run.ROOT/p)==h for p,h in saved['protected_sources'].items())
    assert run.classify_stage(saved['points'])==(saved['reduced_status'],saved['full_status'])


def test_missing_only_zero_calls_bytes_unchanged(saved,monkeypatch):
    def forbidden(*args,**kwargs):pytest.fail('completed controls must be reused')
    for obj,name in ((solver,'solve_mode'),(solver.FullProvider,'matrices'),(sd.AnalyticHalfProvider,'matrices'),
                     (kv,'correct'),(kv,'recover'),(sd,'recover_closed')):
        monkeypatch.setattr(obj,name,forbidden)
    before={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    assert run.compute()==dict(missing_only=True,new_control_rows=0,matrix_calls=0,root_calls=0,form_recoveries=0)
    assert before=={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
