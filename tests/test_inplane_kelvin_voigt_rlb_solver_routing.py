"""D20 focused production regressions; saved controls, no repeated root solve."""
from dataclasses import replace
import json
from unittest.mock import patch
import numpy as np
import pytest
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_solver as solver
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.analysis.laminated_beams import check_inplane_kelvin_voigt_rlb_solver_architecture as run


@pytest.fixture(scope='module')
def properties():
    return kv.section()[1]


@pytest.fixture(scope='module')
def cases():
    return {c['case_id']:c for c in run.inputs()}


@pytest.fixture(scope='module')
def saved():
    if not (run.OUTPUT/'diagnostics.json').exists():
        pytest.skip('local D20 artifacts unavailable; do not recreate old source data')
    return run.screen.read_json(run.OUTPUT/'diagnostics.json')


@pytest.fixture(scope='module',autouse=True)
def counters():
    original=sd.Calls; ledger=[]
    def factory(**kwargs):
        obj=original(**kwargs); ledger.append(obj); return obj
    with patch.object(sd,'Calls',factory):
        yield
    keys=('B','B_z','full_B','half_B','expm','frechet','analytic_transfer','recoveries','corrections')
    print('\nD20_TEST_CALLS',json.dumps({k:sum(getattr(v,k) for v in ledger) for k in keys}),
          'BUILD_EQUIVALENTS',sum(v.cost() for v in ledger))


@pytest.mark.parametrize('field',['L','A','D','m','invS','J'])
def test_exact_rlb_routing(properties,cases,field):
    arm=kv.Arm.reduced('RLB',properties)
    cfg=run.configuration(properties,cases['K12_RLB_ACTIVE_d1'])
    assert cfg.model=='RLB' and solver.route(cfg)=='SYMMETRY_REDUCED'
    assert solver.route(cfg,'full')=='FULL_TWO_ARM'
    other=replace(arm,**{field:np.nextafter(getattr(arm,field),np.inf)})
    unequal=replace(cfg,arms=(arm,other))
    assert solver.route(unequal)=='FULL_TWO_ARM'
    with pytest.raises(ValueError,match='structurally identical'): solver.route(unequal,'reduced')
    assert solver.route(replace(cfg,mu=1e-30))=='FULL_TWO_ARM'


def test_direct_transfer_frechet_only_and_no_eb_provider(properties,cases,monkeypatch):
    case=cases['K12_RLB_ACTIVE_d1']; cfg=run.configuration(properties,case)
    f=solver.FullProvider(cfg,sd.Calls()); arm=cfg.arms[0]; z=case['reference']
    def forbidden(*a,**k): pytest.fail('RLB must not call EB analytic or historical transfer')
    monkeypatch.setattr(kv.Provider,'transfer',forbidden)
    monkeypatch.setattr(sd,'closed_transfer',forbidden)
    monkeypatch.setattr(sd,'recover_closed',forbidden)
    seen=[]; original=solver.expm_frechet
    def frechet(*a,**k):
        seen.append(k.get('compute_expm')); return original(*a,**k)
    monkeypatch.setattr(solver,'expm_frechet',frechet)
    T,Tz=f.transfer(z,arm,derivative=True)
    np.testing.assert_array_equal(T,f.transfer(z,arm))
    assert seen==[False] and np.isfinite(Tz).all()
    assert type(solver.reduced_provider(f,1)) is sd.HalfProvider


@pytest.mark.parametrize('eta',[1,-1])
def test_half_conditions_signs_and_complex_joint(properties,cases,eta):
    c=cases['K12_RLB_ACTIVE_d1']; cfg=run.configuration(properties,c)
    f=solver.FullProvider(cfg,sd.Calls()); p=c['reference']/kv.T_REF
    ch,sh=np.cos(f.beta/2),np.sin(f.beta/2)
    expected=np.zeros((3,6),complex)
    expected[0,[0,1]]=([ch,-sh] if eta==1 else [sh,ch])
    expected[1,[3,4]]=([sh,ch] if eta==1 else [ch,-sh])
    expected[2,5]=1
    if eta==1: expected[2,2]=2*(f.k+f.c*p)
    np.testing.assert_array_equal(sd.conditions(p,f.beta,f.k,f.c,eta),expected)
    if eta==-1:
        np.testing.assert_array_equal(expected,sd.conditions(p,f.beta,0.,0.,eta))
        assert not np.any(sd.conditions(p,f.beta,f.k,f.c,eta,derivative=True))


@pytest.mark.parametrize('index',[0,1,2])
def test_matrix_gate_before_roots_and_derivative_checks(saved,index):
    check=saved['matrix_validation']['checks'][index]
    assert check['accepted'] and check['H_p_exact'] and check['inactive_independent']
    assert check['FH_HF_norm']==0 and check['row_transform_rank']==6
    for key in ('block_relative','derivative_relative','off_block_relative','joint_reflection_relative'):
        assert check[key]<=kv.CRITERIA['H_rtol']
    assert check['B_z_finite_difference']['accepted'] and check['conjugacy']['accepted']
    assert check['T_reflection']['accepted'] and check['T_authoritative']['accepted']
    assert saved['matrix_validation']['calls']['corrections']==0
    assert saved['matrix_validation']['calls']['recoveries']==0


@pytest.mark.parametrize('key',[f'K12_RLB_ACTIVE_d{i}' for i in (1,2,3)])
def test_active_saved_reproduction_full_control_energy_and_form(saved,cases,properties,key):
    point=saved['points'][key]; row=point['row']; cfg=run.configuration(properties,cases[key])
    assert row['accepted'] and row['solver_path']=='SYMMETRY_REDUCED'
    assert point['root_agreement']['accepted'] and row['reference_MAC']>=kv.CRITERIA['MAC']
    assert row['reduced_physical_residual']<=kv.CRITERIA['physical_residual']
    assert row['lifted_full_physical_residual']<=kv.CRITERIA['physical_residual']
    assert row['root_residual']<=kv.CRITERIA['null_residual'] and row['sigma_ratio']<=kv.CRITERIA['sigma_ratio']
    assert row['energy_residual']<=kv.CRITERIA['energy_residual']
    assert row['conjugate_residual']<=kv.CRITERIA['null_residual']
    assert abs(row['alpha']-row['alpha_energy'])<=kv.CRITERIA['a_atol']/kv.T_REF+kv.CRITERIA['a_rtol']*abs(row['alpha'])
    full=saved['points'][key+'_full_control']
    assert full['row']['accepted'] and full['row']['root_origin']=='EVALUATION_ONLY_AT_REDUCED_ROOT'
    assert full['row']['complex_newton_calls']==0 and full['full_reduced_form']['MAC']>=kv.CRITERIA['MAC']
    assert full['row']['p_re']==row['p_re'] and full['row']['p_im']==row['p_im']
    with np.load(run.OUTPUT/'control_shapes.npz',allow_pickle=False) as archive:
        y=archive[key+'__states']; a=archive[key+'__a']; reactions=archive[key+'__reactions']
    np.testing.assert_array_equal(y[1],y[0]@sd.F)
    np.testing.assert_allclose(reactions,a*np.tile([kv.F_REF,kv.F_REF,kv.M_REF],2),rtol=1e-14)
    np.testing.assert_allclose(y[:,0,3:].ravel(),reactions,rtol=1e-14,atol=1e-15)
    np.testing.assert_array_equal(y[:,0,:3],np.zeros((2,3)))
    _,w=kv.quadrature()
    m=sum(arm.L*np.dot(w,arm.m*(abs(v[:,0])**2+abs(v[:,1])**2)+arm.J*abs(v[:,2])**2)
          for arm,v in zip(cfg.arms,y))
    assert m==pytest.approx(1.,abs=1e-12)
    rotary=sum(arm.J*arm.L*np.dot(w,abs(v[:,2])**2) for arm,v in zip(cfg.arms,y))
    assert rotary>0 and abs(m-rotary-1)>1e-6  # catches EB-only mass
    shear=sum(arm.invS*arm.L*np.dot(w,abs(v[:,4])**2) for arm,v in zip(cfg.arms,y))
    assert shear>0
    delta=y[0,-1,2]-y[1,-1,2]
    stiffness=sum(arm.L*np.dot(w,abs(v[:,3])**2/arm.A+abs(v[:,5])**2/arm.D)
                  for arm,v in zip(cfg.arms,y))+shear+kv.M_REF*abs(delta)**2
    assert stiffness==pytest.approx(row['K_phi'],rel=1e-13)
    p=complex(row['p_re'],row['p_im'])
    energy=p*p*m+p*row['C_phi']+stiffness
    assert abs(energy)/(abs(p)**2*m+abs(p)*row['C_phi']+stiffness)<=kv.CRITERIA['energy_residual']
    # RLB section rotation is independent: w'+psi = Q/S, generally nonzero.
    derivative=y[0]@kv.state_matrix(p,cfg.arms[0]).T
    np.testing.assert_allclose(derivative[:,1]+y[0,:,2],cfg.arms[0].invS*y[0,:,4],atol=1e-12)
    assert np.max(abs(derivative[:,1]+y[0,:,2]))>1e-5


def test_inactive_bypasses_newton_and_eb_recovery(properties,cases,monkeypatch):
    case=cases['K12_RLB_INACTIVE_d1']; cfg=run.configuration(properties,case)
    def forbidden(*a,**k): pytest.fail('inactive RLB must not call Newton or EB recovery')
    monkeypatch.setattr(kv,'correct',forbidden); monkeypatch.setattr(sd,'recover_closed',forbidden)
    class ReachedRecovery(Exception): pass
    def stop(half,z,a):
        assert half.eta==-1 and type(half) is sd.HalfProvider
        assert z==case['elastic_z'] and z.real==0
        raise ReachedRecovery
    monkeypatch.setattr(sd,'recover_half',stop)
    with pytest.raises(ReachedRecovery):
        solver.solve_mode(cfg,case['reference'],eta=-1,elastic_z=case['elastic_z'])


@pytest.mark.parametrize('i',[1,2,3])
def test_inactive_saved_exact_reuse(saved,i):
    key=f'K12_RLB_INACTIVE_d{i}'; point=saved['points'][key]; r=point['row']
    assert r['accepted'] and r['inactive_by_symmetry']
    assert point['activity_status']=='EXACT_INACTIVE_BY_SYMMETRY'
    assert r['complex_newton_calls']==r['newton_updates']==r['z_re']==0
    with np.load(run.OUTPUT/'control_shapes.npz') as a:
        y=a[key+'__states']
    np.testing.assert_array_equal(y[1],-y[0]@sd.F)
    assert r['Delta_psi_re']==r['Delta_psi_im']==r['C_phi']==0


def test_exact_limit_matrices_root_and_form(saved,properties):
    arm=replace(kv.Arm.reduced('RLB',properties),invS=0.,J=0.)
    eb=kv.Arm.reduced('EB',properties)
    assert arm.model=='RLB' and arm.invS==arm.J==0
    for derivative in (False,True):
        np.testing.assert_array_equal(kv.state_matrix(-.03+.1j,arm,derivative=derivative),
                                      kv.state_matrix(-.03+.1j,eb,derivative=derivative))
    assert all(v['accepted'] for v in saved['matrix_validation']['exact_limit']+saved['limit_root_form'])
    assert {v['quantity'] for v in saved['matrix_validation']['exact_limit']}=={
        'H','H_p','T','T_p','B_plus','B_plus_p','B_minus','B_minus_p','B_full','B_full_p'}
    point=saved['points']['EXACT_EB_LIMIT']
    assert point['row']['accepted'] and point['row']['reference_MAC']>=kv.CRITERIA['MAC']


def test_asymmetric_full_path_has_no_projection(saved,cases,properties):
    case=cases['K11_RLB_ASYMMETRIC']; cfg=run.configuration(properties,case)
    point=saved['points'][case['case_id']]
    assert [a.L for a in cfg.arms]==[.99,1.01] and cfg.d_theta==0
    assert point['accepted'] and point['eta'] is None and not point['symmetry_gate_applicable']
    assert point['solver_path']=='FULL_TWO_ARM' and not point['symmetry_reduced']
    with np.load(run.OUTPUT/'control_shapes.npz') as a: y=a[case['case_id']+'__states']
    assert np.linalg.norm(y[1]-y[0]@sd.F)/np.linalg.norm(y)>1e-4
    assert point['row']['reference_MAC']>=kv.CRITERIA['MAC']


def test_asymmetric_positive_d_matrix_only(properties,cases):
    # Allowed matrix probe: existing K11 lengths and a K12 d, never a root solve.
    active=cases['K12_RLB_ACTIVE_d1']
    cfg=replace(run.configuration(properties,cases['K11_RLB_ASYMMETRIC']),d_theta=active['d'])
    f=solver.FullProvider(cfg,sd.Calls());z=active['reference']
    B,Bz=f.matrices(z,derivative=True)
    h=1e-4
    fd=(f.matrices(z+h)[0]-f.matrices(z-h)[0])/(2*h)
    error=np.linalg.norm(Bz-fd)/np.linalg.norm(Bz)
    assert error<=kv.CRITERIA['derivative_rtol']
    one,two=(f.transfer(z,arm)[:,3:] for arm in cfg.arms)
    from scipy.linalg import block_diag
    expected=kv.joint_matrix(z/kv.T_REF,f.beta,f.k,f.c)@block_diag(one,two)
    expected*=f.reaction_scales[None,:]/f.row_units[:,None]
    np.testing.assert_array_equal(B,expected)
    assert np.linalg.norm(one-two)/np.linalg.norm(one)>1e-3
    assert solver.route(cfg)=='FULL_TWO_ARM'
    print('\nD20_UNEQUAL_MATRIX_DERIVATIVE_RELATIVE',float(error))


def test_fixed_scope_preservation_and_budget(saved,cases,properties):
    assert tuple(cases)==run.CASE_IDS and len(saved['points'])==11
    assert saved['status']=='RLB_KV_PRODUCTION_PASS'
    assert saved['rlb_active_corrector_calls']==3 and saved['full_evaluation_controls']==3
    for name in ('inactive_complex_newton','asymmetric_positive_d_roots','new_beta',
                 'new_d_outside_K12','new_physical_parameter_study','FEM','high_precision'):
        assert saved[name]==0
    assert saved['calls']['total_build_equivalents']<=run.BUDGET
    assert saved['criteria']==kv.CRITERIA
    assert all(run.screen.sha(run.ROOT/p)==h for p,h in saved['protected_sources'].items())
    for case in cases.values():
        with pytest.raises(ValueError): run.configuration(properties,dict(case,d=.001))
        with pytest.raises(ValueError): run.configuration(properties,dict(case,beta_deg=45.))


def test_manifest_scope_matches_actual_controls(saved):
    config=saved['configuration']; controls=saved['controls']
    assert config==run.configuration_record(config,controls)
    assert config['model']=='RLB' and config['beta_deg']==[5.]
    assert {c['case_id'] for c in config['cases']}==set(run.CASE_IDS)
    for item,control in zip(config['cases'],controls):
        assert (item['mu'],item['d_theta'])==(control['mu'],control['d'])
        assert item['L1']==1-control['mu'] and item['L2']==1+control['mu']


def test_missing_only_has_no_matrix_root_or_shape_calls(saved,monkeypatch):
    def forbidden(*a,**k): pytest.fail('completed regression must perform zero scientific work')
    for obj,name in ((solver,'solve_mode'),(solver.FullProvider,'matrices'),(kv,'correct'),
                     (kv,'recover'),(sd,'recover_half'),(run,'matrix_checks'),(run,'inputs')):
        monkeypatch.setattr(obj,name,forbidden)
    before={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    result=run.compute()
    assert result==dict(missing_only=True,new_rows=0,root_calls=0,matrix_calls=0,form_recoveries=0,status=saved['status'])
    assert before=={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
