"""D16: algebra, scaling and saved diagnostics; no repeated A/C root solve."""
import json
from dataclasses import replace
from unittest.mock import patch
import numpy as np
import pytest
from scipy.linalg import expm
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_kelvin_voigt_symmetry_diagnostics as sd
from scripts.analysis.laminated_beams import diagnose_inplane_kelvin_voigt_ac as run


@pytest.fixture(scope='module', autouse=True)
def counters():
    original = sd.Calls; ledger = []
    keys = ('B','B_z','full_B','half_B','expm','frechet','shape_expm','analytic_transfer','recoveries','corrections')
    analytic = sd.closed_transfer; direct = expm
    extra = dict(analytic=0, direct=0)
    def analytic_counted(*args, **kwargs):
        extra['analytic'] += 1
        return analytic(*args, **kwargs)
    def direct_counted(*args, **kwargs):
        extra['direct'] += 1
        return direct(*args, **kwargs)
    def factory(**kwargs):
        value = original(**kwargs)
        ledger.append((value,{k:getattr(value,k) for k in keys},value.cost()))
        return value
    with patch.object(sd, 'Calls', factory), patch.object(sd, 'closed_transfer', analytic_counted), patch(__name__+'.expm', direct_counted):
        yield
    totals={k:sum(getattr(v,k)-before[k] for v,before,cost in ledger) for k in keys}
    print('\nD16_TEST_CALLS', json.dumps(totals),
          'BUILD_EQUIVALENTS', sum(v.cost()-cost for v,before,cost in ledger),
          'STANDALONE_ANALYTIC_TRANSFERS', extra['analytic']-totals['analytic_transfer'],
          'STANDALONE_DIRECT_EXPM', extra['direct'])


@pytest.fixture(scope='module')
def properties():
    return kv.section()[1]


@pytest.fixture(scope='module')
def saved():
    if not (run.OUTPUT/'diagnostics.json').exists():
        pytest.skip('local D16 diagnostic artifacts unavailable')
    return run.screen.read_json(run.OUTPUT/'diagnostics.json')


@pytest.mark.parametrize('eta', [1, -1])
@pytest.mark.parametrize('beta', [0., 75.])
def test_full_joint_reduction_from_independent_scalar_equations(eta, beta):
    p = -.002+1.1j; k = kv.M_REF; c = .001*k*kv.T_REF
    rng = np.random.default_rng(173)
    y = rng.normal(size=6)+1j*rng.normal(size=6)
    ends = np.r_[y, eta*sd.F@y]
    b = np.deg2rad(beta); ch, sh = np.cos(b/2), np.sin(b/2)
    u,w,psi,N,Q,M = y
    if eta == 1:
        expected = np.array([ch*u-sh*w, sh*N+ch*Q, M+2*(k+c*p)*psi])
        full_expected = np.array([2*ch*expected[0],-2*sh*expected[0],expected[2],
                                  2*sh*expected[1],2*ch*expected[1],0])
    else:
        expected = np.array([sh*u+ch*w, ch*N-sh*Q, M])
        full_expected = np.array([2*sh*expected[0],2*ch*expected[0],M,
                                  2*ch*expected[1],-2*sh*expected[1],2*M])
    np.testing.assert_allclose(sd.conditions(p,b,k,c,eta)@y, expected, atol=1e-14)
    np.testing.assert_allclose(kv.scalar_conditions(ends,p,b,k,c), full_expected, atol=1e-14)
    np.testing.assert_allclose(kv.joint_matrix(p,b,k,c)@ends, full_expected, atol=1e-14)
    mirror = np.r_[sd.F@ends[6:],sd.F@ends[:6]]
    np.testing.assert_allclose(mirror, eta*ends)


@pytest.mark.parametrize('beta', [0., 75.])
def test_complex_block_factorization_and_conjugation(properties, beta):
    sid = 'A_STRONG' if beta == 0 else 'C_WEAK_ACTIVE'
    f = run.provider(properties, sid, beta, .001, sd.Calls())
    assert sd.algebra_check(f, -.1+32j)['accepted']
    for eta in (1,-1):
        half = sd.HalfProvider(f,eta)
        B,_ = half.matrices(-.1+32j)
        conjugate,_ = half.matrices(-.1-32j)
        np.testing.assert_allclose(conjugate, B.conj(), atol=1e-12,rtol=1e-12)
    H = kv.state_matrix(-.1+1j,f.arms[0])
    np.testing.assert_array_equal(sd.F@H,H@sd.F)


def test_minus_independent_plus_full_complex_stiffness(properties):
    arm = kv.Arm.reduced('EB',properties); calls=sd.Calls()
    f = sd.FullProvider((arm,arm),.5,1,.001,calls)
    g = sd.FullProvider((arm,arm),.5,0,0,calls)
    z = -.2+80j
    np.testing.assert_array_equal(sd.HalfProvider(f,-1).matrices(z)[0],sd.HalfProvider(g,-1).matrices(z)[0])
    C = sd.conditions(z/kv.T_REF,f.beta,f.k,f.c,1)
    assert C[2,2] == pytest.approx(2*(f.k+f.c*z/kv.T_REF), rel=1e-15, abs=0)
    assert C[2,2].real != 2*f.k
    np.testing.assert_array_equal(sd.conditions(z/kv.T_REF,f.beta,f.k,f.c,-1,derivative=True),np.zeros((3,6)))


def test_reaction_lift_projection_and_no_scale_mixing():
    a=np.array([1+2j,3,4j,5,6j,7],complex)
    plus,minus=sd.project(a)
    np.testing.assert_allclose(plus+minus,a)
    assert abs(np.vdot(plus,minus))<1e-13
    np.testing.assert_array_equal(sd.R,np.diag([1,-1,-1]))
    for eta,v in ((1,plus),(-1,minus)):
        np.testing.assert_allclose(v[3:],eta*sd.R@v[:3])
        np.testing.assert_allclose(sd.lift(v[:3],eta),v)


@pytest.mark.parametrize('sid',['A_STRONG','C_WEAK_ACTIVE'])
def test_saved_lifted_solution_full_conditions_and_scaling(saved, properties, sid):
    q=saved['points'][sid]; r=q['row']
    f=run.provider(properties,sid,r['beta_deg'],.001,sd.Calls())
    z=complex(r['z_half_re'],r['z_half_im'])
    with np.load(run.OUTPUT/'ac_reduced_shapes.npz') as archive:
        y=archive[sid+'__states']; a=archive[sid+'__a']
    np.testing.assert_array_equal(y[1],y[0]@sd.F)
    np.testing.assert_allclose(a[3:],sd.R@a[:3])
    half=sd.HalfProvider(f,1)
    p=sd.physical_details(half,z,{'states':y})
    assert max(p['full_normalized'])<=kv.CRITERIA['physical_residual']
    assert max(p['full_normalized'][:2])<=kv.CRITERIA['compatibility']
    p2=sd.physical_details(half,z,{'states':y*(-3+2j)})
    np.testing.assert_allclose(p2['full_normalized'],p['full_normalized'],atol=2e-15,rtol=1e-3)
    _,w=kv.quadrature()
    assert np.linalg.norm(kv.mass_vector(y,f.arms,w))**2 == pytest.approx(1.,abs=1e-13)
    assert q['closed_attempt']['full_K12_gate_failures'] == ([] if sid=='A_STRONG' else ['POSSIBLE_MULTIPLICITY'])


@pytest.mark.parametrize('sid',['A_STRONG','C_WEAK_ACTIVE'])
def test_conditional_closed_transfer_derivative_and_real_limit(saved, properties, sid):
    r=saved['points'][sid]['row']; arm=kv.Arm.reduced('EB',properties)
    z=complex(r['z_half_re'],r['z_half_im']); p=z/kv.T_REF; units=arm.scale()
    T,Tp=sd.closed_transfer(p,arm,derivative=True)
    scaled=kv.state_matrix(p,arm)*units[None,:]/units[:,None]
    expected=expm(scaled*arm.L)
    np.testing.assert_allclose(T*units[None,:]/units[:,None],expected,atol=1e-10,rtol=1e-11)
    h=1e-6
    fd=(sd.closed_transfer(p+h,arm)-sd.closed_transfer(p-h,arm))/(2*h)
    assert np.linalg.norm((Tp-fd)*units[None,:]/units[:,None])/np.linalg.norm(Tp*units[None,:]/units[:,None])<1e-7
    # Independent pre-existing real EB formulas, including their physical signs.
    from scripts.lib import inplane_rotational_spring_eb as eb
    e=eb.EBArm(arm.A,arm.D,arm.m,arm.L)
    reactions=np.array([.2,.3,-.4]); xi=np.array([0.,.5,1.])
    old=sd.modes.arm_states(p.imag,e,reactions,xi)
    # Only the endpoint requires another analytic transfer here.
    value=sd.closed_transfer(1j*p.imag,arm)@np.r_[np.zeros(3),reactions]
    np.testing.assert_allclose(value,old[-1],rtol=1e-11,atol=1e-8)
    np.testing.assert_allclose(sd.closed_transfer(p.conjugate(),arm),T.conj(),atol=1e-10,rtol=1e-12)


def test_fixed_equilibration_derivative_no_adaptive_scaling(properties):
    f=run.provider(properties,'C_WEAK_ACTIVE',75.,.001,sd.Calls())
    half=sd.ClosedHalfProvider(f,1,triggered=True)
    b=sd.FrozenBalanced(half,-.006+100j)
    rows,cols=b.rows.copy(),b.cols.copy()
    z=-.007+100.01j; h=1e-5
    B,Bz=b.matrices(z,derivative=True)
    fd=(b.matrices(z+h)[0]-b.matrices(z-h)[0])/(2*h)
    assert np.linalg.norm(Bz-fd)/np.linalg.norm(Bz)<1e-7
    np.testing.assert_array_equal(rows,b.rows); np.testing.assert_array_equal(cols,b.cols)
    r=np.array([1j,2,3+4j])
    np.testing.assert_allclose(B@r,half.matrices(z)[0]@b.reactions(r)/rows)


def test_scope_selection_no_new_d_no_B_no_rlb(saved, properties):
    old=run.screen.read_json(run.K16/'diagnostics.json')
    selected=run.select_inputs(old)
    assert [(sid,s['beta_deg'],s['elastic_sorted_mode']) for sid,p,s in selected]==list(sd.SELECTION)
    assert set(saved['points'])=={'A_STRONG','C_WEAK_ACTIVE'}
    for d in (.005,0.,.002):
        with pytest.raises(ValueError,match='only A/C'):
            run.provider(properties,'A_STRONG',0.,d,sd.Calls())
    with pytest.raises(ValueError):
        run.provider(properties,'B_INTERMEDIATE',45.,.001,sd.Calls())
    arm=kv.Arm.reduced('EB',properties)
    with pytest.raises(ValueError,match='identical EB'):
        sd.FullProvider((arm,replace(arm,L=1.01)),0.,1.,.001,sd.Calls())
    with pytest.raises(ValueError,match='identical EB'):
        sd.FullProvider((replace(arm,model='RLB'),)*2,0.,1.,.001,sd.Calls())
    with pytest.raises(ValueError,match='not triggered'):
        sd.ClosedHalfProvider(run.provider(properties,'A_STRONG',0.,.001,sd.Calls()),1,triggered=False)


def test_conditional_beta_trigger_and_fixed_plan(saved):
    assert sd.LOCAL_BETAS==(70.,72.5,75.,77.5,80.)
    flags=dict(opposite_block_near=False,same_block_suspect=False,second_sigma_unexplained=False)
    assert not sd.beta_trigger(**flags) and sd.local_beta_plan(False)==()
    for key in flags:
        assert sd.beta_trigger(**dict(flags,**{key:True}))
    assert sd.local_beta_plan(True)==sd.LOCAL_BETAS
    assert saved['local_beta_roots']==0 and not saved['local_beta_ran']
    assert not saved['points']['C_WEAK_ACTIVE']['rank_interpretation']['local_beta_trigger']
    assert saved['new_physical_targets']==saved['new_d_values']==saved['B_resolves']==saved['RLB_roots']==0


def test_source_statuses_and_attempts_preserved(saved):
    old=run.screen.read_json(run.K16/'diagnostics.json')
    assert old['stage_status']=='PARTIAL_NUMERICAL_QUALIFICATIONS'
    for sid,q in saved['points'].items():
        assert q['original_row']==old['points'][sid+'_d1']['row']
        assert q['original_attempts']==old['points'][sid+'_d1']['attempts']
        assert q['attempts']<=2
        assert q['correction']['steps']<=20 and q['closed_attempt']['correction']['steps']<=20
        assert q['analytic_transfer_trigger'] and q['analytic_transfer_ran']
        assert q['row']['diagnostic_status']=='FULL_TRANSFER_RECOVERY_CONDITIONING'
    assert saved['source_hashes']==run.protected_hashes()
    assert saved['criteria']==kv.CRITERIA
    assert saved['calls']['total_build_equivalents']<1000


def test_missing_only_zero_calls_unchanged_bytes(saved,monkeypatch):
    def forbidden(*args,**kwargs):
        pytest.fail('missing-only must not recalculate completed diagnostics')
    for obj,name in ((kv,'correct'),(kv,'recover'),(kv,'diagnose'),(sd,'closed_transfer'),
                     (sd.FullProvider,'matrices'),(sd.HalfProvider,'matrices')):
        monkeypatch.setattr(obj,name,forbidden)
    before={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}
    r=run.compute()
    assert r['missing_only'] and all(r[k]==0 for k in ('half_resolves','full_B','half_B','B_z','expm','recoveries'))
    assert before=={p.name:run.screen.sha(p) for p in run.OUTPUT.iterdir() if p.is_file()}


def test_failed_algebra_checkpoint_cannot_resume_into_root_search(saved,monkeypatch,tmp_path):
    import copy
    import shutil
    data=copy.deepcopy(saved)
    data['finished']=False
    data['algebra'][0]['accepted']=False
    run.write_json(tmp_path/'diagnostics.json',data)
    shutil.copyfile(run.OUTPUT/'ac_reduced_shapes.npz',tmp_path/'ac_reduced_shapes.npz')
    monkeypatch.setattr(run,'OUTPUT',tmp_path)
    def forbidden(*args,**kwargs):
        pytest.fail('failed algebra gate must stop before any root audit')
    monkeypatch.setattr(run,'audit_point',forbidden)
    monkeypatch.setattr(kv,'correct',forbidden)
    with pytest.raises(RuntimeError,match='REDUCTION_ALGEBRA_GATE'):
        run.compute()
