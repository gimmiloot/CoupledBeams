"""Targeted complex arithmetic/assembly tests; no repeated spectral pilot."""
from dataclasses import replace
import numpy as np
import pytest
from scipy.linalg import expm
from unittest.mock import patch
import json

from scripts.lib import inplane_kelvin_voigt as kv
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import inplane_rotational_spring_rlb as rlb
from scripts.lib import inplane_spring_modes as elastic

SCALAR_EVALUATIONS = 0


@pytest.fixture(scope='module',autouse=True)
def test_only_counters():
    original, ledger = kv.Calls, []
    def factory(**kwargs):
        value=original(**kwargs)
        ledger.append(value)
        return value
    with patch.object(kv,'Calls',factory):
        yield
    keys=[key for key in original().snapshot() if key!='limit']
    print('\nTEST_ONLY_PROVIDER_CALLS',json.dumps({key:sum(getattr(c,key) for c in ledger) for key in keys}),
          'SCALAR_PENCIL_EVALUATIONS',SCALAR_EVALUATIONS)


@pytest.fixture(scope='module')
def properties():
    return kv.section()[1]


@pytest.mark.parametrize('model', ['EB', 'RLB'])
def test_state_sign_and_complex_coefficients(properties, model):
    arm = kv.Arm.reduced(model, properties)
    p = .17+.31j
    H = kv.state_matrix(p, arm)
    assert H[3, 0] == properties.m*p*p
    assert H[4, 1] == properties.m*p*p
    assert H[5, 2] == arm.J*p*p
    assert H[1, 4] == arm.invS
    old = elastic.Arm(model, properties, 1.).matrix(.31)
    np.testing.assert_allclose(kv.state_matrix(.31j, arm), old, atol=1e-12, rtol=1e-12)
    assert H[3, 0].imag != 0


def test_complex_joint_against_independent_scalar_conditions():
    rng = np.random.default_rng(173)
    p, beta, k, c = -.3+.8j, .32, .43, .17
    u,w,f,n,q,m,U,W,F,N,Q,M = rng.normal(size=12)+1j*rng.normal(size=12)
    expected = [u+np.cos(beta)*U+np.sin(beta)*W,
        w-np.sin(beta)*U+np.cos(beta)*W, m+(k+c*p)*(f-F),
        n-np.cos(beta)*N-np.sin(beta)*Q, q+np.sin(beta)*N-np.cos(beta)*Q, m+M]
    state = [u,w,f,n,q,m,U,W,F,N,Q,M]
    np.testing.assert_allclose(kv.joint_matrix(p,beta,k,c)@state, expected)
    # K_J=0 is regular and retains all six constraints.
    assert np.linalg.matrix_rank(kv.joint_matrix(-k/c,beta,k,c)) == 6
    old = eb.joint_matrix(beta, eb.Joint('SPRING', k))
    np.testing.assert_array_equal(kv.joint_matrix(p,beta,k,0), old)


def test_svd_uses_conjugate_right_vector():
    B = np.array([[1, 1j], [2, 2j]])
    assert np.linalg.norm(B@kv.right_null(B)) < 1e-14
    assert np.linalg.norm(B@np.linalg.svd(B)[2][-1]) > 1


@pytest.mark.parametrize('model', ['EB', 'RLB'])
def test_fixed_scaling_derivative_and_conjugacy(properties, model):
    arm = kv.Arm.reduced(model, properties)
    provider = kv.Provider((arm, arm), .23, 1., .006)
    z = -.08+12.1j
    B, Bz = provider.matrices(z, derivative=True)
    h = 1e-5
    fd = (provider.matrices(z+h)[0]-provider.matrices(z-h)[0])/(2*h)
    np.testing.assert_allclose(Bz, fd, atol=1e-8, rtol=1e-7)
    np.testing.assert_allclose(provider.matrices(z.conjugate())[0], B.conj(), atol=1e-12)
    joint_p = kv.joint_matrix(z/kv.T_REF,.23,provider.k,provider.c,derivative=True)
    assert joint_p[2,2]/kv.T_REF == pytest.approx(provider.d*kv.M_REF, rel=1e-12, abs=0)
    assert provider.c == provider.d*kv.M_REF*kv.T_REF
    physical = kv.joint_matrix(z/kv.T_REF,.23,provider.k,provider.c)
    from scipy.linalg import block_diag
    physical = physical@block_diag(*(provider.transfer(z,a)[:,3:] for a in (arm,arm)))
    np.testing.assert_allclose(B, physical*provider.reaction_scales[None,:]/provider.row_units[:,None])
    assert provider.calls.B_z == 1


def test_limit_follows_rlb_coefficients(properties):
    r = replace(kv.Arm.reduced('RLB', properties), invS=0., J=0.)
    e = kv.Arm.reduced('EB', properties)
    z = -.2+7j
    np.testing.assert_array_equal(kv.state_matrix(z,r), kv.state_matrix(z,e))
    np.testing.assert_array_equal(kv.Provider((r,r),.3,1,.1).matrices(z)[0],
                                  kv.Provider((e,e),.3,1,.1).matrices(z)[0])


@pytest.mark.parametrize('model', ['EB', 'RLB'])
def test_complex_reactions_lengths_mass_and_direct_transfer(properties, model):
    a = kv.Arm.reduced(model, properties, .9)
    b = kv.Arm.reduced(model, properties, 1.1)
    provider = kv.Provider((a,b), .2, 1, .001)
    reactions_hat = np.array([1+2j,.1j,.2,3-.4j,.7j,-.3j])
    z = -.3+7j
    shape = kv.recover(provider,z,reactions_hat,direct_check=True)
    assert max(shape['direct_errors']) < 1e-9
    assert abs(np.vdot(shape['vector'],shape['vector'])-1) < 1e-12
    np.testing.assert_allclose(shape['states'][:,0,3:].ravel(),shape['reactions'])
    np.testing.assert_allclose(shape['reactions']/shape['reactions'][0],
        reactions_hat*provider.reaction_scales/(reactions_hat[0]*provider.reaction_scales[0]))
    for i, arm in enumerate((a,b)):
        direct = expm(kv.state_matrix(z/kv.T_REF,arm)*arm.L)@np.r_[np.zeros(3),shape['reactions'].reshape(2,3)[i]]
        np.testing.assert_allclose(shape['states'][i,-1],direct,rtol=1e-9,atol=1e-10)
    _, weights = kv.quadrature()
    states = np.zeros((2,129,6),complex)
    states[0,:,0], states[1,:,1], states[:,:,2] = 2j,3,4j
    mass = np.linalg.norm(kv.mass_vector(states,(a,b),weights))**2
    assert mass == pytest.approx(4*a.m*a.L+9*b.m*b.L+16*(a.J*a.L+b.J*b.L))
    # Arbitrary reactions must fail physical conditions; moments are not filled from the joint law.
    assert max(kv.diagnose(provider,z,shape)['physical_residuals']) > 1e-5


def test_complex_mac_amplitude_and_phase():
    v = np.array([1,2j,.3+4j])
    assert kv.mac_matrix([v],[(-2+5j)*v])[0,0] == pytest.approx(1)


@pytest.mark.parametrize('inertia,c,k', [(2.,0.,3.),(2.,.2,3.),(1.,3.,1.)])
def test_augmented_corrector_scalar_convention(inertia,c,k):
    expected = (-c+np.sqrt(complex(c*c-4*inertia*k)))/(2*inertia)
    def pencil(p, *, derivative=False):
        global SCALAR_EVALUATIONS
        SCALAR_EVALUATIONS += 1
        return np.array([[inertia*p*p+c*p+k]]), np.array([[2*inertia*p+c]])
    result = kv.correct(pencil, expected+.01+.02j, np.array([1+0j]))
    assert result['status'] == 'CONVERGED'
    assert result['z'] == pytest.approx(expected,abs=1e-10)
    assert result['z'].real <= 1e-10


@pytest.mark.parametrize('value', [-1, np.nan, np.inf, 1j])
def test_invalid_passive_parameters(properties,value):
    arm=kv.Arm.reduced('EB',properties)
    with pytest.raises((ValueError,TypeError)):
        kv.Provider((arm,arm),.2,1,value)


def test_unstable_trial_points_and_cost_cap(properties):
    arm=kv.Arm.reduced('RLB',properties)
    calls=kv.Calls(limit=2)
    provider=kv.Provider((arm,arm),.2,1,.1,calls)
    assert np.isfinite(provider.matrices(1+2j,derivative=True)[0]).all()
    with pytest.raises(RuntimeError,match='COST_LIMIT'):
        provider.matrices(1+2j)
    assert calls.B+calls.B_z == 2


def test_independent_protected_root_not_forced():
    def pencil(z, *, derivative=False):
        return np.diag([z*z+.1*z+1,z*z+4]), np.diag([2*z+.1,2*z])
    original=2j+1e-4
    result=kv.correct(pencil,original,np.array([0.,1.]))
    assert result['steps']>0
    assert abs(result['z']-2j)<1e-10
    assert original.real==1e-4


def test_json_and_complex_shape_serialization(tmp_path):
    from scripts.analysis.laminated_beams import pilot_inplane_kelvin_voigt as run
    text=json.dumps(run.clean(dict(z=-.1+2j,array=np.array([1+3j]))),allow_nan=False)
    assert json.loads(text)['z']==dict(real=-.1,imag=2.)
    with pytest.raises(ValueError):
        run.clean(dict(invalid=np.nan))
    data=np.array([[1+2j,3-4j]])
    np.savez_compressed(tmp_path/'shapes.npz',mode_a=data,mode_b=2j*data)
    with np.load(tmp_path/'shapes.npz') as saved:
        np.testing.assert_array_equal(saved['mode_a'],data)
        np.testing.assert_array_equal(saved['mode_b'],2j*data)


def test_old_output_directory_cannot_be_overwritten():
    from scripts.analysis.laminated_beams import pilot_inplane_kelvin_voigt as run
    with pytest.raises(ValueError,match='read-only'):
        run.Run(run.SOURCE)


def test_inactive_gate_keeps_actual_sign():
    diag=dict(null_residual=0,sigma_ratio=0,r_E=0,symmetry_defect=0,
        conjugate_residual=0,clamp_residual=0,omega_identity_residual=0,
        physical_residuals=[0]*6,next_sigma_ratio=.1,a_energy=0.,Delta_psi=0j)
    z=2e-7+5j
    before=z
    assert 'INACTIVE_GATE' in kv.failures(diag,z,'INACTIVE',5,1.)
    assert z==before and z.real>0
