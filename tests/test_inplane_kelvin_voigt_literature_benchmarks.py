"""Targeted source/sign/scaling tests; no repeated literature spectrum."""
import numpy as np
import pytest
from scripts.lib import inplane_kelvin_voigt_literature_benchmarks as lit
from scripts.lib import inplane_kelvin_voigt as kv


def test_failla_time_and_dimensional_mapping():
    L,m,D=2.,3.,5.
    time=lit.reference_time(L,m,D)
    omega=21.9+.263j
    z=lit.failla_to_project(omega)
    assert z==-0.263+21.9j
    assert time==pytest.approx(4*np.sqrt(3/5))
    assert np.exp(1j*(omega/time)*.12)==pytest.approx(np.exp(z*.12/time))
    cu,cr=.7,.8
    assert cu*L/np.sqrt(m*D)==pytest.approx(cu*L**3/(D*time))
    assert cr/(L*np.sqrt(m*D))==pytest.approx(cr*L/(D*time))


@pytest.mark.parametrize('ku,gu,kr,gr',[(3.,.2,None,0.),(0.,0.,4.,.1),(3.,.2,4.,.1)])
def test_failla_interfaces(ku,gu,kr,gr):
    z=-.2+3j;minus=np.array([.1+.2j,.4,-.3j,2.])
    plus=lit.failla_jump(z,ku,gu,kr,gr)@minus
    expected=minus.copy()
    if kr is not None:expected[1]-=minus[2]/(kr+gr*z)
    expected[3]+=(ku+gu*z)*minus[0]
    np.testing.assert_allclose(plus,expected)
    L,R=lit.failla_interface(z,ku,gu,kr,gr)
    np.testing.assert_allclose(L@minus+R@plus,0,atol=1e-14)


def test_inactive_bare_mode_is_not_assigned_damping():
    beta=4*np.pi;z=1j*beta**2
    for xi in (.25,.5,.75):
        state=np.array([np.sin(beta*xi),beta*np.cos(beta*xi),
                        beta**2*np.sin(beta*xi),beta**3*np.cos(beta*xi)],complex)
        np.testing.assert_allclose(lit.failla_jump(z)@state,state,atol=1e-11)
    # A positive damping coefficient remains present in the operator.
    Lz,Rz=lit.failla_interface(z,derivative=True)
    assert Lz[1,1]==-.1 and Rz[1,1]==.1 and Lz[3,0]==-.1


def test_hong_geometry_shape_factor_and_independent_moduli():
    h=lit.Hong()
    assert h.A==pytest.approx(.000625)
    assert h.I==pytest.approx(3.255208333333334e-8)
    assert h.K==pytest.approx(13/15.3)
    assert h.G==80e9 and h.G!=h.E/(2*(1+h.nu))


@pytest.mark.parametrize('s',[-.3+20j,1.+13j,-2.+800j])
def test_hong_source_matrix_and_project_mapping(s):
    h=lit.Hong()
    expected=np.array([[0,1,-1/(h.K*h.A*h.G),0],
       [0,0,0,1/(h.E*h.I)],[-h.rho*h.A*s*s,0,0,0],
       [0,h.rho*h.I*s*s,1,0]],complex)
    np.testing.assert_allclose(lit.hong_state(s),expected,rtol=kv.CRITERIA['H_rtol'],atol=0.)
    arm=kv.Arm('RLB',h.E*h.A,h.D,h.m,h.L,1/h.S,h.J)
    H=kv.state_matrix(s,arm)[np.ix_([1,2,4,5],[1,2,4,5])]
    S=np.diag([1,-1,-1,-1])
    np.testing.assert_array_equal(lit.hong_state(s),S@H@S)
    assert lit.hong_state(s)[2,0].imag!=0


def test_failla_project_mapping():
    z=-.1+12j
    arm=kv.Arm('EB',1,1,1,1,0,0)
    H=kv.state_matrix(z,arm)[np.ix_([1,2,4,5],[1,2,4,5])]
    P=np.array([[1,0,0,0],[0,-1,0,0],[0,0,0,1],[0,0,1,0]])
    np.testing.assert_array_equal(lit.failla_state(z),P@H@P.T)


@pytest.mark.parametrize('case',['failla','hong_hh','hong_ff','hong_damped'])
def test_boundary_scalar_equations_derivative_and_conjugacy(case):
    beam=lit.Beam(case);z=-.15+24j
    B,Bz=beam.matrices(z,derivative=True)
    step=1e-5
    fd=(beam.matrices(z+step)[0]-beam.matrices(z-step)[0])/(2*step)
    np.testing.assert_allclose(Bz,fd,rtol=1e-7,atol=1e-8)
    np.testing.assert_allclose(beam.matrices(z.conjugate())[0],B.conj(),rtol=1e-12,atol=1e-12)
    rng=np.random.default_rng(316)
    a=rng.normal(size=16)+1j*rng.normal(size=16)
    y,norm_a=beam.recover(z,a)
    scale=np.max(abs(y/beam.scale))
    np.testing.assert_allclose(beam.physical_conditions(z,y),B@norm_a/scale,atol=1e-12,rtol=1e-10)


def test_hong_boundaries_and_support_stiffness():
    s=-3+12j;h=lit.Hong();z=s*h.time
    hh=lit.Beam('hong_hh');ff=lit.Beam('hong_ff');support=lit.Beam('hong_damped')
    L,R=hh.boundary(z)
    np.testing.assert_array_equal(L,np.eye(4)[[0,3]])
    np.testing.assert_array_equal(R,L)
    np.testing.assert_array_equal(ff.boundary(z)[0],np.eye(4)[[2,3]])
    L,R=support.boundary(z)
    K=(2e6+20*s)*h.L**3/h.D
    assert L[0,0]==pytest.approx(K/1000)
    assert R[0,0]==pytest.approx(-K/1000)
    assert L[1,3]==R[1,3]==1


@pytest.mark.parametrize('printed,half',[('21.9037',5e-5),('238.023',5e-4),
    ('3.55778498e+002',5e-7),('-6.6651e-002',5e-7),('0.0120',5e-5)])
def test_published_rounding(printed,half):
    v=float(printed)
    assert lit.rounding(printed,v)['rounding_tolerance']==pytest.approx(half)
    assert lit.rounding(printed,v+.4*half)['rounding_pass']
    assert not lit.rounding(printed,v+3*half)['rounding_pass']
    assert not lit.rounding(printed,np.nan)['rounding_pass']


def test_literal_missing_zero_and_complex_dtype():
    assert lit.rounding(None,1e-11)['rounding_pass']
    assert not lit.rounding(None,1e-4)['rounding_pass']
    assert np.iscomplexobj(lit.Beam('failla').matrices(2j)[0])
