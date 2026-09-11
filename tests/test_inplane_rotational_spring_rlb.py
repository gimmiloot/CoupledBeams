"""Small algebra/matrix checks; no spectral searches or pilot reruns."""
from dataclasses import replace

import numpy as np
import pytest
from scipy.linalg import block_diag, expm

from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import inplane_rotational_spring_rlb as rlb
from scripts.lib import reddy_symmetric_laminated_beam as native

REF = eb.EBArm(.20*.05, .20*.05**3/12, .20*.05, 1.)
SPRING = eb.Joint("SPRING", REF.D/REF.L)


def test_constitutive_chain():
    section, p = rlb.benchmark_section()
    assert len(section.plies) == 4
    assert all(ply.angle_deg == 0 and ply.thickness == .05/4 for ply in section.plies)
    h, nu, g = .05, .3, 1/2.6
    q = np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1-nu)/2]])/(1-nu**2)
    np.testing.assert_allclose(section.A, q*h, rtol=1e-12, atol=0)
    np.testing.assert_allclose(section.D, q*h**3/12, rtol=1e-12, atol=0)
    np.testing.assert_allclose(section.shear, np.eye(2)*g*h, rtol=1e-12, atol=0)
    assert np.linalg.norm(section.B) <= 1e-12*np.linalg.norm(section.A)*h
    assert abs(section.I1) <= 1e-12*section.I0*h
    np.testing.assert_allclose([section.I0, section.I2], [h, h**3/12], rtol=1e-12)
    np.testing.assert_allclose([p.A, p.D, p.m, p.S, p.J],
                               [REF.A, REF.D, REF.m, 5/6*g*.20*h, .20*h**3/12], rtol=1e-12)


@pytest.mark.parametrize("epsilon", [0., .01, .1, 1.])
def test_coefficients_and_frozen_section(epsilon):
    _, p = rlb.benchmark_section()
    arm = rlb.LimitArm(p, 1., epsilon)
    omega = .37
    expected = np.zeros((6, 6))
    expected[0, 3], expected[1, 2], expected[2, 5] = 1/p.A, -1, 1/p.D
    expected[3, 0] = expected[4, 1] = -p.m*omega**2
    expected[5, 4] = 1
    expected[1, 4], expected[5, 2] = epsilon/p.S, -epsilon*p.J*omega**2
    np.testing.assert_array_equal(rlb.state_matrix(omega, arm), expected)
    assert arm.properties is p and (p.A, p.D, p.m) == (arm.properties.A, arm.properties.D, arm.properties.m)
    assert arm.invS == epsilon/p.S and arm.J == epsilon*p.J
    assert SPRING.k_theta == REF.D/REF.L


@pytest.mark.parametrize("beta", [0., np.pi/6])
@pytest.mark.parametrize("Omega", [2., 20., 80.])
def test_native_and_exact_limit(beta, Omega):
    _, p = rlb.benchmark_section()
    omega = Omega / np.sqrt(REF.m/REF.D)
    one, zero = rlb.LimitArm(p, 1., 1.), rlb.LimitArm(p, 1., 0.)
    np.testing.assert_array_equal(rlb.state_matrix(omega, one), native.combined_state_matrix(omega, p))
    z = eb.state_scale(REF)
    for left, right in [(rlb.transfer_matrix(omega, one), native.combined_transfer_matrix(omega, 1., p)),
                        (rlb.transfer_matrix(omega, zero), eb.transfer_matrix(omega, REF))]:
        left, right = left*z[None, :]/z[:, None], right*z[None, :]/z[:, None]
        assert np.linalg.norm(left-right) <= 1e-12 + 1e-9*np.linalg.norm(right)
    left = rlb.boundary_assembly(omega, zero, zero, beta, SPRING, REF)
    right = eb.boundary_assembly(omega, REF, REF, beta, SPRING, REF)
    assert np.linalg.norm(left.dimensionless-right.dimensionless) <= 1e-12+1e-9*np.linalg.norm(right.dimensionless)
    rigid = eb.boundary_assembly(omega, REF, REF, beta, eb.Joint("RIGID"), REF)
    assert not np.allclose(left.dimensionless[2], rigid.dimensionless[2])


def test_rlb_endpoints_scaling_and_no_eb_dispatch(monkeypatch):
    _, p = rlb.benchmark_section()
    one, two = rlb.LimitArm(p, .8, .1), rlb.LimitArm(replace(p, A=1.3*p.A, D=.7*p.D), 1.1, .01)
    def forbidden(*args, **kwargs):
        raise AssertionError("unexpected EB assembly")
    for name in ("state_matrix", "transfer_matrix", "boundary_assembly", "_scaled_transfer"):
        monkeypatch.setattr(eb, name, forbidden)
    for epsilon in (0., .1):
        one = replace(one, epsilon_limit=epsilon)
        assembly = rlb.boundary_assembly(.31, one, two, .4, SPRING, REF)
        physical_maps = block_diag(*(expm(rlb.state_matrix(.31, a)*a.L)[:, 3:] for a in (one, two)))
        np.testing.assert_allclose(assembly.endpoint_map, physical_maps*assembly.reaction_scales[None, :], rtol=1e-9, atol=1e-12)
        np.testing.assert_allclose(assembly.dimensionless, assembly.row_factors[:, None]*assembly.physical*assembly.reaction_scales[None, :])
        diagnostic = eb.endpoint_diagnostics(assembly, .4, SPRING, REF)
        for vector in diagnostic["vectors"]:
            states = np.array(vector["endpoint_states"])
            reactions = np.array(vector["physical_clamp_reactions"])
            np.testing.assert_allclose(physical_maps@reactions, states, rtol=1e-9, atol=1e-12)
            u,w,psi,n,q,m,U,W,P,N,Q,M = states
            c,s = np.cos(.4),np.sin(.4)
            expected = [u+c*U+s*W,w-s*U+c*W,m+SPRING.k_theta*(psi-P),n-c*N-s*Q,q+s*N-c*Q,m+M]
            np.testing.assert_allclose(vector["normalized_physical_residuals"], np.array(expected)/assembly.row_units)


@pytest.mark.parametrize("epsilon", [-1., 1.01, np.nan, np.inf])
def test_invalid_epsilon(epsilon):
    with pytest.raises(ValueError):
        rlb.LimitArm(rlb.benchmark_section()[1], 1., epsilon)


@pytest.mark.parametrize("length", [0., -1., np.nan, np.inf])
def test_invalid_length(length):
    with pytest.raises(ValueError):
        rlb.LimitArm(rlb.benchmark_section()[1], length, 1.)


@pytest.mark.parametrize("frequency", [-1., np.nan, np.inf])
def test_invalid_frequency(frequency):
    with pytest.raises(ValueError):
        rlb.state_matrix(frequency, rlb.LimitArm(rlb.benchmark_section()[1], 1., 0.))


def test_invalid_angle():
    arm = rlb.LimitArm(rlb.benchmark_section()[1], 1., 1.)
    with pytest.raises(ValueError):
        rlb.boundary_assembly(.1, arm, arm, np.nan, SPRING, REF)


def test_provider_counts_builds_and_enforces_shared_limit(monkeypatch):
    from types import SimpleNamespace
    from scripts.analysis.laminated_beams import check_inplane_rotational_spring_rlb_eb_limit as pilot
    made = []
    def assembly(*args):
        made.append(args[0])
        return SimpleNamespace(dimensionless=np.eye(6))
    monkeypatch.setattr(rlb, "boundary_assembly", assembly)
    provider = pilot.Provider(rlb.benchmark_section()[1], pilot.cases()[0])
    provider.builds = pilot.LIMITS["max_builds_per_group"]-1
    first = provider(.31)
    assert provider(.31) is first
    with pytest.raises(pilot.CostLimit, match="COST_LIMIT"):
        provider(.32)
    assert made == [.31] and provider.builds == 6000


def test_rejected_guard_stops_search_without_tail_or_third_recovery(monkeypatch):
    from types import SimpleNamespace
    from scripts.analysis.laminated_beams import check_inplane_rotational_spring_rlb_eb_limit as pilot
    def candidate(value, accepted):
        return SimpleNamespace(omega_bar=value, accepted=accepted,
            rejection_reason="" if accepted else "NULLITY_UNRESOLVED_AT_1E-12",
            detection_sources=("synthetic_test",), interval_left_bar=value-.01, interval_right_bar=value+.01,
            diagnostics=SimpleNamespace(detected_nullity=1 if accepted else 0,
                root_gate_nullity=1, scaled_sigma_ratio=1e-15 if accepted else 1e-11))
    accepted = [candidate(x/4, True) for x in range(1,7)]
    guard = candidate(3., False)
    calls = []
    def scan(*args, **kwargs):
        calls.append((args[2:4],kwargs))
        return accepted+[guard]
    monkeypatch.setattr(pilot.workflow, "scan", scan)
    monkeypatch.setattr(pilot.workflow, "reconcile_local_detections", lambda pool,provider:(pool,[]))
    monkeypatch.setattr(pilot.workflow, "consolidate", lambda pool:(accepted,[]))
    budget = [2]
    result = pilot.solve_group(rlb.benchmark_section()[1], pilot.cases()[0], [], budget)
    assert result["status"] == "INCOMPLETE" and "e_max" not in result
    assert len(calls) == len(result["windows"]) == 1
    assert budget == [2] and result["boundary_builds"] == 0
