"""Targeted algebra/scaling checks; no repetition of the 12-case pilot."""
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.linalg import block_diag
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib.reddy_symmetric_coupled_beams import joint_matrix_closed_form


ARM = eb.EBArm(A=0.01, D=0.2*0.05**3/12, m=0.01, L=1)


@pytest.mark.parametrize("beta", [0.0, np.pi/6, -0.8, np.pi])
@pytest.mark.parametrize("stiffness", [0.0, 0.1, 1e4])
def test_six_rows_and_rank(beta, stiffness):
    joint = eb.Joint("SPRING", stiffness)
    matrix = eb.joint_matrix(beta, joint)
    y = np.array([.2, -.3, .4, 2., -3., 4., -.5, .6, -.7, -5., 6., -7.])
    u1, w1, p1, n1, q1, m1, u2, w2, p2, n2, q2, m2 = y
    c, s = np.cos(beta), np.sin(beta)
    expected = [u1+c*u2+s*w2, w1-s*u2+c*w2, m1+stiffness*(p1-p2),
                n1-c*n2-s*q2, q1+s*n2-c*q2, m1+m2]
    assert matrix.shape == (6, 12)
    assert_allclose(matrix@y, expected, rtol=1e-12, atol=1e-12)
    assert_allclose(matrix[[0, 1, 3, 4, 5]],
                    joint_matrix_closed_form(beta)[[0, 1, 3, 4, 5]], atol=1e-14)
    assert_allclose(np.linalg.det(matrix[:, [0, 1, 3, 4, 5, 11]]), 1, atol=1e-12)
    assert np.linalg.matrix_rank(matrix, tol=1e-12) == 6


@pytest.mark.parametrize("k", [0.0, 0.3, 8.0])
def test_energy_and_hinge(k):
    beta, p1, p2, dp1, dp2 = .6, .4, -.7, -.2, .9
    c, s = np.cos(beta), np.sin(beta)
    du2, dw2, n2, q2 = .3, -.5, 2., -4.
    du1, dw1 = -c*du2-s*dw2, s*du2-c*dw2
    n1, q1 = c*n2+s*q2, -s*n2+c*q2
    m1, m2 = -k*(p1-p2), k*(p1-p2)
    d_energy = k*(p1-p2)*(dp1-dp2)
    work = n1*du1+q1*dw1+m1*dp1+n2*du2+q2*dw2+m2*dp2
    assert_allclose(work+d_energy, 0, atol=1e-12)
    hessian = k*np.array([[1., -1.], [-1., 1.]])
    assert_allclose(hessian@np.ones(2), 0, atol=1e-12)
    assert np.min(np.linalg.eigvalsh(hessian)) >= -1e-12
    if k == 0:
        assert m1 == m2 == 0
        # A hinge still rejects incompatible translations.
        y = np.zeros(12)
        y[0] = 1
        assert (eb.joint_matrix(beta, eb.Joint("SPRING", 0))@y)[0] == 1


def test_exact_rigid_and_state_equations():
    assert eb.STATE_ORDER == ("u", "w", "psi", "N", "Q", "M")
    assert_allclose(eb.joint_matrix(.3, eb.Joint("RIGID")), joint_matrix_closed_form(.3), atol=1e-14)
    omega, y = .7, np.arange(1., 7.)
    expected = [y[3]/ARM.A, -y[2], y[5]/ARM.D, -ARM.m*omega**2*y[0],
                -ARM.m*omega**2*y[1], y[4]]
    assert_allclose(eb.state_matrix(omega, ARM)@y, expected)
    # Exact integration of the static EB equations, in dimensionless units.
    # A second dimensional expm can introduce noise into analytically zero entries.
    length, axial, bending = ARM.L, ARM.A, ARM.D
    static = np.eye(6)
    static[0, 3] = length/axial
    static[1, [2, 4, 5]] = [-length, -length**3/(6*bending), -length**2/(2*bending)]
    static[2, [4, 5]] = [length**2/(2*bending), length/bending]
    static[5, 4] = length
    scale = eb.state_scale(ARM)
    assert_allclose(eb.transfer_matrix(0, ARM)*scale[None, :]/scale[:, None],
                    static*scale[None, :]/scale[:, None], rtol=1e-12, atol=1e-12)
    assert_allclose(eb.clamp_to_joint_map(0, ARM), eb.transfer_matrix(0, ARM)[:, 3:])


@pytest.mark.parametrize("joint", [eb.Joint("RIGID"), eb.Joint("SPRING", 0),
                                   eb.Joint("SPRING", 10000*ARM.D/ARM.L)])
def test_physical_scaling_and_axial_separation(joint):
    other = eb.EBArm(ARM.A*1.1, ARM.D*.8, ARM.m*1.2, .7)
    assembly = eb.boundary_assembly(.2, ARM, other, .3, joint, ARM)
    raw = eb.joint_matrix(.3, joint) @ block_diag(
        eb.clamp_to_joint_map(.2, ARM), eb.clamp_to_joint_map(.2, other))
    assert_allclose(assembly.physical, raw, atol=1e-10, rtol=1e-12)
    assert_allclose(assembly.dimensionless, assembly.row_factors[:, None]*raw*
                    assembly.reaction_scales[None, :], atol=1e-10, rtol=1e-12)
    assert np.all(assembly.row_factors > 0)
    straight = eb.boundary_assembly(.2, ARM, ARM, 0, joint, ARM).physical
    assert_allclose(straight[np.ix_([0, 3], [1, 2, 4, 5])], 0, atol=1e-12)
    assert_allclose(straight[np.ix_([1, 2, 4, 5], [0, 3])], 0, atol=1e-12)


def test_endpoint_recovery_at_one_independent_rigid_root():
    # cos(z)*cosh(z)=1, full length 2L; avoids the coupled determinant.
    z = brentq(lambda value: np.cos(value)-1/np.cosh(value), 4.5, 5.0, xtol=1e-14)
    omega = z**2/(2*ARM.L)**2*np.sqrt(ARM.D/ARM.m)
    joint = eb.Joint("RIGID")
    assembly = eb.boundary_assembly(omega, ARM, ARM, 0, joint, ARM)
    diagnostic = eb.endpoint_diagnostics(assembly, 0, joint, ARM)
    assert diagnostic["sigma_ratio"] <= 1e-9
    assert diagnostic["nullity"] == 1
    for vector in diagnostic["vectors"]:
        assert max(abs(x) for x in vector["normalized_physical_residuals"]) <= 1e-9
        reactions = np.array(vector["physical_clamp_reactions"])
        states = block_diag(eb.clamp_to_joint_map(omega, ARM), eb.clamp_to_joint_map(omega, ARM))@reactions
        assert_allclose(states, vector["endpoint_states"], atol=1e-10)
    factors = np.array([2., .5, 3., 1/3, 5., .2])
    changed = eb.BoundaryAssembly(assembly.physical, factors[:, None]*assembly.dimensionless,
                                 assembly.endpoint_map, assembly.reaction_scales, assembly.row_units,
                                 factors*assembly.row_factors)
    rescaled = eb.endpoint_diagnostics(changed, 0, joint, ARM)
    assert rescaled["nullity"] == 1
    assert rescaled["sigma_ratio"] <= 1e-9
    assert max(abs(x) for x in rescaled["vectors"][0]["normalized_physical_residuals"]) <= 1e-9


def test_close_distinct_candidates_are_not_merged():
    from scripts.analysis.laminated_beams.pilot_inplane_rotational_spring_eb import consolidate
    def candidate(value, lo, hi):
        return SimpleNamespace(accepted=True, omega_bar=value, interval_left_bar=lo, interval_right_bar=hi,
                               diagnostics=SimpleNamespace(detected_nullity=1, root_gate_nullity=1))
    candidates = [candidate(10., 9.9, 10.000000000001), candidate(10.+1e-10, 10.00000000001, 10.1)]
    events, unresolved = consolidate(candidates)
    assert len(events) == len(unresolved) == 2


@pytest.mark.parametrize("doublet", [False, True])
def test_local_reconciliation_requires_one_separated_null_direction(doublet):
    from scripts.analysis.laminated_beams import pilot_inplane_rotational_spring_eb as pilot
    def provider(omega):
        value = omega*pilot.FREQUENCY_SCALE
        return np.diag([1e-8*(value-10.), 1e-8*(value-10.-1e-10) if doublet else 1., 1., 1., 1., 1.])
    candidates = []
    for value in (10., 10.+1e-10):
        diagnostic = pilot.roots.boundary_matrix_diagnostics(value, provider, pilot.FREQUENCY_SCALE)
        candidates.append(pilot.roots.RootCandidate(
            "test", "diagonal", "local", value, ("determinant_bracket",),
            9.99, 10.01, True, diagnostic, True, ""))
    resolved, evidence = pilot.reconcile_local_detections(candidates, provider)
    if doublet:
        assert len(resolved) == 2 and not evidence
    else:
        assert len(resolved) == 1 and len(evidence) == 1
        assert resolved[0].omega_bar == 10.


@pytest.mark.parametrize("name", ["A", "D", "m", "L"])
@pytest.mark.parametrize("bad", [0, -1, np.nan, np.inf])
def test_arm_validation(name, bad):
    values = dict(A=1., D=1., m=1., L=1.)
    values[name] = bad
    with pytest.raises(ValueError):
        eb.EBArm(**values)


@pytest.mark.parametrize("mode,k", [("SPRING", None), ("SPRING", -1), ("SPRING", np.nan),
                                   ("SPRING", np.inf), ("RIGID", 0), ("RIGID", np.inf), ("", None)])
def test_joint_validation(mode, k):
    with pytest.raises(ValueError):
        eb.Joint(mode, k)


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_finite_angle_and_frequency(bad):
    with pytest.raises(ValueError):
        eb.joint_matrix(bad, eb.Joint("RIGID"))
    with pytest.raises(ValueError):
        eb.state_matrix(bad, ARM)


def test_negative_frequency():
    with pytest.raises(ValueError):
        eb.state_matrix(-1, ARM)
