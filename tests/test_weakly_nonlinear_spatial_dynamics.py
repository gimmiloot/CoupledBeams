"""Addressed algebraic tests; no ODE, static Newton, or native FEM jobs."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin
from scripts.lib.weakly_nonlinear_spatial_dynamics import FIELDS, SpatialGalerkin

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def model():
    path = ROOT / "results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json"
    if not path.exists():
        pytest.skip("Saved frozen quartic action is not locally available")
    pol = json.loads(path.read_text(encoding="utf8"))["polynomials"]
    return SimpleNamespace(T4=rod.Polynomial.deserialize(pol["T4"]),
        V4=rod.Polynomial.deserialize(pol["V4"]),
        residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol["residuals_A"]),
        symbols={f: rod.Polynomial.symbol(f) for f in rod.SYMBOL_ORDER})


@pytest.fixture(scope="module")
def coefficients():
    return rod.RodCoefficients.rectangular(1., 1., .3, .2, .1, 1.759089824002232e-5)


@pytest.fixture
def disc(model, coefficients):
    return SpatialGalerkin(coefficients, 6, model=model)


def coordinates(disc, seed=41):
    rng = np.random.default_rng(seed)
    scales = np.array((2e-5, 5e-4, 3e-4, 3e-3, 3e-3, 3e-3, 3e-5))
    raw = np.concatenate([rng.normal(size=disc.n)*a/np.arange(1, disc.n+1)**2 for a in scales])
    velocity_raw = np.concatenate([rng.normal(size=disc.n)*a/np.arange(1, disc.n+1)**2 for a in scales])
    return disc.from_raw_coefficients(raw), disc.from_raw_coefficients(velocity_raw)


def embed(spatial, planar, value):
    result = np.zeros(spatial.ndof)
    for f in planar.fields:
        result[spatial.slices[f]] = value[planar.slices[f]]
    return result


def ids(spatial, planar):
    return np.concatenate([np.arange(spatial.slices[f].start, spatial.slices[f].stop) for f in planar.fields])


def test_contract_and_no_hidden_derivation(coefficients, model, monkeypatch):
    monkeypatch.setattr(rod, "derive_polynomials", lambda: pytest.fail("New derivation forbidden"))
    spatial = SpatialGalerkin(coefficients, 5, model=model)
    assert spatial.fields == ("u", "w", "v", "Phi", "psi", "theta", "c")
    assert spatial.ndof == 28 and spatial.nq == 11
    with pytest.raises(ValueError, match="serialized"):
        SpatialGalerkin(coefficients, 5)
    with pytest.raises(ValueError, match="reduced integration"):
        SpatialGalerkin(coefficients, 5, nq=10, model=model)


@pytest.mark.parametrize("whiten", [False, True])
def test_basis_endpoint_values_only_and_raw_roundtrip(coefficients, model, whiten):
    spatial = SpatialGalerkin(coefficients, 6, model=model, whiten=whiten)
    q, _ = coordinates(spatial)
    assert np.max(abs(spatial.reconstruct(q, [0., 1.]))) == 0.
    assert np.max(abs(spatial.reconstruct(q, [0., 1.], derivative=1))) > 0.
    np.testing.assert_allclose(spatial.from_raw_coefficients(spatial.raw_coefficients(q)), q, rtol=2e-14, atol=2e-20)
    values = spatial.reconstruct(q)
    np.testing.assert_allclose(spatial.project(values), q, rtol=5e-13, atol=2e-19)


def test_resting_and_variable_mass_kinetic_identity(disc):
    q, velocity = coordinates(disc)
    np.testing.assert_allclose(disc.mass_matrix(np.zeros(disc.ndof)), disc.M0, rtol=3e-13, atol=3e-13)
    matrix = disc.mass_matrix(q)
    np.testing.assert_allclose(matrix, matrix.T, rtol=3e-14, atol=3e-14)
    assert np.linalg.eigvalsh(matrix).min() > 0.
    assert np.max(abs(matrix-disc.M0)) > 1e-5
    assert disc.kinetic(q, velocity) == pytest.approx(.5*velocity@matrix@velocity, rel=2e-13)
    assert disc.diagnostics(q)["mass_positive"]


def test_gradient_hessian_directional_derivatives(disc):
    q, direction = coordinates(disc)
    expected = disc.potential(q, hessian=True)
    errors = []
    for step in (1e-2, 1e-3, 1e-4):
        plus = disc.potential(q+step*direction, hessian=True)
        minus = disc.potential(q-step*direction, hessian=True)
        energy_derivative = (plus["V"]-minus["V"])/(2*step)
        gradient_derivative = (plus["gradient"]-minus["gradient"])/(2*step)
        errors.append(np.linalg.norm(gradient_derivative-expected["hessian"]@direction)/np.linalg.norm(expected["hessian"]@direction))
        assert energy_derivative == pytest.approx(expected["gradient"]@direction, rel=2e-5, abs=1e-14)
    assert min(errors) < 2e-8
    np.testing.assert_allclose(expected["hessian"], expected["hessian"].T, rtol=3e-13, atol=3e-13)


def test_full_rhs_jacobian_independent_directional_check(disc):
    q, velocity = coordinates(disc)
    dq, dv = coordinates(disc, seed=97)
    state, direction = np.r_[q, velocity], np.r_[dq, dv]
    analytic = disc.jacobian(0., state)@direction
    errors = [np.linalg.norm((disc.rhs(0., state+h*direction)-disc.rhs(0., state-h*direction))/(2*h)-analytic)/np.linalg.norm(analytic)
              for h in (1e-2, 1e-3, 1e-4)]
    assert min(errors) < 2e-8


def test_variational_strong_weak_and_energy_identity(disc):
    q, velocity = coordinates(disc)
    acceleration = disc.acceleration(q, velocity)
    algebraic = disc.mass_matrix(q)@acceleration+disc.inertial_terms(q, velocity)+disc.potential(q)["gradient"]
    strong_projected = disc.weak_residual(q, velocity, acceleration)
    scale = np.linalg.norm(disc.potential(q)["gradient"])+np.linalg.norm(disc.inertial_terms(q, velocity))
    assert np.linalg.norm(strong_projected-algebraic)/scale < 2e-12
    power_scale = abs(velocity@disc.potential(q)["gradient"])+abs(velocity@disc.mass_matrix(q)@acceleration)
    assert abs(disc.energy_rate(q, velocity, acceleration))/power_scale < 2e-12


def test_planar_action_mass_inertia_rhs_and_jacobian(disc, model, coefficients):
    planar = PlanarGalerkin(coefficients, disc.p, nq=disc.nq, model=model)
    rng = np.random.default_rng(701)
    q4, v4 = (planar.from_raw_coefficients(rng.normal(size=planar.ndof)*1e-4) for _ in range(2))
    q7, v7 = embed(disc, planar, q4), embed(disc, planar, v4)
    active = ids(disc, planar)
    inactive = np.setdiff1d(np.arange(disc.ndof), active)
    pot4, pot7 = planar.potential(q4, hessian=True), disc.potential(q7, hessian=True)
    assert pot7["V"] == pytest.approx(pot4["V"], rel=2e-13)
    np.testing.assert_allclose(pot7["gradient"][active], pot4["gradient"], rtol=5e-13, atol=1e-15)
    np.testing.assert_allclose(pot7["hessian"][np.ix_(active, active)], pot4["hessian"], rtol=5e-13, atol=1e-11)
    np.testing.assert_allclose(disc.mass_matrix(q7)[np.ix_(active, active)], planar.mass_matrix(q4), rtol=5e-13, atol=5e-13)
    np.testing.assert_allclose(disc.inertial_terms(q7, v7)[active], planar.inertial_terms(q4, v4), rtol=1e-12, atol=1e-20)
    assert disc.energy(q7, v7) == pytest.approx(planar.energy(q4, v4), rel=2e-13)
    rhs = disc.rhs(0., np.r_[q7, v7])
    active_state = np.r_[active, disc.ndof+active]
    np.testing.assert_allclose(rhs[active_state], planar.rhs(0., np.r_[q4, v4]), rtol=1e-12, atol=2e-11)
    assert np.max(abs(rhs[np.r_[inactive, disc.ndof+inactive]])) == 0.
    np.testing.assert_allclose(disc.jacobian(0., np.r_[q7, v7])[np.ix_(active_state, active_state)],
                               planar.jacobian(0., np.r_[q4, v4]), rtol=5e-12, atol=2e-9)


def test_reflection_symmetry(disc):
    q, velocity = coordinates(disc)
    signs = np.repeat((1., 1., -1., -1., -1., 1., 1.), disc.n)
    reflected_q, reflected_v = signs*q, signs*velocity
    assert disc.energy(reflected_q, reflected_v) == pytest.approx(disc.energy(q, velocity), rel=3e-13)
    np.testing.assert_allclose(disc.acceleration(reflected_q, reflected_v), signs*disc.acceleration(q, velocity), rtol=2e-12, atol=1e-11)
    np.testing.assert_allclose(disc.mass_matrix(reflected_q), signs[:, None]*disc.mass_matrix(q)*signs[None, :], rtol=2e-12, atol=2e-13)


def test_linear_block_structure_and_frequency_units(disc):
    owners = {"u": 0, "c": 0, "w": 1, "theta": 1, "v": 2, "psi": 2, "Phi": 3}
    for first in FIELDS:
        for second in FIELDS:
            if owners[first] != owners[second]:
                assert np.max(abs(disc.K[disc.slices[first], disc.slices[second]])) == 0.
    modes = disc.linear_eigenpairs("torsion")
    np.testing.assert_allclose(modes["frequency_hz"]*2*np.pi, modes["omega"], rtol=2e-15)
    q, v = coordinates(disc)
    linear = disc.linear_reference(q, v, np.array((0., .03)))
    np.testing.assert_allclose(linear["q"][0], q, rtol=2e-11, atol=1e-15)
    np.testing.assert_allclose(linear["velocity"][0], v, rtol=2e-11, atol=1e-15)
    assert np.isfinite(linear["q"]).all()


def test_mass_cache_and_invalid_inputs(disc):
    q, v = coordinates(disc)
    disc.reset_counters()
    disc.mass_matrix(q); disc.mass_matrix(q.copy())
    assert disc.mass_factorizations == 1
    with pytest.raises(FloatingPointError, match="finite"):
        disc.rhs(0., np.full(2*disc.ndof, np.nan))
    with pytest.raises(ValueError, match="shape"):
        disc.project(np.zeros((disc.nq, 4)))
    with pytest.raises(ValueError, match="coordinates"):
        disc.reconstruct(np.zeros(4*disc.n))

