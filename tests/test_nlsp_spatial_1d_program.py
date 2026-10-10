"""Synthetic orchestration contracts; no real IVP or static equilibrium solve."""
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import nlsp_spatial_1d_program as program
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_spatial_dynamics import SpatialGalerkin

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def disc():
    path = ROOT/"results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json"
    if not path.exists():
        pytest.skip("Saved frozen quartic action unavailable")
    pol = json.loads(path.read_text(encoding="utf8"))["polynomials"]
    model = SimpleNamespace(T4=rod.Polynomial.deserialize(pol["T4"]), V4=rod.Polynomial.deserialize(pol["V4"]),
                            residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol["residuals_A"]))
    coefficients = rod.RodCoefficients.rectangular(1., 1., .3, .2, .1, 1.759089824002232e-5)
    return SpatialGalerkin(coefficients, 6, model=model)


def test_exact_inherited_policies():
    old = json.loads((ROOT/"data/input/nlsp_nonlinear_static_3d_fem.json").read_text(encoding="utf8"))
    assert program.STATIC_POLICY == {k: old["one_d"][k] for k in program.STATIC_POLICY}
    pilot = json.loads((ROOT/"data/input/weakly_nonlinear_planar_time_pilot.json").read_text(encoding="utf8"))
    assert program.TIME_POLICY == pilot["time_levels"]["tight"]
    assert program.SAFETY_POLICY == pilot["safety"]
    with pytest.raises(ValueError, match="unchanged"):
        program._policies({"time_policy": {**program.TIME_POLICY, "rtol": 1e-8}})


def test_canonical_component_vector_tolerances(disc):
    settings = program.time_settings(disc, .004)
    assert settings["atol"].shape == (2*disc.ndof,)
    assert settings["rtol"] == 1e-10 and np.unique(settings["atol"]).size > 1
    p = disc.coefficients
    masses = np.array((p.m, p.m, p.m, p.jp+p.jb, p.jb, p.jp, p.jp))
    expected = np.array((.004, .004, .004, .004, .004, .004, .004))*np.sqrt(masses/disc.n)
    np.testing.assert_allclose(settings["coordinate_scales"], expected, rtol=1e-15)
    np.testing.assert_allclose(settings["atol"][disc.ndof:], settings["atol"][:disc.ndof]*np.sqrt(p.C/p.jp), rtol=1e-15)


def test_two_axis_dead_load_no_torque(disc):
    force = program.line_load(disc, 2e-5, 2.5e-5)
    np.testing.assert_array_equal(force[disc.slices["w"]], 2e-5*(disc.B[1].T@disc.weights))
    np.testing.assert_array_equal(force[disc.slices["v"]], 2.5e-5*(disc.B[2].T@disc.weights))
    for f in ("u", "Phi", "psi", "theta", "c"):
        assert np.max(abs(force[disc.slices[f]])) == 0.


def test_static_solver_reuse_without_newton_execution(disc, monkeypatch):
    calls = []
    def synthetic_static(d, force, policy):
        calls.append((d, force.copy(), policy))
        return {"status": "PASS", "coordinate": np.zeros(d.ndof), "history": [], "load_factor": 1.}
    monkeypatch.setattr(program, "fem2_static_newton", synthetic_static)
    result = program.solve_static_pair(disc, (2e-5, 2.5e-5))
    assert len(calls) == 1 and calls[0][0] is disc and calls[0][2] == program.STATIC_POLICY
    assert result["linear_residual_relative"] < 1e-12
    np.testing.assert_array_equal(result["q_nonlinear"], np.zeros(disc.ndof))


def test_seven_field_safety_indices(disc):
    q = np.zeros(disc.ndof)
    q[disc.slices["Phi"].start] = disc.from_raw_coefficients(
        np.r_[np.zeros(3*disc.n), [1e-3], np.zeros(4*disc.n-1)])[disc.slices["Phi"].start]
    safe = program.safety_check(disc, q)
    assert safe["max_abs_theta"] > 0 and safe["max_abs_c"] == 0
    assert safe["mass_positive"]
    raw = np.zeros(disc.ndof); raw[disc.slices["c"].start] = -.2
    with pytest.raises(ArithmeticError, match="safety gate"):
        program.safety_check(disc, disc.from_raw_coefficients(raw))


class SyntheticDense:
    def __init__(self, old, end, initial):
        self.t_old, self.t, self.y_old = old, end, initial.copy()
        self.Q = np.zeros((len(initial), 3))

    def __call__(self, times):
        return np.repeat(self.y_old[:, None], len(np.atleast_1d(times)), axis=1)


class SyntheticRadau:
    """A synthetic accepted-prefix provider, never scipy's integrator."""
    def __init__(self, rhs, t0, y0, tbound, **settings):
        self.y, self.t, self.tbound, self.status = y0.copy(), t0, tbound, "running"
        self.nfev = self.njev = self.nlu = 0
        self.previous = t0

    def step(self):
        self.previous = self.t
        self.t = min(self.t+.05, self.tbound)
        self.status = "finished" if self.t == self.tbound else "running"

    def dense_output(self):
        return SyntheticDense(self.previous, self.t, self.y)


def test_accepted_prefix_dense_output_and_zero_velocity(disc, monkeypatch, tmp_path):
    monkeypatch.setattr(program, "Radau", SyntheticRadau)
    monkeypatch.setattr(program, "safety_check", lambda *a: {})
    q0 = np.zeros(disc.ndof); q0[disc.slices["w"].start] = 1e-6
    times = np.array((0., .03, .07, .1))
    history, records, stats = program.integrate_prefix(disc, q0, times, program.time_settings(disc, .004), {}, float("inf"))
    assert stats["status"] == "PASS" and stats["accepted_internal_steps"] == 2
    assert stats["actual_accepted_end"] == .1 and stats["time_end"] == .1
    np.testing.assert_array_equal(history[:, disc.ndof:], 0.)
    path = tmp_path/"dense.npz"
    program.save_dense_records(path, records)
    np.testing.assert_array_equal(program.evaluate_dense_records(path, times), history)


def test_no_last_value_padding_after_failure(disc, monkeypatch):
    class FailedSyntheticRadau(SyntheticRadau):
        def step(self):
            if self.t > 0:
                raise ArithmeticError("Synthetic failure after accepted prefix")
            super().step()
    monkeypatch.setattr(program, "Radau", FailedSyntheticRadau)
    monkeypatch.setattr(program, "safety_check", lambda *a: {})
    q0 = np.zeros(disc.ndof)
    history, records, stats = program.integrate_prefix(disc, q0, np.array((0., .03, .07, .1)),
        program.time_settings(disc, .004), {}, float("inf"))
    assert stats["status"] == "PARTIAL" and len(history) == 2 and len(records) == 1
    assert stats["time_end"] == .03 and stats["actual_accepted_end"] == .05
    assert "Synthetic failure" in stats["failure"]


def test_budget_before_radau_is_not_hidden_run(disc, monkeypatch):
    monkeypatch.setattr(program, "Radau", lambda *a, **kw: pytest.fail("Budget-exhausted attempt must not construct Radau"))
    history, records, stats = program.integrate_prefix(disc, np.zeros(disc.ndof), np.array((0., .1)),
        program.time_settings(disc, .004), {}, -1.)
    assert len(history) == 1 and not records and stats["nfev"] == 0 and stats["status"] == "PARTIAL"


def test_authorization_and_interrupted_attempt_guard(disc, tmp_path):
    with pytest.raises(ValueError, match="authorization"):
        program.run_case(disc, (1e-5, 1e-5), np.array((0., .1)), float("inf"), tmp_path/"case", authorization={})
    directory = tmp_path/"interrupted"; directory.mkdir()
    (directory/"attempt_ledger.json").write_text("{}")
    with pytest.raises(RuntimeError, match="no automatic retry"):
        program.cached_case(directory)


def test_failed_static_cache_never_reexecutes(disc, tmp_path, monkeypatch):
    calls = []
    def failed_static(d, load, policy):
        calls.append(1)
        return {"status": "FAIL", "force": program.line_load(d, *load), "q_linear": np.zeros(d.ndof),
                "q_nonlinear": np.zeros(d.ndof), "linear_residual_relative": 0.,
                "nonlinear": {"status": "FAIL", "reason": "Synthetic bounded failure"}}
    monkeypatch.setattr(program, "solve_static_pair", failed_static)
    monkeypatch.setattr(program, "Radau", lambda *a, **kw: pytest.fail("No dynamic solve after static failure"))
    auth = {"user_authorized_spatial_stage_b": True, "id": "synthetic_test"}
    directory = tmp_path/"case"
    first = program.run_case(disc, (1e-5, 1e-5), np.array((0., .1)), float("inf"), directory, authorization=auth)
    second = program.run_case(disc, (1e-5, 1e-5), np.array((0., .1)), float("inf"), directory, authorization=auth)
    assert len(calls) == 1 and first == second and first["calls"]["nonlinear_ODE"] == 0
    assert (directory/"static_states.npz").exists()
    with pytest.raises(ValueError, match="request differs"):
        program.run_case(disc, (2e-5, 1e-5), np.array((0., .1)), float("inf"), directory, authorization=auth)
    (directory/"static.json").write_text("{}")
    with pytest.raises(ValueError, match="artifact changed"):
        program.cached_case(directory)

