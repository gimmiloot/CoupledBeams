"""FEM-3A 1D source/contract tests: no ODE or modal jobs."""
from pathlib import Path
import json
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import nlsp_fem3a_1d_reference as helper

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def reference():
    return helper.load_reference(ROOT,
        "results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7",
        "results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad",
        "results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295")


@pytest.fixture
def pilot():
    return json.loads((ROOT / "data/input/weakly_nonlinear_planar_time_pilot.json").read_text(encoding="utf8"))


def test_saved_coordinates_not_reprojected(reference, monkeypatch):
    monkeypatch.setattr(helper.dynamics.PlanarGalerkin, "project", lambda *a: pytest.fail("no projection"))
    assert reference["saved"]["q_linear"].shape == (252,)
    assert reference["saved"]["q_nonlinear"].shape == (252,)
    assert not reference["saved"]["q_linear"].flags.writeable
    assert not reference["saved"]["q_nonlinear"].flags.writeable
    for kind in ("linear", "nonlinear"):
        np.testing.assert_array_equal(reference["disc"].raw_coefficients(reference["saved"]["q_" + kind]),
            reference["saved"]["raw_" + kind])


def test_saved_profiles_and_correct_field_order(reference):
    assert helper.FIELDS == ("u", "w", "theta", "c")
    assert reference["disc"].nq == 129
    for kind in ("linear", "nonlinear"):
        np.testing.assert_array_equal(reference["disc"].reconstruct(reference["saved"]["q_" + kind],
            reference["saved"]["s"]), reference["saved"][kind])


def test_frozen_action_loaded_without_derivation(monkeypatch):
    monkeypatch.setattr(helper.rod, "derive_polynomials", lambda *a: pytest.fail("no derivation"))
    monkeypatch.setattr(helper.dynamics.PlanarGalerkin, "linear_eigenpairs", lambda *a: pytest.fail("no eigensolve"))
    r = helper.load_reference(ROOT,
        "results/nlsp_nonlinear_static_3d_fem/6714bd9f2778e6d7",
        "results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad",
        "results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295")
    assert r["source"]["symbolic_derivations"] == 0
    assert r["source"]["static_equilibrium_solves"] == 0
    assert r["disc"].linear_eigendecompositions == 0


def test_release_acceleration_and_initial_velocity(reference, pilot, monkeypatch):
    monkeypatch.setattr(helper.dynamics.PlanarGalerkin, "linear_eigenpairs", lambda *a: pytest.fail("no eigensolve"))
    monkeypatch.setattr(helper.runner, "integrate_case", lambda *a, **k: pytest.fail("no ODE"))
    monkeypatch.setattr(helper.static, "fem2_static_newton", lambda *a, **k: pytest.fail("no static solve"))
    monkeypatch.setattr(helper.rod, "derive_polynomials", lambda *a: pytest.fail("no derivation"))
    import scipy.integrate
    monkeypatch.setattr(scipy.integrate, "Radau", lambda *a, **k: pytest.fail("no Radau"))
    result = helper.preflight_reference(reference, pilot)
    assert result["status"] == "PASS"
    assert result["external_force_after_release"] == 0
    assert result["admitted"] is False
    assert result["strict_float64_strong_weak"] == "PARTIAL"
    for kind in ("linear", "nonlinear"):
        state = result["initial_states"][kind]
        assert state["midspan_w_acceleration"] < 0
        assert state["mass_cholesky_positive"]
        assert np.count_nonzero(state["v0"]) == 0
        assert state["released_action_relative_residual"] < 2e-12
        assert state["loaded_equilibrium_relative_residual"] < 1e-10
        assert state["essential_endpoint_max_abs"] == 0
    assert result["slope_constraints"] is False


def test_correct_h_reference_period(reference):
    assert reference["omega1"] == 0.6054167303477958
    assert reference["T1"] == pytest.approx(10.37828159055014, rel=1e-15)
    assert reference["T1"] != 19.791162590151373


def test_tight_atol_uses_current_geometry(reference, pilot):
    helper.runner.load_runtime()
    config = helper.runtime_config(reference, pilot)
    assert config["material_geometry"]["h"] == .1
    assert pilot["material_geometry"]["h"] == .05
    settings = helper.runner.time_settings(reference["disc"], .005, "tight", config)
    assert settings["atol"].shape == (504,)
    assert settings["rtol"] == 1e-10
    assert settings["max_step"] == pytest.approx(.0072093929877006385, rel=1e-14)
    np.testing.assert_allclose(settings["atol"][::63],
        [8.908708063747483e-15, 8.908708063747483e-15,
         2.5717224993681997e-16, 2.5717224993681997e-16,
         3.235077240413132e-13, 3.235077240413132e-13,
         9.338863578008769e-15, 9.338863578008769e-15], rtol=1e-14)


def test_changed_tight_policy_rejected(reference, pilot):
    pilot["time_levels"]["tight"]["rtol"] = 1e-8
    with pytest.raises(ValueError, match="tight time policy"):
        helper.runtime_config(reference, pilot)


@pytest.mark.parametrize("times", [[.01, .1], [0., .1, .1], [0., np.nan], [0., .1, .09], [0., 1.]])
def test_actual_times_rejected(times):
    with pytest.raises(ValueError, match="Actual comparison times"):
        helper._times(times, .5)


def test_exact_linear_reuses_initial_state_without_ode(monkeypatch):
    q0 = np.array([.1, .2])
    calls = []
    def exact(q, v, times):
        calls.append((q.copy(), v.copy(), times.copy()))
        return {"times": times, "q": np.tile(q + 1e-16, (len(times), 1)),
            "velocity": np.zeros((len(times), 2))}
    disc = SimpleNamespace(ndof=2, linear_reference=exact, linear_eigendecompositions=1)
    ref = {"disc": disc, "saved": {"q_linear": q0}, "T1": 2.}
    trajectory = helper.exact_linear_reference(ref, [0., .05, .1])
    np.testing.assert_array_equal(trajectory["q"][0], q0)
    np.testing.assert_array_equal(calls[0][0], q0)
    assert np.count_nonzero(calls[0][1]) == 0
    assert trajectory["metadata"]["ODE_integrations"] == 0
    assert trajectory["metadata"]["new_1D_root_searches"] == 0
    assert trajectory["metadata"]["zero_time_factorization_roundoff_max_abs"] > 0


def test_unauthorized_nonlinear_run_rejected(reference, pilot):
    with pytest.raises(ValueError, match="authorization"):
        helper.integrate_nonlinear_reference(reference, [0., .05 * reference["T1"]], pilot,
            100., authorization={})


def test_only_full_authorized_horizon(reference, pilot):
    with pytest.raises(ValueError, match="authorized 0.05T1"):
        helper.integrate_nonlinear_reference(reference, [0., .01 * reference["T1"]], pilot,
            100., authorization={"user_authorized_FEM3A": True})


def test_mocked_radau_call_uses_exact_q0_and_saves_actual_prefix(reference, pilot, monkeypatch):
    calls = []
    def integrate(disc, shape, initial, config, amplitude_ratio, level, times, deadline, **kw):
        calls.append((shape, initial, config, amplitude_ratio, level, kw))
        hist = np.vstack([np.r_[kw["initial_coordinates"], np.zeros(disc.ndof)],
            np.r_[kw["initial_coordinates"] * .999, np.zeros(disc.ndof)]])
        return hist, {"status": "PARTIAL", "time_end": times[1]}
    monkeypatch.setattr(helper.runner, "integrate_case", integrate)
    horizon = .05 * reference["T1"]
    times = np.array([0., horizon / 2, horizon])
    result = helper.integrate_nonlinear_reference(reference, times, pilot, 100.,
        authorization={"user_authorized_FEM3A": True})
    assert len(calls) == 1
    shape, initial, config, ratio, level, kw = calls[0]
    assert shape is None and level == "tight"
    assert config["material_geometry"]["h"] == .1
    assert ratio == pytest.approx(.05)
    np.testing.assert_array_equal(kw["initial_coordinates"], reference["saved"]["q_nonlinear"])
    np.testing.assert_array_equal(result["times"], times[:2])
    np.testing.assert_array_equal(result["q"][0], reference["saved"]["q_nonlinear"])
    assert result["stats"]["status"] == "PARTIAL"
    assert result["stats"]["admitted"] is False
    assert result["stats"]["execution_mode"] == "EXPLORATORY_NOT_CERTIFIED"
    assert result["stats"]["no_dynamic_derivative_constraints"] is True


def test_artifact_integrity_without_automatic_recovery(tmp_path):
    (tmp_path / "manifest.json").write_text(json.dumps({"artifact_hashes": {"state.npz": "wrong"}}))
    (tmp_path / "state.npz").write_bytes(b"corrupted")
    with pytest.raises(ValueError, match="corrupted"):
        helper._checked_artifact(tmp_path, "state.npz")
