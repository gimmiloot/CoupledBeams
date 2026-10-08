"""Prepared IC checks using saved spectra; no ODE/eigensolves in tests.

Synthetic Hermite/endpoint states test algebra only. They are not admitted
production initial data and do not alter the historical zero-axial IVP.
"""
from __future__ import annotations
import ast
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from numpy.polynomial import legendre as leg
from numpy.polynomial import polynomial as power
import pytest

from scripts.lib import planar_prepared_initial_state as prep
from scripts.lib import planar_second_order_axial_response as leading
from scripts.lib import weakly_nonlinear_planar_dynamics as planar
from scripts.lib import weakly_nonlinear_spatial_rod as rod

ROOT = Path(__file__).resolve().parents[1]
SPECTRAL = ROOT/"results/planar_second_order_axial_response/b3ea4eb6ac95d6e1"
PILOT = ROOT/"results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef"
RECOVERY = ROOT/"results/weakly_nonlinear_planar_recovery/054874a4a4c9c9ff"
PREPARED = ROOT/"results/planar_prepared_initial_state/5ea8d41faf8ede54"
BASELINE_HEAD = "f60f14370713f84b9ada09ce83dae2d1357ec24f"
FROZEN = {
    "scripts/lib/weakly_nonlinear_spatial_rod.py": "aabc5a8657e56061df3d1c86f70801ad8fc0f1f355bc1f659f24ee62950a71f2",
    "scripts/lib/weakly_nonlinear_planar_dynamics.py": "eea98b77babcb4ed840e325efb270f37a68c32cff1d20ededf43c9989defc548",
}
ROUND_OFF = 2e-11


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _scaled_close(actual, expected, tolerance=ROUND_OFF):
    expected = np.asarray(expected)
    scale = max(float(np.max(abs(expected), initial=0)), 1e-30)
    np.testing.assert_allclose(actual, expected, rtol=tolerance, atol=tolerance*scale)


@pytest.fixture(autouse=True)
def no_hidden_integrators_or_eigensolves(monkeypatch):
    import scipy.integrate
    import scipy.linalg
    def forbidden(*args, **kwargs):
        raise AssertionError("Prepared-state tests must not integrate or solve a new eigensystem")
    for module, names in ((scipy.integrate, ("solve_ivp", "Radau")),
                          (scipy.linalg, ("eigh", "eig")),
                          (np.linalg, ("eigh", "eig", "eigvals")),
                          (leading, ("eigh",)), (planar, ("eigh",))):
        for name in names:
            monkeypatch.setattr(module, name, forbidden)


@pytest.fixture(scope="module")
def pilot_config():
    return _read(ROOT/"data/input/weakly_nonlinear_planar_time_pilot.json")


@pytest.fixture(scope="module")
def model():
    return rod.derive_polynomials()


@pytest.fixture(scope="module")
def coefficients(pilot_config):
    path = ROOT/pilot_config["audit_bundle"]/"result.json"
    manifest = _read(path.parent/"manifest.json")
    assert _sha(path) == manifest.get("artifact_hashes", manifest.get("artifacts", {}))["result.json"]
    return rod.RodCoefficients(**_read(path)["coefficients"])


@pytest.fixture(scope="module")
def background(pilot_config):
    return leading.background_from_pilot(pilot_config, ROOT)


@pytest.fixture(scope="module", params=(16, 96))
def cached_model(request, coefficients, background, model):
    p = request.param
    path = SPECTRAL/"models"/f"p{p}.npz"
    if not path.is_file():
        pytest.skip("Historical spectral arrays unavailable; no reproduction is launched")
    manifest = _read(SPECTRAL/"manifest.json")
    hashes = {key.replace(chr(92), "/"): value for key, value in manifest["artifact_hashes"].items()}
    assert _sha(path) == hashes[f"models/p{p}.npz"]
    return leading.SecondOrderAxial.from_saved(coefficients, p, background, path, model=model)


@pytest.fixture(scope="module")
def parts(cached_model):
    return prep.periodic_parts(cached_model)


def test_saved_restoration_retains_every_coordinate_without_new_eigh(cached_model):
    assert cached_model.eigen_decompositions == 0
    assert cached_model.vectors.shape == (2*(cached_model.p-1),)*2
    assert cached_model.b0.shape == cached_model.b2.shape == (cached_model.ndof,)
    assert cached_model.counters()["modal_reduction"] is False


def test_static_profile_satisfies_independent_matrix_equation(cached_model, parts):
    m = cached_model
    residual = m.K@parts["stat"]-m.f0
    scale = np.linalg.norm(m.K@parts["stat"])+np.linalg.norm(m.f0)
    assert np.linalg.norm(residual) <= ROUND_OFF*scale


def test_harmonic_profile_keeps_contraction_inertia(cached_model, parts):
    m = cached_model
    dynamic = m.K-m.driving_omega**2*m.M
    residual = dynamic@parts["harm"]-m.f2
    scale = np.linalg.norm(dynamic@parts["harm"])+np.linalg.norm(m.f2)
    assert np.linalg.norm(residual) <= ROUND_OFF*scale
    assert np.linalg.norm(m.K@parts["harm"]-m.f2) > 1e-8*np.linalg.norm(m.f2)


def test_spectral_profiles_match_independent_direct_linear_solves(cached_model, parts):
    m = cached_model
    _scaled_close(parts["stat"], np.linalg.solve(m.K, m.f0))
    _scaled_close(parts["harm"], np.linalg.solve(m.K-m.driving_omega**2*m.M, m.f2))


@pytest.mark.parametrize("fraction", (0., .001, .1, .25, 1., 5.))
def test_periodic_plus_free_reconstructs_historical_zero_ic_solution(cached_model, parts, fraction):
    m = cached_model
    t = fraction*m.background.T1
    periodic = parts["stat"]+parts["harm"]*np.cos(m.driving_omega*t)
    amplitudes = m.b0/m.omega**2+m.b2/(m.omega**2-m.driving_omega**2)
    free = -m.vectors@(np.cos(m.omega*t)*amplitudes)
    scale = max(np.linalg.norm(periodic)+np.linalg.norm(free), 1e-30)
    assert np.linalg.norm(periodic+free-m.evaluate([t])[0]) <= ROUND_OFF*scale


def test_free_part_is_required_for_old_zero_initial_conditions(cached_model, parts):
    m = cached_model
    amplitudes = m.b0/m.omega**2+m.b2/(m.omega**2-m.driving_omega**2)
    _scaled_close(parts["stat"]+parts["harm"], m.vectors@amplitudes)
    assert np.linalg.norm(parts["stat"]+parts["harm"]) > 0
    np.testing.assert_array_equal(m.evaluate([0.]), 0.)


@pytest.mark.parametrize("part", ("stat", "harm"))
@pytest.mark.parametrize("derivative", (0, 1, 2))
def test_physical_legendre_derivatives_match_reversible_shen_basis(cached_model, parts, part, derivative):
    m = cached_model
    physical = prep.physical_legendre_coefficients(m, parts[part])
    assert physical.shape == (2, m.p+1)
    profiles = prep.LegendreProfiles(physical, m.length)
    points = np.linspace(0., m.length, 71)
    expected = m.reconstruct(parts[part], points, derivative)
    actual = profiles.evaluate(points, derivative)
    scale = max(float(np.max(abs(expected))), 1e-30)
    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=2e-10*scale)


@pytest.mark.parametrize("part", ("stat", "harm"))
def test_periodic_profiles_enforce_values_not_derivative_clamps(cached_model, parts, part):
    m = cached_model
    profiles = prep.LegendreProfiles(prep.physical_legendre_coefficients(m, parts[part]), m.length)
    scale = np.max(abs(profiles.evaluate(np.linspace(0., m.length, 101))))
    assert np.max(abs(profiles.evaluate([0., m.length]))) <= ROUND_OFF*scale
    assert np.linalg.norm(profiles.evaluate([0., m.length], 1)) > 0


def test_common_amplitude_and_mode_normalization_are_not_refitted(background):
    values = background.evaluate([0., .25, .5, .75, 1.])
    np.testing.assert_allclose(values[2, 0], background.h0, rtol=2e-13, atol=0)
    np.testing.assert_array_equal(.05**2*values, 4*(.025**2*values))
    np.testing.assert_allclose(.05*values[:, 0], .0025*(values[:, 0]/background.h0), rtol=2e-15, atol=0)


def test_old_physics_and_initial_case_are_preserved():
    for relative, digest in FROZEN.items():
        assert _sha(ROOT/relative) == digest
    config = _read(ROOT/"data/input/weakly_nonlinear_planar_time_pilot.json")
    assert config["semantics"]["static_initial_correction"] is False
    assert config["boundary_conditions"].endswith("no slope constraints")


def test_actual_old_short_times_and_distinct_time_prescriptions():
    lo, hi = RECOVERY/"controls/new_p32", RECOVERY/"controls/p48_strict_short"
    if not (lo/"trajectory.npz").exists():
        pytest.skip("Historical short data unavailable")
    ma, mb = _read(lo/"case.json"), _read(hi/"case.json")
    with np.load(lo/"trajectory.npz") as a, np.load(hi/"trajectory.npz") as b:
        np.testing.assert_array_equal(a["time"], b["time"])
        assert a["time"][-1] == ma["time_end"] == mb["time_end"]
        assert len(a["time"]) == ma["samples"] == mb["samples"]
        assert np.all(np.diff(a["time"]) > 0)
    assert ma["p"] == 32 and mb["p"] == 48
    assert ma["time_level"] == "tight" and mb["time_level"] == "allowed_extra"


def test_old_partial_case_is_not_extended_by_duplicate_snapshots():
    case = PILOT/"cases/p24_Aoverh0p025_tight"
    if not (case/"trajectory.npz").exists():
        pytest.skip("Historical partial data unavailable")
    metadata = _read(case/"case.json")
    with np.load(case/"trajectory.npz") as data:
        assert data["time"][-1] == metadata["time_end"] < metadata["target_time_end"]
        assert data["time"][data["snapshot_indices"]].max() == data["time"][-1]
        assert len(np.unique(data["snapshot_indices"])) < len(data["snapshot_indices"])
    assert metadata["status"] == "PARTIAL"



def test_profile_gate_distinguishes_small_values_from_unresolved_endpoint_jets():
    from scripts.analysis import prepare_planar_initial_state as cli
    low = np.zeros((2, 97))
    low[:, 0] = [.02, .10]
    low[:, 2] = -low[:, 0]
    high = low.copy()
    high[0, 94] += 1e-12
    high[0, 96] -= 1e-12
    policy = {"relative_tolerance": 1e-6, "endpoint_tolerance": 1e-6,
              "spatial_derivatives": [0, 1, 2]}
    result = cli.profile_comparison(low, high, 1., .05, policy)
    values = next(row for row in result["rows"] if row["field"] == "u" and row["derivative"] == 0)
    jets = next(row for row in result["rows"] if row["field"] == "u" and row["derivative"] == 2)
    assert values["relative_L2"] < policy["relative_tolerance"]
    assert jets["endpoint_fixed_scaled_difference"] > policy["endpoint_tolerance"]
    assert result["pass"] is False
    assert result["derivatives_from_Legendre_not_PDE"] is True


def _synthetic_cache(folder, identity):
    """Small data-cache fixture, never an accepted prepared-state result."""
    from scripts.analysis import prepare_planar_initial_state as cli
    folder.mkdir(parents=True)
    summary = {"synthetic_fixture": True,
               "statuses": {"NLSP_COMMON_INITIAL_PROJECTION": "PARTIAL",
                            "NLSP_PREPARED_SHORT_SPATIAL_CHECK": "NOT_RUN",
                            "NLSP_PREPARED_SHORT_TEMPORAL_CHECK": "NOT_RUN"},
               "new_ODE_integrations": 0}
    cli.write_json(folder/"summary.json", summary)
    points = np.linspace(0., 1., 13)
    values = np.column_stack((points*(1-points)*1e-3,
                              points*(1-points)*2e-3))
    cli.save_npz(folder/"profiles.npz", s=points, p96_stat=values, p96_harm=-values*.3)
    rows = [{"field": field, "derivative": d, "relative_L2": 1e-8,
             "endpoint_fixed_scaled_difference": 1e-8}
            for field in ("u", "c") for d in (0, 1, 2)]
    cli.write_json(folder/"profile_convergence.json",
                   {"pairs": [{"low_p": 64, "high_p": 96, "stat": {"rows": rows},
                               "harm": {"rows": rows}}]})
    files = {str(path.relative_to(folder)): _sha(path) for path in folder.rglob("*") if path.is_file()}
    cli.write_json(folder/"manifest.json", {"identity": identity, "artifact_hashes": files})
    return summary


@pytest.mark.parametrize("action", ("compute", "report-only", "plot-only"))
def test_cached_entrypoints_perform_zero_preparation_and_zero_eigen_or_time_solves(
        tmp_path, monkeypatch, capsys, action):
    from scripts.analysis import prepare_planar_initial_state as cli
    identity = {"synthetic_cache": "prepared-state-v1"}
    output = tmp_path/"results"
    bundle = output/"fixture"
    expected = _synthetic_cache(bundle, identity)

    def forbidden(*args, **kwargs):
        raise AssertionError("Cached workflow attempted preparation, BVP, history or integration")

    monkeypatch.setattr(cli, "run_compute", forbidden)
    monkeypatch.setattr(prep, "periodic_parts", forbidden)
    monkeypatch.setattr(prep, "quintic_theta3", forbidden)
    monkeypatch.setattr(prep.PreparedInitialState, "evaluate", forbidden)
    monkeypatch.setattr(leading.SecondOrderAxial, "from_saved", forbidden)
    monkeypatch.setattr(rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(cli, "identity", lambda *args: ("fixture", identity))
    args = (["--compute", "--output-dir", str(output)] if action == "compute"
            else ["--"+action, str(bundle)])
    returned = cli.main(args)
    reported = json.loads(capsys.readouterr().out)
    assert returned == expected
    assert all(value == 0 for value in reported["this_run_counters"].values())
    assert cli.validate_cache(bundle) == expected


def test_repeat_plot_preserves_figure_hashes_and_immutable_cache(tmp_path):
    from scripts.analysis import prepare_planar_initial_state as cli
    identity = {"synthetic_cache": "deterministic-figures"}
    bundle = tmp_path/"bundle"
    expected = _synthetic_cache(bundle, identity)
    cli.plot_only(bundle)
    cli.write_json(bundle/"manifest.json", cli.manifest_for(bundle, identity))
    before = {str(path.relative_to(bundle)): _sha(path) for path in (bundle/"figures").iterdir()}
    cli.plot_only(bundle)
    after = {str(path.relative_to(bundle)): _sha(path) for path in (bundle/"figures").iterdir()}
    assert before == after
    assert cli.validate_cache(bundle, identity) == expected


def test_cache_rejects_corrupt_data_and_wrong_identity(tmp_path):
    from scripts.analysis import prepare_planar_initial_state as cli
    bundle = tmp_path/"bundle"
    identity = {"synthetic_cache": "immutable"}
    _synthetic_cache(bundle, identity)
    with pytest.raises(ValueError, match="identity mismatch"):
        cli.validate_cache(bundle, {"synthetic_cache": "changed"})
    (bundle/"summary.json").write_text('{"synthetic_fixture":"corrupt"}', encoding="utf8")
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        cli.validate_cache(bundle)


@pytest.mark.parametrize("changed", ("amplitude", "theta3", "reference", "projection"))
def test_cache_identity_changes_with_common_state_and_numerical_policy(tmp_path, changed):
    from scripts.analysis import prepare_planar_initial_state as cli
    original = _read(cli.CONFIG)
    source = tmp_path/"config.json"
    cli.write_json(source, original)
    base, _ = cli.identity(source)
    config = json.loads(json.dumps(original))
    if changed == "amplitude":
        config["amplitude_over_h"] = .025
    elif changed == "theta3":
        config["theta3_rule"] += " synthetic test change"
    elif changed == "reference":
        config["profile_policy"]["final_pair"] = [48, 64]
    else:
        config["nonlinear_policy"]["primary_pair"] = [48, 64]
    cli.write_json(source, config)
    altered, _ = cli.identity(source)
    assert altered != base


def test_profile_policy_preserves_old_trajectory_acceptance(pilot_config):
    from scripts.analysis import prepare_planar_initial_state as cli
    config = _read(cli.CONFIG)
    assert config["profile_policy"]["relative_tolerance"] == 1e-6
    assert config["profile_policy"]["endpoint_tolerance"] == 1e-6
    assert pilot_config["gates"]["u_c_relative"] == 1e-3
    assert pilot_config["gates"]["w_theta_relative"] == 1e-4
    assert pilot_config["gates"]["energy_relative_drift"] == 1e-6
    assert config["nonlinear_policy"]["primary_pair"] == [32, 48]
    assert config["nonlinear_policy"]["allowed_pre_run_replacement"] == [48, 64]
    assert config["nonlinear_policy"]["maximum_integrations"] == 3


def test_test_module_has_only_documented_zero_duration_mock_initializer_call():
    tree = ast.parse(Path(__file__).read_text(encoding="utf8"))
    forbidden = {"solve_ivp", "Radau", "eigh", "eig", "eigvals"}
    calls = []
    for outer in tree.body:
        for node in ast.walk(outer):
            if isinstance(node, ast.Call):
                name = node.func.id if isinstance(node.func, ast.Name) else (
                    node.func.attr if isinstance(node.func, ast.Attribute) else None)
                assert name not in forbidden
                if name == "integrate_case":
                    calls.append(getattr(outer, "name", None))
    # One forwarding unit test adapter calls a patched, already-finished
    # initializer at t0=t_end=0. It never advances or evaluates an ODE.
    assert calls == ["_mock_initializer_call"]



@pytest.fixture(scope="module")
def candidate(background, coefficients, model):
    """ONE common p96 source profile; p96 is not a nonlinear pilot space."""
    path = SPECTRAL/"models/p96.npz"
    if not path.is_file():
        pytest.skip("Common historical source profile unavailable")
    hashes = {key.replace(chr(92), "/"): value
              for key, value in _read(SPECTRAL/"manifest.json")["artifact_hashes"].items()}
    assert _sha(path) == hashes["models/p96.npz"]
    source = leading.SecondOrderAxial.from_saved(coefficients, 96, background, path, model=model)
    split = prep.periodic_parts(source)
    profiles = prep.LegendreProfiles(
        prep.physical_legendre_coefficients(source, split["stat"])
        +prep.physical_legendre_coefficients(source, split["harm"]),
        background.length)
    endpoint = np.array([0., background.length])
    uc_first = profiles.evaluate(endpoint, 1)
    bending_first = background.evaluate(endpoint, 1)
    correction = prep.quintic_theta3(uc_first[:, 0], bending_first[:, 1],
                                    bending_first[:, 0], coefficients, background.length)
    return prep.PreparedInitialState(background, profiles, correction, admitted=False)


@pytest.fixture(scope="module")
def allowed_nonlinear_spaces(coefficients, background, model):
    """Allowed audit spaces only; no time trajectories or eigenproblems."""
    return {p: planar.PlanarGalerkin(coefficients, p, length=background.length, model=model)
            for p in (32, 48, 64)}


def test_endpoint_trace_identities_match_independent_frozen_a_and_b_residuals(model):
    p = model.symbols
    expected = {
        "u": p["C"]*p["u_ss"]+p["nu"]*p["C"]*p["c_s"]+(p["C"]-p["S"])*p["theta_s"]*p["w_s"],
        "w": p["S"]*(p["w_ss"]-p["theta_s"])+(p["C"]-p["S"])*p["u_s"]*p["theta_s"],
        "theta": p["Bp"]*p["theta_ss"]+p["S"]*p["w_s"]-(p["C"]-p["S"])*p["u_s"]*p["w_s"],
        "c": p["H"]*p["c_ss"]-p["nu"]*p["C"]*p["u_s"],
    }
    zero = {name+suffix: 0 for name in rod.FIELD_ORDER
            for suffix in ("", "_t", "_st", "_tt")}
    zero.update({name+suffix: 0 for name in ("v", "Phi", "psi")
                 for suffix in ("_s", "_ss")})
    audit = prep.endpoint_trace_audit(model)
    assert audit["status"] == "PASS"
    for field, index in zip(("u", "w", "theta", "c"), (0, 1, 5, 6)):
        assert rod.Polynomial.deserialize(audit["trace_polynomials"][field]) == expected[field]
        assert -model.residual_a[index].substitute(zero) == expected[field]
        assert -model.residual_b[index].substitute(zero) == expected[field]
    assert all(row["difference_terms"] == [] for row in audit["checks"].values())


def _synthetic_compatible_jets(coefficients):
    """Endpoint algebra fixture, not a continuous production BVP solution."""
    p = coefficients
    ws, ts = np.array([.08, -.08]), np.array([.11, .11])
    us, cs = np.array([.012, .018]), np.array([.002, -.004])
    background_jets = {"value": np.zeros((2, 2)),
                       "first": np.column_stack((ws, ts)),
                       "second": np.column_stack((ts, -p.S/p.Bp*ws))}
    profiles_jets = {"value": np.zeros((2, 2)),
                     "first": np.column_stack((us, cs)),
                     "second": np.column_stack((-p.nu*cs-(p.C-p.S)/p.C*ts*ws,
                                                p.nu*p.C/p.H*us))}
    correction = prep.quintic_theta3(us, ts, ws, p)
    return background_jets, profiles_jets, correction


def test_o2_preparation_cancels_axial_and_contraction_but_leaves_cubic_bending(coefficients):
    b, uc, correction = _synthetic_compatible_jets(coefficients)
    no_correction = {name: np.zeros(2) for name in ("value", "first", "second")}
    formal = prep.formal_endpoint_coefficients(b, uc, no_correction, coefficients)
    for field in ("u", "c"):
        assert np.max(abs(formal[2][field])) < 2e-16
    assert np.linalg.norm(formal[3]["w"]) > 0
    assert np.linalg.norm(formal[3]["theta"]) > 0
    original_uc = {name: np.zeros((2, 2)) for name in ("value", "first", "second")}
    old = prep.formal_endpoint_coefficients(b, original_uc, no_correction, coefficients)
    assert np.linalg.norm(old[2]["u"]) > 0  # Historical mismatch is retained.


@pytest.mark.parametrize("length", (1., .73, 2.))
def test_quintic_meets_six_independent_endpoint_jet_conditions(coefficients, length):
    us, ts, ws = np.array([.012, .018]), np.array([.11, .14]), np.array([.08, -.06])
    q = prep.quintic_theta3(us, ts, ws, coefficients, length)
    target_first = (coefficients.C-coefficients.S)/coefficients.S*us*ts
    target_second = (coefficients.C-coefficients.S)/coefficients.Bp*us*ws
    endpoints = np.array([0., length])
    np.testing.assert_allclose(q.evaluate(endpoints), 0., atol=2e-14, rtol=0)
    _scaled_close(q.evaluate(endpoints, 1), target_first)
    _scaled_close(q.evaluate(endpoints, 2), target_second)
    # Direct six-by-six Hermite solve in s/L provides an independent polynomial.
    rows = []
    for location in (0., 1.):
        for derivative in (0, 1, 2):
            rows.append([power.polyval(location, power.polyder(np.eye(6)[k], m=derivative))
                         for k in range(6)])
    direct = np.linalg.solve(np.array(rows), np.array(
        [0., length*target_first[0], length**2*target_second[0],
         0., length*target_first[1], length**2*target_second[1]]))
    _scaled_close(q.eta_power_coefficients, direct)
    for derivative in (0, 1, 2):
        x = np.linspace(0, length, 43)
        expected = power.polyval(x/length, power.polyder(direct, m=derivative))/length**derivative
        _scaled_close(q.evaluate(x, derivative), expected)


def test_quintic_hermite_basis_is_unique_and_exact_in_rational_arithmetic():
    result = prep.hermite_basis_audit()
    assert result["status"] == "PASS" and result["identities"] == 24
    assert result["full_six_condition_system_determinant"] != 0


def test_symmetric_first_mode_theta3_has_signed_reflection_parity(coefficients):
    q = prep.quintic_theta3([.012, .012], [.11, .11], [.08, -.08], coefficients)
    points = np.linspace(0., 1., 57)
    scale = max(float(np.max(abs(q.evaluate(points)))), 1e-30)
    np.testing.assert_allclose(q.evaluate(1-points), -q.evaluate(points),
                               rtol=ROUND_OFF, atol=ROUND_OFF*scale)
    assert q.parity_metrics()["coefficient_scaled_defect"] < ROUND_OFF


def test_theta3_cancels_through_cubic_order_without_erasing_orders_four_and_five(coefficients):
    b, uc, correction = _synthetic_compatible_jets(coefficients)
    formal = prep.formal_endpoint_coefficients(b, uc, correction.endpoint_jets(), coefficients)
    for order in (1, 2, 3):
        assert max(np.max(abs(value)) for value in formal[order].values()) < 2e-16
    assert np.linalg.norm(formal[4]["u"]) > 0
    assert np.linalg.norm(formal[5]["w"]) > 0
    np.testing.assert_array_equal(formal[4]["theta"], 0.)
    np.testing.assert_array_equal(formal[5]["c"], 0.)


@pytest.mark.parametrize("epsilon", (.05, .025, .0125))
def test_finite_amplitude_endpoint_residual_matches_full_cubic_action_not_o3_truncation(
        coefficients, model, epsilon):
    b, uc, correction = _synthetic_compatible_jets(coefficients)
    jets = correction.endpoint_jets()
    formal = prep.formal_endpoint_coefficients(b, uc, jets, coefficients)
    force = []
    for endpoint in range(2):
        values = {name: 0. for name in rod.SYMBOL_ORDER}
        values.update(coefficients.values())
        for field, i in (("u", 0), ("c", 1)):
            values[field+"_s"] = epsilon**2*uc["first"][endpoint, i]
            values[field+"_ss"] = epsilon**2*uc["second"][endpoint, i]
        values["w_s"], values["w_ss"] = epsilon*b["first"][endpoint, 0], epsilon*b["second"][endpoint, 0]
        values["theta_s"] = epsilon*b["first"][endpoint, 1]+epsilon**3*jets["first"][endpoint]
        values["theta_ss"] = epsilon*b["second"][endpoint, 1]+epsilon**3*jets["second"][endpoint]
        force.append([-model.residual_a[index].evaluate(values) for index in (0, 1, 5, 6)])
    force = np.asarray(force)
    expected = np.column_stack([sum(epsilon**order*formal[order][field] for order in formal)
                                for field in ("u", "w", "theta", "c")])
    np.testing.assert_allclose(force, expected, rtol=ROUND_OFF, atol=2e-18)
    assert np.max(abs(force[:, 0])) > 1e-14
    assert np.max(abs(expected[:, 1])) > 1e-18
    assert np.linalg.norm(force) > 0


def test_unadmitted_candidate_cannot_be_used_as_dynamical_initial_state(candidate):
    with pytest.raises(RuntimeError, match="not admitted"):
        candidate.evaluate(np.array([.25, .5]), .05)
    diagnostic = candidate.evaluate(np.array([.25, .5]), .05, require_admitted=False)
    assert diagnostic.shape == (2, 4) and np.all(np.isfinite(diagnostic))


@pytest.mark.parametrize("derivative", (0, 1, 2))
def test_one_common_initial_evaluator_retains_independent_fields_and_amplitude_orders(candidate, derivative):
    points = np.linspace(0., candidate.length, 47)
    eps = .05
    fields = candidate.evaluate(points, eps, derivative, require_admitted=False)
    uc, bending = candidate.profiles.evaluate(points, derivative), candidate.background.evaluate(points, derivative)
    theta3 = candidate.correction.evaluate(points, derivative)
    _scaled_close(fields[:, 0], eps**2*uc[:, 0])
    _scaled_close(fields[:, 3], eps**2*uc[:, 1])
    _scaled_close(fields[:, 1], eps*bending[:, 0])
    _scaled_close(fields[:, 2], eps*bending[:, 1]+eps**3*theta3)
    np.testing.assert_array_equal(candidate.initial_velocities(points), 0.)
    np.testing.assert_array_equal(candidate.evaluate(points, 0., derivative, require_admitted=False), 0.)


def test_physical_profiles_are_copied_and_readonly_not_mutated_to_fit_projection():
    values = np.array([[.1, 0., -.1], [.2, 0., -.2]])
    profiles = prep.LegendreProfiles(values)
    values[0, 0] = 10.
    assert profiles.coefficients[0, 0] == .1
    with pytest.raises(ValueError):
        profiles.coefficients[0, 0] = 10.


def test_current_result_blocks_both_allowed_pairs_when_any_field_jet_is_unresolved():
    path = PREPARED/"summary.json"
    if not path.is_file():
        pytest.skip("Historical prepared admission result is not locally available")
    assert _sha(path) == _read(PREPARED/"manifest.json")["artifact_hashes"]["summary.json"]
    summary = _read(path)
    if "projection" not in summary:
        pytest.skip("Profile preparation stopped before projection stage")
    projection, decision = summary["projection"], summary["pre_run_decision"]
    chosen = decision["selected_pair"]
    for pair in (decision["primary_pair"], decision["allowed_replacement"]):
        if not all(projection[str(p)]["pass"] for p in pair):
            assert chosen != pair
    if chosen is None:
        assert summary["new_ODE_integrations"] == 0
        assert summary["actual_short_runs"] == []
        assert summary["statuses"]["NLSP_PREPARED_SHORT_TEMPORAL_CHECK"] == "NOT_RUN"
        assert summary["statuses"]["NLSP_PREPARED_SHORT_SPATIAL_CHECK"] == "NOT_RUN"
        assert summary["old_task_status"] == summary["old_recovery_status"] == "PARTIAL"
    # Passing O2 u/c traces does not override a failed bending derivative gate.
    for p, result in projection.items():
        if not all(row["pass"] for row in result["rows"]):
            assert result["pass"] is False


@pytest.mark.parametrize("p", (32, 48, 64))
def test_initial_energy_and_mass_bounds_use_preserved_quartic_model(candidate, allowed_nonlinear_spaces, pilot_config, monkeypatch, p):
    from scripts.analysis import prepare_planar_initial_state as cli
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as old
    disc = allowed_nonlinear_spaces[p]
    fields = candidate.evaluate(disc.x, .05, require_admitted=False)
    q = disc.project(fields)
    velocity = np.zeros(disc.ndof)
    energy = disc.energy(q, velocity)
    bounds = cli.mass_safety_bounds(disc, q)
    assert energy > 0 and np.isfinite(energy)
    assert bounds["mass_positive"] and bounds["relative_mass_eigenvalue_lower_bound"] > 0
    assert bounds["relative_mass_condition_upper_bound"] >= 1
    monkeypatch.setattr(old, "np", np, raising=False)
    old.safety_check(disc, q, pilot_config["safety"])
    independent = 0.
    values = disc.reconstruct(q)
    first = disc.reconstruct(q, derivative=1)
    second = disc.reconstruct(q, derivative=2)
    for i, weight in enumerate(disc.weights):
        data = np.zeros((6, 7))
        data[0, [0, 1, 5, 6]] = values[i]
        data[1, [0, 1, 5, 6]] = first[i]
        data[3, [0, 1, 5, 6]] = second[i]
        evaluated = rod.polynomial_evaluate(rod.FieldJet(*data), disc.coefficients, disc.model)
        independent += weight*evaluated["V"]
    _scaled_close(energy, independent)



@pytest.fixture(scope="module")
def common_weak_diagnostics(candidate, allowed_nonlinear_spaces):
    """Keep absolute/scaled evidence separate on one shared initial evaluator."""
    rows = {}
    for p, d in allowed_nonlinear_spaces.items():
        q = d.project(candidate.evaluate(d.x, .05, require_admitted=False))
        # Synthetic velocity tests inertia algebra; the real initial velocity is zero.
        velocity = d.project(np.column_stack((d.x*(1-d.x)*1e-8,
                                              d.x*(1-d.x)*1e-7,
                                              d.x*(1-d.x)*1e-8,
                                              d.x*(1-d.x)*1e-9)))
        acceleration = d.acceleration(q, velocity)
        gradient = d.potential(q)["gradient"]
        mass_term = d.mass_matrix(q)@acceleration
        inertia = d.inertial_terms(q, velocity)
        action = mass_term+inertia+gradient
        weak = d.weak_residual(q, velocity, acceleration)
        local = d._local_potential(q)[1]
        work = (sum(np.linalg.norm(matrix.T@(d.weights*values))
                    for matrix, values in zip(d._potential_matrices, local))
                +np.linalg.norm(mass_term)+np.linalg.norm(inertia))
        rows[p] = {"difference": weak-action, "work": work,
                   "energy_rate": d.energy_rate(q, velocity, acceleration),
                   "power_scale": np.linalg.norm(velocity)*work}
    return rows


@pytest.mark.parametrize("p", (32, 48, 64))
def test_common_state_absolute_action_weak_and_energy_power_gates(
        common_weak_diagnostics, pilot_config, p):
    row = common_weak_diagnostics[p]
    assert np.max(abs(row["difference"])) <= 2e-12
    assert abs(row["energy_rate"]) <= pilot_config["gates"]["identity_scaled"]*max(row["power_scale"], 1e-30)


@pytest.mark.parametrize("p", (
    32,
    pytest.param(48, marks=pytest.mark.xfail(
        strict=True, raises=AssertionError,
        reason="Recorded common-p96 projected p48 relative weak/action numerical check unresolved; unchanged 2e-12 gate")),
    pytest.param(64, marks=pytest.mark.xfail(
        strict=True, raises=AssertionError,
        reason="Recorded common-p96 projected p64 relative weak/action numerical check unresolved; unchanged 2e-12 gate")),
))
def test_common_state_relative_action_weak_gate_remains_unresolved_when_recorded(
        common_weak_diagnostics, p):
    row = common_weak_diagnostics[p]
    assert np.linalg.norm(row["difference"]) <= 2e-12*max(row["work"], 1e-30)


def test_recorded_auxiliary_p96_failure_is_preserved_as_a_qualification():
    path = PREPARED/"auxiliary_weak_identity.json"
    if not path.is_file():
        pytest.skip("Historical auxiliary numerical check not available")
    assert _sha(path) == _read(PREPARED/"manifest.json")["artifact_hashes"]["auxiliary_weak_identity.json"]
    record = _read(path)
    assert record["same_common_p96_initial_evaluator"] is True
    assert record["ODE_integrations"] == record["mh_timoshenko_eigensolves"] == 0
    rows = {row["p"]: row for row in record["rows"]}
    for p in (32, 48, 64, 96):
        row = rows[p]
        assert row["absolute_gate"] == row["relative_gate"] == 2e-12
        assert row["absolute_pass"] is True
    assert rows[32]["relative_pass"] is True
    for p in (48, 64, 96):
        assert rows[p]["relative_pass"] is False
        assert rows[p]["relative_work_residual"] > rows[p]["relative_gate"]
    assert record["nonlinear_p96_not_authorized"] is True


@pytest.mark.parametrize("frequencies", ([1.], [1., 3.]))
def test_resonant_periodic_preparation_is_explicitly_unresolved_and_json_safe(
        frequencies, monkeypatch):
    from types import SimpleNamespace
    omega = np.array(frequencies)
    saved = SimpleNamespace(omega=omega.copy(), vectors=np.eye(len(omega)),
                            driving_omega=1., M=np.eye(len(omega)), K=np.diag(omega**2),
                            b0=np.ones(len(omega)), b2=np.ones(len(omega)),
                            f0=np.ones(len(omega)), f2=np.ones(len(omega)))

    def forbidden(*args, **kwargs):
        raise AssertionError("A singular preparation must not be regularized or solved")

    monkeypatch.setattr(np.linalg, "solve", forbidden)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        result = prep.periodic_parts(saved, direct_check=True)
    assert result["status"] == "PREPARATION_UNRESOLVED"
    assert "stat" not in result and "harm" not in result
    assert result["checks"]["dynamic_operator_singular_or_unresolved"] is True
    assert result["checks"]["modal_dynamic_condition"] is None
    assert result["checks"]["minimum_absolute_detuning"] == 0.
    np.testing.assert_array_equal(saved.omega, omega)
    json.dumps(result, allow_nan=False)


def test_missing_mandatory_source_stops_before_any_preparation_or_history(
        tmp_path, monkeypatch):
    from scripts.analysis import prepare_planar_initial_state as cli
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as old
    config = _read(cli.CONFIG)
    config["spectral_bundle"] = str(tmp_path/"missing_spectral")

    def forbidden(*args, **kwargs):
        raise AssertionError("Missing mandatory source must stop before numerical preparation")

    monkeypatch.setattr(old, "load_runtime", lambda: None)
    monkeypatch.setattr(leading.SecondOrderAxial, "from_saved", forbidden)
    monkeypatch.setattr(rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(cli.previous, "historical_provenance", forbidden)
    summary = cli.run_compute(config, tmp_path/"stopped")
    assert summary["stop_reason"].startswith("DATA_UNAVAILABLE: mandatory")
    assert summary["new_ODE_integrations"] == summary["new_eigendecompositions"] == 0
    assert summary["statuses"]["NLSP_PREPARED_SHORT_SPATIAL_CHECK"] == "NOT_RUN"
    assert summary["statuses"]["NLSP_PREPARED_SHORT_TEMPORAL_CHECK"] == "NOT_RUN"
    assert _read(tmp_path/"stopped/summary.json") == summary


def test_expired_preparation_budget_saves_partial_before_restoring_models(
        tmp_path, monkeypatch, model, background):
    from types import SimpleNamespace
    from scripts.analysis import prepare_planar_initial_state as cli
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as old
    config = _read(cli.CONFIG)
    values = iter([0., 1000., 1000.])

    def forbidden(*args, **kwargs):
        raise AssertionError("Expired budget must not restore spectra or evaluate histories")

    monkeypatch.setattr(cli, "time", SimpleNamespace(perf_counter=lambda: next(values)))
    monkeypatch.setattr(old, "load_runtime", lambda: None)
    monkeypatch.setattr(cli.previous, "historical_provenance", lambda *args: {})
    monkeypatch.setattr(cli.previous, "validate_cache", lambda *args: {})
    monkeypatch.setattr(leading, "background_from_pilot", lambda *args: background)
    monkeypatch.setattr(rod, "derive_polynomials", lambda: model)
    monkeypatch.setattr(leading.SecondOrderAxial, "from_saved", forbidden)
    summary = cli.run_compute(config, tmp_path/"budget_stop")
    assert summary["stop_reason"] == "PREPARATION_BUDGET_EXHAUSTED"
    assert summary["details"]["completed_degrees"] == []
    assert summary["new_ODE_integrations"] == summary["new_eigendecompositions"] == 0
    assert summary["statuses"]["NLSP_PREPARED_INITIAL_STATE_PILOT"] == "PARTIAL"
    assert summary["statuses"]["NLSP_PREPARED_SHORT_TEMPORAL_CHECK"] == "NOT_RUN"
    assert (tmp_path/"budget_stop/periodic_checks.json").is_file()



# Bounded precision / exploratory short-control continuation.  The historical
# 5ea relative weak/action XFAILs above remain tied to that earlier DP projection.
@pytest.fixture(scope="module")
def frozen_prepared_state():
    if not (PREPARED/"manifest.json").is_file():
        pytest.skip("Frozen prepared target is absent; no target is regenerated")
    return prep.load_frozen_prepared_state(PREPARED)


def test_frozen_target_loader_does_not_refit_profiles_theta3_or_source_eigenpair(monkeypatch):
    if not (PREPARED/"manifest.json").is_file():
        pytest.skip("Frozen physical target absent")
    def forbidden(*args, **kwargs):
        raise AssertionError("Target loading must not solve or refit anything")
    monkeypatch.setattr(prep, "periodic_parts", forbidden)
    monkeypatch.setattr(prep, "quintic_theta3", forbidden)
    monkeypatch.setattr(leading, "background_from_pilot", forbidden)
    state, coeff, provenance = prep.load_frozen_prepared_state(PREPARED)
    summary = _read(PREPARED/"summary.json")
    with np.load(PREPARED/"common_initial_state.npz", allow_pickle=False) as saved:
        np.testing.assert_array_equal(state.profiles.coefficients, saved["U_C_legendre"])
        np.testing.assert_array_equal(state.correction.eta_power_coefficients, saved["theta3_eta_power"])
    np.testing.assert_array_equal(state.background.normalized_coefficients,
                                  summary["background"]["normalized_analytic_coefficients"])
    assert state.background.omega == summary["background"]["omega1"]
    assert coeff.values() == summary["coefficients"]
    assert state.admitted is False
    assert provenance["regenerated_profiles"] is False
    assert provenance["regenerated_Theta3"] is False
    assert provenance["new_eigendecompositions"] == 0
    assert provenance["manifest_sha256"] == _sha(PREPARED/"manifest.json")


@pytest.mark.parametrize("artifact", ("summary.json", "common_initial_state.npz"))
def test_frozen_target_loader_rejects_each_corrupt_physical_input(tmp_path, artifact):
    if not PREPARED.is_dir():
        pytest.skip("Frozen physical target absent")
    for name in ("manifest.json", "summary.json", "common_initial_state.npz"):
        (tmp_path/name).write_bytes((PREPARED/name).read_bytes())
    (tmp_path/artifact).write_bytes((tmp_path/artifact).read_bytes()+b"corrupt")
    with pytest.raises(ValueError, match="artifact hash mismatch"):
        prep.load_frozen_prepared_state(tmp_path)


@pytest.fixture(scope="module")
def dyadic_projection_checks():
    # Independent finite polynomial fixture, not a refit of the common state.
    source = np.array([.125, -.25, -.0625, .28125, -.0625, -.03125])
    return source, {policy: prep.project_saved_legendre(source, 8, length=.75, policy=policy, dps=40)
                    for policy in (prep.UNCONSTRAINED_PROJECTION, prep.CONSTRAINED_PROJECTION)}


@pytest.mark.parametrize("policy", (prep.UNCONSTRAINED_PROJECTION, prep.CONSTRAINED_PROJECTION))
def test_exact_essential_zero_polynomial_is_reproduced_in_full_shen_space(dyadic_projection_checks, policy):
    source, records = dyadic_projection_checks
    projected = records[policy]
    assert projected["raw"].shape == (7,)  # all p-1 spatial functions retained
    assert projected["legendre"].shape == (9,)
    np.testing.assert_allclose(projected["legendre"][:len(source)], source, atol=2e-16, rtol=0)
    np.testing.assert_allclose(projected["legendre"][len(source):], 0., atol=2e-38, rtol=0)
    assert np.max(abs(projected["mp_endpoint_error"])) < 2e-36
    np.testing.assert_array_equal(projected["source_jets"][0], 0.)


@pytest.mark.parametrize("derivative", (0, 1, 2))
def test_analytic_projection_jets_match_independent_legendre_derivatives(dyadic_projection_checks, derivative):
    source, records = dyadic_projection_checks
    length = .75
    independently = leg.legval([-1., 1.], leg.legder(source, m=derivative))*(2/length)**derivative
    for projected in records.values():
        np.testing.assert_allclose(projected["source_jets"][derivative], independently,
                                   rtol=2e-15, atol=2e-15)


def test_initial_projection_endpoint_rule_preserves_full_target_tail_not_filtering():
    source = np.zeros(15)
    source[:6] = [.125, -.25, -.0625, .28125, -.0625, -.03125]
    source[12], source[14] = 1e-9, -1e-9
    saved = source.copy()
    ordinary = prep.project_saved_legendre(source, 8, policy=prep.UNCONSTRAINED_PROJECTION, dps=40)
    common = prep.project_saved_legendre(source, 8, policy=prep.CONSTRAINED_PROJECTION, dps=40)
    np.testing.assert_array_equal(source, saved)
    np.testing.assert_array_equal(common["source_jets"], ordinary["source_jets"])
    # Source tail is beyond p=8, yet its endpoint derivative remains in target.
    assert np.max(abs(ordinary["mp_endpoint_error"][2])) > 1e-6
    assert np.max(abs(common["mp_endpoint_error"][1:])) < 2e-34
    assert common["raw"].shape == ordinary["raw"].shape == (7,)
    assert np.linalg.norm(common["raw"]-ordinary["raw"]) > 0


def test_constraint_applies_only_to_initial_representation_not_future_test_space(dyadic_projection_checks):
    _, records = dyadic_projection_checks
    common = records[prep.CONSTRAINED_PROJECTION]
    future = common["legendre"].copy()
    # This is one of the original essential-zero Shen test functions B_6.
    # Its endpoint first/second derivatives are not excluded from later motion.
    future[6] += 1e-5
    future[8] -= 1e-5
    np.testing.assert_allclose(leg.legval([-1., 1.], future), 0., rtol=0, atol=2e-16)
    original = leg.legval([-1., 1.], leg.legder(common["legendre"]))
    assert np.linalg.norm(leg.legval([-1., 1.], leg.legder(future))-original) > 1e-5


def test_representation_does_not_erase_saved_essential_value_roundoff():
    source = np.array([.125, -.25, -.0625, .28125, -.0625, -.03125])
    source[0] += 2**-40
    result = prep.project_saved_legendre(source, 8, policy=prep.CONSTRAINED_PROJECTION, dps=40)
    assert np.all(result["source_jets"][0] != 0)
    np.testing.assert_allclose(result["projected_jets"][0], 0., atol=2e-38, rtol=0)
    np.testing.assert_allclose(result["mp_endpoint_error"][0], -result["source_jets"][0],
                               rtol=0, atol=2e-38)
    assert np.max(abs(result["mp_endpoint_error"][1:])) < 2e-35


@pytest.mark.parametrize("invalid", (np.array([np.nan]), np.array([np.inf]), np.zeros((2, 3))))
def test_nonfinite_or_nonscalar_source_projection_is_rejected(invalid):
    with pytest.raises(ValueError, match="Finite 1D"):
        prep.project_saved_legendre(invalid, 8, dps=40)


def test_projection_requires_explicit_recognized_policy():
    with pytest.raises(ValueError, match="Unknown projection policy"):
        prep.project_saved_legendre([1., 0., -1.], 8, policy="filter_tail", dps=40)
    with pytest.raises(ValueError, match="p>=5"):
        prep.project_saved_legendre([1., 0., -1.], 4, policy=prep.CONSTRAINED_PROJECTION, dps=40)


def test_preserved_theta3_enters_common_source_directly_without_per_p_reconstruction(frozen_prepared_state):
    state, _, _ = frozen_prepared_state
    x = np.linspace(0., state.length, 53)
    eps = .05
    for derivative in (0, 1, 2):
        actual = state.evaluate(x, eps, derivative, require_admitted=False)
        background_part = eps*state.background.evaluate(x, derivative)[:, 1]
        expected = background_part+eps**3*state.correction.evaluate(x, derivative)
        np.testing.assert_array_equal(actual[:, 2], expected)
    assert state.admitted is False
    with pytest.raises(RuntimeError, match="not admitted"):
        state.evaluate(x, eps)


@pytest.fixture(scope="module")
def historical_runner_ast():
    # Read-only Git object access; no checkout and no working-tree mutation.
    source = subprocess.run(
        ["git", "show", BASELINE_HEAD+":scripts/analysis/simulate_weakly_nonlinear_planar_rod.py"],
        cwd=ROOT, check=True, stdout=subprocess.PIPE).stdout.decode("utf8")
    return ast.parse(source)


def test_optional_initializer_keeps_all_other_old_runner_functions(historical_runner_ast):
    current = ast.parse((ROOT/"scripts/analysis/simulate_weakly_nonlinear_planar_rod.py").read_text(encoding="utf8"))
    before = {x.name:x for x in historical_runner_ast.body if isinstance(x, ast.FunctionDef)}
    after = {x.name:x for x in current.body if isinstance(x, ast.FunctionDef)}
    assert before.keys() == after.keys()
    for name in before.keys()-{"integrate_case"}:
        assert ast.dump(before[name], include_attributes=False) == ast.dump(after[name], include_attributes=False), name
    manifest = _read(PREPARED/"manifest.json")
    assert manifest["identity"]["code_hashes"]["scripts/analysis/simulate_weakly_nonlinear_planar_rod.py"] == (
        "333beed99948d8336dc4bdb5683d9991de684f9f1ba53ed9a1484a523b453759")


def test_default_initializer_preserves_old_operations_and_solver_after_failure_guards(historical_runner_ast):
    current = ast.parse((ROOT/"scripts/analysis/simulate_weakly_nonlinear_planar_rod.py").read_text(encoding="utf8"))
    old = next(x for x in historical_runner_ast.body if isinstance(x, ast.FunctionDef) and x.name=="integrate_case")
    new = next(x for x in current.body if isinstance(x, ast.FunctionDef) and x.name=="integrate_case")
    assert [x.arg for x in new.args.kwonlyargs] == ["initial_coordinates", "history_buffer"]
    assert all(ast.dump(value)==ast.dump(ast.Constant(None)) for value in new.args.kw_defaults)
    # Normalization erases only the explicitly authorized q0 alternate path
    # and failure-preservation guards. Every prior mathematical operation,
    # tolerance expression, solver argument and accepted-step expression stays.
    branch = next(x for x in new.body if isinstance(x, ast.If))
    assert ast.dump(branch.test) == ast.dump(ast.parse("initial_coordinates is None").body[0].value)
    assert ast.dump(branch.body[0]) == ast.dump(old.body[1])
    branch_index = new.body.index(branch)
    new.body[branch_index:branch_index+1] = branch.body
    allocation = next(x for x in new.body if isinstance(x, ast.If)
                      and ast.dump(x.test)==ast.dump(ast.parse("history_buffer is None").body[0].value))
    old_allocation = next(x for x in old.body if isinstance(x, ast.Assign)
                          and any(isinstance(t,ast.Name) and t.id=="history" for t in x.targets))
    assert ast.dump(allocation.body[0])==ast.dump(old_allocation)
    index = new.body.index(allocation)
    new.body[index:index+1] = allocation.body  # caller buffer is storage only
    rhs = next(x for x in new.body if isinstance(x, ast.FunctionDef) and x.name=="rhs")
    guard = rhs.body.pop(0)
    assert isinstance(guard, ast.If)
    assert ast.dump(guard.test) == ast.dump(ast.parse("not np.all(np.isfinite(y))").body[0].value)
    assert ast.dump(guard.body[0]) == ast.dump(ast.parse('raise ArithmeticError("NONFINITE_STATE")').body[0])
    loop = next(x for x in new.body if isinstance(x, ast.While))
    controlled = next(x for x in loop.body if isinstance(x, ast.Try))
    assert not controlled.orelse and not controlled.finalbody
    assert len(controlled.handlers) == 1 and isinstance(controlled.handlers[0].body[-1], ast.Break)
    old_loop = next(x for x in old.body if isinstance(x, ast.While))
    assert ast.dump(controlled.body[0]) == ast.dump(old_loop.body[2])
    at = loop.body.index(controlled)
    loop.body[at:at+1] = controlled.body
    stats = next(x.value for x in new.body if isinstance(x, ast.Assign)
                 and any(isinstance(t, ast.Name) and t.id=="stats" for t in x.targets))
    extra = next(i for i,key in enumerate(stats.keys)
                 if isinstance(key, ast.Constant) and key.value=="internal_time_steps")
    assert ast.dump(stats.values[extra]) == ast.dump(ast.Name(id="steps",ctx=ast.Load()))
    del stats.keys[extra]; del stats.values[extra]  # additional saved metadata only
    new.args = old.args
    assert ast.dump(new, include_attributes=False) == ast.dump(old, include_attributes=False)


@pytest.fixture
def mocked_initializer(monkeypatch):
    from types import SimpleNamespace
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    observed = {"solver_calls":0, "projected":[], "shape_calls":0, "reset_calls":0}
    def forbidden(*args, **kwargs):
        raise AssertionError("An initializer-forwarding test must not step, sample or evaluate an ODE")
    def project(values):
        observed["projected"].append(np.asarray(values).copy())
        return np.array([.11, .22, .33, .44])
    disc = SimpleNamespace(ndof=4, p=2, nq=5, x=np.array([.25, .75]),
        project=project, rhs=forbidden, jacobian=forbidden,
        reset_counters=lambda:observed.update(reset_calls=observed["reset_calls"]+1),
        counters=lambda:{"rhs":0,"jacobian":0})
    settings = {"rtol":2e-10, "atol":np.ones(8)*1e-13,
                "max_step":.01, "coordinate_scales":[1.,1.,1.,1.],
                "velocity_scale_multiplier":1.}
    def create_finished(fun,t0,y0,t_end,**kwargs):
        assert t0 == t_end == 0.
        observed["solver_calls"] += 1
        observed["y0"] = y0.copy()
        observed["settings"] = kwargs
        return SimpleNamespace(status="finished",nfev=0,njev=0,nlu=0,
                               step=forbidden,dense_output=forbidden)
    def shape(points):
        observed["shape_calls"] += 1
        return np.ones((len(points),4))*3
    monkeypatch.setattr(runner,"np",np,raising=False)
    monkeypatch.setattr(runner,"Radau",create_finished,raising=False)
    monkeypatch.setattr(runner,"time_settings",lambda *args:settings)
    return runner, disc, shape, observed, settings


def _mock_initializer_call(runner,disc,shape,**kwargs):
    """Only initialization forwarding: patched finished solver, t_end=0."""
    config = {"material_geometry":{"h":.05}, "safety":{}}
    return runner.integrate_case(disc,shape,{},config,.05,"tight",
                                 np.array([0.]),float("inf"),**kwargs)


def test_old_initialization_default_still_projects_exact_old_fields_once(mocked_initializer):
    runner,disc,shape,observed,settings = mocked_initializer
    history,stats = _mock_initializer_call(runner,disc,shape)
    assert observed["shape_calls"] == len(observed["projected"]) == observed["solver_calls"] == 1
    np.testing.assert_array_equal(observed["projected"][0], .05*.05*np.ones((2,4))*3)
    np.testing.assert_array_equal(history[0,:4], [.11,.22,.33,.44])
    np.testing.assert_array_equal(history[0,4:], 0.)
    assert stats["accepted_internal_steps"] == stats["nfev"] == stats["njev"] == stats["nlu"] == 0
    np.testing.assert_array_equal(observed["settings"]["atol"],settings["atol"])


def test_explicit_whitened_q0_is_copied_without_reprojection_or_readmission(mocked_initializer):
    runner,disc,_,observed,_ = mocked_initializer
    def forbidden(*args, **kwargs):
        raise AssertionError("Already prepared full coefficients must not be reconstructed/reprojected")
    disc.project = forbidden
    q0 = np.array([-.01,.02,.03,-.04])
    history,stats = _mock_initializer_call(runner,disc,forbidden,initial_coordinates=q0)
    np.testing.assert_array_equal(history[0,:4],q0)
    np.testing.assert_array_equal(history[0,4:],0.)
    q0[:] = 100.
    np.testing.assert_array_equal(history[0,:4],[-.01,.02,.03,-.04])
    np.testing.assert_array_equal(observed["y0"][:4],history[0,:4])
    assert not observed["projected"]
    assert stats["time_end"] == stats["target_time_end"] == 0.
    assert stats["accepted_internal_steps"] == 0


@pytest.mark.parametrize("invalid", (
    np.zeros(3), np.zeros((1,4)), np.array([0.,0.,np.nan,0.]),np.array([0.,0.,np.inf,0.])))
def test_explicit_initial_coordinates_reject_bad_vector_before_solver(mocked_initializer,invalid):
    runner,disc,shape,observed,_ = mocked_initializer
    with pytest.raises(ValueError,match="finite full-sized"):
        _mock_initializer_call(runner,disc,shape,initial_coordinates=invalid)
    assert observed["solver_calls"] == observed["shape_calls"] == len(observed["projected"]) == 0

@pytest.fixture(scope="module")
def saved_precision_evidence():
    config_path = ROOT/"data/input/planar_prepared_feasibility.json"
    if not config_path.is_file():
        pytest.skip("Bounded precision configuration absent")
    folder = ROOT/_read(config_path)["precision_evidence"]
    if not (folder/"manifest.json").is_file():
        pytest.skip("Saved precision evidence absent; it is not recomputed by tests")
    manifest = _read(folder/"manifest.json")
    for name,digest in manifest["artifact_hashes"].items():
        assert _sha(folder/name) == digest
    assert manifest["new_ODE_integrations"] == 0
    return folder


def test_saved_precision_proof_keeps_state_and_gate_fixed_without_new_bvp(saved_precision_evidence):
    proof = _read(saved_precision_evidence/"nlsp_strong_weak_precision_20261008.json")
    assert proof["input_hashes"]["source_manifest"] == _sha(PREPARED/"manifest.json")
    assert proof["policy"]["same_float64_coefficient_state"] is True
    assert proof["policy"]["whitening_matrices_held_as_exact_binary_float64"] is True
    assert proof["policy"]["original_relative_gate"] == 2e-12
    assert proof["ODE_integrations"] == proof["eigensolves"] == 0
    for p in (48,64):
        for state in ("initial_zero_velocity","previous_synthetic_velocity"):
            rows = [row for row in proof["rows"] if row["p"]==p and row["state"]==state]
            by_stage = {(row["stage"],str(row["precision"])):row for row in rows}
            for stage,precision in (("frozen_runtime","float64"),
                                    ("stored_arrays","45"),
                                    ("rebuilt_basis_float_gauss","45")):
                assert float(by_stage[(stage,precision)]["relative_residual"]) > 2e-12
            r45 = float(by_stage[("refined_gauss","45")]["relative_residual"])
            r70 = float(by_stage[("refined_gauss","70")]["relative_residual"])
            assert r70 < r45 < 2e-12
    # This is arithmetic localization, not acceptance of the old failed DP gate.
    old = _read(PREPARED/"auxiliary_weak_identity.json")
    assert all(not row["relative_pass"] for row in old["rows"] if row["p"] in (48,64))


@pytest.mark.parametrize("p", (48,64))
def test_saved_initial_representation_precision_reproducible_with_one_common_target(saved_precision_evidence,p):
    proof = _read(saved_precision_evidence/"prepared_precision_projection_20261008.json")
    assert proof["provenance"]["manifest_sha256"] == _sha(PREPARED/"manifest.json")
    assert proof["provenance"]["regenerated_profiles"] is False
    assert proof["provenance"]["regenerated_Theta3"] is False
    options = proof["cases"][str(p)]
    ordinary,common = options[prep.UNCONSTRAINED_PROJECTION],options[prep.CONSTRAINED_PROJECTION]
    np.testing.assert_array_equal(ordinary["source_jets"],common["source_jets"])
    for variant in (ordinary,common):
        assert variant["new_ODE_integrations"] == variant["new_eigendecompositions"] == 0
        assert variant["same_trial_test_space"] is True
        assert variant["40_70_raw_relative_difference"] < 2e-12
        assert variant["40_70_endpoint_absolute_difference"] < 2e-12
    assert all(max(row["relative_L2"],row["relative_max"],row["endpoint_fixed_scaled_error"]) <= 1e-6
               for row in common["rows"])
    assert len(common["rows"]) == 12
    assert {row["field"] for row in common["rows"]} == {"u","w","theta","c"}
    assert np.max(abs(np.asarray(common["source_essential_values"]))) > 0.


def test_saved_precision_mechanism_distinguishes_arithmetic_and_source_roundoff(saved_precision_evidence):
    proof = _read(saved_precision_evidence/"prepared_precision_synthetic_mechanism_20261008.json")
    for p in ("48","64"):
        ordinary = proof["synthetic"][p][prep.UNCONSTRAINED_PROJECTION]
        common = proof["synthetic"][p][prep.CONSTRAINED_PROJECTION]
        assert ordinary["quintic_endpoint_error"] < 2e-12
        assert common["quintic_endpoint_error"] < 2e-12
        assert ordinary["perturbed_second_endpoint_unresolved"] > 1e-6
        assert common["perturbed_second_endpoint_unresolved"] < 2e-12
        assert ordinary["perturbed_low_profile_change"] == 0.
        assert common["perturbed_low_profile_change"] > 0.  # constrained L2, no tail filter
    for row in proof["essential_roundoff_mechanism"]:
        assert np.max(abs(np.asarray(row["original_source_essential_values"]))) > 0.
        assert "diagnostic subtraction only" in row["qualification"]
        assert "actual frozen state" in row["qualification"]


@pytest.mark.parametrize("failed", (False,True))
def test_feasibility_authorization_retains_failed_strict_rows_and_does_not_admit_state(frozen_prepared_state,failed):
    from scripts.analysis import prepare_planar_initial_state as cli
    state,_,_ = frozen_prepared_state
    strict = [{"p":48,"check":"projection","pass":True,"tolerance":1e-6},
              {"p":64,"check":"strong_weak","pass":not failed,"tolerance":2e-12}]
    preserved = json.loads(json.dumps(strict))
    mode = cli.authorize_feasibility(strict,{"positive_mass":True,"finite_RHS":True},True)
    assert mode == ("EXPLORATORY_NOT_CERTIFIED" if failed else "STRICT_ADMITTED")
    assert strict == preserved
    assert state.admitted is False


@pytest.mark.parametrize("basic,evidence", (
    ({"finite_RHS":False,"positive_mass":True},True),
    ({"finite_RHS":True,"positive_mass":False},True),
    ({"finite_RHS":True,"positive_mass":True},False),
))
def test_exploratory_authorization_cannot_cover_unexplained_basic_inconsistency(basic,evidence):
    from scripts.analysis import prepare_planar_initial_state as cli
    with pytest.raises(ArithmeticError,match="BLOCKED_BY_UNEXPLAINED_INCONSISTENCY"):
        cli.authorize_feasibility([{"pass":False}],basic,evidence)


def test_feasibility_configuration_is_bounded_and_preserves_thresholds():
    from scripts.analysis import prepare_planar_initial_state as cli
    config = _read(cli.FEASIBILITY_CONFIG)
    assert config["prepared_bundle"] == str(PREPARED.relative_to(ROOT)).replace(chr(92),"/")
    assert config["degrees"] == [48,64]
    assert config["cases"] == [[48,"tight"],[64,"tight"],[64,"allowed_extra"]]
    assert config["short_periods"] == .1
    assert config["amplitude_over_h"] == .05
    assert config["projection_policy"] == prep.CONSTRAINED_PROJECTION
    assert config["projection_dps"] == [40,70]
    assert config["strict_identity_relative"] == 2e-12
    assert config["profile_policy"]["relative_tolerance"] == config["profile_policy"]["endpoint_tolerance"] == 1e-6
    assert config["budget"]["numerical_wall_seconds"] == 900
    assert config["budget"]["local_precision_seconds"] == 180
    assert config["budget"]["maximum_integrations"] == 3
    assert config["semantics"]["initial_rule_is_not_boundary_condition"] is True
    assert config["semantics"]["initial_endpoint_derivatives_not_dynamic_constraints"] is True
    assert config["semantics"]["all_coordinates_retained"] is True


@pytest.mark.parametrize("changed", ("projection","precision","time_level","source"))
def test_feasibility_identity_invalidates_each_representation_time_or_physical_target_change(tmp_path,changed):
    from scripts.analysis import prepare_planar_initial_state as cli
    config = _read(cli.FEASIBILITY_CONFIG)
    source = tmp_path/"config.json"
    cli.write_json(source,config)
    before,_ = cli.feasibility_identity(source)
    if changed=="projection":config["projection_policy"]=prep.UNCONSTRAINED_PROJECTION
    elif changed=="precision":config["projection_dps"]=[45,70]
    elif changed=="time_level":config["cases"][-1][1]="tight"
    else:
        replacement = tmp_path/"other_frozen_target"
        replacement.mkdir()
        manifest = _read(PREPARED/"manifest.json")
        manifest["test_fixture_changed_target"] = True
        cli.write_json(replacement/"manifest.json",manifest)
        config["prepared_bundle"]=str(replacement)
    cli.write_json(source,config)
    after,item = cli.feasibility_identity(source)
    assert before != after
    assert item["config"] == config


def _synthetic_feasibility_cache(folder,item):
    from scripts.analysis import prepare_planar_initial_state as cli
    folder.mkdir(parents=True)
    rows=[{"field":field,"derivative":d,"relative_L2":1e-8,"relative_max":1e-8,
           "endpoint_fixed_scaled_error":1e-8}
          for d in (0,1,2) for field in ("u","w","theta","c")]
    summary={"schema":"nlsp-prepared-feasibility-v1","synthetic_fixture":True,
        "statuses":{"NLSP_STRICT_INITIAL_VERIFICATION":"PARTIAL",
                    "NLSP_PREPARED_FEASIBILITY_RUN":"COMPLETED_EXPLORATORY_NOT_CERTIFIED"},
        "new_ODE_integrations":3,"new_eigendecompositions":0,
        "execution_mode":"EXPLORATORY_NOT_CERTIFIED","state_admitted_flag":False,
        "original_projection":{str(p):{"rows":rows} for p in (48,64)},
        "projection":{str(p):{"rows":rows} for p in (48,64)}}
    cli.write_json(folder/"summary.json",summary)
    cli.write_json(folder/"manifest.json",cli.manifest_for(folder,item))
    return summary


@pytest.mark.parametrize("action", ("compute","report-only","plot-only"))
def test_cached_feasibility_entrypoints_never_repeat_preparation_precision_or_three_controls(
        tmp_path,monkeypatch,capsys,action):
    from scripts.analysis import prepare_planar_initial_state as cli
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
    item={"synthetic_cache":"feasibility-v1"}
    bundle=tmp_path/"output"/"fixture"
    expected=_synthetic_feasibility_cache(bundle,item)
    def forbidden(*args,**kwargs):
        raise AssertionError("Cached feasibility may only read saved data")
    for name in ("run_feasibility","run_compute","initial_coordinate_metrics","strong_weak_metrics"):
        monkeypatch.setattr(cli,name,forbidden)
    for name in ("stable_initial_projection","high_precision_source_jets","load_frozen_prepared_state"):
        monkeypatch.setattr(prep,name,forbidden)
    monkeypatch.setattr(runner,"integrate_case",forbidden)
    monkeypatch.setattr(cli,"feasibility_identity",lambda *args:("fixture",item))
    args=(["--compute","--feasibility","--output-dir",str(bundle.parent)] if action=="compute"
          else ["--"+action,str(bundle)])
    returned=cli.main(args)
    output=json.loads(capsys.readouterr().out)
    assert returned == expected
    assert all(value==0 for value in output["this_run_counters"].values())
    assert returned["new_ODE_integrations"] == 3  # historical counts, not this run
    assert returned["state_admitted_flag"] is False


def test_feasibility_plot_repeat_preserves_cache_and_runs_zero_precision(tmp_path,monkeypatch):
    from scripts.analysis import prepare_planar_initial_state as cli
    item={"synthetic_cache":"feasibility-figures"}
    bundle=tmp_path/"bundle"
    expected=_synthetic_feasibility_cache(bundle,item)
    def forbidden(*args,**kwargs):
        raise AssertionError("Plotting cannot prepare or integrate")
    monkeypatch.setattr(cli,"run_feasibility",forbidden)
    monkeypatch.setattr(prep,"stable_initial_projection",forbidden)
    cli.plot_only(bundle)
    cli.write_json(bundle/"manifest.json",cli.manifest_for(bundle,item))
    before={p.name:_sha(p) for p in (bundle/"figures").iterdir()}
    cli.plot_only(bundle)
    assert before == {p.name:_sha(p) for p in (bundle/"figures").iterdir()}
    assert cli.validate_cache(bundle,item) == expected


FEASIBILITY_RESULT = ROOT/"results/planar_prepared_feasibility/284a4039177391d1"


@pytest.fixture(scope="module")
def completed_feasibility_evidence():
    from scripts.analysis import prepare_planar_initial_state as cli
    if not (FEASIBILITY_RESULT/"manifest.json").is_file():
        pytest.skip("Completed short controls absent; tests never reproduce them")
    summary=cli.validate_cache(FEASIBILITY_RESULT)
    runs={}
    for stats in summary["actual_short_runs"]:
        name=f'p{stats["p"]}_{stats["time_level"]}'
        folder=FEASIBILITY_RESULT/"short_controls"/name
        assert _read(folder/"case.json")==stats
        with np.load(folder/"trajectory.npz",allow_pickle=False) as data:
            runs[name]={"time":data["time"].copy(),"q0":data["q"][0].copy(),
                        "v0":data["velocity"][0].copy(),"energy":data["energy"].copy(),
                        "q_shape":data["q"].shape,"v_shape":data["velocity"].shape,
                        "stats":stats}
        with np.load(folder/"internal_steps.npz",allow_pickle=False) as data:
            runs[name]["steps"]=data["time_step"].copy()
    return summary,runs


def test_actual_feasibility_does_not_promote_failed_strict_gates_or_source_admission(completed_feasibility_evidence):
    summary,_=completed_feasibility_evidence
    assert summary["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
    assert summary["state_admitted_flag"] is False
    assert summary["independent_precision_evidence"] is True
    assert summary["numerical_policy_chosen_before_ODE"] is True
    assert all(summary["projection"][str(p)]["pass"] for p in (48,64))
    assert any(not row["pass"] for row in summary["strict_table"])
    assert summary["statuses"]["NLSP_STRICT_INITIAL_VERIFICATION"]=="PARTIAL"
    assert summary["statuses"]["NLSP_PREPARED_FEASIBILITY_RUN"]=="COMPLETED_EXPLORATORY_NOT_CERTIFIED"
    assert summary["historical_statuses"]==_read(PREPARED/"summary.json")["statuses"]
    assert summary["source_provenance"]["manifest_sha256"]==_sha(PREPARED/"manifest.json")
    assert summary["new_ODE_integrations"]==3
    assert summary["new_eigendecompositions"]==summary["new_BVP_solves"]==0


def test_actual_three_controls_have_same_real_short_interval_and_full_independent_fields(completed_feasibility_evidence):
    summary,runs=completed_feasibility_evidence
    assert set(runs)=={"p48_tight","p64_tight","p64_allowed_extra"}
    expected=None
    for name,row in runs.items():
        t,stats=row["time"],row["stats"]
        assert t[0]==0.
        assert np.all(np.diff(t)>0)
        assert t[-1]==stats["time_end"]==stats["target_time_end"]==summary["sampling"]["horizon"]
        assert len(t)==stats["samples"]==summary["sampling"]["samples"]
        assert row["q_shape"]==row["v_shape"]==(len(t),4*(stats["p"]-1))
        assert stats["status"]=="PASS"
        assert stats["execution_mode"]=="EXPLORATORY_NOT_CERTIFIED"
        assert stats["projection_policy"]==prep.CONSTRAINED_PROJECTION
        np.testing.assert_array_equal(row["v0"],0.)
        if expected is None:expected=t
        else:np.testing.assert_array_equal(t,expected)
    assert summary["sampling"]["horizon"]==.1*summary["background"]["T1"]
    assert summary["sampling"]["samples_per_bound_period"]==16
    assert summary["sampling"]["output_grid_not_time_accuracy_control"] is True
    with np.load(RECOVERY/"controls/new_p32/trajectory.npz") as old:
        old_times=old["time"][old["time"]<=expected[-1]]
    assert np.all(np.isin(old_times,expected))


def test_actual_runner_uses_exact_projected_q0_and_same_p64_q0_at_both_time_levels(completed_feasibility_evidence):
    _,runs=completed_feasibility_evidence
    for name,row in runs.items():
        p=row["stats"]["p"]
        with np.load(FEASIBILITY_RESULT/"initial_projection"/f"p{p}.npz",allow_pickle=False) as projected:
            np.testing.assert_array_equal(row["q0"],projected["q"])
            np.testing.assert_array_equal(row["v0"],projected["velocity"])
        n=p-1
        assert np.linalg.norm(row["q0"][:n])>0.
        assert np.linalg.norm(row["q0"][3*n:])>0.
    np.testing.assert_array_equal(runs["p64_tight"]["q0"],runs["p64_allowed_extra"]["q0"])


def test_actual_saved_steps_preserve_budget_step_statistics_without_invented_rejections(completed_feasibility_evidence):
    summary,runs=completed_feasibility_evidence
    for row in runs.values():
        steps,stats=row["steps"],row["stats"]
        np.testing.assert_array_equal(steps,stats["internal_time_steps"])
        assert len(steps)==stats["accepted_internal_steps"]
        assert np.all(steps>0.)
        assert np.min(steps)==stats["min_internal_step"]
        assert np.max(steps)==stats["max_internal_step"]
        assert np.max(steps)<=stats["max_step"]*(1+2e-12)
        np.testing.assert_allclose(np.sum(steps),stats["time_end"],rtol=2e-12,atol=0)
        assert "rejected_steps" not in stats
    assert summary["runtime"]["ODE_integrations"]==3
    assert summary["runtime"]["numerical_wall_seconds"]<=900.
    assert summary["runtime"]["local_precision_seconds"]<=180.


def test_actual_saved_mass_safety_and_energy_gates_pass_for_each_short_run(completed_feasibility_evidence,pilot_config):
    _,runs=completed_feasibility_evidence
    for row in runs.values():
        stats,energy=row["stats"],row["energy"]
        assert np.all(np.isfinite(energy)) and np.all(energy>0.)
        drift=float(np.max(abs((energy-energy[0])/energy[0])))
        assert drift==stats["relative_energy_drift_max"]
        assert drift<=pilot_config["gates"]["energy_relative_drift"]
        assert stats["mass_lower_bound_min"]>=pilot_config["safety"]["min_relative_mass_eigenvalue"]
        assert stats["mass_condition_bound_max"]>=1.
        for key,value in stats["safety_extrema"].items():
            if key=="min_one_plus_c":assert value>pilot_config["safety"][key]
            else:assert value<=pilot_config["safety"][key]


@pytest.mark.parametrize("comparison", ("new_spatial_comparison","new_temporal_comparison"))
def test_actual_comparisons_report_all_eight_unchanged_norms_without_phase_or_floor_redefinition(
        completed_feasibility_evidence,pilot_config,comparison):
    summary,_=completed_feasibility_evidence
    report=summary[comparison]
    assert set(report["fields"])=={part+"_"+field for part in ("q","velocity") for field in ("u","w","theta","c")}
    for name,row in report["fields"].items():
        field=name.split("_")[-1]
        tolerance=pilot_config["gates"]["w_theta_relative" if field in ("w","theta") else "u_c_relative"]
        assert row["tolerance"]==tolerance
        assert row["relative_L2"]==row["max_time_L2_difference"]/max(row["reference_max_time_L2"],row["numerical_floor"])
        assert row["relative_max"]==row["max_space_time_difference"]/max(row["reference_max_space_time"],row["numerical_floor"])
        assert row["pass"]==(row["relative_L2"]<=tolerance and row["relative_max"]<=tolerance)
        assert row["fixed_scaled_L2"]==row["max_time_L2_difference"]/(row["fixed_physical_scale"]*np.sqrt(summary["background"]["L"]))
        assert row["fixed_scaled_max"]==row["max_space_time_difference"]/row["fixed_physical_scale"]
    assert "no phase alignment" in report["sampling_qualification"]


def test_actual_spatial_failure_is_not_hidden_by_temporal_pass_or_bending_displacement(completed_feasibility_evidence):
    summary,_=completed_feasibility_evidence
    spatial,temporal=summary["new_spatial_comparison"],summary["new_temporal_comparison"]
    assert temporal["status"]=="PASS" and all(row["pass"] for row in temporal["fields"].values())
    assert spatial["status"]=="PARTIAL"
    assert [name for name,row in spatial["fields"].items() if not row["pass"]]==["velocity_theta"]
    assert spatial["fields"]["velocity_theta"]["relative_max"]>1e-4
    assert spatial["fields"]["velocity_theta"]["relative_L2"]<=1e-4
    assert summary["statuses"]["NLSP_PREPARED_SHORT_SPATIAL_CHECK"]=="PARTIAL"
    assert summary["statuses"]["NLSP_PREPARED_SHORT_TEMPORAL_CHECK"]=="PASS"
    assert "continuous-PDE convergence" in summary["qualification"]
