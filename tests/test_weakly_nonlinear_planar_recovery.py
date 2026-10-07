"""Bounded planar recovery controls; these tests never integrate an ODE.

Historical trajectories and the historical helper are read-only evidence.
No test modifies the initial shape, boundary data, action or output bundle.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
import pytest

from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import weakly_nonlinear_planar_dynamics as planar
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section


ROOT = Path(__file__).resolve().parents[1]
HISTORICAL = ROOT / "results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef"
BASELINE_SHA256 = "fa03369ea8aed8f6475679c319b4026eec169d7cad5ea39237f07354d6f37e90"
BASELINE_HEAD = "1510d75c106a28a4da899c7eea1a337f11791ce6"
FIELDS = ("u", "w", "theta", "c")
ROUND_OFF_GATE = 2e-12


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _scaled_allclose(actual, reference, scale=None):
    reference = np.asarray(reference)
    scale = max(float(np.max(np.abs(reference), initial=0)), 1e-30) if scale is None else scale
    np.testing.assert_allclose(actual, reference, rtol=ROUND_OFF_GATE,
                               atol=ROUND_OFF_GATE * scale)


@pytest.fixture(scope="module")
def model():
    return rod.derive_polynomials()


@pytest.fixture(scope="module")
def historical():
    if not HISTORICAL.is_dir():
        pytest.skip("Read-only historical pilot bundle unavailable")
    return json.loads((HISTORICAL / "summary.json").read_text(encoding="utf8"))


@pytest.fixture(scope="module")
def coefficients(historical):
    return rod.RodCoefficients(**historical["coefficients"])


@pytest.fixture(scope="module")
def small_disc(coefficients, model):
    return planar.PlanarGalerkin(coefficients, p=6, model=model)


def _physical_state(disc):
    def fields(points):
        bubble = points * (1 - points)
        return np.column_stack((2e-4 * bubble * (2 * points - 1),
                                3e-3 * bubble, 4e-3 * bubble * (2 * points - 1),
                                3e-4 * bubble))

    def velocities(points):
        bubble = points * (1 - points)
        return np.column_stack((3e-4 * bubble, 2e-3 * bubble,
                                3e-3 * bubble * (1 + points),
                                -2e-4 * bubble * (1 - .2 * points)))

    return disc.project(fields), disc.project(velocities)


def test_quadratic_axial_flux_is_derived_from_protected_quartic_action(model):
    zero_inactive = {name + suffix: 0 for name in ("v", "Phi", "psi")
                     for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")}
    potential = model.V4.substitute(zero_inactive)
    actual = potential.derivative("u_s").homogeneous(2).substitute({"u_s": 0, "c": 0})
    p = model.symbols
    expected = ((p["C"] - p["S"]) * p["theta"] * p["w_s"]
                + (p["S"] - p["C"] / 2) * p["theta"] ** 2)
    assert actual == expected
    endpoint = actual.total_derivative("s").substitute({"theta": 0})
    assert endpoint == (p["C"] - p["S"]) * p["theta_s"] * p["w_s"]


def test_all_four_initial_endpoint_residuals_before_linear_eigenpair_identities(model):
    zero = {name + suffix: 0 for name in rod.FIELD_ORDER
            for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")
            if not (name in ("w", "theta") and suffix in ("_s", "_ss"))}
    residual = [model.residual_a[index].substitute(zero) for index in (0, 1, 5, 6)]
    p = model.symbols
    assert residual[0] == -(p["C"] - p["S"]) * p["theta_s"] * p["w_s"]
    assert residual[1] == p["S"] * (p["theta_s"] - p["w_ss"])
    assert residual[2] == -p["S"] * p["w_s"] - p["Bp"] * p["theta_ss"]
    assert residual[3] == rod.Polynomial()
    assert residual[0].homogeneous(1) == rod.Polynomial()
    assert residual[0].homogeneous(2) == residual[0]


def test_continuous_initial_trace_mismatch_is_quadratic_not_a_galerkin_trace(historical, coefficients):
    source = mh.project_jang_reduced_rectangular(rectangular_section(
        E=1., nu=.3, rho=1., width=.2, thickness=.05, K=5/6))
    omega = historical["initial_eigenpair"]["omega"]
    modes = np.asarray(historical["initial_eigenpair"]["analytic_coefficients"])
    endpoints = mh.finite_state_basis(source, 1., omega, [0., 1.], "timoshenko") @ modes
    first = mh.finite_state_basis(source, 1., omega, [0., 1.], "timoshenko", 1) @ modes
    second = first @ mh.harmonic_state_matrix(source, omega, "timoshenko").T
    assert np.max(np.abs(endpoints[:, :2])) < 2e-12
    amplitude = .0025
    trace = (coefficients.C - coefficients.S) / coefficients.m * amplitude**2 * first[:, 0] * first[:, 1]
    half = (coefficients.C - coefficients.S) / coefficients.m * (amplitude / 2)**2 * first[:, 0] * first[:, 1]
    assert trace[0] > 1e-5 and trace[1] < -1e-5
    np.testing.assert_allclose(trace, 4 * half, rtol=1e-15, atol=0)
    # Source eigenpair's linear acceleration is zero at the essential ends.
    linear_w = coefficients.S / coefficients.m * (second[:, 0] - first[:, 1])
    linear_theta = (coefficients.Bp * second[:, 1]
                    + coefficients.S * (first[:, 0] - endpoints[:, 1])) / coefficients.jp
    assert np.max(np.abs(amplitude * linear_w)) < 1e-13
    assert np.max(np.abs(amplitude * linear_theta)) < 1e-12


@pytest.mark.parametrize("field", range(4))
@pytest.mark.parametrize("part", ("position", "velocity"))
def test_physical_l2_projection_has_orthogonal_tail_and_common_motion(field, part):
    xi, weights = leggauss(100)
    weights = weights / 2
    vandermonde = legvander(xi, 32)
    low = vandermonde[:, :23] - vandermonde[:, 2:25]
    high = vandermonde[:, :31] - vandermonde[:, 2:33]
    scale = (1e-6, 2e-3, 7e-3, 2e-5)[field] * (3 if part == "velocity" else 1)
    high_coefficients = scale * np.cos(np.arange(31) + .31 * field) / (1 + np.arange(31))**2
    low_coefficients = high_coefficients[:23].copy()
    low_coefficients[3] += scale * .02
    high_profile, low_profile = high @ high_coefficients, low @ low_coefficients
    gram = low.T @ (weights[:, None] * low)
    projection = low @ np.linalg.solve(gram, low.T @ (weights * high_profile))
    tail, common = high_profile - projection, projection - low_profile
    difference = high_profile - low_profile
    assert np.max(np.abs(low.T @ (weights * tail))) < ROUND_OFF_GATE * scale
    _scaled_allclose(difference, tail + common, scale)
    lhs = float(weights @ difference**2)
    rhs = float(weights @ tail**2 + weights @ common**2)
    assert abs(lhs - rhs) <= ROUND_OFF_GATE * max(lhs, rhs)
    assert float(weights @ tail**2) > 0
    assert float(weights @ common**2) > 0


def test_energy_gradient_hessian_cache_request_order_and_state_change(small_disc, monkeypatch):
    disc = small_disc
    q, velocity = _physical_state(disc)
    calls = {"energy": 0, "gradient": 0, "hessian": 0}
    for label, compiler in (("energy", disc._potential_energy),
                            ("gradient", disc._potential_gradient),
                            ("hessian", disc._potential_hessian)):
        original = compiler.evaluate

        def wrapped(variables, original=original, label=label):
            calls[label] += 1
            return original(variables)

        monkeypatch.setattr(compiler, "evaluate", wrapped)
    energy = disc.potential(q, gradient=False)["V"]
    assert calls == {"energy": 1, "gradient": 0, "hessian": 0}
    gradient = disc.potential(q)["gradient"]
    assert calls == {"energy": 1, "gradient": 1, "hessian": 0}
    complete = disc.potential(q, hessian=True)
    assert calls == {"energy": 1, "gradient": 1, "hessian": 1}
    _scaled_allclose(complete["V"], energy)
    np.testing.assert_array_equal(complete["gradient"], gradient)
    np.testing.assert_array_equal(disc.potential(q, hessian=True)["hessian"], complete["hessian"])
    disc.rhs(0., np.r_[q, velocity])
    assert calls == {"energy": 1, "gradient": 1, "hessian": 1}
    changed = q.copy()
    changed[disc.slices["u"].start] += 1e-7
    changed_gradient = disc.potential(changed)["gradient"]
    assert calls == {"energy": 1, "gradient": 2, "hessian": 1}
    assert np.linalg.norm(changed_gradient - gradient) > 0
    disc.potential(changed, hessian=True)
    assert calls["hessian"] == 2


def test_velocity_change_reuses_only_valid_mass_and_changes_inertial_force(small_disc):
    disc = small_disc
    q, velocity = _physical_state(disc)
    disc.reset_counters()
    first_mass = disc.mass_matrix(q)
    count = disc.mass_factorizations
    first_force = disc.inertial_terms(q, velocity)
    first_acceleration = disc.acceleration(q, velocity)
    altered_velocity = velocity.copy()
    altered_velocity[disc.slices["theta"]] *= 1.7
    altered_velocity[disc.slices["c"]] *= .6
    altered_force = disc.inertial_terms(q, altered_velocity)
    altered_acceleration = disc.acceleration(q, altered_velocity)
    np.testing.assert_array_equal(first_mass, disc.mass_matrix(q))
    assert disc.mass_factorizations == count
    assert np.linalg.norm(altered_force - first_force) > 0
    assert np.linalg.norm(altered_acceleration - first_acceleration) > 0
    changed_q = q.copy()
    changed_q[disc.slices["c"].start] += 1e-7
    changed_mass = disc.mass_matrix(changed_q)
    assert disc.mass_factorizations == count + 1
    assert np.linalg.norm(changed_mass - first_mass) > 0


@pytest.fixture(scope="module")
def baseline_module(tmp_path_factory):
    candidates = []
    if os.environ.get("NLSP_RECOVERY_BASELINE_HELPER"):
        candidates.append(Path(os.environ["NLSP_RECOVERY_BASELINE_HELPER"]))
    candidates.append(Path(tempfile.gettempdir()) / "coupledbeams_planar_rhs_baseline_20261007"
                      / "weakly_nonlinear_planar_dynamics_baseline.py")
    if (ROOT / "results/weakly_nonlinear_planar_recovery").is_dir():
        candidates.extend(sorted((ROOT / "results/weakly_nonlinear_planar_recovery").glob(
            "*/baseline/weakly_nonlinear_planar_dynamics.py")))
    path = next((candidate for candidate in candidates if candidate.is_file()
                 and _sha(candidate) == BASELINE_SHA256), None)
    if path is None:
        # Read the specifically authorized published baseline; no Git mutation.
        completed = subprocess.run(["git", "show", BASELINE_HEAD
                                    + ":scripts/lib/weakly_nonlinear_planar_dynamics.py"],
                                   cwd=ROOT, capture_output=True, check=True)
        assert hashlib.sha256(completed.stdout).hexdigest() == BASELINE_SHA256
        path = tmp_path_factory.mktemp("planar_baseline") / "historical_helper.py"
        path.write_bytes(completed.stdout)
    spec = importlib.util.spec_from_file_location("_planar_recovery_historical_baseline", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def real_states(historical):
    name = "p32_Aoverh0p05_tight"
    path = HISTORICAL / "cases" / name / "trajectory.npz"
    expected = historical["cases"][name]["artifact_hashes"]["trajectory.npz"]
    assert _sha(path) == expected
    with np.load(path) as archive:
        count = len(archive["time"])
        indexes = np.unique(np.minimum([0, 1, 123, 10000, 35000, 70984], count - 1))
        positions = archive["q"][indexes]
        velocities = archive["velocity"][indexes]
    assert len(indexes) >= 5
    return list(zip(positions, velocities))


def test_old_new_numerical_path_equivalence_on_distinct_real_and_small_states(
        coefficients, model, baseline_module, real_states):
    old = baseline_module.PlanarGalerkin(coefficients, p=32, model=model)
    new = planar.PlanarGalerkin(coefficients, p=32, model=model)
    states = real_states + [_physical_state(new)]
    tiny_q, tiny_v = _physical_state(new)
    states.append((tiny_q * .013, tiny_v * .021))
    for q, velocity in states:
        old_potential = old.potential(q, hessian=True)
        new_potential = new.potential(q, hessian=True)
        for label in ("V", "gradient", "hessian"):
            _scaled_allclose(new_potential[label], old_potential[label])
        _scaled_allclose(new.mass_matrix(q), old.mass_matrix(q))
        _scaled_allclose(new.inertial_terms(q, velocity), old.inertial_terms(q, velocity))
        acceleration = new.acceleration(q, velocity)
        _scaled_allclose(acceleration, old.acceleration(q, velocity))
        state = np.r_[q, velocity]
        _scaled_allclose(new.rhs(0., state), old.rhs(0., state))
        _scaled_allclose(new.jacobian(0., state), old.jacobian(0., state))
        _scaled_allclose(new.energy(q, velocity), old.energy(q, velocity))
        # Integration by parts and the six potential contributions can cancel
        # strongly for a bending eigenpair.  Scale roundoff by their actual
        # uncancelled work, not by the much smaller final summed force.
        local_gradient = new._local_potential(q)[1]
        work_scale = max(sum(np.linalg.norm(matrix.T @ (new.weights * values))
                             for matrix, values in zip(new._potential_matrices, local_gradient))
                         + np.linalg.norm(new.inertial_terms(q, velocity))
                         + np.linalg.norm(new.mass_matrix(q) @ acceleration), 1e-30)
        new_weak = new.weak_residual(q, velocity, acceleration)
        old_weak = old.weak_residual(q, velocity, acceleration)
        # The independent residual evaluator is unchanged by the optimization.
        np.testing.assert_array_equal(new_weak, old_weak)
        assert np.max(np.abs(new_weak)) <= ROUND_OFF_GATE  # Original absolute gate.
        assert np.linalg.norm(new_weak) <= ROUND_OFF_GATE * work_scale
        power_scale = max(work_scale * np.linalg.norm(velocity), 1e-30)
        _scaled_allclose(new.energy_rate(q, velocity, acceleration),
                         old.energy_rate(q, velocity, acceleration), power_scale)
        assert abs(new.energy_rate(q, velocity, acceleration)) <= ROUND_OFF_GATE * power_scale


def test_recovery_tests_do_not_call_an_integrator():
    """Guard against spending the integration budget inside this test module."""
    import ast
    tree = ast.parse(Path(__file__).read_text(encoding="utf8"))
    # Cached entrypoints are exercised below with the numerical operations
    # replaced by raising guards. Only actual integration calls are forbidden.
    forbidden = {"solve_ivp", "Radau", "integrate_case", "run_controls"}
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = node.func.id if isinstance(node.func, ast.Name) else (
                node.func.attr if isinstance(node.func, ast.Attribute) else None)
            assert name not in forbidden


@pytest.mark.parametrize("length", (1., .73))
def test_diagnostic_physical_basis_matches_action_helper(coefficients, model, length):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    disc = planar.PlanarGalerkin(coefficients, p=6, length=length, model=model)
    points = np.linspace(0, length, 29)
    independent = diagnostic.physical_bases(6, points, coefficients, length)
    action = disc.basis_at(points)
    for name, matrix in zip(FIELDS, independent):
        np.testing.assert_allclose(matrix, action[name], rtol=ROUND_OFF_GATE, atol=ROUND_OFF_GATE)
        np.testing.assert_array_equal(matrix[[0, -1]], 0.)


@pytest.mark.parametrize("field", range(4))
def test_diagnostic_projection_operator_is_physical_and_orthogonal(coefficients, field):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    xi, weights = leggauss(81)
    points, weights = (xi+1)/2, weights/2
    low = diagnostic.physical_bases(6, points, coefficients)[field]
    high = diagnostic.physical_bases(10, points, coefficients)[field]
    projector = diagnostic.projection_operator(6, 10, coefficients, nq=33)[field]
    coordinates = np.cos(np.arange(9)+.37*field)/(np.arange(9)+1)**2 * 1e-5
    reference = high @ coordinates
    projected = low @ (projector @ coordinates)
    tail = reference-projected
    assert np.max(abs(low.T @ (weights*tail))) < ROUND_OFF_GATE * max(np.linalg.norm(reference), 1e-30)
    np.testing.assert_allclose(weights @ reference**2,
                               weights @ projected**2+weights @ tail**2,
                               rtol=ROUND_OFF_GATE, atol=0)
    independently_solved = np.linalg.solve(low.T @ (weights[:, None]*low),
                                          low.T @ (weights*reference))
    _scaled_allclose(projector @ coordinates, independently_solved)


def test_diagnostic_projection_rejects_insufficient_quadrature(coefficients):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    with pytest.raises(ValueError, match="quadrature"):
        diagnostic.projection_operator(24, 48, coefficients, nq=48)


def _synthetic_historical_bundle(folder):
    """Small own-manifest fixture, independent of ignored historical results."""
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    folder.mkdir(parents=True)
    coefficients = {"m": .01, "jp": .00001}
    summary = {"coefficients": coefficients,
               "config": {"material_geometry": {"L": 1.}}}
    diagnostic.write_json(folder / "summary.json", summary)
    case = folder / "cases" / "p4_Aoverh0p05_tight"
    case.mkdir(parents=True)
    times = np.array([0., .5, 1.])
    coordinate = np.arange(36, dtype=float).reshape(3, 12)*1e-8
    indices = np.array([0, 1, 2])
    raw = np.column_stack([coordinate[:, i*3:(i+1)*3] @ transform.T
                           for i, transform in enumerate(diagnostic.physical_transforms(4, coefficients))])
    points = np.linspace(0, 1, 13)
    profiles = np.stack([raw[:, i*3:(i+1)*3] @ diagnostic.raw_basis(4, points).T
                         for i in range(4)], axis=2)
    np.savez_compressed(case / "trajectory.npz", time=times, q=coordinate,
                        snapshot_indices=indices, raw_snapshot_coefficients=raw,
                        snapshot_points=points, snapshots=profiles)
    diagnostic.write_json(case / "case.json", {
        "status": "PASS", "failure": None, "p": 4, "ndof": 12, "nq": 9,
        "amplitude_over_h": .05, "amplitude": .0025, "time_level": "tight",
        "rtol": 1e-10, "max_step": .1, "time_end": 1., "target_time_end": 1.,
        "samples": 3, "accepted_internal_steps": 10,
        "min_internal_step": .1, "max_internal_step": .1,
        "nfev": 70, "njev": 1, "nlu": 2, "integration_seconds": .1, "counters": {}})
    artifacts = {str(path.relative_to(folder)): diagnostic.sha(path)
                 for path in folder.rglob("*") if path.is_file()}
    diagnostic.write_json(folder / "manifest.json", {
        "identity": {"hashes": {"scripts/lib/weakly_nonlinear_planar_dynamics.py":
                                 "historical-code-identity-is-intentionally-different"}},
        "artifact_hashes": artifacts})
    return folder


def test_historical_manifest_validation_does_not_compare_current_code(tmp_path):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    source = _synthetic_historical_bundle(tmp_path / "historical")
    inventory = diagnostic.validate_historical(source)
    assert all(item["match"] for item in inventory["artifact_hash_validation"].values())
    assert inventory["integrations"] == inventory["symbolic_derivations"] == 0
    assert inventory["cases"]["p4_Aoverh0p05_tight"]["actual_array_time_end"] == 1.
    assert inventory["cases"]["p4_Aoverh0p05_tight"]["independent_reconstruction_snapshot_max_difference"] < 1e-18


def test_historical_manifest_validation_rejects_changed_artifact(tmp_path):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    source = _synthetic_historical_bundle(tmp_path / "historical")
    (source / "summary.json").write_text("{}", encoding="utf8")
    with pytest.raises(ValueError, match="Historical artifact hash mismatch"):
        diagnostic.validate_historical(source)


def _synthetic_recovery_bundle(folder, identity):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    folder.mkdir(parents=True)
    summary = {"statuses": {"NLSP_PLANAR_SOLVER_RECOVERY": "PARTIAL"},
               "compatibility": {"T1": 2.}, "integrations": 3,
               "budget_used_seconds": 40.}
    diagnostic.write_json(folder / "summary.json", summary)
    time_values = np.linspace(0, 10, 21)
    errors = np.tile(np.linspace(0, 1e-5, 21)[:, None], (1, 4))
    data = {"time": time_values}
    for part in ("q", "velocity"):
        data.update({part+"_difference_L2": errors,
                     part+"_tail_L2": errors*.3, part+"_common_L2": errors*.7})
        key = part+"_c"
        points = np.linspace(0, 1, 9)
        data.update({key+"_times": np.array([0., 4., 10.]), key+"_x": points,
                     key+"_difference": np.stack([points*(1-points)*scale
                                                    for scale in (0., 1e-6, 3e-6)])})
    np.savez_compressed(folder / "p24_p32_arrays.npz", **data)
    artifacts = {str(path.relative_to(folder)): diagnostic.sha(path)
                 for path in folder.iterdir() if path.is_file()}
    diagnostic.write_json(folder / "manifest.json", {"identity": identity,
                                                      "artifact_hashes": artifacts})
    return folder


@pytest.mark.parametrize("action", ("compute", "report-only", "plot-only"))
def test_recovery_cached_entrypoints_do_zero_new_numerics(tmp_path, monkeypatch, capsys, action):
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    import scipy.integrate
    identity = {"synthetic": "stable-recovery-cache"}
    output = tmp_path / "results"
    bundle = _synthetic_recovery_bundle(output / "fixture", identity)

    def forbidden(*args, **kwargs):
        raise AssertionError("Cached command attempted new integration, profiling, or derivation")

    monkeypatch.setattr(diagnostic, "run_recovery", forbidden)
    monkeypatch.setattr(diagnostic, "profile_equivalence", forbidden)
    monkeypatch.setattr(diagnostic, "initial_compatibility", forbidden)
    monkeypatch.setattr(diagnostic, "diagnostic_pair", forbidden)
    monkeypatch.setattr(pilot, "integrate_case", forbidden)
    monkeypatch.setattr(rod, "derive_polynomials", forbidden)
    monkeypatch.setattr(scipy.integrate, "solve_ivp", forbidden)
    monkeypatch.setattr(scipy.integrate, "Radau", forbidden)
    monkeypatch.setattr(diagnostic, "recovery_identity", lambda *args: ("fixture", identity))
    argv = ["diagnose_weakly_nonlinear_planar_rod.py"]
    argv += (["--compute", "--output-dir", str(output)] if action == "compute"
             else ["--"+action, str(bundle)])
    monkeypatch.setattr(sys, "argv", argv)
    diagnostic.main()
    result = json.loads(capsys.readouterr().out)
    assert result["integrations"] == 0
    assert result["profiling_seconds"] == 0
    assert result["symbolic_derivations"] == 0
    assert result["root_solves"] == 0
    if action == "plot-only":
        assert (bundle / "figures/contraction_projection.pdf").is_file()
        assert (bundle / "figures/contraction_localization.png").is_file()


def test_recovery_identity_locks_budget_and_keeps_source_and_production_separate():
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    if not HISTORICAL.is_dir():
        pytest.skip("Historical bundle unavailable")
    _, identity = diagnostic.recovery_identity(HISTORICAL, None, "compute")
    policy = identity["policy"]
    assert policy["budget_seconds"] == 900
    assert policy["max_short_integrations"] == 3
    assert policy["max_full_integrations"] == 1
    assert policy["p"] == 48
    assert policy["amplitude_over_h"] == .05
    assert policy["full_time_level"] == "tight"
    assert policy["short_p48_time_level"] == "allowed_extra"
    assert policy["forecast_factor"] == 1.25
    assert policy["historical_data_recomputed"] is False
    assert policy["previous_small_amplitude_unchanged"] is True


def test_saved_p48_deferral_obeys_forecast_inequality_without_more_integrations():
    base = ROOT / "results/weakly_nonlinear_planar_recovery"
    if not (base / "current.json").is_file():
        pytest.skip("Recovery result bundle unavailable")
    location = Path(json.loads((base / "current.json").read_text(encoding="utf8"))["bundle"])
    summary = json.loads((location / "summary.json").read_text(encoding="utf8"))
    decision = summary["p48_decision"]
    short = summary["controls"]["p48_strict_short"]["integration_seconds"]
    assert summary["short_integrations"] == 3
    assert summary["full_integrations"] <= 1
    assert summary["budget_used_seconds"] <= 900
    assert decision["forecast_full_seconds"] == short*50*1.25
    assert decision["remaining_integration_budget_seconds"] == 900-summary["budget_used_seconds"]
    if decision["forecast_full_seconds"] > decision["remaining_integration_budget_seconds"]:
        assert summary["full_integrations"] == 0
        assert summary["statuses"]["NLSP_PLANAR_P48_SPATIAL_CHECK"] == "REFINEMENT_DEFERRED_BY_BUDGET"


def test_production_forecast_is_unreachable_if_any_short_control_failed():
    """Execute only the actual guard AST, never its numerical surroundings."""
    import ast
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    tree = ast.parse(Path(diagnostic.__file__).read_text(encoding="utf8"))
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == "run_recovery")
    guard = next(node for node in function.body if isinstance(node, ast.If)
                 and any(isinstance(value, ast.Constant)
                         and value.value == "NOT_RUN_FAILED_SHORT_CONTROL"
                         for value in ast.walk(node)))
    forecast = next(node for node in function.body if isinstance(node, ast.Assign)
                    and any(isinstance(target, ast.Name) and target.id == "forecast"
                            for target in node.targets))
    assert guard.lineno < forecast.lineno
    wrapper = ast.parse("def isolated_guard(controls, statuses, result):\n    pass\n")
    wrapper.body[0].body = [guard, ast.Return(value=ast.Constant(value="CONTINUE_TO_FORECAST"))]
    ast.fix_missing_locations(wrapper)
    namespace = {}
    exec(compile(wrapper, "<actual_short_control_guard>", "exec"), namespace)
    for failed_name in ("old_p32", "new_p32", "p48_strict_short"):
        controls = {name: {"status": "PARTIAL" if name == failed_name else "PASS"}
                    for name in ("old_p32", "new_p32", "p48_strict_short")}
        statuses, result = {}, {"full_integrations": 0}
        actual = namespace["isolated_guard"](controls, statuses, result)
        assert actual is result
        assert statuses["NLSP_PLANAR_P48_SPATIAL_CHECK"] == "NOT_RUN_FAILED_SHORT_CONTROL"
    controls = {name: {"status": "PASS"}
                for name in ("old_p32", "new_p32", "p48_strict_short")}
    assert namespace["isolated_guard"](controls, {}, {}) == "CONTINUE_TO_FORECAST"


def test_partial_p48_metadata_keeps_actual_comparison_end_and_master_end_distinct():
    """Execute the two production metadata assignments without any ODE call."""
    import ast
    from scripts.analysis import diagnose_weakly_nonlinear_planar_rod as diagnostic
    tree = ast.parse(Path(diagnostic.__file__).read_text(encoding="utf8"))
    function = next(node for node in tree.body
                    if isinstance(node, ast.FunctionDef) and node.name == "run_recovery")
    assignments = {}
    for node in function.body:
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if (isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name)
                    and target.value.id == "stats" and isinstance(target.slice, ast.Constant)
                    and target.slice.value in ("master_time_end", "time_end")):
                assignments[target.slice.value] = node
    assert set(assignments) == {"master_time_end", "time_end"}
    metadata_body = ast.Module(body=[assignments["master_time_end"], assignments["time_end"]],
                               type_ignores=[])
    ast.fix_missing_locations(metadata_body)
    present = np.array([0., .05, .10, .15, .20, .25])
    comparison_indices = np.array([0, 2, 4])
    stats = {"status": "PARTIAL", "time_end": .25, "target_time_end": 99.}
    exec(compile(metadata_body, "<actual_partial_p48_metadata>", "exec"),
         {"stats": stats, "present": present, "ids": comparison_indices})
    assert stats["master_time_end"] == .25
    assert stats["time_end"] == .20
    assert stats["time_end"] == present[comparison_indices][-1]
    assert stats["target_time_end"] == 99.
