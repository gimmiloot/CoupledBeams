"""Saved-action Stage A checks, with no ODE, static or FEM jobs."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from scripts.lib import nlsp_spatial_verification_checks as checks
from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin
from scripts.lib.weakly_nonlinear_spatial_dynamics import SpatialGalerkin

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def frozen_model():
    path = ROOT / "results/weakly_nonlinear_spatial_rod/a9cedd4b6de99295/result.json"
    polynomials = json.loads(path.read_text(encoding="utf8"))["polynomials"]
    return SimpleNamespace(T4=rod.Polynomial.deserialize(polynomials["T4"]),
                           V4=rod.Polynomial.deserialize(polynomials["V4"]),
                           residual_a=tuple(rod.Polynomial.deserialize(p) for p in polynomials["residuals_A"]),
                           symbols={name: rod.Polynomial.symbol(name) for name in rod.SYMBOL_ORDER})


@pytest.fixture(scope="module")
def coefficients():
    preflight = json.loads((ROOT / "results/nlsp_linear_rectangular_3d_fem/4262efa427b03dad/preflight.json").read_text(encoding="utf8"))
    return rod.RodCoefficients(**preflight["coefficients"])


@pytest.fixture
def small_disc(coefficients, frozen_model):
    return SpatialGalerkin(coefficients, 8, model=frozen_model)


def _write_manifest(folder, filenames):
    artifacts = {name: hashlib.sha256((folder / name).read_bytes()).hexdigest() for name in filenames}
    (folder / "manifest.json").write_text(json.dumps({"artifacts": artifacts}), encoding="utf8")


def test_exact_frozen_action_planar_parity_variation_and_energy(frozen_model):
    result = checks.exact_action_checks(frozen_model)
    assert result["status"] == "PASS"
    assert (result["T_terms"], result["V_terms"]) == (36, 90)
    assert result["kinetic_velocity_degrees"] == [2]
    assert result["new_symbolic_derivations"] == 0
    assert all(row["difference_terms"] == 0 for row in result["checks"].values())
    assert len([k for k in result["checks"] if k.startswith("frozen_variational_residual_")]) == 7


def test_action_check_detects_bad_residual_without_modifying_original(frozen_model):
    damaged = SimpleNamespace(**vars(frozen_model))
    damaged.residual_a = (frozen_model.residual_a[0] + rod.Polynomial.symbol("w"),) + frozen_model.residual_a[1:]
    result = checks.exact_action_checks(damaged)
    assert result["status"] == "FAIL"
    assert result["checks"]["frozen_variational_residual_u"]["status"] == "FAIL"
    assert checks.exact_action_checks(frozen_model)["status"] == "PASS"


def test_streamed_saved_rows_keep_requested_order(tmp_path):
    data = np.arange(350, dtype=float).reshape(50, 7)
    path = tmp_path / "rows.npz"
    np.savez_compressed(path, q=data)
    np.testing.assert_array_equal(checks.read_npz_rows(path, "q", [49, 0, 23]), data[[49, 0, 23]])


@pytest.mark.parametrize("indices", ([0, 0], [-1], [50]))
def test_saved_rows_reject_fake_or_duplicate_coverage(tmp_path, indices):
    path = tmp_path / "rows.npz"
    np.savez_compressed(path, q=np.zeros((50, 7)))
    with pytest.raises(ValueError):
        checks.read_npz_rows(path, "q", indices)


def test_streamed_saved_rows_reject_fortran_order(tmp_path):
    path = tmp_path / "rows.npz"
    np.savez_compressed(path, q=np.asfortranarray(np.zeros((10, 7))))
    with pytest.raises(ValueError, match="C-order"):
        checks.read_npz_rows(path, "q", [1])


def test_source_manifest_and_hash_rejection_preserve_manifest(tmp_path):
    artifact = tmp_path / "source.json"
    artifact.write_text('{"x":1}', encoding="utf8")
    _write_manifest(tmp_path, [artifact.name])
    before = (tmp_path / "manifest.json").read_bytes()
    result = checks.source_evidence({"a": artifact})
    assert result["a"]["sha256"] == hashlib.sha256(artifact.read_bytes()).hexdigest()
    artifact.write_text('{"x":2}', encoding="utf8")
    with pytest.raises(ValueError, match="hash mismatch"):
        checks.source_evidence({"a": artifact})
    assert (tmp_path / "manifest.json").read_bytes() == before


def test_missing_source_is_not_replaced_by_new_calculation(tmp_path):
    with pytest.raises(FileNotFoundError, match="absent"):
        checks.source_evidence({"a": tmp_path / "missing.npz"})


def test_explicit_zero_reference_policy_has_no_silent_floor():
    assert checks.difference(np.zeros(3), np.zeros(3))["status"] == "PASS"
    result = checks.difference(np.array([0., 1e-30]), np.zeros(2))
    assert result["status"] == "FAIL"
    assert result["relative_max"] is None
    assert result["reference_max"] == 0


def test_planar_embedding_has_exact_order_and_endpoint_values(small_disc):
    q = np.arange(4 * small_disc.n, dtype=float) * 1e-7
    embedded = checks.embed_planar(small_disc, q)
    ids = checks.planar_ids(small_disc)
    np.testing.assert_array_equal(embedded[ids], q)
    assert not embedded[np.setdiff1d(np.arange(small_disc.ndof), ids)].any()
    np.testing.assert_array_equal(small_disc.reconstruct(embedded, points=np.array([0., 1.])), 0.)
    with pytest.raises(ValueError, match="same p"):
        checks.embed_planar(small_disc, q[:-1])


def test_planar_energy_mass_inertia_rhs_jacobian_all_match(small_disc, frozen_model):
    old = PlanarGalerkin(small_disc.coefficients, 8, model=frozen_model)
    q, velocity, _ = checks.smooth_coordinates(small_disc)
    ids = checks.planar_ids(small_disc)
    result = checks.compare_planar_state(small_disc, old, q[ids], velocity[ids])
    assert result["status"] == "PASS"
    assert result["checks"]["inactive_rhs"]["absolute_max"] == 0
    assert result["checks"]["inactive_active_jacobian"]["absolute_max"] == 0


def test_variational_energy_and_jacobian_preserve_declared_thresholds(small_disc):
    result = checks.variational_checks(small_disc)
    assert result["weak_flux_projection"]["status"] == "PASS"
    assert result["weak_flux_projection"]["relative_threshold"] == 2e-12
    assert result["weak_flux_projection"]["absolute_threshold"] == 2e-12
    assert result["weak_flux_projection"]["normalization"] == "unchanged FEM-2 uncancelled potential-work policy"
    assert result["strong_projection_strict"]["status"] == "PASS"
    assert result["energy_rate"]["status"] == "PASS"
    assert result["energy_rate"]["threshold"] == 2e-12
    assert result["jacobian_directional"]["status"] == "PASS"
    assert result["jacobian_directional"]["threshold"] == 2e-7
    assert result["boundary_work"]["essential_endpoint_max"] == 0
    assert result["boundary_work"]["explicit_boundary_work_max"] == 0


def test_weak_flux_continuum_assembly_matches_discrete_action(small_disc):
    q, velocity, acceleration = checks.smooth_coordinates(small_disc)
    weak = checks.weak_flux_projection(small_disc, q, velocity, acceleration)
    action = (small_disc.mass_matrix(q) @ acceleration + small_disc.inertial_terms(q, velocity) +
              small_disc.potential(q)["gradient"])
    np.testing.assert_allclose(weak["residual"], action, rtol=2e-12, atol=2e-18)
    assert weak["boundary_work_max"] == 0
    assert "- [test*flux]_ends" in weak["representation"]


def test_nonessential_test_functions_keep_nonzero_boundary_flux(small_disc):
    q, velocity, acceleration = checks.smooth_coordinates(small_disc)
    constant = tuple(np.ones((small_disc.nq, 1)) for _ in range(7))
    derivative = tuple(np.zeros((small_disc.nq, 1)) for _ in range(7))
    endpoints = tuple(np.ones((2, 1)) for _ in range(7))
    weak = checks.weak_flux_projection(small_disc, q, velocity, acceleration,
                                      test_matrices=constant, test_derivatives=derivative,
                                      endpoint_test_values=endpoints)
    assert weak["boundary_work_max"] > 1e-8
    np.testing.assert_array_equal(weak["residual"], weak["volume_terms"] - weak["boundary_work"])
    variables = []
    for state, order in ((q, 0), (q, 1), (velocity, 0), (q, 2), (velocity, 1), (acceleration, 0)):
        variables.extend(small_disc.reconstruct(state, derivative=order).T)
    strong = small_disc._residual.evaluate(np.asarray(variables)) @ small_disc.weights
    np.testing.assert_allclose(weak["residual"], strong, rtol=2e-12, atol=2e-18)


def test_stable_weak_check_does_not_rename_failed_strong_projection(small_disc, monkeypatch):
    original = small_disc.weak_residual
    monkeypatch.setattr(small_disc, "weak_residual", lambda *args: original(*args) + 1e-10)
    result = checks.variational_checks(small_disc)
    assert result["strong_projection_strict"]["status"] == "PARTIAL"
    assert result["strong_projection_strict"]["relative_threshold"] == 2e-12
    assert result["weak_flux_projection"]["status"] == "PASS"


def test_no_phase_fitting_in_saved_timestamp_selection(tmp_path):
    times = np.array([0., .11, .52, .77, 1.])
    q = np.arange(15, dtype=float).reshape(5, 3)
    path = tmp_path / "saved.npz"
    np.savez_compressed(path, times=times, q=q, velocity=2 * q)
    selected, rows, result, velocity = checks._saved_states(path, (0., .25, .5, .75, 1.))
    np.testing.assert_array_equal(selected, times[rows])
    np.testing.assert_array_equal(result, q[rows])
    np.testing.assert_array_equal(velocity, 2 * q[rows])


def test_run_checks_has_zero_ode_static_fem_and_derivation_calls(tmp_path, small_disc, frozen_model, monkeypatch):
    import scipy.integrate
    import subprocess

    def forbidden(*args, **kwargs):
        raise AssertionError("Stage A matrix/reference checks cannot start a physical solve")

    monkeypatch.setattr(scipy.integrate, "solve_ivp", forbidden)
    monkeypatch.setattr(subprocess, "Popen", forbidden)
    monkeypatch.setattr(rod, "derive_polynomials", forbidden)
    reference = []
    for family, block in (("inplane_bending", "timoshenko"), ("outplane_bending", "outplane"),
                          ("torsion", "torsion"), ("axial_mh", "mh")):
        omega = small_disc.linear_eigenpairs(block)["omega"][0]
        reference.append({"family": family, "local_mode": 1, "omega": float(omega)})
    source = tmp_path / "preflight.json"
    source.write_text(json.dumps({"merged_spectrum": reference}), encoding="utf8")
    _write_manifest(tmp_path, [source.name])
    before = source.read_bytes()
    result = checks.run_stage_a_checks({8: small_disc}, frozen_model, {"fem1_preflight": source}, tmp_path / "checks.json")
    assert result["execution_gate"] == "PASS"
    assert all(value == 0 for value in result["scientific_calls"].values())
    assert result["strict_float64_qualification"]["status"] == "PARTIAL"
    assert result["coordinate_contract"]["derivative_endpoint_constraints"] is False
    assert result["linear_references"]["frequency_units"] == "radians per normalized time"
    assert result["saved_planar_state_status"] == "NOT_RUN"
    assert source.read_bytes() == before
    assert json.loads((tmp_path / "checks.json").read_text())["execution_gate"] == "PASS"


def test_wrong_frequency_preserves_fail_status(small_disc):
    omega = float(small_disc.linear_eigenpairs("torsion")["omega"][0])
    result = checks.linear_reference_checks({8: small_disc}, {"merged_spectrum": [
        {"family": "torsion", "local_mode": 1, "omega": omega * 1.01}]})
    assert result["status"] == "FAIL"
    assert result["threshold"] == 2e-8
