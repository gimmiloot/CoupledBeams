"""FEM-3A 1D orchestration over frozen static coordinates and existing dynamics.

No physics, root search, static solve or symbolic derivation is implemented here.
The full p64 Shen space, audited quartic action and variable mass are retained.
Only the explicitly authorized caller may start the existing Radau runner.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from scripts.analysis import simulate_weakly_nonlinear_planar_rod as runner
from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as static
from scripts.lib import weakly_nonlinear_planar_dynamics as dynamics
from scripts.lib import weakly_nonlinear_spatial_rod as rod

VERSION = "fem3a-frozen-p64-release-reference-v1"
FIELDS = ("u", "w", "theta", "c")


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _checked_artifact(bundle, name):
    bundle = Path(bundle)
    manifest = _read(bundle / "manifest.json")
    artifacts = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
    if name not in artifacts or _sha(bundle / name) != artifacts[name]:
        raise ValueError("Frozen 1D source artifact absent or corrupted: " + str(bundle / name))
    return bundle / name


def load_reference(root, static_bundle, fem1_bundle, action_bundle):
    """Restore precisely the old p64 coordinates, mass scaling and action."""
    root = Path(root)
    sb, fb, ab = (root / Path(p) for p in (static_bundle, fem1_bundle, action_bundle))
    pre_path = _checked_artifact(sb, "one_d_preflight.json")
    coordinates_path = _checked_artifact(sb, "one_d_p64.npz")
    first_path = _checked_artifact(fb, "preflight.json")
    action_path = _checked_artifact(ab, "result.json")
    pre, first, action = _read(pre_path), _read(first_path), _read(action_path)
    if pre["coefficients"] != first["coefficients"]:
        raise ValueError("Static and FEM-1 coefficient references disagree")
    if first["geometry"] != {"L": 1., "b": .2, "h": .1}:
        raise ValueError("FEM-3A reference geometry must remain L=1,b=.20,h=.10")
    if first["material"] != {"E": 1., "rho": 1., "nu": .3, "kappa": 5 / 6}:
        raise ValueError("FEM-3A frozen material mismatch")
    p64 = next(row for row in pre["cases"] if row["p"] == 64)
    if p64["fields"] != list(FIELDS) or p64["ndof"] != 252 or p64["nq"] != 129:
        raise ValueError("Frozen p64 field ordering/dimension/quadrature mismatch")
    pol = action["polynomials"]
    model = SimpleNamespace(T4=rod.Polynomial.deserialize(pol["T4"]),
        V4=rod.Polynomial.deserialize(pol["V4"]),
        residual_a=tuple(rod.Polynomial.deserialize(a) for a in pol["residuals_A"]),
        symbols={a: rod.Polynomial.symbol(a) for a in rod.SYMBOL_ORDER})
    coefficients = rod.RodCoefficients(**pre["coefficients"])
    disc = dynamics.PlanarGalerkin(coefficients, 64, length=1., nq=129, model=model, whiten=True)
    with np.load(coordinates_path, allow_pickle=False) as data:
        saved = {name: data[name].copy() for name in
            ("s", "q_linear", "q_nonlinear", "raw_linear", "raw_nonlinear", "linear", "nonlinear")}
    for name in ("q_linear", "q_nonlinear", "raw_linear", "raw_nonlinear"):
        if saved[name].shape != (252,) or not np.isfinite(saved[name]).all():
            raise ValueError("Invalid frozen coordinates: " + name)
        saved[name].flags.writeable = False
    omega = next(row["omega"] for row in first["merged_spectrum"]
        if row["family"] == "inplane_bending" and row["local_mode"] == 1)
    return {"disc": disc, "saved": saved, "preflight": pre,
        "omega1": float(omega), "T1": 2 * math.pi / omega,
        "source": {"version": VERSION, "static_bundle": str(sb.relative_to(root)),
            "static_manifest_sha256": _sha(sb / "manifest.json"),
            "coordinates_sha256": _sha(coordinates_path),
            "static_preflight_sha256": _sha(pre_path),
            "fem1_bundle": str(fb.relative_to(root)), "fem1_preflight_sha256": _sha(first_path),
            "action_bundle": str(ab.relative_to(root)), "action_result_sha256": _sha(action_path),
            "p": 64, "nq": 129, "ndof": 252, "fields": list(FIELDS),
            "initial_policy": "reuse_saved_static_coordinates_exactly_no_projection",
            "coordinate_scaling": "same resting-mass Cholesky whitening",
            "symbolic_derivations": 0, "static_equilibrium_solves": 0,
            "new_1D_root_searches": 0}}


def runtime_config(reference, pilot_config):
    """Keep the established policy but evaluate scales on the current section."""
    config = copy.deepcopy(pilot_config)
    config["material_geometry"] = {"E": 1., "rho": 1., "nu": .3,
        "b": .2, "h": .1, "L": 1., "kappa": "5/6"}
    config["spatial"]["degrees"] = [64]
    if config["time_levels"]["tight"] != {
        "rtol": 1e-10, "atol_relative": 1e-10,
        "max_step_cutoff_period_fraction": 1 / 24}:
        raise ValueError("Existing tight time policy has changed")
    if reference["disc"].length != config["material_geometry"]["L"]:
        raise ValueError("1D length and time scaling disagree")
    return config


def preflight_reference(reference, pilot_config):
    """Check reuse and restoring acceleration without any time/eigen solve."""
    runner.load_runtime()
    disc, saved = reference["disc"], reference["saved"]
    config = runtime_config(reference, pilot_config)
    force = static.fem2_line_load(disc, reference["preflight"]["load"]["q"])
    result = {"source": reference["source"], "omega1": reference["omega1"],
        "T1": reference["T1"], "horizon": .05 * reference["T1"],
        "external_force_after_release": 0., "slope_constraints": False,
        "initial_velocities": "identically zero in all 252 coordinates",
        "initial_states": {}, "execution_mode": "EXPLORATORY_NOT_CERTIFIED",
        "admitted": False, "strict_float64_strong_weak": "PARTIAL",
        "strict_relative_threshold": 2e-12}
    for kind, linear in (("linear", True), ("nonlinear", False)):
        q0 = saved["q_" + kind]; v0 = np.zeros(disc.ndof)
        mass = disc.M0 if linear else disc.mass_matrix(q0)
        np.linalg.cholesky(mass)
        gradient = disc.K @ q0 if linear else disc.potential(q0)["gradient"]
        acceleration = disc.acceleration(q0, v0, linear=linear)
        physical = disc.reconstruct(q0, saved["s"])
        raw_error = float(np.max(abs(disc.raw_coefficients(q0) - saved["raw_" + kind])))
        physical_error = float(np.max(abs(physical - saved[kind])))
        released = mass @ acceleration + gradient
        scale = float(np.linalg.norm(gradient) + np.linalg.norm(mass @ acceleration))
        diagnostics = _sample_safety(disc, q0)
        runner.safety_check(disc, q0, config["safety"])
        if diagnostics["relative_mass_lower_bound"] < config["safety"]["min_relative_mass_eigenvalue"]:
            raise ValueError("Initial variable mass below existing safety bound")
        endpoint = float(np.max(abs(disc.reconstruct(q0, [0., 1.]))))
        midpoint_acceleration = float(disc.reconstruct(acceleration, [.5])[0, 1])
        equilibrium = float(np.linalg.norm(gradient - force) / max(np.linalg.norm(force), 1e-30))
        released_relative = float(np.linalg.norm(released) / max(scale, 1e-30))
        if (raw_error > 1e-12 or physical_error > 1e-12 or endpoint != 0.
            or equilibrium > 1e-10 or released_relative > 2e-12 or midpoint_acceleration >= 0):
            raise ValueError("Frozen state/acceleration preflight failed for " + kind)
        result["initial_states"][kind] = {"q0": q0.tolist(), "v0": v0.tolist(),
            "acceleration0": acceleration.tolist(), "mass_cholesky_positive": True,
            "raw_coordinate_reproduction_max_abs": raw_error,
            "physical_profile_reproduction_max_abs": physical_error,
            "essential_endpoint_max_abs": endpoint,
            "loaded_equilibrium_relative_residual": equilibrium,
            "released_action_absolute_residual": float(np.max(abs(released))),
            "released_action_relative_residual": released_relative,
            "midspan_w_acceleration": midpoint_acceleration,
            "initial_w_midspan": float(disc.reconstruct(q0, [.5])[0, 1]),
            "safety": diagnostics}
    amplitude = reference["preflight"]["load"]["linear_w_max"]
    settings = runner.time_settings(disc, amplitude, "tight", config)
    result["time_settings"] = {k: v.tolist() if isinstance(v, np.ndarray) else v for k, v in settings.items()}
    result["time_settings"].update({"characteristic_amplitude": amplitude,
        "mass_scale_source": "current frozen m,jp,L,n; no old dimensional atol reused",
        "method": "Radau", "jacobian": "unchanged analytic variable-mass Jacobian"})
    result["status"] = "PASS"
    result["counters"] = disc.counters()
    return result


def _times(times, horizon):
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or not np.isfinite(times).all()
        or times[0] != 0 or np.any(np.diff(times) <= 0)
        or times[-1] > horizon * (1 + 1e-12)):
        raise ValueError("Actual comparison times must increase from zero within 0.05T1")
    return times


def _sample_safety(disc, coordinate):
    values = disc.reconstruct(coordinate)
    gradients = disc.reconstruct(coordinate, derivative=1)
    us, ws, theta, c = gradients[:, 0], gradients[:, 1], values[:, 2], values[:, 3]
    gamma1 = us + theta * ws - theta**2 / 2 - us * theta**2 / 2
    gamma2 = ws - theta - theta * us - ws * theta**2 / 2 + theta**3 / 6
    low = min(1., float(np.min((1 + c)**2)))
    high = max(1., float(np.max((1 + c)**2)))
    return {"min_one_plus_c": float(np.min(1 + c)), "max_abs_c": float(np.max(abs(c))),
        "max_abs_theta": float(np.max(abs(theta))), "max_abs_u_s": float(np.max(abs(us))),
        "max_abs_w_s": float(np.max(abs(ws))),
        "max_L_abs_theta_s": float(disc.length * np.max(abs(gradients[:, 2]))),
        "max_retained_axial_strain": float(np.max(abs(gamma1))),
        "max_retained_shear_strain": float(np.max(abs(gamma2))),
        "relative_mass_lower_bound": low, "relative_mass_upper_bound": high,
        "relative_mass_condition_upper_bound": high / low,
        "mass_positive": low > 0., "mass_bound_method": "positive weighted theta-Gram Loewner bounds"}


def summarize_reference(reference, trajectory, points=None, *, linear=False):
    """Measure actual samples with the appropriate model's own energy."""
    disc = reference["disc"]
    points = np.linspace(0., disc.length, 41) if points is None else np.asarray(points, dtype=float)
    q, v = np.asarray(trajectory["q"]), np.asarray(trajectory["velocity"])
    if q.shape != v.shape or q.ndim != 2 or q.shape[1] != disc.ndof:
        raise ValueError("1D trajectory dimensions disagree")
    if not np.isfinite(q).all() or not np.isfinite(v).all():
        raise ValueError("Nonfinite 1D trajectory")
    energy, safety, norms, speed_norms = [], [], [], []
    for row, speed in zip(q, v):
        energy.append(float((speed @ disc.M0 @ speed + row @ disc.K @ row) / 2)
            if linear else disc.energy(row, speed))
        safety.append(_sample_safety(disc, row))
        values, speeds = disc.reconstruct(row), disc.reconstruct(speed)
        norms.append(np.sqrt(disc.weights @ (values**2)))
        speed_norms.append(np.sqrt(disc.weights @ (speeds**2)))
    energy = np.asarray(energy)
    return {"times": np.asarray(trajectory["times"]), "q": q, "velocity": v,
        "x": points, "fields": disc.reconstruct_series(q, points),
        "physical_velocities": disc.reconstruct_series(v, points),
        "L2_fields": np.asarray(norms), "L2_velocities": np.asarray(speed_norms),
        "energy": energy, "energy_relative_drift": (energy - energy[0]) / energy[0],
        "diagnostics": {"safety_samples": safety,
            "max_relative_energy_drift": float(np.max(abs((energy - energy[0]) / energy[0]))),
            "energy_definition": "0.5*v.T*M0*v+0.5*q.T*K*q" if linear else "0.5*v.T*M(q)*v+V4(q)",
            "removed_GRAV_potential_included": False,
            "temporal_convergence_claimed": False, "spatial_convergence_claimed": False}}


def exact_linear_reference(reference, times):
    times = _times(times, .05 * reference["T1"])
    disc = reference["disc"]; q0 = reference["saved"]["q_linear"]
    trajectory = disc.linear_reference(q0, np.zeros(disc.ndof), times)
    difference = float(np.max(abs(trajectory["q"][0] - q0)))
    # t=0 is the defined saved state; remove only eigenspace reconstruction roundoff.
    trajectory["q"][0] = q0
    trajectory["velocity"][0] = 0.
    trajectory["metadata"] = {"exact_in_time": True,
        "semidiscrete_spectral_factorization": "full 252-coordinate K,M0; no modal reduction",
        "zero_time_factorization_roundoff_max_abs": difference,
        "new_1D_root_searches": 0, "ODE_integrations": 0,
        "linear_eigendecompositions": disc.linear_eigendecompositions}
    return trajectory


def integrate_nonlinear_reference(reference, times, pilot_config, deadline, *, authorization):
    if not authorization or authorization.get("user_authorized_FEM3A") is not True:
        raise ValueError("Explicit FEM-3A nonlinear 1D execution authorization required")
    times = _times(times, .05 * reference["T1"])
    if not math.isclose(times[-1], .05 * reference["T1"], rel_tol=1e-12, abs_tol=0):
        raise ValueError("The sole nonlinear run must target the authorized 0.05T1 horizon")
    config = runtime_config(reference, pilot_config)
    runner.load_runtime()
    disc = reference["disc"]
    amplitude = reference["preflight"]["load"]["linear_w_max"]
    history, stats = runner.integrate_case(disc, None,
        {"omega": reference["omega1"], "T1": reference["T1"]}, config,
        amplitude / config["material_geometry"]["h"], "tight", times, deadline,
        initial_coordinates=reference["saved"]["q_nonlinear"])
    actual_times = times[:len(history)]
    stats.update({"execution_mode": "EXPLORATORY_NOT_CERTIFIED", "admitted": False,
        "authorization": authorization, "initial_coordinates_reused_exactly": True,
        "external_force_after_release": 0., "new_ODE_integrations": 1,
        "no_dynamic_derivative_constraints": True, "strict_float64_strong_weak": "PARTIAL"})
    return {"times": actual_times, "q": history[:, :disc.ndof],
        "velocity": history[:, disc.ndof:], "stats": stats}
