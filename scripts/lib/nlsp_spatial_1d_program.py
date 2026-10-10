"""Scoped seven-field static/release programme over the frozen action helper.

No native FEM operation is implemented. Static continuation reuses the
existing FEM-2 Newton policy and solver. The accepted-prefix Radau loop follows
the existing planar runner, with canonical seven-field tolerances and safety.
Each scientific attempt is recorded before execution and never retried by a
cache lookup, including after a real failure or an interrupted process.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import shutil
import time

import numpy as np
from scipy.integrate import Radau
from scipy.linalg import cho_factor, cho_solve

from scripts.analysis.verify_nlsp_nonlinear_static_3d_fem import fem2_static_newton
from scripts.lib.nlsp_fem3b_continuation import save_dense_records, evaluate_dense_records
from scripts.lib.weakly_nonlinear_spatial_dynamics import FIELDS

ROOT = Path(__file__).resolve().parents[2]
VERSION = "scoped-seven-field-static-release-radau-v1"
STATIC_POLICY = {"residual_relative_tolerance": 1e-10, "increment_relative_tolerance": 1e-12,
                 "load_steps": 10, "newton_max_iterations": 20, "max_load_subdivisions": 6}
TIME_POLICY = {"rtol": 1e-10, "atol_relative": 1e-10, "max_step_cutoff_period_fraction": 1/24}
SAFETY_POLICY = {"min_one_plus_c": .9, "max_abs_c": .1, "max_abs_theta": .1,
                 "max_abs_axial_gradient": .1, "max_abs_transverse_gradient": .1,
                 "max_L_abs_curvature": .2, "min_relative_mass_eigenvalue": .8}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(v) for v in value]
    return value


def _write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(_plain(value), ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf8")
    temp.replace(path)


def _policies(config=None):
    config = {} if config is None else config
    static = config.get("static_policy", STATIC_POLICY)
    timing = config.get("time_policy", TIME_POLICY)
    safety = config.get("safety", SAFETY_POLICY)
    if static != STATIC_POLICY or timing != TIME_POLICY or safety != SAFETY_POLICY:
        raise ValueError("The established static, tight-time, and safety policies must remain unchanged")
    return static, timing, safety


def line_load(disc, q_w, q_v):
    if not all(math.isfinite(float(q)) for q in (q_w, q_v)):
        raise ValueError("Finite frozen line-load components required")
    force = np.zeros(disc.ndof)
    for field, component in (("w", q_w), ("v", q_v)):
        i = FIELDS.index(field)
        force[disc.slices[field]] = component*(disc.B[i].T@disc.weights)
    return force


def time_settings(disc, amplitude, policy=None):
    policy = TIME_POLICY if policy is None else policy
    if policy != TIME_POLICY or not math.isfinite(amplitude) or amplitude <= 0:
        raise ValueError("Positive characteristic amplitude and unchanged tight policy required")
    p, length, n = disc.coefficients, disc.length, disc.n
    cutoff = math.sqrt(p.C/p.jp)
    scales = np.array((amplitude, amplitude, amplitude, amplitude/length,
                       amplitude/length, amplitude/length, amplitude/length))
    masses = np.array((p.m, p.m, p.m, p.jp+p.jb, p.jb, p.jp, p.jp))
    coordinate_scales = scales*np.sqrt(masses*length/n)
    atol_q = np.repeat(coordinate_scales, n)*policy["atol_relative"]
    return {"rtol": policy["rtol"], "atol": np.r_[atol_q, atol_q*cutoff],
            "max_step": 2*np.pi/cutoff*policy["max_step_cutoff_period_fraction"],
            "coordinate_scales": coordinate_scales, "field_scales": scales,
            "resting_field_masses": masses, "velocity_scale_multiplier": cutoff,
            "atol_relative": policy["atol_relative"]}


def safety_check(disc, coordinate, policy=None):
    policy = SAFETY_POLICY if policy is None else policy
    values, gradients = disc.reconstruct(coordinate), disc.reconstruct(coordinate, derivative=1)
    if not np.isfinite(values).all() or not np.isfinite(gradients).all():
        raise ArithmeticError("NONFINITE_RECONSTRUCTED_FIELD")
    result = {"min_one_plus_c": float(np.min(1+values[:, 6])),
              "max_abs_c": float(np.max(abs(values[:, 6]))),
              "max_abs_theta": float(np.max(abs(values[:, 3:6]))),
              "max_rotation_vector_norm": float(np.max(np.linalg.norm(values[:, 3:6], axis=1))),
              "max_abs_axial_gradient": float(np.max(abs(gradients[:, 0]))),
              "max_abs_transverse_gradient": float(np.max(abs(gradients[:, 1:3]))),
              "max_L_abs_curvature": float(disc.length*np.max(abs(gradients[:, 3:6])))}
    if result["min_one_plus_c"] <= policy["min_one_plus_c"] or any(
            result[key] > policy[key] for key in policy if key in result and key != "min_one_plus_c"):
        raise ArithmeticError("Declared small-neighborhood safety gate: "+str(result))
    mass = disc.mass_spectral_bounds(coordinate)
    if not mass["mass_positive"] or mass["relative_mass_lower_bound"] < policy["min_relative_mass_eigenvalue"]:
        raise ArithmeticError("Variable rotational mass violates the existing positivity bound")
    return {**result, **mass}


def solve_static_pair(disc, loads, policy=None):
    policy = STATIC_POLICY if policy is None else policy
    if policy != STATIC_POLICY:
        raise ValueError("Frozen static Newton controls changed")
    force = line_load(disc, *loads)
    if np.linalg.norm(force) == 0:
        raise ValueError("At least one nonzero transverse dead load required")
    linear = cho_solve(cho_factor(disc.K, lower=True, check_finite=True), force, check_finite=True)
    result = fem2_static_newton(disc, force, policy)
    return {"force": force, "q_linear": linear, "q_nonlinear": result["coordinate"],
            "nonlinear": {k: v for k, v in result.items() if k != "coordinate"},
            "linear_residual_relative": float(np.linalg.norm(disc.K@linear-force)/np.linalg.norm(force)),
            "status": result["status"]}


def _snapshot(directory, phase):
    from scripts.lib import weakly_nonlinear_spatial_dynamics as spatial
    from scripts.lib import weakly_nonlinear_planar_dynamics as planar
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as static
    paths = [Path(__file__), Path(spatial.__file__), Path(planar.__file__), Path(rod.__file__), Path(static.__file__)]
    folder = Path(directory)/"execution_code"/phase
    folder.mkdir(parents=True, exist_ok=False)
    for path in paths:
        shutil.copyfile(path, folder/path.name)
    return {path.relative_to(ROOT).as_posix(): _sha(path) for path in paths}


def integrate_prefix(disc, initial, times, settings, safety, deadline):
    """One Radau attempt; accepted prefix and exact existing dense polynomials."""
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0 or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0)):
        raise ValueError("Output timestamps must increase from zero")
    initial = np.asarray(initial, dtype=float)
    if initial.shape != (disc.ndof,) or not np.isfinite(initial).all():
        raise ValueError("Finite full seven-field initial coordinates required")
    history = np.empty((len(times), 2*disc.ndof))
    history[0] = np.r_[initial, np.zeros(disc.ndof)]
    cursor, records, steps = 1, [], []
    started = time.perf_counter()
    disc.reset_counters()
    failure = None
    solver = None

    def rhs(t, state):
        if not np.isfinite(state).all():
            raise ArithmeticError("NONFINITE_STATE")
        safety_check(disc, state[:disc.ndof], safety)
        result = disc.rhs(t, state)
        if not np.isfinite(result).all():
            raise ArithmeticError("NONFINITE_RHS")
        return result

    try:
        if time.perf_counter() > deadline:
            raise TimeoutError("PREDECLARED_COMPUTATIONAL_BUDGET_EXHAUSTED")
        solver = Radau(rhs, 0., history[0], float(times[-1]), jac=disc.jacobian,
                       rtol=settings["rtol"], atol=settings["atol"], max_step=settings["max_step"])
        while solver.status == "running":
            if time.perf_counter() > deadline:
                raise TimeoutError("PREDECLARED_COMPUTATIONAL_BUDGET_EXHAUSTED")
            old = solver.t
            message = solver.step()
            if solver.status == "failed":
                raise ArithmeticError(message or "RADAU_FAILED")
            steps.append(solver.t-old)
            dense = solver.dense_output()
            records.append({"t_old": float(dense.t_old), "t": float(dense.t),
                            "y_old": dense.y_old.copy(), "Q": dense.Q.copy()})
            end = int(np.searchsorted(times, solver.t, side="right"))
            if end > cursor:
                history[cursor:end] = dense(times[cursor:end]).T
                cursor = end
    except (ArithmeticError, ValueError, np.linalg.LinAlgError, TimeoutError) as error:
        failure = type(error).__name__+": "+str(error)
    stats = {"status": "PASS" if cursor == len(times) and failure is None else "PARTIAL",
             "failure": failure, "p": disc.p, "ndof": disc.ndof, "nq": disc.nq,
             "time_end": float(times[cursor-1]), "actual_accepted_end": float(records[-1]["t"] if records else 0.),
             "target_time_end": float(times[-1]), "samples": cursor,
             "accepted_internal_steps": len(steps), "internal_time_steps": steps,
             "min_internal_step": min(steps, default=0.), "max_internal_step": max(steps, default=0.),
             "nfev": 0 if solver is None else solver.nfev, "njev": 0 if solver is None else solver.njev,
             "nlu": 0 if solver is None else solver.nlu, "counters": disc.counters(),
             "integration_seconds": time.perf_counter()-started, "time_settings": _plain(settings),
             "external_force_after_release": 0., "initial_velocities": "all zero",
             "execution_mode": "EXPLORATORY_NOT_CERTIFIED", "admitted": False,
             "strict_float64_strong_weak": "PARTIAL", "automatic_retry": False}
    return history[:cursor], records, stats


def summarize_trajectory(disc, history, safety):
    energies, diagnostics = [], []
    for state in history:
        q, v = state[:disc.ndof], state[disc.ndof:]
        diagnostics.append(safety_check(disc, q, safety))
        energies.append(disc.energy(q, v))
    energies = np.asarray(energies)
    if energies[0] == 0 or not np.isfinite(energies).all():
        raise ArithmeticError("Invalid own mechanical energy reference")
    drift = (energies-energies[0])/energies[0]
    return energies, {"max_relative_energy_drift": float(np.max(abs(drift))),
        "energy_gate": 1e-6, "energy_status": "PASS" if np.max(abs(drift)) <= 1e-6 else "PARTIAL",
        "energy_definition": "T4(q,v)+V4(q); removed load potential excluded",
        "relative_mass_lower_bound": min(r["relative_mass_lower_bound"] for r in diagnostics),
        "relative_mass_upper_bound": max(r["relative_mass_upper_bound"] for r in diagnostics),
        "max_mass_condition_bound": max(r["relative_mass_condition_upper_bound"] for r in diagnostics),
        "min_one_plus_c": min(r["min_one_plus_c"] for r in diagnostics),
        "max_abs_c": max(r["max_abs_c"] for r in diagnostics),
        "max_rotation_vector_norm": max(r["max_rotation_vector_norm"] for r in diagnostics),
        "max_abs_axial_gradient": max(r["max_abs_axial_gradient"] for r in diagnostics),
        "max_abs_transverse_gradient": max(r["max_abs_transverse_gradient"] for r in diagnostics),
        "max_L_abs_curvature": max(r["max_L_abs_curvature"] for r in diagnostics)}


def _finish(directory, result):
    result["artifact_hashes"] = {p.relative_to(directory).as_posix(): _sha(p) for p in sorted(directory.rglob("*"))
                                  if p.is_file() and p.name != "case.json"}
    _write(directory/"case.json", result)
    return result


def cached_case(directory):
    directory = Path(directory)
    if (directory/"case.json").exists():
        result = json.loads((directory/"case.json").read_text(encoding="utf8"))
        for name, digest in result["artifact_hashes"].items():
            if _sha(directory/name) != digest:
                raise ValueError("Immutable seven-field case artifact changed: "+name)
        return result
    if directory.exists() and any(directory.iterdir()):
        raise RuntimeError("Interrupted seven-field scientific attempt: no automatic retry")
    return None


def run_case(disc, loads, times, deadline, bundleprefix, *, authorization, config=None):
    """Static L/NL pair and one nonlinear IVP, or immutable cached evidence."""
    if authorization.get("user_authorized_spatial_stage_b") is not True:
        raise ValueError("Explicit seven-field Stage B scientific authorization required")
    static_policy, timing, safety = _policies(config)
    times = np.asarray(times, dtype=float)
    if (times.ndim != 1 or len(times) < 2 or times[0] != 0 or not np.isfinite(times).all()
            or np.any(np.diff(times) <= 0)):
        raise ValueError("Output timestamps must increase from zero")
    request = {"loads": list(map(float, loads)), "p": disc.p, "nq": disc.nq,
               "length": disc.length, "whiten": disc.whiten,
               "coefficients": disc.coefficients.values(), "authorization": authorization,
               "output_times_sha256": hashlib.sha256(times.tobytes()).hexdigest(),
               "static_policy": static_policy, "time_policy": timing, "safety": safety}
    directory = Path(bundleprefix)
    old = cached_case(directory)
    if old is not None:
        if old.get("request") != request:
            raise ValueError("Immutable case request differs from cached scientific attempt")
        return old
    directory.mkdir(parents=True, exist_ok=True)
    result = {"version": VERSION, "status": "PARTIAL", "authorization": authorization,
              "request": request,
              "p": disc.p, "ndof": disc.ndof, "nq": disc.nq, "fields": list(FIELDS),
              "loads": list(map(float, loads)), "requested_end": float(times[-1]),
              "calls": {"nonlinear_static": 0, "linear_static": 0, "nonlinear_ODE": 0,
                        "exact_time_linear_reference": 0, "native_FEM": 0, "Gmsh": 0},
              "automatic_retry": False, "attempts": []}
    static_hashes = _snapshot(directory, "static")
    result["attempts"].append({"kind": "linear_and_nonlinear_static", "status": "STARTED", "code_hashes": static_hashes})
    _write(directory/"attempt_ledger.json", result)
    try:
        if time.perf_counter() > deadline:
            raise TimeoutError("PREDECLARED_COMPUTATIONAL_BUDGET_EXHAUSTED")
        started = time.perf_counter()
        result["calls"]["linear_static"] = result["calls"]["nonlinear_static"] = 1
        initial = solve_static_pair(disc, loads, static_policy)
        result["static_seconds"] = time.perf_counter()-started
        result["static"] = {k: v for k, v in initial.items() if k not in ("force", "q_linear", "q_nonlinear")}
        np.savez_compressed(directory/"static_states.npz", force=initial["force"],
                            q_linear=initial["q_linear"], q_nonlinear=initial["q_nonlinear"],
                            raw_linear=disc.raw_coefficients(initial["q_linear"]),
                            raw_nonlinear=disc.raw_coefficients(initial["q_nonlinear"]),
                            initial_velocity=np.zeros(disc.ndof))
        _write(directory/"static.json", result["static"])
        result["attempts"][-1]["status"] = initial["status"]
        if initial["status"] != "PASS":
            raise ArithmeticError("Bounded nonlinear static continuation failed")
        safety_check(disc, initial["q_linear"], safety)
        safety_check(disc, initial["q_nonlinear"], safety)
        # Own linear equilibrium is retained, never amplitude-aligned to NL.
        amplitude = float(np.max(abs(disc.reconstruct(initial["q_linear"])[:, (1, 2)])))
        settings = time_settings(disc, amplitude, timing)
        result["initial_safety"] = {kind: disc.diagnostics(initial["q_"+kind], np.zeros(disc.ndof))
                                    for kind in ("linear", "nonlinear")}
        result["initial_release"] = {}
        for kind in ("linear", "nonlinear"):
            q0 = initial["q_"+kind]
            a0 = disc.acceleration(q0, np.zeros(disc.ndof), linear=kind == "linear")
            mass = disc.M0 if kind == "linear" else disc.mass_matrix(q0)
            gradient = disc.K@q0 if kind == "linear" else disc.potential(q0)["gradient"]
            balance = mass@a0+gradient
            result["initial_release"][kind] = {"absolute_residual": float(np.linalg.norm(balance)),
                "relative_residual": float(np.linalg.norm(balance)/(np.linalg.norm(mass@a0)+np.linalg.norm(gradient))),
                "zero_velocities": True, "removed_dead_force": True}
        # Prepare the exact full semidiscrete linear reference and save its
        # existing complete operator before the nonlinear physical attempt.
        # Later comparison at actual FEM times reads these saved arrays only.
        linear_prepared = disc.linear_reference(initial["q_linear"], np.zeros(disc.ndof), times)
        linear_prepared["q"][0], linear_prepared["velocity"][0] = initial["q_linear"], 0.
        np.savez_compressed(directory/"linear_reference_prepared.npz", **linear_prepared)
        np.savez_compressed(directory/"linear_modes.npz", **disc._linear_modes[None], M0=disc.M0)
        result["calls"]["exact_time_linear_reference"] = 1
        ode_hashes = _snapshot(directory, "nonlinear_ODE")
        result["attempts"].append({"kind": "nonlinear_ODE", "status": "STARTED", "code_hashes": ode_hashes})
        result["calls"]["nonlinear_ODE"] = 1
        _write(directory/"attempt_ledger.json", result)
        history, records, stats = integrate_prefix(disc, initial["q_nonlinear"], times, settings, safety, deadline)
        actual_times = times[:len(history)]
        np.savez_compressed(directory/"trajectory.npz", times=actual_times,
                            q=history[:, :disc.ndof], velocity=history[:, disc.ndof:])
        if records:
            result["dense_output"] = save_dense_records(directory/"accepted_dense.npz", records)
        _write(directory/"integration.json", stats)
        result["integration"] = stats
        result["attempts"][-1]["status"] = stats["status"]
        _write(directory/"attempt_ledger.json", result)
        # This linear reference is the complete semidiscrete operator, no
        # analytic root search or modal truncation. Only the actual prefix.
        linear = {name: value[:len(history)] for name, value in linear_prepared.items()}
        np.savez_compressed(directory/"linear_trajectory.npz", **linear)
        energy, diagnostics = summarize_trajectory(disc, history, safety)
        np.savez_compressed(directory/"energy.npz", times=actual_times, energy=energy,
                            relative_drift=(energy-energy[0])/energy[0])
        result["diagnostics"] = diagnostics
        result["status"] = stats["status"]
    except (ArithmeticError, ValueError, np.linalg.LinAlgError, TimeoutError) as error:
        result["failure"] = type(error).__name__+": "+str(error)
        result["status"] = "PARTIAL"
        if result["attempts"][-1]["status"] == "STARTED":
            result["attempts"][-1]["status"] = "FAIL"
    _write(directory/"attempt_ledger.json", result)
    return _finish(directory, result)
