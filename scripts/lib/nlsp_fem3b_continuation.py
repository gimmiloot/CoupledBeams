"""Scoped FEM-3B orchestration for the existing FEM-3AR command.

No constitutive, element, mass, RHS or Jacobian implementation lives here.
The parent remains immutable. A predeclared planning rule selects one horizon
from saved 1D evidence before either of the two native jobs may start.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.metadata
import json
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as resume

base = resume.base
ROOT = resume.ROOT
CONFIG = ROOT / "data/input/nlsp_nonlinear_dynamic_long_horizon.json"
OUTPUT = ROOT / "results/nlsp_nonlinear_dynamic_long_horizon"
read_json, write_json, sha = base.read_json, base.write_json, base.sha
AUTHORIZATION = "explicit_user_FEM3B_2026_10_09"
STATUS_NAMES = ("SOURCE_PRESERVATION", "OLD_SIGNAL_DIAGNOSTIC", "1D_PRELIMINARY",
    "HORIZON_SELECTION", "LINEAR_3D", "NONLINEAR_3D", "PREFIX_REPRODUCTION",
    "RESPONSE_COMPARISON", "DYNAMIC_NONLINEAR_SIGNAL", "ENERGY_DIAGNOSTICS")


def validate_config(c):
    if c.get("schema") != "nlsp-fem3b-long-horizon-v1":
        raise ValueError("FEM-3B schema mismatch")
    auth = c["authorization"]
    if (auth["id"] != AUTHORIZATION or auth["maximum_production_CCX_jobs"] != 2
        or auth["maximum_nonlinear_1D_integrations"] != 2
        or auth["case_order"] != ["linear", "nonlinear"] or auth["automatic_retry"]):
        raise ValueError("Separate bounded FEM-3B authorization required")
    h = c["horizon_selection"]
    if (h["candidate_T1_fractions"] != [.25, .5] or h["preferred_T1_fraction"] != .25
        or h["planning_multiplier"] != 10.
        or h["historical_DAT_FRD_displacement_difference"] != 5e-9
        or h["historical_41_81_correction_profile_difference"] != 7.78842e-9
        or h["rule"] != "prefer_0p25_if_midspan_evolution_signal_ge_10_times_each_historical_indicator_and_observed_p48_p64_evolution_difference_otherwise_0p5"
        or h["physical_validation_threshold"] is not False):
        raise ValueError("Predeclared horizon selection heuristic changed")
    if (c["threads"] != 1 or c["job_timeout_seconds_by_horizon"] != {"0.25": 1200, "0.5": 2400}
        or c["job_memory_limit_bytes"] != 4 * 1024**3 or c["numerical_budget_seconds"] != 6000):
        raise ValueError("FEM-3B bounded resource policy changed")
    if c["one_d"] != {"main_degree": 64, "optional_spatial_degree": 48,
        "target_T1_fraction": .5, "initial_policy": "reuse_saved_static_coordinates_exactly_no_projection",
        "time_level": "tight", "preserve_accepted_Radau_dense_polynomials": True}:
        raise ValueError("Frozen 1D representation or time policy changed")
    if (c["execution_mode"] != "EXPLORATORY_NOT_CERTIFIED" or c["admitted"] is not False
        or any(c[k] for k in ("new_meshes", "new_modal_jobs", "new_static_only_jobs", "new_time_or_space_levels"))):
        raise ValueError("Unauthorized model or scientific extension")
    pairing = c["three_d_pairing"]
    if pairing != {"preferred": "native_equal_timestamps",
        "fallback": "fixed_grid_inside_actual_overlap_linear_and_shape_preserving_cubic_diagnostic",
        "fixed_grid": "T1_over_2000_multiples_inside_overlap", "phase_amplitude_alignment": False}:
        raise ValueError("Predeclared time pairing changed")
    return c


def load_parent(c):
    """Verify the historical manifest and every artifact, without replay."""
    parent = ROOT / c["parent_completed"]["bundle"]
    if sha(parent / "manifest.json") != c["parent_completed"]["manifest_sha256"]:
        raise ValueError("Completed FEM-3AR parent manifest changed")
    summary = resume.validate_cache(parent)
    if (summary["overall"] != "PILOT_COMPLETE_WITH_QUALIFICATIONS"
        or summary["job_calls"]["CCX_production"] != 2
        or any(summary["cases"][k]["status"] != "PASS" for k in ("linear", "nonlinear"))):
        raise ValueError("Wrong or incomplete FEM-3AR parent")
    item = read_json(parent / "provenance.json")
    science = item["config"]
    static_summary = read_json(ROOT / science["source_resume"]["bundle"] / "summary.json")
    binary = Path(static_summary["science_config"]["ccx_exe"])
    if not binary.is_file() or sha(binary) != item["solver_sha256"]:
        raise ValueError("Previously verified CalculiX binary changed or unavailable")
    runtime = {path.name: sha(path) for path in sorted(binary.parent.glob("*.dll"))}
    if runtime != item["runtime_DLLs"]:
        raise ValueError("Previously verified CalculiX runtime libraries changed")
    return parent, item, summary, copy.deepcopy(item["config"])


def validate_cache(bundle):
    summary = base.validate_cache(Path(bundle))
    item = read_json(Path(bundle) / "provenance.json")
    load_parent(validate_config(item["long_horizon_config"]))
    if (Path(bundle) / "horizon_decision.json").exists():
        if sha(Path(bundle) / "horizon_decision.json") != summary["horizon_decision_sha256"]:
            raise ValueError("Frozen pre-FEM horizon decision changed")
    return summary


def save(bundle, item, summary):
    base.finalize(Path(bundle), item, summary)


def existing_attempt(c):
    if not OUTPUT.exists():
        return None
    for bundle in sorted(OUTPUT.iterdir()):
        if not bundle.is_dir() or not (bundle / "provenance.json").exists():
            continue
        item = read_json(bundle / "provenance.json")
        if item.get("authorization", {}).get("id") != c["authorization"]["id"]:
            continue
        if item["long_horizon_config"] != c:
            raise ValueError("FEM-3B authorization already used with different settings")
        if not (bundle / "manifest.json").exists():
            raise RuntimeError("Interrupted unmanifested FEM-3B attempt; no hidden retry")
        return bundle, item, validate_cache(bundle)
    return None


def prepare_stage(config_path=CONFIG):
    """Persist the choice rule and provenance before any new integration."""
    c = validate_config(read_json(config_path))
    found = existing_attempt(c)
    if found:
        return found
    parent, old_item, old_summary, science = load_parent(c)
    helpers = (Path(__file__), Path(resume.__file__), Path(base.__file__),
        Path(base.one.__file__), Path(base.io.__file__), Path(base.one.runner.__file__),
        Path(base.one.dynamics.__file__), Path(base.one.rod.__file__))
    item = {"config": science, "long_horizon_config": c,
        "authorization": c["authorization"], "parent_completed": c["parent_completed"],
        "parent_identity_sha256": hashlib.sha256(json.dumps(old_item, sort_keys=True).encode()).hexdigest(),
        "helper_sha256": {p.relative_to(ROOT).as_posix(): sha(p) for p in helpers},
        "source_mesh_sha256": old_item["source_mesh_sha256"],
        "solver_sha256": old_item["solver_sha256"], "runtime_DLLs": old_item["runtime_DLLs"],
        "python": sys.version, "dependencies": {name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "matplotlib")},
        "HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()}
    key = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()[:16]
    bundle = OUTPUT / key
    if bundle.exists() and any(bundle.iterdir()):
        raise RuntimeError("Existing FEM-3B artifacts cannot be overwritten")
    bundle.mkdir(parents=True)
    write_json(bundle / "provenance.json", item)
    write_json(bundle / "continuation_config.json", c)
    write_json(bundle / "config.json", science)
    write_json(bundle / "predeclared_horizon_rule.json", c["horizon_selection"])
    (bundle / "execution_code").mkdir()
    for path in helpers:
        shutil.copyfile(path, bundle / "execution_code" / path.name)
    for name in ("protocol_evidence.json", "source_input_gate.json"):
        shutil.copyfile(parent / name, bundle / name)
    shutil.copyfile(parent / "one_d_preflight.json", bundle / "historical_p64_preflight.json")
    summary = {"authorization": c["authorization"], "parent_completed": c["parent_completed"],
        "statuses": {"NLSP_FEM3A_" + n: "NOT_RUN" for n in base.STATUS_NAMES},
        "long_horizon_statuses": {"NLSP_FEM3B_" + n: "NOT_RUN" for n in STATUS_NAMES},
        "cases": {}, "attempts": [], "job_calls": {"CCX_production": 0,
            "CCX_fixture": 0, "Gmsh": 0, "1D_nonlinear_ODE": 0, "1D_static": 0,
            "physical_root_search": 0, "symbolic_derivations": 0},
        "one_d_attempts": [], "numerical_seconds": 0., "overall": "NOT_RUN",
        "execution_mode": "EXPLORATORY_NOT_CERTIFIED", "admitted": False,
        "strict_float64_qualification": "PARTIAL",
        "preflight": {"omega1": science["omega1"],
            "T1": base.dynamic_settings(science)["T1"], "target_dynamic_end": None},
        "historical_output_indicators_not_error_bounds": True}
    summary["long_horizon_statuses"]["NLSP_FEM3B_SOURCE_PRESERVATION"] = "PASS"
    save(bundle, item, summary)
    return bundle, item, summary


def load_reference(c, p=64):
    """Restore p64 or the independently saved p48 space, without a static solve."""
    if p not in (48, 64):
        raise ValueError("FEM-3B only permits saved p48/p64 references")
    _, _, _, science = load_parent(c)
    reference = base.one.load_reference(ROOT, science["source_static"]["bundle"],
        science["source_fem1"]["bundle"], science["source_action"]["bundle"])
    if p == 64:
        return reference
    source = ROOT / science["source_static"]["bundle"]
    path = base.one._checked_artifact(source, "one_d_p48.npz")
    row = next(r for r in reference["preflight"]["cases"] if r["p"] == p)
    if row["ndof"] != 4 * (p - 1) or row["nq"] != 2 * p + 1 or row["fields"] != list(base.one.FIELDS):
        raise ValueError("Saved p48 ordering, quadrature or dimension mismatch")
    old_disc = reference["disc"]
    disc = base.one.dynamics.PlanarGalerkin(old_disc.coefficients, p,
        length=old_disc.length, nq=row["nq"], model=old_disc.model, whiten=True)
    with np.load(path, allow_pickle=False) as data:
        saved = {name: data[name].copy() for name in ("s", "q_linear", "q_nonlinear",
            "raw_linear", "raw_nonlinear", "linear", "nonlinear")}
    for kind in ("linear", "nonlinear"):
        q = saved["q_" + kind]
        if q.shape != (disc.ndof,) or not np.isfinite(q).all():
            raise ValueError("Invalid frozen p48 state")
        if (np.max(abs(disc.raw_coefficients(q) - saved["raw_" + kind])) > 1e-12
            or np.max(abs(disc.reconstruct(q, saved["s"]) - saved[kind])) > 1e-12):
            raise ValueError("Frozen p48 raw/physical state reproduction failed")
        if np.max(abs(disc.reconstruct(q, [0., 1.]))) != 0.:
            raise ValueError("Frozen p48 essential clamps changed")
        np.linalg.cholesky(disc.mass_matrix(q))
        q.flags.writeable = False
        saved["raw_" + kind].flags.writeable = False
    reference.update(disc=disc, saved=saved)
    reference["source"] = {**reference["source"], "coordinates_sha256": sha(path),
        "p": p, "nq": disc.nq, "ndof": disc.ndof,
        "initial_policy": "reuse_saved_static_coordinates_exactly_no_projection"}
    return reference


def correction_decomposition(nonlinear, linear_from_linear_ic, linear_from_nonlinear_ic):
    nonlinear, linear, common = map(np.asarray, (nonlinear, linear_from_linear_ic, linear_from_nonlinear_ic))
    if nonlinear.shape != linear.shape or nonlinear.shape != common.shape or nonlinear.ndim < 2:
        raise ValueError("Physical trajectory shapes disagree")
    total = nonlinear - linear
    initial_state = common - linear
    same_ic = nonlinear - common
    identity_error = total - initial_state - same_ic
    return {"total_correction": total, "initial_state_component": initial_state,
        "same_ic_nonlinear_component": same_ic, "evolving_correction": total - total[0],
        "identity_max_abs": float(np.max(abs(identity_error))),
        "physical_initial_states_not_amplitude_aligned": True}


def choose_horizon(candidate_metrics, c):
    """Apply the previously stored planning heuristic, never a validation gate."""
    validate_config(c)
    if len(candidate_metrics) != 2:
        raise ValueError("Exactly the two predeclared horizon candidates required")
    rows = {float(r["horizon_T1"]): r for r in candidate_metrics}
    if set(rows) != {.25, .5}:
        raise ValueError("Exactly the two predeclared horizon candidates required")
    if not all(r.get("safety_passed") is True and r.get("resource_preflight_passed") is True for r in rows.values()):
        raise ValueError("Candidate safety or resource preflight incomplete")
    heuristic = c["horizon_selection"]
    decisions = {}
    for fraction, row in rows.items():
        signal = float(row["midspan_evolution_signal"])
        p_difference = row.get("p48_p64_evolution_difference")
        indicators = [heuristic["historical_DAT_FRD_displacement_difference"],
            heuristic["historical_41_81_correction_profile_difference"]]
        if p_difference is not None:
            indicators.append(float(p_difference))
        if signal < 0 or not math.isfinite(signal) or any(v < 0 or not math.isfinite(v) for v in indicators):
            raise ValueError("Invalid signal or uncertainty scale")
        limit = heuristic["planning_multiplier"] * max(indicators)
        decisions[str(fraction)] = {"midspan_evolution_signal": signal,
            "comparison_indicators": indicators, "planning_requirement": limit,
            "signal_exceeds_planning_requirement": signal >= limit,
            "p48_p64_available": p_difference is not None,
            "physical_validation_threshold": False}
    selected = .25 if decisions["0.25"]["signal_exceeds_planning_requirement"] else .5
    return {"horizon_T1": selected, "candidate_evidence": candidate_metrics,
        "planning_diagnostics": decisions, "selected_before_any_new_3D_result": True,
        "chosen_by_future_1D_3D_agreement": False,
        "exploratory_signal_resolution_qualification": not decisions[str(selected)]["signal_exceeds_planning_requirement"],
        "reason": "Preferred 0.25T1 exceeds each historical indicator and available p-sensitivity by factor10"
            if selected == .25 else "Preferred 0.25T1 insufficient by predeclared heuristic; use predeclared exploratory 0.5T1 candidate"}


def freeze_horizon(bundle, item, summary, candidate_metrics):
    if summary["job_calls"]["CCX_production"] or summary["attempts"]:
        raise ValueError("Cannot choose a horizon after native execution")
    path = Path(bundle) / "horizon_decision.json"
    if path.exists():
        raise ValueError("Pre-FEM horizon is already frozen")
    decision = choose_horizon(candidate_metrics, item["long_horizon_config"])
    horizon = decision["horizon_T1"]
    T = summary["preflight"]["T1"]
    c = copy.deepcopy(item["config"])
    c.update(horizon_T1=horizon, numerical_budget_seconds=6000,
        job_timeout_seconds=item["long_horizon_config"]["job_timeout_seconds_by_horizon"][str(horizon)])
    # The original science configuration remains the historical identity. This
    # distinct execution configuration contains only authorized scope changes.
    c["authorization"] = item["authorization"]
    parent, _, historical, _ = load_parent(item["long_horizon_config"])
    estimates = {kind: {"solver_seconds": historical["cases"][kind]["job"]["seconds"] * horizon / .05,
        "accepted_dynamic_increments": historical["cases"][kind]["dynamic_increments"] * horizon / .05}
        for kind in ("linear", "nonlinear")}
    decision.update(target_end=T * horizon, job_timeout_seconds=c["job_timeout_seconds"],
        frozen_dynamic_policy=c["dynamic"], created_before_native_attempts=0,
        estimated_execution=estimates, estimates_from_parent=str(parent.relative_to(ROOT)),
        estimate_qualification="Runtime scaling is a planning estimate, not a guarantee")
    write_json(path, decision)
    write_json(Path(bundle) / "frozen_3d_config.json", c)
    summary["horizon_decision_sha256"] = sha(path)
    summary["selected_horizon_T1"] = horizon
    summary["preflight"]["target_dynamic_end"] = T * horizon
    summary["long_horizon_statuses"]["NLSP_FEM3B_HORIZON_SELECTION"] = "PASS"
    save(bundle, item, summary)
    return decision


def save_dense_records(path, records):
    if not records:
        raise ValueError("No accepted Radau dense records")
    payload = {name: np.array([record[name] for record in records])
        for name in ("t_old", "t", "y_old", "Q")}
    if payload["Q"].ndim != 3 or payload["Q"].shape[-1] != 3:
        raise ValueError("Expected the existing Radau cubic dense coefficients")
    if any(not np.isfinite(a).all() for a in payload.values()):
        raise ValueError("Nonfinite accepted Radau output")
    if (payload["t_old"][0] != 0. or np.any(payload["t"] <= payload["t_old"])
        or not np.array_equal(payload["t_old"][1:], payload["t"][:-1])):
        raise ValueError("Accepted Radau intervals are not a contiguous actual prefix")
    np.savez_compressed(path, **payload)
    return {"accepted_steps": len(records), "actual_end": float(payload["t"][-1]),
        "formula": "y_old+Q@[x,x^2,x^3];x=(t-t_old)/(t_end-t_old);no_extra_h_factor",
        "new_integrations": 0, "no_trajectory_interpolation": True}


def evaluate_dense_records(path, times):
    """Evaluate the already accepted Radau polynomial at new physical times."""
    times = np.asarray(times, dtype=float)
    with np.load(path, allow_pickle=False) as z:
        old, ends, initial, coefficients = (z[k] for k in ("t_old", "t", "y_old", "Q"))
        if (times.ndim != 1 or not np.isfinite(times).all() or np.any(np.diff(times) < 0)
            or np.any(times < old[0]) or np.any(times > ends[-1])):
            raise ValueError("Requested time is outside the accepted actual Radau prefix")
        indices = np.searchsorted(ends, times, side="left")
        x = (times - old[indices]) / (ends[indices] - old[indices])
        powers = np.stack((x, x*x, x*x*x), axis=-1)
        return initial[indices] + np.einsum("nij,nj->ni", coefficients[indices], powers)


def pair_native_histories(times_linear, fields_linear, times_nonlinear, fields_nonlinear, *, common_grid=None):
    """Pair actual equal samples; explicitly diagnose optional interpolation."""
    tl, tn = np.asarray(times_linear), np.asarray(times_nonlinear)
    fl, fn = np.asarray(fields_linear), np.asarray(fields_nonlinear)
    for t, f in ((tl, fl), (tn, fn)):
        if t.ndim != 1 or len(t) < 2 or not np.isfinite(t).all() or np.any(np.diff(t) <= 0) or f.shape[0] != len(t) or not np.isfinite(f).all():
            raise ValueError("Invalid actual native trajectory")
    if fl.shape[1:] != fn.shape[1:]:
        raise ValueError("Native section coordinates or fields disagree")
    shared, il, inn = np.intersect1d(tl, tn, return_indices=True)
    if np.array_equal(tl, tn):
        return {"time": tl.copy(), "linear": fl.copy(), "nonlinear": fn.copy(),
            "correction": fn-fl, "sample_origin": "native", "interpolation_used": False,
            "interpolation_correction_difference_max_abs": 0.}
    if common_grid is None:
        raise ValueError("Divergent native schedules require the predeclared explicit common grid")
    grid = np.asarray(common_grid, dtype=float)
    low, high = max(tl[0], tn[0]), min(tl[-1], tn[-1])
    if (grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all() or np.any(np.diff(grid) <= 0)
        or grid[0] < low or grid[-1] > high):
        raise ValueError("Interpolation grid must lie inside actual native overlap")
    from scipy.interpolate import PchipInterpolator
    def linear(t, f):
        flat = f.reshape(len(t), -1)
        return np.stack([np.interp(grid, t, column) for column in flat.T], axis=-1).reshape((len(grid),) + f.shape[1:])
    l_linear, n_linear = linear(tl, fl), linear(tn, fn)
    l_cubic, n_cubic = PchipInterpolator(tl, fl, axis=0)(grid), PchipInterpolator(tn, fn, axis=0)(grid)
    delta_linear, delta_cubic = n_linear-l_linear, n_cubic-l_cubic
    return {"time": grid, "linear": l_linear, "nonlinear": n_linear,
        "correction": delta_linear, "linear_PCHIP": l_cubic, "nonlinear_PCHIP": n_cubic,
        "correction_PCHIP": delta_cubic, "sample_origin": "interpolated",
        "interpolation_used": True, "primary_interpolation": "piecewise_linear",
        "diagnostic_interpolation": "shape_preserving_cubic_PCHIP",
        "interpolation_correction_difference_max_abs": float(np.max(abs(delta_cubic-delta_linear))),
        "native_values_not_relabelled": True, "no_extrapolation": True,
        "exact_shared_native_time": shared, "exact_shared_native_correction": fn[inn]-fl[il]}


def run_3d_case(bundle, item, summary, kind):
    """Only the native-job ledger is new; generation and recovery are reused."""
    bundle = Path(bundle)
    if kind not in ("linear", "nonlinear"):
        raise ValueError("Only the authorized linear/nonlinear pair is available")
    if summary.get("hard_stop") or summary["overall"] == "BLOCKED_BY_SOLVER":
        return False
    if not (bundle / "horizon_decision.json").exists():
        raise ValueError("Freeze the pre-FEM horizon before any native job")
    if sha(bundle / "horizon_decision.json") != summary["horizon_decision_sha256"]:
        raise ValueError("Pre-FEM horizon decision changed")
    if kind == "nonlinear" and (summary["cases"].get("linear", {}).get("status") != "PASS"
        or summary["cases"]["linear"].get("continuation_audit", {}).get("status") != "PASS"):
        raise ValueError("Actual linear execution and audit must precede nonlinear")
    science = read_json(bundle / "frozen_3d_config.json")
    previous = summary["cases"].get(kind)
    if previous and previous["status"] == "PASS":
        return True
    if previous and previous["status"] != "OUTPUT_RECOVERY_PENDING":
        return False
    old, source, mesh, audit = base.verify_sources(science)
    case = bundle / "cases" / kind
    if not previous:
        if summary["job_calls"]["CCX_production"] >= 2:
            raise RuntimeError("Two explicitly authorized native attempts exhausted")
        case.mkdir(parents=True)
        contract = base.write_input(case / "motion.inp", science, old, source, mesh, audit, kind == "nonlinear")
        resume.output_safety((case / "motion.inp").read_text(encoding="utf8"))
        write_json(case / "input_contract.json", contract)
        summary["attempts"].append({"kind": kind, "ordinal": len(summary["attempts"])+1,
            "status": "STARTED", "input_sha256": sha(case / "motion.inp"),
            "authorization_id": AUTHORIZATION, "horizon_decision_sha256": summary["horizon_decision_sha256"]})
        summary["cases"][kind] = {"status": "STARTED"}
        remaining = 6000 - summary["numerical_seconds"]
        if remaining <= 0:
            summary["hard_stop"] = True
            summary["cases"][kind].update(status="NOT_RUN", failure="Predeclared numerical budget exhausted before job")
            save(bundle, item, summary)
            return False
        env = dict(os.environ)
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMBER_OF_CPUS="1")
        summary["job_calls"]["CCX_production"] += 1
        save(bundle, item, summary)
        started = time.perf_counter()
        try:
            print("FEM3B production static+dynamic job: " + kind, flush=True)
            command = [old["science_config"]["ccx_exe"], "motion"]
            write_json(case / "execution_environment.json", {"command": command,
                "cwd": str(case), "threads": 1, "PATH_unchanged": True,
                "thread_environment": {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS")}})
            result, stats = base.base.fem1.run_job(command, case,
                min(science["job_timeout_seconds"], remaining), 4*1024**3, case / "motion", env)
            write_json(case / "job.json", stats)
            summary["cases"][kind]["job"] = stats
            stdout = (case / "motion.stdout.txt").read_text(encoding="utf8", errors="replace")
            if result.returncode != 0 or stats["failure"] or "JOB FINISHED" not in stdout.upper() or "*ERROR" in stdout.upper():
                raise RuntimeError("Actual solver failure: " + str(stats))
            summary["cases"][kind]["status"] = "OUTPUT_RECOVERY_PENDING"
            summary["attempts"][-1]["status"] = "SOLVER_FINISHED"
        except Exception as error:
            summary["cases"][kind].update(status="FAIL", failure=str(error))
            summary["attempts"][-1].update(status="FAIL", failure=str(error))
            summary["hard_stop"] = True
            summary["overall"] = "BLOCKED_BY_SOLVER"
        finally:
            summary["numerical_seconds"] += time.perf_counter() - started
            save(bundle, item, summary)
        if summary["cases"][kind]["status"] == "FAIL":
            summary["long_horizon_statuses"]["NLSP_FEM3B_" + kind.upper() + "_3D"] = "FAIL"
            save(bundle, item, summary)
            return False
    try:
        job = summary["cases"][kind]["job"]
        record = base.recover_case(bundle, science, summary, kind, old, mesh, audit)
        record["job"] = job
        summary["cases"][kind] = record
        resume.audit_case(bundle, science, summary, kind, native_time_rounding=True)
        next(a for a in summary["attempts"] if a["kind"] == kind)["status"] = "PASS"
        summary["long_horizon_statuses"]["NLSP_FEM3B_" + kind.upper() + "_3D"] = "PASS"
        summary["long_horizon_statuses"]["NLSP_FEM3B_ENERGY_DIAGNOSTICS"] = "PARTIAL"
        summary["overall"] = "PARTIAL"
        save(bundle, item, summary)
        return True
    except Exception as error:
        summary["cases"][kind].update(status="OUTPUT_RECOVERY_PENDING", recovery_failure=str(error))
        summary["overall"] = "PARTIAL"
        save(bundle, item, summary)
        return False


def _save_trajectory(path, reference, history, times, *, linear, execution):
    """Measure existing trajectories in blocks; retain their own energy origin."""
    disc = reference["disc"]
    points = np.linspace(0., 1., 41)
    count = len(history)
    fields = np.empty((count, len(points), 4))
    velocities = np.empty_like(fields)
    energy = np.empty(count)
    norms, speed_norms = np.empty((count, 4)), np.empty((count, 4))
    safety = {}
    for start in range(0, count, 256):
        stop = min(count, start + 256)
        measured = base.one.summarize_reference(reference,
            {"times": times[start:stop], "q": history[start:stop, :disc.ndof],
             "velocity": history[start:stop, disc.ndof:]}, points, linear=linear)
        fields[start:stop], velocities[start:stop] = measured["fields"], measured["physical_velocities"]
        energy[start:stop] = measured["energy"]
        norms[start:stop], speed_norms[start:stop] = measured["L2_fields"], measured["L2_velocities"]
        for row in measured["diagnostics"]["safety_samples"]:
            for name, value in row.items():
                if isinstance(value, bool) or not isinstance(value, (int, float)):
                    continue
                reducer = min if name in ("min_one_plus_c", "relative_mass_lower_bound") else max
                safety[name] = value if name not in safety else reducer(safety[name], value)
    drift = (energy-energy[0])/energy[0]
    payload = {"times": np.asarray(times), "q": history[:, :disc.ndof],
        "velocity": history[:, disc.ndof:], "x": points, "fields": fields,
        "physical_velocities": velocities, "energy": energy, "energy_relative_drift": drift,
        "L2_fields": norms, "L2_velocities": speed_norms}
    np.savez_compressed(path, **payload)
    diagnostics = {"sampled_safety": safety, "max_relative_energy_drift": float(np.max(abs(drift))),
        "energy_definition": "0.5*v.T*M0*v+0.5*q.T*K*q" if linear else "0.5*v.T*M(q)*v+V4(q)",
        "removed_GRAV_potential_included": False, "temporal_convergence_claimed": False,
        "spatial_convergence_claimed": False, "sampled_maxima_only": True}
    write_json(Path(path).with_suffix(".json"), {"execution": execution,
        "diagnostics": diagnostics, "source": reference["source"]})
    return diagnostics


def _working_history(bundle, name, count, ndof):
    path = Path(bundle) / (name + "_working.npy")
    history = np.lib.format.open_memmap(path, mode="w+", dtype=np.float64, shape=(count, 2*ndof))
    return path, history


def _linear_history(bundle, name, reference, times, initial_kind):
    disc = reference["disc"]
    q0 = reference["saved"]["q_" + initial_kind]
    path, history = _working_history(bundle, name, len(times), disc.ndof)
    roundoff = 0.
    for start in range(0, len(times), 256):
        stop = min(len(times), start + 256)
        current = disc.linear_reference(q0, np.zeros(disc.ndof), times[start:stop])
        if start == 0 and times[0] == 0.:
            roundoff = float(np.max(abs(current["q"][0]-q0)))
            current["q"][0] = q0
            current["velocity"][0] = 0.
        history[start:stop, :disc.ndof] = current["q"]
        history[start:stop, disc.ndof:] = current["velocity"]
    execution = {"exact_in_time": True, "ODE_integrations": 0,
        "semidiscrete_spectral_factorization": "all independent Shen coordinates; no modal reduction",
        "initial_kind": initial_kind, "zero_time_factorization_roundoff_max_abs": roundoff,
        "linear_eigendecompositions": disc.linear_eigendecompositions,
        "authorization_id": AUTHORIZATION}
    _save_trajectory(Path(bundle) / (name + ".npz"), reference, history, times,
        linear=True, execution=execution)
    del history
    path.unlink()


def _safety_passed(diagnostics, policy):
    s = diagnostics["sampled_safety"]
    checks = (s["min_one_plus_c"] >= policy["min_one_plus_c"],
        s["max_abs_c"] <= policy["max_abs_c"], s["max_abs_theta"] <= policy["max_abs_theta"],
        s["max_abs_u_s"] <= policy["max_abs_axial_gradient"],
        s["max_abs_w_s"] <= policy["max_abs_transverse_gradient"],
        s["max_L_abs_theta_s"] <= policy["max_L_abs_curvature"],
        s["relative_mass_lower_bound"] >= policy["min_relative_mass_eigenvalue"])
    return bool(all(checks))


def _complete_degree_references(bundle, prefix, reference, times):
    """Finish missing read-only artifacts after a saved successful NL solve."""
    bundle = Path(bundle)
    mode_path = bundle / ("linear_modes_p" + str(reference["disc"].p) + ".npz")
    if not mode_path.exists():
        raise ValueError("Saved complete linear factors missing; no eigenanalysis retry")
    with np.load(mode_path, allow_pickle=False) as modes:
        reference["disc"]._linear_modes[None] = {k: modes[k].copy() for k in modes.files}
    for suffix, initial_kind in (("linear", "linear"), ("linear_nonlinear_initial", "nonlinear")):
        name = prefix + "_" + suffix
        if not (bundle / (name + ".npz")).exists() or not (bundle / (name + ".json")).exists():
            _linear_history(bundle, name, reference, times, initial_kind)
    path = bundle / (prefix + "_decomposition.npz")
    if path.exists() and path.with_suffix(".json").exists():
        return
    with np.load(bundle / (prefix + "_nonlinear.npz")) as nonlinear, \
         np.load(bundle / (prefix + "_linear.npz")) as linear, \
         np.load(bundle / (prefix + "_linear_nonlinear_initial.npz")) as common:
        decomposition = correction_decomposition(nonlinear["fields"], linear["fields"], common["fields"])
        np.savez_compressed(path, times=times, x=nonlinear["x"],
            **{k: v for k, v in decomposition.items() if isinstance(v, np.ndarray)},
            identity_max_abs=decomposition["identity_max_abs"])
        write_json(path.with_suffix(".json"), {"identity_max_abs": decomposition["identity_max_abs"],
            "same_nonlinear_initial_state_used_for_auxiliary_linear_reference": True,
            "no_additional_nonlinear_integration": True,
            "saved_complete_linear_factors_reused": True})


def run_1d_stage(bundle, item, summary):
    """Two bounded independent-degree runs, each ledgered before execution."""
    bundle = Path(bundle)
    if summary.get("completed") or summary.get("hard_stop"):
        return summary
    if summary["long_horizon_statuses"]["NLSP_FEM3B_OLD_SIGNAL_DIAGNOSTIC"] != "PASS":
        raise ValueError("Read-only old-signal diagnostic must precede 1D integration")
    if (bundle / "candidate_signals.json").exists():
        return summary
    if any(row["status"] != "PASS" for row in summary["one_d_attempts"]):
        summary["hard_stop"] = True
        save(bundle, item, summary)
        return summary
    stage_started = time.perf_counter()
    numerical_at_entry = summary["numerical_seconds"]
    c, science = item["long_horizon_config"], item["config"]
    snapshot = bundle / "execution_code" / "fem3b_stage_b_implementation.py"
    if not snapshot.exists():
        shutil.copyfile(__file__, snapshot)
        write_json(bundle / "stage_b_implementation.json", {"sha256": sha(snapshot),
            "preparation_identity_preserved": True, "snapshot_before_1D_calls": True})
        save(bundle, item, summary)
    refs = {64: load_reference(c, 64)}
    try:
        refs[48] = load_reference(c, 48)
    except (FileNotFoundError, StopIteration) as error:
        summary["p48_spatial_control"] = {"status": "NOT_RUN", "reason": str(error)}
    # Full semidiscrete eigensystems are required for exact-time references;
    # their highest retained frequency controls output sampling, not timestep.
    omega_max = 0.
    for p, reference in refs.items():
        modes_path = bundle / ("linear_modes_p" + str(p) + ".npz")
        if modes_path.exists():
            with np.load(modes_path) as cached_modes:
                reference["disc"]._linear_modes[None] = {k: cached_modes[k].copy() for k in cached_modes.files}
        modes = reference["disc"].linear_eigenpairs()
        omega_max = max(omega_max, float(modes["omega"].max()))
        if not modes_path.exists():
            np.savez_compressed(modes_path, **modes)
    T = summary["preflight"]["T1"]
    horizon = .5*T
    count = int(np.ceil(horizon*omega_max/(2*np.pi)*12))
    parent = ROOT / c["parent_completed"]["bundle"]
    with np.load(parent / "one_d_nonlinear.npz") as old:
        old_times = old["times"].copy()
    times = np.unique(np.r_[np.linspace(0., horizon, count+1), 0., .25*T, horizon, old_times])
    np.save(bundle / "one_d_common_times.npy", times)
    write_json(bundle / "one_d_sampling.json", {"samples": len(times),
        "omega_max_full_retained_spectrum": omega_max, "samples_per_highest_retained_period": 12,
        "old_timestamps_preserved": True, "required_exact_times": [0., .25*T, horizon],
        "no_output_frequency_cap_or_filter": True, "output_sampling_not_time_accuracy_control": True})
    for p in (64, 48):
        if p not in refs:
            continue
        reference, prefix = refs[p], "one_d_p" + str(p)
        disc = reference["disc"]
        existing = next((a for a in summary["one_d_attempts"] if a["p"] == p), None)
        if existing:
            if existing["status"] != "PASS":
                summary["hard_stop"] = True
                save(bundle, item, summary)
                return summary
            _complete_degree_references(bundle, prefix, reference, times)
            continue
        if summary["job_calls"]["1D_nonlinear_ODE"] >= 2:
            raise RuntimeError("Two FEM-3B 1D nonlinear attempts exhausted")
        config = base.one.runtime_config(reference, science)
        config["spatial"]["degrees"] = [p]
        row = {"p": p, "status": "STARTED", "authorization_id": AUTHORIZATION,
            "source_coordinates_sha256": reference["source"]["coordinates_sha256"],
            "target_end": horizon, "automatic_retry": False}
        summary["one_d_attempts"].append(row)
        summary["job_calls"]["1D_nonlinear_ODE"] += 1
        save(bundle, item, summary)
        path, history_buffer = _working_history(bundle, prefix + "_nonlinear", len(times), disc.ndof)
        dense_records = []
        started = time.perf_counter()
        base.one.runner.load_runtime()
        try:
            history, stats = base.one.runner.integrate_case(disc, None,
                {"omega": reference["omega1"], "T1": T}, config,
                reference["preflight"]["load"]["linear_w_max"] / .1, "tight", times,
                time.perf_counter() + max(0., 6000-summary["numerical_seconds"]),
                initial_coordinates=reference["saved"]["q_nonlinear"],
                history_buffer=history_buffer, dense_output_observer=dense_records.append)
            stats.update(authorization_id=AUTHORIZATION, execution_mode="EXPLORATORY_NOT_CERTIFIED",
                admitted=False, initial_coordinates_reused_exactly=True,
                no_dynamic_derivative_constraints=True, external_force_after_release=0.,
                strict_float64_strong_weak="PARTIAL", new_ODE_integrations=1)
            actual = times[:len(history)]
            dense_meta = save_dense_records(bundle / ("dense_p" + str(p) + ".npz"), dense_records)
            check_indices = np.unique(np.linspace(0, len(actual)-1, min(257, len(actual)), dtype=int))
            dense_check = evaluate_dense_records(bundle / ("dense_p" + str(p) + ".npz"), actual[check_indices])
            stats["saved_dense_reproduction_max_abs"] = float(np.max(abs(dense_check-history[check_indices])))
            stats["saved_dense_metadata"] = dense_meta
            diagnostics = _save_trajectory(bundle / (prefix + "_nonlinear.npz"), reference,
                history, actual, linear=False, execution=stats)
            row.update(status=stats["status"], actual_end=float(actual[-1]),
                safety_passed=_safety_passed(diagnostics, config["safety"]),
                max_relative_energy_drift=diagnostics["max_relative_energy_drift"])
            if stats["status"] != "PASS" or not row["safety_passed"]:
                summary["hard_stop"] = True
                summary["overall"] = "PARTIAL"
        except Exception as error:
            row.update(status="FAIL", failure=str(error))
            summary["hard_stop"] = True
            summary["overall"] = "PARTIAL"
            raise
        finally:
            summary["numerical_seconds"] += time.perf_counter()-started
            save(bundle, item, summary)
        del history, history_buffer, dense_records
        path.unlink()
        if summary.get("hard_stop"):
            return summary
        _complete_degree_references(bundle, prefix, reference, times)
        save(bundle, item, summary)
    # Historical all-eight norms retain the old gates and reporting floor.
    if 48 in refs:
        compare_config = copy.deepcopy(science)
        compare_config["gates"] = read_json(base.one.runner.CONFIG)["gates"]
        with np.load(bundle / "one_d_p48_nonlinear.npz") as a, \
             np.load(bundle / "one_d_p64_nonlinear.npz") as b:
            spatial = base.one.runner.compare_histories(refs[48]["disc"], a, refs[64]["disc"], b, compare_config)
        write_json(bundle / "one_d_all8_spatial.json", spatial)
        summary["p48_spatial_control"] = {"status": spatial["status"], "artifact": "one_d_all8_spatial.json"}
    candidates = []
    summary["numerical_seconds"] = max(summary["numerical_seconds"],
        numerical_at_entry+time.perf_counter()-stage_started)
    parent_summary = read_json(ROOT / c["parent_completed"]["bundle"] / "summary.json")
    with np.load(bundle / "one_d_p64_decomposition.npz") as main:
        evolution = main["evolving_correction"][:, :, 1]
        initial_correction = main["total_correction"][0, :, 1]
        for fraction in (.25, .5):
            stop = np.searchsorted(times, fraction*T, side="right")
            p_difference = None
            if 48 in refs:
                with np.load(bundle / "one_d_p48_decomposition.npz") as other:
                    p_difference = float(np.max(abs(evolution[:stop]-other["evolving_correction"][:stop, :, 1])))
            nonlinear_meta = read_json(bundle / "one_d_p64_nonlinear.json")
            estimated_ccx = sum(parent_summary["cases"][kind]["job"]["seconds"] for kind in ("linear", "nonlinear"))*fraction/.05
            candidates.append({"horizon_T1": fraction,
                "midspan_evolution_signal": float(np.max(abs(evolution[:stop, 20]))),
                "spatial_evolution_signal": float(np.max(abs(evolution[:stop]))),
                "max_time_L2_evolution": float(np.max(np.sqrt(np.trapezoid(evolution[:stop]**2, main["x"], axis=1)))),
                "initial_correction_common_scale": float(np.max(abs(initial_correction))),
                "p48_p64_evolution_difference": p_difference,
                "p_sensitivity_policy": "conservative entire41-section evolving-w max on same prefix; fixed before3D",
                "safety_passed": _safety_passed(nonlinear_meta["diagnostics"], science["safety"]),
                "resource_preflight_passed": shutil.disk_usage(bundle).free > 4*1024**3 and estimated_ccx < 6000-summary["numerical_seconds"],
                "estimated_total_CCX_seconds_from_parent": estimated_ccx,
                "available_disk_bytes": shutil.disk_usage(bundle).free,
                "remaining_numerical_budget_seconds": 6000-summary["numerical_seconds"]})
    write_json(bundle / "candidate_signals.json", candidates)
    summary["long_horizon_statuses"]["NLSP_FEM3B_1D_PRELIMINARY"] = "PASS"
    save(bundle, item, summary)
    freeze_horizon(bundle, item, summary, candidates)
    return summary


def evaluate_saved_one_d(bundle, item, times):
    """Saved dense NL and complete linear factors; no new ODE or eigensolve."""
    bundle = Path(bundle)
    times = np.asarray(times, dtype=float)
    if times.ndim != 1 or times[0] != 0. or np.any(np.diff(times) <= 0):
        raise ValueError("Actual comparison times must increase from zero")
    reference = load_reference(item["long_horizon_config"], 64)
    disc = reference["disc"]
    with np.load(bundle / "linear_modes_p64.npz") as saved_modes:
        disc._linear_modes[None] = {k: saved_modes[k].copy() for k in saved_modes.files}
    qv = evaluate_dense_records(bundle / "dense_p64.npz", times)
    points = np.linspace(0., 1., 41)
    linear_fields, linear_velocities = [], []
    nonlinear_fields, nonlinear_velocities = [], []
    for start in range(0, len(times), 256):
        stop = min(len(times), start+256)
        linear = disc.linear_reference(reference["saved"]["q_linear"], np.zeros(disc.ndof), times[start:stop])
        if start == 0:
            linear["q"][0] = reference["saved"]["q_linear"]
            linear["velocity"][0] = 0.
        linear_fields.append(disc.reconstruct_series(linear["q"], points))
        linear_velocities.append(disc.reconstruct_series(linear["velocity"], points))
        nonlinear_fields.append(disc.reconstruct_series(qv[start:stop, :disc.ndof], points))
        nonlinear_velocities.append(disc.reconstruct_series(qv[start:stop, disc.ndof:], points))
    if disc.linear_eigendecompositions:
        raise ValueError("Saved linear factorization was not restored; new eigensolve forbidden")
    return {"times": times, "x": points, "linear_fields": np.concatenate(linear_fields),
        "nonlinear_fields": np.concatenate(nonlinear_fields),
        "linear_velocities": np.concatenate(linear_velocities),
        "nonlinear_velocities": np.concatenate(nonlinear_velocities),
        "initial_linear_fields": disc.reconstruct(reference["saved"]["q_linear"], points),
        "initial_nonlinear_fields": disc.reconstruct(reference["saved"]["q_nonlinear"], points)}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--run-3d", action="store_true")
    mode.add_argument("--run-1d", action="store_true")
    mode.add_argument("--compute", action="store_true")
    mode.add_argument("--report-only", type=Path)
    mode.add_argument("--plot-only", type=Path)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--through-case", choices=("linear", "nonlinear"), default="nonlinear")
    args = parser.parse_args(argv)
    if args.report_only or args.plot_only:
        bundle = args.report_only or args.plot_only
        summary = validate_cache(bundle)
        if args.plot_only:
            # Stage-specific cached rendering is registered by the final report
            # module. It may never call the old .05T1 reference computation.
            from scripts.lib import nlsp_fem3b_diagnostics
            nlsp_fem3b_diagnostics.plot_bundle(bundle)
        print(json.dumps({"bundle": str(bundle), "overall": summary["overall"],
            "statuses": summary["long_horizon_statuses"], "new_scientific_calls": 0}, indent=2))
        return summary
    bundle, item, summary = prepare_stage(args.config)
    calls_before = dict(summary["job_calls"])
    if summary.get("completed") or summary.get("hard_stop"):
        print(json.dumps({"bundle": str(bundle), "overall": summary["overall"], "new_scientific_calls": 0}, indent=2))
        return summary
    if args.run_1d or args.compute:
        run_1d_stage(bundle, item, summary)
    if (args.run_3d or args.compute) and not summary.get("hard_stop"):
        for kind in ("linear", "nonlinear"):
            if not run_3d_case(bundle, item, summary, kind) or args.through_case == kind:
                break
    if args.compute and all(summary["cases"].get(k, {}).get("status") == "PASS" for k in ("linear", "nonlinear")):
        from scripts.lib import nlsp_fem3b_diagnostics
        started = time.perf_counter()
        times = nlsp_fem3b_diagnostics.actual_pairing_times(bundle, item)
        np.savez_compressed(bundle / "comparison_one_d.npz", **evaluate_saved_one_d(bundle, item, times))
        nlsp_fem3b_diagnostics.complete_comparison(bundle, item, summary)
        summary["numerical_seconds"] += time.perf_counter()-started
        summary["completed"] = True
        summary["overall"] = "FEM3B_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
        save(bundle, item, summary)
        nlsp_fem3b_diagnostics.plot_bundle(bundle)
        save(bundle, item, summary)
    print(json.dumps({"bundle": str(bundle), "overall": summary["overall"],
        "statuses": summary["long_horizon_statuses"], "job_calls": summary["job_calls"],
        "new_scientific_calls": {k: value-calls_before.get(k, 0) for k, value in summary["job_calls"].items()}}, indent=2))
    return summary
