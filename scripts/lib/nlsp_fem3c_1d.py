"""Bounded FEM-3C full-period references using the unchanged planar runner.

Internal helper, never a second physics solver. Static coordinates and complete
linear factors are restored from FEM-3B; each nonlinear call is ledgered once.
"""
from __future__ import annotations

import copy
import time
from pathlib import Path

import numpy as np

from scripts.lib import nlsp_fem3b_continuation as previous


def _reference(item, p):
    parent = previous.ROOT / item["validation_config"]["parent_completed"]["bundle"]
    old_item = previous.read_json(parent / "provenance.json")
    reference = previous.load_reference(old_item["long_horizon_config"], p)
    mode_path = parent / f"linear_modes_p{p}.npz"
    with np.load(mode_path, allow_pickle=False) as saved:
        modes = {k: saved[k].copy() for k in saved.files}
    n = reference["disc"].ndof
    if (modes.get("vectors", np.empty(0)).shape != (n, n)
            or any(modes.get(k, np.empty(0)).shape != (n,)
                   for k in ("omega", "eigenvalues"))
            or any(not np.isfinite(v).all() for v in modes.values())
            or np.any(modes["omega"] <= 0) or np.any(modes["eigenvalues"] <= 0)):
        raise ValueError("Saved complete linear factors are invalid or truncated")
    reference["disc"]._linear_modes[None] = modes
    if reference["disc"].linear_eigendecompositions:
        raise ValueError("FEM-3C restores complete saved factors without eigenanalysis")
    return reference, parent, modes


def _linear(bundle, name, reference, times, initial_kind, authorization):
    disc = reference["disc"]
    path, history = previous._working_history(bundle, name, len(times), disc.ndof)
    q0 = reference["saved"]["q_" + initial_kind]
    for start in range(0, len(times), 256):
        stop = min(start + 256, len(times))
        current = disc.linear_reference(q0, np.zeros(disc.ndof), times[start:stop])
        if start == 0:
            current["q"][0], current["velocity"][0] = q0, 0.
        history[start:stop, :disc.ndof] = current["q"]
        history[start:stop, disc.ndof:] = current["velocity"]
    previous._save_trajectory(bundle / (name + ".npz"), reference, history,
        times, linear=True, execution={"authorization_id": authorization,
        "exact_in_time": True, "ODE_integrations": 0,
        "saved_complete_linear_factors_reused": True,
        "no_modal_truncation": True, "initial_kind": initial_kind})
    del history
    path.unlink()


def _complete_read_only(bundle, reference, times, authorization):
    prefix = f"one_d_p{reference['disc'].p}"
    for suffix, initial in (("linear", "linear"), ("linear_nonlinear_initial", "nonlinear")):
        if not (bundle / f"{prefix}_{suffix}.npz").exists():
            _linear(bundle, f"{prefix}_{suffix}", reference, times, initial, authorization)
    with np.load(bundle / f"{prefix}_nonlinear.npz") as nonlinear, \
         np.load(bundle / f"{prefix}_linear.npz") as linear, \
         np.load(bundle / f"{prefix}_linear_nonlinear_initial.npz") as common:
        decomposition = previous.correction_decomposition(nonlinear["fields"],
            linear["fields"], common["fields"])
        np.savez_compressed(bundle / f"{prefix}_decomposition.npz", times=times,
            x=nonlinear["x"], **{k: v for k, v in decomposition.items()
            if isinstance(v, np.ndarray)}, identity_max_abs=decomposition["identity_max_abs"])
        previous.write_json(bundle / f"{prefix}_decomposition.json", {
            "identity_max_abs": decomposition["identity_max_abs"],
            "same_nonlinear_initial_state_used_for_auxiliary_linear_reference": True,
            "additional_nonlinear_integration": False})


def run_full_period(bundle, item, summary):
    """At most one p64 and one p48 nonlinear solve, each from saved q0."""
    from scripts.lib import nlsp_fem3c_validation as stage
    bundle = Path(bundle)
    if summary.get("one_d_completed") or summary.get("hard_stop"):
        return summary
    c, science = item["validation_config"], item["config"]
    authorization = c["authorization"]["id"]
    refs = {p: _reference(item, p)[0] for p in (64, 48)}
    T = 2 * np.pi / science["omega1"]
    omega_max = max(float(r["disc"]._linear_modes[None]["omega"].max()) for r in refs.values())
    count = int(np.ceil(T * omega_max / (2 * np.pi) * 12))
    parent = previous.ROOT / c["parent_completed"]["bundle"]
    old_times = np.load(parent / "one_d_common_times.npy")
    times = np.unique(np.r_[np.linspace(0., T, count + 1), old_times,
        0., .05*T, .25*T, .5*T, .75*T, T])
    np.save(bundle / "one_d_common_times.npy", times)
    previous.write_json(bundle / "one_d_sampling.json", {
        "target_end": T, "samples": len(times),
        "omega_max_full_retained_spectrum": omega_max,
        "samples_per_highest_retained_period": 12, "old_timestamps_preserved": True,
        "required_exact_times": [0., .05*T, .25*T, .5*T, .75*T, T],
        "sampling_is_not_temporal_certification": True})
    summary.setdefault("one_d_attempts", [])
    for p in (64, 48):
        reference, prefix = refs[p], f"one_d_p{p}"
        disc = reference["disc"]
        existing = next((r for r in summary["one_d_attempts"] if r["p"] == p), None)
        if existing:
            if existing["status"] != "PASS":
                summary["hard_stop"] = True
                stage.save(bundle, item, summary)
                return summary
            _complete_read_only(bundle, reference, times, authorization)
            continue
        if summary["job_calls"]["1D_nonlinear_ODE"] >= 2:
            raise RuntimeError("FEM-3C nonlinear 1D attempt limit exhausted")
        runtime = previous.base.one.runtime_config(reference, science)
        runtime["spatial"]["degrees"] = [p]
        row = {"p": p, "status": "STARTED", "authorization_id": authorization,
            "target_end": T, "automatic_retry": False,
            "source_coordinates_sha256": reference["source"]["coordinates_sha256"]}
        summary["one_d_attempts"].append(row)
        summary["job_calls"]["1D_nonlinear_ODE"] += 1
        np.savez_compressed(bundle / f"linear_modes_p{p}.npz", **disc._linear_modes[None])
        stage.save(bundle, item, summary)
        path, buffer = previous._working_history(bundle, prefix + "_nonlinear", len(times), disc.ndof)
        records = []
        started = time.perf_counter()
        previous.base.one.runner.load_runtime()
        try:
            history, stats = previous.base.one.runner.integrate_case(disc, None,
                {"omega": reference["omega1"], "T1": T}, runtime,
                reference["preflight"]["load"]["linear_w_max"] / .1, "tight", times,
                time.perf_counter() + max(0., c["numerical_budget_seconds"] - summary["numerical_seconds"]),
                initial_coordinates=reference["saved"]["q_nonlinear"],
                history_buffer=buffer, dense_output_observer=records.append)
            actual = times[:len(history)]
            stats.update(authorization_id=authorization, execution_mode="EXPLORATORY_NOT_CERTIFIED",
                admitted=False, initial_coordinates_reused_exactly=True,
                no_dynamic_derivative_constraints=True, strict_float64_strong_weak="PARTIAL",
                external_force_after_release=0., new_ODE_integrations=1)
            stats["saved_dense_metadata"] = previous.save_dense_records(bundle / f"dense_p{p}.npz", records)
            indices = np.unique(np.linspace(0, len(actual)-1, min(257, len(actual)), dtype=int))
            check = previous.evaluate_dense_records(bundle / f"dense_p{p}.npz", actual[indices])
            stats["saved_dense_reproduction_max_abs"] = float(np.max(abs(check-history[indices])))
            diagnostics = previous._save_trajectory(bundle / f"{prefix}_nonlinear.npz",
                reference, history, actual, linear=False, execution=stats)
            row.update(status=stats["status"], actual_end=float(actual[-1]),
                safety_passed=previous._safety_passed(diagnostics, runtime["safety"]),
                max_relative_energy_drift=diagnostics["max_relative_energy_drift"])
            if stats["status"] != "PASS" or not row["safety_passed"]:
                summary.update(hard_stop=True, overall="PARTIAL")
        except Exception as error:
            row.update(status="FAIL", failure=str(error))
            summary.update(hard_stop=True, overall="PARTIAL")
            raise
        finally:
            summary["numerical_seconds"] += time.perf_counter() - started
            stage.save(bundle, item, summary)
        del history, buffer, records
        path.unlink()
        if summary.get("hard_stop"):
            return summary
        _complete_read_only(bundle, reference, times, authorization)
        stage.save(bundle, item, summary)
    compare_config = copy.deepcopy(science)
    compare_config["gates"] = previous.read_json(previous.base.one.runner.CONFIG)["gates"]
    with np.load(bundle / "one_d_p48_nonlinear.npz") as a, np.load(bundle / "one_d_p64_nonlinear.npz") as b:
        spatial = previous.base.one.runner.compare_histories(refs[48]["disc"], a,
            refs[64]["disc"], b, compare_config)
    previous.write_json(bundle / "one_d_all8_spatial.json", spatial)
    summary["p48_spatial_control"] = {"status": spatial["status"], "artifact": "one_d_all8_spatial.json"}
    summary["one_d_completed"] = True
    stage.save(bundle, item, summary)
    return summary


def evaluate_saved(bundle, item, times, p=64):
    """Restore dense nonlinear states and all linear factors without solves."""
    bundle, times = Path(bundle), np.asarray(times, dtype=float)
    reference, _, _ = _reference(item, p)
    disc = reference["disc"]
    if times.ndim != 1 or not len(times) or times[0] != 0. or np.any(np.diff(times) <= 0):
        raise ValueError("Comparison times must increase from zero")
    qv = previous.evaluate_dense_records(bundle / f"dense_p{p}.npz", times)
    x = np.linspace(0., 1., 41)
    result = {"times": times, "x": x}
    for kind in ("linear", "nonlinear"):
        fields, velocities = [], []
        for start in range(0, len(times), 256):
            stop = min(start + 256, len(times))
            if kind == "linear":
                state = disc.linear_reference(reference["saved"]["q_linear"], np.zeros(disc.ndof), times[start:stop])
                q, v = state["q"], state["velocity"]
                if start == 0:
                    q[0], v[0] = reference["saved"]["q_linear"], 0.
            else:
                q, v = qv[start:stop, :disc.ndof], qv[start:stop, disc.ndof:]
            fields.append(disc.reconstruct_series(q, x)); velocities.append(disc.reconstruct_series(v, x))
        result[kind + "_fields"] = np.concatenate(fields)
        result[kind + "_velocities"] = np.concatenate(velocities)
    result.update(initial_linear_fields=disc.reconstruct(reference["saved"]["q_linear"], x),
        initial_nonlinear_fields=disc.reconstruct(reference["saved"]["q_nonlinear"], x))
    return result
