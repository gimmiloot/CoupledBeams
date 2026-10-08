"""Analytical-time leading axial diagnostic; never advances an ODE.

The prescribed continuous first bending mode is common to every Shen space.
All discrete coordinates are retained. Old nonlinear histories are immutable
comparison data, not an exact reference for this asymptotic specialization.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import math
import os
import shutil
from pathlib import Path
import subprocess
import sys
import time

if __name__ == "__main__":
    for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
        os.environ[name] = "1"
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
CONFIG = ROOT/"data/input/planar_second_order_axial_response.json"
OUTPUT = ROOT/"results/planar_second_order_axial_response"
VERSION = "nlsp-second-order-axial-workflow-v1"
COMPONENTS = ("u2", "c2", "u2_t", "c2_t")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False, ensure_ascii=False,
                        default=lambda item: item.item() if isinstance(item, np.generic) else item.tolist())+"\n", encoding="utf8")
    temporary.replace(path)


def save_npz(path, **arrays):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    temporary.replace(path)


def identity(config_path=CONFIG, parent=None, checkpoint=None):
    config = read_json(config_path)
    files = ("scripts/analysis/verify_planar_second_order_axial_response.py",
             "scripts/lib/planar_second_order_axial_response.py",
             "scripts/lib/weakly_nonlinear_spatial_rod.py",
             "scripts/lib/weakly_nonlinear_planar_dynamics.py",
             "scripts/lib/mindlin_herrmann_longitudinal.py",
             "scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
             "scripts/analysis/diagnose_weakly_nonlinear_planar_rod.py", config["pilot_config"])
    value = {"version": VERSION, "config": config, "config_sha256": sha(config_path),
             "code_hashes": {name: sha(ROOT/name) for name in files},
             "historical_manifests": {name: sha(ROOT/path/"manifest.json") if (ROOT/path/"manifest.json").exists() else None
                                       for name, path in config["historical_bundles"].items()},
             "dependencies": {name: importlib.metadata.version(name) for name in ("numpy", "scipy", "matplotlib")}}
    pilot = read_json(ROOT/config["pilot_config"])
    value["reference_source_hashes"] = {name:{"manifest":sha(ROOT/pilot[name]/"manifest.json"),
                                             "result":sha(ROOT/pilot[name]/"result.json")}
                                      for name in ("audit_bundle", "linear_reference_bundle")}
    if parent is not None:
        value["conditional_p96_parent"] = {"bundle": str(Path(parent).resolve()), "manifest_sha256": sha(Path(parent)/"manifest.json")}
    if checkpoint is not None:
        value["p96_spectral_checkpoint"]={"path":str(Path(checkpoint).resolve()),"sha256":sha(checkpoint)}
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:16], value


def validate_cache(bundle, expected=None):
    bundle = Path(bundle)
    manifest = read_json(bundle/"manifest.json")
    if expected is not None and manifest["identity"] != expected:
        raise ValueError("Diagnostic cache identity mismatch")
    for name, digest in manifest["artifact_hashes"].items():
        if sha(bundle/name) != digest:
            raise ValueError(f"Diagnostic artifact hash mismatch: {name}")
    return read_json(bundle/"summary.json")


def historical_provenance(config):
    sources = {}
    for name, relative in config["historical_bundles"].items():
        folder = ROOT/relative
        if not (folder/"manifest.json").exists():
            sources[name] = {"status": "DATA_UNAVAILABLE", "bundle": relative}
            continue
        manifest = read_json(folder/"manifest.json")
        checked = {}
        for filename, digest in manifest["artifact_hashes"].items():
            path = folder/filename
            if not path.exists():
                sources[name] = {"status": "DATA_UNAVAILABLE", "bundle": relative, "missing": str(path)}
                break
            checked[filename] = sha(path) == digest
            if not checked[filename]:
                raise ValueError(f"Historical hash mismatch: {path}")
        else:
            cases = {}
            parent = folder/("cases" if name == "pilot" else "controls")
            for case in sorted(parent.iterdir()):
                if not (case/"trajectory.npz").exists():
                    continue
                metadata = read_json(case/"case.json")
                with np.load(case/"trajectory.npz") as saved:
                    times = saved["time"]
                    if np.any(np.diff(times) <= 0) or len(times) != metadata["samples"] or times[-1] != metadata["time_end"]:
                        raise ValueError(f"Historical timestamps mismatch: {case}")
                    row = {"samples": len(times), "actual_time_end": float(times[-1]), "p": metadata["p"],
                           "status": metadata["status"], "amplitude_over_h": metadata["amplitude_over_h"],
                           "time_level": metadata["time_level"]}
                    if "snapshot_indices" in saved:
                        row["actual_snapshot_times"] = times[saved["snapshot_indices"]].tolist()
                    cases[case.name] = row
            sources[name] = {"status": "PASS", "bundle": relative, "manifest_sha256": sha(folder/"manifest.json"),
                             "artifact_hashes_verified": checked, "cases": cases,
                             "historical_code_hash_not_compared_to_current": True}
    pilot = read_json(ROOT/config["pilot_config"])
    for name in ("audit_bundle", "linear_reference_bundle"):
        folder = ROOT/pilot[name]
        manifest = read_json(folder/"manifest.json")
        hashes = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
        assert hashes and all((folder/file).is_file() and sha(folder/file) == digest for file, digest in hashes.items())
        sources[name] = {"status": "PASS", "bundle": pilot[name], "manifest_sha256": sha(folder/"manifest.json")}
    return sources


class Budget:
    def __init__(self, config):
        self.started = time.perf_counter()
        self.limit = config["budget"]["numerical_wall_seconds"]
        self.prior = config["budget"]["prior_numerical_work_charge_seconds"]

    def check(self):
        if self.used() >= self.limit:
            raise TimeoutError("PARTIAL_COVERAGE: fixed actual numerical wall budget exhausted")

    def used(self):
        return self.prior+time.perf_counter()-self.started


def legendre_coefficients(model, coordinates, degree=None):
    """Physical polynomial coefficients, not raw whitening coefficients."""
    degree = model.p if degree is None else degree
    if degree < model.p:
        raise ValueError("Cannot truncate physical polynomial coordinates")
    result = np.zeros((len(coordinates), 2, degree+1))
    for index, field in enumerate(("u", "c")):
        raw = coordinates[:, model.slices[field]]@model.transforms[index].T
        result[:, index, :model.n] += raw
        result[:, index, 2:model.n+2] -= raw
    return result


def squared_l2(coefficients, length):
    weights = length/(2*np.arange(coefficients.shape[-1])+1)
    return np.einsum("tfn,n,tfn->tf", coefficients, weights, coefficients)


class Metric:
    def __init__(self, length, degree, points, config):
        self.length, self.points, self.config = length, points, config
        self.vander = legvander(2*points/length-1, degree)
        self.max_l2 = np.zeros(4)
        self.max_abs = np.zeros(4)
        self.ref_l2 = np.zeros(4)
        self.ref_max = np.zeros(4)
        self.l2_time = np.zeros(4)
        self.peak_time = np.zeros(4)
        self.peak_x = np.zeros(4)
        self.integrals = np.zeros(4)
        self.previous = {}

    def update(self, difference, reference, times, derivative):
        sl = slice(2*derivative, 2*derivative+2)
        norm = np.sqrt(squared_l2(difference, self.length))
        refnorm = np.sqrt(squared_l2(reference, self.length))
        values = (difference.reshape(-1, difference.shape[-1])@self.vander.T).reshape(len(times), 2, len(self.points))
        refvalues = (reference.reshape(-1, reference.shape[-1])@self.vander.T).reshape(len(times), 2, len(self.points))
        for field in range(2):
            k = 2*derivative+field
            row = int(np.argmax(norm[:, field]))
            if norm[row, field] > self.max_l2[k]:
                self.max_l2[k], self.l2_time[k] = norm[row, field], times[row]
            row, point = np.unravel_index(np.argmax(abs(values[:, field])), values[:, field].shape)
            if abs(values[row, field, point]) > self.max_abs[k]:
                self.max_abs[k], self.peak_time[k], self.peak_x[k] = abs(values[row, field, point]), times[row], self.points[point]
        self.ref_l2[sl] = np.maximum(self.ref_l2[sl], refnorm.max(axis=0))
        self.ref_max[sl] = np.maximum(self.ref_max[sl], abs(refvalues).max(axis=(0, 2)))
        sq = norm**2
        self.integrals[sl] += np.trapezoid(sq, times, axis=0)
        if derivative in self.previous:
            previous_time, previous_sq = self.previous[derivative]
            self.integrals[sl] += (times[0]-previous_time)*(previous_sq+sq[0])/2
        self.previous[derivative] = (times[-1], sq[-1])
        return norm

    def result(self):
        rows = {}
        for k, name in enumerate(COMPONENTS):
            floor = self.config["gates"]["relative_numerical_floor"]*max(self.ref_l2[k//2*2:k//2*2+2].max(), 1e-30)
            l2 = float(self.max_l2[k]/max(self.ref_l2[k], floor))
            maximum = float(self.max_abs[k]/max(self.ref_max[k], floor))
            rows[name] = {"absolute_L2": float(self.max_l2[k]), "absolute_max": float(self.max_abs[k]),
                          "reference_L2": float(self.ref_l2[k]), "reference_max": float(self.ref_max[k]),
                          "relative_L2": l2, "relative_max": maximum, "floor": float(floor),
                          "tolerance": self.config["gates"]["u_c_relative"],
                          "pass": l2 <= self.config["gates"]["u_c_relative"] and maximum <= self.config["gates"]["u_c_relative"],
                          "L2_peak_time": float(self.l2_time[k]), "max_peak_time": float(self.peak_time[k]),
                          "max_peak_x": float(self.peak_x[k]), "time_integral_squared_L2": float(self.integrals[k])}
        return {"fields": rows, "status": "PASS" if all(row["pass"] for row in rows.values()) else "PARTIAL"}


def pair_grid(low, high, config, factor=1):
    end = config["periods"]*high.background.T1
    omega = max(low.omega[-1], high.omega[-1], high.driving_omega)
    count = int(math.ceil(end*omega/(2*np.pi)*config["sampling"]["samples_per_fastest_retained_period"]*factor))
    return np.unique(np.r_[np.linspace(0, end, count+1), config["short_periods"]*high.background.T1])


def convergence_pair(low, high, config, bundle, budget, factor=1, projection=False):
    from scripts.analysis.diagnose_weakly_nonlinear_planar_rod import projection_operator
    times = pair_grid(low, high, config, factor)
    nodes, _ = leggauss(config["sampling"]["comparison_Gauss_points"])
    points = (nodes+1)*high.length/2
    metric = Metric(high.length, high.p, points, config)
    short = Metric(high.length, high.p, points, config)
    projectors = projection_operator(low.p, high.p, high.coefficients, length=high.length) if projection else None
    projection_integrals = {part: np.zeros(4) for part in ("total", "tail", "common")}
    projection_residual = np.zeros(4)
    projection_previous = {}
    rows, norm_samples = [], []
    stride = max(1, int(math.ceil(len(times)/5000)))
    started = time.perf_counter()
    block = config["sampling"]["block_size"]
    for start in range(0, len(times), block):
        budget.check()
        t = times[start:start+block]
        local_norms = []
        for derivative in (0, 1):
            xlow, xhigh = low.evaluate(t, derivative), high.evaluate(t, derivative)
            a, b = legendre_coefficients(low, xlow, high.p), legendre_coefficients(high, xhigh)
            difference = b-a
            norm = metric.update(difference, b, t, derivative)
            mask = t <= config["short_periods"]*high.background.T1
            if mask.any():
                short.update(difference[mask], b[mask], t[mask], derivative)
            local_norms.append(norm)
            if projection:
                projected = np.column_stack([xhigh[:, high.slices[name]]@projectors[index].T
                                             for name, index in (("u", 0), ("c", 3))])
                plow = legendre_coefficients(low, projected, high.p)
                tail, common = b-plow, plow-a
                weights = high.length/(2*np.arange(high.p+1)+1)
                cross = np.einsum("tfn,n,tfn->tf", tail, weights, common)
                sq = {"total": squared_l2(difference, high.length),
                      "tail": squared_l2(tail, high.length), "common": squared_l2(common, high.length)}
                psl = slice(2*derivative, 2*derivative+2)
                projection_residual[psl] = np.maximum(projection_residual[psl], np.maximum(np.max(abs(cross), axis=0),
                                                np.max(abs(sq["total"]-sq["tail"]-sq["common"]), axis=0)))
                for name in projection_integrals:
                    key = (derivative, name)
                    values = sq[name]
                    projection_integrals[name][psl] += np.trapezoid(values, t, axis=0)
                    if key in projection_previous:
                        pt, pv = projection_previous[key]
                        projection_integrals[name][psl] += (t[0]-pt)*(pv+values[0])/2
                    projection_previous[key] = (t[-1], values[-1])
        chosen = np.flatnonzero((np.arange(start, start+len(t)) % stride) == 0)
        rows.extend(t[chosen].tolist())
        norm_samples.extend(np.column_stack(local_norms)[chosen].tolist())
    result = metric.result()
    result.update({"pair": [low.p, high.p], "samples": len(times), "time_end": float(times[-1]),
                   "sampling_factor": factor, "max_dt": float(np.diff(times).max()),
                   "samples_per_fastest_retained_period": config["sampling"]["samples_per_fastest_retained_period"]*factor,
                   "short_0p1T1": short.result(), "runtime_seconds": time.perf_counter()-started,
                   "maximum_qualification": config["sampling"]["maximum_qualification"]})
    if projection:
        result["projection"] = {"absolute_identity_residual": projection_residual.tolist(),
                                 "identity_scaled_by_own_reference_L2_squared": (projection_residual/np.maximum(metric.ref_l2**2, 1e-30)).tolist(),
                                 "time_integrals_fields_and_velocities": {key: value.tolist() for key, value in projection_integrals.items()},
                                 "tail_fractions_fields_and_velocities": (projection_integrals["tail"]/np.maximum(projection_integrals["total"], 1e-30)).tolist(),
                                 "qualification": "physical L2 split; common-space part is not proven exclusively phase error"}
    tag = f"p{low.p}_p{high.p}_sampling{factor}"
    write_json(bundle/"convergence"/(tag+".json"), result)
    save_npz(bundle/"convergence"/(tag+".npz"), time=np.asarray(rows), difference_L2=np.asarray(norm_samples))
    selected = np.unique(np.r_[np.asarray(config["sampling"]["snapshots_periods"])*high.background.T1,
                               metric.l2_time, metric.peak_time])
    dense = np.linspace(0, high.length, 1001)
    payload = {"time": selected, "s": dense}
    for order, name in ((0, "q"), (1, "velocity")):
        payload[name+"_high_minus_low"] = high.reconstruct_series(high.evaluate(selected, order), dense)-low.reconstruct_series(low.evaluate(selected, order), dense)
    save_npz(bundle/"convergence"/(tag+"_profiles.npz"), **payload)
    print(json.dumps({"pair": [low.p, high.p], "sampling_factor": factor, "status": result["status"],
                      "seconds": result["runtime_seconds"], "budget_used": budget.used()}), flush=True)
    return result



def nonlinear_coordinates(saved, name, model):
    value = saved[name]
    if value.shape[1] != 4*model.n:
        raise ValueError("Historical coordinate dimension mismatch")
    return np.column_stack((value[:, :model.n], value[:, 3*model.n:4*model.n]))


def compare_saved_case(model, case, config, bundle, budget, label):
    metadata = read_json(case/"case.json")
    epsilon = metadata["amplitude_over_h"]
    nodes, _ = leggauss(config["sampling"]["comparison_Gauss_points"])
    points = (nodes+1)*model.length/2
    metric = Metric(model.length, model.p, points, config)
    traces = {"time": [], "nonlinear_normalized": [], "leading": []}
    observation_points = np.array(config["sampling"]["observation_points"])*model.length
    with np.load(case/"trajectory.npz") as saved:
        times = saved["time"]
        if times[-1] != metadata["time_end"] or len(times) != metadata["samples"]:
            raise ValueError("Cannot substitute nominal time for a historical partial case")
        q, velocity = nonlinear_coordinates(saved, "q", model), nonlinear_coordinates(saved, "velocity", model)
        stride = max(1, len(times)//2000)
        for start in range(0, len(times), config["sampling"]["block_size"]):
            budget.check()
            t = times[start:start+config["sampling"]["block_size"]]
            collected = []
            for derivative, history in ((0, q), (1, velocity)):
                old = history[start:start+len(t)]/epsilon**2
                leading = model.evaluate(t, derivative)
                a, b = legendre_coefficients(model, old), legendre_coefficients(model, leading)
                metric.update(a-b, b, t, derivative)
                chosen = np.flatnonzero(((np.arange(start, start+len(t)) % stride) == 0) | (t <= config["short_periods"]*model.background.T1))
                if derivative == 0:
                    traces["time"].extend(t[chosen].tolist())
                collected.append((model.reconstruct_series(old[chosen], observation_points),
                                  model.reconstruct_series(leading[chosen], observation_points)))
            traces["nonlinear_normalized"].extend(np.concatenate([entry[0] for entry in collected], axis=2).tolist())
            traces["leading"].extend(np.concatenate([entry[1] for entry in collected], axis=2).tolist())
    result = metric.result()
    result.update({"p": model.p, "epsilon_a": epsilon, "time_end": float(times[-1]),
                   "periods_covered": float(times[-1]/model.background.T1), "samples": len(times),
                   "historical_case": str(case.relative_to(ROOT)), "case_sha256": sha(case/"case.json"),
                   "trajectory_sha256": sha(case/"trajectory.npz"), "exact_saved_times": True,
                   "qualification": "Same p but common continuous prescribed bending vs historical projected semidiscrete initial mode; discrepancy includes higher amplitude orders, not solely numerical error",
                   "coverage": "FULL_5T1" if times[-1] == config["periods"]*model.background.T1 else "ACTUAL_SAVED_PREFIX"})
    write_json(bundle/"nonlinear"/(label+".json"), result)
    save_npz(bundle/"nonlinear"/(label+"_observations.npz"),
             **{key: np.asarray(value) for key, value in traces.items()}, s=observation_points)
    return result


def compare_old_discretization(models, config, bundle, budget):
    pilot = ROOT/config["historical_bundles"]["pilot"]
    low, high = models[24], models[32]
    low_case, high_case = pilot/"cases/p24_Aoverh0p05_tight", pilot/"cases/p32_Aoverh0p05_tight"
    nodes, _ = leggauss(config["sampling"]["comparison_Gauss_points"])
    points = (nodes+1)*high.length/2
    metrics = {name: Metric(high.length, high.p, points, config) for name in ("nonlinear_difference", "leading_difference", "remainder")}
    integrals = {name: np.zeros(4) for name in ("old", "leading", "dot")}
    previous = {}
    with np.load(low_case/"trajectory.npz") as a, np.load(high_case/"trajectory.npz") as b:
        times = a["time"]
        if not np.array_equal(times, b["time"]):
            raise ValueError("Discretization comparison cannot interpolate historical times")
        histories = [(nonlinear_coordinates(a, name, low), nonlinear_coordinates(b, name, high)) for name in ("q", "velocity")]
        for start in range(0, len(times), config["sampling"]["block_size"]):
            budget.check()
            t = times[start:start+config["sampling"]["block_size"]]
            for derivative, (old_low, old_high) in enumerate(histories):
                old = legendre_coefficients(high, old_high[start:start+len(t)])-legendre_coefficients(low, old_low[start:start+len(t)], high.p)
                leading = .05**2*(legendre_coefficients(high, high.evaluate(t, derivative))-legendre_coefficients(low, low.evaluate(t, derivative), high.p))
                metrics["nonlinear_difference"].update(old, old, t, derivative)
                metrics["leading_difference"].update(leading, old, t, derivative)
                metrics["remainder"].update(old-leading, old, t, derivative)
                weights = high.length/(2*np.arange(high.p+1)+1)
                samples = {"old": squared_l2(old, high.length), "leading": squared_l2(leading, high.length),
                           "dot": np.einsum("tfn,n,tfn->tf", old, weights, leading)}
                sl = slice(2*derivative, 2*derivative+2)
                for name, values in samples.items():
                    integrals[name][sl] += np.trapezoid(values, t, axis=0)
                    if (name, derivative) in previous:
                        pt, pv = previous[name, derivative]
                        integrals[name][sl] += (t[0]-pt)*(pv+values[0])/2
                    previous[name, derivative] = t[-1], values[-1]
    correlation = integrals["dot"]/np.sqrt(np.maximum(integrals["old"]*integrals["leading"], 1e-60))
    result = {name: metric.result() for name, metric in metrics.items()}
    result.update({"p": [24, 32], "epsilon_a": .05, "samples": len(times), "time_end": float(times[-1]),
                   "space_time_L2_correlations": dict(zip(COMPONENTS, correlation.tolist())),
                   "qualification": "No phase alignment; correlation is diagnostic similarity, not proof that all nonlinear error is explained"})
    write_json(bundle/"nonlinear/discretization_difference_p24_p32.json", result)
    selected = np.array([.001, .01, .1, 1., 4., 5.])*high.background.T1
    indices = np.array([int(np.argmin(abs(times-t))) for t in selected])
    with np.load(low_case/"trajectory.npz") as a, np.load(high_case/"trajectory.npz") as b:
        dense = np.linspace(0, high.length, 1001)
        payload = {"time": times[indices], "s": dense}
        for derivative, name in enumerate(("q", "velocity")):
            ol = nonlinear_coordinates(a, name, low)[indices]
            oh = nonlinear_coordinates(b, name, high)[indices]
            payload[name+"_nonlinear_difference"] = high.reconstruct_series(oh, dense)-low.reconstruct_series(ol, dense)
            payload[name+"_leading_difference"] = .05**2*(high.reconstruct_series(high.evaluate(times[indices], derivative), dense)-low.reconstruct_series(low.evaluate(times[indices], derivative), dense))
        save_npz(bundle/"nonlinear/discretization_difference_profiles.npz", **payload)
    return result


def compare_amplitudes(model, config, bundle, budget):
    parent = ROOT/config["historical_bundles"]["pilot"]/"cases"
    nodes, _ = leggauss(config["sampling"]["comparison_Gauss_points"])
    points = (nodes+1)*model.length/2
    metric = Metric(model.length, model.p, points, config)
    with np.load(parent/"p32_Aoverh0p05_tight/trajectory.npz") as large, np.load(parent/"p32_Aoverh0p025_tight/trajectory.npz") as small:
        times = large["time"]
        if not np.array_equal(times, small["time"]):
            raise ValueError("Amplitude comparison cannot interpolate old times")
        histories = [(nonlinear_coordinates(large, name, model)/.05**2,
                      nonlinear_coordinates(small, name, model)/.025**2) for name in ("q", "velocity")]
        for start in range(0, len(times), config["sampling"]["block_size"]):
            budget.check()
            t = times[start:start+config["sampling"]["block_size"]]
            for derivative, (a, b) in enumerate(histories):
                normalized_difference = legendre_coefficients(model, a[start:start+len(t)]-b[start:start+len(t)])
                reference = legendre_coefficients(model, model.evaluate(t, derivative))
                metric.update(normalized_difference, reference, t, derivative)
    result = metric.result()
    result.update({"p": 32, "samples": len(times), "exact_saved_times": True,
                   "qualification": "Difference of two existing normalized cubic histories; not a fitted amplitude power law",
                   "asymptotic_quarter_rescaling": "(.025/.05)^2=1/4 exactly; arithmetic identity only"})
    write_json(bundle/"nonlinear/amplitude_normalization.json", result)
    return result


def augmented_exponential_check(model, times):
    from scipy.linalg import expm
    n = model.ndof
    A = np.zeros((2*n+3, 2*n+3))
    A[:n, n:2*n] = np.eye(n)
    A[n:2*n, :n] = -np.linalg.solve(model.M, model.K)
    A[n:2*n, 2*n] = np.linalg.solve(model.M, model.f0)
    A[n:2*n, 2*n+1] = np.linalg.solve(model.M, model.f2)
    A[2*n+1, 2*n+2], A[2*n+2, 2*n+1] = -model.driving_omega, model.driving_omega
    initial = np.zeros(2*n+3)
    initial[2*n:2*n+2] = 1.
    reference = np.array([expm(A*t)@initial for t in times])
    exact = np.column_stack((model.evaluate(times), model.evaluate(times, 1)))
    relative = np.linalg.norm(exact-reference[:, :2*n], axis=1)/np.maximum(np.linalg.norm(exact, axis=1), 1e-30)
    return {"times": np.asarray(times).tolist(), "relative": relative.tolist(), "maximum_relative": float(relative.max()),
            "matrix_exponentials": len(times), "time_integrations": 0}


def run_compute(config, bundle):
    from scripts.lib import planar_second_order_axial_response as axial
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    budget = Budget(config)
    pilot = read_json(ROOT/config["pilot_config"])
    provenance = historical_provenance(config)
    write_json(bundle/"source_provenance.json", provenance)
    accepted = read_json(ROOT/pilot["audit_bundle"]/"result.json")
    coefficients = rod.RodCoefficients(**accepted["coefficients"])
    background = axial.background_from_pilot(pilot, ROOT)
    model = rod.derive_polynomials()
    derivation = axial.derive_second_order(model)
    write_json(bundle/"derivation.json", {key: value if key in ("checks", "model_version") else str(value) for key, value in derivation.items()})
    statuses = {"NLSP_SECOND_ORDER_DERIVATION": "PASS", "NLSP_SECOND_ORDER_FORCING_ASSEMBLY": "NOT_RUN",
                "NLSP_SECOND_ORDER_EXACT_TIME_EVALUATOR": "NOT_RUN", "NLSP_SECOND_ORDER_SPATIAL_CONVERGENCE": "PARTIAL_COVERAGE",
                "NLSP_SECOND_ORDER_NONLINEAR_COMPARISON": "NOT_RUN", "NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC": "PARTIAL_COVERAGE"}
    result = {"statuses": statuses, "config": config, "background": background.as_dict(), "coefficients": coefficients.values(),
              "models": {}, "convergence": {}, "nonlinear_comparisons": {}, "time_integrations": 0,
              "previous_nonlinear_statuses": {"pilot": "PARTIAL", "recovery": "PARTIAL"}, "full_nonlinear_spatial_convergence": "UNRESOLVED",
              "all_coordinates_retained": True, "energy_classification": False, "phase_alignment": False}
    models = {}
    try:
        observations = np.linspace(0, config["periods"]*background.T1, config["sampling"]["observation_samples"])
        for p in config["degrees"]:
            budget.check()
            disc = axial.SecondOrderAxial(coefficients, p, background, model=model)
            models[p] = disc
            quadrature = disc.forcing_checks((2*p+1, 3*p+3, 4*p+5))
            identities = disc.equation_and_power_checks(background.T1*np.array([0., 1e-8, .001, .1, .25, 1., 5.]))
            assert max(disc.matrix_checks[key] for key in ("mass_relative_difference", "stiffness_relative_difference")) <= config["gates"]["matrix_match_relative"]
            assert disc.eigen_checks["eigenpair_relative_residual"] <= config["gates"]["eigenpair_scaled"]
            assert disc.eigen_checks["mass_orthogonality_max_absolute"] <= config["gates"]["eigenpair_scaled"]
            assert identities["equation_scaled_max"] <= config["gates"]["exact_time_equation_scaled"]
            assert identities["power_identity_scaled_max"] <= config["gates"]["power_scaled"]
            assert all(max(row["strong_weak_relative"].values()) <= config["quadrature"]["strong_weak_relative_tolerance"] for row in quadrature["rows"])
            assert all(max(row["relative_changes"].values()) <= config["quadrature"]["forcing_relative_tolerance"] for row in quadrature["rows"][1:])
            row = {"p": p, "coordinates": disc.ndof, "matrix_nq": disc.matrix_nq, "forcing_nq": disc.nq,
                   "matrix_checks": disc.matrix_checks, "eigen_checks": disc.eigen_checks,
                   "forcing_checks": quadrature, "exact_time_checks": identities}
            result["models"][str(p)] = row
            save_npz(bundle/"models"/f"p{p}.npz", M=disc.M, K=disc.K, f0=disc.f0, f2=disc.f2,
                     omega=disc.omega, vectors=disc.vectors, b0=disc.b0, b2=disc.b2,
                     detuning=disc.omega-disc.driving_omega, transform_u=disc.transforms[0], transform_c=disc.transforms[1])
            write_json(bundle/"models"/f"p{p}.json", row)
            observation_points = np.asarray(config["sampling"]["observation_points"])*background.length
            save_npz(bundle/"models"/f"p{p}_observations.npz", time=observations, s=observation_points,
                     q=disc.reconstruct_series(disc.evaluate(observations), observation_points),
                     velocity=disc.reconstruct_series(disc.evaluate(observations, 1), observation_points))
            short_end = config["short_periods"]*background.T1
            short_count = int(math.ceil(short_end*disc.omega[-1]/(2*np.pi)*config["sampling"]["samples_per_fastest_retained_period"]))
            short_observations = np.linspace(0, short_end, short_count+1)
            save_npz(bundle/"models"/f"p{p}_short_observations.npz", time=short_observations, s=observation_points,
                     q=disc.reconstruct_series(disc.evaluate(short_observations), observation_points),
                     velocity=disc.reconstruct_series(disc.evaluate(short_observations, 1), observation_points))
        statuses["NLSP_SECOND_ORDER_FORCING_ASSEMBLY"] = "PASS"
        check = augmented_exponential_check(models[16], background.T1*np.array([0., 1e-5, .01, .1, 1., 5.]))
        write_json(bundle/"augmented_exponential.json", check)
        if check["maximum_relative"] > config["gates"]["expm_relative"]:
            raise ArithmeticError("Independent matrix exponential check failed")
        statuses["NLSP_SECOND_ORDER_EXACT_TIME_EVALUATOR"] = "PASS"
        for lo, hi in zip(config["degrees"][:-1], config["degrees"][1:]):
            result["convergence"][f"p{lo}_p{hi}"] = convergence_pair(models[lo], models[hi], config, bundle, budget, projection=(lo, hi) in ((24, 32), (48, 64)))
            write_json(bundle/"progress.json", result)
        last_lo, last_hi = config["degrees"][-2:]
        dense = convergence_pair(models[last_lo], models[last_hi], config, bundle, budget,
                                 factor=config["sampling"]["final_pair_refinement_factor"], projection=True)
        coarse = result["convergence"][f"p{last_lo}_p{last_hi}"]
        changes = {name: {key: abs(dense["fields"][name][key]/max(coarse["fields"][name][key], 1e-30)-1)
                          for key in ("absolute_L2", "absolute_max", "reference_L2", "reference_max")} for name in COMPONENTS}
        result["sampling_refinement"] = {"pair": [last_lo, last_hi], "coarse_samples": coarse["samples"],
                 "fine_samples": dense["samples"], "relative_changes": changes,
                 "status": "PASS" if max(value for fields in changes.values() for value in fields.values()) <= config["sampling"]["maxima_refinement_relative_target"] else "PARTIAL",
                 "qualification": "sampling convergence diagnostic only; not a continuous supremum certificate"}
        result["final_refined_pair"] = dense
        statuses["NLSP_SECOND_ORDER_SPATIAL_CONVERGENCE"] = "PASS" if dense["status"] == "PASS" and result["sampling_refinement"]["status"] == "PASS" else "PARTIAL"
        result["optional_p96"] = {"performed": False, "qualification": "No automatic escalation; primary p16...64 evidence reported before choosing any further diagnostic level"}
        if provenance["pilot"]["status"] == "PASS":
            parent = ROOT/config["historical_bundles"]["pilot"]/"cases"
            for p, amplitude, label in ((24, "0p05", "p24_large"), (32, "0p05", "p32_large"),
                                        (24, "0p025", "p24_small_prefix"), (32, "0p025", "p32_small")):
                case = parent/f"p{p}_Aoverh{amplitude}_tight"
                if (case/"trajectory.npz").exists():
                    result["nonlinear_comparisons"][label] = compare_saved_case(models[p], case, config, bundle, budget, label)
            result["discretization_difference_comparison"] = compare_old_discretization(models, config, bundle, budget)
            result["amplitude_normalization"] = compare_amplitudes(models[32], config, bundle, budget)
            preflight = read_json(ROOT/config["historical_bundles"]["pilot"]/"preflight.json")
            result["historical_background_qualification"] = {"source": "validated historical preflight, no new background roots",
                "rows": [{key: row[key] for key in ("p", "omega", "linear_relative_errors", "projection_L2")} for row in preflight["spatial_linear_controls"] if row["p"] in (24, 32)],
                "common_continuous_background_not_exact_historical_semidiscrete_asymptotic_coefficient": True}
            statuses["NLSP_SECOND_ORDER_NONLINEAR_COMPARISON"] = "COMPLETE_WITH_ACTUAL_PREFIX_QUALIFICATIONS"
        else:
            statuses["NLSP_SECOND_ORDER_NONLINEAR_COMPARISON"] = "DATA_UNAVAILABLE"
        if provenance["recovery"]["status"] == "PASS":
            case = ROOT/config["historical_bundles"]["recovery"]/"controls/p48_strict_short"
            result["nonlinear_comparisons"]["p48_short"] = compare_saved_case(models[48], case, config, bundle, budget, "p48_short")
        statuses["NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC"] = "COMPLETE" if statuses["NLSP_SECOND_ORDER_NONLINEAR_COMPARISON"] != "DATA_UNAVAILABLE" else "PARTIAL_COVERAGE"
    except TimeoutError as exception:
        result["stop_reason"] = str(exception)
    result["runtime"] = {"charged_numerical_seconds": budget.used(), "limit_seconds": budget.limit,
                         "primary_MH_eigendecompositions": sum(disc.eigen_decompositions for disc in models.values()),
                         "exact_time_evaluations": sum(disc.exact_time_evaluations for disc in models.values()),
                         "new_time_integrations": 0, "continuum_root_solves": 0,
                         "model_derivations": 1, "augmented_matrix_exponentials": 6 if "16" in result["models"] else 0,
                         "threads": {name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS")}}
    return result




def run_extend_p96(config, bundle, parent, checkpoint=None):
    """One logged conditional level; restore p64, never repeat p16...64 solves."""
    from scripts.lib import planar_second_order_axial_response as axial
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    parent=Path(parent).resolve()
    result = validate_cache(parent)
    recovered = read_json(checkpoint) if checkpoint is not None else None
    if recovered is not None and (recovered["new_time_integrations"]!=0 or sha(recovered["p96_file"])!=recovered["p96_sha256"]):
        raise ValueError("Invalid p96 spectral checkpoint")
    if result["statuses"]["NLSP_SECOND_ORDER_SPATIAL_CONVERGENCE"] == "PASS":
        raise ValueError("No unresolved leading-response gate justifies optional p96")
    if "96" in result["models"] or config["optional_degree"] != 96:
        raise ValueError("Only one optional p96 level is authorized")
    if result["config"] != config:
        raise ValueError("Conditional refinement must preserve the complete primary config")
    budget = Budget(config)
    budget.prior = result["runtime"]["charged_numerical_seconds"]+(recovered["charged_failed_attempt_seconds"] if recovered else 0.)
    budget.check()
    shutil.copytree(parent, bundle, dirs_exist_ok=True)
    pilot = read_json(ROOT/config["pilot_config"])
    coefficients = rod.RodCoefficients(**result["coefficients"])
    background = axial.background_from_pilot(pilot, ROOT)
    model = rod.derive_polynomials()
    low = axial.SecondOrderAxial.from_saved(coefficients, 64, background, Path(parent)/"models/p64.npz", model)
    high = (axial.SecondOrderAxial.from_saved(coefficients, 96, background, recovered["p96_file"], model)
            if recovered else axial.SecondOrderAxial(coefficients, 96, background, model))
    quadrature = high.forcing_checks((193, 291, 389))
    checks = high.equation_and_power_checks(background.T1*np.array([0., 1e-8, .001, .1, .25, 1., 5.]))
    assert max(high.matrix_checks[key] for key in ("mass_relative_difference", "stiffness_relative_difference")) <= config["gates"]["matrix_match_relative"]
    assert high.eigen_checks["eigenpair_relative_residual"] <= config["gates"]["eigenpair_scaled"]
    assert high.eigen_checks["mass_orthogonality_max_absolute"] <= config["gates"]["eigenpair_scaled"]
    assert checks["equation_scaled_max"] <= config["gates"]["exact_time_equation_scaled"] and checks["power_identity_scaled_max"] <= config["gates"]["power_scaled"]
    assert all(max(row["relative_changes"].values()) <= config["quadrature"]["forcing_relative_tolerance"] for row in quadrature["rows"][1:])
    assert all(max(row["strong_weak_relative"].values()) <= config["quadrature"]["strong_weak_relative_tolerance"] for row in quadrature["rows"])
    row = {"p":96,"coordinates":high.ndof,"matrix_nq":high.matrix_nq,"forcing_nq":high.nq,
           "matrix_checks":high.matrix_checks,"eigen_checks":high.eigen_checks,"forcing_checks":quadrature,"exact_time_checks":checks}
    result["models"]["96"] = row
    save_npz(bundle/"models/p96.npz",M=high.M,K=high.K,f0=high.f0,f2=high.f2,omega=high.omega,vectors=high.vectors,
             b0=high.b0,b2=high.b2,detuning=high.omega-high.driving_omega,transform_u=high.transforms[0],transform_c=high.transforms[1])
    write_json(bundle/"models/p96.json",row)
    short_end=config["short_periods"]*background.T1
    short_count=int(math.ceil(short_end*high.omega[-1]/(2*np.pi)*config["sampling"]["samples_per_fastest_retained_period"]))
    times=np.linspace(0,short_end,short_count+1);points=np.array(config["sampling"]["observation_points"])*background.length
    save_npz(bundle/"models/p96_short_observations.npz",time=times,s=points,
             q=high.reconstruct_series(high.evaluate(times),points),velocity=high.reconstruct_series(high.evaluate(times,1),points))
    result["parent_sampling_refinement"] = result["sampling_refinement"]
    result["parent_final_refined_pair"] = result["final_refined_pair"]
    result["optional_p96"] = {"performed":True,"decision":"Conditional diagnostic: p32-to48 c difference grows again; p48-to64 passes u only, so adequacy of p64 for the leading transient is unresolved",
         "parent_bundle":str(Path(parent).relative_to(ROOT)),"parent_manifest_sha256":sha(Path(parent)/"manifest.json"),
         "old_spatial_cases_recomputed":False,"old_nonlinear_comparisons_recomputed":False,"p64_spectral_restore_eigendecompositions":low.eigen_decompositions,
         "new_MH_eigendecompositions":high.eigen_decompositions+(recovered["new_MH_eigendecompositions"] if recovered else 0),"restore_checks":low.restore_checks,
         "spectral_checkpoint_recovery":recovered, "new_eigh_on_recovery":high.eigen_decompositions}
    result["statuses"]["NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC"]="PARTIAL_COVERAGE"
    try:
        coarse=convergence_pair(low,high,config,bundle,budget,projection=True)
        result["convergence"]["p64_p96"]=coarse
        dense=convergence_pair(low,high,config,bundle,budget,factor=2,projection=True)
        changes={name:{key:abs(dense["fields"][name][key]/max(coarse["fields"][name][key],1e-30)-1)
                      for key in ("absolute_L2","absolute_max","reference_L2","reference_max")} for name in COMPONENTS}
        result["sampling_refinement"]={"pair":[64,96],"coarse_samples":coarse["samples"],"fine_samples":dense["samples"],"relative_changes":changes,
            "status":"PASS" if max(value for fields in changes.values() for value in fields.values())<=config["sampling"]["maxima_refinement_relative_target"] else "PARTIAL",
            "qualification":"sampled maximum diagnostic, not a continuous supremum certificate"}
        result["final_refined_pair"]=dense
        result["statuses"]["NLSP_SECOND_ORDER_SPATIAL_CONVERGENCE"]="PASS" if dense["status"]=="PASS" and result["sampling_refinement"]["status"]=="PASS" else "PARTIAL"
        result["statuses"]["NLSP_SECOND_ORDER_AXIAL_DIAGNOSTIC"]="COMPLETE"
    except TimeoutError as exception:
        result["stop_reason"]=str(exception)
    previous=result["runtime"]
    result["runtime"]={**previous,"charged_numerical_seconds":budget.used(),"primary_MH_eigendecompositions":previous["primary_MH_eigendecompositions"]+high.eigen_decompositions+(recovered["new_MH_eigendecompositions"] if recovered else 0),
       "exact_time_evaluations":previous["exact_time_evaluations"]+low.exact_time_evaluations+high.exact_time_evaluations,
       "model_derivations":previous["model_derivations"]+1+(1 if recovered else 0),"cached_MH_spectral_restorations":1+(1 if recovered else 0),
       "prior_primary_charged_seconds":previous["charged_numerical_seconds"],"new_time_integrations":0}
    return result


def plot_only(bundle):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bundle = Path(bundle)
    summary = validate_cache(bundle)
    T1 = summary["background"]["T1"]
    folder = bundle/"figures"
    folder.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": .2})
    fig, axes = plt.subplots(2, 2, figsize=(8, 5.2), layout="constrained", sharex=True)
    for p in (24, 32, 48, 64, 96):
        path = bundle/"models"/f"p{p}_short_observations.npz"
        if not path.exists():
            continue
        with np.load(path) as data:
            for index, axis in enumerate(axes.flat):
                values = data["q" if index < 2 else "velocity"][:, 0, index%2]
                axis.plot(data["time"]/T1, values, lw=.7, label=f"p={p}")
                axis.set_ylabel(COMPONENTS[index]); axis.set_xlabel(r"$t/T_1$")
    axes[0, 0].legend(fontsize=8)
    fig.savefig(folder/"leading_response.pdf"); fig.savefig(folder/"leading_response.png", dpi=220); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), layout="constrained")
    for name in COMPONENTS:
        entries = list(summary["convergence"].values())
        for axis, key in zip(axes, ("relative_L2", "relative_max")):
            axis.semilogy([row["pair"][1] for row in entries], [row["fields"][name][key] for row in entries], "o-", label=name, lw=.8, ms=3)
            axis.axhline(.001, color="grey", lw=.6, ls="--"); axis.set_xlabel("higher polynomial degree p"); axis.set_ylabel(key)
    axes[0].legend(fontsize=8)
    fig.savefig(folder/"spatial_convergence.pdf"); fig.savefig(folder/"spatial_convergence.png", dpi=220); plt.close(fig)
    path = bundle/"nonlinear/p32_large_observations.npz"
    if path.exists():
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.5), layout="constrained", sharex=True)
        for label, style in (("p32_large", "-"), ("p32_small", ":")):
            with np.load(bundle/"nonlinear"/(label+"_observations.npz")) as data:
                for index, axis in enumerate(axes):
                    mask = data["time"] <= summary["config"]["short_periods"]*T1
                    axis.plot(data["time"][mask]/T1, data["nonlinear_normalized"][mask, 0, index], style, lw=.7,
                              label="cubic / epsilon_a^2, "+label.replace("p32_", ""))
                if label == "p32_large":
                    for index, axis in enumerate(axes):
                        axis.plot(data["time"][mask]/T1, data["leading"][mask, 0, index], "--", lw=.7, label="leading second order")
        for index, axis in enumerate(axes):
            axis.set_ylabel(COMPONENTS[index]); axis.set_xlabel(r"$t/T_1$")
        axes[0].legend(fontsize=7)
        fig.savefig(folder/"nonlinear_comparison.pdf"); fig.savefig(folder/"nonlinear_comparison.png", dpi=220); plt.close(fig)
    return [str(path) for path in sorted(folder.glob("*.pdf"))]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--compute", action="store_true")
    actions.add_argument("--extend-p96", type=Path)
    actions.add_argument("--report-only", type=Path)
    actions.add_argument("--plot-only", type=Path)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--p96-spectral-checkpoint", type=Path)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if args.report_only or args.plot_only:
        bundle = args.report_only or args.plot_only
        summary = validate_cache(bundle)
        figures = plot_only(bundle) if args.plot_only else []
        print(json.dumps({"bundle": str(bundle), "statuses": summary["statuses"], "figures": figures,
                          "new_time_integrations": 0, "eigendecompositions": 0, "analytic_evaluations": 0, "derivations": 0}))
        return 0
    key, expected = identity(args.config, args.extend_p96, args.p96_spectral_checkpoint) if args.extend_p96 else identity(args.config)
    bundle = args.output_dir/key
    if (bundle/"manifest.json").exists():
        summary = validate_cache(bundle, expected)
        print(json.dumps({"bundle": str(bundle), "statuses": summary["statuses"], "cache": "reused",
                          "new_time_integrations": 0, "eigendecompositions": 0, "analytic_evaluations": 0, "derivations": 0}))
        return 0
    # Completed conditional coverage contains the primary grid; reuse it only
    # under the same current code/config/source identity, never by filename.
    current_path = args.output_dir/"current.json"
    if args.compute and current_path.is_file():
        candidate = Path(read_json(current_path)["bundle"])
        candidate_manifest = read_json(candidate/"manifest.json")
        base_identity = dict(candidate_manifest["identity"])
        conditional = base_identity.pop("conditional_p96_parent", None)
        base_identity.pop("p96_spectral_checkpoint", None)
        if conditional is not None and base_identity == expected:
            summary = validate_cache(candidate)
            if all(str(p) in summary["models"] for p in expected["config"]["degrees"]):
                print(json.dumps({"bundle":str(candidate),"statuses":summary["statuses"],"cache":"reused_conditional_coverage",
                                  "new_time_integrations":0,"eigendecompositions":0,"analytic_evaluations":0,"derivations":0}))
                return 0
    bundle.mkdir(parents=True, exist_ok=True)
    summary = (run_extend_p96(read_json(args.config), bundle, args.extend_p96, args.p96_spectral_checkpoint) if args.extend_p96
               else run_compute(read_json(args.config), bundle))
    write_json(bundle/"summary.json", summary)
    manifest = {"identity": expected, "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                "command": sys.argv, "artifact_hashes": {str(path.relative_to(bundle)): sha(path) for path in bundle.rglob("*") if path.is_file() and path.name not in ("manifest.json",) and not path.name.endswith(".tmp")}}
    write_json(bundle/"manifest.json", manifest)
    write_json(args.output_dir/"current.json", {"fingerprint": key, "bundle": str(bundle)})
    print(json.dumps({"bundle": str(bundle), "statuses": summary["statuses"], "runtime": summary["runtime"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
