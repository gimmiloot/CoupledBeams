"""Targeted continuation diagnostics of the preserved four-field time pilot.

This module reads historical manifests without requiring their code identity
to match a later numerical implementation. Physical reconstruction and L2
projection use Shen polynomials independently of the dynamics evaluator.
Importing it performs no file reads, environment changes, symbolic work or
integration. Orchestration calls the diagnostic functions explicitly.
"""
from __future__ import annotations

import hashlib
import json
import argparse
import importlib.metadata
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

if __name__ == "__main__":
    for _variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS"):
        os.environ[_variable] = "1"

import numpy as np
from numpy.polynomial.legendre import leggauss, legvander
from scipy.linalg import cho_factor, cho_solve, solve_triangular


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
SOURCE = ROOT / "results/weakly_nonlinear_planar_time_pilot/c97287772bc461ef"
OUTPUT = ROOT / "results/weakly_nonlinear_planar_recovery"
FIELDS = ("u", "w", "theta", "c")
VERSION = "nlsp-planar-targeted-recovery-v1"


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False,
                                   allow_nan=False) + "\n", encoding="utf8")
    temporary.replace(path)


def _read_json(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def raw_basis(p, points, length=1.):
    """Unscaled essential-only Shen basis, n=0,...,p-2."""
    if p < 2 or length <= 0:
        raise ValueError("A positive length and p>=2 are required")
    vandermonde = legvander(2 * np.asarray(points) / length - 1, p)
    return vandermonde[:, :p-1] - vandermonde[:, 2:p+1]


def physical_transforms(p, coefficients, length=1.):
    """Invert the original constant-resting-mass coordinate whitening."""
    coefficients = coefficients.values() if hasattr(coefficients, "values") and not isinstance(coefficients, dict) else coefficients
    nodes, weights = leggauss(2 * p + 1)
    weights = weights * length / 2
    basis = raw_basis(p, (nodes+1) * length / 2, length)
    gram = basis.T @ (weights[:, None] * basis)
    masses = (coefficients["m"], coefficients["m"], coefficients["jp"], coefficients["jp"])
    return tuple(solve_triangular(np.linalg.cholesky(mass * gram).T,
                                 np.eye(p-1), lower=False)
                 for mass in masses)


def physical_bases(p, points, coefficients, length=1.):
    raw = raw_basis(p, points, length)
    return tuple(raw @ transform for transform in physical_transforms(p, coefficients, length))


def projection_operator(lo, hi, coefficients, length=1., nq=100):
    """Four physical L2 orthogonal projection maps from high to low space."""
    if nq < max(lo, hi) + 1:
        raise ValueError("Projection quadrature must integrate the polynomial products")
    nodes, weights = leggauss(nq)
    points = (nodes+1) * length / 2
    weights = weights * length / 2
    low, high = physical_bases(lo, points, coefficients, length), physical_bases(hi, points, coefficients, length)
    return tuple(cho_solve(cho_factor(left.T @ (weights[:, None] * left), lower=True),
                           left.T @ (weights[:, None] * right))
                 for left, right in zip(low, high))


def _case_path(source, p):
    return Path(source) / "cases" / f"p{p}_Aoverh0p05_tight"


def validate_historical(source=SOURCE):
    """Check the historical artifact manifest and actual, not nominal, times."""
    source = Path(source)
    if not (source / "manifest.json").is_file():
        raise FileNotFoundError(f"MISSING_TRAJECTORY_DATA: {source}")
    manifest = _read_json(source / "manifest.json")
    checks = {}
    for relative, expected in manifest["artifact_hashes"].items():
        path = source / relative
        if not path.is_file():
            raise FileNotFoundError(f"MISSING_TRAJECTORY_DATA: {path}")
        actual = sha(path)
        checks[relative] = {"expected": expected, "actual": actual, "match": actual == expected}
        if actual != expected:
            raise ValueError(f"Historical artifact hash mismatch: {path}")
    summary = _read_json(source / "summary.json")
    coefficients = summary["coefficients"]
    length = summary["config"]["material_geometry"]["L"]
    cases = {}
    keys = ("status", "failure", "p", "ndof", "nq", "amplitude_over_h", "amplitude",
            "time_level", "rtol", "atol", "max_step", "time_end", "target_time_end", "samples",
            "accepted_internal_steps", "min_internal_step", "max_internal_step", "nfev", "njev", "nlu",
            "integration_seconds", "counters")
    for folder in sorted((source / "cases").iterdir()):
        metadata = _read_json(folder / "case.json")
        with np.load(folder / "trajectory.npz") as data:
            times = data["time"]
            intervals = np.diff(times)
            if len(times) < 2 or not np.all(np.isfinite(times)) or np.any(intervals <= 0):
                raise ValueError(f"Historical time grid is invalid: {folder}")
            if len(times) != metadata["samples"] or times[-1] != metadata["time_end"]:
                raise ValueError(f"Actual historical times disagree with metadata: {folder}")
            indices = data["snapshot_indices"]
            if np.any(indices < 0) or np.any(indices >= len(times)):
                raise ValueError(f"Snapshot indices exceed actual history: {folder}")
            row = {key: metadata[key] for key in keys if key in metadata}
            row.update({"array_samples": len(times), "actual_array_time_end": float(times[-1]),
                        "monotone_time": True, "saved_output_dt_min": float(intervals.min()),
                        "saved_output_dt_max": float(intervals.max()),
                        "rhs_per_accepted_step": metadata["nfev"] / metadata["accepted_internal_steps"],
                        "internal_max_step_fraction_of_max_step": metadata["max_internal_step"] / metadata["max_step"],
                        "full_internal_step_distribution_available": False,
                        "actual_snapshot_times": times[indices].tolist(), "actual_snapshot_indices": indices.tolist(),
                        "duplicated_final_snapshots": len(np.unique(indices)) != len(indices)})
            p = metadata["p"]
            coordinate = data["q"][indices]
            raw = np.column_stack([coordinate[:, i*(p-1):(i+1)*(p-1)] @ transform.T
                                   for i, transform in enumerate(physical_transforms(p, coefficients, length))])
            row["whitening_inverse_raw_snapshot_max_difference"] = float(np.max(abs(raw-data["raw_snapshot_coefficients"])))
            profiles = np.stack([raw[:, i*(p-1):(i+1)*(p-1)] @ raw_basis(p, data["snapshot_points"], length).T
                                 for i in range(4)], axis=2)
            row["independent_reconstruction_snapshot_max_difference"] = float(np.max(abs(profiles-data["snapshots"])))
            cases[folder.name] = row
    return {"source_manifest_sha256": sha(source / "manifest.json"), "source_bundle": str(source),
            "artifact_hash_validation": checks, "cases": cases,
            "integrations": 0, "symbolic_derivations": 0,
            "qualification": "Historical manifest verified on its own terms; current numerical code hash need not match historical identity"}


def diagnostic_pair(source, output, lo=24, hi=32, high_case=None):
    """Reproduce full-interval norms and orthogonally split the physical error.

    high_case may be a new p48 case folder. Its common trajectory must have
    exactly the old time grid; no interpolation, phase shift or scaling is
    allowed. Extra high-resolution samples belong in a separate artifact.
    """
    started = time.perf_counter()
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    summary = _read_json(source / "summary.json")
    config, coefficients = summary["config"], summary["coefficients"]
    length, T1 = config["material_geometry"]["L"], summary["initial_eigenpair"]["T1"]
    tag = f"p{lo}_p{hi}"
    low_folder, high_folder = _case_path(source, lo), Path(high_case) if high_case is not None else _case_path(source, hi)
    low_data, high_data = np.load(low_folder / "trajectory.npz"), np.load(high_folder / "trajectory.npz")
    try:
        times = low_data["time"]
        if not np.array_equal(times, high_data["time"]):
            raise ValueError("Unequal actual time arrays: comparison cannot interpolate or extend a partial history")
        if times[-1] != 5 * T1:
            raise ValueError("Full 0...5T1 spatial diagnostic requires completed histories")
        nq = config["spatial"]["comparison_quadrature"]
        nodes, weights = leggauss(nq)
        points, weights = (nodes+1)*length/2, weights*length/2
        low_basis, high_basis = physical_bases(lo, points, coefficients, length), physical_bases(hi, points, coefficients, length)
        projectors = projection_operator(lo, hi, coefficients, length, nq)
        zone_order = max(40, hi+1)
        zones = ((0., .1*length), (.1*length, .9*length), (.9*length, length))
        zone_points, zone_weights = [], []
        for left, right in zones:
            nodes, weights_zone = leggauss(zone_order)
            zone_points.extend((left+(nodes+1)*(right-left)/2).tolist())
            zone_weights.append(weights_zone*(right-left)/2)
        zone_points = np.asarray(zone_points)
        zone_low, zone_high = physical_bases(lo, zone_points, coefficients, length), physical_bases(hi, zone_points, coefficients, length)
        arrays = {"time": times, "comparison_points": points, "zone_bounds": np.asarray(zones)}
        rows, profile_payload = {}, {}
        dense = np.linspace(0, length, 1001)
        dense_low, dense_high = physical_bases(lo, dense, coefficients, length), physical_bases(hi, dense, coefficients, length)
        for part in ("q", "velocity"):
            low_coordinates, high_coordinates = low_data[part], high_data[part]
            if low_coordinates.shape != (len(times), 4*(lo-1)) or high_coordinates.shape != (len(times), 4*(hi-1)):
                raise ValueError("Coefficient history dimension differs from four-field Shen degree")
            shape = (len(times), 4)
            norms, tail_norms, common_norms, refs, max_difference, reference_max, peak_x, cross, identity = [np.empty(shape) for _ in range(9)]
            zone_squared = np.empty((len(times), 4, 3))
            for start in range(0, len(times), 1024):
                stop = min(start+1024, len(times))
                for field in range(4):
                    sl, sh = slice(field*(lo-1), (field+1)*(lo-1)), slice(field*(hi-1), (field+1)*(hi-1))
                    a, b = low_coordinates[start:stop, sl], high_coordinates[start:stop, sh]
                    aa, bb = a @ low_basis[field].T, b @ high_basis[field].T
                    projected_coordinates = b @ projectors[field].T
                    difference, tail = bb-aa, bb-projected_coordinates @ low_basis[field].T
                    common = (projected_coordinates-a) @ low_basis[field].T
                    def squared(values):
                        return np.einsum("ti,i,ti->t", values, weights, values)
                    error_sq, tail_sq, common_sq = squared(difference), squared(tail), squared(common)
                    norms[start:stop, field], tail_norms[start:stop, field], common_norms[start:stop, field] = np.sqrt(error_sq), np.sqrt(tail_sq), np.sqrt(common_sq)
                    refs[start:stop, field], reference_max[start:stop, field] = np.sqrt(squared(bb)), np.max(abs(bb), axis=1)
                    indices = np.argmax(abs(difference), axis=1)
                    max_difference[start:stop, field] = abs(difference[np.arange(stop-start), indices])
                    peak_x[start:stop, field] = points[indices]
                    cross[start:stop, field] = np.einsum("ti,i,ti->t", tail, weights, common)
                    identity[start:stop, field] = error_sq-tail_sq-common_sq
                    difference_zone = b @ zone_high[field].T-a @ zone_low[field].T
                    for zone in range(3):
                        values = difference_zone[:, zone_order*zone:zone_order*(zone+1)]
                        zone_squared[start:stop, field, zone] = np.einsum("ti,i,ti->t", values, zone_weights[zone], values)
            for field, name in enumerate(FIELDS):
                key = part+"_"+name
                scale_l2, scale_max = float(refs[:, field].max()), float(reference_max[:, field].max())
                floor = config["gates"]["relative_numerical_floor"]*max(float(refs.max()), 1e-30)
                peak_l2, peak_max = int(np.argmax(norms[:, field])), int(np.argmax(max_difference[:, field]))
                tolerance = config["gates"]["w_theta_relative" if name in ("w", "theta") else "u_c_relative"]
                former = summary["comparisons"].get("spatial_"+tag, {}).get("fields", {}).get(key)
                temporal = summary["comparisons"].get("temporal_medium_tight", {}).get("fields", {}).get(key)
                total_integral = float(np.trapezoid(norms[:, field]**2, times))
                tail_integral = float(np.trapezoid(tail_norms[:, field]**2, times))
                common_integral = float(np.trapezoid(common_norms[:, field]**2, times))
                relative_l2 = float(norms[:, field].max()/max(scale_l2, floor))
                relative_max = float(max_difference[:, field].max()/max(scale_max, floor))
                row = {"relative_L2": relative_l2, "relative_max": relative_max,
                       "absolute_L2": float(norms[:, field].max()), "absolute_max": float(max_difference[:, field].max()),
                       "reference_scale_L2": scale_l2, "reference_scale_max": scale_max,
                       "reporting_floor": floor, "tolerance": tolerance,
                       "pass": relative_l2 <= tolerance and relative_max <= tolerance,
                       "L2_peak_time": float(times[peak_l2]), "L2_peak_t_over_T1": float(times[peak_l2]/T1),
                       "max_peak_time": float(times[peak_max]), "max_peak_t_over_T1": float(times[peak_max]/T1),
                       "max_peak_x_quadrature_node": float(peak_x[peak_max, field]),
                       "integrated_tail_fraction": tail_integral/max(total_integral, 1e-60),
                       "integrated_common_fraction": common_integral/max(total_integral, 1e-60),
                       "max_tail_L2": float(tail_norms[:, field].max()), "max_common_L2": float(common_norms[:, field].max()),
                       "pythagoras_max_absolute": float(np.max(abs(identity[:, field]))),
                       "pythagoras_scaled_by_field_characteristic_sq": float(np.max(abs(identity[:, field]))/max(scale_l2**2, 1e-60)),
                       "orthogonality_max_absolute": float(np.max(abs(cross[:, field]))),
                       "orthogonality_scaled_by_field_characteristic_sq": float(np.max(abs(cross[:, field]))/max(scale_l2**2, 1e-60)),
                       "zone_time_integrated_fractions": (np.trapezoid(zone_squared[:, field, :], times, axis=0)/max(total_integral, 1e-60)).tolist(),
                       "zone_fractions_at_L2_peak": (zone_squared[peak_l2, field, :]/max(norms[peak_l2, field]**2, 1e-60)).tolist(),
                       "time_windows": {}}
                if former is not None:
                    row["historical_relative_L2_absolute_difference"] = abs(relative_l2-former["relative_L2"])
                    row["historical_relative_max_absolute_difference"] = abs(relative_max-former["relative_max"])
                if temporal is not None:
                    row["temporal_medium_tight_L2"] = temporal["relative_L2"]
                    row["spatial_to_temporal_L2_ratio"] = relative_l2/max(temporal["relative_L2"], 1e-30)
                    row["temporal_comparison_qualification"] = "historical p32 temporal evidence only; not full-horizon temporal evidence for p48"
                for label, left, right in (("initial_0_0p01T1", 0., .01), ("early_0_0p1T1", 0., .1), ("first_period", 0., 1.), ("last_period", 4., 5.)):
                    mask = (times >= left*T1) & (times <= right*T1)
                    row["time_windows"][label] = {"max_relative_L2": float(norms[mask, field].max()/max(scale_l2, floor)),
                        "max_L2_fraction_of_whole_max": float(norms[mask, field].max()/max(norms[:, field].max(), 1e-60)),
                        "tail_fraction_time_integrated": float(np.trapezoid(tail_norms[mask, field]**2, times[mask])/max(np.trapezoid(norms[mask, field]**2, times[mask]), 1e-60))}
                excess = np.flatnonzero(norms[:, field]/max(scale_l2, floor) > tolerance)
                row["first_tolerance_exceedance_time"] = None if not len(excess) else float(times[excess[0]])
                row["first_tolerance_exceedance_t_over_T1"] = None if not len(excess) else float(times[excess[0]]/T1)
                rows[key] = row
                selected = np.unique(np.r_[0, np.argmin(abs(times-.001*T1)), np.argmin(abs(times-.01*T1)), np.argmin(abs(times-.1*T1)), np.argmin(abs(times-T1)), peak_l2, peak_max, len(times)-1]).astype(int)
                sl, sh = slice(field*(lo-1), (field+1)*(lo-1)), slice(field*(hi-1), (field+1)*(hi-1))
                a, b = low_coordinates[selected, sl], high_coordinates[selected, sh]
                projected = b @ projectors[field].T
                profile_payload.update({key+"_times": times[selected], key+"_x": dense,
                    key+"_difference": b @ dense_high[field].T-a @ dense_low[field].T,
                    key+"_tail": b @ dense_high[field].T-projected @ dense_low[field].T,
                    key+"_common": (projected-a) @ dense_low[field].T})
            arrays.update({part+"_difference_L2": norms, part+"_difference_max": max_difference,
                           part+"_max_x": peak_x, part+"_tail_L2": tail_norms, part+"_common_L2": common_norms,
                           part+"_reference_L2": refs, part+"_zone_squared_L2": zone_squared,
                           part+"_orthogonality": cross, part+"_pythagoras": identity})
            del low_coordinates, high_coordinates
        arrays.update(profile_payload)
        temporary = output / (tag+"_arrays.tmp.npz")
        np.savez_compressed(temporary, **arrays)
        temporary.replace(output / (tag+"_arrays.npz"))
        result = {"status": "PASS" if all(row["pass"] for row in rows.values()) else "PARTIAL",
                  "pair": [lo, hi], "source_bundle": str(source), "high_case": str(high_folder),
                  "common_time_array_identical": True, "samples": len(times),
                  "time_start": float(times[0]), "time_end": float(times[-1]), "L2_quadrature": nq,
                  "zones_quadrature_each": zone_order, "zone_bounds": zones,
                  "maximum_definition": "historical maximum over fixed100-point Gauss nodes and common saved times; additional1001-node profiles are diagnostic only",
                  "tail_fraction_definition": "trapezoidal time integral of squared physical L2 tail divided by the same integral of total squared difference; not energy classification",
                  "fields": rows, "processing_seconds": time.perf_counter()-started,
                  "integrations": 0, "symbolic_derivations": 0}
        write_json(output / (tag+"_diagnostic.json"), result)
        return result
    finally:
        low_data.close()
        high_data.close()


def initial_compatibility(source=SOURCE):
    """Continuous initial acceleration traces from the protected cubic action."""
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.lib import mindlin_herrmann_longitudinal as mh
    from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section

    summary = _read_json(Path(source) / "summary.json")
    config, geometry = summary["config"], summary["config"]["material_geometry"]
    audit = _read_json(ROOT / config["audit_bundle"] / "result.json")
    coefficients = rod.RodCoefficients(**audit["coefficients"])
    reference = _read_json(ROOT / config["linear_reference_bundle"] / "result.json")
    omega = reference["timoshenko"]["roots"][0]["omega"]
    length = geometry["L"]
    section = rectangular_section(E=geometry["E"], rho=geometry["rho"], nu=geometry["nu"],
                                  width=geometry["b"], thickness=geometry["h"], K=5/6)
    model = mh.project_jang_reduced_rectangular(section)
    mode = mh.finite_mode(model, length, omega, "timoshenko")
    peak = (mh.finite_state_basis(model, length, omega, [length/2], "timoshenko") @ mode["coefficients"])[0, 0]
    normalized = mode["coefficients"]/peak
    points = np.array([0., length])
    values = mh.finite_state_basis(model, length, omega, points, "timoshenko") @ normalized
    first = mh.finite_state_basis(model, length, omega, points, "timoshenko", 1) @ normalized
    state_matrix = mh.harmonic_state_matrix(model, omega, "timoshenko")
    second = first @ state_matrix.T
    symbolic = rod.derive_polynomials()
    inactive = {name+suffix: 0 for name in ("v", "Phi", "psi")
                for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")}
    potential = symbolic.V4.substitute(inactive)
    flux = potential.derivative("u_s")
    symbols = symbolic.symbols
    expected = ((symbols["C"]-symbols["S"])*symbols["theta"]*symbols["w_s"]
                +(symbols["S"]-symbols["C"]/2)*symbols["theta"]**2)
    second_order_flux = flux.homogeneous(2).substitute({"u_s": 0, "c": 0})
    body = [symbolic.residual_a[index].substitute(inactive) for index in (0, 1, 5, 6)]
    endpoint = {name+suffix: 0 for name in FIELDS for suffix in ("", "_s", "_t", "_ss", "_st", "_tt")
                if not (name in ("w", "theta") and suffix in ("_s", "_ss"))}
    endpoint_body = [residual.substitute(endpoint) for residual in body]
    checks = {"Fu2_exact_match": second_order_flux == expected,
              "endpoint_u_body_exact_match": endpoint_body[0] == -(symbols["C"]-symbols["S"])*symbols["theta_s"]*symbols["w_s"],
              "endpoint_w_body": str(endpoint_body[1]), "endpoint_theta_body": str(endpoint_body[2]),
              "endpoint_c_body": str(endpoint_body[3])}
    if not checks["Fu2_exact_match"] or not checks["endpoint_u_body_exact_match"]:
        raise ArithmeticError("Initial compatibility expression does not match protected action")
    records = []
    active = np.array([0, 1, 5, 6])
    for amplitude in (.0025, .00125):
        for side, coordinate in enumerate(points):
            q, qs, qss = np.zeros(7), np.zeros(7), np.zeros(7)
            q[[1, 5]], qs[[1, 5]], qss[[1, 5]] = amplitude*values[side, :2], amplitude*first[side, :2], amplitude*second[side, :2]
            jet = rod.FieldJet(q, qs, np.zeros(7), qss, np.zeros(7), np.zeros(7))
            variables = jet.values() | coefficients.values()
            residual = np.array([expression.evaluate(variables) for expression in body])
            linear = np.array([expression.homogeneous(1).evaluate(variables) for expression in body])
            mass = np.array([coefficients.m, coefficients.m, coefficients.jp*(1+q[6])**2, coefficients.jp])
            acceleration, linear_acceleration = -residual/mass, -linear/mass
            leading = (coefficients.C-coefficients.S)/coefficients.m*qs[5]*qs[1]
            records.append({"A": amplitude, "s": float(coordinate), "q_end_actual": q[active].tolist(),
                "q_s_end": qs[active].tolist(), "q_ss_end": qss[active].tolist(),
                "linear_acceleration_trace": linear_acceleration.tolist(),
                "linear_exact_mode_acceleration": (-omega**2*q[active]).tolist(),
                "cubic_acceleration_trace": acceleration.tolist(),
                "analytical_trace_with_exact_clamp_and_continuous_linear_mode_identities": [leading, 0., 0., 0.],
                "axial_A2_expression": leading, "nonlinear_minus_linear": (acceleration-linear_acceleration).tolist()})
    provenance = {}
    for relative in (config["audit_bundle"], config["linear_reference_bundle"]):
        bundle = ROOT / relative
        manifest = _read_json(bundle / "manifest.json")
        hashes = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
        verified = {name: (bundle/name).is_file() and sha(bundle/name) == digest for name, digest in hashes.items()}
        if not verified or not all(verified.values()):
            raise ValueError(f"Protected initial reference hash mismatch: {bundle}")
        provenance[relative] = {"manifest_sha256": sha(bundle/"manifest.json"), "artifact_hashes_verified": verified}
    return {"status": "CONFIRMED_LOW_ORDER_MISMATCH", "fields": list(FIELDS), "model_version": rod.MODEL_VERSION,
        "omega": omega, "T1": 2*np.pi/omega, "normalization": "same common signed midpoint factor as original pilot",
        "coefficients": coefficients.values(), "Vp": str(potential), "Fu": str(flux),
        "Fu2_at_u_s_c_zero": str(second_order_flux), "expected_Fu2": str(expected),
        "endpoint_residuals": {field: str(expression) for field, expression in zip(FIELDS, endpoint_body)},
        "symbolic_checks": checks,
        "shape_derivative_source": "continuous analytical finite_state_basis values/first derivative; second derivative from verified harmonic state ODE, no Galerkin reconstruction",
        "shape_linear_state_consistency_max": float(np.max(abs(first-values@state_matrix.T))),
        "unit_amplitude_endpoint_derivatives": {"values": values[:, :2].tolist(), "first": first[:, :2].tolist(), "second": second[:, :2].tolist()},
        "records": records, "hash_provenance": provenance,
        "qualification": "Nonzero internal strong axial acceleration trace is A^2 while fixed Dirichlet requires zero for time-C2 smoothness up to boundary; not invalidity of weak IVP, not a missing slope constraint, and not the sole demonstrated cause of spatial differences",
        "amplitude_halving_axial_trace_ratio": records[0]["axial_A2_expression"]/records[2]["axial_A2_expression"],
        "runtime_integration_seconds": 0, "root_solves": 0}


def recovery_identity(source, baseline, mode):
    source = Path(source)
    historical = _read_json(source/"manifest.json")
    historical_hashes = {name.replace(chr(92), "/"): digest
                         for name, digest in historical["identity"]["hashes"].items()}
    original_hash = historical_hashes["scripts/lib/weakly_nonlinear_planar_dynamics.py"]
    if baseline is not None and sha(baseline) != original_hash:
        raise ValueError("Baseline helper does not match historical execution SHA256")
    item = {"version": VERSION, "mode": mode, "historical_manifest_sha256": sha(source/"manifest.json"),
            "historical_source": str(source.resolve()), "baseline_sha256": original_hash,
            "files": {p:sha(ROOT/p) for p in ("scripts/analysis/diagnose_weakly_nonlinear_planar_rod.py",
                      "scripts/lib/weakly_nonlinear_planar_dynamics.py", "scripts/lib/weakly_nonlinear_spatial_rod.py",
                      "scripts/analysis/simulate_weakly_nonlinear_planar_rod.py", "data/input/weakly_nonlinear_planar_time_pilot.json")},
            "policy": {"budget_seconds":900, "max_short_integrations":3, "max_full_integrations":1,
                       "p":48,"amplitude_over_h":.05,"full_time_level":"tight","short_p48_time_level":"allowed_extra",
                       "forecast_factor":1.25,"historical_data_recomputed":False,"previous_small_amplitude_unchanged":True},
            "dependencies": {n:importlib.metadata.version(n) for n in ("numpy","scipy","matplotlib")}}
    key=hashlib.sha256(json.dumps(item,sort_keys=True).encode()).hexdigest()[:16]
    return key,item


def validate_recovery(bundle, identity=None):
    bundle=Path(bundle);manifest=_read_json(bundle/"manifest.json")
    if identity is not None and manifest["identity"]!=identity:
        raise ValueError("Recovery cache identity mismatch")
    for name,digest in manifest["artifact_hashes"].items():
        if sha(bundle/name)!=digest: raise ValueError(f"Recovery artifact hash mismatch: {name}")
    return _read_json(bundle/"summary.json")


def load_baseline(path):
    spec=importlib.util.spec_from_file_location("nlsp_frozen_numerical_baseline",path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def profile_equivalence(old, new, source):
    """Different real states, warmup, median; no repeated one-state cache trick."""
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    with np.load(_case_path(source,32)/"trajectory.npz") as data:
        ids=np.linspace(0,len(data["time"])-1,16,dtype=int)
        q,v=data["q"][ids],data["velocity"][ids]
    started=time.perf_counter();timings={};equivalence={}
    for label,disc in (("old",old),("new",new)):
        operations={"energy_only":lambda a,b:disc.potential(a,gradient=False),
                    "potential_gradient":lambda a,b:disc.potential(a),
                    "potential_hessian":lambda a,b:disc.potential(a,hessian=True),
                    "mass_assembly":lambda a,b:disc._theta_mass(a),
                    "mass_solve":lambda a,b:disc._solve_mass(a,b),
                    "inertial_terms":disc.inertial_terms,
                    "safety":lambda a,b:pilot.safety_check(disc,a,_read_json(source/"summary.json")["config"]["safety"]),
                    "rhs":lambda a,b:disc.rhs(0,np.r_[a,b]),
                    "jacobian":lambda a,b:disc.jacobian(0,np.r_[a,b]),
                    "reconstruction":lambda a,b:disc.reconstruct(a)}
        # Parse immutable configuration once rather than inside a timed safety call.
        safety=_read_json(source/"summary.json")["config"]["safety"]
        operations["safety"]=lambda a,b:pilot.safety_check(disc,a,safety)
        timings[label]={}
        for name,operation in operations.items():
            for a,b in zip(q,v):operation(a,b)
            repeats=[]
            for _ in range(3):
                begin=time.perf_counter()
                for a,b in zip(q,v):operation(a,b)
                repeats.append((time.perf_counter()-begin)/len(q))
            timings[label][name]={"median_seconds":float(np.median(repeats)),"repeats_seconds":repeats,"different_states":len(q)}
    for name,operation in {"V":lambda d,a,b:d.potential(a,gradient=False)["V"],
                           "gradient":lambda d,a,b:d.potential(a)["gradient"],
                           "hessian":lambda d,a,b:d.potential(a,hessian=True)["hessian"],
                           "mass":lambda d,a,b:d.mass_matrix(a),"inertia":lambda d,a,b:d.inertial_terms(a,b),
                           "acceleration":lambda d,a,b:d.acceleration(a,b),"rhs":lambda d,a,b:d.rhs(0,np.r_[a,b]),
                           "jacobian":lambda d,a,b:d.jacobian(0,np.r_[a,b]),"energy":lambda d,a,b:d.energy(a,b),
                           "energy_rate":lambda d,a,b:d.energy_rate(a,b),
                           "weak_residual":lambda d,a,b:d.weak_residual(a,b,d.acceleration(a,b))}.items():
        absolute=relative=0.;bitwise=True
        for a,b in zip(q,v):
            x,y=np.asarray(operation(old,a,b)),np.asarray(operation(new,a,b))
            absolute=max(absolute,float(np.max(np.abs(x-y))));relative=max(relative,float(np.linalg.norm(x-y)/max(np.linalg.norm(x),1e-30)))
            bitwise=bitwise and np.array_equal(x,y)
        equivalence[name]={"absolute":absolute,"relative":relative,"bitwise_equal":bitwise}
    return {"status":"PASS" if all(x["relative"]<=2e-12 for x in equivalence.values()) else "FAIL",
            "equivalence":equivalence,"timings":timings,"profiling_wall_seconds":time.perf_counter()-started,
            "threads":{n:os.environ.get(n) for n in ("OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","OMP_NUM_THREADS")},
            "policy":"16 different real states; warmup then3 repeats, median; no ODE in profiling"}


def save_control(folder,history,stats,times):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    n=history.shape[1]//2;temporary=folder/"trajectory.tmp.npz"
    np.savez_compressed(temporary,time=times[:len(history)],q=history[:,:n],velocity=history[:,n:])
    temporary.replace(folder/"trajectory.npz");write_json(folder/"case.json",stats)


def run_recovery(source,bundle,baseline_path,item,compute):
    """Exactly two p32 short runs + one stricter p48 short + at most one full."""
    from scripts.analysis import simulate_weakly_nonlinear_planar_rod as pilot
    from scripts.lib import weakly_nonlinear_planar_dynamics as dynamics
    pilot.load_runtime()
    inventory=validate_historical(source);write_json(bundle/"source_inventory.json",inventory)
    pair16=diagnostic_pair(source,bundle,16,24);pair24=diagnostic_pair(source,bundle,24,32)
    compatibility=initial_compatibility(source);write_json(bundle/"initial_compatibility.json",compatibility)
    statuses={"NLSP_PLANAR_ERROR_DIAGNOSTIC":"PASS","NLSP_PLANAR_INITIAL_COMPATIBILITY":compatibility["status"],
              "NLSP_PLANAR_RHS_EQUIVALENCE":"NOT_RUN","NLSP_PLANAR_PERFORMANCE":"NOT_RUN",
              "NLSP_PLANAR_P48_SPATIAL_CHECK":"NOT_RUN","NLSP_PLANAR_SOLVER_RECOVERY":"PARTIAL"}
    result={"statuses":statuses,"source":str(source),"source_inventory":inventory,"diagnostics":{"p16_p24":pair16,"p24_p32":pair24},
            "compatibility":compatibility,"integrations":0,"short_integrations":0,"full_integrations":0,"budget_limit_seconds":900,
            "budget_used_seconds":0.,"no_small_amplitude_rerun":True,"previous_pilot_status":"PARTIAL"}
    if not compute:return result
    if baseline_path is None:raise ValueError("--compute requires the preserved historical --baseline-helper")
    baseline_folder=bundle/"baseline";baseline_folder.mkdir(exist_ok=True)
    copied=baseline_folder/"weakly_nonlinear_planar_dynamics.py";shutil.copyfile(baseline_path,copied)
    baseline=load_baseline(copied)
    historical=_read_json(source/"summary.json");config=historical["config"]
    coefficients,shape,initial,_=pilot.setup(config)
    old=baseline.PlanarGalerkin(coefficients,32,length=1.)
    new=dynamics.PlanarGalerkin(coefficients,32,length=1.,model=old.model)
    prof=profile_equivalence(old,new,source);write_json(bundle/"performance_equivalence.json",prof)
    # Includes earlier isolated profiler measurements, charged conservatively.
    spent=10.+prof["profiling_wall_seconds"];result["budget_used_seconds"]=spent
    statuses["NLSP_PLANAR_RHS_EQUIVALENCE"]=prof["status"]
    statuses["NLSP_PLANAR_PERFORMANCE"]="PASS" if prof["timings"]["new"]["rhs"]["median_seconds"]<prof["timings"]["old"]["rhs"]["median_seconds"] else "NO_GAIN"
    result["performance_equivalence"]=prof
    if prof["status"]!="PASS":return result
    p48=dynamics.PlanarGalerkin(coefficients,48,length=1.,nq=97,model=old.model)
    with np.load(_case_path(source,32)/"trajectory.npz") as data:comparison_times=data["time"]
    short_end=.1*initial["T1"];end=5*initial["T1"]
    omega_max=float(p48.linear_eigenpairs()["omega"][-1]);dt=2*np.pi/omega_max/12
    fine=np.linspace(0,end,int(np.ceil(end/dt))+1)
    master=np.unique(np.r_[comparison_times,fine,short_end]);short_times=master[master<=short_end]
    controls={};histories={}
    for name,disc,level in (("old_p32",old,"tight"),("new_p32",new,"tight"),("p48_strict_short",p48,"allowed_extra")):
        if spent>=900:break
        history,stats=pilot.integrate_case(disc,shape,initial,config,.05,level,short_times,time.perf_counter()+900-spent)
        spent+=stats["integration_seconds"];result["integrations"]+=1;result["short_integrations"]+=1
        stats["initial_energy"]=disc.energy(history[0,:disc.ndof],history[0,disc.ndof:])
        energies=np.array([disc.energy(row[:disc.ndof],row[disc.ndof:]) for row in history])
        stats["max_relative_energy_drift"]=float(np.max(abs(energies/energies[0]-1)))
        save_control(bundle/"controls"/name,history,stats,short_times)
        controls[name]=stats;histories[name]=history
        write_json(bundle/"budget_progress.json",{"spent_seconds":spent,"short_integrations":result["short_integrations"],"controls":controls})
        print(json.dumps({"control":name,"status":stats["status"],"wall_seconds":stats["integration_seconds"],"budget_used":spent}),flush=True)
        if stats["status"]!="PASS":break
    result.update({"controls":controls,"budget_used_seconds":spent})
    if set(controls)!={"old_p32","new_p32","p48_strict_short"}:return result
    if any(control["status"]!="PASS" for control in controls.values()):
        statuses["NLSP_PLANAR_P48_SPATIAL_CHECK"]="NOT_RUN_FAILED_SHORT_CONTROL"
        return result
    def series(history,disc):return {"q":history[:,:disc.ndof],"velocity":history[:,disc.ndof:]}
    short_equivalence=pilot.compare_histories(old,series(histories["old_p32"],old),new,series(histories["new_p32"],new),config)
    result["short_implementation_equivalence"]=short_equivalence;write_json(bundle/"short_equivalence.json",short_equivalence)
    if short_equivalence["status"]!="PASS":return result
    forecast=controls["p48_strict_short"]["integration_seconds"]*50*1.25
    decision={"forecast_full_seconds":forecast,"remaining_integration_budget_seconds":900-spent,"forecast_from":"stricter p48 0.1T1 short, same initial fields; factor50*horizon and1.25 cost guard",
              "p48_omega_max":omega_max,"old_comparison_dt":float(np.max(np.diff(comparison_times))),"fine_output_dt":float(fine[1]),
              "master_samples":len(master),"old_samples_per_p48_fastest_linear_period":2*np.pi/omega_max/float(np.max(np.diff(comparison_times))),
              "sampling_qualification":"p32 saved times are exact common comparison nodes; p48 extra uniform nodes resolve fastest retained linear period with12 samples; no old trajectory interpolation"}
    result["p48_decision"]=decision;write_json(bundle/"p48_decision.json",decision)
    if forecast>900-spent:
        statuses["NLSP_PLANAR_P48_SPATIAL_CHECK"]="REFINEMENT_DEFERRED_BY_BUDGET";return result
    print(json.dumps({"p48":"full integration authorized by gates and fixed budget","forecast":forecast,"remaining":900-spent}),flush=True)
    history,stats=pilot.integrate_case(p48,shape,initial,config,.05,"tight",master,time.perf_counter()+900-spent)
    spent+=stats["integration_seconds"];result["integrations"]+=1;result["full_integrations"]+=1
    # Preserve both the original comparison times and the fine p48 observation grid.
    p48folder=bundle/"p48_tight";p48folder.mkdir(exist_ok=True)
    present=master[:len(history)];ids=np.searchsorted(present,comparison_times[comparison_times<=present[-1]])
    comparison_history=history[ids]
    save_control(p48folder,comparison_history,{**stats,"samples":len(ids),"time_end":float(present[ids[-1]])},present[ids])
    temporary=p48folder/"fine_trajectory.tmp.npz"
    np.savez_compressed(temporary,time=present,q=history[:,:p48.ndof],velocity=history[:,p48.ndof:],comparison_indices=ids)
    temporary.replace(p48folder/"fine_trajectory.npz")
    short_count=min(len(short_times),len(history));reference=histories["p48_strict_short"][:short_count]
    short_time_check=pilot.compare_histories(p48,series(reference,p48),p48,series(history[:short_count],p48),config)
    result["p48_short_temporal_check"]={**short_time_check,"time_end":float(present[short_count-1]),"qualification":"stricter-vs-tight on0..0.1T1 only; full p48 temporal convergence is not independently established"}
    diagnostic,energy,drift,_,_=pilot.series_measure(p48,comparison_history[:,:p48.ndof],comparison_history[:,p48.ndof:],config)
    stats["diagnostics"]=diagnostic;stats["comparison_samples"]=len(ids);stats["master_samples"]=len(history)
    stats["integration_sampling_samples"]=stats["samples"];stats["samples"]=len(ids)
    stats["master_time_end"]=float(present[-1]);stats["time_end"]=float(present[ids[-1]])
    np.savez_compressed(p48folder/"observations.npz",time=present[ids],energy=energy,energy_drift=drift,
                        observations=p48.reconstruct_series(comparison_history[:,:p48.ndof],np.array([.25,.5])))
    write_json(p48folder/"case.json",stats)
    result.update({"p48":stats,"budget_used_seconds":spent})
    if stats["status"]=="PASS":
        last=diagnostic_pair(source,bundle,32,48,high_case=p48folder);result["p32_p48"]=last
        passed=all(row["relative_L2"]<=row["tolerance"] and row["relative_max"]<=row["tolerance"] for row in last["fields"].values())
        statuses["NLSP_PLANAR_P48_SPATIAL_CHECK"]="PASS" if passed and diagnostic["max_relative_energy_drift"]<=config["gates"]["energy_relative_drift"] else "PARTIAL"
    else:statuses["NLSP_PLANAR_P48_SPATIAL_CHECK"]="PARTIAL_BUDGET"
    # The incomplete second-amplitude neighboring-p control is deliberately retained.
    statuses["NLSP_PLANAR_SOLVER_RECOVERY"]="PARTIAL"
    return result


def plot_recovery(bundle):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bundle=Path(bundle);data=np.load(bundle/"p24_p32_arrays.npz");summary=_read_json(bundle/"summary.json")
    T1=summary["compatibility"]["T1"];time_values=data["time"]/T1
    fig,axes=plt.subplots(2,1,figsize=(7,5),layout="constrained",sharex=True)
    for axis,part in zip(axes,("q","velocity")):
        axis.plot(time_values,data[part+"_difference_L2"][:,3],label="total",lw=.8)
        axis.plot(time_values,data[part+"_tail_L2"][:,3],label="outside p24",lw=.7)
        axis.plot(time_values,data[part+"_common_L2"][:,3],label="common space",lw=.7)
        axis.set_ylabel(r"$\|\Delta c\|_{L^2}$" if part=="q" else r"$\|\Delta c_t\|_{L^2}$")
        axis.grid(alpha=.2);axis.legend(fontsize=8)
    axes[-1].set_xlabel(r"$t/T_1$")
    folder=bundle/"figures";folder.mkdir(exist_ok=True)
    fig.savefig(folder/"contraction_projection.pdf");fig.savefig(folder/"contraction_projection.png",dpi=200);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(8,3.4),layout="constrained")
    for axis,key in zip(axes,("q_c","velocity_c")):
        times=data[key+"_times"];index=int(np.argmax(np.max(abs(data[key+"_difference"]),axis=1)))
        axis.plot(data[key+"_x"],data[key+"_difference"][index],lw=.8)
        axis.axvline(.1,color="grey",lw=.6);axis.axvline(.9,color="grey",lw=.6)
        axis.set_xlabel(r"$s/L$");axis.set_ylabel(r"$\Delta c$" if key=="q_c" else r"$\Delta c_t$")
        axis.set_title(f"t/T1={times[index]/T1:.3f}",fontsize=9);axis.grid(alpha=.2)
    fig.savefig(folder/"contraction_localization.pdf");fig.savefig(folder/"contraction_localization.png",dpi=200);plt.close(fig);data.close()
    return [str(folder/name) for name in ("contraction_projection.pdf","contraction_localization.pdf")]


def main():
    parser=argparse.ArgumentParser(description=__doc__);actions=parser.add_mutually_exclusive_group(required=True)
    actions.add_argument("--diagnose",action="store_true");actions.add_argument("--compute",action="store_true")
    actions.add_argument("--report-only",type=Path);actions.add_argument("--plot-only",type=Path)
    parser.add_argument("--historical-bundle",type=Path,default=SOURCE);parser.add_argument("--baseline-helper",type=Path)
    parser.add_argument("--output-dir",type=Path,default=OUTPUT);args=parser.parse_args()
    if args.report_only or args.plot_only:
        bundle=args.report_only or args.plot_only;summary=validate_recovery(bundle)
        figures=plot_recovery(bundle) if args.plot_only else []
        print(json.dumps({"bundle":str(bundle),"statuses":summary["statuses"],"figures":figures,"integrations":0,"profiling_seconds":0,"symbolic_derivations":0,"root_solves":0}));return
    key,item=recovery_identity(args.historical_bundle,args.baseline_helper,"compute" if args.compute else "diagnose")
    bundle=args.output_dir/key
    if (bundle/"manifest.json").is_file():
        summary=validate_recovery(bundle,item);print(json.dumps({"bundle":str(bundle),"statuses":summary["statuses"],"cache":"reused","integrations":0,"profiling_seconds":0,"symbolic_derivations":0,"root_solves":0}));return
    bundle.mkdir(parents=True,exist_ok=True)
    summary=run_recovery(args.historical_bundle,bundle,args.baseline_helper,item,args.compute)
    write_json(bundle/"summary.json",summary)
    hashes={str(p.relative_to(bundle)):sha(p) for p in bundle.rglob("*") if p.is_file() and not p.name.endswith(".tmp")}
    write_json(bundle/"manifest.json",{"identity":item,"artifact_hashes":hashes,"git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),"command":sys.argv})
    write_json(args.output_dir/"current.json",{"bundle":str(bundle),"fingerprint":key})
    print(json.dumps({"bundle":str(bundle),"statuses":summary["statuses"],"integrations":summary["integrations"],"budget_used_seconds":summary["budget_used_seconds"]}));return


if __name__=="__main__":main()
