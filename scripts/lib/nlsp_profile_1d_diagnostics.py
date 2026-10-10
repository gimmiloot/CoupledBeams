"""Saved-state contraction and bending-profile diagnostics, with no solver calls.

The full trigonometric strain is a postprocessing measure from the accepted
kinematics, not a replacement for the retained cubic strain or quartic dynamics.
All coordinates are the saved independently evolving Shen coefficients.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial import Legendre
from numpy.polynomial.legendre import legder, legmul, legval

FIELDS = ("u", "w", "theta", "c")
TIME_FRACTIONS = (0., .25, .5, .75, 1.)
BOUNDARY_LENGTHS = 3.
VERSION = "nlsp-saved-profile-one-d-v1"


def _sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def exact_saved_indices(actual_times, requested_times):
    """Accept only exact saved timestamps; no invented or interpolated time."""
    actual, requested = np.asarray(actual_times), np.asarray(requested_times)
    if actual.ndim != 1 or np.any(np.diff(actual) <= 0) or not np.isfinite(actual).all():
        raise ValueError("Saved physical timestamps must be finite and increasing")
    indices = np.searchsorted(actual, requested)
    if np.any(indices >= len(actual)) or not np.array_equal(actual[indices], requested):
        raise ValueError("NOT_AVAILABLE: a requested physical time was not saved exactly")
    return indices


def planar_measures(fields, gradients, coefficients):
    """Canonical u,w,theta,c ordering and exact/project-retained strain measures."""
    fields, gradients = np.asarray(fields), np.asarray(gradients)
    if fields.shape != gradients.shape or fields.shape[-1] != 4:
        raise ValueError("Expected matching four-field values and gradients")
    theta, c = fields[..., 2], fields[..., 3]
    us, ws = gradients[..., 0], gradients[..., 1]
    gamma = (1. + us) * np.cos(theta) + ws * np.sin(theta) - 1.
    shear = -(1. + us) * np.sin(theta) + ws * np.cos(theta)
    retained = us + theta * ws - theta**2 / 2. - us * theta**2 / 2.
    retained_shear = ws - theta - theta * us - ws * theta**2 / 2. + theta**3 / 6.
    C, nu = coefficients.C, coefficients.nu
    N, Q = C * (gamma + nu * c), coefficients.S * shear
    return {"Gamma1": gamma, "Gamma2": shear,
        "Gamma1_retained": retained, "Gamma2_retained": retained_shear,
        "N": N,
        "N_retained": C * (retained + nu * c),
        "Q": Q, "F_axial_global": N * np.cos(theta) - Q * np.sin(theta),
        "c_plus_nu_Gamma1": c + nu * gamma,
        "minus_nu_Gamma1": -nu * gamma, "theta_minus_w_s": theta - ws}


def shen_to_legendre(raw_field):
    """Exact algebraic P_n-P_(n+2) conversion, no tail deletion or filtering."""
    raw = np.asarray(raw_field, dtype=float)
    if raw.ndim != 1 or not np.isfinite(raw).all():
        raise ValueError("Invalid physical Shen coefficients")
    result = np.zeros(len(raw) + 2)
    result[:len(raw)] += raw
    result[2:] -= raw
    return result


def polynomial_l2(coefficients, length):
    degrees = np.arange(len(coefficients))
    return float(np.sqrt(np.sum(length * np.asarray(coefficients)**2 / (2 * degrees + 1))))


def legendre_tail_diagnostic(coefficients, length, x):
    """Report an explicit upper-degree contribution; never change the profile."""
    coefficients = np.asarray(coefficients)
    p = len(coefficients) - 1
    cutoff = int(np.ceil(.75 * p))
    tail = np.zeros_like(coefficients)
    tail[cutoff:] = coefficients[cutoff:]
    full_l2 = polynomial_l2(coefficients, length)
    tail_l2 = polynomial_l2(tail, length)
    full_derivative = legder(coefficients) * (2 / length)
    tail_derivative = legder(tail) * (2 / length)
    full_d_l2 = polynomial_l2(full_derivative, length)
    tail_d_l2 = polynomial_l2(tail_derivative, length)
    xi = 2 * np.asarray(x) / length - 1
    return {"degree_cutoff": cutoff, "criterion": "Legendre degrees >= ceil(3*p/4)",
        "full_L2": full_l2, "tail_L2": tail_l2,
        "tail_L2_fraction": tail_l2 / full_l2 if full_l2 else None,
        "full_derivative_L2": full_d_l2, "tail_derivative_L2": tail_d_l2,
        "tail_derivative_L2_fraction": tail_d_l2 / full_d_l2 if full_d_l2 else None,
        "tail_max_sampled": float(np.max(abs(legval(xi, tail)))),
        "tail_derivative_max_sampled": float(np.max(abs(legval(xi, tail_derivative)))),
        "tail_profile": legval(xi, tail), "tail_derivative": legval(xi, tail_derivative),
        "qualification": "Unfiltered spectral-content diagnostic; no new convergence gate"}


def classical_axial_benchmark(w_legendre, length, x):
    """Integrate bar(epsilon)-w_s^2/2 exactly as a polynomial diagnostic."""
    slope = legder(np.asarray(w_legendre)) * (2 / length)
    slope2 = legmul(slope, slope)
    slope2_polynomial = Legendre(slope2)
    antiderivative = slope2_polynomial.integ()
    integral = length / 2 * (antiderivative(1.) - antiderivative(-1.))
    average = float(integral / (2 * length))
    us_coefficients = -.5 * slope2
    us_coefficients[0] += average
    u_polynomial = Legendre(us_coefficients).integ() * (length / 2)
    u_polynomial = u_polynomial - u_polynomial(-1.)
    xi = 2 * np.asarray(x) / length - 1
    return {"mean_bending_extension": average,
        "u_s": legval(xi, us_coefficients), "u": u_polynomial(xi),
        "u_endpoint_error": float(max(abs(u_polynomial(-1.)), abs(u_polynomial(1.)))),
        "qualification": "Classical no-shear quasi-static explanatory benchmark, not the dynamic M-H equation"}


def zero_crossings(x, values):
    """Locate sign changes by linear postprocessing of samples, no root solver."""
    x, values = np.asarray(x), np.asarray(values)
    crossing = np.flatnonzero(values[:-1] * values[1:] < 0)
    points = [float(x[i] - values[i] * (x[i+1] - x[i]) / (values[i+1] - values[i])) for i in crossing]
    zero = np.flatnonzero(values == 0)
    points.extend(float(x[i]) for i in zero if 0 < i < len(x)-1)
    return sorted(set(points))


def symmetry_metrics(x, fields):
    x, fields = np.asarray(x), np.asarray(fields)
    if not np.allclose(x + x[::-1], x[0] + x[-1], rtol=0., atol=1e-14):
        raise ValueError("Symmetry sampling coordinates must be reflected pairs")
    parity = np.array((-1., 1., -1., 1.))
    errors = fields - fields[::-1] * parity
    return {f: {"parity": "even" if parity[i] == 1 else "odd",
        "max_absolute_error": float(np.max(abs(errors[:, i]))),
        "full_profile_scale": float(np.max(abs(fields[:, i]))),
        "relative_error": float(np.max(abs(errors[:, i])) / np.max(abs(fields[:, i])))
            if np.max(abs(fields[:, i])) else None,
        "endpoint_max": float(np.max(abs(fields[[0, -1], i])))} for i, f in enumerate(FIELDS)}


def _norm(value, weights):
    return {"max": float(np.max(abs(value))), "L2": float(np.sqrt(weights @ (value * value)))}


def _regional_measure(value, x, weights, ell):
    inner = (x >= BOUNDARY_LENGTHS * ell) & (x <= x[-1] - BOUNDARY_LENGTHS * ell)
    result = {}
    for name, mask in (("full", np.ones(len(x), dtype=bool)), ("interior", inner), ("boundary", ~inner)):
        if not np.any(mask):
            result[name] = {"status": "NOT_AVAILABLE"}
            continue
        ww, vv = weights[mask], value[mask]
        result[name] = {"mean": float(ww @ vv / ww.sum()), "min": float(vv.min()),
            "max": float(vv.max()), "max_abs": float(np.max(abs(vv))),
            "L2": float(np.sqrt(ww @ (vv * vv))),
            "sample_coordinate_range": [float(x[mask][0]), float(x[mask][-1])]}
    return result


def _grid_weights(x):
    """Trapezoidal physical weights for plot/region diagnostics, not old gates."""
    dx = np.diff(x)
    return np.r_[dx[0] / 2, (dx[:-1] + dx[1:]) / 2, dx[-1] / 2]


def _load_discretizations(source):
    """Restore frozen serialized action and mass scaling, without derivation/eigen."""
    from scripts.lib import nlsp_fem3a_1d_reference as one
    from scripts.lib import weakly_nonlinear_planar_dynamics as dynamics
    root = Path(__file__).resolve().parents[2]
    item = json.loads((source / "provenance.json").read_text(encoding="utf8"))
    science = item["config"]
    reference = one.load_reference(root, science["source_static"]["bundle"],
        science["source_fem1"]["bundle"], science["source_action"]["bundle"])
    disc64 = reference["disc"]
    disc48 = dynamics.PlanarGalerkin(disc64.coefficients, 48, length=disc64.length,
        nq=97, model=disc64.model, whiten=True)
    return {64: disc64, 48: disc48}, reference


def audit_saved_one_d(source_bundle, output_bundle, *, spatial_points=2001):
    """Write Stage E/F evidence from five exact saved coefficient rows only."""
    from scripts.analysis.verify_nlsp_nonlinear_static_3d_fem import fem2_tim_uniform
    source, output = Path(source_bundle), Path(output_bundle)
    if output.resolve().is_relative_to(source.resolve()):
        raise ValueError("Diagnostic output must not overwrite the source bundle")
    output.mkdir(parents=True, exist_ok=True)
    discs, reference = _load_discretizations(source)
    assembly_counters = {str(p): disc.counters() for p, disc in discs.items()}
    for counters in assembly_counters.values():
        if any(counters[name] != 0 for name in
               ("rhs_calls", "jacobian_calls", "mass_factorizations", "linear_eigendecompositions")):
            raise ValueError("Source reconstruction unexpectedly called a scientific solver path")
    coeff = discs[64].coefficients
    L, T = discs[64].length, reference["T1"]
    times = np.array(TIME_FRACTIONS) * T
    x = np.linspace(0., L, spatial_points)
    weights = _grid_weights(x)
    ell = float(np.sqrt(coeff.H / coeff.C))
    nodes, gauss_weights = np.polynomial.legendre.leggauss(257)
    gx, gw = (nodes + 1) * L / 2, gauss_weights * L / 2
    arrays = {"x": x, "times": times, "time_fractions": np.array(TIME_FRACTIONS)}
    historical = json.loads((source / "one_d_all8_spatial.json").read_text(encoding="utf8"))
    summary = {"version": VERSION, "source_bundle": str(source),
        "source_manifest_sha256": _sha(source / "manifest.json"),
        "source_hashes": {}, "physical_coefficients": coeff.values(), "ell_c": ell,
        "regional_policy": {"boundary_width": BOUNDARY_LENGTHS * ell,
            "definition": "interior=[3*ell_c,L-3*ell_c]; fixed before profile metrics",
            "no_new_boundary_conditions": True},
        "spatial_policy": {"sample_points": spatial_points, "sensitivity_gauss_points": 257,
            "norms": "physical max and L2; no independently rescaled profiles",
            "five_time_metrics_not_replacement_for_full_horizon_gate": True},
        "historical_full_horizon_c_gate": historical["fields"]["q_c"],
        "historical_full_horizon_c_velocity_gate": historical["fields"]["velocity_c"],
        "contraction_certification_status": "PARTIAL",
        "restored_discretization_assembly_counters": assembly_counters,
        "cases": {}, "five_time_spatial_sensitivity": {},
        "new_scientific_calls": {"CCX": 0, "Gmsh": 0, "ODE": 0,
            "static_Newton": 0, "eigenanalysis": 0, "root_searches": 0, "BVP": 0,
            "symbolic_derivations": 0}}
    states = {}
    for p in (48, 64):
        disc = discs[p]
        path = source / f"one_d_p{p}_nonlinear.npz"
        summary["source_hashes"][path.name] = _sha(path)
        with np.load(path, allow_pickle=False) as saved:
            actual = saved["times"]
            indices = exact_saved_indices(actual, times)
            q, velocity = saved["q"][indices], saved["velocity"][indices]
        if q.shape != (5, disc.ndof) or not np.isfinite(q).all() or not np.isfinite(velocity).all():
            raise ValueError("Saved state dimension or finite-value check failed")
        fields = disc.reconstruct_series(q, x)
        gradients = disc.reconstruct_series(q, x, derivative=1)
        derivatives2 = disc.reconstruct_series(q, x, derivative=2)
        speeds = disc.reconstruct_series(velocity, x)
        measures = planar_measures(fields, gradients, coeff)
        gamma_s = (derivatives2[..., 0] * np.cos(fields[..., 2])
            + derivatives2[..., 1] * np.sin(fields[..., 2])
            + gradients[..., 2] * measures["Gamma2"])
        measures["Gamma1_s"] = gamma_s
        measures["N_s"] = coeff.C * (gamma_s + coeff.nu * gradients[..., 3])
        arrays.update({f"p{p}_fields": fields, f"p{p}_gradients": gradients,
            f"p{p}_derivatives2": derivatives2, f"p{p}_velocities": speeds,
            f"p{p}_q": q, f"p{p}_velocity_coordinates": velocity})
        arrays.update({f"p{p}_{name}": value for name, value in measures.items()})
        states[p] = (q, velocity)
        case = {"saved_indices": indices.tolist(), "actual_times": actual[indices].tolist(),
            "no_time_interpolation": True, "p": p, "nq": disc.nq, "ndof": disc.ndof,
            "length_over_p": L / p, "quadrature_first_endpoint_distance": float(disc.x[0]),
            "quadrature_nodes_in_left_3ell_region": int(np.sum(disc.x <= 3 * ell)),
            "global_polynomial_space_not_uniform_mesh": True, "snapshots": []}
        tail_fields, tail_gradients, classical_fields, classical_gradients = [], [], [], []
        for index, tau in enumerate(TIME_FRACTIONS):
            raw = disc.raw_coefficients(q[index])
            c_leg = shen_to_legendre(raw[disc.slices["c"]])
            w_leg = shen_to_legendre(raw[disc.slices["w"]])
            tail = legendre_tail_diagnostic(c_leg, L, x)
            classic = classical_axial_benchmark(w_leg, L, x)
            tail_fields.append(tail.pop("tail_profile")); tail_gradients.append(tail.pop("tail_derivative"))
            classical_fields.append(classic.pop("u")); classical_gradients.append(classic.pop("u_s"))
            benchmark_u, benchmark_us = classical_fields[-1], classical_gradients[-1]
            norms = {name: _norm(measures[name][index], weights) for name in
                ("Gamma1", "Gamma2", "N", "c_plus_nu_Gamma1", "theta_minus_w_s")}
            classic_error = _norm(fields[index, :, 0] - benchmark_u, weights)
            mean_N = float(weights @ measures["N"][index] / L)
            extrema_locations = zero_crossings(x, gradients[index, :, 0])
            inner = (x >= 3 * ell) & (x <= L - 3 * ell)
            regional_relative = {}
            for label, mask in (("full", np.ones(len(x), bool)), ("interior", inner), ("boundary", ~inner)):
                c_l2 = float(np.sqrt(weights[mask] @ fields[index, mask, 3]**2))
                relation_l2 = float(np.sqrt(weights[mask] @ measures["c_plus_nu_Gamma1"][index, mask]**2))
                regional_relative[label] = {"c_L2": c_l2, "c_plus_nu_Gamma1_L2": relation_l2,
                    "ratio_to_regional_c_L2": relation_l2 / c_l2 if c_l2 else None,
                    "no_exact_dynamic_closure_claimed": True}
            case["snapshots"].append({"tau": tau, "time": float(times[index]),
                "symmetry": symmetry_metrics(x, fields[index]), "norms": norms,
                "c_mean": float(weights @ fields[index, :, 3] / L),
                "Gamma1_mean": float(weights @ measures["Gamma1"][index] / L),
                "contraction_relation": _regional_measure(measures["c_plus_nu_Gamma1"][index], x, weights, ell),
                "contraction_relation_regional_L2": regional_relative,
                "c_regions": _regional_measure(fields[index, :, 3], x, weights, ell),
                "minus_nu_Gamma1_regions": _regional_measure(measures["minus_nu_Gamma1"][index], x, weights, ell),
                "N_mean": mean_N, "N_min": float(measures["N"][index].min()),
                "N_max": float(measures["N"][index].max()),
                "N_nonconstant_max_deviation_from_mean": float(np.max(abs(measures["N"][index] - mean_N))),
                "N_s_max_abs": float(np.max(abs(measures["N_s"][index]))),
                "full_kinematic_global_axial_flux_min": float(measures["F_axial_global"][index].min()),
                "full_kinematic_global_axial_flux_max": float(measures["F_axial_global"][index].max()),
                "u_extremum_locations_sampled": extrema_locations,
                "classical_u_extremum_locations_sampled": zero_crossings(x, benchmark_us),
                "u_zero_crossings": zero_crossings(x, fields[index, :, 0]),
                "u_s_zero_crossings": extrema_locations,
                "Gamma1_zero_crossings": zero_crossings(x, measures["Gamma1"][index]),
                "N_zero_crossings": zero_crossings(x, measures["N"][index]),
                "u_classical_max_difference": classic_error["max"],
                "u_classical_L2_difference": classic_error["L2"],
                "u_classical_max": float(np.max(abs(benchmark_u))),
                "u_actual_max": float(np.max(abs(fields[index, :, 0]))),
                "classical": classic, "c_spectral_tail": tail,
                "full_vs_retained_Gamma1_max": float(np.max(abs(measures["Gamma1"][index] - measures["Gamma1_retained"][index]))),
                "full_vs_retained_N_max": float(np.max(abs(measures["N"][index] - measures["N_retained"][index])))})
        arrays.update({f"p{p}_c_tail": np.array(tail_fields), f"p{p}_c_tail_gradient": np.array(tail_gradients),
            f"p{p}_classical_u": np.array(classical_fields), f"p{p}_classical_u_s": np.array(classical_gradients)})
        summary["cases"][str(p)] = case
    for part, derivative in (("fields", 0), ("gradients", 1), ("velocities", 0)):
        coords = {p: states[p][1 if part == "velocities" else 0] for p in (48, 64)}
        evaluated = {p: discs[p].reconstruct_series(coords[p], gx, derivative=derivative) for p in (48, 64)}
        for field_index, field in enumerate(FIELDS):
            diff = evaluated[48][..., field_index] - evaluated[64][..., field_index]
            abs_l2 = np.sqrt(np.einsum("tx,x,tx->t", diff, gw, diff))
            abs_max = np.max(abs(diff), axis=1)
            ref_l2 = np.sqrt(np.einsum("tx,x,tx->t", evaluated[64][..., field_index], gw, evaluated[64][..., field_index]))
            ref_max = np.max(abs(evaluated[64][..., field_index]), axis=1)
            summary["five_time_spatial_sensitivity"][part + "_" + field] = {
                "absolute_L2": abs_l2.tolist(), "absolute_max": abs_max.tolist(),
                "full_five_time_reference_L2": float(ref_l2.max()),
                "full_five_time_reference_max": float(ref_max.max()),
                "relative_to_fixed_five_time_max": (abs_max / ref_max.max()).tolist() if ref_max.max() else None,
                "relative_to_fixed_five_time_L2": (abs_l2 / ref_l2.max()).tolist() if ref_l2.max() else None,
                "no_new_threshold": True}
    linear_path = source / "one_d_p64_linear.npz"
    summary["source_hashes"][linear_path.name] = _sha(linear_path)
    with np.load(linear_path, allow_pickle=False) as saved:
        linear_q0 = saved["q"][0]
    analytic, analytic_s = fem2_tim_uniform(x, reference["preflight"]["load"]["q"], L, coeff.Bp, coeff.S)
    linear_fields = discs[64].reconstruct(linear_q0, x)
    linear_gradients = discs[64].reconstruct(linear_q0, x, derivative=1)
    arrays.update(analytic_initial_linear_fields=analytic, analytic_initial_linear_gradients=analytic_s,
        saved_initial_linear_fields=linear_fields, saved_initial_linear_gradients=linear_gradients)
    summary["initial_linear_Timoshenko_reference"] = {
        "saved_linear_vs_analytic_field_max": {f: float(np.max(abs(linear_fields[:, i] - analytic[:, i]))) for i, f in enumerate(FIELDS)},
        "saved_linear_vs_analytic_gradient_max": {f: float(np.max(abs(linear_gradients[:, i] - analytic_s[:, i]))) for i, f in enumerate(FIELDS)},
        "nonlinear_initial_vs_analytic_field_max": {f: float(np.max(abs(arrays["p64_fields"][0, :, i] - analytic[:, i]))) for i, f in enumerate(FIELDS)},
        "theta_axis_contract": "positive 1D w is negative global Y; positive theta is rotation about negative global Z",
        "no_new_static_solution": True}
    increments = {str(p): {name: value - assembly_counters[str(p)][name]
        for name, value in disc.counters().items()} for p, disc in discs.items()}
    if any(value != 0 for current in increments.values() for value in current.values()):
        raise ValueError("Postprocessing unexpectedly changed scientific model counters")
    summary["postprocessing_model_counter_increments"] = increments
    np.savez_compressed(output / "one_d_profile_audit.npz", **arrays)
    (output / "one_d_profile_audit.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf8")
    with (output / "one_d_profile_audit.csv").open("w", newline="", encoding="utf8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("p", "tau", "s", "u", "w", "theta", "c", "u_s", "w_s", "theta_s", "c_s", "Gamma1", "N", "c_plus_nu_Gamma1", "classical_u"))
        for p in (48, 64):
            for index, tau in enumerate(TIME_FRACTIONS):
                for j, point in enumerate(x):
                    writer.writerow((p, tau, point, *arrays[f"p{p}_fields"][index, j],
                        *arrays[f"p{p}_gradients"][index, j], arrays[f"p{p}_Gamma1"][index, j],
                        arrays[f"p{p}_N"][index, j], arrays[f"p{p}_c_plus_nu_Gamma1"][index, j],
                        arrays[f"p{p}_classical_u"][index, j]))
    return summary


def audit_actual_one_d_states(source_bundle, output_bundle, selected_states, *, spatial_points=801):
    """Read accepted dense polynomials at the selected FEM physical times.

    No evolution equation is evaluated. Each distinct saved physical time is
    evaluated once, and repeated mesh/state labels retain their own provenance.
    The original five-time exact-row audit and its files are never overwritten.
    """
    import time
    from scripts.lib import nlsp_fem3b_continuation as saved_dense
    source, output = Path(source_bundle), Path(output_bundle)
    if output.resolve().is_relative_to(source.resolve()):
        raise ValueError("Diagnostic output must not overwrite the source bundle")
    states = [dict(state) for state in selected_states]
    if not states or len({state["name"] for state in states}) != len(states):
        raise ValueError("Selected physical state names must be nonempty and unique")
    times = np.array([state["actual_time"] for state in states], dtype=float)
    if not np.isfinite(times).all() or np.any(times < 0):
        raise ValueError("Selected physical times must be finite and nonnegative")
    unique_times, inverse = np.unique(times, return_inverse=True)
    manifest_path = source / "manifest.json"
    manifest_sha = _sha(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf8"))
    expected = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
    source_hashes = {}
    for name in ("dense_p64.npz", "one_d_p64_nonlinear.npz"):
        actual_sha = _sha(source / name)
        if expected.get(name) != actual_sha:
            raise ValueError("Saved actual-time source absent/corrupted: " + name)
        source_hashes[name] = actual_sha
    output.mkdir(parents=True, exist_ok=True)
    code_root = output / "execution_code"
    code_root.mkdir(exist_ok=True)
    helper_source = Path(__file__).read_bytes()
    helper_sha = hashlib.sha256(helper_source).hexdigest()
    snapshot = code_root / "one_d_actual_stage.py"
    if snapshot.exists():
        if snapshot.read_bytes() != helper_source:
            raise ValueError("A different actual-time execution snapshot already exists")
    else:
        snapshot.write_bytes(helper_source)
    started = time.perf_counter()
    discs, reference = _load_discretizations(source)
    disc = discs[64]
    before = disc.counters()
    if any(before[name] for name in
           ("rhs_calls", "jacobian_calls", "mass_factorizations", "linear_eigendecompositions")):
        raise ValueError("Saved actual-time reconstruction invoked a scientific solver")
    qv = saved_dense.evaluate_dense_records(source / "dense_p64.npz", unique_times)
    if qv.shape != (len(unique_times), 2 * disc.ndof) or not np.isfinite(qv).all():
        raise ValueError("Dense coefficient dimensions or finite-value check failed")
    with np.load(source / "one_d_p64_nonlinear.npz", allow_pickle=False) as saved:
        if saved["times"][0] != 0:
            raise ValueError("Saved initial coefficient history does not begin at zero")
        q0, v0 = saved["q"][0], saved["velocity"][0]
    if np.any(v0 != 0):
        raise ValueError("Saved physical initial velocities are not zero")
    if len(unique_times) and unique_times[0] == 0:
        qv[0, :disc.ndof], qv[0, disc.ndof:] = q0, v0
    x = np.linspace(0., disc.length, spatial_points)
    q, velocity = qv[inverse, :disc.ndof], qv[inverse, disc.ndof:]
    fields = disc.reconstruct_series(q, x)
    gradients = disc.reconstruct_series(q, x, derivative=1)
    measures = planar_measures(fields, gradients, disc.coefficients)
    arrays = {"names": np.array([state["name"] for state in states]), "times": times,
        "requested_times": np.array([state["requested_time"] for state in states]),
        "requested_tau": np.array([state["requested_tau"] for state in states]),
        "actual_tau": times / reference["T1"], "x": x, "fields": fields,
        "gradients": gradients, "q": q, "velocity_coordinates": velocity,
        "unique_times": unique_times, "unique_q": qv[:, :disc.ndof],
        "unique_velocity": qv[:, disc.ndof:], "state_unique_index": inverse,
        **measures}
    increments = {name: value - before[name] for name, value in disc.counters().items()}
    if any(increments.values()):
        raise ValueError("Actual-time postprocessing changed scientific model counters")
    if _sha(manifest_path) != manifest_sha or any(_sha(source / name) != value for name, value in source_hashes.items()):
        raise ValueError("Source changed during actual-time postprocessing")
    summary = {"version": VERSION, "phase": "actual-FEM-times-saved-dense-polynomials",
        "source_bundle": str(source), "source_manifest_sha256": manifest_sha,
        "source_hashes": source_hashes, "source_hashes_unchanged_after_processing": True,
        "helper_sha256": helper_sha, "execution_code_snapshot": str(snapshot.relative_to(output)),
        "selected_states": states, "states": len(states), "distinct_physical_times": len(unique_times),
        "unique_times": unique_times.tolist(), "state_unique_index": inverse.tolist(),
        "field_order": list(FIELDS), "p": 64, "nq": disc.nq, "ndof": disc.ndof,
        "T1": reference["T1"], "spatial_points": spatial_points,
        "restored_discretization_assembly_counters": before,
        "postprocessing_model_counter_increments": increments,
        "postprocessing_seconds": time.perf_counter() - started,
        "saved_dense_evaluations": 1, "new_scientific_solver_calls": 0,
        "time_policy": "Evaluate retained accepted Radau polynomials at actual FEM times; no new solve, fitting or time integration",
        "initial_state_policy": "exact saved q0/v0 at t=0; STATIC comparison datum, not a native DYNAMIC frame",
        "historical_five_time_audit_overwritten": False}
    np.savez_compressed(output / "actual_time_one_d.npz", **arrays)
    (output / "actual_time_one_d.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf8")
    return summary
