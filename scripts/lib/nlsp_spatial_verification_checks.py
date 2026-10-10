"""Bounded checks of the frozen seven-field action and its discretization.

This module reads saved references and evaluates algebra, matrices and RHS
values. It never integrates an IVP, solves an equilibrium, or starts FEM.
The limited block eigendecompositions compare already identified linear
frequencies; they are not a new continuum root search.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time
from types import SimpleNamespace
import zipfile

import numpy as np

from scripts.lib import weakly_nonlinear_spatial_rod as rod
from scripts.lib.weakly_nonlinear_planar_dynamics import PlanarGalerkin, _CompiledPolynomials

VERSION = "frozen-seven-field-stage-a-checks-v1"
STRICT_RELATIVE = 2e-12
JACOBIAN_DIRECTIONAL_RELATIVE = 2e-7  # frozen planar pilot config gate
LINEAR_FREQUENCY_RELATIVE = 2e-8  # existing seven-field linear reference check
PLANAR_FIELDS = ("u", "w", "theta", "c")
SUFFIXES = ("", "_s", "_t", "_ss", "_st", "_tt")


def _sha(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def source_evidence(paths):
    """Check each used source against its own nearest immutable manifest."""
    evidence = {}
    for name, supplied in paths.items():
        path = Path(supplied).resolve()
        if not path.is_file():
            raise FileNotFoundError("Required saved Stage A source is absent: " + str(path))
        parent = next((p for p in path.parents if (p / "manifest.json").is_file()), None)
        if parent is None:
            raise ValueError("Saved Stage A source has no parent manifest: " + str(path))
        manifest = json.loads((parent / "manifest.json").read_text(encoding="utf8"))
        artifacts = manifest.get("artifact_hashes", manifest.get("artifacts", {}))
        relative = path.relative_to(parent).as_posix()
        expected = artifacts.get(relative, artifacts.get(relative.replace("/", "\\")))
        actual = _sha(path)
        if isinstance(expected, dict):
            expected = expected.get("sha256")
        if actual != expected:
            raise ValueError("Saved Stage A source hash mismatch: " + str(path))
        evidence[name] = {"path": str(path), "sha256": actual,
                          "manifest": str(parent / "manifest.json"),
                          "manifest_sha256": _sha(parent / "manifest.json")}
    return evidence


def read_npz_rows(path, name, indices):
    """Stream selected C-order rows without retaining the complete history."""
    indices = np.asarray(indices, dtype=int)
    if indices.ndim != 1 or len(set(indices.tolist())) != len(indices):
        raise ValueError("Selected saved-history indices must be unique")
    with zipfile.ZipFile(path) as archive, archive.open(name + ".npy") as stream:
        version = np.lib.format.read_magic(stream)
        shape, fortran, dtype = np.lib.format._read_array_header(stream, version)
        if fortran or len(shape) != 2 or dtype.hasobject:
            raise ValueError("Saved history must be a two-dimensional C-order numeric array")
        if np.any(indices < 0) or np.any(indices >= shape[0]):
            raise ValueError("Selected saved-history row lies outside actual coverage")
        row_bytes = shape[1] * dtype.itemsize
        result = np.empty((len(indices), shape[1]), dtype=dtype)
        position = 0
        for row in np.sort(indices):
            remaining = (int(row) - position) * row_bytes
            while remaining:
                read = stream.read(min(remaining, 1024 * 1024))
                if not read:
                    raise ValueError("Truncated saved history before requested row")
                remaining -= len(read)
            data = stream.read(row_bytes)
            if len(data) != row_bytes:
                raise ValueError("Incomplete saved-history row")
            result[np.flatnonzero(indices == row)[0]] = np.frombuffer(data, dtype=dtype)
            position = int(row) + 1
        return result


def difference(actual, reference, threshold=STRICT_RELATIVE):
    """Absolute and relative differences, with explicit unmodified scales."""
    actual, reference = np.asarray(actual), np.asarray(reference)
    if actual.shape != reference.shape or not np.isfinite(actual).all() or not np.isfinite(reference).all():
        return {"status": "FAIL", "reason": "shape mismatch or nonfinite values"}
    delta = actual - reference
    norm = float(np.linalg.norm(reference.ravel()))
    maximum = float(np.max(abs(reference), initial=0))
    absolute = float(np.max(abs(delta), initial=0))
    norm_delta = float(np.linalg.norm(delta.ravel()))
    relative_l2 = norm_delta / norm if norm else (0. if not norm_delta else None)
    relative_max = absolute / maximum if maximum else (0. if not absolute else None)
    passed = all(v is not None and v <= threshold for v in (relative_l2, relative_max))
    return {"status": "PASS" if passed else "FAIL", "absolute_max": absolute,
            "absolute_l2": norm_delta, "relative_max": relative_max,
            "relative_l2": relative_l2, "reference_max": maximum,
            "reference_l2": norm, "relative_threshold": threshold,
            "zero_reference_policy": "require exact zero difference"}


def exact_action_checks(model):
    """Independent exact restrictions of the loaded serialized polynomials."""
    s = {name: rod.Polynomial.symbol(name) for name in rod.SYMBOL_ORDER}
    zero = {field + suffix: 0 for field in ("v", "Phi", "psi") for suffix in SUFFIXES}
    theta, c = s["theta"], s["c"]
    cosine, sine = 1 - theta**2 / 2 + theta**4 / 24, theta - theta**3 / 6
    gamma1 = ((1 + s["u_s"]) * cosine + s["w_s"] * sine - 1).truncate(3)
    gamma2 = (-(1 + s["u_s"]) * sine + s["w_s"] * cosine).truncate(3)
    planar_t = (s["m"] * (s["u_t"]**2 + s["w_t"]**2) + s["jp"] * s["c_t"]**2 +
                s["jp"] * (1 + c)**2 * s["theta_t"]**2) / 2
    planar_v = (s["C"] * (gamma1**2 + c**2) + s["S"] * gamma2**2 +
                s["H"] * s["c_s"]**2 + s["Bp"] * s["theta_s"]**2) / 2 + s["nu"] * s["C"] * c * gamma1
    checks = {"planar_T": model.T4.substitute(zero) - planar_t,
              "planar_V": model.V4.substitute(zero) - planar_v}
    lagrangian = model.T4 - model.V4
    for field, residual in zip(rod.FIELD_ORDER, model.residual_a):
        varied = (lagrangian.derivative(field + "_t").total_derivative("t") +
                  lagrangian.derivative(field + "_s").total_derivative("s") - lagrangian.derivative(field))
        checks["frozen_variational_residual_" + field] = varied - residual
    for i in (2, 3, 4):
        checks["planar_inactive_" + rod.FIELD_ORDER[i]] = model.residual_a[i].substitute(zero)
    for label, signs, reversed_s in (
            ("outplane_reflection", (1, 1, -1, -1, -1, 1, 1), False),
            ("second_transverse_reflection", (1, -1, 1, -1, 1, -1, 1), False),
            ("midpoint_reflection", (-1, 1, 1, 1, -1, -1, 1), True)):
        reflection = {f + suffix: sign * (-1 if reversed_s and suffix in ("_s", "_st") else 1) * s[f + suffix]
                      for f, sign in zip(rod.FIELD_ORDER, signs) for suffix in SUFFIXES}
        checks[label + "_T"] = model.T4.substitute(reflection) - model.T4
        checks[label + "_V"] = model.V4.substitute(reflection) - model.V4
        for f, sign, residual in zip(rod.FIELD_ORDER, signs, model.residual_a):
            checks[label + "_" + f] = residual.substitute(reflection) - sign * residual
    fluxes = [model.V4.derivative(f + "_s") for f in rod.FIELD_ORDER]
    power = sum((s[f + "_t"] * flux for f, flux in zip(rod.FIELD_ORDER, fluxes)), rod.Polynomial())
    work = sum((s[f + "_t"] * residual for f, residual in zip(rod.FIELD_ORDER, model.residual_a)), rod.Polynomial())
    checks["energy_identity"] = (model.T4 + model.V4).total_derivative("t") - power.total_derivative("s") - work
    velocity_ids = {rod.SYMBOL_ORDER.index(f + "_t") for f in rod.FIELD_ORDER}
    velocity_degrees = sorted({sum(i in velocity_ids for i in key) for key in model.T4.terms})
    rows = {name: {"status": "PASS" if not polynomial else "FAIL",
                   "difference_terms": len(polynomial.terms)} for name, polynomial in checks.items()}
    passed = all(row["status"] == "PASS" for row in rows.values()) and velocity_degrees == [2]
    return {"status": "PASS" if passed else "FAIL", "checks": rows,
            "T_terms": len(model.T4.terms), "V_terms": len(model.V4.terms),
            "kinetic_velocity_degrees": velocity_degrees,
            "new_symbolic_derivations": 0,
            "operation": "exact differentiation/substitution of frozen serialized action"}


def planar_ids(disc):
    return np.concatenate([np.arange(disc.slices[f].start, disc.slices[f].stop) for f in PLANAR_FIELDS])


def embed_planar(disc, values):
    values = np.asarray(values)
    ids = planar_ids(disc)
    if values.shape != ids.shape:
        raise ValueError("Planar embedding requires the same p and coordinate scaling")
    result = np.zeros(disc.ndof)
    result[ids] = values
    return result


def compare_planar_state(disc, planar, q, velocity, *, jacobian=True):
    ids = planar_ids(disc)
    inactive = np.setdiff1d(np.arange(disc.ndof), ids)
    seven_q, seven_v = embed_planar(disc, q), embed_planar(disc, velocity)
    old, new = planar.potential(q, hessian=True), disc.potential(seven_q, hessian=True)
    old_rhs, new_rhs = planar.rhs(0., np.r_[q, velocity]), disc.rhs(0., np.r_[seven_q, seven_v])
    state_ids = np.r_[ids, disc.ndof + ids]
    rows = {
        "potential": difference(new["V"], old["V"]),
        "gradient": difference(new["gradient"][ids], old["gradient"]),
        "hessian": difference(new["hessian"][np.ix_(ids, ids)], old["hessian"]),
        "mass": difference(disc.mass_matrix(seven_q)[np.ix_(ids, ids)], planar.mass_matrix(q)),
        "inertia": difference(disc.inertial_terms(seven_q, seven_v)[ids], planar.inertial_terms(q, velocity)),
        "rhs": difference(new_rhs[state_ids], old_rhs),
        "inactive_force": difference(new["gradient"][inactive], np.zeros(len(inactive))),
        "inactive_rhs": difference(new_rhs[np.r_[inactive, disc.ndof + inactive]], np.zeros(2 * len(inactive))),
        "energy": difference(disc.energy(seven_q, seven_v), planar.energy(q, velocity)),
        "physical_fields": difference(disc.reconstruct(seven_q)[:, (0, 1, 5, 6)], planar.reconstruct(q)),
    }
    if jacobian:
        new_jacobian, old_jacobian = disc.jacobian(0., np.r_[seven_q, seven_v]), planar.jacobian(0., np.r_[q, velocity])
        rows["jacobian"] = difference(new_jacobian[np.ix_(state_ids, state_ids)], old_jacobian)
        rows["inactive_active_jacobian"] = difference(
            new_jacobian[np.ix_(np.r_[inactive, disc.ndof + inactive], state_ids)],
            np.zeros((2 * len(inactive), len(state_ids))))
    return {"status": "PASS" if all(r["status"] == "PASS" for r in rows.values()) else "FAIL", "checks": rows}


def smooth_coordinates(disc, phase=0):
    """Predeclared low-degree all-field states, not prepared scientific IC."""
    raw = np.zeros((7, disc.n))
    amplitudes = np.array((2e-5, 3e-4, -2e-4, 7e-4, -8e-4, 9e-4, 2e-5))
    count = min(4, disc.n)
    for field in range(7):
        raw[field, :count] = amplitudes[field] * np.cos(np.arange(count) + .3 * field + phase) / (1 + np.arange(count))**2
    q = disc.from_raw_coefficients(raw.ravel())
    raw_v = np.roll(raw, 1, axis=1) * .3
    velocity = disc.from_raw_coefficients(raw_v.ravel())
    acceleration = disc.from_raw_coefficients((-.7 * raw).ravel())
    return q, velocity, acceleration


def weak_flux_projection(disc, coordinate, velocity, acceleration, *, test_matrices=None,
                         test_derivatives=None, endpoint_test_values=None):
    """Independent continuum weak projection before spatial flux differentiation.

    The time/body terms and potential flux are compiled directly from the
    frozen local action. This is a separate stable weak representation; it
    does not replace the historical strong-residual projection or any RHS.
    Nonzero test functions retain the explicit minus boundary-flux work.
    """
    fields = rod.FIELD_ORDER
    names = tuple(f + suffix for suffix in ("", "_s", "_t", "_tt") for f in fields)
    kinetic_body = []
    for field in fields:
        momentum = disc.model.T4.derivative(field + "_t")
        derivative = sum((momentum.derivative(f) * rod.Polynomial.symbol(f + "_t") +
                          momentum.derivative(f + "_t") * rod.Polynomial.symbol(f + "_tt")
                          for f in fields), rod.Polynomial())
        kinetic_body.append(derivative - disc.model.T4.derivative(field))
    body = [term + disc.model.V4.derivative(field) for term, field in zip(kinetic_body, fields)]
    flux = [disc.model.V4.derivative(field + "_s") for field in fields]
    evaluator = _CompiledPolynomials(body + flux, names, disc.coefficients)
    values = np.vstack((disc.reconstruct(coordinate).T,
                        disc.reconstruct(coordinate, derivative=1).T,
                        disc.reconstruct(velocity).T, disc.reconstruct(acceleration).T))
    local = evaluator.evaluate(values)
    if test_matrices is None:
        test_matrices = disc.B
    if test_derivatives is None:
        test_derivatives = disc.D
    endpoint = np.array((0., disc.length))
    if endpoint_test_values is None:
        endpoint_test_values = tuple(disc.basis_at(endpoint).values())
    endpoint_variables = np.vstack((disc.reconstruct(coordinate, points=endpoint).T,
                                   disc.reconstruct(coordinate, points=endpoint, derivative=1).T,
                                   disc.reconstruct(velocity, points=endpoint).T,
                                   disc.reconstruct(acceleration, points=endpoint).T))
    endpoint_flux = evaluator.evaluate(endpoint_variables)[7:]
    volume, boundary = [], []
    for i, (basis, derivative, end_values) in enumerate(zip(test_matrices, test_derivatives, endpoint_test_values)):
        if basis.shape[0] != disc.nq or derivative.shape != basis.shape or end_values.shape != (2, basis.shape[1]):
            raise ValueError("Continuum weak test functions have incompatible volume/endpoint shapes")
        volume.append(basis.T @ (disc.weights * local[i]) + derivative.T @ (disc.weights * local[7 + i]))
        boundary.append(end_values[1] * endpoint_flux[i, 1] - end_values[0] * endpoint_flux[i, 0])
    volume, boundary = np.concatenate(volume), np.concatenate(boundary)
    return {"residual": volume - boundary, "volume_terms": volume,
            "boundary_work": boundary, "boundary_work_max": float(np.max(abs(boundary), initial=0)),
            "endpoint_flux": endpoint_flux, "local_body": local[:7], "local_flux": local[7:],
            "representation": "B*continuum body + D*potential flux - [test*flux]_ends"}


def variational_checks(disc, planar_reference=None):
    q, velocity, acceleration = smooth_coordinates(disc)
    gradient = disc.potential(q)["gradient"]
    assembled = disc.mass_matrix(q) @ acceleration + disc.inertial_terms(q, velocity) + gradient
    projected = disc.weak_residual(q, velocity, acceleration)
    weak = difference(projected, assembled)
    # Restore precisely the historical uncancelled potential-work denominator.
    # Norms of the already summed residual are kept only as extra diagnostics.
    _, local_gradient, _ = disc._local_potential(q)
    uncancelled = float(sum(np.linalg.norm(matrix.T @ (disc.weights * values))
                           for matrix, values in zip(disc._potential_matrices, local_gradient)))
    weak["summed_residual_relative_l2_diagnostic"] = weak.pop("relative_l2")
    weak["summed_residual_relative_max_diagnostic"] = weak.pop("relative_max")
    weak["uncancelled_potential_work_scale"] = uncancelled
    weak["relative_l2"] = weak["absolute_l2"] / max(uncancelled, 1e-30)
    weak["floor"] = 1e-30
    weak["absolute_threshold"] = STRICT_RELATIVE
    weak["normalization"] = "unchanged FEM-2 uncancelled potential-work policy"
    weak["status"] = "PASS" if weak["relative_l2"] <= STRICT_RELATIVE and weak["absolute_max"] <= STRICT_RELATIVE else "PARTIAL"
    # Use exactly the old unsigned-work normalization for analytical energy rate.
    physical_acceleration = disc.acceleration(q, velocity)
    rate = disc.energy_rate(q, velocity, physical_acceleration)
    energy_scale = float(np.sum(abs(velocity * gradient)) +
                         np.sum(abs(velocity * (disc.mass_matrix(q) @ physical_acceleration))))
    energy_relative = abs(rate) / max(energy_scale, 1e-30)
    endpoint = float(np.max(abs(disc.reconstruct(q, points=np.array((0., disc.length))))))
    state = np.r_[q, velocity]
    jacobian = disc.jacobian(0., state)
    # Independent central columns follow the existing planar check, including
    # the same step and norm denominator. No derivative is used inside the RHS.
    columns = [disc.slices[f].start + min(1, disc.n - 1) for f in rod.FIELD_ORDER]
    columns += [disc.ndof + disc.slices[f].start for f in ("Phi", "psi", "theta", "c")]
    directional = []
    step = 1e-7
    for column in columns:
        plus, minus = state.copy(), state.copy()
        plus[column] += step
        minus[column] -= step
        observed = (disc.rhs(0., plus) - disc.rhs(0., minus)) / (2 * step)
        absolute = float(np.linalg.norm(observed - jacobian[:, column]))
        scale = max(float(np.linalg.norm(jacobian[:, column])), 1.)
        relative = absolute / scale
        directional.append({"column": column, "absolute_l2": absolute,
                            "denominator": scale, "relative_l2": relative,
                            "status": "PASS" if relative < JACOBIAN_DIRECTIONAL_RELATIVE else "FAIL"})
    stable_projection = weak_flux_projection(disc, q, velocity, acceleration)
    stable_weak = difference(stable_projection["residual"], assembled)
    stable_weak["summed_residual_relative_l2_diagnostic"] = stable_weak.pop("relative_l2")
    stable_weak["summed_residual_relative_max_diagnostic"] = stable_weak.pop("relative_max")
    stable_weak["uncancelled_potential_work_scale"] = uncancelled
    stable_weak["relative_l2"] = stable_weak["absolute_l2"] / max(uncancelled, 1e-30)
    stable_weak["floor"] = 1e-30
    stable_weak["absolute_threshold"] = STRICT_RELATIVE
    stable_weak["normalization"] = "unchanged FEM-2 uncancelled potential-work policy"
    stable_weak["representation"] = stable_projection["representation"]
    stable_weak["status"] = "PASS" if stable_weak["relative_l2"] <= STRICT_RELATIVE and stable_weak["absolute_max"] <= STRICT_RELATIVE else "FAIL"
    result = {"strong_projection_strict": weak,
            "weak_flux_projection": stable_weak,
            "boundary_work": {"status": "PASS" if endpoint == 0. else "FAIL",
                              "essential_endpoint_max": endpoint,
                              "explicit_boundary_work_max": stable_projection["boundary_work_max"],
                              "reason": "test functions vanish at both ends; no derivative BC"},
            "energy_rate": {"status": "PASS" if energy_relative < STRICT_RELATIVE else "FAIL",
                            "rate": rate, "unsigned_work_scale": energy_scale,
                            "relative": energy_relative, "threshold": STRICT_RELATIVE,
                            "floor": 1e-30, "normalization": "unchanged existing planar unsigned-work policy"},
            "jacobian_directional": {"status": "PASS" if all(r["status"] == "PASS" for r in directional) else "FAIL",
                                     "rows": directional, "central_step": step,
                                     "threshold": JACOBIAN_DIRECTIONAL_RELATIVE},
            "state_policy": "smooth deterministic first four raw Shen coefficients only"}
    if planar_reference is not None:
        ids = planar_ids(disc)
        old_q, old_v, old_a = q[ids], velocity[ids], acceleration[ids]
        old_assembled = (planar_reference.mass_matrix(old_q) @ old_a +
                         planar_reference.inertial_terms(old_q, old_v) +
                         planar_reference.potential(old_q)["gradient"])
        old_delta = planar_reference.weak_residual(old_q, old_v, old_a) - old_assembled
        _, old_local_gradient, _ = planar_reference._local_potential(old_q)
        old_scale = float(sum(np.linalg.norm(matrix.T @ (planar_reference.weights * values))
                              for matrix, values in zip(planar_reference._potential_matrices, old_local_gradient)))
        plane_q, plane_v, plane_a = (embed_planar(disc, values) for values in (old_q, old_v, old_a))
        seven_plane_assembled = (disc.mass_matrix(plane_q) @ plane_a + disc.inertial_terms(plane_q, plane_v) +
                                 disc.potential(plane_q)["gradient"])
        seven_plane_delta = (disc.weak_residual(plane_q, plane_v, plane_a) - seven_plane_assembled)[ids]
        result["existing_planar_float64_control"] = {
            "old_absolute_l2": float(np.linalg.norm(old_delta)),
            "old_absolute_max": float(np.max(abs(old_delta))),
            "old_uncancelled_work_scale": old_scale,
            "old_relative_l2": float(np.linalg.norm(old_delta) / max(old_scale, 1e-30)),
            "same_state_seven_absolute_l2": float(np.linalg.norm(seven_plane_delta)),
            "same_state_seven_relative_l2": float(np.linalg.norm(seven_plane_delta) / max(old_scale, 1e-30)),
            "difference_between_old_and_seven_delta_max": float(np.max(abs(old_delta - seven_plane_delta))),
            "threshold": STRICT_RELATIVE,
            "reason": "existing planar implementation evaluated at identical coordinates; no new integration"}
    return result


def linear_reference_checks(discs, reference):
    blocks = {"inplane_bending": "timoshenko", "outplane_bending": "outplane", "torsion": "torsion", "axial_mh": "mh"}
    rows = []
    for p, disc in sorted(discs.items()):
        for row in reference["merged_spectrum"]:
            omega = float(disc.linear_eigenpairs(blocks[row["family"]])["omega"][row["local_mode"] - 1])
            target = float(row["omega"])
            relative = abs(omega - target) / target
            rows.append({"p": p, "family": row["family"], "local_mode": row["local_mode"],
                         "omega_semidiscrete": omega, "omega_frozen_reference": target,
                         "signed_difference": omega - target, "relative_difference": relative,
                         "status": "PASS" if relative < LINEAR_FREQUENCY_RELATIVE else "FAIL"})
    return {"status": "PASS" if all(row["status"] == "PASS" for row in rows) else "FAIL",
            "rows": rows, "threshold": LINEAR_FREQUENCY_RELATIVE,
            "frequency_units": "radians per normalized time", "continuum_root_searches": 0,
            "scope": "only saved identified modes, four finite Galerkin block eigendecompositions per p"}


def _saved_states(path, fractions):
    with np.load(path, allow_pickle=False) as data:
        times = data["times"].copy()
    indices = np.unique([int(np.argmin(abs(times - fraction * times[-1]))) for fraction in fractions])
    q = read_npz_rows(path, "q", indices)
    velocity = read_npz_rows(path, "velocity", indices)
    return times[indices], indices, q, velocity


def run_stage_a_checks(discs, model, source_paths, output=None):
    """Run bounded checks; no new physical solves or historical mutations.

    source_paths keys: fem1_preflight, short_planar, full_planar_p48,
    full_planar_p64, static_p64 and static_preflight. Paths must be actual
    hash-verified artifacts with their own parent manifests. Optional saved
    histories may be omitted explicitly and are reported NOT_RUN.
    """
    started = time.perf_counter()
    evidence = source_evidence(source_paths)
    if "action_result" in source_paths:
        frozen = json.loads(Path(source_paths["action_result"]).read_text(encoding="utf8"))["polynomials"]
        if (model.T4.serialize() != frozen["T4"] or model.V4.serialize() != frozen["V4"] or
                [residual.serialize() for residual in model.residual_a] != frozen["residuals_A"]):
            raise ValueError("Supplied Stage A model differs from the frozen source action")
    exact = exact_action_checks(model)
    reference = json.loads(Path(source_paths["fem1_preflight"]).read_text(encoding="utf8"))
    if "coefficients" in reference and any(disc.coefficients.values() != reference["coefficients"] for disc in discs.values()):
        raise ValueError("Stage A coefficients differ from the frozen FEM-1 reference")
    if "geometry" in reference and any(disc.length != reference["geometry"]["L"] for disc in discs.values()):
        raise ValueError("Stage A length differs from the frozen FEM-1 reference")
    linear = linear_reference_checks(discs, reference)
    variational, saved = {}, []
    planar_model = SimpleNamespace(T4=model.T4, V4=model.V4, residual_a=model.residual_a,
                                   symbols={name: rod.Polynomial.symbol(name) for name in rod.SYMBOL_ORDER})
    old_discs = {p: PlanarGalerkin(d.coefficients, p, length=d.length, nq=d.nq, model=planar_model, whiten=d.whiten)
                 for p, d in discs.items()}
    for p, disc in sorted(discs.items()):
        variational[str(p)] = variational_checks(disc, old_discs[p])
        for source_name, fractions in (("short_planar", (0., .5, 1.)),
                                       ("full_planar_p" + str(p), (0., .25, .5, .75, 1.))):
            if source_name not in source_paths or (source_name == "short_planar" and p != 64):
                continue
            times, indices, coordinates, velocities = _saved_states(source_paths[source_name], fractions)
            if coordinates.shape[1] != old_discs[p].ndof:
                raise ValueError("Saved planar state differs from the requested p/field order")
            for index, actual_time, q, velocity in zip(indices, times, coordinates, velocities):
                check = compare_planar_state(disc, old_discs[p], q, velocity)
                saved.append({"p": p, "source": source_name, "row": int(index),
                              "actual_time": float(actual_time), **check})
    static = {"status": "NOT_RUN", "reason": "saved static source not supplied"}
    historical_strict = {"status": "PARTIAL", "threshold": STRICT_RELATIVE,
                         "reason": "historical explained float64 strong/weak qualification is not promoted"}
    if "static_p64" in source_paths and "static_preflight" in source_paths and 64 in discs:
        preflight = json.loads(Path(source_paths["static_preflight"]).read_text(encoding="utf8"))
        disc, old = discs[64], old_discs[64]
        load = np.zeros(old.ndof)
        load[old.slices["w"]] = preflight["load"]["q"] * old.B[1].T @ old.weights
        rows = []
        with np.load(source_paths["static_p64"], allow_pickle=False) as data:
            for kind in ("linear", "nonlinear"):
                q = data["q_" + kind].copy()
                q7 = embed_planar(disc, q)
                gradient = disc.K @ q7 if kind == "linear" else disc.potential(q7)["gradient"]
                old_gradient = old.K @ q if kind == "linear" else old.potential(q)["gradient"]
                old_equilibrium = old_gradient - load
                new_equilibrium = gradient[planar_ids(disc)] - load
                rows.append({"case": kind, "gradient_equivalence": difference(gradient[planar_ids(disc)], old_gradient),
                             "new_absolute_residual": float(np.linalg.norm(new_equilibrium)),
                             "old_absolute_residual": float(np.linalg.norm(old_equilibrium)),
                             "new_relative_residual": float(np.linalg.norm(new_equilibrium) / np.linalg.norm(load)),
                             "old_relative_residual": float(np.linalg.norm(old_equilibrium) / np.linalg.norm(load)),
                             "new_equilibrium_solves": 0})
        static = {"status": "PASS" if all(r["gradient_equivalence"]["status"] == "PASS" for r in rows) else "FAIL", "rows": rows}
        historical_strict["saved_evidence"] = [
            {"p": row["p"], "case": kind,
             "relative_difference": row[kind]["strong_action_relative_difference"],
             "absolute_difference": row[kind]["strong_action_absolute_difference"],
             "unsigned_work_scale": row[kind]["strong_action_uncancelled_work_scale"],
             "threshold": row[kind]["strong_action_strict_threshold"],
             "status": row[kind]["strong_action_strict_status"]}
            for row in preflight["cases"] for kind in ("linear", "nonlinear")]
    safety = {}
    for p, disc in sorted(discs.items()):
        q, velocity, _ = smooth_coordinates(disc)
        diagonal = np.linalg.eigvalsh(disc.mass_matrix(np.zeros(disc.ndof)))
        info = disc.diagnostics(q, velocity)
        safety[str(p)] = {"status": "PASS" if diagonal.min() > 0 and info["mass_positive"] and
                         info["min_one_plus_c"] > 0 and all(np.isfinite(x) for x in info.values()) else "FAIL",
                         "resting_mass_min_eigenvalue": float(diagonal.min()), "diagnostics": info}
    after = source_evidence(source_paths)
    preservation = {"status": "PASS" if after == evidence else "FAIL", "sources": evidence}
    checked_rows = [exact, linear, preservation, static]
    checked_rows += saved
    checked_rows += [part for row in variational.values() for key, part in row.items()
                     if key != "strong_projection_strict" and isinstance(part, dict) and "status" in part]
    checked_rows += list(safety.values())
    gate = "PASS" if all(row["status"] in ("PASS", "NOT_RUN") for row in checked_rows) else "FAIL"
    qualification = any(row["strong_projection_strict"]["status"] != "PASS" for row in variational.values())
    result = {"version": VERSION, "execution_gate": gate,
              "stage_a_status": "PASS_WITH_QUALIFICATIONS" if gate == "PASS" and qualification else gate,
              "implementation_sha256": {
                  filename: _sha(Path(__file__).with_name(filename)) for filename in
                  ("nlsp_spatial_verification_checks.py", "weakly_nonlinear_spatial_dynamics.py",
                   "weakly_nonlinear_planar_dynamics.py", "weakly_nonlinear_spatial_rod.py")},
              "exact_action": exact, "linear_references": linear, "saved_planar_state_checks": saved,
              "saved_planar_state_status": "PASS" if saved and all(r["status"] == "PASS" for r in saved) else "NOT_RUN" if not saved else "FAIL",
              "saved_static_equivalence": static, "variational_and_jacobian": variational,
              "numerical_mass_safety": safety, "source_preservation": preservation,
              "strict_float64_qualification": historical_strict,
              "STRICT_FLOAT64_STRONG_PROJECTION": {
                  "status": "PARTIAL" if qualification else "PASS",
                  "threshold": STRICT_RELATIVE,
                  "rows": {p: row["strong_projection_strict"] for p, row in variational.items()},
                  "interpretation": "historical strong differentiation/projection retained separately; no failed result promoted"},
              "coordinate_contract": {"fields": list(rod.FIELD_ORDER), "local_basis": [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
                                      "rotation_vector": ["Phi", "-psi", "theta"],
                                      "planar_indices": [0, 1, 5, 6], "derivative_endpoint_constraints": False},
              "scientific_calls": {"ODE_integrations": 0, "static_equilibrium_solves": 0, "FEM_jobs": 0,
                                   "Gmsh_jobs": 0, "BVP_solves": 0, "continuum_root_searches": 0,
                                   "new_symbolic_derivations": 0},
              "elapsed_seconds": time.perf_counter() - started,
              "qualification": "saved-state evaluation is not a new seven-field trajectory or temporal/spatial certificate"}
    if output is not None:
        path = Path(output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf8")
    return result
