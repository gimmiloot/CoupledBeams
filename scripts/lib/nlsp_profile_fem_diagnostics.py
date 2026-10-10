"""Saved-field spatial recovery diagnostics; no FEM or 1D solves.

The historical FEM-2 recovery is called unchanged.  Independent strains below
differentiate the existing C3D10 displacement interpolant at its audited
14-point positive reference-volume quadrature.  A local scalar stretch is a
diagnostic of this 3D field, not an identification of the M-H coordinate c.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import CubicSpline

from scripts.analysis import verify_nlsp_linear_rectangular_3d_fem as fem1
from scripts.analysis import verify_nlsp_nonlinear_static_3d_fem as fem2

LOCAL_SIGNS = np.array((1., -1., -1.))
DEFINITIONS = {
    "axes": "local (s,eta,zeta)=(global X,-global Y,-global Z); eta thickness, zeta width",
    "c_eff": "historical director polar stretch[0,0]-1; thickness eta; not M-H DOF",
    "c_small": "historical fitted local displacement gradient[1,1]",
    "width_effective": "historical director polar stretch[1,1]-1; width zeta",
    "width_small": "fitted local displacement gradient[2,2]",
    "native_small": "local gradient[1,1] from differentiated C3D10 displacements",
    "native_green": "E=0.5*(gradU+gradU.T+gradU.T@gradU), thickness E[1,1]",
    "weights": "positive reference mass weights/rho = undeformed volume; constant rho",
    "finite_strain_relation": "E_eta_eta=c_polar+0.5*c_polar**2+0.5*U_polar[0,1]**2 at same gradient",
    "rigid_rotation": "polar stretch-1 and Green strain vanish; small strain can be O(rotation**2)",
    "roughness": "physical-coordinate differences; TV=sum|dy|; second derivatives use unequal dx",
    "constraints": "historical audited end FACE values only; no derivative conditions",
    "slab_quadrature_qualification": "hard-binned element quadrature estimate for cut slabs, not exact clipped-tetrahedron integral; transverse sampling moments may bias means",
}


def prepare_mesh(mesh, rho=1.):
    """Reuse audited mesh/shape/quadrature and cache reference Jacobian inverses."""
    if not np.isfinite(rho) or rho <= 0:
        raise ValueError("Positive finite density required")
    ids, xyz, eids, conn = fem1.mesh_arrays(mesh)
    quad = fem1.quadrature_arrays(mesh, float(rho))
    _, dN = fem1.tet10_shape(quad["barycentric"])
    jac = np.einsum("eic,qij->eqcj", xyz[conn], dN)
    length, thickness, width = map(float, np.ptp(xyz, axis=0))
    if min(length, thickness, width) <= 0:
        raise ValueError("Three-dimensional rod geometry required")
    return {"quad": quad, "dN": dN, "jacobian_inverse": np.linalg.inv(jac),
            "node_ids": ids, "nodes": xyz, "element_ids": eids,
            "length": length, "thickness": thickness, "width": width,
            "rho": float(rho), "reference_volume": float(np.sum(quad["weights"])/rho)}


def native_quadrature_fields(prepared, nodal_displacement):
    """Independent local gradients/strains from existing nodal U, no polar fit."""
    U = np.asarray(nodal_displacement, float)
    if U.shape != (len(prepared["node_ids"]), 3) or not np.all(np.isfinite(U)):
        raise ValueError("Complete finite nodal displacement field required")
    quad = prepared["quad"]
    displacement = fem1.nlsp_evaluate_tet10_displacements(U, quad["conn"], quad["N"])
    derivative = np.einsum("eic,qij->eqcj", U[quad["conn"]], prepared["dN"])
    global_gradient = derivative @ prepared["jacobian_inverse"]
    gradient = global_gradient * LOCAL_SIGNS[:, None] * LOCAL_SIGNS[None, :]
    linear = .5*(gradient+gradient.swapaxes(-1, -2))
    green = linear+.5*(gradient.swapaxes(-1, -2)@gradient)
    deformation = np.eye(3)+gradient
    directors = deformation[..., :, 1:3]
    gram = directors.swapaxes(-1, -2)@directors
    eig, vec = np.linalg.eigh(gram)
    if np.any(eig <= 0):
        raise ValueError("Nonpositive native transverse director metric")
    stretch = (vec*np.sqrt(eig)[..., None, :])@vec.swapaxes(-1, -2)
    c_polar = stretch[..., 0, 0]-1.
    relation = c_polar+.5*c_polar*c_polar+.5*stretch[..., 0, 1]**2
    return {"global_points": quad["xyz"], "global_displacement": displacement,
        "local_points": fem1.nlsp_local_vectors(quad["xyz"]),
        "local_displacement": fem1.nlsp_local_vectors(displacement),
        "mass_weights": quad["weights"], "volume_weights": quad["weights"]/prepared["rho"],
        "gradient_local": gradient, "green_local": green, "linear_local": linear,
        "point_polar_thickness": c_polar, "point_polar_width": stretch[..., 1, 1]-1.,
        "point_polar_shear": stretch[..., 0, 1],
        "minimum_det_F": float(np.min(np.linalg.det(deformation))),
        "maximum_det_F": float(np.max(np.linalg.det(deformation))),
        "same_gradient_green_polar_relation_max_error": float(np.max(np.abs(green[..., 1, 1]-relation)))}


def profile_roughness(x, values):
    """Sampled roughness with physical length factors for unequal sample spacing."""
    x, y = np.asarray(x, float), np.asarray(values, float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 3 or np.any(np.diff(x) <= 0):
        raise ValueError("Ordered scalar profile with at least three samples required")
    dx = np.diff(x)
    slopes = np.diff(y)/dx
    second = 2*np.diff(slopes)/(dx[:-1]+dx[1:])
    second_weights = .5*(dx[:-1]+dx[1:])
    eps = 64*np.finfo(float).eps*max(float(np.max(np.abs(y))), np.finfo(float).tiny)
    signs = np.where(np.abs(np.diff(y)) > eps, np.sign(np.diff(y)), 0.)
    extrema = np.flatnonzero(signs[:-1]*signs[1:] < 0)+1
    return {"sampled_min": float(np.min(y)), "sampled_max": float(np.max(y)),
        "mean_over_sample_span": float(np.trapezoid(y, x)/(x[-1]-x[0])),
        "total_variation": float(np.sum(np.abs(np.diff(y)))),
        "first_derivative_L2": float(np.sqrt(np.sum(slopes**2*dx))),
        "second_derivative_RMS": float(np.sqrt(np.sum(second**2*second_weights)/sum(second_weights))),
        "second_difference_locations": x[1:-1], "second_derivative_samples": second,
        "raw_extrema_count": int(len(extrema)), "raw_extrema_x": x[extrema],
        "span": float(x[-1]-x[0]), "boundary_values_not_added": True}


def interpolation_diagnostics(raw_x, raw_values, length, dense_x=None, endpoint_values=(0., 0.)):
    """Exact CubicSpline derivative roots and unsmoothed linear comparison.

    End values reproduce the historical plot's audited clamps; extrema in the
    raw-center interval and the entire span are reported separately.
    """
    raw_x, raw = np.asarray(raw_x, float), np.asarray(raw_values, float)
    length = float(length)
    if raw.shape != raw_x.shape or len(raw) < 3 or np.any(np.diff(raw_x) <= 0):
        raise ValueError("Invalid ordered raw scalar profile")
    if not (0 < raw_x[0] < raw_x[-1] < length):
        raise ValueError("Raw section centers must be inside the material span")
    x = np.r_[0., raw_x, length]
    y = np.r_[float(endpoint_values[0]), raw, float(endpoint_values[1])]
    target = np.linspace(0., length, 801) if dense_x is None else np.asarray(dense_x, float)
    if target.ndim != 1 or np.any(np.diff(target) <= 0) or target[0] < 0 or target[-1] > length:
        raise ValueError("Ordered interpolation grid without extrapolation required")
    cubic = CubicSpline(x, y, extrapolate=False)
    roots = np.asarray(cubic.derivative().roots(extrapolate=False), float)
    roots = np.unique(roots[np.isfinite(roots) & (roots > 0) & (roots < length)])
    stationary = []
    for root in roots:
        k = min(len(x)-2, max(0, int(np.searchsorted(x, root)-1)))
        probe = max(32*np.finfo(float).eps*length, min(root-x[k], x[k+1]-root, x[k+1]-x[k])*.001)
        if probe == 0:
            probe = (x[k+1]-x[k])*.00001
        if cubic.derivative()(root-probe)*cubic.derivative()(root+probe) >= 0:
            continue
        val = float(cubic(root))
        lo, hi = sorted((float(y[k]), float(y[k+1])))
        stationary.append({"x": float(root), "value": val, "interval": k,
            "interval_left": float(x[k]), "interval_right": float(x[k+1]),
            "outside_neighbor_range": float(max(lo-val, val-hi, 0.)),
            "inside_raw_center_span": bool(raw_x[0] <= root <= raw_x[-1])})
    inside = [row for row in stationary if row["inside_raw_center_span"]]
    raw_metrics = profile_roughness(raw_x, raw)
    return {"dense_x": target, "cubic": cubic(target), "linear": np.interp(target, x, y),
        "raw_knot_reproduction_max_error": float(np.max(np.abs(cubic(raw_x)-raw))),
        "cubic_minus_linear_max": float(np.max(np.abs(cubic(target)-np.interp(target, x, y)))),
        "stationary_extrema": stationary, "interior_extrema_count": len(inside),
        "additional_interior_extrema_count": max(0, len(inside)-raw_metrics["raw_extrema_count"]),
        "maximum_neighbor_range_overshoot": max((row["outside_neighbor_range"] for row in stationary), default=0.),
        "raw_roughness": raw_metrics, "clamped_endpoint_values": list(endpoint_values),
        "smoothed": False, "interpolated_values_are_native_samples": False}


def _average_statistics(value, weight):
    mean = np.average(value, axis=0, weights=weight)
    variance = np.average((value-mean)**2, axis=0, weights=weight)
    return mean, np.sqrt(variance), np.min(value, axis=0), np.max(value, axis=0)


def _correlation(first, second):
    a, b = np.asarray(first, float), np.asarray(second, float)
    aa, bb = a-a.mean(), b-b.mean()
    denominator = np.linalg.norm(aa)*np.linalg.norm(bb)
    return float(np.dot(aa, bb)/denominator) if denominator > 0 else None


def recover_level(prepared, native, section_count=41, dense_x=None, enforce_clamped_faces=True):
    """Historical fit plus independent volume-strain means in the same slabs."""
    n = int(section_count)
    profile = fem2.fem2_recover_reference_samples(native["global_points"], native["global_displacement"],
        native["mass_weights"], prepared["length"], prepared["thickness"], prepared["width"], n,
        enforce_clamped_faces=enforce_clamped_faces)
    xyz = np.asarray(native["local_points"]).reshape(-1, 3)
    disp = np.asarray(native["local_displacement"]).reshape(-1, 3)
    gradient = native["gradient_local"].reshape(-1, 3, 3)
    green = native["green_local"].reshape(-1, 3, 3)
    linear = native["linear_local"].reshape(-1, 3, 3)
    volume = np.asarray(native["volume_weights"]).reshape(-1)
    mass = np.asarray(native["mass_weights"]).reshape(-1)
    values = np.column_stack((gradient[:, 1, 1], gradient[:, 2, 2], green[:, 1, 1], green[:, 2, 2],
        linear[:, 1, 2], green[:, 1, 2], green[:, 0, 0], gradient[:, 0, 0],
        np.asarray(native["point_polar_thickness"]).reshape(-1),
        np.asarray(native["point_polar_width"]).reshape(-1)))
    names = ("thickness_small", "width_small", "thickness_green", "width_green",
             "transverse_small_shear_tensor", "transverse_green_shear_tensor", "axial_green", "axial_small",
             "point_polar_thickness", "point_polar_width")
    bins = np.clip((xyz[:, 0]/prepared["length"]*n).astype(int), 0, n-1)
    rows = profile["section_rows"]
    means, stds, mins, maxs, displacement_rms, gradient_residuals = [], [], [], [], [], []
    component_displacement_rms, transverse_affine_residual_rms = [], []
    geometric_centroid, geometric_covariance, geometric_third_central = [], [], []
    for k, row in enumerate(rows):
        use = bins == k
        wt = volume[use]
        stats = _average_statistics(values[use], wt)
        for out, stat in zip((means, stds, mins, maxs), stats): out.append(stat)
        p, u = xyz[use], disp[use]
        pc = np.average(p, axis=0, weights=wt)
        geometric_centroid.append(pc)
        geometric_covariance.append(np.einsum("n,ni,nj->ij", wt, p-pc, p-pc)/wt.sum())
        geometric_third_central.append(np.average((p-pc)**3, axis=0, weights=wt))
        span = prepared["length"]/n
        d = (p[:, 0]-row["x"])/span
        eta, zeta = p[:, 1]/prepared["thickness"], p[:, 2]/prepared["width"]
        design = np.column_stack((np.ones(len(d)), d, d*d, d*d*d, eta, zeta, d*eta, d*zeta))
        sw = np.sqrt(mass[use]/mass[use].sum())
        fitted = np.linalg.lstsq(design*sw[:, None], u*sw[:, None], rcond=1e-12)[0]
        residual = u-design@fitted
        component_displacement_rms.append(np.sqrt(np.average(residual**2, axis=0, weights=wt)))
        displacement_rms.append(float(np.sqrt(np.average(np.sum(residual**2, axis=1), weights=wt))))
        fitted_gradient = np.empty((len(d), 3, 3))
        fitted_gradient[:, :, 0] = (fitted[1]+2*d[:, None]*fitted[2]+3*d[:, None]**2*fitted[3]
            +eta[:, None]*fitted[6]+zeta[:, None]*fitted[7])/span
        fitted_gradient[:, :, 1] = (fitted[4]+d[:, None]*fitted[6])/prepared["thickness"]
        fitted_gradient[:, :, 2] = (fitted[5]+d[:, None]*fitted[7])/prepared["width"]
        gres = gradient[use]-fitted_gradient
        gradient_residuals.append(np.sqrt(np.average(gres**2, axis=0, weights=wt)))
        transverse_affine_residual_rms.append(float(np.sqrt(np.average(np.sum(gres[:, :, 1:3]**2, axis=(1, 2)), weights=wt))))
    raw_x = np.asarray([row["x"] for row in rows])
    raw_fields = profile["fields"][1:-1] if enforce_clamped_faces else profile["fields"]
    raw_c = raw_fields[:, 6]
    small = np.asarray([row["affine_gradient_local"][1, 1] for row in rows])
    width = np.asarray([row["effective_width_strain"] for row in rows])
    width_small = np.asarray([row["affine_gradient_local"][2, 2] for row in rows])
    native_mean = np.asarray(means)
    columns = {name: native_mean[:, i] for i, name in enumerate(names)}
    interpolation = interpolation_diagnostics(raw_x, raw_c, prepared["length"], dense_x)
    all_roughness = {name: profile_roughness(raw_x, data) for name, data in {
        "c_eff": raw_c, "c_small": small, "width_effective": width,
        "native_thickness_small": columns["thickness_small"],
        "native_thickness_green": columns["thickness_green"],
        "native_point_polar_thickness": columns["point_polar_thickness"]}.items()}
    mean_difference = raw_c-columns["point_polar_thickness"]
    scalar_volume = np.asarray([float(volume[bins == k].sum()) for k in range(n)])
    fit_grad = np.asarray([row["affine_gradient_local"] for row in rows])
    stretch = np.asarray([row["transverse_stretch"] for row in rows])
    fit_green = .5*(fit_grad+fit_grad.swapaxes(-1, -2)+fit_grad.swapaxes(-1, -2)@fit_grad)
    fit_relation = raw_c+.5*raw_c**2+.5*stretch[:, 0, 1]**2
    geometric_centroid = np.asarray(geometric_centroid)
    nominal_volume = prepared["reference_volume"]/n
    volume_ratio = scalar_volume/nominal_volume
    correlations = {"c_eff_minus_native_point_polar_mean_vs_"+name: _correlation(mean_difference, value)
        for name, value in {"fit_condition": np.asarray([r["fit_condition"] for r in rows]),
            "fit_residual_RMS": displacement_rms, "slab_volume_ratio": volume_ratio,
            "thickness_sample_centroid": geometric_centroid[:, 1],
            "width_sample_centroid": geometric_centroid[:, 2],
            "native_small_std": np.asarray(stds)[:, 0]}.items()}
    return {"section_count": n, "historical_profile": profile, "raw_x": raw_x, "raw_fields": raw_fields,
        "c_eff": raw_c, "c_small": small, "width_effective": width, "width_small": width_small,
        "fit_sample_count": np.asarray([row["samples"] for row in rows]),
        "fit_mass": np.asarray([row["reference_mass"] for row in rows]), "fit_volume": scalar_volume,
        "geometric_centroid_local": geometric_centroid,
        "geometric_covariance_local": np.asarray(geometric_covariance),
        "geometric_third_central_moment_local": np.asarray(geometric_third_central),
        "slab_volume_to_nominal_ratio": volume_ratio, "nominal_slab_volume": nominal_volume,
        "fit_rank": np.asarray([row["fit_rank"] for row in rows]),
        "fit_condition": np.asarray([row["fit_condition"] for row in rows]),
        "fit_residual_mass_L2": np.asarray([row["section_residual_L2"] for row in rows]),
        "fit_residual_displacement_RMS": np.asarray(displacement_rms),
        "fit_residual_component_displacement_RMS": np.asarray(component_displacement_rms),
        "fit_gradient_component_residual_RMS": np.asarray(gradient_residuals),
        "fit_transverse_gradient_residual_RMS": np.asarray(transverse_affine_residual_rms),
        "native_measure_names": list(names), "native_mean": native_mean,
        "native_std": np.asarray(stds), "native_min": np.asarray(mins), "native_max": np.asarray(maxs),
        "native_columns": columns, "interpolation": interpolation, "roughness": all_roughness,
        "associations_not_causal_proof": correlations,
        "summary": {"maximum_fit_condition": profile["maximum_fit_condition"],
            "minimum_fit_rank": int(min(row["fit_rank"] for row in rows)),
            "section_residual_relative_L2": profile["section_residual_relative_L2"],
            "mean_contraction_volume_weighted": float(np.average(raw_c, weights=scalar_volume)),
            "mean_native_small_volume_weighted": float(np.average(columns["thickness_small"], weights=scalar_volume)),
            "mean_native_green_volume_weighted": float(np.average(columns["thickness_green"], weights=scalar_volume)),
            "c_eff_minus_c_small_max": float(np.max(np.abs(raw_c-small))),
            "c_eff_minus_native_point_polar_mean_max": float(np.max(np.abs(mean_difference))),
            "c_eff_minus_native_point_polar_mean_L2": float(np.sqrt(np.trapezoid(mean_difference**2, raw_x))),
            "c_small_minus_native_gradient_mean_max": float(np.max(np.abs(small-columns["thickness_small"]))),
            "fit_green_polar_relation_error_max": float(np.max(np.abs(fit_green[:, 1, 1]-fit_relation))),
            "native_transverse_heterogeneity_std_max": float(np.max(np.asarray(stds)[:, 0])),
            "native_transverse_heterogeneity_range_max": float(np.max(np.asarray(maxs)[:, 0]-np.asarray(mins)[:, 0])),
            "continuous_supremum_claimed": False}}


def compare_levels(first, second, length=1., dense_x=None):
    """Common-coordinate unsmoothed linear and original cubic comparisons."""
    length = float(length)
    x = np.linspace(0., length, 801) if dense_x is None else np.asarray(dense_x, float)
    result = {"x": x, "metrics": {}, "raw_roughness": {}, "alignment_used": False,
              "different_sample_centers": not np.array_equal(first["raw_x"], second["raw_x"])}
    for name in ("c_eff", "c_small", "width_effective", "width_small"):
        a, b = first[name], second[name]
        ac = interpolation_diagnostics(first["raw_x"], a, length, x)
        bc = interpolation_diagnostics(second["raw_x"], b, length, x)
        result["metrics"][name] = {"linear": fem2.fem2_curve_difference(x, ac["linear"], bc["linear"]),
            "cubic": fem2.fem2_curve_difference(x, ac["cubic"], bc["cubic"])}
        result["raw_roughness"][name] = {"first": profile_roughness(first["raw_x"], a),
            "second": profile_roughness(second["raw_x"], b)}
    pa, pb = first["historical_profile"], second["historical_profile"]
    a, b = fem2.fem2_static_sample(pa, x), fem2.fem2_static_sample(pb, x)
    result["control_fields"] = {name: fem2.fem2_curve_difference(x, a[:, k], b[:, k])
        for name, k in (("u", 0), ("w", 1), ("theta", 5))}
    # Direct strain averages have no prescribed zero FACE strain.  Compare only
    # within the intersection of raw-center coverage, never extrapolate them.
    lo = max(first["raw_x"][0], second["raw_x"][0])
    hi = min(first["raw_x"][-1], second["raw_x"][-1])
    interior_x = x[(x >= lo) & (x <= hi)]
    result["native_strain_comparison_span"] = [float(lo), float(hi)]
    result["native_measures"] = {}
    for name in first["native_measure_names"]:
        aa = np.interp(interior_x, first["raw_x"], first["native_columns"][name])
        bb = np.interp(interior_x, second["raw_x"], second["native_columns"][name])
        result["native_measures"][name] = fem2.fem2_curve_difference(interior_x, aa, bb)
    return result


def quadratic_transverse_fit_probe(prepared, native, section_count=41, dense_x=None):
    """ONE separate 11-column local-fit probe, not a replacement recovery.

    Add eta², eta*zeta, zeta² to the historical eight columns.  Preserve the
    same samples, weights, slab centers, longitudinal terms and director polar
    extraction.  The probe tests omitted transverse P2 displacement structure;
    it supplies neither a physical model nor new FE/dynamic constraints.
    """
    n = int(section_count)
    if n not in (21, 41, 81):
        raise ValueError("Probe uses only the declared 21/41/81 windows")
    xyz = np.asarray(native["local_points"]).reshape(-1, 3)
    disp = np.asarray(native["local_displacement"]).reshape(-1, 3)
    wtall = np.asarray(native["mass_weights"]).reshape(-1)
    volume = np.asarray(native["volume_weights"]).reshape(-1)
    gradient = np.asarray(native["gradient_local"]).reshape(-1, 3, 3)
    bins = np.clip((xyz[:, 0]/prepared["length"]*n).astype(int), 0, n-1)
    rows, fields = [], []
    for k in range(n):
        selected = bins == k
        p, u, wt = xyz[selected], disp[selected], wtall[selected]
        if len(wt) < 11:
            return {"status": "PARTIAL", "reason": "Insufficient samples for 11-column probe", "slab": k}
        xc = float(np.average(p[:, 0], weights=wt))
        span = prepared["length"]/n
        d = (p[:, 0]-xc)/span
        eta, zeta = p[:, 1]/prepared["thickness"], p[:, 2]/prepared["width"]
        design = np.column_stack((np.ones(len(d)), d, d*d, d*d*d,
            eta, zeta, d*eta, d*zeta, eta*eta, eta*zeta, zeta*zeta))
        sw = np.sqrt(wt/wt.sum())
        fitted, _, rank, singular = np.linalg.lstsq(design*sw[:, None], u*sw[:, None], rcond=1e-12)
        if rank != 11:
            return {"status": "PARTIAL", "reason": "Rank-deficient 11-column probe", "slab": k, "rank": int(rank)}
        condition = float(singular[0]/singular[-1])
        residual = u-design@fitted
        # Nested weighted regression identity: the 11-column residual is
        # orthogonal to the old eight columns.  Their coefficient change is
        # the aliasing of the three omitted transverse quadratic columns.
        old_design = design[:, :8]
        quadratic_component = design[:, 8:11]@fitted[8:11]
        alias = np.linalg.lstsq(old_design*sw[:, None], quadratic_component*sw[:, None], rcond=1e-12)[0]
        old_fitted = np.linalg.lstsq(old_design*sw[:, None], u*sw[:, None], rcond=1e-12)[0]
        alias_identity_error = np.max(abs(old_fitted-fitted[:8]-alias))
        alias_center_gradient = np.column_stack((alias[1]/span,
            alias[4]/prepared["thickness"], alias[5]/prepared["width"]))
        center_gradient = np.column_stack((fitted[1]/span,
            fitted[4]/prepared["thickness"], fitted[5]/prepared["width"]))
        orient = fem2.fem2_polar_section_orientation((np.eye(3)+center_gradient)[:, 1:3])
        point_gradient = np.empty((len(wt), 3, 3))
        point_gradient[:, :, 0] = (fitted[1]+2*d[:, None]*fitted[2]+3*d[:, None]**2*fitted[3]
            +eta[:, None]*fitted[6]+zeta[:, None]*fitted[7])/span
        point_gradient[:, :, 1] = (fitted[4]+d[:, None]*fitted[6]+2*eta[:, None]*fitted[8]
            +zeta[:, None]*fitted[9])/prepared["thickness"]
        point_gradient[:, :, 2] = (fitted[5]+d[:, None]*fitted[7]+eta[:, None]*fitted[9]
            +2*zeta[:, None]*fitted[10])/prepared["width"]
        gradient_residual = gradient[selected]-point_gradient
        fields.append(np.r_[fitted[0], orient["Phi"], orient["psi"], orient["theta"], orient["stretch"][0, 0]-1.])
        rows.append({"x": xc, "samples": len(wt), "rank": int(rank), "condition": condition,
            "mass": float(wt.sum()), "volume": float(volume[selected].sum()),
            "displacement_residual_mass_L2": float(np.sqrt(np.sum(wt*np.sum(residual**2, axis=1)))),
            "displacement_residual_RMS": float(np.sqrt(np.average(np.sum(residual**2, axis=1), weights=wt))),
            "transverse_gradient_residual_RMS": float(np.sqrt(np.average(
                np.sum(gradient_residual[:, :, 1:3]**2, axis=(1, 2)), weights=wt))),
            "center_gradient_local": center_gradient,
            "fitted_quadratic_transverse_coefficients": fitted[8:11],
            "quadratic_component_displacement_RMS": float(np.sqrt(np.average(np.sum(quadratic_component**2, axis=1), weights=wt))),
            "omitted_quadratic_alias_center_gradient_local": alias_center_gradient,
            "omitted_quadratic_alias_c_small": float(alias_center_gradient[1, 1]),
            "nested_fit_alias_coefficient_identity_max_error": float(alias_identity_error),
            "c_small": float(center_gradient[1, 1]), "width_small": float(center_gradient[2, 2]),
            "width_effective": float(orient["stretch"][1, 1]-1.)})
    x, fields = np.asarray([r["x"] for r in rows]), np.asarray(fields)
    c, small = fields[:, 6], np.asarray([r["c_small"] for r in rows])
    return {"status": "COMPLETE_DIAGNOSTIC_PROBE", "section_count": n, "raw_x": x,
        "raw_fields": fields, "c_eff": c, "c_small": small,
        "width_effective": np.asarray([r["width_effective"] for r in rows]),
        "width_small": np.asarray([r["width_small"] for r in rows]),
        "fit_rows": rows, "fit_rank": np.asarray([r["rank"] for r in rows]),
        "fit_condition": np.asarray([r["condition"] for r in rows]),
        "fit_residual_mass_L2": np.asarray([r["displacement_residual_mass_L2"] for r in rows]),
        "fit_residual_displacement_RMS": np.asarray([r["displacement_residual_RMS"] for r in rows]),
        "fit_transverse_gradient_residual_RMS": np.asarray([r["transverse_gradient_residual_RMS"] for r in rows]),
        "omitted_quadratic_alias_center_gradient_local": np.asarray([r["omitted_quadratic_alias_center_gradient_local"] for r in rows]),
        "omitted_quadratic_alias_c_small": np.asarray([r["omitted_quadratic_alias_c_small"] for r in rows]),
        "nested_fit_alias_coefficient_identity_max_error": float(max(r["nested_fit_alias_coefficient_identity_max_error"] for r in rows)),
        "roughness": {"c_eff": profile_roughness(x, c), "c_small": profile_roughness(x, small)},
        "interpolation": interpolation_diagnostics(x, c, prepared["length"], dense_x),
        "policy": "separate_11col_fit_add_eta2_etazeta_zeta2_to_same_historical_samples",
        "replaces_historical_recovery": False, "new_physical_assumption": False,
        "added_dynamic_constraints": False, "smoothed": False, "scientific_calls": 0}


def analyze_state(prepared, nodal_displacement, section_counts=(21, 41, 81), dense_x=None,
                  keep_sample_arrays=False, quadratic_fit_probe=False):
    """Bounded postprocessing of one supplied actual saved nodal state."""
    counts = tuple(map(int, section_counts))
    if not counts or len(set(counts)) != len(counts) or any(n not in (21, 41, 81) for n in counts):
        raise ValueError("Audit allows only distinct 21/41/81 recovery levels")
    native = native_quadrature_fields(prepared, nodal_displacement)
    levels = {str(n): recover_level(prepared, native, n, dense_x) for n in counts}
    comparisons = {f"{a}_vs_{b}": compare_levels(levels[str(a)], levels[str(b)], prepared["length"], dense_x)
        for a, b in zip(counts[:-1], counts[1:])}
    result = {"definitions": dict(DEFINITIONS), "levels": levels, "window_comparisons": comparisons,
        "native_diagnostics": {k: native[k] for k in ("minimum_det_F", "maximum_det_F",
            "same_gradient_green_polar_relation_max_error")},
        "reference_nodes": int(len(prepared["node_ids"])), "elements": int(len(prepared["element_ids"])),
        "reference_volume": prepared["reference_volume"], "scientific_calls": 0}
    if 41 in counts and 81 in counts:
        result["window_comparisons"]["41_vs_81"] = compare_levels(levels["41"], levels["81"], prepared["length"], dense_x)
    if keep_sample_arrays:
        result["quadrature_samples"] = {k: native[k] for k in ("local_points", "volume_weights",
            "gradient_local", "linear_local", "green_local", "point_polar_thickness", "point_polar_width")}
    if quadratic_fit_probe:
        result["quadratic_transverse_fit_probe"] = {str(n): quadratic_transverse_fit_probe(
            prepared, native, n, dense_x) for n in counts}
        for n in counts:
            probe = result["quadratic_transverse_fit_probe"][str(n)]
            if probe["status"] == "COMPLETE_DIAGNOSTIC_PROBE":
                primary = levels[str(n)]
                difference = primary["c_eff"]-probe["c_eff"]
                probe["historical_minus_probe_raw_max"] = float(np.max(abs(difference)))
                probe["historical_minus_probe_raw_L2"] = float(np.sqrt(np.trapezoid(difference**2, primary["raw_x"])))
    return result


def quadratic_transverse_sampling_control(prepared, section_counts=(21, 41, 81), coefficient=None,
                                         include_quadratic_fit_probe=True):
    """Exactly represented P2 kinematic probe on the supplied existing mesh.

    U_eta=a*eta**2, all other components zero.  This is an analytic synthetic
    displacement, not a new equilibrium or a scientific FEM job.  It separates
    exact C3D10 differentiation from hard-bin moment bias and omission of eta²
    in the historical director fit.  a defaults to 1/L and has units 1/length.
    """
    a = 1./prepared["length"] if coefficient is None else float(coefficient)
    local = fem1.nlsp_local_vectors(prepared["nodes"])
    local_U = np.zeros_like(local)
    local_U[:, 1] = a*local[:, 1]**2
    native = native_quadrature_fields(prepared, local_U*LOCAL_SIGNS)
    exact_gradient = 2*a*native["local_points"][..., 1]
    levels = {}
    for n in section_counts:
        out = recover_level(prepared, native, int(n), enforce_clamped_faces=False)
        expected = 2*a*out["geometric_centroid_local"][:, 1]
        covariance = out["geometric_covariance_local"][:, 1, 1]
        third = out["geometric_third_central_moment_local"][:, 1]
        # The scalar regression identity applies if other historical columns
        # do not correlate with eta.  Actual multi-column LS remains primary.
        scalar_prediction = expected+a*third/covariance
        levels[str(n)] = {"raw_x": out["raw_x"], "fit_c_small": out["c_small"],
            "fit_c_eff": out["c_eff"], "native_small_mean": out["native_columns"]["thickness_small"],
            "exact_bin_mean_from_centroid": expected,
            "scalar_affine_slope_if_no_other_correlated_columns": scalar_prediction,
            "fit_condition": out["fit_condition"],
            "geometric_centroid_local": out["geometric_centroid_local"],
            "geometric_covariance_local": out["geometric_covariance_local"],
            "geometric_third_central_moment_local": out["geometric_third_central_moment_local"],
            "fit_minus_direct_mean": out["c_small"]-expected,
            "fit_minus_direct_mean_max": float(np.max(abs(out["c_small"]-expected))),
            "native_vs_exact_bin_mean_max_error": float(np.max(abs(out["native_columns"]["thickness_small"]-expected))),
            "fit_c_roughness": out["roughness"]["c_small"],
            "direct_mean_roughness": out["roughness"]["native_thickness_small"]}
        if include_quadratic_fit_probe:
            probe = quadratic_transverse_fit_probe(prepared, native, int(n))
            levels[str(n)]["quadratic_fit_probe"] = probe
            if probe["status"] == "COMPLETE_DIAGNOSTIC_PROBE":
                levels[str(n)]["quadratic_probe_center_gradient_max_error"] = float(np.max(abs(probe["c_small"])))
    return {"definition": "synthetic U_eta=a*eta^2 on original nodes; no solver, no mesh change",
        "coefficient_inverse_length": a, "levels": levels,
        "C3D10_exact_gradient_max_error": float(np.max(abs(native["gradient_local"][..., 1, 1]-exact_gradient))),
        "exact_symmetric_rectangular_section_average": 0.,
        "quadrature_clipped_slab_averages_can_be_nonzero": True,
        "no_scientific_equilibrium_claim": True, "scientific_calls": 0}
