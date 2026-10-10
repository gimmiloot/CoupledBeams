"""FEM-3C saved-data robustness comparisons and presentation.

No native, ODE, equilibrium or eigenvalue solver is called here. Positive-time
native observations and the confirmed STATIC initial state remain distinct.
The only temporal transfers are the preregistered linear/PCHIP postprocessing
operators; they do not fit phase, amplitude, period or material coefficients.
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d, PchipInterpolator

from scripts.lib import nlsp_fem3b_diagnostics as old

VERSION = "fem3c-saved-robustness-diagnostics-v1"
CANONICAL_FIELDS = ("u", "w", "v", "Phi", "psi", "theta", "c_eff_diagnostic")
ONE_D_FIELD_ORDER = ("u", "w", "v", "Phi", "psi", "theta", "c")
ACTIVE_FIELDS = {"u": (0, 0), "w": (1, 1), "theta": (2, 5),
                 "c_eff_diagnostic": (3, 6)}
BASELINE_MODEL_DISCREPANCY = 2.824717e-7
HISTORICAL_EVOLUTION_SCALE = 3.86809158318442e-6
read_json, write_json, load_arrays = old.read_json, old.write_json, old.load_arrays


def comparison_grid(policy, T1):
    """The common physical grid is fixed before the new native results exist."""
    if (policy["T1_fraction"] != .25 or policy["points"] != 201
            or policy["primary_interpolation"] != "linear"
            or policy["diagnostic_interpolation"] != "PCHIP"
            or policy["baseline_model_discrepancy"] != BASELINE_MODEL_DISCREPANCY
            or policy["temporal_ratio_limit"] != .25
            or policy["spatial_ratio_limit"] != .25
            or policy["interpolation_ratio_to_effect_limit"] != .25
            or policy["interpolation_ratio_to_baseline_limit"] != .25
            or policy["phase_amplitude_fitting"] is not False
            or not np.isfinite(T1) or T1 <= 0):
        raise ValueError("Preregistered FEM-3C comparison policy changed")
    return np.linspace(0., float(T1) * .25, 201)


def seven_field_one_d(fields):
    """Embed four active fields in the canonical planar invariant subspace.

    The three zeros are the assumed 1D planar subspace, not independently
    validated spatial dynamics. The final coordinate is actual 1D c; its 3D
    counterpart must remain labelled as an effective contraction proxy.
    """
    fields = np.asarray(fields)
    if fields.shape[-1] != 4 or not np.isfinite(fields).all():
        raise ValueError("Expected finite fields ordered (u,w,theta,c)")
    result = np.zeros(fields.shape[:-1] + (7,), dtype=fields.dtype)
    result[..., [0, 1, 5, 6]] = fields
    return result


def transfer_case(case_root, grid, method="linear"):
    """Transfer recorded U/V profiles, with true STATIC provenance at t=0."""
    case_root, grid = Path(case_root), np.asarray(grid, dtype=float)
    history = load_arrays(case_root / "section_history.npz")
    initial = load_arrays(case_root / "initial_sections.npz")
    time, x, fields = history["time"], history["x"], history["fields"]
    if (grid.ndim != 1 or len(grid) < 2 or grid[0] != 0.
            or np.any(np.diff(grid) <= 0) or np.any(np.diff(time) <= 0)
            or time[0] <= 0 or grid[-1] > time[-1]
            or fields.shape != (len(time), len(x), 7)
            or not np.array_equal(x, initial["x"])
            or initial["fields"].shape != (len(x), 7)
            or not bool(initial["not_a_native_dynamic_zero_frame"])
            or not np.isfinite(fields).all() or not np.isfinite(initial["fields"]).all()
            or not np.isfinite(grid).all() or not np.isfinite(time).all()
            or np.any(np.diff(x) <= 0)):
        raise ValueError("Saved native/static profiles or their actual coverage are invalid")
    native_time = np.r_[0., time]
    values = np.concatenate([initial["fields"][None, ...], fields])
    if method == "linear":
        operation = lambda value: interp1d(native_time, value, axis=0, bounds_error=True)(grid)
    elif method == "PCHIP":
        operation = lambda value: PchipInterpolator(native_time, value, axis=0, extrapolate=False)(grid)
    else:
        raise ValueError("Only preregistered linear and PCHIP temporal transfers are permitted")
    transferred = operation(values)
    # Enforce the actual static datum exactly: interpolation must not invent t=0.
    transferred[0] = initial["fields"]
    result = {"time": grid.copy(), "x": x.copy(), "fields": transferred,
        "initial_fields": initial["fields"], "native_time": time.copy(),
        "metadata": {"case_path": str(case_root), "interpolation": method,
            "actual_native_interval": [float(time[0]), float(time[-1])],
            "native_samples": len(time), "comparison_samples": len(grid),
            "zero_source": "confirmed_STATIC_preload_not_native_DYNAMIC_frame",
            "positive_comparison_values": "postprocessing_interpolated",
            "source_static_end_time": float(initial["source_static_end_time"]),
            "physical_dynamic_time_origin": 0., "extrapolation": False,
            "phase_amplitude_alignment": False}}
    if "translation_velocities" in history:
        velocity = history["translation_velocities"]
        if velocity.shape != (len(time), len(x), 3) or not np.isfinite(velocity).all():
            raise ValueError("Incomplete/nonfinite recorded translation velocities")
        result["translation_velocities"] = operation(np.concatenate([np.zeros_like(velocity[:1]), velocity]))
        result["translation_velocities"][0] = 0.
    return result


def correction_pair(linear, nonlinear):
    if (not np.array_equal(linear["time"], nonlinear["time"])
            or not np.array_equal(linear["x"], nonlinear["x"])):
        raise ValueError("Corrections require the same physical time/material section grids")
    correction = nonlinear["fields"] - linear["fields"]
    initial = nonlinear["initial_fields"] - linear["initial_fields"]
    result = {"time": linear["time"], "x": linear["x"],
        "linear": linear["fields"], "nonlinear": nonlinear["fields"],
        "correction": correction, "initial_correction": initial,
        "evolution": correction - initial[None, ...]}
    for key in ("linear", "nonlinear"):
        case = linear if key == "linear" else nonlinear
        if "translation_velocities" in case:
            result[key + "_translation_velocities"] = case["translation_velocities"]
    return result


def _quantities(pair, fields=CANONICAL_FIELDS):
    table, arrays = {}, {}
    for quantity in ("linear", "nonlinear", "correction", "evolution"):
        table[quantity] = {}
        for j, field in enumerate(fields):
            metric, curves = old.sampled_norms(pair[quantity][..., j], pair["x"], pair["time"])
            table[quantity][field] = metric
            for name, value in curves.items():
                arrays[quantity + "_" + field + "_" + name] = value
    return table, arrays


def _difference_pair(first, second):
    if (not np.array_equal(first["time"], second["time"])
            or not np.array_equal(first["x"], second["x"])):
        raise ValueError("Comparison coordinates disagree")
    result = {"time": first["time"], "x": first["x"]}
    for key in ("linear", "nonlinear", "correction", "evolution"):
        result[key] = first[key] - second[key]
    return result


def interpolation_comparability(differences, Dt, Dh, policy):
    """Preregistered comparability heuristic; not an FEM continuum error bound."""
    largest = max(float(value) for value in differences.values())
    if (not np.isfinite([largest, Dt, Dh]).all() or min(largest, Dt, Dh) < 0):
        raise ValueError("Refinement and interpolation differences must be nonnegative finite")
    effect = min(Dt, Dh)
    limit_effect = policy["interpolation_ratio_to_effect_limit"] * effect
    limit_baseline = policy["interpolation_ratio_to_baseline_limit"] * policy["baseline_model_discrepancy"]
    passed = largest <= limit_effect and largest <= limit_baseline
    return {"status": "PASS" if passed else "INTERPOLATION_UNRESOLVED",
        "evolution_w_max_by_resolution": differences,
        "largest_evolution_w_max": largest, "smallest_refinement_effect": effect,
        "ratio_to_smallest_refinement_effect": largest / effect if effect else None,
        "ratio_to_baseline_model_discrepancy": largest / policy["baseline_model_discrepancy"],
        "limit_from_effect": limit_effect, "limit_from_baseline": limit_baseline,
        "zero_effect_requires_zero_difference": True,
        "diagnostic_comparability_only": True, "continuum_error_bound": False}


def _model_comparisons(pair, one_d, common_scales):
    result = {}
    for key in ("linear", "nonlinear", "correction", "evolution"):
        result[key] = {}
        for field, (i, j) in ACTIVE_FIELDS.items():
            a, b = one_d[key][..., i], pair[key][..., j]
            row = old.difference_metrics(a, b, pair["x"], pair["time"])
            row["fixed_all_resolutions_characteristic_scale"] = common_scales[key][field]
            row["relative_max_on_fixed_all_resolutions_scale"] = row["absolute_max"] / common_scales[key][field] if common_scales[key][field] else None
            row["relative_L2_on_fixed_all_resolutions_scale"] = row["max_time_L2"] / common_scales[key][field] if common_scales[key][field] else None
            if key == "evolution" and field == "w":
                row["historical_characteristic_scale"] = HISTORICAL_EVOLUTION_SCALE
                row["relative_max_on_historical_characteristic_scale"] = row["absolute_max"] / HISTORICAL_EVOLUTION_SCALE
                row["relative_L2_on_historical_characteristic_scale"] = row["max_time_L2"] / HISTORICAL_EVOLUTION_SCALE
                row["ratio_to_fixed_baseline_model_difference"] = row["absolute_max"] / BASELINE_MODEL_DISCREPANCY
            row["physical_equivalence"] = "effective_contraction_proxy_only" if field == "c_eff_diagnostic" else "section_projection_coordinate_contract"
            result[key][field] = row
    return result


def robustness(bundle, config, *, case_roots=None, parent_bundle=None, one_d_evaluator=None):
    """Read existing histories, save the bounded C1 comparison, never run jobs.

    ``one_d_evaluator(times)`` may supply cached exact-linear/dense-nonlinear
    data. By default the immutable FEM-3B saved factors and Radau polynomials
    are used, with the helper's explicit no-new-eigensolve check.
    """
    bundle = Path(bundle)
    policy = config["comparison"]
    T1 = config.get("T1", 10.37828159055014)
    grid = comparison_grid(policy, T1)
    if parent_bundle is None:
        parent_bundle = old.pilot.ROOT / config["parent_completed"]["bundle"]
    parent_bundle = Path(parent_bundle)
    if case_roots is None:
        case_roots = {"existing_medium": parent_bundle / "cases",
            "medium_refined_time": bundle / "cases" / "medium_refined_time",
            "fine_refined_time": bundle / "cases" / "fine_refined_time"}
    pairs, metadata = {}, {}
    arrays = {"times": grid}
    for label, root in case_roots.items():
        pairs[label], metadata[label] = {}, {}
        for method in ("linear", "PCHIP"):
            cases = {kind: transfer_case(Path(root) / kind, grid, method) for kind in ("linear", "nonlinear")}
            pairs[label][method] = correction_pair(cases["linear"], cases["nonlinear"])
            metadata[label][method] = {kind: value["metadata"] for kind, value in cases.items()}
    labels = ("existing_medium", "medium_refined_time", "fine_refined_time")
    if tuple(case_roots) != labels:
        raise ValueError("Robustness must contain old medium, time-refined medium and time-refined fine")
    x = pairs[labels[0]]["linear"]["x"]
    if any(not np.array_equal(x, pairs[label][method]["x"]) for label in labels for method in ("linear", "PCHIP")):
        raise ValueError("Fixed material section grid changed between resolutions")
    arrays["x"] = x
    if one_d_evaluator is None:
        from scripts.lib.nlsp_fem3b_continuation import evaluate_saved_one_d
        old_item = read_json(parent_bundle / "provenance.json")
        one_d_evaluator = lambda time: evaluate_saved_one_d(parent_bundle, old_item, time)
    reference = one_d_evaluator(grid)
    if not np.array_equal(reference["times"], grid) or not np.array_equal(reference["x"], x):
        raise ValueError("Saved 1D evaluator returned different physical comparison coordinates")
    one_d = {"linear": reference["linear_fields"], "nonlinear": reference["nonlinear_fields"]}
    one_d["correction"] = one_d["nonlinear"] - one_d["linear"]
    one_d["evolution"] = one_d["correction"] - (reference["initial_nonlinear_fields"] - reference["initial_linear_fields"])[None, ...]
    common_scales = {quantity: {field: max(float(abs(one_d[quantity][..., i]).max()),
        *(float(abs(pairs[label]["linear"][quantity][..., j]).max()) for label in labels))
        for field, (i, j) in ACTIVE_FIELDS.items()} for quantity in one_d}
    result = {"version": VERSION, "comparison_policy": policy,
        "physical_time_interval": [0., float(grid[-1])], "comparison_points": len(grid),
        "section_points": len(x), "sampled_maxima_only": True,
        "temporal_transfer": metadata, "fixed_all_resolutions_characteristic_scales": common_scales,
        "resolutions": {}, "refinement": {}, "scientific_calls": {"CCX": 0, "Gmsh": 0,
            "nonlinear_ODE": 0, "eigenanalysis": 0, "static_equilibrium": 0},
        "qualification": "observed changes of two FEM resolutions, not rigorous continuum error bounds; no phase/amplitude fitting; c_eff is diagnostic only"}
    interpolation_differences = {}
    for label in labels:
        pair = pairs[label]["linear"]
        table, curves = _quantities(pair)
        alternative_table, _ = _quantities(pairs[label]["PCHIP"])
        interpolation_table, interpolation_curves = _quantities(_difference_pair(pair, pairs[label]["PCHIP"]))
        interpolation_differences[label] = interpolation_table["evolution"]["w"]["absolute_max"]
        result["resolutions"][label] = {"quantities": table,
            "PCHIP_diagnostic_quantities": alternative_table,
            "linear_PCHIP_differences": interpolation_table,
            "model_comparisons": _model_comparisons(pair, one_d, common_scales),
            "PCHIP_model_comparisons": _model_comparisons(pairs[label]["PCHIP"], one_d, common_scales)}
        for key, value in curves.items(): arrays[label + "_" + key] = value
        for key, value in interpolation_curves.items(): arrays[label + "_interpolation_" + key] = value
        for quantity in one_d:
            arrays[label + "_" + quantity + "_fields"] = pair[quantity]
            arrays[label + "_PCHIP_" + quantity + "_fields"] = pairs[label]["PCHIP"][quantity]
    for name, first, second in (("temporal", labels[1], labels[0]), ("spatial", labels[2], labels[1])):
        table, curves = _quantities(_difference_pair(pairs[first]["linear"], pairs[second]["linear"]))
        pchip_table, _ = _quantities(_difference_pair(pairs[first]["PCHIP"], pairs[second]["PCHIP"]))
        effect = table["evolution"]["w"]["absolute_max"]
        ratio = effect / policy["baseline_model_discrepancy"]
        baseline_model = result["resolutions"][labels[0]]["model_comparisons"]["evolution"]["w"]["absolute_max"]
        updated_model = result["resolutions"][first]["model_comparisons"]["evolution"]["w"]["absolute_max"]
        result["refinement"][name] = {"first": first, "second": second,
            "quantities": table, "PCHIP_diagnostic_quantities": pchip_table,
            "evolution_w_max": effect, "ratio_to_fixed_baseline_model_discrepancy": ratio,
            "ratio_limit": policy[name + "_ratio_limit"],
            "status": "PASS" if ratio <= policy[name + "_ratio_limit"] else "PARTIAL",
            "updated_absolute_model_discrepancy": updated_model,
            "signed_absolute_model_discrepancy_change_from_existing": updated_model - baseline_model,
            "error_upper_bound_claimed": False}
        for key, value in curves.items(): arrays[name + "_" + key] = value
    result["interpolation_comparability"] = interpolation_comparability(interpolation_differences,
        result["refinement"]["temporal"]["evolution_w_max"],
        result["refinement"]["spatial"]["evolution_w_max"], policy)
    passed = all(result["refinement"][name]["status"] == "PASS" for name in ("temporal", "spatial")) and result["interpolation_comparability"]["status"] == "PASS"
    result["numerical_robustness_status"] = "PASS" if passed else "PARTIAL"
    result["full_period_numerical_robustness_gate"] = passed
    result["full_period_authorization_requires_other_gates_and_resource_preflight"] = True
    for key, value in one_d.items(): arrays["one_d_" + key + "_fields"] = value
    np.savez_compressed(bundle / "robustness_comparison.npz", **arrays)
    write_json(bundle / "robustness_comparison.json", result)
    _summary_csv(bundle / "robustness_summary.csv", result)
    return result


def _summary_csv(path, result):
    with Path(path).open("w", encoding="utf8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("resolution", "quantity", "field", "absolute_max", "max_time_L2",
            "signed_at_max", "x_at_max", "time_at_max", "signed_final_midspan",
            "model_absolute_max", "model_max_time_L2", "model_own_characteristic_scale",
            "model_relative_max_own_scale", "model_fixed_all_resolutions_scale",
            "model_relative_max_fixed_scale", "model_relative_max_historical_evolution_scale"))
        for label, record in result["resolutions"].items():
            for quantity, fields in record["quantities"].items():
                for field, metric in fields.items():
                    model = record["model_comparisons"][quantity].get(field, {})
                    writer.writerow((label, quantity, field, metric["absolute_max"], metric["max_time_L2"],
                        metric["signed_at_max"], metric["x_at_max"], metric["time_at_max"], metric["midspan_final"],
                        model.get("absolute_max"), model.get("max_time_L2"), model.get("characteristic_scale"),
                        model.get("relative_max"), model.get("fixed_all_resolutions_characteristic_scale"),
                        model.get("relative_max_on_fixed_all_resolutions_scale"),
                        model.get("relative_max_on_historical_characteristic_scale")))


def plot_robustness(bundle, T1=10.37828159055014):
    """One compact control figure from saved arrays, without solver calls."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bundle = Path(bundle)
    data = load_arrays(bundle / "robustness_comparison.npz")
    result = read_json(bundle / "robustness_comparison.json")
    fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    tau = data["times"] / T1
    for label, title in (("existing_medium", "medium, original time"),
            ("medium_refined_time", "medium, refined time"),
            ("fine_refined_time", "fine, refined time")):
        ax[0, 0].plot(tau, data[label + "_evolution_fields"][:, len(data["x"])//2, 1] * 1e6, label=title)
        ax[1, 0].plot(tau, data[label + "_interpolation_evolution_w_max"] * 1e9, label=title)
    ax[0, 0].plot(tau, data["one_d_evolution_fields"][:, len(data["x"])//2, 1] * 1e6, "k--", label="1D p64")
    for label in ("temporal", "spatial"):
        ax[0, 1].plot(tau, data[label + "_evolution_w_max"] / BASELINE_MODEL_DISCREPANCY, label=label)
    ax[0, 1].axhline(.25, color="k", linestyle=":", label="preregistered guide")
    labels = list(result["resolutions"])
    vals = [result["resolutions"][name]["model_comparisons"]["evolution"]["w"]["relative_max_on_historical_characteristic_scale"]*100 for name in labels]
    ax[1, 1].bar(range(3), vals)
    ax[1, 1].set_xticks(range(3), ["medium old", "medium time", "fine time"])
    for i, value in enumerate(vals): ax[1, 1].text(i, value, f"{value:.3f}%", ha="center", va="bottom")
    ax[0, 0].set_ylabel("evolving midspan correction × 10⁶")
    ax[0, 1].set_ylabel("observed change / fixed model discrepancy")
    ax[1, 0].set_ylabel("linear/PCHIP evolution difference × 10⁹")
    ax[1, 1].set_ylabel("model discrepancy / historical evolution scale (%)")
    for panel in ax.flat:
        panel.grid(alpha=.25)
        if panel is not ax[1, 1]:
            panel.set_xlabel("t / T₁")
            panel.legend(fontsize=8)
    fig.suptitle("FEM-3C1: bounded temporal and spatial robustness (sampled maxima)")
    paths = []
    for suffix in ("png", "pdf"):
        path = bundle / ("numerical_robustness_controls." + suffix)
        fig.savefig(path, dpi=180)
        paths.append(path)
    plt.close(fig)
    return paths


def _full_grid(config, T1):
    policy = config["full_period_comparison"]
    if (policy["points"] != 401 or policy["primary_interpolation"] != "linear"
            or policy["diagnostic_interpolation"] != "PCHIP"
            or policy["snapshot_T1_fractions"] != [0., .25, .5, .75, 1.]
            or policy["phase_amplitude_fitting"] is not False):
        raise ValueError("Preregistered illustrative full-period comparison changed")
    observations = config["seven_field_observations"]
    if ({name: observations[name] for name in ("u", "w", "v", "Phi", "psi", "theta", "c")}
            != {"u": .25, "w": .5, "v": .5, "Phi": .25, "psi": .25, "theta": .25, "c": .25}
            or observations["inactive_one_d_fields"] != ["v", "Phi", "psi"]
            or observations["three_d_c"] != "effective_contraction_proxy_not_identical_generalized_coordinate"):
        raise ValueError("Fixed seven-field observations/qualifications changed")
    return np.linspace(0., T1, 401)


def _prefix_observations(new_cases, old_cases):
    """Compare only exact native samples on the old interval, without fitting."""
    result = {"scope": "exact shared native times; STATIC datum compared separately",
        "old_interval": [0., .25 * 10.37828159055014], "cases": {}}
    for kind in ("linear", "nonlinear"):
        old_h = load_arrays(Path(old_cases)/kind/"section_history.npz")
        new_h = load_arrays(Path(new_cases)/kind/"section_history.npz")
        common, old_index, new_index = np.intersect1d(old_h["time"], new_h["time"], return_indices=True)
        old_i = load_arrays(Path(old_cases)/kind/"initial_sections.npz")
        new_i = load_arrays(Path(new_cases)/kind/"initial_sections.npz")
        initial = float(abs(old_i["fields"] - new_i["fields"]).max())
        row = {"shared_native_samples": len(common), "initial_static_fields_max_abs": initial,
            "new_output_sampling_differs": True, "unmatched_samples_not_interpolated": True,
            "fields": {}}
        if len(common):
            row["shared_interval"] = [float(common[0]), float(common[-1])]
            for j, field in enumerate(CANONICAL_FIELDS):
                row["fields"][field] = old.difference_metrics(new_h["fields"][new_index, :, j],
                    old_h["fields"][old_index, :, j], old_h["x"], common)
            row["exact_shared_profile_reproduction"] = initial == 0 and all(v["absolute_max"] == 0 for v in row["fields"].values())
        else:
            row["status"] = "NO_EXACT_NATIVE_OVERLAP_SEPARATE_FROM_INTERPOLATED_COMPARISON"
        result["cases"][kind] = row
    return result


def _saved_kinetic_comparison(case):
    """Read the existing C3D10 mass/quadrature result, without reintegration."""
    case = Path(case)
    path = case/"independent_kinetic_energy.npz"
    if not path.exists() or not (case/"energy.json").exists():
        return {"status": "NOT_RUN", "reason": "independent kinetic output unavailable"}
    data = load_arrays(path)
    history = load_arrays(case/"section_history.npz")
    records = read_json(case/"energy.json")["records"]
    native = {row["increment"]: row["kinetic_energy"] for row in records if row["step"] == 2}
    if (not np.array_equal(data["time"], history["time"])
            or len(data["kinetic_energy"]) != len(history["increments"])):
        raise ValueError("Saved independent kinetic energy/native increment coverage mismatch")
    try:
        printed = np.array([native[int(increment)] for increment in history["increments"]])
    except KeyError as error:
        raise ValueError("Native kinetic energy missing at a saved displacement frame") from error
    independent = data["kinetic_energy"]
    if not np.isfinite(independent).all() or not np.isfinite(printed).all():
        raise ValueError("Nonfinite saved independent/native kinetic energy")
    difference = independent-printed
    scale = float(max(abs(independent).max(), abs(printed).max()))
    return {"status": "COMPLETED_DIAGNOSTIC", "samples": len(independent),
        "definition": "0.5*integral_reference_rho*|N*actual_nodal_V|^2",
        "quadrature": "existing affine C3D10 positive14point degree5 consistent mass integral",
        "first_independent_K": float(independent[0]), "first_native_K": float(printed[0]),
        "final_independent_K": float(independent[-1]), "final_native_K": float(printed[-1]),
        "maximum_absolute_difference": float(abs(difference).max()),
        "fixed_kinetic_scale": scale,
        "maximum_relative_difference_fixed_K_scale": float(abs(difference).max()/scale) if scale else None,
        "internal_StVK_reconstruction": "NOT_RUN",
        "native_internal_energy_reference_corrected": False}


def full_period(bundle, item, summary, *, one_d_evaluator=None):
    """Illustrative full-period fields from saved CCX output and 1D histories.

    This extends neither the quarter-period convergence evidence nor its
    practical acceptance guide to T1. T1 remains the fixed first *linear*
    period; returning to the starting state is not a completion criterion.
    """
    bundle = Path(bundle)
    config = item["validation_config"]
    T1 = 2*np.pi/item["config"]["omega1"]
    grid = _full_grid(config, T1)
    cases = bundle/"cases"/"full_period_medium"
    pairs, metadata = {}, {}
    for method in ("linear", "PCHIP"):
        histories = {kind: transfer_case(cases/kind, grid, method) for kind in ("linear", "nonlinear")}
        pairs[method] = correction_pair(histories["linear"], histories["nonlinear"])
        metadata[method] = {kind: value["metadata"] for kind, value in histories.items()}
    if one_d_evaluator is None:
        from scripts.lib.nlsp_fem3c_1d import evaluate_saved
        one_d_evaluator = lambda time: evaluate_saved(bundle, item, time, p=64)
    reference = one_d_evaluator(grid)
    x = pairs["linear"]["x"]
    if not np.array_equal(reference["times"], grid) or not np.array_equal(reference["x"], x):
        raise ValueError("Full-period saved 1D/3D physical coordinates disagree")
    active = {"linear": reference["linear_fields"], "nonlinear": reference["nonlinear_fields"]}
    active["correction"] = active["nonlinear"]-active["linear"]
    active["evolution"] = active["correction"]-(reference["initial_nonlinear_fields"]-reference["initial_linear_fields"])[None, ...]
    arrays = {"times": grid, "x": x}
    common_scales = {quantity: {field: max(float(abs(active[quantity][..., i]).max()),
        float(abs(pairs["linear"][quantity][..., j]).max()))
        for field, (i, j) in ACTIVE_FIELDS.items()} for quantity in active}
    table, curves = _quantities(pairs["linear"])
    interpolation_table, interpolation_curves = _quantities(_difference_pair(pairs["linear"], pairs["PCHIP"]))
    model = _model_comparisons(pairs["linear"], active, common_scales)
    observations = config["seven_field_observations"]
    observation_data, inactive = {}, {}
    active_displacement_scale = max(float(abs(active["nonlinear"][..., 1]).max()),
        float(abs(pairs["linear"]["nonlinear"][..., 1]).max()))
    active_rotation_scale = max(float(abs(active["nonlinear"][..., 2]).max()),
        float(abs(pairs["linear"]["nonlinear"][..., 5]).max()))
    for j, field in enumerate(CANONICAL_FIELDS):
        key = "c" if field == "c_eff_diagnostic" else field
        position = observations[key]
        index = np.flatnonzero(x == position)
        if len(index) != 1: raise ValueError("Preregistered observation is absent from material section grid")
        index = int(index[0])
        one_field = seven_field_one_d(active["nonlinear"])[..., j]
        three_field = pairs["linear"]["nonlinear"][..., j]
        observation_data[field] = {"material_x": position,
            "one_d_initial": float(one_field[0, index]), "one_d_final": float(one_field[-1, index]),
            "three_d_initial": float(three_field[0, index]), "three_d_final": float(three_field[-1, index]),
            "one_d_max_at_observation": float(abs(one_field[:, index]).max()),
            "three_d_max_at_observation": float(abs(three_field[:, index]).max()),
            "qualification": "3D effective contraction proxy, not identical generalized coordinate" if key == "c" else
                "one-dimensional planar assumption; 3D remainder diagnostic, not stability evidence" if key in ("v", "Phi", "psi") else "active planar coordinate"}
        if key in ("v", "Phi", "psi"):
            scale = active_displacement_scale if key == "v" else active_rotation_scale
            remainder, _ = old.sampled_norms(three_field, x, grid)
            inactive[field] = {"one_d_identically_zero_by_planar_subspace": bool(np.all(one_field == 0.)),
                "three_d_observed_max": remainder["absolute_max"],
                "physical_scale": scale, "scale_field": "w" if key == "v" else "theta",
                "remainder_over_active_physical_scale": remainder["absolute_max"]/scale if scale else None,
                "out_of_plane_stability_verified": False}
    for quantity, value in active.items():
        arrays["one_d_"+quantity+"_fields"] = seven_field_one_d(value)
        arrays["three_d_"+quantity+"_fields"] = pairs["linear"][quantity]
        arrays["three_d_PCHIP_"+quantity+"_fields"] = pairs["PCHIP"][quantity]
        for name, curve in curves.items(): arrays["three_d_"+name] = curve
        for name, curve in interpolation_curves.items(): arrays["interpolation_"+name] = curve
    for kind in ("linear", "nonlinear"):
        if kind+"_velocities" in reference:
            arrays["one_d_"+kind+"_velocities"] = seven_field_one_d(reference[kind+"_velocities"])
        if kind+"_translation_velocities" in pairs["linear"]:
            arrays["three_d_"+kind+"_translation_velocities"] = pairs["linear"][kind+"_translation_velocities"]
    parent = old.pilot.ROOT/config["parent_completed"]["bundle"]
    native_endpoints, energy = {}, {"status": "PARTIAL", "raw_native_data_unchanged": True,
        "independent_internal_StVK_energy": "NOT_RUN", "cases": {}, "one_d": {}}
    for kind in ("linear", "nonlinear"):
        h = load_arrays(cases/kind/"section_history.npz")
        ini = load_arrays(cases/kind/"initial_sections.npz")
        native_endpoints[kind] = {"initial_static_w": float(ini["fields"][20, 1]),
            "final_native_dynamic_w": float(h["fields"][-1, 20, 1]),
            "native_first_time": float(h["time"][0]), "native_final_time": float(h["time"][-1]),
            "actual_dynamic_output_frames": len(h["time"]),
            "zero_is_confirmed_static_not_native_dynamic": True}
        recovery = read_json(cases/kind/"recovery.json") if (cases/kind/"recovery.json").exists() else {}
        energy["cases"][kind] = {name: recovery.get(name) for name in (
            "native_initial_internal_energy", "native_dynamic_bookkeeping_initial_energy",
            "native_energy_static_to_dynamic_reference_jump_relative",
            "max_native_dynamic_bookkeeping_reference_relative_drift",
            "maximum_native_external_work_after_release", "maximum_native_damping_work_after_release")}
        energy["cases"][kind]["independent_kinetic_energy"] = _saved_kinetic_comparison(cases/kind)
        one_meta = read_json(bundle/("one_d_p64_"+kind+".json")) if (bundle/("one_d_p64_"+kind+".json")).exists() else {}
        energy["one_d"][kind] = one_meta.get("diagnostics", one_meta)
    result = {"version": VERSION, "status": "ILLUSTRATIVE_FULL_LINEAR_PERIOD_COMPLETE_WITH_QUALIFICATIONS",
        "T1": T1, "comparison_policy": config["full_period_comparison"],
        "physical_interval": [0., float(grid[-1])], "actual_native_endpoints": native_endpoints,
        "comparison_points": len(grid), "section_points": len(x), "seven_field_order": CANONICAL_FIELDS,
        "canonical_one_d_field_order": ONE_D_FIELD_ORDER,
        "recovered_three_d_field_order": CANONICAL_FIELDS,
        "observations": observation_data, "inactive_one_d_fields": inactive,
        "three_d_quantities": table, "model_comparisons": model,
        "interpolation_diagnostic_differences": interpolation_table,
        "temporal_transfer": metadata, "quarter_period_prefix": _prefix_observations(cases, parent/"cases"),
        "energy": energy, "one_d_spatial_status": summary.get("p48_spatial_control", {}).get("status", "NOT_RUN"),
        "quarter_period_robustness_not_extended_to_full_period": True,
        "no_phase_amplitude_frequency_time_fitting": True,
        "nonlinear_periodic_orbit_assumed_or_found": False,
        "p64_exact_continuum_truth": False, "experimental_validation": False,
        "general_seven_field_nonlinear_validation": False,
        "coupled_rod_joint_validation": False, "scientific_calls": 0}
    result["full_period_interpolation_uncertainty"] = {
        "primary_interpolation": "linear", "diagnostic_interpolation": "PCHIP",
        "scope": "fixed401-point full-period grid; native timestamps retained separately",
        "status": "OBSERVED_DIAGNOSTIC_NOT_FULL_PERIOD_CONVERGENCE_CERTIFICATION",
        "quarter_period_robustness_not_extended": True,
        "continuum_error_bound_claimed": False, "w": {}}
    for quantity in ("correction", "evolution"):
        uncertainty = interpolation_table[quantity]["w"].copy()
        scale = common_scales[quantity]["w"]
        model_difference = model[quantity]["w"]["absolute_max"]
        uncertainty.update(fixed_full_period_characteristic_scale=scale,
            relative_max_on_full_period_signal_scale=uncertainty["absolute_max"]/scale if scale else None,
            relative_L2_on_full_period_signal_scale=uncertainty["max_time_L2"]/scale if scale else None,
            full_period_model_discrepancy=model_difference,
            ratio_to_full_period_model_discrepancy=uncertainty["absolute_max"]/model_difference if model_difference else None,
            scientific_acceptance_threshold_assigned=False)
        result["full_period_interpolation_uncertainty"]["w"][quantity] = uncertainty
    result["nonlinear_correction_midspan"] = {}
    for source in ("one_d", "three_d"):
        delta = arrays[source+"_correction_fields"][:, 20, 1]
        evol = arrays[source+"_evolution_fields"][:, 20, 1]
        result["nonlinear_correction_midspan"][source] = {
            "signed_initial": float(delta[0]), "signed_final": float(delta[-1]),
            "signed_final_evolution": float(evol[-1]), "maximum_evolution_abs": float(abs(evol).max())}
    np.savez_compressed(bundle/"full_period_comparison.npz", **arrays)
    write_json(bundle/"full_period_comparison.json", result)
    _full_csv(bundle, arrays, observation_data, cases)
    return result


def _full_csv(bundle, data, observations, cases):
    with (bundle/"full_period_observations.csv").open("w", encoding="utf8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("physical_time", "tau", "source", "kind", "field", "material_x", "value", "provenance"))
        T1 = data["times"][-1]
        for source in ("one_d", "three_d"):
            for kind in ("linear", "nonlinear", "correction", "evolution"):
                for j, field in enumerate(CANONICAL_FIELDS):
                    position = observations[field]["material_x"]
                    index = int(np.flatnonzero(data["x"] == position)[0])
                    values = data[source+"_"+kind+"_fields"][:, index, j]
                    for t, value in zip(data["times"], values):
                        writer.writerow((t, t/T1, source, kind, field, position, value,
                            "static_preload_at_zero" if t == 0 else "saved_dense_or_exact_1D" if source == "one_d" else "postprocessing_linear_interpolation"))
    with (bundle/"full_period_native_observations.csv").open("w", encoding="utf8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("native_physical_time", "kind", "field", "material_x", "value", "provenance"))
        for kind in ("linear", "nonlinear"):
            h = load_arrays(cases/kind/"section_history.npz")
            ini = load_arrays(cases/kind/"initial_sections.npz")
            for j, field in enumerate(CANONICAL_FIELDS):
                position = observations[field]["material_x"]
                index = int(np.flatnonzero(h["x"] == position)[0])
                writer.writerow((0., kind, field, position, ini["fields"][index, j], "confirmed_STATIC_preload_not_native_DYNAMIC_frame"))
                for t, value in zip(h["time"], h["fields"][:, index, j]):
                    writer.writerow((t, kind, field, position, value, "actual_native_DYNAMIC_sample"))


def plot_full_period(bundle):
    """Three full-period illustrations, entirely from the saved comparison."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bundle = Path(bundle)
    data = load_arrays(bundle/"full_period_comparison.npz")
    result = read_json(bundle/"full_period_comparison.json")
    tau = data["times"]/result["T1"]
    figures, paths = [], []
    plt.rcParams.update({"font.size": 10, "pdf.fonttype": 42})
    fig, axes = plt.subplots(2, 4, figsize=(13, 6.5), constrained_layout=True)
    wscale = max(abs(data["one_d_nonlinear_fields"][..., 1]).max(), abs(data["three_d_nonlinear_fields"][..., 1]).max())
    rscale = max(abs(data["one_d_nonlinear_fields"][..., 5]).max(), abs(data["three_d_nonlinear_fields"][..., 5]).max())
    for j, field in enumerate(CANONICAL_FIELDS):
        ax = axes.flat[j]
        position = result["observations"][field]["material_x"]
        index = int(np.flatnonzero(data["x"] == position)[0])
        ax.plot(tau, data["one_d_nonlinear_fields"][:, index, j], label="1D nonlinear")
        ax.plot(tau, data["three_d_nonlinear_fields"][:, index, j], "--", label="3D section recovery, nonlinear")
        title = {"Phi": r"$\Phi$", "psi": r"$\psi$", "theta": r"$\theta$",
            "c_eff_diagnostic": "c\n3D effective contraction proxy"}.get(field, field)
        ax.set(title=(f"c (x/L={position:g})\n3D effective contraction proxy"
            if field == "c_eff_diagnostic" else f"{title} (x/L={position:g})"),
            xlabel="t/T₁", xlim=(0, 1))
        if field in ("v", "Phi", "psi"):
            scale = wscale if field == "v" else rscale
            ax.set_ylim(-1.08*scale, 1.08*scale)
            row = result["inactive_one_d_fields"][field]
            ax.text(.03, .97, "1D zero by planar assumption\n3D max/"+row["scale_field"]+f" scale = {row['remainder_over_active_physical_scale']:.2e}",
                transform=ax.transAxes, va="top", fontsize=8)
        ax.grid(alpha=.25)
    axes.flat[7].axis("off")
    axes.flat[7].text(0., .9, "Four active planar fields: u, w, θ, c.\n\n3D rotations: section polar recovery.\n3D c: effective contraction proxy.\n\nThree 1D zeros are the planar subspace.\n3D remainders do not establish stability.\n\nT₁ is the first linear period.\nFull-period comparison is illustrative.", va="top", fontsize=9)
    axes.flat[7].legend(*axes.flat[0].get_legend_handles_labels(), loc="lower left", fontsize=8)
    figures.append((fig, "full_period_seven_fields"))
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), constrained_layout=True)
    for source, title in (("one_d", "1D"), ("three_d", "3D")):
        for kind, style in (("linear", "--"), ("nonlinear", "-")):
            axes[0].plot(tau, data[source+"_"+kind+"_fields"][:, 20, 1], style, label=title+" "+kind)
        axes[1].plot(tau, data[source+"_correction_fields"][:, 20, 1]*1e6, label=title)
        axes[2].plot(tau, data[source+"_evolution_fields"][:, 20, 1]*1e6, label=title)
    for ax, label in zip(axes, ("Midspan w", "Midspan NL−L correction × 10⁶", "Midspan evolving correction × 10⁶")):
        ax.set(xlabel="t/T₁", ylabel=label, xlim=(0, 1)); ax.grid(alpha=.25); ax.legend(fontsize=8)
    figures.append((fig, "full_period_free_motion_and_correction"))
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    for ax, (field, j) in zip(axes.flat, (("w", 1), ("u", 0), ("theta", 5), ("c_eff_diagnostic", 6))):
        for fraction, color in zip((0., .25, .5, .75, 1.), plt.cm.viridis(np.linspace(0, 1, 5))):
            index = int(round(fraction*(len(tau)-1)))
            ax.plot(data["x"], data["one_d_nonlinear_fields"][index, :, j], color=color, label=f"t/T₁={fraction:g}")
            ax.plot(data["x"], data["three_d_nonlinear_fields"][index, :, j], "--", color=color)
        ax.set(xlabel="Original material x/L", ylabel="c; 3D effective contraction proxy" if j == 6 else field)
        ax.grid(alpha=.25)
    axes.flat[0].legend(fontsize=8)
    fig.suptitle("Nonlinear spatial profiles: 1D solid; 3D dashed (fixed physical times)")
    figures.append((fig, "full_period_representative_profiles"))
    for fig, name in figures:
        for suffix in ("png", "pdf"):
            path = bundle/(name+"."+suffix)
            fig.savefig(path, dpi=180, metadata={"CreationDate": None, "ModDate": None} if suffix == "pdf" else None)
            paths.append(path)
        plt.close(fig)
    return paths
