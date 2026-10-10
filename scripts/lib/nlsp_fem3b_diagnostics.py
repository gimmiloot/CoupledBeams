"""FEM-3B saved-data comparison, recovery diagnostics and cached figures.

This module never executes native solvers, time integration, eigenanalysis,
equilibria or symbolic derivation. The existing physical-field and C3D10
recovery helpers are reused. Initial STATIC data are explicitly distinguished
from positive-time native DYNAMIC observations.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import time

import numpy as np

from scripts.analysis import pilot_nlsp_nonlinear_dynamic_3d_fem as pilot

VERSION = "fem3b-saved-transient-diagnostics-v1"
FIELDS = {"u": (0, 0), "w": (1, 1), "theta": (2, 5),
          "c_eff_diagnostic": (3, 6)}


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf8")


def load_arrays(path):
    with np.load(path, allow_pickle=False) as data:
        return {name: data[name].copy() for name in data.files}


def _exact_indices(available, requested):
    available, requested = np.asarray(available), np.asarray(requested)
    indices = np.searchsorted(available, requested)
    if np.any(indices >= len(available)) or not np.array_equal(available[indices], requested):
        raise ValueError("Saved 1D evaluation must use precisely requested physical timestamps")
    return indices


def _native_data(bundle):
    bundle = Path(bundle)
    histories = {kind: load_arrays(bundle / "cases" / kind / "section_history.npz")
                 for kind in ("linear", "nonlinear")}
    initial = {kind: load_arrays(bundle / "cases" / kind / "initial_sections.npz")
               for kind in histories}
    x = histories["linear"]["x"]
    if (not np.array_equal(x, histories["nonlinear"]["x"])
            or any(not np.array_equal(x, initial[k]["x"]) for k in histories)):
        raise ValueError("Material section coordinates differ between saved cases")
    if not all(bool(initial[k]["not_a_native_dynamic_zero_frame"]) for k in histories):
        raise ValueError("Initial sections must retain their STATIC provenance")
    return histories, initial


def _pairing(bundle, item):
    from scripts.lib.nlsp_fem3b_continuation import pair_native_histories
    histories, initial = _native_data(bundle)
    first, second = (histories[k] for k in ("linear", "nonlinear"))
    tl, tn = first["time"], second["time"]
    grid = None
    if not np.array_equal(tl, tn):
        T = 2 * np.pi / item["config"]["omega1"]
        step = T / 2000
        low, high = max(tl[0], tn[0]), min(tl[-1], tn[-1])
        start, stop = int(np.ceil(low / step)), int(np.floor(high / step))
        grid = np.arange(start, stop + 1) * step
        grid = grid[(grid >= low) & (grid <= high)]
        if len(grid) < 2:
            raise ValueError("Predetermined interpolation grid has insufficient actual overlap")
    paired = pair_native_histories(tl, first["fields"], tn, second["fields"],
                                   common_grid=grid)
    return histories, initial, paired


def actual_pairing_times(bundle, item):
    """Collect exact timestamps needed by the existing saved 1D evaluator.

    Native schedules are retained for separate L/NL comparisons. The approved
    fallback grid is added only when schedules differ. Historical 1D times
    are included solely for the read-only prefix reproduction diagnostic.
    """
    histories, _, paired = _pairing(bundle, item)
    parent = pilot.ROOT / item["parent_completed"]["bundle"]
    old_times = load_arrays(parent / "one_d_linear.npz")["times"]
    horizon = min(histories[k]["time"][-1] for k in histories)
    return np.unique(np.r_[0., histories["linear"]["time"], histories["nonlinear"]["time"],
                            paired["time"], old_times[old_times <= horizon]])


def sampled_norms(values, x, times):
    values, x, times = np.asarray(values), np.asarray(x), np.asarray(times)
    if (values.shape != (len(times), len(x)) or not np.isfinite(values).all()
            or np.any(np.diff(x) <= 0) or np.any(np.diff(times) < 0)):
        raise ValueError("Invalid sampled physical profile history")
    L2 = np.sqrt(np.trapezoid(values**2, x, axis=1))
    maxima = np.max(abs(values), axis=1)
    index = np.unravel_index(np.argmax(abs(values)), values.shape)
    return {"absolute_max": float(maxima.max()), "max_time_L2": float(L2.max()),
        "signed_at_max": float(values[index]), "x_at_max": float(x[index[1]]),
        "time_at_max": float(times[index[0]]), "midspan_max_abs": float(abs(values[:, len(x)//2]).max()),
        "midspan_final": float(values[-1, len(x)//2]), "sampled_maxima_only": True}, {
        "L2": L2, "max": maxima, "cumulative_L2": np.maximum.accumulate(L2),
        "cumulative_max": np.maximum.accumulate(maxima)}


def difference_metrics(first, second, x, times, scale=None, fixed_scale=None):
    result = pilot.field_difference(first, second, x, scale=scale)
    result["time_at_max"] = float(np.asarray(times)[result["time_index_at_max"]])
    result["fixed_physical_scale"] = fixed_scale
    result["relative_max_on_fixed_physical_scale"] = result["absolute_max"] / fixed_scale if fixed_scale else None
    result["relative_L2_on_fixed_physical_scale"] = result["max_time_L2"] / fixed_scale if fixed_scale else None
    return result


def evolving_correction(correction, initial):
    correction, initial = np.asarray(correction), np.asarray(initial)
    if correction.shape[1:] != initial.shape:
        raise ValueError("Initial correction and transient profiles disagree")
    return correction - initial[None, ...]


def signal_resolution(signal, differences):
    """Compare observed scales without manufacturing a validation threshold."""
    signal = float(signal)
    measured = {key: float(value) for key, value in differences.items()
                if value is not None}
    if signal < 0 or not np.isfinite(signal) or any(v < 0 or not np.isfinite(v) for v in measured.values()):
        raise ValueError("Signals and diagnostic differences must be finite nonnegative")
    largest = max(measured.values(), default=0.)
    return {"signal": signal, "observed_diagnostic_differences": measured,
        "largest_observed_difference": largest,
        "signal_to_largest_observed_difference": signal/largest if largest else None,
        "status": "SIGNAL_EXCEEDS_OBSERVED_DIAGNOSTIC_DIFFERENCES" if signal > largest else
                  "DYNAMIC_NONLINEAR_EVOLUTION_UNRESOLVED",
        "physical_validation_threshold": False, "continuum_error_bound_claimed": False,
        "independent_3D_temporal_spatial_certification": False}


def _write_csv(path, labels, columns):
    with Path(path).open("w", newline="", encoding="utf8") as stream:
        writer = csv.writer(stream)
        writer.writerow(labels)
        writer.writerows(zip(*columns))


def one_d_qualifications(bundle, horizon):
    """Read saved degree controls and additive decomposition on fixed scales."""
    bundle=Path(bundle)
    decomp=load_arrays(bundle/"one_d_p64_decomposition.npz")
    mask=decomp["times"]<=horizon
    t,x=decomp["times"][mask],decomp["x"]
    terms={key:decomp[key][mask] for key in
        ("total_correction","initial_state_component","same_ic_nonlinear_component","evolving_correction")}
    identity=terms["total_correction"]-terms["initial_state_component"]-terms["same_ic_nonlinear_component"]
    decomposition={"identity_max_abs":float(abs(identity).max()),
        "qualification":"additive signed motion components, not percentages of an instantaneous correction or modal/energy fractions",
        "fixed_characteristic_scales":{},"fields":{}}
    for name,(i,_) in FIELDS.items():
        scale=float(abs(terms["total_correction"][:,:,i]).max())
        decomposition["fixed_characteristic_scales"][name]=scale
        row={}
        for key,values in terms.items():
            metrics,_=sampled_norms(values[:,:,i],x,t)
            metrics["max_over_fixed_total_correction_scale"]=metrics["absolute_max"]/scale if scale else None
            metrics["L2_over_fixed_total_correction_scale"]=metrics["max_time_L2"]/scale if scale else None
            metrics["signed_final_midspan_over_fixed_total_correction_scale"]=metrics["midspan_final"]/scale if scale else None
            row[key]=metrics
        decomposition["fields"][name]=row
    result={"decomposition":decomposition,"spatial_degree_control":{"status":"NOT_RUN"},
        "temporal_control":{"status":"NOT_CERTIFIED","reason":"one unchanged tight Radau policy; no new temporal refinement authorized or executed"}}
    if (bundle/"one_d_p48_nonlinear.npz").exists() and (bundle/"one_d_p48_linear.npz").exists():
        trajectories={p:{kind:load_arrays(bundle/("one_d_p"+str(p)+"_"+kind+".npz"))
            for kind in ("linear","nonlinear")} for p in (48,64)}
        shared=np.intersect1d(trajectories[48]["nonlinear"]["times"],trajectories[64]["nonlinear"]["times"])
        shared=shared[shared<=horizon]
        if len(shared)<2:raise ValueError("Saved p48/p64 controls have no usable exact common grid")
        states={p:{} for p in trajectories}
        initial={p:{} for p in trajectories}
        for p,data in trajectories.items():
            for kind,d in data.items():
                indices=_exact_indices(d["times"],shared)
                states[p][kind]=d["fields"][indices]
                states[p][kind+"_velocity"]=d["physical_velocities"][indices]
                initial[p][kind]=d["fields"][0]
            states[p]["correction"]=states[p]["nonlinear"]-states[p]["linear"]
            states[p]["evolution"]=states[p]["correction"]-(initial[p]["nonlinear"]-initial[p]["linear"])[None,:,:]
        control={"status":"COMPLETED_DIAGNOSTIC","shared_samples":len(shared),
            "actual_interval":[float(shared[0]),float(shared[-1])],"fields":{},"correction":{},"evolution":{},
            "qualification":"observed physical differences of saved spaces; p64 is not continuum truth"}
        for name,(i,_) in FIELDS.items():
            control["fields"][name]=difference_metrics(states[48]["nonlinear"][:,:,i],states[64]["nonlinear"][:,:,i],x,shared)
            control["fields"][name+"_t"]=difference_metrics(states[48]["nonlinear_velocity"][:,:,i],states[64]["nonlinear_velocity"][:,:,i],x,shared)
            for key in ("correction","evolution"):
                control[key][name]=difference_metrics(states[48][key][:,:,i],states[64][key][:,:,i],x,shared)
        result["spatial_degree_control"]=control
    for kind in ("linear","nonlinear"):
        path=bundle/("one_d_p64_"+kind+".npz")
        data=load_arrays(path); m=data["times"]<=horizon
        result[kind+"_mechanical_energy"]={"initial":float(data["energy"][0]),
            "max_relative_drift_on_selected_horizon":float(abs(data["energy_relative_drift"][m]).max()),
            "removed_gravity_potential_included":False}
    write_json(bundle/"one_d_comparison_qualifications.json",result)
    return result


def _prefix_reproduction(bundle, item, one_d):
    parent = pilot.ROOT / item["parent_completed"]["bundle"]
    result = {"scope": "exact shared native times only; no temporal interpolation",
        "parent": item["parent_completed"], "three_d": {}, "one_d": {}}
    source_coordinates_path = pilot.ROOT/item["config"]["source_static"]["bundle"]/"one_d_p64.npz"
    source_coordinates = load_arrays(source_coordinates_path)
    all_exact = True
    for kind in ("linear", "nonlinear"):
        old = load_arrays(parent / "cases" / kind / "section_history.npz")
        new = load_arrays(bundle / "cases" / kind / "section_history.npz")
        shared, old_index, new_index = np.intersect1d(old["time"], new["time"], return_indices=True)
        if len(shared) < 2:
            raise ValueError("No usable exact native prefix overlap with FEM-3AR")
        fields, velocities = {}, {}
        for name, (_, j) in FIELDS.items():
            fields[name] = difference_metrics(new["fields"][new_index, :, j],
                old["fields"][old_index, :, j], old["x"], shared)
        for j, name in enumerate(("u_t", "w_t", "v_t")):
            velocities[name] = difference_metrics(new["translation_velocities"][new_index, :, j],
                old["translation_velocities"][old_index, :, j], old["x"], shared)
        a, b = (load_arrays(z / "cases" / kind / "initial_sections.npz") for z in (parent, Path(bundle)))
        initial = float(np.max(abs(a["fields"] - b["fields"])))
        exact = initial == 0. and all(v["absolute_max"] == 0. for v in (*fields.values(), *velocities.values()))
        all_exact &= exact
        result["three_d"][kind] = {"shared_native_samples": len(shared), "shared_first_time": float(shared[0]),
            "shared_last_time": float(shared[-1]), "old_native_samples": len(old["time"]),
            "old_times_not_exactly_shared": np.setdiff1d(old["time"], shared).tolist(),
            "initial_static_fields_max_abs": initial, "fields": fields, "translation_velocities": velocities,
            "exact_saved_profile_reproduction": exact, "unmatched_endpoint_not_interpolated": True}
        old_one = load_arrays(parent / ("one_d_" + kind + ".npz"))
        indices = _exact_indices(one_d["times"], old_one["times"])
        fields_1d = {name: difference_metrics(one_d[kind + "_fields"][indices, :, i],
            old_one["fields"][:, :, i], old_one["x"], old_one["times"])
            for name, (i, _) in FIELDS.items()}
        result["one_d"][kind] = {"shared_samples": len(indices), "fields": fields_1d,
            "same_physical_initial_state": bool(np.array_equal(one_d["initial_"+kind+"_fields"], old_one["fields"][0])),
            "initial_reconstructed_fields_max_abs_difference":float(np.max(abs(one_d["initial_"+kind+"_fields"]-old_one["fields"][0]))),
            "physical_bitwise_flag_definition":"bitwise equality of reconstructed float64 profiles only; not a test that target functions or frozen initial coordinates were changed",
            "new_saved_q0_equals_source_bitwise":bool(np.array_equal(load_arrays(Path(bundle)/("one_d_p64_"+kind+".npz"))["q"][0],source_coordinates["q_"+kind])),
            "old_saved_q0_equals_source_bitwise":bool(np.array_equal(old_one["q"][0],source_coordinates["q_"+kind])),
            "source_coordinates_path":str(source_coordinates_path.relative_to(pilot.ROOT)),
            "source_coordinates_sha256":pilot.sha(source_coordinates_path),
            "qualification": "new dense/exact-time evaluation of the same IVP; longer t_bound can alter final adaptive steps. A few-ULP reconstructed-profile difference can arise from eigenspace zero-time arithmetic and is reported separately from exact frozen-q0 identity."}
    result["three_d_exact_saved_reproduction"] = all_exact
    result["status"] = "PASS" if all_exact else "PARTIAL"
    write_json(Path(bundle) / "prefix_reproduction.json", result)
    return result


def _recovery_sensitivity(bundle, item, histories):
    """Limited paired initial/final41/81 recovery from saved native fields."""
    fem2, fem1 = pilot.base, pilot.base.fem1
    mesh_path = pilot.ROOT / item["config"]["source_fem1"]["bundle"] / "meshes/medium/rod.inp"
    mesh = fem1.single.read_gmsh_inp_mesh_data(mesh_path)
    quad = fem1.quadrature_arrays(mesh, item["config"]["material"]["rho"])
    x = np.linspace(0., 1., 801)
    recovered = {}; arrays = {"x": x}
    for kind in ("linear", "nonlinear"):
        case = Path(bundle) / "cases" / kind
        static_frame = sorted((case / "frames").glob("step1_inc*.npz"))[-1]
        final_meta = read_json(case / "frame_metadata.json")[-1]
        final_frame = Path(bundle) / final_meta["frame"]
        recovered[kind] = {}
        for label, path in (("initial", static_frame), ("final", final_frame)):
            with np.load(path, allow_pickle=False) as native:
                U = native["U"]
                if not np.array_equal(native["node_ids"], quad["node_ids"]):
                    raise ValueError("Saved nodal recovery order mismatch")
                Uq = fem1.nlsp_evaluate_tet10_displacements(U, quad["conn"], quad["N"])
                recovered[kind][label] = {}
                for count in (41, 81):
                    profile = fem2.fem2_recover_reference_samples(quad["xyz"], Uq, quad["weights"], 1., .1, .2, count)
                    fields = fem2.fem2_static_sample(profile, x)
                    recovered[kind][label][count] = fields
                    arrays[kind + "_" + label + "_" + str(count)] = fields
    result = {"policy": "paired same-material initial/final recovery41_vs81",
        "qualification": "limited initial/final observations, not an entire-history error bound",
        "initial": {}, "final": {}, "evolving_part": {}}
    for label in ("initial", "final"):
        for count in (41, 81):
            arrays["correction_" + label + "_" + str(count)] = recovered["nonlinear"][label][count] - recovered["linear"][label][count]
        result[label] = {name: fem2.fem2_curve_difference(x,
            arrays["correction_"+label+"_41"][:, j], arrays["correction_"+label+"_81"][:, j])
            for name, (_, j) in FIELDS.items()}
    for count in (41, 81):
        arrays["evolution_" + str(count)] = arrays["correction_final_"+str(count)] - arrays["correction_initial_"+str(count)]
    result["evolving_part"] = {name: fem2.fem2_curve_difference(x,
        arrays["evolution_41"][:, j], arrays["evolution_81"][:, j]) for name, (_, j) in FIELDS.items()}
    np.savez_compressed(Path(bundle) / "paired_recovery_sensitivity.npz", **arrays)
    write_json(Path(bundle) / "paired_recovery_sensitivity.json", result)
    return mesh, quad, result


def independent_kinetic_history(bundle, kind, quad):
    """Integrate actual nodal velocities with existing consistent mass helper."""
    case = Path(bundle) / "cases" / kind
    frames = read_json(case / "frame_metadata.json")
    records = read_json(case / "energy.json")["records"]
    native = {r["increment"]: r for r in records if r["step"] == 2}
    times, independent, printed = [], [], []
    for meta in frames:
        with np.load(Path(bundle) / meta["frame"], allow_pickle=False) as data:
            if not np.array_equal(data["node_ids"], quad["node_ids"]):
                raise ValueError("Saved velocity and mesh node ordering mismatch")
            Vq = pilot.base.fem1.nlsp_evaluate_tet10_displacements(data["VELO"], quad["conn"], quad["N"])
            kinetic = float(.5 * np.sum(quad["weights"] * np.sum(Vq*Vq, axis=2)))
        times.append(meta["time"]); independent.append(kinetic)
        printed.append(native[meta["increment"]]["kinetic_energy"])
    times, independent, printed = map(np.asarray, (times, independent, printed))
    difference = independent - printed
    result = {"status": "COMPLETED_DIAGNOSTIC", "definition": "0.5*integral_reference_rho*|N*actual_nodal_V|^2",
        "quadrature": "existing affine C3D10 positive14point degree5 consistent mass integral",
        "first_independent_K": float(independent[0]), "first_native_K": float(printed[0]),
        "final_independent_K": float(independent[-1]), "final_native_K": float(printed[-1]),
        "maximum_absolute_difference": float(abs(difference).max()),
        "maximum_relative_difference_fixed_K_scale": float(abs(difference).max()/max(independent.max(), printed.max())),
        "internal_StVK_reconstruction": "NOT_RUN", "internal_reason": "no existing checked full-element energy evaluator; no new3D implementation",
        "qualification": "rounded actual velocities and native energy quadrature differ; independent K does not repair native internal-energy bookkeeping"}
    np.savez_compressed(case / "independent_kinetic.npz", time=times, independent_K=independent, native_K=printed, difference=difference)
    _write_csv(case / "independent_kinetic.csv", ["time", "independent_K", "native_K", "difference"], [times, independent, printed, difference])
    write_json(case / "independent_kinetic.json", result)
    return result


def complete_comparison(bundle, item, summary):
    """Compare finished native jobs with pre-evaluated saved1D states only."""
    bundle = Path(bundle); started = time.perf_counter()
    if any(summary["cases"].get(k, {}).get("status") != "PASS" for k in ("linear", "nonlinear")):
        raise ValueError("Both actual3D execution gates must pass before comparison")
    histories, initial, paired = _pairing(bundle, item)
    one_d = load_arrays(bundle / "comparison_one_d.npz")
    x = histories["linear"]["x"]
    if not np.array_equal(one_d["x"], x):
        raise ValueError("1D/3D physical comparison coordinates disagree")
    times = paired["time"]; indices = _exact_indices(one_d["times"], times)
    a, b = one_d["linear_fields"][indices], one_d["nonlinear_fields"][indices]
    delta1 = b-a; delta3 = paired["correction"]
    initial1 = one_d["initial_nonlinear_fields"]-one_d["initial_linear_fields"]
    initial3 = initial["nonlinear"]["fields"]-initial["linear"]["fields"]
    evolution1, evolution3 = evolving_correction(delta1, initial1), evolving_correction(delta3, initial3)
    arrays = {"time": times, "x": x, "linear_1D": a, "nonlinear_1D": b,
        "linear_3D": paired["linear"], "nonlinear_3D": paired["nonlinear"],
        "one_d_correction": delta1, "three_d_correction": delta3,
        "one_d_initial_correction": initial1, "three_d_initial_correction": initial3,
        "one_d_evolution": evolution1, "three_d_evolution": evolution3,
        "paired_sample_origin": np.asarray(paired["sample_origin"])}
    comparison = {"version": VERSION, "field_maps": FIELDS, "pairing": {key: value for key, value in paired.items() if not isinstance(value, np.ndarray)},
        "scope": "one medium C3D10 mesh and original single timestep level; no independent3D convergence certification",
        "fixed_scale_policy": "one full selected-horizon characteristic scale per pair, never instantaneous zeros",
        "initial_STATIC_observations_are_not_native_DYNAMIC_zero_frames": True,
        "absolute_physical_differences": {}, "nonlinear_corrections": {}, "evolving_corrections": {},
        "evolution_norms": {}, "phase_amplitude_time_material_alignment": False,
        "c_eff_not_identical_to_MH_coordinate": True}
    for kind in ("linear", "nonlinear"):
        native = histories[kind]; match = _exact_indices(one_d["times"], native["time"])
        comparison["absolute_physical_differences"][kind] = {name: difference_metrics(
            one_d[kind+"_fields"][match, :, i], native["fields"][:, :, j], x, native["time"],fixed_scale=.1 if i<2 else 1.)
            for name, (i, j) in FIELDS.items()}
    for name, (i, j) in FIELDS.items():
        comparison["nonlinear_corrections"][name] = difference_metrics(delta1[:, :, i], delta3[:, :, j], x, times,fixed_scale=.1 if i<2 else 1.)
        comparison["evolving_corrections"][name] = difference_metrics(evolution1[:, :, i], evolution3[:, :, j], x, times,fixed_scale=.1 if i<2 else 1.)
        for model, ev, first, index in (("1D", evolution1, initial1, i), ("3D", evolution3, initial3, j)):
            values = np.vstack((np.zeros(len(x)), ev[:, :, index]))
            metrics, history = sampled_norms(values, x, np.r_[0., times])
            scale = float(abs(first[:, index]).max())
            metrics["initial_correction_scale"] = scale
            metrics["evolution_over_initial_correction_scale"] = metrics["absolute_max"]/scale if scale else None
            comparison["evolution_norms"][model+"_"+name] = metrics
            for key, value in history.items(): arrays[model+"_"+name+"_evolution_"+key] = value
    if paired["interpolation_used"]:
        arrays.update({key: value for key, value in paired.items() if isinstance(value, np.ndarray)})
    np.savez_compressed(bundle / "dynamic_comparison.npz", **arrays)
    T = 2*np.pi/item["config"]["omega1"]
    _write_csv(bundle / "midspan_response.csv",
        ["time", "tau", "sample_origin", "1D_linear_w", "1D_nonlinear_w", "3D_linear_w", "3D_nonlinear_w", "1D_delta_w", "3D_delta_w", "1D_evolving_delta_w", "3D_evolving_delta_w"],
        [times, times/T, [paired["sample_origin"]]*len(times), a[:,20,1], b[:,20,1], paired["linear"][:,20,1], paired["nonlinear"][:,20,1], delta1[:,20,1], delta3[:,20,1], evolution1[:,20,1], evolution3[:,20,1]])
    rows = []
    for category in ("absolute_physical_differences", "nonlinear_corrections", "evolving_corrections"):
        source = comparison[category]
        if category == "absolute_physical_differences": source = {k+"_"+f:m for k,v in source.items() for f,m in v.items()}
        for field, metric in source.items(): rows.append([category, field, metric["absolute_max"], metric["max_time_L2"], metric["characteristic_scale"], metric["relative_max"], metric["relative_max_L2"], metric["signed_at_max"], metric["x_at_max"], metric["time_at_max"]])
    with (bundle / "comparison_metrics.csv").open("w", newline="", encoding="utf8") as stream:
        writer=csv.writer(stream); writer.writerow(["comparison", "field", "absolute_max", "max_time_L2", "fixed_characteristic_scale", "relative_max", "relative_L2", "signed_at_max", "x_at_max", "time_at_max"]); writer.writerows(rows)
    prefix = _prefix_reproduction(bundle, item, one_d)
    mesh, quad, sensitivity = _recovery_sensitivity(bundle, item, histories)
    kinetic = {kind: independent_kinetic_history(bundle, kind, quad) for kind in ("linear", "nonlinear")}
    one_d_quality=one_d_qualifications(bundle,float(times[-1]))
    output = max(summary["cases"][k]["max_DAT_FRD_displacement_difference"] for k in ("linear", "nonlinear"))
    recovery = sensitivity["evolving_part"]["w"]["absolute_max"]
    interpolation = float(np.max(abs(paired["correction_PCHIP"][:,:,1]-delta3[:,:,1]))) if paired["interpolation_used"] else 0.
    comparison["pairing"]["w_interpolation_correction_difference_max_abs"]=interpolation
    signal = float(np.max(abs(evolution3[:, :, 1])))
    resolution = signal_resolution(signal, {"DAT_FRD_displacement": output,
        "paired_initial_final_41_81_evolving_correction": recovery,
        "time_pairing_linear_PCHIP_correction": interpolation})
    comparison["signal_resolution"] = resolution
    comparison["one_d_qualifications"]=one_d_quality
    comparison["uncertainty_separation"]={
        "1D_spatial":one_d_quality["spatial_degree_control"]["status"],
        "1D_temporal":"one tight level; energy drift observed, independent time convergence NOT_CERTIFIED",
        "3D_temporal":"one original timestep level; independent convergence NOT_CERTIFIED",
        "3D_spatial_mesh":"one medium mesh; independent dynamic mesh convergence NOT_CERTIFIED",
        "3D_recovery_output":"DAT_FRD and paired initial/final41_81 measured; not strict error bounds"}
    comparison["energy"] = {"status": "PARTIAL", "native_STATIC_DYNAMIC_reference_discontinuity_retained": True,
        "independent_kinetic": kinetic, "native": {kind: {key: summary["cases"][kind].get(key)
            for key in ("native_initial_internal_energy", "native_dynamic_bookkeeping_initial_energy", "native_energy_static_to_dynamic_reference_jump_relative", "max_native_dynamic_bookkeeping_reference_relative_drift", "maximum_native_external_work_after_release", "maximum_native_damping_work_after_release")}
            for kind in ("linear", "nonlinear")}, "independent_internal_StVK": "NOT_RUN",
        "no_energy_offset_correction": True}
    comparison["native_midspan"] = {kind: {"initial_static_w": float(initial[kind]["fields"][20,1]),
        "final_dynamic_w": float(histories[kind]["fields"][-1,20,1]), "first_native_time": float(histories[kind]["time"][0]),
        "final_native_time": float(histories[kind]["time"][-1]), "actual_native_frames": len(histories[kind]["time"])} for kind in ("linear", "nonlinear")}
    comparison["postprocessing_seconds"] = time.perf_counter()-started
    write_json(bundle / "dynamic_comparison.json", comparison)
    summary["comparison"] = comparison
    statuses = summary["long_horizon_statuses"]
    statuses["NLSP_FEM3B_PREFIX_REPRODUCTION"] = prefix["status"]
    statuses["NLSP_FEM3B_RESPONSE_COMPARISON"] = "PASS"
    statuses["NLSP_FEM3B_ENERGY_DIAGNOSTICS"] = "PARTIAL"
    statuses["NLSP_FEM3B_DYNAMIC_NONLINEAR_SIGNAL"] = "PASS" if resolution["status"] == "SIGNAL_EXCEEDS_OBSERVED_DIAGNOSTIC_DIFFERENCES" else "PARTIAL"
    summary["overall"] = "PILOT_COMPLETE_WITH_QUALIFICATIONS"
    return comparison


def plot_bundle(bundle):
    """Render exactly three figures from validated saved arrays, zero science."""
    from scripts.lib import nlsp_fem3b_continuation as continuation
    bundle = Path(bundle)
    summary = continuation.validate_cache(bundle)
    if not (bundle / "dynamic_comparison.npz").exists():
        return {"figures": 0, "new_scientific_calls": 0}
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = load_arrays(bundle / "dynamic_comparison.npz")
    metrics = read_json(bundle / "dynamic_comparison.json")
    T = summary["preflight"]["T1"]; tau = data["time"]/T
    horizon = summary["selected_horizon_T1"]
    figures = []; destination = bundle / "figures"; destination.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "savefig.dpi": 240, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5), constrained_layout=True)
    for kind, style in (("linear", "--"), ("nonlinear", "-")):
        for ax, (field, position, label) in zip(axes, (("w", 20, "w(L/2)"), ("theta", 10, "theta(L/4)"), ("u", 10, "u(L/4)"))):
            i, j = FIELDS[field]
            ax.plot(tau, data[kind+"_1D"][:,position,i], style, label="1D "+kind, linewidth=1.3)
            ax.plot(tau, data[kind+"_3D"][:,position,j], style, label="3D "+kind, linewidth=1.1)
            ax.set(xlabel="t/T1", ylabel=label, xlim=(0,horizon)); ax.grid(alpha=.25)
    axes[0].legend(fontsize=8); figures.append((fig,"linear_nonlinear_free_motion"))
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.5), constrained_layout=True)
    for prefix, label in (("one_d", "1D"), ("three_d", "3D")):
        initial = data[prefix+"_initial_correction"][20,1]
        axes[0].plot(np.r_[0.,tau], 1e6*np.r_[initial,data[prefix+"_correction"][:,20,1]], label=label)
        axes[1].plot(np.r_[0.,tau], 1e6*np.r_[0.,data[prefix+"_evolution"][:,20,1]], label=label)
        axes[2].plot(data["x"], 1e6*data[prefix+"_evolution"][-1,:,1], label=label)
    axes[0].set(xlabel="t/T1",ylabel="Midspan Delta w (1e-6)",xlim=(0,horizon))
    axes[1].set(xlabel="t/T1",ylabel="Midspan evolving Delta w (1e-6)",xlim=(0,horizon))
    axes[2].set(xlabel="Material x/L",ylabel="Final evolving Delta w (1e-6)")
    for ax in axes: ax.grid(alpha=.25); ax.legend(fontsize=8)
    if metrics["pairing"]["interpolation_used"]: axes[0].set_title("3D paired values interpolated")
    figures.append((fig,"nonlinear_correction_and_evolution"))
    fig, axes = plt.subplots(1, 2, figsize=(9,3.5), constrained_layout=True)
    decomposition = load_arrays(bundle / "one_d_p64_decomposition.npz")
    mask = decomposition["times"] <= horizon*T
    for key,label in (("total_correction","Total NL-L"),("initial_state_component","Linear initial-state difference"),("same_ic_nonlinear_component","Nonlinear evolution, same initial state")):
        axes[0].plot(decomposition["times"][mask]/T,1e6*decomposition[key][mask,20,1],label=label)
    axes[0].set(xlabel="t/T1",ylabel="1D midspan contribution (1e-6)",xlim=(0,horizon)); axes[0].legend(fontsize=8)
    for key,label in (("initial_state_component","Initial-state component"),("same_ic_nonlinear_component","Same-IC nonlinear component")):
        index=np.flatnonzero(mask)[-1]
        axes[1].plot(decomposition["x"],1e6*decomposition[key][index,:,1],label=label)
    axes[1].set(xlabel="Material x/L",ylabel="1D final contribution (1e-6)"); axes[1].legend(fontsize=8)
    for ax in axes: ax.grid(alpha=.25)
    figures.append((fig,"one_d_initial_state_evolution_decomposition"))
    for fig,name in figures:
        fig.savefig(destination/(name+".png"))
        fig.savefig(destination/(name+".pdf"),metadata={"CreationDate":None,"ModDate":None})
        plt.close(fig)
    continuation.save(bundle,read_json(bundle/"provenance.json"),summary)
    return {"figures":len(figures),"new_scientific_calls":0}
