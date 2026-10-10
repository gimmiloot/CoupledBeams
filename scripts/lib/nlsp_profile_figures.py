"""Read-only saved-data profile metrics and four diagnostic figures.

All curves preserve their dimensional scales.  No smoothing, alignment,
solver, meshing, eigenanalysis or new physical trajectory is performed here.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


def _read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _state_records(bundle):
    return {p.parent.name: _read(p) for p in sorted((Path(bundle)/"states").glob("*/recovery.json"))}


def _slim_roughness(row):
    return {k: row[k] for k in ("sampled_min", "sampled_max", "mean_over_sample_span",
        "total_variation", "first_derivative_L2", "second_derivative_RMS", "raw_extrema_count")}


def aggregate_metrics(bundle):
    """Write scalar diagnostics only; never modify scientific summary/manifest."""
    bundle = Path(bundle)
    states = _state_records(bundle)
    if not states:
        raise ValueError("No completed saved recovery states")
    rows, by_state = [], {}
    for name, state in states.items():
        source = state["actual_source_state"]
        levels = {}
        for level, data in state["levels"].items():
            probe = state.get("quadratic_transverse_fit_probe", {}).get(level, {})
            interp = data["interpolation"]
            c_eff, small = np.asarray(data["c_eff"]), np.asarray(data["c_small"])
            grad = np.asarray([r["affine_gradient_local"] for r in data["historical_profile"]["section_rows"]])
            stretch = np.asarray([r["transverse_stretch"] for r in data["historical_profile"]["section_rows"]])
            finite_difference_prediction = .5*np.sum(grad[:, :, 1]**2, axis=1)-.5*c_eff**2-.5*stretch[:, 0, 1]**2
            record = {"state": name, "mesh": source["mesh"], "requested_tau": source["requested_tau"],
                "actual_tau": source["actual_tau"], "actual_time": source["actual_time"],
                "time_offset": source["time_offset"], "sections": int(level),
                **data["summary"],
                "fit_condition_min": float(min(data["fit_condition"])),
                "fit_condition_max": float(max(data["fit_condition"])),
                "fit_samples_min": int(min(data["fit_sample_count"])),
                "fit_samples_max": int(max(data["fit_sample_count"])),
                "slab_volume_ratio_min": float(min(data["slab_volume_to_nominal_ratio"])),
                "slab_volume_ratio_max": float(max(data["slab_volume_to_nominal_ratio"])),
                "transverse_centroid_thickness_max": float(np.max(abs(np.asarray(data["geometric_centroid_local"])[:, 1]))),
                "transverse_centroid_width_max": float(np.max(abs(np.asarray(data["geometric_centroid_local"])[:, 2]))),
                "interior_cubic_extrema_count": interp["interior_extrema_count"],
                "additional_cubic_extrema_count": interp["additional_interior_extrema_count"],
                "cubic_neighbor_range_overshoot_max": interp["maximum_neighbor_range_overshoot"],
                "cubic_minus_linear_max": interp["cubic_minus_linear_max"],
                "raw_knot_reproduction_max_error": interp["raw_knot_reproduction_max_error"],
                "finite_difference_identity_error": float(np.max(abs(c_eff-small-finite_difference_prediction))),
                "fit_displacement_residual_RMS_max": float(max(data["fit_residual_displacement_RMS"])),
                "fit_transverse_gradient_residual_RMS_max": float(max(data["fit_transverse_gradient_residual_RMS"])),
                "c_eff_roughness": _slim_roughness(data["roughness"]["c_eff"]),
                "native_small_roughness": _slim_roughness(data["roughness"]["native_thickness_small"]),
                "native_green_roughness": _slim_roughness(data["roughness"]["native_thickness_green"]),
                "native_polar_roughness": _slim_roughness(data["roughness"]["native_point_polar_thickness"]),
                "sampling_correlations_not_causal_proof": data["associations_not_causal_proof"]}
            if probe.get("status") == "COMPLETE_DIAGNOSTIC_PROBE":
                condition = np.asarray(probe["fit_condition"])
                record["separate_quadratic_probe"] = {
                    "fit_condition_min": float(condition.min()), "fit_condition_max": float(condition.max()),
                    "minimum_fit_rank": int(min(probe["fit_rank"])),
                    "raw_c_eff_roughness": _slim_roughness(probe["roughness"]["c_eff"]),
                    "historical_minus_probe_raw_max": probe["historical_minus_probe_raw_max"],
                    "historical_minus_probe_raw_L2": probe["historical_minus_probe_raw_L2"],
                    "fit_displacement_residual_RMS_max": float(max(probe["fit_residual_displacement_RMS"])),
                    "fit_transverse_gradient_residual_RMS_max": float(max(probe["fit_transverse_gradient_residual_RMS"])),
                    "alias_coefficient_identity_max_error": probe["nested_fit_alias_coefficient_identity_max_error"],
                    "alias_c_small_max": float(np.max(abs(np.asarray(probe["omitted_quadratic_alias_c_small"])))),
                    "TV_ratio_to_historical": probe["roughness"]["c_eff"]["total_variation"]/data["roughness"]["c_eff"]["total_variation"] if data["roughness"]["c_eff"]["total_variation"] else None,
                    "replaces_original": False}
            levels[level] = record
            flat = {k: v for k, v in record.items() if not isinstance(v, dict)}
            flat.update({"c_eff_"+k: v for k, v in record["c_eff_roughness"].items()})
            rows.append(flat)
        by_state[name] = {"source": source, "levels": levels,
            "window_41_vs_81": {"curve_metrics": state["window_comparisons"]["41_vs_81"]["metrics"],
                "control_fields": state["window_comparisons"]["41_vs_81"]["control_fields"],
                "native_measures": state["window_comparisons"]["41_vs_81"]["native_measures"]}}
    meshes = _read(bundle/"mesh_comparison.json") if (bundle/"mesh_comparison.json").exists() else {}
    mesh_slim = {tau: {"actual_time": data["actual_time"], "levels": {
        level: {key: compared[key] for key in ("metrics", "control_fields", "native_measures", "native_strain_comparison_span")}
        for level, compared in data["levels"].items()}} for tau, data in meshes.items()}
    controls = _read(bundle/"quadratic_sampling_control.json") if (bundle/"quadratic_sampling_control.json").exists() else {}
    probe_slim = {mesh: {"C3D10_exact_gradient_max_error": data["C3D10_exact_gradient_max_error"],
        "coefficient_inverse_length": data["coefficient_inverse_length"], "levels": {
            level: {key: row[key] for key in ("fit_minus_direct_mean_max", "native_vs_exact_bin_mean_max_error",
                "quadratic_probe_center_gradient_max_error")}
            for level, row in data["levels"].items()}} for mesh, data in controls.items()}
    result = {"version": "saved-spatial-profile-scalar-audit-v1", "states": by_state,
        "medium_fine_same_policy": mesh_slim, "synthetic_quadratic_sampling": probe_slim,
        "scientific_calls": 0, "smoothing": False,
        "qualification": "cut-slab quadrature means include geometric sampling bias; richer fit is diagnostic, not continuum truth"}
    (bundle/"numbers.json").write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with (bundle/"FEM_profile_metrics.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)
    return result


def _save(fig, bundle, name):
    fig.savefig(bundle/(name+".png"), dpi=180, bbox_inches="tight")
    fig.savefig(bundle/(name+".pdf"), bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})


def _style(ax, ylabel, title=None):
    ax.set_xlabel(r"material coordinate $s/L$")
    ax.set_ylabel(ylabel)
    if title: ax.set_title(title, fontsize=10)
    ax.grid(alpha=.22)
    ax.set_xlim(0, 1)


def _loaded_arrays(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key] for key in archive.files}


def render_bundle(bundle):
    """Produce at most four PNG/PDF pairs from completed saved diagnostics."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    bundle = Path(bundle)
    aggregate_metrics(bundle)
    states = _state_records(bundle)
    one = _loaded_arrays(bundle/"one_d_profile_audit.npz")
    one_metrics = _read(bundle/"one_d_profile_audit.json")
    actual_path = bundle/"actual_time_one_d.npz"
    actual_one = _loaded_arrays(actual_path) if actual_path.exists() else None
    colors = plt.get_cmap("tab10")(np.arange(5))
    tau = one["time_fractions"]
    x = one["x"]

    # Figure 1: all five actual FEM fields; original 41/81 points remain visible.
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey=True)
    full_states = sorted((state for state in states.values() if state["actual_source_state"]["stage"] == "full_period_medium"),
                         key=lambda state: state["actual_source_state"]["requested_tau"])
    for k, state in enumerate(full_states):
        ax = axes.flat[k]
        source = state["actual_source_state"]
        level = state["levels"]["41"]
        other = state["levels"]["81"]
        interp = level["interpolation"]
        ax.plot(interp["dense_x"], np.asarray(interp["cubic"])*1e5, color="0.55", lw=1, label="original cubic policy (41)")
        ax.plot(level["raw_x"], np.asarray(level["c_eff"])*1e5, "o", ms=3, color="#1565c0", label="raw 41")
        ax.plot(other["raw_x"], np.asarray(other["c_eff"])*1e5, ".", ms=3, color="#ef6c00", label="raw 81")
        probe = state["quadratic_transverse_fit_probe"]["41"]
        ax.plot(probe["raw_x"], np.asarray(probe["c_eff"])*1e5, "x", ms=3, color="#388e3c", label="separate P2-fit probe (41)")
        # Exact actual-time companion, when supplied by root, avoids disguising
        # requested-vs-native offsets in an apparently same-time overlay.
        if actual_one is not None:
            names = actual_one["names"].astype(str).tolist()
            index = names.index(source["name"])
            if abs(actual_one["times"][index]-source["actual_time"]) > 1e-10:
                raise ValueError("Actual-time 1D companion does not match the native FEM state")
            ax.plot(actual_one["x"], actual_one["fields"][index, :, 3]*1e5, color="black", lw=1.3, label="1D c, same actual time")
        else:
            index = int(np.argmin(abs(tau-source["requested_tau"])))
            ax.plot(x, one["p64_fields"][index, :, 3]*1e5, color="black", lw=1.3, label="1D c, requested time")
        title = f"requested τ={source['requested_tau']:g}; actual τ={source['actual_tau']:.6f}"
        _style(ax, r"$c_{eff},c$ [dimensionless] $\times10^5$", title)
    legend_ax = axes.flat[5]
    legend_ax.axis("off")
    handles, labels = axes.flat[0].get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="upper left", frameon=False, fontsize=9)
    legend_ax.text(.03, .43, "Raw sections are not smoothed.\nP2 fit is a separate diagnostic.\n3D c_eff is not the M–H coordinate c.\nNative states differ slightly from requested times.",
                   transform=legend_ax.transAxes, va="top", fontsize=9)
    fig.suptitle("Raw effective contraction, original cubic policy and fit sensitivity", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .96))
    _save(fig, bundle, "raw_vs_interpolated_contraction"); plt.close(fig)

    # Figure 2: only genuinely available same-policy medium/fine quarter states.
    fig, axes = plt.subplots(3, 2, figsize=(11, 9), sharex=True)
    for col, fraction in enumerate((0., .25)):
        selected = {}
        for mesh, stage in (("medium", "medium_refined_time"), ("fine", "fine_refined_time")):
            selected[mesh] = next(state for state in states.values()
                if state["actual_source_state"]["stage"] == stage and state["actual_source_state"]["requested_tau"] == fraction)
        actual_tau = selected["medium"]["actual_source_state"]["actual_tau"]
        for mesh, color in (("medium", "#1565c0"), ("fine", "#ef6c00")):
            data = selected[mesh]["levels"]["41"]
            xx = data["raw_x"]; native = data["native_columns"]
            axes[0, col].plot(xx, np.asarray(data["c_eff"])*1e5, "o-", ms=2.8, color=color, lw=.8, label=mesh+" fitted polar")
            axes[0, col].plot(xx, np.asarray(native["thickness_green"])*1e5, "--", color=color, lw=1.1, label=mesh+" volume GL thickness")
            axes[0, col].plot(xx, np.asarray(native["point_polar_thickness"])*1e5, ":", color=color, lw=1.1, label=mesh+" volume point-polar")
            axes[1, col].plot(xx, np.asarray(data["c_small"])*1e5, "o-", ms=2.8, color=color, lw=.8, label=mesh+" fitted gradient 22")
            axes[1, col].plot(xx, np.asarray(native["thickness_small"])*1e5, "--", color=color, lw=1.1, label=mesh+" volume gradient 22")
            axes[2, col].plot(xx, np.asarray(data["width_effective"])*1e5, "o-", ms=2.8, color=color, lw=.8, label=mesh+" fitted width polar")
            axes[2, col].plot(xx, np.asarray(native["width_green"])*1e5, "--", color=color, lw=1.1, label=mesh+" volume GL width")
        _style(axes[0, col], r"finite transverse measures $\times10^5$", f"actual τ={actual_tau:.6f}")
        _style(axes[1, col], r"small transverse gradient $\times10^5$")
        _style(axes[2, col], r"width measures $\times10^5$")
        for row in range(3): axes[row, col].legend(fontsize=6.8, loc="best")
    fig.suptitle("Independent C3D10 strains and section-fit measures: same-time medium/fine", fontsize=12)
    fig.text(.5, .012, "Volume averages are hard-binned quadrature estimates; finite measures differ. No boundary strains imposed.", ha="center", fontsize=8)
    fig.tight_layout(rect=(0, .025, 1, .965))
    _save(fig, bundle, "native_transverse_strain_reconstruction"); plt.close(fig)

    # Figure 3: fixed dimensional scales, including visible numerical tails.
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True)
    for k, (fraction, color) in enumerate(zip(tau, colors)):
        label = f"τ={fraction:g}"
        axes[0, 0].plot(x, one["p64_fields"][k, :, 3]*1e5, color=color, label=label)
        axes[0, 0].plot(x, one["p64_minus_nu_Gamma1"][k]*1e5, "--", color=color, lw=.85)
        axes[0, 1].plot(x, one["p64_Gamma1"][k]*1e5, color=color, label=label)
        axes[0, 2].plot(x, one["p64_N"][k]*1e6, color=color, label=label)
        axes[1, 0].plot(x, one["p64_c_plus_nu_Gamma1"][k]*1e6, color=color, label=label)
        axes[1, 1].plot(x, (one["p48_fields"][k, :, 3]-one["p64_fields"][k, :, 3])*1e8, color=color, label=label)
        axes[1, 2].plot(x, one["p64_c_tail"][k]*1e8, color=color, label=label)
    descriptions = ((r"$c,-\nu\Gamma_1$ $\times10^5$", "solid c; dashed −νΓ₁"),
        (r"$\Gamma_1$ $\times10^5$", "full kinematic axial strain"),
        (r"$N$ [normalized force] $\times10^6$", "N=C(Γ₁+νc), diagnostic"),
        (r"$c+\nu\Gamma_1$ $\times10^6$", "dynamic departure from static approximation"),
        (r"$c_{48}-c_{64}$ $\times10^8$", "observed p sensitivity; historical gate PARTIAL"),
        (r"$c_{64}$ high-order tail $\times10^8$", "tail component; absolute enlarged scale"))
    for ax, (ylabel, title) in zip(axes.flat, descriptions):
        _style(ax, ylabel, title)
    axes[0, 0].legend(fontsize=7)
    for ax in (axes[0, 0], axes[1, 0]):
        width = one_metrics["regional_policy"]["boundary_width"]
        ax.axvspan(0, width, color="0.85", alpha=.3)
        ax.axvspan(1-width, 1, color="0.85", alpha=.3)
    fig.suptitle("Saved 1D contraction dynamics: profiles, finite strain and spatial sensitivity", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .96))
    _save(fig, bundle, "one_d_contraction_physics"); plt.close(fig)

    # Figure 4: displacement sign is not strain/stress sign; N is not constant.
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=False)
    for k, (fraction, color) in enumerate(zip(tau, colors)):
        label = f"τ={fraction:g}"
        axes[0, 0].plot(x, one["p64_fields"][k, :, 0]*1e6, color=color, label=label)
        axes[0, 0].plot(x, one["p64_classical_u"][k]*1e6, "--", color=color, lw=.8)
        axes[0, 1].plot(x, one["p64_gradients"][k, :, 0]*1e5, color=color, label=label)
        axes[0, 1].plot(x, one["p64_classical_u_s"][k]*1e5, "--", color=color, lw=.8)
        axes[0, 2].plot(x, one["p64_Gamma1"][k]*1e5, color=color, label=label)
        axes[1, 0].plot(x, one["p64_N"][k]*1e6, color=color, label=label)
        axes[1, 1].plot(x, one["p64_fields"][k, :, 2]*1e2, color=color, label=label)
        axes[1, 1].plot(x, one["p64_gradients"][k, :, 1]*1e2, "--", color=color, lw=.8)
    axes[1, 1].plot(x, one["analytic_initial_linear_fields"][:, 2]*1e2, ":", color="black", lw=1., label="analytic linear STATIC θ at τ=0")
    descriptions = ((r"$u$ [normalized length] $\times10^6$", "solid 1D; dashed classical bending benchmark"),
        (r"$u_s$ $\times10^5$", "solid 1D; dashed mean−½wₛ²"),
        (r"$\Gamma_1$ $\times10^5$", "actual axial strain; signs retained"),
        (r"$N$ [normalized force] $\times10^6$", "dynamic axial force; signs retained"),
        (r"$\theta,w_s$ [rad/slope] $\times10^2$", "solid θ; dashed wₛ; shear permits difference"))
    for ax, (ylabel, title) in zip(axes.flat, descriptions):
        _style(ax, ylabel, title)
    axes[0, 0].legend(fontsize=7)
    ax = axes[1, 2]
    snapshots = one_metrics["cases"]["64"]["snapshots"]
    for name, marker in (("u", "o"), ("w", "s"), ("theta", "^"), ("c", "d")):
        ax.plot(tau, [snapshot["symmetry"][name]["relative_error"] for snapshot in snapshots], marker+"-", label=name)
    ax.set_xlabel(r"$t/T_1$"); ax.set_ylabel("relative parity residual")
    ax.set_title("sampled symmetry; own full-profile scales", fontsize=10)
    ax.grid(alpha=.22); ax.legend(fontsize=7)
    fig.suptitle("Mechanical meaning of axial displacement and independent section rotation", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .96))
    _save(fig, bundle, "mechanical_consistency_u_theta"); plt.close(fig)
    return {"figures": ["raw_vs_interpolated_contraction", "native_transverse_strain_reconstruction",
        "one_d_contraction_physics", "mechanical_consistency_u_theta"],
        "scientific_calls": 0, "source": "saved diagnostics only", "smoothed": False}
