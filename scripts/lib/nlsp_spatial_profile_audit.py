"""Scoped saved-data profile audit, dispatched by the existing FEM continuation CLI.

No native, ODE, BVP, static-equilibrium or modal solver is invoked. The historical
recovery remains unchanged; separate diagnostic probes do not replace its data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from scripts.lib import nlsp_profile_fem_diagnostics as fem

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "data/input/nlsp_spatial_profile_audit.json"
OUTPUT = ROOT / "results/nlsp_spatial_profile_audit"
SOURCE_BUNDLES = {
    "FEM3C": {"path": "results/nlsp_nonlinear_dynamic_validation/c6256269eb8143ef",
        "manifest_sha256": "0fa3488d30c1b36de2061894e2f8811443ab80b802f47fddb99fcb30e5677450"},
    "FEM3B": {"path": "results/nlsp_nonlinear_dynamic_long_horizon/7d2b499e6a1eb990",
        "manifest_sha256": "e4bcc291fab5f04a2fb73103c35a2fa6481aecf5f058b8ce6ae4a5661efa0c4e"},
    "FEM3AR": {"path": "results/nlsp_nonlinear_dynamic_3d_fem_resume/f6ac2f2f38510893",
        "manifest_sha256": "187870ece572dfe1b999b83d99016d01c1fdbe10e8835ed45543d1119d6f36ea"}}
POLICIES = {
    "authorization": "explicit_user_NLSP_spatial_profile_postprocessing_2026_10_10",
    "native_time_selection": "nearest_actual_saved_frame; requested and actual times separate; no nodal temporal interpolation",
    "original_figure_time_policy": "preserve historical time-interpolated curves separately",
    "independent_strain_quadrature": "existing positive14point C3D10 reference quadrature",
    "slab_integration_qualification": "hard-binned quadrature estimate, not exact clipped-tetrahedron volume integration",
    "quadratic_sampling_control": "U_eta=eta^2/L on original saved mesh; exact P2 nodal kinematics; no equilibrium solution",
    "quadratic_fit_probe": "one separate 11-column WLS probe adding eta^2,eta*zeta,zeta^2; historical eight-column recovery unchanged"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    return value


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write(path, value):
    Path(path).write_text(json.dumps(plain(value), ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def validate_config(c):
    if (c["schema"] != "nlsp-spatial-profile-audit-v1" or c["section_counts"] != [21, 41, 81]
            or c["snapshot_T1_fractions"] != [0., .25, .5, .75, 1.]
            or c["spatial_comparison_T1_fractions"] != [0., .125, .25]
            or c["FEM_spatial_grid_points"] != 801 or c["one_d_spatial_grid_points"] != 2001
            or c["maximum_main_figures"] != 4 or c["one_d_boundary_partition_ell_c"] != 3.
            or c["primary_section_count"] != 41 or c["source_bundles"] != SOURCE_BUNDLES
            or any(c.get(k) != v for k, v in POLICIES.items())
            or any(c[k] for k in ("smoothing", "phase_amplitude_fitting", "new_physical_solves", "new_meshes"))
            or c["historical_statuses_unchanged"] is not True):
        raise ValueError("Only the bounded saved-profile diagnostic policy is authorized")
    return c


def checked(path, manifest, source, registry):
    path = Path(path).resolve()
    name = path.relative_to(source.resolve()).as_posix()
    expected = manifest.get("artifact_hashes", manifest.get("artifacts", {})).get(name)
    if expected is None or sha(path) != expected:
        raise ValueError("Missing/corrupt immutable source artifact: "+str(path))
    registry[path.relative_to(ROOT).as_posix()] = expected
    return path


def selected_states(source, provenance):
    """Actual frames only; keep original figure's interpolated times distinct."""
    T = 2*np.pi/provenance["config"]["omega1"]
    states = []
    for stage in ("full_period_medium", "medium_refined_time", "fine_refined_time"):
        case = source / "cases" / stage / "nonlinear"
        with np.load(case / "section_history.npz", allow_pickle=False) as z:
            times, increments = z["time"].copy(), z["increments"].copy()
        desired = [0., .25, .5, .75, 1.] if stage == "full_period_medium" else [0., .125, .25]
        if (times.ndim != 1 or not len(times) or increments.shape != times.shape
                or not np.all(np.isfinite(times)) or not np.all(np.isfinite(increments))
                or np.any(np.diff(times) <= 0) or times[0] < 0
                or np.any(increments != increments.astype(int)) or np.any(increments <= 0)
                or np.any(np.diff(increments) <= 0) or times[-1] < max(desired)*T-1e-10):
            raise ValueError("Missing/invalid saved-time coverage: "+stage)
        for tau in desired:
            if tau == 0.:
                increment, actual, index = 10, 0., None
                file = case / "frames" / "step1_inc00010.npz"
            else:
                index = int(np.argmin(abs(times-tau*T)))
                increment, actual = int(increments[index]), float(times[index])
                file = case / "frames" / f"step2_inc{increment:05d}.npz"
            name = stage + "_tau" + str(tau).replace(".", "p")
            states.append({"name": name, "stage": stage,
                "mesh": "fine" if stage == "fine_refined_time" else "medium",
                "requested_tau": tau, "requested_time": tau*T, "actual_time": actual,
                "actual_tau": actual/T, "time_offset": actual-tau*T,
                "time_policy": "actual_STATIC_preload" if tau == 0 else "nearest_actual_native_frame_no_time_interpolation",
                "frame": file.relative_to(ROOT).as_posix(), "increment": increment,
                "step": 1 if tau == 0 else 2, "history_index": index,
                "exact_requested_time_available": actual == tau*T})
    # The mesh comparison uses the same refined-time policy and identical times.
    for tau in (0., .125, .25):
        pair = [r for r in states if r["stage"] in ("medium_refined_time", "fine_refined_time") and r["requested_tau"] == tau]
        if len(pair) != 2 or pair[0]["actual_time"] != pair[1]["actual_time"]:
            raise ValueError("No exact native time pairing for medium/fine diagnostic")
    return T, states


def prepare(config=CONFIG):
    c = validate_config(read(config))
    sources, manifests, registry = {}, {}, {}
    for name, evidence in c["source_bundles"].items():
        source = ROOT/evidence["path"]
        if sha(source/"manifest.json") != evidence["manifest_sha256"]:
            raise ValueError("Immutable parent manifest mismatch: "+name)
        sources[name], manifests[name] = source, read(source/"manifest.json")
    source, manifest = sources["FEM3C"], manifests["FEM3C"]
    provenance = read(checked(source/"provenance.json", manifest, source, registry))
    T, states = selected_states(source, provenance)
    for state in states:
        checked(ROOT/state["frame"], manifest, source, registry)
    for stage in ("full_period_medium", "medium_refined_time", "fine_refined_time"):
        for name in ("section_history.npz", "initial_sections.npz"):
            checked(source/"cases"/stage/"nonlinear"/name, manifest, source, registry)
    for name in ("full_period_comparison.npz", "one_d_all8_spatial.json", "one_d_p48_nonlinear.npz",
                 "one_d_p64_nonlinear.npz", "one_d_p64_linear.npz"):
        checked(source/name, manifest, source, registry)
    for mesh in provenance["source_meshes"].values():
        if sha(ROOT/mesh["include"]) != mesh["include_sha256"]:
            raise ValueError("Frozen mesh include hash mismatch")
        registry[mesh["include"]] = mesh["include_sha256"]
    frozen = ["scripts/analysis/verify_nlsp_nonlinear_static_3d_fem.py",
        "scripts/lib/weakly_nonlinear_spatial_rod.py", "scripts/lib/weakly_nonlinear_planar_dynamics.py",
        "scripts/analysis/simulate_weakly_nonlinear_planar_rod.py"]
    frozen_hashes = {p: sha(ROOT/p) for p in frozen}
    identity = {"config": c, "source_artifacts": registry, "frozen_helpers": frozen_hashes,
        "HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "python": sys.version, "dependencies": {p: importlib.metadata.version(p) for p in ("numpy", "scipy", "matplotlib")}}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
    bundle = OUTPUT/fingerprint
    if (bundle/"manifest.json").exists():
        return bundle, read(bundle/"provenance.json"), validate_cache(bundle)
    if bundle.exists():
        raise ValueError("Interrupted diagnostic bundle: preserve prefix; do not silently overwrite")
    bundle.mkdir(parents=True)
    identity.update(source_geometry=provenance["config"].get("geometry"),
        source_meshes=provenance["source_meshes"], T1=T, selected_states=states,
        initial_checkout_snapshot="results/_smoke/nlsp_spatial_profile_audit/initial_state.json",
        no_scientific_solver_calls=True)
    write(bundle/"provenance.json", identity)
    write(bundle/"config.json", c)
    summary = {"status": "PREPARED", "scientific_calls": 0, "states": {},
        "historical_FEM3C_status_unchanged": True, "strict_float64": "PARTIAL"}
    save(bundle, identity, summary)
    return bundle, identity, summary


def save(bundle, identity, summary):
    write(bundle/"summary.json", summary)
    files = {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(bundle.rglob("*"))
             if p.is_file() and p.name != "manifest.json"}
    write(bundle/"manifest.json", {"schema": "nlsp-profile-diagnostic-artifacts-v1",
        "artifact_hashes": files, "source_manifest_hashes": identity["config"]["source_bundles"],
        "scientific_calls": 0})


def validate_cache(bundle):
    b = Path(bundle)
    for name, digest in read(b/"manifest.json")["artifact_hashes"].items():
        if sha(b/name) != digest:
            raise ValueError("Diagnostic cache artifact mismatch: "+name)
    item = read(b/"provenance.json")
    for evidence in item["config"]["source_bundles"].values():
        if sha(ROOT/evidence["path"]/"manifest.json") != evidence["manifest_sha256"]:
            raise ValueError("Historical manifest changed")
    for registry in ("source_artifacts", "additional_postprocessing_sources", "frozen_helpers"):
        for name, digest in item.get(registry, {}).items():
            if sha(ROOT/name) != digest:
                raise ValueError("Immutable source artifact changed: "+name)
    return read(b/"summary.json")


def arrays_flat(value, prefix="", out=None):
    out = {} if out is None else out
    if isinstance(value, np.ndarray):
        out[prefix] = value
    elif isinstance(value, dict):
        for name, child in value.items():
            arrays_flat(child, prefix+"__"+str(name), out)
    elif isinstance(value, (list, tuple)):
        for i, child in enumerate(value):
            arrays_flat(child, prefix+"__"+str(i), out)
    return out


def store_state(path, result):
    path.mkdir(parents=True, exist_ok=True)
    samples = result.pop("quadrature_samples", None)
    if samples is not None:
        np.savez_compressed(path/"native_quadrature.npz", **samples)
        result["native_quadrature_file"] = "native_quadrature.npz"
    np.savez_compressed(path/"recovery.npz", **arrays_flat(result))
    write(path/"recovery.json", result)
    with (path/"raw_sections.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("sections", "x", "c_eff", "c_small", "width_effective", "width_small",
            "native_E22", "native_grad22", "samples", "mass", "volume", "rank", "condition", "fit_residual_mass_L2"))
        for n, level in result["levels"].items():
            for i, x in enumerate(level["raw_x"]):
                writer.writerow((n, x, *(level[k][i] for k in ("c_eff", "c_small", "width_effective", "width_small")),
                    level["native_columns"]["thickness_green"][i], level["native_columns"]["thickness_small"][i],
                    *(level[k][i] for k in ("fit_sample_count", "fit_mass", "fit_volume", "fit_rank", "fit_condition", "fit_residual_mass_L2"))))


def compute(bundle, item, summary):
    if summary.get("completed"):
        return summary
    prepared = {}
    for name, evidence in item["source_meshes"].items():
        mesh = fem.fem1.single.read_gmsh_inp_mesh_data(ROOT/evidence["include"])
        prepared[name] = fem.prepare_mesh(mesh, 1.)
    start = time.perf_counter()
    code = bundle/"execution_code"; code.mkdir(exist_ok=True)
    for path in (Path(__file__), Path(fem.__file__), ROOT/"scripts/lib/nlsp_profile_1d_diagnostics.py"):
        target = code/path.name
        if target.exists() and sha(target) != sha(path):
            raise ValueError("Interrupted diagnostic code phase changed; preserve its executed snapshot")
        if not target.exists():
            shutil.copyfile(path, target)
    results = {}
    for state in item["selected_states"]:
        name = state["name"]
        if name in summary["states"]:
            results[name] = read(bundle/"states"/name/"recovery.json")
            continue
        with np.load(ROOT/state["frame"], allow_pickle=False) as z:
            if not np.array_equal(z["node_ids"], prepared[state["mesh"]]["node_ids"]):
                raise ValueError("Saved field ordering differs from source mesh")
            U = z["U"].copy()
        result = fem.analyze_state(prepared[state["mesh"]], U,
            section_counts=item["config"]["section_counts"], dense_x=np.linspace(0., 1., 801),
            keep_sample_arrays=True, quadratic_fit_probe=True)
        result["actual_source_state"] = state
        store_state(bundle/"states"/name, result)
        results[name] = result
        summary["states"][name] = {"status": "COMPLETE_DIAGNOSTIC", "actual_source_state": state,
            "result": "states/"+name+"/recovery.json"}
        save(bundle, item, summary)
        print("Saved profile diagnostic: "+name, flush=True)
    comparisons = {}
    for tau in (0., .125, .25):
        names = [s["name"] for s in item["selected_states"] if s["stage"] in ("medium_refined_time", "fine_refined_time") and s["requested_tau"] == tau]
        comparisons[str(tau)] = {"actual_time": summary["states"][names[0]]["actual_source_state"]["actual_time"],
            "levels": {str(n): fem.compare_levels(results[names[0]]["levels"][str(n)],
                results[names[1]]["levels"][str(n)], dense_x=np.linspace(0., 1., 801)) for n in (21, 41, 81)}}
    write(bundle/"mesh_comparison.json", comparisons)
    controls = {name: fem.quadratic_transverse_sampling_control(p) for name, p in prepared.items()}
    np.savez_compressed(bundle/"quadratic_sampling_control.npz", **arrays_flat(controls))
    write(bundle/"quadratic_sampling_control.json", controls)
    staged = ROOT/"results/_smoke/nlsp_spatial_profile_audit/staged_one_d"
    if (staged/"stage_provenance.json").exists() and sha(ROOT/"scripts/lib/nlsp_profile_1d_diagnostics.py") == read(staged/"stage_provenance.json").get("current_helper_sha256", read(staged/"stage_provenance.json").get("helper_sha256")):
        for name in ("one_d_profile_audit.json", "one_d_profile_audit.npz", "one_d_profile_audit.csv", "stage_provenance.json"):
            shutil.copyfile(staged/name, bundle/name)
    else:
        from scripts.lib.nlsp_profile_1d_diagnostics import audit_saved_one_d
        audit_saved_one_d(ROOT/item["config"]["source_bundles"]["FEM3C"]["path"], bundle)
    summary.update(status="DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS", completed=True,
        postprocessing_seconds=time.perf_counter()-start, one_d_c_spatial_status="PARTIAL",
        source_preservation={p: sha(ROOT/p) == h for p, h in item["source_artifacts"].items()})
    if not all(summary["source_preservation"].values()):
        raise ValueError("A source artifact changed during postprocessing")
    save(bundle, item, summary)
    return summary


def complete_saved_evidence(bundle, item, summary):
    """Preserve the old plot's time policy and reproduce its native recovery."""
    source = ROOT/item["config"]["source_bundles"]["FEM3C"]["path"]
    source_manifest = read(source/"manifest.json")
    registry = item.setdefault("additional_postprocessing_sources", {})
    checked(source/"dense_p64.npz", source_manifest, source, registry)
    science = read(source/"provenance.json")["config"]
    nested = item.setdefault("one_d_source_chain", {})
    for key, names in (("source_static", ("one_d_preflight.json", "one_d_p64.npz")),
            ("source_fem1", ("preflight.json",)), ("source_action", ("result.json",))):
        evidence = science[key]
        parent = ROOT/evidence["bundle"]
        expected = evidence["manifest_sha256"]
        if sha(parent/"manifest.json") != expected:
            raise ValueError("Frozen 1D source-chain manifest changed: "+key)
        parent_manifest = read(parent/"manifest.json")
        for name in names:
            checked(parent/name, parent_manifest, parent, registry)
        registry[(parent/"manifest.json").relative_to(ROOT).as_posix()] = expected
        nested[key] = evidence
    phase = bundle/"execution_code"/"postprocessing_extension"
    phase.mkdir(parents=True, exist_ok=True)
    phase_files = (Path(__file__), ROOT/"scripts/lib/nlsp_profile_1d_diagnostics.py",
        ROOT/"scripts/lib/nlsp_profile_figures.py")
    # The original executed numerical code remains in execution_code/ unchanged.
    for path in phase_files:
        target = phase/path.name
        if not target.exists():
            shutil.copyfile(path, target)
    if not (bundle/"actual_time_one_d.json").exists():
        from scripts.lib.nlsp_profile_1d_diagnostics import audit_actual_one_d_states
        audit_actual_one_d_states(source, bundle, item["selected_states"], spatial_points=801)
    if not (bundle/"historical_figure_profiles.json").exists():
        with np.load(source/"full_period_comparison.npz", allow_pickle=False) as z:
            tau = np.asarray(item["config"]["snapshot_T1_fractions"])
            ids = np.array([int(np.argmin(abs(z["times"]-t*item["T1"]))) for t in tau])
            if not np.allclose(z["times"][ids], tau*item["T1"], rtol=0, atol=1e-12):
                raise ValueError("Original representative figure timestamps unavailable")
            np.savez_compressed(bundle/"historical_figure_profiles.npz", x=z["x"],
                times=z["times"][ids], requested_tau=tau,
                three_d_nonlinear_fields=z["three_d_nonlinear_fields"][ids],
                one_d_nonlinear_fields=z["one_d_nonlinear_fields"][ids])
        evidence = {"source": str((source/"full_period_comparison.npz").relative_to(ROOT).as_posix()),
            "time_policy": "historical time-interpolated section curves, not newly observed native nodal frames",
            "spatial_policy": "historical unsmoothed 41-section cubic interpolation",
            "requested_tau": tau, "native_recovery_reproduction": {}}
        for state in item["selected_states"]:
            case = source/"cases"/state["stage"]/"nonlinear"
            filename = "initial_sections.npz" if state["step"] == 1 else "section_history.npz"
            with np.load(case/filename, allow_pickle=False) as z:
                fields = z["fields"] if state["step"] == 1 else z["fields"][state["history_index"]]
                profile = read(bundle/"states"/state["name"]/"recovery.json")["levels"]["41"]["historical_profile"]
                actual = fem.fem2.fem2_static_sample(profile, z["x"])
            differences = np.max(abs(actual-fields), axis=0)
            evidence["native_recovery_reproduction"][state["name"]] = {
                "actual_time": state["actual_time"], "field_order": profile["field_order"],
                "maximum_absolute_difference": differences,
                "within_float64_roundoff": bool(np.max(differences) < 1e-13)}
        if not all(row["within_float64_roundoff"] for row in evidence["native_recovery_reproduction"].values()):
            raise ValueError("Historical native spatial recovery failed reproduction")
        write(bundle/"historical_figure_profiles.json", evidence)
    if not (bundle/"historical_time_policy.json").exists():
        case = source/"cases"/"full_period_medium"/"nonlinear"
        with np.load(case/"initial_sections.npz", allow_pickle=False) as initial, np.load(case/"section_history.npz", allow_pickle=False) as history:
            source_times = np.r_[0., history["time"]]
            source_fields = np.concatenate((initial["fields"][None], history["fields"]), axis=0)
        rows, nearest_fields, left_fields, right_fields = [], [], [], []
        with np.load(bundle/"historical_figure_profiles.npz", allow_pickle=False) as old:
            for i, state in enumerate(s for s in item["selected_states"] if s["stage"] == "full_period_medium"):
                t = float(old["times"][i])
                j = int(np.clip(np.searchsorted(source_times, t, side="right")-1, 0, len(source_times)-2))
                fraction = (t-source_times[j])/(source_times[j+1]-source_times[j])
                replay = (1-fraction)*source_fields[j]+fraction*source_fields[j+1]
                nearest = source_fields[0 if state["step"] == 1 else state["history_index"]+1]
                rows.append({"requested_tau": state["requested_tau"], "requested_time": t,
                    "nearest_actual_time": state["actual_time"],
                    "bracket_times": source_times[j:j+2], "right_fraction": fraction,
                    "linear_time_reproduction_max": float(np.max(abs(replay-old["three_d_nonlinear_fields"][i]))),
                    "historical_minus_nearest_actual_field_max": np.max(abs(old["three_d_nonlinear_fields"][i]-nearest), axis=0)})
                nearest_fields.append(nearest); left_fields.append(source_fields[j]); right_fields.append(source_fields[j+1])
            np.savez_compressed(bundle/"historical_time_policy.npz", x=old["x"],
                requested_times=old["times"], nearest_actual_fields=nearest_fields,
                left_native_fields=left_fields, right_native_fields=right_fields)
        write(bundle/"historical_time_policy.json", {"rows": rows,
            "qualification": "Historical curves are linear time interpolation of already recovered fields; nearest-native differences compare different physical times, not interpolation error bounds"})
    item.setdefault("postprocessing_code_phases", {}).update({
        "original": {p.relative_to(bundle).as_posix(): sha(p)
            for p in sorted((bundle/"execution_code").glob("*.py"))},
        "extension": {p.relative_to(bundle).as_posix(): sha(p) for p in sorted(phase.glob("*.py"))}})
    write(bundle/"provenance.json", item)
    summary["additional_source_preservation"] = {p: sha(ROOT/p) == h for p, h in registry.items()}
    summary["native_recovery_reproduction"] = "PASS"
    summary["actual_time_one_d_scientific_calls"] = 0
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--compute", action="store_true")
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--report-only", type=Path)
    mode.add_argument("--plot-only", type=Path)
    parser.add_argument("--config", type=Path, default=CONFIG)
    args = parser.parse_args(argv)
    if args.report_only or args.plot_only:
        bundle = (args.report_only or args.plot_only).resolve()
        summary = validate_cache(bundle)
        if args.plot_only:
            plot(bundle)
            save(bundle, read(bundle/"provenance.json"), summary)
    else:
        bundle, item, summary = prepare(args.config)
        if args.compute:
            summary = compute(bundle, item, summary)
            summary = complete_saved_evidence(bundle, item, summary)
            plot(bundle)
            save(bundle, item, summary)
    print(json.dumps({"bundle": str(bundle), "status": summary["status"], "scientific_calls": 0}, indent=2))
    return summary


def plot(bundle):
    """Render only completed saved diagnostic arrays, with no solver calls."""
    from scripts.lib.nlsp_profile_figures import render_bundle
    return render_bundle(bundle)
