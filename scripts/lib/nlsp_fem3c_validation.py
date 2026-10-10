"""Bounded FEM-3C preset of the existing continuation entry point.

This module owns authorization, source selection and attempt accounting. Native
generation, element quadrature, parsers, section recovery and the 1D solver are
reused; no new mechanical or finite-element formulation is implemented here.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.metadata
import json
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as resume
from scripts.lib import nlsp_fem3b_continuation as previous

base = resume.base
ROOT = resume.ROOT
CONFIG = ROOT / "data/input/nlsp_nonlinear_dynamic_validation.json"
OUTPUT = ROOT / "results/nlsp_nonlinear_dynamic_validation"
read_json, write_json, sha = base.read_json, base.write_json, base.sha
AUTHORIZATION = "explicit_user_FEM3C_2026_10_09"
CONTROL_STAGES = ("medium_refined_time", "fine_refined_time")
FULL_STAGE = "full_period_medium"
CASE_ORDER = tuple((stage, kind) for stage in (*CONTROL_STAGES, FULL_STAGE)
                   for kind in ("linear", "nonlinear"))
STATUS_NAMES = ("SOURCE_PRESERVATION", "INPUT_PROTOCOL", "TEMPORAL_CONTROL",
    "SPATIAL_CONTROL", "ROBUSTNESS", "FULL_PERIOD_1D", "FULL_PERIOD_3D",
    "SEVEN_FIELD_RECOVERY", "ENERGY_DIAGNOSTICS", "VERIFICATION_SUMMARY")


def validate_config(c):
    if c.get("schema") != "nlsp-fem3c-validation-v1":
        raise ValueError("FEM-3C schema mismatch")
    auth = c["authorization"]
    if (auth["id"] != AUTHORIZATION or auth["maximum_production_CCX_jobs"] != 6
        or auth["maximum_nonlinear_1D_integrations"] != 2 or auth["automatic_retry"]):
        raise ValueError("Separate, bounded FEM-3C authorization required")
    expected = {
        "medium_refined_time": {"mesh_level": "medium", "horizon_T1": .25,
            "initial_T1_fraction": 1/8000, "maximum_T1_fraction": 1/4000,
            "output_frequency": 2, "maximum_increments": 2000},
        "fine_refined_time": {"mesh_level": "fine", "horizon_T1": .25,
            "initial_T1_fraction": 1/8000, "maximum_T1_fraction": 1/4000,
            "output_frequency": 2, "maximum_increments": 2000},
        FULL_STAGE: {"mesh_level": "medium", "horizon_T1": 1.,
            "initial_T1_fraction": 1/4000, "maximum_T1_fraction": 1/2000,
            "output_frequency": 5, "maximum_increments": 3000,
            "conditional_on_robustness": True}}
    if c["stages"] != expected:
        raise ValueError("Frozen FEM-3C time/mesh/output policies changed")
    if c["dynamic"] != {"alpha": 0, "minimum_initial_fraction": 1e-4,
        "release": "OP=NEW plus zero GRAV; STEP AMPLITUDE=STEP"}:
        raise ValueError("Release/integrator policy changed")
    comparison = c["comparison"]
    for name, value in {"T1_fraction": .25, "points": 201,
        "primary_interpolation": "linear", "diagnostic_interpolation": "PCHIP",
        "baseline_model_discrepancy": 2.824717e-7, "temporal_ratio_limit": .25,
        "spatial_ratio_limit": .25, "interpolation_ratio_to_effect_limit": .25,
        "interpolation_ratio_to_baseline_limit": .25,
        "phase_amplitude_fitting": False}.items():
        if comparison.get(name) != value:
            raise ValueError("Predeclared robustness/interpolation policy changed: " + name)
    if (c["threads"] != 1 or c["job_timeout_seconds"] != 5400
        or c["job_memory_limit_bytes"] != 4*1024**3
        or c["numerical_budget_seconds"] != 24000):
        raise ValueError("FEM-3C bounded resources changed")
    if (c["execution_mode"] != "EXPLORATORY_NOT_CERTIFIED" or c["admitted"] is not False
        or any(c[k] for k in ("new_meshes", "new_modal_jobs", "new_static_only_jobs"))):
        raise ValueError("Unauthorized mechanical/scientific extension")
    if c["one_d"] != {"main_degree": 64, "optional_spatial_degree": 48,
        "target_T1_fraction": 1., "initial_policy": "reuse_saved_static_coordinates_exactly_no_projection",
        "time_level": "tight", "preserve_accepted_Radau_dense_polynomials": True}:
        raise ValueError("Frozen full-period 1D policy changed")
    obs = c["seven_field_observations"]
    if obs != {"u": .25, "w": .5, "v": .5, "Phi": .25, "psi": .25,
        "theta": .25, "c": .25, "inactive_one_d_fields": ["v", "Phi", "psi"],
        "three_d_c": "effective_contraction_proxy_not_identical_generalized_coordinate"}:
        raise ValueError("Preselected seven-field observations changed")
    if c["full_period_comparison"] != {"points": 401, "primary_interpolation": "linear",
        "diagnostic_interpolation": "PCHIP", "snapshot_T1_fractions": [0., .25, .5, .75, 1.],
        "phase_amplitude_fitting": False}:
        raise ValueError("Predeclared full-period display/interpolation grid changed")
    return c


def load_parent(c):
    parent = ROOT / c["parent_completed"]["bundle"]
    if sha(parent / "manifest.json") != c["parent_completed"]["manifest_sha256"]:
        raise ValueError("Immutable FEM-3B parent manifest changed")
    summary = previous.validate_cache(parent)
    if (not summary.get("completed") or summary["overall"] != "FEM3B_DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS"
        or summary["job_calls"]["CCX_production"] != 2
        or any(summary["cases"][k]["status"] != "PASS" for k in ("linear", "nonlinear"))):
        raise ValueError("Incomplete or wrong FEM-3B parent")
    item = read_json(parent / "provenance.json")
    return parent, item, summary, copy.deepcopy(item["config"])


def source_for_level(science, level):
    if level not in ("medium", "fine"):
        raise ValueError("FEM-3C only reuses saved medium/fine meshes")
    old, medium_source, medium_mesh, medium_audit = base.verify_sources(science)
    if level == "medium":
        source, mesh, audit = medium_source, medium_mesh, medium_audit
    else:
        sources = base.base.load_fem2_sources(old["science_config"])
        source, mesh, audit = base.base.fem2_source_mesh(old["science_config"], level, sources)
    expected = {"medium": (5649, 3120), "fine": (11553, 6670)}[level]
    if (len(mesh.nodes), len(mesh.solid_elements)) != expected:
        raise ValueError("Saved mesh count mismatch: " + level)
    for kind in ("linear", "nonlinear"):
        path = ROOT / science["source_resume"]["bundle"] / "cases" / level / kind
        for name in ("static_nodal_results.npz", "recovered_sections.json"):
            if not (path / name).is_file():
                raise ValueError("Missing immutable static preload: " + str(path / name))
    return old, source, mesh, audit


def case_science(science, c, stage):
    if stage not in c["stages"]:
        raise ValueError("Unauthorized FEM-3C stage")
    config = copy.deepcopy(science)
    policy = c["stages"][stage]
    config.update(mesh_level=policy["mesh_level"], horizon_T1=policy["horizon_T1"],
        job_timeout_seconds=c["job_timeout_seconds"],
        numerical_budget_seconds=c["numerical_budget_seconds"])
    config["dynamic"] = {**c["dynamic"], **{k: policy[k] for k in
        ("initial_T1_fraction", "maximum_T1_fraction", "output_frequency", "maximum_increments")}}
    return config


def write_case_input(path, science, old, source, mesh, audit, nonlinear):
    """Use the historical generator, then adapt only DYNAMIC output cadence."""
    gate = base.write_input(path, science, old, source, mesh, audit, nonlinear)
    path = Path(path)
    text = path.read_text(encoding="utf8")
    static, dynamic = text.split("*END STEP\n", 1)
    frequency = science["dynamic"]["output_frequency"]
    dynamic = re.sub(r"FREQUENCY=1\b", "FREQUENCY=" + str(frequency), dynamic)
    # Repeat the unchanged inherited element outputs explicitly, so every FRD
    # field follows the declared DYNAMIC cadence rather than STATIC frequency1.
    marker = f"*NODE FILE, GLOBAL=YES, FREQUENCY={frequency}"
    dynamic = dynamic.replace(marker, f"*EL FILE, GLOBAL=YES, FREQUENCY={frequency}\nS,E,ENER\n" + marker, 1)
    text = static + "*END STEP\n" + dynamic
    path.write_text(text, encoding="utf8")
    safety = resume.output_safety(text)
    base.input_contract(text, science)
    if any(f"FREQUENCY={frequency}" not in line for line in dynamic.splitlines()
           if line.startswith(("*NODE FILE", "*EL FILE", "*NODE PRINT", "*EL PRINT"))):
        raise ValueError("Dynamic output cadence was not applied consistently")
    return {**gate, **safety, "static_output_frequency": 1,
        "dynamic_output_frequency": frequency, "final_frame_required": True,
        "physical_generation_unchanged": True, "output_only_adapter": True}


def save(bundle, item, summary):
    base.finalize(Path(bundle), item, summary)


def validate_cache(bundle):
    bundle = Path(bundle)
    manifest = read_json(bundle / "manifest.json")
    for name, digest in manifest["artifact_hashes"].items():
        if sha(bundle / name) != digest:
            raise ValueError("FEM-3C artifact hash mismatch: " + name)
    item = read_json(bundle / "provenance.json")
    load_parent(validate_config(item["validation_config"]))
    for level, evidence in item["source_meshes"].items():
        if sha(ROOT / evidence["include"]) != evidence["include_sha256"]:
            raise ValueError("Frozen source mesh changed: " + level)
    summary = read_json(bundle / "summary.json")
    for name in ("pre_fem_decision", "full_period_decision"):
        key = name + "_sha256"
        if key in summary and sha(bundle / (name + ".json")) != summary[key]:
            raise ValueError("Frozen pre-result decision changed: " + name)
    return summary


def existing_attempt(c):
    if not OUTPUT.exists():
        return None
    for bundle in sorted(OUTPUT.iterdir()):
        if not bundle.is_dir() or not (bundle / "provenance.json").exists():
            continue
        item = read_json(bundle / "provenance.json")
        if item.get("authorization", {}).get("id") != AUTHORIZATION:
            continue
        if item["validation_config"] != c:
            raise ValueError("FEM-3C authorization already used with another config")
        if not (bundle / "manifest.json").exists():
            raise RuntimeError("Interrupted unmanifested FEM-3C attempt: no automatic retry")
        return bundle, item, validate_cache(bundle)
    return None


def prepare_stage(config_path=CONFIG):
    c = validate_config(read_json(config_path))
    found = existing_attempt(c)
    if found:
        return found
    parent, parent_item, parent_summary, science = load_parent(c)
    meshes = {}
    for level in ("medium", "fine"):
        old, source, mesh, audit = source_for_level(science, level)
        meshes[level] = {"include": (source / "solid_mesh.inp").relative_to(ROOT).as_posix(),
            "include_sha256": sha(source / "solid_mesh.inp"), "audit": audit,
            "static_reference": (ROOT / science["source_resume"]["bundle"] /
                "cases" / level).relative_to(ROOT).as_posix()}
    helpers = [Path(__file__), Path(resume.__file__), Path(base.__file__),
        Path(base.one.__file__), Path(base.io.__file__), Path(base.base.__file__),
        Path(base.base.fem1.__file__), Path(base.one.dynamics.__file__),
        Path(base.one.rod.__file__), Path(base.one.runner.__file__), Path(previous.__file__),
        ROOT / "scripts/lib/nlsp_fem3b_diagnostics.py"]
    for name in ("nlsp_fem3c_diagnostics.py", "nlsp_fem3c_1d.py"):
        candidate = Path(__file__).with_name(name)
        if candidate.exists():
            helpers.append(candidate)
    binary = Path(old["science_config"]["ccx_exe"])
    item = {"config": science, "validation_config": c, "authorization": c["authorization"],
        "parent_completed": c["parent_completed"], "source_meshes": meshes,
        "source_mesh_sha256": parent_item["source_mesh_sha256"],
        "config_sha256": sha(config_path), "solver_sha256": sha(binary),
        "runtime_DLLs": {p.name: sha(p) for p in sorted(binary.parent.glob("*.dll"))},
        "helper_sha256": {p.relative_to(ROOT).as_posix(): sha(p) for p in helpers},
        "python": sys.version, "dependencies": {p: importlib.metadata.version(p)
            for p in ("numpy", "scipy", "matplotlib")},
        "HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()}
    initial_snapshot = ROOT / "results/_smoke/fem3c_source/initial_state.json"
    if initial_snapshot.is_file():
        item["initial_checkout_snapshot"] = {"path": initial_snapshot.relative_to(ROOT).as_posix(),
            "sha256": sha(initial_snapshot)}
    if item["solver_sha256"] != parent_item["solver_sha256"] or item["runtime_DLLs"] != parent_item["runtime_DLLs"]:
        raise ValueError("Previously confirmed native binary/runtime identity changed")
    key = hashlib.sha256(json.dumps(item, sort_keys=True).encode()).hexdigest()[:16]
    bundle = OUTPUT / key
    if bundle.exists() and any(bundle.iterdir()):
        raise RuntimeError("Existing scientific artifacts cannot be overwritten")
    runtime = {k: parent_summary["cases"][k]["job"]["seconds"] for k in ("linear", "nonlinear")}
    estimates = {"medium_refined_time": sum(runtime.values()) * 2,
        "fine_refined_time": sum(runtime.values()) * 4.5,
        FULL_STAGE: sum(runtime.values()) * 4}
    estimate_total = sum(estimates.values()) + c["planning"]["postprocessing_allowance_seconds"]
    if estimate_total > c["numerical_budget_seconds"]:
        raise ValueError("Declared resource estimate exceeds bounded program budget")
    free = shutil.disk_usage(ROOT).free
    if free < 30*1024**3:
        raise ValueError("Insufficient free disk for bounded native-output series")
    bundle.mkdir(parents=True)
    write_json(bundle / "provenance.json", item)
    write_json(bundle / "validation_config.json", c)
    write_json(bundle / "config.json", science)
    sources = {"FEM3B": c["parent_completed"], **{k: science[k] for k in
        ("source_resume", "source_static", "source_fem1", "source_action")}}
    sources["FEM3AR"] = parent_item["parent_completed"]
    preservation = {}
    for name, evidence in sources.items():
        source_manifest = read_json(ROOT / evidence["bundle"] / "manifest.json")
        preservation[name] = {**evidence, "actual_manifest_sha256": sha(ROOT / evidence["bundle"] / "manifest.json"),
            "artifacts": len(source_manifest.get("artifact_hashes", source_manifest.get("artifacts", {}))),
            "verification": "Every inherited artifact SHA verified by existing source/cache loaders; no source replay"}
    write_json(bundle / "source_preservation.json", preservation)
    (bundle / "execution_code").mkdir()
    for path in helpers:
        shutil.copyfile(path, bundle / "execution_code" / path.name)
    shutil.copyfile(parent / "protocol_evidence.json", bundle / "protocol_evidence.json")
    decision = {"before_any_FEM3C_result": True, "authorization": AUTHORIZATION,
        "comparison": c["comparison"], "stages": c["stages"],
        "old_actual_CCX_seconds": runtime, "planning_CCX_estimates_seconds": estimates,
        "postprocessing_allowance_seconds": c["planning"]["postprocessing_allowance_seconds"],
        "planning_total_seconds": estimate_total, "estimate_is_not_guarantee": True,
        "fine_factor_rationale": "2x dynamic increments, 2.14x elements/2.05x nodes, conservative 4.5x old-medium native runtime; output frequency2 reduces I/O but no savings assumed",
        "per_job_timeout_seconds": c["job_timeout_seconds"],
        "total_numerical_budget_seconds": c["numerical_budget_seconds"],
        "disk_free_bytes": free, "full_period_jobs_conditionally_authorized_only": True}
    write_json(bundle / "pre_fem_decision.json", decision)
    gates = {}
    for stage in (*CONTROL_STAGES, FULL_STAGE):
        sc = case_science(science, c, stage)
        old, source, mesh, audit = source_for_level(science, sc["mesh_level"])
        for kind in ("linear", "nonlinear"):
            path = bundle / "input_gate" / stage / (kind + ".inp")
            path.parent.mkdir(parents=True, exist_ok=True)
            gates[stage + "/" + kind] = write_case_input(path, sc, old, source, mesh, audit, kind == "nonlinear")
    write_json(bundle / "input_gate.json", gates)
    summary = {"authorization": c["authorization"], "parent_completed": c["parent_completed"],
        "cases": {}, "attempts": [], "one_d_attempts": [],
        "statuses": {"NLSP_FEM3C_" + n: "NOT_RUN" for n in STATUS_NAMES},
        "job_calls": {"CCX_production": 0, "Gmsh": 0, "1D_nonlinear_ODE": 0,
            "1D_static": 0, "physical_root_search": 0, "symbolic_derivations": 0},
        "numerical_seconds": 0., "overall": "NOT_RUN", "execution_mode": "EXPLORATORY_NOT_CERTIFIED",
        "admitted": False, "strict_float64_qualification": "PARTIAL",
        "pre_fem_decision_sha256": sha(bundle / "pre_fem_decision.json"),
        "preflight": {"omega1": science["omega1"], "T1": base.dynamic_settings(science)["T1"],
            "control_horizon": .25 * base.dynamic_settings(science)["T1"], "status": "PASS"}}
    summary["statuses"].update(NLSP_FEM3C_SOURCE_PRESERVATION="PASS", NLSP_FEM3C_INPUT_PROTOCOL="PASS")
    save(bundle, item, summary)
    return bundle, item, summary


def recover_saved_case(bundle, case, science, old, mesh, audit, kind):
    """Read one successful native job with an explicit medium/fine preload.

    This is the scoped I/O adapter needed because historical pilot recovery
    fixes the medium source path. All mechanical reconstruction calls below
    are the same existing FEM-2/FEM-3A functions.
    """
    started = time.perf_counter()
    case = Path(case)
    ids, xyz, _, _ = base.base.fem1.mesh_arrays(mesh)
    sta = base.io.read_transient_sta(case / "motion.sta")
    write_json(case / "increments.json", sta)
    accepted = sta["accepted_increments"]
    static_rows = [r for r in accepted if r["step"] == 1]
    dynamic_rows = [r for r in accepted if r["step"] == 2]
    end = base.dynamic_settings(science)["duration"]
    if not static_rows or not dynamic_rows or abs(static_rows[-1]["step_time"] - 1) > 1e-6:
        raise ValueError("Missing full-load STATIC or actual DYNAMIC increments")
    last_sta = dynamic_rows[-1]
    bound = last_sta["time_rounding_bounds"]["step_time"] + 8*np.finfo(float).eps*max(1., end)
    if abs(last_sta["step_time"] - end) > bound:
        raise ValueError("Actual accepted prefix does not reach requested dynamic horizon")
    static_end = static_rows[-1]["total_time"]
    sets = {"ALL_NODES": ids, "LEFT_FIXED": audit["fixed_left_ids"],
            "RIGHT_FIXED": audit["fixed_right_ids"]}
    datdir = case / "dat_fields"
    datdir.mkdir(exist_ok=True)
    dat = {}
    for block in base.io.iter_transient_dat(case / "motion.dat", sets,
            static_end_time=static_end, increments=accepted):
        key = (block["step"], block["increment"], block["set"], block["name"])
        dest = datdir / ("_".join(map(str, key)) + ".npz")
        np.savez_compressed(dest, values=block["values"])
        dat[key] = dest
    quad = base.base.fem1.quadrature_arrays(mesh, science["material"]["rho"])
    fixed = np.r_[audit["fixed_left_ids"], audit["fixed_right_ids"]]
    framesdir = case / "frames"
    framesdir.mkdir(exist_ok=True)
    x = np.linspace(0., 1., 41)
    rows, static_final, kinetic = [], None, []
    max_round, min_det, max_strain = 0., math.inf, 0.
    for frame in base.io.iter_transient_frd(case / "motion.frd", ids,
            static_end_time=static_end, fixed_node_ids=fixed, increments=accepted):
        key = (frame["step"], frame["increment"], "ALL_NODES", "DISP")
        if key not in dat:
            raise ValueError("Missing complete DAT displacement matching FRD frame")
        with np.load(dat[key], allow_pickle=False) as saved:
            U = saved["values"].copy()
        difference = float(np.max(abs(U - frame["fields"]["DISP"])))
        max_round = max(max_round, difference)
        if (difference > 1e-8 or frame["fixed_displacement_max"] > 1e-12
            or (frame.get("fixed_velocity_max") or 0.) > 1e-12):
            raise ValueError("Unchanged DAT/FRD or fixed-face gate failed")
        meta = {k: v for k, v in frame.items() if k != "fields"}
        payload = {"node_ids": ids, "U": U,
                   **{k: v for k, v in frame["fields"].items() if k != "DISP"}}
        dest = framesdir / f"step{frame['step']}_inc{frame['increment']:05d}.npz"
        np.savez_compressed(dest, **payload)
        if frame["step"] == 1:
            static_final = frame, U, meta
            continue
        if frame["step"] != 2:
            raise ValueError("Unexpected native step")
        strain = base.base.fem2_fe_strain_diagnostics(mesh, U)
        if not strain["finite_values"] or strain["minimum_det_deformation_gradient"] <= 0:
            raise ValueError("Nonfinite or inverted actual dynamic deformation")
        min_det = min(min_det, strain["minimum_det_deformation_gradient"])
        max_strain = max(max_strain, strain["green_lagrange_max_abs"])
        displacement = base.base.fem1.nlsp_evaluate_tet10_displacements(U, quad["conn"], quad["N"])
        profile = base.base.fem2_recover_reference_samples(quad["xyz"], displacement,
            quad["weights"], 1., .1, .2, 41)
        velocity = base.base.fem1.nlsp_evaluate_tet10_displacements(
            frame["fields"]["VELO"], quad["conn"], quad["N"])
        vprofile = base.base.fem2_recover_reference_samples(quad["xyz"], velocity,
            quad["weights"], 1., .1, .2, 41)
        kinetic.append(.5*float(np.sum(quad["weights"] * np.sum(velocity*velocity, axis=2))))
        actual = frame["dynamic_time"]
        time_value = actual
        if frame["increment"] == last_sta["increment"]:
            native_bound = frame.get("total_time_rounding_bound", 1e-8) + 8*np.finfo(float).eps*max(1., end)
            if abs(actual-end) > native_bound:
                raise ValueError("Final native output does not contain target-time rounding interval")
            time_value = end
        rows.append({"time": time_value, "printed_dynamic_time": actual,
            "increment": frame["increment"], "fields": base.base.fem2_static_sample(profile, x),
            "translation_velocities": base.base.fem2_static_sample(vprofile, x)[:, :3],
            "frame": dest.relative_to(bundle).as_posix(), "native_metadata": meta})
    if static_final is None or not rows:
        raise ValueError("Missing actual preload or dynamic output")
    if rows[-1]["increment"] != last_sta["increment"] or rows[-1]["time"] != end:
        raise ValueError("Output cadence omitted final accepted state")
    frame, U, meta = static_final
    reference = ROOT / science["source_resume"]["bundle"] / "cases" / science["mesh_level"] / kind
    with np.load(reference / "static_nodal_results.npz", allow_pickle=False) as saved:
        old_arrays = {k: saved[k].copy() for k in saved.files}
    differences = {"node_displacement_max_difference": float(np.max(abs(U-old_arrays["U"]))),
        "stress_difference": float(np.max(abs(frame["fields"]["STRESS"]-old_arrays["S"]))),
        "strain_difference": float(np.max(abs(frame["fields"]["TOSTRAIN"]-old_arrays["E"])))}
    RF = []
    for name in ("LEFT_FIXED", "RIGHT_FIXED"):
        with np.load(dat[(1, frame["increment"], name, "FORC")], allow_pickle=False) as saved:
            RF.append(saved["values"].copy())
    differences["support_RF_max_difference"] = float(np.max(abs(np.vstack(RF)-old_arrays["RF_support_DAT"])))
    gates = science["preload_reproduction"]
    if (differences["node_displacement_max_difference"] > gates["U_absolute"]
        or differences["support_RF_max_difference"] > gates["RF_absolute"]
        or differences["stress_difference"] > gates["relative_S_E"]*np.max(abs(old_arrays["S"]))
        or differences["strain_difference"] > gates["relative_S_E"]*np.max(abs(old_arrays["E"]))):
        raise ValueError("Static preload did not reproduce its own immutable mesh-level reference")
    static_displacement = base.base.fem1.nlsp_evaluate_tet10_displacements(U, quad["conn"], quad["N"])
    static_profile = base.base.fem2_recover_reference_samples(quad["xyz"], static_displacement,
        quad["weights"], 1., .1, .2, 41)
    fields0 = base.base.fem2_static_sample(static_profile, x)
    differences["section_profile_max_difference"] = float(np.max(abs(
        np.asarray(static_profile["fields"])-np.asarray(read_json(reference/"recovered_sections.json")["fields"]))))
    preload = {"status": "PASS", **differences, "source_static_case": reference.relative_to(ROOT).as_posix(),
        "source_static_nodal_sha256": sha(reference / "static_nodal_results.npz"),
        "static_end_total_time": static_end, "actual_metadata": meta,
        "velocity_initialization": "explicit zero IC and verified source-zeroed static-to-dynamic state; STATIC is not native DYNAMIC t0"}
    write_json(case / "preload_transfer.json", preload)
    np.savez_compressed(case / "initial_sections.npz", x=x, fields=fields0,
        source_static_end_time=static_end, not_a_native_dynamic_zero_frame=np.array(True))
    np.savez_compressed(case / "section_history.npz", time=np.array([r["time"] for r in rows]),
        printed_dynamic_time=np.array([r["printed_dynamic_time"] for r in rows]), x=x,
        fields=np.stack([r["fields"] for r in rows]),
        translation_velocities=np.stack([r["translation_velocities"] for r in rows]),
        increments=np.array([r["increment"] for r in rows]))
    write_json(case / "frame_metadata.json", [{k: v for k, v in r.items()
        if k not in ("fields", "translation_velocities")} for r in rows])
    np.savez_compressed(case / "independent_kinetic_energy.npz",
        time=np.array([r["time"] for r in rows]), kinetic_energy=np.asarray(kinetic))
    energy = base.io.parse_transient_dat_energies(case / "motion.dat", static_end_time=static_end,
        increments=accepted, element_set="SOLID")
    stdout_energy = base.io.parse_transient_stdout_energies(case / "motion.stdout.txt",
        static_end_time=static_end, increments=accepted)
    write_json(case / "energy.json", energy)
    write_json(case / "stdout_energy.json", stdout_energy)
    initial = [r for r in energy["records"] if r["step"] == 1 and abs(r["total_time"]-static_end) < 1e-7]
    dynamic_energy = [r for r in energy["records"] if r["step"] == 2]
    E0 = initial[-1]["internal_energy"] if initial else None
    stdout_dynamic = [r for r in stdout_energy["records"] if r.get("step") == 2]
    work_complete = bool(stdout_dynamic) and all("external_work" in r and "damping_work" in r for r in stdout_dynamic)
    external = max(abs(r["external_work"]) for r in stdout_dynamic) if work_complete else None
    damping = max(abs(r["damping_work"]) for r in stdout_dynamic) if work_complete else None
    if external is None or damping is None or external != 0 or damping != 0:
        raise ValueError("Missing or nonzero actual external/damping work after release")
    if rows[0]["fields"][20, 1] >= fields0[20, 1]:
        raise ValueError("Early free movement is not restoring")
    record = {"status": "PASS", "preload_transfer": preload,
        "dynamic_time_start": rows[0]["time"], "dynamic_time_end": rows[-1]["time"],
        "dynamic_output_frames": len(rows), "static_increments": len(static_rows),
        "dynamic_increments": len(dynamic_rows), "accepted_increments": len(accepted),
        "cutbacks": sta["reported_cutbacks"], "output_frequency": science["dynamic"]["output_frequency"],
        "final_native_frame_reached": True, "max_DAT_FRD_displacement_difference": max_round,
        "initial_midspan_w": fields0[20, 1], "first_midspan_w": rows[0]["fields"][20, 1],
        "final_midspan_w": rows[-1]["fields"][20, 1],
        "first_midspan_w_velocity": rows[0]["translation_velocities"][20, 1],
        "maximum_native_external_work_after_release": external,
        "maximum_native_damping_work_after_release": damping,
        "native_initial_internal_energy": E0, "energy_status": "PARTIAL",
        "max_native_relative_mechanical_energy_drift": max((abs(r["mechanical_energy"]/E0-1)
            for r in dynamic_energy), default=None) if E0 else None,
        "all_saved_frames_strain_diagnostics": {"finite_values": True,
            "minimum_det_deformation_gradient": min_det, "max_abs_green_lagrange_strain": max_strain,
            "sampling": "existing 14-point element quadrature at every actual saved native frame"},
        "final_strain_diagnostics": strain, "independent_internal_energy": "NOT_RUN",
        "recovery_seconds": time.perf_counter()-started}
    resume.qualify_native_energy(case, record)
    write_json(case / "recovery.json", record)
    return record


def audit_saved_case(case, science, old, mesh, audit, kind, record):
    lines = (case / "motion.stdout.txt").read_text(encoding="utf8").splitlines()
    warnings = [line for line in lines if "*WARNING" in line.upper()]
    if warnings:
        raise ValueError("Unexplained solver warnings: " + str(warnings))
    sta = read_json(case / "increments.json")
    static = [r for r in sta["accepted_increments"] if r["step"] == 1]
    inc = static[-1]["increment"]
    ids, xyz, gravity, volume = base.base.consistent_gravity_loads(mesh,
        science["material"]["rho"], science["g"])
    index = {int(node): i for i, node in enumerate(ids)}
    with np.load(case / "frames" / f"step1_inc{inc:05d}.npz", allow_pickle=False) as saved:
        U = saved["U"].copy()
    supports, positions, resultants = [], [], {}
    for name, key in (("LEFT_FIXED", "fixed_left_ids"), ("RIGHT_FIXED", "fixed_right_ids")):
        rows = np.array([index[int(node)] for node in audit[key]])
        with np.load(case / "dat_fields" / f"1_{inc}_{name}_FORC.npz", allow_pickle=False) as saved:
            RF = saved["values"].copy()
        support = RF-gravity[rows]
        current = xyz[rows]+U[rows] if kind == "nonlinear" else xyz[rows]
        supports.append(support)
        positions.append(current)
        resultants[name] = {"force": support.sum(axis=0), "moment_about_face_centroid":
            np.cross(current-current.mean(axis=0), support).sum(axis=0)}
    applied = gravity.sum(axis=0)
    scale = np.linalg.norm(applied)
    current = xyz+U if kind == "nonlinear" else xyz
    force_imbalance = np.linalg.norm(np.vstack(supports).sum(axis=0)+applied)/scale
    moment_imbalance = np.linalg.norm(np.cross(np.vstack(positions), np.vstack(supports)).sum(axis=0)
        + np.cross(current, gravity).sum(axis=0))/scale
    gate = old["science_config"]["gates"]["equilibrium_relative"]
    if force_imbalance > gate or moment_imbalance > gate:
        raise ValueError("Unchanged independent preload equilibrium gate failed")
    result = {"status": "PASS", "warning_lines": warnings,
        "total_applied_force": applied, "reference_volume": volume,
        "support_resultants": resultants, "preload_force_imbalance_relative": float(force_imbalance),
        "preload_moment_imbalance_relative": float(moment_imbalance), "equilibrium_gate_unchanged": gate,
        "bodyload_recovered_independently": True, "source_mesh_level": science["mesh_level"],
        "static_end_total_time": static[-1]["total_time"],
        "release_actual_external_work": record["maximum_native_external_work_after_release"],
        "release_actual_damping_work": record["maximum_native_damping_work_after_release"],
        "zero_physical_initial_velocities": "explicit native IC plus installed verified static-to-dynamic zeroing; first actual frame at positive time"}
    write_json(case / "continuation_audit.json", result)
    record["continuation_audit"] = result
    return result


def run_case(bundle, item, summary, stage, kind):
    """Execute at most once; parser repair only rereads a finished attempt."""
    bundle = Path(bundle)
    c = validate_config(item["validation_config"])
    if (stage, kind) not in CASE_ORDER:
        raise ValueError("Unauthorized production case")
    if summary.get("hard_stop"):
        return False
    if sha(bundle / "pre_fem_decision.json") != summary["pre_fem_decision_sha256"]:
        raise ValueError("Pre-FEM policies changed after authorization")
    name = stage + "/" + kind
    previous_record = summary["cases"].get(name)
    if previous_record and previous_record["status"] == "PASS":
        return True
    if previous_record and previous_record["status"] != "OUTPUT_RECOVERY_PENDING":
        return False
    ordinal = CASE_ORDER.index((stage, kind))
    for prior_stage, prior_kind in CASE_ORDER[:ordinal]:
        prior = summary["cases"].get(prior_stage + "/" + prior_kind, {})
        if prior.get("status") != "PASS" or prior.get("continuation_audit", {}).get("status") != "PASS":
            raise ValueError("Sequential production gate missing: " + prior_stage + "/" + prior_kind)
    if stage == FULL_STAGE:
        path = bundle / "full_period_decision.json"
        if (not path.is_file() or sha(path) != summary.get("full_period_decision_sha256")
            or not read_json(path).get("full_period_3D_authorized_by_actual_gates")):
            raise ValueError("Full-period 3D requires frozen successful numerical-robustness decision")
    science = case_science(item["config"], c, stage)
    old, source, mesh, audit = source_for_level(item["config"], science["mesh_level"])
    case = bundle / "cases" / stage / kind
    if not previous_record:
        if summary["job_calls"]["CCX_production"] >= 6:
            raise RuntimeError("Six authorized production attempts exhausted")
        remaining = c["numerical_budget_seconds"]-summary["numerical_seconds"]
        if remaining <= 0:
            summary.update(hard_stop=True, overall="PARTIAL", stop_reason="Bounded numerical budget exhausted")
            save(bundle, item, summary)
            return False
        if case.exists() and any((case / filename).exists() for filename in
                ("job.json", "motion.stdout.txt", "motion.dat", "motion.frd", "motion.sta")):
            raise RuntimeError("Unledgered native outputs cannot be overwritten or retried")
        case.mkdir(parents=True, exist_ok=True)
        actual_helpers = {name: sha(ROOT / name) for name in item["helper_sha256"]}
        frozen = (Path(base.__file__), Path(base.base.__file__), Path(base.base.fem1.__file__),
            Path(base.one.dynamics.__file__), Path(base.one.rod.__file__))
        for path in frozen:
            helper_name = path.relative_to(ROOT).as_posix()
            if actual_helpers[helper_name] != item["helper_sha256"][helper_name]:
                raise ValueError("Frozen physics/generator helper changed before production: " + helper_name)
        phase = bundle / "execution_code" / ("attempt_" + str(len(summary["attempts"])+1))
        suffix = 1
        while phase.exists():
            suffix += 1
            phase = bundle / "execution_code" / ("attempt_" + str(len(summary["attempts"])+1)
                + "_prelaunch_" + str(suffix))
        phase.mkdir()
        for helper_name in actual_helpers:
            shutil.copyfile(ROOT / helper_name, phase / Path(helper_name).name)
        gate = write_case_input(case / "motion.inp", science, old, source, mesh, audit, kind == "nonlinear")
        preview = bundle / "input_gate" / stage / (kind + ".inp")
        if preview.is_file():
            before = preview.read_text(encoding="utf8")
            after = (case / "motion.inp").read_text(encoding="utf8")
            include_evidence = {}
            def physical_text(text, input_path, label):
                def normalize_include(match):
                    raw = match.group(1).strip().strip('"')
                    resolved = (input_path.parent / raw).resolve()
                    expected = (source / "solid_mesh.inp").resolve()
                    if resolved != expected or not resolved.is_file() or sha(resolved) != sha(expected):
                        raise ValueError("Preflight/production INCLUDE does not resolve to the frozen source mesh")
                    include_evidence[label] = {"written": raw, "resolved": str(resolved), "sha256": sha(resolved)}
                    return "*INCLUDE, INPUT=<FROZEN_SOURCE_MESH>"
                text = re.sub(r"(?im)^\*INCLUDE,\s*INPUT=([^\n\r]+)", normalize_include, text)
                static, dynamic = text.split("*END STEP\n", 1)
                dynamic = re.sub(r"\*EL FILE[^\n]*\nS,E,ENER\n", "", dynamic)
                return static + "*END STEP\n" + dynamic
            if physical_text(before, preview, "preflight") != physical_text(after, case / "motion.inp", "production"):
                raise ValueError("Actual production deck differs physically from saved preflight")
            write_json(case / "input_adapter_revision.json", {
                "preflight_input_sha256": sha(preview), "actual_input_sha256": sha(case / "motion.inp"),
                "different": before != after, "physical_preflight_equivalence": True,
                "resolved_include_evidence": include_evidence,
                "permitted_output_only_change": "Repeat inherited S,E,ENER EL FILE in DYNAMIC at declared frequency; STATIC unchanged"})
        write_json(case / "input_contract.json", gate)
        write_json(case / "science_config.json", science)
        summary["attempts"].append({"case": name, "ordinal": len(summary["attempts"])+1,
            "status": "STARTED", "authorization_id": AUTHORIZATION,
            "input_sha256": sha(case / "motion.inp"), "source_mesh_include_sha256": sha(source / "solid_mesh.inp"),
            "actual_execution_helper_sha256": actual_helpers,
            "pre_fem_decision_sha256": summary["pre_fem_decision_sha256"], "automatic_retry": False})
        summary["cases"][name] = {"status": "STARTED", "mesh_level": science["mesh_level"]}
        env = dict(os.environ)
        env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMBER_OF_CPUS="1")
        command = [old["science_config"]["ccx_exe"], "motion"]
        write_json(case / "execution_environment.json", {"command": command, "cwd": str(case),
            "threads": 1, "PATH_unchanged": True,
            "thread_environment": {k: env[k] for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBER_OF_CPUS")}})
        summary["job_calls"]["CCX_production"] += 1
        save(bundle, item, summary)
        started = time.perf_counter()
        try:
            print("FEM3C production STATIC+DYNAMIC: " + name, flush=True)
            result, stats = base.base.fem1.run_job(command, case,
                min(c["job_timeout_seconds"], remaining), c["job_memory_limit_bytes"], case / "motion", env)
            write_json(case / "job.json", stats)
            summary["cases"][name]["job"] = stats
            stdout = (case / "motion.stdout.txt").read_text(encoding="utf8", errors="replace")
            if result.returncode != 0 or stats["failure"] or "JOB FINISHED" not in stdout.upper() or "*ERROR" in stdout.upper():
                raise RuntimeError("Actual solver failure: " + str(stats))
            summary["cases"][name]["status"] = "OUTPUT_RECOVERY_PENDING"
            summary["attempts"][-1]["status"] = "SOLVER_FINISHED"
        except Exception as error:
            summary["cases"][name].update(status="FAIL", failure=str(error))
            summary["attempts"][-1].update(status="FAIL", failure=str(error))
            summary.update(hard_stop=True, overall="BLOCKED_BY_SOLVER", stop_reason=str(error))
        finally:
            summary["numerical_seconds"] += time.perf_counter()-started
            save(bundle, item, summary)
        if summary["cases"][name]["status"] == "FAIL":
            return False
    started = time.perf_counter()
    try:
        job = summary["cases"][name]["job"]
        record = recover_saved_case(bundle, case, science, old, mesh, audit, kind)
        record["job"] = job
        audit_saved_case(case, science, old, mesh, audit, kind, record)
        summary["cases"][name] = record
        next(row for row in summary["attempts"] if row["case"] == name)["status"] = "PASS"
        summary["statuses"]["NLSP_FEM3C_ENERGY_DIAGNOSTICS"] = "PARTIAL"
        summary["overall"] = "PARTIAL"
        pair_passed = all(summary["cases"].get(stage + "/" + k, {}).get("status") == "PASS"
            for k in ("linear", "nonlinear"))
        if pair_passed:
            status = {"medium_refined_time": "TEMPORAL_CONTROL", "fine_refined_time": "SPATIAL_CONTROL",
                FULL_STAGE: "FULL_PERIOD_3D"}[stage]
            summary["statuses"]["NLSP_FEM3C_" + status] = "PASS"
        return True
    except Exception as error:
        summary["cases"][name].update(status="OUTPUT_RECOVERY_PENDING", recovery_failure=str(error))
        summary["overall"] = "PARTIAL"
        summary["stop_reason"] = "Finished native output needs verified read-only recovery: " + str(error)
        return False
    finally:
        summary["numerical_seconds"] += time.perf_counter()-started
        save(bundle, item, summary)


def run_controls(bundle, item, summary, through_case=None):
    for stage in CONTROL_STAGES:
        for kind in ("linear", "nonlinear"):
            if not run_case(bundle, item, summary, stage, kind):
                return summary
            if through_case == stage + "/" + kind:
                return summary
    return summary


def complete_robustness(bundle, item, summary):
    from scripts.lib import nlsp_fem3c_diagnostics
    if not controls_ready(summary):
        return None
    record_postprocessing_phase(bundle, item, summary, "quarter_period_robustness")
    started = time.perf_counter()
    result = nlsp_fem3c_diagnostics.robustness(bundle, item["validation_config"])
    summary["numerical_seconds"] += time.perf_counter()-started
    summary["statuses"]["NLSP_FEM3C_ROBUSTNESS"] = result["numerical_robustness_status"]
    summary["robustness_gate"] = result["full_period_numerical_robustness_gate"]
    summary["robustness_artifact"] = "robustness_comparison.json"
    save(bundle, item, summary)
    return result


def controls_ready(summary):
    return all(summary["cases"].get(stage + "/" + kind, {}).get("status") == "PASS"
        and summary["cases"][stage + "/" + kind].get("continuation_audit", {}).get("status") == "PASS"
        for stage in CONTROL_STAGES for kind in ("linear", "nonlinear"))


def one_d_stage_allowed(summary):
    """No full-period integration before the bounded C1 outcome is known."""
    if summary.get("hard_stop") or any(case.get("status") == "STARTED"
            for case in summary.get("cases", {}).values()):
        return False
    return controls_ready(summary) or bool(summary.get("c1_terminated_without_solver_failure"))


def record_postprocessing_phase(bundle, item, summary, phase):
    """Capture the actual implementation before analysis or permitted 1D work."""
    bundle = Path(bundle)
    files = [Path(__file__), Path(__file__).with_name("nlsp_fem3c_diagnostics.py"),
        Path(__file__).with_name("nlsp_fem3c_1d.py"), Path(previous.__file__),
        Path(__file__).with_name("nlsp_fem3b_diagnostics.py")]
    rows = summary.setdefault("postprocessing_phases", [])
    ordinal = len(rows)+1
    destination = bundle / "execution_code" / f"postprocessing_{ordinal:03d}_{phase}"
    destination.mkdir()
    hashes = {}
    for path in files:
        hashes[path.relative_to(ROOT).as_posix()] = sha(path)
        shutil.copyfile(path, destination / path.name)
    rows.append({"ordinal": ordinal, "phase": phase, "actual_helper_sha256": hashes,
        "scientific_calls": "at most two authorized nonlinear 1D calls" if phase == "full_period_1d" else 0,
        "no_native_job_from_postprocessing": True})
    save(bundle, item, summary)


def plot_saved(bundle):
    """Render only figures whose underlying comparison arrays already exist."""
    from scripts.lib import nlsp_fem3c_diagnostics
    bundle = Path(bundle)
    paths = []
    if (bundle / "robustness_comparison.npz").exists():
        paths.extend(nlsp_fem3c_diagnostics.plot_robustness(bundle))
    if (bundle / "full_period_comparison.npz").exists():
        paths.extend(nlsp_fem3c_diagnostics.plot_full_period(bundle))
    return [str(path) for path in paths]


def complete_validation(bundle, item, summary):
    """Close only the actually completed, bounded program, keeping qualifications."""
    from scripts.lib import nlsp_fem3c_diagnostics
    bundle = Path(bundle)
    if summary.get("hard_stop") or not controls_ready(summary):
        return summary
    full_cases = all(summary["cases"].get(FULL_STAGE + "/" + kind, {}).get("status") == "PASS"
        and summary["cases"][FULL_STAGE + "/" + kind].get("continuation_audit", {}).get("status") == "PASS"
        for kind in ("linear", "nonlinear"))
    if full_cases and summary.get("one_d_completed"):
        record_postprocessing_phase(bundle, item, summary, "full_period_comparison")
        started = time.perf_counter()
        result = nlsp_fem3c_diagnostics.full_period(bundle, item, summary)
        summary["numerical_seconds"] += time.perf_counter()-started
        summary["full_period_comparison_artifact"] = "full_period_comparison.json"
        summary["statuses"]["NLSP_FEM3C_SEVEN_FIELD_RECOVERY"] = "PASS"
        summary["statuses"]["NLSP_FEM3C_ENERGY_DIAGNOSTICS"] = "PARTIAL"
        summary["full_period_1d_spatial_qualification"] = result["one_d_spatial_status"]
        summary["full_period_illustrative_only"] = True
        summary["statuses"]["NLSP_FEM3C_VERIFICATION_SUMMARY"] = "PASS" if summary.get("robustness_gate") else "PARTIAL"
        summary["overall"] = ("STRAIGHT_ROD_NONLINEAR_3D_FEM_VERIFICATION_COMPLETE_WITH_QUALIFICATIONS"
            if summary.get("robustness_gate") else "PARTIAL")
        summary["completed"] = True
    elif summary.get("one_d_completed") and summary.get("full_period_3D_allowed") is False:
        summary.update(overall="PARTIAL", completed=True)
        summary["statuses"]["NLSP_FEM3C_VERIFICATION_SUMMARY"] = "PARTIAL"
        summary["scientific_qualification"] = summary["full_period_stop_reason"]
    if summary.get("completed"):
        summary["universal_nonlinear_validation"] = False
        summary["experimental_validation"] = False
        save(bundle, item, summary)
        record_postprocessing_phase(bundle, item, summary, "saved_figures")
        summary["figures"] = plot_saved(bundle)
        save(bundle, item, summary)
    return summary


def freeze_full_period(bundle, item, summary):
    """Select the already declared illustration policy before full-period data."""
    bundle = Path(bundle)
    path = bundle / "full_period_decision.json"
    if path.exists():
        if sha(path) != summary.get("full_period_decision_sha256"):
            raise ValueError("Frozen full-period decision changed")
        return read_json(path)["full_period_3D_authorized_by_actual_gates"]
    result = read_json(bundle / "robustness_comparison.json")
    all_cases = all(summary["cases"].get(s + "/" + k, {}).get("status") == "PASS"
        and summary["cases"][s + "/" + k].get("continuation_audit", {}).get("status") == "PASS"
        for s in CONTROL_STAGES for k in ("linear", "nonlinear"))
    old_estimate = read_json(bundle / "pre_fem_decision.json")["planning_CCX_estimates_seconds"][FULL_STAGE]
    observed_refined = sum(summary["cases"][CONTROL_STAGES[0] + "/" + k]["job"]["seconds"]
        for k in ("linear", "nonlinear"))
    estimate = max(old_estimate, 2*observed_refined)
    remaining = item["validation_config"]["numerical_budget_seconds"]-summary["numerical_seconds"]
    resource_ok = estimate + 1800 <= remaining and estimate/2 <= item["validation_config"]["job_timeout_seconds"]
    allowed = bool(result["full_period_numerical_robustness_gate"] and all_cases and resource_ok and not summary.get("hard_stop"))
    decision = {"before_full_period_results": True, "full_period_3D_authorized_by_actual_gates": allowed,
        "robustness_comparison_sha256": sha(bundle / "robustness_comparison.json"),
        "control_solver_preload_release_recovery_gates": all_cases,
        "actual_numerical_robustness_gate": result["full_period_numerical_robustness_gate"],
        "resource_preflight": resource_ok, "planning_CCX_seconds": estimate,
        "remaining_numerical_budget_seconds": remaining,
        "illustration_policy": item["validation_config"]["stages"][FULL_STAGE],
        "full_period_convergence_not_certified_by_quarter_period_control": True,
        "no_choice_by_future_1D_3D_agreement": True}
    write_json(path, decision)
    summary["full_period_decision_sha256"] = sha(path)
    summary["full_period_3D_allowed"] = allowed
    if not allowed:
        summary["full_period_stop_reason"] = "Quarter-period numerical robustness or resource gate is not satisfied; no full-period 3D jobs"
    save(bundle, item, summary)
    return allowed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--run-controls", action="store_true")
    mode.add_argument("--run-full-period", action="store_true")
    mode.add_argument("--run-1d", action="store_true")
    mode.add_argument("--compute", action="store_true")
    mode.add_argument("--report-only", type=Path)
    mode.add_argument("--plot-only", type=Path)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--through-case", choices=[s + "/" + k for s, k in CASE_ORDER])
    args = parser.parse_args(argv)
    if args.report_only or args.plot_only:
        bundle = args.report_only or args.plot_only
        summary = validate_cache(bundle)
        if args.plot_only:
            plot_saved(bundle)
            save(Path(bundle), read_json(Path(bundle) / "provenance.json"), summary)
        print(json.dumps({"bundle": str(bundle), "overall": summary["overall"],
            "statuses": summary["statuses"], "new_scientific_calls": 0}, indent=2))
        return summary
    bundle, item, summary = prepare_stage(args.config)
    before = dict(summary["job_calls"])
    if summary.get("completed") or summary.get("hard_stop") or args.preflight:
        print(json.dumps({"bundle": str(bundle), "overall": summary["overall"], "new_scientific_calls": 0}, indent=2))
        return summary
    if args.run_controls or args.compute:
        run_controls(bundle, item, summary, args.through_case)
        complete_robustness(bundle, item, summary)
    if args.run_1d and not one_d_stage_allowed(summary):
        raise ValueError("Full-period 1D requires the bounded C1 outcome; no integration before controls or during native execution")
    if (args.run_1d or args.compute or args.run_full_period) and one_d_stage_allowed(summary):
        if not (bundle / "robustness_comparison.json").exists() and controls_ready(summary):
            complete_robustness(bundle, item, summary)
        record_postprocessing_phase(bundle, item, summary, "full_period_1d")
        from scripts.lib import nlsp_fem3c_1d
        nlsp_fem3c_1d.run_full_period(bundle, item, summary)
        if summary.get("one_d_completed"):
            summary["statuses"]["NLSP_FEM3C_FULL_PERIOD_1D"] = "PASS"
            save(bundle, item, summary)
    if (args.run_full_period or args.compute) and not summary.get("hard_stop"):
        if (bundle / "robustness_comparison.json").exists() and freeze_full_period(bundle, item, summary):
            for kind in ("linear", "nonlinear"):
                if not run_case(bundle, item, summary, FULL_STAGE, kind) or args.through_case == FULL_STAGE + "/" + kind:
                    break
    if args.compute or args.run_full_period:
        complete_validation(bundle, item, summary)
    print(json.dumps({"bundle": str(bundle), "overall": summary["overall"],
        "statuses": summary["statuses"], "job_calls": summary["job_calls"],
        "new_scientific_calls": {k: v-before.get(k, 0) for k, v in summary["job_calls"].items()}}, indent=2))
    return summary
