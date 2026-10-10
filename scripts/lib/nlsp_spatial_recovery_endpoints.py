"""Four saved endpoint recoveries per existing mesh; no scientific solves.

The 41-section baseline is read from the native program's immutable outputs.
Only the same four saved nodal states are reconstructed with 81 sections.
"""
from __future__ import annotations

import csv
import shutil
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from scripts.lib import nlsp_spatial_fem_protocol as protocol
from scripts.lib import nlsp_spatial_native_program as native

ROOT = Path(__file__).resolve().parents[2]
FIELDS = protocol.FIELD_ORDER


def difference_metrics(first, second, x, *, physical_orientation=True):
    difference = np.asarray(second)-np.asarray(first)
    if difference.shape != (len(x), 7) or not np.isfinite(difference).all():
        raise ValueError("Seven finite fields on the identical material grid required")
    result = {}
    for column, name in enumerate(FIELDS):
        values = difference[:, column]; index = int(np.argmax(abs(values)))
        result[name] = {"absolute_max": float(abs(values[index])),
            "L2_on_41_point_trapezoidal_grid": float(np.sqrt(np.trapezoid(values*values, x))),
            "signed_81_minus_41_at_max": float(values[index]), "maximum_x": float(x[index]),
            "sampled_maximum_not_continuous_supremum": True}
    if physical_orientation:
        first_R = Rotation.from_rotvec(np.asarray(first)[:, 3:6]*[1., -1., 1.]).as_matrix()
        second_R = Rotation.from_rotvec(np.asarray(second)[:, 3:6]*[1., -1., 1.]).as_matrix()
        angles = np.linalg.norm(Rotation.from_matrix(first_R.swapaxes(1, 2)@second_R).as_rotvec(), axis=1)
        result["physical_orientation"] = {"maximum_principal_angle": float(angles.max()),
            "L2_principal_angle_on_41_point_grid": float(np.sqrt(np.trapezoid(angles*angles, x))),
            "definition": "principal angle between exp(a41) and exp(a81) after identical existing spatial cubic transfer; a=(Phi,-psi,theta)"}
    else:
        result["rotation_qualification"] = "Canonical coordinate differences only; adding/subtracting these is not asserted to be an invariant rotation operation"
    return result


def compare_endpoint_fields(fields41, fields81, x):
    metrics = {"full": {}, "correction": {}}
    for kind in ("linear", "nonlinear"):
        for state in ("static", "final"):
            key = kind+"_"+state
            metrics["full"][key] = difference_metrics(fields41[key], fields81[key], x)
    delta41 = {}; delta81 = {}
    for state in ("static", "final"):
        delta41[state] = fields41["nonlinear_"+state]-fields41["linear_"+state]
        delta81[state] = fields81["nonlinear_"+state]-fields81["linear_"+state]
        metrics["correction"][state] = difference_metrics(delta41[state], delta81[state], x, physical_orientation=False)
    evolution41 = delta41["final"]-delta41["static"]
    evolution81 = delta81["final"]-delta81["static"]
    metrics["evolution"] = difference_metrics(evolution41, evolution81, x, physical_orientation=False)
    return metrics, {"delta41_static": delta41["static"], "delta41_final": delta41["final"],
        "delta81_static": delta81["static"], "delta81_final": delta81["final"],
        "evolution41": evolution41, "evolution81": evolution81}


def _record(path, sources):
    sources[path.resolve().as_posix()] = native.sha(path)


def _verify_source(filename, digest):
    path = Path(filename)
    if path.is_file() and native.sha(path) == digest:
        return
    # Finalizing this same scientific attempt may add read-only diagnostics to
    # its comparison JSON. The original context used for the ratio is then kept
    # verbatim in the parent's explicit first-processing archive. This exception
    # applies only to that derived context, never to nodal fields or histories.
    if path.name == "one_d_three_d_comparison.json" and path.parent.name in ("comparison_medium", "comparison_fine"):
        archived = path.parent/"first_processing_evidence"/path.name
        if archived.is_file() and native.sha(archived) == digest:
            return
    raise ValueError("Saved endpoint recovery source changed")


def audit_recovery_endpoints(bundle, mesh_level="medium"):
    """Postprocess exactly L/NL x STATIC/final-DYNAMIC saved states once."""
    start = time.perf_counter(); bundle = Path(bundle)
    if mesh_level not in ("medium", "fine"):
        raise ValueError("Only the two existing meshes are within audit scope")
    output = bundle/"recovery_sensitivity_endpoints"/mesh_level
    if (output/"manifest.json").exists():
        manifest = protocol.static.read_json(output/"manifest.json")
        for filename, digest in manifest["artifact_hashes"].items():
            if native.sha(output/filename) != digest:
                raise ValueError("Endpoint recovery artifact changed")
        for filename, digest in manifest["source_hashes"].items():
            _verify_source(filename, digest)
        return protocol.static.read_json(output/"metrics.json")
    if output.exists() and any(output.iterdir()):
        raise ValueError("Existing unmanifested endpoint diagnostics must be preserved")
    config = protocol.static.read_json(bundle/"config.json")
    if config["geometry"] != {"L": 1., "b": .2, "h": .1}:
        raise ValueError("Frozen reference geometry changed")
    sources = {}; _record(bundle/"config.json", sources)
    for kind in ("linear", "nonlinear"):
        case = bundle/"FEM"/(mesh_level+"_"+kind)
        if protocol.static.read_json(case/"native_attempt.json")["status"] != "PASS":
            raise ValueError("Endpoint recovery requires a completed valid native pair")
    meshdir = ROOT/config["source_bundles"]["FEM1"]["path"]/"meshes"/mesh_level
    source_manifest = ROOT/config["source_bundles"]["FEM1"]["path"]/"manifest.json"
    if native.sha(source_manifest) != config["source_bundles"]["FEM1"]["manifest_sha256"]:
        raise ValueError("Historical source manifest changed")
    _record(source_manifest, sources); _record(meshdir/"solid_mesh.inp", sources)
    source_manifest_data = protocol.static.read_json(source_manifest)
    source_hashes = source_manifest_data.get("artifact_hashes", source_manifest_data.get("artifacts", {}))
    if native.sha(meshdir/"solid_mesh.inp") != source_hashes["meshes/"+mesh_level+"/solid_mesh.inp"]:
        raise ValueError("Saved mesh include hash mismatch")
    mesh = protocol.fem1.single.read_gmsh_inp_mesh_data(meshdir/"solid_mesh.inp")
    quadrature = protocol.fem1.quadrature_arrays(mesh, config["material"]["rho"])
    fields41 = {}; fields81 = {}; payload = {}; diagnostics = {}; endpoint_metadata = {}
    output.mkdir(parents=True)
    phase = output/"execution_code"; phase.mkdir()
    code_sources = [Path(__file__), Path(protocol.__file__), Path(protocol.static.__file__), Path(protocol.fem1.__file__)]
    code_hashes = {}
    for source in code_sources:
        shutil.copyfile(source, phase/source.name); code_hashes[source.as_posix()] = native.sha(source)
    common_x = None
    recovery_calls = 0
    for kind in ("linear", "nonlinear"):
        case = bundle/"FEM"/(mesh_level+"_"+kind)
        attempt = protocol.static.read_json(case/"native_attempt.json")
        if attempt["status"] != "PASS":
            raise ValueError("Endpoint recovery requires a completed valid native pair")
        for filename in ("native_attempt.json", "increments.json", "initial_sections.npz", "section_history.npz"):
            _record(case/filename, sources)
        accepted = protocol.static.read_json(case/"increments.json")["accepted_increments"]
        endpoint_rows = {"static": [row for row in accepted if row["step"] == 1][-1],
                         "final": [row for row in accepted if row["step"] == 2][-1]}
        with np.load(case/"initial_sections.npz", allow_pickle=False) as saved:
            x = saved["x"].copy(); fields41[kind+"_static"] = saved["fields"].copy()
            if not saved["not_a_native_dynamic_zero_frame"]:
                raise ValueError("STATIC anchor must remain distinct from actual dynamic output")
        with np.load(case/"section_history.npz", allow_pickle=False) as saved:
            if not np.array_equal(saved["x"], x) or saved["increments"][-1] != endpoint_rows["final"]["increment"]:
                raise ValueError("Saved section endpoint differs from actual final native increment")
            fields41[kind+"_final"] = saved["fields"][-1].copy()
            endpoint_metadata[kind] = {"dynamic_time": float(saved["time"][-1]),
                "printed_dynamic_time": float(saved["printed_dynamic_time"][-1]),
                "static_source": "confirmed STATIC preload, not native DYNAMIC t=0", "increments": endpoint_rows}
        if common_x is None:
            common_x = x
        elif not np.array_equal(common_x, x):
            raise ValueError("Linear/nonlinear endpoint material grids differ")
        for state, row in endpoint_rows.items():
            key = kind+"_"+state
            frame = case/"frames"/f"step{row['step']}_inc{row['increment']:05d}.npz"
            _record(frame, sources)
            with np.load(frame, allow_pickle=False) as saved:
                U = saved["U"].copy()
                if not np.array_equal(saved["node_ids"], quadrature["node_ids"]):
                    raise ValueError("Native displacement ordering differs from source mesh")
            profile = protocol.recover_spatial_sections(mesh, U, section_count=81, quadrature=quadrature)
            recovery_calls += 1
            fields81[key] = protocol.static.fem2_static_sample(profile, common_x)
            payload[key+"_raw_x81"] = profile["x"]
            payload[key+"_raw_fields81"] = profile["fields"]
            payload[key+"_raw_rotation_matrices81"] = profile["raw_rotation_matrices"]
            protocol.static.write_json(output/(key+"_raw81.json"), profile)
            rows = profile["section_rows"]
            diagnostics[key] = {"raw_sections": len(rows), "rank_min": min(row["fit_rank"] for row in rows),
                "condition_max": max(row["fit_condition"] for row in rows),
                "residual_mass_L2": profile["section_residual_mass_L2"],
                "residual_relative_L2": profile["section_residual_relative_L2"],
                "raw_c_eff_absolute_max": float(np.max(abs(profile["fields"][1:-1, 6]))),
                "c_eff_qualification": "effective thickness stretch proxy, not M-H coordinate",
                "smallest_reference_mass": min(row["reference_mass"] for row in rows),
                "minimum_quadrature_samples": min(row["samples"] for row in rows)}
    if endpoint_metadata["linear"]["dynamic_time"] != endpoint_metadata["nonlinear"]["dynamic_time"]:
        raise ValueError("Actual nonlinear and linear final physical times differ")
    metrics, corrections = compare_endpoint_fields(fields41, fields81, common_x)
    model = protocol.static.read_json(bundle/("comparison_"+mesh_level)/"one_d_three_d_comparison.json")
    _record(bundle/("comparison_"+mesh_level)/"one_d_three_d_comparison.json", sources)
    v_discrepancy = model["metrics"]["evolution"]["v"]["absolute_max"]
    metrics.update(mesh_level=mesh_level, status="DIAGNOSTIC_COMPLETE_WITH_QUALIFICATIONS",
        source_hashes=sources, source_code_sha256=code_hashes, actual_endpoint_metadata=endpoint_metadata,
        raw81_diagnostics=diagnostics, recovery_calls_81=recovery_calls,
        scientific_calls={"CCX": 0, "Gmsh": 0, "Radau": 0, "static_Newton": 0, "eigen": 0, "BVP": 0},
        fixed_grid_policy="The original same 41 material x points; unsmoothed existing cubic spatial transfer only",
        L2_policy="sqrt(trapezoid(error_squared,x)) on the identical41 grid; separate diagnostic, no historical gate redefinition",
        v_model_discrepancy_full_horizon_max=float(v_discrepancy),
        v_endpoint_evolution_recovery_change_over_full_horizon_model_discrepancy=
            metrics["evolution"]["v"]["absolute_max"]/v_discrepancy if v_discrepancy else None,
        limitations="Endpoint-only observed recovery sensitivity, not a full-horizon maximum or strict error bound; no smoothing; c_eff remains a qualified proxy",
        seconds=time.perf_counter()-start)
    for filename, digest in sources.items():
        if native.sha(filename) != digest:
            raise ValueError("A saved source changed during endpoint postprocessing")
    np.savez_compressed(output/"profiles.npz", x=common_x,
        **{key+"_41": value for key,value in fields41.items()},
        **{key+"_81": value for key,value in fields81.items()}, **payload, **corrections)
    protocol.static.write_json(output/"metrics.json", metrics)
    with (output/"metrics.csv").open("w", newline="", encoding="utf8") as stream:
        writer = csv.writer(stream); writer.writerow(("response", "state", "field", "absolute_max", "L2_41_grid", "signed_81_minus_41_at_max", "maximum_x"))
        for response in ("full", "correction", "evolution"):
            states = metrics[response] if response != "evolution" else {"final_minus_static": metrics[response]}
            for state, values in states.items():
                for field in FIELDS:
                    row = values[field]
                    writer.writerow((response, state, field, row["absolute_max"], row["L2_on_41_point_trapezoidal_grid"], row["signed_81_minus_41_at_max"], row["maximum_x"]))
    protocol.static.write_json(output/"manifest.json", {"schema": "bounded-spatial-recovery-endpoints-v1",
        "source_hashes": sources, "artifact_hashes": {path.relative_to(output).as_posix(): native.sha(path)
            for path in sorted(output.rglob("*")) if path.is_file()}})
    return metrics
