"""Scoped seven-field dead-load/recovery adapter; no native execution entry point.

Historical FEM generators, C3D10 quadrature, field interpolation and streaming
readers remain unchanged. This adapter supplies the new two-component load and
uses the complete polar orientation for spatial exponential coordinates.
"""
from __future__ import annotations

import math
import re
import time
from pathlib import Path

import numpy as np
from scipy.spatial.transform import Rotation

from scripts.analysis import resume_nlsp_nonlinear_dynamic_3d_fem as resume

pilot = resume.base
static = pilot.base
fem1 = static.fem1
io = pilot.io
FIELD_ORDER = ("u", "w", "v", "Phi", "psi", "theta", "c_eff")
LOCAL_SIGNS = np.array((1., -1., -1.))
DEFAULT_DYNAMIC = {"alpha": 0, "initial_T1_fraction": 1/8000,
    "maximum_T1_fraction": 1/4000, "minimum_initial_fraction": 1e-4,
    "maximum_increments": 2000, "output_frequency": 2,
    "release": "OP=NEW plus zero GRAV; STEP AMPLITUDE=STEP"}
RESOURCE_POLICY = {"threads": 1, "memory_limit_bytes": 4*1024**3,
    "job_timeouts_seconds": {"medium": 1800, "fine": 4200},
    "total_CCX_budget_seconds": 14400, "maximum_production_jobs": 4,
    "automatic_retry": False}


def load_contract(g_n, g_k, *, rho=1., length=1., width=.2, thickness=.1):
    """Fixed global acceleration, through the undeformed section centroid."""
    values = np.asarray((g_n, g_k, rho, length, width, thickness), float)
    if not np.all(np.isfinite(values)) or np.any(values[2:] <= 0):
        raise ValueError("Finite loads and positive reference geometry/material required")
    acceleration = np.array((0., -float(g_n), -float(g_k)))
    magnitude = float(np.linalg.norm(acceleration))
    if magnitude <= 0:
        raise ValueError("Nonzero dead preload required")
    area = float(width*thickness)
    return {"local_acceleration": [0., float(g_n), float(g_k)],
        "global_acceleration": acceleration.tolist(), "magnitude": magnitude,
        "global_direction": (acceleration/magnitude).tolist(),
        "line_load_local": [0., float(rho*area*g_n), float(rho*area*g_k)],
        "total_force_global": (rho*area*length*acceleration).tolist(),
        "reference_area": area, "reference_volume": area*length,
        "follower_load": False, "distributed_torque": 0.,
        "local_basis_global": np.diag(LOCAL_SIGNS).tolist()}


def _numeric_cards(text):
    """Audit every nonempty numerical field in native width-limited cards."""
    widths = []
    card = ""
    for line in text.splitlines():
        if line.startswith("**") or not line.strip():
            continue
        if line.startswith("*"):
            card = line.upper().split(",")[0]
            continue
        if card not in ("*STATIC", "*DYNAMIC", "*DLOAD", "*CONTROLS", "*ELASTIC", "*DENSITY"):
            continue
        tokens = line.split(",")[2:] if card == "*DLOAD" else line.split(",")
        for token in tokens:
            token = token.strip()
            if token:
                if len(token) > 20 or not math.isfinite(float(token)):
                    raise ValueError("Native numerical field exceeds finite 20-character contract")
                widths.append(len(token))
    return max(widths, default=0)


def input_contract(text, science, *, nonlinear):
    """Check the new physical load without borrowing old planar preload gates."""
    width = _numeric_cards(text)
    safety = resume.output_safety(text)
    parts = text.split("*END STEP")
    if len(parts) != 3 or parts[2].strip():
        raise ValueError("Exactly one STATIC and one DYNAMIC step required")
    first, dynamic = parts[:2]
    lines0 = first.splitlines()
    for card, names in (("*ELASTIC", ("E", "nu")), ("*DENSITY", ("rho",))):
        position = next(i for i, line in enumerate(lines0) if line == card)
        values = tuple(map(float, lines0[position+1].split(",")))
        if values != tuple(science["material"][name] for name in names):
            raise ValueError("Frozen preload material/density changed")
    position = next(i for i, line in enumerate(lines0) if line == "*BOUNDARY")
    boundaries = []
    for line in lines0[position+1:]:
        if line.startswith("*"):
            break
        boundaries.append(line.replace(" ", ""))
    if boundaries != ["LEFT_FIXED,1,3,0", "RIGHT_FIXED,1,3,0"]:
        raise ValueError("Only the two original full translational end-face clamps are allowed")
    if "*BOUNDARY" in dynamic:
        raise ValueError("Dynamic step may not introduce additional constraints")
    if "*STATIC" not in first or "*DYNAMIC, ALPHA=0" not in dynamic:
        raise ValueError("Same-job STATIC/direct DYNAMIC protocol required")
    for block in (first, dynamic):
        step = next(line for line in block.splitlines() if line.startswith("*STEP"))
        actual_nl = "NLGEOM" in step and "NLGEOM=NO" not in step
        if actual_nl != bool(nonlinear):
            raise ValueError("Linear/nonlinear routing changed")
    for forbidden in ("*DAMPING", "*CONTACT", "*SPRING", "*MPC", "*FREQUENCY", "*MODAL", "*CLOAD"):
        if forbidden in text.upper():
            raise ValueError("Unauthorized additional physics: " + forbidden)
    loads = []
    for block in (first, dynamic):
        rows = [line.split(",") for line in block.splitlines() if line.startswith("SOLID,GRAV,")]
        if len(rows) != 1:
            raise ValueError("One resultant GRAV row per step required")
        magnitude = float(rows[0][2]); direction = np.array(list(map(float, rows[0][3:])))
        if direction.shape != (3,) or np.linalg.norm(direction) == 0:
            raise ValueError("Nonzero three-component GRAV direction required even at zero magnitude")
        loads.append(magnitude*direction/np.linalg.norm(direction))
    expected = np.asarray(science["load"]["global_acceleration"])
    if not np.allclose(loads[0], expected, rtol=1e-12, atol=1e-18):
        raise ValueError("Serialized load differs from frozen two-component contract")
    if np.any(loads[1] != 0) or "*DLOAD, OP=NEW" not in dynamic or "AMPLITUDE=STEP" not in dynamic:
        raise ValueError("All previous body loads must be removed instantaneously")
    initial = "*INITIAL CONDITIONS, TYPE=VELOCITY\nALL_NODES,1,0.\nALL_NODES,2,0.\nALL_NODES,3,0."
    if initial not in text:
        raise ValueError("Zero physical initial velocities required")
    expected_time = pilot.dynamic_settings(science)
    lines = dynamic.splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith("*DYNAMIC"))
    actual_time = tuple(map(float, lines[index+1].split(",")))
    target_time = tuple(expected_time[name] for name in
        ("initial_increment", "duration", "minimum_increment", "maximum_increment"))
    if not np.allclose(actual_time, target_time, rtol=1e-12, atol=0):
        raise ValueError("Frozen dynamic time policy changed")
    frequency = science["dynamic"]["output_frequency"]
    if any(f"FREQUENCY={frequency}" not in line for line in lines
           if line.startswith(("*NODE FILE", "*EL FILE", "*NODE PRINT", "*EL PRINT"))):
        raise ValueError("Inconsistent dynamic output cadence")
    return {**safety, "status": "PASS", "numeric_width_max": width,
        "load": science["load"], "dynamic_settings": expected_time,
        "zero_velocity_initial_condition": True, "zero_bodyloads_after_release": True,
        "same_job_state_transfer": True, "historical_preload_equality_required": False,
        "contraction_qualification": "effective proxy, not identical M-H coordinate"}


def write_spatial_input(path, source_mesh_directory, mesh, audit, *, material,
                        g_n, g_k, omega1, static_settings, nonlinear,
                        dynamic=None, horizon_T1=.25):
    """Reuse the tested generator; adapt only requested load and output cadence."""
    if material != {"E": 1., "rho": 1., "nu": .3, "kappa": 5/6}:
        raise ValueError("Frozen FEM material changed")
    policy = dict(DEFAULT_DYNAMIC if dynamic is None else dynamic)
    if policy != DEFAULT_DYNAMIC or horizon_T1 != .25 or float(omega1) != .6054167303477958:
        raise ValueError("Only the predeclared quarter-period refined-time policy is allowed")
    load = load_contract(g_n, g_k, rho=material["rho"])
    science = {"material": dict(material), "g": load["magnitude"], "load": load,
        "omega1": float(omega1), "horizon_T1": float(horizon_T1), "dynamic": policy}
    old = {"science_config": {"static_settings": dict(static_settings)}}
    pilot.write_input(path, science, old, Path(source_mesh_directory), mesh, audit, nonlinear)
    path = Path(path)
    text = path.read_text(encoding="utf8")
    first, dynamic_text = text.split("*END STEP\n", 1)
    fmt = fem1.single.ccx_float
    old_row = f"SOLID,GRAV,{fmt(load['magnitude'])},0,-1,0"
    new_row = "SOLID,GRAV," + ",".join(fmt(value) for value in
        [load["magnitude"], *load["global_direction"]])
    if first.count(old_row) != 1:
        raise ValueError("Historical preload generator no longer matches audited adapter")
    first = first.replace(old_row, new_row)
    frequency = policy["output_frequency"]
    dynamic_text = re.sub(r"FREQUENCY=1\b", f"FREQUENCY={frequency}", dynamic_text)
    marker = f"*NODE FILE, GLOBAL=YES, FREQUENCY={frequency}"
    dynamic_text = dynamic_text.replace(marker,
        f"*EL FILE, GLOBAL=YES, FREQUENCY={frequency}\nS,E,ENER\n" + marker, 1)
    text = first + "*END STEP\n" + dynamic_text
    path.write_text(text, encoding="utf8")
    return science, input_contract(text, science, nonlinear=nonlinear)


def consistent_bodyloads(mesh, rho, acceleration_global, *, quadrature=None):
    """Integrate the vector dead load with existing reference C3D10 quadrature."""
    acceleration = np.asarray(acceleration_global, float)
    if acceleration.shape != (3,) or not np.isfinite(acceleration).all() or rho <= 0:
        raise ValueError("Finite global acceleration and positive density required")
    quad = fem1.quadrature_arrays(mesh, rho) if quadrature is None else quadrature
    ids, xyz, _, _ = fem1.mesh_arrays(mesh)
    if not np.array_equal(ids, quad["node_ids"]):
        raise ValueError("Bodyload quadrature node order mismatch")
    weights = np.einsum("eq,qi->ei", quad["weights"], quad["N"])
    loads = np.zeros((len(ids), 3))
    np.add.at(loads, np.asarray(quad["conn"]).ravel(),
              (weights[:, :, None]*acceleration).reshape(-1, 3))
    return ids, xyz, loads, float(np.sum(quad["weights"])/rho)


def support_equilibrium(mesh, audit, U, support_RF, *, rho, acceleration_global,
                        nonlinear, equilibrium_relative=1e-5, quadrature=None):
    """Recover supports independently of balance; RF includes bodyload at supports."""
    if equilibrium_relative != 1e-5:
        raise ValueError("Historical equilibrium_relative gate must remain unchanged")
    ids, xyz, loads, volume = consistent_bodyloads(mesh, rho, acceleration_global,
                                                 quadrature=quadrature)
    U = np.asarray(U, float)
    if U.shape != xyz.shape or not np.isfinite(U).all():
        raise ValueError("Complete finite preload nodal displacement required")
    lookup = {int(node): i for i, node in enumerate(ids)}
    force = loads.sum(axis=0); scale = float(np.linalg.norm(force))
    if not scale:
        raise ValueError("Nonzero applied force required for equilibrium normalization")
    current = xyz+U if nonlinear else xyz
    recovered, positions, results = [], [], {}
    for name, key in (("LEFT_FIXED", "fixed_left_ids"), ("RIGHT_FIXED", "fixed_right_ids")):
        rows = np.array([lookup[int(node)] for node in audit[key]])
        RF = np.asarray(support_RF[name], float)
        if RF.shape != (len(rows), 3) or not np.isfinite(RF).all():
            raise ValueError("Complete finite support RF required")
        reaction = RF-loads[rows]
        recovered.append(reaction); positions.append(current[rows])
        results[name] = {"force": reaction.sum(axis=0), "raw_RF": RF.sum(axis=0),
            "consistent_bodyload": loads[rows].sum(axis=0), "moment_about_face_centroid":
            np.cross(current[rows]-current[rows].mean(axis=0), reaction).sum(axis=0)}
        if np.max(abs(U[rows])) > 1e-12:
            raise ValueError("Fixed face has nonzero preload displacement")
    reactions = np.vstack(recovered)
    force_residual = reactions.sum(axis=0)+force
    moment_residual = np.cross(np.vstack(positions), reactions).sum(axis=0) + np.cross(current, loads).sum(axis=0)
    length = float(np.ptp(xyz[:, 0]))
    force_relative = float(np.linalg.norm(force_residual)/scale)
    moment_relative = float(np.linalg.norm(moment_residual)/(scale*length))
    return {"status": "PASS" if max(force_relative, moment_relative) <= equilibrium_relative else "FAIL",
        "total_applied_force_global": force, "reference_volume": volume,
        "supports": results, "force_imbalance_global": force_residual,
        "moment_imbalance_global": moment_residual,
        "force_imbalance_relative": force_relative, "moment_imbalance_relative": moment_relative,
        "equilibrium_relative_gate": equilibrium_relative,
        "reaction_contract": "support RF minus independently integrated consistent reference bodyload"}


def canonical_section_profile(profile):
    """Use log(R) for every spatial rotation coordinate, retain old projected theta."""
    fields = np.asarray(profile["fields"], float).copy()
    rows = profile["section_rows"]
    rotations = np.asarray([row["finite_rotation_matrix"] for row in rows])
    a = Rotation.from_matrix(rotations).as_rotvec()
    qrot = a*np.array((1., -1., 1.))
    clamped = bool(profile["clamped_face_endpoint_values_used"])
    target = slice(1, -1) if clamped else slice(None)
    fields[target, 3:6] = qrot
    result = dict(profile)
    result.update(fields=fields, historical_projected_theta=np.asarray(profile["fields"])[:, 5].copy(),
        raw_rotation_matrices=rotations, raw_rotation_coordinates=qrot,
        theta_policy="canonical exponential theta=log(R)_3; historical projected angle retained separately",
        rotation_contract="local a=(Phi,-psi,theta), B_global=diag(1,-1,-1)",
        contraction_status="DIAGNOSTIC_EFFECTIVE_THICKNESS_STRETCH_NOT_MH_DOF")
    return result


def recover_spatial_sections(mesh, nodal_U, *, section_count=41, quadrature=None):
    quad = fem1.quadrature_arrays(mesh, 1.) if quadrature is None else quadrature
    displacement = fem1.nlsp_evaluate_tet10_displacements(nodal_U, quad["conn"], quad["N"])
    original = static.fem2_recover_reference_samples(quad["xyz"], displacement,
        quad["weights"], 1., .1, .2, section_count)
    return canonical_section_profile(original)


def recover_saved_outputs(case, mesh, audit, science, *, nonlinear,
                          equilibrium_relative=1e-5, section_count=41):
    """Read a finished native job; never launch, replay or invent a missing frame.

All raw files stay intact if a parser or quality gate rejects the output. The
new STATIC state is checked on its own load; an old planar equilibrium is not
used as an equality reference for this distinct physical excitation.
"""
    started = time.perf_counter()
    case = Path(case)
    log_lines = (case/"motion.stdout.txt").read_text(encoding="utf8", errors="strict").splitlines()
    if (not any("JOB FINISHED" in line.upper() for line in log_lines)
        or any("*ERROR" in line.upper() or "*WARNING" in line.upper() for line in log_lines)):
        raise ValueError("Native completion missing or unexplained warning/error")
    actual_nl = any("Nonlinear geometric effects are taken into account" in line for line in log_lines)
    if actual_nl != bool(nonlinear):
        raise ValueError("Actual native nonlinear routing mismatch")
    sta = io.read_transient_sta(case/"motion.sta")
    static.write_json(case/"increments.json", sta)
    increments = sta["accepted_increments"]
    preloads = [row for row in increments if row["step"] == 1]
    dynamics = [row for row in increments if row["step"] == 2]
    if not preloads or not dynamics or abs(preloads[-1]["step_time"]-1.) > 1e-6:
        raise ValueError("Full STATIC preload or actual DYNAMIC increments missing")
    end = pilot.dynamic_settings(science)["duration"]
    last = dynamics[-1]
    bound = last["time_rounding_bounds"]["step_time"] + 8*np.finfo(float).eps*max(1., end)
    if abs(last["step_time"]-end) > bound:
        raise ValueError("Accepted dynamic prefix does not reach frozen target horizon")
    offset = preloads[-1]["total_time"]
    ids, _, _, _ = fem1.mesh_arrays(mesh)
    sets = {"ALL_NODES": ids, "LEFT_FIXED": audit["fixed_left_ids"], "RIGHT_FIXED": audit["fixed_right_ids"]}
    datdir = case/"dat_fields"; datdir.mkdir(exist_ok=True)
    blocks = {}
    for block in io.iter_transient_dat(case/"motion.dat", sets, static_end_time=offset, increments=increments):
        key = (block["step"], block["increment"], block["set"], block["name"])
        destination = datdir/("_".join(map(str, key))+".npz")
        np.savez_compressed(destination, values=block["values"])
        blocks[key] = destination
    quad = fem1.quadrature_arrays(mesh, science["material"]["rho"])
    fixed = np.r_[audit["fixed_left_ids"], audit["fixed_right_ids"]]
    frames = case/"frames"; frames.mkdir(exist_ok=True)
    rows = []; static_final = None; maximum_rounding = 0.; minimum_det = math.inf; maximum_strain = 0.
    x = np.linspace(0., 1., section_count)
    for frame in io.iter_transient_frd(case/"motion.frd", ids, static_end_time=offset,
                                     fixed_node_ids=fixed, increments=increments):
        key = (frame["step"], frame["increment"], "ALL_NODES", "DISP")
        if key not in blocks:
            raise ValueError("Native FRD frame has no complete matching DAT displacement")
        with np.load(blocks[key], allow_pickle=False) as saved:
            U = saved["values"].copy()
        difference = float(np.max(abs(U-frame["fields"]["DISP"])))
        maximum_rounding = max(maximum_rounding, difference)
        if (difference > 1e-8 or frame["fixed_displacement_max"] > 1e-12
            or (frame.get("fixed_velocity_max") or 0.) > 1e-12):
            raise ValueError("Unchanged DAT/FRD or fixed-face gate failed")
        destination = frames/f"step{frame['step']}_inc{frame['increment']:05d}.npz"
        np.savez_compressed(destination, node_ids=ids, U=U,
            **{key: value for key, value in frame["fields"].items() if key != "DISP"})
        for name in ("STRESS", "TOSTRAIN"):
            if name not in frame["fields"] or not np.isfinite(frame["fields"][name]).all():
                raise ValueError("Complete finite native stress/strain field required")
        if frame["step"] == 1:
            static_final = frame, U
            continue
        strain = static.fem2_fe_strain_diagnostics(mesh, U)
        if not strain["finite_values"] or strain["minimum_det_deformation_gradient"] <= 0:
            raise ValueError("Nonfinite or inverted saved dynamic deformation")
        minimum_det = min(minimum_det, strain["minimum_det_deformation_gradient"])
        maximum_strain = max(maximum_strain, strain["green_lagrange_max_abs"])
        profile = recover_spatial_sections(mesh, U, section_count=section_count, quadrature=quad)
        velocity = fem1.nlsp_evaluate_tet10_displacements(frame["fields"]["VELO"], quad["conn"], quad["N"])
        velocity_profile = static.fem2_recover_reference_samples(quad["xyz"], velocity,
            quad["weights"], 1., .1, .2, section_count)
        actual = frame["dynamic_time"]; canonical_time = actual
        if frame["increment"] == last["increment"]:
            rounding = frame.get("total_time_rounding_bound", 1e-8) + 8*np.finfo(float).eps*max(1., end)
            if abs(actual-end) > rounding:
                raise ValueError("Final native output misses target-time rounding interval")
            canonical_time = end
        rows.append({"time": canonical_time, "printed_dynamic_time": actual,
            "increment": frame["increment"], "fields": static.fem2_static_sample(profile, x),
            "raw_fields": profile["fields"],
            "raw_rotation_matrices": profile["raw_rotation_matrices"],
            "historical_projected_theta": profile["historical_projected_theta"],
            "translation_velocities": static.fem2_static_sample(velocity_profile, x)[:, :3],
            "kinetic_energy": .5*float(np.sum(quad["weights"]*np.sum(velocity*velocity, axis=2))),
            "frame": str(destination.relative_to(case)),
            "native_metadata": {key: value for key, value in frame.items() if key != "fields"}})
    if static_final is None or not rows or rows[-1]["increment"] != last["increment"] or rows[-1]["time"] != end:
        raise ValueError("Actual STATIC or final DYNAMIC output missing")
    frame, U = static_final
    preload_strain = static.fem2_fe_strain_diagnostics(mesh, U)
    static.write_json(case/"preload_strain_diagnostics.json", preload_strain)
    if not preload_strain["finite_values"] or preload_strain["minimum_det_deformation_gradient"] <= 0:
        raise ValueError("Nonfinite or inverted saved STATIC deformation")
    profile0 = recover_spatial_sections(mesh, U, section_count=section_count, quadrature=quad)
    RF = {}
    for name in ("LEFT_FIXED", "RIGHT_FIXED"):
        with np.load(blocks[(1, frame["increment"], name, "FORC")], allow_pickle=False) as saved:
            RF[name] = saved["values"].copy()
    equilibrium = support_equilibrium(mesh, audit, U, RF, rho=science["material"]["rho"],
        acceleration_global=science["load"]["global_acceleration"], nonlinear=nonlinear,
        equilibrium_relative=equilibrium_relative, quadrature=quad)
    static.write_json(case/"preload_equilibrium.json", equilibrium)
    if equilibrium["status"] != "PASS":
        raise ValueError("New two-component preload fails unchanged equilibrium gate")
    fields0 = static.fem2_static_sample(profile0, x)
    for component in (1, 2):
        if fields0[section_count//2, component] != 0 and (
            (rows[0]["fields"][section_count//2, component]-fields0[section_count//2, component])
            *fields0[section_count//2, component] >= 0):
            raise ValueError("Early free transverse movement is not restoring")
    np.savez_compressed(case/"initial_sections.npz", x=x, fields=fields0,
        raw_x=profile0["x"], raw_fields=profile0["fields"],
        raw_rotation_x=np.asarray(profile0["x"])[1:-1],
        raw_rotation_matrices=profile0["raw_rotation_matrices"],
        historical_projected_theta=profile0["historical_projected_theta"],
        source_static_end_time=offset, not_a_native_dynamic_zero_frame=np.array(True))
    np.savez_compressed(case/"section_history.npz", x=x,
        raw_x=np.asarray(profile0["x"]), raw_rotation_x=np.asarray(profile0["x"])[1:-1],
        time=np.array([row["time"] for row in rows]),
        printed_dynamic_time=np.array([row["printed_dynamic_time"] for row in rows]),
        fields=np.stack([row["fields"] for row in rows]),
        raw_fields=np.stack([row["raw_fields"] for row in rows]),
        raw_rotation_matrices=np.stack([row["raw_rotation_matrices"] for row in rows]),
        historical_projected_theta=np.stack([row["historical_projected_theta"] for row in rows]),
        translation_velocities=np.stack([row["translation_velocities"] for row in rows]),
        increments=np.array([row["increment"] for row in rows]))
    static.write_json(case/"frame_metadata.json", [{key: value for key, value in row.items()
        if key not in ("fields", "raw_fields", "translation_velocities", "raw_rotation_matrices", "historical_projected_theta")}
        for row in rows])
    np.savez_compressed(case/"independent_kinetic_energy.npz",
        time=np.array([row["time"] for row in rows]), kinetic_energy=np.array([row["kinetic_energy"] for row in rows]))
    energy = io.parse_transient_dat_energies(case/"motion.dat", static_end_time=offset,
                                          increments=increments, element_set="SOLID")
    stdout = io.parse_transient_stdout_energies(case/"motion.stdout.txt", static_end_time=offset, increments=increments)
    static.write_json(case/"energy.json", energy); static.write_json(case/"stdout_energy.json", stdout)
    dynamic_energy = [row for row in stdout["records"] if row.get("step") == 2]
    if not dynamic_energy or any("external_work" not in row or "damping_work" not in row for row in dynamic_energy):
        raise ValueError("Actual free-motion external/damping work evidence missing")
    external = max(abs(row["external_work"]) for row in dynamic_energy)
    damping = max(abs(row["damping_work"]) for row in dynamic_energy)
    if external != 0 or damping != 0:
        raise ValueError("Nonzero external or damping work after release")
    initial_energy = [row for row in energy["records"] if row["step"] == 1 and abs(row["total_time"]-offset) < 1e-7]
    record = {"status": "PASS", "dynamic_time_end": end, "dynamic_time_start": rows[0]["time"],
        "dynamic_output_frames": len(rows), "static_increments": len(preloads),
        "dynamic_increments": len(dynamics), "accepted_increments": len(increments),
        "cutbacks": sta["reported_cutbacks"], "final_native_frame_reached": True,
        "preload_equilibrium": equilibrium, "preload_strain_diagnostics": preload_strain,
        "max_DAT_FRD_displacement_difference": maximum_rounding,
        "maximum_native_external_work_after_release": external, "maximum_native_damping_work_after_release": damping,
        "native_initial_internal_energy": initial_energy[-1]["internal_energy"] if initial_energy else None,
        "energy_status": "PARTIAL", "independent_internal_energy": "NOT_RUN",
        "all_saved_frames_strain_diagnostics": {"finite_values": True,
            "minimum_det_deformation_gradient": minimum_det, "max_abs_green_lagrange_strain": maximum_strain},
        "initial_midspan": fields0[section_count//2].tolist(), "final_midspan": rows[-1]["fields"][section_count//2].tolist(),
        "static_to_dynamic_transfer": "consecutive same-job steps; zero velocity IC; STATIC t0 anchor distinguished from native positive-time DYNAMIC frames",
        "historical_preload_equality_required": False, "canonical_rotation_coordinates": True,
        "contraction_qualification": "effective proxy, not identical M-H coordinate", "recovery_seconds": time.perf_counter()-started}
    resume.qualify_native_energy(case, record)
    static.write_json(case/"recovery.json", record)
    return record
