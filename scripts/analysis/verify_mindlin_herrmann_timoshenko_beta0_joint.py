"""Reduced rigid-joint audit: beta=0, three splits of one fixed G20 rod only.

Different contract from direct single-rod/source CLIs: 16 boundary/joint
rows, transparency, interface work, modes and reflection. Reuses their arm
physics, immutable reference bundle and artifact writers; no angle runner.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np

from scripts.analysis import verify_mindlin_herrmann_timoshenko_single_rod as single
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint

OUTPUT = ROOT/"results/mindlin_herrmann_timoshenko_beta0_joint"
SPLITS = (.5, .35, .65)
# Defined before first root calculation; unchanged single-rod root policy.
POLICY = {"frequency_relative_tol": 2e-8, "MAC_loss_tol": 1e-8,
    "component_L2_relative_tol": 1e-6, "interface_reference_tol": 1e-6,
    "interface_scaled_residual_tol": 1e-9, "clamp_scaled_residual_tol": 1e-9,
    "virtual_work_relative_tol": 2e-14, "segmented_transfer_tol": 1e-9,
    "cluster_relative_gap": 1e-6, "quadrature_order": 200,
    "reported_combined": 12, "guard_roots": 1, "beta_deg": 0.,
    "splits": list(SPLITS), "defined_before_computation": True}
SOURCE_CLOSURE = {
    "qualification": "Contraction compatibility is adopted from published common-DOF reduced frame assembly, not derived from finite 3D joint elasticity.",
    "rucka": {"citation_key": "rucka_2010_l_joint_guided_waves",
        "directly_printed": "PDF4/p1763 (10)-(15): q=(u,psi,v,phi); PDF5/p1764 (26)-(28): T_i, Kbar=T^T K T, fbar=T^T f, standard aggregation; section5 L-joint application",
        "transcription": "T_i=[[cos(alpha),0,sin(alpha),0],[0,1,0,0],[-sin(alpha),0,cos(alpha),0],[0,0,0,1]]; q_local=T_i q_global from (27)",
        "inference": "Shared psi/phi nodal DOFs imply common c/theta; dual nodal efforts sum. Scalar c1=c2 is not a separately printed joint equation."},
    "jang": {"citation_key": "jang_2014_timoshenko_composite_patch_guided_waves",
        "directly_printed": "PDF2/p249 (1),(5); PDF3/p250 (6),(12): lateral psi variation and conjugate R; PDF5/p252 (28)-(30),(36): nodal displacement/contraction DOFs and boundary work",
        "inference": "Bare base reduction retains four fields and pair (psi_b,R_b); no attribution of Rucka frame aggregation to Jang."}}


def check_reference(*, allow_compute=False):
    config, model, length, checked = single.check_inputs()
    pointer = single.OUTPUT/"current.json"
    if not pointer.exists():
        if not allow_compute:
            raise ValueError("Direct reference absent; run single-rod --compute or joint --compute first")
        single.main(["--compute"])
    directory = single.OUTPUT/json.loads(pointer.read_text(encoding="utf-8"))["fingerprint"]
    manifest = json.loads((directory/"manifest.json").read_text(encoding="utf-8"))
    _, current_identity = single.identity(config, checked)
    # A new commit alone does not invalidate immutable physical/code evidence.
    historical, current = dict(manifest["identity"]), dict(current_identity)
    historical.pop("head")
    current.pop("head")
    if historical != current or any(not (directory/p).exists() or sha(directory/p) != h
        for p, h in manifest["artifact_hashes"].items()):
        raise ValueError("Direct reference code/input/version or artifacts changed; run unchanged single-rod CLI for a matching bundle")
    result = json.loads((directory/"result.json").read_text(encoding="utf-8"))
    if result["finite_spectrum_status"] != "PASS" or result["length"] != length or result["coefficients"] != model.coefficients:
        raise ValueError("Reference is not the same verified single rod")
    return config, model, length, checked, directory, manifest, result


def virtual_work_check():
    # General frames test the duality, not nonzero-beta spectra.
    frames = (joint.Frame((.6, .8), (.8, -.6)), joint.Frame((-.8, .6), (.6, .8)))
    rng = np.random.default_rng(1783)
    errors = []
    for ends in (("left", "right"), ("right", "right")):
        for _ in range(24):
            states, variation = rng.normal(size=(2, 8)), rng.normal(size=4)
            local, global_value = joint.virtual_work(states, variation, frames, ends)
            errors.append(abs(local-global_value)/max(1., abs(local), abs(global_value)))
    error = max(errors)
    if error > POLICY["virtual_work_relative_tol"] or np.linalg.matrix_rank(joint.joint_matrix()) != 8:
        raise ArithmeticError("Virtual-work gate failed; spectrum must not be calculated")
    return {"status": "PASS", "trials": len(errors), "max_scaled_error": error,
            "joint_rank": 8, "endpoint_signs": {"left": -1, "right": 1}}


def solve_split(model, length, split, config, reference):
    policy = config["policy"]
    record, profiles = {"split": split, "beta_deg": 0.}, []
    for block in ("mh", "timoshenko"):
        direct = reference[block]
        lower, upper = direct["search_attempts"][0]["range_omega"]
        count = mh.finite_count_upper_bound(model, length, upper, block, policy["young_eta"])
        start = mh.finite_count_upper_bound(model, length, lower, block, policy["young_eta"])
        attempts = []
        for attempt in range(2):
            roots, search = joint.roots(model, length, split, block, lower, upper,
                {**policy, "scan_intervals": policy["scan_intervals"]*2**attempt})
            attempts.append(search)
            if len(roots) == count["upper_count"] and not search["failed_intervals"]:
                break
        if len(roots) != count["upper_count"] or start["upper_count"] != 0 or search["failed_intervals"]:
            raise ArithmeticError(f"Count/interval gate failed: split={split}, {block}; attempts={attempts}")
        if len(roots) != len(direct["roots"]):
            raise ArithmeticError("Independent family inventories differ")
        for n, (root, ref) in enumerate(zip(roots, direct["roots"]), 1):
            root.update({"family_index": n, "family": ref["family"], "direct_frequency_hz": ref["frequency_hz"],
                "absolute_difference": abs(root["frequency_hz"]-ref["frequency_hz"]),
                "relative_difference": abs(root["omega"]/ref["omega"]-1)})
            coefficients = joint.mode_coefficients(model, length, split, root["omega"], block, POLICY["quadrature_order"])
            refmode = mh.finite_mode(model, length, ref["omega"], block, POLICY["quadrature_order"])
            diagnostic = joint.mode_diagnostics(model, length, split, root["omega"], block, coefficients,
                ref["omega"], refmode["coefficients"], POLICY["quadrature_order"])
            propagated = joint.segmented_transfer(model, length, split, root["omega"], block,
                policy["qr_step_exponent_cap"], policy["qr_max_steps"])
            root.update({"coefficients": coefficients.tolist(), "mode_comparison": diagnostic,
                         "segmented_transfer": propagated})
            if (not root["converged"] or root["relative_difference"] > POLICY["frequency_relative_tol"] or
                root["singular_ratio"] > policy["boundary_scaled_residual_tol"] or
                root["nonzero_singular_condition"] > policy["nonzero_singular_condition_max"] or
                1-diagnostic["mass_MAC"] > POLICY["MAC_loss_tol"] or
                max(diagnostic["component_L2_relative"].values()) > POLICY["component_L2_relative_tol"] or
                max(diagnostic["joint_residual_scaled"].values()) > POLICY["interface_scaled_residual_tol"] or
                diagnostic["clamp_scaled_residual"] > POLICY["clamp_scaled_residual_tol"] or
                diagnostic["interface_reference_scaled_error"] > POLICY["interface_reference_tol"] or
                propagated["direct_projected_boundary_difference"] > POLICY["segmented_transfer_tol"]):
                raise ArithmeticError(f"Joint mode gate failed: split={split}, block={block}, n={n}: {root}")
            points = np.unique(np.concatenate((np.linspace(0., length, 201), [split*length])))
            states = joint.full_state(joint.mode_state(model, length, split, root["omega"], block,
                coefficients, points)*diagnostic["sign_alignment"], block)
            direct_states = joint.full_state(mh.finite_state_basis(model, length, ref["omega"], points, block)@refmode["coefficients"], block)
            profiles.extend({"split": split, "family": root["family"], "family_index": n, "x": float(x),
                **{name: float(v) for name, v in zip(joint.STATE_ORDER, state)},
                **{"direct_"+name: float(v) for name, v in zip(joint.STATE_ORDER, reference_state)}}
                for x, state, reference_state in zip(points, states, direct_states))
        record[block] = {"roots": roots, "search_attempts": attempts, "completeness": count,
            "below_search_count_bound": start, "certificate_qualification": "Exact beta0 transparency of energy form/domain proved before search; saturated direct min-max bound applies", "status": "PASS"}
    combined = sorted([{**r, "block": block} for block in ("mh", "timoshenko") for r in record[block]["roots"]], key=lambda r: r["omega"])
    if any((b["omega"]-a["omega"])/b["omega"] < POLICY["cluster_relative_gap"] for a, b in zip(combined[:-1], combined[1:])):
        raise ArithmeticError("Cluster detected: individual-vector comparison not authorized; subspace audit needed")
    prefix = combined[:POLICY["reported_combined"]+POLICY["guard_roots"]]
    if prefix[-1]["omega"] >= min(record[b]["roots"][-1]["omega"] for b in ("mh", "timoshenko")):
        raise ArithmeticError("Combined guard exceeds a family inventory")
    record["combined"] = [{"sorted_position": n, "guard": n > POLICY["reported_combined"],
        **{k: r[k] for k in ("block", "family", "family_index", "frequency_hz", "direct_frequency_hz", "absolute_difference", "relative_difference")}}
        for n, r in enumerate(prefix, 1)]
    return record, profiles


def arm_swap(model, length, a, b):
    rows = []
    p = model.coefficients
    for block in ("mh", "timoshenko"):
        nodes, weights = joint.quadrature(length, .35, POLICY["quadrature_order"])
        for left, right in zip(a[block]["roots"], b[block]["roots"]):
            first = joint.mode_state(model, length, .35, left["omega"], block, np.array(left["coefficients"]), nodes)
            reflected = joint.mode_state(model, length, .65, right["omega"], block, np.array(right["coefficients"]), length-nodes)*[-1., 1., 1., -1.]
            mass = np.array([p["m"], p["j"] if block == "mh" else p["r"]])
            overlap = float(weights@((first[:, :2]*reflected[:, :2])@mass))
            overlap /= math.sqrt(float(weights@(first[:, :2]**2@mass))*float(weights@(reflected[:, :2]**2@mass)))
            reflected *= 1. if overlap >= 0 else -1.
            errors = [math.sqrt(float(weights@(first[:, i]-reflected[:, i])**2)/float(weights@first[:, i]**2)) for i in range(2)]
            row = {"block": block, "family_index": left["family_index"],
                "frequency_relative_difference": abs(left["omega"]/right["omega"]-1),
                "mass_MAC": min(1., overlap**2), "component_L2_relative": errors}
            if row["frequency_relative_difference"] > POLICY["frequency_relative_tol"] or 1-row["mass_MAC"] > POLICY["MAC_loss_tol"] or max(errors) > POLICY["component_L2_relative_tol"]:
                raise ArithmeticError(f"Arm-swap/reflection gate failed: {row}")
            rows.append(row)
    return {"status": "PASS", "mapping": "X->L-X, (u,c,w,theta,N,R,Q,M)->(-u,c,-w,theta,N,-R,Q,-M)", "rows": rows}


def compute(config, model, length, reference):
    work = virtual_work_check()  # HARD gate before first determinant evaluation
    # Independent permutations of equations/columns; exact zero mixed blocks.
    row_order = (0, 1, 4, 5, 8, 10, 12, 14, 2, 3, 6, 7, 9, 11, 13, 15)
    column_order = (*range(4), *range(8, 12), *range(4, 8), *range(12, 16))
    matrix_checks = []
    for split in SPLITS:
        matrix = joint.boundary_matrix(model, length, split, 3.)
        permuted = matrix[np.ix_(row_order, column_order)]
        if np.any(permuted[:8, 8:]) or np.any(permuted[8:, :8]):
            raise ArithmeticError("beta0 assembly creates artificial mixed blocks")
        matrix_checks.append({"split": split, "omega": 3., "raw_matrix": matrix.tolist(),
                              "mixed_entries": 0})
    records, profiles = [], []
    for split in SPLITS:
        record, fields = solve_split(model, length, split, config, reference)
        records.append(record)
        profiles.extend(fields)
    swapped = arm_swap(model, length, records[1], records[2])
    return {"statuses": {name: "PASS" for name in ("MHTIM_JOINT_SOURCE_CLOSURE", "MHTIM_JOINT_VARIATIONAL_FORM",
        "MHTIM_BETA0_HOMOGENEOUS_SPECTRUM", "MHTIM_BETA0_SPLIT_INVARIANCE", "MHTIM_BETA0_MODE_SHAPES",
        "MHTIM_BETA0_ARM_SWAP", "MHTIM_BETA0_JOINT_GATE")},
        "source_closure": SOURCE_CLOSURE, "virtual_work": work, "length": length,
        "matrix_checks": {"joint_matrix": joint.joint_matrix().tolist(),
            "joint_row_order": joint.JOINT_ROWS, "coefficient_order": "arm1(MH4,Tim4),arm2(MH4,Tim4)",
            "block_row_permutation": row_order, "block_column_permutation": column_order,
            "samples": matrix_checks},
        "coefficients": model.coefficients, "limits": mh.limits(model), "splits": records,
        "arm_swap": swapped, "policy": POLICY,
        "coordinate_contract": {"x1": "outer left clamp 0 -> joint L1; X=x1",
            "x2": "outer right clamp 0 -> joint L2; X=L-x2", "t1": [1, 0], "n1": [0, -1],
            "t2": [-1, 0], "n2": [0, 1], "rotation_axis": "k=t cross n=-EZ; theta scalar",
            "state_order": joint.STATE_ORDER, "joint_outward_signs": [1, 1],
            "external_clamp_outward_signs": [-1, -1], "beta_deg": 0.}}, profiles


def identity(config, checked, directory, manifest):
    _, base = single.identity(config, checked)
    inputs = {"version": joint.VERSION, "single_arm_identity": base, "policy": POLICY,
        "source_closure": SOURCE_CLOSURE, "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in
            (Path(__file__), ROOT/"scripts/lib/mindlin_herrmann_timoshenko_joint.py")},
        "reference_bundle": directory.relative_to(ROOT).as_posix(),
        "reference_manifest_sha256": sha(directory/"manifest.json"),
        "reference_artifact_hashes": manifest["artifact_hashes"]}
    fingerprint = hashlib.sha256(json.dumps(inputs, sort_keys=True).encode()).hexdigest()[:16]
    return fingerprint, inputs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if not args.check_sources and not args.compute:
        parser.error("Select --check-sources or --compute")
    # Also precedes optional generation of a missing direct reference.
    virtual_work_check()
    config, model, length, checked, directory, manifest, reference = check_reference(allow_compute=args.compute)
    if args.check_sources:
        print("Source closure qualified; dual virtual work PASS; reference hashes verified; beta=0 only")
    if not args.compute:
        return 0
    fingerprint, inputs = identity(config, checked, directory, manifest)
    out = args.output_dir/fingerprint
    manifest_path = out/"manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous["identity"] != inputs or any(not (out/p).exists() or sha(out/p) != h for p, h in previous["artifact_hashes"].items()):
            raise ValueError("Changed/stale beta0 joint bundle")
        print("Verified beta0 joint bundle; zero root evaluations:", out)
        return 0
    out.mkdir(parents=True, exist_ok=True)
    try:
        result, profiles = compute(config, model, length, reference)
    except (ValueError, ArithmeticError) as exc:
        write_json(out/"failure.json", {"status": "FAIL", "reason": str(exc), "identity": inputs,
            "bounded_audit": "coordinates -> endpoint work -> transform -> c/R -> theta/M -> N/Q -> count -> basis; no coefficient/tolerance/sign fitting or nonzero-angle calculation"})
        raise
    write_json(out/"result.json", result)
    write_json(out/"parameters.json", {"single_rod_config": config, "policy": POLICY,
        "material_section": vars(model.section), "variant": model.variant})
    write_csv(out/"mode_profiles.csv", profiles)
    write_csv(out/"frequencies.csv", [{"split": split["split"], **r} for split in result["splits"] for r in split["combined"]])
    artifacts = {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"}
    audits = {p.relative_to(ROOT).as_posix(): sha(p) for p in (OUTPUT/"source_audit").glob("*.png")}
    write_json(manifest_path, {"fingerprint": fingerprint, "identity": inputs, "artifact_hashes": artifacts,
        "source_page_images": audits, "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True, encoding="utf-8"),
        "reference_reused": True, "reference_root_evaluations": 0,
        "spectrum_semantics": "sorted positions and separate families at beta0 only; no tracking",
        "nonzero_angle_validated": False})
    write_json(args.output_dir/"current.json", {"fingerprint": fingerprint, "directory": str(out.resolve())})
    print("MHTIM_BETA0_JOINT_GATE PASS:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
