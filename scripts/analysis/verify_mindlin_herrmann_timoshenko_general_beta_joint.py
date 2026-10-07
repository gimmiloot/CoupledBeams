"""General-frame structural gates, then fixed 5/45/90-degree spectral pilot.

One new contract: geometry/duality/symmetry and count-certified coupled
inventory. Reuses the common joint helper, frozen beta0 CLI/reference and
unchanged single-rod arm physics. No angle map, hierarchy or tracking.
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

from scripts.analysis import verify_mindlin_herrmann_timoshenko_beta0_joint as zero
from scripts.analysis import verify_mindlin_herrmann_timoshenko_single_rod as single
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint

OUTPUT = ROOT/"results/mindlin_herrmann_timoshenko_general_beta_joint"
FROZEN = zero.OUTPUT/"3059d70b1b50ea2e"
ANGLES = (5., 45., 90.)
SMALL_ANGLES = (1e-6, 1e-4, 1e-2)
POLICY = {"defined_before_computation": True, "root_xtol": 1e-11, "root_rtol": 1e-12,
    "scan_intervals": 400, "max_subdivision_depth": 20, "max_subdivisions": 2000,
    "omega_min_pi_c0_L": .001, "omega_max_pi_c0_L": 7.5,
    "pole_mh_ceiling_pi_c0_arm": 4.5, "pole_tim_ss_index": 9.9,
    "pole_exclusion_relative": 1e-7, "schur_symmetry_tol": 1e-10,
    "count_inertia_margin_min": 4e-14, "pole_condition_max": 1e10,
    "coordinate_duality_tol": 2e-14, "matrix_symmetry_map_tol": 5e-14,
    "beta0_frequency_relative_tol": 1e-10, "beta0_component_L2_tol": 5e-9,
    "symmetry_frequency_relative_tol": 5e-10, "symmetry_MAC_loss_tol": 1e-10,
    "symmetry_component_L2_tol": 5e-8, "cluster_relative_gap": 1e-6,
    "zero_limit_relative_allowances": [1e-9, 1e-7, 1e-3],
    "boundary_joint_scaled_tol": 1e-9, "equation_scaled_tol": 1e-9,
    "energy_relative_tol": 5e-8, "mass_gram_tol": 5e-7,
    "nonzero_singular_condition_max": 1e8, "quadrature_order": 200,
    "reported_positions": 12, "guard_roots": 1,
    "symmetry_gate_ceiling_fstar": .7,
    "failure_policy": "bounded count subdivision only; unresolved count/pole/multiplicity or hard gate stops; no tolerance/physics/sign fitting"}
PREFIX = "MHTIM_GENERAL_BETA_"


def check_inputs():
    config, model, length, checked, _, _, direct = zero.check_reference()
    manifest = json.loads((FROZEN/"manifest.json").read_text(encoding="utf-8"))
    if any(sha(FROZEN/p) != h for p, h in manifest["artifact_hashes"].items()):
        raise ValueError("Frozen beta0 artifacts changed")
    subject = "scripts/lib/mindlin_herrmann_timoshenko_joint.py"
    for name, digest in manifest["identity"]["files"].items():
        if name != subject and sha(ROOT/name) != digest:
            raise ValueError(f"Frozen beta0 orchestration changed: {name}")
    _, current = single.identity(config, checked)
    previous = dict(manifest["identity"]["single_arm_identity"])
    previous.pop("head")
    current.pop("head")
    if previous != current:
        raise ValueError("Arm/source physics or accepted inputs changed")
    frozen = json.loads((FROZEN/"result.json").read_text(encoding="utf-8"))
    if frozen["statuses"]["MHTIM_BETA0_JOINT_GATE"] != "PASS":
        raise ValueError("No accepted beta0 reference")
    return config, model, length, checked, direct, frozen


def structural_checks(model, length):
    rng, samples = np.random.default_rng(2751), []
    for beta in (0., *ANGLES):
        frames = joint.frames(beta)
        operator = joint.joint_matrix(frames)
        work_error = 0.
        for _ in range(32):
            state, variation = rng.normal(size=(2, 8)), rng.normal(size=4)
            a, b = joint.virtual_work(state, variation, frames)
            work_error = max(work_error, abs(a-b)/max(1., abs(a), abs(b)))
        orthogonal = max(float(np.max(np.abs(f.nodal_transform.T@f.nodal_transform-np.eye(4)))) for f in frames)
        rank = int(np.linalg.matrix_rank(operator))
        rank_gram = float(np.max(np.abs(operator@operator.T-2*np.eye(8))))
        if orthogonal > POLICY["coordinate_duality_tol"] or work_error > POLICY["coordinate_duality_tol"] or rank != 8 or rank_gram > POLICY["coordinate_duality_tol"]:
            raise ArithmeticError("Coordinate/duality/rank hard gate failed")
        samples.append({"beta_deg": beta, "frames": [{"t": f.t, "n": f.n, "local_to_global": f.nodal_transform.tolist(),
            "global_to_local": f.nodal_transform.T.tolist()} for f in frames],
            "joint_operator": operator.tolist(), "rank": rank, "rank_gram_error": rank_gram,
            "orthogonality_error": orthogonal, "virtual_work_max_error": work_error})
    for frame, t, n in zip(joint.frames(0.), ((1, 0), (-1, 0)), ((0, -1), (0, 1))):
        if not np.array_equal(frame.t, t) or not np.array_equal(frame.n, n):
            raise ArithmeticError("Exact zero geometry recovery failed")
    right = joint.frames(90.)
    right_error = max(float(np.max(np.abs(f.nodal_transform-np.round(f.nodal_transform)))) for f in right)
    if right_error > POLICY["coordinate_duality_tol"]:
        raise ArithmeticError("Right-angle axis permutation failed")
    frames = joint.frames(45.)
    operator = joint.joint_matrix(frames)
    swapped = joint.joint_matrix(frames[::-1])[:, (*range(8, 16), *range(8))]
    swap_error = float(np.max(np.abs(swapped-np.array([-1.]*4+[1.]*4)[:, None]*operator)))
    reflected = joint.reflected_frames(frames)
    reflection_rows = np.array([1., -1., 1., -1.]*2)
    reflection_columns = np.tile(joint.MIRROR_STATE, 2)
    mirror_error = float(np.max(np.abs(joint.joint_matrix(reflected)*reflection_columns-reflection_rows[:, None]*operator)))
    if max(swap_error, mirror_error) > POLICY["matrix_symmetry_map_tol"]:
        raise ArithmeticError("Algebraic swap/reflection gate failed")
    # Local state stays split; only translations mix at the nonzero-angle joint.
    local = mh.full_harmonic_state_matrix(model, 3.)
    if np.any(local[np.ix_(joint.BLOCK_INDICES["mh"], joint.BLOCK_INDICES["timoshenko"])]) or not operator[0, 10]:
        raise ArithmeticError("Local/global coupling structure wrong")
    return {"geometry": samples, "right_angle_max_error": right_error,
        "arm_swap_operator_error": swap_error, "reflection_operator_error": mirror_error,
        "local_mixed_state_entries": 0, "scalar_c_theta_rotated": False}


def beta0_regression(config, model, length, direct, frozen):
    # Frozen CLI uses block extraction of the SAME general assembly at beta0.
    fresh, _ = zero.compute(config, model, length, direct)
    rows, matrix_error = [], 0.
    for sample in frozen["matrix_checks"]["samples"]:
        matrix = joint.frame_boundary_matrix(model, (length*sample["split"], length*(1-sample["split"])), sample["omega"], beta_deg=0.)
        matrix_error = max(matrix_error, float(np.max(np.abs(matrix-np.array(sample["raw_matrix"])))))
    for a, b in zip(fresh["splits"], frozen["splits"]):
        for block in ("mh", "timoshenko"):
            if len(a[block]["roots"]) != len(b[block]["roots"]):
                raise ArithmeticError("Frozen beta0 inventory changed")
            x, weights = joint.quadrature(length, a["split"], POLICY["quadrature_order"])
            for new, old in zip(a[block]["roots"], b[block]["roots"]):
                v = joint.mode_state(model, length, a["split"], new["omega"], block, np.array(new["coefficients"]), x)
                ref = joint.mode_state(model, length, b["split"], old["omega"], block, np.array(old["coefficients"]), x)
                p = model.coefficients
                masses = np.array([p["m"], p["j"] if block == "mh" else p["r"]])
                overlap = float(weights@((v[:, :2]*ref[:, :2])@masses))
                v *= 1. if overlap >= 0 else -1.
                errors = [math.sqrt(float(weights@(v[:, i]-ref[:, i])**2)/float(weights@ref[:, i]**2)) for i in range(2)]
                diff = abs(new["omega"]/old["omega"]-1)
                if diff > POLICY["beta0_frequency_relative_tol"] or max(errors) > POLICY["beta0_component_L2_tol"]:
                    raise ArithmeticError("Frozen beta0 frequency/profile regression failed")
                rows.append({"split": a["split"], "block": block, "family_index": new["family_index"],
                    "frequency_relative_difference": diff, "component_L2_relative": errors,
                    "mass_MAC": min(1., overlap**2), "joint_residuals": new["mode_comparison"]["joint_residual_scaled"]})
    if matrix_error > 1e-14:
        raise ArithmeticError("Frozen physical beta0 matrix changed")
    return {"status": "PASS", "matrix_max_absolute_error": matrix_error, "rows": rows,
        "fresh_arm_swap": fresh["arm_swap"], "split_statuses": fresh["statuses"]}


def pole_catalog(config, model, length):
    """One equal arm, complete fixed--fixed pole inventories, unchanged solver."""
    arm, c0 = length/2, math.sqrt(model.section.E/model.section.rho)
    result = {}
    for block in ("mh", "timoshenko"):
        lower = POLICY["omega_min_pi_c0_L"]*math.pi*c0/length
        upper = POLICY["pole_mh_ceiling_pi_c0_arm"]*math.pi*c0/arm if block == "mh" else math.sqrt(
            mh.blocks(model)[1].temporal(POLICY["pole_tim_ss_index"]*math.pi/arm)[0]["omega_squared"])
        bound = mh.finite_count_upper_bound(model, arm, upper, block)
        roots, search = mh.finite_roots(model, arm, block, lower, upper, config["policy"])
        if len(roots) != bound["upper_count"] or search["failed_intervals"]:
            raise ArithmeticError("Fixed-arm pole catalog not certified")
        for root in roots:
            independent = single.independent_root(model, arm, block, root, config["policy"])
            if independent["relative_difference"] > POLICY["beta0_frequency_relative_tol"]:
                raise ArithmeticError("Independent fixed-arm pole check failed")
            root["independent"] = independent
        result[block] = {"roots": roots, "search": search, "count_bound": bound}
    return result


def count_function(model, length, catalog):
    poles = sorted(r["omega"] for b in catalog.values() for r in b["roots"])
    ceiling = min(b["search"]["range_omega"][1] for b in catalog.values())
    def raw(w, beta, frames):
        if w >= ceiling:
            raise ArithmeticError("Count exceeds certified pole coverage")
        _, diagnostic = joint.nodal_schur_matrix(model, (length/2, length/2), w, beta_deg=beta, arm_frames=frames)
        if diagnostic["scaled_skew_residual"] > POLICY["schur_symmetry_tol"]:
            raise ArithmeticError("Energy Schur symmetry failed")
        j0 = 2*sum(p < w for p in poles)
        return j0+diagnostic["negative_inertia"], {"J0": j0, **diagnostic}
    def count(w, beta, frames=None):
        exclusion = POLICY["pole_exclusion_relative"]*max(1., w)
        if min(abs(w-p) for p in poles) <= exclusion:
            a, b = w-2*exclusion, w+2*exclusion
            ca, da = raw(a, beta, frames)
            cb, db = raw(b, beta, frames)
            if ca != cb or min(da["inertia_relative_margin"], db["inertia_relative_margin"]) < POLICY["count_inertia_margin_min"]:
                raise ArithmeticError("Unresolved global root/pole coincidence")
            return ca, {"pole_exclusion": [a, b], "both_counts": [ca, cb], "sides": [da, db]}
        value, diagnostic = raw(w, beta, frames)
        if diagnostic["inertia_relative_margin"] < POLICY["count_inertia_margin_min"] or max(diagnostic["arm_boundary_conditions"]) > POLICY["pole_condition_max"]:
            raise ArithmeticError("Unresolved inertia sign/Dirichlet conditioning")
        return value, diagnostic
    return count


def solve_case(model, length, beta, count, upper=None, arm_frames=None):
    c0 = math.sqrt(model.section.E/model.section.rho)
    lower = POLICY["omega_min_pi_c0_L"]*math.pi*c0/length
    upper = POLICY["omega_max_pi_c0_L"]*math.pi*c0/length if upper is None else upper
    roots, search = joint.frame_roots(model, (length/2, length/2), beta, lower, upper, POLICY, count, arm_frames)
    if search["lower_count"] != 0:
        raise ArithmeticError("Hidden roots below search start")
    modes = []
    for n, root in enumerate(roots, 1):
        mode = joint.frame_mode(model, (length/2, length/2), beta, root["omega"], POLICY["quadrature_order"], arm_frames)
        d = mode["diagnostics"]
        if max(d["joint_residual_scaled"].values()) > POLICY["boundary_joint_scaled_tol"] or d["clamp_scaled_residual"] > POLICY["boundary_joint_scaled_tol"] or d["equation_scaled_residual"] > POLICY["equation_scaled_tol"] or d["energy_relative_error"] > POLICY["energy_relative_tol"] or d["nonzero_singular_condition"] > POLICY["nonzero_singular_condition_max"] or root["singular_ratio"] > POLICY["boundary_joint_scaled_tol"]:
            raise ArithmeticError(f"Mode residual gate failed: beta={beta}, n={n}, diagnostics={d}")
        root.update({"sorted_position": n, "diagnostics": d, "coefficients": mode["coefficients"].tolist()})
        modes.append(mode)
    nodes, weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    x, weights = (nodes+1)*length/4, weights*length/4
    p = model.coefficients
    masses = np.array([p["m"], p["j"], p["m"], p["r"]])
    fields = [[joint.arm_state(model, length/2, r["omega"], np.array(r["coefficients"])[i], x)[:, :4]
               for i in range(2)] for r in roots]
    gram = np.array([[sum(float(weights@((a*b)@masses)) for a, b in zip(first, second))
                      for second in fields] for first in fields])
    gram_error = float(np.max(np.abs(gram-np.eye(len(roots)))))
    if gram_error > POLICY["mass_gram_tol"]:
        raise ArithmeticError("Normalization/orthogonality gate failed")
    active_frames = joint.frames(beta) if arm_frames is None else arm_frames
    return {"beta_deg": beta, "lengths": [length/2, length/2],
        "arm_frames": [{"t": f.t, "n": f.n, "local_to_global": f.nodal_transform.tolist(),
                        "global_to_local": f.nodal_transform.T.tolist()} for f in active_frames],
        "joint_endpoint_signs": [1, 1], "state_order": list(joint.STATE_ORDER),
        "roots": roots, "search": search, "mass_gram": gram.tolist(),
        "mass_gram_max_error": gram_error, "status": "PASS"}


def symmetry_comparison(model, length, first, second, *, swapped=False, mirror=False):
    rows = []
    nodes, weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    x, weights = (nodes+1)*length/4, weights*length/4
    p = model.coefficients
    masses = np.array([p["m"], p["j"], p["m"], p["r"]])
    if len(first["roots"]) != len(second["roots"]):
        raise ArithmeticError("Symmetry inventories differ")
    for case in (first, second):
        if any((b["omega"]-a["omega"])/b["omega"] < POLICY["cluster_relative_gap"]
               for a, b in zip(case["roots"][:-1], case["roots"][1:])):
            raise ArithmeticError("Unresolved symmetry cluster: subspace comparison required, no forced vector match")
    for a, b in zip(first["roots"], second["roots"]):
        left, right = np.array(a["coefficients"]), np.array(b["coefficients"])
        afields = [joint.arm_state(model, length/2, a["omega"], left[i], x) for i in range(2)]
        bfields = [joint.arm_state(model, length/2, b["omega"], right[1-i if swapped else i], x)*
                   (joint.MIRROR_STATE if mirror else 1.) for i in range(2)]
        overlap = sum(float(weights@((v[:, :4]*r[:, :4])@masses)) for v, r in zip(afields, bfields))
        sign = 1. if overlap >= 0 else -1.
        errors = [math.sqrt(sum(float(weights@(v[:, i]-sign*r[:, i])**2) for v, r in zip(afields, bfields)) /
            max(sum(float(weights@v[:, i]**2) for v in afields), 1e-30)) for i in range(4)]
        diff = abs(a["omega"]/b["omega"]-1)
        if diff > POLICY["symmetry_frequency_relative_tol"] or 1-overlap**2 > POLICY["symmetry_MAC_loss_tol"] or max(errors) > POLICY["symmetry_component_L2_tol"]:
            raise ArithmeticError("Symmetry frequency/profile remapping failed")
        rows.append({"sorted_position": a["sorted_position"], "frequency_relative_difference": diff,
            "same_geometry_mass_MAC": min(1., overlap**2), "component_L2_relative": errors})
    return {"status": "PASS", "rows": rows}


def compute(config, model, length, direct, frozen):
    # These gates precede ANY nonzero-angle root calculation.
    structure = structural_checks(model, length)
    regression = beta0_regression(config, model, length, direct, frozen)
    catalog = pole_catalog(config, model, length)
    count = count_function(model, length, catalog)
    ceiling = 2*math.pi*POLICY["symmetry_gate_ceiling_fstar"]*math.sqrt(model.section.E/model.section.rho)/length
    baseline = solve_case(model, length, 0., count, ceiling)
    zero_limit = []
    for beta, allowance in zip(SMALL_ANGLES, POLICY["zero_limit_relative_allowances"]):
        # Authorized numerical continuity control, not angle-map production.
        case = solve_case(model, length, beta, count, ceiling)
        errors = [abs(a["omega"]/b["omega"]-1) for a, b in zip(case["roots"][:3], baseline["roots"][:3])]
        j0, jb = joint.joint_matrix(joint.frames(0.)), joint.joint_matrix(joint.frames(beta))
        matrix_difference = float(np.linalg.norm(jb-j0, 2))
        if len(errors) != 3 or max(errors) > allowance or matrix_difference > math.radians(beta)*(1+1e-10):
            raise ArithmeticError("Zero-limit hard gate failed")
        zero_limit.append({"beta_deg": beta, "first_three_relative_differences": errors,
            "joint_operator_norm_difference": matrix_difference, "allowance": allowance, "case": case})
    # §12/13 eigenpair controls are required BEFORE the full spectral pilot.
    canonical = solve_case(model, length, 45., count, ceiling)
    swapped = solve_case(model, length, 45., count, ceiling, joint.frames(45.)[::-1])
    mirrored = solve_case(model, length, -45., count, ceiling, joint.reflected_frames(joint.frames(45.)))
    swap = symmetry_comparison(model, length, canonical, swapped, swapped=True)
    reflection = symmetry_comparison(model, length, canonical, mirrored, mirror=True)
    statuses = {PREFIX+n: "PASS" for n in ("COORDINATES", "VIRTUAL_WORK_DUALITY", "BETA0_REGRESSION",
        "ZERO_LIMIT", "JOINT_RANK", "RIGHT_ANGLE_GEOMETRY", "ARM_SWAP", "REFLECTION")}
    # HARD barrier: all structural/symmetry gates pass before this point.
    pilot = []
    for beta in ANGLES:
        case = solve_case(model, length, beta, count)
        if len(case["roots"]) < POLICY["reported_positions"]+POLICY["guard_roots"]:
            raise ArithmeticError("Insufficient pilot guard coverage")
        pilot.append(case)
    statuses.update({PREFIX+"SPECTRAL_PILOT": "PASS", PREFIX+"JOINT_GATE": "PASS"})
    return {"statuses": statuses, "structure": structure, "beta0_regression": regression,
        "geometry_convention": {"source_contract": "docs/laminated_beams/reddy_inplane_coordinate_contract.md section3",
            "beta": "signed joint-to-right-clamp ray angle from +EX, positive toward +EY",
            "local_x": "both outer clamp0 -> jointLi", "rotation_axis": "k=-EZ",
            "joint_origin": [0, 0], "coordinate_fields": "u along t; w along n; c scalar; theta signed rotation about k"},
        "fixed_arm_poles": catalog, "zero_limit": zero_limit,
        "symmetry_gate": {"canonical": canonical, "swapped": swapped, "mirrored": mirrored,
            "arm_swap": swap, "reflection": reflection}, "pilot": pilot,
        "policy": POLICY, "coefficients": model.coefficients, "total_length": length,
        "scope": "sorted positions at three fixed angles; no tracking/hierarchy/applicability study",
        "closure_qualification": zero.SOURCE_CLOSURE["qualification"]}


def identity(config, checked):
    _, base = single.identity(config, checked)
    paths = (Path(__file__), ROOT/"scripts/lib/mindlin_herrmann_timoshenko_joint.py",
             ROOT/"scripts/lib/reddy_inplane_geometry.py")
    record = {"version": joint.GENERAL_VERSION, "arm_identity": base, "policy": POLICY,
        "pilot_angles": list(ANGLES), "small_angles": list(SMALL_ANGLES),
        "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        "frozen_beta0_manifest_hash": sha(FROZEN/"manifest.json"),
        "frozen_beta0_result_hash": sha(FROZEN/"result.json")}
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()[:16], record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if not args.check_sources and not args.compute:
        parser.error("Select --check-sources or --compute")
    config, model, length, checked, direct, frozen = check_inputs()
    structural_checks(model, length)
    if args.check_sources:
        print("Source/arm/frozen hashes and geometry-duality-rank checked; no roots for --check-sources")
    if not args.compute:
        return 0
    fingerprint, inputs = identity(config, checked)
    out = args.output_dir/fingerprint
    manifest_path = out/"manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text(encoding="utf-8"))
        if previous["identity"] != inputs or any(sha(out/p) != h for p, h in previous["artifact_hashes"].items()):
            raise ValueError("Changed/stale general-beta bundle")
        print("Verified general-beta bundle; zero root evaluations:", out)
        return 0
    out.mkdir(parents=True, exist_ok=True)
    write_json(out/"parameters.json", {"single_rod": config, "policy": POLICY,
        "section": vars(model.section), "lengths": [length/2, length/2], "angles": ANGLES})
    try:
        result = compute(config, model, length, direct, frozen)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
        write_json(out/"failure.json", {"status": "FAIL", "reason": str(exc), "identity": inputs,
            "failure_policy": POLICY["failure_policy"]})
        raise
    write_json(out/"result.json", result)
    frequency_rows, profiles = [], []
    for case in result["pilot"]:
        for r in case["roots"]:
            frequency_rows.append({"beta_deg": case["beta_deg"], "sorted_position": r["sorted_position"],
                "frequency_hz": r["frequency_hz"], "guard": r["sorted_position"] == 13,
                "in_reported_prefix": r["sorted_position"] <= 13})
            for arm, a in enumerate(r["coefficients"], 1):
                x = np.linspace(0., length/2, 201)
                fields = joint.arm_state(model, length/2, r["omega"], np.array(a), x)
                profiles.extend({"beta_deg": case["beta_deg"], "sorted_position": r["sorted_position"],
                    "arm": arm, "x_local": float(position), **dict(zip(joint.STATE_ORDER, state.tolist()))}
                    for position, state in zip(x, fields))
    write_csv(out/"frequencies.csv", frequency_rows)
    write_csv(out/"mode_profiles.csv", profiles)
    hashes = {p.name: sha(p) for p in out.iterdir() if p.name != "manifest.json" and p.is_file()}
    write_json(manifest_path, {"identity": inputs, "artifact_hashes": hashes,
        "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True, encoding="utf-8"),
        "spectrum_semantics": "independently sorted positions; same-geometry symmetry overlaps only; no across-beta MAC tracking",
        "source_qualification": zero.SOURCE_CLOSURE})
    write_json(args.output_dir/"current.json", {"fingerprint": fingerprint, "directory": str(out.resolve())})
    print("MHTIM_GENERAL_BETA_JOINT_GATE PASS:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
