"""Finite straight-rod Jang project gate; no frame, joint or parameter map.

New contract: finite essential boundaries, eigenvalues and completeness.
The source-reproduction CLI keeps its explicit Jang input semantics.
Reuses the M-H energy/dispersion helper, existing rectangular Timoshenko
section and Bishop's exact H=0 boundary system / atomic artifact writers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.analysis.reproduce_mindlin_herrmann_timoshenko_literature import source_check
from scripts.lib import bishop_longitudinal as bishop
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section

CONFIG = ROOT/"data/input/mindlin_herrmann_timoshenko_single_rod.json"
OUTPUT = ROOT/"results/mindlin_herrmann_timoshenko_single_rod"


def check_inputs():
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    source_path = ROOT/config["source_contract"]
    if sha(source_path) != config["source_contract_sha256"]:
        raise ValueError("Accepted rectangular contract changed; project kappa/geometry re-audit required")
    contract = json.loads(source_path.read_text(encoding="utf-8"))
    material, geometry = contract["material"], contract["geometries"][config["geometry_id"]]
    section = rectangular_section(E=material["E"], nu=material["nu"], rho=material["rho"],
        K=material["K"], width=geometry["width"], thickness=geometry["thickness"])
    model = mh.project_jang_reduced_rectangular(section)
    fixture, checked = source_check()
    if "fernandes" not in checked:
        raise ValueError("Fernandes source registration missing")
    return config, model, contract["lengths"]["L_ref"], checked


def independent_root(model, length, block, primary, policy):
    from scipy.optimize import brentq
    import numpy as np
    evaluations = 0
    def determinant(omega):
        nonlocal evaluations
        evaluations += 1
        matrix, _ = mh.transfer_boundary_matrix(model, length, omega, block,
            policy["qr_step_exponent_cap"], policy["qr_max_steps"])
        return float(np.linalg.det(matrix))
    left, right = primary["bracket_omega"]
    signs = [determinant(left), determinant(right)]
    if signs[0]*signs[1] > 0:
        raise ArithmeticError("Independent state determinant did not bracket the primary root")
    omega = brentq(determinant, left, right, xtol=policy["root_xtol"], rtol=policy["root_rtol"])
    matrix, transfer = mh.transfer_boundary_matrix(model, length, omega, block,
        policy["qr_step_exponent_cap"], policy["qr_max_steps"])
    singular = np.linalg.svd(matrix, compute_uv=False)
    return {"omega": omega, "frequency_hz": omega/(2*math.pi), "evaluations": evaluations,
        "relative_difference": abs(omega/primary["omega"]-1),
        "projected_boundary_singular_ratio": float(singular[-1]/singular[0]),
        "bracket_determinants": signs, **transfer}


def solve_block(model, length, block, config):
    import numpy as np
    policy = config["policy"]
    c0 = math.sqrt(model.section.E/model.section.rho)
    lower = policy["omega_min_over_pi_c0_L"]*math.pi*c0/length
    upper = policy["mh_omega_max_over_pi_c0_L"]*math.pi*c0/length if block == "mh" else math.sqrt(
        mh.blocks(model)[1].temporal(policy["timo_upper_ss_index"]*math.pi/length)[0]["omega_squared"])
    certificate = mh.finite_count_upper_bound(model, length, upper, block, policy["young_eta"])
    start_count = mh.finite_count_upper_bound(model, length, lower, block, policy["young_eta"])
    if start_count["upper_count"] != 0:
        raise ArithmeticError("Search lower boundary does not exclude hidden roots")
    attempts = []
    for attempt in range(2):
        active = {**policy, "scan_intervals": policy["scan_intervals"]*2**attempt}
        roots, search = mh.finite_roots(model, length, block, lower, upper, active)
        attempts.append(search)
        if len(roots) == certificate["upper_count"]:
            break
    if len(roots) != certificate["upper_count"]:
        raise ArithmeticError(f"Root count not certified: found {len(roots)}, bound {certificate['upper_count']}; attempts={attempts}")
    requested = config["reported_axial" if block == "mh" else "reported_bending"]
    if len(roots) < requested+config["guard_roots"]:
        raise ArithmeticError("Certified interval has insufficient guard roots")
    modes, profiles = [], []
    for index, root in enumerate(roots, 1):
        independent = independent_root(model, length, block, root, policy)
        mode = mh.finite_mode(model, length, root["omega"], block, policy["quadrature_order"])
        diagnostic = mode["diagnostics"]
        if (independent["relative_difference"] > policy["independent_frequency_relative_tol"] or
            root["singular_ratio"] > policy["boundary_scaled_residual_tol"] or
            root["nonzero_singular_condition"] > policy["nonzero_singular_condition_max"] or
            diagnostic["boundary_scaled_residual"] > policy["boundary_scaled_residual_tol"] or
            diagnostic["equation_scaled_residual"] > policy["equation_scaled_residual_tol"] or
            diagnostic["energy_relative_error"] > policy["energy_relative_tol"]):
            raise ArithmeticError(f"Finite mode gate failed: {block}/{index}: {root}, {independent}, {diagnostic}")
        root.update({"family_index": index, "family": "axial_acoustic" if block == "mh" else "bending",
            "guard": index > requested, "independent": independent, "diagnostics": diagnostic})
        modes.append(mode)
        points = np.linspace(0., length, 201)
        state = mh.finite_state_basis(model, length, root["omega"], points, block)@mode["coefficients"]
        profiles.extend({"family_index": index, "x": float(x), "first_displacement": float(y[0]),
            "second_coordinate": float(y[1]), "first_resultant": float(y[2]), "second_resultant": float(y[3])}
            for x, y in zip(points, state))
    p = model.coefficients
    mass2 = p["j"] if block == "mh" else p["r"]
    gram = [[float(modes[0]["weights"]@(p["m"]*a["values"][:, 0]*b["values"][:, 0]+
        mass2*a["values"][:, 1]*b["values"][:, 1])) for b in modes] for a in modes]
    error = float(np.max(np.abs(np.array(gram)-np.eye(len(modes)))))
    if error > policy["mass_orthogonality_tol"]:
        raise ArithmeticError("Finite mass orthogonality failed")
    return {"roots": roots, "search_attempts": attempts, "completeness": certificate,
        "below_search_count_bound": start_count, "mass_gram": gram,
        "mass_orthogonality_max_error": error, "status": "PASS"}, profiles


def hierarchy(model, length, axial, bending, config):
    import numpy as np
    s = model.section
    n = np.arange(1, config["reported_axial"]+config["guard_roots"]+1)
    k = n*math.pi/length
    frequencies = {"elementary": np.sqrt(s.EA*k*k/s.rhoA)/(2*math.pi),
        "rayleigh_love_planar": np.sqrt(s.EA*k*k/(s.rhoA+s.nu**2*s.rhoI*k*k))/(2*math.pi),
        "mindlin_herrmann": np.array([r["frequency_hz"] for r in axial["roots"]])}
    reduced_checks = []
    for name in ("elementary", "rayleigh_love_planar"):
        segment = bishop.Segment(length, s.EA, s.rhoA, H=0.,
            J=0. if name == "elementary" else s.nu**2*s.rhoI)
        for f in frequencies[name]:
            matrix, _ = bishop.boundary_matrix((segment,), f, ("U", "U"))
            ratio = np.linalg.svd(matrix, compute_uv=False)
            reduced_checks.append({"variant": name, "frequency_hz": float(f),
                "H": segment.H, "J": segment.J, "singular_ratio": float(ratio[-1]/ratio[0])})
    if max(r["singular_ratio"] for r in reduced_checks) > config["policy"]["boundary_scaled_residual_tol"]:
        raise ArithmeticError("Exact reduced-order axial boundary check failed")
    spectra = {}
    for name, f_axial in frequencies.items():
        roots = [{"frequency_hz": float(f), "family": "axial_acoustic", "family_index": i}
                 for i, f in enumerate(f_axial, 1)]
        roots += [{k: r[k] for k in ("frequency_hz", "family", "family_index")} for r in bending["roots"]]
        roots.sort(key=lambda r: r["frequency_hz"])
        prefix = roots[:config["reported_combined"]]
        cutoff = prefix[-1]["frequency_hz"]
        if cutoff >= min(f_axial[-1], bending["roots"][-1]["frequency_hz"]):
            raise ArithmeticError("Combined prefix exceeds a block's guard coverage")
        spectra[name] = [{"sorted_position": i, **r} for i, r in enumerate(prefix, 1)]
    return {"axial_frequencies_hz": {k: v[:config["reported_axial"]].tolist() for k, v in frequencies.items()},
        "bending_frequencies_hz": [r["frequency_hz"] for r in bending["roots"][:config["reported_bending"]]],
        "combined_sorted_prefix": spectra, "reduced_boundary_checks": reduced_checks,
        "numerical_status": "PASS", "status": "PARTIAL_PASS",
        "qualification": config["hierarchy_bc_qualification"],
        "rayleigh_love_scope": config["rayleigh_love_scope"]}


def compute(config, model, length):
    axial, axial_profiles = solve_block(model, length, "mh", config)
    bending, bending_profiles = solve_block(model, length, "timoshenko", config)
    p = model.coefficients
    result = {"statuses": {"PRODUCTION_MHTIM_FORMULATION_SELECTED": "PASS",
        "PRODUCTION_MHTIM_KAPPA": "RESOLVED", "MHTIM_SINGLE_ROD_FINITE_SPECTRUM": "PASS",
        "HIERARCHY_SINGLE_ROD": "PARTIAL_PASS"}, "formulation_status": "PRODUCTION_MHTIM_FORMULATION_SELECTED",
        "formulation": config["formulation"], "variant": model.variant,
        "kappa_status": "PRODUCTION_MHTIM_KAPPA_RESOLVED", "kappa": model.section.K,
        "finite_spectrum_status": "PASS", "hierarchy_status": "PARTIAL_PASS",
        "length": length, "normalization": "G20 normalized material E=rho=1; omega=2*pi*f; no experimental frequencies",
        "coefficients": p, "limits": mh.limits(model), "mh": axial, "timoshenko": bending,
        "hierarchy": hierarchy(model, length, axial, bending, config)}
    return result, axial_profiles, bending_profiles


def identity(config, checked):
    from importlib.metadata import version
    paths = [CONFIG, Path(__file__), ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py",
        ROOT/config["source_contract"], ROOT/"data/input/mindlin_herrmann_timoshenko_sources.json",
        ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
        ROOT/"scripts/lib/bishop_longitudinal.py"]
    record = {"version": mh.FINITE_ROD_VERSION, "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        "sources": checked, "python": platform.python_version(),
        "dependencies": {n: version(n) for n in ("numpy", "scipy")},
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()}
    fingerprint = hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()[:16]
    return fingerprint, record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if not args.check_sources and not args.compute:
        parser.error("Select --check-sources or --compute")
    config, model, length, checked = check_inputs()
    if args.check_sources:
        print("Sources verified; Jang closure selected; PROJECT rectangular kappa=5/6; source Jang remains explicit")
    if not args.compute:
        return 0
    fingerprint, inputs = identity(config, checked)
    directory = args.output_dir/fingerprint
    manifest_path = directory/"manifest.json"
    if manifest_path.exists():
        previous = json.loads(manifest_path.read_text())
        if previous["identity"] != inputs or any(sha(directory/p) != h for p, h in previous["artifact_hashes"].items()):
            raise ValueError("Stale or changed finite-rod bundle")
        print("Verified existing finite-rod bundle; zero root evaluations:", directory)
        return 0
    directory.mkdir(parents=True, exist_ok=True)
    try:
        result, axial_profiles, bending_profiles = compute(config, model, length)
    except (ValueError, ArithmeticError) as exc:
        write_json(directory/"failure.json", {"status": "FAIL", "reason": str(exc), "identity": inputs})
        raise
    write_json(directory/"result.json", result)
    write_json(directory/"parameters.json", config)
    for name, records in (("mh_modes.csv", axial_profiles), ("timoshenko_modes.csv", bending_profiles)):
        write_csv(directory/name, records)
    rows = [{"n": i+1, **{name: f[i] for name, f in result["hierarchy"]["axial_frequencies_hz"].items()},
             "bending": result["hierarchy"]["bending_frequencies_hz"][i]}
            for i in range(config["reported_axial"])]
    write_csv(directory/"hierarchy_frequencies.csv", rows)
    artifacts = {p.name: sha(p) for p in directory.iterdir() if p.is_file() and p.name != "manifest.json"}
    write_json(manifest_path, {"fingerprint": fingerprint, "identity": inputs, "artifact_hashes": artifacts,
        "command": subprocess.list2cmdline([sys.executable, *sys.argv]), "executable": sys.executable,
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True, encoding="utf-8"),
        "spectrum_semantics": "sorted family modes / combined positions at one fixed straight-rod input; no tracking"})
    write_json(args.output_dir/"current.json", {"fingerprint": fingerprint, "directory": str(directory.resolve())})
    print("Finite spectrum PASS; hierarchy PARTIAL_PASS (resolved-contraction clamp qualification):", directory)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
