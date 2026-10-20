"""One-axis KV spatial linear split: exact quadratic audit and one OUT prefix.

New diagnostic contract: fourteen vector joint conditions and operator split,
not a preset of a nonlinear evolution or planar frequency-map workflow.
Historical data, linear arms, c-continuous rigid joint and KV roots are read-only.
"""
from __future__ import annotations

import argparse
import hashlib
import inspect
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scripts.analysis import compare_longitudinal_kelvin_voigt as kv_run
from scripts.analysis import verify_weakly_nonlinear_spatial_rod as old_spatial
from scripts.analysis import check_contraction_joint_sensitivity as old_planar
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import uniaxial_kelvin_voigt_split as study
from scripts.lib.weakly_nonlinear_coupled_spatial_dynamics import rectangular_coefficients

CONFIG = ROOT/"data/input/uniaxial_kelvin_voigt_split.json"
OUTPUT = ROOT/"results/uniaxial_kelvin_voigt_split"


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def sources(config):
    if (config["structural_beta_deg"] != [0, 45, 90] or config["spectral_beta_deg"] != 45
        or config["lengths"] != [.5, .5] or config["contraction_closure"] != "FREE"
        or config["outer_clamp"] != "independent_section_rotation"):
        raise ValueError("Outside the declared linear diagnostic")
    protection = {}
    directory = ROOT/config["spatial_bundle"]
    manifest, spatial = read(directory/"manifest.json"), read(directory/"result.json")
    _, identity, fingerprint = old_spatial.identity()
    if identity != manifest["identity"] or fingerprint != directory.name:
        raise ValueError("Spatial source identity mismatch; no recomputation")
    if any(v != "PASS" for v in spatial["statuses"].values() if v != "MATCH"):
        raise ValueError("Unaccepted spatial model audit")
    for name, digest in manifest["artifacts"].items():
        if sha(directory/name) != digest:
            raise ValueError("Spatial artifact hash mismatch: "+name)
        protection[(directory/name).relative_to(ROOT).as_posix()] = digest
    protection[(directory/"manifest.json").relative_to(ROOT).as_posix()] = sha(directory/"manifest.json")
    protection.update({Path(n).as_posix(): h for n, h in identity["hashes"].items()})
    config_kv = read(kv_run.CONFIG)
    mh_model, _, old_hashes = kv_run.check_source(config_kv)
    kv_fingerprint, kv_identity = kv_run.make_identity(config_kv, old_hashes)
    directory = ROOT/config["kv_bundle"]
    manifest, kv = read(directory/"manifest.json"), read(directory/"result.json")
    if kv_fingerprint != directory.name or manifest["identity"] != kv_identity or kv["status"] != "COMPLETE":
        raise ValueError("KV source identity/status mismatch")
    for name, digest in manifest["artifact_hashes"].items():
        if sha(directory/name) != digest:
            raise ValueError("KV artifact hash mismatch: "+name)
        protection[(directory/name).relative_to(ROOT).as_posix()] = digest
    protection[(directory/"manifest.json").relative_to(ROOT).as_posix()] = sha(directory/"manifest.json")
    protection.update(old_hashes)
    protection.update(kv_identity["files"])
    for name in ("docs/theory/weakly_nonlinear_spatial_rod.md", "docs/theory/out_of_plane_eb_torsion.md",
                 "docs/anisotropic_rods/yartsev_ch2_rigid_angular_joint.md", "docs/theory/mindlin_herrmann_timoshenko_rigid_joint.md",
                 "docs/theory/longitudinal_kelvin_voigt_comparison.md", "docs/laminated_beams/inplane_kelvin_voigt_joint_theory.tex",
                 "scripts/lib/weakly_nonlinear_coupled_spatial_dynamics.py", "scripts/lib/yartsev_ch2_coupled_rods.py"):
        protection[name] = sha(ROOT/name)
    coefficient, torsion = rectangular_coefficients(.20, .05)
    if coefficient.values() != spatial["coefficients"]:
        raise ValueError("Current generalized rectangular torsion differs from saved coefficients")
    return spatial, kv, mh_model, protection, torsion


def identity(config, protection):
    paths = (CONFIG, Path(__file__), Path(study.__file__))
    record = {"config": config, "source_hashes": protection,
        "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}}
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()[:16], record


def spectral(outplane, kv, config, result, out, reuse_prefix=None):
    policy = config["policy"]
    c = outplane.c
    trial = 4*np.pi*np.sqrt(c.CT/(c.jp+c.jb))/.5
    cutoff = 2*np.pi*study.mh.blocks(outplane.tim)[1].cutoff_hz
    if not trial < policy["omega_upper"] < cutoff:
        raise ArithmeticError("Preset ceiling is not supported by trial bound/bounded basis")
    result["spectral"] = {"status": "RUNNING", "seventh_trial_bound": trial, "optical_cutoff": cutoff}
    reused = None
    if reuse_prefix is not None:
        directory = Path(reuse_prefix).resolve()
        manifest = read(directory/"manifest.json")
        reused = read(directory/"result.json")
        certificate = read(directory/"prefix_reuse_certificate.json")
        if any(sha(directory/n) != h for n, h in manifest["artifact_hashes"].items()):
            raise ValueError("Prefix source artifact mismatch")
        if reused["parameters"] != config or reused["coefficients"] != result["coefficients"] or reused["spectral"]["search"]["status"] != "PASS":
            raise ValueError("Prefix source physical contract/count mismatch")
        if sha(directory/"result.json") != certificate["source_result_sha256"]:
            raise ValueError("Prefix result differs from its reuse certificate")
        for name, digest in certificate["method_hashes"].items():
            actual = hashlib.sha256(inspect.getsource(getattr(study.OutOfPlane, name)).encode()).hexdigest()
            if actual != digest:
                raise ValueError("Root-search implementation changed: "+name)
        if any(sha(ROOT/n) != h for n, h in manifest["identity"]["source_hashes"].items()):
            raise ValueError("Prefix source/model hashes changed")
        result["prefix_reuse"] = {"directory": str(directory), "reason": certificate["reason"],
            "source_files": {str(directory/n): sha(directory/n) for n in ("manifest.json", "result.json", "prefix_reuse_certificate.json")},
            "unchanged_method_hashes": certificate["method_hashes"], "original_solve_ledger": reused["solve_ledger"]}
        catalog = reused["spectral"]["catalog"]
    else:
        catalog = outplane.pole_catalog(policy)
        result["solve_ledger"]["fixed_arm_bending_catalog_scans"] += 1
        result["solve_ledger"]["fixed_arm_bending_poles"] = len(catalog["tim_roots"])
    result["spectral"]["catalog"] = catalog
    if catalog["status"] != "PASS":
        raise ArithmeticError("Fixed-arm pole completeness is not certified at the declared ceiling")
    green = [{"omega": w, "relative": outplane.green_check(w)} for w in (.317, 2.3)]
    result["spectral"]["schur_green_identity"] = green
    if max(g["relative"] for g in green) > policy["energy_tol"]:
        raise ArithmeticError("Schur Green identity failed")
    if reused is None:
        result["solve_ledger"]["global_prefix_searches"] += 1
        search = result["spectral"]["search"] = {}
        roots = outplane.prefix(catalog, policy, search)
        result["solve_ledger"]["global_out_roots"] = len(roots)
    else:
        result["spectral"]["search"] = reused["spectral"]["search"]
        roots = reused["spectral"]["search"]["roots"]
        result["solve_ledger"]["reused_out_roots"] = len(roots)
    result["spectral"]["modes"] = []
    rows, profiles, vectors = [], [], []
    for k, root in enumerate(roots, 1):
        omega = root["omega"]
        coeff, d = outplane.mode(omega, kv["k_theta"]+1j*omega*kv["c_theta"], policy)
        vectors.append((d["values"][:, :, :3]*np.sqrt(d["weights"][:, None]*[c.m, c.jb, c.jp+c.jb])).ravel())
        diag = {n: v for n, v in d.items() if n not in ("x", "weights", "values", "gradients")}
        result["spectral"]["modes"].append({"sorted_position": k, "guard": k == 7, "omega": omega,
            "p": 1j*omega, "coefficients": coeff, "diagnostics": diag})
        rows.append({"k": k, "guard": k == 7, "omega": omega, "alpha": 0., "mass": d["mass"],
            "joint_scaled_max": max(d["joint_scaled"]), "clamp_scaled": d["clamp_scaled"],
            "ode_scaled": d["ode_scaled"], "energy_relative": d["energy_relative"],
            "rigid_KV_matrix_max": d["rigid_kv_boundary_max"], "delta_theta": 0.})
        x = np.linspace(0., .5, 201)
        for arm in range(2):
            values = outplane.basis(omega, x)@coeff[arm]
            for xx, yy in zip(x, values):
                profiles.append({"k": k, "arm": arm+1, "x": xx,
                    **dict(zip(("v", "psi", "Phi", "Qv", "Mpsi", "MT"), yy))})
    gram = np.array(vectors)@np.array(vectors).T
    result["spectral"].update(mass_gram=gram, mass_gram_error=float(np.max(abs(gram-np.eye(7)))))
    write_csv(out/"frequencies.csv", rows)
    write_csv(out/"profiles.csv", profiles)
    if result["spectral"]["mass_gram_error"] > policy["mass_gram_tol"]:
        raise ArithmeticError("Mass orthogonality gate failed")
    result["spectral"]["status"] = "PASS"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--structure-only", action="store_true")
    parser.add_argument("--reuse-prefix", help="Hash/method-certified local prefix; physical checks only, zero root searches")
    args = parser.parse_args(argv)
    config = read(CONFIG)
    spatial, kv, mh_model, protection, torsion = sources(config)
    fingerprint, ident = identity(config, protection)
    out = OUTPUT/fingerprint
    if (out/"manifest.json").exists():
        manifest = read(out/"manifest.json")
        if manifest["identity"] != ident or any(sha(out/n) != h for n, h in manifest["artifact_hashes"].items()):
            raise ValueError("Diagnostic cache hash mismatch; no overwrite or retry")
        result = read(out/"result.json")
        if any(sha(Path(n)) != h for n, h in result.get("prefix_reuse", {}).get("source_files", {}).items()):
            raise ValueError("Reused prefix provenance changed")
        if result["status"] != "STRUCTURE_PASS" or not args.compute or args.structure_only:
            print(json.dumps({"status": result["status"], "cache": "VALIDATED", "fingerprint": fingerprint,
                              "new_spectral_solves_this_invocation": 0}))
            return
    if not (args.compute or args.structure_only):
        print(json.dumps({"sources": "PASS", "fingerprint": fingerprint, "new_spectral_solves_this_invocation": 0}))
        return
    out.mkdir(parents=True, exist_ok=True)
    result = {"status": "RUNNING", "parameters": config, "torsion_provenance": torsion,
              "k_theta": kv["k_theta"], "c_theta": kv["c_theta"], "t_ref": kv["t_ref"],
              "solve_ledger": {"fixed_arm_bending_catalog_scans": 0, "fixed_arm_bending_poles": 0,
                  "global_prefix_searches": 0, "global_out_roots": 0,
                  "repeated_KV_spectrum_searches": 0, "inplane_root_solves": 0, "complex_correctors": 0}}
    try:
        coefficients, operator, action = study.quadratic_action(spatial)
        result.update(coefficients=coefficients.values(), quadratic_action=action)
        result["structure"] = study.structural_checks(operator, coefficients, mh_model, kv, config["policy"])
        result["status"] = "STRUCTURE_PASS"
        write_json(out/"result.json", old_planar.plain(result))
        if args.compute and not args.structure_only:
            spectral(study.OutOfPlane(coefficients, operator), kv, config, result, out, args.reuse_prefix)
            result["status"] = "COMPLETE"
    except Exception as error:
        result["status"] = "SPECTRAL_UNRESOLVED" if result.get("structure", {}).get("status") == "PASS" else "STRUCTURE_UNRESOLVED"
        result["stop_reason"] = f"{type(error).__name__}: {error}"
        if "spectral" in result:
            result["spectral"]["status"] = "SPECTRAL_UNRESOLVED"
    result["source_preservation"] = "PASS" if all(sha(ROOT/n) == h for n, h in protection.items()) else "FAIL"
    if result["source_preservation"] != "PASS":
        result["status"] = "PRESERVATION_FAILED"
    write_json(out/"result.json", old_planar.plain(result))
    manifest = {"identity": ident, "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "artifact_hashes": {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"}}
    write_json(out/"manifest.json", manifest)
    write_json(OUTPUT/"current.json", {"fingerprint": fingerprint, "status": result["status"]})
    print(json.dumps({"status": result["status"], "fingerprint": fingerprint, "solve_ledger": result["solve_ledger"],
                      "stop_reason": result.get("stop_reason")}))
    if result["status"] not in ("COMPLETE", "STRUCTURE_PASS"):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
