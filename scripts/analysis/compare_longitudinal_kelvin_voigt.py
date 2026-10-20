"""Cache-only elastic comparison, then exactly four local elementary/RL KV roots.

Diagnostic-only new contract: the completed c-closure workflow cannot safely
serve as a preset because it searches full elastic prefixes and solves MH.
This entry point never calls those searches or modifies the source bundle.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scripts.analysis import check_contraction_joint_sensitivity as source
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import contraction_joint_sensitivity as old
from scripts.lib import longitudinal_kelvin_voigt_comparison as study

CONFIG = ROOT/"data/input/longitudinal_kelvin_voigt_comparison.json"
OUTPUT = ROOT/"results/longitudinal_kelvin_voigt_comparison"
VARIANTS = ("elementary", "rayleigh_love_planar", "mindlin_herrmann")


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def check_source(config):
    """Strict read-only provenance check; there is deliberately no fallback solve."""
    if (config["beta_deg"] != 45 or config["kappa_theta"] != 1 or config["d_theta"] != .001
        or config["variants"] != list(VARIANTS) or config["mh_closure"] != "FREE"
        or config["complex_positions"] != [1, 6] or config["elastic_positions"] != list(range(1, 7))):
        raise ValueError("Outside the authorized fixed G20 comparison")
    bundle = ROOT/config["source_bundle"]
    manifest, result, previous = read(bundle/"manifest.json"), read(bundle/"result.json"), read(source.CONFIG)
    model, _, sources, refs = source.references(previous)
    fingerprint, identity = source.identity(previous, sources, refs)
    if fingerprint != bundle.name or identity != manifest["identity"]:
        raise ValueError("Source identity/versions/hash mismatch; no automatic recomputation")
    for name, digest in manifest["artifact_hashes"].items():
        if sha(bundle/name) != digest:
            raise ValueError(f"Source artifact mismatch: {name}")
    if result["parameters"] != previous or (result["status"], result["elastic_status"], result["complex_status"]) != ("COMPLETE", "PASS", "PASS"):
        raise ValueError("Source parameters or acceptance do not match")
    if previous["policy"]["d_theta"] != config["d_theta"] or result["coefficients"] != model.coefficients:
        raise ValueError("Damper/coefficient mismatch")
    c = model.coefficients
    k = c["B"]/.5
    tref = .5**2*np.sqrt(c["m"]/c["B"])
    ctheta = .001*c["B"]*tref/.5
    if any(not np.isclose(result[name], value, rtol=1e-14, atol=0.)
           for name, value in (("k_theta", k), ("t_ref", tref), ("c_theta", ctheta))):
        raise ValueError("Physical reference normalizations differ")
    protection = {str((bundle/name).relative_to(ROOT)).replace("\\", "/"): sha(bundle/name)
                  for name in ["manifest.json", *manifest["artifact_hashes"]]}
    protection.update(identity["files"])
    for name in ("docs/theory/contraction_joint_sensitivity.md", "results/contraction_joint_sensitivity/README.md",
                 "results/contraction_joint_sensitivity/current.json"):
        protection[name] = sha(ROOT/name)
    return model, result, protection


def make_identity(config, protection):
    paths = (CONFIG, Path(__file__), Path(study.__file__))
    identity = {"config": config, "source_hashes": protection,
        "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}}
    return hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16], identity


def records(previous):
    return {v: next(c for c in previous["cases"] if c["case"]["beta"] == 45
                   and c["case"]["variant"] == v and c["case"]["rotation"] == "SPRING"
                   and (v != "mindlin_herrmann" or c["case"]["closure"] == "FREE")) for v in VARIANTS}


def stage_a(model, previous, config):
    policy, lengths = config["policy"], previous["parameters"]["geometry"]["lengths"]
    catalog, diags, table, failures = records(previous), {}, [], []
    for variant, item in catalog.items():
        case = old.Case(**item["case"])
        diags[variant] = []
        for root in item["roots"][:6]:
            d = study.saved_elastic_diagnostics(model, lengths, case, root, previous["k_theta"], policy["quadrature_order"])
            diags[variant].append(d)
            energy_error = abs(d["elastic_energy"]/(root["omega"]**2*d["mass"])-1)
            mass_error = abs(d["mass"]-1)
            saved_mass_error = abs(d["mass"]-root["diagnostics"]["mass"])
            if (mass_error > policy["mass_normalization_tol"] or saved_mass_error > policy["mass_normalization_tol"]
                or energy_error > policy["energy_tol"] or max(d["joint_scaled"]) > policy["boundary_tol"]
                or d["clamp_scaled"] > policy["boundary_tol"] or d["ode_scaled"] > policy["ode_tol"]):
                failures.append(f"elastic diagnostics {variant} k{root['sorted_position']}")
            omega, mass, jump = root["omega"], d["mass"], float(d["delta_theta"])
            alpha = previous["c_theta"]*abs(jump)**2/(2*mass)
            fraction = d["spring_energy"]/d["elastic_energy"]
            common = mass-d["contraction_or_gradient_mass"]
            table.append({"variant": variant, "k": root["sorted_position"], "omega_0": omega,
                "delta_theta": jump, "mass": mass, "common_mass": common,
                "additional_mass": d["contraction_or_gradient_mass"],
                "spring_potential_energy": .5*d["spring_energy"], "spring_energy_fraction": fraction,
                "active": fraction > policy["active_spring_energy_fraction_min"],
                "alpha_pred": alpha, "zeta_pred": alpha/omega,
                "rotation_participation_common": abs(jump)**2/common, "mass_correction": common/mass,
                "mass_error": mass_error, "energy_relative": energy_error,
                "joint_scaled_max": float(max(d["joint_scaled"])), "clamp_scaled": d["clamp_scaled"], "ode_scaled": d["ode_scaled"]})
    correspondence = {}
    ref_case = old.Case(**catalog["mindlin_herrmann"]["case"])
    for variant in VARIANTS[:2]:
        correspondence[variant] = study.projection_correspondence(model, old.Case(**catalog[variant]["case"]),
            ref_case, diags[variant], diags["mindlin_herrmann"], policy)
        if correspondence[variant]["status"] != "PASS":
            failures.append("correspondence "+variant)
    for row in table:
        ref = next(r for r in table if r["variant"] == "mindlin_herrmann" and r["k"] == row["k"])
        row["omega_vs_MH_percent"] = 100*(row["omega_0"]/ref["omega_0"]-1)
        active = row["active"] and ref["active"]
        for name in ("alpha_pred", "zeta_pred", "rotation_participation_common", "mass_correction"):
            row[name+"_vs_MH_percent"] = 100*(row[name]/ref[name]-1) if active else None
    for variant in VARIANTS:
        frequencies = np.array([r["omega"] for r in catalog[variant]["roots"]])
        for k in config["complex_positions"]:
            w = frequencies[k-1]
            row = next(r for r in table if r["variant"] == variant and r["k"] == k)
            if not row["active"] or min(abs(np.delete(frequencies, k-1)-w))/w <= policy["simple_relative_gap_min"]:
                failures.append(f"inactive or unresolved target {variant} k{k}")
    return {"status": "PASS" if not failures else "NUMERICAL_UNRESOLVED", "failures": failures,
            "table": table, "correspondence": correspondence}, diags


def scalar_complex_table(result):
    table = []
    for value in result["complex_modes"]:
        if value["status"] != "PASS":
            continue
        p = value["p"]
        p = complex(p["real"], p["imag"]) if isinstance(p, dict) else p
        row = {"variant": value["variant"], "k": value["k"], "provenance": value["provenance"],
               "omega_0": value["omega_elastic"], "omega_d": p.imag, "alpha": -p.real,
               "zeta": -p.real/abs(p), "alpha_pred": value["alpha_predictor"],
               "predictor_error_percent": 100*value["predictor_relative_error"]}
        table.append(row)
    for row in table:
        ref = next((r for r in table if r["variant"] == "mindlin_herrmann" and r["k"] == row["k"]), None)
        for name in ("omega_0", "omega_d", "alpha", "zeta", "alpha_pred"):
            row[name+"_vs_MH_percent"] = 100*(row[name]/ref[name]-1) if ref else None
    return table


def stage_b(model, previous, config, result, diags, out):
    if result["stage_a"]["status"] != "PASS":
        raise ArithmeticError("Stage A hard gate; complex stage prohibited")
    policy, catalog = config["policy"], records(previous)
    lengths = previous["parameters"]["geometry"]["lengths"]
    providers, result["preflight"] = {}, []
    for variant in VARIANTS[:2]:
        for k in config["complex_positions"]:
            omega = catalog[variant]["roots"][k-1]["omega"]
            provider = study.ReducedComplexBasis(model, lengths, omega, old.Case(**catalog[variant]["case"]), previous["k_theta"], previous["t_ref"])
            providers[variant, k] = provider
            check = study.provider_preflight(provider, previous["c_theta"], policy)
            result["preflight"].append({"variant": variant, "k": k, **check})
    if any(c["status"] != "PASS" for c in result["preflight"]):
        raise ArithmeticError("Complex provider preflight failed before all correctors")
    for k in config["complex_positions"]:
        value = next(v for v in previous["complex_modes"] if v["case"] == catalog["mindlin_herrmann"]["key"] and v["k"] == k)
        if value["status"] != "PASS" or value["energy_identity_relative"] > policy["energy_tol"] or value["alpha_identity_relative"] > policy["energy_tol"]:
            raise ArithmeticError("Unaccepted source MH complex root")
        result["complex_modes"].append({**value, "variant": "mindlin_herrmann", "provenance": "REUSED_HASH_VERIFIED"})
        result["solve_ledger"]["reused_MH_complex_roots"] += 1
    profiles = []
    for (variant, k), provider in providers.items():
        print(f"local complex target {variant} k{k}", flush=True)
        row = next(r for r in result["stage_a"]["table"] if r["variant"] == variant and r["k"] == k)
        frequencies = np.array([r["omega"] for r in catalog[variant]["roots"]])
        gap = min(abs(np.delete(frequencies, k-1)-provider.omega))
        value, d = study.complex_target(provider, row["alpha_pred"], diags[variant], k-1, gap,
            previous["c_theta"], policy, result["solve_ledger"])
        value.update(variant=variant, k=k, omega_elastic=provider.omega, provenance="NEW_LOCAL_COMPLEX")
        if d is not None:
            value["diagnostics"] = source.diagnostics_record(d)
            for arm, length in enumerate(lengths):
                x = np.linspace(0, length, 201)
                fields = provider.values(value["p"], arm, x)@value["coefficients"][arm]
                for xx, yy in zip(x, fields):
                    profiles.append({"variant": variant, "k": k, "arm": arm+1, "x": xx,
                        **{name+suffix: float(getattr(y, component)) for name, y in zip(("u", "w", "theta", "N", "Q", "M"), yy)
                           for suffix, component in (("_real", "real"), ("_imag", "imag"))}})
        result["complex_modes"].append(value)
        write_json(out/"result.json", source.plain(result))
        if profiles:
            write_csv(out/"new_complex_profiles.csv", profiles)
        if value["status"] != "PASS":
            raise ArithmeticError(f"NUMERICAL_UNRESOLVED: {variant} k{k}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--stage-a", action="store_true", help="saved-form evaluation only; zero correctors")
    args = parser.parse_args(argv)
    config = read(CONFIG)
    model, previous, protection = check_source(config)
    fingerprint, identity = make_identity(config, protection)
    out = OUTPUT/fingerprint
    if (out/"manifest.json").exists():
        manifest = read(out/"manifest.json")
        if manifest["identity"] != identity or any(sha(out/n) != h for n, h in manifest["artifact_hashes"].items()):
            raise ValueError("New bundle hash mismatch; will not overwrite or re-solve")
        cached = read(out/"result.json")
        if cached["status"] == "COMPLETE" or not args.compute or args.stage_a:
            print(json.dumps({"status": cached["status"], "fingerprint": fingerprint, "cache": "VALIDATED", "new_root_correctors_this_invocation": 0}))
            return
        if cached["status"] == "NUMERICAL_UNRESOLVED":
            raise RuntimeError("Unresolved bundle retained; no automatic retry")
    if not (args.compute or args.stage_a):
        print(json.dumps({"source": "PASS", "fingerprint": fingerprint, "new_root_correctors_this_invocation": 0}))
        return
    out.mkdir(parents=True, exist_ok=True)
    result = {"status": "RUNNING", "parameters": config, "coefficients": model.coefficients,
        **{k: previous[k] for k in ("k_theta", "c_theta", "t_ref")}, "complex_modes": [],
        "solve_ledger": {"elastic_spectrum_searches": 0, "elastic_shape_solves": 0,
            "reused_elastic_modes": 18, "reused_MH_complex_roots": 0,
            "elastic_provider_correctors": 0, "new_damped_correctors": 0}}
    try:
        result["stage_a"], diags = stage_a(model, previous, config)
        write_csv(out/"elastic_predictors.csv", result["stage_a"]["table"])
        write_json(out/"result.json", source.plain(result))
        if result["stage_a"]["status"] != "PASS":
            raise ArithmeticError("Stage A failed: "+str(result["stage_a"]["failures"]))
        result["status"] = "STAGE_A_PASS"
        if args.compute and not args.stage_a:
            stage_b(model, previous, config, result, diags, out)
            result["complex_comparison"] = scalar_complex_table(result)
            write_csv(out/"complex_comparison.csv", result["complex_comparison"])
            result["status"] = "COMPLETE"
    except Exception as error:
        result.update(status="NUMERICAL_UNRESOLVED", stop_reason=f"{type(error).__name__}: {error}")
    result["source_preservation"] = "PASS" if all(sha(ROOT/n) == h for n, h in protection.items()) else "FAIL"
    if result["source_preservation"] != "PASS":
        result["status"] = "NUMERICAL_UNRESOLVED"
    write_json(out/"result.json", source.plain(result))
    manifest = {"identity": identity, "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
        "head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "artifact_hashes": {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"}}
    write_json(out/"manifest.json", manifest)
    write_json(OUTPUT/"current.json", {"fingerprint": fingerprint, "status": result["status"]})
    print(json.dumps({"status": result["status"], "fingerprint": fingerprint,
                      "ledger": result["solve_ledger"], "stop_reason": result.get("stop_reason")}))
    if result["status"] == "NUMERICAL_UNRESOLVED":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
