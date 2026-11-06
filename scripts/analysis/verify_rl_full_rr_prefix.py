"""Bounded RL_FULL/RR prefix verification: stages A and B only.

Script Proliferation Control: event/slot/shape output and immutable old source
fingerprints require a separate diagnostic entry point, not another old preset.
Reuses frozen physics, count_at, tim_catalog and Ritz; no physical duplication.
--tables-only and --verify-cache never call a root solver.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scripts.lib import rl_full_rr_prefix as engine
from scripts.analysis.screen_coupled_longitudinal_theory_hierarchy import tim_catalog

p = engine.physics
CONFIG = ROOT/"data/input/rl_full_rr_prefix_ab.json"
OUTPUT = ROOT/"results/rl_full_rr_prefix"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_new(path, record):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(jsonable(record), stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def validate_scope(config):
    expected = dict(E=1., rho=1., l=1., nu=.3, epsilon=.005, aspect=4., ks=5/6)
    stages = {"A": dict(K_plot=2, K_guard=3, angles=[0., 6.1194, 30., 90.]),
              "B": dict(K_plot=8, K_guard=9, angles=[0., 5., 6.1194, 7., 30., 90.])}
    if (config["material"] != expected or config["stages"] != stages or
        (config["mu"], config["tau"], config["k_joint"], config["k_end"]) != (.01, 0., "RIGID", "RIGID")):
        raise ValueError("This CLI is restricted to the approved A-B configuration")


def verify_point_record(record):
    body = {k: v for k, v in record.items() if k != "record_digest"}
    if record.get("record_digest") != engine.digest(body):
        raise engine.PrefixFailure("POINT_RECORD_DIGEST_MISMATCH")


def verify_reference(folder):
    manifest = json.loads((folder/"manifest.json").read_text())
    failed = [name for name, h in manifest["artifacts"].items() if sha(folder/name) != h]
    if manifest["status"] != "PASS" or failed:
        raise engine.PrefixFailure("REFERENCE_CHANGED", failed=failed)
    return {"manifest_sha256": sha(folder/"manifest.json"), "verified_artifacts": len(manifest["artifacts"])}


def identity(config, reference):
    files = {CONFIG, Path(__file__), Path(engine.__file__),
             ROOT/"tests/test_rl_full_rr_prefix.py"}
    # Full loaded local dependency closure, not merely a Git revision.
    for mod in tuple(sys.modules.values()):
        name = getattr(mod, "__file__", None)
        if name:
            path = Path(name).resolve()
            if path.is_relative_to(ROOT) and path.suffix == ".py":
                files.add(path)
    record = dict(schema=config["schema"], engine=engine.VERSION, config=config,
                  sources={x.relative_to(ROOT).as_posix(): sha(x) for x in sorted(files)},
                  reference=reference, head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    return engine.digest(record)[:16], record


def load_modes(folder, record):
    file = folder/record["shape_file"]
    if sha(file) != record["shape_sha256"]:
        raise engine.PrefixFailure("SHAPE_HASH_MISMATCH", path=str(file))
    with np.load(file, allow_pickle=False) as arrays:
        return [dict(arms=[{field: arrays[f"mode{i+1}_arm{j+1}_{field}"]
                    for field in ("x", "weights", "y", "dy", "ends")} for j in range(2)])
                for i in range(len(record["roots"]))]


def catalogues(model, case, upper, policy):
    result, attempts = {}, []
    for L in case.lengths(model):
        if L in result:
            continue
        for intervals in (policy["pole_scan_intervals"], 2*policy["pole_scan_intervals"]):
            # Inspect the existing catalog's proposed ceiling before basis calls.
            dispersion = p.tim.blocks(model.carrier)[1]
            if upper >= 2*math.pi*dispersion.cutoff_hz:
                raise engine.PrefixFailure("TIM_BASIS_DOMAIN_BLOCKED", upper=upper)
            k = dispersion.spatial(upper/(2*math.pi))[0]["wavenumber_per_m"]
            index = math.ceil(k*L/math.pi-.9)+.9
            ceiling = math.sqrt(dispersion.temporal(index*math.pi/L)[0]["omega_squared"])
            limit = min(math.sqrt(model.section.KGA/model.section.rhoI), math.sqrt(model.section.EA/model.J))
            if ceiling >= limit:
                raise engine.PrefixFailure("CATALOG_CEILING_OUTSIDE_DOMAIN", ceiling=ceiling, limit=limit, length=L)
            try:
                cat = tim_catalog(model.carrier, L, upper, {"policy": {
                    "scan_intervals": intervals, "root_xtol": policy["root_xtol"], "root_rtol": policy["root_rtol"]}})
                cat["arm_contract"] = dict(length=L, material=asdict(model),
                    EA=model.section.EA, EI=model.section.EI, KGA=model.section.KGA,
                    rhoA=model.section.rhoA, rhoI=model.section.rhoI, J=model.J)
                result[L] = cat
                attempts.append(dict(length=L, intervals=intervals, status="PASS"))
                break
            except (ArithmeticError, ValueError) as error:
                attempts.append(dict(length=L, intervals=intervals, status="FAILED", reason=str(error)))
        else:
            raise engine.PrefixFailure("POLE_CATALOG_UNRESOLVED", attempts=attempts)
    return result, attempts


def ritz_checks(model, case, roots, initial, policy):
    target = np.array([r["omega"] for r in roots])
    trials = list(initial)
    for order in [*policy["ritz_orders"], policy["ritz_extra_order"]]:
        prior = next((t for t in trials if t["order"] == order), None)
        if prior is None:
            values, _, diag = p.ritz(model, case, order, count=len(target))
            trials.append(dict(order=order, omega=values.tolist(), diagnostics=diag))
        trials.sort(key=lambda t: t["order"])
        current = next(t for t in trials if t["order"] == order)
        current["relative_frequency_error"] = float(max(abs(np.array(current["omega"])/target-1)))
        predecessors = [t for t in trials if t["order"] < order]
        change = float(max(abs(np.array(current["omega"])/np.array(predecessors[-1]["omega"])-1))) if predecessors else None
        current["relative_p_change"] = change
        d = current["diagnostics"]
        current["accepted"] = bool(change is not None and change <= policy["ritz_frequency_tol"] and
            current["relative_frequency_error"] <= policy["ritz_frequency_tol"] and
            d["mass_gram_error"] <= policy["mass_gram_tol"] and d["energy_error"] <= policy["energy_tol"])
        if current["accepted"]:
            return dict(status="PASS", accepted_order=order, trials=trials)
    return dict(status="RITZ_NOT_CONVERGED", trials=trials)


def quadrature_check(model, case, roots, modes, policy):
    high, changes = [], []
    adapter = SimpleNamespace(arms=(model, model), section=model.section, l=model.l)
    for root, mode in zip(roots, modes):
        arms = engine.fields(model, case, root["omega"], mode["coefficients"], policy["quadrature_check_order"])
        mass, elastic, springs = engine.energies(adapter, case, arms)
        mass0, elastic0, _ = engine.energies(adapter, case, mode["arms"])
        error = dict(mass_change=float(abs(mass.sum()-mass0.sum())),
                     energy_relative_change=float(abs(elastic.sum()/elastic0.sum()-1)),
                     high_energy_error=float(abs((elastic.sum()+sum(springs))/mass.sum()/root["omega"]**2-1)))
        high.append(dict(arms=arms)); changes.append(error)
    gram = np.array([[p.inner(model, a, b) for b in high] for a in high])
    err = float(np.max(abs(gram-np.eye(len(high)))))
    passed = err <= policy["mass_gram_tol"] and all(
        c["mass_change"] <= policy["mass_gram_tol"] and c["energy_relative_change"] <= policy["energy_tol"] and
        c["high_energy_error"] <= policy["energy_tol"] for c in changes)
    return dict(status="PASS" if passed else "QUADRATURE_FAILED", low_order=policy["quadrature_order"],
                high_order=policy["quadrature_check_order"], gram_max_error=err, modes=changes)


def regression(model, case, roots, modes, counter, reference, policy):
    contract = {**asdict(case), "variant": "RL_FULL"}
    key = hashlib.sha256(json.dumps(contract, sort_keys=True).encode()).hexdigest()[:16]
    path = reference/"points"/(key+".json")
    old = json.loads(path.read_text())
    if old["case"] != contract or old["status"] != "PASS" or old["count"]["status"] != "PASS":
        raise engine.PrefixFailure("REFERENCE_CONTRACT_MISMATCH")
    old_modes = load_modes(reference, old)
    comparisons = []
    for root, mode, prior, prior_mode in zip(roots, modes, old["roots"], old_modes):
        overlap = p.inner(model, mode, prior_mode)
        mac = overlap**2/(p.inner(model, mode, mode)*p.inner(model, prior_mode, prior_mode))
        sign = 1 if overlap >= 0 else -1
        errors = []
        for field in ("y", "dy", "ends"):
            new = np.vstack([a[field] for a in mode["arms"]])
            ref = np.vstack([a[field] for a in prior_mode["arms"]])
            for cols, scale in ((slice(0, 3), [1, 1, model.l]), (slice(3, 6), [1, 1, 1/model.l])):
                errors.append(float(np.max(abs((new[:, cols]*sign-ref[:, cols])*scale))/max(np.max(abs(ref[:, cols]*scale)), 1e-30)))
        comparisons.append(dict(relative_omega=abs(root["omega"]/prior["omega"]-1),
            relative_Lambda=abs(root["Lambda"]/prior["Lambda"]-1), mass_MAC=mac,
            max_scaled_field_error=max(errors)))
    checks = []
    for sample in old["count"]["samples"]:
        actual, count, _ = counter(sample["requested"])
        checks.append(dict(requested=sample["requested"], actual=actual, old_effective=sample["effective"],
                           count=count, old_count=sample["count"]))
    passed = len(roots) == len(old["roots"]) == 3 and all(
        max(c["relative_omega"], c["relative_Lambda"]) <= policy["reference_frequency_tol"] and
        abs(c["mass_MAC"]-1) <= policy["reference_shape_tol"] and
        c["max_scaled_field_error"] <= policy["reference_shape_tol"] for c in comparisons) and all(
        c["count"] == c["old_count"] and c["actual"] == c["old_effective"] for c in checks)
    return dict(status="PASS" if passed else "REGRESSION_FAILED", reference=path.relative_to(ROOT).as_posix(),
                reference_sha256=sha(path), modes=comparisons, count_checks=checks)


def straight_split(model, case, roots, modes, upper, policy):
    """Single CC rod of total length, independently separated axial/bending.

    Same bending functions, different boundary assembly; Ritz is the independent
    energy discretization. Transform copied from the accepted beta0 shape audit.
    """
    length = sum(case.lengths(model))
    cat = tim_catalog(model.carrier, length, upper, {"policy": dict(
        scan_intervals=policy["pole_scan_intervals"], root_xtol=policy["root_xtol"], root_rtol=policy["root_rtol"])})
    segment = p.axial.Segment(length, model.section.EA, model.section.rhoA, H=0., J=model.J)
    nmax = math.floor(cat["upper"]*math.sqrt(segment.m/(segment.EA-segment.J*cat["upper"]**2))*length/math.pi)
    axial = p.axial.explicit_frequencies(segment, nmax)*2*math.pi
    direct = sorted([(w, "axial", n+1) for n, w in enumerate(axial)] +
                    [(r["omega"], "bending", n+1) for n, r in enumerate(cat["roots"])])
    checks = []
    for root, mode, (w, family, index) in zip(roots, modes, direct):
        if family == "bending":
            b, _, col = p.balance(p.tim.finite_boundary_matrix(model.carrier, length, w, "timoshenko"))
            coef = col*np.linalg.svd(b)[2][-1]
        reference = dict(arms=[])
        for i, arm in enumerate(mode["arms"]):
            x = arm["x"] if i == 0 else length-arm["x"]
            y, dy = np.zeros((len(x), 6)), np.zeros((len(x), 6))
            if family == "bending":
                y[:, [1, 2, 4, 5]] = p.tim.finite_state_basis(model.carrier, length, w, x, "timoshenko")@coef
                dy[:, [1, 2, 4, 5]] = p.tim.finite_state_basis(model.carrier, length, w, x, "timoshenko", 1)@coef
            else:
                y[:, 0] = p.axial.basis(segment, w, x)[:, 1]
                dy[:, 0] = p.axial.basis(segment, w, x, 1)[:, 1]
                y[:, 3] = (segment.EA-segment.J*w*w)*dy[:, 0]
                dy[:, 3] = (segment.EA-segment.J*w*w)*p.axial.basis(segment, w, x, 2)[:, 1]
            if i == 1:
                transform = np.array([-1, -1, 1, 1, 1, -1])
                y *= transform; dy *= -transform
            reference["arms"].append(dict(y=y, dy=dy, weights=arm["weights"]))
        mac = p.inner(model, mode, reference)**2/(p.inner(model, mode, mode)*p.inner(model, reference, reference))
        checks.append(dict(family=family, family_index=index, omega=float(w),
                           relative_omega=abs(root["omega"]/w-1), mass_MAC=mac))
    ok = len(checks) == len(roots) and all(c["relative_omega"] <= policy["same_physics_frequency_tol"] and
                                         abs(c["mass_MAC"]-1) <= policy["reference_shape_tol"] for c in checks)
    return dict(status="PASS" if ok else "STRAIGHT_SPLIT_FAILED", modes=checks,
                bending_catalog=cat, axial_poles=axial.tolist(), lower_count=0,
                first_axial_omega=float(p.axial.explicit_frequencies(segment, 1)[0]*2*math.pi),
                count_below_catalog_ceiling=len(direct), reference="single CC rod; axial exact + saturated Tim catalog")


def solve_point(config, stage, beta, folder, source_id, neighbour=None):
    model = p.Model(**config["material"])
    case = p.Case(beta, config["mu"], config["k_joint"], config["k_end"])
    policy, scope = config["policy"], config["stages"][stage]
    contract = engine.point_contract(model, case, scope["K_plot"], scope["K_guard"], policy)
    key = engine.digest(dict(contract=contract, source_id=source_id, stage=stage))[:16]
    path = folder/"points"/(key+".json")
    if path.exists():
        saved = json.loads(path.read_text())
        verify_point_record(saved)
        if saved["contract"] != contract or saved["source_id"] != source_id:
            raise engine.PrefixFailure("POINT_CACHE_CONTRACT_MISMATCH")
        load_modes(folder, saved) if "shape_file" in saved else None
        return saved
    record = dict(stage=stage, beta_deg=beta, contract=contract, source_id=source_id, status="RUNNING")
    start = time.perf_counter()
    try:
        initial = []
        lower = float(model.omega(policy["Lambda_lower"]))
        upper = float(model.omega(policy["Lambda_upper"]))
        if stage == "B":
            values, _, diag = p.ritz(model, case, policy["ritz_orders"][-1], count=scope["K_guard"])
            initial = [dict(order=policy["ritz_orders"][-1], omega=values.tolist(), diagnostics=diag)]
            record["initial_ritz"] = initial
            upper = max(upper, float(values[-1])*policy["ritz_upper_padding"])
        attempts = []
        for attempt in range(policy["max_bound_retries"]+1):
            cats, catalog_attempts = catalogues(model, case, upper, policy)
            record["catalogs"] = cats
            counter = engine.PhysicalCounter(model, case, cats, policy)
            try:
                events, search = engine.counted_prefix(lambda w: p.balance(p.boundary(model, case, w))[0],
                    counter, lower, upper, scope["K_guard"], policy, neighbour)
                break
            except engine.PrefixFailure as error:
                attempts.append(dict(attempt=attempt, lower=lower, upper=upper,
                    catalog_attempts=catalog_attempts, code=error.code, details=error.details))
                record["bound_attempts"] = attempts
                if error.code == "INSUFFICIENT_RIGHT_COUNT":
                    upper *= policy["upper_growth"]
                elif error.code == "NONZERO_LOWER_COUNT":
                    lower /= 2
                else:
                    raise
        else:
            raise engine.PrefixFailure("BOUND_BUDGET_EXHAUSTED", attempts=attempts)
        record.update(search=search, catalogs=cats, catalog_attempts=catalog_attempts, bound_attempts=attempts)
        modes, roots = [], []
        for i, event in enumerate(events):
            local = engine.recovered_modes(model, case, event, policy["quadrature_order"], policy)
            for j, mode in enumerate(local):
                rank = event["first_rank"]+j
                roots.append(dict(omega=event["omega"], Omega=event["omega"]*model.time,
                    Lambda=float(model.Lambda(event["omega"])), sorted_position=rank, event_id=i+1,
                    multiplicity=event["multiplicity"], role="research" if rank <= scope["K_plot"] else (
                        "guard" if rank == scope["K_guard"] else "cluster_completion"),
                    coefficients=mode["coefficients"], diagnostics=mode["diagnostics"]))
                modes.append(mode)
        record["roots"] = roots
        shape = path.with_suffix(".npz")
        shape.parent.mkdir(parents=True, exist_ok=True)
        with shape.open("xb") as stream:
            np.savez_compressed(stream, **{f"mode{i+1}_arm{j+1}_{field}": arm[field]
                for i, mode in enumerate(modes) for j, arm in enumerate(mode["arms"])
                for field in ("x", "weights", "y", "dy", "ends")})
        record.update(shape_file=shape.relative_to(folder).as_posix(), shape_sha256=sha(shape),
                      quality=engine.mode_gates(model, modes, policy))
        record["quadrature"] = quadrature_check(model, case, roots, modes, policy)
        record["ritz"] = ritz_checks(model, case, roots, initial, policy)
        record["domain"] = dict(status="PASS", minimum_basis_limit=min(counter.rl_limit, counter.tim_limit),
            RL_limit=counter.rl_limit, Tim_limit=counter.tim_limit, catalog_limit=counter.domain,
            guard_omega=roots[scope["K_guard"]-1]["omega"], right_omega=search["guard_right_omega"],
            guard_to_catalog_margin=counter.domain-roots[scope["K_guard"]-1]["omega"],
            guard_to_Tim_ratio=roots[scope["K_guard"]-1]["omega"]/counter.tim_limit)
        record["close_gaps"] = [dict(ranks=[a["sorted_position"], b["sorted_position"]],
            gap_omega=b["omega"]-a["omega"], relative_gap=(b["omega"]-a["omega"])/b["omega"])
            for a, b in zip(roots[:-1], roots[1:])]
        checks = [record[k]["status"] for k in ("quality", "quadrature", "ritz", "domain")]
        if stage == "A":
            record["regression"] = regression(model, case, roots, modes, counter, ROOT/config["reference"], policy)
            checks.append(record["regression"]["status"])
        if stage == "B" and beta == 0:
            record["straight_split"] = straight_split(model, case, roots, modes, upper, policy)
            checks.append(record["straight_split"]["status"])
        record["status"] = "PASS" if all(c == "PASS" for c in checks) else "CHECK_FAILED"
        record["mode_identification"] = "SORTED_SUBSPACES" if any(e["multiplicity"] > 1 for e in events) else "SIMPLE_SORTED_MODES_NO_TRACKING"
    except (engine.PrefixFailure, ArithmeticError, ValueError, RuntimeError, OSError, np.linalg.LinAlgError) as error:
        record.update(status=getattr(error, "code", "NUMERICAL_FAILURE"),
                      failure=getattr(error, "details", {"reason": str(error), "type": type(error).__name__}))
    record["elapsed_seconds"] = time.perf_counter()-start
    record = jsonable(record)
    record["record_digest"] = engine.digest(record)
    write_new(path, record)
    return record


def verify_bundle(folder):
    manifest = json.loads((folder/"manifest.json").read_text())
    failures = [name for name, h in manifest["artifacts"].items() if sha(folder/name) != h]
    if failures:
        raise engine.PrefixFailure("ARTIFACT_HASH_MISMATCH", files=failures)
    points = [json.loads((folder/name).read_text()) for name in manifest["point_files"]]
    return manifest, points


def table(points):
    stream = io.StringIO()
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(["stage", "beta_deg", "rank", "role", "omega", "Omega", "Lambda", "multiplicity", "point_status"])
    for rec in points:
        for r in rec.get("roots", []):
            writer.writerow([rec["stage"], rec["beta_deg"], r["sorted_position"], r["role"],
                             r["omega"], r["Omega"], r["Lambda"], r["multiplicity"], rec["status"]])
    return stream.getvalue()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--stage", choices=["A", "B", "AB"])
    group.add_argument("--verify-cache", type=Path)
    group.add_argument("--tables-only", type=Path)
    args = parser.parse_args(argv)
    if args.verify_cache or args.tables_only:
        manifest, points = verify_bundle(args.verify_cache or args.tables_only)
        print(table(points) if args.tables_only else json.dumps(dict(status=manifest["status"], verified=len(manifest["artifacts"]))))
        return 0
    config = json.loads(CONFIG.read_text())
    validate_scope(config)
    reference = verify_reference(ROOT/config["reference"])
    _, provenance = identity(config, reference)
    provenance["execution_scope"] = "A" if args.stage == "A" else "AB"
    key = engine.digest(provenance)[:16]
    folder = OUTPUT/key
    folder.mkdir(parents=True, exist_ok=True)
    if (folder/"manifest.json").exists():
        manifest, points = verify_bundle(folder)
        print(json.dumps(dict(folder=str(folder), status=manifest["status"], mode="verified_cache")))
        return 0 if manifest["status"] == "PASS" else 2
    if not (folder/"identity.json").exists():
        write_new(folder/"identity.json", provenance)
    points, stages = [], {}
    # B always verifies A first, even when invoked directly.
    for stage in (["A"] if args.stage == "A" else ["A", "B"]):
        previous = None
        for beta in config["stages"][stage]["angles"]:
            rec = solve_point(config, stage, beta, folder, key, previous)
            points.append(rec)
            print(json.dumps(dict(stage=stage, beta=beta, status=rec["status"], seconds=rec["elapsed_seconds"])), flush=True)
            if rec["status"] == "PASS":
                previous = [r["omega"] for r in rec["roots"]]
            elif stage == "A":
                break
            else:
                previous = None
        stages[stage] = "PASS" if all(r["status"] == "PASS" for r in points if r["stage"] == stage) else "FAILED"
        if stages[stage] != "PASS" and stage == "A":
            break
    reference_after = verify_reference(ROOT/config["reference"])
    if reference_after != reference:
        raise engine.PrefixFailure("REFERENCE_CHANGED_DURING_RUN")
    csvpath = folder/"spectrum.csv"
    with csvpath.open("x", encoding="utf-8", newline="") as stream:
        stream.write(table(points))
    manifest = dict(status="PASS" if all(s == "PASS" for s in stages.values()) else "FAILED", stages=stages,
        command=sys.argv, identity=provenance, reference_after=reference_after,
        point_files=[x.relative_to(folder).as_posix() for x in sorted((folder/"points").glob("*.json"))],
        artifacts={x.relative_to(folder).as_posix(): sha(x) for x in sorted(folder.rglob("*")) if x.is_file()},
        excluded_stages=["C", "D", "E"])
    write_new(folder/"manifest.json", manifest)
    print(json.dumps(dict(folder=str(folder), stages=stages, status=manifest["status"])), flush=True)
    return 0 if manifest["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
