"""Configurable, sorted RL_FULL/RR Lambda maps; no physical equations here.

Separate from immutable A-B/C0-C1 certification CLIs: fast_plot, resumable
point transactions and read-only CSV rendering are a different I/O contract.
The scientific engine is the unchanged individual-section counted prefix.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time
import uuid

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
OUTPUT = ROOT / "results/rl_full_rr_maps"
CONFIG = ROOT / "data/input/rl_full_rr_maps.json"
AXES = {"BETA": "beta_deg", "MU": "mu", "TAU": "tau"}
# Reviewed interface-only migration; a future unreviewed wrapper edit cannot
# silently reuse these older source snapshots. Physical dependencies stay exact.
PREFIX_INTERFACE_AST_SHA256 = "7687bfb9dc930cfa1ffbe9851fe07ba44315ad4dd2851d0cfda3df0c2732fc5c"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def publish(path, obj):
    """Atomic, non-overwriting point/manifest transaction (also on Windows)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + "." + uuid.uuid4().hex + ".tmp")
    with temporary.open("x", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False, allow_nan=False, default=lambda value:value.item())
        f.write("\n")
    if path.exists():
        raise FileExistsError(path)
    temporary.rename(path)


def configuration(path):
    c = read(path)
    reference = read(ROOT / "data/input/rl_full_rr_geometry_c.json")
    fixed = dict(c["material"], aspect=4., epsilon=.005)
    if (fixed != reference["material"] or not math.isfinite(c["material"]["aspect"])
        or c["material"]["aspect"] <= 0 or c["material"]["epsilon"] not in (.005, .01)
        or c["numerical_policy"] != "data/input/rl_full_rr_geometry_c.json"):
        raise ValueError("Positive aspect, epsilon=.005/.01 and the accepted material/numerical policy required")
    prefix_ranks(c)
    names = set()
    for spec in c["maps"]:
        if spec["name"] in names or not spec["name"].replace("_", "").isalnum():
            raise ValueError("Unique alphanumeric/underscore map names required")
        names.add(spec["name"])
        axis = AXES[spec["mode"]]
        if set(spec["fixed"]) != set(AXES.values()) - {axis}:
            raise ValueError("Exactly two fixed geometry parameters required")
        grid = spec["grid"]
        if not grid or not all(math.isfinite(x) for x in grid) or any(a >= b for a, b in zip(grid, grid[1:])):
            raise ValueError("A finite, strictly increasing grid is required")
        for value in grid:
            q = dict(spec["fixed"], **{axis: value})
            if not (0 <= q["beta_deg"] <= 180 and abs(q["mu"]) < 1 and abs(q["tau"]) < 1):
                raise ValueError("Require beta in [0,180], |mu|<1 and |tau|<1")
    if not 0 < c["initial_upper_omega"] < 10:
        raise ValueError("Finite positive initial search bound below 10 required")
    return c, reference["policy"]


def prefix_ranks(config):
    pair = config.get("K_plot", 8), config.get("K_guard", 9)
    if any(type(k) is not int for k in pair) or pair not in ((6, 7), (8, 9)):
        raise ValueError("Supported prefixes: K6/guard7 and K8/guard9")
    return pair


def local_policy(spec, k_plot=8, k_guard=9):
    prefix_ranks(dict(K_plot=k_plot, K_guard=k_guard))
    return dict(frequency_map_policy="frequency-map-v1", calculation_mode="fast_plot",
        spectrum_semantics="sorted_positions", sweep_parameter=AXES[spec["mode"]],
        parameter_grid=spec["grid"], K_plot=k_plot, K_guard=k_guard, guard_root_role="completeness_only",
        neighbour_audit=f"max |omega_next/omega_previous-1| over ranks 1..{k_plot} > 0.2",
        local_repair_policy="one count-only rerun of affected point without predictor; no added parameter nodes",
        strict_audit_default=False,
        shape_policy="recover full guard prefix solely for required BC/ODE/energy/mass gates; no MAC",
        unresolved_policy="retain failure and gap; no automatic strict audit")


def physical_key(contract):
    """Comparable A-B and individual-section contracts, excluding audit extras."""
    model = dict(contract["model"])
    mu = contract["case"]["mu"]
    tau = contract["tau"]
    model.update(mu=mu, tau=tau)
    case = dict(contract["case"], tau=tau)
    def canonical(value):
        if isinstance(value, dict): return {k:canonical(v) for k,v in value.items()}
        if isinstance(value, list): return [canonical(v) for v in value]
        if isinstance(value, (int, float)) and not isinstance(value, bool): return float(value)
        return value
    return digest(canonical(dict(model=model, case=case, K_plot=contract["K_plot"],
        K_guard=contract["K_guard"], policy=contract["policy"],
        time_reference=contract["time_reference"], state_order=contract["state_order"],
        coordinates=contract["coordinates"], semantics=contract["spectrum_semantics"])))


def reusable_prefix(source, target):
    """A certified higher prefix can supply a lower one, never conversely."""
    sk, sg = prefix_ranks(source)
    tk, tg = prefix_ranks(target)
    if sk < tk or sg < tg or physical_key(dict(source, K_plot=tk, K_guard=tg)) != physical_key(target):
        return False
    normalizations = {"Omega=omega*t_ref=Lambda**2", "Omega=omega*t_ref=Lambda**2; common reference section"}
    if source.get("normalization") not in normalizations or target.get("normalization") not in normalizations:
        return False
    # Same model parameters imply the same canonical coefficients. Also reject
    # internally inconsistent stored derived coefficients, not just parameter tags.
    from scripts.lib import rl_full_rr_sections_prefix as e
    from scripts.lib.rl_full_rr_aspect import make_model
    for contract in (source, target):
        model = make_model(**dict(contract["model"], mu=contract["case"]["mu"], tau=contract["tau"]))
        case = e.p.Case(**dict(contract["case"], tau=contract["tau"]))
        expected = e.point_contract(model, case, contract["K_plot"], contract["K_guard"], contract["policy"])
        for field in ("arms", "reference_section"):
            if field in contract and contract[field] != expected[field]:
                return False
        if "lengths" in contract and contract["lengths"] != list(case.lengths(model)):
            return False
        if "section" in contract and contract["section"] != {k:expected["arms"][0][k] for k in contract["section"]}:
            return False
    return True


def find_reusable(contract, known):
    for kp, kg in (prefix_ranks(contract), (8, 9)):
        found = known.get(physical_key(dict(contract, K_plot=kp, K_guard=kg)))
        if found and reusable_prefix(found[1]["contract"], contract):
            return found
    return None


def prefix_view(source_path, target, folder, source_id):
    """Immutable derived record/NPZ, certified by a saved larger prefix.

    No matrix, basis, pole catalogue, count or eigenvalue call. Complete a
    cluster through the guard using the already certified source interval.
    """
    import copy
    import numpy as np
    from scripts.lib import rl_full_rr_sections_prefix as e
    from scripts.lib.rl_full_rr_aspect import make_model
    from scripts.analysis import verify_rl_full_rr_prefix as ab
    source = verify_record(source_path)
    if source["status"] != "PASS" or not reusable_prefix(source["contract"], target):
        raise ValueError("Source does not certify the requested physical prefix")
    kp, kg = prefix_ranks(target)
    search = source["search"]
    if search["status"] != "PASS" or search["lower_count"] != 0 or search["unresolved_intervals"]:
        raise ValueError("Source completeness is unresolved")
    if len(source["roots"]) < source["contract"]["K_guard"] or source["quality"]["status"] != "PASS":
        raise ValueError("Incomplete certified source")
    retained = kg
    while True:
        last = source["roots"][retained-1]
        right = last["omega"] + target["policy"]["cluster_probe_relative"]*max(1., abs(last["omega"]))
        if retained < len(source["roots"]) and source["roots"][retained]["omega"] <= right:
            retained += 1
        else:
            break
    if retained-kg > target["policy"]["max_cluster_extra_ranks"] or right > search["guard_right_omega"]:
        raise ValueError("Source does not cover the complete requested guard cluster")
    if search["guard_right_count"] != len(source["roots"]):
        raise ValueError("Source right count and saved positions disagree")
    proof = dict(source_record=source_path.relative_to(ROOT).as_posix(), source_record_sha256=sha(source_path),
                 source_npz=(source_path.parent.parent/source["shape_file"]).relative_to(ROOT).as_posix(),
                 source_npz_sha256=source["shape_sha256"], source_K_plot=source["contract"]["K_plot"],
                 source_K_guard=source["contract"]["K_guard"], selected_positions=list(range(1,retained+1)),
                 source_search=search, spectral_calls=0)
    key = digest(dict(contract=target, source_id=source_id, reuse=proof))[:16]
    path = folder/"points"/(key+".json")
    if path.exists():
        return path, verify_record(path)
    modes = ab.load_modes(source_path.parent.parent,source)[:retained]
    roots = copy.deepcopy(source["roots"][:retained])
    for rank, (root, mode) in enumerate(zip(roots,modes), 1):
        root["role"] = "plot" if rank<=kp else "guard" if rank==kg else "cluster_completion"
        mode["diagnostics"] = root["diagnostics"]
    model = make_model(**target["model"])
    quality = e.mode_gates(model,modes,target["policy"])
    if quality["status"] != "PASS":
        raise ValueError("Saved prefix failed physical/mass gates")
    shape = path.with_suffix(".npz"); shape.parent.mkdir(parents=True,exist_ok=True)
    with shape.open("xb") as f:
        np.savez_compressed(f, **{f"mode{i+1}_arm{j+1}_{field}":arm[field]
            for i, mode in enumerate(modes) for j, arm in enumerate(mode["arms"])
            for field in ("x","weights","y","dy","ends")})
    record = dict(contract=target, source_id=source_id, calculation_mode="reuse_certified_prefix",
        status="PASS", beta_deg=source["beta_deg"], mu=source["mu"], tau=source["tau"],
        geometry=source.get("geometry"), roots=roots, quality=quality, reuse=proof,
        search=dict(status="PASS", lower_count=0, requested_guard=kg, computed_ranks=retained,
            guard_right_omega=right, guard_right_count=retained, unresolved_intervals=[],
            proof_type="subprefix_of_certified_counted_interval", new_count_queries=0),
        domain=dict(source.get("domain",{}),guard_omega=roots[kg-1]["omega"],right_omega=right),
        shape_file=shape.relative_to(folder).as_posix(), shape_sha256=sha(shape), elapsed_seconds=0.)
    record=ab.jsonable(record);record["record_digest"]=digest(record);publish(path,record)
    return path,record


def saved_control_reference(c, policy):
    """Explicit archived square controls; never run new Ritz as a side effect."""
    folder = ROOT/c["control_reference"]
    if sha(folder/"manifest.json") != c["control_reference_sha256"]:
        raise ValueError("Changed control reference")
    saved = verify_manifest(folder)
    material = saved["identity"]["material"]
    if (saved["status"] != "PASS" or material["aspect"] != 1 or c["material"]["aspect"] != 1
            or dict(material,epsilon=c["material"]["epsilon"]) != c["material"]
            or saved["identity"]["policy"] != policy):
        raise ValueError("Incompatible archived square controls")
    for file,h in saved["identity"]["sources"].items():
        if file in ("scripts/analysis/map_rl_full_rr.py", "scripts/lib/rl_full_rr_square_checks.py"):
            continue
        if sha(ROOT/file)!=h:raise ValueError("Changed physical control dependency: "+file)
    return dict(saved,folder=folder.relative_to(ROOT).as_posix(),manifest_sha256=sha(folder/"manifest.json"),
                coverage_epsilon=material["epsilon"],new_ritz_calls=0)


def verify_record(path):
    rec = read(path)
    saved = rec.pop("record_digest")
    if digest(rec) != saved:
        raise ValueError(f"Record digest mismatch: {path}")
    rec["record_digest"] = saved
    if "shape_file" in rec and sha(path.parent.parent / rec["shape_file"]) != rec["shape_sha256"]:
        raise ValueError(f"Shape checksum mismatch: {path}")
    return rec


def verify_manifest(folder):
    m = read(folder / "manifest.json")
    for name, h in m["artifacts"].items():
        if sha(folder / name) != h:
            raise ValueError(f"Artifact checksum mismatch: {folder/name}")
    return m


def accepted_sources(c):
    """Only named accepted RR bundles; never the old spring cache."""
    known, manifests, source_hashes = {}, {}, {}
    for name in c["references"]:
        folder = ROOT / name
        m = verify_manifest(folder)
        if m["status"] != "PASS":
            raise ValueError("Reference bundle did not pass")
        for file, h in m["identity"]["sources"].items():
            if sha(ROOT/file) != h:
                raise ValueError(f"Historical source changed: {file}")
        source_hashes.update(m["identity"]["sources"])
        manifests[name] = sha(folder/"manifest.json")
        paths = ([folder/n for n in m["point_files"]] if "point_files" in m else
                 [ROOT/x["record"] for x in m["C0_points"] + m["C1_points"]])
        for path in paths:
            rec = verify_record(path)
            if rec["status"] == "PASS" and rec["contract"]["K_guard"] == 9:
                known.setdefault(physical_key(rec["contract"]), (path, rec, "accepted_reference"))
    return known, manifests, source_hashes


def geometry(model, case):
    arms = []
    for a, length in zip(model.arms, case.lengths(model)):
        s = a.section
        h = math.sqrt(12)*a.l*a.epsilon
        b = a.aspect*h
        arms.append(dict(L=length, h=h, b=b, aspect=a.aspect, A=s.rhoA/a.rho, I=s.EI/a.E,
            I_P=b*h*(b*b+h*h)/12, EA=s.EA, EI=s.EI, KGA=s.KGA,
            rhoA=s.rhoA, rhoI=s.rhoI, J=a.J, L_over_h=length/h, L_over_b=length/b,
            slenderness_geometric_pass=length/h>=10))
    error = abs(sum(a["A"]*a["L"] for a in arms)/(2*model.section.rhoA/model.rho*model.l)-1)
    if error > 1e-12:
        raise ValueError("Volume conservation failed")
    from scripts.lib.frequency_geometry_style import classify
    return dict(arms=arms, volume_relative_error=error, reference_time=model.time,
                **classify(*(a["L_over_h"] for a in arms)),
                slenderness_geometric_pass=all(a["slenderness_geometric_pass"] for a in arms),
                slenderness_rule="L_i/h_i >= 10; preliminary geometry only, not a frequency cutoff")


def solve(contract, folder, source_id, initial_upper, neighbour=None, repair=False):
    # Lazy import: plot-only/tables-only cannot even assemble the physical engine.
    import numpy as np
    from scripts.lib import rl_full_rr_sections_prefix as e
    from scripts.lib.rl_full_rr_aspect import make_model
    from scripts.analysis import verify_rl_full_rr_prefix as ab
    p = e.p
    policy = contract["policy"]
    k_plot, k_guard = prefix_ranks(contract)
    model, case = make_model(**contract["model"]), p.Case(**contract["case"])
    key = digest(dict(contract=contract, source_id=source_id, repair=repair))[:16]
    path = folder/"points"/(key+".json")
    if path.exists():
        return path, verify_record(path), False
    start = time.perf_counter()
    rec = dict(contract=contract, source_id=source_id, calculation_mode="fast_plot",
        beta_deg=case.beta, mu=case.mu, tau=case.tau, status="RUNNING",
        repair=repair, predictor_omega=neighbour, independent_ritz="accepted reference controls only")
    try:
        rec["geometry"] = geometry(model, case)
        lo = float(model.omega(policy["Lambda_lower"]))
        hi = initial_upper if neighbour is None else max(float(neighbour[k_guard-1])*1.1, float(model.omega(policy["Lambda_upper"])))
        if repair:
            hi *= 1.125  # different partition and upper bound, unchanged physics/tolerances
        rec["bound_attempts"] = []
        for attempt in range(policy["max_bound_retries"]+1):
            cats, trials = e.catalogues(model, case, hi, policy)
            rec.update(catalogs=cats, catalog_attempts=trials)
            counter = e.PhysicalCounter(model, case, cats, policy)
            try:
                events, search = e.counted_prefix(lambda w:p.balance(p.boundary(model, case, w))[0],
                    counter, lo, hi, k_guard, policy, neighbour)
                break
            except e.PrefixFailure as ex:
                rec["bound_attempts"].append(dict(lower=lo, upper=hi, code=ex.code, details=ex.details))
                if ex.code == "INSUFFICIENT_RIGHT_COUNT": hi *= policy["upper_growth"]
                elif ex.code == "NONZERO_LOWER_COUNT": lo /= 2
                else: raise
        else:
            raise e.PrefixFailure("BOUND_BUDGET_EXHAUSTED")
        rec["search"] = search
        modes, roots = [], []
        for event_id, event in enumerate(events, 1):
            for j, mode in enumerate(e.recovered_modes(model, case, event, policy["quadrature_order"], policy)):
                rank = event["first_rank"]+j
                roots.append(dict(omega=event["omega"], Omega=event["omega"]*model.time,
                    Lambda=float(model.Lambda(event["omega"])), sorted_position=rank,
                    multiplicity=event["multiplicity"], event_id=event_id,
                    role="plot" if rank<=k_plot else "guard" if rank==k_guard else "cluster_completion",
                    coefficients=mode["coefficients"], diagnostics=mode["diagnostics"]))
                modes.append(mode)
        rec.update(roots=roots, quality=e.mode_gates(model, modes, policy),
            domain=dict(arms=counter.limits, catalog_limit=counter.domain,
                        guard_omega=roots[k_guard-1]["omega"], right_omega=search["guard_right_omega"]))
        # Unique shape transaction also permits resuming after an interrupted JSON publication.
        shape = path.with_name(key+"_"+uuid.uuid4().hex[:8]+".npz")
        shape.parent.mkdir(parents=True, exist_ok=True)
        tmp = shape.with_suffix(".tmp")
        with tmp.open("xb") as f:
            np.savez_compressed(f, **{f"mode{i+1}_arm{j+1}_{field}":arm[field]
                for i, mode in enumerate(modes) for j, arm in enumerate(mode["arms"])
                for field in ("x", "weights", "y", "dy", "ends")})
        tmp.rename(shape)
        rec.update(shape_file=shape.relative_to(folder).as_posix(), shape_sha256=sha(shape), status=rec["quality"]["status"])
    except (ArithmeticError, ValueError, RuntimeError, np.linalg.LinAlgError) as ex:
        rec.update(status=getattr(ex, "code", "NUMERICAL_FAILURE"),
                   failure=getattr(ex, "details", dict(reason=str(ex), type=type(ex).__name__)))
    rec["elapsed_seconds"] = time.perf_counter()-start
    rec = ab.jsonable(rec)
    rec["record_digest"] = digest(rec)
    publish(path, rec)
    return path, rec, True


def jump(a, b):
    k = min(a.get("contract", {}).get("K_plot", 8), b.get("contract", {}).get("K_plot", 8))
    return max(abs(y["omega"]/x["omega"]-1) for x,y in zip(a["roots"][:k], b["roots"][:k]))


def agreement(a, b, tolerance):
    ra, rb = a["roots"], b["roots"]
    return (len(ra) == len(rb) and all(x["multiplicity"] == y["multiplicity"] and
        abs(x["omega"]/y["omega"]-1) <= tolerance for x,y in zip(ra, rb)))


def solver_sources_compatible(previous, current, legacy_a4=False):
    """Explicit source snapshots permit metadata/control-only migrations.

    Gate/report code is not a spectral dependency; all physical kernels and
    aspect coefficients must match. Unknown source revisions are rejected.
    """
    import ast
    name = Path(__file__).relative_to(ROOT).as_posix()
    helper = "scripts/lib/rl_full_rr_square_checks.py"
    previous = {k:v for k,v in previous.items() if k!=helper}
    current = {k:v for k,v in current.items() if k!=helper}
    if previous == current:return True
    snapshots = {
        "4e564b482def6708b26a6c16c3e46469aa3c62d829ca8663cd1e93a60fc18381":"private/dissertation/_history/20261108_slenderness_rule/scripts/analysis/map_rl_full_rr.py",
        "240b1591abc2a948691ab4f1cedb362d97731a19c05c354d9c0c420d386ae120":"private/dissertation/_history/20261109_square_pilot/scripts/analysis/map_rl_full_rr.py",
        "9dbce4239c0748f9a0956157744779c8ddbf1c27c91ce0ee4c0cf1ea184f0605":"results/rl_full_rr_maps/_audit/map_rl_full_rr_v1.py",
        "edc4f24f75f3e7c670c98ff83af7cda386fa099aa0416611d6648487a8de5587":"results/rl_full_rr_maps/_audit/map_rl_full_rr_v2.py",
        "61371b890eba4c6f549f50d52190afcd9350c31ed8ed74b81255df03fbe4d794":"results/rl_full_rr_square/_audit/map_before_json_fix.py"}
    legacy_sha = previous.get(name)
    if legacy_sha not in snapshots: return False
    legacy = ROOT/snapshots[legacy_sha]
    if not legacy.exists() or sha(legacy) != legacy_sha:
        return False
    sources = dict(previous); sources[name] = current[name]
    extension = "scripts/lib/rl_full_rr_aspect.py"
    if legacy_a4 and extension in current: sources.setdefault(extension,current[extension])
    if sources != current:return False
    def functions(path):
        return {n.name:ast.dump(n, include_attributes=False) for n in ast.parse(path.read_text(encoding="utf-8")).body if isinstance(n, ast.FunctionDef)}
    before, after = functions(legacy), functions(Path(__file__))
    reviewed = ("configuration", "prefix_ranks", "physical_key", "reusable_prefix", "find_reusable",
                "prefix_view", "saved_control_reference", "local_policy", "digest", "verify_record",
                "accepted_sources", "solve", "geometry", "jump", "agreement", "write_tables", "run", "cache_compatible")
    if (current.get(name) == sha(__file__) and
            digest({k:after[k] for k in reviewed}) == PREFIX_INTERFACE_AST_SHA256):
        return True
    functions_to_check = ["physical_key", "local_policy", "digest", "verify_record"]
    if not legacy_a4 or legacy_sha.startswith("61371"):
        functions_to_check += ["solve", "geometry"]
    return all(before[k] == after[k] for k in functions_to_check)


def cache_compatible(old, current):
    if any(old.get(k) != current[k] for k in ("versions", "policy", "initial_upper_omega", "material")):return False
    if old.get("engine_id") == current["engine_id"]:return True
    return solver_sources_compatible(old.get("sources",{}), current["sources"], legacy_a4=current["material"]["aspect"]==4)


def write_tables(folder, entries):
    fields = ["map", "mode", "x", "beta_deg", "mu", "tau", "rank", "Lambda", "omega", "Omega",
              "multiplicity", "status", "origin", "record", "record_sha256", "aspect",
              "L1_over_h1", "L2_over_h2", "L1_over_b1", "L2_over_b2", "slenderness_geometric_pass",
              "K_plot", "K_guard", "role", "epsilon", "chi", "slenderness_criterion_pass", "slenderness_threshold", "line_style"]
    with (folder/"spectrum.csv").open("x", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for item in entries:
            rec = verify_record(ROOT/item["record"])
            kp, kg = prefix_ranks(rec["contract"])
            for rank in range(1, max(kg, len(rec.get("roots", [])))+1):
                row = {k:item[k] for k in fields if k in item}
                row["rank"] = rank
                row.update(K_plot=kp, K_guard=kg, epsilon=rec["contract"]["model"]["epsilon"],
                           role="plot" if rank<=kp else "guard" if rank==kg else "cluster_completion")
                geo = item["geometry"]
                row.update(aspect=rec["contract"]["model"]["aspect"],
                    L1_over_h1=geo["arms"][0]["L_over_h"], L2_over_h2=geo["arms"][1]["L_over_h"],
                    L1_over_b1=geo["arms"][0]["L_over_b"], L2_over_b2=geo["arms"][1]["L_over_b"],
                    slenderness_geometric_pass=geo["slenderness_geometric_pass"])
                from scripts.lib.frequency_geometry_style import classify
                row.update({k:v for k,v in classify(row["L1_over_h1"],row["L2_over_h2"]).items() if k in fields})
                if item["status"] == "PASS":
                    root = rec["roots"][rank-1]
                    row.update({k:root[k] for k in ("Lambda", "omega", "Omega", "multiplicity")})
                writer.writerow(row)


def run(c, policy, specs):
    import numpy as np
    import scipy
    from scripts.lib import rl_full_rr_sections_prefix as e
    from scripts.lib.rl_full_rr_aspect import make_model
    k_plot, k_guard = prefix_ranks(c)
    known, refs, sources = accepted_sources(c)
    sources[Path(__file__).relative_to(ROOT).as_posix()] = sha(__file__)
    sources["scripts/lib/rl_full_rr_aspect.py"] = sha(ROOT/"scripts/lib/rl_full_rr_aspect.py")
    gate = None
    if c["material"]["aspect"] != 4:
        from scripts.lib.rl_full_rr_square_checks import control_gate
        gate = saved_control_reference(c, policy) if "control_reference" in c else control_gate(c, policy)
        if gate["status"] != "PASS":
            raise RuntimeError("Aspect control gate failed; maps not started: "+gate["folder"])
        sources["scripts/lib/rl_full_rr_square_checks.py"] = sha(ROOT/"scripts/lib/rl_full_rr_square_checks.py")
    versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
    engine_id = digest(dict(sources=sources, versions=versions, initial_upper=c["initial_upper_omega"],
        policy=policy, mode="fast_plot", repair_policy=local_policy(specs[0],k_plot,k_guard)["local_repair_policy"]))
    identity = dict(schema="rl-full-rr-map-prefix-v3", material=c["material"], maps=specs,
        K_plot=k_plot, K_guard=k_guard,
        local_policies={s["name"]:local_policy(s,k_plot,k_guard) for s in specs}, policy=policy,
        initial_upper_omega=c["initial_upper_omega"], sources=sources, references=refs,
        engine_id=engine_id, versions=versions,
        head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip())
    if gate:
        identity["aspect_gate"] = {k:gate[k] for k in ("folder", "manifest_sha256")}
        for item in gate["points"]:
            rec = verify_record(ROOT/item["record"])
            known.setdefault(physical_key(rec["contract"]), (ROOT/item["record"], rec, "accepted_aspect_gate"))
    folder = OUTPUT/digest(identity)[:16]
    if (folder/"manifest.json").exists():
        m = verify_manifest(folder)
        print(json.dumps(dict(bundle=str(folder), status=m["status"], root_calculations=0, reused_bundle=True)), flush=True)
        return folder
    # Reuse completed fast_plot bundles only under identical engine/numerical identity.
    known_repairs = {}
    reused_bundles = []
    for old in sorted(OUTPUT.glob("*/manifest.json")):
        previous = read(old)
        if not cache_compatible(previous.get("identity", {}), identity): continue
        previous = verify_manifest(old.parent)
        reused_bundles.append(dict(manifest=old.relative_to(ROOT).as_posix(), sha256=sha(old),
            compatibility="same_engine" if previous["identity"]["engine_id"]==engine_id else "verified_snapshot_same_physics"))
        for item in previous["points"]:
            if item["status"] != "PASS": continue
            path = ROOT/item["record"]
            if sha(path) != item["record_sha256"]: raise ValueError("Referenced point changed")
            rec = verify_record(path)
            known.setdefault(physical_key(rec["contract"]), (path, rec, "previous_fast_plot"))
        for item in previous["repairs"]:
            if item["status"] != "REPRODUCED": continue
            path = ROOT/item["rerun"]
            if sha(path) != item["rerun_sha256"]: raise ValueError("Referenced repair changed")
            rec = verify_record(path)
            if rec["status"] == "PASS": known_repairs.setdefault(physical_key(rec["contract"]), (path, rec))
    folder.mkdir(parents=True, exist_ok=True)
    if not (folder/"identity.json").exists(): publish(folder/"identity.json", identity)
    start = time.perf_counter(); entries = []; repairs = []; new_count = 0; checkpoints = 0
    for spec in specs:
        previous = None
        axis = AXES[spec["mode"]]
        for value in spec["grid"]:
            q = dict(spec["fixed"], **{axis:value})
            model = make_model(**c["material"], mu=q["mu"], tau=q["tau"])
            case = e.p.Case(q["beta_deg"], q["mu"], "RIGID", "RIGID", q["tau"])
            contract = e.point_contract(model, case, k_plot, k_guard, policy)
            geometry_check = geometry(model, case)
            key = physical_key(contract)
            reusable = find_reusable(contract, known)
            if reusable:
                path, rec, origin = reusable
                if prefix_ranks(rec["contract"]) != (k_plot,k_guard):
                    path, rec = prefix_view(path,contract,folder,engine_id)
                    origin = "higher_certified_prefix"
                    known[key]=(path,rec,origin)
            else:
                seed = [r["omega"] for r in previous["roots"]] if previous else None
                path, rec, calculated = solve(contract, folder, engine_id, c["initial_upper_omega"], seed)
                new_count += int(calculated); checkpoints += int(not calculated)
                origin = "new_fast_plot" if calculated else "checkpoint"
                if rec["status"] == "PASS": known[key] = (path, rec, "same_run")
            item = dict(map=spec["name"], mode=spec["mode"], x=value, **q, status=rec["status"],
                origin=origin, record=path.relative_to(ROOT).as_posix(), record_sha256=sha(path), geometry=geometry_check)
            if "reuse" in rec:
                item["reuse"] = {k:v for k,v in rec["reuse"].items() if k!="source_search"}
            if rec["status"] == "PASS":
                if previous is not None:
                    item["neighbour_relative_jump"] = jump(previous, rec)
                previous = rec
            else:
                previous = None
            entries.append(item)
            print(json.dumps({k:item[k] for k in ("map", "x", "status", "origin")}), flush=True)
    # Audit only actual adjacent requested nodes; no bridging a failed point.
    for left, right in zip(entries, entries[1:]):
        if left["map"] != right["map"] or left["status"] != "PASS" or right["status"] != "PASS": continue
        a, b = verify_record(ROOT/left["record"]), verify_record(ROOT/right["record"])
        delta = jump(a, b)
        if delta <= policy["neighbour_relative_jump"]: continue
        key = physical_key(b["contract"])
        if key in known_repairs:
            rp, rr = known_repairs[key]; calculated = False
            repair_origin = "reused_verified_repair"
        else:
            rp, rr, calculated = solve(b["contract"], folder, engine_id, c["initial_upper_omega"], repair=True)
            repair_origin = "new_repair" if calculated else "checkpoint"
        new_count += int(calculated); checkpoints += int(repair_origin == "checkpoint")
        ok = rr["status"] == "PASS" and agreement(b, rr, policy["same_physics_frequency_tol"])
        repairs.append(dict(map=right["map"], x=right["x"], trigger="NEIGHBOUR_JUMP", jump=delta,
            original=right["record"], rerun=rp.relative_to(ROOT).as_posix(), rerun_sha256=sha(rp),
            origin=repair_origin, status="REPRODUCED" if ok else "UNRESOLVED"))
        if not ok: right["status"] = "NEIGHBOUR_UNRESOLVED"
    write_tables(folder, entries)
    artifacts = {p.relative_to(folder).as_posix():sha(p) for p in sorted(folder.rglob("*")) if p.is_file() and p.suffix != ".tmp"}
    m = dict(status="PASS" if all(p["status"] == "PASS" for p in entries) else "PARTIAL",
        identity=identity, points=entries, repairs=repairs, artifacts=artifacts, reused_bundles=reused_bundles,
        statistics=dict(grid_points=len(entries), new_grid_points=sum(p["origin"]=="new_fast_plot" for p in entries),
            reused_grid_points=sum(p["origin"]!="new_fast_plot" for p in entries),
            root_calculations=new_count, resumed_transactions=checkpoints, local_repairs=len(repairs),
            reused_repairs=sum(r["origin"]=="reused_verified_repair" for r in repairs),
            independent_ritz_calls=0, spectrum_seconds=time.perf_counter()-start))
    publish(folder/"manifest.json", m)
    print(json.dumps(dict(bundle=str(folder), status=m["status"], statistics=m["statistics"])), flush=True)
    return folder


def saved_table(folder):
    """Presentation uses only manifest + CSV, never point arrays or solver imports."""
    m = read(folder/"manifest.json")
    if sha(folder/"spectrum.csv") != m["artifacts"]["spectrum.csv"]:
        raise ValueError("Spectrum table checksum mismatch")
    with (folder/"spectrum.csv").open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    return m, rows


def render(folder, style, tables_only=False):
    from scripts.lib import frequency_geometry_style as line_geometry
    m, rows = saved_table(folder)
    if tables_only:
        result = dict(bundle=str(folder), csv=str(folder/"spectrum.csv"), rows=len(rows),
                      root_calculations=0, calculation_mode="tables_only")
        print(json.dumps(result)); return result
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    start = time.perf_counter()
    spec = dict(calculation_mode="plot_only", source_manifest=sha(folder/"manifest.json"),
        source_csv=sha(folder/"spectrum.csv"), style=style, renderer=sha(__file__), matplotlib=matplotlib.__version__,
        geometry_sources=line_geometry.provenance())
    out = OUTPUT/"figures"/(folder.name+"_"+digest(spec)[:12])
    if (out/"manifest.json").exists():
        verify_manifest(out); print(json.dumps(dict(figures=str(out), root_calculations=0, reused=True))); return out
    out.mkdir(parents=True, exist_ok=True)
    material = m["identity"].get("material", dict(epsilon=.005, aspect=4.))
    k_plot, _ = prefix_ranks(m["identity"])
    geometries = {entry["name"]: line_geometry.FrequencyGeometry.from_spec(entry, material,
                  entry.get("grid") or [float(r["x"]) for r in rows if r["map"] == entry["name"] and int(r["rank"]) == 1])
                  for entry in m["identity"]["maps"]}
    geometry_diagnostics = line_geometry.write_geometry_diagnostics(out, geometries)
    captions = {}
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#332288", "#AA4499", "#008899", "#886600"]
    markers = ["o", "s", "^", "D", "v", "P", "X", "h"]
    layouts = ["all", "split"] if style.get("layout", "both") == "both" else [style["layout"]]
    if any(v not in ("all", "split") for v in layouts): raise ValueError("layout: all, split or both")
    gallery = ["# RL_FULL / Timoshenko: Lambda maps", "", "Sorted positions; markers are computed nodes. Lines are visual guides, not tracked modal branches.", ""]
    with plt.rc_context({"font.size":style.get("font_size",11), "pdf.fonttype":42}):
        for entry in m["identity"]["maps"]:
            subset = [r for r in rows if r["map"]==entry["name"]]
            fixed = ", ".join((r"$\beta="+f"{v:g}"+r"^\circ$" if k=="beta_deg" else "$\\"+k+f"={v:g}$") for k,v in entry["fixed"].items())
            aspect = m["identity"].get("material", {}).get("aspect", 4)
            thickness = "" if entry["mode"] in ("EPS", "EPSILON") else fr", $\varepsilon={material.get('epsilon', .005):g}$"
            title = style.get("title") or (fixed + thickness + fr", $a={aspect:g}$; RL_FULL / Timoshenko, RR")
            for layout in layouts:
                half = k_plot//2
                groups = [range(1,k_plot+1)] if layout=="all" else [range(1,half+1), range(half+1,k_plot+1)]
                fig, axes = plt.subplots(len(groups), 1, figsize=style.get("figsize", [10,6] if layout=="all" else [10,8]), squeeze=False, sharex=True)
                for ax, ranks in zip(axes[:,0], groups):
                    for rank in ranks:
                        data = [r for r in subset if int(r["rank"])==rank]
                        if [float(r["x"]) for r in data] != list(geometries[entry["name"]].x):
                            raise ValueError("Saved frequency grid differs from the declared geometry grid")
                        line_geometry.plot_frequency_curve(ax, geometries[entry["name"]],
                            [float(r["Lambda"]) if r["status"]=="PASS" else float("nan") for r in data],
                            color=colors[rank-1], marker=markers[rank-1], markersize=4,
                            linewidth=1.3, label=fr"$\Lambda_{rank}$")
                    for x in sorted({float(r["x"]) for r in subset if r["status"]!="PASS"}):
                        ax.axvline(x, color="red", alpha=.4, linestyle=":")
                    ax.set_ylabel(style.get("ylabel", r"$\Lambda$")); ax.grid(alpha=.23)
                    ax.set_yscale(style.get("yscale", "linear"))
                    if "ylim" in style: ax.set_ylim(style["ylim"])
                    if "xlim" in style: ax.set_xlim(style["xlim"])
                    ax.legend(ncol=4, loc=style.get("legend_loc", "best"), fontsize=style.get("legend_font_size",10))
                axes[-1,0].set_xlabel(style.get("xlabel", {"BETA":r"$\beta$ (deg)","MU":r"$\mu$","TAU":r"$\tau$","EPS":r"$\varepsilon$","EPSILON":r"$\varepsilon$"}[entry["mode"]]))
                fig.suptitle(title, fontsize=style.get("title_font_size",11))
                fig.text(.5,.01,f"Sorted positions 1–{k_plot}; markers: computed nodes; segments: visual guides", ha="center", fontsize=9)
                fig.tight_layout(rect=(0,.025,1,.96))
                stem = entry["name"]+"_"+layout
                captions[stem] = line_geometry.add_geometry_caption(fig)
                for ext in ("png", "pdf"): fig.savefig(out/(stem+"."+ext), dpi=style.get("dpi",220), bbox_inches="tight")
                plt.close(fig)
                gallery.extend([f"## {entry['name']} / {layout}", "", f"![{stem}]({stem}.png)", "", f"[Vector PDF]({stem}.pdf)", ""])
    (out/"README.md").write_text("\n".join(gallery), encoding="utf-8")
    publish(out/"manifest.json", dict(identity=spec, root_calculations=0,
        geometry_diagnostics=geometry_diagnostics, geometry_captions=captions,
        rendering_seconds=time.perf_counter()-start,
        artifacts={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file()}))
    print(json.dumps(dict(figures=str(out), root_calculations=0)), flush=True)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=CONFIG)
    parser.add_argument("--mode", choices=["BETA", "MU", "TAU", "ALL"], default="ALL")
    parser.add_argument("--name", help="Optional exact map name within the selected mode")
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--plot-only", type=Path, metavar="BUNDLE")
    action.add_argument("--tables-only", type=Path, metavar="BUNDLE")
    parser.add_argument("--style", type=Path, help="JSON rendering overrides; never changes the spectrum identity")
    parser.add_argument("--no-plots", action="store_true")
    parser.add_argument("--check-aspect", action="store_true", help="Run the bounded section controls only, without maps")
    args = parser.parse_args(argv)
    if args.tables_only:
        render(args.tables_only, {}, tables_only=True); return 0
    style = read(args.config).get("plot", {})
    if args.style: style.update(read(args.style))
    if args.plot_only:
        render(args.plot_only, style); return 0
    c, policy = configuration(args.config)
    if args.check_aspect:
        from scripts.lib.rl_full_rr_square_checks import control_gate
        gate = control_gate(c, policy)
        print(json.dumps({k:gate[k] for k in ("folder", "status")}))
        return 0 if gate["status"] == "PASS" else 2
    specs = [s for s in c["maps"] if (args.mode=="ALL" or s["mode"]==args.mode) and (args.name is None or s["name"]==args.name)]
    if not specs: parser.error("No matching maps")
    folder = run(c, policy, specs)
    if not args.no_plots: render(folder, style)
    return 0 if read(folder/"manifest.json")["status"]=="PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
