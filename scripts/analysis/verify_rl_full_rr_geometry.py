"""Bounded C0/C1 extension of RL_FULL/RR K8; no epsilon or tau sweep.

Separate CLI: individual-section contracts, reflection gate and explicit
external archive references differ from immutable A-B output. Reuses its
catalogue builder, counted-prefix algorithm and serialization. Diagnostic only.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import csv
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
import scipy
from scripts.lib import rl_full_rr_sections_prefix as e
from scripts.analysis import verify_rl_full_rr_prefix as ab

p=e.p
CONFIG=ROOT/"data/input/rl_full_rr_geometry_c.json"
OUTPUT=ROOT/"results/rl_full_rr_geometry"
sha,write_new,jsonable=ab.sha,ab.write_new,ab.jsonable


def read_config():
    c=json.loads(CONFIG.read_text())
    original=json.loads(ab.CONFIG.read_text())
    if (c["material"]!=original["material"] or c["policy"]!=original["policy"] or
        c["K_plot"]!=8 or c["K_guard"]!=9 or c["reflection_geometries"]!=[[.1,-.1],[-.1,.1]] or
        c["reflection_angles"]!=[0.,30.] or c["pilot_mu"]!=[0.,.01,.1] or
        c["pilot_angles"]!=[0.,5.,6.1194,7.,30.]):
        raise ValueError("Only the authorized C0/C1 geometry and unchanged A-B policy")
    return c


def reference(c):
    folder=ROOT/c["reference"]
    m,points=ab.verify_bundle(folder)
    if m["status"]!="PASS" or m["stages"]!={"A":"PASS","B":"PASS"}:
        raise e.PrefixFailure("REFERENCE_NOT_ACCEPTED")
    for name,h in m["identity"]["sources"].items():
        if sha(ROOT/name)!=h: raise e.PrefixFailure("HISTORICAL_SOURCE_CHANGED",path=name)
    return folder,m,points


def provenance(c):
    folder,m,_=reference(c)
    paths={CONFIG,Path(__file__),ROOT/"tests/test_rl_full_rr_sections_prefix.py"}
    for mod in tuple(sys.modules.values()):
        name=getattr(mod,"__file__",None)
        if name:
            path=Path(name).resolve()
            if path.is_relative_to(ROOT) and path.suffix==".py":paths.add(path)
    return dict(schema=c["schema"],engine=e.VERSION,config=c,
        reference_manifest_sha256=sha(folder/"manifest.json"),
        reference_artifacts=m["artifacts"],
        sources={x.relative_to(ROOT).as_posix():sha(x) for x in sorted(paths)},
        head=subprocess.check_output(["git","rev-parse","HEAD"],text=True,cwd=ROOT).strip(),
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))


def ritz_check(model,case,roots,initial,policy):
    target=np.array([r["omega"] for r in roots]); trials=list(initial)
    for order in [*policy["ritz_orders"],policy["ritz_extra_order"]]:
        if not any(t["order"]==order for t in trials):
            w,_,d=p.ritz(model,case,order,count=len(target))
            trials.append(dict(order=order,omega=w.tolist(),diagnostics=d))
        trials.sort(key=lambda t:t["order"])
        cur=next(t for t in trials if t["order"]==order)
        prev=[t for t in trials if t["order"]<order]
        cur["relative_frequency_error"]=float(max(abs(np.array(cur["omega"])/target-1)))
        cur["relative_p_change"]=float(max(abs(np.array(cur["omega"])/np.array(prev[-1]["omega"])-1))) if prev else None
        d=cur["diagnostics"]
        cur["accepted"]=bool(prev and cur["relative_p_change"]<=policy["ritz_frequency_tol"] and
            cur["relative_frequency_error"]<=policy["ritz_frequency_tol"] and
            d["mass_gram_error"]<=policy["mass_gram_tol"] and d["energy_error"]<=policy["energy_tol"] and
            d["constraint_residual"]<=policy["boundary_tol"] and d["positive_K_M"])
        if cur["accepted"]:return dict(status="PASS",accepted_order=order,trials=trials)
    return dict(status="RITZ_NOT_CONVERGED",trials=trials)


def solve_point(c,stage,mu,tau,beta,folder,key,reference_records=(),neighbour=None):
    model=p.Model(**c["material"],mu=mu,tau=tau)
    case=p.Case(beta,mu,p.base.RIGID,p.base.RIGID,tau)
    policy=c["policy"]; contract=e.point_contract(model,case,8,9,policy)
    point_key=e.digest(dict(contract=contract,source_id=key,stage=stage))[:16]
    path=folder/"points"/(point_key+".json")
    if path.exists():
        rec=json.loads(path.read_text()); ab.verify_point_record(rec)
        if rec["contract"]!=contract or rec["source_id"]!=key:raise e.PrefixFailure("POINT_CACHE_MISMATCH")
        if "shape_file" in rec:ab.load_modes(folder,rec)
        return path,rec
    rec=dict(stage=stage,mu=mu,tau=tau,beta_deg=beta,contract=contract,source_id=key,status="RUNNING")
    start=time.perf_counter()
    try:
        w0,_,d0=p.ritz(model,case,policy["ritz_orders"][-1],count=9)
        initial=[dict(order=policy["ritz_orders"][-1],omega=w0.tolist(),diagnostics=d0)]
        rec["initial_ritz"]=initial
        lo=float(model.omega(policy["Lambda_lower"]))
        hi=max(float(model.omega(policy["Lambda_upper"])),float(w0[-1])*policy["ritz_upper_padding"])
        attempts=[]
        for attempt in range(policy["max_bound_retries"]+1):
            cats,cat_attempts=e.catalogues(model,case,hi,policy)
            rec.update(catalogs=cats,catalog_attempts=cat_attempts)
            counter=e.PhysicalCounter(model,case,cats,policy)
            try:
                events,search=e.counted_prefix(lambda w:p.balance(p.boundary(model,case,w))[0],
                                               counter,lo,hi,9,policy,neighbour)
                break
            except e.PrefixFailure as ex:
                attempts.append(dict(attempt=attempt,lower=lo,upper=hi,code=ex.code,details=ex.details))
                rec["bound_attempts"]=attempts
                if ex.code=="INSUFFICIENT_RIGHT_COUNT":hi*=policy["upper_growth"]
                elif ex.code=="NONZERO_LOWER_COUNT":lo/=2
                else:raise
        else:raise e.PrefixFailure("BOUND_BUDGET_EXHAUSTED",attempts=attempts)
        rec.update(search=search,bound_attempts=attempts)
        modes,roots=[],[]
        for event_id,event in enumerate(events,1):
            for j,mode in enumerate(e.recovered_modes(model,case,event,policy["quadrature_order"],policy)):
                rank=event["first_rank"]+j
                roots.append(dict(omega=event["omega"],Omega=event["omega"]*model.time,
                    Lambda=float(model.Lambda(event["omega"])),sorted_position=rank,event_id=event_id,
                    multiplicity=event["multiplicity"],role="research" if rank<=8 else ("guard" if rank==9 else "cluster_completion"),
                    coefficients=mode["coefficients"],diagnostics=mode["diagnostics"]))
                modes.append(mode)
        rec["roots"]=roots
        shape=path.with_suffix(".npz");shape.parent.mkdir(parents=True,exist_ok=True)
        with shape.open("xb") as stream:
            np.savez_compressed(stream,**{f"mode{i+1}_arm{j+1}_{field}":arm[field]
                for i,mode in enumerate(modes) for j,arm in enumerate(mode["arms"])
                for field in ("x","weights","y","dy","ends")})
        rec.update(shape_file=shape.relative_to(folder).as_posix(),shape_sha256=sha(shape),
            quality=e.mode_gates(model,modes,policy),
            quadrature=e.quadrature_check(model,case,roots,modes,policy),
            ritz=ritz_check(model,case,roots,initial,policy),
            domain=dict(status="PASS",arms=counter.limits,catalog_limit=counter.domain,
                        guard_omega=roots[8]["omega"],right_omega=search["guard_right_omega"]))
        checks=[rec[k]["status"] for k in ("quality","quadrature","ritz","domain")]
        if beta==0:
            rec["straight_interface"]=e.straight_interface(model,case,roots,modes,policy)
            checks.append(rec["straight_interface"]["status"])
            if tau==0:
                rec["straight_uniform"]=ab.straight_split(model.arms[0],case,roots,modes,hi,policy)
                checks.append(rec["straight_uniform"]["status"])
        regressions=[]
        for old in reference_records:
            ref_folder=ROOT/c["reference"]; ref_modes=ab.load_modes(ref_folder,old)
            n=len(old["roots"])
            cmp=e.compare_modes(model,roots[:n],modes[:n],old["roots"],ref_modes,policy)
            samples=[]
            for query in old["search"]["samples"]:
                if query["status"]!="PASS":continue
                actual,count,_=counter(query["requested"])
                samples.append(dict(requested=query["requested"],effective=actual,count=count,
                                    old_effective=query["effective"],old_count=query["count"]))
            cmp.update(reference_stage=old["stage"],count_checks=samples)
            if any(x["count"]!=x["old_count"] or x["effective"]!=x["old_effective"] for x in samples):
                cmp["status"]="COUNT_REGRESSION_FAILED"
            regressions.append(cmp);checks.append(cmp["status"])
        if regressions:rec["regressions"]=regressions
        rec["status"]="PASS" if all(x=="PASS" for x in checks) else "CHECK_FAILED"
        rec["mode_identification"]="SORTED_SUBSPACES" if any(x["multiplicity"]>1 for x in roots) else "SIMPLE_SORTED_NO_BRANCH_IDS"
    except (ArithmeticError,ValueError,RuntimeError,OSError,np.linalg.LinAlgError) as ex:
        rec.update(status=getattr(ex,"code","NUMERICAL_FAILURE"),
                   failure=getattr(ex,"details",dict(reason=str(ex),type=type(ex).__name__)))
    rec["elapsed_seconds"]=time.perf_counter()-start
    rec=jsonable(rec);rec["record_digest"]=e.digest(rec);write_new(path,rec)
    return path,rec


def resolve(entry):
    path=ROOT/entry["record"]
    if sha(path)!=entry["sha256"]:raise e.PrefixFailure("REFERENCED_POINT_CHANGED")
    rec=json.loads(path.read_text());ab.verify_point_record(rec)
    return rec,ab.load_modes(path.parent.parent,rec)


def entry(path,rec,origin):
    return dict(mu=rec.get("mu",rec["contract"]["case"]["mu"]),tau=rec["contract"]["tau"],
        beta_deg=rec["beta_deg"],record=path.relative_to(ROOT).as_posix(),sha256=sha(path),origin=origin)


def verify_bundle(folder):
    m=json.loads((folder/"manifest.json").read_text())
    for name,h in m["artifacts"].items():
        if sha(folder/name)!=h:raise e.PrefixFailure("ARTIFACT_HASH_MISMATCH",path=name)
    for item in m["C1_points"]:resolve(item)
    return m


def tables(folder,entries):
    rows=[]
    for item in entries:
        rec,_=resolve(item)
        for root in rec["roots"]:
            rows.append(dict(mu=item["mu"],tau=item["tau"],beta_deg=item["beta_deg"],
                rank=root["sorted_position"],omega=root["omega"],Omega=root["Omega"],Lambda=root["Lambda"],
                multiplicity=root["multiplicity"],role=root["role"],origin=item["origin"],source=item["record"]))
    with (folder/"spectrum_C1.csv").open("x",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def run(c,only_c0=False):
    provenance_record=provenance(c)
    # Both CLI stages share identity; C0 is sealed before C1 and always verified.
    key=e.digest(provenance_record)[:16];folder=OUTPUT/key
    if (folder/"manifest.json").exists():
        m=verify_bundle(folder);print(json.dumps(dict(folder=str(folder),status=m["status"],mode="verified_cache")));return
    folder.mkdir(parents=True,exist_ok=True)
    if not (folder/"identity.json").exists():write_new(folder/"identity.json",provenance_record)
    _,ref_manifest,ref_points=reference(c)
    points=[];known={};status="PASS";reflections=[]
    c0_cases=[(.01,0.,b) for b in [0.,5.,6.1194,7.,30.,90.]]
    c0_cases += [(mu,tau,b) for mu,tau in c["reflection_geometries"] for b in c["reflection_angles"]]
    c0_cases += [(0.,0.,0.)]
    gate_path=folder/"C0_gate.json"
    if gate_path.exists():
        gate=json.loads(gate_path.read_text())
        for item in gate["points"]:
            rec,_=resolve(item);path=ROOT/item["record"]
            points.append(item);known[item["mu"],item["tau"],item["beta_deg"]]=(path,rec)
        status=gate["status"]
    else:
        for mu,tau,beta in c0_cases:
            refs=[r for r in ref_points if r["beta_deg"]==beta] if (mu,tau)==(.01,0.) else []
            path,rec=solve_point(c,"C0",mu,tau,beta,folder,key,refs)
            item=entry(path,rec,"new_C0")
            points.append(item);known[mu,tau,beta]=(path,rec)
            print(json.dumps(dict(stage="C0",mu=mu,tau=tau,beta=beta,status=rec["status"],seconds=rec["elapsed_seconds"])),flush=True)
            if rec["status"]!="PASS":status="STOP_C0";break
        if status=="PASS":
            for beta in c["reflection_angles"]:
                pa,a=known[.1,-.1,beta];pb,b=known[-.1,.1,beta]
                ma=ab.load_modes(folder,a);mb=ab.load_modes(folder,b)
                model=p.Model(**b["contract"]["model"])
                cmp=e.compare_modes(model,b["roots"],mb,a["roots"],[e.reflect_mode(x) for x in ma],c["policy"])
                cmp.update(beta_deg=beta,from_record=pa.name,to_record=pb.name,
                           state_transform=e.REFLECTION_STATE.tolist(),arms_permuted=True)
                reflections.append(cmp)
                if cmp["status"]!="PASS":status="STOP_C0";break
        gate=dict(status=status,points=points,reflections=reflections)
        write_new(gate_path,gate)
    if only_c0 and status=="PASS":
        print(json.dumps(dict(folder=str(folder),C0="PASS",C1="NOT_RUN")));return
    pilot=[]
    if status=="PASS":
        for mu in c["pilot_mu"]:
            previous=None
            for beta in c["pilot_angles"]:
                if mu==.01:
                    rec=next(x for x in ref_points if x["stage"]=="B" and x["beta_deg"]==beta)
                    name=next(n for n in ref_manifest["point_files"] if json.loads((ROOT/c["reference"]/n).read_text())==rec)
                    item=entry(ROOT/c["reference"]/name,rec,"accepted_AB")
                elif (mu,0.,beta) in known:
                    path,rec=known[mu,0.,beta];item=entry(path,rec,"verified_C0")
                else:
                    path,rec=solve_point(c,"C1",mu,0.,beta,folder,key,neighbour=previous)
                    item=entry(path,rec,"new_C1")
                pilot.append(item)
                print(json.dumps(dict(stage="C1",mu=mu,beta=beta,status=rec["status"],origin=item["origin"])),flush=True)
                if rec["status"]!="PASS":status="STOP_C1";break
                previous=[r["omega"] for r in rec["roots"]]
            if status!="PASS":break
    reference(c)
    if pilot:tables(folder,pilot)
    manifest=dict(status=status,C0=gate["status"],C1="PASS" if status=="PASS" else ("NOT_RUN" if status=="STOP_C0" else "STOP"),
        identity=provenance_record,C0_points=points,C1_points=pilot,
        artifacts={x.relative_to(folder).as_posix():sha(x) for x in sorted(folder.rglob("*")) if x.is_file()},
        excluded=["epsilon_stage_D","tau_sweep","dense_angles","FEM"])
    write_new(folder/"manifest.json",manifest)
    print(json.dumps(dict(folder=str(folder),status=status,C0=manifest["C0"],C1=manifest["C1"])),flush=True)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--stage",choices=["C0","C"])
    group.add_argument("--verify-cache",type=Path)
    group.add_argument("--tables-only",type=Path)
    args=parser.parse_args(argv)
    if args.verify_cache or args.tables_only:
        folder=args.verify_cache or args.tables_only;m=verify_bundle(folder)
        print((folder/"spectrum_C1.csv").read_text() if args.tables_only else json.dumps(dict(status=m["status"],artifacts=len(m["artifacts"]))))
        return
    run(read_config(),args.stage=="C0")


if __name__=="__main__":main()

