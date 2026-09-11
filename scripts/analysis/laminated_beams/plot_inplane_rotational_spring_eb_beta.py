"""Five fixed EB joint states, one common angular table, six figures.

Diagnostic frequency-map-v1/fast_plot continuation, not the pilot's full scan
in a loop. Physics and detector/refiner are reused unchanged. plot-only never
calls a matrix provider, detector, SVD, or refiner and never writes spectra.
"""
from __future__ import annotations

import argparse
from collections import OrderedDict, Counter
from dataclasses import asdict
from datetime import datetime, timezone
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

for _thread_variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_thread_variable] = "1"
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/"src"))

import numpy as np
import scipy
from scipy.linalg import block_diag
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.analysis.laminated_beams import pilot_inplane_rotational_spring_eb as pilot

ARM = pilot.ARM
FS = pilot.FREQUENCY_SCALE
OUTPUT = ROOT/"results/laminated_beams/inplane_rotational_spring_eb_beta"
OLD = ROOT/"results/laminated_beams/inplane_rotational_spring_eb_pilot"
STATES = (("k0",0.),("k0.1",.1),("k1",1.),("k100",100.),("RIGID",None))
GOOD = ("COMPLETED","TARGET_CONFIRMED_GUARD_QUALIFIED")
LIMITS = dict(sigma_ratio=1e-9,rank_rtol=1e-12,physical_residual=1e-9,
    compatibility=1e-10,null_residual=1e-9,frequency_relative=1e-6,
    guard_margin_Omega=.02,max_B_per_point=6000,max_recoveries=30,
    repeated_failure_recoveries=3,max_extra_angles=20,transfer_cache_size=512,
    window_points=9,gap_points=5,window_floor=.04,window_relative=.002,
    neighbour_relative_defect=.002,neighbour_slope_jump_per_degree=.01)
CALLS = Counter()


def grid_tenths():
    return [*range(0,101),*range(105,301,5),*range(310,901,10)]


def case(state, beta_deg, role="BASE"):
    kappa = dict(STATES)[state]
    return dict(point_id=f"{state}_b{beta_deg:g}",state=state,beta_deg=float(beta_deg),
                beta_rad=math.radians(beta_deg),mode="RIGID" if kappa is None else "SPRING",
                kappa_theta=kappa,k_theta=None if kappa is None else kappa*ARM.D/ARM.L,
                grid_role=role)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def contract():
    return dict(frequency_map_policy="frequency-map-v1",calculation_mode="fast_plot",
        spectrum_semantics="sorted_positions",sweep_parameter="beta_deg",
        parameter_grid_tenths=grid_tenths(),states=list(STATES),K_plot=6,K_guard=7,
        guard_root_role="completeness_only",neighbour_audit="nonuniform_interpolation_defect_and_slope_jump",
        local_repair_policy="one_triggered_attempt_per_point_max30",strict_audit_default=False,
        geometry=dict(l=1,b=.20,h=.05,E=1,rho=1),arm=asdict(ARM),
        normalization="Omega=omega*l^2*sqrt(rho*Ag/(E*Ig)); Lambda=sqrt(Omega)",limits=LIMITS,
        physics_sha256={p:sha(ROOT/p) for p in ("scripts/lib/inplane_rotational_spring_eb.py",
            "scripts/lib/reddy_inplane_geometry.py","scripts/lib/reddy_symmetric_coupled_beams.py")})


class CostLimit(RuntimeError):
    pass


class Transfers:
    """Bounded exact-key cache: only the fixed identical EB arms of this map."""
    def __init__(self):
        self.cache = OrderedDict()
        self.expm_calls = self.hits = 0

    def endpoints(self, omega):
        key = (float(omega), ARM.A, ARM.D, ARM.m, ARM.L)
        if key in self.cache:
            self.hits += 1
            self.cache.move_to_end(key)
            return self.cache[key]
        self.expm_calls += 1
        CALLS["transfer_expm"] += 1
        block = eb.state_scale(ARM)[:,None]*eb._scaled_transfer(omega,ARM)[:,3:]
        value = block_diag(block,block)
        self.cache[key] = value
        if len(self.cache)>LIMITS["transfer_cache_size"]:
            self.cache.popitem(last=False)
        return value


class Provider:
    def __init__(self, point, transfers):
        self.point,self.transfers = point,transfers
        self.joint = eb.Joint(point["mode"],point["k_theta"])
        self.joint_matrix = eb.joint_matrix(point["beta_rad"],self.joint)
        self.reactions = np.tile(eb.state_scale(ARM)[3:],2)
        moment,force = ARM.D/ARM.L,ARM.D/ARM.L**2
        self.units = np.array([ARM.L,ARM.L,1. if point["mode"]=="RIGID" else moment,force,force,moment])
        self.factors = 1/self.units
        if point["mode"]=="SPRING":
            self.factors[2] /= max(1.,point["kappa_theta"])
        self.builds,self.cache = 0,{}

    def assembly(self, omega):
        key = float(omega)
        if key not in self.cache:
            if self.builds>=LIMITS["max_B_per_point"]:
                raise CostLimit("COST_LIMIT")
            self.builds += 1
            CALLS["matrix_builds"] += 1
            endpoints = self.transfers.endpoints(key)
            reacted = self.joint_matrix@endpoints
            self.cache[key] = eb.BoundaryAssembly(reacted/self.reactions[None,:],
                self.factors[:,None]*reacted,endpoints,self.reactions,self.units,self.factors)
        return self.cache[key]

    def __call__(self,omega):
        return self.assembly(omega).dimensionless


def diagnostics(provider, Omega):
    d = eb.endpoint_diagnostics(provider.assembly(Omega/FS),provider.point["beta_rad"],provider.joint,ARM)
    physical = np.max(np.abs([v["normalized_physical_residuals"] for v in d["vectors"]]),axis=0)
    residual = max(max(v["boundary_residual"],v["scaled_residual"]) for v in d["vectors"])
    failures = []
    if d["nullity"]<1 or d["sigma_ratio"]>LIMITS["sigma_ratio"]:
        failures.append("ENDPOINT_SINGULARITY_FAIL")
    if max(physical)>LIMITS["physical_residual"]:
        failures.append("PHYSICAL_RESIDUAL_FAIL")
    if max(physical[:2])>LIMITS["compatibility"]:
        failures.append("COMPATIBILITY_FAIL")
    if residual>LIMITS["null_residual"]:
        failures.append("NULL_RESIDUAL_FAIL")
    return dict(Omega=Omega,nullity=d["nullity"],sigma_ratio=d["sigma_ratio"],
                physical_residuals=physical.tolist(),null_residual=residual,failures=failures)


def slots_from_events(events):
    return [event for event in events for _ in range(event.diagnostics.detected_nullity)]


def merge_windows(windows):
    merged=[]
    for lo,hi in sorted(windows):
        if merged and lo<=merged[-1][1]:
            merged[-1][1]=max(merged[-1][1],hi)
        else:
            merged.append([lo,hi])
    return merged


def continuation_windows(history, beta):
    last=history[-1]
    values=np.array([r["Omega"] for r in last["rows"] if r["role"] in ("ROOT","GUARD","GUARD_CANDIDATE")])
    shift=np.zeros(len(values))
    if len(history)>1:
        previous=history[-2]
        old=np.array([r["Omega"] for r in previous["rows"]])
        if len(old)==len(values):
            shift=(values-old)*(beta-last["case"]["beta_deg"])/(last["case"]["beta_deg"]-previous["case"]["beta_deg"])
    forecast=values+shift
    half=np.maximum(LIMITS["window_floor"],LIMITS["window_relative"]*values)+2*np.abs(shift)
    windows=merge_windows([(max(1e-8,min(a,b)-w),max(a,b)+w) for a,b,w in zip(values,forecast,half)])
    upper=windows[-1][1]+max(.05,.002*values[-1])
    return windows,upper


def scan_interval(provider,lo,hi,points,phase="PRIMARY"):
    CALLS["detector_calls"]+=1
    return pilot.roots._scan_candidates(provider,FS,pilot.policy(lo,hi),
        case_id=provider.point["point_id"],builder_id="physical_EB_spring_map",
        scan_id=phase,points=points,phases=(0.,))[0]


def search_intervals(windows,upper):
    # Sparse gap checks include the lower spectrum; no dense global rescan.
    intervals=[];left=1e-8
    for lo,hi in windows:
        if lo>left+1e-10:
            intervals.append((left,lo,LIMITS["gap_points"]))
        intervals.append((lo,hi,LIMITS["window_points"]))
        left=hi
    if left<upper:
        intervals.append((left,upper,LIMITS["gap_points"]))
    return intervals


def assess(pool,provider,upper):
    pool,evidence=pilot.reconcile_local_detections(pool,provider)
    events,ambiguous=pilot.consolidate(pool)
    slots=slots_from_events(events)
    suspect=[c for c in pool if pilot.suspicious(c)]+ambiguous
    rows=[];checks=[];warnings=[]
    if len(slots)<6:
        return dict(status="MISSING_TARGET",rows=rows,endpoints=checks),pool,evidence,suspect
    sixth=slots[5].omega_bar
    bad_target=[c for c in suspect if c.interval_left_bar<=sixth+LIMITS["guard_margin_Omega"]]
    if len(slots)>=7:
        guard=slots[6].omega_bar
        kept=[c for c in slots if c.omega_bar<=guard]
    else:
        guards=[c for c in suspect if c.interval_left_bar>sixth+LIMITS["guard_margin_Omega"]]
        if not guards:
            return dict(status="MISSING_GUARD",rows=rows,endpoints=checks),pool,evidence,suspect
        candidate=min(guards,key=lambda c:c.omega_bar)
        guard=candidate.omega_bar
        kept=slots+[candidate]
    for index,event in enumerate(kept,1):
        d=diagnostics(provider,event.omega_bar)
        if event.accepted and d["nullity"]!=event.diagnostics.detected_nullity:
            d["failures"].append("MULTIPLICITY_MISMATCH")
        checks.append(d)
        rows.append(dict(sorted_position=index,role="ROOT" if index<=6 else ("GUARD" if event.accepted else "GUARD_CANDIDATE"),
            Omega=event.omega_bar,omega=event.omega_bar/FS,Lambda=math.sqrt(event.omega_bar),
            multiplicity=event.diagnostics.detected_nullity,endpoint_accepted=not d["failures"]))
    if bad_target or any(d["failures"] for d in checks[:6]):
        status="TARGET_UNCONFIRMED"
    elif guard-sixth<=LIMITS["guard_margin_Omega"] or upper-guard<=LIMITS["guard_margin_Omega"]:
        status="GUARD_NOT_SEPARATED"
    elif any(c.omega_bar<=guard for c in suspect) or any(d["failures"] for d in checks[6:]):
        status="TARGET_CONFIRMED_GUARD_QUALIFIED"
    else:
        status="COMPLETED"
    for c in suspect:
        if c.omega_bar<=guard:
            warnings.append(pilot.candidate_record(c))
    return dict(status=status,rows=rows,endpoints=checks,warnings=warnings,
                guard_gap_Omega=upper-guard,target_guard_gap_Omega=guard-sixth),pool,evidence,bad_target


def recovery_allowed(state,point,signature):
    return (point["point_id"] not in state["recovery_points"] and
            len(state["recovery_points"])<LIMITS["max_recoveries"] and
            state["failure_recovery_counts"].get(point["state"]+":"+signature,0)<LIMITS["repeated_failure_recoveries"])


def solve_point(point,history,transfers,state):
    started=time.perf_counter();start_expm=transfers.expm_calls;start_hits=transfers.hits
    provider=Provider(point,transfers);pool=[];evidence=[];recovery=None
    result=dict(status="NOT_COMPLETED",rows=[],endpoints=[])
    if history:
        windows,upper=continuation_windows(history,point["beta_deg"])
        intervals=search_intervals(windows,upper)
    else:
        # Only a missing initial anchor, not the regular angular continuation.
        windows=[];upper=120.;intervals=[(1e-8,upper,1201)]
    try:
        for lo,hi,n in intervals:
            pool.extend(scan_interval(provider,lo,hi,n))
        result,pool,evidence,suspects=assess(pool,provider,upper)
        if result["status"] not in GOOD:
            signature=result["status"]
            if recovery_allowed(state,point,signature):
                state["recovery_points"].append(point["point_id"])
                key=point["state"]+":"+signature
                state["failure_recovery_counts"][key]=state["failure_recovery_counts"].get(key,0)+1
                recovery=dict(trigger=signature,original_warnings=result.get("warnings",[]),
                              original_candidates=[pilot.candidate_record(c) for c in pool])
                # One attempt, restricted to suspect/missing-root windows.
                if suspects:
                    repairs=merge_windows([(max(1e-8,c.interval_left_bar-.02),min(upper,c.interval_right_bar+.02)) for c in suspects])
                elif windows:
                    repairs=merge_windows([(max(1e-8,lo-(hi-lo)),hi+(hi-lo)) for lo,hi in windows])
                    upper=max(upper,repairs[-1][1]+.05)
                else:
                    repairs=[(1e-8,upper)]
                recovery["intervals"]=repairs
                recovery["builds_before"]=provider.builds
                for lo,hi in repairs:
                    local=scan_interval(provider,lo,hi,25,"RECOVERY")
                    pool=[c for c in pool if not lo<c.omega_bar<hi]+local
                result,pool,proof,_=assess(pool,provider,upper)
                evidence.extend(proof)
                recovery["builds"]=provider.builds-recovery["builds_before"]
    except CostLimit:
        result["status"]="COST_LIMIT"
    except (ValueError,RuntimeError,FloatingPointError,np.linalg.LinAlgError) as error:
        result.update(status="NUMERICAL_FAILURE",error=str(error))
    result.update(case=point,origin="NEW_COMPUTATION",matrix_builds=provider.builds,
        expm_calls=transfers.expm_calls-start_expm,transfer_cache_hits=transfers.hits-start_hits,
        seconds=time.perf_counter()-started,search_intervals=intervals,recovery=recovery,
        detector_reconciliations=evidence,candidates=[pilot.candidate_record(c) for c in pool])
    return result


def atomic(path,text):
    tmp=path.with_suffix(path.suffix+".tmp")
    tmp.write_text(text,encoding="utf-8",newline="\n");tmp.replace(path)


def table_rows(state):
    rows=[]
    for group in state["points"].values():
        source_rows=group["rows"]
        if group["status"] not in GOOD:
            # Failed frequencies stay in diagnostics. Explicit blank rows also
            # preserve gaps at ADDED angles, which are absent from the BASE grid.
            source_rows=[dict(sorted_position=j,role="ROOT",Omega=None,omega=None,Lambda=None,
                              multiplicity=None,endpoint_accepted=False) for j in range(1,7)]
        for r in source_rows:
            rows.append(dict(**group["case"],**r,status=group["status"],origin=group["origin"],
                target_confirmed=group["status"] in GOOD,source_version=group.get("source_version","EB_PILOT")))
    return sorted(rows,key=lambda r:(r["state"],r["beta_deg"],r["sorted_position"]))


def save(state):
    OUTPUT.mkdir(parents=True,exist_ok=True)
    atomic(OUTPUT/"diagnostics.json",json.dumps(state,ensure_ascii=False,indent=1,allow_nan=False)+"\n")
    rows=table_rows(state)
    if rows:
        stream=io.StringIO(newline="");writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        atomic(OUTPUT/"spectrum_roots.csv",stream.getvalue())
    counts=Counter(g["status"] for g in state["points"].values())
    manifest={k:v for k,v in state.items() if k!="points"}
    manifest.update(status_counts=dict(counts),BASE_processed=sum(g["case"]["grid_role"]=="BASE" for g in state["points"].values()),
        BASE_confirmed=sum(g["case"]["grid_role"]=="BASE" and g["status"] in GOOD for g in state["points"].values()),
        reused_points=sum(g["origin"]=="REUSED_EB_REFERENCE" for g in state["points"].values()),
        newly_computed_points=sum(g["origin"]=="NEW_COMPUTATION" for g in state["points"].values()),
        matrix_builds=sum(g["matrix_builds"] for g in state["points"].values()),
        expm_calls=sum(g["expm_calls"] for g in state["points"].values()),
        spectrum_seconds=sum(g["seconds"] for g in state["points"].values()))
    path=OUTPUT/"run_manifest.json"
    if path.exists():
        old=json.loads(path.read_text(encoding="utf-8"))
        if "render" in old:manifest["render"]=old["render"]
    atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+"\n")


def reuse(state):
    paths=[OLD/name for name in ("spectrum_roots.csv","diagnostics.json","run_manifest.json")]
    if not all(p.exists() for p in paths):
        state["reuse_note"]="Local EB files unavailable; only missing anchors will be computed"
        return
    old=json.loads(paths[1].read_text(encoding="utf-8"));manifest=json.loads(paths[2].read_text(encoding="utf-8"))
    c=old["contract"]
    if c["geometry"]!=state["contract"]["geometry"] or c["arm"]!=asdict(ARM):
        raise ValueError("EB reuse geometry/properties mismatch")
    if c["normalization"]!="Omega=omega*l^2*sqrt(m/D); Lambda=sqrt(Omega)":
        raise ValueError("EB reuse normalization mismatch")
    with paths[0].open(encoding="utf-8",newline="") as stream:csv_rows=list(csv.DictReader(stream))
    state["reuse_origin"]=dict(files_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        source_HEAD=manifest["source_HEAD"],overall_status=manifest["status"],
        qualification="beta30_RIGID BASE accepted; additional legacy guard regression unresolved; kappa10000 not included")
    for label,kappa in STATES:
        for beta in (0.,30.):
            point=case(label,beta);name=f"beta{beta:g}_{label}"
            if point["point_id"] in state["points"]:continue
            group=old["groups"].get(name)
            if not group or group["status"]!="COMPLETED":continue
            physics=state["contract"]["physics_sha256"]
            if any(group["source_code_sha256"].get(p)!=digest for p,digest in physics.items()):
                raise ValueError("EB saved physical assembly differs")
            source=[r for r in csv_rows if r["case_id"]==name]
            if len(source)!=len(group["rows"]) or group["unresolved_below_guard"]!=0:
                raise ValueError("EB BASE row/count qualification mismatch")
            rows=[]
            for saved,row in zip(group["rows"],source):
                if (any(float(row[k])!=saved[k] for k in ("Omega","omega","Lambda")) or
                    saved["mode"]!=point["mode"] or saved["kappa_theta"]!=kappa or
                    saved["k_theta"]!=point["k_theta"] or saved["beta_rad"]!=point["beta_rad"] or
                    not math.isclose(saved["Omega"],saved["omega"]*FS,rel_tol=1e-14) or
                    not math.isclose(saved["Lambda"]**2,saved["Omega"],rel_tol=1e-14)):
                    raise ValueError("EB saved scalar mismatch")
                rows.append({**{k:saved[k] for k in ("sorted_position","role","Omega","omega","Lambda","multiplicity")},"endpoint_accepted":True})
            state["points"][point["point_id"]]=dict(case=point,status="COMPLETED",origin="REUSED_EB_REFERENCE",rows=rows,
                matrix_builds=0,expm_calls=0,seconds=0.,source_group=name,source_code_sha256=group["source_code_sha256"],
                source_guard_gap_Omega=group["guard_gap_Omega"],
                qualification="legacy guard regression unresolved (BASE accepted)" if name=="beta30_RIGID" else "")


def load_compute_state():
    path=OUTPUT/"diagnostics.json";expected=json.loads(json.dumps(contract()))
    if path.exists():
        state=json.loads(path.read_text(encoding="utf-8"))
        if state["contract"]!=expected:raise ValueError("Spectral contract changed; preserve existing checkpoint")
    else:
        state=dict(contract=expected,points={},recovery_points=[],failure_recovery_counts={},source_versions=[],extra_points=[])
        reuse(state)
    version=dict(HEAD=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        working_tree_status=subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True),
        runner_sha256=sha(Path(__file__)),detector_sha256=sha(ROOT/pilot.CODE_FILES[4]),
        executable=sys.executable,versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        started_utc=datetime.now(timezone.utc).isoformat())
    state["source_versions"].append(version)
    return state


def compute_point(state,point,transfers):
    if point["point_id"] in state["points"]:return False
    history=sorted([g for g in state["points"].values() if g["case"]["state"]==point["state"] and
        g["case"]["beta_deg"]<point["beta_deg"] and g["status"] in GOOD],key=lambda g:g["case"]["beta_deg"])[-2:]
    result=solve_point(point,history,transfers,state)
    result["source_version"]=len(state["source_versions"])-1
    state["points"][point["point_id"]]=result
    save(state)
    return True


def audit(state):
    flags=[];stiffness=[]
    for label,_ in STATES:
        groups=sorted([g for g in state["points"].values() if g["case"]["state"]==label and g["case"]["grid_role"]=="BASE"],key=lambda g:g["case"]["beta_deg"])
        for a,b,c in zip(groups,groups[1:],groups[2:]):
            if any(g["status"] not in GOOD for g in (a,b,c)):continue
            x0,x1,x2=[g["case"]["beta_deg"] for g in (a,b,c)]
            values=[np.log([r["Omega"] for r in g["rows"][:6]]) for g in (a,b,c)]
            interpolation=values[0]+(values[2]-values[0])*(x1-x0)/(x2-x0)
            defect=np.abs(values[1]-interpolation)
            jump=np.abs((values[2]-values[1])/(x2-x1)-(values[1]-values[0])/(x1-x0))
            indices=np.flatnonzero((defect>LIMITS["neighbour_relative_defect"]) & (jump>LIMITS["neighbour_slope_jump_per_degree"]))
            if len(indices):
                interval=(x0,x1) if x1-x0>=x2-x1 else (x1,x2)
                flags.append(dict(state=label,interval=interval,beta_mid=sum(interval)/2,
                    positions=(indices+1).tolist(),max_log_defect=float(max(defect)),max_slope_jump=float(max(jump))))
    for tick in grid_tenths():
        groups=[state["points"].get(case(label,tick/10)["point_id"]) for label,_ in STATES]
        for left,right in zip(groups,groups[1:]):
            if not left or not right or left["status"] not in GOOD or right["status"] not in GOOD:continue
            lo=np.array([r["Omega"] for r in left["rows"][:6]]);hi=np.array([r["Omega"] for r in right["rows"][:6]])
            if np.any(hi<lo*(1-LIMITS["frequency_relative"])):
                stiffness.append(dict(beta_deg=tick/10,left=left["case"]["state"],right=right["case"]["state"],relative_increment=((hi-lo)/lo).tolist()))
    return dict(neighbour_flags=flags,stiffness_violations=stiffness)


def repair_audit_point(state,point_id,transfers):
    """One existing-refiner attempt at a recorded stiffness-order trigger.

Search only the affected prediction window(s). Other saved events are
evaluated at their existing frequencies, never globally searched again.
"""
    group=state["points"][point_id];point=group["case"]
    triggers=[f for f in state.get("audit",{}).get("stiffness_violations",[])
              if f["beta_deg"]==point["beta_deg"] and f["left"]==point["state"]]
    if not triggers:raise ValueError("No recorded stiffness-order trigger for this point")
    if not recovery_allowed(state,point,"AUDIT_STIFFNESS_ORDER"):
        raise ValueError("Point/global recovery budget already used")
    state["recovery_points"].append(point_id)
    history=sorted([g for g in state["points"].values() if g["case"]["state"]==point["state"] and
        g["case"]["beta_deg"]<point["beta_deg"] and g["status"] in GOOD],key=lambda g:g["case"]["beta_deg"])[-2:]
    windows,_=continuation_windows(history,point["beta_deg"])
    indices=sorted({i for f in triggers for i,v in enumerate(f["relative_increment"]) if v < -LIMITS["frequency_relative"]})
    centers=[history[-1]["rows"][i]["Omega"] for i in indices]
    repairs=[(lo,hi) for lo,hi in windows if any(lo<=x<=hi for x in centers)]
    if not repairs:raise ValueError("Affected predictor window not located")
    started=time.perf_counter();start_expm=transfers.expm_calls
    provider=Provider(point,transfers);provider.builds=group["matrix_builds"]
    pool=[]
    for row in group["rows"]:
        Omega=row["Omega"]
        if any(lo<=Omega<=hi for lo,hi in repairs):continue
        saved=next(c for c in group["candidates"] if c["Omega"]==Omega and c["accepted"])
        diag=pilot.roots.boundary_matrix_diagnostics(Omega,provider,FS,
             rank_relative_tolerance=LIMITS["rank_rtol"],root_ratio_tolerance=LIMITS["sigma_ratio"])
        ok,reason=pilot.roots._candidate_quality(diag,pilot.policy(*saved["interval"]))
        pool.append(pilot.roots.RootCandidate(point_id,"EB_saved_event","AUDIT_EVALUATION",Omega,
            tuple(saved["sources"]),*saved["interval"],True,diag,ok,reason))
    for lo,hi in repairs:pool.extend(scan_interval(provider,lo,hi,129,"AUDIT_RECOVERY"))
    upper=group["search_intervals"][-1][1]
    result,pool,evidence,_=assess(pool,provider,upper)
    result.update(case=point,origin="NEW_COMPUTATION",source_version=len(state["source_versions"])-1,
        matrix_builds=provider.builds,expm_calls=group["expm_calls"]+transfers.expm_calls-start_expm,
        seconds=group["seconds"]+time.perf_counter()-started,search_intervals=group["search_intervals"],
        candidates=[pilot.candidate_record(c) for c in pool],detector_reconciliations=evidence,
        recovery=dict(trigger="AUDIT_STIFFNESS_ORDER",intervals=repairs,points=129,
                      builds=provider.builds-group["matrix_builds"],seconds=time.perf_counter()-started),
        previous_attempt=group)
    state["points"][point_id]=result
    save(state)
    print(point_id,result["status"],result["recovery"],flush=True)


def curve_arrays(rows,label,position):
    relevant=[r for r in rows if r["state"]==label and int(r["sorted_position"])==position]
    mapping={float(r["beta_deg"]):r for r in relevant}
    x=sorted(set(t/10 for t in grid_tenths())|set(mapping))
    y=[float(mapping[a]["Lambda"]) if a in mapping and mapping[a]["status"] in GOOD else math.nan for a in x]
    return np.array(x),np.array(y)


def render(output=OUTPUT):
    started=time.perf_counter();before=dict(CALLS)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with (output/"spectrum_roots.csv").open(encoding="utf-8",newline="") as stream:rows=list(csv.DictReader(stream))
    colors=["#0072B2","#D55E00","#009E73","#CC79A7","#222222"]
    styles=["-",(0,(2,2)),(0,(6,2,1,2)),(0,(1,1)),(0,(7,3))]
    markers=["o","s","^","D",None]
    labels=["κθ=0","κθ=0.1","κθ=1","κθ=100","жёсткое соединение"]
    with plt.rc_context({"font.family":"DejaVu Sans","font.size":12,"text.usetex":False,"pdf.fonttype":42}):
        for position in range(1,7):
            fig,ax=plt.subplots(figsize=(8.6,5.4));fig.subplots_adjust(left=.10,right=.98,bottom=.22,top=.97)
            for index,(label,_) in enumerate(STATES):
                x,y=curve_arrays(rows,label,position)
                ax.plot(x,y,color=colors[index],linestyle=styles[index],lw=1.65,
                    marker=markers[index],markersize=3.5,markerfacecolor="none",markevery=19,label=labels[index])
            ax.set(xlim=(0,90),xlabel="β, °",ylabel="Λ"+"₀₁₂₃₄₅₆"[position])
            ax.xaxis.label.set_size(15);ax.yaxis.label.set_size(15)
            ax.set_xticks(np.arange(0,91,10));ax.margins(y=.07);ax.grid(alpha=.23)
            fig.legend(*ax.get_legend_handles_labels(),loc="lower center",ncol=3,frameon=False,fontsize=10.5,bbox_to_anchor=(.5,.015))
            name=f"eb_spring_lambda{position:02d}_vs_beta"
            fig.savefig(output/(name+".png"),dpi=300)
            fig.savefig(output/(name+".pdf"))
            plt.close(fig)
    if dict(CALLS)!=before:raise AssertionError("plot-only called numerical code")
    result=dict(seconds=time.perf_counter()-started,calculation_mode="plot_only",renderer_sha256=sha(Path(__file__)),
                matplotlib=matplotlib.__version__,matrix_calls=0,root_calls=0,
                spectral_csv_sha256=sha(output/"spectrum_roots.csv"),figures=6,formats=["PDF","PNG 300 dpi"])
    path=output/"run_manifest.json"
    manifest=json.loads(path.read_text(encoding="utf-8"));manifest["render"]=result
    atomic(path,json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("plot-only",result,flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode",choices=("compute","benchmark","audit","plot-only"),required=True)
    parser.add_argument("--state",choices=[s for s,_ in STATES])
    parser.add_argument("--beta-min",type=float,default=0.)
    parser.add_argument("--beta-max",type=float,default=90.)
    parser.add_argument("--add-neighbours",action="store_true",help="audit only: one midpoint level, at most 20 triggered combinations")
    parser.add_argument("--repair-audit-point",help="audit only: one recorded stiffness-order trigger, without a global rescan")
    args=parser.parse_args()
    if args.mode=="plot-only":render();return
    if not 0<=args.beta_min<=args.beta_max<=90:parser.error("beta range must be within [0,90] degrees")
    state=load_compute_state();transfers=Transfers();save(state)
    if args.mode=="audit":
        if args.repair_audit_point:
            repair_audit_point(state,args.repair_audit_point,transfers)
        if "audit" in state:state.setdefault("audit_history",[]).append(state["audit"])
        state["audit"]=audit(state)
        if args.add_neighbours:
            for flag in sorted(state["audit"]["neighbour_flags"],key=lambda f:-f["max_log_defect"]):
                p=case(flag["state"],flag["beta_mid"],"ADDED")
                if len(state["extra_points"])>=LIMITS["max_extra_angles"]:break
                if p["point_id"] in state["points"]:continue
                state["extra_points"].append(dict(**flag,point_id=p["point_id"]))
                compute_point(state,p,transfers)
        save(state);print(json.dumps(state["audit"],ensure_ascii=False),flush=True);return
    labels=[args.state] if args.state else [s for s,_ in STATES]
    ticks=grid_tenths()
    if args.mode=="benchmark":labels=["k1"];ticks=[1,2,3]
    completed=0;started=time.perf_counter()
    for label in labels:
        for tick in ticks:
            if not args.beta_min<=tick/10<=args.beta_max:continue
            p=case(label,tick/10)
            if compute_point(state,p,transfers):
                completed+=1;g=state["points"][p["point_id"]]
                if args.mode=="benchmark" or completed%10==0 or g["status"] not in GOOD:
                    print(p["point_id"],g["status"],"B",g["matrix_builds"],"expm",g["expm_calls"],"s",round(g["seconds"],3),
                        "processed",len(state["points"]),"recoveries",len(state["recovery_points"]),flush=True)
    state["last_compute_invocation"]=dict(new_points=completed,wall_seconds=time.perf_counter()-started,
                                        calls=dict(CALLS),transfer_cache_size=len(transfers.cache))
    save(state)
    print("finished",state["last_compute_invocation"],flush=True)


if __name__=="__main__":main()
