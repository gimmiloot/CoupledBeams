"""Bounded seven-field verification, dispatched by the existing FEM CLI.

Stage A verifies the frozen action. Stage B performs four authorized 1D
cases. Only a frozen successful nonlinear-signal decision can enable the
sequential medium and conditional fine native pairs. Historical workflows
and failed-attempt guards are not reused as permission for this experiment.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CONFIG = ROOT / "data/input/nlsp_spatial_nonlinear_3d_fem_verification.json"
OUTPUT = ROOT / "results/nlsp_spatial_nonlinear_3d_fem_verification"
AUTHORIZATION = "explicit_user_NLSP_spatial_seven_field_verification_2026_10_10"
VERSION = "seven-field-straight-rod-bounded-verification-v1"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def plain(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(v) for v in value]
    return value


def read(path):
    return json.loads(Path(path).read_text(encoding="utf8"))


def write(path, value):
    Path(path).write_text(json.dumps(plain(value), ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf8")


def validate_config(c):
    if (c["schema"] != "nlsp-spatial-nonlinear-verification-v1"
        or c["authorization"] != AUTHORIZATION or c["degrees"] != [48, 64]
        or c["geometry"] != {"L": 1., "b": .2, "h": .1}
        or c["material"] != {"E": 1., "rho": 1., "nu": .3, "kappa": 5/6}
        or c["fields"] != ["u", "w", "v", "Phi", "psi", "theta", "c"]
        or c["local_basis_global"] != [[1,0,0],[0,-1,0],[0,0,-1]]
        or c["rotation_vector_local"] != ["Phi", "-psi", "theta"]
        or c["essential_BC"] != "all seven field values zero at s=0,L; no derivative constraints"
        or c["quadrature"] != "2*p+1 positive Gauss-Legendre points"
        or c["horizon_T1"] != .25 or c["omega1"] != .6054167303477958
        or c["one_d_cases"] != [["joint_p48",48,"joint"],["joint_p64",64,"joint"],
            ["isolated_w_p64",64,"w"],["isolated_v_p64",64,"v"]]
        or c["one_d_maximum_nonlinear_ODE_calls"] != 4
        or c["execution_mode"] != "EXPLORATORY_NOT_CERTIFIED" or c["admitted"] is not False
        or any(c[k] for k in ("new_meshes", "new_modal_FEM_jobs", "full_period_FEM"))):
        raise ValueError("Only the explicitly bounded seven-field contract is authorized")
    from scripts.lib import nlsp_spatial_fem_protocol as protocol
    if c["FEM_dynamic"] != protocol.DEFAULT_DYNAMIC or c["FEM_budget"] != protocol.RESOURCE_POLICY:
        raise ValueError("Frozen native time/resource policies changed")
    if c["load_policy"] != {"primary_w_over_h": .04, "fallback_w_over_h": .03,
        "g_k_over_g_n": 1.25, "combined_corner_bending_strain_limit": .01,
        "kind": "uniform dead reference body acceleration; one global resultant GRAV vector", "distributed_torque": 0.}:
        raise ValueError("The two predeclared load candidates must remain unchanged")
    if (c["comparison"]["planning_signal_ratio"] != 10.
        or c["comparison"]["physical_grid_points"] != 801
        or c["comparison"]["L2_gauss_points"] != 100
        or c["comparison"]["common_FEM_time_points"] != 201
        or c["comparison"]["primary_time_interpolation"] != "linear"
        or c["comparison"]["diagnostic_time_interpolation"] != "PCHIP"
        or c["comparison"]["relative_numerical_floor"] != 1e-10
        or c["comparison"]["phase_amplitude_time_fitting"] is not False
        or c["comparison"]["u_c_relative"] != 1e-3
        or c["comparison"]["bending_rotation_relative"] != 1e-4
        or c["comparison"]["energy_relative_drift"] != 1e-6):
        raise ValueError("Predeclared signal/physical-comparison policy changed")
    return c


def checked(source, name, registry):
    manifest = read(source/"manifest.json")
    expected = manifest.get("artifact_hashes", manifest.get("artifacts", {})).get(name)
    path = source/name
    if expected is None or sha(path) != expected:
        raise ValueError("Missing/corrupt immutable source: "+str(path))
    registry[path.relative_to(ROOT).as_posix()] = expected
    return path


def source_registry(c):
    registry, sources = {}, {}
    for name, evidence in c["source_bundles"].items():
        path = ROOT/evidence["path"]
        if sha(path/"manifest.json") != evidence["manifest_sha256"]:
            raise ValueError("Immutable source manifest changed: "+name)
        sources[name] = path
    selected = {"action": ["result.json"], "FEM1": ["preflight.json",
        "meshes/medium/solid_mesh.inp", "meshes/medium/mesh_audit.json",
        "meshes/fine/solid_mesh.inp", "meshes/fine/mesh_audit.json"],
        "FEM2": ["one_d_p64.npz", "one_d_preflight.json"],
        "FEM3AR": ["one_d_nonlinear.npz"],
        "FEM3C": ["provenance.json", "one_d_p48_nonlinear.npz", "one_d_p64_nonlinear.npz",
            "cases/medium_refined_time/nonlinear/section_history.npz",
            "cases/fine_refined_time/nonlinear/section_history.npz"],
        "profile_audit": ["summary.json", "numbers.json"]}
    for name, files in selected.items():
        for file in files:
            checked(sources[name], file, registry)
    if sha(Path(c["ccx_exe"])) != c["ccx_sha256"]:
        raise ValueError("Actual CalculiX binary differs from the verified build")
    return sources, registry


def load_discretizations(c):
    from scripts.lib import weakly_nonlinear_spatial_rod as rod
    from scripts.lib.weakly_nonlinear_spatial_dynamics import SpatialGalerkin
    action = read(ROOT/c["source_bundles"]["action"]["path"]/"result.json")["polynomials"]
    model = SimpleNamespace(T4=rod.Polynomial.deserialize(action["T4"]),
        V4=rod.Polynomial.deserialize(action["V4"]),
        residual_a=tuple(rod.Polynomial.deserialize(p) for p in action["residuals_A"]),
        symbols={n: rod.Polynomial.symbol(n) for n in rod.SYMBOL_ORDER})
    reference = read(ROOT/c["source_bundles"]["FEM1"]["path"]/"preflight.json")
    if reference["geometry"] != c["geometry"] or reference["material"] != c["material"]:
        raise ValueError("Frozen geometry/material references disagree")
    coefficients = rod.RodCoefficients(**reference["coefficients"])
    discs = {p: SpatialGalerkin(coefficients,p,length=1.,nq=2*p+1,model=model,whiten=True) for p in c["degrees"]}
    return discs, model, coefficients


def select_load(c, coefficients):
    from scripts.analysis.verify_nlsp_nonlinear_static_3d_fem import fem2_tim_uniform
    from scripts.lib.nlsp_spatial_fem_protocol import load_contract
    candidates = []
    for fraction in (.04, .03):
        unit_w = float(fem2_tim_uniform(np.array([.5]),1.,1.,coefficients.Bp,coefficients.S)[0][0,1])
        q_w = fraction*.1/unit_w
        q_v = 1.25*q_w
        corner = .05*abs(q_w/(12*coefficients.Bp))+.1*abs(q_v/(12*coefficients.Bb))
        candidates.append({"w_linear_over_h":fraction,"q_w":q_w,"q_v":q_v,
            "estimated_combined_corner_bending_strain":corner,"passed":corner<=.01})
        if corner<=.01:
            contract = load_contract(q_w/.02,q_v/.02)
            contract.update(q_w=q_w,q_v=q_v,w_linear_midspan=fraction*.1,
                v_linear_midspan=float(fem2_tim_uniform(np.array([.5]),q_v,1.,coefficients.Bb,coefficients.S)[0][0,1]),
                estimated_combined_corner_bending_strain=corner,candidate_history=candidates,
                selected_before_new_3D_results=True)
            return contract
    raise ValueError("Both authorized load candidates exceed the predeclared corner-strain guide")


def save(bundle, item, summary):
    write(bundle/"summary.json",summary)
    write(bundle/"manifest.json",{"schema":VERSION,"artifact_hashes":{
        p.relative_to(bundle).as_posix():sha(p) for p in sorted(bundle.rglob("*"))
        if p.is_file() and p != bundle/"manifest.json"},"source_manifests":item["config"]["source_bundles"]})


def validate_cache(bundle):
    bundle=Path(bundle)
    for name,digest in read(bundle/"manifest.json")["artifact_hashes"].items():
        if sha(bundle/name)!=digest:
            raise ValueError("Immutable verification cache artifact changed: "+name)
    item=read(bundle/"provenance.json")
    for evidence in item["config"]["source_bundles"].values():
        if sha(ROOT/evidence["path"]/"manifest.json")!=evidence["manifest_sha256"]:
            raise ValueError("Historical source manifest changed")
    for name,digest in item["source_artifacts"].items():
        if sha(ROOT/name)!=digest:
            raise ValueError("Selected historical source changed: "+name)
    return read(bundle/"summary.json")


def prepare(config=CONFIG):
    c=validate_config(read(config));sources,registry=source_registry(c)
    identity={"version":VERSION,"config":c,"source_artifacts":registry,
        "HEAD":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "python":sys.version,"dependencies":{n:importlib.metadata.version(n) for n in ("numpy","scipy","matplotlib")}}
    fingerprint=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
    bundle=OUTPUT/fingerprint
    # Authorization ledger is invariant under implementation-hash changes.
    for previous in OUTPUT.glob("*/provenance.json") if OUTPUT.exists() else []:
        old=read(previous)
        if old["config"]["authorization"]==AUTHORIZATION:
            if old["config"]!=c or previous.parent!=bundle:
                raise ValueError("A scientific attempt already exists; do not bypass its identity or failure guard")
    if (bundle/"manifest.json").exists():
        return bundle,read(bundle/"provenance.json"),validate_cache(bundle)
    if bundle.exists():
        raise ValueError("Interrupted unmanifested attempt; preserve evidence instead of automatic retry")
    bundle.mkdir(parents=True)
    identity["initial_checkout_snapshot"]="results/_smoke/nlsp_spatial_nonlinear_verification/initial_state.json"
    write(bundle/"provenance.json",identity);write(bundle/"config.json",c)
    summary={"overall":"PREPARED","completed":False,"hard_stop":False,
        "stage_A":"NOT_RUN","stage_B":"NOT_RUN","stage_C":"NOT_RUN",
        "one_d_cases":{},"FEM_cases":{},"attempts":[],"one_d_seconds":0.,"CCX_seconds":0.,
        "calls":{"nonlinear_1D_ODE":0,"CCX_production":0,"Gmsh":0,"FEM_modal":0},
        "historical_strict_float64":"PARTIAL","historical_statuses_unchanged":True}
    save(bundle,identity,summary);return bundle,identity,summary


def compute_stage_ab(bundle,item,summary,through_stage="B"):
    if summary["completed"] or summary["hard_stop"] or summary["stage_B"]!="NOT_RUN":
        return summary
    from scripts.lib import nlsp_spatial_verification_checks as checks
    from scripts.lib import nlsp_spatial_1d_program as program
    c=item["config"];discs,model,coefficients=load_discretizations(c)
    source=lambda key:ROOT/c["source_bundles"][key]["path"]
    paths={"fem1_preflight":source("FEM1")/"preflight.json","action_result":source("action")/"result.json",
        "short_planar":source("FEM3AR")/"one_d_nonlinear.npz",
        "full_planar_p48":source("FEM3C")/"one_d_p48_nonlinear.npz",
        "full_planar_p64":source("FEM3C")/"one_d_p64_nonlinear.npz",
        "static_p64":source("FEM2")/"one_d_p64.npz","static_preflight":source("FEM2")/"one_d_preflight.json"}
    if summary["stage_A"]=="NOT_RUN":
        phase=bundle/"execution_code"/"stage_A";phase.mkdir(parents=True)
        for module in (checks,program):shutil.copyfile(module.__file__,phase/Path(module.__file__).name)
        shutil.copyfile(ROOT/"scripts/lib/weakly_nonlinear_spatial_dynamics.py",phase/"weakly_nonlinear_spatial_dynamics.py")
        report=checks.run_stage_a_checks(discs,model,paths,output=bundle/"stage_a_checks.json")
        summary["stage_A"]=report["stage_a_status"]
        if report["execution_gate"]!="PASS":
            summary.update(hard_stop=True,overall="BLOCKED",stop_reason="Required new Stage A implementation/weak-form checks did not pass")
            save(bundle,item,summary);return summary
        save(bundle,item,summary)
    if through_stage=="A":
        return summary
    load=select_load(c,coefficients);write(bundle/"frozen_load.json",load)
    T=2*np.pi/c["omega1"];H=.25*T
    if (bundle/"one_d_common_times.npy").exists():times=np.load(bundle/"one_d_common_times.npy")
    else:
        omega_max=0.
        for p,disc in discs.items():
            modes=disc.linear_eigenpairs();omega_max=max(omega_max,float(modes["omega"].max()))
            np.savez_compressed(bundle/f"linear_modes_p{p}.npz",**modes,M0=disc.M0)
        count=int(np.ceil(H*omega_max/(2*np.pi)*12))
        times=np.unique(np.r_[np.linspace(0.,H,count+1),np.linspace(0.,H,201)])
        np.save(bundle/"one_d_common_times.npy",times)
        write(bundle/"one_d_sampling.json",{"omega_max":omega_max,"samples_per_fastest_retained_period":12,
            "samples":len(times),"T1":T,"H":H,"no_frequency_filter":True})
    for name,p,direction in c["one_d_cases"]:
        if summary["one_d_cases"].get(name,{}).get("status")=="PASS":continue
        if name in summary["one_d_cases"]:
            summary.update(hard_stop=True,overall="NUMERICAL_PARTIAL",stop_reason="An earlier 1D attempt cannot be automatically retried")
            save(bundle,item,summary);return summary
        qw=load["q_w"] if direction in ("joint","w") else 0.
        qv=load["q_v"] if direction in ("joint","v") else 0.
        summary["one_d_cases"][name]={"status":"STARTED","p":p,"loads":[qw,qv]}
        save(bundle,item,summary)
        started=time.perf_counter();print("Stage B seven-field 1D: "+name,flush=True)
        try:
            remaining=c["one_d_total_seconds"]-summary["one_d_seconds"]
            result=program.run_case(discs[p],(qw,qv),times,time.perf_counter()+min(c["one_d_per_case_seconds"],remaining),
                bundle/"one_d"/name,authorization={"user_authorized_spatial_stage_b":True,"id":AUTHORIZATION})
            summary["one_d_cases"][name]={"status":result["status"],"p":p,"case_json":"one_d/"+name+"/case.json"}
            if result["status"]!="PASS":
                summary.update(hard_stop=True,overall="NUMERICAL_PARTIAL",stop_reason=result.get("failure","Incomplete Stage B 1D case"))
        except Exception as error:
            summary["one_d_cases"][name].update(status="FAIL",failure=str(error))
            summary.update(hard_stop=True,overall="BLOCKED",stop_reason=str(error))
        finally:
            summary["one_d_seconds"]+=time.perf_counter()-started
            summary["calls"]["nonlinear_1D_ODE"]=sum(read(p)["calls"]["nonlinear_ODE"]
                for p in (bundle/"one_d").glob("*/case.json"))
            save(bundle,item,summary)
        if summary["hard_stop"]:return summary
    from scripts.lib import nlsp_spatial_comparison as comparison
    result=comparison.analyze_stage_b({n:bundle/"one_d"/n for n,_,_ in c["one_d_cases"]},discs,bundle,T1=T)
    summary["stage_B"]=result.get("status",result.get("overall","PARTIAL"))
    summary["overall"]="ONE_D_PREFLIGHT_COMPLETE"
    save(bundle,item,summary)
    return summary


def compute_additional_controls(bundle,item,summary):
    """Same excitation types at the second p; keep original four-case evidence.

    The user's Stage B permits isolated w/v controls and two resolutions. The
    original four-call cap was our primary-run organization, not a user cap.
    This separately recorded amendment permits exactly two additional calls
    within the unchanged 3600 s 1D budget; it does not certify full velocities.
    """
    if summary["completed"] or summary["hard_stop"] or summary["stage_B"]=="NOT_RUN":
        return summary
    if summary.get("additional_controls")=="PASS":return summary
    c=item["config"]
    names=(("isolated_w_p48",48,"w"),("isolated_v_p48",48,"v"))
    amendment=bundle/"additional_1d_controls_config.json"
    if not amendment.exists():
        baseline=bundle/"primary_four_case_evidence";baseline.mkdir()
        shutil.copyfile(bundle/"manifest.json",baseline/"phase_manifest.json")
        shutil.copyfile(bundle/"summary.json",baseline/"summary.json")
        for path in bundle.glob("stage_b*"):
            if path.is_file():shutil.copyfile(path,baseline/path.name)
        phase=bundle/"execution_code"/"additional_controls";phase.mkdir(parents=True)
        for name in ("nlsp_spatial_nonlinear_verification","nlsp_spatial_1d_program",
                     "weakly_nonlinear_spatial_dynamics","nlsp_spatial_comparison"):
            path=ROOT/"scripts/lib"/(name+".py");shutil.copyfile(path,phase/path.name)
        write(amendment,{"schema":"bounded-existing-excitation-p48-controls-v1",
            "authorization":AUTHORIZATION,"permission_origin":"Original user Stage B section12 allows isolated w/v/joint controls; section19 requires two p. This is an agent implementation decision, not a new user answer.",
            "primary_manifest_snapshot":"primary_four_case_evidence/phase_manifest.json",
            "primary_manifest_sha256":sha(baseline/"phase_manifest.json"),
            "original_config_sha256":sha(bundle/"config.json"),"frozen_load_sha256":sha(bundle/"frozen_load.json"),
            "cases":names,"maximum_additional_nonlinear_ODE_calls":2,
            "maximum_total_nonlinear_ODE_calls":6,"total_1D_seconds":c["one_d_total_seconds"],
            "unchanged_primary_cases":c["one_d_cases"],"unchanged_load_horizon_time_policy":True,
            "purpose":"Measure actual mixed nonlinear response sensitivity, including both isolated controls, before any FEM results.",
            "full_field_and_velocity_thresholds_unchanged":True,"automatic_retry":False})
        save(bundle,item,summary)
    evidence=read(amendment)
    if (evidence["cases"]!=[list(x) for x in names]
        or evidence["maximum_additional_nonlinear_ODE_calls"]!=2
        or evidence["original_config_sha256"]!=sha(bundle/"config.json")
        or evidence["primary_manifest_sha256"]!=sha(bundle/evidence["primary_manifest_snapshot"])
        or evidence["frozen_load_sha256"]!=sha(bundle/"frozen_load.json")):
        raise ValueError("Supplemental scope/source identity changed")
    from scripts.lib import nlsp_spatial_1d_program as program
    discs,_,_=load_discretizations(c);load=read(bundle/"frozen_load.json")
    times=np.load(bundle/"one_d_common_times.npy")
    for name,p,direction in names:
        if summary["one_d_cases"].get(name,{}).get("status")=="PASS":continue
        if name in summary["one_d_cases"]:
            summary.update(hard_stop=True,overall="NUMERICAL_PARTIAL",stop_reason="Supplemental attempt cannot be automatically retried")
            save(bundle,item,summary);return summary
        remaining=c["one_d_total_seconds"]-summary["one_d_seconds"]
        if remaining<=0:
            summary.update(hard_stop=True,overall="NUMERICAL_PARTIAL",stop_reason="Original 1D budget exhausted before additional control")
            save(bundle,item,summary);return summary
        loads=(load["q_w"] if direction=="w" else 0.,load["q_v"] if direction=="v" else 0.)
        summary["one_d_cases"][name]={"status":"STARTED","p":p,"loads":loads,
            "amendment_sha256":sha(amendment)}
        save(bundle,item,summary);started=time.perf_counter()
        print("Additional same-excitation 1D control: "+name,flush=True)
        try:
            result=program.run_case(discs[p],loads,times,
                time.perf_counter()+min(c["one_d_per_case_seconds"],remaining),bundle/"one_d"/name,
                authorization={"user_authorized_spatial_stage_b":True,"id":AUTHORIZATION,
                    "additional_existing_excitation_controls_sha256":sha(amendment)})
            summary["one_d_cases"][name].update(status=result["status"],case_json="one_d/"+name+"/case.json")
            if result["status"]!="PASS":summary.update(hard_stop=True,overall="NUMERICAL_PARTIAL",stop_reason=result.get("failure","Additional control incomplete"))
        except Exception as error:
            summary["one_d_cases"][name].update(status="FAIL",failure=str(error))
            summary.update(hard_stop=True,overall="NUMERICAL_PARTIAL",stop_reason=str(error))
        finally:
            summary["one_d_seconds"]+=time.perf_counter()-started
            summary["calls"]["nonlinear_1D_ODE"]=sum(read(p)["calls"]["nonlinear_ODE"] for p in (bundle/"one_d").glob("*/case.json"))
            save(bundle,item,summary)
        if summary["hard_stop"]:return summary
    from scripts.lib import nlsp_spatial_comparison as comparison
    all_cases={n:bundle/"one_d"/n for n,_,_ in (*c["one_d_cases"],*names)}
    result=comparison.analyze_stage_b(all_cases,discs,bundle,T1=2*np.pi/c["omega1"])
    summary["additional_controls"]="PASS"
    summary["stage_B"]=result["status"]
    save(bundle,item,summary);return summary


def freeze_fem_decision(bundle,item,summary):
    """Decide on primary nonlinear bending signals before any native result."""
    decision=bundle/"pre_fem_decision.json"
    if decision.exists():return read(decision)
    stage_a=read(bundle/"stage_a_checks.json");stage_b=read(bundle/"stage_b_comparison.json")
    c=item["config"]
    eligible=(stage_a["execution_gate"]=="PASS" and stage_b["stage_c_allowed"] is True
        and summary.get("additional_controls")=="PASS" and not summary["hard_stop"])
    estimates={"medium_linear":1122.,"medium_nonlinear":1185.,"fine_linear":2922.,"fine_nonlinear":3071.}
    # Historical counters are planning estimates for comparable quarter-period
    # meshes and steps, not promises or measurements of this new excitation.
    free=shutil.disk_usage(bundle).free
    result={"authorization":AUTHORIZATION,"selected_before_any_new_3D_result":True,
        "stage_C_allowed":eligible,"stage_A_execution_gate":stage_a["execution_gate"],
        "stage_B_primary_signal_gate":stage_b["stage_c_decision"],
        "full14_spatial_status":stage_b["all14_spatial"]["status"],
        "full_velocity_spatial_status":stage_b["all14_spatial"]["velocity"]["status"],
        "stage_B_overall_status_unchanged":summary["stage_B"],
        "qualification":"Limited primary displacement/evolving mixed-coupling permission; c and full velocities PARTIAL, single 1D time level, torsion not yet physically verified.",
        "original_primary_decision":"primary_four_case_evidence/stage_b_comparison.json",
        "original_primary_decision_sha256":sha(bundle/"primary_four_case_evidence/stage_b_comparison.json"),
        "amended_signal_analysis_sha256":sha(bundle/"stage_b_comparison.json"),
        "additional_controls_sha256":sha(bundle/"additional_1d_controls_config.json"),
        "frozen_load_sha256":sha(bundle/"frozen_load.json"),"config_sha256":sha(bundle/"config.json"),
        "source_manifests":c["source_bundles"],"source_artifacts":item["source_artifacts"],
        "dynamic_policy":c["FEM_dynamic"],"H":.25*2*np.pi/c["omega1"],
        "resource_policy":c["FEM_budget"],"planning_runtime_seconds":estimates,
        "planning_total_seconds":sum(estimates.values()),"free_disk_bytes":free,
        "code_identity":{name:sha(ROOT/"scripts/lib"/(name+".py")) for name in
            ("nlsp_spatial_nonlinear_verification","nlsp_spatial_native_program","nlsp_spatial_fem_protocol","nlsp_spatial_comparison")},
        "binary_sha256":c["ccx_sha256"],"no_physics_or_gate_denominator_changes":True,
        "maximum_new_production_CCX_jobs":4,"fine_pair_requires_actual_medium_signal":True}
    if free<8*1024**3:
        result.update(stage_C_allowed=False,stop_reason="Insufficient free disk for bounded native output (8 GiB planning reserve)")
    write(decision,result);return result


def compute_stage_c(bundle,item,summary):
    if summary["completed"] or summary["hard_stop"]:return summary
    if summary.get("additional_controls")!="PASS":return summary
    from scripts.lib import nlsp_spatial_fem_protocol as protocol
    from scripts.lib import nlsp_spatial_native_program as native
    from scripts.lib import nlsp_spatial_comparison as comparison
    decision=freeze_fem_decision(bundle,item,summary)
    if not decision["stage_C_allowed"]:
        summary.update(completed=True,stage_C="NOT_RUN",overall="SPATIAL_NONLINEAR_SIGNAL_NOT_RESOLVED",
            stop_reason=decision.get("stop_reason","Pre-FEM nonlinear bending/mixed signal gate not passed"))
        save(bundle,item,summary);return summary
    c=item["config"];load=read(bundle/"frozen_load.json");stage_b=read(bundle/"stage_b_comparison.json")
    sources,_=source_registry(c);discs,_,_=load_discretizations(c)
    remaining=c["FEM_budget"]["total_CCX_budget_seconds"]-summary["CCX_seconds"]
    for ordinal,(level,nonlinear) in enumerate((("medium",False),("medium",True),("fine",False),("fine",True)),1):
        name=level+("_nonlinear" if nonlinear else "_linear")
        if summary["FEM_cases"].get(name,{}).get("status")=="PASS":continue
        if name in summary["FEM_cases"]:
            summary.update(hard_stop=True,stage_C="PARTIAL",overall="NUMERICAL_PARTIAL",stop_reason="Native attempt already exists; no automatic retry")
            save(bundle,item,summary);return summary
        if level=="fine":
            medium=read(bundle/"comparison_medium/one_d_three_d_comparison.json")
            informative=all(medium["signal_planning"][f]["resolved_for_pre_FEM_planning"] for f in ("w","v"))
            if not informative:
                summary.update(completed=True,stage_C="PARTIAL",overall="SPATIAL_NONLINEAR_SIGNAL_NOT_RESOLVED",
                    stop_reason="Actual medium evolving bend corrections do not clear predeclared signal/interpolation indicators; conditional fine pair NOT_RUN")
                save(bundle,item,summary);return summary
            required=sum(decision["planning_runtime_seconds"][k] for k in ("fine_linear","fine_nonlinear") if k not in summary["FEM_cases"])
            if remaining<required:
                summary.update(completed=True,stage_C="PARTIAL",overall="NUMERICAL_PARTIAL",stop_reason="Completed-prefix cost leaves insufficient budget for estimated conditional fine pair")
                save(bundle,item,summary);return summary
        source=sources["FEM1"]/"meshes"/level
        audit=read(source/"mesh_audit.json")
        mesh=protocol.fem1.single.read_gmsh_inp_mesh_data(source/"solid_mesh.inp")
        auth={"id":AUTHORIZATION,"explicit_user_authorization":True,"stage_A_pass":True,"stage_B_pass":True,
            "frozen_decision_path":str(bundle/"pre_fem_decision.json"),"frozen_decision_sha256":sha(bundle/"pre_fem_decision.json"),
            "source_manifest_verified":True,"source_mesh_sha256":sha(source/"solid_mesh.inp"),
            "mesh_level":level,"nonlinear":nonlinear,"ordinal":ordinal,
            "completed_medium_pair":level=="fine","medium_signal_resolved":level=="fine"}
        summary["FEM_cases"][name]={"status":"STARTED","ordinal":ordinal}
        summary["attempts"].append({"case":name,"ordinal":ordinal,"status":"STARTED","automatic_retry":False})
        save(bundle,item,summary);print("New bounded native job: "+name,flush=True)
        try:
            result=native.run_attempt(bundle/"FEM"/name,source,mesh,audit,c,load,authorization=auth,remainingBudget=remaining)
            summary["FEM_cases"][name]={"status":result["status"],"ordinal":ordinal,
                "native_seconds":result["native_seconds"],"recovery_seconds":result["recovery_seconds"],
                "native_attempt":"FEM/"+name+"/native_attempt.json"}
            summary["attempts"][-1].update(status=result["status"],new_solver_calls=result["new_solver_calls"])
            summary["CCX_seconds"]+=result["native_seconds"]
            remaining=c["FEM_budget"]["total_CCX_budget_seconds"]-summary["CCX_seconds"]
            if result["status"]!="PASS":
                summary.update(hard_stop=True,stage_C="PARTIAL",overall="NUMERICAL_PARTIAL",
                    stop_reason=result.get("native_failure",result.get("recovery_failure","Native job did not pass actual output gates")))
        except Exception as error:
            summary["FEM_cases"][name].update(status="FAIL",failure=str(error))
            summary["attempts"][-1].update(status="FAIL",failure=str(error))
            summary.update(hard_stop=True,stage_C="PARTIAL",overall="NUMERICAL_PARTIAL",stop_reason=str(error))
        finally:
            summary["calls"]["CCX_production"]=sum(read(p)["solver_calls"] for p in (bundle/"FEM").glob("*/native_attempt.json"))
            save(bundle,item,summary)
        if summary["hard_stop"]:return summary
        if nonlinear:
            comparison.compare_fem_pair(bundle/"FEM"/(level+"_linear"),bundle/"FEM"/name,
                bundle/"one_d/joint_p64",discs[64],bundle/("comparison_"+level),
                T1=2*np.pi/c["omega1"],stage_b=stage_b)
            summary["stage_C"]="PARTIAL";save(bundle,item,summary)
    summary.update(completed=True,stage_C="COMPLETE_WITH_QUALIFICATIONS",overall="NUMERICAL_PARTIAL",
        limited_spatial_bending_comparison="COMPLETE_WITH_QUALIFICATIONS",energy_diagnostics="PARTIAL")
    save(bundle,item,summary);return summary


def write_output_tables(bundle,summary):
    """Export compact tables and compare already integrated kinetic energies."""
    energy={}
    for name,case in summary["FEM_cases"].items():
        if case["status"]!="PASS":continue
        path=bundle/"FEM"/name;recovery=read(path/"recovery.json")
        native={row["increment"]:row for row in read(path/"energy.json")["records"]
            if row["step"]==2 and "kinetic_energy" in row}
        with np.load(path/"independent_kinetic_energy.npz") as z:
            times,independent=z["time"].copy(),z["kinetic_energy"].copy()
        with np.load(path/"section_history.npz") as z:
            increments=z["increments"].copy()
            if not np.array_equal(times,z["time"]):raise ValueError("Saved kinetic/section timestamps differ")
        missing=[int(i) for i in increments if int(i) not in native]
        if missing:
            energy[name]={"energy_status":"PARTIAL","missing_native_kinetic_increment_ids":missing,
                "independent_internal_energy":"NOT_RUN","native_reference_jump_not_corrected":True}
            continue
        observed=np.array([native[int(i)]["kinetic_energy"] for i in increments])
        difference=independent-observed;scale=max(float(abs(independent).max()),float(abs(observed).max()))
        initial=recovery.get("native_initial_internal_energy")
        energy[name]={"max_absolute_kinetic_difference":float(abs(difference).max()),
            "common_full_horizon_kinetic_scale":scale,
            "relative_kinetic_difference_on_common_scale":float(abs(difference).max()/scale) if scale>0 else None,
            "relative_difference_to_own_STATIC_internal_energy":float(abs(difference).max()/initial) if initial is not None and initial>0 else None,
            "source":"already saved positive C3D10 reference-volume quadrature of actual FRD nodal velocities; matched native DAT increment IDs",
            "independent_internal_energy":"NOT_RUN","native_reference_jump_not_corrected":True,"energy_status":"PARTIAL"}
        with (path/"kinetic_comparison.csv").open("w",newline="",encoding="utf8") as stream:
            writer=csv.writer(stream);writer.writerow(("increment","physical_time","native_K","integrated_K","integrated_minus_native_K"))
            writer.writerows(zip(increments,times,observed,independent,difference))
    write(bundle/"energy_diagnostics.json",energy)
    for level in ("medium","fine"):
        path=bundle/("comparison_"+level)
        if not (path/"one_d_three_d_comparison.json").exists():continue
        data=read(path/"one_d_three_d_comparison.json")
        with (path/"comparison_metrics.csv").open("w",newline="",encoding="utf8") as stream:
            writer=csv.writer(stream);writer.writerow(("quantity","field","absolute_max","max_time_L2","relative_max_common","relative_L2_common","qualification"))
            for quantity,rows in data["metrics"].items():
                for field,row in rows.items():writer.writerow((quantity,field,*[row[key] for key in ("absolute_max","max_time_L2","relative_max_common","relative_L2_common")],row.get("qualification","")))
        with np.load(path/"one_d_three_d_comparison.npz") as z:
            mid=int(np.argmin(abs(z["x"]-.5)))
            with (path/"midspan_motion.csv").open("w",newline="",encoding="utf8") as stream:
                writer=csv.writer(stream);keys=("one_d_linear","one_d_nonlinear","three_d_linear","three_d_nonlinear","one_d_correction","three_d_correction","one_d_evolution","three_d_evolution")
                writer.writerow(("physical_time",*[key+"_"+field for key in keys for field in ("w","v")]))
                writer.writerows((t,*[z[key][i,mid,j] for key in keys for j in (1,2)]) for i,t in enumerate(z["time"]))
    return energy


def postprocess_saved(bundle,item,summary):
    """Finalize saved comparisons; this route contains no scientific solves."""
    if summary.get("postprocessing_complete"):return summary
    if not (summary["completed"] or summary["hard_stop"]):
        raise ValueError("Postprocess the finished/accepted prefix after native execution stops")
    from scripts.lib import nlsp_spatial_comparison as comparison
    phase=bundle/"execution_code/final_postprocessing";phase.mkdir(parents=True,exist_ok=True)
    for name in ("nlsp_spatial_nonlinear_verification","nlsp_spatial_comparison"):
        path=ROOT/"scripts/lib"/(name+".py");shutil.copyfile(path,phase/path.name)
    discs,_,_=load_discretizations(item["config"])
    for level in ("medium","fine"):
        if not all(summary["FEM_cases"].get(level+suffix,{}).get("status")=="PASS" for suffix in ("_linear","_nonlinear")):continue
        output=bundle/("comparison_"+level)
        previous=output/"first_processing_evidence"
        if output.exists() and not previous.exists():
            previous.mkdir()
            for path in output.glob("*"):
                if path.is_file():shutil.copyfile(path,previous/path.name)
        comparison.compare_fem_pair(bundle/"FEM"/(level+"_linear"),bundle/"FEM"/(level+"_nonlinear"),
            bundle/"one_d/joint_p64",discs[64],output,T1=2*np.pi/item["config"]["omega1"],
            stage_b=read(bundle/"stage_b_comparison.json"))
        # This audit uses only this completed pair's saved STATIC/final nodal
        # fields. Existing endpoint diagnostics replay from their own cache;
        # incomplete fine cases never enter this branch.
        from scripts.lib.nlsp_spatial_recovery_endpoints import audit_recovery_endpoints
        endpoint = audit_recovery_endpoints(bundle,level)
        endpoint_path = bundle/"recovery_sensitivity_endpoints"/level
        summary.setdefault("recovery_endpoint_sensitivity",{})[level] = {
            "status":endpoint["status"],
            "metrics":"recovery_sensitivity_endpoints/"+level+"/metrics.json",
            "manifest_sha256":sha(endpoint_path/"manifest.json") if (endpoint_path/"manifest.json").exists() else None,
            "recorded_81_section_recoveries":endpoint.get("recovery_calls_81",0),
            "scientific_calls":endpoint.get("scientific_calls",{}),
            "coverage":"confirmed STATIC anchor and actual final DYNAMIC frame only"}
    comparison.analyze_mesh_comparison(bundle)
    if not (bundle/"initial_state_decomposition.json").exists():
        cases={name:bundle/"one_d"/name for name in summary["one_d_cases"]
            if summary["one_d_cases"][name]["status"]=="PASS"}
        comparison.write_initial_state_decomposition(cases,bundle,discs=discs,
            T1=2*np.pi/item["config"]["omega1"])
    cases={name:bundle/"one_d"/name for name in summary["one_d_cases"]
        if summary["one_d_cases"][name]["status"]=="PASS"}
    comparison.audit_static_torsional_supports(bundle,cases,discs=discs,
        action_result=ROOT/item["config"]["source_bundles"]["action"]["path"]/"result.json")
    write_output_tables(bundle,summary)
    summary["figures"]=comparison.render_figures(bundle)
    summary["postprocessing_complete"]=True
    summary["postprocessing_scientific_calls"]={"CCX":0,"Gmsh":0,"Radau":0,"static_Newton":0,"eigen":0,"BVP":0}
    save(bundle,item,summary);return summary


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    mode=parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--preflight",action="store_true");mode.add_argument("--compute",action="store_true")
    mode.add_argument("--report-only",type=Path);mode.add_argument("--plot-only",type=Path)
    mode.add_argument("--postprocess-only",type=Path)
    parser.add_argument("--config",type=Path,default=CONFIG)
    parser.add_argument("--through-stage",choices=("A","B","C"),default="C")
    parser.add_argument("--additional-controls",action="store_true",
        help="Bounded p48 repeats of the two already authorized isolated excitations; no new physical scenario")
    args=parser.parse_args(argv)
    if args.report_only or args.plot_only or args.postprocess_only:
        bundle=(args.report_only or args.plot_only or args.postprocess_only).resolve();summary=validate_cache(bundle)
        if args.postprocess_only:
            summary=postprocess_saved(bundle,read(bundle/"provenance.json"),summary)
        if args.plot_only:
            # Matching plot replay leaves the completed manifest and images
            # unchanged. First rendering is part of finalizing this bundle.
            if not (bundle/"figure_data_provenance.json").exists():
                from scripts.lib.nlsp_spatial_comparison import render_figures
                render_figures(bundle);save(bundle,read(bundle/"provenance.json"),summary)
    else:
        bundle,item,summary=prepare(args.config)
        if args.compute:
            summary=compute_stage_ab(bundle,item,summary,through_stage=args.through_stage)
            if args.additional_controls or args.through_stage=="C":
                summary=compute_additional_controls(bundle,item,summary)
            if args.through_stage=="C":
                summary=compute_stage_c(bundle,item,summary)
    print(json.dumps({"bundle":str(bundle),"overall":summary["overall"],"stage_A":summary["stage_A"],
        "stage_B":summary["stage_B"],"stage_C":summary["stage_C"],"calls":summary["calls"]},indent=2))
    return summary
