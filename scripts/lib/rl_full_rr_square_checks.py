"""Bounded aspect gate for the existing map CLI; no separate runnable solver."""
from dataclasses import asdict
from pathlib import Path
import csv
import platform
import subprocess
import time
import numpy as np
import scipy
from scripts.analysis import map_rl_full_rr as m
from scripts.analysis import verify_rl_full_rr_prefix as ab
from scripts.analysis import verify_rl_full_rr_geometry as cg
from scripts.lib import rl_full_rr_sections_prefix as e
from scripts.lib.rl_full_rr_aspect import make_model


def check_coefficients(model, case):
    """Direct rectangular identities, including full polar-gradient inertia."""
    checks = []
    for arm in (make_model(**{**asdict(model), "mu":0., "tau":0.}), *model.arms):
        h = np.sqrt(12)*arm.l*arm.epsilon; b = arm.aspect*h
        expected = dict(A=b*h, I=b*h**3/12, I_P=b*h*(b*b+h*h)/12)
        expected["J"] = arm.rho*arm.nu**2*expected["I_P"]
        actual = dict(A=arm.section.rhoA/arm.rho, I=arm.section.EI/arm.E,
                      I_P=(1+arm.aspect**2)*arm.section.EI/arm.E, J=arm.J)
        checks.append(dict(expected=expected, actual=actual,
                           error=max(abs(actual[k]/v-1) for k,v in expected.items())))
    geometry = m.geometry(model,case)
    expected_time = model.l/model.epsilon*np.sqrt(model.rho/model.E)
    ok = max(x["error"] for x in checks)<1e-13 and abs(model.time-expected_time)<1e-12
    return dict(status="PASS" if ok else "SECTION_CONTRACT_FAILED", sections=checks, geometry=geometry)


def control_gate(c, policy):
    known, refs, sources = m.accepted_sources(c)
    for name in ["scripts/analysis/map_rl_full_rr.py", "scripts/lib/rl_full_rr_aspect.py", __file__]:
        p=Path(name).resolve(); sources[p.relative_to(m.ROOT).as_posix()]=m.sha(p)
    aspect=c["material"]["aspect"]
    identity=dict(schema="rl-full-rr-aspect-gate-v1", material=c["material"],policy=policy,
        sources=sources,references=refs,initial_upper_omega=c["initial_upper_omega"],
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=m.ROOT,text=True).strip(),
        scope="four a4 regressions; four aspect comparisons; reflected pair at beta0,30; Ritz in four new aspect controls")
    folder=m.ROOT/"results/rl_full_rr_square"/m.digest(identity)[:16]
    if (folder/"manifest.json").exists():
        result=m.verify_manifest(folder)
        return dict(result,folder=folder.relative_to(m.ROOT).as_posix(),manifest_sha256=m.sha(folder/"manifest.json"))
    folder.mkdir(parents=True,exist_ok=True)
    if not (folder/"identity.json").exists():m.publish(folder/"identity.json",identity)
    # A changed control/report helper need not repeat already accepted solves.
    # Re-run ALL checks below, reusing only same-contract points whose complete
    # solver dependency hashes (excluding this non-solver helper) are unchanged.
    reusable={}
    helper=Path(__file__).relative_to(m.ROOT).as_posix()
    numerical_sources={k:v for k,v in sources.items() if k!=helper}
    for old_identity in folder.parent.glob("*/identity.json"):
        previous=m.read(old_identity)
        if (previous.get("policy")!=policy or previous.get("versions")!=identity["versions"] or
            not m.solver_sources_compatible(previous.get("sources",{}),numerical_sources)):continue
        for point in (old_identity.parent/"points").glob("*.json"):
            saved=m.verify_record(point)
            if saved["status"]=="PASS":reusable.setdefault(m.physical_key(saved["contract"]),(point,saved))
    key=m.digest(identity); points=[]; checks=[]; data={}; status="PASS"; start=time.perf_counter()
    # Fixed finite controls: do not expand into section-shape or symmetry maps.
    cases=[(4.,.1,0.,b) for b in [0.,5.,30.,90.]]+[(4.,.1,-.1,30.)]
    cases += [(aspect,.1,0.,b) for b in [0.,5.,30.,90.]]
    cases += [(aspect,mu,tau,b) for mu,tau in [(.1,-.1),(-.1,.1)] for b in [0.,30.]]
    try:
        for a,mu,tau,beta in dict.fromkeys(cases):
            model=make_model(**dict(c["material"],aspect=a),mu=mu,tau=tau)
            case=e.p.Case(beta,mu,"RIGID","RIGID",tau)
            contract=e.point_contract(model,case,8,9,policy)
            if m.physical_key(contract) in reusable:
                path,rec=reusable[m.physical_key(contract)];calculated=False
            else:
                path,rec,calculated=m.solve(contract,folder,key,c["initial_upper_omega"])
            item=dict(aspect=a,mu=mu,tau=tau,beta_deg=beta,record=path.relative_to(m.ROOT).as_posix(),sha256=m.sha(path))
            points.append(item)
            check=dict(**item,status=rec["status"],coefficients=check_coefficients(model,case))
            if rec["status"]!="PASS":
                checks.append(check);status="STOP_POINT";break
            modes=ab.load_modes(path.parent.parent,rec)
            for mode,root in zip(modes,rec["roots"]):mode["coefficients"]=root["coefficients"]
            data[a,mu,tau,beta]=(rec,modes,model,case)
            if a==4 and m.physical_key(contract) in known:
                oldpath,old,_=known[m.physical_key(contract)]
                check["regression"]=e.compare_modes(model,rec["roots"],modes,old["roots"],ab.load_modes(oldpath.parent.parent,old),policy)
                check["regression"]["reference"]=oldpath.relative_to(m.ROOT).as_posix()
            if a==aspect and beta in (0.,30.):
                # Four Ritz controls: uniform straight, uniform angled, and two reflected angled sections.
                if tau==0 or beta==30.:
                    check["ritz"]=cg.ritz_check(model,case,rec["roots"],[],policy)
                check["quadrature"]=e.quadrature_check(model,case,rec["roots"],modes,policy)
                if beta==0:
                    check["straight_interface"]=e.straight_interface(model,case,rec["roots"],modes,policy)
                    if tau==0:
                        check["straight_uniform"]=ab.straight_split(model.arms[0],case,rec["roots"],modes,c["initial_upper_omega"],policy)
            if any(v.get("status")!="PASS" for v in check.values() if isinstance(v,dict) and "status" in v):
                check["status"]="CONTROL_FAILED"
            checks.append(check)
            print(f"aspect gate a={a:g} mu={mu:g} tau={tau:g} beta={beta:g}: {check['status']}",flush=True)
            if check["status"]!="PASS":status="STOP_CONTROL";break
        reflections=[]
        if status=="PASS":
            for beta in (0.,30.):
                ra,ma,_,_=data[aspect,.1,-.1,beta]
                rb,mb,model,_=data[aspect,-.1,.1,beta]
                compare=e.compare_modes(model,rb["roots"],mb,ra["roots"],[e.reflect_mode(x) for x in ma],policy)
                reflections.append(dict(compare,beta_deg=beta,arm_permutation=[2,1],state_transform=e.REFLECTION_STATE.tolist()))
                if compare["status"]!="PASS":status="STOP_REFLECTION"
        comparisons=[]
        if status=="PASS":
            for beta in (0.,5.,30.,90.):
                ra=data[aspect,.1,0.,beta][0]; r4=data[4.,.1,0.,beta][0]
                for i,(x,y) in enumerate(zip(ra["roots"],r4["roots"]),1):
                    comparisons.append(dict(beta_deg=beta,rank=i,aspect=aspect,omega_a=x["omega"],omega_4=y["omega"],
                        Lambda_a=x["Lambda"],Lambda_4=y["Lambda"],omega_relative=x["omega"]/y["omega"]-1,
                        Lambda_relative=x["Lambda"]/y["Lambda"]-1))
            with (folder/"aspect_comparison.csv").open("x",newline="",encoding="utf-8") as f:
                w=csv.DictWriter(f,fieldnames=list(comparisons[0]));w.writeheader();w.writerows(comparisons)
        result=dict(status=status,identity=identity,points=points,checks=checks,reflections=reflections,
            elapsed_seconds=time.perf_counter()-start,
            artifacts={p.relative_to(folder).as_posix():m.sha(p) for p in sorted(folder.rglob("*")) if p.is_file()})
    except (ArithmeticError,ValueError,RuntimeError,KeyError,np.linalg.LinAlgError) as ex:
        result=dict(status="STOP_CONTROL_EXCEPTION",identity=identity,points=points,checks=checks,
                    failure=dict(type=type(ex).__name__,reason=str(ex)),
                    artifacts={p.relative_to(folder).as_posix():m.sha(p) for p in sorted(folder.rglob("*")) if p.is_file()})
    result=ab.jsonable(result)
    m.publish(folder/"manifest.json",result)
    return dict(result,folder=folder.relative_to(m.ROOT).as_posix(),manifest_sha256=m.sha(folder/"manifest.json"))
