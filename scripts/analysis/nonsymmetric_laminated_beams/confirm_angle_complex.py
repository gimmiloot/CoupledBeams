"""ANGLE-C fixed preset: reuse NSL-3A corrector for two roots, then three figures.

Script Proliferation Control: parameterize the EXISTING NSL provider and
evaluate_target, without copying Newton/physics. This diagnostic entry has a
new validated ANGLE input plus a three-source figure-table contract. Editing
the historical runners would invalidate their immutable source hashes.
The separate fixed renderer has no spectral code and is not a general app.
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
import time

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
import numpy as np
import scipy
from scripts.analysis.nonsymmetric_laminated_beams import screen_elastic_angles as angle
from scripts.analysis.nonsymmetric_laminated_beams import verify_kelvin_voigt_pilot as pilot
from scripts.analysis.check_contraction_joint_sensitivity import plain
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import nonsymmetric_laminated_joint as elastic
from scripts.lib import nsl_angle_research_figures as figures

CONFIG=ROOT/"data/input/nonsymmetric_laminated_beam_angle_complex_confirmation.json"
OUTPUT=ROOT/"results/nonsymmetric_laminated_beams/nsl3a_angle_complex_confirmation"
read=angle.read


def load_context():
    config=read(CONFIG)
    old,source,props,ref,manifests,kv=angle.context(read(angle.CONFIG))
    fp,ident=angle.identity(read(angle.CONFIG),manifests)
    discovered=list((ROOT/"results/nonsymmetric_laminated_beams").rglob("manifest.json"))
    paths=[p for p in discovered if p.parent.name==fp]
    assert len(paths)==1 and fp==config["source_bundles"][-1]
    path=paths[0];screen=angle.validate_cache(path.parent,ident)
    assert screen["status"]=="PASS"
    manifests[fp]=dict(path=path.relative_to(ROOT).as_posix(),sha256=sha(path))
    assert set(manifests)==set(config["source_bundles"])
    assert config["beta_deg"]==15 and config["kappa_theta"]==1 and config["d_theta"]==.001
    assert config["variant"]==old["variant"]=="PROJECT_REDDY_COUPLED"
    assert ref==screen["reference"]
    assert config["states"]==[
        dict(id="C1",case="MIRROR",position=5,role="ACTIVE"),
        dict(id="C2",case="SAME_LOCAL",position=5,role="ACTIVE")]
    # All numerical tolerances and universal solver functions are inherited.
    policy={**read(ROOT/config["policy_from"])["policy"],**config["policy_overrides"]}
    config={**config,"policy":policy}
    pair=next(r for r in screen["identity_diagnostics"] if
        r["kind"]=="SAME_ANGLE_CROSS_LAYOUT" and r["reference"]==["MIRROR",15])
    match=pair["assignments"][4]
    assert match["reference_position"]==match["target_position"]==5 and match["status"]=="RESOLVED"
    states=[]
    for specification in config["states"]:
        point=next(p for p in screen["points"] if p["case"]==specification["case"] and p["beta_deg"]==15)
        mode=point["modes"][4]
        assert point["status"]==mode["status"]=="PASS" and not mode["structural_zero"]
        assert mode["alpha_pred"]==ref["c_theta"]*mode["delta_psi"]**2/(2*mode["mass"])
        states.append(dict(**specification,omega=mode["omega"],saved=mode,
            count=point["count"]["roots"][4],
            gap=min(abs(m["omega"]-mode["omega"]) for m in point["modes"] if m["position"]!=5),
            source_manifest=manifests[fp],
            source_pointer=f'/points/{screen["points"].index(point)}/modes/4',
            accepted_cross_layout_pair=match))
    return config,old,source,props,ref,manifests,kv,screen,states


def preservation(manifests):
    path=OUTPUT/"baseline.json"
    if not path.exists():
        hashes=dict(read(angle.OUTPUT/"baseline.json")["sha256"])
        for item in manifests.values():
            for p in (ROOT/item["path"]).parent.rglob("*"):
                if p.is_file():hashes[p.relative_to(ROOT).as_posix()]=sha(p)
        for p in (angle.CONFIG,Path(angle.__file__),
                  ROOT/"tests/test_nsl3a_elastic_angle_screening.py",
                  ROOT/"docs/nonsymmetric_laminated_beams/nsl3a_elastic_angle_screening.md",
                  angle.OUTPUT/"README.md",angle.OUTPUT/"baseline.json",
                  angle.OUTPUT/"current.json",angle.OUTPUT/"verification.json"):
            hashes[p.relative_to(ROOT).as_posix()]=sha(p)
        memory={}
        for name in ("decisions.md","knowledge.md"):
            p=ROOT/"docs/memory"/name
            memory[p.relative_to(ROOT).as_posix()]=dict(size=p.stat().st_size,sha256=sha(p))
        write_json(path,dict(sha256=hashes,memory_prefix=memory))
    b=read(path)
    assert all(sha(ROOT/n)==h for n,h in b["sha256"].items())
    for n,v in b["memory_prefix"].items():
        assert hashlib.sha256((ROOT/n).read_bytes()[:v["size"]]).hexdigest()==v["sha256"]
    return b


def identity(config,manifests):
    preservation(manifests)
    paths=[CONFIG,Path(__file__),Path(pilot.__file__),Path(pilot.model.__file__)]
    ident=dict(config=config,source_manifests=manifests,
        baseline_sha256=sha(OUTPUT/"baseline.json"),
        files={p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16],ident


def prepare_saved_tables(config,old,source,screen,kv,manifests,out):
    """CSV/profile processing only; no matrix assembly, SVD or elastic solve."""
    maxima=[]
    for beta in (15,30,45,60,75):
        for case in ("MIRROR","SAME_LOCAL"):
            selected=[r for r in screen["table"] if r["case"]==case and r["beta_deg"]==beta]
            assert len(selected)==6
            row=max(selected,key=lambda r:r["zeta_pred"])
            maxima.append(dict(beta_deg=beta,case=case,position=row["position"],
                zeta_pred_max=row["zeta_pred"],alpha_pred=row["alpha_pred"],omega0=row["omega"]))
    write_csv(out/"figure1_maxima.csv",maxima)
    profiles_path=(ROOT/manifests["61842957b5e4d8d1"]["path"]).parent/"profiles.csv"
    with profiles_path.open(encoding="utf-8",newline="") as f:
        profiles=list(csv.DictReader(f))
    n=old["policy"]["mode_quadrature"]
    x,w=np.polynomial.legendre.leggauss(n)
    nodes=(x+1)/2
    pair_records=[];rotations=[]
    for pos in (4,6):
        accepted=source["cross_case"]["MIRROR_SPRING__SAME_LOCAL_SPRING"]["assignments"][pos-1]
        assert accepted["identity"]=="RESOLVED" and accepted["right_position"]==pos
        fields={};masses={};scales={}
        # Here all arms have equal m,j2 and j1=0: the common physical pairing
        # at the SAME beta equals this local pairing. Verify these assumptions.
        props={name:source["materials"][name]["properties"] for name in ("HHLL","LLHH")}
        for name in ("HHLL","LLHH"):
            assert props[name]["j1"]==0
            assert abs(props[name]["m"]/props["HHLL"]["m"]-1)<1e-14
            assert abs(props[name]["j2"]/props["HHLL"]["j2"]-1)<1e-14
        m,j2=props["HHLL"]["m"],props["HHLL"]["j2"]
        for case in ("MIRROR","SAME_LOCAL"):
            arms=[]
            for arm in (1,2):
                rows=[r for r in profiles if r["case"]==case+"_SPRING" and
                      int(r["position"])==pos and int(r["arm"])==arm]
                assert len(rows)==n+2
                xx=np.array([float(r["x"]) for r in rows]);L=old["lengths"][arm-1]
                assert xx[0]==0 and xx[-1]==L
                assert np.allclose(xx[1:-1]/L,nodes,rtol=0,atol=1e-15)
                q=np.array([[float(r[k]) for k in ("u","w","psi")] for r in rows])
                arms.append((xx/L,q,w*L/2))
            mass=sum(weight@(m*np.sum(q[1:-1,:2]**2,axis=1)+j2*q[1:-1,2]**2)
                     for _,q,weight in arms)
            saved=source["spectral"][case+"_SPRING"]["modes"][pos-1]
            assert abs(mass-saved["mass"])<config["policy"]["mass_tol"]
            masses[case]=float(mass);scales[case]=1/np.sqrt(mass);fields[case]=arms
        overlap=sum(weight@(m*np.sum(a[1:-1,:2]*b[1:-1,:2],axis=1)+j2*a[1:-1,2]*b[1:-1,2])
                    for (_,a,weight),(_,b,_) in zip(fields["MIRROR"],fields["SAME_LOCAL"]))
        mac=float(overlap**2/(masses["MIRROR"]*masses["SAME_LOCAL"]))
        assert abs(mac-accepted["MAC"])<1e-10
        # One common display orientation, then the SIGNED mass pairing.
        qm=fields["MIRROR"][0][1][:,2]
        sign_m=float(np.sign(qm[np.argmax(abs(qm))]))
        signs={"MIRROR":sign_m,"SAME_LOCAL":sign_m*(1. if overlap>=0 else -1.)}
        for case in ("MIRROR","SAME_LOCAL"):
            factor=signs[case]*scales[case]
            state=next(r for r in kv["states"] if r["case"]==case and r["elastic_position"]==pos)
            for arm,(xi,q,weight) in enumerate(fields[case],1):
                for xx,qq in zip(xi,q):
                    rotations.append(dict(position=pos,case=case,arm=arm,xi=float(xx),
                        psi=float(factor*qq[2]),source_psi=float(qq[2]),display_factor=float(factor)))
            delta=factor*(fields[case][0][1][-1,2]-fields[case][1][1][-1,2])
            pair_records.append(dict(position=pos,case=case,mass_original=masses[case],mass_display=1.,
                display_factor=float(factor),signed_mass_overlap=float(overlap),pair_MAC=mac,
                structural_zero=case=="MIRROR",delta_plot=float(delta),
                actual_zeta=state["zeta"],complex_source_id=state["id"]))
    write_csv(out/"figure2_rotations.csv",rotations)
    write_csv(out/"figure2_pairs.csv",pair_records)
    return dict(status="PASS",figure1_rows=len(maxima),figure2_rows=len(rotations),
        pairs=pair_records,profile_source=dict(path=profiles_path.relative_to(ROOT).as_posix(),sha256=sha(profiles_path)),
        heatmap=dict(path=(ROOT/manifests[config["source_bundles"][-1]]["path"]).parent
                         .joinpath("angle_screening.png").relative_to(ROOT).as_posix(),
                     reused_without_rendering=True))


def beta15_power(ref):
    maps,jump=elastic.trace_maps(np.deg2rad(15),elastic.Joint("SPRING",ref["k_theta"]))
    q=np.array([.2,-.3,.4,-.1]);vel=np.array([-.1,.5,.7,-.2])
    moment=ref["k_theta"]*(jump@q)+ref["c_theta"]*(jump@vel)
    force=np.array([.3,-.7])
    forces=[np.r_[e[:2,:2]@force*sign,-moment*sign] for e,sign in zip(maps,(1,-1))]
    power=sum(f@(e@vel) for f,e in zip(forces,maps))+ref["k_theta"]*(jump@q)*(jump@vel)
    error=abs(power+ref["c_theta"]*(jump@vel)**2)
    return dict(beta_deg=15,power=float(power),residual=float(error),
                status="PASS" if power<=0 and error<1e-14 else "FAIL")


def comparison(states):
    a,b=states
    result={}
    for metric,pred in (("alpha","alpha_pred"),("zeta","zeta_pred")):
        ratio=b[metric]/a[metric];prediction=b[pred]/a[pred]
        result[metric]=dict(MIRROR=a[metric],SAME_LOCAL=b[metric],
            ratio=ratio,absolute_difference=b[metric]-a[metric],
            predicted_ratio=prediction,predicted_absolute_difference=b[pred]-a[pred],
            ratio_relative_difference=ratio/prediction-1,
            contrast_percent=100*(ratio-1),predicted_contrast_percent=100*(prediction-1))
    return result


def envelope_tables(states,config,out):
    modes=[];rows=[]
    for r in states:
        assert r["status"]=="PASS"
        ne=r["omega_d"]/(2*np.pi*r["alpha"])
        modes.append(dict(id=r["id"],case=r["case"],alpha=r["alpha"],omega_d=r["omega_d"],
            zeta=r["zeta"],N_1e=ne,amplitude_one_period=np.exp(-1/ne),
            amplitude_N100=np.exp(-100/ne),amplitude_N200=np.exp(-200/ne)))
        for n in np.linspace(0,config["figures"]["envelope_max_periods"],config["figures"]["envelope_samples"]):
            rows.append(dict(case=r["case"],N=float(n),amplitude=float(np.exp(-n/ne)),
                time=float(2*np.pi*n/r["omega_d"])))
    write_csv(out/"figure3_modes.csv",modes);write_csv(out/"figure3_envelopes.csv",rows)
    return plain(modes)


def validate_cache(out,ident=None):
    manifest=read(out/"manifest.json")
    if ident is not None:assert manifest["identity"]==ident
    assert all(sha(out/n)==h for n,h in manifest["artifact_hashes"].items())
    return read(out/"result.json")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--plot-only",action="store_true")
    args=parser.parse_args()
    if args.plot_only:
        out=OUTPUT/read(OUTPUT/"current.json")["fingerprint"]
        result=validate_cache(out)
        assert result["status"]=="PASS"
        destination=OUTPUT/"plot_only"
        figure=figures.render(out,destination,read(CONFIG)["figures"]["dpi"])
        print(json.dumps(dict(status=figure["status"],mode="plot_only",figures=3,new_eigenvalue_solves=0)))
        return
    config,old,source,props,ref,manifests,kv,screen,states=load_context()
    fp,ident=identity(config,manifests);out=OUTPUT/fp
    if (out/"manifest.json").exists():
        result=validate_cache(out,ident)
        if result["status"]=="PASS":
            f=read(out/"figures/figure_manifest.json")
            assert all(sha(out/"figures"/n)==h for n,h in f["artifact_hashes"].items())
        print(json.dumps(dict(status=result["status"],fingerprint=fp,cache="VALIDATED",new_complex_solves=0)))
        return
    out.mkdir(parents=True,exist_ok=True)
    prepared=prepare_saved_tables(config,old,source,screen,kv,manifests,out)
    checks=pilot.structural(config,old,props,ref,states)
    checks["beta15_power"]=beta15_power(ref)
    write_json(out/"preflight_contract.json",plain(dict(identity=ident,reference=ref,
        states=states,structural=checks,figure_sources=prepared)))
    assert checks["status"]==checks["beta15_power"]["status"]=="PASS"
    if not args.compute:
        print(json.dumps(dict(status="PASS",stage="PREFLIGHT",fingerprint=fp,new_complex_solves=0)))
        return
    if (out/"ledger.json").exists():
        raise RuntimeError("Prior corrector attempt exists; no automatic replay or extra solve")
    ledger=dict(primary_complex_calls=0,auxiliary_complex_calls=0,newton_updates=0,
                protected_newton_calls=0,elastic_root_solves=0,prefix_searches=0)
    result=dict(status="IN_PROGRESS",states=[],reference=ref,structural=checks,
                figure_sources=prepared,solve_ledger=ledger)
    started=time.perf_counter()
    try:
        for state in states:
            assert ledger["primary_complex_calls"]<config["max_complex_solves"]
            row=pilot.evaluate_target(config,old,props,ref,state,out,ledger)
            row.update(zeta_pred=row["alpha_pred"]/row["omega0"],
                zeta_predictor_relative_difference=row["zeta"]/(row["alpha_pred"]/row["omega0"])-1)
            result["states"].append(row)
            write_json(out/(state["id"]+".json"),plain(row))
            write_json(out/"progress.json",plain(result))
            print(state["id"],row["status"],row["p"],"steps",row["correction"]["steps"],flush=True)
            if row["status"]!="PASS":
                raise ArithmeticError(f'{state["id"]}: physical/numerical gate failed')
        assert ledger["primary_complex_calls"]==2
        result.update(status="PASS",comparison=comparison(result["states"]),
                      envelopes=envelope_tables(result["states"],config,out))
        write_csv(out/"comparison.csv",[
            {n:r[n] for n in ("id","case","omega0","omega_d","alpha_pred","alpha","zeta_pred","zeta",
                              "predictor_relative_difference","zeta_predictor_relative_difference","MAC",
                              "frequency_shift","alpha_energy_absolute_error","alpha_energy_relative_error","status")}
            | {"p_real":r["p"].real,"p_imag":r["p"].imag,"newton_updates":r["correction"]["steps"]}
            for r in result["states"]])
    except Exception as error:
        import traceback
        result.update(status="PARTIAL_PASS",stop_reason=f"{type(error).__name__}: {error}",traceback=traceback.format_exc())
    result["spectrum_seconds"]=time.perf_counter()-started
    result["protected_files_verified"]=len(preservation(manifests)["sha256"])
    write_json(out/"result.json",plain(result))
    write_json(out/"manifest.json",dict(identity=ident,
        git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        command=subprocess.list2cmdline([sys.executable,*sys.argv]),
        artifact_hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}))
    write_json(OUTPUT/"current.json",dict(fingerprint=fp,status=result["status"]))
    if result["status"]=="PASS":
        figures.render(out,out/"figures",config["figures"]["dpi"])
    print(json.dumps(dict(status=result["status"],fingerprint=fp,ledger=ledger,
                          stop_reason=result.get("stop_reason"))))


if __name__=="__main__":
    main()
