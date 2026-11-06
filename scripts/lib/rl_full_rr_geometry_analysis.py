"""Saved C0/C1 diagnostics. Run: python -c "from scripts.lib.rl_full_rr_geometry_analysis import run; run()"

Read-only physical postprocessing, separated to preserve the sealed calculation
source fingerprint. No new solver, no parameter points. Imports accepted energy,
global mass and sample-node functions; writes a distinct derived result family.
"""
from pathlib import Path
import csv
import json
from unittest.mock import patch
from contextlib import ExitStack
import numpy as np
from scripts.analysis import verify_rl_full_rr_geometry as calc
from scripts.analysis import postprocess_rl_full_rr_modes as saved

ROOT=calc.ROOT
SOURCE=ROOT/"results/rl_full_rr_geometry/a79da9b64ef44143"
OUTPUT=ROOT/"results/rl_full_rr_geometry_analysis"


def csv_new(path,rows):
    with path.open("x",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def analyze(manifest):
    c=manifest["identity"]["config"];by={}
    for entry in manifest["C1_points"]:
        rec,modes=calc.resolve(entry)
        by[entry["mu"],entry["beta_deg"]]=(entry,rec,modes)
    rows,nodes,macs,matches,packets,gaps=[],[],[],[],[],[]
    for (mu,beta),(entry,rec,modes) in by.items():
        model=calc.p.Model(**c["material"],mu=mu,tau=0.)
        case=calc.p.Case(beta,mu)
        initial=by[mu,0.][1]
        for k,(root,mode) in enumerate(zip(rec["roots"],modes),1):
            mc,kc,sp=calc.p.energies(model,case,mode["arms"])
            mass,elastic=float(mc.sum()),float(kc.sum())
            assert abs(elastic/mass/root["omega"]**2-1)<=c["policy"]["energy_tol"]
            assert abs(mass-1)<=c["policy"]["mass_gram_tol"]
            if k>8:continue
            mf,kf=mc.sum(0)/mass,kc.sum(0)/elastic
            peak=max(np.linalg.norm(a["y"][:,:2],axis=1).max() for a in mode["arms"])
            row=dict(mu=mu,beta_deg=beta,rank=k,omega=root["omega"],Lambda=root["Lambda"],
                F_omega=root["omega"]/initial["roots"][k-1]["omega"],K=elastic,M=mass,
                K_axial=float(kf[0]),K_bending=float(kf[1]),K_shear=float(kf[2]),
                M_axial=float(mf[0]),M_transverse=float(mf[1]),M_rotation=float(mf[2]),M_RL=float(mf[3]),
                K_arm1=float(kc[0].sum()/elastic),M_arm1=float(mc[0].sum()/mass),
                joint_w_over_peak=float(abs(mode["arms"][0]["ends"][1,1])/peak),
                source_point=entry["record"],source_shape=str((ROOT/entry["record"]).parent.parent/rec["shape_file"]),
                origin=entry["origin"])
            rows.append(row)
            if beta==0:
                aa,bb=mode["arms"];L1,L2=case.lengths(model)
                s=np.r_[0.,aa["x"],L1,2*model.l-bb["x"][::-1],2*model.l]
                w=np.r_[-aa["ends"][0,1],-aa["y"][:,1],-aa["ends"][1,1],
                         bb["y"][::-1,1],bb["ends"][0,1]]
                node=saved.node_brackets(s,w,peak)
                brackets=node["brackets"]
                nearest=min(brackets,key=lambda z:abs((z[0]+z[1])/2-L1)) if brackets else None
                distance=([max(0.,nearest[0]-L1,L1-nearest[1]),max(abs(L1-nearest[0]),abs(L1-nearest[1]))]
                          if nearest else None)
                nodes.append(dict(mu=mu,rank=k,joint_s_over_2l=L1/(2*model.l),
                    joint_w_over_peak=row["joint_w_over_peak"],node_count=node["count"],
                    brackets_s=brackets,brackets_s_over_2l=[[v/(2*model.l) for v in z] for z in brackets],
                    nearest_node_distance_interval_over_2l=[v/(2*model.l) for v in distance] if distance else None,
                    amplitude_floor=saved.AMPLITUDE_FLOOR,source=entry["record"]))
        w=np.array([x["omega"] for x in rec["roots"]])
        for k in range(8):
            gaps.append(dict(mu=mu,beta_deg=beta,rank_low=k+1,rank_high=k+2,
                gap_omega=float(w[k+1]-w[k]),relative_gap=float((w[k+1]-w[k])/w[k]),
                includes_guard=k==7,source=entry["record"]))
    for mu in c["pilot_mu"]:
        model=calc.p.Model(**c["material"],mu=mu,tau=0.)
        for a,b in zip(c["pilot_angles"][:-1],c["pilot_angles"][1:]):
            _,ra,ma=by[mu,a];_,rb,mb=by[mu,b]
            na=np.array([calc.p.inner(model,x,x) for x in ma])
            nb=np.array([calc.p.inner(model,x,x) for x in mb])
            overlap=np.array([[calc.e.global_inner(model,x,y,a,b) for y in mb] for x in ma])/np.sqrt(na[:,None]*nb)
            mac=overlap**2
            assert mac.max()<=1+1e-10
            for i in range(9):
                for j in range(9):
                    macs.append(dict(mu=mu,beta_from=a,beta_to=b,rank_from=i+1,rank_to=j+1,
                                     overlap=float(overlap[i,j]),MAC=float(mac[i,j]),includes_guard=(i==8 or j==8)))
                if i<8:
                    ids=np.argsort(mac[i])[::-1];maximum=float(mac[i,ids[0]]);margin=maximum-float(mac[i,ids[1]])
                    ambiguous=maximum<c["policy"]["mac_min"] or margin<c["policy"]["assignment_margin_min"]
                    matches.append(dict(mu=mu,beta_from=a,beta_to=b,rank_from=i+1,best_rank=int(ids[0])+1,
                        best_MAC=maximum,margin=margin,diagonal_MAC=float(mac[i,i]),
                        status="UNRESOLVED_DISCRETE_MATCH" if ambiguous else "STRONG_DISCRETE_MATCH_NOT_BRANCH_ID"))
            for pair in [(1,2),(3,4),(5,6),(7,8)]:
                ids=np.array(pair)-1;sv=np.linalg.svd(overlap[np.ix_(ids,ids)],compute_uv=False)
                packets.append(dict(mu=mu,beta_from=a,beta_to=b,ranks=list(pair),singular_values=sv.tolist(),
                                    min_projection=float(sv[-1]**2)))
    ref=np.array([x["omega"] for x in by[.01,0.][1]["roots"]])
    invariant=[dict(mu=mu,relative_omega_errors=(np.array([x["omega"] for x in by[mu,0.][1]["roots"]])/ref-1).tolist())
               for mu in c["pilot_mu"]]
    symmetry=[]
    model=calc.p.Model(**c["material"],mu=0.,tau=0.)
    for beta in c["pilot_angles"]:
        entry,rec,modes=by[0.,beta]
        for k,mode in enumerate(modes,1):
            reflected=calc.e.reflect_mode(mode)
            parity=calc.p.inner(model,mode,reflected)/calc.p.inner(model,mode,mode)
            symmetry.append(dict(beta_deg=beta,rank=k,reflection_overlap=float(parity),
                parity=1 if parity>=0 else -1,error=float(abs(abs(parity)-1)),source=entry["record"]))
    return rows,nodes,macs,matches,packets,gaps,invariant,symmetry


def technical(source,manifest):
    records=[json.loads(x.read_text()) for x in sorted((source/"points").glob("*.json"))]
    modes=[r for x in records for r in x["roots"]]
    gate=json.loads((source/"C0_gate.json").read_text())
    summary=dict(new_points=len(records),new_modes=len(modes),C0_new=sum(x["stage"]=="C0" for x in records),
        C1_new=sum(x["stage"]=="C1" for x in records),pilot_points=len(manifest["C1_points"]),
        reused_AB=sum(x["origin"]=="accepted_AB" for x in manifest["C1_points"]),
        reused_C0=sum(x["origin"]=="verified_C0" for x in manifest["C1_points"]),
        max_boundary_error=float(max(max(np.max(r["diagnostics"]["joint_scaled"]),np.max(r["diagnostics"]["external_scaled"])) for r in modes)),
        max_ode_error=max(r["diagnostics"]["ode_error"] for r in modes),
        max_energy_error=max(r["diagnostics"]["energy_error"] for r in modes),
        max_gram_error=max(x["quality"]["gram_max_error"] for x in records),
        max_Ritz_error=max(t["relative_frequency_error"] for x in records for t in x["ritz"]["trials"] if t["accepted"]),
        max_quadrature_energy_change=max(r["energy_relative_change"] for x in records for r in x["quadrature"]["modes"]),
        max_reference_frequency_error=max(r["relative_frequency_error"] for x in records for r in x.get("regressions",[])),
        max_reference_field_error=max(g["max_field_error"] for x in records for r in x.get("regressions",[]) for g in r["groups"]),
        reflection=gate["reflections"],max_multiplicity=max(r["multiplicity"] for r in modes),
        count_queries=[len(x["search"]["samples"]) for x in records],
        min_guard_catalog_margin=min(x["domain"]["catalog_limit"]-x["roots"][8]["omega"] for x in records),
        max_guard_Tim_ratio=max(x["roots"][8]["omega"]/a["Tim"] for x in records for a in x["domain"]["arms"]),
        max_RL_ratio=max(x["roots"][8]["omega"]/a["RL"] for x in records for a in x["domain"]["arms"]),
        all_lower_counts_zero=all(x["search"]["lower_count"]==0 for x in records),
        all_guard_counts_nine=all(x["search"]["guard_right_count"]==9 for x in records))
    interfaces=[]
    for rec in records:
        if rec["beta_deg"]!=0 or rec["tau"]==0:continue
        model=calc.p.Model(**rec["contract"]["model"]);modes=calc.ab.load_modes(source,rec)
        for root,mode in zip(rec["roots"],modes):
            a,b=mode["arms"];T=np.array([-1.,-1.,1.,1.,1.,-1.])
            left=a["ends"][1];right=T*b["ends"][1]
            dl=calc.p.base.state_matrix(model.arms[0],root["omega"])@left
            dr=-T*(calc.p.base.state_matrix(model.arms[1],root["omega"])@b["ends"][1])
            interfaces.append(dict(mu=rec["mu"],tau=rec["tau"],rank=root["sorted_position"],
                M_left=float(left[5]),M_right=float(right[5]),Q_left=float(left[4]),Q_right=float(right[4]),
                psi_prime_left=float(dl[2]),psi_prime_right=float(dr[2]),
                EI_left=model.arms[0].section.EI,EI_right=model.arms[1].section.EI))
    return summary,interfaces


def figures(folder,rows):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for k,ax in enumerate(axes.flat,1):
        for mu,color,marker in [(0.,"C0","o"),(.01,"C1","s"),(.1,"C2","^")]:
            r=[x for x in rows if x["mu"]==mu and x["rank"]==k]
            ax.scatter([x["beta_deg"] for x in r],[x["F_omega"] for x in r],label=f"mu={mu:g}",color=color,marker=marker)
        ax.set(title=f"Sorted position {k}",xlabel="beta (degrees)",ylabel="omega_k(beta) / omega_k(0)")
        ax.grid(alpha=.2)
    axes[0,0].legend()
    fig.suptitle("Five saved angles, three cuts of the same straight rod; no branch tracking")
    fig.savefig(folder/"low_frequency_ratios.png",dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(12,6.5),constrained_layout=True)
    for col,mu in enumerate([0.,.01,.1]):
        for k,color,marker in [(1,"C0","o"),(2,"C1","s")]:
            r=[x for x in rows if x["mu"]==mu and x["rank"]==k]
            axes[0,col].scatter([x["beta_deg"] for x in r],[100*x["K_axial"] for x in r],label=f"rank {k}",color=color,marker=marker)
            axes[1,col].scatter([x["beta_deg"] for x in r],[100*x["M_arm1"] for x in r],color=color,marker=marker)
        axes[0,col].set(title=f"mu={mu:g}",ylabel="Axial stiffness share (%)",ylim=(-2,70))
        axes[1,col].set(xlabel="beta (degrees)",ylabel="Mass form in arm 1 (%)",ylim=(-2,102))
        axes[1,col].axhline(50,color="k",lw=.7,ls=":")
        for ax in axes[:,col]:ax.grid(alpha=.2)
    axes[0,0].legend()
    fig.suptitle("First pair: energy redistribution and arm mass are separate diagnostics")
    fig.savefig(folder/"first_pair_energy_mass.png",dpi=170);plt.close(fig)


def run(source=SOURCE):
    manifest=calc.verify_bundle(source)
    if manifest["status"]!="PASS":raise ValueError("Completed C0/C1 PASS required")
    import matplotlib
    identity=dict(version="saved-c1-mu-audit-v1",source=source.relative_to(ROOT).as_posix(),
        source_manifest_sha256=calc.sha(source/"manifest.json"),source_artifacts=manifest["artifacts"],
        code_sha256=calc.sha(Path(__file__)),node_helper_sha256=calc.sha(Path(saved.__file__)),
        versions=dict(numpy=np.__version__,matplotlib=matplotlib.__version__),
        method="stored quadrature; global mass MAC within each mu, no cross-mu field identity",
        amplitude_floor=saved.AMPLITUDE_FLOOR,energy_reporting_floor=saved.FLOOR,
        frequency_policy=dict(manifest["identity"]["config"],calculation_mode="plot_only"))
    key=calc.e.digest(identity)[:16];folder=OUTPUT/key
    if (folder/"manifest.json").exists():
        old=json.loads((folder/"manifest.json").read_text())
        assert all(calc.sha(folder/n)==h for n,h in old["artifacts"].items())
        print(json.dumps(dict(folder=str(folder),status="PASS",mode="verified_cache")));return folder
    folder.mkdir(parents=True,exist_ok=False)
    with saved.no_solves(),ExitStack() as stack:
        def forbidden(*a,**kw):raise AssertionError("No new spectra in saved C1 audit")
        for name in ("ritz","count_at","mode","boundary"):
            stack.enter_context(patch.object(calc.p,name,forbidden))
        rows,nodes,macs,matches,packets,gaps,invariant,symmetry=analyze(manifest)
        summary,interfaces=technical(source,manifest)
        calc.write_new(folder/"identity.json",identity)
        for name,data in [("modal_components.csv",rows),("global_mass_MAC.csv",macs),
                          ("discrete_matches.csv",matches),("spectral_gaps.csv",gaps),
                          ("straight_interface_components.csv",interfaces)]:
            csv_new(folder/name,data)
        for name,data in [("baseline_nodes.json",nodes),("pair_subspaces.json",packets),
                          ("straight_mu_invariance.json",invariant),("technical_summary.json",summary),
                          ("symmetric_reflection_parity.json",symmetry)]:
            calc.write_new(folder/name,data)
        events=[dict(kind="LOW_OR_AMBIGUOUS_MAC",**x) for x in matches if x["status"].startswith("UNRESOLVED")]
        # A 1% gap is an explicitly diagnostic alert, not a root-merge tolerance.
        events += [dict(kind="SAMPLED_GAP_BELOW_ONE_PERCENT",**x) for x in gaps if x["relative_gap"]<.01]
        calc.write_new(folder/"events.json",dict(events=events,action="record_only; no new angles",gap_alert_fraction=.01,
            MAC_min=manifest["identity"]["config"]["policy"]["mac_min"],assignment_margin_min=.1))
        figures(folder,rows)
    calc.verify_bundle(source)
    calc.write_new(folder/"manifest.json",dict(status="PASS",identity=identity,spectral_solves=0,basis_evaluations=0,
        artifacts={x.name:calc.sha(x) for x in sorted(folder.iterdir()) if x.is_file()}))
    print(json.dumps(dict(folder=str(folder),status="PASS",primary_modes=len(rows),events=len(events),spectral_solves=0)))
    return folder

