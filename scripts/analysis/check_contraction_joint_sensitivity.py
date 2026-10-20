"""Fixed G20 c-closure study; elastic hard gates precede local KV checks.

New diagnostic contract (not a rigid-model preset): admissible trace spaces,
variational index, nested-domain ordering, conditional complex confirmation.
Reuses unchanged arm helpers, source/reference validation and atomic writers.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy

from scripts.analysis import verify_mindlin_herrmann_timoshenko_single_rod as single
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import contraction_joint_sensitivity as study
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as rigid

CONFIG = ROOT/"data/input/contraction_joint_sensitivity.json"
OUTPUT = ROOT/"results/contraction_joint_sensitivity"


def plain(value):
    if isinstance(value, np.ndarray):
        return plain(value.tolist())
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    if isinstance(value, (complex, np.complexfloating)):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, np.generic):
        return value.item()
    return value


def references(config):
    base, model, length, sources = single.check_inputs()
    refs, manifests = {}, {}
    for key in ("direct_bundle", "rigid_bundle", "pole_bundle"):
        directory = ROOT/config[key]
        manifest = json.loads((directory/"manifest.json").read_text(encoding="utf-8"))
        for name, digest in manifest["artifact_hashes"].items():
            if sha(directory/name) != digest:
                raise ValueError(f"Historical artifact hash mismatch: {directory/name}")
        refs[key] = json.loads((directory/"result.json").read_text(encoding="utf-8"))
        if refs[key]["coefficients"] != model.coefficients:
            raise ValueError("Historical reference coefficients differ from the accepted model")
        manifests[config[key]+"/manifest.json"] = sha(directory/"manifest.json")
    s = model.section
    if (length != 1 or config["material"] != dict(E=s.E, rho=s.rho, nu=s.nu, kappa=s.K)
        or config["geometry"] != dict(b=s.width, h=s.thickness, lengths=[.5, .5])
        or config["beta_deg"] != [0., 45.] or config["kappa_theta"] != 1
        or config["l_ref"] != .5 or config["policy"]["guard_position"] != 7):
        raise ValueError("This study is restricted to the requested G20 cases")
    if refs["pole_bundle"]["statuses"]["MHTIM_GENERAL_BETA_JOINT_GATE"] != "PASS":
        raise ValueError("Unaccepted historical pole/reference bundle")
    return model, refs, sources, manifests


def identity(config, sources, manifests):
    paths = [CONFIG, Path(__file__), Path(study.__file__), Path(single.__file__),
        ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py",
        ROOT/"scripts/lib/mindlin_herrmann_timoshenko_joint.py",
        ROOT/"scripts/lib/coupled_longitudinal_comparators.py",
        ROOT/"scripts/lib/bishop_longitudinal.py",
        ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
        ROOT/"scripts/lib/reddy_inplane_geometry.py", ROOT/"scripts/lib/inplane_kelvin_voigt.py",
        single.CONFIG, ROOT/"tests/data/reddy_four_ply_isotropic_limit_cases.json"]
    data = {"config": config, "files": {p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        "sources": sources, "reference_manifests": manifests,
        "versions": {"python":platform.python_version(), "numpy":np.__version__, "scipy":scipy.__version__}}
    return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()[:16], data


def diagnostics_record(diag):
    return {k:plain(v) for k,v in diag.items() if k not in ("values", "gradients", "weights")}


def solve(model, lengths, case, k, catalog, policy):
    frequencies, search = study.prefix_roots(model, lengths, case, k, catalog, policy)
    roots, all_diag, vectors = [], [], []
    for index, omega in enumerate(frequencies, 1):
        coeff, diag = study.elastic_mode(model, lengths, omega, case, k, policy["quadrature_order"])
        study.validate_mode(diag, policy)
        operator = study.state_operator(model, omega, case)
        def evaluator(arm, x):
            return (study.basis(model,lengths[arm],omega,x,case)@coeff[arm],
                study.basis(model,lengths[arm],omega,x,case,1)@coeff[arm], operator)
        refined = study.field_diagnostics(model,lengths,case,k,evaluator,policy["quadrature_check_order"])
        diag["quadrature_relative"] = max(abs(refined[name]/diag[name]-1) for name in ("mass","elastic_energy"))
        if diag["quadrature_relative"] > policy["energy_tol"]:
            raise ArithmeticError("Elastic quadrature gate failed")
        vectors.append(study.mass_vectors(model, case, diag))
        all_diag.append(diag)
        roots.append({"sorted_position":index, "guard":index == 7, "omega":omega,
            "fstar":omega/(2*math.pi), "coefficients":coeff.tolist(),
            "diagnostics":diagnostics_record(diag)})
    gram = np.array(vectors)@np.array(vectors).T
    error = float(np.max(abs(gram-np.eye(7))))
    if error > policy["mass_gram_tol"]:
        raise ArithmeticError(f"Mass orthogonality failed: {case.key}: {error}")
    return {"key":case.key, "case":asdict(case), "roots":roots,
        "search":search, "mass_gram":gram.tolist(), "mass_gram_error":error}, all_diag


def rigid_regression(model, config, result, refs, diags):
    length = sum(config["geometry"]["lengths"])
    inventory = sorted([(r["omega"], block, r) for block in ("mh", "timoshenko")
        for r in refs["direct_bundle"][block]["roots"]])[:7]
    frozen = refs["rigid_bundle"]["splits"]
    half = next(s for s in frozen if s["split"] == .5)
    frozen_freq = sorted(r["omega"] for block in ("mh", "timoshenko") for r in half[block]["roots"])[:7]
    nodes, weights = np.polynomial.legendre.leggauss(config["policy"]["quadrature_order"])
    checks = []
    for root, diag, (old_w, block, saved), frozen_w in zip(result["roots"], diags, inventory, frozen_freq):
        old_mode = mh.finite_mode(model, length, old_w, block, config["policy"]["quadrature_order"])
        states = []
        for arm in range(2):
            x = (nodes+1)*.25
            global_x = x if arm == 0 else 1-x
            local = np.zeros((len(x), 8))
            local[:, rigid.BLOCK_INDICES[block]] = mh.finite_state_basis(model, 1., old_w, global_x, block)@old_mode["coefficients"]
            if arm == 1:
                local *= rigid.REFLECTION
            states.append(local)
        ref_diag = {**diag, "values":np.array(states)}
        case = study.Case(0., rotation="RIGID")
        a, b = study.mass_vectors(model, case, diag), study.mass_vectors(model, case, ref_diag)
        sign = 1 if a@b >= 0 else -1
        mac = float((a@b)**2/((a@a)*(b@b)))
        l2 = float(np.linalg.norm(sign*a-b)/np.linalg.norm(b))
        frequency_error = max(abs(root["omega"]/old_w-1), abs(root["omega"]/frozen_w-1))
        if frequency_error > config["policy"]["rigid_relative_tol"] or l2 > 5e-8:
            raise ArithmeticError("Exact rigid/direct G20 regression failed")
        checks.append({"k":root["sorted_position"], "direct_family":block,
            "relative_frequency_error":frequency_error, "mass_MAC":mac, "mass_L2":l2})
    return {"status":"PASS", "exact_constraint":True, "checks":checks}


def comparisons(model, cases, all_diags, policy):
    ordering, overlaps, metrics = [], [], []
    for beta in (0., 45.):
        group = [c for c in cases if c["case"]["beta"] == beta and c["case"]["rotation"] == "SPRING"]
        mh_cases = {c["case"]["closure"]:c for c in group if c["case"]["variant"] == "mindlin_herrmann"}
        for k in range(7):
            w = {cl:c["roots"][k]["omega"] for cl,c in mh_cases.items()}
            violation = max(0., w["FREE"]-w["CONTINUOUS"], w["CONTINUOUS"]-w["BLOCKED"])/w["FREE"]
            row = {"beta":beta, "k":k+1, "relative_violation":violation,
                "FREE":w["FREE"], "CONTINUOUS":w["CONTINUOUS"], "BLOCKED":w["BLOCKED"],
                "blocked_minus_free_percent":100*(w["BLOCKED"]/w["FREE"]-1)}
            ordering.append(row)
            if violation > policy["ordering_relative_tol"]:
                raise ArithmeticError(f"HARD GATE: min-max ordering failure: {row}")
            others = {c["case"]["variant"]:c["roots"][k]["omega"] for c in group if c["case"]["variant"] != "mindlin_herrmann"}
            metrics.append({**row, **others,
                "continuous_minus_free_percent":100*(w["CONTINUOUS"]/w["FREE"]-1),
                "blocked_minus_continuous_percent":100*(w["BLOCKED"]/w["CONTINUOUS"]-1),
                "elementary_vs_continuous_percent":100*(others["elementary"]/w["CONTINUOUS"]-1),
                "love_vs_continuous_percent":100*(others["rayleigh_love_planar"]/w["CONTINUOUS"]-1),
                "love_vs_elementary_percent":100*(others["rayleigh_love_planar"]/others["elementary"]-1)})
        for first, second in itertools.combinations(group, 2):
            ca, cb = study.Case(**first["case"]), study.Case(**second["case"])
            common = not (ca.is_mh and cb.is_mh)
            va = np.array([study.mass_vectors(model, ca, d, common) for d in all_diags[first["key"]]])
            vb = np.array([study.mass_vectors(model, cb, d, common) for d in all_diags[second["key"]]])
            va /= np.linalg.norm(va, axis=1)[:,None]
            vb /= np.linalg.norm(vb, axis=1)[:,None]
            matrix = np.minimum(1., (va@vb.T)**2)
            rows = []
            for i in range(7):
                row_other = max(np.delete(matrix[i], i))
                col_other = max(np.delete(matrix[:,i], i))
                accepted = (matrix[i,i] >= policy["shape_mac_min"] and
                            matrix[i,i]-max(row_other,col_other) >= policy["shape_margin_min"])
                rows.append({"k":i+1,"diagonal":float(matrix[i,i]),
                    "row_best":int(np.argmax(matrix[i]))+1, "column_best":int(np.argmax(matrix[:,i]))+1,
                    "margin":float(matrix[i,i]-max(row_other,col_other)), "accepted_fixed_case_correspondence":bool(accepted)})
            overlaps.append({"beta":beta,"first":first["key"],"second":second["key"],
                "metric":"common u,w,theta mass-weighted projection" if common else "full MH physical mass MAC",
                "matrix":matrix.tolist(), "rows":rows,
                "qualification":"fixed-case correspondence only; no across-beta tracking, no root reassignment"})
    return {"status":"PASS", "ordering":ordering, "metrics":metrics, "overlaps":overlaps}


def component_comparison(model, result, all_diags):
    rows = []
    for beta in (0.,45.):
        group = [c for c in result["cases"] if c["case"]["beta"]==beta
            and c["case"]["variant"]=="mindlin_herrmann" and c["case"]["rotation"]=="SPRING"]
        for first,second in itertools.combinations(group,2):
            for index,(da,db) in enumerate(zip(all_diags[first["key"]],all_diags[second["key"]]),1):
                case = study.Case(**first["case"])
                va,vb = study.mass_vectors(model,case,da),study.mass_vectors(model,case,db)
                sign = 1 if va@vb>=0 else -1
                a,b,weights = da["values"],sign*db["values"],da["weights"]
                cond = max(da["nonzero_condition"],db["nonzero_condition"])
                qscale = max(np.max(abs(a[:,:,:4])),np.max(abs(b[:,:,:4])))
                pscale = max(np.max(abs(a[:,:,4:])),np.max(abs(b[:,:,4:])))
                record = {"beta":beta,"k":index,"first":first["case"]["closure"],"second":second["case"]["closure"]}
                for j,name in enumerate(rigid.STATE_ORDER):
                    norm_a = float(np.sqrt(np.sum(weights*a[:,:,j]**2)))
                    norm_b = float(np.sqrt(np.sum(weights*b[:,:,j]**2)))
                    diff = float(np.sqrt(np.sum(weights*(a[:,:,j]-b[:,:,j])**2)))
                    floor = 128*np.finfo(float).eps*cond*(qscale if j<4 else pscale)
                    record[name+"_relative_L2"] = diff/max(norm_a,norm_b) if max(norm_a,norm_b)>floor else None
                    record[name+"_norm_first"] = norm_a
                    record[name+"_norm_second"] = norm_b
                    record[name+"_difference_L2"] = diff
                for arm in range(2):
                    record[f"c{arm+1}_first"] = float(da["end_states"][arm,1])
                    record[f"c{arm+1}_second_aligned"] = float(sign*db["end_states"][arm,1])
                    record[f"R{arm+1}_first"] = float(da["end_states"][arm,5])
                    record[f"R{arm+1}_second_aligned"] = float(sign*db["end_states"][arm,5])
                rows.append(record)
    return rows


def complex_compute(model, config, result, all_diags, out):
    if result.get("elastic_status")!="PASS":
        raise ArithmeticError("Complex stage requires all elastic gates")
    policy = config["policy"]
    group = [r for r in result["cases"] if r["case"]["beta"]==45
        and r["case"]["variant"]=="mindlin_herrmann" and r["case"]["rotation"]=="SPRING"]
    candidates = []
    for index in range(6):
        k = index+1
        relevant = [o for o in result["comparison"]["overlaps"] if o["beta"]==45
            and "mindlin_herrmann" in o["first"] and "mindlin_herrmann" in o["second"]]
        correspondence = all(o["rows"][index]["accepted_fixed_case_correspondence"] for o in relevant)
        gaps = []
        for r in group:
            frequencies = np.array([v["omega"] for v in r["roots"]])
            gaps.append(float(min(abs(np.delete(frequencies,index)-frequencies[index]))))
        relative_gap = min(g/r["roots"][index]["omega"] for g,r in zip(gaps,group))
        activity = min(r["roots"][index]["predictor"]["sensitivity_s"] for r in group)
        candidates.append({"k":k,"correspondence":correspondence,"min_relative_gap":relative_gap,
            "min_activity":activity,"eligible":bool(correspondence and relative_gap>policy["simple_relative_gap_min"]
                and activity>policy["active_spring_energy_fraction_min"]),"gaps":gaps})
    selected = sorted([r for r in candidates if r["eligible"]],key=lambda r:r["min_activity"],reverse=True)[:policy["complex_mode_limit"]]
    result["complex_selection"] = {"rule":"largest minimum spring energy fraction among eligible first-six positions",
        "candidates":candidates,"selected_positions":[r["k"] for r in selected]}
    result["complex_modes"] = []
    profiles, table = [], []
    for selected_mode in selected:
        index = selected_mode["k"]-1
        for item,gap in zip(group,selected_mode["gaps"]):
            case = study.Case(**item["case"])
            root = item["roots"][index]
            print("complex",case.key,"k",index+1,flush=True)
            value,provider = study.complex_mode(model,config["geometry"]["lengths"],case,root["omega"],
                root["predictor"]["alpha"],result["k_theta"],result["c_theta"],result["t_ref"],
                all_diags[item["key"]][index],gap,policy)
            for arm,length in enumerate(config["geometry"]["lengths"]):
                x = np.linspace(0.,length,201)
                states = provider.values(value["p"],arm,x)@value["coefficients"][arm]
                for xx,yy in zip(x,states):
                    profiles.append({"case":case.key,"k":index+1,"arm":arm+1,"x":float(xx),
                        **{name+suffix:float(getattr(y,component)) for name,y in zip(rigid.STATE_ORDER,yy)
                           for suffix,component in (("_real","real"),("_imag","imag"))}})
            value["diagnostics"] = diagnostics_record(value["diagnostics"])
            value.update({"case":case.key,"k":index+1,"omega_elastic":root["omega"],
                "basis_reference_rates":provider.rates,"basis_pivots":provider.pivots,"basis_scales":provider.scales})
            result["complex_modes"].append(plain(value))
            table.append({"case":case.key,"k":index+1,"p_real":value["p"].real,"p_imag":value["p"].imag,
                "alpha_predictor":value["alpha_predictor"],"predictor_relative_error":value["predictor_relative_error"],
                "energy_identity_relative":value["energy_identity_relative"],
                "alpha_identity_relative":value["alpha_identity_relative"],"mass_MAC":value["mass_MAC_with_elastic"]})
            write_json(out/"result.json",plain(result))
    if profiles:
        write_csv(out/"complex_profiles.csv",profiles)
        write_csv(out/"complex_modes.csv",table)
    result["complex_status"] = "PASS" if selected else "NOT_RUN_NO_ELIGIBLE_MODE"
    result["status"] = "COMPLETE"


def elastic_compute(model, config, refs, out):
    policy, lengths = config["policy"], config["geometry"]["lengths"]
    k = config["kappa_theta"]*model.coefficients["B"]/config["l_ref"]
    tref = config["l_ref"]**2*math.sqrt(model.coefficients["m"]/model.coefficients["B"])
    ctheta = policy["d_theta"]*model.coefficients["B"]*tref/config["l_ref"]
    result = {"status":"RUNNING", "complex_status":"NOT_RUN", "parameters":config,
        "coefficients":model.coefficients, "k_theta":k,"t_ref":tref,"c_theta":ctheta,
        "structure":study.structural_checks(model, lengths, k), "cases":[]}
    write_json(out/"result.json", plain(result))
    all_diags = {}
    catalog = refs["pole_bundle"]["fixed_arm_poles"]
    configurations = [study.Case(0., rotation="RIGID")]
    configurations += [study.Case(beta, closure) for beta in config["beta_deg"] for closure in config["closures"]]
    configurations += [study.Case(beta, variant=variant) for beta in config["beta_deg"] for variant in config["comparators"]]
    for case in configurations:
        print("elastic", case.key, flush=True)
        record, diags = solve(model, lengths, case, k, catalog, policy)
        result["cases"].append(record)
        all_diags[case.key] = diags
        if case.rotation == "RIGID":
            result["rigid_regression"] = rigid_regression(model, config, record, refs, diags)
        for root in record["roots"]:
            d = root["diagnostics"]
            jump = d["delta_theta"] if case.rotation=="SPRING" else 0.
            root["predictor"] = {"alpha":ctheta*jump**2/2,
                "alpha_t_ref":ctheta*jump**2*tref/2,
                "zeta":ctheta*jump**2/(2*root["omega"]),
                "sensitivity_s":k*jump**2/root["omega"]**2,
                "rotation":"constraint, no finite damper" if case.rotation == "RIGID" else "SPRING"}
        write_json(out/"result.json", plain(result))
    result["comparison"] = comparisons(model, result["cases"], all_diags, policy)
    result["component_comparison"] = component_comparison(model,result,all_diags)
    result["elastic_status"] = "PASS"
    result["status"] = "ELASTIC_PASS"
    write_json(out/"result.json", plain(result))
    return result, all_diags


def save_tables(model, config, result, out):
    rows, profiles = [], []
    for record in result["cases"]:
        case = study.Case(**record["case"])
        for root in record["roots"]:
            d = root["diagnostics"]
            rows.append({"case":case.key,"beta":case.beta,"variant":case.variant,
                "closure":case.closure if case.is_mh else "NA", "rotation":case.rotation,
                "k":root["sorted_position"],"guard":root["guard"],"omega":root["omega"],
                "fstar":root["fstar"],"mass":d["mass"],"delta_theta":d["delta_theta"],
                "c_or_gradient_mass_fraction":d["contraction_or_gradient_mass"]/d["mass"],
                "axial_energy_fraction":d["axial_energy"]/d["elastic_energy"],
                "spring_energy_fraction":d["spring_energy"]/d["elastic_energy"],
                "joint_scaled_max":max(d["joint_scaled"]),"energy_relative":d["energy_relative"],
                "alpha_predictor":root["predictor"]["alpha"],"sensitivity_s":root["predictor"]["sensitivity_s"]})
            for arm, length in enumerate(config["geometry"]["lengths"]):
                x = np.linspace(0, length, 201)
                y = study.basis(model, length, root["omega"], x, case)@np.array(root["coefficients"])[arm]
                names = rigid.STATE_ORDER if case.is_mh else ("u","w","theta","N","Q","M")
                for xx, yy in zip(x,y):
                    state = dict(zip(names, map(float, yy)))
                    profiles.append({"case":case.key,"k":root["sorted_position"],"arm":arm+1,
                        "x":float(xx), **{name:state.get(name, "") for name in rigid.STATE_ORDER}})
    write_csv(out/"frequencies.csv", rows)
    write_csv(out/"profiles.csv", profiles)
    write_csv(out/"comparison.csv", result["comparison"]["metrics"])
    write_csv(out/"component_comparison.csv", result["component_comparison"])


def plot_saved(out):
    """One local c/R view from saved mass-normalized fields; zero eigensolves."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = list(csv.DictReader((out/"profiles.csv").open(encoding="utf-8")))
    result = json.loads((out/"result.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(2,2,figsize=(10,6),sharex=True,layout="constrained")
    colors = {"FREE":"#176b9d","CONTINUOUS":"#dc8b20","BLOCKED":"#9e315b"}
    styles = {"FREE":"-","CONTINUOUS":"--","BLOCKED":"-."}
    for row,k in enumerate((1,3)):
        for closure in ("FREE","CONTINUOUS","BLOCKED"):
            key = study.Case(45.,closure).key
            values = [v for v in data if v["case"]==key and int(v["k"])==k and int(v["arm"])==1]
            x = np.array([float(v["x"]) for v in values])
            sign = 1
            if closure!="FREE":
                # BLOCKED end c is zero: align by a saved full-mass dot product.
                reference = [v for v in data if v["case"]==study.Case(45.,"FREE").key and int(v["k"])==k and int(v["arm"])==1]
                dot = sum(float(a["u"])*float(b["u"])+float(a["w"])*float(b["w"])
                    +result["coefficients"]["r"]/result["coefficients"]["m"]*float(a["theta"])*float(b["theta"])
                    +result["coefficients"]["j"]/result["coefficients"]["m"]*float(a["c"])*float(b["c"])
                    for a,b in zip(reference,values))
                sign = 1 if dot>=0 else -1
            for col,field in enumerate(("c","R")):
                scale = 1 if field=="c" else 1e6
                visible = x>=.35
                axes[row,col].plot(x[visible],sign*scale*np.array([float(v[field]) for v in values])[visible],
                    styles[closure],color=colors[closure],label=closure,lw=1.8)
                axes[row,col].set_xlim(.35,.5)
                axes[row,col].grid(alpha=.25)
                axes[row,col].set_ylabel(f"k = {k}: "+("c" if field=="c" else r"$10^6 R$"))
    axes[0,0].legend(frameon=False,fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel(r"$x_1$ (clamp $0$, joint $0.5$)")
    fig.suptitle(r"G20, $\beta=45^\circ$, $\kappa_\theta=1$: contraction near the joint"+"\n"
        +"Sorted positions with verified fixed-case correspondence; physical modal mass = 1",fontsize=12)
    fig.savefig(out/"contraction_near_joint.png",dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--elastic-only", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args(argv)
    if args.plot_only:
        pointer = json.loads((OUTPUT/"current.json").read_text(encoding="utf-8"))
        out = OUTPUT/pointer["fingerprint"]
        manifest = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
        if any(sha(out/p)!=h for p,h in manifest["artifact_hashes"].items()):
            raise ValueError("Saved-data artifact hash mismatch")
        plot_saved(out)
        manifest["artifact_hashes"]["contraction_near_joint.png"] = sha(out/"contraction_near_joint.png")
        write_json(out/"manifest.json",manifest)
        print("Plot from saved profiles, zero root evaluations:",out)
        return
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    model, refs, sources, manifests = references(config)
    fingerprint, ident = identity(config, sources, manifests)
    if not args.compute:
        print(json.dumps({"inputs":"PASS", "fingerprint":fingerprint, "mode":"read-only check"}))
        return
    out = OUTPUT/fingerprint
    out.mkdir(parents=True, exist_ok=True)
    manifest = {"identity":ident, "head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "command":subprocess.list2cmdline([sys.executable,*sys.argv]), "executable":sys.executable}
    if (out/"manifest.json").exists():
        old = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
        if old["identity"] == ident and all(sha(out/p)==h for p,h in old.get("artifact_hashes",{}).items()):
            saved = json.loads((out/"result.json").read_text(encoding="utf-8"))
            if saved.get("status") == "COMPLETE" or (args.elastic_only and saved.get("elastic_status") == "PASS"):
                print("Validated cache, zero root evaluations:",out)
                return
    try:
        result, diags = elastic_compute(model, config, refs, out)
        save_tables(model, config, result, out)
        if not args.elastic_only:
            complex_compute(model,config,result,diags,out)
        write_json(out/"result.json", plain(result))
        plot_saved(out)
    except Exception as error:
        checkpoint = json.loads((out/"result.json").read_text(encoding="utf-8")) if (out/"result.json").exists() else {}
        write_json(out/"failure.json", {"status":"HARD_GATE_FAILED", "error":str(error),
            "complex_status":"STOPPED" if checkpoint.get("complex_modes") else "NOT_RUN",
            "completed_complex_modes":len(checkpoint.get("complex_modes",[]))})
        raise
    manifest["artifact_hashes"] = {p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"}
    write_json(out/"manifest.json", manifest)
    write_json(OUTPUT/"current.json", {"fingerprint":fingerprint,"directory":str(out),"status":result["status"]})
    print(result["status"], out)


if __name__ == "__main__":
    main()
