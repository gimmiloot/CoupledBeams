"""Fixed NSL-3A: five saved seeds, two protected substitutions, three roots.

Script Proliferation Control: new contract = validated NSL bundles -> complex
NSL boundary provider -> bounded local KV corrections. This is not a parameter
preset of the uncoupled RLB provider. Reuse NSL materials/geometry, NSL-2R
column scaling and the existing augmented Newton; no prefix/Ritz solver.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scipy.linalg import block_diag, expm
from scripts.analysis.nonsymmetric_laminated_beams import diagnose_mode_recovery as prior
from scripts.lib import nonsymmetric_laminated_kelvin_voigt as model
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as elastic
from scripts.lib import inplane_kelvin_voigt as kv
from scripts.analysis.check_contraction_joint_sensitivity import plain
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv

CONFIG = ROOT/"data/input/nonsymmetric_laminated_beam_kelvin_voigt_pilot.json"
OUTPUT = ROOT/"results/nonsymmetric_laminated_beams/nsl3a_kelvin_voigt_pilot"
read = prior.read


def context(config):
    old, data, properties, _ = prior.source(read(prior.CONFIG))
    fp, ident = prior.identity(read(prior.CONFIG))
    assert fp == "de35ae8dcf5ac831"
    discovered = list((ROOT/"results/nonsymmetric_laminated_beams").rglob("manifest.json"))
    manifests = {}
    for name in config["source_bundles"]:
        paths = [p for p in discovered if p.parent.name == name]
        if len(paths) != 1:
            raise ValueError(f"Expected exactly one source manifest for {name}")
        path = paths[0]; saved = read(path)
        assert all(sha(path.parent/n) == h for n, h in saved["artifact_hashes"].items())
        if name == fp:
            assert saved["identity"] == ident
        manifests[name] = dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path))
    assert old["main_beta_deg"] == config["beta_deg"] == 30.
    assert old["main_kappa"] == config["kappa_theta"] == 1.
    assert old["variant"] == config["variant"] == "PROJECT_REDDY_COUPLED"
    baseline = read(OUTPUT/"baseline.json")
    assert all(sha(ROOT/n) == h for n, h in baseline["sha256"].items())
    r = old["reference"]
    time = r["length"]**2*np.sqrt(r["m"]/r["D"])
    reference = dict(**r, t_ref=float(time), k_theta=config["kappa_theta"]*r["D"]/r["length"],
                     c_theta=config["d_theta"]*r["D"]*time/r["length"])
    return old, data, properties, reference, manifests


def provider(config, old, properties, reference, state, d=None):
    case = old["cases"][state["case"]]
    return model.Provider([properties[n] for n in case], old["lengths"],
        np.deg2rad(config["beta_deg"]), reference["k_theta"],
        reference["c_theta"] if d is None else d*reference["D"]*reference["t_ref"]/reference["length"],
        reference["t_ref"], state["omega"])


def seeds(config, data):
    result = []
    for specification in config["states"]:
        s = data["spectral"][specification["case"]+"_SPRING"]
        mode = s["modes"][specification["position"]-1]
        result.append(dict(**specification, omega=mode["omega"], saved=mode,
            count=s["count"]["roots"][specification["position"]-1],
            gap=min(abs(m["omega"]-mode["omega"]) for m in s["modes"] if m is not mode)))
    return result


def structural(config, old, properties, ref, states):
    policy = config["policy"]; rows = []
    reference_errors = dict(D=abs(ref["D"]/kv.REFERENCE.D-1),
        m=abs(ref["m"]/kv.REFERENCE.m-1), length=abs(ref["length"]/kv.REFERENCE.L-1),
        time=abs(ref["t_ref"]/kv.T_REF-1))
    for state in states:
        pr = provider(config, old, properties, ref, state)
        w = state["omega"]; c0 = np.array(state["saved"]["coefficients"], complex)
        er = elastic.CoupledBeams(pr.properties, pr.lengths, pr.beta,
                                 elastic.Joint("SPRING", pr.k))
        b0 = er.boundary(w, normalized=False)
        raw0, _ = pr.raw(1j*w, c_override=0.)
        alpha0 = pr.c*state["saved"]["delta_psi"]**2/(2*state["saved"]["mass"])
        p = -alpha0+1j*w
        b, bp = pr.raw(p, True)
        z = p*ref["t_ref"]; h = policy["derivative_step_z"]
        bhat, bz = pr.matrices(z, derivative=True)
        fd = (pr.matrices(z+h)[0]-pr.matrices(z-h)[0])/(2*h)
        # An independent imaginary direction also checks holomorphic scaling.
        fdi = (pr.matrices(z+1j*h)[0]-pr.matrices(z-1j*h)[0])/(2j*h)
        offsets = [a.B/a.A for a in pr.properties]
        shifted, shifted_p = pr.raw(p, True, offsets=offsets)
        scale = block_diag(*(np.diag(a.scale) for a in pr.arms))
        sc = np.linalg.solve(scale, elastic.REFLECTION@scale)
        wr = model.reflection_rows(pr.beta)
        hd = block_diag(*(single.state_matrix(p, a) for a in pr.properties))
        eq = pr.evaluate(1j*w, c0, policy["quadrature"])
        basepr = provider(config, old, properties, ref, state, d=0.)
        eq0 = basepr.evaluate(1j*w, c0, policy["quadrature"])
        sigma = np.linalg.svd(basepr.matrices(1j*w*ref["t_ref"])[0], compute_uv=False)
        dd = kv.joint_matrix(p, pr.beta, pr.k, pr.c)-kv.joint_matrix(p, pr.beta, pr.k, 0.)
        only_moment = np.zeros_like(dd); only_moment[2] = dd[2]
        minus = np.vstack([np.eye(6), -np.linalg.solve(np.diag(pr.arms[1].scale),
                              single.REVERSE@np.diag(pr.arms[0].scale))])/np.sqrt(2)
        undamped, _ = pr.raw(p, c_override=0.)
        row = dict(id=state["id"], elastic_boundary=model.relative(raw0, b0),
            derivative_real=model.relative(bz, fd), derivative_imag=model.relative(bz, fdi),
            conjugate=model.relative(pr.raw(p.conjugate())[0], b.conj()),
            axis_boundary=model.relative(b, shifted), axis_derivative=model.relative(bp, shifted_p),
            reflection_boundary=model.relative(b@sc, wr@b),
            reflection_operator=model.relative(hd@elastic.REFLECTION, elastic.REFLECTION@hd),
            reflection_square=float(np.max(abs(sc@sc-np.eye(12)))),
            protected_block_d_independence=model.relative(b@minus, undamped@minus),
            viscous_rows_error=float(np.max(abs(dd-only_moment))),
            joint_rank=int(np.linalg.matrix_rank(kv.joint_matrix(p,pr.beta,pr.k,pr.c))),
            elastic_singular_ratio=float(sigma[-1]/sigma[0]),
            second_singular_ratio=float(sigma[-2]/sigma[0]),
            nearest_gap_relative=state["gap"]/w,
            seed_gates=model.physical_gates(eq0, policy),
            positive_d_saved_form_gates=model.physical_gates(eq, policy),
            seed_mass=eq0["mass"], seed_delta=eq0["delta_psi"],
            seed_mass_difference=abs(eq0["mass"]-state["saved"]["mass"]),
            seed_delta_difference=abs(eq0["delta_psi"]-state["saved"]["delta_psi"]))
        checks = dict(matrix=max(row[n] for n in ("elastic_boundary","conjugate",
            "axis_boundary","axis_derivative","reflection_square")) <= policy["matrix_tol"],
            derivative=max(row["derivative_real"],row["derivative_imag"]) <= policy["derivative_rtol"],
            viscous_rows=row["viscous_rows_error"] == 0., rank=row["joint_rank"] == 6,
            elastic_seed=all(row["seed_gates"].values()),
            elastic_singularity=row["elastic_singular_ratio"] <= policy["sigma_tol"],
            simple=row["second_singular_ratio"] >= policy["simple_sigma_min"] and
                   row["nearest_gap_relative"] > 1e-5 and state["count"]["counts"] ==
                   [state["position"]-1,state["position"]])
        if state["case"] == "MIRROR":
            checks["reflection"] = max(row["reflection_boundary"],row["reflection_operator"],
                row["protected_block_d_independence"]) <= policy["matrix_tol"]
        else:
            checks["not_mirror"] = row["reflection_operator"] > 1e-5
        if state["role"] == "PROTECTED":
            checks["protected_saved_form"] = all(row["positive_d_saved_form_gates"].values())
        row["gates"] = checks; rows.append(row)
    # Same-coefficient symmetric limit only, not an EB/RLB scientific comparison.
    a = properties["BASE"]; arm = kv.Arm("RLB",a.A,a.D,a.m,1.,1/a.S,a.j2)
    p = -.0001+.7j
    hlim = model.relative(single.state_matrix(p,a), kv.state_matrix(p,arm))
    dlim = model.relative(model.state_derivative(p,a),kv.state_matrix(p,arm,derivative=True))
    scale = np.diag(single.SingleBeam(a,1.).scale)
    transfer = model.relative(expm(np.linalg.solve(scale,single.state_matrix(p,a)@scale)*.5),
        expm(np.linalg.solve(scale,kv.state_matrix(p,arm)@scale)*.5))
    # Time-domain virtual power, with arbitrary compatible translations/rotations.
    maps, jump = elastic.trace_maps(np.pi/6,elastic.Joint("SPRING",ref["k_theta"]))
    q = np.array([.2,-.3,.4,-.1]); vel = np.array([-.1,.5,.7,-.2])
    delta, deltav = jump@q, jump@vel
    moment = ref["k_theta"]*delta+ref["c_theta"]*deltav
    force = np.array([.3,-.7])
    forces = [np.r_[e[:2,:2]@force*sign, -moment*sign]
              for e,sign in zip(maps,(1,-1))]
    power = sum(f@(e@vel) for f,e in zip(forces,maps))+ref["k_theta"]*delta*deltav
    balance_error = abs(power+ref["c_theta"]*deltav**2)
    final_gates = dict(seeds=all(all(r["gates"].values()) for r in rows),
        reference=max(reference_errors.values()) < 1e-14,
        rlb_limit=max(hlim,dlim,transfer) <= policy["matrix_tol"],
        passive=ref["c_theta"] >= 0 and power <= 0 and balance_error < 1e-14)
    return dict(status="PASS" if all(final_gates.values()) else "FAIL",
        gates=final_gates, states=rows, reference_errors=reference_errors,
        symmetric_limit=dict(H=hlim,H_p=dlim,transfer=transfer),
        virtual_power=dict(value=power,expected=-ref["c_theta"]*deltav**2,error=balance_error))


def compact(data):
    return {k:v for k,v in data.items() if k != "arms"}


def evaluate_target(config, old, properties, ref, state, out, ledger):
    policy = config["policy"]
    pr = provider(config,old,properties,ref,state)
    w = state["omega"]; elastic_p = 1j*w
    seed = pr.evaluate(elastic_p,np.array(state["saved"]["coefficients"],complex),policy["quadrature"])
    prediction = pr.c*state["saved"]["delta_psi"]**2/(2*state["saved"]["mass"])
    if state["role"] == "PROTECTED":
        p = elastic_p
        coef = np.array(state["saved"]["coefficients"],complex)
        fields = seed
        fine = pr.evaluate(p,coef,policy["quadrature_check"])
        fields["quadrature_relative"] = max(abs(fine[n]/fields[n]-1) for n in ("mass","stiffness"))
        correction = dict(status="EXACT_INACTIVE_BY_SYMMETRY",steps=0,history=[],z=p*pr.t_ref)
    else:
        z0 = (-prediction+1j*w)*pr.t_ref
        _, v = pr.null(z0)
        ledger["primary_complex_calls"] += 1
        # Persist before entering the corrector: an interruption must not hide a solve.
        write_json(out/"ledger.json",ledger)
        criteria = {**kv.CRITERIA,**{n:policy[n] for n in ("newton_residual","newton_step","max_steps")}}
        correction = kv.correct(pr.matrices,z0,v,criteria)
        ledger["newton_updates"] += correction["steps"]
        write_json(out/(state["id"]+"_newton.json"),plain(correction))
        write_json(out/"ledger.json",ledger)
        if correction["status"] != "CONVERGED":
            raise ArithmeticError(f'{state["id"]}: Newton did not converge; retain diagnostic, no automatic retry')
        p = correction["z"]/pr.t_ref
        coef, fields = pr.recover(correction["z"],policy["quadrature"],policy["quadrature_check"])
    b,_ = pr.matrices(p*pr.t_ref)
    vector = coef*pr.columns
    sigma = np.linalg.svd(b,compute_uv=False)
    null = float(np.linalg.norm(b@vector)/(np.linalg.norm(b)*np.linalg.norm(vector)))
    partner = pr.evaluate(p.conjugate(),coef.conjugate(),policy["quadrature"])
    overlap = model.mac(fields,seed,pr.properties)
    gates = model.physical_gates(fields,policy)
    gates.update(null=null <= policy["null_tol"],
        sigma=sigma[-1]/sigma[0] <= policy["sigma_tol"],MAC=overlap >= policy["mac_min"],
        conjugate_matrix=model.relative(pr.matrices(p.conjugate()*pr.t_ref)[0],b.conj()) <= policy["matrix_tol"],
        conjugate_form=all(model.physical_gates(partner,policy).values()),frequency=p.imag > 0)
    actual_alpha = -p.real
    alpha_difference = abs(actual_alpha-fields["alpha_energy"])
    if state["role"] == "PROTECTED":
        gates["protected_delta"] = abs(fields["delta_psi"]) <= policy["inactive_delta_atol"]
        alpha,pred_report,zeta,pred_error = 0.,0.,0.,None
        gates["structural_zero"] = actual_alpha == 0.
        # Numerical residuals are preserved separately, never overwritten with zero.
    else:
        alpha = actual_alpha; pred_report = prediction
        zeta = alpha/abs(p); pred_error = (alpha/prediction-1)
        gates["damping"] = alpha > 0
        gates["exact_alpha"] = alpha_difference <= policy["alpha_identity_atol"]+policy["alpha_identity_rtol"]*abs(alpha)
    row = dict(id=state["id"],case=state["case"],elastic_position=state["position"],
        role=state["role"],status="PASS" if all(gates.values()) else "NUMERICAL_UNRESOLVED",
        omega0=w,abs_delta_psi0=abs(state["saved"]["delta_psi"]),mass0=state["saved"]["mass"],
        d_theta=config["d_theta"],alpha_pred=pred_report,elastic_predictor_raw=prediction,
        p=p,z=p*pr.t_ref,alpha=alpha,zeta=zeta,omega_d=p.imag,frequency_shift=p.imag-w,
        frequency_relative_shift=p.imag/w-1, predictor_relative_difference=pred_error,
        MAC=overlap,alpha_energy_absolute_error=alpha_difference,
        alpha_energy_relative_error=alpha_difference/alpha if alpha else None,
        null_residual=null,singular_ratio=float(sigma[-1]/sigma[0]),
        singular_values=sigma,metrics=compact(fields),gates=gates,
        coefficients=coef,correction=correction,calls=pr.calls,
        fixed_rows=pr.rows,fixed_columns=pr.columns,seed_provenance=state)
    profile = []
    for i,(arm,c) in enumerate(zip(pr.arms,coef.reshape(2,6)),1):
        x=np.linspace(0,arm.length,policy["profile_nodes"])
        y=pr.basis(p,arm,x)[0]@c
        for xx,yy in zip(x,y):
            record=dict(id=state["id"],arm=i,x=xx)
            for name,value in zip(("u","w","psi","N","Q","M"),yy):
                record[name+"_real"]=value.real;record[name+"_imag"]=value.imag
            profile.append(record)
    write_json(out/(state["id"]+".json"),plain(row))
    write_csv(out/(state["id"]+"_profiles.csv"),profile)
    write_json(out/"ledger.json",ledger)
    return row


def identity(config, manifests):
    baseline=read(OUTPUT/"baseline.json")
    assert all(sha(ROOT/n)==h for n,h in baseline["sha256"].items())
    files=[CONFIG,Path(__file__),Path(model.__file__)]
    ident=dict(config=config,baseline_sha256=sha(OUTPUT/"baseline.json"),
        source_manifests=manifests,files={p.relative_to(ROOT).as_posix():sha(p) for p in files},
        versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16],ident


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--structural-only",action="store_true")
    args=parser.parse_args()
    config=read(CONFIG);old,data,properties,ref,manifests=context(config)
    fingerprint,ident=identity(config,manifests);out=OUTPUT/fingerprint
    if (out/"manifest.json").exists():
        saved=read(out/"manifest.json")
        assert saved["identity"]==ident
        assert all(sha(out/n)==h for n,h in saved["artifact_hashes"].items())
        result=read(out/"result.json")
        print(json.dumps(dict(status=result["status"],cache="VALIDATED",fingerprint=fingerprint,
                              new_complex_solves=0)));return
    if (out/"ledger.json").exists() and read(out/"ledger.json")["primary_complex_calls"]:
        raise RuntimeError("Unfinished solve ledger exists; no automatic replay of roots")
    out.mkdir(parents=True,exist_ok=True)
    states=seeds(config,data)
    checks=structural(config,old,properties,ref,states)
    write_json(out/"structural.json",plain(checks))
    write_json(out/"preflight_contract.json",dict(config=config,reference=ref,identity=ident))
    if not args.compute or checks["status"]!="PASS":
        print(json.dumps(dict(status=checks["status"],stage="STRUCTURAL",fingerprint=fingerprint,new_complex_solves=0)))
        if checks["status"]!="PASS":raise SystemExit(1)
        return
    ledger=dict(primary_complex_calls=0,auxiliary_complex_calls=0,newton_updates=0,
                protected_newton_calls=0,elastic_root_solves=0,prefix_searches=0)
    result=dict(status="IN_PROGRESS",reference=ref,structural=checks,states=[],solve_ledger=ledger)
    try:
        for state in states:
            row=evaluate_target(config,old,properties,ref,state,out,ledger)
            result["states"].append(row)
            write_json(out/"result.json",plain(result))
            print(state["id"],row["status"],row["p"],"MAC",row["MAC"],flush=True)
        result["status"]="PASS" if all(s["status"]=="PASS" for s in result["states"]) else "PARTIAL_PASS"
    except Exception as error:
        import traceback
        result.update(status="PARTIAL_PASS",stop_reason=f"{type(error).__name__}: {error}",
                      traceback=traceback.format_exc())
    result["protected_sources_unchanged"]=all(sha(ROOT/n)==h for n,h in read(OUTPUT/"baseline.json")["sha256"].items())
    write_json(out/"result.json",plain(result))
    comparison=[]
    for row in result["states"]:
        names=("id","case","elastic_position","omega0","abs_delta_psi0","d_theta",
               "alpha_pred","alpha","omega_d","zeta","predictor_relative_difference",
               "frequency_shift","frequency_relative_shift","MAC","status")
        comparison.append({**{n:row[n] for n in names},"p_real":row["p"].real,"p_imag":row["p"].imag})
    if comparison:write_csv(out/"comparison.csv",comparison)
    write_json(out/"manifest.json",dict(identity=ident,
        git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        command=subprocess.list2cmdline([sys.executable,*sys.argv]),
        artifact_hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}))
    write_json(OUTPUT/"current.json",dict(fingerprint=fingerprint,status=result["status"]))
    print(json.dumps(dict(status=result["status"],fingerprint=fingerprint,ledger=ledger,
                          stop_reason=result.get("stop_reason"))))


if __name__=="__main__":
    main()
