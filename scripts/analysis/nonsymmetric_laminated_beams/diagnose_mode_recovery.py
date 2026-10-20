"""NSL-2R: saved-data recovery diagnosis and TWO fixed local derivative roots.

Script Proliferation Control: this input is a frozen prior bundle, not a new
spectral case. The workflow never invokes prefix/Ritz searches. Keeping this
diagnostic entry separate preserves the historical NSL-1/2 source hashes.
It reuses their H, boundary assembly, physical residuals and mass integration.
"""
from __future__ import annotations
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT))
import numpy as np
import scipy
from scipy.linalg import block_diag, eig, expm
from scipy.optimize import brentq
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as joint
from scripts.analysis.nonsymmetric_laminated_beams import verify_coupled_joint as previous
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.analysis.check_contraction_joint_sensitivity import plain

CONFIG=ROOT/"data/input/nsl2r_numerical_diagnostics.json"
OUTPUT=ROOT/"results/nonsymmetric_laminated_beams/nsl2r_numerical_diagnostics"
read=previous.read


class BoundedBeam(single.SingleBeam):
    """Local basis change for the six existing spatial exponential solutions.

    Growing real columns start at x=L, decaying columns at x=0; conjugate
    columns are represented by real/imaginary parts about L/2. No entries of H
    are changed and no conditions are imposed on the columns themselves.
    """
    @lru_cache(maxsize=64)
    def decomposition(self,omega):
        h=self.balanced(omega)
        values,vectors=eig(h)
        columns=[]
        for i,value in enumerate(values):
            v=vectors[:,i]
            if abs(value.imag)<1e-10*max(1,abs(value)):
                if abs(value.imag)>1e-13:
                    raise ArithmeticError("Unexpected nonreal spatial exponent")
                v=v.real
                if v[np.argmax(abs(v))]<0:v=-v
                columns.append((float(value.real),v,self.length if value.real>0 else 0.,"REAL"))
            elif value.imag>0:
                # Constant phase choice in each local interval. Retain the
                # numerical real part of lambda; no altered dispersion law.
                v=v*np.exp(-1j*np.angle(v[np.argmax(abs(v))]))
                columns.append((complex(value),v,self.length/2,"PAIR"))
        columns.sort(key=lambda z:(z[3],z[0].imag if z[3]=="PAIR" else z[0]))
        if sum(2 if z[3]=="PAIR" else 1 for z in columns)!=6:
            raise ArithmeticError("Spatial eigensystem is not a complete real basis")
        residual=np.linalg.norm(h@vectors-vectors*values)/max(np.linalg.norm(h)*np.linalg.norm(vectors),1e-300)
        if residual>1e-12 or np.linalg.cond(vectors)>1e10:
            raise ArithmeticError("Unsupported ill-conditioned spatial diagonalization")
        return columns,{"spatial_eigenvalues":values,"spatial_vectors_condition":float(np.linalg.cond(vectors)),
            "spatial_eigen_residual":float(residual)}

    def basis(self,omega,x,derivative=False):
        columns,_=self.decomposition(float(omega));points=np.atleast_1d(x);result=[]
        for value,vector,origin,kind in columns:
            solution=np.exp(value*(points-origin))[:,None]*vector[None,:]
            if derivative:solution=solution*value
            if kind=="PAIR":
                result.extend([solution.real,solution.imag])
            else:
                result.append(solution.real)
        return self.scale[None,:,None]*np.stack(result,axis=2)


def bounded_model(model):
    if isinstance(model,joint.CoupledBeams):
        result=joint.CoupledBeams(model.properties,model.lengths,model.beta,model.joint)
        result.arms=[BoundedBeam(a.a,a.length) for a in model.arms]
        return result
    return BoundedBeam(model.a,model.length)


def raw_boundary(model,omega,origin="CENTER"):
    arms=model.arms if isinstance(model,joint.CoupledBeams) else [model]
    matrices=[]
    for arm in arms:
        if origin=="CENTER":
            b=arm.basis(omega,[0,arm.length])
        else:
            h=arm.balanced(omega)
            b=np.array([arm.scale[:,None]*expm(h*x) for x in (0,arm.length)])
        matrices.append(b)
    if len(arms)==1:return matrices[0][:,:3].reshape(6,6)
    clamps=block_diag(*(b[0,:3] for b in matrices))
    end=block_diag(*(b[1] for b in matrices))
    return np.vstack([clamps,joint.joint_matrix(model.beta,model.joint)@end])


def null_diagnostics(matrix):
    rows=np.linalg.norm(matrix,axis=1)
    b=matrix/rows[:,None]
    _,s,vh=np.linalg.svd(b)
    columns=np.linalg.norm(b,axis=0)
    _,sc,vc=np.linalg.svd(b/columns)
    return {"singular_values":s,"nonzero_condition":float(s[0]/s[-2]),
        "full_condition":float(s[0]/s[-1]),"row_norm_ratio":float(max(rows)/min(rows)),
        "column_norm_ratio_after_rows":float(max(columns)/min(columns)),
        "column_balanced_singular_values":sc,"column_balanced_nonzero_condition":float(sc[0]/sc[-2])},vh[-1],vc[-1]/columns


def metrics(d):
    names=("mass","stiffness","boundary_scaled","joint_scaled","joint_physical_residual",
           "ode_scaled","energy_relative","quadrature_relative","singular_ratio","nonzero_condition","delta_psi")
    answer={n:d[n] for n in names if n in d}
    arms=d["arms"] if "arms" in d else [d]
    answer["clamp_physical_values"]=[v["end_states"][[0] if "arms" in d else [0,1],:3] for v in arms]
    answer["component_peaks"]=[np.max(abs(v["states"][:,:3]),axis=0) for v in arms]
    return answer


def recover_mode(model,omega,policy):
    """Minimal accepted change: column equilibration before the existing SVD.

    Bc=B diag(1/norm(B[:,j])); c=diag(1/norm(B[:,j])) vc.
    Physical evaluation, row scaling, mass normalization and tolerances remain
    the historical routines. The bounded basis above is diagnostic only.
    """
    b=model.boundary(omega)
    scales=np.linalg.norm(b,axis=0)
    if np.any(scales<=0) or not np.isfinite(scales).all():
        raise ArithmeticError("Zero/nonfinite fundamental column")
    _,sigma,vh=np.linalg.svd(b/scales)
    coefficient=vh[-1]/scales
    before=model.evaluate(omega,coefficient,policy["mode_quadrature"])
    coefficient=coefficient/np.sqrt(before["mass"])
    d=model.evaluate(omega,coefficient,policy["mode_quadrature"])
    fine=model.evaluate(omega,coefficient,policy["mode_quadrature_check"])
    d.update(quadrature_relative=max(abs(fine[n]/d[n]-1) for n in ("mass","stiffness")),
             singular_ratio=float(sigma[-1]/sigma[0]),
             nonzero_condition=float(sigma[0]/sigma[-2]),
             unnormalized_mass=before["mass"],column_scales=scales)
    return coefficient,d


def physical_gates(d,policy):
    return {
        "clamps":d["boundary_scaled"]<=policy["boundary_tol"],
        "joint":d.get("joint_scaled",0.)<=policy["boundary_tol"],
        "ODE":d["ode_scaled"]<=policy["ode_tol"],
        "energy":d["energy_relative"]<=policy["energy_tol"],
        "mass":abs(d["mass"]-1)<=policy["mass_gram_tol"],
        "quadrature":d["quadrature_relative"]<=policy["energy_tol"]}


def cancellation(model,omega,coefficient):
    arms=model.arms if isinstance(model,joint.CoupledBeams) else [model]
    result=[]
    for arm,c in zip(arms,np.asarray(coefficient).reshape(-1,6)):
        b=arm.basis(omega,[0,arm.length]);y=b@c
        # Sum of individual column contributions relative to an interior field
        # amplitude, so exact endpoint zeros do not cause divisions by zero.
        grid=arm.basis(omega,np.linspace(0,arm.length,17))@c
        scale=np.maximum(np.max(abs(grid),axis=0),1e-30)
        result.append({"end_column_cancellation":np.einsum("xij,j->xi",abs(b),abs(c))/scale,
            "balanced_fundamental_condition_at_end":float(np.linalg.cond(b[1]/arm.scale[:,None])),
            "balanced_fundamental_max":float(np.max(abs(b/arm.scale[None,:,None]))),
            "endpoint_values":y})
    return result


def high_precision_diagnostic(model,omega,coefficient,policy,dps):
    """Same saved omega/H, 50-digit expm/SVD; no root solve.

    Boundary values are evaluated at high precision. The unchanged Gauss grid
    supplies component peaks in double precision; this is a diagnostic only,
    not a replacement acceptance metric for the bounded physical forms.
    """
    sys.path.insert(0,str(ROOT/"tmp/nsl_tools"))
    import mpmath as mp
    arms=model.arms if isinstance(model,joint.CoupledBeams) else [model]
    with mp.workdps(dps):
        def matrix(a):return mp.matrix(np.asarray(a).tolist())
        ends=[];errors=[]
        for arm in arms:
            h=matrix(arm.balanced(omega));scale=mp.diag(arm.scale.tolist())
            bs=[scale*mp.expm(h*mp.mpf(float(x-arm.length/2))) for x in (0.,arm.length)]
            ends.append(bs)
            old=arm.basis(omega,[0,arm.length])
            ref=np.array([np.array(b.tolist(),float) for b in bs])
            errors.append(previous.relative(old,ref))
        dim=6*len(arms);raw=mp.zeros(dim)
        if len(arms)==1:
            for i in range(6):
                for j in range(6):raw[i,j]=ends[0][i//3][i%3,j]
        else:
            endpoint=mp.zeros(12)
            for arm in range(2):
                for i in range(6):
                    for j in range(6):
                        endpoint[arm*6+i,arm*6+j]=ends[arm][1][i,j]
                        if i<3:raw[3*arm+i,6*arm+j]=ends[arm][0][i,j]
            rows=matrix(joint.joint_matrix(model.beta,model.joint))*endpoint
            for i in range(6):
                for j in range(12):raw[6+i,j]=rows[i,j]
        balanced=raw.copy()
        for i in range(dim):
            norm=mp.sqrt(sum(raw[i,j]**2 for j in range(dim)))
            for j in range(dim):balanced[i,j]/=norm
        _,sigma,v=mp.svd_r(balanced)
        hp=mp.matrix([v[dim-1,j] for j in range(dim)])
        old=mp.matrix(np.asarray(coefficient).tolist())
        hp_float=np.array(list(hp),float)
        peaks=model.evaluate(omega,hp_float,policy["mode_quadrature"])
        legacy=model.evaluate(omega,np.asarray(coefficient),policy["mode_quadrature"])
        def physical_clamp(c,field):
            ds=field["arms"] if "arms" in field else [field]
            maxima=[]
            for i,(arm,d) in enumerate(zip(arms,ds)):
                cc=c[6*i:6*i+6,0]
                for end in ([0] if len(arms)==2 else [0,1]):
                    yy=ends[i][end]*cc
                    scale=np.maximum(np.max(abs(d["states"][:,:3]),axis=0),1e-30)
                    maxima.extend(float(abs(yy[j]))/scale[j] for j in range(3))
            return max(maxima)
        return {"dps":dps,"expm_relative_errors":errors,
            "hp_null_clamp_same_peak_definition":physical_clamp(hp,peaks),
            "legacy_null_hp_propagation_clamp":physical_clamp(old,legacy),
            "hp_singular_values":[str(s) for s in sigma],
            "hp_nonzero_condition":str(sigma[0]/sigma[dim-2]),
            "frequency_unchanged":repr(omega),"new_root_solves":0}


def source(config):
    oldconfig=read(ROOT/"data/input/nonsymmetric_laminated_beam_coupled_joint.json")
    fp,ident=previous.identity(oldconfig)
    directory=ROOT/config["source_bundle"];manifest=read(directory/"manifest.json")
    assert directory==previous.OUTPUT/fp and manifest["identity"]==ident
    assert all(sha(directory/n)==h for n,h in manifest["artifact_hashes"].items())
    baseline=read(OUTPUT/"baseline.json")
    assert all(sha(ROOT/n)==h for n,h in baseline["sha256"].items())
    data=read(directory/"result.json")
    props,_=previous.materials(oldconfig)
    single_model=single.SingleBeam(props["HHLL"],2.)
    coupled=joint.CoupledBeams([props["HHLL"],props["LLHH"]],[1.,1.],0.,joint.Joint("RIGID"))
    return oldconfig,data,props,{"COUPLED":coupled,"SINGLE":single_model}


def diagnose(config,oldconfig,data,models):
    policy=oldconfig["policy"];result=[]
    for kind,position in config["problem_cases"]:
        model=models[kind]
        row=data["spectral"]["STRAIGHT_RIGID"]["modes"][position-1] if kind=="COUPLED" else data["straight_single"]["count"]["roots"][position-1]
        omega=row["omega"]
        null,coefficient,cc=null_diagnostics(raw_boundary(model,omega))
        d0=model.evaluate(omega,coefficient,policy["mode_quadrature"])
        c,d=model.mode(omega,policy)
        balanced=model.evaluate(omega,cc,policy["mode_quadrature"])
        bounded=bounded_model(model);bc,bd=bounded.mode(omega,policy)
        bnull,_,_=null_diagnostics(raw_boundary(bounded,omega))
        upstream,_,_=null_diagnostics(raw_boundary(model,omega,"LEFT"))
        record={"model":kind,"position":position,"omega":omega,
            "original":metrics(d),"before_mass_normalization":metrics(d0),
            "column_normalization_only":metrics(balanced),
            "original_boundary":null,"left_origin_boundary":upstream,
            "bounded":metrics(bd),"bounded_boundary":bnull,
            "original_cancellation":cancellation(model,omega,c),
            "bounded_cancellation":cancellation(bounded,omega,bc),
            "high_precision":high_precision_diagnostic(model,omega,c,policy,config["high_precision_dps"])}
        result.append(record)
        print(kind,position,"BC",d["boundary_scaled"],"bounded",bd["boundary_scaled"],flush=True)
    return result


def direct_control(oldconfig,data,models,out):
    """Seven SAVED roots, original physical bases, two representations."""
    policy=oldconfig["policy"];rows=[];profiles=[]
    runtime={};result={"status":"IN_PROGRESS","new_root_solves":0,"modes":[]}
    for kind,model in models.items():
        saved=data["spectral"]["STRAIGHT_RIGID"]["modes"] if kind=="COUPLED" else data["straight_single"]["count"]["roots"]
        modes=[]
        for k,row in enumerate(saved):
            w=row["omega"];coefficient,d=recover_mode(model,w,policy)
            record={"model":kind,"position":k+1,"omega_old":w,"omega_used":w,
                    "frequency_change":0.,"coefficients":coefficient,
                    "metrics":metrics(d),"gates":physical_gates(d,policy)}
            result["modes"].append(record);modes.append((coefficient,d))
            rows.append({"model":kind,"position":k+1,"omega":w,
                **{n:d.get(n) for n in ("boundary_scaled","joint_scaled","ode_scaled","mass","energy_relative","quadrature_relative")}})
            arms=model.arms if kind=="COUPLED" else [model]
            for armidx,(arm,c) in enumerate(zip(arms,coefficient.reshape(-1,6))):
                x=np.r_[0.,np.linspace(0,arm.length,101)[1:-1],arm.length]
                y=arm.basis(w,x)@c
                for xx,yy in zip(x,y):
                    profiles.append({"model":kind,"position":k+1,"arm":armidx+1,"x":float(xx),
                        **dict(zip(("u","w","psi","N","Q","M"),map(float,yy)))})
        runtime[kind]=modes
    comparisons=[]
    coupled=models["COUPLED"];single_model=models["SINGLE"]
    for k in range(7):
        cc,cm=runtime["COUPLED"][k];sc,sm=runtime["SINGLE"][k]
        wc=data["spectral"]["STRAIGHT_RIGID"]["modes"][k]["omega"]
        ws=data["straight_single"]["count"]["roots"][k]["omega"]
        mapped={"arms":[]};field_errors=[]
        for i,arm in enumerate(cm["arms"]):
            x=arm["x"] if i==0 else single_model.length-arm["x"]
            y=single_model.basis(ws,x)@sc
            if i==1:y=y@joint.CONTINUOUS
            mapped["arms"].append({"states":y,"weights":arm["weights"]})
        overlap=joint.inner(cm,mapped,coupled.properties)
        mapped_mass=joint.inner(mapped,mapped,coupled.properties)
        sign=1 if overlap>=0 else -1
        error=0.
        for i,(arm,reference) in enumerate(zip(cm["arms"],mapped["arms"])):
            diff=arm["states"]-sign*reference["states"]
            peaks=np.maximum(np.max(abs(arm["states"]),axis=0),1e-30)
            field_errors.append(np.max(abs(diff),axis=0)/peaks)
            error+=np.sum(arm["weights"]*np.einsum("ni,ij,nj->n",diff[:,:3],coupled.properties[i].mass_matrix(),diff[:,:3]))
        center=single_model.basis(ws,[1.])[0]@sc
        interface=np.r_[center,joint.CONTINUOUS@center]
        interface_residual=joint.joint_matrix(0.,joint.Joint("RIGID"))@interface
        comp={"position":k+1,"frequency_relative":abs(wc/ws-1),
            "one_minus_mac":max(0.,1-overlap**2/(cm["mass"]*mapped_mass)),
            "signed_shape_mass_error":float(np.sqrt(error)),
            "fieldwise_relative_u_w_psi_N_Q_M":field_errors,
            "single_mapped_interface_residual":interface_residual}
        comp["gates"]={"frequency":comp["frequency_relative"]<=policy["frequency_relative_tol"],
            "MAC":comp["one_minus_mac"]<=policy["one_minus_mac_tol"],
            "physical_shape":comp["signed_shape_mass_error"]<=policy["signed_shape_mass_error_tol"],
            "interface":float(np.max(abs(interface_residual)))<=policy["boundary_tol"]}
        comparisons.append(comp)
    # Same physical mass metric used in the historical checks.
    cg=np.array([[joint.inner(a[1],b[1],coupled.properties) for b in runtime["COUPLED"]] for a in runtime["COUPLED"]])
    fields=np.array([d["states"][:,:3] for c,d in runtime["SINGLE"]])
    weights=runtime["SINGLE"][0][1]["weights"]
    sg=np.einsum("kni,ij,lnj,n->kl",fields,single_model.a.mass_matrix(),fields,weights)
    result.update(comparisons=comparisons,
        mass_gram_errors={"COUPLED":float(np.max(abs(cg-np.eye(7)))),"SINGLE":float(np.max(abs(sg-np.eye(7))))})
    result["gates"]={
        "all_physical_modes":all(all(m["gates"].values()) for m in result["modes"]),
        "representations":all(all(m["gates"].values()) for m in comparisons),
        "mass_orthogonality":max(result["mass_gram_errors"].values())<=policy["mass_gram_tol"]}
    result["status"]="PASS" if all(result["gates"].values()) else "PARTIAL_PASS"
    write_csv(out/"straight_modes.csv",rows);write_csv(out/"straight_profiles.csv",profiles)
    return result


def active_derivative(config,oldconfig,data,props,ledger,out):
    policy=oldconfig["policy"];settings=config["derivative"]
    saved=data["spectral"][settings["case"]];index=settings["position"]-1
    central=saved["modes"][index];w0=central["omega"]
    kscale=oldconfig["reference"]["D"]/oldconfig["reference"]["length"]
    model=joint.CoupledBeams([props["HHLL"],props["HHLL"]],oldconfig["lengths"],
        np.deg2rad(saved["beta_deg"]),joint.Joint("SPRING",kscale))
    baseline=model.evaluate(w0,np.asarray(central["coefficients"]),policy["mode_quadrature"])
    assert abs(baseline["mass"]-central["mass"])<1e-13
    assert abs(baseline["delta_psi"]-central["delta_psi"])<1e-10
    predictor=kscale*central["delta_psi"]**2/central["mass"]
    gap=min(abs(w0-m["omega"]) for i,m in enumerate(saved["modes"]) if i!=index)
    singular=np.linalg.svd(model.boundary(w0),compute_uv=False)
    isolation=saved["count"]["roots"][index]
    if predictor<=0 or gap/w0<1e-5 or isolation["counts"]!=[index,index+1] or singular[-2]/singular[0]<1e-10:
        raise ArithmeticError("Preset active mode is not demonstrably isolated/simple")
    delta=settings["delta_kappa"]
    shift=delta*predictor/(2*w0)
    radius=min(w0-isolation["bracket"][0],isolation["bracket"][1]-w0,gap/4,max(10*abs(shift),1e-8))
    result={"status":"IN_PROGRESS","position":index+1,"omega0":w0,
        "delta_psi0":central["delta_psi"],"mass0":central["mass"],"D_ref_over_l_ref":kscale,
        "energy_prediction_dlambda_dkappa":predictor,"simple_gap_relative":gap/w0,
        "central_singular_values":singular,"saved_isolation":isolation,"local_halfwidth":radius,
        "settings":settings,"local_roots":[]}
    for sign in (-1,1):
        kappa=1+sign*delta
        current=joint.CoupledBeams(model.properties,model.lengths,model.beta,joint.Joint("SPRING",kscale*kappa))
        lo,hi=w0-radius,w0+radius
        function=lambda w:float(np.linalg.det(current.boundary(w)))
        if function(lo)*function(hi)>=0:
            raise ArithmeticError("Preset local interval has no bracket; no enlargement allowed")
        if ledger["active_derivative_roots"]>=2:
            raise ArithmeticError("Two-root diagnostic budget exhausted")
        ledger["active_derivative_roots"]+=1
        w,info=brentq(function,lo,hi,xtol=settings["root_xtol"],rtol=settings["root_rtol"],full_output=True)
        c,d=recover_mode(current,w,policy)
        overlap=joint.inner(d,baseline,current.properties)
        mac=overlap**2/(d["mass"]*baseline["mass"])
        row={"kappa":kappa,"omega":w,"local_bracket":[lo,hi],"iterations":info.iterations,
             "converged":bool(info.converged),"coefficients":c,"metrics":metrics(d),
             "MAC_to_saved_mode":mac,"gates":physical_gates(d,policy)}
        row["gates"]["MAC"]=1-mac<=settings["one_minus_mac_tol"]
        result["local_roots"].append(row)
        write_json(out/"active_derivative_checkpoint.json",plain(result))
    minus,plus=[r["omega"] for r in result["local_roots"]]
    derivative=(plus**2-minus**2)/(2*delta)
    # Conservative solver-tolerance bound, distinct from rounding-only error.
    e1=settings["root_xtol"]+settings["root_rtol"]*abs(minus)
    e2=settings["root_xtol"]+settings["root_rtol"]*abs(plus)
    uncertainty=(2*abs(minus)*e1+e1**2+2*abs(plus)*e2+e2**2)/(2*delta)
    rounding=[]
    for digits in (15,12,9):
        m,p=round(minus,digits),round(plus,digits)
        value=(p*p-m*m)/(2*delta)
        rounding.append({"decimal_places":digits,"derivative":value,
            "relative_to_full_precision":abs(value/derivative-1),
            "relative_to_predictor":abs(value/predictor-1)})
    result.update(central_difference=derivative,relative_error=abs(derivative/predictor-1),
        solver_uncertainty_bound=uncertainty,solver_uncertainty_relative=uncertainty/predictor,
        decimal_rounding=rounding,
        subtractive_vs_factored_difference=abs(derivative-(plus-minus)*(plus+minus)/(2*delta)))
    result["gates"]={"positive":derivative>0,
        "derivative":result["relative_error"]<=settings["relative_tol"],
        "root_error_budget":result["solver_uncertainty_relative"]<=settings["relative_tol"]/10,
        "physical_forms":all(all(r["gates"].values()) for r in result["local_roots"])}
    result["status"]="PASS" if all(result["gates"].values()) else "PARTIAL_PASS"
    write_csv(out/"active_derivative_roots.csv",[{"kappa":r["kappa"],"omega":r["omega"],
        "MAC":r["MAC_to_saved_mode"],**{n:r["metrics"][n] for n in ("mass","delta_psi","boundary_scaled","joint_scaled","ode_scaled","energy_relative")}} for r in result["local_roots"]])
    return result


def identity(config):
    baseline=read(OUTPUT/"baseline.json")
    assert all(sha(ROOT/n)==h for n,h in baseline["sha256"].items())
    paths=[CONFIG,Path(__file__)]
    ident={"config":config,"baseline_sha256":sha(OUTPUT/"baseline.json"),
        "source_hashes":baseline["sha256"],"files":{p.relative_to(ROOT).as_posix():sha(p) for p in paths},
        "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}}
    return hashlib.sha256(json.dumps(ident,sort_keys=True).encode()).hexdigest()[:16],ident


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diagnose-only",action="store_true")
    parser.add_argument("--compute",action="store_true")
    args=parser.parse_args()
    config=read(CONFIG);oldconfig,data,props,models=source(config)
    fingerprint,ident=identity(config);out=OUTPUT/fingerprint
    if (out/"manifest.json").exists():
        manifest=read(out/"manifest.json")
        assert manifest["identity"]==ident
        assert all(sha(out/n)==h for n,h in manifest["artifact_hashes"].items())
        saved=read(out/"result.json")
        if saved["status"]!="IN_PROGRESS" or not args.compute:
            print(json.dumps({"status":saved["status"],"cache":"VALIDATED","fingerprint":fingerprint,"new_eigenvalue_solves":0}));return
    if not(args.compute or args.diagnose_only):
        print(json.dumps({"source_identity":"PASS","fingerprint":fingerprint,"new_eigenvalue_solves":0}));return
    out.mkdir(parents=True,exist_ok=True)
    result={"status":"IN_PROGRESS","config":config,"solve_ledger":{"straight_root_refinements":0,"active_derivative_roots":0,"prefix_searches":0,"ritz_eigenproblems":0,"complex_damped_roots":0}}
    try:
        result["diagnosis"]=diagnose(config,oldconfig,data,models)
        write_json(out/"result.json",plain(result))
        if args.compute:
            result["direct_control"]=direct_control(oldconfig,data,models,out)
            # No straight refinement is needed if unchanged-root recovery passes.
            # An unresolved result is retained, not automatically retried.
            result["active_derivative"]=active_derivative(config,oldconfig,data,props,result["solve_ledger"],out)
            result["status"]="PASS" if result["direct_control"]["status"]==result["active_derivative"]["status"]=="PASS" else "PARTIAL_PASS"
    except Exception as error:
        result["status"]="PARTIAL_PASS" if "diagnosis" in result else "FAIL"
        result["stop_reason"]=f"{type(error).__name__}: {error}"
        import traceback
        result["traceback"]=traceback.format_exc()
    result["source_preservation"]="PASS" if all(sha(ROOT/n)==h for n,h in ident["source_hashes"].items()) else "FAIL"
    write_json(out/"result.json",plain(result))
    write_json(out/"manifest.json",{"identity":ident,
        "command":subprocess.list2cmdline([sys.executable,*sys.argv]),
        "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "artifact_hashes":{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}})
    write_json(OUTPUT/"current.json",{"fingerprint":fingerprint,"status":result["status"]})
    print(json.dumps({"status":result["status"],"fingerprint":fingerprint,"solve_ledger":result["solve_ledger"],"stop_reason":result.get("stop_reason")}))
    if "stop_reason" in result:raise SystemExit(1)


if __name__=="__main__":
    main()
