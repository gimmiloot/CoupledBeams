"""NSL-1 only: two reductions, source audit, one CC spectrum by propagation/Ritz.

Script Proliferation Control: nonzero B/j1 cannot be a symmetric-RLB preset.
Reuses layer integration and the independent RLB single-beam essential basis.
No angular assembly, joint parameters, damped roots or parameter sweep.
"""
from __future__ import annotations
from dataclasses import asdict, replace
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scripts.lib import nonsymmetric_laminated_beam as nsl
from scripts.lib import reddy_symmetric_laminated_beam as rlb
from scripts.analysis.reproduce_bishop_literature import sha, write_csv, write_json
from scripts.analysis.check_contraction_joint_sensitivity import plain

CONFIG = ROOT/"data/input/nonsymmetric_laminated_beam_single_rod.json"
OUTPUT = ROOT/"results/nonsymmetric_laminated_beams/single_beam_verification"


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def symbolic_audit():
    if importlib.util.find_spec("sympy") is None:
        sys.path.insert(0, str(ROOT/"tmp/nsl_tools"))
    import sympy as s
    A, B, D, S, m, j, j1, k, w, p = s.symbols("A B D S m j j1 k omega p")
    # Independently substitute partial_x=-i*k, partial_t^2=-omega^2
    # into energy balances, order axial/transverse/moment.
    own = s.Matrix([[m*w*w-A*k*k, 0, j1*w*w-B*k*k],
                    [0, m*w*w-S*k*k, -s.I*k*S],
                    [j1*w*w-B*k*k, s.I*k*S, j*w*w-D*k*k-S]])
    source = s.Matrix([[-k*k*B, s.I*k*S, j*w*w-k*k*D-S],
                       [0, m*w*w-k*k*S, -s.I*k*S],
                       [m*w*w-k*k*A, 0, -k*k*B]])
    assert own.subs(j1, 0)[[2, 1, 0], :] == source
    polynomial = s.Poly(s.expand(source.det()/S), k)
    expected = [
        A*D-B*B,
        w*w*(m*B*B-m*D*S-A*(m*D+j*S))/S,
        m*w*w*(A*(j*w*w-S)+w*w*(m*D+j*S))/S,
        m*m*w**4*(1-j*w*w/S)]
    coefficients = [polynomial.coeff_monomial(k**degree) for degree in (6, 4, 2, 0)]
    assert all(s.simplify(a-b) == 0 for a, b in zip(coefficients, expected))
    g = s.Matrix([[A, 0, B], [0, S, 0], [B, 0, D]])
    c = s.zeros(3); c[1, 2] = S
    v = s.diag(0, 0, S)
    mass = s.Matrix([[m, 0, j1], [0, m, 0], [j1, 0, j]])
    gi, gc = g.inv(), g.inv()*c  # symbolic, not numerical full-ABD inversion
    h = (-gc).row_join(gi).col_join((v+p*p*mass-c.T*gc).row_join(c.T*gi))
    delta = A*D-B*B
    displayed = s.Matrix([[0,0,0,D/delta,0,-B/delta],[0,0,-1,0,1/S,0],
        [0,0,0,-B/delta,0,A/delta],[m*p*p,0,j1*p*p,0,0,0],
        [0,m*p*p,0,0,0,0],[j1*p*p,0,j*p*p,0,1,0]])
    assert s.simplify(h-displayed) == s.zeros(6)
    t = s.diag(1,-1,-1,1,-1,-1)
    assert s.simplify(h.subs({B:-B,j1:-j1})-t*h*t) == s.zeros(6)
    canonical = s.zeros(6); canonical[:3, 3:] = s.eye(3); canonical[3:, :3] = -s.eye(3)
    assert s.simplify(h.T*canonical+canonical*h) == s.zeros(6)
    rho, area, inertia = s.symbols("rho area I")
    printed_d = rho*area**2*w**4*(1-rho*inertia*w*w/S)
    derived_d = expected[-1].subs({m:rho*area,j:rho*inertia})
    difference = s.factor(derived_d-printed_d)
    assert difference != 0
    return {"status": "PASS_LOCAL_EQUATIONS_WITH_SOURCE_PRINT_DISCREPANCIES",
        "sympy_version": s.__version__, "Lean": "not available; direct variation plus exact SymPy verification",
        "energy_state_equals_control": True, "work_dual_Hamilton_structure": True,
        "stack_reversal_similarity": True, "source_eq9_row_permutation": [2,1,0],
        "det_source_divisor": "S", "coefficients_k6_k4_k2_k0": [str(s.factor(a)) for a in coefficients],
        "source_eq11_abc": "EXACT_MATCH", "source_eq11d": "PRINTED_IDENTITY_FAIL",
        "eq11d_printed": str(printed_d), "eq11d_derived": str(derived_d),
        "eq11d_difference": str(difference),
        "eq4b_warning": "Printed nu12*E1 contradicts reciprocity/4e; existing Reddy nu12*E2 retained. Q12 drops out of Qbar11 at 0/90.",
        "source_no_silent_repair": True}


def stack(material, angles, thickness):
    return rlb.integrate_laminate([rlb.Ply(material, angle, thickness/len(angles)) for angle in angles])


def hhll(config, layout=None):
    c = config["spectral_case"]
    m0 = rlb.OrthotropicLamina(**c["base_material"])
    materials = {}
    for label, factor in (("H", 1+c["chi"]), ("L", 1-c["chi"])):
        materials[label] = replace(m0, **{n:factor*getattr(m0,n) for n in ("E1","E2","G12","G13","G23")}, name=label)
    return rlb.integrate_laminate([rlb.Ply(materials[label], c["angle_deg"], c["thickness"]/4)
                                  for label in (layout or c["stack"])]), m0


def relative(a, b):
    a,b=np.asarray(a),np.asarray(b)
    return float(np.max(abs(a-b))/max(np.max(abs(a)), np.max(abs(b)), 1e-300))


def coefficient_audit(config):
    c = config["bachoo"]
    mat = rlb.OrthotropicLamina(**c["material"])
    b, h = c["width"], c["thickness"]
    den = mat.E1-mat.E2*mat.nu12**2
    source_rows = []
    for n in c["repetitions"]:
        section = stack(mat, c["project_stack_bottom_to_top"]*n, h)
        prop, audit = nsl.source_bachoo_properties(section, width=b, K=c["K"], provenance=mat.name)
        formula = [b*mat.E1*(mat.E1+mat.E2)*h/(2*den),
                   b*mat.E1*(mat.E1-mat.E2)*h*h/(8*n*den),
                   b*mat.E1*(mat.E1+mat.E2)*h**3/(24*den),
                   c["K"]*b*h*(mat.G13+mat.G23)/2]
        obtained = [prop.A,prop.B,prop.D,prop.S]
        errors = [abs(a/z-1) for a,z in zip(obtained,formula)]
        reverse = rlb.integrate_laminate(tuple(reversed(section.plies)))
        rev, _ = nsl.source_bachoo_properties(reverse,width=b,K=c["K"],provenance="reversed "+mat.name)
        source_rows.append({"n":n,"formula_6":formula,"computed":obtained,"relative_errors":errors,
            "B_times_n":prop.B*n,"j1_scaled":abs(prop.j1)/(prop.m*h),
            "reverse_A":relative(section.A,reverse.A),"reverse_B":relative(section.B,-reverse.B),
            "reverse_D":relative(section.D,reverse.D),"reverse_S":abs(rev.S/prop.S-1),
            "reverse_mass":relative([prop.m,prop.j1,prop.j2],[rev.m,-rev.j1,rev.j2]),"reduction":audit})
    section, m0 = hhll(config)
    cc = config["spectral_case"]
    h, b, chi = cc["thickness"], cc["width"], cc["chi"]
    q0 = rlb.transformed_reduced_stiffness(m0, 0.)
    matrix_errors = [relative(section.A,h*q0),relative(section.B,-chi*h*h*q0/4),
                     relative(section.D,h**3*q0/12)]
    props, reductions, checks = {}, {}, []
    for variant in (nsl.SOURCE_BACHOO,nsl.PROJECT_REDDY_COUPLED):
        prop,audit = nsl.reduce_laminate(section,width=b,K=cc["K"],variant=variant,provenance="RLB-2E M0 H/H/L/L chi=.4")
        props[variant],reductions[variant] = prop,audit
        checks.append({"variant":variant,"B2_over_AD":prop.B**2/(prop.A*prop.D),
                       "Deff_over_D":(prop.D-prop.B**2/prop.A)/prop.D,
                       "j1_scaled":abs(prop.j1)/(prop.m*h)})
    # Previously verified symmetric cross-ply, isotropic and H/L control sections.
    historical = rlb.OrthotropicLamina(25.,1.,.25,.5,.5,.2,1.,"Reddy verified ratios")
    iso = rlb.OrthotropicLamina(1.,1.,.3,1/2.6,1/2.6,1/2.6,1.,"RLB-1C-ISO")
    sym_sections = [(name,stack(material,angles,.05)) for name,material,angles in (
        ("Reddy0",historical,[0]),("Reddy90",historical,[90]),
        ("Reddy0090s",historical,[0,90,90,0]),("isotropic4",iso,[0,90,90,0]))]
    sym_sections += [(name,hhll(config,list(name))[0]) for name in ("HLLH","LHHL")]
    symmetric = []
    for name, sec in sym_sections:
        old = rlb.reduce_to_beam_properties(sec,width=.2,K=cc["K"])
        new, aud = nsl.project_reddy_coupled_properties(sec,width=.2,K=cc["K"],provenance=name)
        errors = [abs(a/b-1) for a,b in zip([new.A,new.D,new.S,new.m,new.j2],[old.A,old.D,old.S,old.m,old.J])]
        op = nsl.state_matrix(1.7j,new).real
        ref = rlb.combined_state_matrix(1.7,old)
        scale = np.array([1,1,1,new.A,new.D,new.D])
        symmetric.append({"case":name,"property_errors":errors,
            "B_relative":abs(new.B)/np.sqrt(new.A*new.D),"j1_relative":abs(new.j1)/np.sqrt(new.m*new.j2),
            "state_relative":relative(op*scale[None,:]/scale[:,None],ref*scale[None,:]/scale[:,None]),
            "reduction_error":aud["reduction_relative_error"]})
    metrics = matrix_errors+[abs(x["B2_over_AD"]-.12) for x in checks]+[abs(x["Deff_over_D"]-.88) for x in checks]
    for x in source_rows:
        metrics += x["relative_errors"]+[x[k] for k in ("j1_scaled","reverse_A","reverse_B","reverse_D","reverse_S","reverse_mass")]
    for x in symmetric:
        metrics += x["property_errors"]+[x["B_relative"],x["j1_relative"],x["state_relative"],x["reduction_error"]]
    if max(metrics) > config["policy"]["coefficient_tol"]:
        raise ArithmeticError("Constitutive/symmetric/source-6 gate failed")
    return props[nsl.PROJECT_REDDY_COUPLED], section, {"status":"PASS","source_eq6":source_rows,
        "HHLL_matrix_relative":matrix_errors,"HHLL_reductions":reductions,"HHLL_checks":checks,
        "symmetric_limits":symmetric,"maximum_error":max(metrics)}


def spectral(config, prop, section, result, out):
    policy,c = config["policy"],config["spectral_case"]
    beam=nsl.SingleBeam(prop,c["length"])
    spectral=result["spectral"]={"status":"IN_PROGRESS","ritz":[],"search":{},"modes":[]}
    spectral["segment_green_relative"]=beam.green_identity(.7,policy["count_segments"])
    if spectral["segment_green_relative"]>policy["energy_tol"]:
        raise ArithmeticError("Exact-segment Green identity failed")
    # Freeze independent Ritz results before propagator roots are computed.
    spectra=[]
    for order in policy["ritz_orders"]:
        w,v,d=nsl.ritz_spectrum(prop,c["length"],order,policy["ritz_quadrature"])
        spectra.append(w)
        spectral["ritz"].append({"order":order,"omega":w,"diagnostics":d})
        result["solve_ledger"]["ritz_eigenproblems"]+=1
    convergence=np.abs(spectra[-1]/spectra[-2]-1)
    monotone=float(np.max(np.diff(np.array(spectra),axis=0)/np.array(spectra[:-1])))
    spectral.update(ritz_convergence=convergence,ritz_monotonicity_max=monotone)
    if max(convergence)>policy["ritz_convergence_tol"] or monotone>policy["ritz_monotonicity_tol"]:
        raise ArithmeticError("Fixed Ritz sequence did not converge")
    result["solve_ledger"]["propagation_prefix_searches"]+=1
    roots=beam.prefix(policy,spectral["search"])
    result["solve_ledger"]["propagation_roots"]=len(roots)
    rev_section=rlb.integrate_laminate(tuple(reversed(section.plies)))
    rev,_=nsl.project_reddy_coupled_properties(rev_section,width=c["width"],K=c["K"],provenance="physical stack reversal")
    reversed_beam=nsl.SingleBeam(rev,c["length"])
    mass_matrix=prop.mass_matrix()
    modes,rows,profiles,allmetrics=[],[],[],[]
    for i,root in enumerate(roots):
        omega=root["omega"]
        coeff,d=beam.mode(omega,policy)
        x,weights=d["x"],d["weights"]
        q=d["states"][:,:3]
        rq,rdq=nsl.ritz_fields(v[:,i],c["length"],x)
        inner=float(np.sum(weights*np.einsum("ni,ij,nj->n",q,mass_matrix,rq)))
        sign=1 if inner>=0 else -1
        rq*=sign;rdq*=sign
        ritzmass=float(np.sum(weights*np.einsum("ni,ij,nj->n",rq,mass_matrix,rq)))
        mac=inner**2/(d["mass"]*ritzmass)
        shape_error=float(np.sqrt(np.sum(weights*np.einsum("ni,ij,nj->n",q-rq,mass_matrix,q-rq))))
        # Reflection exact similarity, transformed actual shape and same omega.
        h=nsl.state_matrix(1j*omega,prop)
        hr=nsl.state_matrix(1j*omega,rev)
        reflection_operator=relative(hr,nsl.REVERSE@h@nsl.REVERSE)
        reflected_coeff=nsl.REVERSE@coeff
        reflected=reversed_beam.evaluate(omega,reflected_coeff,policy["mode_quadrature"])
        reflection_fields=relative(reflected["states"],d["states"]@nsl.REVERSE)
        reflection_boundary=float(np.linalg.norm(reversed_beam.boundary(omega)@reflected_coeff)/np.linalg.norm(reflected_coeff))
        diag={k:z for k,z in d.items() if k not in ("x","weights","states","derivatives")}
        diag.update(frequency_relative=abs(omega/w[i]-1),one_minus_mac=max(0.,1-mac),
            signed_shape_mass_error=shape_error,ritz_mass=ritzmass,ritz_sign=sign,
            reversal_operator=reflection_operator,reversal_fields=reflection_fields,
            reversal_boundary=reflection_boundary,reversal_mass=abs(reflected["mass"]/d["mass"]-1),
            reversal_energy=reflected["energy_relative"])
        gates={"frequency":diag["frequency_relative"]<=policy["frequency_relative_tol"],
            "boundary":diag["boundary_scaled"]<=policy["boundary_tol"],"ode":diag["ode_scaled"]<=policy["ode_tol"],
            "energy":diag["energy_relative"]<=policy["energy_tol"],
            "quadrature":diag["quadrature_relative"]<=policy["energy_tol"],
            "shape":diag["one_minus_mac"]<=policy["one_minus_mac_tol"] and shape_error<=policy["signed_shape_mass_error_tol"],
            "reversal":max(reflection_operator,reflection_fields,reflection_boundary,diag["reversal_mass"])<=policy["reverse_tol"]}
        modes.append(q);allmetrics.append(gates)
        spectral["modes"].append({"position":i+1,"guard":i==6,"omega":omega,"ritz_omega":w[i],
                                  "coefficients":coeff,"ritz_coefficients":v[:,i]*sign,"diagnostics":diag,"gates":gates})
        rows.append({"position":i+1,"guard":i==6,"omega_propagation":omega,"omega_ritz":w[i],
                     **{k:diag[k] for k in ("frequency_relative","mass","energy_relative","boundary_scaled","ode_scaled","one_minus_mac","signed_shape_mass_error")}})
        xx=np.linspace(0,c["length"],201); yy=beam.basis(omega,xx)@coeff
        rqq,_=nsl.ritz_fields(v[:,i]*sign,c["length"],xx)
        for xxx,yyy,rr in zip(xx,yy,rqq):
            profiles.append({"position":i+1,"x":xxx,**dict(zip(("u","w","psi","N","Q","M"),yyy)),
                             **dict(zip(("ritz_u","ritz_w","ritz_psi"),rr))})
    gram=np.einsum("kni,ij,lnj,n->kl",np.array(modes),mass_matrix,np.array(modes),weights)
    spectral.update(mass_gram=gram,mass_gram_error=float(np.max(abs(gram-np.eye(7)))))
    write_csv(out/"frequencies.csv",rows);write_csv(out/"profiles.csv",profiles)
    if spectral["mass_gram_error"]>policy["mass_gram_tol"] or not all(all(g.values()) for g in allmetrics):
        raise ArithmeticError("Physical mode/comparison gate failed; inspect saved modes")
    spectral["status"]="PASS"


def identity(config):
    baseline=read(OUTPUT/"baseline.json")
    if any(sha(ROOT/n)!=h for n,h in baseline["sha256"].items()):
        raise ValueError("A protected historical source changed")
    if sha(ROOT/config["source_pdf"])!=config["source_sha256"]:
        raise ValueError("Bachoo PDF identity mismatch")
    paths=[CONFIG,Path(__file__),Path(nsl.__file__)]
    data={"config":config,"source_hashes":baseline["sha256"],"baseline_sha256":sha(OUTPUT/"baseline.json"),
          "files":{p.relative_to(ROOT).as_posix():sha(p) for p in paths},
          "versions":{"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__}}
    return hashlib.sha256(json.dumps(data,sort_keys=True).encode()).hexdigest()[:16],data


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--algebra-only",action="store_true")
    args=parser.parse_args()
    config=read(CONFIG);fp,ident=identity(config);out=OUTPUT/fp
    if (out/"manifest.json").exists():
        manifest=read(out/"manifest.json")
        if manifest["identity"]!=ident or any(sha(out/n)!=h for n,h in manifest["artifact_hashes"].items()):
            raise ValueError("NSL cache mismatch; no overwrite")
        saved=read(out/"result.json")
        if saved["status"]!="IN_PROGRESS" or not args.compute:
            print(json.dumps({"status":saved["status"],"cache":"VALIDATED","fingerprint":fp,"new_eigenvalue_solves":0}))
            return
    if not(args.compute or args.algebra_only):
        print(json.dumps({"sources":"PASS","fingerprint":fp,"new_eigenvalue_solves":0}));return
    out.mkdir(parents=True,exist_ok=True)
    result={"status":"IN_PROGRESS","parameters":config,"solve_ledger":{"ritz_eigenproblems":0,
            "propagation_prefix_searches":0,"propagation_roots":0,"reverse_root_searches":0,"complex_root_solves":0}}
    try:
        result["symbolic"]=symbolic_audit()
        prop,section,result["constitutive"]=coefficient_audit(config)
        result["properties"]=asdict(prop)
        write_json(out/"result.json",plain(result))
        if args.compute:
            spectral(config,prop,section,result,out)
            # Numerical/model gates pass, but printed Eq.(11d) is not reproduced.
            result["status"]="PARTIAL_PASS"
            result["qualification"]="All NSL model/numerical gates PASS; literal published Eq.(11d) FAIL. Local Eq.(9) determinant retained, no silent correction."
    except Exception as error:
        result["status"]="PARTIAL_PASS" if result.get("constitutive",{}).get("status")=="PASS" else "FAIL"
        result["stop_reason"]=f"{type(error).__name__}: {error}"
        if "spectral" in result: result["spectral"]["status"]="NUMERICAL_UNRESOLVED"
    result["source_preservation"]="PASS" if all(sha(ROOT/n)==h for n,h in ident["source_hashes"].items()) else "FAIL"
    write_json(out/"result.json",plain(result))
    write_json(out/"manifest.json",{"identity":ident,"command":subprocess.list2cmdline([sys.executable,*sys.argv]),
        "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "artifact_hashes":{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}})
    write_json(OUTPUT/"current.json",{"fingerprint":fp,"status":result["status"]})
    print(json.dumps({"status":result["status"],"fingerprint":fp,"solve_ledger":result["solve_ledger"],"stop_reason":result.get("stop_reason")}))
    if "stop_reason" in result:raise SystemExit(1)


if __name__=="__main__":
    main()
