"""Production-only geometry checks and sorted Lambda(beta) plot sweeps.

Per-arm heterogeneous composition calls verified arm basis/Dirichlet maps
and the SAME general joint operator. No wave/constitutive equation is copied.
Numerical predictor uses sorted frequencies only, followed by exact counts;
it assigns no modal identities. Previous physics helpers/bundles are immutable.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.linalg import block_diag, expm
from scipy.optimize import brentq
from scripts.analysis import screen_coupled_longitudinal_theory_hierarchy as hierarchy
from scripts.analysis import verify_mindlin_herrmann_timoshenko_general_beta_joint as general
from scripts.analysis import verify_mindlin_herrmann_timoshenko_single_rod as single
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section

CONFIG = ROOT/"data/input/mindlin_herrmann_timoshenko_lambda_beta_large_checks.json"
OUTPUT = ROOT/"results/mindlin_herrmann_timoshenko_lambda_beta_large_checks"
VERSION = "paired-production-arms-count-predictor-fixed-Lambda-v1"
POLICY = general.POLICY
PREFIX = "MHTIM_"


def lambda_factor(config):
    r, p = config["reference"], config["material"]
    A, I = r["b_ref"]*r["h_ref"], r["b_ref"]*r["h_ref"]**3/12
    return p["rho"]*A*r["l_ref"]**4/(p["E"]*I)


def lambda_from_frequency(frequency, config):
    f = np.asarray(frequency)
    if np.any(~np.isfinite(f)) or np.any(f < 0):
        raise ValueError("Finite nonnegative frequency required")
    return np.sqrt(2*math.pi*f)*lambda_factor(config)**.25


def geometry(config, family, parameter):
    r, p = config["reference"], config["material"]
    if family == "length_asymmetry":
        lengths = (r["l_ref"]*(1-parameter), r["l_ref"]*(1+parameter))
        heights = (r["h_ref"], r["h_ref"])
    elif family == "thickness_contrast":
        lengths = (r["l_ref"], r["l_ref"])
        heights = (r["h_ref"]*(1-parameter), r["h_ref"]*(1+parameter))
    else:
        raise ValueError("Explicit length_asymmetry or thickness_contrast required")
    models = tuple(mh.project_jang_reduced_rectangular(rectangular_section(E=p["E"],rho=p["rho"],nu=p["nu"],
        K=p["kappa"],width=r["b_ref"],thickness=h)) for h in heights)
    mass = sum(m.section.rhoA*l for m, l in zip(models, lengths))
    return models, lengths, {"family": family, "geometry_parameter": parameter,
        "lengths": lengths, "heights": heights, "sections": [vars(m.section) for m in models],
        "coefficients": [m.coefficients for m in models], "total_mass": mass}


def check_inputs():
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    base, single_config, model, length, checked, frozen, reference = hierarchy.check_inputs()
    if (config["material"] != base["material"] or config["reference"] != {"l_ref": .5, "b_ref": .2, "h_ref": .05}
        or config["beta_deg"] != np.arange(0, 90.1, 2.5).tolist() or config["length_mu"] != [0,.25,.5]
        or config["thickness_contrast"] != [0,.2,.4] or config["Lambda_definition"] != "Lambda^4 = rho*A_ref*omega^2*l_ref^4/(E*I_ref)"):
        raise ValueError("Declared geometry/reference/grid contract changed")
    # Canonical equation identity, with exact rational geometry.
    assert Fraction(1, 100)*Fraction(1, 2)**4/(Fraction(1, 5)*Fraction(1, 20)**3/12) == 300
    if not math.isclose(lambda_factor(config), 300, rel_tol=3e-15):
        raise ArithmeticError("Lambda reference scale failed")
    previous = ROOT/config["baseline_bundle"]
    manifest = json.loads((previous/"manifest.json").read_text(encoding="utf-8"))
    _, identity = hierarchy.identity(base, checked)
    if manifest["identity"] != identity or any(sha(previous/n) != h for n,h in manifest["artifact_hashes"].items()):
        raise ValueError("Preserved hierarchy baseline changed")
    return config, single_config, checked, previous


def identity(config, checked):
    paths = [CONFIG, Path(__file__), Path(mh.__file__), Path(joint.__file__),
        ROOT/"scripts/lib/reddy_inplane_geometry.py", Path(hierarchy.__file__), Path(general.__file__),
        Path(single.__file__), ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py"]
    record = {"version": VERSION, "config": config, "policy": POLICY,
        "files": {p.relative_to(ROOT).as_posix():sha(p) for p in paths}, "source_hashes": checked,
        "normalization_files": {n:sha(ROOT/n) for n in config["normalization_sources"]},
        "baseline_manifest": sha(ROOT/config["baseline_bundle"]/"manifest.json"),
        "versions": hierarchy.identity(json.loads(hierarchy.CONFIG.read_text(encoding="utf-8")), checked)[1]["versions"]}
    return hashlib.sha256(json.dumps(record,sort_keys=True).encode()).hexdigest()[:16], record


def boundary(models, lengths, omega, beta, straight=False):
    """Composition of verified state bases, no new local physics/sign rules."""
    matrix = np.zeros((16,16))
    if straight:
        # Both reference segment axes point in positive global X; no angles/maps.
        a0 = joint.arm_basis(models[0],lengths[0],omega,0.)
        a1 = joint.arm_basis(models[0],lengths[0],omega,lengths[0])
        b0 = joint.arm_basis(models[1],lengths[1],omega,0.)
        b1 = joint.arm_basis(models[1],lengths[1],omega,lengths[1])
        matrix[:4,:8], matrix[4:8,8:] = a0[:4], b1[:4]
        matrix[8:,:8], matrix[8:,8:] = a1, -b0
    else:
        ends = []
        for i,(model,length) in enumerate(zip(models,lengths)):
            matrix[i*4:(i+1)*4,i*8:(i+1)*8] = joint.arm_basis(model,length,omega,0.)[:4]
            ends.append(joint.arm_basis(model,length,omega,length))
        matrix[8:] = joint.joint_matrix(joint.frames(beta))@block_diag(*ends)
    return matrix


def scaled_boundary(models,lengths,omega,beta,straight=False):
    matrix = boundary(models,lengths,omega,beta,straight)
    return matrix/np.linalg.norm(matrix,axis=1)[:,None]


def pole_catalog(model,length,ceiling,base,perf):
    """Existing bounded basis/min-max certificates, independently per arm."""
    arm = model.section
    cutoff = mh.blocks(model)[1].cutoff_hz
    fraction = .9
    k = mh.blocks(model)[1].spatial(min(ceiling,.99*cutoff))[0]["wavenumber_per_m"]
    index = math.ceil(k*length/math.pi-fraction)+fraction
    tim_upper = mh.blocks(model)[1].temporal(index*math.pi/length)[0]["frequency_hz"]
    if tim_upper >= .99*cutoff:
        k = mh.blocks(model)[1].spatial(.99*cutoff)[0]["wavenumber_per_m"]
        index = math.floor(k*length/math.pi-fraction)+fraction
        tim_upper = mh.blocks(model)[1].temporal(index*math.pi/length)[0]["frequency_hz"]
    lower_cut = mh.finite_count_upper_bound(model,length,0.,"mh")["lower_contraction_cutoff_hz"]
    spacing = math.sqrt(arm.E/arm.rho)/(2*length)
    mh_upper = (math.ceil(ceiling/spacing-.5)+.5)*spacing
    if mh_upper >= .99*lower_cut:
        mh_upper = (math.floor(.99*lower_cut/spacing-.5)+.5)*spacing
    catalog = {}
    for block,upper in (("mh",mh_upper),("timoshenko",tim_upper)):
        count = mh.finite_count_upper_bound(model,length,2*math.pi*upper,block)
        attempts = []
        for retry in range(2):
            policy = {**base["policy"],"scan_intervals":400*2**retry}
            roots,search = mh.finite_roots(model,length,block,.001*math.pi,2*math.pi*upper,policy)
            attempts.append(search); perf["catalog_evaluations"] += search["evaluations"]
            if len(roots) == count["upper_count"]:
                break
        if len(roots) != count["upper_count"]:
            raise ArithmeticError(f"Arm catalog not saturated: {length},h={arm.thickness},{block}, {count}, attempts={attempts}")
        catalog[block] = {"roots":roots,"search":search,"count_bound":count,"attempts":attempts}
    return catalog


def count_query(models,lengths,omega,beta,catalogs,refmodel,perf):
    perf["count_evaluations"] += 1
    poles = sorted(r["omega"] for c in catalogs for block in c.values() for r in block["roots"])
    w = float(omega)
    exclusion = POLICY["pole_exclusion_relative"]*max(1.,w)
    excluded = []
    for pole in poles:
        if abs(w-pole) <= exclusion:
            excluded.append(pole); w = pole+2*exclusion
    schur, conditions = np.zeros((4,4)), []
    for model,length,frame in zip(models,lengths,joint.frames(beta)):
        stiffness,condition = joint.arm_dynamic_stiffness(model,length,w)
        schur += frame.nodal_transform@stiffness@frame.nodal_transform.T
        conditions.append(condition)
    p, total = refmodel.coefficients,sum(lengths)
    scales = np.array([math.sqrt(total/p["C"]),(p["C"]*p["H"])**(-.25),math.sqrt(total**3/p["B"]),math.sqrt(total/p["B"])])
    balanced = schur*scales[:,None]*scales[None,:]
    skew = float(np.max(abs(balanced-balanced.T))/max(np.linalg.norm(balanced,2),1e-30))
    eig = np.linalg.eigvalsh((balanced+balanced.T)/2)
    margin = float(min(abs(eig))/max(abs(eig)))
    if skew > POLICY["schur_symmetry_tol"] or margin < POLICY["count_inertia_margin_min"] or max(conditions)>POLICY["pole_condition_max"]:
        raise ArithmeticError(f"Count conditioning unresolved: skew={skew},margin={margin},arm={conditions}")
    j0,negative = sum(p<w for p in poles),int(np.count_nonzero(eig<0))
    return w,j0+negative,{"requested_omega":omega,"effective_omega":w,"excluded_poles":excluded,
        "J0":j0,"negative_inertia":negative,"scaled_skew_residual":skew,
        "inertia_relative_margin":margin,"arm_boundary_conditions":conditions}


def localized_roots(models,lengths,beta,catalogs,refmodel,upper,previous,config,perf):
    """Existing count-guided architecture with a numerical sorted-root seed."""
    lower = .001*math.pi
    attempts = []
    def run(seed):
        samples,brackets,failed = [],[],[]
        subdivisions = 0
        def det(w):
            perf["determinant_evaluations"] += 1
            return float(np.linalg.det(scaled_boundary(models,lengths,w,beta)))
        def query(w):
            effective,value,diagnostic = count_query(models,lengths,w,beta,catalogs,refmodel,perf)
            samples.append({"count":value,**diagnostic})
            return effective,value
        def interval(a,b,ca,cb,depth):
            nonlocal subdivisions
            if cb<ca:
                raise ArithmeticError("Count not monotone")
            if cb==ca or ca>=13:
                return
            if cb-ca==1:
                fa,fb = det(a),det(b)
                if fa*fb<=0:
                    brackets.append((a,b,fa,fb,ca,cb));return
            if depth>=POLICY["max_subdivision_depth"] or subdivisions>=POLICY["max_subdivisions"]:
                failed.append({"bracket_omega":[a,b],"counts":[ca,cb],"reason":"bounded subdivision exhausted"});return
            subdivisions+=1
            middle,cm = query((a+b)/2)
            if not a<middle<b:
                raise ArithmeticError("Pole exclusion crosses local interval")
            interval(a,middle,ca,cm,depth+1);interval(middle,b,cm,cb,depth+1)
        if seed is None:
            nodes = np.linspace(lower,upper,401)
        else:
            roots = [r["omega"] for r in seed["roots"][:13]]
            end = min(upper,roots[-1]+config["predictor_last_gap_padding"]*(roots[-1]-roots[-2]))
            nodes = np.array([lower,*[(a+b)/2 for a,b in zip(roots[:-1],roots[1:])],end])
        queries = [query(float(w)) for w in nodes]
        if queries[-1][1]<13:
            raise ArithmeticError("INCOMPLETE_ROOT_INVENTORY: predictor/ceiling does not cover guard13")
        if queries[0][1]!=0 or any(b[0]<=a[0] for a,b in zip(queries[:-1],queries[1:])):
            raise ArithmeticError("Invalid lower count or shifted query ordering")
        for (a,ca),(b,cb) in zip(queries[:-1],queries[1:]):
            interval(a,b,ca,cb,0)
        records = []
        for a,b,fa,fb,ca,cb in brackets:
            w,info = brentq(det,a,b,xtol=POLICY["root_xtol"],rtol=POLICY["root_rtol"],full_output=True)
            sv = np.linalg.svd(scaled_boundary(models,lengths,w,beta),compute_uv=False)
            records.append({"omega":w,"frequency_hz":w/(2*math.pi),"sorted_position":ca+1,
                "bracket_omega":[a,b],"bracket_counts":[ca,cb],"bracket_determinants":[fa,fb],
                "iterations":info.iterations,"sigma_min":float(sv[-1]),"singular_ratio":float(sv[-1]/sv[0]),
                "nonzero_singular_condition":float(sv[0]/sv[-2])})
        records.sort(key=lambda r:r["omega"])
        if failed or [r["sorted_position"] for r in records]!=list(range(1,14)):
            raise ArithmeticError(f"Prefix not certified: failed={failed},counts={[r['sorted_position'] for r in records]}")
        if any(b["omega"]-a["omega"]<=POLICY["root_xtol"] for a,b in zip(records[:-1],records[1:])):
            raise ArithmeticError("Duplicate/multiple root unresolved")
        perf["subdivisions"] += subdivisions
        return records,{"method":"full_seed_scan" if seed is None else "sorted_frequency_bracket_predictor",
            "count_samples":samples,"subdivisions":subdivisions,"failed_intervals":failed,
            "lower_count":0,"guard_right_count":13,"prefix_certificate":"PASS",
            "certified_range_omega":[lower,records[-1]["bracket_omega"][1]],"attempts":attempts}
    if previous is not None:
        perf["local_continuations"] += 1
        try:
            return run(previous)
        except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
            attempts.append({"method":"local_predictor","status":"FAILED_LOCALIZATION","reason":str(exc)})
            perf["fallback_full_scans"] += 1
    perf["full_scans"] += 1
    return run(None)


def mode(models,lengths,beta,omega,straight=False):
    matrix = scaled_boundary(models,lengths,omega,beta,straight)
    _,sv,right = np.linalg.svd(matrix)
    coefficients = right[-1].reshape(2,8)
    nodes,weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    mass,energy,clamp,qpeak,ppeak,pde = 0.,0.,0.,0.,0.,0.
    ends = []
    for i,(model,length,a) in enumerate(zip(models,lengths,coefficients)):
        x,weight = (nodes+1)*length/2,weights*length/2
        v = joint.arm_state(model,length,omega,a,x)
        dx = joint.arm_state(model,length,omega,a,x,1)
        end = joint.arm_state(model,length,omega,a,[0.,length])
        p,nu = model.coefficients,model.section.nu
        mass += float(weight@(v[:,:4]**2@np.array([p["m"],p["j"],p["m"],p["r"]])))
        energy += float(weight@(p["C"]*(dx[:,0]**2+2*nu*dx[:,0]*v[:,1]+v[:,1]**2)+p["H"]*dx[:,1]**2+p["B"]*dx[:,3]**2+p["S"]*(dx[:,2]-v[:,3])**2))
        qpeak = max(qpeak,float(np.max(abs(v[:,:4]))))
        ppeak = max(ppeak,float(np.max(abs(v[:,4:]))))  # total L=1, same work scale as old gate
        clamp = max(clamp,float(np.max(abs(end[1 if straight and i==1 else 0,:4]))))
        ends.append(end[0 if straight and i==1 else 1])
        rhs = v@mh.full_harmonic_state_matrix(model,omega).T
        scale = np.maximum(np.max(abs(rhs),axis=0)+np.max(abs(dx),axis=0),1e-30)
        pde = max(pde,float(np.max(abs(dx-rhs)/scale)))
    residual = ends[0]-ends[1] if straight else joint.joint_residual(*ends,joint.frames(beta))
    row_names = joint.STATE_ORDER if straight else joint.JOINT_ROWS
    scaled = abs(residual)/np.array([qpeak]*4+[ppeak]*4)
    diagnostics = {"mass_norm":mass/math.sqrt(mass)**2,"energy_relative_error":abs(energy/mass/omega**2-1),
        "clamp_scaled_residual":clamp/qpeak,"equation_scaled_residual":pde,
        "joint_residual_scaled":dict(zip(row_names,scaled.tolist())),"singular_ratio":float(sv[-1]/sv[0]),
        "nonzero_singular_condition":float(sv[0]/sv[-2])}
    if (max(scaled)>POLICY["boundary_joint_scaled_tol"] or diagnostics["clamp_scaled_residual"]>POLICY["boundary_joint_scaled_tol"]
        or pde>POLICY["equation_scaled_tol"] or diagnostics["energy_relative_error"]>POLICY["energy_relative_tol"]
        or diagnostics["singular_ratio"]>POLICY["boundary_joint_scaled_tol"] or diagnostics["nonzero_singular_condition"]>POLICY["nonzero_singular_condition_max"]):
        raise ArithmeticError(f"Root quality failed: {beta},omega={omega},diagnostics={diagnostics}")
    return {"coefficients":(coefficients/math.sqrt(mass)).tolist(),"diagnostics":diagnostics}


def solve_case(models,lengths,beta,catalogs,refmodel,usable,previous,config,perf):
    roots,search = localized_roots(models,lengths,beta,catalogs,refmodel,2*math.pi*usable,previous,config,perf)
    for r in roots:
        r.update(mode(models,lengths,beta,r["omega"]))
        r["f_star"] = r["frequency_hz"]*config["total_length"]/math.sqrt(config["material"]["E"]/config["material"]["rho"])
        r["Lambda"] = float(lambda_from_frequency(r["frequency_hz"],config))
    nodes,weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    gram = np.zeros((13,13))
    for i,(model,length) in enumerate(zip(models,lengths)):
        x,weight = (nodes+1)*length/2,weights*length/2
        fields = np.array([joint.arm_state(model,length,r["omega"],np.array(r["coefficients"])[i],x)[:,:4] for r in roots])
        masses = np.array([model.coefficients[k] for k in ("m","j","m","r")])
        balanced = (fields*np.sqrt(weight[:,None]*masses)).reshape(13,-1)
        gram += balanced@balanced.T
    error = float(np.max(abs(gram-np.eye(13))))
    if error>POLICY["mass_gram_tol"]:
        raise ArithmeticError(f"Frame mass normalization/orthogonality failed: {error}")
    return {"beta_deg":beta,"lengths":lengths,"roots":roots,"search":search,"mass_gram_max_error":error,"status":"PASS"}


def independent_endpoint(model,length,omega,perf):
    """Existing harmonic state/positive QR method; retain forces at end."""
    result = np.zeros((8,4));steps_record = []
    for block,cols in (("mh",slice(0,2)),("timoshenko",slice(2,4))):
        operator = mh.harmonic_state_matrix(model,omega,block)
        rate = float(max(abs(np.linalg.eigvals(operator))))
        p = model.coefficients
        elastic,gradient = (p["C"],p["H"]) if block=="mh" else (p["S"],p["B"])
        scales = np.array([1.,length,1/(elastic*rate),length/(gradient*rate)])
        balanced = operator*scales[:,None]/scales[None,:]
        steps = max(1,math.ceil(rate*length))
        if steps>512:
            raise ArithmeticError("Independent QR step budget exceeded")
        step = expm(balanced*length/steps)
        frame = np.vstack((np.zeros((2,2)),np.eye(2)))
        for _ in range(steps):
            frame,triangular = np.linalg.qr(step@frame,mode="reduced")
            frame *= np.where(np.diag(triangular)>=0,1.,-1.)[None,:]
        result[np.ix_(joint.BLOCK_INDICES[block],range(cols.start,cols.stop))] = frame/scales[:,None]
        steps_record.append(steps)
    return result,steps_record


def review_anomaly(models,lengths,case,position,perf):
    root = case["roots"][position-1]
    def matrix(w):
        perf["anomaly_independent_evaluations"] += 1
        ends = [independent_endpoint(m,l,w,perf) for m,l in zip(models,lengths)]
        raw = joint.joint_matrix(joint.frames(case["beta_deg"]))@block_diag(*(e[0] for e in ends))
        return raw/np.linalg.norm(raw,axis=1)[:,None],[e[1] for e in ends]
    def det(w):
        return float(np.linalg.det(matrix(w)[0]))
    a,b = root["bracket_omega"]
    if det(a)*det(b)>0:
        raise ArithmeticError("Independent anomaly determinant does not bracket root")
    omega = brentq(det,a,b,xtol=POLICY["root_xtol"],rtol=POLICY["root_rtol"])
    independent,steps = matrix(omega)
    sv = np.linalg.svd(independent,compute_uv=False)
    error = abs(omega/root["omega"]-1)
    if error>POLICY["beta0_frequency_relative_tol"] or sv[-1]/sv[0]>POLICY["boundary_joint_scaled_tol"]:
        raise ArithmeticError("Independent anomaly/root audit failed")
    return {"relative_frequency_difference":error,"boundary_singular_ratio":float(sv[-1]/sv[0]),
        "steps_by_arm_block":steps,"method":"same-angle independent state-expm/positive QR; no transfer-matrix product"}


def profiles_compare(models,lengths,first,second,swapped=False,straight=False):
    rows = []
    nodes,weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    for a,b in zip(first["roots"],second["roots"]):
        vfields,rfields,weightfields = [],[],[]
        for i,(model,length) in enumerate(zip(models,lengths)):
            x,weight = (nodes+1)*length/2,weights*length/2
            v = joint.arm_state(model,length,a["omega"],np.array(a["coefficients"])[i],x)
            j = 1-i if swapped else i
            r = joint.arm_state(model,length,b["omega"],np.array(b["coefficients"])[j],length-x if straight and i==1 else x)
            if swapped:
                r *= joint.MIRROR_STATE
            if straight and i==1:
                r *= joint.REFLECTION
            vfields.append(v);rfields.append(r);weightfields.append(weight)
        masses = [np.array([m.coefficients[k] for k in ("m","j","m","r")]) for m in models]
        overlap = sum(float(w@((v[:,:4]*r[:,:4])@mass)) for w,v,r,mass in zip(weightfields,vfields,rfields,masses))
        sign = 1. if overlap>=0 else -1.
        errors = []
        for ids in (slice(0,4),slice(4,8)):
            numerator = sum(float(w@np.sum((v[:,ids]-sign*r[:,ids])**2,axis=1)) for w,v,r in zip(weightfields,vfields,rfields))
            denominator = sum(float(w@np.sum(v[:,ids]**2,axis=1)) for w,v in zip(weightfields,vfields))
            errors.append(math.sqrt(numerator/denominator))
        frequency = abs(a["omega"]/b["omega"]-1)
        if frequency>POLICY["symmetry_frequency_relative_tol"] or max(errors)>POLICY["symmetry_component_L2_tol"]:
            raise ArithmeticError(f"Reference/swap profile gate failed: k={a['sorted_position']},freq={frequency},L2={errors}")
        rows.append({"position":a["sorted_position"],"frequency_relative_difference":frequency,"kinematic_L2_relative":errors[0],"resultant_L2_relative":errors[1]})
    return {"status":"PASS","rows":rows}


def straight_reference(models,lengths,primary,perf):
    roots = []
    def determinant(w):
        perf["reference_determinant_evaluations"] += 1
        return float(np.linalg.det(scaled_boundary(models,lengths,w,0,True)))
    for r in primary["roots"]:
        a,b = r["bracket_omega"]
        if determinant(a)*determinant(b)>0:
            raise ArithmeticError("Straight matching determinant does not bracket root")
        w = brentq(determinant,a,b,xtol=POLICY["root_xtol"],rtol=POLICY["root_rtol"])
        roots.append({"omega":w,"sorted_position":r["sorted_position"],**mode(models,lengths,0,w,True)})
    reference = {"roots":roots,"method":"two positive-X segment states; eight literal state-continuity rows; no angle transformation"}
    return reference,profiles_compare(models,lengths,primary,reference,straight=True)


def anomaly_flags(cases,config):
    flags = []
    values = np.array([[r["Lambda"] for r in c["roots"][:12]] for c in cases])
    for k in range(12):
        curvature = np.diff(values[:,k],n=2)
        absolute = abs(curvature)
        median = float(np.median(absolute));mad = float(np.median(abs(absolute-median)))
        threshold = median+config["anomaly_MAD_multiplier"]*mad+128*np.finfo(float).eps*max(abs(values[:,k]))
        for i in np.flatnonzero(absolute>threshold):
            r = cases[i+1]["roots"][k]
            flags.append({"beta_deg":cases[i+1]["beta_deg"],"position":k+1,"second_difference":float(curvature[i]),
                "diagnostic_threshold":threshold,"status":"REVIEWED_COUNT_SVD_RESIDUAL_PASS",
                "bracket_omega":r["bracket_omega"],"bracket_counts":r["bracket_counts"],"diagnostics":r["diagnostics"],
                "interpretation":"numerical candidate only; certified root retained without smoothing or physical-extremum claim"})
    return flags


def compute(config,base,previous,out):
    start = time.perf_counter()
    perf = {k:0 for k in ("determinant_evaluations","count_evaluations","catalog_evaluations","reference_determinant_evaluations",
        "full_scans","local_continuations","fallback_full_scans","subdivisions","anomaly_independent_evaluations","ceiling_expansions")}
    refmodels,_,_ = geometry(config,"length_asymmetry",0)
    cache = {}
    def catalogs_for(models,lengths,ceiling=None):
        ceiling = config["frequency_ceiling_fstar"] if ceiling is None else ceiling
        records = []
        for m,l in zip(models,lengths):
            key = (m.section.thickness,l,ceiling)
            if key not in cache:
                cache[key] = pole_catalog(m,l,ceiling,base,perf)
            records.append(cache[key])
        ceiling = min(ceiling,.98*min(c["search"]["range_omega"][1]/(2*math.pi) for arm in records for c in arm.values()))
        return records,ceiling
    def solve_controlled(models,lengths,beta,catalogs,usable,previous):
        attempts = [{"configured_ceiling_fstar":config["frequency_ceiling_fstar"],"usable_ceiling_fstar":usable}]
        try:
            case = solve_case(models,lengths,beta,catalogs,refmodels[0],usable,previous,config,perf)
        except ArithmeticError as exc:
            if "INCOMPLETE_ROOT_INVENTORY" not in str(exc):
                raise
            expanded = config["frequency_ceiling_fstar"]*config["one_ceiling_expansion_factor"]
            next_catalogs,next_usable = catalogs_for(models,lengths,expanded)
            attempts.append({"configured_ceiling_fstar":expanded,"usable_ceiling_fstar":next_usable,"reason":str(exc)})
            perf["ceiling_expansions"] += 1
            if next_usable<=usable:
                raise ArithmeticError(f"INCOMPLETE_ROOT_INVENTORY after one domain-limited expansion: {attempts}")
            case = solve_case(models,lengths,beta,next_catalogs,refmodels[0],next_usable,None,config,perf)
        case["ceiling_attempts"] = attempts
        return case
    result = {"geometries":[],"series":{},"reference_checks":[],"swap_checks":[],"baseline_regression":[]}
    geometries = [("length_asymmetry",p) for p in config["length_mu"]]+[("thickness_contrast",p) for p in config["thickness_contrast"] if p!=0]
    for family,parameter in geometries:
        models,lengths,description = geometry(config,family,parameter)
        catalogs,usable = catalogs_for(models,lengths)
        key = f"{family}_{parameter:g}"
        result["geometries"].append(description)
        series = []
        for beta in config["beta_deg"]:
            case = solve_controlled(models,lengths,beta,catalogs,usable,series[-1] if series else None)
            series.append(case)
            if beta==0:
                ref,check = straight_reference(models,lengths,case,perf)
                if family=="length_asymmetry":
                    direct = json.loads((previous/"inventory_mindlin_herrmann_0.json").read_text(encoding="utf-8"))
                    diffs = [abs(a["omega"]/b["omega"]-1) for a,b in zip(case["roots"],direct["roots"])]
                    if max(diffs)>POLICY["beta0_frequency_relative_tol"]:
                        raise ArithmeticError("Length beta0 collapse failed")
                    check["direct_homogeneous_max_frequency_relative_difference"] = max(diffs)
                result["reference_checks"].append({"family":family,"geometry_parameter":parameter,"reference":ref,"comparison":check})
            if parameter==0 and beta in (0,5,15,30,45,60,75,90):
                preserved = json.loads((previous/f"inventory_mindlin_herrmann_{int(beta)}.json").read_text(encoding="utf-8"))
                differences = [abs(a["omega"]/b["omega"]-1) for a,b in zip(case["roots"],preserved["roots"])]
                if max(differences)>POLICY["beta0_frequency_relative_tol"]:
                    raise ArithmeticError("Preserved baseline frequency regression failed")
                result["baseline_regression"].append({"beta_deg":beta,"frequency_relative_differences":differences,
                    "Lambda_relative_differences":[abs(a["Lambda"]/float(lambda_from_frequency(b["frequency_hz"],config))-1) for a,b in zip(case["roots"],preserved["roots"])]})
            if beta%15==0:
                print(f"{key}: beta={beta:g}, prefix13 certified, predictor={case['search']['method']}",flush=True)
        flags = anomaly_flags(series,config)
        for flag in flags:
            case = next(c for c in series if c["beta_deg"]==flag["beta_deg"])
            flag["independent_same_point_review"] = review_anomaly(models,lengths,case,flag["position"],perf)
            flag["status"] = "REVIEWED_INDEPENDENT_QR_COUNT_SVD_RESIDUAL_PASS"
        result["series"][key] = {"geometry":description,"cases":series,"anomaly_flags":flags,"usable_ceiling_fstar":usable}
        write_json(out/f"series_{key}.json",result["series"][key])
        # Same physical system: canonical negative parameter swaps via reflection
        # in the angle bisector. Local w/theta/Q/M reverse; c/R do not.
        if parameter:
            swapped_models,swapped_lengths,_ = geometry(config,family,-parameter)
            swapped_catalogs,swapped_usable = catalogs_for(swapped_models,swapped_lengths)
            for beta in config["symmetry_beta_deg"]:
                canonical = next(c for c in series if c["beta_deg"]==beta)
                swapped = solve_controlled(swapped_models,swapped_lengths,beta,swapped_catalogs,swapped_usable,canonical)
                comparison = profiles_compare(models,lengths,canonical,swapped,swapped=True)
                result["swap_checks"].append({"family":family,"geometry_parameter":parameter,"beta_deg":beta,"comparison":comparison,"negative_parameter_case":swapped})
        write_json(out/"gate_checks.json",{k:result[k] for k in ("reference_checks","swap_checks","baseline_regression")})
    result["pole_catalogs"] = [{"h":key[0],"length":key[1],"configured_ceiling_fstar":key[2],"catalog":value} for key,value in cache.items()]
    statuses = {PREFIX+name:"PASS" for name in ("LAMBDA_NORMALIZATION","LAMBDA_MAP_BASELINE_REGRESSION","LENGTH_BETA0_COLLAPSE","LENGTH_ARM_SWAP",
        "LENGTH_LAMBDA_BETA_MAP","THICKNESS_STEPPED_BETA0","THICKNESS_ARM_SWAP","THICKNESS_LAMBDA_BETA_MAP","LAMBDA_BETA_ROOT_QUALITY","LAMBDA_BETA_LARGE_CHECKS")}
    result.update({"statuses":statuses,"normalization":{"definition":config["Lambda_definition"],"factor_Lambda4_over_omega2":lambda_factor(config),
        "reference":config["reference"],"mapping":"Lambda=sqrt(omega)*300^(1/4)=sqrt(2*pi*f_star)*300^(1/4) for E=rho=L=1; RLB-2B local Lambda is canonical Lambda^2"},
        "performance":{**perf,"compute_runtime_seconds":time.perf_counter()-start},"spectrum_semantics":"independently sorted positions; no shape continuation"})
    return result


def plot_saved(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    result = json.loads((out/"result.json").read_text(encoding="utf-8"))
    config = json.loads((out/"parameters.json").read_text(encoding="utf-8"))["config"]
    with plt.rc_context({"font.family":"DejaVu Serif","font.size":9,"axes.labelsize":10,"pdf.fonttype":42}):
        values = [r["Lambda"] for s in result["series"].values() for c in s["cases"] for r in c["roots"][:12]]
        ylim = (0,max(values)*1.04)
        colors = [plt.get_cmap("tab20")(i) for i in (0,2,4,6,8,10,12,14,16,18,1,7)]
        styles = ["-","--","-."]
        for family,parameters,symbol in (("length_asymmetry",config["length_mu"],r"\mu"),("thickness_contrast",config["thickness_contrast"],r"\delta_h")):
            fig,axes = plt.subplots(1,3,figsize=(10,4),sharex=True,sharey=True)
            for n,(ax,parameter) in enumerate(zip(axes,parameters)):
                key = f"{family}_{parameter:g}" if parameter else "length_asymmetry_0"
                series = result["series"][key]["cases"]
                for k in range(12):
                    ax.plot([c["beta_deg"] for c in series],[c["roots"][k]["Lambda"] for c in series],color=colors[k],ls=styles[k//4],lw=1.05)
                ax.set(xlim=(0,90),ylim=ylim,xticks=(0,30,60,90),title=rf"({chr(97+n)}) ${symbol}={parameter:.2f}$")
                ax.grid(alpha=.16)
            axes[0].set_ylabel(r"$\Lambda$")
            fig.supxlabel(r"$\beta$, deg",y=.16)
            handles = [Line2D([0],[0],color=colors[k],ls=styles[k//4],lw=1.05,label=rf"$k={k+1}$") for k in range(12)]
            fig.legend(handles=handles,loc="lower center",ncol=6,frameon=False,bbox_to_anchor=(.5,-.01),handlelength=2.2,columnspacing=1.2)
            fig.subplots_adjust(left=.065,right=.995,top=.91,bottom=.26,wspace=.10)
            name = family+"_Lambda_beta"
            fig.savefig(out/(name+".png"),dpi=300)
            fig.savefig(out/(name+".pdf"),metadata={"Creator":"CoupledBeams","CreationDate":None,"ModDate":None})
            plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-sources",action="store_true")
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--plot-only",type=Path,metavar="BUNDLE")
    parser.add_argument("--output-dir",type=Path,default=OUTPUT)
    args = parser.parse_args(argv)
    if not(args.check_sources or args.compute or args.plot_only):
        parser.error("Select --check-sources, --compute or --plot-only BUNDLE")
    if args.plot_only:
        if args.compute:
            parser.error("plot-only cannot compute")
        manifest = json.loads((args.plot_only/"manifest.json").read_text(encoding="utf-8"))
        if any(sha(args.plot_only/n)!=h for n,h in manifest["artifact_hashes"].items()):
            raise ValueError("Saved map changed")
        plot_saved(args.plot_only)
        if any(sha(args.plot_only/n)!=h for n,h in manifest["artifact_hashes"].items()):
            raise ValueError("Plot-only bytes changed")
        print("Plot-only: zero root evaluations",args.plot_only);return 0
    config,base,checked,previous = check_inputs()
    fingerprint,inputs = identity(config,checked)
    if args.check_sources:
        print("Canonical Lambda and source/baseline hashes PASS; zero roots")
    if not args.compute:
        return 0
    out = args.output_dir/fingerprint
    if (out/"manifest.json").exists():
        manifest = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
        if manifest["identity"]!=inputs or any(sha(out/n)!=h for n,h in manifest["artifact_hashes"].items()):
            raise ValueError("Stale geometry map cache")
        print("Verified map cache: zero root evaluations",out);return 0
    out.mkdir(parents=True,exist_ok=True)
    write_json(out/"parameters.json",{"config":config,"policy":POLICY})
    try:
        result = compute(config,base,previous,out)
    except (ValueError,ArithmeticError,np.linalg.LinAlgError) as exc:
        write_json(out/"failure.json",{"status":"FAIL","reason":str(exc),"identity":inputs});raise
    write_json(out/"result.json",result)
    rows,summary = [],[]
    for family,parameters in (("length_asymmetry",config["length_mu"]),("thickness_contrast",config["thickness_contrast"])):
        for parameter in parameters:
            series = result["series"][f"{family}_{parameter:g}" if parameter else "length_asymmetry_0"]
            ds = [r["diagnostics"] for c in series["cases"] for r in c["roots"]]
            summary.append({"family":family,"geometry_parameter":parameter,"cases_passed":37,
                "max_residual":max(max(d["joint_residual_scaled"].values()) for d in ds),
                "max_condition":max(d["nonzero_singular_condition"] for d in ds),
                "anomaly_flags":len(series["anomaly_flags"]),"subdivisions":sum(c["search"]["subdivisions"] for c in series["cases"]),
                "fallback_scans":sum(bool(c["search"]["attempts"]) for c in series["cases"]),
                "mass_gram_max_error":max(c["mass_gram_max_error"] for c in series["cases"]),
                "beta0_reference_difference":max(r["frequency_relative_difference"] for check in result["reference_checks"] if check["family"]==("length_asymmetry" if parameter==0 else family) and check["geometry_parameter"]==parameter for r in check["comparison"]["rows"]),
                "arm_swap_difference":max((r["frequency_relative_difference"] for check in result["swap_checks"] if check["family"]==family and check["geometry_parameter"]==parameter for r in check["comparison"]["rows"]),default=0.)})
            for c in series["cases"]:
                for r in c["roots"]:
                    d = r["diagnostics"]
                    rows.append({"family":family,"geometry_parameter":parameter,"beta_deg":c["beta_deg"],"position":r["sorted_position"],
                        "frequency":r["frequency_hz"],"f_star":r["f_star"],"Lambda":r["Lambda"],"guard":r["sorted_position"]==13,
                        "bracket_left":r["bracket_omega"][0],"bracket_right":r["bracket_omega"][1],"sigma_min":r["sigma_min"],
                        "singular_ratio":r["singular_ratio"],"condition":r["nonzero_singular_condition"],
                        "clamp_residual":d["clamp_scaled_residual"],"PDE_residual":d["equation_scaled_residual"],
                        **d["joint_residual_scaled"],"completeness_status":"PASS"})
    write_csv(out/"Lambda_spectra.csv",rows);write_csv(out/"summary.csv",summary);write_json(out/"summary.json",summary)
    plot_saved(out)
    write_json(out/"manifest.json",{"identity":inputs,"statuses":result["statuses"],
        "artifact_hashes":{p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="manifest.json"},
        "git_branch":subprocess.check_output(["git","branch","--show-current"],text=True).strip(),
        "git_head":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
        "git_status":subprocess.check_output(["git","status","--short"],text=True,encoding="utf-8"),
        "command":subprocess.list2cmdline([sys.executable,*sys.argv])})
    write_json(args.output_dir/"current.json",{"fingerprint":fingerprint,"directory":str(out.resolve())})
    print("MHTIM_LAMBDA_BETA_LARGE_CHECKS PASS",out);return 0


if __name__=="__main__":
    raise SystemExit(main())
