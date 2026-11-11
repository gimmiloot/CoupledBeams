"""Read-only continuation of saved square RL/Timoshenko physical modes.

New workflow: unlike the a=4 modal audit and reflection-restricted renderer,
this compares unequal-length states using a documented common-domain inertia
product. Reuses analytic_branch_tracking.assignment_diagnostics only; no
production basis, spectrum, Ritz or root routine can be called. Diagnostic
results, not a modification of the archive contract or an accepted global ID.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import hashlib
import json
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
for _path in (ROOT, ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))
from scripts.lib import analytic_branch_tracking as tracker
from scripts.lib import frequency_geometry_style as line_geometry
from scripts.lib import rl_full_rr_sections_prefix as sections_engine
from scripts.lib.rl_full_rr_aspect import make_model

DEFAULT = ROOT / "private/dissertation/research/chapter1/chapter1_existing_square_postprocess.json"
MAPS = ROOT / "results/rl_full_rr_maps"
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def table(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def write_csv(path, rows):
    if not rows:
        Path(path).write_text("", encoding="utf-8")
        return
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


@contextmanager
def saved_only():
    """Fail closed on any accidental call into spectral/basis implementations."""
    prior = sys.getprofile()
    allowed_tracker = {"assignment_diagnostics", "mac_matrix", "shape_vector_from_components", "validate_shape_metric", "__init__"}
    # Canonical geometry/contract checks only. State matrices, bases, recovery,
    # root counters and Ritz remain forbidden, including in these modules.
    geometry_only = {
        "scripts.lib.rl_full_rr_aspect": {"make_model", "__init__", "__post_init__"},
        "scripts.lib.rl_full_rr_sections_prefix": {"point_contract", "arm_contract", "<dictcomp>", "<listcomp>"},
        "scripts.lib.isotropic_two_springs_thickness_asymmetry": {"__init__", "__post_init__", "ratios", "arms", "compatible", "<genexpr>"},
        "scripts.lib.isotropic_two_springs_rl_timoshenko": {"__init__", "__post_init__", "section", "J", "time", "lengths"},
        "scripts.lib.isotropic_rectangular_timoshenko_coupled_beams": set("rectangular_section __init__ __post_init__ _positive _finite _poisson_ratio EA EI KGA rhoA rhoI axial_stiffness bending_stiffness shear_stiffness shear_modulus mass_per_length rotary_inertia_per_length".split()),
    }
    def guard(frame, event, arg):
        if event != "call":
            return
        module = frame.f_globals.get("__name__", "")
        name = frame.f_code.co_name
        if module == "scripts.lib.analytic_branch_tracking":
            if name not in allowed_tracker and name not in {"<listcomp>", "<genexpr>"}:
                raise RuntimeError("Forbidden tracking/solver entry: "+name)
        elif module.startswith(("scripts.lib.isotropic_", "scripts.lib.rl_full_", "scripts.lib.mindlin_", "my_project.analytic.solvers", "my_project.analytic.formulas")):
            if name not in geometry_only.get(module, set()):
                raise RuntimeError("Forbidden physical evaluation: "+module+":"+name)
    sys.setprofile(guard)
    try:
        yield
    finally:
        sys.setprofile(prior)


def frames(beta):
    b = np.deg2rad(beta)
    return ((np.array([1., 0.]), np.array([0., -1.])),
            (np.array([-np.cos(b), -np.sin(b)]), np.array([-np.sin(b), np.cos(b)])))


def metric_between(left, right, fraction=.5):
    # Fixed-geometry beta paths have identical coefficients on each arm at
    # both ends, including tau != 0. Legacy mu paths are restricted to tau=0.
    return [{k: (1-fraction)*a[k]+fraction*b[k] for k in ("length", "rhoA", "rhoI", "J")}
            for a, b in zip(left["sections"], right["sections"])]


def features(point, metric=None):
    """Full inertia embedding after pullback x=L_source*xi.

    Integral on xi: Lbar*rhoA*d.d + Lbar*rhoI*psi*psi +
    Jbar/Lbar*U_xi*U_xi. Translations use each source's global frame;
    section rotations and material axial derivatives use corresponding arms.
    No interpolation: saved normalized Gauss nodes/weights must coincide.
    """
    metric = point["sections"] if metric is None else metric
    blocks = []
    for arm, sec, met, (t, n) in zip(point["arms"], point["sections"], metric, frames(point["beta"])):
        y, dy = arm["y"], arm["dy"]
        w = arm["weights"] / sec["length"]
        if np.any(w <= 0):
            raise ValueError("Nonpositive quadrature weight")
        d = y[:, :, 0, None]*t + y[:, :, 1, None]*n
        translations = (np.sqrt(w*met["length"]*met["rhoA"])[None, :, None]*d).reshape(len(y), -1)
        rotation = np.sqrt(w*met["length"]*met["rhoI"])[None, :]*y[:, :, 2]
        axial = np.sqrt(w*met["J"]/met["length"])[None, :]*sec["length"]*dy[:, :, 0]
        blocks.append(np.concatenate((translations, rotation, axial), axis=1))
    return np.concatenate(blocks, axis=1)


def normalized(a):
    norm = np.linalg.norm(a, axis=1)
    if np.any(~np.isfinite(norm)) or np.any(norm <= 0):
        raise ValueError("Invalid inertia norm")
    return a/norm[:, None]


def assignment(a, b, fa, fb, penalty=0.):
    # Four legacy component names are storage slots, NOT the old EB metric.
    prev = [SimpleNamespace(Lambda=f, components={"u_left": v, "w_left": np.empty(0),
            "u_right": np.empty(0), "w_right": np.empty(0)}) for v, f in zip(a, fa)]
    current = [SimpleNamespace(Lambda=f, shape_vector=v) for v, f in zip(normalized(b), fb)]
    return tracker.assignment_diagnostics(previous_points=prev, current_states=current,
        sorted_lambdas=fb, freq_weight=penalty, shape_metric="full")


def principal_cosines(a, b):
    # QR removes the effect of nonorthogonality in a cross-geometry metric.
    qa, ra = np.linalg.qr(a.T)
    qb, rb = np.linalg.qr(b.T)
    if np.linalg.matrix_rank(ra) != len(a) or np.linalg.matrix_rank(rb) != len(b):
        raise ValueError("Rank-deficient comparison subspace")
    return np.linalg.svd(qa.T @ qb, compute_uv=False)


def accept_transition(overlap, margin, stable, jump, threshold=.95, min_margin=.1):
    return bool(overlap >= threshold and margin >= min_margin and stable and jump <= 1)


class Archive:
    def __init__(self):
        self.hashes, self.cache, self.manifests = {}, {}, {}
        self.provenance_path = MAPS / "_audit_beta_mu_refined_20261107/point_provenance.csv"
        self.pin(self.provenance_path)
        self.provenance = {(r["map"], float(r["beta_deg"])): r for r in table(self.provenance_path)}

    def pin(self, path, expected=None):
        path = Path(path).resolve()
        key = path.relative_to(ROOT).as_posix()
        actual = sha(path)
        if expected is not None and actual != expected:
            raise ValueError("SHA256 mismatch: "+key)
        if key in self.hashes and self.hashes[key] != actual:
            raise ValueError("Input changed during run: "+key)
        self.hashes[key] = actual
        return path

    def manifest(self, path):
        path = self.pin(path)
        if path not in self.manifests:
            self.manifests[path] = read(path)
        m = self.manifests[path]
        if m.get("status") != "PASS":
            raise ValueError("Unaccepted source manifest")
        return m

    def point(self, row):
        path = self.pin(ROOT/row["record"], row["record_sha256"])
        if path in self.cache:
            return self.cache[path]
        r = read(path)
        c = r["contract"]
        expected = {"epsilon": .005, "aspect": 1., "E": 1., "rho": 1., "l": 1., "nu": .3, "ks": 5/6, "variant": "RL_FULL"}
        if any(c["model"].get(k) != v for k, v in expected.items()):
            raise ValueError("Different model contract")
        kp, kg = c["K_plot"], c["K_guard"]
        count = len(r["roots"])
        if (r["status"] != "PASS" or r["quality"]["status"] != "PASS" or (kp, kg) not in ((6,7),(8,9))
                or not kg <= count <= kg+c["policy"]["max_cluster_extra_ranks"]
                or c["case"]["k_joint"] != "RIGID" or c["case"]["k_end"] != "RIGID"
                or c["coordinates"] != "outer clamp to joint"
                or c["state_order"] != ["u", "w", "psi", "N", "Q", "M"]):
            raise ValueError("Incomplete/incompatible point")
        model = make_model(**c["model"])
        case = sections_engine.p.Case(**c["case"])
        canonical = sections_engine.point_contract(model, case, kp, kg, c["policy"])
        if (c["arms"] != canonical["arms"] or c["reference_section"] != canonical["reference_section"]
                or c["time_reference"] != canonical["time_reference"]
                or c["policy"] != read(ROOT/"data/input/rl_full_rr_geometry_c.json")["policy"]
                or r["tau"] != c["tau"] or r["tau"] != case.tau or r["mu"] != case.mu
                or r["beta_deg"] != case.beta
                or c["normalization"] not in ("Omega=omega*t_ref=Lambda**2", "Omega=omega*t_ref=Lambda**2; common reference section")):
            raise ValueError("Inconsistent canonical physical contract")
        if any(not all(g.values()) for g in r["quality"]["gates"]):
            raise ValueError("Failed saved physical check")
        if (r["search"].get("unresolved_intervals") or r["search"]["lower_count"] != 0
                or r["search"]["guard_right_count"] != count):
            raise ValueError("Unresolved saved spectrum")
        owner = path.parent.parent
        npz = self.pin(owner/r["shape_file"], r["shape_sha256"])
        # Some accepted points were salvaged from an interrupted producer.
        # The final map manifest pins the JSON, which pins the complete NPZ;
        # a missing historical producer manifest does not break that chain.
        if (owner/"manifest.json").exists():
            owner_manifest = self.manifest(owner/"manifest.json")
            for file in (path, npz):
                key = file.relative_to(owner).as_posix()
                if key in owner_manifest["artifacts"]:
                    self.pin(file, owner_manifest["artifacts"][key])
        arms = []
        with np.load(npz, allow_pickle=False) as z:
            for i in (1, 2):
                arms.append({k: np.stack([z[f"mode{j}_arm{i}_{k}"] for j in range(1, count+1)])
                             for k in ("y", "dy", "ends")})
                for key in ("x", "weights"):
                    arrays = [z[f"mode{j}_arm{i}_{key}"] for j in range(1, count+1)]
                    if any(not np.array_equal(arrays[0], q) for q in arrays[1:]):
                        raise ValueError("Inconsistent modal quadrature")
                    arms[-1][key] = arrays[0].copy()
        p = dict(beta=float(r["beta_deg"]), mu=float(r["mu"]), tau=float(r["tau"]), epsilon=model.epsilon,
                 K_plot=kp, K_guard=kg, sections=c["arms"], arms=arms,
                 frequencies=np.array([q["Lambda"] for q in r["roots"]]), record=row["record"],
                 npz=npz.relative_to(ROOT).as_posix(), json=r)
        for i, sec in enumerate(p["sections"]):
            np.testing.assert_allclose(sec["length"], 1+(-1 if i == 0 else 1)*p["mu"], atol=1e-14)
        f = features(p)
        gram = f @ f.T
        np.testing.assert_allclose(gram, r["quality"]["gram"], atol=2e-12, rtol=2e-12)
        if np.max(np.abs(gram-np.eye(count))) > c["policy"]["mass_gram_tol"]:
            raise ValueError("Saved full-inertia orthonormality failed")
        for root in r["roots"]:
            np.testing.assert_allclose(root["omega"], .005*root["Lambda"]**2, rtol=5e-15)
        p["mass_gram_error"] = float(np.max(np.abs(gram-np.eye(count))))
        for k, root in enumerate(r["roots"]):
            mass = []
            for arm, sec in zip(arms, c["arms"]):
                y, dy, w = arm["y"][k], arm["dy"][k], arm["weights"]
                mass.append([np.dot(w,sec["rhoA"]*y[:,0]**2),np.dot(w,sec["rhoA"]*y[:,1]**2),
                             np.dot(w,sec["rhoI"]*y[:,2]**2),np.dot(w,sec["J"]*dy[:,0]**2)])
            np.testing.assert_allclose(mass,root["diagnostics"]["mass_components"],atol=2e-13,rtol=2e-12)
        self.cache[path] = p
        return p

    def series(self, spec):
        cfg = read(self.pin(ROOT/spec["config"]))
        config_map = next(a for a in cfg["maps"] if a["name"] == spec["map"])
        man = self.manifest(ROOT/spec["manifest"])
        if spec["axis"] != {"BETA":"beta_deg","MU":"mu"}.get(config_map["mode"]):
            raise ValueError("Sweep axis differs from the physical map")
        if cfg["material"] != man["identity"]["material"]:
            raise ValueError("Configuration/manifest material mismatch")
        ownman = self.manifest((ROOT/spec["csv"]).parent/"manifest.json")
        self.pin(ROOT/spec["csv"], ownman["artifacts"]["spectrum.csv"])
        rows = [r for r in table(ROOT/spec["csv"]) if r["map"] == spec["map"]]
        points_meta = {(q["map"], float(q["x"])): q for q in man["points"]}
        if spec["axis"] == "mu":
            dense = (ROOT/spec["manifest"]).parent/"spectrum.csv"
            self.pin(dense, man["artifacts"]["spectrum.csv"])
            dense_rows = [r for r in table(dense) if r["map"] == spec["map"]]
            # The derivative individual bundle reports reuse where the parent
            # reports its original solve. This origin label is bundle-relative.
            strip_origin = lambda r: {k:v for k,v in r.items() if k != "origin"}
            if list(map(strip_origin, rows)) != list(map(strip_origin, dense_rows)):
                raise ValueError("Individual and common map CSV disagree")
        grid = sorted({float(r["x"]) for r in rows})
        if grid != config_map["grid"]:
            raise ValueError("Incomplete requested grid")
        data = []
        for x in grid:
            group = sorted([r for r in rows if float(r["x"]) == x], key=lambda r:int(r["rank"]))
            if [int(r["rank"]) for r in group] != list(range(1,len(group)+1)):
                raise ValueError("Incomplete guard prefix")
            p = self.point(group[0])
            if len(group) != len(p["frequencies"]) or p["K_guard"] != spec.get("K_guard",9):
                raise ValueError("Unexpected saved prefix")
            meta = points_meta[(spec["map"],x)]
            if meta["record"] != p["record"] or meta["record_sha256"] != group[0]["record_sha256"]:
                raise ValueError("Manifest/CSV provenance mismatch")
            for r in group:
                k = int(r["rank"])-1
                if (r["record"] != p["record"] or float(r["Lambda"]) != p["frequencies"][k]
                        or float(r["mu"]) != p["mu"] or float(r["beta_deg"]) != p["beta"]
                        or float(r["tau"]) != p["tau"] or r["status"] != "PASS"):
                    raise ValueError("Point/CSV mismatch")
            if spec["axis"] == "beta_deg" and spec.get("refinement_provenance",True):
                pr = self.provenance[(spec["map"], x)]
                if pr["record"] != p["record"] or pr["shape_file"] != p["npz"] or pr["shape_sha256"] != sha(ROOT/p["npz"]):
                    raise ValueError("Refinement provenance mismatch")
            for k, v in config_map["fixed"].items():
                if float(group[0][k]) != v:
                    raise ValueError("Fixed parameter mismatch")
            data.append(p)
        for p in data[1:]:
            if spec["axis"] == "beta_deg" and (p["sections"] != data[0]["sections"] or p["mu"] != data[0]["mu"] or p["tau"] != data[0]["tau"]):
                raise ValueError("Angular path changes its physical sections")
            if spec["axis"] == "mu" and (p["tau"] != 0 or data[0]["tau"] != 0):
                raise ValueError("Variable-section mu transport is outside this postprocessor contract")
            for a,b,sa,sb in zip(data[0]["arms"],p["arms"],data[0]["sections"],p["sections"]):
                np.testing.assert_allclose(a["x"]/sa["length"], b["x"]/sb["length"], atol=3e-15, rtol=0)
                np.testing.assert_allclose(a["weights"]/sa["length"], b["weights"]/sb["length"], atol=3e-15, rtol=0)
        return data


def seed_ids(spec, data, archive):
    n = spec["branches"]
    if spec["axis"] == "beta_deg":
        return [(f"B{i+1}", "", 0) for i in range(n)]
    folder = MAPS/"_audit_reflection_classification_20261107"
    man = read(archive.pin(folder/"manifest.json"))
    path = folder/"class_spectrum.csv"
    expected = man.get("outputs", man.get("artifacts", {})).get("class_spectrum.csv")
    archive.pin(path, expected)
    rows = [r for r in table(path) if float(r["beta_deg"]) == data[0]["beta"]]
    result = []
    for k in range(1,n+1):
        r = next(r for r in rows if int(r["original_sorted_position"]) == k)
        if float(r["Lambda"]) != data[0]["frequencies"][k-1] or r["classification_status"] not in ("PURE_PLUS","PURE_MINUS"):
            raise ValueError("Unconfirmed seed class")
        rank, cls = int(r["class_rank"]), r["symmetry_class"]
        result.append((f"B{2*rank-(1 if cls == '+' else 0)}", cls, rank))
    return result


def characterize(point, k, floors):
    mass, stiffness, nodes = [], [], []
    peak = max(float(np.max(np.linalg.norm(a["y"][k,:,:2],axis=1))) for a in point["arms"])
    for a,s in zip(point["arms"],point["sections"]):
        y,dy,w = a["y"][k],a["dy"][k],a["weights"]
        mass.append([float(np.dot(w,s["rhoA"]*y[:,0]**2)),float(np.dot(w,s["rhoA"]*y[:,1]**2)),
                     float(np.dot(w,s["rhoI"]*y[:,2]**2)),float(np.dot(w,s["J"]*dy[:,0]**2))])
        stiffness.append(float(np.dot(w,s["EA"]*dy[:,0]**2+s["EI"]*dy[:,2]**2+s["KGA"]*(dy[:,1]-y[:,2])**2)))
        counts, brackets = [], []
        for floor in floors:
            selected=np.flatnonzero(np.abs(y[:,1])>floor*peak)
            cuts=[(int(i),int(j)) for i,j in zip(selected[:-1],selected[1:]) if y[i,1]*y[j,1]<0]
            counts.append(len(cuts) if len(selected) else -1)
            brackets.append([[float(a["x"][i]/s["length"]),float(a["x"][j]/s["length"])] for i,j in cuts])
        nodes.append(dict(counts=counts, stable=len(set(counts))==1 and counts[0]>=0, brackets=brackets[1]))
    axial_nodes=[]
    for a,s in zip(point["arms"],point["sections"]):
        counts=[];brackets=[]
        for floor in floors:
            selected=np.flatnonzero(np.abs(a["y"][k,:,0])>floor*peak)
            cuts=[(int(i),int(j)) for i,j in zip(selected[:-1],selected[1:]) if a["y"][k,i,0]*a["y"][k,j,0]<0]
            counts.append(len(cuts) if len(selected) else -1)
            brackets.append([[float(a["x"][i]/s["length"]),float(a["x"][j]/s["length"])] for i,j in cuts])
        axial_nodes.append(dict(counts=counts,stable=len(set(counts))==1 and counts[0]>=0,brackets=brackets[1]))
    m=np.array(mass); total=float(m.sum());omega=point.get("epsilon",.005)*point["frequencies"][k]**2
    return dict(inertia_arm1=float(m[0].sum()/total), inertia_arm2=float(m[1].sum()/total),
                inertia_axial=float(m[:,[0,3]].sum()/total), inertia_rotation=float(m[:,2].sum()/total),
                inertia_transverse=float(m[:,1].sum()/total),inertia_components=json.dumps(m.tolist()),
                mass=total, energy_error=abs(sum(stiffness)/(total*omega**2)-1),
                nodes1=nodes[0]["counts"][1],nodes2=nodes[1]["counts"][1],
                nodes_stable=all(v["stable"] for v in nodes), node_details=json.dumps(nodes),
                axial_nodes1=axial_nodes[0]["counts"][1],axial_nodes2=axial_nodes[1]["counts"][1],
                axial_nodes_stable=all(v["stable"] for v in axial_nodes),axial_node_details=json.dumps(axial_nodes))


def continue_series(spec, data, ids, cfg):
    n=spec["branches"]; positions=np.arange(n); confirmed=np.ones(n,dtype=bool); signs=np.ones(n)
    rows, transitions, subspaces, allmac=[],[],[],[]
    for step,p in enumerate(data):
        if step:
            old=data[step-1]; f=old["frequencies"][positions]; scores=[]; assignments=[]
            for fraction in (0.,.5,1.):
                metric=metric_between(old,p,fraction)
                a,b=features(old,metric)[positions],features(p,metric)
                diag=assignment(a,b,f,p["frequencies"])
                assignments.append(diag.assignment);scores.append(diag.mac)
            mid=scores[1]; new=assignments[1]
            stable=np.all(np.stack(assignments)==new,axis=0)
            for penalty in (.03,.1):
                d=assignment(features(old,metric_between(old,p))[positions],features(p,metric_between(old,p)),f,p["frequencies"],penalty)
                stable &= d.assignment==new
            a,b=features(old,metric_between(old,p))[positions],features(p,metric_between(old,p))
            reverse=assignment(b[new],features(old,metric_between(old,p)),p["frequencies"][new],old["frequencies"])
            reverse_stable=reverse.assignment==positions
            stable &= reverse_stable
            # Adjacent seed pairs: principal angles, not a forced identification
            # of the individual vectors inside a nearly multiple subspace.
            for i in range(n-1):
                c=principal_cosines(a[[i,i+1]],b[new[[i,i+1]]])
                subspaces.append(dict(series=spec["id"],step=step,x_from=old["beta"] if spec["axis"]=="beta_deg" else old["mu"],x_to=p["beta"] if spec["axis"]=="beta_deg" else p["mu"],branch_a=ids[i][0],branch_b=ids[i+1][0],minimum_principal_cosine=float(c.min()),maximum_principal_cosine=float(c.max()),principal_angles_deg=json.dumps(np.rad2deg(np.arccos(np.clip(c,0,1))).tolist()),positions_from=json.dumps((positions[[i,i+1]]+1).tolist()),positions_to=json.dumps((new[[i,i+1]]+1).tolist())))
            for i,k in enumerate(new):
                ov=float(mid[i,k]);margin=ov-float(np.max(np.delete(mid[i],k)))
                jump=abs(int(k)-int(positions[i])); accepted=accept_transition(ov,margin,stable[i],jump,cfg["overlap_min"],cfg["margin_min"])
                status="ACCEPTED" if accepted else ("PREFIX_OR_ASSIGNMENT_UNRESOLVED" if k==len(p["frequencies"])-1 or positions[i]==len(old["frequencies"])-1 else "AMBIGUOUS_ASSIGNMENT")
                transitions.append(dict(series=spec["id"],step=step,branch_id=ids[i][0],x_from=old["beta"] if spec["axis"]=="beta_deg" else old["mu"],x_to=p["beta"] if spec["axis"]=="beta_deg" else p["mu"],previous_position=int(positions[i]+1),candidate_position=int(k+1),overlap=ov,MAC=ov**2,margin=margin,overlap_left=float(scores[0][i,k]),overlap_right=float(scores[2][i,k]),metric_penalty_stable=bool(stable[i]),reverse_stable=bool(reverse_stable[i]),position_jump=jump,local_status=status,previous_canonical_position=int(positions[i]+1) if confirmed[i] else None,origin_confirmed=bool(confirmed[i] and accepted)))
                confirmed[i] &= accepted
                signs[i] *= 1 if np.dot(a[i],b[k])>=0 else -1
                for col in range(len(p["frequencies"])):
                    allmac.append(dict(series=spec["id"],step=step,branch_id=ids[i][0],candidate_position=col+1,overlap_mid=float(mid[i,col]),MAC_mid=float(mid[i,col]**2),overlap_left=float(scores[0][i,col]),overlap_right=float(scores[2][i,col]),selected=col==k))
            positions=new.copy()
        for i,k in enumerate(positions):
            detail=characterize(p,k,cfg["node_relative_floors"])
            rows.append(dict(series=spec["id"],step=step,x=p["beta"] if spec["axis"]=="beta_deg" else p["mu"],beta_deg=p["beta"],mu=p["mu"],branch_id=ids[i][0],seed_position=i+1,seed_class=ids[i][1],seed_class_rank=ids[i][2],candidate_sorted_position=int(k+1),current_sorted_position=int(k+1) if confirmed[i] else None,Lambda=float(p["frequencies"][k]),status="CONFIRMED_FROM_SEED" if confirmed[i] else "UNRESOLVED_FROM_SEED",sign=float(signs[i]),shape_key=f"mode{k+1}",source_record=p["record"],source_npz=p["npz"],**detail))
    return rows,transitions,subspaces,allmac


def summary(spec,data,rows,transitions,cfg):
    gaps=[float(p["frequencies"][1]-p["frequencies"][0]) for p in data];j=int(np.argmin(gaps));p=data[j]
    sensitivity={}
    for ov in cfg["overlap_sensitivity"]:
        for margin in cfg["margin_sensitivity"]:
            rejected=[r for r in transitions if not accept_transition(r["overlap"],r["margin"],r["metric_penalty_stable"],r["position_jump"],ov,margin)]
            sensitivity[f"overlap={ov};margin={margin}"]=len(rejected)
    result=dict(nodes=len(data),saved_forms=sum(len(p["frequencies"]) for p in data),seed=rows[:spec["branches"]],minimum_sorted_gap=dict(step=j,beta_deg=p["beta"],mu=p["mu"],Lambda1=float(p["frequencies"][0]),Lambda2=float(p["frequencies"][1]),gap=gaps[j]),
                confirmed_rows=sum(r["status"]=="CONFIRMED_FROM_SEED" for r in rows),total_rows=len(rows),
                unresolved_steps=[r for r in transitions if r["local_status"]!="ACCEPTED"],
                overlap_quantiles=np.quantile([r["overlap"] for r in transitions],[0,.01,.5,1]).tolist(),
                threshold_sensitivity=sensitivity,maximum_mass_gram_error=max(p["mass_gram_error"] for p in data),
                maximum_energy_error=max(r["energy_error"] for r in rows))
    pair=[]
    for step in range(len(data)):
        rs=[r for r in rows if r["step"]==step and r["seed_position"]<=2]
        if len(rs)==2 and all(r["status"]=="CONFIRMED_FROM_SEED" for r in rs):
            pair.append((abs(rs[0]["Lambda"]-rs[1]["Lambda"]),step,rs))
    best=min(pair,key=lambda t:t[0]);result["minimum_confirmed_seed_pair_gap"]=dict(gap=best[0],step=best[1],x=best[2][0]["x"],branches=[r["branch_id"] for r in best[2]],frequencies=[r["Lambda"] for r in best[2]],complete_interval=len(pair)==len(data))
    if "pair" in spec:
        a,b=np.array(spec["pair"])-1
        gaps=[float(p["frequencies"][b]-p["frequencies"][a]) for p in data]
        j=int(np.argmin(gaps));p=data[j]
        result["minimum_requested_sorted_pair_gap"]=dict(sorted_positions=spec["pair"],step=j,
            beta_deg=p["beta"],gap=gaps[j],frequencies=p["frequencies"][[a,b]].tolist(),
            interpretation="minimum over saved nodes only")
    return result


def candidate_scenes(spec,data,rows,stats):
    # Select only existing nodes: angular minimum +/- nearest one-degree nodes;
    # length sequences use first confirmed changes in the sampled node count.
    grid=np.array([p["beta"] if spec["axis"]=="beta_deg" else p["mu"] for p in data])
    chosen=[]
    if spec["axis"]=="beta_deg":
        if "scene_angles" in spec:
            indices=[int(np.flatnonzero(np.isclose(grid,x,atol=1e-12,rtol=0))[0]) for x in spec["scene_angles"]]
            return sorted([r for r in rows if r["seed_position"] in spec["pair"] and r["step"] in indices],key=lambda r:(r["seed_position"],r["step"]))
        j=stats["minimum_sorted_gap"]["step"]
        indices=[int(np.argmin(abs(grid-(grid[j]-1)))),j,int(np.argmin(abs(grid-(grid[j]+1))))]
        for seed in (1,2):
            chosen += [r for r in rows if r["seed_position"]==seed and r["step"] in indices]
    else:
        for seed in (1,2,3):
            rr=[r for r in rows if r["seed_position"]==seed and r["status"]=="CONFIRMED_FROM_SEED"]
            changes=[r for a,r in zip(rr[:-1],rr[1:]) if r["nodes_stable"] and (a["nodes1"],a["nodes2"])!=(r["nodes1"],r["nodes2"]) and r["x"]>.03]
            selected=[rr[0]]
            for change in changes:
                target=min(change["step"]+3,rr[-1]["step"])
                r=next((v for v in rr if v["step"]==target),change)
                if r["step"]-selected[-1]["step"]>=8:
                    selected.append(r)
                if len(selected)==3:break
            if len(selected)<3:
                for target in (.4,.8):
                    r=min(rr,key=lambda r:abs(r["x"]-target))
                    if r not in selected:selected.append(r)
                    if len(selected)==3:break
            chosen += selected
    return sorted(chosen,key=lambda r:(r["seed_position"],r["step"]))


def character_comparisons(config,data_by_id,scenes):
    result=[]
    for spec in config["series"][:2]:
        pair=spec.get("pair",[1,2]);positions=np.array(pair)-1
        rows=[r for r in scenes if r["series"]==spec["id"] and r["seed_position"]==pair[0]]
        a,b=[data_by_id[spec["id"]][r["step"]] for r in (rows[0],rows[-1])]
        matrices=[]
        for t in (0.,.5,1.):
            metric=metric_between(a,b,t)
            fa,fb=features(a,metric)[positions],features(b,metric)[positions]
            matrices.append((normalized(fa)@normalized(fb).T)**2)
        c=principal_cosines(fa,fb)
        result.append(dict(series=spec["id"],beta_before=a["beta"],beta_after=b["beta"],
            MAC_squared_mid=matrices[1].tolist(),MAC_squared_left=matrices[0].tolist(),MAC_squared_right=matrices[2].tolist(),
            principal_cosines=c.tolist(),cross_similarity_dominates=bool(matrices[1][0,1]>matrices[1][0,0] and matrices[1][1,0]>matrices[1][1,1]),
            sorted_positions=pair,principal_angles_deg=np.rad2deg(np.arccos(np.clip(c,0,1))).tolist(),
            interpretation="finite-window comparison only; not an exact identity or a proof about the continuous minimum"))
    return result


def symmetric_reference(config, archive):
    """Reuse the approved six descendants, without reclassification/tracking."""
    path=archive.pin(ROOT/config["symmetric_reference"]); man=read(path);folder=path.parent
    expected=dict(aspect=1.,mu=0.,tau=0.,epsilon=.005,E=1.,rho=1.,l=1.,nu=.3,ks=5/6,
                  joint="RIGID",external="two clamped",time_reference=200.)
    if man["status"]!="COMPLETE" or any(man["physical_contract"].get(k)!=v for k,v in expected.items()):
        raise ValueError("Unaccepted symmetric reference")
    for name,h in man["artifacts"].items():archive.pin(folder/name,h)
    for name,h in man["identity"]["inputs"].items():archive.pin(ROOT/name,h)
    rows=[r for r in table(folder/"plotted_frequencies.csv") if int(r["base_sorted_index"])<=6]
    classification=table(MAPS/"_audit_reflection_classification_20261107/class_spectrum.csv")
    byclass={(float(r["beta_deg"]),r["symmetry_class"],int(r["class_rank"])):r for r in classification}
    local=table(MAPS/"_audit_reflection_intersection_20261107/block_roots.csv")
    bylocal={(float(r["beta_deg"]),r["symmetry_class"]):r for r in local}
    crossing=read(MAPS/"_audit_reflection_intersection_20261107/crossing_solution.json")
    source=archive.manifest(MAPS/"b775d0980d7bad77/manifest.json")
    records={float(r["x"]):r for r in source["points"] if r["map"]==man["identity"]["map_name"]}
    for r in rows:
        beta=float(r["beta_deg"]);lam=float(r["Lambda"]);cls=r["symmetry_class"]
        if r["plot_source"]=="class_spectrum.csv":
            cr=byclass[(beta,cls,int(r["class_rank"]))]
            if cr["classification_status"] not in ("PURE_PLUS","PURE_MINUS") or float(cr["Lambda"])!=lam:
                raise ValueError("Symmetric classification mismatch")
            item=records[beta];rp=archive.pin(ROOT/item["record"],item["record_sha256"]);rec=read(rp)
            rank=int(r["original_sorted_position"])
            if rec["roots"][rank-1]["Lambda"]!=lam:raise ValueError("Symmetric JSON/CSV mismatch")
            npz=archive.pin(rp.parent.parent/rec["shape_file"],rec["shape_sha256"])
            if npz!=ROOT/r["source_npz"]:raise ValueError("Symmetric NPZ mismatch")
        elif r["plot_source"]=="block_roots.csv":
            q=bylocal[(beta,cls)]
            if float(q["Lambda"])!=lam or q["first_class_certified"]!="True":raise ValueError("Local reference mismatch")
        elif r["plot_source"]=="crossing_solution.json":
            if beta!=crossing["beta_deg"] or lam!=crossing["Lambda"] or not crossing["first_class_roots_certified"]:
                raise ValueError("Crossing reference mismatch")
        else:raise ValueError("Unknown reference source")
    return rows


def pair_diagnostics(spec,data,rows):
    """Raw sorted pair observations remain usable when seed identity is lost."""
    pair=np.array(spec["pair"])-1; records=[];comparisons=[]
    selected=[p for p in data if spec["zoom"][0]<=p["beta"]<=spec["zoom"][1]]
    for step,p in enumerate(data):
        if not spec["zoom"][0]<=p["beta"]<=spec["zoom"][1]:continue
        for k in pair:
            canonical=next((r["branch_id"] for r in rows if r["step"]==step and r["current_sorted_position"]==k+1),"")
            records.append(dict(series=spec["id"],step=step,beta_deg=p["beta"],sorted_position=int(k+1),
                canonical_branch=canonical,Lambda=float(p["frequencies"][k]),source_record=p["record"],
                **characterize(p,int(k),[1e-6,1e-5,1e-4])))
    def compare(a,b,kind):
        metric=metric_between(a,b);fa,fb=features(a,metric)[pair],features(b,metric)[pair]
        overlap=normalized(fa)@normalized(fb).T;cs=principal_cosines(fa,fb)
        return dict(series=spec["id"],kind=kind,beta_from=a["beta"],beta_to=b["beta"],positions=spec["pair"],
                    overlap_absolute=abs(overlap).tolist(),MAC_squared=(overlap**2).tolist(),
                    principal_cosines=cs.tolist(),principal_angles_deg=np.rad2deg(np.arccos(np.clip(cs,0,1))).tolist())
    for a,b in zip(selected[:-1],selected[1:]):comparisons.append(compare(a,b,"adjacent"))
    scenes=[next(p for p in data if abs(p["beta"]-x)<1e-12) for x in spec["scene_angles"]]
    for i,j in ((0,1),(1,2),(0,2)):comparisons.append(compare(scenes[i],scenes[j],"selected_scenes"))
    return records,comparisons


def refinement_diagnostics(spec,data,rows,cfg):
    """Compare retained saved-node paths; no new spectra or invented midpoints."""
    lo,hi=spec["refinement_window"]
    window=[p for p in data if lo<=p["beta"]<=hi]
    local_rows=[r for r in rows if lo<=r["beta_deg"]<=hi]
    output=[];summaries=[]
    for stride in (10,5,2,1):
        indices=sorted(set(range(0,len(window),stride))|{len(window)-1})
        sample=[window[j] for j in indices]
        # Local IDs explicitly denote sorted positions at the window entrance.
        local_ids=[(f"local{k}","",0) for k in range(1,spec["branches"]+1)]
        rr,tr,_,_=continue_series(spec,sample,local_ids,cfg)
        key=f"retained_step_{.05*stride:g}"
        for t in tr:output.append(dict(path=key,**t))
        last=[r for r in rr if r["step"]==len(sample)-1]
        summaries.append(dict(path=key,beta_nodes=[p["beta"] for p in sample],
            minimum_overlap=min(t["overlap"] for t in tr),
            rejected_transitions=sum(t["local_status"]!="ACCEPTED" for t in tr),
            final_candidate_positions=[r["candidate_sorted_position"] for r in last],
            final_confirmed_positions=[r["current_sorted_position"] for r in last]))
    ids=[(f"B{k}","",0) for k in range(1,spec["branches"]+1)]
    reverse,rt,_,_=continue_series(spec,list(reversed(data)),ids,cfg)
    forward={(r["beta_deg"],r["branch_id"]):r for r in rows}
    matches=all(r["current_sorted_position"] is not None and
        r["current_sorted_position"]==forward[(r["beta_deg"],r["branch_id"])]["current_sorted_position"] for r in reverse)
    result=dict(series=spec["id"],window=[lo,hi],saved_node_paths=summaries,
        reverse_start="sorted positions at beta=90, independently propagated to beta=0",
        reverse_minimum_overlap=min(t["overlap"] for t in rt),
        reverse_rejected_transitions=sum(t["local_status"]!="ACCEPTED" for t in rt),
        forward_reverse_all_nodes_agree=matches,
        limitation="Step-halving compares saved subsets only; no step below 0.05 was computed")
    return result,output,reverse,rt,local_rows


def displayed_mode(rows,series,beta,seed):
    """Select by confirmed descendant, falling back only to an explicit raw position."""
    record=next((r for r in rows if r["series"]==series and r["beta_deg"]==beta and
                 r["seed_position"]==seed and r["status"]=="CONFIRMED_FROM_SEED"),None)
    return (int(record["current_sorted_position"])-1,record) if record else (seed-1,None)


def draw_publication(output,config,allrows,data_by_id,symmetric):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from datetime import datetime, timezone
    plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","DejaVu Serif"],
        "mathtext.fontset":"stix","font.size":10,"axes.linewidth":.65,"pdf.fonttype":42,"axes.grid":False})
    geometry={};captions={};scales=[];plotted_modes=[]
    for spec in config["series"]:
        data=data_by_id[spec["id"]];p=data[0]
        geometry[spec["id"]]=line_geometry.FrequencyGeometry([a["beta"] for a in data],axis="beta_deg",
            fixed=dict(mu=p["mu"],tau=p["tau"]),material=p["json"]["contract"]["model"])
    symgeometry={}
    for branch in (range(1,7) if symmetric else ()):
        rr=[r for r in symmetric if int(r["base_sorted_index"])==branch]
        symgeometry[branch]=line_geometry.FrequencyGeometry([float(r["beta_deg"]) for r in rr],axis="beta_deg",
            fixed=dict(mu=0,tau=0),material=dict(epsilon=.005,aspect=1))
    diagnostics=line_geometry.write_geometry_diagnostics(output,dict(geometry,**{f"symmetric_B{k}":v for k,v in symgeometry.items()}))
    def requested(name):
        return config.get("render_figures") is None or name in config["render_figures"]
    def save(fig,name):
        if not requested(name):
            plt.close(fig);return
        captions[name]=line_geometry.add_geometry_caption(fig)
        fig.savefig(output/(name+".pdf"),bbox_inches="tight",metadata={"CreationDate":datetime(2026,11,9,tzinfo=timezone.utc),"ModDate":None})
        fig.savefig(output/(name+".png"),bbox_inches="tight",dpi=config["dpi"]);plt.close(fig)
    def format_axis(ax):
        ax.set_xlabel(r"$\beta,\ ^\circ$");ax.set_ylabel(r"$\Lambda$");ax.tick_params(direction="in",width=.65)
    def panel(ax,spec=None,zoom=False):
        if spec is None:
            for branch in range(1,7):
                rr=[r for r in symmetric if int(r["base_sorted_index"])==branch]
                line_geometry.plot_frequency_curve(ax,symgeometry[branch],[float(r["Lambda"]) for r in rr],color=COLORS[branch-1],lw=1.)
        else:
            rr=[r for r in allrows if r["series"]==spec["id"]];geom=geometry[spec["id"]]
            branches=spec["pair"] if zoom else range(1,7)
            for branch in branches:
                br=[r for r in rr if r["seed_position"]==branch]
                yy=[r["Lambda"] if r["status"]=="CONFIRMED_FROM_SEED" else np.nan for r in br]
                line_geometry.plot_frequency_curve(ax,geom,yy,color=COLORS[branch-1],lw=1.)
            for step,p in enumerate(data_by_id[spec["id"]]):
                if zoom and not spec["zoom"][0]<=p["beta"]<=spec["zoom"][1]:continue
                confirmed={r["current_sorted_position"] for r in rr if r["step"]==step and r["current_sorted_position"] is not None}
                observed=spec["pair"] if zoom else range(1,7)
                missing=[k for k in observed if k not in confirmed]
                if missing:ax.plot([p["beta"]]*len(missing),p["frequencies"][np.array(missing)-1],ls="none",marker="o",mfc="white",mec=".35",mew=.5,ms=2.5 if zoom else 1.5)
            if zoom:
                ax.set_xlim(*spec["zoom"])
                values=[p["frequencies"][np.array(spec["pair"])-1] for p in data_by_id[spec["id"]] if spec["zoom"][0]<=p["beta"]<=spec["zoom"][1]]
                low,high=np.min(values),np.max(values);ax.set_ylim(low-.06*(high-low),high+.06*(high-low))
        if not zoom:ax.set_xlim(0,90);ax.set_ylim(*config.get("frequency_ylim",[2,11.5]));ax.set_xticks([0,30,60,90])
        format_axis(ax)
    layouts=[(3,"frequency_three_panels",[None,*config["series"]]),(2,"frequency_two_panels",config["series"])]
    if "frequency_panel_series" in config:
        selected=[next(s for s in config["series"] if s["id"]==name) for name in config["frequency_panel_series"]]
        layouts=[(len(selected),"frequency_three_panels",selected)]
    for count,name,specs in layouts:
        if not requested(name):continue
        fig,axes=plt.subplots(1,count,figsize=(3.45*count,3.15))
        for i,(ax,spec) in enumerate(zip(axes,specs)):
            panel(ax,spec);p=dict(mu=0,tau=0) if spec is None else data_by_id[spec["id"]][0]
            ax.text(.5,-.27,f"({'абв'[i]}) "+fr"$\mu={p['mu']:g},\ \tau={p['tau']:g}$",transform=ax.transAxes,ha="center")
        fig.subplots_adjust(wspace=.3,bottom=.24);save(fig,name)
    if requested("frequency_closeups"):
        fig,axes=plt.subplots(1,2,figsize=(8,3.25))
        for i,(ax,spec) in enumerate(zip(axes,config["series"])):
            panel(ax,spec,True);p=data_by_id[spec["id"]][0]
            ax.text(.5,-.27,f"({'аб'[i]}) "+fr"$\mu={p['mu']:g},\ \tau={p['tau']:g}$",transform=ax.transAxes,ha="center")
        fig.subplots_adjust(wspace=.28,bottom=.25);save(fig,"frequency_closeups")
    for spec in config["series"]:
        if not spec.get("scene_angles"):continue
        components=requested(spec["id"]+"_components")
        if not requested(spec["id"]+"_modes") and not components:continue
        data=data_by_id[spec["id"]];pair=np.array(spec["pair"])-1
        scenes=[next(p for p in data if abs(p["beta"]-x)<1e-12) for x in spec["scene_angles"]]
        selected={(seed,p["beta"]):displayed_mode(allrows,spec["id"],p["beta"],seed) for seed in spec["pair"] for p in scenes}
        peak=max(np.max(np.linalg.norm(a["y"][selected[(seed,p["beta"])][0],:,:2],axis=1)) for seed in spec["pair"] for p in scenes for a in p["arms"])
        scale=.22/float(peak);scales.append(dict(series=spec["id"],displacement_multiplier=scale,normalization="unit full inertia, one multiplier per sheet",labels="confirmed descendant; raw position if unresolved" if config.get("mode_label_style")=="branch_only" else "confirmed descendant and sorted position; raw position if unresolved",deformed_axis_color=".15",deformed_axis_linestyle="solid"))
        fig,axes=plt.subplots(2,3,figsize=(10,2.8))
        fig2,axes2=plt.subplots(2,3,figsize=(10,4.7),sharex=True,sharey=True) if components else (None,None)
        for row,seed in enumerate(spec["pair"]):
            for col,p in enumerate(scenes):
                k,record=selected[(seed,p["beta"])]
                first_w=p["arms"][0]["y"][k,:,1]
                sign=record["sign"] if record else (1. if first_w[np.argmax(abs(first_w))]>=0 else -1.)
                color=COLORS[record["seed_position"]-1] if record else ".35"
                plotted_modes.append(dict(series=spec["id"],beta_deg=p["beta"],sorted_position=int(k+1),
                    Lambda=float(p["frequencies"][k]),canonical_branch=record["branch_id"] if record else "",
                    sign=sign,phase_rule="confirmed_transport" if record else "positive_largest_arm1_W_display_only",
                    displacement_multiplier=scale,deformed_axis_color=".15",deformed_axis_linestyle="solid",source_record=p["record"],source_npz=p["npz"]))
                ax,ac=axes[row,col],axes2[row,col] if components else None
                for arm,sec,(t,n) in zip(p["arms"],p["sections"],frames(p["beta"])):
                    x=np.r_[0,arm["x"],sec["length"]];y=np.vstack((arm["ends"][k,0],arm["y"][k],arm["ends"][k,1]))*sign
                    base=(x-sec["length"])[:,None]*t;curve=base+scale*(y[:,0,None]*t+y[:,1,None]*n)
                    ax.plot(base[:,0],base[:,1],color=".75",lw=.55);ax.plot(curve[:,0],curve[:,1],color=".15",ls="-",lw=1.1)
                    coordinate=x/sec["length"]-1 if sec is p["sections"][0] else 1-x/sec["length"]
                    if ac is not None:
                        ac.plot(coordinate,y[:,0],color=".1",lw=.8,label="$U$" if sec is p["sections"][0] else None)
                        ac.plot(coordinate,y[:,1],color=color if record else ".6",lw=1.1,label="$W$" if sec is p["sections"][0] else None)
                ax.set_aspect("equal");ax.set_xlim(*config.get("mode_xlim",[-1.25,1.5]));ax.set_ylim(*config.get("mode_ylim",[-.30,.50]));ax.axis("off")
                if config.get("mode_label_style")=="ordinal_rows":
                    if record is None or int(record["current_sorted_position"])!=seed:
                        raise ValueError("Ordinal mode label requires a confirmed matching sorted position")
                    label=f"{seed}-я форма"
                    if col==0:ax.text(-.07,.5,label,transform=ax.transAxes,ha="right",va="center",fontsize=config.get("mode_label_fontsize",14))
                    if row==0:ax.set_title(fr"$\beta={p['beta']:g}^\circ$",fontsize=config.get("mode_label_fontsize",14),pad=10)
                elif config.get("mode_label_style")=="branch_only" and record:
                    branch=record["branch_id"][1:]
                    label=fr"$\mathrm{{B}}_{{{branch}}},\ \beta={p['beta']:g}^\circ$"
                else:
                    label=(record["branch_id"]+", " if record else "")+fr"$k={k+1},\ \beta={p['beta']:g}^\circ$"
                if config.get("mode_label_style")!="ordinal_rows":
                    ax.text(.5,-.06,label,transform=ax.transAxes,ha="center",fontsize=config.get("mode_label_fontsize"),va="top" if config.get("mode_label_style")=="branch_only" else "baseline")
                if ac is not None:
                    ac.set_title(label,fontsize=10);ac.tick_params(direction="in");ac.axvline(0,color=".8",lw=.5)
                    if row==1:ac.set_xlabel(r"$s$")
                    if col==0:ac.set_ylabel(r"$U,\ W$")
        fig.subplots_adjust(wspace=.04,hspace=.22)
        if config.get("mode_label_style")=="ordinal_rows":
            scales[-1]["labels"]="ordinal rows; confirmed matching positions"
            fig.subplots_adjust(left=.14,top=.82,bottom=.06,hspace=.35,wspace=.12)
        save(fig,spec["id"]+"_modes")
        if components:
            axes2[0,0].legend(frameon=False,fontsize=8);fig2.tight_layout();save(fig2,spec["id"]+"_components")
    write_json(output/"geometry_presentation.json",dict(geometry_diagnostics=diagnostics,geometry_captions=captions))
    write_csv(output/"plotted_mode_panels.csv",plotted_modes)
    return scales


def draw(output,config,allrows,candidates,data_by_id):
    if config.get("presentation")=="p26_figure5":
        symmetric=table(output/"symmetric_reference.csv") if (output/"symmetric_reference.csv").exists() else []
        return draw_publication(output,config,allrows,data_by_id,symmetric)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family":"serif","font.serif":["Times New Roman","DejaVu Serif"],"mathtext.fontset":"stix","font.size":11,"axes.linewidth":.7,"pdf.fonttype":42,"axes.grid":False})
    geometries, captions = {}, {}
    for spec in config["series"]:
        points = data_by_id[spec["id"]]
        model = points[0]["json"]["contract"]["model"]
        geometries[spec["id"]] = line_geometry.FrequencyGeometry(
            [p["beta"] if spec["axis"] == "beta_deg" else p["mu"] for p in points],
            axis=spec["axis"], fixed=dict(beta_deg=points[0]["beta"]), material=model)
    geometry_diagnostics = line_geometry.write_geometry_diagnostics(output, geometries)
    def save(fig,name):
        captions[name] = line_geometry.add_geometry_caption(fig)
        fig.savefig(output/(name+".pdf"),bbox_inches="tight")
        fig.savefig(output/(name+".png"),bbox_inches="tight",dpi=config["dpi"])
        plt.close(fig)
    def panel(ax,spec):
        rr=[r for r in allrows if r["series"]==spec["id"]]
        for branch in sorted({r["branch_id"] for r in rr},key=lambda x:int(x[1:])):
            values=[r for r in rr if r["branch_id"]==branch]
            ys=[r["Lambda"] if r["status"]=="CONFIRMED_FROM_SEED" else np.nan for r in values]
            geometry = geometries[spec["id"]]
            if tuple(float(r["x"]) for r in values) != geometry.x:
                raise ValueError("Branch observations differ from the saved geometry grid")
            line_geometry.plot_frequency_curve(ax,geometry,ys,lw=1.2,color=COLORS[(int(branch[1:])-1)%8],label=branch)
        ax.set_xlabel(r"$\beta,\ ^\circ$" if spec["axis"]=="beta_deg" else r"$\mu$")
        ax.set_ylabel(r"$\Lambda$");ax.set_xlim(min(geometries[spec["id"]].x),max(geometries[spec["id"]].x))
        ax.tick_params(direction="in",width=.7);ax.legend(frameon=False,ncol=4,fontsize=8,loc="upper center",bbox_to_anchor=(.5,1.14))
    for spec in config["series"]:
        fig,ax=plt.subplots(figsize=(6.6,4.2));panel(ax,spec);save(fig,spec["id"]+"_branches")
    for name,specs in (("beta_comparison",config["series"][:2]),("mu_comparison",config["series"][2:])):
        fig,axes=plt.subplots(1,2,figsize=(11.4,3.8))
        for i,(ax,spec) in enumerate(zip(axes,specs)):
            panel(ax,spec);ax.text(.5,-.24,"(а)" if i==0 else "(б)",transform=ax.transAxes,ha="center")
        fig.subplots_adjust(wspace=.22);save(fig,name)
    scales=[]
    for spec in config["series"]:
        scenes=[r for r in candidates if r["series"]==spec["id"]]
        if not scenes:
            continue
        nrows=2 if spec["axis"]=="beta_deg" else 3
        fig,axes=plt.subplots(nrows,3,figsize=(10,2*nrows),squeeze=False)
        # Normalize with the existing full-inertia norm; ONE multiplier for all panels.
        peak=max(np.max(np.linalg.norm(a["y"][int(r["candidate_sorted_position"])-1,:,:2],axis=1))/np.sqrt(r["mass"])
                 for r in scenes for a in data_by_id[spec["id"]][r["step"]]["arms"])
        scale=.22/float(peak)
        scales.append(dict(series=spec["id"],displacement_multiplier=scale,normalization="unit full inertia; shared multiplier for entire sheet"))
        for ax in axes.flat:ax.set_visible(False)
        for ax,r in zip(axes.flat,scenes):
            ax.set_visible(True);p=data_by_id[spec["id"]][r["step"]];k=int(r["candidate_sorted_position"])-1
            for arm,sec,(t,n),color in zip(p["arms"],p["sections"],frames(p["beta"]),COLORS[:2]):
                x=np.r_[0,arm["x"],sec["length"]];y=np.vstack((arm["ends"][k,0],arm["y"][k],arm["ends"][k,1]))
                base=(x-sec["length"])[:,None]*t
                displacement=(y[:,0,None]*t+y[:,1,None]*n)*r["sign"]/np.sqrt(r["mass"])
                ax.plot(base[:,0],base[:,1],color=".65",lw=.6,ls="--")
                curve=base+scale*displacement;ax.plot(curve[:,0],curve[:,1],color=color,lw=1.1)
            ax.set_aspect("equal");ax.set_xlim(-1.08,1.88);ax.set_ylim(-.28,1.13)
            ax.axis("off")
            label=(r"$\beta="+f"{r['beta_deg']:g}"+r"^\circ$" if spec["axis"]=="beta_deg" else r"$\mu="+f"{r['mu']:g}"+"$")
            ax.text(.5,-.04,r["branch_id"]+", "+label,transform=ax.transAxes,ha="center",fontsize=10)
        fig.subplots_adjust(hspace=.05,wspace=.04);save(fig,spec["id"]+"_candidate_modes")
    write_json(output/"geometry_presentation.json",dict(geometry_diagnostics=geometry_diagnostics,geometry_captions=captions))
    return scales


def verified_redraw_source(source):
    """Verify a historical derived bundle; only the current renderer may differ."""
    source=Path(source);old=read(source/"manifest.json")
    for path,h in old["inputs"].items():
        if (ROOT/path).resolve()==Path(__file__).resolve():continue
        if sha(ROOT/path)!=h:raise ValueError("Changed redraw input: "+path)
    for name,h in old["outputs"].items():
        if sha(source/name)!=h:raise ValueError("Changed derived data: "+name)
    return old


def confirmed_redraw_rows(rows,spec,points):
    """Select saved descendants and reconcile their provenance with physical fields."""
    int_fields={"step","seed_position","seed_class_rank","candidate_sorted_position","current_sorted_position"}
    float_fields={"x","beta_deg","mu","Lambda","sign","mass"}
    selected=[]
    for r in rows:
        if r["series"]!=spec["id"] or int(r["seed_position"])>spec["branches"]:continue
        q={k:int(v) if k in int_fields and v else float(v) if k in float_fields and v else v for k,v in r.items()}
        if q["status"]!="CONFIRMED_FROM_SEED":raise ValueError("Unconfirmed saved descendant in requested series")
        selected.append(q)
    for step,p in enumerate(points):
        rr=[r for r in selected if r["step"]==step]
        if sorted(r["seed_position"] for r in rr)!=list(range(1,spec["branches"]+1)):
            raise ValueError("Missing or duplicate saved descendants")
        if sorted(r["current_sorted_position"] for r in rr)!=list(range(1,spec["branches"]+1)):
            raise ValueError("Selected descendants do not cover the requested lowest frequencies")
        for r in rr:
            k=r["current_sorted_position"]-1
            if (r["beta_deg"]!=p["beta"] or r["mu"]!=p["mu"] or r["Lambda"]!=p["frequencies"][k]
                    or r["source_record"]!=p["record"] or r["source_npz"]!=p["npz"]):
                raise ValueError("Saved continuation and physical archive disagree")
            if r["sign"] not in (-1.,1.):raise ValueError("Invalid transported phase")
    if len(selected)!=len(points)*spec["branches"]:raise ValueError("Unexpected saved grid points")
    return sorted(selected,key=lambda r:(r["step"],r["seed_position"]))


def redraw_selected(output,config):
    """Compose specified confirmed series from existing bundles, without tracking."""
    output=Path(output)
    if output.exists():raise FileExistsError("Refusing to overwrite results: "+str(output))
    sources={name:ROOT/path for name,path in config["continuation_sources"].items()}
    manifests={path:verified_redraw_source(path) for path in set(sources.values())}
    archive=Archive();data={};rows=[];provenance={};inputs={}
    for path,manifest in manifests.items():
        inputs.update({p:h for p,h in manifest["inputs"].items() if (ROOT/p).resolve()!=Path(__file__).resolve()})
        inputs[(path/"manifest.json").relative_to(ROOT).as_posix()]=sha(path/"manifest.json")
    with saved_only():
        for spec in config["series"]:
            name=spec["id"];source=sources[name];data[name]=archive.series(spec)
            selected=confirmed_redraw_rows(table(source/"branch_records.csv"),spec,data[name])
            for r in selected:r["continuation_source"]=source.relative_to(ROOT).as_posix()
            rows.extend(selected)
            provenance[name]=dict(source=source.relative_to(ROOT).as_posix(),manifest_sha256=sha(source/"manifest.json"),
                records_sha256=sha(source/"branch_records.csv"),selected_records=len(selected),K_plot=spec["branches"],K_guard=spec.get("K_guard",9))
    inputs.update(archive.hashes);inputs[Path(__file__).relative_to(ROOT).as_posix()]=sha(__file__)
    output.mkdir(parents=True)
    with saved_only():scales=draw(output,config,rows,[],data)
    fields=list(dict.fromkeys(k for row in rows for k in row))
    write_csv(output/"branch_records.csv",[{k:r.get(k,"") for k in fields} for r in rows])
    write_json(output/"presentation.json",scales)
    write_json(output/"source_series.json",provenance)
    for path,h in inputs.items():
        if sha(ROOT/path)!=h:raise ValueError("Input changed during rendering: "+path)
    write_json(output/"manifest.json",dict(schema="saved-square-presentation-v1",status="COMPLETE",configuration=config,
        inputs=inputs,source_series=provenance,source_unique_points=len(archive.cache),
        source_unique_forms=sum(len(p["frequencies"]) for p in archive.cache.values()),
        assignment_repeated=False,spectral_calls=0,basis_calls=0,ritz_calls=0,fem_calls=0,
        **read(output/"geometry_presentation.json"),outputs={p.name:sha(p) for p in output.iterdir() if p.is_file()}))
    return output


def redraw(source,output,config):
    """Render derived tables and original NPZ without repeating assignment."""
    if config.get("continuation_sources"):
        return redraw_selected(output,config)
    source,output=Path(source),Path(output)
    old=read(source/"manifest.json")
    for path,h in old["inputs"].items():
        # A redraw uses the current presenter, preserving the old presenter hash
        # in redraw_source. All scientific inputs must still match verbatim.
        if (ROOT/path).resolve() == Path(__file__).resolve():continue
        if sha(ROOT/path)!=h:raise ValueError("Changed redraw input: "+path)
    for name,h in old["outputs"].items():
        if sha(source/name)!=h:raise ValueError("Changed derived data: "+name)
    output.mkdir(parents=True,exist_ok=False)
    archive=Archive();data={}
    with saved_only():
        for spec in config["series"]:data[spec["id"]]=archive.series(spec)
    numeric_int={"step","seed_position","seed_class_rank","candidate_sorted_position","current_sorted_position","nodes1","nodes2"}
    numeric_float={"x","beta_deg","mu","Lambda","sign","inertia_arm1","inertia_arm2","inertia_axial","inertia_rotation","mass","energy_error"}
    def convert(rows):
        return [{k:(int(v) if k in numeric_int and v else float(v) if k in numeric_float else v=="True" if k=="nodes_stable" else v) for k,v in r.items()} for r in rows]
    rows=convert(table(source/"branch_records.csv"));scenes=convert(table(source/"candidate_modes.csv"))
    if config.get("presentation")=="p26_figure5":shutil.copyfile(source/"symmetric_reference.csv",output/"symmetric_reference.csv")
    scales=draw(output,config,rows,scenes,data)
    write_json(output/"presentation.json",scales)
    for f in source.iterdir():
        if f.suffix in (".csv",".json") and f.name not in ("manifest.json","geometry_nodes.csv","geometry_diagnostics.json","geometry_presentation.json","presentation.json","plotted_mode_panels.csv"):shutil.copyfile(f,output/f.name)
    write_json(output/"manifest.json",dict(old,**read(output/"geometry_presentation.json"),configuration=config,renderer_sha256=sha(__file__),redraw_source=str(source),redraw_source_sha256=sha(source/"manifest.json"),outputs={p.name:sha(p) for p in output.iterdir()},assignment_repeated=False))
    return output


def run(config_path,output_override=None):
    config=read(config_path);output=Path(output_override) if output_override else ROOT/config["output"]
    if output.exists() and any(output.iterdir()):
        raise FileExistsError("Refusing to overwrite results: "+str(output))
    output.mkdir(parents=True,exist_ok=True)
    archive=Archive();archive.pin(config_path);archive.pin(Path(__file__))
    archive.pin(ROOT/"scripts/lib/analytic_branch_tracking.py")
    for source in ("scripts/lib/isotropic_two_springs_rl_timoshenko.py", "scripts/lib/isotropic_two_springs_thickness_asymmetry.py", "scripts/lib/rl_full_rr_aspect.py", "scripts/lib/rl_full_rr_sections_prefix.py"):
        archive.pin(ROOT/source)
    allrows,transitions,subspaces,macs,scenes,stats,data_by_id=[],[],[],[],[],{},{}
    with saved_only():
        if config.get("presentation")=="p26_figure5":
            archive.pin(ROOT/config["publication_pdf"])
            write_csv(output/"symmetric_reference.csv",symmetric_reference(config,archive))
        for spec in config["series"]:
            data=archive.series(spec);data_by_id[spec["id"]]=data
            ids=seed_ids(spec,data,archive)
            rows,tr,sub,mac=continue_series(spec,data,ids,config)
            stats[spec["id"]]=summary(spec,data,rows,tr,config)
            scenes+=candidate_scenes(spec,data,rows,stats[spec["id"]])
            allrows+=rows;transitions+=tr;subspaces+=sub;macs+=mac
            print(spec["id"],len(data),"nodes;",stats[spec["id"]]["confirmed_rows"],"confirmed branch records;",len(stats[spec["id"]]["unresolved_steps"]),"unresolved transitions",flush=True)
    for name,rows in [("branch_records",allrows),("transitions",transitions),("pair_subspaces",subspaces),("MAC_candidates",macs),("candidate_modes",scenes)]:
        write_csv(output/(name+".csv"),rows)
    write_json(output/"quality_checks.json",stats)
    write_json(output/"character_comparisons.json",character_comparisons(config,data_by_id,scenes))
    if config.get("presentation")=="p26_figure5":
        pairrows,pairchecks=[],[]
        with saved_only():
            for spec in config["series"]:
                rr=[r for r in allrows if r["series"]==spec["id"]]
                pr,pc=pair_diagnostics(spec,data_by_id[spec["id"]],rr)
                pairrows+=pr;pairchecks+=pc
        write_csv(output/"pair_observations.csv",pairrows)
        write_json(output/"pair_comparisons.json",pairchecks)
    for spec in config["series"]:
        if "refinement_window" not in spec:continue
        with saved_only():
            result,steps,rr,rt,local=refinement_diagnostics(spec,data_by_id[spec["id"]],
                [r for r in allrows if r["series"]==spec["id"]],config)
        write_json(output/"step_refinement.json",result)
        write_csv(output/"step_refinement_transitions.csv",steps)
        write_csv(output/"reverse_branch_records.csv",rr)
        write_csv(output/"reverse_transitions.csv",rt)
        write_csv(output/"refined_window_branches.csv",local)
    scales=draw(output,config,allrows,scenes,data_by_id)
    write_json(output/"presentation.json",scales)
    for p,h in archive.hashes.items():
        if sha(ROOT/p)!=h:raise ValueError("Input modified during postprocessing")
    write_json(output/"manifest.json",dict(schema="saved-square-continuation-v1",status="COMPLETE_WITH_EXPLICIT_TRACKING_LIMITS" if any(v["unresolved_steps"] for v in stats.values()) else "COMPLETE",inputs=archive.hashes,outputs={p.name:sha(p) for p in output.iterdir() if p.is_file()},spectral_calls=0,basis_calls=0,ritz_calls=0,fem_calls=0,source_unique_points=len(archive.cache),source_unique_forms=sum(len(p["frequencies"]) for p in archive.cache.values()),configuration=config,coordinate_metric="global translation; common xi; midpoint lengths and physical inertia; derivative pulled back to xi",canonical_after_ambiguity=False,**read(output/"geometry_presentation.json")))
    return output


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=DEFAULT)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--redraw-from",type=Path,help="Reuse derived branch tables; requires a new --output directory")
    args=parser.parse_args()
    if args.redraw_from:
        if not args.output:parser.error("--redraw-from requires --output")
        print(redraw(args.redraw_from,args.output,read(args.config)))
    else:
        print(run(args.config,args.output))


if __name__=="__main__":
    main()
