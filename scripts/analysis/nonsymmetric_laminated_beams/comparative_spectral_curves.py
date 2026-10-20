"""NSL-VIS-1: two fixed Ritz maps, benchmark first, no damped roots.

Script Proliferation Control: this finite-dimensional map/cache contract is
different from the immutable prefix and complex verification runners. Reuse
the NSL-2 Ritz solver and basis; do not create another eigensolver.
The plot-only branch imports no physical provider and reads completed CSVs.
"""
from __future__ import annotations

import argparse
import csv
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
CONFIG = ROOT/"data/input/nonsymmetric_laminated_beam_comparative_spectral_curves.json"
OUTPUT = ROOT/"results/nonsymmetric_laminated_beams/comparative_spectral_curves"
CASES = ("BASE", "MIRROR", "SAME_LOCAL")


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plain(x):
    import numpy as np
    if isinstance(x, dict):
        return {k: plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def write(path, value):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(plain(value), indent=2, ensure_ascii=False)+"\n", encoding="utf-8")
    temporary.replace(path)


def table(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def grids(config):
    import numpy as np
    spring = config["parameter_grid"]["spring"]
    values = np.expm1(np.linspace(0, np.log1p(spring["maximum"]), spring["uniform_log1p_count"]))
    values[0], values[-1] = 0., float(spring["maximum"])
    values = sorted(set(map(float, values)) | set(map(float, spring["include_exact"])))
    return list(range(91)), values


def key(case, beta, kappa):
    # Full binary float identity, not rounded display labels.
    token = json.dumps([case, float(beta), kappa], separators=(",", ":"))
    return hashlib.sha256(token.encode()).hexdigest()[:16]


def context():
    from scripts.analysis.nonsymmetric_laminated_beams import confirm_angle_complex as previous
    config = read(CONFIG)
    _, old, source, props, ref, manifests, kv, screen, _ = previous.load_context()
    oldconfig, *_ = previous.load_context()
    fp, ident = previous.identity(oldconfig, manifests)
    assert fp == config["source_bundle"]
    previous.validate_cache(previous.OUTPUT/fp, ident)
    manifests[fp] = dict(path=(previous.OUTPUT/fp/"manifest.json").relative_to(ROOT).as_posix(),
                         sha256=sha(previous.OUTPUT/fp/"manifest.json"))
    baseline_path = OUTPUT/"baseline.json"
    if not baseline_path.exists():
        hashes = dict(previous.preservation({k: v for k, v in manifests.items() if k != fp})["sha256"])
        for item in manifests.values():
            for p in (ROOT/item["path"]).parent.rglob("*"):
                if p.is_file():
                    hashes[p.relative_to(ROOT).as_posix()] = sha(p)
        for p in [previous.CONFIG, Path(previous.__file__), Path(previous.figures.__file__),
                  ROOT/"docs/nonsymmetric_laminated_beams/nsl3a_angle_complex_confirmation.md",
                  previous.OUTPUT/"README.md", previous.OUTPUT/"verification.json",
                  ROOT/"tests/test_nsl3a_angle_complex_confirmation.py"]:
            hashes[p.relative_to(ROOT).as_posix()] = sha(p)
        memory = {}
        for name in ("decisions.md", "knowledge.md"):
            p = ROOT/"docs/memory"/name
            memory[p.relative_to(ROOT).as_posix()] = dict(size=p.stat().st_size, sha256=sha(p))
        write(baseline_path, dict(sha256=hashes, memory_prefix=memory))
    baseline = read(baseline_path)
    assert all(sha(ROOT/n) == h for n, h in baseline["sha256"].items())
    for name, value in baseline["memory_prefix"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()[:value["size"]]).hexdigest() == value["sha256"]
    assert config["variant"] == old["variant"] == "PROJECT_REDDY_COUPLED"
    assert config["c_theta"] == 0 and config["ritz_orders"] == [36, 44]
    import numpy as np
    import scipy
    from scripts.lib import nonsymmetric_laminated_joint as joint
    paths = [CONFIG, Path(__file__), Path(joint.__file__), Path(joint.nsl.__file__),
             ROOT/"scripts/lib/reddy_symmetric_coupled_beams_ritz.py"]
    identity = dict(config=config, source_manifests=manifests, baseline_sha256=sha(baseline_path),
                    source_hashes={p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
                    versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:16]
    return config, old, source, props, ref, screen, identity, fingerprint


@lru_cache(maxsize=16)
def basis(order, quadrature):
    import numpy as np
    from scripts.lib.nonsymmetric_laminated_joint import shifted_legendre_clamped_basis
    x, w = np.polynomial.legendre.leggauss(quadrature)
    b = shifted_legendre_clamped_basis(order, (x+1)/2)
    return (x+1)/2, w/2, b.values, b.derivatives


def physical(coefficients, properties, lengths, law, quadrature):
    """Vectorized evaluation of the unchanged NSL-2 physical Rayleigh forms."""
    import numpy as np
    order = coefficients.shape[0]//6
    xi, weights, v, d = basis(order, quadrature)
    cc = coefficients.reshape(2, 3, order, -1)
    q = np.einsum("nx,afnk->axfk", v, cc)
    dq = np.einsum("nx,afnk->axfk", d, cc)/np.array(lengths)[:, None, None, None]
    gram = np.zeros((coefficients.shape[1],)*2)
    components = {n: np.zeros(coefficients.shape[1]) for n in ("axial", "mixed", "bending", "shear", "spring")}
    for i, (a, length) in enumerate(zip(properties, lengths)):
        ww = weights*length
        gram += np.einsum("x,xik,ij,xjl->kl", ww, q[i], a.mass_matrix(), q[i])
        components["axial"] += ww@(a.A*dq[i, :, 0]**2)
        components["mixed"] += ww@(2*a.B*dq[i, :, 0]*dq[i, :, 2])
        components["bending"] += ww@(a.D*dq[i, :, 2]**2)
        components["shear"] += ww@(a.S*(dq[i, :, 1]+q[i, :, 2])**2)
    delta = cc[0, 2].sum(axis=0)-cc[1, 2].sum(axis=0)
    if law.mode != "RIGID":
        components["spring"] = law.k_theta*delta**2
    mass = np.diag(gram)
    energy = sum(components.values())
    return dict(q=q, dq=dq, gram=gram, mass=mass, components=components,
                energy=energy, omega=np.sqrt(energy/mass), delta=delta, xi=xi, weights=weights)


def overlap(q, v, properties, lengths, weights):
    import numpy as np
    return sum(np.einsum("x,xik,ij,xjl->kl", weights*length, a, p.mass_matrix(), b)
               for a, b, p, length in zip(q, v, properties, lengths))


def clusters(frequencies, tolerance):
    groups = [[0]]
    for i in range(1, len(frequencies)):
        if (frequencies[i]-frequencies[i-1])/frequencies[i] <= tolerance:
            groups[-1].append(i)
        else:
            groups.append([i])
    return groups


def structural(properties, lengths, beta, law, order, case):
    """No extra eigensolve: old assembly, reflection and allowed trace space."""
    import numpy as np
    from scripts.lib import nonsymmetric_laminated_joint as joint
    k, m, z, c, jump = joint.ritz_matrices(properties, lengths, beta, law, order)
    signs = np.repeat([1., -1., -1.], order)
    def reflected(a):
        b = a.reshape(2, 3*order, -1)[::-1]*signs[None, :, None]
        return b.reshape(a.shape)
    sr = z.T@reflected(z)
    kr, mr = z.T@k@z, z.T@m@z
    rel = lambda a, b: float(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.linalg.norm(b), 1e-300))
    ks, ms = rel(sr.T@kr@sr, kr), rel(sr.T@mr@sr, mr)
    subspace = float(np.linalg.norm(c@reflected(z))/np.linalg.norm(c))
    jump_minus = float(np.linalg.norm(jump@(z-reflected(z))/2))
    symmetry = case in ("BASE", "MIRROR")
    return dict(constraint_rows=len(c), expected_rows=3 if law.mode == "RIGID" else 2,
                reflection_K=ks, reflection_M=ms, reflection_constraint=subspace,
                protected_jump=jump_minus, reflection_claim=symmetry,
                pass_structure=(len(c) == (3 if law.mode == "RIGID" else 2) and
                                ((max(ks, ms, subspace, jump_minus) < 1e-11) if symmetry else ks > 1e-6)))


def solve_order(properties, lengths, beta, law, order, case, policy):
    import numpy as np
    from scripts.lib import nonsymmetric_laminated_joint as joint
    start = time.perf_counter()
    raw, coeff, diagnostics = joint.ritz_spectrum(properties, lengths, beta, law, order, 64)
    d = physical(coeff, properties, lengths, law, 64)
    coeff = coeff/np.sqrt(d["mass"])[None, :]
    d = physical(coeff, properties, lengths, law, 64)
    check = physical(coeff, properties, lengths, law, 96)
    permutation = np.argsort(d["omega"], kind="stable")
    coeff, raw = coeff[:, permutation], raw[permutation]
    if not np.array_equal(permutation, np.arange(7)):
        d = physical(coeff, properties, lengths, law, 64)
        check = physical(coeff, properties, lengths, law, 96)
    groups = clusters(d["omega"], policy["cluster_relative_gap"])
    cluster_of = {j: g for g in groups for j in g}
    parity = [None]*7
    if case in ("BASE", "MIRROR"):
        reflected = d["q"][::-1]*np.array([1, -1, -1])[None, None, :, None]
        pair = np.diag(overlap(d["q"], reflected, properties, lengths, d["weights"]))
        for j in range(7):
            eta = 1 if pair[j] >= 0 else -1
            error = reflected[..., j:j+1]-eta*d["q"][..., j:j+1]
            norm = float(np.sqrt(max(0, overlap(error, error, properties, lengths, d["weights"])[0, 0])))
            parity[j] = dict(eta=eta if norm <= policy["reflection_mass_tol"] else None,
                             residual=norm, pairing=float(pair[j]),
                             status="RESOLVED" if norm <= policy["reflection_mass_tol"] else "CLUSTER_UNRESOLVED")
    metrics = dict(
        mass_gram=float(np.max(abs(d["gram"]-np.eye(7)))),
        energy_quadrature=float(np.max(abs(check["energy"]/d["energy"]-1))),
        mass_quadrature=float(np.max(abs(check["mass"]/d["mass"]-1))),
        raw_rayleigh_frequency=float(np.max(abs(raw/d["omega"]-1))),
        algebraic_residual=float(max(diagnostics["eigen_residual"])),
        kinematic=float(diagnostics["constraint_residual"]/max(1., np.max(abs(coeff)))),
        energy_raw_eigenvalue=float(np.max(abs(d["energy"]-raw**2*d["mass"])/d["energy"])))
    gates = dict(
        positive=bool(np.all(d["omega"] > 0)),
        mass=metrics["mass_gram"] <= policy["mass_orthogonality_tol"],
        energy=metrics["energy_quadrature"] <= policy["quadrature_energy_rtol"],
        mass_quadrature=metrics["mass_quadrature"] <= policy["quadrature_energy_rtol"],
        raw_rayleigh=metrics["raw_rayleigh_frequency"] <= policy["raw_rayleigh_frequency_rtol"],
        algebraic=metrics["algebraic_residual"] <= policy["algebraic_residual_tol"],
        kinematic=metrics["kinematic"] <= policy["kinematic_tol"],
        constraints=diagnostics["constraint_rank"] == (3 if law.mode == "RIGID" else 2),
        parity=all(p is None or p["eta"] is not None or len(cluster_of[j]) > 1 for j, p in enumerate(parity)))
    record = dict(order=order, omega=d["omega"], raw_lapack_omega=raw,
                  coefficients=coeff, diagnostics=diagnostics, metrics=metrics,
                  gates=gates, parity=parity, clusters=[[j+1 for j in g] for g in groups],
                  mass=d["mass"], delta_psi_raw=d["delta"], stiffness=d["energy"],
                  components=d["components"], seconds=time.perf_counter()-start)
    return record, d


def compare_spaces(first, second, properties, lengths, policy):
    import numpy as np
    cross = overlap(first["q"], second["q"], properties, lengths, first["weights"])
    groups = clusters(second["omega"], policy["cluster_relative_gap"])
    rows = []
    for g in groups:
        small = cross[np.ix_(g, g)]
        sv = np.linalg.svd(small, compute_uv=False)
        rows.append(dict(positions=[i+1 for i in g], dimension=len(g),
                         individual_identity=len(g) == 1,
                         minimum_MAC=float(min(sv)**2),
                         one_minus_MAC=float(max(0, 1-min(sv)**2))))
    return dict(relative_frequency=abs(first["omega"]/second["omega"]-1),
                mass_overlap=cross, groups=rows,
                maximum_one_minus_MAC=max(r["one_minus_MAC"] for r in rows))


def material_contract(props):
    import numpy as np
    from dataclasses import asdict
    base = props["BASE"]
    differences = {n: max(abs(getattr(a, n)/getattr(base, n)-1) for a in props.values())
                   for n in ("A", "D", "S", "m", "j2")}
    assert max(differences.values()) < 1e-13
    assert base.B == 0 and np.isclose(props["HHLL"].B, -props["LLHH"].B, rtol=1e-13, atol=0)
    assert props["HHLL"].B != 0 and all(a.j1 == 0 for a in props.values())
    return dict(status="PASS",properties={n: asdict(a) for n,a in props.items()},
                relative_equalities=differences,BASE_is_homogeneous=True)


def get_point(case, beta, kappa, env, out):
    import numpy as np
    from scripts.lib import nonsymmetric_laminated_joint as joint
    config, old, props, ref = env
    point_id = key(case, beta, kappa)
    path = out/"points"/(point_id+".json")
    if path.exists():
        point = read(path)
        assert point["identity"] == [case, beta, kappa]
        assert sha(out/"points"/(point_id+".npz")) == point["coefficients_sha256"]
        return point, False
    properties = [props[n] for n in old["cases"][case]]
    law = joint.Joint("RIGID") if kappa is None else joint.Joint("SPRING", kappa*ref["D"]/ref["length"])
    records=[]; fields=[]
    ledger_path=out/"ledger.json"
    for order in config["ritz_orders"]:
        ledger=read(ledger_path)
        ledger["ritz_eigenproblems"]+=1
        write(ledger_path,ledger)
        record, d = solve_order(properties, old["lengths"], np.deg2rad(beta), law, order, case, config["policy"])
        records.append(record); fields.append(d)
    convergence=compare_spaces(*fields,properties,old["lengths"],config["policy"])
    structure=structural(properties,old["lengths"],np.deg2rad(beta),law,44,case)
    gates=dict(orders=all(all(r["gates"].values()) for r in records),
               convergence=max(convergence["relative_frequency"]) <= config["policy"]["frequency_convergence_rtol"],
               forms=convergence["maximum_one_minus_MAC"] <= config["policy"]["one_minus_mac_tol"],
               structure=structure["pass_structure"])
    path.parent.mkdir(parents=True,exist_ok=True)
    coeff_path=path.with_suffix(".npz")
    with coeff_path.open("wb") as stream:
        np.savez_compressed(stream,**{f'N{r["order"]}':r.pop("coefficients") for r in records})
    point=dict(identity=[case,beta,kappa],status="PASS" if all(gates.values()) else "NUMERICAL_UNRESOLVED",
               orders=records,convergence=convergence,structure=structure,gates=gates,
               coefficients_sha256=sha(coeff_path),solver="unchanged NSL-2 ritz_spectrum + physical Rayleigh forms")
    write(path,point)
    return plain(point),True


def reference_check(point, saved, env, out):
    """Evaluate saved analytic fields; no nullspace recovery and no roots."""
    import numpy as np
    from scripts.lib import nonsymmetric_laminated_joint as joint
    config, old, props, ref=env
    case,beta,kappa=point["identity"]
    properties=[props[n] for n in old["cases"][case]]
    path=out/"points"/(key(case,beta,kappa)+".npz")
    coeff=np.load(path)["N44"]
    law=joint.Joint("SPRING",ref["k_theta"])
    d=physical(coeff,properties,old["lengths"],law,96)
    model=joint.CoupledBeams(properties,old["lengths"],np.deg2rad(beta),law)
    q=np.zeros_like(d["q"])
    frequencies=[]
    for j, mode in enumerate(saved):
        frequencies.append(mode["omega"])
        c=np.asarray(mode["coefficients"])
        for i,arm in enumerate(model.arms):
            q[i,:,:,j]=(arm.basis(mode["omega"],d["xi"]*arm.length)@c[i*6:i*6+6])[:,:3]
    grams=overlap(q,q,properties,old["lengths"],d["weights"])
    q=q/np.sqrt(np.diag(grams))[None,None,None,:]
    analytical=dict(q=q,omega=np.array(frequencies),weights=d["weights"])
    comparison=compare_spaces(analytical,d,properties,old["lengths"],config["policy"])
    frequency=float(max(comparison["relative_frequency"]))
    passed=(frequency<=config["policy"]["reference_frequency_rtol"] and
            comparison["maximum_one_minus_MAC"]<=config["policy"]["one_minus_mac_tol"])
    return dict(case=case,beta_deg=beta,status="PASS" if passed else "NUMERICAL_UNRESOLVED",
                max_relative_frequency=frequency,comparison=comparison,
                new_boundary_root_solves=0,source="saved full-precision determinant roots and coefficients")


def benchmark(env, source, screen, out):
    config, old, props, ref=env
    path=out/"benchmark.json"
    if path.exists():
        return read(path)
    jobs=[(c,b,1.) for b in (15,30,45,60,75) for c in ("MIRROR","SAME_LOCAL")]
    jobs+=[("BASE",30,1.)]
    jobs+=[(c,0,1.) for c in CASES]+[(c,30,100.) for c in CASES]+[(c,30,None) for c in CASES]
    report=dict(status="IN_PROGRESS",points=[],references=[])
    for c,b,k in jobs:
        p,new=get_point(c,b,k,env,out)
        report["points"].append(dict(identity=p["identity"],status=p["status"],
                                    seconds=sum(r["seconds"] for r in p["orders"])))
        if p["status"]!="PASS":
            report.update(status="PARTIAL_PASS",stop_reason=f"Benchmark {c,b,k} failed gates")
            write(path,report);return report
        if k==1. and b in (15,30,45,60,75):
            saved=(source["spectral"]["BASE_SPRING"]["modes"] if c=="BASE" else
                   next(v["modes"] for v in screen["points"] if v["case"]==c and v["beta_deg"]==b))
            comparison=reference_check(p,saved,env,out)
            report["references"].append(comparison)
            if comparison["status"]!="PASS":
                report.update(status="PARTIAL_PASS",stop_reason=f"Reference {c,b,k} failed")
                write(path,report);return report
        write(out/"benchmark_progress.json",report)
        print("benchmark",c,b,k,"PASS",flush=True)
    import numpy as np
    times=[p["seconds"] for p in report["points"]]
    angles,springs=grids(config)
    unique=3*(len(angles)+len(springs)-1)+3
    estimate=float(np.median(times)*(unique-len(jobs)))
    report.update(status="PASS" if estimate<=config["policy"]["runtime_review_seconds"] else "RUNTIME_REVIEW",
                  median_seconds_per_pair=float(np.median(times)),remaining_pairs=unique-len(jobs),
                  estimated_remaining_seconds=estimate,unique_total_points=unique,
                  no_additional_high_frequency_inventory=True)
    write(path,report)
    return plain(report)


def rows_for(point, family, ref):
    import numpy as np
    case,beta,kappa=point["identity"];r=point["orders"][-1]
    rows=[]
    for j in range(7):
        p=r["parity"][j];eta=p["eta"] if p else None
        protected=(eta==-1 and case in ("BASE","MIRROR"))
        delta=r["delta_psi_raw"][j];mass=r["mass"][j];w=r["omega"][j]
        raw=delta**2/mass;participation=0. if protected else raw
        # RIGID suppresses delta kinematically; it is not a finite-spring predictor.
        if kappa is None:participation=0.
        omega=w*ref["t_ref"]
        g=.5*omega*(ref["D"]/ref["length"])*participation/w**2
        rows.append(dict(family=family,case=case,beta_deg=beta,kappa_theta=kappa,
                         position=j+1,guard=j==6,omega=w,Omega=omega,Lambda=np.sqrt(omega),
                         mass=mass,delta_psi_raw=delta,P_raw=raw,P=participation,G=g,
                         eta=eta,structural_zero=protected or kappa is None,
                         zero_reason="RIGID" if kappa is None else ("reflection_eta_minus" if protected else ""),
                         parity_residual=p["residual"] if p else None,
                         cluster=next("-".join(map(str,gp)) for gp in r["clusters"] if j+1 in gp),
                         ritz_frequency_change=point["convergence"]["relative_frequency"][j],
                         status=point["status"],point_id=key(case,beta,kappa)))
    return rows


def audit_tables(angle_rows,spring_rows,rigid_rows,config):
    import numpy as np
    findings=[];protected=[];jumps=[]
    for case in CASES:
        data=[r for r in spring_rows if r["case"]==case]
        kappas=sorted({r["kappa_theta"] for r in data})
        w=np.array([[next(r["omega"] for r in data if r["kappa_theta"]==k and r["position"]==j)
                     for j in range(1,8)] for k in kappas])
        wr=np.array([r["omega"] for r in rigid_rows if r["case"]==case])
        violation=float(max(0.,np.max(-np.diff(w,axis=0)/w[:-1]),np.max(w/wr-1)))
        findings.append(dict(case=case,monotonicity_and_rigid_violation=violation,
                             status="PASS" if violation<=config["policy"]["monotonicity_rtol"] else "NUMERICAL_UNRESOLVED"))
        if case in ("BASE","MIRROR"):
            initial=[r for r in data if r["kappa_theta"]==0 and r["eta"]==-1]
            for item in initial:
                matches=[]
                for k in kappas:
                    candidates=[r for r in data if r["kappa_theta"]==k and r["eta"]==-1]
                    if candidates:
                        match=min(candidates,key=lambda r:abs(r["omega"]/item["omega"]-1))
                        if abs(match["omega"]/item["omega"]-1)<=config["policy"]["protected_frequency_rtol"]:
                            matches.append(match)
                protected.append(dict(case=case,omega=item["omega"],Lambda=item["Lambda"],eta=-1,
                    count=len(matches),expected_count=len(kappas),all_points_resolved=len(matches)==len(kappas),
                    sorted_positions=sorted({r["position"] for r in matches}),
                    maximum_relative_variation=max([abs(r["omega"]/item["omega"]-1) for r in matches],default=0),
                    plotted_somewhere=any(r["position"]<=6 for r in matches)))
        for family,rows,param in (("angle",angle_rows,"beta_deg"),("spring",spring_rows,"kappa_theta")):
            for j in range(1,7):
                seq=sorted([r for r in rows if r["case"]==case and r["position"]==j],key=lambda r:r[param])
                for a,b in zip(seq,seq[1:]):
                    change=abs(b["omega"]/a["omega"]-1)
                    if change>config["neighbour_audit"]["relative_jump_trigger"]:
                        jumps.append(dict(family=family,case=case,position=j,left=a[param],right=b[param],
                                          relative_change=change,
                                          classification="Ritz 36/44 independently agree at both real points; retained, no interpolation"))
    return dict(status="PASS" if all(r["status"]=="PASS" for r in findings) else "PARTIAL_PASS",
                spring_order=findings,protected_families=protected,neighbour_triggers=jumps,
                continuum_completeness="NOT_CLAIMED")


def render(out,destination):
    """Only saved tables; exactly two figures, with no solver or model access."""
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FixedLocator, FixedFormatter
    destination.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    protected=read(out/"audit.json")["protected_families"]
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"pdf.fonttype":42,
                         "axes.facecolor":"white","figure.facecolor":"white"})
    colors=["#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9"]
    for family,filename in (("angle","angle_spectrum.csv"),("spring","spring_spectrum.csv")):
        with (out/filename).open(encoding="utf-8") as stream:rows=list(csv.DictReader(stream))
        fig,axs=plt.subplots(1,3,figsize=(12,4.1),sharey=True,layout="constrained")
        ymax=max(float(r["Lambda"]) for r in rows if int(r["position"])<=6)
        for ax,case in zip(axs,CASES):
            for j,color in enumerate(colors,1):
                data=[r for r in rows if r["case"]==case and int(r["position"])==j]
                parameter="beta_deg" if family=="angle" else "kappa_theta"
                x=np.array([float(r[parameter]) for r in data])
                y=np.array([float(r["Lambda"]) if r["status"]=="PASS" else np.nan for r in data])
                ax.plot(x if family=="angle" else np.log10(1+x),y,color=color,lw=1.25,label=rf"$\Lambda_{j}$")
            ax.set_title(case.replace("_"," "),fontsize=11)
            ax.set_ylim(0,ymax*1.045)
            ax.grid(alpha=.18,lw=.5)
            ax.spines[["top","right"]].set_visible(False)
            if family=="angle":
                ax.set_xlim(0,90);ax.set_xticks(range(0,91,15));ax.set_xlabel(r"$\beta,\ ^\circ$")
            else:
                ticks=np.array([0,1,3,10,30,100])
                ax.set_xlim(0,np.log10(101));ax.xaxis.set_major_locator(FixedLocator(np.log10(1+ticks)))
                ax.xaxis.set_major_formatter(FixedFormatter([str(x) for x in ticks]))
                ax.set_xlabel(r"$\kappa_\theta$")
                if case in ("BASE","MIRROR"):
                    lines=[p for p in protected if p["case"]==case and p["all_points_resolved"] and p["plotted_somewhere"]]
                    for p in lines:
                        ax.axhline(p["Lambda"],color=".25",ls=(0,(3,3)),lw=.7,alpha=.75)
                    if lines:ax.text(.97,.035,r"$\eta=-1$  --",transform=ax.transAxes,ha="right",fontsize=9,color=".25")
        axs[0].set_ylabel(r"$\Lambda=\sqrt{\omega t_{\mathrm{ref}}}$")
        handles,labels=axs[0].get_legend_handles_labels()
        fig.legend(handles,labels,ncol=6,loc="outside lower center",frameon=False,handlelength=2)
        for ext in ("png","pdf"):
            fig.savefig(destination/f"{family}_spectral_curves.{ext}",dpi=240)
        plt.close(fig)
    files=list(destination.glob("*.*"))
    manifest=dict(status="PASS",new_eigenvalue_solves=0,figure_count=2,render_seconds=time.perf_counter()-started,
                  spectrum_semantics="sorted_positions",source_tables={n:sha(out/n) for n in
                  ("angle_spectrum.csv","spring_spectrum.csv","audit.json")},
                  renderer_sha256=sha(__file__),matplotlib=matplotlib.__version__,
                  artifacts={p.name:sha(p) for p in files if p.suffix in (".png",".pdf")})
    write(destination/"figure_manifest.json",manifest)
    return manifest


def validate(out,identity=None):
    m=read(out/"manifest.json")
    if identity is not None:assert m["identity"]==identity
    assert all(sha(out/n)==h for n,h in m["artifacts"].items())
    return read(out/"result.json")


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark",action="store_true")
    parser.add_argument("--compute",action="store_true")
    parser.add_argument("--plot-only",action="store_true")
    args=parser.parse_args()
    if args.plot_only:
        out=OUTPUT/read(OUTPUT/"current.json")["fingerprint"]
        validate(out)
        result=render(out,OUTPUT/"plot_only")
        print(json.dumps(dict(status=result["status"],mode="plot_only",new_eigenvalue_solves=0)))
        return
    config,old,source,props,ref,screen,identity,fp=context()
    out=OUTPUT/fp;out.mkdir(parents=True,exist_ok=True)
    if (out/"manifest.json").exists():
        result=validate(out,identity)
        print(json.dumps(dict(status=result["status"],cache="VALIDATED",fingerprint=fp,new_eigenvalue_solves=0)))
        return
    write(out/"contract.json",identity)
    write(out/"materials.json",material_contract(props))
    if not (out/"ledger.json").exists():
        write(out/"ledger.json",dict(ritz_eigenproblems=0,boundary_root_solves=0,complex_solves=0))
    if not(args.benchmark or args.compute):
        print(json.dumps(dict(status="PREFLIGHT",fingerprint=fp,new_eigenvalue_solves=0)));return
    env=(config,old,props,ref)
    started=time.perf_counter()
    bench=benchmark(env,source,screen,out)
    print("benchmark summary",bench["status"],bench.get("estimated_remaining_seconds"),flush=True)
    if bench["status"]!="PASS" or not args.compute:return
    angles,springs=grids(config)
    angle_rows=[];spring_rows=[];rigid_rows=[]
    reused=0;new=0;stop=None
    for family,values in (("angle",angles),("spring",springs),("rigid",[None])):
        for case in CASES:
            for value in values:
                b,k=(value,1.) if family=="angle" else (30,value)
                p,created=get_point(case,b,k,env,out)
                new+=int(created);reused+=int(not created)
                if p["status"]!="PASS":
                    stop=f"{family} {case} {value}: {p['gates']}";break
                rows=rows_for(p,family,ref)
                {"angle":angle_rows,"spring":spring_rows,"rigid":rigid_rows}[family].extend(rows)
                if created and new%20==0:print("completed new points",new,flush=True)
            if stop:break
        if stop:break
    if stop:
        write(out/"stop.json",dict(status="PARTIAL_PASS",reason=stop,ledger=read(out/"ledger.json")))
        print(stop);return
    table(out/"angle_spectrum.csv",angle_rows)
    table(out/"spring_spectrum.csv",spring_rows)
    table(out/"rigid_spectrum.csv",rigid_rows)
    audit=audit_tables(angle_rows,spring_rows,rigid_rows,config);write(out/"audit.json",audit)
    points=[read(p) for p in (out/"points").glob("*.json")]
    result=dict(status=audit["status"],reference=ref,point_counts=dict(angle=273,spring=183,rigid=3,
                unique=len(points),shared_angle_spring=3,benchmark_points=20),
                production_created_points=new,production_reused_points=reused,
                benchmark=bench,solve_ledger=read(out/"ledger.json"),audit=audit,
                maximum_frequency_convergence=max(max(p["convergence"]["relative_frequency"]) for p in points),
                maximum_one_minus_MAC=max(p["convergence"]["maximum_one_minus_MAC"] for p in points),
                maximum_metrics={n:max(r["metrics"][n] for p in points for r in p["orders"]) for n in points[0]["orders"][0]["metrics"]},
                elapsed_this_invocation_seconds=time.perf_counter()-started,
                total_point_solver_seconds=sum(r["seconds"] for p in points for r in p["orders"]),
                historical_statuses_unchanged=True,continuum_completeness="NOT_CLAIMED")
    write(out/"result.json",result)
    write(out/"manifest.json",dict(identity=identity,git_head=subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
         artifacts={p.relative_to(out).as_posix():sha(p) for p in out.rglob("*") if p.is_file() and p.name!="manifest.json"}))
    write(OUTPUT/"current.json",dict(fingerprint=fp,status=result["status"]))
    render(out,out/"figures")
    print(json.dumps(dict(status=result["status"],fingerprint=fp,ledger=result["solve_ledger"],points=result["point_counts"])))


if __name__=="__main__":
    main()
