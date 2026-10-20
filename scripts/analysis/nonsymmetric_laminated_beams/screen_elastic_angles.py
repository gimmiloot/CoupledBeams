"""NSL-3A-ANGLE: fixed elastic screening, saved centre, eight missing points.

Script Proliferation Control: diagnostic-only saved/new-data workflow, not a
new solver or a parameter fork of production physics. The historical NSL-2
verification runner also runs Ritz/straight/axis controls and cannot serve
this missing-only cache contract without changing its frozen source hash.
Reuse its prefix/count/energies, NSL-2R recovery and NSL-3A provenance.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import numpy as np
import scipy
from scipy.linalg import block_diag
from scripts.analysis.nonsymmetric_laminated_beams import diagnose_mode_recovery as recovery
from scripts.analysis.nonsymmetric_laminated_beams import verify_kelvin_voigt_pilot as pilot
from scripts.lib import nonsymmetric_laminated_beam as single
from scripts.lib import nonsymmetric_laminated_joint as joint
from scripts.lib.nonsymmetric_laminated_kelvin_voigt import reflection_rows
from scripts.analysis.check_contraction_joint_sensitivity import plain
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv

CONFIG = ROOT/"data/input/nonsymmetric_laminated_beam_elastic_angle_screening.json"
OUTPUT = ROOT/"results/nonsymmetric_laminated_beams/nsl3a_elastic_angle_screening"
read = pilot.read


def context(config):
    old, data, properties, reference, manifests = pilot.context(read(pilot.CONFIG))
    fp, ident = pilot.identity(read(pilot.CONFIG), manifests)
    paths = list((ROOT/"results/nonsymmetric_laminated_beams").rglob("manifest.json"))
    matches = [p for p in paths if p.parent.name == fp]
    assert fp == config["source_bundles"][-1] and len(matches) == 1
    path = matches[0]; manifest = read(path)
    assert manifest["identity"] == ident
    assert all(sha(path.parent/n) == h for n, h in manifest["artifact_hashes"].items())
    manifests[fp] = dict(path=path.relative_to(ROOT).as_posix(), sha256=sha(path))
    assert set(manifests) == set(config["source_bundles"])
    assert config["variant"] == old["variant"] == "PROJECT_REDDY_COUPLED"
    assert config["kappa_theta"] == old["main_kappa"] == 1
    assert config["reuse_beta_deg"] == old["main_beta_deg"] == 30
    assert config["c_theta_in_operator"] == 0
    assert config["parameter_grid"] == [15, 30, 45, 60, 75]
    assert config["cases"] == ["MIRROR", "SAME_LOCAL"]
    assert config["K_plot"] == 6 and config["K_guard"] == 7
    assert config["illustrative_d_theta"] == read(pilot.CONFIG)["d_theta"]
    assert data["spectral"]["MIRROR_SPRING"]["status"] == "PASS"
    assert data["spectral"]["SAME_LOCAL_SPRING"]["status"] == "PASS"
    return old, data, properties, reference, manifests, read(path.parent/"result.json")


def baseline(manifests):
    path = OUTPUT/"baseline.json"
    if not path.exists():
        hashes = dict(read(pilot.OUTPUT/"baseline.json")["sha256"])
        for item in manifests.values():
            directory = (ROOT/item["path"]).parent
            for p in directory.rglob("*"):
                if p.is_file():
                    hashes[p.relative_to(ROOT).as_posix()] = sha(p)
        for p in (pilot.CONFIG, Path(pilot.__file__),
                  ROOT/"scripts/lib/nonsymmetric_laminated_kelvin_voigt.py",
                  ROOT/"tests/test_nonsymmetric_laminated_kelvin_voigt.py",
                  ROOT/"docs/nonsymmetric_laminated_beams/nsl3a_kelvin_voigt_pilot.md"):
            hashes[p.relative_to(ROOT).as_posix()] = sha(p)
        # Append-only scientific records: preserve all bytes preceding new records.
        memory = {}
        for name in ("decisions.md", "knowledge.md"):
            p = ROOT/"docs/memory"/name
            memory[p.relative_to(ROOT).as_posix()] = dict(size=p.stat().st_size, sha256=sha(p))
        write_json(path, dict(sha256=hashes, memory_prefix=memory))
    b = read(path)
    assert all(sha(ROOT/n) == h for n, h in b["sha256"].items())
    for n, v in b["memory_prefix"].items():
        assert hashlib.sha256((ROOT/n).read_bytes()[:v["size"]]).hexdigest() == v["sha256"]
    return b


def identity(config, manifests):
    baseline(manifests)
    files = [CONFIG, Path(__file__), Path(recovery.__file__),
             Path(single.__file__), Path(joint.__file__)]
    ident = dict(config=config, baseline_sha256=sha(OUTPUT/"baseline.json"),
        source_manifests=manifests,
        files={p.relative_to(ROOT).as_posix(): sha(p) for p in files},
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    return hashlib.sha256(json.dumps(ident, sort_keys=True).encode()).hexdigest()[:16], ident


def make_model(old, properties, ref, case, angle):
    return joint.CoupledBeams([properties[n] for n in old["cases"][case]],
        old["lengths"], np.deg2rad(angle), joint.Joint("SPRING", ref["k_theta"]))


def relative(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(a), np.linalg.norm(b), 1e-300))


def structural(model, case, policy):
    s = joint.REFLECTION
    h = block_diag(*(single.state_matrix(.37j, a) for a in model.properties))
    operator_error = relative(h@s, s@h)
    errors = []
    for omega in (.37, .83):
        b = model.boundary(omega, normalized=False)
        errors.append(relative(b@s, reflection_rows(model.beta)@b))
    j = joint.joint_matrix(model.beta, model.joint)
    jump = np.zeros(12); jump[[2, 8]] = [1, -1]
    rank = int(np.linalg.matrix_rank(j/np.linalg.norm(j, axis=1)[:, None]))
    # This is an exact algebraic projector statement, not rounding a mode.
    protected_jump = float(np.linalg.norm(jump@(np.eye(12)-s)/2))
    passed = rank == 6 and protected_jump == 0
    if case == "MIRROR":
        passed &= max(operator_error, *errors) <= policy["operator_tol"]
    else:
        passed &= operator_error > policy["operator_tol"]
    return dict(status="PASS" if passed else "FAIL", joint_rank=rank,
        reflection_operator_error=operator_error, reflection_boundary_errors=errors,
        protected_subspace_jump=protected_jump, symmetry_claim=case == "MIRROR")


def parity(mode, properties):
    transformed = []
    for arm in reversed(mode["arms"]):
        transformed.append({**arm, "states": arm["states"]@single.REVERSE})
    reflected = dict(arms=transformed)
    pair = joint.inner(mode, reflected, properties)/mode["mass"]
    eta = 1 if pair >= 0 else -1
    difference = dict(arms=[{**a, "states": b["states"]-eta*a["states"]}
                           for a, b in zip(mode["arms"], transformed)])
    error = np.sqrt(max(0., joint.inner(difference, difference, properties)/mode["mass"]))
    return dict(eta=eta, mass_error=float(error), mass_pairing=float(pair))


def indicators(omega, mass, delta, ref, d, protected=False):
    p_raw = delta**2/mass
    p = 0. if protected else p_raw
    s = ref["D"]/ref["length"]*p/omega**2
    om = omega*ref["t_ref"]
    g = .5*om*s
    alpha = ref["c_theta"]*p/2
    zeta = d*g
    error = abs(zeta-alpha/omega)/max(abs(zeta), abs(alpha/omega), 1e-300)
    return dict(P=p, P_numeric=p_raw, s=s, G=g, alpha_pred=alpha,
                zeta_pred=zeta, Omega=om, predictor_algebra_relative=error)


def mode_checks(model, omega, coefficient, data, oldpolicy, policy, case):
    b = model.boundary(omega)
    null = np.linalg.norm(b@coefficient)/(np.linalg.norm(b)*np.linalg.norm(coefficient))
    ends = np.concatenate([a["end_states"][1] for a in data["arms"]])
    peaks = np.concatenate([np.max(abs(a["states"]), axis=0) for a in data["arms"]])
    j = joint.joint_matrix(model.beta, model.joint)
    rows = abs(j@ends)/np.maximum(abs(j)@peaks, 1e-30)
    second = ends[11]-model.joint.k_theta*(ends[2]-ends[8])
    second_scaled = abs(second)/(peaks[11]+model.joint.k_theta*(peaks[2]+peaks[8]))
    reflect = parity(data, model.properties) if case == "MIRROR" else None
    gates = recovery.physical_gates(data, oldpolicy)
    gates.update(null=float(null) <= policy["null_tol"],
        singular=data["singular_ratio"] <= policy["null_tol"],
        second_moment=float(second_scaled) <= policy["physical_tol"],
        simple_kernel=1/data["nonzero_condition"] >= policy["simple_sigma_min"],
        reflection=reflect is None or reflect["mass_error"] <= policy["symmetry_shape_tol"])
    return dict(null_relative=float(null), physical_joint_rows=rows,
        second_moment_physical=float(second), second_moment_scaled=float(second_scaled),
        reflection=reflect, gates=gates)


def save_profiles(path, model, modes, nodes):
    rows = []
    for mode in modes:
        for i, (arm, c) in enumerate(zip(model.arms, np.array(mode["coefficients"]).reshape(2, 6)), 1):
            x = np.linspace(0, arm.length, nodes)
            for xx, y in zip(x, arm.basis(mode["omega"], x)@c):
                rows.append(dict(position=mode["position"], arm=i, x=float(xx),
                                 **dict(zip(("u", "w", "psi", "N", "Q", "M"), y))))
    write_csv(path, rows)


def point(config, old, source, model, ref, case, angle, out, ident, neighbour):
    name = f"{case}_{angle:02d}"
    path = out/(name+".json")
    count_path = out/(name+"_roots.json")
    if path.exists():
        p = read(path)
        assert p["identity"] == ident
        if p["status"] == "PASS":
            return p, True
        assert p["status"] == "IN_PROGRESS" and count_path.exists(), \
            "Unresolved/incomplete root search cannot be automatically replayed"
    reused = angle == config["reuse_beta_deg"]
    count = {}
    started = time.perf_counter()
    p = dict(identity=ident, case=case, beta_deg=angle, status="IN_PROGRESS",
        source="NSL2_SAVED" if reused else "NEW_ELASTIC", modes=[],
        neighbour_predictor=neighbour,
        new_roots=0, new_prefix_calls=0)
    try:
        if reused:
            saved = source["spectral"][case+"_SPRING"]
            roots = saved["count"]["roots"]; count = saved["count"]
            p["source_form_reference"] = dict(bundle="61842957b5e4d8d1",
                json_pointer=f"/spectral/{case}_SPRING/modes")
        elif count_path.exists():
            record = read(count_path)
            assert record["identity"] == ident and record["count"]["status"] == "PASS"
            count = record["count"]; roots = count["roots"]
            p.update(new_roots=len(roots), new_prefix_calls=1)
        else:
            # Existing NSL prefix/count only; deliberately no Ritz, mode tracking
            # or extra root above the completeness guard.
            p["new_prefix_calls"] = 1
            write_json(path, plain(p))
            roots = model.prefix(old["policy"], count)
            p["new_roots"] = len(roots)
            write_json(count_path, dict(identity=ident, count=plain(count)))
        p["count"] = count
        fields = []
        for i, root in enumerate(roots):
            omega = root["omega"]
            if reused:
                original = saved["modes"][i]
                c = np.array(original["coefficients"])
                # Evaluate the SAVED function; no kernel extraction, new
                # normalization or new mode shape. Preserve every coefficient.
                data = model.evaluate(omega, c, old["policy"]["mode_quadrature"])
                for key in ("quadrature_relative", "singular_ratio", "nonzero_condition"):
                    data[key] = original[key]
                assert abs(data["mass"]-original["mass"]) <= config["policy"]["mass_gram_tol"]
                assert data["delta_psi"] == original["delta_psi"]
            else:
                c, data = recovery.recover_mode(model, omega, old["policy"])
            checks = mode_checks(model, omega, c, data, old["policy"], config["policy"], case)
            gap = min(abs(r["omega"]-omega) for k, r in enumerate(roots) if k != i)/omega
            checks["gates"]["isolated"] = gap >= config["policy"]["simple_relative_gap_min"]
            protected = case == "MIRROR" and checks["reflection"]["eta"] == -1
            metrics = indicators(omega, data["mass"], data["delta_psi"], ref,
                                 config["illustrative_d_theta"], protected)
            checks["gates"]["predictor"] = metrics["predictor_algebra_relative"] <= config["policy"]["predictor_algebra_rtol"]
            row = dict(position=i+1, omega=omega, coefficients=c,
                mass=data["mass"], delta_psi=data["delta_psi"],
                structural_zero=protected, **metrics, **checks,
                physical=recovery.metrics(data), stiffness_components=data["stiffness_components"],
                spring_energy=data["spring_energy"], relative_neighbour_gap=gap,
                status="PASS" if all(checks["gates"].values()) else "NUMERICAL_UNRESOLVED")
            p["modes"].append(plain(row)); fields.append(data)
            if row["status"] != "PASS":
                raise ArithmeticError(f"{name} position {i+1}: {checks['gates']}")
        gram = np.array([[joint.inner(a, b, model.properties) for b in fields] for a in fields])
        p["mass_gram"] = gram
        p["mass_gram_error"] = float(np.max(abs(gram-np.eye(7))))
        assert p["mass_gram_error"] <= config["policy"]["mass_gram_tol"], "Mass orthogonality"
        assert count["status"] == "PASS" and [r["counts"] for r in roots] == [[i, i+1] for i in range(7)]
        p.update(status="PASS", elapsed_seconds=time.perf_counter()-started)
        if not reused:
            save_profiles(out/(name+"_profiles.csv"), model, p["modes"], config["policy"]["profile_nodes"])
    except Exception as error:
        p.update(status="NUMERICAL_UNRESOLVED", stop_reason=f"{type(error).__name__}: {error}",
                 count=count, new_roots=0 if reused else len(count.get("roots", [])))
        write_json(path, plain(p))
        raise
    write_json(path, plain(p))
    return plain(p), False


def global_pair(left, right, beta_left, beta_right, properties):
    """Common physical metric for THIS equal-density family (j1=0).

    Coordinates x_i/L_i and arm labels coincide. Rotations are about the same
    global normal. This is a geometric overlap, not angle continuation.
    """
    total = 0.
    for i, (a, b, p) in enumerate(zip(left["arms"], right["arms"], properties), 1):
        if p.j1 != 0:
            raise ValueError("This diagnostic physical metric requires the fixed j1=0 family")
        assert np.array_equal(a["x"], b["x"])
        ga = joint.global_kinematic_map(beta_left, i)
        gb = joint.global_kinematic_map(beta_right, i)
        qa = a["states"][:, :3]@ga.T; qb = b["states"][:, :3]@gb.T
        total += np.sum(a["weights"]*(p.m*(qa[:, 0]*qb[:, 0]+qa[:, 1]*qb[:, 1])
                                     +p.j2*qa[:, 2]*qb[:, 2]))
    return float(total)


def assignments(mac, policy, selected=None):
    rows = []
    for i in range(len(mac)) if selected is None else selected:
        j = int(np.argmax(mac[i]))
        row_gap = mac[i, j]-max(np.delete(mac[i], j))
        col_gap = mac[i, j]-max(np.delete(mac[:, j], i))
        passed = (mac[i, j] >= policy["mac_min"] and row_gap >= policy["mac_margin"]
                  and col_gap >= policy["mac_margin"] and np.argmax(mac[:, j]) == i)
        rows.append(dict(reference_position=i+1, target_position=j+1,
            MAC=float(mac[i, j]), row_margin=float(row_gap), column_margin=float(col_gap),
            status="RESOLVED" if passed else "UNRESOLVED"))
    return rows


def postprocess(config, old, properties, ref, points, kvsource):
    table = []; guard = []; fields = {}; models = {}
    for p in points:
        key = (p["case"], p["beta_deg"])
        model = make_model(old, properties, ref, *key); models[key] = model
        fields[key] = [model.evaluate(m["omega"], np.array(m["coefficients"]), old["policy"]["mode_quadrature"])
                       for m in p["modes"][:6]]
        for m in p["modes"]:
            row = dict(beta_deg=p["beta_deg"], case=p["case"], position=m["position"],
                **{n:m[n] for n in ("omega", "Omega", "mass", "delta_psi", "P", "P_numeric", "s", "G",
                                    "alpha_pred", "zeta_pred", "structural_zero")},
                eta=m["reflection"]["eta"] if m["reflection"] else None,
                spectral_status=p["count"]["status"], modal_status=m["status"], source=p["source"])
            (table if m["position"] <= 6 else guard).append(row)
    reuse_checks = []
    for state in kvsource["states"]:
        row = next(r for r in table if r["case"] == state["case"] and
                   r["beta_deg"] == 30 and r["position"] == state["elastic_position"])
        expected = state["alpha_pred"]
        if row["structural_zero"]:
            error = 0.
            assert row["alpha_pred"] == expected == 0
        else:
            error = abs(row["alpha_pred"]/expected-1)
        reuse_checks.append(dict(id=state["id"], alpha_pred=row["alpha_pred"],
            expected=expected, relative_error=error,
            status="PASS" if error <= config["policy"]["saved_predictor_rtol"] else "FAIL"))
    macs = []
    def match(left_key, right_key, kind, selected=None):
        left, right = fields[left_key], fields[right_key]
        ma, mb = models[left_key], models[right_key]
        values = np.array([[global_pair(a, b, ma.beta, mb.beta, ma.properties)**2/(a["mass"]*b["mass"])
                           for b in right] for a in left])
        macs.append(dict(kind=kind, reference=list(left_key), target=list(right_key),
            MAC=values, assignments=assignments(values, config["policy"], selected)))
    for angle in config["parameter_grid"]:
        match(("MIRROR", angle), ("SAME_LOCAL", angle), "SAME_ANGLE_CROSS_LAYOUT")
        if angle != 30:
            for case in config["cases"]:
                match((case, 30), (case, angle), "BETA30_REFERENCE_ONLY", [3, 5])
    neighbours = []
    for case in config["cases"]:
        for lo, hi in zip(config["parameter_grid"][:-1], config["parameter_grid"][1:]):
            a = next(p for p in points if (p["case"], p["beta_deg"]) == (case, lo))
            b = next(p for p in points if (p["case"], p["beta_deg"]) == (case, hi))
            ratios = [max(x["omega"]/y["omega"], y["omega"]/x["omega"])
                      for x, y in zip(a["modes"][:6], b["modes"][:6])]
            neighbours.append(dict(case=case, angles=[lo, hi], frequency_ratios=ratios,
                status="PASS" if max(ratios) <= config["neighbour_audit"]["max_adjacent_frequency_ratio"] else "UNRESOLVED"))
    maxima = []
    for case in config["cases"]:
        rows = [r for r in table if r["case"] == case]
        for name in ("P", "alpha_pred", "zeta_pred"):
            maximum = max(rows, key=lambda r:r[name])
            centre = max((r for r in rows if r["beta_deg"] == 30), key=lambda r:r[name])
            same_position = next(r for r in rows if r["beta_deg"] == 30 and r["position"] == maximum["position"])
            maxima.append(dict(case=case, metric=name, maximum=maximum,
                beta30_max=centre[name], ratio_to_beta30_max=maximum[name]/centre[name],
                beta30_same_position=same_position[name],
                ratio_to_beta30_same_position=maximum[name]/same_position[name] if same_position[name] else None,
                semantics="grid maximum / sorted-position comparison, not branch identity"))
    return plain(dict(table=table, guard=guard, beta30_checks=reuse_checks,
                      identity_diagnostics=macs, neighbour_audit=neighbours, maxima=maxima))


def plot_table(table, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle
    positive = [r["zeta_pred"] for r in table if not r["structural_zero"]]
    norm = LogNorm(min(positive), max(positive))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout="constrained")
    angles = [15, 30, 45, 60, 75]
    for ax, case in zip(axes, ("MIRROR", "SAME_LOCAL")):
        a = np.array([[next(r["zeta_pred"] for r in table if
            (r["case"], r["beta_deg"], r["position"]) == (case, angle, j))
            for j in range(1, 7)] for angle in angles])
        im = ax.imshow(np.ma.masked_equal(a, 0), norm=norm, cmap="viridis", aspect="auto")
        for i in range(5):
            for j in range(6):
                if a[i,j] == 0:
                    ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor="#eeeeee",hatch="//",edgecolor="#bbbbbb"))
                    ax.text(j,i,"0*",ha="center",va="center",fontsize=9,color="black")
                else:
                    ax.text(j,i,f"{a[i,j]:.1e}",ha="center",va="center",fontsize=8,
                            color="black" if norm(a[i,j])>.65 else "white")
        ax.set(xticks=range(6), xticklabels=range(1,7), yticks=range(5), yticklabels=angles,
               xlabel="Sorted elastic position", ylabel="Angle beta (degrees)", title=case)
    fig.colorbar(im, ax=axes, shrink=.85, label="Elastic predictor zeta (illustrative d = 0.001)")
    fig.suptitle("NSL-3A-ANGLE: elastic predictors only\n0* = exact protected subspace; no complex roots")
    fig.savefig(out/"angle_screening.png", dpi=180)
    fig.savefig(out/"angle_screening.pdf")
    plt.close(fig)


def validate_cache(out, identity_value=None):
    manifest = read(out/"manifest.json")
    if identity_value is not None:
        assert manifest["identity"] == identity_value
    assert all(sha(out/n) == h for n, h in manifest["artifact_hashes"].items())
    return read(out/"result.json")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true", help="Compute only missing fixed points")
    parser.add_argument("--plot-only", action="store_true", help="No matrix assembly, SVD or root calculation")
    args = parser.parse_args()
    if args.plot_only:
        out = OUTPUT/read(OUTPUT/"current.json")["fingerprint"]
        result = validate_cache(out)
        # Render outside the immutable bundle; do not rewrite hashed figures.
        render = OUTPUT/"plot_only"; render.mkdir(exist_ok=True)
        plot_table(result["table"], render)
        print(json.dumps(dict(status=result["status"], mode="plot_only", new_eigenvalue_solves=0)))
        return
    config = read(CONFIG)
    old, data, properties, ref, manifests, kvsource = context(config)
    fp, ident = identity(config, manifests); out = OUTPUT/fp
    if (out/"manifest.json").exists():
        result = validate_cache(out, ident)
        print(json.dumps(dict(status=result["status"], fingerprint=fp, cache="VALIDATED",
                              new_eigenvalue_solves=0, reused_points=len(result["points"]))))
        return
    out.mkdir(parents=True, exist_ok=True)
    checks = []
    for angle in config["parameter_grid"]:
        for case in config["cases"]:
            model = make_model(old, properties, ref, case, angle)
            checks.append(dict(beta_deg=angle, case=case, **structural(model, case, config["policy"])))
    write_json(out/"preflight_contract.json", plain(dict(identity=ident, reference=ref,
        numerical_policy=old["policy"], structural=checks)))
    assert all(c["status"] == "PASS" for c in checks), "Structural gate"
    if not args.compute:
        print(json.dumps(dict(status="PASS", stage="PREFLIGHT", fingerprint=fp, new_eigenvalue_solves=0)))
        return
    points = []; new_calls = 0
    result = dict(status="IN_PROGRESS", points=points, structural=checks, reference=ref)
    try:
        # Load the already completed centre first; then traverse the specified grid.
        order = [30, 15, 45, 60, 75]
        for angle in order:
            for case in config["cases"]:
                model = make_model(old, properties, ref, case, angle)
                candidates = [p for p in points if p["case"] == case]
                near = min(candidates, key=lambda p:abs(p["beta_deg"]-angle)) if candidates else None
                neighbour = None if near is None else dict(beta_deg=near["beta_deg"],
                    frequencies=[m["omega"] for m in near["modes"]])
                p, cached = point(config, old, data, model, ref, case, angle, out, ident, neighbour)
                points.append(p)
                if not cached and angle != 30:
                    new_calls += 1
                assert sum(p["new_roots"] for p in points) <= config["max_new_roots"]
                write_json(out/"progress.json", plain(dict(completed=[(p["case"],p["beta_deg"]) for p in points],
                    total_new_roots=sum(p["new_roots"] for p in points))))
                print(case, angle, p["status"], "roots", p["new_roots"],
                      "max BC", max(m["physical"]["boundary_scaled"] for m in p["modes"]), flush=True)
        result.update(postprocess(config, old, properties, ref, points, kvsource))
        result["status"] = "PASS" if all(x["status"] == "PASS" for x in
            result["beta30_checks"]+result["neighbour_audit"]) else "PARTIAL_PASS"
        write_csv(out/"participation.csv", result["table"])
        write_csv(out/"guard_diagnostics.csv", result["guard"])
        write_json(out/"identity_diagnostics.json", result["identity_diagnostics"])
        plot_table(result["table"], out)
    except Exception as error:
        import traceback
        result.update(status="PARTIAL_PASS", stop_reason=f"{type(error).__name__}: {error}",
                      traceback=traceback.format_exc())
    saved_points = [read(p) for p in out.glob("*.json") if p.stem.startswith(("MIRROR_", "SAME_LOCAL_"))
                    and not p.stem.endswith("_roots")]
    result["solve_ledger"] = dict(new_spectral_points=sum(p["source"] == "NEW_ELASTIC" for p in saved_points),
        accepted_new_points=sum(p["source"] == "NEW_ELASTIC" and p["status"] == "PASS" for p in saved_points),
        new_elastic_roots=sum(p["new_roots"] for p in saved_points),
        reused_beta30_points=sum(p["source"] == "NSL2_SAVED" for p in saved_points),
        reused_beta30_roots=sum(len(p["modes"]) for p in saved_points if p["source"] == "NSL2_SAVED"),
        complex_roots=0, new_ritz_eigenproblems=0, completed_new_prefix_calls_this_invocation=new_calls)
    result["protected_file_count"] = len(baseline(manifests)["sha256"])
    write_json(out/"result.json", plain(result))
    write_json(out/"manifest.json", dict(identity=ident,
        git_head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        command=subprocess.list2cmdline([sys.executable,*sys.argv]),
        artifact_hashes={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"}))
    write_json(OUTPUT/"current.json", dict(fingerprint=fp, status=result["status"]))
    print(json.dumps(dict(status=result["status"], fingerprint=fp, ledger=result["solve_ledger"],
                          stop_reason=result.get("stop_reason"))))


if __name__ == "__main__":
    main()
