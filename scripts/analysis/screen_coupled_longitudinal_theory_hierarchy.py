"""Fixed eight-angle screening; sorted positions and fixed-beta geometry overlap.

New workflow contract: compare three theories, with comparator gates first.
Reuses verified MH/Timoshenko arms, geometry, independent energy root counts
and immutable MH references. No energy classification or modal continuation.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
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
import matplotlib
from scripts.analysis import verify_mindlin_herrmann_timoshenko_general_beta_joint as general
from scripts.analysis import verify_mindlin_herrmann_timoshenko_single_rod as single
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import coupled_longitudinal_comparators as reduced
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint

CONFIG = ROOT/"data/input/coupled_longitudinal_theory_hierarchy_screening.json"
OUTPUT = ROOT/"results/coupled_longitudinal_theory_hierarchy_screening"
VERSION = "fixed-beta-sorted-hierarchy-geometric-overlap-v1"
POLICY = dict(general.POLICY)  # unchanged accepted residual/search tolerances


def check_inputs(config_path=CONFIG):
    config = json.loads(config_path.read_text(encoding="utf-8"))
    base, model, length, checked, direct, frozen = general.check_inputs()
    s = model.section
    if (config["material"] != {"E": s.E, "rho": s.rho, "nu": s.nu, "kappa": s.K} or
        config["geometry"] != {"b": s.width, "h": s.thickness, "L": length, "split": .5} or
        config["beta_deg"] != [0, 5, 15, 30, 45, 60, 75, 90] or
        config["models"] != [*reduced.NAMES, "mindlin_herrmann"] or
        (config["K_plot"], config["K_guard"]) != (12, 13)):
        raise ValueError("This bounded screening accepts the declared G20/grid/prefix only")
    reference = ROOT/config["reference_bundle"]
    manifest = json.loads((reference/"manifest.json").read_text(encoding="utf-8"))
    if any(sha(reference/p) != h for p, h in manifest["artifact_hashes"].items()):
        raise ValueError("MH reference artifact hash mismatch")
    _, current = single.identity(base, checked)
    previous = dict(manifest["identity"]["arm_identity"])
    previous.pop("head", None)
    current.pop("head", None)
    if current != previous or any(sha(ROOT/p) != h for p, h in manifest["identity"]["files"].items()):
        raise ValueError("MH reference code/source/arm identity mismatch")
    reference_data = json.loads((reference/"result.json").read_text(encoding="utf-8"))
    if reference_data["statuses"]["MHTIM_GENERAL_BETA_JOINT_GATE"] != "PASS":
        raise ValueError("MH reference gate not accepted")
    return config, base, model, length, checked, frozen, reference_data


def identity(config, checked):
    paths = [Path(__file__), CONFIG, ROOT/"scripts/lib/coupled_longitudinal_comparators.py",
        ROOT/"scripts/lib/bishop_longitudinal.py", ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py",
        ROOT/"scripts/lib/mindlin_herrmann_timoshenko_joint.py",
        ROOT/"scripts/lib/reddy_inplane_geometry.py", single.CONFIG,
        ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
        Path(general.__file__), Path(single.__file__)]
    reference = ROOT/config["reference_bundle"]
    record = {"version": VERSION, "config": config, "policy": POLICY,
        "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        "source_hashes": checked, "reference_manifest": sha(reference/"manifest.json"),
        "beta0_manifest": sha(general.FROZEN/"manifest.json"),
        "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__}}
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()[:16], record


def tim_catalog(model, length, upper, base):
    """Min-max saturated fixed-arm catalog; choose SS half-gap before search."""
    k = mh.blocks(model)[1].spatial(upper/(2*math.pi))[0]["wavenumber_per_m"]
    index = math.ceil(k*length/math.pi-.9)+.9
    ceiling = math.sqrt(mh.blocks(model)[1].temporal(index*math.pi/length)[0]["omega_squared"])
    lower = POLICY["omega_min_pi_c0_L"]*math.pi
    roots, search = mh.finite_roots(model, length, "timoshenko", lower, ceiling, base["policy"])
    bound = mh.finite_count_upper_bound(model, length, ceiling, "timoshenko")
    if len(roots) != bound["upper_count"] or search["failed_intervals"]:
        raise ArithmeticError(f"Fixed Tim pole count not saturated: L={length}, {len(roots)}, {bound}")
    return {"poles": [r["omega"] for r in roots], "roots": roots, "upper": ceiling,
        "search": search, "count_bound": bound, "ss_ceiling_index": index, "status": "PASS"}


def solve_reduced(model, lengths, beta, variant, catalogs, upper, frames=None):
    active = joint.frames(beta) if frames is None else frames
    roots, search = reduced.roots(model, lengths, variant, active, .001*math.pi, upper, catalogs, POLICY)
    for n, root in enumerate(roots, 1):
        mode = reduced.mode(model, lengths, root["omega"], variant, active, POLICY["quadrature_order"])
        root.update(mode)
        root["sorted_position"] = n
        reduced.validate_mode(root, POLICY)
    # Gradient inertia is part of the Love mass inner product. No fractions.
    nodes, weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    gram = np.zeros((min(13, len(roots)),)*2)
    for arm, length in enumerate(lengths):
        x, weight = (nodes+1)*length/2, weights*length/2
        fields = np.array([reduced.state(model, length, r["omega"], np.array(r["coefficients"])[arm], x, variant) for r in roots[:13]])
        gradients = np.array([reduced.state(model, length, r["omega"], np.array(r["coefficients"])[arm], x, variant, 1)[:, 0] for r in roots[:13]])
        s = reduced.segment(model, length, variant)
        for i, first in enumerate(fields):
            for j, second in enumerate(fields):
                gram[i, j] += float(weight@(s.m*(first[:, 0]*second[:, 0]+first[:, 1]*second[:, 1])+model.coefficients["r"]*first[:, 2]*second[:, 2]+s.J*gradients[i]*gradients[j]))
    gram_error = float(np.max(abs(gram-np.eye(len(gram)))))
    if gram_error > POLICY["mass_gram_tol"]:
        raise ArithmeticError("Comparator gradient-mass orthogonality failed")
    return {"model": variant, "beta_deg": beta, "lengths": list(lengths),
        "roots": roots, "search": search, "mass_gram": gram.tolist(), "mass_gram_max_error": gram_error, "status": "PASS"}


def geometry_gates():
    rng = np.random.default_rng(1948)
    rows = []
    for beta in (0, 5, 45, 90):
        frames = joint.frames(beta)
        error = 0.
        for frame in frames:
            g = reduced.transform(frame)
            q, p = rng.normal(size=(2, 3))
            error = max(error, abs(p@(g.T@q)-(g@p)@q), float(np.max(abs(g.T@g-np.eye(3)))))
        rank = int(np.linalg.matrix_rank(reduced.joint_matrix(frames)))
        if error > POLICY["coordinate_duality_tol"] or rank != 6:
            raise ArithmeticError("Comparator geometry/duality/rank gate failed")
        rows.append({"beta_deg": beta, "rank": rank, "virtual_work_orthogonality_error": error,
            "transforms": [reduced.transform(f).tolist() for f in frames]})
    return rows


def sample_case(model, case, frames=None, order=None):
    """Common quadrature; physical d and theta kept dimensionally separate."""
    nodes, weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"] if order is None else order)
    active = joint.frames(case["beta_deg"]) if frames is None else frames
    common_weights, values, gradients = [], [], []
    for arm, (length, frame) in enumerate(zip(case["lengths"], active)):
        x, weight = (nodes+1)*length/2, weights*length/2
        common_weights.append(weight)
        arm_values, arm_gradients = [], []
        for r in case["roots"][:13]:
            coefficients = np.array(r["coefficients"])[arm]
            if case["model"] == "mindlin_herrmann":
                v = joint.arm_state(model, length, r["omega"], coefficients, x)
                derivative = joint.arm_state(model, length, r["omega"], coefficients, x, 1)
                fields = v[:, [0, 2, 3, 1]]  # u,w,theta,c
            else:
                v = reduced.state(model, length, r["omega"], coefficients, x, case["model"])
                derivative = reduced.state(model, length, r["omega"], coefficients, x, case["model"], 1)
                fields = np.column_stack((v[:, :3], np.zeros(len(x))))
            displacement = fields[:, :2]@frame.translation.T
            arm_values.append(np.column_stack((displacement, fields[:, 2:])))
            arm_gradients.append(derivative[:, 0])
        values.append(np.array(arm_values))
        gradients.append(np.array(arm_gradients))
    return np.concatenate(values, axis=1), np.concatenate(gradients, axis=1), np.concatenate(common_weights)


def symmetry_error(model, first, second, operation):
    nodes, weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    rows = []
    for a, b in zip(first["roots"][:13], second["roots"][:13]):
        vectors = [[], []]
        for arm, length in enumerate(first["lengths"]):
            x, weight = (nodes+1)*length/2, weights*length/2
            v = reduced.state(model, length, a["omega"], np.array(a["coefficients"])[arm], x, first["model"])
            refarm = 1-arm if operation == "swap" else arm
            r = reduced.state(model, length, b["omega"], np.array(b["coefficients"])[refarm], x, first["model"])
            if operation == "reflection":
                r *= reduced.MIRROR
            vectors[0].append((v[:, :3]*[1., 1., sum(first["lengths"])])*np.sqrt(weight)[:, None])
            vectors[1].append((r[:, :3]*[1., 1., sum(first["lengths"])])*np.sqrt(weight)[:, None])
        v, r = (np.concatenate(vector).ravel() for vector in vectors)
        sign = 1. if v@r >= 0 else -1.
        error = float(np.linalg.norm(v-sign*r)/np.linalg.norm(v))
        freq = abs(a["omega"]/b["omega"]-1)
        if freq > POLICY["symmetry_frequency_relative_tol"] or error > POLICY["symmetry_component_L2_tol"]:
            raise ArithmeticError("Comparator swap/reflection gate failed")
        rows.append({"position": a["sorted_position"], "relative_frequency_difference": freq, "kinematic_L2_relative": error})
    return {"status": "PASS", "rows": rows}


def comparator_gates(config, base, model, length, catalogs, upper):
    geometry = geometry_gates()
    direct_tim = tim_catalog(model, length, upper, base)
    report, reusable = {}, {}
    for variant in reduced.NAMES:
        splits = []
        s = reduced.segment(model, length, variant)
        direct_axial = reduced.axial_poles(model, length, variant, upper)
        direct = sorted(direct_axial+[w for w in direct_tim["poles"] if w < upper])
        for split in config["comparator_splits"]:
            lengths = (length*split, length*(1-split))
            case = solve_reduced(model, lengths, 0., variant, catalogs, upper)
            errors = [abs(r["omega"]/w-1) for r, w in zip(case["roots"], direct)]
            if len(case["roots"]) != len(direct) or max(errors) > POLICY["beta0_frequency_relative_tol"]:
                raise ArithmeticError(f"Direct/split recovery failed: {variant}, split={split}, {errors}")
            splits.append({"split": split, "max_frequency_relative_difference": max(errors), "case": case})
            if split == .5:
                reusable[(variant, 0)] = case
        canonical = solve_reduced(model, (.5, .5), 45., variant, catalogs, upper)
        swapped = solve_reduced(model, (.5, .5), 45., variant, catalogs, upper, joint.frames(45)[::-1])
        mirrored = solve_reduced(model, (.5, .5), 45., variant, catalogs, upper,
            joint.reflected_frames(joint.frames(45)))
        swap = symmetry_error(model, canonical, swapped, "swap")
        reflection = symmetry_error(model, canonical, mirrored, "reflection")
        reusable[(variant, 45)] = canonical
        report[variant] = {"status": "PASS", "geometry": geometry, "J_planar": s.J,
            "splits": splits, "swap": swap, "reflection": reflection,
            "symmetry_cases": [swapped, mirrored]}
        print("Comparator gates PASS:", variant, flush=True)
    return report, reusable


def mh_case(model, length, beta, reference, frozen, count, upper):
    if beta in (5, 45, 90):
        case = copy.deepcopy(next(c for c in reference["pilot"] if c["beta_deg"] == beta))
        case["origin"] = "immutable general-beta reference"
    elif beta == 0:
        split = next(s for s in frozen["splits"] if s["split"] == .5)
        records = sorted([r for block in ("mh", "timoshenko") for r in split[block]["roots"]], key=lambda r: r["omega"])[:13]
        roots = []
        for n, r in enumerate(records, 1):
            mode = joint.frame_mode(model, (.5, .5), 0., r["omega"], POLICY["quadrature_order"])
            record = {key: r[key] for key in ("omega", "frequency_hz", "bracket_omega", "singular_ratio", "nonzero_singular_condition")}
            record.update({"sorted_position": n, "coefficients": mode["coefficients"].tolist(), "diagnostics": mode["diagnostics"]})
            reduced.validate_mode(record, POLICY)
            roots.append(record)
        end = roots[-1]["bracket_omega"][1]
        nc, diagnostic = count(end, 0)
        start, sd = count(.001*math.pi, 0)
        if nc != 13 or start != 0:
            raise ArithmeticError("Cached MH beta0 prefix not count-certified")
        case = {"beta_deg": 0, "lengths": [.5, .5], "roots": roots, "origin": "immutable beta0 roots; fresh unchanged full-frame profile/residual check",
            "search": {"upper_count": nc, "lower_count": start, "range_omega": [.001*math.pi, end],
                "count_samples": [{"omega": end, "count": nc, **diagnostic}, {"omega": .001*math.pi, "count": start, **sd}],
                "failed_intervals": [], "status": "PASS", "root_evaluations": 0}, "status": "PASS"}
    else:
        case = general.solve_case(model, length, beta, count, upper)
        case["origin"] = "unchanged production MH frame solver"
    case["model"] = "mindlin_herrmann"
    if case["search"]["status"] != "PASS":
        raise ArithmeticError("MH reference search not verified")
    for r in case["roots"][:13]:
        reduced.validate_mode(r, POLICY)
    return case


def ranking(matrix):
    def direction(rows):
        records = []
        for n, row in enumerate(rows, 1):
            finite = sorted(((float(v), i+1) for i, v in enumerate(row) if v is not None), reverse=True)
            if not finite:
                records.append({"position": n, "status": "NOT_INFORMATIVE", "argmax": None})
                continue
            best, second = finite[0], finite[1] if len(finite) > 1 else (0., None)
            records.append({"position": n, "diagonal": row[n-1], "argmax": best[1],
                "best": best[0], "second_best": second[0], "margin": best[0]-second[0],
                "status": "POSITION_CORRESPONDENCE_NONDIAGONAL" if best[1] != n else "DIAGONAL"})
        return records
    return {"rows": direction(matrix), "columns": direction(list(map(list, zip(*matrix))))}


def diagnostics(config, model, cases):
    table, overlaps, gaps, contraction, summary, non_diagonal = [], [], [], [], [], []
    epsilon = np.finfo(float).eps
    for beta in config["beta_deg"]:
        active = {name: cases[(name, beta)] for name in config["models"]}
        fields = {}
        zeros = {}
        zero_scales = {}
        dcrows = []
        for name, case in active.items():
            values, ux, weight = sample_case(model, case)
            displacement = values[:12, :, :2].reshape(12, -1)
            theta = values[:12, :, 2]
            normd = np.sqrt(np.sum(displacement**2*np.repeat(weight, 2), axis=1))
            normtheta = np.sqrt(np.sum(theta**2*weight, axis=1))
            scales = np.array([config["numerical_zero_roundoff_multiplier"]*epsilon*
                r["diagnostics"]["nonzero_singular_condition"]*normd[k]/sum(case["lengths"])
                for k, r in enumerate(case["roots"][:12])])
            fields[name] = (displacement, theta, weight)
            zeros[name] = normtheta <= scales
            zero_scales[name] = scales.tolist()
            frequencies = [r["frequency_hz"] for r in case["roots"][:12]]
            gaps.extend({"model": name, "beta_deg": beta, "position": n+1, "next_position": n+2,
                "relative_adjacent_gap": b/a-1} for n, (a, b) in enumerate(zip(frequencies[:-1], frequencies[1:])))
            if name == "mindlin_herrmann":
                for k, r in enumerate(case["roots"][:12]):
                    spatial = mh.blocks(model)[0].spatial(r["frequency_hz"])
                    rate = max(math.sqrt(abs(row["k_squared_per_m2"])) for row in spatial)
                    scale = scales[k]*max(1., rate*sum(case["lengths"]))
                    dc = reduced.contraction_diagnostic(values[k, :, 3], model.section.nu*ux[k], weight, scale)
                    if dc["D_c"] is not None and not -1e-14 <= dc["D_c"] <= 1+1e-14:
                        raise ArithmeticError("Contraction triangle inequality failed")
                    dc.update({"beta_deg": beta, "position": k+1})
                    dcrows.append(dc)
                contraction.extend(dcrows)
        pair_records = {}
        for a, b in (("elementary", "mindlin_herrmann"), ("rayleigh_love_planar", "mindlin_herrmann"), ("elementary", "rayleigh_love_planar")):
            va, ta, weight = fields[a]
            vb, tb, _ = fields[b]
            od, dna, dnb = reduced.overlap(va, vb, np.repeat(weight, 2))
            ot, tna, tnb = reduced.overlap(ta, tb, weight, zeros[a], zeros[b])
            record = {"beta_deg": beta, "first": a, "second": b, "O_d": od, "O_theta": ot,
                "d_norm_first": dna, "d_norm_second": dnb, "theta_norm_first": tna, "theta_norm_second": tnb,
                "theta_numerical_zero_scale_first": zero_scales[a], "theta_numerical_zero_scale_second": zero_scales[b],
                "theta_norm_status_first": ["SMALL_NORM" if v else "INFORMATIVE" for v in zeros[a]],
                "theta_norm_status_second": ["SMALL_NORM" if v else "INFORMATIVE" for v in zeros[b]],
                "displacement_correspondence": ranking(od), "theta_correspondence": ranking(ot)}
            overlaps.append(record)
            pair_records[a, b] = record
            for direction, rows in record["displacement_correspondence"].items():
                non_diagonal.extend({"beta_deg": beta, "first": a, "second": b, "direction": direction, **row}
                    for row in rows if row["status"] == "POSITION_CORRESPONDENCE_NONDIAGONAL")
        beta_table = []
        for k in range(12):
            fe, fr, fm = (active[name]["roots"][k]["frequency_hz"] for name in config["models"])
            de, dr = fe/fm-1, fr/fm-1
            e = pair_records["elementary", "mindlin_herrmann"]["displacement_correspondence"]["rows"][k]
            r = pair_records["rayleigh_love_planar", "mindlin_herrmann"]["displacement_correspondence"]["rows"][k]
            row = {"beta_deg": beta, "position": k+1, "f_E": fe, "f_RL": fr, "f_MH": fm,
                "delta_E_signed": de, "abs_delta_E": abs(de), "delta_RL_signed": dr, "abs_delta_RL": abs(dr),
                "delta_RL_vs_E": fr/fe-1, "O_d_E_MH_diagonal": e["diagonal"],
                "O_d_RL_MH_diagonal": r["diagonal"], "best_match_E_to_MH": e["argmax"],
                "best_match_RL_to_MH": r["argmax"], "position_correspondence_E": e["status"],
                "position_correspondence_RL": r["status"], "D_c": dcrows[k]["D_c"], "D_c_status": dcrows[k]["status"]}
            beta_table.append(row)
        table.extend(beta_table)
        er, rr = (max(beta_table, key=lambda row: row[key]) for key in ("abs_delta_E", "abs_delta_RL"))
        defined = [d for d in dcrows if d["D_c"] is not None]
        maximum = max(defined, key=lambda d: d["D_c"]) if defined else None
        summary.append({"beta_deg": beta, "max_abs_delta_E": er["abs_delta_E"], "position_E": er["position"],
            "max_abs_delta_RL": rr["abs_delta_RL"], "position_RL": rr["position"],
            "median_abs_delta_E": float(np.median([r["abs_delta_E"] for r in beta_table])),
            "median_abs_delta_RL": float(np.median([r["abs_delta_RL"] for r in beta_table])),
            "non_diagonal_E_to_MH": sum(r["best_match_E_to_MH"] != r["position"] for r in beta_table),
            "non_diagonal_RL_to_MH": sum(r["best_match_RL_to_MH"] != r["position"] for r in beta_table),
            "max_D_c": maximum["D_c"] if maximum else None, "position_D_c": maximum["position"] if maximum else None})
    maxima = {name: max(table, key=lambda r: r[key]) for name, key in (("E_vs_MH", "abs_delta_E"), ("RL_vs_MH", "abs_delta_RL"))}
    return {"frequency_table": table, "overlaps": overlaps, "adjacent_gaps": gaps, "contraction": contraction,
        "summary": summary, "global_maxima": maxima, "non_diagonal_correspondence": non_diagonal}


def compute(config, base, model, length, frozen, reference, out):
    upper = 2*math.pi*config["frequency_ceiling_fstar"]
    catalogs = {l: tim_catalog(model, l, upper, base) for l in (.5, .35, .65)}
    gates, cases = comparator_gates(config, base, model, length, catalogs, upper)
    write_json(out/"comparator_gates.json", gates)
    statuses = {"ELEMENTARY_TIM_COUPLED_COMPARATOR": "PASS", "RAYLEIGH_LOVE_TIM_COUPLED_COMPARATOR": "PASS",
        "MHTIM_REFERENCE_REUSE": "PASS"}
    mh_count = general.count_function(model, length, reference["fixed_arm_poles"])
    failures = []
    for beta in config["beta_deg"]:
        for variant in config["models"]:
            try:
                if (variant, beta) not in cases:
                    if variant == "mindlin_herrmann":
                        case = mh_case(model, length, beta, reference, frozen, mh_count, upper)
                    else:
                        case = solve_reduced(model, (.5, .5), beta, variant, catalogs, upper)
                    cases[variant, beta] = case
                case = cases[variant, beta]
                case["ceiling_attempts"] = [{"upper_fstar": upper/(2*math.pi), "root_count": len(case["roots"]), "reason": "initial bounded inventory"}]
                if len(case["roots"]) < 13:
                    expanded = upper*config["one_ceiling_expansion_factor"]
                    if variant == "mindlin_herrmann":
                        raise ArithmeticError("MH guard missing; expanded MH pole coverage unavailable")
                    expanded_catalog = {.5: tim_catalog(model, .5, expanded, base)}
                    fresh = solve_reduced(model, (.5, .5), beta, variant, expanded_catalog, expanded)
                    fresh["ceiling_attempts"] = case["ceiling_attempts"]+[{"upper_fstar": expanded/(2*math.pi), "root_count": len(fresh["roots"]), "reason": "one allowed guard expansion"}]
                    case = cases[variant, beta] = fresh
                if len(case["roots"]) < 13:
                    raise ArithmeticError("INCOMPLETE_ROOT_INVENTORY after one ceiling expansion")
                case["reported_positions"] = 12
                case["guard_position"] = 13
                case["prefix_certificate"] = {"sorted_positions": list(range(1, 14)),
                    "positive_root_count": len(case["roots"]), "guard_omega": case["roots"][12]["omega"],
                    "lower_count": case["search"]["lower_count"], "upper_count": case["search"]["upper_count"],
                    "count_range_omega": case["search"]["range_omega"], "status": "PASS"}
                write_json(out/f"inventory_{variant}_{beta}.json", case)
                print(f"Inventory PASS: {variant}, beta={beta}, count={len(case['roots'])}, guard13={case['roots'][12]['frequency_hz']:.12g}", flush=True)
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                failures.append({"model": variant, "beta_deg": beta, "status": "UNRESOLVED", "reason": str(exc)})
                write_json(out/"unresolved_cases.json", failures)
    if failures:
        statuses.update({"COUPLED_HIERARCHY_ROOT_INVENTORY": "PARTIAL_PASS", "COUPLED_HIERARCHY_FIXED_BETA_SHAPE_OVERLAP": "PARTIAL", "COUPLED_HIERARCHY_SCREENING": "PARTIAL"})
        return {"statuses": statuses, "failures": failures}, cases
    data = diagnostics(config, model, cases)
    residual_maxima = {}
    for case in cases.values():
        for root in case["roots"]:
            d = root["diagnostics"]
            for key in ("clamp_scaled_residual", "equation_scaled_residual", "singular_ratio", "nonzero_singular_condition", "energy_relative_error"):
                residual_maxima[key] = max(residual_maxima.get(key, 0.), d[key])
            for key, value in d["joint_residual_scaled"].items():
                residual_maxima[key] = max(residual_maxima.get(key, 0.), value)
        residual_maxima["mass_gram_max_error"] = max(residual_maxima.get("mass_gram_max_error", 0.), case.get("mass_gram_max_error", 0.))
    statuses.update({"COUPLED_HIERARCHY_ROOT_INVENTORY": "PASS", "COUPLED_HIERARCHY_FIXED_BETA_SHAPE_OVERLAP": "COMPLETE", "COUPLED_HIERARCHY_SCREENING": "COMPLETE"})
    data.update({"statuses": statuses, "failures": [], "policy": POLICY, "residual_maxima": residual_maxima,
        "model_definitions": {"elementary": "H=J=0, N=EA U_x; no c/R",
            "rayleigh_love_planar": "H=0,J=nu^2 rho Iy,N=(EA-J omega^2) U_x; no c/R",
            "mindlin_herrmann": "unchanged project_jang_reduced_rectangular; c/R closure retained"},
        "tim_pole_catalogs": {str(l): c for l, c in catalogs.items()},
        "spectrum_semantics": "independently sorted positions at each fixed beta; no reassignment"})
    return data, cases


def profiles(model, cases):
    for (variant, beta), case in cases.items():
        for r in case["roots"][:13]:
            for arm, length in enumerate(case["lengths"]):
                x = np.linspace(0., length, 201)
                a = np.array(r["coefficients"])[arm]
                if variant == "mindlin_herrmann":
                    fields = joint.arm_state(model, length, r["omega"], a, x)
                    names = joint.STATE_ORDER
                else:
                    fields = reduced.state(model, length, r["omega"], a, x, variant)
                    names = reduced.STATE_ORDER
                for point, value in zip(x, fields):
                    yield {"model": variant, "beta_deg": beta, "position": r["sorted_position"], "arm": arm+1, "x_local": float(point),
                        **{name: float(value[names.index(name)]) if name in names else None for name in joint.STATE_ORDER}}


def plot_saved(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = json.loads((out/"summary.json").read_text(encoding="utf-8"))
    rows = data["summary"]
    fig, ax = plt.subplots(figsize=(6.5, 4))
    for name, label in (("E", "Elementary / MH"), ("RL", "Planar Love / MH")):
        ax.plot([r["beta_deg"] for r in rows], [100*r[f"max_abs_delta_{name}"] for r in rows], "o-", label=label)
    ax.set(xlabel="beta (deg), fixed grid", ylabel="max model difference, first 12 (%)")
    ax.legend(); ax.grid(alpha=.25); fig.tight_layout(); fig.savefig(out/"max_differences.png", dpi=180); plt.close(fig)
    table = data["frequency_table"]
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    for ax, name, title in zip(axes, ("E", "RL"), ("Elementary / MH", "Planar Love / MH")):
        matrix = np.array([r[f"abs_delta_{name}"]*100 for r in table]).reshape(8, 12)
        im = ax.imshow(matrix, origin="lower", aspect="auto", extent=(.5, 12.5, -.5, 7.5))
        ax.set(xlabel="sorted position", title=title, yticks=range(8), yticklabels=[r["beta_deg"] for r in rows])
        fig.colorbar(im, ax=ax, label="absolute model difference (%)")
    axes[0].set_ylabel("beta (deg)"); fig.tight_layout(); fig.savefig(out/"position_differences.png", dpi=180); plt.close(fig)
    matrix = np.array([r["D_c"] if r["D_c"] is not None else np.nan for r in table]).reshape(8, 12)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.imshow(matrix, origin="lower", aspect="auto", vmin=0, vmax=1, extent=(.5, 12.5, -.5, 7.5))
    ax.set(xlabel="MH sorted position", ylabel="beta (deg)", yticks=range(8), yticklabels=[r["beta_deg"] for r in rows])
    fig.colorbar(im, ax=ax, label="D_c (undefined fields masked)"); fig.tight_layout(); fig.savefig(out/"contraction_diagnostic.png", dpi=180); plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check-sources", action="store_true")
    parser.add_argument("--compute", action="store_true")
    parser.add_argument("--plot-only", type=Path, metavar="BUNDLE")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    if not (args.check_sources or args.compute or args.plot_only):
        parser.error("Select --check-sources, --compute or --plot-only BUNDLE")
    if args.plot_only:
        if args.compute:
            parser.error("--plot-only and --compute are separate workflows")
        old = json.loads((args.plot_only/"manifest.json").read_text(encoding="utf-8"))
        if any(sha(args.plot_only/p) != h for p, h in old["artifact_hashes"].items()):
            raise ValueError("Saved artifact hashes changed")
        plot_saved(args.plot_only)
        if any(sha(args.plot_only/p) != h for p, h in old["artifact_hashes"].items()):
            raise ValueError("Plot regeneration differs from accepted artifact bytes; retain current files for review")
        # Reproducible plot-only leaves identical images and manifest intact.
        print("Plots from saved tables; zero root evaluations:", args.plot_only)
        return 0 if old["statuses"]["COUPLED_HIERARCHY_SCREENING"] == "COMPLETE" else 1
    config, base, model, length, checked, frozen, reference = check_inputs()
    geometry_gates()
    if args.check_sources:
        print("Accepted inputs, source/reference hashes, common geometry and duality PASS; no roots")
    if not args.compute:
        return 0
    fingerprint, inputs = identity(config, checked)
    out = args.output_dir/fingerprint
    if (out/"manifest.json").exists():
        old = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
        if old["identity"] != inputs or any(sha(out/p) != h for p, h in old["artifact_hashes"].items()):
            raise ValueError("Stale screening cache")
        print("Verified screening cache; zero root evaluations:", out)
        return 0 if old["statuses"]["COUPLED_HIERARCHY_SCREENING"] == "COMPLETE" else 1
    out.mkdir(parents=True, exist_ok=True)
    write_json(out/"parameters.json", {"config": config, "policy": POLICY, "section": vars(model.section), "coefficients_MH": model.coefficients})
    try:
        data, cases = compute(config, base, model, length, frozen, reference, out)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
        write_json(out/"failure.json", {"status": "FAIL", "reason": str(exc), "identity": inputs})
        raise
    write_json(out/"summary.json", data)
    if data["statuses"]["COUPLED_HIERARCHY_SCREENING"] == "COMPLETE":
        write_csv(out/"frequencies.csv", data["frequency_table"])
        write_csv(out/"adjacent_gaps.csv", data["adjacent_gaps"])
        write_csv(out/"contraction.csv", data["contraction"])
        write_json(out/"overlaps.json", data["overlaps"])
        write_csv(out/"mode_profiles.csv", list(profiles(model, cases)))
        plot_saved(out)
    write_json(out/"manifest.json", {"identity": inputs,
        "artifact_hashes": {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"},
        "git_branch": subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True, encoding="utf-8"),
        "command": subprocess.list2cmdline([sys.executable, *sys.argv]), "statuses": data["statuses"]})
    write_json(args.output_dir/"current.json", {"fingerprint": fingerprint, "directory": str(out.resolve())})
    print(data["statuses"], out)
    return 0 if data["statuses"]["COUPLED_HIERARCHY_SCREENING"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
