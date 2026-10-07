"""Five-thickness/three-angle screening using unchanged hierarchy solvers.

New orchestration/output contract: coefficient scaling, thickness trends and
beta0 direct profiles. No new physics solver, classification or continuation.
Only sub-cutoff prefixes are needed; unsupported tails are never simulated.
"""
from __future__ import annotations

import argparse
import copy
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
import numpy as np
from scripts.analysis import screen_coupled_longitudinal_theory_hierarchy as hierarchy
from scripts.analysis import verify_mindlin_herrmann_timoshenko_general_beta_joint as general
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import coupled_longitudinal_comparators as reduced
from scripts.lib import mindlin_herrmann_longitudinal as mh
from scripts.lib import mindlin_herrmann_timoshenko_joint as joint
from scripts.lib.isotropic_rectangular_timoshenko_coupled_beams import rectangular_section

CONFIG = ROOT/"data/input/coupled_longitudinal_theory_thickness_screening.json"
OUTPUT = ROOT/"results/coupled_longitudinal_theory_thickness_screening"
VERSION = "bounded-thickness-prefix-direct-profile-audit-v1"
POLICY = hierarchy.POLICY


def make_model(config, h):
    p, g = config["material"], config["geometry"]
    section = rectangular_section(E=p["E"], rho=p["rho"], nu=p["nu"], K=p["kappa"], width=g["b"], thickness=h)
    return mh.project_jang_reduced_rectangular(section)


def coefficient_ratios(E, rho, nu, kappa, b, h):
    """Also exact Fraction arithmetic; same coefficient definitions as source."""
    A, I = b*h, b*h**3/12
    G = E/(2*(1+nu))
    m, C = rho*A, E*A/(1-nu**2)
    J, j, H = nu**2*rho*I, rho*I, kappa*G*I
    B, S, r = E*I, kappa*G*A, rho*I
    return {"EA_over_m": E*A/m, "J_over_m": J/m, "C_over_m": C/m,
        "j_over_m": j/m, "H_over_C": H/C, "H_over_j": H/j, "C_over_j": C/j,
        "B_over_m": B/m, "r_over_m": r/m, "S_over_m": S/m, "B_over_S": B/S,
        "S_over_r": S/r, "static_MH_layer_squared": H/(E*A)}


def scaling_audit(config):
    records = []
    for sh, h in zip(config["s_h"], config["h"]):
        model = make_model(config, h)
        s, p = model.section, model.coefficients
        arm = config["geometry"]["L"]*.5
        ratios = coefficient_ratios(s.E, s.rho, s.nu, s.K, s.width, s.thickness)
        doubled = coefficient_ratios(s.E, s.rho, s.nu, s.K, s.width*2, s.thickness)
        if not all(math.isclose(ratios[key], doubled[key], rel_tol=3e-15) for key in ratios):
            raise ArithmeticError("Width cancellation audit failed")
        # Same object, hence exact shared bending inputs for all three models.
        common = {key: p[key] for key in ("m", "B", "S", "r")}
        record = {"s_h": sh, "h": h, "A": s.area, "I": s.inertia, "I_over_A": s.inertia/s.area,
            "q_h": s.inertia/s.area/arm**2, "q_h_over_q_h0": sh**2,
            "h_over_L_arm": h/arm, "L_arm_over_h": arm/h, "b_over_h": s.width/h,
            "coefficient_ratios": ratios, "MH_coefficients": p, "RL_J": s.nu**2*s.rhoI,
            "bending_inputs_by_model": {name: common for name in config["models"]},
            "MH_cutoff_fstar": mh.blocks(model)[0].cutoff_hz,
            "Tim_cutoff_fstar": mh.blocks(model)[1].cutoff_hz}
        records.append(record)
    # Exact rational checks, without a new symbolic dependency.
    for sh in (F(1), F(5, 4), F(3, 2), F(7, 4), F(2)):
        b, h0 = F(1, 5), F(1, 20)
        h = sh*h0
        assert (b*h**3/12)/(b*h) == h**2/12
        assert (h**2/12)/(h0**2/12) == sh**2
        a = coefficient_ratios(F(1), F(1), F(3, 10), F(5, 6), b, h)
        assert a == coefficient_ratios(F(1), F(1), F(3, 10), F(5, 6), 3*b, h)
    return {"status": "PASS", "raw_h_powers": {"A,m,EA,C,S": 1, "I,J,j,H,B,r": 3},
        "normalized_scale": "I/A=h^2/12; common b cancels; cutoffs~h^-1, static MH layer~h",
        "not_a_spectral_power_law": True, "records": records}


def check_inputs():
    config = json.loads(CONFIG.read_text(encoding="utf-8"))
    baseline, base, model, length, checked, frozen, reference = hierarchy.check_inputs()
    if config["s_h"] != [1., 1.25, 1.5, 1.75, 2.] or config["h"] != [.05, .0625, .075, .0875, .1] or config["beta_deg"] != [0, 45, 90]:
        raise ValueError("Only the five predeclared thicknesses and three angles are authorized")
    if config["material"] != baseline["material"] or config["geometry"] != {"b": .2, "h0": .05, "L": 1., "split": .5}:
        raise ValueError("Material/width/length/kappa control changed")
    path = ROOT/config["baseline_bundle"]
    manifest = json.loads((path/"manifest.json").read_text(encoding="utf-8"))
    _, current = hierarchy.identity(baseline, checked)
    if manifest["identity"] != current or any(sha(path/name) != digest for name, digest in manifest["artifact_hashes"].items()):
        raise ValueError("Previous hierarchy bundle stale or changed")
    if manifest["statuses"]["COUPLED_HIERARCHY_SCREENING"] != "COMPLETE":
        raise ValueError("No accepted hierarchy baseline")
    scaling_audit(config)
    return config, base, checked, path


def identity(config, base, checked):
    paths = [CONFIG, Path(__file__), Path(hierarchy.__file__), Path(general.__file__),
        Path(reduced.__file__), Path(mh.__file__), Path(joint.__file__),
        ROOT/"scripts/lib/bishop_longitudinal.py", ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
        ROOT/"scripts/lib/reddy_inplane_geometry.py"]
    _, inherited = hierarchy.identity(json.loads((ROOT/config["baseline_config"]).read_text(encoding="utf-8")), checked)
    record = {"version": VERSION, "config": config, "policy": POLICY, "single_config": base,
        "inherited_identity": inherited, "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in paths},
        "baseline_manifest": sha(ROOT/config["baseline_bundle"]/"manifest.json")}
    return hashlib.sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()[:16], record


def catalogs(model, base, config, ceiling):
    """Same bounded basis/counts; no use of unsupported above-cutoff tail."""
    policy, arm = config["catalog_policy"], .5
    tim = mh.blocks(model)[1]
    k = tim.spatial(min(ceiling, policy["cutoff_margin"]*tim.cutoff_hz))[0]["wavenumber_per_m"]
    fraction = policy["tim_ss_fraction"]
    index = math.ceil(k*arm/math.pi-fraction)+fraction
    upper_tim = tim.temporal(index*math.pi/arm)[0]["frequency_hz"]
    if upper_tim >= policy["cutoff_margin"]*tim.cutoff_hz:
        klim = tim.spatial(policy["cutoff_margin"]*tim.cutoff_hz)[0]["wavenumber_per_m"]
        index = math.floor(klim*arm/math.pi-fraction)+fraction
        upper_tim = tim.temporal(index*math.pi/arm)[0]["frequency_hz"]
    # Keep MH Young lower-form optical count zero for a sharp finite catalog.
    eta = base["policy"]["young_eta"]
    bound_cut = mh.finite_count_upper_bound(model, arm, 0., "mh", eta)["lower_contraction_cutoff_hz"]
    limit = min(ceiling, policy["cutoff_margin"]*bound_cut)
    spacing = math.sqrt(model.section.E/model.section.rho)/(2*arm)
    frac = policy["mh_scalar_fraction"]
    upper_mh = (math.ceil(limit/spacing-frac)+frac)*spacing
    if upper_mh >= policy["cutoff_margin"]*bound_cut:
        upper_mh = (math.floor(policy["cutoff_margin"]*bound_cut/spacing-frac)+frac)*spacing
    result = {}
    for block, upper in (("mh", upper_mh), ("timoshenko", upper_tim)):
        attempts = []
        bound = mh.finite_count_upper_bound(model, arm, 2*math.pi*upper, block)
        for attempt in range(2):
            active = {**base["policy"], "scan_intervals": base["policy"]["scan_intervals"]*2**attempt}
            roots, search = mh.finite_roots(model, arm, block, .001*math.pi, 2*math.pi*upper, active)
            attempts.append(search)
            if len(roots) == bound["upper_count"]:
                break
        if len(roots) != bound["upper_count"] or search["failed_intervals"]:
            raise ArithmeticError(f"Pole catalog not count-certified: {block}, {upper}, found={len(roots)}, bound={bound}, attempts={attempts}")
        result[block] = {"roots": roots, "search": search, "count_bound": bound, "attempts": attempts}
    usable = min(ceiling, policy["query_coverage_margin"]*upper_tim, policy["query_coverage_margin"]*upper_mh)
    return result, {.5: {"roots": result["timoshenko"]["roots"], "poles": [r["omega"] for r in result["timoshenko"]["roots"]],
        "upper": 2*math.pi*upper_tim, "count_bound": result["timoshenko"]["count_bound"]}}, usable


def direct_spot(model, case):
    """Independent direct boundary systems + full physical beta0 profiles."""
    upper = case["roots"][12]["bracket_omega"][1]
    variant = case["model"]
    candidates, attempts = [], []
    if variant != "mindlin_herrmann":
        candidates.extend({"omega": w, "operator_block": "reduced_axial"}
            for w in reduced.axial_poles(model, 1., variant, upper))
    for block in (("mh", "timoshenko") if variant == "mindlin_herrmann" else ("timoshenko",)):
        roots, search = mh.finite_roots(model, 1., block, .001*math.pi, upper, POLICY)
        candidates.extend({**r, "operator_block": block} for r in roots)
        attempts.append({"block": block, **search})
    candidates.sort(key=lambda r: r["omega"])
    if len(candidates) != 13:
        raise ArithmeticError(f"Direct spot inventory mismatch: {variant}, {len(candidates)}")
    nodes, weights = np.polynomial.legendre.leggauss(POLICY["quadrature_order"])
    x, weight = (nodes+1)/4, weights/4
    rows = []
    for root, direct in zip(case["roots"][:13], candidates):
        difference = abs(root["omega"]/direct["omega"]-1)
        vectors, references = [], []
        if direct["operator_block"] != "reduced_axial":
            pure = mh.finite_mode(model, 1., direct["omega"], direct["operator_block"], POLICY["quadrature_order"])
        for arm in range(2):
            positions = x if arm == 0 else 1-x
            a = np.array(root["coefficients"])[arm]
            if variant == "mindlin_herrmann":
                value = joint.arm_state(model, .5, root["omega"], a, x)[:, :4]
                reference = np.zeros_like(value)
                ids = (0, 1) if direct["operator_block"] == "mh" else (2, 3)
                reference[:, ids] = mh.finite_state_basis(model, 1., direct["omega"], positions, direct["operator_block"])[:, :2]@pure["coefficients"]
                if arm == 1:
                    reference *= [-1., 1., -1., 1.]
                scale = [1., 1., 1., 1.]  # L=1, dimensionally scaled q
            else:
                value = reduced.state(model, .5, root["omega"], a, x, variant)[:, :3]
                reference = np.zeros_like(value)
                if direct["operator_block"] == "reduced_axial":
                    s = reduced.segment(model, 1., variant)
                    k = direct["omega"]*math.sqrt(s.m/(s.EA-s.J*direct["omega"]**2))
                    reference[:, 0] = math.sqrt(2/(s.m+s.J*k*k))*np.sin(k*positions)
                else:
                    reference[:, (1, 2)] = mh.finite_state_basis(model, 1., direct["omega"], positions, "timoshenko")[:, :2]@pure["coefficients"]
                if arm == 1:
                    reference *= [-1., -1., 1.]
                scale = [1., 1., 1.]
            vectors.append(value*np.array(scale)*np.sqrt(weight)[:, None])
            references.append(reference*np.array(scale)*np.sqrt(weight)[:, None])
        v, ref = np.concatenate(vectors).ravel(), np.concatenate(references).ravel()
        sign = 1. if v@ref >= 0 else -1.
        error = float(np.linalg.norm(v-sign*ref)/np.linalg.norm(ref))
        if difference > POLICY["beta0_frequency_relative_tol"] or error > POLICY["symmetry_component_L2_tol"]:
            raise ArithmeticError(f"Direct beta0 frequency/profile recovery failed: {variant}, {difference}, {error}")
        rows.append({"position": root["sorted_position"], "relative_frequency_difference": difference,
            "kinematic_profile_L2_relative": error, "joint_residuals": root["diagnostics"]["joint_residual_scaled"]})
    return {"model": variant, "status": "PASS", "upper_omega": upper, "direct_roots": candidates,
        "direct_search": attempts, "rows": rows, "count_qualification": "13 independent direct eigenpairs match the count-certified artificial-interface prefix; no cross-case tracking"}


def compute(config, base, baseline_path, out):
    scaling = scaling_audit(config)
    results, all_cases, spot_checks, failures = [], {}, [], []
    for sh, h in zip(config["s_h"], config["h"]):
        model = make_model(config, h)
        cases = {}
        catalog_record = None
        if sh != 1:
            try:
                catalog_record, tim_catalog, usable = catalogs(model, base, config, config["frequency_ceiling_fstar"])
                count = general.count_function(model, 1., catalog_record)
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                failures.append({"s_h": sh, "h": h, "status": "UNRESOLVED", "reason": str(exc)})
                continue
        for beta in config["beta_deg"]:
            for name in config["models"]:
                try:
                    if sh == 1:
                        case = json.loads((baseline_path/f"inventory_{name}_{beta}.json").read_text(encoding="utf-8"))
                        case["origin_thickness"] = "immutable hierarchy baseline"
                    else:
                        if name == "mindlin_herrmann":
                            case = general.solve_case(model, 1., beta, count, 2*math.pi*usable)
                            case["model"] = name
                        else:
                            case = hierarchy.solve_reduced(model, (.5, .5), beta, name, tim_catalog, 2*math.pi*usable)
                        case["ceiling_attempts"] = [{"configured_ceiling_fstar": config["frequency_ceiling_fstar"],
                            "usable_count_window_fstar": usable, "reason": "count-certified prefix inside unchanged sub-cutoff analytic domain"}]
                        if len(case["roots"]) < 13:
                            expanded = config["frequency_ceiling_fstar"]*config["one_ceiling_expansion_factor"]
                            cat2, tim2, use2 = catalogs(model, base, config, expanded)
                            if use2 <= usable:
                                raise ArithmeticError("INCOMPLETE_ROOT_INVENTORY: one expanded ceiling still limited by unchanged sub-cutoff domain")
                            fresh = general.solve_case(model, 1., beta, general.count_function(model, 1., cat2), 2*math.pi*use2) if name == "mindlin_herrmann" else hierarchy.solve_reduced(model, (.5, .5), beta, name, tim2, 2*math.pi*use2)
                            fresh["model"] = name
                            fresh["ceiling_attempts"] = case["ceiling_attempts"]+[{"configured_ceiling_fstar": expanded, "usable_count_window_fstar": use2, "reason": "one allowed guard expansion"}]
                            case = fresh
                    if len(case["roots"]) < 13:
                        raise ArithmeticError("INCOMPLETE_ROOT_INVENTORY after bounded attempts")
                    for r in case["roots"]:
                        reduced.validate_mode(r, POLICY)
                    case.update({"s_h": sh, "h": h, "guard_position": 13, "reported_positions": 12,
                        "prefix_certificate": {"status": "PASS", "lower_count": case["search"]["lower_count"],
                            "upper_count": case["search"]["upper_count"], "guard_omega": case["roots"][12]["omega"],
                            "certified_range_omega": case["search"]["range_omega"]}})
                    cases[name, beta] = case
                    all_cases[sh, name, beta] = case
                    write_json(out/f"inventory_{sh:g}_{name}_{beta}.json", case)
                    print(f"Inventory PASS h={h:g}, beta={beta}, {name}, count={len(case['roots'])}", flush=True)
                except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                    failures.append({"s_h": sh, "h": h, "beta_deg": beta, "model": name, "status": "UNRESOLVED", "reason": str(exc)})
        if len(cases) != 9:
            write_json(out/"unresolved_cases.json", failures)
            continue
        spots = [direct_spot(model, cases[name, 0]) for name in config["models"]]
        spot_checks.extend({"s_h": sh, "h": h, **spot} for spot in spots)
        data = hierarchy.diagnostics(config, model, cases)
        for key in ("frequency_table", "overlaps", "adjacent_gaps", "contraction", "non_diagonal_correspondence", "summary"):
            for row in data[key]:
                row.update({"s_h": sh, "h": h})
        for row in data["summary"]:
            beta = row["beta_deg"]
            primary = [r for r in data["frequency_table"] if r["beta_deg"] == beta]
            row["min_diagonal_O_d_E_MH"] = min(r["O_d_E_MH_diagonal"] for r in primary)
            row["min_diagonal_O_d_RL_MH"] = min(r["O_d_RL_MH_diagonal"] for r in primary)
            row["min_gap_by_model"] = {name: min((r for r in data["adjacent_gaps"] if r["beta_deg"] == beta and r["model"] == name), key=lambda r: r["relative_adjacent_gap"]) for name in config["models"]}
        data["pole_catalogs"] = catalog_record
        results.append(data)
        write_json(out/f"diagnostics_{sh:g}.json", data)
        write_json(out/"direct_beta0_checks.json", spot_checks)
        profile_rows = list(hierarchy.profiles(model, cases))
        write_csv(out/f"mode_profiles_{sh:g}.csv", [{"s_h": sh, "h": h, **r} for r in profile_rows])
    final = {key: [r for data in results for r in data[key]] for key in ("frequency_table", "overlaps", "adjacent_gaps", "contraction", "non_diagonal_correspondence", "summary")}
    for row in final["summary"]:
        baseline = next(r for r in final["summary"] if r["s_h"] == 1 and r["beta_deg"] == row["beta_deg"])
        row.update({"max_delta_E_over_baseline": row["max_abs_delta_E"]/baseline["max_abs_delta_E"],
            "max_delta_RL_over_baseline": row["max_abs_delta_RL"]/baseline["max_abs_delta_RL"], "q_h_over_q_h0": row["s_h"]**2})
    final.update({"scaling_audit": scaling, "failures": failures,
        "global_maxima": {name: max(final["frequency_table"], key=lambda r: r[key]) for name, key in (("E_vs_MH", "abs_delta_E"), ("RL_vs_MH", "abs_delta_RL"))} if results else {},
        "statuses": {"COUPLED_THICKNESS_SCALING_AUDIT": "PASS", "COUPLED_THICKNESS_ROOT_INVENTORY": "PASS" if not failures else "PARTIAL_PASS",
            "COUPLED_THICKNESS_FIXED_CASE_OVERLAP": "COMPLETE" if not failures else "PARTIAL", "COUPLED_THICKNESS_SCREENING": "COMPLETE" if not failures else "PARTIAL"}})
    quality = {}
    for case in all_cases.values():
        for r in case["roots"]:
            d = r["diagnostics"]
            for key in ("clamp_scaled_residual", "equation_scaled_residual", "singular_ratio", "nonzero_singular_condition", "energy_relative_error"):
                quality[key] = max(quality.get(key, 0), d[key])
            for key, value in d["joint_residual_scaled"].items():
                quality[key] = max(quality.get(key, 0), value)
        quality["mass_gram_max_error"] = max(quality.get("mass_gram_max_error", 0), case.get("mass_gram_max_error", 0))
    final["residual_maxima"] = quality
    return final


def plot_saved(out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    data = json.loads((out/"summary.json").read_text(encoding="utf-8"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
    for ax, name in zip(axes, ("E", "RL")):
        for beta in (0, 45, 90):
            rows = [r for r in data["summary"] if r["beta_deg"] == beta]
            ax.plot([r["s_h"] for r in rows], [100*r[f"max_abs_delta_{name}"] for r in rows], "o-", label=f"beta={beta} deg")
        ax.set(xlabel="s_h = h / h0", ylabel=f"max |delta {name} / MH| (%)")
        ax.grid(alpha=.25); ax.legend()
    fig.tight_layout(); fig.savefig(out/"thickness_max_differences.png", dpi=180); plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), sharey=True)
    for ax, name in zip(axes, ("E", "RL")):
        matrix = np.array([100*r[f"abs_delta_{name}"] for r in data["frequency_table"] if r["beta_deg"] == 45]).reshape(5, 12)
        im = ax.imshow(matrix, origin="lower", aspect="auto", extent=(.5, 12.5, .875, 2.125))
        ax.set(xlabel="sorted position at beta45", ylabel="s_h", title=f"{name} / MH")
        fig.colorbar(im, ax=ax, label="model difference (%)")
    fig.tight_layout(); fig.savefig(out/"thickness_position_differences.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4))
    for beta in (0, 45, 90):
        rows = [r for r in data["summary"] if r["beta_deg"] == beta]
        ax.plot([r["s_h"] for r in rows], [r["max_D_c"] for r in rows], "o-", label=f"beta={beta} deg")
    ax.set(xlabel="s_h", ylabel="max defined D_c, first12")
    ax.grid(alpha=.25); ax.legend(); fig.tight_layout(); fig.savefig(out/"thickness_contraction_diagnostic.png", dpi=180); plt.close(fig)


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
            parser.error("plot-only and compute are separate")
        manifest = json.loads((args.plot_only/"manifest.json").read_text(encoding="utf-8"))
        if any(sha(args.plot_only/name) != digest for name, digest in manifest["artifact_hashes"].items()):
            raise ValueError("Saved thickness artifacts changed")
        plot_saved(args.plot_only)
        if any(sha(args.plot_only/name) != digest for name, digest in manifest["artifact_hashes"].items()):
            raise ValueError("Plot-only bytes differ")
        print("Plot-only, zero root evaluations:", args.plot_only)
        return 0
    config, base, checked, baseline = check_inputs()
    fingerprint, inputs = identity(config, base, checked)
    if args.check_sources:
        print("Sources/baseline hashes and exact thickness scaling PASS; zero roots")
    if not args.compute:
        return 0
    out = args.output_dir/fingerprint
    if (out/"manifest.json").exists():
        manifest = json.loads((out/"manifest.json").read_text(encoding="utf-8"))
        if manifest["identity"] != inputs or any(sha(out/name) != digest for name, digest in manifest["artifact_hashes"].items()):
            raise ValueError("Changed/stale thickness cache")
        print("Verified thickness cache, zero roots:", out)
        return 0 if manifest["statuses"]["COUPLED_THICKNESS_SCREENING"] == "COMPLETE" else 1
    out.mkdir(parents=True, exist_ok=True)
    write_json(out/"parameters.json", {"config": config, "policy": POLICY})
    write_json(out/"scaling_audit.json", scaling_audit(config))
    try:
        data = compute(config, base, baseline, out)
    except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
        write_json(out/"failure.json", {"status": "FAIL", "reason": str(exc), "identity": inputs})
        raise
    write_json(out/"summary.json", data)
    for key in ("frequency_table", "adjacent_gaps", "contraction"):
        if data[key]:
            write_csv(out/("frequencies.csv" if key == "frequency_table" else key+".csv"), data[key])
    write_json(out/"overlaps.json", data["overlaps"])
    if data["statuses"]["COUPLED_THICKNESS_SCREENING"] == "COMPLETE":
        plot_saved(out)
    write_json(out/"manifest.json", {"identity": inputs, "statuses": data["statuses"],
        "artifact_hashes": {p.name: sha(p) for p in out.iterdir() if p.is_file() and p.name != "manifest.json"},
        "git_branch": subprocess.check_output(["git", "branch", "--show-current"], text=True).strip(),
        "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        "git_status": subprocess.check_output(["git", "status", "--short"], text=True, encoding="utf-8"),
        "command": subprocess.list2cmdline([sys.executable, *sys.argv])})
    write_json(args.output_dir/"current.json", {"directory": str(out.resolve()), "fingerprint": fingerprint})
    print(data["statuses"], out)
    return 0 if data["statuses"]["COUPLED_THICKNESS_SCREENING"] == "COMPLETE" else 1


if __name__ == "__main__":
    raise SystemExit(main())
