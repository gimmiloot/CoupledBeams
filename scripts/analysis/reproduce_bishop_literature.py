"""One bounded literature workflow, separate source/compute/plot/precision actions.

Diagnostic-only. New fourth-order longitudinal/source-precision contract cannot
be a preset of the project's angled bending scripts. Reuses bishop_longitudinal.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
FIXTURE = ROOT / "data/input/bishop_literature_sources.json"
MODULE = ROOT / "scripts/lib/bishop_longitudinal.py"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temp, path)


def write_csv(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def source_check():
    fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
    checked = {}
    for name, source in fixture["sources"].items():
        actual = sha(ROOT / source["path"])
        if actual != source["sha256"]:
            raise ValueError("Source changed; repeat transcription audit: " + source["path"])
        checked[name] = {"path": source["path"], "sha256": actual, "pages": source["pages"]}
    if len(fixture["sources"]["popov"]["ratio_strings"]) != 30:
        raise ValueError("Expected 30 source ratios")
    return fixture, checked


def provenance(checked):
    from importlib.metadata import version
    def git(*args):
        result = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
        return result.stdout.strip() if result.returncode == 0 else "unavailable"
    identity = {
        "schema": "bishop-literature-results-v1", "sources": checked,
        "files": {str(p.relative_to(ROOT)).replace("\\", "/"): sha(p) for p in (FIXTURE, MODULE, Path(__file__))},
        "python": platform.python_version(),
        "versions": {name: version(name) for name in ("numpy", "scipy", "matplotlib", "mpmath")},
        "head": git("rev-parse", "HEAD"),
    }
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    return {"fingerprint": fingerprint, "identity": identity, "git_branch": git("branch", "--show-current"),
            "git_status": git("status", "--short"), "executable": sys.executable,
            "cwd": str(ROOT), "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
            "spectrum_semantics": "sorted_positions", "scope": "diagnostic-only fixed literature rods"}


def print_comparison(value, printed, digits=4):
    import math
    p = float(printed)
    step = 10.**(math.floor(math.log10(p))-digits+1)
    nearest = p-step/2 <= value < p+step/2
    truncation = p <= value < p+step
    return {"printed_hz": printed, "calculated_hz": value, "signed_difference_hz": value-p,
            "absolute_difference_hz": abs(value-p), "relative_difference": (value-p)/p,
            "absolute_relative_difference": abs(value-p)/p, "source_significant_digits": digits,
            "print_step_hz": step, "nearest_interval_hz": [p-step/2, p+step/2],
            "truncation_interval_hz": [p, p+step],
            "nearest_status": "PRINT_MATCH" if nearest else "PRINT_MISMATCH",
            "truncation_status": "PRINT_MATCH" if truncation else "PRINT_MISMATCH"}


def numerical_status(v, contract):
    passed = (v["mass_orthogonality_max"] <= contract["mass_orthogonality_tol"] and
              all(m["physical_bc_max"] <= contract["bc_residual_tol"] and
                  m["ode_relative"] <= contract["ode_residual_tol"] and
                  m["energy_relative"] <= contract["energy_relative_tol"] and
                  m["nonnull_condition"] <= contract["nonnull_condition_max"] for m in v["modes"]))
    return "PASS" if passed else "UNRESOLVED"


def solve(segments, ends, count, fixture, label):
    from scripts.lib import bishop_longitudinal as b
    contract = fixture["numerical_contract"]
    case = "marais" if label == "marais" else "popov"
    roots, search = b.bounded_roots(lambda f: b.characteristic(segments, f, ends),
                                  contract[case+"_search_hz"], contract[case+"_scan_intervals"], count, contract)
    if search["status"] != "PASS":
        return roots, search, None, None
    used = roots if label == "marais" else roots[:30]
    checks, profiles = b.verify_modes(segments, used, ends)
    checks["status"] = numerical_status(checks, contract)
    return roots, search, checks, profiles


def compute_marais(fixture, directory):
    import numpy as np
    from scripts.lib import bishop_longitudinal as b
    source = fixture["sources"]["marais"]
    p = source["parameters_si"]
    segments = [b.circular_segment(L, p["E"], p["rho"], p["nu"], r)
                for L, r in zip(p["lengths"], p["radii"])]
    roots, search, checks, profiles = solve(segments, ("C", "F"), 5, fixture, "marais")
    write_json(directory/"search.json", search)
    if checks is None:
        return {"NUMERICAL_VERIFICATION": "UNRESOLVED"}
    contours = [b.marais_argument_count(segments, horizontal_samples=n) for n in (2048, 4096)]
    # U(0)=0 implies integral U² <= (2L/pi)² integral U'²; discard positive H energy.
    length = sum(s.L for s in segments)
    lower = np.sqrt(min(s.EA for s in segments)/(max(s.m for s in segments)*(2*length/np.pi)**2
                                                    + max(s.J for s in segments)))/(2*np.pi)
    complete = all(c["count"] == 5 and c["max_phase_step_rad"] < np.pi/2 and
                   c["branch_points_outside_contour"] for c in contours) and lower > 1
    checks["completeness"] = {"status": "PASS" if complete else "UNRESOLVED", "contours": contours,
                              "lower_bound_hz": float(lower), "claim": "finite numerical count in [1,31000] Hz; no root below 1 by energy bound"}
    rows = [{"mode": n, **print_comparison(float(f), printed, source["significant_digits"])}
            for n, (f, printed) in enumerate(zip(roots, source["frequency_strings_hz"]), 1)]
    write_json(directory/"frequencies.json", rows)
    write_csv(directory/"frequencies.csv", rows)
    write_json(directory/"verification.json", checks)
    write_csv(directory/"profiles.csv", profiles)
    return {"SOURCE_TRANSCRIPTION": "VERIFIED_LOCAL_PDF", "EQUATION_AND_BC_CONSISTENCY": "PASS_WITH_DOCUMENTED_RECONSTRUCTIONS",
            "NUMERICAL_VERIFICATION": checks["status"] if complete else "UNRESOLVED",
            "INDEPENDENT_VERIFICATION": "PENDING_HIGH_PRECISION_ACTION",
            "SOURCE_PRINT_MATCH": {"nearest": [r["nearest_status"] for r in rows], "truncation": [r["truncation_status"] for r in rows]},
            "MODE_SHAPE_COMPARISON": "VISUAL_ONLY_SEE_TRACKED_REPORT", "EXPERIMENT_REPRODUCTION": "NOT_APPLICABLE"}


def fit_rayleigh(p, experimental):
    import numpy as np
    from scipy.optimize import minimize_scalar
    n = np.arange(1, 31)
    wave = n*p["c_printed"]/(2*p["L"])
    a = np.pi**2*n*n*p["d"]**2/(8*p["L"]**2)
    calls = 0
    def objective(nu):
        nonlocal calls
        calls += 1
        return float(np.mean(np.abs(wave/np.sqrt(1+a*nu**2)-experimental)))
    lo, hi = p["nu_interval"]
    crossings = np.sqrt(np.maximum(0., ((wave/experimental)**2-1)/a))
    breaks = sorted(set([lo, hi, *[float(v) for v in crossings if lo < v < hi]]))
    candidates = [(v, objective(v)) for v in breaks]
    attempts = []
    # The L1 objective has kinks; include all exact residual-zero candidates,
    # and bounded minimization in each smooth interval. No Bishop fit.
    for left, right in zip(breaks[:-1], breaks[1:]):
        result = minimize_scalar(objective, bounds=(left, right), method="bounded",
                                 options={"xatol": 1e-12, "maxiter": 100})
        attempts.append({"interval": [left, right], "success": bool(result.success), "nfev": result.nfev})
        candidates.append((float(result.x), float(result.fun)))
    nu, error = min(candidates, key=lambda pair: pair[1])
    return {"nu": nu, "mae_hz": error, "interval": [lo, hi], "evaluations": calls,
            "method": "piecewise bounded MAE minimization including every zero-residual kink",
            "attempts": attempts, "source_nu_string": "0.337", "rounded_3_decimals": f"{nu:.3f}",
            "status": "PASS" if all(r["success"] for r in attempts) else "UNRESOLVED"}


def compute_popov(fixture, directory):
    import numpy as np
    from scripts.lib import bishop_longitudinal as b
    source, contract = fixture["sources"]["popov"], fixture["numerical_contract"]
    p = source["parameters_si"]
    experimental = np.array([float(r)*p["f1_exp"] for r in source["ratio_strings"]])
    fit = fit_rayleigh(p, experimental)
    write_json(directory/"fit.json", fit)
    rows, statistics, all_searches, all_checks = [], [], {}, {}
    nus = [*p["nu_fig5"], fit["nu"]]
    for nu in nus:
        s = b.speed_segment(p["L"], p["c_printed"], nu, p["d"])
        spectra = {}
        for end in ("F", "C", "UP"):
            roots, search, checks, profiles = solve([s], (end, end), contract["popov_expected_roots_in_range"], fixture, "popov")
            key = f"nu={nu:.12g}_{end}"
            all_searches[key] = search
            if checks is None:
                write_json(directory/"search.json", all_searches)
                return {"NUMERICAL_VERIFICATION": "UNRESOLVED", "case": key}
            spectra[end] = roots[:30]
            checks["guard_frequency_hz"] = float(roots[30])
            if end in ("F", "C"):
                independent, eq_search = b.bounded_roots(lambda f: b.popov_characteristic(s, f, end),
                    contract["popov_search_hz"], contract["popov_scan_intervals"], 31, contract)
                all_searches[key+"_printed_eq13"] = eq_search
                checks["eq13_relative_difference"] = float(max(abs(independent-roots)/roots)) if len(independent) == len(roots) else None
                if (eq_search["status"] != "PASS" or checks["eq13_relative_difference"] is None or
                        checks["eq13_relative_difference"] > contract["independent_frequency_relative_tol"]):
                    checks["status"] = "UNRESOLVED"
            else:
                checks["eq15_relative_difference"] = float(max(abs(b.explicit_frequencies(s, 30)-roots[:30])/roots[:30]))
                if checks["eq15_relative_difference"] > contract["independent_frequency_relative_tol"]:
                    checks["status"] = "UNRESOLVED"
            all_checks[key] = checks
            # F-F profiles at the two Fig.5 parameters; C/UP coefficients retained in JSON.
            if end == "F" and nu in p["nu_fig5"]:
                write_csv(directory/f"profiles_nu_{nu:.2f}.csv", profiles)
        wave = b.explicit_frequencies(b.speed_segment(p["L"], p["c_printed"], nu, p["d"], "wave"), 30)
        rayleigh = b.explicit_frequencies(b.speed_segment(p["L"], p["c_printed"], nu, p["d"], "rayleigh"), 30)
        for n in range(30):
            row = {"nu": nu, "mode": n+1, "ratio_printed": source["ratio_strings"][n],
                   "ratio_decimal_places": 3, "frequency_reconstructed_hz": float(experimental[n]),
                   "ratio_rounding_halfwidth_hz_if_nearest_fixed_f1": 0. if n == 0 else p["f1_exp"]*.0005,
                   "f_separately_printed_hz": "1297.812" if n == 0 else ("38557.99" if n == 29 else ""),
                   "wave_hz": float(wave[n]), "rayleigh_hz": float(rayleigh[n]),
                   "bishop_free_hz": float(spectra["F"][n]), "bishop_up_hz": float(spectra["UP"][n]),
                   "bishop_clamped_hz": float(spectra["C"][n])}
            for theory in ("wave", "rayleigh", "bishop_free"):
                delta = row[theory+"_hz"]-experimental[n]
                row[theory+"_signed_delta_hz"] = float(delta)
                row[theory+"_absolute_error_hz"] = float(abs(delta))
            rows.append(row)
        stats = {"nu": nu, "mae_hz": {"wave": float(np.mean(abs(wave-experimental))),
                 "rayleigh": float(np.mean(abs(rayleigh-experimental))),
                 "bishop_free": float(np.mean(abs(spectra["F"]-experimental)))},
                 "free_vs_up_max_hz": float(max(abs(spectra["F"]-spectra["UP"]))),
                 "free_vs_up_max_relative": float(max(abs(spectra["F"]-spectra["UP"])/spectra["UP"])),
                 "clamped_vs_up_min_hz": float(min(abs(spectra["C"]-spectra["UP"]))),
                 "clamped_vs_up_max_hz": float(max(abs(spectra["C"]-spectra["UP"]))),
                 "clamped_vs_up_max_relative": float(max(abs(spectra["C"]-spectra["UP"])/spectra["UP"]))}
        stats["mae_ranking"] = sorted(stats["mae_hz"], key=stats["mae_hz"].get)
        statistics.append(stats)
    write_csv(directory/"comparison.csv", rows)
    write_json(directory/"comparison.json", rows)
    write_json(directory/"summary.json", {"statistics": statistics, "fit": fit,
        "c_printed_m_s": p["c_printed"], "c_from_2Lf1_m_s": 2*p["L"]*p["f1_exp"],
        "f30_reconstructed_hz": float(experimental[-1]), "f30_separately_printed_hz": "38557.99",
        "rigid_mode": {"frequency_hz": 0, "shape": "constant", "positive_mode_index": None},
        "caution": "c derived from experimental f1; nu fitted to same rounded table; no independent material validation"})
    write_json(directory/"search.json", all_searches)
    write_json(directory/"verification.json", all_checks)
    passed = all(v["status"] == "PASS" for v in all_checks.values()) and fit["status"] == "PASS"
    return {"SOURCE_TRANSCRIPTION": "VERIFIED_LOCAL_PDF", "EQUATION_AND_BC_CONSISTENCY": "PASS_EQUATIONS_9_TO_15_WITH_REFERENCE_WARNINGS",
            "NUMERICAL_VERIFICATION": "PASS" if passed else "UNRESOLVED",
            "SOURCE_PRINT_MATCH": {"nu_rounded": fit["rounded_3_decimals"], "nu_printed": "0.337",
                                   "raw_measurements": "UNAVAILABLE", "figure5": "VISUAL_DIFFERENCES_SEE_REPORT"},
            "MODE_SHAPE_COMPARISON": "NOT_PUBLISHED_FOR_THIS_COMPARISON",
            "EXPERIMENT_REPRODUCTION": "ROUNDED_TABLE_REPRODUCTION_SEE_REPORT"}


def compute_uniform(fixture, directory):
    import numpy as np
    from scripts.lib import bishop_longitudinal as b
    p, contract = fixture["sources"]["popov"]["parameters_si"], fixture["numerical_contract"]
    rows, checks = [], {}
    for model in ("wave", "rayleigh", "bishop"):
        s = b.speed_segment(p["L"], p["c_printed"], .34, p["d"], model)
        end = "UP" if model == "bishop" else "U"
        exact = b.explicit_frequencies(s, 5)
        roots, search = b.bounded_roots(lambda f: b.characteristic([s], f, (end, end)),
            [1., 7000.], 140, 5, contract)
        if search["status"] != "PASS":
            checks[model] = search
            continue
        v, _ = b.verify_modes([s], roots, (end, end))
        parts = [b.Segment(s.L*.4, s.EA, s.m, s.H, s.J), b.Segment(s.L*.6, s.EA, s.m, s.H, s.J)]
        split, split_search = b.bounded_roots(lambda f: b.characteristic(parts, f, (end, end)),
                                            [1., 7000.], 140, 5, contract)
        if len(split) != 5:
            checks[model] = {"status": "UNRESOLVED", "search": search, "split_search": split_search}
            continue
        v["formula_max_relative"] = float(max(abs(roots-exact)/exact))
        v["split_max_relative"] = float(max(abs(split-roots)/roots))
        # Analytical constant represented exactly at omega=0; independently apply BC matrix.
        rigid = np.array([1., 0., 0., 0.] if s.H else [1., 0.])
        zero_matrix, _ = b.boundary_matrix([s], 0., ("F", "F"), balanced=False)
        v["rigid_bc_max"] = float(max(abs(zero_matrix@rigid)))
        v["rigid_energy"] = 0.0
        v["status"] = numerical_status(v, contract)
        if max(v["formula_max_relative"], v["split_max_relative"], v["rigid_bc_max"]) > contract["independent_frequency_relative_tol"]:
            v["status"] = "UNRESOLVED"
        v.update(search=search, split_search=split_search)
        checks[model] = v
        rows.extend({"model": model, "mode": n+1, "formula_hz": float(exact[n]),
                     "boundary_hz": float(roots[n]), "split_boundary_hz": float(split[n])} for n in range(5))
    write_json(directory/"verification.json", checks)
    if rows:
        write_csv(directory/"frequencies.csv", rows)
    return {"NUMERICAL_VERIFICATION": "PASS" if all(v["status"] == "PASS" for v in checks.values()) else "UNRESOLVED",
            "EQUATION_AND_BC_CONSISTENCY": "SEPARATE_SECOND_AND_FOURTH_ORDER_CONTROLS"}


def render(case, directory):
    # Plot-only reads completed data and never imports the solver module.
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    if case == "marais":
        with (directory/"profiles.csv").open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        fig, ax = plt.subplots(figsize=(8, 5))
        for n, style in zip(range(1, 6), ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]):
            selected = [r for r in rows if int(r["mode"]) == n]
            ax.plot([float(r["x_m"]) for r in selected], [float(r["Y"]) for r in selected],
                    linestyle=style, label=f"n={n}")
        ax.set(xlabel="x, m", ylabel="Y = U / max|U|", title="Marais: calculated first five modes (Fig. 2 comparison)")
        ax.grid(alpha=.3); ax.legend()
        fig.tight_layout()
        for extension in ("png", "pdf"):
            fig.savefig(directory/("marais_fig2."+extension), dpi=180)
        plt.close(fig)
    elif case == "popov":
        import numpy as np
        rows = json.loads((directory/"comparison.json").read_text(encoding="utf-8"))
        fig, axes = plt.subplots(2, 1, figsize=(11, 8), sharex=True)
        for ax, nu in zip(axes, [.31, .34]):
            selected = [r for r in rows if r["nu"] == nu]
            n = np.arange(1, 31)
            for offset, model in zip([-.25, 0, .25], ["wave", "rayleigh", "bishop_free"]):
                ax.bar(n+offset, [r[model+"_signed_delta_hz"] for r in selected], width=.25, label=model)
            ax.axhline(0, color="black", linewidth=.5)
            ax.set(ylabel="f_cal - f_reconstructed, Hz", title=f"Popov-Sadovsky Fig. 5 comparison; nu={nu}")
            ax.legend(); ax.grid(axis="y", alpha=.2)
        axes[-1].set(xlabel="positive sorted mode n", xticks=np.arange(1, 31))
        fig.tight_layout()
        for extension in ("png", "pdf"):
            fig.savefig(directory/("popov_fig5."+extension), dpi=180)
        plt.close(fig)


def validate_bundle(directory, fingerprint):
    manifest = json.loads((directory/"manifest.json").read_text(encoding="utf-8"))
    if manifest["fingerprint"] != fingerprint:
        raise ValueError("Stale provenance: compute this configuration first")
    for name, digest in manifest["artifacts"].items():
        if sha(directory/name) != digest:
            raise ValueError("Saved artifact changed: " + name)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=("all", "uniform", "marais", "popov"), default="all")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check-sources", action="store_true")
    action.add_argument("--compute", action="store_true")
    action.add_argument("--plot-only", action="store_true")
    action.add_argument("--high-precision", action="store_true")
    parser.add_argument("--output-dir", type=Path, default=ROOT/"results/bishop_literature")
    args = parser.parse_args(argv)
    fixture, checked = source_check()
    prov = provenance(checked)
    destination = args.output_dir/prov["fingerprint"][:16]
    if args.check_sources:
        write_json(args.output_dir/"source_check.json", {**prov, "status": "VERIFIED_LOCAL_PDF_HASHES", "fixtures": fixture})
        print("SOURCE_TRANSCRIPTION: verified PDF hashes and fixture structure; manual audit recorded in tracked report")
        return 0
    cases = ("uniform", "marais", "popov") if args.case == "all" else (args.case,)
    if args.high_precision:
        if args.case not in ("all", "marais"):
            parser.error("High precision is scoped to the five Marais targets")
        cases = ("marais",)
    failed = False
    for case in cases:
        directory = destination/case
        if args.compute and (directory/"manifest.json").exists():
            manifest = validate_bundle(directory, prov["fingerprint"])
            print(f"{case}: reused verified bundle (zero root calls): {directory}")
            failed |= manifest["statuses"].get("NUMERICAL_VERIFICATION") == "UNRESOLVED"
            continue
        start = time.perf_counter()
        if args.compute:
            directory.mkdir(parents=True, exist_ok=True)
            status = {"uniform": compute_uniform, "marais": compute_marais, "popov": compute_popov}[case](fixture, directory)
            manifest = {**prov, "case": case, "parameters_and_source_precision": fixture,
                        "statuses": status, "compute_seconds": time.perf_counter()-start, "actions": []}
        else:
            manifest = validate_bundle(directory, prov["fingerprint"])
            if args.plot_only:
                render(case, directory)
            elif args.high_precision:
                from scripts.lib import bishop_longitudinal as b
                import mpmath as mp
                if (directory/"independent.json").exists():
                    print(f"{case}: reused high-precision results: {directory}")
                    failed |= manifest["statuses"]["INDEPENDENT_VERIFICATION"] != "PASS"
                    continue
                seeds = [r["calculated_hz"] for r in json.loads((directory/"frequencies.json").read_text(encoding="utf-8"))]
                runs = [b.independent_marais_transfer(fixture["sources"]["marais"]["parameters_si"], seeds, dps)
                        for dps in fixture["numerical_contract"]["high_precision_dps"]]
                stability = []
                for a, z in zip(runs[0]["roots"], runs[1]["roots"]):
                    if a["status"] == z["status"] == "PASS":
                        with mp.workdps(75):
                            stability.append(float(abs(mp.mpf(a["frequency_hz"])-mp.mpf(z["frequency_hz"]))/mp.mpf(z["frequency_hz"])))
                contract = fixture["numerical_contract"]
                passed = (len(stability) == 5 and max(stability) <= contract["high_precision_stability_relative_tol"]
                          and all(r["relative_to_double"] <= contract["independent_frequency_relative_tol"]
                                  for run in runs for r in run["roots"] if r["status"] == "PASS"))
                write_json(directory/"independent.json", {"runs": runs, "precision_stability_relative": stability,
                           "status": "PASS" if passed else "UNRESOLVED"})
                manifest["statuses"]["INDEPENDENT_VERIFICATION"] = "PASS" if passed else "UNRESOLVED"
                failed |= not passed
        manifest["actions"].append({"command": prov["command"], "seconds": time.perf_counter()-start})
        manifest["artifacts"] = {p.name: sha(p) for p in sorted(directory.iterdir()) if p.is_file() and p.name != "manifest.json"}
        write_json(directory/"manifest.json", manifest)
        failed |= manifest["statuses"].get("NUMERICAL_VERIFICATION") == "UNRESOLVED"
        print(f"{case}: {manifest['statuses']}; {directory}")
    write_json(args.output_dir/"current.json", {"fingerprint": prov["fingerprint"], "directory": str(destination), "cases_requested": cases})
    return 2 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
