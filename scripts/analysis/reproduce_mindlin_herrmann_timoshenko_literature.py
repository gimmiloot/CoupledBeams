"""Bounded one-rectangle M-H/Timoshenko source audit; no production defaults.

New independent-contraction/dispersion contract cannot be a Bishop preset.
Reuses rectangular section/basis and Bishop workflow's atomic artifact writers.
Jang's unstated numeric kappa requires an explicit conditional-control input.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from importlib.metadata import version
import json
import math
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.analysis.reproduce_bishop_literature import sha, write_json, write_csv
from scripts.lib import mindlin_herrmann_longitudinal as mh

FIXTURE = ROOT / "data/input/mindlin_herrmann_timoshenko_sources.json"
OUTPUT = ROOT / "results/mindlin_herrmann_timoshenko_literature"


def source_check():
    fixture = json.loads(FIXTURE.read_text(encoding="utf-8"))
    if fixture["equations_version"] != mh.EQUATIONS_VERSION:
        raise ValueError("Equation/config version mismatch")
    checked = {}
    for name, entry in fixture["sources"].items():
        path = ROOT/entry["path"]
        if sha(path) != entry["sha256"]:
            raise ValueError("Source hash changed; re-audit " + entry["path"])
        checked[name] = {"path": entry["path"], "sha256": entry["sha256"], "pages": entry["pages"]}
    # Changing a fixture must not silently retain a source-reproduction label.
    for name in ("rucka", "jang"):
        printed = fixture["sources"][name]["printed_parameters"]
        parameters = fixture["cases"][name]["parameters_si"]
        conversion = {"E": ("E_GPa", 10**9), "rho": ("rho_kg_m3", 1),
                      "nu": ("nu", 1), "b": ("b_mm", Fraction(1, 1000)),
                      "h": ("h_mm", Fraction(1, 1000))}
        for field, (source_field, factor) in conversion.items():
            if float(parameters[field]) != float(Fraction(printed[source_field])*factor):
                raise ValueError(f"Source/SI parameter conflict: {name}.{field}")
    r = fixture["cases"]["rucka"]
    for field, source_field in (("mh_shear_factor", "K_MH1"), ("mh_inertia_factor", "K_MH2"),
                                ("tim_shear_factor", "K_Tim1")):
        if r[field] != float(Fraction(fixture["sources"]["rucka"]["printed_parameters"][source_field])):
            raise ValueError("Source correction-factor conflict: " + field)
    return fixture, checked


def provenance(checked, kappa_input):
    import hashlib
    def git(*args):
        r = subprocess.run(["git", *args], cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
        return r.stdout.strip() if r.returncode == 0 else "unavailable"
    files = [FIXTURE, Path(__file__), ROOT/"scripts/lib/mindlin_herrmann_longitudinal.py",
             ROOT/"scripts/lib/isotropic_rectangular_timoshenko_coupled_beams.py",
             ROOT/"scripts/analysis/reproduce_bishop_literature.py"]
    identity = {"schema": "mh-tim-results-v1", "sources": checked,
        "files": {p.relative_to(ROOT).as_posix(): sha(p) for p in files},
        "equations_version": mh.EQUATIONS_VERSION, "jang_kappa_input": kappa_input,
        "versions": {n: version(n) for n in ("numpy", "scipy", "matplotlib")},
        "python": platform.python_version(), "head": git("rev-parse", "HEAD")}
    fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
    return {"fingerprint": fingerprint, "identity": identity, "git_branch": git("branch", "--show-current"),
            "git_status": git("status", "--short"), "executable": sys.executable,
            "cwd": str(ROOT), "command": subprocess.list2cmdline([sys.executable, *sys.argv]),
            "spectrum_semantics": "source_acoustic_optical_dispersion_branches; not descendant tracking",
            "frequency_map_policy": "frequency-map-v1", "calculation_mode": "fast_plot",
            "policy_scope": "fixed literature dispersion curves; no geometry/frequency-eigenvalue map",
            "root_method": "analytic quadratic; fixed bounded frequency grids; zero iterative root searches"}


def make_model(fixture, variant, kappa=None, parameters=None):
    case = fixture["cases"]["rucka" if variant == "rucka_2010" else "jang"]
    params = parameters or case["parameters_si"]
    if variant == "rucka_2010":
        factors = (case["mh_shear_factor"], case["mh_inertia_factor"],
                   case["tim_shear_factor"], 12*case["tim_shear_factor"]/math.pi**2)
    elif variant == "jang_2014_bare_isotropic":
        if kappa is None:
            raise ValueError("Jang numeric kappa_b is not established; pass --jang-kappa explicitly")
        factors = (kappa, 1., kappa, 1.)
    else:
        raise ValueError("Unknown source variant: " + variant)
    return mh.source_model(params, mh_shear_factor=factors[0], mh_inertia_factor=factors[1],
        tim_shear_factor=factors[2], tim_rotary_factor=factors[3], variant=variant)


def model_record(model):
    s = model.section
    return {"variant": model.variant, "parameters_si": {"E": s.E, "rho": s.rho, "nu": s.nu,
        "b": s.width, "h": s.thickness}, "moments": mh.rectangle_moments(s.width, s.thickness),
        "correction_factors": {"K_MH1": model.mh_shear_factor, "K_MH2": model.mh_inertia_factor,
            "K_Tim1": s.K, "K_Tim2": model.tim_rotary_factor},
        "coefficients": model.coefficients, "coefficient_units_kg_m_s": mh.COEFFICIENT_UNITS,
        "limits": mh.limits(model)}


def verify_model(model, contract):
    """Independent D*E*D energy eigenproblem and Hellmann--Feynman derivative."""
    import numpy as np
    from scipy.linalg import eigh
    from scripts.lib import isotropic_rectangular_timoshenko_coupled_beams as timo
    elastic, mass = mh.energy_matrices(model)
    inv = 1/np.sqrt(np.diag(mass))
    records = []
    dp = mh.strain_operator(1.)-mh.strain_operator(0.)
    for k in contract["verification_k_per_m"]:
        stiffness, _ = mh.fourier_matrices(model, k)
        d = mh.strain_operator(k)
        dk = dp.conj().T@elastic@d + d.conj().T@elastic@dp
        normal = stiffness*inv[:, None]*inv[None, :]
        derivative = dk*inv[:, None]*inv[None, :]
        for indices, block in zip(((0, 1), (2, 3)), mh.blocks(model)):
            sub = normal[np.ix_(indices, indices)]
            deriv = derivative[np.ix_(indices, indices)]
            values, vectors = eigh(sub)
            scale = float(np.linalg.norm(sub, 2))
            for i, result in enumerate(block.temporal(k)):
                lam = result["omega_squared"]
                v = vectors[:, i]
                group = float(np.real(v.conj()@deriv@v))/(2*math.sqrt(lam))
                residual = float(np.linalg.norm(sub@v-lam*v))/(scale+abs(lam))
                spatial = block.spatial(result["frequency_hz"])[i]
                records.append({"k_per_m": k, "branch": result["branch"],
                    "equation_scaled_residual": residual,
                    "eigenvalue_scaled_error": abs(lam-values[i])/scale,
                    "eigenvalue_relative_error": abs(lam-values[i])/lam,
                    "spectral_scale_over_eigenvalue": scale/lam,
                    "group_hf_relative_error": abs(group-result["group_velocity_m_s"])/max(abs(group), 1e-30),
                    "spatial_roundtrip_relative_error": abs(spatial["wavenumber_per_m"]-k)/k,
                    "polynomial_scaled_residual": spatial["polynomial_scaled_residual"]})
        # One matrix built from the four-field energy, compared with both polynomials.
        full_values = np.linalg.eigvalsh(normal)
        union = sorted(r["omega_squared"] for b in mh.blocks(model) for r in b.temporal(k))
        if np.max(np.abs(full_values-union))/np.linalg.norm(normal, 2) > contract["eigenvalue_relative_tol"]:
            raise ArithmeticError("Four-field variational roots do not match source blocks")
        if np.any(normal[:2, 2:] != 0) or np.any(mass[:2, 2:] != 0):
            raise ArithmeticError("Source axial/bending coupling present; block implementation rejected")
    maximum = {field: max(r[field] for r in records) for field in (
        "equation_scaled_residual", "eigenvalue_scaled_error", "group_hf_relative_error",
        "spatial_roundtrip_relative_error", "polynomial_scaled_residual")}
    passed = (maximum["equation_scaled_residual"] <= contract["mass_normalized_equation_residual_tol"] and
        maximum["eigenvalue_scaled_error"] <= contract["eigenvalue_relative_tol"] and
        maximum["group_hf_relative_error"] <= contract["group_velocity_relative_tol"] and
        maximum["spatial_roundtrip_relative_error"] <= contract["spatial_roundtrip_relative_tol"] and
        maximum["polynomial_scaled_residual"] <= contract["polynomial_scaled_residual_tol"])
    project_checks = []
    if model.tim_rotary_factor == 1:
        block = mh.blocks(model)[1]
        for ratio in (.1, .9, 1., 1.1, 2.):
            f = ratio*block.cutoff_hz
            basis = timo.timoshenko_spatial_basis(2*math.pi*f, model.section)
            actual = sorted(-r["k_squared_per_m2"] for r in block.spatial(f))
            expected = sorted((basis.z_a, basis.z_b))
            error = max(abs(a-b)/max(1., abs(b)) for a, b in zip(actual, expected))
            project_checks.append({"frequency_hz": f, "relative_spatial_root_error": error})
        passed &= all(r["relative_spatial_root_error"] <= contract["project_spatial_root_relative_tol"] for r in project_checks)
    return {"status": "PASS" if passed else "UNRESOLVED", "maxima": maximum, "records": records,
        "project_bending_checks": project_checks, "mass_metric_scaling": inv.tolist(),
        "conditioning_note": "Eigenvalue error is normalized by spectral scale; tiny acoustic eigenvalues also retain raw relative errors and scale/lambda conditioning diagnostics."}


def dispersion_rows(model, grid):
    import numpy as np
    rows = []
    frequencies = np.linspace(grid["start"], grid["stop"], grid["count"]).tolist()
    frequencies += [b.cutoff_hz for b in mh.blocks(model) if grid["start"] <= b.cutoff_hz <= grid["stop"]]
    for f in sorted(set(frequencies)):
        for block in mh.blocks(model):
            try:
                roots = block.spatial(f)
            except ArithmeticError as error:
                raise ArithmeticError(f"{model.variant}, {block.labels}, f={f}: {error}") from error
            for row in roots:
                rows.append({"variant": model.variant, **row})
    return rows


def compute(case, fixture, kappa, directory):
    if case in ("jang", "comparison") and kappa is None:
        p = fixture["cases"]["rucka" if case == "comparison" else "jang"]["parameters_si"]
        write_json(directory/"summary.json", {"status": "SOURCE_NUMERIC_CONFIG_UNRESOLVED",
            "reason": "Jang kappa_b is not numerically specified in the audited source; explicit input required",
            "unconditional_axial_speed_m_s": math.sqrt(p["E"]/p["rho"]),
            "unconditional_contraction_cutoff_hz": math.sqrt(12*p["E"]/(p["rho"]*(1-p["nu"]**2)*p["h"]**2))/(2*math.pi),
            "production_status": "PRODUCTION_MH_COEFFICIENTS_UNRESOLVED"})
        return {"SOURCE_REPRODUCTION": "SOURCE_NUMERIC_CONFIG_UNRESOLVED",
                "PRODUCTION_COEFFICIENTS": "PRODUCTION_MH_COEFFICIENTS_UNRESOLVED"}
    names = ("rucka_2010", "jang_2014_bare_isotropic") if case == "comparison" else (
        "rucka_2010" if case == "rucka" else "jang_2014_bare_isotropic",)
    common = fixture["cases"]["rucka"]["parameters_si"] if case == "comparison" else None
    models = [make_model(fixture, n, kappa, common) for n in names]
    grid = fixture["cases"][case]["frequency_grid_hz"]
    rows = [r for m in models for r in dispersion_rows(m, grid)]
    write_csv(directory/"dispersion.csv", rows)
    diagnostics = {m.variant: verify_model(m, fixture["numerical_contract"]) for m in models}
    polynomial_max = max(r["polynomial_scaled_residual"] for r in rows)
    diagnostics["grid_polynomial_scaled_residual_max"] = polynomial_max
    write_json(directory/"diagnostics.json", diagnostics)
    summary = {"models": [model_record(m) for m in models], "source_grid_hz": grid,
        "production_status": "PRODUCTION_MH_COEFFICIENTS_UNRESOLVED",
        "mapping_status": "MH_SOURCE_VARIANTS_NOT_EQUIVALENT",
        "mapping": "Exact family mapping: K_MH1=kappa_b,K_MH2=1,K_Tim1=kappa_b,K_Tim2=1; Rucka published fitted factors do not satisfy it.",
        "jang_kappa": kappa, "jang_kappa_status": "EXPLICIT_CONDITIONAL_CONTROL; not source-recovered or production",
        "figure_status": "QUALITATIVE_REPRODUCTION" if case == "rucka" else "CONDITIONAL_QUALITATIVE_REPRODUCTION",
        "scientific_status": "MHTIM_VARIANT_DEPENDENT", "attempts": [],
        "spatial_roots_stored": len(rows), "quadratic_evaluations": len(rows)//2,
        "grid_rule": "declared bounded uniform grid plus exact in-range analytic cutoffs",
        "fourier_verification_points": len(fixture["numerical_contract"]["verification_k_per_m"])*len(models),
        "iterative_root_searches": 0, "retries": 0}
    if case == "rucka":
        m = models[0]
        import numpy as np
        p = m.coefficients
        speed = math.sqrt(m.section.E/m.section.rho)
        references = []
        for f in np.linspace(grid["start"], grid["stop"], grid["count"]):
            omega = 2*math.pi*f
            references.extend([
                {"branch": "elementary", "frequency_hz": f, "wavenumber_per_m": omega/speed,
                 "group_velocity_m_s": speed},
                {"branch": "Euler-Bernoulli", "frequency_hz": f,
                 "wavenumber_per_m": (p["m"]/p["B"])**.25*math.sqrt(omega),
                 "group_velocity_m_s": 2*(p["B"]/p["m"])**.25*math.sqrt(omega)}])
        write_csv(directory/"reference_dispersion.csv", references)
        checks = []
        for statement in fixture["sources"]["rucka"]["statements"]:
            block = mh.blocks(m)[0 if statement["block"] == "mh" else 1]
            interval = statement["frequency_interval_hz"]
            counts = [sum(r["state"] == "PROPAGATING" for r in block.spatial(f)) for f in interval]
            # Cutoff above the entire interval proves the count between endpoints as well.
            passed = counts == [statement["propagating_count"]]*2 and block.cutoff_hz > interval[1]
            checks.append({"source": statement, "endpoint_counts": counts,
                "cutoff_hz": block.cutoff_hz, "status": "PASS" if passed else "FAIL"})
        summary["source_statement_checks"] = checks
        samples = []
        for f in (100000, 120000):
            samples.extend(r for b in mh.blocks(m) for r in b.spatial(f) if r["state"] == "PROPAGATING")
        summary["calibration_frequency_samples"] = samples
        wave_speed = math.sqrt(m.section.E/m.section.rho)
        summary["elementary_velocity_m_s"] = wave_speed
        summary["eb_group_velocity_formula"] = "2*(EI/rhoA)^0.25*sqrt(omega)"
    if case == "comparison":
        comparison = []
        for k in fixture["cases"]["comparison"]["k_control_per_m"]:
            for m in models:
                for block in mh.blocks(m):
                    for r in block.temporal(k):
                        comparison.append({"variant": m.variant, "k_per_m": k, **r,
                            "phase_velocity_m_s": 2*math.pi*r["frequency_hz"]/k if k else None})
        write_csv(directory/"comparison.csv", comparison)
    numeric = (all(diagnostics[n]["status"] == "PASS" for n in names) and
               polynomial_max <= fixture["numerical_contract"]["polynomial_scaled_residual_tol"])
    source_pass = all(c["status"] == "PASS" for c in summary.get("source_statement_checks", []))
    summary["scientific_status"] = "MHTIM_VARIANT_DEPENDENT" if numeric and source_pass else "MHTIM_SOURCE_REPRODUCTION_FAILED"
    write_json(directory/"summary.json", summary)
    return {"SOURCE_TRANSCRIPTION": "VERIFIED_LOCAL_PDFS", "NUMERICAL_VERIFICATION": "PASS" if numeric else "UNRESOLVED",
        "SOURCE_STATEMENTS": ("PASS" if source_pass else "FAIL") if case == "rucka" else "NOT_NUMERICALLY_PRINTED",
        "SOURCE_FIGURE": summary["figure_status"], "SOURCE_REPRODUCTION": summary["scientific_status"],
        "MH_MAPPING": summary["mapping_status"], "PRODUCTION_COEFFICIENTS": summary["production_status"]}


def validate_bundle(directory, fingerprint):
    manifest = json.loads((directory/"manifest.json").read_text(encoding="utf-8"))
    if manifest["fingerprint"] != fingerprint:
        raise ValueError("Stale provenance: recompute with current inputs/code/versions")
    for name, digest in manifest["artifacts"].items():
        if sha(directory/name) != digest:
            raise ValueError("Artifact changed: " + name)
    return manifest


def render(case, directory):
    """Saved data only: no call to source_model, dispersion or verification."""
    import csv
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    summary = json.loads((directory/"summary.json").read_text(encoding="utf-8"))
    if "models" not in summary:
        return
    with (directory/"dispersion.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    fig, axes = plt.subplots(1, 2 if case == "jang" else 1, figsize=(11 if case == "jang" else 7, 4.5), squeeze=False)
    ax = axes[0, -1]
    colors = {"axial": "black", "contraction": "green", "bending": "red", "shear": "blue"}
    for m in summary["models"]:
        for label in colors:
            group = [r for r in rows if r["variant"] == m["variant"] and r["branch"] == label]
            frequencies = np.array([float(r["frequency_hz"]) for r in group])
            vg = np.array([float(r["group_velocity_m_s"]) if r["group_velocity_m_s"] else np.nan for r in group])
            style = "--" if m["variant"].startswith("jang") and case == "comparison" else "-"
            ax.plot(frequencies/(1e6 if case == "jang" else 1e3), vg/(1e3 if case == "jang" else 1), style,
                    color=colors[label], label=label+(" / "+m["variant"] if case == "comparison" else ""))
            if case == "jang":
                h = m["parameters_si"]["h"]
                axes[0, 0].plot(frequencies/1e6, [float(r["wavenumber_per_m"])*h for r in group], color=colors[label], label=label)
                axes[0, 0].plot(frequencies/1e6, [-float(r["attenuation_per_m"])*h for r in group], "--", color=colors[label])
    if case == "rucka":
        with (directory/"reference_dispersion.csv").open(encoding="utf-8") as f:
            reference = list(csv.DictReader(f))
        for label, color in (("elementary", "gray"), ("Euler-Bernoulli", "orange")):
            data = [r for r in reference if r["branch"] == label]
            ax.plot([float(r["frequency_hz"])/1e3 for r in data],
                    [float(r["group_velocity_m_s"]) for r in data], "--", color=color, label=label)
    for a in axes[0]:
        a.set_xlabel("Frequency (MHz)" if case == "jang" else "Frequency (kHz)")
        a.grid(alpha=.25)
        a.legend(fontsize=7)
    ax.set_ylabel("Group velocity (km/s)" if case == "jang" else "Group velocity (m/s)")
    if case == "jang":
        axes[0, 0].set_ylabel("Re(k h) / -attenuation h")
    if case == "rucka":
        title = "Rucka Fig.4: source factors; visual comparison"
    elif case == "jang":
        title = f"Jang Fig.9(a): conditional kappa={summary['jang_kappa']:.6g}; source value unstated"
    else:
        title = f"One common Rucka geometry; explicit Jang kappa={summary['jang_kappa']:.6g}"
        fig.set_size_inches(8, 5.6)
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.16), ncol=2, fontsize=7)
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(directory/(case+".png"), dpi=160)
    fig.savefig(directory/(case+".pdf"))
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check-sources", action="store_true")
    action.add_argument("--compute", action="store_true")
    action.add_argument("--plot-only", action="store_true")
    parser.add_argument("--case", choices=("rucka", "jang", "comparison", "all"), default="all")
    parser.add_argument("--variant", choices=("rucka_2010", "jang_2014_bare_isotropic"))
    parser.add_argument("--jang-kappa", help="Explicit conditional value, e.g. 5/6; never a recovered source/default")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args(argv)
    try:
        kappa = None if args.jang_kappa is None else float(Fraction(args.jang_kappa))
    except (ValueError, ZeroDivisionError):
        parser.error("--jang-kappa must be a number or fraction, e.g. 5/6")
    if kappa is not None and (not math.isfinite(kappa) or kappa <= 0):
        parser.error("--jang-kappa must be finite and positive")
    if args.variant:
        selected = "rucka" if args.variant == "rucka_2010" else "jang"
        if args.case not in (selected, "all"):
            parser.error("--variant conflicts with --case")
        args.case = selected
    fixture, checked = source_check()
    prov = provenance(checked, args.jang_kappa)
    if args.check_sources:
        print("VERIFIED_LOCAL_PDF_HASHES; explicit Jang kappa required for numerical control")
        return 0
    destination = args.output_dir/prov["fingerprint"][:16]
    failed = False
    cases = ("rucka", "jang", "comparison") if args.case == "all" else (args.case,)
    for case in cases:
        directory = destination/case
        if args.compute and not (directory/"manifest.json").exists():
            directory.mkdir(parents=True, exist_ok=True)
            try:
                statuses = compute(case, fixture, kappa, directory)
            except (ArithmeticError, ValueError) as error:
                write_json(directory/"failure.json", {"exception": type(error).__name__, "reason": str(error),
                    "attempt": 1, "retries": 0, "status": "UNRESOLVED"})
                statuses = {"NUMERICAL_VERIFICATION": "UNRESOLVED"}
            artifacts = {p.name: sha(p) for p in directory.iterdir() if p.is_file() and p.name != "manifest.json"}
            summary = json.loads((directory/"summary.json").read_text(encoding="utf-8")) if (directory/"summary.json").exists() else {}
            write_json(directory/"manifest.json", {**prov, "case": case, "config": fixture["cases"][case],
                "exact_model_parameter_sets": summary.get("models", []),
                "source_printed_data": {n: s.get("printed_parameters", {}) for n, s in fixture["sources"].items()},
                "numerical_contract": fixture["numerical_contract"], "statuses": statuses, "artifacts": artifacts})
        else:
            manifest = validate_bundle(directory, prov["fingerprint"])
            statuses = manifest["statuses"]
            if args.plot_only:
                render(case, directory)
            else:
                print(case+": reused validated data; zero root evaluations")
        failed |= (statuses.get("NUMERICAL_VERIFICATION") == "UNRESOLVED" or
                   statuses.get("SOURCE_STATEMENTS") == "FAIL" or
                   statuses.get("SOURCE_REPRODUCTION") == "SOURCE_NUMERIC_CONFIG_UNRESOLVED")
        print(case, statuses, directory)
    write_json(args.output_dir/"current.json", {"fingerprint": prov["fingerprint"], "directory": str(destination), "cases": cases})
    return 2 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
