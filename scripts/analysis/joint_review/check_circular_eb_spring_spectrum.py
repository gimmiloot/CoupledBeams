"""Bounded circular-EB diagnostic: existing physics and public spectrum solver.

Composes existing spectrum and physical-mode APIs; no new physical formulas.
--mode shapes reads saved roots and cannot call a root search.
--mode kappa-continuation adds only triggered, bounded stiffness midpoints.
"""
from __future__ import annotations

import argparse
import csv
from copy import deepcopy
from dataclasses import asdict
import hashlib
import inspect
import json
import math
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[3]
for directory in (ROOT, ROOT / "src"):
    if str(directory) not in sys.path:
        sys.path.insert(0, str(directory))

import numpy as np

from my_project.analytic.formulas import (
    BeamParams, assemble_clamped_coupled_matrix, lambdas_to_frequencies,
    segment_lengths,
)
from scripts.lib.inplane_rotational_spring_eb import (
    EBArm, Joint, boundary_assembly, endpoint_diagnostics, scalar_joint_residuals,
)
from scripts.lib import inplane_rotational_spring_eb_modes as modes
from scripts.lib.reddy_symmetric_coupled_beams import boundary_matrix_diagnostics
from scripts.lib.general_spectrum_completeness import (
    Geometry, SearchSettings, resolve_matrix_spectrum,
)

OUTPUT = ROOT / "results/joint_review/circular_eb_spring_general_spectrum"
K_TARGET, K_GUARD = 6, 7
# Existing EB spring-mode gates; the requested min/next gate is explicit.
# Sources: eb_modes.recover, track_inplane_rotational_spring_eb.CRITERIA,
# check_inplane_spring_robustness.CRITERIA (quadrature_rtol).
SHAPE_GATES = dict(sigma_ratio=1e-9, physical_residual=1e-9,
    compatibility=1e-10, boundary_residual=1e-9, mass_error=1e-6,
    MAC=.95, margin=.20)
SOURCE_PATHS = (
    "src/my_project/analytic/formulas.py",
    "scripts/lib/inplane_rotational_spring_eb.py",
    "scripts/lib/general_spectrum_completeness.py",
    "scripts/analysis/joint_review/check_circular_eb_spring_spectrum.py",
)
REUSE = """# Circular EB spring spectrum: diagnostic smoke-test

Existing APIs reused directly:

- `formulas.BeamParams`, `segment_lengths`, `lambdas_to_frequencies`,
  `assemble_clamped_coupled_matrix`: circular baseline and normalization.
- `inplane_rotational_spring_eb.EBArm`, `Joint`, `boundary_assembly`:
  existing section-independent EB mechanics; only A=E*S, D=E*I, m=rho*S
  are passed from BeamParams to EBArm.
- `general_spectrum_completeness.SearchSettings()`, `Geometry`,
  `resolve_matrix_spectrum`: unchanged public solver and all default settings,
  with no seeds. Twelve requested sorted positions; internal reserve 20/24.

This separate diagnostic output contract composes existing components; it is
not a new physical/helper module or a preset of the old pilot's root search.
No kernel, solver, threshold, old pilot or old result is changed. This bounded
smoke-test is not a frequency map or article result.

Parameters: E=2.1e11 Pa, rho=7800 kg/m^3, r=.005 m, L_total=2 m,
l=BeamParams.L_base=1 m, epsilon=BeamParams.eps=.0025. Lengths come from
segment_lengths. Lambda is the baseline frequency variable;
f=lambdas_to_frequencies([Lambda], params)[0], omega=2*pi*f, Omega=Lambda^2.
Finite stiffness is k_theta=kappa_theta*reference.D/reference.L;
RIGID is exactly Joint("RIGID"). No mass renormalization is performed.

Only sorted eigenvalues are compared: sorted_position is not a descendant
branch. No shapes, cross-state modal MAC, Delta psi or sensitivity are computed.
Native algebraic self_MAC/null_vector fields from the generic solver are saved
unchanged; they are internal search checks, not physical mode reconstruction.
Native sigma_ratio means sigma_1/sigma_2, not sigma_min/sigma_max.

Scientific scope: K_target=6, K_guard=7, only (mu,beta)=(.30,15 degrees).
Gate A compares baseline/exact RIGID only over positions 1..7; the full
12-root comparison and its historical failure remain explicit qualifications.
Gate B: kappa=1; Gate C: kappa=10,100. Only a target-prefix failure stops this
pilot. Native spectrum_status is never overwritten or reinterpreted as a pass.
Numerical comparisons use unchanged root_match_tol=2e-4. The prefix audit
uses the native margin max(seed_half_width,2*scan_step)=.03 above guard 7.
Saved results are reused after checking parameters, settings, library hashes
and the unchanged provider source. No old RIGID recalculation is required.
"""


def source_hashes():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCE_PATHS}


def providers(params, mu, beta_deg, kappa_or_mode="RIGID"):
    """Compose baseline and EB-kernel matrix providers; no new physics."""
    lengths = segment_lengths(params, mu)
    resultants = dict(A=params.E * params.S, D=params.E * params.I, m=params.rho * params.S)
    arm1, arm2 = (EBArm(**resultants, L=length) for length in lengths)
    reference = EBArm(**resultants, L=params.L_base)
    joint = (Joint("RIGID") if kappa_or_mode == "RIGID" else
             Joint("SPRING", float(kappa_or_mode) * reference.D / reference.L))
    beta_rad = math.radians(beta_deg)

    def baseline(value):
        return assemble_clamped_coupled_matrix(value, beta_rad, mu, params.eps)

    def kernel(value):
        frequency_hz = lambdas_to_frequencies(np.array([value]), params)[0]
        omega = 2 * np.pi * frequency_hz
        return boundary_assembly(omega, arm1, arm2, beta_rad, joint, reference).dimensionless

    kernel.case = (arm1, arm2, beta_rad, joint, reference)
    return baseline, kernel


def provider_source_hash():
    """Preserve the saved physics fingerprint; ignore only metadata exposure.

    The sole added line attaches the already constructed arms/joint for shapes.
    All original parameter, normalization and matrix-provider code is unchanged.
    Full entry-point hashes still include that line and all postprocessing code.
    """
    source = inspect.getsource(providers).replace(
        "    kernel.case = (arm1, arm2, beta_rad, joint, reference)\n", "")
    return hashlib.sha256(source.encode()).hexdigest()


def finite_json(value):
    """Keep nonfinite native diagnostics explicit in standards-compliant JSON."""
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [finite_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    return value


def write_csv(path, rows, fields):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def target_prefix(result):
    """Read native records only; no matrix calls, new root gates or recovery."""
    settings = result["settings"]
    margin = max(settings["seed_half_width"], 2 * settings["scan_step"])
    configurations = {name: result[name] for name in ("primary", "verification")}
    groups = [result["roots"], *(c["roots"] for c in configurations.values())]
    enough = all(len(roots) >= K_GUARD for roots in groups)
    limit = max(roots[K_GUARD - 1]["Lambda"] for roots in groups) + margin if enough else math.inf
    reasons, higher, checks = [], [], {}
    if not enough:
        reasons.append("missing_target_or_guard")
    for name, config in configurations.items():
        roots = config["roots"]
        prefix = roots[:K_GUARD]
        ordered = all(math.isfinite(r["Lambda"]) for r in roots) and all(
            a["Lambda"] <= b["Lambda"] for a, b in zip(prefix, prefix[1:]))
        ordered = ordered and [r["sorted_index"] for r in prefix] == list(range(1, len(prefix) + 1))
        if enough:
            ordered = ordered and all(r["Lambda"] >= prefix[-1]["Lambda"] for r in roots[K_GUARD:])
        if not ordered:
            reasons.append(name + ":incorrect_order")
        quality = []
        for root in prefix:
            evidence = [c for c in config["candidates"] if c["Lambda"] == root["Lambda"]]
            accepted = root["acceptance_status"] == "accepted_full_matrix_svd" and any(
                c["acceptance_status"] == "accepted_full_matrix_svd"
                and c["diagnostics"]["finite_matrix_status"] == "finite" for c in evidence)
            # Multiple-root quality is the native candidate acceptance decision;
            # the SVD/nullity calculation itself is not duplicated here.
            quality.append(accepted and math.isfinite(root["sigma_1"])
                and root["sigma_1"] <= settings["sigma_accept"]
                and ((root["detected_nullity"] == 1 and root["multiplicity_status"] == "simple_root"
                      and math.isfinite(root["sigma_ratio"]) and root["sigma_ratio"] <= settings["sigma_ratio_accept"])
                     or (root["detected_nullity"] == 2 and root["multiplicity_status"] == "verified_nullity_2")))
        if not all(quality):
            reasons.append(name + ":native_root_quality_failure")
        # A saved accepted root above the native audit margin witnesses actual
        # evaluation there, even if an incomplete run's lambda_upper grew after
        # its last scan. This does not infer a new search interval or root.
        witnesses = [r["Lambda"] for r in roots if limit < r["Lambda"] <= config["lambda_upper"]]
        covered = enough and config["lambda_upper"] > limit and bool(witnesses)
        if not covered:
            reasons.append(name + ":guard_coverage_not_established")
        unresolved = []
        for entry in config["unresolved_intervals"]:
            try:
                lower = float(entry.split(":", 1)[0])
            except (TypeError, ValueError):
                lower = math.nan
            unresolved.append((lower, entry))
        for row in config["interval_rows"]:
            status = row["resolution_status"]
            if "unresolved" in status or status == "pending":
                unresolved.append((float(row["Lambda_left"]), row))
        for candidate in config["candidates"]:
            if candidate["acceptance_status"] != "accepted_full_matrix_svd":
                unresolved.append((float(candidate["interval_left"]), candidate))
        blocking = [entry for lower, entry in unresolved if not math.isfinite(lower) or lower <= limit]
        above = [entry for lower, entry in unresolved if math.isfinite(lower) and lower > limit]
        if blocking:
            reasons.append(name + ":unresolved_below_guard_audit_limit")
        unusual_tail = [dict(sorted_position=r["sorted_index"], Lambda=r["Lambda"],
                            detected_nullity=r["detected_nullity"], cluster_id=r["root_cluster_id"],
                            cluster_size=r["cluster_size"], multiplicity_status=r["multiplicity_status"])
                        for r in roots[K_GUARD:] if r["detected_nullity"] > 1 or r["root_cluster_id"]]
        if above or unusual_tail:
            higher.append(dict(configuration=name, unresolved_above_guard_margin=above,
                               native_tail_multiplicity_or_cluster=unusual_tail))
        checks[name] = dict(has_seven=len(roots) >= K_GUARD, ordered=ordered,
            native_root_quality=quality, blocking_unresolved=blocking, lambda_upper=config["lambda_upper"],
            guard_coverage_pass=covered, accepted_coverage_witness=min(witnesses) if witnesses else None)
    comparisons = result["primary_vs_verification"][:K_GUARD]
    agreement = len(comparisons) == K_GUARD and all(
        row["status"] == "pass" and row["multiplicity_agreement"] for row in comparisons)
    pairs = zip(configurations["primary"]["roots"][:K_GUARD], configurations["verification"]["roots"][:K_GUARD])
    multiplicity = enough and all(all(a[k] == b[k] for k in (
        "detected_nullity", "track_multiplicity", "multiplicity_status")) for a, b in pairs)
    if not agreement:
        reasons.append("target_primary_verification_disagreement")
    if not multiplicity:
        reasons.append("target_multiplicity_disagreement")
    higher_comparisons = [row for row in result["primary_vs_verification"][K_GUARD:] if row["status"] != "pass"]
    if higher_comparisons or result["exclusion_reason"] or result["root12_boundary_warning"]:
        higher.append(dict(native_solver_exclusion_reason=result["exclusion_reason"],
            root12_boundary_warning=result["root12_boundary_warning"], higher_comparisons=higher_comparisons))
    return dict(solver_spectrum_status=result["spectrum_status"],
        target_prefix_status="TARGET_PREFIX_FAIL" if reasons else "TARGET_PREFIX_PASS",
        K_target=K_TARGET, K_guard=K_GUARD, native_audit_margin=margin, guard_audit_limit=limit,
        primary_verification_agreement=agreement, multiplicity_agreement=multiplicity,
        checks=checks, failure_reasons=reasons,
        higher_spectrum_status="HIGHER_SPECTRUM_QUALIFICATION" if higher else "NONE",
        higher_spectrum_qualifications=higher,
        candidate_audit_scope="native retained candidate union and native unresolved interval records; raw rejected detections are not exported by this public API")


def rigid_comparison(baseline, rigid):
    rows = []
    tolerance = baseline["settings"]["root_match_tol"]
    for j, (a, b) in enumerate(zip(baseline["roots"], rigid["roots"]), 1):
        error = abs(a["Lambda"] - b["Lambda"])
        rows.append(dict(case="main", mu=.30, beta_deg=15., sorted_position=j,
            baseline_Lambda=a["Lambda"], rigid_Lambda=b["Lambda"], absolute_difference=error,
            relative_difference=error / a["Lambda"], within_root_match_tol=error <= tolerance,
            baseline_nullity=a["detected_nullity"], rigid_nullity=b["detected_nullity"],
            baseline_cluster=a["root_cluster_id"], rigid_cluster=b["root_cluster_id"],
            baseline_spectrum_status=baseline["spectrum_status"], rigid_spectrum_status=rigid["spectrum_status"],
            scope="TARGET_PREFIX" if j <= K_GUARD else "HIGHER_SPECTRUM"))
    prefix_pass = len(rows) >= K_GUARD and all(r["within_root_match_tol"]
        and r["baseline_nullity"] == r["rigid_nullity"] for r in rows[:K_GUARD])
    prefix_pass = prefix_pass and all(target_prefix(r)["target_prefix_status"] == "TARGET_PREFIX_PASS" for r in (baseline, rigid))
    summary = dict(target_prefix_status="RIGID_TARGET_PREFIX_PASS" if prefix_pass else "RIGID_TARGET_PREFIX_FAIL",
        full12_equivalent=len(rows) == 12 and all(r["within_root_match_tol"] and r["baseline_nullity"] == r["rigid_nullity"] for r in rows),
        higher_mismatched_positions=[r["sorted_position"] for r in rows[K_GUARD:]
            if not r["within_root_match_tol"] or r["baseline_nullity"] != r["rigid_nullity"]])
    for count in (K_GUARD, 12):
        for name, field in (("absolute", "absolute_difference"), ("relative", "relative_difference")):
            summary[f"first{count}_max_{name}"] = max(r[field] for r in rows[:count]) if len(rows) >= count else None
    return summary, rows


def historical_event(result, params):
    """Locate the old Hz value in native event records, without a local search."""
    selected = []
    for name in ("primary", "verification"):
        candidates = result[name]["candidates"]
        if not candidates:
            return dict(status="NO_NATIVE_CANDIDATES")
        selected.append(min(candidates, key=lambda c: abs(float(
            lambdas_to_frequencies(np.array([c["Lambda"]]), params)[0]) - 204.287)))
    lo = min(c["interval_left"] for c in selected)
    hi = max(c["interval_right"] for c in selected)
    output = dict(locator_hz=204.287, native_bracket_union=[lo, hi], configurations={})
    for name in ("primary", "verification"):
        events = []
        for candidate in result[name]["candidates"]:
            if lo <= candidate["Lambda"] <= hi:
                records = [r for r in result[name]["roots"] if r["Lambda"] == candidate["Lambda"]]
                events.append(dict(Lambda=candidate["Lambda"], frequency_hz=float(lambdas_to_frequencies(
                    np.array([candidate["Lambda"]]), params)[0]), sorted_positions=[r["sorted_index"] for r in records],
                    root_records=records, interval_left=candidate["interval_left"], interval_right=candidate["interval_right"],
                    detection_sources=candidate["detection_sources"], acceptance_status=candidate["acceptance_status"]))
        output["configurations"][name] = dict(native_root_event_count=len(events), events=events)
    indices = {i for c in output["configurations"].values() for e in c["events"] for i in e["sorted_positions"]}
    output["in_target_prefix"] = bool(indices) and all(i <= K_GUARD for i in indices)
    output["native_comparisons"] = [r for r in result["primary_vs_verification"] if r["sorted_index"] in indices]
    return output


def trend_rows(spectra, tolerance):
    rows = []
    for j in range(K_TARGET):
        values = [spectra[key]["roots"][j]["Lambda"] for key in ("1", "10", "100", "RIGID")]
        distances = [abs(value - values[-1]) for value in values[:-1]]
        monotone = all(b >= a - tolerance for a, b in zip(values, values[1:]))
        if max(values) - min(values) <= tolerance:
            convergence = "LOW_SENSITIVITY_SORTED_POSITION"
        elif all(a - b > tolerance for a, b in zip(distances, distances[1:])):
            convergence = "DECREASING_RESOLVED"
        elif distances[-1] < distances[0] - tolerance and all(b <= a + tolerance for a, b in zip(distances, distances[1:])):
            convergence = "DECREASING_WITH_UNRESOLVED_STEPS"
        else:
            convergence = "CONVERGENCE_NOT_CONFIRMED"
        rows.append(dict(sorted_position=j + 1, **dict(zip(
            ("Lambda_k1", "Lambda_k10", "Lambda_k100", "Lambda_rigid"), values)), **dict(zip(
            ("delta_vs_rigid_k1", "delta_vs_rigid_k10", "delta_vs_rigid_k100"), [x / values[-1] for x in distances])),
            monotonic_status="NONDECREASING_WITHIN_ROOT_MATCH_TOL" if monotone else "VIOLATION",
            convergence_status=convergence))
    return rows


def save(output, data):
    output.mkdir(parents=True, exist_ok=True)
    (output / "diagnostics.json").write_text(json.dumps(finite_json(data), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_csv(output / "rigid_equivalence.csv", data.get("rigid_rows", []), (
        "case", "mu", "beta_deg", "sorted_position", "baseline_Lambda", "rigid_Lambda",
        "absolute_difference", "relative_difference", "within_root_match_tol", "baseline_nullity",
        "rigid_nullity", "baseline_cluster", "rigid_cluster", "baseline_spectrum_status", "rigid_spectrum_status", "scope"))
    spring_rows = []
    params = BeamParams(**data["params"])
    for state, key in (("1", "main_kappa_1"), ("10", "main_kappa_10"), ("100", "main_kappa_100"), ("RIGID", "main_RIGID")):
        if key not in data["runs"]:
            continue
        entry = data["runs"][key]
        result = entry["result"]
        for root, comparison in zip(result["roots"], result["primary_vs_verification"]):
            spring_rows.append(dict(joint_mode="RIGID" if state == "RIGID" else "SPRING",
                kappa_theta="" if state == "RIGID" else int(state), sorted_position=root["sorted_index"],
                Lambda=root["Lambda"], frequency_hz=float(lambdas_to_frequencies(np.array([root["Lambda"]]), params)[0]),
                detected_nullity=root["detected_nullity"], cluster_id=root["root_cluster_id"], cluster_size=root["cluster_size"],
                sigma_1=root["sigma_1"], sigma_2=root["sigma_2"], sigma_ratio=root["sigma_ratio"],
                primary_verification_difference=comparison["absolute_difference"],
                detection_sources=";".join(root["detection_sources"]), solver_spectrum_status=result["spectrum_status"],
                target_prefix_status=entry["prefix"]["target_prefix_status"]))
    write_csv(output / "spring_spectrum.csv", spring_rows, (
        "joint_mode", "kappa_theta", "sorted_position", "Lambda", "frequency_hz", "detected_nullity",
        "cluster_id", "cluster_size", "sigma_1", "sigma_2", "sigma_ratio", "primary_verification_difference",
        "detection_sources", "solver_spectrum_status", "target_prefix_status"))
    trends = data.get("sorted_trends", [])
    write_csv(output / "trend_summary.csv", trends, ("sorted_position", "Lambda_k1", "Lambda_k10", "Lambda_k100",
        "Lambda_rigid", "delta_vs_rigid_k1", "delta_vs_rigid_k10", "delta_vs_rigid_k100", "monotonic_status", "convergence_status"))
    lines = [REUSE, "\n## Status\n", f"Workflow: `{data['status']}`. Gates: `{json.dumps(data['gates'])}`.\n",
        "| Provider | Native solver status | Target prefix | P/V agree 1..7 | Blocking reasons | Higher spectrum |",
        "|---|---|---|---|---|---|"]
    for key, entry in data["runs"].items():
        prefix = entry["prefix"]
        lines.append(f"| {key} | {prefix['solver_spectrum_status']} | {prefix['target_prefix_status']} | {prefix['primary_verification_agreement']} | {prefix['failure_reasons']} | {prefix['higher_spectrum_status']} |")
    if "rigid_comparison" in data:
        comparison = data["rigid_comparison"]
        lines += ["\n## RIGID equivalence and retained qualification\n", f"`{json.dumps(comparison)}`\n"]
        if comparison["target_prefix_status"] == "RIGID_TARGET_PREFIX_PASS":
            lines.append("Exact RIGID equivalence is confirmed for the scientific target prefix consisting of the first six sorted eigenvalues plus guard root 7. A separate higher-spectrum SVD-nullity classification discrepancy remains at positions 11?12 and is outside the scope of the present pilot.\n")
        lines += ["The full 12-root comparison is NOT accepted as physical equivalence. Native RIGID",
            "nullity=2 at Lambda=18.13952880970414 occupies positions 11?12, whereas baseline",
            "has a simple root there and Lambda=19.65271588480334 at position 12. These are",
            "native numerical classifications, not a new proof of physical multiplicity.",
            "The previous STOP_GATE_A, full CSV/report and 69-pass/1-fail regression record",
            "are preserved under previous_full_spectrum_gate in diagnostics.json.\n"]
    lines += ["\n## Prefix evidence\n", "| Provider/configuration | Guard audit limit | Reported upper bound | Accepted root beyond audit limit | Blocking unresolved records |",
        "|---|---|---|---|---|"]
    for key, entry in data["runs"].items():
        for name, checks in entry["prefix"]["checks"].items():
            lines.append(f"| {key}/{name} | {entry['prefix']['guard_audit_limit']:.10g} | {checks['lambda_upper']:.10g} | {checks['accepted_coverage_witness']} | {len(checks['blocking_unresolved'])} |")
    lines += ["\nPrefix quality uses native accepted_full_matrix_svd records and native SVD",
        "limits, native P/V comparison rows, ordering, matching multiplicity, and native",
        "unresolved records below guard + .03. An accepted higher root witnesses actual",
        "coverage beyond that margin. Nonlocalized unresolved records fail the prefix.",
        "Higher warnings and all native search records remain in JSON. Candidate completeness",
        "is assessed within the public API: its retained union contains accepted candidates;",
        "it does not export every raw rejected detection. No new raw-candidate audit is claimed.\n"]
    if "historical_204Hz_event" in data:
        event = data["historical_204Hz_event"]
        lines += ["\n## Historical 204.287 Hz area\n", f"Native bracket union in Lambda: `{event.get('native_bracket_union')}`; in prefix: `{event.get('in_target_prefix')}`."]
        for name, item in event.get("configurations", {}).items():
            for found in item["events"]:
                classifications = [(r["multiplicity_status"], r["cluster_size"]) for r in found["root_records"]]
                lines.append(f"{name}: {item['native_root_event_count']} native event(s); sorted positions {found['sorted_positions']}; f={found['frequency_hz']:.12g} Hz; Lambda={found['Lambda']:.12g}; {classifications}.")
        lines.append(f"Native P/V comparison: `{json.dumps(event.get('native_comparisons', []))}`\n")
    if trends:
        lines += ["\n## First six sorted positions\n", "| j | Lambda(1) | Lambda(10) | Lambda(100) | Lambda(RIGID) | delta(1) | delta(10) | delta(100) | Monotonicity | Convergence |",
            "|---|---|---|---|---|---|---|---|---|---|"]
        for row in trends:
            numbers = " | ".join(f"{row[k]:.10g}" for k in ("Lambda_k1", "Lambda_k10", "Lambda_k100", "Lambda_rigid", "delta_vs_rigid_k1", "delta_vs_rigid_k10", "delta_vs_rigid_k100"))
            lines.append(f"| {row['sorted_position']} | {numbers} | {row['monotonic_status']} | {row['convergence_status']} |")
        lines += ["\nMonotonicity is a post-check only: no root was chosen, changed or repaired using it.",
            "LOW_SENSITIVITY_SORTED_POSITION, if present, means the four Lambda values",
            "span no more than the existing root_match_tol; it is not a modal sensitivity",
            "calculation and gives no information about Delta psi. Cluster members, if any,",
            "are interpreted as sorted cluster positions, not identified physical modes.\n"]
    else:
        lines.append("\nFour-state frequency trend: NOT_EVALUATED. Stop reasons are retained above.\n")
    lines += ["\n## Reproduction and limits\n", "`python -B scripts/analysis/joint_review/check_circular_eb_spring_spectrum.py`",
        "\nCheckpoint reuse is missing-only; unchanged native results and historical provenance are preserved.",
        "Only the frequency tendency is tested. No relative rotation, real monolithic-joint",
        "stiffness, r/L relationship or applicability to a real construction is established.",
        "No new physical/helper module, root solver, matrix, threshold or shape analysis."]
    (output / "report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run(output=OUTPUT):
    output = Path(output)
    params = BeamParams(E=2.1e11, rho=7800., r=.005, L_total=2.)
    settings = SearchSettings()
    sources = source_hashes()
    provider_hash = provider_source_hash()
    checkpoint = output / "diagnostics.json"
    if checkpoint.exists():
        data = json.loads(checkpoint.read_text(encoding="utf-8"))
        if (data["params"] != asdict(params) or data["settings"] != asdict(settings)
            or data.get("provider_source_sha256") != provider_hash
            or any(data["source_hashes"][p] != sources[p] for p in SOURCE_PATHS[:-1])):
            raise ValueError("Saved physics/provider/settings provenance changed; no automatic recomputation")
        # The prior full-spectrum audit remains historical, not a current gate.
        data.pop("audit", None)
    else:
        data = dict(params=asdict(params), settings=asdict(settings), source_hashes=sources,
                    provider_source_sha256=provider_hash, runs={})
    data.update(workflow_source_hashes=sources, spectrum_semantics="sorted_positions",
        scope=dict(K_target=K_TARGET, K_guard=K_GUARD, mu=.30, beta_deg=15., states=[1, 10, 100, "RIGID"]),
        status="RUNNING_GATE_A", gates=dict(A="PENDING", B="NOT_RUN", C="NOT_RUN"), solver_calls_this_invocation=0)
    for key, entry in data["runs"].items():
        if entry["result"]["geometry"] != asdict(Geometry(params.eps, 15., .30, 0.)):
            raise ValueError(f"Unexpected saved geometry in {key}")
        entry.setdefault("source_hashes", data["source_hashes"])
        entry["prefix"] = target_prefix(entry["result"])
    data.pop("sorted_trends", None)
    save(output, data)

    def solve(key, provider):
        if key in data["runs"]:
            print(f"Reusing {key}: no solver call", flush=True)
            return data["runs"][key]["result"]
        print(f"Starting {key}: unchanged public solver defaults", flush=True)
        started = time.perf_counter()
        result = asdict(resolve_matrix_spectrum(provider, settings=settings,
            geometry=Geometry(params.eps, 15., .30, 0.), model=key))
        data["solver_calls_this_invocation"] += 1
        data["runs"][key] = dict(elapsed_seconds=time.perf_counter() - started, result=result,
            source_hashes=sources, prefix=target_prefix(result))
        save(output, data)
        print(f"{key}: solver={result['spectrum_status']}; target={data['runs'][key]['prefix']['target_prefix_status']}", flush=True)
        return result

    baseline_provider, rigid_provider = providers(params, .30, 15.)
    baseline = solve("main_baseline", baseline_provider)
    rigid = solve("main_RIGID", rigid_provider)
    data["rigid_comparison"], data["rigid_rows"] = rigid_comparison(baseline, rigid)
    data["gates"]["A"] = data["rigid_comparison"]["target_prefix_status"]
    if data["gates"]["A"] != "RIGID_TARGET_PREFIX_PASS":
        data["status"] = "STOP_GATE_A_TARGET_PREFIX"
        save(output, data)
        return data
    spectra = {"RIGID": rigid}
    for kappa in (1, 10, 100):
        data["status"] = "RUNNING_GATE_B" if kappa == 1 else "RUNNING_GATE_C"
        _, provider = providers(params, .30, 15., kappa)
        key = f"main_kappa_{kappa}"
        result = solve(key, provider)
        spectra[str(kappa)] = result
        status = data["runs"][key]["prefix"]["target_prefix_status"]
        data["gates"]["B" if kappa == 1 else "C"] = status
        if kappa == 1:
            data["historical_204Hz_event"] = historical_event(result, params)
        if status != "TARGET_PREFIX_PASS":
            data["status"] = f"STOP_TARGET_PREFIX_kappa_{kappa}"
            break
    else:
        data["status"] = "COMPLETED_TARGET_PREFIX"
        data["sorted_trends"] = trend_rows(spectra, settings.root_match_tol)
    save(output, data)
    return data


def saved_shape_inputs(output):
    """Fail closed on missing/provenance-invalid data; never search for roots."""
    path = Path(output) / "diagnostics.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    params = BeamParams(**data["params"])
    if (data["params"] != asdict(BeamParams(2.1e11, 7800., .005, 2.))
        or data["provider_source_sha256"] != provider_source_hash()
        or any(data["source_hashes"][p] != source_hashes()[p] for p in SOURCE_PATHS[:-1])):
        raise ValueError("SAVED_FREQUENCY_PROVENANCE_MISMATCH")
    results = {}
    for state, key in (("1", "main_kappa_1"), ("10", "main_kappa_10"),
                       ("100", "main_kappa_100"), ("RIGID", "main_RIGID")):
        result = data["runs"][key]["result"]
        if (result["geometry"] != asdict(Geometry(params.eps, 15., .30, 0.))
            or data["runs"][key]["prefix"]["target_prefix_status"] != "TARGET_PREFIX_PASS"
            or len(result["roots"]) < K_GUARD):
            raise ValueError(f"MISSING_SAVED_TARGET_PREFIX:{state}")
        roots = result["roots"][:K_GUARD]
        if any(r["detected_nullity"] != 1 or r["sorted_index"] != j
               or not math.isfinite(r["Lambda"]) for j, r in enumerate(roots, 1)):
            raise ValueError(f"INCONSISTENT_SAVED_SIMPLE_ROOTS:{state}")
        results[state] = roots
    return params, results, data


def reconstruct_saved_mode(params, state, root, case):
    """Compose existing EB reaction recovery, arm fields and mass metric."""
    arm1, arm2, beta, joint, reference = case
    frequency = float(lambdas_to_frequencies(np.array([root["Lambda"]]), params)[0])
    omega = 2 * np.pi * frequency
    assembly = boundary_assembly(omega, arm1, arm2, beta, joint, reference)
    endpoint = endpoint_diagnostics(assembly, beta, joint, reference)
    singular = boundary_matrix_diagnostics(assembly.dimensionless).scaled_singular_values
    ratio = float(singular[-1] / singular[-2])
    row = dict(joint_mode=joint.mode, kappa_theta="" if state == "RIGID" else float(state),
        sorted_position=root["sorted_index"], Lambda=root["Lambda"], frequency_hz=frequency,
        detected_nullity=endpoint["nullity"], sigma_ratio=ratio,
        endpoint_sigma_min_over_max=endpoint["sigma_ratio"], reconstruction_status="FAILED", shape_key="")
    diagnostic = dict(endpoint=endpoint, scaled_singular_values=singular.tolist(), failures=[])
    failures = diagnostic["failures"]
    if endpoint["nullity"] != 1:
        failures.append("RECONSTRUCTION_NULLITY_NOT_ONE")
    if not math.isfinite(ratio) or ratio > SHAPE_GATES["sigma_ratio"]:
        failures.append("FULL_MATRIX_MIN_OVER_NEXT_GATE")
    if endpoint["sigma_ratio"] > SHAPE_GATES["sigma_ratio"]:
        failures.append("ENDPOINT_ROOT_GATE")
    if failures:
        return row, None, diagnostic
    record = endpoint["vectors"][0]
    reactions = np.asarray(record["physical_clamp_reactions"]).reshape(2, 3)
    xi, weights = modes.quadrature(129)
    states = np.array([modes.arm_states(omega, arm, reaction, xi)
                       for arm, reaction in zip((arm1, arm2), reactions)])
    components = [modes.mass_vector(y[None, ...], arm, weights)
                  for y, arm in zip(states, (arm1, arm2))]
    vector = np.concatenate(components)
    mass = float(np.vdot(vector, vector).real)
    diagnostic["arm_masses_before_normalization"] = [float(np.vdot(v, v).real) for v in components]
    if not math.isfinite(mass) or mass <= 0:
        failures.append("NONPOSITIVE_MODAL_MASS")
        return row, None, diagnostic
    factor = np.sqrt(mass)
    states, reactions, vector = states / factor, reactions / factor, vector / factor
    mass_after = float(np.vdot(vector, vector).real)
    ends = states[:, -1, :].ravel()
    # Same residual normalization as eb_modes.recover, with one fixed reference.
    units = np.tile([reference.L, reference.L, 1., reference.D/reference.L**2,
                     reference.D/reference.L**2, reference.D/reference.L], 2)
    amplitude = float(np.max(abs(ends / units)))
    residuals = scalar_joint_residuals(ends / amplitude, beta, joint) / assembly.row_units
    hat = reactions.ravel() / assembly.reaction_scales
    boundary = float(np.linalg.norm(assembly.dimensionless @ hat) /
                     (np.linalg.norm(assembly.dimensionless) * np.linalg.norm(hat)))
    delta = float(states[0, -1, 2] - states[1, -1, 2])
    diagnostic.update(normalized_reconstructed_physical_residuals=residuals.tolist(),
        mass_normalized_physical_residuals=scalar_joint_residuals(ends, beta, joint).tolist(),
        endpoint_amplitude_in_reference_units=amplitude,
        normalized_endpoint_residual_max=max(abs(np.asarray(record["normalized_physical_residuals"]))),
        endpoint_boundary_residual=record["boundary_residual"],
        endpoint_scaled_residual=record["scaled_residual"])
    physical = max(float(max(abs(residuals))), diagnostic["normalized_endpoint_residual_max"])
    compatible = max(float(max(abs(residuals[:2]))), max(abs(np.asarray(record["normalized_physical_residuals"])[:2])))
    null = max(boundary, record["boundary_residual"], record["scaled_residual"])
    if physical > SHAPE_GATES["physical_residual"]:
        failures.append("RECONSTRUCTED_PHYSICAL_GATE")
    if compatible > SHAPE_GATES["compatibility"]:
        failures.append("TRANSLATIONAL_COMPATIBILITY_GATE")
    if null > SHAPE_GATES["boundary_residual"]:
        failures.append("RECONSTRUCTED_NULL_GATE")
    if abs(mass_after - 1) > SHAPE_GATES["mass_error"]:
        failures.append("MASS_NORMALIZATION_GATE")
    key = f"k{state}_p{root['sorted_index']:02d}"
    row.update(mass_before_normalization=mass, mass_after_normalization=mass_after,
        mass_normalization_error=abs(mass_after - 1), psi1_mass_normalized=float(states[0, -1, 2]),
        psi2_mass_normalized=float(states[1, -1, 2]), Delta_psi_mass_normalized=delta,
        abs_Delta_psi_mass_normalized=abs(delta),
        s="NOT_APPLICABLE_RIGID" if state == "RIGID" else reference.D/reference.L * delta**2/(omega**2 * mass_after),
        max_physical_residual=physical, max_compatibility_residual=compatible, boundary_residual=null,
        reconstruction_status="FAILED" if failures else "CONFIRMED", shape_key=key)
    return row, dict(states=states, reactions=reactions, physical_vector=vector), diagnostic


def assign_physical_modes(left, right):
    """Delegate the full rectangular assignment to the existing mass MAC API."""
    return modes.assign([r["physical_vector"] for r in left],
                        [r["physical_vector"] for r in right])


def track_reconstructed_modes(by_state, arrays, state_status):
    """Local seeds only; unresolved links are never propagated as descendants."""
    seeds = [f"k1_seed_{j:02d}" for j in range(1, K_TARGET + 1)]
    paths = {seed: {} for seed in seeds}
    mappings, assignments = [], []
    if state_status["1"] == "CONFIRMED":
        for seed, row in zip(seeds, by_state["1"][:K_TARGET]):
            paths[seed]["1"] = row
    active = [seed for seed in seeds if paths[seed]]
    for source, target in (("1", "10"), ("10", "100"), ("100", "RIGID")):
        if not active or state_status[target] != "CONFIRMED":
            break
        previous = [paths[seed][source] for seed in active]
        candidates = by_state[target]  # All seven, including the guard.
        left = [arrays[row["shape_key"]] for row in previous]
        right = [arrays[row["shape_key"]] for row in candidates]
        columns, mac, margins = assign_physical_modes(left, right)
        assignments.append(dict(source=source, target=target, seeds=list(active),
            source_sorted_positions=[r["sorted_position"] for r in previous],
            candidate_sorted_positions=[r["sorted_position"] for r in candidates],
            columns=columns.tolist(), MAC_matrix=mac.tolist(), margins=margins.tolist(),
            metric="physical_mass_vectors", frequency_in_cost=False))
        next_active = []
        for i, seed in enumerate(active):
            j = int(columns[i])
            row = candidates[j]
            confirmed = mac[i, j] >= SHAPE_GATES["MAC"] and margins[i] >= SHAPE_GATES["margin"]
            phase = 1
            if confirmed:
                shape = right[j]
                if np.vdot(left[i]["physical_vector"], shape["physical_vector"]).real < 0:
                    phase = -1
                    for field in ("states", "reactions", "physical_vector"):
                        shape[field] *= -1
                    for field in ("psi1_mass_normalized", "psi2_mass_normalized", "Delta_psi_mass_normalized"):
                        row[field] *= -1
                row["phase_sign"] = phase
                paths[seed][target] = row
                next_active.append(seed)
            mappings.append(dict(branch_id=seed, source_kappa=source, target_kappa=target,
                source_sorted_position=previous[i]["sorted_position"], target_sorted_position=row["sorted_position"],
                MAC=float(mac[i, j]), margin=float(margins[i]), phase_sign=phase,
                mapping_status="CONFIRMED" if confirmed else "UNRESOLVED"))
        active = next_active
    trends, descendants = [], []
    for seed in seeds:
        path = paths[seed]
        confirmed = len(path) == 4
        row = dict(branch_id=seed, tracking_status="CONFIRMED" if confirmed else "UNRESOLVED",
                   rotation_trend_status="UNRESOLVED", sensitivity_trend_status="UNRESOLVED")
        for state, suffix in (("1", "k1"), ("10", "k10"), ("100", "k100"), ("RIGID", "rigid")):
            if state not in path:
                continue
            found = path[state]
            row["sorted_position_" + suffix] = found["sorted_position"]
            row["abs_Delta_psi_" + suffix] = found["abs_Delta_psi_mass_normalized"]
            if state != "RIGID":
                row["s_" + suffix] = found["s"]
            descendants.append(dict(branch_id=seed, state=state, Lambda=found["Lambda"],
                frequency_hz=found["frequency_hz"], sorted_position=found["sorted_position"],
                abs_Delta_psi_mass_normalized=found["abs_Delta_psi_mass_normalized"], s=found["s"],
                full_path_status=row["tracking_status"]))
        for mapping in mappings:
            if mapping["branch_id"] == seed:
                row[f"MAC_{mapping['source_kappa']}_to_{mapping['target_kappa'].lower()}"] = mapping["MAC"]
        if confirmed:
            rotation = [path[state]["abs_Delta_psi_mass_normalized"] for state in ("1", "10", "100", "RIGID")]
            sensitivity = [path[state]["s"] for state in ("1", "10", "100")]
            row["rotation_trend_status"] = ("DECREASING" if all(b <= a for a, b in zip(rotation, rotation[1:]))
                else "NONMONOTONE_APPROACH_TO_RIGID" if rotation[-1] < rotation[0] else "NO_DECREASE")
            row["sensitivity_trend_status"] = ("DECREASING" if all(b <= a for a, b in zip(sensitivity, sensitivity[1:])) else "NONMONOTONE")
        trends.append(row)
    return mappings, trends, assignments, descendants


def shapes_run(output=OUTPUT):
    """No call path to the spectral solver: saved-root postprocessing only."""
    output = Path(output)
    params, saved, spectrum = saved_shape_inputs(output)
    frequency_files = ("diagnostics.json", "rigid_equivalence.csv", "spring_spectrum.csv", "trend_summary.csv")
    input_hashes = {name: hashlib.sha256((output / name).read_bytes()).hexdigest() for name in frequency_files}
    rows, arrays, details, by_state, status = [], {}, {}, {}, {}
    for state, roots in saved.items():
        _, provider = providers(params, .30, 15., state)
        by_state[state] = []
        status[state] = "CONFIRMED"
        for root in roots:
            row, shape, diagnostic = reconstruct_saved_mode(params, state, root, provider.case)
            row["phase_sign"] = 1
            rows.append(row)
            by_state[state].append(row)
            details[f"k{state}_p{root['sorted_index']:02d}"] = diagnostic
            if shape is not None:
                arrays[row["shape_key"]] = shape
            if row["reconstruction_status"] != "CONFIRMED":
                status[state] = "FAILED"
                break  # Stop this state; never correct a saved root.
        print(f"shapes {state}: {status[state]}, {len(by_state[state])}/7 evaluated", flush=True)
    mappings, trends, assignments, descendants = track_reconstructed_modes(by_state, arrays, status)
    confirmed = [r for r in rows if r["reconstruction_status"] == "CONFIRMED"]
    rigid_rows = [r for r in confirmed if r["joint_mode"] == "RIGID"]
    summary = dict(state_status=status, attempted_modes=len(rows), confirmed_modes=len(confirmed),
        root_search_calls=0, boundary_assemblies=len(rows),
        max_physical_residual=max((r["max_physical_residual"] for r in confirmed), default=None),
        max_compatibility_residual=max((r["max_compatibility_residual"] for r in confirmed), default=None),
        max_boundary_residual=max((r["boundary_residual"] for r in confirmed), default=None),
        max_sigma_min_over_next=max((r["sigma_ratio"] for r in rows), default=None),
        max_mass_normalization_error=max((r["mass_normalization_error"] for r in confirmed), default=None),
        rigid_confirmed_count=len(rigid_rows),
        max_abs_Delta_psi_rigid=max((r["abs_Delta_psi_mass_normalized"] for r in rigid_rows), default=None),
        confirmed_descendants=[r["branch_id"] for r in trends if r["tracking_status"] == "CONFIRMED"],
        unresolved_descendants=[r["branch_id"] for r in trends if r["tracking_status"] != "CONFIRMED"],
        min_attempted_MAC=min((r["MAC"] for r in mappings), default=None),
        min_attempted_margin=min((r["margin"] for r in mappings), default=None),
        sorted_position_changes=[r for r in mappings if r["mapping_status"] == "CONFIRMED"
                                 and r["source_sorted_position"] != r["target_sorted_position"]])
    for name, digest in input_hashes.items():
        if hashlib.sha256((output / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"FREQUENCY_INPUT_CHANGED:{name}")
    payload = dict(summary=summary, parameters=spectrum["params"], geometry=spectrum["scope"],
        frequency_input_sha256=input_hashes, current_source_sha256=source_hashes(),
        provider_physics_sha256=provider_source_hash(), gates=SHAPE_GATES,
        gate_sources=["inplane_rotational_spring_eb_modes.recover", "track_inplane_rotational_spring_eb.CRITERIA",
                      "check_inplane_spring_robustness.CRITERIA", "explicit user sigma_min/sigma_next <= 1e-9"],
        shape_nodes=129, mass_metric="concatenate mass_vector of each actual arm; u,w only; M=1",
        sigma_ratio_definition="min/next from existing positively-equilibrated boundary_matrix_diagnostics",
        endpoint_sigma_ratio_definition="min/max; retained separately, never substituted for min/next",
        sensitivity="(reference.D/reference.L)*abs(Delta_psi)^2/(omega^2*M); RIGID not applicable",
        phase_convention="raw endpoint diagnostics retain SVD phase; normalized CSV/NPZ fields use applied phase_sign",
        branch_id_scope="local k1_seed identifiers only; not canonical project branch_id",
        reconstruction=rows, endpoint_and_field_diagnostics=details,
        assignments=assignments, mapping=mappings, trends=trends, descendant_records=descendants)
    extra = ("scripts/lib/inplane_rotational_spring_eb_modes.py", "scripts/lib/reddy_symmetric_coupled_beams.py")
    payload["current_source_sha256"].update({p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in extra})
    (output / "shape_diagnostics.json").write_text(json.dumps(finite_json(payload), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    fields = ("joint_mode", "kappa_theta", "sorted_position", "Lambda", "frequency_hz",
        "mass_before_normalization", "mass_after_normalization", "mass_normalization_error",
        "psi1_mass_normalized", "psi2_mass_normalized", "Delta_psi_mass_normalized", "abs_Delta_psi_mass_normalized",
        "s", "detected_nullity", "sigma_ratio", "endpoint_sigma_min_over_max", "max_physical_residual",
        "max_compatibility_residual", "boundary_residual", "reconstruction_status", "shape_key", "phase_sign")
    write_csv(output / "mode_reconstruction.csv", [{k: r.get(k, "") for k in fields} for r in rows], fields)
    write_csv(output / "mode_mapping.csv", mappings, ("branch_id", "source_kappa", "target_kappa",
        "source_sorted_position", "target_sorted_position", "MAC", "margin", "phase_sign", "mapping_status"))
    fields = ("branch_id", "sorted_position_k1", "sorted_position_k10", "sorted_position_k100", "sorted_position_rigid",
        "MAC_1_to_10", "MAC_10_to_100", "MAC_100_to_rigid", "abs_Delta_psi_k1", "abs_Delta_psi_k10",
        "abs_Delta_psi_k100", "abs_Delta_psi_rigid", "s_k1", "s_k10", "s_k100", "tracking_status",
        "rotation_trend_status", "sensitivity_trend_status")
    write_csv(output / "joint_rotation_trend.csv", [{k: r.get(k, "") for k in fields} for r in trends], fields)
    xi, weights = modes.quadrature(129)
    archive = {f"{key}__{field}": value for key, shape in arrays.items() for field, value in shape.items()}
    np.savez_compressed(output / "shapes.npz", xi=xi, quadrature_weights=weights, **archive)
    marker = "\n## Physical mode post-processing\n"
    original = (output / "report.md").read_text(encoding="utf-8").split(marker)[0]
    lines = [marker, "This section follows the preceding frequency-only stage. The prior sorted-spectrum",
        "result and higher-spectrum qualification remain unchanged. Command: `--mode shapes`.\n",
        "Existing endpoint_diagnostics recovers physical reactions; existing arm_states and",
        "mass_vector use each actual arm (L1=.7 m, L2=1.3 m). Existing quadrature(129)",
        "sets the mass metric. No parity projection or generic-solver coefficient vector is used.",
        "M=1 for the whole structure; Delta psi values below depend on this normalization.",
        "Physical residual normalization follows eb_modes.recover. Algebraic diagnostics reuse",
        "boundary_matrix_diagnostics; this is not a Reddy/RLB physical calculation.\n",
        f"Summary: `{json.dumps(summary)}`\n", f"Unchanged gates: `{json.dumps(SHAPE_GATES)}`\n",
        "endpoint_diagnostics.sigma_ratio is min/max. The separate user-required min/next",
        "ratio is obtained from existing equilibrated SVD diagnostics and must also pass 1e-9.",
        "Endpoint nullity must equal 1 with the existing rank_rtol=1e-12.\n",
        "MAC uses only physical mass vectors. Candidate pool: positions 1..7.",
        "Local k1_seed IDs are not canonical project branch IDs. An unresolved link stops",
        "that descendant; no intermediate stiffness or frequency cost is introduced.\n",
        "| Seed | Link | Positions | MAC | Margin | Status |", "|---|---|---|---|---|---|"]
    for r in mappings:
        lines.append(f"| {r['branch_id']} | {r['source_kappa']} -> {r['target_kappa']} | {r['source_sorted_position']} -> {r['target_sorted_position']} | {r['MAC']:.8g} | {r['margin']:.8g} | {r['mapping_status']} |")
    lines += ["\nOnly full-path CONFIRMED descendants enter the following trend table.",
        "| Seed | abs Delta psi (1) | (10) | (100) | RIGID | s(1) | s(10) | s(100) | Rotation | s trend |",
        "|---|---|---|---|---|---|---|---|---|---|"]
    for r in trends:
        if r["tracking_status"] == "CONFIRMED":
            numbers = " | ".join(f"{r[k]:.10g}" for k in ("abs_Delta_psi_k1", "abs_Delta_psi_k10", "abs_Delta_psi_k100", "abs_Delta_psi_rigid", "s_k1", "s_k10", "s_k100"))
            lines.append(f"| {r['branch_id']} | {numbers} | {r['rotation_trend_status']} | {r['sensitivity_trend_status']} |")
    lines += ["\nRotation/sensitivity monotonicity is descriptive, never an acceptance gate.",
        "RIGID s is NOT_APPLICABLE_RIGID. Failed reconstruction details, partial confirmed",
        "paths and all assignment matrices remain in shape_diagnostics.json. No failed",
        "root is refined or replaced. Source frequency files are unchanged byte for byte.",
        "No inference about real monolithic-joint stiffness or r/L follows from this 1D check."]
    (output / "report.md").write_text(original + "\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary), flush=True)
    return payload


def geometric_midpoint(left, right):
    """Finite positive stiffness only; exact RIGID has no logarithmic midpoint."""
    if left == "RIGID" or right == "RIGID":
        raise ValueError("RIGID_IS_NOT_A_FINITE_KAPPA")
    if not 0 < float(left) < float(right) < math.inf:
        raise ValueError("INVALID_FINITE_INTERVAL")
    return math.sqrt(float(left) * float(right))


def continuation_assignment(source, candidates, source_state, target_state, depth, attempts):
    """Existing global mass-MAC assignment; metadata never enter its cost."""
    if len(candidates) != K_GUARD:
        raise ValueError("CONTINUATION_REQUIRES_GUARD_7")
    columns, mac, margins = assign_physical_modes(source, candidates)
    selected, rows = [], []
    for i, column in enumerate(columns):
        target = deepcopy(candidates[column])
        accepted = bool(mac[i, column] >= SHAPE_GATES["MAC"] and margins[i] >= SHAPE_GATES["margin"])
        sign = -1 if accepted and np.vdot(source[i]["physical_vector"], target["physical_vector"]).real < 0 else 1
        for key in ("states", "reactions", "physical_vector"):
            target[key] *= sign
        for key in ("psi1_mass_normalized", "psi2_mass_normalized", "Delta_psi_mass_normalized"):
            target[key] *= sign
        target["branch_id"] = source[i]["branch_id"]
        selected.append(target)
        rows.append(dict(branch_id=target["branch_id"], source_kappa=source_state, target_kappa=target_state,
            source_sorted_position=source[i]["sorted_position"], target_sorted_position=target["sorted_position"],
            MAC=float(mac[i, column]), margin=float(margins[i]), phase_sign=sign,
            step_status="CONFIRMED" if accepted else "UNRESOLVED", refinement_depth=depth,
            direct_endpoint_agreement="", used_in_path=True))
    attempts.append(dict(source_kappa=source_state, target_kappa=target_state,
        candidate_sorted_positions=[r["sorted_position"] for r in candidates],
        source_branch_ids=[r["branch_id"] for r in source], MAC_matrix=mac.tolist(), rows=rows))
    return selected, rows


def continue_finite_interval(source, left, right, point, attempts, *, depth=0, force_midpoint=False):
    """At most midpoint plus failed quarter intervals, never a regular sweep.

    A terminal failed best candidate is only tentative diagnostic scaffolding
    for the remaining global assignments. Its path stays unresolved and is
    never admitted as a confirmed descendant to the next stiffness decade.
    """
    if left == "RIGID" or right == "RIGID":
        raise ValueError("RIGID_IS_NOT_A_FINITE_KAPPA")
    if depth not in (0, 1, 2):
        raise ValueError("CONTINUATION_DEPTH_BUDGET_EXCEEDED")
    if force_midpoint and depth != 0:
        raise ValueError("MIDPOINT_MUST_START_AT_DEPTH_ZERO")
    if not force_midpoint:
        selected, records = continuation_assignment(source, point(right), left, right, depth, attempts)
        if all(r["step_status"] == "CONFIRMED" for r in records) or depth == 2:
            return selected, [records], [left, right]
        for record in records:
            record["used_in_path"] = False
    middle = geometric_midpoint(left, right)
    selected, first, path1 = continue_finite_interval(source, left, middle, point, attempts, depth=depth+1)
    selected, second, path2 = continue_finite_interval(selected, middle, right, point, attempts, depth=depth+1)
    return selected, first + second, path1 + path2[1:]


def continue_rigid_endpoint(source, point, attempts):
    """One permitted bridge, 1000, only after the direct 100 -> RIGID fails."""
    selected, direct = continuation_assignment(source, point("RIGID"), 100., "RIGID", 0, attempts)
    if all(r["step_status"] == "CONFIRMED" for r in direct):
        return selected, [direct], [100., "RIGID"]
    for record in direct:
        record["used_in_path"] = False
    middle, first = continuation_assignment(source, point(1000.), 100., 1000., 1, attempts)
    selected, second = continuation_assignment(middle, point("RIGID"), 1000., "RIGID", 1, attempts)
    return selected, [first, second], [100., 1000., "RIGID"]


def continuation_status(direct, final, records):
    """Keep the historical direct verdict, including low-MAC best candidates."""
    agreement = ("DIRECT_ENDPOINT_AGREEMENT" if final["sorted_position"] == direct["target_sorted_position"]
                 else "DIRECT_ENDPOINT_DISAGREEMENT")
    if any(r["step_status"] != "CONFIRMED" for r in records):
        return "CONTINUATION_UNRESOLVED", agreement
    if agreement == "DIRECT_ENDPOINT_DISAGREEMENT":
        return "PATH_MAPPING_CONFLICT", agreement
    return ("DIRECT_UNRESOLVED_CONTINUATION_CONFIRMED" if direct["mapping_status"] == "UNRESOLVED"
            else "CONTINUATION_CONFIRMED"), agreement


CONTINUATION_INPUTS = ("diagnostics.json", "rigid_equivalence.csv", "spring_spectrum.csv", "trend_summary.csv",
    "mode_reconstruction.csv", "mode_mapping.csv", "joint_rotation_trend.csv", "shape_diagnostics.json", "shapes.npz")


def saved_continuation_inputs(output):
    params, roots, _ = saved_shape_inputs(output)
    payload = json.loads((output / "shape_diagnostics.json").read_text(encoding="utf-8"))
    if payload["gates"] != SHAPE_GATES or payload["provider_physics_sha256"] != provider_source_hash():
        raise ValueError("SAVED_SHAPE_PROVENANCE_MISMATCH")
    for name, digest in payload["current_source_sha256"].items():
        if name != SOURCE_PATHS[-1] and hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"SAVED_SHAPE_SOURCE_CHANGED:{name}")
    for name, digest in payload["frequency_input_sha256"].items():
        if hashlib.sha256((output / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"SAVED_SHAPE_FREQUENCY_CHANGED:{name}")
    points = {state: [] for state in roots}
    with np.load(output / "shapes.npz", allow_pickle=False) as archive:
        for row in payload["reconstruction"]:
            state = "RIGID" if row["joint_mode"] == "RIGID" else f"{row['kappa_theta']:g}"
            if row["reconstruction_status"] != "CONFIRMED":
                raise ValueError("SAVED_RECONSTRUCTION_NOT_CONFIRMED")
            mode = deepcopy(row)
            for field in ("states", "reactions", "physical_vector"):
                mode[field] = archive[row["shape_key"] + "__" + field].copy()
            points[state].append(mode)
    for state, pool in points.items():
        pool.sort(key=lambda r: r["sorted_position"])
        if len(pool) != K_GUARD or any(r["sorted_position"] != j or r["Lambda"] != roots[state][j-1]["Lambda"]
                                     for j, r in enumerate(pool, 1)):
            raise ValueError(f"SAVED_SHAPES_ROOTS_INCONSISTENT:{state}")
    direct = {r["branch_id"]: r for r in payload["mapping"] if r["source_kappa"] == "1" and r["target_kappa"] == "10"}
    if len(direct) != K_TARGET:
        raise ValueError("MISSING_DIRECT_1_TO_10_RECORDS")
    return params, points, direct


def continuation_run(output=OUTPUT):
    """Bounded continuation, separate artifacts; original endpoint files immutable."""
    output = Path(output)
    params, points, direct = saved_continuation_inputs(output)
    hashes = {n: hashlib.sha256((output / n).read_bytes()).hexdigest() for n in CONTINUATION_INPUTS}
    checkpoint = output / "kappa_continuation_diagnostics.json"
    settings = SearchSettings()
    data = json.loads(checkpoint.read_text(encoding="utf-8")) if checkpoint.exists() else dict(points={})
    if data.get("input_sha256", hashes) != hashes or data.get("settings", asdict(settings)) != asdict(settings):
        raise ValueError("CONTINUATION_CACHE_PROVENANCE_MISMATCH")
    data.setdefault("initial_source_sha256", data.get("source_sha256", source_hashes()))
    data.pop("stop_reason", None)
    data.update(input_sha256=hashes, settings=asdict(settings), gates=SHAPE_GATES,
        source_sha256=source_hashes(), provider_physics_sha256=provider_source_hash(),
        DIRECT_1_TO_10=direct, status="RUNNING", attempts=[], trends=[], solver_calls_this_invocation=0,
        scope=dict(mu=.30, beta_deg=15., K_target=6, K_guard=7, maximum_depth=2,
            finite_decade_budget=3, rigid_bridge=1000., symmetry_classes=False, frequency_in_cost=False),
        precedent="docs/laminated_beams/inplane_rotational_spring_eb_hinge_comparison.md; prepare_seed_mapping/continue_seed_mapping",
        phase_convention="point arrays retain reconstruction phase; each mapping phase_sign multiplies the target point arrays",
        branch_id_scope="local k1_seed IDs, not canonical project branches")

    def checkpoint_save():
        checkpoint.write_text(json.dumps(finite_json(data), indent=2, allow_nan=False) + "\n", encoding="utf-8")

    def point(state):
        key = "RIGID" if state == "RIGID" else f"{float(state):g}"
        if key in points:
            return points[key]
        # Whitelist only the authorized midpoints/quarters and the sole RIGID bridge.
        allowed = [10**q for q in (.25, .5, .75, 1.25, 1.5, 1.75)] + [1000.]
        canonical = next((x for x in allowed if math.isclose(float(state), x, rel_tol=2e-15)), None)
        if canonical is None:
            raise ValueError(f"UNAUTHORIZED_CONTINUATION_POINT:{state}")
        key = repr(canonical)
        if key in points:
            return points[key]
        if key not in data["points"]:
            _, provider = providers(params, .30, 15., canonical)
            checkpoint_save()  # Retain the triggering MAC failure before a long solve.
            print(f"continuation: new kappa={canonical:.17g}; unchanged public solver", flush=True)
            started = time.perf_counter()
            result = asdict(resolve_matrix_spectrum(provider, settings=settings,
                geometry=Geometry(params.eps, 15., .30, 0.), model=f"circular_EB_continuation_k{key}"))
            data["solver_calls_this_invocation"] += 1
            entry = dict(kappa_theta=canonical, result=result, prefix=target_prefix(result),
                source_sha256=source_hashes(), elapsed_seconds=time.perf_counter()-started,
                reconstruction=[], shape_arrays={}, diagnostics={})
            data["points"][key] = entry
            checkpoint_save()
        entry = data["points"][key]
        entry.setdefault("source_sha256", data["initial_source_sha256"])
        if entry["prefix"]["target_prefix_status"] != "TARGET_PREFIX_PASS":
            raise ValueError(f"STOP_TARGET_PREFIX:{key}")
        if any(r["reconstruction_status"] != "CONFIRMED" for r in entry["reconstruction"]):
            raise ValueError(f"STOP_RECONSTRUCTION:{key}")
        _, provider = providers(params, .30, 15., canonical)
        for root in entry["result"]["roots"][len(entry["reconstruction"]):K_GUARD]:
            row, shape, details = reconstruct_saved_mode(params, key, root, provider.case)
            entry["reconstruction"].append(row)
            entry["diagnostics"][str(root["sorted_index"])] = details
            if shape is not None:
                entry["shape_arrays"][row["shape_key"]] = {k: v.tolist() for k, v in shape.items()}
            checkpoint_save()
            if row["reconstruction_status"] != "CONFIRMED":
                raise ValueError(f"STOP_RECONSTRUCTION:{key}:{root['sorted_index']}")
        if len(entry["reconstruction"]) != K_GUARD or any(r["reconstruction_status"] != "CONFIRMED" for r in entry["reconstruction"]):
            raise ValueError(f"STOP_RECONSTRUCTION:{key}")
        points[key] = [dict(r, **{k: np.array(v) for k, v in entry["shape_arrays"][r["shape_key"]].items()})
                       for r in entry["reconstruction"]]
        print(f"continuation kappa={canonical:.10g}: TARGET_PREFIX_PASS; 7/7 forms confirmed", flush=True)
        return points[key]

    seeds = [dict(deepcopy(r), branch_id=f"k1_seed_{j:02d}") for j, r in enumerate(points["1"][:K_TARGET], 1)]
    paths = {r["branch_id"]: {"1": r} for r in seeds}
    records = {r["branch_id"]: [] for r in seeds}
    trends = {seed: dict(branch_id=seed, direct_1_to_10_status=direct[seed]["mapping_status"],
        continuation_1_to_10_status="NOT_RUN", status_10_to_100="NOT_RUN", status_100_to_rigid="NOT_RUN",
        full_path_status="CONTINUATION_UNRESOLVED") for seed in paths}

    def retain_steps(steps):
        for step in steps:
            for record in step:
                records[record["branch_id"]].append(record)

    try:
        selected, steps, path = continue_finite_interval(seeds, 1., 10., point, data["attempts"], force_midpoint=True)
        retain_steps(steps)
        active = []
        for mode in selected:
            seed = mode["branch_id"]
            status, agreement = continuation_status(direct[seed], mode, records[seed])
            trends[seed].update(continuation_1_to_10_status=status, path_to_10=path, direct_endpoint_agreement=agreement)
            for record in records[seed]:
                if record["target_kappa"] == 10.:
                    record["direct_endpoint_agreement"] = agreement
            if status in ("CONTINUATION_CONFIRMED", "DIRECT_UNRESOLVED_CONTINUATION_CONFIRMED"):
                paths[seed]["10"] = mode
                active.append(mode)
            elif status == "PATH_MAPPING_CONFLICT":
                trends[seed]["full_path_status"] = status
        for right, name in ((100., "10_to_100"), ("RIGID", "100_to_rigid")):
            if not active:
                break
            if right == "RIGID":
                selected, steps, path = continue_rigid_endpoint(active, point, data["attempts"])
            else:
                selected, steps, path = continue_finite_interval(active, 10., 100., point, data["attempts"])
            retain_steps(steps)
            active = []
            for mode in selected:
                seed = mode["branch_id"]
                ok = all(r["step_status"] == "CONFIRMED" for step in steps for r in step if r["branch_id"] == seed)
                suffix = "rigid" if right == "RIGID" else "100"
                trends[seed].update({f"status_{name}": "CONTINUATION_CONFIRMED" if ok else "CONTINUATION_UNRESOLVED",
                                     f"path_to_{suffix}": path})
                if ok:
                    paths[seed]["RIGID" if right == "RIGID" else "100"] = mode
                    active.append(mode)
        data["status"] = "COMPLETED_BOUNDED_CONTINUATION"
    except ValueError as error:
        data.update(status="STOP_CONTINUATION", stop_reason=str(error))
        print(str(error), flush=True)
    data["descendant_path_records"] = []
    for seed, path in paths.items():
        row = trends[seed]
        full = "RIGID" in path
        if full:
            row["full_path_status"] = "FULL_KAPPA_PATH_CONFIRMED"
        for state, mode in path.items():
            suffix = "rigid" if state == "RIGID" else f"k{state}"
            row[f"sorted_position_{suffix}"] = mode["sorted_position"]
            # Partial paths retain metadata, never an unconfirmed physical trend.
            if full:
                row[f"abs_Delta_psi_{suffix}"] = mode["abs_Delta_psi_mass_normalized"]
                if state != "RIGID":
                    row[f"s_{suffix}"] = mode["s"]
        accepted = [r for r in records[seed] if r["step_status"] == "CONFIRMED"]
        if accepted:
            row.update(minimum_step_MAC=min(r["MAC"] for r in accepted), minimum_step_margin=min(r["margin"] for r in accepted))
        if full:
            for field, states, result_name in (("abs_Delta_psi_mass_normalized", ("1", "10", "100", "RIGID"), "rotation"),
                                              ("s", ("1", "10", "100"), "sensitivity")):
                values = [path[s][field] for s in states]
                row[f"{result_name}_trend_status"] = "DECREASING" if all(b <= a for a, b in zip(values, values[1:])) else "NONMONOTONE"
            # Include every accepted intermediate state, without another recovery.
            ordered_modes = [(1., path["1"])] + [(r["target_kappa"],
                next(m for m in point(r["target_kappa"]) if m["sorted_position"] == r["target_sorted_position"]))
                for r in records[seed]]
            for state, mode in ordered_modes:
                data["descendant_path_records"].append(dict(branch_id=seed, kappa_or_mode=state,
                    **{k: mode[k] for k in ("sorted_position", "Lambda", "frequency_hz", "shape_key", "abs_Delta_psi_mass_normalized", "s")}))
            for field, result_name in (("abs_Delta_psi_mass_normalized", "rotation"), ("s", "sensitivity")):
                values = [mode[field] for state, mode in ordered_modes if field != "s" or state != "RIGID"]
                row[f"{result_name}_including_intermediates"] = "DECREASING" if all(b <= a for a, b in zip(values, values[1:])) else "NONMONOTONE"
    data["trends"] = list(trends.values())
    data["mapping"] = [r for attempt in data["attempts"] for r in attempt["rows"]]
    accepted = [r for r in data["mapping"] if r["used_in_path"] and r["step_status"] == "CONFIRMED"]
    data["summary"] = dict(new_kappas=[p["kappa_theta"] for p in data["points"].values()],
        full_confirmed_count=sum(r["full_path_status"] == "FULL_KAPPA_PATH_CONFIRMED" for r in trends.values()),
        minimum_accepted_MAC=min((r["MAC"] for r in accepted), default=None),
        minimum_accepted_margin=min((r["margin"] for r in accepted), default=None),
        endpoint_conflicts=[s for s, r in trends.items() if r["continuation_1_to_10_status"] == "PATH_MAPPING_CONFLICT"],
        sorted_position_changes=[r for r in accepted if r["source_sorted_position"] != r["target_sorted_position"]],
        decreasing_rotation_count=sum(r.get("rotation_trend_status") == "DECREASING" for r in trends.values()),
        decreasing_sensitivity_count=sum(r.get("sensitivity_trend_status") == "DECREASING" for r in trends.values()))
    new_modes = [r for p in data["points"].values() for r in p["reconstruction"] if r["reconstruction_status"] == "CONFIRMED"]
    data["summary"]["new_modes_confirmed"] = len(new_modes)
    data["summary"]["total_new_physical_points"] = len(data["points"])
    data["summary"]["new_point_maxima"] = {key: max((r[key] for r in new_modes), default=None)
        for key in ("sigma_ratio", "max_physical_residual", "max_compatibility_residual", "boundary_residual", "mass_normalization_error")}
    data["summary"]["decreasing_rotation_including_intermediates"] = sum(r.get("rotation_including_intermediates") == "DECREASING" for r in trends.values())
    data["summary"]["decreasing_sensitivity_including_intermediates"] = sum(r.get("sensitivity_including_intermediates") == "DECREASING" for r in trends.values())
    for name, digest in hashes.items():
        if hashlib.sha256((output / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"IMMUTABLE_ENDPOINT_INPUT_CHANGED:{name}")
    checkpoint_save()
    save_continuation_tables(output, data)
    print(json.dumps(data["summary"]), flush=True)
    return data


def save_continuation_tables(output, data):
    points = []
    for entry in data["points"].values():
        for row in entry["reconstruction"]:
            item = {k: row[k] for k in ("kappa_theta", "sorted_position", "Lambda", "frequency_hz", "shape_key", "reconstruction_status")}
            item["target_prefix_status"] = entry["prefix"]["target_prefix_status"]
            for key, source in (("abs_Delta_psi", "abs_Delta_psi_mass_normalized"), ("s", "s"),
                ("mass_error", "mass_normalization_error"), ("physical_residual", "max_physical_residual"), ("boundary_residual", "boundary_residual")):
                if source in row:
                    item[key] = row[source]
            points.append(item)
    for name, rows in (("points", points), ("mapping", data["mapping"]), ("trends", data["trends"])):
        fields = list(dict.fromkeys(k for r in rows for k in r))
        if fields:
            write_csv(output / f"kappa_continuation_{name}.csv", rows, fields)
    marker = "\n## Bounded kappa continuation\n"
    previous = (output / "report.md").read_bytes()
    for boundary in (marker.encode("utf-8"), marker.replace("\n", "\r\n").encode("utf-8")):
        previous = previous.split(boundary)[0]
    lines = [marker, "Command: `--mode kappa-continuation`. Original frequency/shape files are immutable SHA-256 inputs.",
        "DIRECT_1_TO_10 retains its four UNRESOLVED mappings. CONTINUATION_1_TO_10 is a separate result.",
        "Precedent: inplane_rotational_spring_eb_hinge_comparison.md and prepare_seed_mapping/continue_seed_mapping.",
        "Existing circular providers, public generic solver with default SearchSettings, target_prefix,",
        "reconstruct_saved_mode and modes.assign are reused. No new model/helper, symmetry classes or solver.",
        "Only failed finite intervals are subdivided (logarithmic depth <=2); exact RIGID permits only bridge 1000.",
        "MAC>=.95 and margin>=.20 are unchanged; assignment cost uses only physical mass vectors.",
        "The six local seeds are assigned together to seven candidates; failed terminal paths remain unresolved.",
        f"\nStatus: `{data['status']}`; summary: `{json.dumps(data['summary'])}`\n",
        "| kappa | Native spectrum | Target prefix | Reconstructed forms | Higher spectrum |", "|---|---|---|---|---|"]
    for entry in data["points"].values():
        p = entry["prefix"]
        lines.append(f"| {entry['kappa_theta']:.12g} | {p['solver_spectrum_status']} | {p['target_prefix_status']} | {sum(r['reconstruction_status']=='CONFIRMED' for r in entry['reconstruction'])}/7 | {p['higher_spectrum_status']} |")
    lines += ["\nNative full spectra, higher qualifications, reconstruction residuals and new shape arrays are in the separate continuation JSON.",
        "No old root or shape was recalculated; the previous high-spectrum qualification is unchanged.",
        "\n| Seed | Direct 1->10 | Continuation 1->10 | Path to 10 | Path to 100 | Path to RIGID | Endpoint check | Full path |",
        "|---|---|---|---|---|---|---|---|"]
    for r in data["trends"]:
        lines.append(f"| {r['branch_id']} | {r['direct_1_to_10_status']} | {r['continuation_1_to_10_status']} | {r.get('path_to_10','')} | {r.get('path_to_100','')} | {r.get('path_to_rigid','')} | {r.get('direct_endpoint_agreement','')} | {r['full_path_status']} |")
    lines += ["\n| Seed | Local step | Sorted positions | MAC | Margin | Status | Used in final path |",
        "|---|---|---|---|---|---|---|"]
    for r in data["mapping"]:
        lines.append(f"| {r['branch_id']} | {r['source_kappa']} -> {r['target_kappa']} | {r['source_sorted_position']} -> {r['target_sorted_position']} | {r['MAC']:.10g} | {r['margin']:.10g} | {r['step_status']} | {r['used_in_path']} |")
    lines += ["\n| Seed (M=1) | abs Delta psi(1) | (10) | (100) | RIGID | s(1) | s(10) | s(100) | Rotation | s trend |",
        "|---|---|---|---|---|---|---|---|---|---|"]
    for r in data["trends"]:
        if r["full_path_status"] == "FULL_KAPPA_PATH_CONFIRMED":
            numbers = " | ".join(f"{r[k]:.10g}" for k in ("abs_Delta_psi_k1", "abs_Delta_psi_k10", "abs_Delta_psi_k100", "abs_Delta_psi_rigid", "s_k1", "s_k10", "s_k100"))
            lines.append(f"| {r['branch_id']} | {numbers} | {r['rotation_trend_status']} | {r['sensitivity_trend_status']} |")
    confirmed = [r for r in data["trends"] if r["full_path_status"] == "FULL_KAPPA_PATH_CONFIRMED"]
    nonmonotone = [r["branch_id"] for r in confirmed if r["rotation_trend_status"] != "DECREASING"
                  or r["sensitivity_trend_status"] != "DECREASING"]
    lines += [f"\nConfirmed full paths: {len(confirmed)}/6. Nonmonotone endpoint sequences: {nonmonotone}.",
        "The table states actual endpoint changes; including intermediate points gives the counts recorded in the summary."]
    for r in confirmed:
        if r["branch_id"] in nonmonotone:
            lines.append(f"{r['branch_id']}: 1->10 changes abs Delta psi by {100*(r['abs_Delta_psi_k10']/r['abs_Delta_psi_k1']-1):.6g}% "
                f"and s by {100*(r['s_k10']/r['s_k1']-1):.6g}%. "
                "This is retained despite the confirmed MAC path. Sorted-frequency monotonicity does not impose monotonicity on either quantity.")
    if len(confirmed) == 6 and all(r["abs_Delta_psi_k100"] < r["abs_Delta_psi_k1"] for r in confirmed):
        lines.append("For all six tracked low modes of this EB example, the recorded stiffness path connects the finite-spring shapes to exact RIGID shapes; "
            "the rotation mismatch at 100 is smaller than at 1 and the exact RIGID mismatch is within the existing residual gate. "
            "Suppression is supported, but strict monotonic suppression for every mode is not.")
    lines += ["\nDelta psi and s never enter matching or acceptance. Nonmonotone sequences remain physical outputs.",
        "Full confirmed paths, including intermediate rotation and sensitivity values, are in descendant_path_records.",
        "s(kappa)=d ln(omega^2)/d kappa is local sensitivity; finite distance to RIGID accumulates changes",
        "over later stiffnesses. Different mode rankings by these quantities are not a contradiction.",
        "RIGID is exact, with s=NOT_APPLICABLE_RIGID. IDs are local diagnostic descendants, not canonical branches.",
        "This continuation establishes only the recorded path, not path independence in a parameter plane.",
        "No claim about real monolithic joints, their effective stiffness or an r/L-to-kappa relation follows."]
    if data.get("validation"):
        lines.append(f"\nValidation: `{json.dumps(data['validation'])}`")
    if data.get("stop_reason"):
        lines.append(f"\nStopped: `{data['stop_reason']}`. No dependent continuation is attempted.")
    (output / "report.md").write_bytes(previous + ("\n".join(lines) + "\n").encode("utf-8"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    parser.add_argument("--mode", choices=("spectrum", "shapes", "kappa-continuation"), default="spectrum")
    args = parser.parse_args()
    if args.mode == "kappa-continuation":
        result = continuation_run(args.output_dir)
        sys.exit(0 if result["status"] == "COMPLETED_BOUNDED_CONTINUATION" and result["summary"]["full_confirmed_count"] == 6 else 2)
    elif args.mode == "shapes":
        result = shapes_run(args.output_dir)
        sys.exit(0 if all(s == "CONFIRMED" for s in result["summary"]["state_status"].values()) else 2)
    else:
        result = run(args.output_dir)
        print(result["status"], flush=True)
        sys.exit(0 if result["status"] == "COMPLETED_TARGET_PREFIX" else 2)
