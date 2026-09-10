"""Fixed, diagnostic-only EB spring pilot; missing-only, sequential, no plots.

Uses the existing matrix determinant AND sigma-min detector/local refiner,
but never its full RLB inventory runner. See the tracked pilot report before use.
"""
from __future__ import annotations

import argparse
import copy
import csv
from dataclasses import asdict, replace
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

for _variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_variable] = "1"
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import scipy
from scipy.optimize import brentq

from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.analysis.laminated_beams import pilot_reddy_symmetric_coupled_beams_beta0 as roots
from my_project.analytic.formulas import assemble_clamped_coupled_matrix

ARM = eb.EBArm(A=.20*.05, D=.20*.05**3/12, m=.20*.05, L=1.)
FREQUENCY_SCALE = ARM.L**2 * math.sqrt(ARM.m/ARM.D)
KAPPAS = (0., .1, 1., 100., 10000.)
LIMITS = dict(sigma_ratio=1e-9, null_residual=1e-9, physical_residual=1e-9,
              compatibility=1e-10, rank_rtol=1e-12, spectrum_relative=1e-9,
              monotonic_relative=1e-9, root_xtol_Omega=1e-11,
              scan_step_Omega=4/64, window_Omega=4., overlap_Omega=.125,
              guard_margin_Omega=.02, max_Omega=240., group_seconds=120.)
ROW_CONTROL = np.array([2., .5, 3., 1/3, 5., .2])
OUTPUT = ROOT / "results/laminated_beams/inplane_rotational_spring_eb_pilot"
CODE_FILES = [Path(__file__).relative_to(ROOT).as_posix(),
              "scripts/lib/inplane_rotational_spring_eb.py",
              "scripts/lib/reddy_symmetric_coupled_beams.py",
              "scripts/lib/reddy_inplane_geometry.py",
              "scripts/analysis/laminated_beams/pilot_reddy_symmetric_coupled_beams_beta0.py",
              "src/my_project/analytic/formulas.py"]


def cases() -> list[dict]:
    return [dict(case_id=f"beta{beta}_{label}", beta_deg=beta, beta_rad=math.radians(beta),
                 mode=mode, kappa_theta=kappa,
                 k_theta=None if kappa is None else kappa*ARM.D/ARM.L)
            for beta in (0, 30)
            for mode, label, kappa in [*(('SPRING', f'k{k:g}', k) for k in KAPPAS),
                                      ('RIGID', 'RIGID', None)]]


def policy(left: float, right: float):
    # As in the existing bounded pilots: adapt a COPY, never shared defaults.
    result = copy.copy(roots.SearchPolicy())
    for name, value in dict(requested_roots=6, guard_roots=1, omega_bar_min=left,
                            omega_bar_max=right, post_guard_tail_bar=0.).items():
        object.__setattr__(result, name, value)
    return result


def scan(provider, case_id, left, right, *, repair=False):
    return roots._scan_candidates(
        provider, FREQUENCY_SCALE, policy(left, right), case_id=case_id,
        builder_id="physical_EB_rotational_spring", scan_id="LOCAL" if repair else "BASE",
        points=max(33, int(math.ceil((right-left)/LIMITS["scan_step_Omega"]))+1)*(2 if repair else 1),
        phases=(0., .5) if repair else (0.,))[0]


def candidate_record(candidate):
    return dict(Omega=candidate.omega_bar, accepted=candidate.accepted,
                reason=candidate.rejection_reason, sources=candidate.detection_sources,
                interval=[candidate.interval_left_bar, candidate.interval_right_bar],
                sigma_ratio=float(candidate.diagnostics.scaled_sigma_ratio),
                nullity=candidate.diagnostics.detected_nullity)


def consolidate(candidates):
    """Only reconcile duplicate *detections* at floating-point resolution.

    An overlap, matching strict/root-gate nullity and <=64 ULP separation are
    all required. Other close candidates remain unresolved, not one averaged
    root. Multiplicity is the matrix nullity, never the number of detections.
    """
    accepted, ambiguous = [], []
    for candidate in sorted((c for c in candidates if c.accepted), key=lambda c: c.omega_bar):
        if accepted:
            previous = accepted[-1]
            gap = candidate.omega_bar - previous.omega_bar
            near = 5e-10 + 5e-12*abs(candidate.omega_bar)
            if gap <= near:
                overlap = max(candidate.interval_left_bar, previous.interval_left_bar) <= min(
                    candidate.interval_right_bar, previous.interval_right_bar)
                counts = [c.diagnostics.detected_nullity for c in (candidate, previous)]
                gates = [c.diagnostics.root_gate_nullity for c in (candidate, previous)]
                if gap <= 64*abs(np.spacing(candidate.omega_bar)) and overlap and counts == gates and counts[0] == counts[1]:
                    best = min((candidate, previous), key=lambda c: c.diagnostics.scaled_sigma_ratio)
                    accepted[-1] = replace(best, detection_sources=tuple(sorted(set(
                        candidate.detection_sources + previous.detection_sources))))
                    continue
                ambiguous.extend([previous, candidate])
        accepted.append(candidate)
    return accepted, ambiguous


def suspicious(candidate):
    if candidate.accepted:
        return candidate.diagnostics.detected_nullity != candidate.diagnostics.root_gate_nullity
    return (candidate.rejection_reason != "FALSE_SIGMA_VALLEY" or
            candidate.diagnostics.scaled_sigma_ratio <= roots.SearchPolicy().sigma_prefilter)


def checked_endpoints(Omega, case):
    joint = eb.Joint(case["mode"], case["k_theta"])
    assembly = eb.boundary_assembly(Omega/FREQUENCY_SCALE, ARM, ARM, case["beta_rad"], joint, ARM)
    result = eb.endpoint_diagnostics(assembly, case["beta_rad"], joint, ARM, LIMITS["rank_rtol"])
    if result["nullity"] < 1 or result["sigma_ratio"] > LIMITS["sigma_ratio"]:
        raise RuntimeError("ENDPOINT_SINGULARITY_FAIL")
    for vector in result["vectors"]:
        residual = np.abs(vector["normalized_physical_residuals"])
        if (np.max(residual[:2]) > LIMITS["compatibility"] or
                np.max(residual) > LIMITS["physical_residual"] or
                max(vector["boundary_residual"], vector["scaled_residual"]) > LIMITS["null_residual"]):
            raise RuntimeError("PHYSICAL_ENDPOINT_RESIDUAL_FAIL")
        if case["kappa_theta"] == 0 and max(abs(x) for x in vector["moments_in_reference_units"]) > LIMITS["physical_residual"]:
            raise RuntimeError("HINGE_MOMENT_FAIL")
    return result


def solve_group(case, predictors=()):
    started = time.perf_counter()
    joint = eb.Joint(case["mode"], case["k_theta"])
    provider = lambda omega: eb.boundary_assembly(omega, ARM, ARM, case["beta_rad"], joint, ARM).dimensionless
    left, pool, windows, repairs = 1e-8, [], [], []
    while left < LIMITS["max_Omega"]:
        if time.perf_counter()-started > LIMITS["group_seconds"]:
            raise RuntimeError("GROUP_COST_LIMIT: completed groups remain saved")
        right = min(left+LIMITS["window_Omega"], LIMITS["max_Omega"])
        # Neighbours only position a window edge; all intervals below guard
        # still pass through both detectors, no predictor is an output root.
        for estimate in predictors:
            if left+1 < estimate+.25 < right:
                right = estimate+.25
                break
        found = scan(provider, case["case_id"], left, right)
        pool.extend(found)
        windows.append([left, right])
        events, ambiguous = consolidate(pool)
        provisional = [event for event in events for _ in range(event.diagnostics.detected_nullity)]
        cutoff = provisional[6].omega_bar if len(provisional) >= 7 else right
        suspects = [c for c in pool if c.omega_bar <= cutoff and suspicious(c)]
        suspects += [c for c in ambiguous if c.omega_bar <= cutoff]
        if suspects:
            # One bounded local recovery per affected interval; no global retry.
            for candidate in suspects:
                lo = max(1e-8, candidate.interval_left_bar-.01)
                hi = min(right, candidate.interval_right_bar+.01)
                if any(a <= candidate.omega_bar <= b for a, b in repairs):
                    continue
                local = scan(provider, case["case_id"], lo, hi, repair=True)
                pool = [c for c in pool if not lo < c.omega_bar < hi] + local
                repairs.append([lo, hi])
            events, ambiguous = consolidate(pool)
            provisional = [event for event in events for _ in range(event.diagnostics.detected_nullity)]
            cutoff = provisional[6].omega_bar if len(provisional) >= 7 else right
            if any(c.omega_bar <= cutoff for c in ambiguous) or any(
                    c.omega_bar <= cutoff and suspicious(c) for c in pool):
                raise RuntimeError("UNRESOLVED_CANDIDATE_AFTER_LOCAL_RECOVERY: " +
                                   json.dumps([candidate_record(c) for c in pool if c.omega_bar <= cutoff]))
        if len(provisional) >= 7:
            guard = provisional[6].omega_bar
            if len(provisional) > 7 and provisional[7].omega_bar == guard:
                raise RuntimeError("MULTIPLICITY_CROSSES_GUARD: guard event must not be truncated")
            if right-guard <= LIMITS["guard_margin_Omega"]:
                raise RuntimeError("GUARD_AT_WINDOW_EDGE: no automatic spectral tail")
            rows, endpoints = [], []
            event_index = {id(event): index for index, event in enumerate(events, 1)}
            repeats = {}
            for slot, event in enumerate(provisional[:7], 1):
                diagnostic = checked_endpoints(event.omega_bar, case)
                if diagnostic["nullity"] != event.diagnostics.detected_nullity:
                    raise RuntimeError("MULTIPLICITY_DIAGNOSTIC_MISMATCH")
                repeats[id(event)] = repeats.get(id(event), 0)+1
                rows.append(dict(**case, group_role="BASE", sorted_position=slot,
                                 role="GUARD" if slot == 7 else "ROOT",
                                 event_id=f"event{event_index[id(event)]}", repeated_root_slot=repeats[id(event)],
                                 multiplicity=diagnostic["nullity"], Omega=event.omega_bar,
                                 omega=event.omega_bar/FREQUENCY_SCALE, Lambda=math.sqrt(event.omega_bar)))
                endpoints.append(diagnostic)
            return dict(status="COMPLETED", rows=rows, endpoints=endpoints, windows=windows,
                        local_recoveries=repairs, candidates=[candidate_record(c) for c in pool if c.omega_bar <= guard],
                        guard_gap_Omega=right-guard, unresolved_below_guard=0,
                        seconds=time.perf_counter()-started)
        left = right-LIMITS["overlap_Omega"]
        if right == LIMITS["max_Omega"]:
            break
    raise RuntimeError("GUARD_NOT_FOUND_WITHIN_BOUNDED_RANGE")


def local_control(provider, group, case):
    """Only seven saved-root neighbourhoods, counted separately from BASE."""
    differences, physical_maxima = [], []
    for row in group["rows"]:
        center = row["Omega"]
        found = scan(provider, "CONTROL", center-1e-4, center+1e-4, repair=True)
        events, ambiguous = consolidate(found)
        if len(events) != 1 or ambiguous or any(suspicious(c) for c in found):
            raise RuntimeError("LOCAL_CONTROL_UNRESOLVED")
        if events[0].diagnostics.detected_nullity != row["multiplicity"]:
            raise RuntimeError("LOCAL_CONTROL_MULTIPLICITY_FAIL")
        differences.append(abs(events[0].omega_bar-center)/center)
        endpoint = checked_endpoints(events[0].omega_bar, case)
        physical_maxima.append(max(abs(value) for vector in endpoint["vectors"]
                                   for value in vector["normalized_physical_residuals"]))
    return dict(relative_differences=differences,
                normalized_physical_maxima=physical_maxima, local_root_checks=len(differences),
                status="WITHIN_TOLERANCE" if max(differences) <= LIMITS["spectrum_relative"] else "MISMATCH")


def axial_reference(guard):
    """Exact axial family, restricted to the saved group's guard."""
    axial = []
    n = 1
    while True:
        value = n*math.pi/(2*ARM.L)*math.sqrt(ARM.A/ARM.m)*FREQUENCY_SCALE
        if value > guard*(1+LIMITS["spectrum_relative"]):
            break
        axial.append(value)
        n += 1
    return axial


def straight_reference(guard):
    """Independent single-beam scalar equations; no roots above guard."""
    values, axial = [], axial_reference(guard)
    zmax = 2*math.sqrt(guard*(1+LIMITS["spectrum_relative"]))
    n = 1
    while (n+.25)*math.pi < zmax:
        lo, hi = (n+.25)*math.pi, min((n+.75)*math.pi, zmax)
        function = lambda z: math.cos(z)-1/math.cosh(z)
        if function(lo)*function(hi) <= 0:
            z = brentq(function, lo, hi, xtol=1e-13)
            values.append(z*z/4)
        n += 1
    return sorted(values+axial), axial


def controls_for_group(case, group):
    controls = {}
    spectrum = np.array([row["Omega"] for row in group["rows"]])
    if case["beta_deg"] == 0:
        axial = axial_reference(spectrum[-1])
        errors = [float(np.min(np.abs(spectrum-value))/value) for value in axial]
        controls["axial_in_range"] = dict(reference_Omega=axial, relative_differences=errors,
            status=("NOT_IN_RANGE" if not errors else "WITHIN_TOLERANCE" if max(errors) <= LIMITS["spectrum_relative"] else "MISMATCH"))
        if case["mode"] == "RIGID":
            reference, _ = straight_reference(spectrum[-1])
            errors = [abs(a-b)/b for a, b in zip(spectrum, reference)]
            controls["independent_straight_rigid"] = dict(reference_Omega=reference, relative_differences=errors,
                status="WITHIN_TOLERANCE" if len(reference) == 7 and max(errors) <= LIMITS["spectrum_relative"] else "MISMATCH")
    if case["mode"] == "RIGID":
        # The rectangle has eps=h/(sqrt(12)*l); no circular BeamParams/defaults.
        provider = lambda omega: assemble_clamped_coupled_matrix(
            math.sqrt(omega*FREQUENCY_SCALE), case["beta_rad"], 0., .05/math.sqrt(12))
        controls["legacy_rigid_regression"] = local_control(provider, group, case)
    if case["case_id"] == "beta30_k1":
        joint = eb.Joint(case["mode"], case["k_theta"])
        provider = lambda omega: ROW_CONTROL[:, None]*eb.boundary_assembly(
            omega, ARM, ARM, case["beta_rad"], joint, ARM).dimensionless
        controls["positive_row_scaling"] = local_control(provider, group, case)
    return controls


def cross_case_controls(groups):
    result = {}
    for beta in (0, 30):
        selected = [c for c in cases() if c["beta_deg"] == beta]
        if any(c["case_id"] not in groups for c in selected):
            result[str(beta)] = dict(status="PENDING_GROUPS")
            continue
        spectra = np.array([[r["Omega"] for r in groups[c["case_id"]]["rows"][:6]] for c in selected])
        increments = np.diff(spectra, axis=0)/spectra[:-1]
        errors = np.abs(spectra[-3:-1]-spectra[-1])/spectra[-1]
        result[str(beta)] = dict(
            relative_increments=increments.tolist(), finite_to_rigid_relative=errors.tolist(),
            monotonic=bool(np.min(increments) >= -LIMITS["monotonic_relative"]),
            approaches_rigid=bool(np.all(errors[1] <= errors[0]+LIMITS["spectrum_relative"])),
            status="WITHIN_TOLERANCE" if np.min(increments) >= -LIMITS["monotonic_relative"] and
            np.all(errors[1] <= errors[0]+LIMITS["spectrum_relative"]) else "MISMATCH")
    return result


def atomic_text(path, content):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(content, encoding="utf-8")
    os.replace(temporary, path)


def json_text(value):
    return json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmark-only", action="store_true", help="only beta=30 deg, kappa=1; saves this BASE group")
    args = parser.parse_args()
    started = time.perf_counter()
    versions = dict(python=sys.version, executable=sys.executable, numpy=np.__version__,
                    scipy=scipy.__version__, platform=platform.platform())
    hashes = {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in CODE_FILES}
    contract = dict(schema="eb-rotational-pilot-v1", policy="frequency-map-v1/fast_plot",
                    spectrum_semantics="sorted_positions", geometry=dict(l=1, b=.20, h=.05, E=1, rho=1),
                    arm=asdict(ARM), normalization="Omega=omega*l^2*sqrt(m/D); Lambda=sqrt(Omega)",
                    cases=cases(), requested_roots=6, guard_roots=1, limits=LIMITS,
                    search_policy_seed=asdict(roots.SearchPolicy()), code_sha256=hashes,
                    dependency_versions=versions)
    contract = json.loads(json_text(contract))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    path = OUTPUT/"diagnostics.json"
    if not path.exists() and any((OUTPUT/name).exists() for name in ("spectrum_roots.csv", "run_manifest.json")):
        raise SystemExit("Orphan output without authoritative diagnostics; no existing data overwritten")
    state = json.loads(path.read_text(encoding="utf-8")) if path.exists() else dict(contract=contract, groups={}, attempts=[])
    if state["contract"] != contract:
        raise SystemExit("Existing output contract/version differs; no data overwritten")
    groups = state["groups"]
    reused, new, controls_completed = list(groups), [], []
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    git_status = subprocess.check_output(["git", "status", "--short"], cwd=ROOT, text=True)
    manifest = dict(contract=contract, source_HEAD=head, source_working_tree_status=git_status,
                    run_started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                    reused_groups=reused, new_groups=new, controls_completed=controls_completed, failure=None)

    def save():
        state["cross_case_controls"] = cross_case_controls(groups)
        atomic_text(path, json_text(state))  # authoritative completed-group checkpoint
        rows = [row for case in cases() if case["case_id"] in groups for row in groups[case["case_id"]]["rows"]]
        if rows:
            output = io.StringIO(newline="")
            writer = csv.DictWriter(output, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
            atomic_text(OUTPUT/"spectrum_roots.csv", output.getvalue())
        manifest.update(completed_BASE_groups=len(groups), BASE_rows=len(rows),
                        elapsed_seconds=time.perf_counter()-started,
                        group_seconds={key: value["seconds"] for key, value in groups.items()},
                        controls_seconds={key: value.get("controls_seconds") for key, value in groups.items()},
                        local_recoveries={key: value["local_recoveries"] for key, value in groups.items()})
        checks = [c for group in groups.values() for c in group.get("controls", {}).values()]
        checks += list(state["cross_case_controls"].values())
        manifest["status"] = (
            "STOPPED" if manifest["failure"] else
            "CONTROL_MISMATCH" if any(c.get("status") == "MISMATCH" for c in checks) else
            "COMPLETED_FINITE_CHECKS" if len(groups) == 12 and all(g.get("controls_done") for g in groups.values()) else
            "PARTIAL")
        atomic_text(OUTPUT/"run_manifest.json", json_text(manifest))

    ordered = sorted(cases(), key=lambda c: c["case_id"] != "beta30_k1")
    if args.benchmark_only:
        ordered = ordered[:1]
    save()  # repairs CSV/manifest after an interruption, without recalculation
    for case in ordered:
        if case["case_id"] in groups and groups[case["case_id"]].get("controls_done"):
            if any(c.get("status") == "MISMATCH" for c in groups[case["case_id"]]["controls"].values()):
                raise SystemExit("Saved control mismatch requires local diagnosis; no BASE rerun")
            continue
        try:
            if case["case_id"] not in groups:
                neighbours = [g for key, g in groups.items() if key.startswith(f"beta{case['beta_deg']}_")]
                predictors = [r["Omega"] for r in neighbours[-1]["rows"]] if neighbours else []
                groups[case["case_id"]] = solve_group(case, predictors)
                new.append(case["case_id"])
                save()  # BASE is durable even if a subsequent control fails
            group = groups[case["case_id"]]
            control_started = time.perf_counter()
            group["controls"] = controls_for_group(case, group)
            group["controls_seconds"] = time.perf_counter()-control_started
            group["controls_done"] = True
            controls_completed.append(case["case_id"])
            save()
            print(f"{case['case_id']}: 7 BASE rows, {group['seconds']:.2f} s", flush=True)
            controls = list(group["controls"].values())+list(state["cross_case_controls"].values())
            if any(c.get("status") == "MISMATCH" for c in controls):
                raise RuntimeError("CONTROL_MISMATCH: saved result retained for local diagnosis")
        except (RuntimeError, ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
            manifest["failure"] = dict(case_id=case["case_id"], reason=str(error))
            state["attempts"].append(manifest["failure"])
            save()
            raise SystemExit(str(error)) from error
    save()
    if any(c.get("status") == "MISMATCH" for c in state["cross_case_controls"].values()):
        raise SystemExit("Saved cross-case mismatch requires local diagnosis; no BASE rerun")


if __name__ == "__main__":
    main()
