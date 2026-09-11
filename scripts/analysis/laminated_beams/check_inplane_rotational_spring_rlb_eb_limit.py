"""Bounded, sequential RLB -> EB check at fixed finite kappa_theta=1.

New I/O contract: reduced plies + epsilon_limit, RLB endpoints, reused EB
references. This is not an EB runner preset. No full inventory or legacy
controls are called. --matrices-only has no root search or root evaluation.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import csv
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
from scripts.lib import inplane_rotational_spring_rlb as rlb
from scripts.lib import inplane_rotational_spring_eb as eb
from scripts.lib import reddy_symmetric_laminated_beam as native
from scripts.analysis.laminated_beams import pilot_inplane_rotational_spring_eb as workflow

REF = workflow.ARM
FS = workflow.FREQUENCY_SCALE
SPRING = eb.Joint("SPRING", REF.D/REF.L)
OUTPUT = ROOT / "results/laminated_beams/inplane_rotational_spring_rlb_eb_limit"
EB_OUTPUT = ROOT / "results/laminated_beams/inplane_rotational_spring_eb_pilot"
# Frozen before evaluation; matrix tolerances apply in common fixed units.
LIMITS = dict(matrix_atol=1e-12, constitutive_rtol=1e-12, H_rtol=1e-12,
              transfer_boundary_rtol=1e-9, frequency_interpretation_rtol=1e-6,
              sigma_ratio=1e-9, rank_rtol=1e-12, physical_residual=1e-9,
              compatibility=1e-10, null_residual=1e-9,
              max_builds_per_group=6000, max_local_recoveries_total=2,
              scan_step_Omega=4/64, window_Omega=4., overlap_Omega=.125,
              guard_margin_Omega=.02, max_Omega=114.)
CODE_FILES = [Path(__file__).relative_to(ROOT).as_posix(),
              "scripts/lib/inplane_rotational_spring_rlb.py",
              "scripts/lib/reddy_symmetric_laminated_beam.py", *workflow.CODE_FILES]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cases():
    return [dict(case_id=f"beta{beta}_eps{epsilon:g}", beta_deg=beta,
                 beta_rad=math.radians(beta), epsilon_limit=epsilon,
                 kappa_theta=1., k_theta=SPRING.k_theta)
            for beta in (0, 30) for epsilon in (.01, .1, 1.)]


def matrix_checks(properties, section):
    started = time.perf_counter()
    checks = []
    def compare(label, actual, expected, rtol, **context):
        error = float(np.linalg.norm(np.asarray(actual)-expected))
        norm = float(np.linalg.norm(expected))
        checks.append(dict(label=label, **context, absolute_error=error,
                           relative_error=error/norm if norm else None,
                           reference_norm=norm, tolerance=LIMITS["matrix_atol"]+rtol*norm,
                           passed=error <= LIMITS["matrix_atol"]+rtol*norm))
    expected = [REF.A, REF.D, REF.m, 5/6/2.6*.20*.05, .20*.05**3/12]
    for name, value in zip(("A", "D", "m", "S", "J"), expected):
        # Each scalar is normalized independently; absolute term is dimensionless.
        compare(name, getattr(properties, name)/value, 1., LIMITS["constitutive_rtol"])
    compare("B_symmetry", section.B/(np.linalg.norm(section.A)*.05), np.zeros((3,3)), 0.)
    compare("I1_symmetry", section.I1/(section.I0*.05), 0., 0.)
    z = eb.state_scale(REF)
    common = lambda matrix: matrix*z[None, :]/z[:, None]
    for beta in (0, 30):
        for Omega in (2., 20., 80.):
            omega = Omega/FS
            one, zero = rlb.LimitArm(properties, 1., 1.), rlb.LimitArm(properties, 1., 0.)
            context = dict(beta_deg=beta, Omega=Omega)
            for label, left, right, rtol in (
                ("H_native_eps1", rlb.state_matrix(omega, one), native.combined_state_matrix(omega, properties), LIMITS["H_rtol"]),
                ("T_native_eps1", rlb.transfer_matrix(omega, one), native.combined_transfer_matrix(omega, 1., properties), LIMITS["transfer_boundary_rtol"]),
                ("H_EB_eps0", rlb.state_matrix(omega, zero), eb.state_matrix(omega, REF), LIMITS["H_rtol"]),
                ("T_EB_eps0", rlb.transfer_matrix(omega, zero), eb.transfer_matrix(omega, REF), LIMITS["transfer_boundary_rtol"]),
            ):
                compare(label, common(left), common(right), rtol, **context)
            boundary = rlb.boundary_assembly(omega, zero, zero, math.radians(beta), SPRING, REF)
            reference = eb.boundary_assembly(omega, REF, REF, math.radians(beta), SPRING, REF)
            # Fixed reaction/row scales, BEFORE adaptive equilibration.
            left = boundary.row_factors[:, None]*boundary.physical*reference.reaction_scales[None, :]
            compare("B_EB_eps0", left, reference.dimensionless, LIMITS["transfer_boundary_rtol"], **context)
    law = eb.joint_matrix(0., SPRING)[2]
    expected_law = np.zeros(12)
    expected_law[[2,5,8]] = [REF.D, 1., -REF.D]
    compare("finite_spring_row", law, expected_law, 0.)
    return dict(status="MATRICES_CONFIRMED" if all(c["passed"] for c in checks) else "MATRIX_MISMATCH",
                norm="Frobenius (scalars: absolute); fixed common dimensionless units",
                checks=checks, seconds=time.perf_counter()-started,
                properties={name: getattr(properties, name) for name in ("A","D","S","m","J","K","width")},
                laminate={name: np.asarray(getattr(section,name)).tolist() for name in ("z_interfaces","A","B","D","shear","I0","I1","I2")})


def read_eb_references():
    paths = [EB_OUTPUT/name for name in ("spectrum_roots.csv", "diagnostics.json", "run_manifest.json")]
    data, manifest = [json.loads(p.read_text(encoding="utf-8")) for p in paths[1:]]
    contract = data["contract"]
    if contract["geometry"] != dict(l=1,b=.20,h=.05,E=1,rho=1):
        raise ValueError("EB reference geometry mismatch")
    if contract["normalization"] != "Omega=omega*l^2*sqrt(m/D); Lambda=sqrt(Omega)":
        raise ValueError("EB reference normalization mismatch")
    if not np.allclose(list(contract["arm"].values()), list(asdict(REF).values()), rtol=1e-12, atol=0):
        raise ValueError("EB reference arm mismatch")
    with paths[0].open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    references = {}
    for beta in (0,30):
        name = f"beta{beta}_k1"
        group = data["groups"][name]
        if group["status"] != "COMPLETED" or group["unresolved_below_guard"] != 0:
            raise ValueError("EB reference group is not accepted")
        selected = [row for row in rows if row["case_id"] == name]
        if len(selected) != len(group["rows"]) or len(selected) < 7:
            raise ValueError("EB CSV/diagnostics count mismatch")
        for row, saved in zip(selected, group["rows"]):
            if (float(row["Omega"]) != saved["Omega"] or float(row["k_theta"]) != SPRING.k_theta or
                    float(row["kappa_theta"]) != 1 or row["mode"] != "SPRING" or
                    float(row["beta_deg"]) != beta or
                    not math.isclose(float(row["Omega"]), float(row["omega"])*FS, rel_tol=1e-14)):
                raise ValueError("EB reference row mismatch")
        references[str(beta)] = dict(source_case=name, origin="REUSED_EB_REFERENCE",
            status=group["status"], rows=group["rows"], guard_gap_Omega=group["guard_gap_Omega"],
            source_code_sha256=group["source_code_sha256"])
    return references, dict(files_sha256={p.relative_to(ROOT).as_posix(): digest(p) for p in paths},
        source_HEAD=manifest["source_HEAD"], source_working_tree_status=manifest["source_working_tree_status"],
        original_contract=contract, original_overall_status=manifest["status"],
        original_unaccepted_BASE_cases=manifest["unaccepted_BASE_cases"],
        original_unresolved_control=manifest["unresolved_control"])


class CostLimit(RuntimeError):
    pass


class Provider:
    """Count actual RLB boundary builds, including endpoint diagnostics.

Exact float-key reuse is per group; cache hits do not construct a matrix.
The cached endpoint map is always the RLB one, including epsilon_limit=0.
"""
    def __init__(self, properties, case):
        self.arm = rlb.LimitArm(properties, 1., case["epsilon_limit"])
        self.case, self.builds, self.cache = case, 0, {}

    def assembly(self, omega):
        key = float(omega)
        if key not in self.cache:
            if self.builds >= LIMITS["max_builds_per_group"]:
                raise CostLimit("COST_LIMIT")
            self.builds += 1
            self.cache[key] = rlb.boundary_assembly(key, self.arm, self.arm, self.case["beta_rad"], SPRING, REF)
        return self.cache[key]

    def __call__(self, omega):
        return self.assembly(omega).dimensionless


def endpoints(provider, Omega):
    result = eb.endpoint_diagnostics(provider.assembly(Omega/FS), provider.case["beta_rad"], SPRING, REF, LIMITS["rank_rtol"])
    failures = []
    if result["nullity"] < 1 or result["sigma_ratio"] > LIMITS["sigma_ratio"]:
        failures.append("ENDPOINT_SINGULARITY_FAIL")
    for vector in result["vectors"]:
        residual = np.abs(vector["normalized_physical_residuals"])
        if max(residual[:2]) > LIMITS["compatibility"]:
            failures.append("COMPATIBILITY_FAIL")
        if max(residual) > LIMITS["physical_residual"]:
            failures.append("PHYSICAL_RESIDUAL_FAIL")
        if max(vector["boundary_residual"],vector["scaled_residual"]) > LIMITS["null_residual"]:
            failures.append("NULL_RESIDUAL_FAIL")
    result.update(Omega=Omega, failures=sorted(set(failures)))
    return result


def exact_limit(properties, beta, reference):
    started = time.perf_counter()
    case = dict(case_id=f"beta{beta}_eps0", beta_deg=beta, beta_rad=math.radians(beta), epsilon_limit=0.)
    provider = Provider(properties, case)
    diagnostics = [endpoints(provider, row["Omega"]) for row in reference["rows"]]
    return dict(origin="EXACT_LIMIT_EVALUATION", frequencies_origin="REUSED_EB_REFERENCE",
                status="CONFIRMED_AT_REUSED_ROOTS" if not any(d["failures"] for d in diagnostics) else "QUALIFIED",
                endpoints=diagnostics, boundary_builds=provider.builds, seconds=time.perf_counter()-started)


def solve_group(properties, case, predictor, recovery_budget):
    started = time.perf_counter()
    provider = Provider(properties, case)
    pool, windows, reconciliations, recoveries = [], [], [], []
    result = dict(case=case, origin="NEW_RLB_SPECTRUM", status="INCOMPLETE", rows=[], endpoints=[])
    left = 1e-8
    try:
        while left < LIMITS["max_Omega"]:
            right = min(left+LIMITS["window_Omega"], LIMITS["max_Omega"])
            for estimate in predictor:
                if left+1 < estimate+.25 < right:
                    right = estimate+.25
                    break
            found = workflow.scan(provider, case["case_id"], left, right)
            pool.extend(found)
            windows.append([left,right])
            # Existing finite isolation diagnostics, not a new search grid:
            # intersecting brackets + one monotone sign change + second-sigma gap.
            pool, evidence = workflow.reconcile_local_detections(pool, provider)
            reconciliations.extend(evidence)
            events, ambiguous = workflow.consolidate(pool)
            slots = [event for event in events for _ in range(event.diagnostics.detected_nullity)]
            cutoff = slots[6].omega_bar if len(slots) >= 7 else right
            suspects = [c for c in pool if c.omega_bar <= cutoff and workflow.suspicious(c)]
            suspects += [c for c in ambiguous if c.omega_bar <= cutoff]
            if suspects and recovery_budget[0] < LIMITS["max_local_recoveries_total"]:
                candidate = suspects[0]
                lo, hi = max(1e-8,candidate.interval_left_bar-.01), min(right,candidate.interval_right_bar+.01)
                if not any(a <= candidate.omega_bar <= b for a,b in recoveries):
                    recovery_budget[0] += 1
                    recoveries.append([lo,hi])
                    local = workflow.scan(provider, case["case_id"], lo, hi, repair=True)
                    pool = [c for c in pool if not lo < c.omega_bar < hi]+local
                    pool, evidence = workflow.reconcile_local_detections(pool, provider)
                    reconciliations.extend(evidence)
                    events, ambiguous = workflow.consolidate(pool)
                    slots = [event for event in events for _ in range(event.diagnostics.detected_nullity)]
            if len(slots) >= 7:
                guard = slots[6].omega_bar
                # Never cut a matrix-nullity multiplicity at the seventh slot.
                kept = [c for c in slots if c.omega_bar <= guard]
                suspects = [c for c in pool if c.omega_bar <= guard and workflow.suspicious(c)]
                suspects += [c for c in ambiguous if c.omega_bar <= guard]
                sixth = slots[5].omega_bar
                unresolved_target = [c for c in suspects if c.interval_left_bar <= sixth+LIMITS["guard_margin_Omega"]]
                result.update(guard_gap_Omega=right-guard, target_guard_gap_Omega=guard-sixth,
                              unresolved=[workflow.candidate_record(c) for c in suspects])
                for position, event in enumerate(kept,1):
                    diagnostic = endpoints(provider,event.omega_bar)
                    if diagnostic["nullity"] != event.diagnostics.detected_nullity:
                        diagnostic["failures"].append("MULTIPLICITY_MISMATCH")
                    result["endpoints"].append(diagnostic)
                    result["rows"].append(dict(**case, origin="NEW_RLB_SPECTRUM", sorted_position=position,
                        role="ROOT" if position <= 6 else "GUARD", multiplicity=event.diagnostics.detected_nullity,
                        Omega=event.omega_bar, omega=event.omega_bar/FS, Lambda=math.sqrt(event.omega_bar),
                        accepted_endpoint=not diagnostic["failures"]))
                target_fail = any(d["failures"] for d in result["endpoints"][:6])
                guard_fail = any(d["failures"] for d in result["endpoints"][6:])
                if unresolved_target or target_fail:
                    result["status"] = "TARGET_QUALIFIED"
                elif result["guard_gap_Omega"] <= LIMITS["guard_margin_Omega"] or guard == sixth:
                    result["status"] = "GUARD_NOT_SEPARATED"
                elif suspects or guard_fail:
                    result["status"] = "TARGET_CONFIRMED_GUARD_QUALIFIED"
                else:
                    result["status"] = "COMPLETED"
                break
            if len(slots) == 6:
                unresolved_above = [c for c in pool if workflow.suspicious(c) and
                                    c.interval_left_bar > slots[5].omega_bar+LIMITS["guard_margin_Omega"]]
                if unresolved_above:
                    # Do not continue scanning a tail just because a detected
                    # guard failed the strict gate. Preserve it for evaluation.
                    result["status"] = "INCOMPLETE"
                    break
            left = right-LIMITS["overlap_Omega"]
            if right == LIMITS["max_Omega"]:
                break
    except CostLimit:
        result["status"] = "COST_LIMIT"
    except (RuntimeError, ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
        result.update(status="NUMERICAL_FAILURE", failure=str(error))
    finally:
        result.update(boundary_builds=provider.builds, seconds=time.perf_counter()-started,
            windows=windows, local_recoveries=recoveries, detector_reconciliations=reconciliations,
            candidates=[workflow.candidate_record(c) for c in pool])
    return result


def evaluate_saved_target(properties, group, reference):
    """Evaluate saved frequencies only; no detector, refiner or root search.

This narrow completion records a guard qualification, never repairs its
frequency or discards its original rejection. Intended for six accepted
events followed by an isolated, rejected seventh event.
"""
    started = time.perf_counter()
    provider = Provider(properties, group["case"])
    # The provider limit covers original search AND this endpoint-only pass.
    original_builds = group["boundary_builds"]
    provider.builds = original_builds
    candidates = []
    for saved in group["candidates"]:
        if not saved["accepted"]:
            continue
        diag = workflow.roots.boundary_matrix_diagnostics(saved["Omega"], provider, FS,
            rank_relative_tolerance=LIMITS["rank_rtol"],root_ratio_tolerance=LIMITS["sigma_ratio"])
        accepted, reason = workflow.roots._candidate_quality(diag, workflow.policy(*saved["interval"]))
        candidates.append(workflow.roots.RootCandidate(group["case"]["case_id"],"RLB_saved_endpoint_evaluation",
            "NO_ROOT_SEARCH",saved["Omega"],tuple(saved["sources"]),*saved["interval"],True,diag,accepted,reason))
    events, ambiguous = workflow.consolidate(candidates)
    slots = [event for event in events for _ in range(event.diagnostics.detected_nullity)]
    if ambiguous or len(slots) != 6 or any(not c.accepted for c in candidates):
        raise ValueError("Saved target is not six unambiguous accepted positions")
    sixth = slots[-1].omega_bar
    rejected = [c for c in group["candidates"] if not c["accepted"] and
                (c["reason"] != "FALSE_SIGMA_VALLEY" or c["sigma_ratio"] <= workflow.roots.SearchPolicy().sigma_prefilter)]
    if len(rejected) != 1 or rejected[0]["interval"][0] <= sixth+LIMITS["guard_margin_Omega"]:
        raise ValueError("Unresolved event is not isolated above target range")
    guard = rejected[0]
    endpoint_data = [endpoints(provider,c.omega_bar) for c in slots]+[endpoints(provider,guard["Omega"])]
    if any(d["failures"] for d in endpoint_data[:6]):
        raise ValueError("Saved target physical gate failed")
    rows = []
    for position, diagnostic in enumerate(endpoint_data,1):
        Omega = diagnostic["Omega"]
        rows.append(dict(**group["case"],origin="NEW_RLB_SPECTRUM",sorted_position=position,
            role="ROOT" if position<=6 else "GUARD_CANDIDATE",multiplicity=diagnostic["nullity"],
            Omega=Omega,omega=Omega/FS,Lambda=math.sqrt(Omega),accepted_endpoint=not diagnostic["failures"]))
    group.update(original_search_status=group["status"],status="TARGET_CONFIRMED_GUARD_QUALIFIED",
        rows=rows,endpoints=endpoint_data,unresolved=[guard],
        target_guard_gap_Omega=guard["interval"][0]-sixth,
        guard_gap_Omega=group["windows"][-1][1]-guard["Omega"],
        qualification="Six target positions checked; rejected guard retains original nullity failure; no root refinement",
        saved_endpoint_evaluation=dict(boundary_builds=provider.builds-original_builds,
            seconds=time.perf_counter()-started,original_guard=guard),boundary_builds=provider.builds)
    group["seconds"] += group["saved_endpoint_evaluation"]["seconds"]
    fixed = [r["Omega"] for r in reference["rows"][:6]]
    group["e_j"] = [abs(r["Omega"]-v)/v for r,v in zip(rows[:6],fixed)]
    group["e_max"] = max(group["e_j"])


def atomic_text(path, value):
    temporary = path.with_suffix(path.suffix+".tmp")
    temporary.write_text(value, encoding="utf-8", newline="\n")
    temporary.replace(path)


def save(state):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    atomic_text(OUTPUT/"diagnostics.json",json.dumps(state,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    rows = []
    for beta, reference in state.get("eb_references",{}).items():
        for row in reference["rows"]:
            rows.append(dict(beta_deg=int(beta), epsilon_limit=0., origin="REUSED_EB_REFERENCE",
                check_origin="EXACT_LIMIT_EVALUATION", sorted_position=row["sorted_position"], role=row["role"],
                multiplicity=row["multiplicity"], Omega=row["Omega"],omega=row["omega"],Lambda=row["Lambda"],
                kappa_theta=1.,k_theta=SPRING.k_theta,
                status=state.get("exact_limit",{}).get(beta,{}).get("status","NOT_EVALUATED"),e_j=None))
    for group in state["groups"].values():
        for row in group["rows"]:
            pos = row["sorted_position"]
            rows.append(dict(beta_deg=row["beta_deg"], epsilon_limit=row["epsilon_limit"], origin=row["origin"],
                check_origin="NEW_RLB_SPECTRUM", sorted_position=pos,role=row["role"],multiplicity=row["multiplicity"],
                Omega=row["Omega"],omega=row["omega"],Lambda=row["Lambda"],kappa_theta=1.,k_theta=SPRING.k_theta,
                status=group["status"], e_j=group.get("e_j",[None]*6)[pos-1] if pos<=6 else None))
    if rows:
        stream = io.StringIO(newline="")
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        atomic_text(OUTPUT/"spectrum_comparison.csv",stream.getvalue())
    manifest = {key:value for key,value in state.items() if key not in ("groups","matrices","exact_limit","eb_references")}
    manifest.update(matrix_status=state["matrices"]["status"],
        groups={key:{k:g.get(k) for k in ("status","boundary_builds","seconds","local_recoveries","e_max")} for key,g in state["groups"].items()},
        exact_limit={key:{k:g[k] for k in ("status","boundary_builds","seconds")} for key,g in state.get("exact_limit",{}).items()},
        matrix_seconds=state["matrices"]["seconds"],
        total_measured_compute_seconds=state["matrices"]["seconds"]+sum(g["seconds"] for g in state["groups"].values())+sum(g["seconds"] for g in state.get("exact_limit",{}).values()))
    atomic_text(OUTPUT/"run_manifest.json",json.dumps(manifest,ensure_ascii=False,indent=2,allow_nan=False)+"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrices-only",action="store_true")
    parser.add_argument("--case",choices=[c["case_id"] for c in cases()],help="one declared group; default: remaining groups sequentially")
    parser.add_argument("--evaluate-saved-target",choices=[c["case_id"] for c in cases()],
                        help="endpoint-only evaluation of six saved events with an isolated rejected guard; no root search")
    args = parser.parse_args()
    section, properties = rlb.benchmark_section()
    code = {name:digest(ROOT/name) for name in dict.fromkeys(CODE_FILES)}
    path = OUTPUT/"diagnostics.json"
    if path.exists():
        state = json.loads(path.read_text(encoding="utf-8"))
        known_code = state.get("endpoint_evaluation_code_sha256",state["code_sha256"])
        if known_code != code:
            runner = Path(__file__).relative_to(ROOT).as_posix()
            only_runner_changed = all(value == code[name] for name,value in known_code.items() if name != runner)
            if not args.evaluate_saved_target or not only_runner_changed:
                raise ValueError("Code changed since saved run; inspect provenance before reuse")
            # Explicit endpoint-only phase keeps original search code hashes.
            state["endpoint_evaluation_code_sha256"] = code
            state["code_change_note"] = "Added saved-target endpoint evaluation and stop at rejected guard; existing spectra not rerun"
    else:
        state = dict(schema="rlb-eb-spring-limit-v1",policy="frequency-map-v1/fast_plot",
            spectrum_semantics="sorted_positions", requested_roots=6,guard_roots=1,
            geometry=dict(l=1,b=.2,h=.05,E=1,rho=1,nu=.3,ply_count=4,K=5/6),
            normalization="Omega=omega*l^2*sqrt(rho*Ag/(E*Ig)); Lambda=sqrt(Omega)",
            coefficient_path="invS=epsilon_limit/S_star; J=epsilon_limit*J_star; fixed A,D,m,L,beta,k_theta",
            cases=cases(),limits=LIMITS, code_sha256=code, source_HEAD=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
            source_working_tree_status=subprocess.check_output(["git","status","--short"],cwd=ROOT,text=True),
            executable=sys.executable, versions=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
            started_utc=datetime.now(timezone.utc).isoformat(),groups={},local_recovery_attempts=0,
            matrices=matrix_checks(properties,section))
        save(state)
    print("matrices:",state["matrices"]["status"],flush=True)
    if args.matrices_only or state["matrices"]["status"] != "MATRICES_CONFIRMED":
        return
    references, origin = read_eb_references()
    if "eb_origin" in state and state["eb_origin"] != origin:
        raise ValueError("Original EB files changed during this pilot")
    state.update(eb_references=references,eb_origin=origin)
    if args.evaluate_saved_target:
        group = state["groups"][args.evaluate_saved_target]
        if "saved_endpoint_evaluation" not in group:
            evaluate_saved_target(properties,group,references[str(group["case"]["beta_deg"])])
            save(state)
        print(args.evaluate_saved_target,group["status"],"no root search",flush=True)
        return
    state.setdefault("exact_limit",{})
    for beta in (0,30):
        if str(beta) not in state["exact_limit"]:
            state["exact_limit"][str(beta)] = exact_limit(properties,beta,references[str(beta)])
            save(state)
    budget = [state["local_recovery_attempts"]]
    for case in cases():
        if args.case and case["case_id"] != args.case:
            continue
        if case["case_id"] in state["groups"]:
            print(case["case_id"],"preserved (including qualifications)",flush=True)
            continue
        fixed = [r["Omega"] for r in references[str(case["beta_deg"])]["rows"]]
        previous = [g for g in state["groups"].values() if g["case"]["beta_deg"]==case["beta_deg"] and g["status"]=="COMPLETED"]
        predictor = [r["Omega"] for r in previous[-1]["rows"]] if previous else fixed
        group = solve_group(properties,case,predictor,budget)
        if group["status"] in ("COMPLETED","TARGET_CONFIRMED_GUARD_QUALIFIED"):
            group["e_j"] = [abs(r["Omega"]-value)/value for r,value in zip(group["rows"][:6],fixed[:6])]
            group["e_max"] = max(group["e_j"])
        state["groups"][case["case_id"]] = group
        state["local_recovery_attempts"] = budget[0]
        save(state)
        print(case["case_id"],group["status"],"builds",group["boundary_builds"],"seconds",round(group["seconds"],3),"e_max",group.get("e_max"),flush=True)


if __name__ == "__main__":
    main()
