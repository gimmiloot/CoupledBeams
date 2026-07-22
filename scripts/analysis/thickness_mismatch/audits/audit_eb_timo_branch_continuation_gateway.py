from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import csv
import json
import math
from pathlib import Path
import sys
import time
from typing import Mapping, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.analysis.thickness_mismatch.audits import run_eb_epsilon_apriori_pilot as pilot  # noqa: E402
from scripts.analysis.thickness_mismatch.postprocess import analyze_eb_epsilon_apriori_pilot as postprocess  # noqa: E402
from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402


DEFAULT_OUTPUT_DIR = Path("results") / "eb_timo_branch_continuation_gateway"
DEFAULT_PILOT_OUTPUT_DIR = Path("results") / "eb_epsilon_apriori_pilot_branch_continuation_v1"
DEFAULT_CACHE_DIR = DEFAULT_OUTPUT_DIR / "cache"
STEP3_MANIFEST_PATH = (
    REPO_ROOT
    / "scripts"
    / "analysis"
    / "thickness_mismatch"
    / "audits"
    / "data"
    / "eb_epsilon_lower_envelope_step3_cases.csv"
)
PREVIOUS_STRICT_MATRIX_EVALUATIONS = 2_290_411
PREVIOUS_STRICT_SVD_CALLS = 2_106_161
R_CASES = {
    "R1": 0.0228017578125,
    "R2": 0.0159306640625,
    "R3": 0.015625,
}
BLOCKER_IDS = ("B07", "G01", "G02", "M02")

OUTPUT_NAMES = (
    "branch_parent_spectrum.csv",
    "branch_continuation_steps.csv",
    "branch_seed_refinement_audit.csv",
    "branch_cluster_continuation_audit.csv",
    "branch_global_guard_audit.csv",
    "branch_strict_fallback_audit.csv",
    "branch_k10_guard_summary.csv",
    "branch_primary_vs_verification.csv",
    "branch_close_pair_regression.csv",
    "branch_pilot_legacy_comparison.csv",
    "branch_rule_metrics_comparison.csv",
    "branch_operation_counts.csv",
    "eb_timo_branch_continuation_gateway_report.md",
)


@dataclass(frozen=True)
class Args:
    output_dir: Path
    pilot_output_dir: Path
    cache_dir: Path
    manifest: Path
    baseline_dir: Path
    strict_audit_dir: Path
    legacy_pilot_dir: Path
    auto_complete_pilot_dir: Path
    k_max: int
    n_spectrum_roots: int
    n_candidate_roots: int
    verification_candidate_roots: int
    initial_beta_step: float
    min_beta_step: float
    max_beta_step: float
    seed_half_width: float
    sigma_accept: float
    sigma_ratio_accept: float
    mac_accept: float
    subspace_mac_accept: float
    cluster_gap_absolute: float
    cluster_gap_relative: float
    guard_scan_step: float
    run_global_guard: bool
    allow_strict_fallback: bool
    run_close_pair_regressions: bool
    run_pilot: bool
    reuse_cache: bool
    force_recompute: bool
    force_strict_verification: bool
    pilot_spectrum_method: str
    force: bool
    smoke: bool
    plot_only: bool
    write_step3_manifest_if_ready: bool


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Audit the branch-informed EB/Timoshenko K=10 spectrum gateway.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--pilot-output-dir", type=Path, default=DEFAULT_PILOT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--manifest", type=Path, default=pilot.DEFAULT_MANIFEST)
    parser.add_argument("--baseline-dir", type=Path, default=Path("results") / "eb_epsilon_baseline_thresholds")
    parser.add_argument("--strict-audit-dir", type=Path, default=Path("results") / "eb_timo_general_spectrum_completeness")
    parser.add_argument("--legacy-pilot-dir", type=Path, default=Path("results") / "eb_epsilon_apriori_pilot")
    parser.add_argument("--auto-complete-pilot-dir", type=Path, default=Path("results") / "eb_epsilon_apriori_pilot_complete_spectrum_v1")
    parser.add_argument("--corrected-pilot-dir", type=Path, dest="pilot_output_dir")
    parser.add_argument("--k-max", type=int, default=10)
    parser.add_argument("--n-spectrum-roots", type=int, default=12)
    parser.add_argument("--n-candidate-roots", type=int, default=20)
    parser.add_argument("--verification-candidate-roots", type=int, default=24)
    parser.add_argument("--initial-beta-step", type=float, default=0.5)
    parser.add_argument("--min-beta-step", type=float, default=0.0125)
    parser.add_argument("--max-beta-step", type=float, default=5.0)
    parser.add_argument("--seed-half-width", type=float, default=0.055)
    parser.add_argument("--sigma-accept", type=float, default=5.0e-6)
    parser.add_argument("--sigma-ratio-accept", type=float, default=5.0e-4)
    parser.add_argument("--mac-accept", type=float, default=0.25)
    parser.add_argument("--subspace-mac-accept", type=float, default=0.5)
    parser.add_argument("--cluster-gap-absolute", type=float, default=0.04)
    parser.add_argument("--cluster-gap-relative", type=float, default=3.0e-3)
    parser.add_argument("--guard-scan-step", type=float, default=0.05)
    parser.add_argument("--global-guard", dest="run_global_guard", action="store_true", default=True)
    parser.add_argument("--no-global-guard", dest="run_global_guard", action="store_false")
    parser.add_argument("--strict-fallback", dest="allow_strict_fallback", action="store_true", default=True)
    parser.add_argument("--no-strict-fallback", dest="allow_strict_fallback", action="store_false")
    parser.add_argument("--run-close-pair-regressions", action="store_true", default=True)
    parser.add_argument("--no-close-pair-regressions", dest="run_close_pair_regressions", action="store_false")
    parser.add_argument("--run-pilot", action="store_true", default=True)
    parser.add_argument("--no-run-pilot", dest="run_pilot", action="store_false")
    parser.add_argument(
        "--pilot-spectrum-method",
        choices=(pilot.SPECTRUM_METHOD_LEGACY, pilot.SPECTRUM_METHOD_AUTO, pilot.SPECTRUM_METHOD_BRANCH),
        default=pilot.SPECTRUM_METHOD_BRANCH,
    )
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--force-strict-verification", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--write-step3-manifest-if-ready", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    args = Args(
        output_dir=repo_path(Path(ns.output_dir)),
        pilot_output_dir=repo_path(Path(ns.pilot_output_dir or DEFAULT_PILOT_OUTPUT_DIR)),
        cache_dir=repo_path(Path(ns.cache_dir)),
        manifest=repo_path(Path(ns.manifest)),
        baseline_dir=repo_path(Path(ns.baseline_dir)),
        strict_audit_dir=repo_path(Path(ns.strict_audit_dir)),
        legacy_pilot_dir=repo_path(Path(ns.legacy_pilot_dir)),
        auto_complete_pilot_dir=repo_path(Path(ns.auto_complete_pilot_dir)),
        k_max=int(ns.k_max),
        n_spectrum_roots=int(ns.n_spectrum_roots),
        n_candidate_roots=int(ns.n_candidate_roots),
        verification_candidate_roots=int(ns.verification_candidate_roots),
        initial_beta_step=float(ns.initial_beta_step),
        min_beta_step=float(ns.min_beta_step),
        max_beta_step=float(ns.max_beta_step),
        seed_half_width=float(ns.seed_half_width),
        sigma_accept=float(ns.sigma_accept),
        sigma_ratio_accept=float(ns.sigma_ratio_accept),
        mac_accept=float(ns.mac_accept),
        subspace_mac_accept=float(ns.subspace_mac_accept),
        cluster_gap_absolute=float(ns.cluster_gap_absolute),
        cluster_gap_relative=float(ns.cluster_gap_relative),
        guard_scan_step=float(ns.guard_scan_step),
        run_global_guard=bool(ns.run_global_guard),
        allow_strict_fallback=bool(ns.allow_strict_fallback),
        run_close_pair_regressions=bool(ns.run_close_pair_regressions),
        run_pilot=bool(ns.run_pilot),
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        force_strict_verification=bool(ns.force_strict_verification),
        pilot_spectrum_method=str(ns.pilot_spectrum_method),
        force=bool(ns.force),
        smoke=bool(ns.smoke),
        plot_only=bool(ns.plot_only),
        write_step3_manifest_if_ready=bool(ns.write_step3_manifest_if_ready),
    )
    if args.pilot_spectrum_method != pilot.SPECTRUM_METHOD_BRANCH:
        raise ValueError("the branch gateway requires --pilot-spectrum-method branch_informed_continuation_v1")
    if args.output_dir == args.pilot_output_dir:
        raise ValueError("gateway and pilot output directories must be distinct")
    if args.output_dir in {args.legacy_pilot_dir, args.auto_complete_pilot_dir} or args.pilot_output_dir in {
        args.legacy_pilot_dir,
        args.auto_complete_pilot_dir,
    }:
        raise ValueError("new output directories must not replace legacy/auto-complete pilot directories")
    if args.k_max != 10 or args.n_spectrum_roots < args.k_max + 2:
        raise ValueError("the gateway requires K=10 and n_spectrum_roots >= K+2")
    if args.n_candidate_roots < args.n_spectrum_roots or args.verification_candidate_roots <= args.n_candidate_roots:
        raise ValueError("candidate and verification root margins are invalid")
    if not 0.0 < args.min_beta_step <= args.initial_beta_step <= args.max_beta_step:
        raise ValueError("beta steps must satisfy 0 < min <= initial <= max")
    positive = (
        args.seed_half_width,
        args.sigma_accept,
        args.sigma_ratio_accept,
        args.mac_accept,
        args.subspace_mac_accept,
        args.cluster_gap_absolute,
        args.cluster_gap_relative,
        args.guard_scan_step,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in positive):
        raise ValueError("window, sigma, MAC, cluster, and guard settings must be positive")
    if args.mac_accept > 1.0 or args.subspace_mac_accept > 1.0:
        raise ValueError("MAC thresholds must not exceed one")
    return args


def fmt(value: object) -> object:
    if isinstance(value, float):
        if math.isnan(value):
            return "nan"
        if math.isinf(value):
            return "inf" if value > 0.0 else "-inf"
        return f"{value:.16e}"
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (tuple, list, dict)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    if not fields:
        fields = ["status"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key, "")) for key in fields})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def bool_text(value: object) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def ensure_output_policy(args: Args) -> None:
    existing = [args.output_dir / name for name in OUTPUT_NAMES if (args.output_dir / name).exists()]
    if existing and not args.force:
        raise FileExistsError("gateway outputs already exist; pass --force to replace them")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def case_key(label: str, model: str, geometry: complete.Geometry) -> dict[str, object]:
    return {
        "audit_case": label,
        "model": model,
        "epsilon_0": geometry.epsilon_0,
        "beta_deg": geometry.beta_deg,
        "mu": geometry.mu,
        "eta": geometry.eta,
    }


def collect_result(
    label: str,
    result: branch.BranchContinuationResult,
    tables: dict[str, list[dict[str, object]]],
) -> None:
    key = case_key(label, result.model, result.geometry)
    parents = sorted(result.parent_branches, key=lambda item: (item.Lambda, item.family, item.family_index))
    multiplicity_groups: dict[str, str] = {}
    group_index = 0
    group_value = float("nan")
    for parent_item in parents:
        if not math.isfinite(group_value) or abs(parent_item.Lambda - group_value) > result.settings.root_match_tolerance:
            group_index += 1
            group_value = parent_item.Lambda
        multiplicity_groups[parent_item.branch_id] = f"M{group_index:03d}"
    beta0_geometry = complete.Geometry(
        result.geometry.epsilon_0, 0.0, result.geometry.mu, result.geometry.eta
    )
    parent_provider = complete.model_matrix_provider(result.model, beta0_geometry)
    block_map = branch.BETA0_BLOCK_MAPS[result.model]
    for parent_sorted_index, parent_item in enumerate(parents, start=1):
        rows = block_map.axial_rows if parent_item.family == "axial" else block_map.bending_rows
        columns = block_map.axial_columns if parent_item.family == "axial" else block_map.bending_columns
        matrix = parent_provider(parent_item.Lambda)
        _u, _v, _singular, block_sigma_1, block_sigma_ratio, _nullity = branch._canonical_block_svd(
            matrix[np.ix_(rows, columns)]
        )
        result.operations.parent_matrix_evaluations += 1
        result.operations.parent_block_SVD_calls += 1
        quality = (
            block_sigma_1 <= result.settings.sigma_accept
            and block_sigma_ratio <= result.settings.sigma_ratio_accept
            and parent_item.sigma_1 <= result.settings.sigma_accept
            and parent_item.sigma_ratio <= result.settings.sigma_ratio_accept
        )
        tables["parent"].append(
            {
                **key,
                "epsilon": result.geometry.epsilon_0,
                **asdict(parent_item),
                "parent_family": parent_item.family,
                "parent_family_index": parent_item.family_index,
                "Lambda_beta0": parent_item.Lambda,
                "parent_sorted_index": parent_sorted_index,
                "multiplicity_group": multiplicity_groups[parent_item.branch_id],
                "block_sigma_1": block_sigma_1,
                "block_sigma_ratio": block_sigma_ratio,
                "full_sigma_1": parent_item.sigma_1,
                "full_sigma_ratio": parent_item.sigma_ratio,
                "quality_status": "resolved_parent" if quality else "parent_quality_failed",
            }
        )
    for step_index, step in enumerate(result.steps, start=1):
        step_base = {
            **key,
            "step_index": step_index,
            "beta_from_deg": step.beta_from_deg,
            "beta_to_deg": step.beta_to_deg,
            "beta_previous": step.beta_from_deg,
            "beta_target": step.beta_to_deg,
            "attempted_step_deg": step.attempted_step_deg,
            "accepted": step.accepted,
            "accepted_step": step.accepted,
            "status": step.status,
            "branch_count": len(step.branches),
            "cluster_count": len(step.clusters),
            "step_reduction_reason": "" if step.accepted else step.status,
        }
        tables["steps"].append({**step_base, "record_type": "step_attempt"})
        ordered_branches = sorted(step.branches, key=lambda item: item.Lambda)
        for branch_index, branch_item in enumerate(ordered_branches):
            neighbor_gaps = []
            if branch_index:
                neighbor_gaps.append(branch_item.Lambda - ordered_branches[branch_index - 1].Lambda)
            if branch_index + 1 < len(ordered_branches):
                neighbor_gaps.append(ordered_branches[branch_index + 1].Lambda - branch_item.Lambda)
            tables["steps"].append(
                {
                    **step_base,
                    "record_type": "branch",
                    "branch_or_cluster_id": branch_item.branch_id,
                    "branch_id": branch_item.branch_id,
                    "cluster_id": branch_item.cluster_id,
                    "predicted_Lambda": branch_item.predicted_Lambda,
                    "refined_Lambda": branch_item.Lambda,
                    "MAC": branch_item.self_MAC,
                    "MAC_margin": branch_item.self_MAC - result.settings.mac_accept,
                    "gap": min(neighbor_gaps) if neighbor_gaps else float("nan"),
                    "refinement_status": branch_item.refinement_status,
                }
            )
        for cluster in step.clusters:
            tables["steps"].append(
                {
                    **step_base,
                    "record_type": "cluster",
                    "branch_or_cluster_id": cluster.cluster_id,
                    "cluster_id": cluster.cluster_id,
                    "MAC": cluster.subspace_MAC,
                    "MAC_margin": cluster.subspace_MAC - result.settings.subspace_mac_accept,
                    "refinement_status": cluster.status,
                }
            )
        cluster_by_branch = {
            branch_id: cluster.cluster_id for cluster in step.clusters for branch_id in cluster.branch_ids
        }
        new_refinements = [
            item for item in step.refinements if item.acceptance_status == branch.SEED_REFINED_TO_NEW_ROOT
        ]
        for refinement in step.refinements:
            related = ""
            if refinement.acceptance_status == branch.SEED_REFINED_TO_NEW_ROOT:
                related = refinement.Lambda
            elif refinement.acceptance_status == branch.SEED_REFINED_TO_EXISTING_ROOT and new_refinements:
                related = min(new_refinements, key=lambda item: abs(item.Lambda - refinement.Lambda)).Lambda
            tables["seed"].append(
                {
                    **key,
                    "step_index": step_index,
                    "beta_to_deg": step.beta_to_deg,
                    "seed_source": (
                        "cluster_reduced_window"
                        if refinement.branch_id in cluster_by_branch
                        else "isolated_projected_window"
                    ),
                    **asdict(refinement),
                    "seed_Lambda": refinement.seed,
                    "window": (refinement.window_left, refinement.window_right),
                    "stationary_minimum": refinement.primary_candidate,
                    "refined_Lambda": refinement.Lambda,
                    "sigma_values": (refinement.sigma_1, refinement.sigma_2, refinement.sigma_ratio),
                    "merged_new_rejected_status": refinement.acceptance_status,
                    "related_accepted_root": related,
                }
            )
            tables["verification"].append(
                {
                    **key,
                    "verification_scope": "local_independent_refinement",
                    "step_index": step_index,
                    "beta_deg": step.beta_to_deg,
                    "branch_id": refinement.branch_id,
                    "primary_Lambda": refinement.Lambda,
                    "verification_Lambda": refinement.independent_candidate,
                    "absolute_difference": abs(refinement.Lambda - refinement.independent_candidate),
                    "within_tolerance": refinement.independent_agreement,
                }
            )
        for cluster in step.clusters:
            cluster_refinements = [
                item for item in step.refinements if item.branch_id in cluster.branch_ids
            ]
            full_roots = tuple(
                sorted(
                    item.Lambda
                    for item in cluster_refinements
                    if item.acceptance_status == branch.SEED_REFINED_TO_NEW_ROOT
                )
            )
            tables["cluster"].append(
                {
                    **key,
                    "step_index": step_index,
                    **asdict(cluster),
                    "found_root_count": len(full_roots),
                    "full_roots": full_roots,
                    "cluster_status": cluster.status,
                }
            )
    tables["guard"].append({**key, **asdict(result.guard)})
    for row in result.strict_fallback_audit:
        tables["fallback"].append({**key, **row})
    for row in result.primary_vs_verification:
        tables["verification"].append({**key, "verification_scope": "force_global_strict", **row})
    first_disagreement = next(
        (
            int(row["sorted_index"])
            for row in result.primary_vs_verification
            if not bool(row.get("within_tolerance")) and row.get("sorted_index") is not None
        ),
        "",
    )
    first11_ids = {item.branch_id for item in sorted(result.branches, key=lambda item: item.Lambda)[:11]}
    cluster_ambiguity_count = sum(
        not cluster.resolved
        for step in result.steps
        if step.accepted
        for cluster in step.clusters
        if first11_ids.intersection(cluster.branch_ids)
    )
    seed_only_accepted_count = sum(
        item.refinement_status == "seed_only"
        for item in sorted(result.branches, key=lambda item: item.Lambda)[:11]
    )
    tables["summary"].append(
        {
            **key,
            "algorithm_version": result.algorithm_version,
            "K10_guard_resolved": result.k10_guard_resolved,
            "full12_resolved": result.full12_resolved,
            "spectrum_status": result.spectrum_status,
            "exclusion_reason": result.exclusion_reason,
            "oracle_agreement": result.oracle_agreement,
            "force_verification_agreement": result.force_verification_agreement,
            "root_count": len(result.values),
            "root11": result.values[10] if len(result.values) >= 11 else float("nan"),
            "root12": result.values[11] if len(result.values) >= 12 else float("nan"),
            "roots_1_10_status": "resolved" if result.k10_guard_resolved and len(result.values) >= 10 else "unresolved",
            "root11_status": "resolved_guard" if result.k10_guard_resolved and len(result.values) >= 11 else "unresolved",
            "root12_status": "resolved" if result.full12_resolved else "not_resolved_optional_margin",
            "first_disagreement_index": first_disagreement,
            "unresolved_interval_below_guard": result.guard.unresolved_intervals,
            "seed_only_accepted_count": seed_only_accepted_count,
            "cluster_ambiguity_count": cluster_ambiguity_count,
            "strict_fallback_used": result.operations.strict_fallback_runs > 0,
            "cache_status": result.cache_status,
        }
    )
    tables["operations"].append(
        {
            **key,
            "cost_scope": "branch_primary_plus_triggered_fallbacks",
            **asdict(result.operations),
            "previous_strict_matrix_evaluations_metadata": PREVIOUS_STRICT_MATRIX_EVALUATIONS,
            "previous_strict_full_6x6_SVD_calls_metadata": PREVIOUS_STRICT_SVD_CALLS,
        }
    )


def audit_geometries(cases: Sequence[pilot.PilotCase], *, smoke: bool) -> list[tuple[str, complete.Geometry]]:
    rows: list[tuple[str, complete.Geometry]] = []
    beta_values = (0.0, 0.1) if smoke else (0.0, 0.1, 1.0, 5.0)
    r_items = (("R3", R_CASES["R3"]),) if smoke else tuple(R_CASES.items())
    perturbations = (0.0,) if smoke else (-1.0e-6, 0.0, 1.0e-6)
    for label, epsilon in r_items:
        for perturbation in perturbations:
            suffix = "base" if perturbation == 0.0 else ("minus_1e-6" if perturbation < 0.0 else "plus_1e-6")
            for beta in beta_values:
                rows.append(
                    (
                        f"{label}_{suffix}_beta{beta:g}",
                        complete.Geometry(epsilon + perturbation, beta, 0.0, 0.0),
                    )
                )
    by_id = {case.case_id: case for case in cases}
    blocker_ids = ("B07", "M02") if smoke else BLOCKER_IDS
    for case_id in blocker_ids:
        item = by_id[case_id]
        rows.append((case_id, complete.Geometry(item.epsilon_0, item.beta_deg, item.mu, item.eta)))
    unique: dict[tuple[float, float, float, float, str], tuple[str, complete.Geometry]] = {}
    for label, geometry in rows:
        for model in complete.SUPPORTED_MODELS:
            unique[(geometry.epsilon_0, geometry.beta_deg, geometry.mu, geometry.eta, model)] = (label, geometry)
    return [(label, geometry) for label, geometry in rows]


def close_pair_rows(
    label: str, result: branch.BranchContinuationResult
) -> list[dict[str, object]]:
    roots = result.values[:12]
    rows: list[dict[str, object]] = []
    for index in range(len(roots) - 1):
        gap = roots[index + 1] - roots[index]
        if gap <= max(result.settings.cluster_gap_absolute, result.settings.cluster_gap_relative * roots[index + 1]):
            rows.append(
                {
                    **case_key(label, result.model, result.geometry),
                    "left_sorted_index": index + 1,
                    "right_sorted_index": index + 2,
                    "left_Lambda": roots[index],
                    "right_Lambda": roots[index + 1],
                    "gap": gap,
                    "positive_gap_or_verified_nullity": gap > result.settings.root_match_tolerance,
                    "K10_guard_resolved": result.k10_guard_resolved,
                }
            )
    return rows


def run_pilot(args: Args) -> dict[str, object]:
    pilot_argv = [
        "--manifest", str(args.manifest),
        "--output-dir", str(args.pilot_output_dir),
        "--cache-dir", str(args.cache_dir),
        "--spectrum-method", pilot.SPECTRUM_METHOD_BRANCH,
        "--force",
    ]
    if args.reuse_cache:
        pilot_argv.append("--reuse-cache")
    else:
        pilot_argv.append("--no-reuse-cache")
    if args.force_recompute:
        pilot_argv.append("--force-recompute")
    if args.smoke:
        pilot_argv.append("--smoke")
    generated = pilot.main(pilot_argv)
    analysis_dir = args.pilot_output_dir / "analysis"
    post_argv = [
        "--pilot-dir", str(args.pilot_output_dir),
        "--output-dir", str(analysis_dir),
        "--force",
    ]
    if args.smoke:
        post_argv.extend(("--smoke", "--candidate-grid-sizes", "8"))
    analyzed = postprocess.main(post_argv)
    return {"generated": generated, "analyzed": analyzed}


def pilot_comparison_rows(args: Args) -> list[dict[str, object]]:
    sources = {
        "legacy": args.legacy_pilot_dir,
        "auto_complete_spectrum_v1": args.auto_complete_pilot_dir,
        branch.BRANCH_CONTINUATION_ALGORITHM_VERSION: args.pilot_output_dir,
    }
    geometry = {
        name: {row["case_id"]: row for row in read_csv(path / "epsilon_pilot_geometry_metrics.csv")}
        for name, path in sources.items()
    }
    modes = {
        name: {
            (row["case_id"], int(row["sorted_index"])): row
            for row in read_csv(path / "epsilon_pilot_mode_metrics.csv")
            if row.get("sorted_index", "").isdigit()
        }
        for name, path in sources.items()
    }
    case_ids = sorted(set().union(*(set(rows) for rows in geometry.values())))
    out: list[dict[str, object]] = []
    for case_id in case_ids:
        row: dict[str, object] = {"case_id": case_id}
        for name, records in geometry.items():
            item = records.get(case_id, {})
            row[f"{name}_quality_status"] = item.get("quality_status", "missing")
            row[f"{name}_N_true"] = item.get("N_true", "")
            row[f"{name}_first_true_failed_mode"] = item.get("first_true_failed_mode", "")
            row[f"{name}_exclusion_reason"] = item.get("exclusion_reason", "")
            row[f"{name}_root12_available"] = item.get("root12_available", "")
        row["branch_vs_auto_N_true_equal"] = (
            row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_N_true")
            == row.get("auto_complete_spectrum_v1_N_true")
        )
        branch_modes = {
            index: modes[branch.BRANCH_CONTINUATION_ALGORITHM_VERSION].get((case_id, index), {})
            for index in range(1, 11)
        }
        for reference in ("legacy", "auto_complete_spectrum_v1"):
            differences: list[float] = []
            changed = 0
            for index in range(1, 11):
                reference_mode = modes[reference].get((case_id, index), {})
                branch_mode = branch_modes[index]
                for field in ("Lambda_EB", "Lambda_Timo"):
                    try:
                        difference = abs(float(branch_mode[field]) - float(reference_mode[field]))
                    except (KeyError, TypeError, ValueError):
                        continue
                    differences.append(difference)
                    changed += int(difference > 2.0e-4)
            row[f"branch_vs_{reference}_changed_root_count"] = changed
            row[f"branch_vs_{reference}_max_abs_root_difference"] = (
                max(differences) if differences else float("nan")
            )
        out.append(row)
    return out


def rule_comparison_rows(args: Args) -> list[dict[str, object]]:
    sources = {
        "legacy": args.legacy_pilot_dir / "analysis",
        "auto_complete_spectrum_v1": args.auto_complete_pilot_dir / "analysis",
        branch.BRANCH_CONTINUATION_ALGORITHM_VERSION: args.pilot_output_dir / "analysis",
    }
    out: list[dict[str, object]] = []
    for name, path in sources.items():
        for row in read_csv(path / "epsilon_rule_validation_summary.csv"):
            out.append({"spectrum_source": name, **row})
    return out


def write_step3_manifest(path: Path, baseline_dir: Path | None = None) -> Path:
    threshold_path = (baseline_dir or (REPO_ROOT / "results" / "eb_epsilon_baseline_thresholds")) / (
        "baseline_critical_prefix_thresholds.csv"
    )
    source_rows = [
        row
        for row in read_csv(threshold_path)
        if row.get("threshold_status") == "resolved"
        and row.get("epsilon_near_n")
        and row.get("epsilon_buffer_n")
    ]
    threshold_groups: list[dict[str, object]] = []
    for source in source_rows:
        near = float(source["epsilon_near_n"])
        buffer = float(source["epsilon_buffer_n"])
        matched = next(
            (
                group
                for group in threshold_groups
                if abs(float(group["near"]) - near) <= 1.0e-14
                and abs(float(group["buffer"]) - buffer) <= 1.0e-14
            ),
            None,
        )
        if matched is None:
            threshold_groups.append(
                {"prefixes": [int(source["prefix_n"])], "near": near, "buffer": buffer}
            )
        else:
            prefixes = matched["prefixes"]
            assert isinstance(prefixes, list)
            prefixes.append(int(source["prefix_n"]))
    if not threshold_groups:
        raise ValueError(f"no resolved epsilon_near_n/epsilon_buffer_n rows in {threshold_path}")

    adversarial = (
        (0.1, 0.0, 0.0, "small-beta avoided-crossing probe"),
        (1.0, 0.0, 0.0, "small-beta continuation probe"),
        (45.0, 0.0, 0.0, "mid-angle spectrum probe"),
        (90.0, 0.0, 0.0, "right-angle spectrum probe"),
        (45.0, 0.7, 0.0, "high-mu mid-angle probe"),
        (90.0, 0.7, 0.0, "high-mu right-angle probe"),
        (45.0, 0.5, -0.1, "negative-eta mixed probe"),
        (45.0, 0.5, 0.1, "positive-eta mixed probe"),
        (0.1, 0.7, -0.1, "small-beta high-mu negative-eta probe"),
        (1.0, 0.7, 0.1, "small-beta high-mu positive-eta probe"),
        (45.0, 0.7, -0.1, "high-mu negative-eta mixed probe"),
        (90.0, 0.7, 0.1, "high-mu positive-eta right-angle probe"),
        (5.0, 0.5, -0.1, "low-angle negative-eta mixed probe"),
        (5.0, 0.5, 0.1, "low-angle positive-eta mixed probe"),
    )
    rows: list[dict[str, object]] = []
    source_index = 0
    for group in threshold_groups:
        prefixes = group["prefixes"]
        assert isinstance(prefixes, list)
        prefix_group = "prefixes_" + "_".join(str(value) for value in prefixes)
        for epsilon_source, field in (("epsilon_near_n", "near"), ("epsilon_buffer_n", "buffer")):
            epsilon = float(group[field])
            designs = (
                (0.0, 0.0, 0.0, "baseline control at corrected threshold"),
                adversarial[source_index % len(adversarial)],
            )
            source_index += 1
            for beta, mu, eta, rationale in designs:
                rows.append(
                    {
                        "case_id": f"S3_{len(rows) + 1:02d}",
                        "prefix_group": prefix_group,
                        "epsilon_source": epsilon_source,
                        "epsilon": epsilon,
                        "beta_deg": beta,
                        "mu": mu,
                        "eta": eta,
                        "case_group": "baseline_control" if beta == mu == eta == 0.0 else "adversarial",
                        "adversarial_rationale": rationale,
                        "required_spectrum_method": branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
                        "required_K10_guard": True,
                        "notes": "future step3 proposal only; root search not executed",
                    }
                )
    if not 25 <= len(rows) <= 40:
        raise ValueError(f"step-3 manifest must contain 25..40 compact cases, got {len(rows)}")
    write_csv(path, rows)
    return path


def write_report(
    args: Args,
    tables: Mapping[str, Sequence[Mapping[str, object]]],
    pilot_comparison: Sequence[Mapping[str, object]],
    rule_comparison: Sequence[Mapping[str, object]],
    readiness: Mapping[str, bool],
    *,
    elapsed: float,
    manifest_path: Path | None,
) -> Path:
    def md_value(value: object) -> str:
        if isinstance(value, float):
            return "nan" if not math.isfinite(value) else f"{value:.6g}"
        return str(value).replace("|", "\\|")

    def md_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> list[str]:
        if not rows:
            return ["No rows were produced."]
        return [
            "| " + " | ".join(headers) + " |",
            "| " + " | ".join("---" for _ in headers) + " |",
            *("| " + " | ".join(md_value(value) for value in row) + " |" for row in rows),
        ]

    ready = all(readiness.values())
    summary = tables["summary"]
    k10_count = sum(bool(row.get("K10_guard_resolved")) for row in summary)
    full12_count = sum(bool(row.get("full12_resolved")) for row in summary)
    operation_rows = tables["operations"]
    operation_fields = tuple(branch.BranchOperationCounts.__dataclass_fields__)
    operation_totals = {
        name: sum(int(row.get(name, 0)) for row in operation_rows) for name in operation_fields
    }
    primary_svd_by_result = [
        int(row.get("parent_block_SVD_calls", 0))
        + int(row.get("local_full_6x6_SVD_calls", 0))
        + int(row.get("local_reduced_matrix_SVD_calls", 0))
        + int(row.get("guard_full_6x6_SVD_calls", 0))
        + int(row.get("strict_fallback_full_6x6_SVD_calls", 0))
        for row in operation_rows
    ]
    sorted_primary_svd = sorted(primary_svd_by_result)
    if sorted_primary_svd:
        middle = len(sorted_primary_svd) // 2
        median_primary_svd = (
            float(sorted_primary_svd[middle])
            if len(sorted_primary_svd) % 2
            else 0.5 * (sorted_primary_svd[middle - 1] + sorted_primary_svd[middle])
        )
    else:
        median_primary_svd = float("nan")
    total_primary_svd = sum(primary_svd_by_result)
    mean_primary_svd = total_primary_svd / len(primary_svd_by_result) if primary_svd_by_result else float("nan")
    max_primary_svd = max(primary_svd_by_result, default=0)
    strict_svd_reduction = 1.0 - total_primary_svd / PREVIOUS_STRICT_SVD_CALLS
    k10_per_1000_svd = 1000.0 * k10_count / total_primary_svd if total_primary_svd else float("nan")
    fallback_case_count = sum(int(row.get("strict_fallback_runs", 0)) > 0 for row in operation_rows)
    fallback_fraction = fallback_case_count / len(operation_rows) if operation_rows else 0.0
    seed_counts: dict[str, int] = {}
    for row in tables["seed"]:
        status = str(row.get("acceptance_status", "unknown"))
        seed_counts[status] = seed_counts.get(status, 0) + 1
    rejected_seed_count = sum(
        count for status, count in seed_counts.items() if status.startswith("seed_rejected_")
    )
    maximum_off_block = max(
        (float(row.get("off_block_max_abs", 0.0)) for row in tables["parent"]),
        default=float("nan"),
    )
    oracle_failures = sum(row.get("oracle_agreement") is False for row in summary)
    family_roots: dict[tuple[str, str], list[float]] = {}
    for row in tables["parent"]:
        if row.get("audit_case") != "R3_base_beta0":
            continue
        key = (str(row.get("model")), str(row.get("family")))
        values = family_roots.setdefault(key, [])
        value = float(row.get("Lambda", float("nan")))
        if math.isfinite(value) and all(abs(value - old) > 1.0e-8 for old in values):
            values.append(value)
    parent_rows = []
    for model in complete.SUPPORTED_MODELS:
        block = branch.BETA0_BLOCK_MAPS[model]
        for family in ("bending", "axial"):
            roots = sorted(family_roots.get((model, family), ()))[:6]
            parent_rows.append(
                (
                    model,
                    "R3_base_beta0",
                    family,
                    block.bending_rows if family == "bending" else block.axial_rows,
                    block.bending_columns if family == "bending" else block.axial_columns,
                    ", ".join(f"{value:.7g}" for value in roots),
                )
            )
    close_gap = {
        (str(row.get("audit_case")), str(row.get("model"))): min(
            float(candidate.get("gap", float("inf")))
            for candidate in tables["close"]
            if candidate.get("audit_case") == row.get("audit_case")
            and candidate.get("model") == row.get("model")
        )
        for row in tables["summary"]
        if any(
            candidate.get("audit_case") == row.get("audit_case")
            and candidate.get("model") == row.get("model")
            for candidate in tables["close"]
        )
    }
    fallback_keys = {
        (str(row.get("audit_case")), str(row.get("model"))) for row in tables["fallback"]
    }
    operation_by_key = {
        (str(row.get("audit_case")), str(row.get("model"))): row for row in operation_rows
    }
    close_pair_report_rows = []
    for row in summary:
        label = str(row.get("audit_case"))
        if not (label.startswith(("R1_base_", "R2_base_", "R3_base_"))):
            continue
        key = (label, str(row.get("model")))
        operations = operation_by_key.get(key, {})
        close_pair_report_rows.append(
            (
                label,
                row.get("model"),
                11,
                row.get("root_count"),
                close_gap.get(key, "not detected"),
                row.get("K10_guard_resolved"),
                row.get("full12_resolved"),
                key in fallback_keys,
                operations.get("local_matrix_evaluations", 0),
            )
        )
    blocker_report_rows = [
        (
            row.get("audit_case"),
            row.get("model"),
            row.get("root_count"),
            row.get("root11"),
            row.get("root12"),
            row.get("K10_guard_resolved"),
            row.get("full12_resolved"),
            row.get("exclusion_reason", ""),
        )
        for row in summary
        if row.get("audit_case") in BLOCKER_IDS
    ]
    pilot_summary = [row for row in summary if str(row.get("audit_case", "")).startswith("pilot_")]
    pilot_included = sum(
        row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_quality_status") == "included"
        for row in pilot_comparison
    )
    pilot_excluded = [
        str(row.get("case_id"))
        for row in pilot_comparison
        if row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_quality_status") != "included"
    ]
    changed_vs_legacy = [
        str(row.get("case_id"))
        for row in pilot_comparison
        if int(row.get("branch_vs_legacy_changed_root_count", 0)) > 0
        or row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_N_true")
        != row.get("legacy_N_true")
        or row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_first_true_failed_mode")
        != row.get("legacy_first_true_failed_mode")
    ]
    changed_vs_auto = [
        str(row.get("case_id"))
        for row in pilot_comparison
        if int(row.get("branch_vs_auto_complete_spectrum_v1_changed_root_count", 0)) > 0
        or not bool(row.get("branch_vs_auto_N_true_equal"))
        or row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_first_true_failed_mode")
        != row.get("auto_complete_spectrum_v1_first_true_failed_mode")
    ]
    rule_rows = [
        (
            row.get("spectrum_source"),
            row.get("rule"),
            row.get("false_safe_geometry_count"),
            row.get("false_safe_frequency_count"),
            row.get("mean_N_true"),
            row.get("mean_N_hat"),
            row.get("usable_frequency_retention"),
        )
        for row in rule_comparison
        if row.get("split_kind") == "all_held_out"
    ]
    failed_reasons = [name for name, passed in readiness.items() if not passed]
    report = args.output_dir / "eb_timo_branch_continuation_gateway_report.md"
    lines = [
        "# EB/Timoshenko Branch-Informed Continuation Gateway",
        "",
        "## Scope",
        "",
        "- This is step 2.5b only: `K=10` plus the root-11 guard.",
        "- No step-3 lower-envelope/root search was run.",
        "- No determinant, physical model, coefficient ordering, or FEM model was modified.",
        "",
        "## Executive Verdict",
        "",
        f"- Step-3 readiness: `{'ready_for_targeted_step3' if ready else 'not_ready_for_step3'}`.",
        f"- K10 guard resolved: `{k10_count}/{len(summary)}` audited model/geometries.",
        f"- Full 12-root status: `{full12_count}/{len(summary)}` (reported separately; not substituted for K10 readiness).",
        f"- Triggered strict fallbacks: `{fallback_case_count}/{len(operation_rows)}` result records.",
        f"- Elapsed wall-clock time: `{elapsed:.3f}` seconds (auxiliary only).",
        "",
        "## Relation to Previous Tracking",
        "",
        "The older MAC tracker assumed that a complete candidate-root list already existed. This gateway first establishes the local root records, then assigns branch identity. MAC remains an assignment diagnostic for isolated branches, while subspace MAC is used for close-cluster continuation.",
        "",
        "## Beta=0 Parent Spectra",
        "",
        *md_table(("model", "representative", "family", "row indices", "column indices", "first roots"), parent_rows),
        "",
        f"Maximum measured off-block absolute entry: `{maximum_off_block:.3e}`. Straight-family oracle failures: `{oracle_failures}`. Indices above are zero-based and preserve each production matrix's existing coefficient ordering.",
        "",
        "## Seed Correction",
        "",
        "Seeds only define local windows. Direct terminal acceptance is prohibited; only `seed_refined_to_new_root` creates a new root record.",
        "",
        f"- Seeds evaluated: `{len(tables['seed'])}`.",
        f"- Seeds refined to genuinely new roots: `{seed_counts.get(branch.SEED_REFINED_TO_NEW_ROOT, 0)}`.",
        f"- Seeds merged into an existing stationary root: `{seed_counts.get(branch.SEED_REFINED_TO_EXISTING_ROOT, 0)}`.",
        f"- Rejected seeds: `{rejected_seed_count}`.",
        f"- False duplicate root records removed: `{seed_counts.get(branch.SEED_REFINED_TO_EXISTING_ROOT, 0)}`.",
        "",
        "## Continuation Algorithm",
        "",
        "Isolated branches use a beta predictor only to place a projected local window; every accepted point is independently narrowed and checked by the unchanged full 6x6 matrix. Close branches are continued as a left/right null-subspace cluster, with reduced candidates used only as seeds. Beta steps shrink adaptively on failure. A cheap global guard checks below root 11; strict global search is allowed only after minimum-step failure, and force-global strict verification is kept in a separate cache scope.",
        "",
        f"Accepted beta steps: `{operation_totals['beta_steps_accepted']}/{operation_totals['beta_steps_attempted']}`; step reductions: `{operation_totals['beta_step_reductions']}`; cluster resolutions: `{operation_totals['cluster_steps_resolved']}/{operation_totals['cluster_steps_attempted']}`.",
        "",
        "## Close-Pair Regressions",
        "",
        *md_table(
            ("case", "model", "expected", "found", "minimum gap", "K10", "full12", "fallback", "local matrix evals"),
            close_pair_report_rows,
        ),
        "",
        "## Pilot Blockers",
        "",
        *md_table(("case", "model", "roots", "root11", "root12", "K10", "full12", "reason"), blocker_report_rows),
        "",
        "## 21-Case Pilot",
        "",
        f"- Included pilot geometries: `{pilot_included}/{len(pilot_comparison)}`; excluded: `{', '.join(pilot_excluded) if pilot_excluded else 'none'}`.",
        f"- Model spectra with K10 guard resolved: `{sum(bool(row.get('K10_guard_resolved')) for row in pilot_summary)}/{len(pilot_summary)}`.",
        f"- Model spectra with full12 resolved: `{sum(bool(row.get('full12_resolved')) for row in pilot_summary)}/{len(pilot_summary)}`.",
        f"- Cases with roots/N_true/first-failure changed versus legacy: `{', '.join(changed_vs_legacy) if changed_vs_legacy else 'none'}`.",
        f"- Cases changed versus auto-complete-v1 on the same metrics: `{', '.join(changed_vs_auto) if changed_vs_auto else 'none'}`.",
        "",
        "## K10 vs Full12",
        "",
        "`K10_guard_resolved` certifies roots 1..11 and no unresolved candidate below the root-11 guard. `full12_resolved` additionally requires root 12. A full12 failure is reported but cannot replace or erase a valid K10 certificate.",
        "",
        "## Rule Comparison",
        "",
        *md_table(("source", "rule", "false-safe geom.", "false-safe freq.", "mean N_true", "mean N_hat", "frequency retention"), rule_rows),
        "",
        "The historical pilot still defaults to the legacy spectrum source; this gateway writes its branch-informed pilot to a separate directory.",
        "",
        "## Cost",
        "",
        f"- Local matrix evaluations: `{operation_totals['local_matrix_evaluations']}`; local full 6x6 SVD calls: `{operation_totals['local_full_6x6_SVD_calls']}`.",
        f"- Guard matrix evaluations: `{operation_totals['guard_matrix_evaluations']}`; guard full 6x6 SVD calls: `{operation_totals['guard_full_6x6_SVD_calls']}`; recovered roots: `{operation_totals['guard_recovered_roots']}`.",
        f"- Strict-fallback fraction: `{fallback_fraction:.3%}`; force-verification runs: `{operation_totals['force_verification_runs']}`.",
        f"- Primary SVD primitives per model/geometry (mean/median/maximum): `{mean_primary_svd:.3f}/{median_primary_svd:.3f}/{max_primary_svd}`; force-verification work is excluded from this comparison.",
        f"- Primary SVD-call reduction versus the previous strict-audit metadata: `{strict_svd_reduction:.3%}`; K10-resolved records per 1000 primary SVD primitives: `{k10_per_1000_svd:.6g}`.",
        f"- Cache hits/non-hits in this aggregation: `{sum(row.get('cache_status') == 'hit' for row in summary)}/{sum(row.get('cache_status') != 'hit' for row in summary)}`.",
        f"- Previous full strict audit metadata: `{PREVIOUS_STRICT_MATRIX_EVALUATIONS}` matrix evaluations and `{PREVIOUS_STRICT_SVD_CALLS}` full SVD calls; these are not reclassified as continuation work.",
        "",
        "## Step-3 Readiness",
        "",
        f"Decision: `{'ready_for_targeted_step3' if ready else 'not_ready_for_step3'}`.",
        f"Exact blocking reasons: `{', '.join(failed_reasons) if failed_reasons else 'none'}`.",
        "",
    ]
    for name, passed in readiness.items():
        lines.append(f"- `{name}`: `{'PASS' if passed else 'FAIL'}`.")
    lines.extend(
        [
            "",
            "## Limitations",
            "",
            "- This is numerical continuation, not a mathematical proof of spectral completeness.",
            "- Coverage is limited to the selected pilot and stress geometries.",
            "- The general production solver default remains unchanged.",
            "- No lower-envelope search has been run.",
            "- Timoshenko remains the existing 1D reference model.",
            "- Future counterexamples still require strict-fallback verification.",
            "- No FEM, 3D FEM, Gmsh, or CalculiX task was run.",
        ]
    )
    if manifest_path is not None:
        lines.extend(["", f"A future-only step-3 manifest was written to `{manifest_path.relative_to(REPO_ROOT)}`; no roots were calculated from it."])
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    if args.plot_only:
        analyzed = postprocess.main(
            [
                "--pilot-dir", str(args.pilot_output_dir),
                "--output-dir", str(args.pilot_output_dir / "analysis"),
                "--force",
                *(["--smoke", "--candidate-grid-sizes", "8"] if args.smoke else []),
            ]
        )
        return {"args": args, "plot_only": True, "analyzed": analyzed}
    ensure_output_policy(args)
    cases = pilot.load_manifest(args.manifest)
    pilot.validate_manifest(cases, require_full_manifest=not args.smoke)
    selected_cases = pilot.select_cases(cases, smoke=args.smoke)
    settings = branch.ContinuationSettings(
        requested_roots=args.n_spectrum_roots,
        candidate_roots=args.n_candidate_roots,
        verification_candidate_roots=args.verification_candidate_roots,
        beta_initial_step_deg=args.initial_beta_step,
        beta_min_step_deg=args.min_beta_step,
        beta_max_step_deg=args.max_beta_step,
        seed_half_width=args.seed_half_width,
        sigma_accept=args.sigma_accept,
        sigma_ratio_accept=args.sigma_ratio_accept,
        mac_accept=args.mac_accept,
        subspace_mac_accept=args.subspace_mac_accept,
        cluster_gap_absolute=args.cluster_gap_absolute,
        cluster_gap_relative=args.cluster_gap_relative,
        guard_scan_step=args.guard_scan_step,
        run_global_guard=args.run_global_guard,
        allow_strict_fallback=args.allow_strict_fallback,
    )
    primary_cache = branch.BranchContinuationCache(
        args.cache_dir,
        reuse_cache=args.reuse_cache,
        force_recompute=args.force_recompute,
        verification_scope="primary",
    )
    force_cache = branch.BranchContinuationCache(
        args.cache_dir,
        reuse_cache=args.reuse_cache,
        force_recompute=args.force_recompute,
        verification_scope="force_strict_verification",
    )
    tables: dict[str, list[dict[str, object]]] = {
        name: [] for name in ("parent", "steps", "seed", "cluster", "guard", "fallback", "summary", "verification", "close", "operations")
    }
    audited: dict[tuple[str, str], branch.BranchContinuationResult] = {}
    started = time.perf_counter()
    for label, geometry in audit_geometries(cases, smoke=args.smoke):
        for model in complete.SUPPORTED_MODELS:
            force_sample = args.force_strict_verification and (
                "_base_" in label or label in BLOCKER_IDS
            )
            cache = force_cache if force_sample else primary_cache
            result = cache.resolve(model, geometry, settings)
            audited[(label, model)] = result
            collect_result(label, result, tables)
            if args.run_close_pair_regressions:
                tables["close"].extend(close_pair_rows(label, result))
            print(f"audited {label} {model}: {result.spectrum_status}")

    pilot_products = run_pilot(args) if args.run_pilot else None
    for case in selected_cases if args.run_pilot else ():
        geometry = complete.Geometry(case.epsilon_0, case.beta_deg, case.mu, case.eta)
        for model in complete.SUPPORTED_MODELS:
            label = f"pilot_{case.case_id}"
            result = primary_cache.resolve(model, geometry, settings)
            audited[(label, model)] = result
            collect_result(label, result, tables)
            if args.run_close_pair_regressions:
                tables["close"].extend(close_pair_rows(label, result))

    pilot_compare = pilot_comparison_rows(args) if args.run_pilot else []
    rule_compare = rule_comparison_rows(args) if args.run_pilot else []
    write_csv(args.output_dir / "branch_parent_spectrum.csv", tables["parent"])
    write_csv(args.output_dir / "branch_continuation_steps.csv", tables["steps"])
    write_csv(args.output_dir / "branch_seed_refinement_audit.csv", tables["seed"])
    write_csv(args.output_dir / "branch_cluster_continuation_audit.csv", tables["cluster"])
    write_csv(args.output_dir / "branch_global_guard_audit.csv", tables["guard"])
    write_csv(args.output_dir / "branch_strict_fallback_audit.csv", tables["fallback"])
    write_csv(args.output_dir / "branch_k10_guard_summary.csv", tables["summary"])
    write_csv(args.output_dir / "branch_primary_vs_verification.csv", tables["verification"])
    write_csv(args.output_dir / "branch_close_pair_regression.csv", tables["close"])
    write_csv(args.output_dir / "branch_pilot_legacy_comparison.csv", pilot_compare)
    write_csv(args.output_dir / "branch_rule_metrics_comparison.csv", rule_compare)
    write_csv(args.output_dir / "branch_operation_counts.csv", tables["operations"])

    pilot_geometry = (
        pilot_products["generated"]["geometry_rows"]  # type: ignore[index]
        if pilot_products is not None
        else []
    )
    readiness = {
        "all_audited_K10_guard_resolved": all(item.k10_guard_resolved for item in audited.values()),
        "all_pilot_geometries_included": bool(args.run_pilot) and all(
            row.get("quality_status") == "included" for row in pilot_geometry
        ),
        "no_seed_only_candidates": all(
            row.get("acceptance_status") != "seed_only" for row in tables["seed"]
        ),
        "no_unresolved_interval_below_guard": all(item.guard.passed for item in audited.values()),
        "all_accepted_clusters_resolved": all(
            cluster.resolved
            for item in audited.values()
            for step in item.steps if step.accepted
            for cluster in step.clusters
            if {value.branch_id for value in item.branches[:11]}.intersection(cluster.branch_ids)
        ),
        "straight_oracle_regression": all(item.oracle_agreement is not False for item in audited.values()),
        "local_independent_refinement": all(
            bool(row.get("within_tolerance"))
            for row in tables["verification"]
            if row.get("verification_scope") == "local_independent_refinement"
            and row.get("branch_id") in {
                value.branch_id
                for item in audited.values()
                for value in item.branches[:11]
            }
            and any(
                seed.get("audit_case") == row.get("audit_case")
                and seed.get("model") == row.get("model")
                and seed.get("step_index") == row.get("step_index")
                and seed.get("branch_id") == row.get("branch_id")
                and seed.get("acceptance_status") == branch.SEED_REFINED_TO_NEW_ROOT
                for seed in tables["seed"]
            )
        ),
        "force_global_strict_verification": args.force_strict_verification and all(
            item.force_verification_agreement is True
            for (label, _model), item in audited.items()
            if "_base_" in label or label in BLOCKER_IDS
        ),
        "no_direct_seed_acceptance_in_fallback": all(
            not bool(row.get("direct_seed_acceptance_present")) for row in tables["fallback"]
        ),
        "pilot_branch_vs_auto_N_true_consistent": all(
            bool(row.get("branch_vs_auto_N_true_equal"))
            for row in pilot_compare
            if row.get(f"{branch.BRANCH_CONTINUATION_ALGORITHM_VERSION}_quality_status") != "missing"
            and row.get("auto_complete_spectrum_v1_quality_status") != "missing"
        ),
    }
    manifest_path: Path | None = None
    if all(readiness.values()) and not args.smoke and args.write_step3_manifest_if_ready:
        manifest_path = write_step3_manifest(STEP3_MANIFEST_PATH, args.baseline_dir)
    elapsed = time.perf_counter() - started
    report = write_report(
        args,
        tables,
        pilot_compare,
        rule_compare,
        readiness,
        elapsed=elapsed,
        manifest_path=manifest_path,
    )
    print(f"gateway readiness: {'READY' if all(readiness.values()) else 'NOT READY'}")
    print(f"output: {args.output_dir}")
    return {
        "args": args,
        "tables": tables,
        "pilot": pilot_products,
        "pilot_comparison": pilot_compare,
        "rule_comparison": rule_compare,
        "readiness": readiness,
        "manifest_path": manifest_path,
        "report": report,
        "elapsed_seconds": elapsed,
    }


if __name__ == "__main__":
    main()
