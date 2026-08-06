"""Run the finite article-facing EB/Timoshenko epsilon upper-envelope grid.

The script is an orchestration and persistence layer.  It does not implement
or alter either scientific model.  Sorted spectra come from the existing
general completeness solver; suspicious cases are checked with the existing
branch-informed strict gateway.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import socket
import statistics
import subprocess
import sys
import time
from typing import Mapping, Sequence


for _thread_variable in (
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_thread_variable] = "1"

import numpy as np
import scipy


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.lib import article_epsilon_upper_envelope as workflow  # noqa: E402
from scripts.lib import article_epsilon_prefix_optimization as prefix  # noqa: E402
from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import article_epsilon_family_inventory_integration as family_integration  # noqa: E402
from scripts.lib import article_epsilon_family_reconciliation as family_reconciliation  # noqa: E402
from scripts.lib import article_epsilon_compact_certificates as compact_certificates  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402
from scripts.lib import variable_length_timoshenko as timo  # noqa: E402


COARSE_GRID_RUNNER_VERSION = "article_epsilon_coarse_grid_runner_v4_compact_certificates"
COARSE_GRID_OUTPUT_DIR = (
    Path("results") / "article_epsilon_upper_envelope" / "coarse_grid_v1"
)
CONTROL_SAMPLE_SEED = 20260803
SCIENTIFIC_SOLVER_FILES = (
    "scripts/lib/variable_length_timoshenko.py",
    "scripts/lib/general_spectrum_completeness.py",
    "src/my_project/analytic/formulas_thickness_mismatch.py",
    "scripts/lib/branch_informed_spectrum_continuation.py",
    "scripts/lib/article_epsilon_prefix_optimization.py",
    "scripts/lib/article_epsilon_upper_envelope.py",
)


SPECTRA_FIELDS = (
    "case_id",
    "case_identity",
    "grid_group",
    "regression_label",
    "claim_eligible",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "model",
    "sorted_index",
    "root_role",
    "Lambda",
    "Lambda_squared",
    "solver_status",
    "mode_status",
    "root_quality",
    "guard_status",
    "strict_verification_status",
    "strict_trigger_reasons",
    "strict_policy_effective",
    "cache_provenance",
    "unresolved_reason",
    "cluster_id",
    "cluster_member_index",
    "cluster_size",
    "detected_nullity",
    "track_multiplicity",
    "multiplicity_status",
    "branch_id",
    "parent_family",
    "branch_reordered",
    "detection_sources",
    "sigma_1",
    "sigma_ratio",
    "source_path",
    "cache_source_path",
    "omega",
    "omega_over_cutoff_1",
    "omega_over_cutoff_2",
    "max_cutoff_ratio",
    "timo_basis_regime_arm1",
    "timo_basis_regime_arm2",
    "execution_status",
    "N_true_cached",
    "first_failed_mode",
    "first_failed_delta_f",
    "prefix_guard_status",
    "prefix_guard_resolved_through",
    "full_K10_guard_status",
    "early_stop_used",
    "early_stop_reason",
)

CASE_EXECUTION_FIELDS = (
    "case_id",
    "execution_status",
    "N_true_cached",
    "first_failed_mode",
    "first_failed_delta_f",
    "prefix_guard_status",
    "prefix_guard_resolved_through",
    "full_K10_guard_status",
    "early_stop_used",
    "early_stop_reason",
    "strict_verification_status",
    "strict_trigger_reasons",
    "strict_policy_effective",
    "cache_provenance",
    "unresolved_reason",
)

RUNTIME_FIELDS = (
    "phase",
    "ordinal",
    "phase_total",
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "execution_status",
    "N_true",
    "first_failed_mode",
    "wall_seconds",
    "primary_seconds_cumulative",
    "verification_seconds_cumulative",
    "strict_seconds_cumulative",
    "partial_cache_load_status",
    "cache_hit_current_invocation",
)

SANITY_FIELDS = (
    "sanity_role",
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "execution_status",
    "N_true",
    "delta_f_5",
    "strict_verification_status",
    "sanity_status",
)

CONTROL_MANIFEST_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "N_true_prefix",
    "first_failed_mode_prefix",
    "first_failed_delta_f_prefix",
    "thin_0p1_flag",
    "basis_regimes",
    "first_failure_stratum",
    "selection_reasons",
)

CONTROL_RESULT_FIELDS = (
    *CONTROL_MANIFEST_FIELDS,
    "full_execution_status",
    "full_N_true",
    "full_first_failed_mode",
    "full_strict_verification_status",
    "full_K10_guard_status",
    "wall_seconds",
    "cache_hit_current_invocation",
    *(f"full_EB_root_{index}" for index in range(1, 12)),
    *(f"full_Timo_root_{index}" for index in range(1, 12)),
    *(f"full_delta_f_{index}" for index in range(1, 11)),
)

CONTROL_COMPARISON_FIELDS = (
    "case_id",
    "prefix_execution_status",
    "full_execution_status",
    "prefix_N_true",
    "full_N_true",
    "N_true_agreement",
    "prefix_first_failed_mode",
    "full_first_failed_mode",
    "first_failed_mode_agreement",
    "compared_root_count_EB",
    "compared_root_count_Timo",
    "max_root_absolute_difference",
    "root_agreement",
    "cluster_identity_agreement",
    "prefix_guard_agreement",
    "execution_status_agreement",
    "full_strict_pass",
    "comparison_status",
    "disagreement_reason",
)

DEFERRED_CURRENT_FIELDS = (
    "case_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "execution_status",
    "first_apparent_failed_mode",
    "required_guard",
    "defer_reason",
    "expensive_escalation_kind",
    "runtime_before_defer",
    "available_EB_roots",
    "available_Timoshenko_roots",
    "diagnostic_cache_path",
)

KNOWN_EXPENSIVE_STRICT_TAIL_CASE_IDS = (
    "AUE_d892e8483e7ab865bee4",
    "AUE_25852c53792842ff9227",
    "AUE_a461708b1bae16fe60fc",
    "AUE_114d1dbd85d90ba98c49",
)


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def _sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _git_text(*args: str) -> str:
    completed = subprocess.run(
        ["git", *args],
        cwd=REPO_ROOT,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    return completed.stdout.strip()


def _scientific_solver_hashes() -> dict[str, str]:
    return {name: _sha256_path(REPO_ROOT / name) for name in SCIENTIFIC_SOLVER_FILES}


def _command_line(argv: Sequence[str] | None) -> str:
    values = list(sys.argv[1:] if argv is None else argv)
    return subprocess.list2cmdline([sys.executable, str(SCRIPT_PATH), *values])


def _pre_run_metadata(
    *,
    args: argparse.Namespace,
    argv: Sequence[str] | None,
    manifest_path: Path,
    manifest: Sequence[workflow.GridPoint],
    started_utc: str,
) -> dict[str, object]:
    configuration = workflow.solver_configuration()
    dirty = _git_text("status", "--short")
    return {
        "runner_version": COARSE_GRID_RUNNER_VERSION,
        "workflow_version": workflow.WORKFLOW_VERSION,
        "run_started_utc": started_utc,
        "repository": {
            "cwd": str(Path.cwd()),
            "git_root": _git_text("rev-parse", "--show-toplevel"),
            "branch": _git_text("branch", "--show-current"),
            "HEAD": _git_text("rev-parse", "HEAD"),
            "dirty_git_status_short": dirty.splitlines() if dirty else [],
        },
        "scientific_solver_file_sha256": _scientific_solver_hashes(),
        "timoshenko_basis_evaluator_version": timo.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
        "cache_schema": prefix.PREFIX_CACHE_SCHEMA_VERSION,
        "solver_configuration": configuration,
        "solver_configuration_hash": workflow.cache_digest(configuration),
        "exact_command": _command_line(argv),
        "environment": {
            "python": platform.python_version(),
            "python_executable": sys.executable,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "hostname": socket.gethostname(),
            "cpu_count": os.cpu_count(),
            "platform": platform.platform(),
        },
        "manifest_sha256": _sha256_path(manifest_path),
        "manifest_unique_count": len({point.case_id for point in manifest}),
        "manifest_identity_unique_count": len({point.case_identity for point in manifest}),
        "manifest_counts": workflow.group_counts(manifest),
        "manifest_max_absolute_mass_residual": max(abs(point.mass_residual) for point in manifest),
        "run_scope": {
            "smoke": bool(args.smoke),
            "base_only": bool(args.base_only),
            "low_angle_only": bool(args.low_angle_only),
            "regressions_only": bool(args.regressions_only),
            "reuse_cache": bool(args.reuse_cache),
            "force": bool(args.force),
            "solver_mode": "prefix-until-failure" if args.prefix_until_failure else "full-k10",
            "prefix_strategy": args.prefix_strategy,
            "strict_policy": args.strict_policy,
            "family_inventory_policy": args.family_inventory_policy,
            "defer_expensive_strict": bool(args.defer_expensive_strict),
            "model_scientific_scope": family_integration.SCIENTIFIC_SCOPE,
            "workers": int(args.workers),
            "envelope_only": bool(args.envelope_only),
            "readiness_sanity_required": not bool(args.smoke or args.main_pass_only),
            "full_k10_controls_required": bool(
                args.prefix_until_failure and not args.smoke and not args.main_pass_only
            ),
            "full_k10_controls_prepared_only": bool(args.main_pass_only),
            "readiness_sanity_only": bool(args.readiness_sanity_only),
            "main_pass_only": bool(args.main_pass_only),
            "skip_existing_unresolved": bool(args.skip_existing_unresolved),
            "skip_deferred": bool(args.skip_deferred),
            "skip_interrupted": bool(args.skip_interrupted),
            "defer_case_list": str(args.defer_case_list or ""),
        },
        "scientific_scope": (
            "empirical upper envelope on the declared finite grid; no continuous-domain proof; "
            "s_max and Timoshenko cutoff are diagnostics only"
        ),
    }


def _validate_full_manifest(manifest: Sequence[workflow.GridPoint]) -> None:
    expected = {"base": 1400, "low_angle": 144, "s3_14_sweep": 8, "regression": 2}
    if len(manifest) != 1554 or workflow.group_counts(manifest) != expected:
        raise RuntimeError("the approved 1554-point manifest contract is not satisfied")
    if len({point.case_id for point in manifest}) != 1554:
        raise RuntimeError("manifest case IDs are not unique")
    if len({point.case_identity for point in manifest}) != 1554:
        raise RuntimeError("manifest full-precision identities are not unique")
    if max(abs(point.mass_residual) for point in manifest) > 1.0e-12:
        raise RuntimeError("manifest mass preservation check failed")
    if not any(not point.thin_0p1_flag for point in manifest):
        raise RuntimeError("manifest appears to have filtered s_max > 0.1 cases")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute the finite article epsilon upper-envelope grid with K=10 and root-11 guards."
    )
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--force", action="store_true")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--base-only", action="store_true")
    selection.add_argument("--low-angle-only", action="store_true")
    selection.add_argument("--regressions-only", action="store_true")
    parser.add_argument("--postprocess-only", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--workers", type=int, choices=(1, 2, 4), default=1)
    parser.add_argument("--readiness-sanity-only", action="store_true")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--full-k10", action="store_true")
    mode.add_argument("--prefix-until-failure", action="store_true")
    parser.add_argument("--prefix-strategy", choices=prefix.PREFIX_STRATEGIES, default="paired")
    parser.add_argument("--strict-policy", choices=prefix.STRICT_POLICIES, default="auto")
    parser.add_argument(
        "--family-inventory-policy",
        choices=family_integration.FAMILY_POLICIES,
        default="off",
        help="Opt-in parent-process sorted-family post-stage; default 'off' preserves legacy behavior.",
    )
    parser.add_argument(
        "--defer-expensive-strict",
        action="store_true",
        help="With family local-repair, defer an unresolved required prefix instead of executing expensive strict.",
    )
    parser.add_argument("--envelope-only", action="store_true")
    parser.add_argument(
        "--main-pass-only",
        action="store_true",
        help="Run only previously not-attempted regular points and prepare controls without solving them.",
    )
    parser.add_argument("--skip-existing-unresolved", action="store_true")
    parser.add_argument("--skip-deferred", action="store_true")
    parser.add_argument("--skip-interrupted", action="store_true")
    parser.add_argument("--defer-case-list", type=Path)
    parser.add_argument(
        "--reconcile-family-local-repair-shadow",
        action="store_true",
        help="Explicitly promote verified shadow results without running scientific code.",
    )
    parser.add_argument("--shadow-dir", type=Path)
    parser.add_argument(
        "--promotion-policy",
        choices=(family_reconciliation.PROMOTION_POLICY,),
        default=family_reconciliation.PROMOTION_POLICY,
    )
    parser.add_argument("--reconcile-only", action="store_true")
    parser.add_argument("--no-new-point-solves", action="store_true")
    parser.add_argument("--build-compact-point-certificates", action="store_true")
    parser.add_argument("--compact-certificate-dir", type=Path)
    compact_mode = parser.add_mutually_exclusive_group()
    compact_mode.add_argument("--compact-only", action="store_true")
    compact_mode.add_argument("--family-post-stage-only", action="store_true")
    compact_mode.add_argument("--epsilon-005-target-diagnostics-only", action="store_true")
    compact_mode.add_argument("--epsilon-005-targeted-resolution", action="store_true")
    parser.add_argument("--use-compact-point-certificates", action="store_true")
    args = parser.parse_args(argv)
    reconciliation_auxiliary = (
        args.shadow_dir is not None or args.reconcile_only
    )
    if reconciliation_auxiliary and not args.reconcile_family_local_repair_shadow:
        parser.error(
            "--shadow-dir/--reconcile-only/--no-new-point-solves require "
            "--reconcile-family-local-repair-shadow"
        )
    if args.reconcile_family_local_repair_shadow:
        if not args.reconcile_only or not args.no_new_point_solves:
            parser.error(
                "reconciliation requires --reconcile-only --no-new-point-solves"
            )
        incompatible = (
            args.smoke or args.reuse_cache or args.force or args.base_only
            or args.low_angle_only or args.regressions_only or args.postprocess_only
            or args.workers != 1 or args.readiness_sanity_only or args.full_k10
            or args.prefix_until_failure or args.envelope_only or args.main_pass_only
            or args.skip_existing_unresolved or args.skip_deferred
            or args.skip_interrupted or args.defer_case_list is not None
            or args.family_inventory_policy != "off" or args.defer_expensive_strict
        )
        if incompatible:
            parser.error("--reconcile-only cannot be combined with solve or resume options")
        return args
    if args.compact_only:
        if not args.build_compact_point_certificates or not args.no_new_point_solves:
            parser.error("--compact-only requires --build-compact-point-certificates --no-new-point-solves")
        incompatible = (
            args.smoke or args.reuse_cache or args.force or args.base_only
            or args.low_angle_only or args.regressions_only or args.postprocess_only
            or args.workers != 1 or args.readiness_sanity_only or args.full_k10
            or args.prefix_until_failure or args.envelope_only or args.main_pass_only
            or args.skip_existing_unresolved or args.skip_deferred
            or args.skip_interrupted or args.defer_case_list is not None
            or args.family_inventory_policy != "off" or args.defer_expensive_strict
            or args.use_compact_point_certificates
        )
        if incompatible:
            parser.error("--compact-only cannot be combined with solve, resume, or postprocess options")
        return args
    if args.family_post_stage_only:
        if not args.use_compact_point_certificates or not args.no_new_point_solves:
            parser.error("--family-post-stage-only requires --use-compact-point-certificates --no-new-point-solves")
        if not args.defer_expensive_strict:
            parser.error("--family-post-stage-only requires --defer-expensive-strict")
        incompatible = (
            args.smoke or args.reuse_cache or args.force or args.base_only
            or args.low_angle_only or args.regressions_only or args.postprocess_only
            or args.workers != 1 or args.readiness_sanity_only or args.full_k10
            or args.prefix_until_failure or args.envelope_only or args.main_pass_only
            or args.skip_existing_unresolved or args.skip_deferred
            or args.skip_interrupted or args.defer_case_list is not None
            or args.family_inventory_policy != "off" or args.build_compact_point_certificates
        )
        if incompatible:
            parser.error("--family-post-stage-only cannot be combined with primary solve or resume options")
        return args
    if args.epsilon_005_target_diagnostics_only or args.epsilon_005_targeted_resolution:
        incompatible = (
            args.smoke or args.force or args.base_only or args.low_angle_only
            or args.regressions_only or args.postprocess_only or args.workers != 1
            or args.readiness_sanity_only or args.full_k10 or args.prefix_until_failure
            or args.envelope_only or args.main_pass_only
            or args.skip_existing_unresolved or args.skip_deferred
            or args.skip_interrupted or args.defer_case_list is not None
            or args.family_inventory_policy != "off" or args.defer_expensive_strict
            or args.build_compact_point_certificates or args.compact_certificate_dir
            or args.use_compact_point_certificates or args.no_new_point_solves
            or args.reconcile_family_local_repair_shadow
        )
        if incompatible:
            parser.error("epsilon_0=0.050 target modes are isolated workers=1 workflows")
        return args
    if args.no_new_point_solves:
        parser.error("--no-new-point-solves requires reconciliation, --compact-only, or --family-post-stage-only")
    if args.build_compact_point_certificates or args.compact_certificate_dir or args.use_compact_point_certificates:
        parser.error("compact-certificate options require --compact-only or --family-post-stage-only")
    if args.force and args.postprocess_only:
        parser.error("--force cannot be combined with --postprocess-only")
    if args.smoke and (args.base_only or args.low_angle_only or args.regressions_only):
        parser.error("--smoke uses its own fixed 16-point manifest")
    if args.envelope_only and not args.prefix_until_failure:
        parser.error("--envelope-only requires --prefix-until-failure")
    if args.workers > 1 and not (args.prefix_until_failure or args.postprocess_only):
        parser.error("workers>1 is supported only for prefix-until-failure or postprocess-only")
    if args.readiness_sanity_only and args.postprocess_only:
        parser.error("--readiness-sanity-only cannot be combined with --postprocess-only")
    if args.family_inventory_policy == "shadow" and not args.postprocess_only:
        parser.error("--family-inventory-policy shadow requires --postprocess-only")
    if args.family_inventory_policy == "local-repair" and not args.main_pass_only:
        parser.error("--family-inventory-policy local-repair requires --main-pass-only")
    if args.family_inventory_policy == "local-repair" and not args.defer_expensive_strict:
        parser.error("--family-inventory-policy local-repair requires --defer-expensive-strict")
    if args.defer_expensive_strict and args.family_inventory_policy != "local-repair":
        parser.error("--defer-expensive-strict requires --family-inventory-policy local-repair")
    if args.main_pass_only and not args.prefix_until_failure and not args.postprocess_only:
        parser.error("--main-pass-only requires --prefix-until-failure")
    if args.main_pass_only and (args.force or args.envelope_only or args.full_k10):
        parser.error("--main-pass-only forbids --force, --envelope-only, and --full-k10")
    if (
        args.skip_existing_unresolved
        or args.skip_deferred
        or args.skip_interrupted
        or args.defer_case_list
    ) and not args.main_pass_only:
        parser.error("main-pass resume filters require --main-pass-only")
    return args


def _deferred_case_rows(
    path: Path | None,
    manifest: Sequence[workflow.GridPoint],
) -> list[dict[str, str]]:
    if path is None:
        return []
    resolved = repo_path(path)
    if not resolved.exists():
        raise FileNotFoundError(f"deferred case list does not exist: {resolved}")
    rows = workflow.read_csv(resolved)
    if not rows or "case_id" not in rows[0]:
        raise RuntimeError("deferred case list must contain a case_id column")
    manifest_ids = {point.case_id for point in manifest}
    identifiers = [str(row.get("case_id", "")) for row in rows]
    if any(not value for value in identifiers):
        raise RuntimeError("deferred case list contains an empty case_id")
    if len(set(identifiers)) != len(identifiers):
        raise RuntimeError("deferred case list contains duplicate case IDs")
    unknown = sorted(set(identifiers) - manifest_ids)
    if unknown:
        raise RuntimeError(f"deferred case IDs are outside the approved manifest: {unknown}")
    return rows


def _prefix_runtime_case_ids(output_dir: Path) -> set[str]:
    path = output_dir / "runtime_by_case.csv"
    if not path.exists():
        return set()
    return {
        str(row.get("case_id", ""))
        for row in workflow.read_csv(path)
        if row.get("phase") == "prefix_sweep" and row.get("case_id")
    }


def _existing_prefix_records(
    manifest: Sequence[workflow.GridPoint],
    output_dir: Path,
    *,
    strategy: str,
    strict_policy: str,
) -> dict[str, str]:
    """Read scalar statuses without retaining decompressed solver payloads."""

    runtime_case_ids = _prefix_runtime_case_ids(output_dir)
    compact_index = output_dir / "compact_point_certificates_v1" / "compact_index.csv"
    compact_statuses: dict[str, str] = {}
    if compact_index.exists():
        compact_statuses = {
            str(row.get("case_id", "")): str(row.get("execution_status", "attempted_unresolved"))
            for row in workflow.read_csv(compact_index)
            if row.get("case_id")
        }
    point_cache = prefix.PartialPointCache(output_dir / "cache" / "prefix", reuse_cache=True, force=False)
    policies = tuple(dict.fromkeys((strict_policy, "auto")))
    statuses: dict[str, str] = {}
    for point in manifest:
        if point.case_id in compact_statuses:
            raw_status = compact_statuses[point.case_id]
            if raw_status == "deferred_expensive_strict":
                statuses[point.case_id] = raw_status
            elif point.case_id not in runtime_case_ids and raw_status not in {
                "resolved_prefix_early_stop", "resolved_full_K10",
            }:
                statuses[point.case_id] = "interrupted_incomplete"
            else:
                statuses[point.case_id] = raw_status
            continue
        raw = None
        for policy in policies:
            raw = point_cache.load(point, strategy=strategy, strict_policy=policy)
            if raw is not None:
                break
        if raw is None:
            statuses[point.case_id] = "not_attempted"
            continue
        raw_status = str(raw.get("execution_status", "attempted_unresolved"))
        if raw_status == "deferred_expensive_strict":
            statuses[point.case_id] = raw_status
        elif point.case_id not in runtime_case_ids and raw_status not in {
            "resolved_prefix_early_stop",
            "resolved_full_K10",
        }:
            statuses[point.case_id] = "interrupted_incomplete"
        else:
            statuses[point.case_id] = raw_status
        del raw
    return statuses


def _resolved_root_values(raw: Mapping[str, object], model: str) -> list[float]:
    explicit = raw.get("roots_already_resolved", {})
    if isinstance(explicit, Mapping) and isinstance(explicit.get(model), Sequence):
        return [float(value) for value in explicit[model]]  # type: ignore[index]
    models = raw.get("models", {})
    state = models.get(model, {}) if isinstance(models, Mapping) else {}
    values = state.get("resolved_roots", ()) if isinstance(state, Mapping) else ()
    return [float(value) for value in values] if isinstance(values, Sequence) else []


def repair_expensive_strict_tail_statuses(
    output_dir: Path,
    case_ids: Sequence[str] = KNOWN_EXPENSIVE_STRICT_TAIL_CASE_IDS,
) -> list[dict[str, object]]:
    """Repair interrupted auto-policy tails using saved payloads only.

    This routine deliberately never constructs an evaluator or a strict
    callback.  It records that the already-observed escalation belongs to the
    later complex-case pass and preserves all model/diagnostic state in place.
    """

    output_dir = repo_path(output_dir)
    manifest = {point.case_id: point for point in workflow.build_manifest()}
    point_cache = prefix.PartialPointCache(
        output_dir / "cache" / "prefix", reuse_cache=True, force=False
    )
    repaired: list[dict[str, object]] = []
    for case_id in case_ids:
        point = manifest.get(str(case_id))
        if point is None:
            raise RuntimeError(f"tail repair case is outside the manifest: {case_id}")
        raw = point_cache.load(point, strategy="paired", strict_policy="auto")
        if raw is None:
            raise RuntimeError(f"tail repair requires the saved auto payload: {case_id}")
        deltas = {
            int(key): float(value)
            for key, value in dict(raw.get("deltas", {})).items()
        }
        first_failed = next(
            (
                mode
                for mode in range(1, workflow.K_MAX + 1)
                if mode in deltas and deltas[mode] > workflow.DELTA_TOL
            ),
            raw.get("first_failed_mode"),
        )
        models = raw.get("models", {})
        guard = max(
            (
                int(state.get("highest_guard_mode", 0))
                for state in models.values()
                if isinstance(state, Mapping)
            ),
            default=0,
        )
        stage_rows = raw.get("stage_timings", ())
        elapsed_before = sum(
            float(row.get("primary_seconds", 0.0))
            + float(row.get("verification_seconds", 0.0))
            + float(row.get("preparation_seconds", 0.0))
            for row in stage_rows
            if isinstance(row, Mapping)
        ) if isinstance(stage_rows, Sequence) else 0.0
        raw.update(
            {
                "execution_status": "deferred_expensive_strict",
                "n_true_status": "unresolved_pending_complex_pass",
                "N_true": None,
                "expensive_escalation_requested": True,
                "expensive_escalation_kind": "force_strict_verification",
                "trigger_reason": raw.get(
                    "strict_trigger_reasons", ["saved_force_strict_entry"]
                ),
                "first_apparent_failed_mode": first_failed,
                "first_failed_mode": first_failed,
                "required_guard_mode": guard,
                "roots_already_resolved": {
                    model: _resolved_root_values(raw, model)
                    for model in complete.SUPPORTED_MODELS
                },
                "primary_status": "saved_primary_and_local_state_available",
                "local_status": "local_prefix_count_and_guard_pass",
                "elapsed_time_before_escalation": elapsed_before,
                "runtime_before_defer": elapsed_before,
                "strict_policy_effective": "deferred_before_full",
                "strict_verification_status": "not_executed_deferred",
                "required_prefix_strict_status": "unresolved_without_expensive_strict",
                "prefix_guard_status": "unresolved_without_expensive_strict",
                "required_prefix_guard_status": "unresolved_without_expensive_strict",
                "upper_spectrum_audit_status": "not_attempted",
                "optional_upper_spectrum_full_audit_status": "not_attempted",
                "full_K10_control_status": "not_attempted",
                "full_K10_guard_status": "not_attempted",
                "defer_reason": "force_strict_verification_required",
                "unresolved_reason": "force_strict_verification_required",
                "force_strict_requested": max(
                    1, int(raw.get("force_strict_requested", 0))
                ),
                "force_strict_executed": 0,
                "deferred_before_force_strict": max(
                    1, int(raw.get("deferred_before_force_strict", 0))
                ),
                "status_repair_zero_solve": True,
                "status_repair_utc": datetime.now(timezone.utc).isoformat(),
            }
        )
        path = point_cache.save(
            point, raw, strategy="paired", strict_policy="auto"
        )
        repaired.append(
            {
                "case_id": case_id,
                "cache_path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "execution_status": raw["execution_status"],
                "root_calculations": 0,
            }
        )
    return repaired


def build_deferred_complex_cases_current(
    output_dir: Path,
    *,
    seed_path: Path | None = None,
) -> list[dict[str, object]]:
    """Merge all saved unresolved/deferred case IDs without any solves."""

    output_dir = repo_path(output_dir)
    manifest = workflow.build_manifest()
    by_id = {point.case_id: point for point in manifest}
    seed = seed_path or output_dir / "deferred_complex_cases_pre_run.csv"
    seed_rows = workflow.read_csv(seed) if seed.exists() else []
    seed_by_id = {str(row["case_id"]): row for row in seed_rows}
    statuses = _existing_prefix_records(
        manifest, output_dir, strategy="paired", strict_policy="main-pass"
    )
    included = set(seed_by_id)
    included.update(
        case_id
        for case_id, status in statuses.items()
        if status
        in {
            "attempted_unresolved",
            "attempted_error",
            "interrupted_incomplete",
            "deferred_expensive_strict",
            "deferred_complex",
        }
    )
    cache = prefix.PartialPointCache(
        output_dir / "cache" / "prefix", reuse_cache=True, force=False
    )
    compact_index_path = output_dir / "compact_point_certificates_v1" / "compact_index.csv"
    compact_paths = {
        str(row.get("case_id", "")): repo_path(Path(str(row.get("certificate_path", ""))))
        for row in (workflow.read_csv(compact_index_path) if compact_index_path.exists() else [])
        if row.get("case_id") and row.get("certificate_path")
    }
    rows: list[dict[str, object]] = []
    for case_id in sorted(included):
        point = by_id.get(case_id)
        if point is None:
            raise RuntimeError(f"deferred case is outside the manifest: {case_id}")
        raw: Mapping[str, object] = {}
        compact_path = compact_paths.get(case_id)
        if compact_path is not None and compact_path.exists():
            certificate = compact_certificates.load_certificate(compact_path)
            raw = compact_certificates.compact_pseudo_payload(certificate)
        else:
            for policy in ("main-pass", "auto"):
                loaded = cache.load(point, strategy="paired", strict_policy=policy)
                if loaded is not None:
                    raw = loaded
                    break
        seed_row = seed_by_id.get(case_id, {})
        selected_policy = "main-pass"
        cache_path = cache.path(point, strategy="paired", strict_policy=selected_policy)
        if not cache_path.exists():
            selected_policy = "auto"
            cache_path = cache.path(point, strategy="paired", strict_policy=selected_policy)
        status = statuses.get(case_id, "not_attempted")
        if status == "not_attempted":
            status = str(seed_row.get("execution_status", "deferred_complex"))
        first_failed = raw.get(
            "first_apparent_failed_mode", raw.get("first_failed_mode", "")
        )
        guard = raw.get(
            "required_guard_mode", raw.get("prefix_guard_resolved_through", "")
        )
        rows.append(
            {
                "case_id": case_id,
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "execution_status": status,
                "first_apparent_failed_mode": first_failed,
                "required_guard": guard,
                "defer_reason": raw.get(
                    "defer_reason",
                    seed_row.get("exact_reason", seed_row.get("reason", "known_complex_case")),
                ),
                "expensive_escalation_kind": raw.get(
                    "expensive_escalation_kind", ""
                ),
                "runtime_before_defer": raw.get(
                    "runtime_before_defer",
                    raw.get("elapsed_time_before_escalation", ""),
                ),
                "available_EB_roots": json.dumps(
                    _resolved_root_values(raw, complete.MODEL_EB), separators=(",", ":")
                ),
                "available_Timoshenko_roots": json.dumps(
                    _resolved_root_values(raw, complete.MODEL_TIMO), separators=(",", ":")
                ),
                "diagnostic_cache_path": (
                    str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/")
                    if cache_path.exists()
                    else str(seed_row.get("diagnostic_cache_path", ""))
                ),
            }
        )
        del raw
    rows.sort(
        key=lambda row: (
            float(row["epsilon_0"]),
            float(row["beta_deg"]),
            float(row["mu"]),
            float(row["eta"]),
            str(row["case_id"]),
        )
    )
    workflow.write_csv(
        output_dir / "deferred_complex_cases_current.csv",
        rows,
        DEFERRED_CURRENT_FIELDS,
    )
    return rows


def _main_pass_selection(
    selected: Sequence[workflow.GridPoint],
    statuses: Mapping[str, str],
    deferred_case_ids: set[str],
) -> list[workflow.GridPoint]:
    """Select only never-attempted, non-deferred points for the regular pass."""

    return [
        point
        for point in selected
        if statuses.get(point.case_id, "not_attempted") == "not_attempted"
        and point.case_id not in deferred_case_ids
    ]


def _all_available_prefix_payloads(
    manifest: Sequence[workflow.GridPoint],
    output_dir: Path,
    *,
    strategy: str,
    strict_policy: str,
) -> dict[str, dict[str, object]]:
    """Load every compatible prefix cache into export form, with zero solves."""

    point_cache = prefix.PartialPointCache(output_dir / "cache" / "prefix", reuse_cache=True, force=False)
    policies = tuple(dict.fromkeys((strict_policy, "auto")))
    payloads: dict[str, dict[str, object]] = {}
    for point in manifest:
        raw = None
        selected_policy = strict_policy
        for policy in policies:
            raw = point_cache.load(point, strategy=strategy, strict_policy=policy)
            if raw is not None:
                selected_policy = policy
                break
        if raw is None:
            continue
        cache_path = point_cache.path(point, strategy=strategy, strict_policy=selected_policy)
        payloads[point.case_id] = _prefix_case_payload(point, raw, cache_path)
    return payloads


def _cached_s3_regression_gate(output_dir: Path) -> dict[str, object]:
    """Validate the saved S3 evidence without a scientific evaluation."""

    expected = {
        "exact_regression_S3_12": (4, 0.11739469908796035),
        "exact_regression_S3_14": (4, 0.10050934855181458),
    }
    path = output_dir / "readiness_sanity.csv"
    if not path.exists():
        raise FileNotFoundError("main pass requires saved readiness_sanity.csv")
    rows = {str(row.get("sanity_role", "")): row for row in workflow.read_csv(path)}
    observed: dict[str, object] = {}
    for role, (expected_n, expected_delta) in expected.items():
        row = rows.get(role)
        if row is None:
            raise RuntimeError(f"saved readiness sanity is missing {role}")
        n_true = int(row["N_true"])
        delta = float(row["delta_f_5"])
        if (
            row.get("sanity_status") != "PASS"
            or n_true != expected_n
            or not math.isclose(delta, expected_delta, rel_tol=0.0, abs_tol=1.0e-14)
        ):
            raise RuntimeError(f"cached regression gate failed for {role}")
        observed[role] = {"N_true": n_true, "delta_f_5": delta, "status": "PASS"}
    gate_path = output_dir / "resume_readiness_gate.csv"
    gate_rows = workflow.read_csv(gate_path) if gate_path.exists() else []
    if not gate_rows or gate_rows[0].get("resume_readiness_gate") != "PASS":
        raise RuntimeError("main pass requires resume_readiness_gate=PASS")
    return {"status": "PASS", "root_calculations": 0, "cases": observed}


def _model_token(model: str) -> str:
    return "EB" if model == complete.MODEL_EB else "Timo"


def _case_cache_path(output_dir: Path, point: workflow.GridPoint, config_digest: str) -> Path:
    return output_dir / "cache" / "cases" / f"case_{point.case_id}_{config_digest}.json"


def _load_case_cache(
    path: Path,
    point: workflow.GridPoint,
    configuration: Mapping[str, object],
) -> dict[str, object] | None:
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if payload.get("identity") != workflow.case_cache_identity(point, configuration):
        return None
    return payload


def _general_source_path(point: workflow.GridPoint) -> str:
    paths = ["scripts/lib/general_spectrum_completeness.py"]
    if (
        abs(point.beta_deg) <= 1.0e-14
        and abs(point.mu) <= 1.0e-14
        and abs(point.eta) <= 1.0e-14
    ):
        paths.append("scripts/lib/straight_rod_factorized_spectrum.py")
    return ";".join(paths)


def _general_model_payload(
    result: complete.CompleteSpectrumResult,
    *,
    cache_path: Path,
    point: workflow.GridPoint,
) -> dict[str, object]:
    roots: list[dict[str, object]] = []
    for item in result.roots[: workflow.REQUESTED_ROOTS]:
        roots.append(
            {
                "sorted_index": int(item.sorted_index),
                "Lambda": float(item.Lambda),
                "cluster_id": item.root_cluster_id,
                "cluster_member_index": int(item.cluster_member_index),
                "cluster_size": int(item.cluster_size),
                "detected_nullity": int(item.detected_nullity),
                "track_multiplicity": int(item.track_multiplicity),
                "multiplicity_status": item.multiplicity_status,
                "branch_id": "",
                "parent_family": "",
                "branch_reordered": False,
                "detection_sources": ";".join(item.detection_sources),
                "sigma_1": float(item.sigma_1),
                "sigma_ratio": float(item.sigma_ratio),
                "root_quality": (
                    "pass"
                    if item.sigma_1 <= result.settings.sigma_accept
                    and item.sigma_ratio <= result.settings.sigma_ratio_accept
                    else "fail"
                ),
            }
        )
    first11_comparison = tuple(result.primary_vs_verification[: workflow.ROOT_GUARD_INDEX])
    first11_agree = len(first11_comparison) == workflow.ROOT_GUARD_INDEX and all(
        row.get("status") == "pass" for row in first11_comparison
    )
    root11_value = roots[workflow.ROOT_GUARD_INDEX - 1]["Lambda"] if len(roots) >= 11 else float("nan")

    def unresolved_below_guard(entries: Sequence[str]) -> bool:
        if not math.isfinite(float(root11_value)):
            return True
        for entry in entries:
            try:
                left = float(str(entry).split(":", 1)[0])
            except (TypeError, ValueError):
                return True
            if left <= float(root11_value) + max(result.settings.seed_half_width, 2.0 * result.settings.scan_step):
                return True
        return False

    relevant_unresolved = unresolved_below_guard(
        (*result.primary.unresolved_intervals, *result.verification.unresolved_intervals)
    )
    first11_quality = len(roots) >= 11 and all(root["root_quality"] == "pass" for root in roots[:11])
    # K=10 requires root 11 as the right completeness guard.  Failure to fill
    # the more distant 20/24 candidate margins is retained as diagnostics but
    # is not a K10 exclusion when both independent searches agree through 11.
    guard_pass = (
        result.root11_available
        and len(roots) >= workflow.ROOT_GUARD_INDEX
        and first11_quality
        and first11_agree
        and not relevant_unresolved
    )
    return {
        "model": result.model,
        "algorithm_version": result.algorithm_version,
        "spectrum_status": result.spectrum_status,
        "exclusion_reason": result.exclusion_reason,
        "guard_passed": bool(guard_pass),
        "guard_status": "root11_resolved" if guard_pass else "root11_unresolved",
        "independent_agreement": bool(first11_agree),
        "full_requested_agreement": bool(result.independent_agreement),
        "unresolved_interval_below_root11_guard": bool(relevant_unresolved),
        "root12_available": bool(result.root12_available),
        "root12_boundary_warning": bool(result.root12_boundary_warning),
        "cache_status": result.cache_status,
        "cache_source_path": str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "source_path": _general_source_path(point),
        "operations": asdict(result.operations),
        "roots": roots,
    }


def _branch_reordered(result: branch.BranchContinuationResult, item: branch.ContinuedBranch, index: int) -> bool:
    parent_ids = [
        parent.branch_id
        for parent in sorted(
            result.parent_branches,
            key=lambda parent: (parent.Lambda, parent.family, parent.family_index),
        )
    ]
    try:
        return parent_ids.index(item.branch_id) + 1 != int(index)
    except ValueError:
        return True


def _strict_model_payload(
    result: branch.BranchContinuationResult,
    *,
    cache_path: Path,
) -> dict[str, object]:
    ordered = tuple(sorted(result.branches, key=lambda item: item.Lambda))
    roots: list[dict[str, object]] = []
    for index, item in enumerate(ordered[: workflow.REQUESTED_ROOTS], start=1):
        cluster_size = sum(other.cluster_id == item.cluster_id for other in ordered) if item.cluster_id else 1
        roots.append(
            {
                "sorted_index": index,
                "Lambda": float(item.Lambda),
                "cluster_id": item.cluster_id,
                "cluster_member_index": 1,
                "cluster_size": int(cluster_size),
                "detected_nullity": int(item.nullity),
                "track_multiplicity": 1,
                "multiplicity_status": "branch_cluster" if item.cluster_id else "simple_root",
                "branch_id": item.branch_id,
                "parent_family": item.parent_family,
                "branch_reordered": _branch_reordered(result, item, index),
                "detection_sources": item.detection_source,
                "sigma_1": float(item.sigma_1),
                "sigma_ratio": float(item.sigma_ratio),
                "root_quality": (
                    "pass"
                    if item.sigma_1 <= result.settings.sigma_accept
                    and item.sigma_ratio <= result.settings.sigma_ratio_accept
                    else "fail"
                ),
            }
        )
    guard_pass = (
        result.k10_guard_resolved
        and result.guard.passed
        and len(roots) >= workflow.ROOT_GUARD_INDEX
        and result.force_verification_agreement is not False
    )
    return {
        "model": result.model,
        "algorithm_version": result.algorithm_version,
        "spectrum_status": result.spectrum_status,
        "exclusion_reason": result.exclusion_reason,
        "guard_passed": bool(guard_pass),
        "guard_status": result.guard.status,
        "independent_agreement": result.force_verification_agreement is not False,
        "root12_available": len(roots) >= 12,
        "root12_boundary_warning": not result.full12_resolved,
        "cache_status": result.cache_status,
        "cache_source_path": str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "source_path": (
            "scripts/lib/branch_informed_spectrum_continuation.py;"
            "scripts/lib/general_spectrum_completeness.py"
        ),
        "operations": asdict(result.operations),
        "roots": roots,
    }


def _model_values(payload: Mapping[str, object]) -> tuple[float, ...]:
    roots = payload.get("roots", ())
    return tuple(float(item["Lambda"]) for item in roots)  # type: ignore[index]


def _primary_case_status(payload: Mapping[str, object]) -> str:
    models = payload.get("models", {})
    if not isinstance(models, Mapping):
        return "unresolved"
    if all(
        isinstance(models.get(model), Mapping)
        and bool(models[model].get("guard_passed"))  # type: ignore[index,union-attr]
        for model in complete.SUPPORTED_MODELS
    ):
        return "resolved"
    return "unresolved"


def _accepted_model(payload: Mapping[str, object], model: str) -> Mapping[str, object] | None:
    accepted = payload.get("accepted_models", {})
    if isinstance(accepted, Mapping) and isinstance(accepted.get(model), Mapping):
        return accepted[model]  # type: ignore[return-value]
    models = payload.get("models", {})
    if isinstance(models, Mapping) and isinstance(models.get(model), Mapping):
        return models[model]  # type: ignore[return-value]
    return None


def _provisional_n_true(payload: Mapping[str, object]) -> tuple[int | None, tuple[float, ...]]:
    eb = _accepted_model(payload, complete.MODEL_EB)
    tm = _accepted_model(payload, complete.MODEL_TIMO)
    if eb is None or tm is None or not bool(eb.get("guard_passed")) or not bool(tm.get("guard_passed")):
        return None, ()
    eb_values = _model_values(eb)
    tm_values = _model_values(tm)
    if len(eb_values) < 11 or len(tm_values) < 11:
        return None, ()
    deltas = tuple(
        workflow.squared_frequency_delta(eb_values[index], tm_values[index])
        for index in range(workflow.K_MAX)
    )
    return workflow.true_safe_prefix(deltas), deltas


def _case_payload(
    point: workflow.GridPoint,
    models: Mapping[str, Mapping[str, object]],
    configuration: Mapping[str, object],
) -> dict[str, object]:
    payload: dict[str, object] = {
        "identity": workflow.case_cache_identity(point, configuration),
        "point": point.manifest_row(),
        "models": dict(models),
        "accepted_models": dict(models),
        "strict": {
            "status": "not_evaluated",
            "trigger_reasons": [],
            "models": {},
        },
    }
    payload["case_status"] = _primary_case_status(payload)
    return payload


def _root_order_event(previous: Mapping[str, object], current: Mapping[str, object]) -> bool:
    previous_values = _model_values(previous)
    current_values = _model_values(current)
    count = min(len(previous_values), len(current_values), workflow.ROOT_GUARD_INDEX)
    if count < workflow.ROOT_GUARD_INDEX:
        return False
    for index, value in enumerate(previous_values[:count]):
        nearest = min(range(count), key=lambda candidate: abs(current_values[candidate] - value))
        if nearest != index:
            local_gap = abs(current_values[nearest] - current_values[index])
            if local_gap <= max(0.04, 3.0e-3 * max(current_values[nearest], 1.0)):
                return True
    return False


def strict_trigger_map(
    points: Sequence[workflow.GridPoint],
    payloads: Mapping[str, Mapping[str, object]],
) -> dict[str, set[str]]:
    triggers: dict[str, set[str]] = defaultdict(set)
    by_geometry: dict[tuple[float, float, float], list[workflow.GridPoint]] = defaultdict(list)
    provisional: dict[str, tuple[int | None, tuple[float, ...]]] = {}
    for point in points:
        payload = payloads.get(point.case_id)
        if payload is None:
            continue
        provisional[point.case_id] = _provisional_n_true(payload)
        by_geometry[(point.beta_deg, point.mu, point.eta)].append(point)
        if _primary_case_status(payload) != "resolved":
            triggers[point.case_id].add("primary_root11_or_quality_unresolved")
        if point.regression_label:
            triggers[point.case_id].add("exact_regression_point")
        n_true, deltas = provisional[point.case_id]
        if n_true is not None:
            relevant: list[float] = []
            if n_true > 0:
                relevant.append(deltas[n_true - 1])
            if n_true < workflow.K_MAX:
                relevant.append(deltas[n_true])
            if any(0.095 <= value <= 0.105 for value in relevant):
                triggers[point.case_id].add("near_delta_threshold")
        for model in complete.SUPPORTED_MODELS:
            item = _accepted_model(payload, model)
            if item is None:
                continue
            for root in item.get("roots", ())[: workflow.ROOT_GUARD_INDEX]:  # type: ignore[index]
                if (
                    int(root.get("cluster_size", 1)) > 1
                    or int(root.get("track_multiplicity", 1)) > 1
                    or str(root.get("multiplicity_status", "simple_root")) != "simple_root"
                ):
                    triggers[point.case_id].add(f"{_model_token(model)}_cluster_or_multiplicity")

    for geometry_points in by_geometry.values():
        ordered = sorted(geometry_points, key=lambda item: item.epsilon_0)
        for previous_point, current_point in zip(ordered, ordered[1:]):
            previous_payload = payloads[previous_point.case_id]
            current_payload = payloads[current_point.case_id]
            previous_n, _ = provisional.get(previous_point.case_id, (None, ()))
            current_n, _ = provisional.get(current_point.case_id, (None, ()))
            if previous_n is not None and current_n is not None and abs(current_n - previous_n) > 1:
                triggers[previous_point.case_id].add("abrupt_N_true_change")
                triggers[current_point.case_id].add("abrupt_N_true_change")
            for model in complete.SUPPORTED_MODELS:
                previous_model = _accepted_model(previous_payload, model)
                current_model = _accepted_model(current_payload, model)
                if previous_model is None or current_model is None:
                    continue
                if _root_order_event(previous_model, current_model):
                    triggers[previous_point.case_id].add(f"{_model_token(model)}_close_root_order_event")
                    triggers[current_point.case_id].add(f"{_model_token(model)}_close_root_order_event")
            previous_timo = _accepted_model(previous_payload, complete.MODEL_TIMO)
            current_timo = _accepted_model(current_payload, complete.MODEL_TIMO)
            if previous_timo is not None and current_timo is not None:
                left = _model_values(previous_timo)[: workflow.K_MAX]
                right = _model_values(current_timo)[: workflow.K_MAX]
                if len(left) == workflow.K_MAX and len(right) == workflow.K_MAX and any(
                    new > old * (1.0 + 1.0e-4) for old, new in zip(left, right)
                ):
                    triggers[current_point.case_id].add("unusual_Timoshenko_frequency_increase")

    by_epsilon: dict[float, list[workflow.GridPoint]] = defaultdict(list)
    for point in points:
        if point.claim_eligible and point.case_id in provisional:
            by_epsilon[point.epsilon_0].append(point)
    for epsilon_points in by_epsilon.values():
        values = [
            provisional[point.case_id][0]
            for point in epsilon_points
            if provisional[point.case_id][0] is not None
        ]
        if not values:
            continue
        n_up = max(int(value) for value in values if value is not None)
        for point in epsilon_points:
            n_true = provisional[point.case_id][0]
            if n_true == n_up:
                triggers[point.case_id].add("defines_N_up_raw")
            elif n_true == n_up - 1:
                triggers[point.case_id].add("one_mode_below_N_up_raw")
    return triggers


def _compare_first11(primary: Mapping[str, object], strict: Mapping[str, object], tolerance: float = 2.0e-4) -> bool:
    left = _model_values(primary)
    right = _model_values(strict)
    return len(left) >= 11 and len(right) >= 11 and all(
        abs(left[index] - right[index]) <= tolerance for index in range(11)
    )


def _resolve_strict_case(
    point: workflow.GridPoint,
    payload: dict[str, object],
    reasons: Sequence[str],
    cache: branch.BranchContinuationCache,
    settings: branch.ContinuationSettings,
    *,
    force: bool,
    counters: dict[str, int],
) -> None:
    existing = payload.get("strict", {})
    if (
        not force
        and isinstance(existing, Mapping)
        and existing.get("status") in {"pass", "repaired_by_strict", "fail"}
        and isinstance(existing.get("models"), Mapping)
        and len(existing.get("models", {})) == 2
    ):
        counters["strict_case_cache_hits"] += 1
        existing = dict(existing)
        existing["trigger_reasons"] = list(sorted(set(reasons)))
        payload["strict"] = existing
        return

    strict_models: dict[str, object] = {}
    for model in complete.SUPPORTED_MODELS:
        path = cache.path(model, point.geometry, settings)
        try:
            result = cache.resolve(model, point.geometry, settings)
        except Exception as exc:  # preserve the sweep and report the scientific gap
            strict_models[model] = {
                "model": model,
                "algorithm_version": branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
                "spectrum_status": "unresolved_exception",
                "exclusion_reason": f"{type(exc).__name__}: {exc}",
                "guard_passed": False,
                "guard_status": "strict_exception",
                "independent_agreement": False,
                "root12_available": False,
                "root12_boundary_warning": True,
                "cache_status": "exception",
                "cache_source_path": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "source_path": (
                    "scripts/lib/branch_informed_spectrum_continuation.py;"
                    "scripts/lib/general_spectrum_completeness.py"
                ),
                "operations": {},
                "roots": [],
            }
            counters[f"{_model_token(model)}_strict_errors"] += 1
        else:
            strict_models[model] = _strict_model_payload(result, cache_path=path)
            if result.cache_status == "hit":
                counters[f"{_model_token(model)}_strict_cache_hits"] += 1
            else:
                counters[f"{_model_token(model)}_strict_solves"] += 1

    primary_models = payload["models"]
    accepted: dict[str, object] = {}
    all_pass = True
    repaired = False
    comparisons: dict[str, bool] = {}
    for model in complete.SUPPORTED_MODELS:
        primary = primary_models[model]  # type: ignore[index]
        strict_model = strict_models[model]  # type: ignore[index]
        strict_pass = bool(strict_model["guard_passed"])
        primary_values = _model_values(primary)
        comparison_required = len(primary_values) >= workflow.ROOT_GUARD_INDEX
        agrees = _compare_first11(primary, strict_model) if comparison_required else True
        comparisons[model] = agrees
        if not strict_pass or not agrees:
            all_pass = False
        if bool(primary.get("guard_passed")):
            accepted[model] = primary
        elif strict_pass and agrees:
            accepted[model] = strict_model
            repaired = True
        else:
            accepted[model] = primary
    status = "repaired_by_strict" if all_pass and repaired else ("pass" if all_pass else "fail")
    payload["strict"] = {
        "status": status,
        "trigger_reasons": list(sorted(set(reasons))),
        "first11_agreement": comparisons,
        "models": strict_models,
    }
    payload["accepted_models"] = accepted
    payload["case_status"] = "resolved" if all_pass and all(
        bool(accepted[model].get("guard_passed")) for model in complete.SUPPORTED_MODELS  # type: ignore[union-attr]
    ) else "unresolved"
    counters["strict_case_checks"] += 1


def _make_prefix_full_strict_callback(
    strict_cache: branch.BranchContinuationCache,
    settings: branch.ContinuationSettings,
    counters: dict[str, int],
):
    def callback(
        point: workflow.GridPoint,
        sessions: Mapping[str, prefix.IncrementalModelSession],
        required_guard: int,
    ) -> Mapping[str, object]:
        models: dict[str, object] = {}
        all_prefix_pass = True
        all_upper_pass = True
        for model in complete.SUPPORTED_MODELS:
            token = _model_token(model)
            started = time.perf_counter()
            try:
                result = strict_cache.resolve(model, point.geometry, settings)
            except Exception as exc:
                elapsed = time.perf_counter() - started
                models[model] = {
                    "status": "exception",
                    "reason": f"{type(exc).__name__}: {exc}",
                    "seconds": elapsed,
                    "agreement_through_guard": False,
                    "required_prefix_guard_pass": False,
                    "optional_upper_spectrum_full_audit_pass": False,
                }
                counters[f"{token}_strict_errors"] += 1
                all_prefix_pass = False
                all_upper_pass = False
                continue
            elapsed = time.perf_counter() - started
            staged = sessions[model].latest_result
            staged_values = tuple(root.Lambda for root in staged.primary.roots) if staged else ()
            tolerance = staged.settings.root_match_tol if staged else complete.DEFAULT_ROOT_MATCH_TOL
            required_assessment = prefix.assess_branch_strict_through_guard(
                result,
                staged_values,
                required_guard,
                tolerance,
            )
            prefix_pass = bool(required_assessment["agreement"])
            upper_pass = (
                result.k10_guard_resolved
                and result.guard.passed
                and result.force_verification_agreement is not False
            )
            all_prefix_pass = all_prefix_pass and prefix_pass
            all_upper_pass = all_upper_pass and bool(upper_pass)
            models[model] = {
                "status": "pass" if prefix_pass else "fail",
                "seconds": elapsed,
                "cache_status": result.cache_status,
                "agreement_through_guard": prefix_pass,
                "required_prefix_guard_pass": prefix_pass,
                "required_prefix_first_disagreement_root": required_assessment.get(
                    "first_disagreement"
                ),
                "required_prefix_reasons": list(required_assessment.get("reasons", ())),
                "optional_upper_spectrum_full_audit_pass": bool(upper_pass),
                "strict_K10_guard": bool(upper_pass),
                "required_guard": required_guard,
            }
            if result.cache_status == "hit":
                counters[f"{token}_strict_cache_hits"] += 1
            else:
                counters[f"{token}_strict_solves"] += 1
        counters["strict_case_checks"] += 1
        if all_prefix_pass and all_upper_pass:
            status = "full_strict_pass"
        elif all_prefix_pass:
            status = "required_prefix_strict_pass_upper_audit_incomplete"
        else:
            status = "full_strict_failed"
        return {
            "status": status,
            "required_prefix_strict_status": "pass" if all_prefix_pass else "fail",
            "optional_upper_spectrum_full_audit_status": (
                "pass" if all_upper_pass else "incomplete"
            ),
            "required_guard": required_guard,
            "models": models,
        }

    return callback


def _prefix_model_payload(
    point: workflow.GridPoint,
    model: str,
    model_state: Mapping[str, object],
    point_payload: Mapping[str, object],
    cache_path: Path,
) -> dict[str, object]:
    latest = model_state.get("latest_result")
    if not isinstance(latest, Mapping) or not latest:
        return {
            "model": model,
            "guard_passed": False,
            "guard_status": "not_attempted",
            "roots": [],
            "source_path": _general_source_path(point),
            "cache_source_path": str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        }
    result = complete._result_from_payload(latest, "partial_cache")  # type: ignore[attr-defined]
    roots = []
    for item in result.roots:
        roots.append(
            {
                "sorted_index": item.sorted_index,
                "Lambda": item.Lambda,
                "cluster_id": item.root_cluster_id,
                "cluster_member_index": item.cluster_member_index,
                "cluster_size": item.cluster_size,
                "detected_nullity": item.detected_nullity,
                "track_multiplicity": item.track_multiplicity,
                "multiplicity_status": item.multiplicity_status,
                "branch_id": "",
                "parent_family": "",
                "branch_reordered": False,
                "detection_sources": ";".join(item.detection_sources),
                "sigma_1": item.sigma_1,
                "sigma_ratio": item.sigma_ratio,
                "root_quality": (
                    "pass"
                    if item.acceptance_status == "accepted_full_matrix_svd"
                    and item.sigma_1 <= result.settings.sigma_accept
                    else "fail"
                ),
            }
        )
    execution = str(point_payload.get("execution_status", "attempted_unresolved"))
    resolved = execution in {"resolved_prefix_early_stop", "resolved_full_K10"}
    return {
        "model": model,
        "algorithm_version": result.algorithm_version,
        "spectrum_status": result.spectrum_status,
        "exclusion_reason": result.exclusion_reason,
        "guard_passed": resolved,
        "guard_status": point_payload.get("prefix_guard_status", "prefix_guard_unresolved"),
        "independent_agreement": result.independent_agreement,
        "root12_available": result.root12_available,
        "root12_boundary_warning": result.root12_boundary_warning,
        "cache_status": "partial_cache",
        "cache_source_path": str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "source_path": _general_source_path(point),
        "operations": asdict(result.operations),
        "roots": roots,
    }


def _prefix_case_payload(
    point: workflow.GridPoint,
    raw: Mapping[str, object],
    cache_path: Path,
) -> dict[str, object]:
    states = raw.get("models", {})
    models = {
        model: _prefix_model_payload(point, model, states.get(model, {}), raw, cache_path)  # type: ignore[union-attr,arg-type]
        for model in complete.SUPPORTED_MODELS
    }
    resolved = raw.get("execution_status") in {"resolved_prefix_early_stop", "resolved_full_K10"}
    strict = {
        "status": raw.get("strict_verification_status", "not_attempted"),
        "trigger_reasons": raw.get("strict_trigger_reasons", ()),
        "models": raw.get("strict_details", {}),
    }
    return {
        "point": point.manifest_row(),
        "models": models,
        "accepted_models": models,
        "strict": strict,
        "case_status": "resolved" if resolved else "unresolved",
        "execution_status": raw.get("execution_status", "not_attempted"),
        "N_true": raw.get("N_true"),
        "first_failed_mode": raw.get("first_failed_mode"),
        "first_failed_delta_f": raw.get("first_failed_delta_f"),
        "prefix_guard_status": raw.get("prefix_guard_status", "not_attempted"),
        "prefix_guard_resolved_through": raw.get("prefix_guard_resolved_through", 0),
        "full_K10_guard_status": raw.get("full_K10_guard_status", "not_attempted"),
        "early_stop_used": raw.get("early_stop_used", False),
        "early_stop_reason": raw.get("early_stop_reason", ""),
        "deltas": raw.get("deltas", {}),
        "unresolved_reason": raw.get("unresolved_reason", ""),
        "strict_policy_effective": raw.get("strict_policy_effective", "not_attempted"),
        "n_true_status": raw.get("n_true_status", ""),
        "expensive_escalation_requested": raw.get("expensive_escalation_requested", False),
        "expensive_escalation_kind": raw.get("expensive_escalation_kind", ""),
        "trigger_reason": raw.get("trigger_reason", ()),
        "required_guard_mode": raw.get("required_guard_mode", ""),
        "roots_already_resolved": raw.get("roots_already_resolved", {}),
        "clusters_already_resolved": raw.get("clusters_already_resolved", {}),
        "primary_status": raw.get("primary_status", {}),
        "local_status": raw.get("local_status", ""),
        "elapsed_time_before_escalation": raw.get("elapsed_time_before_escalation", ""),
        "defer_reason": raw.get("defer_reason", ""),
        "force_strict_requested": raw.get("force_strict_requested", 0),
        "force_strict_executed": raw.get("force_strict_executed", 0),
        "deferred_before_force_strict": raw.get("deferred_before_force_strict", 0),
        "cache_provenance": str(cache_path.relative_to(REPO_ROOT)).replace("\\", "/"),
    }


def _cumulative_stage_seconds(raw: Mapping[str, object], key: str) -> float:
    rows = raw.get("stage_timings", ())
    if not isinstance(rows, Sequence):
        return 0.0
    return sum(
        float(row.get(key, 0.0))
        for row in rows
        if isinstance(row, Mapping)
    )


def _cumulative_strict_seconds(raw: Mapping[str, object]) -> float:
    details = raw.get("strict_details", {})
    if not isinstance(details, Mapping):
        return 0.0
    models = details.get("models", {})
    if not isinstance(models, Mapping):
        return 0.0
    return sum(
        float(item.get("seconds", 0.0))
        for item in models.values()
        if isinstance(item, Mapping)
    )


def _runtime_record(
    *,
    phase: str,
    ordinal: int,
    phase_total: int,
    point: workflow.GridPoint,
    raw: Mapping[str, object],
    wall_seconds: float,
    cache_status: str,
) -> dict[str, object]:
    return {
        "phase": phase,
        "ordinal": ordinal,
        "phase_total": phase_total,
        "case_id": point.case_id,
        "epsilon_0": point.epsilon_0,
        "beta_deg": point.beta_deg,
        "mu": point.mu,
        "eta": point.eta,
        "execution_status": raw.get("execution_status", "attempted_unresolved"),
        "N_true": raw.get("N_true", ""),
        "first_failed_mode": raw.get("first_failed_mode", ""),
        "wall_seconds": wall_seconds,
        "primary_seconds_cumulative": _cumulative_stage_seconds(raw, "primary_seconds"),
        "verification_seconds_cumulative": _cumulative_stage_seconds(raw, "verification_seconds"),
        "strict_seconds_cumulative": _cumulative_strict_seconds(raw),
        "partial_cache_load_status": cache_status,
        "cache_hit_current_invocation": cache_status == "hit",
    }


def _persist_progress(
    output_dir: Path,
    runtime_records: Mapping[tuple[str, str], Mapping[str, object]],
    *,
    phase: str,
    completed: int,
    total: int,
    payloads: Mapping[str, Mapping[str, object]],
    current: workflow.GridPoint,
    started: float,
    recent_seconds: Sequence[float],
) -> None:
    ordered_runtime = sorted(
        runtime_records.values(),
        key=lambda row: (str(row.get("phase", "")), int(row.get("ordinal", 0)), str(row.get("case_id", ""))),
    )
    workflow.write_csv(output_dir / "runtime_by_case.csv", ordered_runtime, RUNTIME_FIELDS)
    resolved = sum(
        row.get("case_status") == "resolved" for row in payloads.values()
    )
    attempted_unresolved = sum(
        row.get("execution_status") in {"attempted_unresolved", "attempted_error"}
        for row in payloads.values()
    )
    cache_hits = sum(
        bool(row.get("cache_hit_current_invocation"))
        for (row_phase, _case_id), row in runtime_records.items()
        if row_phase == phase
    )
    rolling = statistics.median(recent_seconds[-25:]) if recent_seconds else float("nan")
    eta_seconds = rolling * max(0, total - completed) if math.isfinite(rolling) else float("nan")
    progress = {
        "phase": phase,
        "completed": completed,
        "total": total,
        "resolved": resolved,
        "attempted_unresolved": attempted_unresolved,
        "not_attempted": max(0, total - completed),
        "cache_hits_current_invocation": cache_hits,
        "current_geometry": {
            "beta_deg": current.beta_deg,
            "mu": current.mu,
            "eta": current.eta,
        },
        "current_epsilon_0": current.epsilon_0,
        "elapsed_seconds": time.perf_counter() - started,
        "rolling_median_seconds": rolling,
        "eta_seconds": eta_seconds,
        "updated_utc": datetime.now(timezone.utc).isoformat(),
    }
    workflow.atomic_write_json(output_dir / "logs" / "progress.json", progress)
    print(
        f"[{phase}] {completed}/{total} resolved={resolved} "
        f"attempted_unresolved={attempted_unresolved} not_attempted={max(0, total-completed)} "
        f"cache_hits={cache_hits} geometry=({current.beta_deg:g},{current.mu:g},{current.eta:g}) "
        f"epsilon_0={current.epsilon_0:.17g} elapsed={progress['elapsed_seconds']:.1f}s "
        f"ETA={eta_seconds:.1f}s",
        flush=True,
    )


def _cache_size_bytes(output_dir: Path) -> int:
    return sum(
        path.stat().st_size
        for path in (output_dir / "cache").rglob("*")
        if path.is_file()
    )


def _persist_main_pass_progress(
    output_dir: Path,
    runtime_records: Mapping[tuple[str, str], Mapping[str, object]],
    *,
    payloads: Mapping[str, Mapping[str, object]],
    progress_context: Mapping[str, int],
    new_completed: int,
    selected_total: int,
    cache_hits: int,
    active_geometries: Sequence[tuple[float, float, float]],
    started: float,
    recent_seconds: Sequence[float],
    completed_chains: int,
    failed_chains: int,
    write_checkpoint: bool,
) -> None:
    ordered_runtime = sorted(
        runtime_records.values(),
        key=lambda row: (
            str(row.get("phase", "")),
            int(row.get("ordinal", 0)),
            str(row.get("case_id", "")),
        ),
    )
    workflow.write_csv(output_dir / "runtime_by_case.csv", ordered_runtime, RUNTIME_FIELDS)
    elapsed = time.perf_counter() - started
    new_resolved = sum(row.get("case_status") == "resolved" for row in payloads.values())
    new_unresolved = sum(
        row.get("execution_status") in {"attempted_unresolved", "attempted_error"}
        for row in payloads.values()
    )
    new_deferred = sum(
        row.get("execution_status") == "deferred_expensive_strict"
        for row in payloads.values()
    )
    force_requested = sum(
        int(row.get("force_strict_requested", 0)) for row in payloads.values()
    )
    force_executed = sum(
        int(row.get("force_strict_executed", 0)) for row in payloads.values()
    )
    remaining = max(0, selected_total - new_completed)
    rolling_point = statistics.median(recent_seconds[-50:]) if recent_seconds else float("nan")
    workers = max(1, int(progress_context.get("workers", 1)))
    rolling_eta = (
        rolling_point * remaining / workers if math.isfinite(rolling_point) else float("nan")
    )
    active_workers = len(active_geometries)
    progress = {
        "phase": "main_regular_pass",
        "completed_total": int(progress_context.get("initial_completed", 0)) + new_completed,
        "manifest_total": int(progress_context.get("manifest_total", 1554)),
        "newly_completed_this_run": new_completed,
        "resolved_exact": int(progress_context.get("initial_resolved", 0)) + new_resolved,
        "new_unresolved": new_unresolved,
        "unresolved_total": int(progress_context.get("initial_unresolved", 0)) + new_unresolved,
        "deferred_expensive_strict": int(
            progress_context.get("deferred_expensive_strict", 0)
        ) + new_deferred,
        "deferred": int(progress_context.get("deferred", 0)) + new_deferred,
        "interrupted": int(progress_context.get("interrupted", 0)),
        "remaining_not_attempted": remaining,
        "cache_hits_current_invocation": cache_hits,
        "active_workers": active_workers,
        "configured_workers": workers,
        "current_geometry_chains": [
            {"beta_deg": beta, "mu": mu, "eta": eta}
            for beta, mu, eta in active_geometries
        ],
        "elapsed_seconds": elapsed,
        "rolling_median_point_seconds": rolling_point,
        "rolling_eta_seconds": rolling_eta,
        "force_strict_requested": force_requested,
        "force_strict_executed": force_executed,
        "updated_utc": datetime.now(timezone.utc).isoformat(),
    }
    workflow.atomic_write_json(output_dir / "logs" / "progress.json", progress)
    print(
        "[main_regular_pass] "
        f"{progress['completed_total']}/{progress['manifest_total']} "
        f"new={new_completed} resolved={progress['resolved_exact']} "
        f"new_unresolved={new_unresolved} deferred_strict={progress['deferred_expensive_strict']} "
        f"interrupted={progress['interrupted']} remaining={remaining} "
        f"force_requested={force_requested} force_executed={force_executed} "
        f"cache_hits={cache_hits} active={active_workers} "
        f"chains={active_geometries} elapsed={elapsed:.1f}s ETA={rolling_eta:.1f}s",
        flush=True,
    )
    if write_checkpoint or remaining == 0:
        checkpoint = {
            "timestamp": progress["updated_utc"],
            "completed": progress["completed_total"],
            "resolved": progress["resolved_exact"],
            "unresolved": progress["unresolved_total"],
            "deferred": progress["deferred"],
            "deferred_expensive_strict": progress["deferred_expensive_strict"],
            "not_attempted": remaining,
            "wall_time": elapsed,
            "cache_size": _cache_size_bytes(output_dir),
            "worker_health": {
                "configured": workers,
                "active": active_workers,
                "completed_chains": completed_chains,
                "failed_chains": failed_chains,
                "status": "PASS" if failed_chains == 0 else "DEGRADED",
            },
            "rolling_ETA": rolling_eta,
        }
        workflow.atomic_write_json(output_dir / "main_pass_checkpoint.json", checkpoint)


def _solve_prefix_chain_worker(
    task: tuple[
        int,
        str,
        tuple[workflow.GridPoint, ...],
        bool,
        bool,
        str,
        str,
        bool,
        bool,
        str,
        str,
        str,
    ],
) -> dict[str, object]:
    """Solve one complete geometry chain; return data for parent-only output."""

    (
        chain_ordinal,
        output_text,
        points,
        reuse_cache,
        force,
        strategy,
        strict_policy,
        full_k10,
        force_audit_all,
        phase,
        cache_namespace,
        strict_cache_namespace,
    ) = task
    for variable in (
        "OMP_NUM_THREADS",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        os.environ[variable] = "1"
    output_dir = Path(output_text)
    partial_cache = prefix.PartialPointCache(
        output_dir / "cache" / cache_namespace,
        reuse_cache=reuse_cache,
        force=force,
    )
    strict_cache = branch.BranchContinuationCache(
        output_dir / "cache" / strict_cache_namespace / strategy / strict_policy,
        reuse_cache=reuse_cache,
        force_recompute=force,
        verification_scope="force_strict_verification",
    )
    counters: dict[str, int] = defaultdict(int)
    callback = _make_prefix_full_strict_callback(
        strict_cache,
        workflow.strict_settings(),
        counters,
    )
    continuation: dict[str, tuple[float, ...]] = {}
    payloads: dict[str, dict[str, object]] = {}
    runtime_rows: list[dict[str, object]] = []
    log_rows: list[str] = []
    chain_started = time.perf_counter()
    ordered = tuple(sorted(points, key=lambda point: point.epsilon_0))
    for epsilon_ordinal, point in enumerate(ordered, start=1):
        case_started = time.perf_counter()
        raw = prefix.run_staged_point(
            point,
            cache=partial_cache,
            strategy=strategy,
            strict_policy=strict_policy,
            full_k10=full_k10,
            continuation_roots=continuation,
            strict_callback=callback,
            force_audit=(
                force_audit_all
                or bool(point.regression_label)
                or prefix.is_known_full_path_unresolved_smoke(point)
            ),
        )
        case_seconds = time.perf_counter() - case_started
        cache_path = partial_cache.path(point, strategy=strategy, strict_policy=strict_policy)
        payloads[point.case_id] = _prefix_case_payload(point, raw, cache_path)
        status = str(raw.get("execution_status"))
        counters[status] += 1
        counters["force_strict_requested"] += int(raw.get("force_strict_requested", 0))
        counters["force_strict_executed"] += int(raw.get("force_strict_executed", 0))
        counters["deferred_before_force_strict"] += int(
            raw.get("deferred_before_force_strict", 0)
        )
        current_cache_hit = partial_cache.last_load_status == "hit"
        counters["prefix_case_cache_hits"] += int(current_cache_hit)
        if not current_cache_hit:
            states = raw.get("models", {})
            if isinstance(states, Mapping):
                for model in complete.SUPPORTED_MODELS:
                    if states.get(model):
                        counters[f"{_model_token(model)}_spectrum_solves"] += 1
        if status in {"resolved_prefix_early_stop", "resolved_full_K10"}:
            states = raw.get("models", {})
            for model in complete.SUPPORTED_MODELS:
                state = states.get(model, {}) if isinstance(states, Mapping) else {}
                values = tuple(
                    float(value) for value in state.get("resolved_roots", ())  # type: ignore[union-attr]
                )
                if values:
                    continuation[model] = values
        runtime_rows.append(
            _runtime_record(
                phase=phase,
                ordinal=epsilon_ordinal,
                phase_total=len(ordered),
                point=point,
                raw=raw,
                wall_seconds=case_seconds,
                cache_status=partial_cache.last_load_status,
            )
        )
        log_rows.append(
            f"{phase} chain={chain_ordinal} epsilon={epsilon_ordinal}/{len(ordered)} "
            f"{point.case_id} status={status} N_true={raw.get('N_true')} "
            f"strategy={strategy} strict={strict_policy}"
        )
    return {
        "chain_ordinal": chain_ordinal,
        "payloads": payloads,
        "runtime_rows": runtime_rows,
        "log_rows": log_rows,
        "counters": dict(counters),
        "wall_seconds": time.perf_counter() - chain_started,
        "worker_pid": os.getpid(),
    }


def _solve_selected_points_prefix_parallel(
    points: Sequence[workflow.GridPoint],
    *,
    workers: int,
    output_dir: Path,
    reuse_cache: bool,
    force: bool,
    strategy: str,
    strict_policy: str,
    log_lines: list[str],
    full_k10: bool,
    force_audit_all: bool,
    phase: str,
    runtime_records: dict[tuple[str, str], dict[str, object]] | None,
    cache_namespace: str,
    strict_cache_namespace: str,
    progress_context: Mapping[str, int] | None = None,
) -> tuple[dict[str, dict[str, object]], dict[str, int], dict[str, float]]:
    grouped: dict[tuple[float, float, float], list[workflow.GridPoint]] = defaultdict(list)
    for point in points:
        grouped[(point.beta_deg, point.mu, point.eta)].append(point)
    chains = [
        tuple(sorted(items, key=lambda point: point.epsilon_0))
        for _geometry, items in sorted(grouped.items())
    ]
    tasks = [
        (
            chain_ordinal,
            str(output_dir),
            chain,
            reuse_cache,
            force,
            strategy,
            strict_policy,
            full_k10,
            force_audit_all,
            phase,
            cache_namespace,
            strict_cache_namespace,
        )
        for chain_ordinal, chain in enumerate(chains, start=1)
    ]
    payloads: dict[str, dict[str, object]] = {}
    counters: dict[str, int] = defaultdict(int)
    started = time.perf_counter()
    recent_seconds: list[float] = []
    completed = 0
    completed_chains = 0
    failed_chains = 0
    last_checkpoint_bucket = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(_solve_prefix_chain_worker, task): task for task in tasks}
        for future in as_completed(futures):
            task = futures[future]
            try:
                result = future.result()
            except Exception as exc:  # pragma: no cover - process failure safety net
                failed_chains += 1
                point_cache = prefix.PartialPointCache(
                    output_dir / "cache" / cache_namespace,
                    reuse_cache=True,
                    force=False,
                )
                failed_payloads: dict[str, dict[str, object]] = {}
                failed_rows: list[dict[str, object]] = []
                for epsilon_ordinal, point in enumerate(task[2], start=1):
                    raw = point_cache.load(
                        point, strategy=strategy, strict_policy=strict_policy
                    )
                    if raw is None:
                        raw = prefix.fresh_point_payload(
                            point, strategy=strategy, strict_policy=strict_policy
                        )
                        raw.update(
                            {
                                "execution_status": "attempted_error",
                                "unresolved_reason": f"worker_exception:{type(exc).__name__}:{exc}",
                            }
                        )
                        point_cache.save(
                            point,
                            raw,
                            strategy=strategy,
                            strict_policy=strict_policy,
                        )
                    cache_path = point_cache.path(
                        point, strategy=strategy, strict_policy=strict_policy
                    )
                    failed_payloads[point.case_id] = _prefix_case_payload(point, raw, cache_path)
                    failed_rows.append(
                        _runtime_record(
                            phase=phase,
                            ordinal=epsilon_ordinal,
                            phase_total=len(task[2]),
                            point=point,
                            raw=raw,
                            wall_seconds=0.0,
                            cache_status=point_cache.last_load_status,
                        )
                    )
                result = {
                    "payloads": failed_payloads,
                    "runtime_rows": failed_rows,
                    "log_rows": [
                        f"{phase} chain={task[0]} worker_exception={type(exc).__name__}:{exc}"
                    ],
                    "counters": {"worker_chain_exceptions": 1},
                }
            chain_payloads = result["payloads"]
            if isinstance(chain_payloads, Mapping):
                payloads.update(chain_payloads)  # type: ignore[arg-type]
            for key, value in dict(result["counters"]).items():  # type: ignore[arg-type]
                counters[str(key)] += int(value)
            rows = list(result["runtime_rows"])  # type: ignore[arg-type]
            completed += len(rows)
            completed_chains += 1
            recent_seconds.extend(float(row["wall_seconds"]) for row in rows)
            if runtime_records is not None:
                for row in rows:
                    key = (str(row["phase"]), str(row["case_id"]))
                    if force or key not in runtime_records:
                        runtime_records[key] = dict(row)
                if progress_context is None:
                    current_point = next(
                        point for point in points if point.case_id == str(rows[-1]["case_id"])
                    )
                    _persist_progress(
                        output_dir,
                        runtime_records,
                        phase=phase,
                        completed=completed,
                        total=len(points),
                        payloads=payloads,
                        current=current_point,
                        started=started,
                        recent_seconds=recent_seconds,
                    )
                else:
                    active = []
                    for pending, pending_task in futures.items():
                        if pending is future or pending.done():
                            continue
                        chain = pending_task[2]
                        if chain:
                            active.append(
                                (chain[0].beta_deg, chain[0].mu, chain[0].eta)
                            )
                    bucket = completed // 100
                    _persist_main_pass_progress(
                        output_dir,
                        runtime_records,
                        payloads=payloads,
                        progress_context=progress_context,
                        new_completed=completed,
                        selected_total=len(points),
                        cache_hits=int(counters.get("prefix_case_cache_hits", 0)),
                        active_geometries=active[:workers],
                        started=started,
                        recent_seconds=recent_seconds,
                        completed_chains=completed_chains,
                        failed_chains=failed_chains,
                        write_checkpoint=bucket > last_checkpoint_bucket,
                    )
                    last_checkpoint_bucket = max(last_checkpoint_bucket, bucket)
            log_lines.extend(str(value) for value in result["log_rows"])  # type: ignore[arg-type]
            workflow.atomic_write_text(output_dir / "logs" / "live.log", "\n".join(log_lines) + "\n")
    return payloads, dict(counters), {f"{phase}_seconds": time.perf_counter() - started}


def solve_selected_points_prefix(
    points: Sequence[workflow.GridPoint],
    *,
    output_dir: Path,
    reuse_cache: bool,
    force: bool,
    strategy: str,
    strict_policy: str,
    envelope_only: bool,
    log_lines: list[str],
    full_k10: bool = False,
    force_audit_all: bool = False,
    phase: str = "prefix_sweep",
    runtime_records: dict[tuple[str, str], dict[str, object]] | None = None,
    cache_namespace: str = "prefix",
    strict_cache_namespace: str = "prefix_strict_branch",
    workers: int = 1,
    progress_context: Mapping[str, int] | None = None,
) -> tuple[dict[str, dict[str, object]], dict[str, int], dict[str, float]]:
    if workers not in {1, 2, 4}:
        raise ValueError("workers must be one of 1, 2, or 4")
    if workers > 1:
        if envelope_only:
            raise ValueError("parallel prefix execution does not support envelope-only mode")
        return _solve_selected_points_prefix_parallel(
            points,
            workers=workers,
            output_dir=output_dir,
            reuse_cache=reuse_cache,
            force=force,
            strategy=strategy,
            strict_policy=strict_policy,
            log_lines=log_lines,
            full_k10=full_k10,
            force_audit_all=force_audit_all,
            phase=phase,
            runtime_records=runtime_records,
            cache_namespace=cache_namespace,
            strict_cache_namespace=strict_cache_namespace,
            progress_context=progress_context,
        )
    partial_cache = prefix.PartialPointCache(
        output_dir / "cache" / cache_namespace,
        reuse_cache=reuse_cache,
        force=force,
    )
    strict_cache = branch.BranchContinuationCache(
        output_dir / "cache" / strict_cache_namespace / strategy / strict_policy,
        reuse_cache=reuse_cache,
        force_recompute=force,
        verification_scope="force_strict_verification",
    )
    counters: dict[str, int] = defaultdict(int)
    callback = _make_prefix_full_strict_callback(strict_cache, workflow.strict_settings(), counters)
    payloads: dict[str, dict[str, object]] = {}
    continuation: dict[tuple[float, float, float], dict[str, tuple[float, ...]]] = defaultdict(dict)
    ordered = sorted(points, key=lambda item: (item.beta_deg, item.mu, item.eta, item.epsilon_0))
    saturated_epsilon: set[float] = set()
    started = time.perf_counter()
    recent_seconds: list[float] = []
    for ordinal, point in enumerate(ordered, start=1):
        if envelope_only and point.epsilon_0 in saturated_epsilon:
            continue
        geometry_key = (point.beta_deg, point.mu, point.eta)
        case_started = time.perf_counter()
        raw = prefix.run_staged_point(
            point,
            cache=partial_cache,
            strategy=strategy,
            strict_policy=strict_policy,
            full_k10=full_k10,
            continuation_roots=continuation[geometry_key],
            strict_callback=callback,
            force_audit=(
                force_audit_all
                or bool(point.regression_label)
                or prefix.is_known_full_path_unresolved_smoke(point)
            ),
        )
        case_seconds = time.perf_counter() - case_started
        recent_seconds.append(case_seconds)
        path = partial_cache.path(point, strategy=strategy, strict_policy=strict_policy)
        converted = _prefix_case_payload(point, raw, path)
        payloads[point.case_id] = converted
        status = str(raw.get("execution_status"))
        counters[status] += 1
        current_cache_hit = partial_cache.last_load_status == "hit"
        counters["prefix_case_cache_hits"] += int(current_cache_hit)
        if not current_cache_hit:
            states = raw.get("models", {})
            if isinstance(states, Mapping):
                for model in complete.SUPPORTED_MODELS:
                    if states.get(model):
                        counters[f"{_model_token(model)}_spectrum_solves"] += 1
        if status in {"resolved_prefix_early_stop", "resolved_full_K10"}:
            states = raw.get("models", {})
            for model in complete.SUPPORTED_MODELS:
                state = states.get(model, {}) if isinstance(states, Mapping) else {}
                values = tuple(float(value) for value in state.get("resolved_roots", ()))  # type: ignore[union-attr]
                if values:
                    continuation[geometry_key][model] = values
        if envelope_only and raw.get("N_true") == workflow.K_MAX:
            saturated_epsilon.add(point.epsilon_0)
        log_lines.append(
            f"{phase} {ordinal}/{len(ordered)} {point.case_id} status={status} "
            f"N_true={raw.get('N_true')} strategy={strategy} strict={strict_policy}"
        )
        if runtime_records is not None:
            record = _runtime_record(
                phase=phase,
                ordinal=ordinal,
                phase_total=len(ordered),
                point=point,
                raw=raw,
                wall_seconds=case_seconds,
                cache_status=partial_cache.last_load_status,
            )
            runtime_records[(phase, point.case_id)] = record
            _persist_progress(
                output_dir,
                runtime_records,
                phase=phase,
                completed=ordinal,
                total=len(ordered),
                payloads=payloads,
                current=point,
                started=started,
                recent_seconds=recent_seconds,
            )
        workflow.atomic_write_text(output_dir / "logs" / "live.log", "\n".join(log_lines) + "\n")
    return payloads, dict(counters), {f"{phase}_seconds": time.perf_counter() - started}


def solve_selected_points(
    points: Sequence[workflow.GridPoint],
    *,
    output_dir: Path,
    reuse_cache: bool,
    force: bool,
    configuration: Mapping[str, object],
    log_lines: list[str],
    strict_policy: str = "auto",
) -> tuple[dict[str, dict[str, object]], dict[str, int], dict[str, float]]:
    settings = workflow.primary_settings()
    config_token = workflow.cache_digest(configuration)
    general_cache = complete.GeneralSpectrumCache(
        output_dir / "cache" / "spectra",
        reuse_cache=reuse_cache,
        force_recompute=force,
    )
    counters: dict[str, int] = defaultdict(int)
    timings: dict[str, float] = {}
    payloads: dict[str, dict[str, object]] = {}
    continuation: dict[tuple[float, float, float], dict[str, tuple[float, ...]]] = defaultdict(dict)

    ordered_points = sorted(points, key=lambda item: (item.beta_deg, item.mu, item.eta, item.epsilon_0))
    primary_started = time.perf_counter()
    for ordinal, point in enumerate(ordered_points, start=1):
        case_path = _case_cache_path(output_dir, point, config_token)
        cached = _load_case_cache(case_path, point, configuration) if reuse_cache and not force else None
        geometry_key = (point.beta_deg, point.mu, point.eta)
        if cached is not None:
            payloads[point.case_id] = cached
            counters["case_cache_hits"] += 1
            counters["EB_cache_hits"] += 1
            counters["Timo_cache_hits"] += 1
            for model in complete.SUPPORTED_MODELS:
                item = _accepted_model(cached, model)
                if item is not None and bool(item.get("guard_passed")):
                    continuation[geometry_key][model] = _model_values(item)
            continue

        models: dict[str, Mapping[str, object]] = {}
        eb_seeds = continuation[geometry_key].get(complete.MODEL_EB, ())
        eb_result = general_cache.resolve(
            complete.MODEL_EB,
            point.geometry,
            settings,
            continuation_seeds=eb_seeds,
        )
        eb_path = general_cache.path(complete.MODEL_EB, point.geometry, settings)
        models[complete.MODEL_EB] = _general_model_payload(eb_result, cache_path=eb_path, point=point)
        if eb_result.cache_status == "hit":
            counters["EB_cache_hits"] += 1
        else:
            counters["EB_spectrum_solves"] += 1
        if bool(models[complete.MODEL_EB]["guard_passed"]):
            continuation[geometry_key][complete.MODEL_EB] = _model_values(models[complete.MODEL_EB])

        timo_seeds = continuation[geometry_key].get(complete.MODEL_TIMO, ())
        timo_result = general_cache.resolve(
            complete.MODEL_TIMO,
            point.geometry,
            settings,
            continuation_seeds=timo_seeds,
            eb_seed_roots=_model_values(models[complete.MODEL_EB]),
        )
        timo_path = general_cache.path(complete.MODEL_TIMO, point.geometry, settings)
        models[complete.MODEL_TIMO] = _general_model_payload(timo_result, cache_path=timo_path, point=point)
        if timo_result.cache_status == "hit":
            counters["Timo_cache_hits"] += 1
        else:
            counters["Timo_spectrum_solves"] += 1
        if bool(models[complete.MODEL_TIMO]["guard_passed"]):
            continuation[geometry_key][complete.MODEL_TIMO] = _model_values(models[complete.MODEL_TIMO])

        payload = _case_payload(point, models, configuration)
        payloads[point.case_id] = payload
        workflow.atomic_write_json(case_path, payload)
        log_lines.append(
            f"primary {ordinal}/{len(ordered_points)} {point.case_id} "
            f"eps={point.epsilon_0:.17g} beta={point.beta_deg:g} mu={point.mu:g} eta={point.eta:g} "
            f"status={payload['case_status']}"
        )
    timings["primary_spectrum_seconds"] = time.perf_counter() - primary_started

    trigger_map = strict_trigger_map(points, payloads)
    # Every general-spectrum solve already contains a second, independent
    # half-step verification search with a larger candidate target.  That
    # existing verification satisfies the explicit recheck for ordinary
    # envelope and N_up-1 rows.  The much more expensive branch-informed
    # gateway is reserved for triggers that need branch/cluster diagnostics
    # or repair of a failed root-11 contract.
    primary_independent_reasons = {
        "defines_N_up_raw",
        "one_mode_below_N_up_raw",
        "near_delta_threshold",
        "abrupt_N_true_change",
        "unusual_Timoshenko_frequency_increase",
    }
    branch_trigger_map = {
        case_id: reasons
        for case_id, reasons in trigger_map.items()
        if set(reasons) - primary_independent_reasons
    }
    if strict_policy == "full":
        branch_trigger_map = {
            point.case_id: set(trigger_map.get(point.case_id, ())) | {"explicit_full_strict_policy"}
            for point in points
        }
    strict_cache = branch.BranchContinuationCache(
        output_dir / "cache" / "strict_branch",
        reuse_cache=reuse_cache,
        force_recompute=force,
        verification_scope="force_strict_verification",
    )
    strict = workflow.strict_settings()
    strict_started = time.perf_counter()
    point_by_id = {point.case_id: point for point in points}
    for ordinal, case_id in enumerate(sorted(branch_trigger_map), start=1):
        point = point_by_id[case_id]
        payload = payloads[case_id]
        _resolve_strict_case(
            point,
            payload,
            sorted(branch_trigger_map[case_id]),
            strict_cache,
            strict,
            force=force,
            counters=counters,
        )
        workflow.atomic_write_json(_case_cache_path(output_dir, point, config_token), payload)
        log_lines.append(
            f"strict {ordinal}/{len(branch_trigger_map)} {case_id} "
            f"status={payload.get('strict', {}).get('status')} "  # type: ignore[union-attr]
            f"reasons={';'.join(sorted(branch_trigger_map[case_id]))}"
        )
    for case_id, payload in payloads.items():
        if case_id in branch_trigger_map:
            continue
        reasons = sorted(trigger_map.get(case_id, ()))
        primary_status = _primary_case_status(payload)
        payload["accepted_models"] = payload["models"]
        payload["case_status"] = primary_status
        payload["strict"] = {
            "status": (
                "primary_independent_pass"
                if reasons and primary_status == "resolved"
                else (
                    "not_triggered_primary_independent_pass"
                    if primary_status == "resolved"
                    else "not_triggered_primary_unresolved"
                )
            ),
            "trigger_reasons": reasons,
            "models": {},
        }
        point = point_by_id[case_id]
        workflow.atomic_write_json(_case_cache_path(output_dir, point, config_token), payload)
    timings["strict_verification_seconds"] = time.perf_counter() - strict_started
    counters["strict_triggered_cases"] = len(trigger_map)
    counters["branch_strict_triggered_cases"] = len(branch_trigger_map)
    counters["primary_independent_envelope_rechecks"] = len(trigger_map) - len(branch_trigger_map)
    return payloads, dict(counters), timings


def _all_available_payloads(
    manifest: Sequence[workflow.GridPoint],
    output_dir: Path,
    configuration: Mapping[str, object],
) -> dict[str, dict[str, object]]:
    token = workflow.cache_digest(configuration)
    result: dict[str, dict[str, object]] = {}
    for point in manifest:
        payload = _load_case_cache(_case_cache_path(output_dir, point, token), point, configuration)
        if payload is not None:
            result[point.case_id] = payload
    return result


def export_spectra_long(
    manifest: Sequence[workflow.GridPoint],
    payloads: Mapping[str, Mapping[str, object]],
    output_dir: Path,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for point in manifest:
        payload = payloads.get(point.case_id)
        if payload is None:
            continue
        strict = payload.get("strict", {})
        strict_status = str(strict.get("status", "not_evaluated")) if isinstance(strict, Mapping) else "not_evaluated"
        strict_reasons = (
            ";".join(str(value) for value in strict.get("trigger_reasons", ()))
            if isinstance(strict, Mapping)
            else ""
        )
        execution_status = str(
            payload.get(
                "execution_status",
                "resolved_full_K10" if payload.get("case_status") == "resolved" else "attempted_unresolved",
            )
        )
        for model in complete.SUPPORTED_MODELS:
            model_payload = _accepted_model(payload, model)
            if model_payload is None:
                continue
            model_resolved = bool(model_payload.get("guard_passed"))
            section1 = timo.section_from_epsilon_tau(point.epsilon_0, point.tau1)
            section2 = timo.section_from_epsilon_tau(point.epsilon_0, point.tau2)
            for root in model_payload.get("roots", ()):  # type: ignore[union-attr]
                index = int(root["sorted_index"])
                value = float(root["Lambda"])
                omega = float("nan")
                ratio1 = ratio2 = ratio_max = float("nan")
                regime1 = regime2 = ""
                if model == complete.MODEL_TIMO:
                    omega = timo.project_omega(value, point.epsilon_0)
                    ratio1 = omega / timo.omega_cutoff(section1)
                    ratio2 = omega / timo.omega_cutoff(section2)
                    ratio_max = max(ratio1, ratio2)
                    regime1 = timo.timo_basis(value, point.epsilon_0, section1).regime
                    regime2 = timo.timo_basis(value, point.epsilon_0, section2).regime
                rows.append(
                    {
                        "case_id": point.case_id,
                        "case_identity": point.case_identity,
                        "grid_group": point.grid_group,
                        "regression_label": point.regression_label,
                        "claim_eligible": point.claim_eligible,
                        "epsilon_0": point.epsilon_0,
                        "beta_deg": point.beta_deg,
                        "mu": point.mu,
                        "eta": point.eta,
                        "model": model,
                        "sorted_index": index,
                        "root_role": (
                            "reported_K10"
                            if index <= workflow.K_MAX
                            else ("K10_right_guard" if index == 11 else "optional_right_margin")
                        ),
                        "Lambda": value,
                        "Lambda_squared": value**2,
                        "solver_status": "resolved" if model_resolved else "unresolved",
                        "mode_status": (
                            "not_needed_after_first_failure"
                            if payload.get("first_failed_mode") is not None
                            and index > int(payload["first_failed_mode"]) + 1
                            else "resolved"
                        ),
                        "root_quality": root.get("root_quality", "unknown"),
                        "guard_status": model_payload.get("guard_status", "unresolved"),
                        "strict_verification_status": strict_status,
                        "strict_trigger_reasons": strict_reasons,
                        "cluster_id": root.get("cluster_id", ""),
                        "cluster_member_index": root.get("cluster_member_index", 1),
                        "cluster_size": root.get("cluster_size", 1),
                        "detected_nullity": root.get("detected_nullity", 1),
                        "track_multiplicity": root.get("track_multiplicity", 1),
                        "multiplicity_status": root.get("multiplicity_status", "simple_root"),
                        "branch_id": root.get("branch_id", ""),
                        "parent_family": root.get("parent_family", ""),
                        "branch_reordered": root.get("branch_reordered", False),
                        "detection_sources": root.get("detection_sources", ""),
                        "sigma_1": root.get("sigma_1", float("nan")),
                        "sigma_ratio": root.get("sigma_ratio", float("nan")),
                        "source_path": model_payload.get("source_path", ""),
                        "cache_source_path": model_payload.get("cache_source_path", ""),
                        "omega": omega,
                        "omega_over_cutoff_1": ratio1,
                        "omega_over_cutoff_2": ratio2,
                        "max_cutoff_ratio": ratio_max,
                        "timo_basis_regime_arm1": regime1,
                        "timo_basis_regime_arm2": regime2,
                        "execution_status": execution_status,
                        "N_true_cached": payload.get("N_true", ""),
                        "first_failed_mode": payload.get("first_failed_mode", ""),
                        "first_failed_delta_f": payload.get("first_failed_delta_f", ""),
                        "prefix_guard_status": payload.get("prefix_guard_status", model_payload.get("guard_status", "")),
                        "prefix_guard_resolved_through": payload.get("prefix_guard_resolved_through", ""),
                        "full_K10_guard_status": payload.get("full_K10_guard_status", "full_K10_guard_resolved"),
                        "early_stop_used": payload.get("early_stop_used", False),
                        "early_stop_reason": payload.get("early_stop_reason", ""),
                        "strict_policy_effective": payload.get("strict_policy_effective", "not_attempted"),
                        "cache_provenance": payload.get("cache_provenance", ""),
                        "unresolved_reason": payload.get("unresolved_reason", ""),
                    }
                )
    rows.sort(key=lambda row: (str(row["case_id"]), str(row["model"]), int(row["sorted_index"])))
    workflow.write_csv(output_dir / "spectra_long.csv", rows, SPECTRA_FIELDS)
    return rows


def export_case_execution(
    manifest: Sequence[workflow.GridPoint],
    payloads: Mapping[str, Mapping[str, object]],
    output_dir: Path,
    *,
    envelope_only: bool,
) -> list[dict[str, object]]:
    saturated = {
        float(payload["point"]["epsilon_0"])  # type: ignore[index]
        for payload in payloads.values()
        if payload.get("N_true") == workflow.K_MAX and isinstance(payload.get("point"), Mapping)
    }
    rows: list[dict[str, object]] = []
    for point in manifest:
        payload = payloads.get(point.case_id)
        if payload is None:
            status = (
                "not_needed_after_envelope_saturation"
                if envelope_only and point.epsilon_0 in saturated
                else "not_attempted"
            )
            rows.append({"case_id": point.case_id, "execution_status": status})
            continue
        strict = payload.get("strict", {})
        rows.append(
            {
                "case_id": point.case_id,
                "execution_status": payload.get("execution_status", "attempted_unresolved"),
                "N_true_cached": payload.get("N_true", ""),
                "first_failed_mode": payload.get("first_failed_mode", ""),
                "first_failed_delta_f": payload.get("first_failed_delta_f", ""),
                "prefix_guard_status": payload.get("prefix_guard_status", ""),
                "prefix_guard_resolved_through": payload.get("prefix_guard_resolved_through", ""),
                "full_K10_guard_status": payload.get("full_K10_guard_status", ""),
                "early_stop_used": payload.get("early_stop_used", False),
                "early_stop_reason": payload.get("early_stop_reason", ""),
                "strict_verification_status": (
                    strict.get("status", "") if isinstance(strict, Mapping) else ""
                ),
                "strict_trigger_reasons": (
                    ";".join(str(value) for value in strict.get("trigger_reasons", ()))
                    if isinstance(strict, Mapping)
                    else ""
                ),
                "strict_policy_effective": payload.get("strict_policy_effective", "not_attempted"),
                "cache_provenance": payload.get("cache_provenance", ""),
                "unresolved_reason": payload.get("unresolved_reason", ""),
            }
        )
    workflow.write_csv(output_dir / "case_execution.csv", rows, CASE_EXECUTION_FIELDS)
    return rows


def _point_lookup(
    manifest: Sequence[workflow.GridPoint],
    epsilon_0: float,
    beta_deg: float,
    mu: float,
    eta: float,
    *,
    regression_label: str = "",
) -> workflow.GridPoint:
    matches = [
        point
        for point in manifest
        if point.epsilon_0 == float(epsilon_0)
        and point.beta_deg == float(beta_deg)
        and point.mu == float(mu)
        and point.eta == float(eta)
        and (not regression_label or point.regression_label == regression_label)
    ]
    if len(matches) != 1:
        raise RuntimeError(f"sanity point lookup returned {len(matches)} rows")
    return matches[0]


def _sanity_manifest(manifest: Sequence[workflow.GridPoint]) -> list[tuple[str, workflow.GridPoint]]:
    return [
        ("exact_regression_S3_12", _point_lookup(manifest, 0.029408510742187498, 90, 0.7, 0, regression_label="S3_12")),
        ("exact_regression_S3_14", _point_lookup(manifest, 0.024798906738281248, 45, 0.5, -0.1, regression_label="S3_14")),
        ("mixed_regime", _point_lookup(manifest, 0.02, 45, 0.5, 0.0)),
        ("cutoff_limit", _point_lookup(manifest, 0.05, 90, 0.0, 0.0)),
        ("two_trigonometric", _point_lookup(manifest, 0.06, 90, 0.9, -0.5)),
        ("EB_close_root", _point_lookup(manifest, 0.05, 0, 0.7, -0.5)),
    ]


def _basis_regimes(point: workflow.GridPoint, payload: Mapping[str, object]) -> tuple[str, ...]:
    model = _accepted_model(payload, complete.MODEL_TIMO)
    if model is None:
        return ()
    sections = (
        timo.section_from_epsilon_tau(point.epsilon_0, point.tau1),
        timo.section_from_epsilon_tau(point.epsilon_0, point.tau2),
    )
    regimes: set[str] = set()
    for root in model.get("roots", ()):  # type: ignore[union-attr]
        value = float(root["Lambda"])
        regimes.update(timo.timo_basis(value, point.epsilon_0, section).regime for section in sections)
    return tuple(sorted(regimes))


def run_readiness_sanity(
    manifest: Sequence[workflow.GridPoint],
    *,
    output_dir: Path,
    strategy: str,
    strict_policy: str,
    force: bool,
    log_lines: list[str],
    runtime_records: dict[tuple[str, str], dict[str, object]],
) -> tuple[list[dict[str, object]], dict[str, int], dict[str, float]]:
    declared = _sanity_manifest(manifest)
    points = [point for _role, point in declared]
    payloads, counters, timings = solve_selected_points_prefix(
        points,
        output_dir=output_dir,
        reuse_cache=True,
        force=force,
        strategy=strategy,
        strict_policy=strict_policy,
        envelope_only=False,
        log_lines=log_lines,
        full_k10=True,
        force_audit_all=True,
        phase="readiness_sanity",
        runtime_records=runtime_records,
        cache_namespace="sanity",
        strict_cache_namespace="sanity_strict_branch",
    )
    expected = {
        "S3_12": (4, workflow.REGRESSION_POINTS[0][6]),
        "S3_14": (4, workflow.REGRESSION_POINTS[1][6]),
    }
    rows: list[dict[str, object]] = []
    for role, point in declared:
        payload = payloads[point.case_id]
        deltas = payload.get("deltas", {})
        delta5 = float(deltas.get("5", float("nan"))) if isinstance(deltas, Mapping) else float("nan")
        strict = payload.get("strict", {})
        strict_status = str(strict.get("status", "")) if isinstance(strict, Mapping) else ""
        passed = payload.get("execution_status") == "resolved_full_K10" and strict_status == "full_strict_pass"
        if point.regression_label:
            expected_n, expected_delta = expected[point.regression_label]
            passed = passed and payload.get("N_true") == expected_n and abs(delta5 - expected_delta) <= 1.0e-9
        regimes = _basis_regimes(point, payload)
        if role == "mixed_regime":
            passed = passed and timo.TIMO_REGIME_MIXED in regimes
        if role == "two_trigonometric":
            passed = passed and timo.TIMO_REGIME_TWO_TRIG in regimes
        if role == "cutoff_limit":
            sections = (
                timo.section_from_epsilon_tau(point.epsilon_0, point.tau1),
                timo.section_from_epsilon_tau(point.epsilon_0, point.tau2),
            )
            cutoff_regimes = []
            for section in sections:
                lambda_cutoff = math.sqrt(timo.omega_cutoff(section) / point.epsilon_0)
                cutoff_regimes.append(timo.timo_basis(lambda_cutoff, point.epsilon_0, section).regime)
            passed = passed and all(value == timo.TIMO_REGIME_CUTOFF for value in cutoff_regimes)
        rows.append(
            {
                "sanity_role": role,
                "case_id": point.case_id,
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "execution_status": payload.get("execution_status"),
                "N_true": payload.get("N_true", ""),
                "delta_f_5": delta5,
                "strict_verification_status": strict_status,
                "sanity_status": "PASS" if passed else "FAIL",
            }
        )
    workflow.write_csv(output_dir / "readiness_sanity.csv", rows, SANITY_FIELDS)
    if any(row["sanity_status"] != "PASS" for row in rows):
        raise RuntimeError("readiness sanity or exact regression failed; coarse grid was not started")
    return rows, counters, timings


def _failure_stratum(payload: Mapping[str, object]) -> str:
    value = payload.get("first_failed_mode")
    if value in {None, ""}:
        return "N_true_10"
    mode = int(value)
    if mode <= 4:
        return "early_failure"
    if mode <= 7:
        return "middle_prefix"
    return "late_prefix"


def _cutoff_warning(point: workflow.GridPoint, payload: Mapping[str, object]) -> bool:
    model = _accepted_model(payload, complete.MODEL_TIMO)
    if model is None:
        return False
    sections = (
        timo.section_from_epsilon_tau(point.epsilon_0, point.tau1),
        timo.section_from_epsilon_tau(point.epsilon_0, point.tau2),
    )
    ratios: list[float] = []
    regimes: set[str] = set()
    for root in model.get("roots", ()):  # type: ignore[union-attr]
        value = float(root["Lambda"])
        omega = timo.project_omega(value, point.epsilon_0)
        for section in sections:
            ratios.append(omega / timo.omega_cutoff(section))
            regimes.add(timo.timo_basis(value, point.epsilon_0, section).regime)
    return (
        timo.TIMO_REGIME_CUTOFF in regimes
        or (timo.TIMO_REGIME_MIXED in regimes and timo.TIMO_REGIME_TWO_TRIG in regimes)
        or any(0.95 <= value <= 1.05 for value in ratios)
        or bool(model.get("root12_boundary_warning"))
    )


def _cluster_or_order_warning(payload: Mapping[str, object]) -> bool:
    reasons = ";".join(str(value) for value in payload.get("strict", {}).get("trigger_reasons", ())) if isinstance(payload.get("strict"), Mapping) else ""
    if "cluster" in reasons or "order" in reasons:
        return True
    for model in complete.SUPPORTED_MODELS:
        item = _accepted_model(payload, model)
        if item is None:
            continue
        for root in item.get("roots", ()):  # type: ignore[union-attr]
            if int(root.get("cluster_size", 1)) > 1 or str(root.get("multiplicity_status", "simple_root")) != "simple_root":
                return True
    return False


def _sample_controls(
    rows: Sequence[Mapping[str, object]],
    *,
    fraction: float = 0.05,
    seed: int = CONTROL_SAMPLE_SEED,
    population_size: int | None = None,
) -> set[str]:
    """Choose >=5% deterministically while covering every observed diagnostic level."""

    target = max(1, math.ceil(float(fraction) * int(population_size or len(rows))))
    target = min(target, len(rows))
    fields = (
        "epsilon_0",
        "beta_deg",
        "mu",
        "eta",
        "N_true",
        "thin_0p1_flag",
        "basis_regimes",
        "first_failure_stratum",
    )

    def rank(row: Mapping[str, object]) -> str:
        return hashlib.sha256(f"{seed}|{row['case_id']}".encode("utf-8")).hexdigest()

    selected: set[str] = set()
    for field in fields:
        values = sorted({str(row.get(field, "")) for row in rows})
        for value in values:
            eligible = [row for row in rows if str(row.get(field, "")) == value]
            selected.add(str(min(eligible, key=rank)["case_id"]))
    for row in sorted(rows, key=rank):
        if len(selected) >= target:
            break
        selected.add(str(row["case_id"]))
    return selected


def build_full_control_manifest(
    manifest: Sequence[workflow.GridPoint],
    payloads: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    by_epsilon: dict[float, int] = {}
    for point in manifest:
        payload = payloads.get(point.case_id)
        if point.claim_eligible and payload is not None and payload.get("case_status") == "resolved":
            value = int(payload["N_true"])
            by_epsilon[point.epsilon_0] = max(value, by_epsilon.get(point.epsilon_0, -1))
    sample_rows: list[dict[str, object]] = []
    for point in manifest:
        payload = payloads.get(point.case_id)
        if not point.claim_eligible or payload is None or payload.get("case_status") != "resolved":
            continue
        sample_rows.append(
            {
                "case_id": point.case_id,
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "N_true": payload.get("N_true"),
                "thin_0p1_flag": point.thin_0p1_flag,
                "basis_regimes": ";".join(_basis_regimes(point, payload)),
                "first_failure_stratum": _failure_stratum(payload),
            }
        )
    sampled = _sample_controls(sample_rows, population_size=sum(point.claim_eligible for point in manifest))
    rows: list[dict[str, object]] = []
    for point in manifest:
        payload = payloads.get(point.case_id)
        if payload is None:
            continue
        reasons: list[str] = []
        if point.regression_label in {"S3_12", "S3_14"}:
            reasons.append("exact_regression_point")
        if point.claim_eligible and payload.get("case_status") == "resolved" and int(payload["N_true"]) == by_epsilon.get(point.epsilon_0):
            reasons.extend(("defines_N_up_raw", "argmax_case"))
        first_delta = payload.get("first_failed_delta_f")
        if first_delta not in {None, ""} and 0.095 <= float(first_delta) <= 0.105:
            reasons.append("near_first_failed_threshold")
        if payload.get("strict_policy_effective") == "full":
            reasons.append("auto_escalated_full_strict")
        if _cluster_or_order_warning(payload):
            reasons.append("cluster_or_order_warning")
        if _cutoff_warning(point, payload):
            reasons.append("cutoff_warning")
        if payload.get("case_status") != "resolved":
            reasons.append("prefix_unresolved_or_error")
        if point.case_id in sampled:
            reasons.append("stratified_five_percent_control")
        if not reasons:
            continue
        rows.append(
            {
                "case_id": point.case_id,
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "N_true_prefix": payload.get("N_true", ""),
                "first_failed_mode_prefix": payload.get("first_failed_mode", ""),
                "first_failed_delta_f_prefix": payload.get("first_failed_delta_f", ""),
                "thin_0p1_flag": point.thin_0p1_flag,
                "basis_regimes": ";".join(_basis_regimes(point, payload)),
                "first_failure_stratum": _failure_stratum(payload),
                "selection_reasons": ";".join(dict.fromkeys(reasons)),
            }
        )
    rows.sort(key=lambda row: (float(row["epsilon_0"]), float(row["beta_deg"]), float(row["mu"]), float(row["eta"]), str(row["case_id"])))
    return rows


def _root_rows(payload: Mapping[str, object], model: str) -> dict[int, Mapping[str, object]]:
    item = _accepted_model(payload, model)
    if item is None:
        return {}
    return {int(root["sorted_index"]): root for root in item.get("roots", ())}  # type: ignore[union-attr]


def _compare_prefix_full(
    prefix_payload: Mapping[str, object],
    full_payload: Mapping[str, object],
) -> dict[str, object]:
    prefix_resolved = prefix_payload.get("case_status") == "resolved"
    full_resolved = full_payload.get("case_status") == "resolved"
    first = prefix_payload.get("first_failed_mode")
    required_guard = workflow.ROOT_GUARD_INDEX if first in {None, ""} else min(workflow.ROOT_GUARD_INDEX, int(first) + 1)
    tolerance = workflow.primary_settings().root_match_tol
    differences: list[float] = []
    root_agreement = True
    cluster_agreement = True
    compared: dict[str, int] = {}
    reasons: list[str] = []
    for model in complete.SUPPORTED_MODELS:
        left = _root_rows(prefix_payload, model)
        right = _root_rows(full_payload, model)
        token = _model_token(model)
        compared[token] = 0
        for index in range(1, required_guard + 1):
            if index not in left or index not in right:
                root_agreement = False
                reasons.append(f"{token}_missing_root_{index}")
                continue
            compared[token] += 1
            difference = abs(float(left[index]["Lambda"]) - float(right[index]["Lambda"]))
            differences.append(difference)
            if difference > tolerance:
                root_agreement = False
                reasons.append(f"{token}_root_{index}_difference")
            identity_fields = ("cluster_id", "cluster_member_index", "cluster_size", "detected_nullity", "multiplicity_status")
            if any(left[index].get(field) != right[index].get(field) for field in identity_fields):
                cluster_agreement = False
                reasons.append(f"{token}_cluster_{index}_difference")
    n_agreement = prefix_resolved and full_resolved and prefix_payload.get("N_true") == full_payload.get("N_true")
    first_agreement = prefix_resolved and full_resolved and prefix_payload.get("first_failed_mode") == full_payload.get("first_failed_mode")
    execution_agreement = prefix_resolved == full_resolved
    guard_agreement = prefix_payload.get("prefix_guard_status") == "prefix_guard_resolved" and full_payload.get("full_K10_guard_status") == "full_K10_guard_resolved"
    strict = full_payload.get("strict", {})
    strict_pass = isinstance(strict, Mapping) and strict.get("status") == "full_strict_pass"
    if not n_agreement:
        reasons.append("N_true_disagreement")
    if not first_agreement:
        reasons.append("first_failed_mode_disagreement")
    if not execution_agreement:
        reasons.append("execution_status_disagreement")
    if not guard_agreement:
        reasons.append("guard_disagreement")
    if not strict_pass:
        reasons.append("full_strict_not_passed")
    passed = all((prefix_resolved, full_resolved, n_agreement, first_agreement, root_agreement, cluster_agreement, guard_agreement, execution_agreement, strict_pass))
    return {
        "prefix_execution_status": prefix_payload.get("execution_status", ""),
        "full_execution_status": full_payload.get("execution_status", ""),
        "prefix_N_true": prefix_payload.get("N_true", ""),
        "full_N_true": full_payload.get("N_true", ""),
        "N_true_agreement": n_agreement,
        "prefix_first_failed_mode": prefix_payload.get("first_failed_mode", ""),
        "full_first_failed_mode": full_payload.get("first_failed_mode", ""),
        "first_failed_mode_agreement": first_agreement,
        "compared_root_count_EB": compared.get("EB", 0),
        "compared_root_count_Timo": compared.get("Timo", 0),
        "max_root_absolute_difference": max(differences) if differences else "",
        "root_agreement": root_agreement,
        "cluster_identity_agreement": cluster_agreement,
        "prefix_guard_agreement": guard_agreement,
        "execution_status_agreement": execution_agreement,
        "full_strict_pass": strict_pass,
        "comparison_status": "PASS" if passed else "FAIL",
        "disagreement_reason": ";".join(dict.fromkeys(reasons)),
    }


def run_full_k10_controls(
    manifest: Sequence[workflow.GridPoint],
    prefix_payloads: Mapping[str, Mapping[str, object]],
    control_manifest: Sequence[Mapping[str, object]],
    *,
    output_dir: Path,
    strategy: str,
    strict_policy: str,
    force: bool,
    log_lines: list[str],
    runtime_records: dict[tuple[str, str], dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, int], dict[str, float]]:
    by_id = {point.case_id: point for point in manifest}
    points = [by_id[str(row["case_id"])] for row in control_manifest]
    full_payloads, counters, timings = solve_selected_points_prefix(
        points,
        output_dir=output_dir,
        reuse_cache=True,
        force=force,
        strategy=strategy,
        strict_policy=strict_policy,
        envelope_only=False,
        log_lines=log_lines,
        full_k10=True,
        force_audit_all=True,
        phase="full_k10_controls",
        runtime_records=runtime_records,
    )
    partial_cache = prefix.PartialPointCache(output_dir / "cache" / "prefix", reuse_cache=True, force=False)
    strict_cache = branch.BranchContinuationCache(
        output_dir / "cache" / "prefix_strict_branch" / strategy / strict_policy,
        reuse_cache=True,
        force_recompute=False,
        verification_scope="force_strict_verification",
    )
    callback = _make_prefix_full_strict_callback(strict_cache, workflow.strict_settings(), counters)
    results: list[dict[str, object]] = []
    comparisons: list[dict[str, object]] = []
    manual_strict_seconds_total = 0.0
    manifest_by_id = {str(row["case_id"]): row for row in control_manifest}
    for point in points:
        full_payload = full_payloads[point.case_id]
        strict = full_payload.get("strict", {})
        if not (isinstance(strict, Mapping) and strict.get("status") == "full_strict_pass"):
            raw = partial_cache.load(point, strategy=strategy, strict_policy=strict_policy)
            if raw is not None and raw.get("execution_status") == "resolved_full_K10":
                states = raw.get("models", {})
                sessions = {
                    model: prefix.IncrementalModelSession(
                        model,
                        point.geometry,
                        state=states.get(model) if isinstance(states, Mapping) else None,  # type: ignore[arg-type,union-attr]
                    )
                    for model in complete.SUPPORTED_MODELS
                }
                manual_started = time.perf_counter()
                audit = callback(point, sessions, workflow.ROOT_GUARD_INDEX)
                manual_seconds = time.perf_counter() - manual_started
                manual_strict_seconds_total += manual_seconds
                raw["strict_policy_effective"] = "full"
                raw["strict_trigger_reasons"] = list(dict.fromkeys([*raw.get("strict_trigger_reasons", ()), "deterministic_control_case"]))
                raw["strict_verification_status"] = audit.get("status", "full_strict_failed")
                raw["strict_details"] = dict(audit)
                if audit.get("status") != "full_strict_pass":
                    raw["execution_status"] = "attempted_unresolved"
                    raw["unresolved_reason"] = "full_control_strict_failure_or_disagreement"
                partial_cache.save(point, raw, strategy=strategy, strict_policy=strict_policy)
                full_payload = _prefix_case_payload(
                    point,
                    raw,
                    partial_cache.path(point, strategy=strategy, strict_policy=strict_policy),
                )
                full_payloads[point.case_id] = full_payload
                runtime = runtime_records.get(("full_k10_controls", point.case_id))
                if runtime is not None:
                    runtime["wall_seconds"] = float(runtime.get("wall_seconds", 0.0)) + manual_seconds
                    runtime["strict_seconds_cumulative"] = _cumulative_strict_seconds(raw)
        comparison = {"case_id": point.case_id, **_compare_prefix_full(prefix_payloads[point.case_id], full_payload)}
        comparisons.append(comparison)
        declared = dict(manifest_by_id[point.case_id])
        eb = _root_rows(full_payload, complete.MODEL_EB)
        tm = _root_rows(full_payload, complete.MODEL_TIMO)
        strict = full_payload.get("strict", {})
        runtime = runtime_records.get(("full_k10_controls", point.case_id), {})
        result: dict[str, object] = {
            **declared,
            "full_execution_status": full_payload.get("execution_status", ""),
            "full_N_true": full_payload.get("N_true", ""),
            "full_first_failed_mode": full_payload.get("first_failed_mode", ""),
            "full_strict_verification_status": strict.get("status", "") if isinstance(strict, Mapping) else "",
            "full_K10_guard_status": full_payload.get("full_K10_guard_status", ""),
            "wall_seconds": runtime.get("wall_seconds", ""),
            "cache_hit_current_invocation": runtime.get("cache_hit_current_invocation", ""),
        }
        for index in range(1, 12):
            result[f"full_EB_root_{index}"] = eb.get(index, {}).get("Lambda", "")
            result[f"full_Timo_root_{index}"] = tm.get(index, {}).get("Lambda", "")
        for index in range(1, 11):
            if index in eb and index in tm:
                result[f"full_delta_f_{index}"] = workflow.squared_frequency_delta(float(eb[index]["Lambda"]), float(tm[index]["Lambda"]))
            else:
                result[f"full_delta_f_{index}"] = ""
        results.append(result)
    results.sort(key=lambda row: str(row["case_id"]))
    comparisons.sort(key=lambda row: str(row["case_id"]))
    workflow.write_csv(output_dir / "full_k10_control_results.csv", results, CONTROL_RESULT_FIELDS)
    workflow.write_csv(output_dir / "prefix_full_comparison.csv", comparisons, CONTROL_COMPARISON_FIELDS)
    workflow.write_csv(
        output_dir / "runtime_by_case.csv",
        sorted(
            runtime_records.values(),
            key=lambda row: (str(row.get("phase", "")), int(row.get("ordinal", 0)), str(row.get("case_id", ""))),
        ),
        RUNTIME_FIELDS,
    )
    timings["full_k10_controls_manual_strict_seconds"] = manual_strict_seconds_total
    return results, comparisons, counters, timings


def postprocess_only(output_dir: Path) -> dict[str, object]:
    spectra = output_dir / "spectra_long.csv"
    manifest = output_dir / "grid_manifest.csv"
    if not manifest.exists() or not spectra.exists():
        raise FileNotFoundError("postprocess-only requires existing grid_manifest.csv and spectra_long.csv")
    from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: WPS433
        analyze_article_epsilon_upper_envelope as analyzer,
    )

    return analyzer.analyze(output_dir)


def main_pass_postprocess_only(
    output_dir: Path,
    *,
    defer_case_list: Path | None = None,
) -> dict[str, object]:
    from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: WPS433
        analyze_article_epsilon_upper_envelope as analyzer,
    )

    return analyzer.analyze_main_pass(
        output_dir,
        defer_case_list=repo_path(defer_case_list) if defer_case_list else None,
    )


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    initial_dirty = _git_text("status", "--short")
    initial_repository = {
        "cwd": str(Path.cwd()),
        "git_root": _git_text("rev-parse", "--show-toplevel"),
        "branch": _git_text("branch", "--show-current"),
        "HEAD": _git_text("rev-parse", "HEAD"),
        "dirty_git_status_short": initial_dirty.splitlines() if initial_dirty else [],
    }
    default_output = workflow.SMOKE_OUTPUT_DIR if args.smoke else COARSE_GRID_OUTPUT_DIR
    output_dir = repo_path(args.output_dir or default_output)
    if args.reconcile_family_local_repair_shadow:
        result = family_reconciliation.reconcile(
            output_dir,
            shadow_dir=(repo_path(args.shadow_dir) if args.shadow_dir else None),
            promotion_policy=args.promotion_policy,
        )
        result["root_calculations"] = 0
        return result
    if args.compact_only:
        compact_dir = (
            repo_path(args.compact_certificate_dir)
            if args.compact_certificate_dir
            else output_dir / "compact_point_certificates_v1"
        )
        run = compact_certificates.build_compact_certificates(
            output_dir, compact_dir=compact_dir,
        )
        return {
            "output_dir": run.output_dir,
            "certificate_count": run.certificate_count,
            "converted_count": run.converted_count,
            "cache_hit_count": run.cache_hit_count,
            "failure_count": run.failure_count,
            "peak_rss_bytes": run.peak_rss,
            "wall_seconds": run.wall_seconds,
            "root_calculations": 0,
            "matrix_evaluator_calls": 0,
            "local_repair_calls": 0,
        }
    if args.family_post_stage_only:
        from scripts.lib import article_epsilon_compact_poststage as compact_poststage  # noqa: WPS433

        compact_dir = (
            repo_path(args.compact_certificate_dir)
            if args.compact_certificate_dir
            else output_dir / "compact_point_certificates_v1"
        )
        result = compact_poststage.run_compact_family_poststage(
            output_dir, compact_dir=compact_dir,
        )
        result["root_calculations"] = 0
        return result
    if args.epsilon_005_target_diagnostics_only:
        from scripts.lib import article_epsilon_targeted_resolution as targeted  # noqa: WPS433

        result = targeted.prepare_target_diagnostics(output_dir)
        result["root_calculations"] = 0
        return result
    if args.epsilon_005_targeted_resolution:
        from scripts.lib import article_epsilon_targeted_resolution as targeted  # noqa: WPS433

        result = targeted.run_targeted_resolution(output_dir)
        result["root_calculations"] = 0
        return result
    if args.postprocess_only:
        if args.family_inventory_policy == "shadow":
            result = family_integration.run_shadow(output_dir)
            result["root_calculations"] = int(result.get("root_calculations", 0))
            return result
        result = (
            main_pass_postprocess_only(
                output_dir, defer_case_list=args.defer_case_list
            )
            if args.main_pass_only
            else postprocess_only(output_dir)
        )
        result["root_calculations"] = 0
        return result
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(parents=True, exist_ok=True)
    (output_dir / "cache").mkdir(parents=True, exist_ok=True)

    manifest_started = time.perf_counter()
    manifest = workflow.build_manifest(smoke=args.smoke)
    if not args.smoke:
        _validate_full_manifest(manifest)
    manifest_path = output_dir / "grid_manifest.csv"
    workflow.write_csv(
        manifest_path,
        [point.manifest_row() for point in manifest],
        workflow.MANIFEST_FIELDS,
    )
    manifest_seconds = time.perf_counter() - manifest_started

    selected = workflow.select_points(
        manifest,
        base_only=args.base_only,
        low_angle_only=args.low_angle_only,
        regressions_only=args.regressions_only,
    )
    full_contract_run = (
        not args.smoke
        and not args.base_only
        and not args.low_angle_only
        and not args.regressions_only
    )
    if args.readiness_sanity_only and not full_contract_run:
        raise RuntimeError("--readiness-sanity-only requires the full non-smoke manifest")
    if full_contract_run and (
        not args.prefix_until_failure
        or args.prefix_strategy != "paired"
        or args.envelope_only
    ):
        raise RuntimeError(
            "coarse_grid_v1 requires --prefix-until-failure --prefix-strategy paired "
            "and forbids --envelope-only"
        )
    if full_contract_run and args.main_pass_only and args.strict_policy != "main-pass":
        raise RuntimeError("cheap --main-pass-only requires --strict-policy main-pass")
    if full_contract_run and not args.main_pass_only and args.strict_policy != "auto":
        raise RuntimeError("the legacy full-contract workflow requires --strict-policy auto")
    if args.main_pass_only and not full_contract_run:
        raise RuntimeError("--main-pass-only requires the complete 1554-point manifest")
    configuration = workflow.solver_configuration()
    metadata = _pre_run_metadata(
        args=args,
        argv=argv,
        manifest_path=manifest_path,
        manifest=manifest,
        started_utc=started_utc,
    )
    metadata["repository"] = initial_repository
    metadata["selected_point_count"] = len(selected)
    metadata["output_dir"] = str(output_dir.relative_to(REPO_ROOT)).replace("\\", "/")
    metadata["run_state"] = "pre_run_ready"
    metadata["timings_seconds"] = {"manifest": manifest_seconds}
    workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
    log_lines = [
        f"started_utc={started_utc}",
        f"runner_version={COARSE_GRID_RUNNER_VERSION}",
        f"workflow_version={workflow.WORKFLOW_VERSION}",
        f"output_dir={output_dir}",
        f"manifest_count={len(manifest)} selected_count={len(selected)}",
        f"selection=base_only:{args.base_only},low_angle_only:{args.low_angle_only},regressions_only:{args.regressions_only}",
    ]
    runtime_records: dict[tuple[str, str], dict[str, object]] = {}
    runtime_path = output_dir / "runtime_by_case.csv"
    if runtime_path.exists() and args.reuse_cache and not args.force:
        for row in workflow.read_csv(runtime_path):
            runtime_records[(row.get("phase", ""), row.get("case_id", ""))] = dict(row)
    aggregate_counters: dict[str, int] = defaultdict(int)
    solve_timings: dict[str, float] = {}

    deferred_rows: list[dict[str, str]] = []
    deferred_case_ids: set[str] = set()
    progress_context: dict[str, int] | None = None
    if args.main_pass_only:
        effective_defer_path = args.defer_case_list
        default_current = output_dir / "deferred_complex_cases_current.csv"
        if effective_defer_path is None and default_current.exists():
            effective_defer_path = default_current
        deferred_rows = _deferred_case_rows(effective_defer_path, manifest)
        deferred_case_ids = {str(row["case_id"]) for row in deferred_rows}
        existing_statuses = _existing_prefix_records(
            manifest,
            output_dir,
            strategy=args.prefix_strategy,
            strict_policy=args.strict_policy,
        )
        initial_resolved = sum(
            status in {"resolved_prefix_early_stop", "resolved_full_K10"}
            for status in existing_statuses.values()
        )
        initial_unresolved = sum(
            status in {"attempted_unresolved", "attempted_error"}
            for status in existing_statuses.values()
        )
        initial_interrupted = sum(
            status == "interrupted_incomplete" for status in existing_statuses.values()
        )
        initial_deferred_expensive = sum(
            status == "deferred_expensive_strict"
            for status in existing_statuses.values()
        )
        deferred_status_count = len(
            {
                case_id
                for case_id in deferred_case_ids
                if existing_statuses.get(case_id) == "not_attempted"
            }
        ) + initial_deferred_expensive
        initial_completed = initial_resolved + initial_unresolved + initial_deferred_expensive
        selected = _main_pass_selection(selected, existing_statuses, deferred_case_ids)
        reconciliation_resume_ids: set[str] | None = None
        if args.family_inventory_policy == "local-repair":
            reconciliation_dir = output_dir / family_reconciliation.OUTPUT_DIRECTORY_NAME
            if reconciliation_dir.exists():
                reconciliation_resume_ids = family_reconciliation.load_ready_resume_ids(
                    output_dir
                )
                selected = [
                    point for point in selected if point.case_id in reconciliation_resume_ids
                ]
        progress_context = {
            "manifest_total": len(manifest),
            "initial_completed": initial_completed,
            "initial_resolved": initial_resolved,
            "initial_unresolved": initial_unresolved,
            "interrupted": initial_interrupted,
            "deferred": deferred_status_count,
            "deferred_expensive_strict": initial_deferred_expensive,
            "workers": int(args.workers),
        }
        metadata["cached_regression_gate"] = _cached_s3_regression_gate(output_dir)
        metadata["main_pass_initial_status_counts"] = {
            "completed": initial_completed,
            "resolved": initial_resolved,
            "unresolved": initial_unresolved,
            "interrupted": initial_interrupted,
            "not_attempted": sum(
                status == "not_attempted" for status in existing_statuses.values()
            ),
            "deferred_complex": deferred_status_count,
            "deferred_expensive_strict": initial_deferred_expensive,
        }
        metadata["main_pass_initial_complex_case_ids"] = sorted(
            case_id
            for case_id, status in existing_statuses.items()
            if status in {
                "attempted_unresolved",
                "attempted_error",
                "interrupted_incomplete",
                "deferred_expensive_strict",
            }
        )
        metadata["deferred_case_ids"] = sorted(deferred_case_ids)
        metadata["main_pass_selected_count"] = len(selected)
        metadata["reconciliation_resume_filter"] = {
            "used": reconciliation_resume_ids is not None,
            "ready_case_count": (
                len(reconciliation_resume_ids)
                if reconciliation_resume_ids is not None else None
            ),
            "source": (
                str(
                    (output_dir / family_reconciliation.OUTPUT_DIRECTORY_NAME / "resume_plan.csv")
                    .relative_to(REPO_ROOT)
                ).replace("\\", "/")
                if reconciliation_resume_ids is not None else None
            ),
        }
        metadata["selected_point_count"] = len(selected)
        metadata["run_state"] = "main_pass_preflight_passed"
        workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)

    if full_contract_run and not args.main_pass_only:
        sanity_rows, sanity_counters, sanity_timings = run_readiness_sanity(
            manifest,
            output_dir=output_dir,
            strategy=args.prefix_strategy,
            strict_policy=args.strict_policy,
            force=args.force,
            log_lines=log_lines,
            runtime_records=runtime_records,
        )
        for key, value in sanity_counters.items():
            aggregate_counters[key] += value
        solve_timings.update(sanity_timings)
        metadata["readiness_sanity"] = {
            "status": "PASS",
            "case_count": len(sanity_rows),
            "cache_namespace": "cache/sanity",
        }
        metadata["run_state"] = "readiness_sanity_passed"
        workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
        if args.readiness_sanity_only:
            metadata["run_state"] = "readiness_sanity_only_complete"
            metadata["run_finished_utc"] = datetime.now(timezone.utc).isoformat()
            metadata["timings_seconds"] = {"manifest": manifest_seconds, **solve_timings}
            metadata["scientific_solver_file_sha256_at_finish"] = _scientific_solver_hashes()
            metadata["scientific_solver_files_unchanged_during_run"] = (
                metadata["scientific_solver_file_sha256_at_finish"]
                == metadata["scientific_solver_file_sha256"]
            )
            workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
            print(f"readiness sanity PASS: {len(sanity_rows)} cases; output={output_dir}")
            return {
                "output_dir": output_dir,
                "manifest_count": len(manifest),
                "sanity_count": len(sanity_rows),
                "readiness_sanity": "PASS",
                "root_calculations": sum(
                    sanity_counters.get(key, 0)
                    for key in ("EB_spectrum_solves", "Timo_spectrum_solves")
                ),
            }

    if args.prefix_until_failure:
        _selected_payloads, counters, prefix_timings = solve_selected_points_prefix(
            selected,
            output_dir=output_dir,
            reuse_cache=args.reuse_cache,
            force=args.force,
            strategy=args.prefix_strategy,
            strict_policy=args.strict_policy,
            envelope_only=args.envelope_only,
            log_lines=log_lines,
            phase="prefix_sweep",
            runtime_records=runtime_records,
            workers=args.workers,
            progress_context=progress_context,
        )
        for key, value in counters.items():
            aggregate_counters[key] += value
        solve_timings.update(prefix_timings)
    else:
        _selected_payloads, counters, full_timings = solve_selected_points(
            selected,
            output_dir=output_dir,
            reuse_cache=args.reuse_cache,
            force=args.force,
            configuration=configuration,
            log_lines=log_lines,
            strict_policy=args.strict_policy,
        )
        for key, value in counters.items():
            aggregate_counters[key] += value
        solve_timings.update(full_timings)

    export_started = time.perf_counter()
    available = (
        _all_available_prefix_payloads(
            manifest,
            output_dir,
            strategy=args.prefix_strategy,
            strict_policy=args.strict_policy,
        )
        if args.prefix_until_failure and args.main_pass_only
        else (
            _selected_payloads
            if args.prefix_until_failure
            else _all_available_payloads(manifest, output_dir, configuration)
        )
    )
    if args.main_pass_only:
        build_deferred_complex_cases_current(
            output_dir,
            seed_path=args.defer_case_list,
        )
    if args.prefix_until_failure:
        export_case_execution(
            manifest,
            available,
            output_dir,
            envelope_only=args.envelope_only,
        )
    spectra_rows = export_spectra_long(manifest, available, output_dir)
    export_seconds = time.perf_counter() - export_started
    resolved_available = sum(payload.get("case_status") == "resolved" for payload in available.values())
    unresolved_available = len(available) - resolved_available

    control_manifest: list[dict[str, object]] = []
    control_results: list[dict[str, object]] = []
    control_comparisons: list[dict[str, object]] = []
    if full_contract_run:
        control_manifest = build_full_control_manifest(manifest, available)
        workflow.write_csv(
            output_dir / "full_k10_control_manifest.csv",
            control_manifest,
            CONTROL_MANIFEST_FIELDS,
        )
        if not args.main_pass_only:
            control_results, control_comparisons, control_counters, control_timings = run_full_k10_controls(
                manifest,
                available,
                control_manifest,
                output_dir=output_dir,
                strategy=args.prefix_strategy,
                strict_policy=args.strict_policy,
                force=args.force,
                log_lines=log_lines,
                runtime_records=runtime_records,
            )
            for key, value in control_counters.items():
                aggregate_counters[key] += value
            solve_timings.update(control_timings)

    comparisons_pass = bool(control_comparisons) and all(
        row["comparison_status"] == "PASS" for row in control_comparisons
    )
    final_validation_gate = (
        "NOT_RUN_MAIN_PASS"
        if args.main_pass_only
        else
        "PASS"
        if (not full_contract_run or (resolved_available == len(manifest) and comparisons_pass))
        else "FAIL"
    )

    elapsed_before_postprocess = time.perf_counter() - started
    metadata.update(
        {
            "run_state": "root_solves_complete",
            "run_finished_utc": datetime.now(timezone.utc).isoformat(),
            "available_cached_point_count": len(available),
            "resolved_cached_point_count": resolved_available,
            "unresolved_cached_point_count": unresolved_available,
            "not_attempted_point_count": len(manifest) - len(available),
            "spectra_long_row_count": len(spectra_rows),
            "operation_counts_current_invocation": dict(aggregate_counters),
            "full_k10_control_count": len(control_manifest),
            "stratified_control_count": sum(
                "stratified_five_percent_control" in str(row["selection_reasons"])
                for row in control_manifest
            ),
            "full_k10_control_pass_count": sum(
                row["comparison_status"] == "PASS" for row in control_comparisons
            ),
            "final_validation_gate": final_validation_gate,
            "max_prefix_full_root_absolute_difference": max(
                (float(row["max_root_absolute_difference"]) for row in control_comparisons if row["max_root_absolute_difference"] != ""),
                default=None,
            ),
        }
    )
    metadata["timings_seconds"] = {
        "manifest": manifest_seconds,
        **solve_timings,
        "spectra_export": export_seconds,
        "pre_postprocess_total": elapsed_before_postprocess,
    }
    workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
    post_started = time.perf_counter()
    post_result = (
        main_pass_postprocess_only(
            output_dir,
            defer_case_list=args.defer_case_list,
        )
        if args.main_pass_only
        else postprocess_only(output_dir)
    )
    family_result: dict[str, object] | None = None
    if args.family_inventory_policy == "local-repair":
        # Never hand the family stage a process-wide collection of decompressed
        # solver traces.  Stream the immutable point caches into compact
        # certificates, then load only one small beta-family at a time.
        from scripts.lib import article_epsilon_compact_poststage as compact_poststage  # noqa: WPS433

        compact_run = compact_certificates.build_compact_certificates(output_dir)
        family_result = compact_poststage.run_compact_family_poststage(
            output_dir, compact_dir=compact_run.output_dir,
        )
        family_result["compaction"] = {
            "certificate_count": compact_run.certificate_count,
            "converted_count": compact_run.converted_count,
            "cache_hit_count": compact_run.cache_hit_count,
            "failure_count": compact_run.failure_count,
            "peak_rss_bytes": compact_run.peak_rss,
            "point_solver_calls": 0,
        }
    post_seconds = time.perf_counter() - post_started
    metadata["timings_seconds"]["postprocess"] = post_seconds  # type: ignore[index]
    metadata["timings_seconds"]["total"] = time.perf_counter() - started  # type: ignore[index]
    metadata["run_finished_utc"] = datetime.now(timezone.utc).isoformat()
    final_hashes = _scientific_solver_hashes()
    metadata["scientific_solver_file_sha256_at_finish"] = final_hashes
    metadata["scientific_solver_files_unchanged_during_run"] = (
        final_hashes == metadata["scientific_solver_file_sha256"]
    )
    final_dirty = _git_text("status", "--short")
    metadata["repository_final"] = {
        **{key: initial_repository[key] for key in ("cwd", "git_root", "branch", "HEAD")},
        "dirty_git_status_short": final_dirty.splitlines() if final_dirty else [],
    }
    metadata["run_state"] = "complete" if final_validation_gate == "PASS" else "complete_with_failed_validation"

    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ")
    log_lines.extend(
        [
            f"available={len(available)} resolved={resolved_available} unresolved={unresolved_available}",
            f"controls={len(control_manifest)} control_pass={sum(row['comparison_status'] == 'PASS' for row in control_comparisons)}",
            f"final_validation_gate={final_validation_gate}",
            f"counters={json.dumps(aggregate_counters, sort_keys=True)}",
            f"timings={json.dumps(metadata['timings_seconds'], sort_keys=True)}",
        ]
    )
    log_path = output_dir / "logs" / f"run_{run_id}.log"
    workflow.atomic_write_text(log_path, "\n".join(log_lines) + "\n")
    metadata["log_path"] = str(log_path.relative_to(REPO_ROOT)).replace("\\", "/")
    workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
    print(
        f"manifest={len(manifest)} selected={len(selected)} available={len(available)} "
        f"resolved={resolved_available} unresolved={unresolved_available} controls={len(control_manifest)} "
        f"validation={final_validation_gate}"
    )
    print(json.dumps(aggregate_counters, sort_keys=True))
    print(f"output={output_dir}")
    return {
        "output_dir": output_dir,
        "manifest_count": len(manifest),
        "selected_count": len(selected),
        "available_count": len(available),
        "resolved_count": resolved_available,
        "unresolved_count": unresolved_available,
        "control_count": len(control_manifest),
        "final_validation_gate": final_validation_gate,
        "counters": dict(aggregate_counters),
        "metadata": metadata,
        "postprocess": post_result,
        "family_inventory": family_result,
    }


if __name__ == "__main__":
    main()
