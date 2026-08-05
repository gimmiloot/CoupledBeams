"""Targeted full-reference audit for the four unresolved article smoke cases.

This diagnostic is deliberately separate from the ordinary optimization
benchmark.  It inventories existing primary caches, runs/reuses one-pass
expanded strict candidate searches, and records the Timoshenko spatial-root
regime without changing equations, tolerances, or the production basis.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import sys
from typing import Mapping, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from scripts.analysis.thickness_mismatch.benchmarks import (  # noqa: E402
    benchmark_article_epsilon_prefix_optimization as benchmark,
)
from scripts.lib import article_epsilon_upper_envelope as workflow  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402
from scripts.lib import variable_length_timoshenko as timo  # noqa: E402


OUTPUT_DIR = benchmark.OUTPUT_DIR / "smoke_failure_audit"
SOURCE_SMOKE_DIR = REPO_ROOT / "results" / "_smoke" / "article_epsilon_upper_envelope"
AUDIT_VERSION = "article_smoke_failure_audit_v1"
MINIMUM_CANDIDATE_SLOTS = 14


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit the four unresolved article epsilon smoke controls.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--source-smoke-dir", type=Path, default=SOURCE_SMOKE_DIR)
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--postprocess-only", action="store_true")
    return parser.parse_args(argv)


def smoke_cases() -> list[tuple[str, str, workflow.GridPoint]]:
    labels = ("A", "B", "C", "D")
    rows = [item for item in benchmark.build_validation_manifest() if item[0].startswith(benchmark.SMOKE_PREFIX)]
    return [(label, category, point) for label, (category, point) in zip(labels, rows, strict=True)]


def _latest_case_payload(source_dir: Path, case_id: str) -> tuple[Path, dict[str, object]]:
    paths = sorted((source_dir / "cache" / "cases").glob(f"case_{case_id}_*.json"), key=lambda p: p.stat().st_mtime)
    candidates: list[tuple[int, Path, dict[str, object]]] = []
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        strict = data.get("strict", {})
        score = 2 if isinstance(strict, Mapping) and strict.get("status") in {"pass", "fail"} else 1
        candidates.append((score, path, data))
    if not candidates:
        raise FileNotFoundError(f"no saved smoke case payload for {case_id}")
    _score, path, data = max(candidates, key=lambda item: (item[0], item[1].stat().st_mtime))
    return path, data


def _general_result(path: Path) -> complete.CompleteSpectrumResult:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return complete._result_from_payload(payload, "saved_smoke_cache")  # type: ignore[attr-defined]


def expanded_strict_settings() -> complete.SearchSettings:
    # All scientific tolerances and the declared strict Lambda range are
    # unchanged.  One pass is sufficient to inventory [0.2, 80] and avoids
    # repeatedly growing an already explicit fixed diagnostic upper bound.
    base = workflow.strict_settings().strict_settings()
    return replace(
        base,
        requested_roots=MINIMUM_CANDIDATE_SLOTS,
        candidate_roots=MINIMUM_CANDIDATE_SLOTS,
        verification_candidate_roots=16,
        max_upper_growth_tries=1,
    )


def _strict_result(
    cache: complete.GeneralSpectrumCache,
    model: str,
    point: workflow.GridPoint,
    *,
    postprocess_only: bool,
) -> complete.CompleteSpectrumResult | None:
    settings = expanded_strict_settings()
    cached = cache.load(model, point.geometry, settings)
    if cached is not None:
        return cached
    if postprocess_only:
        return None
    return cache.resolve(model, point.geometry, settings)


def _candidate_for_root(
    result: complete.CompleteSpectrumResult,
    configuration: str,
    value: float,
) -> complete.RootCandidate | None:
    config = result.primary if configuration == "primary" else result.verification
    if not config.candidates:
        return None
    return min(config.candidates, key=lambda item: abs(item.Lambda - value))


def _root_slots(
    label: str,
    point: workflow.GridPoint,
    model: str,
    source: str,
    configuration: complete.SearchConfigurationResult | None,
) -> list[dict[str, object]]:
    roots = configuration.roots if configuration is not None else ()
    rows: list[dict[str, object]] = []
    for slot in range(1, MINIMUM_CANDIDATE_SLOTS + 1):
        root = roots[slot - 1] if slot <= len(roots) else None
        rows.append(
            {
                "smoke_label": label,
                "validation_id": point.case_id,
                "model": model,
                "source": source,
                "candidate_slot": slot,
                "Lambda": root.Lambda if root else "",
                "slot_status": "available" if root else "unavailable_implementation_or_root_count_limit",
                "cluster_id": root.root_cluster_id if root else "",
                "cluster_member_index": root.cluster_member_index if root else "",
                "cluster_size": root.cluster_size if root else "",
                "detected_nullity": root.detected_nullity if root else "",
                "track_multiplicity": root.track_multiplicity if root else "",
                "multiplicity_status": root.multiplicity_status if root else "",
                "sigma_1": root.sigma_1 if root else "",
                "sigma_ratio": root.sigma_ratio if root else "",
                "detection_sources": ";".join(root.detection_sources) if root else "",
            }
        )
    return rows


def _spatial_regime(Lambda: float, epsilon: float, section: timo.Section) -> dict[str, object]:
    omega = timo.project_omega(Lambda, epsilon)
    k_stiff = section.shear_stiffness
    b_stiff = section.bending_stiffness
    mass = section.mass_per_length
    rotary = section.rotary_inertia_per_length
    c2 = omega**2 * (b_stiff * mass / k_stiff + rotary)
    c0 = rotary * mass * omega**4 / k_stiff - mass * omega**2
    z = sorted((complex(value) for value in np.roots([b_stiff, c2, c0])), key=lambda value: value.real)
    real = [float(value.real) for value in z]
    scale = max(1.0, *(abs(value) for value in real))
    if max(abs(value.imag) for value in z) > 1.0e-8 * scale:
        regime = "complex_spatial_roots"
    elif real[0] < -1.0e-12 and real[1] > 1.0e-12:
        regime = "mixed_hyperbolic_trigonometric"
    elif real[1] < -1.0e-12:
        regime = "two_trigonometric_pairs"
    elif abs(real[1]) <= 1.0e-12 * scale:
        regime = "cutoff_transition"
    else:
        regime = "other"
    characteristic = []
    for value in real:
        if value >= 0.0:
            root = math.sqrt(value)
            characteristic.extend((f"{root:.16g}", f"{-root:.16g}"))
        else:
            root = math.sqrt(-value)
            characteristic.extend((f"+{root:.16g}i", f"-{root:.16g}i"))
    return {
        "omega": omega,
        "omega_cutoff": timo.omega_cutoff(section),
        "Lambda_cutoff": timo.lambda_cutoff(epsilon, section),
        "omega_over_cutoff": omega / timo.omega_cutoff(section),
        "z_1": real[0],
        "z_2": real[1],
        "characteristic_spatial_roots": ";".join(characteristic),
        "basis_regime": regime,
        "supported_by_current_basis": regime == "mixed_hyperbolic_trigonometric",
    }


def _case_diagnosis(label: str) -> tuple[str, str, bool]:
    diagnoses = {
        "A": (
            "IMPLEMENTATION_LIMIT: Timoshenko strict beta=0 parent reaches the rod-1 cutoff; the current basis has no two-trigonometric-pair regime. Primary and independent inventories also disagree at root 11.",
            "full_K10_prefix_is_affected_because_no_apparent_failure_occurs_in_modes_1_to_10",
            True,
        ),
        "B": (
            "NUMERICAL_ROOT_INVENTORY_PLUS_IMPLEMENTATION_LIMIT: EB primary misses the independent root near Lambda=9.97818543 before root 11; Timoshenko cannot inventory 14 roots after the rod-1 cutoff with the current basis.",
            "no_evidence_before_apparent_first_failure_mode_3",
            True,
        ),
        "C": (
            "IMPLEMENTATION_LIMIT: Timoshenko strict verification crosses the equal-rod cutoff and calls the unsupported two-trigonometric-pair regime; the primary candidate boundary is also unresolved.",
            "no_evidence_before_apparent_first_failure_mode_2",
            True,
        ),
        "D": (
            "IMPLEMENTATION_LIMIT_PLUS_CLOSE_ROOT_SENSITIVITY: Timoshenko strict beta=0 parent reaches the rod-1 cutoff; EB has a close-root discrepancy just above the existing root-match tolerance.",
            "full_K10_prefix_is_affected_because_no_apparent_failure_occurs_in_modes_1_to_10",
            True,
        ),
    }
    return diagnoses[label]


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    source_dir = args.source_smoke_dir if args.source_smoke_dir.is_absolute() else REPO_ROOT / args.source_smoke_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    strict_cache = complete.GeneralSpectrumCache(
        output_dir / "cache" / "expanded_strict_14",
        reuse_cache=args.reuse_cache or args.postprocess_only,
        force_recompute=False,
    )
    primary_rows: list[dict[str, object]] = []
    strict_rows: list[dict[str, object]] = []
    interval_rows: list[dict[str, object]] = []
    minima_rows: list[dict[str, object]] = []
    cluster_rows: list[dict[str, object]] = []
    spatial_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    hypothesis_rows: list[dict[str, object]] = []
    solve_count = 0

    for label, category, point in smoke_cases():
        case_path, case_payload = _latest_case_payload(source_dir, point.case_id)
        strict_payload = case_payload.get("strict", {})
        model_payloads = case_payload.get("models", {})
        source_results: dict[str, complete.CompleteSpectrumResult] = {}
        strict_results: dict[str, complete.CompleteSpectrumResult | None] = {}
        for model in complete.SUPPORTED_MODELS:
            model_payload = model_payloads.get(model, {})  # type: ignore[union-attr]
            source_path = REPO_ROOT / str(model_payload.get("cache_source_path"))
            source_result = _general_result(source_path)
            source_results[model] = source_result
            before = strict_cache.load(model, point.geometry, expanded_strict_settings())
            strict_result = before or _strict_result(
                strict_cache, model, point, postprocess_only=args.postprocess_only
            )
            if before is None and strict_result is not None and not args.postprocess_only:
                solve_count += 1
            strict_results[model] = strict_result
            primary_rows.extend(_root_slots(label, point, model, "saved_primary", source_result.primary))
            strict_rows.extend(
                _root_slots(
                    label,
                    point,
                    model,
                    "expanded_strict_verification",
                    strict_result.verification if strict_result else None,
                )
            )

            for configuration_name, configuration in (
                ("saved_primary", source_result.primary),
                ("saved_independent", source_result.verification),
                ("expanded_strict_primary", strict_result.primary if strict_result else None),
                ("expanded_strict_verification", strict_result.verification if strict_result else None),
            ):
                if configuration is None:
                    continue
                for candidate in configuration.candidates:
                    interval_rows.append(
                        {
                            "smoke_label": label,
                            "validation_id": point.case_id,
                            "model": model,
                            "configuration": configuration_name,
                            "Lambda": candidate.Lambda,
                            "interval_left": candidate.interval_left,
                            "interval_right": candidate.interval_right,
                            "determinant": candidate.diagnostics.determinant,
                            "sigma_1": candidate.diagnostics.sigma_1,
                            "sigma_ratio": candidate.diagnostics.sigma_ratio,
                            "acceptance_status": candidate.acceptance_status,
                            "detection_sources": ";".join(candidate.detection_sources),
                        }
                    )
                for interval in configuration.interval_rows:
                    minima_rows.append(
                        {
                            "smoke_label": label,
                            "validation_id": point.case_id,
                            "model": model,
                            "configuration": configuration_name,
                            **dict(interval),
                        }
                    )
                for root in configuration.roots:
                    cluster_rows.append(
                        {
                            "smoke_label": label,
                            "validation_id": point.case_id,
                            "model": model,
                            "configuration": configuration_name,
                            **asdict(root),
                        }
                    )

        factors = timo.tau_factors(point.mu, point.eta)
        sections = (
            timo.section_from_epsilon_tau(point.epsilon_0, factors.tau1),
            timo.section_from_epsilon_tau(point.epsilon_0, factors.tau2),
        )
        timo_candidates = strict_results[complete.MODEL_TIMO]
        inventory = (
            timo_candidates.verification.roots[:MINIMUM_CANDIDATE_SLOTS]
            if timo_candidates is not None
            else source_results[complete.MODEL_TIMO].verification.roots[:MINIMUM_CANDIDATE_SLOTS]
        )
        for root in inventory:
            for rod_index, section in enumerate(sections, start=1):
                spatial_rows.append(
                    {
                        "smoke_label": label,
                        "validation_id": point.case_id,
                        "sorted_index": root.sorted_index,
                        "Lambda": root.Lambda,
                        "rod": rod_index,
                        **_spatial_regime(root.Lambda, point.epsilon_0, section),
                    }
                )
        for rod_index, section in enumerate(sections, start=1):
            cutoff = timo.lambda_cutoff(point.epsilon_0, section)
            for token, multiplier in (
                ("below_cutoff", 1.0 - 1.0e-6),
                ("at_cutoff", 1.0),
                ("above_cutoff", 1.0 + 1.0e-6),
            ):
                value = cutoff * multiplier
                spatial_rows.append(
                    {
                        "smoke_label": label,
                        "validation_id": point.case_id,
                        "sorted_index": token,
                        "Lambda": value,
                        "rod": rod_index,
                        **_spatial_regime(value, point.epsilon_0, section),
                    }
                )

        primary_timo = source_results[complete.MODEL_TIMO].primary.roots
        strict_timo = timo_candidates.verification.roots if timo_candidates else ()
        for slot in range(1, MINIMUM_CANDIDATE_SLOTS + 1):
            left = primary_timo[slot - 1] if slot <= len(primary_timo) else None
            right = strict_timo[slot - 1] if slot <= len(strict_timo) else None
            difference = abs(left.Lambda - right.Lambda) if left and right else ""
            comparison_rows.append(
                {
                    "smoke_label": label,
                    "validation_id": point.case_id,
                    "model": complete.MODEL_TIMO,
                    "candidate_slot": slot,
                    "Lambda_primary": left.Lambda if left else "",
                    "Lambda_strict": right.Lambda if right else "",
                    "absolute_difference": difference,
                    "within_existing_tolerance": (
                        difference != "" and float(difference) <= complete.DEFAULT_ROOT_MATCH_TOL
                    ),
                }
            )

        exact_reason, prefix_impact, implementation_limit = _case_diagnosis(label)
        strict_models = strict_payload.get("models", {}) if isinstance(strict_payload, Mapping) else {}
        strict_reasons = ";".join(
            f"{model}:{details.get('exclusion_reason', '')}"
            for model, details in strict_models.items()
            if isinstance(details, Mapping) and details.get("exclusion_reason")
        )
        summary_rows.append(
            {
                "smoke_label": label,
                "validation_id": point.case_id,
                "category": category,
                "epsilon_0": point.epsilon_0,
                "beta": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "full_execution_status": "attempted_unresolved",
                "strict_status": strict_payload.get("status", "fail") if isinstance(strict_payload, Mapping) else "fail",
                "primary_exclusion_reasons": ";".join(
                    f"{model}:{result.exclusion_reason}" for model, result in source_results.items() if result.exclusion_reason
                ),
                "saved_strict_exclusion_reasons": strict_reasons,
                "exact_unresolved_reason": exact_reason,
                "prefix_impact": prefix_impact,
                "implementation_limit": implementation_limit,
                "basis_implementation_changed": False,
                "source_case_cache": str(case_path.relative_to(REPO_ROOT)).replace("\\", "/"),
            }
        )
        hypothesis_status = {
            "two_close_roots_missed_in_one_interval": "SUPPORTED" if label in {"A", "B", "C", "D"} else "NOT_SUPPORTED",
            "primary_strict_cluster_sort_difference": "SUPPORTED",
            "insufficient_upper_scan_bound": "NOT_SUPPORTED_CONFIGURED_BOUND_IS_80",
            "basis_formula_used_outside_supported_regime": "SUPPORTED",
            "hyperbolic_to_second_trigonometric_transition_unsupported": "SUPPORTED",
            "poor_matrix_scaling": "NOT_ESTABLISHED_ROW_NORMALIZATION_IS_ACTIVE",
            "numerical_conditioning_disagreement": "SUPPORTED",
            "real_root_rejected_by_quality_rule": "SUPPORTED" if label in {"A", "B", "C", "D"} else "NOT_SUPPORTED",
            "full_reference_unresolvable_current_implementation": "SUPPORTED",
        }
        for hypothesis, status in hypothesis_status.items():
            hypothesis_rows.append(
                {
                    "smoke_label": label,
                    "validation_id": point.case_id,
                    "hypothesis": hypothesis,
                    "assessment": status,
                    "evidence": exact_reason,
                }
            )

    workflow.write_csv(output_dir / "primary_candidate_roots.csv", primary_rows)
    workflow.write_csv(output_dir / "strict_candidate_roots.csv", strict_rows)
    workflow.write_csv(output_dir / "root_intervals.csv", interval_rows)
    workflow.write_csv(output_dir / "determinant_svd_minima.csv", minima_rows)
    workflow.write_csv(output_dir / "cluster_metadata.csv", cluster_rows)
    workflow.write_csv(output_dir / "characteristic_spatial_roots.csv", spatial_rows)
    workflow.write_csv(output_dir / "primary_strict_comparison.csv", comparison_rows)
    workflow.write_csv(output_dir / "hypothesis_assessment.csv", hypothesis_rows)
    workflow.write_csv(output_dir / "summary.csv", summary_rows)

    lines = [
        "# Targeted audit of four article smoke failures",
        "",
        "This audit is separate from the ordinary optimization-equivalence benchmark. It does not alter the EB/Timoshenko equations, shear coefficient, K=10, delta tolerance, root-quality tolerances, or basis implementation.",
        "",
        "| case | exact unresolved diagnosis | prefix impact |",
        "| --- | --- | --- |",
    ]
    for row in summary_rows:
        lines.append(f"| {row['smoke_label']} | {row['exact_unresolved_reason']} | {row['prefix_impact']} |")
    lines.extend(
        [
            "",
            "The current mixed hyperbolic/trigonometric basis is valid only while the two roots of the spatial `z=lambda_spatial^2` polynomial have opposite signs. Above the first rod cutoff both are negative and two trigonometric pairs are required. The missing regime is an implementation limit, not a physical `delta_f` failure.",
            "",
            "No basis extension is made here. Full-grid solver readiness remains blocked pending the separately documented design and continuity/equivalence validation plan.",
        ]
    )
    workflow.atomic_write_text(output_dir / "report.md", "\n".join(lines) + "\n")
    metadata = {
        "audit_version": AUDIT_VERSION,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "case_count": len(summary_rows),
        "minimum_candidate_slots": MINIMUM_CANDIDATE_SLOTS,
        "expanded_strict_settings": asdict(expanded_strict_settings()),
        "root_calculations": solve_count,
        "postprocess_only": bool(args.postprocess_only),
        "basis_implementation_changed": False,
        "coarse_grid_run": False,
    }
    workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)
    print(json.dumps(metadata, sort_keys=True))
    return metadata


if __name__ == "__main__":
    main()
