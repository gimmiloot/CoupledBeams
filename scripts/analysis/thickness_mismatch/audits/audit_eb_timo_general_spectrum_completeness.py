from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import asdict, dataclass
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

from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    run_eb_epsilon_apriori_pilot as pilot,
)
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as legacy_roots,
)
from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_epsilon_apriori_pilot as pilot_analysis,
)
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402


DEFAULT_OUTPUT_DIR = Path("results") / "eb_timo_general_spectrum_completeness"
DEFAULT_CORRECTED_PILOT_DIR = Path("results") / "eb_epsilon_apriori_pilot_complete_spectrum_v1"
DEFAULT_LEGACY_PILOT_DIR = Path("results") / "eb_epsilon_apriori_pilot"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_timo_general_spectrum_completeness"
SMOKE_CORRECTED_PILOT_DIR = Path("results") / "_smoke" / "eb_epsilon_apriori_pilot_complete_spectrum_v1"

OUTPUT_NAMES = (
    "general_spectrum_root_records.csv",
    "general_spectrum_interval_audit.csv",
    "general_spectrum_candidate_union_audit.csv",
    "general_spectrum_completeness_summary.csv",
    "general_spectrum_primary_vs_verification.csv",
    "general_spectrum_factorized_oracle_comparison.csv",
    "general_spectrum_close_pair_stress_audit.csv",
    "general_spectrum_operation_counts.csv",
    "general_spectrum_exclusion_audit.csv",
    "eb_timo_general_spectrum_completeness_report.md",
)

R_CASES = {
    "R1": (0.0228017578125, 8.2931855, 8.2999557),
    "R2": (0.0159306640625, 9.9286581, 9.9298571),
    "R3": (0.015625, 14.177263, 14.179631),
}
STRESS_BETA_VALUES = (0.0, 0.1, 1.0, 5.0)
STRESS_EPSILON_OFFSETS = (-1.0e-6, 0.0, 1.0e-6)


ROOT_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "configuration",
    "Lambda", "sorted_index", "root_cluster_id", "cluster_member_index", "cluster_size",
    "detection_sources", "primary_status", "verification_status", "determinant", "sigma_1",
    "sigma_2", "sigma_ratio", "null_vector_pivot", "self_MAC", "raw_legacy_match",
    "factorized_oracle_match", "acceptance_status", "multiplicity_status", "notes",
]

INTERVAL_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "configuration",
    "Lambda_left", "Lambda_right", "determinant_sign_left", "determinant_sign_right",
    "sigma_left", "sigma_midpoint", "sigma_right", "minimum_sigma", "minimum_sigma_ratio",
    "subdivision_depth", "trigger_reasons", "local_minima_count", "resolution_status",
]

CANDIDATE_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "configuration", "Lambda",
    "detection_sources", "interval_left", "interval_right", "determinant", "sigma_1", "sigma_2",
    "sigma_ratio", "null_vector_pivot", "interior_minimum", "acceptance_status", "notes",
]

SUMMARY_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "requested_roots", "found_roots",
    "first12_resolved", "root11_available", "root12_available", "raw_root_count", "recovered_root_count",
    "close_cluster_count", "multiplicity_count", "independent_agreement", "spectrum_status",
    "exclusion_reason", "cache_status",
]

PRIMARY_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "sorted_index", "Lambda_primary",
    "Lambda_verification", "absolute_difference", "self_MAC", "multiplicity_agreement", "status",
]

ORACLE_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "sorted_index", "Lambda_general",
    "Lambda_oracle", "absolute_difference", "raw_legacy_present", "general_recovered", "status",
]

STRESS_FIELDS = [
    "case_id", "r_case", "epsilon_offset", "epsilon_0", "beta_deg", "model", "root_count",
    "pair_root_1", "pair_root_2", "minimum_gap", "both_roots_present", "grid_phase_independent",
    "primary_verification_agreement", "spectrum_status", "notes",
]

EXCLUSION_FIELDS = [
    "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "spectrum_status", "exclusion_reason", "notes",
]


@dataclass(frozen=True)
class Args:
    manifest: Path
    legacy_pilot_dir: Path
    output_dir: Path
    corrected_pilot_dir: Path
    k_max: int
    n_spectrum_roots: int
    n_candidate_roots: int
    verification_candidate_roots: int
    lambda_min: float
    lambda_max: float | None
    scan_step: float
    local_points: int
    adaptive_depth: int
    sigma_prefilter: float
    sigma_accept: float
    sigma_ratio_accept: float
    root_dedup_tol: float
    root_match_tol: float
    reuse_cache: bool
    force_recompute: bool
    plot_only: bool
    smoke: bool
    run_stress_cases: bool
    skip_corrected_pilot_analysis: bool


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        value_f = float(value)
        if math.isnan(value_f):
            return "nan"
        if math.isinf(value_f):
            return "inf" if value_f > 0.0 else "-inf"
        return f"{value_f:.16e}"
    if isinstance(value, (bool, np.bool_)):
        return "true" if bool(value) else "false"
    if isinstance(value, (list, tuple, dict)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def bool_value(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def finite_float(row: Mapping[str, object], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))
    except (TypeError, ValueError):
        return float("nan")


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Audit first-12 completeness for general coupled EB/Timoshenko 6x6 spectra.",
    )
    parser.add_argument("--manifest", type=Path, default=pilot.DEFAULT_MANIFEST)
    parser.add_argument("--legacy-pilot-dir", type=Path, default=DEFAULT_LEGACY_PILOT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--corrected-pilot-dir", type=Path, default=DEFAULT_CORRECTED_PILOT_DIR)
    parser.add_argument("--k-max", type=int, default=10)
    parser.add_argument("--n-spectrum-roots", type=int, default=12)
    parser.add_argument("--n-candidate-roots", type=int, default=20)
    parser.add_argument("--verification-candidate-roots", type=int, default=24)
    parser.add_argument("--lambda-min", type=float, default=complete.DEFAULT_SCAN_START)
    parser.add_argument("--lambda-max", type=float, default=None)
    parser.add_argument("--scan-step", type=float, default=complete.DEFAULT_SCAN_STEP)
    parser.add_argument("--local-points", type=int, default=complete.DEFAULT_LOCAL_POINTS)
    parser.add_argument("--adaptive-depth", type=int, default=complete.DEFAULT_ADAPTIVE_DEPTH)
    parser.add_argument("--sigma-prefilter", type=float, default=complete.DEFAULT_SIGMA_PREFILTER)
    parser.add_argument("--sigma-accept", type=float, default=complete.DEFAULT_SIGMA_ACCEPT)
    parser.add_argument("--sigma-ratio-accept", type=float, default=complete.DEFAULT_SIGMA_RATIO_ACCEPT)
    parser.add_argument("--root-dedup-tol", type=float, default=complete.DEFAULT_ROOT_DEDUP_TOL)
    parser.add_argument("--root-match-tol", type=float, default=complete.DEFAULT_ROOT_MATCH_TOL)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--run-stress-cases", action="store_true")
    parser.add_argument("--skip-corrected-pilot-analysis", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    output_dir = repo_path(Path(ns.output_dir))
    corrected_dir = repo_path(Path(ns.corrected_pilot_dir))
    if ns.smoke:
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            output_dir = repo_path(SMOKE_OUTPUT_DIR)
        if Path(ns.corrected_pilot_dir) == DEFAULT_CORRECTED_PILOT_DIR:
            corrected_dir = repo_path(SMOKE_CORRECTED_PILOT_DIR)
    args = Args(
        manifest=repo_path(Path(ns.manifest)),
        legacy_pilot_dir=repo_path(Path(ns.legacy_pilot_dir)),
        output_dir=output_dir,
        corrected_pilot_dir=corrected_dir,
        k_max=int(ns.k_max),
        n_spectrum_roots=int(ns.n_spectrum_roots),
        n_candidate_roots=int(ns.n_candidate_roots),
        verification_candidate_roots=int(ns.verification_candidate_roots),
        lambda_min=float(ns.lambda_min),
        lambda_max=float(ns.lambda_max) if ns.lambda_max is not None else None,
        scan_step=float(ns.scan_step),
        local_points=int(ns.local_points),
        adaptive_depth=int(ns.adaptive_depth),
        sigma_prefilter=float(ns.sigma_prefilter),
        sigma_accept=float(ns.sigma_accept),
        sigma_ratio_accept=float(ns.sigma_ratio_accept),
        root_dedup_tol=float(ns.root_dedup_tol),
        root_match_tol=float(ns.root_match_tol),
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        plot_only=bool(ns.plot_only),
        smoke=bool(ns.smoke),
        run_stress_cases=bool(ns.run_stress_cases),
        skip_corrected_pilot_analysis=bool(ns.skip_corrected_pilot_analysis),
    )
    settings_from_args(args).validate(k_max=args.k_max)
    if args.output_dir.resolve() == args.legacy_pilot_dir.resolve():
        raise ValueError("general audit output directory must differ from legacy pilot directory")
    if args.corrected_pilot_dir.resolve() == args.legacy_pilot_dir.resolve():
        raise ValueError("corrected and legacy pilot directories must differ")
    return args


def settings_from_args(args: Args) -> complete.SearchSettings:
    return complete.SearchSettings(
        requested_roots=args.n_spectrum_roots,
        candidate_roots=args.n_candidate_roots,
        verification_candidate_roots=args.verification_candidate_roots,
        lambda_min=args.lambda_min,
        lambda_max=args.lambda_max,
        scan_step=args.scan_step,
        local_points=args.local_points,
        adaptive_depth=args.adaptive_depth,
        sigma_prefilter=args.sigma_prefilter,
        sigma_accept=args.sigma_accept,
        sigma_ratio_accept=args.sigma_ratio_accept,
        root_dedup_tol=args.root_dedup_tol,
        root_match_tol=args.root_match_tol,
    )


def deterministic_cases(cases: Sequence[pilot.PilotCase]) -> list[pilot.PilotCase]:
    return sorted(
        cases,
        key=lambda case: (
            0 if abs(case.beta_deg) <= 1.0e-14 and abs(case.eta) <= 1.0e-14 else 1,
            abs(case.beta_deg), abs(case.mu), abs(case.eta), case.epsilon_0, case.case_id,
        ),
    )


def stress_cases(smoke: bool) -> list[pilot.PilotCase]:
    rows: list[pilot.PilotCase] = []
    for r_case, (epsilon, _left, _right) in R_CASES.items():
        for offset in STRESS_EPSILON_OFFSETS:
            for beta in STRESS_BETA_VALUES:
                case_id = f"{r_case}_e{offset:+.0e}_b{beta:g}".replace(".", "p").replace("+", "p").replace("-", "m")
                rows.append(
                    pilot.PilotCase(
                        case_id=case_id,
                        epsilon_0=epsilon + offset,
                        beta_deg=beta,
                        mu=0.0,
                        eta=0.0,
                        case_groups=("close_pair_stress", r_case),
                        notes=f"{r_case}; epsilon offset {offset:+.1e}; diagnostic stress only",
                    )
                )
    if smoke:
        return [row for row in rows if row.case_id.startswith("R1_") and row.beta_deg in {0.0, 0.1} and abs(row.epsilon_0 - R_CASES["R1"][0]) <= 1.0e-15]
    return rows


def geometry(case: pilot.PilotCase) -> complete.Geometry:
    return complete.Geometry(case.epsilon_0, case.beta_deg, case.mu, case.eta)


def geometry_distance(left: complete.Geometry, right: complete.Geometry) -> float:
    return (
        abs(left.epsilon_0 - right.epsilon_0) / 0.05
        + abs(left.beta_deg - right.beta_deg) / 90.0
        + abs(left.mu - right.mu)
        + abs(left.eta - right.eta)
    )


def nearest_roots(
    target: complete.Geometry,
    bank: Sequence[tuple[complete.Geometry, tuple[float, ...]]],
) -> tuple[float, ...]:
    if not bank:
        return ()
    return min(bank, key=lambda item: (geometry_distance(target, item[0]), item[0].beta_deg, item[0].mu, item[0].eta))[1]


def legacy_result(model: str, case: pilot.PilotCase, count: int) -> legacy_roots.RootResult:
    spec = legacy_roots.CaseSpec(mu=case.mu, eta=case.eta, epsilon=case.epsilon_0)
    return legacy_roots.solve_model(spec, case.beta_deg, count, model)


def nearest_match(value: float, candidates: Sequence[float], tolerance: float) -> tuple[bool, float]:
    finite = [float(item) for item in candidates if math.isfinite(float(item))]
    if not finite:
        return False, float("nan")
    difference = min(abs(float(value) - item) for item in finite)
    return difference <= float(tolerance), difference


def case_rows(
    case: pilot.PilotCase,
    model: str,
    result: complete.CompleteSpectrumResult,
    raw: legacy_roots.RootResult,
    settings: complete.SearchSettings,
) -> dict[str, list[dict[str, object]]]:
    base = {
        "case_id": case.case_id,
        "epsilon_0": case.epsilon_0,
        "beta_deg": case.beta_deg,
        "mu": case.mu,
        "eta": case.eta,
        "model": model,
    }
    raw_values = [float(value) for value in raw.roots if math.isfinite(float(value))]
    oracle = (
        complete.straight_oracle_values(model, result.geometry, settings.requested_roots)
        if abs(case.beta_deg) <= 1.0e-14 and abs(case.eta) <= 1.0e-14
        else ()
    )
    verification_by_index = {root.sorted_index: root for root in result.verification.roots}
    root_rows: list[dict[str, object]] = []
    for root in result.primary.roots:
        raw_match, _raw_diff = nearest_match(root.Lambda, raw_values, settings.root_match_tol)
        oracle_match, _oracle_diff = nearest_match(root.Lambda, oracle, settings.root_match_tol)
        verify = verification_by_index.get(root.sorted_index)
        root_rows.append(
            {
                **base,
                "configuration": "primary",
                **asdict(root),
                "detection_sources": "+".join(root.detection_sources),
                "primary_status": root.acceptance_status,
                "verification_status": "matched" if verify is not None and abs(root.Lambda - verify.Lambda) <= settings.root_match_tol else "unmatched",
                "raw_legacy_match": raw_match,
                "factorized_oracle_match": oracle_match if oracle else "not_applicable",
            }
        )
    interval_rows = [
        {**base, "configuration": config.configuration, **row}
        for config in (result.primary, result.verification)
        for row in config.interval_rows
    ]
    candidate_rows = [
        {
            **base,
            "configuration": config.configuration,
            "Lambda": candidate.Lambda,
            "detection_sources": "+".join(candidate.detection_sources),
            "interval_left": candidate.interval_left,
            "interval_right": candidate.interval_right,
            "determinant": candidate.diagnostics.determinant,
            "sigma_1": candidate.diagnostics.sigma_1,
            "sigma_2": candidate.diagnostics.sigma_2,
            "sigma_ratio": candidate.diagnostics.sigma_ratio,
            "null_vector_pivot": candidate.diagnostics.null_vector_pivot,
            "interior_minimum": candidate.interior_minimum,
            "acceptance_status": candidate.acceptance_status,
            "notes": candidate.notes,
        }
        for config in (result.primary, result.verification)
        for candidate in config.candidates
    ]
    recovered = [root for root in result.roots if not nearest_match(root.Lambda, raw_values, settings.root_match_tol)[0]]
    summary = {
        **base,
        "requested_roots": settings.requested_roots,
        "found_roots": len(result.primary.roots),
        "first12_resolved": result.spectrum_status == "resolved_complete",
        "root11_available": result.root11_available,
        "root12_available": result.root12_available,
        "raw_root_count": raw.root_count_found,
        "recovered_root_count": len(recovered),
        "close_cluster_count": len({root.root_cluster_id for root in result.roots if root.root_cluster_id}),
        "multiplicity_count": sum(root.detected_nullity > 1 for root in result.roots),
        "independent_agreement": result.independent_agreement,
        "spectrum_status": result.spectrum_status,
        "exclusion_reason": result.exclusion_reason,
        "cache_status": result.cache_status,
    }
    primary_rows = [{**base, **row} for row in result.primary_vs_verification]
    oracle_rows: list[dict[str, object]] = []
    for index, oracle_value in enumerate(oracle, start=1):
        general_value = result.roots[index - 1].Lambda if index <= len(result.roots) else float("nan")
        raw_present, _raw_diff = nearest_match(oracle_value, raw_values, settings.root_match_tol)
        difference = abs(general_value - oracle_value) if math.isfinite(general_value) else float("nan")
        oracle_rows.append(
            {
                **base,
                "sorted_index": index,
                "Lambda_general": general_value,
                "Lambda_oracle": oracle_value,
                "absolute_difference": difference,
                "raw_legacy_present": raw_present,
                "general_recovered": not raw_present and difference <= settings.root_match_tol,
                "status": "pass" if difference <= settings.root_match_tol else "mismatch",
            }
        )
    operation_rows = [
        {
            **base,
            "cost_scope": "offline_general_completeness",
            "algorithm_version": result.algorithm_version,
            **asdict(result.operations),
            "wall_clock_seconds": "auxiliary_not_stored_per_cached_result",
        }
    ]
    exclusion_rows = [] if result.spectrum_status == "resolved_complete" else [
        {
            **base,
            "spectrum_status": result.spectrum_status,
            "exclusion_reason": result.exclusion_reason,
            "notes": "unresolved corrected spectrum blocks N_true and step 3",
        }
    ]
    return {
        "root": root_rows,
        "interval": interval_rows,
        "candidate": candidate_rows,
        "summary": [summary],
        "primary": primary_rows,
        "oracle": oracle_rows,
        "operation": operation_rows,
        "exclusion": exclusion_rows,
    }


def stress_row(case: pilot.PilotCase, model: str, result: complete.CompleteSpectrumResult) -> dict[str, object]:
    r_case = next((name for name in R_CASES if name in case.case_groups), "")
    epsilon_center, left_ref, right_ref = R_CASES[r_case]
    center = 0.5 * (left_ref + right_ref)
    window = [root.Lambda for root in result.roots if abs(root.Lambda - center) <= 0.18]
    if len(window) >= 2:
        pair = min(zip(window[:-1], window[1:]), key=lambda values: abs(values[1] - values[0]))
        first, second = pair
        gap = abs(second - first)
    else:
        first = window[0] if window else float("nan")
        second = gap = float("nan")
    pair_required = model == complete.MODEL_TIMO
    return {
        "case_id": case.case_id,
        "r_case": r_case,
        "epsilon_offset": case.epsilon_0 - epsilon_center,
        "epsilon_0": case.epsilon_0,
        "beta_deg": case.beta_deg,
        "model": model,
        "root_count": len(result.roots),
        "pair_root_1": first,
        "pair_root_2": second,
        "minimum_gap": gap,
        "both_roots_present": (len(window) >= 2) if pair_required else result.spectrum_status == "resolved_complete",
        "grid_phase_independent": result.independent_agreement,
        "primary_verification_agreement": result.independent_agreement,
        "spectrum_status": result.spectrum_status,
        "notes": (
            "Timoshenko close-pair perturbation tracks; roots do not replace the central geometry roots"
            if pair_required
            else "EB first-12 completeness control; the R-window Timoshenko close pair is not imposed on EB"
        ),
    }


def legacy_close_case_ids(legacy_pilot_dir: Path) -> set[str]:
    selected: set[str] = set()
    for row in read_csv(legacy_pilot_dir / "epsilon_pilot_matching_audit.csv"):
        clustered = bool(str(row.get("cluster_EB", "")).strip() or str(row.get("cluster_Timo", "")).strip())
        gaps = [finite_float(row, "gap_EB"), finite_float(row, "gap_Timo")]
        small_gap = any(math.isfinite(value) and value <= complete.DEFAULT_CLOSE_GAP for value in gaps)
        if clustered or small_gap:
            selected.add(str(row.get("case_id", "")))
    return selected


def pilot_close_stress_row(
    case: pilot.PilotCase,
    model: str,
    result: complete.CompleteSpectrumResult,
) -> dict[str, object]:
    values = [root.Lambda for root in result.roots]
    if len(values) >= 2:
        first, second = min(zip(values[:-1], values[1:]), key=lambda pair: pair[1] - pair[0])
        gap = second - first
    else:
        first = second = gap = float("nan")
    return {
        "case_id": case.case_id,
        "r_case": "legacy_pilot_close_cluster",
        "epsilon_offset": "not_applicable",
        "epsilon_0": case.epsilon_0,
        "beta_deg": case.beta_deg,
        "model": model,
        "root_count": len(result.roots),
        "pair_root_1": first,
        "pair_root_2": second,
        "minimum_gap": gap,
        "both_roots_present": len(values) >= 12,
        "grid_phase_independent": result.independent_agreement,
        "primary_verification_agreement": result.independent_agreement,
        "spectrum_status": result.spectrum_status,
        "notes": "legacy pilot close-cluster/small-gap control selected independently of N_true",
    }


def write_legacy_comparisons(legacy_dir: Path, corrected_dir: Path) -> tuple[Path, Path, Path]:
    legacy_geometry = {row["case_id"]: row for row in read_csv(legacy_dir / "epsilon_pilot_geometry_metrics.csv")}
    corrected_geometry = {row["case_id"]: row for row in read_csv(corrected_dir / "epsilon_pilot_geometry_metrics.csv")}
    legacy_modes = defaultdict(list)
    corrected_modes = defaultdict(list)
    for row in read_csv(legacy_dir / "epsilon_pilot_mode_metrics.csv"):
        legacy_modes[row["case_id"]].append(row)
    for row in read_csv(corrected_dir / "epsilon_pilot_mode_metrics.csv"):
        corrected_modes[row["case_id"]].append(row)
    legacy_predictions = {row["case_id"]: row for row in read_csv(legacy_dir / "analysis" / "epsilon_vs_existing_rules_geometry_comparison.csv")}
    corrected_predictions = {row["case_id"]: row for row in read_csv(corrected_dir / "analysis" / "epsilon_vs_existing_rules_geometry_comparison.csv")}
    comparison_rows: list[dict[str, object]] = []
    root_rows: list[dict[str, object]] = []
    for case_id in sorted(corrected_geometry):
        legacy = legacy_geometry.get(case_id, {})
        corrected = corrected_geometry.get(case_id, {})
        old_n = finite_float(legacy, "N_true")
        new_n = finite_float(corrected, "N_true")
        old_modes = sorted(legacy_modes.get(case_id, []), key=lambda row: int(float(row["sorted_index"])))
        new_modes = sorted(corrected_modes.get(case_id, []), key=lambda row: int(float(row["sorted_index"])))
        changed_eb = sum(
            abs(finite_float(left, "Lambda_EB") - finite_float(right, "Lambda_EB")) > complete.DEFAULT_ROOT_MATCH_TOL
            for left, right in zip(old_modes, new_modes)
        ) + abs(len(old_modes) - len(new_modes))
        changed_timo = sum(
            abs(finite_float(left, "Lambda_Timo") - finite_float(right, "Lambda_Timo")) > complete.DEFAULT_ROOT_MATCH_TOL
            for left, right in zip(old_modes, new_modes)
        ) + abs(len(old_modes) - len(new_modes))
        old_pred = legacy_predictions.get(case_id, {})
        new_pred = corrected_predictions.get(case_id, {})
        prediction_fields = sorted(field for field in set(old_pred) | set(new_pred) if field.startswith("N_hat_"))
        changed_predictions = [field for field in prediction_fields if old_pred.get(field, "") != new_pred.get(field, "")]
        comparison_rows.append(
            {
                "case_id": case_id,
                "legacy_N_true": old_n,
                "corrected_N_true": new_n,
                "N_true_difference": new_n - old_n if math.isfinite(old_n) and math.isfinite(new_n) else "",
                "legacy_first_failed_mode": legacy.get("first_true_failed_mode", ""),
                "corrected_first_failed_mode": corrected.get("first_true_failed_mode", ""),
                "changed_EB_root_count": changed_eb,
                "changed_Timo_root_count": changed_timo,
                "missing_legacy_roots": max(0, changed_eb + changed_timo),
                "recovered_roots": max(0, changed_eb + changed_timo),
                "changed_close_cluster": legacy.get("close_cluster_count", "") != corrected.get("close_cluster_count", ""),
                "changed_rule_predictions": ";".join(changed_predictions),
            }
        )
        for index, (left, right) in enumerate(zip(old_modes, new_modes), start=1):
            root_rows.append(
                {
                    "case_id": case_id,
                    "sorted_index": index,
                    "legacy_Lambda_EB": left.get("Lambda_EB", ""),
                    "corrected_Lambda_EB": right.get("Lambda_EB", ""),
                    "legacy_Lambda_Timo": left.get("Lambda_Timo", ""),
                    "corrected_Lambda_Timo": right.get("Lambda_Timo", ""),
                    "EB_changed": abs(finite_float(left, "Lambda_EB") - finite_float(right, "Lambda_EB")) > complete.DEFAULT_ROOT_MATCH_TOL,
                    "Timo_changed": abs(finite_float(left, "Lambda_Timo") - finite_float(right, "Lambda_Timo")) > complete.DEFAULT_ROOT_MATCH_TOL,
                }
            )
    comparison_path = write_csv(
        corrected_dir / "pilot_legacy_vs_complete_spectrum.csv",
        comparison_rows,
        [
            "case_id", "legacy_N_true", "corrected_N_true", "N_true_difference", "legacy_first_failed_mode",
            "corrected_first_failed_mode", "changed_EB_root_count", "changed_Timo_root_count", "missing_legacy_roots",
            "recovered_roots", "changed_close_cluster", "changed_rule_predictions",
        ],
    )
    root_path = write_csv(
        corrected_dir / "pilot_root_changes_by_case.csv",
        root_rows,
        [
            "case_id", "sorted_index", "legacy_Lambda_EB", "corrected_Lambda_EB", "legacy_Lambda_Timo",
            "corrected_Lambda_Timo", "EB_changed", "Timo_changed",
        ],
    )
    legacy_metrics = {
        (row.get("split_kind", ""), row.get("rule", "")): row
        for row in read_csv(legacy_dir / "analysis" / "epsilon_rule_validation_summary.csv")
    }
    corrected_metrics = {
        (row.get("split_kind", ""), row.get("rule", "")): row
        for row in read_csv(corrected_dir / "analysis" / "epsilon_rule_validation_summary.csv")
    }
    metric_rows: list[dict[str, object]] = []
    for key in sorted(set(legacy_metrics) | set(corrected_metrics)):
        old = legacy_metrics.get(key, {})
        new = corrected_metrics.get(key, {})
        metric_rows.append(
            {
                "split_kind": key[0],
                "rule": key[1],
                "legacy_false_safe_geometry_count": old.get("false_safe_geometry_count", ""),
                "corrected_false_safe_geometry_count": new.get("false_safe_geometry_count", ""),
                "legacy_retention": old.get("usable_frequency_retention", ""),
                "corrected_retention": new.get("usable_frequency_retention", ""),
                "legacy_mean_conservative_loss": old.get("mean_conservative_loss", ""),
                "corrected_mean_conservative_loss": new.get("mean_conservative_loss", ""),
                "legacy_unique_false_safe_geometries": old.get("false_safe_geometry_count", ""),
                "corrected_unique_false_safe_geometries": new.get("false_safe_geometry_count", ""),
            }
        )
    metric_path = write_csv(
        corrected_dir / "pilot_rule_metrics_legacy_vs_complete_spectrum.csv",
        metric_rows,
        [
            "split_kind", "rule", "legacy_false_safe_geometry_count", "corrected_false_safe_geometry_count",
            "legacy_retention", "corrected_retention", "legacy_mean_conservative_loss",
            "corrected_mean_conservative_loss", "legacy_unique_false_safe_geometries",
            "corrected_unique_false_safe_geometries",
        ],
    )
    return comparison_path, metric_path, root_path


def report_from_rows(
    args: Args,
    summaries: Sequence[Mapping[str, object]],
    oracle_rows: Sequence[Mapping[str, object]],
    stress_rows: Sequence[Mapping[str, object]],
    operation_rows: Sequence[Mapping[str, object]],
    *,
    pilot_comparison_rows: Sequence[Mapping[str, object]] = (),
    rule_metric_rows: Sequence[Mapping[str, object]] = (),
) -> Path:
    pilot_summaries = [row for row in summaries if not str(row.get("case_id", "")).startswith("R")]
    all_pilot_cases = {str(row.get("case_id")) for row in pilot_summaries}
    resolved_pilot_cases = {
        case_id
        for case_id in all_pilot_cases
        if {
            str(row.get("model"))
            for row in pilot_summaries
            if str(row.get("case_id")) == case_id
            and str(row.get("spectrum_status")) == "resolved_complete"
        } == {complete.MODEL_EB, complete.MODEL_TIMO}
    }
    unresolved_pilot = [
        row for row in pilot_summaries if str(row.get("spectrum_status")) != "resolved_complete"
    ]
    stress_summaries = [row for row in summaries if str(row.get("case_id", "")).startswith("R")]
    unresolved_stress = [
        row for row in stress_summaries if str(row.get("spectrum_status")) != "resolved_complete"
    ]
    corrected_exclusions = read_csv(args.corrected_pilot_dir / "epsilon_pilot_exclusion_audit.csv")
    corrected_excluded_cases = {str(row.get("case_id")) for row in corrected_exclusions}
    corrected_pilot_total = len({str(row.get("case_id")) for row in pilot_comparison_rows})
    corrected_pilot_included = corrected_pilot_total - len(corrected_excluded_cases)
    oracle_failures = [row for row in oracle_rows if str(row.get("status")) != "pass"]
    stress_failures = [
        row for row in stress_rows
        if not bool_value(row.get("both_roots_present"))
        or str(row.get("spectrum_status")) != "resolved_complete"
        or not bool_value(row.get("primary_verification_agreement"))
    ]
    ready = (
        len(all_pilot_cases) == (3 if args.smoke else 21)
        and resolved_pilot_cases == all_pilot_cases
        and not oracle_failures
        and (not args.run_stress_cases or not stress_failures)
    )
    decision = "ready_for_targeted_step3" if ready and not args.smoke else "not_ready_for_step3"
    raw_missing = sum(
        int(finite_float(row, "recovered_root_count"))
        for row in summaries
        if math.isfinite(finite_float(row, "recovered_root_count"))
    )
    recovered = raw_missing
    minimum_gap = min(
        (finite_float(row, "minimum_gap") for row in stress_rows if math.isfinite(finite_float(row, "minimum_gap"))),
        default=float("nan"),
    )
    changed_n = [
        row for row in pilot_comparison_rows
        if math.isfinite(finite_float(row, "N_true_difference")) and abs(finite_float(row, "N_true_difference")) > 0.0
    ]
    changed_first_failed = [
        row for row in pilot_comparison_rows
        if str(row.get("legacy_first_failed_mode", "")) != str(row.get("corrected_first_failed_mode", ""))
    ]
    changed_root_cases = [
        row for row in pilot_comparison_rows
        if int(finite_float(row, "changed_EB_root_count")) > 0
        or int(finite_float(row, "changed_Timo_root_count")) > 0
    ]
    changed_close_clusters = [
        row for row in pilot_comparison_rows if bool_value(row.get("changed_close_cluster"))
    ]
    changed_rule_prediction_cases = [
        row for row in pilot_comparison_rows if str(row.get("changed_rule_predictions", "")).strip()
    ]
    false_safe_metric_changes = [
        row for row in rule_metric_rows
        if finite_float(row, "legacy_false_safe_geometry_count")
        != finite_float(row, "corrected_false_safe_geometry_count")
        or finite_float(row, "legacy_unique_false_safe_geometries")
        != finite_float(row, "corrected_unique_false_safe_geometries")
    ]
    retention_changes = [
        abs(finite_float(row, "legacy_retention") - finite_float(row, "corrected_retention"))
        for row in rule_metric_rows
        if math.isfinite(finite_float(row, "legacy_retention"))
        and math.isfinite(finite_float(row, "corrected_retention"))
    ]
    conservative_loss_changes = [
        abs(
            finite_float(row, "legacy_mean_conservative_loss")
            - finite_float(row, "corrected_mean_conservative_loss")
        )
        for row in rule_metric_rows
        if math.isfinite(finite_float(row, "legacy_mean_conservative_loss"))
        and math.isfinite(finite_float(row, "corrected_mean_conservative_loss"))
    ]
    all_held_out_metrics = [
        row for row in rule_metric_rows if str(row.get("split_kind")) == "all_held_out"
    ]

    def retention_ranking(field: str) -> str:
        ranked = sorted(
            (
                (str(row.get("rule")), finite_float(row, field))
                for row in all_held_out_metrics
                if math.isfinite(finite_float(row, field))
            ),
            key=lambda item: (-item[1], item[0]),
        )
        return " > ".join(name for name, _value in ranked) or "unavailable"

    all_held_out_false_safe = ", ".join(
        f"{row.get('rule')} {int(finite_float(row, 'legacy_unique_false_safe_geometries'))}"
        f"->{int(finite_float(row, 'corrected_unique_false_safe_geometries'))}"
        for row in all_held_out_metrics
        if math.isfinite(finite_float(row, "legacy_unique_false_safe_geometries"))
        and math.isfinite(finite_float(row, "corrected_unique_false_safe_geometries"))
    ) or "unavailable"
    close_cluster_total = sum(
        int(finite_float(row, "close_cluster_count"))
        for row in summaries
        if math.isfinite(finite_float(row, "close_cluster_count"))
    )
    multiplicity_total = sum(
        int(finite_float(row, "multiplicity_count"))
        for row in summaries
        if math.isfinite(finite_float(row, "multiplicity_count"))
    )
    unresolved_pilot_description = "; ".join(
        f"{row.get('case_id')} {row.get('model')} ({row.get('exclusion_reason')})"
        for row in unresolved_pilot
    ) or "none"
    primitive_names = (
        "characteristic_matrix_evaluations",
        "full_6x6_SVD_calls",
        "adaptive_interval_subdivisions",
        "Brent_or_bisection_evaluations",
        "cache_hits",
        "cache_misses",
    )
    primitive_totals = {
        name: sum(
            int(finite_float(row, name))
            for row in operation_rows
            if math.isfinite(finite_float(row, name))
        )
        for name in primitive_names
    }
    lines = [
        "# EB/Timoshenko General Spectrum Completeness Report",
        "",
        "## Scope",
        "",
        "- General coupled Euler--Bernoulli and Timoshenko 6x6 matrices over selected beta/mu/eta geometries.",
        f"- Target `K={args.k_max}`; first `{args.n_spectrum_roots}` roots audited, with candidate margins `{args.n_candidate_roots}/{args.verification_candidate_roots}`.",
        "- This is a finite numerical audit, not a mathematical root-count proof.",
        "- No physical model, formula, matrix, production solver, FEM, or article workflow changed.",
        "- Research step 3 was not implemented or run.",
        "",
        "## Why Sign Scan Is Insufficient",
        "",
        "Two simple roots can lie inside one determinant interval, produce two sign changes, and leave equal endpoint signs. The corrected straight baseline demonstrated this defect. Determinant brackets remain useful candidates, but every accepted root here is checked by the row-normalized full-matrix SVD.",
        "",
        "## Recovery Algorithm",
        "",
        "Candidate paths are legacy and shifted sign scans, an independent half-step shifted scan, global sigma/sigma-ratio valleys, adaptive local intervals, EB seed windows for Timoshenko, deterministic continuation seeds, and straight-oracle seeds used only for full-matrix verification. Close roots are retained unless Lambda, self-MAC, and detection history all support deduplication. Exact nullity and coalesced-track metadata remain explicit.",
        "",
        "## Straight Oracle Validation",
        "",
        f"- Oracle rows: `{len(oracle_rows)}`; failures: `{len(oracle_failures)}`.",
        f"- Raw sign-scan missing roots within audited prefixes: `{raw_missing}`; recovered: `{recovered}`.",
        "- R1--R3 straight roots are accepted only after the unchanged full 6x6 matrix passes the SVD criterion.",
        "",
        "## Small-Angle Stress Cases",
        "",
        f"- Stress rows: `{len(stress_rows)}`; failures/unresolved: `{len(stress_failures)}`.",
        f"- Minimum recovered pair gap: `{minimum_gap:.12g}`." if math.isfinite(minimum_gap) else "- Minimum recovered pair gap: unavailable.",
        "",
        "## 21-Case Pilot",
        "",
        f"- General 6x6 audit resolved both models for `{len(resolved_pilot_cases)}/{len(all_pilot_cases)}` pilot cases; unresolved pilot model rows: `{len(unresolved_pilot)}`.",
        (
        f"- Corrected auto-spectrum pilot included `{corrected_pilot_included}/{corrected_pilot_total}` cases; "
            f"explicitly excluded: `{', '.join(sorted(corrected_excluded_cases)) or 'none'}`."
            if corrected_pilot_total
            else "- Corrected auto-spectrum pilot comparison: unavailable."
        ),
        f"- Unresolved requested stress model rows: `{len(unresolved_stress)}` (reported separately from pilot rows).",
        f"- Cases with changed `N_true`: `{len(changed_n)}` ({', '.join(str(row.get('case_id')) for row in changed_n) or 'none'}).",
        f"- Cases with changed first-failed mode: `{len(changed_first_failed)}`; with changed first-ten EB/Timoshenko roots: `{len(changed_root_cases)}`.",
        f"- Recorded close clusters across audit summaries: `{close_cluster_total}`; asserted multiplicities: `{multiplicity_total}`; corrected-pilot cases with changed close-cluster status: `{len(changed_close_clusters)}`.",
        f"- Cases whose rule predictions changed after explicit exclusion/reanalysis: `{len(changed_rule_prediction_cases)}`.",
        f"- General unresolved quality rows: {unresolved_pilot_description}.",
        "- Root 11 is the right gap guard for mode 10; root 12 is candidate-boundary margin only.",
        "- An unresolved corrected spectrum is excluded and is never replaced by legacy roots.",
        "",
        "## Rule Comparison",
        "",
        f"- Legacy/corrected rule metric rows: `{len(rule_metric_rows)}`.",
        "- Compared E0-ref, Emax-ref, E0-cal, Emax-cal, Rule A, Rule A-gap, Rule B, Rule C, and Rule D.",
        f"- False-safe and unique-false-safe geometry counts changed in `{len(false_safe_metric_changes)}/{len(rule_metric_rows)}` rows.",
        f"- All-held-out unique false-safe geometries, legacy->corrected: {all_held_out_false_safe}.",
        f"- All-held-out retention ranking, legacy: {retention_ranking('legacy_retention')}.",
        f"- All-held-out retention ranking, corrected: {retention_ranking('corrected_retention')}.",
        f"- Maximum absolute retention change: `{max(retention_changes, default=0.0):.12g}`; maximum absolute mean conservative-loss change: `{max(conservative_loss_changes, default=0.0):.12g}`.",
        "- Rules and thresholds were not tuned in response to corrected roots; the existing CSV-only postprocessor was reused.",
        "",
        "## Operation Counts",
        "",
        f"- Per-case/model operation rows: `{len(operation_rows)}`. Primitive counts are in `general_spectrum_operation_counts.csv`; wall time is auxiliary only.",
        (
            "- Totals: "
            f"matrix evaluations `{primitive_totals['characteristic_matrix_evaluations']}`, "
            f"full 6x6 SVD calls `{primitive_totals['full_6x6_SVD_calls']}`, "
            f"adaptive subdivisions `{primitive_totals['adaptive_interval_subdivisions']}`, "
            f"Brent/bisection evaluations `{primitive_totals['Brent_or_bisection_evaluations']}`, "
            f"cache hits/misses `{primitive_totals['cache_hits']}/{primitive_totals['cache_misses']}`."
        ),
        "",
        "## Limitations",
        "",
        "- Audit-resolved under independent search configurations does not mean mathematically certified.",
        "- Only the selected pilot and requested stress geometries were checked.",
        "- The general production root solver is unchanged; Timoshenko remains a one-dimensional reference.",
        "- No lower-envelope conclusion is available yet; unresolved spectra must block step 3.",
        "",
        "## Step-3 Readiness Decision",
        "",
        f"`{decision}`",
        "",
    ]
    if decision == "not_ready_for_step3":
        reasons = []
        if args.smoke:
            reasons.append("smoke run is not a full readiness audit")
        if resolved_pilot_cases != all_pilot_cases:
            reasons.append(
                "not all pilot geometries are audit-resolved in both models "
                f"({', '.join(sorted(all_pilot_cases - resolved_pilot_cases))})"
            )
        if oracle_failures:
            reasons.append("straight oracle mismatch remains")
        if args.run_stress_cases and stress_failures:
            reasons.append("small-angle stress failures remain")
        lines.append("Reasons: " + "; ".join(reasons or ["full required stress set was not run"]) + ".")
    report = args.output_dir / "eb_timo_general_spectrum_completeness_report.md"
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def plot_only(args: Args) -> dict[str, object]:
    summary_rows = read_csv(args.output_dir / "general_spectrum_completeness_summary.csv")
    if not summary_rows:
        raise FileNotFoundError("plot-only requires existing general spectrum CSV outputs")
    report = report_from_rows(
        args,
        summary_rows,
        read_csv(args.output_dir / "general_spectrum_factorized_oracle_comparison.csv"),
        read_csv(args.output_dir / "general_spectrum_close_pair_stress_audit.csv"),
        read_csv(args.output_dir / "general_spectrum_operation_counts.csv"),
        pilot_comparison_rows=read_csv(args.corrected_pilot_dir / "pilot_legacy_vs_complete_spectrum.csv"),
        rule_metric_rows=read_csv(args.corrected_pilot_dir / "pilot_rule_metrics_legacy_vs_complete_spectrum.csv"),
    )
    print("plot-only: root calculations performed: 0")
    return {"report": report, "root_calculations": 0}


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    if args.plot_only:
        return plot_only(args)
    settings = settings_from_args(args)
    manifest_cases = pilot.load_manifest(args.manifest)
    pilot.validate_manifest(manifest_cases, require_full_manifest=not args.smoke)
    selected_pilot = pilot.select_cases(manifest_cases, smoke=args.smoke)
    audit_cases = deterministic_cases(selected_pilot)
    stress = stress_cases(args.smoke) if args.run_stress_cases else []
    cache = complete.GeneralSpectrumCache(
        args.corrected_pilot_dir / "cache",
        reuse_cache=args.reuse_cache,
        force_recompute=args.force_recompute,
    )
    banks: dict[str, list[tuple[complete.Geometry, tuple[float, ...]]]] = {
        complete.MODEL_EB: [], complete.MODEL_TIMO: []
    }
    collections: dict[str, list[dict[str, object]]] = defaultdict(list)
    results_by_case_model: dict[tuple[str, str], complete.CompleteSpectrumResult] = {}
    close_pilot_ids = legacy_close_case_ids(args.legacy_pilot_dir) if args.run_stress_cases else set()
    started = time.perf_counter()
    for case in [*audit_cases, *deterministic_cases(stress)]:
        current_geometry = geometry(case)
        eb_result = cache.resolve(
            complete.MODEL_EB,
            current_geometry,
            settings,
            continuation_seeds=nearest_roots(current_geometry, banks[complete.MODEL_EB]),
        )
        if eb_result.spectrum_status == "resolved_complete":
            banks[complete.MODEL_EB].append((current_geometry, tuple(root.Lambda for root in eb_result.primary.roots)))
        timo_result = cache.resolve(
            complete.MODEL_TIMO,
            current_geometry,
            settings,
            continuation_seeds=nearest_roots(current_geometry, banks[complete.MODEL_TIMO]),
            eb_seed_roots=eb_result.values,
        )
        if timo_result.spectrum_status == "resolved_complete":
            banks[complete.MODEL_TIMO].append((current_geometry, tuple(root.Lambda for root in timo_result.primary.roots)))
        for model, result in ((complete.MODEL_EB, eb_result), (complete.MODEL_TIMO, timo_result)):
            raw = legacy_result(model, case, args.n_candidate_roots)
            products = case_rows(case, model, result, raw, settings)
            for key, rows in products.items():
                collections[key].extend(rows)
            results_by_case_model[(case.case_id, model)] = result
            if case in stress:
                collections["stress"].append(stress_row(case, model, result))
            elif args.run_stress_cases and case.case_id in close_pilot_ids:
                collections["stress"].append(pilot_close_stress_row(case, model, result))
        print(
            f"completed {case.case_id}: EB={eb_result.spectrum_status}, "
            f"Timo={timo_result.spectrum_status}"
        )
    elapsed = time.perf_counter() - started
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / OUTPUT_NAMES[0], collections["root"], ROOT_FIELDS)
    write_csv(args.output_dir / OUTPUT_NAMES[1], collections["interval"], INTERVAL_FIELDS)
    write_csv(args.output_dir / OUTPUT_NAMES[2], collections["candidate"], CANDIDATE_FIELDS)
    write_csv(args.output_dir / OUTPUT_NAMES[3], collections["summary"], SUMMARY_FIELDS)
    write_csv(args.output_dir / OUTPUT_NAMES[4], collections["primary"], PRIMARY_FIELDS)
    write_csv(args.output_dir / OUTPUT_NAMES[5], collections["oracle"], ORACLE_FIELDS)
    write_csv(args.output_dir / OUTPUT_NAMES[6], collections["stress"], STRESS_FIELDS)
    operation_fields = [
        "case_id", "epsilon_0", "beta_deg", "mu", "eta", "model", "cost_scope", "algorithm_version",
        *complete.OperationCounts.__dataclass_fields__.keys(), "wall_clock_seconds",
    ]
    write_csv(args.output_dir / OUTPUT_NAMES[7], collections["operation"], operation_fields)
    write_csv(args.output_dir / OUTPUT_NAMES[8], collections["exclusion"], EXCLUSION_FIELDS)

    pilot_comparison_rows: list[dict[str, str]] = []
    rule_metric_rows: list[dict[str, str]] = []
    pilot_result: dict[str, object] | None = None
    analysis_result: dict[str, object] | None = None
    if not args.skip_corrected_pilot_analysis:
        pilot_argv = [
            "--manifest", str(args.manifest),
            "--output-dir", str(args.corrected_pilot_dir),
            "--cache-dir", str(args.corrected_pilot_dir / "cache"),
            "--k-max", str(args.k_max),
            "--n-spectrum-roots", str(args.n_spectrum_roots),
            "--n-candidate-roots", str(args.n_candidate_roots),
            "--verification-candidate-roots", str(args.verification_candidate_roots),
            "--spectrum-method", pilot.SPECTRUM_METHOD_AUTO,
            "--reuse-cache", "--force",
        ]
        if args.smoke:
            pilot_argv.append("--smoke")
        pilot_result = pilot.main(pilot_argv)
        analysis_argv = [
            "--pilot-dir", str(args.corrected_pilot_dir),
            "--output-dir", str(args.corrected_pilot_dir / "analysis"),
            "--k-max", str(args.k_max),
            "--force",
        ]
        if args.smoke:
            analysis_argv.append("--smoke")
        analysis_result = pilot_analysis.main(analysis_argv)
        write_legacy_comparisons(args.legacy_pilot_dir, args.corrected_pilot_dir)
        pilot_comparison_rows = read_csv(args.corrected_pilot_dir / "pilot_legacy_vs_complete_spectrum.csv")
        rule_metric_rows = read_csv(args.corrected_pilot_dir / "pilot_rule_metrics_legacy_vs_complete_spectrum.csv")
    report = report_from_rows(
        args,
        collections["summary"],
        collections["oracle"],
        collections["stress"],
        collections["operation"],
        pilot_comparison_rows=pilot_comparison_rows,
        rule_metric_rows=rule_metric_rows,
    )
    print(f"general audit wall-clock seconds (auxiliary): {elapsed:.6f}")
    print(f"output: {args.output_dir}")
    return {
        "args": args,
        "collections": collections,
        "report": report,
        "pilot_result": pilot_result,
        "analysis_result": analysis_result,
        "elapsed_seconds": elapsed,
    }


if __name__ == "__main__":
    main()
