from __future__ import annotations

import argparse
from collections import Counter
import csv
from dataclasses import dataclass
import json
import math
from math import isfinite
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
    audit_eb_validity_fixed_epsilon_geometry_scan as fixed_scan,
)
from scripts.analysis.thickness_mismatch.audits import (  # noqa: E402
    audit_eb_validity_vs_timoshenko_stage1 as stage1,
)
from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_safe_prefix_certification as certification,
)


DEFAULT_MANIFEST = SCRIPT_PATH.parent / "data" / "eb_epsilon_apriori_pilot_cases.csv"
DEFAULT_OUTPUT_DIR = Path("results") / "eb_epsilon_apriori_pilot"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "eb_epsilon_apriori_pilot"
DEFAULT_K_MAX = 10
DEFAULT_N_CANDIDATE_ROOTS = 20
DEFAULT_N_SHAPE_POINTS = fixed_scan.DEFAULT_SHAPE_POINTS
SMOKE_N_SHAPE_POINTS = 51
DEFAULT_CLUSTER_GAP_THRESHOLD = fixed_scan.DEFAULT_CLUSTER_GAP_THRESHOLD
FREQUENCY_ERROR_THRESHOLD = 0.10
SMOKE_CASE_IDS = ("B01", "G04", "M04")

SOURCE_OUTPUT_NAMES = (
    "epsilon_pilot_case_manifest_resolved.csv",
    "epsilon_pilot_mode_metrics.csv",
    "epsilon_pilot_geometry_metrics.csv",
    "epsilon_pilot_matching_audit.csv",
    "epsilon_pilot_exclusion_audit.csv",
    "epsilon_pilot_root_operation_counts.csv",
    "epsilon_pilot_generation_report.md",
)

RESOLVED_FIELDS = [
    "case_id",
    "geometry_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "tau_1",
    "tau_2",
    "epsilon_1",
    "epsilon_2",
    "epsilon_max",
    "radius_to_length_1",
    "radius_to_length_2",
    "diameter_to_length_1",
    "diameter_to_length_2",
    "case_groups",
    "notes",
]

MODE_FIELDS = [
    "case_id",
    "geometry_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "case_groups",
    "sorted_index",
    "Lambda_EB",
    "Lambda_Timo",
    "delta_f",
    "stored_delta_f",
    "delta_f_consistent",
    "true_pass_10",
    "tau_1",
    "tau_2",
    "epsilon_1",
    "epsilon_2",
    "epsilon_max",
    "chi_max_EB",
    "chi_eff_EB",
    "Theta_max_EB",
    "Pi_shear_EB",
    "Pi_rotary_EB",
    "Pi_EB",
    "EB_axial_fraction",
    "EB_bending_fraction",
    "EB_classification",
    "Timo_axial_fraction",
    "Timo_bending_fraction",
    "Timo_shear_fraction",
    "Timo_bending_shear_fraction",
    "Timo_classification",
    "gap_left_EB",
    "gap_right_EB",
    "gap_EB",
    "gap_sidedness",
    "Lambda_EB_guard_11",
    "mode10_upper_guard_source",
    "matched_timo_sorted_index",
    "shape_MAC_EB_Timo",
    "matching_status",
    "cluster_id",
    "cluster_size_EB",
    "cluster_size_Timo",
    "root_warning",
    "candidate_boundary",
    "included_in_certification",
    "quality_status",
    "source_warning",
    "notes",
]

GEOMETRY_FIELDS = [
    "case_id",
    "geometry_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "epsilon_1",
    "epsilon_2",
    "epsilon_max",
    "case_groups",
    "N_true",
    "first_true_failed_mode",
    "pass_pattern_10",
    "late_pass_after_failure",
    "complete_K10",
    "quality_status",
    "root_warning_count",
    "candidate_boundary_count",
    "ambiguous_matching_count",
    "close_cluster_count",
    "root11_guard_available",
    "exclusion_reason",
    "notes",
]

MATCHING_FIELDS = [
    "case_id",
    "geometry_id",
    "epsilon_0",
    "case_groups",
    *fixed_scan.MATCHING_FIELDS,
]

EXCLUSION_FIELDS = [
    "case_id",
    "geometry_id",
    "epsilon_0",
    "beta_deg",
    "mu",
    "eta",
    "complete_K10",
    "root_warning_count",
    "candidate_boundary_count",
    "exclusion_reason",
    "notes",
]

ROOT_OPERATION_FIELDS = [
    "case_id",
    "cost_scope",
    "measurement_status",
    "n_candidate_roots",
    "n_shape_points",
    "EB_root_solver_calls",
    "Timoshenko_root_solver_calls",
    "root_cache_hits",
    "root_cache_misses",
    "boundary_retries",
    "EB_root_seconds",
    "Timoshenko_root_seconds",
    "shape_and_predictor_seconds",
    "matching_seconds",
    "matching_pair_rows",
    "EB_characteristic_matrix_evaluations",
    "EB_determinant_evaluations",
    "Timoshenko_characteristic_matrix_evaluations",
    "Timoshenko_determinant_evaluations",
    "root_bracketing_intervals",
    "Brent_function_evaluations",
    "notes",
]


@dataclass(frozen=True)
class PilotCase:
    case_id: str
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float
    case_groups: tuple[str, ...]
    notes: str


@dataclass(frozen=True)
class Args:
    manifest: Path
    output_dir: Path
    cache_dir: Path
    k_max: int
    n_candidate_roots: int
    n_shape_points: int
    cluster_gap_threshold: float
    reuse_cache: bool
    force_recompute: bool
    force: bool
    smoke: bool


def repo_path(path: Path) -> Path:
    path_obj = Path(path)
    return path_obj if path_obj.is_absolute() else REPO_ROOT / path_obj


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
    if isinstance(value, (dict, list, tuple)):
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


def finite_float(row: Mapping[str, object], key: str) -> float:
    try:
        return float(row.get(key, float("nan")))
    except (TypeError, ValueError):
        return float("nan")


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Generate the manifest-driven K=10 EB epsilon a-priori pilot source CSV files.",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--k-max", type=int, default=DEFAULT_K_MAX)
    parser.add_argument("--n-candidate-roots", type=int, default=DEFAULT_N_CANDIDATE_ROOTS)
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_N_SHAPE_POINTS)
    parser.add_argument("--cluster-gap-threshold", type=float, default=DEFAULT_CLUSTER_GAP_THRESHOLD)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--force", action="store_true", help="Allow replacement of existing pilot source outputs.")
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    output_dir = repo_path(Path(ns.output_dir))
    n_shape_points = int(ns.n_shape_points)
    if bool(ns.smoke):
        n_shape_points = min(n_shape_points, SMOKE_N_SHAPE_POINTS)
        if Path(ns.output_dir) == DEFAULT_OUTPUT_DIR:
            output_dir = repo_path(SMOKE_OUTPUT_DIR)
    cache_dir = repo_path(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / "cache"
    args = Args(
        manifest=repo_path(Path(ns.manifest)),
        output_dir=output_dir,
        cache_dir=cache_dir,
        k_max=int(ns.k_max),
        n_candidate_roots=int(ns.n_candidate_roots),
        n_shape_points=n_shape_points,
        cluster_gap_threshold=float(ns.cluster_gap_threshold),
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        force=bool(ns.force),
        smoke=bool(ns.smoke),
    )
    if args.k_max != 10:
        raise ValueError("this pilot contract requires --k-max 10")
    if args.n_candidate_roots < args.k_max + 1:
        raise ValueError("--n-candidate-roots must include the K+1 EB guard root")
    if args.n_shape_points < 51:
        raise ValueError("--n-shape-points must be at least 51")
    if args.cluster_gap_threshold <= 0.0:
        raise ValueError("--cluster-gap-threshold must be positive")
    return args


def parse_groups(value: object) -> tuple[str, ...]:
    return tuple(part.strip() for part in str(value).split(";") if part.strip())


def load_manifest(path: Path) -> list[PilotCase]:
    if not path.exists():
        raise FileNotFoundError(f"pilot manifest does not exist: {path}")
    with path.open("r", newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    required = {"case_id", "epsilon_0", "beta_deg", "mu", "eta", "case_groups", "notes"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"pilot manifest must contain columns: {sorted(required)}")
    cases: list[PilotCase] = []
    for row in rows:
        cases.append(
            PilotCase(
                case_id=str(row["case_id"]).strip(),
                epsilon_0=float(row["epsilon_0"]),
                beta_deg=float(row["beta_deg"]),
                mu=float(row["mu"]),
                eta=float(row["eta"]),
                case_groups=parse_groups(row["case_groups"]),
                notes=str(row.get("notes", "")).strip(),
            )
        )
    return cases


def local_parameters(case: PilotCase) -> dict[str, float]:
    return stage1.local_thickness_parameters(case.epsilon_0, case.mu, case.eta)


def group_members(cases: Sequence[PilotCase], group: str) -> set[str]:
    return {case.case_id for case in cases if group in case.case_groups}


def validate_manifest(cases: Sequence[PilotCase], *, require_full_manifest: bool = True) -> None:
    if require_full_manifest and len(cases) != 21:
        raise ValueError(f"the production pilot manifest must contain exactly 21 rows, found {len(cases)}")
    if not cases:
        raise ValueError("pilot manifest is empty")
    case_ids = [case.case_id for case in cases]
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("pilot case_id values must be unique")
    keys = [certification.geometry_key(case.epsilon_0, case.beta_deg, case.mu, case.eta) for case in cases]
    if len(keys) != len(set(keys)):
        raise ValueError("pilot geometry keys must be unique")
    for case in cases:
        values = (case.epsilon_0, case.beta_deg, case.mu, case.eta)
        if not all(isfinite(value) for value in values):
            raise ValueError(f"non-finite parameter in case {case.case_id}")
        if case.epsilon_0 <= 0.0:
            raise ValueError(f"epsilon_0 must be positive in case {case.case_id}")
        if not 0.0 <= case.beta_deg <= 90.0:
            raise ValueError(f"beta_deg must lie in [0, 90] in case {case.case_id}")
        local_parameters(case)

    if not require_full_manifest:
        return
    expected_ids = {
        *(f"B{index:02d}" for index in range(1, 8)),
        *(f"G{index:02d}" for index in range(1, 9)),
        "H01",
        "H02",
        *(f"M{index:02d}" for index in range(1, 5)),
    }
    if set(case_ids) != expected_ids:
        raise ValueError("pilot manifest case IDs differ from the specified 21-case design")
    expected_groups = {
        "baseline_reference": {f"B{index:02d}" for index in range(1, 8)},
        "fixed_epsilon_0p02_block": {"B03", *(f"G{index:02d}" for index in range(1, 9))},
        "matched_epsilon_max_0p02": {"G03", "M01", "M02"},
        "matched_epsilon_max_0p04": {"M03", "G04", "M04"},
        "eta_sign_triplet": {"H01", "G04", "H02"},
    }
    for group, expected in expected_groups.items():
        observed = group_members(cases, group)
        if observed != expected:
            raise ValueError(f"group {group} has {sorted(observed)}, expected {sorted(expected)}")
    if group_members(cases, "eta_sign_triplet_center") != {"G04"}:
        raise ValueError("G04 must be the unique eta_sign_triplet_center case")
    by_id = {case.case_id: case for case in cases}
    for group, target in (("matched_epsilon_max_0p02", 0.02), ("matched_epsilon_max_0p04", 0.04)):
        for case_id in expected_groups[group]:
            observed = local_parameters(by_id[case_id])["epsilon_max"]
            if not math.isclose(observed, target, rel_tol=0.0, abs_tol=1.0e-12):
                raise ValueError(f"{case_id} has epsilon_max={observed}, expected {target}")


def select_cases(cases: Sequence[PilotCase], *, smoke: bool) -> list[PilotCase]:
    if not smoke:
        return list(cases)
    by_id = {case.case_id: case for case in cases}
    selected = [by_id[case_id] for case_id in SMOKE_CASE_IDS if case_id in by_id]
    return selected if selected else list(cases[: min(3, len(cases))])


def ensure_output_policy(args: Args) -> None:
    existing = [args.output_dir / name for name in SOURCE_OUTPUT_NAMES if (args.output_dir / name).exists()]
    if existing and not args.force:
        names = ", ".join(path.name for path in existing)
        raise FileExistsError(f"pilot source outputs already exist ({names}); pass --force to replace them")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def fixed_args_for_case(args: Args, case: PilotCase) -> fixed_scan.Args:
    return fixed_scan.Args(
        epsilon=case.epsilon_0,
        beta_values=(case.beta_deg,),
        mu_values=(case.mu,),
        eta_values=(case.eta,),
        n_reported_modes=args.k_max,
        n_candidate_roots=args.n_candidate_roots,
        n_shape_points=args.n_shape_points,
        cluster_gap_threshold=args.cluster_gap_threshold,
        benchmark_local_timo=False,
        workers=1,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        reuse_cache=args.reuse_cache,
        force_recompute=args.force_recompute,
        plot_only=False,
        smoke=args.smoke,
    )


def runtime_snapshot(runtime: fixed_scan.RuntimeTracker) -> dict[str, float | int]:
    return {
        name: getattr(runtime, name)
        for name in runtime.__dataclass_fields__
    }


def runtime_difference(before: Mapping[str, float | int], after: Mapping[str, float | int]) -> dict[str, float | int]:
    return {name: after[name] - before[name] for name in before}


def gap_from_pair(left: float, right: float) -> float:
    if not (isfinite(left) and isfinite(right)):
        return float("nan")
    return fixed_scan.sorted_gap((left, right), 1)


def root_quality_reasons(
    case: PilotCase,
    sorted_rows: Sequence[Mapping[str, object]],
    matching_rows: Sequence[Mapping[str, object]],
    eb_entry: fixed_scan.RootCacheEntry,
    timo_entry: fixed_scan.RootCacheEntry,
    *,
    k_max: int,
    n_candidate_roots: int,
) -> list[str]:
    reasons: list[str] = []
    if len(sorted_rows) != int(k_max):
        reasons.append(f"missing_sorted_rows={len(sorted_rows)}/{int(k_max)}")
    indices = [int(float(row.get("eb_sorted_index", 0))) for row in sorted_rows]
    if len(indices) != len(set(indices)):
        reasons.append("duplicate_sorted_indices")
    if indices != list(range(1, int(k_max) + 1)):
        reasons.append("noncanonical_sorted_indices")
    for label, entry in (("EB", eb_entry), ("Timoshenko", timo_entry)):
        if entry.root_count_found < int(n_candidate_roots):
            reasons.append(f"missing_{label}_candidate_roots={entry.root_count_found}/{int(n_candidate_roots)}")
        roots = list(entry.roots[: int(k_max)])
        if len(roots) < int(k_max):
            reasons.append(f"missing_{label}_roots")
        elif not all(isfinite(value) for value in roots):
            reasons.append(f"nonfinite_{label}_roots")
        elif any(right <= left for left, right in zip(roots, roots[1:])):
            reasons.append(f"nonincreasing_{label}_roots")
        if entry.warnings:
            reasons.append(f"{label}_root_warning")
    selected = [row for row in matching_rows if str(row.get("assignment_selected", "")).lower() == "true"]
    if any(str(row.get("matching_status", "")) == "candidate_boundary" for row in selected):
        reasons.append("candidate_boundary_warning")
    inconsistent: list[int] = []
    for row in sorted_rows:
        recomputed = certification.squared_lambda_delta(
            finite_float(row, "Lambda_EB"), finite_float(row, "Lambda_Timo")
        )
        stored = finite_float(row, "delta_f")
        tolerance = certification.DELTA_ABS_TOLERANCE + certification.DELTA_REL_TOLERANCE * max(
            abs(recomputed), abs(stored), 1.0e-12
        )
        if not (isfinite(recomputed) and isfinite(stored)) or abs(recomputed - stored) > tolerance:
            inconsistent.append(int(float(row.get("eb_sorted_index", 0))))
    if inconsistent:
        reasons.append("delta_f_mismatch_at=" + ",".join(str(value) for value in inconsistent))
    del case
    return list(dict.fromkeys(reasons))


def solve_case(
    args: Args,
    case: PilotCase,
    runtime: fixed_scan.RuntimeTracker,
) -> dict[str, object]:
    point_args = fixed_args_for_case(args, case)
    cache = fixed_scan.RootCache(point_args, runtime)
    before = runtime_snapshot(runtime)
    products = fixed_scan.point_products(
        point_args,
        cache,
        runtime,
        beta_deg=case.beta_deg,
        mu=case.mu,
        eta=case.eta,
    )
    after = runtime_snapshot(runtime)
    operation_delta = runtime_difference(before, after)

    audit_cache = fixed_scan.RootCache(point_args, None)
    eb_entry, timo_entry = fixed_scan.solve_point_roots(
        point_args,
        audit_cache,
        beta_deg=case.beta_deg,
        mu=case.mu,
        eta=case.eta,
        n_roots=args.n_candidate_roots,
    )
    sorted_rows = sorted(
        (row for row in products.mode_rows if row.get("comparison_type") == "sorted_index"),
        key=lambda row: int(row["eb_sorted_index"]),
    )
    homologous = {
        int(row["eb_sorted_index"]): row
        for row in products.mode_rows
        if row.get("comparison_type") == "homologous_mode"
        and isinstance(row.get("eb_sorted_index"), int)
    }
    local = local_parameters(case)
    geometry_key = certification.geometry_key(case.epsilon_0, case.beta_deg, case.mu, case.eta)
    geometry_id = certification.geometry_identifier(geometry_key)
    case_groups = ";".join(case.case_groups)
    reasons = root_quality_reasons(
        case,
        sorted_rows,
        products.matching_rows,
        eb_entry,
        timo_entry,
        k_max=args.k_max,
        n_candidate_roots=args.n_candidate_roots,
    )
    included = not reasons
    eb_roots = list(eb_entry.roots)
    guard_11 = eb_roots[args.k_max] if len(eb_roots) > args.k_max and isfinite(eb_roots[args.k_max]) else float("nan")
    mode_rows: list[dict[str, object]] = []
    deltas: list[float] = []
    for source in sorted_rows:
        index = int(source["eb_sorted_index"])
        lambda_eb = finite_float(source, "Lambda_EB")
        lambda_timo = finite_float(source, "Lambda_Timo")
        delta = certification.squared_lambda_delta(lambda_eb, lambda_timo)
        deltas.append(delta)
        left_gap = gap_from_pair(eb_roots[index - 2], eb_roots[index - 1]) if index > 1 and len(eb_roots) >= index else float("nan")
        right_gap = gap_from_pair(eb_roots[index - 1], eb_roots[index]) if len(eb_roots) > index else float("nan")
        gaps = [value for value in (left_gap, right_gap) if isfinite(value)]
        matched = homologous.get(index, {})
        stored = finite_float(source, "delta_f")
        tolerance = certification.DELTA_ABS_TOLERANCE + certification.DELTA_REL_TOLERANCE * max(
            abs(delta), abs(stored), 1.0e-12
        )
        source_warning = ";".join(
            item
            for item in (
                "root_warning" if str(source.get("root_warning", "")).lower() == "true" else "",
                str(matched.get("matching_status", "")) if str(matched.get("matching_status", "")) not in {"", "reliable_individual", "reliable_cluster"} else "",
            )
            if item
        )
        mode_rows.append(
            {
                "case_id": case.case_id,
                "geometry_id": geometry_id,
                "epsilon_0": case.epsilon_0,
                "beta_deg": case.beta_deg,
                "mu": case.mu,
                "eta": case.eta,
                "case_groups": case_groups,
                "sorted_index": index,
                "Lambda_EB": lambda_eb,
                "Lambda_Timo": lambda_timo,
                "delta_f": delta,
                "stored_delta_f": stored,
                "delta_f_consistent": isfinite(delta) and isfinite(stored) and abs(delta - stored) <= tolerance,
                "true_pass_10": delta <= FREQUENCY_ERROR_THRESHOLD,
                **local,
                **{field: source.get(field, "") for field in (
                    "chi_max_EB",
                    "chi_eff_EB",
                    "Theta_max_EB",
                    "Pi_shear_EB",
                    "Pi_rotary_EB",
                    "Pi_EB",
                    "EB_axial_fraction",
                    "EB_bending_fraction",
                    "EB_classification",
                    "Timo_axial_fraction",
                    "Timo_bending_fraction",
                    "Timo_shear_fraction",
                    "Timo_bending_shear_fraction",
                    "Timo_classification",
                )},
                "gap_left_EB": left_gap,
                "gap_right_EB": right_gap,
                "gap_EB": min(gaps) if gaps else float("nan"),
                "gap_sidedness": "right_only" if index == 1 else ("two_sided_with_guard" if index == args.k_max and isfinite(right_gap) else ("left_only_guard_unavailable" if index == args.k_max else "two_sided")),
                "Lambda_EB_guard_11": guard_11 if index == args.k_max else "",
                "mode10_upper_guard_source": ("candidate_root_11" if isfinite(guard_11) else "unavailable") if index == args.k_max else "not_applicable",
                "matched_timo_sorted_index": matched.get("timo_sorted_index", ""),
                "shape_MAC_EB_Timo": matched.get("MAC_uw", ""),
                "matching_status": matched.get("matching_status", ""),
                "cluster_id": matched.get("cluster_id", ""),
                "cluster_size_EB": matched.get("cluster_size_EB", ""),
                "cluster_size_Timo": matched.get("cluster_size_Timo", ""),
                "root_warning": source.get("root_warning", ""),
                "candidate_boundary": matched.get("candidate_boundary", ""),
                "included_in_certification": included,
                "quality_status": "included" if included else "excluded",
                "source_warning": source_warning,
                "notes": "sorted-spectrum target; homologous matching fields are quality metadata only",
            }
        )

    complete = len(deltas) == args.k_max and all(isfinite(value) for value in deltas)
    n_true = certification.true_safe_prefix(deltas, k_max=args.k_max, strict=True) if complete else None
    passes = [value <= FREQUENCY_ERROR_THRESHOLD for value in deltas]
    first_failure_position = next((position for position, value in enumerate(passes) if not value), None)
    late_pass = bool(first_failure_position is not None and any(passes[first_failure_position + 1 :]))
    selected_matching = [
        row for row in products.matching_rows if str(row.get("assignment_selected", "")).lower() == "true"
    ]
    root_warning_count = int(bool(eb_entry.warnings)) + int(bool(timo_entry.warnings))
    boundary_count = sum(str(row.get("matching_status", "")) == "candidate_boundary" for row in selected_matching)
    ambiguous_count = sum(str(row.get("matching_status", "")) == "ambiguous" for row in selected_matching)
    cluster_ids = {
        str(row.get("cluster_EB", ""))
        for row in selected_matching
        if str(row.get("cluster_EB", ""))
    }
    geometry_row = {
        "case_id": case.case_id,
        "geometry_id": geometry_id,
        "epsilon_0": case.epsilon_0,
        "beta_deg": case.beta_deg,
        "mu": case.mu,
        "eta": case.eta,
        "epsilon_1": local["epsilon_1"],
        "epsilon_2": local["epsilon_2"],
        "epsilon_max": local["epsilon_max"],
        "case_groups": case_groups,
        "N_true": n_true if n_true is not None else "",
        "first_true_failed_mode": n_true + 1 if n_true is not None and n_true < args.k_max else "",
        "pass_pattern_10": ",".join("Y" if value else "N" for value in passes),
        "late_pass_after_failure": late_pass,
        "complete_K10": complete,
        "quality_status": "included" if included else "excluded",
        "root_warning_count": root_warning_count,
        "candidate_boundary_count": boundary_count,
        "ambiguous_matching_count": ambiguous_count,
        "close_cluster_count": len(cluster_ids),
        "root11_guard_available": isfinite(guard_11),
        "exclusion_reason": ";".join(reasons),
        "notes": "N_true uses the first continuous sorted-frequency prefix",
    }
    matching_rows = [
        {
            "case_id": case.case_id,
            "geometry_id": geometry_id,
            "epsilon_0": case.epsilon_0,
            "case_groups": case_groups,
            **row,
        }
        for row in products.matching_rows
    ]
    exclusion_row = None if included else {
        "case_id": case.case_id,
        "geometry_id": geometry_id,
        "epsilon_0": case.epsilon_0,
        "beta_deg": case.beta_deg,
        "mu": case.mu,
        "eta": case.eta,
        "complete_K10": complete,
        "root_warning_count": root_warning_count,
        "candidate_boundary_count": boundary_count,
        "exclusion_reason": ";".join(reasons),
        "notes": "excluded explicitly; no silent removal",
    }
    operation_row = {
        "case_id": case.case_id,
        "cost_scope": "offline_reference_generation",
        "measurement_status": "measured_wrapper_counters; determinant_primitives_unavailable",
        "n_candidate_roots": args.n_candidate_roots,
        "n_shape_points": args.n_shape_points,
        "EB_root_solver_calls": operation_delta["eb_root_calls"],
        "Timoshenko_root_solver_calls": operation_delta["timo_root_calls"],
        "root_cache_hits": operation_delta["root_cache_hits"],
        "root_cache_misses": operation_delta["root_cache_misses"],
        "boundary_retries": operation_delta["boundary_retries"],
        "EB_root_seconds": operation_delta["eb_roots_seconds"],
        "Timoshenko_root_seconds": operation_delta["timo_roots_seconds"],
        "shape_and_predictor_seconds": operation_delta["eb_shape_pi_seconds"],
        "matching_seconds": operation_delta["matching_seconds"],
        "matching_pair_rows": len(products.matching_rows),
        "EB_characteristic_matrix_evaluations": "unavailable",
        "EB_determinant_evaluations": "unavailable",
        "Timoshenko_characteristic_matrix_evaluations": "unavailable",
        "Timoshenko_determinant_evaluations": "unavailable",
        "root_bracketing_intervals": "unavailable",
        "Brent_function_evaluations": "unavailable",
        "notes": "actual existing-wrapper counters only; no counts inferred from wall-clock time",
    }
    resolved_row = {
        "case_id": case.case_id,
        "geometry_id": geometry_id,
        "epsilon_0": case.epsilon_0,
        "beta_deg": case.beta_deg,
        "mu": case.mu,
        "eta": case.eta,
        **local,
        "case_groups": case_groups,
        "notes": case.notes,
    }
    return {
        "resolved_row": resolved_row,
        "mode_rows": mode_rows,
        "geometry_row": geometry_row,
        "matching_rows": matching_rows,
        "exclusion_row": exclusion_row,
        "operation_row": operation_row,
    }


def write_generation_report(
    args: Args,
    cases: Sequence[PilotCase],
    geometry_rows: Sequence[Mapping[str, object]],
    operation_rows: Sequence[Mapping[str, object]],
    *,
    elapsed_seconds: float,
) -> Path:
    included = [row for row in geometry_rows if row.get("quality_status") == "included"]
    excluded = [row for row in geometry_rows if row.get("quality_status") != "included"]
    root_warnings = sum(int(row.get("root_warning_count", 0)) for row in geometry_rows)
    boundary = sum(int(row.get("candidate_boundary_count", 0)) for row in geometry_rows)
    cache_hits = sum(int(row.get("root_cache_hits", 0)) for row in operation_rows)
    cache_misses = sum(int(row.get("root_cache_misses", 0)) for row in operation_rows)
    report = args.output_dir / "epsilon_pilot_generation_report.md"
    lines = [
        "# EB Epsilon A-Priori Pilot Generation Report",
        "",
        "## Scope",
        "",
        f"- Manifest geometries solved: `{len(cases)}`.",
        f"- Sorted target modes: `K = {args.k_max}`.",
        f"- Candidate roots per theory: `{args.n_candidate_roots}`.",
        f"- Shape samples per rod: `{args.n_shape_points}`.",
        "- The dimensional-frequency discrepancy is recomputed from squared `Lambda`.",
        "- Timoshenko is used only as the one-dimensional reference target and matching audit.",
        "- No FEM, 3D FEM, Gmsh, or CalculiX workflow was run.",
        "",
        "## Quality",
        "",
        f"- Included geometries: `{len(included)}`.",
        f"- Excluded geometries: `{len(excluded)}`.",
        f"- Root-warning count: `{root_warnings}`.",
        f"- Candidate-boundary count: `{boundary}`.",
        f"- Root cache hits/misses: `{cache_hits}/{cache_misses}`.",
        f"- Generation wall-clock time: `{elapsed_seconds:.6f}` seconds (auxiliary only).",
        "",
        "Excluded cases are retained in `epsilon_pilot_exclusion_audit.csv`; they are never removed silently.",
        "The operation CSV reports only counters exposed by the existing wrappers. Missing determinant-level counters are marked unavailable rather than reconstructed from elapsed time.",
    ]
    if excluded:
        lines.extend(["", "## Exclusions", ""])
        for row in excluded:
            lines.append(f"- `{row.get('case_id')}`: {row.get('exclusion_reason')}")
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    all_cases = load_manifest(args.manifest)
    validate_manifest(all_cases, require_full_manifest=not args.smoke or args.manifest == DEFAULT_MANIFEST)
    cases = select_cases(all_cases, smoke=args.smoke)
    ensure_output_policy(args)

    runtime = fixed_scan.RuntimeTracker()
    resolved_rows: list[dict[str, object]] = []
    mode_rows: list[dict[str, object]] = []
    geometry_rows: list[dict[str, object]] = []
    matching_rows: list[dict[str, object]] = []
    exclusion_rows: list[dict[str, object]] = []
    operation_rows: list[dict[str, object]] = []
    started = time.perf_counter()
    for case in cases:
        product = solve_case(args, case, runtime)
        resolved_rows.append(product["resolved_row"])  # type: ignore[arg-type]
        mode_rows.extend(product["mode_rows"])  # type: ignore[arg-type]
        geometry_rows.append(product["geometry_row"])  # type: ignore[arg-type]
        matching_rows.extend(product["matching_rows"])  # type: ignore[arg-type]
        if product["exclusion_row"] is not None:
            exclusion_rows.append(product["exclusion_row"])  # type: ignore[arg-type]
        operation_rows.append(product["operation_row"])  # type: ignore[arg-type]
        print(f"completed {case.case_id}: N_true={product['geometry_row']['N_true']}, status={product['geometry_row']['quality_status']}")  # type: ignore[index]
    elapsed = time.perf_counter() - started

    write_csv(args.output_dir / "epsilon_pilot_case_manifest_resolved.csv", resolved_rows, RESOLVED_FIELDS)
    write_csv(args.output_dir / "epsilon_pilot_mode_metrics.csv", mode_rows, MODE_FIELDS)
    write_csv(args.output_dir / "epsilon_pilot_geometry_metrics.csv", geometry_rows, GEOMETRY_FIELDS)
    write_csv(args.output_dir / "epsilon_pilot_matching_audit.csv", matching_rows, MATCHING_FIELDS)
    write_csv(args.output_dir / "epsilon_pilot_exclusion_audit.csv", exclusion_rows, EXCLUSION_FIELDS)
    write_csv(args.output_dir / "epsilon_pilot_root_operation_counts.csv", operation_rows, ROOT_OPERATION_FIELDS)
    report = write_generation_report(
        args,
        cases,
        geometry_rows,
        operation_rows,
        elapsed_seconds=elapsed,
    )
    counts = Counter(str(row["quality_status"]) for row in geometry_rows)
    print(f"pilot geometries: {len(cases)}; included/excluded: {counts['included']}/{counts['excluded']}")
    print(f"root cache hits/misses: {runtime.root_cache_hits}/{runtime.root_cache_misses}")
    print(f"output: {args.output_dir}")
    return {
        "args": args,
        "cases": cases,
        "resolved_rows": resolved_rows,
        "mode_rows": mode_rows,
        "geometry_rows": geometry_rows,
        "matching_rows": matching_rows,
        "exclusion_rows": exclusion_rows,
        "operation_rows": operation_rows,
        "report": report,
        "elapsed_seconds": elapsed,
    }


if __name__ == "__main__":
    main()
