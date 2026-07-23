from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import statistics
import sys
import time
from typing import Iterator, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
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
from scripts.analysis.thickness_mismatch.maps import (  # noqa: E402
    plot_eb_vs_timoshenko_lambda_beta_cases as beta_workflow,
)
from scripts.analysis.thickness_mismatch.postprocess import (  # noqa: E402
    analyze_eb_safe_prefix_certification as certification,
)
from scripts.lib import variable_length_timoshenko as TIMO  # noqa: E402


ALGORITHM_VERSION = "rule_s_cost_break_even_v2"
K_MAX = 10
EXPECTED_CASE_IDS = ("B01", "G04", "S3_06", "S3_12", "S3_14")
DEFAULT_PHASE1_DIR = REPO_ROOT / "results" / "eb_rule_ab_exact_pareto"
DEFAULT_PROPOSAL_CSV = DEFAULT_PHASE1_DIR / "rule_B_cost_break_even_proposal.csv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "eb_rule_s_cost_break_even"
DEFAULT_RUNTIME_REPEATS = 3
DEFAULT_N_SHAPE_POINTS = 401
REFERENCE_ABS_TOL = 1.0e-6
REFERENCE_REL_TOL = 1.0e-8
CACHE_NAME = "benchmark_cache.json"
REPORT_NAME = "rule_S_cost_break_even_report.md"
PLOT_NAMES = ("runtime_comparison.png", "timo_operation_comparison.png")

OUTPUT_FIELDS: dict[str, tuple[str, ...]] = {
    "cost_benchmark_manifest.csv": (
        "algorithm_version", "proposal_status", "execution_status", "case_id", "source", "partition",
        "epsilon_0", "beta_deg", "mu", "eta", "frozen_T_s", "N_true", "N_S", "suffix_frequency_count",
        "eb_root_source", "timo_reference_source", "online_timo_reference_policy", "purpose",
    ),
    "cost_benchmark_input_audit.csv": (
        "algorithm_version", "case_id", "proposal_csv", "prefix_csv", "rule_S_predictions_csv",
        "proposal_sha256", "prefix_sha256", "rule_S_predictions_sha256", "proposal_geometry_match",
        "sorted_EB_modes_1_to_10", "saved_Timo_reference_count", "online_case_contains_Timo_references",
        "online_bracketing_uses_saved_Timo_references", "input_quality_status", "notes",
    ),
    "direct_timo_predictions.csv": (
        "algorithm_version", "case_id", "mode_index", "Lambda_Timo_direct", "Lambda_Timo_reference",
        "absolute_error", "relative_error", "verification_passed", "root_source", "warnings",
    ),
    "hybrid_rule_S_predictions.csv": (
        "algorithm_version", "case_id", "mode_index", "N_S", "accepted_by_rule_S", "frequency_source",
        "Lambda_hybrid", "Lambda_reference", "absolute_error", "relative_error", "verification_passed",
        "fallback_used", "local_attempt_notes",
    ),
    "suffix_root_verification.csv": (
        "algorithm_version", "case_id", "mode_index", "online_root_source", "Lambda_online",
        "Lambda_saved_reference", "absolute_error", "relative_error", "verification_tolerance",
        "verification_passed", "verification_performed_after_online_solve", "reference_used_for_bracketing",
    ),
    "per_case_operation_counts.csv": (
        "algorithm_version", "case_id", "workflow", "canonical_repeat", "EB_mode_reconstructions",
        "EB_characteristic_matrix_evaluations", "EB_SVD_6x6_calls", "EB_shape_grid_point_evaluations",
        "EB_quadrature_calls", "EB_quadrature_point_evaluations", "selector_scalar_comparisons",
        "Timo_root_solver_calls", "Timo_roots_requested", "Timo_roots_returned_by_solver_calls",
        "Timo_characteristic_determinant_evaluations", "Timo_characteristic_matrix_evaluations",
        "Timo_scan_grid_evaluations", "Timo_brent_calls", "Timo_brent_callback_evaluations",
        "Timo_minimizer_callback_evaluations", "Timo_SVD_calls", "strict_root_quality_checks",
        "full_spectrum_fallback_count", "suffix_frequency_count", "Timoshenko_frequency_count", "warnings",
    ),
    "aggregate_operation_counts.csv": (
        "algorithm_version", "workflow", "case_count", "suffix_case_count", "EB_mode_reconstructions",
        "EB_characteristic_matrix_evaluations", "EB_SVD_6x6_calls", "EB_shape_grid_point_evaluations",
        "EB_quadrature_calls", "EB_quadrature_point_evaluations", "selector_scalar_comparisons",
        "Timo_root_solver_calls", "Timo_roots_requested", "Timo_roots_returned_by_solver_calls",
        "Timo_characteristic_determinant_evaluations", "Timo_characteristic_matrix_evaluations",
        "Timo_scan_grid_evaluations", "Timo_brent_calls", "Timo_brent_callback_evaluations",
        "Timo_minimizer_callback_evaluations", "Timo_SVD_calls", "strict_root_quality_checks",
        "full_spectrum_fallback_count", "suffix_frequency_count", "Timoshenko_frequency_count",
    ),
    "runtime_measurements.csv": (
        "algorithm_version", "case_id", "workflow", "repeat_index", "temperature", "elapsed_seconds",
        "result_hash", "N_S", "suffix_frequency_count", "full_spectrum_fallback_count",
        "Timo_root_solver_calls", "Timo_roots_returned_by_solver_calls",
        "Timo_characteristic_determinant_evaluations", "Timo_SVD_calls", "Timoshenko_frequency_count",
        "EB_mode_reconstructions", "EB_SVD_6x6_calls", "EB_quadrature_calls",
    ),
    "cost_break_even_summary.csv": (
        "algorithm_version", "status", "case_count", "runtime_repeats", "all_reference_verifications_passed",
        "suffix_case_count", "suffix_cases_with_full_fallback", "systematic_full_fallback",
        "direct_total_Timo_determinant_evaluations", "hybrid_total_Timo_determinant_evaluations",
        "direct_total_Timo_root_solver_calls", "hybrid_total_Timo_root_solver_calls",
        "hybrid_total_EB_mode_reconstructions", "direct_median_runtime_seconds",
        "hybrid_median_runtime_seconds", "common_scalar_cost_model_used", "decision_basis",
        "new_EB_root_calculations", "new_Timo_root_solver_calls_executed",
        "new_Timo_roots_returned_during_all_repeats", "new_Timo_determinant_evaluations_executed",
        "new_EB_mode_reconstructions_executed", "new_EB_SVD_6x6_calls_executed",
        "new_EB_quadrature_calls_executed", "cache_reused",
    ),
}


@dataclass(frozen=True)
class Args:
    proposal_csv: Path
    output_dir: Path
    force: bool
    reuse_cache: bool
    plot_only: bool
    runtime_repeats: int
    n_shape_points: int


@dataclass(frozen=True)
class CaseRecord:
    case_id: str
    source: str
    partition: str
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float
    frozen_T_s: float
    N_true: int
    expected_N_S: int
    eb_roots: tuple[float, ...]
    timo_references: tuple[float, ...]
    eb_root_source: str


@dataclass(frozen=True)
class OnlineCase:
    """Inputs permitted during an online workflow; intentionally has no Timoshenko references."""

    case_id: str
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float
    eb_roots: tuple[float, ...]


@dataclass
class OperationCounts:
    EB_mode_reconstructions: int = 0
    EB_characteristic_matrix_evaluations: int = 0
    EB_SVD_6x6_calls: int = 0
    EB_shape_grid_point_evaluations: int = 0
    EB_quadrature_calls: int = 0
    EB_quadrature_point_evaluations: int = 0
    selector_scalar_comparisons: int = 0
    Timo_root_solver_calls: int = 0
    Timo_roots_requested: int = 0
    Timo_roots_returned_by_solver_calls: int = 0
    Timo_characteristic_determinant_evaluations: int = 0
    Timo_characteristic_matrix_evaluations: int = 0
    Timo_scan_grid_evaluations: int = 0
    Timo_brent_calls: int = 0
    Timo_brent_callback_evaluations: int = 0
    Timo_minimizer_callback_evaluations: int = 0
    Timo_SVD_calls: int = 0
    strict_root_quality_checks: int = 0
    full_spectrum_fallback_count: int = 0
    suffix_frequency_count: int = 0
    Timoshenko_frequency_count: int = 0
    warnings: list[str] | None = None

    def __post_init__(self) -> None:
        if self.warnings is None:
            self.warnings = []


@dataclass(frozen=True)
class HybridResult:
    N_S: int
    pi_shear: tuple[float, ...]
    suffix_roots: tuple[float, ...]
    suffix_sources: tuple[str, ...]
    fallback_used: bool
    local_notes: tuple[str, ...]


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Benchmark frozen Rule S against a direct Timoshenko K=10 solve on five proposal cases.",
    )
    parser.add_argument("--proposal-csv", type=Path, default=DEFAULT_PROPOSAL_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--runtime-repeats", type=int, default=DEFAULT_RUNTIME_REPEATS)
    parser.add_argument("--n-shape-points", type=int, default=DEFAULT_N_SHAPE_POINTS)
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if int(ns.runtime_repeats) < 2:
        raise ValueError("--runtime-repeats must be at least 2")
    if int(ns.n_shape_points) < 51:
        raise ValueError("--n-shape-points must be at least 51")
    return Args(
        proposal_csv=repo_path(Path(ns.proposal_csv)),
        output_dir=repo_path(Path(ns.output_dir)),
        force=bool(ns.force),
        reuse_cache=bool(ns.reuse_cache),
        plot_only=bool(ns.plot_only),
        runtime_repeats=int(ns.runtime_repeats),
        n_shape_points=int(ns.n_shape_points),
    )


def fmt(value: object) -> object:
    return certification.fmt(value)


def read_csv(path: Path) -> list[dict[str, str]]:
    return certification.read_csv(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stable_hash(values: Sequence[float]) -> str:
    payload = "|".join(f"{float(value):.17e}" for value in values)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})


def _float(row: Mapping[str, object], key: str) -> float:
    return certification.finite_float(row, key)


def _int(row: Mapping[str, object], key: str) -> int:
    value = _float(row, key)
    if not float(value).is_integer():
        raise ValueError(f"{key} must be an integer")
    return int(value)


def load_cases(proposal_csv: Path) -> tuple[list[CaseRecord], list[dict[str, object]], list[dict[str, object]]]:
    phase1_dir = proposal_csv.parent
    prefix_csv = phase1_dir / "prefix_predictor_metrics.csv"
    predictions_csv = phase1_dir / "rule_S_predictions.csv"
    proposal = read_csv(proposal_csv)
    prefix = read_csv(prefix_csv)
    predictions = read_csv(predictions_csv)
    if tuple(row.get("case_id", "") for row in proposal) != EXPECTED_CASE_IDS:
        raise ValueError(f"proposal must contain exactly {EXPECTED_CASE_IDS} in that order")
    if any(row.get("proposal_status") != "proposal_only" or row.get("execution_status") != "not_executed" for row in proposal):
        raise ValueError("proposal rows must be proposal_only/not_executed")

    prefix_by_case: dict[str, list[dict[str, str]]] = {}
    for row in prefix:
        prefix_by_case.setdefault(row["case_id"], []).append(row)
    prediction_by_case = {row["case_id"]: row for row in predictions}
    proposal_hash = sha256_file(proposal_csv)
    prefix_hash = sha256_file(prefix_csv)
    prediction_hash = sha256_file(predictions_csv)
    cases: list[CaseRecord] = []
    manifest: list[dict[str, object]] = []
    input_audit: list[dict[str, object]] = []
    thresholds: set[float] = set()
    for proposal_row in proposal:
        case_id = proposal_row["case_id"]
        rows = sorted(prefix_by_case.get(case_id, []), key=lambda row: _int(row, "prefix_n"))
        if [int(row["prefix_n"]) for row in rows] != list(range(1, K_MAX + 1)):
            raise ValueError(f"{case_id}: prefix CSV does not contain sorted modes 1..10")
        pred = prediction_by_case.get(case_id)
        if pred is None:
            raise ValueError(f"{case_id}: missing Rule S prediction")
        geometry = tuple(_float(rows[0], key) for key in ("epsilon_0", "beta_deg", "mu", "eta"))
        proposal_geometry = tuple(_float(proposal_row, key) for key in ("epsilon_0", "beta_deg", "mu", "eta"))
        geometry_match = all(abs(left - right) <= 1.0e-12 for left, right in zip(geometry, proposal_geometry, strict=True))
        if not geometry_match:
            raise ValueError(f"{case_id}: proposal geometry differs from Phase-I prefix data")
        eb_roots = tuple(_float(row, "Lambda_EB") for row in rows)
        refs = tuple(_float(row, "Lambda_Timo") for row in rows)
        if any(b <= a for a, b in zip(eb_roots, eb_roots[1:])):
            raise ValueError(f"{case_id}: EB roots are not strictly sorted")
        threshold = _float(pred, "T_s")
        thresholds.add(threshold)
        record = CaseRecord(
            case_id=case_id,
            source=rows[0]["source"],
            partition=rows[0]["assigned_partition"],
            epsilon_0=geometry[0],
            beta_deg=geometry[1],
            mu=geometry[2],
            eta=geometry[3],
            frozen_T_s=threshold,
            N_true=_int(pred, "N_true"),
            expected_N_S=_int(pred, "N_hat"),
            eb_roots=eb_roots,
            timo_references=refs,
            eb_root_source=rows[0]["root_source"],
        )
        cases.append(record)
        manifest.append({
            "algorithm_version": ALGORITHM_VERSION,
            "proposal_status": proposal_row["proposal_status"],
            "execution_status": "executed_targeted_benchmark",
            "case_id": case_id,
            "source": record.source,
            "partition": record.partition,
            "epsilon_0": record.epsilon_0,
            "beta_deg": record.beta_deg,
            "mu": record.mu,
            "eta": record.eta,
            "frozen_T_s": record.frozen_T_s,
            "N_true": record.N_true,
            "N_S": record.expected_N_S,
            "suffix_frequency_count": K_MAX - record.expected_N_S,
            "eb_root_source": record.eb_root_source,
            "timo_reference_source": record.eb_root_source,
            "online_timo_reference_policy": "withheld_until_post_solve_verification",
            "purpose": proposal_row["purpose"],
        })
        input_audit.append({
            "algorithm_version": ALGORITHM_VERSION,
            "case_id": case_id,
            "proposal_csv": str(proposal_csv.relative_to(REPO_ROOT)),
            "prefix_csv": str(prefix_csv.relative_to(REPO_ROOT)),
            "rule_S_predictions_csv": str(predictions_csv.relative_to(REPO_ROOT)),
            "proposal_sha256": proposal_hash,
            "prefix_sha256": prefix_hash,
            "rule_S_predictions_sha256": prediction_hash,
            "proposal_geometry_match": geometry_match,
            "sorted_EB_modes_1_to_10": True,
            "saved_Timo_reference_count": len(refs),
            "online_case_contains_Timo_references": False,
            "online_bracketing_uses_saved_Timo_references": False,
            "input_quality_status": "quality_approved_existing_phase1_data",
            "notes": "saved Timoshenko roots are attached only to CaseRecord for post-solve verification",
        })
    if len(thresholds) != 1:
        raise ValueError("the benchmark must use one frozen Rule S threshold")
    return cases, manifest, input_audit


def online_case(record: CaseRecord) -> OnlineCase:
    return OnlineCase(
        case_id=record.case_id,
        epsilon_0=record.epsilon_0,
        beta_deg=record.beta_deg,
        mu=record.mu,
        eta=record.eta,
        eb_roots=record.eb_roots,
    )


@contextmanager
def instrument_timo(counts: OperationCounts) -> Iterator[None]:
    original_det = TIMO.timo_det_for_scan
    original_timo_brentq = TIMO.brentq
    original_fixed_brentq = fixed_scan.brentq

    def counted_det(*args: object, **kwargs: object) -> float:
        counts.Timo_characteristic_determinant_evaluations += 1
        counts.Timo_characteristic_matrix_evaluations += 1
        return float(original_det(*args, **kwargs))

    def counted_brent_factory(original: object):
        def counted_brent(func: object, *args: object, **kwargs: object) -> float:
            counts.Timo_brent_calls += 1

            def callback(value: float) -> float:
                counts.Timo_brent_callback_evaluations += 1
                return float(func(value))  # type: ignore[operator]

            return float(original(callback, *args, **kwargs))  # type: ignore[operator]

        return counted_brent

    TIMO.timo_det_for_scan = counted_det  # type: ignore[assignment]
    TIMO.brentq = counted_brent_factory(original_timo_brentq)  # type: ignore[assignment]
    fixed_scan.brentq = counted_brent_factory(original_fixed_brentq)  # type: ignore[assignment]
    try:
        yield
    finally:
        TIMO.timo_det_for_scan = original_det  # type: ignore[assignment]
        TIMO.brentq = original_timo_brentq  # type: ignore[assignment]
        fixed_scan.brentq = original_fixed_brentq  # type: ignore[assignment]


def run_rule_s_selector(case: OnlineCase, threshold: float, n_points: int, counts: OperationCounts) -> tuple[int, tuple[float, ...]]:
    values: list[float] = []
    accepted = 0
    for mode_index, root in enumerate(case.eb_roots, start=1):
        mode = fixed_scan.eb_mode_result(
            epsilon=case.epsilon_0,
            beta_deg=case.beta_deg,
            mu=case.mu,
            eta=case.eta,
            sorted_index=mode_index,
            Lambda=root,
            n_points=n_points,
            root_warnings=(),
        )
        pi_shear = float(fixed_scan.eb_shear_predictor(mode)["Pi_shear_EB"])
        values.append(pi_shear)
        counts.EB_mode_reconstructions += 1
        counts.EB_characteristic_matrix_evaluations += 1
        counts.EB_SVD_6x6_calls += 1
        counts.EB_shape_grid_point_evaluations += 2 * n_points
        counts.EB_quadrature_calls += 6
        counts.EB_quadrature_point_evaluations += 6 * n_points
        counts.selector_scalar_comparisons += 1
        if pi_shear <= threshold:
            accepted = mode_index
        else:
            break
    return accepted, tuple(values)


def solve_direct_timo(case: OnlineCase, counts: OperationCounts) -> tuple[float, ...]:
    spec = beta_workflow.CaseSpec(mu=case.mu, eta=case.eta, epsilon=case.epsilon_0)
    before_det = counts.Timo_characteristic_determinant_evaluations
    before_brent = counts.Timo_brent_callback_evaluations
    counts.Timo_root_solver_calls += 1
    counts.Timo_roots_requested += K_MAX
    result = beta_workflow.solve_timo(spec, case.beta_deg, K_MAX)
    counts.Timo_roots_returned_by_solver_calls += result.root_count_found
    if result.root_count_found < K_MAX:
        counts.Timo_root_solver_calls += 1
        counts.Timo_roots_requested += K_MAX
        retry = beta_workflow.solve_timo(spec, case.beta_deg, K_MAX, retry=True)
        counts.Timo_roots_returned_by_solver_calls += retry.root_count_found
        if retry.root_count_found > result.root_count_found:
            result = retry
    counts.Timo_scan_grid_evaluations += (
        counts.Timo_characteristic_determinant_evaluations - before_det
        - (counts.Timo_brent_callback_evaluations - before_brent)
    )
    counts.warnings.extend(result.warnings)
    roots = tuple(float(value) for value in result.roots[:K_MAX] if math.isfinite(float(value)))
    if len(roots) != K_MAX or any(b <= a for a, b in zip(roots, roots[1:])):
        raise RuntimeError(f"{case.case_id}: direct Timoshenko K=10 solve did not return ten sorted roots")
    return roots


def local_interval(roots: Sequence[float], mode_index: int) -> tuple[float, float]:
    """Use the existing fixed-scan audit's 45%-of-neighbour-gap local interval policy."""

    index = mode_index - 1
    root = float(roots[index])
    left_neighbor = float(roots[index - 1]) if index > 0 else beta_workflow.ROOT_SCAN_START
    right_neighbor = float(roots[index + 1]) if index + 1 < len(roots) else root + max(0.4, root - left_neighbor)
    left = max(beta_workflow.ROOT_SCAN_START, root - 0.45 * (root - left_neighbor))
    right = root + 0.45 * (right_neighbor - root)
    return float(left), float(right)


def _local_root(case: OnlineCase, mode_index: int, counts: OperationCounts) -> tuple[float, bool, str]:
    left, right = local_interval(case.eb_roots, mode_index)
    before_det = counts.Timo_characteristic_determinant_evaluations
    before_brent = counts.Timo_brent_callback_evaluations
    counts.Timo_root_solver_calls += 1
    counts.Timo_roots_requested += 1
    root, _reported_evals, failure, notes = fixed_scan.local_timo_root(
        beta_deg=case.beta_deg,
        mu=case.mu,
        eta=case.eta,
        epsilon=case.epsilon_0,
        lambda_eb=case.eb_roots[mode_index - 1],
        left=left,
        right=right,
    )
    det_delta = counts.Timo_characteristic_determinant_evaluations - before_det
    brent_delta = counts.Timo_brent_callback_evaluations - before_brent
    grid_count = min(80, max(0, det_delta - brent_delta))
    counts.Timo_scan_grid_evaluations += grid_count
    counts.Timo_minimizer_callback_evaluations += max(0, det_delta - brent_delta - grid_count)
    if failure or not math.isfinite(float(root)):
        return float("nan"), False, f"mode {mode_index}: {notes or 'local root failure'}"
    counts.Timo_roots_returned_by_solver_calls += 1
    spec = beta_workflow.CaseSpec(mu=case.mu, eta=case.eta, epsilon=case.epsilon_0)
    counts.strict_root_quality_checks += 1
    counts.Timo_characteristic_matrix_evaluations += 1
    counts.Timo_SVD_calls += 1
    sv_min, _sv_second = beta_workflow.model_singular_values(spec, case.beta_deg, float(root), beta_workflow.MODEL_TIMO)
    if sv_min > beta_workflow.CONTINUATION_SVD_ACCEPT:
        return float(root), False, f"mode {mode_index}: local residual SVD {sv_min:.6e} exceeds existing acceptance"
    if not (left <= float(root) <= right):
        return float(root), False, f"mode {mode_index}: local root escaped online interval"
    return float(root), True, f"mode {mode_index}: {notes}; sv_min={sv_min:.6e}"


def suffix_has_cluster_risk(case: OnlineCase, n_safe: int) -> bool:
    """Flag suffix roots in an EB cluster using the existing continuation cluster tolerance."""

    for mode_index in range(n_safe + 1, K_MAX + 1):
        index = mode_index - 1
        neighbour_gaps = []
        if index > 0:
            neighbour_gaps.append(case.eb_roots[index] - case.eb_roots[index - 1])
        if index + 1 < K_MAX:
            neighbour_gaps.append(case.eb_roots[index + 1] - case.eb_roots[index])
        if any(gap <= beta_workflow.SVD_RECOVERY_CLUSTER_TOL for gap in neighbour_gaps):
            return True
    return False


def solve_hybrid_timo(case: OnlineCase, threshold: float, n_points: int, counts: OperationCounts) -> HybridResult:
    n_safe, pi_values = run_rule_s_selector(case, threshold, n_points, counts)
    counts.suffix_frequency_count = K_MAX - n_safe
    counts.Timoshenko_frequency_count = K_MAX - n_safe
    if n_safe == K_MAX:
        return HybridResult(n_safe, pi_values, (), (), False, ())
    local_roots: list[float] = []
    local_notes: list[str] = []
    reliable = not suffix_has_cluster_risk(case, n_safe)
    if not reliable:
        local_notes.append("online reliability trigger: EB-neighbour gap is within existing cluster tolerance")
    else:
        for mode_index in range(n_safe + 1, K_MAX + 1):
            root, passed, notes = _local_root(case, mode_index, counts)
            local_notes.append(notes)
            if not passed or (local_roots and root <= local_roots[-1] + 1.0e-8):
                reliable = False
                break
            local_roots.append(root)
    if reliable and len(local_roots) == K_MAX - n_safe:
        return HybridResult(
            n_safe,
            pi_values,
            tuple(local_roots),
            tuple("local_existing_EB_interval" for _ in local_roots),
            False,
            tuple(local_notes),
        )
    counts.full_spectrum_fallback_count += 1
    full_roots = solve_direct_timo(case, counts)
    suffix = full_roots[n_safe:]
    local_notes.append("online reliability trigger: full K=10 fallback; saved references were not consulted")
    return HybridResult(
        n_safe,
        pi_values,
        tuple(suffix),
        tuple("full_K10_fallback" for _ in suffix),
        True,
        tuple(local_notes),
    )


def reference_error(value: float, reference: float) -> tuple[float, float, float, bool]:
    absolute = abs(float(value) - float(reference))
    relative = absolute / max(abs(float(reference)), np.finfo(float).tiny)
    tolerance = REFERENCE_ABS_TOL + REFERENCE_REL_TOL * abs(float(reference))
    return absolute, relative, tolerance, bool(absolute <= tolerance)


def apply_postsolve_verification_fallback(
    case: OnlineCase,
    references: Sequence[float],
    hybrid: HybridResult,
    counts: OperationCounts,
) -> HybridResult:
    """Apply the specified verification-mismatch fallback after local solving, never during bracketing."""

    if hybrid.fallback_used or hybrid.N_S == K_MAX:
        return hybrid
    mismatches: list[int] = []
    for mode_index, value in enumerate(hybrid.suffix_roots, start=hybrid.N_S + 1):
        _absolute, _relative, _tolerance, passed = reference_error(value, references[mode_index - 1])
        if not passed:
            mismatches.append(mode_index)
    if not mismatches:
        return hybrid
    counts.full_spectrum_fallback_count += 1
    full_roots = solve_direct_timo(case, counts)
    notes = list(hybrid.local_notes)
    notes.append(
        "post-solve verification mismatch at modes "
        + ",".join(str(index) for index in mismatches)
        + "; independent global K=10 fallback (references not used for bracketing)"
    )
    suffix = full_roots[hybrid.N_S:]
    return HybridResult(
        hybrid.N_S,
        hybrid.pi_shear,
        tuple(suffix),
        tuple("full_K10_fallback_after_postsolve_verification" for _ in suffix),
        True,
        tuple(notes),
    )


def operation_row(case_id: str, workflow: str, counts: OperationCounts) -> dict[str, object]:
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "case_id": case_id,
        "workflow": workflow,
        "canonical_repeat": 0,
        **{key: value for key, value in asdict(counts).items() if key != "warnings"},
        "warnings": " | ".join(counts.warnings or []),
    }


def aggregate_operation_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    numeric_fields = [field for field in OUTPUT_FIELDS["aggregate_operation_counts.csv"] if field not in {
        "algorithm_version", "workflow", "case_count", "suffix_case_count"
    }]
    output: list[dict[str, object]] = []
    for workflow in ("direct_Timoshenko_K10", "hybrid_Rule_S"):
        subset = [row for row in rows if row["workflow"] == workflow]
        output.append({
            "algorithm_version": ALGORITHM_VERSION,
            "workflow": workflow,
            "case_count": len(subset),
            "suffix_case_count": sum(int(row["suffix_frequency_count"]) > 0 for row in subset),
            **{field: sum(int(row.get(field, 0)) for row in subset) for field in numeric_fields},
        })
    return output


def run_benchmark(args: Args) -> dict[str, list[dict[str, object]]]:
    cases, manifest, input_audit = load_cases(args.proposal_csv)
    direct_rows: list[dict[str, object]] = []
    hybrid_rows: list[dict[str, object]] = []
    verification_rows: list[dict[str, object]] = []
    operation_rows: list[dict[str, object]] = []
    runtime_rows: list[dict[str, object]] = []
    canonical_direct: dict[str, tuple[float, ...]] = {}
    canonical_hybrid: dict[str, HybridResult] = {}

    for repeat_index in range(args.runtime_repeats):
        for record in cases:
            permitted = online_case(record)
            direct_counts = OperationCounts(
                suffix_frequency_count=K_MAX - record.expected_N_S,
                Timoshenko_frequency_count=K_MAX,
            )
            start = time.perf_counter()
            with instrument_timo(direct_counts):
                direct_roots = solve_direct_timo(permitted, direct_counts)
            direct_elapsed = time.perf_counter() - start
            runtime_rows.append({
                "algorithm_version": ALGORITHM_VERSION,
                "case_id": record.case_id,
                "workflow": "direct_Timoshenko_K10",
                "repeat_index": repeat_index,
                "temperature": "cold" if repeat_index == 0 else "warm",
                "elapsed_seconds": direct_elapsed,
                "result_hash": stable_hash(direct_roots),
                "N_S": "",
                "suffix_frequency_count": K_MAX - record.expected_N_S,
                "full_spectrum_fallback_count": direct_counts.full_spectrum_fallback_count,
                "Timo_root_solver_calls": direct_counts.Timo_root_solver_calls,
                "Timo_roots_returned_by_solver_calls": direct_counts.Timo_roots_returned_by_solver_calls,
                "Timo_characteristic_determinant_evaluations": direct_counts.Timo_characteristic_determinant_evaluations,
                "Timo_SVD_calls": direct_counts.Timo_SVD_calls,
                "Timoshenko_frequency_count": direct_counts.Timoshenko_frequency_count,
                "EB_mode_reconstructions": direct_counts.EB_mode_reconstructions,
                "EB_SVD_6x6_calls": direct_counts.EB_SVD_6x6_calls,
                "EB_quadrature_calls": direct_counts.EB_quadrature_calls,
            })
            if repeat_index == 0:
                canonical_direct[record.case_id] = direct_roots
                operation_rows.append(operation_row(record.case_id, "direct_Timoshenko_K10", direct_counts))

            hybrid_counts = OperationCounts()
            start = time.perf_counter()
            with instrument_timo(hybrid_counts):
                hybrid = solve_hybrid_timo(permitted, record.frozen_T_s, args.n_shape_points, hybrid_counts)
                hybrid = apply_postsolve_verification_fallback(
                    permitted,
                    record.timo_references,
                    hybrid,
                    hybrid_counts,
                )
            hybrid_elapsed = time.perf_counter() - start
            if hybrid.N_S != record.expected_N_S:
                raise RuntimeError(
                    f"{record.case_id}: online Rule S produced N={hybrid.N_S}, expected frozen Phase-I N={record.expected_N_S}"
                )
            runtime_rows.append({
                "algorithm_version": ALGORITHM_VERSION,
                "case_id": record.case_id,
                "workflow": "hybrid_Rule_S",
                "repeat_index": repeat_index,
                "temperature": "cold" if repeat_index == 0 else "warm",
                "elapsed_seconds": hybrid_elapsed,
                "result_hash": stable_hash(hybrid.suffix_roots),
                "N_S": hybrid.N_S,
                "suffix_frequency_count": K_MAX - hybrid.N_S,
                "full_spectrum_fallback_count": hybrid_counts.full_spectrum_fallback_count,
                "Timo_root_solver_calls": hybrid_counts.Timo_root_solver_calls,
                "Timo_roots_returned_by_solver_calls": hybrid_counts.Timo_roots_returned_by_solver_calls,
                "Timo_characteristic_determinant_evaluations": hybrid_counts.Timo_characteristic_determinant_evaluations,
                "Timo_SVD_calls": hybrid_counts.Timo_SVD_calls,
                "Timoshenko_frequency_count": hybrid_counts.Timoshenko_frequency_count,
                "EB_mode_reconstructions": hybrid_counts.EB_mode_reconstructions,
                "EB_SVD_6x6_calls": hybrid_counts.EB_SVD_6x6_calls,
                "EB_quadrature_calls": hybrid_counts.EB_quadrature_calls,
            })
            if repeat_index == 0:
                canonical_hybrid[record.case_id] = hybrid
                operation_rows.append(operation_row(record.case_id, "hybrid_Rule_S", hybrid_counts))

    for record in cases:
        direct = canonical_direct[record.case_id]
        hybrid = canonical_hybrid[record.case_id]
        for mode_index, (root, reference) in enumerate(zip(direct, record.timo_references, strict=True), start=1):
            absolute, relative, _tolerance, passed = reference_error(root, reference)
            direct_rows.append({
                "algorithm_version": ALGORITHM_VERSION,
                "case_id": record.case_id,
                "mode_index": mode_index,
                "Lambda_Timo_direct": root,
                "Lambda_Timo_reference": reference,
                "absolute_error": absolute,
                "relative_error": relative,
                "verification_passed": passed,
                "root_source": "existing_fast_direct_Timoshenko_sign_scan",
                "warnings": next(
                    str(row["warnings"]) for row in operation_rows
                    if row["case_id"] == record.case_id and row["workflow"] == "direct_Timoshenko_K10"
                ),
            })
        suffix_by_mode = {
            mode_index: (root, source)
            for mode_index, (root, source) in enumerate(
                zip(hybrid.suffix_roots, hybrid.suffix_sources, strict=True), start=hybrid.N_S + 1
            )
        }
        for mode_index in range(1, K_MAX + 1):
            reference = record.timo_references[mode_index - 1]
            if mode_index <= hybrid.N_S:
                value = record.eb_roots[mode_index - 1]
                source = "accepted_EB_frequency"
            else:
                value, source = suffix_by_mode[mode_index]
            absolute, relative, tolerance, passed = reference_error(value, reference)
            if mode_index <= hybrid.N_S:
                delta_squared = abs(value * value - reference * reference) / (reference * reference)
                passed = bool(delta_squared <= 0.10)
            hybrid_rows.append({
                "algorithm_version": ALGORITHM_VERSION,
                "case_id": record.case_id,
                "mode_index": mode_index,
                "N_S": hybrid.N_S,
                "accepted_by_rule_S": mode_index <= hybrid.N_S,
                "frequency_source": source,
                "Lambda_hybrid": value,
                "Lambda_reference": reference,
                "absolute_error": absolute,
                "relative_error": relative,
                "verification_passed": passed,
                "fallback_used": hybrid.fallback_used,
                "local_attempt_notes": " | ".join(hybrid.local_notes),
            })
            if mode_index > hybrid.N_S:
                verification_rows.append({
                    "algorithm_version": ALGORITHM_VERSION,
                    "case_id": record.case_id,
                    "mode_index": mode_index,
                    "online_root_source": source,
                    "Lambda_online": value,
                    "Lambda_saved_reference": reference,
                    "absolute_error": absolute,
                    "relative_error": relative,
                    "verification_tolerance": tolerance,
                    "verification_passed": passed,
                    "verification_performed_after_online_solve": True,
                    "reference_used_for_bracketing": False,
                })

    aggregate_rows = aggregate_operation_rows(operation_rows)
    aggregates = {row["workflow"]: row for row in aggregate_rows}
    all_verified = all(bool(row["verification_passed"]) for row in direct_rows + hybrid_rows)
    suffix_cases = [record for record in cases if record.expected_N_S < K_MAX]
    fallback_cases = sum(canonical_hybrid[record.case_id].fallback_used for record in suffix_cases)
    systematic_fallback = bool(suffix_cases and fallback_cases == len(suffix_cases))
    direct_runtime = [float(row["elapsed_seconds"]) for row in runtime_rows if row["workflow"] == "direct_Timoshenko_K10"]
    hybrid_runtime = [float(row["elapsed_seconds"]) for row in runtime_rows if row["workflow"] == "hybrid_Rule_S"]
    direct_aggregate = aggregates["direct_Timoshenko_K10"]
    hybrid_aggregate = aggregates["hybrid_Rule_S"]
    if not all_verified:
        status = "rule_S_cost_benchmark_inconclusive_data_quality"
        basis = "at least one direct or hybrid post-solve reference verification failed"
    elif systematic_fallback:
        status = "rule_S_cost_not_beneficial"
        basis = "every nonempty suffix case triggered a full K=10 fallback after local or post-solve verification checks"
    elif (
        int(hybrid_aggregate["Timo_characteristic_determinant_evaluations"])
        >= int(direct_aggregate["Timo_characteristic_determinant_evaluations"])
        and int(hybrid_aggregate["EB_mode_reconstructions"]) > 0
    ):
        status = "rule_S_cost_not_beneficial"
        basis = "hybrid used no fewer Timoshenko determinant evaluations and added EB reconstruction work"
    elif (
        int(hybrid_aggregate["Timo_characteristic_determinant_evaluations"])
        < int(direct_aggregate["Timo_characteristic_determinant_evaluations"])
        and statistics.median(hybrid_runtime) < statistics.median(direct_runtime)
        and fallback_cases == 0
    ):
        status = "rule_S_cost_beneficial_on_benchmark"
        basis = "hybrid reduced dominant Timoshenko evaluations and median runtime without fallback"
    else:
        status = "rule_S_cost_inconclusive_without_common_operation_model"
        basis = "componentwise operation and runtime evidence do not establish dominance without arbitrary weights"
    summary_rows = [{
        "algorithm_version": ALGORITHM_VERSION,
        "status": status,
        "case_count": len(cases),
        "runtime_repeats": args.runtime_repeats,
        "all_reference_verifications_passed": all_verified,
        "suffix_case_count": len(suffix_cases),
        "suffix_cases_with_full_fallback": fallback_cases,
        "systematic_full_fallback": systematic_fallback,
        "direct_total_Timo_determinant_evaluations": direct_aggregate["Timo_characteristic_determinant_evaluations"],
        "hybrid_total_Timo_determinant_evaluations": hybrid_aggregate["Timo_characteristic_determinant_evaluations"],
        "direct_total_Timo_root_solver_calls": direct_aggregate["Timo_root_solver_calls"],
        "hybrid_total_Timo_root_solver_calls": hybrid_aggregate["Timo_root_solver_calls"],
        "hybrid_total_EB_mode_reconstructions": hybrid_aggregate["EB_mode_reconstructions"],
        "direct_median_runtime_seconds": statistics.median(direct_runtime),
        "hybrid_median_runtime_seconds": statistics.median(hybrid_runtime),
        "common_scalar_cost_model_used": False,
        "decision_basis": basis,
        "new_EB_root_calculations": 0,
        "new_Timo_root_solver_calls_executed": sum(int(row["Timo_root_solver_calls"]) for row in runtime_rows),
        "new_Timo_roots_returned_during_all_repeats": sum(
            int(row["Timo_roots_returned_by_solver_calls"]) for row in runtime_rows
        ),
        "new_Timo_determinant_evaluations_executed": sum(
            int(row["Timo_characteristic_determinant_evaluations"]) for row in runtime_rows
        ),
        "new_EB_mode_reconstructions_executed": sum(int(row["EB_mode_reconstructions"]) for row in runtime_rows),
        "new_EB_SVD_6x6_calls_executed": sum(int(row["EB_SVD_6x6_calls"]) for row in runtime_rows),
        "new_EB_quadrature_calls_executed": sum(int(row["EB_quadrature_calls"]) for row in runtime_rows),
        "cache_reused": False,
    }]
    return {
        "cost_benchmark_manifest.csv": manifest,
        "cost_benchmark_input_audit.csv": input_audit,
        "direct_timo_predictions.csv": direct_rows,
        "hybrid_rule_S_predictions.csv": hybrid_rows,
        "suffix_root_verification.csv": verification_rows,
        "per_case_operation_counts.csv": operation_rows,
        "aggregate_operation_counts.csv": aggregate_rows,
        "runtime_measurements.csv": runtime_rows,
        "cost_break_even_summary.csv": summary_rows,
    }


def save_cache(path: Path, outputs: Mapping[str, Sequence[Mapping[str, object]]], args: Args) -> None:
    phase1_dir = args.proposal_csv.parent
    payload = {
        "algorithm_version": ALGORITHM_VERSION,
        "input_sha256": {
            "proposal": sha256_file(args.proposal_csv),
            "prefix": sha256_file(phase1_dir / "prefix_predictor_metrics.csv"),
            "rule_S_predictions": sha256_file(phase1_dir / "rule_S_predictions.csv"),
        },
        "runtime_repeats": args.runtime_repeats,
        "n_shape_points": args.n_shape_points,
        "outputs": outputs,
    }
    path.write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def load_cache(path: Path, args: Args) -> dict[str, list[dict[str, object]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    phase1_dir = args.proposal_csv.parent
    expected_hashes = {
        "proposal": sha256_file(args.proposal_csv),
        "prefix": sha256_file(phase1_dir / "prefix_predictor_metrics.csv"),
        "rule_S_predictions": sha256_file(phase1_dir / "rule_S_predictions.csv"),
    }
    cached_hashes = payload.get("input_sha256")
    if cached_hashes is None:
        legacy_audit = payload["outputs"]["cost_benchmark_input_audit.csv"][0]
        cached_hashes = {
            "proposal": payload.get("proposal_sha256"),
            "prefix": legacy_audit.get("prefix_sha256"),
            "rule_S_predictions": legacy_audit.get("rule_S_predictions_sha256"),
        }
    expected = (ALGORITHM_VERSION, expected_hashes, args.runtime_repeats, args.n_shape_points)
    actual = (
        payload.get("algorithm_version"), cached_hashes, payload.get("runtime_repeats"),
        payload.get("n_shape_points"),
    )
    if actual != expected:
        raise ValueError("benchmark cache metadata does not match current inputs/options")
    outputs = {name: list(rows) for name, rows in payload["outputs"].items()}
    summary = outputs["cost_break_even_summary.csv"][0]
    summary["cache_reused"] = True
    return outputs


def normalize_output_accounting(outputs: dict[str, list[dict[str, object]]]) -> dict[str, list[dict[str, object]]]:
    """Upgrade cached bookkeeping without recomputing any spectrum."""

    n_s_by_case = {
        str(row["case_id"]): int(row["N_S"])
        for row in outputs["cost_benchmark_manifest.csv"]
    }
    for row in outputs["per_case_operation_counts.csv"]:
        n_s = n_s_by_case[str(row["case_id"])]
        row["suffix_frequency_count"] = K_MAX - n_s
        row["Timoshenko_frequency_count"] = K_MAX if row["workflow"] == "direct_Timoshenko_K10" else K_MAX - n_s
    outputs["aggregate_operation_counts.csv"] = aggregate_operation_rows(outputs["per_case_operation_counts.csv"])
    operation_by_key = {
        (str(row["case_id"]), str(row["workflow"])): row
        for row in outputs["per_case_operation_counts.csv"]
    }
    for row in outputs["runtime_measurements.csv"]:
        n_s = n_s_by_case[str(row["case_id"])]
        row["suffix_frequency_count"] = K_MAX - n_s
        row["Timoshenko_frequency_count"] = K_MAX if row["workflow"] == "direct_Timoshenko_K10" else K_MAX - n_s
        reconstructions = int(row.get("EB_mode_reconstructions", 0))
        canonical = operation_by_key[(str(row["case_id"]), str(row["workflow"]))]
        row.setdefault("Timo_SVD_calls", int(canonical.get("Timo_SVD_calls", 0)))
        row.setdefault("EB_SVD_6x6_calls", reconstructions)
        row.setdefault("EB_quadrature_calls", 6 * reconstructions)
    summary = outputs["cost_break_even_summary.csv"][0]
    runtime = outputs["runtime_measurements.csv"]
    summary["new_Timo_determinant_evaluations_executed"] = sum(
        int(row["Timo_characteristic_determinant_evaluations"]) for row in runtime
    )
    summary["new_EB_mode_reconstructions_executed"] = sum(int(row["EB_mode_reconstructions"]) for row in runtime)
    summary["new_EB_SVD_6x6_calls_executed"] = sum(int(row["EB_SVD_6x6_calls"]) for row in runtime)
    summary["new_EB_quadrature_calls_executed"] = sum(int(row["EB_quadrature_calls"]) for row in runtime)
    if bool(summary.get("systematic_full_fallback")):
        summary["decision_basis"] = (
            "every nonempty suffix case triggered a full K=10 fallback after local or post-solve verification checks"
        )
    return outputs


def create_plots(output_dir: Path) -> list[Path]:
    runtime = read_csv(output_dir / "runtime_measurements.csv")
    operations = read_csv(output_dir / "per_case_operation_counts.csv")
    labels = list(EXPECTED_CASE_IDS)
    x = np.arange(len(labels), dtype=float)
    width = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    for offset, workflow, color, label in (
        (-width / 2, "direct_Timoshenko_K10", "#4c78a8", "direct Timoshenko K=10"),
        (width / 2, "hybrid_Rule_S", "#f58518", "hybrid Rule S"),
    ):
        values = [statistics.median(
            float(row["elapsed_seconds"]) for row in runtime
            if row["case_id"] == case_id and row["workflow"] == workflow
        ) for case_id in labels]
        ax.bar(x + offset, values, width, color=color, label=label)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Median elapsed time, s")
    ax.set_title("Targeted cost benchmark (all recorded repeats)")
    ax.legend()
    runtime_path = output_dir / PLOT_NAMES[0]
    fig.savefig(runtime_path, dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.8), constrained_layout=True)
    for offset, workflow, color, label in (
        (-width / 2, "direct_Timoshenko_K10", "#4c78a8", "direct"),
        (width / 2, "hybrid_Rule_S", "#f58518", "hybrid"),
    ):
        values = [int(next(
            row["Timo_characteristic_determinant_evaluations"] for row in operations
            if row["case_id"] == case_id and row["workflow"] == workflow
        )) for case_id in labels]
        ax.bar(x + offset, values, width, color=color, label=label)
    fallback_cases = {
        row["case_id"] for row in operations
        if row["workflow"] == "hybrid_Rule_S" and int(row["full_spectrum_fallback_count"]) > 0
    }
    for index, case_id in enumerate(labels):
        if case_id in fallback_cases:
            ax.text(index + width / 2, 0, "fallback", rotation=90, va="bottom", ha="center", fontsize=8)
    ax.set_xticks(x, labels)
    ax.set_ylabel("Timoshenko determinant evaluations")
    ax.set_title("Canonical-repeat operation counts")
    ax.legend()
    operation_path = output_dir / PLOT_NAMES[1]
    fig.savefig(operation_path, dpi=180)
    plt.close(fig)
    return [runtime_path, operation_path]


def write_report(output_dir: Path) -> None:
    summary = read_csv(output_dir / "cost_break_even_summary.csv")[0]
    manifest = read_csv(output_dir / "cost_benchmark_manifest.csv")
    operations = read_csv(output_dir / "aggregate_operation_counts.csv")
    report = [
        "# Frozen Rule S cost break-even benchmark",
        "",
        f"Algorithm: `{ALGORITHM_VERSION}`.",
        "",
        "This targeted benchmark compares the existing fast direct Timoshenko K=10 workflow with the frozen Rule S "
        "selector followed by suffix-only local Timoshenko attempts and the documented online full-spectrum fallback. "
        "No EB root was solved. Saved Timoshenko roots were withheld from online bracketing and used only after each solve.",
        "",
        "## Locked inputs",
        "",
        f"The cases are `{', '.join(row['case_id'] for row in manifest)}` and the frozen threshold is "
        f"`T_s={manifest[0]['frozen_T_s']}`. No geometry, threshold, or rule definition was changed after inspecting outcomes.",
        "",
        "## Result",
        "",
        f"Primary cost status: `{summary['status']}`.",
        "",
        f"Decision basis: {summary['decision_basis']}.",
        "",
        f"All post-solve checks against saved references passed: `{summary['all_reference_verifications_passed']}`. "
        f"Full K=10 fallback occurred in `{summary['suffix_cases_with_full_fallback']}` of "
        f"`{summary['suffix_case_count']}` nonempty-suffix cases.",
        "",
        "## Componentwise operation accounting",
        "",
        "| workflow | Timo determinant evaluations | Timo solver calls | EB reconstructions | full fallbacks |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in operations:
        report.append(
            f"| {row['workflow']} | {row['Timo_characteristic_determinant_evaluations']} | "
            f"{row['Timo_root_solver_calls']} | {row['EB_mode_reconstructions']} | "
            f"{row['full_spectrum_fallback_count']} |"
        )
    report.extend([
        "",
        f"Median measured runtime over all case/repeat observations was `{summary['direct_median_runtime_seconds']}` s "
        f"for direct and `{summary['hybrid_median_runtime_seconds']}` s for hybrid. Runtime is supporting evidence only.",
        "",
        "No arbitrary scalar weighting of EB reconstructions, determinant evaluations, SVDs, quadratures, and runtime was used.",
        "",
        "## Scope",
        "",
        "This is a finite five-geometry operation-cost experiment, not a guarantee on the continuous parameter domain. "
        "It does not alter the separate finite-sample safety result for Rule S/Rule B.",
        "",
    ])
    (output_dir / REPORT_NAME).write_text("\n".join(report), encoding="utf-8")


def ensure_write_policy(args: Args) -> None:
    existing = [args.output_dir / name for name in (*OUTPUT_FIELDS.keys(), REPORT_NAME, CACHE_NAME) if (args.output_dir / name).exists()]
    if existing and not args.force:
        raise FileExistsError("outputs already exist; pass --force or use --plot-only")


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.plot_only:
        missing = [args.output_dir / name for name in ("runtime_measurements.csv", "per_case_operation_counts.csv") if not (args.output_dir / name).exists()]
        if missing:
            raise FileNotFoundError(f"plot-only inputs are missing: {missing}")
        create_plots(args.output_dir)
        print(f"algorithm {ALGORITHM_VERSION}")
        print("plot-only: numerical CSV unchanged; root calculations 0")
        return 0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ensure_write_policy(args)
    cache_path = args.output_dir / CACHE_NAME
    if args.reuse_cache and cache_path.exists():
        outputs = normalize_output_accounting(load_cache(cache_path, args))
        save_cache(cache_path, outputs, args)
    else:
        outputs = normalize_output_accounting(run_benchmark(args))
        save_cache(cache_path, outputs, args)
    for name, rows in outputs.items():
        write_csv(args.output_dir / name, rows, OUTPUT_FIELDS[name])
    create_plots(args.output_dir)
    write_report(args.output_dir)
    summary = outputs["cost_break_even_summary.csv"][0]
    print(f"algorithm {ALGORITHM_VERSION}")
    print(f"status {summary['status']}")
    print("new EB root calculations 0")
    calls_this_invocation = 0 if bool(summary.get("cache_reused")) else summary["new_Timo_root_solver_calls_executed"]
    print(f"new Timoshenko root solver calls {calls_this_invocation}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
