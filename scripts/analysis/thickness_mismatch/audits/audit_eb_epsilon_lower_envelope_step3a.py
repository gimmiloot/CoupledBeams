from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path
import sys
import time
from typing import Callable, Mapping, Sequence

import numpy as np


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[4]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import epsilon_lower_envelope_metrics as metrics  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402
from scripts.lib import straight_rod_factorized_spectrum as factorized  # noqa: E402


STEP3A_ALGORITHM_VERSION = "epsilon_lower_envelope_step3a_v1"
DEFAULT_MANIFEST = (
    REPO_ROOT
    / "scripts"
    / "analysis"
    / "thickness_mismatch"
    / "audits"
    / "data"
    / "eb_epsilon_lower_envelope_step3_cases.csv"
)
DEFAULT_BASELINE_DIR = REPO_ROOT / "results" / "eb_epsilon_baseline_thresholds"
DEFAULT_GATEWAY_DIR = REPO_ROOT / "results" / "eb_timo_branch_continuation_gateway"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "eb_epsilon_lower_envelope_step3a"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "eb_epsilon_lower_envelope_step3a"
THRESHOLD_FILE = "baseline_critical_prefix_thresholds.csv"
FACTOR_FILE = "baseline_factorized_spectrum_audit.csv"
BASELINE_MODE_FILE = "baseline_mode_metrics.csv"
BASELINE_REPORT = "eb_epsilon_baseline_thresholds_report.md"
GATEWAY_SUMMARY = "branch_k10_guard_summary.csv"
GATEWAY_REPORT = "eb_timo_branch_continuation_gateway_report.md"

OUTPUT_NAMES = (
    "step3a_manifest_resolved.csv",
    "step3a_case_spectrum_summary.csv",
    "step3a_mode_metrics.csv",
    "step3a_case_summary.csv",
    "step3a_prefix_screening.csv",
    "step3a_prefix_group_summary.csv",
    "step3a_baseline_control_audit.csv",
    "step3a_strict_verification_audit.csv",
    "step3a_counterexample_audit.csv",
    "step3a_exclusion_audit.csv",
    "step3a_operation_counts.csv",
    "step3b_followup_case_manifest.csv",
    "eb_epsilon_lower_envelope_step3a_report.md",
)
PLOT_NAMES = (
    "step3a_margin_by_prefix.png",
    "step3a_violation_by_case.png",
    "step3a_N_true_vs_baseline_certificate.png",
    "step3a_near_buffer_screen.png",
    "step3a_trigger_mode_map.png",
    "step3a_operation_costs.png",
)
SMOKE_CASE_IDS = (
    "S3_01",
    "S3_02",
    "S3_03",
    "S3_04",
    "S3_09",
    "S3_10",
    "S3_11",
    "S3_12",
    "S3_13",
    "S3_14",
    "S3_15",
    "S3_16",
)


@dataclass(frozen=True)
class ManifestCase:
    case_id: str
    prefix_group: str
    epsilon_source: str
    epsilon: float
    beta_deg: float
    mu: float
    eta: float
    case_group: str
    adversarial_rationale: str
    required_spectrum_method: str
    required_K10_guard: bool
    notes: str

    @property
    def prefixes(self) -> tuple[int, ...]:
        return metrics.EXPECTED_PREFIX_GROUPS[self.prefix_group]

    @property
    def geometry(self) -> complete.Geometry:
        return complete.Geometry(self.epsilon, self.beta_deg, self.mu, self.eta)

    @property
    def geometry_id(self) -> str:
        return (
            f"eps={self.epsilon:.16e}|beta={self.beta_deg:.16e}|"
            f"mu={self.mu:.16e}|eta={self.eta:.16e}"
        )


@dataclass(frozen=True)
class BaselineThreshold:
    prefix_n: int
    threshold_status: str
    epsilon_certified_n: float
    epsilon_safe_lower: float
    epsilon_unsafe_upper: float
    epsilon_star_estimate: float
    epsilon_near_n: float
    epsilon_buffer_n: float
    numerical_verification_status: str
    verification_bracket_agreement: bool
    verification_root_order_agreement: bool


@dataclass(frozen=True)
class BaselineData:
    thresholds: Mapping[int, BaselineThreshold]
    factorized_roots: Mapping[tuple[float, str], tuple[float, ...]]
    mode_deltas: Mapping[float, tuple[float, ...]]
    algorithm_version: str


@dataclass(frozen=True)
class Args:
    manifest: Path
    baseline_dir: Path
    gateway_dir: Path
    output_dir: Path
    cache_dir: Path
    verification_cache_dir: Path
    spectrum_method: str
    k_max: int
    n_spectrum_roots: int
    n_candidate_roots: int
    verification_candidate_roots: int
    threshold: float
    strict_verification_margin: float
    violation_tolerance: float
    baseline_provenance_tolerance: float
    initial_beta_step: float
    min_beta_step: float
    max_beta_step: float
    guard_scan_step: float
    seed_half_width: float
    sigma_accept: float
    sigma_ratio_accept: float
    mac_accept: float
    subspace_mac_accept: float
    cluster_gap_absolute: float
    cluster_gap_relative: float
    run_global_guard: bool
    allow_strict_fallback: bool
    verification_initial_beta_step: float
    verification_min_beta_step: float
    verification_max_beta_step: float
    verification_guard_scan_step: float
    reuse_cache: bool
    force_recompute: bool
    force: bool
    plot_only: bool
    smoke: bool
    write_step3b_followup_manifest: bool


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_bool(value: object) -> bool:
    token = str(value).strip().lower()
    if token in {"true", "1", "yes"}:
        return True
    if token in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def finite_or_nan(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


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
    if isinstance(value, (tuple, list, dict)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(fields or ())
    if not fieldnames:
        for row in rows:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
    if not fieldnames:
        fieldnames = ["status"]
        rows = ({"status": "no_rows"},)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fieldnames})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Run the fixed 28-case EB epsilon lower-envelope Step-3A screen.",
    )
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--baseline-dir", type=Path, default=DEFAULT_BASELINE_DIR)
    parser.add_argument("--gateway-dir", type=Path, default=DEFAULT_GATEWAY_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--verification-cache-dir", type=Path, default=None)
    parser.add_argument(
        "--spectrum-method",
        choices=(branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,),
        default=branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
    )
    parser.add_argument("--k-max", type=int, default=10)
    parser.add_argument("--n-spectrum-roots", type=int, default=12)
    parser.add_argument("--n-candidate-roots", type=int, default=20)
    parser.add_argument("--verification-candidate-roots", type=int, default=24)
    parser.add_argument("--threshold", type=float, default=0.10)
    parser.add_argument("--strict-verification-margin", type=float, default=0.005)
    parser.add_argument("--violation-tolerance", type=float, default=1.0e-5)
    parser.add_argument("--baseline-provenance-tolerance", type=float, default=1.0e-14)
    parser.add_argument("--initial-beta-step", type=float, default=0.5)
    parser.add_argument("--min-beta-step", type=float, default=0.005)
    parser.add_argument("--max-beta-step", type=float, default=5.0)
    parser.add_argument("--guard-scan-step", type=float, default=0.05)
    parser.add_argument("--seed-half-width", type=float, default=0.055)
    parser.add_argument("--sigma-accept", type=float, default=5.0e-6)
    parser.add_argument("--sigma-ratio-accept", type=float, default=5.0e-4)
    parser.add_argument("--mac-accept", type=float, default=0.25)
    parser.add_argument("--subspace-mac-accept", type=float, default=0.5)
    parser.add_argument("--cluster-gap-absolute", type=float, default=0.04)
    parser.add_argument("--cluster-gap-relative", type=float, default=3.0e-3)
    parser.add_argument("--global-guard", dest="run_global_guard", action="store_true", default=True)
    parser.add_argument("--no-global-guard", dest="run_global_guard", action="store_false")
    parser.add_argument("--strict-fallback", dest="allow_strict_fallback", action="store_true", default=True)
    parser.add_argument("--no-strict-fallback", dest="allow_strict_fallback", action="store_false")
    parser.add_argument("--verification-initial-beta-step", type=float, default=0.5)
    parser.add_argument("--verification-min-beta-step", type=float, default=0.005)
    parser.add_argument("--verification-max-beta-step", type=float, default=4.0)
    parser.add_argument("--verification-guard-scan-step", type=float, default=0.025)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--write-step3b-followup-manifest", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    output_dir = repo_path(Path(ns.output_dir))
    if bool(ns.smoke):
        output_dir = SMOKE_OUTPUT_DIR
    cache_dir = repo_path(Path(ns.cache_dir)) if ns.cache_dir is not None else output_dir / "cache"
    verification_cache_dir = (
        repo_path(Path(ns.verification_cache_dir))
        if ns.verification_cache_dir is not None
        else output_dir / "cache_verification"
    )
    if bool(ns.smoke):
        cache_dir = output_dir / "cache"
        verification_cache_dir = output_dir / "cache_verification"
    args = Args(
        manifest=repo_path(Path(ns.manifest)),
        baseline_dir=repo_path(Path(ns.baseline_dir)),
        gateway_dir=repo_path(Path(ns.gateway_dir)),
        output_dir=output_dir,
        cache_dir=cache_dir,
        verification_cache_dir=verification_cache_dir,
        spectrum_method=str(ns.spectrum_method),
        k_max=int(ns.k_max),
        n_spectrum_roots=int(ns.n_spectrum_roots),
        n_candidate_roots=int(ns.n_candidate_roots),
        verification_candidate_roots=int(ns.verification_candidate_roots),
        threshold=float(ns.threshold),
        strict_verification_margin=float(ns.strict_verification_margin),
        violation_tolerance=float(ns.violation_tolerance),
        baseline_provenance_tolerance=float(ns.baseline_provenance_tolerance),
        initial_beta_step=float(ns.initial_beta_step),
        min_beta_step=float(ns.min_beta_step),
        max_beta_step=float(ns.max_beta_step),
        guard_scan_step=float(ns.guard_scan_step),
        seed_half_width=float(ns.seed_half_width),
        sigma_accept=float(ns.sigma_accept),
        sigma_ratio_accept=float(ns.sigma_ratio_accept),
        mac_accept=float(ns.mac_accept),
        subspace_mac_accept=float(ns.subspace_mac_accept),
        cluster_gap_absolute=float(ns.cluster_gap_absolute),
        cluster_gap_relative=float(ns.cluster_gap_relative),
        run_global_guard=bool(ns.run_global_guard),
        allow_strict_fallback=bool(ns.allow_strict_fallback),
        verification_initial_beta_step=float(ns.verification_initial_beta_step),
        verification_min_beta_step=float(ns.verification_min_beta_step),
        verification_max_beta_step=float(ns.verification_max_beta_step),
        verification_guard_scan_step=float(ns.verification_guard_scan_step),
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        force=bool(ns.force),
        plot_only=bool(ns.plot_only),
        smoke=bool(ns.smoke),
        write_step3b_followup_manifest=bool(ns.write_step3b_followup_manifest),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if args.k_max != 10 or args.n_spectrum_roots < args.k_max + 2:
        raise ValueError("Step 3A requires --k-max 10 and --n-spectrum-roots >= K+2")
    if args.n_candidate_roots < args.n_spectrum_roots:
        raise ValueError("--n-candidate-roots must include all requested roots")
    if args.verification_candidate_roots <= args.n_candidate_roots:
        raise ValueError("--verification-candidate-roots must exceed --n-candidate-roots")
    if not 0.0 < args.threshold < 1.0:
        raise ValueError("--threshold must lie in (0, 1)")
    positives = (
        args.strict_verification_margin,
        args.violation_tolerance,
        args.baseline_provenance_tolerance,
        args.initial_beta_step,
        args.min_beta_step,
        args.max_beta_step,
        args.guard_scan_step,
        args.seed_half_width,
        args.sigma_accept,
        args.sigma_ratio_accept,
        args.mac_accept,
        args.subspace_mac_accept,
        args.cluster_gap_absolute,
        args.cluster_gap_relative,
        args.verification_initial_beta_step,
        args.verification_min_beta_step,
        args.verification_max_beta_step,
        args.verification_guard_scan_step,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in positives):
        raise ValueError("margins, tolerances, and continuation settings must be finite and positive")
    if not args.min_beta_step <= args.initial_beta_step <= args.max_beta_step:
        raise ValueError("primary beta steps must satisfy min <= initial <= max")
    if not args.verification_min_beta_step <= args.verification_initial_beta_step <= args.verification_max_beta_step:
        raise ValueError("verification beta steps must satisfy min <= initial <= max")
    if args.output_dir.resolve() in {args.baseline_dir.resolve(), args.gateway_dir.resolve()}:
        raise ValueError("Step-3A output directory must differ from baseline and gateway directories")
    if args.cache_dir.resolve() == args.verification_cache_dir.resolve():
        raise ValueError("primary and verification cache directories must be distinct")
    if not args.manifest.exists():
        raise FileNotFoundError(f"manifest does not exist: {args.manifest}")
    if not (args.baseline_dir / THRESHOLD_FILE).exists():
        raise FileNotFoundError(f"baseline threshold CSV does not exist: {args.baseline_dir / THRESHOLD_FILE}")
    if not (args.gateway_dir / GATEWAY_SUMMARY).exists() or not (args.gateway_dir / GATEWAY_REPORT).exists():
        raise FileNotFoundError("verified branch-gateway summary/report is missing")


def load_manifest(path: Path) -> list[ManifestCase]:
    rows = read_csv(path)
    required = {
        "case_id",
        "prefix_group",
        "epsilon_source",
        "epsilon",
        "beta_deg",
        "mu",
        "eta",
        "case_group",
        "adversarial_rationale",
        "required_spectrum_method",
        "required_K10_guard",
        "notes",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Step-3A manifest must contain columns: {sorted(required)}")
    cases: list[ManifestCase] = []
    for row in rows:
        cases.append(
            ManifestCase(
                case_id=str(row["case_id"]).strip(),
                prefix_group=str(row["prefix_group"]).strip(),
                epsilon_source=str(row["epsilon_source"]).strip(),
                epsilon=float(row["epsilon"]),
                beta_deg=float(row["beta_deg"]),
                mu=float(row["mu"]),
                eta=float(row["eta"]),
                case_group=str(row["case_group"]).strip(),
                adversarial_rationale=str(row["adversarial_rationale"]).strip(),
                required_spectrum_method=str(row["required_spectrum_method"]).strip(),
                required_K10_guard=parse_bool(row["required_K10_guard"]),
                notes=str(row.get("notes", "")).strip(),
            )
        )
    return cases


def validate_manifest_contract(cases: Sequence[ManifestCase]) -> None:
    if len(cases) != 28:
        raise ValueError(f"Step-3A manifest must contain exactly 28 rows, found {len(cases)}")
    case_ids = [case.case_id for case in cases]
    if len(case_ids) != len(set(case_ids)):
        raise ValueError("Step-3A manifest case IDs must be unique")
    geometry_ids = [case.geometry_id for case in cases]
    if len(geometry_ids) != len(set(geometry_ids)):
        raise ValueError("Step-3A manifest geometry keys must be unique")
    if {case.prefix_group for case in cases} != set(metrics.EXPECTED_PREFIX_GROUPS):
        raise ValueError("Step-3A manifest prefix groups differ from the explicit seven-group mapping")
    for case in cases:
        case.geometry.validate()
        if case.required_spectrum_method != branch.BRANCH_CONTINUATION_ALGORITHM_VERSION:
            raise ValueError(f"{case.case_id} does not require branch_informed_continuation_v1")
        if not case.required_K10_guard:
            raise ValueError(f"{case.case_id} does not require the K10 guard")
        if case.epsilon_source not in {"epsilon_near_n", "epsilon_buffer_n"}:
            raise ValueError(f"invalid epsilon source in {case.case_id}")
        if case.case_group not in {"baseline_control", "adversarial"}:
            raise ValueError(f"invalid case group in {case.case_id}")
        if case.case_group == "baseline_control" and any(
            abs(value) > 1.0e-14 for value in (case.beta_deg, case.mu, case.eta)
        ):
            raise ValueError(f"baseline control {case.case_id} must have beta=mu=eta=0")
    for group in metrics.EXPECTED_PREFIX_GROUPS:
        observed = Counter(
            (case.epsilon_source, case.case_group) for case in cases if case.prefix_group == group
        )
        expected = Counter(
            {
                ("epsilon_near_n", "baseline_control"): 1,
                ("epsilon_near_n", "adversarial"): 1,
                ("epsilon_buffer_n", "baseline_control"): 1,
                ("epsilon_buffer_n", "adversarial"): 1,
            }
        )
        if observed != expected:
            raise ValueError(f"prefix group {group} lacks the required near/buffer control/adversarial design")


def load_thresholds(baseline_dir: Path) -> dict[int, BaselineThreshold]:
    rows = read_csv(baseline_dir / THRESHOLD_FILE)
    thresholds: dict[int, BaselineThreshold] = {}
    for row in rows:
        prefix_n = int(row["prefix_n"])
        thresholds[prefix_n] = BaselineThreshold(
            prefix_n=prefix_n,
            threshold_status=str(row["threshold_status"]),
            epsilon_certified_n=float(row["epsilon_certified_n"]),
            epsilon_safe_lower=float(row["epsilon_safe_lower"]),
            epsilon_unsafe_upper=finite_or_nan(row.get("epsilon_unsafe_upper")),
            epsilon_star_estimate=finite_or_nan(row.get("epsilon_star_estimate")),
            epsilon_near_n=finite_or_nan(row.get("epsilon_near_n")),
            epsilon_buffer_n=finite_or_nan(row.get("epsilon_buffer_n")),
            numerical_verification_status=str(row.get("numerical_verification_status", "")),
            verification_bracket_agreement=parse_bool(row["verification_bracket_agreement"])
            if row.get("verification_bracket_agreement")
            else prefix_n == 1,
            verification_root_order_agreement=parse_bool(row["verification_root_order_agreement"])
            if row.get("verification_root_order_agreement")
            else prefix_n == 1,
        )
    if set(thresholds) != set(range(1, 11)):
        raise ValueError("corrected baseline threshold table must contain prefixes 1..10")
    for prefix_n, item in thresholds.items():
        if prefix_n == 1:
            continue
        if item.threshold_status != "resolved" or item.numerical_verification_status != "pass":
            raise ValueError(f"baseline prefix {prefix_n} is not a corrected verified threshold")
        if not item.verification_bracket_agreement or not item.verification_root_order_agreement:
            raise ValueError(f"baseline prefix {prefix_n} failed corrected verification gates")
    return thresholds


def validate_manifest_provenance(
    cases: Sequence[ManifestCase],
    thresholds: Mapping[int, BaselineThreshold],
    *,
    tolerance: float,
    baseline_dir: Path = DEFAULT_BASELINE_DIR,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in cases:
        differences: list[float] = []
        for prefix_n in case.prefixes:
            baseline = thresholds[prefix_n]
            expected = getattr(baseline, case.epsilon_source)
            difference = abs(case.epsilon - expected)
            differences.append(difference)
            if not math.isfinite(expected) or difference > tolerance:
                raise ValueError(
                    "manifest_baseline_provenance_failure: "
                    f"{case.case_id} {case.epsilon_source} differs from corrected prefix {prefix_n} by {difference}"
                )
        rows.append(
            {
                **asdict(case),
                "parsed_prefix_tuple": case.prefixes,
                "corrected_baseline_threshold_source": display_path(baseline_dir / THRESHOLD_FILE),
                "epsilon_provenance_difference": max(differences, default=0.0),
                "manifest_validation_status": "pass",
                "geometry_id": case.geometry_id,
                "step3a_algorithm_version": STEP3A_ALGORITHM_VERSION,
            }
        )
    return rows


def _matching_epsilon(value: float, targets: Sequence[float], tolerance: float) -> float | None:
    matches = [target for target in targets if abs(float(value) - float(target)) <= tolerance]
    return min(matches, key=lambda target: abs(float(value) - float(target))) if matches else None


def load_baseline_data(
    baseline_dir: Path,
    cases: Sequence[ManifestCase],
    thresholds: Mapping[int, BaselineThreshold],
    *,
    tolerance: float,
) -> BaselineData:
    report_path = baseline_dir / BASELINE_REPORT
    factor_path = baseline_dir / FACTOR_FILE
    mode_path = baseline_dir / BASELINE_MODE_FILE
    for path in (report_path, factor_path, mode_path):
        if not path.exists():
            raise FileNotFoundError(f"corrected baseline provenance file is missing: {path}")
    if factorized.ALGORITHM_VERSION not in report_path.read_text(encoding="utf-8"):
        raise ValueError("baseline report does not identify factorized_straight_spectrum_v2")
    target_eps = sorted({case.epsilon for case in cases if case.case_group == "baseline_control"})
    roots_lists: dict[tuple[float, str], dict[int, float]] = defaultdict(dict)
    with factor_path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            epsilon = _matching_epsilon(float(row["epsilon"]), target_eps, tolerance)
            if epsilon is None or abs(float(row["mu"])) > tolerance:
                continue
            if str(row.get("algorithm_version")) != factorized.ALGORITHM_VERSION:
                raise ValueError("legacy factorized baseline row encountered for a Step-3A control")
            if str(row.get("quality_status")) != "pass":
                raise ValueError("corrected factorized baseline root failed its quality gate")
            index = int(row["sorted_index"])
            if index <= 11:
                roots_lists[(epsilon, str(row["model"]))][index] = float(row["Lambda_factorized"])
    factor_roots: dict[tuple[float, str], tuple[float, ...]] = {}
    for epsilon in target_eps:
        for model in complete.SUPPORTED_MODELS:
            indexed = roots_lists[(epsilon, model)]
            if set(indexed) == set(range(1, 12)):
                factor_roots[(epsilon, model)] = tuple(indexed[index] for index in range(1, 12))
            else:
                # Near/buffer points are full-precision derivatives of the
                # corrected thresholds and are not necessarily persisted as
                # baseline scan rows.  Evaluate the same corrected factorized
                # oracle directly; never substitute a legacy general scan.
                geometry = complete.Geometry(epsilon, 0.0, 0.0, 0.0)
                factor_roots[(epsilon, model)] = complete.straight_oracle_values(model, geometry, 11)
    persisted_delta_lists: dict[float, dict[int, float]] = defaultdict(dict)
    with mode_path.open("r", newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            epsilon = _matching_epsilon(float(row["epsilon"]), target_eps, tolerance)
            if epsilon is None:
                continue
            index = int(row["sorted_index"])
            if index <= 10:
                persisted_delta_lists[epsilon][index] = float(row["delta_f"])
    mode_deltas: dict[float, tuple[float, ...]] = {}
    for epsilon in target_eps:
        eb_roots = factor_roots[(epsilon, complete.MODEL_EB)]
        timo_roots = factor_roots[(epsilon, complete.MODEL_TIMO)]
        direct = tuple(
            metrics.squared_frequency_delta(eb_roots[index], timo_roots[index])
            for index in range(10)
        )
        indexed = persisted_delta_lists[epsilon]
        if set(indexed) == set(range(1, 11)):
            persisted = tuple(indexed[index] for index in range(1, 11))
            if max(abs(left - right) for left, right in zip(direct, persisted)) > 1.0e-10:
                raise ValueError(f"direct corrected oracle disagrees with persisted baseline modes at epsilon={epsilon}")
        mode_deltas[epsilon] = direct
    return BaselineData(thresholds, factor_roots, mode_deltas, factorized.ALGORITHM_VERSION)


def validate_gateway_provenance(gateway_dir: Path) -> None:
    rows = read_csv(gateway_dir / GATEWAY_SUMMARY)
    if not rows or any(
        row.get("algorithm_version") != branch.BRANCH_CONTINUATION_ALGORITHM_VERSION for row in rows
    ):
        raise ValueError("gateway summary is not uniformly branch_informed_continuation_v1")
    report = (gateway_dir / GATEWAY_REPORT).read_text(encoding="utf-8")
    if "ready_for_targeted_step3" not in report:
        raise ValueError("branch gateway report is not marked READY")


def primary_settings(args: Args) -> branch.ContinuationSettings:
    return branch.ContinuationSettings(
        requested_roots=args.n_spectrum_roots,
        candidate_roots=args.n_candidate_roots,
        verification_candidate_roots=args.verification_candidate_roots,
        beta_initial_step_deg=args.initial_beta_step,
        beta_min_step_deg=args.min_beta_step,
        beta_max_step_deg=args.max_beta_step,
        guard_scan_step=args.guard_scan_step,
        seed_half_width=args.seed_half_width,
        sigma_accept=args.sigma_accept,
        sigma_ratio_accept=args.sigma_ratio_accept,
        mac_accept=args.mac_accept,
        subspace_mac_accept=args.subspace_mac_accept,
        cluster_gap_absolute=args.cluster_gap_absolute,
        cluster_gap_relative=args.cluster_gap_relative,
        run_global_guard=args.run_global_guard,
        allow_strict_fallback=args.allow_strict_fallback,
    )


def verification_settings(args: Args) -> branch.ContinuationSettings:
    return replace(
        primary_settings(args),
        beta_initial_step_deg=args.verification_initial_beta_step,
        beta_min_step_deg=args.verification_min_beta_step,
        beta_max_step_deg=args.verification_max_beta_step,
        guard_scan_step=args.verification_guard_scan_step,
        allow_strict_fallback=True,
    )


def resolve_case_models(
    case: ManifestCase,
    settings: branch.ContinuationSettings,
    cache: branch.BranchContinuationCache,
) -> dict[str, branch.BranchContinuationResult]:
    """Resolve roots from geometry only; screening metrics never enter this boundary."""

    geometry = case.geometry
    return {model: cache.resolve(model, geometry, settings) for model in complete.SUPPORTED_MODELS}


def normalized_gaps(values: Sequence[float]) -> tuple[float, ...]:
    gaps: list[float] = []
    for left, right in zip(values, values[1:]):
        denominator = float(left) + float(right)
        gaps.append(2.0 * abs(float(right) - float(left)) / denominator if denominator > 0.0 else float("nan"))
    return tuple(gaps)


def sided_min_gap(values: Sequence[float], index_zero_based: int) -> float:
    gaps = normalized_gaps(values)
    candidates: list[float] = []
    if index_zero_based > 0 and index_zero_based - 1 < len(gaps):
        candidates.append(gaps[index_zero_based - 1])
    if index_zero_based < len(gaps):
        candidates.append(gaps[index_zero_based])
    finite = [value for value in candidates if math.isfinite(value)]
    return min(finite) if finite else float("nan")


def ordered_branches(result: branch.BranchContinuationResult) -> tuple[branch.ContinuedBranch, ...]:
    return tuple(sorted(result.branches, key=lambda item: item.Lambda))


def cluster_size(item: branch.ContinuedBranch, ordered: Sequence[branch.ContinuedBranch]) -> int:
    if not item.cluster_id:
        return 1
    return sum(other.cluster_id == item.cluster_id for other in ordered)


def branch_reordered(result: branch.BranchContinuationResult, item: branch.ContinuedBranch, sorted_index: int) -> bool:
    parent_ids = [
        parent.branch_id
        for parent in sorted(result.parent_branches, key=lambda parent: (parent.Lambda, parent.family, parent.family_index))
    ]
    try:
        return parent_ids.index(item.branch_id) + 1 != int(sorted_index)
    except ValueError:
        return True


def operation_mapping(result: branch.BranchContinuationResult) -> dict[str, int]:
    operations = result.operations
    return {name: int(getattr(operations, name)) for name in branch.BranchOperationCounts.__dataclass_fields__}


def case_quality_triggers(
    case: ManifestCase,
    spectra: Mapping[str, branch.BranchContinuationResult],
    prefix_metrics: Sequence[metrics.PrefixMetric],
    settings: branch.ContinuationSettings,
) -> tuple[str, ...]:
    reasons: set[str] = set()
    target_triggers = {
        index
        for prefix_n in case.prefixes
        for index in prefix_metrics[prefix_n - 1].triggering_indices
    }
    for model, result in spectra.items():
        ordered = ordered_branches(result)
        values = result.values
        if result.operations.strict_fallback_runs > 0:
            reasons.add(f"{model}:strict_fallback_used")
        if result.exclusion_reason or result.guard.unresolved_intervals:
            reasons.add(f"{model}:continuation_warning")
        for trigger_index in target_triggers:
            if trigger_index <= len(ordered):
                item = ordered[trigger_index - 1]
                if cluster_size(item, ordered) > 1:
                    reasons.add(f"{model}:triggering_mode_close_cluster")
                if branch_reordered(result, item, trigger_index):
                    reasons.add(f"{model}:root_order_event")
                if sided_min_gap(values[:11], trigger_index - 1) <= settings.cluster_gap_relative:
                    reasons.add(f"{model}:unusually_small_root_gap")
        for prefix_n in case.prefixes:
            if prefix_n < len(ordered):
                left = ordered[prefix_n - 1]
                right = ordered[prefix_n]
                if left.cluster_id and left.cluster_id == right.cluster_id:
                    reasons.add(f"{model}:cluster_crosses_prefix_boundary")
        if not result.full12_resolved:
            root12_close = len(values) < 12
            if len(values) >= 12:
                root12_close = (
                    abs(values[11] - values[10]) <= settings.cluster_gap_absolute
                    or normalized_gaps(values[10:12])[0] <= settings.cluster_gap_relative
                )
            if root12_close:
                reasons.add(f"{model}:root11_depends_on_unresolved_full12_feature")
    return tuple(sorted(reasons))


def build_spectrum_summary_row(
    case: ManifestCase,
    result: branch.BranchContinuationResult,
    *,
    run_scope: str,
) -> dict[str, object]:
    ordered = ordered_branches(result)
    values = result.values
    cluster_count = len({item.cluster_id for item in ordered[:11] if item.cluster_id})
    min_gap = min(normalized_gaps(values[:11]), default=float("nan"))
    return {
        "case_id": case.case_id,
        "epsilon": case.epsilon,
        "beta_deg": case.beta_deg,
        "mu": case.mu,
        "eta": case.eta,
        "model": result.model,
        "run_scope": run_scope,
        "algorithm_version": result.algorithm_version,
        "root_count": len(values),
        "roots_1_11_status": "resolved" if result.k10_guard_resolved and len(values) >= 11 else "unresolved",
        "root12_status": "resolved" if result.full12_resolved else "optional_not_resolved",
        "K10_guard_resolved": result.k10_guard_resolved,
        "full12_resolved": result.full12_resolved,
        "strict_fallback_used": result.operations.strict_fallback_runs > 0,
        "cluster_count_first11": cluster_count,
        "minimum_normalized_gap_first11": min_gap,
        "root11": values[10] if len(values) >= 11 else float("nan"),
        "root12": values[11] if len(values) >= 12 else float("nan"),
        "guard_status": result.guard.status,
        "guard_passed": result.guard.passed,
        "unresolved_intervals_below_guard": result.guard.unresolved_intervals,
        "oracle_agreement": result.oracle_agreement,
        "force_verification_agreement": result.force_verification_agreement,
        "cache_status": result.cache_status,
        "warnings": result.exclusion_reason,
        **operation_mapping(result),
    }


def build_primary_products(
    cases: Sequence[ManifestCase],
    spectra_by_case: Mapping[str, Mapping[str, branch.BranchContinuationResult]],
    baseline: BaselineData,
    args: Args,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, dict[str, object]],
    list[dict[str, object]],
]:
    spectrum_rows: list[dict[str, object]] = []
    mode_rows: list[dict[str, object]] = []
    prefix_rows: list[dict[str, object]] = []
    operation_rows: list[dict[str, object]] = []
    calculations: dict[str, dict[str, object]] = {}
    safe_lower = {prefix_n: item.epsilon_certified_n for prefix_n, item in baseline.thresholds.items()}

    pair_lookup: dict[tuple[str, str, str], ManifestCase] = {}
    for case in cases:
        pair_lookup[(case.prefix_group, case.case_group, case.epsilon_source)] = case

    for case in cases:
        spectra = spectra_by_case[case.case_id]
        for result in spectra.values():
            spectrum_rows.append(build_spectrum_summary_row(case, result, run_scope="primary"))
            operation_rows.append(
                {
                    "case_id": case.case_id,
                    "model": result.model,
                    "cost_scope": "primary_screening",
                    **operation_mapping(result),
                }
            )
        eb = spectra[complete.MODEL_EB]
        timo = spectra[complete.MODEL_TIMO]
        included = (
            eb.k10_guard_resolved
            and timo.k10_guard_resolved
            and len(eb.values) >= 11
            and len(timo.values) >= 11
        )
        if included:
            deltas = tuple(
                metrics.squared_frequency_delta(eb.values[index], timo.values[index])
                for index in range(args.k_max)
            )
            prefix_metrics = metrics.running_prefix_metrics(deltas, threshold=args.threshold)
            n_true = metrics.true_safe_prefix(deltas, threshold=args.threshold)
            late = metrics.late_pass_indices(deltas, threshold=args.threshold)
        else:
            deltas = tuple(float("nan") for _ in range(args.k_max))
            prefix_metrics = metrics.running_prefix_metrics(deltas, threshold=args.threshold)
            n_true = 0
            late = ()
        n_certified = metrics.certified_prefix(case.epsilon, safe_lower, k_max=args.k_max)
        opposite = "epsilon_buffer_n" if case.epsilon_source == "epsilon_near_n" else "epsilon_near_n"
        paired = pair_lookup[(case.prefix_group, case.case_group, opposite)]
        same_pair = metrics.same_geometry(asdict(case), asdict(paired))
        quality_triggers = case_quality_triggers(case, spectra, prefix_metrics, primary_settings(args)) if included else ()
        calculations[case.case_id] = {
            "case": case,
            "deltas": deltas,
            "prefix_metrics": prefix_metrics,
            "N_true": n_true,
            "N_certified_0": n_certified,
            "late_pass_indices": late,
            "included": included,
            "quality_triggers": quality_triggers,
            "same_geometry_near_buffer_pair": same_pair,
        }

        eb_ordered = ordered_branches(eb)
        timo_ordered = ordered_branches(timo)
        for index in range(1, args.k_max + 1):
            eb_item = eb_ordered[index - 1] if index <= len(eb_ordered) else None
            timo_item = timo_ordered[index - 1] if index <= len(timo_ordered) else None
            delta = deltas[index - 1]
            mode_rows.append(
                {
                    "case_id": case.case_id,
                    "prefix_group": case.prefix_group,
                    "epsilon_source": case.epsilon_source,
                    "case_group": case.case_group,
                    "epsilon": case.epsilon,
                    "beta_deg": case.beta_deg,
                    "mu": case.mu,
                    "eta": case.eta,
                    "sorted_index": index,
                    "Lambda_EB_primary": eb.values[index - 1] if index <= len(eb.values) else float("nan"),
                    "Lambda_Timo_primary": timo.values[index - 1] if index <= len(timo.values) else float("nan"),
                    "delta_f_primary": delta,
                    "pass_10_primary": math.isfinite(delta) and delta <= args.threshold,
                    "EB_branch_id": eb_item.branch_id if eb_item else "",
                    "Timo_branch_id": timo_item.branch_id if timo_item else "",
                    "EB_parent_family": eb_item.parent_family if eb_item else "",
                    "Timo_parent_family": timo_item.parent_family if timo_item else "",
                    "EB_cluster_id": eb_item.cluster_id if eb_item else "",
                    "Timo_cluster_id": timo_item.cluster_id if timo_item else "",
                    "EB_cluster_size": cluster_size(eb_item, eb_ordered) if eb_item else 0,
                    "Timo_cluster_size": cluster_size(timo_item, timo_ordered) if timo_item else 0,
                    "EB_adjacent_gap": sided_min_gap(eb.values[:11], index - 1),
                    "Timo_adjacent_gap": sided_min_gap(timo.values[:11], index - 1),
                    "EB_sigma_1": eb_item.sigma_1 if eb_item else float("nan"),
                    "Timo_sigma_1": timo_item.sigma_1 if timo_item else float("nan"),
                    "EB_sigma_ratio": eb_item.sigma_ratio if eb_item else float("nan"),
                    "Timo_sigma_ratio": timo_item.sigma_ratio if timo_item else float("nan"),
                    "EB_root_quality_status": "pass" if eb_item and eb_item.sigma_1 <= eb.settings.sigma_accept and eb_item.sigma_ratio <= eb.settings.sigma_ratio_accept else "fail",
                    "Timo_root_quality_status": "pass" if timo_item and timo_item.sigma_1 <= timo.settings.sigma_accept and timo_item.sigma_ratio <= timo.settings.sigma_ratio_accept else "fail",
                    "EB_branch_reordered": branch_reordered(eb, eb_item, index) if eb_item else "",
                    "Timo_branch_reordered": branch_reordered(timo, timo_item, index) if timo_item else "",
                    "EB_root11_guard": eb.values[10] if len(eb.values) >= 11 else float("nan"),
                    "Timo_root11_guard": timo.values[10] if len(timo.values) >= 11 else float("nan"),
                    "strict_fallback_used": eb.operations.strict_fallback_runs > 0 or timo.operations.strict_fallback_runs > 0,
                    "Lambda_EB_verification": float("nan"),
                    "Lambda_Timo_verification": float("nan"),
                    "delta_f_verification": float("nan"),
                    "verification_residual_status": "not_run",
                }
            )
        target_set = set(case.prefixes)
        for item in prefix_metrics:
            primary = metrics.primary_status(
                V_n=item.V_n,
                N_true=n_true,
                required_prefix_n=item.prefix_n,
                N_certified_0=n_certified,
                strict_margin=args.strict_verification_margin,
                k10_resolved=included,
            )
            prefix_rows.append(
                {
                    "case_id": case.case_id,
                    "prefix_group": case.prefix_group,
                    "epsilon_source": case.epsilon_source,
                    "case_group": case.case_group,
                    "epsilon": case.epsilon,
                    "beta_deg": case.beta_deg,
                    "mu": case.mu,
                    "eta": case.eta,
                    "prefix_n": item.prefix_n,
                    "is_target_prefix": item.prefix_n in target_set,
                    "is_full_certificate_probe": item.prefix_n == n_certified,
                    "N_certified_0": n_certified,
                    "N_true_primary": n_true,
                    "Delta_n": item.Delta_n,
                    "V_n": item.V_n,
                    "M_n": item.M_n,
                    "triggering_sorted_indices": item.triggering_indices,
                    "target_prefix_pass": item.Delta_n <= args.threshold if math.isfinite(item.Delta_n) else False,
                    "full_certificate_pass": n_true >= n_certified,
                    "certificate_overprediction": max(n_certified - n_true, 0),
                    "primary_status": primary,
                    "quality_triggers": quality_triggers,
                    "same_geometry_near_buffer_pair": same_pair,
                    "Delta_n_verification": float("nan"),
                    "V_n_verification": float("nan"),
                    "M_n_verification": float("nan"),
                    "N_true_verification": "",
                    "verification_spread": float("nan"),
                    "roots_max_abs_difference": float("nan"),
                    "roots_max_rel_difference": float("nan"),
                    "primary_K10_resolved": included,
                    "verification_K10_resolved": "",
                    "cluster_agreement": "",
                    "strict_verification_required": False,
                    "strict_verification_status": "not_triggered",
                    "final_status": "pending",
                }
            )
    return spectrum_rows, mode_rows, prefix_rows, calculations, operation_rows


def audit_baseline_controls(
    cases: Sequence[ManifestCase],
    spectra_by_case: Mapping[str, Mapping[str, branch.BranchContinuationResult]],
    calculations: Mapping[str, Mapping[str, object]],
    baseline: BaselineData,
    args: Args,
) -> tuple[list[dict[str, object]], set[str]]:
    rows: list[dict[str, object]] = []
    failed_cases: set[str] = set()
    controls = [case for case in cases if case.case_group == "baseline_control"]
    by_group_source = {(case.prefix_group, case.epsilon_source): case for case in controls}
    for case in controls:
        calculation = calculations[case.case_id]
        spectra = spectra_by_case[case.case_id]
        baseline_deltas = baseline.mode_deltas[case.epsilon]
        baseline_prefix = metrics.running_prefix_metrics(baseline_deltas, threshold=args.threshold)
        baseline_n_true = metrics.true_safe_prefix(baseline_deltas, threshold=args.threshold)
        observed_deltas = calculation["deltas"]
        for prefix_n in case.prefixes:
            root_differences: dict[str, float] = {}
            model_root_pass: dict[str, bool] = {}
            for model in complete.SUPPORTED_MODELS:
                observed = spectra[model].values[:11]
                oracle = baseline.factorized_roots[(case.epsilon, model)]
                maximum = max((abs(left - right) for left, right in zip(observed, oracle)), default=float("inf"))
                root_differences[model] = maximum
                model_root_pass[model] = len(observed) == 11 and maximum <= spectra[model].settings.root_match_tolerance
            delta_max_abs = max(
                (abs(float(left) - float(right)) for left, right in zip(observed_deltas, baseline_deltas)),
                default=float("inf"),
            )
            near_case = by_group_source[(case.prefix_group, "epsilon_near_n")]
            buffer_case = by_group_source[(case.prefix_group, "epsilon_buffer_n")]
            near_metric = calculations[near_case.case_id]["prefix_metrics"][prefix_n - 1]  # type: ignore[index]
            buffer_metric = calculations[buffer_case.case_id]["prefix_metrics"][prefix_n - 1]  # type: ignore[index]
            near_closer = abs(near_metric.V_n) <= abs(buffer_metric.V_n) + args.violation_tolerance
            observed_metric = calculation["prefix_metrics"][prefix_n - 1]  # type: ignore[index]
            passed = (
                all(model_root_pass.values())
                and delta_max_abs <= max(1.0e-10, args.violation_tolerance)
                and int(calculation["N_true"]) == baseline_n_true
                and observed_metric.Delta_n <= args.threshold
                and baseline_prefix[prefix_n - 1].Delta_n <= args.threshold
                and near_closer
            )
            if not passed:
                failed_cases.add(case.case_id)
            rows.append(
                {
                    "case_id": case.case_id,
                    "prefix_group": case.prefix_group,
                    "prefix_n": prefix_n,
                    "epsilon_source": case.epsilon_source,
                    "epsilon": case.epsilon,
                    "oracle_algorithm_version": baseline.algorithm_version,
                    "EB_roots_1_11_max_abs_difference": root_differences[complete.MODEL_EB],
                    "Timo_roots_1_11_max_abs_difference": root_differences[complete.MODEL_TIMO],
                    "EB_roots_1_11_agree": model_root_pass[complete.MODEL_EB],
                    "Timo_roots_1_11_agree": model_root_pass[complete.MODEL_TIMO],
                    "Delta_n_observed": observed_metric.Delta_n,
                    "Delta_n_corrected_baseline": baseline_prefix[prefix_n - 1].Delta_n,
                    "delta_first10_max_abs_difference": delta_max_abs,
                    "N_true_observed": calculation["N_true"],
                    "N_true_corrected_baseline": baseline_n_true,
                    "target_prefix_safe": observed_metric.Delta_n <= args.threshold,
                    "near_closer_than_buffer_same_prefix": near_closer,
                    "same_geometry_near_buffer_pair": True,
                    "baseline_control_status": "pass" if passed else "baseline_control_pipeline_failure",
                }
            )
    return rows, failed_cases


def root_comparison(
    primary: Mapping[str, branch.BranchContinuationResult],
    verification: Mapping[str, branch.BranchContinuationResult],
) -> tuple[float, float, bool, bool]:
    maximum_abs = 0.0
    maximum_rel = 0.0
    roots_agree = True
    clusters_agree = True
    for model in complete.SUPPORTED_MODELS:
        left = primary[model]
        right = verification[model]
        left_values = left.values[:11]
        right_values = right.values[:11]
        if len(left_values) != 11 or len(right_values) != 11:
            return float("inf"), float("inf"), False, False
        differences = [abs(a - b) for a, b in zip(left_values, right_values)]
        relatives = [difference / max(abs(b), 1.0e-30) for difference, b in zip(differences, right_values)]
        maximum_abs = max(maximum_abs, max(differences, default=0.0))
        maximum_rel = max(maximum_rel, max(relatives, default=0.0))
        roots_agree = roots_agree and max(differences, default=float("inf")) <= left.settings.root_match_tolerance
        left_ordered = ordered_branches(left)[:11]
        right_ordered = ordered_branches(right)[:11]
        left_signature = tuple(cluster_size(item, left_ordered) for item in left_ordered)
        right_signature = tuple(cluster_size(item, right_ordered) for item in right_ordered)
        clusters_agree = clusters_agree and left_signature == right_signature
    return maximum_abs, maximum_rel, roots_agree, clusters_agree


def select_strict_cases(
    cases: Sequence[ManifestCase],
    prefix_rows: Sequence[Mapping[str, object]],
    calculations: Mapping[str, Mapping[str, object]],
    baseline_failures: set[str],
) -> dict[str, tuple[str, ...]]:
    reasons: dict[str, set[str]] = defaultdict(set)
    for row in prefix_rows:
        if not (bool(row["is_target_prefix"]) or bool(row["is_full_certificate_probe"])):
            continue
        status = str(row["primary_status"])
        if status == "provisional_counterexample":
            reasons[str(row["case_id"])].add("provisional_counterexample")
        elif status == "primary_near_boundary":
            reasons[str(row["case_id"])].add("primary_near_boundary")
    for case in cases:
        for reason in calculations[case.case_id]["quality_triggers"]:  # type: ignore[union-attr]
            reasons[case.case_id].add(str(reason))
        if case.case_id in baseline_failures:
            reasons[case.case_id].add("baseline_control_pipeline_failure")
    for group in metrics.EXPECTED_PREFIX_GROUPS:
        candidates = [
            row
            for row in prefix_rows
            if row["prefix_group"] == group
            and bool(row["is_target_prefix"])
            and math.isfinite(float(row["V_n"]))
        ]
        worst = metrics.worst_row(candidates)
        if worst is not None:
            reasons[str(worst["case_id"])].add("force_sample_worst_prefix_group")
    return {case_id: tuple(sorted(values)) for case_id, values in reasons.items()}


def apply_strict_verification(
    cases: Sequence[ManifestCase],
    primary_spectra: Mapping[str, Mapping[str, branch.BranchContinuationResult]],
    prefix_rows: list[dict[str, object]],
    mode_rows: list[dict[str, object]],
    calculations: dict[str, dict[str, object]],
    strict_reasons: Mapping[str, tuple[str, ...]],
    baseline_failures: set[str],
    args: Args,
    resolver: Callable[[ManifestCase, branch.ContinuationSettings, branch.BranchContinuationCache], Mapping[str, branch.BranchContinuationResult]],
) -> tuple[
    dict[str, Mapping[str, branch.BranchContinuationResult]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
]:
    strict_spectra: dict[str, Mapping[str, branch.BranchContinuationResult]] = {}
    strict_rows: list[dict[str, object]] = []
    spectrum_rows: list[dict[str, object]] = []
    operation_rows: list[dict[str, object]] = []
    case_lookup = {case.case_id: case for case in cases}
    cache = branch.BranchContinuationCache(
        args.verification_cache_dir,
        reuse_cache=False,
        force_recompute=True,
        verification_scope="primary",
    )
    settings = verification_settings(args)
    for case_id in sorted(strict_reasons):
        case = case_lookup[case_id]
        resolved = dict(resolver(case, settings, cache))
        if resolved is primary_spectra[case_id]:
            raise RuntimeError("strict verification reused the primary spectrum mapping")
        if any(resolved[model] is primary_spectra[case_id][model] for model in complete.SUPPORTED_MODELS):
            raise RuntimeError("strict verification reused a primary result object")
        strict_spectra[case_id] = resolved
        for result in resolved.values():
            spectrum_rows.append(build_spectrum_summary_row(case, result, run_scope="strict_force_recompute"))
            operation_rows.append(
                {
                    "case_id": case_id,
                    "model": result.model,
                    "cost_scope": "strict_verification",
                    **operation_mapping(result),
                }
            )
        eb = resolved[complete.MODEL_EB]
        timo = resolved[complete.MODEL_TIMO]
        k10 = eb.k10_guard_resolved and timo.k10_guard_resolved and len(eb.values) >= 11 and len(timo.values) >= 11
        if k10:
            deltas = tuple(
                metrics.squared_frequency_delta(eb.values[index], timo.values[index])
                for index in range(args.k_max)
            )
            prefix_metrics = metrics.running_prefix_metrics(deltas, threshold=args.threshold)
            n_true = metrics.true_safe_prefix(deltas, threshold=args.threshold)
        else:
            deltas = tuple(float("nan") for _ in range(args.k_max))
            prefix_metrics = metrics.running_prefix_metrics(deltas, threshold=args.threshold)
            n_true = 0
        maximum_abs, maximum_rel, roots_agree, clusters_agree = root_comparison(
            primary_spectra[case_id], resolved
        )
        calculations[case_id]["strict_deltas"] = deltas
        calculations[case_id]["strict_prefix_metrics"] = prefix_metrics
        calculations[case_id]["strict_N_true"] = n_true
        calculations[case_id]["strict_K10"] = k10
        calculations[case_id]["roots_max_abs_difference"] = maximum_abs
        calculations[case_id]["roots_max_rel_difference"] = maximum_rel
        calculations[case_id]["roots_agree"] = roots_agree
        calculations[case_id]["clusters_agree"] = clusters_agree
        for item in prefix_metrics:
            primary_item = calculations[case_id]["prefix_metrics"][item.prefix_n - 1]  # type: ignore[index]
            strict_rows.append(
                {
                    "case_id": case_id,
                    "prefix_n": item.prefix_n,
                    "verification_reasons": strict_reasons[case_id],
                    "independent_settings": asdict(settings),
                    "separate_cache_directory": display_path(args.verification_cache_dir),
                    "force_recompute": True,
                    "Delta_n_primary": primary_item.Delta_n,
                    "Delta_n_verification": item.Delta_n,
                    "V_n_primary": primary_item.V_n,
                    "V_n_verification": item.V_n,
                    "verification_spread": abs(primary_item.V_n - item.V_n),
                    "N_true_primary": calculations[case_id]["N_true"],
                    "N_true_verification": n_true,
                    "roots_max_abs_difference": maximum_abs,
                    "roots_max_rel_difference": maximum_rel,
                    "primary_K10_resolved": calculations[case_id]["included"],
                    "verification_K10_resolved": k10,
                    "roots_agree": roots_agree,
                    "cluster_agreement": clusters_agree,
                    "EB_force_global_full_SVD_agreement": eb.force_verification_agreement,
                    "Timo_force_global_full_SVD_agreement": timo.force_verification_agreement,
                    "verification_status": "pass" if k10 and roots_agree and clusters_agree else "fail",
                }
            )

    mode_lookup = {(str(row["case_id"]), int(row["sorted_index"])): row for row in mode_rows}
    for case_id, spectra in strict_spectra.items():
        eb = spectra[complete.MODEL_EB]
        timo = spectra[complete.MODEL_TIMO]
        strict_deltas = calculations[case_id]["strict_deltas"]
        for index in range(1, args.k_max + 1):
            row = mode_lookup[(case_id, index)]
            row["Lambda_EB_verification"] = eb.values[index - 1] if index <= len(eb.values) else float("nan")
            row["Lambda_Timo_verification"] = timo.values[index - 1] if index <= len(timo.values) else float("nan")
            row["delta_f_verification"] = strict_deltas[index - 1]  # type: ignore[index]
            row["verification_residual_status"] = "pass" if calculations[case_id]["strict_K10"] else "fail"

    for row in prefix_rows:
        case_id = str(row["case_id"])
        strict_required = case_id in strict_reasons
        row["strict_verification_required"] = strict_required
        if case_id in baseline_failures:
            row["primary_status"] = "baseline_control_failure"
        if strict_required and case_id in strict_spectra:
            item = calculations[case_id]["strict_prefix_metrics"][int(row["prefix_n"]) - 1]  # type: ignore[index]
            row["Delta_n_verification"] = item.Delta_n
            row["V_n_verification"] = item.V_n
            row["M_n_verification"] = item.M_n
            row["N_true_verification"] = calculations[case_id]["strict_N_true"]
            row["verification_spread"] = abs(float(row["V_n"]) - item.V_n)
            row["roots_max_abs_difference"] = calculations[case_id]["roots_max_abs_difference"]
            row["roots_max_rel_difference"] = calculations[case_id]["roots_max_rel_difference"]
            row["verification_K10_resolved"] = calculations[case_id]["strict_K10"]
            row["cluster_agreement"] = calculations[case_id]["clusters_agree"]
            row["strict_verification_status"] = (
                "pass"
                if calculations[case_id]["strict_K10"]
                and calculations[case_id]["roots_agree"]
                and calculations[case_id]["clusters_agree"]
                else "fail"
            )
            verification_v: float | None = item.V_n
            verification_n: int | None = int(calculations[case_id]["strict_N_true"])
            verification_k10: bool | None = bool(calculations[case_id]["strict_K10"])
            roots_agree: bool | None = bool(calculations[case_id]["roots_agree"])
            clusters_agree: bool | None = bool(calculations[case_id]["clusters_agree"])
        else:
            verification_v = None
            verification_n = None
            verification_k10 = None
            roots_agree = None
            clusters_agree = None
        row["final_status"] = metrics.final_status(
            primary=str(row["primary_status"]),
            V_primary=float(row["V_n"]),
            V_verification=verification_v,
            N_true_primary=int(row["N_true_primary"]),
            N_true_verification=verification_n,
            required_prefix_n=int(row["prefix_n"]),
            N_certified_0=int(row["N_certified_0"]),
            primary_k10=bool(row["primary_K10_resolved"]),
            verification_k10=verification_k10,
            roots_agree=roots_agree,
            clusters_agree=clusters_agree,
            violation_tolerance=args.violation_tolerance,
            strict_required=strict_required,
        )
    return strict_spectra, strict_rows, spectrum_rows, operation_rows


def build_case_summary(
    cases: Sequence[ManifestCase],
    prefix_rows: Sequence[Mapping[str, object]],
    calculations: Mapping[str, Mapping[str, object]],
) -> list[dict[str, object]]:
    by_case: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in prefix_rows:
        by_case[str(row["case_id"])].append(row)
    rows: list[dict[str, object]] = []
    priority = {
        "confirmed_counterexample": 5,
        "numerically_indeterminate_near_threshold": 4,
        "unresolved_spectrum": 3,
        "confirmed_safe_at_screen_point": 2,
        "screen_safe_not_strictly_recomputed": 1,
    }
    for case in cases:
        relevant = [
            row
            for row in by_case[case.case_id]
            if bool(row["is_target_prefix"]) or bool(row["is_full_certificate_probe"])
        ]
        final = max((str(row["final_status"]) for row in relevant), key=lambda status: priority.get(status, 0))
        n_true = int(calculations[case.case_id]["N_true"])
        deltas = calculations[case.case_id]["deltas"]
        rows.append(
            {
                "case_id": case.case_id,
                "prefix_group": case.prefix_group,
                "epsilon_source": case.epsilon_source,
                "case_group": case.case_group,
                "epsilon": case.epsilon,
                "beta_deg": case.beta_deg,
                "mu": case.mu,
                "eta": case.eta,
                "N_certified_0": calculations[case.case_id]["N_certified_0"],
                "N_true": n_true,
                "first_failed_mode": n_true + 1 if n_true < 10 else "",
                "late_pass_indices": calculations[case.case_id]["late_pass_indices"],
                "maximum_delta_first10": max(deltas) if all(math.isfinite(float(value)) for value in deltas) else float("nan"),
                "K10_status": "included" if calculations[case.case_id]["included"] else "unresolved",
                "strict_verification_triggered": case.case_id in {
                    str(row["case_id"]) for row in relevant if bool(row["strict_verification_required"])
                },
                "same_geometry_near_buffer_pair": calculations[case.case_id]["same_geometry_near_buffer_pair"],
                "final_scientific_status": final,
            }
        )
    return rows


def build_group_summary(prefix_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    grouped: dict[tuple[str, int, str], list[Mapping[str, object]]] = defaultdict(list)
    for row in prefix_rows:
        if row["case_group"] == "adversarial" and bool(row["is_target_prefix"]):
            grouped[(str(row["prefix_group"]), int(row["prefix_n"]), str(row["epsilon_source"]))].append(row)
    output: list[dict[str, object]] = []
    for (group, prefix_n, source), rows in sorted(grouped.items()):
        worst = metrics.worst_row(rows)
        assert worst is not None
        output.append(
            {
                "prefix_group": group,
                "prefix_n": prefix_n,
                "epsilon_source": source,
                "adversarial_case_count": len(rows),
                "resolved_count": sum(row["final_status"] != "unresolved_spectrum" for row in rows),
                "confirmed_counterexample_count": sum(row["final_status"] == "confirmed_counterexample" for row in rows),
                "indeterminate_count": sum(row["final_status"] == "numerically_indeterminate_near_threshold" for row in rows),
                "minimum_margin": min(float(row["M_n"]) for row in rows),
                "maximum_violation": max(float(row["V_n"]) for row in rows),
                "worst_case_id": worst["case_id"],
                "worst_beta_deg": worst["beta_deg"],
                "worst_mu": worst["mu"],
                "worst_eta": worst["eta"],
                "worst_N_true": worst["N_true_primary"],
                "triggering_sorted_indices": worst["triggering_sorted_indices"],
                "quality_triggers": worst["quality_triggers"],
                "worst_final_status": worst["final_status"],
            }
        )
    return output


def build_counterexample_audit(prefix_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for row in prefix_rows:
        if not (bool(row["is_target_prefix"]) or bool(row["is_full_certificate_probe"])):
            continue
        primary = str(row["primary_status"])
        final = str(row["final_status"])
        if primary != "provisional_counterexample" and final not in {
            "confirmed_counterexample",
            "numerically_indeterminate_near_threshold",
        }:
            continue
        if final == "confirmed_counterexample":
            disposition = "confirmed"
        elif final == "numerically_indeterminate_near_threshold":
            disposition = "indeterminate"
        elif primary == "provisional_counterexample":
            disposition = "rejected" if final.startswith("confirmed_safe") else "provisional"
        else:
            disposition = "provisional"
        rows.append(
            {
                "case_id": row["case_id"],
                "prefix_n": row["prefix_n"],
                "candidate_disposition": disposition,
                "primary_status": primary,
                "final_status": final,
                "V_n_primary": row["V_n"],
                "V_n_verification": row["V_n_verification"],
                "N_true_primary": row["N_true_primary"],
                "N_true_verification": row["N_true_verification"],
                "N_certified_0": row["N_certified_0"],
                "triggering_sorted_indices": row["triggering_sorted_indices"],
                "verification_status": row["strict_verification_status"],
                "reason": "target-prefix failure or corrected baseline-certificate overprediction candidate",
            }
        )
    return rows


def build_step3b_proposal(
    cases: Sequence[ManifestCase],
    prefix_rows: Sequence[Mapping[str, object]],
    thresholds: Mapping[int, BaselineThreshold],
) -> list[dict[str, object]]:
    del cases  # Proposal selection is CSV-metric-only and never calls a root solver.
    output: list[dict[str, object]] = []
    seen: set[tuple[int, str, float, float, float]] = set()
    for prefix_n in range(2, 11):
        candidates = [
            row
            for row in prefix_rows
            if int(row["prefix_n"]) == prefix_n
            and bool(row["is_target_prefix"])
            and row["case_group"] == "adversarial"
            and row["final_status"] != "unresolved_spectrum"
            and math.isfinite(float(row.get("V_n_verification", float("nan"))))
        ]
        if not candidates:
            continue
        worst = max(candidates, key=lambda row: float(row["V_n_verification"]))
        purpose = (
            "refine_counterexample_epsilon_star"
            if worst["final_status"] == "confirmed_counterexample"
            else "paired_worst_margin_check"
        )
        threshold = thresholds[prefix_n]
        for source in ("epsilon_near_n", "epsilon_buffer_n"):
            epsilon = getattr(threshold, source)
            key = (prefix_n, source, float(worst["beta_deg"]), float(worst["mu"]), float(worst["eta"]))
            if key in seen:
                continue
            seen.add(key)
            output.append(
                {
                    "case_id": f"S3B_n{prefix_n:02d}_{'near' if source == 'epsilon_near_n' else 'buffer'}_01",
                    "source_step3a_case_id": worst["case_id"],
                    "prefix_n": prefix_n,
                    "prefix_group": worst["prefix_group"],
                    "epsilon_source": source,
                    "epsilon": epsilon,
                    "beta_deg": worst["beta_deg"],
                    "mu": worst["mu"],
                    "eta": worst["eta"],
                    "followup_purpose": purpose,
                    "observed_Delta_n": worst["Delta_n_verification"],
                    "observed_margin": -float(worst["V_n_verification"]),
                    "observed_violation": worst["V_n_verification"],
                    "strict_verification_status": worst["strict_verification_status"],
                    "required_spectrum_method": branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
                    "required_K10_guard": True,
                    "notes": "Step-3B proposal only; paired epsilon row has not been executed",
                }
            )
    return output


def verified_value(row: Mapping[str, object], key: str, fallback: str) -> float:
    value = finite_or_nan(row.get(key))
    return value if math.isfinite(value) else finite_or_nan(row.get(fallback))


def prefix_worst_rows(prefix_rows: Sequence[Mapping[str, object]]) -> list[Mapping[str, object]]:
    output: list[Mapping[str, object]] = []
    for prefix_n in range(1, 11):
        candidates = [
            row
            for row in prefix_rows
            if int(row["prefix_n"]) == prefix_n
            and bool(row["is_target_prefix"])
            and row["final_status"] != "unresolved_spectrum"
        ]
        if candidates:
            output.append(
                max(
                    candidates,
                    key=lambda row: verified_value(row, "V_n_verification", "V_n"),
                )
            )
    return output


def create_plots(output_dir: Path) -> tuple[Path, ...]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prefix_rows = read_csv(output_dir / "step3a_prefix_screening.csv")
    case_rows = read_csv(output_dir / "step3a_case_summary.csv")
    mode_rows = read_csv(output_dir / "step3a_mode_metrics.csv")
    operation_rows = read_csv(output_dir / "step3a_operation_counts.csv")
    target = [row for row in prefix_rows if parse_bool(row["is_target_prefix"])]
    worst = prefix_worst_rows(target)
    paths: list[Path] = []

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    for group, marker in (("baseline_control", "o"), ("adversarial", "x")):
        rows = [row for row in target if row["case_group"] == group]
        ax.scatter(
            [int(row["prefix_n"]) for row in rows],
            [verified_value(row, "M_n_verification", "M_n") for row in rows],
            marker=marker,
            label=group,
            alpha=0.75,
        )
    for row in worst:
        ax.annotate(str(row["case_id"]), (int(row["prefix_n"]), verified_value(row, "M_n_verification", "M_n")), fontsize=7)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set(xlabel="prefix n", ylabel=r"verified margin $M_n$", title="Step 3A margin by prefix")
    ax.legend()
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[0]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(11.0, 4.8))
    colors = {
        "confirmed_counterexample": "tab:red",
        "numerically_indeterminate_near_threshold": "tab:orange",
        "unresolved_spectrum": "tab:gray",
        "confirmed_safe_at_screen_point": "tab:green",
        "screen_safe_not_strictly_recomputed": "tab:blue",
    }
    ax.scatter(
        range(len(target)),
        [verified_value(row, "V_n_verification", "V_n") for row in target],
        c=[colors.get(row["final_status"], "tab:blue") for row in target],
        s=24,
    )
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xticks(range(len(target)))
    ax.set_xticklabels([f"{row['case_id']}:n{row['prefix_n']}" for row in target], rotation=90, fontsize=6)
    ax.set(ylabel=r"verified violation $V_n$", title="Target-prefix violations")
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[1]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(6.0, 5.5))
    n_cert = [int(row["N_certified_0"]) for row in case_rows]
    n_true = [int(row["N_true"]) for row in case_rows]
    ax.scatter(n_cert, n_true, c=["tab:blue" if row["case_group"] == "adversarial" else "tab:green" for row in case_rows])
    ax.plot([0, 10], [0, 10], color="black", linewidth=0.8)
    ax.set(xlabel=r"baseline certificate $N^{(0)}_{certified}$", ylabel=r"observed $N_{true}$", xlim=(0, 10.5), ylim=(0, 10.5))
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[2]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    fig, ax = plt.subplots(figsize=(10.0, 4.8))
    for source, marker in (("epsilon_near_n", "o"), ("epsilon_buffer_n", "s")):
        rows = [row for row in target if row["epsilon_source"] == source]
        ax.scatter(
            range(len(rows)),
            [verified_value(row, "M_n_verification", "M_n") for row in rows],
            marker=marker,
            label=source,
        )
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set(title="Near/buffer screen (adversarial rows are not generally geometry-paired)", ylabel=r"$M_n$")
    ax.legend()
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[3]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    mode_lookup = {(row["case_id"], int(row["sorted_index"])): row for row in mode_rows}
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    for row in target:
        triggers = json.loads(row["triggering_sorted_indices"])
        for trigger in triggers:
            mode = mode_lookup.get((row["case_id"], int(trigger)), {})
            clustered = max(int(mode.get("EB_cluster_size", 1)), int(mode.get("Timo_cluster_size", 1)))
            ax.scatter(int(row["prefix_n"]), int(trigger), s=22 + 12 * clustered, alpha=0.65)
    ax.set(xlabel="prefix n", ylabel="triggering sorted mode", title="Trigger-mode map; marker size reflects cluster size")
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[4]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)

    scopes = ("primary_screening", "strict_verification", "baseline_control_oracle", "report_postprocessing")
    matrix_fields = (
        "parent_matrix_evaluations",
        "local_matrix_evaluations",
        "guard_matrix_evaluations",
        "strict_fallback_matrix_evaluations",
        "force_verification_matrix_evaluations",
    )
    svd_fields = (
        "parent_block_SVD_calls",
        "local_full_6x6_SVD_calls",
        "guard_full_6x6_SVD_calls",
        "strict_fallback_full_6x6_SVD_calls",
        "force_verification_full_6x6_SVD_calls",
    )
    matrix_totals = [
        sum(int(float(row.get(field, 0) or 0)) for row in operation_rows if row["cost_scope"] == scope for field in matrix_fields)
        for scope in scopes
    ]
    svd_totals = [
        sum(int(float(row.get(field, 0) or 0)) for row in operation_rows if row["cost_scope"] == scope for field in svd_fields)
        for scope in scopes
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
    axes[0].bar(range(len(scopes)), matrix_totals)
    axes[0].set_title("Matrix evaluations")
    axes[1].bar(range(len(scopes)), svd_totals)
    axes[1].set_title("SVD calls")
    for axis in axes:
        axis.set_xticks(range(len(scopes)))
        axis.set_xticklabels(scopes, rotation=35, ha="right", fontsize=7)
    fig.suptitle("Operation costs by scope (separate primitives; no combined score)")
    fig.tight_layout()
    path = output_dir / PLOT_NAMES[5]
    fig.savefig(path, dpi=180)
    plt.close(fig)
    paths.append(path)
    return tuple(paths)


def operation_totals(rows: Sequence[Mapping[str, object]], scope: str) -> dict[str, int]:
    fields = tuple(branch.BranchOperationCounts.__dataclass_fields__) + (
        "delta_evaluations",
        "prefix_maximum_operations",
        "margin_violation_comparisons",
        "aggregation_comparisons",
        "root_comparison_operations",
    )
    return {
        field: sum(int(float(row.get(field, 0) or 0)) for row in rows if row.get("cost_scope") == scope)
        for field in fields
    }


def write_report(
    args: Args,
    *,
    cases: Sequence[ManifestCase],
    prefix_rows: Sequence[Mapping[str, object]],
    case_rows: Sequence[Mapping[str, object]],
    spectrum_rows: Sequence[Mapping[str, object]],
    baseline_rows: Sequence[Mapping[str, object]],
    counterexample_rows: Sequence[Mapping[str, object]],
    exclusion_rows: Sequence[Mapping[str, object]],
    operation_rows: Sequence[Mapping[str, object]],
    proposal_rows: Sequence[Mapping[str, object]],
    strict_reasons: Mapping[str, tuple[str, ...]],
    decision: str,
    elapsed: float,
) -> Path:
    target = [row for row in prefix_rows if bool(row["is_target_prefix"])]
    confirmed = [row for row in counterexample_rows if row["candidate_disposition"] == "confirmed"]
    provisional = [row for row in counterexample_rows if row["candidate_disposition"] == "provisional"]
    rejected = [row for row in counterexample_rows if row["candidate_disposition"] == "rejected"]
    indeterminate = [row for row in counterexample_rows if row["candidate_disposition"] == "indeterminate"]
    k10_geometries = sum(row["K10_status"] == "included" for row in case_rows)
    full12_model_rows = sum(bool(row["full12_resolved"]) for row in spectrum_rows)
    fallback_model_rows = sum(bool(row["strict_fallback_used"]) for row in spectrum_rows)
    baseline_pass = sum(row["baseline_control_status"] == "pass" for row in baseline_rows)
    worst = prefix_worst_rows(prefix_rows)
    primary_ops = operation_totals(operation_rows, "primary_screening")
    strict_ops = operation_totals(operation_rows, "strict_verification")
    lines = [
        "# EB epsilon lower-envelope Step-3A report",
        "",
        "## Scope",
        "",
        f"This run is Step 3A only: `{len(cases)}` fixed manifest cases, corrected straight-baseline thresholds, K=10 plus the mandatory root-11 guard, and Timoshenko as the existing 1D reference. No Step 3B calculation, FEM, formula, matrix, determinant, shared-solver, or article workflow was run or modified.",
        "",
        "## Manifest Provenance",
        "",
        f"The source manifest contains 28 validated rows and exactly seven explicit prefix groups. All epsilon values matched the full-precision corrected `{factorized.ALGORITHM_VERSION}` threshold source within `{args.baseline_provenance_tolerance:.3e}` before any root solve. The selected run contains `{sum(case.case_group == 'baseline_control' for case in cases)}` baseline and `{sum(case.case_group == 'adversarial' for case in cases)}` adversarial geometries.",
        "",
        "## Spectrum Quality",
        "",
        f"K10/root-11 inclusion passed for `{k10_geometries}/{len(cases)}` geometries. `full12_resolved` passed for `{full12_model_rows}/{2 * len(cases)}` primary case/model rows and remained diagnostic only. Strict fallback was used in `{fallback_model_rows}` primary case/model rows. Exclusions: `{len(exclusion_rows)}`. Strict verification was triggered for `{len(strict_reasons)}` geometries.",
        "",
        "## Baseline Controls",
        "",
        f"Corrected factorized-oracle checks passed for `{baseline_pass}/{len(baseline_rows)}` control/prefix rows.",
        "",
        "| case | n | source | Delta_n | N_true | oracle status |",
        "|---|---:|---|---:|---:|---|",
    ]
    for row in baseline_rows:
        lines.append(
            f"| {row['case_id']} | {row['prefix_n']} | {row['epsilon_source']} | {float(row['Delta_n_observed']):.8e} | {row['N_true_observed']} | {row['baseline_control_status']} |"
        )
    lines.extend(
        [
            "",
            "## Targeted Screening",
            "",
            "| case | n | epsilon source | epsilon | beta | mu | eta | Ncert0 | Ntrue | Delta_n | V_n | M_n | trigger | strict | final |",
            "|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|",
        ]
    )
    for row in target:
        lines.append(
            "| {case_id} | {prefix_n} | {epsilon_source} | {epsilon:.8e} | {beta_deg:g} | {mu:g} | {eta:g} | {N_certified_0} | {N_true_primary} | {Delta_n:.8e} | {V_n:.8e} | {M_n:.8e} | {triggering_sorted_indices} | {strict_verification_status} | {final_status} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Prefix-Wise Worst Cases",
            "",
            "| n | case | verified V_n | verified M_n | beta | mu | eta | trigger | status |",
            "|---:|---|---:|---:|---:|---:|---:|---|---|",
        ]
    )
    for row in worst:
        violation = verified_value(row, "V_n_verification", "V_n")
        lines.append(
            f"| {row['prefix_n']} | {row['case_id']} | {violation:.8e} | {-violation:.8e} | {float(row['beta_deg']):g} | {float(row['mu']):g} | {float(row['eta']):g} | {row['triggering_sorted_indices']} | {row['final_status']} |"
        )
    lines.extend(
        [
            "",
            "## Counterexample Assessment",
            "",
            f"Provisional: `{len(provisional)}`; confirmed: `{len(confirmed)}`; rejected after verification: `{len(rejected)}`; numerically indeterminate: `{len(indeterminate)}`.",
            "",
            "No local epsilon refinement was executed, including when a candidate was detected.",
            "",
            "## Near/Buffer Design Caveat",
            "",
            "Baseline near/buffer controls are paired at identical beta, mu, and eta. Adversarial near/buffer cases in the fixed manifest generally use different geometries, so they are independent screen points, not a pure-thickness trajectory. `same_geometry_near_buffer_pair` records this distinction.",
            "",
            "## Step-3A Decision",
            "",
            f"Decision: `{decision}`.",
            "",
        ]
    )
    if not confirmed:
        lines.append("No counterexample was found in the executed 28-case targeted screen." if len(cases) == 28 else "No counterexample was found in this isolated smoke subset.")
    lines.extend(
        [
            "",
            "## Step-3B Proposal",
            "",
            f"The proposal contains `{len(proposal_rows)}` paired near/buffer rows selected from verified worst margins. It is a proposal only and was not executed.",
            "",
            "## Operation Counts",
            "",
            f"Primary matrix evaluations (local + guard + fallbacks): `{primary_ops['local_matrix_evaluations'] + primary_ops['guard_matrix_evaluations'] + primary_ops['strict_fallback_matrix_evaluations']}`; primary full 6x6 SVD calls: `{primary_ops['local_full_6x6_SVD_calls'] + primary_ops['guard_full_6x6_SVD_calls'] + primary_ops['strict_fallback_full_6x6_SVD_calls']}`.",
            f"Strict matrix evaluations (including force-global checks): `{strict_ops['local_matrix_evaluations'] + strict_ops['guard_matrix_evaluations'] + strict_ops['strict_fallback_matrix_evaluations'] + strict_ops['force_verification_matrix_evaluations']}`; strict full 6x6 SVD calls: `{strict_ops['local_full_6x6_SVD_calls'] + strict_ops['guard_full_6x6_SVD_calls'] + strict_ops['strict_fallback_full_6x6_SVD_calls'] + strict_ops['force_verification_full_6x6_SVD_calls']}`. Wall time `{elapsed:.3f}` s is auxiliary metadata only.",
            "",
            "## Limitations",
            "",
            "This is a finite targeted screen, not a continuous-domain proof. It uses selected adversarial geometries, no optimizer, no non-baseline epsilon-star refinement, and the existing 1D Timoshenko reference. Root completeness is numerically audited rather than mathematically proved. No 3D FEM rerun was performed.",
            "",
            "## Scientific Wording",
            "",
            "The result must not be described as proof that the straight baseline is a global lower envelope. It supports only the finite executed design stated above.",
        ]
    )
    path = args.output_dir / "eb_epsilon_lower_envelope_step3a_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def ensure_output_policy(args: Args) -> None:
    existing = [args.output_dir / name for name in (*OUTPUT_NAMES, *PLOT_NAMES) if (args.output_dir / name).exists()]
    if existing and not args.force:
        raise FileExistsError("Step-3A outputs already exist; pass --force to replace generated products")
    args.output_dir.mkdir(parents=True, exist_ok=True)


def plot_only(args: Args) -> dict[str, object]:
    required = [args.output_dir / name for name in OUTPUT_NAMES[:11]]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(f"plot-only requires saved Step-3A CSVs; missing: {missing}")
    paths = create_plots(args.output_dir)
    return {"args": args, "plot_only": True, "root_calculations": 0, "plots": paths}


def execute(
    args: Args,
    *,
    resolver: Callable[[ManifestCase, branch.ContinuationSettings, branch.BranchContinuationCache], Mapping[str, branch.BranchContinuationResult]] = resolve_case_models,
) -> dict[str, object]:
    if args.plot_only:
        return plot_only(args)
    started = time.perf_counter()

    # The complete contract and corrected provenance are validated before the
    # first resolver call.  Smoke selection happens only after these checks.
    all_cases = load_manifest(args.manifest)
    validate_manifest_contract(all_cases)
    thresholds = load_thresholds(args.baseline_dir)
    manifest_rows = validate_manifest_provenance(
        all_cases,
        thresholds,
        tolerance=args.baseline_provenance_tolerance,
        baseline_dir=args.baseline_dir,
    )
    validate_gateway_provenance(args.gateway_dir)
    baseline = load_baseline_data(
        args.baseline_dir,
        all_cases,
        thresholds,
        tolerance=args.baseline_provenance_tolerance,
    )
    cases = (
        [case for case in all_cases if case.case_id in set(SMOKE_CASE_IDS)]
        if args.smoke
        else list(all_cases)
    )
    ensure_output_policy(args)

    settings = primary_settings(args)
    primary_cache = branch.BranchContinuationCache(
        args.cache_dir,
        reuse_cache=args.reuse_cache,
        force_recompute=args.force_recompute,
        verification_scope="primary",
    )
    primary_spectra: dict[str, Mapping[str, branch.BranchContinuationResult]] = {}
    for case in cases:
        primary_spectra[case.case_id] = dict(resolver(case, settings, primary_cache))
        print(
            f"primary {case.case_id}: "
            + ", ".join(
                f"{model}={primary_spectra[case.case_id][model].spectrum_status}"
                for model in complete.SUPPORTED_MODELS
            )
        )

    spectrum_rows, mode_rows, prefix_rows, calculations, operation_rows = build_primary_products(
        cases, primary_spectra, baseline, args
    )
    baseline_rows, baseline_failures = audit_baseline_controls(
        cases, primary_spectra, calculations, baseline, args
    )
    strict_reasons = select_strict_cases(cases, prefix_rows, calculations, baseline_failures)
    strict_spectra, strict_rows, _strict_spectrum_rows, strict_operation_rows = apply_strict_verification(
        cases,
        primary_spectra,
        prefix_rows,
        mode_rows,
        calculations,
        strict_reasons,
        baseline_failures,
        args,
        resolver,
    )
    operation_rows.extend(strict_operation_rows)
    for case_id in strict_spectra:
        print(f"strict {case_id}: {'pass' if calculations[case_id]['strict_K10'] else 'fail'}")

    case_rows = build_case_summary(cases, prefix_rows, calculations)
    group_rows = build_group_summary(prefix_rows)
    counterexample_rows = build_counterexample_audit(prefix_rows)
    exclusion_rows: list[dict[str, object]] = []
    for case in cases:
        reasons: list[str] = []
        if not calculations[case.case_id]["included"]:
            for model in complete.SUPPORTED_MODELS:
                result = primary_spectra[case.case_id][model]
                if not result.k10_guard_resolved:
                    reasons.append(f"{model}:{result.exclusion_reason or result.spectrum_status}")
        if case.case_id in baseline_failures:
            reasons.append("baseline_control_pipeline_failure")
        if case.case_id in strict_spectra and not calculations[case.case_id]["strict_K10"]:
            reasons.append("strict_verification_K10_failure")
        if reasons:
            exclusion_rows.append(
                {
                    "case_id": case.case_id,
                    "geometry_id": case.geometry_id,
                    "exclusion_status": "unresolved_or_invalid",
                    "reasons": tuple(reasons),
                    "legacy_fallback_used": False,
                }
            )

    proposal_rows = (
        build_step3b_proposal(cases, prefix_rows, thresholds)
        if args.write_step3b_followup_manifest
        else []
    )
    baseline_controls_pass = bool(baseline_rows) and all(
        row["baseline_control_status"] == "pass" for row in baseline_rows
    )
    relevant_statuses = [
        str(row["final_status"])
        for row in prefix_rows
        if bool(row["is_target_prefix"]) or bool(row["is_full_certificate_probe"])
    ]
    strict_verified_count = sum(
        bool(calculations[case_id].get("strict_K10"))
        and bool(calculations[case_id].get("roots_agree"))
        and bool(calculations[case_id].get("clusters_agree"))
        for case_id in strict_reasons
    )
    decision = metrics.decide_step3a(
        geometry_count=len(cases),
        resolved_geometry_count=sum(bool(calculations[case.case_id]["included"]) for case in cases),
        baseline_controls_pass=baseline_controls_pass,
        final_target_statuses=relevant_statuses,
        strict_trigger_count=len(strict_reasons),
        strict_verified_count=strict_verified_count,
    )

    operation_rows.append(
        {
            "case_id": "baseline_controls",
            "model": "both",
            "cost_scope": "baseline_control_oracle",
            "root_comparison_operations": len(baseline_rows) * 22,
            "delta_evaluations": len(baseline_rows) * args.k_max,
            "prefix_maximum_operations": len(baseline_rows) * (args.k_max - 1),
            "margin_violation_comparisons": len(baseline_rows),
        }
    )
    operation_rows.append(
        {
            "case_id": "step3a_metrics",
            "model": "both",
            "cost_scope": "report_postprocessing",
            "delta_evaluations": len(cases) * args.k_max + len(strict_spectra) * args.k_max,
            "prefix_maximum_operations": (len(cases) + len(strict_spectra)) * (args.k_max - 1),
            "margin_violation_comparisons": len(prefix_rows) * 2,
            "aggregation_comparisons": len(prefix_rows) + len(group_rows),
            "root_comparison_operations": len(strict_spectra) * 2 * 11,
            "wall_seconds": time.perf_counter() - started,
        }
    )

    write_csv(args.output_dir / "step3a_manifest_resolved.csv", manifest_rows)
    write_csv(args.output_dir / "step3a_case_spectrum_summary.csv", spectrum_rows)
    write_csv(args.output_dir / "step3a_mode_metrics.csv", mode_rows)
    write_csv(args.output_dir / "step3a_case_summary.csv", case_rows)
    write_csv(args.output_dir / "step3a_prefix_screening.csv", prefix_rows)
    write_csv(args.output_dir / "step3a_prefix_group_summary.csv", group_rows)
    write_csv(args.output_dir / "step3a_baseline_control_audit.csv", baseline_rows)
    write_csv(args.output_dir / "step3a_strict_verification_audit.csv", strict_rows)
    write_csv(
        args.output_dir / "step3a_counterexample_audit.csv",
        counterexample_rows,
        fields=(
            "case_id",
            "prefix_n",
            "candidate_disposition",
            "primary_status",
            "final_status",
            "V_n_primary",
            "V_n_verification",
            "N_true_primary",
            "N_true_verification",
            "N_certified_0",
            "triggering_sorted_indices",
            "verification_status",
            "reason",
        ),
    )
    write_csv(
        args.output_dir / "step3a_exclusion_audit.csv",
        exclusion_rows,
        fields=("case_id", "geometry_id", "exclusion_status", "reasons", "legacy_fallback_used"),
    )
    write_csv(args.output_dir / "step3a_operation_counts.csv", operation_rows)
    proposal_fields = (
        "case_id",
        "source_step3a_case_id",
        "prefix_n",
        "prefix_group",
        "epsilon_source",
        "epsilon",
        "beta_deg",
        "mu",
        "eta",
        "followup_purpose",
        "observed_Delta_n",
        "observed_margin",
        "observed_violation",
        "strict_verification_status",
        "required_spectrum_method",
        "required_K10_guard",
        "notes",
    )
    write_csv(args.output_dir / "step3b_followup_case_manifest.csv", proposal_rows, fields=proposal_fields)
    elapsed = time.perf_counter() - started
    plots = create_plots(args.output_dir)
    report = write_report(
        args,
        cases=cases,
        prefix_rows=prefix_rows,
        case_rows=case_rows,
        spectrum_rows=spectrum_rows,
        baseline_rows=baseline_rows,
        counterexample_rows=counterexample_rows,
        exclusion_rows=exclusion_rows,
        operation_rows=operation_rows,
        proposal_rows=proposal_rows,
        strict_reasons=strict_reasons,
        decision=decision,
        elapsed=elapsed,
    )
    print(f"Step-3A decision: {decision}")
    print(f"output: {args.output_dir}")
    return {
        "args": args,
        "all_cases": all_cases,
        "cases": cases,
        "manifest_rows": manifest_rows,
        "primary_spectra": primary_spectra,
        "strict_spectra": strict_spectra,
        "prefix_rows": prefix_rows,
        "case_rows": case_rows,
        "baseline_rows": baseline_rows,
        "strict_rows": strict_rows,
        "counterexample_rows": counterexample_rows,
        "exclusion_rows": exclusion_rows,
        "operation_rows": operation_rows,
        "proposal_rows": proposal_rows,
        "strict_reasons": strict_reasons,
        "decision": decision,
        "plots": plots,
        "report": report,
        "elapsed_seconds": elapsed,
    }


def main(
    argv: Sequence[str] | None = None,
    *,
    resolver: Callable[[ManifestCase, branch.ContinuationSettings, branch.BranchContinuationCache], Mapping[str, branch.BranchContinuationResult]] = resolve_case_models,
) -> dict[str, object]:
    return execute(parse_args(argv), resolver=resolver)


if __name__ == "__main__":
    main()
