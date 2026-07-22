from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path
import sys
from typing import Callable, Mapping, Sequence

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

from my_project.analytic.formulas import (  # noqa: E402
    BeamParams,
    frequency_scale,
    lambdas_to_frequencies,
)
from scripts.lib import branch_informed_spectrum_continuation as branch  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402


DEFAULT_STEP3A_DIR = REPO_ROOT / "results" / "eb_epsilon_lower_envelope_step3a"
DEFAULT_GATEWAY_DIR = REPO_ROOT / "results" / "eb_timo_branch_continuation_gateway"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "results" / "eb_timo_counterexample_dimensional_frequency_beta"
SMOKE_OUTPUT_DIR = REPO_ROOT / "results" / "_smoke" / "eb_timo_counterexample_dimensional_frequency_beta"

MAIN_CSV = "counterexample_dimensional_frequency_beta.csv"
SCALE_CSV = "counterexample_frequency_scale_audit.csv"
QUALITY_CSV = "counterexample_frequency_beta_quality_audit.csv"
REPORT_MD = "counterexample_dimensional_frequency_beta_report.md"
PDF_NAMES = {
    "S3_12": "S3_12_dimensional_frequency_vs_beta.pdf",
    "S3_14": "S3_14_dimensional_frequency_vs_beta.pdf",
}
OUTPUT_NAMES = (MAIN_CSV, SCALE_CSV, QUALITY_CSV, REPORT_MD, *PDF_NAMES.values())

MODEL_EB = complete.MODEL_EB
MODEL_TIMO = complete.MODEL_TIMO
MODELS = (MODEL_EB, MODEL_TIMO)
FREQUENCY_KIND = "cyclic_frequency"
FREQUENCY_UNIT = "Hz"
FREQUENCY_FORMULA = (
    "f = Lambda**2 * sqrt(E*I/(rho*S)) / (2*pi*L_base**2)"
)
SCALING_SOURCE = "src/my_project/analytic/formulas.py"
SCALING_HELPERS = "BeamParams; frequency_scale; lambdas_to_frequencies"

# Canonical analytic diagnostic parameters already used by FreqMuNet.py,
# FreqFromAngle.py, and the epsilon-to-radius check in
# scripts/analysis/check_mu_to_one_single_rod_limit.py.  The radius itself is
# never chosen independently: BeamParams.eps = r/(2*L_base) is inverted below.
CANONICAL_E_PA = 2.1e11
CANONICAL_RHO_KG_M3 = 7800.0
CANONICAL_L_TOTAL_M = 2.0

EXPECTED_CASE_GEOMETRY = {
    "S3_12": {"mu": 0.7, "eta": 0.0, "counterexample_beta_deg": 90.0},
    "S3_14": {"mu": 0.5, "eta": -0.1, "counterexample_beta_deg": 45.0},
}

MAIN_FIELDS = (
    "case_id",
    "epsilon_0",
    "mu",
    "eta",
    "beta_deg",
    "sorted_index",
    "model",
    "branch_id",
    "parent_family",
    "Lambda",
    "dimensional_frequency",
    "frequency_kind",
    "frequency_unit",
    "K10_guard_resolved",
    "full12_resolved",
    "cluster_id",
    "adjacent_gap",
    "strict_fallback_used",
    "warnings",
)

QUALITY_FIELDS = (
    "case_id",
    "epsilon_0",
    "mu",
    "eta",
    "model",
    "beta_deg",
    "root_count",
    "K10_status",
    "K10_guard_resolved",
    "root11_status",
    "root11",
    "full12_status",
    "full12_resolved",
    "minimum_gap",
    "cluster_count",
    "fallback",
    "strict_fallback_used",
    "strict_fallback_count",
    "cache_status",
    "warnings",
)

SCALE_FIELDS = (
    "case_id",
    "epsilon_0",
    "scaling_formula",
    "material_constants",
    "length_radius_parameters",
    "E_Pa",
    "rho_kg_m3",
    "L_total_m",
    "L_base_m",
    "r0_m",
    "S0_m2",
    "I0_m4",
    "common_scale_factor",
    "frequency_kind",
    "unit",
    "source_helper_file",
    "source_helper",
    "EB_scale_factor",
    "Timoshenko_scale_factor",
    "EB_Timoshenko_scale_equality_check",
)


@dataclass(frozen=True)
class CaseSpec:
    case_id: str
    epsilon_0: float
    mu: float
    eta: float
    counterexample_beta_deg: float


@dataclass(frozen=True)
class Args:
    step3a_dir: Path
    gateway_dir: Path
    output_dir: Path
    beta_min: float
    beta_max: float
    beta_step: float
    k_max: int
    n_spectrum_roots: int
    n_candidate_roots: int
    verification_candidate_roots: int
    spectrum_method: str
    initial_beta_step: float
    min_beta_step: float
    max_beta_step: float
    beta_growth_factor: float
    beta_shrink_factor: float
    guard_scan_step: float
    seed_half_width: float
    sigma_accept: float
    sigma_ratio_accept: float
    mac_accept: float
    subspace_mac_accept: float
    cluster_gap_absolute: float
    cluster_gap_relative: float
    workers: int
    reuse_cache: bool
    force_recompute: bool
    plot_only: bool
    smoke: bool


@dataclass(frozen=True)
class SpectrumPoint:
    case_id: str
    model: str
    beta_deg: float
    branches: tuple[branch.ContinuedBranch, ...]
    k10_guard_resolved: bool
    full12_resolved: bool
    guard_status: str
    strict_fallback_used: bool
    strict_fallback_count: int
    cache_status: str
    warnings: tuple[str, ...]

    @property
    def values(self) -> tuple[float, ...]:
        return tuple(item.Lambda for item in self.branches)


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def fmt(value: object) -> object:
    if isinstance(value, (float, np.floating)):
        number = float(value)
        if math.isnan(number):
            return "nan"
        if math.isinf(number):
            return "inf" if number > 0.0 else "-inf"
        return f"{number:.16e}"
    if isinstance(value, (bool, np.bool_)):
        return "true" if bool(value) else "false"
    if isinstance(value, (tuple, list, dict)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def write_csv(
    path: Path,
    rows: Sequence[Mapping[str, object]],
    fields: Sequence[str],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in fields})
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: object) -> bool:
    token = str(value).strip().lower()
    if token in {"true", "1", "yes"}:
        return True
    if token in {"false", "0", "no", ""}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description=(
            "Plot the first ten sorted dimensional EB/Timoshenko frequencies "
            "for the two confirmed Step-3A counterexamples."
        ),
    )
    parser.add_argument("--step3a-dir", type=Path, default=DEFAULT_STEP3A_DIR)
    parser.add_argument("--gateway-dir", type=Path, default=DEFAULT_GATEWAY_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--beta-min", type=float, default=0.0)
    parser.add_argument("--beta-max", type=float, default=90.0)
    parser.add_argument("--beta-step", type=float, default=0.5)
    parser.add_argument("--k-max", type=int, default=10)
    parser.add_argument("--n-spectrum-roots", type=int, default=12)
    parser.add_argument("--n-candidate-roots", type=int, default=20)
    parser.add_argument("--verification-candidate-roots", type=int, default=24)
    parser.add_argument(
        "--spectrum-method",
        choices=(branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,),
        default=branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
    )
    parser.add_argument("--initial-beta-step", type=float, default=0.5)
    parser.add_argument("--min-beta-step", type=float, default=0.005)
    parser.add_argument("--max-beta-step", type=float, default=5.0)
    parser.add_argument("--beta-growth-factor", type=float, default=1.6)
    parser.add_argument("--beta-shrink-factor", type=float, default=0.5)
    parser.add_argument("--guard-scan-step", type=float, default=0.05)
    parser.add_argument("--seed-half-width", type=float, default=0.055)
    parser.add_argument("--sigma-accept", type=float, default=5.0e-6)
    parser.add_argument("--sigma-ratio-accept", type=float, default=5.0e-4)
    parser.add_argument("--mac-accept", type=float, default=0.25)
    parser.add_argument("--subspace-mac-accept", type=float, default=0.5)
    parser.add_argument("--cluster-gap-absolute", type=float, default=0.04)
    parser.add_argument("--cluster-gap-relative", type=float, default=3.0e-3)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--reuse-cache", dest="reuse_cache", action="store_true", default=True)
    parser.add_argument("--no-reuse-cache", dest="reuse_cache", action="store_false")
    parser.add_argument("--force-recompute", action="store_true")
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    ns = parser.parse_args(list(sys.argv[1:] if argv is None else argv))

    output_dir = repo_path(Path(ns.output_dir))
    beta_step = float(ns.beta_step)
    if bool(ns.smoke):
        output_dir = SMOKE_OUTPUT_DIR
        beta_step = 15.0
    args = Args(
        step3a_dir=repo_path(Path(ns.step3a_dir)),
        gateway_dir=repo_path(Path(ns.gateway_dir)),
        output_dir=output_dir,
        beta_min=float(ns.beta_min),
        beta_max=float(ns.beta_max),
        beta_step=beta_step,
        k_max=int(ns.k_max),
        n_spectrum_roots=int(ns.n_spectrum_roots),
        n_candidate_roots=int(ns.n_candidate_roots),
        verification_candidate_roots=int(ns.verification_candidate_roots),
        spectrum_method=str(ns.spectrum_method),
        initial_beta_step=float(ns.initial_beta_step),
        min_beta_step=float(ns.min_beta_step),
        max_beta_step=float(ns.max_beta_step),
        beta_growth_factor=float(ns.beta_growth_factor),
        beta_shrink_factor=float(ns.beta_shrink_factor),
        guard_scan_step=float(ns.guard_scan_step),
        seed_half_width=float(ns.seed_half_width),
        sigma_accept=float(ns.sigma_accept),
        sigma_ratio_accept=float(ns.sigma_ratio_accept),
        mac_accept=float(ns.mac_accept),
        subspace_mac_accept=float(ns.subspace_mac_accept),
        cluster_gap_absolute=float(ns.cluster_gap_absolute),
        cluster_gap_relative=float(ns.cluster_gap_relative),
        workers=int(ns.workers),
        reuse_cache=bool(ns.reuse_cache),
        force_recompute=bool(ns.force_recompute),
        plot_only=bool(ns.plot_only),
        smoke=bool(ns.smoke),
    )
    validate_args(args)
    return args


def validate_args(args: Args) -> None:
    if args.k_max != 10:
        raise ValueError("this diagnostic requires exactly --k-max 10")
    if args.n_spectrum_roots < args.k_max + 2:
        raise ValueError("--n-spectrum-roots must retain roots 11 and 12")
    if args.n_candidate_roots < args.n_spectrum_roots:
        raise ValueError("--n-candidate-roots must include all requested roots")
    if args.verification_candidate_roots <= args.n_candidate_roots:
        raise ValueError("--verification-candidate-roots must exceed --n-candidate-roots")
    if not (
        math.isfinite(args.beta_min)
        and math.isfinite(args.beta_max)
        and math.isfinite(args.beta_step)
        and 0.0 <= args.beta_min <= args.beta_max <= 90.0
        and args.beta_step > 0.0
    ):
        raise ValueError("beta grid must satisfy 0 <= min <= max <= 90 and step > 0")
    positives = (
        args.initial_beta_step,
        args.min_beta_step,
        args.max_beta_step,
        args.beta_growth_factor,
        args.beta_shrink_factor,
        args.guard_scan_step,
        args.seed_half_width,
        args.sigma_accept,
        args.sigma_ratio_accept,
        args.mac_accept,
        args.subspace_mac_accept,
        args.cluster_gap_absolute,
        args.cluster_gap_relative,
    )
    if any(not math.isfinite(value) or value <= 0.0 for value in positives):
        raise ValueError("continuation settings must be finite and positive")
    if not args.min_beta_step <= args.initial_beta_step <= args.max_beta_step:
        raise ValueError("beta continuation steps must satisfy min <= initial <= max")
    if args.beta_growth_factor <= 1.0 or not 0.0 < args.beta_shrink_factor < 1.0:
        raise ValueError("invalid beta growth/shrink factors")
    if args.plot_only and args.force_recompute:
        raise ValueError("--plot-only and --force-recompute are mutually exclusive")
    if args.workers < 1 or args.workers > 4:
        raise ValueError("--workers must lie in 1..4")


def load_step3a_cases(step3a_dir: Path) -> tuple[CaseSpec, ...]:
    manifest_path = step3a_dir / "step3a_manifest_resolved.csv"
    summary_path = step3a_dir / "step3a_case_summary.csv"
    source_path = manifest_path if manifest_path.exists() else summary_path
    if not source_path.exists():
        raise FileNotFoundError(
            "Step-3A full-precision manifest/summary is missing: "
            f"{manifest_path} or {summary_path}"
        )
    rows = {row["case_id"]: row for row in read_csv(source_path)}
    summary_rows = (
        {row["case_id"]: row for row in read_csv(summary_path)}
        if summary_path.exists()
        else {}
    )
    output: list[CaseSpec] = []
    for case_id, expected in EXPECTED_CASE_GEOMETRY.items():
        if case_id not in rows:
            raise ValueError(f"{case_id} is missing from {source_path}")
        row = rows[case_id]
        epsilon_text = row.get("epsilon", row.get("epsilon_0", ""))
        epsilon = float(epsilon_text)
        mu = float(row["mu"])
        eta = float(row["eta"])
        beta = float(row["beta_deg"])
        if not math.isfinite(epsilon) or epsilon <= 0.0:
            raise ValueError(f"invalid full-precision epsilon for {case_id}: {epsilon_text!r}")
        if (
            abs(mu - float(expected["mu"])) > 1.0e-14
            or abs(eta - float(expected["eta"])) > 1.0e-14
            or abs(beta - float(expected["counterexample_beta_deg"])) > 1.0e-12
        ):
            raise ValueError(f"Step-3A geometry mismatch for {case_id}")
        if summary_rows:
            status = summary_rows[case_id].get("final_scientific_status", "")
            if status != "confirmed_counterexample":
                raise ValueError(f"{case_id} is not confirmed in Step-3A: {status!r}")
        output.append(
            CaseSpec(
                case_id=case_id,
                epsilon_0=epsilon,
                mu=mu,
                eta=eta,
                counterexample_beta_deg=beta,
            )
        )
    return tuple(output)


def validate_gateway_provenance(gateway_dir: Path) -> None:
    summary_path = gateway_dir / "branch_k10_guard_summary.csv"
    report_path = gateway_dir / "eb_timo_branch_continuation_gateway_report.md"
    if not summary_path.exists() or not report_path.exists():
        raise FileNotFoundError("verified branch-continuation gateway outputs are missing")
    rows = read_csv(summary_path)
    if not rows or any(
        row.get("algorithm_version") != branch.BRANCH_CONTINUATION_ALGORITHM_VERSION
        for row in rows
    ):
        raise ValueError("gateway provenance is not uniformly branch_informed_continuation_v1")
    if "ready_for_targeted_step3" not in report_path.read_text(encoding="utf-8"):
        raise ValueError("branch-continuation gateway is not marked ready")


def beam_params_from_epsilon(epsilon_0: float) -> BeamParams:
    """Invert the project definition epsilon=sqrt(I/S)/L_base=r/(2L_base)."""

    l_base = 0.5 * CANONICAL_L_TOTAL_M
    radius = 2.0 * float(epsilon_0) * l_base
    params = BeamParams(
        E=CANONICAL_E_PA,
        rho=CANONICAL_RHO_KG_M3,
        r=radius,
        L_total=CANONICAL_L_TOTAL_M,
    )
    if not math.isclose(params.eps, float(epsilon_0), rel_tol=0.0, abs_tol=2.0e-16):
        raise RuntimeError("epsilon-to-radius conversion disagrees with BeamParams.eps")
    return params


def dimensional_frequencies(lambdas: Sequence[float], epsilon_0: float) -> np.ndarray:
    return lambdas_to_frequencies(
        np.asarray(lambdas, dtype=float), beam_params_from_epsilon(epsilon_0)
    )


def build_scale_rows(cases: Sequence[CaseSpec]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in cases:
        params = beam_params_from_epsilon(case.epsilon_0)
        scale = float(frequency_scale(params))
        eb_scale = float(dimensional_frequencies((1.0,), case.epsilon_0)[0])
        timo_scale = float(dimensional_frequencies((1.0,), case.epsilon_0)[0])
        rows.append(
            {
                "case_id": case.case_id,
                "epsilon_0": case.epsilon_0,
                "scaling_formula": FREQUENCY_FORMULA,
                "material_constants": f"E={params.E:.16e} Pa; rho={params.rho:.16e} kg/m^3",
                "length_radius_parameters": (
                    f"L_total={params.L_total:.16e} m; L_base={params.L_base:.16e} m; "
                    f"r0=2*epsilon_0*L_base={params.r:.16e} m"
                ),
                "E_Pa": params.E,
                "rho_kg_m3": params.rho,
                "L_total_m": params.L_total,
                "L_base_m": params.L_base,
                "r0_m": params.r,
                "S0_m2": params.S,
                "I0_m4": params.I,
                "common_scale_factor": scale,
                "frequency_kind": FREQUENCY_KIND,
                "unit": FREQUENCY_UNIT,
                "source_helper_file": SCALING_SOURCE,
                "source_helper": SCALING_HELPERS,
                "EB_scale_factor": eb_scale,
                "Timoshenko_scale_factor": timo_scale,
                "EB_Timoshenko_scale_equality_check": math.isclose(
                    eb_scale, timo_scale, rel_tol=0.0, abs_tol=0.0
                ),
            }
        )
    return rows


def beta_grid(beta_min: float, beta_max: float, beta_step: float) -> np.ndarray:
    count = int(math.floor((beta_max - beta_min) / beta_step + 1.0e-12)) + 1
    values = beta_min + beta_step * np.arange(max(count, 1), dtype=float)
    values = values[values <= beta_max + 1.0e-10]
    values = np.append(values, [beta_min, beta_max])
    for required in (45.0, 90.0):
        if beta_min - 1.0e-12 <= required <= beta_max + 1.0e-12:
            values = np.append(values, required)
    return np.unique(np.round(values, 12))


def continuation_settings(args: Args) -> branch.ContinuationSettings:
    return branch.ContinuationSettings(
        requested_roots=args.n_spectrum_roots,
        candidate_roots=args.n_candidate_roots,
        verification_candidate_roots=args.verification_candidate_roots,
        beta_initial_step_deg=args.initial_beta_step,
        beta_min_step_deg=args.min_beta_step,
        beta_max_step_deg=args.max_beta_step,
        beta_growth_factor=args.beta_growth_factor,
        beta_shrink_factor=args.beta_shrink_factor,
        guard_scan_step=args.guard_scan_step,
        seed_half_width=args.seed_half_width,
        sigma_accept=args.sigma_accept,
        sigma_ratio_accept=args.sigma_ratio_accept,
        mac_accept=args.mac_accept,
        subspace_mac_accept=args.subspace_mac_accept,
        cluster_gap_absolute=args.cluster_gap_absolute,
        cluster_gap_relative=args.cluster_gap_relative,
        run_global_guard=True,
        allow_strict_fallback=True,
    )


def minimum_gap(values: Sequence[float]) -> float:
    gaps = [float(right) - float(left) for left, right in zip(values, values[1:])]
    return min(gaps) if gaps else float("nan")


def adjacent_gap(values: Sequence[float], index_zero: int) -> float:
    candidates: list[float] = []
    if index_zero > 0:
        candidates.append(float(values[index_zero]) - float(values[index_zero - 1]))
    if index_zero + 1 < len(values):
        candidates.append(float(values[index_zero + 1]) - float(values[index_zero]))
    return min(candidates) if candidates else float("nan")


def _strict_recover(
    *,
    model: str,
    case: CaseSpec,
    beta_deg: float,
    previous: Sequence[branch.ContinuedBranch],
    settings: branch.ContinuationSettings,
    strict_cache: complete.GeneralSpectrumCache,
) -> tuple[tuple[branch.ContinuedBranch, ...], bool, str, str]:
    geometry = complete.Geometry(case.epsilon_0, beta_deg, case.mu, case.eta)
    strict = strict_cache.resolve(
        model,
        geometry,
        settings.strict_settings(),
        continuation_seeds=[item.Lambda for item in previous],
    )
    recovered = branch._branches_from_strict(model, beta_deg, previous, strict)
    direct_seed = any(
        "direct_full_matrix_SVD" in source
        for root in strict.primary.roots[:11]
        for source in root.detection_sources
    )
    warning = "" if strict.spectrum_status == "resolved_complete" else strict.exclusion_reason
    return recovered, direct_seed, strict.cache_status, warning


def _point_from_branch_result(
    case: CaseSpec,
    result: branch.BranchContinuationResult,
) -> SpectrumPoint:
    ordered = tuple(sorted(result.branches, key=lambda item: item.Lambda))
    warnings = tuple(
        item
        for item in (
            result.exclusion_reason,
            ";".join(result.guard.unresolved_intervals),
        )
        if item
    )
    return SpectrumPoint(
        case_id=case.case_id,
        model=result.model,
        beta_deg=result.geometry.beta_deg,
        branches=ordered,
        k10_guard_resolved=result.k10_guard_resolved,
        full12_resolved=result.full12_resolved,
        guard_status=result.guard.status,
        strict_fallback_used=result.operations.strict_fallback_runs > 0,
        strict_fallback_count=result.operations.strict_fallback_runs,
        cache_status=result.cache_status,
        warnings=warnings,
    )


def compute_case_model_sweep(
    case: CaseSpec,
    model: str,
    beta_values: Sequence[float],
    settings: branch.ContinuationSettings,
    branch_cache: branch.BranchContinuationCache,
    strict_cache: complete.GeneralSpectrumCache,
) -> tuple[SpectrumPoint, ...]:
    """Continue existing branch records sequentially to every requested beta node."""

    targets = sorted({0.0, *(float(value) for value in beta_values)})
    base_geometry = complete.Geometry(case.epsilon_0, 0.0, case.mu, case.eta)
    base = branch_cache.resolve(model, base_geometry, settings)
    current = tuple(sorted(base.branches, key=lambda item: item.Lambda))
    base_point = _point_from_branch_result(case, base)
    if not base_point.k10_guard_resolved:
        base_operations = branch.BranchOperationCounts()
        recovered, direct_seed, strict_status, strict_warning = _strict_recover(
            model=model,
            case=case,
            beta_deg=0.0,
            previous=current,
            settings=settings,
            strict_cache=strict_cache,
        )
        base_warnings = list(base_point.warnings)
        if strict_warning:
            base_warnings.append(strict_warning)
        if len(recovered) == len(current):
            current = tuple(sorted(recovered, key=lambda item: item.Lambda))
            base_guard = branch._guard(model, base_geometry, current, settings, base_operations)
            first11 = current[:11]
            base_k10 = bool(
                len(first11) == 11
                and all(
                    item.sigma_1 <= settings.sigma_accept
                    and item.sigma_ratio <= settings.sigma_ratio_accept
                    and item.refinement_status == branch.SEED_REFINED_TO_NEW_ROOT
                    for item in first11
                )
                and base_guard.passed
                and not direct_seed
            )
            first12 = current[:12]
            base_full12 = bool(
                base_k10
                and len(first12) == 12
                and all(
                    item.sigma_1 <= settings.sigma_accept
                    and item.sigma_ratio <= settings.sigma_ratio_accept
                    for item in first12
                )
            )
            if not base_guard.passed:
                base_warnings.extend(base_guard.unresolved_intervals)
            if direct_seed:
                base_warnings.append("strict_direct_seed_acceptance_present")
            base_point = SpectrumPoint(
                case_id=case.case_id,
                model=model,
                beta_deg=0.0,
                branches=current,
                k10_guard_resolved=base_k10,
                full12_resolved=base_full12,
                guard_status=base_guard.status,
                strict_fallback_used=True,
                strict_fallback_count=base_point.strict_fallback_count + 1,
                cache_status=f"strict_{strict_status}",
                warnings=tuple(dict.fromkeys(item for item in base_warnings if item)),
            )
        else:
            base_warnings.append("strict_fallback_branch_count_mismatch")
            base_point = replace(
                base_point,
                strict_fallback_used=True,
                strict_fallback_count=base_point.strict_fallback_count + 1,
                cache_status=f"strict_{strict_status}",
                warnings=tuple(dict.fromkeys(base_warnings)),
            )
    points: dict[float, SpectrumPoint] = {0.0: base_point}
    earlier: tuple[branch.ContinuedBranch, ...] = ()
    current_beta = 0.0

    for grid_target in targets[1:]:
        interval_warnings: list[str] = []
        interval_fallback = False
        interval_fallback_count = 0
        cache_status = "sequential_continuation"
        interval_steps: list[branch.ContinuationStep] = []
        reached = current_beta >= grid_target - 1.0e-12
        step_size = min(settings.beta_initial_step_deg, max(grid_target - current_beta, 0.0))
        operations = branch.BranchOperationCounts()
        direct_seed_acceptance = False
        while current_beta < grid_target - 1.0e-12:
            step_target = min(grid_target, current_beta + step_size)
            operations.beta_steps_attempted += 1
            attempted = branch._attempt_step(
                model,
                complete.Geometry(case.epsilon_0, grid_target, case.mu, case.eta),
                current_beta,
                step_target,
                current,
                {item.branch_id: item for item in earlier},
                settings,
                operations,
            )
            interval_steps.append(attempted)
            if attempted.accepted:
                earlier = current
                current = tuple(sorted(attempted.branches, key=lambda item: item.Lambda))
                current_beta = step_target
                operations.beta_steps_accepted += 1
                step_size = min(
                    settings.beta_max_step_deg,
                    max(settings.beta_min_step_deg, step_size * settings.beta_growth_factor),
                )
                continue
            if step_size * settings.beta_shrink_factor >= settings.beta_min_step_deg:
                step_size *= settings.beta_shrink_factor
                operations.beta_step_reductions += 1
                continue
            interval_fallback = True
            interval_fallback_count += 1
            recovered, direct_seed, strict_status, strict_warning = _strict_recover(
                model=model,
                case=case,
                beta_deg=grid_target,
                previous=current,
                settings=settings,
                strict_cache=strict_cache,
            )
            cache_status = f"strict_{strict_status}"
            direct_seed_acceptance = direct_seed_acceptance or direct_seed
            if strict_warning:
                interval_warnings.append(strict_warning)
            if len(recovered) != len(current):
                interval_warnings.append("strict_fallback_branch_count_mismatch")
                break
            earlier = current
            current = tuple(sorted(recovered, key=lambda item: item.Lambda))
            current_beta = grid_target
            operations.beta_steps_accepted += 1
            break
        reached = current_beta >= grid_target - 1.0e-12

        ordered = tuple(sorted(current, key=lambda item: item.Lambda))
        geometry = complete.Geometry(case.epsilon_0, grid_target, case.mu, case.eta)
        guard = (
            branch._guard(model, geometry, ordered, settings, operations)
            if reached
            else branch.GlobalGuardResult(
                0.0,
                (),
                (),
                (),
                ("continuation_target_not_reached",),
                False,
                "guard_unresolved",
            )
        )
        if reached and not guard.passed:
            interval_fallback = True
            interval_fallback_count += 1
            recovered, direct_seed, strict_status, strict_warning = _strict_recover(
                model=model,
                case=case,
                beta_deg=grid_target,
                previous=ordered,
                settings=settings,
                strict_cache=strict_cache,
            )
            cache_status = f"strict_{strict_status}"
            direct_seed_acceptance = direct_seed_acceptance or direct_seed
            if strict_warning:
                interval_warnings.append(strict_warning)
            if len(recovered) == len(ordered):
                earlier = ordered
                current = tuple(sorted(recovered, key=lambda item: item.Lambda))
                current_beta = grid_target
                ordered = current
                guard = branch._guard(model, geometry, ordered, settings, operations)
            else:
                interval_warnings.append("guard_strict_fallback_branch_count_mismatch")

        first11 = ordered[:11]
        final_ids = {item.branch_id for item in first11}
        clusters_resolved = all(
            cluster.resolved
            for step in interval_steps
            if step.accepted
            for cluster in step.clusters
            if final_ids.intersection(cluster.branch_ids)
        )
        root_quality = all(
            item.sigma_1 <= settings.sigma_accept
            and item.sigma_ratio <= settings.sigma_ratio_accept
            and item.refinement_status == branch.SEED_REFINED_TO_NEW_ROOT
            for item in first11
        )
        k10 = bool(
            reached
            and len(first11) == 11
            and root_quality
            and clusters_resolved
            and guard.passed
            and not direct_seed_acceptance
        )
        first12 = ordered[:12]
        full12 = bool(
            k10
            and len(first12) == 12
            and all(
                item.sigma_1 <= settings.sigma_accept
                and item.sigma_ratio <= settings.sigma_ratio_accept
                for item in first12
            )
        )
        if operations.beta_step_reductions:
            interval_warnings.append(
                f"adaptive_step_reductions={operations.beta_step_reductions}"
            )
        if not reached:
            interval_warnings.append("continuation_target_not_reached")
        if not guard.passed:
            interval_warnings.extend(guard.unresolved_intervals)
        if not clusters_resolved:
            interval_warnings.append("cluster_continuation_unresolved")
        if direct_seed_acceptance:
            interval_warnings.append("strict_direct_seed_acceptance_present")
        points[grid_target] = SpectrumPoint(
            case_id=case.case_id,
            model=model,
            beta_deg=grid_target,
            branches=ordered,
            k10_guard_resolved=k10,
            full12_resolved=full12,
            guard_status=guard.status,
            strict_fallback_used=interval_fallback,
            strict_fallback_count=interval_fallback_count,
            cache_status=cache_status,
            warnings=tuple(dict.fromkeys(item for item in interval_warnings if item)),
        )
        print(
            f"{case.case_id} {model} beta={grid_target:g}: "
            f"{'K10_guard_resolved' if k10 else 'unresolved'}"
        )
    return tuple(points[float(value)] for value in beta_values)


def _compute_sweep_process(
    case: CaseSpec,
    model: str,
    beta_values: tuple[float, ...],
    settings: branch.ContinuationSettings,
    gateway_dir: Path,
    reuse_cache: bool,
    force_recompute: bool,
) -> tuple[str, str, tuple[SpectrumPoint, ...]]:
    """Process-safe wrapper; the four case/model sweeps share no mutable state."""

    branch_cache = branch.BranchContinuationCache(
        gateway_dir / "cache",
        reuse_cache=reuse_cache,
        force_recompute=force_recompute,
        verification_scope="primary",
    )
    strict_cache = complete.GeneralSpectrumCache(
        gateway_dir / "cache" / "counterexample_dimensional_frequency_beta_strict",
        reuse_cache=reuse_cache,
        force_recompute=force_recompute,
    )
    values = compute_case_model_sweep(
        case,
        model,
        beta_values,
        settings,
        branch_cache,
        strict_cache,
    )
    return case.case_id, model, values


def build_frequency_rows(
    cases: Sequence[CaseSpec],
    points: Mapping[tuple[str, str], Sequence[SpectrumPoint]],
    k_max: int,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in cases:
        for model in MODELS:
            for point in points[(case.case_id, model)]:
                values = point.values
                frequencies = dimensional_frequencies(values[:k_max], case.epsilon_0)
                for index in range(1, k_max + 1):
                    resolved = point.k10_guard_resolved and index <= len(point.branches)
                    item = point.branches[index - 1] if resolved else None
                    rows.append(
                        {
                            "case_id": case.case_id,
                            "epsilon_0": case.epsilon_0,
                            "mu": case.mu,
                            "eta": case.eta,
                            "beta_deg": point.beta_deg,
                            "sorted_index": index,
                            "model": model,
                            "branch_id": item.branch_id if item else "",
                            "parent_family": item.parent_family if item else "",
                            "Lambda": item.Lambda if item else float("nan"),
                            "dimensional_frequency": (
                                float(frequencies[index - 1])
                                if item is not None and index <= len(frequencies)
                                else float("nan")
                            ),
                            "frequency_kind": FREQUENCY_KIND,
                            "frequency_unit": FREQUENCY_UNIT,
                            "K10_guard_resolved": point.k10_guard_resolved,
                            "full12_resolved": point.full12_resolved,
                            "cluster_id": item.cluster_id if item else "",
                            "adjacent_gap": (
                                adjacent_gap(values[:11], index - 1)
                                if item is not None
                                else float("nan")
                            ),
                            "strict_fallback_used": point.strict_fallback_used,
                            "warnings": ";".join(point.warnings),
                        }
                    )
    return sorted(
        rows,
        key=lambda row: (
            str(row["case_id"]),
            float(row["beta_deg"]),
            int(row["sorted_index"]),
            str(row["model"]),
        ),
    )


def build_quality_rows(
    cases: Sequence[CaseSpec],
    points: Mapping[tuple[str, str], Sequence[SpectrumPoint]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for case in cases:
        for model in MODELS:
            for point in points[(case.case_id, model)]:
                values = point.values
                clusters = {item.cluster_id for item in point.branches[:11] if item.cluster_id}
                rows.append(
                    {
                        "case_id": case.case_id,
                        "epsilon_0": case.epsilon_0,
                        "mu": case.mu,
                        "eta": case.eta,
                        "model": model,
                        "beta_deg": point.beta_deg,
                        "root_count": len(values),
                        "K10_status": (
                            "K10_guard_resolved" if point.k10_guard_resolved else "unresolved"
                        ),
                        "K10_guard_resolved": point.k10_guard_resolved,
                        "root11_status": (
                            "resolved_guard" if point.k10_guard_resolved and len(values) >= 11 else "unresolved"
                        ),
                        "root11": values[10] if len(values) >= 11 else float("nan"),
                        "full12_status": (
                            "resolved" if point.full12_resolved else "optional_not_resolved"
                        ),
                        "full12_resolved": point.full12_resolved,
                        "minimum_gap": minimum_gap(values[:11]),
                        "cluster_count": len(clusters),
                        "fallback": (
                            "strict_fallback_used" if point.strict_fallback_used else "not_used"
                        ),
                        "strict_fallback_used": point.strict_fallback_used,
                        "strict_fallback_count": point.strict_fallback_count,
                        "cache_status": point.cache_status,
                        "warnings": ";".join(point.warnings),
                    }
                )
    return sorted(rows, key=lambda row: (str(row["case_id"]), str(row["model"]), float(row["beta_deg"])))


def create_case_figure(
    rows: Sequence[Mapping[str, object]],
    case_id: str,
    *,
    k_max: int = 10,
) -> tuple[plt.Figure, plt.Axes]:
    selected = [row for row in rows if str(row["case_id"]) == case_id]
    if not selected:
        raise ValueError(f"no dimensional-frequency rows for {case_id}")
    fig, ax = plt.subplots(figsize=(8.0, 5.2), constrained_layout=True)
    color_map = plt.get_cmap("tab10")
    for sorted_index in range(1, k_max + 1):
        color = color_map(sorted_index - 1)
        for model, linestyle in ((MODEL_EB, "--"), (MODEL_TIMO, "-")):
            series = sorted(
                (
                    row
                    for row in selected
                    if int(row["sorted_index"]) == sorted_index
                    and str(row["model"]) == model
                ),
                key=lambda row: float(row["beta_deg"]),
            )
            ax.plot(
                [float(row["beta_deg"]) for row in series],
                [float(row["dimensional_frequency"]) for row in series],
                color=color,
                linestyle=linestyle,
                linewidth=1.45,
            )
    ax.set_xlabel(r"$\beta,\ ^\circ$")
    ax.set_ylabel(r"$f,\ \mathrm{Hz}$")
    ax.margins(x=0.0)
    return fig, ax


def create_plots(
    rows: Sequence[Mapping[str, object]],
    output_dir: Path,
    *,
    k_max: int = 10,
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for case_id in ("S3_12", "S3_14"):
        fig, _ax = create_case_figure(rows, case_id, k_max=k_max)
        path = output_dir / PDF_NAMES[case_id]
        fig.savefig(
            path,
            format="pdf",
            metadata={
                "Creator": "CoupledBeams diagnostic",
                "CreationDate": None,
                "ModDate": None,
            },
        )
        plt.close(fig)
        paths.append(path)
    graphical = sorted(
        path.name
        for path in output_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".pdf", ".png", ".svg"}
    )
    expected = sorted(PDF_NAMES.values())
    if graphical != expected:
        raise RuntimeError(
            f"graphical output contract requires exactly {expected}, found {graphical}"
        )
    return paths[0], paths[1]


def validate_plot_csv_rows(
    rows: Sequence[Mapping[str, object]],
    *,
    k_max: int,
) -> None:
    if not rows:
        raise ValueError("dimensional-frequency CSV is empty")
    case_ids = {str(row["case_id"]) for row in rows}
    if case_ids != set(PDF_NAMES):
        raise ValueError(f"plot CSV must contain exactly {sorted(PDF_NAMES)}, got {sorted(case_ids)}")
    if {int(row["sorted_index"]) for row in rows} != set(range(1, k_max + 1)):
        raise ValueError("plot CSV does not contain exactly sorted indices 1..10")
    if {str(row["model"]) for row in rows} != set(MODELS):
        raise ValueError("plot CSV must contain both EB and Timoshenko")
    if any(str(row["frequency_kind"]) != FREQUENCY_KIND for row in rows):
        raise ValueError("plot CSV frequency kind mismatch")
    if any(str(row["frequency_unit"]) != FREQUENCY_UNIT for row in rows):
        raise ValueError("plot CSV frequency unit mismatch")


def write_report(
    path: Path,
    cases: Sequence[CaseSpec],
    beta_values: Sequence[float],
    scale_rows: Sequence[Mapping[str, object]],
    quality_rows: Sequence[Mapping[str, object]],
) -> Path:
    k10_count = sum(bool(row["K10_guard_resolved"]) for row in quality_rows)
    full12_count = sum(bool(row["full12_resolved"]) for row in quality_rows)
    fallback_point_count = sum(bool(row["strict_fallback_used"]) for row in quality_rows)
    fallback_count = sum(int(row["strict_fallback_count"]) for row in quality_rows)
    unresolved_count = len(quality_rows) - k10_count
    lines = [
        "# Counterexample dimensional frequency versus beta",
        "",
        "## Dimensional scaling",
        "",
        f"- Formula: `{FREQUENCY_FORMULA}`.",
        f"- Quantity and unit: `{FREQUENCY_KIND}`, `{FREQUENCY_UNIT}`.",
        f"- Shared helper: `{SCALING_HELPERS}` from `{SCALING_SOURCE}`.",
        f"- Material convention: `E={CANONICAL_E_PA:.16e} Pa`, `rho={CANONICAL_RHO_KG_M3:.16e} kg/m^3`.",
        f"- Length convention: `L_total={CANONICAL_L_TOTAL_M:.16e} m`, `L_base={0.5 * CANONICAL_L_TOTAL_M:.16e} m`; `r0=2*epsilon_0*L_base`.",
        "- The EB and Timoshenko curves use the same dimensional scale for each case.",
        "",
        "## Cases",
        "",
    ]
    scale_by_case = {str(row["case_id"]): row for row in scale_rows}
    for case in cases:
        scale = scale_by_case[case.case_id]
        lines.append(
            f"- `{case.case_id}`: `epsilon_0={case.epsilon_0:.16e}`, `mu={case.mu:.16e}`, "
            f"`eta={case.eta:.16e}`, `r0={float(scale['r0_m']):.16e} m`, "
            f"`f/Lambda^2={float(scale['common_scale_factor']):.16e} Hz`."
        )
    lines.extend(
        [
            "",
            "## Beta grid and spectrum quality",
            "",
            f"- Grid: `{float(beta_values[0]):g}..{float(beta_values[-1]):g} deg`, "
            f"`{len(beta_values)}` points; exact `beta=45 deg` and `beta=90 deg` are included.",
            f"- K10/root-11 resolved model-points: `{k10_count}/{len(quality_rows)}`.",
            f"- Optional full-12 resolved model-points: `{full12_count}/{len(quality_rows)}`.",
            f"- Strict fallback invocations: `{fallback_count}` across `{fallback_point_count}` model-points.",
            f"- Unresolved model-points: `{unresolved_count}`; unresolved target values remain NaN and are not connected across gaps.",
            "",
            "## Outputs and figure contract",
            "",
            f"- `{PDF_NAMES['S3_12']}`",
            f"- `{PDF_NAMES['S3_14']}`",
            f"- `{MAIN_CSV}`",
            f"- `{SCALE_CSV}`",
            f"- `{QUALITY_CSV}`",
            "- Each PDF contains sorted positions 1..10 only: Euler--Bernoulli dashed, Timoshenko solid, and a shared color per sorted index.",
            "- Figures have no title, suptitle, legend, annotation, marker, grid, inset, or case/parameter text; only axes, ticks, axis labels, and curves are present.",
            "- Root 11 is retained only as the right K10 gap guard; root 12 is diagnostic and neither is plotted.",
            "",
            "## Scope",
            "",
            "No physical formula, matrix entry, unknown ordering, shared/global solver default, tolerance, Timoshenko shear coefficient, FEM workflow, or article workspace was changed.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def plot_only(args: Args) -> dict[str, object]:
    csv_path = args.output_dir / MAIN_CSV
    if not csv_path.exists():
        raise FileNotFoundError(f"--plot-only requires {csv_path}")
    rows = read_csv(csv_path)
    validate_plot_csv_rows(rows, k_max=args.k_max)
    pdfs = create_plots(rows, args.output_dir, k_max=args.k_max)
    return {
        "args": args,
        "plot_only": True,
        "root_calculations": 0,
        "pdfs": pdfs,
    }


def execute(
    args: Args,
    *,
    sweep_solver: Callable[
        [
            CaseSpec,
            str,
            Sequence[float],
            branch.ContinuationSettings,
            branch.BranchContinuationCache,
            complete.GeneralSpectrumCache,
        ],
        tuple[SpectrumPoint, ...],
    ] = compute_case_model_sweep,
) -> dict[str, object]:
    if args.plot_only:
        return plot_only(args)
    cases = load_step3a_cases(args.step3a_dir)
    validate_gateway_provenance(args.gateway_dir)
    beta_values = beta_grid(args.beta_min, args.beta_max, args.beta_step)
    if not any(abs(value - 45.0) <= 1.0e-12 for value in beta_values):
        raise ValueError("the selected beta range must include beta=45 deg")
    if not any(abs(value - 90.0) <= 1.0e-12 for value in beta_values):
        raise ValueError("the selected beta range must include beta=90 deg")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    existing_extra_graphics = [
        path
        for path in args.output_dir.iterdir()
        if path.is_file()
        and path.suffix.lower() in {".pdf", ".png", ".svg"}
        and path.name not in PDF_NAMES.values()
    ]
    if existing_extra_graphics:
        raise RuntimeError(
            "refusing to violate the exactly-two-graphics contract; unexpected files: "
            + ", ".join(path.name for path in existing_extra_graphics)
        )

    settings = continuation_settings(args)
    settings.validate(k_max=args.k_max)
    points: dict[tuple[str, str], tuple[SpectrumPoint, ...]] = {}
    root_calculations = 0
    tasks = [(case, model) for case in cases for model in MODELS]
    if args.workers > 1 and sweep_solver is compute_case_model_sweep:
        with ProcessPoolExecutor(max_workers=min(args.workers, len(tasks))) as executor:
            futures = [
                executor.submit(
                    _compute_sweep_process,
                    case,
                    model,
                    tuple(float(value) for value in beta_values),
                    settings,
                    args.gateway_dir,
                    args.reuse_cache,
                    args.force_recompute,
                )
                for case, model in tasks
            ]
            for future in as_completed(futures):
                case_id, model, values = future.result()
                points[(case_id, model)] = values
                root_calculations += 1
    else:
        branch_cache = branch.BranchContinuationCache(
            args.gateway_dir / "cache",
            reuse_cache=args.reuse_cache,
            force_recompute=args.force_recompute,
            verification_scope="primary",
        )
        strict_cache = complete.GeneralSpectrumCache(
            args.gateway_dir / "cache" / "counterexample_dimensional_frequency_beta_strict",
            reuse_cache=args.reuse_cache,
            force_recompute=args.force_recompute,
        )
        for case, model in tasks:
            points[(case.case_id, model)] = sweep_solver(
                case,
                model,
                beta_values,
                settings,
                branch_cache,
                strict_cache,
            )
            root_calculations += 1

    frequency_rows = build_frequency_rows(cases, points, args.k_max)
    quality_rows = build_quality_rows(cases, points)
    scale_rows = build_scale_rows(cases)
    write_csv(args.output_dir / MAIN_CSV, frequency_rows, MAIN_FIELDS)
    write_csv(args.output_dir / SCALE_CSV, scale_rows, SCALE_FIELDS)
    write_csv(args.output_dir / QUALITY_CSV, quality_rows, QUALITY_FIELDS)
    pdfs = create_plots(frequency_rows, args.output_dir, k_max=args.k_max)
    report = write_report(
        args.output_dir / REPORT_MD,
        cases,
        beta_values,
        scale_rows,
        quality_rows,
    )
    return {
        "args": args,
        "cases": cases,
        "beta_values": beta_values,
        "frequency_rows": frequency_rows,
        "quality_rows": quality_rows,
        "scale_rows": scale_rows,
        "root_calculations": root_calculations,
        "pdfs": pdfs,
        "report": report,
    }


def main(
    argv: Sequence[str] | None = None,
    *,
    sweep_solver: Callable[
        [
            CaseSpec,
            str,
            Sequence[float],
            branch.ContinuationSettings,
            branch.BranchContinuationCache,
            complete.GeneralSpectrumCache,
        ],
        tuple[SpectrumPoint, ...],
    ] = compute_case_model_sweep,
) -> dict[str, object]:
    return execute(parse_args(argv), sweep_solver=sweep_solver)


if __name__ == "__main__":
    main()
