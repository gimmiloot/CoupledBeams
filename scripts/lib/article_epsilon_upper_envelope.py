"""Shared value-only helpers for the article epsilon upper-envelope workflow.

This module deliberately contains no characteristic equation or root solver.
It builds the finite manifest, defines cache identities, and implements the
sorted-frequency prefix statistics consumed by the runner and CSV-only
postprocessor.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import csv
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from scripts.lib import branch_informed_spectrum_continuation as branch
from scripts.lib import general_spectrum_completeness as complete
from src.my_project.analytic.formulas_thickness_mismatch import (
    thickness_mismatch_factors,
    thickness_to_length_ratios,
)


WORKFLOW_VERSION = "article_epsilon_upper_envelope_v8_solver_readiness_v2"
K_MAX = 10
ROOT_GUARD_INDEX = 11
REQUESTED_ROOTS = 12
DELTA_TOL = 0.10
THIN_DIAGNOSTIC_LIMIT = 0.10

BASE_BETA_DEG = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0)
BASE_MU = (0.0, 0.3, 0.5, 0.7, 0.9)
BASE_ETA = (-0.5, -0.25, 0.0, 0.25, 0.5)
EPSILON_VALUES = (0.010, 0.015, 0.020, 0.025, 0.030, 0.040, 0.050, 0.060)

LOW_ANGLE_BETA_DEG = (5.0, 10.0)
LOW_ANGLE_MU = (0.0, 0.5, 0.9)
LOW_ANGLE_ETA = (-0.5, 0.0, 0.5)

S3_14_SWEEP_GEOMETRY = (45.0, 0.5, -0.1)
REGRESSION_POINTS = (
    (
        "S3_12",
        0.029408510742187498,
        90.0,
        0.7,
        0.0,
        4,
        0.11739469909177262,
    ),
    (
        "S3_14",
        0.024798906738281248,
        45.0,
        0.5,
        -0.1,
        4,
        0.10050934854803291,
    ),
)

SMOKE_BETA_DEG = (0.0, 90.0)
SMOKE_MU = (0.0, 0.7)
SMOKE_ETA = (0.0, -0.5)
SMOKE_EPSILON = (0.02, 0.05)

MAIN_OUTPUT_DIR = Path("results") / "article_epsilon_upper_envelope"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "article_epsilon_upper_envelope"

MANIFEST_FIELDS = (
    "case_id",
    "case_identity",
    "grid_group",
    "regression_label",
    "claim_eligible",
    "epsilon_0",
    "epsilon_0_hex",
    "beta_deg",
    "beta_deg_hex",
    "mu",
    "mu_hex",
    "eta",
    "eta_hex",
    "tau1",
    "tau2",
    "mass_factor",
    "mass_residual",
    "l1_over_l",
    "l2_over_l",
    "r1_over_r0",
    "r2_over_r0",
    "s1",
    "s2",
    "s_max",
    "thin_0p1_flag",
)


@dataclass(frozen=True)
class GridPoint:
    case_id: str
    case_identity: str
    grid_group: str
    regression_label: str
    claim_eligible: bool
    epsilon_0: float
    beta_deg: float
    mu: float
    eta: float
    tau1: float
    tau2: float
    mass_factor: float
    mass_residual: float
    l1_over_l: float
    l2_over_l: float
    r1_over_r0: float
    r2_over_r0: float
    s1: float
    s2: float
    s_max: float
    thin_0p1_flag: bool

    @property
    def geometry(self) -> complete.Geometry:
        return complete.Geometry(self.epsilon_0, self.beta_deg, self.mu, self.eta)

    def manifest_row(self) -> dict[str, object]:
        row = asdict(self)
        row.update(
            {
                "epsilon_0_hex": float(self.epsilon_0).hex(),
                "beta_deg_hex": float(self.beta_deg).hex(),
                "mu_hex": float(self.mu).hex(),
                "eta_hex": float(self.eta).hex(),
            }
        )
        return row


def geometry_identity(epsilon_0: float, beta_deg: float, mu: float, eta: float) -> str:
    """Return a deterministic full-binary-precision geometry identity."""

    payload = {
        "epsilon_0": float(epsilon_0).hex(),
        "beta_deg": float(beta_deg).hex(),
        "mu": float(mu).hex(),
        "eta": float(eta).hex(),
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _point(
    grid_group: str,
    epsilon_0: float,
    beta_deg: float,
    mu: float,
    eta: float,
    *,
    regression_label: str = "",
) -> GridPoint:
    factors = thickness_mismatch_factors(float(mu), float(eta))
    s1, s2 = thickness_to_length_ratios(float(epsilon_0), float(mu), float(eta))
    identity = geometry_identity(epsilon_0, beta_deg, mu, eta)
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:20]
    mass_residual = float(factors.mass_factor - 2.0)
    if abs(mass_residual) > 5.0e-13:
        raise ValueError(f"mass preservation failed for {identity}: residual={mass_residual:.17g}")
    return GridPoint(
        case_id=f"AUE_{digest}",
        case_identity=identity,
        grid_group=str(grid_group),
        regression_label=str(regression_label),
        claim_eligible=grid_group != "regression",
        epsilon_0=float(epsilon_0),
        beta_deg=float(beta_deg),
        mu=float(mu),
        eta=float(eta),
        tau1=float(factors.tau1),
        tau2=float(factors.tau2),
        mass_factor=float(factors.mass_factor),
        mass_residual=mass_residual,
        l1_over_l=1.0 - float(mu),
        l2_over_l=1.0 + float(mu),
        r1_over_r0=float(factors.tau1),
        r2_over_r0=float(factors.tau2),
        s1=float(s1),
        s2=float(s2),
        s_max=max(float(s1), float(s2)),
        thin_0p1_flag=max(float(s1), float(s2)) <= THIN_DIAGNOSTIC_LIMIT,
    )


def _deduplicate(points: Iterable[GridPoint]) -> list[GridPoint]:
    by_identity: dict[str, GridPoint] = {}
    for point in points:
        previous = by_identity.get(point.case_identity)
        if previous is not None:
            raise ValueError(
                "unexpected full-precision geometry collision: "
                f"{previous.grid_group}/{previous.regression_label} and "
                f"{point.grid_group}/{point.regression_label}: {point.case_identity}"
            )
        by_identity[point.case_identity] = point
    return sorted(
        by_identity.values(),
        key=lambda item: (
            {"base": 0, "low_angle": 1, "s3_14_sweep": 2, "regression": 3, "smoke": 4}.get(
                item.grid_group, 99
            ),
            item.beta_deg,
            item.mu,
            item.eta,
            item.epsilon_0,
            item.regression_label,
        ),
    )


def build_manifest(*, smoke: bool = False) -> list[GridPoint]:
    points: list[GridPoint] = []
    if smoke:
        for beta_deg in SMOKE_BETA_DEG:
            for mu in SMOKE_MU:
                for eta in SMOKE_ETA:
                    for epsilon_0 in SMOKE_EPSILON:
                        points.append(_point("smoke", epsilon_0, beta_deg, mu, eta))
        result = _deduplicate(points)
        if len(result) != 16:
            raise ValueError(f"smoke manifest must contain 16 unique points, found {len(result)}")
        return result

    for beta_deg in BASE_BETA_DEG:
        for mu in BASE_MU:
            for eta in BASE_ETA:
                for epsilon_0 in EPSILON_VALUES:
                    points.append(_point("base", epsilon_0, beta_deg, mu, eta))
    for beta_deg in LOW_ANGLE_BETA_DEG:
        for mu in LOW_ANGLE_MU:
            for eta in LOW_ANGLE_ETA:
                for epsilon_0 in EPSILON_VALUES:
                    points.append(_point("low_angle", epsilon_0, beta_deg, mu, eta))
    beta_deg, mu, eta = S3_14_SWEEP_GEOMETRY
    for epsilon_0 in EPSILON_VALUES:
        points.append(_point("s3_14_sweep", epsilon_0, beta_deg, mu, eta))
    for label, epsilon_0, beta_deg, mu, eta, _expected_n, _expected_delta in REGRESSION_POINTS:
        points.append(
            _point(
                "regression",
                epsilon_0,
                beta_deg,
                mu,
                eta,
                regression_label=label,
            )
        )
    result = _deduplicate(points)
    counts = group_counts(result)
    expected = {"base": 1400, "low_angle": 144, "s3_14_sweep": 8, "regression": 2}
    if counts != expected or len(result) != 1554:
        raise ValueError(f"manifest contract mismatch: counts={counts}, total={len(result)}")
    return result


def group_counts(points: Sequence[GridPoint]) -> dict[str, int]:
    result: dict[str, int] = {}
    for point in points:
        result[point.grid_group] = result.get(point.grid_group, 0) + 1
    return result


def select_points(
    points: Sequence[GridPoint],
    *,
    base_only: bool = False,
    low_angle_only: bool = False,
    regressions_only: bool = False,
) -> list[GridPoint]:
    selected_groups: set[str] | None = None
    if base_only:
        selected_groups = {"base"}
    elif low_angle_only:
        selected_groups = {"low_angle"}
    elif regressions_only:
        selected_groups = {"regression"}
    if selected_groups is None:
        return list(points)
    return [point for point in points if point.grid_group in selected_groups]


def primary_settings() -> complete.SearchSettings:
    return complete.SearchSettings(
        requested_roots=REQUESTED_ROOTS,
        candidate_roots=20,
        verification_candidate_roots=24,
    )


def strict_settings() -> branch.ContinuationSettings:
    return branch.ContinuationSettings(
        requested_roots=REQUESTED_ROOTS,
        # K=10 needs roots 1..11, while root 12 is an optional right margin.
        # Asking the beta=0 Timoshenko parent for 20 roots can cross into the
        # high-frequency basis-regime change before the relevant guard.  Keep
        # the branch candidate target at the complete first-12 spectrum and
        # retain a larger independent verification target.
        candidate_roots=REQUESTED_ROOTS,
        verification_candidate_roots=16,
        # The declared grid reaches mu=0.9 and |eta|=0.5.  The unchanged
        # beta=0 block solver needs a wider technical search interval there
        # to collect its 20 candidate roots before continuation starts.
        lambda_max=80.0,
    )


def solver_configuration() -> dict[str, object]:
    primary = primary_settings()
    strict = strict_settings()
    return {
        "workflow_version": WORKFLOW_VERSION,
        "general_algorithm_version": complete.GENERAL_SPECTRUM_ALGORITHM_VERSION,
        "eb_matrix_evaluator_version": complete.EB_MATRIX_EVALUATOR_VERSION,
        "branch_algorithm_version": branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
        "timoshenko_basis_evaluator_version": complete.TIMO.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
        "primary_settings": asdict(primary),
        "strict_settings": asdict(strict),
        "k_max": K_MAX,
        "root_guard_index": ROOT_GUARD_INDEX,
        "delta_tol": DELTA_TOL,
        "sorted_frequency_definition": True,
    }


def case_cache_identity(point: GridPoint, configuration: Mapping[str, object] | None = None) -> dict[str, object]:
    return {
        "case_identity": point.case_identity,
        "solver_configuration": dict(configuration or solver_configuration()),
    }


def cache_digest(configuration: Mapping[str, object] | None = None) -> str:
    encoded = json.dumps(dict(configuration or solver_configuration()), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()[:16]


def squared_frequency_delta(lambda_eb: float, lambda_timo: float) -> float:
    denominator = float(lambda_timo) ** 2
    if not math.isfinite(denominator) or denominator <= 0.0:
        return float("nan")
    return abs(float(lambda_eb) ** 2 - denominator) / denominator


def true_safe_prefix(deltas: Sequence[float], *, threshold: float = DELTA_TOL, k_max: int = K_MAX) -> int:
    if len(deltas) < int(k_max):
        raise ValueError(f"N_true requires at least {k_max} ordered deltas")
    safe = 0
    for value in deltas[: int(k_max)]:
        if not math.isfinite(float(value)) or float(value) > float(threshold):
            break
        safe += 1
    return safe


def suffix_max(values: Sequence[int | float]) -> list[int | float]:
    result: list[int | float] = [float("nan")] * len(values)
    running: int | float = float("nan")
    for index in range(len(values) - 1, -1, -1):
        value = values[index]
        if math.isfinite(float(value)):
            running = value if not math.isfinite(float(running)) else max(value, running)
        result[index] = running
    return result


def atomic_write_text(path: Path, text: str) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(text, encoding="utf-8", newline="\n")
    os.replace(temporary, path)
    return path


def atomic_write_json(path: Path, payload: Mapping[str, object]) -> Path:
    return atomic_write_text(
        Path(path),
        json.dumps(payload, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=True) + "\n",
    )


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str] | None = None) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    chosen = list(fields or (list(rows[0].keys()) if rows else []))
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=chosen, extrasaction="ignore")
        if chosen:
            writer.writeheader()
            writer.writerows(rows)
    os.replace(temporary, path)
    return path


def read_csv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def parse_bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


__all__ = [
    "BASE_BETA_DEG",
    "BASE_ETA",
    "BASE_MU",
    "DELTA_TOL",
    "EPSILON_VALUES",
    "GridPoint",
    "K_MAX",
    "LOW_ANGLE_BETA_DEG",
    "LOW_ANGLE_ETA",
    "LOW_ANGLE_MU",
    "MAIN_OUTPUT_DIR",
    "MANIFEST_FIELDS",
    "REGRESSION_POINTS",
    "REQUESTED_ROOTS",
    "ROOT_GUARD_INDEX",
    "SMOKE_OUTPUT_DIR",
    "THIN_DIAGNOSTIC_LIMIT",
    "WORKFLOW_VERSION",
    "atomic_write_json",
    "atomic_write_text",
    "build_manifest",
    "cache_digest",
    "case_cache_identity",
    "geometry_identity",
    "group_counts",
    "parse_bool",
    "primary_settings",
    "read_csv",
    "select_points",
    "solver_configuration",
    "squared_frequency_delta",
    "strict_settings",
    "suffix_max",
    "true_safe_prefix",
    "write_csv",
]
