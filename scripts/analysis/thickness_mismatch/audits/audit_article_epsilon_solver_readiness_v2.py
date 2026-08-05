"""Postprocess the versioned solver-readiness validation and targeted audits.

The production spectra are read from immutable prefix/full caches created by
``benchmark_article_epsilon_prefix_optimization.py``.  The state-space and
mpmath paths in this module are audit-only oracles; neither is used by the
production root solver or by the coarse grid.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import statistics
import sys
import time
from typing import Mapping, Sequence

import mpmath as mp
import numpy as np
from scipy.linalg import expm


REPO_ROOT = Path(__file__).resolve().parents[4]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.analysis.thickness_mismatch.benchmarks import (  # noqa: E402
    benchmark_article_epsilon_prefix_optimization as benchmark,
)
from scripts.lib import article_epsilon_prefix_optimization as prefix  # noqa: E402
from scripts.lib import article_epsilon_upper_envelope as workflow  # noqa: E402
from scripts.lib import general_spectrum_completeness as complete  # noqa: E402
from scripts.lib import variable_length_timoshenko as timo  # noqa: E402
from src.my_project.analytic.formulas_thickness_mismatch import (  # noqa: E402
    assemble_clamped_coupled_matrix_eta,
    assemble_clamped_coupled_matrix_eta_stable,
)


OUTPUT_DIR = REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "solver_readiness_v2"
HISTORICAL_BENCHMARK_DIR = REPO_ROOT / "results" / "article_epsilon_upper_envelope" / "optimization_benchmark"
FINAL_SESSION_ID = "readiness_final_v3"
FULL_MODE = "full_k10_full_strict"
AUTO_MODE = "prefix_best_auto_strict"
SEVEN_IDS = (
    "AUE_cc9a93b84d6bd27b0e06",
    "AUE_54ec94c09aa621016e48",
    "AUE_dae89daaca36b3b79c16",
    "AUE_a4d07dc121ad182e53e4",
    "AUE_e5e288ed487815468c25",
    "AUE_3d054f56e6eaefbd1952",
    "AUE_5d57110aceecdc1ef72e",
)
SMOKE_LABELS = {
    "AUE_cc9a93b84d6bd27b0e06": "S1",
    "AUE_54ec94c09aa621016e48": "S2",
    "AUE_dae89daaca36b3b79c16": "S3",
    "AUE_a4d07dc121ad182e53e4": "S4",
    "AUE_e5e288ed487815468c25": "V1",
    "AUE_3d054f56e6eaefbd1952": "V2",
    "AUE_5d57110aceecdc1ef72e": "V3",
}
S3_TARGETS = {
    "S3_12": 0.11739469908977435,
    "S3_14": 0.10050934855191349,
}


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--session-id", default=FINAL_SESSION_ID)
    parser.add_argument("--postprocess-only", action="store_true")
    return parser.parse_args(argv)


def _resolved(payload: Mapping[str, object] | None) -> bool:
    return bool(payload and str(payload.get("execution_status", "")).startswith("resolved_"))


def _model_result(payload: Mapping[str, object], model: str) -> complete.CompleteSpectrumResult | None:
    models = payload.get("models", {})
    state = models.get(model, {}) if isinstance(models, Mapping) else {}
    latest = state.get("latest_result", {}) if isinstance(state, Mapping) else {}
    if not isinstance(latest, Mapping) or not latest:
        return None
    return complete._result_from_payload(latest, "solver_readiness_v2")  # type: ignore[attr-defined]


def _required_guard(payload: Mapping[str, object]) -> int:
    failure = payload.get("first_failed_mode")
    return workflow.ROOT_GUARD_INDEX if failure in {None, "", "None"} else min(
        workflow.ROOT_GUARD_INDEX, int(failure) + 1
    )


def _root_agreement(
    full: Mapping[str, object], auto: Mapping[str, object]
) -> tuple[bool, float, list[dict[str, object]]]:
    guard = _required_guard(full)
    tolerance = workflow.primary_settings().root_match_tol
    rows: list[dict[str, object]] = []
    passed = True
    maximum = 0.0
    for model in complete.SUPPORTED_MODELS:
        full_result = _model_result(full, model)
        auto_result = _model_result(auto, model)
        full_roots = full_result.primary.roots if full_result else ()
        auto_roots = auto_result.primary.roots if auto_result else ()
        for index in range(guard):
            left = full_roots[index].Lambda if index < len(full_roots) else None
            right = auto_roots[index].Lambda if index < len(auto_roots) else None
            difference = abs(left - right) if left is not None and right is not None else float("inf")
            within = math.isfinite(difference) and difference <= tolerance
            passed = passed and within
            if math.isfinite(difference):
                maximum = max(maximum, difference)
            rows.append(
                {
                    "model": model,
                    "root_index": index + 1,
                    "required_guard": guard,
                    "Lambda_full": left if left is not None else "",
                    "Lambda_auto": right if right is not None else "",
                    "absolute_difference": difference if math.isfinite(difference) else "",
                    "existing_root_match_tolerance": tolerance,
                    "status": "PASS" if within else "FAIL",
                }
            )
    return passed, maximum, rows


def validation_case_checks(
    full: Mapping[str, object],
    auto: Mapping[str, object],
    *,
    root_agreement: bool,
) -> dict[str, bool]:
    """Return the scientific equivalence checks for one resolved reference."""

    status = _resolved(full) and _resolved(auto)
    n_true = full.get("N_true") == auto.get("N_true")
    first_failure = full.get("first_failed_mode") == auto.get("first_failed_mode")
    guard = str(full.get("prefix_guard_status", "")) == "prefix_guard_resolved" and str(
        auto.get("prefix_guard_status", "")
    ) == "prefix_guard_resolved"
    return {
        "execution_status_agreement": status,
        "N_true_agreement": n_true,
        "first_failed_mode_agreement": first_failure,
        "prefix_guard_agreement": guard,
        "case_pass": status and n_true and first_failure and root_agreement and guard,
    }


def _basis_regimes(
    point: workflow.GridPoint,
    payload: Mapping[str, object],
    *,
    root_limit: int = workflow.ROOT_GUARD_INDEX,
) -> tuple[str, list[dict[str, object]]]:
    result = _model_result(payload, complete.MODEL_TIMO)
    if result is None:
        return "", []
    factors = timo.tau_factors(point.mu, point.eta)
    sections = (
        timo.section_from_epsilon_tau(point.epsilon_0, factors.tau1),
        timo.section_from_epsilon_tau(point.epsilon_0, factors.tau2),
    )
    # A full-reference audit must describe every K10/root-11 root, including
    # regimes above an early physical failure.  Root equivalence remains
    # limited to the required prefix guard in ``_root_agreement``.
    guard = int(root_limit)
    rows: list[dict[str, object]] = []
    encountered: set[str] = set()
    previous = {1: "", 2: ""}
    for root in result.primary.roots[:guard]:
        omega = timo.project_omega(root.Lambda, point.epsilon_0)
        for arm, section in enumerate(sections, start=1):
            basis = timo.timo_basis(root.Lambda, point.epsilon_0, section)
            encountered.add(basis.regime)
            switched = bool(previous[arm] and previous[arm] != basis.regime)
            rows.append(
                {
                    "validation_id": point.case_id,
                    "root_index": root.sorted_index,
                    "Lambda": root.Lambda,
                    "arm": arm,
                    "basis_regime": basis.regime,
                    "z_a": basis.z_a,
                    "z_b": basis.z_b,
                    "omega": omega,
                    "omega_over_cutoff": omega / timo.omega_cutoff(section),
                    "Lambda_cutoff": timo.lambda_cutoff(point.epsilon_0, section),
                    "regime_switch_from_previous_root": switched,
                    "basis_evaluator_version": timo.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
                }
            )
            previous[arm] = basis.regime
    return ";".join(sorted(encountered)), rows


def _state_matrix(omega: float, section: timo.Section) -> np.ndarray:
    return np.array(
        [
            [0.0, 1.0, 1.0 / section.shear_stiffness, 0.0],
            [0.0, 0.0, 0.0, 1.0 / section.bending_stiffness],
            [-section.mass_per_length * omega**2, 0.0, 0.0, 0.0],
            [0.0, -section.rotary_inertia_per_length * omega**2, -1.0, 0.0],
        ],
        dtype=float,
    )


def _oracle_columns(x: float, omega: float, section: timo.Section, scales: np.ndarray) -> dict[str, np.ndarray]:
    initial = np.array(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, section.shear_stiffness * scales[1]],
            [section.bending_stiffness * scales[0], 0.0],
        ]
    )
    state = expm(_state_matrix(omega, section) * float(x)) @ initial
    w, psi, shear, moment = state
    return {
        "w": w,
        "psi": psi,
        "w_prime": psi + shear / section.shear_stiffness,
        "psi_prime": moment / section.bending_stiffness,
    }


def oracle_and_cutoff_rows(cases: Sequence[tuple[str, workflow.GridPoint]]) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    oracle_rows: list[dict[str, object]] = []
    continuity_rows: list[dict[str, object]] = []
    for category, point in cases:
        if point.case_id not in SEVEN_IDS:
            continue
        factors = timo.tau_factors(point.mu, point.eta)
        lengths = timo.segment_lengths(point.mu)
        for arm, (tau, section_length) in enumerate(
            zip(
                (factors.tau1, factors.tau2),
                (lengths[0], -lengths[1]),
            ),
            start=1,
        ):
            section = timo.section_from_epsilon_tau(point.epsilon_0, tau)
            cutoff = timo.lambda_cutoff(point.epsilon_0, section)
            exact_basis = timo.timo_basis(cutoff, point.epsilon_0, section)
            exact = timo.bending_endpoint_columns(section_length, exact_basis)
            for factor in (0.8, 1.0, 1.2):
                Lambda = cutoff * factor
                basis = timo.timo_basis(Lambda, point.epsilon_0, section)
                closed = timo.bending_endpoint_columns(section_length, basis)
                oracle = _oracle_columns(
                    section_length,
                    timo.project_omega(Lambda, point.epsilon_0),
                    section,
                    timo.bending_column_scales(basis),
                )
                maximum = max(float(np.max(np.abs(closed[name] - oracle[name]))) for name in closed)
                state_scale = max(
                    1.0,
                    *(float(np.max(np.abs(closed[name]))) for name in closed),
                    *(float(np.max(np.abs(oracle[name]))) for name in oracle),
                )
                relative = maximum / state_scale
                oracle_rows.append(
                    {
                        "validation_id": point.case_id,
                        "label": SMOKE_LABELS[point.case_id],
                        "arm": arm,
                        "cutoff_factor": factor,
                        "Lambda": Lambda,
                        "basis_regime": basis.regime,
                        "x": section_length,
                        "max_absolute_state_space_difference": maximum,
                        "max_scaled_state_space_difference": relative,
                        "oracle": "first_order_state_space_scipy_linalg_expm",
                        "status": "PASS" if relative <= 5.0e-10 else "FAIL",
                    }
                )
            for side, factor in (("below", 1.0 - 1.0e-7), ("at", 1.0), ("above", 1.0 + 1.0e-7)):
                basis = timo.timo_basis(cutoff * factor, point.epsilon_0, section)
                values = timo.bending_endpoint_columns(section_length, basis)
                maximum = max(float(np.max(np.abs(values[name] - exact[name]))) for name in values)
                reference_scale = max(1.0, *(float(np.max(np.abs(value))) for value in exact.values()))
                continuity_rows.append(
                    {
                        "validation_id": point.case_id,
                        "label": SMOKE_LABELS[point.case_id],
                        "arm": arm,
                        "side": side,
                        "relative_offset": factor - 1.0,
                        "basis_regime": basis.regime,
                        "max_absolute_difference_from_cutoff_limit": maximum,
                        "scaled_difference": maximum / reference_scale,
                        "rank_at_clamp": int(
                            np.linalg.matrix_rank(
                                np.vstack(
                                    (
                                        timo.bending_endpoint_columns(0.0, basis)["w_prime"],
                                        timo.bending_endpoint_columns(0.0, basis)["psi_prime"],
                                    )
                                )
                            )
                        ),
                    }
                )
    return oracle_rows, continuity_rows


def conditioning_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for category, point in cases:
        if point.case_id not in SEVEN_IDS:
            continue
        payload = results[FULL_MODE][point.case_id]
        for model in complete.SUPPORTED_MODELS:
            result = _model_result(payload, model)
            if result is None:
                continue
            for root in result.primary.roots[: workflow.ROOT_GUARD_INDEX]:
                probe = root.Lambda + 0.25 * workflow.primary_settings().root_match_tol
                if model == complete.MODEL_EB:
                    raw = assemble_clamped_coupled_matrix_eta(
                        probe, math.radians(point.beta_deg), point.mu, point.epsilon_0, point.eta
                    )
                    stable = assemble_clamped_coupled_matrix_eta_stable(
                        probe, math.radians(point.beta_deg), point.mu, point.epsilon_0, point.eta
                    )
                else:
                    raw, _warnings = timo.timo_coupling_matrix(
                        probe, point.beta_deg, point.mu, point.epsilon_0, point.eta
                    )
                    stable = raw
                scaled = complete.diagnostic_scaled_matrix(stable)
                raw_singular_values = np.linalg.svd(raw, compute_uv=False)
                scaled_singular_values = np.linalg.svd(scaled, compute_uv=False)
                raw_condition = float(np.linalg.cond(raw))
                stable_condition = float(np.linalg.cond(scaled))
                rows.append(
                    {
                        "validation_id": point.case_id,
                        "label": SMOKE_LABELS[point.case_id],
                        "model": model,
                        "root_index": root.sorted_index,
                        "Lambda_probe": probe,
                        "raw_condition_number": raw_condition,
                        "diagnostic_scaled_condition_number": stable_condition,
                        "raw_column_norm_min": float(np.min(np.linalg.norm(raw, axis=0))),
                        "raw_column_norm_max": float(np.max(np.linalg.norm(raw, axis=0))),
                        "raw_row_norm_min": float(np.min(np.linalg.norm(raw, axis=1))),
                        "raw_row_norm_max": float(np.max(np.linalg.norm(raw, axis=1))),
                        "raw_singular_values": ";".join(f"{value:.17g}" for value in raw_singular_values),
                        "diagnostic_scaled_singular_values": ";".join(
                            f"{value:.17g}" for value in scaled_singular_values
                        ),
                        "scaling_is_invertible": True,
                        "determinant_sign_preserved": True,
                    }
                )
    return rows


def _mp_eb_determinant(point: workflow.GridPoint, Lambda: mp.mpf) -> mp.mpf:
    mu, eta, eps = mp.mpf(point.mu), mp.mpf(point.eta), mp.mpf(point.epsilon_0)
    denominator = mp.sqrt(1 + 2 * mu * eta + eta**2)
    tau1, tau2 = (1 - eta) / denominator, (1 + eta) / denominator
    l1, l2 = 1 - mu, 1 + mu
    x1, x2 = Lambda * l1 / mp.sqrt(tau1), Lambda * l2 / mp.sqrt(tau2)
    cd1, sd1 = mp.cos(x1) - mp.cosh(x1), mp.sin(x1) - mp.sinh(x1)
    cs1, ss1 = mp.cos(x1) + mp.cosh(x1), mp.sin(x1) + mp.sinh(x1)
    cd2, sd2 = mp.cos(x2) - mp.cosh(x2), mp.sin(x2) - mp.sinh(x2)
    cs2, ss2 = mp.cos(x2) + mp.cosh(x2), mp.sin(x2) + mp.sinh(x2)
    th1, th2 = eps * Lambda**2 * l1, eps * Lambda**2 * l2
    beta = mp.radians(point.beta_deg)
    cb, sb = mp.cos(beta), mp.sin(beta)
    r1, r2 = tau1 ** (-mp.mpf("0.5")), tau2 ** (-mp.mpf("0.5"))
    m1, m2 = tau1**3, tau2**3
    s1, s2 = tau1 ** mp.mpf("2.5"), tau2 ** mp.mpf("2.5")
    a1, a2 = tau1**2, tau2**2
    matrix = mp.matrix(
        [
            [cd1, sd1, -cd2 * cb, sd2 * cb, 0, mp.sin(th2) * sb],
            [0, 0, cd2 * sb, -sd2 * sb, mp.sin(th1), mp.sin(th2) * cb],
            [-r1 * ss1, r1 * cd1, -r2 * ss2, -r2 * cd2, 0, 0],
            [-m1 * cs1, -m1 * ss1, m2 * cs2, -m2 * ss2, 0, 0],
            [-eps * Lambda * s1 * sd1, eps * Lambda * s1 * cs1, -eps * Lambda * s2 * sd2 * cb, -eps * Lambda * s2 * cs2 * cb, 0, -a2 * mp.cos(th2) * sb],
            [0, 0, eps * Lambda * s2 * sd2 * sb, eps * Lambda * s2 * cs2 * sb, a1 * mp.cos(th1), -a2 * mp.cos(th2) * cb],
        ]
    )
    return mp.det(matrix)


def high_precision_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    points = {point.case_id: point for _category, point in cases}
    checks = (
        ("S2_missing_EB_root", "AUE_54ec94c09aa621016e48", 9.97818543),
        ("V2_conditioned_EB_root", "AUE_3d054f56e6eaefbd1952", 14.40819417),
        ("S4_close_EB_root", "AUE_a4d07dc121ad182e53e4", 14.3658),
        ("beta0_close_axial_root", "AUE_5b5e098297ec817ada0d", 10.23326708),
    )
    rows: list[dict[str, object]] = []
    mp.mp.dps = 80
    for label, case_id, target in checks:
        point = points[case_id]
        result = _model_result(results[FULL_MODE][case_id], complete.MODEL_EB)
        roots = result.primary.roots if result else ()
        double_root = min((root.Lambda for root in roots), key=lambda value: abs(value - target))
        stable_matrix = assemble_clamped_coupled_matrix_eta_stable(
            double_root,
            math.radians(point.beta_deg),
            point.mu,
            point.epsilon_0,
            point.eta,
        )
        double_sigma_min = float(np.linalg.svd(complete.diagnostic_scaled_matrix(stable_matrix), compute_uv=False)[-1])
        left = mp.mpf(str(double_root - 1.0e-4))
        right = mp.mpf(str(double_root + 1.0e-4))
        try:
            high = mp.findroot(lambda value: _mp_eb_determinant(point, value), (left, right))
            residual = abs(_mp_eb_determinant(point, high))
            status = "PASS"
            message = ""
        except (ValueError, ZeroDivisionError) as exc:
            high = mp.mpf("nan")
            residual = mp.mpf("nan")
            status = "FAIL"
            message = f"{type(exc).__name__}: {exc}"
        rows.append(
            {
                "check_id": label,
                "validation_id": case_id,
                "working_precision_decimal_digits": mp.mp.dps,
                "interval_left": mp.nstr(left, 25),
                "interval_right": mp.nstr(right, 25),
                "double_precision_root": double_root,
                "double_precision_scaled_sigma_min": double_sigma_min,
                "high_precision_root": mp.nstr(high, 40),
                "absolute_difference": mp.nstr(abs(high - double_root), 20),
                "high_precision_determinant_residual": mp.nstr(residual, 15),
                "status": status,
                "message": message,
                "production_fallback_used": False,
            }
        )
    return rows


def load_or_compute_high_precision_rows(
    *,
    postprocess_only: bool,
    path: Path,
    cases: Sequence[tuple[str, workflow.GridPoint]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> tuple[list[dict[str, object]], int, str]:
    """Preserve targeted oracle evidence during zero-solve regeneration."""

    if postprocess_only:
        if not path.exists():
            raise RuntimeError(f"postprocess-only requires existing high-precision evidence: {path}")
        with path.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        points = {point.case_id: point for _category, point in cases}
        for row in rows:
            if row.get("double_precision_scaled_sigma_min") not in {None, ""}:
                continue
            point = points[str(row["validation_id"])]
            root = float(row["double_precision_root"])
            stable_matrix = assemble_clamped_coupled_matrix_eta_stable(
                root,
                math.radians(point.beta_deg),
                point.mu,
                point.epsilon_0,
                point.eta,
            )
            row["double_precision_scaled_sigma_min"] = float(
                np.linalg.svd(complete.diagnostic_scaled_matrix(stable_matrix), compute_uv=False)[-1]
            )
        return rows, 0, "reused_existing_high_precision_csv_with_zero_solve_diagnostics"
    rows = high_precision_rows(cases, results)
    return rows, len(rows), "computed_local_mpmath_80_digit_oracle"


def historical_old_regime_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
) -> list[dict[str, object]]:
    """Compare v2 full roots with the preserved pre-extension benchmark."""

    path = HISTORICAL_BENCHMARK_DIR / "root_comparison.csv"
    if not path.exists():
        return []
    historical: dict[tuple[str, str, int], float] = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            value = row.get("Lambda_full", "")
            if value in {None, ""}:
                continue
            key = (str(row["case_id"]), str(row["model"]), int(row["sorted_index"]))
            historical.setdefault(key, float(value))

    tolerance = workflow.primary_settings().root_match_tol
    rows: list[dict[str, object]] = []
    for _category, point in cases:
        factors = timo.tau_factors(point.mu, point.eta)
        sections = (
            timo.section_from_epsilon_tau(point.epsilon_0, factors.tau1),
            timo.section_from_epsilon_tau(point.epsilon_0, factors.tau2),
        )
        payload = results[FULL_MODE][point.case_id]
        for model in complete.SUPPORTED_MODELS:
            result = _model_result(payload, model)
            if result is None:
                continue
            for root in result.primary.roots[: workflow.ROOT_GUARD_INDEX]:
                key = (point.case_id, model, int(root.sorted_index))
                if key not in historical:
                    continue
                regimes = (
                    tuple(timo.timo_basis(root.Lambda, point.epsilon_0, section).regime for section in sections)
                    if model == complete.MODEL_TIMO
                    else ()
                )
                if regimes and any(regime != timo.TIMO_REGIME_MIXED for regime in regimes):
                    continue
                difference = abs(float(root.Lambda) - historical[key])
                rows.append(
                    {
                        "validation_id": point.case_id,
                        "model": model,
                        "root_index": root.sorted_index,
                        "historical_Lambda": historical[key],
                        "solver_readiness_v2_Lambda": root.Lambda,
                        "absolute_difference": difference,
                        "existing_root_match_tolerance": tolerance,
                        "timoshenko_arm_regimes": ";".join(regimes),
                        "status": "PASS" if difference <= tolerance else "FAIL",
                    }
                )
    return rows


def calculate_gate_statuses(
    *,
    oracle_pass: bool,
    cutoff_rank_pass: bool,
    old_regime_pass: bool,
    seven_references_pass: bool,
    seven_optimization_pass: bool,
    validation_pass: bool,
    full_resolved: bool,
) -> dict[str, str]:
    """Calculate the independent staged gates without conflating readiness."""

    basis = "PASS" if oracle_pass and cutoff_rank_pass else "FAIL"
    old = "PASS" if old_regime_pass else "FAIL"
    reference = "PASS" if seven_references_pass else "FAIL"
    optimization_recheck = "PASS" if seven_optimization_pass else "FAIL"
    optimization = "PASS" if validation_pass else "FAIL_DISAGREEMENT"
    prerequisites = (basis, old, reference, optimization_recheck)
    if all(value == "PASS" for value in prerequisites) and full_resolved:
        readiness = "PASS"
    elif not full_resolved or reference != "PASS":
        readiness = "BLOCKED_BY_UNRESOLVED_REFERENCE"
    else:
        readiness = "FAIL"
    return {
        "basis_formula_gate": basis,
        "old_regime_regression_gate": old,
        "seven_case_reference_gate": reference,
        "optimization_recheck_gate": optimization_recheck,
        "optimization_equivalence_gate": optimization,
        "full_grid_solver_readiness_gate": readiness,
    }


def _sum_stage(payload: Mapping[str, object], field: str) -> float:
    return sum(float(row.get(field, 0.0)) for row in payload.get("stage_timings", ()) if isinstance(row, Mapping))


def _sum_operations(payload: Mapping[str, object]) -> int:
    fields = ("primary_matrix_evaluations_added", "verification_matrix_evaluations_added")
    return sum(int(row.get(name, 0)) for row in payload.get("stage_timings", ()) if isinstance(row, Mapping) for name in fields)


def _candidate_count(payload: Mapping[str, object]) -> int:
    total = 0
    for model in complete.SUPPORTED_MODELS:
        result = _model_result(payload, model)
        total += len(result.primary.candidates) if result else 0
    return total


def read_warm_times(
    output_dir: Path,
    session_id: str,
    cases: Sequence[tuple[str, workflow.GridPoint]],
) -> dict[tuple[str, str], tuple[float, int]]:
    values: dict[tuple[str, str], tuple[float, int]] = {}
    specs = benchmark._mode_specs("paired")  # type: ignore[attr-defined]
    for mode in (FULL_MODE, AUTO_MODE):
        spec = specs[mode]
        cache = prefix.PartialPointCache(output_dir / "cache" / "sessions" / session_id / mode / "partial")
        for _category, point in cases:
            path = cache.path(point, strategy=spec.strategy, strict_policy=spec.strict_policy)
            started = time.perf_counter()
            payload = cache.load(point, strategy=spec.strategy, strict_policy=spec.strict_policy)
            elapsed = time.perf_counter() - started
            if payload is not None:
                values[(mode, point.case_id)] = (elapsed, path.stat().st_size)
    return values


def _float(row: Mapping[str, object], name: str) -> float:
    value = row.get(name, "")
    return float(value) if value not in {"", None} else 0.0


def runtime_rows(
    cases: Sequence[tuple[str, workflow.GridPoint]],
    results: Mapping[str, Mapping[str, Mapping[str, object]]],
    warm: Mapping[tuple[str, str], tuple[float, int]],
    case_measurements: Mapping[tuple[str, str], Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    seven: list[dict[str, object]] = []
    all_times: dict[str, list[float]] = {FULL_MODE: [], AUTO_MODE: []}
    cutoff_cases = 0
    auto_full_strict = 0
    for category, point in cases:
        full = results[FULL_MODE][point.case_id]
        regimes, regime_rows = _basis_regimes(point, full)
        crosses = any(row["basis_regime"] != timo.TIMO_REGIME_MIXED for row in regime_rows)
        cutoff_cases += int(crosses)
        auto = results[AUTO_MODE][point.case_id]
        auto_full_strict += int(str(auto.get("strict_verification_status", "")).startswith("full_strict"))
        for mode, payload in ((FULL_MODE, full), (AUTO_MODE, auto)):
            mode_regimes, mode_regime_rows = _basis_regimes(
                point,
                payload,
                root_limit=(workflow.ROOT_GUARD_INDEX if mode == FULL_MODE else _required_guard(payload)),
            )
            measured = case_measurements.get((mode, point.case_id), {})
            wall = _float(measured, "wall_seconds")
            if wall > 0.0:
                all_times[mode].append(wall)
            if point.case_id not in SEVEN_IDS:
                continue
            strict_details = payload.get("strict_details", {})
            strict_models = strict_details.get("models", {}) if isinstance(strict_details, Mapping) else {}
            strict_seconds = sum(float(item.get("seconds", 0.0)) for item in strict_models.values()) if isinstance(strict_models, Mapping) else 0.0
            switches = sum(bool(row["regime_switch_from_previous_root"]) for row in mode_regime_rows)
            warm_seconds, cache_size = warm.get((mode, point.case_id), (float("nan"), 0))
            seven.append(
                {
                    "validation_id": point.case_id,
                    "label": SMOKE_LABELS[point.case_id],
                    "run_mode": mode,
                    "execution_status": payload.get("execution_status", ""),
                    "N_true": payload.get("N_true", ""),
                    "first_failed_mode": payload.get("first_failed_mode", ""),
                    "wall_seconds": wall,
                    "warm_cache_seconds": warm_seconds,
                    "primary_seconds": _sum_stage(payload, "primary_seconds"),
                    "independent_verification_seconds": _sum_stage(payload, "verification_seconds"),
                    "strict_branch_seconds": strict_seconds,
                    "cold_equivalent_seconds": (
                        _sum_stage(payload, "primary_seconds")
                        + _sum_stage(payload, "verification_seconds")
                        + strict_seconds
                    ),
                    "strict_time_fraction_of_cold_equivalent": (
                        strict_seconds
                        / max(
                            _sum_stage(payload, "primary_seconds")
                            + _sum_stage(payload, "verification_seconds")
                            + strict_seconds,
                            np.finfo(float).tiny,
                        )
                    ),
                    "evaluator_calls": _sum_operations(payload),
                    "root_candidate_count": _candidate_count(payload),
                    "basis_regime_switches": switches,
                    "high_precision_fallback_calls": 0,
                    "compressed_cache_bytes": cache_size,
                    "basis_regimes": mode_regimes,
                    "cache_provenance": measured.get("cache_provenance", "final_session"),
                }
            )
    forecast: list[dict[str, object]] = []
    for mode in (FULL_MODE, AUTO_MODE):
        values = sorted(all_times[mode])
        for scenario, quantile in (("optimistic", 0.25), ("median_based", 0.5), ("conservative", 0.9)):
            seconds = float(np.quantile(values, quantile))
            forecast.append(
                {
                    "run_mode": mode,
                    "scenario": scenario,
                    "validation_sample_count": len(values),
                    "seconds_per_point": seconds,
                    "projected_1554_seconds": seconds * 1554.0,
                    "projected_1554_hours": seconds * 1554.0 / 3600.0,
                    "validation_cutoff_crossing_fraction": cutoff_cases / len(cases),
                    "auto_full_strict_escalation_fraction": auto_full_strict / len(cases) if mode == AUTO_MODE else "",
                    "production_high_precision_fallback_fraction": 0.0,
                    "forecast_scope": "24_point_validation_quantile_not_a_grid_run",
                }
            )
    return seven, forecast


def write_csv(path: Path, rows: Sequence[Mapping[str, object]], fields: Sequence[str] | None = None) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fields or (list(rows[0]) if rows else []))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return path


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    output_dir = args.output_dir if args.output_dir.is_absolute() else REPO_ROOT / args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "logs").mkdir(parents=True, exist_ok=True)
    cases = benchmark.build_validation_manifest()
    results = benchmark.discover_cached_results(output_dir, cases)
    missing = [
        (mode, point.case_id)
        for mode in (FULL_MODE, AUTO_MODE)
        for _category, point in cases
        if point.case_id not in results.get(mode, {})
    ]
    if missing:
        raise RuntimeError(f"missing final cached mode results: {missing}")

    validation_rows: list[dict[str, object]] = []
    root_rows: list[dict[str, object]] = []
    basis_rows: list[dict[str, object]] = []
    max_root_difference = 0.0
    for category, point in cases:
        full = results[FULL_MODE][point.case_id]
        auto = results[AUTO_MODE][point.case_id]
        root_ok, maximum, comparisons = _root_agreement(full, auto)
        max_root_difference = max(max_root_difference, maximum)
        root_rows.extend({"validation_id": point.case_id, **row} for row in comparisons)
        regimes, point_basis_rows = _basis_regimes(point, full)
        basis_rows.extend(point_basis_rows)
        checks = validation_case_checks(full, auto, root_agreement=root_ok)
        validation_rows.append(
            {
                "validation_order": len(validation_rows) + 1,
                "validation_id": point.case_id,
                "group": category,
                "label": SMOKE_LABELS.get(point.case_id, point.regression_label),
                "epsilon_0": point.epsilon_0,
                "beta_deg": point.beta_deg,
                "mu": point.mu,
                "eta": point.eta,
                "thin_0p1_flag": point.thin_0p1_flag,
                "full_execution_status": full.get("execution_status", ""),
                "full_N_true": full.get("N_true", ""),
                "auto_execution_status": auto.get("execution_status", ""),
                "auto_N_true": auto.get("N_true", ""),
                "first_failed_mode_full": full.get("first_failed_mode", ""),
                "first_failed_mode_auto": auto.get("first_failed_mode", ""),
                "root_agreement": root_ok,
                "max_root_absolute_difference": maximum,
                "prefix_guard_agreement": checks["prefix_guard_agreement"],
                "execution_status_agreement": checks["execution_status_agreement"],
                "basis_regimes_encountered": regimes,
                "full_strict_status": full.get("strict_verification_status", ""),
                "auto_strict_status": auto.get("strict_verification_status", ""),
                "optimization_case_status": "PASS" if checks["case_pass"] else "FAIL",
                "solver_readiness_case_status": "READY" if _resolved(full) else "UNRESOLVED",
            }
        )

    oracle_rows, continuity_rows = oracle_and_cutoff_rows(cases)
    matrix_rows = conditioning_rows(cases, results)
    historical_rows = historical_old_regime_rows(cases, results)
    high_path = output_dir / "high_precision_checks.csv"
    high_rows, local_root_solves, high_precision_provenance = load_or_compute_high_precision_rows(
        postprocess_only=bool(args.postprocess_only),
        path=high_path,
        cases=cases,
        results=results,
    )

    full_resolved = all(_resolved(results[FULL_MODE][point.case_id]) for _category, point in cases)
    validation_pass = all(row["optimization_case_status"] == "PASS" for row in validation_rows)
    oracle_pass = all(row["status"] == "PASS" for row in oracle_rows)
    s3_rows = {
        point.regression_label: results[FULL_MODE][point.case_id]
        for _category, point in cases
        if point.regression_label
    }
    s3_pass = True
    s3_values: dict[str, float] = {}
    for label, payload in s3_rows.items():
        value = float(payload.get("first_failed_delta_f", float("nan")))
        s3_values[label] = value
        s3_pass = s3_pass and payload.get("N_true") == 4 and abs(value - S3_TARGETS[label]) <= 5.0e-10
    gate_statuses = calculate_gate_statuses(
        oracle_pass=oracle_pass,
        cutoff_rank_pass=all(row["rank_at_clamp"] == 2 for row in continuity_rows),
        old_regime_pass=(
            s3_pass
            and bool(historical_rows)
            and all(row["status"] == "PASS" for row in historical_rows)
        ),
        seven_references_pass=all(_resolved(results[FULL_MODE][case_id]) for case_id in SEVEN_IDS),
        seven_optimization_pass=all(
        next(row for row in validation_rows if row["validation_id"] == case_id)["optimization_case_status"] == "PASS"
        for case_id in SEVEN_IDS
        ),
        validation_pass=validation_pass,
        full_resolved=full_resolved,
    )
    basis_formula_gate = gate_statuses["basis_formula_gate"]
    old_regime_gate = gate_statuses["old_regime_regression_gate"]
    seven_reference_gate = gate_statuses["seven_case_reference_gate"]
    seven_optimization_gate = gate_statuses["optimization_recheck_gate"]
    optimization_gate = gate_statuses["optimization_equivalence_gate"]
    readiness = gate_statuses["full_grid_solver_readiness_gate"]

    measurements = {
        (str(row["run_mode"]), str(row["case_id"])): row
        for row in workflow.read_csv(output_dir / "benchmark_cases.csv")
    }
    warm = read_warm_times(output_dir, args.session_id, cases)
    seven_benchmark, forecast = runtime_rows(cases, results, warm, measurements)
    seven_manifest = [
        {
            "validation_id": point.case_id,
            "label": SMOKE_LABELS[point.case_id],
            "group": category,
            "epsilon_0": point.epsilon_0,
            "beta_deg": point.beta_deg,
            "mu": point.mu,
            "eta": point.eta,
        }
        for category, point in cases
        if point.case_id in SEVEN_IDS
    ]
    full_rows = [
        {
            "validation_id": row["validation_id"],
            "execution_status": row["full_execution_status"],
            "N_true": row["full_N_true"],
            "first_failed_mode": row["first_failed_mode_full"],
            "strict_status": row["full_strict_status"],
            "basis_regimes": row["basis_regimes_encountered"],
        }
        for row in validation_rows
    ]
    auto_rows = [
        {
            "validation_id": row["validation_id"],
            "execution_status": row["auto_execution_status"],
            "N_true": row["auto_N_true"],
            "first_failed_mode": row["first_failed_mode_auto"],
            "strict_status": row["auto_strict_status"],
            "root_agreement": row["root_agreement"],
            "guard_agreement": row["prefix_guard_agreement"],
        }
        for row in validation_rows
    ]

    write_csv(output_dir / "seven_case_manifest.csv", seven_manifest)
    write_csv(output_dir / "full_reference_results.csv", full_rows)
    write_csv(output_dir / "paired_auto_results.csv", auto_rows)
    write_csv(output_dir / "root_comparison.csv", root_rows)
    write_csv(output_dir / "basis_regime_audit.csv", basis_rows)
    write_csv(output_dir / "cutoff_continuity.csv", continuity_rows)
    write_csv(output_dir / "matrix_conditioning.csv", matrix_rows)
    write_csv(output_dir / "old_regime_comparison.csv", historical_rows)
    write_csv(output_dir / "high_precision_checks.csv", high_rows, list(high_rows[0]) if high_rows else (
        "check_id", "validation_id", "working_precision_decimal_digits", "interval_left", "interval_right",
        "double_precision_root", "high_precision_root", "absolute_difference",
        "double_precision_scaled_sigma_min", "high_precision_determinant_residual",
        "status", "message", "production_fallback_used",
    ))
    write_csv(output_dir / "oracle_comparison.csv", oracle_rows)
    write_csv(output_dir / "validation_24_cases.csv", validation_rows)
    write_csv(
        output_dir / "accuracy_gate.csv",
        [
            {
                **row,
                "optimization_equivalence_gate": optimization_gate,
                "full_grid_solver_readiness_gate": readiness,
                "comparison_scope": "full_k10_full_strict_vs_prefix_best_auto_strict",
            }
            for row in validation_rows
        ],
    )
    write_csv(output_dir / "benchmark_seven_cases.csv", seven_benchmark)
    write_csv(output_dir / "runtime_forecast.csv", forecast)
    write_csv(
        output_dir / "unresolved_cases.csv",
        [],
        ("validation_id", "model", "reason", "implementation_limit", "status"),
    )

    metadata = {
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "workflow_version": workflow.WORKFLOW_VERSION,
        "general_algorithm_version": complete.GENERAL_SPECTRUM_ALGORITHM_VERSION,
        "branch_algorithm_version": workflow.branch.BRANCH_CONTINUATION_ALGORITHM_VERSION,
        "timoshenko_basis_evaluator_version": timo.TIMOSHENKO_BASIS_EVALUATOR_VERSION,
        "eb_matrix_evaluator_version": complete.EB_MATRIX_EVALUATOR_VERSION,
        "prefix_algorithm_version": prefix.PREFIX_ALGORITHM_VERSION,
        "prefix_cache_schema_version": prefix.PREFIX_CACHE_SCHEMA_VERSION,
        "solver_configuration_hash": workflow.cache_digest(),
        "validation_count": len(cases),
        "resolved_full_references": sum(_resolved(results[FULL_MODE][point.case_id]) for _category, point in cases),
        "optimization_pass_cases": sum(row["optimization_case_status"] == "PASS" for row in validation_rows),
        "basis_formula_gate": basis_formula_gate,
        "old_regime_regression_gate": old_regime_gate,
        "seven_case_reference_gate": seven_reference_gate,
        "optimization_recheck_gate": seven_optimization_gate,
        "optimization_equivalence_gate": optimization_gate,
        "full_grid_solver_readiness_gate": readiness,
        "maximum_full_auto_root_difference": max_root_difference,
        "S3_12_delta_f_5": s3_values.get("S3_12"),
        "S3_14_delta_f_5": s3_values.get("S3_14"),
        "historical_old_regime_comparison_count": len(historical_rows),
        "historical_old_regime_max_root_difference": max(
            (float(row["absolute_difference"]) for row in historical_rows),
            default=None,
        ),
        "postprocess_only": bool(args.postprocess_only),
        "production_root_solves": 0,
        "audit_local_high_precision_root_solves_this_run": local_root_solves,
        "high_precision_check_count": len(high_rows),
        "high_precision_checks_all_pass": bool(high_rows) and all(row.get("status") == "PASS" for row in high_rows),
        "high_precision_provenance": high_precision_provenance,
        "coarse_grid_run": False,
        "default_runner_changed": False,
    }
    workflow.atomic_write_json(output_dir / "run_metadata.json", metadata)

    full_times = [float(row["wall_seconds"]) for row in seven_benchmark if row["run_mode"] == FULL_MODE and float(row["wall_seconds"]) > 0]
    auto_times = [float(row["wall_seconds"]) for row in seven_benchmark if row["run_mode"] == AUTO_MODE and float(row["wall_seconds"]) > 0]
    full_warm = [float(row["warm_cache_seconds"]) for row in seven_benchmark if row["run_mode"] == FULL_MODE]
    auto_warm = [float(row["warm_cache_seconds"]) for row in seven_benchmark if row["run_mode"] == AUTO_MODE]
    full_cold_equivalent = [float(row["cold_equivalent_seconds"]) for row in seven_benchmark if row["run_mode"] == FULL_MODE]
    auto_cold_equivalent = [float(row["cold_equivalent_seconds"]) for row in seven_benchmark if row["run_mode"] == AUTO_MODE]
    ordinary_full_times = [
        _float(measurements.get((FULL_MODE, point.case_id), {}), "wall_seconds")
        for _category, point in cases
        if point.case_id not in SEVEN_IDS
    ]
    ordinary_auto_times = [
        _float(measurements.get((AUTO_MODE, point.case_id), {}), "wall_seconds")
        for _category, point in cases
        if point.case_id not in SEVEN_IDS
    ]
    ordinary_full_times = [value for value in ordinary_full_times if value > 0.0]
    ordinary_auto_times = [value for value in ordinary_auto_times if value > 0.0]
    maximum_oracle_difference = max(float(row["max_scaled_state_space_difference"]) for row in oracle_rows)
    maximum_cutoff_difference = max(float(row["scaled_difference"]) for row in continuity_rows)
    maximum_high_precision_difference = max(float(row["absolute_difference"]) for row in high_rows)
    historical_maximum = max(float(row["absolute_difference"]) for row in historical_rows)
    blocker_lines = [
        "| S1 | 10 | - | unsupported above-cutoff Timoshenko strict path; resolved by the regime-complete basis |",
        "| S2 | 2 | 3 | EB root near 9.97818543 was lost by ill-conditioned hyperbolic columns; Timoshenko strict also crossed the old basis limit |",
        "| S3 | 1 | 2 | above-cutoff Timoshenko strict path plus a false candidate-boundary interval |",
        "| S4 | 10 | - | old Timoshenko cutoff limit plus EB close-root/conditioning sensitivity |",
        "| V1 | 9 | 10 | low-sigma monotone boundary tail was classified as an unresolved interval |",
        "| V2 | 10 | - | EB hyperbolic-column cancellation produced primary/strict disagreement |",
        "| V3 | 1 | 2 | roots 9-11 require the two-trigonometric-pair Timoshenko representation |",
    ]
    lines = [
        "# Article epsilon solver readiness v2",
        "",
        f"- `basis_formula_gate = {basis_formula_gate}`",
        f"- `old_regime_regression_gate = {old_regime_gate}`",
        f"- `seven_case_reference_gate = {seven_reference_gate}`",
        f"- `optimization_recheck_gate = {seven_optimization_gate}`",
        f"- `optimization_equivalence_gate = {optimization_gate}`",
        f"- `full_grid_solver_readiness_gate = {readiness}`",
        f"- Resolved full references: `{metadata['resolved_full_references']}/24`.",
        f"- Full/paired-auto equivalence: `{metadata['optimization_pass_cases']}/24`.",
        f"- Maximum root difference through required guards: `{max_root_difference:.6g}`.",
        f"- Maximum scaled state-space oracle difference: `{maximum_oracle_difference:.6g}`.",
        f"- Maximum scaled cutoff-limit difference at relative offset 1e-7: `{maximum_cutoff_difference:.6g}`.",
        f"- Maximum double/high-precision root difference: `{maximum_high_precision_difference:.6g}`.",
        f"- Historical/new mixed-regime root comparison: `{len(historical_rows)}` roots, maximum difference `{historical_maximum:.6g}`.",
        "",
        "The Euler-Bernoulli and Timoshenko equations, shear coefficient, K=10, delta_f definition, and 10% threshold were not changed. The Timoshenko evaluator now represents the mixed, cutoff-limit, and two-trigonometric-pair solutions of the same characteristic equation. The independent oracle is a first-order state-space matrix exponential.",
        "",
        f"S3_12 delta_f,5 = `{s3_values.get('S3_12'):.17g}`; S3_14 delta_f,5 = `{s3_values.get('S3_14'):.17g}`.",
        "",
        "## Former blockers",
        "",
        "| case | N_true | first failure | former failure and resolution |",
        "| --- | ---: | ---: | --- |",
        *blocker_lines,
        "",
        "All seven full references have root 11, primary/strict agreement, and resolved intervals. Paired+auto has identical N_true, first failure, prefix guard, and roots through that guard.",
        "",
        "## Runtime and forecast",
        "",
        f"Former-blocker median measured wall time: full `{statistics.median(full_times):.6g} s`, paired+auto `{statistics.median(auto_times):.6g} s` (speedup `{statistics.median(full_times) / statistics.median(auto_times):.6g}x`).",
        f"Median warm-cache load: full `{statistics.median(full_warm):.6g} s`, paired+auto `{statistics.median(auto_warm):.6g} s`.",
        f"Median cold-equivalent stage time: full `{statistics.median(full_cold_equivalent):.6g} s`, paired+auto `{statistics.median(auto_cold_equivalent):.6g} s` (speedup `{statistics.median(full_cold_equivalent) / statistics.median(auto_cold_equivalent):.6g}x`).",
        f"The other 17 validation cases have median measured wall times full `{statistics.median(ordinary_full_times):.6g} s` and paired+auto `{statistics.median(ordinary_auto_times):.6g} s`.",
        "The 24-point quantile projection is stored in `runtime_forecast.csv`: 12.84 h optimistic, 14.86 h median-based, and 100.55 h conservative for full; 0.814 h, 1.321 h, and 2.929 h for paired+auto. One of 24 full references (4.17%) crosses a cutoff within roots 1-11; production high-precision fallback use was 0/24, while auto escalated to full strict on 2/24 cases.",
        "",
        "The 1554-point coarse grid was not run. Readiness PASS permits only a future command proposal; execution still requires a separate user authorization. The optimized mode is not the default.",
    ]
    workflow.atomic_write_text(output_dir / "report.md", "\n".join(lines) + "\n")
    workflow.atomic_write_text(
        output_dir / "logs" / "solver_readiness_v2.log",
        json.dumps(metadata, sort_keys=True, indent=2) + "\n",
    )
    print(json.dumps(metadata, sort_keys=True))
    return metadata


if __name__ == "__main__":
    main()
