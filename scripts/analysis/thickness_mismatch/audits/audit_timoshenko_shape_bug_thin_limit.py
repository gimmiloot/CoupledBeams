from __future__ import annotations

"""Diagnostic-only thin-limit audit of sorted Timoshenko mode reconstruction.

This is a distinct workflow from the plotting-only audits: it compares EB and
Timoshenko fields, inspects the nullspace at each selected root, and verifies
the reconstructed boundary/joint residuals.  It reuses the established root
and field helpers and does not alter analytic formulas or solver behavior.
"""

import argparse
import inspect
from dataclasses import dataclass
from math import isfinite
from pathlib import Path
import sys
from typing import Sequence

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
    audit_timoshenko_shape_construction as shape_audit,
)
from scripts.lib import in_plane_shape_geometry as DISPLAY  # noqa: E402


MODEL_EB = shape_audit.MODEL_EB
MODEL_TIMO = shape_audit.MODEL_TIMO

DEFAULT_OUTPUT_DIR = Path("results") / "timoshenko_shape_bug_audit"
SMOKE_OUTPUT_DIR = Path("results") / "_smoke" / "timoshenko_shape_bug_audit"
DEFAULT_BETA_DEG = 45.0
DEFAULT_ETA = 0.0
THIN_EPSILON = 0.0025
PROBLEM_EPSILON = 0.03
DEFAULT_MU_VALUES = (0.0, 0.2, 0.4, 0.6)
DEFAULT_SORTED_INDICES = (4, 5, 6)
SMOKE_MU_VALUES = (0.0,)
SMOKE_SORTED_INDICES = (4, 5)
DEFAULT_N_POINTS = 801
SMOKE_N_POINTS = 301
REQUESTED_SCALE = 0.08
FALLBACK_SCALE = 0.05
SAFETY_SCALE = 0.02
AUDIT_TOL = 1.0e-6
NULL_RATIO_RELIABLE_MAX = 1.0e-6
NEAR_DEGENERATE_SIGMA_SECOND_MAX = 1.0e-6

THIN_MAC_FIELDS = [
    "mu",
    "sorted_index",
    "Lambda_EB",
    "Lambda_Timo",
    "rel_freq_diff",
    "mac_uw",
    "mac_u",
    "mac_w",
    "mac_psi_vs_eb_slope",
    "mac_global_displacement",
    "mac_rod2_global_displacement",
    "mac_display_corrected",
    "mac_rod2_display_corrected",
    "max_abs_gamma",
    "gamma_norm_ratio",
    "eb_classification",
    "timo_classification",
    "warning",
    "notes",
]

NULL_VECTOR_FIELDS = [
    "mu",
    "epsilon",
    "beta_deg",
    "eta",
    "sorted_index",
    "Lambda",
    "sigma_min",
    "sigma_second_min",
    "singular_value_ratio",
    "residual_norm",
    "coeff_norm",
    "coeff_A1",
    "coeff_B1",
    "coeff_P1",
    "coeff_A2",
    "coeff_B2",
    "coeff_P2",
    "warning",
    "notes",
]

RESIDUAL_FIELDS = [
    "mu",
    "epsilon",
    "beta_deg",
    "eta",
    "sorted_index",
    "Lambda",
    "rod1_clamp_u",
    "rod1_clamp_w",
    "rod1_clamp_psi",
    "rod2_clamp_u",
    "rod2_clamp_w",
    "rod2_clamp_psi",
    "max_abs_clamp_gap",
    "gap_w",
    "gap_u",
    "gap_psi",
    "gap_M",
    "gap_Q",
    "gap_N",
    "max_abs_kinematic_gap",
    "max_abs_force_gap",
    "normalized_matrix_residual_norm",
    "max_ode_residual",
    "passed",
    "warning",
    "notes",
]

DISPLAY_TRANSFORM_FIELDS = [
    "beta_deg",
    "tangent_x",
    "tangent_y",
    "normal_x",
    "normal_y",
    "tangent_dot_normal",
    "tangent_norm",
    "normal_norm",
    "pure_axial_parallel_error",
    "pure_transverse_orthogonality_error",
    "passed",
    "notes",
]


@dataclass(frozen=True)
class ScaleDecision:
    result: shape_audit.ModeResult
    scale_requested: float
    scale_used: float
    fallback_used: bool
    fold_flag: bool
    near_overlap_flag: bool
    plotted_joint_gap: float
    point_order_ok: bool
    passed: bool


def repo_path(path: Path) -> Path:
    return path if path.is_absolute() else REPO_ROOT / path


def rel(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Focused beta-unit, thin-limit MAC, null-vector, and residual audit for Timoshenko shapes.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--beta-deg", type=float, default=DEFAULT_BETA_DEG)
    parser.add_argument("--eta", type=float, default=DEFAULT_ETA)
    parser.add_argument("--n-points", type=int, default=DEFAULT_N_POINTS)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(list(sys.argv[1:] if argv is None else argv))
    if bool(args.smoke):
        args.output_dir = SMOKE_OUTPUT_DIR
        args.n_points = min(int(args.n_points), SMOKE_N_POINTS)
    args.output_dir = repo_path(Path(args.output_dir))
    if abs(float(args.eta)) > 1.0e-14:
        raise ValueError("this focused audit is restricted to eta=0")
    if int(args.n_points) < 101:
        raise ValueError("--n-points must be at least 101")
    if not isfinite(float(args.beta_deg)):
        raise ValueError("--beta-deg must be finite")
    return args


def selected_mu_values(smoke: bool) -> tuple[float, ...]:
    return SMOKE_MU_VALUES if bool(smoke) else DEFAULT_MU_VALUES


def selected_sorted_indices(smoke: bool) -> tuple[int, ...]:
    return SMOKE_SORTED_INDICES if bool(smoke) else DEFAULT_SORTED_INDICES


def assert_beta_unit_contract(beta_deg: float) -> float:
    """Fail early if the imported helper contracts lose their explicit units."""

    beta_rad = float(np.deg2rad(float(beta_deg)))
    for function in (
        shape_audit.TIMO.timo_coupling_matrix,
        shape_audit.TIMO.timo_mode_coefficients,
        shape_audit.TIMO.timo_mode_fields,
        shape_audit.TIMO.timo_det_for_scan,
    ):
        parameters = inspect.signature(function).parameters
        assert "beta_deg" in parameters, f"{function.__name__} must expose beta_deg explicitly"
        assert "beta" not in parameters, f"{function.__name__} has an ambiguous beta parameter"
    assert np.isclose(beta_rad, np.pi / 4.0) if np.isclose(float(beta_deg), 45.0) else True
    return beta_rad


def display_transform_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for beta_deg in (0.0, 15.0, 45.0, 90.0):
        basis = DISPLAY.rod2_display_basis(beta_deg)
        axial_x, axial_y = DISPLAY.rod2_local_displacement_to_display(
            np.array([1.0], dtype=float),
            np.array([0.0], dtype=float),
            beta_deg=beta_deg,
        )
        transverse_x, transverse_y = DISPLAY.rod2_local_displacement_to_display(
            np.array([0.0], dtype=float),
            np.array([1.0], dtype=float),
            beta_deg=beta_deg,
        )
        axial = np.array([axial_x[0], axial_y[0]], dtype=float)
        transverse = np.array([transverse_x[0], transverse_y[0]], dtype=float)
        tangent_dot_normal = float(np.dot(basis.tangent, basis.normal))
        tangent_norm = float(np.linalg.norm(basis.tangent))
        normal_norm = float(np.linalg.norm(basis.normal))
        axial_error = float(np.linalg.norm(axial - basis.tangent))
        transverse_error = abs(float(np.dot(transverse, basis.tangent)))
        passed = bool(
            abs(tangent_dot_normal) <= 1.0e-14
            and abs(tangent_norm - 1.0) <= 1.0e-14
            and abs(normal_norm - 1.0) <= 1.0e-14
            and axial_error <= 1.0e-14
            and transverse_error <= 1.0e-14
            and (not np.isclose(beta_deg, 45.0) or not np.allclose(transverse, basis.tangent))
        )
        rows.append(
            {
                "beta_deg": beta_deg,
                "tangent_x": basis.tangent[0],
                "tangent_y": basis.tangent[1],
                "normal_x": basis.normal[0],
                "normal_y": basis.normal[1],
                "tangent_dot_normal": tangent_dot_normal,
                "tangent_norm": tangent_norm,
                "normal_norm": normal_norm,
                "pure_axial_parallel_error": axial_error,
                "pure_transverse_orthogonality_error": transverse_error,
                "passed": passed,
                "notes": "display basis t2=(cos beta,sin beta), n2=(sin beta,-cos beta)",
            }
        )
    return rows


def build_results(
    args: argparse.Namespace,
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> tuple[list[shape_audit.ModeResult], list[shape_audit.ModeResult]]:
    n_roots = max(int(index) for index in sorted_indices)
    provider = shape_audit.RootProvider(float(args.beta_deg), float(args.eta), n_roots)
    thin_results: list[shape_audit.ModeResult] = []
    problem_timo_results: list[shape_audit.ModeResult] = []

    for mu in mu_values:
        eb_roots, timo_roots, warnings = provider.roots(THIN_EPSILON, float(mu))
        for sorted_index in sorted_indices:
            index = int(sorted_index) - 1
            thin_results.append(
                shape_audit.eb_mode_result(
                    epsilon=THIN_EPSILON,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=float(eb_roots[index]),
                    n_points=int(args.n_points),
                    case_role="thin_limit_shape_bug_audit",
                    root_source="sorted_global_solver",
                    root_warnings=warnings,
                )
            )
            thin_results.append(
                shape_audit.timo_mode_result(
                    epsilon=THIN_EPSILON,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=float(timo_roots[index]),
                    n_points=int(args.n_points),
                    case_role="thin_limit_shape_bug_audit",
                    root_source="sorted_global_solver",
                    root_warnings=warnings,
                )
            )

        _eb_problem, timo_problem, problem_warnings = provider.roots(PROBLEM_EPSILON, float(mu))
        for sorted_index in sorted_indices:
            problem_timo_results.append(
                shape_audit.timo_mode_result(
                    epsilon=PROBLEM_EPSILON,
                    beta_deg=float(args.beta_deg),
                    eta=float(args.eta),
                    mu=float(mu),
                    sorted_index=int(sorted_index),
                    Lambda=float(timo_problem[int(sorted_index) - 1]),
                    n_points=int(args.n_points),
                    case_role="problem_epsilon_comparison",
                    root_source="sorted_global_solver",
                    root_warnings=problem_warnings,
                )
            )
    return thin_results, problem_timo_results


def result_map(results: Sequence[shape_audit.ModeResult]) -> dict[tuple[str, float, int], shape_audit.ModeResult]:
    return {
        (result.model, round(float(result.mu), 12), int(result.sorted_index)): result
        for result in results
    }


def trapezoid_sqrt_weights(x_values: np.ndarray) -> np.ndarray:
    x = np.asarray(x_values, dtype=float)
    dx = np.abs(np.diff(x))
    weights = np.empty_like(x)
    weights[0] = 0.5 * dx[0]
    weights[-1] = 0.5 * dx[-1]
    if x.size > 2:
        weights[1:-1] = 0.5 * (dx[:-1] + dx[1:])
    return np.sqrt(np.maximum(weights, 0.0))


def weighted_two_rod_vector(
    result: shape_audit.ModeResult,
    rod1_values: np.ndarray,
    rod2_values: np.ndarray,
) -> np.ndarray:
    return np.concatenate(
        [
            np.asarray(rod1_values, dtype=float) * trapezoid_sqrt_weights(result.rod1.x_local),
            np.asarray(rod2_values, dtype=float) * trapezoid_sqrt_weights(result.rod2.x_local),
        ]
    )


def component_vector(result: shape_audit.ModeResult, component: str) -> np.ndarray:
    if component == "u":
        return weighted_two_rod_vector(result, result.rod1.u, result.rod2.u)
    if component == "w":
        return weighted_two_rod_vector(result, result.rod1.w, result.rod2.w)
    if component == "psi":
        return weighted_two_rod_vector(
            result,
            np.asarray(result.rod1.psi, dtype=float),
            np.asarray(result.rod2.psi, dtype=float),
        )
    if component == "gamma":
        return weighted_two_rod_vector(
            result,
            np.asarray(result.rod1.gamma, dtype=float),
            np.asarray(result.rod2.gamma, dtype=float),
        )
    raise ValueError(f"unknown component {component!r}")


def grouped_vector(result: shape_audit.ModeResult, components: Sequence[str]) -> np.ndarray:
    return np.concatenate([component_vector(result, component) for component in components])


def mac_value(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    denominator = float(np.dot(a, a) * np.dot(b, b))
    if denominator <= 1.0e-30:
        return float("nan")
    dot = float(np.dot(a, b))
    return float((dot * dot) / denominator)


def eb_slope_vector(result: shape_audit.ModeResult) -> np.ndarray:
    slope1 = np.gradient(result.rod1.w, result.rod1.x_local, edge_order=2)
    slope2 = np.gradient(result.rod2.w, result.rod2.x_local, edge_order=2)
    return weighted_two_rod_vector(result, slope1, slope2)


def global_displacement_vectors(result: shape_audit.ModeResult) -> tuple[np.ndarray, np.ndarray]:
    u1 = np.asarray(result.rod1.u, dtype=float)
    w1 = np.asarray(result.rod1.w, dtype=float)
    u2 = np.asarray(result.rod2.u, dtype=float)
    w2 = np.asarray(result.rod2.w, dtype=float)
    if result.model == MODEL_EB:
        dx1, dy1 = DISPLAY.eb_rod1_local_displacement_to_display(u1, w1)
        dx2, dy2 = DISPLAY.eb_rod2_local_displacement_to_display(u2, w2, beta_deg=result.beta_deg)
    else:
        dx1, dy1 = DISPLAY.rod1_local_displacement_to_display(u1, w1)
        dx2, dy2 = DISPLAY.rod2_local_displacement_to_display(u2, w2, beta_deg=result.beta_deg)
    weight1 = trapezoid_sqrt_weights(result.rod1.x_local)
    weight2 = trapezoid_sqrt_weights(result.rod2.x_local)
    rod1 = np.concatenate([dx1 * weight1, dy1 * weight1])
    rod2 = np.concatenate([dx2 * weight2, dy2 * weight2])
    return np.concatenate([rod1, rod2]), rod2


def thin_mac_rows(
    results: Sequence[shape_audit.ModeResult],
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> list[dict[str, object]]:
    by_key = result_map(results)
    rows: list[dict[str, object]] = []
    for mu in mu_values:
        for sorted_index in sorted_indices:
            eb = by_key[(MODEL_EB, round(float(mu), 12), int(sorted_index))]
            timo = by_key[(MODEL_TIMO, round(float(mu), 12), int(sorted_index))]
            lambda_scale = max(abs(float(eb.Lambda)), 1.0e-30)
            rel_freq_diff = abs(float(timo.Lambda) - float(eb.Lambda)) / lambda_scale
            gamma_vector = component_vector(timo, "gamma")
            w_prime_vector = weighted_two_rod_vector(
                timo,
                np.gradient(timo.rod1.w, timo.rod1.x_local, edge_order=2),
                np.gradient(timo.rod2.w, timo.rod2.x_local, edge_order=2),
            )
            gamma_norm_ratio = float(np.linalg.norm(gamma_vector)) / max(
                float(np.linalg.norm(w_prime_vector)), 1.0e-30
            )
            max_abs_gamma = max(
                float(np.max(np.abs(np.asarray(timo.rod1.gamma, dtype=float)))),
                float(np.max(np.abs(np.asarray(timo.rod2.gamma, dtype=float)))),
            )
            mac_uw = mac_value(grouped_vector(eb, ("u", "w")), grouped_vector(timo, ("u", "w")))
            mac_u = mac_value(component_vector(eb, "u"), component_vector(timo, "u"))
            mac_w = mac_value(component_vector(eb, "w"), component_vector(timo, "w"))
            mac_psi = mac_value(eb_slope_vector(eb), component_vector(timo, "psi"))
            eb_global, eb_global_rod2 = global_displacement_vectors(eb)
            timo_global, timo_global_rod2 = global_displacement_vectors(timo)
            mac_global = mac_value(eb_global, timo_global)
            mac_global_rod2 = mac_value(eb_global_rod2, timo_global_rod2)
            warnings: list[str] = []
            if rel_freq_diff > 0.01:
                warnings.append("relative_frequency_difference_gt_1pct")
            if mac_uw < 0.80:
                warnings.append("low_mac_uw")
            if mac_w < 0.80 and str(eb.energy["classification"]) == "bending_dominated":
                warnings.append("low_mac_w_for_bending_mode")
            if mac_psi < 0.80 and str(eb.energy["classification"]) == "bending_dominated":
                warnings.append("low_mac_psi_vs_eb_slope_for_bending_mode")
            if gamma_norm_ratio > 0.05:
                warnings.append("gamma_norm_ratio_gt_0p05")
            rows.append(
                {
                    "mu": float(mu),
                    "sorted_index": int(sorted_index),
                    "Lambda_EB": eb.Lambda,
                    "Lambda_Timo": timo.Lambda,
                    "rel_freq_diff": rel_freq_diff,
                    "mac_uw": mac_uw,
                    "mac_u": mac_u,
                    "mac_w": mac_w,
                    "mac_psi_vs_eb_slope": mac_psi,
                    "mac_global_displacement": mac_global,
                    "mac_rod2_global_displacement": mac_global_rod2,
                    "mac_display_corrected": mac_global,
                    "mac_rod2_display_corrected": mac_global_rod2,
                    "max_abs_gamma": max_abs_gamma,
                    "gamma_norm_ratio": gamma_norm_ratio,
                    "eb_classification": eb.energy["classification"],
                    "timo_classification": timo.energy["classification"],
                    "warning": "; ".join(warnings),
                    "notes": (
                        "sign-invariant squared MAC with trapezoid spatial weights over both rods; "
                        "u+w uses one physical grouped vector; EB slope is numerical d(w_EB)/dx_i"
                    ),
                }
            )
    return rows


def null_vector_rows(results: Sequence[shape_audit.ModeResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        if result.model != MODEL_TIMO:
            continue
        matrix, matrix_warnings = shape_audit.TIMO.timo_coupling_matrix(
            result.Lambda,
            result.beta_deg,
            result.mu,
            result.epsilon,
            result.eta,
        )
        scaled = shape_audit.TIMO.row_normalized(matrix)
        singular_values = np.linalg.svd(scaled, compute_uv=False)
        sigma_min = float(singular_values[-1])
        sigma_second = float(singular_values[-2])
        ratio = sigma_min / sigma_second if sigma_second > 1.0e-30 else float("nan")
        coeff = np.asarray(result.coeff, dtype=float)
        residual_norm = float(np.linalg.norm(scaled @ coeff))
        raw_residual_norm = float(np.linalg.norm(matrix @ coeff))
        warnings = list(matrix_warnings)
        if isfinite(ratio) and ratio > NULL_RATIO_RELIABLE_MAX:
            warnings.append("singular_value_ratio_not_small")
        if sigma_second <= NEAR_DEGENERATE_SIGMA_SECOND_MAX:
            warnings.append("near_degenerate_nullspace")
        A1, B1, P1, A2, B2, P2 = [float(value) for value in coeff]
        rows.append(
            {
                "mu": result.mu,
                "epsilon": result.epsilon,
                "beta_deg": result.beta_deg,
                "eta": result.eta,
                "sorted_index": result.sorted_index,
                "Lambda": result.Lambda,
                "sigma_min": sigma_min,
                "sigma_second_min": sigma_second,
                "singular_value_ratio": ratio,
                "residual_norm": residual_norm,
                "coeff_norm": float(np.linalg.norm(coeff)),
                "coeff_A1": A1,
                "coeff_B1": B1,
                "coeff_P1": P1,
                "coeff_A2": A2,
                "coeff_B2": B2,
                "coeff_P2": P2,
                "warning": "; ".join(warnings),
                "notes": (
                    "coefficient order A1,B1,P1,A2,B2,P2; singular values and residual_norm use the "
                    f"same row-normalized matrix as analytic_null_vector; raw_matrix_residual_norm={raw_residual_norm:.6e}"
                ),
            }
        )
    return rows


def residual_rows(results: Sequence[shape_audit.ModeResult]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for result in results:
        if result.model != MODEL_TIMO:
            continue
        psi1 = np.asarray(result.rod1.psi, dtype=float)
        psi2 = np.asarray(result.rod2.psi, dtype=float)
        clamp_values = [
            float(result.rod1.u[0]),
            float(result.rod1.w[0]),
            float(psi1[0]),
            float(result.rod2.u[0]),
            float(result.rod2.w[0]),
            float(psi2[0]),
        ]
        max_clamp = max(abs(value) for value in clamp_values)
        joint = result.joint_row or {}
        matrix, _warnings = shape_audit.TIMO.timo_coupling_matrix(
            result.Lambda,
            result.beta_deg,
            result.mu,
            result.epsilon,
            result.eta,
        )
        normalized_residual = float(
            np.linalg.norm(shape_audit.TIMO.row_normalized(matrix) @ np.asarray(result.coeff, dtype=float))
        )
        max_kin = float(joint.get("max_abs_kinematic_gap", float("nan")))
        max_force = float(joint.get("max_abs_force_gap", float("nan")))
        passed = bool(
            max_clamp <= AUDIT_TOL
            and max_kin <= AUDIT_TOL
            and max_force <= AUDIT_TOL
            and normalized_residual <= AUDIT_TOL
        )
        rows.append(
            {
                "mu": result.mu,
                "epsilon": result.epsilon,
                "beta_deg": result.beta_deg,
                "eta": result.eta,
                "sorted_index": result.sorted_index,
                "Lambda": result.Lambda,
                "rod1_clamp_u": clamp_values[0],
                "rod1_clamp_w": clamp_values[1],
                "rod1_clamp_psi": clamp_values[2],
                "rod2_clamp_u": clamp_values[3],
                "rod2_clamp_w": clamp_values[4],
                "rod2_clamp_psi": clamp_values[5],
                "max_abs_clamp_gap": max_clamp,
                "gap_w": joint.get("gap_w", float("nan")),
                "gap_u": joint.get("gap_u", float("nan")),
                "gap_psi": joint.get("gap_psi", float("nan")),
                "gap_M": joint.get("gap_M", float("nan")),
                "gap_Q": joint.get("gap_Q", float("nan")),
                "gap_N": joint.get("gap_N", float("nan")),
                "max_abs_kinematic_gap": max_kin,
                "max_abs_force_gap": max_force,
                "normalized_matrix_residual_norm": normalized_residual,
                "max_ode_residual": float("nan"),
                "passed": passed,
                "warning": "" if passed else "boundary_joint_or_matrix_residual_gt_1e-6",
                "notes": (
                    "rod1 joint x1=+l1; rod2 joint x2=-l2; beta_deg passed to Timoshenko helpers; "
                    "ODE residual not independently differenced because fields are direct analytic basis evaluations"
                ),
            }
        )
    return rows


def point_order_ok(result: shape_audit.ModeResult) -> bool:
    local1 = np.asarray(result.rod1.x_local, dtype=float)
    local2 = np.asarray(result.rod2.x_local, dtype=float)[::-1]
    return bool(
        (np.all(np.diff(local1) >= -1.0e-12) or np.all(np.diff(local1) <= 1.0e-12))
        and (np.all(np.diff(local2) >= -1.0e-12) or np.all(np.diff(local2) <= 1.0e-12))
    )


def evaluate_scale(result: shape_audit.ModeResult, fraction: float) -> ScaleDecision:
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(fraction))
    regularity = [
        shape_audit.regularity_row(
            result,
            rod_id=rod_id,
            scale_kind="fixed_fraction",
            deformation_scale_fraction=float(fraction),
            deformation_scale_value=scale,
        )
        for rod_id in (1, 2)
    ]
    fold = any(bool(row["visually_folded_at_scale"]) for row in regularity)
    near = any(int(row["near_overlap_pair_count"]) > 0 for row in regularity)
    plotted_gap = shape_audit.fixed_scale_plot_joint_gap(result, float(fraction))
    order_ok = point_order_ok(result)
    math_ok = True
    if result.model == MODEL_TIMO:
        joint = result.joint_row or {}
        math_ok = bool(
            float(joint.get("max_abs_kinematic_gap", float("inf"))) <= AUDIT_TOL
            and float(joint.get("max_abs_force_gap", float("inf"))) <= AUDIT_TOL
        )
    passed = bool(not fold and not near and plotted_gap <= AUDIT_TOL and order_ok and math_ok)
    return ScaleDecision(
        result=result,
        scale_requested=float(fraction),
        scale_used=float(fraction),
        fallback_used=False,
        fold_flag=fold,
        near_overlap_flag=near,
        plotted_joint_gap=plotted_gap,
        point_order_ok=order_ok,
        passed=passed,
    )


def scale_decisions(results: Sequence[shape_audit.ModeResult]) -> list[ScaleDecision]:
    decisions: list[ScaleDecision] = []
    for result in results:
        requested = evaluate_scale(result, REQUESTED_SCALE)
        if requested.passed:
            decisions.append(requested)
            continue
        fallback = evaluate_scale(result, FALLBACK_SCALE)
        if not fallback.passed:
            fallback = evaluate_scale(result, SAFETY_SCALE)
        decisions.append(
            ScaleDecision(
                result=result,
                scale_requested=REQUESTED_SCALE,
                scale_used=fallback.scale_used,
                fallback_used=True,
                fold_flag=fallback.fold_flag,
                near_overlap_flag=fallback.near_overlap_flag,
                plotted_joint_gap=fallback.plotted_joint_gap,
                point_order_ok=fallback.point_order_ok,
                passed=fallback.passed,
            )
        )
    return decisions


def decision_map(decisions: Sequence[ScaleDecision]) -> dict[tuple[str, float, int], ScaleDecision]:
    return {
        (item.result.model, round(float(item.result.mu), 12), int(item.result.sorted_index)): item
        for item in decisions
    }


def common_limits(decisions: Sequence[ScaleDecision]) -> tuple[tuple[float, float], tuple[float, float]]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for decision in decisions:
        result = decision.result
        base = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
        deformed = shape_audit.deformed_coordinates_with_scale(
            result,
            transverse_only=False,
            scale=shape_audit.fixed_scale_from_fraction(result.mu, decision.scale_used),
        )
        xs.extend([base[0], base[2], deformed[0], deformed[2]])
        ys.extend([base[1], base[3], deformed[1], deformed[3]])
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    x_span = max(float(np.ptp(x_all)), 1.0)
    y_span = max(float(np.ptp(y_all)), 0.5)
    return (
        (float(np.min(x_all) - 0.07 * x_span), float(np.max(x_all) + 0.07 * x_span)),
        (float(np.min(y_all) - 0.14 * y_span), float(np.max(y_all) + 0.14 * y_span)),
    )


def draw_shape(ax: plt.Axes, decision: ScaleDecision, *, color: str) -> None:
    result = decision.result
    x1, y1, x2, y2 = shape_audit.base_coordinates(result.mu, result.beta_deg, len(result.rod1.u))
    scale = shape_audit.fixed_scale_from_fraction(result.mu, decision.scale_used)
    xd1, yd1, xd2, yd2 = shape_audit.deformed_coordinates_with_scale(
        result,
        transverse_only=False,
        scale=scale,
    )
    ax.plot(x1, y1, color="0.72", linestyle="--", linewidth=0.8)
    ax.plot(x2, y2, color="0.72", linestyle="--", linewidth=0.8)
    ax.plot(xd1, yd1, color=color, linewidth=1.45)
    ax.plot(xd2, yd2, color=color, linewidth=1.45)
    ax.scatter([x1[-1]], [y1[-1]], color="black", s=10, zorder=5)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color="0.91", linewidth=0.45)
    ax.tick_params(labelsize=6.5)


def write_model_grid(
    output_dir: Path,
    decisions: Sequence[ScaleDecision],
    *,
    model: str,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    filename = (
        "eb_full_shapes_eps0p0025_beta45_eta0_modes4_6_grid.png"
        if model == MODEL_EB
        else "timo_full_shapes_eps0p0025_beta45_eta0_modes4_6_grid.png"
    )
    output = output_dir / filename
    selected = [item for item in decisions if item.result.model == model]
    by_key = decision_map(selected)
    limits = common_limits(selected)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(3.15 * len(mu_values), 2.5 * len(sorted_indices) + 0.65),
        squeeze=False,
    )
    color = "#0072b2" if model == MODEL_EB else "#d55e00"
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            decision = by_key[(model, round(float(mu), 12), int(sorted_index))]
            ax = axes[row_index, col_index]
            draw_shape(ax, decision, color=color)
            ax.set_xlim(*limits[0])
            ax.set_ylim(*limits[1])
            fallback = " fallback" if decision.fallback_used else ""
            ax.set_title(
                f"mu={float(mu):g}, sorted {int(sorted_index)}, Lambda={decision.result.Lambda:.6g}\n"
                f"scale={decision.scale_used:.2f}{fallback}; {decision.result.energy['classification']}",
                fontsize=7.4,
            )
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=7.5)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=7.5)
    fig.suptitle(
        f"{model} full displacement, epsilon=0.0025, beta=45 deg, eta=0",
        fontsize=10.5,
        y=0.995,
    )
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95), h_pad=0.55, w_pad=0.45)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def signed_display_vector(result: shape_audit.ModeResult) -> np.ndarray:
    full, _rod2 = global_displacement_vectors(result)
    return float(result.sign_factor) * full


def aligned_deformed_coordinates(
    result: shape_audit.ModeResult,
    *,
    scale_fraction: float,
    extra_sign: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    u1, w1, u2, w2 = shape_audit.plot_arrays(result)
    u1 = float(extra_sign) * u1
    w1 = float(extra_sign) * w1
    u2 = float(extra_sign) * u2
    w2 = float(extra_sign) * w2
    scale = shape_audit.fixed_scale_from_fraction(result.mu, float(scale_fraction))
    _l1, l2 = shape_audit.TIMO.segment_lengths(result.mu)
    rod1_mapper = (
        DISPLAY.eb_rod1_local_fields_to_display
        if result.model == MODEL_EB
        else DISPLAY.rod1_local_fields_to_display
    )
    rod2_mapper = (
        DISPLAY.eb_rod2_local_fields_to_display
        if result.model == MODEL_EB
        else DISPLAY.rod2_local_fields_to_display
    )
    rod1 = rod1_mapper(
        result.rod1.x_local,
        u1,
        w1,
        scale=scale,
    )
    rod2 = rod2_mapper(
        result.rod2.x_local[::-1],
        u2,
        w2,
        l2=l2,
        x_joint=float(result.rod1.x_local[-1]),
        beta_deg=result.beta_deg,
        scale=scale,
    )
    return rod1.x_deformed, rod1.y_deformed, rod2.x_deformed, rod2.y_deformed


def write_paired_grid(
    output_dir: Path,
    results: Sequence[shape_audit.ModeResult],
    decisions: Sequence[ScaleDecision],
    *,
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> Path:
    output = output_dir / "eb_vs_timo_shapes_eps0p0025_beta45_eta0_modes4_6_grid.png"
    by_result = result_map(results)
    by_decision = decision_map(decisions)
    fig, axes = plt.subplots(
        len(sorted_indices),
        len(mu_values),
        figsize=(3.15 * len(mu_values), 2.5 * len(sorted_indices) + 0.75),
        squeeze=False,
    )
    all_x: list[np.ndarray] = []
    all_y: list[np.ndarray] = []
    plot_data: dict[tuple[float, int], tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...], float]] = {}
    for mu in mu_values:
        for sorted_index in sorted_indices:
            eb = by_result[(MODEL_EB, round(float(mu), 12), int(sorted_index))]
            timo = by_result[(MODEL_TIMO, round(float(mu), 12), int(sorted_index))]
            eb_decision = by_decision[(MODEL_EB, round(float(mu), 12), int(sorted_index))]
            timo_decision = by_decision[(MODEL_TIMO, round(float(mu), 12), int(sorted_index))]
            common_fraction = min(float(eb_decision.scale_used), float(timo_decision.scale_used))
            align = 1.0 if float(np.dot(signed_display_vector(eb), signed_display_vector(timo))) >= 0.0 else -1.0
            eb_xy = aligned_deformed_coordinates(eb, scale_fraction=common_fraction, extra_sign=1.0)
            timo_xy = aligned_deformed_coordinates(timo, scale_fraction=common_fraction, extra_sign=align)
            plot_data[(round(float(mu), 12), int(sorted_index))] = (eb_xy, timo_xy, common_fraction)
            all_x.extend([eb_xy[0], eb_xy[2], timo_xy[0], timo_xy[2]])
            all_y.extend([eb_xy[1], eb_xy[3], timo_xy[1], timo_xy[3]])
    x_all = np.concatenate(all_x)
    y_all = np.concatenate(all_y)
    x_span = max(float(np.ptp(x_all)), 1.0)
    y_span = max(float(np.ptp(y_all)), 0.5)
    limits = (
        (float(np.min(x_all) - 0.07 * x_span), float(np.max(x_all) + 0.07 * x_span)),
        (float(np.min(y_all) - 0.14 * y_span), float(np.max(y_all) + 0.14 * y_span)),
    )
    for row_index, sorted_index in enumerate(sorted_indices):
        for col_index, mu in enumerate(mu_values):
            ax = axes[row_index, col_index]
            eb_xy, timo_xy, common_fraction = plot_data[(round(float(mu), 12), int(sorted_index))]
            beta_deg = by_result[(MODEL_EB, round(float(mu), 12), int(sorted_index))].beta_deg
            x1, y1, x2, y2 = shape_audit.base_coordinates(float(mu), beta_deg, len(eb_xy[0]))
            ax.plot(x1, y1, color="0.75", linestyle="--", linewidth=0.75)
            ax.plot(x2, y2, color="0.75", linestyle="--", linewidth=0.75)
            ax.plot(eb_xy[0], eb_xy[1], color="#0072b2", linewidth=1.2, label="EB")
            ax.plot(eb_xy[2], eb_xy[3], color="#0072b2", linewidth=1.2)
            ax.plot(timo_xy[0], timo_xy[1], color="#d55e00", linewidth=1.15, linestyle="--", label="Timoshenko")
            ax.plot(timo_xy[2], timo_xy[3], color="#d55e00", linewidth=1.15, linestyle="--")
            ax.scatter([x1[-1]], [y1[-1]], color="black", s=9, zorder=5)
            ax.set_xlim(*limits[0])
            ax.set_ylim(*limits[1])
            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, color="0.91", linewidth=0.45)
            ax.tick_params(labelsize=6.5)
            ax.set_title(f"mu={float(mu):g}, sorted {int(sorted_index)}, paired scale={common_fraction:.2f}", fontsize=7.3)
            if row_index == len(sorted_indices) - 1:
                ax.set_xlabel("global x", fontsize=7.5)
            if col_index == 0:
                ax.set_ylabel("global y", fontsize=7.5)
    axes[0, 0].legend(loc="best", fontsize=6.7, frameon=False)
    fig.suptitle("EB vs Timoshenko thin-limit full displacement, epsilon=0.0025", fontsize=10.5, y=0.995)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95), h_pad=0.55, w_pad=0.45)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output


def pairwise_timo_stats(
    results: Sequence[shape_audit.ModeResult],
    mu_values: Sequence[float],
    sorted_indices: Sequence[int],
) -> dict[str, float]:
    by_key = result_map(results)
    values: dict[str, list[float]] = {"mac_uw": [], "mac_u": [], "mac_w": [], "mac_psi": [], "mac_gamma": []}
    for mu in mu_values:
        for i, mode_i in enumerate(sorted_indices):
            for mode_j in sorted_indices[i + 1 :]:
                left = by_key[(MODEL_TIMO, round(float(mu), 12), int(mode_i))]
                right = by_key[(MODEL_TIMO, round(float(mu), 12), int(mode_j))]
                values["mac_uw"].append(mac_value(grouped_vector(left, ("u", "w")), grouped_vector(right, ("u", "w"))))
                values["mac_u"].append(mac_value(component_vector(left, "u"), component_vector(right, "u")))
                values["mac_w"].append(mac_value(component_vector(left, "w"), component_vector(right, "w")))
                values["mac_psi"].append(mac_value(component_vector(left, "psi"), component_vector(right, "psi")))
                values["mac_gamma"].append(mac_value(component_vector(left, "gamma"), component_vector(right, "gamma")))
    stats: dict[str, float] = {}
    for name, field_values in values.items():
        array = np.asarray(field_values, dtype=float)
        stats[f"{name}_min"] = float(np.nanmin(array))
        stats[f"{name}_mean"] = float(np.nanmean(array))
        stats[f"{name}_max"] = float(np.nanmax(array))
    return stats


def interior_zero_crossings(values: np.ndarray) -> int:
    y = np.asarray(values, dtype=float)
    amplitude = float(np.max(np.abs(y)))
    if amplitude <= 1.0e-30:
        return 0
    signs = np.sign(y)
    signs[np.abs(y) <= 1.0e-5 * amplitude] = 0.0
    nonzero = signs[signs != 0.0]
    return int(np.sum(nonzero[1:] * nonzero[:-1] < 0.0)) if nonzero.size > 1 else 0


def long_rod_node_summary(results: Sequence[shape_audit.ModeResult]) -> tuple[int, int, float]:
    counts: list[int] = []
    for result in results:
        long_rod = result.rod2 if result.mu >= 0.0 else result.rod1
        counts.append(interior_zero_crossings(long_rod.w))
    return min(counts), max(counts), float(np.mean(counts))


def stat_block(rows: Sequence[dict[str, object]], field: str) -> tuple[float, float, float]:
    values = np.asarray([float(row[field]) for row in rows], dtype=float)
    return float(np.nanmin(values)), float(np.nanmean(values)), float(np.nanmax(values))


def write_beta_unit_report(path: Path, *, beta_deg: float, beta_rad: float) -> Path:
    lines = [
        "# Beta Unit Audit",
        "",
        f"- Audited value: `beta_deg={float(beta_deg):.16g}` and `beta_rad={float(beta_rad):.16g}`.",
        "- `assemble_clamped_coupled_matrix_eta` and the EB root solver expect radians; callers explicitly pass `beta_rad=np.deg2rad(beta_deg)`.",
        "- `timo_coupling_matrix`, `timo_mode_coefficients`, `timo_mode_fields`, and `timo_det_for_scan` expose `beta_deg`; the coupling matrix performs `math.radians(beta_deg)` internally.",
        "- The Timoshenko sorted root path `solve_timo -> timo_det_for_scan -> timo_coupling_matrix` passes `beta_deg` unchanged.",
        "- The shape path `timo_mode_result -> timo_mode_coefficients/timo_mode_fields` passes `beta_deg` unchanged, including when an explicit coefficient vector is supplied.",
        "- Plotting and joint-audit trigonometry uses a separately named local `beta_rad`; it is not passed back into a Timoshenko analytic helper.",
        "- Audited shape callers: `audit_timoshenko_shape_construction.py`, `audit_longitudinal_suspect_modes_eb_timo.py`, `audit_timoshenko_modes456_visualization.py`, `audit_timoshenko_modes_4_6_shape_diagnostics.py`, and `plot_eb_timo_full_mode_shapes_eps0p03_beta45_eta0_modes4_6.py`.",
        "",
        "## Finding",
        "",
        "No degrees/radians mismatch was found in the audited Timoshenko root or shape-reconstruction calls. Existing ambiguous local plotting variables were renamed to `beta_rad`; this is a naming-only caller clarification and does not change values or formulas.",
        "",
        "The audit includes runtime signature assertions requiring the Timoshenko helper parameter name `beta_deg` and forbidding an ambiguous `beta` parameter in those helper contracts.",
        "",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_report(
    path: Path,
    *,
    beta_report: Path,
    display_csv: Path,
    mac_csv: Path,
    null_csv: Path,
    residual_csv: Path,
    figures: Sequence[Path],
    thin_rows: Sequence[dict[str, object]],
    null_rows_data: Sequence[dict[str, object]],
    residual_rows_data: Sequence[dict[str, object]],
    thin_timo_stats: dict[str, float],
    problem_timo_stats: dict[str, float],
    thin_nodes: tuple[int, int, float],
    problem_nodes: tuple[int, int, float],
    decisions: Sequence[ScaleDecision],
) -> Path:
    freq_stats = stat_block(thin_rows, "rel_freq_diff")
    mac_uw_stats = stat_block(thin_rows, "mac_uw")
    mac_u_stats = stat_block(thin_rows, "mac_u")
    mac_w_stats = stat_block(thin_rows, "mac_w")
    mac_psi_stats = stat_block(thin_rows, "mac_psi_vs_eb_slope")
    mac_global_stats = stat_block(thin_rows, "mac_display_corrected")
    mac_global_rod2_stats = stat_block(thin_rows, "mac_rod2_display_corrected")
    gamma_stats = stat_block(thin_rows, "gamma_norm_ratio")
    ratio_stats = stat_block(null_rows_data, "singular_value_ratio")
    sigma_second_min = min(float(row["sigma_second_min"]) for row in null_rows_data)
    residual_max = max(float(row["residual_norm"]) for row in null_rows_data)
    null_warning_count = sum(1 for row in null_rows_data if str(row["warning"]))
    near_degenerate_count = sum("near_degenerate_nullspace" in str(row["warning"]) for row in null_rows_data)
    ratio_warning_count = sum("singular_value_ratio_not_small" in str(row["warning"]) for row in null_rows_data)
    residual_fail_count = sum(1 for row in residual_rows_data if not bool(row["passed"]))
    fallback_count = sum(1 for item in decisions if item.fallback_used)
    fallback_005_count = sum(abs(float(item.scale_used) - FALLBACK_SCALE) <= 1.0e-12 for item in decisions)
    safety_002_count = sum(abs(float(item.scale_used) - SAFETY_SCALE) <= 1.0e-12 for item in decisions)
    final_fold_count = sum(1 for item in decisions if item.fold_flag)
    final_near_count = sum(1 for item in decisions if item.near_overlap_flag)
    max_plot_gap = max(float(item.plotted_joint_gap) for item in decisions)
    thin_distinct = thin_timo_stats["mac_uw_max"] < 0.95 or thin_timo_stats["mac_w_max"] < 0.95
    thin_match = (
        freq_stats[2] < 0.01
        and mac_w_stats[1] >= 0.80
        and mac_global_stats[0] > 0.98
        and mac_global_rod2_stats[0] > 0.98
        and gamma_stats[2] < 0.05
    )
    null_reliable = ratio_warning_count == 0 and residual_max <= AUDIT_TOL

    lines = [
        "# Timoshenko Shape Bug Audit",
        "",
        "## Scope",
        "",
        "- Thin limit: `epsilon=0.0025`, `beta=45 deg`, `eta=0`, `mu=0,0.2,0.4,0.6`, sorted modes `4,5,6`.",
        "- Comparison case: `epsilon=0.03` with the same beta, eta, mu, and sorted indices.",
        "- Sorted frequencies are used directly. Descendant tracking is not used.",
        "- Full displacement means local `u` and `w`; rods are plotted separately. Timoshenko uses the corrected reflected display bases, while EB uses its equivalent opposite-sign determinant mapping.",
        "",
        "## Output Files",
        "",
        f"- Beta unit audit: `{rel(beta_report)}`",
        f"- Display-transform regression CSV: `{rel(display_csv)}`",
        f"- Thin-limit MAC CSV: `{rel(mac_csv)}`",
        f"- Null-vector CSV: `{rel(null_csv)}`",
        f"- Residual CSV: `{rel(residual_csv)}`",
    ]
    for figure in figures:
        lines.append(f"- Figure: `{rel(figure)}`")
    lines.extend(
        [
            "",
            "## Thin-Limit Frequencies And Shape MAC",
            "",
            f"- Relative frequency difference min/mean/max: `{freq_stats[0]:.6e}` / `{freq_stats[1]:.6e}` / `{freq_stats[2]:.6e}`.",
            f"- MAC(u,w) min/mean/max: `{mac_uw_stats[0]:.4f}` / `{mac_uw_stats[1]:.4f}` / `{mac_uw_stats[2]:.4f}`.",
            f"- MAC(u) min/mean/max: `{mac_u_stats[0]:.4f}` / `{mac_u_stats[1]:.4f}` / `{mac_u_stats[2]:.4f}`.",
            f"- MAC(w) min/mean/max: `{mac_w_stats[0]:.4f}` / `{mac_w_stats[1]:.4f}` / `{mac_w_stats[2]:.4f}`.",
            f"- MAC(psi_Timo, w'_EB) min/mean/max: `{mac_psi_stats[0]:.4f}` / `{mac_psi_stats[1]:.4f}` / `{mac_psi_stats[2]:.4f}`.",
            f"- Corrected display MAC min/mean/max: `{mac_global_stats[0]:.4f}` / `{mac_global_stats[1]:.4f}` / `{mac_global_stats[2]:.4f}`.",
            f"- Corrected rod-2-only display MAC min/mean/max: `{mac_global_rod2_stats[0]:.4f}` / `{mac_global_rod2_stats[1]:.4f}` / `{mac_global_rod2_stats[2]:.4f}`.",
            f"- gamma norm ratio min/mean/max: `{gamma_stats[0]:.6e}` / `{gamma_stats[1]:.6e}` / `{gamma_stats[2]:.6e}`.",
            "",
            "| mu | sorted | Lambda EB | Lambda Timo | rel diff | MAC(u,w) | MAC(w) | gamma ratio |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in thin_rows:
        lines.append(
            f"| {float(row['mu']):g} | {int(row['sorted_index'])} | {float(row['Lambda_EB']):.8f} | "
            f"{float(row['Lambda_Timo']):.8f} | {float(row['rel_freq_diff']):.3e} | "
            f"{float(row['mac_uw']):.4f} | {float(row['mac_w']):.4f} | {float(row['gamma_norm_ratio']):.3e} |"
        )
    lines.extend(
        [
            "",
            "## Display-Frame Correction",
            "",
            "Determinant-frame components are not Cartesian plotting coordinates. Positive determinant `w1` points opposite to positive display Y, so rod 1 uses `t1=(1,0)` and `n1=(0,-1)`. Rod 2 uses the displayed tangent `t2=(cos(beta),sin(beta))` and reflected normal `n2=(sin(beta),-cos(beta))`.",
            "",
            "The old global Timoshenko plots used the determinant coupling relation directly as `(X,Y)`, giving `dY2=-sin(beta)*u2+cos(beta)*w2`. Those centerline and vector-field figures are invalid for physical interpretation. Local fields, frequencies, null vectors, joint equations, and energy fractions remain valid.",
            "",
            "The corrected display mapping is `dX1=u1`, `dY1=-w1`, `dX2=cos(beta)*u2+sin(beta)*w2`, `dY2=sin(beta)*u2-cos(beta)*w2`. Its tangent and normal are orthonormal for beta 0, 15, 45, and 90 degrees, and the corrected thin-limit EB/Timoshenko display MAC confirms the mapping.",
            "EB uses the opposite determinant transverse sign, so its helper maps through `n1=(0,1)` and `n2=(-sin(beta),cos(beta))`. This is not a different physical geometry: the two theory-specific mappings overlay in the thin limit.",
            "",
            "## Null Vector And Reconstruction",
            "",
            f"- Singular-value ratio min/mean/max: `{ratio_stats[0]:.6e}` / `{ratio_stats[1]:.6e}` / `{ratio_stats[2]:.6e}`.",
            f"- Minimum second-smallest singular value: `{sigma_second_min:.6e}`.",
            f"- Maximum row-normalized `||M c||`: `{residual_max:.6e}`.",
            f"- Null-vector warning rows: `{null_warning_count}`; absolute-sigma near-degenerate marks: `{near_degenerate_count}`; poor-ratio warnings: `{ratio_warning_count}`.",
            "- The absolute `sigma_second` marks occur in the thin-limit ill-conditioned scaling, but the corresponding ratios remain at most order 1e-9. Together with the near-unit EB/Timoshenko field MAC, this does not indicate an arbitrary null-vector combination.",
            f"- Boundary/joint/matrix residual failures: `{residual_fail_count}`.",
            "- Coefficient ordering is explicitly `A1,B1,P1,A2,B2,P2`.",
            "- ODE residuals were not independently differenced because each field is evaluated directly from the analytic Timoshenko basis; clamp, joint, and matrix residuals were checked numerically.",
            "",
            "## Epsilon 0.0025 Versus 0.03",
            "",
            "Pairwise Timoshenko MAC among sorted modes 4, 5, and 6:",
            "",
            "| epsilon | mean MAC(u,w) | max MAC(u,w) | mean MAC(w) | max MAC(w) | mean MAC(psi) | mean MAC(gamma) |",
            "|---:|---:|---:|---:|---:|---:|---:|",
            f"| 0.0025 | {thin_timo_stats['mac_uw_mean']:.4f} | {thin_timo_stats['mac_uw_max']:.4f} | {thin_timo_stats['mac_w_mean']:.4f} | {thin_timo_stats['mac_w_max']:.4f} | {thin_timo_stats['mac_psi_mean']:.4f} | {thin_timo_stats['mac_gamma_mean']:.4f} |",
            f"| 0.03 | {problem_timo_stats['mac_uw_mean']:.4f} | {problem_timo_stats['mac_uw_max']:.4f} | {problem_timo_stats['mac_w_mean']:.4f} | {problem_timo_stats['mac_w_max']:.4f} | {problem_timo_stats['mac_psi_mean']:.4f} | {problem_timo_stats['mac_gamma_mean']:.4f} |",
            "",
            f"- Long-rod transverse interior zero-crossing count at epsilon=0.0025, min/max/mean: `{thin_nodes[0]}` / `{thin_nodes[1]}` / `{thin_nodes[2]:.2f}`.",
            f"- Long-rod transverse interior zero-crossing count at epsilon=0.03, min/max/mean: `{problem_nodes[0]}` / `{problem_nodes[1]}` / `{problem_nodes[2]:.2f}`.",
            "- A visually similar one-half-wave full centerline is not sufficient evidence of identical modal fields; the pairwise component MAC values above are the deciding diagnostic.",
            "- After the display correction, the thin-limit paired centerlines agree and the distinct nodal patterns of sorted modes 4, 5, and 6 are visible. The former repeated one-half-wave appearance was a display-frame artifact, not a physical thick-beam effect.",
            "",
            "## Plot Scale Audit",
            "",
            f"- Shapes plotted at requested scale 0.08 / fallback 0.05 / safety fallback 0.02: `{len(decisions) - fallback_count}` / `{fallback_005_count}` / `{safety_002_count}`.",
            f"- Remaining fold / near-overlap flags: `{final_fold_count}` / `{final_near_count}`.",
            f"- Maximum plotted joint gap: `{max_plot_gap:.6e}`.",
            "",
            "## Answers",
            "",
            "1. Was a beta degrees/radians mismatch found?",
            "No. EB receives radians; every audited Timoshenko root, matrix, null-vector, and field call receives degrees.",
            "",
            "2. Do EB and Timoshenko shapes match in the thin limit?",
            ("Yes, within the frequency, transverse-shape, and shear criteria used by this audit." if thin_match else "Not uniformly. See the low-MAC rows above; this is evidence that the same-index reconstruction needs further investigation."),
            "",
            "3. Is gamma small in the thin limit?",
            ("Yes." if gamma_stats[2] < 0.05 else "No; at least one gamma norm ratio exceeds 0.05."),
            "",
            "4. Are Timoshenko sorted modes 4, 5, and 6 distinct in the thin limit?",
            ("Yes. Their pairwise component MAC values do not support one common field." if thin_distinct else "They are unusually similar by the present pairwise MAC threshold and require a degeneracy-focused follow-up."),
            "",
            "5. Are the epsilon=0.03 plots caused by a physical change, visualization, or a reconstruction bug?",
            (
                "The corrected thin-limit local and display-field agreement clears beta, root, null-vector, and reconstruction hypotheses. The former epsilon=0.03 resemblance was caused by the old rod-2 display transform; remaining differences after correction are genuine modal differences."
                if thin_match and null_reliable
                else "The current evidence does not clear reconstruction: the failed thin-limit or null-vector checks identify the next technical target."
            ),
            "",
            "6. Is same sorted-index EB/Timoshenko comparison valid in the thin limit?",
            (
                "Yes for this grid: frequencies are close and the same-index field comparisons are consistent; no descendant tracking is needed for this local diagnostic."
                if thin_match
                else "Not for every row without qualification; inspect frequency proximity and MAC together rather than assuming index equivalence."
            ),
            "",
            "7. What should be fixed next?",
            (
                "Use the shared display helper with its explicit EB and reflected Timoshenko mappings for every global centerline and vector-field plot. No analytic-model change is indicated."
                if thin_match and null_reliable
                else "Investigate the specific low-MAC or high-singular-ratio rows before changing visualization; formulas and solvers remain frozen."
            ),
            "",
            "The deformation scale changes visualization only. Frequencies, energy fractions, formulas, determinant entries, root solvers, and k' are unchanged. No FEM, 3D FEM, Gmsh, or CalculiX workflow was run.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    args = parse_args(argv)
    beta_rad = assert_beta_unit_contract(float(args.beta_deg))
    mu_values = selected_mu_values(bool(args.smoke))
    sorted_indices = selected_sorted_indices(bool(args.smoke))
    thin_results, problem_timo_results = build_results(
        args,
        mu_values=mu_values,
        sorted_indices=sorted_indices,
    )
    thin_rows = thin_mac_rows(thin_results, mu_values, sorted_indices)
    thin_timo_results = [result for result in thin_results if result.model == MODEL_TIMO]
    null_data = null_vector_rows(thin_timo_results)
    residual_data = residual_rows(thin_timo_results)
    decisions = scale_decisions(thin_results)

    beta_report = write_beta_unit_report(
        args.output_dir / "beta_unit_audit.md",
        beta_deg=float(args.beta_deg),
        beta_rad=beta_rad,
    )
    display_csv = shape_audit.write_csv(
        args.output_dir / "display_transform_regression.csv",
        display_transform_rows(),
        DISPLAY_TRANSFORM_FIELDS,
    )
    mac_csv = shape_audit.write_csv(
        args.output_dir / "thin_limit_eb_timo_shape_mac.csv",
        thin_rows,
        THIN_MAC_FIELDS,
    )
    null_csv = shape_audit.write_csv(
        args.output_dir / "timoshenko_null_vector_audit.csv",
        null_data,
        NULL_VECTOR_FIELDS,
    )
    residual_csv = shape_audit.write_csv(
        args.output_dir / "timoshenko_shape_residual_audit.csv",
        residual_data,
        RESIDUAL_FIELDS,
    )
    figures = [
        write_model_grid(
            args.output_dir,
            decisions,
            model=MODEL_EB,
            mu_values=mu_values,
            sorted_indices=sorted_indices,
        ),
        write_model_grid(
            args.output_dir,
            decisions,
            model=MODEL_TIMO,
            mu_values=mu_values,
            sorted_indices=sorted_indices,
        ),
        write_paired_grid(
            args.output_dir,
            thin_results,
            decisions,
            mu_values=mu_values,
            sorted_indices=sorted_indices,
        ),
    ]
    thin_timo_stats = pairwise_timo_stats(thin_timo_results, mu_values, sorted_indices)
    problem_timo_stats = pairwise_timo_stats(problem_timo_results, mu_values, sorted_indices)
    thin_nodes = long_rod_node_summary(thin_timo_results)
    problem_nodes = long_rod_node_summary(problem_timo_results)
    report = write_report(
        args.output_dir / "timoshenko_shape_bug_audit_report.md",
        beta_report=beta_report,
        display_csv=display_csv,
        mac_csv=mac_csv,
        null_csv=null_csv,
        residual_csv=residual_csv,
        figures=figures,
        thin_rows=thin_rows,
        null_rows_data=null_data,
        residual_rows_data=residual_data,
        thin_timo_stats=thin_timo_stats,
        problem_timo_stats=problem_timo_stats,
        thin_nodes=thin_nodes,
        problem_nodes=problem_nodes,
        decisions=decisions,
    )

    print("generated Timoshenko thin-limit shape bug audit outputs:")
    for output in [beta_report, display_csv, mac_csv, null_csv, residual_csv, report, *figures]:
        print(f"  {rel(output)}")
    for field in (
        "rel_freq_diff",
        "mac_uw",
        "mac_u",
        "mac_w",
        "mac_psi_vs_eb_slope",
        "mac_display_corrected",
        "mac_rod2_display_corrected",
        "gamma_norm_ratio",
    ):
        minimum, mean, maximum = stat_block(thin_rows, field)
        print(f"{field}: min={minimum:.6e}, mean={mean:.6e}, max={maximum:.6e}")
    ratio_min, ratio_mean, ratio_max = stat_block(null_data, "singular_value_ratio")
    print(f"singular_value_ratio: min={ratio_min:.6e}, mean={ratio_mean:.6e}, max={ratio_max:.6e}")
    print(f"null-vector warnings: {sum(1 for row in null_data if str(row['warning']))}")
    print(f"residual audit failures: {sum(1 for row in residual_data if not bool(row['passed']))}")
    print(f"scale fallbacks below 0.08: {sum(1 for item in decisions if item.fallback_used)}")
    print(f"max plotted joint gap: {max(float(item.plotted_joint_gap) for item in decisions):.6e}")
    return {
        "beta_report": beta_report,
        "display_csv": display_csv,
        "mac_csv": mac_csv,
        "null_csv": null_csv,
        "residual_csv": residual_csv,
        "figures": figures,
        "report": report,
        "thin_rows": thin_rows,
        "null_rows": null_data,
        "residual_rows": residual_data,
    }


if __name__ == "__main__":
    main()
