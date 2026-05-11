"""Compare analytic and FEM joint dynamic stiffness operators.

This diagnostic condenses each clamped arm to the coupled joint and compares
the resulting dynamic stiffness matrices at a tracked analytic Lambda.  The FEM
side uses the same element matrices and right-arm rotation convention as
``python_fem.py``.  The analytic side uses full cos/sin/cosh/sinh bending
solutions plus sin/cos axial solutions; it does not modify the determinant or
baseline FEM model.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from typing import Sequence

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from my_project.analytic.formulas import assemble_clamped_coupled_matrix  # noqa: E402
from my_project.fem import python_fem as fem  # noqa: E402
from scripts.analysis.check_analytic_shape_in_fem_residual import (  # noqa: E402
    DEFAULT_NUM_ANALYTIC_ROOTS,
    analytic_local_fields,
    assemble_fem_matrices,
    build_analytic_global_q,
    build_params_for_case,
    residual_metrics,
    resolve_analytic_branch,
    track_fem_branch,
)
from scripts.compare_beta0_analytic_vs_fem import fem_parameter_override  # noqa: E402
from scripts.lib.analytic_coupled_rods_shapes import (  # noqa: E402
    NEAR_ZERO_NORM,
    ROW_LABELS,
    analytic_null_vector,
)


DEFAULT_BRANCH_ID = "bending_desc_05"
DEFAULT_BETA = 15.0
DEFAULT_MU = 0.2
DEFAULT_EPSILON = 0.0025
DEFAULT_L_TOTAL = 2.0
DEFAULT_OUTPUT_PREFIX = None

DOF_LABELS = ("Fx", "Fy", "M")
JOINT_Q_LABELS = ("ux", "uy", "theta")

MATRIX_FIELDNAMES = [
    "branch_id",
    "beta",
    "mu",
    "epsilon",
    "analytic_lambda",
    "omega_sq",
    "matrix_name",
    "row",
    "col",
    "fem_value",
    "analytic_value",
    "abs_diff",
    "rel_diff",
]

SUMMARY_FIELDNAMES = [
    "branch_id",
    "beta",
    "mu",
    "epsilon",
    "analytic_lambda",
    "omega_sq",
    "current_sorted_index",
    "branch_tracking_warning_flag",
    "norm_D_left_diff",
    "rel_norm_D_left_diff",
    "norm_D_right_local_diff",
    "rel_norm_D_right_local_diff",
    "norm_D_right_global_diff",
    "rel_norm_D_right_global_diff",
    "norm_D_total_diff",
    "rel_norm_D_total_diff",
    "largest_entry_diff",
    "largest_entry_matrix",
    "largest_entry_row",
    "largest_entry_col",
    "largest_row_family",
    "determinant_dynamic_scaled_max_abs_diff",
    "determinant_dynamic_scaled_norm_diff",
    "q_joint_ux",
    "q_joint_uy",
    "q_joint_theta",
]

CROSSCHECK_FIELDNAMES = [
    "branch_id",
    "beta",
    "mu",
    "epsilon",
    "analytic_lambda",
    "row_index",
    "row_label",
    "determinant_matrix_value",
    "dynamic_total_scaled_value",
    "dynamic_total_fem_scaled_value",
    "abs_diff",
    "rel_diff",
    "scale_note",
]

ROW_AUDIT_FIELDNAMES = [
    "case_id",
    "matrix_pair",
    "row_index",
    "row_label",
    "col_index",
    "col_label",
    "current_value",
    "comparison_value",
    "abs_diff",
    "rel_diff",
    "force_order_variant",
    "sign_variant",
    "notes",
]

ROW_SCALING_AUDIT_FIELDNAMES = [
    "case_id",
    "matrix_pair",
    "row_index",
    "row_label",
    "force_order_variant",
    "sign_variant",
    "best_scalar",
    "row_l2_relative_diff_before_scaling",
    "row_l2_relative_diff_after_best_scalar",
    "row_correlation",
    "notes",
]

FORCE_ROW_CONTRIBUTION_FIELDNAMES = [
    "case_id",
    "row_label",
    "contribution_label",
    "col_label",
    "current_value",
    "independent_endpoint_value",
    "dynamic_stiffness_value",
    "abs_diff_current_vs_endpoint",
    "abs_diff_current_vs_dynamic",
    "rel_diff_current_vs_endpoint",
    "rel_diff_current_vs_dynamic",
    "inferred_missing_factor",
    "notes",
]

NULL_VECTOR_ROW_AUDIT_FIELDNAMES = [
    "case_id",
    "row_label",
    "current_residual",
    "independent_endpoint_residual",
    "dynamic_stiffness_residual",
    "current_vs_endpoint_abs_diff",
    "current_vs_dynamic_abs_diff",
    "notes",
]

CONSTRAINED_ROW_AUDIT_FIELDNAMES = [
    "case_id",
    "variant_name",
    "subspace_dim",
    "row_label",
    "row_index",
    "subspace_col",
    "determinant_value",
    "dynamic_value",
    "abs_diff",
    "rel_diff",
    "row_best_scalar",
    "row_error_after_scaling",
    "notes",
]

CONSTRAINED_SUMMARY_FIELDNAMES = [
    "case_id",
    "variant_name",
    "subspace_dim",
    "Rk_singular_values",
    "Rk_rank",
    "Rk_nullspace_tol",
    "Rk_condition_number",
    "kinematic_nullspace_residual_norm",
    "q_left_right_subspace_mismatch",
    "frobenius_rel_diff_before_scaling",
    "frobenius_rel_diff_after_row_scaling",
    "max_row_error_after_scaling",
    "best_row_scalars",
    "row_correlations",
    "conclusion_flag",
    "notes",
]

CONSTRAINED_NULL_VECTOR_FIELDNAMES = [
    "case_id",
    "variant_name",
    "row_label",
    "row_index",
    "determinant_c0",
    "dynamic_c0",
    "determinant_projected",
    "dynamic_projected",
    "det_vs_dynamic_c0_abs_diff",
    "det_vs_dynamic_projected_abs_diff",
    "c0_projection_abs_error",
    "c0_projection_relative_error",
    "kinematic_residual_norm",
    "notes",
]

FORCE_ROW_SCALING_AUDIT_FIELDNAMES = [
    "case_id",
    "row_label",
    "contribution_label",
    "col_label",
    "determinant_current_value",
    "fem_compatible_value",
    "ratio_det_to_fem",
    "candidate_factor_match",
    "abs_diff",
    "rel_diff",
    "notes",
]

FORCE_ROW_SCALING_SUMMARY_FIELDNAMES = [
    "case_id",
    "row_label",
    "contribution_label",
    "ratio_pattern",
    "likely_missing_factor",
    "supports_scaling_mismatch",
    "nonzero_ratio_count",
    "ratio_min",
    "ratio_max",
    "notes",
]

FORCE_ROW_SCALING_DIAGNOSTIC_FIELDNAMES = [
    "case_id",
    "diagnostic_variant",
    "current_null_current_force_norm",
    "diagnostic_force_on_current_null_norm",
    "diagnostic_new_null_residual_norm",
    "diagnostic_new_null_smallest_singular_value",
    "diagnostic_new_null_analytic_in_fem_relative_residual",
    "diagnostic_new_null_analytic_in_fem_residual_l2",
    "notes",
]

JOINT_NULLSPACE_MATRIX_FIELDNAMES = [
    "case_id",
    "matrix_name",
    "row",
    "col",
    "fem_value",
    "analytic_value",
    "abs_diff",
    "rel_diff",
]

JOINT_NULLSPACE_SUMMARY_FIELDNAMES = [
    "case_id",
    "lambda",
    "fem_singular_values",
    "analytic_singular_values",
    "fem_condition_number",
    "analytic_condition_number",
    "frobenius_rel_diff",
    "max_abs_diff",
    "smallest_singular_value_fem",
    "smallest_singular_value_analytic",
]

JOINT_NULLSPACE_VECTOR_FIELDNAMES = [
    "case_id",
    "vector_name",
    "ux",
    "uy",
    "theta",
    "normalized_ux",
    "normalized_uy",
    "normalized_theta",
]

JOINT_NULLSPACE_PAIRWISE_FIELDNAMES = [
    "case_id",
    "vector_a",
    "vector_b",
    "mac",
    "relative_l2_after_sign_alignment",
    "dot_after_sign_alignment",
    "notes",
]

JOINT_RECONSTRUCTION_RESIDUAL_FIELDNAMES = [
    "case_id",
    "reconstruction_name",
    "lambda",
    "omega_sq",
    "residual_l2",
    "relative_residual",
    "residual_inf",
    "relative_residual_inf",
    "q_joint_ux",
    "q_joint_uy",
    "q_joint_theta",
    "notes",
]

DETERMINANT_NULLSPACE_CONDITIONING_FIELDNAMES = [
    "method",
    "row_scaling",
    "column_scaling",
    "smallest_singular_value",
    "second_singular_value",
    "singular_ratio",
    "raw_matrix_residual",
    "row_normalized_matrix_residual",
    "q_joint_mac_vs_joint_null",
    "q_joint_rel_l2_vs_joint_null",
    "fem_relative_residual",
    "rayleigh_lambda",
    "notes",
]

DETERMINANT_AS_JOINT_MATRIX_FIELDNAMES = [
    "row",
    "col",
    "D_det_joint_value",
    "D_analytic_joint_value",
    "D_fem_joint_value",
    "abs_diff_det_vs_analytic",
    "rel_diff_det_vs_analytic",
    "abs_diff_det_vs_fem",
    "rel_diff_det_vs_fem",
]

DETERMINANT_AS_JOINT_SUMMARY_FIELDNAMES = [
    "case_id",
    "det_joint_singular_values",
    "analytic_joint_singular_values",
    "fem_joint_singular_values",
    "frobenius_rel_diff_det_vs_analytic",
    "frobenius_rel_diff_det_vs_fem",
    "smallest_singular_vector_MAC_det_vs_analytic",
    "smallest_singular_vector_MAC_det_vs_fem",
    "condition_QN",
    "notes",
]

DETERMINANT_AS_JOINT_VARIANT_FIELDNAMES = [
    "case_id",
    "variant_name",
    "row",
    "col",
    "D_det_joint_variant_value",
    "D_analytic_joint_value",
    "abs_diff_det_vs_analytic",
    "rel_diff_det_vs_analytic",
    "frobenius_rel_diff_det_vs_analytic",
    "smallest_singular_vector_MAC_det_vs_analytic",
    "notes",
]

BASIS_LABELS = ("A1", "B1", "A2", "B2", "P1", "P2")


def filename_number_token(value: float) -> str:
    return f"{float(value):g}".replace("-", "m").replace(".", "p")


def default_output_prefix(branch_id: str, beta: float, mu: float, epsilon: float) -> Path:
    return (
        REPO_ROOT
        / "results"
        / (
            f"joint_dynamic_stiffness_compare_beta{filename_number_token(beta)}_{branch_id}"
            f"_mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}"
        )
    )


def joint_nullspace_output_prefix(branch_id: str, beta: float, mu: float, epsilon: float) -> Path:
    return (
        REPO_ROOT
        / "results"
        / (
            f"joint_nullspace_compare_beta{filename_number_token(beta)}_{branch_id}"
            f"_mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}"
        )
    )


def resolve_output_prefix(value: str | None, branch_id: str, beta: float, mu: float, epsilon: float) -> Path:
    if value is None:
        return default_output_prefix(branch_id, beta, mu, epsilon)
    path = Path(value)
    return path if path.is_absolute() else REPO_ROOT / path


def suffixed_path(prefix: Path, suffix: str) -> Path:
    return Path(f"{prefix}{suffix}")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        allow_abbrev=False,
        description="Compare analytic and FEM Schur-complement dynamic stiffness at the coupled joint.",
    )
    parser.add_argument("--branch-id", default=DEFAULT_BRANCH_ID)
    parser.add_argument("--beta", type=float, default=DEFAULT_BETA)
    parser.add_argument("--mu", type=float, default=DEFAULT_MU)
    parser.add_argument("--epsilon", type=float, default=DEFAULT_EPSILON)
    parser.add_argument("--l-total", type=float, default=DEFAULT_L_TOTAL)
    parser.add_argument("--num-analytic-roots", type=int, default=DEFAULT_NUM_ANALYTIC_ROOTS)
    parser.add_argument("--output-prefix", default=DEFAULT_OUTPUT_PREFIX)
    parser.add_argument(
        "--row-audit",
        "--audit-determinant-rows",
        dest="row_audit",
        action="store_true",
        help="Write row-level determinant force-row audit CSVs without modifying formulas.",
    )
    parser.add_argument(
        "--constrained-row-audit",
        action="store_true",
        help=(
            "Compare determinant and dynamic-stiffness force rows only on the "
            "kinematic-compatibility nullspace."
        ),
    )
    parser.add_argument(
        "--audit-force-row-scaling",
        action="store_true",
        help="Audit determinant force-row scaling against FEM-compatible N/Q/M endpoint quantities.",
    )
    parser.add_argument(
        "--compare-joint-nullspace",
        action="store_true",
        default=False,
        help=(
            "Compare D_total_FEM and D_total_analytic through the physical "
            "joint displacement nullspace q_joint=[ux,uy,theta]."
        ),
    )
    parser.add_argument(
        "--audit-determinant-nullspace-conditioning",
        action="store_true",
        default=False,
        help=(
            "Audit whether the determinant coefficient null vector is sensitive "
            "to row/column scaling, without changing the determinant."
        ),
    )
    parser.add_argument(
        "--derive-determinant-as-joint-system",
        action="store_true",
        default=False,
        help=(
            "Restrict determinant force rows to the kinematic nullspace and "
            "express them as a 3x3 q_joint=[ux,uy,theta] system."
        ),
    )
    args = parser.parse_args(argv)
    if args.epsilon <= 0.0:
        parser.error("--epsilon must be positive.")
    if args.l_total <= 0.0:
        parser.error("--l-total must be positive.")
    if args.mu < -1e-12:
        parser.error("--mu must be non-negative for canonical branch tracking.")
    return args


def safe_relative_diff(value: float, reference: float) -> float:
    return float(abs(float(value) - float(reference)) / abs(float(reference))) if abs(float(reference)) > NEAR_ZERO_NORM else np.nan


def relative_matrix_norm(diff: np.ndarray, reference: np.ndarray) -> float:
    denom = float(np.linalg.norm(reference))
    return float(np.linalg.norm(diff) / denom) if denom > NEAR_ZERO_NORM else np.nan


def rotation3(beta_deg: float) -> np.ndarray:
    beta_rad = float(np.deg2rad(beta_deg))
    c, s = float(np.cos(beta_rad)), float(np.sin(beta_rad))
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)


def row_mapping_note() -> str:
    return (
        "q_joint=[ux,uy,theta_fem=dw/ds]; f_joint=[Fx,Fy,M]; "
        "determinant force rows assumed [M/Lambda^2, eps*Fy/Lambda^2, eps*Fx/Lambda^2]"
    )


def assemble_fem_arm_matrices(
    *,
    length: float,
    beta_deg: float,
    rotate_to_global: bool,
) -> tuple[np.ndarray, np.ndarray]:
    n = fem.N_ELEM
    ndof = 3 * (n + 1)
    le = float(length) / n
    K = np.zeros((ndof, ndof), dtype=float)
    M = np.zeros((ndof, ndof), dtype=float)
    Ke = fem.elem_K(le)
    Me = fem.elem_M(le)
    if rotate_to_global:
        T = fem.rotation_matrix_6x6(np.deg2rad(float(beta_deg)))
        Ke = T @ Ke @ T.T
        Me = T @ Me @ T.T

    for elem in range(n):
        dofs = np.array(
            [
                3 * elem,
                3 * elem + 1,
                3 * elem + 2,
                3 * (elem + 1),
                3 * (elem + 1) + 1,
                3 * (elem + 1) + 2,
            ],
            dtype=int,
        )
        K[np.ix_(dofs, dofs)] += Ke
        M[np.ix_(dofs, dofs)] += Me
    return K, M


def schur_joint_dynamic_stiffness(
    K: np.ndarray,
    M: np.ndarray,
    *,
    omega_sq: float,
    joint_node: int,
    clamp_node: int,
) -> np.ndarray:
    A = np.asarray(K, dtype=float) - float(omega_sq) * np.asarray(M, dtype=float)
    n_nodes = A.shape[0] // 3
    joint = np.array([3 * joint_node, 3 * joint_node + 1, 3 * joint_node + 2], dtype=int)
    clamp = {3 * clamp_node, 3 * clamp_node + 1, 3 * clamp_node + 2}
    joint_set = set(int(idx) for idx in joint)
    internal = np.array(
        [idx for idx in range(3 * n_nodes) if idx not in clamp and idx not in joint_set],
        dtype=int,
    )
    Ajj = A[np.ix_(joint, joint)]
    if internal.size == 0:
        return Ajj
    Aji = A[np.ix_(joint, internal)]
    Aij = A[np.ix_(internal, joint)]
    Aii = A[np.ix_(internal, internal)]
    return Ajj - Aji @ np.linalg.solve(Aii, Aij)


def full_basis_rows(k: float, x: float) -> dict[str, np.ndarray]:
    z = float(k) * float(x)
    return {
        "w": np.array([np.cos(z), np.sin(z), np.cosh(z), np.sinh(z)], dtype=float),
        "theta": float(k)
        * np.array([-np.sin(z), np.cos(z), np.sinh(z), np.cosh(z)], dtype=float),
        "curvature": float(k) ** 2
        * np.array([-np.cos(z), -np.sin(z), np.cosh(z), np.sinh(z)], dtype=float),
        "third": float(k) ** 3
        * np.array([np.sin(z), -np.cos(z), np.sinh(z), np.cosh(z)], dtype=float),
    }


def analytic_bending_dynamic_stiffness(
    *,
    length: float,
    Lambda: float,
    joint_at_right: bool,
) -> np.ndarray:
    k = float(Lambda)
    rows_0 = full_basis_rows(k, 0.0)
    rows_L = full_basis_rows(k, float(length))
    if joint_at_right:
        system = np.vstack([rows_0["w"], rows_0["theta"], rows_L["w"], rows_L["theta"]])
        bc_template = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)
        joint_positions = (2, 3)
        force_rows = np.vstack([-rows_L["third"], rows_L["curvature"]])
    else:
        system = np.vstack([rows_0["w"], rows_0["theta"], rows_L["w"], rows_L["theta"]])
        bc_template = np.array([0.0, 0.0, 0.0, 0.0], dtype=float)
        joint_positions = (0, 1)
        force_rows = np.vstack([rows_0["third"], -rows_0["curvature"]])

    D = np.zeros((2, 2), dtype=float)
    for col, pos in enumerate(joint_positions):
        rhs = bc_template.copy()
        rhs[pos] = 1.0
        coeff = np.linalg.solve(system, rhs)
        D[:, col] = force_rows @ coeff
    return D


def axial_value_rows(gamma: float, x: float) -> tuple[np.ndarray, np.ndarray]:
    z = float(gamma) * float(x)
    u = np.array([np.cos(z), np.sin(z)], dtype=float)
    du = float(gamma) * np.array([-np.sin(z), np.cos(z)], dtype=float)
    return u, du


def analytic_axial_dynamic_stiffness(
    *,
    length: float,
    Lambda: float,
    epsilon: float,
    joint_at_right: bool,
) -> float:
    omega_sq = float(Lambda) ** 4
    EA = 1.0 / float(epsilon) ** 2
    gamma = np.sqrt(max(omega_sq / EA, 0.0))
    u0, du0 = axial_value_rows(gamma, 0.0)
    uL, duL = axial_value_rows(gamma, float(length))
    if joint_at_right:
        system = np.vstack([u0, uL])
        rhs = np.array([0.0, 1.0], dtype=float)
        force_row = EA * duL
    else:
        system = np.vstack([u0, uL])
        rhs = np.array([1.0, 0.0], dtype=float)
        force_row = -EA * du0
    coeff = np.linalg.solve(system, rhs)
    return float(force_row @ coeff)


def analytic_arm_dynamic_stiffness_local(
    *,
    length: float,
    Lambda: float,
    epsilon: float,
    joint_at_right: bool,
) -> np.ndarray:
    D = np.zeros((3, 3), dtype=float)
    D[0, 0] = analytic_axial_dynamic_stiffness(
        length=float(length),
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        joint_at_right=joint_at_right,
    )
    Db = analytic_bending_dynamic_stiffness(
        length=float(length),
        Lambda=float(Lambda),
        joint_at_right=joint_at_right,
    )
    D[np.ix_([1, 2], [1, 2])] = Db
    return D


def matrix_entry_rows(
    *,
    branch_id: str,
    beta: float,
    mu: float,
    epsilon: float,
    Lambda: float,
    omega_sq: float,
    matrices: dict[str, tuple[np.ndarray, np.ndarray]],
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for matrix_name, (fem_matrix, analytic_matrix) in matrices.items():
        diff = np.asarray(fem_matrix, dtype=float) - np.asarray(analytic_matrix, dtype=float)
        for row_idx, row_label in enumerate(DOF_LABELS):
            for col_idx, col_label in enumerate(JOINT_Q_LABELS):
                fem_value = float(fem_matrix[row_idx, col_idx])
                analytic_value = float(analytic_matrix[row_idx, col_idx])
                rows.append(
                    {
                        "branch_id": branch_id,
                        "beta": float(beta),
                        "mu": float(mu),
                        "epsilon": float(epsilon),
                        "analytic_lambda": float(Lambda),
                        "omega_sq": float(omega_sq),
                        "matrix_name": matrix_name,
                        "row": row_label,
                        "col": col_label,
                        "fem_value": fem_value,
                        "analytic_value": analytic_value,
                        "abs_diff": float(abs(diff[row_idx, col_idx])),
                        "rel_diff": safe_relative_diff(fem_value, analytic_value),
                    }
                )
    return rows


def scaled_force_rows_for_determinant(force_global: np.ndarray, *, Lambda: float, epsilon: float) -> np.ndarray:
    force = np.asarray(force_global, dtype=float)
    return np.array(
        [
            force[2] / (float(Lambda) ** 2),
            float(epsilon) * force[1] / (float(Lambda) ** 2),
            float(epsilon) * force[0] / (float(Lambda) ** 2),
        ],
        dtype=float,
    )


def determinant_crosscheck_rows(
    *,
    branch_id: str,
    beta: float,
    mu: float,
    epsilon: float,
    Lambda: float,
    matrix_residual: np.ndarray,
    dynamic_force_analytic: np.ndarray,
    dynamic_force_fem: np.ndarray,
) -> list[dict[str, float | str]]:
    determinant_values = np.asarray(matrix_residual, dtype=float)[3:6]
    dynamic_scaled = scaled_force_rows_for_determinant(dynamic_force_analytic, Lambda=Lambda, epsilon=epsilon)
    dynamic_fem_scaled = scaled_force_rows_for_determinant(dynamic_force_fem, Lambda=Lambda, epsilon=epsilon)
    rows: list[dict[str, float | str]] = []
    for offset, label in enumerate(ROW_LABELS[3:6]):
        det_value = float(determinant_values[offset])
        dyn_value = float(dynamic_scaled[offset])
        rows.append(
            {
                "branch_id": branch_id,
                "beta": float(beta),
                "mu": float(mu),
                "epsilon": float(epsilon),
                "analytic_lambda": float(Lambda),
                "row_index": int(offset + 3),
                "row_label": str(label),
                "determinant_matrix_value": det_value,
                "dynamic_total_scaled_value": dyn_value,
                "dynamic_total_fem_scaled_value": float(dynamic_fem_scaled[offset]),
                "abs_diff": float(abs(det_value - dyn_value)),
                "rel_diff": safe_relative_diff(dyn_value, det_value),
                "scale_note": "assumed determinant order [M/Lambda^2, eps*Fy/Lambda^2, eps*Fx/Lambda^2]",
            }
        )
    return rows


def reduced_bending_values(a_value: float, b_value: float, z: float) -> dict[str, float]:
    a = float(a_value)
    b = float(b_value)
    zz = float(z)
    return {
        "w": float(a * (np.cos(zz) - np.cosh(zz)) + b * (np.sin(zz) - np.sinh(zz))),
        "slope_z": float(a * (-np.sin(zz) - np.sinh(zz)) + b * (np.cos(zz) - np.cosh(zz))),
        "curvature_z": float(a * (-np.cos(zz) - np.cosh(zz)) + b * (-np.sin(zz) - np.sinh(zz))),
        "third_z": float(a * (np.sin(zz) - np.sinh(zz)) + b * (-np.cos(zz) - np.cosh(zz))),
    }


def endpoint_residual_from_coeff(
    coeff: np.ndarray,
    *,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
) -> np.ndarray:
    A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
    x1 = float(Lambda) * (1.0 - float(mu))
    x2 = float(Lambda) * (1.0 + float(mu))
    z2 = -x2
    th1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu))
    th2 = -float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu))
    cb, sb = float(np.cos(beta_rad)), float(np.sin(beta_rad))
    left = reduced_bending_values(A1, B1, x1)
    right = reduced_bending_values(A2, B2, z2)
    u1 = P1 * np.sin(th1)
    u2 = P2 * np.sin(th2)
    axial1 = P1 * np.cos(th1)
    axial2 = P2 * np.cos(th2)
    shear1 = -float(epsilon) * float(Lambda) * left["third_z"]
    shear2 = -float(epsilon) * float(Lambda) * right["third_z"]
    return np.array(
        [
            left["w"] - right["w"] * cb - u2 * sb,
            u1 + right["w"] * sb - u2 * cb,
            left["slope_z"] - right["slope_z"],
            left["curvature_z"] - right["curvature_z"],
            shear1 - shear2 * cb - axial2 * sb,
            axial1 + shear2 * sb - axial2 * cb,
        ],
        dtype=float,
    )


def independent_endpoint_matrix(
    *,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
) -> np.ndarray:
    rows = []
    for col in range(6):
        basis = np.zeros(6, dtype=float)
        basis[col] = 1.0
        rows.append(
            endpoint_residual_from_coeff(
                basis,
                Lambda=Lambda,
                beta_rad=beta_rad,
                mu=mu,
                epsilon=epsilon,
            )
        )
    return np.column_stack(rows)


def joint_q_matrix(
    *,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
    q_source_variant: str,
) -> np.ndarray:
    left_columns = []
    right_columns = []
    cb, sb = float(np.cos(beta_rad)), float(np.sin(beta_rad))
    for col in range(6):
        coeff = np.zeros(6, dtype=float)
        coeff[col] = 1.0
        A1, B1, A2, B2, P1, P2 = [float(value) for value in coeff]
        x1 = float(Lambda) * (1.0 - float(mu))
        x2 = float(Lambda) * (1.0 + float(mu))
        z2 = -x2
        th1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu))
        th2 = -float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu))
        left = reduced_bending_values(A1, B1, x1)
        right = reduced_bending_values(A2, B2, z2)
        theta_scale = float(Lambda) if "fem_theta" in q_source_variant else 1.0
        q_left = np.array(
            [
                P1 * np.sin(th1),
                left["w"],
                theta_scale * left["slope_z"],
            ],
            dtype=float,
        )
        u2 = P2 * np.sin(th2)
        q_right = np.array(
            [
                cb * u2 - sb * right["w"],
                sb * u2 + cb * right["w"],
                theta_scale * right["slope_z"],
            ],
            dtype=float,
        )
        left_columns.append(q_left)
        right_columns.append(q_right)
    q_left_matrix = np.column_stack(left_columns)
    q_right_matrix = np.column_stack(right_columns)
    if q_source_variant.startswith("left_joint"):
        return q_left_matrix
    if q_source_variant.startswith("right_joint_global"):
        return q_right_matrix
    if q_source_variant.startswith("average_joint"):
        return 0.5 * (q_left_matrix + q_right_matrix)
    raise ValueError(f"Unknown q_source_variant: {q_source_variant}")


def q_source_variants() -> tuple[str, ...]:
    return (
        "left_joint_fem_theta",
        "right_joint_global_fem_theta",
        "average_joint_fem_theta",
        "left_joint_slope_theta",
        "right_joint_global_slope_theta",
        "average_joint_slope_theta",
    )


def force_order_variants() -> dict[str, tuple[tuple[int, float | str], ...]]:
    return {
        "M_epsFy_epsFx": ((2, "moment"), (1, "eps_force"), (0, "eps_force")),
        "M_epsFx_epsFy": ((2, "moment"), (0, "eps_force"), (1, "eps_force")),
        "M_Fy_Fx_unscaled_force": ((2, "moment"), (1, "force"), (0, "force")),
        "M_Fx_Fy_unscaled_force": ((2, "moment"), (0, "force"), (1, "force")),
    }


def sign_variants() -> tuple[tuple[str, np.ndarray], ...]:
    variants = []
    for m in (-1.0, 1.0):
        for t in (-1.0, 1.0):
            for a in (-1.0, 1.0):
                label = f"{'+' if m > 0 else '-'}{'+' if t > 0 else '-'}{'+' if a > 0 else '-'}"
                variants.append((label, np.array([m, t, a], dtype=float)))
    return tuple(variants)


def force_scale(scale_kind: float | str, *, Lambda: float, epsilon: float) -> float:
    if scale_kind == "moment":
        return 1.0 / (float(Lambda) ** 2)
    if scale_kind == "eps_force":
        return float(epsilon) / (float(Lambda) ** 2)
    if scale_kind == "force":
        return 1.0 / (float(Lambda) ** 2)
    return float(scale_kind)


def dynamic_stiffness_matrix_variant(
    *,
    endpoint_matrix: np.ndarray,
    D_total: np.ndarray,
    Lambda: float,
    epsilon: float,
    q_matrix: np.ndarray,
    force_order_variant: str,
    sign_values: np.ndarray,
) -> np.ndarray:
    force_matrix = np.asarray(D_total, dtype=float) @ np.asarray(q_matrix, dtype=float)
    out = np.asarray(endpoint_matrix, dtype=float).copy()
    order = force_order_variants()[force_order_variant]
    for row_offset, (force_index, scale_kind) in enumerate(order):
        out[3 + row_offset, :] = (
            float(sign_values[row_offset])
            * force_scale(scale_kind, Lambda=Lambda, epsilon=epsilon)
            * force_matrix[int(force_index), :]
        )
    return out


def row_audit_matrix_rows(
    *,
    case_id: str,
    current: np.ndarray,
    comparison: np.ndarray,
    matrix_pair: str,
    force_order_variant: str,
    sign_variant: str,
    notes: str,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for row_idx, row_label in enumerate(ROW_LABELS):
        for col_idx, col_label in enumerate(BASIS_LABELS):
            current_value = float(current[row_idx, col_idx])
            comparison_value = float(comparison[row_idx, col_idx])
            rows.append(
                {
                    "case_id": case_id,
                    "matrix_pair": matrix_pair,
                    "row_index": int(row_idx),
                    "row_label": str(row_label),
                    "col_index": int(col_idx),
                    "col_label": str(col_label),
                    "current_value": current_value,
                    "comparison_value": comparison_value,
                    "abs_diff": float(abs(current_value - comparison_value)),
                    "rel_diff": safe_relative_diff(current_value, comparison_value),
                    "force_order_variant": force_order_variant,
                    "sign_variant": sign_variant,
                    "notes": notes,
                }
            )
    return rows


def best_row_scalar(current_row: np.ndarray, comparison_row: np.ndarray) -> tuple[float, float, float, float]:
    current = np.asarray(current_row, dtype=float)
    comparison = np.asarray(comparison_row, dtype=float)
    denom = float(current @ current)
    best_scalar = float((current @ comparison) / denom) if denom > NEAR_ZERO_NORM else np.nan
    reference_norm = float(np.linalg.norm(comparison))
    before = (
        float(np.linalg.norm(current - comparison) / reference_norm)
        if reference_norm > NEAR_ZERO_NORM
        else np.nan
    )
    after = (
        float(np.linalg.norm(best_scalar * current - comparison) / reference_norm)
        if reference_norm > NEAR_ZERO_NORM and np.isfinite(best_scalar)
        else np.nan
    )
    corr_denom = float(np.linalg.norm(current) * np.linalg.norm(comparison))
    correlation = float((current @ comparison) / corr_denom) if corr_denom > NEAR_ZERO_NORM else np.nan
    return best_scalar, before, after, correlation


def row_scaling_audit_rows(
    *,
    case_id: str,
    current: np.ndarray,
    comparison: np.ndarray,
    matrix_pair: str,
    force_order_variant: str,
    sign_variant: str,
    notes: str,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for row_idx, row_label in enumerate(ROW_LABELS):
        best_scalar, before, after, correlation = best_row_scalar(current[row_idx], comparison[row_idx])
        rows.append(
            {
                "case_id": case_id,
                "matrix_pair": matrix_pair,
                "row_index": int(row_idx),
                "row_label": str(row_label),
                "force_order_variant": force_order_variant,
                "sign_variant": sign_variant,
                "best_scalar": float(best_scalar),
                "row_l2_relative_diff_before_scaling": float(before),
                "row_l2_relative_diff_after_best_scalar": float(after),
                "row_correlation": float(correlation),
                "notes": notes,
            }
        )
    return rows


def force_contribution_matrices_current(
    *,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
) -> dict[tuple[str, str], np.ndarray]:
    cb, sb = float(np.cos(beta_rad)), float(np.sin(beta_rad))
    x1 = float(Lambda) * (1.0 - float(mu))
    x2 = float(Lambda) * (1.0 + float(mu))
    th1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu))
    th2 = float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu))
    Sd1 = np.sin(x1) - np.sinh(x1)
    Cs1 = np.cos(x1) + np.cosh(x1)
    Sd2 = np.sin(x2) - np.sinh(x2)
    Cs2 = np.cos(x2) + np.cosh(x2)
    Ss1 = np.sin(x1) + np.sinh(x1)
    Cs_curv1 = np.cos(x1) + np.cosh(x1)
    Ss2 = np.sin(x2) + np.sinh(x2)
    Cs_curv2 = np.cos(x2) + np.cosh(x2)
    out: dict[tuple[str, str], np.ndarray] = {}

    def vec(values: dict[int, float]) -> np.ndarray:
        row = np.zeros(6, dtype=float)
        for idx, value in values.items():
            row[int(idx)] = float(value)
        return row

    out[("moment_equilibrium", "left_moment_scaled")] = vec({0: -Cs_curv1, 1: -Ss1})
    out[("moment_equilibrium", "right_moment_signed_scaled")] = vec({2: Cs_curv2, 3: -Ss2})
    out[("transverse_global_force_equilibrium", "left_shear_scaled")] = vec(
        {0: -float(epsilon) * float(Lambda) * Sd1, 1: float(epsilon) * float(Lambda) * Cs1}
    )
    out[("transverse_global_force_equilibrium", "right_shear_cos_signed_scaled")] = vec(
        {2: -float(epsilon) * float(Lambda) * Sd2 * cb, 3: -float(epsilon) * float(Lambda) * Cs2 * cb}
    )
    out[("transverse_global_force_equilibrium", "right_axial_sin_signed_scaled")] = vec(
        {5: -np.cos(th2) * sb}
    )
    out[("axial_global_force_equilibrium", "left_axial_scaled")] = vec({4: np.cos(th1)})
    out[("axial_global_force_equilibrium", "right_shear_sin_scaled")] = vec(
        {2: float(epsilon) * float(Lambda) * Sd2 * sb, 3: float(epsilon) * float(Lambda) * Cs2 * sb}
    )
    out[("axial_global_force_equilibrium", "right_axial_cos_signed_scaled")] = vec(
        {5: -np.cos(th2) * cb}
    )
    return out


def force_contribution_matrices_dynamic(
    *,
    D_left: np.ndarray,
    D_right_local: np.ndarray,
    q_matrix: np.ndarray,
    beta_deg: float,
    Lambda: float,
    epsilon: float,
) -> dict[tuple[str, str], np.ndarray]:
    T = rotation3(beta_deg)
    q_global = np.asarray(q_matrix, dtype=float)
    q_right_local = T.T @ q_global
    f_left = np.asarray(D_left, dtype=float) @ q_global
    f_right_local = np.asarray(D_right_local, dtype=float) @ q_right_local
    right_axial_global = T @ np.vstack(
        [f_right_local[0, :], np.zeros(q_global.shape[1]), np.zeros(q_global.shape[1])]
    )
    right_shear_global = T @ np.vstack(
        [np.zeros(q_global.shape[1]), f_right_local[1, :], np.zeros(q_global.shape[1])]
    )
    right_moment_global = T @ np.vstack(
        [np.zeros(q_global.shape[1]), np.zeros(q_global.shape[1]), f_right_local[2, :]]
    )
    force_scale_eps = float(epsilon) / (float(Lambda) ** 2)
    moment_scale = 1.0 / (float(Lambda) ** 2)
    return {
        ("moment_equilibrium", "left_moment_scaled"): moment_scale * f_left[2, :],
        ("moment_equilibrium", "right_moment_signed_scaled"): moment_scale * right_moment_global[2, :],
        ("transverse_global_force_equilibrium", "left_shear_scaled"): force_scale_eps * f_left[1, :],
        (
            "transverse_global_force_equilibrium",
            "right_shear_cos_signed_scaled",
        ): force_scale_eps * right_shear_global[1, :],
        (
            "transverse_global_force_equilibrium",
            "right_axial_sin_signed_scaled",
        ): force_scale_eps * right_axial_global[1, :],
        ("axial_global_force_equilibrium", "left_axial_scaled"): force_scale_eps * f_left[0, :],
        ("axial_global_force_equilibrium", "right_shear_sin_scaled"): force_scale_eps * right_shear_global[0, :],
        (
            "axial_global_force_equilibrium",
            "right_axial_cos_signed_scaled",
        ): force_scale_eps * right_axial_global[0, :],
    }


def force_contribution_matrices_fem_physical(
    *,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
) -> dict[tuple[str, str], np.ndarray]:
    """FEM-compatible physical endpoint force contributions.

    The values are unscaled generalized forces in the current nondimensional
    FEM convention: EI=1, EA=1/epsilon^2, and the right-arm local joint is the
    first node of an element whose positive local x direction points from joint
    to the external clamp.
    """
    cb, sb = float(np.cos(beta_rad)), float(np.sin(beta_rad))
    x1 = float(Lambda) * (1.0 - float(mu))
    x2 = float(Lambda) * (1.0 + float(mu))
    th1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu))
    th2 = float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu))
    Sd1 = np.sin(x1) - np.sinh(x1)
    Cs1 = np.cos(x1) + np.cosh(x1)
    Ss1 = np.sin(x1) + np.sinh(x1)
    Sd2 = np.sin(x2) - np.sinh(x2)
    Cs2 = np.cos(x2) + np.cosh(x2)
    Ss2 = np.sin(x2) + np.sinh(x2)
    lam2 = float(Lambda) ** 2
    lam3 = float(Lambda) ** 3
    axial_scale = lam2 / float(epsilon)

    def vec(values: dict[int, float]) -> np.ndarray:
        row = np.zeros(6, dtype=float)
        for idx, value in values.items():
            row[int(idx)] = float(value)
        return row

    left_moment = vec({0: -lam2 * Cs1, 1: -lam2 * Ss1})
    right_moment_global = vec({2: lam2 * Cs2, 3: -lam2 * Ss2})
    left_shear_global_y = vec({0: -lam3 * Sd1, 1: lam3 * Cs1})
    right_shear_local = vec({2: -lam3 * Sd2, 3: -lam3 * Cs2})
    right_axial_local = vec({5: -axial_scale * np.cos(th2)})
    left_axial_global_x = vec({4: axial_scale * np.cos(th1)})

    return {
        ("moment_equilibrium", "left_moment_scaled"): left_moment,
        ("moment_equilibrium", "right_moment_signed_scaled"): right_moment_global,
        ("transverse_global_force_equilibrium", "left_shear_scaled"): left_shear_global_y,
        (
            "transverse_global_force_equilibrium",
            "right_shear_cos_signed_scaled",
        ): cb * right_shear_local,
        (
            "transverse_global_force_equilibrium",
            "right_axial_sin_signed_scaled",
        ): -sb * right_axial_local,
        ("axial_global_force_equilibrium", "left_axial_scaled"): left_axial_global_x,
        ("axial_global_force_equilibrium", "right_shear_sin_scaled"): sb * right_shear_local,
        (
            "axial_global_force_equilibrium",
            "right_axial_cos_signed_scaled",
        ): cb * right_axial_local,
    }


def candidate_factor_values(
    *,
    Lambda: float,
    epsilon: float,
    mu: float,
) -> dict[str, float]:
    lam = float(Lambda)
    eps = float(epsilon)
    l1 = 1.0 - float(mu)
    l2 = 1.0 + float(mu)
    raw = {
        "1": 1.0,
        "epsilon": eps,
        "1/epsilon": 1.0 / eps,
        "Lambda": lam,
        "Lambda^2": lam**2,
        "1/Lambda": 1.0 / lam,
        "1/Lambda^2": 1.0 / (lam**2),
        "epsilon/Lambda": eps / lam,
        "epsilon/Lambda^2": eps / (lam**2),
        "epsilon*Lambda": eps * lam,
        "epsilon*Lambda^2": eps * (lam**2),
        "1/(epsilon*Lambda)": 1.0 / (eps * lam),
        "1/(epsilon*Lambda^2)": 1.0 / (eps * lam**2),
        "L1": l1,
        "L2": l2,
        "1/L1": 1.0 / l1 if abs(l1) > NEAR_ZERO_NORM else np.nan,
        "1/L2": 1.0 / l2 if abs(l2) > NEAR_ZERO_NORM else np.nan,
        "L1^2": l1**2,
        "L2^2": l2**2,
        "1/L1^2": 1.0 / (l1**2) if abs(l1) > NEAR_ZERO_NORM else np.nan,
        "1/L2^2": 1.0 / (l2**2) if abs(l2) > NEAR_ZERO_NORM else np.nan,
        "L1^3": l1**3,
        "L2^3": l2**3,
        "1/L1^3": 1.0 / (l1**3) if abs(l1) > NEAR_ZERO_NORM else np.nan,
        "1/L2^3": 1.0 / (l2**3) if abs(l2) > NEAR_ZERO_NORM else np.nan,
    }
    out: dict[str, float] = {}
    for label, value in raw.items():
        if np.isfinite(float(value)):
            out[label] = float(value)
            out[f"-{label}"] = -float(value)
    return out


def candidate_factor_match(
    ratio: float,
    *,
    Lambda: float,
    epsilon: float,
    mu: float,
    tolerance: float = 1e-7,
) -> str:
    if not np.isfinite(float(ratio)):
        return ""
    candidates = candidate_factor_values(Lambda=Lambda, epsilon=epsilon, mu=mu)
    best_label = ""
    best_error = np.inf
    for label, value in candidates.items():
        denom = max(abs(float(value)), NEAR_ZERO_NORM)
        error = abs(float(ratio) - float(value)) / denom
        if error < best_error:
            best_error = error
            best_label = label
    return best_label if best_error <= tolerance else f"unmatched(best={best_label}, relerr={best_error:.3e})"


def force_row_scaling_audit_rows(
    *,
    case_id: str,
    current_contrib: dict[tuple[str, str], np.ndarray],
    fem_contrib: dict[tuple[str, str], np.ndarray],
    Lambda: float,
    epsilon: float,
    mu: float,
    notes: str,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for key in current_contrib:
        row_label, contribution_label = key
        det = np.asarray(current_contrib[key], dtype=float)
        fem_values = np.asarray(fem_contrib.get(key, np.zeros_like(det)), dtype=float)
        for col_idx, col_label in enumerate(BASIS_LABELS):
            det_value = float(det[col_idx])
            fem_value = float(fem_values[col_idx])
            ratio = det_value / fem_value if abs(fem_value) > NEAR_ZERO_NORM else np.nan
            rows.append(
                {
                    "case_id": case_id,
                    "row_label": str(row_label),
                    "contribution_label": str(contribution_label),
                    "col_label": str(col_label),
                    "determinant_current_value": det_value,
                    "fem_compatible_value": fem_value,
                    "ratio_det_to_fem": float(ratio),
                    "candidate_factor_match": candidate_factor_match(
                        ratio,
                        Lambda=Lambda,
                        epsilon=epsilon,
                        mu=mu,
                    ),
                    "abs_diff": float(abs(det_value - fem_value)),
                    "rel_diff": safe_relative_diff(det_value, fem_value),
                    "notes": notes,
                }
            )
    return rows


def summarize_force_row_scaling(
    rows: Sequence[dict[str, float | str]],
    *,
    case_id: str,
) -> list[dict[str, float | str]]:
    grouped: dict[tuple[str, str], list[dict[str, float | str]]] = {}
    for row in rows:
        grouped.setdefault((str(row["row_label"]), str(row["contribution_label"])), []).append(row)
    summary: list[dict[str, float | str]] = []
    for (row_label, contribution_label), group in grouped.items():
        finite = [
            row
            for row in group
            if np.isfinite(float(row["ratio_det_to_fem"]))
            and abs(float(row["determinant_current_value"])) > NEAR_ZERO_NORM
        ]
        ratios = np.array([float(row["ratio_det_to_fem"]) for row in finite], dtype=float)
        matches = sorted({str(row["candidate_factor_match"]) for row in finite if str(row["candidate_factor_match"])})
        ratio_pattern = "no_nonzero_overlap"
        likely = ""
        supports = "no"
        if ratios.size:
            spread = float(np.max(ratios) - np.min(ratios))
            scale = max(float(np.max(np.abs(ratios))), NEAR_ZERO_NORM)
            ratio_pattern = "constant" if spread / scale < 1e-7 else "varies_by_column"
            likely = ";".join(matches)
            supports = "yes" if ratio_pattern == "constant" and likely and not likely.startswith("unmatched") else "needs_review"
        summary.append(
            {
                "case_id": case_id,
                "row_label": row_label,
                "contribution_label": contribution_label,
                "ratio_pattern": ratio_pattern,
                "likely_missing_factor": likely,
                "supports_scaling_mismatch": supports,
                "nonzero_ratio_count": int(ratios.size),
                "ratio_min": float(np.min(ratios)) if ratios.size else np.nan,
                "ratio_max": float(np.max(ratios)) if ratios.size else np.nan,
                "notes": (
                    "constant det/fem ratio identifies the multiplier that maps FEM-compatible "
                    "physical contribution into the current determinant contribution"
                ),
            }
        )
    return summary


def diagnostic_matrix_from_fem_physical_contributions(
    *,
    matrix: np.ndarray,
    fem_contrib: dict[tuple[str, str], np.ndarray],
    Lambda: float,
    epsilon: float,
) -> np.ndarray:
    diagnostic = np.asarray(matrix, dtype=float).copy()
    moment_scale = 1.0 / (float(Lambda) ** 2)
    force_scale_eps = float(epsilon) / (float(Lambda) ** 2)
    for row_offset, row_label in enumerate(ROW_LABELS[3:6]):
        scale = moment_scale if row_label == "moment_equilibrium" else force_scale_eps
        total = np.zeros(6, dtype=float)
        for (key_row, _), contribution in fem_contrib.items():
            if key_row == row_label:
                total += np.asarray(contribution, dtype=float)
        diagnostic[3 + row_offset, :] = scale * total
    return diagnostic


def diagnostic_force_matrix_rows(
    *,
    case_id: str,
    matrix: np.ndarray,
    diagnostic_matrix: np.ndarray,
    coeff: np.ndarray,
    Lambda: float,
    mu: float,
    epsilon: float,
    beta: float,
    params,
) -> list[dict[str, float | str]]:
    current_force = np.asarray(matrix[3:6, :], dtype=float) @ np.asarray(coeff, dtype=float)
    diagnostic_force = np.asarray(diagnostic_matrix[3:6, :], dtype=float) @ np.asarray(coeff, dtype=float)
    diagnostic_coeff, diagnostic_smin, _ = analytic_null_vector(diagnostic_matrix)
    diagnostic_matrix_residual = np.asarray(diagnostic_matrix, dtype=float) @ diagnostic_coeff
    fields = analytic_local_fields(
        Lambda,
        mu=mu,
        epsilon=epsilon,
        coeff=diagnostic_coeff,
        s_norm=np.linspace(0.0, 1.0, fem.N_ELEM + 1),
        right_theta_sign=1.0,
    )
    q_full, _ = build_analytic_global_q(fields, beta_deg=beta)
    k_free, m_free, _, _, free = assemble_fem_matrices(params, mu=mu, beta_deg=beta)
    residual = residual_metrics(k_free, m_free, q_full, free, omega_sq=float(Lambda) ** 4)
    return [
        {
            "case_id": case_id,
            "diagnostic_variant": "fem_physical_common_row_scaling",
            "current_null_current_force_norm": float(np.linalg.norm(current_force)),
            "diagnostic_force_on_current_null_norm": float(np.linalg.norm(diagnostic_force)),
            "diagnostic_new_null_residual_norm": float(np.linalg.norm(diagnostic_matrix_residual)),
            "diagnostic_new_null_smallest_singular_value": float(diagnostic_smin),
            "diagnostic_new_null_analytic_in_fem_relative_residual": float(residual["relative_residual"]),
            "diagnostic_new_null_analytic_in_fem_residual_l2": float(residual["residual_l2"]),
            "notes": (
                "diagnostic-only matrix keeps kinematic rows and replaces force rows by "
                "[M/Lambda^2, epsilon*Fy/Lambda^2, epsilon*Fx/Lambda^2] from FEM-compatible endpoint forces"
            ),
        }
    ]


def contribution_audit_rows(
    *,
    case_id: str,
    current_contrib: dict[tuple[str, str], np.ndarray],
    independent_contrib: dict[tuple[str, str], np.ndarray],
    dynamic_contrib: dict[tuple[str, str], np.ndarray],
    notes: str,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for key in current_contrib:
        row_label, contribution_label = key
        cur = np.asarray(current_contrib[key], dtype=float)
        endp = np.asarray(independent_contrib.get(key, np.zeros_like(cur)), dtype=float)
        dyn = np.asarray(dynamic_contrib.get(key, np.zeros_like(cur)), dtype=float)
        for col_idx, col_label in enumerate(BASIS_LABELS):
            current_value = float(cur[col_idx])
            endpoint_value = float(endp[col_idx])
            dynamic_value = float(dyn[col_idx])
            inferred = (
                float(dynamic_value / current_value)
                if abs(current_value) > NEAR_ZERO_NORM
                else np.nan
            )
            rows.append(
                {
                    "case_id": case_id,
                    "row_label": row_label,
                    "contribution_label": contribution_label,
                    "col_label": col_label,
                    "current_value": current_value,
                    "independent_endpoint_value": endpoint_value,
                    "dynamic_stiffness_value": dynamic_value,
                    "abs_diff_current_vs_endpoint": float(abs(current_value - endpoint_value)),
                    "abs_diff_current_vs_dynamic": float(abs(current_value - dynamic_value)),
                    "rel_diff_current_vs_endpoint": safe_relative_diff(current_value, endpoint_value),
                    "rel_diff_current_vs_dynamic": safe_relative_diff(current_value, dynamic_value),
                    "inferred_missing_factor": inferred,
                    "notes": notes,
                }
            )
    return rows


def null_vector_audit_rows(
    *,
    case_id: str,
    current: np.ndarray,
    endpoint: np.ndarray,
    dynamic: np.ndarray,
    coeff: np.ndarray,
    notes: str,
) -> list[dict[str, float | str]]:
    cur = np.asarray(current, dtype=float) @ np.asarray(coeff, dtype=float)
    endp = np.asarray(endpoint, dtype=float) @ np.asarray(coeff, dtype=float)
    dyn = np.asarray(dynamic, dtype=float) @ np.asarray(coeff, dtype=float)
    rows: list[dict[str, float | str]] = []
    for row_idx, row_label in enumerate(ROW_LABELS):
        rows.append(
            {
                "case_id": case_id,
                "row_label": str(row_label),
                "current_residual": float(cur[row_idx]),
                "independent_endpoint_residual": float(endp[row_idx]),
                "dynamic_stiffness_residual": float(dyn[row_idx]),
                "current_vs_endpoint_abs_diff": float(abs(cur[row_idx] - endp[row_idx])),
                "current_vs_dynamic_abs_diff": float(abs(cur[row_idx] - dyn[row_idx])),
                "notes": notes,
            }
        )
    return rows


def best_scaling_rows(
    rows: Sequence[dict[str, float | str]],
    *,
    matrix_pair: str,
    force_only: bool = False,
) -> list[dict[str, float | str]]:
    selected = [row for row in rows if str(row["matrix_pair"]) == matrix_pair]
    if force_only:
        selected = [row for row in selected if int(row["row_index"]) >= 3]
    grouped: dict[int, list[dict[str, float | str]]] = {}
    for row in selected:
        grouped.setdefault(int(row["row_index"]), []).append(row)
    best = []
    for row_idx in sorted(grouped):
        finite = [
            row
            for row in grouped[row_idx]
            if np.isfinite(float(row["row_l2_relative_diff_after_best_scalar"]))
        ]
        if finite:
            best.append(min(finite, key=lambda row: float(row["row_l2_relative_diff_after_best_scalar"])))
    return best


def format_float_sequence(values: Sequence[float]) -> str:
    return ";".join(f"{float(value):.16g}" for value in values)


def kinematic_nullspace_basis(Rk: np.ndarray) -> tuple[np.ndarray, np.ndarray, int, float, float, float]:
    matrix = np.asarray(Rk, dtype=float)
    _, singular_values, vh = np.linalg.svd(matrix, full_matrices=True)
    largest = float(singular_values[0]) if singular_values.size else 0.0
    tol = float(max(matrix.shape) * np.finfo(float).eps * largest)
    rank = int(np.sum(singular_values > tol))
    basis = vh[rank:, :].T.copy()
    condition = (
        float(singular_values[0] / singular_values[-1])
        if singular_values.size and abs(float(singular_values[-1])) > NEAR_ZERO_NORM
        else np.inf
    )
    residual_norm = float(np.linalg.norm(matrix @ basis)) if basis.size else np.nan
    return basis, singular_values, rank, tol, condition, residual_norm


def dynamic_force_row_operator(
    *,
    D_total: np.ndarray,
    q_matrix: np.ndarray,
    Lambda: float,
    epsilon: float,
    force_order_variant: str,
    sign_values: np.ndarray,
) -> np.ndarray:
    force_matrix = np.asarray(D_total, dtype=float) @ np.asarray(q_matrix, dtype=float)
    out = np.zeros((3, force_matrix.shape[1]), dtype=float)
    order = force_order_variants()[force_order_variant]
    for row_offset, (force_index, scale_kind) in enumerate(order):
        out[row_offset, :] = (
            float(sign_values[row_offset])
            * force_scale(scale_kind, Lambda=Lambda, epsilon=epsilon)
            * force_matrix[int(force_index), :]
        )
    return out


def row_scalars_for_force_subspace(
    determinant: np.ndarray,
    dynamic: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    determinant = np.asarray(determinant, dtype=float)
    dynamic = np.asarray(dynamic, dtype=float)
    scalars = np.zeros(determinant.shape[0], dtype=float)
    errors = np.zeros(determinant.shape[0], dtype=float)
    correlations = np.zeros(determinant.shape[0], dtype=float)
    for row_idx in range(determinant.shape[0]):
        best_scalar, _, after, correlation = best_row_scalar(determinant[row_idx], dynamic[row_idx])
        scalars[row_idx] = float(best_scalar)
        errors[row_idx] = float(after)
        correlations[row_idx] = float(correlation)
    return scalars, errors, correlations


def force_subspace_summary_flag(
    *,
    subspace_dim: int,
    fro_before: float,
    fro_after: float,
    max_row_after: float,
) -> str:
    if int(subspace_dim) != 3:
        return "unexpected_kinematic_nullspace_dimension"
    if np.isfinite(fro_before) and fro_before < 1e-8:
        return "matches_without_scaling_on_kinematic_subspace"
    if np.isfinite(fro_after) and fro_after < 1e-8:
        return "matches_after_row_scaling_on_kinematic_subspace"
    if np.isfinite(max_row_after) and max_row_after < 1e-2:
        return "near_match_after_row_scaling_on_kinematic_subspace"
    return "mismatch_on_kinematic_subspace"


def constrained_force_audit_for_variant(
    *,
    case_id: str,
    variant_name: str,
    Rf: np.ndarray,
    nullspace_basis: np.ndarray,
    dynamic_operator: np.ndarray,
    singular_values: np.ndarray,
    rank: int,
    tol: float,
    condition: float,
    kinematic_residual_norm: float,
    q_left_right_mismatch: float,
    notes: str,
) -> tuple[list[dict[str, float | str]], dict[str, float | str]]:
    N = np.asarray(nullspace_basis, dtype=float)
    determinant_subspace = np.asarray(Rf, dtype=float) @ N
    dynamic_subspace = np.asarray(dynamic_operator, dtype=float) @ N
    denom = float(np.linalg.norm(dynamic_subspace))
    fro_before = (
        float(np.linalg.norm(determinant_subspace - dynamic_subspace) / denom)
        if denom > NEAR_ZERO_NORM
        else np.nan
    )
    scalars, row_errors, row_correlations = row_scalars_for_force_subspace(
        determinant_subspace,
        dynamic_subspace,
    )
    scaled = determinant_subspace.copy()
    for row_idx in range(scaled.shape[0]):
        scaled[row_idx, :] *= scalars[row_idx]
    fro_after = (
        float(np.linalg.norm(scaled - dynamic_subspace) / denom)
        if denom > NEAR_ZERO_NORM
        else np.nan
    )
    max_row_error = float(np.nanmax(row_errors)) if row_errors.size else np.nan

    detail_rows: list[dict[str, float | str]] = []
    for row_offset, row_label in enumerate(ROW_LABELS[3:6]):
        for subspace_col in range(N.shape[1]):
            determinant_value = float(determinant_subspace[row_offset, subspace_col])
            dynamic_value = float(dynamic_subspace[row_offset, subspace_col])
            detail_rows.append(
                {
                    "case_id": case_id,
                    "variant_name": variant_name,
                    "subspace_dim": int(N.shape[1]),
                    "row_label": str(row_label),
                    "row_index": int(row_offset + 3),
                    "subspace_col": int(subspace_col),
                    "determinant_value": determinant_value,
                    "dynamic_value": dynamic_value,
                    "abs_diff": float(abs(determinant_value - dynamic_value)),
                    "rel_diff": safe_relative_diff(determinant_value, dynamic_value),
                    "row_best_scalar": float(scalars[row_offset]),
                    "row_error_after_scaling": float(row_errors[row_offset]),
                    "notes": notes,
                }
            )

    summary = {
        "case_id": case_id,
        "variant_name": variant_name,
        "subspace_dim": int(N.shape[1]),
        "Rk_singular_values": format_float_sequence(singular_values),
        "Rk_rank": int(rank),
        "Rk_nullspace_tol": float(tol),
        "Rk_condition_number": float(condition),
        "kinematic_nullspace_residual_norm": float(kinematic_residual_norm),
        "q_left_right_subspace_mismatch": float(q_left_right_mismatch),
        "frobenius_rel_diff_before_scaling": float(fro_before),
        "frobenius_rel_diff_after_row_scaling": float(fro_after),
        "max_row_error_after_scaling": max_row_error,
        "best_row_scalars": format_float_sequence(scalars),
        "row_correlations": format_float_sequence(row_correlations),
        "conclusion_flag": force_subspace_summary_flag(
            subspace_dim=N.shape[1],
            fro_before=fro_before,
            fro_after=fro_after,
            max_row_after=max_row_error,
        ),
        "notes": notes,
    }
    return detail_rows, summary


def constrained_null_vector_rows(
    *,
    case_id: str,
    variant_name: str,
    Rk: np.ndarray,
    Rf: np.ndarray,
    nullspace_basis: np.ndarray,
    dynamic_operator: np.ndarray,
    coeff: np.ndarray,
    notes: str,
) -> list[dict[str, float | str]]:
    N = np.asarray(nullspace_basis, dtype=float)
    coeff_vector = np.asarray(coeff, dtype=float)
    if N.size:
        coordinates, *_ = np.linalg.lstsq(N, coeff_vector, rcond=None)
        projected = N @ coordinates
    else:
        projected = np.zeros_like(coeff_vector)
    projection_abs = float(np.linalg.norm(coeff_vector - projected))
    coeff_norm = float(np.linalg.norm(coeff_vector))
    projection_rel = projection_abs / coeff_norm if coeff_norm > NEAR_ZERO_NORM else np.nan
    kinematic_residual = float(np.linalg.norm(np.asarray(Rk, dtype=float) @ coeff_vector))
    det_c0 = np.asarray(Rf, dtype=float) @ coeff_vector
    dyn_c0 = np.asarray(dynamic_operator, dtype=float) @ coeff_vector
    det_projected = np.asarray(Rf, dtype=float) @ projected
    dyn_projected = np.asarray(dynamic_operator, dtype=float) @ projected
    rows: list[dict[str, float | str]] = []
    for row_offset, row_label in enumerate(ROW_LABELS[3:6]):
        rows.append(
            {
                "case_id": case_id,
                "variant_name": variant_name,
                "row_label": str(row_label),
                "row_index": int(row_offset + 3),
                "determinant_c0": float(det_c0[row_offset]),
                "dynamic_c0": float(dyn_c0[row_offset]),
                "determinant_projected": float(det_projected[row_offset]),
                "dynamic_projected": float(dyn_projected[row_offset]),
                "det_vs_dynamic_c0_abs_diff": float(abs(det_c0[row_offset] - dyn_c0[row_offset])),
                "det_vs_dynamic_projected_abs_diff": float(
                    abs(det_projected[row_offset] - dyn_projected[row_offset])
                ),
                "c0_projection_abs_error": projection_abs,
                "c0_projection_relative_error": projection_rel,
                "kinematic_residual_norm": kinematic_residual,
                "notes": notes,
            }
        )
    return rows


def safe_scalar_ratio(numerator: float, denominator: float) -> float:
    return float(numerator) / float(denominator) if abs(float(denominator)) > NEAR_ZERO_NORM else np.nan


def condition_number_from_singular_values(singular_values: np.ndarray) -> float:
    s = np.asarray(singular_values, dtype=float)
    if s.size == 0:
        return np.nan
    return float(s[0] / s[-1]) if abs(float(s[-1])) > NEAR_ZERO_NORM else float("inf")


def normalized_vector(vector: np.ndarray) -> np.ndarray:
    values = np.asarray(vector, dtype=float)
    norm = float(np.linalg.norm(values))
    return values / norm if norm > NEAR_ZERO_NORM else np.full(values.shape, np.nan, dtype=float)


def right_singular_null_vector(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    _, singular_values, vh = np.linalg.svd(np.asarray(matrix, dtype=float), full_matrices=True)
    vector = normalized_vector(vh[-1, :].copy())
    smallest = float(singular_values[-1]) if singular_values.size else np.nan
    return vector, singular_values, smallest


def component_ratio_note(vector_name: str, vector: np.ndarray) -> str:
    values = np.asarray(vector, dtype=float)
    ux, uy, theta = [float(value) for value in values]
    ratios = {
        "ux/theta": safe_scalar_ratio(ux, theta),
        "uy/theta": safe_scalar_ratio(uy, theta),
        "ux/uy": safe_scalar_ratio(ux, uy),
    }
    formatted = ", ".join(f"{name}={value:.6e}" for name, value in ratios.items())
    return f"{vector_name}: {formatted}"


def joint_nullspace_matrix_rows(
    *,
    case_id: str,
    D_total_fem: np.ndarray,
    D_total_analytic: np.ndarray,
) -> list[dict[str, float | str]]:
    diff = np.asarray(D_total_fem, dtype=float) - np.asarray(D_total_analytic, dtype=float)
    rows: list[dict[str, float | str]] = []
    for row_idx, row_label in enumerate(DOF_LABELS):
        for col_idx, col_label in enumerate(JOINT_Q_LABELS):
            fem_value = float(D_total_fem[row_idx, col_idx])
            analytic_value = float(D_total_analytic[row_idx, col_idx])
            rows.append(
                {
                    "case_id": case_id,
                    "matrix_name": "D_total",
                    "row": row_label,
                    "col": col_label,
                    "fem_value": fem_value,
                    "analytic_value": analytic_value,
                    "abs_diff": float(abs(diff[row_idx, col_idx])),
                    "rel_diff": safe_relative_diff(fem_value, analytic_value),
                }
            )
    return rows


def joint_nullspace_summary_row(
    *,
    case_id: str,
    Lambda: float,
    D_total_fem: np.ndarray,
    D_total_analytic: np.ndarray,
    fem_singular_values: np.ndarray,
    analytic_singular_values: np.ndarray,
) -> dict[str, float | str]:
    diff = np.asarray(D_total_fem, dtype=float) - np.asarray(D_total_analytic, dtype=float)
    return {
        "case_id": case_id,
        "lambda": float(Lambda),
        "fem_singular_values": format_float_sequence(fem_singular_values),
        "analytic_singular_values": format_float_sequence(analytic_singular_values),
        "fem_condition_number": condition_number_from_singular_values(fem_singular_values),
        "analytic_condition_number": condition_number_from_singular_values(analytic_singular_values),
        "frobenius_rel_diff": relative_matrix_norm(diff, D_total_analytic),
        "max_abs_diff": float(np.max(np.abs(diff))),
        "smallest_singular_value_fem": float(fem_singular_values[-1]) if fem_singular_values.size else np.nan,
        "smallest_singular_value_analytic": float(analytic_singular_values[-1])
        if analytic_singular_values.size
        else np.nan,
    }


def joint_vector_rows(
    *,
    case_id: str,
    vectors: dict[str, np.ndarray],
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for vector_name, vector in vectors.items():
        values = np.asarray(vector, dtype=float)
        unit = normalized_vector(values)
        rows.append(
            {
                "case_id": case_id,
                "vector_name": vector_name,
                "ux": float(values[0]),
                "uy": float(values[1]),
                "theta": float(values[2]),
                "normalized_ux": float(unit[0]),
                "normalized_uy": float(unit[1]),
                "normalized_theta": float(unit[2]),
            }
        )
    return rows


def aligned_unit_vector_metrics(reference: np.ndarray, candidate: np.ndarray) -> tuple[float, float, float]:
    a = normalized_vector(reference)
    b = normalized_vector(candidate)
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return np.nan, np.nan, np.nan
    dot = float(np.dot(a, b))
    if dot < 0.0:
        b = -b
        dot = -dot
    relative_l2 = float(np.linalg.norm(a - b) / np.linalg.norm(a))
    return float(dot * dot), relative_l2, dot


def joint_pairwise_rows(
    *,
    case_id: str,
    vectors: dict[str, np.ndarray],
) -> list[dict[str, float | str]]:
    names = list(vectors.keys())
    rows: list[dict[str, float | str]] = []
    for i, name_a in enumerate(names):
        for name_b in names[i + 1 :]:
            mac, relative_l2, dot = aligned_unit_vector_metrics(vectors[name_a], vectors[name_b])
            rows.append(
                {
                    "case_id": case_id,
                    "vector_a": name_a,
                    "vector_b": name_b,
                    "mac": mac,
                    "relative_l2_after_sign_alignment": relative_l2,
                    "dot_after_sign_alignment": dot,
                    "notes": (
                        "unit-vector comparison; "
                        f"{component_ratio_note(name_a, vectors[name_a])}; "
                        f"{component_ratio_note(name_b, vectors[name_b])}"
                    ),
                }
            )
    return rows


def solve_bending_coefficients_from_joint_values(
    *,
    length: float,
    Lambda: float,
    joint_at_right: bool,
    w_joint: float,
    theta_joint: float,
) -> np.ndarray:
    rows_0 = full_basis_rows(float(Lambda), 0.0)
    rows_L = full_basis_rows(float(Lambda), float(length))
    system = np.vstack([rows_0["w"], rows_0["theta"], rows_L["w"], rows_L["theta"]])
    rhs = (
        np.array([0.0, 0.0, float(w_joint), float(theta_joint)], dtype=float)
        if joint_at_right
        else np.array([float(w_joint), float(theta_joint), 0.0, 0.0], dtype=float)
    )
    return np.linalg.solve(system, rhs)


def solve_axial_coefficients_from_joint_value(
    *,
    length: float,
    Lambda: float,
    epsilon: float,
    joint_at_right: bool,
    u_joint: float,
) -> np.ndarray:
    omega_sq = float(Lambda) ** 4
    EA = 1.0 / float(epsilon) ** 2
    gamma = np.sqrt(max(omega_sq / EA, 0.0))
    u0, _ = axial_value_rows(gamma, 0.0)
    uL, _ = axial_value_rows(gamma, float(length))
    system = np.vstack([u0, uL])
    rhs = (
        np.array([0.0, float(u_joint)], dtype=float)
        if joint_at_right
        else np.array([float(u_joint), 0.0], dtype=float)
    )
    return np.linalg.solve(system, rhs)


def sample_analytic_arm_from_joint_q(
    *,
    length: float,
    Lambda: float,
    epsilon: float,
    joint_at_right: bool,
    q_local: np.ndarray,
    n_nodes: int,
) -> dict[str, np.ndarray]:
    q = np.asarray(q_local, dtype=float)
    axial_coeff = solve_axial_coefficients_from_joint_value(
        length=float(length),
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        joint_at_right=joint_at_right,
        u_joint=float(q[0]),
    )
    bending_coeff = solve_bending_coefficients_from_joint_values(
        length=float(length),
        Lambda=float(Lambda),
        joint_at_right=joint_at_right,
        w_joint=float(q[1]),
        theta_joint=float(q[2]),
    )
    x_values = np.linspace(0.0, float(length), int(n_nodes))
    omega_sq = float(Lambda) ** 4
    gamma = np.sqrt(max(omega_sq / (1.0 / float(epsilon) ** 2), 0.0))
    u_values = []
    w_values = []
    theta_values = []
    for x_value in x_values:
        u_row, _ = axial_value_rows(gamma, float(x_value))
        rows = full_basis_rows(float(Lambda), float(x_value))
        u_values.append(float(u_row @ axial_coeff))
        w_values.append(float(rows["w"] @ bending_coeff))
        theta_values.append(float(rows["theta"] @ bending_coeff))
    return {
        "u": np.asarray(u_values, dtype=float),
        "w": np.asarray(w_values, dtype=float),
        "theta": np.asarray(theta_values, dtype=float),
    }


def build_q_global_from_joint_null_vector(
    *,
    q_joint_global: np.ndarray,
    Lambda: float,
    mu: float,
    epsilon: float,
    beta: float,
    params,
) -> np.ndarray:
    n = fem.N_ELEM
    n_nodes = 2 * n + 1
    q_joint = np.asarray(q_joint_global, dtype=float)
    length_left = float(params.L_base) * (1.0 - float(mu))
    length_right = float(params.L_base) * (1.0 + float(mu))
    T = rotation3(float(beta))
    left = sample_analytic_arm_from_joint_q(
        length=length_left,
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        joint_at_right=True,
        q_local=q_joint,
        n_nodes=n + 1,
    )
    right_local_joint_q = T.T @ q_joint
    right = sample_analytic_arm_from_joint_q(
        length=length_right,
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        joint_at_right=False,
        q_local=right_local_joint_q,
        n_nodes=n + 1,
    )
    ux = np.zeros(n_nodes, dtype=float)
    uy = np.zeros(n_nodes, dtype=float)
    theta = np.zeros(n_nodes, dtype=float)
    ux[: n + 1] = left["u"]
    uy[: n + 1] = left["w"]
    theta[: n + 1] = left["theta"]
    right_global = T @ np.vstack([right["u"], right["w"], right["theta"]])
    ux[n + 1 :] = right_global[0, 1:]
    uy[n + 1 :] = right_global[1, 1:]
    theta[n + 1 :] = right_global[2, 1:]
    q_full = np.empty(3 * n_nodes, dtype=float)
    q_full[0::3] = ux
    q_full[1::3] = uy
    q_full[2::3] = theta
    return q_full


def joint_reconstruction_residual_rows(
    *,
    case_id: str,
    Lambda: float,
    omega_sq: float,
    mu: float,
    epsilon: float,
    beta: float,
    params,
    q_full_determinant: np.ndarray,
    q_joint_determinant: np.ndarray,
    q_joint_null_analytic: np.ndarray,
) -> list[dict[str, float | str]]:
    k_free, m_free, _, _, free = assemble_fem_matrices(params, mu=float(mu), beta_deg=float(beta))
    q_full_joint_null = build_q_global_from_joint_null_vector(
        q_joint_global=q_joint_null_analytic,
        Lambda=float(Lambda),
        mu=float(mu),
        epsilon=float(epsilon),
        beta=float(beta),
        params=params,
    )
    cases = [
        (
            "determinant_null_vector_reconstruction",
            np.asarray(q_full_determinant, dtype=float),
            np.asarray(q_joint_determinant, dtype=float),
            "current determinant-null coefficient embedding from build_analytic_global_q",
        ),
        (
            "joint_null_analytic_D_reconstruction",
            q_full_joint_null,
            np.asarray(q_joint_null_analytic, dtype=float),
            "arm coefficients solved from q_joint_null_analytic_D using the D_global=T*D_local*T.T convention",
        ),
    ]
    rows: list[dict[str, float | str]] = []
    for name, q_full, q_joint, notes in cases:
        metrics = residual_metrics(k_free, m_free, q_full, free, omega_sq=float(omega_sq))
        rows.append(
            {
                "case_id": case_id,
                "reconstruction_name": name,
                "lambda": float(Lambda),
                "omega_sq": float(omega_sq),
                "residual_l2": float(metrics["residual_l2"]),
                "relative_residual": float(metrics["relative_residual"]),
                "residual_inf": float(metrics["residual_inf"]),
                "relative_residual_inf": float(metrics["relative_residual_inf"]),
                "q_joint_ux": float(q_joint[0]),
                "q_joint_uy": float(q_joint[1]),
                "q_joint_theta": float(q_joint[2]),
                "notes": notes,
            }
        )
    return rows


def pairwise_value(
    rows: Sequence[dict[str, float | str]],
    name_a: str,
    name_b: str,
    key: str,
) -> float:
    wanted = {name_a, name_b}
    for row in rows:
        if {str(row["vector_a"]), str(row["vector_b"])} == wanted:
            return float(row[key])
    return np.nan


def joint_nullspace_interpretation_lines(
    *,
    summary_row: dict[str, float | str],
    pairwise_rows: Sequence[dict[str, float | str]],
) -> list[str]:
    operator_rel = float(summary_row["frobenius_rel_diff"])
    operator_close = np.isfinite(operator_rel) and operator_rel <= 1e-3
    fem_null_mac = pairwise_value(
        pairwise_rows,
        "q_joint_from_FEM_eigenvector",
        "q_joint_null_FEM_D",
        "mac",
    )
    dynamic_null_mac = pairwise_value(
        pairwise_rows,
        "q_joint_null_FEM_D",
        "q_joint_null_analytic_D",
        "mac",
    )
    det_vs_dynamic_mac = pairwise_value(
        pairwise_rows,
        "q_joint_from_analytic_determinant_coeffs",
        "q_joint_null_analytic_D",
        "mac",
    )
    lines = []
    if operator_close and np.isfinite(dynamic_null_mac) and dynamic_null_mac >= 0.999:
        lines.append(
            "A: D_total_FEM and D_total_analytic agree at the joint-displacement operator level."
        )
    elif operator_close:
        lines.append(
            "B-like: D_total operators are close, but the FEM/analytic joint null vectors need review."
        )
    else:
        lines.append(
            "C: D_total_FEM and D_total_analytic are not close enough under the current 1e-3 Frobenius threshold."
        )
    if operator_close and np.isfinite(det_vs_dynamic_mac) and det_vs_dynamic_mac < 0.999:
        lines.append(
            "B: determinant-null q_joint differs from the analytic D_total nullspace; inspect coefficient-to-joint mapping or determinant force rows."
        )
    if np.isfinite(fem_null_mac) and fem_null_mac < 0.999:
        lines.append(
            "D: FEM eigenvector q_joint differs from FEM Schur D_total nullspace; inspect Schur assembly/order or frequency mismatch."
        )
    return lines


def row_unit_scale(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(np.asarray(matrix, dtype=float), axis=1)
    return np.array([1.0 / value if value > NEAR_ZERO_NORM else 1.0 for value in norms], dtype=float)


def column_unit_scale(matrix: np.ndarray, row_scale: np.ndarray | None = None) -> np.ndarray:
    base = np.asarray(matrix, dtype=float)
    if row_scale is not None:
        base = np.asarray(row_scale, dtype=float)[:, None] * base
    norms = np.linalg.norm(base, axis=0)
    return np.array([1.0 / value if value > NEAR_ZERO_NORM else 1.0 for value in norms], dtype=float)


def row_normalized_matrix(matrix: np.ndarray) -> np.ndarray:
    return row_unit_scale(matrix)[:, None] * np.asarray(matrix, dtype=float)


def determinant_coefficients_from_joint_q(
    *,
    q_joint_global: np.ndarray,
    Lambda: float,
    beta: float,
    mu: float,
    epsilon: float,
) -> np.ndarray:
    q_left = np.asarray(q_joint_global, dtype=float)
    q_right = rotation3(float(beta)) @ q_left
    x1 = float(Lambda) * (1.0 - float(mu))
    x2 = float(Lambda) * (1.0 + float(mu))
    z2 = -x2
    left_rows = reduced_bending_values(1.0, 0.0, x1), reduced_bending_values(0.0, 1.0, x1)
    left_system = np.array(
        [
            [left_rows[0]["w"], left_rows[1]["w"]],
            [left_rows[0]["slope_z"], left_rows[1]["slope_z"]],
        ],
        dtype=float,
    )
    right_rows = reduced_bending_values(1.0, 0.0, z2), reduced_bending_values(0.0, 1.0, z2)
    right_system = np.array(
        [
            [right_rows[0]["w"], right_rows[1]["w"]],
            [right_rows[0]["slope_z"], right_rows[1]["slope_z"]],
        ],
        dtype=float,
    )
    A1, B1 = np.linalg.solve(left_system, np.array([q_left[1], q_left[2] / float(Lambda)], dtype=float))
    A2, B2 = np.linalg.solve(right_system, np.array([q_right[1], q_right[2] / float(Lambda)], dtype=float))
    th1 = float(epsilon) * float(Lambda) ** 2 * (1.0 - float(mu))
    th2 = -float(epsilon) * float(Lambda) ** 2 * (1.0 + float(mu))
    P1 = safe_scalar_ratio(float(q_left[0]), float(np.sin(th1)))
    P2 = safe_scalar_ratio(float(q_right[0]), float(np.sin(th2)))
    return np.array([A1, B1, A2, B2, P1, P2], dtype=float)


def normalize_coefficients(coeff: np.ndarray) -> np.ndarray:
    values = np.asarray(coeff, dtype=float)
    norm = float(np.linalg.norm(values))
    out = values / norm if norm > NEAR_ZERO_NORM else values.copy()
    if out.size:
        pivot = int(np.argmax(np.abs(out)))
        if out[pivot] < 0.0:
            out = -out
    return out


def scaled_svd_null_coeff(
    matrix: np.ndarray,
    *,
    row_scale: np.ndarray,
    column_scale: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    scaled = np.asarray(row_scale, dtype=float)[:, None] * np.asarray(matrix, dtype=float)
    scaled = scaled * np.asarray(column_scale, dtype=float)[None, :]
    _, singular_values, vh = np.linalg.svd(scaled, full_matrices=True)
    scaled_coeff = vh[-1, :].astype(float)
    coeff = np.asarray(column_scale, dtype=float) * scaled_coeff
    return normalize_coefficients(coeff), singular_values


def determinant_q_joint_from_coeff(
    *,
    coeff: np.ndarray,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
) -> np.ndarray:
    q_matrix = joint_q_matrix(
        Lambda=float(Lambda),
        beta_rad=float(beta_rad),
        mu=float(mu),
        epsilon=float(epsilon),
        q_source_variant="left_joint_fem_theta",
    )
    return q_matrix @ np.asarray(coeff, dtype=float)


def rayleigh_lambda_for_q(
    *,
    k_free: np.ndarray,
    m_free: np.ndarray,
    q_full: np.ndarray,
    free: np.ndarray,
) -> float:
    q_free = np.asarray(q_full, dtype=float)[np.asarray(free, dtype=int)]
    mass_q = np.asarray(m_free, dtype=float) @ q_free
    denom = float(q_free @ mass_q)
    if abs(denom) <= NEAR_ZERO_NORM:
        return np.nan
    omega_sq_rayleigh = float((q_free @ (np.asarray(k_free, dtype=float) @ q_free)) / denom)
    return float(omega_sq_rayleigh ** 0.25) if omega_sq_rayleigh > 0.0 else np.nan


def determinant_conditioning_method_row(
    *,
    method: str,
    row_scaling_label: str,
    column_scaling_label: str,
    matrix: np.ndarray,
    coeff: np.ndarray,
    singular_values: np.ndarray,
    q_joint_null: np.ndarray,
    Lambda: float,
    beta: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
    params,
    notes: str,
) -> dict[str, float | str]:
    c = normalize_coefficients(coeff)
    raw_residual = np.asarray(matrix, dtype=float) @ c
    row_matrix = row_normalized_matrix(matrix)
    row_residual = row_matrix @ c
    matrix_norm = float(np.linalg.norm(matrix))
    row_matrix_norm = float(np.linalg.norm(row_matrix))
    coeff_norm = float(np.linalg.norm(c))
    raw_rel = (
        float(np.linalg.norm(raw_residual) / (matrix_norm * coeff_norm))
        if matrix_norm > NEAR_ZERO_NORM and coeff_norm > NEAR_ZERO_NORM
        else np.nan
    )
    row_rel = (
        float(np.linalg.norm(row_residual) / (row_matrix_norm * coeff_norm))
        if row_matrix_norm > NEAR_ZERO_NORM and coeff_norm > NEAR_ZERO_NORM
        else np.nan
    )
    q_joint = determinant_q_joint_from_coeff(
        coeff=c,
        Lambda=float(Lambda),
        beta_rad=float(beta_rad),
        mu=float(mu),
        epsilon=float(epsilon),
    )
    mac, q_rel_l2, _ = aligned_unit_vector_metrics(q_joint_null, q_joint)
    fields = analytic_local_fields(
        Lambda,
        mu=mu,
        epsilon=epsilon,
        coeff=c,
        s_norm=np.linspace(0.0, 1.0, fem.N_ELEM + 1),
        right_theta_sign=1.0,
    )
    q_full, _ = build_analytic_global_q(fields, beta_deg=beta)
    k_free, m_free, _, _, free = assemble_fem_matrices(params, mu=mu, beta_deg=beta)
    fem_residual = residual_metrics(k_free, m_free, q_full, free, omega_sq=float(Lambda) ** 4)
    rayleigh_lambda = rayleigh_lambda_for_q(k_free=k_free, m_free=m_free, q_full=q_full, free=free)
    s = np.asarray(singular_values, dtype=float)
    smallest = float(s[-1]) if s.size else np.nan
    second = float(s[-2]) if s.size >= 2 else np.nan
    ratio = safe_scalar_ratio(smallest, second)
    row_components = format_float_sequence(row_residual)
    raw_components = format_float_sequence(raw_residual)
    raw_labeled = ";".join(
        f"{label}={float(value):.16g}" for label, value in zip(ROW_LABELS, raw_residual)
    )
    row_labeled = ";".join(
        f"{label}={float(value):.16g}" for label, value in zip(ROW_LABELS, row_residual)
    )
    return {
        "method": method,
        "row_scaling": row_scaling_label,
        "column_scaling": column_scaling_label,
        "smallest_singular_value": smallest,
        "second_singular_value": second,
        "singular_ratio": ratio,
        "raw_matrix_residual": raw_rel,
        "row_normalized_matrix_residual": row_rel,
        "q_joint_mac_vs_joint_null": mac,
        "q_joint_rel_l2_vs_joint_null": q_rel_l2,
        "fem_relative_residual": float(fem_residual["relative_residual"]),
        "rayleigh_lambda": rayleigh_lambda,
        "notes": (
            f"{notes}; raw_abs={float(np.linalg.norm(raw_residual)):.16g}; "
            f"row_abs={float(np.linalg.norm(row_residual)):.16g}; "
            f"raw_row_components={raw_components}; row_normalized_components={row_components}; "
            f"raw_labeled_residual={raw_labeled}; row_normalized_labeled_residual={row_labeled}"
        ),
    }


def determinant_conditioning_rows(
    *,
    matrix: np.ndarray,
    c_det_raw: np.ndarray,
    c_joint: np.ndarray,
    q_joint_null: np.ndarray,
    Lambda: float,
    beta: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
    params,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    raw_singular_values = np.linalg.svd(np.asarray(matrix, dtype=float), compute_uv=False)
    rows.append(
        determinant_conditioning_method_row(
            method="reference_raw_svd",
            row_scaling_label="none",
            column_scaling_label="none",
            matrix=matrix,
            coeff=c_det_raw,
            singular_values=raw_singular_values,
            q_joint_null=q_joint_null,
            Lambda=Lambda,
            beta=beta,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            params=params,
            notes="current analytic_null_vector(raw determinant matrix)",
        )
    )
    rows.append(
        determinant_conditioning_method_row(
            method="reference_joint_null_coefficients",
            row_scaling_label="none",
            column_scaling_label="none",
            matrix=matrix,
            coeff=c_joint,
            singular_values=raw_singular_values,
            q_joint_null=q_joint_null,
            Lambda=Lambda,
            beta=beta,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            params=params,
            notes="coefficients reconstructed from q_joint_null_analytic_D in determinant basis",
        )
    )
    ones_row = np.ones(np.asarray(matrix, dtype=float).shape[0], dtype=float)
    ones_col = np.ones(np.asarray(matrix, dtype=float).shape[1], dtype=float)
    row_scale = row_unit_scale(matrix)
    method_specs: list[tuple[str, str, str, np.ndarray, np.ndarray]] = [
        ("scaled_raw_svd", "none", "none", ones_row, ones_col),
        ("scaled_row_normalized_svd", "unit_row_norm", "none", row_scale, ones_col),
        (
            "scaled_column_normalized_svd",
            "none",
            "unit_column_norm",
            ones_row,
            column_unit_scale(matrix),
        ),
        (
            "scaled_row_column_normalized_svd",
            "unit_row_norm",
            "unit_column_norm_after_row_scaling",
            row_scale,
            column_unit_scale(matrix, row_scale=row_scale),
        ),
    ]
    physical_scales = {
        "P_epsilon": float(epsilon),
        "P_inv_epsilon": 1.0 / float(epsilon),
        "P_Lambda": float(Lambda),
        "P_inv_Lambda": 1.0 / float(Lambda),
        "P_Lambda2": float(Lambda) ** 2,
        "P_inv_Lambda2": 1.0 / (float(Lambda) ** 2),
    }
    for label, p_scale in physical_scales.items():
        col_scale = np.array([1.0, 1.0, 1.0, 1.0, p_scale, p_scale], dtype=float)
        method_specs.append((f"scaled_physical_{label}", "none", label, ones_row, col_scale))
        method_specs.append((f"scaled_row_physical_{label}", "unit_row_norm", label, row_scale, col_scale))
    for method, row_label, col_label, method_row_scale, method_col_scale in method_specs:
        coeff_scaled, singular_values = scaled_svd_null_coeff(
            matrix,
            row_scale=method_row_scale,
            column_scale=method_col_scale,
        )
        rows.append(
            determinant_conditioning_method_row(
                method=method,
                row_scaling_label=row_label,
                column_scaling_label=col_label,
                matrix=matrix,
                coeff=coeff_scaled,
                singular_values=singular_values,
                q_joint_null=q_joint_null,
                Lambda=Lambda,
                beta=beta,
                beta_rad=beta_rad,
                mu=mu,
                epsilon=epsilon,
                params=params,
                notes="SVD null vector of scaled determinant matrix, then mapped back to raw coefficient basis",
            )
        )
    return rows


def determinant_conditioning_conclusion(rows: Sequence[dict[str, float | str]]) -> list[str]:
    joint_row = next((row for row in rows if str(row["method"]) == "reference_joint_null_coefficients"), None)
    raw_row = next((row for row in rows if str(row["method"]) == "reference_raw_svd"), None)
    lines: list[str] = []
    if joint_row is not None:
        joint_row_res = float(joint_row["row_normalized_matrix_residual"])
        if np.isfinite(joint_row_res) and joint_row_res <= 1e-8:
            lines.append("c_joint satisfies the current determinant matrix in row-normalized residual.")
        else:
            lines.append("c_joint does not satisfy the current determinant matrix closely; equation mismatch remains.")
    if raw_row is not None and joint_row is not None:
        raw_q = float(raw_row["q_joint_mac_vs_joint_null"])
        joint_q = float(joint_row["q_joint_mac_vs_joint_null"])
        if np.isfinite(joint_q) and joint_q >= 0.999999 and (not np.isfinite(raw_q) or raw_q < 0.999):
            lines.append("raw determinant SVD selects a different q_joint direction than the joint-null reference.")
    finite = [row for row in rows if np.isfinite(float(row["fem_relative_residual"]))]
    if finite:
        best_fem = min(finite, key=lambda row: float(row["fem_relative_residual"]))
        lines.append(
            "best FEM residual method: "
            f"{best_fem['method']} ({float(best_fem['fem_relative_residual']):.6e})"
        )
    finite_q = [row for row in rows if np.isfinite(float(row["q_joint_rel_l2_vs_joint_null"]))]
    if finite_q:
        best_q = min(finite_q, key=lambda row: float(row["q_joint_rel_l2_vs_joint_null"]))
        lines.append(
            "best q_joint match method: "
            f"{best_q['method']} (rel_l2={float(best_q['q_joint_rel_l2_vs_joint_null']):.6e})"
        )
    return lines


def print_top_conditioning_methods(rows: Sequence[dict[str, float | str]], *, limit: int = 5) -> None:
    finite_fem = [row for row in rows if np.isfinite(float(row["fem_relative_residual"]))]
    finite_q = [row for row in rows if np.isfinite(float(row["q_joint_rel_l2_vs_joint_null"]))]
    print("determinant nullspace conditioning audit:")
    print("  top by FEM relative residual:")
    for row in sorted(finite_fem, key=lambda item: float(item["fem_relative_residual"]))[:limit]:
        print(
            f"    {row['method']}: fem={float(row['fem_relative_residual']):.6e}, "
            f"q_rel_l2={float(row['q_joint_rel_l2_vs_joint_null']):.6e}, "
            f"row_res={float(row['row_normalized_matrix_residual']):.6e}"
        )
    print("  top by q_joint match:")
    for row in sorted(finite_q, key=lambda item: float(item["q_joint_rel_l2_vs_joint_null"]))[:limit]:
        print(
            f"    {row['method']}: q_rel_l2={float(row['q_joint_rel_l2_vs_joint_null']):.6e}, "
            f"mac={float(row['q_joint_mac_vs_joint_null']):.6e}, "
            f"fem={float(row['fem_relative_residual']):.6e}"
        )
    for line in determinant_conditioning_conclusion(rows):
        print(f"  conclusion: {line}")


def determinant_force_raw_to_physical_transform(
    *,
    Lambda: float,
    epsilon: float,
    mapping_variant: str,
) -> np.ndarray:
    lambda_sq = float(Lambda) ** 2
    force_scale = lambda_sq / float(epsilon)
    if mapping_variant == "current_scaled_rows_to_Fx_Fy_M":
        return np.array(
            [
                [0.0, 0.0, force_scale],
                [0.0, force_scale, 0.0],
                [lambda_sq, 0.0, 0.0],
            ],
            dtype=float,
        )
    if mapping_variant == "swap_Fx_Fy_scaled_rows":
        return np.array(
            [
                [0.0, force_scale, 0.0],
                [0.0, 0.0, force_scale],
                [lambda_sq, 0.0, 0.0],
            ],
            dtype=float,
        )
    raise ValueError(f"Unknown determinant force mapping variant: {mapping_variant}")


def determinant_as_joint_raw_matrix(
    *,
    Rf: np.ndarray,
    N: np.ndarray,
    QN_inverse: np.ndarray,
) -> np.ndarray:
    return np.asarray(Rf, dtype=float) @ np.asarray(N, dtype=float) @ np.asarray(QN_inverse, dtype=float)


def determinant_as_joint_physical_matrix(
    *,
    Rf: np.ndarray,
    N: np.ndarray,
    QN_inverse: np.ndarray,
    Lambda: float,
    epsilon: float,
    mapping_variant: str,
) -> np.ndarray:
    D_raw = determinant_as_joint_raw_matrix(Rf=Rf, N=N, QN_inverse=QN_inverse)
    transform = determinant_force_raw_to_physical_transform(
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        mapping_variant=mapping_variant,
    )
    return transform @ D_raw


def determinant_force_row_variants(matrix: np.ndarray) -> list[tuple[str, np.ndarray, str, str]]:
    Rf = np.asarray(matrix, dtype=float)[3:6, :]
    variants: list[tuple[str, np.ndarray, str, str]] = [
        (
            "current",
            Rf.copy(),
            "current_scaled_rows_to_Fx_Fy_M",
            "current determinant force rows mapped from [M/Lambda^2, eps*Fy/Lambda^2, eps*Fx/Lambda^2] to [Fx,Fy,M]",
        )
    ]
    mixed = Rf.copy()
    mixed[1, 5] *= -1.0
    mixed[2, 2:4] *= -1.0
    variants.append(
        (
            "flip_mixed_right_arm_signs_candidate",
            mixed,
            "current_scaled_rows_to_Fx_Fy_M",
            "flips candidate right-arm N2*sin(beta) term in transverse row and Q2*sin(beta) terms in axial row",
        )
    )
    axial = Rf.copy()
    axial[2, :] *= -1.0
    variants.append(
        (
            "flip_axial_row_sign",
            axial,
            "current_scaled_rows_to_Fx_Fy_M",
            "flips the determinant axial/global-x force row before physical row mapping",
        )
    )
    transverse = Rf.copy()
    transverse[1, :] *= -1.0
    variants.append(
        (
            "flip_transverse_row_sign",
            transverse,
            "current_scaled_rows_to_Fx_Fy_M",
            "flips the determinant transverse/global-y force row before physical row mapping",
        )
    )
    variants.append(
        (
            "swap_Fx_Fy_row_mapping",
            Rf.copy(),
            "swap_Fx_Fy_scaled_rows",
            "interprets determinant force rows as [M/Lambda^2, eps*Fx/Lambda^2, eps*Fy/Lambda^2]",
        )
    )
    return variants


def determinant_as_joint_matrix_rows(
    *,
    D_det_joint: np.ndarray,
    D_analytic_joint: np.ndarray,
    D_fem_joint: np.ndarray,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for row_idx, row_label in enumerate(DOF_LABELS):
        for col_idx, col_label in enumerate(JOINT_Q_LABELS):
            det_value = float(D_det_joint[row_idx, col_idx])
            analytic_value = float(D_analytic_joint[row_idx, col_idx])
            fem_value = float(D_fem_joint[row_idx, col_idx])
            rows.append(
                {
                    "row": row_label,
                    "col": col_label,
                    "D_det_joint_value": det_value,
                    "D_analytic_joint_value": analytic_value,
                    "D_fem_joint_value": fem_value,
                    "abs_diff_det_vs_analytic": float(abs(det_value - analytic_value)),
                    "rel_diff_det_vs_analytic": safe_relative_diff(det_value, analytic_value),
                    "abs_diff_det_vs_fem": float(abs(det_value - fem_value)),
                    "rel_diff_det_vs_fem": safe_relative_diff(det_value, fem_value),
                }
            )
    return rows


def smallest_right_singular_vector(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    _, singular_values, vh = np.linalg.svd(np.asarray(matrix, dtype=float), full_matrices=True)
    return normalized_vector(vh[-1, :].copy()), singular_values


def determinant_as_joint_summary_row(
    *,
    case_id: str,
    D_det_joint: np.ndarray,
    D_analytic_joint: np.ndarray,
    D_fem_joint: np.ndarray,
    condition_QN: float,
    notes: str,
) -> dict[str, float | str]:
    det_vec, det_s = smallest_right_singular_vector(D_det_joint)
    analytic_vec, analytic_s = smallest_right_singular_vector(D_analytic_joint)
    fem_vec, fem_s = smallest_right_singular_vector(D_fem_joint)
    mac_det_analytic, _, _ = aligned_unit_vector_metrics(det_vec, analytic_vec)
    mac_det_fem, _, _ = aligned_unit_vector_metrics(det_vec, fem_vec)
    return {
        "case_id": case_id,
        "det_joint_singular_values": format_float_sequence(det_s),
        "analytic_joint_singular_values": format_float_sequence(analytic_s),
        "fem_joint_singular_values": format_float_sequence(fem_s),
        "frobenius_rel_diff_det_vs_analytic": relative_matrix_norm(D_det_joint - D_analytic_joint, D_analytic_joint),
        "frobenius_rel_diff_det_vs_fem": relative_matrix_norm(D_det_joint - D_fem_joint, D_fem_joint),
        "smallest_singular_vector_MAC_det_vs_analytic": mac_det_analytic,
        "smallest_singular_vector_MAC_det_vs_fem": mac_det_fem,
        "condition_QN": float(condition_QN),
        "notes": notes,
    }


def determinant_as_joint_variant_rows(
    *,
    case_id: str,
    matrix: np.ndarray,
    N: np.ndarray,
    QN_inverse: np.ndarray,
    Lambda: float,
    epsilon: float,
    D_analytic_joint: np.ndarray,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    analytic_vec, _ = smallest_right_singular_vector(D_analytic_joint)
    for variant_name, Rf_variant, mapping_variant, note in determinant_force_row_variants(matrix):
        D_variant = determinant_as_joint_physical_matrix(
            Rf=Rf_variant,
            N=N,
            QN_inverse=QN_inverse,
            Lambda=float(Lambda),
            epsilon=float(epsilon),
            mapping_variant=mapping_variant,
        )
        diff = D_variant - np.asarray(D_analytic_joint, dtype=float)
        rel_fro = relative_matrix_norm(diff, D_analytic_joint)
        variant_vec, _ = smallest_right_singular_vector(D_variant)
        mac, _, _ = aligned_unit_vector_metrics(variant_vec, analytic_vec)
        for row_idx, row_label in enumerate(DOF_LABELS):
            for col_idx, col_label in enumerate(JOINT_Q_LABELS):
                det_value = float(D_variant[row_idx, col_idx])
                analytic_value = float(D_analytic_joint[row_idx, col_idx])
                rows.append(
                    {
                        "case_id": case_id,
                        "variant_name": variant_name,
                        "row": row_label,
                        "col": col_label,
                        "D_det_joint_variant_value": det_value,
                        "D_analytic_joint_value": analytic_value,
                        "abs_diff_det_vs_analytic": float(abs(det_value - analytic_value)),
                        "rel_diff_det_vs_analytic": safe_relative_diff(det_value, analytic_value),
                        "frobenius_rel_diff_det_vs_analytic": rel_fro,
                        "smallest_singular_vector_MAC_det_vs_analytic": mac,
                        "notes": f"{note}; physical rows=[Fx,Fy,M], cols=[ux,uy,theta]",
                    }
                )
    return rows


def determinant_as_joint_system_rows(
    *,
    case_id: str,
    matrix: np.ndarray,
    Lambda: float,
    beta_rad: float,
    mu: float,
    epsilon: float,
    D_analytic_joint: np.ndarray,
    D_fem_joint: np.ndarray,
) -> tuple[list[dict[str, float | str]], list[dict[str, float | str]], list[dict[str, float | str]]]:
    Rk = np.asarray(matrix, dtype=float)[:3, :]
    Rf = np.asarray(matrix, dtype=float)[3:6, :]
    N, singular_values, rank, tol, condition_Rk, kinematic_residual = kinematic_nullspace_basis(Rk)
    Q = joint_q_matrix(
        Lambda=float(Lambda),
        beta_rad=float(beta_rad),
        mu=float(mu),
        epsilon=float(epsilon),
        q_source_variant="left_joint_fem_theta",
    )
    QN = Q @ N
    if QN.shape == (3, 3):
        condition_QN = float(np.linalg.cond(QN))
        QN_inverse = np.linalg.inv(QN)
        qn_note = "QN inverted directly"
    else:
        condition_QN = float(np.linalg.cond(QN)) if QN.size else np.nan
        QN_inverse = np.linalg.pinv(QN)
        qn_note = "QN was not 3x3; Moore-Penrose pseudoinverse used"
    D_det_joint = determinant_as_joint_physical_matrix(
        Rf=Rf,
        N=N,
        QN_inverse=QN_inverse,
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        mapping_variant="current_scaled_rows_to_Fx_Fy_M",
    )
    notes = (
        "determinant rows restricted to null(Rk), then mapped from raw determinant force order "
        "[M/Lambda^2, eps*Fy/Lambda^2, eps*Fx/Lambda^2] to physical [Fx,Fy,M]; "
        f"Rk_rank={rank}; Rk_singular_values={format_float_sequence(singular_values)}; "
        f"Rk_tol={tol:.16g}; Rk_condition={condition_Rk:.16g}; "
        f"kinematic_residual_norm={kinematic_residual:.16g}; {qn_note}"
    )
    matrix_rows = determinant_as_joint_matrix_rows(
        D_det_joint=D_det_joint,
        D_analytic_joint=D_analytic_joint,
        D_fem_joint=D_fem_joint,
    )
    summary_rows = [
        determinant_as_joint_summary_row(
            case_id=case_id,
            D_det_joint=D_det_joint,
            D_analytic_joint=D_analytic_joint,
            D_fem_joint=D_fem_joint,
            condition_QN=condition_QN,
            notes=notes,
        )
    ]
    variant_rows = determinant_as_joint_variant_rows(
        case_id=case_id,
        matrix=matrix,
        N=N,
        QN_inverse=QN_inverse,
        Lambda=float(Lambda),
        epsilon=float(epsilon),
        D_analytic_joint=D_analytic_joint,
    )
    return matrix_rows, summary_rows, variant_rows


def print_determinant_as_joint_top_entries(
    matrix_rows: Sequence[dict[str, float | str]],
    variant_rows: Sequence[dict[str, float | str]],
    summary_rows: Sequence[dict[str, float | str]],
    *,
    limit: int = 5,
) -> None:
    print("determinant as q_joint system:")
    if summary_rows:
        summary = summary_rows[0]
        print(
            "  current mapping: "
            f"rel_fro_det_vs_analytic={float(summary['frobenius_rel_diff_det_vs_analytic']):.6e}, "
            f"rel_fro_det_vs_fem={float(summary['frobenius_rel_diff_det_vs_fem']):.6e}, "
            f"condition_QN={float(summary['condition_QN']):.6e}"
        )
        print(
            "  null-vector MAC: "
            f"det-vs-analytic={float(summary['smallest_singular_vector_MAC_det_vs_analytic']):.6e}, "
            f"det-vs-fem={float(summary['smallest_singular_vector_MAC_det_vs_fem']):.6e}"
        )
    print("  largest current entry differences vs analytic:")
    for row in sorted(matrix_rows, key=lambda item: float(item["abs_diff_det_vs_analytic"]), reverse=True)[:limit]:
        print(
            f"    {row['row']}<-{row['col']}: "
            f"det={float(row['D_det_joint_value']):.6e}, "
            f"analytic={float(row['D_analytic_joint_value']):.6e}, "
            f"abs={float(row['abs_diff_det_vs_analytic']):.6e}"
        )
    variant_summaries: dict[str, dict[str, float | str]] = {}
    for row in variant_rows:
        variant_summaries.setdefault(str(row["variant_name"]), row)
    print("  variants by Frobenius difference vs analytic:")
    for row in sorted(
        variant_summaries.values(),
        key=lambda item: float(item["frobenius_rel_diff_det_vs_analytic"]),
    )[:limit]:
        print(
            f"    {row['variant_name']}: "
            f"rel_fro={float(row['frobenius_rel_diff_det_vs_analytic']):.6e}, "
            f"null_MAC={float(row['smallest_singular_vector_MAC_det_vs_analytic']):.6e}"
        )


def largest_matrix_entry(rows: Sequence[dict[str, float | str]]) -> dict[str, float | str]:
    return max(rows, key=lambda row: float(row["abs_diff"]))


def row_family(row_label: str) -> str:
    if row_label == "Fx":
        return "axial/global_x"
    if row_label == "Fy":
        return "shear/global_y"
    if row_label == "M":
        return "moment"
    return "unknown"


def write_csv(path: Path, fieldnames: Sequence[str], rows: Sequence[dict[str, float | str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: Sequence[str] | None = None) -> dict[str, float | str]:
    args = parse_args(argv)
    branch_id = str(args.branch_id)
    beta = float(args.beta)
    mu = float(args.mu)
    epsilon = float(args.epsilon)
    params = build_params_for_case(epsilon=epsilon, l_total=float(args.l_total))

    tracking = resolve_analytic_branch(
        branch_id=branch_id,
        beta_deg=beta,
        mu=mu,
        epsilon=epsilon,
        num_analytic_roots=int(args.num_analytic_roots),
    )
    point = tracking["point"]
    Lambda = float(point.Lambda)
    omega_sq = Lambda**4
    beta_rad = float(np.deg2rad(beta))
    matrix = assemble_clamped_coupled_matrix(Lambda, beta_rad, mu, epsilon)
    coeff, _, _ = analytic_null_vector(matrix)
    matrix_residual = matrix @ coeff

    length_left = float(params.L_base) * (1.0 - mu)
    length_right = float(params.L_base) * (1.0 + mu)
    with fem_parameter_override(params):
        K_left, M_left = assemble_fem_arm_matrices(length=length_left, beta_deg=beta, rotate_to_global=False)
        K_right_local, M_right_local = assemble_fem_arm_matrices(
            length=length_right,
            beta_deg=beta,
            rotate_to_global=False,
        )
        K_right_global, M_right_global = assemble_fem_arm_matrices(
            length=length_right,
            beta_deg=beta,
            rotate_to_global=True,
        )

    n = fem.N_ELEM
    D_left_fem = schur_joint_dynamic_stiffness(
        K_left,
        M_left,
        omega_sq=omega_sq,
        joint_node=n,
        clamp_node=0,
    )
    D_right_fem_local = schur_joint_dynamic_stiffness(
        K_right_local,
        M_right_local,
        omega_sq=omega_sq,
        joint_node=0,
        clamp_node=n,
    )
    D_right_fem_global = schur_joint_dynamic_stiffness(
        K_right_global,
        M_right_global,
        omega_sq=omega_sq,
        joint_node=0,
        clamp_node=n,
    )
    D_total_fem = D_left_fem + D_right_fem_global

    D_left_analytic = analytic_arm_dynamic_stiffness_local(
        length=length_left,
        Lambda=Lambda,
        epsilon=epsilon,
        joint_at_right=True,
    )
    D_right_analytic_local = analytic_arm_dynamic_stiffness_local(
        length=length_right,
        Lambda=Lambda,
        epsilon=epsilon,
        joint_at_right=False,
    )
    T3 = rotation3(beta)
    D_right_analytic_global = T3 @ D_right_analytic_local @ T3.T
    D_total_analytic = D_left_analytic + D_right_analytic_global

    fields = analytic_local_fields(
        Lambda,
        mu=mu,
        epsilon=epsilon,
        coeff=coeff,
        s_norm=np.linspace(0.0, 1.0, fem.N_ELEM + 1),
        right_theta_sign=1.0,
    )
    q_full, _ = build_analytic_global_q(fields, beta_deg=beta)
    q_joint = q_full[3 * fem.N_ELEM : 3 * fem.N_ELEM + 3]
    dynamic_force_analytic = D_total_analytic @ q_joint
    dynamic_force_fem = D_total_fem @ q_joint

    matrices = {
        "left_arm": (D_left_fem, D_left_analytic),
        "right_arm_local": (D_right_fem_local, D_right_analytic_local),
        "right_arm_global": (D_right_fem_global, D_right_analytic_global),
        "total_joint": (D_total_fem, D_total_analytic),
    }
    matrix_rows = matrix_entry_rows(
        branch_id=branch_id,
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        Lambda=Lambda,
        omega_sq=omega_sq,
        matrices=matrices,
    )
    crosscheck_rows = determinant_crosscheck_rows(
        branch_id=branch_id,
        beta=beta,
        mu=mu,
        epsilon=epsilon,
        Lambda=Lambda,
        matrix_residual=matrix_residual,
        dynamic_force_analytic=dynamic_force_analytic,
        dynamic_force_fem=dynamic_force_fem,
    )
    case_id = (
        f"beta{filename_number_token(beta)}_{branch_id}_"
        f"mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}"
    )
    endpoint_matrix = independent_endpoint_matrix(
        Lambda=Lambda,
        beta_rad=beta_rad,
        mu=mu,
        epsilon=epsilon,
    )
    row_audit_rows: list[dict[str, float | str]] = []
    row_scaling_rows: list[dict[str, float | str]] = []
    contribution_rows: list[dict[str, float | str]] = []
    null_rows: list[dict[str, float | str]] = []
    row_audit_best_force_rows: list[dict[str, float | str]] = []
    constrained_rows: list[dict[str, float | str]] = []
    constrained_summary_rows: list[dict[str, float | str]] = []
    constrained_null_rows: list[dict[str, float | str]] = []
    force_scaling_rows: list[dict[str, float | str]] = []
    force_scaling_summary_rows: list[dict[str, float | str]] = []
    force_scaling_diagnostic_rows: list[dict[str, float | str]] = []
    baseline_dynamic_matrix = endpoint_matrix.copy()
    if bool(args.row_audit):
        endpoint_note = "independent explicit endpoint rows from cos/sin/cosh/sinh and axial sin/cos"
        row_audit_rows.extend(
            row_audit_matrix_rows(
                case_id=case_id,
                current=matrix,
                comparison=endpoint_matrix,
                matrix_pair="current_vs_independent_endpoint",
                force_order_variant="endpoint_explicit",
                sign_variant="as_written",
                notes=endpoint_note,
            )
        )
        row_scaling_rows.extend(
            row_scaling_audit_rows(
                case_id=case_id,
                current=matrix,
                comparison=endpoint_matrix,
                matrix_pair="current_vs_independent_endpoint",
                force_order_variant="endpoint_explicit",
                sign_variant="as_written",
                notes=endpoint_note,
            )
        )
        for q_variant in q_source_variants():
            q_matrix = joint_q_matrix(
                Lambda=Lambda,
                beta_rad=beta_rad,
                mu=mu,
                epsilon=epsilon,
                q_source_variant=q_variant,
            )
            for order_name in force_order_variants():
                for sign_label, signs in sign_variants():
                    dynamic_matrix = dynamic_stiffness_matrix_variant(
                        endpoint_matrix=endpoint_matrix,
                        D_total=D_total_analytic,
                        Lambda=Lambda,
                        epsilon=epsilon,
                        q_matrix=q_matrix,
                        force_order_variant=order_name,
                        sign_values=signs,
                    )
                    variant = f"q_source={q_variant};order={order_name}"
                    note = f"dynamic stiffness induced rows; {row_mapping_note()}"
                    row_audit_rows.extend(
                        row_audit_matrix_rows(
                            case_id=case_id,
                            current=matrix,
                            comparison=dynamic_matrix,
                            matrix_pair="current_vs_dynamic_stiffness",
                            force_order_variant=variant,
                            sign_variant=sign_label,
                            notes=note,
                        )
                    )
                    row_scaling_rows.extend(
                        row_scaling_audit_rows(
                            case_id=case_id,
                            current=matrix,
                            comparison=dynamic_matrix,
                            matrix_pair="current_vs_dynamic_stiffness",
                            force_order_variant=variant,
                            sign_variant=sign_label,
                            notes=note,
                        )
                    )
                    if (
                        q_variant == "left_joint_fem_theta"
                        and order_name == "M_epsFy_epsFx"
                        and sign_label == "+++"
                    ):
                        baseline_dynamic_matrix = dynamic_matrix
        row_audit_best_force_rows = best_scaling_rows(
            row_scaling_rows,
            matrix_pair="current_vs_dynamic_stiffness",
            force_only=True,
        )
        q_left_matrix = joint_q_matrix(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            q_source_variant="left_joint_fem_theta",
        )
        current_contrib = force_contribution_matrices_current(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
        )
        independent_contrib = force_contribution_matrices_current(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
        )
        dynamic_contrib = force_contribution_matrices_dynamic(
            D_left=D_left_analytic,
            D_right_local=D_right_analytic_local,
            q_matrix=q_left_matrix,
            beta_deg=beta,
            Lambda=Lambda,
            epsilon=epsilon,
        )
        contribution_rows = contribution_audit_rows(
            case_id=case_id,
            current_contrib=current_contrib,
            independent_contrib=independent_contrib,
            dynamic_contrib=dynamic_contrib,
            notes=f"baseline contribution audit uses q_source=left_joint_fem_theta; {row_mapping_note()}",
        )
        null_rows = null_vector_audit_rows(
            case_id=case_id,
            current=matrix,
            endpoint=endpoint_matrix,
            dynamic=baseline_dynamic_matrix,
            coeff=coeff,
            notes=f"dynamic column uses q_source=left_joint_fem_theta/order=M_epsFy_epsFx/sign=+++; {row_mapping_note()}",
        )
    constrained_best_summary: dict[str, float | str] | None = None
    if bool(args.constrained_row_audit):
        Rk = np.asarray(matrix[:3, :], dtype=float)
        Rf = np.asarray(matrix[3:6, :], dtype=float)
        nullspace_basis, singular_values, rank, tol, condition, kinematic_residual_norm = (
            kinematic_nullspace_basis(Rk)
        )
        q_matrix = joint_q_matrix(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            q_source_variant="average_joint_fem_theta",
        )
        q_left_matrix = joint_q_matrix(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            q_source_variant="left_joint_fem_theta",
        )
        q_right_matrix = joint_q_matrix(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            q_source_variant="right_joint_global_fem_theta",
        )
        q_subspace = q_matrix @ nullspace_basis
        q_mismatch_denom = float(np.linalg.norm(q_subspace))
        q_left_right_mismatch = (
            float(np.linalg.norm((q_left_matrix - q_right_matrix) @ nullspace_basis) / q_mismatch_denom)
            if q_mismatch_denom > NEAR_ZERO_NORM
            else np.nan
        )
        constrained_note = (
            "force rows compared after restricting coefficients to null(Rk); "
            "Q uses q_source=average_joint_fem_theta; "
            f"{row_mapping_note()}"
        )
        for order_name in force_order_variants():
            for sign_label, signs in sign_variants():
                variant_name = f"q_source=average_joint_fem_theta;order={order_name};sign={sign_label}"
                dynamic_operator = dynamic_force_row_operator(
                    D_total=D_total_analytic,
                    q_matrix=q_matrix,
                    Lambda=Lambda,
                    epsilon=epsilon,
                    force_order_variant=order_name,
                    sign_values=signs,
                )
                detail_rows, summary_row = constrained_force_audit_for_variant(
                    case_id=case_id,
                    variant_name=variant_name,
                    Rf=Rf,
                    nullspace_basis=nullspace_basis,
                    dynamic_operator=dynamic_operator,
                    singular_values=singular_values,
                    rank=rank,
                    tol=tol,
                    condition=condition,
                    kinematic_residual_norm=kinematic_residual_norm,
                    q_left_right_mismatch=q_left_right_mismatch,
                    notes=constrained_note,
                )
                constrained_rows.extend(detail_rows)
                constrained_summary_rows.append(summary_row)
                constrained_null_rows.extend(
                    constrained_null_vector_rows(
                        case_id=case_id,
                        variant_name=variant_name,
                        Rk=Rk,
                        Rf=Rf,
                        nullspace_basis=nullspace_basis,
                        dynamic_operator=dynamic_operator,
                        coeff=coeff,
                        notes=constrained_note,
                    )
                )
        finite_constrained = [
            row
            for row in constrained_summary_rows
            if np.isfinite(float(row["frobenius_rel_diff_after_row_scaling"]))
        ]
        if finite_constrained:
            constrained_best_summary = min(
                finite_constrained,
                key=lambda row: (
                    float(row["frobenius_rel_diff_after_row_scaling"]),
                    float(row["frobenius_rel_diff_before_scaling"]),
                ),
            )
    if bool(args.audit_force_row_scaling):
        current_contrib = force_contribution_matrices_current(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
        )
        fem_physical_contrib = force_contribution_matrices_fem_physical(
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
        )
        force_scaling_note = (
            "determinant_current_value is the current force-row contribution; "
            "fem_compatible_value is the unscaled physical generalized-force contribution "
            "with EI=1, EA=1/epsilon^2 in the FEM local/global force convention"
        )
        force_scaling_rows = force_row_scaling_audit_rows(
            case_id=case_id,
            current_contrib=current_contrib,
            fem_contrib=fem_physical_contrib,
            Lambda=Lambda,
            epsilon=epsilon,
            mu=mu,
            notes=force_scaling_note,
        )
        force_scaling_summary_rows = summarize_force_row_scaling(
            force_scaling_rows,
            case_id=case_id,
        )
        diagnostic_force_matrix = diagnostic_matrix_from_fem_physical_contributions(
            matrix=matrix,
            fem_contrib=fem_physical_contrib,
            Lambda=Lambda,
            epsilon=epsilon,
        )
        force_scaling_diagnostic_rows = diagnostic_force_matrix_rows(
            case_id=case_id,
            matrix=matrix,
            diagnostic_matrix=diagnostic_force_matrix,
            coeff=coeff,
            Lambda=Lambda,
            mu=mu,
            epsilon=epsilon,
            beta=beta,
            params=params,
        )
    joint_nullspace_matrix_data: list[dict[str, float | str]] = []
    joint_nullspace_summary_data: list[dict[str, float | str]] = []
    joint_nullspace_vector_data: list[dict[str, float | str]] = []
    joint_nullspace_pairwise_data: list[dict[str, float | str]] = []
    joint_reconstruction_residual_data: list[dict[str, float | str]] = []
    joint_nullspace_prefix = joint_nullspace_output_prefix(branch_id, beta, mu, epsilon)
    joint_nullspace_matrix_csv = Path(f"{joint_nullspace_prefix}_matrix.csv")
    joint_nullspace_summary_csv = Path(f"{joint_nullspace_prefix}_summary.csv")
    joint_nullspace_vectors_csv = Path(f"{joint_nullspace_prefix}_vectors.csv")
    joint_nullspace_pairwise_csv = Path(f"{joint_nullspace_prefix}_pairwise.csv")
    joint_reconstruction_residual_csv = REPO_ROOT / "results" / (
        f"joint_nullspace_reconstruction_residual_beta{filename_number_token(beta)}_{branch_id}"
        f"_mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}.csv"
    )
    joint_interpretation_lines: list[str] = []
    q_joint_null_analytic_for_audit: np.ndarray | None = None
    if bool(args.compare_joint_nullspace):
        q_joint_null_fem, fem_singular_values, _ = right_singular_null_vector(D_total_fem)
        q_joint_null_analytic, analytic_singular_values, _ = right_singular_null_vector(D_total_analytic)
        q_joint_null_analytic_for_audit = q_joint_null_analytic
        fem_branch = track_fem_branch(
            params,
            branch_id=branch_id,
            beta_deg=beta,
            mu=mu,
        )
        q_joint_from_fem_eigenvector = fem_branch.full_vec[3 * fem.N_ELEM : 3 * fem.N_ELEM + 3]
        joint_vectors = {
            "q_joint_from_FEM_eigenvector": q_joint_from_fem_eigenvector,
            "q_joint_from_analytic_determinant_coeffs": q_joint,
            "q_joint_null_FEM_D": q_joint_null_fem,
            "q_joint_null_analytic_D": q_joint_null_analytic,
        }
        joint_nullspace_matrix_data = joint_nullspace_matrix_rows(
            case_id=case_id,
            D_total_fem=D_total_fem,
            D_total_analytic=D_total_analytic,
        )
        joint_summary = joint_nullspace_summary_row(
            case_id=case_id,
            Lambda=Lambda,
            D_total_fem=D_total_fem,
            D_total_analytic=D_total_analytic,
            fem_singular_values=fem_singular_values,
            analytic_singular_values=analytic_singular_values,
        )
        joint_nullspace_summary_data = [joint_summary]
        joint_nullspace_vector_data = joint_vector_rows(case_id=case_id, vectors=joint_vectors)
        joint_nullspace_pairwise_data = joint_pairwise_rows(case_id=case_id, vectors=joint_vectors)
        joint_reconstruction_residual_data = joint_reconstruction_residual_rows(
            case_id=case_id,
            Lambda=Lambda,
            omega_sq=omega_sq,
            mu=mu,
            epsilon=epsilon,
            beta=beta,
            params=params,
            q_full_determinant=q_full,
            q_joint_determinant=q_joint,
            q_joint_null_analytic=q_joint_null_analytic,
        )
        joint_interpretation_lines = joint_nullspace_interpretation_lines(
            summary_row=joint_summary,
            pairwise_rows=joint_nullspace_pairwise_data,
        )
    determinant_conditioning_rows_data: list[dict[str, float | str]] = []
    determinant_conditioning_csv = REPO_ROOT / "results" / (
        f"determinant_nullspace_conditioning_beta{filename_number_token(beta)}_{branch_id}"
        f"_mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}.csv"
    )
    if bool(args.audit_determinant_nullspace_conditioning):
        if q_joint_null_analytic_for_audit is None:
            q_joint_null_analytic_for_audit, _, _ = right_singular_null_vector(D_total_analytic)
        c_joint = determinant_coefficients_from_joint_q(
            q_joint_global=q_joint_null_analytic_for_audit,
            Lambda=Lambda,
            beta=beta,
            mu=mu,
            epsilon=epsilon,
        )
        determinant_conditioning_rows_data = determinant_conditioning_rows(
            matrix=matrix,
            c_det_raw=coeff,
            c_joint=c_joint,
            q_joint_null=q_joint_null_analytic_for_audit,
            Lambda=Lambda,
            beta=beta,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            params=params,
        )
    determinant_as_joint_matrix_data: list[dict[str, float | str]] = []
    determinant_as_joint_summary_data: list[dict[str, float | str]] = []
    determinant_as_joint_variant_data: list[dict[str, float | str]] = []
    determinant_as_joint_prefix = REPO_ROOT / "results" / (
        f"determinant_as_joint_system_beta{filename_number_token(beta)}_{branch_id}"
        f"_mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}"
    )
    determinant_as_joint_matrix_csv = Path(f"{determinant_as_joint_prefix}_matrix.csv")
    determinant_as_joint_summary_csv = Path(f"{determinant_as_joint_prefix}_summary.csv")
    determinant_as_joint_variants_csv = REPO_ROOT / "results" / (
        f"determinant_as_joint_system_variants_beta{filename_number_token(beta)}_{branch_id}"
        f"_mu{filename_number_token(mu)}_eps{filename_number_token(epsilon)}.csv"
    )
    if bool(args.derive_determinant_as_joint_system):
        (
            determinant_as_joint_matrix_data,
            determinant_as_joint_summary_data,
            determinant_as_joint_variant_data,
        ) = determinant_as_joint_system_rows(
            case_id=case_id,
            matrix=matrix,
            Lambda=Lambda,
            beta_rad=beta_rad,
            mu=mu,
            epsilon=epsilon,
            D_analytic_joint=D_total_analytic,
            D_fem_joint=D_total_fem,
        )
    largest = largest_matrix_entry(matrix_rows)
    determinant_diff_norm = float(
        np.linalg.norm([float(row["abs_diff"]) for row in crosscheck_rows])
    )
    determinant_diff_max = float(max(float(row["abs_diff"]) for row in crosscheck_rows))

    summary = {
        "branch_id": branch_id,
        "beta": float(beta),
        "mu": float(mu),
        "epsilon": float(epsilon),
        "analytic_lambda": float(Lambda),
        "omega_sq": float(omega_sq),
        "current_sorted_index": int(point.current_sorted_index),
        "branch_tracking_warning_flag": str(tracking["warning_flag"]),
        "norm_D_left_diff": float(np.linalg.norm(D_left_fem - D_left_analytic)),
        "rel_norm_D_left_diff": relative_matrix_norm(D_left_fem - D_left_analytic, D_left_analytic),
        "norm_D_right_local_diff": float(np.linalg.norm(D_right_fem_local - D_right_analytic_local)),
        "rel_norm_D_right_local_diff": relative_matrix_norm(
            D_right_fem_local - D_right_analytic_local,
            D_right_analytic_local,
        ),
        "norm_D_right_global_diff": float(np.linalg.norm(D_right_fem_global - D_right_analytic_global)),
        "rel_norm_D_right_global_diff": relative_matrix_norm(
            D_right_fem_global - D_right_analytic_global,
            D_right_analytic_global,
        ),
        "norm_D_total_diff": float(np.linalg.norm(D_total_fem - D_total_analytic)),
        "rel_norm_D_total_diff": relative_matrix_norm(D_total_fem - D_total_analytic, D_total_analytic),
        "largest_entry_diff": float(largest["abs_diff"]),
        "largest_entry_matrix": str(largest["matrix_name"]),
        "largest_entry_row": str(largest["row"]),
        "largest_entry_col": str(largest["col"]),
        "largest_row_family": row_family(str(largest["row"])),
        "determinant_dynamic_scaled_max_abs_diff": determinant_diff_max,
        "determinant_dynamic_scaled_norm_diff": determinant_diff_norm,
        "q_joint_ux": float(q_joint[0]),
        "q_joint_uy": float(q_joint[1]),
        "q_joint_theta": float(q_joint[2]),
    }

    output_prefix = resolve_output_prefix(args.output_prefix, branch_id, beta, mu, epsilon)
    matrix_csv = suffixed_path(output_prefix, ".csv")
    summary_csv = suffixed_path(output_prefix, "_summary.csv")
    crosscheck_csv = suffixed_path(output_prefix, "_determinant_crosscheck.csv")
    row_audit_csv = suffixed_path(output_prefix, "_determinant_row_audit.csv")
    row_scaling_csv = suffixed_path(output_prefix, "_determinant_row_scaling_audit.csv")
    contribution_csv = suffixed_path(output_prefix, "_determinant_force_row_contribution_audit.csv")
    null_vector_csv = suffixed_path(output_prefix, "_determinant_null_vector_row_audit.csv")
    constrained_prefix = REPO_ROOT / "results" / f"constrained_determinant_force_row_audit_{case_id}"
    constrained_csv = Path(f"{constrained_prefix}.csv")
    constrained_summary_csv = Path(f"{constrained_prefix}_summary.csv")
    constrained_null_csv = REPO_ROOT / "results" / f"constrained_null_vector_force_check_{case_id}.csv"
    force_scaling_prefix = REPO_ROOT / "results" / f"force_row_scaling_audit_{case_id}"
    force_scaling_csv = Path(f"{force_scaling_prefix}.csv")
    force_scaling_summary_csv = Path(f"{force_scaling_prefix}_summary.csv")
    force_scaling_diagnostic_csv = Path(f"{force_scaling_prefix}_diagnostic_matrix.csv")
    write_csv(matrix_csv, MATRIX_FIELDNAMES, matrix_rows)
    write_csv(summary_csv, SUMMARY_FIELDNAMES, [summary])
    write_csv(crosscheck_csv, CROSSCHECK_FIELDNAMES, crosscheck_rows)
    if bool(args.row_audit):
        write_csv(row_audit_csv, ROW_AUDIT_FIELDNAMES, row_audit_rows)
        write_csv(row_scaling_csv, ROW_SCALING_AUDIT_FIELDNAMES, row_scaling_rows)
        write_csv(contribution_csv, FORCE_ROW_CONTRIBUTION_FIELDNAMES, contribution_rows)
        write_csv(null_vector_csv, NULL_VECTOR_ROW_AUDIT_FIELDNAMES, null_rows)
    if bool(args.constrained_row_audit):
        write_csv(constrained_csv, CONSTRAINED_ROW_AUDIT_FIELDNAMES, constrained_rows)
        write_csv(constrained_summary_csv, CONSTRAINED_SUMMARY_FIELDNAMES, constrained_summary_rows)
        write_csv(constrained_null_csv, CONSTRAINED_NULL_VECTOR_FIELDNAMES, constrained_null_rows)
    if bool(args.audit_force_row_scaling):
        write_csv(force_scaling_csv, FORCE_ROW_SCALING_AUDIT_FIELDNAMES, force_scaling_rows)
        write_csv(force_scaling_summary_csv, FORCE_ROW_SCALING_SUMMARY_FIELDNAMES, force_scaling_summary_rows)
        write_csv(
            force_scaling_diagnostic_csv,
            FORCE_ROW_SCALING_DIAGNOSTIC_FIELDNAMES,
            force_scaling_diagnostic_rows,
        )
    if bool(args.compare_joint_nullspace):
        write_csv(joint_nullspace_matrix_csv, JOINT_NULLSPACE_MATRIX_FIELDNAMES, joint_nullspace_matrix_data)
        write_csv(joint_nullspace_summary_csv, JOINT_NULLSPACE_SUMMARY_FIELDNAMES, joint_nullspace_summary_data)
        write_csv(joint_nullspace_vectors_csv, JOINT_NULLSPACE_VECTOR_FIELDNAMES, joint_nullspace_vector_data)
        write_csv(joint_nullspace_pairwise_csv, JOINT_NULLSPACE_PAIRWISE_FIELDNAMES, joint_nullspace_pairwise_data)
        write_csv(
            joint_reconstruction_residual_csv,
            JOINT_RECONSTRUCTION_RESIDUAL_FIELDNAMES,
            joint_reconstruction_residual_data,
        )
    if bool(args.audit_determinant_nullspace_conditioning):
        write_csv(
            determinant_conditioning_csv,
            DETERMINANT_NULLSPACE_CONDITIONING_FIELDNAMES,
            determinant_conditioning_rows_data,
        )
    if bool(args.derive_determinant_as_joint_system):
        write_csv(
            determinant_as_joint_matrix_csv,
            DETERMINANT_AS_JOINT_MATRIX_FIELDNAMES,
            determinant_as_joint_matrix_data,
        )
        write_csv(
            determinant_as_joint_summary_csv,
            DETERMINANT_AS_JOINT_SUMMARY_FIELDNAMES,
            determinant_as_joint_summary_data,
        )
        write_csv(
            determinant_as_joint_variants_csv,
            DETERMINANT_AS_JOINT_VARIANT_FIELDNAMES,
            determinant_as_joint_variant_data,
        )

    print("Joint dynamic stiffness analytic/FEM comparison")
    print(f"branch_id: {branch_id}")
    print(f"beta={beta:g} deg, mu={mu:g}, epsilon={epsilon:g}")
    print(f"Lambda={Lambda:.12g}, omega_sq=Lambda^4={omega_sq:.12g}")
    print(
        "relative Frobenius differences: "
        f"left={float(summary['rel_norm_D_left_diff']):.6e}, "
        f"right_local={float(summary['rel_norm_D_right_local_diff']):.6e}, "
        f"right_global={float(summary['rel_norm_D_right_global_diff']):.6e}, "
        f"total={float(summary['rel_norm_D_total_diff']):.6e}"
    )
    print(
        "largest entry diff: "
        f"{summary['largest_entry_matrix']}[{summary['largest_entry_row']},"
        f"{summary['largest_entry_col']}]={float(summary['largest_entry_diff']):.6e} "
        f"({summary['largest_row_family']})"
    )
    print(
        "determinant scaled force-row cross-check: "
        f"max_abs_diff={determinant_diff_max:.6e}, norm_diff={determinant_diff_norm:.6e}"
    )
    print(f"joint displacement order: q_joint = [ux, uy, theta_fem=dw/ds]")
    print(f"joint force order: f_joint = [Fx, Fy, M]")
    print(
        "determinant force-row mapping audited as: "
        "[moment_equilibrium, transverse_global_force_equilibrium, axial_global_force_equilibrium]"
    )
    if bool(args.row_audit):
        print("Best dynamic-stiffness row matches after per-row scalar fit:")
        for row in row_audit_best_force_rows:
            print(
                f"  {row['row_label']}: after={float(row['row_l2_relative_diff_after_best_scalar']):.6e}, "
                f"before={float(row['row_l2_relative_diff_before_scaling']):.6e}, "
                f"scalar={float(row['best_scalar']):.6e}, "
                f"corr={float(row['row_correlation']):.6e}, "
                f"{row['force_order_variant']}, sign={row['sign_variant']}"
            )
        print(f"determinant row audit CSV: {row_audit_csv}")
        print(f"row scaling audit CSV: {row_scaling_csv}")
        print(f"force-row contribution audit CSV: {contribution_csv}")
        print(f"null-vector row audit CSV: {null_vector_csv}")
    if bool(args.constrained_row_audit):
        if constrained_summary_rows:
            first = constrained_summary_rows[0]
            print(
                "constrained kinematic nullspace: "
                f"dim={int(first['subspace_dim'])}, "
                f"singular_values={first['Rk_singular_values']}, "
                f"q_left_right_mismatch={float(first['q_left_right_subspace_mismatch']):.6e}"
            )
        if constrained_best_summary is not None:
            print(
                "best constrained force-row variant after row scaling: "
                f"{constrained_best_summary['variant_name']}; "
                f"before={float(constrained_best_summary['frobenius_rel_diff_before_scaling']):.6e}, "
                f"after={float(constrained_best_summary['frobenius_rel_diff_after_row_scaling']):.6e}, "
                f"max_row_after={float(constrained_best_summary['max_row_error_after_scaling']):.6e}, "
                f"scalars={constrained_best_summary['best_row_scalars']}, "
                f"flag={constrained_best_summary['conclusion_flag']}"
            )
        print(f"constrained row audit CSV: {constrained_csv}")
        print(f"constrained summary CSV: {constrained_summary_csv}")
        print(f"constrained null-vector force check CSV: {constrained_null_csv}")
    if bool(args.audit_force_row_scaling):
        print("force-row scaling audit ratios det/current over FEM-compatible physical contribution:")
        for row in force_scaling_summary_rows:
            if int(row["nonzero_ratio_count"]) <= 0:
                continue
            print(
                f"  {row['row_label']} / {row['contribution_label']}: "
                f"pattern={row['ratio_pattern']}, factor={row['likely_missing_factor']}, "
                f"range=[{float(row['ratio_min']):.6e}, {float(row['ratio_max']):.6e}]"
            )
        for row in force_scaling_diagnostic_rows:
            print(
                "diagnostic common physical-force rows: "
                f"old-null force norm={float(row['diagnostic_force_on_current_null_norm']):.6e}, "
                f"new-null matrix residual={float(row['diagnostic_new_null_residual_norm']):.6e}, "
                f"new-null analytic-in-FEM relative residual="
                f"{float(row['diagnostic_new_null_analytic_in_fem_relative_residual']):.6e}"
            )
        print(f"force-row scaling audit CSV: {force_scaling_csv}")
        print(f"force-row scaling summary CSV: {force_scaling_summary_csv}")
        print(f"force-row scaling diagnostic matrix CSV: {force_scaling_diagnostic_csv}")
    if bool(args.compare_joint_nullspace):
        joint_summary_print = joint_nullspace_summary_data[0]
        print("joint displacement nullspace comparison:")
        print(
            "  D_total FEM/analytic: "
            f"frobenius_rel_diff={float(joint_summary_print['frobenius_rel_diff']):.6e}, "
            f"max_abs_diff={float(joint_summary_print['max_abs_diff']):.6e}"
        )
        print(
            "  singular values: "
            f"FEM=[{joint_summary_print['fem_singular_values']}], "
            f"analytic=[{joint_summary_print['analytic_singular_values']}]"
        )
        print(
            "  q_null FEM_D vs analytic_D: "
            f"MAC={pairwise_value(joint_nullspace_pairwise_data, 'q_joint_null_FEM_D', 'q_joint_null_analytic_D', 'mac'):.6e}, "
            "rel_l2="
            f"{pairwise_value(joint_nullspace_pairwise_data, 'q_joint_null_FEM_D', 'q_joint_null_analytic_D', 'relative_l2_after_sign_alignment'):.6e}"
        )
        print(
            "  q FEM eigenvector vs FEM Schur null: "
            f"MAC={pairwise_value(joint_nullspace_pairwise_data, 'q_joint_from_FEM_eigenvector', 'q_joint_null_FEM_D', 'mac'):.6e}, "
            "rel_l2="
            f"{pairwise_value(joint_nullspace_pairwise_data, 'q_joint_from_FEM_eigenvector', 'q_joint_null_FEM_D', 'relative_l2_after_sign_alignment'):.6e}"
        )
        print(
            "  q determinant coeffs vs analytic D null: "
            "MAC="
            f"{pairwise_value(joint_nullspace_pairwise_data, 'q_joint_from_analytic_determinant_coeffs', 'q_joint_null_analytic_D', 'mac'):.6e}, "
            "rel_l2="
            f"{pairwise_value(joint_nullspace_pairwise_data, 'q_joint_from_analytic_determinant_coeffs', 'q_joint_null_analytic_D', 'relative_l2_after_sign_alignment'):.6e}"
        )
        for line in joint_interpretation_lines:
            print(f"  interpretation: {line}")
        print(f"joint nullspace matrix CSV: {joint_nullspace_matrix_csv}")
        print(f"joint nullspace summary CSV: {joint_nullspace_summary_csv}")
        print(f"joint nullspace vectors CSV: {joint_nullspace_vectors_csv}")
        print(f"joint nullspace pairwise CSV: {joint_nullspace_pairwise_csv}")
        print(f"joint nullspace reconstruction residual CSV: {joint_reconstruction_residual_csv}")
    if bool(args.audit_determinant_nullspace_conditioning):
        print_top_conditioning_methods(determinant_conditioning_rows_data)
        print(f"determinant nullspace conditioning CSV: {determinant_conditioning_csv}")
    if bool(args.derive_determinant_as_joint_system):
        print_determinant_as_joint_top_entries(
            determinant_as_joint_matrix_data,
            determinant_as_joint_variant_data,
            determinant_as_joint_summary_data,
        )
        print(f"determinant as joint matrix CSV: {determinant_as_joint_matrix_csv}")
        print(f"determinant as joint summary CSV: {determinant_as_joint_summary_csv}")
        print(f"determinant as joint variants CSV: {determinant_as_joint_variants_csv}")
    print(f"matrix CSV: {matrix_csv}")
    print(f"summary CSV: {summary_csv}")
    print(f"determinant cross-check CSV: {crosscheck_csv}")
    return summary


if __name__ == "__main__":
    main()
